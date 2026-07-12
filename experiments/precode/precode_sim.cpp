// Standalone precode / dense-matrix-replacement simulator.
//
// This intentionally does not include or modify WirehairCodec.cpp.  It reuses
// WirehairTools row-weight generation and models the FULL intermediate-block
// linear system, unlike experiments/peeling/peel_sweep.cpp which models only
// the peel graph and reports residual sizes.
//
// System model per trial:
//   L = K source columns + D binary precode columns + H heavy columns.
//   Constraint rows: D binary precode rows (scheme-dependent pattern) plus
//   H heavy GF(256) rows.  Heavy rows are excluded from the binary system and
//   modeled as an MDS rank patch of up to H deficiencies (full-coverage
//   assumption; production covers only the last 18 GE columns, so production
//   failure rates are >= the ones reported here for the same def counts).
//   Received rows: K + OH rows, each with peel degree d over the K source
//   columns plus `mix` distinct columns over the D+H precode columns.
//
// After peeling, the simulator replays peeled-column substitution in solve
// order to build the exact GF(2) projection of every unused row onto the
// inactivated column set, then computes the binary rank of that residual.
// Decode success requires  peeled + residual_rank + heavy_patch == L.
//
// The self-test proves the projection replay correct via rank invariance:
// peeled + residual_rank must equal the brute-force GE rank of the entire
// binary system for any instance, because peeling is just a solving strategy.
//
// N is the number of blocks; payload bytes and memory pressure scale as
// N * block_bytes.  This harness counts block-XOR work in block units, so its
// cost columns are block-size independent; convert with xor_bench timings.

#include "WirehairTools.h"

#include "codec/WirehairV2Precode.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using namespace wirehair;

namespace {

//------------------------------------------------------------------------------
// RNG (splitmix64, same family as the v2 scaffold)

struct Rng
{
    uint64_t State;

    explicit Rng(uint64_t seed) : State(seed) {}

    uint64_t Next()
    {
        uint64_t z = (State += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }

    uint32_t U32()
    {
        return (uint32_t)(Next() >> 32);
    }

    // Unbiased uniform in [0, bound) via Lemire multiply-shift rejection
    uint32_t Below(uint32_t bound)
    {
        if (bound <= 1u) {
            return 0;
        }
        // 2^32 mod bound, computed in 32-bit arithmetic: 0u - bound wraps to
        // 2^32 - bound.  (Widening to 64 bits first would compute 2^64 mod
        // bound, which is the wrong rejection threshold for a 32-bit lane.)
        const uint32_t threshold = (0u - bound) % bound;
        for (;;)
        {
            const uint32_t x = U32();
            const uint64_t m = (uint64_t)x * bound;
            if ((uint32_t)m >= threshold) {
                return (uint32_t)(m >> 32);
            }
        }
    }

    double Unit()
    {
        return (Next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

static uint64_t Mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

static uint64_t HashString64(const std::string& s)
{
    uint64_t h = UINT64_C(0xcbf29ce484222325);
    for (char c : s)
    {
        h ^= (uint8_t)c;
        h *= UINT64_C(0x100000001b3);
    }
    return h;
}

//------------------------------------------------------------------------------
// Row degree distributions for received rows

enum class RowDist
{
    Wirehair,     // production GeneratePeelRowWeight
    LtM1C64,      // soliton-ish truncated LT, min degree 1, cap 64
    LtM1C16,      // same LT family, min degree 1, cap 16
    LtM2C1024,    // min degree 2, cap 1024 (residual-quality frontier family)
    RsC001D50C128, // robust soliton c=0.01 delta=0.50, min degree 1, cap 128
                  // (peel_sweep robust_soliton_weight law, renormalized)
    Fixed32,
    Fixed36,
    Fixed40,
    Fixed44,
    Fixed48,
    Fixed64,
    Fixed80,
    Fixed96,
    Fixed128,
    Fixed256,
    Mix32_128,
    Mix36_80P10,
    Mix36_80P25,
    Mix36_128P10,
    Mix8_64P50,
    Mix8_48P75,
    Mix12_64P50,
    Mix16_64P50,
    Mix16_96P25,
    Staged44_16P125,
    Staged44_16P167,
    Staged44_16P25,
    Staged44_12P167,
    Log6,
    Log8
};

static const char* RowDistName(RowDist dist)
{
    switch (dist)
    {
    case RowDist::Wirehair: return "wirehair";
    case RowDist::LtM1C64: return "lt_m1_c64";
    case RowDist::LtM1C16: return "lt_m1_c16";
    case RowDist::LtM2C1024: return "lt_m2_c1024";
    case RowDist::RsC001D50C128: return "rs_c001_d50_c128";
    case RowDist::Fixed32: return "fixed32";
    case RowDist::Fixed36: return "fixed36";
    case RowDist::Fixed40: return "fixed40";
    case RowDist::Fixed44: return "fixed44";
    case RowDist::Fixed48: return "fixed48";
    case RowDist::Fixed64: return "fixed64";
    case RowDist::Fixed80: return "fixed80";
    case RowDist::Fixed96: return "fixed96";
    case RowDist::Fixed128: return "fixed128";
    case RowDist::Fixed256: return "fixed256";
    case RowDist::Mix32_128: return "mix32_128";
    case RowDist::Mix36_80P10: return "mix36_80_p10";
    case RowDist::Mix36_80P25: return "mix36_80_p25";
    case RowDist::Mix36_128P10: return "mix36_128_p10";
    case RowDist::Mix8_64P50: return "mix8_64_p50";
    case RowDist::Mix8_48P75: return "mix8_48_p75";
    case RowDist::Mix12_64P50: return "mix12_64_p50";
    case RowDist::Mix16_64P50: return "mix16_64_p50";
    case RowDist::Mix16_96P25: return "mix16_96_p25";
    case RowDist::Staged44_16P125: return "staged44_16_p125";
    case RowDist::Staged44_16P167: return "staged44_16_p167";
    case RowDist::Staged44_16P25: return "staged44_16_p25";
    case RowDist::Staged44_12P167: return "staged44_12_p167";
    case RowDist::Log6: return "log6";
    case RowDist::Log8: return "log8";
    }
    return "unknown";
}

struct DegreeSampler
{
    RowDist Dist;
    unsigned K;
    std::vector<double> Cumulative; // for LT families
    double Total;
    unsigned MinDegree;
    unsigned MaxDegree;

    DegreeSampler(RowDist dist, unsigned k)
        : Dist(dist), K(k), Total(0.0), MinDegree(1), MaxDegree(1)
    {
        if (dist == RowDist::Wirehair) {
            return;
        }
        unsigned fixed = 0u;
        switch (dist)
        {
        case RowDist::Fixed32: fixed = 32u; break;
        case RowDist::Fixed36: fixed = 36u; break;
        case RowDist::Fixed40: fixed = 40u; break;
        case RowDist::Fixed44: fixed = 44u; break;
        case RowDist::Fixed48: fixed = 48u; break;
        case RowDist::Fixed64: fixed = 64u; break;
        case RowDist::Fixed80: fixed = 80u; break;
        case RowDist::Fixed96: fixed = 96u; break;
        case RowDist::Fixed128: fixed = 128u; break;
        case RowDist::Fixed256: fixed = 256u; break;
        case RowDist::Log6:
            fixed = (unsigned)std::ceil(6.0 * std::log2((double)k));
            break;
        case RowDist::Log8:
            fixed = (unsigned)std::ceil(8.0 * std::log2((double)k));
            break;
        default: break;
        }
        if (fixed != 0u)
        {
            MinDegree = std::min(fixed, k);
            MaxDegree = MinDegree;
            return;
        }
        if (dist == RowDist::Mix32_128)
        {
            MinDegree = std::min(32u, k);
            MaxDegree = std::min(128u, k);
            return;
        }
        if (dist == RowDist::Mix36_80P10 ||
            dist == RowDist::Mix36_80P25 ||
            dist == RowDist::Mix36_128P10)
        {
            MinDegree = std::min(36u, k);
            MaxDegree = std::min(
                dist == RowDist::Mix36_128P10 ? 128u : 80u, k);
            return;
        }
        if (dist == RowDist::Mix8_64P50 ||
            dist == RowDist::Mix8_48P75 ||
            dist == RowDist::Mix12_64P50 ||
            dist == RowDist::Mix16_64P50 ||
            dist == RowDist::Mix16_96P25)
        {
            unsigned low = 16u, high = 64u;
            if (dist == RowDist::Mix8_64P50 ||
                dist == RowDist::Mix8_48P75) {
                low = 8u;
            }
            else if (dist == RowDist::Mix12_64P50) {
                low = 12u;
            }
            if (dist == RowDist::Mix8_48P75) {
                high = 48u;
            }
            else if (dist == RowDist::Mix16_96P25) {
                high = 96u;
            }
            MinDegree = std::min(low, k);
            MaxDegree = std::min(high, k);
            return;
        }
        if (dist == RowDist::Staged44_16P125 ||
            dist == RowDist::Staged44_16P167 ||
            dist == RowDist::Staged44_16P25 ||
            dist == RowDist::Staged44_12P167)
        {
            MinDegree = std::min(
                dist == RowDist::Staged44_12P167 ? 12u : 16u, k);
            MaxDegree = std::min(44u, k);
            return;
        }
        unsigned min_d = 1, cap = 64;
        bool robust = false;
        if (dist == RowDist::LtM2C1024) {
            min_d = 2;
            cap = 1024;
        }
        else if (dist == RowDist::LtM1C16) {
            cap = 16;
        }
        else if (dist == RowDist::RsC001D50C128) {
            cap = 128;
            robust = true;
        }
        if (min_d > k) {
            min_d = k;
        }
        unsigned max_d = cap;
        if (max_d > k) {
            max_d = k;
        }
        if (max_d < min_d) {
            max_d = min_d;
        }
        MinDegree = min_d;
        MaxDegree = max_d;
        // Robust soliton tau-spike parameters; this mirrors peel_sweep's
        // robust_soliton_weight EXACTLY: R = max(1, c*ln(N/delta)*sqrt(N)),
        // spike = clamp(floor(N/R), 1, N), tau(d) = R/(d*N) for d < spike,
        // R*ln(R/delta)/N at d == spike, 0 above; the cap truncates the law
        // (the spike itself is truncated away when spike > cap) and the
        // cumulative sum renormalizes.
        const double rs_c = 0.01;
        const double rs_delta = 0.50;
        double rs_R = 0.0;
        unsigned rs_spike = 0;
        if (robust)
        {
            const double n = (double)k;
            rs_R = rs_c * std::log(n / rs_delta) * std::sqrt(n);
            if (rs_R < 1.0) {
                rs_R = 1.0;
            }
            rs_spike = (unsigned)std::floor(n / rs_R);
            if (rs_spike < 1u) {
                rs_spike = 1u;
            }
            if (rs_spike > k) {
                rs_spike = k;
            }
        }
        Cumulative.reserve(max_d - min_d + 1u);
        for (unsigned d = min_d; d <= max_d; ++d)
        {
            double w;
            if (d <= 1u) {
                w = 1.0 / (double)k;
            }
            else {
                w = 1.0 / ((double)d * (double)(d - 1u));
            }
            if (robust)
            {
                double tau = 0.0;
                if (d < rs_spike) {
                    tau = rs_R / ((double)d * (double)k);
                }
                else if (d == rs_spike) {
                    tau = rs_R * std::log(rs_R / rs_delta) / (double)k;
                }
                if (tau < 0.0) {
                    tau = 0.0;
                }
                w += tau;
            }
            Total += w;
            Cumulative.push_back(Total);
        }
    }

    unsigned Sample(Rng& rng) const
    {
        if (Dist == RowDist::Wirehair)
        {
            unsigned degree = GeneratePeelRowWeight(rng.U32(), (uint16_t)K);
            unsigned max_weight = K / 2u;
            if (max_weight < 1u) {
                max_weight = 1u;
            }
            if (degree > max_weight) {
                degree = max_weight;
            }
            if (degree < 1u) {
                degree = 1u;
            }
            return degree;
        }
        if (MinDegree == MaxDegree) {
            return MinDegree;
        }
        if (Dist == RowDist::Mix32_128) {
            return rng.Unit() < 0.5 ? MinDegree : MaxDegree;
        }
        if (Dist == RowDist::Mix36_80P10 ||
            Dist == RowDist::Mix36_80P25 ||
            Dist == RowDist::Mix36_128P10)
        {
            const double high_mass =
                Dist == RowDist::Mix36_80P25 ? 0.25 : 0.10;
            return rng.Unit() < high_mass ? MaxDegree : MinDegree;
        }
        if (Dist == RowDist::Mix8_64P50 ||
            Dist == RowDist::Mix8_48P75 ||
            Dist == RowDist::Mix12_64P50 ||
            Dist == RowDist::Mix16_64P50 ||
            Dist == RowDist::Mix16_96P25)
        {
            double high_mass = 0.50;
            if (Dist == RowDist::Mix8_48P75) {
                high_mass = 0.75;
            }
            else if (Dist == RowDist::Mix16_96P25) {
                high_mass = 0.25;
            }
            return rng.Unit() < high_mass ? MaxDegree : MinDegree;
        }
        if (Dist == RowDist::Staged44_16P125 ||
            Dist == RowDist::Staged44_16P167 ||
            Dist == RowDist::Staged44_16P25 ||
            Dist == RowDist::Staged44_12P167) {
            return MaxDegree;
        }
        const double target = rng.Unit() * Total;
        for (size_t i = 0; i < Cumulative.size(); ++i)
        {
            if (target <= Cumulative[i]) {
                return MinDegree + (unsigned)i;
            }
        }
        return MaxDegree;
    }
};

enum class PacketSchedule
{
    None,
    Iid,
    Burst,
    Permutation,
    RepairOnly
};

static const double kMaximumSupportedLoss = 0.99;

static const char* PacketScheduleName(PacketSchedule schedule)
{
    switch (schedule)
    {
    case PacketSchedule::None: return "none";
    case PacketSchedule::Iid: return "iid";
    case PacketSchedule::Burst: return "burst";
    case PacketSchedule::Permutation: return "permutation";
    case PacketSchedule::RepairOnly: return "repair-only";
    }
    return "unknown";
}

static std::vector<uint32_t> BuildPacketIds(
    uint32_t K,
    uint32_t delivered_count,
    double loss,
    uint64_t seed,
    PacketSchedule schedule)
{
    std::vector<uint32_t> output;
    output.reserve(delivered_count);
    Rng rng(seed ^ UINT64_C(0x10fade));
    uint64_t candidate = 0u;
    uint32_t burst_remaining = 0u;
    std::vector<uint32_t> permutation;
    size_t permutation_index = 0u;
    uint32_t permutation_base = 0u;
    const uint64_t candidate_limit =
        (uint64_t)delivered_count * 256u + 65536u;
    for (uint64_t attempts = 0u;
         output.size() < delivered_count && attempts < candidate_limit;
         ++attempts)
    {
        uint32_t id;
        if (schedule == PacketSchedule::RepairOnly) {
            id = K + (uint32_t)candidate++;
        }
        else if (schedule == PacketSchedule::Permutation)
        {
            if (permutation_index >= permutation.size())
            {
                const uint32_t count = std::min<uint32_t>(K + 512u, 65536u);
                permutation.resize(count);
                for (uint32_t i = 0u; i < count; ++i) {
                    permutation[i] = permutation_base + i;
                }
                for (uint32_t i = count; i > 1u; --i) {
                    std::swap(permutation[i - 1u],
                        permutation[rng.U32() % i]);
                }
                permutation_base += count;
                permutation_index = 0u;
            }
            id = permutation[permutation_index++];
        }
        else {
            id = (uint32_t)candidate++;
        }

        bool drop = false;
        if (schedule == PacketSchedule::Burst)
        {
            static const uint32_t kBurstLength = 8u;
            if (burst_remaining != 0u)
            {
                --burst_remaining;
                drop = true;
            }
            else
            {
                const double start_probability = loss /
                    (kBurstLength - (kBurstLength - 1u) * loss);
                if (rng.Unit() < start_probability)
                {
                    burst_remaining = kBurstLength - 1u;
                    drop = true;
                }
            }
        }
        else {
            drop = rng.Unit() < loss;
        }
        if (!drop) {
            output.push_back(id);
        }
    }
    return output;
}

//------------------------------------------------------------------------------
// Precode schemes

enum class SchemeKind
{
    None,        // no precode at all (pure LT control)
    DenseP50,    // D rows, iid p=0.5 over [0, K+D) -- production-dense idealization
    DenseSparse, // D rows, fixed weight W random columns over [0, K+D)
    LdpcStair,   // D parity columns; each source column connects to 3 parities;
                 // staircase double-diagonal among parity columns
    LdpcTri,     // RaptorQ-flavored: like stair but each source connects to 3
                 // parities chosen as {a, a+b, a+2b} mod D (circulant-ish)
    LdpcDense,   // S staircase parity columns + D2 dense p=0.5 rows over all
                 // K + S + D2 binary columns: staircase peels cheaply, dense
                 // rows pin the residual tail
    LdpcDenseShuffle2, // like LdpcDense but the D2 rows are generated
                 // Shuffle-2-style from a seeded column deck (production
                 // MultiplyDenseRows mechanics; see README for the exact
                 // rule): first row = half the deck, each next row = previous
                 // row with exactly 2 deck-driven bit flips
    CodecPort    // constraint rows imported from the actual codec-side
                 // construction (codec/WirehairV2Precode.cpp: PCGRandom
                 // staircase + unbiased Shuffle-2 deck), so the certified
                 // rank/failure machinery here validates the real port
};

struct Scheme
{
    std::string Name;
    SchemeKind Kind;
    unsigned D;       // binary precode rows/columns (staircase S for LdpcDense)
    unsigned H;       // heavy GF(256) rows/columns (MDS patch model)
    unsigned Weight;  // DenseSparse row weight
    unsigned Dense2;  // LdpcDense extra dense row/column count
    unsigned N1;      // staircase parities each SOURCE column connects to
                      // (LdpcStair/LdpcDense family; default 3 = historical)
    bool IdentCorner; // CodecPort only: identity-corner dense variant
                      // (deck excludes dense columns; row r owns column
                      // K+S+r) -- the encoder-feasible restructuring
};

//------------------------------------------------------------------------------
// Trial result

struct TrialResult
{
    bool Ok;                 // generation/solve consistency
    bool SuccessNoHeavy;     // def == 0
    bool Success;            // def <= H
    bool Runaway;            // --max-inact guard fired: peeling abandoned;
                             // counts as a decode failure, excluded from
                             // inact/cost means
    unsigned Def;            // residual column deficiency before heavy patch
    unsigned Inactivated;    // |I|
    unsigned ResidualRank;
    unsigned Peeled;
    unsigned ResidualRows;   // unused binary rows
    uint64_t RecvRowGenXors; // encoder per-received-row block XORs, summed
    uint64_t PrecodeGenXors; // encoder precode value generation block XORs
    uint64_t SparseSolveXors;// decoder peel-side block XOR proxy
    uint64_t BacksubXors;    // sum of popcount(proj) over peeled columns
    uint64_t GeBlockXors;    // ~R^2/2 dense block-op estimate
    uint64_t GeBitOps;       // residual_rows * R^2 / 64 word-op estimate
    uint64_t HeavyValueMuladds; // H*Inactivated + H*H block-unit GF(256) muladds
                                // (heavy value substitution + heavy GE elim)
    uint64_t HeavyDivs;      // H block-unit GF(256) divides
    // --ge-replay only (zero otherwise):
    uint64_t GeRealWordXors; // 64-bit words actually XORed during eliminations
    uint64_t GeRealRowOps;   // row eliminations performed
    uint64_t GeFillIn;       // pivot-row popcount growth above initial weight
    unsigned DefBand;        // deepest deficient col from-the-end position + 1
    bool DefOutsideW18;      // any deficient col with from-the-end pos >= 18
};

//------------------------------------------------------------------------------
// Sparse binary row container

struct SparseRow
{
    std::vector<uint32_t> Columns;
    bool IsConstraint;
};

//------------------------------------------------------------------------------
// Row generation

static void AddDistinctColumns(
    std::vector<uint32_t>& row,
    unsigned base,
    unsigned space,
    unsigned want,
    Rng& rng)
{
    if (want > space) {
        want = space;
    }
    if (want == 0u) {
        return;
    }
    if (want * 4u >= space)
    {
        // Partial Fisher-Yates over the whole space
        std::vector<uint32_t> deck(space);
        for (unsigned i = 0; i < space; ++i) {
            deck[i] = base + i;
        }
        for (unsigned i = 0; i < want; ++i)
        {
            const unsigned j = i + rng.Below(space - i);
            std::swap(deck[i], deck[j]);
            row.push_back(deck[i]);
        }
        return;
    }
    const size_t start = row.size();
    while (row.size() - start < want)
    {
        const uint32_t column = base + rng.Below(space);
        bool duplicate = false;
        for (size_t i = start; i < row.size(); ++i)
        {
            if (row[i] == column) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            row.push_back(column);
        }
    }
}

// Sattolo-style inside-out shuffle mirroring production ShuffleDeck16
// (WirehairTools.cpp:398): deck[0] = 0; for ii in [1, n): jj = rand % ii;
// deck[ii] = deck[jj]; deck[jj] = ii.  Production consumes its PCGRandom in
// 8-/16-bit chunks; the simulator draws jj as an unbiased full-width uniform
// in [0, ii) from the same Rng stream instead (documented deviation -- the
// codec implementation uses ShuffleDeck16 + PCGRandom directly).
static void Shuffle2Deck(std::vector<uint32_t>& deck, Rng& rng)
{
    if (deck.empty()) {
        return;
    }
    deck[0] = 0;
    for (uint32_t ii = 1; ii < (uint32_t)deck.size(); ++ii)
    {
        const uint32_t jj = rng.Below(ii);
        deck[ii] = deck[jj];
        deck[jj] = ii;
    }
}

struct GeneratedSystem
{
    unsigned L;
    unsigned K;
    unsigned D;
    unsigned H;
    std::vector<SparseRow> Rows;     // binary rows only (constraints first)
    uint64_t RecvRowGenXors;
    uint64_t PrecodeGenXors;
    bool PacketScheduleExhausted;
};

// paired_recv_seed (optional, --paired): seed for a SECOND RNG dedicated to
// received-row generation, whose value omits the scheme token, so different
// schemes at the same (K, oh, trial, --seed) draw an identical received-row
// stream.  Constraint rows always keep the scheme-dependent stream (`seed`).
static GeneratedSystem GenerateSystem(
    const Scheme& scheme,
    unsigned K,
    unsigned received,
    RowDist dist,
    unsigned mix,
    uint64_t seed,
    const uint64_t* paired_recv_seed = nullptr,
    PacketSchedule packet_schedule = PacketSchedule::None,
    double loss = 0.0,
    bool identity_systematic = false)
{
    GeneratedSystem sys;
    sys.K = K;
    sys.D = scheme.D + scheme.Dense2;
    sys.H = scheme.H;
    sys.L = K + sys.D + scheme.H;
    sys.RecvRowGenXors = 0;
    sys.PrecodeGenXors = 0;
    sys.PacketScheduleExhausted = false;

    Rng rng(seed);
    const unsigned precode_space = sys.D + scheme.H;

    // --- Binary precode constraint rows ---
    switch (scheme.Kind)
    {
    case SchemeKind::None:
        break;
    case SchemeKind::DenseP50:
    {
        // Each row: iid p=0.5 over the K+D source+dense columns.
        // Production uses the Shuffle-2 code, whose VALUE generation is
        // incremental: approximately 2.5*(K+D) block XORs total.
        for (unsigned r = 0; r < scheme.D; ++r)
        {
            SparseRow row;
            row.IsConstraint = true;
            row.Columns.reserve((K + scheme.D) / 2u + 8u);
            // Walk 64-bit words of random bits
            const unsigned span = K + scheme.D;
            for (unsigned base = 0; base < span; base += 64u)
            {
                uint64_t bits = rng.Next();
                const unsigned limit = (span - base < 64u) ? (span - base) : 64u;
                for (unsigned b = 0; b < limit; ++b)
                {
                    if (bits & 1u) {
                        row.Columns.push_back(base + b);
                    }
                    bits >>= 1;
                }
            }
            if (row.Columns.empty()) {
                row.Columns.push_back(rng.Below(span));
            }
            sys.Rows.push_back(std::move(row));
        }
        sys.PrecodeGenXors = (uint64_t)((K + scheme.D) * 5u) / 2u;
        break;
    }
    case SchemeKind::DenseSparse:
    {
        for (unsigned r = 0; r < scheme.D; ++r)
        {
            SparseRow row;
            row.IsConstraint = true;
            AddDistinctColumns(row.Columns, 0, K + scheme.D, scheme.Weight, rng);
            sys.Rows.push_back(std::move(row));
            if (scheme.Weight > 0u) {
                sys.PrecodeGenXors += scheme.Weight - 1u;
            }
        }
        break;
    }
    case SchemeKind::LdpcStair:
    case SchemeKind::LdpcTri:
    case SchemeKind::LdpcDense:
    case SchemeKind::LdpcDenseShuffle2:
    {
        const unsigned S = scheme.D;
        unsigned n1 = scheme.N1;
        if (n1 < 1u) {
            n1 = 1u;
        }
        if (n1 > 8u) {
            n1 = 8u;
        }
        std::vector<std::vector<uint32_t> > checks(S);
        uint32_t picks[8];
        for (unsigned c = 0; c < K; ++c)
        {
            if (scheme.Kind != SchemeKind::LdpcTri)
            {
                // n1 distinct random parities per source column (as many
                // distinct as S allows).  For n1 == 3 this draws EXACTLY the
                // same RNG stream as the historical hardcoded 3-hit loop, so
                // pre-n1 tokens keep generating identical systems.
                for (unsigned hit = 0; hit < n1; ++hit)
                {
                    uint32_t p = rng.Below(S);
                    if (S > hit)
                    {
                        for (bool collide = true; collide;)
                        {
                            collide = false;
                            for (unsigned j = 0; j < hit; ++j)
                            {
                                if (picks[j] == p)
                                {
                                    collide = true;
                                    p = rng.Below(S);
                                    break;
                                }
                            }
                        }
                    }
                    picks[hit] = p;
                    if (hit == 0u || S > hit) {
                        checks[p].push_back(c);
                    }
                }
            }
            else
            {
                // RFC-6330-flavored circulant triple {a, a+b, a+2b} mod S
                const uint32_t a = rng.Below(S);
                uint32_t b = 0;
                if (S > 1u) {
                    b = 1u + rng.Below(S - 1u);
                }
                checks[a].push_back(c);
                if (S > 1u) {
                    checks[(a + b) % S].push_back(c);
                    checks[(a + 2u * b) % S].push_back(c);
                }
            }
        }
        for (unsigned j = 0; j < S; ++j)
        {
            SparseRow row;
            row.IsConstraint = true;
            row.Columns = checks[j];
            row.Columns.push_back(K + j);          // own parity column
            if (j > 0u) {
                row.Columns.push_back(K + j - 1u); // staircase link
            }
            // Deduplicate source columns that landed twice on this check
            std::sort(row.Columns.begin(), row.Columns.end());
            std::vector<uint32_t> dedup;
            dedup.reserve(row.Columns.size());
            for (size_t i = 0; i < row.Columns.size();)
            {
                size_t j2 = i + 1u;
                while (j2 < row.Columns.size() &&
                       row.Columns[j2] == row.Columns[i]) {
                    ++j2;
                }
                if (((j2 - i) & 1u) != 0u) {
                    dedup.push_back(row.Columns[i]); // odd multiplicity survives GF(2)
                }
                i = j2;
            }
            row.Columns = std::move(dedup);
            if (row.Columns.size() > 1u) {
                sys.PrecodeGenXors += row.Columns.size() - 1u;
            }
            sys.Rows.push_back(std::move(row));
        }

        // LdpcDense: D2 extra dense p=0.5 rows over every binary column
        // (source + parity + dense), pinning the residual rank tail the
        // staircase alone cannot.
        if (scheme.Kind == SchemeKind::LdpcDense)
        {
            const unsigned span = K + S + scheme.Dense2;
            for (unsigned r = 0; r < scheme.Dense2; ++r)
            {
                SparseRow row;
                row.IsConstraint = true;
                row.Columns.reserve(span / 2u + 8u);
                for (unsigned base = 0; base < span; base += 64u)
                {
                    uint64_t bits = rng.Next();
                    const unsigned limit =
                        (span - base < 64u) ? (span - base) : 64u;
                    for (unsigned b = 0; b < limit; ++b)
                    {
                        if (bits & 1u) {
                            row.Columns.push_back(base + b);
                        }
                        bits >>= 1;
                    }
                }
                if (row.Columns.empty()) {
                    row.Columns.push_back(rng.Below(span));
                }
                sys.Rows.push_back(std::move(row));
            }
            // Shuffle-2-style incremental value generation estimate
            sys.PrecodeGenXors += (uint64_t)(span * 5u) / 2u;
        }
        // LdpcDenseShuffle2: D2 rows generated Shuffle-2-style from a single
        // column deck over all span = K + S + D2 binary columns, mirroring
        // production MultiplyDenseRows (WirehairCodec.cpp:786) as closely as
        // the bitset representation allows.  Exact rule (see README, must
        // stay codec-portable):
        //   1. deck = Shuffle2Deck over span columns (constraint-row Rng);
        //      set_count = ceil(span/2); first row = deck[0 .. set_count).
        //   2. Reshuffle deck; set half = deck[0 .. set_count), clear half =
        //      deck[set_count .. span).  Next floor(D2/2) rows: row[i+1] =
        //      row[i] XOR {set_half[ii], clear_half[ii]}, ii = 0,1,...
        //   3. Reshuffle deck again; remaining floor(D2/2) - 1 + (D2 & 1)
        //      rows by the same flip rule, ii restarting at 0.
        // Every row after the first differs from its predecessor in exactly
        // 2 columns (deck entries at distinct positions are distinct).  No
        // destination-row permutation deck: row order does not change the
        // linear system (production MAY add it for stream parity).
        else if (scheme.Kind == SchemeKind::LdpcDenseShuffle2 &&
                 scheme.Dense2 > 0u)
        {
            const unsigned span = K + S + scheme.Dense2;
            const unsigned set_count = (span + 1u) >> 1;
            std::vector<uint32_t> deck(span);
            Shuffle2Deck(deck, rng);
            std::vector<uint8_t> bitmap(span, 0);
            for (unsigned i = 0; i < set_count; ++i) {
                bitmap[deck[i]] = 1;
            }
            auto emit_row = [&]() {
                SparseRow row;
                row.IsConstraint = true;
                row.Columns.reserve(set_count + 8u);
                for (unsigned col = 0; col < span; ++col)
                {
                    if (bitmap[col]) {
                        row.Columns.push_back(col);
                    }
                }
                sys.Rows.push_back(std::move(row));
            };
            emit_row();
            // Incremental value generation: first-row accumulation, then two
            // block XORs per subsequent row
            sys.PrecodeGenXors += set_count;
            const unsigned halves[2] = {
                scheme.Dense2 >> 1,
                (scheme.Dense2 >> 1) + (scheme.Dense2 & 1u) - 1u
            };
            for (unsigned half = 0; half < 2u; ++half)
            {
                Shuffle2Deck(deck, rng); // Shuffle-2 reshuffle
                for (unsigned ii = 0; ii < halves[half]; ++ii)
                {
                    bitmap[deck[ii]] ^= 1;             // set-half flip
                    bitmap[deck[set_count + ii]] ^= 1; // clear-half flip
                    emit_row();
                    sys.PrecodeGenXors += 2u;
                }
            }
        }
        break;
    }
    case SchemeKind::CodecPort:
    {
        // Constraint rows from the real codec-side construction so the
        // rank/failure machinery here certifies the port itself.  The
        // codec uses its own PCGRandom streams seeded from `seed`; the
        // sim's constraint rng is deliberately unused in this branch
        // (received rows keep their own stream, paired or not).
        wirehair_v2::PrecodeParams params;
        params.BlockCount = K;
        params.Staircase = scheme.D;
        params.DenseRows = scheme.Dense2;
        params.HeavyRows = scheme.H;
        params.SourceHits = scheme.N1;
        params.DenseIdentityCorner = scheme.IdentCorner;
        params.Seed = seed;

        wirehair_v2::PrecodeSystem codec_sys;
        if (!wirehair_v2::BuildPrecodeSystem(params, codec_sys))
        {
            fprintf(stderr,
                "codecport: BuildPrecodeSystem rejected K=%u S=%u D2=%u N1=%u\n",
                K, scheme.D, scheme.Dense2, scheme.N1);
            exit(1);
        }

        for (const std::vector<uint32_t>& columns : codec_sys.StaircaseRows)
        {
            SparseRow row;
            row.IsConstraint = true;
            row.Columns = columns;
            if (row.Columns.size() > 1u) {
                sys.PrecodeGenXors += row.Columns.size() - 1u;
            }
            sys.Rows.push_back(std::move(row));
        }

        // Same incremental-generation accounting as LdpcDenseShuffle2:
        // first-row accumulation over known deck columns, then two block
        // XORs per subsequent row.  Identity-corner rows own one dense
        // output column each, which is not part of the known accumulation.
        const unsigned deck_span = scheme.IdentCorner ?
            K + scheme.D : K + scheme.D + scheme.Dense2;
        bool first_dense = true;
        for (const std::vector<uint32_t>& columns : codec_sys.DenseRowColumns)
        {
            SparseRow row;
            row.IsConstraint = true;
            row.Columns = columns;
            sys.PrecodeGenXors +=
                first_dense ? ((deck_span + 1u) >> 1) : 2u;
            first_dense = false;
            sys.Rows.push_back(std::move(row));
        }
        break;
    }
    }

    // --- Received rows ---
    // Paired mode: every RNG draw below is scheme-independent (degree sample,
    // source columns over [0, K), and exactly `mix` raw 64-bit canonical mix
    // draws per row -- consumed even when the precode space is empty), so two
    // schemes at the same (K, oh, trial, --seed) stay in stream lockstep.  The
    // canonical mix draws are mapped into the scheme-local D+H precode space
    // by modulo + deterministic linear probing (distinctness without extra
    // draws); the mapped columns coincide across schemes iff D+H coincides.
    const DegreeSampler degrees(dist, K);
    Rng paired_rng(paired_recv_seed ? *paired_recv_seed : 0);
    Rng& recv_rng = paired_recv_seed ? paired_rng : rng;
    const uint64_t packet_seed = paired_recv_seed ? *paired_recv_seed : seed;
    const std::vector<uint32_t> packet_ids =
        packet_schedule == PacketSchedule::None ? std::vector<uint32_t>() :
        BuildPacketIds(K, received, loss, packet_seed, packet_schedule);
    if (packet_schedule != PacketSchedule::None &&
        packet_ids.size() != received)
    {
        sys.PacketScheduleExhausted = true;
        return sys;
    }
    for (unsigned r = 0; r < received; ++r)
    {
        SparseRow row;
        row.IsConstraint = false;
        const uint32_t packet_id = packet_ids.empty() ? r : packet_ids[r];
        if (identity_systematic && packet_id < K)
        {
            row.Columns.push_back(packet_id);
            sys.Rows.push_back(std::move(row));
            continue;
        }
        Rng packet_rng(Mix64(
            packet_seed ^
            (uint64_t)packet_id * UINT64_C(0xd6e8feb86659fd93)));
        Rng& row_rng = packet_ids.empty() ? recv_rng : packet_rng;
        unsigned degree = degrees.Sample(row_rng);
        if (dist == RowDist::Staged44_16P125 ||
            dist == RowDist::Staged44_16P167 ||
            dist == RowDist::Staged44_16P25 ||
            dist == RowDist::Staged44_12P167)
        {
            const uint32_t repair_ordinal = packet_id >= K ? packet_id - K : 0u;
            uint32_t high_count = (K + 5u) / 6u;
            if (dist == RowDist::Staged44_16P125) {
                high_count = (K + 7u) / 8u;
            }
            else if (dist == RowDist::Staged44_16P25) {
                high_count = (K + 3u) / 4u;
            }
            degree = repair_ordinal < high_count ?
                degrees.MaxDegree : degrees.MinDegree;
        }
        AddDistinctColumns(row.Columns, 0, K, degree, row_rng);
        if (paired_recv_seed)
        {
            const size_t mix_start = row.Columns.size();
            for (unsigned m = 0; m < mix; ++m)
            {
                const uint64_t draw = row_rng.Next();
                if (precode_space == 0u ||
                    row.Columns.size() - mix_start >= precode_space) {
                    continue; // draw consumed anyway: keep streams in lockstep
                }
                const uint32_t offset = (uint32_t)(draw % precode_space);
                for (unsigned probe = 0; probe < precode_space; ++probe)
                {
                    const uint32_t column =
                        K + (offset + probe) % precode_space;
                    bool duplicate = false;
                    for (size_t i = mix_start; i < row.Columns.size(); ++i)
                    {
                        if (row.Columns[i] == column)
                        {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate)
                    {
                        row.Columns.push_back(column);
                        break;
                    }
                }
            }
        }
        else if (precode_space > 0u) {
            AddDistinctColumns(row.Columns, K, precode_space, mix, row_rng);
        }
        if (row.Columns.size() > 1u) {
            sys.RecvRowGenXors += row.Columns.size() - 1u;
        }
        sys.Rows.push_back(std::move(row));
    }
    return sys;
}

//------------------------------------------------------------------------------
// Peel solver with peel-order recording

struct PeelOutcome
{
    std::vector<uint8_t> ColumnState;       // 0=peeled 1=inactivated
    std::vector<uint32_t> SolveRow;         // for peeled columns
    std::vector<uint32_t> PeelOrder;        // columns in chronological solve order
    std::vector<uint32_t> InactOrder;       // columns in inactivation order
    uint64_t SparseSolveXors;
    bool Aborted;                           // --max-inact guard fired
};

static const uint8_t kStatePeeled = 0;
static const uint8_t kStateInactivated = 1;

class PeelSolver
{
public:
    PeelSolver(unsigned column_count, const std::vector<SparseRow>& rows)
        : ColumnCount(column_count),
          Rows(rows),
          Live(rows.size()),
          Used(rows.size(), 0),
          ColRows(column_count),
          Resolved(column_count, 0)
    {
        for (uint32_t r = 0; r < rows.size(); ++r)
        {
            Live[r] = (uint32_t)rows[r].Columns.size();
            for (uint32_t c : rows[r].Columns) {
                ColRows[c].push_back(r);
            }
            if (Live[r] == 1u) {
                Queue.push_back(r);
            }
        }
    }

    // max_inact == 0 means unlimited; otherwise the solve aborts (Aborted set)
    // as soon as the inactivated column count EXCEEDS max_inact.
    PeelOutcome Run(unsigned max_inact = 0)
    {
        PeelOutcome out;
        out.ColumnState.assign(ColumnCount, kStatePeeled);
        out.SolveRow.assign(ColumnCount, UINT32_MAX);
        out.SparseSolveXors = 0;
        out.Aborted = false;

        unsigned todo = ColumnCount;
        size_t queue_head = 0;
        while (todo > 0u)
        {
            // Exhaust singleton rows
            while (queue_head < Queue.size())
            {
                const uint32_t row_i = Queue[queue_head++];
                if (Live[row_i] != 1u || Used[row_i]) {
                    continue;
                }
                uint32_t column = UINT32_MAX;
                for (uint32_t c : Rows[row_i].Columns)
                {
                    if (!Resolved[c]) {
                        column = c;
                        break;
                    }
                }
                if (column == UINT32_MAX) {
                    continue;
                }
                Used[row_i] = 1;
                out.SolveRow[column] = row_i;
                out.PeelOrder.push_back(column);
                if (Rows[row_i].Columns.size() > 1u) {
                    out.SparseSolveXors += Rows[row_i].Columns.size() - 1u;
                }
                ResolveColumn(column);
                --todo;
            }
            if (todo == 0u) {
                break;
            }

            const uint32_t column = SelectInactivationColumn();
            out.ColumnState[column] = kStateInactivated;
            out.InactOrder.push_back(column);
            if (max_inact != 0u && out.InactOrder.size() > (size_t)max_inact)
            {
                out.Aborted = true;
                return out;
            }
            ResolveColumn(column);
            --todo;
        }
        return out;
    }

private:
    void ResolveColumn(uint32_t column)
    {
        Resolved[column] = 1;
        for (uint32_t row_i : ColRows[column])
        {
            if (Live[row_i] == 0u) {
                continue;
            }
            --Live[row_i];
            if (Live[row_i] == 1u && !Used[row_i]) {
                Queue.push_back(row_i);
            }
        }
    }

    // Production-flavored greedy: prefer the column covered by the most live
    // degree-2 rows (max avalanche), tie-break by total live references, then
    // lowest index.  Fallback: max live references.
    uint32_t SelectInactivationColumn()
    {
        std::vector<uint32_t> deg2(ColumnCount, 0);
        bool any_deg2 = false;
        for (uint32_t r = 0; r < Rows.size(); ++r)
        {
            if (Live[r] != 2u || Used[r]) {
                continue;
            }
            for (uint32_t c : Rows[r].Columns)
            {
                if (!Resolved[c])
                {
                    ++deg2[c];
                    any_deg2 = true;
                }
            }
        }

        uint32_t best = UINT32_MAX;
        uint64_t best_score = 0;
        for (uint32_t c = 0; c < ColumnCount; ++c)
        {
            if (Resolved[c]) {
                continue;
            }
            uint32_t live_refs = 0;
            for (uint32_t row_i : ColRows[c])
            {
                if (Live[row_i] > 0u && !Used[row_i]) {
                    ++live_refs;
                }
            }
            const uint64_t primary = any_deg2 ? deg2[c] : live_refs;
            const uint64_t score = (primary << 32) | live_refs;
            if (best == UINT32_MAX || score > best_score)
            {
                best = c;
                best_score = score;
            }
        }
        return best;
    }

    unsigned ColumnCount;
    const std::vector<SparseRow>& Rows;
    std::vector<uint32_t> Live;
    std::vector<uint8_t> Used;
    std::vector<std::vector<uint32_t> > ColRows;
    std::vector<uint8_t> Resolved;
    std::vector<uint32_t> Queue;
};

//------------------------------------------------------------------------------
// Residual projection + rank

struct ResidualAnalysis
{
    unsigned Inactivated;
    unsigned ResidualRows;
    unsigned Rank;
    unsigned Peeled;
    uint64_t BacksubXors;
};

class BitMatrix
{
public:
    explicit BitMatrix(unsigned bits)
        : Words((bits + 63u) / 64u), Bits(bits) {}

    unsigned Words;

    std::vector<uint64_t> NewRow() const
    {
        return std::vector<uint64_t>(Words, 0);
    }

    // Returns rank gained (0 or 1); reduces row in place against pivots.
    unsigned Insert(std::vector<uint64_t>& row)
    {
        for (;;)
        {
            int hi = -1;
            for (unsigned w = 0; w < Words; ++w)
            {
                if (row[w] != 0)
                {
                    hi = (int)(w * 64u + (unsigned)__builtin_ctzll(row[w]));
                    break;
                }
            }
            if (hi < 0) {
                return 0;
            }
            auto it = Pivots.find_slot((unsigned)hi);
            if (it == UINT32_MAX)
            {
                Pivots.add((unsigned)hi, (uint32_t)Stored.size());
                Stored.push_back(row);
                return 1;
            }
            const std::vector<uint64_t>& pivot = Stored[it];
            for (unsigned w = 0; w < Words; ++w) {
                row[w] ^= pivot[w];
            }
        }
    }

private:
    struct PivotMap
    {
        std::vector<uint32_t> Slot; // bit -> stored index

        uint32_t find_slot(unsigned bit) const
        {
            if (bit >= Slot.size()) {
                return UINT32_MAX;
            }
            return Slot[bit];
        }

        void add(unsigned bit, uint32_t index)
        {
            if (bit >= Slot.size()) {
                Slot.resize(bit + 1u, UINT32_MAX);
            }
            Slot[bit] = index;
        }
    };

    unsigned Bits;
    PivotMap Pivots;
    std::vector<std::vector<uint64_t> > Stored;
};

// residual_bits_out (optional): the exact per-unused-row projections onto the
// inactivated columns (bit i = i-th column in INACTIVATION order) that the
// rank accumulator consumes, captured before in-place reduction, for replay.
static ResidualAnalysis AnalyzeResidual(
    const GeneratedSystem& sys,
    const PeelOutcome& peel,
    std::vector<std::vector<uint64_t> >* residual_bits_out = nullptr)
{
    ResidualAnalysis res;
    res.Inactivated = (unsigned)peel.InactOrder.size();
    res.Peeled = (unsigned)peel.PeelOrder.size();
    res.BacksubXors = 0;

    const unsigned R = res.Inactivated;
    std::vector<uint32_t> inact_index(sys.L, UINT32_MAX);
    for (unsigned i = 0; i < R; ++i) {
        inact_index[peel.InactOrder[i]] = i;
    }

    const unsigned words = (R + 63u) / 64u;
    // Projection of each peeled column onto the inactivated set
    std::vector<std::vector<uint64_t> > proj(sys.L);

    std::vector<uint64_t> acc(words == 0 ? 1 : words);
    for (uint32_t column : peel.PeelOrder)
    {
        std::fill(acc.begin(), acc.end(), 0);
        const uint32_t row_i = peel.SolveRow[column];
        bool nonzero = false;
        for (uint32_t c : sys.Rows[row_i].Columns)
        {
            if (c == column) {
                continue;
            }
            if (inact_index[c] != UINT32_MAX)
            {
                acc[inact_index[c] >> 6] ^= UINT64_C(1) << (inact_index[c] & 63u);
                nonzero = true;
            }
            else if (!proj[c].empty())
            {
                for (unsigned w = 0; w < words; ++w) {
                    acc[w] ^= proj[c][w];
                }
                nonzero = true;
            }
        }
        if (nonzero && words > 0u)
        {
            bool all_zero = true;
            for (unsigned w = 0; w < words; ++w)
            {
                if (acc[w] != 0) {
                    all_zero = false;
                    break;
                }
            }
            if (!all_zero)
            {
                proj[column].assign(acc.begin(), acc.begin() + words);
                unsigned pop = 0;
                for (unsigned w = 0; w < words; ++w) {
                    pop += (unsigned)__builtin_popcountll(proj[column][w]);
                }
                res.BacksubXors += pop;
            }
        }
    }

    // Residual rows: every binary row never used as a peel pivot
    std::vector<uint8_t> used_row(sys.Rows.size(), 0);
    for (uint32_t column : peel.PeelOrder) {
        used_row[peel.SolveRow[column]] = 1;
    }

    BitMatrix rank_matrix(R == 0 ? 1 : R);
    unsigned residual_rows = 0;
    unsigned rank = 0;
    std::vector<uint64_t> rowbits(words == 0 ? 1 : words);
    for (uint32_t r = 0; r < sys.Rows.size(); ++r)
    {
        if (used_row[r]) {
            continue;
        }
        ++residual_rows;
        if (R == 0u) {
            continue;
        }
        std::fill(rowbits.begin(), rowbits.end(), 0);
        for (uint32_t c : sys.Rows[r].Columns)
        {
            if (inact_index[c] != UINT32_MAX)
            {
                rowbits[inact_index[c] >> 6] ^= UINT64_C(1) << (inact_index[c] & 63u);
            }
            else if (!proj[c].empty())
            {
                for (unsigned w = 0; w < words; ++w) {
                    rowbits[w] ^= proj[c][w];
                }
            }
        }
        if (residual_bits_out) {
            residual_bits_out->push_back(rowbits); // copy: Insert reduces in place
        }
        rank += rank_matrix.Insert(rowbits);
    }

    res.ResidualRows = residual_rows;
    res.Rank = rank;
    return res;
}

//------------------------------------------------------------------------------
// --ge-replay: explicit plain GE on the residual binary system
//
// Rows are the unused-row projections onto the inactivated columns (the same
// bitsets the rank accumulator consumes); columns are in INACTIVATION order.
// Rows are processed in ascending initial-popcount order (cheapest-pivot-first,
// a production-like heuristic; --ge-replay-reverse flips to descending as a
// sensitivity check) and the order is fixed up front, not re-sorted during
// elimination.  For each column j in order, the first remaining row with bit j
// normally becomes the pivot.  --ge-pivot-window=N instead scans the next N
// remaining rows and picks the candidate with the lowest remaining tail
// popcount, falling back to the first later candidate if the window has none.
// Counted per elimination: one row op, and (words - j/64) word XORs -- the
// 64-bit words a triangular implementation actually touches.  Fill-in is the
// popcount growth of a pivot row above its initial weight, measured at the
// moment it is selected as pivot.  Columns that end pivotless (deficient) are
// recorded by their from-the-end position (last inactivated column = 0).

struct GeReplayResult
{
    unsigned Rank;
    unsigned Def;
    uint64_t WordXors;       // 64-bit words XORed during eliminations
    uint64_t RowOps;         // row eliminations
    uint64_t FillIn;         // pivot popcount growth above initial weight
    unsigned DefBand;        // deepest deficient from-the-end position + 1 (0 = none)
    bool DefOutsideW18;      // any deficient from-the-end position >= 18
};

static unsigned TailPopcount(
    const std::vector<uint64_t>& row,
    unsigned word,
    uint64_t mask,
    unsigned words)
{
    unsigned pop = (unsigned)__builtin_popcountll(row[word] & ~(mask - 1u));
    for (unsigned w = word + 1u; w < words; ++w) {
        pop += (unsigned)__builtin_popcountll(row[w]);
    }
    return pop;
}

static GeReplayResult ReplayGe(
    std::vector<std::vector<uint64_t> > rows, // by value: consumed as workspace
    unsigned R,
    bool descending_weight,
    unsigned pivot_window = 0)
{
    GeReplayResult out;
    out.Rank = 0;
    out.Def = 0;
    out.WordXors = 0;
    out.RowOps = 0;
    out.FillIn = 0;
    out.DefBand = 0;
    out.DefOutsideW18 = false;
    if (R == 0u) {
        return out;
    }
    const unsigned words = (R + 63u) / 64u;
    const size_t n = rows.size();

    std::vector<uint32_t> init_pop(n, 0);
    for (size_t i = 0; i < n; ++i)
    {
        unsigned pop = 0;
        for (unsigned w = 0; w < words; ++w) {
            pop += (unsigned)__builtin_popcountll(rows[i][w]);
        }
        init_pop[i] = pop;
    }

    // Stable sort by initial weight: deterministic given the system row order
    std::vector<uint32_t> order(n);
    for (size_t i = 0; i < n; ++i) {
        order[i] = (uint32_t)i;
    }
    std::stable_sort(order.begin(), order.end(),
        [&](uint32_t a, uint32_t b) {
            return descending_weight ? init_pop[a] > init_pop[b]
                                     : init_pop[a] < init_pop[b];
        });
    std::vector<std::vector<uint64_t> > work;
    std::vector<uint32_t> pops;
    work.reserve(n);
    pops.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        work.push_back(std::move(rows[order[i]]));
        pops.push_back(init_pop[order[i]]);
    }

    size_t pivot_count = 0;
    for (unsigned j = 0; j < R; ++j)
    {
        const unsigned word = j >> 6;
        const uint64_t mask = UINT64_C(1) << (j & 63u);
        size_t found = SIZE_MAX;
        if (pivot_window > 0u)
        {
            const size_t window_end = std::min(
                n, pivot_count + (size_t)pivot_window);
            unsigned best_tail = UINT32_MAX;
            for (size_t r = pivot_count; r < window_end; ++r)
            {
                if (work[r][word] & mask)
                {
                    const unsigned tail = TailPopcount(work[r], word, mask, words);
                    if (tail < best_tail)
                    {
                        best_tail = tail;
                        found = r;
                        if (tail <= 1u) {
                            break;
                        }
                    }
                }
            }
        }
        const size_t search_begin =
            found == SIZE_MAX ? pivot_count : pivot_count + (size_t)pivot_window;
        for (size_t r = search_begin; found == SIZE_MAX && r < n; ++r)
        {
            if (work[r][word] & mask)
            {
                found = r;
                break;
            }
        }
        if (found == SIZE_MAX)
        {
            // Deficient column: record its from-the-end position
            const unsigned from_end = R - 1u - j;
            if (from_end + 1u > out.DefBand) {
                out.DefBand = from_end + 1u;
            }
            if (from_end >= 18u) {
                out.DefOutsideW18 = true;
            }
            ++out.Def;
            continue;
        }
        std::swap(work[found], work[pivot_count]);
        std::swap(pops[found], pops[pivot_count]);
        const std::vector<uint64_t>& pivot = work[pivot_count];

        unsigned pop_now = 0;
        for (unsigned w = 0; w < words; ++w) {
            pop_now += (unsigned)__builtin_popcountll(pivot[w]);
        }
        if (pop_now > pops[pivot_count]) {
            out.FillIn += pop_now - pops[pivot_count];
        }

        for (size_t r = pivot_count + 1u; r < n; ++r)
        {
            if (work[r][word] & mask)
            {
                for (unsigned w = word; w < words; ++w) {
                    work[r][w] ^= pivot[w];
                }
                ++out.RowOps;
                out.WordXors += words - word;
            }
        }
        ++pivot_count;
    }
    out.Rank = (unsigned)pivot_count;
    return out;
}

//------------------------------------------------------------------------------
// Brute-force binary rank of the entire system (self-test reference)

static unsigned BruteForceRank(const GeneratedSystem& sys)
{
    BitMatrix m(sys.L);
    unsigned rank = 0;
    for (const SparseRow& row : sys.Rows)
    {
        std::vector<uint64_t> bits = m.NewRow();
        for (uint32_t c : row.Columns) {
            bits[c >> 6] ^= UINT64_C(1) << (c & 63u);
        }
        rank += m.Insert(bits);
    }
    return rank;
}

//------------------------------------------------------------------------------
// One trial

static TrialResult RunTrial(
    const Scheme& scheme,
    unsigned K,
    unsigned overhead,
    RowDist dist,
    unsigned mix,
    uint64_t seed,
    const uint64_t* paired_recv_seed = nullptr,
    bool ge_replay = false,
    bool ge_replay_reverse = false,
    unsigned ge_pivot_window = 0,
    unsigned max_inact = 0,
    PacketSchedule packet_schedule = PacketSchedule::None,
    double loss = 0.0,
    bool identity_systematic = false)
{
    TrialResult t;
    std::memset(&t, 0, sizeof(t));

    // Received budget: the decoder collects K + overhead packets.  The D + H
    // constraint rows are free (known structure), so the binary system is
    // (K + overhead + D) rows over L = K + D + H columns, and L - K of the
    // unknowns are precode blocks the decoder does not return to the caller.
    const GeneratedSystem sys = GenerateSystem(
        scheme, K, K + overhead, dist, mix, seed, paired_recv_seed,
        packet_schedule, loss, identity_systematic);
    if (sys.PacketScheduleExhausted) {
        return t;
    }

    const PeelOutcome peel = PeelSolver(sys.L, sys.Rows).Run(max_inact);
    if (peel.Aborted)
    {
        // --max-inact runaway: peeling abandoned before completion, so no
        // residual rank (or cost ledger) exists.  Recorded as a failed trial;
        // the aggregator excludes it from inact/cost means.
        t.Ok = true;
        t.Runaway = true;
        t.Success = false;
        t.SuccessNoHeavy = false;
        t.Inactivated = (unsigned)peel.InactOrder.size();
        return t;
    }
    std::vector<std::vector<uint64_t> > residual_bits;
    const ResidualAnalysis res = AnalyzeResidual(
        sys, peel, ge_replay ? &residual_bits : nullptr);

    const unsigned R = res.Inactivated;
    const unsigned def = R - res.Rank;

    t.Ok = true;
    t.Def = def;
    t.SuccessNoHeavy = (def == 0u);
    t.Success = (def <= scheme.H);
    t.Inactivated = R;
    t.ResidualRank = res.Rank;
    t.Peeled = res.Peeled;
    t.ResidualRows = res.ResidualRows;
    t.RecvRowGenXors = sys.RecvRowGenXors;
    t.PrecodeGenXors = sys.PrecodeGenXors;
    t.SparseSolveXors = peel.SparseSolveXors;
    t.BacksubXors = res.BacksubXors;
    t.GeBlockXors = (uint64_t)R * R / 2u;
    t.GeBitOps = (uint64_t)res.ResidualRows * R * R / 64u;
    t.HeavyValueMuladds =
        (uint64_t)scheme.H * R + (uint64_t)scheme.H * scheme.H;
    t.HeavyDivs = scheme.H;

    if (ge_replay)
    {
        const GeReplayResult rep =
            ReplayGe(std::move(residual_bits), R, ge_replay_reverse,
                     ge_pivot_window);
        if (rep.Rank != res.Rank)
        {
            fprintf(stderr,
                "FATAL: --ge-replay rank %u != accumulator rank %u "
                "(scheme=%s K=%u oh=%u seed=0x%llx)\n",
                rep.Rank, res.Rank, scheme.Name.c_str(), K, overhead,
                (unsigned long long)seed);
            std::abort();
        }
        t.GeRealWordXors = rep.WordXors;
        t.GeRealRowOps = rep.RowOps;
        t.GeFillIn = rep.FillIn;
        t.DefBand = rep.DefBand;
        t.DefOutsideW18 = rep.DefOutsideW18;
    }
    return t;
}

//------------------------------------------------------------------------------
// Scheme parsing

static bool ParseUnsigned(const std::string& s, unsigned& out)
{
    // Require a leading digit: strtoul silently accepts whitespace and
    // sign prefixes, so "-1" would wrap to 4294967295
    if (s.empty() || s[0] < '0' || s[0] > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long v = strtoull(s.c_str(), &end, 10);
    if (!end || *end != '\0' || errno == ERANGE ||
        v > 0xffffffffull)
    {
        // Out-of-range values used to truncate to 32 bits, so e.g.
        // dense_d4294967296 ran with D=0 while the CSV kept the label
        return false;
    }
    out = (unsigned)v;
    return true;
}

static bool ParseU64(const std::string& s, uint64_t& out)
{
    if (s.empty() || s[0] < '0' || s[0] > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long v = strtoull(s.c_str(), &end, 0);
    if (!end || *end != '\0' || errno == ERANGE) {
        return false;
    }
    out = (uint64_t)v;
    return true;
}

static bool ParseNonNegativeDouble(const std::string& s, double& out)
{
    if (s.empty() ||
        ((s[0] < '0' || s[0] > '9') && s[0] != '.'))
    {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const double v = strtod(s.c_str(), &end);
    if (!end || *end != '\0' || errno == ERANGE ||
        !std::isfinite(v) || v < 0.0)
    {
        return false;
    }
    out = v;
    return true;
}

// Token grammar:
//   none
//   dense            (D = GetDenseCount(K), H = 6)
//   dense_d<D>       (explicit D)
//   densesparse_w<W> (D = GetDenseCount(K), weight W rows)
//   densesparse_w<W>_d<D>
//   ldpc             (S = GetDenseCount(K))
//   ldpc_s<S>
//   ldpc2x           (S = 2 * GetDenseCount(K))
//   ldpctri / ldpctri_s<S>
//   ldpcdense[_n1<X>]_s<S>_d<D2>[_s2]
//                    (X in {2,3,4}: staircase parities per source column,
//                     absent = 3 = historical behavior; _s2 = Shuffle-2
//                     structured D2 rows instead of iid p=0.5)
//   heavyonly        (D = 0, H = 16)
// Optional suffix on any token: _h<H> to override heavy row count
// (after _s2 when both are present: ..._s2_h<H>).
static bool MakeScheme(const std::string& token, unsigned K, Scheme& out)
{
    std::string body = token;
    out.H = 6;
    out.Weight = 0;
    out.Dense2 = 0;
    out.N1 = 3;
    out.IdentCorner = false;

    const size_t hpos = body.rfind("_h");
    if (hpos != std::string::npos)
    {
        unsigned h = 0;
        if (ParseUnsigned(body.substr(hpos + 2u), h))
        {
            out.H = h;
            body = body.substr(0, hpos);
        }
    }

    // _ic suffix (codecport family only): identity-corner dense variant
    // (deck excludes the dense columns; row r owns column K+S+r)
    if (body.rfind("codecport", 0) == 0 &&
        body.size() > 3u &&
        body.compare(body.size() - 3u, 3u, "_ic") == 0)
    {
        out.IdentCorner = true;
        body = body.substr(0, body.size() - 3u);
    }

    // _s2 suffix: Shuffle-2 structured D2 rows; only defined on the
    // ldpcdense family (the prefix check keeps ldpc_s2 meaning S = 2)
    bool shuffle2 = false;
    if (body.rfind("ldpcdense", 0) == 0 &&
        body.size() > 3u &&
        body.compare(body.size() - 3u, 3u, "_s2") == 0)
    {
        shuffle2 = true;
        body = body.substr(0, body.size() - 3u);
    }

    const unsigned dense_count = GetDenseCount((uint32_t)K);

    if (body == "none")
    {
        out.Kind = SchemeKind::None;
        out.D = 0;
        if (token.find("_h") == std::string::npos) {
            out.H = 0;
        }
    }
    else if (body == "dense")
    {
        out.Kind = SchemeKind::DenseP50;
        out.D = dense_count;
    }
    else if (body.rfind("dense_d", 0) == 0)
    {
        out.Kind = SchemeKind::DenseP50;
        if (!ParseUnsigned(body.substr(7), out.D)) {
            return false;
        }
    }
    else if (body.rfind("densesparse_w", 0) == 0)
    {
        out.Kind = SchemeKind::DenseSparse;
        std::string rest = body.substr(13);
        const size_t dpos = rest.find("_d");
        if (dpos == std::string::npos)
        {
            if (!ParseUnsigned(rest, out.Weight)) {
                return false;
            }
            out.D = dense_count;
        }
        else
        {
            if (!ParseUnsigned(rest.substr(0, dpos), out.Weight)) {
                return false;
            }
            if (!ParseUnsigned(rest.substr(dpos + 2u), out.D)) {
                return false;
            }
        }
    }
    else if (body == "ldpc")
    {
        out.Kind = SchemeKind::LdpcStair;
        out.D = dense_count;
    }
    else if (body == "ldpc2x")
    {
        out.Kind = SchemeKind::LdpcStair;
        out.D = dense_count * 2u;
    }
    else if (body.rfind("ldpc_s", 0) == 0)
    {
        out.Kind = SchemeKind::LdpcStair;
        if (!ParseUnsigned(body.substr(6), out.D)) {
            return false;
        }
    }
    else if (body.rfind("ldpcdense_n1", 0) == 0)
    {
        // ldpcdense_n1<X>_s<S>_d<D2>
        std::string rest = body.substr(12);
        const size_t spos = rest.find("_s");
        if (spos == std::string::npos) {
            return false;
        }
        unsigned n1 = 0;
        if (!ParseUnsigned(rest.substr(0, spos), n1)) {
            return false;
        }
        if (n1 < 2u || n1 > 4u) {
            return false; // E5 sweeps X in {2,3,4} only
        }
        out.N1 = n1;
        rest = rest.substr(spos + 2u);
        const size_t dpos = rest.find("_d");
        if (dpos == std::string::npos) {
            return false;
        }
        if (!ParseUnsigned(rest.substr(0, dpos), out.D)) {
            return false;
        }
        if (!ParseUnsigned(rest.substr(dpos + 2u), out.Dense2)) {
            return false;
        }
        out.Kind = shuffle2 ? SchemeKind::LdpcDenseShuffle2
                            : SchemeKind::LdpcDense;
    }
    else if (body.rfind("ldpcdense_s", 0) == 0)
    {
        out.Kind = shuffle2 ? SchemeKind::LdpcDenseShuffle2
                            : SchemeKind::LdpcDense;
        std::string rest = body.substr(11);
        const size_t dpos = rest.find("_d");
        if (dpos == std::string::npos) {
            return false;
        }
        if (!ParseUnsigned(rest.substr(0, dpos), out.D)) {
            return false;
        }
        if (!ParseUnsigned(rest.substr(dpos + 2u), out.Dense2)) {
            return false;
        }
    }
    else if (body == "ldpctri")
    {
        out.Kind = SchemeKind::LdpcTri;
        out.D = dense_count;
    }
    else if (body.rfind("ldpctri_s", 0) == 0)
    {
        out.Kind = SchemeKind::LdpcTri;
        if (!ParseUnsigned(body.substr(9), out.D)) {
            return false;
        }
    }
    else if (body == "heavyonly")
    {
        out.Kind = SchemeKind::None;
        out.D = 0;
        if (token.find("_h") == std::string::npos) {
            out.H = 16;
        }
    }
    else if (body == "codecport" || body.rfind("codecport_n1", 0) == 0)
    {
        // Certified rule as ported to codec/WirehairV2Precode.cpp:
        // S = GetDenseCount(K), D2 = 12, default H = 12.  N1 defaults to
        // the codec's banded verdict (2 below K=10000, 3 from K=10000
        // upward); codecport_n1<X> overrides it so old n12/n13 controls
        // remain reproducible on the real port.
        out.Kind = SchemeKind::CodecPort;
        out.D = dense_count;
        out.Dense2 = 12;
        out.N1 = K >= 10000u ? 3u : 2u;
        if (body != "codecport")
        {
            unsigned n1 = 0;
            if (!ParseUnsigned(body.substr(12), n1) || n1 < 1u || n1 > 8u) {
                return false;
            }
            out.N1 = n1;
        }
        if (token.find("_h") == std::string::npos) {
            out.H = 12;
        }
    }
    else
    {
        return false;
    }

    // Keep the simulator inside the production 16-bit column domain.  Apart
    // from preventing arithmetic wrap, these caps stop a malformed CLI token
    // from turning into an unbounded allocation before any CSV row is emitted.
    const uint64_t binary_columns =
        (uint64_t)out.D + (uint64_t)out.Dense2;
    const uint64_t total_columns =
        (uint64_t)K + binary_columns + (uint64_t)out.H;
    static const uint64_t kMaxModelPrecodeColumns = 4096;
    if (binary_columns > kMaxModelPrecodeColumns ||
        out.H > 128u ||
        total_columns > UINT16_MAX)
    {
        return false;
    }

    if (out.Kind == SchemeKind::DenseSparse &&
        (out.Weight == 0u || (uint64_t)out.Weight > (uint64_t)K + out.D))
    {
        return false;
    }

    // The LDPC family indexes parity buckets with rng.Below(S); S == 0 would
    // write through an empty checks vector in GenerateSystem
    if (out.D == 0u &&
        (out.Kind == SchemeKind::LdpcStair ||
         out.Kind == SchemeKind::LdpcTri ||
         out.Kind == SchemeKind::LdpcDense ||
         out.Kind == SchemeKind::LdpcDenseShuffle2 ||
         out.Kind == SchemeKind::CodecPort)) {
        return false;
    }

    out.Name = token;
    return true;
}

//------------------------------------------------------------------------------
// Aggregation

struct Aggregate
{
    uint64_t Trials = 0;
    uint64_t Fails = 0;
    uint64_t FailsNoHeavy = 0;
    uint64_t Runaways = 0;  // --max-inact aborts: in fail_rate, out of means
    uint64_t ScheduleExhaustions = 0; // bounded loss generator, out of means
    double DefSum = 0, DefMax = 0;
    double InactSum = 0, InactSqSum = 0, InactMax = 0;
    double RankSum = 0;
    double PeeledSum = 0;
    double ResidualRowsSum = 0;
    double RecvXorSum = 0;
    double PrecodeXorSum = 0;
    double SparseXorSum = 0;
    double BacksubXorSum = 0;
    double GeBlockXorSum = 0;
    double GeBitOpSum = 0;
    double HeavyMuladdSum = 0;
    double HeavyDivSum = 0;
    double GeRealWordXorSum = 0;
    double GeRealRowOpSum = 0;
    double FillInSum = 0;
    uint64_t DefOutsideW18Count = 0;
    std::vector<uint32_t> DefBands; // per-trial deepest def band (--ge-replay)
    std::vector<uint32_t> DefHist;  // def occurrence counts

    void Add(const TrialResult& t)
    {
        ++Trials;
        if (!t.Ok)
        {
            ++Fails;
            ++FailsNoHeavy;
            ++ScheduleExhaustions;
            return;
        }
        if (!t.Success) {
            ++Fails;
        }
        if (!t.SuccessNoHeavy) {
            ++FailsNoHeavy;
        }
        if (t.Runaway)
        {
            // Counts as a decode failure (above), but peeling was abandoned:
            // no def/inact/cost statistics exist for this trial.
            ++Runaways;
            return;
        }
        DefSum += t.Def;
        if (t.Def > DefMax) {
            DefMax = t.Def;
        }
        if (t.Def >= DefHist.size()) {
            DefHist.resize(t.Def + 1u, 0);
        }
        ++DefHist[t.Def];
        InactSum += t.Inactivated;
        InactSqSum += (double)t.Inactivated * t.Inactivated;
        if (t.Inactivated > InactMax) {
            InactMax = t.Inactivated;
        }
        RankSum += t.ResidualRank;
        PeeledSum += t.Peeled;
        ResidualRowsSum += t.ResidualRows;
        RecvXorSum += (double)t.RecvRowGenXors;
        PrecodeXorSum += (double)t.PrecodeGenXors;
        SparseXorSum += (double)t.SparseSolveXors;
        BacksubXorSum += (double)t.BacksubXors;
        GeBlockXorSum += (double)t.GeBlockXors;
        GeBitOpSum += (double)t.GeBitOps;
        HeavyMuladdSum += (double)t.HeavyValueMuladds;
        HeavyDivSum += (double)t.HeavyDivs;
        GeRealWordXorSum += (double)t.GeRealWordXors;
        GeRealRowOpSum += (double)t.GeRealRowOps;
        FillInSum += (double)t.GeFillIn;
        if (t.DefOutsideW18) {
            ++DefOutsideW18Count;
        }
        DefBands.push_back(t.DefBand);
    }
};

// Nearest-rank percentile over a SORTED vector; 0 for empty input
static unsigned PercentileSorted(const std::vector<uint32_t>& sorted, double p)
{
    if (sorted.empty()) {
        return 0;
    }
    const double rank = std::ceil(p * (double)sorted.size());
    size_t idx = (rank < 1.0) ? 0 : (size_t)rank - 1u;
    if (idx >= sorted.size()) {
        idx = sorted.size() - 1u;
    }
    return sorted[idx];
}

//------------------------------------------------------------------------------
// Self-test

// Self-test only: an INDEPENDENT transcription of the analytic degree laws
// (peel_sweep's lt_weight + robust_soliton_weight) so the chi-square clause
// can catch a mis-ported formula in DegreeSampler's cumulative table.
static double SelfTestDegreeWeight(RowDist dist, unsigned d, unsigned K)
{
    double w;
    if (d <= 1u) {
        w = 1.0 / (double)K;
    }
    else {
        w = 1.0 / ((double)d * (double)(d - 1u));
    }
    if (dist == RowDist::RsC001D50C128)
    {
        const double c = 0.01;
        const double delta = 0.50;
        const double n = (double)K;
        double R = c * std::log(n / delta) * std::sqrt(n);
        if (R < 1.0) {
            R = 1.0;
        }
        unsigned spike = (unsigned)std::floor(n / R);
        if (spike < 1u) {
            spike = 1u;
        }
        if (spike > K) {
            spike = K;
        }
        double tau = 0.0;
        if (d < spike) {
            tau = R / ((double)d * n);
        }
        else if (d == spike) {
            tau = R * std::log(R / delta) / n;
        }
        if (tau < 0.0) {
            tau = 0.0;
        }
        w += tau;
    }
    return w;
}

static bool SelfTest()
{
    // Rank invariance: peeled + residual rank must equal brute-force rank of
    // the whole binary system, for every scheme and instance.
    const char* tokens[] = {
        "none", "dense_d8", "densesparse_w6_d8", "ldpc_s5", "ldpctri_s5",
        "heavyonly_h4", "ldpcdense_s5_d4", "ldpcdense_s6_d3_h8",
        "ldpcdense_s5_d4_s2", "ldpcdense_n12_s5_d4", "ldpcdense_n14_s5_d4",
        "ldpcdense_n14_s5_d4_s2", "ldpcdense_n12_s6_d3_s2_h8",
        "codecport", "codecport_ic", "codecport_n13_ic"
    };
    struct ParseCase
    {
        const char* Token;
        unsigned K;
        bool Valid;
    };
    const ParseCase parse_cases[] = {
        { "none_h0", 32, true },
        { "none_h128", 32, true },
        { "none_h129", 32, false },
        { "none_h4294967295", 32, false },
        { "none_h4294967296", 32, false },
        { "dense_d0", 32, true },
        { "dense_d4096", 32, true },
        { "dense_d4097", 32, false },
        { "dense_d4294967295", 32, false },
        { "dense_d4294967296", 32, false },
        { "dense_d4096", 61433, true },
        { "dense_d4096", 61434, false },
        { "dense_d1529", 64000, true },
        { "dense_d1530", 64000, false },
        { "densesparse_w0_d8", 32, false },
        { "densesparse_w1_d0", 32, true },
        { "densesparse_w40_d8", 32, true },
        { "densesparse_w41_d8", 32, false },
        { "densesparse_w4294967295_d8", 32, false },
        { "ldpc_s0", 32, false },
        { "ldpc_s1", 32, true },
        { "ldpc_s4096", 32, true },
        { "ldpc_s4097", 32, false },
        { "ldpctri_s0", 32, false },
        { "ldpctri_s1", 32, true },
        { "ldpctri_s4096", 32, true },
        { "ldpctri_s4097", 32, false },
        { "ldpcdense_s0_d4", 32, false },
        { "ldpcdense_s1_d4095", 32, true },
        { "ldpcdense_s1_d4096", 32, false },
        { "ldpcdense_s4096_d0", 32, true },
        { "ldpcdense_s4096_d1", 32, false },
        { "ldpcdense_s0_d4_s2", 32, false },
        { "ldpcdense_s1_d4095_s2", 32, true },
        { "ldpcdense_n11_s5_d4", 32, false },
        { "ldpcdense_n12_s1_d4095", 32, true },
        { "ldpcdense_n14_s1_d4096", 32, false },
        { "ldpcdense_n15_s5_d4", 32, false },
        { "ldpcdense_n1_s5_d4", 32, false },
        { "ldpcdense_n13_d4", 32, false },
        { "heavyonly_h0", 32, true },
        { "heavyonly_h128", 32, true },
        { "heavyonly_h129", 32, false },
        { "codecport_n11", 32, true },
        { "codecport_n18_ic", 32, true },
        { "codecport_n10", 32, false },
        { "codecport_n19", 32, false },
        { "codecport_h128", 32, true },
        { "codecport_h129", 32, false },
    };
    for (const ParseCase& test : parse_cases)
    {
        Scheme scheme;
        const bool valid = MakeScheme(test.Token, test.K, scheme);
        if (valid != test.Valid) {
            fprintf(stderr,
                "self-test: scheme %s at K=%u was unexpectedly %s\n",
                test.Token,
                test.K,
                valid ? "accepted" : "rejected");
            return false;
        }
    }
    unsigned checked = 0;
    for (unsigned K = 8; K <= 48; K += 8)
    {
        for (const char* token : tokens)
        {
            Scheme scheme;
            if (!MakeScheme(token, K, scheme)) {
                fprintf(stderr, "self-test: bad scheme %s\n", token);
                return false;
            }
            for (unsigned oh = 0; oh <= 4; oh += 2)
            {
                for (unsigned trial = 0; trial < 40; ++trial)
                {
                    const uint64_t seed = Mix64(
                        HashString64(token) ^ Mix64(K * 977u + oh * 131u + trial));
                    const GeneratedSystem sys = GenerateSystem(
                        scheme, K, K + oh, RowDist::LtM1C64, 3, seed);

                    // No duplicate columns within any row
                    for (const SparseRow& row : sys.Rows)
                    {
                        std::vector<uint32_t> sorted = row.Columns;
                        std::sort(sorted.begin(), sorted.end());
                        for (size_t i = 1; i < sorted.size(); ++i)
                        {
                            if (sorted[i] == sorted[i - 1u])
                            {
                                fprintf(stderr,
                                    "self-test: duplicate column in %s K=%u\n",
                                    token, K);
                                return false;
                            }
                        }
                        for (uint32_t c : sorted)
                        {
                            if (c >= sys.L)
                            {
                                fprintf(stderr,
                                    "self-test: column out of range in %s\n",
                                    token);
                                return false;
                            }
                        }
                    }

                    const PeelOutcome peel = PeelSolver(sys.L, sys.Rows).Run();
                    if (peel.PeelOrder.size() + peel.InactOrder.size() != sys.L)
                    {
                        fprintf(stderr,
                            "self-test: column accounting broke in %s K=%u\n",
                            token, K);
                        return false;
                    }
                    std::vector<std::vector<uint64_t> > residual_bits;
                    const ResidualAnalysis res =
                        AnalyzeResidual(sys, peel, &residual_bits);
                    const unsigned via_peel =
                        (unsigned)peel.PeelOrder.size() + res.Rank;
                    const unsigned direct = BruteForceRank(sys);
                    if (via_peel != direct)
                    {
                        fprintf(stderr,
                            "self-test: rank invariance broke: %s K=%u oh=%u "
                            "trial=%u peel+residual=%u direct=%u\n",
                            token, K, oh, trial, via_peel, direct);
                        return false;
                    }

                    // GE replay (both row orders) must reproduce the
                    // accumulator rank on the same residual system
                    const unsigned R = (unsigned)peel.InactOrder.size();
                    const GeReplayResult fwd = ReplayGe(residual_bits, R, false);
                    const GeReplayResult rev = ReplayGe(residual_bits, R, true);
                    if (fwd.Rank != res.Rank || rev.Rank != res.Rank ||
                        fwd.Def != R - res.Rank || rev.Def != R - res.Rank)
                    {
                        fprintf(stderr,
                            "self-test: ge-replay rank broke: %s K=%u oh=%u "
                            "trial=%u accum=%u fwd=%u rev=%u\n",
                            token, K, oh, trial, res.Rank, fwd.Rank, rev.Rank);
                        return false;
                    }
                    ++checked;
                }
            }
        }
    }

    // Determinism + heavy-cost identity
    {
        Scheme scheme;
        MakeScheme("ldpc_s7", 64, scheme);
        const TrialResult a = RunTrial(scheme, 64, 2, RowDist::Wirehair, 3, 12345);
        const TrialResult b = RunTrial(scheme, 64, 2, RowDist::Wirehair, 3, 12345);
        if (a.Def != b.Def || a.Inactivated != b.Inactivated ||
            a.ResidualRank != b.ResidualRank ||
            a.SparseSolveXors != b.SparseSolveXors)
        {
            fprintf(stderr, "self-test: determinism broke\n");
            return false;
        }
        if (a.HeavyValueMuladds !=
                (uint64_t)scheme.H * a.Inactivated +
                (uint64_t)scheme.H * scheme.H ||
            a.HeavyDivs != scheme.H)
        {
            fprintf(stderr, "self-test: heavy cost identity broke\n");
            return false;
        }

        // Same determinism guarantee for the new kinds (Shuffle-2 D2 rows,
        // n1 variant, robust-soliton rowdist)
        Scheme s2scheme;
        if (!MakeScheme("ldpcdense_n12_s7_d4_s2", 64, s2scheme)) {
            fprintf(stderr, "self-test: s2 determinism scheme broke\n");
            return false;
        }
        const TrialResult c = RunTrial(
            s2scheme, 64, 2, RowDist::RsC001D50C128, 3, 777);
        const TrialResult d = RunTrial(
            s2scheme, 64, 2, RowDist::RsC001D50C128, 3, 777);
        if (c.Def != d.Def || c.Inactivated != d.Inactivated ||
            c.ResidualRank != d.ResidualRank ||
            c.SparseSolveXors != d.SparseSolveXors ||
            c.PrecodeGenXors != d.PrecodeGenXors)
        {
            fprintf(stderr, "self-test: s2/rowdist determinism broke\n");
            return false;
        }
    }

    // --paired: different schemes at equal sizes must generate identical
    // received-row degree sequences + source columns (and identical mix
    // columns when their D+H precode spaces coincide); without --paired the
    // received rows must differ for at least one instance.
    {
        const unsigned K = 32, received = 34;
        Scheme a, b, c;
        if (!MakeScheme("dense_d8", K, a) ||
            !MakeScheme("ldpc_s8", K, b) ||
            !MakeScheme("none", K, c))
        {
            fprintf(stderr, "self-test: paired scheme setup broke\n");
            return false;
        }
        auto source_cols = [&](const GeneratedSystem& sys,
                               std::vector<std::vector<uint32_t> >& out,
                               std::vector<std::vector<uint32_t> >& out_mix) {
            out.clear();
            out_mix.clear();
            for (const SparseRow& row : sys.Rows)
            {
                if (row.IsConstraint) {
                    continue;
                }
                std::vector<uint32_t> src, mixc;
                for (uint32_t col : row.Columns)
                {
                    if (col < K) {
                        src.push_back(col);
                    }
                    else {
                        mixc.push_back(col - K);
                    }
                }
                out.push_back(src);
                out_mix.push_back(mixc);
            }
        };
        bool any_unpaired_diff = false;
        for (unsigned trial = 0; trial < 8; ++trial)
        {
            const uint64_t seed_a =
                Mix64(HashString64("dense_d8") ^ Mix64(1000u + trial));
            const uint64_t seed_b =
                Mix64(HashString64("ldpc_s8") ^ Mix64(1000u + trial));
            const uint64_t seed_c =
                Mix64(HashString64("none") ^ Mix64(1000u + trial));
            const uint64_t recv_seed = Mix64(UINT64_C(0x9a17ed) + trial);

            const GeneratedSystem pa = GenerateSystem(
                a, K, received, RowDist::Wirehair, 3, seed_a, &recv_seed);
            const GeneratedSystem pb = GenerateSystem(
                b, K, received, RowDist::Wirehair, 3, seed_b, &recv_seed);
            const GeneratedSystem pc = GenerateSystem(
                c, K, received, RowDist::Wirehair, 3, seed_c, &recv_seed);
            std::vector<std::vector<uint32_t> > sa, sb, sc, ma, mb, mc;
            source_cols(pa, sa, ma);
            source_cols(pb, sb, mb);
            source_cols(pc, sc, mc);
            if (sa != sb || sa != sc)
            {
                fprintf(stderr,
                    "self-test: paired source columns differ (trial %u)\n",
                    trial);
                return false;
            }
            // dense_d8 (D+H = 8+6) and ldpc_s8 (8+6) share the precode space
            // width, so the canonical mix mapping must coincide exactly
            if (pa.D + pa.H == pb.D + pb.H && ma != mb)
            {
                fprintf(stderr,
                    "self-test: paired mix columns differ at equal D+H "
                    "(trial %u)\n", trial);
                return false;
            }

            const GeneratedSystem ua = GenerateSystem(
                a, K, received, RowDist::Wirehair, 3, seed_a);
            const GeneratedSystem ub = GenerateSystem(
                b, K, received, RowDist::Wirehair, 3, seed_b);
            std::vector<std::vector<uint32_t> > usa, usb, uma, umb;
            source_cols(ua, usa, uma);
            source_cols(ub, usb, umb);
            if (usa != usb) {
                any_unpaired_diff = true;
            }
        }
        if (!any_unpaired_diff)
        {
            fprintf(stderr,
                "self-test: unpaired schemes never diverged (suspicious)\n");
            return false;
        }
    }

    // Packet-ID schedules are part of the paired experiment contract: an
    // identical seed must reproduce the exact delivered-ID trace.  Every
    // schedule is duplicate-free here, IID with zero loss is sequential, and
    // repair-only must never expose a systematic ID.
    {
        const unsigned schedule_k = 1000;
        const unsigned delivered = 1100;
        const uint64_t schedule_seed = UINT64_C(0x1d5ced01);
        const PacketSchedule schedules[] = {
            PacketSchedule::Iid,
            PacketSchedule::Burst,
            PacketSchedule::Permutation,
            PacketSchedule::RepairOnly,
        };
        for (PacketSchedule schedule : schedules)
        {
            const double schedule_loss =
                schedule == PacketSchedule::RepairOnly ? 0.0 : 0.3;
            const std::vector<uint32_t> a = BuildPacketIds(
                schedule_k, delivered, schedule_loss, schedule_seed,
                schedule);
            const std::vector<uint32_t> b = BuildPacketIds(
                schedule_k, delivered, schedule_loss, schedule_seed,
                schedule);
            std::vector<uint32_t> sorted = a;
            std::sort(sorted.begin(), sorted.end());
            if (a != b || a.size() != delivered ||
                std::adjacent_find(sorted.begin(), sorted.end()) !=
                    sorted.end())
            {
                fprintf(stderr,
                    "self-test: packet schedule determinism/uniqueness broke\n");
                return false;
            }
            if (schedule == PacketSchedule::RepairOnly &&
                sorted.front() < schedule_k)
            {
                fprintf(stderr,
                    "self-test: repair-only schedule emitted source ID\n");
                return false;
            }
        }
        const std::vector<uint32_t> sequential = BuildPacketIds(
            schedule_k, delivered, 0.0, schedule_seed, PacketSchedule::Iid);
        for (unsigned i = 0; i < delivered; ++i)
        {
            if (sequential[i] != i)
            {
                fprintf(stderr,
                    "self-test: zero-loss IID schedule is not sequential\n");
                return false;
            }
        }
        const std::vector<uint32_t> extended = BuildPacketIds(
            schedule_k, delivered + 17u, 0.3, schedule_seed,
            PacketSchedule::Iid);
        const std::vector<uint32_t> prefix = BuildPacketIds(
            schedule_k, delivered, 0.3, schedule_seed,
            PacketSchedule::Iid);
        if (prefix.size() != delivered || extended.size() != delivered + 17u ||
            !std::equal(prefix.begin(), prefix.end(), extended.begin()))
        {
            fprintf(stderr,
                "self-test: overhead extension changed packet-ID prefix\n");
            return false;
        }

        Scheme prefix_scheme;
        if (!MakeScheme("codecport_ic", schedule_k, prefix_scheme)) {
            fprintf(stderr, "self-test: prefix scheme setup broke\n");
            return false;
        }
        const uint64_t system_seed = UINT64_C(0xc011ec7);
        const uint64_t receive_seed = UINT64_C(0xfacefeed);
        const GeneratedSystem short_system = GenerateSystem(
            prefix_scheme, schedule_k, delivered, RowDist::Fixed44, 3,
            system_seed, &receive_seed, PacketSchedule::Iid, 0.3, true);
        const GeneratedSystem long_system = GenerateSystem(
            prefix_scheme, schedule_k, delivered + 17u, RowDist::Fixed44, 3,
            system_seed, &receive_seed, PacketSchedule::Iid, 0.3, true);
        if (short_system.PacketScheduleExhausted ||
            long_system.PacketScheduleExhausted ||
            long_system.Rows.size() != short_system.Rows.size() + 17u)
        {
            fprintf(stderr, "self-test: overhead-prefix system setup broke\n");
            return false;
        }
        for (size_t i = 0; i < short_system.Rows.size(); ++i)
        {
            if (short_system.Rows[i].IsConstraint !=
                    long_system.Rows[i].IsConstraint ||
                short_system.Rows[i].Columns != long_system.Rows[i].Columns)
            {
                fprintf(stderr,
                    "self-test: overhead extension changed system prefix\n");
                return false;
            }
        }

        const GeneratedSystem exhausted = GenerateSystem(
            prefix_scheme, schedule_k, delivered, RowDist::Fixed44, 3,
            system_seed, &receive_seed, PacketSchedule::Iid, 1.0, true);
        const TrialResult exhausted_trial = RunTrial(
            prefix_scheme, schedule_k, 0u, RowDist::Fixed44, 3,
            system_seed, &receive_seed, false, false, 0u, 0u,
            PacketSchedule::Iid, 1.0, true);
        if (!exhausted.PacketScheduleExhausted || exhausted_trial.Ok)
        {
            fprintf(stderr,
                "self-test: exhausted packet schedule was not structured\n");
            return false;
        }
    }

    // Identity-systematic packet rows must be unit rows for source IDs and
    // must use the requested repair law for repair IDs.  Also exercise the
    // production PeelRowIterator with an overridden fixed degree: it must be
    // deterministic and visit distinct in-range columns for every packet ID.
    {
        const unsigned identity_k = 1000;
        const unsigned received = 1100;
        const uint64_t recv_seed = UINT64_C(0x1d3a71c0);
        Scheme scheme;
        if (!MakeScheme("codecport_ic", identity_k, scheme))
        {
            fprintf(stderr, "self-test: identity scheme setup broke\n");
            return false;
        }
        const GeneratedSystem sys = GenerateSystem(
            scheme, identity_k, received, RowDist::Fixed44, 3,
            UINT64_C(0x51757e), &recv_seed, PacketSchedule::Iid, 0.0, true);
        unsigned packet_row = 0;
        for (const SparseRow& row : sys.Rows)
        {
            if (row.IsConstraint) {
                continue;
            }
            unsigned source_refs = 0;
            std::vector<uint32_t> source_columns;
            for (uint32_t column : row.Columns)
            {
                if (column < identity_k)
                {
                    ++source_refs;
                    source_columns.push_back(column);
                }
            }
            std::sort(source_columns.begin(), source_columns.end());
            if (packet_row < identity_k)
            {
                if (row.Columns.size() != 1u ||
                    row.Columns.front() != packet_row)
                {
                    fprintf(stderr,
                        "self-test: identity source row %u is not singleton\n",
                        packet_row);
                    return false;
                }
            }
            else if (source_refs != 44u ||
                std::adjacent_find(source_columns.begin(),
                    source_columns.end()) != source_columns.end())
            {
                fprintf(stderr,
                    "self-test: fixed44 repair row %u is malformed\n",
                    packet_row);
                return false;
            }
            ++packet_row;
        }
        if (packet_row != received)
        {
            fprintf(stderr, "self-test: identity packet row count broke\n");
            return false;
        }

        const GeneratedSystem staged = GenerateSystem(
            scheme, identity_k, identity_k, RowDist::Staged44_16P167, 3,
            UINT64_C(0x51757e), &recv_seed,
            PacketSchedule::RepairOnly, 0.0, true);
        packet_row = 0;
        const unsigned staged_high_count = (identity_k + 5u) / 6u;
        for (const SparseRow& row : staged.Rows)
        {
            if (row.IsConstraint) {
                continue;
            }
            unsigned source_refs = 0;
            for (uint32_t column : row.Columns) {
                source_refs += column < identity_k;
            }
            const unsigned expected =
                packet_row < staged_high_count ? 44u : 16u;
            if (source_refs != expected)
            {
                fprintf(stderr,
                    "self-test: staged repair row %u has degree %u, want %u\n",
                    packet_row, source_refs, expected);
                return false;
            }
            ++packet_row;
        }

        const uint16_t source_prime = wirehair::NextPrime16(identity_k);
        std::vector<uint64_t> marginal_hits(identity_k, 0u);
        for (uint32_t packet_id = 0; packet_id < 20000u; ++packet_id)
        {
            wirehair::PeelRowParameters a, b;
            a.Initialize(packet_id, 0x9e3779b9u, identity_k, 18u);
            b.Initialize(packet_id, 0x9e3779b9u, identity_k, 18u);
            a.PeelCount = b.PeelCount = 44u;
            wirehair::PeelRowIterator ia(a, identity_k, source_prime);
            wirehair::PeelRowIterator ib(b, identity_k, source_prime);
            std::vector<uint16_t> columns;
            do {
                if (ia.GetColumn() != ib.GetColumn() ||
                    ia.GetColumn() >= identity_k)
                {
                    fprintf(stderr,
                        "self-test: production fixed44 iterator compatibility "
                        "broke at packet %u\n", packet_id);
                    return false;
                }
                columns.push_back(ia.GetColumn());
                ++marginal_hits[ia.GetColumn()];
                const bool more_a = ia.Iterate();
                const bool more_b = ib.Iterate();
                if (more_a != more_b)
                {
                    fprintf(stderr,
                        "self-test: production iterator length mismatch\n");
                    return false;
                }
                if (!more_a) {
                    break;
                }
            } while (true);
            std::sort(columns.begin(), columns.end());
            if (columns.size() != 44u ||
                std::adjacent_find(columns.begin(), columns.end()) !=
                    columns.end())
            {
                fprintf(stderr,
                    "self-test: production fixed44 iterator repeated a column\n");
                return false;
            }
        }
        const double expected_hits = 20000.0 * 44.0 / identity_k;
        double marginal_chi2 = 0.0;
        for (uint64_t hits : marginal_hits)
        {
            const double delta = (double)hits - expected_hits;
            marginal_chi2 += delta * delta / expected_hits;
        }
        // df=999: this generous cutoff is over 20 standard deviations above
        // the uniform expectation, so it catches column-selection bias without
        // making the deterministic self-test sensitive to normal variance.
        if (marginal_chi2 > 2000.0)
        {
            fprintf(stderr,
                "self-test: production fixed44 marginal chi2 %.1f is biased\n",
                marginal_chi2);
            return false;
        }
    }

    // Fixed and mixture degree samplers must stay within their declared
    // supports.  Mixture mass checks are intentionally loose enough to avoid
    // flaky sampling failures while detecting parser/profile wiring errors.
    {
        struct DegreeCase
        {
            RowDist Dist;
            unsigned Low;
            unsigned High;
            double HighMass;
        };
        const DegreeCase cases[] = {
            { RowDist::Fixed44, 44, 44, 0.0 },
            { RowDist::Fixed48, 48, 48, 0.0 },
            { RowDist::Mix36_80P10, 36, 80, 0.10 },
            { RowDist::Mix36_80P25, 36, 80, 0.25 },
            { RowDist::Mix36_128P10, 36, 128, 0.10 },
            { RowDist::Mix8_64P50, 8, 64, 0.50 },
            { RowDist::Mix8_48P75, 8, 48, 0.75 },
            { RowDist::Mix12_64P50, 12, 64, 0.50 },
            { RowDist::Mix16_64P50, 16, 64, 0.50 },
            { RowDist::Mix16_96P25, 16, 96, 0.25 },
        };
        const unsigned samples = 100000;
        for (const DegreeCase& dc : cases)
        {
            const DegreeSampler sampler(dc.Dist, 3200);
            Rng rng(UINT64_C(0xd36e5eed) + (unsigned)dc.Dist);
            unsigned high_count = 0;
            for (unsigned i = 0; i < samples; ++i)
            {
                const unsigned degree = sampler.Sample(rng);
                if (degree != dc.Low && degree != dc.High)
                {
                    fprintf(stderr,
                        "self-test: fixed/mixture sampler emitted degree %u\n",
                        degree);
                    return false;
                }
                high_count += degree == dc.High && dc.High != dc.Low;
            }
            if (dc.High != dc.Low)
            {
                const double observed = (double)high_count / samples;
                if (std::fabs(observed - dc.HighMass) > 0.01)
                {
                    fprintf(stderr,
                        "self-test: mixture mass %.4f differs from %.4f\n",
                        observed, dc.HighMass);
                    return false;
                }
            }
        }
    }

    // _s2 schemes: the D2 dense rows must follow the documented Shuffle-2
    // rule exactly -- first row has ceil(span/2) columns; every subsequent
    // row differs from its predecessor in EXACTLY 2 columns.
    {
        auto symdiff = [](const std::vector<uint32_t>& a,
                          const std::vector<uint32_t>& b) -> size_t {
            size_t i = 0, j = 0, diff = 0;
            while (i < a.size() && j < b.size())
            {
                if (a[i] == b[j]) {
                    ++i;
                    ++j;
                }
                else if (a[i] < b[j]) {
                    ++diff;
                    ++i;
                }
                else {
                    ++diff;
                    ++j;
                }
            }
            return diff + (a.size() - i) + (b.size() - j);
        };
        const char* s2_tokens[] = {
            "ldpcdense_s5_d4_s2", "ldpcdense_n12_s6_d5_s2",
            "ldpcdense_s7_d1_s2", "ldpcdense_s6_d2_s2",
            "ldpcdense_n14_s8_d12_s2"
        };
        for (const char* token : s2_tokens)
        {
            for (unsigned k2 = 16; k2 <= 48; k2 += 16)
            {
                Scheme scheme;
                if (!MakeScheme(token, k2, scheme) ||
                    scheme.Kind != SchemeKind::LdpcDenseShuffle2)
                {
                    fprintf(stderr, "self-test: bad _s2 scheme %s\n", token);
                    return false;
                }
                const uint64_t seed =
                    Mix64(HashString64(token) ^ (uint64_t)(k2 * 7919u));
                const GeneratedSystem sys = GenerateSystem(
                    scheme, k2, k2 + 2u, RowDist::LtM1C64, 3, seed);
                const unsigned S = scheme.D;
                const unsigned D2 = scheme.Dense2;
                const unsigned span = k2 + S + D2;
                for (unsigned r = 0; r < S + D2; ++r)
                {
                    if (!sys.Rows[r].IsConstraint)
                    {
                        fprintf(stderr,
                            "self-test: %s constraint layout broke\n", token);
                        return false;
                    }
                }
                if ((unsigned)sys.Rows[S].Columns.size() != (span + 1u) >> 1)
                {
                    fprintf(stderr,
                        "self-test: %s K=%u first _s2 row weight %zu != "
                        "ceil(span/2) = %u\n",
                        token, k2, sys.Rows[S].Columns.size(),
                        (span + 1u) >> 1);
                    return false;
                }
                for (unsigned r = 1; r < D2; ++r)
                {
                    const size_t diff = symdiff(
                        sys.Rows[S + r - 1u].Columns,
                        sys.Rows[S + r].Columns);
                    if (diff != 2u)
                    {
                        fprintf(stderr,
                            "self-test: %s K=%u _s2 row %u differs from "
                            "previous in %zu columns (want 2)\n",
                            token, k2, r, diff);
                        return false;
                    }
                }
            }
        }
    }

    // codecport default follows the codec's K-banded N1 verdict, while
    // explicit n1 suffixes remain fixed controls.
    {
        struct CodecPortCase
        {
            const char* Token;
            unsigned K;
            unsigned N1;
            bool IdentCorner;
        };
        const CodecPortCase cases[] = {
            { "codecport", 9999, 2, false },
            { "codecport", 10000, 3, false },
            { "codecport_n12", 10000, 2, false },
            { "codecport_n13_ic", 10000, 3, true },
        };
        for (const CodecPortCase& c : cases)
        {
            Scheme scheme;
            if (!MakeScheme(c.Token, c.K, scheme) ||
                scheme.Kind != SchemeKind::CodecPort ||
                scheme.N1 != c.N1 ||
                scheme.IdentCorner != c.IdentCorner)
            {
                fprintf(stderr,
                    "self-test: codecport default/override broke for %s K=%u\n",
                    c.Token, c.K);
                return false;
            }
        }
    }

    // codecport identity-corner value generation uses a K+S deck and an
    // owned dense output column, so its first dense row should not be
    // charged for the excluded D2 dense-column half.
    {
        const unsigned cost_k = 128;
        const unsigned S = GetDenseCount((uint32_t)cost_k);
        const unsigned D2 = 12;
        Scheme non_ic, ic;
        if (!MakeScheme("codecport_n13", cost_k, non_ic) ||
            !MakeScheme("codecport_n13_ic", cost_k, ic))
        {
            fprintf(stderr, "self-test: bad codecport cost schemes\n");
            return false;
        }
        const GeneratedSystem a = GenerateSystem(
            non_ic, cost_k, 0, RowDist::LtM1C64, 3,
            UINT64_C(0x1dca7c05));
        const GeneratedSystem b = GenerateSystem(
            ic, cost_k, 0, RowDist::LtM1C64, 3,
            UINT64_C(0x1dca7c05));
        const uint64_t staircase = 3ull * cost_k + S - 1u;
        const uint64_t non_ic_model =
            staircase + (((cost_k + S + D2) + 1u) >> 1) +
            2ull * (D2 - 1u);
        const uint64_t ic_model =
            staircase + (((cost_k + S) + 1u) >> 1) +
            2ull * (D2 - 1u);
        if (a.PrecodeGenXors != non_ic_model ||
            b.PrecodeGenXors != ic_model)
        {
            fprintf(stderr,
                "self-test: codecport generation cost model broke "
                "(non-ic=%llu/%llu ic=%llu/%llu)\n",
                (unsigned long long)a.PrecodeGenXors,
                (unsigned long long)non_ic_model,
                (unsigned long long)b.PrecodeGenXors,
                (unsigned long long)ic_model);
            return false;
        }
    }

    // n1 variants: each SOURCE column must hit exactly X staircase parities
    // (tokens use S > X so distinct parities are always achievable); tokens
    // without _n1 must keep the historical X = 3.
    {
        struct N1Case
        {
            const char* Token;
            unsigned X;
        };
        const N1Case n1_cases[] = {
            { "ldpcdense_n12_s8_d3", 2 },
            { "ldpcdense_n13_s8_d3", 3 },
            { "ldpcdense_n14_s8_d3", 4 },
            { "ldpcdense_s8_d3", 3 },        // n1 absent = 3 (compatibility)
            { "ldpcdense_n12_s8_d3_s2", 2 }, // _s2 must not disturb the hits
            { "ldpc_s8", 3 },
        };
        const unsigned k3 = 24;
        for (const N1Case& nc : n1_cases)
        {
            Scheme scheme;
            if (!MakeScheme(nc.Token, k3, scheme))
            {
                fprintf(stderr, "self-test: bad n1 scheme %s\n", nc.Token);
                return false;
            }
            const GeneratedSystem sys = GenerateSystem(
                scheme, k3, k3, RowDist::LtM1C64, 3,
                Mix64(HashString64(nc.Token) ^ UINT64_C(0x71a7)));
            std::vector<unsigned> hits(k3, 0);
            for (unsigned r = 0; r < scheme.D; ++r)
            {
                for (uint32_t col : sys.Rows[r].Columns)
                {
                    if (col < k3) {
                        ++hits[col];
                    }
                }
            }
            for (unsigned col = 0; col < k3; ++col)
            {
                if (hits[col] != nc.X)
                {
                    fprintf(stderr,
                        "self-test: %s source column %u has %u parity hits "
                        "(want %u)\n",
                        nc.Token, col, hits[col], nc.X);
                    return false;
                }
            }
        }
    }

    // Degree-law chi-square: for every non-Wirehair rowdist, ~200k sampled
    // degrees at K=3200 must match the analytic law transcribed independently
    // in SelfTestDegreeWeight (guards against a silently mis-ported weight
    // formula).  Bins are pooled to expected count >= 10; the bound
    // df + 8*sqrt(2*df) + 30 is loose against sampling noise (typical chi2
    // is df +- 2*sqrt(2*df)) yet far below the chi2 a wrong spike/tau/cap
    // produces (orders of magnitude larger).
    {
        struct DistCase
        {
            RowDist Dist;
            const char* Name;
            unsigned MinD;
            unsigned Cap;
        };
        const DistCase dist_cases[] = {
            { RowDist::LtM1C64, "lt_m1_c64", 1, 64 },
            { RowDist::LtM1C16, "lt_m1_c16", 1, 16 },
            { RowDist::LtM2C1024, "lt_m2_c1024", 2, 1024 },
            { RowDist::RsC001D50C128, "rs_c001_d50_c128", 1, 128 },
        };
        const unsigned chi_k = 3200;
        const unsigned samples = 200000;
        for (const DistCase& dc : dist_cases)
        {
            const DegreeSampler sampler(dc.Dist, chi_k);
            Rng rng(Mix64(HashString64(dc.Name) ^ UINT64_C(0xc41a55a)));
            std::vector<uint64_t> hist(dc.Cap + 1u, 0);
            for (unsigned s = 0; s < samples; ++s)
            {
                const unsigned deg = sampler.Sample(rng);
                if (deg < dc.MinD || deg > dc.Cap)
                {
                    fprintf(stderr,
                        "self-test: %s degree %u outside [%u, %u]\n",
                        dc.Name, deg, dc.MinD, dc.Cap);
                    return false;
                }
                ++hist[deg];
            }
            double total = 0.0;
            std::vector<double> weight(dc.Cap + 1u, 0.0);
            for (unsigned deg = dc.MinD; deg <= dc.Cap; ++deg)
            {
                weight[deg] = SelfTestDegreeWeight(dc.Dist, deg, chi_k);
                total += weight[deg];
            }
            // Pool ascending-degree bins to expected >= 10; the leftover
            // tail merges into the last pooled bin
            std::vector<double> exp_bins, obs_bins;
            double e = 0.0, o = 0.0;
            for (unsigned deg = dc.MinD; deg <= dc.Cap; ++deg)
            {
                e += (double)samples * weight[deg] / total;
                o += (double)hist[deg];
                if (e >= 10.0)
                {
                    exp_bins.push_back(e);
                    obs_bins.push_back(o);
                    e = 0.0;
                    o = 0.0;
                }
            }
            if ((e > 0.0 || o > 0.0) && !exp_bins.empty())
            {
                exp_bins.back() += e;
                obs_bins.back() += o;
            }
            if (exp_bins.size() < 2u)
            {
                fprintf(stderr,
                    "self-test: %s chi-square degenerate\n", dc.Name);
                return false;
            }
            double chi2 = 0.0;
            for (size_t b = 0; b < exp_bins.size(); ++b)
            {
                const double diff = obs_bins[b] - exp_bins[b];
                chi2 += diff * diff / exp_bins[b];
            }
            const double df = (double)(exp_bins.size() - 1u);
            const double bound = df + 8.0 * std::sqrt(2.0 * df) + 30.0;
            if (chi2 > bound)
            {
                fprintf(stderr,
                    "self-test: %s degree chi-square %.1f > bound %.1f "
                    "(df=%.0f): weight-law port mismatch?\n",
                    dc.Name, chi2, bound, df);
                return false;
            }
        }
    }

    // --max-inact runaway guard: capping below the trial's natural
    // inactivation count must mark the trial as a runaway failure; capping
    // at exactly the natural count must reproduce the unlimited result.
    {
        Scheme scheme;
        if (!MakeScheme("dense", 64, scheme))
        {
            fprintf(stderr, "self-test: max-inact scheme broke\n");
            return false;
        }
        const TrialResult full = RunTrial(
            scheme, 64, 0, RowDist::Wirehair, 3, 24681357);
        if (full.Runaway || full.Inactivated < 2u)
        {
            fprintf(stderr,
                "self-test: max-inact reference trial unusable "
                "(runaway=%d inact=%u)\n",
                (int)full.Runaway, full.Inactivated);
            return false;
        }
        const TrialResult capped = RunTrial(
            scheme, 64, 0, RowDist::Wirehair, 3, 24681357,
            nullptr, false, false, 0, full.Inactivated - 1u);
        if (!capped.Runaway || capped.Success || capped.SuccessNoHeavy)
        {
            fprintf(stderr, "self-test: max-inact cap did not abort\n");
            return false;
        }
        const TrialResult exact = RunTrial(
            scheme, 64, 0, RowDist::Wirehair, 3, 24681357,
            nullptr, false, false, 0, full.Inactivated);
        if (exact.Runaway || exact.Def != full.Def ||
            exact.Inactivated != full.Inactivated ||
            exact.ResidualRank != full.ResidualRank)
        {
            fprintf(stderr,
                "self-test: max-inact == natural inact changed the trial\n");
            return false;
        }
    }

    printf("self-test passed (%u rank-invariance instances)\n", checked);
    return true;
}

//------------------------------------------------------------------------------
// CLI

static std::vector<std::string> SplitCsv(const std::string& s)
{
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= s.size())
    {
        size_t comma = s.find(',', start);
        if (comma == std::string::npos) {
            comma = s.size();
        }
        out.push_back(s.substr(start, comma - start));
        start = comma + 1u;
    }
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    std::vector<unsigned> k_list;
    std::vector<unsigned> oh_list;
    std::vector<std::string> scheme_tokens;
    unsigned trials = 100;
    unsigned threads = std::thread::hardware_concurrency();
    uint64_t base_seed = UINT64_C(0x5eed5eed5eed5eed);
    RowDist dist = RowDist::Wirehair;
    unsigned mix = 3;
    unsigned max_inact = 0; // 0 = unlimited
    double max_row_seconds = 0.0; // 0 = unlimited
    unsigned ge_pivot_window = 0; // 0 = first eligible pivot
    bool self_test = false;
    bool ge_replay = false;
    bool ge_replay_reverse = false;
    bool paired = false;
    PacketSchedule packet_schedule = PacketSchedule::None;
    double loss = 0.0;
    bool identity_systematic = false;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                fprintf(stderr, "missing value for %s\n", arg.c_str());
                exit(1);
            }
            return argv[++i];
        };
        if (arg == "--K") {
            for (const std::string& t : SplitCsv(next())) {
                unsigned v = 0;
                if (!ParseUnsigned(t, v)) {
                    fprintf(stderr, "bad --K token %s\n", t.c_str());
                    return 1;
                }
                if (v == 0u) {
                    fprintf(stderr, "bad --K token %s\n", t.c_str());
                    return 1;
                }
                k_list.push_back(v);
            }
        }
        else if (arg == "--oh") {
            for (const std::string& t : SplitCsv(next())) {
                unsigned v = 0;
                if (!ParseUnsigned(t, v)) {
                    fprintf(stderr, "bad --oh token %s\n", t.c_str());
                    return 1;
                }
                oh_list.push_back(v);
            }
        }
        else if (arg == "--schemes") {
            scheme_tokens = SplitCsv(next());
        }
        else if (arg == "--trials") {
            const std::string v = next();
            if (!ParseUnsigned(v, trials)) {
                fprintf(stderr, "bad --trials value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--threads") {
            const std::string v = next();
            if (!ParseUnsigned(v, threads)) {
                fprintf(stderr, "bad --threads value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--seed") {
            const std::string v = next();
            if (!ParseU64(v, base_seed)) {
                fprintf(stderr, "bad --seed value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--mix") {
            const std::string v = next();
            if (!ParseUnsigned(v, mix)) {
                fprintf(stderr, "bad --mix value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--rowdist") {
            const std::string v = next();
            if (v == "wirehair") {
                dist = RowDist::Wirehair;
            }
            else if (v == "lt_m1_c64") {
                dist = RowDist::LtM1C64;
            }
            else if (v == "lt_m1_c16") {
                dist = RowDist::LtM1C16;
            }
            else if (v == "lt_m2_c1024") {
                dist = RowDist::LtM2C1024;
            }
            else if (v == "rs_c001_d50_c128") {
                dist = RowDist::RsC001D50C128;
            }
            else if (v == "fixed32") { dist = RowDist::Fixed32; }
            else if (v == "fixed36") { dist = RowDist::Fixed36; }
            else if (v == "fixed40") { dist = RowDist::Fixed40; }
            else if (v == "fixed44") { dist = RowDist::Fixed44; }
            else if (v == "fixed48") { dist = RowDist::Fixed48; }
            else if (v == "fixed64") { dist = RowDist::Fixed64; }
            else if (v == "fixed80") { dist = RowDist::Fixed80; }
            else if (v == "fixed96") { dist = RowDist::Fixed96; }
            else if (v == "fixed128") { dist = RowDist::Fixed128; }
            else if (v == "fixed256") { dist = RowDist::Fixed256; }
            else if (v == "mix32_128") { dist = RowDist::Mix32_128; }
            else if (v == "mix36_80_p10") {
                dist = RowDist::Mix36_80P10;
            }
            else if (v == "mix36_80_p25") {
                dist = RowDist::Mix36_80P25;
            }
            else if (v == "mix36_128_p10") {
                dist = RowDist::Mix36_128P10;
            }
            else if (v == "mix8_64_p50") {
                dist = RowDist::Mix8_64P50;
            }
            else if (v == "mix8_48_p75") {
                dist = RowDist::Mix8_48P75;
            }
            else if (v == "mix12_64_p50") {
                dist = RowDist::Mix12_64P50;
            }
            else if (v == "mix16_64_p50") {
                dist = RowDist::Mix16_64P50;
            }
            else if (v == "mix16_96_p25") {
                dist = RowDist::Mix16_96P25;
            }
            else if (v == "staged44_16_p125") {
                dist = RowDist::Staged44_16P125;
            }
            else if (v == "staged44_16_p167") {
                dist = RowDist::Staged44_16P167;
            }
            else if (v == "staged44_16_p25") {
                dist = RowDist::Staged44_16P25;
            }
            else if (v == "staged44_12_p167") {
                dist = RowDist::Staged44_12P167;
            }
            else if (v == "log6") { dist = RowDist::Log6; }
            else if (v == "log8") { dist = RowDist::Log8; }
            else {
                fprintf(stderr, "unknown rowdist %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--max-inact") {
            const std::string v = next();
            if (!ParseUnsigned(v, max_inact)) {
                fprintf(stderr, "bad --max-inact value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--max-row-seconds") {
            const std::string v = next();
            if (!ParseNonNegativeDouble(v, max_row_seconds)) {
                fprintf(stderr, "bad --max-row-seconds value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--ge-pivot-window") {
            ge_replay = true; // pivot-window selection is a replay mode
            const std::string v = next();
            if (!ParseUnsigned(v, ge_pivot_window)) {
                fprintf(stderr, "bad --ge-pivot-window value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--ge-replay") {
            ge_replay = true;
        }
        else if (arg == "--ge-replay-reverse") {
            ge_replay = true; // reverse implies replay
            ge_replay_reverse = true;
        }
        else if (arg == "--paired") {
            paired = true;
        }
        else if (arg == "--packet-schedule")
        {
            const std::string v = next();
            if (v == "iid") { packet_schedule = PacketSchedule::Iid; }
            else if (v == "burst") {
                packet_schedule = PacketSchedule::Burst;
            }
            else if (v == "permutation") {
                packet_schedule = PacketSchedule::Permutation;
            }
            else if (v == "repair-only") {
                packet_schedule = PacketSchedule::RepairOnly;
            }
            else {
                fprintf(stderr, "unknown packet schedule %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--loss")
        {
            const std::string v = next();
            if (!ParseNonNegativeDouble(v, loss) ||
                loss > kMaximumSupportedLoss) {
                fprintf(stderr, "bad --loss value %s\n", v.c_str());
                return 1;
            }
        }
        else if (arg == "--identity-systematic") {
            identity_systematic = true;
        }
        else if (arg == "--self-test") {
            self_test = true;
        }
        else {
            fprintf(stderr, "unknown arg %s\n", arg.c_str());
            return 1;
        }
    }

    if (self_test) {
        return SelfTest() ? 0 : 1;
    }

    if (k_list.empty()) {
        k_list.push_back(1000);
    }
    if (oh_list.empty()) {
        oh_list = {0, 1, 2};
    }
    if (scheme_tokens.empty()) {
        scheme_tokens = {"none", "dense", "ldpc"};
    }
    if (threads < 1u) {
        threads = 1;
    }
    if (trials == 0u) {
        fprintf(stderr, "--trials must be positive\n");
        return 1;
    }
    if (packet_schedule == PacketSchedule::None && loss != 0.0)
    {
        fprintf(stderr,
            "--loss requires --packet-schedule so it is not silently ignored\n");
        return 1;
    }
    for (unsigned K : k_list)
    {
        if (K < CAT_WIREHAIR_MIN_N || K > CAT_WIREHAIR_MAX_N)
        {
            fprintf(stderr,
                "--K value %u is outside Wirehair production range [%u, %u]\n",
                K,
                (unsigned)CAT_WIREHAIR_MIN_N,
                (unsigned)CAT_WIREHAIR_MAX_N);
            return 1;
        }
    }

    // Validate every model shape before printing the CSV header or starting
    // workers.  Result pipelines can therefore treat any stdout header as a
    // promise that the complete requested experiment was syntactically valid.
    for (unsigned K : k_list)
    {
        for (const std::string& token : scheme_tokens)
        {
            Scheme scheme;
            if (!MakeScheme(token, K, scheme))
            {
                fprintf(stderr, "bad or out-of-range scheme token %s for K=%u\n",
                    token.c_str(), K);
                return 1;
            }
        }
        for (unsigned oh : oh_list)
        {
            if ((uint64_t)K + oh > UINT16_MAX)
            {
                fprintf(stderr,
                    "bad --oh value %u for K=%u: received-row count exceeds %u\n",
                    oh, K, (unsigned)UINT16_MAX);
                return 1;
            }
        }
    }

    printf("K,scheme,D,H,oh,trials,fail_rate,fail_rate_noheavy,def_mu,def_max,"
           "def_pdf,inact_mu,inact_sd,inact_max,rank_mu,peeled_mu,"
           "residual_rows_mu,recv_xors_per_packet,precode_gen_xors_mu,"
           "sparse_solve_xors_mu,backsub_xors_mu,ge_block_xors_mu,"
           "ge_bitops_mu,heavy_muladds_mu,heavy_divs_mu");
    if (ge_replay) {
        printf(",ge_real_bitops_mu,ge_real_rowops_mu,fill_in_mu,"
               "def_outside_w18_rate,def_band_w95,def_band_w99");
    }
    // Preserve every historical column position and append run metadata so
    // positional consumers continue to work while new CSVs self-describe.
    printf(",runaway_rate,packet_schedule_exhausted_rate,rowdist,"
           "packet_schedule,loss,identity_systematic,mix,base_seed,paired,"
           "max_inact,max_row_seconds,requested_trials,threads,ge_replay,"
           "ge_replay_reverse,ge_pivot_window");
    printf("\n");
    fflush(stdout);

    for (unsigned K : k_list)
    {
        for (const std::string& token : scheme_tokens)
        {
            Scheme scheme;
            if (!MakeScheme(token, K, scheme))
            {
                fprintf(stderr, "bad scheme token %s\n", token.c_str());
                return 1;
            }
            for (unsigned oh : oh_list)
            {
                Aggregate agg;
                std::mutex agg_mutex;
                std::atomic<unsigned> next_trial(0);
                std::atomic<bool> stop_row(false);
                const auto row_start = std::chrono::steady_clock::now();

                auto worker = [&]() {
                    for (;;)
                    {
                        if (stop_row.load(std::memory_order_relaxed)) {
                            return;
                        }
                        const unsigned trial = next_trial.fetch_add(1u);
                        if (trial >= trials) {
                            return;
                        }
                        if (max_row_seconds > 0.0)
                        {
                            // Always allow the first trial to start,
                            // so tiny bounds do not emit zero-trial rows.
                            const double elapsed =
                                std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() -
                                    row_start).count();
                            if (trial > 0u && elapsed >= max_row_seconds)
                            {
                                stop_row.store(
                                    true, std::memory_order_relaxed);
                                return;
                            }
                        }
                        const uint64_t seed = Mix64(
                            base_seed ^
                            Mix64((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)) ^
                            Mix64(HashString64(token)) ^
                            Mix64((uint64_t)trial * UINT64_C(0xbf58476d1ce4e5b9)));
                        // --paired: received-row seed omits both scheme and
                        // overhead.  K+OH rows are therefore a prefix of the
                        // same trial trace at every overhead point, while the
                        // system seed above keeps the scheme-specific precode
                        // fixed throughout that sweep.
                        const uint64_t recv_seed = Mix64(
                            base_seed ^
                            Mix64((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)) ^
                            Mix64((uint64_t)trial * UINT64_C(0xbf58476d1ce4e5b9)));
                        const TrialResult t = RunTrial(
                            scheme, K, oh, dist, mix, seed,
                            paired ? &recv_seed : nullptr,
                            ge_replay, ge_replay_reverse, ge_pivot_window,
                            max_inact, packet_schedule, loss,
                            identity_systematic);
                        std::lock_guard<std::mutex> lock(agg_mutex);
                        agg.Add(t);
                    }
                };

                std::vector<std::thread> pool;
                const unsigned use_threads =
                    threads < trials ? threads : (trials > 0 ? trials : 1u);
                for (unsigned w = 0; w < use_threads; ++w) {
                    pool.emplace_back(worker);
                }
                for (std::thread& th : pool) {
                    th.join();
                }
                if (max_row_seconds > 0.0 && agg.Trials < trials)
                {
                    fprintf(stderr,
                        "# max-row-seconds stopped K=%u scheme=%s oh=%u "
                        "after %llu/%u trials\n",
                        K,
                        token.c_str(),
                        oh,
                        (unsigned long long)agg.Trials,
                        trials);
                }

                const double n = (double)(agg.Trials ? agg.Trials : 1u);
                // Runaway trials (--max-inact) are in fail_rate but have no
                // def/inact/cost statistics: means use completed trials only
                const uint64_t completed = agg.Trials - agg.Runaways -
                    agg.ScheduleExhaustions;
                const double m = (double)(completed ? completed : 1u);
                const double inact_mu = agg.InactSum / m;
                const double inact_var =
                    agg.InactSqSum / m - inact_mu * inact_mu;
                // def_pdf keeps its contract (re-scoring it at any H must
                // reproduce the fail rate): it is normalized over ALL trials
                // and runaways appear as a sentinel def=999999 bucket that
                // fails at every H.  With no runaways this is the historical
                // per-completed-trial pdf unchanged.
                std::string def_pdf;
                for (size_t d = 0; d < agg.DefHist.size(); ++d)
                {
                    if (agg.DefHist[d] == 0u) {
                        continue;
                    }
                    char buf[64];
                    snprintf(buf, sizeof(buf), "%s%zu:%.6f",
                        def_pdf.empty() ? "" : "|",
                        d, agg.DefHist[d] / n);
                    def_pdf += buf;
                }
                if (agg.Runaways > 0u)
                {
                    char buf[64];
                    snprintf(buf, sizeof(buf), "%s999999:%.6f",
                        def_pdf.empty() ? "" : "|",
                        agg.Runaways / n);
                    def_pdf += buf;
                }
                if (agg.ScheduleExhaustions > 0u)
                {
                    char buf[64];
                    snprintf(buf, sizeof(buf), "%s999998:%.6f",
                        def_pdf.empty() ? "" : "|",
                        agg.ScheduleExhaustions / n);
                    def_pdf += buf;
                }
                printf("%u,%s,%u,%u,%u,%llu,%.6f,%.6f,%.4f,%.0f,%s,"
                       "%.2f,%.2f,%.0f,%.2f,%.2f,%.2f,%.4f,%.1f,%.1f,%.1f,"
                       "%.1f,%.1f,%.1f,%.1f",
                    K, token.c_str(), scheme.D, scheme.H, oh,
                    (unsigned long long)agg.Trials,
                    agg.Fails / n,
                    agg.FailsNoHeavy / n,
                    agg.DefSum / m,
                    agg.DefMax,
                    def_pdf.empty() ? "0:1.000000" : def_pdf.c_str(),
                    inact_mu,
                    inact_var > 0.0 ? std::sqrt(inact_var) : 0.0,
                    agg.InactMax,
                    agg.RankSum / m,
                    agg.PeeledSum / m,
                    agg.ResidualRowsSum / m,
                    agg.RecvXorSum / m / (double)(K + oh),
                    agg.PrecodeXorSum / m,
                    agg.SparseXorSum / m,
                    agg.BacksubXorSum / m,
                    agg.GeBlockXorSum / m,
                    agg.GeBitOpSum / m,
                    agg.HeavyMuladdSum / m,
                    agg.HeavyDivSum / m);
                if (ge_replay)
                {
                    std::sort(agg.DefBands.begin(), agg.DefBands.end());
                    printf(",%.1f,%.1f,%.2f,%.6f,%u,%u",
                        agg.GeRealWordXorSum / m,
                        agg.GeRealRowOpSum / m,
                        agg.FillInSum / m,
                        (double)agg.DefOutsideW18Count / m,
                        PercentileSorted(agg.DefBands, 0.95),
                        PercentileSorted(agg.DefBands, 0.99));
                }
                printf(",%.6f,%.6f,%s,%s,%.17g,%u,%u,%llu,%u,%u,%.17g,%u,%u,"
                       "%u,%u,%u",
                    agg.Runaways / n,
                    agg.ScheduleExhaustions / n,
                    RowDistName(dist), PacketScheduleName(packet_schedule),
                    loss, identity_systematic ? 1u : 0u, mix,
                    (unsigned long long)base_seed, paired ? 1u : 0u,
                    max_inact, max_row_seconds, trials, threads,
                    ge_replay ? 1u : 0u, ge_replay_reverse ? 1u : 0u,
                    ge_pivot_window);
                printf("\n");
                fflush(stdout);
            }
        }
    }
    return 0;
}
