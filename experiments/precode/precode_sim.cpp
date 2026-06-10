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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
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
        const uint32_t threshold = (uint32_t)((0u - (uint64_t)bound) % bound);
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
    Wirehair,   // production GeneratePeelRowWeight
    LtM1C64,    // soliton-ish truncated LT, min degree 1, cap 64
    LtM2C1024   // min degree 2, cap 1024 (residual-quality frontier family)
};

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
        unsigned min_d = 1, cap = 64;
        if (dist == RowDist::LtM2C1024) {
            min_d = 2;
            cap = 1024;
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
    LdpcDense    // S staircase parity columns + D2 dense p=0.5 rows over all
                 // K + S + D2 binary columns: staircase peels cheaply, dense
                 // rows pin the residual tail
};

struct Scheme
{
    std::string Name;
    SchemeKind Kind;
    unsigned D;       // binary precode rows/columns (staircase S for LdpcDense)
    unsigned H;       // heavy GF(256) rows/columns (MDS patch model)
    unsigned Weight;  // DenseSparse row weight
    unsigned Dense2;  // LdpcDense extra dense row/column count
};

//------------------------------------------------------------------------------
// Trial result

struct TrialResult
{
    bool Ok;                 // generation/solve consistency
    bool SuccessNoHeavy;     // def == 0
    bool Success;            // def <= H
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

struct GeneratedSystem
{
    unsigned L;
    unsigned K;
    unsigned D;
    unsigned H;
    std::vector<SparseRow> Rows;     // binary rows only (constraints first)
    uint64_t RecvRowGenXors;
    uint64_t PrecodeGenXors;
};

static GeneratedSystem GenerateSystem(
    const Scheme& scheme,
    unsigned K,
    unsigned received,
    RowDist dist,
    unsigned mix,
    uint64_t seed)
{
    GeneratedSystem sys;
    sys.K = K;
    sys.D = scheme.D + scheme.Dense2;
    sys.H = scheme.H;
    sys.L = K + sys.D + scheme.H;
    sys.RecvRowGenXors = 0;
    sys.PrecodeGenXors = 0;

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
    {
        const unsigned S = scheme.D;
        std::vector<std::vector<uint32_t> > checks(S);
        for (unsigned c = 0; c < K; ++c)
        {
            if (scheme.Kind != SchemeKind::LdpcTri)
            {
                // 3 distinct random parities per source column
                uint32_t p1 = rng.Below(S);
                uint32_t p2 = rng.Below(S);
                while (S > 1u && p2 == p1) {
                    p2 = rng.Below(S);
                }
                uint32_t p3 = rng.Below(S);
                while (S > 2u && (p3 == p1 || p3 == p2)) {
                    p3 = rng.Below(S);
                }
                checks[p1].push_back(c);
                if (S > 1u) {
                    checks[p2].push_back(c);
                }
                if (S > 2u) {
                    checks[p3].push_back(c);
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
        break;
    }
    }

    // --- Received rows ---
    const DegreeSampler degrees(dist, K);
    for (unsigned r = 0; r < received; ++r)
    {
        SparseRow row;
        row.IsConstraint = false;
        const unsigned degree = degrees.Sample(rng);
        AddDistinctColumns(row.Columns, 0, K, degree, rng);
        if (precode_space > 0u) {
            AddDistinctColumns(row.Columns, K, precode_space, mix, rng);
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

    PeelOutcome Run()
    {
        PeelOutcome out;
        out.ColumnState.assign(ColumnCount, kStatePeeled);
        out.SolveRow.assign(ColumnCount, UINT32_MAX);
        out.SparseSolveXors = 0;

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
        : Bits(bits), Words((bits + 63u) / 64u) {}

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

static ResidualAnalysis AnalyzeResidual(
    const GeneratedSystem& sys,
    const PeelOutcome& peel)
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
        rank += rank_matrix.Insert(rowbits);
    }

    res.ResidualRows = residual_rows;
    res.Rank = rank;
    return res;
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
    uint64_t seed)
{
    TrialResult t;
    std::memset(&t, 0, sizeof(t));

    // Received budget: the decoder collects K + overhead packets.  The D + H
    // constraint rows are free (known structure), so the binary system is
    // (K + overhead + D) rows over L = K + D + H columns, and L - K of the
    // unknowns are precode blocks the decoder does not return to the caller.
    const GeneratedSystem sys =
        GenerateSystem(scheme, K, K + overhead, dist, mix, seed);

    const PeelOutcome peel = PeelSolver(sys.L, sys.Rows).Run();
    const ResidualAnalysis res = AnalyzeResidual(sys, peel);

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
    return t;
}

//------------------------------------------------------------------------------
// Scheme parsing

static bool ParseUnsigned(const std::string& s, unsigned& out)
{
    if (s.empty()) {
        return false;
    }
    char* end = nullptr;
    const unsigned long v = strtoul(s.c_str(), &end, 10);
    if (!end || *end != '\0') {
        return false;
    }
    out = (unsigned)v;
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
//   heavyonly        (D = 0, H = 16)
// Optional suffix on any token: _h<H> to override heavy row count.
static bool MakeScheme(const std::string& token, unsigned K, Scheme& out)
{
    std::string body = token;
    out.H = 6;
    out.Weight = 0;
    out.Dense2 = 0;

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
    else if (body.rfind("ldpcdense_s", 0) == 0)
    {
        out.Kind = SchemeKind::LdpcDense;
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
    else
    {
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
    std::vector<uint32_t> DefHist; // def occurrence counts

    void Add(const TrialResult& t)
    {
        ++Trials;
        if (!t.Success) {
            ++Fails;
        }
        if (!t.SuccessNoHeavy) {
            ++FailsNoHeavy;
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
    }
};

//------------------------------------------------------------------------------
// Self-test

static bool SelfTest()
{
    // Rank invariance: peeled + residual rank must equal brute-force rank of
    // the whole binary system, for every scheme and instance.
    const char* tokens[] = {
        "none", "dense_d8", "densesparse_w6_d8", "ldpc_s5", "ldpctri_s5",
        "heavyonly_h4", "ldpcdense_s5_d4", "ldpcdense_s6_d3_h8"
    };
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
                    const ResidualAnalysis res = AnalyzeResidual(sys, peel);
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
                    ++checked;
                }
            }
        }
    }

    // Determinism
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
        if (comma > start) {
            out.push_back(s.substr(start, comma - start));
        }
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
    bool self_test = false;

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
                if (ParseUnsigned(t, v)) {
                    k_list.push_back(v);
                }
            }
        }
        else if (arg == "--oh") {
            for (const std::string& t : SplitCsv(next())) {
                unsigned v = 0;
                if (ParseUnsigned(t, v)) {
                    oh_list.push_back(v);
                }
            }
        }
        else if (arg == "--schemes") {
            scheme_tokens = SplitCsv(next());
        }
        else if (arg == "--trials") {
            ParseUnsigned(next(), trials);
        }
        else if (arg == "--threads") {
            ParseUnsigned(next(), threads);
        }
        else if (arg == "--seed") {
            base_seed = strtoull(next().c_str(), nullptr, 0);
        }
        else if (arg == "--mix") {
            ParseUnsigned(next(), mix);
        }
        else if (arg == "--rowdist") {
            const std::string v = next();
            if (v == "wirehair") {
                dist = RowDist::Wirehair;
            }
            else if (v == "lt_m1_c64") {
                dist = RowDist::LtM1C64;
            }
            else if (v == "lt_m2_c1024") {
                dist = RowDist::LtM2C1024;
            }
            else {
                fprintf(stderr, "unknown rowdist %s\n", v.c_str());
                return 1;
            }
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

    printf("K,scheme,D,H,oh,trials,fail_rate,fail_rate_noheavy,def_mu,def_max,"
           "def_pdf,inact_mu,inact_sd,inact_max,rank_mu,peeled_mu,"
           "residual_rows_mu,recv_xors_per_packet,precode_gen_xors_mu,"
           "sparse_solve_xors_mu,backsub_xors_mu,ge_block_xors_mu,"
           "ge_bitops_mu\n");
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

                auto worker = [&]() {
                    for (;;)
                    {
                        const unsigned trial = next_trial.fetch_add(1u);
                        if (trial >= trials) {
                            return;
                        }
                        const uint64_t seed = Mix64(
                            base_seed ^
                            Mix64((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)) ^
                            Mix64(HashString64(token)) ^
                            Mix64((uint64_t)oh * UINT64_C(0xd6e8feb86659fd93)) ^
                            Mix64((uint64_t)trial * UINT64_C(0xbf58476d1ce4e5b9)));
                        const TrialResult t =
                            RunTrial(scheme, K, oh, dist, mix, seed);
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

                const double n = (double)(agg.Trials ? agg.Trials : 1u);
                const double inact_mu = agg.InactSum / n;
                const double inact_var =
                    agg.InactSqSum / n - inact_mu * inact_mu;
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
                printf("%u,%s,%u,%u,%u,%llu,%.6f,%.6f,%.4f,%.0f,%s,"
                       "%.2f,%.2f,%.0f,%.2f,%.2f,%.2f,%.4f,%.1f,%.1f,%.1f,"
                       "%.1f,%.1f\n",
                    K, token.c_str(), scheme.D, scheme.H, oh,
                    (unsigned long long)agg.Trials,
                    agg.Fails / n,
                    agg.FailsNoHeavy / n,
                    agg.DefSum / n,
                    agg.DefMax,
                    def_pdf.empty() ? "0:1.000000" : def_pdf.c_str(),
                    inact_mu,
                    inact_var > 0.0 ? std::sqrt(inact_var) : 0.0,
                    agg.InactMax,
                    agg.RankSum / n,
                    agg.PeeledSum / n,
                    agg.ResidualRowsSum / n,
                    agg.RecvXorSum / n / (double)(K + oh),
                    agg.PrecodeXorSum / n,
                    agg.SparseXorSum / n,
                    agg.BacksubXorSum / n,
                    agg.GeBlockXorSum / n,
                    agg.GeBitOpSum / n);
                fflush(stdout);
            }
        }
    }
    return 0;
}
