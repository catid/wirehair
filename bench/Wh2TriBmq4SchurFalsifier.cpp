#include "Wh2FrozenTrace.h"

#include "../WirehairTools.h"
#include "../codec/WirehairV2PrecodeEncode.h"
#include "../codec/WirehairV2RawSeed.h"
#include "../codec/WirehairV2Seeds.h"
#include "../codec/WirehairV2Solve.h"
#include "../gf256.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#define WIREHAIR_WH2_SOURCE_GIT_COMMIT \
    "0000000000000000000000000000000000000000"
#endif

#if !defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "tri-bmq4-schur-v0 falsifier requires the test-hook policy library"
#endif

#if defined(WH_SEED_KNOBS)
#error "tri-bmq4-schur-v0 falsifier forbids seed-knob overrides"
#endif

namespace {

using wirehair::wh2_benchmark::CanonicalRecoveryCellJson;
using wirehair::wh2_benchmark::CopyNestedPrefix;
using wirehair::wh2_benchmark::EnumerateDevelopmentRecoveryCells;
using wirehair::wh2_benchmark::FrozenPacketTrace;
using wirehair::wh2_benchmark::FrozenRecoveryCell;
using wirehair::wh2_benchmark::FrozenSchedule;
using wirehair::wh2_benchmark::FrozenTraceStatus;
using wirehair::wh2_benchmark::GenerateFrozenPacketTrace;
using wirehair::wh2_benchmark::PacketIdsSha256;
using wirehair::wh2_benchmark::RecoveryCellSha256;
using wirehair::wh2_benchmark::Sha256Hex;

static const uint32_t kK = 64000u;
static const uint32_t kBlockBytes = 2u;
static const uint32_t kLossPpm = 500000u;
static const uint32_t kAttempt = 0u;
static const uint32_t kDenseRows = 12u;
static const uint32_t kHeavyRows = 12u;
static const uint32_t kCoreLimit = 1024u;
static const uint32_t kDeadlineSeconds = 120u;
static const uint32_t kExpectedPrecodeContractVersion = 2u;
static const uint32_t kExpectedPacketContractVersion = 4u;
static const int kExitPass = 0;
static const int kExitInvalid = 1;
static const int kExitReject = 2;
static const int kExitUsage = 64;
static const uint64_t kLossRoot = UINT64_C(0xd1b54a32d192ed03);
static const uint64_t kFreeSalt = UINT64_C(0x465245452d424d34);
static const uint64_t kRowSalt = UINT64_C(0x524f572d424d5134);
static const uint64_t kDegreeSalt = UINT64_C(0x4445472d424d5134);
static const uint64_t kExcessSalt = UINT64_C(0x4558542d424d5134);
static const char kCandidate[] = "tri-bmq4-schur-v0";
static const char kSchema[] =
    "wh2-tri-bmq4-schur-v0-falsifier-receipt-v1";
static const char kBead[] = "wirehair-sxvz.16.1.20.24";
static const char kDesignCommit[] =
    "04d8114e908e2c3038da8e45c4ab24a214db711a";
static const char kExpectedCellSha[] =
    "88873aae0ae86c8e623f9a9a87755990d20ba5cb32e4c1ee20ee5e5a50c97f49";
static const char kExpectedTraceSha[] =
    "7321bfedca2201737d3a2da7dd5fb397d1a10a93f404837ec46a58a6927f236d";
static const char kExpectedPrefixSha[] =
    "7ac4dc7fd943f192030178a7faf7b9acf2deb96efbbc5cdf6780b235733eb551";
static const char kSeedScheduleSha[] =
    "90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7";
static const uint64_t kExpectedTraceSeed = UINT64_C(0x83dae2a65577a471);
static const uint64_t kExpectedBurstSeed = UINT64_C(0x83dae2a655675eaf);
static const uint64_t kExpectedCandidateLimit = UINT64_C(16450560);
static const uint64_t kExpectedAttemptedCandidates = UINT64_C(128308);
static const uint64_t kExpectedPrecodeSeed = UINT64_C(0x487468302aad7105);
static const uint32_t kExpectedPacketSeed = UINT32_C(0x4ec72102);

static_assert(
    wirehair_v2::kPrecodeContractVersion ==
        kExpectedPrecodeContractVersion,
    "tri-bmq4-schur-v0 requires precode contract v2");
static_assert(
    wirehair_v2::kPacketRowContractVersion ==
        kExpectedPacketContractVersion,
    "tri-bmq4-schur-v0 requires packet-row contract v4");

enum class EvaluationOutcome : uint8_t
{
    Pass,
    Reject,
    Invalid,
    Deadline
};

struct EvaluationMapping
{
    const char* Status;
    const char* Stage;
    int ExitCode;
};

EvaluationMapping MapEvaluationOutcome(
    EvaluationOutcome outcome,
    const char* reject_stage)
{
    switch (outcome)
    {
    case EvaluationOutcome::Pass:
        return EvaluationMapping{ "pass", "none", kExitPass };
    case EvaluationOutcome::Reject:
        return EvaluationMapping{ "reject", reject_stage, kExitReject };
    case EvaluationOutcome::Deadline:
        return EvaluationMapping{ "invalid", "deadline", kExitInvalid };
    case EvaluationOutcome::Invalid:
    default:
        return EvaluationMapping{ "invalid", "oracle", kExitInvalid };
    }
}

class Deadline
{
public:
    explicit Deadline(uint32_t seconds)
        : End(std::chrono::steady_clock::now() +
              std::chrono::seconds(seconds))
    {
    }

    bool Expired() const
    {
        return std::chrono::steady_clock::now() >= End;
    }

private:
    std::chrono::steady_clock::time_point End;
};

std::string Hex64(uint64_t value)
{
    std::ostringstream out;
    out << "0x" << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

std::string Hex32(uint32_t value)
{
    std::ostringstream out;
    out << "0x" << std::hex << std::setfill('0') << std::setw(8) << value;
    return out.str();
}

class DigestBuilder
{
public:
    void U8(uint8_t value) { Bytes.push_back(value); }

    void U32(uint32_t value)
    {
        for (unsigned i = 0u; i < 4u; ++i) {
            Bytes.push_back(static_cast<uint8_t>(value >> (i * 8u)));
        }
    }

    void U64(uint64_t value)
    {
        for (unsigned i = 0u; i < 8u; ++i) {
            Bytes.push_back(static_cast<uint8_t>(value >> (i * 8u)));
        }
    }

    void String(const std::string& value)
    {
        U64(static_cast<uint64_t>(value.size()));
        Bytes.insert(Bytes.end(), value.begin(), value.end());
    }

    void Tag(const char* value) { String(value); }

    std::string Finish() const
    {
        return Sha256Hex(
            Bytes.empty() ? NULL : Bytes.data(), Bytes.size());
    }

private:
    std::vector<uint8_t> Bytes;
};

struct RowView
{
    const uint32_t* First;
    const uint32_t* Last;

    const uint32_t* begin() const { return First; }
    const uint32_t* end() const { return Last; }
    size_t size() const { return static_cast<size_t>(Last - First); }
    bool empty() const { return First == Last; }
};

class ExplicitRows
{
public:
    ExplicitRows() { Offsets.push_back(0u); }

    void Reserve(size_t rows, size_t references)
    {
        Offsets.reserve(rows + 1u);
        Columns.reserve(references);
    }

    void Append(const std::vector<uint32_t>& row)
    {
        Columns.insert(Columns.end(), row.begin(), row.end());
        Offsets.push_back(Columns.size());
    }

    RowView operator[](size_t row) const
    {
        static const uint32_t empty_storage = 0u;
        const uint32_t* const base = Columns.empty() ?
            &empty_storage : Columns.data();
        RowView view;
        view.First = base + Offsets[row];
        view.Last = base + Offsets[row + 1u];
        return view;
    }

    size_t size() const { return Offsets.size() - 1u; }
    size_t references() const { return Columns.size(); }

private:
    std::vector<size_t> Offsets;
    std::vector<uint32_t> Columns;
};

bool ValidateRows(
    const ExplicitRows& rows,
    uint32_t column_count,
    std::string& error)
{
    std::vector<uint32_t> last_row(column_count, UINT32_MAX);
    for (uint32_t row = 0u; row < rows.size(); ++row)
    {
        const RowView view = rows[row];
        if (view.empty())
        {
            error = "empty explicit row";
            return false;
        }
        for (uint32_t column : view)
        {
            if (column >= column_count)
            {
                error = "explicit row column out of range";
                return false;
            }
            if (last_row[column] == row)
            {
                error = "duplicate column in explicit row";
                return false;
            }
            last_row[column] = row;
        }
    }
    return true;
}

class MaskBasis12
{
public:
    MaskBasis12() : RankValue(0u) { Basis.fill(0u); }

    bool Insert(uint16_t value)
    {
        uint16_t reduced = value;
        for (uint32_t bit = 0u; bit < kDenseRows; ++bit)
        {
            const uint16_t mask = static_cast<uint16_t>(1u << bit);
            if ((reduced & mask) == 0u) {
                continue;
            }
            if (Basis[bit] != 0u) {
                reduced ^= Basis[bit];
                continue;
            }
            Basis[bit] = reduced;
            ++RankValue;
            return true;
        }
        return false;
    }

    uint32_t Rank() const { return RankValue; }

private:
    std::array<uint16_t, kDenseRows> Basis;
    uint32_t RankValue;
};

bool SameParams(
    const wirehair_v2::PrecodeParams& a,
    const wirehair_v2::PrecodeParams& b)
{
    return a.BlockCount == b.BlockCount &&
        a.Staircase == b.Staircase &&
        a.DenseRows == b.DenseRows &&
        a.HeavyRows == b.HeavyRows &&
        a.SourceHits == b.SourceHits &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.HeavyFamily == b.HeavyFamily &&
        a.DenseAnchors == b.DenseAnchors &&
        a.Seed == b.Seed;
}

bool ValidateRawSeedLaw(std::string& error)
{
    if (std::strcmp(
            wirehair_v2::test::kRawArchitectureSeedBasis,
            "uniform-raw-v1") != 0 ||
        wirehair_v2::test::kRawArchitecturePrecodeSeed !=
            kExpectedPrecodeSeed ||
        wirehair_v2::test::kRawArchitecturePacketSeed !=
            kExpectedPacketSeed ||
        std::strcmp(
            wirehair_v2::test::kRawArchitectureSeedScheduleSha256,
            kSeedScheduleSha) != 0 ||
        Sha256Hex(
            wirehair_v2::test::kRawArchitectureSeedScheduleCanonical) !=
            kSeedScheduleSha)
    {
        error = "uniform-raw-v1 seed basis receipt mismatch";
        return false;
    }
    return true;
}

bool ValidateContractVersions(std::string& error)
{
    if (wirehair_v2::PrecodeContractVersion() !=
            kExpectedPrecodeContractVersion ||
        wirehair_v2::kPacketRowContractVersion !=
            kExpectedPacketContractVersion)
    {
        error = "precode or packet-row contract version mismatch";
        return false;
    }
    return true;
}

uint32_t RankGf256(std::vector<uint8_t>& matrix,
                   uint32_t rows,
                   uint32_t columns,
                   Deadline* deadline,
                   bool& timed_out)
{
    timed_out = false;
    uint32_t rank = 0u;
    for (uint32_t column = 0u; column < columns && rank < rows; ++column)
    {
        if ((column & 15u) == 0u && deadline && deadline->Expired())
        {
            timed_out = true;
            return rank;
        }
        uint32_t pivot = rank;
        while (pivot < rows &&
               matrix[static_cast<size_t>(pivot) * columns + column] == 0u)
        {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        if (pivot != rank)
        {
            for (uint32_t j = 0u; j < columns; ++j) {
                std::swap(matrix[static_cast<size_t>(pivot) * columns + j],
                          matrix[static_cast<size_t>(rank) * columns + j]);
            }
        }
        const uint8_t pivot_value =
            matrix[static_cast<size_t>(rank) * columns + column];
        const uint8_t inverse = gf256_inv(pivot_value);
        const uint8_t* const pivot_row =
            matrix.data() + static_cast<size_t>(rank) * columns;
        for (uint32_t target = rank + 1u; target < rows; ++target)
        {
            uint8_t* const target_row =
                matrix.data() + static_cast<size_t>(target) * columns;
            const uint8_t coefficient = target_row[column];
            if (coefficient == 0u) {
                continue;
            }
            const uint8_t scale = gf256_mul(coefficient, inverse);
            gf256_muladd_mem(
                target_row + column,
                scale,
                pivot_row + column,
                static_cast<int>(columns - column));
        }
        ++rank;
    }
    return rank;
}

bool RankPrecodeRestriction(
    const wirehair_v2::PrecodeSystem& system,
    const std::vector<uint32_t>& selected_columns,
    Deadline* deadline,
    uint32_t& rank,
    std::string& error)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint32_t P = S + D2 + H;
    const uint32_t L = K + P;
    const uint32_t width = static_cast<uint32_t>(selected_columns.size());
    rank = 0u;
    if (selected_columns.size() > UINT32_MAX)
    {
        error = "restriction width overflow";
        return false;
    }
    if (width == 0u) {
        return true;
    }
    std::vector<int32_t> selected_index(L, -1);
    for (uint32_t index = 0u; index < width; ++index)
    {
        const uint32_t column = selected_columns[index];
        if (column >= L || selected_index[column] >= 0)
        {
            error = "invalid restriction column list";
            return false;
        }
        selected_index[column] = static_cast<int32_t>(index);
    }
    if (static_cast<uint64_t>(P) * width >
        static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
    {
        error = "restriction allocation overflow";
        return false;
    }
    std::vector<uint8_t> matrix(static_cast<size_t>(P) * width, 0u);
    uint32_t row_index = 0u;
    const auto add_binary = [&](const std::vector<uint32_t>& row) {
        for (uint32_t column : row)
        {
            const int32_t index = selected_index[column];
            if (index >= 0) {
                matrix[static_cast<size_t>(row_index) * width +
                       static_cast<uint32_t>(index)] ^= 1u;
            }
        }
        ++row_index;
    };
    for (const std::vector<uint32_t>& row : system.StaircaseRows) {
        add_binary(row);
    }
    for (const std::vector<uint32_t>& row : system.DenseBasisRowColumns) {
        add_binary(row);
    }
    for (uint32_t heavy = 0u; heavy < H; ++heavy)
    {
        uint8_t* const row = matrix.data() +
            static_cast<size_t>(row_index++) * width;
        for (uint32_t index = 0u; index < width; ++index) {
            row[index] = wirehair_v2::HeavyCoefficientForParams(
                system.Params, heavy, selected_columns[index]);
        }
    }
    if (row_index != P)
    {
        error = "precode restriction row-count mismatch";
        return false;
    }
    bool timed_out = false;
    rank = RankGf256(matrix, P, width, deadline, timed_out);
    if (timed_out)
    {
        error = "deadline";
        return false;
    }
    return true;
}

struct Partition
{
    std::vector<uint32_t> Qd;
    std::vector<uint16_t> QdMasks;
    std::vector<uint32_t> Q;
    std::vector<uint32_t> F;
    std::vector<uint8_t> IsQ;
    std::vector<uint8_t> IsQd;
    std::vector<uint16_t> DenseIncidence;
    std::vector<uint16_t> SourceHitRows;
    uint32_t QSourceCount = 0u;
    uint32_t SchurRank = 0u;
    uint32_t FullQRank = 0u;
};

EvaluationOutcome BuildPartition(
    const wirehair_v2::PrecodeSystem& system,
    Deadline* deadline,
    Partition& out,
    std::string& error)
{
    out = Partition();
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint32_t N1 = system.Params.SourceHits;
    const uint32_t P = S + D2 + H;
    const uint32_t L = K + P;
    const uint32_t dense_base = K + S;
    const uint32_t heavy_base = dense_base + D2;
    if (K < 2u || S == 0u || S > UINT16_MAX ||
        D2 != kDenseRows || H != kHeavyRows || N1 == 0u || N1 > 8u ||
        system.StaircaseRows.size() != S ||
        system.DenseBasisRowColumns.size() != D2 ||
        static_cast<uint64_t>(K) + P > UINT32_MAX ||
        static_cast<uint64_t>(K) * N1 >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
    {
        error = "partition input shape is invalid";
        return EvaluationOutcome::Invalid;
    }

    out.DenseIncidence.assign(L, 0u);
    for (uint32_t row = 0u; row < D2; ++row)
    {
        for (uint32_t column : system.DenseBasisRowColumns[row])
        {
            if (column >= heavy_base)
            {
                error = "dense basis column outside binary span";
                return EvaluationOutcome::Invalid;
            }
            out.DenseIncidence[column] ^=
                static_cast<uint16_t>(1u << row);
        }
    }

    out.SourceHitRows.assign(
        static_cast<size_t>(K) * N1, UINT16_MAX);
    std::vector<uint8_t> hit_count(K, 0u);
    for (uint32_t row = 0u; row < S; ++row)
    {
        for (uint32_t column : system.StaircaseRows[row])
        {
            if (column >= dense_base)
            {
                error = "staircase column outside staircase span";
                return EvaluationOutcome::Invalid;
            }
            if (column >= K) {
                continue;
            }
            uint8_t& count = hit_count[column];
            if (count >= N1)
            {
                error = "too many staircase hits for source";
                return EvaluationOutcome::Invalid;
            }
            out.SourceHitRows[static_cast<size_t>(column) * N1 + count] =
                static_cast<uint16_t>(row);
            ++count;
        }
    }
    for (uint32_t source = 0u; source < K; ++source)
    {
        if (hit_count[source] != N1)
        {
            error = "source staircase hit-count mismatch";
            return EvaluationOutcome::Invalid;
        }
    }

    std::vector<uint16_t> suffix(S + 1u, 0u);
    for (uint32_t row = S; row-- > 0u;) {
        suffix[row] = static_cast<uint16_t>(
            suffix[row + 1u] ^ out.DenseIncidence[K + row]);
    }
    const auto source_vector = [&](uint32_t source) {
        uint16_t value = out.DenseIncidence[source];
        for (uint32_t hit = 0u; hit < N1; ++hit) {
            value ^= suffix[out.SourceHitRows[
                static_cast<size_t>(source) * N1 + hit]];
        }
        return value;
    };

    MaskBasis12 basis;
    for (uint32_t column = dense_base;
         column < heavy_base && out.Qd.size() < D2;
         ++column)
    {
        const uint16_t value = out.DenseIncidence[column];
        if (basis.Insert(value))
        {
            out.Qd.push_back(column);
            out.QdMasks.push_back(value);
        }
    }
    for (uint32_t source = 0u;
         source < K && out.Qd.size() < D2;
         ++source)
    {
        if ((source & 4095u) == 0u && deadline && deadline->Expired())
        {
            error = "deadline";
            return EvaluationOutcome::Deadline;
        }
        const uint16_t value = source_vector(source);
        if (basis.Insert(value))
        {
            out.Qd.push_back(source);
            out.QdMasks.push_back(value);
        }
    }
    out.SchurRank = basis.Rank();
    if (out.Qd.size() != D2 || out.SchurRank != D2)
    {
        error = "Schur construction did not find 12 pivots";
        return EvaluationOutcome::Reject;
    }

    out.IsQ.assign(L, 0u);
    out.IsQd.assign(L, 0u);
    out.Q.reserve(P);
    for (uint32_t column = K; column < dense_base; ++column)
    {
        out.IsQ[column] = 1u;
        out.Q.push_back(column);
    }
    for (uint32_t column : out.Qd)
    {
        out.IsQ[column] = 1u;
        out.IsQd[column] = 1u;
        out.QSourceCount += column < K ? 1u : 0u;
        out.Q.push_back(column);
    }
    for (uint32_t column = heavy_base; column < L; ++column)
    {
        out.IsQ[column] = 1u;
        out.Q.push_back(column);
    }
    for (uint32_t column = 0u; column < L; ++column) {
        if (!out.IsQ[column]) {
            out.F.push_back(column);
        }
    }
    if (out.Q.size() != P || out.F.size() != K)
    {
        error = "information-set partition size mismatch";
        return EvaluationOutcome::Invalid;
    }
    if (!RankPrecodeRestriction(
            system, out.Q, deadline, out.FullQRank, error))
    {
        return error == "deadline" ?
            EvaluationOutcome::Deadline : EvaluationOutcome::Invalid;
    }
    if (out.FullQRank != P)
    {
        error = "full precode Q submatrix is singular";
        return EvaluationOutcome::Reject;
    }
    return EvaluationOutcome::Pass;
}

std::string DigestPrecode(const wirehair_v2::PrecodeSystem& system)
{
    DigestBuilder digest;
    digest.Tag("tri-bmq4-precode-equations-v1");
    const wirehair_v2::PrecodeParams& p = system.Params;
    digest.U32(p.BlockCount);
    digest.U32(p.Staircase);
    digest.U32(p.DenseRows);
    digest.U32(p.HeavyRows);
    digest.U32(p.SourceHits);
    digest.U8(p.DenseIdentityCorner ? 1u : 0u);
    digest.U32(static_cast<uint32_t>(p.HeavyFamily));
    digest.U32(static_cast<uint32_t>(p.DenseAnchors));
    digest.U64(p.Seed);
    digest.U32(static_cast<uint32_t>(system.StaircaseRows.size()));
    for (const std::vector<uint32_t>& row : system.StaircaseRows)
    {
        digest.U32(static_cast<uint32_t>(row.size()));
        for (uint32_t column : row) {
            digest.U32(column);
        }
    }
    digest.U32(static_cast<uint32_t>(system.DenseBasisRowColumns.size()));
    for (const std::vector<uint32_t>& row : system.DenseBasisRowColumns)
    {
        digest.U32(static_cast<uint32_t>(row.size()));
        for (uint32_t column : row) {
            digest.U32(column);
        }
    }
    const uint32_t L = p.BlockCount + p.Staircase + p.DenseRows + p.HeavyRows;
    for (uint32_t heavy = 0u; heavy < p.HeavyRows; ++heavy) {
        for (uint32_t column = 0u; column < L; ++column) {
            digest.U8(wirehair_v2::HeavyCoefficientForParams(
                p, heavy, column));
        }
    }
    return digest.Finish();
}

uint32_t UniformBelow(wirehair::PCGRandom& prng, uint32_t bound)
{
    if (bound <= 1u) {
        return 0u;
    }
    const uint32_t threshold = (0u - bound) % bound;
    for (;;)
    {
        const uint32_t random = prng.Next();
        const uint64_t product = static_cast<uint64_t>(random) * bound;
        if (static_cast<uint32_t>(product) >= threshold) {
            return static_cast<uint32_t>(product >> 32);
        }
    }
}

wirehair::PCGRandom CandidatePrng(
    uint32_t K,
    uint64_t precode_seed,
    uint32_t packet_seed,
    uint64_t salt)
{
    wirehair::PCGRandom prng;
    const uint64_t sequence = precode_seed ^ salt;
    const uint64_t state =
        (static_cast<uint64_t>(packet_seed) << 32) | K;
    prng.Seed(sequence, state);
    return prng;
}

std::vector<uint32_t> BuildDeck(
    uint32_t count,
    uint64_t precode_seed,
    uint32_t packet_seed,
    uint64_t salt)
{
    std::vector<uint32_t> deck(count);
    if (count == 0u) {
        return deck;
    }
    wirehair::PCGRandom prng = CandidatePrng(
        count, precode_seed, packet_seed, salt);
    deck[0] = 0u;
    for (uint32_t index = 1u; index < count; ++index)
    {
        const uint32_t selected = UniformBelow(prng, index);
        deck[index] = deck[selected];
        deck[selected] = index;
    }
    return deck;
}

bool ValidatePermutation(
    const std::vector<uint32_t>& values,
    std::string& error)
{
    std::vector<uint8_t> seen(values.size(), 0u);
    for (uint32_t value : values)
    {
        if (value >= values.size() || seen[value])
        {
            error = "invalid deck permutation";
            return false;
        }
        seen[value] = 1u;
    }
    return true;
}

struct Topology
{
    std::vector<uint32_t> FreeDeck;
    std::vector<uint32_t> RowDeck;
    std::vector<uint32_t> DegreeDeck;
    std::vector<uint32_t> FreeColumn;
    std::vector<uint32_t> PublicId;
    std::vector<uint32_t> InversePublicId;
    std::vector<uint32_t> Degree;
    std::vector<size_t> ReferenceOffsets;
    std::vector<uint32_t> References;
    uint64_t ReferenceCount = 0u;
    uint64_t FifoReferenceCount = 0u;
    uint64_t ExcessReferenceCount = 0u;
    uint32_t PendingMaximum = 0u;
    uint32_t PendingTerminal = 0u;
    uint32_t PendingNeedTerminal = 0u;
};

struct PendingDemand
{
    uint32_t Position;
};

EvaluationOutcome BuildTopology(
    const Partition& partition,
    uint64_t precode_seed,
    uint32_t packet_seed,
    Deadline* deadline,
    Topology& out,
    std::string& error)
{
    out = Topology();
    if (partition.F.size() > UINT16_MAX)
    {
        error = "topology K outside frozen domain";
        return EvaluationOutcome::Invalid;
    }
    const uint32_t K = static_cast<uint32_t>(partition.F.size());
    if (K < 2u || K > UINT16_MAX)
    {
        error = "topology K outside frozen domain";
        return EvaluationOutcome::Invalid;
    }
    out.FreeDeck = BuildDeck(K, precode_seed, packet_seed, kFreeSalt);
    out.RowDeck = BuildDeck(K, precode_seed, packet_seed, kRowSalt);
    out.DegreeDeck = BuildDeck(K, precode_seed, packet_seed, kDegreeSalt);
    if (!ValidatePermutation(out.FreeDeck, error) ||
        !ValidatePermutation(out.RowDeck, error) ||
        !ValidatePermutation(out.DegreeDeck, error))
    {
        return EvaluationOutcome::Invalid;
    }
    out.FreeColumn.resize(K);
    out.PublicId.resize(K);
    out.InversePublicId.assign(K, UINT32_MAX);
    for (uint32_t position = 0u; position < K; ++position)
    {
        out.FreeColumn[position] = partition.F[out.FreeDeck[position]];
        out.PublicId[position] = out.RowDeck[position];
        if (out.InversePublicId[out.PublicId[position]] != UINT32_MAX)
        {
            error = "duplicate public row id";
            return EvaluationOutcome::Invalid;
        }
        out.InversePublicId[out.PublicId[position]] = position;
    }

    out.Degree.resize(K);
    out.ReferenceOffsets.reserve(static_cast<size_t>(K) + 1u);
    out.ReferenceOffsets.push_back(0u);
    out.References.reserve(static_cast<size_t>(K) * 5u);
    std::deque<PendingDemand> pending;
    std::vector<uint8_t> need(K, 0u);
    std::vector<uint32_t> chosen_stamp(K, 0u);
    wirehair::PCGRandom excess_prng = CandidatePrng(
        K, precode_seed, packet_seed, kExcessSalt);

    for (uint32_t position = 0u; position < K; ++position)
    {
        if ((position & 1023u) == 0u && deadline && deadline->Expired())
        {
            error = "deadline";
            return EvaluationOutcome::Deadline;
        }
        const uint32_t roster = out.DegreeDeck[position];
        const uint64_t numerator =
            (static_cast<uint64_t>(2u * roster + 1u) << 31);
        const uint32_t random_value =
            static_cast<uint32_t>(numerator / K);
        const uint32_t weight = wirehair::GeneratePeelRowWeight(
            random_value, static_cast<uint16_t>(K));
        if (weight == 0u)
        {
            error = "zero production row weight";
            return EvaluationOutcome::Invalid;
        }
        const uint32_t degree = std::min(weight - 1u, position);
        out.Degree[position] = degree;
        const uint32_t stamp = position + 1u;

        const size_t pending_snapshot = pending.size();
        const uint32_t fifo_count = std::min<uint32_t>(
            degree, static_cast<uint32_t>(pending_snapshot));
        for (uint32_t index = 0u; index < fifo_count; ++index)
        {
            const PendingDemand demand = pending.front();
            pending.pop_front();
            if (demand.Position >= position ||
                need[demand.Position] == 0u ||
                chosen_stamp[demand.Position] == stamp)
            {
                error = "FIFO demand invariant failed";
                return EvaluationOutcome::Invalid;
            }
            chosen_stamp[demand.Position] = stamp;
            out.References.push_back(demand.Position);
            --need[demand.Position];
            if (need[demand.Position] != 0u) {
                pending.push_back(demand);
            }
        }
        out.FifoReferenceCount += fifo_count;

        uint32_t remaining = degree - fifo_count;
        if (remaining != 0u)
        {
            if (position == 0u)
            {
                error = "excess demand has empty prior domain";
                return EvaluationOutcome::Invalid;
            }
            const uint32_t start = UniformBelow(excess_prng, position);
            for (uint32_t scanned = 0u;
                 scanned < position && remaining != 0u;
                 ++scanned)
            {
                uint32_t candidate = start + scanned;
                if (candidate >= position) {
                    candidate -= position;
                }
                if (need[candidate] != 0u ||
                    chosen_stamp[candidate] == stamp)
                {
                    continue;
                }
                chosen_stamp[candidate] = stamp;
                out.References.push_back(candidate);
                --remaining;
                ++out.ExcessReferenceCount;
            }
            if (remaining != 0u)
            {
                error = "excess cyclic scan exhausted";
                return EvaluationOutcome::Invalid;
            }
        }

        const size_t row_begin = out.ReferenceOffsets.back();
        if (out.References.size() - row_begin != degree)
        {
            error = "topology row degree mismatch";
            return EvaluationOutcome::Invalid;
        }
        for (size_t index = row_begin;
             index < out.References.size();
             ++index)
        {
            if (out.References[index] >= position)
            {
                error = "non-triangular topology reference";
                return EvaluationOutcome::Invalid;
            }
        }
        out.ReferenceOffsets.push_back(out.References.size());
        out.ReferenceCount += degree;
        need[position] = 4u;
        pending.push_back(PendingDemand{position});
        out.PendingMaximum = std::max<uint32_t>(
            out.PendingMaximum, static_cast<uint32_t>(pending.size()));
    }

    out.PendingTerminal = static_cast<uint32_t>(pending.size());
    uint32_t need_total = 0u;
    for (const PendingDemand& demand : pending)
    {
        if (demand.Position >= K || need[demand.Position] == 0u)
        {
            error = "terminal FIFO demand invariant failed";
            return EvaluationOutcome::Invalid;
        }
        need_total += need[demand.Position];
    }
    out.PendingNeedTerminal = need_total;
    if (out.ReferenceCount != out.References.size() ||
        out.ReferenceCount !=
            out.FifoReferenceCount + out.ExcessReferenceCount)
    {
        error = "topology reference accounting mismatch";
        return EvaluationOutcome::Invalid;
    }
    return EvaluationOutcome::Pass;
}

EvaluationOutcome EvaluateTopologyTerminalGate(
    const Topology& topology,
    std::string& error)
{
    if (topology.PendingTerminal <= 64u) {
        return EvaluationOutcome::Pass;
    }
    error = "terminal pending count exceeds 64";
    return EvaluationOutcome::Reject;
}

std::string DigestVector(
    const char* tag,
    const std::vector<uint32_t>& values)
{
    DigestBuilder digest;
    digest.Tag(tag);
    digest.U64(static_cast<uint64_t>(values.size()));
    for (uint32_t value : values) {
        digest.U32(value);
    }
    return digest.Finish();
}

std::string DigestPartition(const Partition& partition)
{
    if (partition.Qd.size() != partition.QdMasks.size()) {
        throw std::runtime_error("partition digest shape mismatch");
    }
    DigestBuilder digest;
    digest.Tag("tri-bmq4-partition-v1");
    digest.U32(static_cast<uint32_t>(partition.Qd.size()));
    for (size_t i = 0u; i < partition.Qd.size(); ++i)
    {
        digest.U32(partition.Qd[i]);
        digest.U32(partition.QdMasks[i]);
    }
    digest.U32(static_cast<uint32_t>(partition.Q.size()));
    for (uint32_t column : partition.Q) {
        digest.U32(column);
    }
    digest.U32(static_cast<uint32_t>(partition.F.size()));
    for (uint32_t column : partition.F) {
        digest.U32(column);
    }
    return digest.Finish();
}

std::string DigestTopology(const Topology& topology)
{
    DigestBuilder digest;
    digest.Tag("tri-bmq4-topology-v1");
    digest.U32(static_cast<uint32_t>(topology.FreeColumn.size()));
    for (uint32_t position = 0u;
         position < topology.FreeColumn.size();
         ++position)
    {
        digest.U32(position);
        digest.U32(topology.FreeColumn[position]);
        digest.U32(topology.PublicId[position]);
        digest.U32(topology.Degree[position]);
        const size_t first = topology.ReferenceOffsets[position];
        const size_t last = topology.ReferenceOffsets[position + 1u];
        digest.U32(static_cast<uint32_t>(last - first));
        for (size_t index = first; index < last; ++index) {
            digest.U32(topology.References[index]);
        }
    }
    digest.U64(topology.ReferenceCount);
    digest.U32(topology.PendingMaximum);
    digest.U32(topology.PendingTerminal);
    digest.U32(topology.PendingNeedTerminal);
    return digest.Finish();
}

bool BuildSystematicRow(
    const Topology& topology,
    uint32_t public_id,
    std::vector<uint32_t>& row,
    std::string& error)
{
    row.clear();
    if (public_id >= topology.InversePublicId.size())
    {
        error = "systematic public id out of range";
        return false;
    }
    const uint32_t position = topology.InversePublicId[public_id];
    if (position >= topology.FreeColumn.size())
    {
        error = "invalid inverse public row deck";
        return false;
    }
    row.reserve(static_cast<size_t>(topology.Degree[position]) + 1u);
    row.push_back(topology.FreeColumn[position]);
    const size_t first = topology.ReferenceOffsets[position];
    const size_t last = topology.ReferenceOffsets[position + 1u];
    for (size_t index = first; index < last; ++index) {
        row.push_back(topology.FreeColumn[topology.References[index]]);
    }
    std::sort(row.begin(), row.end());
    if (std::adjacent_find(row.begin(), row.end()) != row.end())
    {
        error = "duplicate systematic row column";
        return false;
    }
    return true;
}

struct DeliveredRows
{
    ExplicitRows BinaryRows;
    std::vector<uint8_t> PacketCovered;
    uint32_t SystematicCount = 0u;
    uint32_t RepairCount = 0u;
    std::string SystematicDigest;
    std::string RepairDigest;
    std::string RepairOracleDigest;
    std::string PacketDigest;
};

bool BuildDeliveredRows(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const std::vector<uint32_t>& packet_ids,
    const Topology& topology,
    Deadline* deadline,
    DeliveredRows& out,
    std::string& error)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint32_t L = K + P;
    if (packet_ids.size() != K)
    {
        error = "OH0 packet prefix does not contain K rows";
        return false;
    }
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount))
    {
        error = "packet row runtime initialization failed";
        return false;
    }
    out = DeliveredRows();
    size_t estimated_references = 0u;
    for (const std::vector<uint32_t>& row : system.StaircaseRows) {
        estimated_references += row.size();
    }
    for (const std::vector<uint32_t>& row : system.DenseBasisRowColumns) {
        estimated_references += row.size();
    }
    estimated_references += static_cast<size_t>(K) * 9u;
    out.BinaryRows.Reserve(
        static_cast<size_t>(system.Params.Staircase) +
            system.Params.DenseRows + K,
        estimated_references);
    for (const std::vector<uint32_t>& row : system.StaircaseRows) {
        out.BinaryRows.Append(row);
    }
    for (const std::vector<uint32_t>& row : system.DenseBasisRowColumns) {
        out.BinaryRows.Append(row);
    }
    out.PacketCovered.assign(L, 0u);
    DigestBuilder systematic_digest;
    DigestBuilder repair_digest;
    DigestBuilder repair_oracle_digest;
    DigestBuilder packet_digest;
    systematic_digest.Tag("tri-bmq4-systematic-rows-v1");
    repair_digest.Tag("tri-bmq4-repair-rows-v1");
    repair_oracle_digest.Tag("tri-bmq4-repair-rows-v1");
    packet_digest.Tag("tri-bmq4-delivered-rows-v1");

    std::vector<uint32_t> row;
    for (uint32_t packet_index = 0u; packet_index < K; ++packet_index)
    {
        if ((packet_index & 1023u) == 0u && deadline && deadline->Expired())
        {
            error = "deadline";
            return false;
        }
        const uint32_t block_id = packet_ids[packet_index];
        if (block_id < K)
        {
            if (!BuildSystematicRow(topology, block_id, row, error)) {
                return false;
            }
            ++out.SystematicCount;
            systematic_digest.U32(block_id);
            systematic_digest.U32(static_cast<uint32_t>(row.size()));
            for (uint32_t column : row) {
                systematic_digest.U32(column);
            }
        }
        else
        {
            row = wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                K, P, block_id, config, runtime);
            const std::vector<uint32_t> oracle =
                wirehair_v2::GeneratePacketMatrixRow(
                    K, P, block_id, config);
            if (row.empty() || row != oracle)
            {
                error = "cached versus public repair-row oracle mismatch";
                return false;
            }
            ++out.RepairCount;
            repair_digest.U32(block_id);
            repair_digest.U32(static_cast<uint32_t>(row.size()));
            repair_oracle_digest.U32(block_id);
            repair_oracle_digest.U32(static_cast<uint32_t>(oracle.size()));
            for (uint32_t column : row) {
                repair_digest.U32(column);
            }
            for (uint32_t column : oracle) {
                repair_oracle_digest.U32(column);
            }
        }
        packet_digest.U32(block_id);
        packet_digest.U32(static_cast<uint32_t>(row.size()));
        for (uint32_t column : row)
        {
            if (column >= L)
            {
                error = "delivered row column out of range";
                return false;
            }
            out.PacketCovered[column] = 1u;
            packet_digest.U32(column);
        }
        out.BinaryRows.Append(row);
    }
    if (out.SystematicCount + out.RepairCount != K)
    {
        error = "delivered row classification mismatch";
        return false;
    }
    if (!ValidateRows(out.BinaryRows, L, error)) {
        return false;
    }
    out.SystematicDigest = systematic_digest.Finish();
    out.RepairDigest = repair_digest.Finish();
    out.RepairOracleDigest = repair_oracle_digest.Finish();
    out.PacketDigest = packet_digest.Finish();
    return true;
}

struct WorkReceipt
{
    uint64_t Free = 0u;
    uint64_t Triangular = 0u;
    uint64_t Staircase = 0u;
    uint64_t DenseKnown = 0u;
    uint64_t SchurSolveXors = 0u;
    uint64_t SchurSolveCopies = 0u;
    uint64_t SchurSolve = 0u;
    uint64_t StaircasePatch = 0u;
    uint64_t HeavyKnownXors = 0u;
    uint64_t HeavyKnownMulAdds = 0u;
    uint64_t HeavyKnown = 0u;
    uint64_t HeavySolveXors = 0u;
    uint64_t HeavySolveMulAdds = 0u;
    uint64_t HeavySolveDivides = 0u;
    uint64_t HeavySolveCopies = 0u;
    uint64_t HeavySolve = 0u;
    uint64_t Emit = 0u;
    uint64_t XorCopyTotal = 0u;
    uint64_t Gf256Total = 0u;
    uint64_t Total = 0u;
    uint64_t SourceOwnership = 0u;
    uint64_t XorCopyTotalWithOwnedSource = 0u;
    uint64_t TotalWithOwnedSource = 0u;
};

EvaluationOutcome CountSchurSolveOps(
    const Partition& partition,
    uint64_t& xors,
    uint64_t& copies,
    std::string& error)
{
    xors = 0u;
    copies = 0u;
    if (partition.QdMasks.size() != kDenseRows)
    {
        error = "Schur work matrix shape is invalid";
        return EvaluationOutcome::Invalid;
    }
    std::array<uint16_t, kDenseRows> row_masks;
    row_masks.fill(0u);
    for (uint32_t column = 0u; column < kDenseRows; ++column) {
        for (uint32_t row = 0u; row < kDenseRows; ++row) {
            if ((partition.QdMasks[column] & (1u << row)) != 0u) {
                row_masks[row] ^= static_cast<uint16_t>(1u << column);
            }
        }
    }
    std::array<uint8_t, kDenseRows> used;
    used.fill(0u);
    for (uint32_t column = 0u; column < kDenseRows; ++column)
    {
        const uint16_t bit = static_cast<uint16_t>(1u << column);
        uint32_t pivot = kDenseRows;
        for (uint32_t row = 0u; row < kDenseRows; ++row) {
            if (!used[row] && (row_masks[row] & bit) != 0u) {
                pivot = row;
                break;
            }
        }
        if (pivot == kDenseRows)
        {
            error = "Schur work solve found singular matrix";
            return EvaluationOutcome::Invalid;
        }
        used[pivot] = 1u;
        for (uint32_t row = 0u; row < kDenseRows; ++row) {
            if (row != pivot && (row_masks[row] & bit) != 0u) {
                row_masks[row] ^= row_masks[pivot];
                ++xors;
            }
        }
    }
    copies = kDenseRows;
    if (xors + copies > 156u)
    {
        error = "Schur solve work exceeds frozen bound";
        return EvaluationOutcome::Reject;
    }
    return EvaluationOutcome::Pass;
}

EvaluationOutcome CountHeavySolveOps(
    const wirehair_v2::PrecodeSystem& system,
    uint64_t& xors,
    uint64_t& muladds,
    uint64_t& divides,
    uint64_t& copies,
    std::string& error)
{
    const uint32_t H = system.Params.HeavyRows;
    xors = 0u;
    muladds = 0u;
    divides = 0u;
    copies = 0u;
    if (H != kHeavyRows ||
        system.Params.HeavyFamily !=
            wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy)
    {
        error = "heavy work matrix shape is invalid";
        return EvaluationOutcome::Invalid;
    }
    const uint32_t heavy_base = system.Params.BlockCount +
        system.Params.Staircase + system.Params.DenseRows;
    std::vector<uint8_t> corner(static_cast<size_t>(H) * H);
    for (uint32_t row = 0u; row < H; ++row) {
        for (uint32_t column = 0u; column < H; ++column) {
            corner[static_cast<size_t>(row) * H + column] =
                wirehair_v2::HeavyCoefficientForParams(
                    system.Params, row, heavy_base + column);
        }
    }
    std::vector<uint8_t> used(H, 0u);
    for (uint32_t column = 0u; column < H; ++column)
    {
        uint32_t pivot = H;
        for (uint32_t row = 0u; row < H; ++row) {
            if (!used[row] &&
                corner[static_cast<size_t>(row) * H + column] != 0u)
            {
                pivot = row;
                break;
            }
        }
        if (pivot == H)
        {
            error = "heavy work solve found singular corner";
            return EvaluationOutcome::Invalid;
        }
        used[pivot] = 1u;
        const uint8_t pivot_value =
            corner[static_cast<size_t>(pivot) * H + column];
        if (pivot_value != 1u)
        {
            const uint8_t inverse = gf256_inv(pivot_value);
            for (uint32_t j = 0u; j < H; ++j) {
                corner[static_cast<size_t>(pivot) * H + j] = gf256_mul(
                    corner[static_cast<size_t>(pivot) * H + j], inverse);
            }
            ++divides;
        }
        for (uint32_t row = 0u; row < H; ++row)
        {
            const uint8_t scale =
                corner[static_cast<size_t>(row) * H + column];
            if (row == pivot || scale == 0u) {
                continue;
            }
            for (uint32_t j = 0u; j < H; ++j) {
                corner[static_cast<size_t>(row) * H + j] ^= gf256_mul(
                    scale,
                    corner[static_cast<size_t>(pivot) * H + j]);
            }
            if (scale == 1u) {
                ++xors;
            }
            else {
                ++muladds;
            }
        }
    }
    copies = H;
    if (xors + muladds + divides + copies > 156u)
    {
        error = "heavy solve work exceeds frozen bound";
        return EvaluationOutcome::Reject;
    }
    return EvaluationOutcome::Pass;
}

EvaluationOutcome ComputeWork(
    const wirehair_v2::PrecodeSystem& system,
    const Partition& partition,
    const Topology& topology,
    WorkReceipt& work,
    std::string& error)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint32_t N1 = system.Params.SourceHits;
    const uint32_t heavy_base = K + S + D2;
    const uint32_t L = heavy_base + H;
    work = WorkReceipt();
    const uint64_t source_hit_count = static_cast<uint64_t>(K) * N1;
    if (K < 2u || S == 0u || D2 != kDenseRows || H != kHeavyRows ||
        N1 == 0u || partition.QSourceCount > K ||
        partition.Qd.size() != D2 || partition.QdMasks.size() != D2 ||
        partition.F.size() != K || partition.IsQd.size() != L ||
        source_hit_count != partition.SourceHitRows.size() ||
        system.DenseBasisRowColumns.size() != D2 ||
        topology.ReferenceCount != topology.References.size() ||
        topology.Degree.size() != K ||
        topology.ReferenceOffsets.size() != static_cast<size_t>(K) + 1u ||
        topology.ReferenceOffsets.front() != 0u ||
        topology.ReferenceOffsets.back() != topology.References.size())
    {
        error = "work accounting input shape is invalid";
        return EvaluationOutcome::Invalid;
    }
    uint32_t counted_q_sources = 0u;
    for (uint32_t column : partition.Qd)
    {
        if (column >= heavy_base || !partition.IsQd[column])
        {
            error = "work partition pivot invariant failed";
            return EvaluationOutcome::Invalid;
        }
        counted_q_sources += column < K ? 1u : 0u;
    }
    if (counted_q_sources != partition.QSourceCount)
    {
        error = "work Q-source count invariant failed";
        return EvaluationOutcome::Invalid;
    }
    for (uint32_t position = 0u; position < K; ++position)
    {
        const size_t first = topology.ReferenceOffsets[position];
        const size_t last = topology.ReferenceOffsets[position + 1u];
        if (first > last || last > topology.References.size() ||
            last - first != topology.Degree[position])
        {
            error = "work topology offset invariant failed";
            return EvaluationOutcome::Invalid;
        }
        for (size_t index = first; index < last; ++index) {
            if (topology.References[index] >= position)
            {
                error = "work topology triangular invariant failed";
                return EvaluationOutcome::Invalid;
            }
        }
    }
    work.Free = K;
    work.Triangular = topology.ReferenceCount;
    work.Staircase =
        static_cast<uint64_t>(N1) * (K - partition.QSourceCount) + (S - 1u);
    for (const std::vector<uint32_t>& row : system.DenseBasisRowColumns) {
        for (uint32_t column : row) {
            if (column >= partition.IsQd.size())
            {
                error = "work dense column outside precode span";
                return EvaluationOutcome::Invalid;
            }
            if (!partition.IsQd[column]) {
                ++work.DenseKnown;
            }
        }
    }
    const EvaluationOutcome schur_outcome = CountSchurSolveOps(
            partition,
            work.SchurSolveXors,
            work.SchurSolveCopies,
            error);
    if (schur_outcome != EvaluationOutcome::Pass) {
        return schur_outcome;
    }
    work.SchurSolve = work.SchurSolveXors + work.SchurSolveCopies;
    for (uint32_t column : partition.Qd)
    {
        if (column >= K) {
            continue;
        }
        uint32_t weight = 0u;
        uint8_t parity = 0u;
        uint32_t hit_index = 0u;
        for (uint32_t row = 0u; row < S; ++row)
        {
            while (hit_index < N1 &&
                   partition.SourceHitRows[
                       static_cast<size_t>(column) * N1 + hit_index] == row)
            {
                parity ^= 1u;
                ++hit_index;
            }
            weight += parity;
        }
        if (hit_index != N1)
        {
            error = "source-pivot staircase patch hit mismatch";
            return EvaluationOutcome::Invalid;
        }
        work.StaircasePatch += weight;
    }
    if (heavy_base < 488u) {
        work.HeavyKnownMulAdds = static_cast<uint64_t>(H) * heavy_base;
    }
    else
    {
        work.HeavyKnownXors = heavy_base;
        work.HeavyKnownMulAdds = static_cast<uint64_t>(H) * 244u;
    }
    work.HeavyKnown = work.HeavyKnownXors + work.HeavyKnownMulAdds;
    const EvaluationOutcome heavy_outcome = CountHeavySolveOps(
            system,
            work.HeavySolveXors,
            work.HeavySolveMulAdds,
            work.HeavySolveDivides,
            work.HeavySolveCopies,
            error);
    if (heavy_outcome != EvaluationOutcome::Pass) {
        return heavy_outcome;
    }
    work.HeavySolve = work.HeavySolveXors + work.HeavySolveMulAdds +
        work.HeavySolveDivides + work.HeavySolveCopies;
    work.Emit = K;
    work.XorCopyTotal = work.Free + work.Triangular + work.Staircase +
        work.DenseKnown + work.SchurSolveXors + work.SchurSolveCopies +
        work.StaircasePatch + work.HeavyKnownXors + work.HeavySolveXors +
        work.HeavySolveCopies + work.Emit;
    work.Gf256Total = work.HeavyKnownMulAdds + work.HeavySolveMulAdds +
        work.HeavySolveDivides;
    work.Total = work.XorCopyTotal + work.Gf256Total;
    work.SourceOwnership = K;
    work.XorCopyTotalWithOwnedSource =
        work.XorCopyTotal + work.SourceOwnership;
    work.TotalWithOwnedSource = work.Total + work.SourceOwnership;
    if (work.TotalWithOwnedSource >= static_cast<uint64_t>(13u) * K)
    {
        error = "owned-source exact block work is not below 13K";
        return EvaluationOutcome::Reject;
    }
    return EvaluationOutcome::Pass;
}

struct PeelRowState
{
    uint16_t Live = 0u;
    uint16_t LowDegreeXor = 0u;
};

struct PeelOutcome
{
    std::vector<uint32_t> SolveRow;
    std::vector<uint32_t> PeelOrder;
    std::vector<uint32_t> InactiveOrder;
    std::vector<uint8_t> UsedRows;
    bool Complete = false;
    bool CoreLimitExceeded = false;
    bool TimedOut = false;
};

uint64_t DirectDegreeTwoKey(
    uint32_t live_references,
    uint32_t total_references,
    uint32_t column)
{
    return static_cast<uint64_t>(live_references) << 32 |
        static_cast<uint64_t>(total_references) << 16 |
        (static_cast<uint32_t>(UINT16_MAX) - column);
}

uint32_t DirectDegreeTwoColumn(uint64_t key)
{
    return static_cast<uint32_t>(UINT16_MAX) -
        static_cast<uint16_t>(key);
}

bool PeelProductionOrder(
    uint32_t column_count,
    const ExplicitRows& rows,
    uint32_t inactive_limit,
    Deadline* deadline,
    PeelOutcome& out,
    std::string& error)
{
    out = PeelOutcome();
    if (column_count == 0u || column_count >= UINT16_MAX ||
        rows.size() > UINT32_MAX)
    {
        error = "production-order peeler domain is invalid";
        return false;
    }
    out.SolveRow.assign(column_count, UINT32_MAX);
    out.UsedRows.assign(rows.size(), 0u);
    out.PeelOrder.reserve(column_count);

    std::vector<PeelRowState> row_state(rows.size());
    std::vector<size_t> column_offsets(
        static_cast<size_t>(column_count) + 1u, 0u);
    std::vector<uint8_t> resolved(column_count, 0u);
    std::vector<uint32_t> queue;
    queue.reserve(rows.size());
    std::vector<uint32_t> degree_two_refs(column_count, 0u);
    class DegreeTwoQueue : public std::priority_queue<uint64_t>
    {
    public:
        void Reserve(size_t count) { this->c.reserve(count); }
    };
    DegreeTwoQueue degree_two_queue;
    degree_two_queue.Reserve(column_count);

    for (uint32_t row = 0u; row < rows.size(); ++row)
    {
        const RowView view = rows[row];
        if (view.size() > UINT16_MAX)
        {
            error = "explicit row exceeds uint16 degree";
            return false;
        }
        row_state[row].Live = static_cast<uint16_t>(view.size());
        if (row_state[row].Live == 1u)
        {
            queue.push_back(row);
            row_state[row].LowDegreeXor =
                static_cast<uint16_t>(*view.begin());
        }
        for (uint32_t column : view)
        {
            if (column >= column_count)
            {
                error = "peel column out of range";
                return false;
            }
            ++column_offsets[static_cast<size_t>(column) + 1u];
        }
    }
    for (uint32_t column = 0u; column < column_count; ++column) {
        column_offsets[static_cast<size_t>(column) + 1u] +=
            column_offsets[column];
    }
    std::vector<uint32_t> column_rows(column_offsets[column_count]);
    for (uint32_t row = 0u; row < rows.size(); ++row) {
        for (uint32_t column : rows[row]) {
            const size_t destination =
                column_offsets[column] + degree_two_refs[column]++;
            column_rows[destination] = row;
        }
    }
    std::fill(degree_two_refs.begin(), degree_two_refs.end(), 0u);

    uint32_t maximum_reference_count = 0u;
    for (uint32_t column = 0u; column < column_count; ++column) {
        const size_t count = column_offsets[static_cast<size_t>(column) + 1u] -
            column_offsets[column];
        if (count > UINT16_MAX)
        {
            error = "column reference count exceeds direct-key domain";
            return false;
        }
        maximum_reference_count = std::max<uint32_t>(
            maximum_reference_count, static_cast<uint32_t>(count));
    }
    std::vector<size_t> bucket_offsets(
        static_cast<size_t>(maximum_reference_count) + 2u, 0u);
    for (uint32_t column = 0u; column < column_count; ++column)
    {
        const uint32_t count = static_cast<uint32_t>(
            column_offsets[static_cast<size_t>(column) + 1u] -
            column_offsets[column]);
        ++bucket_offsets[static_cast<size_t>(count) + 1u];
    }
    for (uint32_t count = 0u; count <= maximum_reference_count; ++count) {
        bucket_offsets[static_cast<size_t>(count) + 1u] +=
            bucket_offsets[count];
    }
    std::vector<size_t> bucket_cursor = bucket_offsets;
    std::vector<uint32_t> reference_columns(column_count);
    for (uint32_t column = 0u; column < column_count; ++column)
    {
        const uint32_t count = static_cast<uint32_t>(
            column_offsets[static_cast<size_t>(column) + 1u] -
            column_offsets[column]);
        reference_columns[bucket_cursor[count]++] = column;
    }
    bucket_cursor = bucket_offsets;

    const auto degree_two_key = [&](uint32_t column) {
        const uint32_t total = static_cast<uint32_t>(
            column_offsets[static_cast<size_t>(column) + 1u] -
            column_offsets[column]);
        return DirectDegreeTwoKey(
            degree_two_refs[column], total, column);
    };
    const auto add_degree_two = [&](uint32_t row) {
        PeelRowState& state = row_state[row];
        if (state.Live != 2u) {
            return true;
        }
        uint32_t pair_xor = 0u;
        uint32_t pair_count = 0u;
        for (uint32_t column : rows[row])
        {
            if (resolved[column]) {
                continue;
            }
            pair_xor ^= column;
            ++pair_count;
            ++degree_two_refs[column];
            degree_two_queue.push(degree_two_key(column));
        }
        if (pair_count != 2u || pair_xor > UINT16_MAX) {
            return false;
        }
        state.LowDegreeXor = static_cast<uint16_t>(pair_xor);
        return true;
    };
    const auto remove_degree_two = [&](uint32_t row,
                                       uint32_t resolved_column) {
        PeelRowState& state = row_state[row];
        if (state.Live != 2u) {
            return true;
        }
        const uint32_t other =
            static_cast<uint32_t>(state.LowDegreeXor) ^ resolved_column;
        if (other >= column_count || resolved[other]) {
            return false;
        }
        if (degree_two_refs[other] == 0u) {
            return false;
        }
        --degree_two_refs[other];
        state.LowDegreeXor = static_cast<uint16_t>(other);
        return true;
    };
    for (uint32_t row = 0u; row < rows.size(); ++row) {
        if (!add_degree_two(row))
        {
            error = "initial degree-two invariant failed";
            return false;
        }
    }

    const auto resolve = [&](uint32_t column) {
        if (column >= column_count || resolved[column]) {
            return false;
        }
        resolved[column] = 1u;
        for (size_t reference = column_offsets[column];
             reference < column_offsets[static_cast<size_t>(column) + 1u];
             ++reference)
        {
            const uint32_t row = column_rows[reference];
            if (row_state[row].Live == 0u) {
                continue;
            }
            if (!remove_degree_two(row, column)) {
                return false;
            }
            --row_state[row].Live;
            if (!add_degree_two(row)) {
                return false;
            }
            if (row_state[row].Live == 1u) {
                queue.push_back(row);
            }
        }
        return true;
    };

    uint32_t remaining = column_count;
    size_t queue_head = 0u;
    uint64_t iterations = 0u;
    while (remaining > 0u)
    {
        if (((iterations++) & 1023u) == 0u &&
            deadline && deadline->Expired())
        {
            out.TimedOut = true;
            error = "deadline";
            return false;
        }
        while (queue_head < queue.size())
        {
            if ((queue_head & 1023u) == 0u &&
                deadline && deadline->Expired())
            {
                out.TimedOut = true;
                error = "deadline";
                return false;
            }
            const uint32_t row = queue[queue_head++];
            if (row_state[row].Live != 1u) {
                continue;
            }
            const uint32_t column = row_state[row].LowDegreeXor;
            if (column >= column_count || resolved[column])
            {
                error = "singleton low-degree XOR invariant failed";
                return false;
            }
            out.UsedRows[row] = 1u;
            out.SolveRow[column] = row;
            out.PeelOrder.push_back(column);
            if (!resolve(column) || row_state[row].Live != 0u)
            {
                error = "singleton resolution invariant failed";
                return false;
            }
            --remaining;
        }
        if (remaining == 0u) {
            break;
        }

        uint32_t best = UINT32_MAX;
        while (!degree_two_queue.empty())
        {
            const uint64_t candidate = degree_two_queue.top();
            const uint32_t column = DirectDegreeTwoColumn(candidate);
            if (column >= column_count || resolved[column] ||
                degree_two_refs[column] !=
                    static_cast<uint32_t>(candidate >> 32))
            {
                degree_two_queue.pop();
                continue;
            }
            best = column;
            break;
        }
        while (best == UINT32_MAX)
        {
            size_t& cursor = bucket_cursor[maximum_reference_count];
            const size_t end =
                bucket_offsets[static_cast<size_t>(maximum_reference_count) + 1u];
            while (cursor < end && resolved[reference_columns[cursor]]) {
                ++cursor;
            }
            if (cursor < end)
            {
                best = reference_columns[cursor];
                break;
            }
            if (maximum_reference_count == 0u) {
                break;
            }
            --maximum_reference_count;
        }
        if (best == UINT32_MAX)
        {
            error = "production-order peeler could not select a column";
            return false;
        }
        out.InactiveOrder.push_back(best);
        if (inactive_limit != 0u &&
            out.InactiveOrder.size() > inactive_limit)
        {
            out.CoreLimitExceeded = true;
            return true;
        }
        if (!resolve(best))
        {
            error = "inactive-column resolution invariant failed";
            return false;
        }
        --remaining;
    }
    out.Complete = out.PeelOrder.size() + out.InactiveOrder.size() ==
        column_count;
    if (!out.Complete)
    {
        error = "production-order peeler incomplete";
        return false;
    }
    return true;
}

bool PeelSlowReference(
    uint32_t column_count,
    const ExplicitRows& rows,
    PeelOutcome& out,
    std::string& error)
{
    out = PeelOutcome();
    if (column_count == 0u) {
        error = "slow peeler empty domain";
        return false;
    }
    out.SolveRow.assign(column_count, UINT32_MAX);
    out.UsedRows.assign(rows.size(), 0u);
    std::vector<uint32_t> live(rows.size(), 0u);
    std::vector<uint8_t> resolved(column_count, 0u);
    std::vector<std::vector<uint32_t> > column_rows(column_count);
    std::vector<uint32_t> queue;
    for (uint32_t row = 0u; row < rows.size(); ++row)
    {
        live[row] = static_cast<uint32_t>(rows[row].size());
        if (live[row] == 1u) {
            queue.push_back(row);
        }
        for (uint32_t column : rows[row]) {
            if (column >= column_count)
            {
                error = "slow peeler column out of range";
                return false;
            }
            column_rows[column].push_back(row);
        }
    }
    const auto resolve = [&](uint32_t column) {
        if (column >= column_count || resolved[column]) {
            return false;
        }
        resolved[column] = 1u;
        for (uint32_t row : column_rows[column]) {
            if (live[row] == 0u) {
                continue;
            }
            --live[row];
            if (live[row] == 1u) {
                queue.push_back(row);
            }
        }
        return true;
    };

    uint32_t remaining = column_count;
    size_t queue_head = 0u;
    while (remaining != 0u)
    {
        while (queue_head < queue.size())
        {
            const uint32_t row = queue[queue_head++];
            if (live[row] != 1u) {
                continue;
            }
            uint32_t column = UINT32_MAX;
            for (uint32_t candidate : rows[row]) {
                if (!resolved[candidate]) {
                    column = candidate;
                    break;
                }
            }
            if (column == UINT32_MAX)
            {
                error = "slow singleton scan failed";
                return false;
            }
            out.UsedRows[row] = 1u;
            out.SolveRow[column] = row;
            out.PeelOrder.push_back(column);
            if (!resolve(column))
            {
                error = "slow singleton resolve failed";
                return false;
            }
            --remaining;
        }
        if (remaining == 0u) {
            break;
        }
        uint32_t best = UINT32_MAX;
        uint32_t best_degree_two = 0u;
        uint32_t best_total = 0u;
        bool have_degree_two = false;
        for (uint32_t column = 0u; column < column_count; ++column)
        {
            if (resolved[column]) {
                continue;
            }
            uint32_t degree_two = 0u;
            for (uint32_t row : column_rows[column]) {
                degree_two += live[row] == 2u ? 1u : 0u;
            }
            if (degree_two != 0u) {
                have_degree_two = true;
            }
            const uint32_t total =
                static_cast<uint32_t>(column_rows[column].size());
            if (best == UINT32_MAX ||
                degree_two > best_degree_two ||
                (degree_two == best_degree_two && total > best_total) ||
                (degree_two == best_degree_two && total == best_total &&
                 column < best))
            {
                best = column;
                best_degree_two = degree_two;
                best_total = total;
            }
        }
        if (!have_degree_two)
        {
            best = UINT32_MAX;
            best_total = 0u;
            for (uint32_t column = 0u; column < column_count; ++column) {
                if (resolved[column]) {
                    continue;
                }
                const uint32_t total =
                    static_cast<uint32_t>(column_rows[column].size());
                if (best == UINT32_MAX || total > best_total ||
                    (total == best_total && column < best))
                {
                    best = column;
                    best_total = total;
                }
            }
        }
        if (best == UINT32_MAX)
        {
            error = "slow peeler could not select a column";
            return false;
        }
        out.InactiveOrder.push_back(best);
        if (!resolve(best))
        {
            error = "slow inactive resolve failed";
            return false;
        }
        --remaining;
    }
    out.Complete = true;
    return true;
}

bool SamePeel(const PeelOutcome& a, const PeelOutcome& b)
{
    return a.Complete == b.Complete &&
        a.SolveRow == b.SolveRow &&
        a.PeelOrder == b.PeelOrder &&
        a.InactiveOrder == b.InactiveOrder &&
        a.UsedRows == b.UsedRows;
}

class Gf2Basis
{
public:
    explicit Gf2Basis(uint32_t width)
        : Width(width), Words((width + 63u) / 64u), Pivot(width, UINT32_MAX)
    {
    }

    bool Insert(std::vector<uint64_t>& row)
    {
        for (;;)
        {
            uint32_t pivot_bit = Width;
            for (uint32_t word = 0u; word < Words; ++word) {
                if (row[word] != 0u)
                {
                    pivot_bit = word * 64u +
                        wirehair::NonzeroLowestBitIndex64(row[word]);
                    break;
                }
            }
            if (pivot_bit >= Width) {
                return false;
            }
            const uint32_t existing = Pivot[pivot_bit];
            if (existing == UINT32_MAX)
            {
                Pivot[pivot_bit] = static_cast<uint32_t>(Rows.size());
                Rows.push_back(row);
                return true;
            }
            const std::vector<uint64_t>& pivot = Rows[existing];
            for (uint32_t word = 0u; word < Words; ++word) {
                row[word] ^= pivot[word];
            }
        }
    }

private:
    uint32_t Width;
    uint32_t Words;
    std::vector<uint32_t> Pivot;
    std::vector<std::vector<uint64_t> > Rows;
};

struct ResidualReceipt
{
    uint32_t Rows = 0u;
    uint32_t BinaryRank = 0u;
    uint32_t Deficiency = 0u;
};

bool AnalyzeBinaryResidual(
    uint32_t column_count,
    const ExplicitRows& rows,
    const PeelOutcome& peel,
    Deadline* deadline,
    ResidualReceipt& receipt,
    std::string& error)
{
    receipt = ResidualReceipt();
    if (!peel.Complete ||
        peel.PeelOrder.size() + peel.InactiveOrder.size() != column_count ||
        peel.SolveRow.size() != column_count ||
        peel.UsedRows.size() != rows.size())
    {
        error = "residual analysis requires complete peel";
        return false;
    }
    std::vector<uint8_t> classified(column_count, 0u);
    for (uint32_t column : peel.InactiveOrder)
    {
        if (column >= column_count || classified[column])
        {
            error = "invalid inactive-column order";
            return false;
        }
        classified[column] = 1u;
    }
    for (uint32_t column : peel.PeelOrder)
    {
        if (column >= column_count || classified[column])
        {
            error = "invalid peeled-column order";
            return false;
        }
        classified[column] = 1u;
    }
    const uint32_t R = static_cast<uint32_t>(peel.InactiveOrder.size());
    const uint32_t words = (R + 63u) / 64u;
    if (R == 0u)
    {
        std::vector<uint8_t> projected(column_count, 0u);
        for (uint32_t solve_index = 0u;
             solve_index < peel.PeelOrder.size();
             ++solve_index)
        {
            if ((solve_index & 1023u) == 0u &&
                deadline && deadline->Expired())
            {
                error = "deadline";
                return false;
            }
            const uint32_t column = peel.PeelOrder[solve_index];
            const uint32_t solve_row = peel.SolveRow[column];
            if (solve_row >= rows.size())
            {
                error = "peeled column has invalid solve row";
                return false;
            }
            for (uint32_t dependency : rows[solve_row])
            {
                if (dependency == column) {
                    continue;
                }
                if (dependency >= column_count || !projected[dependency])
                {
                    error = "zero-width peeled dependency is unavailable";
                    return false;
                }
            }
            projected[column] = 1u;
        }
        for (uint32_t row = 0u; row < rows.size(); ++row)
        {
            if ((row & 1023u) == 0u && deadline && deadline->Expired())
            {
                error = "deadline";
                return false;
            }
            if (peel.UsedRows[row]) {
                continue;
            }
            for (uint32_t column : rows[row])
            {
                if (column >= column_count || !projected[column])
                {
                    error = "zero-width residual projection is unavailable";
                    return false;
                }
            }
            ++receipt.Rows;
        }
        return true;
    }
    std::vector<uint32_t> inactive_index(column_count, UINT32_MAX);
    for (uint32_t index = 0u; index < R; ++index) {
        inactive_index[peel.InactiveOrder[index]] = index;
    }
    if (words != 0u &&
        static_cast<uint64_t>(column_count) * words >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
    {
        error = "projection allocation overflow";
        return false;
    }
    std::vector<uint64_t> projection(
        static_cast<size_t>(column_count) * words, 0u);
    std::vector<uint8_t> projected(column_count, 0u);
    std::vector<uint64_t> accumulator(words, 0u);
    // This derives each peeled variable's inactive-column projection from
    // dependencies that were inactive or resolved earlier, so chronological
    // PeelOrder is intentional.  Production's reverse pass performs the
    // algebraically equivalent operation of propagating row coefficients
    // into dependencies instead.
    for (uint32_t solve_index = 0u;
         solve_index < peel.PeelOrder.size();
         ++solve_index)
    {
        if ((solve_index & 1023u) == 0u && deadline && deadline->Expired())
        {
            error = "deadline";
            return false;
        }
        std::fill(accumulator.begin(), accumulator.end(), 0u);
        const uint32_t column = peel.PeelOrder[solve_index];
        const uint32_t solve_row = peel.SolveRow[column];
        if (solve_row >= rows.size())
        {
            error = "peeled column has invalid solve row";
            return false;
        }
        for (uint32_t dependency : rows[solve_row])
        {
            if (dependency == column) {
                continue;
            }
            if (dependency >= column_count)
            {
                error = "peeled dependency column is out of range";
                return false;
            }
            const uint32_t inactive = inactive_index[dependency];
            if (inactive != UINT32_MAX) {
                accumulator[inactive >> 6] ^=
                    UINT64_C(1) << (inactive & 63u);
            }
            else
            {
                if (!projected[dependency])
                {
                    error = "peeled dependency projection is unavailable";
                    return false;
                }
                const uint64_t* source = projection.data() +
                    static_cast<size_t>(dependency) * words;
                for (uint32_t word = 0u; word < words; ++word) {
                    accumulator[word] ^= source[word];
                }
            }
        }
        uint64_t* destination = projection.data() +
            static_cast<size_t>(column) * words;
        std::copy(accumulator.begin(), accumulator.end(), destination);
        projected[column] = 1u;
    }

    Gf2Basis basis(R);
    for (uint32_t row = 0u; row < rows.size(); ++row)
    {
        if ((row & 1023u) == 0u && deadline && deadline->Expired())
        {
            error = "deadline";
            return false;
        }
        if (peel.UsedRows[row]) {
            continue;
        }
        ++receipt.Rows;
        std::fill(accumulator.begin(), accumulator.end(), 0u);
        for (uint32_t column : rows[row])
        {
            if (column >= column_count)
            {
                error = "residual column is out of range";
                return false;
            }
            const uint32_t inactive = inactive_index[column];
            if (inactive != UINT32_MAX) {
                accumulator[inactive >> 6] ^=
                    UINT64_C(1) << (inactive & 63u);
            }
            else
            {
                if (!projected[column])
                {
                    error = "residual projection is unavailable";
                    return false;
                }
                const uint64_t* source = projection.data() +
                    static_cast<size_t>(column) * words;
                for (uint32_t word = 0u; word < words; ++word) {
                    accumulator[word] ^= source[word];
                }
            }
        }
        if (basis.Insert(accumulator)) {
            ++receipt.BinaryRank;
        }
    }
    if (receipt.BinaryRank > R)
    {
        error = "binary residual rank exceeds inactive width";
        return false;
    }
    receipt.Deficiency = R - receipt.BinaryRank;
    return true;
}

uint32_t BruteForceGf2Rank(uint32_t columns, const ExplicitRows& rows)
{
    const uint32_t words = (columns + 63u) / 64u;
    Gf2Basis basis(columns);
    uint32_t rank = 0u;
    std::vector<uint64_t> bits(words, 0u);
    for (uint32_t row = 0u; row < rows.size(); ++row)
    {
        std::fill(bits.begin(), bits.end(), 0u);
        for (uint32_t column : rows[row]) {
            bits[column >> 6] ^= UINT64_C(1) << (column & 63u);
        }
        rank += basis.Insert(bits) ? 1u : 0u;
    }
    return rank;
}

class ReferencePcg
{
public:
    void Seed(uint64_t sequence, uint64_t state = 0u)
    {
        State = 0u;
        Increment = (sequence << 1u) | 1u;
        Next();
        State += state;
        Next();
    }

    uint32_t Next()
    {
        const uint64_t old_state = State;
        State = old_state * UINT64_C(6364136223846793005) + Increment;
        const uint32_t shifted = static_cast<uint32_t>(
            ((old_state >> 18u) ^ old_state) >> 27u);
        const uint32_t rotation = static_cast<uint32_t>(old_state >> 59u);
        return (shifted >> rotation) |
            (shifted << (static_cast<uint32_t>(
                -static_cast<int32_t>(rotation)) & 31u));
    }

private:
    uint64_t State = 0u;
    uint64_t Increment = 0u;
};

uint32_t RankGf256Scalar(
    std::vector<uint8_t> matrix,
    uint32_t rows,
    uint32_t columns)
{
    uint32_t rank = 0u;
    for (uint32_t column = 0u; column < columns && rank < rows; ++column)
    {
        uint32_t pivot = rank;
        while (pivot < rows &&
               matrix[static_cast<size_t>(pivot) * columns + column] == 0u)
        {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        for (uint32_t j = 0u; j < columns; ++j) {
            std::swap(matrix[static_cast<size_t>(rank) * columns + j],
                      matrix[static_cast<size_t>(pivot) * columns + j]);
        }
        const uint8_t inverse = gf256_inv(
            matrix[static_cast<size_t>(rank) * columns + column]);
        for (uint32_t j = column; j < columns; ++j) {
            matrix[static_cast<size_t>(rank) * columns + j] = gf256_mul(
                matrix[static_cast<size_t>(rank) * columns + j], inverse);
        }
        for (uint32_t target = 0u; target < rows; ++target)
        {
            if (target == rank) {
                continue;
            }
            const uint8_t scale =
                matrix[static_cast<size_t>(target) * columns + column];
            if (scale == 0u) {
                continue;
            }
            for (uint32_t j = column; j < columns; ++j) {
                matrix[static_cast<size_t>(target) * columns + j] ^= gf256_mul(
                    scale,
                    matrix[static_cast<size_t>(rank) * columns + j]);
            }
        }
        ++rank;
    }
    return rank;
}

bool SelfCheck(bool condition, const char* message, std::string& error)
{
    if (condition) {
        return true;
    }
    error = message;
    return false;
}

bool SelfTestPcgAndDeck(std::string& error)
{
    static const uint64_t sequences[] = {
        UINT64_C(0), UINT64_C(1), UINT64_C(0x487468302aad7105)
    };
    static const uint64_t states[] = {
        UINT64_C(0), UINT64_C(0x123456789abcdef0),
        UINT64_C(0x4ec7210200000101)
    };
    for (uint64_t sequence : sequences) {
        for (uint64_t state : states)
        {
            wirehair::PCGRandom production;
            ReferencePcg reference;
            production.Seed(sequence, state);
            reference.Seed(sequence, state);
            for (unsigned draw = 0u; draw < 100u; ++draw) {
                if (production.Next() != reference.Next()) {
                    error = "PCGRandom reference mismatch";
                    return false;
                }
            }
        }
    }

    const std::vector<uint32_t> first = BuildDeck(
        257u, UINT64_C(0x123456789abcdef0), UINT32_C(0x89abcdef),
        kFreeSalt);
    const std::vector<uint32_t> second = BuildDeck(
        257u, UINT64_C(0x123456789abcdef0), UINT32_C(0x89abcdef),
        kFreeSalt);
    if (!SelfCheck(first == second, "deck determinism failed", error) ||
        !ValidatePermutation(first, error))
    {
        return false;
    }
    std::vector<uint8_t> visited(first.size(), 0u);
    uint32_t cursor = 0u;
    for (uint32_t count = 0u; count < first.size(); ++count)
    {
        if (cursor >= first.size() || visited[cursor])
        {
            error = "Sattolo deck is not one cycle";
            return false;
        }
        visited[cursor] = 1u;
        cursor = first[cursor];
    }
    return SelfCheck(cursor == 0u, "Sattolo cycle does not close", error);
}

bool SelfTestGf256Rank(std::string& error)
{
    ReferencePcg rng;
    rng.Seed(UINT64_C(0x72616e6b2d6f7263), UINT64_C(0x1234));
    for (uint32_t rows = 1u; rows <= 20u; ++rows) {
        for (uint32_t columns = 1u; columns <= 20u; ++columns)
        {
            std::vector<uint8_t> matrix(
                static_cast<size_t>(rows) * columns);
            for (uint8_t& value : matrix) {
                value = static_cast<uint8_t>(rng.Next());
            }
            std::vector<uint8_t> optimized = matrix;
            bool timed_out = false;
            const uint32_t actual = RankGf256(
                optimized, rows, columns, NULL, timed_out);
            const uint32_t expected =
                RankGf256Scalar(matrix, rows, columns);
            if (timed_out || actual != expected)
            {
                error = "GF256 rank differential failed";
                return false;
            }
        }
    }
    return true;
}

bool SelfTestRawSeedConfiguration(std::string& error)
{
    if (!ValidateRawSeedLaw(error) || !ValidateContractVersions(error)) {
        return false;
    }

    wirehair_v2::PrecodeParams params;
    wirehair_v2::PacketRowConfig packet;
    wirehair_v2::test::MakeRawArchitectureConfiguration(
        127u, params, packet);
    static const uint32_t attempts[] = { 0u, 1u, 2u, 255u };
    static const uint64_t expected_precode[] = {
        UINT64_C(0x487468302aad7105),
        UINT64_C(0xe6abe1e9a9f7ed1a),
        UINT64_C(0x84e35ba32942692f),
        UINT64_C(0xe1b6a7f5f5df09f0)
    };
    static const uint32_t expected_packet[] = {
        UINT32_C(0x4ec72102), UINT32_C(0xecfe9abb),
        UINT32_C(0x8b361474), UINT32_C(0xe8096049)
    };
    if (params.BlockCount != 127u ||
        params.Seed != kExpectedPrecodeSeed ||
        packet.PeelSeed != kExpectedPacketSeed ||
        packet.MixCount != wirehair_v2::kCertifiedPacketMixCount ||
        wirehair::GetDenseCount(kK) != 346u)
    {
        error = "bounded raw configuration mismatch";
        return false;
    }
    for (size_t i = 0u; i < sizeof(attempts) / sizeof(attempts[0]); ++i)
    {
        if (wirehair_v2::PrecodeParamsForAttempt(
                params, attempts[i]).Seed != expected_precode[i] ||
            wirehair_v2::PacketConfigForAttempt(
                packet, attempts[i]).PeelSeed != expected_packet[i])
        {
            error = "uniform-raw-v1 attempt schedule mismatch";
            return false;
        }
    }
    return true;
}

bool SelfTestPartitionAndTopology(std::string& error)
{
    bool tested_partition = false;
    static const uint32_t candidates[] = { 32u, 100u, 128u, 512u };
    for (uint32_t K : candidates)
    {
        const wirehair_v2::SeedProfile profile =
            wirehair_v2::SelectSeedProfile(K, 2u);
        wirehair_v2::MessagePrecodeEncoderOptions options;
        wirehair_v2::PrecodeParams params;
        wirehair_v2::PacketRowConfig config;
        if (!wirehair_v2::ResolveMessagePrecodeConfiguration(
                profile, options, params, config))
        {
            error = "selftest profile configuration failed";
            return false;
        }
        wirehair_v2::PrecodeSystem system;
        if (!wirehair_v2::BuildPrecodeSystem(params, system))
        {
            error = "selftest precode construction failed";
            return false;
        }
        Partition partition;
        std::string local_error;
        const EvaluationOutcome partition_outcome =
            BuildPartition(system, NULL, partition, local_error);
        if (partition_outcome == EvaluationOutcome::Reject) {
            continue;
        }
        if (partition_outcome != EvaluationOutcome::Pass)
        {
            error = local_error.empty() ?
                "bounded partition returned invalid outcome" : local_error;
            return false;
        }
        if (partition.FullQRank !=
                params.Staircase + params.DenseRows + params.HeavyRows ||
            partition.F.size() != K || partition.Qd.size() != kDenseRows ||
            config.MixCount != wirehair_v2::kCertifiedPacketMixCount)
        {
            error = "selftest partition invariant failed";
            return false;
        }
        uint64_t schur_xors = 0u;
        uint64_t schur_copies = 0u;
        uint64_t heavy_xors = 0u;
        uint64_t heavy_muladds = 0u;
        uint64_t heavy_divides = 0u;
        uint64_t heavy_copies = 0u;
        if (CountSchurSolveOps(
                partition, schur_xors, schur_copies, error) !=
                    EvaluationOutcome::Pass ||
            CountHeavySolveOps(
                system, heavy_xors, heavy_muladds,
                heavy_divides, heavy_copies, error) !=
                    EvaluationOutcome::Pass ||
            schur_copies != kDenseRows ||
            schur_xors + schur_copies > 156u ||
            heavy_copies != kHeavyRows ||
            heavy_xors + heavy_muladds + heavy_divides + heavy_copies > 156u)
        {
            if (error.empty()) {
                error = "bounded split work accounting failed";
            }
            return false;
        }
        tested_partition = true;
        break;
    }
    if (!tested_partition)
    {
        error = "no bounded selftest precode admitted a Schur partition";
        return false;
    }

    wirehair_v2::PrecodeParams rejected_params;
    wirehair_v2::PacketRowConfig rejected_config;
    wirehair_v2::test::MakeRawArchitectureConfiguration(
        32u, rejected_params, rejected_config);
    wirehair_v2::PrecodeSystem rejected_system;
    if (!wirehair_v2::BuildPrecodeSystem(
            rejected_params, rejected_system))
    {
        error = "injected Schur-gate precode construction failed";
        return false;
    }
    for (std::vector<uint32_t>& row :
         rejected_system.DenseBasisRowColumns)
    {
        row.clear();
    }
    Partition rejected_partition;
    std::string rejected_error;
    if (BuildPartition(
            rejected_system, NULL, rejected_partition, rejected_error) !=
            EvaluationOutcome::Reject ||
        rejected_partition.SchurRank != 0u ||
        !rejected_partition.Qd.empty() ||
        rejected_error != "Schur construction did not find 12 pivots")
    {
        error = "injected Schur construction-gate mapping failed";
        return false;
    }

    Partition identity_partition;
    identity_partition.F.resize(257u);
    for (uint32_t i = 0u; i < identity_partition.F.size(); ++i) {
        identity_partition.F[i] = i;
    }
    Topology first;
    Topology second;
    if (BuildTopology(
            identity_partition,
            UINT64_C(0x123456789abcdef0), UINT32_C(0x76543210),
            NULL, first, error) != EvaluationOutcome::Pass ||
        BuildTopology(
            identity_partition,
            UINT64_C(0x123456789abcdef0), UINT32_C(0x76543210),
            NULL, second, error) != EvaluationOutcome::Pass)
    {
        return false;
    }
    std::string gate_error;
    if (EvaluateTopologyTerminalGate(first, gate_error) !=
        EvaluationOutcome::Pass)
    {
        error = "bounded topology unexpectedly failed terminal gate";
        return false;
    }
    Deadline expired_topology(0u);
    Topology timed_out_topology;
    std::string timeout_error;
    if (BuildTopology(
            identity_partition,
            UINT64_C(0x123456789abcdef0), UINT32_C(0x76543210),
            &expired_topology, timed_out_topology, timeout_error) !=
            EvaluationOutcome::Deadline ||
        timeout_error != "deadline")
    {
        error = "injected topology deadline mapping failed";
        return false;
    }
    Topology injected_reject = first;
    injected_reject.PendingTerminal = 65u;
    injected_reject.PendingMaximum = std::max<uint32_t>(
        injected_reject.PendingMaximum, injected_reject.PendingTerminal);
    gate_error.clear();
    if (EvaluateTopologyTerminalGate(injected_reject, gate_error) !=
            EvaluationOutcome::Reject ||
        injected_reject.PendingTerminal != 65u ||
        gate_error != "terminal pending count exceeds 64")
    {
        error = "injected topology terminal-gate differential failed";
        return false;
    }
    Partition invalid_topology_partition;
    invalid_topology_partition.F.push_back(0u);
    Topology invalid_topology;
    if (BuildTopology(
            invalid_topology_partition, 0u, 0u, NULL,
            invalid_topology, gate_error) != EvaluationOutcome::Invalid)
    {
        error = "topology invariant outcome mapping failed";
        return false;
    }
    return SelfCheck(
        DigestTopology(first) == DigestTopology(second) &&
            first.ReferenceCount == first.References.size() &&
            first.PendingTerminal <= 64u,
        "bounded topology differential failed", error);
}

bool SelfTestGoldenTopologyAndWork(std::string& error)
{
    Partition identity_partition;
    identity_partition.F.resize(17u);
    for (uint32_t column = 0u; column < identity_partition.F.size(); ++column) {
        identity_partition.F[column] = column;
    }
    static const uint32_t expected_free[] = {
        13u, 2u, 3u, 16u, 6u, 11u, 10u, 8u, 5u,
        15u, 0u, 12u, 4u, 9u, 1u, 14u, 7u
    };
    static const uint32_t expected_row[] = {
        1u, 15u, 0u, 16u, 11u, 8u, 5u, 2u, 4u,
        14u, 13u, 3u, 9u, 6u, 10u, 12u, 7u
    };
    static const uint32_t expected_degree[] = {
        10u, 2u, 13u, 9u, 0u, 7u, 5u, 8u, 12u,
        4u, 16u, 6u, 1u, 14u, 15u, 3u, 11u
    };
    const std::vector<uint32_t> expected_free_deck(
        expected_free, expected_free + 17u);
    const std::vector<uint32_t> expected_row_deck(
        expected_row, expected_row + 17u);
    const std::vector<uint32_t> expected_degree_deck(
        expected_degree, expected_degree + 17u);
    Topology small_topology;
    if (BuildTopology(
            identity_partition,
            kExpectedPrecodeSeed,
            kExpectedPacketSeed,
            NULL,
            small_topology,
            error) != EvaluationOutcome::Pass ||
        small_topology.FreeDeck != expected_free_deck ||
        small_topology.RowDeck != expected_row_deck ||
        small_topology.DegreeDeck != expected_degree_deck ||
        DigestTopology(small_topology) !=
            "879a7c66b8cfbe1464a6cf9429a2fc9ce3ad85a0b29a6b437fe4abce005c14d7" ||
        small_topology.ReferenceCount != 45u ||
        small_topology.FifoReferenceCount != 43u ||
        small_topology.ExcessReferenceCount != 2u ||
        small_topology.PendingMaximum != 11u ||
        small_topology.PendingTerminal != 10u ||
        small_topology.PendingNeedTerminal != 25u)
    {
        if (error.empty()) {
            error = "K17 independent topology golden mismatch";
        }
        return false;
    }

    static const uint32_t work_k = 2048u;
    wirehair_v2::PrecodeParams params;
    wirehair_v2::PacketRowConfig config;
    wirehair_v2::test::MakeRawArchitectureConfiguration(
        work_k, params, config);
    params = wirehair_v2::PrecodeParamsForAttempt(params, 0u);
    config = wirehair_v2::PacketConfigForAttempt(config, 0u);
    if (params.Seed != kExpectedPrecodeSeed ||
        config.PeelSeed != kExpectedPacketSeed ||
        config.MixCount != wirehair_v2::kCertifiedPacketMixCount ||
        params.Staircase != 54u || params.DenseRows != kDenseRows ||
        params.HeavyRows != kHeavyRows || params.SourceHits != 2u)
    {
        error = "K2048 independent parameter golden mismatch";
        return false;
    }
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system) ||
        !wirehair_v2::ValidatePrecodeSystem(system))
    {
        error = "K2048 golden precode construction failed";
        return false;
    }
    Partition partition;
    if (BuildPartition(system, NULL, partition, error) !=
            EvaluationOutcome::Pass ||
        partition.SchurRank != kDenseRows ||
        partition.FullQRank != 78u ||
        partition.QSourceCount != 11u)
    {
        if (error.empty()) {
            error = "K2048 independent partition golden mismatch";
        }
        return false;
    }
    Topology topology;
    if (BuildTopology(
            partition,
            params.Seed,
            config.PeelSeed,
            NULL,
            topology,
            error) != EvaluationOutcome::Pass ||
        EvaluateTopologyTerminalGate(topology, error) !=
            EvaluationOutcome::Pass ||
        topology.ReferenceCount != 8674u ||
        topology.FifoReferenceCount != 8025u ||
        topology.ExcessReferenceCount != 649u ||
        topology.PendingMaximum != 87u ||
        topology.PendingTerminal != 59u ||
        topology.PendingNeedTerminal != 167u)
    {
        if (error.empty()) {
            error = "K2048 independent topology golden mismatch";
        }
        return false;
    }
    WorkReceipt work;
    if (ComputeWork(system, partition, topology, work, error) !=
            EvaluationOutcome::Pass ||
        work.Free != 2048u ||
        work.Triangular != 8674u ||
        work.Staircase != 4127u ||
        work.DenseKnown != 1063u ||
        work.SchurSolveXors != 6u ||
        work.SchurSolveCopies != 12u ||
        work.SchurSolve != 18u ||
        work.StaircasePatch != 146u ||
        work.HeavyKnownXors != 2114u ||
        work.HeavyKnownMulAdds != 2928u ||
        work.HeavyKnown != 5042u ||
        work.HeavySolveXors != 0u ||
        work.HeavySolveMulAdds != 132u ||
        work.HeavySolveDivides != 12u ||
        work.HeavySolveCopies != 12u ||
        work.HeavySolve != 156u ||
        work.Emit != 2048u ||
        work.XorCopyTotal != 20250u ||
        work.Gf256Total != 3072u ||
        work.Total != 23322u ||
        work.SourceOwnership != 2048u ||
        work.XorCopyTotalWithOwnedSource != 22298u ||
        work.TotalWithOwnedSource != 25370u ||
        work.TotalWithOwnedSource >= static_cast<uint64_t>(13u) * work_k)
    {
        if (error.empty()) {
            error = "K2048 independent work golden mismatch";
        }
        return false;
    }
    Topology invalid_work_topology = topology;
    invalid_work_topology.ReferenceOffsets.pop_back();
    WorkReceipt invalid_work;
    std::string invalid_work_error;
    if (ComputeWork(
            system, partition, invalid_work_topology,
            invalid_work, invalid_work_error) != EvaluationOutcome::Invalid ||
        invalid_work_error != "work accounting input shape is invalid")
    {
        error = "injected work-invariant outcome mapping failed";
        return false;
    }
    return true;
}

bool SelfTestPeel(std::string& error)
{
    // This asymmetric chain is deliberately supplied with a valid triangular
    // solve order.  Forward projection makes row {0,2} dependent; reversing
    // the traversal asks for column 1 before its projection is available.
    ExplicitRows chain_rows;
    chain_rows.Append(std::vector<uint32_t>{ 0u, 1u });
    chain_rows.Append(std::vector<uint32_t>{ 1u, 2u });
    chain_rows.Append(std::vector<uint32_t>{ 0u, 2u });
    PeelOutcome chain;
    chain.SolveRow.assign(3u, UINT32_MAX);
    chain.SolveRow[1u] = 0u;
    chain.SolveRow[2u] = 1u;
    chain.PeelOrder.push_back(1u);
    chain.PeelOrder.push_back(2u);
    chain.InactiveOrder.push_back(0u);
    chain.UsedRows.push_back(1u);
    chain.UsedRows.push_back(1u);
    chain.UsedRows.push_back(0u);
    chain.Complete = true;
    ResidualReceipt chain_residual;
    if (!AnalyzeBinaryResidual(
            3u, chain_rows, chain, NULL, chain_residual, error) ||
        chain_residual.Rows != 1u || chain_residual.BinaryRank != 0u ||
        chain_residual.Deficiency != 1u ||
        chain.PeelOrder.size() + chain_residual.BinaryRank !=
            BruteForceGf2Rank(3u, chain_rows))
    {
        if (error.empty()) {
            error = "chronological projection chain differential failed";
        }
        return false;
    }
    PeelOutcome reverse_chain = chain;
    std::reverse(
        reverse_chain.PeelOrder.begin(), reverse_chain.PeelOrder.end());
    ResidualReceipt rejected_projection;
    std::string projected_error;
    if (AnalyzeBinaryResidual(
            3u, chain_rows, reverse_chain, NULL,
            rejected_projection, projected_error) ||
        projected_error != "peeled dependency projection is unavailable")
    {
        error = "projected-state invariant guard failed";
        return false;
    }

    ExplicitRows full_peel_rows;
    full_peel_rows.Append(std::vector<uint32_t>(1u, 0u));
    PeelOutcome full_peel;
    if (!PeelProductionOrder(
            1u, full_peel_rows, 0u, NULL, full_peel, error))
    {
        return false;
    }
    ResidualReceipt full_peel_residual;
    if (!AnalyzeBinaryResidual(
            1u, full_peel_rows, full_peel, NULL,
            full_peel_residual, error) ||
        full_peel_residual.Rows != 0u ||
        full_peel_residual.BinaryRank != 0u ||
        full_peel_residual.Deficiency != 0u)
    {
        if (error.empty()) {
            error = "zero-width residual path failed";
        }
        return false;
    }
    Deadline expired(0u);
    ResidualReceipt expired_residual;
    std::string deadline_error;
    if (AnalyzeBinaryResidual(
            1u, full_peel_rows, full_peel, &expired,
            expired_residual, deadline_error) ||
        deadline_error != "deadline")
    {
        error = "zero-width residual deadline guard failed";
        return false;
    }

    ReferencePcg rng;
    rng.Seed(UINT64_C(0x7065656c2d6f7263), UINT64_C(0x5678));
    for (uint32_t trial = 0u; trial < 80u; ++trial)
    {
        const uint32_t columns = 2u + rng.Next() % 31u;
        const uint32_t row_count = 1u + rng.Next() % (columns + 20u);
        ExplicitRows rows;
        for (uint32_t row_index = 0u; row_index < row_count; ++row_index)
        {
            const uint32_t degree =
                1u + rng.Next() % std::min<uint32_t>(columns, 8u);
            std::vector<uint32_t> row;
            while (row.size() < degree)
            {
                const uint32_t column = rng.Next() % columns;
                if (std::find(row.begin(), row.end(), column) == row.end()) {
                    row.push_back(column);
                }
            }
            rows.Append(row);
        }
        PeelOutcome fast;
        PeelOutcome slow;
        if (!PeelProductionOrder(
                columns, rows, 0u, NULL, fast, error) ||
            !PeelSlowReference(columns, rows, slow, error) ||
            !SamePeel(fast, slow))
        {
            if (error.empty()) {
                error = "fast versus slow peel-order mismatch";
            }
            return false;
        }
        ResidualReceipt residual;
        if (!AnalyzeBinaryResidual(
                columns, rows, fast, NULL, residual, error))
        {
            return false;
        }
        const uint32_t brute_rank = BruteForceGf2Rank(columns, rows);
        if (fast.PeelOrder.size() + residual.BinaryRank != brute_rank)
        {
            error = "peeled plus residual rank identity failed";
            return false;
        }
    }
    return true;
}

bool SelfTestProductionDifferential(std::string& error)
{
    wirehair_v2::ResetBinaryPeelOracleComparisonsForTesting();
    wirehair_v2::SetBinaryPeelOracleForTesting(true);
    uint32_t comparisons = 0u;
    static const uint32_t cases[] = { 8u, 32u, 100u };
    for (uint32_t K : cases)
    {
        const wirehair_v2::SeedProfile profile =
            wirehair_v2::SelectSeedProfile(K, 1u);
        wirehair_v2::MessagePrecodeEncoderOptions options;
        wirehair_v2::PrecodeParams params;
        wirehair_v2::PacketRowConfig config;
        if (!wirehair_v2::ResolveMessagePrecodeConfiguration(
                profile, options, params, config))
        {
            wirehair_v2::SetBinaryPeelOracleForTesting(false);
            error = "production differential configuration failed";
            return false;
        }
        wirehair_v2::PrecodeSystem system;
        if (!wirehair_v2::BuildPrecodeSystem(params, system))
        {
            wirehair_v2::SetBinaryPeelOracleForTesting(false);
            error = "production differential precode failed";
            return false;
        }
        const uint32_t P = params.Staircase + params.DenseRows +
            params.HeavyRows;
        const uint32_t L = K + P;
        ExplicitRows rows;
        for (const std::vector<uint32_t>& row : system.StaircaseRows) {
            rows.Append(row);
        }
        for (const std::vector<uint32_t>& row :
             system.DenseBasisRowColumns)
        {
            rows.Append(row);
        }
        const uint32_t packet_count = K + 4u;
        std::vector<uint8_t> zeros(packet_count, 0u);
        std::vector<wirehair_v2::SolvePacket> packets(packet_count);
        for (uint32_t row = 0u; row < packet_count; ++row)
        {
            packets[row].BlockId = row;
            packets[row].Data = &zeros[row];
            const std::vector<uint32_t> columns =
                wirehair_v2::GeneratePacketMatrixRow(K, P, row, config);
            if (columns.empty())
            {
                wirehair_v2::SetBinaryPeelOracleForTesting(false);
                error = "production differential row generation failed";
                return false;
            }
            rows.Append(columns);
        }
        PeelOutcome peel;
        ResidualReceipt residual;
        if (!PeelProductionOrder(L, rows, 0u, NULL, peel, error) ||
            !AnalyzeBinaryResidual(
                L, rows, peel, NULL, residual, error))
        {
            wirehair_v2::SetBinaryPeelOracleForTesting(false);
            return false;
        }
        std::vector<uint8_t> output;
        wirehair_v2::PrecodeSolveStats stats;
        const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
            system, config, packets, 1u, output, &stats);
        if ((result != Wirehair_Success && result != Wirehair_NeedMore) ||
            stats.PeeledColumns != peel.PeelOrder.size() ||
            stats.InactivatedColumns != peel.InactiveOrder.size() ||
            stats.BinaryResidualRank != residual.BinaryRank ||
            stats.ResidualRows != residual.Rows + params.HeavyRows)
        {
            wirehair_v2::SetBinaryPeelOracleForTesting(false);
            error = "production solve-stat peel differential failed";
            return false;
        }
        ++comparisons;
    }
    wirehair_v2::SetBinaryPeelOracleForTesting(false);
    return SelfCheck(
        wirehair_v2::BinaryPeelOracleComparisonsForTesting() == comparisons,
        "production internal binary-peel oracle did not run", error);
}

bool SelfTestFrozenCellMetadata(std::string& error)
{
    const std::vector<FrozenRecoveryCell> cells =
        EnumerateDevelopmentRecoveryCells();
    if (cells.size() <= 59u)
    {
        error = "frozen development domain lacks ordinal 59";
        return false;
    }
    const FrozenRecoveryCell& cell = cells[59u];
    return SelfCheck(
        cell.ordinal == 59u && cell.K == kK &&
            cell.block_bytes == kBlockBytes &&
            cell.loss_ppm == kLossPpm &&
            cell.schedule == FrozenSchedule::Burst &&
            cell.trial == 0u && cell.base_seed_attempt == kAttempt &&
            cell.loss_seed == kLossRoot &&
            RecoveryCellSha256(cell) == kExpectedCellSha &&
            !CanonicalRecoveryCellJson(cell).empty(),
        "frozen ordinal-59 metadata mismatch", error);
}

bool SelfTestReceiptExitSemantics(std::string& error);

bool RunSelfTests(std::string& error)
{
    if (gf256_init() != 0)
    {
        error = "GF256 initialization failed";
        return false;
    }
    return SelfTestPcgAndDeck(error) &&
        SelfTestGf256Rank(error) &&
        SelfTestRawSeedConfiguration(error) &&
        SelfTestPartitionAndTopology(error) &&
        SelfTestGoldenTopologyAndWork(error) &&
        SelfTestPeel(error) &&
        SelfTestProductionDifferential(error) &&
        SelfTestFrozenCellMetadata(error) &&
        SelfTestReceiptExitSemantics(error);
}

std::string JsonEscape(const std::string& value)
{
    static const char hex[] = "0123456789abcdef";
    std::string output;
    output.push_back('"');
    for (unsigned char ch : value)
    {
        switch (ch)
        {
        case '"': output += "\\\""; break;
        case '\\': output += "\\\\"; break;
        case '\b': output += "\\b"; break;
        case '\f': output += "\\f"; break;
        case '\n': output += "\\n"; break;
        case '\r': output += "\\r"; break;
        case '\t': output += "\\t"; break;
        default:
            if (ch < 0x20u || ch >= 0x7fu)
            {
                output += "\\u00";
                output.push_back(hex[ch >> 4]);
                output.push_back(hex[ch & 15u]);
            }
            else {
                output.push_back(static_cast<char>(ch));
            }
            break;
        }
    }
    output.push_back('"');
    return output;
}

struct RunReceipt
{
    std::string Status = "invalid";
    std::string RejectStage = "oracle";
    std::string Reason = "not_started";
    std::string License = "none";
    std::string CellJson;
    std::string CellSha;
    std::string TraceSha;
    std::string PrefixSha;
    std::string PrecodeSha;
    std::string PartitionSha;
    std::string FreeDeckSha;
    std::string RowDeckSha;
    std::string DegreeDeckSha;
    std::string TopologySha;
    std::string SystematicSha;
    std::string RepairSha;
    std::string RepairOracleSha;
    std::string PacketSha;
    std::string UncoveredAllSha;
    std::string UncoveredFSha;
    uint64_t TraceSeed = 0u;
    uint64_t BurstSeed = 0u;
    uint64_t CandidateLimit = 0u;
    uint64_t AttemptedCandidates = 0u;
    uint64_t PrecodeSeed = 0u;
    uint32_t PacketSeed = 0u;
    uint32_t S = 0u;
    uint32_t P = 0u;
    uint32_t L = 0u;
    uint32_t N1 = 0u;
    uint32_t SchurRank = 0u;
    uint32_t QdCount = 0u;
    uint32_t FullQRank = 0u;
    uint32_t QSourceCount = 0u;
    uint64_t TriangularReferences = 0u;
    uint64_t FifoReferences = 0u;
    uint64_t ExcessReferences = 0u;
    uint32_t PendingMaximum = 0u;
    uint32_t PendingTerminal = 0u;
    uint32_t PendingNeedTerminal = 0u;
    uint32_t SystematicPackets = 0u;
    uint32_t RepairPackets = 0u;
    uint32_t UncoveredAll = 0u;
    uint32_t UncoveredAllRank = 0u;
    bool UncoveredAllRankComputed = false;
    uint32_t UncoveredF = 0u;
    uint32_t UncoveredFRank = 0u;
    bool UncoveredFRankComputed = false;
    WorkReceipt Work;
    uint32_t Peeled = 0u;
    uint32_t Inactivated = 0u;
    uint32_t ResidualRows = 0u;
    uint32_t BinaryResidualRank = 0u;
    uint32_t Deficiency = 0u;
    uint64_t ElapsedMilliseconds = 0u;
};

void RecordPartitionReceipt(
    RunReceipt& receipt,
    const Partition& partition)
{
    receipt.SchurRank = partition.SchurRank;
    receipt.QdCount = static_cast<uint32_t>(partition.Qd.size());
    receipt.FullQRank = partition.FullQRank;
    receipt.QSourceCount = partition.QSourceCount;
    receipt.PartitionSha = DigestPartition(partition);
}

void RecordTopologyReceipt(
    RunReceipt& receipt,
    const Topology& topology,
    bool complete)
{
    // Preserve allocation-free numeric evidence before any digest can throw.
    receipt.TriangularReferences = topology.ReferenceCount;
    receipt.FifoReferences = topology.FifoReferenceCount;
    receipt.ExcessReferences = topology.ExcessReferenceCount;
    receipt.PendingMaximum = topology.PendingMaximum;
    receipt.PendingTerminal = topology.PendingTerminal;
    receipt.PendingNeedTerminal = topology.PendingNeedTerminal;
    if (complete)
    {
        receipt.FreeDeckSha = DigestVector("tri-bmq4-free-deck-v1",
            topology.FreeDeck);
        receipt.RowDeckSha = DigestVector("tri-bmq4-row-deck-v1",
            topology.RowDeck);
        receipt.DegreeDeckSha = DigestVector("tri-bmq4-degree-deck-v1",
            topology.DegreeDeck);
        receipt.TopologySha = DigestTopology(topology);
    }
}

int ReceiptExitCode(const RunReceipt& receipt)
{
    if (receipt.Status == "pass") {
        return kExitPass;
    }
    if (receipt.Status == "reject") {
        return kExitReject;
    }
    return kExitInvalid;
}

std::string ReceiptJson(const RunReceipt& r)
{
    std::ostringstream out;
    out << '{'
        << "\"schema\":" << JsonEscape(kSchema)
        << ",\"candidate\":" << JsonEscape(kCandidate)
        << ",\"bead\":" << JsonEscape(kBead)
        << ",\"design_base_commit\":" << JsonEscape(kDesignCommit)
        << ",\"source_commit\":"
        << JsonEscape(WIREHAIR_WH2_SOURCE_GIT_COMMIT)
        << ",\"status\":" << JsonEscape(r.Status)
        << ",\"exit_code\":" << ReceiptExitCode(r)
        << ",\"reject_stage\":" << JsonEscape(r.RejectStage)
        << ",\"reason\":" << JsonEscape(r.Reason)
        << ",\"license\":" << JsonEscape(r.License)
        << ",\"cell_ordinal\":59"
        << ",\"cell_json\":" << JsonEscape(r.CellJson)
        << ",\"cell_sha256\":" << JsonEscape(r.CellSha)
        << ",\"K\":" << kK
        << ",\"block_bytes\":" << kBlockBytes
        << ",\"loss_ppm\":" << kLossPpm
        << ",\"schedule\":\"burst\""
        << ",\"overhead\":0"
        << ",\"trial\":0"
        << ",\"base_seed_attempt\":0"
        << ",\"loss_root\":" << JsonEscape(Hex64(kLossRoot))
        << ",\"trace_seed\":" << JsonEscape(Hex64(r.TraceSeed))
        << ",\"burst_rng_seed\":" << JsonEscape(Hex64(r.BurstSeed))
        << ",\"candidate_limit\":" << r.CandidateLimit
        << ",\"attempted_candidates\":" << r.AttemptedCandidates
        << ",\"trace_k_plus_4_sha256\":" << JsonEscape(r.TraceSha)
        << ",\"trace_oh0_sha256\":" << JsonEscape(r.PrefixSha)
        << ",\"construction\":\"uniform-raw-v1\""
        << ",\"seed_schedule_sha256\":" << JsonEscape(kSeedScheduleSha)
        << ",\"precode_seed\":" << JsonEscape(Hex64(r.PrecodeSeed))
        << ",\"packet_seed\":" << JsonEscape(Hex32(r.PacketSeed))
        << ",\"precode_contract_version\":"
        << wirehair_v2::PrecodeContractVersion()
        << ",\"packet_contract_version\":"
        << wirehair_v2::kPacketRowContractVersion
        << ",\"mix_count\":3"
        << ",\"S\":" << r.S
        << ",\"D2\":12"
        << ",\"H\":12"
        << ",\"P\":" << r.P
        << ",\"L\":" << r.L
        << ",\"N1\":" << r.N1
        << ",\"dense_identity_corner\":false"
        << ",\"dense_anchors\":\"disabled\""
        << ",\"heavy_family\":\"periodic-cauchy\""
        << ",\"precode_sha256\":" << JsonEscape(r.PrecodeSha)
        << ",\"partition_sha256\":" << JsonEscape(r.PartitionSha)
        << ",\"free_deck_sha256\":" << JsonEscape(r.FreeDeckSha)
        << ",\"row_deck_sha256\":" << JsonEscape(r.RowDeckSha)
        << ",\"degree_deck_sha256\":" << JsonEscape(r.DegreeDeckSha)
        << ",\"topology_sha256\":" << JsonEscape(r.TopologySha)
        << ",\"systematic_rows_sha256\":" << JsonEscape(r.SystematicSha)
        << ",\"repair_rows_sha256\":" << JsonEscape(r.RepairSha)
        << ",\"repair_oracle_sha256\":" << JsonEscape(r.RepairOracleSha)
        << ",\"delivered_rows_sha256\":" << JsonEscape(r.PacketSha)
        << ",\"uncovered_all_sha256\":" << JsonEscape(r.UncoveredAllSha)
        << ",\"uncovered_f_sha256\":" << JsonEscape(r.UncoveredFSha)
        << ",\"schur_rank\":" << r.SchurRank
        << ",\"qd_count\":" << r.QdCount
        << ",\"full_q_rank\":" << r.FullQRank
        << ",\"q_source_count\":" << r.QSourceCount
        << ",\"triangular_references\":" << r.TriangularReferences
        << ",\"fifo_references\":" << r.FifoReferences
        << ",\"excess_references\":" << r.ExcessReferences
        << ",\"pending_maximum\":" << r.PendingMaximum
        << ",\"pending_terminal\":" << r.PendingTerminal
        << ",\"pending_need_terminal\":" << r.PendingNeedTerminal
        << ",\"systematic_packets\":" << r.SystematicPackets
        << ",\"repair_packets\":" << r.RepairPackets
        << ",\"uncovered_all\":" << r.UncoveredAll
        << ",\"uncovered_all_rank\":" << r.UncoveredAllRank
        << ",\"uncovered_all_rank_computed\":"
        << (r.UncoveredAllRankComputed ? "true" : "false")
        << ",\"uncovered_f\":" << r.UncoveredF
        << ",\"uncovered_f_rank\":" << r.UncoveredFRank
        << ",\"uncovered_f_rank_computed\":"
        << (r.UncoveredFRankComputed ? "true" : "false")
        << ",\"uncovered_f_scope\":\"diagnostic-redundant-oracle\""
        << ",\"work_free\":" << r.Work.Free
        << ",\"work_triangular\":" << r.Work.Triangular
        << ",\"work_staircase\":" << r.Work.Staircase
        << ",\"work_dense_known\":" << r.Work.DenseKnown
        << ",\"work_schur_solve_xors\":" << r.Work.SchurSolveXors
        << ",\"work_schur_solve_copies\":" << r.Work.SchurSolveCopies
        << ",\"work_schur_solve\":" << r.Work.SchurSolve
        << ",\"work_staircase_patch\":" << r.Work.StaircasePatch
        << ",\"work_heavy_known_xors\":" << r.Work.HeavyKnownXors
        << ",\"work_heavy_known_muladds\":"
        << r.Work.HeavyKnownMulAdds
        << ",\"work_heavy_known\":" << r.Work.HeavyKnown
        << ",\"work_heavy_solve_xors\":" << r.Work.HeavySolveXors
        << ",\"work_heavy_solve_muladds\":"
        << r.Work.HeavySolveMulAdds
        << ",\"work_heavy_solve_divides\":"
        << r.Work.HeavySolveDivides
        << ",\"work_heavy_solve_copies\":"
        << r.Work.HeavySolveCopies
        << ",\"work_heavy_solve\":" << r.Work.HeavySolve
        << ",\"work_emit\":" << r.Work.Emit
        << ",\"work_xor_copy_total\":" << r.Work.XorCopyTotal
        << ",\"work_gf256_total\":" << r.Work.Gf256Total
        << ",\"work_total\":" << r.Work.Total
        << ",\"work_source_ownership\":" << r.Work.SourceOwnership
        << ",\"work_xor_copy_total_with_owned_source\":"
        << r.Work.XorCopyTotalWithOwnedSource
        << ",\"work_total_with_owned_source\":"
        << r.Work.TotalWithOwnedSource
        << ",\"work_limit\":" << static_cast<uint64_t>(13u) * kK
        << ",\"core_limit\":" << kCoreLimit
        << ",\"peeled\":" << r.Peeled
        << ",\"inactivated\":" << r.Inactivated
        << ",\"residual_rows\":" << r.ResidualRows
        << ",\"binary_residual_rank\":" << r.BinaryResidualRank
        << ",\"deficiency\":" << r.Deficiency
        << ",\"deadline_seconds\":" << kDeadlineSeconds
        << ",\"elapsed_milliseconds\":" << r.ElapsedMilliseconds
        << '}';
    if (!out) {
        throw std::runtime_error("receipt JSON formatting failed");
    }
    return out.str();
}

int EmissionExitCode(int receipt_exit_code, bool emitted)
{
    return emitted ? receipt_exit_code : kExitInvalid;
}

int EmitSerializedReceipt(
    const std::string& json,
    int receipt_exit_code)
{
    bool emitted =
        std::fwrite(json.data(), 1u, json.size(), stdout) == json.size();
    if (emitted) {
        emitted = std::fputc('\n', stdout) != EOF;
    }
    if (std::fflush(stdout) != 0 || std::ferror(stdout) != 0) {
        emitted = false;
    }
    if (!emitted) {
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: failed to emit complete receipt\n");
    }
    return EmissionExitCode(receipt_exit_code, emitted);
}

bool SelfTestReceiptExitSemantics(std::string& error)
{
    const EvaluationMapping pass = MapEvaluationOutcome(
        EvaluationOutcome::Pass, "ignored");
    const EvaluationMapping reject = MapEvaluationOutcome(
        EvaluationOutcome::Reject, "construction_gate");
    const EvaluationMapping invalid = MapEvaluationOutcome(
        EvaluationOutcome::Invalid, "ignored");
    const EvaluationMapping deadline = MapEvaluationOutcome(
        EvaluationOutcome::Deadline, "ignored");
    if (std::strcmp(pass.Status, "pass") != 0 ||
        std::strcmp(pass.Stage, "none") != 0 || pass.ExitCode != kExitPass ||
        std::strcmp(reject.Status, "reject") != 0 ||
        std::strcmp(reject.Stage, "construction_gate") != 0 ||
        reject.ExitCode != kExitReject ||
        std::strcmp(invalid.Status, "invalid") != 0 ||
        std::strcmp(invalid.Stage, "oracle") != 0 ||
        invalid.ExitCode != kExitInvalid ||
        std::strcmp(deadline.Status, "invalid") != 0 ||
        std::strcmp(deadline.Stage, "deadline") != 0 ||
        deadline.ExitCode != kExitInvalid ||
        EmissionExitCode(kExitPass, false) != kExitInvalid ||
        EmissionExitCode(kExitReject, false) != kExitInvalid)
    {
        error = "typed outcome or emission-failure mapping mismatch";
        return false;
    }
    std::string hostile_json_input;
    hostile_json_input.push_back(static_cast<char>(0x7f));
    hostile_json_input.push_back(static_cast<char>(0x80));
    hostile_json_input.push_back(static_cast<char>(0xff));
    if (JsonEscape(hostile_json_input) != "\"\\u007f\\u0080\\u00ff\"")
    {
        error = "canonical JSON byte escaping mismatch";
        return false;
    }
    RunReceipt receipt;
    if (ReceiptExitCode(receipt) != kExitInvalid) {
        error = "default receipt exit semantics mismatch";
        return false;
    }
    receipt.Status = "reject";
    if (ReceiptExitCode(receipt) != kExitReject ||
        ReceiptJson(receipt).find("\"exit_code\":2") == std::string::npos)
    {
        error = "reject receipt exit semantics mismatch";
        return false;
    }
    Partition partial_partition;
    partial_partition.Qd.assign(7u, 0u);
    partial_partition.QdMasks.assign(7u, 0u);
    partial_partition.SchurRank = 7u;
    RecordPartitionReceipt(receipt, partial_partition);
    Topology terminal_reject;
    terminal_reject.ReferenceCount = 123u;
    terminal_reject.FifoReferenceCount = 100u;
    terminal_reject.ExcessReferenceCount = 23u;
    terminal_reject.PendingMaximum = 71u;
    terminal_reject.PendingTerminal = 65u;
    terminal_reject.PendingNeedTerminal = 129u;
    RecordTopologyReceipt(receipt, terminal_reject, true);
    std::string gate_error;
    if (receipt.SchurRank != 7u || receipt.QdCount != 7u ||
        receipt.PartitionSha.empty() || receipt.TriangularReferences != 123u ||
        receipt.PendingMaximum != 71u || receipt.PendingTerminal != 65u ||
        receipt.TopologySha.empty() ||
        EvaluateTopologyTerminalGate(terminal_reject, gate_error) !=
            EvaluationOutcome::Reject)
    {
        error = "construction-reject receipt ordering mismatch";
        return false;
    }
    receipt.Status = "pass";
    return SelfCheck(
        ReceiptExitCode(receipt) == kExitPass &&
            ReceiptJson(receipt).find("\"exit_code\":0") !=
                std::string::npos,
        "pass receipt exit semantics mismatch", error);
}

bool PrepareFrozenTrace(
    RunReceipt& receipt,
    std::vector<uint32_t>& prefix,
    std::string& error)
{
    const std::vector<FrozenRecoveryCell> cells =
        EnumerateDevelopmentRecoveryCells();
    if (cells.size() <= 59u)
    {
        error = "frozen development domain lacks ordinal 59";
        return false;
    }
    const FrozenRecoveryCell& cell = cells[59u];
    receipt.CellJson = CanonicalRecoveryCellJson(cell);
    receipt.CellSha = RecoveryCellSha256(cell);
    if (cell.ordinal != 59u || cell.K != kK ||
        cell.block_bytes != kBlockBytes || cell.loss_ppm != kLossPpm ||
        cell.schedule != FrozenSchedule::Burst || cell.trial != 0u ||
        cell.base_seed_attempt != kAttempt || cell.loss_seed != kLossRoot ||
        receipt.CellJson.empty() || receipt.CellSha != kExpectedCellSha)
    {
        error = "frozen ordinal-59 cell receipt mismatch";
        return false;
    }
    FrozenPacketTrace trace;
    if (GenerateFrozenPacketTrace(cell, trace) != FrozenTraceStatus::Complete ||
        !CopyNestedPrefix(trace, 0u, prefix))
    {
        error = "frozen ordinal-59 trace generation failed";
        return false;
    }
    receipt.TraceSeed = trace.trace_seed;
    receipt.BurstSeed = trace.trace_seed ^ UINT64_C(0x10fade);
    receipt.CandidateLimit = trace.candidate_limit;
    receipt.AttemptedCandidates = trace.attempted_candidates;
    receipt.TraceSha = trace.trace_sha256;
    receipt.PrefixSha = PacketIdsSha256(prefix);
    if (prefix.size() != kK ||
        receipt.TraceSeed != kExpectedTraceSeed ||
        receipt.BurstSeed != kExpectedBurstSeed ||
        receipt.CandidateLimit != kExpectedCandidateLimit ||
        receipt.AttemptedCandidates != kExpectedAttemptedCandidates ||
        receipt.TraceSha != kExpectedTraceSha ||
        receipt.PrefixSha != kExpectedPrefixSha)
    {
        error = "frozen ordinal-59 trace receipt mismatch";
        return false;
    }
    return true;
}

void FinishReceipt(
    RunReceipt& receipt,
    const char* status,
    const char* stage,
    const std::string& reason,
    const std::chrono::steady_clock::time_point& start)
{
    receipt.Status = status;
    receipt.RejectStage = stage;
    receipt.Reason = reason;
    receipt.License = std::strcmp(status, "pass") == 0 ?
        "broader_bounded_structure_screen_only" : "none";
    receipt.ElapsedMilliseconds = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start).count());
}

struct RunContext
{
    RunContext()
        : Start(std::chrono::steady_clock::now())
        , Timer(kDeadlineSeconds)
    {
    }

    std::chrono::steady_clock::time_point Start;
    Deadline Timer;
    RunReceipt Receipt;
};

EvaluationOutcome InvalidOrDeadline(const std::string& error)
{
    return error == "deadline" ?
        EvaluationOutcome::Deadline : EvaluationOutcome::Invalid;
}

int EmitTerminalReceipt(RunContext& context)
{
    const auto expire_scientific_result = [&]() {
        const bool expired = context.Timer.Expired();
        if ((context.Receipt.Status == "pass" ||
             context.Receipt.Status == "reject") &&
            expired)
        {
            FinishReceipt(context.Receipt, "invalid", "deadline", "deadline",
                context.Start);
            return true;
        }
        return false;
    };
    expire_scientific_result();
    std::string json = ReceiptJson(context.Receipt);
    // Receipt construction itself can allocate and format many fields.  Check
    // once more after it completes so a pass/reject is never emitted after the
    // frozen wall-clock deadline.
    if (expire_scientific_result()) {
        json = ReceiptJson(context.Receipt);
    }
    return EmitSerializedReceipt(json, ReceiptExitCode(context.Receipt));
}

int CompleteRun(
    RunContext& context,
    EvaluationOutcome outcome,
    const char* reject_stage,
    const std::string& reason)
{
    const EvaluationMapping mapping =
        MapEvaluationOutcome(outcome, reject_stage);
    FinishReceipt(
        context.Receipt,
        mapping.Status,
        mapping.Stage,
        reason,
        context.Start);
    return EmitTerminalReceipt(context);
}

int RunFrozenCellImpl(RunContext& context)
{
    RunReceipt& receipt = context.Receipt;
    Deadline& deadline = context.Timer;
    std::string error;
    std::vector<uint32_t> packet_prefix;
    if (gf256_init() != 0)
    {
        return CompleteRun(context, EvaluationOutcome::Invalid, "oracle",
            "GF256 initialization failed");
    }
    if (!ValidateRawSeedLaw(error) || !ValidateContractVersions(error))
    {
        return CompleteRun(
            context, EvaluationOutcome::Invalid, "oracle", error);
    }
    if (!PrepareFrozenTrace(receipt, packet_prefix, error))
    {
        return CompleteRun(
            context, EvaluationOutcome::Invalid, "oracle", error);
    }
    if (deadline.Expired())
    {
        return CompleteRun(
            context, EvaluationOutcome::Deadline, "deadline", "deadline");
    }

    wirehair_v2::PrecodeParams params;
    wirehair_v2::PacketRowConfig config;
    wirehair_v2::test::MakeRawArchitectureConfiguration(kK, params, config);
    params = wirehair_v2::PrecodeParamsForAttempt(params, kAttempt);
    config = wirehair_v2::PacketConfigForAttempt(config, kAttempt);
    receipt.PrecodeSeed = params.Seed;
    receipt.PacketSeed = config.PeelSeed;
    receipt.S = params.Staircase;
    receipt.P = params.Staircase + params.DenseRows + params.HeavyRows;
    receipt.L = params.BlockCount + receipt.P;
    receipt.N1 = params.SourceHits;
    const wirehair_v2::PrecodeParams expected =
        wirehair_v2::MakeCertifiedParams(kK, kExpectedPrecodeSeed);
    if (!SameParams(params, expected) ||
        params.Staircase != 346u || params.DenseRows != kDenseRows ||
        params.HeavyRows != kHeavyRows || params.SourceHits != 3u ||
        params.Seed != kExpectedPrecodeSeed ||
        config.PeelSeed != kExpectedPacketSeed || config.MixCount != 3u ||
        receipt.P != 370u || receipt.L != 64370u)
    {
        return CompleteRun(context, EvaluationOutcome::Invalid, "oracle",
            "attempt-zero construction receipt mismatch");
    }
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system) ||
        !wirehair_v2::ValidatePrecodeSystem(system))
    {
        return CompleteRun(context, EvaluationOutcome::Invalid, "oracle",
            "certified precode construction failed");
    }
    receipt.PrecodeSha = DigestPrecode(system);

    Partition partition;
    const EvaluationOutcome partition_outcome =
        BuildPartition(system, &deadline, partition, error);
    RecordPartitionReceipt(receipt, partition);
    if (partition_outcome != EvaluationOutcome::Pass)
    {
        return CompleteRun(
            context, partition_outcome, "construction_gate", error);
    }

    Topology topology;
    const EvaluationOutcome topology_outcome = BuildTopology(
            partition, params.Seed, config.PeelSeed,
            &deadline, topology, error);
    RecordTopologyReceipt(
        receipt, topology, topology_outcome == EvaluationOutcome::Pass);
    if (topology_outcome != EvaluationOutcome::Pass)
    {
        return CompleteRun(
            context, topology_outcome, "construction_gate", error);
    }
    const EvaluationOutcome terminal_outcome =
        EvaluateTopologyTerminalGate(topology, error);
    if (terminal_outcome != EvaluationOutcome::Pass)
    {
        return CompleteRun(
            context, terminal_outcome, "construction_gate", error);
    }

    DeliveredRows delivered;
    if (!BuildDeliveredRows(
            system, config, packet_prefix, topology,
            &deadline, delivered, error))
    {
        return CompleteRun(
            context, InvalidOrDeadline(error), "oracle", error);
    }
    receipt.SystematicPackets = delivered.SystematicCount;
    receipt.RepairPackets = delivered.RepairCount;
    receipt.SystematicSha = delivered.SystematicDigest;
    receipt.RepairSha = delivered.RepairDigest;
    receipt.RepairOracleSha = delivered.RepairOracleDigest;
    receipt.PacketSha = delivered.PacketDigest;
    if (receipt.RepairSha != receipt.RepairOracleSha)
    {
        return CompleteRun(context, EvaluationOutcome::Invalid, "oracle",
            "repair row digest differs from head oracle");
    }

    std::vector<uint32_t> uncovered_all;
    std::vector<uint32_t> uncovered_f;
    uncovered_all.reserve(receipt.P);
    uncovered_f.reserve(receipt.P);
    for (uint32_t column = 0u; column < receipt.L; ++column)
    {
        if ((column & 1023u) == 0u && deadline.Expired())
        {
            return CompleteRun(context, EvaluationOutcome::Deadline,
                "deadline", "deadline");
        }
        if (delivered.PacketCovered[column]) {
            continue;
        }
        uncovered_all.push_back(column);
        if (!partition.IsQ[column]) {
            uncovered_f.push_back(column);
        }
    }
    receipt.UncoveredAll = static_cast<uint32_t>(uncovered_all.size());
    receipt.UncoveredF = static_cast<uint32_t>(uncovered_f.size());
    receipt.UncoveredAllSha = DigestVector(
        "tri-bmq4-uncovered-all-v1", uncovered_all);
    receipt.UncoveredFSha = DigestVector(
        "tri-bmq4-uncovered-f-v1", uncovered_f);

    if (uncovered_all.size() > receipt.P)
    {
        return CompleteRun(context, EvaluationOutcome::Reject,
            "coverage_rank", "uncovered all-global column count exceeds P");
    }
    if (!RankPrecodeRestriction(
            system, uncovered_all, &deadline,
            receipt.UncoveredAllRank, error))
    {
        return CompleteRun(
            context, InvalidOrDeadline(error), "oracle", error);
    }
    receipt.UncoveredAllRankComputed = true;
    if (receipt.UncoveredAllRank != receipt.UncoveredAll)
    {
        return CompleteRun(context, EvaluationOutcome::Reject,
            "coverage_rank",
            "precode restriction on all-global U is rank deficient");
    }
    if (!RankPrecodeRestriction(
            system, uncovered_f, &deadline,
            receipt.UncoveredFRank, error))
    {
        return CompleteRun(
            context, InvalidOrDeadline(error), "oracle", error);
    }
    receipt.UncoveredFRankComputed = true;
    if (receipt.UncoveredFRank != receipt.UncoveredF)
    {
        return CompleteRun(context, EvaluationOutcome::Invalid, "oracle",
            "restricted U_F rank disagrees with full-U oracle");
    }

    const EvaluationOutcome work_outcome =
        ComputeWork(system, partition, topology, receipt.Work, error);
    if (work_outcome != EvaluationOutcome::Pass)
    {
        return CompleteRun(context, work_outcome, "work", error);
    }

    PeelOutcome peel;
    if (!PeelProductionOrder(
            receipt.L, delivered.BinaryRows, kCoreLimit,
            &deadline, peel, error))
    {
        return CompleteRun(
            context, InvalidOrDeadline(error), "oracle", error);
    }
    receipt.Peeled = static_cast<uint32_t>(peel.PeelOrder.size());
    receipt.Inactivated = static_cast<uint32_t>(peel.InactiveOrder.size());
    if (peel.CoreLimitExceeded)
    {
        return CompleteRun(context, EvaluationOutcome::Reject, "peel_core",
            "inactive column 1025 selected");
    }
    ResidualReceipt residual;
    if (!AnalyzeBinaryResidual(
            receipt.L, delivered.BinaryRows, peel,
            &deadline, residual, error))
    {
        return CompleteRun(
            context, InvalidOrDeadline(error), "oracle", error);
    }
    receipt.ResidualRows = residual.Rows;
    receipt.BinaryResidualRank = residual.BinaryRank;
    receipt.Deficiency = residual.Deficiency;
    if (receipt.Deficiency < kHeavyRows)
    {
        return CompleteRun(context, EvaluationOutcome::Invalid, "oracle",
            "OH0 binary deficiency is below H");
    }
    if (receipt.Deficiency > kHeavyRows)
    {
        return CompleteRun(context, EvaluationOutcome::Reject,
            "residual_rank", "binary residual deficiency exceeds H12");
    }
    return CompleteRun(context, EvaluationOutcome::Pass, "none",
        "bounded structure screen passed");
}

int EmitRunException(RunContext& context, const char* reason) noexcept
{
    try {
        return CompleteRun(
            context, EvaluationOutcome::Invalid, "oracle", reason);
    }
    catch (...)
    {
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: failed to preserve exception receipt\n");
        return kExitInvalid;
    }
}

int RunFrozenCell() noexcept
{
    try
    {
        RunContext context;
        try {
            return RunFrozenCellImpl(context);
        }
        catch (const std::bad_alloc&) {
            return EmitRunException(context, "allocation failure");
        }
        catch (const std::length_error&) {
            return EmitRunException(context, "length failure");
        }
        catch (const std::exception& exception) {
            return EmitRunException(context, exception.what());
        }
        catch (...) {
            return EmitRunException(context, "unknown exception");
        }
    }
    catch (...)
    {
        // This residual path is reachable only if the receipt context itself
        // cannot be constructed, so no in-progress state exists to preserve.
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: exception before run receipt context\n");
        return kExitInvalid;
    }
}

} // namespace

int main(int argc, char** argv)
{
#if defined(SIGPIPE)
    if (std::signal(SIGPIPE, SIG_IGN) == SIG_ERR)
    {
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: failed to ignore SIGPIPE\n");
        return kExitInvalid;
    }
#endif
    if (argc != 2)
    {
        std::fprintf(stderr, "usage: %s --selftest|--run\n", argv[0]);
        return kExitUsage;
    }
    try
    {
        const std::string command(argv[1]);
        if (command == "--selftest")
        {
            std::string error;
            if (!RunSelfTests(error))
            {
                std::fprintf(stderr,
                    "tri-bmq4-schur-v0 selftest: FAIL: %s\n",
                    error.c_str());
                return kExitInvalid;
            }
            bool emitted =
                std::puts("tri-bmq4-schur-v0 selftest: PASS") != EOF;
            if (std::fflush(stdout) != 0 || std::ferror(stdout) != 0) {
                emitted = false;
            }
            if (!emitted)
            {
                std::fprintf(stderr,
                    "tri-bmq4-schur-v0 selftest: PASS output failed\n");
                return kExitInvalid;
            }
            return kExitPass;
        }
        if (command == "--run") {
            return RunFrozenCell();
        }
        std::fprintf(stderr, "usage: %s --selftest|--run\n", argv[0]);
        return kExitUsage;
    }
    catch (const std::bad_alloc&)
    {
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: allocation failure before dispatch\n");
        return kExitInvalid;
    }
    catch (const std::length_error&)
    {
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: length failure before dispatch\n");
        return kExitInvalid;
    }
    catch (const std::exception& exception)
    {
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: exception before dispatch: %s\n",
            exception.what());
        return kExitInvalid;
    }
    catch (...)
    {
        std::fprintf(stderr,
            "tri-bmq4-schur-v0: unknown exception before dispatch\n");
        return kExitInvalid;
    }
}
