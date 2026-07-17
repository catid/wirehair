#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Codec.h"
#include "WirehairV2Peel.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

// Encoder value phase tests for the certified precode (wirehair-axd
// Phase 3, first slice):
//
// 1. Correctness on as-built systems: after ComputePrecodeValues succeeds,
//    every staircase and dense constraint row XORs to zero over the full
//    column values and every heavy row's GF(256) weighted sum is zero.
//    Seeds with a singular dense corner are skipped and counted, and the
//    success/failure verdict is cross-checked against the structure-only
//    DenseCornerInvertible on every seed.
// 2. Correctness on doctored systems (dense corner forced to identity) so
//    the dense-solve and heavy phases are exercised at every K even where
//    as-built corners are singular.
// 3. Feasibility measurement: dense-corner invertibility rate and rank
//    distribution over many seeds -- the Phase 4 seed-gating input.
// 4. Cost check against the certified model 2K + S - 1 + ceil(span/2) +
//    2*(D2 - 1) block ops for the staircase + dense known-part path.

namespace {

class MixedCoefficientPeriodScope
{
public:
    explicit MixedCoefficientPeriodScope(uint32_t period)
        : Previous(wirehair_v2::ActiveMixedCoefficientPeriod())
        , Valid(wirehair_v2::SetMixedCoefficientPeriodForTesting(period))
    {
    }

    ~MixedCoefficientPeriodScope()
    {
        (void)wirehair_v2::SetMixedCoefficientPeriodForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    uint32_t Previous;
    bool Valid;
};

class MixedGF16RowsScope
{
public:
    explicit MixedGF16RowsScope(uint32_t rows)
        : Previous(wirehair_v2::ActiveMixedGF16Rows())
        , Valid(wirehair_v2::SetMixedGF16RowsForTesting(rows))
    {
    }

    ~MixedGF16RowsScope()
    {
        (void)wirehair_v2::SetMixedGF16RowsForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    uint32_t Previous;
    bool Valid;
};

class MixedGF256RowsScope
{
public:
    explicit MixedGF256RowsScope(uint32_t rows)
        : Previous(wirehair_v2::ActiveMixedGF256Rows())
        , Valid(wirehair_v2::SetMixedGF256RowsForTesting(rows))
    {
    }

    ~MixedGF256RowsScope()
    {
        (void)wirehair_v2::SetMixedGF256RowsForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    uint32_t Previous;
    bool Valid;
};

class MixedCoefficientGeometryScope
{
public:
    explicit MixedCoefficientGeometryScope(
        wirehair_v2::MixedCoefficientGeometry geometry)
        : Previous(wirehair_v2::ActiveMixedCoefficientGeometry())
        , Valid(wirehair_v2::SetMixedCoefficientGeometryForTesting(geometry))
    {
    }

    ~MixedCoefficientGeometryScope()
    {
        (void)wirehair_v2::SetMixedCoefficientGeometryForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    wirehair_v2::MixedCoefficientGeometry Previous;
    bool Valid;
};

class MixedResidueSkewScope
{
public:
    explicit MixedResidueSkewScope(uint32_t skew)
        : Previous(wirehair_v2::ActiveMixedResidueSkew())
        , Valid(wirehair_v2::SetMixedResidueSkewForTesting(skew))
    {
    }

    ~MixedResidueSkewScope()
    {
        (void)wirehair_v2::SetMixedResidueSkewForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    uint32_t Previous;
    bool Valid;
};

class MixedResidueScheduleScope
{
public:
    explicit MixedResidueScheduleScope(
        wirehair_v2::MixedResidueSchedule schedule)
        : Previous(wirehair_v2::ActiveMixedResidueSchedule())
        , Valid(wirehair_v2::SetMixedResidueScheduleForTesting(schedule))
    {
    }

    ~MixedResidueScheduleScope()
    {
        (void)wirehair_v2::SetMixedResidueScheduleForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    wirehair_v2::MixedResidueSchedule Previous;
    bool Valid;
};

class MixedResidueHashSeedScope
{
public:
    explicit MixedResidueHashSeedScope(uint32_t seed)
        : Previous(wirehair_v2::ActiveMixedResidueHashSeed())
    {
        wirehair_v2::SetMixedResidueHashSeedForTesting(seed);
    }

    ~MixedResidueHashSeedScope()
    {
        wirehair_v2::SetMixedResidueHashSeedForTesting(Previous);
    }

private:
    uint32_t Previous;
};

class MixedIndependentExtensionResiduesScope
{
public:
    explicit MixedIndependentExtensionResiduesScope(bool enabled)
        : Valid(
            wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(
                enabled))
    {
    }

    ~MixedIndependentExtensionResiduesScope()
    {
        (void)wirehair_v2::
            SetMixedIndependentExtensionResiduesForTesting(false);
    }

    bool IsValid() const { return Valid; }

private:
    bool Valid;
};

class MixedResidueBucketModeScope
{
public:
    explicit MixedResidueBucketModeScope(
        wirehair_v2::MixedResidueBucketMode mode)
        : Previous(
            wirehair_v2::ActiveMixedResidueBucketModeForTesting())
        , Valid(wirehair_v2::SetMixedResidueBucketModeForTesting(mode))
    {
    }

    ~MixedResidueBucketModeScope()
    {
        (void)wirehair_v2::SetMixedResidueBucketModeForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    wirehair_v2::MixedResidueBucketMode Previous;
    bool Valid;
};

const char* MixedGeometryName(
    wirehair_v2::MixedCoefficientGeometry geometry)
{
    return geometry == wirehair_v2::MixedCoefficientGeometry::SharedCauchyX ?
        "shared-x" : "frozen";
}

bool ParsePositiveU32(const char* text, uint32_t& out)
{
    if (!text || !*text || *text < '0' || *text > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long value = std::strtoul(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' ||
        value == 0 || value > UINT32_MAX)
    {
        return false;
    }
    out = (uint32_t)value;
    return true;
}

const uint8_t* ColumnValue(
    const wirehair_v2::PrecodeSystem& system,
    const uint8_t* source, const uint8_t* parity,
    uint32_t block_bytes, uint32_t c)
{
    const uint32_t K = system.Params.BlockCount;
    return c < K ?
        source + (size_t)c * block_bytes :
        parity + (size_t)(c - K) * block_bytes;
}

bool ManualPacketBlock(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const uint8_t* intermediate,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint32_t L = K + P;
    const std::vector<uint32_t> row =
        wirehair_v2::GeneratePacketMatrixRow(K, P, block_id, config);
    if (!intermediate || !block_out || row.empty()) {
        return false;
    }
    std::memset(block_out, 0, block_bytes);
    for (const uint32_t column : row)
    {
        if (column >= L) {
            return false;
        }
        const uint8_t* const value =
            intermediate + (size_t)column * block_bytes;
        for (uint32_t b = 0; b < block_bytes; ++b) {
            block_out[b] ^= value[b];
        }
    }
    return true;
}

/// Verify every constraint row of the system sums to zero over the values
bool VerifyValues(
    const wirehair_v2::PrecodeSystem& system,
    const uint8_t* source, const uint8_t* parity,
    uint32_t block_bytes, const char* tag)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint32_t L = K + S + D2 + H;
    std::vector<uint8_t> acc(block_bytes);
    const wirehair_v2::MixedCoefficientRows* mixed_rows =
        system.Params.Field ==
                wirehair_v2::CompletionField::MixedGF256GF16 ?
            wirehair_v2::GetMixedCoefficientRows() : nullptr;
    if (system.Params.Field ==
            wirehair_v2::CompletionField::MixedGF256GF16 && !mixed_rows)
    {
        return false;
    }

    const auto is_zero = [&]() {
        for (uint32_t i = 0; i < block_bytes; ++i) {
            if (acc[i] != 0u) {
                return false;
            }
        }
        return true;
    };

    // Binary rows: staircase then dense
    for (uint32_t family = 0; family < 2u; ++family)
    {
        const std::vector<std::vector<uint32_t>>& rows =
            family == 0u ? system.StaircaseRows : system.DenseRowColumns;
        for (size_t r = 0; r < rows.size(); ++r)
        {
            std::memset(acc.data(), 0, block_bytes);
            for (const uint32_t col : rows[r])
            {
                gf256_add_mem(acc.data(),
                    ColumnValue(system, source, parity, block_bytes, col),
                    (int)block_bytes);
            }
            if (!is_zero())
            {
                std::fprintf(stderr,
                    "%s: K=%u bb=%u: %s row %zu does not sum to zero\n",
                    tag, K, block_bytes,
                    family == 0u ? "staircase" : "dense", r);
                return false;
            }
        }
    }

    // Completion rows: legacy uses GF(256) throughout; the mixed profile
    // embeds its leading rows in the subfield and uses GF(2^16) afterward.
    const uint32_t mixed_subfield_rows =
        wirehair_v2::ActiveMixedGF256Rows();
    for (uint32_t r = 0; r < H; ++r)
    {
        std::memset(acc.data(), 0, block_bytes);
        for (uint32_t c = 0; c < L; ++c)
        {
            uint32_t coefficient_column =
                system.Params.Field ==
                    wirehair_v2::CompletionField::MixedGF256GF16 ?
                wirehair_v2::ActiveMixedCoefficientResidue(c) : c;
            const uint8_t* v =
                ColumnValue(system, source, parity, block_bytes, c);
            if (system.Params.Field ==
                    wirehair_v2::CompletionField::MixedGF256GF16 &&
                r >= mixed_subfield_rows)
            {
                coefficient_column =
                    wirehair_v2::ActiveMixedExtensionCoefficientResidue(c);
                const uint16_t active_y = mixed_rows->Extension[
                    r - mixed_subfield_rows][coefficient_column];
                if (!wirehair_v2::GF16MulAddMem(
                        acc.data(), active_y, v, block_bytes))
                {
                    return false;
                }
            }
            else
            {
                const uint8_t y = system.Params.Field ==
                        wirehair_v2::CompletionField::MixedGF256GF16 ?
                    mixed_rows->Subfield[r][coefficient_column] :
                    wirehair_v2::HeavyCoefficient(r, coefficient_column, H);
                if (y == 1u) {
                    gf256_add_mem(acc.data(), v, (int)block_bytes);
                }
                else {
                    gf256_muladd_mem(acc.data(), y, v, (int)block_bytes);
                }
            }
        }
        if (!is_zero())
        {
            std::fprintf(stderr,
                "%s: K=%u bb=%u: heavy row %u weighted sum nonzero\n",
                tag, K, block_bytes, r);
            return false;
        }
    }
    return true;
}

void FillRandomBlocks(
    uint8_t* blocks, size_t bytes, uint64_t seed)
{
    wirehair::PCGRandom prng;
    prng.Seed(seed, 0xb10cda7au);
    for (size_t i = 0; i < bytes; ++i) {
        blocks[i] = (uint8_t)prng.Next();
    }
}

uint32_t MixedCornerRank(
    const std::vector<uint32_t>& subfield_columns,
    const std::vector<uint32_t>& extension_columns)
{
    const uint32_t subfield_rows = wirehair_v2::ActiveMixedGF256Rows();
    const uint32_t H = subfield_rows +
        wirehair_v2::ActiveMixedGF16Rows();
    if (subfield_columns.size() != H || extension_columns.size() != H) {
        return 0u;
    }
    const wirehair_v2::MixedCoefficientRows* rows =
        wirehair_v2::GetMixedCoefficientRows();
    if (!rows) return 0u;
    uint16_t matrix[
        (wirehair_v2::kMixedGF256RowsMax +
         wirehair_v2::kMixedGF16RowsMax) *
        (wirehair_v2::kMixedGF256RowsMax +
         wirehair_v2::kMixedGF16RowsMax)] = {};
    for (uint32_t r = 0; r < subfield_rows; ++r) {
        for (uint32_t j = 0; j < H; ++j) {
            matrix[(size_t)r * H + j] =
                rows->Subfield[r][subfield_columns[j]];
        }
    }
    for (uint32_t er = 0;
         er < wirehair_v2::ActiveMixedGF16Rows(); ++er) {
        const uint32_t r = subfield_rows + er;
        for (uint32_t j = 0; j < H; ++j) {
            matrix[(size_t)r * H + j] =
                rows->Extension[er][extension_columns[j]];
        }
    }

    uint32_t rank = 0u;
    for (uint32_t col = 0; col < H; ++col)
    {
        uint32_t pivot = rank;
        while (pivot < H && matrix[(size_t)pivot * H + col] == 0u) {
            ++pivot;
        }
        if (pivot >= H) continue;
        for (uint32_t k = 0; k < H; ++k) {
            std::swap(
                matrix[(size_t)rank * H + k],
                matrix[(size_t)pivot * H + k]);
        }
        const uint16_t inverse = wirehair_v2::GF16Inverse(
            matrix[(size_t)rank * H + col]);
        for (uint32_t k = 0; k < H; ++k) {
            matrix[(size_t)rank * H + k] = wirehair_v2::GF16Multiply(
                matrix[(size_t)rank * H + k], inverse);
        }
        for (uint32_t r = 0; r < H; ++r)
        {
            if (r == rank) continue;
            const uint16_t scale = matrix[(size_t)r * H + col];
            for (uint32_t k = 0; k < H; ++k) {
                matrix[(size_t)r * H + k] ^=
                    wirehair_v2::GF16Multiply(
                        scale, matrix[(size_t)rank * H + k]);
            }
        }
        ++rank;
    }
    return rank;
}

uint32_t MixedCornerRank(const std::vector<uint32_t>& columns)
{
    return MixedCornerRank(columns, columns);
}

bool TestIndependentMixedCornerAcrossKForPeriod(uint32_t period)
{
    MixedGF16RowsScope rows_scope(wirehair_v2::kMixedGF16RowsMax);
    MixedCoefficientPeriodScope period_scope(period);
    MixedCoefficientGeometryScope geometry_scope(
        wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
    MixedGF256RowsScope subfield_scope(
        wirehair_v2::kMixedGF256RowsMax);
    MixedResidueScheduleScope schedule_scope(
        wirehair_v2::MixedResidueSchedule::Hashed);
    MixedResidueHashSeedScope hash_seed_scope(68u);
    if (!rows_scope.IsValid() || !period_scope.IsValid() ||
        !geometry_scope.IsValid() || !subfield_scope.IsValid() ||
        !schedule_scope.IsValid())
    {
        return false;
    }

    const uint32_t H = wirehair_v2::ActiveMixedGF256Rows() +
        wirehair_v2::ActiveMixedGF16Rows();
    std::vector<uint32_t> subfield_columns(H);
    std::vector<uint32_t> extension_columns(H);
    for (uint32_t K = 2u; K <= 64000u; ++K)
    {
        uint32_t selected_seed = 0u;
        if (!wirehair_v2::
                SelectFullCycleMixedResidueKeyedSeedForTesting(
                    68u, K, selected_seed) ||
            wirehair_v2::ActiveMixedResidueHashSeed() != selected_seed ||
            !wirehair_v2::
                SetMixedIndependentExtensionResiduesForTesting(true))
        {
            return false;
        }
        const wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeMixedParams(K, 0u);
        const uint32_t heavy_base =
            K + params.Staircase + params.DenseRows;
        for (uint32_t j = 0u; j < H; ++j) {
            subfield_columns[j] =
                wirehair_v2::ActiveMixedCoefficientResidue(heavy_base + j);
            extension_columns[j] =
                wirehair_v2::ActiveMixedExtensionCoefficientResidue(
                    heavy_base + j);
        }
        if (MixedCornerRank(subfield_columns, extension_columns) != H)
        {
            std::fprintf(stderr,
                "mixed independent corner singular at K=%u\n", K);
            return false;
        }
    }
    (void)wirehair_v2::
        SetMixedIndependentExtensionResiduesForTesting(false);
    std::printf(
        "mixed independent P%u H16 keyed corners K=2..64000: PASS\n",
        period);
    return true;
}

bool TestIndependentMixedCornerAcrossK()
{
    return TestIndependentMixedCornerAcrossKForPeriod(31u) &&
        TestIndependentMixedCornerAcrossKForPeriod(32u);
}

bool TestMixedCornerRankForGeometry(
    wirehair_v2::MixedCoefficientGeometry geometry,
    uint32_t extension_rows,
    const uint32_t* periods,
    size_t period_count)
{
    MixedGF16RowsScope rows_scope(extension_rows);
    MixedCoefficientGeometryScope geometry_scope(geometry);
    if (!rows_scope.IsValid() || !geometry_scope.IsValid()) {
        return false;
    }
    const uint32_t H = wirehair_v2::kMixedGF256Rows +
        wirehair_v2::ActiveMixedGF16Rows();
    std::vector<uint32_t> columns(H);
    wirehair::PCGRandom prng;
    prng.Seed(
        UINT64_C(0x16c0a4e7) ^ (uint32_t)geometry,
        UINT64_C(0x244));
    for (size_t period_index = 0;
         period_index < period_count; ++period_index)
    {
        const uint32_t period = periods[period_index];
        for (uint32_t start = 0; start < period; ++start)
        {
            for (uint32_t j = 0; j < H; ++j) {
                columns[j] = (start + j) % period;
            }
            if (MixedCornerRank(columns) != H) {
                std::fprintf(stderr,
                    "mixed corner: period %u consecutive start %u is "
                    "singular geometry=%s\n",
                    period, start, MixedGeometryName(geometry));
                return false;
            }
        }

        std::vector<uint32_t> deck(period);
        for (uint32_t trial = 0; trial < 10000u; ++trial)
        {
            for (uint32_t i = 0; i < period; ++i) deck[i] = i;
            for (uint32_t i = 0; i < H; ++i)
            {
                const uint32_t pick = i + prng.Next() % (period - i);
                std::swap(deck[i], deck[pick]);
                columns[i] = deck[i];
            }
            if (MixedCornerRank(columns) != H) {
                std::fprintf(stderr,
                    "mixed corner: period %u nonconsecutive trial %u is "
                    "singular geometry=%s\n",
                    period, trial, MixedGeometryName(geometry));
                return false;
            }
        }
    }
    std::printf(
        "mixed %ux%u corner rank geometry=%s (%zu periods, all starts + "
        "10000 samples each): PASS\n",
        H, H, MixedGeometryName(geometry), period_count);
    return true;
}

bool TestMixedCornerRank()
{
    // Frozen smaller periods are sweep points, not certified coefficient
    // sets.  Period 64 was rejected after producing a sampled singular
    // arbitrary corner.  The shared-X construction should be Cauchy/MDS at
    // every valid period because it retains distinct X and Y coordinates.
    const uint32_t frozen_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u
    };
    const uint32_t shared_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 16u, 12u
    };
    const uint32_t shared_h13_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 16u, 13u
    };
    const uint32_t shared_h14_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 16u, 14u
    };
    if (!TestMixedCornerRankForGeometry(
            wirehair_v2::MixedCoefficientGeometry::FrozenPowerX,
            wirehair_v2::kMixedGF16Rows,
            frozen_periods,
            sizeof(frozen_periods) / sizeof(frozen_periods[0])) ||
        !TestMixedCornerRankForGeometry(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16Rows,
            shared_periods,
            sizeof(shared_periods) / sizeof(shared_periods[0])) ||
        !TestMixedCornerRankForGeometry(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16Rows + 1u,
            shared_h13_periods,
            sizeof(shared_h13_periods) / sizeof(shared_h13_periods[0])) ||
        !TestMixedCornerRankForGeometry(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            shared_h14_periods,
            sizeof(shared_h14_periods) / sizeof(shared_h14_periods[0])))
    {
        return false;
    }

    const uint32_t original_period =
        wirehair_v2::ActiveMixedCoefficientPeriod();
    const wirehair_v2::MixedCoefficientGeometry original_geometry =
        wirehair_v2::ActiveMixedCoefficientGeometry();
    const uint32_t original_rows = wirehair_v2::ActiveMixedGF16Rows();
    const uint32_t original_skew =
        wirehair_v2::ActiveMixedResidueSkew();
    const wirehair_v2::MixedResidueSchedule original_schedule =
        wirehair_v2::ActiveMixedResidueSchedule();
    const uint32_t original_hash_seed =
        wirehair_v2::ActiveMixedResidueHashSeed();
    if (wirehair_v2::SetMixedCoefficientPeriodForTesting(
            wirehair_v2::kMixedGF256Rows +
                wirehair_v2::kMixedGF16Rows - 1u) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(
            wirehair_v2::kMixedCoefficientPeriod + 1u) ||
        wirehair_v2::SetMixedCoefficientGeometryForTesting(
            static_cast<wirehair_v2::MixedCoefficientGeometry>(2u)) ||
        wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16Rows - 1u) ||
        wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16RowsMax + 1u) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16RowsMax) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(
            wirehair_v2::kMixedGF256Rows +
                wirehair_v2::kMixedGF16Rows) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(original_rows) ||
        wirehair_v2::ActiveMixedCoefficientPeriod() != original_period ||
        wirehair_v2::ActiveMixedCoefficientGeometry() != original_geometry ||
        wirehair_v2::ActiveMixedGF16Rows() != original_rows)
    {
        std::fprintf(stderr,
            "mixed corner: invalid experiment override was accepted\n");
        return false;
    }
    if (wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256Rows - 1u) ||
        wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256RowsMax + 1u) ||
        wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256RowsMax) ||
        !wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX) ||
        wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256RowsMax) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16RowsMax) ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(32u) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256RowsMax) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(28u) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(30u) ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(31u) ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(32u) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(34u) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(29u) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(33u) ||
        wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::FrozenPowerX) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256Rows) ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(original_period) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(original_rows) ||
        !wirehair_v2::SetMixedCoefficientGeometryForTesting(
            original_geometry) ||
        wirehair_v2::ActiveMixedGF256Rows() !=
            wirehair_v2::kMixedGF256Rows)
    {
        std::fprintf(stderr,
            "mixed corner: GF256 row override contract failed\n");
        return false;
    }

    bool skew_corners_ok =
        wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16RowsMax) &&
        wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
    static const struct {
        uint32_t Period;
        uint32_t Skew;
    } kSkewCases[] = {
        { 29u, 1u },
        { 29u, 14u },
        { 29u, 15u },
        { 32u, 1u },
        { 32u, 18u }
    };
    for (const auto& skew_case : kSkewCases)
    {
        skew_corners_ok = skew_corners_ok &&
            wirehair_v2::SetMixedCoefficientPeriodForTesting(
                skew_case.Period) &&
            !wirehair_v2::SetMixedResidueSkewForTesting(
                skew_case.Period - 14u + 1u) &&
            wirehair_v2::SetMixedResidueSkewForTesting(skew_case.Skew) &&
            wirehair_v2::ActiveMixedResidueSkew() == skew_case.Skew;
        for (uint32_t start = 0u;
             start < skew_case.Period * skew_case.Period && skew_corners_ok;
             ++start)
        {
            bool seen[wirehair_v2::kMixedCoefficientPeriod] = {};
            for (uint32_t j = 0u; j < 14u; ++j)
            {
                const uint32_t residue =
                    wirehair_v2::ActiveMixedCoefficientResidue(start + j);
                if (residue >= skew_case.Period || seen[residue]) {
                    skew_corners_ok = false;
                    break;
                }
                seen[residue] = true;
            }
        }
    }
    static const wirehair_v2::MixedResidueSchedule kAperiodicSchedules[] = {
        wirehair_v2::MixedResidueSchedule::Ramp,
        wirehair_v2::MixedResidueSchedule::Hashed
    };
    for (const auto schedule : kAperiodicSchedules)
    {
        skew_corners_ok = skew_corners_ok &&
            wirehair_v2::SetMixedCoefficientPeriodForTesting(28u) &&
            wirehair_v2::SetMixedResidueSkewForTesting(0u) &&
            wirehair_v2::SetMixedResidueScheduleForTesting(schedule) &&
            wirehair_v2::ActiveMixedResidueSchedule() == schedule;
        for (uint32_t start = 0u;
             start < 127u * 28u * 28u && skew_corners_ok;
             ++start)
        {
            bool seen[wirehair_v2::kMixedCoefficientPeriod] = {};
            for (uint32_t j = 0u; j < 14u; ++j)
            {
                const uint32_t residue =
                    wirehair_v2::ActiveMixedCoefficientResidue(start + j);
                if (residue >= 28u || seen[residue]) {
                    skew_corners_ok = false;
                    break;
                }
                seen[residue] = true;
            }
        }
    }
    for (uint32_t period = 28u; period <= 32u && skew_corners_ok; ++period)
    {
        wirehair_v2::SetMixedResidueHashSeedForTesting(0u);
        skew_corners_ok =
            wirehair_v2::SetMixedCoefficientPeriodForTesting(period) &&
            wirehair_v2::SetMixedResidueSkewForTesting(0u) &&
            wirehair_v2::SetMixedResidueScheduleForTesting(
                wirehair_v2::MixedResidueSchedule::Hashed);
        uint32_t previous_shift =
            wirehair_v2::ActiveMixedResidueBlockShift(0u);
        for (uint32_t block = 0u; block < 127u && skew_corners_ok; ++block)
        {
            const uint32_t next_shift =
                wirehair_v2::ActiveMixedResidueBlockShift(block + 1u);
            const uint32_t step =
                (next_shift + period - previous_shift) % period;
            if (step < 1u || step > period - 14u) {
                skew_corners_ok = false;
            }
            previous_shift = next_shift;
        }
        for (uint32_t cycle = 1u; cycle < period && skew_corners_ok; ++cycle)
        {
            if (wirehair_v2::ActiveMixedResidueBlockShift(127u * cycle) ==
                0u)
            {
                skew_corners_ok = false;
            }
        }
        if (wirehair_v2::ActiveMixedResidueBlockShift(127u * period) != 0u) {
            skew_corners_ok = false;
        }
    }
    uint32_t seed_zero_shifts[128] = {};
    skew_corners_ok = skew_corners_ok &&
        wirehair_v2::SetMixedCoefficientPeriodForTesting(28u) &&
        wirehair_v2::SetMixedResidueScheduleForTesting(
            wirehair_v2::MixedResidueSchedule::Hashed);
    wirehair_v2::SetMixedResidueHashSeedForTesting(0u);
    for (uint32_t block = 0u; block <= 127u; ++block) {
        seed_zero_shifts[block] =
            wirehair_v2::ActiveMixedResidueBlockShift(block);
    }
    wirehair_v2::SetMixedResidueHashSeedForTesting(7u);
    bool seeded_sequence_differs = false;
    uint32_t previous_shift =
        wirehair_v2::ActiveMixedResidueBlockShift(0u);
    for (uint32_t block = 0u; block < 127u && skew_corners_ok; ++block)
    {
        const uint32_t next_shift =
            wirehair_v2::ActiveMixedResidueBlockShift(block + 1u);
        const uint32_t step = (next_shift + 28u - previous_shift) % 28u;
        if (step < 1u || step > 14u) {
            skew_corners_ok = false;
        }
        seeded_sequence_differs = seeded_sequence_differs ||
            next_shift != seed_zero_shifts[block + 1u];
        previous_shift = next_shift;
    }
    skew_corners_ok = skew_corners_ok && seeded_sequence_differs;

    bool independent_schedule_ok = skew_corners_ok &&
        wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(true) &&
        wirehair_v2::ActiveMixedIndependentExtensionResidues() &&
        wirehair_v2::ActiveMixedExtensionResidueBlockShift(0u) == 0u;
    bool independent_sequence_differs = false;
    previous_shift =
        wirehair_v2::ActiveMixedExtensionResidueBlockShift(0u);
    for (uint32_t block = 0u;
         block < 127u && independent_schedule_ok; ++block)
    {
        const uint32_t next_shift =
            wirehair_v2::ActiveMixedExtensionResidueBlockShift(block + 1u);
        const uint32_t step = (next_shift + 28u - previous_shift) % 28u;
        if (step < 1u || step > 14u) {
            independent_schedule_ok = false;
        }
        independent_sequence_differs = independent_sequence_differs ||
            next_shift !=
                wirehair_v2::ActiveMixedResidueBlockShift(block + 1u);
        previous_shift = next_shift;
    }
    for (uint32_t cycle = 1u;
         cycle < 28u && independent_schedule_ok; ++cycle)
    {
        if (wirehair_v2::ActiveMixedExtensionResidueBlockShift(
                127u * cycle) == 0u)
        {
            independent_schedule_ok = false;
        }
    }
    independent_schedule_ok = independent_schedule_ok &&
        independent_sequence_differs &&
        wirehair_v2::ActiveMixedExtensionResidueBlockShift(127u * 28u) ==
            0u &&
        wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(false) &&
        !wirehair_v2::ActiveMixedIndependentExtensionResidues();
    for (uint32_t block = 0u;
         block < 128u && independent_schedule_ok; ++block)
    {
        independent_schedule_ok =
            wirehair_v2::ActiveMixedExtensionResidueBlockShift(block) ==
            wirehair_v2::ActiveMixedResidueBlockShift(block);
    }
    independent_schedule_ok = independent_schedule_ok &&
        wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(true);
    wirehair_v2::SetMixedResidueHashSeedForTesting(8u);
    independent_schedule_ok = independent_schedule_ok &&
        !wirehair_v2::ActiveMixedIndependentExtensionResidues() &&
        wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(true) &&
        wirehair_v2::SetMixedResidueScheduleForTesting(
            wirehair_v2::MixedResidueSchedule::Constant) &&
        !wirehair_v2::ActiveMixedIndependentExtensionResidues() &&
        !wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(true);
    skew_corners_ok = skew_corners_ok && independent_schedule_ok;

    wirehair_v2::SetMixedResidueHashSeedForTesting(original_hash_seed);
    const bool skew_restored =
        wirehair_v2::SetMixedCoefficientPeriodForTesting(original_period) &&
        wirehair_v2::SetMixedCoefficientGeometryForTesting(original_geometry) &&
        wirehair_v2::SetMixedGF16RowsForTesting(original_rows) &&
        wirehair_v2::SetMixedResidueSkewForTesting(original_skew) &&
        wirehair_v2::SetMixedResidueScheduleForTesting(original_schedule);
    if (!skew_corners_ok || !skew_restored)
    {
        std::fprintf(stderr,
            "mixed corner: balanced residue skew invariant failed\n");
        return false;
    }
    return true;
}

/// GF(2) rank of the D2 x D2 dense corner (structure only)
uint32_t DenseCornerRank(const wirehair_v2::PrecodeSystem& system)
{
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t dense_base =
        system.Params.BlockCount + system.Params.Staircase;
    std::vector<uint64_t> masks(D2, 0);
    for (uint32_t r = 0; r < D2; ++r)
    {
        for (const uint32_t col : system.DenseRowColumns[r])
        {
            if (col >= dense_base) {
                masks[r] ^= UINT64_C(1) << (col - dense_base);
            }
        }
    }
    uint32_t rank = 0;
    for (uint32_t col = 0; col < D2; ++col)
    {
        const uint64_t bit = UINT64_C(1) << col;
        uint32_t pivot = rank;
        while (pivot < D2 && 0u == (masks[pivot] & bit)) {
            ++pivot;
        }
        if (pivot >= D2) {
            continue;
        }
        std::swap(masks[rank], masks[pivot]);
        for (uint32_t r = 0; r < D2; ++r) {
            if (r != rank && 0u != (masks[r] & bit)) {
                masks[r] ^= masks[rank];
            }
        }
        ++rank;
    }
    return rank;
}

bool TestCorrectnessAsBuilt()
{
    const uint32_t Ks[] = {2u, 16u, 1000u, 3200u};
    const uint32_t bbs[] = {1u, 13u, 1280u};
    const uint32_t seeds = 30u;
    bool ok = true;

    for (const uint32_t K : Ks)
    {
        for (const uint32_t bb : bbs)
        {
            uint32_t succeeded = 0, singular = 0;
            for (uint32_t seed = 0; seed < seeds; ++seed)
            {
                wirehair_v2::PrecodeSystem system;
                if (!BuildPrecodeSystem(
                        wirehair_v2::MakeCertifiedParams(K, seed), system))
                {
                    std::fprintf(stderr, "K=%u: build failed\n", K);
                    return false;
                }
                const uint32_t parity_count = system.Params.Staircase +
                    system.Params.DenseRows + system.Params.HeavyRows;
                std::vector<uint8_t> source((size_t)K * bb);
                std::vector<uint8_t> parity((size_t)parity_count * bb);
                FillRandomBlocks(source.data(), source.size(),
                    seed * 1000003ull + K * 17ull + bb);

                const bool feasible =
                    wirehair_v2::DenseCornerInvertible(system);
                const bool computed = wirehair_v2::ComputePrecodeValues(
                    system, source.data(), bb, parity.data());
                if (computed != feasible)
                {
                    std::fprintf(stderr,
                        "K=%u bb=%u seed=%u: ComputePrecodeValues=%d but "
                        "DenseCornerInvertible=%d\n",
                        K, bb, seed, (int)computed, (int)feasible);
                    return false;
                }
                if (!computed)
                {
                    ++singular;
                    continue;
                }
                ++succeeded;
                if (!VerifyValues(system, source.data(), parity.data(),
                        bb, "as-built"))
                {
                    ok = false;
                }
            }
            std::printf(
                "as-built K=%5u bb=%4u: %2u/%u seeds feasible "
                "(%u singular dense corners)\n",
                K, bb, succeeded, seeds, singular);
        }
    }
    return ok;
}

bool TestCorrectnessDoctored()
{
    const uint32_t Ks[] = {16u, 1000u, 3200u};
    const uint32_t bbs[] = {1u, 13u, 1280u};

    for (const uint32_t K : Ks)
    {
        for (const uint32_t bb : bbs)
        {
            for (uint32_t seed = 100; seed < 103u; ++seed)
            {
                wirehair_v2::PrecodeParams params =
                    wirehair_v2::MakeCertifiedParams(K, seed);
                params.DenseIdentityCorner = true;
                wirehair_v2::PrecodeSystem system;
                if (!BuildPrecodeSystem(params, system))
                {
                    std::fprintf(stderr, "K=%u: build failed\n", K);
                    return false;
                }
                if (!wirehair_v2::DenseCornerInvertible(system))
                {
                    std::fprintf(stderr,
                        "K=%u: doctored corner reported singular\n", K);
                    return false;
                }

                const uint32_t parity_count = system.Params.Staircase +
                    system.Params.DenseRows + system.Params.HeavyRows;
                std::vector<uint8_t> source((size_t)K * bb);
                std::vector<uint8_t> parity((size_t)parity_count * bb);
                FillRandomBlocks(source.data(), source.size(),
                    seed * 999983ull + K * 29ull + bb);

                if (!wirehair_v2::ComputePrecodeValues(
                        system, source.data(), bb, parity.data()))
                {
                    std::fprintf(stderr,
                        "K=%u bb=%u seed=%u: doctored encode failed\n",
                        K, bb, seed);
                    return false;
                }
                if (!VerifyValues(system, source.data(), parity.data(),
                        bb, "doctored"))
                {
                    return false;
                }
            }
        }
    }
    std::printf("doctored-corner correctness: PASS\n");
    return true;
}

bool TestTwoAnchorPhasedEncoder()
{
    // The public phased encoder is not the intended benchmark path for this
    // experiment, but it must still evaluate the deliberately dense row 6->7
    // difference correctly whenever the full-span dense corner is feasible.
    const uint32_t K = 2u;
    const uint32_t block_bytes = 13u;
    for (uint32_t seed = 0u; seed < 10000u; ++seed)
    {
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(K, seed);
        params.DenseTwoAnchor = true;
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(params, system)) {
            return false;
        }
        if (!wirehair_v2::DenseCornerInvertible(system)) {
            continue;
        }

        const uint32_t parity_count = system.Params.Staircase +
            system.Params.DenseRows + system.Params.HeavyRows;
        std::vector<uint8_t> source((size_t)K * block_bytes);
        std::vector<uint8_t> parity((size_t)parity_count * block_bytes);
        FillRandomBlocks(
            source.data(), source.size(), UINT64_C(0x7a12e0c0) + seed);
        if (!wirehair_v2::ComputePrecodeValues(
                system, source.data(), block_bytes, parity.data()) ||
            !VerifyValues(
                system, source.data(), parity.data(), block_bytes,
                "two-anchor"))
        {
            std::fprintf(stderr,
                "two-anchor phased encode failed at seed=%u\n", seed);
            return false;
        }
        std::printf(
            "two-anchor phased encoder seed=%u: PASS\n", seed);
        return true;
    }
    std::fprintf(stderr,
        "two-anchor phased encoder found no feasible tiny corner\n");
    return false;
}

bool TestFeasibility(uint32_t trials)
{
    const uint32_t Ks[] = {1000u, 3200u, 32000u};

    std::printf(
        "\ndense-corner feasibility (%u seeds per K; Phase 4 gating "
        "input):\n", trials);
    for (const uint32_t K : Ks)
    {
        uint32_t invertible = 0;
        uint32_t max_rank = 0;
        uint64_t rank_sum = 0;
        for (uint32_t seed = 0; seed < trials; ++seed)
        {
            wirehair_v2::PrecodeSystem system;
            if (!BuildPrecodeSystem(
                    wirehair_v2::MakeCertifiedParams(K, seed), system))
            {
                std::fprintf(stderr, "K=%u: build failed\n", K);
                return false;
            }
            const uint32_t rank = DenseCornerRank(system);
            rank_sum += rank;
            max_rank = std::max(max_rank, rank);
            const bool inv = rank == system.Params.DenseRows;
            if (inv != wirehair_v2::DenseCornerInvertible(system))
            {
                std::fprintf(stderr,
                    "K=%u seed=%u: DenseCornerInvertible disagrees with "
                    "rank\n", K, seed);
                return false;
            }
            if (inv) {
                ++invertible;
            }
        }
        std::printf(
            "  K=%5u: invertible %u/%u (%.4f%%), corner rank mean %.2f "
            "max %u of %u\n",
            K, invertible, trials, 100.0 * invertible / trials,
            (double)rank_sum / trials, max_rank, 12u);
    }
    return true;
}

bool TestCostModel()
{
    const uint32_t Ks[] = {16u, 1000u, 3200u, 10000u};
    bool ok = true;

    std::printf("\ncost check vs model min(N1,S)*K + S - 1 + ceil(span/2) + "
        "2*(D2-1):\n");
    for (const uint32_t K : Ks)
    {
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(
                wirehair_v2::MakeCertifiedParams(K, 0x5eedu), system))
        {
            std::fprintf(stderr, "K=%u: build failed\n", K);
            return false;
        }
        const uint32_t S = system.Params.Staircase;
        const uint32_t D2 = system.Params.DenseRows;
        const uint32_t span = K + S + D2;
        const uint32_t bb = 64u;

        const uint32_t parity_count =
            S + D2 + system.Params.HeavyRows;
        std::vector<uint8_t> source((size_t)K * bb);
        std::vector<uint8_t> parity((size_t)parity_count * bb);
        FillRandomBlocks(source.data(), source.size(), K);

        // Stats for the staircase + dense known-part path are recorded even
        // when the dense corner is singular (the common case at these K)
        wirehair_v2::PrecodeEncodeStats stats;
        (void)wirehair_v2::ComputePrecodeValues(
            system, source.data(), bb, parity.data(), &stats);

        const uint32_t hits =
            system.Params.SourceHits < S ? system.Params.SourceHits : S;
        const uint64_t staircase_model = (uint64_t)hits * K + S - 1u;
        const uint64_t dense_model =
            ((span + 1u) >> 1) + 2ull * (D2 - 1u);
        const uint64_t measured =
            stats.StaircaseBlockOps + stats.DenseKnownBlockOps;
        const uint64_t model = staircase_model + dense_model;

        std::printf(
            "  K=%5u: staircase %llu (model %llu), dense known %llu "
            "(model %llu), total %llu vs %llu (%+lld)\n",
            K,
            (unsigned long long)stats.StaircaseBlockOps,
            (unsigned long long)staircase_model,
            (unsigned long long)stats.DenseKnownBlockOps,
            (unsigned long long)dense_model,
            (unsigned long long)measured,
            (unsigned long long)model,
            (long long)(measured - model));

        if (stats.StaircaseBlockOps != staircase_model)
        {
            std::fprintf(stderr,
                "K=%u: staircase ops off model\n", K);
            ok = false;
        }
        // The dense known part undershoots the model by exactly the number
        // of dense-column touches (row-0 dense columns + dense-hitting
        // flips), each of which costs a mask toggle instead of a block op
        if (stats.DenseKnownBlockOps > dense_model ||
            dense_model - stats.DenseKnownBlockOps > D2 + 2u * (D2 - 1u))
        {
            std::fprintf(stderr,
                "K=%u: dense known-part ops outside model window\n", K);
            ok = false;
        }
    }
    return ok;
}

bool TestHeavyResidueDispatch()
{
    const uint32_t K = 500u;
    const uint32_t bb = 37u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x57ea4));
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) {
        std::fprintf(stderr, "heavy residue dispatch: build failed\n");
        return false;
    }
    const uint32_t parity_count = params.Staircase + params.DenseRows +
        params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * bb);
    std::vector<uint8_t> full((size_t)parity_count * bb);
    std::vector<uint8_t> streamed(full.size());
    FillRandomBlocks(source.data(), source.size(), UINT64_C(0x57ea45eed));

    wirehair_v2::PrecodeEncodeStats full_stats;
    wirehair_v2::PrecodeEncodeStats streamed_stats;
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(UINT64_MAX);
    const bool full_ok = wirehair_v2::ComputePrecodeValues(
        system, source.data(), bb, full.data(), &full_stats);
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(0u);
    const bool streamed_ok = wirehair_v2::ComputePrecodeValues(
        system, source.data(), bb, streamed.data(), &streamed_stats);
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(UINT64_C(64) << 20);

    const bool stats_equal =
        full_stats.StaircaseBlockOps == streamed_stats.StaircaseBlockOps &&
        full_stats.DenseKnownBlockOps == streamed_stats.DenseKnownBlockOps &&
        full_stats.DenseSolveBlockOps == streamed_stats.DenseSolveBlockOps &&
        full_stats.HeavyBucketXors == streamed_stats.HeavyBucketXors &&
        full_stats.HeavyMulAdds == streamed_stats.HeavyMulAdds &&
        full_stats.HeavySolveBlockOps == streamed_stats.HeavySolveBlockOps;
    const uint32_t heavy_base = K + params.Staircase + params.DenseRows;
    const uint32_t window = 256u - params.HeavyRows;
    if (!full_ok || !streamed_ok || full != streamed || !stats_equal ||
        streamed_stats.HeavyBucketXors != heavy_base ||
        streamed_stats.HeavyMulAdds != (uint64_t)params.HeavyRows * window ||
        !VerifyValues(
            system, source.data(), streamed.data(), bb,
            "heavy residue streaming"))
    {
        std::fprintf(stderr,
            "heavy residue dispatch: full/streaming mismatch\n");
        return false;
    }
    std::printf("heavy residue full/streaming differential: PASS\n");
    return true;
}

bool StatsAreZero(const wirehair_v2::PrecodeEncodeStats& stats);

bool TestMixedCompletionForPeriod(
    uint32_t period,
    wirehair_v2::MixedCoefficientGeometry geometry,
    uint32_t extension_rows,
    uint32_t residue_skew = 0u,
    wirehair_v2::MixedResidueSchedule residue_schedule =
        wirehair_v2::MixedResidueSchedule::Constant,
    bool independent_extension_residues = false,
    uint32_t subfield_rows = wirehair_v2::kMixedGF256Rows,
    wirehair_v2::MixedResidueBucketMode bucket_mode =
        wirehair_v2::MixedResidueBucketMode::Automatic)
{
    MixedCoefficientGeometryScope geometry_scope(geometry);
    MixedGF16RowsScope rows_scope(extension_rows);
    MixedCoefficientPeriodScope period_scope(period);
    MixedGF256RowsScope subfield_scope(subfield_rows);
    MixedResidueSkewScope skew_scope(residue_skew);
    MixedResidueScheduleScope schedule_scope(residue_schedule);
    MixedResidueHashSeedScope hash_seed_scope(
        independent_extension_residues ? 68u :
            wirehair_v2::ActiveMixedResidueHashSeed());
    MixedIndependentExtensionResiduesScope independent_scope(
        independent_extension_residues);
    MixedResidueBucketModeScope bucket_mode_scope(bucket_mode);
    if (!geometry_scope.IsValid() || !subfield_scope.IsValid() ||
        !rows_scope.IsValid() || !period_scope.IsValid() ||
        !skew_scope.IsValid() ||
        !schedule_scope.IsValid() || !independent_scope.IsValid() ||
        !bucket_mode_scope.IsValid())
    {
        std::fprintf(stderr,
            "mixed completion: invalid test scope g=%u s=%u r=%u p=%u "
            "k=%u q=%u i=%u\n",
            geometry_scope.IsValid() ? 1u : 0u,
            subfield_scope.IsValid() ? 1u : 0u,
            rows_scope.IsValid() ? 1u : 0u,
            period_scope.IsValid() ? 1u : 0u,
            skew_scope.IsValid() ? 1u : 0u,
            schedule_scope.IsValid() ? 1u : 0u,
            independent_scope.IsValid() ? 1u : 0u);
        return false;
    }
    // Exercise the exact production geometry (default full-span dense
    // corner).  Most such direct-source corners are singular by design, so
    // select a deterministically feasible tiny seed; the public packet solve
    // covers larger systems independently.
    bool production_geometry_checked = false;
    for (uint64_t seed = 0; seed < 1000u && !production_geometry_checked;
         ++seed)
    {
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeMixedParams(2u, seed);
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(params, system) ||
            !wirehair_v2::DenseCornerInvertible(system))
        {
            continue;
        }
        const uint32_t bb = 16u;
        const uint32_t parity_count = params.Staircase +
            params.DenseRows + params.HeavyRows;
        std::vector<uint8_t> source((size_t)params.BlockCount * bb);
        std::vector<uint8_t> parity((size_t)parity_count * bb);
        FillRandomBlocks(source.data(), source.size(), seed ^ UINT64_C(0xcafe));
        if (!wirehair_v2::ComputePrecodeValues(
                system, source.data(), bb, parity.data()) ||
            !VerifyValues(
                system, source.data(), parity.data(), bb,
                "mixed production geometry"))
        {
            std::fprintf(stderr,
                "mixed completion: production geometry oracle failed\n");
            return false;
        }
        production_geometry_checked = true;
    }
    if (!production_geometry_checked) {
        std::fprintf(stderr,
            "mixed completion: no feasible production geometry seed\n");
        return false;
    }

    const uint32_t Ks[] = {16u, 500u, 1000u};
    const uint32_t bbs[] = {2u, 16u, 1280u};
    for (uint32_t K : Ks)
    {
        for (uint32_t bb : bbs)
        {
            wirehair_v2::PrecodeParams params =
                wirehair_v2::MakeMixedParams(K, UINT64_C(0x16c0de) + K + bb);
            params.DenseIdentityCorner = true;
            wirehair_v2::PrecodeSystem system;
            if (!BuildPrecodeSystem(params, system)) {
                std::fprintf(stderr, "mixed completion: build failed\n");
                return false;
            }
            const uint32_t parity_count = params.Staircase +
                params.DenseRows + params.HeavyRows;
            std::vector<uint8_t> source((size_t)K * bb);
            std::vector<uint8_t> full((size_t)parity_count * bb, 0xa5u);
            std::vector<uint8_t> streamed(full.size(), 0x5au);
            FillRandomBlocks(
                source.data(), source.size(), UINT64_C(0x16f1e1d) ^ K ^ bb);

            wirehair_v2::PrecodeEncodeStats full_stats;
            wirehair_v2::PrecodeEncodeStats streamed_stats;
            wirehair_v2::SetHeavyBucketStorageLimitForTesting(UINT64_MAX);
            const bool full_ok = wirehair_v2::ComputePrecodeValues(
                system, source.data(), bb, full.data(), &full_stats);
            wirehair_v2::SetHeavyBucketStorageLimitForTesting(0u);
            const bool streamed_ok = wirehair_v2::ComputePrecodeValues(
                system, source.data(), bb, streamed.data(), &streamed_stats);
            wirehair_v2::SetHeavyBucketStorageLimitForTesting(
                UINT64_C(64) << 20);

            const uint32_t heavy_base =
                K + params.Staircase + params.DenseRows;
            const uint32_t coefficient_period =
                wirehair_v2::ActiveMixedCoefficientPeriod();
            const bool bucketed =
                heavy_base >= 2u * coefficient_period;
            const bool joint_bucketed =
                bucketed && independent_extension_residues &&
                bucket_mode ==
                    wirehair_v2::MixedResidueBucketMode::JointDelta;
            const bool dual_bucketed =
                bucketed && independent_extension_residues &&
                bucket_mode ==
                    wirehair_v2::MixedResidueBucketMode::Dual;
            const uint64_t residues = bucketed ?
                coefficient_period : heavy_base;
            const uint64_t separate_bucket_xors = bucketed ?
                (uint64_t)heavy_base *
                    (independent_extension_residues ? 2u : 1u) :
                0u;
            const bool joint_stats_ok = !joint_bucketed ||
                (full_stats.MixedJointSourceXors == heavy_base &&
                 full_stats.MixedJointActiveDeltas > 0u &&
                 full_stats.MixedJointMarginalXors ==
                    (uint64_t)2u * coefficient_period *
                        (full_stats.MixedJointActiveDeltas - 1u) &&
                 full_stats.MixedJointMarginalCopies ==
                    (uint64_t)2u * coefficient_period &&
                 full_stats.MixedJointScratchBytes ==
                    (uint64_t)3u * coefficient_period * bb);
            const bool dual_stats_ok = !dual_bucketed ||
                (full_stats.MixedDualSourceColumns == heavy_base &&
                 full_stats.MixedJointSourceXors == 0u &&
                 full_stats.MixedJointMarginalXors == 0u &&
                 full_stats.MixedJointMarginalCopies == 0u &&
                 full_stats.MixedJointScratchBytes == 0u &&
                 full_stats.MixedJointActiveDeltas == 0u);
            const uint64_t expected_full_bucket_xors = joint_bucketed ?
                full_stats.MixedJointSourceXors +
                    full_stats.MixedJointMarginalXors :
                separate_bucket_xors;
            if (!full_ok || !streamed_ok || full != streamed ||
                !VerifyValues(
                    system, source.data(), streamed.data(), bb, "mixed") ||
                !joint_stats_ok ||
                !dual_stats_ok ||
                full_stats.HeavyBucketXors != expected_full_bucket_xors ||
                streamed_stats.HeavyBucketXors != separate_bucket_xors ||
                full_stats.HeavyMulAdds != params.HeavyRows * residues ||
                full_stats.MixedGF16MulAdds !=
                    wirehair_v2::ActiveMixedGF16Rows() * residues ||
                full_stats.MixedPlaneConversions !=
                    residues + wirehair_v2::ActiveMixedGF256Rows() +
                        params.HeavyRows ||
                full_stats.HeavyMulAdds != streamed_stats.HeavyMulAdds ||
                full_stats.HeavySolveBlockOps !=
                    streamed_stats.HeavySolveBlockOps ||
                full_stats.MixedGF16MulAdds !=
                    streamed_stats.MixedGF16MulAdds ||
                full_stats.MixedGF16SolveBlockOps !=
                    streamed_stats.MixedGF16SolveBlockOps ||
                full_stats.MixedPlaneConversions !=
                    streamed_stats.MixedPlaneConversions)
            {
                std::fprintf(stderr,
                    "mixed completion: K=%u bb=%u differential/oracle failed\n",
                    K, bb);
                return false;
            }

            if (K == 16u && bb == 16u)
            {
                const wirehair_v2::PeelingCodec codec =
                    wirehair_v2::MakePeelingCodec(
                        wirehair_v2::PeelStructure::LtM1C32,
                        wirehair_v2::PeelSolver::KsBmaxTop16);
                wirehair_v2::PrecodeEncoder encoder;
                if (encoder.InitializeResult(
                        system, codec, UINT64_C(0x51eed), 5u,
                        source.data(), bb) != Wirehair_Success ||
                    encoder.System().Params.Field !=
                        wirehair_v2::CompletionField::MixedGF256GF16 ||
                    std::memcmp(
                        encoder.ParityBlocks(), full.data(), full.size()) != 0 ||
                    encoder.EncodeStats().MixedGF16MulAdds !=
                        full_stats.MixedGF16MulAdds ||
                    encoder.EncodeStats().MixedGF16SolveBlockOps !=
                        full_stats.MixedGF16SolveBlockOps)
                {
                    std::fprintf(stderr,
                        "mixed completion: stateful encoder mismatch\n");
                    return false;
                }
            }
        }
    }

    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeMixedParams(16u, UINT64_C(0x0dd));
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) return false;
    const uint32_t parity_count =
        params.Staircase + params.DenseRows + params.HeavyRows;
    std::vector<uint8_t> source((size_t)params.BlockCount * 3u, 0x11u);
    std::vector<uint8_t> parity((size_t)parity_count * 3u, 0x5au);
    const std::vector<uint8_t> before = parity;
    wirehair_v2::PrecodeEncodeStats stats;
    if (wirehair_v2::ComputePrecodeValues(
            system, source.data(), 3u, parity.data(), &stats) ||
        parity != before || !StatsAreZero(stats))
    {
        std::fprintf(stderr,
            "mixed completion: odd block rejection was not transactional\n");
        return false;
    }
    const wirehair_v2::PeelingCodec codec = wirehair_v2::MakePeelingCodec(
        wirehair_v2::PeelStructure::LtM1C32,
        wirehair_v2::PeelSolver::KsBmaxTop16);
    wirehair_v2::PrecodeEncoder encoder;
    if (encoder.InitializeResult(
            system, codec, UINT64_C(0x51eed), 5u,
            source.data(), 3u) != Wirehair_InvalidInput ||
        encoder.IsInitialized())
    {
        std::fprintf(stderr,
            "mixed completion: odd stateful encoder classification failed\n");
        return false;
    }

    std::printf(
        "mixed completion scalar/full/streamed oracle period=%u geometry=%s "
        "gf256_rows=%u gf16_rows=%u skew=%u schedule=%u "
        "independent_extension=%u bucket_mode=%u: "
        "PASS\n",
        period, MixedGeometryName(geometry), subfield_rows, extension_rows,
        residue_skew,
        (uint32_t)residue_schedule,
        independent_extension_residues ? 1u : 0u,
        (uint32_t)bucket_mode);
    return true;
}

bool TestMixedCompletion()
{
    const uint32_t frozen_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u
    };
    const uint32_t shared_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 12u
    };
    for (const uint32_t period : frozen_periods) {
        if (!TestMixedCompletionForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::FrozenPowerX,
                wirehair_v2::kMixedGF16Rows))
        {
            return false;
        }
    }
    for (const uint32_t period : shared_periods) {
        if (!TestMixedCompletionForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16Rows))
        {
            return false;
        }
    }
    const uint32_t shared_h13_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 13u
    };
    const uint32_t shared_h14_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 14u
    };
    for (const uint32_t period : shared_h13_periods) {
        if (!TestMixedCompletionForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16Rows + 1u))
        {
            return false;
        }
    }
    for (const uint32_t period : shared_h14_periods) {
        if (!TestMixedCompletionForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16RowsMax))
        {
            return false;
        }
    }
    if (!TestMixedCompletionForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            18u) ||
        !TestMixedCompletionForPeriod(
            28u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Ramp) ||
        !TestMixedCompletionForPeriod(
            28u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed) ||
        !TestMixedCompletionForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true) ||
        !TestMixedCompletionForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true,
            wirehair_v2::kMixedGF256RowsMax) ||
        !TestMixedCompletionForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true,
            wirehair_v2::kMixedGF256Rows + 1u,
            wirehair_v2::MixedResidueBucketMode::Dual) ||
        !TestMixedCompletionForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true,
            wirehair_v2::kMixedGF256Rows + 1u,
            wirehair_v2::MixedResidueBucketMode::JointDelta))
    {
        return false;
    }
    return true;
}

bool TestAutomaticMixedJointBucketPolicy()
{
    if (wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            3199u, 4096u, 32u) ||
        !wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            3200u, 4096u, 32u) ||
        wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            9999u, 1280u, 32u) ||
        !wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            10000u, 1280u, 32u) ||
        wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            64000u, 1279u, 32u) ||
        wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            64000u, 4096u, 31u) ||
        wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            64000u, 4096u, 33u) ||
        wirehair_v2::UseAutomaticMixedJointResidueBuckets(
            64000u, 4096u, 244u) ||
        !wirehair_v2::UseAutomaticMixedJointResidueBucketsForTesting(
            3200u, 4096u, 32u) ||
        !wirehair_v2::MixedJointResidueBucketStorageFits(
            32u, 699050u,
            wirehair_v2::kMixedJointResidueBucketDataByteCap) ||
        wirehair_v2::MixedJointResidueBucketStorageFits(
            32u, 699051u,
            wirehair_v2::kMixedJointResidueBucketDataByteCap) ||
        wirehair_v2::MixedJointResidueBucketStorageFits(
            0u, 4096u,
            wirehair_v2::kMixedJointResidueBucketDataByteCap) ||
        wirehair_v2::MixedJointResidueBucketStorageFits(
            wirehair_v2::kMixedCoefficientPeriod + 1u, 4096u,
            wirehair_v2::kMixedJointResidueBucketDataByteCap) ||
        wirehair_v2::MixedJointResidueBucketStorageFits(
            32u, 0u,
            wirehair_v2::kMixedJointResidueBucketDataByteCap))
    {
        std::fprintf(stderr, "mixed automatic bucket policy boundary failed\n");
        return false;
    }

    MixedCoefficientGeometryScope geometry_scope(
        wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
    MixedGF16RowsScope extension_scope(wirehair_v2::kMixedGF16RowsMax);
    MixedCoefficientPeriodScope period_scope(32u);
    MixedGF256RowsScope subfield_scope(wirehair_v2::kMixedGF256Rows + 1u);
    MixedResidueSkewScope skew_scope(0u);
    MixedResidueScheduleScope schedule_scope(
        wirehair_v2::MixedResidueSchedule::Hashed);
    MixedResidueHashSeedScope hash_scope(68u);
    MixedIndependentExtensionResiduesScope independent_scope(true);
    MixedResidueBucketModeScope bucket_scope(
        wirehair_v2::MixedResidueBucketMode::Automatic);
    if (!geometry_scope.IsValid() || !extension_scope.IsValid() ||
        !period_scope.IsValid() || !subfield_scope.IsValid() ||
        !skew_scope.IsValid() || !schedule_scope.IsValid() ||
        !independent_scope.IsValid() || !bucket_scope.IsValid())
    {
        return false;
    }

    static const uint32_t K = 3200u;
    static const uint32_t block_bytes = 4096u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeMixedParams(K, UINT64_C(0xa070b0c0d0e0f001));
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) return false;
    const uint32_t parity_count =
        params.Staircase + params.DenseRows + params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * block_bytes);
    std::vector<uint8_t> parity((size_t)parity_count * block_bytes);
    std::vector<uint8_t> capped(parity.size());
    FillRandomBlocks(source.data(), source.size(), UINT64_C(0xa070f111));
    wirehair_v2::PrecodeEncodeStats stats;
    const bool joint_ok = wirehair_v2::ComputePrecodeValues(
        system, source.data(), block_bytes, parity.data(), &stats);
    const uint64_t joint_plane_bytes = UINT64_C(32) * block_bytes;
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(
        3u * joint_plane_bytes - 1u);
    wirehair_v2::PrecodeEncodeStats capped_stats;
    const bool capped_ok = wirehair_v2::ComputePrecodeValues(
        system, source.data(), block_bytes, capped.data(), &capped_stats);
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(
        wirehair_v2::kMixedJointResidueBucketDataByteCap);
    if (!joint_ok || !capped_ok || parity != capped ||
        !VerifyValues(
            system, source.data(), parity.data(), block_bytes,
            "mixed automatic joint buckets") ||
        stats.MixedJointSourceXors !=
            (uint64_t)K + params.Staircase + params.DenseRows ||
        stats.MixedJointMarginalCopies != 64u ||
        stats.MixedJointActiveDeltas == 0u ||
        stats.MixedDualSourceColumns != 0u ||
        capped_stats.MixedJointSourceXors != 0u ||
        capped_stats.MixedJointMarginalXors != 0u ||
        capped_stats.MixedJointMarginalCopies != 0u ||
        capped_stats.MixedJointScratchBytes != 0u ||
        capped_stats.MixedJointActiveDeltas != 0u ||
        capped_stats.MixedDualSourceColumns != 0u)
    {
        std::fprintf(stderr,
            "mixed automatic bucket dispatch/accounting failed\n");
        return false;
    }
    return true;
}

bool TestMalformedDenseCorner()
{
    wirehair_v2::PrecodeSystem system;
    system.Params.BlockCount = 2u;
    system.Params.Staircase = 1u;
    system.Params.DenseRows = 1u;
    system.Params.HeavyRows = 0u;
    system.Params.SourceHits = 1u;
    system.Params.DenseIdentityCorner = false;
    system.Params.Seed = 0u;
    system.StaircaseRows.push_back(std::vector<uint32_t>{0u, 2u});

    // dense_base is K + S = 3.  Column 67 is dense_base + 64, which used
    // to shift by 64 in DenseCornerInvertible before row validation.  Keep
    // a valid dense column at the tail so the test also catches validators
    // that only scan the sorted dense tail and miss malformed unsorted rows.
    system.DenseRowColumns.push_back(std::vector<uint32_t>{67u, 0u, 3u});
    if (wirehair_v2::DenseCornerInvertible(system)) {
        std::fprintf(stderr,
            "malformed dense corner should not be reported invertible\n");
        return false;
    }

    system.DenseRowColumns[0] = std::vector<uint32_t>{3u, 3u};
    if (wirehair_v2::DenseCornerInvertible(system)) {
        std::fprintf(stderr,
            "duplicate dense columns should cancel in GF(2)\n");
        return false;
    }
    return true;
}

bool StatsAreZero(const wirehair_v2::PrecodeEncodeStats& stats)
{
    return stats.StaircaseBlockOps == 0u &&
        stats.DenseKnownBlockOps == 0u &&
        stats.DenseSolveBlockOps == 0u &&
        stats.HeavyBucketXors == 0u &&
        stats.HeavyMulAdds == 0u &&
        stats.HeavySolveBlockOps == 0u &&
        stats.MixedGF16MulAdds == 0u &&
        stats.MixedGF16SolveBlockOps == 0u &&
        stats.MixedPlaneConversions == 0u &&
        stats.MixedJointSourceXors == 0u &&
        stats.MixedJointMarginalXors == 0u &&
        stats.MixedJointMarginalCopies == 0u &&
        stats.MixedJointScratchBytes == 0u &&
        stats.MixedJointActiveDeltas == 0u &&
        stats.MixedDualSourceColumns == 0u;
}

bool TestStrictSystemValidation()
{
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(16u, UINT64_C(0x51a1d));
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem valid;
    if (!BuildPrecodeSystem(params, valid) ||
        !wirehair_v2::ValidatePrecodeSystem(valid))
    {
        std::fprintf(stderr, "strict validation: valid builder system failed\n");
        return false;
    }

    std::vector<wirehair_v2::PrecodeSystem> invalid;
    wirehair_v2::PrecodeSystem bad;

    bad = valid;
    bad.Params.BlockCount = 0u;
    invalid.push_back(bad);
    bad = valid;
    bad.Params.Staircase = UINT32_MAX;
    invalid.push_back(bad);
    bad = valid;
    bad.Params.DenseRows = UINT32_MAX;
    invalid.push_back(bad);
    bad = valid;
    bad.Params.HeavyRows = 129u;
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows.pop_back();
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns.pop_back();
    invalid.push_back(bad);

    const uint32_t K = valid.Params.BlockCount;
    const uint32_t S = valid.Params.Staircase;
    const uint32_t dense_base = K + S;
    const uint32_t binary_span = dense_base + valid.Params.DenseRows;

    bad = valid;
    bad.StaircaseRows[0].erase(
        std::find(bad.StaircaseRows[0].begin(),
            bad.StaircaseRows[0].end(), K));
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows[1].erase(
        std::find(bad.StaircaseRows[1].begin(),
            bad.StaircaseRows[1].end(), K));
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows[1].push_back(K + 2u);
    std::sort(bad.StaircaseRows[1].begin(), bad.StaircaseRows[1].end());
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows[0].insert(
        bad.StaircaseRows[0].begin(), bad.StaircaseRows[0][0]);
    invalid.push_back(bad);
    bad = valid;
    std::swap(bad.StaircaseRows[0][0], bad.StaircaseRows[0][1]);
    invalid.push_back(bad);

    bad = valid;
    bad.DenseRowColumns[0].insert(
        bad.DenseRowColumns[0].begin(), bad.DenseRowColumns[0][0]);
    invalid.push_back(bad);
    bad = valid;
    std::swap(bad.DenseRowColumns[0][0], bad.DenseRowColumns[0][1]);
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns[0].push_back(binary_span);
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns[0].back() = dense_base + 1u;
    std::sort(
        bad.DenseRowColumns[0].begin(), bad.DenseRowColumns[0].end());
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns[1].erase(bad.DenseRowColumns[1].begin());
    invalid.push_back(bad);

    wirehair::PCGRandom mutation_rng;
    mutation_rng.Seed(UINT64_C(0xf0225eed), UINT64_C(0x51a1d));
    for (uint32_t trial = 0; trial < 256u; ++trial)
    {
        bad = valid;
        const uint32_t row = mutation_rng.Next() % S;
        switch (mutation_rng.Next() % 8u)
        {
        case 0:
            bad.Params.BlockCount = mutation_rng.Next() & 1u;
            break;
        case 1:
            bad.StaircaseRows[row].erase(std::find(
                bad.StaircaseRows[row].begin(),
                bad.StaircaseRows[row].end(), K + row));
            break;
        case 2:
            bad.StaircaseRows[row].insert(
                bad.StaircaseRows[row].begin(),
                bad.StaircaseRows[row][0]);
            break;
        case 3:
            bad.DenseRowColumns[mutation_rng.Next() % valid.Params.DenseRows]
                .push_back(binary_span);
            break;
        case 4:
        {
            std::vector<uint32_t>& dense = bad.DenseRowColumns[
                mutation_rng.Next() % valid.Params.DenseRows];
            std::swap(dense[0], dense[1]);
            break;
        }
        case 5:
        {
            const uint32_t dense_row =
                mutation_rng.Next() % valid.Params.DenseRows;
            std::vector<uint32_t>& dense = bad.DenseRowColumns[dense_row];
            dense.erase(std::find(
                dense.begin(), dense.end(), dense_base + dense_row));
            break;
        }
        case 6:
            bad.Params.HeavyRows = 129u + mutation_rng.Next() % 100u;
            break;
        default:
            bad.DenseRowColumns.pop_back();
            break;
        }
        invalid.push_back(bad);
    }

    const uint32_t bb = 7u;
    const uint32_t parity_count = S + valid.Params.DenseRows +
        valid.Params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * bb, 0x3cu);
    for (size_t i = 0; i < invalid.size(); ++i)
    {
        std::vector<uint8_t> parity((size_t)parity_count * bb, 0xa5u);
        const std::vector<uint8_t> before = parity;
        wirehair_v2::PrecodeEncodeStats stats;
        stats.StaircaseBlockOps = UINT64_MAX;
        if (wirehair_v2::ValidatePrecodeSystem(invalid[i]) ||
            wirehair_v2::DenseCornerInvertible(invalid[i]) ||
            wirehair_v2::ComputePrecodeValues(
                invalid[i], source.data(), bb, parity.data(), &stats) ||
            parity != before || !StatsAreZero(stats))
        {
            std::fprintf(stderr,
                "strict validation: malformed case %zu was accepted/written\n",
                i);
            return false;
        }
    }

    wirehair_v2::PrecodeSystem singular;
    if (!BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(16u, 0u), singular) ||
        !wirehair_v2::ValidatePrecodeSystem(singular) ||
        wirehair_v2::DenseCornerInvertible(singular))
    {
        std::fprintf(stderr,
            "strict validation: expected structurally valid singular system\n");
        return false;
    }
    std::vector<uint8_t> singular_parity((size_t)parity_count * bb, 0xa5u);
    wirehair_v2::PrecodeEncodeStats singular_stats;
    if (wirehair_v2::ComputePrecodeValues(
            singular, source.data(), bb, singular_parity.data(),
            &singular_stats) ||
        singular_stats.StaircaseBlockOps == 0u ||
        singular_stats.DenseKnownBlockOps == 0u)
    {
        std::fprintf(stderr,
            "strict validation: singular system lost partial cost behavior\n");
        return false;
    }

    wirehair::PCGRandom property_rng;
    property_rng.Seed(UINT64_C(0xc0de51a1d), UINT64_C(0x5eed));
    for (uint32_t trial = 0; trial < 32u; ++trial)
    {
        const uint32_t random_k = 16u + property_rng.Next() % 1009u;
        wirehair_v2::PrecodeParams random_params =
            wirehair_v2::MakeCertifiedParams(random_k, property_rng.Next());
        random_params.DenseIdentityCorner = true;
        wirehair_v2::PrecodeSystem random_system;
        if (!BuildPrecodeSystem(random_params, random_system) ||
            !wirehair_v2::ValidatePrecodeSystem(random_system))
        {
            std::fprintf(stderr,
                "strict validation: random feasible build failed\n");
            return false;
        }
        const uint32_t random_parity_count = random_params.Staircase +
            random_params.DenseRows + random_params.HeavyRows;
        std::vector<uint8_t> random_source((size_t)random_k * 3u);
        std::vector<uint8_t> random_parity(
            (size_t)random_parity_count * 3u);
        FillRandomBlocks(
            random_source.data(), random_source.size(), property_rng.Next());
        if (!wirehair_v2::ComputePrecodeValues(
                random_system, random_source.data(), 3u,
                random_parity.data()) ||
            !VerifyValues(
                random_system, random_source.data(), random_parity.data(),
                3u, "random feasible property"))
        {
            std::fprintf(stderr,
                "strict validation: random values violated constraints\n");
            return false;
        }
    }

    std::printf("strict precode validation: PASS\n");
    return true;
}

void ManualRecoveryBlock(
    const wirehair_v2::PrecodeSystem& system,
    const uint8_t* source, const uint8_t* parity,
    uint32_t block_bytes,
    const std::vector<uint32_t>& row_columns,
    uint8_t* block_out)
{
    std::memset(block_out, 0, block_bytes);
    for (const uint32_t col : row_columns)
    {
        gf256_add_mem(block_out,
            ColumnValue(system, source, parity, block_bytes, col),
            (int)block_bytes);
    }
}

bool EqualBlock(
    const uint8_t* a, const uint8_t* b, uint32_t block_bytes)
{
    return std::memcmp(a, b, block_bytes) == 0;
}

bool IsZeroBlock(const uint8_t* a, uint32_t block_bytes)
{
    for (uint32_t i = 0; i < block_bytes; ++i) {
        if (a[i] != 0u) {
            return false;
        }
    }
    return true;
}

bool TestRecoveryBlockEncoding()
{
    const uint32_t K = 1000u;
    const uint32_t bb = 37u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x5eed1234));
    params.DenseIdentityCorner = true;

    wirehair_v2::PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) {
        std::fprintf(stderr, "recovery encode: build failed\n");
        return false;
    }

    const uint32_t parity_count =
        system.Params.Staircase + system.Params.DenseRows +
        system.Params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * bb);
    std::vector<uint8_t> parity((size_t)parity_count * bb);
    FillRandomBlocks(source.data(), source.size(), UINT64_C(0xbeadfeed));
    wirehair_v2::PrecodeEncodeStats precode_stats;
    if (!wirehair_v2::ComputePrecodeValues(
            system, source.data(), bb, parity.data(), &precode_stats))
    {
        std::fprintf(stderr, "recovery encode: precode values failed\n");
        return false;
    }
    if (!VerifyValues(
            system, source.data(), parity.data(), bb, "recovery encode"))
    {
        return false;
    }

    const wirehair_v2::PeelingCodec codec =
        wirehair_v2::MakePeelingCodec(
            wirehair_v2::PeelStructure::LtM1C32,
            wirehair_v2::PeelSolver::KsBmaxTop16);
    const uint64_t row_seed = UINT64_C(0xdec0de);
    const uint32_t recovery_mix = 5u;
    const std::vector<std::vector<uint32_t> > rows =
        wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, parity_count, 8u, recovery_mix, row_seed);
    if (rows.size() != 8u) {
        std::fprintf(stderr, "recovery encode: row generation failed\n");
        return false;
    }

    std::vector<uint8_t> got(bb), want(bb);
    for (size_t r = 0; r < rows.size(); ++r)
    {
        uint64_t ops = UINT64_MAX;
        if (!wirehair_v2::ComputeRecoveryBlock(
                system, source.data(), parity.data(), bb, rows[r],
                got.data(), &ops))
        {
            std::fprintf(stderr,
                "recovery encode: ComputeRecoveryBlock failed for row %zu\n",
                r);
            return false;
        }
        ManualRecoveryBlock(
            system, source.data(), parity.data(), bb, rows[r], want.data());
        if (!EqualBlock(got.data(), want.data(), bb)) {
            std::fprintf(stderr,
                "recovery encode: row %zu did not match manual XOR\n", r);
            return false;
        }
        if (ops != rows[r].size()) {
            std::fprintf(stderr,
                "recovery encode: row %zu ops %llu != columns %zu\n",
                r, (unsigned long long)ops, rows[r].size());
            return false;
        }
    }

    wirehair_v2::PrecodeEncoder uninitialized_encoder;
    uint64_t ops = UINT64_MAX;
    if (uninitialized_encoder.HasCompleteSystem() ||
        uninitialized_encoder.Encode(0u, got.data(), &ops) ||
        ops != UINT64_MAX)
    {
        std::fprintf(stderr,
            "recovery encode: uninitialized encoder should fail\n");
        return false;
    }

    wirehair_v2::PrecodeEncoder encoder_state;
    if (!encoder_state.Initialize(
            system, codec, row_seed, recovery_mix, source.data(), bb))
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder initialize failed\n");
        return false;
    }
    if (!encoder_state.IsInitialized() ||
        !encoder_state.HasCompleteSystem() ||
        !wirehair_v2::ValidatePrecodeSystem(encoder_state.System()) ||
        encoder_state.System().StaircaseRows != system.StaircaseRows ||
        encoder_state.System().DenseRowColumns != system.DenseRowColumns ||
        encoder_state.SourceBlockCount() != K ||
        encoder_state.ParityBlockCount() != parity_count ||
        encoder_state.BlockBytes() != bb ||
        !encoder_state.ParityBlocks())
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder accessors failed\n");
        return false;
    }
    if (encoder_state.EncodeStats().StaircaseBlockOps !=
            precode_stats.StaircaseBlockOps ||
        encoder_state.EncodeStats().DenseKnownBlockOps !=
            precode_stats.DenseKnownBlockOps ||
        encoder_state.EncodeStats().DenseSolveBlockOps !=
            precode_stats.DenseSolveBlockOps ||
        encoder_state.EncodeStats().HeavyBucketXors !=
            precode_stats.HeavyBucketXors ||
        encoder_state.EncodeStats().HeavyMulAdds !=
            precode_stats.HeavyMulAdds ||
        encoder_state.EncodeStats().HeavySolveBlockOps !=
            precode_stats.HeavySolveBlockOps)
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder stats mismatch\n");
        return false;
    }
    if (std::memcmp(
            encoder_state.ParityBlocks(), parity.data(),
            (size_t)parity_count * bb) != 0)
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder parity mismatch\n");
        return false;
    }
    if (encoder_state.Initialize(
            system, codec, row_seed, recovery_mix, nullptr, bb) ||
        !encoder_state.IsInitialized() ||
        !encoder_state.HasCompleteSystem() ||
        !wirehair_v2::ValidatePrecodeSystem(encoder_state.System()) ||
        !encoder_state.ParityBlocks() ||
        encoder_state.BlockBytes() != bb)
    {
        std::fprintf(stderr,
            "recovery encode: failed initialize should preserve state\n");
        return false;
    }
    if (!encoder_state.Initialize(
            system, codec, row_seed, recovery_mix, source.data(), bb))
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder reinitialize failed\n");
        return false;
    }

    const uint32_t source_id = 17u;
    if (!wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), nullptr, bb, source_id, got.data(), &ops) ||
        !EqualBlock(
            got.data(), source.data() + (size_t)source_id * bb, bb) ||
        ops != 1u)
    {
        std::fprintf(stderr, "recovery encode: encoded source block failed\n");
        return false;
    }
    if (!encoder_state.Encode(source_id, got.data(), &ops) ||
        !EqualBlock(
            got.data(), source.data() + (size_t)source_id * bb, bb) ||
        ops != 1u)
    {
        std::fprintf(stderr,
            "recovery encode: state encoded source block failed\n");
        return false;
    }

    const uint32_t recovery_index = 4u;
    if (!wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), parity.data(), bb, K + recovery_index,
            got.data(), &ops))
    {
        std::fprintf(stderr,
            "recovery encode: encoded recovery block failed\n");
        return false;
    }
    ManualRecoveryBlock(
        system, source.data(), parity.data(), bb,
        rows[recovery_index], want.data());
    if (!EqualBlock(got.data(), want.data(), bb) ||
        ops != rows[recovery_index].size())
    {
        std::fprintf(stderr,
            "recovery encode: encoded recovery block mismatch\n");
        return false;
    }
    if (!encoder_state.Encode(K + recovery_index, got.data(), &ops) ||
        !EqualBlock(got.data(), want.data(), bb) ||
        ops != rows[recovery_index].size())
    {
        std::fprintf(stderr,
            "recovery encode: state encoded recovery block mismatch\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), nullptr, bb, K + recovery_index,
            got.data(), &ops) ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }) ||
        ops != 0u)
    {
        std::fprintf(stderr,
            "recovery encode: recovery block should require parity values\n");
        return false;
    }

    const uint32_t high_block_id = UINT32_MAX;
    const uint32_t high_row_id = high_block_id - K;
    const std::vector<uint32_t> high_row =
        wirehair_v2::GenerateRecoveryMatrixRow(
            codec, K, parity_count, high_row_id, recovery_mix, row_seed);
    ManualRecoveryBlock(
        system, source.data(), parity.data(), bb, high_row, want.data());
    if (!wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), parity.data(), bb,
            high_block_id, got.data(), &ops) ||
        !EqualBlock(got.data(), want.data(), bb) ||
        ops != high_row.size())
    {
        std::fprintf(stderr,
            "recovery encode: full uint32 block id failed\n");
        return false;
    }
    if (!encoder_state.Encode(high_block_id, got.data(), &ops) ||
        !EqualBlock(got.data(), want.data(), bb) ||
        ops != high_row.size())
    {
        std::fprintf(stderr,
            "recovery encode: state full uint32 block id failed\n");
        return false;
    }

    const std::vector<uint32_t> source_only{0u};
    if (!wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, source_only,
            got.data(), &ops) ||
        !EqualBlock(got.data(), source.data(), bb) || ops != 1u)
    {
        std::fprintf(stderr, "recovery encode: source-only row failed\n");
        return false;
    }

    const std::vector<uint32_t> parity_only{K + system.Params.Staircase};
    if (!wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, parity_only,
            got.data(), &ops) ||
        !EqualBlock(
            got.data(), parity.data() + (size_t)system.Params.Staircase * bb,
            bb) ||
        ops != 1u)
    {
        std::fprintf(stderr, "recovery encode: parity-only row failed\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    const std::vector<uint32_t> empty_row;
    if (!wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, empty_row,
            got.data(), &ops) ||
        !IsZeroBlock(got.data(), bb) || ops != 0u)
    {
        std::fprintf(stderr, "recovery encode: empty row failed\n");
        return false;
    }

    const std::vector<uint32_t> bad_row{
        0u,
        K + parity_count
    };
    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, bad_row,
            got.data(), &ops) ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }) ||
        ops != 0u)
    {
        std::fprintf(stderr,
            "recovery encode: invalid row should fail without writes\n");
        return false;
    }

    if (wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), 0u, source_only,
            got.data()))
    {
        std::fprintf(stderr,
            "recovery encode: zero block_bytes should fail\n");
        return false;
    }

    std::printf("recovery block encoding: PASS\n");
    return true;
}

bool TestMessagePrecodeEncoder()
{
    const uint32_t K = 1000u;
    const uint32_t bb = 37u;
    const uint32_t tail = 13u;
    const uint64_t message_bytes = (uint64_t)(K - 1u) * bb + tail;
    std::vector<uint8_t> message((size_t)message_bytes);
    FillRandomBlocks(message.data(), message.size(), UINT64_C(0x1234beef));

    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.DenseIdentityCorner = true;
    options.RecoveryMixCount = wirehair_v2::kCertifiedPacketMixCount;

    wirehair_v2::MessagePrecodeEncoder encoder;
    if (!encoder.Initialize(
            message.data(), message_bytes, bb, nullptr, &options))
    {
        std::fprintf(stderr,
            "message encoder: identity-corner initialize failed\n");
        return false;
    }
    if (!encoder.IsInitialized() ||
        encoder.MessageBytes() != message_bytes ||
        encoder.SourceBlockCount() != K ||
        encoder.BlockBytes() != bb ||
        encoder.Profile().BlockCount != K ||
        encoder.Profile().BlockBytes != bb ||
        encoder.Options().RecoveryMixCount != options.RecoveryMixCount ||
        !encoder.Options().DenseIdentityCorner ||
        !encoder.IntermediateBlocks() ||
        !encoder.BlockEncoder().IsInitialized())
    {
        std::fprintf(stderr, "message encoder: accessors failed\n");
        return false;
    }

    std::vector<uint8_t> got(bb, 0xac), want(bb, 0xbd);
    uint32_t data_bytes = UINT32_MAX;
    uint64_t ops = UINT64_MAX;
    if (!encoder.Encode(0u, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != bb ||
        ops == 0u ||
        std::memcmp(got.data(), message.data(), bb) != 0)
    {
        std::fprintf(stderr, "message encoder: source block failed\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    const uint8_t* tail_src =
        message.data() + (size_t)(K - 1u) * bb;
    if (!encoder.Encode(K - 1u, got.data(), tail, &data_bytes, &ops) ||
        data_bytes != tail ||
        ops == 0u ||
        std::memcmp(got.data(), tail_src, tail) != 0 ||
        !std::all_of(got.begin() + tail, got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr, "message encoder: final partial source failed\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (!encoder.Encode(K - 1u, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != tail ||
        std::memcmp(got.data(), tail_src, tail) != 0 ||
        !std::all_of(got.begin() + tail, got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr,
            "message encoder: full-capacity partial source over-wrote tail\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    data_bytes = UINT32_MAX;
    ops = UINT64_MAX;
    if (encoder.Encode(K - 1u, got.data(), tail - 1u, &data_bytes, &ops) ||
        data_bytes != UINT32_MAX ||
        ops != UINT64_MAX ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr,
            "message encoder: short partial-source buffer should fail\n");
        return false;
    }

    const uint32_t recovery_id = K + 3u;
    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (!encoder.Encode(recovery_id, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != bb ||
        ops == 0u)
    {
        std::fprintf(stderr, "message encoder: recovery block failed\n");
        return false;
    }
    uint64_t want_ops = UINT64_MAX;
    if (!encoder.BlockEncoder().Encode(
            recovery_id, want.data(), &want_ops) ||
        want_ops != ops ||
        !EqualBlock(got.data(), want.data(), bb))
    {
        std::fprintf(stderr,
            "message encoder: recovery block mismatch\n");
        return false;
    }

    const wirehair_v2::PrecodeEncoder& blocks = encoder.BlockEncoder();
    const wirehair_v2::PrecodeSystem& encoder_system = blocks.System();
    wirehair_v2::PacketRowConfig packet_config;
    packet_config.PeelSeed = (uint32_t)blocks.RecoveryRowSeed();
    packet_config.MixCount = blocks.RecoveryMixCount();
    if (!blocks.IntermediateBlocks() || blocks.HasCompleteSystem() ||
        !encoder_system.StaircaseRows.empty() ||
        !encoder_system.DenseRowColumns.empty() ||
        encoder_system.Params.BlockCount != K ||
        encoder_system.Params.Staircase != encoder.Profile().V2StaircaseCount ||
        encoder_system.Params.DenseRows != encoder.Profile().V2DenseRowCount ||
        encoder_system.Params.HeavyRows != encoder.Profile().V2HeavyRowCount)
    {
        std::fprintf(stderr,
            "message encoder: solved state retained or lost system data\n");
        return false;
    }
    const uint32_t reordered_ids[] = {
        K + 9u, 2u, K + 1u, K - 1u, K + 15u, 0u, K + 4u
    };
    for (uint32_t packet_id : reordered_ids)
    {
        std::fill(got.begin(), got.end(), uint8_t{0xac});
        std::fill(want.begin(), want.end(), uint8_t{0xbd});
        data_bytes = UINT32_MAX;
        if (!encoder.Encode(
                packet_id, got.data(), bb, &data_bytes, &ops))
        {
            std::fprintf(stderr,
                "message encoder: reordered packet %u failed\n", packet_id);
            return false;
        }
        if (packet_id < K)
        {
            const uint32_t expected_bytes =
                packet_id + 1u < K ? bb : tail;
            if (data_bytes != expected_bytes ||
                std::memcmp(
                    got.data(),
                    message.data() + (size_t)packet_id * bb,
                    expected_bytes) != 0)
            {
                std::fprintf(stderr,
                    "message encoder: reordered source oracle mismatch\n");
                return false;
            }
        }
        else
        {
            if (!ManualPacketBlock(
                encoder_system,
                packet_config,
                blocks.IntermediateBlocks(),
                bb,
                packet_id,
                want.data()))
            {
                std::fprintf(stderr,
                    "message encoder: packet oracle evaluation failed\n");
                return false;
            }
            if (data_bytes != bb || !EqualBlock(got.data(), want.data(), bb))
            {
                std::fprintf(stderr,
                    "message encoder: reordered repair oracle mismatch\n");
                return false;
            }
        }
    }

    const uint32_t max_packet_id = UINT32_MAX;
    if (!ManualPacketBlock(
            encoder_system, packet_config, blocks.IntermediateBlocks(), bb,
            max_packet_id, want.data()))
    {
        std::fprintf(stderr,
            "message encoder: maximum packet-id oracle failed\n");
        return false;
    }
    if (!encoder.Encode(
            max_packet_id, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != bb || !EqualBlock(got.data(), want.data(), bb))
    {
        std::fprintf(stderr,
            "message encoder: maximum packet-id oracle mismatch\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    data_bytes = UINT32_MAX;
    ops = UINT64_MAX;
    if (encoder.Encode(recovery_id, got.data(), bb - 1u, &data_bytes, &ops) ||
        data_bytes != UINT32_MAX ||
        ops != UINT64_MAX ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr,
            "message encoder: short recovery buffer should fail\n");
        return false;
    }

    wirehair_v2::SeedProfile mismatch =
        wirehair_v2::SelectSeedProfile(K, bb);
    ++mismatch.BlockCount;
    if (encoder.Initialize(
            message.data(), message_bytes, bb, &mismatch, &options) ||
        !encoder.IsInitialized() ||
        !encoder.IntermediateBlocks() ||
        encoder.BlockEncoder().HasCompleteSystem() ||
        !encoder.BlockEncoder().System().StaircaseRows.empty() ||
        !encoder.BlockEncoder().System().DenseRowColumns.empty() ||
        encoder.SourceBlockCount() != K ||
        !encoder.BlockEncoder().IsInitialized())
    {
        std::fprintf(stderr,
            "message encoder: failed reinitialize should preserve state\n");
        return false;
    }

    if (encoder.Initialize(nullptr, message_bytes, bb, nullptr, &options) ||
        !encoder.IsInitialized())
    {
        std::fprintf(stderr,
            "message encoder: null message should fail\n");
        return false;
    }
    const uint8_t dummy = 0u;
    if (encoder.Initialize(
            &dummy,
            UINT64_C(0x100000000),
            UINT32_C(0x80000000),
            nullptr,
            &options) ||
        !encoder.IsInitialized())
    {
        std::fprintf(stderr,
            "message encoder: oversized block should fail before allocation\n");
        return false;
    }

    std::printf("message precode encoder: PASS\n");
    return true;
}

bool TestBorrowedMessageLifetime()
{
    const uint32_t K = 64u;
    const uint32_t bb = 67u;
    const uint32_t final_sizes[] = { bb, 13u };
    const uint32_t packet_ids[] = { 0u, K - 1u, K + 7u };

    for (uint32_t final_bytes : final_sizes)
    {
        const uint64_t message_bytes =
            (uint64_t)(K - 1u) * bb + final_bytes;
        std::vector<uint8_t> message((size_t)message_bytes);
        FillRandomBlocks(
            message.data(), message.size(),
            UINT64_C(0xb0770ed) ^ final_bytes);

        wirehair_v2::MessagePrecodeEncoder encoder;
        if (encoder.InitializeResult(
                message.data(), message_bytes, bb) != Wirehair_Success)
        {
            std::fprintf(stderr,
                "borrowed lifetime: initialize failed for final=%u\n",
                final_bytes);
            return false;
        }

        std::vector<std::vector<uint8_t> > expected;
        std::vector<uint32_t> expected_bytes;
        for (uint32_t packet_id : packet_ids)
        {
            expected.push_back(std::vector<uint8_t>(bb, 0xa5u));
            uint32_t data_bytes = UINT32_MAX;
            if (encoder.EncodeResult(
                    packet_id, expected.back().data(), bb, &data_bytes) !=
                    Wirehair_Success)
            {
                std::fprintf(stderr,
                    "borrowed lifetime: initial encode failed id=%u\n",
                    packet_id);
                return false;
            }
            expected_bytes.push_back(data_bytes);
        }

        if (final_bytes != bb)
        {
            std::vector<uint8_t> padded(bb, 0xa5u);
            uint64_t ops = UINT64_MAX;
            if (!encoder.BlockEncoder().Encode(
                    K - 1u, padded.data(), &ops) ||
                std::memcmp(
                    padded.data(),
                    message.data() + (size_t)(K - 1u) * bb,
                    final_bytes) != 0 ||
                !std::all_of(
                    padded.begin() + final_bytes, padded.end(),
                    [](uint8_t value) { return value == 0u; }))
            {
                std::fprintf(stderr,
                    "borrowed lifetime: partial final padding was not zero\n");
                return false;
            }
        }

        // Release the caller allocation completely.  ASan will catch any
        // post-return read through a stale borrowed pointer; byte comparisons
        // also cover builds without sanitizers.
        std::vector<uint8_t>().swap(message);
        for (size_t i = 0; i <
                sizeof(packet_ids) / sizeof(packet_ids[0]); ++i)
        {
            std::vector<uint8_t> actual(bb, 0x5au);
            uint32_t data_bytes = UINT32_MAX;
            if (encoder.EncodeResult(
                    packet_ids[i], actual.data(), bb, &data_bytes) !=
                    Wirehair_Success ||
                data_bytes != expected_bytes[i] ||
                std::memcmp(
                    actual.data(), expected[i].data(), expected_bytes[i]) != 0 ||
                !std::all_of(
                    actual.begin() + expected_bytes[i], actual.end(),
                    [](uint8_t value) { return value == 0x5au; }))
            {
                std::fprintf(stderr,
                    "borrowed lifetime: encoder retained caller storage "
                    "for final=%u id=%u\n",
                    final_bytes, packet_ids[i]);
                return false;
            }
        }
    }

    std::printf("borrowed message lifetime: PASS\n");
    return true;
}

bool TestSystematicSourceCache()
{
    const uint32_t K = 64u;
    const uint32_t bb = 67u;
    const uint32_t tail = 13u;
    const uint64_t message_bytes = (uint64_t)(K - 1u) * bb + tail;
    std::vector<uint8_t> message((size_t)message_bytes);
    FillRandomBlocks(message.data(), message.size(), UINT64_C(0xca5ced));
    const std::vector<uint8_t> original = message;

    wirehair_v2::MessagePrecodeEncoder uncached;
    if (uncached.InitializeResult(
            message.data(), message_bytes, bb) != Wirehair_Success ||
        uncached.HasSystematicSourceCache() ||
        uncached.SystematicSourceCacheBytes() != 0u)
    {
        std::fprintf(stderr, "systematic cache: default path retained source\n");
        return false;
    }

    wirehair_v2::MessagePrecodeEncoderOptions cached_options;
    cached_options.CacheSystematicSource = true;
    wirehair_v2::MessagePrecodeEncoder cached;
    const wirehair_v2::SeedProfile selected_profile = uncached.Profile();
    if (cached.InitializeResult(
            message.data(), message_bytes, bb,
            &selected_profile, &cached_options) != Wirehair_Success ||
        !cached.Options().CacheSystematicSource ||
        !cached.HasSystematicSourceCache() ||
        cached.SystematicSourceCacheBytes() != message_bytes)
    {
        std::fprintf(stderr,
            "systematic cache: selected-profile opt-in failed\n");
        return false;
    }

    // A local storage policy must not alter any serialized matrix or packet
    // contract field selected by the uncached encoder.
    const wirehair_v2::SeedProfile& cached_profile = cached.Profile();
    if (cached_profile.BlockCount != selected_profile.BlockCount ||
        cached_profile.BlockBytes != selected_profile.BlockBytes ||
        cached_profile.DenseCount != selected_profile.DenseCount ||
        cached_profile.V2SeedAttempt != selected_profile.V2SeedAttempt ||
        cached_profile.V2PrecodeContractVersion !=
            selected_profile.V2PrecodeContractVersion ||
        cached_profile.V2PacketRowContractVersion !=
            selected_profile.V2PacketRowContractVersion ||
        cached_profile.V2StaircaseCount !=
            selected_profile.V2StaircaseCount ||
        cached_profile.V2DenseRowCount != selected_profile.V2DenseRowCount ||
        cached_profile.V2HeavyRowCount != selected_profile.V2HeavyRowCount ||
        cached_profile.V2CompletionField !=
            selected_profile.V2CompletionField ||
        cached_profile.V2SourceHits != selected_profile.V2SourceHits ||
        cached_profile.V2PrecodeSeed != selected_profile.V2PrecodeSeed ||
        cached_profile.V2PacketPeelSeed !=
            selected_profile.V2PacketPeelSeed ||
        cached_profile.V2RecoveryMixCount !=
            selected_profile.V2RecoveryMixCount ||
        cached_profile.V2DenseIdentityCorner !=
            selected_profile.V2DenseIdentityCorner ||
        cached_profile.V2PrecodeSeedSalt !=
            selected_profile.V2PrecodeSeedSalt ||
        cached_profile.V2RecoveryRowSeedSalt !=
            selected_profile.V2RecoveryRowSeedSalt)
    {
        std::fprintf(stderr,
            "systematic cache: local option changed selected profile\n");
        return false;
    }

    const uint32_t packet_ids[] = { 0u, K - 1u, K + 7u };
    std::vector<std::vector<uint8_t> > expected;
    std::vector<uint32_t> expected_bytes;
    for (const uint32_t packet_id : packet_ids)
    {
        expected.push_back(std::vector<uint8_t>(bb, 0xa5u));
        uint32_t bytes = UINT32_MAX;
        uint64_t uncached_ops = UINT64_MAX;
        if (uncached.EncodeResult(
                packet_id, expected.back().data(), bb,
                &bytes, &uncached_ops) != Wirehair_Success ||
            uncached_ops == 0u)
        {
            std::fprintf(stderr,
                "systematic cache: uncached reference encode failed\n");
            return false;
        }
        expected_bytes.push_back(bytes);

        std::vector<uint8_t> actual(bb, 0xa5u);
        uint32_t actual_bytes = UINT32_MAX;
        uint64_t cached_ops = UINT64_MAX;
        if (cached.EncodeResult(
                packet_id, actual.data(), bb,
                &actual_bytes, &cached_ops) != Wirehair_Success ||
            actual_bytes != bytes || actual != expected.back() ||
            (packet_id < K ? cached_ops != 0u : cached_ops != uncached_ops))
        {
            std::fprintf(stderr,
                "systematic cache: attached wire mismatch id=%u\n",
                packet_id);
            return false;
        }
    }

    // The cache is an owned byte-for-byte snapshot, not another borrowed
    // view of the caller allocation.
    std::fill(message.begin(), message.end(), uint8_t{0x5a});
    for (size_t i = 0u; i < 2u; ++i)
    {
        std::vector<uint8_t> actual(bb, 0xa5u);
        uint32_t bytes = UINT32_MAX;
        uint64_t ops = UINT64_MAX;
        if (cached.EncodeResult(
                packet_ids[i], actual.data(), bb, &bytes, &ops) !=
                Wirehair_Success ||
            bytes != expected_bytes[i] || actual != expected[i] || ops != 0u)
        {
            std::fprintf(stderr,
                "systematic cache: caller mutation changed cached bytes\n");
            return false;
        }
    }

    cached.ReleaseSystematicSourceCache();
    cached.ReleaseSystematicSourceCache();
    if (cached.HasSystematicSourceCache() ||
        cached.SystematicSourceCacheBytes() != 0u ||
        cached.Options().CacheSystematicSource)
    {
        std::fprintf(stderr, "systematic cache: release was not idempotent\n");
        return false;
    }
    for (size_t i = 0u; i < 3u; ++i)
    {
        std::vector<uint8_t> actual(bb, 0xa5u);
        uint32_t bytes = UINT32_MAX;
        uint64_t ops = UINT64_MAX;
        if (cached.EncodeResult(
                packet_ids[i], actual.data(), bb, &bytes, &ops) !=
                Wirehair_Success ||
            bytes != expected_bytes[i] || actual != expected[i] || ops == 0u)
        {
            std::fprintf(stderr,
                "systematic cache: detached fallback mismatch id=%u\n",
                packet_ids[i]);
            return false;
        }
    }

    // Failed cached initialization is transactional, including cache state.
    message = original;
    if (cached.InitializeResult(
            message.data(), message_bytes, bb,
            &selected_profile, &cached_options) != Wirehair_Success)
    {
        std::fprintf(stderr, "systematic cache: reattach initialize failed\n");
        return false;
    }
    const uint8_t* const working_intermediate = cached.IntermediateBlocks();
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult cache_oom = cached.InitializeResult(
        message.data(), message_bytes, bb,
        &selected_profile, &cached_options);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (cache_oom != Wirehair_OOM ||
        cached.IntermediateBlocks() != working_intermediate ||
        !cached.HasSystematicSourceCache() ||
        cached.SystematicSourceCacheBytes() != message_bytes)
    {
        std::fprintf(stderr,
            "systematic cache: OOM did not preserve active encoder\n");
        return false;
    }

    // A successful default reinitialization returns to x8rs.4's zero-copy
    // source-lifetime policy and releases the prior opt-in cache.
    if (cached.InitializeResult(
            message.data(), message_bytes, bb,
            &selected_profile, nullptr) != Wirehair_Success ||
        cached.HasSystematicSourceCache() ||
        cached.Options().CacheSystematicSource)
    {
        std::fprintf(stderr,
            "systematic cache: default reinitialize retained cache\n");
        return false;
    }

    wirehair_v2::Codec facade;
    if (facade.ReleasePrecodeEncoderSystematicCache() !=
            Wirehair_InvalidInput ||
        facade.InitializePrecodeEncoder(
            message.data(), message_bytes, bb,
            &selected_profile, &cached_options) != Wirehair_Success ||
        facade.ReleasePrecodeEncoderSystematicCache() != Wirehair_Success ||
        facade.ReleasePrecodeEncoderSystematicCache() != Wirehair_Success)
    {
        std::fprintf(stderr, "systematic cache: facade detach failed\n");
        return false;
    }
    std::vector<uint8_t> facade_block(bb, 0xa5u);
    uint32_t facade_bytes = UINT32_MAX;
    if (facade.Encode(
            K - 1u, facade_block.data(), bb, &facade_bytes) !=
            Wirehair_Success ||
        facade_bytes != tail ||
        std::memcmp(
            facade_block.data(), original.data() + (size_t)(K - 1u) * bb,
            tail) != 0)
    {
        std::fprintf(stderr,
            "systematic cache: facade fallback changed systematic bytes\n");
        return false;
    }

    std::printf("systematic source cache: PASS\n");
    return true;
}

bool TestTypedFailuresAndAllocationContainment()
{
    const uint32_t K = 16u;
    const uint32_t bb = 19u;
    const uint64_t message_bytes = (uint64_t)K * bb;
    std::vector<uint8_t> message((size_t)message_bytes);
    FillRandomBlocks(message.data(), message.size(), UINT64_C(0xa110ca7e));

    wirehair_v2::MessagePrecodeEncoderOptions identity;
    identity.DenseIdentityCorner = true;
    wirehair_v2::MessagePrecodeEncoder encoder;
    if (encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &identity) !=
            Wirehair_Success)
    {
        std::fprintf(stderr, "typed failures: initial encoder failed\n");
        return false;
    }

    const wirehair_v2::SeedProfile working_profile = encoder.Profile();
    const uint8_t* working_source = encoder.IntermediateBlocks();
    const uint8_t dummy = 0u;
    if (encoder.InitializeResult(
            &dummy, UINT64_C(0x100000000), UINT32_C(0x80000000),
            nullptr, &identity) != Wirehair_InvalidInput ||
        !encoder.IsInitialized() ||
        encoder.IntermediateBlocks() != working_source ||
        encoder.Profile().BlockCount != working_profile.BlockCount)
    {
        std::fprintf(stderr,
            "typed failures: invalid input classification/state failed\n");
        return false;
    }

    wirehair_v2::MessagePrecodeEncoderOptions certified;
    if (encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &certified) !=
            Wirehair_Success || !encoder.IsInitialized() ||
        encoder.Options().DenseIdentityCorner)
    {
        std::fprintf(stderr,
            "typed failures: certified global solve failed\n");
        return false;
    }
    working_source = encoder.IntermediateBlocks();

    // The solved encoder no longer copies the nested precode row graph, so
    // only the three real guarded allocations before ownership transfer remain.
    for (int64_t failure = 0; failure < 3; ++failure)
    {
        wirehair_v2::SetAllocationFailureCountdownForTesting(failure);
        const WirehairResult init_oom = encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &identity);
        wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
        if (init_oom != Wirehair_OOM || !encoder.IsInitialized() ||
            encoder.IntermediateBlocks() != working_source)
        {
            std::fprintf(stderr,
                "typed failures: init OOM containment failed at %lld\n",
                (long long)failure);
            return false;
        }
    }

    // A partial input owns exactly one padded tail during the synchronous
    // solve.  Failure of that allocation must preserve the prior encoder.
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult tail_oom = encoder.InitializeResult(
        message.data(), message_bytes - 1u, bb, nullptr, &identity);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (tail_oom != Wirehair_OOM || !encoder.IsInitialized() ||
        encoder.IntermediateBlocks() != working_source)
    {
        std::fprintf(stderr,
            "typed failures: padded-tail OOM containment failed\n");
        return false;
    }

    std::vector<uint8_t> output(bb, 0xa5u);
    const std::vector<uint8_t> before = output;
    uint32_t data_bytes = UINT32_MAX;
    uint64_t ops = UINT64_MAX;
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult encode_oom = encoder.EncodeResult(
        K, output.data(), bb, &data_bytes, &ops);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (encode_oom != Wirehair_OOM || output != before ||
        data_bytes != UINT32_MAX || ops != UINT64_MAX)
    {
        std::fprintf(stderr,
            "typed failures: recovery Encode OOM modified output\n");
        return false;
    }
    if (encoder.EncodeResult(
            K, output.data(), bb, &data_bytes, &ops) != Wirehair_Success ||
        data_bytes != bb || ops == 0u)
    {
        std::fprintf(stderr, "typed failures: post-OOM encode failed\n");
        return false;
    }

    wirehair_v2::Codec facade;
    if (facade.InitializePrecodeEncoder(
            message.data(), message_bytes, bb, nullptr, &identity) !=
            Wirehair_Success)
    {
        std::fprintf(stderr, "typed failures: facade initial encode failed\n");
        return false;
    }
    if (facade.InitializePrecodeEncoder(
            &dummy, UINT64_C(0x100000000), UINT32_C(0x80000000),
            nullptr, &identity) != Wirehair_InvalidInput ||
        facade.Profile().BlockCount != K ||
        facade.InitializePrecodeEncoder(
            message.data(), message_bytes, bb, nullptr, &certified) !=
            Wirehair_Success ||
        facade.Profile().BlockCount != K)
    {
        std::fprintf(stderr,
            "typed failures: facade invalid/singular classification failed\n");
        return false;
    }
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult facade_oom = facade.InitializePrecodeEncoder(
        message.data(), message_bytes, bb, nullptr, &identity);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (facade_oom != Wirehair_OOM ||
        facade.Profile().BlockCount != K)
    {
        std::fprintf(stderr,
            "typed failures: facade OOM classification/state failed\n");
        return false;
    }

    std::fill(output.begin(), output.end(), uint8_t{0x5a});
    const std::vector<uint8_t> facade_before = output;
    data_bytes = UINT32_MAX;
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult facade_encode_oom = facade.Encode(
        K, output.data(), bb, &data_bytes);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (facade_encode_oom != Wirehair_OOM || output != facade_before ||
        data_bytes != UINT32_MAX)
    {
        std::fprintf(stderr,
            "typed failures: facade Encode OOM modified output\n");
        return false;
    }

    std::printf("typed failures and OOM containment: PASS\n");
    return true;
}

int RunMixedBucketEncodeBenchmark(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t iterations,
    wirehair_v2::MixedResidueBucketMode mode)
{
    if (K < 2u || K > 64000u || block_bytes == 0u ||
        (block_bytes & 1u) != 0u || iterations == 0u ||
        (uint64_t)K * block_bytes >
            (uint64_t)std::numeric_limits<size_t>::max())
    {
        return 2;
    }
    if (!wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16RowsMax) ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(32u) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256Rows + 1u) ||
        !wirehair_v2::SetMixedResidueSkewForTesting(0u) ||
        !wirehair_v2::SetMixedResidueScheduleForTesting(
            wirehair_v2::MixedResidueSchedule::Hashed))
    {
        return 2;
    }
    wirehair_v2::SetMixedResidueHashSeedForTesting(68u);
    if (!wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(true) ||
        !wirehair_v2::SetMixedResidueBucketModeForTesting(mode))
    {
        return 2;
    }

    wirehair_v2::PrecodeParams params = wirehair_v2::MakeMixedParams(
        K, UINT64_C(0x6a6f696e7464656c) ^ K);
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        return 1;
    }
    const uint32_t parity_count =
        params.Staircase + params.DenseRows + params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * block_bytes);
    std::vector<uint8_t> parity((size_t)parity_count * block_bytes);
    FillRandomBlocks(
        source.data(), source.size(),
        UINT64_C(0x656e636f6465626d) ^ K ^ block_bytes);
    wirehair_v2::PrecodeEncodeStats stats;
    if (!wirehair_v2::ComputePrecodeValues(
            system, source.data(), block_bytes, parity.data(), &stats))
    {
        return 1;
    }
    const auto start = std::chrono::steady_clock::now();
    for (uint32_t iteration = 0u; iteration < iterations; ++iteration)
    {
        if (!wirehair_v2::ComputePrecodeValues(
                system, source.data(), block_bytes, parity.data(), &stats))
        {
            return 1;
        }
    }
    const auto stop = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(stop - start).count();
    uint64_t checksum = 0u;
    for (size_t i = 0u; i < parity.size(); i += 4096u) {
        checksum = checksum * 257u + parity[i];
    }
    std::printf(
        "mixed_bucket_encode K=%u bb=%u iterations=%u mode=%u "
        "elapsed_ms=%.6f per_iteration_ms=%.6f bucket_xors=%llu "
        "joint_source=%llu joint_marginal=%llu joint_copies=%llu "
        "dual_source_columns=%llu deltas=%u scratch=%llu "
        "checksum=0x%llx\n",
        K, block_bytes, iterations, (uint32_t)mode,
        elapsed_ms, elapsed_ms / iterations,
        (unsigned long long)stats.HeavyBucketXors,
        (unsigned long long)stats.MixedJointSourceXors,
        (unsigned long long)stats.MixedJointMarginalXors,
        (unsigned long long)stats.MixedJointMarginalCopies,
        (unsigned long long)stats.MixedDualSourceColumns,
        stats.MixedJointActiveDeltas,
        (unsigned long long)stats.MixedJointScratchBytes,
        (unsigned long long)checksum);
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    gf256_init();

    if (argc == 6 &&
        std::strcmp(argv[1], "--mixed-bucket-bench") == 0)
    {
        uint32_t K = 0u;
        uint32_t block_bytes = 0u;
        uint32_t iterations = 0u;
        if (!ParsePositiveU32(argv[2], K) ||
            !ParsePositiveU32(argv[3], block_bytes) ||
            !ParsePositiveU32(argv[4], iterations))
        {
            return 2;
        }
        wirehair_v2::MixedResidueBucketMode mode;
        if (std::strcmp(argv[5], "auto") == 0) {
            mode = wirehair_v2::MixedResidueBucketMode::Automatic;
        }
        else if (std::strcmp(argv[5], "separate") == 0) {
            mode = wirehair_v2::MixedResidueBucketMode::Separate;
        }
        else if (std::strcmp(argv[5], "dual") == 0) {
            mode = wirehair_v2::MixedResidueBucketMode::Dual;
        }
        else if (std::strcmp(argv[5], "joint-delta") == 0) {
            mode = wirehair_v2::MixedResidueBucketMode::JointDelta;
        }
        else {
            return 2;
        }
        return RunMixedBucketEncodeBenchmark(
            K, block_bytes, iterations, mode);
    }

    if (argc == 2 && std::strcmp(argv[1], "--oom-probe") == 0)
    {
        const uint8_t dummy = 0u;
        wirehair_v2::MessagePrecodeEncoderOptions options;
        options.DenseIdentityCorner = true;
        wirehair_v2::Codec codec;
        const WirehairResult result = codec.InitializePrecodeEncoder(
            &dummy,
            UINT64_C(512) * 1024u * 1024u,
            UINT32_C(8) * 1024u * 1024u,
            nullptr,
            &options);
        if (result != Wirehair_OOM) {
            std::fprintf(stderr,
                "constrained precode init returned %d, expected OOM\n",
                (int)result);
            return 1;
        }
        std::printf("constrained precode init: Wirehair_OOM\n");
        return 0;
    }

    // Feasibility trials tunable so the sanitizer build stays quick
    uint32_t trials = 2000u;
    if (argc > 2 || (argc == 2 && !ParsePositiveU32(argv[1], trials))) {
        std::fprintf(stderr,
            "usage: %s [positive feasibility-trials]\n", argv[0]);
        return 1;
    }

    bool ok = true;
    ok = TestCorrectnessAsBuilt() && ok;
    ok = TestCorrectnessDoctored() && ok;
    ok = TestTwoAnchorPhasedEncoder() && ok;
    ok = TestMalformedDenseCorner() && ok;
    ok = TestStrictSystemValidation() && ok;
    ok = TestCostModel() && ok;
    ok = TestHeavyResidueDispatch() && ok;
    ok = TestMixedCornerRank() && ok;
    ok = TestIndependentMixedCornerAcrossK() && ok;
    ok = TestMixedCompletion() && ok;
    ok = TestAutomaticMixedJointBucketPolicy() && ok;
    ok = TestRecoveryBlockEncoding() && ok;
    ok = TestMessagePrecodeEncoder() && ok;
    ok = TestBorrowedMessageLifetime() && ok;
    ok = TestSystematicSourceCache() && ok;
    ok = TestTypedFailuresAndAllocationContainment() && ok;
    ok = TestFeasibility(trials) && ok;

    if (!ok) {
        std::fprintf(stderr, "precode_encode_test: FAIL\n");
        return 1;
    }
    std::printf("\nprecode_encode_test: PASS\n");
    return 0;
}
