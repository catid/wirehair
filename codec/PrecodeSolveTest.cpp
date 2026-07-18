#include "WirehairV2Plan.h"
#include "WirehairV2Precode.h"
#include "WirehairV2PrecodeDecode.h"
#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"
#include "V2TinyDenseOracle.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <atomic>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cstring>
#include <thread>
#include <vector>

namespace {

class MixedProjectionOracleScope
{
public:
    MixedProjectionOracleScope()
    {
        wirehair_v2::SetMixedProjectionOracleForTesting(true);
    }

    ~MixedProjectionOracleScope()
    {
        wirehair_v2::SetMixedProjectionOracleForTesting(false);
    }
};

class SolveValueArenaPoisonScope
{
public:
    SolveValueArenaPoisonScope()
    {
        wirehair_v2::SetSolveValueArenaPoisonForTesting(true);
    }

    ~SolveValueArenaPoisonScope()
    {
        wirehair_v2::SetSolveValueArenaPoisonForTesting(false);
    }
};

class FusedBlockInitializationScope
{
public:
    FusedBlockInitializationScope()
    {
        wirehair_v2::SetFusedBlockInitializationForTesting(true);
    }

    ~FusedBlockInitializationScope()
    {
        wirehair_v2::SetFusedBlockInitializationForTesting(false);
    }
};

class SolveValueArenaAllocationFailureScope
{
public:
    SolveValueArenaAllocationFailureScope()
    {
        wirehair_v2::SetSolveValueArenaAllocationFailureForTesting(true);
    }

    ~SolveValueArenaAllocationFailureScope()
    {
        wirehair_v2::SetSolveValueArenaAllocationFailureForTesting(false);
    }
};

bool SameSolveStats(
    const wirehair_v2::PrecodeSolveStats& a,
    const wirehair_v2::PrecodeSolveStats& b)
{
    return a.PacketRows == b.PacketRows &&
        a.PeeledColumns == b.PeeledColumns &&
        a.InactivatedColumns == b.InactivatedColumns &&
        a.ResidualRows == b.ResidualRows &&
        a.ResidualRank == b.ResidualRank &&
        a.BinaryResidualRank == b.BinaryResidualRank &&
        a.BinaryRowReferences == b.BinaryRowReferences &&
        a.BinaryRowStorageBytes == b.BinaryRowStorageBytes &&
        a.BinaryAdjacencyStorageBytes == b.BinaryAdjacencyStorageBytes &&
        a.BinaryRowStorageAllocations == b.BinaryRowStorageAllocations &&
        a.BinaryAdjacencyStorageAllocations ==
            b.BinaryAdjacencyStorageAllocations &&
        a.BlockXors == b.BlockXors &&
        a.BlockMulAdds == b.BlockMulAdds &&
        a.BuildNanoseconds == b.BuildNanoseconds &&
        a.PeelNanoseconds == b.PeelNanoseconds &&
        a.ProjectNanoseconds == b.ProjectNanoseconds &&
        a.ResidualNanoseconds == b.ResidualNanoseconds &&
        a.BackSubNanoseconds == b.BackSubNanoseconds &&
        a.PacketSeedAttempt == b.PacketSeedAttempt &&
        a.MixedJointSourceXors == b.MixedJointSourceXors &&
        a.MixedJointMarginalXors == b.MixedJointMarginalXors &&
        a.MixedJointMarginalCopies == b.MixedJointMarginalCopies &&
        a.MixedJointScratchBytes == b.MixedJointScratchBytes &&
        a.MixedJointActiveDeltas == b.MixedJointActiveDeltas &&
        a.MixedDualSourceColumns == b.MixedDualSourceColumns &&
        a.SolveValueArenaBytes == b.SolveValueArenaBytes &&
        a.SolveValueArenaEagerZeroBytes ==
            b.SolveValueArenaEagerZeroBytes &&
        a.SolveValueArenaCommitCopyBytes ==
            b.SolveValueArenaCommitCopyBytes;
}

bool SameResumeState(
    const wirehair_v2::PrecodeSolveResumeState& a,
    const wirehair_v2::PrecodeSolveResumeState& b)
{
    return a.SourceCount == b.SourceCount &&
        a.PrecodeCount == b.PrecodeCount &&
        a.ColumnCount == b.ColumnCount &&
        a.BlockBytes == b.BlockBytes &&
        a.InactiveCount == b.InactiveCount &&
        a.ProjectionWords == b.ProjectionWords &&
        a.Rank == b.Rank &&
        a.Config.PeelSeed == b.Config.PeelSeed &&
        a.Config.MixCount == b.Config.MixCount &&
        a.Runtime.SourcePrime() == b.Runtime.SourcePrime() &&
        a.Runtime.PrecodePrime() == b.Runtime.PrecodePrime() &&
        SameSolveStats(a.Stats, b.Stats) &&
        a.InactiveIndex == b.InactiveIndex &&
        a.InactiveColumns == b.InactiveColumns &&
        a.Projection == b.Projection &&
        a.Values == b.Values &&
        a.PivotCoefficients == b.PivotCoefficients &&
        a.PivotRhs == b.PivotRhs &&
        a.HavePivot == b.HavePivot &&
        a.CoefficientScratch == b.CoefficientScratch &&
        a.RhsScratch == b.RhsScratch &&
        a.Active == b.Active;
}

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

class MixedGroupedGF256RowsScope
{
public:
    explicit MixedGroupedGF256RowsScope(uint32_t rows)
        : Previous(wirehair_v2::ActiveMixedGroupedGF256Rows())
        , Valid(wirehair_v2::SetMixedGroupedGF256RowsForTesting(rows))
    {
    }

    ~MixedGroupedGF256RowsScope()
    {
        (void)wirehair_v2::SetMixedGroupedGF256RowsForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    uint32_t Previous;
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

class BinaryPeelOracleScope
{
public:
    BinaryPeelOracleScope()
    {
        wirehair_v2::SetBinaryPeelOracleForTesting(true);
    }

    ~BinaryPeelOracleScope()
    {
        wirehair_v2::SetBinaryPeelOracleForTesting(false);
    }
};

bool CheckLowestBitIndex()
{
    for (unsigned bit = 0u; bit < 64u; ++bit)
    {
        const uint64_t word = UINT64_C(1) << bit;
        if (wirehair::NonzeroLowestBitIndex64(word) != bit)
        {
            std::fprintf(stderr,
                "solve: lowest-bit singleton mismatch bit=%u\n", bit);
            return false;
        }
    }
    struct Pattern
    {
        uint64_t Word;
        unsigned Expected;
    };
    static const Pattern kPatterns[] = {
        { UINT64_MAX, 0u },
        { UINT64_C(0x8000000080000000), 31u },
        { UINT64_C(0xffffffff00000000), 32u },
        { UINT64_C(0xc000000000000000), 62u }
    };
    for (const Pattern& pattern : kPatterns)
    {
        if (wirehair::NonzeroLowestBitIndex64(pattern.Word) !=
            pattern.Expected)
        {
            std::fprintf(stderr,
                "solve: lowest-bit pattern mismatch word=%016llx\n",
                (unsigned long long)pattern.Word);
            return false;
        }
    }
    std::printf("portable lowest-bit boundaries: PASS\n");
    return true;
}

std::vector<uint32_t> ReferencePacketRow(
    uint32_t K,
    uint32_t P,
    uint32_t block_id,
    const wirehair_v2::PacketRowConfig& config)
{
    wirehair::PeelRowParameters params;
    params.Initialize(
        block_id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
    std::vector<uint32_t> row;
    row.reserve((size_t)params.PeelCount + config.MixCount);
    wirehair::PeelRowIterator source(
        params, (uint16_t)K, wirehair::NextPrime16((uint16_t)K));
    do {
        row.push_back(source.GetColumn());
    } while (source.Iterate());
    const wirehair::RowMixIterator mix(
        params, (uint16_t)P, wirehair::NextPrime16((uint16_t)P));
    for (uint32_t i = 0u; i < config.MixCount; ++i) {
        row.push_back(K + mix.Columns[i]);
    }
    return row;
}

std::vector<uint8_t> ReferencePacket(
    uint32_t K,
    uint32_t P,
    uint32_t block_id,
    const wirehair_v2::PacketRowConfig& config,
    const uint8_t* intermediate,
    uint32_t block_bytes)
{
    std::vector<uint8_t> expected(block_bytes, 0u);
    const std::vector<uint32_t> row =
        ReferencePacketRow(K, P, block_id, config);
    for (uint32_t column : row) {
        const uint8_t* source =
            intermediate + (size_t)column * block_bytes;
        for (uint32_t i = 0; i < block_bytes; ++i) {
            expected[i] ^= source[i];
        }
    }
    return expected;
}

bool CheckPacketEvaluationCase(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_id,
    uint32_t block_bytes,
    unsigned input_offset,
    unsigned output_offset,
    bool fast_path)
{
    static const size_t kGuardBytes = 64u;
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const size_t intermediate_bytes = (size_t)(K + P) * block_bytes;
    std::vector<uint8_t> input_storage(
        intermediate_bytes + 2u * kGuardBytes + 128u, 0xc7u);
    const uintptr_t input_start = reinterpret_cast<uintptr_t>(
        input_storage.data() + kGuardBytes);
    uint8_t* intermediate = reinterpret_cast<uint8_t*>(
        (input_start + 63u) & ~uintptr_t(63u)) + input_offset;
    for (size_t i = 0; i < intermediate_bytes; ++i) {
        intermediate[i] = (uint8_t)(
            i * 131u + (i >> 7) + block_id * 17u + block_bytes);
    }
    const std::vector<uint8_t> input_before = input_storage;

    std::vector<uint8_t> output_storage(
        block_bytes + 2u * kGuardBytes + 128u, 0xa5u);
    const uintptr_t output_start = reinterpret_cast<uintptr_t>(
        output_storage.data() + kGuardBytes);
    uint8_t* output = reinterpret_cast<uint8_t*>(
        (output_start + 63u) & ~uintptr_t(63u)) + output_offset;
    std::vector<uint8_t> expected_storage = output_storage;
    const std::vector<uint8_t> expected = ReferencePacket(
        K, P, block_id, config, intermediate, block_bytes);
    const std::vector<uint32_t> expected_row =
        ReferencePacketRow(K, P, block_id, config);
    const std::vector<uint32_t> generated_row =
        wirehair_v2::GeneratePacketMatrixRow(K, P, block_id, config);
    if (generated_row != expected_row)
    {
        std::fprintf(stderr,
            "solve: packet row/reference mismatch id=%u mix=%u\n",
            block_id, config.MixCount);
        return false;
    }
    std::memcpy(
        expected_storage.data() + (output - output_storage.data()),
        expected.data(),
        block_bytes);

    uint64_t operations = UINT64_C(0xfedcba9876543210);
    const bool evaluated = fast_path ?
        wirehair_v2::EvaluatePacketBlockForValidatedSystem(
            system, config, intermediate, block_bytes, block_id,
            output, &operations) :
        wirehair_v2::EvaluatePacketBlock(
            system, config, intermediate, block_bytes, block_id,
            output, &operations);
    const size_t row_size = expected_row.size();
    if (!evaluated || operations != row_size ||
        output_storage != expected_storage || input_storage != input_before)
    {
        std::fprintf(stderr,
            "solve: packet evaluation mismatch id=%u bb=%u inoff=%u "
            "outoff=%u mix=%u fast=%u operations=%" PRIu64
            " expected=%zu\n",
            block_id, block_bytes, input_offset, output_offset,
            config.MixCount, fast_path ? 1u : 0u,
            operations, row_size);
        return false;
    }

    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }
    std::fill(output_storage.begin(), output_storage.end(), uint8_t{0xa5u});
    operations = UINT64_C(0xfedcba9876543210);
    if (!wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
            system, config, runtime, intermediate, block_bytes, block_id,
            output, &operations) ||
        operations != row_size || output_storage != expected_storage ||
        input_storage != input_before)
    {
        std::fprintf(stderr,
            "solve: cached packet evaluation mismatch id=%u bb=%u "
            "mix=%u operations=%" PRIu64 " expected=%zu\n",
            block_id, block_bytes, config.MixCount, operations, row_size);
        return false;
    }
    return true;
}

bool CheckPacketEvaluationFusion()
{
    const uint32_t K = 128u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(
                K, UINT64_C(0x786f72667573696f)),
            system))
    {
        return false;
    }
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = UINT32_C(0x4d241359);
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;

    // Pin the source-weight transitions used by packet-tail pairing: exact
    // degrees one through five, a mid-weight row, and the capped heavy tail.
    // A high-bit public id separately guards the full repair-id domain.
    static const unsigned kWeightFixtureCount = 7u;
    static const unsigned kFixtureCount = kWeightFixtureCount + 1u;
    uint32_t ids[kFixtureCount] = {
        UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
        UINT32_MAX, UINT32_MAX,
        UINT32_C(0xf1234567)
    };
    for (uint32_t id = 0u; id < 1000000u; ++id)
    {
        const size_t row_size = wirehair_v2::GeneratePacketMatrixRow(
            K, P, id, config).size();
        const uint32_t degree =
            (uint32_t)row_size - config.MixCount;
        unsigned slot = kWeightFixtureCount;
        if (degree == 1u) slot = 0u;
        else if (degree == 2u) slot = 1u;
        else if (degree == 3u) slot = 2u;
        else if (degree == 4u) slot = 3u;
        else if (degree == 5u) slot = 4u;
        else if (degree >= 8u && degree <= 16u) slot = 5u;
        else if (degree >= 32u) slot = 6u;
        if (slot < kWeightFixtureCount && ids[slot] == UINT32_MAX) {
            ids[slot] = id;
        }
        if (ids[0] != UINT32_MAX && ids[1] != UINT32_MAX &&
            ids[2] != UINT32_MAX && ids[3] != UINT32_MAX &&
            ids[4] != UINT32_MAX && ids[5] != UINT32_MAX &&
            ids[6] != UINT32_MAX)
        {
            break;
        }
    }
    for (unsigned i = 0; i < kWeightFixtureCount; ++i) {
        if (ids[i] == UINT32_MAX) {
            std::fprintf(stderr,
                "solve: packet evaluation source-weight fixture %u missing\n",
                i);
            return false;
        }
    }

    static const uint32_t kLengths[] = {
        1u, 2u, 3u, 7u, 15u, 16u, 17u, 31u, 32u, 33u,
        63u, 64u, 65u, 127u, 128u, 129u, 255u, 256u, 257u, 1280u,
        32767u, 32768u, 32769u
    };
    static const unsigned kOffsets[] = { 0u, 1u, 7u, 15u, 31u, 63u };
    static const uint32_t kFusedMixCounts[] = {
        wirehair_v2::kCertifiedPacketMixCount, 2u, 1u
    };
    for (uint32_t mix_count : kFusedMixCounts)
    {
        wirehair_v2::PacketRowConfig fused = config;
        fused.MixCount = mix_count;
        for (unsigned weight_i = 0; weight_i < kFixtureCount; ++weight_i)
        {
            for (unsigned length_i = 0;
                 length_i < sizeof(kLengths) / sizeof(kLengths[0]);
                 ++length_i)
            {
                const unsigned offset_i =
                    (weight_i * 3u + length_i) %
                    (sizeof(kOffsets) / sizeof(kOffsets[0]));
                if (!CheckPacketEvaluationCase(
                        system, fused, ids[weight_i], kLengths[length_i],
                        kOffsets[offset_i],
                        kOffsets[(offset_i * 5u + 1u) %
                            (sizeof(kOffsets) / sizeof(kOffsets[0]))],
                        ((weight_i + length_i) & 1u) != 0u))
                {
                    return false;
                }
            }
        }
    }

    // K=64000 removes singleton rows and exercises a very different modulus
    // and working-set shape.  Check both sides of each mix-count crossover:
    // total row weight six means d3/m3, d4/m2, and d5/m1 respectively.
    {
        static const uint32_t large_K = 64000u;
        wirehair_v2::PrecodeSystem large_system;
        if (!wirehair_v2::BuildPrecodeSystem(
                wirehair_v2::MakeCertifiedParams(
                    large_K, UINT64_C(0x6b36346b7061636b)),
                large_system))
        {
            return false;
        }
        const uint32_t large_P = large_system.Params.Staircase +
            large_system.Params.DenseRows + large_system.Params.HeavyRows;
        uint32_t degree_ids[5] = {
            UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX
        };
        for (uint32_t id = 0u; id < 1000000u; ++id)
        {
            wirehair::PeelRowParameters params;
            params.Initialize(
                id, config.PeelSeed, (uint16_t)large_K, (uint16_t)large_P);
            if (params.PeelCount >= 2u && params.PeelCount <= 6u) {
                degree_ids[params.PeelCount - 2u] = id;
            }
            bool complete = true;
            for (uint32_t fixture : degree_ids) {
                complete = complete && fixture != UINT32_MAX;
            }
            if (complete) {
                break;
            }
        }
        for (uint32_t fixture : degree_ids) {
            if (fixture == UINT32_MAX) {
                std::fprintf(stderr,
                    "solve: K64k packet degree fixture missing\n");
                return false;
            }
        }
        for (uint32_t mix_count : kFusedMixCounts)
        {
            wirehair_v2::PacketRowConfig fused = config;
            fused.MixCount = mix_count;
            for (unsigned i = 0u; i < 5u; ++i)
            {
                if (!CheckPacketEvaluationCase(
                        large_system, fused, degree_ids[i], 33u,
                        (i * 7u + mix_count) & 31u,
                        (i * 11u + mix_count * 3u) & 31u,
                        ((i + mix_count) & 1u) != 0u))
                {
                    return false;
                }
            }
        }
    }

    // Exercise the one-pass packet set-XOR gate at its exact 16-term boundary.
    // The existing length/offset matrix above covers every tail shape on the
    // paired fallback; this larger-K case selects the fixed-count family
    // with independently unaligned input and output ranges.
    {
        static const uint32_t set_xor_K = 10000u;
        static const uint32_t set_xor_block_bytes = 1280u;
        wirehair_v2::PrecodeSystem set_xor_system;
        if (!wirehair_v2::BuildPrecodeSystem(
                wirehair_v2::MakeCertifiedParams(
                    set_xor_K, UINT64_C(0x736574786f723136)),
                set_xor_system))
        {
            return false;
        }
        const uint32_t set_xor_P =
            set_xor_system.Params.Staircase +
            set_xor_system.Params.DenseRows +
            set_xor_system.Params.HeavyRows;
        wirehair_v2::PacketRowConfig set_xor_config = config;
        set_xor_config.MixCount = 2u;
        uint32_t set_xor_id = UINT32_MAX;
        for (uint32_t id = 0u; id < 1000000u; ++id)
        {
            wirehair::PeelRowParameters params;
            params.Initialize(
                id, set_xor_config.PeelSeed,
                (uint16_t)set_xor_K, (uint16_t)set_xor_P);
            if (params.PeelCount + set_xor_config.MixCount == 16u)
            {
                set_xor_id = id;
                break;
            }
        }
        if (set_xor_id == UINT32_MAX ||
            !CheckPacketEvaluationCase(
                set_xor_system, set_xor_config, set_xor_id,
                set_xor_block_bytes, 7u, 31u, true))
        {
            std::fprintf(stderr,
                "solve: 16-term packet set-XOR fixture failed\n");
            return false;
        }
    }

    // Hard-coded packet bytes guard the shipping equation and fused schedule
    // independently of the row-column golden vectors in PolicyTest.
    {
        static const uint32_t block_bytes = 32u;
        const size_t intermediate_bytes = (size_t)(K + P) * block_bytes;
        std::vector<uint8_t> intermediate(intermediate_bytes);
        for (size_t i = 0; i < intermediate.size(); ++i) {
            intermediate[i] = (uint8_t)(i * 131u + (i >> 7) + 0x5du);
        }
        uint8_t output[block_bytes] = {};
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, intermediate.data(), block_bytes,
                UINT32_C(0xf1234567), output))
        {
            return false;
        }
        static const uint8_t expected[block_bytes] = {
            0xbd, 0xd8, 0xbb, 0x46, 0x29, 0xc4, 0xa7, 0x32,
            0x95, 0x70, 0xf3, 0x5e, 0xc1, 0x5c, 0xff, 0x4a,
            0xed, 0x48, 0xab, 0x36, 0x99, 0x74, 0x97, 0x62,
            0x05, 0xe0, 0x63, 0x0e, 0xf1, 0x4c, 0x2f, 0xba
        };
        if (std::memcmp(output, expected, block_bytes) != 0)
        {
            std::fprintf(stderr, "solve: packet byte golden changed:");
            for (uint8_t byte : output) {
                std::fprintf(stderr, " %02x", (unsigned)byte);
            }
            std::fprintf(stderr, "\n");
            return false;
        }
    }

    // The API rejects every overlap shape before touching output/work count.
    {
        static const uint32_t block_bytes = 17u;
        const size_t intermediate_bytes = (size_t)(K + P) * block_bytes;
        std::vector<uint8_t> storage(
            intermediate_bytes + 2u * block_bytes + 2u, 0x6bu);
        uint8_t* intermediate = storage.data() + block_bytes + 1u;
        for (size_t i = 0; i < intermediate_bytes; ++i) {
            intermediate[i] = (uint8_t)(i * 29u + 7u);
        }
        uint8_t* overlap_outputs[] = {
            intermediate,
            intermediate - 1u,
            intermediate + intermediate_bytes - 1u
        };
        for (uint32_t mix_count : kFusedMixCounts)
        {
            wirehair_v2::PacketRowConfig fused = config;
            fused.MixCount = mix_count;
            for (unsigned i = 0; i < 3u; ++i)
            {
                const std::vector<uint8_t> before = storage;
                uint64_t operations = UINT64_C(0x0123456789abcdef);
                const bool evaluated = i == 1u ?
                    wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                        system, fused, intermediate, block_bytes, 17u,
                        overlap_outputs[i], &operations) :
                    wirehair_v2::EvaluatePacketBlock(
                        system, fused, intermediate, block_bytes, 17u,
                        overlap_outputs[i], &operations);
                if (evaluated || storage != before ||
                    operations != UINT64_C(0x0123456789abcdef))
                {
                    std::fprintf(stderr,
                        "solve: packet overlap shape %u mix=%u was not "
                        "no-write\n",
                        i, mix_count);
                    return false;
                }
            }

            // Exact endpoint adjacency is non-overlap and remains supported.
            uint8_t* adjacent_output = intermediate + intermediate_bytes;
            const std::vector<uint8_t> expected = ReferencePacket(
                K, P, 17u, fused, intermediate, block_bytes);
            if (!wirehair_v2::EvaluatePacketBlock(
                    system, fused, intermediate, block_bytes, 17u,
                    adjacent_output) ||
                std::memcmp(
                    adjacent_output, expected.data(), block_bytes) != 0)
            {
                std::fprintf(stderr,
                    "solve: adjacent packet output failed mix=%u\n",
                    mix_count);
                return false;
            }
        }
    }

    // A cache is accepted only for the exact immutable packet domain.  Failed
    // initialization also clears an earlier valid cache so stale prime values
    // can never influence equations.
    {
        wirehair_v2::PacketRowRuntime runtime;
        if (!runtime.Initialize(K, P, config.MixCount) ||
            runtime.SourcePrime() != wirehair::NextPrime16((uint16_t)K) ||
            runtime.PrecodePrime() != wirehair::NextPrime16((uint16_t)P))
        {
            std::fprintf(stderr, "solve: packet runtime initialization failed\n");
            return false;
        }
        static const uint32_t kIds[] = {
            0u, K - 1u, K, UINT32_C(0xf1234567), UINT32_MAX
        };
        for (uint32_t id : kIds)
        {
            if (wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                    K, P, id, config, runtime) !=
                wirehair_v2::GeneratePacketMatrixRow(K, P, id, config))
            {
                std::fprintf(stderr,
                    "solve: cached packet row mismatch id=%u\n", id);
                return false;
            }
        }

        wirehair_v2::PacketRowRuntime stale;
        if (!stale.Initialize(K + 1u, P, config.MixCount)) {
            return false;
        }
        static const uint32_t block_bytes = 17u;
        std::vector<uint8_t> intermediate(
            (size_t)(K + P) * block_bytes, 0x3cu);
        std::vector<uint8_t> output(block_bytes, 0xa5u);
        const std::vector<uint8_t> output_before = output;
        uint64_t operations = UINT64_C(0x0123456789abcdef);
        if (!wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                K, P, 17u, config, stale).empty() ||
            wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
                system, config, stale, intermediate.data(), block_bytes,
                17u, output.data(), &operations) ||
            output != output_before ||
            operations != UINT64_C(0x0123456789abcdef))
        {
            std::fprintf(stderr,
                "solve: stale packet runtime was not rejected/no-write\n");
            return false;
        }
        std::vector<wirehair_v2::SolvePacket> packets(K);
        for (uint32_t id = 0u; id < K; ++id) {
            packets[id].BlockId = id;
            packets[id].Data = intermediate.data() + (size_t)id * block_bytes;
        }
        std::vector<uint8_t> solved(11u, 0x6du);
        const std::vector<uint8_t> solved_before = solved;
        wirehair_v2::PrecodeSolveStats stats;
        stats.PacketRows = UINT32_C(0x76543210);
        if (wirehair_v2::SolvePrecodeSystemWithRuntime(
                system, config, stale, packets, block_bytes,
                solved, &stats) != Wirehair_InvalidInput ||
            solved != solved_before || stats.PacketRows != UINT32_C(0x76543210))
        {
            std::fprintf(stderr,
                "solve: stale runtime solve was not invalid/no-write\n");
            return false;
        }
        wirehair_v2::PrecodeSystem malformed = system;
        malformed.StaircaseRows[0].clear();
        stats.PacketRows = UINT32_C(0x76543210);
        if (wirehair_v2::SolvePrecodeSystemWithRuntime(
                malformed, config, runtime, packets, block_bytes,
                solved, &stats) != Wirehair_InvalidInput ||
            solved != solved_before || stats.PacketRows != UINT32_C(0x76543210))
        {
            std::fprintf(stderr,
                "solve: public runtime solve skipped system validation\n");
            return false;
        }
        if (runtime.Initialize(1u, P, config.MixCount) ||
            runtime.IsValidFor(K, P, config.MixCount) ||
            runtime.SourcePrime() != 0u || runtime.PrecodePrime() != 0u)
        {
            std::fprintf(stderr,
                "solve: invalid initialization retained packet runtime\n");
            return false;
        }
    }

    std::printf("packet evaluation fusion/cache/alias contract: PASS\n");
    return true;
}

bool CheckOddPacketPeelSeedInterleave()
{
    static const uint32_t K = 251u;
    static const uint32_t even_id = 100u;
    static const uint32_t odd_id = 101u;
    static const uint32_t seed_xor = 19u;
    wirehair_v2::PacketRowConfig base;
    base.PeelSeed = UINT32_C(0x12345678);
    base.MixCount = 2u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(
                K, UINT64_C(0x4f44445041434b54)),
            system))
    {
        return false;
    }
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, base.MixCount)) {
        return false;
    }

    wirehair_v2::SetOddPacketPeelSeedXorForTesting(0u);
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(false);
    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u)) {
        return false;
    }
    const std::vector<uint32_t> base_even =
        wirehair_v2::GeneratePacketMatrixRow(K, P, even_id, base);
    const std::vector<uint32_t> base_odd =
        wirehair_v2::GeneratePacketMatrixRow(K, P, odd_id, base);
    wirehair_v2::PacketRowConfig alternate = base;
    alternate.PeelSeed ^= seed_xor;
    const std::vector<uint32_t> alternate_odd =
        wirehair_v2::GeneratePacketMatrixRow(K, P, odd_id, alternate);

    wirehair_v2::SetOddPacketPeelSeedXorForTesting(seed_xor);
    const std::vector<uint32_t> blended_even =
        wirehair_v2::GeneratePacketMatrixRow(K, P, even_id, base);
    const std::vector<uint32_t> blended_odd =
        wirehair_v2::GeneratePacketMatrixRow(K, P, odd_id, base);
    const std::vector<uint32_t> cached_even =
        wirehair_v2::GeneratePacketMatrixRowWithRuntime(
            K, P, even_id, base, runtime);
    const std::vector<uint32_t> cached_odd =
        wirehair_v2::GeneratePacketMatrixRowWithRuntime(
            K, P, odd_id, base, runtime);
    static const uint32_t block_bytes = 17u;
    std::vector<uint8_t> intermediate(
        (size_t)(K + P) * block_bytes);
    for (size_t i = 0; i < intermediate.size(); ++i) {
        intermediate[i] = (uint8_t)(i * 131u + (i >> 5));
    }
    const std::vector<uint8_t> expected_odd = ReferencePacket(
        K, P, odd_id, alternate, intermediate.data(), block_bytes);
    std::vector<uint8_t> evaluated_odd(block_bytes, 0xa5u);
    uint64_t operations = 0u;
    const bool evaluated =
        wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
            system, base, runtime, intermediate.data(), block_bytes,
            odd_id, evaluated_odd.data(), &operations);
    wirehair_v2::SetOddPacketPeelSeedXorForTesting(0u);

    if (base_even.empty() || base_odd.empty() || alternate_odd.empty() ||
        alternate_odd == base_odd || blended_even != base_even ||
        blended_odd != alternate_odd || cached_even != blended_even ||
        cached_odd != blended_odd || !evaluated ||
        evaluated_odd != expected_odd || operations != alternate_odd.size())
    {
        std::fprintf(stderr,
            "solve: odd packet peel-seed interleave mismatch\n");
        return false;
    }
    std::printf("odd packet peel-seed interleave: PASS\n");
    return true;
}

bool CheckPacketRowSeedPermutation()
{
    static const uint32_t K = 251u;
    static const uint32_t block_id = 101u;
    static const uint32_t multiplier = UINT32_C(0x9e3779b1);
    static const uint32_t block_bytes = 17u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(
                K, UINT64_C(0x524f575045524d55)),
            system))
    {
        return false;
    }
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = UINT32_C(0x13579bdf);
    config.MixCount = 2u;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }

    wirehair_v2::SetOddPacketPeelSeedXorForTesting(0u);
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(false);
    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u)) {
        return false;
    }
    const uint32_t permuted_id = block_id * multiplier;
    const std::vector<uint32_t> expected_row =
        wirehair_v2::GeneratePacketMatrixRow(
            K, P, permuted_id, config);
    std::vector<uint8_t> intermediate(
        (size_t)(K + P) * block_bytes);
    for (size_t i = 0; i < intermediate.size(); ++i) {
        intermediate[i] = (uint8_t)(i * 67u + (i >> 3));
    }
    const std::vector<uint8_t> expected_block = ReferencePacket(
        K, P, permuted_id, config, intermediate.data(), block_bytes);

    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(multiplier)) {
        return false;
    }
    const std::vector<uint32_t> actual_row =
        wirehair_v2::GeneratePacketMatrixRow(K, P, block_id, config);
    const std::vector<uint32_t> cached_row =
        wirehair_v2::GeneratePacketMatrixRowWithRuntime(
            K, P, block_id, config, runtime);
    std::vector<uint8_t> actual_block(block_bytes, 0xa5u);
    uint64_t operations = 0u;
    const bool evaluated =
        wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
            system, config, runtime, intermediate.data(), block_bytes,
            block_id, actual_block.data(), &operations);
    const bool invalid_preserved =
        !wirehair_v2::SetPacketRowSeedMultiplierForTesting(0u) &&
        !wirehair_v2::SetPacketRowSeedMultiplierForTesting(2u) &&
        wirehair_v2::GeneratePacketMatrixRow(K, P, block_id, config) ==
            actual_row;
    (void)wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u);

    uint32_t avalanche_id = permuted_id;
    avalanche_id = (avalanche_id ^ (avalanche_id >> 16)) *
        UINT32_C(0x7feb352d);
    avalanche_id = (avalanche_id ^ (avalanche_id >> 15)) *
        UINT32_C(0x846ca68b);
    avalanche_id ^= avalanche_id >> 16;
    const std::vector<uint32_t> expected_avalanche_row =
        wirehair_v2::GeneratePacketMatrixRow(
            K, P, avalanche_id, config);
    const std::vector<uint8_t> expected_avalanche_block = ReferencePacket(
        K, P, avalanche_id, config, intermediate.data(), block_bytes);
    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(multiplier)) {
        return false;
    }
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(true);
    const std::vector<uint32_t> actual_avalanche_row =
        wirehair_v2::GeneratePacketMatrixRow(K, P, block_id, config);
    std::vector<uint8_t> actual_avalanche_block(block_bytes, 0xa5u);
    uint64_t avalanche_operations = 0u;
    const bool avalanche_evaluated =
        wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
            system, config, runtime, intermediate.data(), block_bytes,
            block_id, actual_avalanche_block.data(), &avalanche_operations);
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(false);
    (void)wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u);

    if (expected_row.empty() || actual_row != expected_row ||
        cached_row != expected_row || !evaluated ||
        actual_block != expected_block || operations != expected_row.size() ||
        !invalid_preserved || expected_avalanche_row.empty() ||
        actual_avalanche_row != expected_avalanche_row ||
        !avalanche_evaluated ||
        actual_avalanche_block != expected_avalanche_block ||
        avalanche_operations != expected_avalanche_row.size())
    {
        std::fprintf(stderr,
            "solve: packet row-seed permutation mismatch\n");
        return false;
    }
    std::printf("packet row-seed permutation: PASS\n");
    return true;
}

bool CheckPacketRuntimeBoundaries()
{
    struct Domain
    {
        uint32_t K;
        uint32_t P;
    };
    static const Domain kDomains[] = {
        { 2u, 2u },
        { 251u, 251u },
        { 252u, 250u },
        { 2u, 65521u },
        { 64000u, 1535u }
    };
    static const uint32_t kIds[] = {
        0u, 1u, UINT32_C(0xf1234567), UINT32_MAX
    };
    for (const Domain& domain : kDomains)
    {
        wirehair_v2::PacketRowConfig config;
        config.PeelSeed = UINT32_C(0x8d12a4f7);
        config.MixCount = std::min<uint32_t>(
            wirehair_v2::kCertifiedPacketMixCount, domain.P);
        wirehair_v2::PacketRowRuntime runtime;
        if (!runtime.Initialize(domain.K, domain.P, config.MixCount) ||
            runtime.SourcePrime() !=
                wirehair::NextPrime16((uint16_t)domain.K) ||
            runtime.PrecodePrime() !=
                wirehair::NextPrime16((uint16_t)domain.P))
        {
            std::fprintf(stderr,
                "solve: runtime boundary init failed K=%u P=%u\n",
                domain.K, domain.P);
            return false;
        }
        for (uint32_t id : kIds)
        {
            if (wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                    domain.K, domain.P, id, config, runtime) !=
                wirehair_v2::GeneratePacketMatrixRow(
                    domain.K, domain.P, id, config))
            {
                std::fprintf(stderr,
                    "solve: runtime boundary row mismatch K=%u P=%u id=%u\n",
                    domain.K, domain.P, id);
                return false;
            }
        }
    }
    std::printf("packet runtime domain boundaries: PASS\n");
    return true;
}

bool CheckTinyDenseOracle()
{
    const uint32_t K = 2u;
    const uint32_t block_bytes = 3u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x54494e594f524143));
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x13579bdf);
    base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    uint32_t attempt = 0u;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config, &attempt) !=
            Wirehair_Success)
    {
        std::fprintf(stderr, "solve: tiny oracle configuration failed\n");
        return false;
    }

    const uint32_t L = K + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint32_t coefficient_period = 256u - system.Params.HeavyRows;
    const uint64_t empty_residue_muladds =
        (uint64_t)(coefficient_period - std::min(coefficient_period, L)) *
        system.Params.HeavyRows;
    static const uint64_t kExpectedExecutedMulAdds = 528u;
    static const uint32_t kBlockBytes[] = { 1u, block_bytes, 1280u };
    uint64_t expected_muladds = UINT64_MAX;
    for (uint32_t bytes : kBlockBytes)
    {
        std::vector<uint8_t> message((size_t)K * bytes);
        for (size_t i = 0; i < message.size(); ++i) {
            message[i] = (uint8_t)(i * 29u + bytes * 7u + 3u);
        }
        std::vector<wirehair_v2::SolvePacket> packets(K);
        for (uint32_t id = 0; id < K; ++id) {
            packets[id].BlockId = id;
            packets[id].Data = message.data() + (size_t)id * bytes;
        }
        std::vector<uint8_t> solved;
        wirehair_v2::PrecodeSolveStats stats;
        if (wirehair_v2::SolvePrecodeSystem(
                system, config, packets, bytes, solved, &stats) !=
            Wirehair_Success)
        {
            std::fprintf(stderr,
                "solve: tiny sparse solve failed bb=%u\n", bytes);
            return false;
        }
        std::vector<uint8_t> oracle;
        const WirehairResult oracle_result =
            wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
                system, config, packets, bytes, oracle);
        if (oracle_result != Wirehair_Success || oracle != solved)
        {
            std::fprintf(stderr,
                "solve: extracted tiny dense oracle mismatch attempt=%u "
                "bb=%u result=%d\n",
                attempt, bytes, (int)oracle_result);
            return false;
        }
        if (stats.BlockMulAdds != kExpectedExecutedMulAdds ||
            stats.BlockMulAdds >= empty_residue_muladds ||
            (expected_muladds != UINT64_MAX &&
             stats.BlockMulAdds != expected_muladds))
        {
            std::fprintf(stderr,
                "solve: tiny heavy residue work mismatch bb=%u muladds=%llu "
                "empty_residue_muladds=%llu expected=%llu\n",
                bytes,
                (unsigned long long)stats.BlockMulAdds,
                (unsigned long long)empty_residue_muladds,
                (unsigned long long)expected_muladds);
            return false;
        }
        expected_muladds = stats.BlockMulAdds;
    }
    std::printf(
        "tiny heavy residues: L=%u period=%u skipped=%llu muladds=%llu: PASS\n",
        L,
        coefficient_period,
        (unsigned long long)empty_residue_muladds,
        (unsigned long long)expected_muladds);

    wirehair_v2::PrecodeParams hashed_params = params;
    hashed_params.HeavyFamily =
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero;
    wirehair_v2::PrecodeSystem hashed_system;
    wirehair_v2::PacketRowConfig hashed_config;
    if (wirehair_v2::SelectSystematicConfiguration(
            hashed_params, base_config, hashed_system, hashed_config) !=
        Wirehair_Success)
    {
        std::fprintf(stderr, "solve: tiny hashed configuration failed\n");
        return false;
    }
    std::vector<uint8_t> hashed_message(K * block_bytes);
    for (size_t i = 0; i < hashed_message.size(); ++i) {
        hashed_message[i] = (uint8_t)(i * 73u + 19u);
    }
    std::vector<wirehair_v2::SolvePacket> hashed_packets(K);
    for (uint32_t id = 0u; id < K; ++id) {
        hashed_packets[id].BlockId = id;
        hashed_packets[id].Data =
            hashed_message.data() + (size_t)id * block_bytes;
    }
    std::vector<uint8_t> hashed_solved;
    std::vector<uint8_t> hashed_oracle;
    if (wirehair_v2::SolvePrecodeSystem(
            hashed_system, hashed_config, hashed_packets, block_bytes,
            hashed_solved) != Wirehair_Success ||
        wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
            hashed_system, hashed_config, hashed_packets, block_bytes,
            hashed_oracle) != Wirehair_Success ||
        hashed_solved != hashed_oracle)
    {
        std::fprintf(stderr, "solve: tiny hashed dense oracle mismatch\n");
        return false;
    }
    std::printf("tiny hashed-family dense oracle: PASS\n");
    return true;
}

bool CheckHeavyCoefficientBoundaryOracle()
{
    static const uint32_t kHeavyRows[] = {0u, 1u, 12u, 128u};
    static const wirehair_v2::HeavyCoefficientFamily kFamilies[] = {
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy,
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero
    };
    const uint32_t K = 8u;
    const uint32_t block_bytes = 5u;
    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 97u + (i >> 2) + 41u);
    }
    std::vector<wirehair_v2::SolvePacket> packets(K);
    for (uint32_t id = 0u; id < K; ++id) {
        packets[id].BlockId = id;
        packets[id].Data = message.data() + (size_t)id * block_bytes;
    }

    for (wirehair_v2::HeavyCoefficientFamily family : kFamilies)
    {
        for (uint32_t heavy_rows : kHeavyRows)
        {
            wirehair_v2::PrecodeParams params =
                wirehair_v2::MakeCertifiedParams(
                    K,
                    UINT64_C(0x4845415659424f55) ^
                        ((uint64_t)heavy_rows << 8) ^ (uint32_t)family);
            params.HeavyRows = heavy_rows;
            params.HeavyFamily = family;
            wirehair_v2::PacketRowConfig base_config;
            base_config.PeelSeed =
                UINT32_C(0xa17e31d9) ^ heavy_rows ^ (uint32_t)family;
            base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
            wirehair_v2::PrecodeSystem system;
            wirehair_v2::PacketRowConfig config;
            if (wirehair_v2::SelectSystematicConfiguration(
                    params, base_config, system, config) != Wirehair_Success)
            {
                std::fprintf(stderr,
                    "solve: heavy boundary configuration failed H=%u "
                    "family=%u\n",
                    heavy_rows, (unsigned)family);
                return false;
            }

            std::vector<uint8_t> solved;
            std::vector<uint8_t> oracle;
            if (wirehair_v2::SolvePrecodeSystem(
                    system, config, packets, block_bytes, solved) !=
                    Wirehair_Success ||
                wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
                    system, config, packets, block_bytes, oracle) !=
                    Wirehair_Success ||
                solved != oracle ||
                !wirehair_v2::VerifyPrecodeSolution(
                    system, config, packets, solved.data(), block_bytes))
            {
                std::fprintf(stderr,
                    "solve: heavy boundary oracle mismatch H=%u family=%u\n",
                    heavy_rows, (unsigned)family);
                return false;
            }
        }
    }
    std::printf("heavy H=0/1/12/128 periodic/hashed dense oracle: PASS\n");
    return true;
}

bool CheckBinaryQuotientBoundary()
{
    const uint32_t K = 2u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(
            K, UINT64_C(0x51554f5449454e54));
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x2468ace0);
    base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config) != Wirehair_Success)
    {
        std::fprintf(stderr, "solve: quotient configuration failed\n");
        return false;
    }

    uint64_t muladds[2] = {};
    const uint32_t block_sizes[2] = {
        wirehair_v2::kBinaryQuotientMinBlockBytes - 1u,
        wirehair_v2::kBinaryQuotientMinBlockBytes
    };
    for (uint32_t case_i = 0; case_i < 2u; ++case_i)
    {
        const uint32_t block_bytes = block_sizes[case_i];
        std::vector<uint8_t> message((size_t)K * block_bytes);
        for (size_t i = 0; i < message.size(); ++i) {
            message[i] = (uint8_t)(i * 149u + (i >> 7) + 23u);
        }
        std::vector<wirehair_v2::SolvePacket> packets(K);
        for (uint32_t id = 0; id < K; ++id) {
            packets[id].BlockId = id;
            packets[id].Data =
                message.data() + (size_t)id * block_bytes;
        }
        std::vector<uint8_t> solved;
        wirehair_v2::PrecodeSolveStats stats;
        if (wirehair_v2::SolvePrecodeSystem(
                system, config, packets, block_bytes, solved, &stats) !=
            Wirehair_Success)
        {
            std::fprintf(stderr,
                "solve: quotient boundary failed bb=%u\n", block_bytes);
            return false;
        }
        std::vector<uint8_t> oracle;
        if (wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
                system, config, packets, block_bytes, oracle) !=
                Wirehair_Success ||
            oracle != solved)
        {
            std::fprintf(stderr,
                "solve: quotient boundary oracle mismatch bb=%u\n",
                block_bytes);
            return false;
        }
        muladds[case_i] = stats.BlockMulAdds;
    }
    std::printf(
        "binary quotient threshold: legacy_muladds=%llu "
        "quotient_muladds=%llu: PASS\n",
        (unsigned long long)muladds[0],
        (unsigned long long)muladds[1]);

    // Exercise the quotient on a nontrivial full-rank systematic system and
    // compare both successful and inconsistent RHS behavior to an independent
    // bounded dense GF(256) oracle.
    const uint32_t oracle_K = 64u;
    const uint32_t oracle_block_bytes =
        wirehair_v2::kBinaryQuotientMinBlockBytes;
    wirehair_v2::PrecodeParams oracle_params =
        wirehair_v2::MakeCertifiedParams(
            oracle_K, UINT64_C(0x4b363451554f5449));
    wirehair_v2::PrecodeSystem oracle_system;
    wirehair_v2::PacketRowConfig oracle_config;
    if (wirehair_v2::SelectSystematicConfiguration(
            oracle_params, base_config, oracle_system, oracle_config) !=
        Wirehair_Success)
    {
        std::fprintf(stderr, "solve: K64 quotient configuration failed\n");
        return false;
    }
    std::vector<uint8_t> oracle_message(
        (size_t)oracle_K * oracle_block_bytes);
    for (size_t i = 0; i < oracle_message.size(); ++i) {
        oracle_message[i] = (uint8_t)(i * 109u + (i >> 9) + 31u);
    }
    std::vector<wirehair_v2::SolvePacket> oracle_packets(oracle_K);
    for (uint32_t id = 0u; id < oracle_K; ++id) {
        oracle_packets[id].BlockId = id;
        oracle_packets[id].Data =
            oracle_message.data() + (size_t)id * oracle_block_bytes;
    }
    std::vector<uint8_t> production;
    std::vector<uint8_t> dense;
    if (wirehair_v2::SolvePrecodeSystem(
            oracle_system, oracle_config, oracle_packets,
            oracle_block_bytes, production) != Wirehair_Success ||
        wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
            oracle_system, oracle_config, oracle_packets,
            oracle_block_bytes, dense) != Wirehair_Success ||
        production != dense)
    {
        std::fprintf(stderr, "solve: K64 quotient dense oracle mismatch\n");
        return false;
    }

    std::vector<uint8_t> conflicting(
        oracle_message.begin(),
        oracle_message.begin() + oracle_block_bytes);
    conflicting[0] ^= 1u;
    wirehair_v2::SolvePacket conflict_packet;
    conflict_packet.BlockId = 0u;
    conflict_packet.Data = conflicting.data();
    oracle_packets.push_back(conflict_packet);
    production.assign(11u, 0xa5u);
    dense.assign(13u, 0x6du);
    const std::vector<uint8_t> production_before = production;
    const std::vector<uint8_t> dense_before = dense;
    if (wirehair_v2::SolvePrecodeSystem(
            oracle_system, oracle_config, oracle_packets,
            oracle_block_bytes, production) != Wirehair_Error ||
        wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
            oracle_system, oracle_config, oracle_packets,
            oracle_block_bytes, dense) != Wirehair_Error ||
        production != production_before || dense != dense_before)
    {
        std::fprintf(stderr,
            "solve: K64 quotient conflicting RHS acceptance/no-write mismatch\n");
        return false;
    }
    std::printf("K64 binary quotient dense oracle/conflict: PASS\n");
    return true;
}

bool CheckIncrementalResumeCase(uint32_t block_bytes)
{
    const uint32_t K = 64u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(
            K, UINT64_C(0x524553554d455354));
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x91e10da5);
    base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config) != Wirehair_Success)
    {
        std::fprintf(stderr, "solve: resume configuration failed\n");
        return false;
    }

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 197u + (i >> 3) + 41u);
    }
    std::vector<wirehair_v2::SolvePacket> systematic(K);
    for (uint32_t id = 0; id < K; ++id) {
        systematic[id].BlockId = id;
        systematic[id].Data =
            message.data() + (size_t)id * block_bytes;
    }
    std::vector<uint8_t> expected;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, systematic, block_bytes, expected) !=
            Wirehair_Success)
    {
        return false;
    }

    const uint32_t L = K + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::SolveValueStorage lazy_expected;
    lazy_expected.assign(11u, 0x3cu);
    wirehair_v2::PrecodeSolveStats lazy_full_stats;
    WirehairResult lazy_full_result = Wirehair_Error;
    {
        SolveValueArenaPoisonScope poison;
        lazy_full_result = wirehair_v2::SolvePrecodeSystem(
            system, config, systematic, block_bytes,
            lazy_expected, &lazy_full_stats);
    }
    if (lazy_full_result != Wirehair_Success ||
        lazy_expected.size() != expected.size() ||
        !std::equal(
            lazy_expected.begin(), lazy_expected.end(), expected.begin()) ||
        lazy_full_stats.SolveValueArenaBytes !=
            (uint64_t)L * block_bytes ||
        lazy_full_stats.SolveValueArenaEagerZeroBytes != 0u ||
        lazy_full_stats.SolveValueArenaCommitCopyBytes != 0u)
    {
        std::fprintf(stderr,
            "solve: poisoned no-init full solve failed bb=%u\n", block_bytes);
        return false;
    }
    std::vector<wirehair_v2::SolvePacket> null_packet = systematic;
    null_packet[K / 2u].Data = nullptr;
    wirehair_v2::SolveValueStorage null_output;
    null_output.assign(5u, 0x72u);
    const std::vector<uint8_t> null_sentinel(
        null_output.begin(), null_output.end());
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, null_packet, block_bytes, null_output) !=
            Wirehair_InvalidInput ||
        null_output.size() != null_sentinel.size() ||
        !std::equal(
            null_output.begin(), null_output.end(), null_sentinel.begin()))
    {
        std::fprintf(stderr,
            "solve: no-init null packet contract failed bb=%u\n", block_bytes);
        return false;
    }
    wirehair_v2::SolveValueStorage oom_output;
    oom_output.assign(5u, 0x39u);
    const std::vector<uint8_t> oom_sentinel(
        oom_output.begin(), oom_output.end());
    wirehair_v2::PrecodeSolveResumeState oom_resume;
    oom_resume.Active = true;
    oom_resume.SourceCount = UINT32_C(0x13579bdf);
    wirehair_v2::PrecodeSolveStats oom_stats;
    oom_stats.PacketRows = UINT32_C(0x2468ace0);
    WirehairResult oom_result = Wirehair_Error;
    {
        SolveValueArenaAllocationFailureScope fail_arena;
        oom_result = wirehair_v2::SolvePrecodeSystem(
            system, config, systematic, block_bytes,
            oom_output, &oom_stats, &oom_resume);
    }
    if (oom_result != Wirehair_OOM ||
        oom_output.size() != oom_sentinel.size() ||
        !std::equal(
            oom_output.begin(), oom_output.end(), oom_sentinel.begin()) ||
        !oom_resume.Active ||
        oom_resume.SourceCount != UINT32_C(0x13579bdf) ||
        oom_stats.PacketRows != UINT32_C(0x2468ace0))
    {
        std::fprintf(stderr,
            "solve: no-init arena OOM was not transactional bb=%u\n",
            block_bytes);
        return false;
    }

    std::vector<wirehair_v2::SolvePacket> deficient(K);
    for (wirehair_v2::SolvePacket& packet : deficient) {
        packet.BlockId = 0u;
        packet.Data = message.data();
    }

    wirehair_v2::SolveValueStorage lazy_output;
    lazy_output.assign(11u, 0x6du);
    const std::vector<uint8_t> lazy_sentinel(
        lazy_output.begin(), lazy_output.end());
    wirehair_v2::PrecodeSolveResumeState lazy_resume;
    wirehair_v2::PrecodeSolveStats lazy_resume_stats;
    WirehairResult lazy_deficient_result = Wirehair_Error;
    {
        SolveValueArenaPoisonScope poison;
        lazy_deficient_result = wirehair_v2::SolvePrecodeSystem(
            system, config, deficient, block_bytes,
            lazy_output, &lazy_resume_stats, &lazy_resume);
    }
    if (lazy_deficient_result != Wirehair_NeedMore ||
        !lazy_resume.Active ||
        lazy_output.size() != lazy_sentinel.size() ||
        !std::equal(
            lazy_output.begin(), lazy_output.end(), lazy_sentinel.begin()) ||
        lazy_resume_stats.SolveValueArenaBytes !=
            (uint64_t)L * block_bytes ||
        lazy_resume_stats.SolveValueArenaEagerZeroBytes != 0u ||
        lazy_resume_stats.SolveValueArenaCommitCopyBytes != 0u)
    {
        std::fprintf(stderr,
            "solve: poisoned no-init checkpoint failed bb=%u\n", block_bytes);
        return false;
    }
    const wirehair_v2::PrecodeSolveResumeState lazy_resume_copy = lazy_resume;
    if (lazy_resume_copy.Values.size() != lazy_resume.Values.size() ||
        !std::equal(
            lazy_resume_copy.Values.begin(), lazy_resume_copy.Values.end(),
            lazy_resume.Values.begin()))
    {
        std::fprintf(stderr,
            "solve: no-init checkpoint copy failed bb=%u\n", block_bytes);
        return false;
    }
    // Checkpoint publication materializes these implicit-zero constants so the
    // state stays copyable.  Re-poison their values numerically to prove resume
    // does not depend on them before the residual variables are solved.
    for (uint32_t column : lazy_resume.InactiveColumns) {
        std::memset(
            lazy_resume.Values.data() + (size_t)column * block_bytes,
            0xa5,
            block_bytes);
    }
    wirehair_v2::SolveValueStorage lazy_ignored;
    lazy_ignored.assign(7u, 0x4bu);
    if (wirehair_v2::ResumePrecodeSystem(
            system, config, 0u, message.data(), block_bytes,
            lazy_resume, lazy_ignored, nullptr, false) != Wirehair_NeedMore)
    {
        std::fprintf(stderr,
            "solve: poisoned checkpoint duplicate failed bb=%u\n", block_bytes);
        return false;
    }
    WirehairResult lazy_resume_result = Wirehair_NeedMore;
    bool checked_lazy_resume_oom = false;
    for (uint32_t id = 1u; id < K; ++id)
    {
        if (!checked_lazy_resume_oom &&
            lazy_resume.Rank + 1u == lazy_resume.InactiveCount)
        {
            const wirehair_v2::PrecodeSolveResumeState state_before =
                lazy_resume;
            const size_t persistent_bytes_before =
                lazy_resume.PersistentBytes();
            const std::vector<uint8_t> output_before(
                lazy_output.begin(), lazy_output.end());
            WirehairResult oom_resume_result = Wirehair_Error;
            {
                SolveValueArenaAllocationFailureScope fail_arena;
                oom_resume_result = wirehair_v2::ResumePrecodeSystem(
                    system, config, id,
                    message.data() + (size_t)id * block_bytes,
                    block_bytes, lazy_resume, lazy_output, nullptr, true);
            }
            if (oom_resume_result != Wirehair_OOM ||
                !SameResumeState(lazy_resume, state_before) ||
                lazy_resume.PersistentBytes() != persistent_bytes_before ||
                lazy_output.size() != output_before.size() ||
                !std::equal(
                    lazy_output.begin(), lazy_output.end(),
                    output_before.begin()))
            {
                std::fprintf(stderr,
                    "solve: no-init resume OOM was not transactional bb=%u\n",
                    block_bytes);
                return false;
            }
            checked_lazy_resume_oom = true;
        }
        lazy_resume_result = wirehair_v2::ResumePrecodeSystem(
            system, config, id,
            message.data() + (size_t)id * block_bytes,
            block_bytes, lazy_resume, lazy_output,
            &lazy_resume_stats, true);
        if (id + 1u < K && lazy_resume_result != Wirehair_NeedMore) {
            return false;
        }
    }
    if (lazy_resume_result != Wirehair_Success || lazy_resume.Active ||
        !checked_lazy_resume_oom ||
        lazy_output.size() != expected.size() ||
        !std::equal(
            lazy_output.begin(), lazy_output.end(), expected.begin()) ||
        lazy_resume_stats.SolveValueArenaCommitCopyBytes !=
            (uint64_t)L * block_bytes)
    {
        std::fprintf(stderr,
            "solve: poisoned no-init resume failed bb=%u\n", block_bytes);
        return false;
    }

    std::vector<uint8_t> output(11u, 0xa5u);
    const std::vector<uint8_t> sentinel = output;
    wirehair_v2::PrecodeSolveResumeState resume;
    wirehair_v2::PrecodeSolveStats stats;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, deficient, block_bytes,
            output, &stats, &resume) != Wirehair_NeedMore ||
        !resume.Active || output != sentinel ||
        resume.Rank >= resume.InactiveCount)
    {
        std::fprintf(stderr, "solve: rank-deficient checkpoint missing\n");
        return false;
    }

    const uint32_t rank_before = resume.Rank;
    const size_t bytes_before = resume.PersistentBytes();
    const std::vector<uint8_t> coefficient_scratch_before =
        resume.CoefficientScratch;
    const std::vector<uint8_t> rhs_scratch_before = resume.RhsScratch;
    const std::vector<uint8_t> pivot_coefficients_before =
        resume.PivotCoefficients;
    const std::vector<uint8_t> pivot_rhs_before = resume.PivotRhs;
    const std::vector<uint8_t> have_pivot_before = resume.HavePivot;
    if (wirehair_v2::ResumePrecodeSystem(
            system, config, 0u, message.data(), block_bytes,
            resume, output, nullptr, false) != Wirehair_NeedMore ||
        resume.Rank != rank_before ||
        resume.PersistentBytes() != bytes_before ||
        resume.CoefficientScratch != coefficient_scratch_before ||
        resume.RhsScratch != rhs_scratch_before ||
        resume.PivotCoefficients != pivot_coefficients_before ||
        resume.PivotRhs != pivot_rhs_before ||
        resume.HavePivot != have_pivot_before ||
        output != sentinel)
    {
        std::fprintf(stderr, "solve: exact checkpoint duplicate changed state\n");
        return false;
    }
    std::vector<uint8_t> corrupt(message.begin(), message.begin() + block_bytes);
    corrupt[0] ^= 1u;
    if (wirehair_v2::ResumePrecodeSystem(
            system, config, 0u, corrupt.data(), block_bytes,
            resume, output, nullptr, false) != Wirehair_Error ||
        resume.Rank != rank_before ||
        resume.CoefficientScratch != coefficient_scratch_before ||
        resume.RhsScratch != rhs_scratch_before ||
        resume.PivotCoefficients != pivot_coefficients_before ||
        resume.PivotRhs != pivot_rhs_before ||
        resume.HavePivot != have_pivot_before ||
        output != sentinel)
    {
        std::fprintf(stderr,
            "solve: conflicting checkpoint duplicate was accepted\n");
        return false;
    }

    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t id = 1u; id < K; ++id)
    {
        result = wirehair_v2::ResumePrecodeSystem(
            system,
            config,
            id,
            message.data() + (size_t)id * block_bytes,
            block_bytes,
            resume,
            output,
            &stats,
            true);
        if (id + 1u < K && result != Wirehair_NeedMore) {
            std::fprintf(stderr,
                "solve: checkpoint completed early id=%u result=%d\n",
                id, (int)result);
            return false;
        }
    }
    if (result != Wirehair_Success || resume.Active || output != expected ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, systematic,
            output.data(), block_bytes))
    {
        std::fprintf(stderr, "solve: resumed solution mismatch\n");
        return false;
    }
    std::printf(
        "incremental rank-deficient resume bb=%u: PASS\n", block_bytes);
    return true;
}

bool CheckIncrementalResume()
{
    return CheckIncrementalResumeCase(17u) &&
        CheckIncrementalResumeCase(
            wirehair_v2::kBinaryQuotientMinBlockBytes);
}

bool SamePrecodeSolveWork(
    const wirehair_v2::PrecodeSolveStats& a,
    const wirehair_v2::PrecodeSolveStats& b)
{
    return a.PacketRows == b.PacketRows &&
        a.PeeledColumns == b.PeeledColumns &&
        a.InactivatedColumns == b.InactivatedColumns &&
        a.ResidualRows == b.ResidualRows &&
        a.ResidualRank == b.ResidualRank &&
        a.BinaryResidualRank == b.BinaryResidualRank &&
        a.BinaryRowReferences == b.BinaryRowReferences &&
        a.BinaryRowStorageBytes == b.BinaryRowStorageBytes &&
        a.BinaryAdjacencyStorageBytes == b.BinaryAdjacencyStorageBytes &&
        a.BinaryRowStorageAllocations == b.BinaryRowStorageAllocations &&
        a.BinaryAdjacencyStorageAllocations ==
            b.BinaryAdjacencyStorageAllocations &&
        a.BlockXors == b.BlockXors &&
        a.BlockMulAdds == b.BlockMulAdds &&
        a.PacketSeedAttempt == b.PacketSeedAttempt
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        && a.MixedJointSourceXors == b.MixedJointSourceXors
        && a.MixedJointMarginalXors == b.MixedJointMarginalXors
        && a.MixedJointMarginalCopies == b.MixedJointMarginalCopies
        && a.MixedJointScratchBytes == b.MixedJointScratchBytes
        && a.MixedJointActiveDeltas == b.MixedJointActiveDeltas
        && a.MixedDualSourceColumns == b.MixedDualSourceColumns
#endif
        ;
}

bool CheckMixedProjectionResidueBucketsOracleForPeriod(
    uint32_t period,
    wirehair_v2::MixedCoefficientGeometry geometry,
    uint32_t extension_rows,
    uint32_t residue_skew = 0u,
    wirehair_v2::MixedResidueSchedule residue_schedule =
        wirehair_v2::MixedResidueSchedule::Constant,
    bool independent_extension_residues = false,
    uint32_t subfield_rows = wirehair_v2::kMixedGF256Rows,
    wirehair_v2::MixedResidueBucketMode bucket_mode =
        wirehair_v2::MixedResidueBucketMode::Automatic,
    bool dense_two_anchor = false,
    uint32_t grouped_gf256_rows = 0u)
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
    // The grouped suffix validates the complete H12/shared-X/constant-A
    // configuration and prerequisite setters clear it, so configure it last.
    MixedGroupedGF256RowsScope grouped_scope(grouped_gf256_rows);
    if (!geometry_scope.IsValid() || !subfield_scope.IsValid() ||
        !rows_scope.IsValid() || !period_scope.IsValid() ||
        !skew_scope.IsValid() ||
        !schedule_scope.IsValid() || !independent_scope.IsValid() ||
        !bucket_mode_scope.IsValid() || !grouped_scope.IsValid())
    {
        std::fprintf(stderr,
            "solve: invalid mixed projection scope g=%u s=%u r=%u p=%u "
            "k=%u q=%u i=%u b=%u c=%u\n",
            geometry_scope.IsValid() ? 1u : 0u,
            subfield_scope.IsValid() ? 1u : 0u,
            rows_scope.IsValid() ? 1u : 0u,
            period_scope.IsValid() ? 1u : 0u,
            skew_scope.IsValid() ? 1u : 0u,
            schedule_scope.IsValid() ? 1u : 0u,
            independent_scope.IsValid() ? 1u : 0u,
            bucket_mode_scope.IsValid() ? 1u : 0u,
            grouped_scope.IsValid() ? 1u : 0u);
        return false;
    }
    wirehair_v2::ResetMixedProjectionOracleComparisonsForTesting();
    MixedProjectionOracleScope oracle_scope;
    static const uint32_t kBlockCounts[] = {
        2u, 3u, 10u, 63u, 64u, 127u, 128u,
        // With H12, S=30 and D2+H=24, so these exercise exact total-column
        // counts L=243, 244, and 245 around the full-period transition.
        189u, 190u, 191u,
        243u, 244u, 245u, 320u, 1000u,
        // Automatic P32 joint-delta crossover (decoder dispatch coverage).
        3200u,
        // The independent case also covers the large-payload dual-bucket
        // mixed-RHS path selected by SolveMixedCompletionQuotient.
        30000u
    };
    // The first pair straddles the block-zero identity boundary L-H <= P;
    // the second pair straddles a later partial-period boundary.  The tiny
    // table is intentionally explicit because its certified S varies by K.
    static const uint32_t kGroupedBlockCountsP32[] = {
        10u, 11u, 89u, 91u
    };
    static const uint32_t kGroupedBlockCountsP48[] = {
        20u, 21u, 60u, 65u
    };
    static const uint32_t kGroupedBlockCountsP64[] = {
        35u, 31u, 89u, 91u
    };
    const uint32_t* grouped_block_counts = kGroupedBlockCountsP64;
    if (period == 32u) {
        grouped_block_counts = kGroupedBlockCountsP32;
    }
    else if (period == 48u) {
        grouped_block_counts = kGroupedBlockCountsP48;
    }
    const uint32_t* const block_counts = grouped_gf256_rows != 0u ?
        grouped_block_counts : kBlockCounts;
    const size_t configured_case_count = grouped_gf256_rows != 0u ?
        sizeof(kGroupedBlockCountsP32) /
            sizeof(kGroupedBlockCountsP32[0]) :
        sizeof(kBlockCounts) / sizeof(kBlockCounts[0]);
    size_t expected_case_count = 0u;
    for (size_t case_index = 0;
         case_index < configured_case_count;
         ++case_index)
    {
        const uint32_t K = block_counts[case_index];
        if (K == 30000u && !independent_extension_residues) continue;
        ++expected_case_count;
        const bool automatic_joint_case = K == 3200u &&
            independent_extension_residues && period == 32u &&
            bucket_mode == wirehair_v2::MixedResidueBucketMode::Automatic;
        const uint32_t block_bytes = automatic_joint_case ? 4096u :
            (K == 30000u ? 1024u :
             ((case_index & 1u) == 0u ? 2u : 6u));
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeMixedParams(
                K,
                UINT64_C(0x70726f6a65637400) ^
                    ((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)));
        params.DenseTwoAnchor = dense_two_anchor;
        wirehair_v2::PacketRowConfig base_config;
        base_config.PeelSeed =
            UINT32_C(0x6f72636c) ^ K * UINT32_C(0x9e3779b9);
        base_config.MixCount = (case_index & 1u) == 0u ? 2u : 3u;
        wirehair_v2::PrecodeSystem system;
        wirehair_v2::PacketRowConfig config;
        if (wirehair_v2::SelectSystematicConfiguration(
                params, base_config, system, config) != Wirehair_Success)
        {
            std::fprintf(stderr,
                "solve: mixed projection oracle selection failed K=%u\n", K);
            return false;
        }
        const uint32_t column_count = K + system.Params.Staircase +
            system.Params.DenseRows + system.Params.HeavyRows;
        const uint32_t first_heavy_column =
            column_count - system.Params.HeavyRows;
        if (grouped_gf256_rows != 0u)
        {
            const bool identity_boundary_case = case_index < 2u;
            const bool boundary_side_ok = identity_boundary_case ?
                (case_index == 0u ?
                    first_heavy_column <= period :
                    first_heavy_column > period) :
                (period == 48u ?
                    (case_index == 2u ?
                        first_heavy_column < 2u * period :
                        first_heavy_column > 2u * period) :
                    first_heavy_column % period ==
                        (K == 89u ? period - 1u : 1u));
            if (!boundary_side_ok)
            {
                std::fprintf(stderr,
                    "solve: grouped mixed projection boundary "
                    "case=%zu K=%u L-H=%u P=%u\n",
                    case_index,
                    K, first_heavy_column, period);
                return false;
            }
        }
        const bool grouped_block_zero_identity =
            grouped_gf256_rows != 0u && first_heavy_column <= period;
        const bool secondary_schedule =
            independent_extension_residues ||
            (grouped_gf256_rows != 0u && !grouped_block_zero_identity);
        const bool automatic_joint_bucketed =
            bucket_mode == wirehair_v2::MixedResidueBucketMode::Automatic &&
            wirehair_v2::UseAutomaticMixedJointResidueBucketsForTesting(
                K, block_bytes, period);
        const bool joint_bucketed = secondary_schedule &&
            (bucket_mode ==
                 wirehair_v2::MixedResidueBucketMode::JointDelta ||
             automatic_joint_bucketed);
        const bool dual_bucketed = secondary_schedule &&
            (bucket_mode == wirehair_v2::MixedResidueBucketMode::Dual ||
             (bucket_mode ==
                  wirehair_v2::MixedResidueBucketMode::Automatic &&
              !automatic_joint_bucketed && column_count >= 30000u &&
              block_bytes >= 1024u &&
              (uint64_t)2u * period * block_bytes <=
                  (UINT64_C(128) << 10)));
        if (K >= 189u && K <= 191u &&
            column_count != 243u + (K - 189u) +
                (extension_rows - wirehair_v2::kMixedGF16Rows) +
                (subfield_rows - wirehair_v2::kMixedGF256Rows))
        {
            std::fprintf(stderr,
                "solve: mixed projection boundary K=%u produced L=%u\n",
                K, column_count);
            return false;
        }

        std::vector<uint8_t> message((size_t)K * block_bytes);
        for (size_t i = 0; i < message.size(); ++i) {
            message[i] = (uint8_t)(
                i * 157u + (i >> 3) + K * 29u + case_index);
        }
        std::vector<wirehair_v2::SolvePacket> systematic(K);
        for (uint32_t id = 0u; id < K; ++id) {
            systematic[id].BlockId = id;
            systematic[id].Data =
                message.data() + (size_t)id * block_bytes;
        }
        std::vector<uint8_t> expected;
        wirehair_v2::PrecodeSolveStats systematic_stats;
        const uint64_t systematic_comparisons_before =
            wirehair_v2::MixedProjectionOracleComparisonsForTesting();
        const WirehairResult systematic_result =
            wirehair_v2::SolvePrecodeSystem(
                system, config, systematic, block_bytes, expected,
                &systematic_stats);
        const uint64_t systematic_comparisons_after =
            wirehair_v2::MixedProjectionOracleComparisonsForTesting();
        const bool systematic_verified =
            systematic_result == Wirehair_Success &&
            wirehair_v2::VerifyPrecodeSolution(
                system, config, systematic, expected.data(), block_bytes);
        if (systematic_result != Wirehair_Success ||
            systematic_comparisons_after !=
                systematic_comparisons_before + 1u ||
            !systematic_verified)
        {
            std::fprintf(stderr,
                "solve: mixed projection systematic oracle failed K=%u "
                "result=%d comparisons=%llu->%llu verified=%u grouped=%u "
                "independent=%u\n",
                K, (int)systematic_result,
                (unsigned long long)systematic_comparisons_before,
                (unsigned long long)systematic_comparisons_after,
                systematic_verified ? 1u : 0u,
                grouped_gf256_rows,
                independent_extension_residues ? 1u : 0u);
            return false;
        }
        const auto bucket_stats_valid = [&](
            const wirehair_v2::PrecodeSolveStats& stats) -> bool
        {
            if (joint_bucketed)
            {
                if (grouped_gf256_rows != 0u)
                {
                    return stats.MixedJointSourceXors > 0u &&
                        stats.MixedJointActiveDeltas > 0u &&
                        stats.MixedJointMarginalXors > 0u &&
                        stats.MixedJointMarginalCopies > 0u &&
                        stats.MixedJointScratchBytes > 0u &&
                        stats.MixedDualSourceColumns == 0u;
                }
                return stats.MixedJointSourceXors +
                            stats.InactivatedColumns == column_count &&
                    stats.MixedJointActiveDeltas > 0u &&
                    stats.MixedJointMarginalXors ==
                        (uint64_t)2u * period *
                            (stats.MixedJointActiveDeltas - 1u) &&
                    stats.MixedJointMarginalCopies ==
                        (uint64_t)2u * period &&
                    stats.MixedJointScratchBytes ==
                        (uint64_t)3u * period * block_bytes &&
                    stats.MixedDualSourceColumns == 0u;
            }
            if (dual_bucketed)
            {
                return stats.MixedDualSourceColumns +
                            stats.InactivatedColumns == column_count &&
                    stats.MixedJointSourceXors == 0u &&
                    stats.MixedJointMarginalXors == 0u &&
                    stats.MixedJointMarginalCopies == 0u &&
                    stats.MixedJointScratchBytes == 0u &&
                    stats.MixedJointActiveDeltas == 0u;
            }
            return stats.MixedJointSourceXors == 0u &&
                stats.MixedJointMarginalXors == 0u &&
                stats.MixedJointMarginalCopies == 0u &&
                stats.MixedJointScratchBytes == 0u &&
                stats.MixedJointActiveDeltas == 0u &&
                stats.MixedDualSourceColumns == 0u;
        };
        if (!bucket_stats_valid(systematic_stats))
        {
            std::fprintf(stderr,
                "solve: mixed bucket systematic accounting failed K=%u\n",
                K);
            return false;
        }
        if (grouped_block_zero_identity)
        {
            std::vector<uint8_t> canonical;
            wirehair_v2::PrecodeSolveStats canonical_stats;
            WirehairResult canonical_result = Wirehair_Error;
            bool canonical_verified = false;
            {
                MixedGroupedGF256RowsScope canonical_scope(0u);
                if (!canonical_scope.IsValid()) return false;
                canonical_result = wirehair_v2::SolvePrecodeSystem(
                    system, config, systematic, block_bytes,
                    canonical, &canonical_stats);
                canonical_verified = canonical_result == Wirehair_Success &&
                    wirehair_v2::VerifyPrecodeSolution(
                        system, config, systematic,
                        canonical.data(), block_bytes);
            }
            if (canonical_result != Wirehair_Success ||
                !canonical_verified || canonical != expected ||
                !SamePrecodeSolveWork(systematic_stats, canonical_stats))
            {
                std::fprintf(stderr,
                    "solve: grouped block-zero systematic identity/work "
                    "mismatch K=%u L-H=%u P=%u result=%d "
                    "xors=%llu/%llu muladds=%llu/%llu\n",
                    K, first_heavy_column, period, (int)canonical_result,
                    (unsigned long long)systematic_stats.BlockXors,
                    (unsigned long long)canonical_stats.BlockXors,
                    (unsigned long long)systematic_stats.BlockMulAdds,
                    (unsigned long long)canonical_stats.BlockMulAdds);
                return false;
            }
        }
        if (joint_bucketed || dual_bucketed)
        {
            std::vector<uint8_t> separate_expected;
            {
                MixedResidueBucketModeScope separate_scope(
                    wirehair_v2::MixedResidueBucketMode::Separate);
                if (!separate_scope.IsValid() ||
                    wirehair_v2::SolvePrecodeSystem(
                        system, config, systematic, block_bytes,
                        separate_expected) != Wirehair_Success)
                {
                    return false;
                }
            }
            if (separate_expected != expected) {
                std::fprintf(stderr,
                    "solve: bucket/separate systematic mismatch K=%u\n", K);
                return false;
            }
        }

        // Replay a deterministic lossy schedule containing repair ids.  This
        // changes the peel graph while preserving one exact expected solution;
        // the enabled hook independently compares every packed inactive
        // coefficient against the original dense projection expansion.
        const size_t delivered_count = (size_t)K + 20u;
        std::vector<uint8_t> delivered_storage(
            delivered_count * block_bytes);
        std::vector<wirehair_v2::SolvePacket> delivered;
        delivered.reserve(delivered_count);
        for (uint32_t id = 0u; delivered.size() < delivered_count; ++id)
        {
            if ((id + (uint32_t)case_index) % 11u == 0u) {
                continue;
            }
            uint8_t* block = delivered_storage.data() +
                delivered.size() * block_bytes;
            if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                    system, config, expected.data(), block_bytes, id, block))
            {
                return false;
            }
            wirehair_v2::SolvePacket packet;
            packet.BlockId = id;
            packet.Data = block;
            delivered.push_back(packet);
        }
        std::vector<uint8_t> recovered;
        wirehair_v2::PrecodeSolveStats repair_stats;
        const uint64_t repair_comparisons_before =
            wirehair_v2::MixedProjectionOracleComparisonsForTesting();
        const WirehairResult repair_result = wirehair_v2::SolvePrecodeSystem(
                system, config, delivered, block_bytes, recovered,
                &repair_stats);
        const uint64_t repair_comparisons_after =
            wirehair_v2::MixedProjectionOracleComparisonsForTesting();
        if (repair_result != Wirehair_Success ||
            repair_comparisons_after != repair_comparisons_before + 1u ||
            recovered != expected ||
            !wirehair_v2::VerifyPrecodeSolution(
                system, config, delivered, recovered.data(), block_bytes))
        {
            std::fprintf(stderr,
                "solve: mixed projection repair oracle failed K=%u\n", K);
            return false;
        }
        if (!bucket_stats_valid(repair_stats))
        {
            std::fprintf(stderr,
                "solve: mixed bucket repair accounting failed K=%u\n", K);
            return false;
        }
        if (grouped_block_zero_identity)
        {
            std::vector<uint8_t> canonical;
            wirehair_v2::PrecodeSolveStats canonical_stats;
            WirehairResult canonical_result = Wirehair_Error;
            bool canonical_verified = false;
            {
                MixedGroupedGF256RowsScope canonical_scope(0u);
                if (!canonical_scope.IsValid()) return false;
                canonical_result = wirehair_v2::SolvePrecodeSystem(
                    system, config, delivered, block_bytes,
                    canonical, &canonical_stats);
                canonical_verified = canonical_result == Wirehair_Success &&
                    wirehair_v2::VerifyPrecodeSolution(
                        system, config, delivered,
                        canonical.data(), block_bytes);
            }
            if (canonical_result != Wirehair_Success ||
                !canonical_verified || canonical != recovered ||
                !SamePrecodeSolveWork(repair_stats, canonical_stats))
            {
                std::fprintf(stderr,
                    "solve: grouped block-zero repair identity/work "
                    "mismatch K=%u L-H=%u P=%u result=%d "
                    "xors=%llu/%llu muladds=%llu/%llu\n",
                    K, first_heavy_column, period, (int)canonical_result,
                    (unsigned long long)repair_stats.BlockXors,
                    (unsigned long long)canonical_stats.BlockXors,
                    (unsigned long long)repair_stats.BlockMulAdds,
                    (unsigned long long)canonical_stats.BlockMulAdds);
                return false;
            }
            std::printf(
                "grouped block-zero identity K=%u L-H=%u P=%u "
                "xors=%llu muladds=%llu: PASS\n",
                K, first_heavy_column, period,
                (unsigned long long)repair_stats.BlockXors,
                (unsigned long long)repair_stats.BlockMulAdds);
        }
        if (joint_bucketed || dual_bucketed)
        {
            std::vector<uint8_t> separate_recovered;
            {
                MixedResidueBucketModeScope separate_scope(
                    wirehair_v2::MixedResidueBucketMode::Separate);
                if (!separate_scope.IsValid() ||
                    wirehair_v2::SolvePrecodeSystem(
                        system, config, delivered, block_bytes,
                        separate_recovered) != Wirehair_Success)
                {
                    return false;
                }
            }
            if (separate_recovered != recovered) {
                std::fprintf(stderr,
                    "solve: bucket/separate repair mismatch K=%u\n", K);
                return false;
            }
        }
    }
    const uint64_t comparisons =
        wirehair_v2::MixedProjectionOracleComparisonsForTesting();
    if (comparisons < 2u * expected_case_count)
    {
        std::fprintf(stderr,
            "solve: mixed projection oracle comparison count=%llu\n",
            (unsigned long long)comparisons);
        return false;
    }
    std::printf(
        "mixed residue-bucket projection oracle period=%u geometry=%u "
        "gf256_rows=%u gf16_rows=%u skew=%u schedule=%u "
        "independent_extension=%u bucket_mode=%u dense_two_anchor=%u "
        "grouped_gf256_rows=%u grouped_hash_seed=0x%x "
        "comparisons=%llu: PASS\n",
        period, (uint32_t)geometry, subfield_rows, extension_rows,
        residue_skew, (uint32_t)residue_schedule,
        independent_extension_residues ? 1u : 0u,
        (uint32_t)bucket_mode,
        dense_two_anchor ? 1u : 0u,
        grouped_gf256_rows,
        wirehair_v2::ActiveMixedGroupedGF256HashSeed(),
        (unsigned long long)comparisons);
    return true;
}

bool CheckMixedProjectionResidueBucketsOracle()
{
    const uint32_t periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u
    };
    const wirehair_v2::MixedCoefficientGeometry geometries[] = {
        wirehair_v2::MixedCoefficientGeometry::FrozenPowerX,
        wirehair_v2::MixedCoefficientGeometry::SharedCauchyX
    };
    for (const wirehair_v2::MixedCoefficientGeometry geometry : geometries) {
        for (const uint32_t period : periods) {
            if (!CheckMixedProjectionResidueBucketsOracleForPeriod(
                    period, geometry, wirehair_v2::kMixedGF16Rows))
            {
                return false;
            }
        }
    }
    const uint32_t h13_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 13u
    };
    const uint32_t h14_periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u, 14u
    };
    for (const uint32_t period : h13_periods) {
        if (!CheckMixedProjectionResidueBucketsOracleForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16Rows + 1u))
        {
            return false;
        }
    }
    for (const uint32_t period : h14_periods) {
        if (!CheckMixedProjectionResidueBucketsOracleForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16RowsMax))
        {
            return false;
        }
    }
    if (!CheckMixedProjectionResidueBucketsOracleForPeriod(
            29u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            14u) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            18u) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
            28u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Ramp) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
            28u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true,
            11u,
            wirehair_v2::MixedResidueBucketMode::Automatic,
            true) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true,
            wirehair_v2::kMixedGF256RowsMax) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
            32u,
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
            wirehair_v2::kMixedGF16RowsMax,
            0u,
            wirehair_v2::MixedResidueSchedule::Hashed,
            true,
            wirehair_v2::kMixedGF256Rows + 1u,
            wirehair_v2::MixedResidueBucketMode::Dual) ||
        !CheckMixedProjectionResidueBucketsOracleForPeriod(
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
    // Keep grouped coverage bounded to the two boundary pairs above: every
    // solve still compares optimized projection against the dense expansion,
    // while explicit separate/joint modes cover P32, P48, and P64 at the
    // minimum, finalist, and maximum useful grouped suffix sizes.
    const uint32_t grouped_periods[] = {32u, 48u, 64u};
    for (const uint32_t period : grouped_periods)
    {
        if (!CheckMixedProjectionResidueBucketsOracleForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16Rows,
                0u,
                wirehair_v2::MixedResidueSchedule::Constant,
                false,
                wirehair_v2::kMixedGF256Rows,
                wirehair_v2::MixedResidueBucketMode::Separate,
                false,
                1u) ||
            !CheckMixedProjectionResidueBucketsOracleForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16Rows,
                0u,
                wirehair_v2::MixedResidueSchedule::Constant,
                false,
                wirehair_v2::kMixedGF256Rows,
                wirehair_v2::MixedResidueBucketMode::Separate,
                false,
                period == 32u ? 7u : 3u) ||
            !CheckMixedProjectionResidueBucketsOracleForPeriod(
                period,
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX,
                wirehair_v2::kMixedGF16Rows,
                0u,
                wirehair_v2::MixedResidueSchedule::Constant,
                false,
                wirehair_v2::kMixedGF256Rows,
                wirehair_v2::MixedResidueBucketMode::JointDelta,
                false,
                9u))
        {
            return false;
        }
    }
    return true;
}

bool CheckMixedMix1EndToEnd()
{
    const uint32_t K = 320u;
    const uint32_t block_bytes = 1280u;
    wirehair_v2::PrecodeParams params = wirehair_v2::MakeMixedParams(
        K, UINT64_C(0x6d697831656e6432));
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x6d697831);
    base_config.MixCount = 1u;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config) != Wirehair_Success ||
        config.MixCount != 1u)
    {
        std::fprintf(stderr, "solve: mixed mix1 configuration failed\n");
        return false;
    }

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 149u + (i >> 5) + 0x31u);
    }
    std::vector<wirehair_v2::SolvePacket> systematic(K);
    for (uint32_t id = 0u; id < K; ++id) {
        systematic[id].BlockId = id;
        systematic[id].Data = message.data() + (size_t)id * block_bytes;
    }
    std::vector<uint8_t> expected;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, systematic, block_bytes, expected) !=
            Wirehair_Success ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, systematic, expected.data(), block_bytes))
    {
        std::fprintf(stderr, "solve: mixed mix1 systematic solve failed\n");
        return false;
    }

    const size_t delivered_count = (size_t)K + 16u;
    std::vector<uint8_t> delivered_storage(
        delivered_count * block_bytes);
    std::vector<wirehair_v2::SolvePacket> delivered;
    delivered.reserve(delivered_count);
    for (uint32_t id = 0u; delivered.size() < delivered_count; ++id)
    {
        if (id % 9u == 4u) {
            continue;
        }
        uint8_t* block = delivered_storage.data() +
            delivered.size() * block_bytes;
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, expected.data(), block_bytes, id, block))
        {
            return false;
        }
        wirehair_v2::SolvePacket packet;
        packet.BlockId = id;
        packet.Data = block;
        delivered.push_back(packet);
    }

    std::vector<uint8_t> recovered;
    wirehair_v2::PrecodeSolveStats stats;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, delivered, block_bytes, recovered, &stats) !=
            Wirehair_Success ||
        recovered != expected ||
        stats.PacketRows != delivered_count ||
        stats.ResidualRank != stats.InactivatedColumns ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, delivered, recovered.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: mixed mix1 repair solve failed rows=%u R=%u rank=%u\n",
            stats.PacketRows, stats.InactivatedColumns, stats.ResidualRank);
        return false;
    }
    std::vector<uint8_t> recovered_message(message.size());
    for (uint32_t id = 0u; id < K; ++id)
    {
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, recovered.data(), block_bytes, id,
                recovered_message.data() + (size_t)id * block_bytes))
        {
            std::fprintf(stderr,
                "solve: mixed mix1 source evaluation failed id=%u\n", id);
            return false;
        }
    }
    if (recovered_message != message)
    {
        std::fprintf(stderr,
            "solve: mixed mix1 recovered source bytes differ\n");
        return false;
    }
    std::printf("mixed mix1 encode/loss/decode/verify: PASS\n");
    return true;
}

bool CheckBinaryPeelLowDegreeXorOracle()
{
    const uint32_t K = 64000u;
    const uint32_t block_bytes = 2u;
    wirehair_v2::PrecodeParams params = wirehair_v2::MakeMixedParams(
        K, UINT64_C(0x7065656c786f7231));
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x786f7231);
    base_config.MixCount = 1u;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "solve: binary peel oracle configuration failed\n");
        return false;
    }

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 173u + (i >> 7) + 0x5bu);
    }
    std::vector<wirehair_v2::SolvePacket> packets(K);
    for (uint32_t id = 0u; id < K; ++id) {
        packets[id].BlockId = id;
        packets[id].Data = message.data() + (size_t)id * block_bytes;
    }

    wirehair_v2::ResetBinaryPeelOracleComparisonsForTesting();
    std::vector<uint8_t> intermediate;
    wirehair_v2::PrecodeSolveStats stats;
    WirehairResult result = Wirehair_Error;
    {
        BinaryPeelOracleScope oracle_scope;
        result = wirehair_v2::SolvePrecodeSystem(
            system, config, packets, block_bytes, intermediate, &stats);
    }
    const uint64_t comparisons =
        wirehair_v2::BinaryPeelOracleComparisonsForTesting();
    const uint32_t L = K + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    if (result != Wirehair_Success || comparisons != 1u ||
        stats.PacketRows != K ||
        stats.PeeledColumns + stats.InactivatedColumns != L ||
        stats.ResidualRank != stats.InactivatedColumns ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, packets, intermediate.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: binary peel oracle failed result=%d comparisons=%llu "
            "peeled=%u inactive=%u rank=%u L=%u\n",
            (int)result, (unsigned long long)comparisons,
            stats.PeeledColumns, stats.InactivatedColumns,
            stats.ResidualRank, L);
        return false;
    }

    std::vector<uint8_t> recovered_message(message.size());
    for (uint32_t id = 0u; id < K; ++id)
    {
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, intermediate.data(), block_bytes, id,
                recovered_message.data() + (size_t)id * block_bytes))
        {
            std::fprintf(stderr,
                "solve: binary peel oracle source evaluation failed id=%u\n",
                id);
            return false;
        }
    }
    if (recovered_message != message)
    {
        std::fprintf(stderr,
            "solve: binary peel oracle source bytes differ\n");
        return false;
    }
    std::printf(
        "K=64000 low-degree-XOR/scan peel oracle comparisons=%llu: PASS\n",
        (unsigned long long)comparisons);
    return true;
}

bool CheckMixedSystematicSolve()
{
    const uint32_t K = 64u;
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x6d697865);
    base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    uint32_t selected_attempt = 0u;
    bool found_nonzero_attempt = false;
    for (uint32_t trial = 0; trial < 128u; ++trial)
    {
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeMixedParams(
                K, UINT64_C(0x6d69786564000000) + trial);
        base_config.PeelSeed = UINT32_C(0x6d697865) + trial;
        if (wirehair_v2::SelectSystematicConfiguration(
                params, base_config, system, config,
                &selected_attempt) == Wirehair_Success &&
            selected_attempt != 0u)
        {
            found_nonzero_attempt = true;
            break;
        }
    }
    if (!found_nonzero_attempt) {
        std::fprintf(stderr,
            "solve: mixed nonzero systematic attempt fixture missing\n");
        return false;
    }

    const uint32_t L = K + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint32_t block_sizes[] = {2u, 16u, 1280u, 2048u, 4096u};
    for (uint32_t block_bytes : block_sizes)
    {
        std::vector<uint8_t> message((size_t)K * block_bytes);
        for (size_t i = 0; i < message.size(); ++i) {
            message[i] = (uint8_t)(i * 113u + block_bytes + selected_attempt);
        }
        std::vector<wirehair_v2::SolvePacket> packets(K);
        for (uint32_t id = 0; id < K; ++id) {
            packets[id].BlockId = id;
            packets[id].Data = message.data() + (size_t)id * block_bytes;
        }
        std::vector<uint8_t> output;
        wirehair_v2::PrecodeSolveStats stats;
        if (wirehair_v2::SolvePrecodeSystem(
                system, config, packets, block_bytes, output, &stats) !=
                Wirehair_Success ||
            stats.ResidualRank != stats.InactivatedColumns ||
            stats.BinaryResidualRank > stats.ResidualRank ||
            stats.ResidualRank - stats.BinaryResidualRank >
                wirehair_v2::kMixedGF256Rows +
                    wirehair_v2::kMixedGF16Rows ||
            stats.SolveValueArenaBytes != (uint64_t)L * block_bytes ||
            stats.SolveValueArenaEagerZeroBytes != 0u ||
            stats.SolveValueArenaCommitCopyBytes !=
                (uint64_t)L * block_bytes ||
            !wirehair_v2::VerifyPrecodeSolution(
                system, config, packets, output.data(), block_bytes))
        {
            std::fprintf(stderr,
                "solve: mixed systematic failed bb=%u q=%u\n",
                block_bytes,
                stats.ResidualRank - stats.BinaryResidualRank);
            return false;
        }

        wirehair_v2::SolveValueStorage lazy_output;
        lazy_output.assign(9u, 0x5au);
        wirehair_v2::PrecodeSolveStats lazy_stats;
        WirehairResult lazy_result = Wirehair_Error;
        {
            SolveValueArenaPoisonScope poison;
            FusedBlockInitializationScope fused;
            lazy_result = wirehair_v2::SolvePrecodeSystem(
                system, config, packets, block_bytes,
                lazy_output, &lazy_stats);
        }
        if (lazy_result != Wirehair_Success ||
            lazy_output.size() != output.size() ||
            !std::equal(
                lazy_output.begin(), lazy_output.end(), output.begin()) ||
            lazy_stats.SolveValueArenaBytes !=
                (uint64_t)L * block_bytes ||
            lazy_stats.SolveValueArenaEagerZeroBytes != 0u ||
            lazy_stats.SolveValueArenaCommitCopyBytes != 0u ||
            lazy_stats.ResidualRank != stats.ResidualRank ||
            lazy_stats.BinaryResidualRank != stats.BinaryResidualRank)
        {
            std::fprintf(stderr,
                "solve: poisoned mixed no-init solve failed bb=%u\n",
                block_bytes);
            return false;
        }
    }

    // Exercise exact mixed quotient widths through the production packet
    // projection.  Consistent repair equations reduce the 12-column base
    // quotient one binary rank at a time.
    const uint32_t boundary_bytes = 16u;
    std::vector<uint8_t> boundary_message((size_t)K * boundary_bytes);
    for (size_t i = 0; i < boundary_message.size(); ++i) {
        boundary_message[i] = (uint8_t)(i * 71u + selected_attempt);
    }
    std::vector<wirehair_v2::SolvePacket> boundary_packets(K);
    for (uint32_t id = 0; id < K; ++id) {
        boundary_packets[id].BlockId = id;
        boundary_packets[id].Data =
            boundary_message.data() + (size_t)id * boundary_bytes;
    }
    std::vector<uint8_t> boundary_intermediate;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, boundary_packets, boundary_bytes,
            boundary_intermediate) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> mixed_oracle_output(7u, 0xa5u);
    const std::vector<uint8_t> mixed_oracle_before = mixed_oracle_output;
    if (wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
            system, config, boundary_packets, boundary_bytes,
            mixed_oracle_output) != Wirehair_InvalidInput ||
        mixed_oracle_output != mixed_oracle_before)
    {
        std::fprintf(stderr,
            "solve: GF256 tiny oracle accepted mixed coefficients\n");
        return false;
    }
    static const uint32_t target_q[] = {0u, 1u, 2u, 10u, 12u};
    bool saw_q[sizeof(target_q) / sizeof(target_q[0])] = {};
    std::vector<std::vector<uint8_t> > repair_blocks(
        20u, std::vector<uint8_t>(boundary_bytes));
    for (uint32_t i = 0; i < repair_blocks.size(); ++i)
    {
        uint64_t operations = 0u;
        const uint32_t id = K + i;
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, boundary_intermediate.data(),
                boundary_bytes, id, repair_blocks[i].data(), &operations) ||
            operations == 0u)
        {
            return false;
        }
    }
    for (uint32_t overhead = 0u; overhead <= repair_blocks.size();
         ++overhead)
    {
        if (overhead != 0u) {
            wirehair_v2::SolvePacket packet;
            packet.BlockId = K + overhead - 1u;
            packet.Data = repair_blocks[overhead - 1u].data();
            boundary_packets.push_back(packet);
        }
        std::vector<uint8_t> output(9u, 0xabu);
        wirehair_v2::PrecodeSolveStats stats;
        const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
            system, config, boundary_packets, boundary_bytes,
            output, &stats);
        const uint32_t q =
            stats.InactivatedColumns - stats.BinaryResidualRank;
        if (result != Wirehair_Success ||
            output != boundary_intermediate)
        {
            std::fprintf(stderr,
                "solve: mixed quotient overhead=%u failed result=%d q=%u\n",
                overhead, (int)result, q);
            return false;
        }
        for (uint32_t t = 0;
             t < sizeof(target_q) / sizeof(target_q[0]); ++t) {
            if (q == target_q[t]) saw_q[t] = true;
        }
    }
    for (uint32_t t = 0;
         t < sizeof(target_q) / sizeof(target_q[0]); ++t) {
        if (!saw_q[t]) {
            std::fprintf(stderr,
                "solve: mixed quotient q=%u fixture missing\n", target_q[t]);
            return false;
        }
    }

    for (uint32_t corrupt_byte = 0u; corrupt_byte < 2u; ++corrupt_byte)
    {
        std::vector<wirehair_v2::SolvePacket> inconsistent(
            boundary_packets.begin(), boundary_packets.begin() + K + 1u);
        std::vector<uint8_t> corrupt_repair = repair_blocks[0];
        corrupt_repair[corrupt_byte] ^= 1u;
        inconsistent.back().Data = corrupt_repair.data();
        std::vector<uint8_t> inconsistent_output(9u, 0xedu);
        const std::vector<uint8_t> inconsistent_before =
            inconsistent_output;
        wirehair_v2::PrecodeSolveStats inconsistent_stats;
        if (wirehair_v2::SolvePrecodeSystem(
                system, config, inconsistent, boundary_bytes,
                inconsistent_output, &inconsistent_stats) !=
                Wirehair_Error ||
            inconsistent_stats.InactivatedColumns -
                inconsistent_stats.BinaryResidualRank != 11u ||
            inconsistent_output != inconsistent_before)
        {
            std::fprintf(stderr,
                "solve: mixed zero-coeff inconsistency byte=%u failed\n",
                corrupt_byte);
            return false;
        }

        wirehair_v2::SolveValueStorage lazy_inconsistent;
        lazy_inconsistent.assign(9u, 0xedu);
        WirehairResult lazy_inconsistent_result = Wirehair_NeedMore;
        {
            SolveValueArenaPoisonScope poison;
            FusedBlockInitializationScope fused;
            lazy_inconsistent_result = wirehair_v2::SolvePrecodeSystem(
                system, config, inconsistent, boundary_bytes,
                lazy_inconsistent);
        }
        if (lazy_inconsistent_result != Wirehair_Error ||
            lazy_inconsistent.size() != inconsistent_before.size() ||
            !std::equal(
                lazy_inconsistent.begin(), lazy_inconsistent.end(),
                inconsistent_before.begin()))
        {
            std::fprintf(stderr,
                "solve: poisoned mixed inconsistency byte=%u failed\n",
                corrupt_byte);
            return false;
        }
    }
    for (uint32_t corrupt_byte = 0u; corrupt_byte < 2u; ++corrupt_byte)
    {
        std::vector<uint8_t> corrupt_intermediate = boundary_intermediate;
        corrupt_intermediate[corrupt_byte] ^= 1u;
        if (wirehair_v2::VerifyPrecodeSolution(
                system, config, boundary_packets,
                corrupt_intermediate.data(), boundary_bytes))
        {
            std::fprintf(stderr,
                "solve: mixed verifier accepted byte=%u corruption\n",
                corrupt_byte);
            return false;
        }
    }

    std::vector<wirehair_v2::SolvePacket> q13_packets(
        boundary_packets.begin(), boundary_packets.begin() + (K - 1u));
    q13_packets.push_back(boundary_packets[0]);
    std::vector<uint8_t> q13_output(9u, 0xcdu);
    const std::vector<uint8_t> q13_before = q13_output;
    wirehair_v2::PrecodeSolveStats q13_stats;
    wirehair_v2::PrecodeSolveResumeState q13_resume;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, q13_packets, boundary_bytes,
            q13_output, &q13_stats, &q13_resume) != Wirehair_NeedMore ||
        q13_stats.InactivatedColumns - q13_stats.BinaryResidualRank != 13u ||
        q13_stats.ResidualRank != q13_stats.BinaryResidualRank ||
        q13_output != q13_before || q13_resume.Active)
    {
        std::fprintf(stderr,
            "solve: mixed q13 boundary failed q=%u\n",
            q13_stats.InactivatedColumns - q13_stats.BinaryResidualRank);
        return false;
    }
    wirehair_v2::PrecodeSolveResumeState prior_resume;
    prior_resume.Active = true;
    prior_resume.SourceCount = UINT32_C(0x12345678);
    prior_resume.CoefficientScratch.assign(3u, 0xa5u);
    std::vector<uint8_t> prior_output = q13_before;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, q13_packets, boundary_bytes,
            prior_output, nullptr, &prior_resume) != Wirehair_NeedMore ||
        !prior_resume.Active ||
        prior_resume.SourceCount != UINT32_C(0x12345678) ||
        prior_resume.CoefficientScratch !=
            std::vector<uint8_t>(3u, 0xa5u) ||
        prior_output != q13_before)
    {
        std::fprintf(stderr,
            "solve: mixed q13 changed caller resume/output state\n");
        return false;
    }

    std::vector<uint8_t> odd_output(17u, 0xa5u);
    const std::vector<uint8_t> odd_before = odd_output;
    uint8_t odd_data[3] = {};
    std::vector<wirehair_v2::SolvePacket> odd_packets(K);
    for (wirehair_v2::SolvePacket& packet : odd_packets) {
        packet.BlockId = 0u;
        packet.Data = odd_data;
    }
    wirehair_v2::PrecodeSolveResumeState resume;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, odd_packets, 3u, odd_output, nullptr, &resume) !=
            Wirehair_InvalidInput ||
        odd_output != odd_before || resume.Active)
    {
        std::fprintf(stderr, "solve: mixed odd rejection failed\n");
        return false;
    }

    std::vector<uint8_t> deficient_output(17u, 0x5au);
    const std::vector<uint8_t> deficient_before = deficient_output;
    uint8_t duplicate[16] = {};
    std::vector<wirehair_v2::SolvePacket> deficient(K);
    for (wirehair_v2::SolvePacket& packet : deficient) {
        packet.BlockId = 0u;
        packet.Data = duplicate;
    }
    wirehair_v2::PrecodeSolveStats deficient_stats;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, deficient, sizeof(duplicate),
            deficient_output, &deficient_stats, &resume) !=
            Wirehair_NeedMore ||
        deficient_output != deficient_before || resume.Active ||
        wirehair_v2::ResumePrecodeSystem(
            system, config, 1u, duplicate, sizeof(duplicate),
            resume, deficient_output) != Wirehair_InvalidInput ||
        deficient_output != deficient_before || resume.Active)
    {
        std::fprintf(stderr,
            "solve: mixed deficient/no-resume contract failed\n");
        return false;
    }
    wirehair_v2::SolveValueStorage lazy_deficient;
    lazy_deficient.assign(17u, 0x5au);
    WirehairResult lazy_deficient_result = Wirehair_Error;
    {
        SolveValueArenaPoisonScope poison;
        FusedBlockInitializationScope fused;
        lazy_deficient_result = wirehair_v2::SolvePrecodeSystem(
            system, config, deficient, sizeof(duplicate), lazy_deficient);
    }
    if (lazy_deficient_result != Wirehair_NeedMore ||
        lazy_deficient.size() != deficient_before.size() ||
        !std::equal(
            lazy_deficient.begin(), lazy_deficient.end(),
            deficient_before.begin()))
    {
        std::fprintf(stderr,
            "solve: poisoned mixed deficient output changed\n");
        return false;
    }
    std::printf(
        "mixed systematic solve/nonzero attempt=%u/no-resume: PASS\n",
        selected_attempt);
    return true;
}

bool CheckPackedBinaryResidualOracle()
{
    if (!wirehair_v2::CheckPackedBinaryResidualOracleForTesting())
    {
        std::fprintf(stderr,
            "solve: packed GF2 residual differential oracle failed\n");
        return false;
    }
    std::printf(
        "packed GF2 residual R=63..193 word-boundary differential: PASS\n");
    return true;
}

bool CheckMixedRhsFusionOracle()
{
    if (!wirehair_v2::CheckMixedRhsFusionOracleForTesting())
    {
        std::fprintf(stderr,
            "solve: mixed RHS initialization fusion oracle failed\n");
        return false;
    }
    std::printf(
        "mixed q=0..16 pivot/residual RHS fusion oracle: PASS\n");
    return true;
}

bool CheckMixedNullWitnessCanonicalization()
{
    if (!wirehair_v2::CheckMixedNullWitnessCanonicalizationForTesting())
    {
        std::fprintf(stderr,
            "solve: mixed null-witness canonicalization oracle failed\n");
        return false;
    }
    std::printf(
        "mixed GF16 null-witness basis invariance: PASS\n");
    return true;
}

bool CheckMixedQuotientRankFirstOracles()
{
    if (!wirehair_v2::CheckMixedQuotientFactorReplayForTesting())
    {
        std::fprintf(stderr,
            "solve: mixed quotient factor/replay oracle failed\n");
        return false;
    }
    {
        MixedCoefficientPeriodScope period(32u);
        MixedCoefficientGeometryScope geometry(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
        MixedResidueScheduleScope schedule(
            wirehair_v2::MixedResidueSchedule::Constant);
        MixedIndependentExtensionResiduesScope independent(false);
        MixedResidueBucketModeScope buckets(
            wirehair_v2::MixedResidueBucketMode::Separate);
        if (!period.IsValid() || !geometry.IsValid() ||
            !schedule.IsValid() || !independent.IsValid() ||
            !buckets.IsValid() ||
            !wirehair_v2::
                CheckMixedQuotientDeficientSyndromeForTesting())
        {
            std::fprintf(stderr,
                "solve: constant mixed deficient-syndrome oracle failed\n");
            return false;
        }
    }
    {
        MixedCoefficientPeriodScope period(32u);
        MixedCoefficientGeometryScope geometry(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
        MixedResidueScheduleScope schedule(
            wirehair_v2::MixedResidueSchedule::Hashed);
        MixedResidueHashSeedScope hash_seed(68u);
        MixedIndependentExtensionResiduesScope independent(true);
        MixedResidueBucketModeScope buckets(
            wirehair_v2::MixedResidueBucketMode::Separate);
        if (!period.IsValid() || !geometry.IsValid() ||
            !schedule.IsValid() || !independent.IsValid() ||
            !buckets.IsValid() ||
            !wirehair_v2::
                CheckMixedQuotientDeficientSyndromeForTesting())
        {
            std::fprintf(stderr,
                "solve: independent mixed deficient-syndrome oracle failed\n");
            return false;
        }
    }
    {
        MixedCoefficientPeriodScope period(32u);
        MixedCoefficientGeometryScope geometry(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
        MixedResidueScheduleScope schedule(
            wirehair_v2::MixedResidueSchedule::Constant);
        MixedIndependentExtensionResiduesScope independent(false);
        MixedResidueBucketModeScope buckets(
            wirehair_v2::MixedResidueBucketMode::Separate);
        // Configure the grouped suffix last: prerequisite setters deliberately
        // clear experiment schedules on the current thread.
        MixedGroupedGF256RowsScope grouped(3u);
        if (!period.IsValid() || !geometry.IsValid() ||
            !schedule.IsValid() || !independent.IsValid() ||
            !buckets.IsValid() || !grouped.IsValid() ||
            !wirehair_v2::
                CheckMixedQuotientDeficientSyndromeForTesting())
        {
            std::fprintf(stderr,
                "solve: grouped GF256 mixed deficient-syndrome oracle "
                "failed\n");
            return false;
        }
    }
    std::printf(
        "mixed quotient factor/replay and constant/independent/grouped "
        "deficient syndromes: PASS\n");
    return true;
}

bool RunCase(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t loss_stride,
    wirehair_v2::HeavyCoefficientFamily heavy_family =
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy)
{
    const wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    wirehair_v2::PrecodeParams params = wirehair_v2::MakeCertifiedParams(
        K,
        wirehair_v2::MatrixSeedFromProfile(
            profile, 0u, wirehair_v2::kMessagePrecodeSeedSalt));
    params.HeavyFamily = heavy_family;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        std::fprintf(stderr, "solve: precode build failed K=%u\n", K);
        return false;
    }

    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
        profile, wirehair_v2::kMessageRecoveryRowSeedSalt);
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 131u + K * 17u + 29u);
    }
    std::vector<wirehair_v2::SolvePacket> packets;
    for (uint32_t id = 0; id < K; ++id)
    {
        wirehair_v2::SolvePacket packet;
        packet.BlockId = id;
        packet.Data = message.data() + (size_t)id * block_bytes;
        packets.push_back(packet);
    }

    std::vector<uint8_t> intermediate;
    wirehair_v2::PrecodeSolveStats stats;
    const WirehairResult encoded = wirehair_v2::SolvePrecodeSystem(
        system, config, packets, block_bytes, intermediate, &stats);
    if (encoded != Wirehair_Success ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, packets, intermediate.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: systematic solve failed K=%u result=%d R=%u rank=%u\n",
            K, (int)encoded, stats.InactivatedColumns, stats.ResidualRank);
        return false;
    }
    const uint64_t row_count = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + packets.size();
    const uint64_t column_count = (uint64_t)K +
        system.Params.Staircase + system.Params.DenseRows +
        system.Params.HeavyRows;
    const uint64_t old_allocation_lower_bound =
        row_count + column_count + 2u;
    const uint64_t pooled_allocations =
        stats.BinaryRowStorageAllocations +
        stats.BinaryAdjacencyStorageAllocations;
    const uint64_t old_peak_bytes_lower_bound =
        row_count * (sizeof(std::vector<uint32_t>) + sizeof(const uint8_t*)) +
        column_count * sizeof(std::vector<uint32_t>) +
        2u * stats.BinaryRowReferences * sizeof(uint32_t);
    const uint64_t pooled_peak_bytes =
        stats.BinaryRowStorageBytes + stats.BinaryAdjacencyStorageBytes;
    if (stats.BinaryRowStorageAllocations != 3u ||
        stats.BinaryAdjacencyStorageAllocations != 2u ||
        pooled_allocations >= old_allocation_lower_bound ||
        pooled_peak_bytes >= old_peak_bytes_lower_bound)
    {
        std::fprintf(stderr,
            "solve: pooled binary storage regression K=%u alloc=%llu/%llu "
            "bytes=%llu/%llu\n",
            K,
            (unsigned long long)pooled_allocations,
            (unsigned long long)old_allocation_lower_bound,
            (unsigned long long)pooled_peak_bytes,
            (unsigned long long)old_peak_bytes_lower_bound);
        return false;
    }

    std::vector<uint8_t> block(block_bytes);
    for (uint32_t id = 0; id < K; ++id)
    {
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, intermediate.data(), block_bytes,
                id, block.data()) ||
            std::memcmp(
                block.data(),
                message.data() + (size_t)id * block_bytes,
                block_bytes) != 0)
        {
            std::fprintf(stderr,
                "solve: systematic row mismatch K=%u id=%u\n", K, id);
            return false;
        }
    }

    std::vector<std::vector<uint8_t> > delivered_data;
    std::vector<wirehair_v2::SolvePacket> delivered;
    delivered_data.reserve(K + 32u);
    delivered.reserve(K + 32u);
    for (uint32_t id = K; id-- > 0u;)
    {
        if (id % loss_stride == 0u) {
            continue;
        }
        delivered_data.push_back(std::vector<uint8_t>(block_bytes));
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, intermediate.data(), block_bytes,
                id, delivered_data.back().data()))
        {
            return false;
        }
        wirehair_v2::SolvePacket packet;
        packet.BlockId = id;
        packet.Data = delivered_data.back().data();
        delivered.push_back(packet);
    }
    WirehairResult decoded = Wirehair_NeedMore;
    std::vector<uint8_t> recovered_intermediate;
    uint32_t repair_id = K;
    while (repair_id < K + 32u && decoded != Wirehair_Success)
    {
        delivered_data.push_back(std::vector<uint8_t>(block_bytes));
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, intermediate.data(), block_bytes,
                repair_id, delivered_data.back().data()))
        {
            return false;
        }
        wirehair_v2::SolvePacket packet;
        packet.BlockId = repair_id++;
        packet.Data = delivered_data.back().data();
        delivered.push_back(packet);
        if (delivered.size() >= K) {
            decoded = wirehair_v2::SolvePrecodeSystem(
                system, config, delivered, block_bytes,
                recovered_intermediate, &stats);
        }
    }
    if (decoded != Wirehair_Success ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, delivered,
            recovered_intermediate.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: lossy solve failed K=%u result=%d delivered=%zu\n",
            K, (int)decoded, delivered.size());
        return false;
    }
    for (uint32_t id = 0; id < K; ++id)
    {
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, recovered_intermediate.data(), block_bytes,
                id, block.data()) ||
            std::memcmp(
                block.data(),
                message.data() + (size_t)id * block_bytes,
                block_bytes) != 0)
        {
            std::fprintf(stderr,
                "solve: recovered message mismatch K=%u id=%u\n", K, id);
            return false;
        }
    }

    recovered_intermediate[0] ^= 1u;
    if (wirehair_v2::VerifyPrecodeSolution(
            system, config, delivered,
            recovered_intermediate.data(), block_bytes))
    {
        std::fprintf(stderr, "solve: corrupted solution verified K=%u\n", K);
        return false;
    }

    std::printf(
        "global solve K=%u bb=%u family=%u delivered=%zu inact=%u "
        "rank=%u binary_storage_alloc=%llu/%llu bytes=%llu/%llu: PASS\n",
        K, block_bytes, (unsigned)heavy_family, delivered.size(),
        stats.InactivatedColumns, stats.ResidualRank,
        (unsigned long long)pooled_allocations,
        (unsigned long long)old_allocation_lower_bound,
        (unsigned long long)pooled_peak_bytes,
        (unsigned long long)old_peak_bytes_lower_bound);
    return true;
}

bool CheckMixDomainValidation()
{
    wirehair_v2::PrecodeParams params;
    params.BlockCount = 2u;
    params.Staircase = 2u;
    params.DenseRows = 0u;
    params.HeavyRows = 0u;
    params.SourceHits = 1u;
    params.Seed = 7u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        return false;
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = 11u;
    config.MixCount = 3u;
    if (!wirehair_v2::GeneratePacketMatrixRow(2u, 2u, 0u, config).empty()) {
        std::fprintf(stderr, "solve: oversized mix generated duplicates\n");
        return false;
    }

    const uint8_t intermediate[4] = {1u, 2u, 3u, 4u};
    uint8_t output = 0xa5u;
    if (wirehair_v2::EvaluatePacketBlock(
            system, config, intermediate, 1u, 0u, &output) ||
        output != 0xa5u)
    {
        std::fprintf(stderr, "solve: invalid mix modified packet output\n");
        return false;
    }
    const uint8_t zero = 0u;
    std::vector<wirehair_v2::SolvePacket> packets(2u);
    packets[0].BlockId = 0u;
    packets[0].Data = &zero;
    packets[1].BlockId = 1u;
    packets[1].Data = &zero;
    std::vector<uint8_t> solved(3u, 0xccu);
    const std::vector<uint8_t> before = solved;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, packets, 1u, solved) != Wirehair_InvalidInput ||
        solved != before)
    {
        std::fprintf(stderr, "solve: invalid mix solve was not no-write\n");
        return false;
    }
    const uint8_t zero_intermediate[4] = {};
    if (wirehair_v2::VerifyPrecodeSolution(
            system, config, packets, zero_intermediate, 1u))
    {
        std::fprintf(stderr,
            "solve: oversized mix verified an empty packet equation\n");
        return false;
    }
    config.MixCount = 2u;
    const std::vector<uint32_t> boundary_row =
        wirehair_v2::GeneratePacketMatrixRow(2u, 2u, 0u, config);
    const std::vector<uint8_t> boundary_expected = ReferencePacket(
        2u, 2u, 0u, config, intermediate, 1u);
    uint64_t boundary_operations = 0u;
    output = 0xa5u;
    if (boundary_row.empty() || boundary_expected.size() != 1u ||
        !wirehair_v2::EvaluatePacketBlock(
            system, config, intermediate, 1u, 0u,
            &output, &boundary_operations) ||
        output != boundary_expected[0] ||
        boundary_operations != boundary_row.size())
    {
        std::fprintf(stderr,
            "solve: P=2 two-mix fused boundary mismatch\n");
        return false;
    }
    config.MixCount = 0u;
    if (wirehair_v2::VerifyPrecodeSolution(
            system, config, packets, zero_intermediate, 1u))
    {
        std::fprintf(stderr,
            "solve: zero mix verified an empty packet equation\n");
        return false;
    }
    return true;
}

bool CheckConcurrentMixedCoefficientTable()
{
    static const uint32_t kThreadCount = 16u;
    const wirehair_v2::MixedPackedCoefficients* tables[kThreadCount] = {};
    std::vector<std::thread> workers;
    workers.reserve(kThreadCount);
    std::atomic<uint32_t> ready(0u);
    std::atomic<bool> start(false);
    try
    {
        for (uint32_t thread = 0u; thread < kThreadCount; ++thread)
        {
            workers.push_back(std::thread([&, thread]() {
                ready.fetch_add(1u, std::memory_order_release);
                while (!start.load(std::memory_order_acquire)) {
                    std::this_thread::yield();
                }
                tables[thread] =
                    wirehair_v2::GetMixedPackedCoefficients();
            }));
        }
    }
    catch (...)
    {
        start.store(true, std::memory_order_release);
        for (std::thread& worker : workers) {
            worker.join();
        }
        std::fprintf(stderr,
            "solve: concurrent mixed table thread launch failed\n");
        return false;
    }
    while (ready.load(std::memory_order_acquire) != kThreadCount) {
        std::this_thread::yield();
    }
    start.store(true, std::memory_order_release);
    for (std::thread& worker : workers) {
        worker.join();
    }
    for (uint32_t thread = 0u; thread < kThreadCount; ++thread) {
        if (!tables[thread] || tables[thread] != tables[0]) {
            std::fprintf(stderr,
                "solve: concurrent mixed table publication mismatch\n");
            return false;
        }
    }

    const wirehair_v2::MixedCoefficientRows* rows =
        wirehair_v2::GetMixedCoefficientRows();
    const uint32_t H =
        wirehair_v2::kMixedGF256Rows + wirehair_v2::kMixedGF16Rows;
    if (!rows) {
        return false;
    }
    for (uint32_t residue = 0u;
         residue < wirehair_v2::kMixedCoefficientPeriod;
         ++residue)
    {
        for (uint32_t row = 0u; row < H; ++row)
        {
            const uint16_t expected =
                row < wirehair_v2::kMixedGF256Rows ?
                    wirehair_v2::HeavyCoefficient(row, residue, H) :
                    wirehair_v2::MixedGF16Coefficient(
                        row - wirehair_v2::kMixedGF256Rows, residue);
            const uint16_t row_value =
                row < wirehair_v2::kMixedGF256Rows ?
                    rows->Subfield[row][residue] :
                    rows->Extension[
                        row - wirehair_v2::kMixedGF256Rows][residue];
            const uint16_t packed_value = (uint16_t)(
                tables[0]->ByResidue[residue][row >> 2] >>
                ((row & 3u) * 16u));
            if (row_value != expected || packed_value != expected)
            {
                std::fprintf(stderr,
                    "solve: mixed coefficient cache mismatch r=%u m=%u\n",
                    row, residue);
                return false;
            }
        }
    }
    return true;
}

bool CheckConcurrentCoefficientCacheCase(bool mixed)
{
    const uint32_t K = 320u;
    const uint32_t block_bytes = mixed ? 38u : 37u;
    const wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    const uint64_t matrix_seed = wirehair_v2::MatrixSeedFromProfile(
        profile, 0u, wirehair_v2::kMessagePrecodeSeedSalt);
    wirehair_v2::PrecodeParams params = mixed ?
        wirehair_v2::MakeMixedParams(K, matrix_seed) :
        wirehair_v2::MakeCertifiedParams(K, matrix_seed);
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        return false;
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
        profile, wirehair_v2::kMessageRecoveryRowSeedSalt);
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 89u + (i >> 3) + 17u);
    }
    std::vector<wirehair_v2::SolvePacket> packets(K);
    for (uint32_t id = 0u; id < K; ++id) {
        packets[id].BlockId = id;
        packets[id].Data = message.data() + (size_t)id * block_bytes;
    }

    static const uint32_t kThreadCount = 16u;
    std::vector<std::vector<uint8_t> > outputs(kThreadCount);
    std::vector<std::thread> workers;
    workers.reserve(kThreadCount);
    std::atomic<uint32_t> ready(0u);
    std::atomic<bool> start(false);
    std::atomic<bool> failed(false);
    try
    {
        for (uint32_t thread = 0u; thread < kThreadCount; ++thread)
        {
            workers.push_back(std::thread([&, thread]() {
                ready.fetch_add(1u, std::memory_order_release);
                while (!start.load(std::memory_order_acquire)) {
                    std::this_thread::yield();
                }
                wirehair_v2::PrecodeSolveStats stats;
                if (wirehair_v2::SolvePrecodeSystem(
                        system, config, packets, block_bytes,
                        outputs[thread], &stats) != Wirehair_Success ||
                    !wirehair_v2::VerifyPrecodeSolution(
                        system, config, packets,
                        outputs[thread].data(), block_bytes))
                {
                    failed.store(true, std::memory_order_relaxed);
                }
            }));
        }
    }
    catch (...)
    {
        start.store(true, std::memory_order_release);
        for (std::thread& worker : workers) {
            worker.join();
        }
        std::fprintf(stderr,
            "solve: concurrent %s cache thread launch failed\n",
            mixed ? "mixed" : "H12");
        return false;
    }
    while (ready.load(std::memory_order_acquire) != kThreadCount) {
        std::this_thread::yield();
    }
    start.store(true, std::memory_order_release);
    for (std::thread& worker : workers) {
        worker.join();
    }
    if (failed.load(std::memory_order_relaxed)) {
        std::fprintf(stderr,
            "solve: concurrent %s coefficient-cache solve failed\n",
            mixed ? "mixed" : "H12");
        return false;
    }
    for (uint32_t thread = 1u; thread < kThreadCount; ++thread)
    {
        if (outputs[thread] != outputs[0]) {
            std::fprintf(stderr,
                "solve: concurrent %s coefficient-cache mismatch\n",
                mixed ? "mixed" : "H12");
            return false;
        }
    }
    return true;
}

bool CheckConcurrentCoefficientCaches()
{
    // Publish both mixed tables under contention, verify every coefficient
    // against its independent generator, then exercise cached mixed solving.
    // The certified case independently covers first use of the H12 table.
    if (!CheckConcurrentMixedCoefficientTable() ||
        !CheckConcurrentCoefficientCacheCase(true) ||
        !CheckConcurrentCoefficientCacheCase(false))
    {
        return false;
    }
    std::printf(
        "concurrent mixed/H12 coefficient-cache first use: PASS\n");
    return true;
}

bool PacketRowHasDistinctMix(
    const std::vector<uint32_t>& row,
    uint32_t K,
    uint32_t P,
    uint32_t mix_count)
{
    if (row.size() <= mix_count) {
        return false;
    }
    const size_t mix_begin = row.size() - mix_count;
    for (size_t i = 0; i < mix_begin; ++i) {
        if (row[i] >= K) {
            return false;
        }
        for (size_t j = 0; j < i; ++j) {
            if (row[j] == row[i]) {
                return false;
            }
        }
    }
    for (size_t i = mix_begin; i < row.size(); ++i)
    {
        if (row[i] < K || (uint64_t)row[i] >= (uint64_t)K + P) {
            return false;
        }
        for (size_t j = mix_begin; j < i; ++j) {
            if (row[j] == row[i]) {
                return false;
            }
        }
    }
    return true;
}

bool CheckPacketRowDomainBoundaries()
{
    static_assert(
        wirehair_v2::kMaxPacketPrecodeCount == 65521u,
        "packet domain must track the last 16-bit prime");
    struct BoundaryCase
    {
        uint32_t K;
        uint32_t P;
        uint32_t MixCount;
        bool PacketValid;
    };
    const BoundaryCase cases[] = {
        { 2u, 1u, 1u, false },
        { 2u, 2u, 2u, true },
        { 2u, 65521u, 3u, true },
        { 2u, 65522u, 3u, false },
        // Both sides of the maximum uint16 structural-span boundary.
        { 14u, 65521u, 3u, true },
        { 2u, 65533u, 3u, false }
    };
    const uint32_t sample_ids[] = {
        0u, 1u, 13u, 65521u, 153334u, UINT32_MAX
    };

    for (const BoundaryCase& test : cases)
    {
        wirehair_v2::PrecodeParams params;
        params.BlockCount = test.K;
        params.Staircase = test.P;
        params.DenseRows = 0u;
        params.HeavyRows = 0u;
        params.SourceHits = 1u;
        params.Seed = UINT64_C(0x5041434b4554444f) ^ test.P;
        wirehair_v2::PrecodeSystem system;
        if (!wirehair_v2::BuildPrecodeSystem(params, system) ||
            !wirehair_v2::ValidatePrecodeSystem(system))
        {
            std::fprintf(stderr,
                "solve: structure boundary rejected K=%u P=%u\n",
                test.K, test.P);
            return false;
        }

        wirehair_v2::PacketRowConfig config;
        // Seed/id pair pins the original P=65522 duplicate-column
        // reproducer in Release builds; rejection now precedes generation.
        config.PeelSeed = 1u;
        config.MixCount = test.MixCount;
        if (wirehair_v2::IsPacketRowDomainValid(
                test.K, test.P, test.MixCount) != test.PacketValid)
        {
            std::fprintf(stderr,
                "solve: packet-domain predicate mismatch K=%u P=%u\n",
                test.K, test.P);
            return false;
        }

        for (uint32_t id : sample_ids)
        {
            const std::vector<uint32_t> row =
                wirehair_v2::GeneratePacketMatrixRow(
                    test.K, test.P, id, config);
            if (test.PacketValid ?
                    !PacketRowHasDistinctMix(
                        row, test.K, test.P, test.MixCount) :
                    !row.empty())
            {
                std::fprintf(stderr,
                    "solve: packet row boundary mismatch K=%u P=%u id=%u\n",
                    test.K, test.P, id);
                return false;
            }
        }

        std::vector<uint8_t> intermediate(
            (size_t)test.K + test.P, 0u);
        for (size_t i = 0; i < intermediate.size(); ++i) {
            intermediate[i] = (uint8_t)(i * 29u + test.P * 7u + 3u);
        }
        uint8_t output = 0xa5u;
        uint64_t operations = UINT64_C(0xfeedfacecafebeef);
        const uint32_t evaluation_id = 153334u;
        const bool evaluated = wirehair_v2::EvaluatePacketBlock(
            system, config, intermediate.data(), 1u, evaluation_id,
            &output, &operations);
        if (evaluated != test.PacketValid)
        {
            std::fprintf(stderr,
                "solve: packet evaluation boundary mismatch K=%u P=%u\n",
                test.K, test.P);
            return false;
        }
        if (test.PacketValid)
        {
            const std::vector<uint32_t> row =
                wirehair_v2::GeneratePacketMatrixRow(
                    test.K, test.P, evaluation_id, config);
            uint8_t expected = 0u;
            for (uint32_t column : row) {
                expected ^= intermediate[column];
            }
            if (output != expected || operations != row.size()) {
                std::fprintf(stderr,
                    "solve: accepted packet evaluation mismatch K=%u P=%u\n",
                    test.K, test.P);
                return false;
            }
        }
        else
        {
            if (output != 0xa5u ||
                operations != UINT64_C(0xfeedfacecafebeef))
            {
                std::fprintf(stderr,
                    "solve: rejected packet evaluation wrote output K=%u "
                    "P=%u\n", test.K, test.P);
                return false;
            }
            output = 0x5au;
            operations = UINT64_C(0x0123456789abcdef);
            if (wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                    system, config, intermediate.data(), 1u, evaluation_id,
                    &output, &operations) ||
                output != 0x5au ||
                operations != UINT64_C(0x0123456789abcdef))
            {
                std::fprintf(stderr,
                    "solve: rejected fast packet evaluation wrote output "
                    "K=%u P=%u\n", test.K, test.P);
                return false;
            }

            const uint8_t zero = 0u;
            std::vector<wirehair_v2::SolvePacket> packets(test.K);
            for (uint32_t id = 0; id < test.K; ++id) {
                packets[id].BlockId = id;
                packets[id].Data = &zero;
            }
            std::vector<uint8_t> solved(7u, 0xccu);
            const std::vector<uint8_t> before = solved;
            wirehair_v2::PrecodeSolveStats stats;
            stats.PacketRows = UINT32_MAX;
            if (wirehair_v2::SolvePrecodeSystem(
                    system, config, packets, 1u, solved, &stats) !=
                    Wirehair_InvalidInput ||
                solved != before || stats.PacketRows != UINT32_MAX ||
                wirehair_v2::VerifyPrecodeSolution(
                    system, config, packets, intermediate.data(), 1u))
            {
                std::fprintf(stderr,
                    "solve: rejected packet solve/verify contract failed "
                    "K=%u P=%u\n", test.K, test.P);
                return false;
            }

            wirehair_v2::PacketRowConfig selected;
            selected.PeelSeed = UINT32_C(0xdecafbad);
            selected.MixCount = UINT32_MAX;
            uint32_t attempt = UINT32_MAX;
            if (wirehair_v2::SelectSystematicPacketConfig(
                    system, config, selected, &attempt) !=
                    Wirehair_InvalidInput ||
                selected.PeelSeed != UINT32_C(0xdecafbad) ||
                selected.MixCount != UINT32_MAX || attempt != UINT32_MAX)
            {
                std::fprintf(stderr,
                    "solve: rejected packet selector wrote output K=%u "
                    "P=%u\n", test.K, test.P);
                return false;
            }

            wirehair_v2::PrecodeSystem selected_system;
            selected_system.Params.BlockCount = UINT32_MAX;
            if (wirehair_v2::SelectSystematicConfiguration(
                    params, config, selected_system, selected, &attempt) !=
                    Wirehair_InvalidInput ||
                selected_system.Params.BlockCount != UINT32_MAX ||
                selected.PeelSeed != UINT32_C(0xdecafbad) ||
                selected.MixCount != UINT32_MAX || attempt != UINT32_MAX)
            {
                std::fprintf(stderr,
                    "solve: rejected joint selector wrote output K=%u P=%u\n",
                    test.K, test.P);
                return false;
            }
        }
    }

    // Exhaust the 16-bit public-id sample domain at a tiny span, then pin
    // high-bit IDs separately above.  Every accepted row must retain exactly
    // two distinct in-range precode columns.
    wirehair_v2::PacketRowConfig small_config;
    small_config.PeelSeed = UINT32_C(0x6f726163);
    small_config.MixCount = 2u;
    for (uint32_t id = 0; id <= UINT16_MAX; ++id)
    {
        const std::vector<uint32_t> row =
            wirehair_v2::GeneratePacketMatrixRow(
                2u, 2u, id, small_config);
        if (!PacketRowHasDistinctMix(row, 2u, 2u, 2u)) {
            std::fprintf(stderr,
                "solve: tiny packet-id sweep failed id=%u\n", id);
            return false;
        }
    }

    // Fixed-seed property samples cover the full accepted P range and both
    // low/high public block-id bits.  This complements the exhaustive tiny
    // sweep without turning the ordinary test lane into a multi-billion-row
    // job over the complete uint32 id space.
    uint64_t random_state = UINT64_C(0x8b79d42f6a1ce503);
    const auto next_random = [&random_state]() -> uint32_t {
        random_state ^= random_state >> 12;
        random_state ^= random_state << 25;
        random_state ^= random_state >> 27;
        return (uint32_t)(
            random_state * UINT64_C(0x2545f4914f6cdd1d) >> 32);
    };
    for (uint32_t trial = 0; trial < 100000u; ++trial)
    {
        const uint32_t P = 2u +
            next_random() % (wirehair_v2::kMaxPacketPrecodeCount - 1u);
        const uint32_t max_k =
            std::min<uint32_t>(64000u, UINT16_MAX - P);
        const uint32_t K = 2u + next_random() % (max_k - 1u);
        wirehair_v2::PacketRowConfig config;
        config.PeelSeed = next_random();
        config.MixCount = 1u + next_random() %
            std::min<uint32_t>(wirehair_v2::kCertifiedPacketMixCount, P);
        const uint32_t id = next_random();
        const std::vector<uint32_t> row =
            wirehair_v2::GeneratePacketMatrixRow(K, P, id, config);
        if (!PacketRowHasDistinctMix(row, K, P, config.MixCount))
        {
            std::fprintf(stderr,
                "solve: fixed property failed trial=%u K=%u P=%u id=%u\n",
                trial, K, P, id);
            return false;
        }
    }

    // The certified profile rule stays wholly inside the narrower packet
    // domain for every supported source count, so message facades cannot
    // select a structurally valid but unevaluable packet profile.
    for (uint32_t K = 2u; K <= 64000u; ++K)
    {
        const wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(K, 0u);
        const uint32_t P = params.Staircase +
            params.DenseRows + params.HeavyRows;
        if (!wirehair_v2::IsPacketRowDomainValid(
                K, P, wirehair_v2::kCertifiedPacketMixCount))
        {
            std::fprintf(stderr,
                "solve: certified packet profile escaped domain K=%u P=%u\n",
                K, P);
            return false;
        }
    }

    std::printf("packet row domain boundaries: PASS\n");
    return true;
}

bool CheckLargePacketEvaluationWork()
{
    const uint32_t K = 64000u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x987654321));
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        return false;
    }
    const uint32_t L = K + params.Staircase +
        params.DenseRows + params.HeavyRows;
    std::vector<uint8_t> intermediate(L);
    for (uint32_t i = 0; i < L; ++i) {
        intermediate[i] = (uint8_t)(i * 17u + 3u);
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = 0x5eedu;
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;

    uint64_t total_work = 0u;
    uint64_t digest = 0u;
    uint8_t output = 0u;
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (uint32_t i = 0; i < 4096u; ++i)
    {
        uint64_t work = 0u;
        const uint32_t id = UINT32_MAX - i * 7919u;
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, intermediate.data(), 1u, id, &output, &work) ||
            work < 4u || work > 67u)
        {
            std::fprintf(stderr,
                "solve: K=64000 packet work invalid id=%u work=%llu\n",
                id, (unsigned long long)work);
            return false;
        }
        total_work += work;
        digest = digest * UINT64_C(0x9e3779b97f4a7c15) + output + work;
    }
    const double milliseconds =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin).count();
    if (total_work > UINT64_C(4096) * 67u || digest == 0u) {
        return false;
    }
    std::printf(
        "K=64000 packet evaluation: rows=4096 work=%llu time=%.3fms "
        "digest=%llu: PASS\n",
        (unsigned long long)total_work,
        milliseconds,
        (unsigned long long)digest);
    return true;
}

bool CheckInactiveResidualCap()
{
    const uint32_t K = 5000u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x524553434150)),
            system))
    {
        return false;
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = 17u;
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    const uint8_t zero = 0u;
    std::vector<wirehair_v2::SolvePacket> packets(K);
    for (wirehair_v2::SolvePacket& packet : packets) {
        packet.BlockId = 0u;
        packet.Data = &zero;
    }
    std::vector<uint8_t> output(7u, 0xccu);
    const std::vector<uint8_t> before = output;
    wirehair_v2::PrecodeSolveStats stats;
    const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
        system, config, packets, 1u, output, &stats);
    if (result != Wirehair_NeedMore || output != before ||
        stats.InactivatedColumns <= wirehair_v2::kMaxInactiveColumns ||
        stats.PeelNanoseconds == 0u)
    {
        std::fprintf(stderr,
            "solve: inactive cap failed result=%d inact=%u\n",
            (int)result, stats.InactivatedColumns);
        return false;
    }
    return true;
}

uint64_t RecoveryMix64(uint64_t x)
{
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

bool RunMixedRecoveryBenchmark(uint32_t trials)
{
    const uint32_t Ks[] = {1000u, 10000u};
    const uint32_t overheads[] = {0u, 1u, 2u};
    for (uint32_t K : Ks)
    {
        for (uint32_t field = 0; field < 2u; ++field)
        {
            const bool mixed = field != 0u;
            const uint32_t block_bytes = 2u;
            std::vector<uint8_t> message((size_t)K * block_bytes, 0u);
            wirehair_v2::MessagePrecodeEncoderOptions options;
            if (mixed) {
                options.Completion =
                    wirehair_v2::CompletionField::MixedGF256GF16;
            }
            wirehair_v2::MessagePrecodeEncoder encoder;
            if (encoder.InitializeResult(
                    message.data(), message.size(), block_bytes,
                    nullptr, &options) != Wirehair_Success)
            {
                return false;
            }
            wirehair_v2::MessagePrecodeDecoder decoder;
            if (decoder.InitializeResult(
                    message.size(), block_bytes, &encoder.Profile()) !=
                    Wirehair_Success)
            {
                return false;
            }
            wirehair_v2::PacketRowConfig config;
            config.PeelSeed = decoder.PacketPeelSeed();
            config.MixCount = options.RecoveryMixCount;
            const wirehair_v2::PrecodeSystem& system = decoder.System();
            for (uint32_t mode = 0; mode < 2u; ++mode)
            {
                for (uint32_t overhead : overheads)
                {
                    std::atomic<uint32_t> next_trial(0u);
                    std::atomic<uint32_t> failures(0u);
                    std::atomic<uint32_t> errors(0u);
                    const uint32_t hardware = std::max(
                        1u, std::thread::hardware_concurrency());
                    const uint32_t worker_count = std::min(trials, hardware);
                    std::vector<std::thread> workers;
                    workers.reserve(worker_count);
                    const std::chrono::steady_clock::time_point begin =
                        std::chrono::steady_clock::now();
                    for (uint32_t worker = 0;
                         worker < worker_count; ++worker)
                    {
                        workers.push_back(std::thread([&]() {
                            const uint8_t zero[2] = {0u, 0u};
                            std::vector<wirehair_v2::SolvePacket> packets;
                            std::vector<uint8_t> output;
                            for (;;)
                            {
                                const uint32_t trial = next_trial.fetch_add(1u);
                                if (trial >= trials) break;
                                packets.clear();
                                packets.reserve((size_t)K + overhead);
                                if (mode == 0u)
                                {
                                    uint32_t id = 0u;
                                    while (packets.size() <
                                           (size_t)K + overhead)
                                    {
                                        const uint64_t random = RecoveryMix64(
                                            ((uint64_t)trial << 32) ^ id);
                                        if (random % 10u != 0u) {
                                            wirehair_v2::SolvePacket packet;
                                            packet.BlockId = id;
                                            packet.Data = zero;
                                            packets.push_back(packet);
                                        }
                                        ++id;
                                    }
                                }
                                else
                                {
                                    const uint32_t start = K +
                                        (uint32_t)RecoveryMix64(trial) * 2u;
                                    for (uint32_t i = 0;
                                         i < K + overhead; ++i)
                                    {
                                        wirehair_v2::SolvePacket packet;
                                        packet.BlockId = start + i;
                                        packet.Data = zero;
                                        packets.push_back(packet);
                                    }
                                }
                                output.clear();
                                const WirehairResult result =
                                    wirehair_v2::SolvePrecodeSystem(
                                        system, config, packets, block_bytes,
                                        output);
                                if (result == Wirehair_NeedMore) {
                                    failures.fetch_add(1u);
                                }
                                else if (result != Wirehair_Success) {
                                    errors.fetch_add(1u);
                                }
                            }
                        }));
                    }
                    for (std::thread& worker : workers) worker.join();
                    const double seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - begin).count();
                    if (errors.load() != 0u) return false;
                    std::printf(
                        "mixed_recovery_bench,profile=%s,K=%u,mode=%s,"
                        "overhead=%u,trials=%u,failures=%u,attempt=%u,"
                        "seconds=%.3f\n",
                        mixed ? "mixed" : "certified", K,
                        mode == 0u ? "loss10" : "repair-only",
                        overhead, trials, failures.load(),
                        encoder.Profile().V2SeedAttempt, seconds);
                }
            }
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    static_assert(
        wirehair_v2::kPacketRowContractVersion == 4u,
        "shipping packet-row contract must be version 4");
    if (argc == 2 &&
        std::strcmp(argv[1], "--mixed-projection-oracle") == 0)
    {
        return CheckMixedProjectionResidueBucketsOracle() ? 0 : 1;
    }
    if (argc == 2 &&
        std::strcmp(argv[1], "--mixed-rhs-fusion-oracle") == 0)
    {
        return CheckMixedRhsFusionOracle() ? 0 : 1;
    }
    if (argc == 2 &&
        std::strcmp(argv[1], "--solve-arena-oracle") == 0)
    {
        return CheckIncrementalResume() && CheckMixedSystematicSolve() ?
            0 : 1;
    }
    if (argc == 3 &&
        std::strcmp(argv[1], "--recovery-benchmark") == 0)
    {
        const uint32_t trials =
            (uint32_t)std::strtoul(argv[2], nullptr, 10);
        return trials == 0u || !RunMixedRecoveryBenchmark(trials) ? 1 : 0;
    }
    if (argc == 3)
    {
        const uint32_t K = (uint32_t)std::strtoul(argv[1], nullptr, 10);
        const uint32_t block_bytes =
            (uint32_t)std::strtoul(argv[2], nullptr, 10);
        if (K < 2u || K > 64000u || block_bytes == 0u) {
            return 2;
        }
        std::vector<uint8_t> message((size_t)K * block_bytes, 0x5au);
        const wirehair_v2::SeedProfile diagnostic_profile =
            wirehair_v2::SelectSeedProfile(K, block_bytes);
        wirehair_v2::PrecodeParams diagnostic_params =
            wirehair_v2::MakeCertifiedParams(
                K,
                wirehair_v2::MatrixSeedFromProfile(
                    diagnostic_profile,
                    0u,
                    wirehair_v2::kMessagePrecodeSeedSalt));
        wirehair_v2::PrecodeSystem diagnostic_system;
        if (!wirehair_v2::BuildPrecodeSystem(
                diagnostic_params, diagnostic_system))
        {
            return 1;
        }
        wirehair_v2::PacketRowConfig diagnostic_config;
        diagnostic_config.PeelSeed =
            wirehair_v2::PacketPeelSeedFromProfile(
                diagnostic_profile,
                wirehair_v2::kMessageRecoveryRowSeedSalt);
        diagnostic_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
        std::printf(
            "seeds profile_peel=%u profile_dense=%u dense_count=%u "
            "precode=0x%llx packet=0x%x\n",
            diagnostic_profile.PeelSeed,
            diagnostic_profile.DenseSeed,
            diagnostic_profile.DenseCount,
            (unsigned long long)diagnostic_params.Seed,
            diagnostic_config.PeelSeed);
        std::vector<wirehair_v2::SolvePacket> diagnostic_packets(K);
        for (uint32_t id = 0; id < K; ++id) {
            diagnostic_packets[id].BlockId = id;
            diagnostic_packets[id].Data =
                message.data() + (size_t)id * block_bytes;
        }
        uint32_t direct_binary_rank = 0u;
        if (K + diagnostic_params.Staircase + diagnostic_params.DenseRows +
                diagnostic_params.HeavyRows <= 64u)
        {
            std::vector<uint64_t> masks;
            const auto add_mask = [&](const std::vector<uint32_t>& row) {
                uint64_t mask = 0u;
                for (uint32_t column : row) {
                    mask ^= UINT64_C(1) << column;
                }
                masks.push_back(mask);
            };
            for (const std::vector<uint32_t>& row :
                    diagnostic_system.StaircaseRows) {
                add_mask(row);
            }
            for (const std::vector<uint32_t>& row :
                    diagnostic_system.DenseRowColumns) {
                add_mask(row);
            }
            for (uint32_t id = 0; id < K; ++id) {
                add_mask(wirehair_v2::GeneratePacketMatrixRow(
                    K,
                    diagnostic_params.Staircase +
                        diagnostic_params.DenseRows +
                        diagnostic_params.HeavyRows,
                    id,
                    diagnostic_config));
            }
            for (uint32_t column = 0; column < 64u; ++column)
            {
                uint32_t pivot = direct_binary_rank;
                while (pivot < masks.size() &&
                       (masks[pivot] & (UINT64_C(1) << column)) == 0u) {
                    ++pivot;
                }
                if (pivot == masks.size()) {
                    continue;
                }
                std::swap(masks[pivot], masks[direct_binary_rank]);
                for (uint32_t r = 0; r < masks.size(); ++r) {
                    if (r != direct_binary_rank &&
                        (masks[r] & (UINT64_C(1) << column)) != 0u) {
                        masks[r] ^= masks[direct_binary_rank];
                    }
                }
                ++direct_binary_rank;
            }
        }
        std::vector<uint8_t> diagnostic_intermediate;
        wirehair_v2::PrecodeSolveStats diagnostic_stats;
        const WirehairResult diagnostic_result =
            wirehair_v2::SolvePrecodeSystem(
                diagnostic_system,
                diagnostic_config,
                diagnostic_packets,
                block_bytes,
                diagnostic_intermediate,
                &diagnostic_stats);
        std::printf(
            "base result=%d peeled=%u inact=%u residual_rows=%u "
            "binary_rank=%u direct_binary=%u rank=%u "
            "row_refs=%llu row_alloc=%u adjacency_alloc=%u "
            "binary_storage_bytes=%llu\n",
            (int)diagnostic_result,
            diagnostic_stats.PeeledColumns,
            diagnostic_stats.InactivatedColumns,
            diagnostic_stats.ResidualRows,
            diagnostic_stats.BinaryResidualRank,
            direct_binary_rank,
            diagnostic_stats.ResidualRank,
            (unsigned long long)diagnostic_stats.BinaryRowReferences,
            diagnostic_stats.BinaryRowStorageAllocations,
            diagnostic_stats.BinaryAdjacencyStorageAllocations,
            (unsigned long long)(diagnostic_stats.BinaryRowStorageBytes +
                diagnostic_stats.BinaryAdjacencyStorageBytes));
        wirehair_v2::MessagePrecodeEncoder encoder;
        const WirehairResult result = encoder.InitializeResult(
            message.data(), message.size(), block_bytes);
        const wirehair_v2::PrecodeSolveStats& stats = encoder.SolveStats();
        std::printf(
            "profile K=%u bb=%u result=%d attempt=%u inact=%u rank=%u "
            "build=%.3fms peel=%.3fms project=%.3fms residual=%.3fms "
            "backsub=%.3fms\n",
            K, block_bytes, (int)result,
            stats.PacketSeedAttempt,
            stats.InactivatedColumns, stats.ResidualRank,
            stats.BuildNanoseconds / 1000000.0,
            stats.PeelNanoseconds / 1000000.0,
            stats.ProjectNanoseconds / 1000000.0,
            stats.ResidualNanoseconds / 1000000.0,
            stats.BackSubNanoseconds / 1000000.0);
        if (result == Wirehair_Success)
        {
            wirehair_v2::MessagePrecodeDecoder decoder;
            if (decoder.InitializeResult(
                    message.size(), block_bytes, &encoder.Profile()) !=
                    Wirehair_Success)
            {
                return 1;
            }
            std::vector<uint8_t> block(block_bytes);
            WirehairResult decode_result = Wirehair_NeedMore;
            uint32_t delivered = 0u;
            for (uint32_t id = 0u;
                 decode_result == Wirehair_NeedMore && id < K * 2u;
                 ++id)
            {
                if (id % 10u == 0u) {
                    continue;
                }
                uint32_t bytes = 0u;
                if (encoder.EncodeResult(
                        id, block.data(), block_bytes, &bytes) !=
                        Wirehair_Success)
                {
                    return 1;
                }
                decode_result = decoder.DecodeResult(
                    id, block.data(), bytes);
                ++delivered;
            }
            const wirehair_v2::PrecodeSolveStats& decode_stats =
                decoder.SolveStats();
            std::printf(
                "decode result=%d delivered=%u inact=%u rank=%u "
                "build=%.3fms peel=%.3fms project=%.3fms residual=%.3fms "
                "backsub=%.3fms\n",
                (int)decode_result, delivered,
                decode_stats.InactivatedColumns,
                decode_stats.ResidualRank,
                decode_stats.BuildNanoseconds / 1000000.0,
                decode_stats.PeelNanoseconds / 1000000.0,
                decode_stats.ProjectNanoseconds / 1000000.0,
                decode_stats.ResidualNanoseconds / 1000000.0,
                decode_stats.BackSubNanoseconds / 1000000.0);
            if (decode_result != Wirehair_Success) {
                return 1;
            }
        }
        return result == Wirehair_Success ? 0 : 1;
    }
    if (argc != 1) {
        return 2;
    }
    bool ok = true;
    ok = CheckLowestBitIndex() && ok;
    ok = CheckPacketEvaluationFusion() && ok;
    ok = CheckOddPacketPeelSeedInterleave() && ok;
    ok = CheckPacketRowSeedPermutation() && ok;
    ok = CheckPacketRuntimeBoundaries() && ok;
    ok = CheckTinyDenseOracle() && ok;
    ok = CheckHeavyCoefficientBoundaryOracle() && ok;
    ok = CheckBinaryQuotientBoundary() && ok;
    ok = CheckConcurrentCoefficientCaches() && ok;
    ok = CheckIncrementalResume() && ok;
    ok = CheckMixedProjectionResidueBucketsOracle() && ok;
    ok = CheckMixedMix1EndToEnd() && ok;
    ok = CheckBinaryPeelLowDegreeXorOracle() && ok;
    ok = CheckMixedSystematicSolve() && ok;
    ok = CheckPackedBinaryResidualOracle() && ok;
    ok = CheckMixedRhsFusionOracle() && ok;
    ok = CheckMixedNullWitnessCanonicalization() && ok;
    ok = CheckMixedQuotientRankFirstOracles() && ok;
    ok = CheckMixDomainValidation() && ok;
    ok = CheckPacketRowDomainBoundaries() && ok;
    ok = CheckInactiveResidualCap() && ok;
    ok = CheckLargePacketEvaluationWork() && ok;
    ok = RunCase(64u, 17u, 7u) && ok;
    ok = RunCase(
        64u, 17u, 7u,
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero) && ok;
    ok = RunCase(
        64u, wirehair_v2::kBinaryQuotientMinBlockBytes, 7u,
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero) && ok;
    ok = RunCase(320u, 37u, 11u) && ok;
    return ok ? 0 : 1;
}
