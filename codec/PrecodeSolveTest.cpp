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

class HeavyProjectionOracleScope
{
public:
    HeavyProjectionOracleScope()
    {
        wirehair_v2::SetHeavyProjectionOracleForTesting(true);
    }

    ~HeavyProjectionOracleScope()
    {
        wirehair_v2::SetHeavyProjectionOracleForTesting(false);
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
        a.PacketSeedAttempt == b.PacketSeedAttempt;
}

void ClearSolveTiming(wirehair_v2::PrecodeSolveStats& stats)
{
    stats.BuildNanoseconds = 0u;
    stats.PeelNanoseconds = 0u;
    stats.ProjectNanoseconds = 0u;
    stats.ResidualNanoseconds = 0u;
    stats.BackSubNanoseconds = 0u;
}

bool SameSolveStatsIgnoringTiming(
    wirehair_v2::PrecodeSolveStats a,
    wirehair_v2::PrecodeSolveStats b)
{
    ClearSolveTiming(a);
    ClearSolveTiming(b);
    return SameSolveStats(a, b);
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

bool SameResumeStateIgnoringTiming(
    wirehair_v2::PrecodeSolveResumeState a,
    wirehair_v2::PrecodeSolveResumeState b)
{
    ClearSolveTiming(a.Stats);
    ClearSolveTiming(b.Stats);
    return SameResumeState(a, b);
}

struct ResumeStorageIdentity
{
    const uint32_t* InactiveIndexData;
    size_t InactiveIndexCapacity;
    const uint32_t* InactiveColumnsData;
    size_t InactiveColumnsCapacity;
    const uint64_t* ProjectionData;
    size_t ProjectionCapacity;
    const uint8_t* ValuesData;
    size_t ValuesCapacity;
    const uint8_t* PivotCoefficientsData;
    size_t PivotCoefficientsCapacity;
    const uint8_t* PivotRhsData;
    size_t PivotRhsCapacity;
    const uint8_t* HavePivotData;
    size_t HavePivotCapacity;
    const uint8_t* CoefficientScratchData;
    size_t CoefficientScratchCapacity;
    const uint8_t* RhsScratchData;
    size_t RhsScratchCapacity;
};

ResumeStorageIdentity CaptureResumeStorageIdentity(
    const wirehair_v2::PrecodeSolveResumeState& state)
{
    ResumeStorageIdentity identity = {};
    identity.InactiveIndexData = state.InactiveIndex.data();
    identity.InactiveIndexCapacity = state.InactiveIndex.capacity();
    identity.InactiveColumnsData = state.InactiveColumns.data();
    identity.InactiveColumnsCapacity = state.InactiveColumns.capacity();
    identity.ProjectionData = state.Projection.data();
    identity.ProjectionCapacity = state.Projection.capacity();
    identity.ValuesData = state.Values.data();
    identity.ValuesCapacity = state.Values.capacity();
    identity.PivotCoefficientsData = state.PivotCoefficients.data();
    identity.PivotCoefficientsCapacity = state.PivotCoefficients.capacity();
    identity.PivotRhsData = state.PivotRhs.data();
    identity.PivotRhsCapacity = state.PivotRhs.capacity();
    identity.HavePivotData = state.HavePivot.data();
    identity.HavePivotCapacity = state.HavePivot.capacity();
    identity.CoefficientScratchData = state.CoefficientScratch.data();
    identity.CoefficientScratchCapacity =
        state.CoefficientScratch.capacity();
    identity.RhsScratchData = state.RhsScratch.data();
    identity.RhsScratchCapacity = state.RhsScratch.capacity();
    return identity;
}

bool SameResumeStorageIdentity(
    const wirehair_v2::PrecodeSolveResumeState& state,
    const ResumeStorageIdentity& identity)
{
    return state.InactiveIndex.data() == identity.InactiveIndexData &&
        state.InactiveIndex.capacity() == identity.InactiveIndexCapacity &&
        state.InactiveColumns.data() == identity.InactiveColumnsData &&
        state.InactiveColumns.capacity() ==
            identity.InactiveColumnsCapacity &&
        state.Projection.data() == identity.ProjectionData &&
        state.Projection.capacity() == identity.ProjectionCapacity &&
        state.Values.data() == identity.ValuesData &&
        state.Values.capacity() == identity.ValuesCapacity &&
        state.PivotCoefficients.data() == identity.PivotCoefficientsData &&
        state.PivotCoefficients.capacity() ==
            identity.PivotCoefficientsCapacity &&
        state.PivotRhs.data() == identity.PivotRhsData &&
        state.PivotRhs.capacity() == identity.PivotRhsCapacity &&
        state.HavePivot.data() == identity.HavePivotData &&
        state.HavePivot.capacity() == identity.HavePivotCapacity &&
        state.CoefficientScratch.data() ==
            identity.CoefficientScratchData &&
        state.CoefficientScratch.capacity() ==
            identity.CoefficientScratchCapacity &&
        state.RhsScratch.data() == identity.RhsScratchData &&
        state.RhsScratch.capacity() == identity.RhsScratchCapacity;
}

enum class ColdSolveEntryPoint
{
    Public,
    Runtime,
    ValidatedRuntime
};

const char* ColdSolveEntryPointName(ColdSolveEntryPoint entry_point)
{
    switch (entry_point)
    {
    case ColdSolveEntryPoint::Public:
        return "public";
    case ColdSolveEntryPoint::Runtime:
        return "runtime";
    case ColdSolveEntryPoint::ValidatedRuntime:
        return "validated-runtime";
    }
    return "unknown";
}

WirehairResult CallColdSolve(
    ColdSolveEntryPoint entry_point,
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    const std::vector<wirehair_v2::SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& output,
    wirehair_v2::PrecodeSolveStats* stats,
    wirehair_v2::PrecodeSolveResumeState* resume_state)
{
    switch (entry_point)
    {
    case ColdSolveEntryPoint::Public:
        return wirehair_v2::SolvePrecodeSystem(
            system, config, packets, block_bytes,
            output, stats, resume_state);
    case ColdSolveEntryPoint::Runtime:
        return wirehair_v2::SolvePrecodeSystemWithRuntime(
            system, config, runtime, packets, block_bytes,
            output, stats, resume_state);
    case ColdSolveEntryPoint::ValidatedRuntime:
        return wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
            system, config, runtime, packets, block_bytes,
            output, stats, resume_state);
    }
    return Wirehair_Error;
}

wirehair_v2::PrecodeSolveResumeState MakeStatsAliasSentinel(uint32_t tag)
{
    wirehair_v2::PrecodeSolveResumeState state;
    state.SourceCount = UINT32_C(0x10203040) ^ tag;
    state.PrecodeCount = UINT32_C(0x21314151) ^ tag;
    state.ColumnCount = UINT32_C(0x32425262) ^ tag;
    state.BlockBytes = UINT32_C(0x43536373) ^ tag;
    state.InactiveCount = UINT32_C(0x54647484) ^ tag;
    state.ProjectionWords = UINT32_C(0x65758595) ^ tag;
    state.Rank = UINT32_C(0x768696a6) ^ tag;
    state.Config.PeelSeed = UINT32_C(0x8797a7b7) ^ tag;
    state.Config.MixCount = 2u;
    state.Stats.PacketRows = UINT32_C(0x98a8b8c8) ^ tag;
    state.Stats.PeeledColumns = UINT32_C(0xa9b9c9d9) ^ tag;
    state.Stats.InactivatedColumns = UINT32_C(0xbacadbea) ^ tag;
    state.Stats.ResidualRows = UINT32_C(0xcbdcedfb) ^ tag;
    state.Stats.BlockXors = UINT64_C(0x123456789abcdef0) ^ tag;
    state.Stats.BlockMulAdds = UINT64_C(0xfedcba9876543210) ^ tag;
    state.Stats.PacketSeedAttempt = UINT32_C(0xdcedfe0f) ^ tag;
    state.InactiveIndex.assign(5u, UINT32_C(0x11111111) ^ tag);
    state.InactiveColumns.assign(6u, UINT32_C(0x22222222) ^ tag);
    state.Projection.assign(7u, UINT64_C(0x3333333333333333) ^ tag);
    state.Values.assign(8u, (uint8_t)(0x40u ^ tag));
    state.PivotCoefficients.assign(9u, (uint8_t)(0x50u ^ tag));
    state.PivotRhs.assign(10u, (uint8_t)(0x60u ^ tag));
    state.HavePivot.assign(11u, (uint8_t)(0x70u ^ tag));
    state.CoefficientScratch.assign(12u, (uint8_t)(0x80u ^ tag));
    state.RhsScratch.assign(13u, (uint8_t)(0x90u ^ tag));
    state.Active = (tag & 1u) != 0u;
    return state;
}

bool ExpectColdSolveStatsAliasRejected(
    const char* scenario,
    ColdSolveEntryPoint entry_point,
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    const std::vector<wirehair_v2::SolvePacket>& packets,
    uint32_t block_bytes,
    uint32_t tag)
{
    wirehair_v2::PrecodeSolveResumeState state =
        MakeStatsAliasSentinel(tag);
    const wirehair_v2::PrecodeSolveResumeState state_before = state;
    const ResumeStorageIdentity storage_before =
        CaptureResumeStorageIdentity(state);
    const size_t persistent_before = state.PersistentBytes();
    std::vector<uint8_t> output(15u, (uint8_t)(0xa0u ^ tag));
    const std::vector<uint8_t> output_before = output;
    const uint8_t* const output_data_before = output.data();
    const size_t output_capacity_before = output.capacity();

    const WirehairResult result = CallColdSolve(
        entry_point, system, config, runtime, packets, block_bytes,
        output, &state.Stats, &state);
    if (result != Wirehair_InvalidInput ||
        !SameResumeState(state, state_before) ||
        !SameResumeStorageIdentity(state, storage_before) ||
        state.PersistentBytes() != persistent_before ||
        output != output_before ||
        output.data() != output_data_before ||
        output.capacity() != output_capacity_before)
    {
        std::fprintf(stderr,
            "solve: cold stats/checkpoint alias changed state "
            "scenario=%s entry=%s result=%d\n",
            scenario, ColdSolveEntryPointName(entry_point), (int)result);
        return false;
    }
    return true;
}

bool ExpectResumeStatsAliasRejected(
    const char* scenario,
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_id,
    const uint8_t* block_data,
    uint32_t block_bytes,
    const wirehair_v2::PrecodeSolveResumeState& prototype,
    bool allow_insert)
{
    wirehair_v2::PrecodeSolveResumeState state = prototype;
    const wirehair_v2::PrecodeSolveResumeState state_before = state;
    const ResumeStorageIdentity storage_before =
        CaptureResumeStorageIdentity(state);
    const size_t persistent_before = state.PersistentBytes();
    std::vector<uint8_t> output(16u, 0xb3u);
    const std::vector<uint8_t> output_before = output;
    const uint8_t* const output_data_before = output.data();
    const size_t output_capacity_before = output.capacity();

    const WirehairResult result = wirehair_v2::ResumePrecodeSystem(
        system, config, block_id, block_data, block_bytes,
        state, output, &state.Stats, allow_insert);
    if (result != Wirehair_InvalidInput ||
        !SameResumeState(state, state_before) ||
        !SameResumeStorageIdentity(state, storage_before) ||
        state.PersistentBytes() != persistent_before ||
        output != output_before ||
        output.data() != output_data_before ||
        output.capacity() != output_capacity_before)
    {
        std::fprintf(stderr,
            "solve: resume stats/checkpoint alias changed state "
            "scenario=%s allow_insert=%u result=%d\n",
            scenario, allow_insert ? 1u : 0u, (int)result);
        return false;
    }
    return true;
}

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

enum PacketEvaluatorKind
{
    PacketEvaluatorGeneric,
    PacketEvaluatorValidated,
    PacketEvaluatorRuntime
};

bool CallPacketEvaluator(
    PacketEvaluatorKind kind,
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    const uint8_t* intermediate,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* output,
    uint64_t* operations)
{
    switch (kind)
    {
    case PacketEvaluatorGeneric:
        return wirehair_v2::EvaluatePacketBlock(
            system, config, intermediate, block_bytes, block_id,
            output, operations);
    case PacketEvaluatorValidated:
        return wirehair_v2::EvaluatePacketBlockForValidatedSystem(
            system, config, intermediate, block_bytes, block_id,
            output, operations);
    case PacketEvaluatorRuntime:
        return wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
            system, config, runtime, intermediate, block_bytes, block_id,
            output, operations);
    }
    return true;
}

bool CheckPacketOperationAliases(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_id,
    uint32_t block_bytes,
    bool all_entry_points)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const size_t intermediate_bytes = (size_t)(K + P) * block_bytes;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }
    const unsigned entry_point_count = all_entry_points ? 3u : 1u;
    for (unsigned entry_point = 0u;
         entry_point < entry_point_count;
         ++entry_point)
    {
        const PacketEvaluatorKind kind =
            (PacketEvaluatorKind)entry_point;

        // Exact counter/output overlap.  Both views begin at an aligned,
        // live uint64_t object so the rejected call itself is sanitizer-safe.
        std::vector<uint64_t> input_words(
            (intermediate_bytes + 7u) / 8u, UINT64_C(0x4b4b4b4b4b4b4b4b));
        std::vector<uint64_t> output_words(
            (block_bytes + 15u) / 8u, UINT64_C(0xa5a5a5a5a5a5a5a5));
        uint8_t* const input =
            reinterpret_cast<uint8_t*>(input_words.data());
        uint8_t* const output =
            reinterpret_cast<uint8_t*>(output_words.data());
        for (size_t i = 0u; i < intermediate_bytes; ++i) {
            input[i] = (uint8_t)(i * 29u + block_id * 7u + 3u);
        }
        const std::vector<uint64_t> input_before = input_words;
        const std::vector<uint64_t> output_before = output_words;
        if (CallPacketEvaluator(
                kind, system, config, runtime, input, block_bytes, block_id,
                output, output_words.data()) ||
            input_words != input_before || output_words != output_before)
        {
            std::fprintf(stderr,
                "solve: exact packet counter/output alias wrote kind=%u "
                "mix=%u bb=%u\n",
                entry_point, config.MixCount, block_bytes);
            return false;
        }

        // Move the byte output one byte past the aligned allocation start.
        // The last aligned counter object before its end then crosses the
        // output boundary by one to seven bytes.
        output_words = std::vector<uint64_t>(
            (block_bytes + 23u) / 8u, UINT64_C(0x9696969696969696));
        uint8_t* const partial_output =
            reinterpret_cast<uint8_t*>(output_words.data()) + 1u;
        const size_t partial_output_end = 1u + block_bytes;
        uint64_t* const partial_output_ops = output_words.data() +
            ((partial_output_end - 1u) / 8u);
        const std::vector<uint64_t> partial_output_before = output_words;
        if (CallPacketEvaluator(
                kind, system, config, runtime, input, block_bytes, block_id,
                partial_output, partial_output_ops) ||
            input_words != input_before ||
            output_words != partial_output_before)
        {
            std::fprintf(stderr,
                "solve: partial packet counter/output alias wrote kind=%u "
                "mix=%u bb=%u\n",
                entry_point, config.MixCount, block_bytes);
            return false;
        }

        output_words = std::vector<uint64_t>(
            (block_bytes + 15u) / 8u, UINT64_C(0x8787878787878787));
        uint8_t* const disjoint_output =
            reinterpret_cast<uint8_t*>(output_words.data());
        const std::vector<uint64_t> disjoint_output_before = output_words;

        // Exact counter/intermediate overlap.
        if (CallPacketEvaluator(
                kind, system, config, runtime, input, block_bytes, block_id,
                disjoint_output, input_words.data()) ||
            input_words != input_before ||
            output_words != disjoint_output_before)
        {
            std::fprintf(stderr,
                "solve: exact packet counter/input alias wrote kind=%u "
                "mix=%u bb=%u\n",
                entry_point, config.MixCount, block_bytes);
            return false;
        }

        // As above, offset the byte input by one so an aligned counter object
        // straddles its final byte without ever forming a misaligned uint64_t*.
        std::vector<uint64_t> partial_input_words(
            (intermediate_bytes + 16u) / 8u,
            UINT64_C(0x7878787878787878));
        uint8_t* const partial_input =
            reinterpret_cast<uint8_t*>(partial_input_words.data()) + 1u;
        for (size_t i = 0u; i < intermediate_bytes; ++i) {
            partial_input[i] = (uint8_t)(i * 31u + block_id * 5u + 9u);
        }
        const size_t partial_input_end = 1u + intermediate_bytes;
        uint64_t* const partial_input_ops = partial_input_words.data() +
            ((partial_input_end - 1u) / 8u);
        const std::vector<uint64_t> partial_input_before = partial_input_words;
        if (CallPacketEvaluator(
                kind, system, config, runtime, partial_input, block_bytes,
                block_id, disjoint_output, partial_input_ops) ||
            partial_input_words != partial_input_before ||
            output_words != disjoint_output_before)
        {
            std::fprintf(stderr,
                "solve: partial packet counter/input alias wrote kind=%u "
                "mix=%u bb=%u\n",
                entry_point, config.MixCount, block_bytes);
            return false;
        }
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

    // The normal schedule uses fewer than six terms; the paired schedule uses
    // the same small-block domain above the crossover.  Test all three entry
    // points and all supported mix counts because they share one preflight.
    for (uint32_t mix_count : kFusedMixCounts)
    {
        wirehair_v2::PacketRowConfig fused = config;
        fused.MixCount = mix_count;
        if (!CheckPacketOperationAliases(
                system, fused, ids[0], 17u, true) ||
            !CheckPacketOperationAliases(
                system, fused, ids[5], 17u, true))
        {
            return false;
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
        if (!CheckPacketOperationAliases(
                set_xor_system, set_xor_config, set_xor_id,
                set_xor_block_bytes, false))
        {
            std::fprintf(stderr,
                "solve: set-XOR packet operation alias guard failed\n");
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
    // The equation-preserving dense basis peels the eleven two-column deltas
    // directly, reducing the tiny residual completion work from 528.
    static const uint64_t kExpectedExecutedMulAdds = 492u;
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool CheckDenseAnchorPayloadOracle()
{
    static const wirehair_v2::DenseAnchorLayout kLayouts[] = {
        wirehair_v2::DenseAnchorLayout::Two07,
        wirehair_v2::DenseAnchorLayout::Four0369
    };
    static const uint32_t kBlockCounts[] = { 2u, 17u, 64u };
    static const uint32_t kBlockBytes[] = { 1u, 13u };

    for (wirehair_v2::DenseAnchorLayout layout : kLayouts)
    {
        for (uint32_t K : kBlockCounts)
        {
            wirehair_v2::PrecodeParams params =
                wirehair_v2::MakeCertifiedParams(
                    K,
                    UINT64_C(0x616e63686f72736f) ^
                        ((uint64_t)(uint32_t)layout << 32) ^ K);
            params.DenseAnchors = layout;
            wirehair_v2::PacketRowConfig base_config;
            base_config.PeelSeed =
                UINT32_C(0x9e3779b9) ^ K ^ (uint32_t)layout;
            base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;

            wirehair_v2::PrecodeSystem system;
            wirehair_v2::PacketRowConfig config;
            uint32_t attempt = 0u;
            if (wirehair_v2::SelectSystematicConfiguration(
                    params, base_config, system, config, &attempt) !=
                    Wirehair_Success ||
                system.Params.DenseAnchors != layout)
            {
                std::fprintf(stderr,
                    "solve: dense-anchor configuration failed "
                    "layout=%u K=%u\n",
                    (unsigned)layout, K);
                return false;
            }

            for (uint32_t block_bytes : kBlockBytes)
            {
                std::vector<uint8_t> message((size_t)K * block_bytes);
                for (size_t i = 0; i < message.size(); ++i) {
                    message[i] = (uint8_t)(
                        i * 43u + K * 11u + block_bytes * 7u +
                        (uint32_t)layout);
                }
                std::vector<wirehair_v2::SolvePacket> packets(K);
                for (uint32_t id = 0; id < K; ++id) {
                    packets[id].BlockId = id;
                    packets[id].Data =
                        message.data() + (size_t)id * block_bytes;
                }

                std::vector<uint8_t> solved;
                const WirehairResult solve_result =
                    wirehair_v2::SolvePrecodeSystem(
                        system, config, packets, block_bytes, solved);
                std::vector<uint8_t> oracle;
                const WirehairResult oracle_result =
                    wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
                        system, config, packets, block_bytes, oracle);
                if (solve_result != Wirehair_Success ||
                    oracle_result != Wirehair_Success || solved != oracle)
                {
                    std::fprintf(stderr,
                        "solve: dense-anchor payload oracle mismatch "
                        "layout=%u K=%u attempt=%u bb=%u solve=%d "
                        "oracle=%d\n",
                        (unsigned)layout, K, attempt, block_bytes,
                        (int)solve_result, (int)oracle_result);
                    return false;
                }

                std::vector<uint8_t> recovered(block_bytes);
                for (uint32_t id = 0; id < K; ++id)
                {
                    if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                            system,
                            config,
                            solved.data(),
                            block_bytes,
                            id,
                            recovered.data()) ||
                        std::memcmp(
                            recovered.data(),
                            message.data() + (size_t)id * block_bytes,
                            block_bytes) != 0)
                    {
                        std::fprintf(stderr,
                            "solve: dense-anchor payload replay mismatch "
                            "layout=%u K=%u attempt=%u bb=%u id=%u\n",
                            (unsigned)layout, K, attempt, block_bytes, id);
                        return false;
                    }
                }
            }
        }
    }
    std::printf("dense-anchor payload/oracle solve: PASS\n");
    return true;
}
#endif

bool CheckZeroRhsAcrossBlockWidths(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    WirehairResult expected)
{
    static const uint32_t kBlockBytes[] = { 1u, 64u, 2048u, 4096u };
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }
    for (uint32_t block_bytes : kBlockBytes)
    {
        std::vector<uint8_t> message((size_t)K * block_bytes, 0u);
        std::vector<wirehair_v2::SolvePacket> packets(K);
        for (uint32_t id = 0u; id < K; ++id)
        {
            packets[id].BlockId = id;
            packets[id].Data =
                message.data() + (size_t)id * block_bytes;
        }
        std::vector<uint8_t> intermediate;
        const WirehairResult actual =
            wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
                system,
                config,
                runtime,
                packets,
                block_bytes,
                intermediate);
        if (actual != expected)
        {
            std::fprintf(stderr,
                "solve: zero-RHS rank changed at bb=%u expected=%d actual=%d\n",
                block_bytes, (int)expected, (int)actual);
            return false;
        }
    }
    return true;
}

bool CheckExactSystematicFailureClassification()
{
    const uint32_t K = 8u;
    const wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(K, 64u);
    wirehair_v2::MessagePrecodeEncoderOptions options;
    wirehair_v2::PrecodeParams params;
    wirehair_v2::PacketRowConfig base_config;
    if (!wirehair_v2::ResolveMessagePrecodeOptions(
            profile, nullptr, options) ||
        !wirehair_v2::ResolveMessagePrecodeConfiguration(
            profile, options, params, base_config))
    {
        std::fprintf(stderr,
            "solve: exact-failure profile setup failed\n");
        return false;
    }

    wirehair_v2::PrecodeSystem full_rank_system;
    wirehair_v2::PacketRowConfig full_rank_config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params,
            base_config,
            full_rank_system,
            full_rank_config) != Wirehair_Success ||
        !CheckZeroRhsAcrossBlockWidths(
            full_rank_system, full_rank_config, Wirehair_Success))
    {
        std::fprintf(stderr,
            "solve: exact-failure full-rank setup failed\n");
        return false;
    }

    wirehair_v2::PrecodeSystem deficient_system;
    const wirehair_v2::PrecodeParams deficient_params =
        wirehair_v2::PrecodeParamsForAttempt(params, 2u);
    const wirehair_v2::PacketRowConfig deficient_config =
        wirehair_v2::PacketConfigForAttempt(base_config, 2u);
    if (!wirehair_v2::BuildPrecodeSystem(
            deficient_params, deficient_system) ||
        !CheckZeroRhsAcrossBlockWidths(
            deficient_system, deficient_config, Wirehair_NeedMore))
    {
        std::fprintf(stderr,
            "solve: exact-failure deficient setup failed\n");
        return false;
    }
    const uint32_t full_rank_P = full_rank_system.Params.Staircase +
        full_rank_system.Params.DenseRows + full_rank_system.Params.HeavyRows;
    const uint32_t deficient_P = deficient_system.Params.Staircase +
        deficient_system.Params.DenseRows + deficient_system.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime full_rank_runtime;
    wirehair_v2::PacketRowRuntime deficient_runtime;
    if (!full_rank_runtime.Initialize(
            K, full_rank_P, full_rank_config.MixCount) ||
        !deficient_runtime.Initialize(
            K, deficient_P, deficient_config.MixCount))
    {
        std::fprintf(stderr,
            "solve: exact-failure runtime setup failed\n");
        return false;
    }
    if (wirehair_v2::ClassifyExactSystematicConstructionFailure(
            full_rank_system,
            full_rank_config,
            full_rank_runtime,
            Wirehair_Error) != Wirehair_Error ||
        wirehair_v2::ClassifyExactSystematicConstructionFailure(
            deficient_system,
            deficient_config,
            deficient_runtime,
            Wirehair_Error) != Wirehair_BadPeelSeed ||
        wirehair_v2::ClassifyExactSystematicConstructionFailure(
            deficient_system,
            deficient_config,
            deficient_runtime,
            Wirehair_NeedMore) !=
                Wirehair_BadPeelSeed ||
        wirehair_v2::ClassifyExactSystematicConstructionFailure(
            full_rank_system,
            full_rank_config,
            full_rank_runtime,
            Wirehair_InvalidInput) !=
                Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "solve: exact-failure classification changed fatal status\n");
        return false;
    }
    std::printf("exact systematic failure classification: PASS\n");
    return true;
}

bool CheckHeavyCoefficientBoundaryOracle()
{
    wirehair_v2::ResetHeavyProjectionOracleCountersForTesting();
    struct Geometry
    {
        uint32_t DenseRows;
        uint32_t HeavyRows;
    };
    static const uint32_t kBoundaryHeavyRows[] = {
        0u, 1u, 11u, 12u, 13u, 128u
    };
    static const uint32_t kCandidateBlockCounts[] = {3u, 4u};
    static const Geometry kCandidateGeometries[] = {
        {12u, 11u},
        {12u, 13u},
        {13u, 12u}
    };
    static const wirehair_v2::HeavyCoefficientFamily kFamilies[] = {
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy,
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero
    };
    const uint32_t block_bytes = 5u;

    const auto check_case = [](
        uint32_t K,
        uint32_t dense_rows,
        uint32_t heavy_rows,
        wirehair_v2::HeavyCoefficientFamily family) -> bool
    {
        std::vector<uint8_t> message((size_t)K * block_bytes);
        for (size_t i = 0; i < message.size(); ++i) {
            message[i] = (uint8_t)(
                i * 97u + (i >> 2) + K * 17u + dense_rows * 7u +
                heavy_rows * 3u + (uint32_t)family);
        }
        std::vector<wirehair_v2::SolvePacket> packets(K);
        for (uint32_t id = 0u; id < K; ++id) {
            packets[id].BlockId = id;
            packets[id].Data = message.data() + (size_t)id * block_bytes;
        }

        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(
                K,
                UINT64_C(0x4845415659424f55) ^
                    ((uint64_t)K << 32) ^
                    ((uint64_t)dense_rows << 16) ^
                    ((uint64_t)heavy_rows << 8) ^ (uint32_t)family);
        params.DenseRows = dense_rows;
        params.HeavyRows = heavy_rows;
        params.HeavyFamily = family;
        wirehair_v2::PacketRowConfig base_config;
        base_config.PeelSeed =
            UINT32_C(0xa17e31d9) ^ K ^ (dense_rows << 16) ^
            heavy_rows ^ (uint32_t)family;
        base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
        wirehair_v2::PrecodeSystem system;
        wirehair_v2::PacketRowConfig config;
        if (wirehair_v2::SelectSystematicConfiguration(
                params, base_config, system, config) != Wirehair_Success ||
            system.Params.DenseRows != dense_rows ||
            system.Params.HeavyRows != heavy_rows ||
            system.Params.HeavyFamily != family)
        {
            std::fprintf(stderr,
                "solve: candidate configuration failed K=%u D=%u H=%u "
                "family=%u\n",
                K, dense_rows, heavy_rows, (unsigned)family);
            return false;
        }

        std::vector<uint8_t> solved;
        std::vector<uint8_t> oracle;
        WirehairResult solve_result = Wirehair_Error;
        {
            HeavyProjectionOracleScope oracle_scope;
            solve_result = wirehair_v2::SolvePrecodeSystem(
                system, config, packets, block_bytes, solved);
        }
        if (solve_result != Wirehair_Success ||
            wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
                system, config, packets, block_bytes, oracle) !=
                Wirehair_Success ||
            solved != oracle ||
            !wirehair_v2::VerifyPrecodeSolution(
                system, config, packets, solved.data(), block_bytes))
        {
            std::fprintf(stderr,
                "solve: dense oracle mismatch K=%u D=%u H=%u family=%u\n",
                K, dense_rows, heavy_rows, (unsigned)family);
            return false;
        }
        return true;
    };

    for (wirehair_v2::HeavyCoefficientFamily family : kFamilies) {
        for (uint32_t heavy_rows : kBoundaryHeavyRows) {
            if (!check_case(8u, 12u, heavy_rows, family)) {
                return false;
            }
        }
    }
    for (wirehair_v2::HeavyCoefficientFamily family : kFamilies)
    {
        for (uint32_t K : kCandidateBlockCounts) {
            for (const Geometry& geometry : kCandidateGeometries) {
                if (!check_case(
                        K, geometry.DenseRows, geometry.HeavyRows, family))
                {
                    return false;
                }
            }
        }
    }
    const uint64_t comparisons =
        wirehair_v2::HeavyProjectionOracleComparisonsForTesting();
    const uint64_t fallbacks =
        wirehair_v2::HeavyProjectionLegacyFallbacksForTesting();
    if (comparisons != 0u || fallbacks == 0u)
    {
        std::fprintf(stderr,
            "solve: heavy projection fallback boundary mismatch "
            "comparisons=%llu fallbacks=%llu\n",
            (unsigned long long)comparisons,
            (unsigned long long)fallbacks);
        return false;
    }
    std::printf(
        "heavy H=0/1/11/12/13/128 and K=3/4 D/H candidates "
        "periodic/hashed dense oracle fallbacks=%llu: PASS\n",
        (unsigned long long)fallbacks);
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

enum class ResumePacketAliasTarget
{
    CoefficientScratch,
    RhsScratch
};

struct ResumePacketAliasCase
{
    const char* Name;
    ResumePacketAliasTarget Target;
    size_t Offset;
};

bool CheckResumePacketAliasCall(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PrecodeSolveResumeState& initial_state,
    uint32_t block_id,
    const uint8_t* packet_data,
    uint32_t block_bytes,
    const ResumePacketAliasCase& alias_case,
    bool allow_insert,
    WirehairResult expected_result,
    const std::vector<uint8_t>& expected_output)
{
    wirehair_v2::PrecodeSolveResumeState control = initial_state;
    wirehair_v2::PrecodeSolveResumeState aliased = initial_state;
    const auto stage_packet = [&](
        wirehair_v2::PrecodeSolveResumeState& state,
        const uint8_t** aliased_data) -> bool
    {
        std::vector<uint8_t>& storage =
            alias_case.Target == ResumePacketAliasTarget::CoefficientScratch ?
                state.CoefficientScratch : state.RhsScratch;
        if (alias_case.Offset > storage.size() ||
            block_bytes > storage.size() - alias_case.Offset)
        {
            return false;
        }
        std::memcpy(
            storage.data() + alias_case.Offset,
            packet_data,
            block_bytes);
        if (aliased_data) {
            *aliased_data = storage.data() + alias_case.Offset;
        }
        return true;
    };

    const uint8_t* aliased_data = nullptr;
    if (!stage_packet(control, nullptr) ||
        !stage_packet(aliased, &aliased_data))
    {
        std::fprintf(stderr,
            "solve: resume packet alias fixture too small bb=%u case=%s\n",
            block_bytes, alias_case.Name);
        return false;
    }
    const wirehair_v2::PrecodeSolveResumeState aliased_before = aliased;
    const ResumeStorageIdentity storage_before =
        CaptureResumeStorageIdentity(aliased);
    const size_t persistent_before = aliased.PersistentBytes();

    std::vector<uint8_t> control_output(11u, 0xa5u);
    std::vector<uint8_t> aliased_output = control_output;
    wirehair_v2::PrecodeSolveStats initial_stats = initial_state.Stats;
    initial_stats.PacketRows ^= UINT32_C(0x86b39a27);
    initial_stats.BlockXors ^= UINT64_C(0xc1d2e3f405162738);
    wirehair_v2::PrecodeSolveStats control_stats = initial_stats;
    wirehair_v2::PrecodeSolveStats aliased_stats = initial_stats;
    const WirehairResult control_result =
        wirehair_v2::ResumePrecodeSystem(
            system, config, block_id, packet_data, block_bytes,
            control, control_output, &control_stats, allow_insert);
    const WirehairResult aliased_result =
        wirehair_v2::ResumePrecodeSystem(
            system, config, block_id, aliased_data, block_bytes,
            aliased, aliased_output, &aliased_stats, allow_insert);

    if (control_result != expected_result ||
        aliased_result != control_result ||
        control_output != expected_output ||
        aliased_output != control_output ||
        !SameResumeStateIgnoringTiming(control, aliased) ||
        !SameSolveStatsIgnoringTiming(control_stats, aliased_stats))
    {
        std::fprintf(stderr,
            "solve: resume packet alias mismatch bb=%u case=%s "
            "insert=%u expected=%d control=%d alias=%d "
            "control_rank=%u alias_rank=%u\n",
            block_bytes, alias_case.Name, allow_insert ? 1u : 0u,
            (int)expected_result, (int)control_result,
            (int)aliased_result, control.Rank, aliased.Rank);
        return false;
    }
    if (!allow_insert &&
        (!SameResumeState(aliased, aliased_before) ||
         !SameResumeStorageIdentity(aliased, storage_before) ||
         aliased.PersistentBytes() != persistent_before ||
         !SameSolveStats(aliased_stats, initial_stats)))
    {
        std::fprintf(stderr,
            "solve: checked resume packet alias changed checkpoint "
            "bb=%u case=%s result=%d\n",
            block_bytes, alias_case.Name, (int)aliased_result);
        return false;
    }
    return true;
}

bool CheckResumePacketInputAliases(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const std::vector<uint8_t>& message,
    uint32_t block_bytes,
    const wirehair_v2::PrecodeSolveResumeState& initial_state,
    const std::vector<uint8_t>& expected)
{
    if (initial_state.RhsScratch.size() != block_bytes)
    {
        std::fprintf(stderr,
            "solve: resume packet alias fixture dimensions bb=%u R=%zu\n",
            block_bytes, initial_state.CoefficientScratch.size());
        return false;
    }
    std::vector<ResumePacketAliasCase> alias_cases;
    if (initial_state.CoefficientScratch.size() >=
        (size_t)block_bytes + 2u)
    {
        const size_t right_boundary =
            initial_state.CoefficientScratch.size() - block_bytes;
        alias_cases.push_back({
            "coeff-exact",
            ResumePacketAliasTarget::CoefficientScratch,
            0u });
        alias_cases.push_back({
            "coeff-left-offset",
            ResumePacketAliasTarget::CoefficientScratch,
            1u });
        alias_cases.push_back({
            "coeff-right-boundary",
            ResumePacketAliasTarget::CoefficientScratch,
            right_boundary });
    }
    alias_cases.push_back({
        "rhs-exact", ResumePacketAliasTarget::RhsScratch, 0u });
    const std::vector<uint8_t> sentinel(11u, 0xa5u);
    std::vector<uint8_t> conflicting(
        message.begin(), message.begin() + block_bytes);
    conflicting[0] ^= 1u;

    for (const ResumePacketAliasCase& alias_case : alias_cases)
    {
        if (!CheckResumePacketAliasCall(
                system, config, initial_state,
                1u, message.data() + block_bytes, block_bytes,
                alias_case, true, Wirehair_NeedMore, sentinel) ||
            !CheckResumePacketAliasCall(
                system, config, initial_state,
                0u, message.data(), block_bytes,
                alias_case, false, Wirehair_NeedMore, sentinel) ||
            !CheckResumePacketAliasCall(
                system, config, initial_state,
                0u, conflicting.data(), block_bytes,
                alias_case, false, Wirehair_Error, sentinel))
        {
            return false;
        }
    }

    wirehair_v2::PrecodeSolveResumeState near_complete = initial_state;
    std::vector<uint8_t> output = sentinel;
    wirehair_v2::PrecodeSolveStats stats = initial_state.Stats;
    for (uint32_t id = 1u; id + 1u < initial_state.SourceCount; ++id)
    {
        if (wirehair_v2::ResumePrecodeSystem(
                system, config, id,
                message.data() + (size_t)id * block_bytes,
                block_bytes, near_complete, output, &stats, true) !=
                Wirehair_NeedMore ||
            output != sentinel)
        {
            std::fprintf(stderr,
                "solve: resume packet alias completion fixture failed "
                "bb=%u id=%u\n",
                block_bytes, id);
            return false;
        }
    }
    const uint32_t final_id = initial_state.SourceCount - 1u;
    for (const ResumePacketAliasCase& alias_case : alias_cases)
    {
        if (!CheckResumePacketAliasCall(
                system, config, near_complete,
                final_id,
                message.data() + (size_t)final_id * block_bytes,
                block_bytes, alias_case, true,
                Wirehair_Success, expected))
        {
            return false;
        }
    }
    std::printf(
        "resume packet scratch aliases bb=%u: PASS\n", block_bytes);
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

    std::vector<wirehair_v2::SolvePacket> deficient(K);
    for (wirehair_v2::SolvePacket& packet : deficient) {
        packet.BlockId = 0u;
        packet.Data = message.data();
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

    typedef std::vector<uint8_t>
        wirehair_v2::PrecodeSolveResumeState::* ResumeByteVectorMember;
    static const ResumeByteVectorMember resume_byte_vectors[] = {
        &wirehair_v2::PrecodeSolveResumeState::Values,
        &wirehair_v2::PrecodeSolveResumeState::PivotCoefficients,
        &wirehair_v2::PrecodeSolveResumeState::PivotRhs,
        &wirehair_v2::PrecodeSolveResumeState::HavePivot,
        &wirehair_v2::PrecodeSolveResumeState::CoefficientScratch,
        &wirehair_v2::PrecodeSolveResumeState::RhsScratch
    };
    static const char* const resume_byte_vector_names[] = {
        "Values",
        "PivotCoefficients",
        "PivotRhs",
        "HavePivot",
        "CoefficientScratch",
        "RhsScratch"
    };
    static_assert(
        sizeof(resume_byte_vectors) / sizeof(resume_byte_vectors[0]) ==
            sizeof(resume_byte_vector_names) /
                sizeof(resume_byte_vector_names[0]),
        "resume byte-vector alias names must stay complete");
    for (size_t alias_i = 0u;
         alias_i < sizeof(resume_byte_vectors) /
             sizeof(resume_byte_vectors[0]);
         ++alias_i)
    {
        // Cold checkpoint publication replaces an existing resume state.  Its
        // output cannot be one of the vectors being replaced.
        wirehair_v2::PrecodeSolveResumeState cold_alias = resume;
        std::vector<uint8_t>& cold_output =
            cold_alias.*resume_byte_vectors[alias_i];
        cold_output.assign(alias_i + 3u, (uint8_t)(0x80u + alias_i));
        const wirehair_v2::PrecodeSolveResumeState cold_before = cold_alias;
        const ResumeStorageIdentity cold_storage_before =
            CaptureResumeStorageIdentity(cold_alias);
        const std::vector<uint8_t> cold_output_before = cold_output;
        const size_t cold_persistent_before = cold_alias.PersistentBytes();
        wirehair_v2::PrecodeSolveStats cold_stats = stats;
        cold_stats.PacketRows ^= UINT32_C(0x13579bdf);
        cold_stats.BlockXors ^= UINT64_C(0x2468ace013579bdf);
        const wirehair_v2::PrecodeSolveStats cold_stats_before = cold_stats;
        if (wirehair_v2::SolvePrecodeSystem(
                system, config, deficient, block_bytes,
                cold_output, &cold_stats, &cold_alias) !=
                Wirehair_InvalidInput ||
            !SameResumeState(cold_alias, cold_before) ||
            !SameResumeStorageIdentity(cold_alias, cold_storage_before) ||
            cold_alias.PersistentBytes() != cold_persistent_before ||
            cold_output != cold_output_before ||
            !SameSolveStats(cold_stats, cold_stats_before))
        {
            std::fprintf(stderr,
                "solve: cold checkpoint/output alias changed state "
                "bb=%u member=%s\n",
                block_bytes, resume_byte_vector_names[alias_i]);
            return false;
        }

        // Resume completion clears the checkpoint after publishing.  Reject
        // aliases before row generation, scratch writes, or output allocation.
        wirehair_v2::PrecodeSolveResumeState resume_alias = resume;
        std::vector<uint8_t>& resume_output =
            resume_alias.*resume_byte_vectors[alias_i];
        const wirehair_v2::PrecodeSolveResumeState resume_before =
            resume_alias;
        const ResumeStorageIdentity resume_storage_before =
            CaptureResumeStorageIdentity(resume_alias);
        const std::vector<uint8_t> resume_output_before = resume_output;
        const size_t resume_persistent_before =
            resume_alias.PersistentBytes();
        wirehair_v2::PrecodeSolveStats resume_stats = stats;
        resume_stats.PacketRows ^= UINT32_C(0x89abcdef);
        resume_stats.BlockMulAdds ^= UINT64_C(0xfedcba9876543210);
        const wirehair_v2::PrecodeSolveStats resume_stats_before =
            resume_stats;
        if (wirehair_v2::ResumePrecodeSystem(
                system, config, 1u,
                message.data() + block_bytes, block_bytes,
                resume_alias, resume_output, &resume_stats, true) !=
                Wirehair_InvalidInput ||
            !SameResumeState(resume_alias, resume_before) ||
            !SameResumeStorageIdentity(
                resume_alias, resume_storage_before) ||
            resume_alias.PersistentBytes() != resume_persistent_before ||
            resume_output != resume_output_before ||
            !SameSolveStats(resume_stats, resume_stats_before))
        {
            std::fprintf(stderr,
                "solve: incremental checkpoint/output alias changed state "
                "bb=%u member=%s\n",
                block_bytes, resume_byte_vector_names[alias_i]);
            return false;
        }
    }

    if (!CheckResumePacketInputAliases(
            system, config, message, block_bytes, resume, expected))
    {
        return false;
    }

    // A diagnostic pointer into the checkpoint itself cannot satisfy the
    // resume output contract.  Exercise paths that would otherwise validate a
    // dependent row, reject a conflicting row, insert an independent row, or
    // fail validation.  Every one must reject before row generation, scratch
    // allocation, or state mutation.
    std::vector<uint8_t> alias_corrupt(
        message.begin(), message.begin() + block_bytes);
    alias_corrupt[0] ^= 1u;
    wirehair_v2::PrecodeSolveResumeState independent_control = resume;
    wirehair_v2::PrecodeSolveStats independent_stats = {};
    std::vector<uint8_t> independent_output(5u, 0xc4u);
    const std::vector<uint8_t> independent_output_before =
        independent_output;
    if (wirehair_v2::ResumePrecodeSystem(
            system, config, 1u, message.data() + block_bytes, block_bytes,
            independent_control, independent_output, &independent_stats,
            true) != Wirehair_NeedMore ||
        independent_control.Rank <= resume.Rank ||
        !SameSolveStats(independent_stats, independent_control.Stats) ||
        independent_output != independent_output_before)
    {
        std::fprintf(stderr,
            "solve: independent resume control did not advance rank bb=%u\n",
            block_bytes);
        return false;
    }
    if (!ExpectResumeStatsAliasRejected(
            "dependent-check", system, config, 0u, message.data(),
            block_bytes, resume, false) ||
        !ExpectResumeStatsAliasRejected(
            "conflicting-check", system, config, 0u,
            alias_corrupt.data(), block_bytes, resume, false) ||
        !ExpectResumeStatsAliasRejected(
            "dependent-insert", system, config, 0u, message.data(),
            block_bytes, resume, true) ||
        !ExpectResumeStatsAliasRejected(
            "independent-insert", system, config, 1u,
            message.data() + block_bytes, block_bytes, resume, true))
    {
        return false;
    }
    wirehair_v2::PrecodeSolveResumeState invalid_resume = resume;
    invalid_resume.Active = false;
    if (!ExpectResumeStatsAliasRejected(
            "invalid-checkpoint", system, config, 1u,
            message.data() + block_bytes, block_bytes,
            invalid_resume, true))
    {
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
    wirehair_v2::PrecodeSolveResumeState success_alias_state;
    uint32_t success_alias_id = UINT32_MAX;
    for (uint32_t id = 1u; id < K; ++id)
    {
        wirehair_v2::PrecodeSolveResumeState state_before_success;
        if (block_bytes == 17u) {
            state_before_success = resume;
        }
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
        if (block_bytes == 17u && result == Wirehair_Success)
        {
            success_alias_state.Swap(state_before_success);
            success_alias_id = id;
        }
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
    if (block_bytes == 17u &&
        (success_alias_id == UINT32_MAX ||
         !ExpectResumeStatsAliasRejected(
             "success", system, config, success_alias_id,
             message.data() + (size_t)success_alias_id * block_bytes,
             block_bytes, success_alias_state, true)))
    {
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

struct HeavyProjectionCase
{
    const char* Name;
    uint32_t K;
    uint32_t HeavyRows;
    wirehair_v2::DenseAnchorLayout DenseAnchors;
};

bool CheckHeavyProjectionCase(const HeavyProjectionCase& test_case)
{
    const uint32_t block_bytes = 1u;
    wirehair_v2::PrecodeParams params = wirehair_v2::MakeCertifiedParams(
        test_case.K,
        UINT64_C(0x4850524f4a454354) ^
            ((uint64_t)test_case.K << 17) ^
            ((uint64_t)test_case.HeavyRows << 8) ^
            (uint32_t)test_case.DenseAnchors);
    params.HeavyRows = test_case.HeavyRows;
    params.DenseAnchors = test_case.DenseAnchors;
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x4850524f) ^ test_case.K ^
        (test_case.HeavyRows << 16) ^ (uint32_t)test_case.DenseAnchors;
    base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config) != Wirehair_Success ||
        system.Params.HeavyRows != test_case.HeavyRows ||
        system.Params.DenseAnchors != test_case.DenseAnchors)
    {
        std::fprintf(stderr,
            "solve: heavy projection configuration failed case=%s\n",
            test_case.Name);
        return false;
    }

    std::vector<uint8_t> message(test_case.K * block_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(
            i * 109u + (i >> 5) + test_case.K * 3u +
            test_case.HeavyRows * 7u + (uint32_t)test_case.DenseAnchors);
    }
    std::vector<wirehair_v2::SolvePacket> packets(test_case.K);
    for (uint32_t id = 0u; id < test_case.K; ++id) {
        packets[id].BlockId = id;
        packets[id].Data = message.data() + (size_t)id * block_bytes;
    }

    std::vector<uint8_t> control;
    std::vector<uint8_t> checked;
    wirehair_v2::PrecodeSolveStats control_stats;
    wirehair_v2::PrecodeSolveStats checked_stats;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, packets, block_bytes, control,
            &control_stats) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "solve: heavy projection control failed case=%s\n",
            test_case.Name);
        return false;
    }

    const uint64_t comparisons_before =
        wirehair_v2::HeavyProjectionOracleComparisonsForTesting();
    const uint64_t fallbacks_before =
        wirehair_v2::HeavyProjectionLegacyFallbacksForTesting();
    WirehairResult checked_result = Wirehair_Error;
    {
        HeavyProjectionOracleScope oracle_scope;
        checked_result = wirehair_v2::SolvePrecodeSystem(
            system, config, packets, block_bytes, checked, &checked_stats);
    }
    const uint64_t comparison_delta =
        wirehair_v2::HeavyProjectionOracleComparisonsForTesting() -
        comparisons_before;
    const uint64_t fallback_delta =
        wirehair_v2::HeavyProjectionLegacyFallbacksForTesting() -
        fallbacks_before;
    const uint32_t L = test_case.K + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const bool optimized =
        system.Params.HeavyFamily ==
            wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy &&
        system.Params.HeavyRows == 12u && L >= 244u &&
        control_stats.InactivatedColumns != 0u;
    const uint64_t expected_comparisons = optimized ? 1u : 0u;
    const uint64_t expected_fallbacks =
        !optimized && control_stats.InactivatedColumns != 0u ? 1u : 0u;
    if (checked_result != Wirehair_Success || checked != control ||
        !SameSolveStatsIgnoringTiming(control_stats, checked_stats) ||
        comparison_delta != expected_comparisons ||
        fallback_delta != expected_fallbacks ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, packets, checked.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: heavy projection oracle mismatch case=%s result=%d "
            "comparisons=%llu/%llu fallbacks=%llu/%llu R=%u L=%u\n",
            test_case.Name,
            (int)checked_result,
            (unsigned long long)comparison_delta,
            (unsigned long long)expected_comparisons,
            (unsigned long long)fallback_delta,
            (unsigned long long)expected_fallbacks,
            control_stats.InactivatedColumns,
            L);
        return false;
    }

    static const uint32_t kSampleNumerators[] = { 0u, 1u, 2u };
    uint8_t recovered = 0u;
    for (uint32_t numerator : kSampleNumerators)
    {
        const uint32_t id = numerator == 0u ? 0u :
            (numerator == 1u ? test_case.K / 2u : test_case.K - 1u);
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, checked.data(), block_bytes, id,
                &recovered) ||
            recovered != message[id])
        {
            std::fprintf(stderr,
                "solve: heavy projection replay mismatch case=%s id=%u\n",
                test_case.Name, id);
            return false;
        }
    }
    std::printf(
        "heavy projection case=%s K=%u H=%u layout=%u "
        "comparisons=%llu fallbacks=%llu: PASS\n",
        test_case.Name,
        test_case.K,
        test_case.HeavyRows,
        (unsigned)test_case.DenseAnchors,
        (unsigned long long)comparison_delta,
        (unsigned long long)fallback_delta);
    return true;
}

bool CheckHeavyProjectionResumeOracle()
{
    const uint32_t K = 320u;
    const uint32_t block_bytes = 17u;
    wirehair_v2::PrecodeParams params = wirehair_v2::MakeCertifiedParams(
        K, UINT64_C(0x4850524f4a52534d));
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x4a52534d);
    base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "solve: heavy projection resume configuration failed\n");
        return false;
    }

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 181u + (i >> 4) + 0x37u);
    }
    std::vector<wirehair_v2::SolvePacket> systematic(K);
    for (uint32_t id = 0u; id < K; ++id) {
        systematic[id].BlockId = id;
        systematic[id].Data =
            message.data() + (size_t)id * block_bytes;
    }
    std::vector<uint8_t> expected;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, systematic, block_bytes, expected) !=
            Wirehair_Success)
    {
        std::fprintf(stderr,
            "solve: heavy projection resume control failed\n");
        return false;
    }

    std::vector<wirehair_v2::SolvePacket> deficient(
        K, systematic.front());
    wirehair_v2::PrecodeSolveResumeState resume;
    wirehair_v2::PrecodeSolveStats stats;
    std::vector<uint8_t> output(7u, uint8_t{0x5a});
    const std::vector<uint8_t> sentinel = output;
    const uint64_t comparisons_before =
        wirehair_v2::HeavyProjectionOracleComparisonsForTesting();
    const uint64_t fallbacks_before =
        wirehair_v2::HeavyProjectionLegacyFallbacksForTesting();
    WirehairResult result = Wirehair_Error;
    {
        HeavyProjectionOracleScope oracle_scope;
        result = wirehair_v2::SolvePrecodeSystem(
            system, config, deficient, block_bytes, output, &stats, &resume);
    }
    const uint64_t comparison_delta =
        wirehair_v2::HeavyProjectionOracleComparisonsForTesting() -
        comparisons_before;
    const uint64_t fallback_delta =
        wirehair_v2::HeavyProjectionLegacyFallbacksForTesting() -
        fallbacks_before;
    if (result != Wirehair_NeedMore || !resume.Active ||
        resume.Rank >= resume.InactiveCount || output != sentinel ||
        comparison_delta != 1u || fallback_delta != 0u)
    {
        std::fprintf(stderr,
            "solve: heavy projection resume checkpoint failed result=%d "
            "rank=%u/%u comparisons=%llu fallbacks=%llu\n",
            (int)result,
            resume.Rank,
            resume.InactiveCount,
            (unsigned long long)comparison_delta,
            (unsigned long long)fallback_delta);
        return false;
    }

    for (uint32_t id = 1u; id < K && result == Wirehair_NeedMore; ++id)
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
    }
    if (result != Wirehair_Success || resume.Active || output != expected ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, systematic, output.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: heavy projection resumed output mismatch result=%d\n",
            (int)result);
        return false;
    }
    std::printf(
        "heavy projection deficient resume K=%u comparisons=%llu: PASS\n",
        K,
        (unsigned long long)comparison_delta);
    return true;
}

bool CheckHeavyProjectionPropagationOracle()
{
    static const HeavyProjectionCase kCases[] = {
        { "tiny-h12", 8u, 12u,
          wirehair_v2::DenseAnchorLayout::Disabled },
        { "middle-production", 1000u, 12u,
          wirehair_v2::DenseAnchorLayout::Disabled },
        { "middle-two07", 1000u, 12u,
          wirehair_v2::DenseAnchorLayout::Two07 },
        { "large-production", 64000u, 12u,
          wirehair_v2::DenseAnchorLayout::Disabled },
        { "fallback-h0", 320u, 0u,
          wirehair_v2::DenseAnchorLayout::Disabled },
        { "fallback-h13", 320u, 13u,
          wirehair_v2::DenseAnchorLayout::Disabled },
        { "fallback-h128", 320u, 128u,
          wirehair_v2::DenseAnchorLayout::Disabled }
    };
    wirehair_v2::ResetHeavyProjectionOracleCountersForTesting();
    for (const HeavyProjectionCase& test_case : kCases) {
        if (!CheckHeavyProjectionCase(test_case)) {
            return false;
        }
    }
    if (!CheckHeavyProjectionResumeOracle()) {
        return false;
    }
    std::printf("packed H=12 peel-schedule heavy projection oracle: PASS\n");
    return true;
}

bool CheckColdSolveStatsAlias()
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 17u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(
            K, UINT64_C(0x5354415453414c49));
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed = UINT32_C(0x5a17c0de);
    base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PrecodeSystem system;
    wirehair_v2::PacketRowConfig config;
    if (wirehair_v2::SelectSystematicConfiguration(
            params, base_config, system, config) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "solve: cold stats-alias configuration failed\n");
        return false;
    }
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 181u + (i >> 4) + 37u);
    }
    std::vector<wirehair_v2::SolvePacket> systematic(K);
    for (uint32_t id = 0u; id < K; ++id) {
        systematic[id].BlockId = id;
        systematic[id].Data =
            message.data() + (size_t)id * block_bytes;
    }
    std::vector<wirehair_v2::SolvePacket> deficient(K);
    for (wirehair_v2::SolvePacket& packet : deficient) {
        packet.BlockId = 0u;
        packet.Data = message.data();
    }
    std::vector<wirehair_v2::SolvePacket> too_few(
        systematic.begin(), systematic.end() - 1);

    static const ColdSolveEntryPoint entry_points[] = {
        ColdSolveEntryPoint::Public,
        ColdSolveEntryPoint::Runtime,
        ColdSolveEntryPoint::ValidatedRuntime
    };
    for (size_t entry_i = 0u;
         entry_i < sizeof(entry_points) / sizeof(entry_points[0]);
         ++entry_i)
    {
        const ColdSolveEntryPoint entry_point = entry_points[entry_i];
        const uint32_t tag = (uint32_t)entry_i + 1u;

        // A disjoint diagnostic object retains the ordinary successful path,
        // while an unused checkpoint remains byte-for-byte and allocation-
        // identity stable.
        wirehair_v2::PrecodeSolveResumeState success_state =
            MakeStatsAliasSentinel(UINT32_C(0x100) + tag);
        const wirehair_v2::PrecodeSolveResumeState success_before =
            success_state;
        const ResumeStorageIdentity success_storage_before =
            CaptureResumeStorageIdentity(success_state);
        const size_t success_persistent_before =
            success_state.PersistentBytes();
        wirehair_v2::PrecodeSolveStats success_stats = {};
        std::vector<uint8_t> success_output(9u, 0x6du);
        if (CallColdSolve(
                entry_point, system, config, runtime, systematic,
                block_bytes, success_output, &success_stats,
                &success_state) != Wirehair_Success ||
            !wirehair_v2::VerifyPrecodeSolution(
                system, config, systematic,
                success_output.data(), block_bytes) ||
            success_stats.PacketRows != K ||
            !SameResumeState(success_state, success_before) ||
            !SameResumeStorageIdentity(
                success_state, success_storage_before) ||
            success_state.PersistentBytes() != success_persistent_before)
        {
            std::fprintf(stderr,
                "solve: disjoint cold success changed checkpoint entry=%s\n",
                ColdSolveEntryPointName(entry_point));
            return false;
        }
        if (!ExpectColdSolveStatsAliasRejected(
                "success", entry_point, system, config, runtime,
                systematic, block_bytes, UINT32_C(0x200) + tag))
        {
            return false;
        }

        // Rank deficiency may publish a new checkpoint and the same counters
        // to two distinct objects.  The output remains exactly untouched.
        wirehair_v2::PrecodeSolveResumeState deficient_state =
            MakeStatsAliasSentinel(UINT32_C(0x300) + tag);
        wirehair_v2::PrecodeSolveStats deficient_stats = {};
        std::vector<uint8_t> deficient_output(10u, 0x7eu);
        const std::vector<uint8_t> deficient_output_before =
            deficient_output;
        const uint8_t* const deficient_output_data_before =
            deficient_output.data();
        const size_t deficient_output_capacity_before =
            deficient_output.capacity();
        if (CallColdSolve(
                entry_point, system, config, runtime, deficient,
                block_bytes, deficient_output, &deficient_stats,
                &deficient_state) != Wirehair_NeedMore ||
            !deficient_state.Active ||
            deficient_state.Rank >= deficient_state.InactiveCount ||
            !SameSolveStats(deficient_stats, deficient_state.Stats) ||
            deficient_output != deficient_output_before ||
            deficient_output.data() != deficient_output_data_before ||
            deficient_output.capacity() !=
                deficient_output_capacity_before)
        {
            std::fprintf(stderr,
                "solve: disjoint deficient checkpoint mismatch entry=%s\n",
                ColdSolveEntryPointName(entry_point));
            return false;
        }
        if (!ExpectColdSolveStatsAliasRejected(
                "deficient", entry_point, system, config, runtime,
                deficient, block_bytes, UINT32_C(0x400) + tag))
        {
            return false;
        }

        // The pre-allocation too-few-packets result leaves disjoint
        // diagnostics and checkpoint state alone; alias rejection has higher
        // precedence because that output contract is impossible to honor.
        wirehair_v2::PrecodeSolveResumeState short_state =
            MakeStatsAliasSentinel(UINT32_C(0x500) + tag);
        const wirehair_v2::PrecodeSolveResumeState short_before = short_state;
        const ResumeStorageIdentity short_storage_before =
            CaptureResumeStorageIdentity(short_state);
        const size_t short_persistent_before = short_state.PersistentBytes();
        wirehair_v2::PrecodeSolveStats short_stats = short_state.Stats;
        const wirehair_v2::PrecodeSolveStats short_stats_before = short_stats;
        std::vector<uint8_t> short_output(11u, 0x8fu);
        const std::vector<uint8_t> short_output_before = short_output;
        if (CallColdSolve(
                entry_point, system, config, runtime, too_few,
                block_bytes, short_output, &short_stats,
                &short_state) != Wirehair_NeedMore ||
            !SameResumeState(short_state, short_before) ||
            !SameResumeStorageIdentity(short_state, short_storage_before) ||
            short_state.PersistentBytes() != short_persistent_before ||
            !SameSolveStats(short_stats, short_stats_before) ||
            short_output != short_output_before)
        {
            std::fprintf(stderr,
                "solve: disjoint pre-allocation NeedMore mutated output "
                "entry=%s\n", ColdSolveEntryPointName(entry_point));
            return false;
        }
        if (!ExpectColdSolveStatsAliasRejected(
                "pre-allocation-need-more", entry_point,
                system, config, runtime, too_few, block_bytes,
                UINT32_C(0x600) + tag))
        {
            return false;
        }
    }

    // Public and validating-runtime entry points preserve all caller outputs
    // on stored-graph validation failures.  The trusted overload deliberately
    // skips this validation and is covered only by its alias fast rejection.
    wirehair_v2::PrecodeSystem invalid_system = system;
    if (invalid_system.StaircaseRows.empty()) {
        return false;
    }
    invalid_system.StaircaseRows.pop_back();
    static const ColdSolveEntryPoint validating_entry_points[] = {
        ColdSolveEntryPoint::Public,
        ColdSolveEntryPoint::Runtime
    };
    for (size_t entry_i = 0u;
         entry_i < sizeof(validating_entry_points) /
             sizeof(validating_entry_points[0]);
         ++entry_i)
    {
        const ColdSolveEntryPoint entry_point =
            validating_entry_points[entry_i];
        wirehair_v2::PrecodeSolveResumeState state =
            MakeStatsAliasSentinel(UINT32_C(0x700) + (uint32_t)entry_i);
        const wirehair_v2::PrecodeSolveResumeState state_before = state;
        const ResumeStorageIdentity storage_before =
            CaptureResumeStorageIdentity(state);
        const size_t persistent_before = state.PersistentBytes();
        wirehair_v2::PrecodeSolveStats stats = state.Stats;
        const wirehair_v2::PrecodeSolveStats stats_before = stats;
        std::vector<uint8_t> output(12u, 0x91u);
        const std::vector<uint8_t> output_before = output;
        if (CallColdSolve(
                entry_point, invalid_system, config, runtime, systematic,
                block_bytes, output, &stats, &state) !=
                Wirehair_InvalidInput ||
            !SameResumeState(state, state_before) ||
            !SameResumeStorageIdentity(state, storage_before) ||
            state.PersistentBytes() != persistent_before ||
            !SameSolveStats(stats, stats_before) || output != output_before)
        {
            std::fprintf(stderr,
                "solve: validation failure mutated outputs entry=%s\n",
                ColdSolveEntryPointName(entry_point));
            return false;
        }
    }
    for (size_t entry_i = 0u;
         entry_i < sizeof(entry_points) / sizeof(entry_points[0]);
         ++entry_i)
    {
        if (!ExpectColdSolveStatsAliasRejected(
                "invalid-system", entry_points[entry_i],
                invalid_system, config, runtime, systematic, block_bytes,
                UINT32_C(0x800) + (uint32_t)entry_i))
        {
            return false;
        }
    }

    std::printf(
        "cold solve stats/checkpoint alias across all entry points: PASS\n");
    return true;
}

bool CheckBinaryPeelLowDegreeXorOracle()
{
    const uint32_t K = 64000u;
    const uint32_t block_bytes = 2u;
    wirehair_v2::PrecodeParams params = wirehair_v2::MakeCertifiedParams(
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

bool CheckConcurrentCoefficientCache()
{
    const uint32_t K = 320u;
    const uint32_t block_bytes = 37u;
    const wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    const uint64_t matrix_seed = wirehair_v2::MatrixSeedFromProfile(
        profile, 0u, wirehair_v2::kMessagePrecodeSeedSalt);
    wirehair_v2::PrecodeParams params =
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
            "solve: concurrent H12 cache thread launch failed\n");
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
            "solve: concurrent H12 coefficient-cache solve failed\n");
        return false;
    }
    for (uint32_t thread = 1u; thread < kThreadCount; ++thread)
    {
        if (outputs[thread] != outputs[0]) {
            std::fprintf(stderr,
                "solve: concurrent H12 coefficient-cache mismatch\n");
            return false;
        }
    }
    return true;
}

bool CheckConcurrentCoefficientCaches()
{
    if (!CheckConcurrentCoefficientCache()) {
        return false;
    }
    std::printf("concurrent H12 coefficient-cache first use: PASS\n");
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
    const uint8_t* const output_data_before = output.data();
    const size_t output_capacity_before = output.capacity();
    wirehair_v2::PrecodeSolveResumeState resume =
        MakeStatsAliasSentinel(UINT32_C(0x900));
    const wirehair_v2::PrecodeSolveResumeState resume_before = resume;
    const ResumeStorageIdentity storage_before =
        CaptureResumeStorageIdentity(resume);
    const size_t persistent_before = resume.PersistentBytes();
    wirehair_v2::PrecodeSolveStats stats = {};
    const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
        system, config, packets, 1u, output, &stats, &resume);
    if (result != Wirehair_NeedMore || output != before ||
        output.data() != output_data_before ||
        output.capacity() != output_capacity_before ||
        !SameResumeState(resume, resume_before) ||
        !SameResumeStorageIdentity(resume, storage_before) ||
        resume.PersistentBytes() != persistent_before ||
        stats.InactivatedColumns <= wirehair_v2::kMaxInactiveColumns ||
        stats.PeelNanoseconds == 0u)
    {
        std::fprintf(stderr,
            "solve: inactive cap failed result=%d inact=%u\n",
            (int)result, stats.InactivatedColumns);
        return false;
    }

    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }
    static const ColdSolveEntryPoint entry_points[] = {
        ColdSolveEntryPoint::Public,
        ColdSolveEntryPoint::Runtime,
        ColdSolveEntryPoint::ValidatedRuntime
    };
    for (size_t entry_i = 0u;
         entry_i < sizeof(entry_points) / sizeof(entry_points[0]);
         ++entry_i)
    {
        if (!ExpectColdSolveStatsAliasRejected(
                "inactive-cap", entry_points[entry_i],
                system, config, runtime, packets, 1u,
                UINT32_C(0xa00) + (uint32_t)entry_i))
        {
            return false;
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
                    diagnostic_system.DenseBasisRowColumns) {
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    ok = CheckDenseAnchorPayloadOracle() && ok;
#endif
    ok = CheckExactSystematicFailureClassification() && ok;
    ok = CheckHeavyCoefficientBoundaryOracle() && ok;
    ok = CheckBinaryQuotientBoundary() && ok;
    ok = CheckConcurrentCoefficientCaches() && ok;
    ok = CheckIncrementalResume() && ok;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    ok = CheckHeavyProjectionPropagationOracle() && ok;
#endif
    ok = CheckColdSolveStatsAlias() && ok;
    ok = CheckBinaryPeelLowDegreeXorOracle() && ok;
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
