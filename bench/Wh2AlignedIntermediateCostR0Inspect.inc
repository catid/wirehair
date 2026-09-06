// Appended verbatim to BOTH complete facade translation units, never standalone.
#if !defined(WIREHAIR_BUILDING) || WIREHAIR_BUILDING != 1 || \
    !defined(WIREHAIR_STATIC) || WIREHAIR_STATIC != 1 || !defined(__PIC__)
#error "Inspection requires the original static PIC facade compile state"
#endif
#if defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS) || defined(WH_SEED_KNOBS) || \
    defined(WIREHAIR_HPP)
#error "Inspection forbids hooks, experimental equations and C++ facade wrappers"
#endif
#if !defined(WH2_ALIGNMENT_INSPECT_FN) || !defined(WH2_ALIGNMENT_COLUMNS_FN) || \
    !defined(WH2_ALIGNMENT_SCALAR_FN) || !defined(WH2_ALIGNMENT_INSPECT_ALIGNED)
#error "Inspection export names and storage type must be explicit"
#endif
#if WH2_ALIGNMENT_INSPECT_ALIGNED != 0 && WH2_ALIGNMENT_INSPECT_ALIGNED != 1
#error "Invalid inspection storage variant"
#endif

#include "Wh2AlignedIntermediateCostR0Bridge.h"
#include "Wh2AlignedIntermediateR0.h"

namespace {

// Local types are essential: the two TUs refer to different V2 namespaces.
struct CostPrecodeTag
{
    typedef std::unique_ptr<wirehair_v2::MessagePrecodeEncoder>
        wirehair_v2::Codec::* Type;
    friend Type CostMember(CostPrecodeTag);
};

struct CostStorageTag
{
#if WH2_ALIGNMENT_INSPECT_ALIGNED
    typedef wh2_aligned_intermediate_r0::ByteVector Storage;
#else
    typedef std::vector<uint8_t> Storage;
#endif
    typedef Storage wirehair_v2::PrecodeEncoder::* Type;
    friend Type CostMember(CostStorageTag);
};

template<class Tag, typename Tag::Type member>
struct CostAccess
{
    friend typename Tag::Type CostMember(Tag) { return member; }
};
template struct CostAccess<CostPrecodeTag, &wirehair_v2::Codec::PrecodeImpl>;
template struct CostAccess<CostStorageTag,
    &wirehair_v2::PrecodeEncoder::SolvedIntermediateStorage>;

struct CostInspection
{
    Wh2AlignedIntermediateCostSnapshot pod = {};
    const PublicCodec* facade = nullptr;
    const wirehair_v2::MessagePrecodeEncoder* message = nullptr;
    const wirehair_v2::PrecodeEncoder* block = nullptr;
    wirehair_v2::PacketRowConfig config = {};
    wirehair_v2::PacketRowRuntime runtime = {};
};

bool CostRangeFits(const void* pointer, size_t bytes)
{
    return pointer && bytes != 0u && bytes <=
        std::numeric_limits<uintptr_t>::max() -
        reinterpret_cast<uintptr_t>(pointer);
}

bool CostWritable(const CostInspection& state, const void* pointer, size_t bytes)
{
    return CostRangeFits(pointer, bytes) &&
        !MemoryRangesOverlap(pointer, bytes, state.pod.source,
            static_cast<size_t>(state.pod.message_bytes)) &&
        !MemoryRangesOverlap(pointer, bytes, state.pod.intermediate,
            static_cast<size_t>(state.pod.intermediate_capacity)) &&
        !MemoryRangesOverlap(pointer, bytes, state.facade, sizeof(*state.facade)) &&
        !MemoryRangesOverlap(pointer, bytes, state.message, sizeof(*state.message)) &&
        !MemoryRangesOverlap(pointer, bytes, state.block, sizeof(*state.block));
}

WirehairV2Result CostCapture(WirehairV2Codec handle, CostInspection& state)
{
    const PublicCodec* const facade = FromHandle(handle);
    if (!facade || facade->Mode != CodecMode::Encoder || facade->Decoded ||
        facade->SourceState != EncoderSourceState::BorrowedImmutable ||
        !facade->BorrowedSource || facade->MessageBytes == 0u ||
        !CompleteMessageRangeRepresentable(facade->BorrowedSource,
            facade->MessageBytes))
    {
        return WirehairV2_InvalidInput;
    }
    const wirehair_v2::MessagePrecodeEncoder* const message =
        (facade->Impl.*CostMember(CostPrecodeTag{})).get();
    if (!message || !message->IsInitialized() ||
        message->MessageBytes() != facade->MessageBytes ||
        message->BlockBytes() != facade->BlockBytes)
    {
        return WirehairV2_InvalidInput;
    }
    const wirehair_v2::PrecodeEncoder& block = message->BlockEncoder();
    const CostStorageTag::Storage& storage = block.*CostMember(CostStorageTag{});
    const uint32_t K = block.SourceBlockCount();
    const uint32_t P = block.ParityBlockCount();
    const uint32_t B = block.BlockBytes();
    const uint64_t columns = static_cast<uint64_t>(K) + P;
    if (!block.IsInitialized() || block.HasCompleteSystem() ||
        !K || !P || !B || B > UINT32_C(0x7fffffff) ||
        K != message->SourceBlockCount() || B != facade->BlockBytes ||
        BlockCountWide(facade->MessageBytes, B) != K ||
        columns > static_cast<uint64_t>(std::numeric_limits<size_t>::max()) / B ||
        storage.size() != static_cast<size_t>(columns) * B ||
        storage.capacity() < storage.size() ||
        storage.data() != block.IntermediateBlocks() ||
        message->IntermediateBlocks() != block.IntermediateBlocks() ||
        !CostRangeFits(storage.data(), storage.capacity()) ||
        block.ParityBlocks() != storage.data() + static_cast<size_t>(K) * B ||
        block.RecoveryRowSeed() > UINT32_MAX)
    {
        return WirehairV2_InvalidInput;
    }
#if WH2_ALIGNMENT_INSPECT_ALIGNED
    if (reinterpret_cast<uintptr_t>(storage.data()) % 32u != 0u) {
        return WirehairV2_Error;
    }
#endif
    const wirehair_v2::PrecodeParams& params = block.System().Params;
    const wirehair_v2::SeedProfile& profile = message->Profile();
    state.config.PeelSeed = static_cast<uint32_t>(block.RecoveryRowSeed());
    state.config.MixCount = block.RecoveryMixCount();
    if (profile.BlockBytes != B || params.BlockCount != profile.BlockCount ||
        params.BlockCount != K ||
        params.Staircase != profile.V2StaircaseCount ||
        params.DenseRows != profile.V2DenseRowCount ||
        params.HeavyRows != profile.V2HeavyRowCount ||
        static_cast<uint64_t>(params.Staircase) + params.DenseRows +
            params.HeavyRows != P ||
        params.SourceHits != profile.V2SourceHits ||
        params.DenseIdentityCorner != profile.V2DenseIdentityCorner ||
        params.HeavyFamily != wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
        static_cast<uint32_t>(params.DenseAnchors) != profile.V2DenseAnchorLayout ||
        params.Seed != profile.V2PrecodeSeed ||
        state.config.PeelSeed != profile.V2PacketPeelSeed ||
        state.config.MixCount != profile.V2RecoveryMixCount ||
        !state.runtime.Initialize(K, P, state.config.MixCount))
    {
        return WirehairV2_InvalidInput;
    }
    WirehairV2Profile inner = {}, outer = {};
    WirehairV2Result result = MakePublicProfile(profile, facade->MessageBytes,
        WIREHAIR_V2_PROFILE_CERTIFIED_2026_07, inner);
    if (result != WirehairV2_Success) { return result; }
    result = MakePublicProfile(facade->Impl.Profile(), facade->MessageBytes,
        WIREHAIR_V2_PROFILE_CERTIFIED_2026_07, outer);
    if (result != WirehairV2_Success) { return result; }
    uint8_t outer_bytes[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t written = 0;
    result = wirehair_v2_profile_serialize(&inner, state.pod.serialized_profile,
        sizeof(state.pod.serialized_profile), &written);
    if (result != WirehairV2_Success || written != sizeof(state.pod.serialized_profile)) {
        return WirehairV2_Error;
    }
    result = wirehair_v2_profile_serialize(&outer, outer_bytes,
        sizeof(outer_bytes), &written);
    if (result != WirehairV2_Success || written != sizeof(outer_bytes) ||
        std::memcmp(outer_bytes, state.pod.serialized_profile, sizeof(outer_bytes)) != 0)
    {
        return WirehairV2_Error;
    }
    state.pod.handle = handle;
    state.pod.source = facade->BorrowedSource;
    state.pod.intermediate = storage.data();
    state.pod.message_bytes = facade->MessageBytes;
    state.pod.intermediate_bytes = storage.size();
    state.pod.intermediate_capacity = storage.capacity();
    state.pod.block_bytes = B;
    state.pod.source_count = K;
    state.pod.precode_count = P;
    state.pod.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    state.facade = facade;
    state.message = message;
    state.block = &block;
    return WirehairV2_Success;
}

std::vector<uint32_t> CostColumns(const CostInspection& state, uint32_t id)
{
    return wirehair_v2::GeneratePacketMatrixRowWithRuntime(
        state.pod.source_count, state.pod.precode_count, id,
        state.config, state.runtime);
}

bool CostColumnsValid(const CostInspection& state, const std::vector<uint32_t>& row)
{
    const uint64_t columns = static_cast<uint64_t>(state.pod.source_count) +
        state.pod.precode_count;
    if (row.empty() || row.size() > columns || row.size() > UINT32_MAX) {
        return false;
    }
    for (uint32_t column : row) {
        if (column >= columns) { return false; }
    }
    return true;
}

} // namespace

extern "C" WirehairV2Result WH2_ALIGNMENT_INSPECT_FN(
    WirehairV2Codec encoder, Wh2AlignedIntermediateCostSnapshot* output)
{
    if (!output) { return WirehairV2_InvalidInput; }
    try {
        CostInspection state;
        const WirehairV2Result result = CostCapture(encoder, state);
        if (result != WirehairV2_Success) { return result; }
        if (!CostWritable(state, output, sizeof(*output))) {
            return WirehairV2_InvalidInput;
        }
        *output = state.pod;
        return WirehairV2_Success;
    }
    catch (const std::bad_alloc&) { return WirehairV2_OOM; }
    catch (const std::length_error&) { return WirehairV2_OOM; }
}

extern "C" WirehairV2Result WH2_ALIGNMENT_COLUMNS_FN(
    WirehairV2Codec encoder, uint32_t id, uint32_t* columns,
    uint32_t capacity, uint32_t* count_out)
{
    if (!count_out) { return WirehairV2_InvalidInput; }
    try {
        CostInspection state;
        const WirehairV2Result result = CostCapture(encoder, state);
        if (result != WirehairV2_Success) { return result; }
        const std::vector<uint32_t> row = CostColumns(state, id);
        if (!CostColumnsValid(state, row)) { return WirehairV2_Error; }
        const size_t bytes = row.size() * sizeof(uint32_t);
        if (!CostWritable(state, count_out, sizeof(*count_out)) ||
            (columns && (!CostWritable(state, columns, bytes) ||
             MemoryRangesOverlap(columns, bytes, count_out, sizeof(*count_out)))))
        {
            return WirehairV2_InvalidInput;
        }
        if (capacity < row.size()) {
            *count_out = static_cast<uint32_t>(row.size());
            return WirehairV2_BufferTooSmall;
        }
        if (!columns) { return WirehairV2_InvalidInput; }
        std::memcpy(columns, row.data(), bytes);
        *count_out = static_cast<uint32_t>(row.size());
        return WirehairV2_Success;
    }
    catch (const std::bad_alloc&) { return WirehairV2_OOM; }
    catch (const std::length_error&) { return WirehairV2_OOM; }
}

extern "C" WirehairV2Result WH2_ALIGNMENT_SCALAR_FN(
    WirehairV2Codec encoder, uint32_t id, void* output,
    uint32_t capacity, uint32_t* bytes_out)
{
    if (!bytes_out) { return WirehairV2_InvalidInput; }
    try {
        CostInspection state;
        const WirehairV2Result result = CostCapture(encoder, state);
        if (result != WirehairV2_Success) { return result; }
        uint32_t required = state.pod.block_bytes;
        if (id == state.pod.source_count - 1u) {
            required = static_cast<uint32_t>(state.pod.message_bytes -
                static_cast<uint64_t>(id) * state.pod.block_bytes);
        }
        if (!CostWritable(state, bytes_out, sizeof(*bytes_out)) ||
            (output && (!CostWritable(state, output, required) ||
             MemoryRangesOverlap(output, required, bytes_out, sizeof(*bytes_out)))))
        {
            return WirehairV2_InvalidInput;
        }
        if (capacity < required) {
            *bytes_out = required;
            return WirehairV2_BufferTooSmall;
        }
        if (!output) { return WirehairV2_InvalidInput; }
        const std::vector<uint32_t> row = CostColumns(state, id);
        if (!CostColumnsValid(state, row)) { return WirehairV2_Error; }
        uint8_t* const destination = static_cast<uint8_t*>(output);
        // Independent scalar XOR, not a production evaluator/GF call and not
        // a borrowed-systematic shortcut.  All possible failures precede writes.
        for (uint32_t j = 0; j < required; ++j) {
            uint8_t value = 0;
            for (uint32_t column : row) {
                value ^= state.pod.intermediate[
                    static_cast<size_t>(column) * state.pod.block_bytes + j];
            }
            destination[j] = value;
        }
        *bytes_out = required;
        return WirehairV2_Success;
    }
    catch (const std::bad_alloc&) { return WirehairV2_OOM; }
    catch (const std::length_error&) { return WirehairV2_OOM; }
}
