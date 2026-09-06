// This translation unit owns the ONLY V2 facade in the diagnostic executable.
// Link the original archive, but never extract its WirehairV2Profile.cpp.o.
// These are the immutable archive facade's flags, not ordinary bench flags.
#if !defined(WIREHAIR_BUILDING) || WIREHAIR_BUILDING != 1 || \
    !defined(WIREHAIR_STATIC) || WIREHAIR_STATIC != 1 || !defined(__PIC__)
#error "alignment bridge requires original BUILDING=1 STATIC=1 PIC facade flags"
#endif
#if defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WH_SEED_KNOBS)
#error "alignment bridge requires the original production equation/runtime domain"
#endif

#include "../codec/WirehairV2Profile.cpp"
#include "Wh2RepairAlignmentR0Bridge.h"

#include <utility>

namespace wh2_repair_alignment_r0 {
namespace {

// Explicit-instantiation access exemption; no altered class definitions,
// layout guesses, untyped offsets, ownership transfer or private/public macro.
struct PrecodeTag
{
    using Type = std::unique_ptr<wirehair_v2::MessagePrecodeEncoder>
        wirehair_v2::Codec::*;
    friend Type Member(PrecodeTag);
};

template<class Tag, typename Tag::Type pointer>
struct Access
{
    friend typename Tag::Type Member(Tag) { return pointer; }
};

template struct Access<PrecodeTag, &wirehair_v2::Codec::PrecodeImpl>;

bool Extent(uint32_t K, uint32_t P, uint32_t B, size_t& bytes)
{
    const uint64_t columns = (uint64_t)K + P;
    if (!K || !P || !B || B > UINT32_C(0x7fffffff) ||
        columns > (uint64_t)std::numeric_limits<size_t>::max() / B)
    {
        return false;
    }
    bytes = (size_t)columns * B;
    return true;
}

bool RangeFits(const void* pointer, size_t bytes)
{
    return pointer && bytes && bytes <=
        std::numeric_limits<uintptr_t>::max() -
        reinterpret_cast<uintptr_t>(pointer);
}

bool ParamsMatch(
    const wirehair_v2::PrecodeParams& params,
    const wirehair_v2::SeedProfile& profile,
    const wirehair_v2::PacketRowConfig& config)
{
    return params.BlockCount == profile.BlockCount &&
        params.Staircase == profile.V2StaircaseCount &&
        params.DenseRows == profile.V2DenseRowCount &&
        params.HeavyRows == profile.V2HeavyRowCount &&
        params.SourceHits == profile.V2SourceHits &&
        params.DenseIdentityCorner == profile.V2DenseIdentityCorner &&
        params.HeavyFamily == wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy &&
        (uint32_t)params.DenseAnchors == profile.V2DenseAnchorLayout &&
        params.Seed == profile.V2PrecodeSeed &&
        config.PeelSeed == profile.V2PacketPeelSeed &&
        config.MixCount == profile.V2RecoveryMixCount;
}

bool FixedDimensions(const Snapshot& snapshot)
{
    size_t bytes = 0;
    return snapshot.message_bytes == 7680u && snapshot.block_bytes == 1280u &&
        snapshot.source_count == 6u && snapshot.precode_count == 30u &&
        Extent(snapshot.source_count, snapshot.precode_count,
            snapshot.block_bytes, bytes) &&
        snapshot.intermediate_bytes == bytes && bytes == 46080u;
}

} // namespace

WirehairV2Result Capture(WirehairV2Codec live_encoder, Snapshot& output)
{
    const PublicCodec* const facade = FromHandle(live_encoder);
    if (!facade || facade->Mode != CodecMode::Encoder || facade->Decoded ||
        facade->SourceState != EncoderSourceState::BorrowedImmutable ||
        !facade->BorrowedSource || facade->MessageBytes != 7680u ||
        facade->BlockBytes != 1280u)
    {
        return WirehairV2_InvalidInput;
    }
    const wirehair_v2::MessagePrecodeEncoder* const message =
        (facade->Impl.*Member(PrecodeTag{})).get();
    if (!message || !message->IsInitialized() ||
        message->MessageBytes() != facade->MessageBytes ||
        message->BlockBytes() != facade->BlockBytes ||
        message->SourceBlockCount() != 6u)
    {
        return WirehairV2_InvalidInput;
    }
    const wirehair_v2::PrecodeEncoder& block = message->BlockEncoder();
    if (!block.IsInitialized() || block.HasCompleteSystem() ||
        !block.IntermediateBlocks() ||
        message->IntermediateBlocks() != block.IntermediateBlocks() ||
        block.SourceBlockCount() != 6u || block.ParityBlockCount() != 30u ||
        block.BlockBytes() != 1280u || block.RecoveryRowSeed() > UINT32_MAX ||
        block.RecoveryMixCount() != 3u)
    {
        return WirehairV2_InvalidInput;
    }

    try
    {
        Snapshot next;
        next.handle = live_encoder;
        next.source = facade->BorrowedSource;
        next.intermediate = block.IntermediateBlocks();
        next.message_bytes = facade->MessageBytes;
        next.block_bytes = block.BlockBytes();
        next.source_count = block.SourceBlockCount();
        next.precode_count = block.ParityBlockCount();
        next.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
        next.expanded_profile = message->Profile();
        next.system.Params = block.System().Params;
        next.config.PeelSeed = (uint32_t)block.RecoveryRowSeed();
        next.config.MixCount = block.RecoveryMixCount();
        if (!Extent(next.source_count, next.precode_count, next.block_bytes,
                next.intermediate_bytes) || !FixedDimensions(next) ||
            !RangeFits(next.source, (size_t)next.message_bytes) ||
            !RangeFits(next.intermediate, next.intermediate_bytes) ||
            !ParamsMatch(next.system.Params, next.expanded_profile, next.config) ||
            !next.runtime.Initialize(next.source_count, next.precode_count,
                next.config.MixCount))
        {
            return WirehairV2_InvalidInput;
        }

        // Validate both original profile copies using the unchanged facade's
        // canonical expansion, not only a compact descriptor/seed-attempt test.
        WirehairV2Result result = MakePublicProfile(next.expanded_profile,
            next.message_bytes, WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            next.profile);
        if (result != WirehairV2_Success) {
            return result;
        }
        WirehairV2Profile outer_profile = {};
        result = MakePublicProfile(facade->Impl.Profile(), next.message_bytes,
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07, outer_profile);
        if (result != WirehairV2_Success) {
            return result;
        }
        uint32_t bytes = 0;
        result = wirehair_v2_profile_serialize(&next.profile,
            next.serialized_profile, sizeof(next.serialized_profile), &bytes);
        if (result != WirehairV2_Success || bytes != sizeof(next.serialized_profile)) {
            return WirehairV2_Error;
        }
        uint8_t outer_serialized[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
        result = wirehair_v2_profile_serialize(&outer_profile,
            outer_serialized, sizeof(outer_serialized), &bytes);
        if (result != WirehairV2_Success || bytes != sizeof(outer_serialized) ||
            std::memcmp(next.serialized_profile, outer_serialized,
                sizeof(outer_serialized)) != 0)
        {
            return WirehairV2_Error;
        }
        output = std::move(next);
        return WirehairV2_Success;
    }
    catch (const std::bad_alloc&) {
        return WirehairV2_OOM;
    }
    catch (const std::length_error&) {
        return WirehairV2_OOM;
    }
}

bool ValidateShadowBuffers(
    const Snapshot& snapshot,
    const uint8_t* intermediate,
    size_t intermediate_capacity,
    uint8_t* output,
    size_t output_capacity)
{
    return FixedDimensions(snapshot) &&
        intermediate_capacity >= snapshot.intermediate_bytes &&
        output_capacity >= snapshot.block_bytes &&
        RangeFits(intermediate, snapshot.intermediate_bytes) &&
        RangeFits(output, snapshot.block_bytes) &&
        RangeFits(snapshot.intermediate, snapshot.intermediate_bytes) &&
        RangeFits(snapshot.source, (size_t)snapshot.message_bytes) &&
        !MemoryRangesOverlap(output, snapshot.block_bytes,
            intermediate, snapshot.intermediate_bytes) &&
        !MemoryRangesOverlap(output, snapshot.block_bytes,
            snapshot.intermediate, snapshot.intermediate_bytes) &&
        !MemoryRangesOverlap(output, snapshot.block_bytes,
            snapshot.source, (size_t)snapshot.message_bytes);
}

Evaluation EvaluateShadow(
    const Snapshot& snapshot,
    const uint8_t* intermediate,
    uint32_t block_id,
    uint8_t* output)
{
    Evaluation result;
    if (!intermediate || !output || block_id < 6u || block_id > 11u) {
        return result;
    }
    try
    {
        // Mirror the real PrecodeEncoder::EncodeResult operation-scratch branch.
        uint64_t local_ops = 0;
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
                snapshot.system, snapshot.config, snapshot.runtime, intermediate,
                snapshot.block_bytes, block_id, output, &local_ops))
        {
            return result;
        }
        result.status = WirehairV2_Success;
        result.bytes = snapshot.block_bytes;
        result.operations = local_ops;
        return result;
    }
    catch (const std::bad_alloc&) {
        result.status = WirehairV2_OOM;
    }
    catch (const std::length_error&) {
        result.status = WirehairV2_OOM;
    }
    return result;
}

bool NeutralSelfTest()
{
    size_t bytes = 123;
    if (!Extent(6, 30, 1280, bytes) || bytes != 46080 ||
        Extent(0, 30, 1280, bytes) || Extent(6, 0, 1280, bytes) ||
        Extent(6, 30, 0, bytes) || Extent(6, 30, UINT32_MAX, bytes))
    {
        return false;
    }
    Snapshot neutral;
    neutral.message_bytes = 7680;
    neutral.block_bytes = 1280;
    neutral.source_count = 6;
    neutral.precode_count = 30;
    neutral.intermediate_bytes = 46080;
    // Integer addresses are never dereferenced; no fixture/codec is evaluated.
    neutral.source = reinterpret_cast<const uint8_t*>(uintptr_t(0x10000));
    neutral.intermediate = reinterpret_cast<const uint8_t*>(uintptr_t(0x20000));
    const uint8_t* const view =
        reinterpret_cast<const uint8_t*>(uintptr_t(0x40000));
    uint8_t* const out = reinterpret_cast<uint8_t*>(uintptr_t(0x60000));
    if (!ValidateShadowBuffers(neutral, view, 46080, out, 1280) ||
        !ValidateShadowBuffers(neutral, neutral.intermediate, 46080, out, 1280) ||
        !ValidateShadowBuffers(neutral, view, 46080,
            reinterpret_cast<uint8_t*>(uintptr_t(0x40000 + 46080)), 1280) ||
        !ValidateShadowBuffers(neutral, view, 46080,
            reinterpret_cast<uint8_t*>(uintptr_t(0x40000 - 1280)), 1280) ||
        ValidateShadowBuffers(neutral, view, 46079, out, 1280) ||
        ValidateShadowBuffers(neutral, view, 46080, out, 1279) ||
        ValidateShadowBuffers(neutral, nullptr, 46080, out, 1280) ||
        ValidateShadowBuffers(neutral, view, 46080, nullptr, 1280) ||
        ValidateShadowBuffers(neutral, view, 46080,
            const_cast<uint8_t*>(view), 1280) ||
        ValidateShadowBuffers(neutral, view, 46080,
            reinterpret_cast<uint8_t*>(uintptr_t(0x40000 + 46079)), 1280) ||
        ValidateShadowBuffers(neutral, view, 46080,
            reinterpret_cast<uint8_t*>(uintptr_t(0x40000 - 1279)), 1280) ||
        ValidateShadowBuffers(neutral, view, 46080,
            const_cast<uint8_t*>(neutral.source), 1280) ||
        ValidateShadowBuffers(neutral, view, 46080,
            const_cast<uint8_t*>(neutral.intermediate), 1280) ||
        RangeFits(reinterpret_cast<const void*>(
            std::numeric_limits<uintptr_t>::max() - 1279u), 1280u))
    {
        return false;
    }
    neutral.intermediate_bytes = 46079;
    if (ValidateShadowBuffers(neutral, view, 46080, out, 1280)) {
        return false;
    }
    wirehair_v2::PrecodeParams params;
    wirehair_v2::SeedProfile profile;
    wirehair_v2::PacketRowConfig config;
    profile.V2RecoveryMixCount = config.MixCount;
    if (!ParamsMatch(params, profile, config)) {
        return false;
    }
    params.Seed = 1;
    if (ParamsMatch(params, profile, config)) {
        return false;
    }
    params.Seed = 0;
    config.PeelSeed = 1;
    return !ParamsMatch(params, profile, config);
}

} // namespace wh2_repair_alignment_r0
