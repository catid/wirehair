#include "WirehairV2Codec.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <new>
#include <utility>

namespace wirehair_v2 {
namespace {

uint64_t BlockCountWideFor(uint64_t message_bytes, uint32_t block_bytes)
{
    // Divide before adding the round-up so the sum cannot wrap uint64
    return message_bytes / block_bytes +
        ((message_bytes % block_bytes != 0u) ? 1u : 0u);
}

uint32_t BlockCountFor(uint64_t message_bytes, uint32_t block_bytes)
{
    // Only valid after ValidateCodecInput has bounded the wide count
    return (uint32_t)BlockCountWideFor(message_bytes, block_bytes);
}

WirehairResult ValidateCodecInput(uint64_t message_bytes, uint32_t block_bytes)
{
    if (message_bytes < 1u || block_bytes < 1u) {
        return Wirehair_InvalidInput;
    }
    // Validate at full width before narrowing: truncating first would wrap
    // huge messages into the accepted range or report the wrong error code
    const uint64_t block_count = BlockCountWideFor(message_bytes, block_bytes);
    if (block_count < CAT_WIREHAIR_MIN_N) {
        return Wirehair_BadInput_SmallN;
    }
    if (block_count > CAT_WIREHAIR_MAX_N) {
        return Wirehair_BadInput_LargeN;
    }
    return Wirehair_Success;
}

WirehairResult EnsureGf256Initialized()
{
    return gf256_init() == 0 ? Wirehair_Success : Wirehair_UnsupportedPlatform;
}

} // namespace

Codec::Codec()
{
}

Codec::~Codec()
{
}

WirehairResult Codec::InitializeEncoder(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override)
{
    if (!message) {
        return Wirehair_InvalidInput;
    }
    const WirehairResult gf_result = EnsureGf256Initialized();
    if (gf_result != Wirehair_Success) {
        return gf_result;
    }
    const WirehairResult input_result =
        ValidateCodecInput(message_bytes, block_bytes);
    if (input_result != Wirehair_Success) {
        return input_result;
    }

    const uint32_t block_count = BlockCountFor(message_bytes, block_bytes);
    if (seed_override &&
        (seed_override->BlockCount != block_count ||
         seed_override->BlockBytes != block_bytes))
    {
        return Wirehair_InvalidInput;
    }
    const SeedProfile profile = seed_override ? *seed_override :
        SelectSeedProfile(block_count, block_bytes);

    Impl.OverrideSeeds(
        profile.DenseCount,
        profile.PeelSeed,
        profile.DenseSeed);

    WirehairResult result = Impl.InitializeEncoder(message_bytes, block_bytes);
    if (result == Wirehair_Success) {
        result = Impl.EncodeFeed(message);
    }
    if (result == Wirehair_Success) {
        // Only publish the profile for initializations that succeeded
        CurrentProfile = profile;
        PrecodeImpl.reset();
        PrecodeDecoderImpl.reset();
    }
    return result;
}

WirehairResult Codec::InitializePrecodeEncoder(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override,
    const MessagePrecodeEncoderOptions* options)
{
    if (!message) {
        return Wirehair_InvalidInput;
    }
    const WirehairResult gf_result = EnsureGf256Initialized();
    if (gf_result != Wirehair_Success) {
        return gf_result;
    }
    const WirehairResult input_result =
        ValidateCodecInput(message_bytes, block_bytes);
    if (input_result != Wirehair_Success) {
        return input_result;
    }

    const uint32_t block_count = BlockCountFor(message_bytes, block_bytes);
    if (seed_override &&
        (seed_override->BlockCount != block_count ||
         seed_override->BlockBytes != block_bytes))
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        std::unique_ptr<MessagePrecodeEncoder> next(
            new MessagePrecodeEncoder());
        const WirehairResult result = next->InitializeResult(
            message, message_bytes, block_bytes, seed_override, options);
        if (result != Wirehair_Success) {
            return result;
        }

        CurrentProfile = next->Profile();
        PrecodeImpl = std::move(next);
        PrecodeDecoderImpl.reset();
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
}

WirehairResult Codec::InitializePrecodeDecoder(
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override,
    const MessagePrecodeEncoderOptions* options)
{
    const WirehairResult gf_result = EnsureGf256Initialized();
    if (gf_result != Wirehair_Success) {
        return gf_result;
    }
    const WirehairResult input_result =
        ValidateCodecInput(message_bytes, block_bytes);
    if (input_result != Wirehair_Success) {
        return input_result;
    }

    const uint32_t block_count = BlockCountFor(message_bytes, block_bytes);
    if (seed_override &&
        (seed_override->BlockCount != block_count ||
         seed_override->BlockBytes != block_bytes))
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        std::unique_ptr<MessagePrecodeDecoder> next(
            new MessagePrecodeDecoder());
        const WirehairResult result = next->InitializeResult(
            message_bytes, block_bytes, seed_override, options);
        if (result != Wirehair_Success) {
            return result;
        }

        CurrentProfile = next->Profile();
        PrecodeDecoderImpl = std::move(next);
        PrecodeImpl.reset();
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
}

WirehairResult Codec::InitializeDecoder(
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override)
{
    const WirehairResult gf_result = EnsureGf256Initialized();
    if (gf_result != Wirehair_Success) {
        return gf_result;
    }
    const WirehairResult input_result =
        ValidateCodecInput(message_bytes, block_bytes);
    if (input_result != Wirehair_Success) {
        return input_result;
    }

    const uint32_t block_count = BlockCountFor(message_bytes, block_bytes);
    if (seed_override &&
        (seed_override->BlockCount != block_count ||
         seed_override->BlockBytes != block_bytes))
    {
        return Wirehair_InvalidInput;
    }
    const SeedProfile profile = seed_override ? *seed_override :
        SelectSeedProfile(block_count, block_bytes);

    Impl.OverrideSeeds(
        profile.DenseCount,
        profile.PeelSeed,
        profile.DenseSeed);

    const WirehairResult result =
        Impl.InitializeDecoder(message_bytes, block_bytes);
    if (result == Wirehair_Success) {
        // Only publish the profile for initializations that succeeded
        CurrentProfile = profile;
        PrecodeImpl.reset();
        PrecodeDecoderImpl.reset();
    }
    return result;
}

WirehairResult Codec::Encode(
    uint32_t block_id,
    void* block_out,
    uint32_t out_bytes,
    uint32_t* data_bytes_out)
{
    if (!block_out || !data_bytes_out) {
        return Wirehair_InvalidInput;
    }
    if (PrecodeDecoderImpl && PrecodeDecoderImpl->IsInitialized()) {
        return Wirehair_InvalidInput;
    }
    if (PrecodeImpl && PrecodeImpl->IsInitialized())
    {
        return PrecodeImpl->EncodeResult(
            block_id,
            static_cast<uint8_t*>(block_out),
            out_bytes,
            data_bytes_out);
    }
    const uint32_t written = Impl.Encode(block_id, block_out, out_bytes);
    *data_bytes_out = written;
    if (written == 0u) {
        return Wirehair_InvalidInput;
    }
    return Wirehair_Success;
}

WirehairResult Codec::Decode(
    uint32_t block_id,
    const void* block_in,
    uint32_t block_bytes)
{
    if (PrecodeImpl && PrecodeImpl->IsInitialized()) {
        return Wirehair_InvalidInput;
    }
    if (!block_in || block_bytes == 0u) {
        return Wirehair_InvalidInput;
    }
    if (PrecodeDecoderImpl && PrecodeDecoderImpl->IsInitialized()) {
        return PrecodeDecoderImpl->DecodeResult(
            block_id,
            static_cast<const uint8_t*>(block_in),
            block_bytes);
    }
    return Impl.DecodeFeed(block_id, block_in, block_bytes);
}

WirehairResult Codec::Recover(void* message_out, uint64_t message_bytes)
{
    if (PrecodeImpl && PrecodeImpl->IsInitialized()) {
        return Wirehair_InvalidInput;
    }
    if (!message_out) {
        return Wirehair_InvalidInput;
    }
    if (PrecodeDecoderImpl && PrecodeDecoderImpl->IsInitialized()) {
        return PrecodeDecoderImpl->RecoverResult(message_out, message_bytes);
    }
    return Impl.ReconstructOutput(message_out, message_bytes);
}

const SeedProfile& Codec::Profile() const
{
    return CurrentProfile;
}

} // namespace wirehair_v2
