#include "WirehairV2Codec.h"

#include "../WirehairTools.h"

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
    }
    return result;
}

WirehairResult Codec::InitializeDecoder(
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override)
{
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
    if (!block_in || block_bytes == 0u) {
        return Wirehair_InvalidInput;
    }
    return Impl.DecodeFeed(block_id, block_in, block_bytes);
}

WirehairResult Codec::Recover(void* message_out, uint64_t message_bytes)
{
    if (!message_out) {
        return Wirehair_InvalidInput;
    }
    return Impl.ReconstructOutput(message_out, message_bytes);
}

const SeedProfile& Codec::Profile() const
{
    return CurrentProfile;
}

} // namespace wirehair_v2
