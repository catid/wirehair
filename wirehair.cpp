/*
    Copyright (c) 2012-2018 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Wirehair nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include <wirehair/wirehair.h>
#include "WirehairCodec.h"
#include "WirehairEnvironment.h"

#include <atomic>
#include <cstdlib>
#include <mutex>
#include <new> // std::nothrow

namespace {

std::once_flag InitOnce;
std::atomic<int> PublishedInitResult(Wirehair_Error);
WirehairResult InitResult = Wirehair_Error;

static_assert(sizeof(WirehairWireProfile) == 16,
    "WirehairWireProfile is a fixed C ABI descriptor");

bool ResolveWireProfileId(uint64_t profile_id, wirehair::WireProfile& profile)
{
    switch (profile_id)
    {
    case WIREHAIR_LEGACY_PROFILE_PRE_FIXUP:
        profile = wirehair::WireProfile::LegacyPreFixup;
        return true;
    case WIREHAIR_LEGACY_PROFILE_CURRENT:
        profile = wirehair::WireProfile::LegacyCurrent;
        return true;
    default:
        return false;
    }
}

bool ResolveWireProfileDescriptor(
    const WirehairWireProfile* descriptor,
    wirehair::WireProfile& profile)
{
    return descriptor &&
        descriptor->struct_bytes == sizeof(WirehairWireProfile) &&
        descriptor->profile_version == WIREHAIR_WIRE_PROFILE_VERSION &&
        ResolveWireProfileId(descriptor->profile_id, profile);
}

bool NamedLegacyWireProfilesAvailable()
{
    // These private experiment knobs alter packet equations.  Raw APIs remain
    // available to experiment binaries, but such binaries must not claim one
    // of the immutable named profiles through the explicit descriptor path.
#if defined(WH_SEED_KNOBS) || defined(CAT_IDENTITY_LOWER_RIGHT) || \
    WH_HEAVY_ROWS != 6 || WH_HEAVY_COLS != 18
    return false;
#else
    return true;
#endif
}

WirehairResult GetPublishedInitResult()
{
    return static_cast<WirehairResult>(
        PublishedInitResult.load(std::memory_order_acquire));
}

bool TestForcesAllocationFailure()
{
#if defined(WIREHAIR_TESTING)
    const wirehair::EnvironmentValue environment(
        "WIREHAIR_TEST_FORCE_OOM");
    const char* value = environment.Get();
    return value && value[0] == '1' && value[1] == '\0';
#else
    return false;
#endif
}

WirehairResult TestForcedCreationResult()
{
#if defined(WIREHAIR_TESTING)
    const wirehair::EnvironmentValue environment(
        "WIREHAIR_TEST_FORCE_CREATE_RESULT");
    const char* value = environment.Get();
    if (value && *value) {
        const int result = std::atoi(value);
        if (result > Wirehair_Success && result < WirehairResult_Count) {
            return static_cast<WirehairResult>(result);
        }
    }
#endif
    return Wirehair_Success;
}

WirehairResult CreateEncoder(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    bool copyInput,
    wirehair::WireProfile profile,
    WirehairCodec* codecOut)
{
    if (!codecOut) {
        return Wirehair_InvalidInput;
    }
    *codecOut = nullptr;

    wirehair::Codec* codec = reinterpret_cast<wirehair::Codec*>(reuseOpt);
    const WirehairResult initResult = GetPublishedInitResult();
    if (initResult != Wirehair_Success) {
        delete codec;
        return initResult;
    }
    if (!message || messageBytes < 1 || blockBytes < 1) {
        delete codec;
        return Wirehair_InvalidInput;
    }
    if (TestForcesAllocationFailure()) {
        delete codec;
        return Wirehair_OOM;
    }
    const WirehairResult forcedResult = TestForcedCreationResult();
    if (forcedResult != Wirehair_Success) {
        delete codec;
        return forcedResult;
    }

    if (!codec) {
        codec = new (std::nothrow) wirehair::Codec;
        if (!codec) {
            return Wirehair_OOM;
        }
    }

    codec->SetWireProfile(profile);
    WirehairResult result = codec->InitializeEncoder(messageBytes, blockBytes);
    if (result == Wirehair_Success) {
        result = codec->EncodeFeed(message, copyInput);
    }
    if (result != Wirehair_Success) {
        delete codec;
        return result;
    }

    *codecOut = reinterpret_cast<WirehairCodec>(codec);
    return Wirehair_Success;
}

WirehairResult CreateDecoder(
    WirehairCodec reuseOpt,
    uint64_t messageBytes,
    uint32_t blockBytes,
    wirehair::WireProfile profile,
    WirehairCodec* codecOut)
{
    if (!codecOut) {
        return Wirehair_InvalidInput;
    }
    *codecOut = nullptr;

    wirehair::Codec* codec = reinterpret_cast<wirehair::Codec*>(reuseOpt);
    const WirehairResult initResult = GetPublishedInitResult();
    if (initResult != Wirehair_Success) {
        delete codec;
        return initResult;
    }
    if (messageBytes < 1 || blockBytes < 1) {
        delete codec;
        return Wirehair_InvalidInput;
    }
    if (TestForcesAllocationFailure()) {
        delete codec;
        return Wirehair_OOM;
    }
    const WirehairResult forcedResult = TestForcedCreationResult();
    if (forcedResult != Wirehair_Success) {
        delete codec;
        return forcedResult;
    }

    if (!codec) {
        codec = new (std::nothrow) wirehair::Codec;
        if (!codec) {
            return Wirehair_OOM;
        }
    }

    codec->SetWireProfile(profile);
    const WirehairResult result =
        codec->InitializeDecoder(messageBytes, blockBytes);
    if (result != Wirehair_Success) {
        delete codec;
        return result;
    }

    *codecOut = reinterpret_cast<WirehairCodec>(codec);
    return Wirehair_Success;
}

} // namespace


extern "C" {


//-----------------------------------------------------------------------------
// Wirehair API

WIREHAIR_EXPORT const char *wirehair_result_string(
    WirehairResult result ///< Result code to convert to string
)
{
    static_assert(WirehairResult_Count == 11, "Update this switch too");

    switch (result)
    {
    case Wirehair_Success:           return "Wirehair_Success";
    case Wirehair_NeedMore:    return "Wirehair_NeedMore";
    case Wirehair_BadDenseSeed:      return "Wirehair_BadDenseSeed";
    case Wirehair_BadPeelSeed:       return "Wirehair_BadPeelSeed";
    case Wirehair_BadInput_SmallN:   return "Wirehair_BadInput_SmallN";
    case Wirehair_BadInput_LargeN:   return "Wirehair_BadInput_LargeN";
    case Wirehair_ExtraInsufficient: return "Wirehair_ExtraInsufficient";
    case Wirehair_InvalidInput:      return "Wirehair_InvalidInput";
    case Wirehair_OOM:               return "Wirehair_OOM";
    case Wirehair_UnsupportedPlatform: return "Wirehair_UnsupportedPlatform";
    case Wirehair_Error:             return "Wirehair_Error";
    default:
        break;
    }

    return "Unknown";
}

WIREHAIR_EXPORT WirehairResult wirehair_init_(int expected_version)
{
    // If version does not match:
    if (expected_version != WIREHAIR_VERSION) {
        return Wirehair_InvalidInput;
    }

    std::call_once(InitOnce, []() {
        InitResult = gf256_init() == 0 ?
            Wirehair_Success : Wirehair_UnsupportedPlatform;
        PublishedInitResult.store(InitResult, std::memory_order_release);
    });

    return InitResult;
}

WIREHAIR_EXPORT WirehairResult wirehair_wire_profile_init(
    uint64_t profileId,
    WirehairWireProfile* profileOut)
{
    if (!profileOut) {
        return Wirehair_InvalidInput;
    }

    profileOut->struct_bytes = 0;
    profileOut->profile_version = 0;
    profileOut->profile_id = 0;

    wirehair::WireProfile ignored;
    if (!ResolveWireProfileId(profileId, ignored)) {
        return Wirehair_InvalidInput;
    }
    if (!NamedLegacyWireProfilesAvailable()) {
        return Wirehair_UnsupportedPlatform;
    }

    profileOut->struct_bytes = sizeof(WirehairWireProfile);
    profileOut->profile_version = WIREHAIR_WIRE_PROFILE_VERSION;
    profileOut->profile_id = profileId;
    return Wirehair_Success;
}

WIREHAIR_EXPORT WirehairCodec wirehair_encoder_create(
    WirehairCodec reuseOpt, ///< [Optional] Pointer to prior codec object
    const void*    message, ///< Pointer to message
    uint64_t  messageBytes, ///< Bytes in the message
    uint32_t    blockBytes  ///< Bytes in an output block
)
{
    WirehairCodec codec = nullptr;
    (void)CreateEncoder(
        reuseOpt, message, messageBytes, blockBytes, false,
        wirehair::WireProfile::LegacyCurrent, &codec);
    return codec;
}

WIREHAIR_EXPORT WirehairResult wirehair_encoder_create_ex(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    WirehairCodec* codecOut)
{
    return CreateEncoder(
        reuseOpt, message, messageBytes, blockBytes, false,
        wirehair::WireProfile::LegacyCurrent, codecOut);
}

WIREHAIR_EXPORT WirehairCodec wirehair_encoder_create_owned(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes)
{
    WirehairCodec codec = nullptr;
    (void)CreateEncoder(
        reuseOpt, message, messageBytes, blockBytes, true,
        wirehair::WireProfile::LegacyCurrent, &codec);
    return codec;
}

WIREHAIR_EXPORT WirehairResult wirehair_encoder_create_owned_ex(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    WirehairCodec* codecOut)
{
    return CreateEncoder(
        reuseOpt, message, messageBytes, blockBytes, true,
        wirehair::WireProfile::LegacyCurrent, codecOut);
}

WIREHAIR_EXPORT WirehairResult wirehair_encoder_create_profile_ex(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    const WirehairWireProfile* profile,
    uint32_t flags,
    WirehairCodec* codecOut)
{
    if (!codecOut) {
        return Wirehair_InvalidInput;
    }
    wirehair::WireProfile resolved;
    if ((flags & ~WIREHAIR_ENCODER_OWN_INPUT) != 0 ||
        !ResolveWireProfileDescriptor(profile, resolved))
    {
        *codecOut = nullptr;
        delete reinterpret_cast<wirehair::Codec*>(reuseOpt);
        return Wirehair_InvalidInput;
    }
    if (!NamedLegacyWireProfilesAvailable()) {
        *codecOut = nullptr;
        delete reinterpret_cast<wirehair::Codec*>(reuseOpt);
        return Wirehair_UnsupportedPlatform;
    }
    return CreateEncoder(
        reuseOpt, message, messageBytes, blockBytes,
        (flags & WIREHAIR_ENCODER_OWN_INPUT) != 0,
        resolved, codecOut);
}

WIREHAIR_EXPORT WirehairResult wirehair_encode(
    WirehairCodec    codec, ///< Pointer to codec from wirehair_encoder_init()
    unsigned       blockId, ///< Identifier of block to generate
    void*     blockDataOut, ///< Pointer to output block data
    uint32_t      outBytes, ///< Bytes in the output buffer
    uint32_t* dataBytesOut  ///< Number of bytes written <= blockBytes
)
{
    if (dataBytesOut) {
        *dataBytesOut = 0;
    }
    if (!codec || !blockDataOut || !dataBytesOut) {
        return Wirehair_InvalidInput;
    }

    wirehair::Codec* session = reinterpret_cast<wirehair::Codec*>(codec);
    if (!session->CanEncode()) {
        return Wirehair_InvalidInput;
    }

    const uint32_t writtenBytes = session->Encode(blockId, blockDataOut, outBytes);
    *dataBytesOut = writtenBytes;

    if (writtenBytes <= 0) {
        return Wirehair_InvalidInput;
    }

    return Wirehair_Success;
}

WIREHAIR_EXPORT WirehairCodec wirehair_decoder_create(
    WirehairCodec reuseOpt, ///< Codec object to reuse
    uint64_t  messageBytes, ///< Bytes in the message to decode
    uint32_t    blockBytes  ///< Bytes in each encoded block
)
{
    WirehairCodec codec = nullptr;
    (void)CreateDecoder(
        reuseOpt, messageBytes, blockBytes,
        wirehair::WireProfile::LegacyCurrent, &codec);
    return codec;
}

WIREHAIR_EXPORT WirehairResult wirehair_decoder_create_ex(
    WirehairCodec reuseOpt,
    uint64_t messageBytes,
    uint32_t blockBytes,
    WirehairCodec* codecOut)
{
    return CreateDecoder(
        reuseOpt, messageBytes, blockBytes,
        wirehair::WireProfile::LegacyCurrent, codecOut);
}

WIREHAIR_EXPORT WirehairResult wirehair_decoder_create_profile_ex(
    WirehairCodec reuseOpt,
    uint64_t messageBytes,
    uint32_t blockBytes,
    const WirehairWireProfile* profile,
    WirehairCodec* codecOut)
{
    if (!codecOut) {
        return Wirehair_InvalidInput;
    }
    wirehair::WireProfile resolved;
    if (!ResolveWireProfileDescriptor(profile, resolved))
    {
        *codecOut = nullptr;
        delete reinterpret_cast<wirehair::Codec*>(reuseOpt);
        return Wirehair_InvalidInput;
    }
    if (!NamedLegacyWireProfilesAvailable()) {
        *codecOut = nullptr;
        delete reinterpret_cast<wirehair::Codec*>(reuseOpt);
        return Wirehair_UnsupportedPlatform;
    }
    return CreateDecoder(
        reuseOpt, messageBytes, blockBytes, resolved, codecOut);
}

WIREHAIR_EXPORT WirehairResult wirehair_decode(
    WirehairCodec   codec, ///< Codec object
    unsigned      blockId, ///< ID number of received block
    const void* blockData, ///< Pointer to block data
    uint32_t    dataBytes  ///< Number of bytes in the data block
)
{
    // If input is invalid:
    if (!codec || !blockData || dataBytes < 1) {
        return Wirehair_InvalidInput;
    }

    wirehair::Codec* decoder = reinterpret_cast<wirehair::Codec*>(codec);
    if (!decoder->CanDecode()) {
        return Wirehair_InvalidInput;
    }

    return decoder->DecodeFeed(blockId, blockData, dataBytes);
}

WIREHAIR_EXPORT WirehairResult wirehair_recover(
    WirehairCodec    codec, ///< Codec object
    void*       messageOut, ///< Buffer where reconstructed message will be written
    uint64_t  messageBytes  ///< Bytes in the message
)
{
    // If input is invalid:
    if (!codec || !messageOut) {
        return Wirehair_InvalidInput;
    }

    wirehair::Codec* decoder = reinterpret_cast<wirehair::Codec*>(codec);

    return decoder->ReconstructOutput(messageOut, messageBytes);
}

WIREHAIR_EXPORT WirehairResult wirehair_recover_block(
    WirehairCodec codec, ///< Codec object
    unsigned    blockId, ///< ID of the block to reconstruct between 0..N-1
    void*  blockDataOut, ///< Pointer to block data
    uint32_t*  bytesOut  ///< Set to the number of data bytes in the block
)
{
    return wirehair_recover_block_ex(
        codec, blockId, blockDataOut, UINT32_MAX, bytesOut);
}

WIREHAIR_EXPORT WirehairResult wirehair_recover_block_ex(
    WirehairCodec codec,
    unsigned blockId,
    void* blockDataOut,
    uint32_t outputCapacity,
    uint32_t* bytesOut)
{
    if (bytesOut) {
        *bytesOut = 0;
    }
    if (!codec || !blockDataOut || !bytesOut || blockId > 0xffffu) {
        return Wirehair_InvalidInput;
    }

    wirehair::Codec* decoder = reinterpret_cast<wirehair::Codec*>(codec);
    return decoder->ReconstructBlock(
        static_cast<uint16_t>(blockId),
        blockDataOut,
        outputCapacity,
        bytesOut);
}

WIREHAIR_EXPORT WirehairResult wirehair_decoder_becomes_encoder(
    WirehairCodec codec ///< Codec to change
)
{
    // If input is invalid:
    if (!codec) {
        return Wirehair_InvalidInput;
    }

    wirehair::Codec* encoder = reinterpret_cast<wirehair::Codec*>(codec);

    return encoder->InitializeEncoderFromDecoder();
}

WIREHAIR_EXPORT void wirehair_free(
    WirehairCodec codec ///< Codec object to free
)
{
    wirehair::Codec* object = reinterpret_cast<wirehair::Codec*>(codec);

    delete object;
}


} // extern "C"
