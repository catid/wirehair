#include "V2FuzzDriver.h"

#include <wirehair/wirehair.h>
#include <wirehair/wirehair.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#ifndef WIREHAIR_FUZZ_TARGET_NAME
#define WIREHAIR_FUZZ_TARGET_NAME "borrowed-facade"
#endif

#ifndef WIREHAIR_FUZZ_CORPUS_MANIFEST
#define WIREHAIR_FUZZ_CORPUS_MANIFEST \
    "codec/fuzz/corpus/stateful/manifest.txt"
#endif

namespace {

bool Fail(std::string& failure, const char* message)
{
    failure = message;
    return false;
}

struct CodecOwner
{
    WirehairV2Codec Handle = nullptr;

    ~CodecOwner()
    {
        wirehair_v2_free(Handle);
    }

    CodecOwner() = default;
    CodecOwner(const CodecOwner&) = delete;
    CodecOwner& operator=(const CodecOwner&) = delete;
};

struct Fixture
{
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    uint64_t MessageBytes = 0u;
    std::vector<uint32_t> Storage;

    uint8_t* Message()
    {
        return reinterpret_cast<uint8_t*>(Storage.data());
    }

    const uint8_t* Message() const
    {
        return reinterpret_cast<const uint8_t*>(Storage.data());
    }
};

Fixture MakeFixture(wirehair_v2::fuzz::Input& input)
{
    Fixture fixture;
    fixture.K = 2u + input.U8() % 15u;
    fixture.BlockBytes = 4u + input.U8() % 61u;
    const uint32_t final_bytes = 1u + input.U8() % fixture.BlockBytes;
    fixture.MessageBytes =
        (uint64_t)(fixture.K - 1u) * fixture.BlockBytes + final_bytes;
    fixture.Storage.assign(
        ((size_t)fixture.MessageBytes + sizeof(uint32_t) - 1u) /
            sizeof(uint32_t),
        0u);
    uint64_t state = input.U64() | 1u;
    for (size_t i = 0; i < (size_t)fixture.MessageBytes; ++i)
    {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        fixture.Message()[i] = (uint8_t)(
            state * UINT64_C(2685821657736338717));
    }
    return fixture;
}

uint32_t RequiredBytes(const Fixture& fixture, uint32_t id)
{
    if (id == fixture.K - 1u) {
        return (uint32_t)(fixture.MessageBytes -
            (uint64_t)(fixture.K - 1u) * fixture.BlockBytes);
    }
    return fixture.BlockBytes;
}

bool EncodeExact(
    WirehairV2Codec codec,
    uint32_t id,
    uint32_t block_bytes,
    std::vector<uint8_t>& output,
    uint32_t& written)
{
    output.assign(block_bytes + 2u, 0xa5u);
    written = UINT32_MAX;
    const WirehairV2Result result = wirehair_v2_encode(
        codec, id, output.data() + 1u, block_bytes, &written);
    if (result != WirehairV2_Success || written == 0u ||
        written > block_bytes || output.front() != 0xa5u ||
        output.back() != 0xa5u)
    {
        return false;
    }
    output.erase(output.begin());
    output.resize(written);
    return true;
}

bool ComparePacket(
    WirehairV2Codec first,
    WirehairV2Codec second,
    uint32_t id,
    uint32_t block_bytes)
{
    std::vector<uint8_t> a;
    std::vector<uint8_t> b;
    uint32_t a_bytes = 0u;
    uint32_t b_bytes = 0u;
    return EncodeExact(first, id, block_bytes, a, a_bytes) &&
        EncodeExact(second, id, block_bytes, b, b_bytes) &&
        a_bytes == b_bytes && a == b;
}

WirehairV2Result CreateBorrowed(
    unsigned route,
    const Fixture& fixture,
    const std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES>& profile,
    const WirehairV2EncoderOptions* options,
    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES>& output_profile,
    uint32_t& output_bytes,
    WirehairV2Codec* codec_out)
{
    if (route == 0u) {
        return wirehair_v2_encoder_create_with_options(
            fixture.Message(), fixture.MessageBytes, fixture.BlockBytes,
            options, output_profile.data(), (uint32_t)output_profile.size(),
            &output_bytes, codec_out);
    }
    if (route == 1u) {
        return wirehair_v2_encoder_create_profile_id_with_options(
            WIREHAIR_V2_PROFILE_CURRENT,
            fixture.Message(), fixture.MessageBytes, fixture.BlockBytes,
            options, output_profile.data(), (uint32_t)output_profile.size(),
            &output_bytes, codec_out);
    }
    return wirehair_v2_encoder_create_profile_with_options(
        fixture.Message(), profile.data(), (uint32_t)profile.size(),
        options, codec_out);
}

WirehairV2Result ExpectedOptionsResult(
    const WirehairV2EncoderOptions* options)
{
    if (!options) return WirehairV2_InvalidInput;
    if (options->struct_bytes != sizeof(WirehairV2EncoderOptions)) {
        return WirehairV2_InvalidSize;
    }
    if (options->options_version !=
        WIREHAIR_V2_ENCODER_OPTIONS_VERSION)
    {
        return WirehairV2_UnsupportedVersion;
    }
    if (options->reserved != 0u) {
        return WirehairV2_ReservedNonzero;
    }
    if (options->source_policy !=
            (uint32_t)WirehairV2EncoderSource_Independent &&
        options->source_policy !=
            (uint32_t)WirehairV2EncoderSource_BorrowedImmutable)
    {
        return WirehairV2_InvalidInput;
    }
    return WirehairV2_Success;
}

bool MalformedOptions(
    wirehair_v2::fuzz::Input& input,
    const Fixture& fixture,
    const std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES>& profile,
    std::string& failure)
{
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    const WirehairV2EncoderOptions* options_pointer = &options;
    switch (input.U8() % 8u)
    {
    case 0: options_pointer = nullptr; break;
    case 1: options.struct_bytes = 0u; break;
    case 2: ++options.options_version; break;
    case 3: options.reserved = input.U32() | 1u; break;
    case 4: options.source_policy = WirehairV2EncoderSource_Invalid; break;
    case 5: options.source_policy = WirehairV2EncoderSource_Count; break;
    case 6: options.source_policy = WirehairV2EncoderSource_Padding; break;
    default: options.source_policy = input.U32() | UINT32_C(0x80000000); break;
    }

    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> output;
    output.fill(0x5au);
    const std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> before =
        output;
    uint32_t output_bytes = UINT32_MAX;
    CodecOwner codec;
    codec.Handle = reinterpret_cast<WirehairV2Codec>(
        static_cast<uintptr_t>(1u));
    const unsigned route = input.U8() % 3u;
    const WirehairV2Result result = CreateBorrowed(
        route,
        fixture,
        profile,
        options_pointer,
        output,
        output_bytes,
        &codec.Handle);
    const WirehairV2Result expected = ExpectedOptionsResult(options_pointer);
    if (result != expected || codec.Handle != nullptr) {
        codec.Handle = nullptr;
        return Fail(failure, "malformed public options result/state mismatch");
    }
    if (route != 2u &&
        (output != before ||
         output_bytes != WIREHAIR_V2_PROFILE_SERIALIZED_BYTES))
    {
        return Fail(failure, "malformed selector options modified output");
    }
    return true;
}

bool BorrowedSequence(
    wirehair_v2::fuzz::Input& input,
    Fixture& fixture,
    WirehairV2Codec baseline,
    const std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES>& profile,
    std::string& failure)
{
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> output_profile{};
    uint32_t output_bytes = 0u;
    CodecOwner borrowed;
    const unsigned route = input.U8() % 3u;
    if (CreateBorrowed(
            route,
            fixture,
            profile,
            &options,
            output_profile,
            output_bytes,
            &borrowed.Handle) != WirehairV2_Success ||
        !borrowed.Handle ||
        (route != 2u &&
            (output_bytes != output_profile.size() ||
             output_profile != profile)))
    {
        return Fail(failure, "valid public borrowed construction failed");
    }

    const uint32_t ids[] = {
        0u,
        fixture.K - 1u,
        fixture.K,
        fixture.K + 7u,
        UINT32_MAX
    };
    for (uint32_t id : ids) {
        if (!ComparePacket(
                baseline, borrowed.Handle, id, fixture.BlockBytes))
        {
            return Fail(failure, "borrowed/default packet identity mismatch");
        }
    }

    const std::vector<uint32_t> source_before = fixture.Storage;
    uint32_t counter = UINT32_C(0xa5a5a5a5);
    if (wirehair_v2_encode(
            borrowed.Handle,
            fixture.K,
            fixture.Message(),
            0u,
            &counter) != WirehairV2_InvalidInput ||
        counter != UINT32_C(0xa5a5a5a5) ||
        fixture.Storage != source_before)
    {
        return Fail(failure, "borrowed packet/source overlap was not no-write");
    }
    uint32_t* const source_counter = fixture.Storage.data();
    const uint32_t source_counter_before = *source_counter;
    if (wirehair_v2_encode(
            borrowed.Handle,
            fixture.K,
            nullptr,
            fixture.BlockBytes,
            source_counter) != WirehairV2_InvalidInput ||
        *source_counter != source_counter_before ||
        fixture.Storage != source_before)
    {
        return Fail(failure, "borrowed counter/source overlap was not no-write");
    }

    std::vector<uint8_t> short_output(fixture.BlockBytes, 0x3cu);
    const std::vector<uint8_t> short_before = short_output;
    uint32_t required = 0u;
    if (wirehair_v2_encode(
            borrowed.Handle,
            fixture.K,
            short_output.data(),
            fixture.BlockBytes - 1u,
            &required) != WirehairV2_BufferTooSmall ||
        required != RequiredBytes(fixture, fixture.K) ||
        short_output != short_before)
    {
        return Fail(failure, "borrowed short output contract mismatch");
    }

    if (wirehair_v2_encoder_detach_input(borrowed.Handle) !=
            WirehairV2_Success ||
        wirehair_v2_encoder_detach_input(borrowed.Handle) !=
            WirehairV2_Success)
    {
        return Fail(failure, "public borrowed detach was not idempotent");
    }
    std::fill(fixture.Storage.begin(), fixture.Storage.end(), UINT32_MAX);
    std::vector<uint32_t>().swap(fixture.Storage);
    for (uint32_t id : ids) {
        if (!ComparePacket(
                baseline, borrowed.Handle, id, fixture.BlockBytes))
        {
            return Fail(failure, "detached public packet identity mismatch");
        }
    }
    return true;
}

bool WrapperTransaction(
    const Fixture& fixture,
    WirehairV2Codec baseline,
    std::string& failure)
{
    std::vector<uint8_t> source(
        fixture.Message(), fixture.Message() + fixture.MessageBytes);
    wirehair::v2::SerializedProfile profile;
    wirehair::v2::Encoder encoder;
    if (encoder.CreateBorrowed(
            source.data(), source.size(), fixture.BlockBytes, profile) !=
            WirehairV2_Success)
    {
        return Fail(failure, "C++ borrowed wrapper construction failed");
    }
    wirehair::v2::SerializedProfile replacement_profile;
    if (encoder.CreateBorrowed(
            UINT64_C(0x0123456789abcdef),
            source.data(), source.size(), fixture.BlockBytes,
            replacement_profile) != WirehairV2_UnsupportedProfile ||
        !encoder)
    {
        return Fail(failure, "C++ failed replacement was not transactional");
    }
    wirehair::v2::Encoder moved(std::move(encoder));
    if (encoder || !moved || moved.DetachInput() != WirehairV2_Success ||
        moved.DetachInput() != WirehairV2_Success)
    {
        return Fail(failure, "C++ move/detach contract mismatch");
    }
    std::fill(source.begin(), source.end(), 0xccu);
    std::vector<uint8_t>().swap(source);
    std::vector<uint8_t> expected;
    uint32_t expected_bytes = 0u;
    if (!EncodeExact(
            baseline,
            fixture.K + 7u,
            fixture.BlockBytes,
            expected,
            expected_bytes))
    {
        return Fail(failure, "C++ baseline packet generation failed");
    }
    std::vector<uint8_t> actual(fixture.BlockBytes, 0u);
    uint32_t actual_bytes = 0u;
    if (moved.Encode(
            fixture.K + 7u,
            actual.data(),
            fixture.BlockBytes,
            actual_bytes) != WirehairV2_Success ||
        actual_bytes != expected_bytes ||
        !std::equal(
            actual.begin(), actual.begin() + actual_bytes,
            expected.begin()))
    {
        return Fail(failure, "C++ detached packet identity mismatch");
    }
    return true;
}

bool FuzzBorrowedFacadeCase(
    const uint8_t* data,
    size_t size,
    std::string& failure)
{
    if (size > wirehair_v2::fuzz::kMaxFuzzInputBytes) return true;
    wirehair_v2::fuzz::Input input(data, size);
    Fixture fixture = MakeFixture(input);

    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> profile{};
    uint32_t profile_bytes = 0u;
    CodecOwner baseline;
    if (wirehair_v2_encoder_create(
            fixture.Message(),
            fixture.MessageBytes,
            fixture.BlockBytes,
            profile.data(),
            (uint32_t)profile.size(),
            &profile_bytes,
            &baseline.Handle) != WirehairV2_Success ||
        profile_bytes != profile.size() || !baseline.Handle)
    {
        return Fail(failure, "public baseline construction failed");
    }

    switch (input.U8() % 3u)
    {
    case 0:
        return MalformedOptions(input, fixture, profile, failure);
    case 1:
        return WrapperTransaction(fixture, baseline.Handle, failure);
    default:
        return BorrowedSequence(
            input, fixture, baseline.Handle, profile, failure);
    }
}

} // namespace

#if defined(WIREHAIR_ENABLE_LIBFUZZER)
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    wirehair_v2::fuzz::RunCoverageGuidedCaseOrAbort(
        WIREHAIR_FUZZ_TARGET_NAME, FuzzBorrowedFacadeCase, data, size);
    return 0;
}
#else
int main(int argc, char** argv)
{
    return wirehair_v2::fuzz::RunDeterministicFuzzer(
        argc,
        argv,
        WIREHAIR_FUZZ_TARGET_NAME,
        WIREHAIR_FUZZ_CORPUS_MANIFEST,
        FuzzBorrowedFacadeCase);
}
#endif
