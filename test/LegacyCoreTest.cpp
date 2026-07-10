#include <wirehair/wirehair.h>

#include "../WirehairCodec.h"
#include "../WirehairEnvironment.h"
#include "../gf256.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <thread>
#include <vector>

namespace {

int Failures = 0;

#define CHECK(expression) do { \
    if (!(expression)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", \
            __FILE__, __LINE__, #expression); \
        ++Failures; \
    } \
} while (0)

bool SetEnv(const char* name, const char* value)
{
#if defined(_WIN32)
    return _putenv_s(name, value ? value : "") == 0;
#else
    if (value) {
        return setenv(name, value, 1) == 0;
    }
    return unsetenv(name) == 0;
#endif
}

void TestEnvironmentValue()
{
    const char* name = "WIREHAIR_ENVIRONMENT_VALUE_TEST";
    const char* other_name = "WIREHAIR_ENVIRONMENT_VALUE_OTHER_TEST";
    CHECK(SetEnv(name, nullptr));
    {
        const wirehair::EnvironmentValue environment(name);
        CHECK(!environment.IsSet());
        CHECK(environment.Get() == nullptr);
    }

#if !defined(_WIN32)
    // The Windows CRT's _putenv_s removes a variable for an empty value.
    CHECK(SetEnv(name, ""));
    {
        const wirehair::EnvironmentValue environment(name);
        CHECK(environment.IsSet());
        CHECK(environment.Get() != nullptr);
        if (environment.Get()) {
            CHECK(environment.Get()[0] == '\0');
        }
    }
#endif

    CHECK(SetEnv(name, "12345"));
    {
        const wirehair::EnvironmentValue environment(name);
        CHECK(environment.IsSet());
        CHECK(environment.Get() != nullptr);
        if (environment.Get()) {
            CHECK(std::strcmp(environment.Get(), "12345") == 0);
        }
    }

    CHECK(SetEnv(name, "first-value"));
    CHECK(SetEnv(other_name, "second-value"));
    {
        const wirehair::EnvironmentValue first(name);
        const wirehair::EnvironmentValue second(other_name);
        CHECK(SetEnv(name, "mutated-value"));
        CHECK(first.IsSet());
        CHECK(second.IsSet());
        if (first.Get()) {
            CHECK(std::strcmp(first.Get(), "first-value") == 0);
        }
        if (second.Get()) {
            CHECK(std::strcmp(second.Get(), "second-value") == 0);
        }
    }
    CHECK(SetEnv(name, nullptr));
    CHECK(SetEnv(other_name, nullptr));
}

std::vector<uint8_t> MakeMessage(size_t bytes, uint32_t salt = 0)
{
    std::vector<uint8_t> message(bytes);
    for (size_t i = 0; i < bytes; ++i) {
        message[i] = static_cast<uint8_t>((i * 73u + salt * 29u + 11u) & 0xffu);
    }
    return message;
}

bool BasicRoundTrip(uint32_t salt)
{
    const uint32_t block_bytes = 37;
    const uint64_t message_bytes = block_bytes * 8u - 9u;
    const std::vector<uint8_t> message = MakeMessage(
        static_cast<size_t>(message_bytes), salt);

    WirehairCodec encoder = wirehair_encoder_create_owned(
        nullptr, message.data(), message_bytes, block_bytes);
    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message_bytes, block_bytes);
    if (!encoder || !decoder) {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return false;
    }

    uint8_t block[37];
    WirehairResult result = Wirehair_NeedMore;
    for (unsigned id = 0; id < 8; ++id)
    {
        uint32_t written = 0;
        if (wirehair_encode(encoder, id, block, sizeof(block), &written) !=
            Wirehair_Success)
        {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return false;
        }
        result = wirehair_decode(decoder, id, block, written);
    }

    std::vector<uint8_t> recovered(static_cast<size_t>(message_bytes));
    const bool ok = result == Wirehair_Success &&
        wirehair_recover(decoder, recovered.data(), message_bytes) ==
            Wirehair_Success &&
        recovered == message;
    wirehair_free(encoder);
    wirehair_free(decoder);
    return ok;
}

int RunPreinit()
{
    uint8_t message[64] = {};
    WirehairCodec codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_encoder_create_ex(
        nullptr, message, sizeof(message), 32, &codec) == Wirehair_Error);
    CHECK(codec == nullptr);
    CHECK(wirehair_encoder_create(nullptr, message, sizeof(message), 32) == nullptr);

    codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_decoder_create_ex(
        nullptr, sizeof(message), 32, &codec) == Wirehair_Error);
    CHECK(codec == nullptr);
    CHECK(wirehair_init_(WIREHAIR_VERSION + 1) == Wirehair_InvalidInput);
    return Failures == 0 ? 0 : 1;
}

int RunInitFailure()
{
    SetEnv("WIREHAIR_GF256_TEST_INIT_RESULT", "-3");
    CHECK(wirehair_init_(WIREHAIR_VERSION + 1) == Wirehair_InvalidInput);
    CHECK(wirehair_init() == Wirehair_UnsupportedPlatform);
    SetEnv("WIREHAIR_GF256_TEST_INIT_RESULT", nullptr);
    CHECK(wirehair_init() == Wirehair_UnsupportedPlatform);
    CHECK(gf256_init_(GF256_VERSION + 1) == -1);
    CHECK(gf256_init() == -3);

    uint8_t message[64] = {};
    WirehairCodec codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_encoder_create_ex(
        nullptr, message, sizeof(message), 32, &codec) ==
        Wirehair_UnsupportedPlatform);
    CHECK(codec == nullptr);
    return Failures == 0 ? 0 : 1;
}

int RunColdStart()
{
    const unsigned thread_count = 24;
    const unsigned iterations = 6;
    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    for (unsigned thread_i = 0; thread_i < thread_count; ++thread_i)
    {
        threads.emplace_back([&, thread_i]() {
            ready.fetch_add(1, std::memory_order_release);
            while (!go.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }

            if ((thread_i & 1u) != 0 &&
                wirehair_init_(WIREHAIR_VERSION + 1) != Wirehair_InvalidInput)
            {
                failures.fetch_add(1, std::memory_order_relaxed);
            }
            if ((thread_i % 3u) == 0 &&
                gf256_init_(GF256_VERSION + 1) != -1)
            {
                failures.fetch_add(1, std::memory_order_relaxed);
            }
            if (wirehair_init() != Wirehair_Success) {
                failures.fetch_add(1, std::memory_order_relaxed);
                return;
            }
            for (unsigned i = 0; i < iterations; ++i) {
                if (!BasicRoundTrip(thread_i * iterations + i)) {
                    failures.fetch_add(1, std::memory_order_relaxed);
                }
            }
        });
    }

    while (ready.load(std::memory_order_acquire) != thread_count) {
        std::this_thread::yield();
    }
    go.store(true, std::memory_order_release);
    for (std::thread& thread : threads) {
        thread.join();
    }

    CHECK(failures.load(std::memory_order_relaxed) == 0);
    return Failures == 0 ? 0 : 1;
}

void TestCpuFeatureSelection()
{
    gf256_x86_cpu_snapshot snapshot = {};
    gf256_x86_cpu_features features;

    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(!features.SSSE3 && !features.AVX2 && !features.GFNI &&
        !features.AVX512);

    snapshot.MaxBasicLeaf = 7;
    snapshot.Leaf1ECX = (1u << 9) | (1u << 28);
    snapshot.Leaf7EBX = (1u << 5) | (1u << 16) | (1u << 30);
    snapshot.Leaf7ECX = (1u << 8);
    snapshot.XCR0 = UINT64_C(0xe6);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.SSSE3 && !features.AVX2 && !features.GFNI &&
        !features.AVX512);

    snapshot.Leaf1ECX |= (1u << 27);
    snapshot.XCR0 = UINT64_C(0x2);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(!features.AVX2 && !features.GFNI && !features.AVX512);

    snapshot.XCR0 = UINT64_C(0x6);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.AVX2 && !features.GFNI && !features.AVX512);

    snapshot.XCR0 = UINT64_C(0xe6);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.AVX2 && features.GFNI && features.AVX512);

    snapshot.Leaf7EBX &= ~(1u << 30);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.AVX2 && !features.GFNI && features.AVX512);

    snapshot.MaxBasicLeaf = 1;
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.SSSE3 && !features.AVX2 && !features.GFNI &&
        !features.AVX512);

    const gf256_x86_cpu_snapshot full = {
        7,
        (1u << 9) | (1u << 27) | (1u << 28),
        (1u << 5) | (1u << 16) | (1u << 30),
        (1u << 8),
        UINT64_C(0xe6)
    };
    snapshot = full;
    snapshot.Leaf1ECX &= ~(1u << 9);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(!features.SSSE3 && features.AVX2 && features.GFNI &&
        features.AVX512);

    snapshot = full;
    snapshot.Leaf1ECX &= ~(1u << 28);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.SSSE3 && !features.AVX2 && !features.GFNI &&
        !features.AVX512);

    for (uint64_t missing : { UINT64_C(0x2), UINT64_C(0x4) })
    {
        snapshot = full;
        snapshot.XCR0 &= ~missing;
        gf256_select_x86_cpu_features(&snapshot, &features);
        CHECK(features.SSSE3 && !features.AVX2 && !features.GFNI &&
            !features.AVX512);
    }
    for (uint64_t missing : {
            UINT64_C(0x20), UINT64_C(0x40), UINT64_C(0x80) })
    {
        snapshot = full;
        snapshot.XCR0 &= ~missing;
        gf256_select_x86_cpu_features(&snapshot, &features);
        CHECK(features.SSSE3 && features.AVX2 && !features.GFNI &&
            !features.AVX512);
    }

    snapshot = full;
    snapshot.Leaf7EBX &= ~(1u << 5);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.SSSE3 && !features.AVX2 && features.GFNI &&
        features.AVX512);

    snapshot = full;
    snapshot.Leaf7EBX &= ~(1u << 16);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.SSSE3 && features.AVX2 && !features.GFNI &&
        !features.AVX512);

    snapshot = full;
    snapshot.Leaf7ECX &= ~(1u << 8);
    gf256_select_x86_cpu_features(&snapshot, &features);
    CHECK(features.SSSE3 && features.AVX2 && !features.GFNI &&
        features.AVX512);

    gf256_select_x86_cpu_features(nullptr, &features);
    CHECK(!features.SSSE3 && !features.AVX2 && !features.GFNI &&
        !features.AVX512);
    gf256_select_x86_cpu_features(&snapshot, nullptr);
}

void TestCreationResultsAndBoundaries()
{
    uint8_t message[128] = {};
    WirehairCodec codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));

    CHECK(wirehair_encoder_create_ex(
        nullptr, nullptr, 64, 32, &codec) == Wirehair_InvalidInput);
    CHECK(codec == nullptr);
    codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_encoder_create_ex(
        nullptr, message, 0, 32, &codec) == Wirehair_InvalidInput);
    CHECK(codec == nullptr);
    CHECK(wirehair_encoder_create_ex(
        nullptr, message, 1, 1, &codec) == Wirehair_BadInput_SmallN);
    CHECK(codec == nullptr);
    CHECK(wirehair_decoder_create_ex(
        nullptr, std::numeric_limits<uint64_t>::max(), 2, &codec) ==
        Wirehair_BadInput_LargeN);
    CHECK(codec == nullptr);
    CHECK(wirehair_decoder_create_ex(
        nullptr, 64, 0x80000000u, &codec) == Wirehair_InvalidInput);
    CHECK(codec == nullptr);
    CHECK(wirehair_decoder_create_ex(
        nullptr, 63, 32, &codec) == Wirehair_Success);
    CHECK(reinterpret_cast<wirehair::Codec*>(codec)->BlockCount() == 2);
    wirehair_free(codec);
    codec = nullptr;
    CHECK(wirehair_decoder_create_ex(
        nullptr, 65, 32, &codec) == Wirehair_Success);
    CHECK(reinterpret_cast<wirehair::Codec*>(codec)->BlockCount() == 3);
    wirehair_free(codec);

    struct BoundaryCase {
        uint64_t message_bytes;
        uint32_t block_bytes;
        WirehairResult expected;
        uint32_t expected_blocks;
    };
    const BoundaryCase boundaries[] = {
        { 0, 0, Wirehair_InvalidInput, 0 },
        { 0, 1, Wirehair_InvalidInput, 0 },
        { 1, 0, Wirehair_InvalidInput, 0 },
        { 1, 1, Wirehair_BadInput_SmallN, 0 },
        { 2, 1, Wirehair_Success, 2 },
        { 63, 32, Wirehair_Success, 2 },
        { 64, 32, Wirehair_Success, 2 },
        { 65, 32, Wirehair_Success, 3 },
        { UINT64_MAX - 3, 1, Wirehair_BadInput_LargeN, 0 },
        { UINT64_MAX - 2, 2, Wirehair_BadInput_LargeN, 0 },
        { UINT64_MAX - 1, 1, Wirehair_BadInput_LargeN, 0 },
        { UINT64_MAX, 2, Wirehair_BadInput_LargeN, 0 },
        { UINT64_MAX, UINT32_MAX, Wirehair_InvalidInput, 0 }
    };
    for (const BoundaryCase& boundary : boundaries)
    {
        wirehair::Codec candidate;
        const WirehairResult result = candidate.InitializeDecoder(
            boundary.message_bytes, boundary.block_bytes);
        CHECK(result == boundary.expected);
        if (result == Wirehair_Success) {
            CHECK(candidate.BlockCount() == boundary.expected_blocks);
        }
    }

    wirehair::Codec internal;
    uint64_t state = UINT64_C(0x9e3779b97f4a7c15);
    for (unsigned i = 0; i < 256; ++i)
    {
        state = state * UINT64_C(6364136223846793005) + 1;
        const uint32_t block_bytes = 100u + static_cast<uint32_t>(state % 901u);
        state = state * UINT64_C(6364136223846793005) + 1;
        const uint64_t message_bytes = 1u + state % 100000u;
        const uint64_t expected = message_bytes / block_bytes +
            (message_bytes % block_bytes != 0 ? 1 : 0);
        const WirehairResult result = internal.InitializeDecoder(
            message_bytes, block_bytes);
        if (expected < 2) {
            CHECK(result == Wirehair_BadInput_SmallN);
        }
        else {
            CHECK(result == Wirehair_Success);
            CHECK(internal.BlockCount() == expected);
        }
    }

    SetEnv("WIREHAIR_TEST_FORCE_OOM", "1");
    codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_encoder_create_owned_ex(
        nullptr, message, 64, 32, &codec) == Wirehair_OOM);
    CHECK(codec == nullptr);
    SetEnv("WIREHAIR_TEST_FORCE_OOM", nullptr);

    const WirehairResult injected_results[] = {
        Wirehair_BadDenseSeed,
        Wirehair_BadPeelSeed
    };
    for (WirehairResult injected : injected_results)
    {
        char value[16];
        std::snprintf(value, sizeof(value), "%d", static_cast<int>(injected));
        SetEnv("WIREHAIR_TEST_FORCE_CREATE_RESULT", value);
        codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
        CHECK(wirehair_encoder_create_ex(
            nullptr, message, 64, 32, &codec) == injected);
        CHECK(codec == nullptr);
    }
    SetEnv("WIREHAIR_TEST_FORCE_CREATE_RESULT", nullptr);

    codec = wirehair_encoder_create(nullptr, message, 64, 32);
    CHECK(codec != nullptr);
    CHECK(wirehair_decoder_create_ex(codec, 64, 32, nullptr) ==
        Wirehair_InvalidInput);
    uint8_t block[32] = {};
    uint32_t written = 0;
    CHECK(wirehair_encode(codec, 0, block, sizeof(block), &written) ==
        Wirehair_Success);
    wirehair_free(codec);
}

WirehairCodec DecodeLossy(
    WirehairCodec encoder,
    const std::vector<uint8_t>& expected,
    uint32_t block_bytes)
{
    const uint64_t message_bytes = expected.size();
    const unsigned block_count = static_cast<unsigned>(
        message_bytes / block_bytes + (message_bytes % block_bytes != 0));
    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message_bytes, block_bytes);
    CHECK(decoder != nullptr);
    if (!decoder) {
        return nullptr;
    }

    std::vector<unsigned> ids;
    for (unsigned id = 0; id < block_count; ++id) {
        if (id % 5u != 1u) {
            ids.push_back(id);
        }
    }
    for (unsigned id = block_count; id < block_count + 96u; ++id) {
        ids.push_back(id);
    }
    for (size_t i = 0; i < ids.size(); ++i) {
        const size_t j = (i * 47u + 13u) % ids.size();
        std::swap(ids[i], ids[j]);
    }

    std::vector<uint8_t> block(block_bytes);
    WirehairResult decode_result = Wirehair_NeedMore;
    for (unsigned id : ids)
    {
        uint32_t written = 0;
        CHECK(wirehair_encode(
            encoder, id, block.data(), block_bytes, &written) ==
            Wirehair_Success);
        decode_result = wirehair_decode(decoder, id, block.data(), written);
        CHECK(decode_result == Wirehair_NeedMore ||
            decode_result == Wirehair_Success);
        if (decode_result == Wirehair_Success) {
            break;
        }
    }
    CHECK(decode_result == Wirehair_Success);

    std::vector<uint8_t> recovered(expected.size());
    CHECK(wirehair_recover(
        decoder, recovered.data(), recovered.size()) == Wirehair_Success);
    CHECK(recovered == expected);
    return decoder;
}

void TestOwnedInputAndRecoverBlock()
{
    const uint32_t block_bytes = 31;
    const size_t message_bytes = block_bytes * 257u - 7u;
    const std::vector<uint8_t> expected = MakeMessage(message_bytes, 17);
    uint8_t* source = new uint8_t[message_bytes];
    std::memcpy(source, expected.data(), message_bytes);

    WirehairCodec encoder = nullptr;
    CHECK(wirehair_encoder_create_owned_ex(
        nullptr, source, message_bytes, block_bytes, &encoder) ==
        Wirehair_Success);
    CHECK(encoder != nullptr);
    std::memset(source, 0xa7, message_bytes);
    delete[] source;

    std::vector<uint8_t> systematic(block_bytes, uint8_t{0});
    uint32_t written = 0;
    CHECK(wirehair_encode(
        encoder, 0, systematic.data(), block_bytes, &written) ==
        Wirehair_Success);
    CHECK(written == block_bytes);
    CHECK(std::memcmp(systematic.data(), expected.data(), block_bytes) == 0);

    WirehairCodec decoder = DecodeLossy(encoder, expected, block_bytes);
    CHECK(decoder != nullptr);
    if (!decoder) {
        wirehair_free(encoder);
        return;
    }

    const unsigned block_count = 257;
    for (unsigned id = 0; id < block_count; ++id)
    {
        const uint32_t expected_bytes = id + 1 == block_count ?
            static_cast<uint32_t>(message_bytes % block_bytes) : block_bytes;
        std::vector<uint8_t> output(block_bytes + 2u, uint8_t{0xa5});
        uint32_t bytes = 1234;
        CHECK(wirehair_recover_block_ex(
            decoder, id, output.data() + 1, expected_bytes - 1, &bytes) ==
            Wirehair_InvalidInput);
        CHECK(bytes == 0);
        CHECK(std::all_of(output.begin(), output.end(),
            [](uint8_t value) { return value == 0xa5; }));

        CHECK(wirehair_recover_block_ex(
            decoder, id, output.data() + 1, expected_bytes, &bytes) ==
            Wirehair_Success);
        CHECK(bytes == expected_bytes);
        CHECK(output.front() == 0xa5 && output.back() == 0xa5);
        CHECK(std::memcmp(
            output.data() + 1,
            expected.data() + static_cast<size_t>(id) * block_bytes,
            expected_bytes) == 0);
    }

    std::vector<uint8_t> output(block_bytes, uint8_t{0xa5});
    uint32_t bytes = 99;
    CHECK(wirehair_recover_block_ex(
        decoder, 0, nullptr, 0, &bytes) == Wirehair_InvalidInput);
    CHECK(bytes == 0);
    CHECK(wirehair_recover_block_ex(
        decoder, 0x10000u, output.data(), block_bytes, &bytes) ==
        Wirehair_InvalidInput);
    CHECK(bytes == 0);
    CHECK(wirehair_recover_block_ex(
        decoder, 0, output.data(), block_bytes, nullptr) ==
        Wirehair_InvalidInput);
    CHECK(std::all_of(output.begin(), output.end(),
        [](uint8_t value) { return value == 0xa5; }));
    CHECK(wirehair_recover_block(
        decoder, 0, output.data(), &bytes) == Wirehair_Success);
    CHECK(bytes == block_bytes);

    CHECK(wirehair_decoder_create_ex(
        decoder, message_bytes, block_bytes, &decoder) == Wirehair_Success);
    WirehairResult reuse_result = Wirehair_NeedMore;
    for (unsigned id = block_count; id-- > 0;)
    {
        const uint32_t input_bytes = id + 1 == block_count ?
            static_cast<uint32_t>(message_bytes % block_bytes) : block_bytes;
        reuse_result = wirehair_decode(
            decoder,
            id,
            expected.data() + static_cast<size_t>(id) * block_bytes,
            input_bytes);
    }
    CHECK(reuse_result == Wirehair_Success);
    CHECK(wirehair_recover_block_ex(
        decoder, 0, output.data(), block_bytes, &bytes) == Wirehair_Success);
    CHECK(bytes == block_bytes);
    CHECK(std::memcmp(output.data(), expected.data(), block_bytes) == 0);

    std::vector<uint8_t> replacement = MakeMessage(message_bytes, 29);
    CHECK(wirehair_encoder_create_owned_ex(
        encoder, replacement.data(), replacement.size(), block_bytes,
        &encoder) == Wirehair_Success);
    std::fill(replacement.begin(), replacement.end(), uint8_t{0});
    CHECK(wirehair_encode(
        encoder, 0, output.data(), block_bytes, &bytes) == Wirehair_Success);
    const std::vector<uint8_t> replacement_expected =
        MakeMessage(message_bytes, 29);
    CHECK(std::memcmp(output.data(), replacement_expected.data(), block_bytes) == 0);

    wirehair_free(decoder);
    wirehair_free(encoder);
}

void TestLifecycleAndBorrowedMutation()
{
    const uint32_t block_bytes = 32;
    std::vector<uint8_t> message = MakeMessage(63, 3);
    WirehairCodec encoder = wirehair_encoder_create(
        nullptr, message.data(), message.size(), block_bytes);
    CHECK(encoder != nullptr);
    message[0] ^= 0xff;

    std::vector<uint8_t> output(block_bytes, uint8_t{0xa5});
    uint32_t bytes = 77;
    CHECK(wirehair_encode(
        encoder, 0, output.data(), block_bytes, &bytes) == Wirehair_Success);
    CHECK(bytes == block_bytes && output[0] == message[0]);
    CHECK(wirehair_decode(
        encoder, 0, output.data(), bytes) == Wirehair_InvalidInput);
    CHECK(wirehair_recover(
        encoder, output.data(), message.size()) == Wirehair_InvalidInput);
    CHECK(wirehair_recover_block_ex(
        encoder, 0, output.data(), block_bytes, &bytes) ==
        Wirehair_InvalidInput);
    CHECK(bytes == 0);
    CHECK(wirehair_decoder_becomes_encoder(encoder) == Wirehair_InvalidInput);

    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message.size(), block_bytes);
    CHECK(decoder != nullptr);
    std::fill(output.begin(), output.end(), uint8_t{0xa5});
    bytes = 77;
    CHECK(wirehair_encode(
        decoder, 0, output.data(), block_bytes, &bytes) ==
        Wirehair_InvalidInput);
    CHECK(bytes == 0);
    CHECK(std::all_of(output.begin(), output.end(),
        [](uint8_t value) { return value == 0xa5; }));
    CHECK(wirehair_recover(
        decoder, message.data(), message.size()) == Wirehair_NeedMore);
    CHECK(wirehair_recover_block_ex(
        decoder, 0, output.data(), block_bytes, &bytes) ==
        Wirehair_NeedMore);
    CHECK(bytes == 0);
    CHECK(wirehair_decoder_becomes_encoder(decoder) == Wirehair_InvalidInput);

    WirehairResult decode_result = Wirehair_NeedMore;
    for (unsigned id = 0; id < 2; ++id)
    {
        CHECK(wirehair_encode(
            encoder, id, output.data(), block_bytes, &bytes) ==
            Wirehair_Success);
        decode_result = wirehair_decode(decoder, id, output.data(), bytes);
    }
    CHECK(decode_result == Wirehair_Success);
    std::vector<uint8_t> recovered(message.size());
    CHECK(wirehair_recover(
        decoder, recovered.data(), recovered.size()) == Wirehair_Success);
    CHECK(recovered == message);

    bytes = 77;
    std::fill(output.begin(), output.end(), uint8_t{0xa5});
    CHECK(wirehair_encode(
        decoder, 0, output.data(), block_bytes, &bytes) ==
        Wirehair_InvalidInput);
    CHECK(bytes == 0);
    CHECK(wirehair_decoder_becomes_encoder(decoder) == Wirehair_Success);
    CHECK(wirehair_decoder_becomes_encoder(decoder) == Wirehair_InvalidInput);
    CHECK(wirehair_decode(
        decoder, 3, output.data(), block_bytes) == Wirehair_InvalidInput);
    CHECK(wirehair_recover(
        decoder, recovered.data(), recovered.size()) == Wirehair_InvalidInput);
    CHECK(wirehair_encode(
        decoder, 3, output.data(), block_bytes, &bytes) == Wirehair_Success);

    WirehairCodec failed = wirehair_decoder_create(
        nullptr, message.size(), block_bytes);
    CHECK(failed != nullptr);
    CHECK(wirehair_encode(
        encoder, 0, output.data(), block_bytes, &bytes) == Wirehair_Success);
    CHECK(wirehair_decode(failed, 0, output.data(), bytes) == Wirehair_NeedMore);
    CHECK(wirehair_decode(failed, 0, output.data(), bytes) ==
        Wirehair_InvalidInput);
    CHECK(wirehair_decode(failed, 1, output.data(), bytes) ==
        Wirehair_InvalidInput);
    CHECK(wirehair_recover(
        failed, recovered.data(), recovered.size()) == Wirehair_InvalidInput);
    CHECK(wirehair_decoder_becomes_encoder(failed) == Wirehair_InvalidInput);
    bytes = 77;
    CHECK(wirehair_encode(
        failed, 0, output.data(), block_bytes, &bytes) ==
        Wirehair_InvalidInput);
    CHECK(bytes == 0);

    WirehairCodec reused = nullptr;
    CHECK(wirehair_decoder_create_ex(
        encoder, message.size(), block_bytes, &reused) == Wirehair_Success);
    CHECK(reused != nullptr);
    bytes = 77;
    CHECK(wirehair_encode(
        reused, 0, output.data(), block_bytes, &bytes) ==
        Wirehair_InvalidInput);
    CHECK(bytes == 0);

    wirehair_free(reused);
    wirehair_free(failed);
    wirehair_free(decoder);
    wirehair_free(nullptr);
}

void TestHighNRecoveryIndex()
{
    const unsigned block_count = 64000;
    std::vector<uint8_t> message(block_count);
    for (unsigned i = 0; i < block_count; ++i) {
        message[i] = static_cast<uint8_t>((i * 17u + 5u) & 0xffu);
    }

    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message.size(), 1);
    CHECK(decoder != nullptr);
    if (!decoder) {
        return;
    }
    WirehairResult result = Wirehair_NeedMore;
    for (unsigned id = block_count; id-- > 0;) {
        result = wirehair_decode(decoder, id, &message[id], 1);
    }
    CHECK(result == Wirehair_Success);

    const auto start = std::chrono::steady_clock::now();
    uint64_t checksum = 0;
    for (unsigned repeat = 0; repeat < 3; ++repeat)
    {
        for (unsigned id = 0; id < block_count; ++id)
        {
            uint8_t value = 0;
            uint32_t bytes = 0;
            CHECK(wirehair_recover_block_ex(
                decoder, id, &value, 1, &bytes) == Wirehair_Success);
            CHECK(bytes == 1 && value == message[id]);
            checksum += value;
        }
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start).count();
    CHECK(checksum != 0);
    CHECK(elapsed < 300);
    std::printf(
        "N=64000 indexed recover_block: %lld ms for three full passes\n",
        static_cast<long long>(elapsed));
    wirehair_free(decoder);
}

void BenchmarkOwnedCreationCost()
{
    const uint32_t block_bytes = 1024;
    const size_t message_bytes = block_bytes * 256u - 13u;
    const std::vector<uint8_t> message = MakeMessage(message_bytes, 41);

    const auto measure = [&](bool owned) {
        const auto start = std::chrono::steady_clock::now();
        for (unsigned i = 0; i < 5; ++i)
        {
            WirehairCodec encoder = owned ?
                wirehair_encoder_create_owned(
                    nullptr, message.data(), message.size(), block_bytes) :
                wirehair_encoder_create(
                    nullptr, message.data(), message.size(), block_bytes);
            CHECK(encoder != nullptr);
            wirehair_free(encoder);
        }
        return std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - start).count();
    };

    const long long borrowed_us = measure(false);
    const long long owned_us = measure(true);
    std::printf(
        "Encoder creation (5 x %zu bytes): borrowed=%lld us owned=%lld us\n",
        message_bytes, borrowed_us, owned_us);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], "--preinit") == 0) {
        return RunPreinit();
    }
    if (argc == 2 && std::strcmp(argv[1], "--init-failure") == 0) {
        return RunInitFailure();
    }
    if (argc == 2 && std::strcmp(argv[1], "--cold-start") == 0) {
        return RunColdStart();
    }

    TestEnvironmentValue();
    TestCpuFeatureSelection();
    CHECK(wirehair_init_(WIREHAIR_VERSION + 1) == Wirehair_InvalidInput);
    CHECK(gf256_init_(GF256_VERSION + 1) == -1);
    CHECK(wirehair_init() == Wirehair_Success);
    CHECK(wirehair_init() == Wirehair_Success);
    gf256_x86_cpu_features active;
    gf256_get_active_x86_cpu_features(&active);
    CHECK((active.SSSE3 == 0 || active.SSSE3 == 1) &&
        (active.AVX2 == 0 || active.AVX2 == 1) &&
        (active.GFNI == 0 || active.GFNI == 1) &&
        (active.AVX512 == 0 || active.AVX512 == 1));
    std::printf("Active x86 kernels: SSSE3=%d AVX2=%d GFNI=%d AVX512=%d\n",
        active.SSSE3, active.AVX2, active.GFNI, active.AVX512);

    TestCreationResultsAndBoundaries();
    TestOwnedInputAndRecoverBlock();
    TestLifecycleAndBorrowedMutation();
    TestHighNRecoveryIndex();
    BenchmarkOwnedCreationCost();

    if (Failures != 0) {
        std::fprintf(stderr, "Legacy core test failed: %d checks\n", Failures);
        return 1;
    }
    std::puts("Legacy core test passed");
    return 0;
}
