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
#include <map>
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

uint64_t Fnv1aByte(uint64_t hash, uint8_t value)
{
    return (hash ^ value) * UINT64_C(1099511628211);
}

uint64_t Fnv1aU16(uint64_t hash, uint16_t value)
{
    hash = Fnv1aByte(hash, static_cast<uint8_t>(value));
    return Fnv1aByte(hash, static_cast<uint8_t>(value >> 8));
}

uint64_t Fnv1aU32(uint64_t hash, uint32_t value)
{
    hash = Fnv1aU16(hash, static_cast<uint16_t>(value));
    return Fnv1aU16(hash, static_cast<uint16_t>(value >> 16));
}

uint64_t TestOnlyMessageFingerprint(const std::vector<uint8_t>& message)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (uint8_t value : message) {
        hash = Fnv1aByte(hash, value);
    }
    return hash;
}

uint64_t DenseCountFingerprint()
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned n = 2; n <= 64000; ++n)
    {
        hash = Fnv1aU16(hash, static_cast<uint16_t>(n));
        hash = Fnv1aU16(hash, wirehair::GetDenseCount(n));
    }
    return hash;
}

void TestDenseCountInterpolation()
{
    using wirehair::detail::InterpolateDenseCount;

    // Ascending and endpoint cases.
    CHECK(InterpolateDenseCount(100, 200, 10, 30, 100) == 10);
    CHECK(InterpolateDenseCount(100, 200, 10, 30, 125) == 15);
    CHECK(InterpolateDenseCount(100, 200, 10, 30, 200) == 30);

    // A flat segment must remain flat throughout its domain.
    CHECK(InterpolateDenseCount(100, 200, 17, 17, 137) == 17);

    // The product and quotient are both negative on this descending segment.
    CHECK(InterpolateDenseCount(100, 200, 30, 10, 126) == 25);

    // Exercise a product outside signed 32-bit range and the uint16_t limit.
    CHECK(InterpolateDenseCount(
        0, 64000, UINT16_MAX, 1, 64000) == 1);

    const uint64_t fingerprint = DenseCountFingerprint();
    CHECK(fingerprint == UINT64_C(0x1a472addbc37c434));
    std::printf("Dense-count N=2..64000 fingerprint: %016llx\n",
        static_cast<unsigned long long>(fingerprint));
}

uint64_t SeedProfileFingerprint(bool base)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned n = 2; n <= 64000; ++n)
    {
        const uint16_t dense_count = wirehair::GetDenseCount(n);
        const uint16_t dense_seed = base ?
            wirehair::GetDenseSeedPreFixup(n, dense_count) :
            wirehair::GetDenseSeed(n, dense_count);
        const uint16_t peel_seed = base ?
            wirehair::GetPeelSeedPreFixup(n) :
            wirehair::GetPeelSeed(n);
        hash = Fnv1aU16(hash, static_cast<uint16_t>(n));
        hash = Fnv1aU16(hash, dense_count);
        hash = Fnv1aU16(hash, dense_seed);
        hash = Fnv1aU16(hash, peel_seed);
    }
    return hash;
}

uint64_t PacketSampleFingerprint(
    const WirehairWireProfile& profile,
    uint32_t block_count,
    uint32_t salt)
{
    const std::vector<uint8_t> message = MakeMessage(block_count, salt);
    WirehairCodec encoder = nullptr;
    CHECK(wirehair_encoder_create_profile_ex(
        nullptr, message.data(), message.size(), 1, &profile,
        WIREHAIR_ENCODER_OWN_INPUT, &encoder) == Wirehair_Success);
    if (!encoder) {
        return 0;
    }

    const uint32_t ids[] = {
        block_count,
        block_count + 1,
        block_count * 2 + 17,
        UINT32_C(0x10203040)
    };
    uint64_t hash = UINT64_C(14695981039346656037);
    for (uint32_t id : ids)
    {
        uint8_t packet = 0;
        uint32_t written = 0;
        CHECK(wirehair_encode(
            encoder, id, &packet, 1, &written) == Wirehair_Success);
        CHECK(written == 1);
        hash = Fnv1aU32(hash, id);
        hash = Fnv1aByte(hash, packet);
    }
    wirehair_free(encoder);
    return hash;
}

uint64_t PacketCorpusFingerprint(const WirehairWireProfile& profile)
{
    std::vector<uint32_t> counts = {
        2, 3, 4, 7, 8, 15, 16, 17, 31, 32, 63, 64, 65,
        127, 128, 129, 255, 256, 257, 511, 512, 513,
        1023, 1024, 1025, 2047, 2048, 2049, 2373, 4095,
        4096, 4097, 5550, 8191, 8192, 8193, 16383, 16384,
        16385, 31999, 32000, 32001, 63998, 63999, 64000
    };
    uint64_t state = UINT64_C(0x853c49e6748fea9b);
    for (unsigned i = 0; i < 256; ++i)
    {
        state = state * UINT64_C(6364136223846793005) +
            UINT64_C(1442695040888963407);
        counts.push_back(2 + static_cast<uint32_t>(state % 63999));
    }
    std::sort(counts.begin(), counts.end());
    counts.erase(std::unique(counts.begin(), counts.end()), counts.end());

    uint64_t hash = UINT64_C(14695981039346656037);
    for (uint32_t count : counts)
    {
        const std::vector<uint8_t> message = MakeMessage(count, count);
        WirehairCodec encoder = nullptr;
        CHECK(wirehair_encoder_create_profile_ex(
            nullptr, message.data(), message.size(), 1, &profile,
            WIREHAIR_ENCODER_OWN_INPUT, &encoder) == Wirehair_Success);
        if (!encoder) {
            return 0;
        }

        const uint32_t ids[] = {
            count, count + 1, count * 2 + 17,
            UINT32_C(0x10203040), UINT32_MAX
        };
        hash = Fnv1aU32(hash, count);
        for (uint32_t id : ids)
        {
            uint8_t packet = 0;
            uint32_t written = 0;
            CHECK(wirehair_encode(
                encoder, id, &packet, 1, &written) == Wirehair_Success);
            CHECK(written == 1);
            hash = Fnv1aU32(hash, id);
            hash = Fnv1aByte(hash, packet);
        }
        wirehair_free(encoder);
    }
    return hash;
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

void TestGF256WireContract()
{
    // These values pin the production GF(256) field used by the heavy rows.
    // They are interoperability vectors, not a replacement for the exhaustive
    // backend self-tests in gf256_init().
    CHECK(GF256Ctx.Polynomial == 0x14d);
    CHECK(gf256_mul(0x53, 0xca) == 0x94);
    CHECK(gf256_mul(0x57, 0x83) == 0x43);
    CHECK(gf256_div(0x53, 0xca) == 0x44);
    CHECK(gf256_inv(0xca) == 0x5d);

    // Pin an externally observable packet too.  This catches changes to the
    // heavy-field contract even when the internal arithmetic vectors still
    // happen to be updated together with the implementation.
    static const uint8_t kExpectedRepair[16] = {
        0x17, 0xed, 0xd1, 0x3d, 0xc5, 0xa4, 0x7a, 0x85,
        0xc4, 0x43, 0xc6, 0x57, 0x64, 0xe8, 0x7e, 0x56
    };
    const std::vector<uint8_t> message = MakeMessage(128);
    WirehairCodec encoder = wirehair_encoder_create(
        nullptr, message.data(), message.size(), 16);
    CHECK(encoder != nullptr);
    if (encoder)
    {
        uint8_t repair[sizeof(kExpectedRepair)] = {};
        uint32_t written = 0;
        CHECK(wirehair_encode(
            encoder, 12345, repair, sizeof(repair), &written) ==
            Wirehair_Success);
        CHECK(written == sizeof(kExpectedRepair));
        CHECK(std::memcmp(repair, kExpectedRepair, sizeof(repair)) == 0);
        wirehair_free(encoder);
    }
}

void TestLegacyWireProfiles()
{
    CHECK(sizeof(WirehairWireProfile) == 16);
    CHECK(WIREHAIR_LEGACY_PROFILE_PRE_FIXUP ==
        UINT64_C(0xe1b9f77f1c90f680));
    CHECK(WIREHAIR_LEGACY_PROFILE_FIXUPS_2026_07 ==
        UINT64_C(0x4d241359db07bb07));
    CHECK(WIREHAIR_LEGACY_PROFILE_CURRENT ==
        WIREHAIR_LEGACY_PROFILE_FIXUPS_2026_07);
    CHECK(WIREHAIR_LEGACY_PROFILE_PRE_FIXUP !=
        WIREHAIR_LEGACY_PROFILE_CURRENT);

    WirehairWireProfile base = {};
    WirehairWireProfile current = {};
    CHECK(wirehair_wire_profile_init(
        WIREHAIR_LEGACY_PROFILE_PRE_FIXUP, &base) == Wirehair_Success);
    CHECK(wirehair_wire_profile_init(
        WIREHAIR_LEGACY_PROFILE_CURRENT, &current) == Wirehair_Success);
    CHECK(base.struct_bytes == sizeof(base));
    CHECK(base.profile_version == WIREHAIR_WIRE_PROFILE_VERSION);
    CHECK(base.profile_id == WIREHAIR_LEGACY_PROFILE_PRE_FIXUP);
    CHECK(current.struct_bytes == sizeof(current));
    CHECK(current.profile_version == WIREHAIR_WIRE_PROFILE_VERSION);
    CHECK(current.profile_id == WIREHAIR_LEGACY_PROFILE_CURRENT);
    CHECK(wirehair_wire_profile_init(0, nullptr) == Wirehair_InvalidInput);

    WirehairWireProfile unknown = { 1, 1, 1 };
    CHECK(wirehair_wire_profile_init(UINT64_C(0x123456789abcdef0), &unknown) ==
        Wirehair_InvalidInput);
    CHECK(unknown.struct_bytes == 0 && unknown.profile_version == 0 &&
        unknown.profile_id == 0);

    uint8_t tiny_message[2] = { 1, 2 };
    WirehairCodec rejected = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    WirehairWireProfile malformed = current;
    --malformed.struct_bytes;
    CHECK(wirehair_encoder_create_profile_ex(
        nullptr, tiny_message, sizeof(tiny_message), 1, &malformed, 0,
        &rejected) == Wirehair_InvalidInput);
    CHECK(rejected == nullptr);
    malformed = current;
    ++malformed.profile_version;
    rejected = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_decoder_create_profile_ex(
        nullptr, sizeof(tiny_message), 1, &malformed, &rejected) ==
        Wirehair_InvalidInput);
    CHECK(rejected == nullptr);
    malformed = current;
    malformed.profile_id = UINT64_C(0x123456789abcdef0);
    rejected = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_decoder_create_profile_ex(
        nullptr, sizeof(tiny_message), 1, &malformed, &rejected) ==
        Wirehair_InvalidInput);
    CHECK(rejected == nullptr);
    rejected = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_encoder_create_profile_ex(
        nullptr, tiny_message, sizeof(tiny_message), 1, &current, 2,
        &rejected) == Wirehair_InvalidInput);
    CHECK(rejected == nullptr);

    const uint64_t base_seed_hash = SeedProfileFingerprint(true);
    const uint64_t current_seed_hash = SeedProfileFingerprint(false);
    CHECK(base_seed_hash == UINT64_C(0xea803c84c587877a));
    CHECK(current_seed_hash == UINT64_C(0x4294485edd1cc4dd));
    std::printf("Legacy profile seed hashes: pre-fixup=%016llx current=%016llx\n",
        static_cast<unsigned long long>(base_seed_hash),
        static_cast<unsigned long long>(current_seed_hash));
    const uint64_t base_n5550 = PacketSampleFingerprint(base, 5550, 91);
    const uint64_t current_n5550 = PacketSampleFingerprint(current, 5550, 91);
    CHECK(base_n5550 == UINT64_C(0x18497889fd2e5fd7));
    CHECK(current_n5550 == UINT64_C(0xa87074dc238666c0));
    std::printf("Legacy N=5550 sample hashes: pre-fixup=%016llx current=%016llx\n",
        static_cast<unsigned long long>(base_n5550),
        static_cast<unsigned long long>(current_n5550));
    const uint64_t base_corpus = PacketCorpusFingerprint(base);
    const uint64_t current_corpus = PacketCorpusFingerprint(current);
    CHECK(base_corpus == UINT64_C(0xd69a6563d81cf247));
    CHECK(current_corpus == UINT64_C(0x26898d9bc72d6f10));
    std::printf("Legacy 301-N packet corpus: pre-fixup=%016llx current=%016llx\n",
        static_cast<unsigned long long>(base_corpus),
        static_cast<unsigned long long>(current_corpus));

    // The historical anchor encoder accidentally compared only the low 16
    // bits of a repair ID when deciding whether to truncate the final block.
    // Both named profiles follow the anchor decoder's full-width contract.
    static const uint8_t kShortFinalRepair[7] = {
        0xa9, 0x4c, 0x4d, 0x5e, 0xcd, 0x52, 0x31
    };
    const std::vector<uint8_t> short_message = MakeMessage(53, 123);
    const WirehairWireProfile* short_profiles[] = { &base, &current };
    for (const WirehairWireProfile* profile : short_profiles)
    {
        WirehairCodec encoder = nullptr;
        CHECK(wirehair_encoder_create_profile_ex(
            nullptr, short_message.data(), short_message.size(), 7, profile,
            WIREHAIR_ENCODER_OWN_INPUT, &encoder) == Wirehair_Success);
        if (!encoder) {
            continue;
        }
        uint8_t packet[7] = {};
        uint32_t written = 0;
        CHECK(wirehair_encode(
            encoder, UINT32_C(65543), packet, sizeof(packet), &written) ==
            Wirehair_Success);
        CHECK(written == sizeof(packet));
        CHECK(std::memcmp(
            packet, kShortFinalRepair, sizeof(kShortFinalRepair)) == 0);

        written = 0;
        CHECK(wirehair_encode(
            encoder, 7, packet, sizeof(packet), &written) == Wirehair_Success);
        CHECK(written == 4);
        CHECK(std::memcmp(packet, short_message.data() + 49, written) == 0);
        wirehair_free(encoder);
    }

    // N=2373 is the first exact-N peel fixup.  A repair-only stream from the
    // pre-fixup profile is solvable by the current raw decoder but reconstructs
    // a different message, demonstrating why the trusted profile id and an
    // application integrity check are both required.
    const uint32_t block_count = 2373;
    const std::vector<uint8_t> message = MakeMessage(block_count, 77);
    WirehairCodec base_encoder = nullptr;
    WirehairCodec base_decoder = nullptr;
    WirehairCodec raw_mismatched_decoder = nullptr;
    CHECK(wirehair_encoder_create_profile_ex(
        nullptr, message.data(), message.size(), 1, &base,
        WIREHAIR_ENCODER_OWN_INPUT, &base_encoder) == Wirehair_Success);
    CHECK(wirehair_decoder_create_profile_ex(
        nullptr, message.size(), 1, &base, &base_decoder) == Wirehair_Success);
    // Model a receiver that did not carry the profile in trusted framing and
    // therefore used the raw API's frozen current profile.
    raw_mismatched_decoder = wirehair_decoder_create(
        nullptr, message.size(), 1);
    CHECK(raw_mismatched_decoder != nullptr);

    uint64_t base_packet_hash = UINT64_C(14695981039346656037);
    WirehairResult base_result = Wirehair_NeedMore;
    WirehairResult mismatch_result = Wirehair_NeedMore;
    for (uint32_t i = 0; i < block_count; ++i)
    {
        const uint32_t id = block_count + i;
        uint8_t packet = 0;
        uint32_t written = 0;
        CHECK(wirehair_encode(
            base_encoder, id, &packet, 1, &written) == Wirehair_Success);
        CHECK(written == 1);
        base_packet_hash = Fnv1aU32(base_packet_hash, id);
        base_packet_hash = Fnv1aByte(base_packet_hash, packet);
        base_result = wirehair_decode(base_decoder, id, &packet, written);
        mismatch_result = wirehair_decode(
            raw_mismatched_decoder, id, &packet, written);
    }
    CHECK(base_result == Wirehair_Success);
    CHECK(mismatch_result == Wirehair_Success);
    std::vector<uint8_t> base_recovered(message.size());
    std::vector<uint8_t> mismatch_recovered(message.size());
    CHECK(wirehair_recover(
        base_decoder, base_recovered.data(), base_recovered.size()) ==
        Wirehair_Success);
    CHECK(wirehair_recover(
        raw_mismatched_decoder, mismatch_recovered.data(),
        mismatch_recovered.size()) == Wirehair_Success);
    CHECK(base_recovered == message);
    CHECK(mismatch_recovered != message);
    // FNV is only a deterministic test stand-in.  Production callers need a
    // cryptographic digest/MAC obtained independently from trusted metadata.
    const uint64_t trusted_fingerprint = TestOnlyMessageFingerprint(message);
    CHECK(TestOnlyMessageFingerprint(base_recovered) == trusted_fingerprint);
    CHECK(TestOnlyMessageFingerprint(mismatch_recovered) !=
        trusted_fingerprint);
    CHECK(base_packet_hash == UINT64_C(0xa6c1a4d75faec4d8));
    std::printf("Legacy pre-fixup N=2373 repair hash: %016llx\n",
        static_cast<unsigned long long>(base_packet_hash));

    // Decoder-to-encoder conversion must retain the selected profile.
    CHECK(wirehair_decoder_becomes_encoder(base_decoder) == Wirehair_Success);
    uint8_t converted_packet = 0, base_packet = 0;
    uint32_t converted_written = 0, base_written = 0;
    CHECK(wirehair_encode(base_decoder, block_count * 2 + 31,
        &converted_packet, 1, &converted_written) == Wirehair_Success);
    CHECK(wirehair_encode(base_encoder, block_count * 2 + 31,
        &base_packet, 1, &base_written) == Wirehair_Success);
    CHECK(converted_written == 1 && base_written == 1 &&
        converted_packet == base_packet);

    // The unframed API is deliberately fixed to the current profile.
    WirehairCodec current_encoder = nullptr;
    CHECK(wirehair_encoder_create_profile_ex(
        nullptr, message.data(), message.size(), 1, &current, 0,
        &current_encoder) == Wirehair_Success);
    WirehairCodec raw_encoder = wirehair_encoder_create(
        nullptr, message.data(), message.size(), 1);
    CHECK(raw_encoder != nullptr);
    WirehairCodec current_decoder = nullptr;
    CHECK(wirehair_decoder_create_profile_ex(
        nullptr, message.size(), 1, &current, &current_decoder) ==
        Wirehair_Success);
    WirehairCodec raw_decoder = wirehair_decoder_create(
        nullptr, message.size(), 1);
    CHECK(raw_decoder != nullptr);
    uint64_t current_packet_hash = UINT64_C(14695981039346656037);
    WirehairResult current_result = Wirehair_NeedMore;
    WirehairResult raw_result = Wirehair_NeedMore;
    for (uint32_t i = 0; i < block_count; ++i)
    {
        const uint32_t id = block_count + i;
        uint8_t explicit_packet = 0, raw_packet = 0;
        uint32_t explicit_written = 0, raw_written = 0;
        CHECK(wirehair_encode(current_encoder, id, &explicit_packet, 1,
            &explicit_written) == Wirehair_Success);
        CHECK(wirehair_encode(raw_encoder, id, &raw_packet, 1,
            &raw_written) == Wirehair_Success);
        CHECK(explicit_written == 1 && raw_written == 1);
        CHECK(explicit_packet == raw_packet);
        current_packet_hash = Fnv1aU32(current_packet_hash, id);
        current_packet_hash = Fnv1aByte(current_packet_hash, explicit_packet);
        current_result = wirehair_decode(
            current_decoder, id, &explicit_packet, explicit_written);
        raw_result = wirehair_decode(
            raw_decoder, id, &raw_packet, raw_written);
    }
    CHECK(current_packet_hash == UINT64_C(0xa9020b085ecd5dfa));
    CHECK(current_result == Wirehair_Success);
    CHECK(raw_result == Wirehair_Success);
    std::vector<uint8_t> current_recovered(message.size());
    std::vector<uint8_t> raw_recovered(message.size());
    CHECK(wirehair_recover(current_decoder, current_recovered.data(),
        current_recovered.size()) == Wirehair_Success);
    CHECK(wirehair_recover(raw_decoder, raw_recovered.data(),
        raw_recovered.size()) == Wirehair_Success);
    CHECK(current_recovered == message);
    CHECK(raw_recovered == message);
    std::printf("Legacy current N=2373 repair hash: %016llx\n",
        static_cast<unsigned long long>(current_packet_hash));

    // Reuse must replace, rather than retain, the prior equation profile.
    WirehairCodec reused_encoder = nullptr;
    CHECK(wirehair_encoder_create_profile_ex(
        base_encoder, message.data(), message.size(), 1, &current,
        WIREHAIR_ENCODER_OWN_INPUT, &reused_encoder) == Wirehair_Success);
    base_encoder = nullptr;
    uint8_t reused_packet = 0, expected_packet = 0;
    uint32_t reused_written = 0, expected_written = 0;
    CHECK(wirehair_encode(reused_encoder, block_count + 19,
        &reused_packet, 1, &reused_written) == Wirehair_Success);
    CHECK(wirehair_encode(current_encoder, block_count + 19,
        &expected_packet, 1, &expected_written) == Wirehair_Success);
    CHECK(reused_written == 1 && expected_written == 1 &&
        reused_packet == expected_packet);

    WirehairCodec reused_decoder = nullptr;
    CHECK(wirehair_decoder_create_profile_ex(
        base_decoder, message.size(), 1, &current, &reused_decoder) ==
        Wirehair_Success);
    base_decoder = nullptr;
    WirehairResult reused_result = Wirehair_NeedMore;
    for (uint32_t i = 0; i < block_count; ++i)
    {
        const uint32_t id = block_count + i;
        uint8_t packet = 0;
        uint32_t written = 0;
        CHECK(wirehair_encode(current_encoder, id, &packet, 1, &written) ==
            Wirehair_Success);
        reused_result = wirehair_decode(reused_decoder, id, &packet, written);
    }
    CHECK(reused_result == Wirehair_Success);
    std::vector<uint8_t> reused_recovered(message.size());
    CHECK(wirehair_recover(reused_decoder, reused_recovered.data(),
        reused_recovered.size()) == Wirehair_Success);
    CHECK(reused_recovered == message);

    wirehair_free(base_encoder);
    wirehair_free(base_decoder);
    wirehair_free(raw_mismatched_decoder);
    wirehair_free(current_encoder);
    wirehair_free(raw_encoder);
    wirehair_free(current_decoder);
    wirehair_free(raw_decoder);
    wirehair_free(reused_encoder);
    wirehair_free(reused_decoder);
}

int RunPreinit()
{
    uint8_t message[64] = {};
    WirehairWireProfile profile = {};
    CHECK(wirehair_wire_profile_init(
        WIREHAIR_LEGACY_PROFILE_CURRENT, &profile) == Wirehair_Success);
    CHECK(profile.profile_id == WIREHAIR_LEGACY_PROFILE_CURRENT);
    WirehairCodec codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_encoder_create_ex(
        nullptr, message, sizeof(message), 32, &codec) == Wirehair_Error);
    CHECK(codec == nullptr);
    codec = reinterpret_cast<WirehairCodec>(uintptr_t(1));
    CHECK(wirehair_encoder_create_profile_ex(
        nullptr, message, sizeof(message), 32, &profile, 0, &codec) ==
        Wirehair_Error);
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

    WirehairCodec duplicate_decoder = wirehair_decoder_create(
        nullptr, message.size(), block_bytes);
    CHECK(duplicate_decoder != nullptr);
    CHECK(wirehair_encode(
        encoder, 0, output.data(), block_bytes, &bytes) == Wirehair_Success);
    CHECK(wirehair_decode(
        duplicate_decoder, 0, output.data(), bytes) == Wirehair_NeedMore);
    CHECK(wirehair_decode(
        duplicate_decoder, 0, output.data(), bytes) == Wirehair_NeedMore);
    CHECK(wirehair_encode(
        encoder, 1, output.data(), block_bytes, &bytes) == Wirehair_Success);
    CHECK(wirehair_decode(
        duplicate_decoder, 1, output.data(), bytes) == Wirehair_Success);
    CHECK(wirehair_recover(
        duplicate_decoder, recovered.data(), recovered.size()) ==
        Wirehair_Success);
    CHECK(recovered == message);
    CHECK(wirehair_decoder_becomes_encoder(duplicate_decoder) ==
        Wirehair_Success);
    bytes = 77;
    CHECK(wirehair_encode(
        duplicate_decoder, 0, output.data(), block_bytes, &bytes) ==
        Wirehair_Success);
    CHECK(bytes == block_bytes);

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
    wirehair_free(duplicate_decoder);
    wirehair_free(decoder);
    wirehair_free(nullptr);
}

void TestPacketFingerprintVectors()
{
    // Independently generated with a compact reference SipHash-2-4
    // implementation.  These pin keys, rounds, byte order, tail packing, and
    // both halves of the fixed duplicate-equality contract.
    struct Vector
    {
        uint32_t Bytes;
        uint64_t Hash0;
        uint64_t Hash1;
    };
    static const Vector vectors[] = {
        {   0, UINT64_C(0x751bf6774078774a), UINT64_C(0xf591d34bdde27228) },
        {   1, UINT64_C(0x2447d1d7f46786aa), UINT64_C(0xd415a66f30deb419) },
        {   7, UINT64_C(0x223c277ec6187d37), UINT64_C(0x31ac659a8c7ebd32) },
        {   8, UINT64_C(0x63a02fa1962d9e59), UINT64_C(0xd160a0d3392fe438) },
        {   9, UINT64_C(0x4be6afd93cfbb34d), UINT64_C(0x1274f83d0918ff39) },
        {  15, UINT64_C(0xed73050006ce018d), UINT64_C(0xfe836b9cc5cb8a44) },
        {  16, UINT64_C(0x846186dbdd674086), UINT64_C(0x062cd0618c4507f9) },
        {  31, UINT64_C(0x2e09bdda950ba278), UINT64_C(0xf4135c3d74122e53) },
        {  63, UINT64_C(0xd8e3e83994fd5d43), UINT64_C(0x521cc2869b7a8186) },
        {  64, UINT64_C(0x28f561b0b4bf7938), UINT64_C(0x61fd47bf6a3e9b8a) },
        {  65, UINT64_C(0x54e3457ba238d8b9), UINT64_C(0x9c481d1302a967fa) },
        { 257, UINT64_C(0x1983e075a1af665e), UINT64_C(0xe8b76cac6ac432f2) }
    };
    uint8_t data[257];
    for (uint32_t i = 0; i < sizeof(data); ++i) {
        data[i] = static_cast<uint8_t>(i * 37u + 11u);
    }
    for (const Vector& vector : vectors)
    {
        uint64_t hash0 = 0;
        uint64_t hash1 = 0;
        wirehair::Codec::TestingPacketFingerprint128(
            data, vector.Bytes, hash0, hash1);
        CHECK(hash0 == vector.Hash0);
        CHECK(hash1 == vector.Hash1);
    }
}

void TestDuplicatePacketIdentity()
{
    const uint32_t block_bytes = 17;
    const uint32_t block_count = 8;
    const size_t message_bytes = block_bytes * block_count - 5;
    const std::vector<uint8_t> message = MakeMessage(message_bytes, 53);
    WirehairCodec encoder = wirehair_encoder_create(
        nullptr, message.data(), message.size(), block_bytes);
    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message.size(), block_bytes);
    CHECK(encoder != nullptr && decoder != nullptr);
    if (!encoder || !decoder) {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return;
    }

    const wirehair::Codec* internal =
        reinterpret_cast<const wirehair::Codec*>(decoder);
    const uint32_t expected_limit = block_count + 1024;
    uint32_t expected_slots = 1;
    while (expected_slots < expected_limit * 2) {
        expected_slots <<= 1;
    }
    CHECK(internal->TestingReceivedPacketLimit() == expected_limit);
    CHECK(internal->TestingReceivedPacketSlots() == expected_slots);
    CHECK(internal->TestingReceivedPacketMemoryBytes() ==
        static_cast<uint64_t>(expected_slots) * 24);

    std::vector<uint8_t> packet(block_bytes, uint8_t{0xa5});
    uint32_t written = 0;
    CHECK(wirehair_encode(
        encoder, 0, packet.data(), block_bytes, &written) ==
        Wirehair_Success);
    CHECK(written == block_bytes);
    const std::vector<uint8_t> systematic_packet = packet;
    CHECK(wirehair_decode(
        decoder, 0, packet.data(), written) == Wirehair_NeedMore);
    CHECK(internal->TestingDecoderRowCount() == 1);
    CHECK(internal->TestingReceivedPacketCount() == 1);

    for (unsigned repeat = 0; repeat < 512; ++repeat) {
        CHECK(wirehair_decode(
            decoder, 0, systematic_packet.data(), written) ==
            Wirehair_NeedMore);
    }
    CHECK(internal->TestingDecoderRowCount() == 1);
    CHECK(internal->TestingReceivedPacketCount() == 1);

    std::vector<uint8_t> conflict = systematic_packet;
    conflict[3] ^= 0x80;
    CHECK(wirehair_decode(
        decoder, 0, conflict.data(), written) == Wirehair_InvalidInput);
    CHECK(wirehair_decode(
        decoder, 0, conflict.data(), written) == Wirehair_InvalidInput);
    CHECK(internal->TestingDecoderRowCount() == 1);
    CHECK(internal->TestingReceivedPacketCount() == 1);

    uint32_t repair_id = block_count + 17;
    std::vector<uint8_t> repair_packet;
    const uint32_t before_repair_count =
        internal->TestingReceivedPacketCount();
    for (unsigned attempt = 0; attempt < 64; ++attempt, ++repair_id)
    {
        packet.assign(block_bytes, uint8_t{0});
        CHECK(wirehair_encode(
            encoder, repair_id, packet.data(), block_bytes, &written) ==
            Wirehair_Success);
        CHECK(wirehair_decode(
            decoder, repair_id, packet.data(), written) == Wirehair_NeedMore);
        if (internal->TestingReceivedPacketCount() > before_repair_count) {
            repair_packet = packet;
            break;
        }
    }
    CHECK(!repair_packet.empty());
    if (!repair_packet.empty())
    {
        CHECK(wirehair_decode(
            decoder, repair_id, repair_packet.data(), block_bytes) ==
            Wirehair_NeedMore);
        conflict = repair_packet;
        conflict[0] ^= 1;
        CHECK(wirehair_decode(
            decoder, repair_id, conflict.data(), block_bytes) ==
            Wirehair_InvalidInput);
    }

    const uint32_t final_id = block_count - 1;
    packet.assign(block_bytes, uint8_t{0x11});
    CHECK(wirehair_encode(
        encoder, final_id, packet.data(), block_bytes, &written) ==
        Wirehair_Success);
    CHECK(written == block_bytes - 5);
    std::vector<uint8_t> final_packet = packet;
    CHECK(wirehair_decode(
        decoder, final_id, final_packet.data(), written) ==
        Wirehair_NeedMore);
    std::fill(
        final_packet.begin() + written, final_packet.end(), uint8_t{0xee});
    CHECK(wirehair_decode(
        decoder, final_id, final_packet.data(), block_bytes) ==
        Wirehair_NeedMore);
    conflict = final_packet;
    conflict[written - 1] ^= 1;
    CHECK(wirehair_decode(
        decoder, final_id, conflict.data(), block_bytes) ==
        Wirehair_InvalidInput);

    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t id = 1; id + 1 < block_count && result != Wirehair_Success;
         ++id)
    {
        packet.assign(block_bytes, uint8_t{0});
        CHECK(wirehair_encode(
            encoder, id, packet.data(), block_bytes, &written) ==
            Wirehair_Success);
        result = wirehair_decode(decoder, id, packet.data(), written);
        CHECK(result == Wirehair_NeedMore || result == Wirehair_Success);
    }
    for (uint32_t id = 1000;
         result == Wirehair_NeedMore && id < 1200;
         ++id)
    {
        packet.assign(block_bytes, uint8_t{0});
        CHECK(wirehair_encode(
            encoder, id, packet.data(), block_bytes, &written) ==
            Wirehair_Success);
        result = wirehair_decode(decoder, id, packet.data(), written);
        CHECK(result == Wirehair_NeedMore || result == Wirehair_Success);
    }
    CHECK(result == Wirehair_Success);

    std::vector<uint8_t> recovered(message.size());
    CHECK(wirehair_recover(
        decoder, recovered.data(), recovered.size()) == Wirehair_Success);
    CHECK(recovered == message);
    const uint32_t completed_rows = internal->TestingDecoderRowCount();
    const uint32_t completed_rank = internal->TestingDecoderRank();

    CHECK(wirehair_decode(
        decoder, 0, systematic_packet.data(), block_bytes) ==
        Wirehair_Success);
    std::vector<uint8_t> completed_conflict = systematic_packet;
    completed_conflict[3] ^= 0x80;
    CHECK(wirehair_decode(
        decoder, 0, completed_conflict.data(), block_bytes) ==
        Wirehair_InvalidInput);
    CHECK(internal->TestingDecoderRowCount() == completed_rows);
    CHECK(internal->TestingDecoderRank() == completed_rank);

    std::vector<uint8_t> max_id_packet(block_bytes, uint8_t{0x5a});
    CHECK(wirehair_decode(
        decoder, UINT32_MAX, max_id_packet.data(), block_bytes) ==
        Wirehair_Success);
    CHECK(wirehair_decode(
        decoder, UINT32_MAX, max_id_packet.data(), block_bytes) ==
        Wirehair_Success);
    max_id_packet[0] ^= 1;
    CHECK(wirehair_decode(
        decoder, UINT32_MAX, max_id_packet.data(), block_bytes) ==
        Wirehair_InvalidInput);

    uint32_t last_id = 0;
    std::vector<uint8_t> last_packet(block_bytes);
    for (uint32_t sequence = 0;
         internal->TestingReceivedPacketCount() < expected_limit;
         ++sequence)
    {
        last_id = UINT32_C(0x40000000) + sequence;
        for (uint32_t i = 0; i < block_bytes; ++i) {
            last_packet[i] = static_cast<uint8_t>(
                last_id * 29u + i * 71u + 3u);
        }
        CHECK(wirehair_decode(
            decoder, last_id, last_packet.data(), block_bytes) ==
            Wirehair_Success);
    }
    CHECK(internal->TestingReceivedPacketCount() == expected_limit);
    CHECK(internal->TestingDecoderRowCount() == completed_rows);
    CHECK(internal->TestingDecoderRank() == completed_rank);

    CHECK(wirehair_decode(
        decoder, last_id, last_packet.data(), block_bytes) ==
        Wirehair_Success);
    conflict = last_packet;
    conflict[5] ^= 1;
    CHECK(wirehair_decode(
        decoder, last_id, conflict.data(), block_bytes) ==
        Wirehair_InvalidInput);
    CHECK(wirehair_decode(
        decoder, UINT32_C(0x70000000), last_packet.data(), block_bytes) ==
        Wirehair_ExtraInsufficient);
    CHECK(wirehair_decode(
        decoder, UINT32_C(0x70000000), last_packet.data(), block_bytes) ==
        Wirehair_ExtraInsufficient);
    CHECK(wirehair_decode(
        decoder, 0, systematic_packet.data(), block_bytes) ==
        Wirehair_Success);
    CHECK(internal->TestingReceivedPacketCount() == expected_limit);
    CHECK(internal->TestingDecoderRowCount() == completed_rows);
    CHECK(internal->TestingDecoderRank() == completed_rank);
    CHECK(wirehair_recover(
        decoder, recovered.data(), recovered.size()) == Wirehair_Success);
    CHECK(recovered == message);

    WirehairCodec reused = nullptr;
    CHECK(wirehair_decoder_create_ex(
        decoder, message.size(), block_bytes, &reused) == Wirehair_Success);
    decoder = nullptr;
    CHECK(reused != nullptr);
    if (reused)
    {
        const wirehair::Codec* reused_internal =
            reinterpret_cast<const wirehair::Codec*>(reused);
        CHECK(reused_internal->TestingReceivedPacketCount() == 0);
        packet.assign(block_bytes, uint8_t{0});
        CHECK(wirehair_encode(
            encoder, 0, packet.data(), block_bytes, &written) ==
            Wirehair_Success);
        CHECK(wirehair_decode(
            reused, 0, packet.data(), written) == Wirehair_NeedMore);
        CHECK(reused_internal->TestingReceivedPacketCount() == 1);
    }

    wirehair_free(reused);
    wirehair_free(decoder);
    wirehair_free(encoder);
}

void TestDuplicatePacketDuringExtraResume()
{
    const uint32_t block_bytes = 23;
    const uint32_t block_count = 32;
    const std::vector<uint8_t> message =
        MakeMessage(block_bytes * block_count, 67);
    WirehairCodec encoder = wirehair_encoder_create(
        nullptr, message.data(), message.size(), block_bytes);
    CHECK(encoder != nullptr);
    if (!encoder) {
        return;
    }

    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message.size(), block_bytes);
    CHECK(decoder != nullptr);
    if (!decoder) {
        wirehair_free(encoder);
        return;
    }
    const wirehair::Codec* internal =
        reinterpret_cast<const wirehair::Codec*>(decoder);

    // Find two distinct repair IDs that generate exactly the same equation.
    // Feeding both makes the first N-row solve rank-deficient by construction,
    // so the next unique packet necessarily exercises extra-row resume.
    typedef std::pair<uint64_t, uint16_t> ParameterKey;
    std::map<ParameterKey, uint32_t> seen;
    uint32_t collision_a = 0;
    uint32_t collision_b = 0;
    const uint16_t mix_count =
        static_cast<uint16_t>(wirehair::GetDenseCount(block_count) + 6);
    for (uint32_t id = block_count; id < 200000 && collision_a == 0; ++id)
    {
        wirehair::PeelRowParameters params;
        params.Initialize(
            id, internal->PSeed(), block_count, mix_count);
        const uint64_t packed =
            static_cast<uint64_t>(params.PeelCount) |
            (static_cast<uint64_t>(params.PeelFirst) << 16) |
            (static_cast<uint64_t>(params.PeelAdd) << 32) |
            (static_cast<uint64_t>(params.MixFirst) << 48);
        const ParameterKey key(packed, params.MixAdd);
        const std::map<ParameterKey, uint32_t>::const_iterator prior =
            seen.find(key);
        if (prior != seen.end()) {
            collision_a = prior->second;
            collision_b = id;
        }
        else {
            seen.insert(std::make_pair(key, id));
        }
    }
    CHECK(collision_a >= block_count && collision_b > collision_a);

    std::vector<uint8_t> packet(block_bytes);
    std::vector<uint8_t> collision_packet_a(block_bytes);
    std::vector<uint8_t> collision_packet_b(block_bytes);
    uint32_t written = 0;
    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t id = 0; id + 2 < block_count; ++id)
    {
        CHECK(wirehair_encode(
            encoder, id, packet.data(), block_bytes, &written) ==
            Wirehair_Success);
        result = wirehair_decode(decoder, id, packet.data(), written);
        CHECK(result == Wirehair_NeedMore);
    }
    CHECK(wirehair_encode(
        encoder, collision_a, collision_packet_a.data(), block_bytes,
        &written) == Wirehair_Success);
    CHECK(wirehair_decode(
        decoder, collision_a, collision_packet_a.data(), written) ==
        Wirehair_NeedMore);
    CHECK(wirehair_encode(
        encoder, collision_b, collision_packet_b.data(), block_bytes,
        &written) == Wirehair_Success);
    CHECK(collision_packet_a == collision_packet_b);
    result = wirehair_decode(
        decoder, collision_b, collision_packet_b.data(), written);
    CHECK(result == Wirehair_NeedMore);
    CHECK(internal->TestingDecoderRowCount() == block_count);

    const uint32_t rows_before = internal->TestingDecoderRowCount();
    const uint32_t rank_before = internal->TestingDecoderRank();
    const uint32_t packets_before = internal->TestingReceivedPacketCount();
    for (unsigned repeat = 0; repeat < 256; ++repeat) {
        CHECK(wirehair_decode(
            decoder, collision_b, collision_packet_b.data(), block_bytes) ==
            Wirehair_NeedMore);
    }
    CHECK(internal->TestingDecoderRowCount() == rows_before);
    CHECK(internal->TestingDecoderRank() == rank_before);
    CHECK(internal->TestingReceivedPacketCount() == packets_before);
    std::vector<uint8_t> conflict_packet = collision_packet_b;
    conflict_packet[7] ^= 1;
    CHECK(wirehair_decode(
        decoder, collision_b, conflict_packet.data(), block_bytes) ==
        Wirehair_InvalidInput);

    for (uint32_t id = block_count - 2;
         id < block_count && result == Wirehair_NeedMore;
         ++id)
    {
        CHECK(wirehair_encode(
            encoder, id, packet.data(), block_bytes, &written) ==
            Wirehair_Success);
        result = wirehair_decode(decoder, id, packet.data(), written);
        CHECK(result == Wirehair_NeedMore || result == Wirehair_Success);
    }
    for (uint32_t id = 300000;
         result == Wirehair_NeedMore && id < 300256;
         ++id)
    {
        CHECK(wirehair_encode(
            encoder, id, packet.data(), block_bytes, &written) ==
            Wirehair_Success);
        result = wirehair_decode(decoder, id, packet.data(), written);
        CHECK(result == Wirehair_NeedMore || result == Wirehair_Success);
    }
    CHECK(result == Wirehair_Success);
    CHECK(internal->TestingDecoderRowCount() > rows_before);
    std::vector<uint8_t> recovered(message.size());
    CHECK(wirehair_recover(
        decoder, recovered.data(), recovered.size()) == Wirehair_Success);
    CHECK(recovered == message);

    wirehair_free(decoder);
    wirehair_free(encoder);
}

void TestDuplicatePacketActiveCapacity()
{
    const uint32_t block_count = 2;
    const uint32_t block_bytes = 1;
    const std::vector<uint8_t> message = MakeMessage(block_count, 79);
    WirehairCodec encoder = wirehair_encoder_create(
        nullptr, message.data(), message.size(), block_bytes);
    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message.size(), block_bytes);
    CHECK(encoder != nullptr && decoder != nullptr);
    if (!encoder || !decoder) {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return;
    }

    const wirehair::Codec* internal =
        reinterpret_cast<const wirehair::Codec*>(decoder);
    const uint32_t limit = block_count + 1024;
    CHECK(internal->TestingReceivedPacketLimit() == limit);
    CHECK(internal->TestingReceivedPacketSlots() == 4096);
    CHECK(internal->TestingReceivedPacketMemoryBytes() ==
        UINT64_C(4096) * 24);

    const uint16_t mix_count =
        static_cast<uint16_t>(wirehair::GetDenseCount(block_count) + 6);
    wirehair::PeelRowParameters target;
    target.Initialize(2, internal->PSeed(), block_count, mix_count);
    std::vector<uint32_t> equivalent_ids;
    equivalent_ids.push_back(2);
    for (uint32_t id = 3;
         id < 4000000 && equivalent_ids.size() < limit + 1;
         ++id)
    {
        wirehair::PeelRowParameters params;
        params.Initialize(
            id, internal->PSeed(), block_count, mix_count);
        if (params.PeelCount == target.PeelCount &&
            params.PeelFirst == target.PeelFirst &&
            params.PeelAdd == target.PeelAdd &&
            params.MixFirst == target.MixFirst &&
            params.MixAdd == target.MixAdd)
        {
            equivalent_ids.push_back(id);
        }
    }
    CHECK(equivalent_ids.size() == limit + 1);
    if (equivalent_ids.size() != limit + 1) {
        wirehair_free(decoder);
        wirehair_free(encoder);
        return;
    }

    uint8_t packet = 0;
    uint8_t first_packet = 0;
    uint32_t written = 0;
    for (uint32_t i = 0; i < limit; ++i)
    {
        CHECK(wirehair_encode(
            encoder, equivalent_ids[i], &packet, sizeof(packet), &written) ==
            Wirehair_Success);
        CHECK(written == 1);
        if (i == 0) {
            first_packet = packet;
        }
        else {
            CHECK(packet == first_packet);
        }
        CHECK(wirehair_decode(
            decoder, equivalent_ids[i], &packet, written) ==
            Wirehair_NeedMore);
    }
    CHECK(internal->TestingReceivedPacketCount() == limit);
    const uint32_t rows = internal->TestingDecoderRowCount();
    const uint32_t rank = internal->TestingDecoderRank();

    for (unsigned repeat = 0; repeat < 256; ++repeat) {
        CHECK(wirehair_decode(
            decoder, equivalent_ids[0], &first_packet, 1) ==
            Wirehair_NeedMore);
    }
    CHECK(internal->TestingDecoderRowCount() == rows);
    CHECK(internal->TestingDecoderRank() == rank);
    CHECK(internal->TestingReceivedPacketCount() == limit);

    uint8_t conflict = static_cast<uint8_t>(first_packet ^ 1);
    CHECK(wirehair_decode(
        decoder, equivalent_ids[0], &conflict, 1) ==
        Wirehair_InvalidInput);
    CHECK(wirehair_decode(
        decoder, equivalent_ids[limit], &first_packet, 1) ==
        Wirehair_ExtraInsufficient);
    CHECK(wirehair_decode(
        decoder, equivalent_ids[limit], &first_packet, 1) ==
        Wirehair_ExtraInsufficient);
    CHECK(wirehair_decode(
        decoder, equivalent_ids[0], &first_packet, 1) ==
        Wirehair_NeedMore);
    CHECK(internal->TestingDecoderRowCount() == rows);
    CHECK(internal->TestingDecoderRank() == rank);
    CHECK(internal->TestingReceivedPacketCount() == limit);

    wirehair_free(decoder);
    wirehair_free(encoder);
}

void TestEncoderDetachCore()
{
    CHECK(wirehair_encoder_detach_input(nullptr) == Wirehair_InvalidInput);

    const uint32_t block_bytes = 31;
    const uint32_t block_count = 65;
    const size_t message_bytes = block_bytes * block_count - 7;
    const std::vector<uint8_t> expected = MakeMessage(message_bytes, 83);

    typedef std::pair<uint32_t, std::vector<uint8_t> > Packet;
    std::vector<uint32_t> ids;
    for (uint32_t id = 0; id < block_count; ++id) {
        ids.push_back(id);
    }
    ids.push_back(block_count);
    ids.push_back(block_count + 1);
    ids.push_back(block_count + 37);
    ids.push_back(UINT32_MAX);

    const auto capture = [&](WirehairCodec encoder) {
        std::vector<Packet> packets;
        for (uint32_t id : ids)
        {
            std::vector<uint8_t> data(block_bytes, uint8_t{0xa5});
            uint32_t written = 0;
            CHECK(wirehair_encode(
                encoder, id, data.data(), data.size(), &written) ==
                Wirehair_Success);
            CHECK(written > 0 && written <= block_bytes);
            data.resize(written);
            packets.push_back(Packet(id, data));
        }
        return packets;
    };
    const auto compare = [&](WirehairCodec encoder,
                             const std::vector<Packet>& packets) {
        for (const Packet& packet : packets)
        {
            std::vector<uint8_t> actual(block_bytes, uint8_t{0x5a});
            uint32_t written = 0;
            CHECK(wirehair_encode(
                encoder, packet.first, actual.data(), actual.size(),
                &written) == Wirehair_Success);
            actual.resize(written);
            CHECK(actual == packet.second);
        }
    };

    std::vector<uint8_t> borrowed_source = expected;
    WirehairCodec borrowed = wirehair_encoder_create(
        nullptr, borrowed_source.data(), borrowed_source.size(), block_bytes);
    CHECK(borrowed != nullptr);
    if (!borrowed) {
        return;
    }
    wirehair::Codec* borrowed_internal =
        reinterpret_cast<wirehair::Codec*>(borrowed);
    CHECK(borrowed_internal->TestingRecoveryReady());
    CHECK(borrowed_internal->TestingInputAttached());
    CHECK(borrowed_internal->TestingInputAllocatedBytes() == 0);
    const std::vector<Packet> borrowed_packets = capture(borrowed);
    CHECK(wirehair_encoder_detach_input(borrowed) == Wirehair_Success);
    CHECK(wirehair_encoder_detach_input(borrowed) == Wirehair_Success);
    CHECK(!borrowed_internal->TestingInputAttached());
    CHECK(borrowed_internal->TestingInputAllocatedBytes() == 0);

    // Release the actual borrowed allocation before asking for any packet.
    // ASan therefore catches even a systematic fast-path use-after-free.
    std::fill(borrowed_source.begin(), borrowed_source.end(), uint8_t{0xcc});
    std::vector<uint8_t>().swap(borrowed_source);
    compare(borrowed, borrowed_packets);

    WirehairCodec owned = wirehair_encoder_create_owned(
        nullptr, expected.data(), expected.size(), block_bytes);
    CHECK(owned != nullptr);
    if (!owned) {
        wirehair_free(borrowed);
        return;
    }
    wirehair::Codec* owned_internal =
        reinterpret_cast<wirehair::Codec*>(owned);
    CHECK(owned_internal->TestingRecoveryReady());
    CHECK(owned_internal->TestingInputAttached());
    CHECK(owned_internal->TestingInputAllocatedBytes() ==
        static_cast<uint64_t>(block_count) * block_bytes);
    const std::vector<Packet> owned_packets = capture(owned);
    CHECK(wirehair_encoder_detach_input(owned) == Wirehair_Success);
    CHECK(!owned_internal->TestingInputAttached());
    CHECK(owned_internal->TestingInputAllocatedBytes() == 0);
    compare(owned, owned_packets);

    wirehair::Codec unmaterialized;
    CHECK(unmaterialized.InitializeEncoder(message_bytes, block_bytes) ==
        Wirehair_Success);
    CHECK(!unmaterialized.TestingRecoveryReady());
    CHECK(unmaterialized.DetachInput() == Wirehair_InvalidInput);

    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, expected.size(), block_bytes);
    CHECK(decoder != nullptr);
    if (decoder)
    {
        wirehair::Codec* decoder_internal =
            reinterpret_cast<wirehair::Codec*>(decoder);
        CHECK(wirehair_encoder_detach_input(decoder) == Wirehair_InvalidInput);
        WirehairResult result = Wirehair_NeedMore;
        std::vector<uint8_t> packet(block_bytes);
        for (uint32_t id = 0; id < block_count; ++id)
        {
            uint32_t written = 0;
            CHECK(wirehair_encode(
                owned, id, packet.data(), packet.size(), &written) ==
                Wirehair_Success);
            result = wirehair_decode(
                decoder, id, packet.data(), written);
        }
        CHECK(result == Wirehair_Success);
        CHECK(decoder_internal->TestingInputAllocatedBytes() ==
            static_cast<uint64_t>(block_count + CAT_MAX_EXTRA_ROWS) *
                block_bytes);
        CHECK(wirehair_encoder_detach_input(decoder) == Wirehair_InvalidInput);
        CHECK(wirehair_decoder_becomes_encoder(decoder) == Wirehair_Success);
        CHECK(decoder_internal->TestingRecoveryReady());
        const std::vector<Packet> converted_packets = capture(decoder);
        CHECK(wirehair_encoder_detach_input(decoder) == Wirehair_Success);
        CHECK(!decoder_internal->TestingInputAttached());
        CHECK(decoder_internal->TestingInputAllocatedBytes() == 0);
        compare(decoder, converted_packets);
    }

    std::vector<uint8_t> replacement = MakeMessage(message_bytes, 97);
    WirehairCodec reused = nullptr;
    CHECK(wirehair_encoder_create_owned_ex(
        owned, replacement.data(), replacement.size(), block_bytes,
        &reused) == Wirehair_Success);
    owned = nullptr;
    CHECK(reused != nullptr);
    if (reused)
    {
        wirehair::Codec* reused_internal =
            reinterpret_cast<wirehair::Codec*>(reused);
        CHECK(reused_internal->TestingInputAllocatedBytes() ==
            static_cast<uint64_t>(block_count) * block_bytes);
        const std::vector<Packet> replacement_packets = capture(reused);
        CHECK(wirehair_encoder_detach_input(reused) == Wirehair_Success);
        std::fill(replacement.begin(), replacement.end(), uint8_t{0x33});
        compare(reused, replacement_packets);
    }

    wirehair_free(reused);
    wirehair_free(decoder);
    wirehair_free(owned);
    wirehair_free(borrowed);
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
    const wirehair::Codec* internal =
        reinterpret_cast<const wirehair::Codec*>(decoder);
    CHECK(internal->TestingReceivedPacketLimit() == block_count + 1024);
    CHECK(internal->TestingReceivedPacketSlots() == 131072);
    CHECK(internal->TestingReceivedPacketMemoryBytes() ==
        UINT64_C(131072) * 24);
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

void BenchmarkDetachTradeoff()
{
    const uint32_t block_bytes = 1200;
    const uint32_t block_count = 1024;
    const std::vector<uint8_t> message =
        MakeMessage(static_cast<size_t>(block_bytes) * block_count, 101);
    WirehairCodec encoder = wirehair_encoder_create_owned(
        nullptr, message.data(), message.size(), block_bytes);
    CHECK(encoder != nullptr);
    if (!encoder) {
        return;
    }
    wirehair::Codec* internal = reinterpret_cast<wirehair::Codec*>(encoder);
    const uint64_t released_bytes = internal->TestingInputAllocatedBytes();
    std::vector<uint8_t> packet(block_bytes);

    const auto measure = [&](uint64_t& checksum) {
        checksum = 0;
        const auto start = std::chrono::steady_clock::now();
        for (unsigned repeat = 0; repeat < 2; ++repeat)
        {
            for (uint32_t id = 0; id < block_count; ++id)
            {
                uint32_t written = 0;
                CHECK(wirehair_encode(
                    encoder, id, packet.data(), packet.size(), &written) ==
                    Wirehair_Success);
                checksum += packet[0];
                checksum += packet[written - 1];
            }
        }
        return std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - start).count();
    };

    uint64_t attached_checksum = 0;
    uint64_t detached_checksum = 0;
    const long long attached_us = measure(attached_checksum);
    CHECK(wirehair_encoder_detach_input(encoder) == Wirehair_Success);
    const long long detached_us = measure(detached_checksum);
    CHECK(attached_checksum == detached_checksum);
    CHECK(released_bytes == message.size());
    CHECK(internal->TestingInputAllocatedBytes() == 0);
    std::printf(
        "Encoder detach tradeoff (2 x %u systematic, %u bytes): "
        "attached=%lld us detached=%lld us released=%llu bytes\n",
        block_count,
        block_bytes,
        attached_us,
        detached_us,
        static_cast<unsigned long long>(released_bytes));
    wirehair_free(encoder);
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
    if (argc == 2 && std::strcmp(argv[1], "--duplicates") == 0)
    {
        CHECK(wirehair_init() == Wirehair_Success);
        TestPacketFingerprintVectors();
        TestDuplicatePacketIdentity();
        TestDuplicatePacketDuringExtraResume();
        TestDuplicatePacketActiveCapacity();
        if (Failures != 0) {
            std::fprintf(stderr,
                "Legacy duplicate packet test failed: %d checks\n", Failures);
            return 1;
        }
        std::puts("Legacy duplicate packet test passed");
        return 0;
    }
    if (argc == 2 && std::strcmp(argv[1], "--detach") == 0)
    {
        CHECK(wirehair_init() == Wirehair_Success);
        TestEncoderDetachCore();
        BenchmarkDetachTradeoff();
        if (Failures != 0) {
            std::fprintf(stderr,
                "Legacy encoder detach test failed: %d checks\n", Failures);
            return 1;
        }
        std::puts("Legacy encoder detach test passed");
        return 0;
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

    TestGF256WireContract();
    TestDenseCountInterpolation();
    TestLegacyWireProfiles();
    TestCreationResultsAndBoundaries();
    TestOwnedInputAndRecoverBlock();
    TestLifecycleAndBorrowedMutation();
    TestPacketFingerprintVectors();
    TestDuplicatePacketIdentity();
    TestDuplicatePacketDuringExtraResume();
    TestDuplicatePacketActiveCapacity();
    TestEncoderDetachCore();
    TestHighNRecoveryIndex();
    BenchmarkOwnedCreationCost();
    BenchmarkDetachTradeoff();

    if (Failures != 0) {
        std::fprintf(stderr, "Legacy core test failed: %d checks\n", Failures);
        return 1;
    }
    std::puts("Legacy core test passed");
    return 0;
}
