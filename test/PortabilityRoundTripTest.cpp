#include <wirehair/wirehair.h>

#include "../gf256.h"

#include <cstdint>
#include <cstdio>
#include <cstring>

namespace {

const char* ArchitectureName()
{
#if defined(__aarch64__) || defined(_M_ARM64)
    return "aarch64";
#elif defined(__arm__) || defined(_M_ARM)
    return "arm";
#elif defined(__s390x__)
    return "s390x";
#elif defined(__powerpc64__)
    return "ppc64";
#elif defined(__riscv) && __riscv_xlen == 64
    return "riscv64";
#elif defined(__i386__) || defined(_M_IX86)
    return "i686";
#elif defined(__x86_64__) || defined(_M_X64) || defined(_M_AMD64)
    return "x86_64";
#else
    return "unknown";
#endif
}

const char* BackendName()
{
#if defined(GF256_TRY_NEON)
    return "neon";
#elif defined(GF256_TARGET_MOBILE)
    return "scalar";
#else
    return "x86-simd";
#endif
}

const char* CompileEndianName()
{
#if defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) && \
    __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return "big";
#else
    return "little";
#endif
}

bool RuntimeIsBigEndian()
{
    const uint16_t value = UINT16_C(0x0102);
    return *reinterpret_cast<const uint8_t*>(&value) == 1;
}

bool Check(bool condition, const char* description)
{
    if (condition) {
        return true;
    }
    std::fprintf(stderr, "portability round trip failed: %s\n", description);
    return false;
}

} // namespace

int main()
{
    const bool compile_big = std::strcmp(CompileEndianName(), "big") == 0;
    if (!Check(
            compile_big == RuntimeIsBigEndian(),
            "compile-time and runtime byte order disagree"))
    {
        return 1;
    }

    std::printf(
        "Portability diagnostics: arch=%s pointer-bits=%u endian=%s "
        "gf256-backend=%s\n",
        ArchitectureName(), static_cast<unsigned>(sizeof(void*) * 8),
        CompileEndianName(), BackendName());

    if (!Check(wirehair_init() == Wirehair_Success, "wirehair_init")) {
        return 1;
    }

    enum : unsigned {
        MessageBytes = 128,
        BlockBytes = 16,
        BlockCount = MessageBytes / BlockBytes
    };
    uint8_t message[MessageBytes];
    for (unsigned i = 0; i < MessageBytes; ++i) {
        message[i] = static_cast<uint8_t>(i * 73u + 11u);
    }

    WirehairCodec encoder = wirehair_encoder_create(
        nullptr, message, MessageBytes, BlockBytes);
    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, MessageBytes, BlockBytes);
    if (!Check(encoder != nullptr && decoder != nullptr, "codec creation")) {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return 1;
    }

    // This externally visible packet pins the GF(256) heavy-field and packet
    // byte-order contract on every architecture in the portability matrix.
    static const uint8_t ExpectedRepair[BlockBytes] = {
        0x17, 0xed, 0xd1, 0x3d, 0xc5, 0xa4, 0x7a, 0x85,
        0xc4, 0x43, 0xc6, 0x57, 0x64, 0xe8, 0x7e, 0x56
    };
    uint8_t block[BlockBytes] = {};
    uint32_t written = 0;
    bool ok = Check(
        wirehair_encode(
            encoder, UINT32_C(12345), block, sizeof(block), &written) ==
            Wirehair_Success,
        "golden repair encode");
    ok = Check(written == sizeof(block), "golden repair length") && ok;
    ok = Check(
        std::memcmp(block, ExpectedRepair, sizeof(block)) == 0,
        "golden repair bytes") && ok;
    if (!ok) {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return 1;
    }

    WirehairResult decode_result = Wirehair_NeedMore;
    for (uint32_t id = 0; id < BlockCount; ++id)
    {
        if (id == 2) {
            continue;
        }
        written = 0;
        ok = Check(
            wirehair_encode(encoder, id, block, sizeof(block), &written) ==
                Wirehair_Success,
            "systematic encode") && ok;
        if (!ok) {
            break;
        }
        decode_result = wirehair_decode(decoder, id, block, written);
        ok = Check(
            decode_result == Wirehair_NeedMore,
            "decoder completed before the missing block was supplied") && ok;
    }

    for (uint32_t id = BlockCount;
        ok && decode_result == Wirehair_NeedMore && id < BlockCount + 64;
        ++id)
    {
        written = 0;
        ok = Check(
            wirehair_encode(encoder, id, block, sizeof(block), &written) ==
                Wirehair_Success,
            "repair encode") && ok;
        if (!ok) {
            break;
        }
        decode_result = wirehair_decode(decoder, id, block, written);
        ok = Check(
            decode_result == Wirehair_NeedMore ||
                decode_result == Wirehair_Success,
            "repair decode") && ok;
    }

    ok = Check(decode_result == Wirehair_Success, "decoder completion") && ok;
    if (decode_result == Wirehair_Success)
    {
        uint8_t recovered[MessageBytes] = {};
        ok = Check(
            wirehair_recover(decoder, recovered, sizeof(recovered)) ==
                Wirehair_Success,
            "message recovery") && ok;
        ok = Check(
            std::memcmp(message, recovered, sizeof(message)) == 0,
            "recovered message bytes") && ok;
    }

    wirehair_free(encoder);
    wirehair_free(decoder);
    if (!ok) {
        return 1;
    }

    // The serialized V2 profile and a recovery packet are separately pinned
    // here so every emulated endian/architecture lane executes the public
    // descriptor parser rather than merely compiling it.
    static const uint8_t ExpectedV2Profile[
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {
        0x57, 0x48, 0x56, 0x32, 0x01, 0x00, 0x20, 0x00,
        0xc9, 0xf9, 0xf4, 0x47, 0xbb, 0x5b, 0x29, 0x4b,
        0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };
    static const uint8_t ExpectedV2Repair[BlockBytes] = {
        0xe0, 0x0c, 0x23, 0x1b, 0x7c, 0x47, 0xff, 0xd4,
        0x80, 0x17, 0x47, 0x48, 0x26, 0xcb, 0x88, 0xb4
    };
    uint8_t v2_profile[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t v2_profile_bytes = 0;
    WirehairV2Codec v2_encoder = nullptr;
    WirehairV2Codec v2_decoder = nullptr;
    ok = Check(wirehair_v2_encoder_create(
        message, MessageBytes, BlockBytes,
        v2_profile, sizeof(v2_profile), &v2_profile_bytes, &v2_encoder) ==
        WirehairV2_Success, "V2 encoder selection") && ok;
    ok = Check(v2_profile_bytes == sizeof(v2_profile) &&
        std::memcmp(v2_profile, ExpectedV2Profile, sizeof(v2_profile)) == 0,
        "V2 cross-endian profile golden") && ok;
    ok = Check(wirehair_v2_decoder_create(
        v2_profile, v2_profile_bytes, &v2_decoder) == WirehairV2_Success,
        "V2 decoder from profile") && ok;
    written = 0;
    ok = Check(wirehair_v2_encode(
        v2_encoder, UINT32_C(12345), block, sizeof(block), &written) ==
        WirehairV2_Success && written == sizeof(block),
        "V2 repair encode") && ok;
    if (ok && std::memcmp(block, ExpectedV2Repair, sizeof(block)) != 0)
    {
        std::fprintf(stderr, "V2 portability repair golden actual:");
        for (uint8_t byte : block) {
            std::fprintf(stderr, " 0x%02x", byte);
        }
        std::fprintf(stderr, "\n");
        ok = false;
    }

    WirehairV2Result v2_decode_result = WirehairV2_NeedMore;
    for (uint32_t id = 0; ok && id < BlockCount; ++id)
    {
        written = 0;
        ok = Check(wirehair_v2_encode(
            v2_encoder, id, block, sizeof(block), &written) ==
            WirehairV2_Success, "V2 systematic encode") && ok;
        if (!ok) {
            break;
        }
        v2_decode_result = wirehair_v2_decode(
            v2_decoder, id, block, written);
        ok = Check(v2_decode_result == (id + 1u == BlockCount ?
            WirehairV2_Success : WirehairV2_NeedMore),
            "V2 systematic decode") && ok;
    }
    uint8_t v2_recovered[MessageBytes] = {};
    uint64_t v2_recovered_bytes = 0;
    ok = Check(v2_decode_result == WirehairV2_Success &&
        wirehair_v2_recover(
            v2_decoder, v2_recovered, sizeof(v2_recovered),
            &v2_recovered_bytes) == WirehairV2_Success,
        "V2 profile recovery") && ok;
    ok = Check(v2_recovered_bytes == MessageBytes &&
        std::memcmp(v2_recovered, message, sizeof(message)) == 0,
        "V2 recovered bytes") && ok;
    wirehair_v2_free(v2_encoder);
    wirehair_v2_free(v2_decoder);
    if (!ok) {
        return 1;
    }
    std::puts("Portability golden round trip passed");
    return 0;
}
