// Metadata-only contract probe: no encoder, GF initialization, huge allocation,
// production validator helper or production serializer builds the wire inputs.
#include <wirehair/wirehair.h>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
using Byte = std::uint8_t;
using Status = WirehairV2Result;
using Host = WirehairV2Profile;
using Record = std::array<Byte, 32>;
std::size_t host_cases = 0, wire_cases = 0;
struct Kind { std::uint64_t id; unsigned k; };
const Kind kinds[] = {
    {UINT64_C(0x4b295bbb47f4f9c9), 0},
    {UINT64_C(0x67c1043ecaa9e184), 3},
    {UINT64_C(0x80070c81bfe375f1), 5},
    {UINT64_C(0x7a9276b85c730ae0), 8}
};
static_assert(sizeof(Host) == 32, "Fixed descriptor ABI");
static_assert(WIREHAIR_V2_PROFILE_CURRENT == UINT64_C(0x4b295bbb47f4f9c9), "CURRENT is not a small-profile dispatch policy");
void Check(bool value, const char* what)
{
    if (!value) {
        std::fprintf(stderr, "FAIL %s at host=%zu wire=%zu\n", what, host_cases, wire_cases);
        std::exit(1);
    }
}
Host Profile(std::uint64_t id, std::uint64_t message, std::uint32_t block, unsigned attempt)
{
    Host p = {};
    p.struct_bytes = 32; p.profile_version = 1; p.profile_id = id;
    p.message_bytes = message; p.block_bytes = block;
    p.seed_attempt = static_cast<Byte>(attempt);
    return p;
}
Status Expected(const Host& p)
{
    if (p.struct_bytes != 32) return WirehairV2_InvalidSize;
    if (p.profile_version != 1) return WirehairV2_UnsupportedVersion;
    for (Byte b : p.reserved) if (b) return WirehairV2_ReservedNonzero;
    const Kind* selected = nullptr;
    for (const Kind& kind : kinds) if (p.profile_id == kind.id) selected = &kind;
    if (!selected) return WirehairV2_UnsupportedProfile;
    // Independent floor formula; no M+B overflow and no production ceil helper.
    if (!p.message_bytes || !p.block_bytes || p.block_bytes > UINT32_C(0x7fffffff))
        return WirehairV2_InvalidDimensions;
    const std::uint64_t last = (p.message_bytes - 1) / p.block_bytes;
    if (last < 1 || last >= 64000) return WirehairV2_InvalidDimensions;
    if (selected->k) {
        if (last != selected->k - 1 ||
            p.block_bytes > UINT32_C(268435456) / (selected->k + 1))
            return WirehairV2_InvalidDimensions;
        if (p.seed_attempt) return WirehairV2_BadSeed;
    }
    return WirehairV2_Success;
}
void Little(Byte* output, std::uint64_t value, unsigned bytes)
{
    for (unsigned i = 0; i < bytes; ++i) output[i] = static_cast<Byte>(value >> (8 * i));
}
Record Pack(const Host& p)
{
    Record r = {{'W', 'H', 'V', '2', 1, 0, 32, 0}};
    Little(r.data() + 8, p.profile_id, 8);
    Little(r.data() + 16, p.message_bytes, 8);
    Little(r.data() + 24, p.block_bytes, 4);
    r[28] = p.seed_attempt;
    for (unsigned i = 0; i < 3; ++i) r[29 + i] = p.reserved[i];
    return r;
}
void Serialize(const Host& p, Status expected)
{
    ++host_cases;
    const Host before = p;
    const Record packed = Pack(p);
    for (std::uint32_t capacity : {0u, 31u, 32u, 40u}) {
        std::array<Byte, 64> output; output.fill(0xa7);
        auto wanted = output;
        const Status result = expected == WirehairV2_Success && capacity < 32 ?
            WirehairV2_BufferTooSmall : expected;
        if (result == WirehairV2_Success) std::memcpy(wanted.data() + 8, packed.data(), 32);
        std::uint32_t bytes = UINT32_MAX;
        Check(wirehair_v2_profile_serialize(&p, output.data() + 8, capacity, &bytes) == result,
            "serialize status/precedence");
        Check(bytes == 32 && output == wanted, "serialize bytes/guards/transaction");
    }
    std::uint32_t bytes = UINT32_MAX;
    Check(wirehair_v2_profile_serialize(&p, nullptr, 32, &bytes) ==
        (expected == WirehairV2_Success ? WirehairV2_BufferTooSmall : expected), "null output precedence");
    Check(bytes == 32 && std::memcmp(&p, &before, 32) == 0, "host input unchanged");
}
void Parse(const Record& record, Status expected, const Host* wanted = nullptr)
{
    ++wire_cases;
    std::array<Byte, 48> input; input.fill(0x5e);
    std::memcpy(input.data() + 8, record.data(), 32);
    const auto before = input;
    struct Guarded { Byte before[8]; Host value; Byte after[8]; } output;
    std::memset(&output, 0xb6, sizeof(output));
    Guarded reference = output;
    if (expected == WirehairV2_Success) {
        Check(wanted != nullptr, "successful parse has independent expected host");
        reference.value = *wanted;
    } else std::memset(&reference.value, 0, sizeof(Host));
    Check(wirehair_v2_profile_deserialize(input.data() + 8, 32, &output.value) == expected,
        "parse status/precedence");
    Check(std::memcmp(&output, &reference, sizeof(output)) == 0, "parse complete output/guards");
    Check(wirehair_v2_profile_validate(input.data() + 8, 32) == expected, "validate status");
    Check(input == before, "wire input and guards unchanged");
}
void Case(const Host& p)
{
    const Status expected = Expected(p);
    Serialize(p, expected);
    Parse(Pack(p), expected, &p);
}
void Matrix()
{
    struct Shape { std::uint64_t message; std::uint32_t block; };
    std::array<Shape, 128> shapes = {};
    std::size_t count = 0;
    const auto add = [&](std::uint64_t message, std::uint32_t block) {
        Check(count < shapes.size(), "geometry capacity");
        shapes[count++] = Shape{message, block};
    };
    for (unsigned k : {2u, 3u, 4u, 5u, 6u, 8u, 9u}) for (unsigned b : {1u, 2u}) {
        add(std::uint64_t(k) * b, b); add(std::uint64_t(k - 1) * b + 1, b);
    }
    add(0, 0); add(0, 2); add(2, 0); add(1, 1); add(2, 2);
    add(64000, 1); add(64001, 1); add(128000, 2); add(128001, 2);
    add(UINT64_MAX, 1); add(UINT64_MAX, UINT32_MAX);
    add(UINT64_C(0x7fffffff) * 2, UINT32_C(0x7fffffff));
    add(UINT64_C(0x80000000) * 2, UINT32_C(0x80000000));
    for (unsigned k : {3u, 5u, 8u}) {
        const std::uint32_t bound = UINT32_C(268435456) / (k + 1);
        for (std::uint32_t b : {bound, bound + 1}) {
            add(std::uint64_t(k - 1) * b, b); add(std::uint64_t(k - 1) * b + 1, b);
            add(std::uint64_t(k) * b, b); add(std::uint64_t(k) * b + 1, b);
        }
    }
    Check(count == 65, "frozen geometry count");
    for (const Kind& kind : kinds) for (unsigned attempt = 0; attempt < 256; ++attempt)
        for (std::size_t i = 0; i < count; ++i)
            Case(Profile(kind.id, shapes[i].message, shapes[i].block, attempt));
}
void Precedence()
{
    // Unknown IDs, permanent tombstones and all separate prototype/API IDs.
    const std::uint64_t unsupported[] = {0, UINT64_MAX,
        UINT64_C(0xe161ce5d456f9bb7), UINT64_C(0x20a4f27a870612a2),
        UINT64_C(0x5748324b32544d31), UINT64_C(0x5748324b33544d31),
        UINT64_C(0x5748324b35544d31), UINT64_C(0x5748324b36544d31),
        UINT64_C(0x5748324b38544d31)};
    for (std::uint64_t id : unsupported) for (unsigned attempt : {0u, 255u}) {
        for (bool valid_dimensions : {false, true}) {
            Host p = Profile(id, valid_dimensions ? 6 : 0, valid_dimensions ? 2 : 0, attempt);
            Check(Expected(p) == WirehairV2_UnsupportedProfile, "unsupported oracle");
            Case(p);
            for (unsigned i = 0; i < 3; ++i) {
                Host bad = p; bad.reserved[i] = 1;
                Serialize(bad, WirehairV2_ReservedNonzero);
                Parse(Pack(bad), WirehairV2_ReservedNonzero);
                bad.profile_version = 2;
                Serialize(bad, WirehairV2_UnsupportedVersion);
                bad.struct_bytes = 31;
                Serialize(bad, WirehairV2_InvalidSize);
            }
            Record wire = Pack(p);
            wire[29] = 1;
            wire[6] = 31;
            Parse(wire, WirehairV2_InvalidSize);
            wire[4] = 2;
            Parse(wire, WirehairV2_UnsupportedVersion);
            wire[0] = '?';
            Parse(wire, WirehairV2_InvalidMagic);
        }
    }
    for (const Kind& kind : kinds) {
        Host p = Profile(kind.id, 0, 0, 255);
        Check(Expected(p) == WirehairV2_InvalidDimensions, "dimensions before seed oracle");
        Case(p);
        p.reserved[2] = 255;
        Serialize(p, WirehairV2_ReservedNonzero);
        Parse(Pack(p), WirehairV2_ReservedNonzero);
    }
}
} // namespace
int main()
{
    Matrix(); Precedence();
    Check(host_cases == 66928 && wire_cases == 66820, "complete frozen case counts");
    std::printf("Profile validation passed: %zu host cases, %zu wire cases; no timing or recovery sample\n",
        host_cases, wire_cases);
}
