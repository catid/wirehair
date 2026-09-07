#include "wirehair/wirehair_k6.h"
#include "wirehair/wirehair.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <vector>

namespace {
bool tracking = false;
size_t allocations = 0, fail_at = SIZE_MAX;
void* pointers[8] = {};
size_t sizes[8] = {};
}
void* operator new(size_t n)
{
    const size_t index = tracking ? allocations++ : SIZE_MAX;
    if (index == fail_at && tracking) throw std::bad_alloc();
    void* p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    if (index < 8) { pointers[index] = p; sizes[index] = n; }
    return p;
}
void* operator new[](size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void* operator new(size_t n, const std::nothrow_t&) noexcept
{
    try { return ::operator new(n); } catch (const std::bad_alloc&) { return nullptr; }
}
void* operator new[](size_t n, const std::nothrow_t&) noexcept
{
    try { return ::operator new[](n); } catch (const std::bad_alloc&) { return nullptr; }
}
void operator delete(void* p, const std::nothrow_t&) noexcept { std::free(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, size_t) noexcept { std::free(p); }
void operator delete[](void* p, size_t) noexcept { std::free(p); }
#endif

namespace {
using Byte = uint8_t;
using Matrix = std::array<Byte, 36>;
using Row = std::array<Byte, 6>;
using Profile = std::array<Byte, WIREHAIR_K6_PROFILE_BYTES>;
void Check(bool ok, const char* why)
{
    if (!ok) { std::cerr << "FAIL: " << why << '\n'; std::exit(1); }
}
void Start(size_t failure = SIZE_MAX)
{
    allocations = 0; fail_at = failure;
    std::fill(pointers, pointers + 8, nullptr); std::fill(sizes, sizes + 8, 0);
    tracking = true;
}
size_t Stop() { tracking = false; return allocations; }
Byte Multiply(Byte a, Byte b)
{
    unsigned p = 0;
    for (unsigned i = 0; i < 8; ++i) if (b & (1u << i)) p ^= unsigned(a) << i;
    for (int i = 14; i >= 8; --i) if (p & (1u << i)) p ^= 0x14du << (i - 8);
    return static_cast<Byte>(p);
}
Matrix Product(const Matrix& a, const Matrix& b)
{
    Matrix p = {};
    for (unsigned r = 0; r < 6; ++r) for (unsigned c = 0; c < 6; ++c)
        for (unsigned k = 0; k < 6; ++k) p[r * 6 + c] ^= Multiply(a[r * 6 + k], b[k * 6 + c]);
    return p;
}
struct Oracle {
    Matrix powers[2][32];
    Oracle()
    {
        const Byte feedback[2][6] = {{124, 127, 152, 84, 241, 63}, {125, 127, 152, 84, 241, 63}};
        for (unsigned phase = 0; phase < 2; ++phase) {
            powers[phase][0].fill(0);
            for (unsigned i = 0; i < 5; ++i) powers[phase][0][(i + 1) * 6 + i] = 1;
            for (unsigned i = 0; i < 6; ++i) powers[phase][0][i * 6 + 5] = feedback[phase][i];
        }
        for (unsigned level = 1; level < 32; ++level) {
            powers[0][level] = Product(powers[0][level - 1], powers[1][level - 1]);
            powers[1][level] = Product(powers[1][level - 1], powers[0][level - 1]);
        }
    }
    Row Coefficients(uint32_t id) const
    {
        Row result = {{1, 0, 0, 0, 0, 0}};
        for (unsigned bit = 0; bit < 32; ++bit) if (id & (uint32_t(1) << bit)) {
            unsigned phase = 0;
            for (unsigned higher = bit + 1; higher < 32; ++higher) phase ^= (id >> higher) & 1u;
            Row next = {};
            for (unsigned r = 0; r < 6; ++r) for (unsigned c = 0; c < 6; ++c)
                next[r] ^= Multiply(powers[phase][bit][r * 6 + c], result[c]);
            result = next;
        }
        return result;
    }
    std::vector<Byte> Packet(const std::vector<Byte>& source, uint32_t block, uint32_t id) const
    {
        std::vector<Byte> result(id == 5 ? source.size() - size_t(block) * 5 : block, 0);
        const Row row = Coefficients(id);
        for (size_t j = 0; j < result.size(); ++j) for (unsigned i = 0; i < 6; ++i) {
            const size_t offset = size_t(block) * i + j;
            if (offset < source.size()) result[j] ^= Multiply(row[i], source[offset]);
        }
        return result;
    }
};
void Put(Byte* p, uint64_t x, unsigned bytes)
{
    for (unsigned i = 0; i < bytes; ++i) { p[i] = static_cast<Byte>(x & 255); x >>= 8; }
}
Profile Descriptor(uint64_t message, uint32_t block)
{
    Profile p = {{'W', 'H', 'K', '6', 1, 0, 32, 0}};
    Put(p.data() + 8, WIREHAIR_K6_PROFILE_ID, 8);
    Put(p.data() + 16, message, 8); Put(p.data() + 24, block, 4);
    return p;
}
std::vector<Byte> Message(size_t bytes)
{
    std::vector<Byte> result(bytes);
    for (size_t i = 0; i < bytes; ++i) result[i] = static_cast<Byte>(37 * i + i / 11);
    return result;
}
const uint32_t ids[] = {0,1,2,3,4,5,6,7,8,9,10,11,31,255,1023,1024,1025,
    65535,131071,131072,0xffffff,0x1000000,0xffffffff,0xfffffffd,0xfffffffb,
    0xfffffff9,0xfffffff7,0xfffffff5};
void PacketCheck(WirehairK6Codec encoder, uint32_t id, const std::vector<Byte>& expected)
{
    std::vector<Byte> packet(expected.size() + 2, 0xa5);
    const WirehairK6Result short_result = wirehair_k6_encode(encoder, id, packet.data() + 1, expected.size() - 1);
    Check(short_result.status == WirehairK6_BufferTooSmall && !short_result.bytes_written &&
        short_result.bytes_required == expected.size() &&
        std::all_of(packet.begin(), packet.end(), [](Byte b) { return b == 0xa5; }), "short encode no-write");
    Start(0);
    const WirehairK6Result result = wirehair_k6_encode(encoder, id, packet.data() + 1, expected.size());
    const size_t count = Stop();
    Check(!count && result.status == WirehairK6_Success && result.bytes_required == expected.size() &&
        result.bytes_written == expected.size() && packet.front() == 0xa5 && packet.back() == 0xa5 &&
        std::equal(expected.begin(), expected.end(), packet.begin() + 1), "exact independent packet oracle, guards, no allocation");
}
void Lifecycle(const Oracle& oracle, uint32_t block, uint32_t tail, uint32_t policy)
{
    const std::vector<Byte> message = Message(size_t(block) * 5 + tail);
    std::vector<Byte> source = message;
    Profile profile;
    profile.fill(0xa5);
    const WirehairK6CreateResult encoder = wirehair_k6_encoder_create(source.data(), source.size(), block,
        policy, profile.data(), profile.size());
    Check(encoder.status == WirehairK6_Success && encoder.codec && profile == Descriptor(message.size(), block), "create and exact serialized bytes");
    Check(wirehair_v2_profile_validate(profile.data(), profile.size()) == WirehairV2_InvalidMagic,
          "existing WH2 API rejects distinct K6 wire format");
    if (policy == WirehairK6_Independent) {
        std::fill(source.begin(), source.end(), 0x5a);
        std::vector<Byte>().swap(source);
    } else {
        const auto alias = wirehair_k6_encode(encoder.codec, 6, source.data(), block);
        Check(alias.status == WirehairK6_InvalidInput && source == message, "borrowed-source write exclusion");
    }
    for (uint32_t id : ids) PacketCheck(encoder.codec, id, oracle.Packet(message, block, id));
    Check(wirehair_k6_encoder_detach_input(encoder.codec) == WirehairK6_Success, "detach succeeds");
    Start(0);
    const auto detached = wirehair_k6_encoder_detach_input(encoder.codec);
    Check(Stop() == 0 && detached == WirehairK6_Success, "detach idempotent allocation free");
    std::fill(source.begin(), source.end(), 0x3c);
    std::vector<Byte>().swap(source);
    for (uint32_t id : ids) PacketCheck(encoder.codec, id, oracle.Packet(message, block, id));
    wirehair_k6_free(encoder.codec);

    // Fresh decoder from bytes after encoder AND original source are gone.
    const auto decoder = wirehair_k6_decoder_create(profile.data(), profile.size());
    Check(decoder.status == WirehairK6_Success && decoder.codec, "standalone decoder");
    profile.fill(0); // Descriptor is not retained either.
    std::vector<Byte> recovered(message.size() + 2, 0xa5);
    auto r = wirehair_k6_recover(decoder.codec, recovered.data() + 1, message.size());
    Check(r.status == WirehairK6_NeedMore && !r.bytes_written && r.bytes_required == message.size(), "initial NeedMore");
    for (unsigned i = 0; i < 6; ++i) {
        const uint32_t id = UINT32_MAX - 2 * i;
        std::vector<Byte> packet = oracle.Packet(message, block, id);
        Start(0);
        const auto status = wirehair_k6_decode(decoder.codec, id, packet.data(), packet.size());
        const auto again = wirehair_k6_decode(decoder.codec, id, packet.data(), packet.size());
        Check(Stop() == 0 && (status == WirehairK6_NeedMore || status == WirehairK6_Success) && status == again,
              "distant feed, duplicates, no receive allocation");
        std::fill(packet.begin(), packet.end(), 0); // Feed owns accepted bytes.
    }
    for (unsigned i = 0; i < 6; ++i) {
        const auto packet = oracle.Packet(message, block, i);
        const auto status = wirehair_k6_decode(decoder.codec, i, packet.data(), packet.size());
        Check(status == WirehairK6_NeedMore || status == WirehairK6_Success, "resume with systematic equations");
    }
    r = wirehair_k6_recover(decoder.codec, recovered.data() + 1, message.size() - 1);
    Check(r.status == WirehairK6_BufferTooSmall && !r.bytes_written &&
        std::all_of(recovered.begin(), recovered.end(), [](Byte b) { return b == 0xa5; }), "failed recover no-write");
    for (unsigned repeat = 0; repeat < 2; ++repeat) {
        Start(0);
        r = wirehair_k6_recover(decoder.codec, recovered.data() + 1, message.size());
        Check(Stop() == 0 && r.status == WirehairK6_Success && r.bytes_written == message.size() &&
            r.bytes_required == message.size() && recovered.front() == 0xa5 && recovered.back() == 0xa5 &&
            std::equal(message.begin(), message.end(), recovered.begin() + 1), "repeat exact recovery, guards, no allocation");
        const auto packet = oracle.Packet(message, block, 6);
        Check(wirehair_k6_decode(decoder.codec, 6, packet.data(), packet.size()) == WirehairK6_Success, "feed after solved basis");
    }
    wirehair_k6_free(decoder.codec);
}

void Malformed()
{
    auto source = Message(384);
    const Profile good = Descriptor(source.size(), 64);
    auto reject = [&](const void* p, size_t bytes, WirehairK6Status expected) {
        Start(0);
        const auto s = wirehair_k6_profile_validate(p, bytes);
        const auto d = wirehair_k6_decoder_create(p, bytes);
        const auto e = wirehair_k6_encoder_create_profile(source.data(), p, bytes, WirehairK6_Independent);
        Check(Stop() == 0 && s == expected && d.status == expected && !d.codec &&
              e.status == expected && !e.codec, "malformed descriptor rejects without allocation");
    };
    for (size_t bytes = 0; bytes < 32; ++bytes) reject(good.data(), bytes, WirehairK6_InvalidInput);
    reject(good.data(), 33, WirehairK6_InvalidInput);
    reject(nullptr, 32, WirehairK6_InvalidInput);
    reject(reinterpret_cast<void*>(UINTPTR_MAX - 15), 32, WirehairK6_InvalidInput);
    for (unsigned i = 0; i < 32; ++i) if (i < 16 || i >= 28)
        for (unsigned bit = 0; bit < 8; ++bit) {
            Profile p = good; p[i] ^= Byte(1u << bit);
            reject(p.data(), p.size(), WirehairK6_UnsupportedProfile);
        }
    for (uint64_t id : {WIREHAIR_V2_PROFILE_CURRENT, UINT64_C(0xe161ce5d456f9bb7),
                        UINT64_C(0x20a4f27a870612a2), UINT64_C(0), UINT64_MAX}) {
        Profile p = good; Put(p.data() + 8, id, 8);
        reject(p.data(), p.size(), WirehairK6_UnsupportedProfile);
        std::memcpy(p.data(), "WHV2", 4);
        reject(p.data(), p.size(), WirehairK6_UnsupportedProfile);
        const auto status = wirehair_v2_profile_validate(p.data(), p.size());
        Check(status == (id == WIREHAIR_V2_PROFILE_CURRENT ? WirehairV2_Success : WirehairV2_UnsupportedProfile),
              "actual public current/retired profile validation unchanged");
    }
    for (uint64_t message : {UINT64_C(0), UINT64_C(320), UINT64_C(385), UINT64_MAX}) {
        const Profile p = Descriptor(message, 64);
        reject(p.data(), p.size(), WirehairK6_InvalidDimensions);
    }
    for (uint32_t block : {0u, WIREHAIR_K6_MAX_BLOCK_BYTES + 1, UINT32_MAX}) {
        const Profile p = Descriptor(uint64_t(block) * 6, block);
        reject(p.data(), p.size(), WirehairK6_InvalidDimensions);
    }
    Profile p;
    p.fill(0xa5);
    for (uint32_t policy : {0u, 3u, UINT32_MAX}) {
        Start(0);
        auto e = wirehair_k6_encoder_create(source.data(), source.size(), 64, policy, p.data(), p.size());
        Check(Stop() == 0 && e.status == WirehairK6_InvalidInput && !e.codec && p[0] == 0xa5, "invalid policy");
    }
    Start(0);
    const auto short_profile = wirehair_k6_encoder_create(source.data(), source.size(), 64,
        WirehairK6_Independent, p.data(), p.size() - 1);
    Check(Stop() == 0 && short_profile.status == WirehairK6_BufferTooSmall && !short_profile.codec &&
          std::all_of(p.begin(), p.end(), [](Byte b) { return b == 0xa5; }), "short profile no work/write");
    const auto original = source;
    Start(0);
    auto e = wirehair_k6_encoder_create(source.data(), source.size(), 64, WirehairK6_BorrowedImmutable, source.data(), 32);
    Check(Stop() == 0 && e.status == WirehairK6_InvalidInput && !e.codec && source == original, "borrowed descriptor alias rejected");
    e = wirehair_k6_encoder_create(source.data(), source.size(), 64, WirehairK6_Independent, source.data(), 32);
    Check(e.status == WirehairK6_Success && std::equal(good.begin(), good.end(), source.begin()), "independent descriptor alias staged");
    Check(wirehair_k6_decode(e.codec, 0, source.data(), 64) == WirehairK6_InvalidInput &&
          wirehair_k6_recover(e.codec, source.data(), source.size()).status == WirehairK6_InvalidInput,
          "encoder cannot decode or recover");
    Check(wirehair_k6_encode(e.codec, 0, source.data(), SIZE_MAX).status == WirehairK6_InvalidInput,
          "encoder rejects capacity wrapping");
    PacketCheck(e.codec, 0, std::vector<Byte>(original.begin(), original.begin() + 64));
    wirehair_k6_free(e.codec);
}

// Exercise every packed table selector and both parity states through the
// public encoder. Six byte-basis source blocks expose the full row directly;
// expected rows come from the polynomial oracle, never the packed table.
void LookupCoverage(const Oracle& oracle)
{
    std::vector<Byte> source(36, 0);
    for (unsigned i = 0; i < 6; ++i) source[i * 6 + i] = 1;
    Profile profile;
    const auto encoder = wirehair_k6_encoder_create(source.data(), source.size(), 6,
        WirehairK6_BorrowedImmutable, profile.data(), profile.size());
    Check(encoder.status == WirehairK6_Success, "lookup coverage encoder");
    const auto check = [&](uint32_t id) {
        const Row expected = oracle.Coefficients(id);
        Byte output[6] = {};
        const auto r = wirehair_k6_encode(encoder.codec, id, output, sizeof(output));
        Check(r.status == WirehairK6_Success && r.bytes_written == 6 &&
            std::equal(expected.begin(), expected.end(), output), "all lookup selectors match polynomial oracle");
    };
    for (uint32_t low = 0; low < 1024; ++low) {
        check(low); check(1024u | low);
    }
    for (uint32_t mid = 0; mid < 128; ++mid) for (uint32_t column = 0; column < 6; ++column) {
        check((mid << 10) | column);
        check((1u << 17) | (mid << 10) | column);
        check((mid << 17) | column);
        check((1u << 24) | (mid << 17) | column);
    }
    for (uint32_t high = 0; high < 256; ++high) for (uint32_t column = 0; column < 6; ++column)
        check((high << 24) | column);
    wirehair_k6_free(encoder.codec);
}

void FailuresAndAliases(const Oracle& oracle)
{
    for (unsigned tail : {1u, 64u}) for (uint32_t policy : {1u, 2u}) {
        auto source = Message(320 + tail);
        Profile profile = Descriptor(source.size(), 64);
        const auto expected = oracle.Packet(source, 64, UINT32_MAX);
        Start();
        auto created = wirehair_k6_encoder_create(source.data(), source.size(), 64, policy, profile.data(), profile.size());
        const size_t count = Stop();
        Check(created.status == WirehairK6_Success && count == 2 + (policy == 1 ? 1u : 0u) + (tail != 64 ? 1u : 0u), "account every encoder allocation");
        for (size_t i = 0; i < count; ++i) {
            std::vector<Byte> snapshot(static_cast<Byte*>(pointers[i]), static_cast<Byte*>(pointers[i]) + sizes[i]);
            const auto alias = wirehair_k6_encode(created.codec, 0, pointers[i], 64);
            Check(alias.status == WirehairK6_InvalidInput && std::equal(snapshot.begin(), snapshot.end(), static_cast<Byte*>(pointers[i])), "all private encoder allocations protected");
        }
        wirehair_k6_free(created.codec);
        for (size_t failure = 0; failure < count; ++failure) {
            profile.fill(0xa5); Start(failure);
            created = wirehair_k6_encoder_create(source.data(), source.size(), 64, policy, profile.data(), profile.size());
            const size_t attempted = Stop();
            Check(attempted == failure + 1 && created.status == WirehairK6_OutOfMemory && !created.codec &&
                std::all_of(profile.begin(), profile.end(), [](Byte b) { return b == 0xa5; }), "every encoder OOM atomic");
        }
        profile = Descriptor(source.size(), 64);
        for (size_t failure = 0; failure < count; ++failure) {
            Start(failure);
            created = wirehair_k6_encoder_create_profile(source.data(), profile.data(), profile.size(), policy);
            Check(Stop() == failure + 1 && created.status == WirehairK6_OutOfMemory && !created.codec &&
                  profile == Descriptor(source.size(), 64), "every serialized encoder OOM atomic");
        }
        created = wirehair_k6_encoder_create_profile(source.data(), profile.data(), profile.size(), policy);
        Check(created.status == WirehairK6_Success, "serialized encoder reconstruction");
        PacketCheck(created.codec, UINT32_MAX, expected);
        wirehair_k6_free(created.codec);
        if (policy == 2) for (size_t failure = 0; failure < (tail == 64 ? 2u : 3u); ++failure) {
            created = wirehair_k6_encoder_create_profile(source.data(), profile.data(), profile.size(), policy);
            Check(created.status == WirehairK6_Success, "detach fixture");
            Start(failure);
            const auto status = wirehair_k6_encoder_detach_input(created.codec);
            Check(Stop() == failure + 1 && status == WirehairK6_OutOfMemory, "every detach OOM");
            PacketCheck(created.codec, UINT32_MAX, expected);
            Check(wirehair_k6_encode(created.codec, 0, source.data(), 64).status == WirehairK6_InvalidInput, "failed detach retains borrowing");
            Check(wirehair_k6_encoder_detach_input(created.codec) == WirehairK6_Success, "detach retry succeeds");
            wirehair_k6_free(created.codec);
        }
    }
    const Profile profile = Descriptor(383, 64);
    for (size_t failure = 0; failure < 3; ++failure) {
        Start(failure);
        const auto d = wirehair_k6_decoder_create(profile.data(), profile.size());
        Check(Stop() == failure + 1 && d.status == WirehairK6_OutOfMemory && !d.codec, "every decoder OOM");
    }
    Start();
    const auto d = wirehair_k6_decoder_create(profile.data(), profile.size());
    Check(Stop() == 3 && d.status == WirehairK6_Success, "decoder wrapper/core/slab accounting");
    for (unsigned i = 0; i < 3; ++i) {
        std::vector<Byte> snapshot(static_cast<Byte*>(pointers[i]), static_cast<Byte*>(pointers[i]) + sizes[i]);
        Check(wirehair_k6_decode(d.codec, 0, pointers[i], 64) == WirehairK6_InvalidInput &&
            wirehair_k6_recover(d.codec, pointers[i], 383).status == WirehairK6_InvalidInput &&
            std::equal(snapshot.begin(), snapshot.end(), static_cast<Byte*>(pointers[i])), "all private decoder allocations protected");
    }
    Check(wirehair_k6_encoder_detach_input(d.codec) == WirehairK6_InvalidInput, "decoder cannot detach");
    Byte packet[64] = {};
    Check(wirehair_k6_encode(d.codec, 0, packet, 64).status == WirehairK6_InvalidInput, "decoder cannot encode");
    Check(wirehair_k6_decode(d.codec, 0, packet, 63) == WirehairK6_InvalidInput, "invalid feed does not poison");
    Check(wirehair_k6_decode(d.codec, 0, packet, 64) == WirehairK6_NeedMore, "rank one");
    packet[0] = 1;
    Check(wirehair_k6_decode(d.codec, 0, packet, 64) == WirehairK6_Conflict, "contradiction detected");
    packet[0] = 0;
    Check(wirehair_k6_decode(d.codec, 0, packet, 64) == WirehairK6_Conflict, "persistent poison");
    std::array<Byte, 383> output;
    output.fill(0xa5);
    const auto r = wirehair_k6_recover(d.codec, output.data(), output.size());
    Check(r.status == WirehairK6_Conflict && !r.bytes_written && r.bytes_required == output.size() &&
        std::all_of(output.begin(), output.end(), [](Byte b) { return b == 0xa5; }), "poison recover no-write");
    wirehair_k6_free(d.codec);
    const auto resumed = wirehair_k6_decoder_create(profile.data(), profile.size());
    Check(resumed.status == WirehairK6_Success, "resume fixture");
    for (unsigned id = 0; id < 5; ++id)
        Check(wirehair_k6_decode(resumed.codec, id, packet, sizeof(packet)) == WirehairK6_NeedMore, "retain deficient basis");
    Check(wirehair_k6_recover(resumed.codec, output.data(), output.size()).status == WirehairK6_NeedMore &&
        std::all_of(output.begin(), output.end(), [](Byte b) { return b == 0xa5; }), "NeedMore preserves output and basis");
    Check(wirehair_k6_decode(resumed.codec, 5, packet, 64) == WirehairK6_InvalidInput &&
          wirehair_k6_decode(resumed.codec, 5, packet, 63) == WirehairK6_Success, "exact partial tail and resume");
    const void* wrapped = reinterpret_cast<const void*>(UINTPTR_MAX - 15);
    Check(wirehair_k6_decode(resumed.codec, 0, wrapped, 64) == WirehairK6_InvalidInput &&
          wirehair_k6_recover(resumed.codec, output.data(), SIZE_MAX).status == WirehairK6_InvalidInput,
          "wrapped input/capacity rejected");
    Check(wirehair_k6_recover(resumed.codec, output.data(), output.size()).status == WirehairK6_Success &&
          std::all_of(output.begin(), output.end(), [](Byte b) { return b == 0; }), "resumed zero-message oracle");
    packet[0] = 1;
    Check(wirehair_k6_decode(resumed.codec, 0, packet, 64) == WirehairK6_Conflict, "contradiction after solved basis");
    wirehair_k6_free(resumed.codec);
    wirehair_k6_free(nullptr);
    Check(wirehair_k6_encoder_detach_input(nullptr) == WirehairK6_InvalidInput &&
        wirehair_k6_encode(nullptr, 0, packet, 64).status == WirehairK6_InvalidInput &&
        wirehair_k6_decode(nullptr, 0, packet, 64) == WirehairK6_InvalidInput &&
        wirehair_k6_recover(nullptr, packet, 64).status == WirehairK6_InvalidInput, "null handles");
}
} // namespace

int main()
{
    Check(wirehair_init() == Wirehair_Success, "shared runtime initialization");
    const Oracle oracle;
    unsigned shapes = 0;
    for (uint32_t block : {1u,2u,3u,7u,16u,31u,32u,63u,64u,65u,127u,128u,129u,1280u}) {
        std::vector<uint32_t> tails = {1, block};
        if (block > 1) tails.push_back(block - 1);
        std::sort(tails.begin(), tails.end()); tails.erase(std::unique(tails.begin(), tails.end()), tails.end());
        for (uint32_t tail : tails) for (uint32_t policy : {1u,2u}) {
            Lifecycle(oracle, block, tail, policy); ++shapes;
        }
    }
    Malformed(); LookupCoverage(oracle); FailuresAndAliases(oracle);
    std::cout << "PASS " << shapes << " serialized lifecycle shapes; independent packet/recovery oracle, C ABI separately, ownership/detach/OOM/alias/profile/poison gates\n";
}
