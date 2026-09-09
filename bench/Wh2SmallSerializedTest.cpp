#include "Wh2SmallSerializedCore.h"
#include "../codec/WirehairK6Payload.h"
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


#ifndef WH2_SMALL_TEST_K
#define WH2_SMALL_TEST_K WH2_SMALL_CODEC_K
#endif
#if WH2_SMALL_TEST_K == 2 || WH2_SMALL_TEST_K == 3 || WH2_SMALL_TEST_K == 5
static_assert(WH2_SMALL_TEST_K == WH2_SMALL_CODEC_K, "test matches external boundary");
struct Api {
    static constexpr uint64_t ProfileId = WH2_SMALL_PROFILE_ID;
    static constexpr uint32_t MaxBlockBytes = WH2_SMALL_MAX_BLOCK_BYTES;
    static Wh2SmallStatus ProfileValidate(const void* p,size_t n) { return wh2_small_profile_validate(p,n); }
    static Wh2SmallCreateResult EncoderCreate(const void* s,uint64_t m,uint32_t b,uint32_t policy,void* p,size_t n)
    { return wh2_small_encoder_create(s,m,b,policy,p,n); }
    static Wh2SmallCreateResult EncoderCreateProfile(const void* s,const void* p,size_t n,uint32_t policy)
    { return wh2_small_encoder_create_profile(s,p,n,policy); }
    static Wh2SmallCreateResult DecoderCreate(const void* p,size_t n) { return wh2_small_decoder_create(p,n); }
    static Wh2SmallStatus Detach(Wh2SmallCodec h) { return wh2_small_encoder_detach_input(h); }
    static Wh2SmallResult Encode(Wh2SmallCodec h,uint32_t id,void* p,size_t n) { return wh2_small_encode(h,id,p,n); }
    static Wh2SmallStatus Decode(Wh2SmallCodec h,uint32_t id,const void* p,size_t n) { return wh2_small_decode(h,id,p,n); }
    static Wh2SmallResult Recover(Wh2SmallCodec h,void* p,size_t n) { return wh2_small_recover(h,p,n); }
    static void Free(Wh2SmallCodec h) { wh2_small_free(h); }
};
#elif WH2_SMALL_TEST_K == 6
namespace {
alignas(64) const uint8_t six_lookup[] = {
#include "../codec/WirehairK6Lookup.inc"
};
struct SixTraits {
    static wirehair_small_core::Lookup View() { return {six_lookup,sizeof(six_lookup)}; }
};
}
using Api = wh2_small_serialized::Facade<6,SixTraits>;
static_assert(Api::ProfileId == WIREHAIR_K6_PROFILE_ID, "old K6 profile identity");
static_assert(Api::MaxBlockBytes == WIREHAIR_K6_MAX_BLOCK_BYTES, "old K6 memory policy");
#else
#error "Unsupported test dimension"
#endif

namespace {
bool tracking = false;
size_t allocations = 0, fail_at = SIZE_MAX;
void* pointers[8] = {};
size_t sizes[8] = {};
}
__attribute__((noinline)) void* operator new(size_t n)
{
    const size_t index = tracking ? allocations++ : SIZE_MAX;
    if (index == fail_at && tracking) throw std::bad_alloc();
    void* p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    if (index < 8) { pointers[index] = p; sizes[index] = n; }
    return p;
}
__attribute__((noinline)) void* operator new[](size_t n) { return ::operator new(n); }
__attribute__((noinline)) void operator delete(void* p) noexcept { std::free(p); }
__attribute__((noinline)) void operator delete[](void* p) noexcept { std::free(p); }
__attribute__((noinline)) void* operator new(size_t n, const std::nothrow_t&) noexcept
{
    try { return ::operator new(n); } catch (const std::bad_alloc&) { return nullptr; }
}
__attribute__((noinline)) void* operator new[](size_t n, const std::nothrow_t&) noexcept
{
    try { return ::operator new[](n); } catch (const std::bad_alloc&) { return nullptr; }
}
__attribute__((noinline)) void operator delete(void* p, const std::nothrow_t&) noexcept { std::free(p); }
__attribute__((noinline)) void operator delete[](void* p, const std::nothrow_t&) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
__attribute__((noinline)) void operator delete(void* p, size_t) noexcept { std::free(p); }
__attribute__((noinline)) void operator delete[](void* p, size_t) noexcept { std::free(p); }
#endif

namespace {
using Byte = uint8_t;
constexpr unsigned K = WH2_SMALL_TEST_K;
using Matrix = std::array<Byte, K * K>;
using Row = std::array<Byte, K>;
using Profile = std::array<Byte, WH2_SMALL_PROFILE_BYTES>;
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
    for (unsigned r = 0; r < K; ++r) for (unsigned c = 0; c < K; ++c)
        for (unsigned k = 0; k < K; ++k) p[r * K + c] ^= Multiply(a[r * K + k], b[k * K + c]);
    return p;
}
struct Oracle {
    Matrix powers[2][32];
    Oracle()
    {
        const Byte six[6] = {124, 127, 152, 84, 241, 63};
        const Byte three[3] = {8, 14, 7};
        const Byte five[5] = {121, 110, 207, 198, 31};
        const Byte two[2] = {2, 3};
        for (unsigned phase = 0; phase < 2; ++phase) {
            powers[phase][0].fill(0);
            for (unsigned i = 0; i < K - 1; ++i) powers[phase][0][(i + 1) * K + i] = 1;
            for (unsigned i = 0; i < K; ++i) powers[phase][0][i * K + K - 1] =
                static_cast<Byte>((K == 2 ? two[i] : K == 3 ? three[i] : K == 5 ? five[i] : six[i]) ^
                    (i == 0 ? phase : 0));
        }
        for (unsigned level = 1; level < 32; ++level) {
            powers[0][level] = Product(powers[0][level - 1], powers[1][level - 1]);
            powers[1][level] = Product(powers[1][level - 1], powers[0][level - 1]);
        }
    }
    Row Coefficients(uint32_t id) const
    {
        Row result = {}; result[0] = 1;
        for (unsigned bit = 0; bit < 32; ++bit) if (id & (uint32_t(1) << bit)) {
            unsigned phase = 0;
            for (unsigned higher = bit + 1; higher < 32; ++higher) phase ^= (id >> higher) & 1u;
            Row next = {};
            for (unsigned r = 0; r < K; ++r) for (unsigned c = 0; c < K; ++c)
                next[r] ^= Multiply(powers[phase][bit][r * K + c], result[c]);
            result = next;
        }
        return result;
    }
    std::vector<Byte> Packet(const std::vector<Byte>& source, uint32_t block, uint32_t id) const
    {
        std::vector<Byte> result(id == K - 1 ? source.size() - size_t(block) * (K - 1) : block, 0);
        const Row row = Coefficients(id);
        for (size_t j = 0; j < result.size(); ++j) for (unsigned i = 0; i < K; ++i) {
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
    Profile p = {{'W', 'H', 'K', Byte('0' + K), 1, 0, 32, 0}};
    Put(p.data() + 8, Api::ProfileId, 8);
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
void PacketCheck(Wh2SmallCodec encoder, uint32_t id, const std::vector<Byte>& expected,
                 WirehairK6Codec legacy = nullptr)
{
    std::vector<Byte> packet(expected.size() + 2, 0xa5);
    const Wh2SmallResult short_result = Api::Encode(encoder, id, packet.data() + 1, expected.size() - 1);
    Check(short_result.status == Wh2Small_BufferTooSmall && !short_result.bytes_written &&
        short_result.bytes_required == expected.size() &&
        std::all_of(packet.begin(), packet.end(), [](Byte b) { return b == 0xa5; }), "short encode no-write");
    Start(0);
    const Wh2SmallResult result = Api::Encode(encoder, id, packet.data() + 1, expected.size());
    const size_t count = Stop();
    Check(!count && result.status == Wh2Small_Success && result.bytes_required == expected.size() &&
        result.bytes_written == expected.size() && packet.front() == 0xa5 && packet.back() == 0xa5 &&
        std::equal(expected.begin(), expected.end(), packet.begin() + 1), "exact independent packet oracle, guards, no allocation");
    if (legacy) {
        std::vector<Byte> old(packet.size(), 0xa5);
        Start(0);
        const auto r = wirehair_k6_encode(legacy, id, old.data() + 1, expected.size());
        Check(Stop() == 0 && int(r.status) == int(result.status) && r.bytes_required == result.bytes_required &&
              r.bytes_written == result.bytes_written && old == packet, "installed K6 packet/result parity");
    }
}
void Lifecycle(const Oracle& oracle, uint32_t block, uint32_t tail, uint32_t policy)
{
    const std::vector<Byte> message = Message(size_t(block) * (K - 1) + tail);
    std::vector<Byte> source = message;
    Profile profile;
    profile.fill(0xa5);
    const Wh2SmallCreateResult encoder = Api::EncoderCreate(source.data(), source.size(), block,
        policy, profile.data(), profile.size());
    Check(encoder.status == Wh2Small_Success && encoder.codec && profile == Descriptor(message.size(), block), "create and exact serialized bytes");
    Check(wirehair_v2_profile_validate(profile.data(), profile.size()) == WirehairV2_InvalidMagic,
          "production rejects experimental wire format");
    WirehairK6Codec legacy_encoder = nullptr;
    if (K == 6) {
        Profile old_profile = {};
        const auto old = wirehair_k6_encoder_create(source.data(), source.size(), block, policy,
            old_profile.data(), old_profile.size());
        Check(old.status == WirehairK6_Success && old.codec && old_profile == profile,
              "installed K6 descriptor/create parity");
        legacy_encoder = old.codec;
    } else {
        Check(wirehair_k6_profile_validate(profile.data(), profile.size()) == WirehairK6_UnsupportedProfile,
              "installed K6 rejects other dimension profile");
    }
    if (policy == Wh2Small_Independent) {
        std::fill(source.begin(), source.end(), 0x5a);
        std::vector<Byte>().swap(source);
    } else {
        const auto alias = Api::Encode(encoder.codec, K, source.data(), block);
        Check(alias.status == Wh2Small_InvalidInput && source == message, "borrowed-source write exclusion");
    }
    for (uint32_t id : ids) PacketCheck(encoder.codec, id, oracle.Packet(message, block, id), legacy_encoder);
    Check(Api::Detach(encoder.codec) == Wh2Small_Success, "detach succeeds");
    if (legacy_encoder) Check(wirehair_k6_encoder_detach_input(legacy_encoder) == WirehairK6_Success,
                              "installed K6 detach parity");
    Start(0);
    const auto detached = Api::Detach(encoder.codec);
    Check(Stop() == 0 && detached == Wh2Small_Success, "detach idempotent allocation free");
    std::fill(source.begin(), source.end(), 0x3c);
    std::vector<Byte>().swap(source);
    for (uint32_t id : ids) PacketCheck(encoder.codec, id, oracle.Packet(message, block, id), legacy_encoder);
    Api::Free(encoder.codec);
    wirehair_k6_free(legacy_encoder);

    // Fresh decoder from bytes after encoder AND original source are gone.
    const auto decoder = Api::DecoderCreate(profile.data(), profile.size());
    Check(decoder.status == Wh2Small_Success && decoder.codec, "standalone decoder");
    WirehairK6Codec legacy_decoder = nullptr;
    if (K == 6) {
        const auto old = wirehair_k6_decoder_create(profile.data(), profile.size());
        Check(old.status == WirehairK6_Success && old.codec, "installed K6 decoder create parity");
        legacy_decoder = old.codec;
    }
    profile.fill(0); // Descriptor is not retained either.
    std::vector<Byte> recovered(message.size() + 2, 0xa5);
    auto r = Api::Recover(decoder.codec, recovered.data() + 1, message.size());
    Check(r.status == Wh2Small_NeedMore && !r.bytes_written && r.bytes_required == message.size(), "initial NeedMore");
    for (unsigned i = 0; i < K; ++i) {
        const uint32_t id = UINT32_MAX - 2 * i;
        std::vector<Byte> packet = oracle.Packet(message, block, id);
        Start(0);
        const auto status = Api::Decode(decoder.codec, id, packet.data(), packet.size());
        const auto again = Api::Decode(decoder.codec, id, packet.data(), packet.size());
        Check(Stop() == 0 && (status == Wh2Small_NeedMore || status == Wh2Small_Success) && status == again,
              "distant feed, duplicates, no receive allocation");
        if (legacy_decoder) {
            Start(0);
            const auto old = wirehair_k6_decode(legacy_decoder, id, packet.data(), packet.size());
            const auto repeat = wirehair_k6_decode(legacy_decoder, id, packet.data(), packet.size());
            Check(Stop() == 0 && int(old) == int(status) && int(repeat) == int(again), "installed K6 feed parity");
        }
        std::fill(packet.begin(), packet.end(), 0); // Feed owns accepted bytes.
    }
    for (unsigned i = 0; i < K; ++i) {
        const auto packet = oracle.Packet(message, block, i);
        const auto status = Api::Decode(decoder.codec, i, packet.data(), packet.size());
        Check(status == Wh2Small_NeedMore || status == Wh2Small_Success, "resume with systematic equations");
        if (legacy_decoder) Check(int(wirehair_k6_decode(legacy_decoder, i, packet.data(), packet.size())) == int(status),
                                  "installed K6 resume parity");
    }
    r = Api::Recover(decoder.codec, recovered.data() + 1, message.size() - 1);
    Check(r.status == Wh2Small_BufferTooSmall && !r.bytes_written &&
        std::all_of(recovered.begin(), recovered.end(), [](Byte b) { return b == 0xa5; }), "failed recover no-write");
    for (unsigned repeat = 0; repeat < 2; ++repeat) {
        std::fill(recovered.begin(), recovered.end(), 0xa5);
        Start(0);
        r = Api::Recover(decoder.codec, recovered.data() + 1, message.size());
        Check(Stop() == 0 && r.status == Wh2Small_Success && r.bytes_written == message.size() &&
            r.bytes_required == message.size() && recovered.front() == 0xa5 && recovered.back() == 0xa5 &&
            std::equal(message.begin(), message.end(), recovered.begin() + 1), "repeat exact recovery, guards, no allocation");
        if (legacy_decoder) {
            std::vector<Byte> old_output(recovered.size(), 0xa5);
            Start(0);
            const auto old = wirehair_k6_recover(legacy_decoder, old_output.data() + 1, message.size());
            Check(Stop() == 0 && int(old.status) == int(r.status) && old.bytes_written == r.bytes_written &&
                old.bytes_required == r.bytes_required && old_output == recovered, "installed K6 recovery parity");
        }
        const auto packet = oracle.Packet(message, block, K);
        Check(Api::Decode(decoder.codec, K, packet.data(), packet.size()) == Wh2Small_Success, "feed after solved basis");
        if (legacy_decoder) Check(wirehair_k6_decode(legacy_decoder, K, packet.data(), packet.size()) == WirehairK6_Success,
                                  "installed K6 solved-feed parity");
    }
    Api::Free(decoder.codec);
    wirehair_k6_free(legacy_decoder);
}

void Malformed()
{
    auto source = Message(K * 64);
    const Profile good = Descriptor(source.size(), 64);
    auto reject = [&](const void* p, size_t bytes, Wh2SmallStatus expected) {
        Start(0);
        const auto s = Api::ProfileValidate(p, bytes);
        const auto d = Api::DecoderCreate(p, bytes);
        const auto e = Api::EncoderCreateProfile(source.data(), p, bytes, Wh2Small_Independent);
        Check(Stop() == 0 && s == expected && d.status == expected && !d.codec &&
              e.status == expected && !e.codec, "malformed descriptor rejects without allocation");
    };
    for (size_t bytes = 0; bytes < 32; ++bytes) reject(good.data(), bytes, Wh2Small_InvalidInput);
    reject(good.data(), 33, Wh2Small_InvalidInput);
    reject(nullptr, 32, Wh2Small_InvalidInput);
    reject(reinterpret_cast<void*>(UINTPTR_MAX - 15), 32, Wh2Small_InvalidInput);
    for (unsigned other_k : {2u, 3u, 5u, 6u}) if (other_k != K) {
        Profile other = good;
        other[3] = static_cast<Byte>('0' + other_k);
        Put(other.data() + 8, UINT64_C(0x5748324b30544d31) + (uint64_t(other_k) << 24), 8);
        reject(other.data(), other.size(), Wh2Small_UnsupportedProfile);
    }
    for (unsigned i = 0; i < 32; ++i) if (i < 16 || i >= 28)
        for (unsigned bit = 0; bit < 8; ++bit) {
            Profile p = good; p[i] ^= Byte(1u << bit);
            reject(p.data(), p.size(), Wh2Small_UnsupportedProfile);
        }
    for (uint64_t id : {WIREHAIR_V2_PROFILE_CURRENT, UINT64_C(0xe161ce5d456f9bb7),
                        UINT64_C(0x20a4f27a870612a2), UINT64_C(0), UINT64_MAX}) {
        Profile p = good; Put(p.data() + 8, id, 8);
        reject(p.data(), p.size(), Wh2Small_UnsupportedProfile);
        std::memcpy(p.data(), "WHV2", 4);
        reject(p.data(), p.size(), Wh2Small_UnsupportedProfile);
        const auto status = wirehair_v2_profile_validate(p.data(), p.size());
        Check(status == (id == WIREHAIR_V2_PROFILE_CURRENT ? WirehairV2_Success : WirehairV2_UnsupportedProfile),
              "actual public current/retired profile validation unchanged");
    }
    for (uint64_t message : {UINT64_C(0), uint64_t((K - 1) * 64), uint64_t(K * 64 + 1), UINT64_MAX}) {
        const Profile p = Descriptor(message, 64);
        reject(p.data(), p.size(), Wh2Small_InvalidDimensions);
    }
    for (uint32_t block : {0u, Api::MaxBlockBytes + 1, UINT32_MAX}) {
        const Profile p = Descriptor(uint64_t(block) * K, block);
        reject(p.data(), p.size(), Wh2Small_InvalidDimensions);
    }
    {
        const Profile maximum = Descriptor(uint64_t(Api::MaxBlockBytes) * K, Api::MaxBlockBytes);
        Start(0);
        const auto s = Api::ProfileValidate(maximum.data(), maximum.size());
        Check(Stop() == 0 && s == Wh2Small_Success, "maximum slab shape accepted without allocation");
    }
    Profile p;
    p.fill(0xa5);
    for (uint32_t policy : {0u, 3u, UINT32_MAX}) {
        Start(0);
        auto e = Api::EncoderCreate(source.data(), source.size(), 64, policy, p.data(), p.size());
        Check(Stop() == 0 && e.status == Wh2Small_InvalidInput && !e.codec && p[0] == 0xa5, "invalid policy");
    }
    Start(0);
    const auto short_profile = Api::EncoderCreate(source.data(), source.size(), 64,
        Wh2Small_Independent, p.data(), p.size() - 1);
    Check(Stop() == 0 && short_profile.status == Wh2Small_BufferTooSmall && !short_profile.codec &&
          std::all_of(p.begin(), p.end(), [](Byte b) { return b == 0xa5; }), "short profile no work/write");
    const auto original = source;
    for (const void* invalid : {static_cast<const void*>(nullptr), reinterpret_cast<const void*>(UINTPTR_MAX - 15)}) {
        Start(0);
        const auto bad_source = Api::EncoderCreate(invalid, source.size(), 64, Wh2Small_Independent, p.data(), p.size());
        const auto bad_output = Api::EncoderCreate(source.data(), source.size(), 64, Wh2Small_Independent,
            const_cast<void*>(invalid), p.size());
        Check(Stop() == 0 && bad_source.status == Wh2Small_InvalidInput && !bad_source.codec &&
            bad_output.status == Wh2Small_InvalidInput && !bad_output.codec, "null/wrapped create buffers rejected");
    }
    Start(0);
    auto e = Api::EncoderCreate(source.data(), source.size(), 64, Wh2Small_BorrowedImmutable, source.data(), 32);
    Check(Stop() == 0 && e.status == Wh2Small_InvalidInput && !e.codec && source == original, "borrowed descriptor alias rejected");
    e = Api::EncoderCreate(source.data(), source.size(), 64, Wh2Small_Independent, source.data(), 32);
    Check(e.status == Wh2Small_Success && std::equal(good.begin(), good.end(), source.begin()), "independent descriptor alias staged");
    Check(Api::Decode(e.codec, 0, source.data(), 64) == Wh2Small_InvalidInput &&
          Api::Recover(e.codec, source.data(), source.size()).status == Wh2Small_InvalidInput,
          "encoder cannot decode or recover");
    Check(Api::Encode(e.codec, 0, source.data(), SIZE_MAX).status == Wh2Small_InvalidInput,
          "encoder rejects capacity wrapping");
    PacketCheck(e.codec, 0, std::vector<Byte>(original.begin(), original.begin() + 64));
    Api::Free(e.codec);
}

void FailuresAndAliases(const Oracle& oracle)
{
    for (unsigned tail : {1u, 64u}) for (uint32_t policy : {1u, 2u}) {
        auto source = Message((K - 1) * 64 + tail);
        Profile profile = Descriptor(source.size(), 64);
        const auto expected = oracle.Packet(source, 64, UINT32_MAX);
        Start();
        auto created = Api::EncoderCreate(source.data(), source.size(), 64, policy, profile.data(), profile.size());
        const size_t count = Stop();
        Check(created.status == Wh2Small_Success && count == 2 + (policy == 1 ? 1u : 0u) + (tail != 64 ? 1u : 0u), "account every encoder allocation");
        for (size_t i = 0; i < count; ++i) {
            std::vector<Byte> snapshot(static_cast<Byte*>(pointers[i]), static_cast<Byte*>(pointers[i]) + sizes[i]);
            const auto alias = Api::Encode(created.codec, 0, pointers[i], 64);
            Check(alias.status == Wh2Small_InvalidInput && std::equal(snapshot.begin(), snapshot.end(), static_cast<Byte*>(pointers[i])), "all private encoder allocations protected");
        }
        Api::Free(created.codec);
        for (size_t failure = 0; failure < count; ++failure) {
            profile.fill(0xa5); Start(failure);
            created = Api::EncoderCreate(source.data(), source.size(), 64, policy, profile.data(), profile.size());
            const size_t attempted = Stop();
            Check(attempted == failure + 1 && created.status == Wh2Small_OutOfMemory && !created.codec &&
                std::all_of(profile.begin(), profile.end(), [](Byte b) { return b == 0xa5; }), "every encoder OOM atomic");
        }
        profile = Descriptor(source.size(), 64);
        for (size_t failure = 0; failure < count; ++failure) {
            Start(failure);
            created = Api::EncoderCreateProfile(source.data(), profile.data(), profile.size(), policy);
            Check(Stop() == failure + 1 && created.status == Wh2Small_OutOfMemory && !created.codec &&
                  profile == Descriptor(source.size(), 64), "every serialized encoder OOM atomic");
        }
        created = Api::EncoderCreateProfile(source.data(), profile.data(), profile.size(), policy);
        Check(created.status == Wh2Small_Success, "serialized encoder reconstruction");
        PacketCheck(created.codec, UINT32_MAX, expected);
        Api::Free(created.codec);
        if (policy == 2) for (size_t failure = 0; failure < (tail == 64 ? 2u : 3u); ++failure) {
            created = Api::EncoderCreateProfile(source.data(), profile.data(), profile.size(), policy);
            Check(created.status == Wh2Small_Success, "detach fixture");
            Start(failure);
            const auto status = Api::Detach(created.codec);
            Check(Stop() == failure + 1 && status == Wh2Small_OutOfMemory, "every detach OOM");
            PacketCheck(created.codec, UINT32_MAX, expected);
            Check(Api::Encode(created.codec, 0, source.data(), 64).status == Wh2Small_InvalidInput, "failed detach retains borrowing");
            Check(Api::Detach(created.codec) == Wh2Small_Success, "detach retry succeeds");
            Api::Free(created.codec);
        }
    }
    const Profile profile = Descriptor((K * 64 - 1), 64);
    for (size_t failure = 0; failure < 3; ++failure) {
        Start(failure);
        const auto d = Api::DecoderCreate(profile.data(), profile.size());
        Check(Stop() == failure + 1 && d.status == Wh2Small_OutOfMemory && !d.codec, "every decoder OOM");
    }
    Start();
    const auto d = Api::DecoderCreate(profile.data(), profile.size());
    Check(Stop() == 3 && d.status == Wh2Small_Success, "decoder wrapper/core/slab accounting");
    for (unsigned i = 0; i < 3; ++i) {
        std::vector<Byte> snapshot(static_cast<Byte*>(pointers[i]), static_cast<Byte*>(pointers[i]) + sizes[i]);
        Check(Api::Decode(d.codec, 0, pointers[i], 64) == Wh2Small_InvalidInput &&
            Api::Recover(d.codec, pointers[i], (K * 64 - 1)).status == Wh2Small_InvalidInput &&
            std::equal(snapshot.begin(), snapshot.end(), static_cast<Byte*>(pointers[i])), "all private decoder allocations protected");
    }
    Check(Api::Detach(d.codec) == Wh2Small_InvalidInput, "decoder cannot detach");
    Byte packet[64] = {};
    Check(Api::Encode(d.codec, 0, packet, 64).status == Wh2Small_InvalidInput, "decoder cannot encode");
    Check(Api::Decode(d.codec, 0, packet, 63) == Wh2Small_InvalidInput, "invalid feed does not poison");
    Check(Api::Decode(d.codec, 0, packet, 64) == Wh2Small_NeedMore, "rank one");
    packet[0] = 1;
    Check(Api::Decode(d.codec, 0, packet, 64) == Wh2Small_Conflict, "contradiction detected");
    packet[0] = 0;
    Check(Api::Decode(d.codec, 0, packet, 64) == Wh2Small_Conflict, "persistent poison");
    std::array<Byte, (K * 64 - 1)> output;
    output.fill(0xa5);
    const auto r = Api::Recover(d.codec, output.data(), output.size());
    Check(r.status == Wh2Small_Conflict && !r.bytes_written && r.bytes_required == output.size() &&
        std::all_of(output.begin(), output.end(), [](Byte b) { return b == 0xa5; }), "poison recover no-write");
    Api::Free(d.codec);
    const auto resumed = Api::DecoderCreate(profile.data(), profile.size());
    Check(resumed.status == Wh2Small_Success, "resume fixture");
    for (unsigned id = 0; id < K - 1; ++id)
        Check(Api::Decode(resumed.codec, id, packet, sizeof(packet)) == Wh2Small_NeedMore, "retain deficient basis");
    Check(Api::Recover(resumed.codec, output.data(), output.size()).status == Wh2Small_NeedMore &&
        std::all_of(output.begin(), output.end(), [](Byte b) { return b == 0xa5; }), "NeedMore preserves output and basis");
    Check(Api::Decode(resumed.codec, K - 1, packet, 64) == Wh2Small_InvalidInput &&
          Api::Decode(resumed.codec, K - 1, packet, 63) == Wh2Small_Success, "exact partial tail and resume");
    const void* wrapped = reinterpret_cast<const void*>(UINTPTR_MAX - 15);
    Check(Api::Decode(resumed.codec, 0, wrapped, 64) == Wh2Small_InvalidInput &&
          Api::Recover(resumed.codec, output.data(), SIZE_MAX).status == Wh2Small_InvalidInput,
          "wrapped input/capacity rejected");
    Check(Api::Recover(resumed.codec, output.data(), output.size()).status == Wh2Small_Success &&
          std::all_of(output.begin(), output.end(), [](Byte b) { return b == 0; }), "resumed zero-message oracle");
    packet[0] = 1;
    Check(Api::Decode(resumed.codec, 0, packet, 64) == Wh2Small_Conflict, "contradiction after solved basis");
    Api::Free(resumed.codec);
    Api::Free(nullptr);
    Check(Api::Detach(nullptr) == Wh2Small_InvalidInput &&
        Api::Encode(nullptr, 0, packet, 64).status == Wh2Small_InvalidInput &&
        Api::Decode(nullptr, 0, packet, 64) == Wh2Small_InvalidInput &&
        Api::Recover(nullptr, packet, 64).status == Wh2Small_InvalidInput, "null handles");
}
} // namespace

int main()
{
    Check(wirehair_init() == Wirehair_Success, "shared runtime initialization");
#ifdef WH2_SMALL_EXPECT_PORTABLE
    Check(!wirehair_k6_payload::Available(), "portable GF backend selected");
#endif
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
    Check(shapes == 78, "exact lifecycle shape roster");
    Malformed(); FailuresAndAliases(oracle);
    std::cout << "PASS K" << K << " " << shapes << " serialized lifecycle shapes; independent packet/recovery oracle, C ABI separately, ownership/detach/OOM/alias/profile/poison gates; GFNI="
              << wirehair_k6_payload::Available() << '\n';
}
