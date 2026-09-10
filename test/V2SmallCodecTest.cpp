#include <wirehair/wirehair.h>
#include <wirehair/wirehair.hpp>
#include "../codec/WirehairV2Codec.h"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <vector>
#if defined(__unix__)
#include <sys/mman.h>
#include <unistd.h>
#endif

#if defined(_MSC_VER)
#define SMALL_TEST_NOINLINE __declspec(noinline)
#else
#define SMALL_TEST_NOINLINE __attribute__((noinline))
#endif

namespace {
bool tracking = false;
size_t allocations = 0, first_allocation_bytes = 0, fail_at = SIZE_MAX;
}
SMALL_TEST_NOINLINE void* operator new(size_t n)
{
    const size_t index = tracking ? allocations++ : SIZE_MAX;
    if (tracking && index == 0) first_allocation_bytes = n;
    if (index == fail_at && tracking) throw std::bad_alloc();
    void* p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    return p;
}
SMALL_TEST_NOINLINE void* operator new[](size_t n) { return ::operator new(n); }
SMALL_TEST_NOINLINE void operator delete(void* p) noexcept { std::free(p); }
SMALL_TEST_NOINLINE void operator delete[](void* p) noexcept { std::free(p); }
SMALL_TEST_NOINLINE void* operator new(size_t n, const std::nothrow_t&) noexcept
{
    try { return ::operator new(n); } catch (const std::bad_alloc&) { return nullptr; }
}
SMALL_TEST_NOINLINE void* operator new[](size_t n, const std::nothrow_t&) noexcept
{
    try { return ::operator new[](n); } catch (const std::bad_alloc&) { return nullptr; }
}
SMALL_TEST_NOINLINE void operator delete(void* p, const std::nothrow_t&) noexcept { std::free(p); }
SMALL_TEST_NOINLINE void operator delete[](void* p, const std::nothrow_t&) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
SMALL_TEST_NOINLINE void operator delete(void* p, size_t) noexcept { std::free(p); }
SMALL_TEST_NOINLINE void operator delete[](void* p, size_t) noexcept { std::free(p); }
#endif


namespace {
using Byte = uint8_t;
#ifndef WIREHAIR_V2_SMALL_TEST_K
#define WIREHAIR_V2_SMALL_TEST_K 3
#endif
constexpr unsigned K = WIREHAIR_V2_SMALL_TEST_K;
static_assert(K == 3 || K == 5 || K == 8, "Installed small WHV2 profiles only");
constexpr uint64_t ProfileId = K == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
    K == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09;
constexpr uint32_t MaxBlockBytes = UINT32_C(268435456) / (K + 1);
// K5/K8 are explicit-only until ordinary-path admission gates pass.
constexpr unsigned FirstRoute = K == 3 ? 0 : 1;
using Matrix = std::array<Byte, K * K>;
using Row = std::array<Byte, K>;
using Profile = std::array<Byte, 32>;
void Check(bool value, const char* what)
{
    if (!value) { std::cerr << "FAIL: " << what << '\n'; std::exit(1); }
}
void Start(size_t failure = SIZE_MAX)
{ allocations = 0; first_allocation_bytes = 0; fail_at = failure; tracking = true; }
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
        const Byte three[3] = {8, 14, 7};
        const Byte five[5] = {121, 110, 207, 198, 31};
        const Byte eight[8] = {96, 19, 186, 153, 85, 252, 7, 255};
        const Byte* feedback = K == 3 ? three : K == 5 ? five : eight;
        const unsigned lambda = K == 8 ? 2 : 1;
        for (unsigned phase = 0; phase < 2; ++phase) {
            powers[phase][0].fill(0);
            for (unsigned i = 0; i < K - 1; ++i) powers[phase][0][(i + 1) * K + i] = 1;
            for (unsigned i = 0; i < K; ++i) powers[phase][0][i * K + K - 1] = static_cast<Byte>(feedback[i] ^ (i == 0 ? lambda * phase : 0));
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

Profile Descriptor(uint64_t message, uint32_t block)
{
    WirehairV2Profile host = {};
    host.struct_bytes = sizeof(host);
    host.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    host.profile_id = ProfileId;
    host.message_bytes = message;
    host.block_bytes = block;
    Profile p = {};
    Check(wirehair_v2_profile_serialize(&host, p.data(), static_cast<uint32_t>(p.size()), nullptr) == WirehairV2_Success,
          "descriptor fixture");
    return p;
}
// policy 0: original source-independent entrypoints; 1: independent options;
// 2: borrowed options. Every public encoder constructor is exercised.
WirehairV2Result Create(unsigned route, unsigned policy, const void* source,
    uint64_t message, uint32_t block, Profile& p, WirehairV2Codec* codec)
{
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = policy == 2 ? WirehairV2EncoderSource_BorrowedImmutable :
        WirehairV2EncoderSource_Independent;
    uint32_t bytes = 0;
    if (route == 2) {
        return policy ? wirehair_v2_encoder_create_profile_with_options(
            source, p.data(), static_cast<uint32_t>(p.size()), &options, codec) :
            wirehair_v2_encoder_create_profile(source, p.data(), static_cast<uint32_t>(p.size()), codec);
    }
    WirehairV2Result result;
    if (route == 1) {
        result = policy ? wirehair_v2_encoder_create_profile_id_with_options(
            ProfileId, source, message, block, &options,
            p.data(), static_cast<uint32_t>(p.size()), &bytes, codec) :
            wirehair_v2_encoder_create_profile_id(ProfileId,
                source, message, block, p.data(), static_cast<uint32_t>(p.size()), &bytes, codec);
    } else {
        result = policy ? wirehair_v2_encoder_create_with_options(
            source, message, block, &options, p.data(), static_cast<uint32_t>(p.size()), &bytes, codec) :
            wirehair_v2_encoder_create(source, message, block, p.data(), static_cast<uint32_t>(p.size()), &bytes, codec);
    }
    Check(bytes == static_cast<uint32_t>(p.size()), "constructor required descriptor bytes");
    return result;
}
std::vector<Byte> Message(size_t n)
{
    std::vector<Byte> v(n);
    for (size_t i = 0; i < n; ++i) v[i] = static_cast<Byte>(37 * i + i / 11);
    return v;
}
void Packet(WirehairV2Codec h, uint32_t id, const std::vector<Byte>& expected)
{
    std::vector<Byte> out(expected.size() + 2, 0xa5);
    uint32_t bytes = 999;
    Start(0);
    const auto result = wirehair_v2_encode(h, id, out.data() + 1,
        static_cast<uint32_t>(expected.size()), &bytes);
    const size_t count = Stop();
    Check(result == WirehairV2_Success && count == 0, "allocation-free encode");
    Check(bytes == expected.size() && out.front() == 0xa5 && out.back() == 0xa5 &&
        std::equal(expected.begin(), expected.end(), out.begin() + 1), "packet oracle and guards");
}
void Recover(WirehairV2Codec h, const std::vector<Byte>& expected)
{
    std::vector<Byte> out(expected.size() + 2, 0xa5);
    uint64_t bytes = 999;
    Start(0);
    auto result = wirehair_v2_recover(h, out.data() + 1, expected.size(), &bytes);
    size_t count = Stop();
    Check(result == WirehairV2_Success && count == 0, "allocation-free recovery");
    Check(bytes == expected.size() && out.front() == 0xa5 && out.back() == 0xa5 &&
        std::equal(expected.begin(), expected.end(), out.begin() + 1), "recovery bytes and guards");
}
void LookupCoverage(const Oracle& oracle)
{
    // Each systematic block is a unit vector, so an encoded packet exposes
    // the complete coefficient row. Cover every packed chunk index/phase,
    // plus mixed chunk values, without using the implementation's lookup.
    std::vector<Byte> identity(K * K, 0);
    for (unsigned i = 0; i < K; ++i) identity[i * K + i] = 1;
    Profile p = Descriptor(identity.size(), K);
    WirehairV2Codec h = nullptr;
    Check(Create(FirstRoute, 0, identity.data(), identity.size(), K, p, &h) ==
        WirehairV2_Success, "lookup coverage encoder");
    size_t count = 0;
    const auto check = [&](uint32_t id) {
        const Row expected = oracle.Coefficients(id);
        Packet(h, id, std::vector<Byte>(expected.begin(), expected.end()));
        ++count;
    };
    for (uint32_t low = 0; low < 1024; ++low) {
        check(low);
        check(low | (1u << 10));
    }
    for (unsigned shift : {10u, 17u}) for (uint32_t index = 0; index < 128; ++index)
        for (unsigned phase = 0; phase < 2; ++phase) for (unsigned column = 0; column < K; ++column)
            check((index << shift) | (phase << (shift + 7)) | column);
    for (uint32_t high = 0; high < 256; ++high) for (unsigned column = 0; column < K; ++column)
        check((high << 24) | column);
    uint32_t state = 0x52d397a1u;
    for (unsigned i = 0; i < 512; ++i) {
        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        check(state);
    }
    wirehair_v2_free(h);
    std::cout << "K" << K << " packed lookup oracle: " << count << " packets\n";
}
void DecoderErrors(const Oracle& oracle)
{
    const auto original = Message(K * 64 - 1);
    const Profile p = Descriptor(original.size(), 64);
    WirehairV2Codec h = nullptr;
    Check(wirehair_v2_decoder_create(p.data(), 32, &h) == WirehairV2_Success, "error fixture decoder");
    auto good = oracle.Packet(original, 64, 0);
    auto bad = good; bad[0] ^= 1;
    Check(wirehair_v2_decode(h, 0, good.data(), 64) == WirehairV2_NeedMore, "first pivot");
    Check(wirehair_v2_decode(h, 0, bad.data(), 64) == WirehairV2_Error, "conflict before full rank");
    std::vector<Byte> untouched(original.size(), 0xa5);
    const auto before = untouched;
    Check(wirehair_v2_recover(h, untouched.data(), untouched.size(), nullptr) == WirehairV2_NeedMore &&
        untouched == before, "conflict preserves incomplete recovery");
    Check(wirehair_v2_decode(h, 1, good.data(), 63) == WirehairV2_InvalidInput &&
        wirehair_v2_decode(h, 1, good.data(), 65) == WirehairV2_InvalidInput &&
        wirehair_v2_decode(h, 1, nullptr, 64) == WirehairV2_InvalidInput, "invalid packets do not advance rank");
    Check(wirehair_v2_decode(h, 0, good.data(), 64) == WirehairV2_NeedMore, "good duplicate after conflict");
    for (unsigned id = 1; id < K; ++id) {
        const auto packet = oracle.Packet(original, 64, id);
        Check(wirehair_v2_decode(h, id, packet.data(), static_cast<uint32_t>(packet.size())) ==
            (id == K - 1 ? WirehairV2_Success : WirehairV2_NeedMore), "resume after conflict");
    }
    Check(wirehair_v2_decode(h, 0, bad.data(), 64) == WirehairV2_Error, "full-rank conflict before recover");
    Recover(h, original);
    Recover(h, original);
    wirehair_v2_free(h);
}
void Lifecycle(const Oracle& oracle)
{
    const uint32_t ids[] = {0,1,2,3,4,5,6,7,8,1023,1024,1025,131071,131072,
        16777215,16777216,0x80000000u,0xffffffffu,0xfffffffdu,0xfffffffbu};
    size_t shapes = 0;
    for (uint32_t block : {1u,2u,3u,7u,8u,15u,16u,17u,31u,32u,33u,63u,64u,65u,
            127u,128u,129u,255u,256u,257u,1279u,1280u,1281u,4096u}) {
        std::vector<uint32_t> tails = {1, (block + 1) / 2, block};
        tails.erase(std::unique(tails.begin(), tails.end()), tails.end());
        for (uint32_t tail : tails) {
            ++shapes;
            const auto original = Message(size_t(block) * (K - 1) + tail);
            const Profile canonical = Descriptor(original.size(), block);
            size_t independent_count = 0;
            for (unsigned route = FirstRoute; route < 3; ++route) for (unsigned policy = 0; policy < 3; ++policy) {
                auto source = original;
                Profile p = canonical;
                WirehairV2Codec h = nullptr;
                Start();
                auto result = Create(route, policy, source.data(), source.size(), block, p, &h);
                size_t count = Stop();
                Check(result == WirehairV2_Success && h && p == canonical, "all constructors / identity");
                if (policy == 0) independent_count = count;
                else Check(count == independent_count, "borrowing adds no allocation");
                if (policy != 2) std::fill(source.begin(), source.end(), 0xcc);
                for (uint32_t id : ids) Packet(h, id, oracle.Packet(original, block, id));
                Start(0);
                result = wirehair_v2_encoder_detach_input(h);
                Check(wirehair_v2_encoder_detach_input(h) == WirehairV2_Success, "idempotent detach");
                count = Stop();
                Check(result == WirehairV2_Success && count == 0, "allocation-free detach");
                std::vector<Byte>().swap(source);
                for (uint32_t id : ids) Packet(h, id, oracle.Packet(original, block, id));
                wirehair_v2_free(h);
            }
            // Literal descriptor and independent oracle packets: no live encoder.
            for (unsigned stream = 0; stream < 3; ++stream) {
                WirehairV2Codec decoder = nullptr;
                Check(wirehair_v2_decoder_create(canonical.data(), canonical.size(), &decoder) ==
                    WirehairV2_Success, "standalone decoder");
                Byte untouched = 0x5a;
                Check(wirehair_v2_recover(decoder, &untouched, 0, nullptr) ==
                    WirehairV2_BufferTooSmall && untouched == 0x5a, "short recovery unchanged");
                for (unsigned i = 0; i < K; ++i) {
                    const uint32_t id = stream == 0 ? i : stream == 1 ? K + i : UINT32_MAX - 2 * i;
                    const auto packet = oracle.Packet(original, block, id);
                    Start(0);
                    auto result = wirehair_v2_decode(decoder, id, packet.data(), static_cast<uint32_t>(packet.size()));
                    Check(wirehair_v2_decode(decoder, id, packet.data(), static_cast<uint32_t>(packet.size())) == result,
                        "duplicate packet");
                    size_t count = Stop();
                    Check(count == 0 && result == (i == K - 1 ? WirehairV2_Success : WirehairV2_NeedMore),
                        "first-success / allocation-free feed");
                }
                Recover(decoder, original);
                Recover(decoder, original);
                auto good = oracle.Packet(original, block, UINT32_MAX);
                auto bad = good; bad[0] ^= 1;
                Check(wirehair_v2_decode(decoder, UINT32_MAX, bad.data(), static_cast<uint32_t>(bad.size())) ==
                    WirehairV2_Error, "contradiction after recovery");
                Check(wirehair_v2_decode(decoder, UINT32_MAX, good.data(), static_cast<uint32_t>(good.size())) ==
                    WirehairV2_Success, "conflict preserves retained basis");
                Recover(decoder, original);
                Check(wirehair_v2_encoder_detach_input(decoder) == WirehairV2_InvalidInput,
                    "decoder detach rejected");
                wirehair_v2_free(decoder);
            }
        }
    }
    std::cout << "K" << K << " WHV2 lifecycles: " << shapes << " width/tail shapes x "
              << (3 - FirstRoute) * 3 << " constructors\n";
}
void PivotOrders(const Oracle& oracle)
{
    // Especially for K8: exercise pivot-mask bit7 before the lower pivots,
    // and rotated/reversed RHS/scratch orderings with full and one-byte tails.
    size_t cases = 0;
    for (uint32_t block : {2u, 64u, 1280u}) for (uint32_t tail : {1u, block})
        for (unsigned reverse = 0; reverse < 2; ++reverse)
            for (unsigned rotation = 0; rotation < K; ++rotation) {
                const auto original = Message(size_t(block) * (K - 1) + tail);
                const auto p = Descriptor(original.size(), block);
                WirehairV2Codec h = nullptr;
                Check(wirehair_v2_decoder_create(p.data(), 32, &h) == WirehairV2_Success,
                    "pivot order receiver");
                std::vector<Byte> untouched(original.size(), 0xa5);
                const auto before = untouched;
                for (unsigned i = 0; i < K; ++i) {
                    const unsigned position = (rotation + i) % K;
                    const uint32_t id = reverse ? K - 1 - position : position;
                    const auto good = oracle.Packet(original, block, id);
                    auto bad = good; bad.back() ^= 1;
                    const auto expected = i == K - 1 ? WirehairV2_Success : WirehairV2_NeedMore;
                    Start(0);
                    const auto fed = wirehair_v2_decode(h, id, good.data(), static_cast<uint32_t>(good.size()));
                    const auto conflict = wirehair_v2_decode(h, id, bad.data(), static_cast<uint32_t>(bad.size()));
                    const auto duplicate = wirehair_v2_decode(h, id, good.data(), static_cast<uint32_t>(good.size()));
                    const size_t count = Stop();
                    Check(fed == expected && conflict == WirehairV2_Error && duplicate == expected && count == 0,
                        "pivot permutation / nonpoisoning conflict / allocation-free feed");
                    if (i != K - 1) {
                        Check(wirehair_v2_recover(h, untouched.data(), untouched.size(), nullptr) ==
                            WirehairV2_NeedMore && untouched == before, "incomplete permuted basis no-write");
                    }
                }
                Recover(h, original);
                Recover(h, original);
                for (uint32_t id : {K - 1, K, UINT32_MAX}) {
                    const auto good = oracle.Packet(original, block, id);
                    auto bad = good; bad.back() ^= 1;
                    Check(wirehair_v2_decode(h, id, bad.data(), static_cast<uint32_t>(bad.size())) ==
                        WirehairV2_Error, "post-recovery highest-pivot/repair conflict");
                    Check(wirehair_v2_decode(h, id, good.data(), static_cast<uint32_t>(good.size())) ==
                        WirehairV2_Success, "post-recovery identity retains basis");
                    Recover(h, original);
                }
                wirehair_v2_free(h);
                ++cases;
            }
    std::cout << "K" << K << " pivot order cases: " << cases << '\n';
}
void AllocationFailures()
{
    for (uint32_t tail : {1u, 64u}) for (unsigned route = FirstRoute; route < 3; ++route)
        for (unsigned policy = 0; policy < 3; ++policy) {
            const auto source = Message(64 * (K - 1) + tail);
            const Profile canonical = Descriptor(source.size(), 64);
            Profile p = canonical;
            WirehairV2Codec h = nullptr;
            Start();
            auto result = Create(route, policy, source.data(), source.size(), 64, p, &h);
            const size_t count = Stop();
            Check(result == WirehairV2_Success && count == (tail == 64 ? 3u : 4u),
                "small isolation adds no encoder allocation");
            wirehair_v2_free(h);
            for (size_t fail = 0; fail < count; ++fail) {
                p = canonical;
                if (route != 2) p.fill(0xa5);
                const Profile before = p;
                h = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
                Start(fail);
                result = Create(route, policy, source.data(), source.size(), 64, p, &h);
                Stop();
                Check(result == WirehairV2_OOM && !h && p == before, "every constructor OOM transactional");
            }
        }
    const Profile p = Descriptor(64 * K - 1, 64);
    WirehairV2Codec h = nullptr;
    Start();
    auto result = wirehair_v2_decoder_create(p.data(), static_cast<uint32_t>(p.size()), &h);
    size_t count = Stop();
    Check(result == WirehairV2_Success && count == 3, "decoder allocations");
    wirehair_v2_free(h);
    for (size_t fail = 0; fail < count; ++fail) {
        h = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        Start(fail);
        result = wirehair_v2_decoder_create(p.data(), static_cast<uint32_t>(p.size()), &h);
        Stop();
        Check(result == WirehairV2_OOM && !h, "every decoder OOM transactional");
    }
}
void HandleAllocationIsolation()
{
    // The pre-admission certified handle's layout, using the current target
    // ABI rather than a hard-coded 64-bit allocation size. Never instantiate
    // or reinterpret a live handle as this test-only model.
    struct PriorCertifiedLayout {
        wirehair_v2::Codec Impl;
        uint64_t MessageBytes;
        uint32_t BlockBytes;
        int Mode;
        bool Decoded;
        int SourceState;
        const uint8_t* BorrowedSource;
    };
    for (unsigned k : {2u, 3u, 4u, 5u, 6u, 8u, 128u}) for (unsigned policy : {1u, 2u}) {
        const auto source = Message(size_t(k) * 64 - 1);
        Profile p = {};
        WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
        options.source_policy = policy;
        WirehairV2Codec h = nullptr;
        uint32_t bytes = 0;
        Start();
        auto result = wirehair_v2_encoder_create_profile_id_with_options(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07, source.data(), source.size(),
            64, &options, p.data(), 32, &bytes, &h);
        const size_t encoder_count = Stop();
        Check(result == WirehairV2_Success && h && bytes == 32 && encoder_count != 0 &&
            first_allocation_bytes == sizeof(PriorCertifiedLayout),
            "certified encoder retains pre-admission handle size");
        wirehair_v2_free(h);
        h = nullptr;
        Start();
        result = wirehair_v2_decoder_create(p.data(), 32, &h);
        const size_t decoder_count = Stop();
        Check(result == WirehairV2_Success && h && decoder_count != 0 &&
            first_allocation_bytes == sizeof(PriorCertifiedLayout),
            "certified decoder retains pre-admission handle size");
        wirehair_v2_free(h);
    }
    for (unsigned k : {3u, 5u, 8u}) {
        const auto source = Message(size_t(k) * 64 - 1);
        Profile p = {};
        WirehairV2Codec h = nullptr;
        uint32_t bytes = 0;
        Start();
        const auto result = wirehair_v2_encoder_create_profile_id(
            k == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
                k == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09,
            source.data(), source.size(), 64, p.data(), 32, &bytes, &h);
        Stop();
        Check(result == WirehairV2_Success && h &&
            first_allocation_bytes == sizeof(PriorCertifiedLayout) + 3 * sizeof(void*),
            "small handles retain isolated K3 size");
        wirehair_v2_free(h);
    }
    wirehair_v2_free(nullptr);
}
void Contracts(const Oracle& oracle)
{
    auto source = Message(64 * K);
    const auto original = source;
    Profile p = Descriptor(source.size(), 64);
    WirehairV2Codec h = nullptr;
    uint32_t bytes = 777;
    // Independent descriptor/message overlap is still supported.
    const auto alias_result = K == 3 ?
        wirehair_v2_encoder_create(source.data(), source.size(), 64, source.data(), 32, &bytes, &h) :
        wirehair_v2_encoder_create_profile_id(ProfileId, source.data(), source.size(), 64,
            source.data(), 32, &bytes, &h);
    Check(alias_result == WirehairV2_Success, "staged descriptor/message alias");
    Check(std::equal(p.begin(), p.end(), source.begin()), "overlapping descriptor output");
    Packet(h, UINT32_MAX, oracle.Packet(original, 64, UINT32_MAX));
    wirehair_v2_free(h);
    source = original;
    WirehairV2EncoderOptions opt = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    opt.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    h = reinterpret_cast<WirehairV2Codec>(uintptr_t(1)); bytes = 777;
    Check(wirehair_v2_encoder_create_with_options(source.data(), source.size(), 64, &opt,
        source.data(), 0, &bytes, &h) == WirehairV2_InvalidInput &&
        source == original && bytes == 777 && h == reinterpret_cast<WirehairV2Codec>(uintptr_t(1)),
        "borrowed descriptor alias before capacity, no writes");
    h = nullptr;
    Check(Create(FirstRoute, 2, source.data(), source.size(), 64, p, &h) == WirehairV2_Success, "borrowed alias fixture");
    Byte out[64] = {};
    for (uint32_t id : {0u, K - 1, UINT32_MAX}) {
        bytes = 777;
        Check(wirehair_v2_encode(h, id, source.data() + 80, 0, &bytes) == WirehairV2_InvalidInput &&
            bytes == 777 && source == original, "packet/source alias before capacity");
        Check(wirehair_v2_encode(h, id, nullptr, 0, reinterpret_cast<uint32_t*>(source.data() + 64)) ==
            WirehairV2_InvalidInput && source == original, "counter/source alias before null output");
        bytes = 777;
        Check(wirehair_v2_encode(h, id, out, 0, &bytes) == WirehairV2_BufferTooSmall &&
            bytes == 64 && out[0] == 0, "short encode no writes");
        std::array<uint32_t, 17> shared = {};
        Check(wirehair_v2_encode(h, id, shared.data(), 0, shared.data() + 1) ==
            WirehairV2_InvalidInput && shared[1] == 0, "packet/counter alias before capacity");
    }
    Check(wirehair_v2_encoder_detach_input(h) == WirehairV2_Success, "detach alias fixture");
    Check(wirehair_v2_encode(h, 0, source.data(), 64, &bytes) == WirehairV2_Success,
        "former source output permitted after detach");
    wirehair_v2_free(h);
    // Validation priority and exact profile geometry/seed domain.
    WirehairV2Profile host = {};
    Check(wirehair_v2_profile_deserialize(p.data(), static_cast<uint32_t>(p.size()), &host) == WirehairV2_Success, "parse small profile");
    host.seed_attempt = 1;
    Check(wirehair_v2_profile_serialize(&host, out, sizeof(out), nullptr) == WirehairV2_BadSeed,
        "small nonzero seed rejected");
    host.block_bytes = 63;
    Check(wirehair_v2_profile_serialize(&host, out, sizeof(out), nullptr) == WirehairV2_InvalidDimensions,
        "geometry precedes seed");
    for (uint64_t retired : {UINT64_C(0xe161ce5d456f9bb7), UINT64_C(0x20a4f27a870612a2),
            UINT64_C(0x5748324b35544d31), UINT64_C(0x5748324b38544d31), UINT64_MAX}) {
        host.profile_id = retired;
        Check(wirehair_v2_profile_serialize(&host, out, sizeof(out), nullptr) ==
            WirehairV2_UnsupportedProfile, "retired profile never reinterpreted");
    }
    p[28] = 1;
    opt.options_version = 99;
    h = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    Check(wirehair_v2_encoder_create_profile_with_options(source.data(), p.data(), static_cast<uint32_t>(p.size()), &opt, &h) ==
        WirehairV2_BadSeed && !h, "descriptor error before options");
    bytes = 777;
    Check(wirehair_v2_encoder_create_with_options(source.data(), source.size(), 64, &opt,
        out, 0, &bytes, &h) == WirehairV2_UnsupportedVersion && bytes == 32 && !h,
        "options error before capacity");
    // Explicit CURRENT remains the certified equation ID at every K.
    Check(WIREHAIR_V2_PROFILE_CURRENT == WIREHAIR_V2_PROFILE_CERTIFIED_2026_07, "stable current alias");
    Check(wirehair_v2_encoder_create_profile_id(WIREHAIR_V2_PROFILE_CURRENT, original.data(),
        original.size(), 64, p.data(), static_cast<uint32_t>(p.size()), &bytes, &h) == WirehairV2_Success, "explicit old K3");
    Check(wirehair_v2_profile_deserialize(p.data(), static_cast<uint32_t>(p.size()), &host) == WirehairV2_Success &&
        host.profile_id == WIREHAIR_V2_PROFILE_CURRENT, "explicit old descriptor unchanged");
    wirehair_v2_free(h);
    // Explicit K5/K8 must not silently alter either ordinary selector.
    for (unsigned k : {5u, 8u}) for (unsigned policy : {0u, 1u, 2u}) {
        const auto ordinary_source = Message(k * 64 - 1);
        Profile ordinary = {};
        h = nullptr;
        Check(Create(0, policy, ordinary_source.data(), ordinary_source.size(), 64,
            ordinary, &h) == WirehairV2_Success, "ordinary K5/K8 still constructs");
        Check(wirehair_v2_profile_deserialize(ordinary.data(), 32, &host) == WirehairV2_Success &&
            host.profile_id == WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            "K5/K8 defaults stay certified until admission gates pass");
        wirehair_v2_free(h);
    }
}
void ProtectedSource(const Oracle& oracle)
{
#if defined(__unix__)
    const long page_size = sysconf(_SC_PAGESIZE);
    Check(page_size > 0, "page size");
    const size_t page = static_cast<size_t>(page_size);
    const size_t n = (K * 1280 + page - 1) / page * page;
    void* memory = mmap(nullptr, n, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    Check(memory != MAP_FAILED, "protected source mmap");
    const auto original = Message(K * 1280 - 1);
    std::memcpy(memory, original.data(), original.size());
    auto p = Descriptor(original.size(), 1280);
    WirehairV2Codec h = nullptr;
    Check(Create(FirstRoute, 2, memory, original.size(), 1280, p, &h) == WirehairV2_Success, "protected source create");
    // Deliberate white-box contract test: public callers must keep input readable
    // until detach. The protection detects any hidden repair or detach read.
    Check(mprotect(memory, n, PROT_NONE) == 0, "protect borrowed source");
    Packet(h, UINT32_MAX, oracle.Packet(original, 1280, UINT32_MAX));
    Start(0);
    auto result = wirehair_v2_encoder_detach_input(h);
    size_t count = Stop();
    Check(result == WirehairV2_Success && count == 0, "unreadable-source allocation-free detach");
    Check(munmap(memory, n) == 0, "release source");
    for (uint32_t id = 0; id <= K; ++id) Packet(h, id, oracle.Packet(original, 1280, id));
    Packet(h, UINT32_MAX, oracle.Packet(original, 1280, UINT32_MAX));
    wirehair_v2_free(h);
#else
    (void)oracle;
#endif
}
void CppOwnership()
{
    using namespace wirehair::v2;
    const auto expected = Message(64 * K - 1);
    auto source = expected;
    SerializedProfile profile;
    Encoder encoder;
    const auto created = K == 3 ? encoder.CreateBorrowed(source.data(), source.size(), 64, profile) :
        encoder.CreateBorrowed(ProfileId, source.data(), source.size(), 64, profile);
    Check(created == WirehairV2_Success, "C++ borrowed create");
    Encoder moved(std::move(encoder));
    Check(!encoder && moved, "C++ move transfers lifetime obligation");
    const auto replacement = Message(64 * K);
    SerializedProfile next;
    Start(0);
    const auto failed = K == 3 ? moved.CreateBorrowed(replacement.data(), replacement.size(), 64, next) :
        moved.CreateBorrowed(ProfileId, replacement.data(), replacement.size(), 64, next);
    Stop();
    Check(failed == WirehairV2_OOM && moved, "C++ failed replacement retains old handle");
    Start(0);
    Check(moved.DetachInput() == WirehairV2_Success, "C++ detach");
    Check(Stop() == 0, "C++ detach allocation-free");
    std::vector<Byte>().swap(source);
    Encoder assigned;
    assigned = std::move(moved);
    Check(assigned && !moved, "C++ move assignment");
    Decoder decoder;
    Check(decoder.Create(profile) == WirehairV2_Success, "C++ decoder");
    Byte packet[64];
    for (uint32_t id = K; id < 2 * K; ++id) {
        uint32_t bytes = 0;
        Check(assigned.Encode(id, packet, sizeof(packet), bytes) == WirehairV2_Success,
            "C++ repair after source release");
        Check(decoder.Decode(id, packet, bytes) == (id == 2 * K - 1 ? WirehairV2_Success : WirehairV2_NeedMore),
            "C++ first success");
    }
    std::vector<Byte> output(expected.size());
    uint64_t bytes = 0;
    Check(decoder.Recover(output.data(), output.size(), bytes) == WirehairV2_Success &&
        bytes == output.size() && output == expected, "C++ recovery");
}
void ProfileBounds()
{
    const Profile valid = Descriptor(2 * K - 1, 2);
    Profile foreign_magic = valid;
    foreign_magic[2] = 'K'; foreign_magic[3] = static_cast<Byte>('0' + K);
    Check(wirehair_v2_profile_validate(foreign_magic.data(), 32) == WirehairV2_InvalidMagic,
        "prototype descriptor magic is not WHV2");
    for (unsigned attempt = 1; attempt < 256; ++attempt) {
        Profile p = valid;
        p[28] = static_cast<Byte>(attempt);
        WirehairV2Codec decoder = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        Start(0);
        const auto validation = wirehair_v2_profile_validate(p.data(), 32);
        const auto create = wirehair_v2_decoder_create(p.data(), 32, &decoder);
        const size_t count = Stop();
        Check(validation == WirehairV2_BadSeed && create == validation && !decoder && count == 0,
            "all nonzero small attempts rejected before allocation");
    }
    for (uint32_t block : {0u, 1u, MaxBlockBytes, MaxBlockBytes + 1, UINT32_MAX}) {
        WirehairV2Profile host = {};
        host.struct_bytes = sizeof(host);
        host.profile_version = WIREHAIR_V2_PROFILE_VERSION;
        host.profile_id = ProfileId;
        host.block_bytes = block;
        host.message_bytes = uint64_t(block) * K;
        Profile p = {};
        Start(0);
        const auto result = wirehair_v2_profile_serialize(&host, p.data(), 32, nullptr);
        const size_t count = Stop();
        Check(result == (block == 1 || block == MaxBlockBytes ? WirehairV2_Success : WirehairV2_InvalidDimensions)
            && count == 0, "exact small block-size bound without allocation");
    }
    for (uint32_t block : {1u, 2u, MaxBlockBytes}) {
        for (uint64_t message : {uint64_t(block) * (K - 1),
                uint64_t(block) * (K - 1) + 1, uint64_t(block) * K,
                uint64_t(block) * K + 1}) {
            WirehairV2Profile host = {};
            host.struct_bytes = sizeof(host);
            host.profile_version = WIREHAIR_V2_PROFILE_VERSION;
            host.profile_id = ProfileId;
            host.block_bytes = block;
            host.message_bytes = message;
            Profile p = {};
            Start(0);
            const auto result = wirehair_v2_profile_serialize(&host, p.data(), 32, nullptr);
            const size_t count = Stop();
            const bool valid_shape = message > uint64_t(block) * (K - 1) && message <= uint64_t(block) * K;
            Check(result == (valid_shape ? WirehairV2_Success : WirehairV2_InvalidDimensions) && count == 0,
                "exact small message geometry without allocation");
        }
    }
    auto p = valid;
    auto* wrapped = reinterpret_cast<const void*>(UINTPTR_MAX - 3);
    WirehairV2Codec encoder = nullptr;
    Start(0);
    const auto result = Create(FirstRoute, 0, wrapped, 2 * K - 1, 2, p, &encoder);
    const size_t count = Stop();
    Check(result == WirehairV2_InvalidInput && !encoder && count == 0,
        "independent wrapped input rejected before allocation/read");
}
}
int main()
{
    Check(wirehair_init() == Wirehair_Success, "GF256 initialization");
    Oracle oracle;
    HandleAllocationIsolation();
    LookupCoverage(oracle);
    DecoderErrors(oracle);
    Lifecycle(oracle);
    PivotOrders(oracle);
    AllocationFailures();
    Contracts(oracle);
    ProtectedSource(oracle);
    CppOwnership();
    ProfileBounds();
    std::cout << "WHV2 K" << K << " correctness passed (not a speed or recovery-rate gate)\n";
}
