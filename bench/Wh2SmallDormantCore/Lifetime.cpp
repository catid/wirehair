// Link-wrap the real out-of-line Codec constructor/destructor only in this
// qualification executable. Real candidate libraries contain no test hooks.
#include <wirehair/wirehair.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#if defined(WH2_BASIS_ASAN)
#include <sanitizer/asan_interface.h>
#define POISON(p,n) __asan_poison_memory_region(p,n)
#define UNPOISON(p,n) __asan_unpoison_memory_region(p,n)
#else
#define POISON(p,n) ((void)0)
#define UNPOISON(p,n) ((void)0)
#endif

namespace {
bool tracking = false, small = false, poisoned = false;
unsigned count = 0, failure = ~0u, constructors = 0, destructors = 0;
void* facade = nullptr;
void* live[256] = {};
constexpr unsigned CoreBytes = 240, SmallBytes = 296;
void Check(bool ok, const char* why)
{ if (!ok) { std::fprintf(stderr,"FAIL: %s\n",why); std::abort(); } }
__attribute__((noinline)) void* Allocate(std::size_t bytes)
{
    const unsigned index = tracking ? count++ : ~0u;
    if (tracking && index == failure) throw std::bad_alloc();
    void* p = std::malloc(bytes ? bytes : 1);
    if (!p) throw std::bad_alloc();
    if (tracking) {
        unsigned slot = index == 0 ? 0 : 1;
        while (slot < 256 && live[slot]) ++slot;
        Check(slot < 256,"simultaneously live allocation registry bound");
        live[slot] = p;
        if (index == 0) {
            Check(bytes == (small ? SmallBytes : 272),"facade size");
            facade = p;
#if defined(DORMANT_CANDIDATE)
            if (small) {
                // Before outer construction; remains poisoned through the
                // complete derived/base destruction and every OOM unwind.
                std::memset(p,0xa5,CoreBytes); POISON(p,CoreBytes); poisoned = true;
            }
#endif
        }
    }
    return p;
}
__attribute__((noinline)) void Release(void* p) noexcept
{
    if (!p) return;
    for (unsigned i = 0; i < 256; ++i) if (live[i] == p) {
        live[i] = nullptr;
        if (i == 0 && poisoned) {
            UNPOISON(p,CoreBytes); poisoned = false;
            const auto* bytes = static_cast<const unsigned char*>(p);
            for (unsigned j = 0; j < CoreBytes; ++j) Check(bytes[j] == 0xa5,"inactive core canary survived entire lifetime");
        }
        break;
    }
    std::free(p);
}
void CoreConstruct(void* p)
{ if (p == facade) { Check(!poisoned,"unexpected dormant-core constructor"); ++constructors; } }
void CoreDestroy(void* p)
{ if (p == facade) { Check(!poisoned,"unexpected dormant-core destructor"); ++destructors; } }
}

#define WRAP(SYMBOL, CHECK) \
extern "C" void __real_##SYMBOL(void*); \
extern "C" void __wrap_##SYMBOL(void* p) { CHECK(p); __real_##SYMBOL(p); }
WRAP(_ZN11wirehair_v25CodecC1Ev, CoreConstruct)
WRAP(_ZN11wirehair_v25CodecC2Ev, CoreConstruct)
WRAP(_ZN11wirehair_v25CodecD1Ev, CoreDestroy)
WRAP(_ZN11wirehair_v25CodecD2Ev, CoreDestroy)
#undef WRAP
void* operator new(std::size_t n) { return Allocate(n); }
void* operator new[](std::size_t n) { return Allocate(n); }
void operator delete(void* p) noexcept { Release(p); }
void operator delete[](void* p) noexcept { Release(p); }
void* operator new(std::size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n); } catch (...) { return nullptr; } }
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n); } catch (...) { return nullptr; } }
void operator delete(void* p, const std::nothrow_t&) noexcept { Release(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { Release(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { Release(p); }
void operator delete[](void* p, std::size_t) noexcept { Release(p); }
#endif

namespace {
void Start(bool is_small, unsigned fail = ~0u)
{
    Check(!tracking && !poisoned && std::all_of(live,live+256,[](void* p){return !p;}),"previous cleanup");
    small = is_small; count = constructors = destructors = 0; failure = fail; facade = nullptr; tracking = true;
}
void Finished()
{
    Check(!tracking && !poisoned && std::all_of(live,live+256,[](void* p){return !p;}),"all constructor allocations freed");
    unsigned expected = facade ? 1u : 0u;
#if defined(DORMANT_CANDIDATE)
    if (small) expected = 0;
#endif
    Check(constructors == expected && destructors == expected,"exact real Codec constructor/destructor counts");
    facade = nullptr;
}
uint64_t Profile(unsigned k, bool is_small)
{
    return !is_small ? WIREHAIR_V2_PROFILE_CERTIFIED_2026_07 :
        k == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
        k == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09;
}
}

int main()
{
    Check(wirehair_init() == Wirehair_Success,"init");
    unsigned char source[16*1280], packet[1280], repairs[48][1280], recovered[16*1280];
    for (unsigned i = 0; i < sizeof(source); ++i) source[i] = static_cast<unsigned char>(i*37+i/11);
    unsigned lifetimes = 0, oom = 0, max_calls = 0;
    for (unsigned flavor = 0; flavor < 5; ++flavor) {
        const unsigned k = flavor == 0 || flavor == 3 ? 3 : flavor == 1 ? 5 : flavor == 2 ? 8 : 16;
        const bool is_small = flavor < 3;
        for (unsigned route = 0; route < 3; ++route) {
            if (route == 0 && (flavor == 1 || flavor == 2 || flavor == 3)) continue;
        for (unsigned width : {2u,64u,1280u}) for (unsigned tail : {1u,width})
            for (unsigned policy = 0; policy < 3; ++policy) {
                const std::size_t message = std::size_t(k-1)*width+tail;
                WirehairV2EncoderOptions opts = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
                opts.source_policy = policy == 2 ? WirehairV2EncoderSource_BorrowedImmutable : WirehairV2EncoderSource_Independent;
                unsigned char descriptor[32], canonical[32];
                unsigned enc_allocs = 0, dec_allocs = 0;
                uint32_t n = 0;
                // Canonical selected descriptor also supplies the serialized
                // constructor route; never guess a certified attempt number.
                WirehairV2Codec oracle = nullptr;
                Check(wirehair_v2_encoder_create_profile_id_with_options(Profile(k,is_small),source,message,width,&opts,
                      canonical,32,&n,&oracle) == WirehairV2_Success && oracle && n == 32,"descriptor oracle");
                wirehair_v2_free(oracle);
                // All small constructor boundaries; certified positive controls
                // cover facade failure and the first nested allocation only.
                for (unsigned trial = 0; trial <= enc_allocs; ++trial) {
                    WirehairV2Codec h = reinterpret_cast<WirehairV2Codec>(std::uintptr_t(1));
                    if (route == 2) std::memcpy(descriptor,canonical,32);
                    else std::memset(descriptor,0xcc,32);
                    unsigned char before[32]; std::memcpy(before,descriptor,32);
                    Start(is_small,trial == 0 ? ~0u : trial-1);
                    WirehairV2Result status;
                    if (route == 0) status = policy ? wirehair_v2_encoder_create_with_options(source,message,width,&opts,descriptor,32,&n,&h) :
                        wirehair_v2_encoder_create(source,message,width,descriptor,32,&n,&h);
                    else if (route == 1) status = policy ?
                        wirehair_v2_encoder_create_profile_id_with_options(Profile(k,is_small),source,message,width,&opts,descriptor,32,&n,&h) :
                        wirehair_v2_encoder_create_profile_id(Profile(k,is_small),source,message,width,descriptor,32,&n,&h);
                    else status = policy ? wirehair_v2_encoder_create_profile_with_options(source,descriptor,32,&opts,&h) :
                        wirehair_v2_encoder_create_profile(source,descriptor,32,&h);
                    tracking = false;
                    if (!trial) {
                        Check(status == WirehairV2_Success && h == facade && (route == 2 || n == 32) &&
                              std::memcmp(canonical,descriptor,32) == 0,"encoder lifetime");
                        Check(count >= 2 && (!is_small || count == (tail == width ? 3u : 4u)),"encoder allocation count");
                        max_calls = std::max(max_calls,count);
                        enc_allocs = is_small ? count : 2;
                        for (unsigned j = 0; j < k; ++j) Check(wirehair_v2_encode(h,j,packet,width,&n) == WirehairV2_Success &&
                            n == (j == k-1 ? tail : width) && std::memcmp(packet,source+std::size_t(j)*width,n) == 0,"systematic");
                        for (unsigned j = 0; j < 48; ++j) Check(wirehair_v2_encode(h,k+j,repairs[j],width,&n) == WirehairV2_Success && n == width,"repair");
                        Check(wirehair_v2_encode(h,UINT32_MAX,packet,width,&n) == WirehairV2_Success && n == width,"distant repair");
                        Check(wirehair_v2_encoder_detach_input(h) == WirehairV2_Success,"detach");
                        wirehair_v2_free(h); ++lifetimes;
                    } else {
                        Check(status == WirehairV2_OOM && !h && count == trial &&
                              std::memcmp(descriptor,before,32) == 0,"encoder OOM transaction");
                        ++oom;
                    }
                    Finished();
                }
                // Retain the successful encoder's exact selected descriptor;
                // certified selection must not be replaced with assumed attempt0.
                std::memcpy(descriptor,canonical,32);
                for (unsigned trial = 0; trial <= dec_allocs; ++trial) {
                    WirehairV2Codec h = reinterpret_cast<WirehairV2Codec>(std::uintptr_t(1));
                    Start(is_small,trial == 0 ? ~0u : trial-1);
                    const auto status = wirehair_v2_decoder_create(descriptor,32,&h);
                    tracking = false;
                    if (!trial) {
                        Check(status == WirehairV2_Success && h == facade,"decoder lifetime");
                        Check(count >= 2 && (!is_small || count == 3),"decoder allocation count");
                        max_calls = std::max(max_calls,count);
                        dec_allocs = is_small ? count : 2;
                        bool done = false;
                        for (unsigned j = 0; j < 48; ++j) {
                            const auto result = wirehair_v2_decode(h,k+j,repairs[j],width);
                            Check(result == WirehairV2_Success || result == WirehairV2_NeedMore,"repair feed");
                            if (result == WirehairV2_Success) { done = true; break; }
                        }
                        Check(done,"bounded first success");
                        for (unsigned repeat = 0; repeat < 2; ++repeat) {
                            uint64_t bytes = 0;
                            Check(wirehair_v2_recover(h,recovered,message,&bytes) == WirehairV2_Success && bytes == message &&
                                  std::memcmp(recovered,source,message) == 0,"recovery");
                        }
                        wirehair_v2_free(h); ++lifetimes;
                    } else {
                        Check(status == WirehairV2_OOM && !h && count == trial,"decoder OOM transaction");
                        ++oom;
                    }
                    Finished();
                }
            }
        }
    }
    Check(lifetimes == 432 && oom == 1179,"complete lifetime roster");
    std::printf("PASS %u lifetimes; %u OOM checks (all small, first two certified); max constructor calls %u\n",lifetimes,oom,max_calls);
}
