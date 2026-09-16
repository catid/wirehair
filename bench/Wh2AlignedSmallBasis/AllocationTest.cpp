// Force every supported raw new[] alignment; verify the owned pointer, exact
// copied view, guards, partial padding, packet parity and every OOM cleanup.
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
constexpr unsigned MaxMessage = 8*1280;
alignas(64) std::uint8_t storage[MaxMessage+256];
bool tracking = false, bad = false;
unsigned calls = 0, fail_at = ~0u, raw_offset = 0;
void* live[4] = {};
std::size_t requested[4] = {}, basis_bytes = 0;
void Check(bool ok, const char* why)
{
    if (!ok) { std::fprintf(stderr,"FAIL: %s\n",why); std::exit(1); }
}
void* Allocate(std::size_t n, bool array)
{
    const unsigned slot = tracking ? calls++ : ~0u;
    if (tracking && slot == fail_at) throw std::bad_alloc();
    if (tracking && (slot >= 4 || (slot == 1 ? !array || n != basis_bytes : slot == 3 ? !array : array))) {
        bad = true; throw std::bad_alloc();
    }
    void* p;
    if (tracking && slot == 1) {
        p = storage+64+raw_offset;
        UNPOISON(p,n);
    } else {
        p = std::malloc(n ? n : 1);
        if (!p) throw std::bad_alloc();
    }
    if (tracking) { live[slot] = p; requested[slot] = n; }
    return p;
}
void Release(void* p) noexcept
{
    if (!p) return;
    for (unsigned i = 0; i < 4; ++i) if (p == live[i]) {
        live[i] = nullptr;
        if (i == 1) { POISON(p,requested[i]); return; }
        std::free(p); return;
    }
    const auto address = reinterpret_cast<std::uintptr_t>(p);
    const auto start = reinterpret_cast<std::uintptr_t>(storage);
    if (address >= start && address-start < sizeof(storage)) { bad = true; return; }
    std::free(p);
}
}

void* operator new(std::size_t n) { return Allocate(n,false); }
void* operator new[](std::size_t n) { return Allocate(n,true); }
void operator delete(void* p) noexcept { Release(p); }
void operator delete[](void* p) noexcept { Release(p); }
void* operator new(std::size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n,false); } catch (...) { return nullptr; } }
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n,true); } catch (...) { return nullptr; } }
void operator delete(void* p, const std::nothrow_t&) noexcept { Release(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { Release(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { Release(p); }
void operator delete[](void* p, std::size_t) noexcept { Release(p); }
#endif

namespace {
std::uint64_t Profile(unsigned k)
{
    return k == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
        k == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09;
}
std::uint32_t Packet(unsigned k, unsigned j)
{ return j < 2*k ? j : UINT32_MAX-2*(j-2*k); }
void Start(unsigned offset, std::size_t bytes, unsigned failure)
{
    Check(!tracking && !bad && !live[0] && !live[1] && !live[2] && !live[3],"previous lifetime");
    Check(offset <= 48 && offset%16 == 0 && bytes <= MaxMessage+63,"arena bounds");
    UNPOISON(storage,sizeof(storage));
    std::memset(storage,173,sizeof(storage));
    POISON(storage,sizeof(storage));
    std::fill(requested,requested+4,0);
    calls = 0; raw_offset = offset; basis_bytes = bytes; fail_at = failure; tracking = true;
}
void NoLive()
{ Check(!bad && !live[0] && !live[1] && !live[2] && !live[3],"all allocations released via original pointer"); }
void Guards(const std::uint8_t* source, std::size_t message, unsigned skip, bool copied)
{
    UNPOISON(storage,sizeof(storage));
    const std::size_t begin = 64+raw_offset+skip;
    for (std::size_t i = 0; i < sizeof(storage); ++i) {
        const auto wanted = copied && i >= begin && i < begin+message ? source[i-begin] : 173;
        Check(storage[i] == wanted,"copy position and exact surrounding guards");
    }
    POISON(storage,sizeof(storage));
    if (live[1]) UNPOISON(live[1],requested[1]);
}
}

int main()
{
    Check(wirehair_init() == Wirehair_Success,"init");
    alignas(64) std::uint8_t source[MaxMessage], expected[24][1280], output[1282];
    for (unsigned i = 0; i < MaxMessage; ++i) source[i] = static_cast<std::uint8_t>(i*37+i/11);
    unsigned constructions = 0, failures = 0, packets = 0;
    for (unsigned k : {3u,5u,8u})
        for (unsigned width : {1u,2u,63u,64u,65u,127u,128u,255u,256u,257u,320u,1024u,1280u})
            for (unsigned tail : {1u,(width+1)/2,width}) {
                const std::size_t message = (k-1)*width+tail;
                for (unsigned policy = 0; policy < 2; ++policy) {
                    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
                    options.source_policy = policy ? WirehairV2EncoderSource_BorrowedImmutable : WirehairV2EncoderSource_Independent;
                    std::uint8_t canonical[32]; std::uint32_t written = 0; WirehairV2Codec oracle = nullptr;
                    Check(wirehair_v2_encoder_create_profile_id_with_options(Profile(k),source,message,width,&options,
                          canonical,32,&written,&oracle) == WirehairV2_Success && oracle && written == 32,"ordinary oracle");
                    for (unsigned j = 0; j < 3*k; ++j)
                        Check(wirehair_v2_encode(oracle,Packet(k,j),expected[j],width,&written) == WirehairV2_Success &&
                              written == (j == k-1 ? tail : width),"oracle packet");
                    wirehair_v2_free(oracle);
                    const unsigned count = tail == width ? 3u : 4u;
                    unsigned extra = 0;
#if defined(ALIGN_CANDIDATE)
                    if (width >= 256 && width%64 == 0) extra = 63;
#endif
                    for (unsigned offset : {0u,16u,32u,48u}) {
                        const unsigned skip = extra && offset ? 64-offset : 0;
                        std::uint8_t descriptor[32]; WirehairV2Codec h = nullptr;
                        Start(offset,message+extra,~0u);
                        const auto result = wirehair_v2_encoder_create_profile_id_with_options(Profile(k),source,message,width,&options,
                            descriptor,32,&written,&h);
                        tracking = false;
                        Check(result == WirehairV2_Success && h && written == 32 && calls == count &&
                              std::memcmp(descriptor,canonical,32) == 0,"controlled constructor");
                        Check(h == live[0] && live[1] == storage+64+offset && live[2] &&
                              bool(live[3]) == (tail != width) &&
                              (tail == width || requested[3] == width),"exact allocator interception");
                        Guards(source,message,skip,true);
                        for (unsigned j = 0; j < 3*k; ++j) {
                            std::memset(output,173,sizeof(output));
                            const unsigned bytes = j == k-1 ? tail : width;
                            calls = 0; fail_at = 0; tracking = true;
                            const auto status = wirehair_v2_encode(h,Packet(k,j),output+1,width,&written);
                            tracking = false;
                            Check(status == WirehairV2_Success && calls == 0 && written == bytes &&
                                  std::memcmp(output+1,expected[j],bytes) == 0 && output[0] == 173 &&
                                  std::all_of(output+1+bytes,output+sizeof(output),[](std::uint8_t x){return x==173;}),"allocation-free packet parity and guards");
                            ++packets;
                        }
                        wirehair_v2_free(h); NoLive(); Guards(source,message,skip,true);
                        ++constructions;
                        for (unsigned failure = 0; failure < count; ++failure) {
                            h = reinterpret_cast<WirehairV2Codec>(std::uintptr_t(1));
                            std::memset(descriptor,211,sizeof(descriptor));
                            Start(offset,message+extra,failure);
                            const auto status = wirehair_v2_encoder_create_profile_id_with_options(Profile(k),source,message,width,&options,
                                descriptor,32,&written,&h);
                            tracking = false;
                            Check(status == WirehairV2_OOM && !h && calls == failure+1 &&
                                  std::all_of(descriptor,descriptor+32,[](std::uint8_t x){return x==211;}),"every OOM boundary is transactional");
                            NoLive(); Guards(source,message,skip,failure > 1);
                            ++failures;
                        }
                    }
                }
            }
    std::printf("PASS %u forced-placement constructions, %u OOM boundaries, %u packets\n",constructions,failures,packets);
}
