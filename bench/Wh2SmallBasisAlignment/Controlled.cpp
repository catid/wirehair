// Baseline-only, hot-repair alignment diagnosis. No candidate or WH1 speed gate.
#include <wirehair/wirehair.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <new>
#include <sched.h>
#include <stdexcept>
#include <vector>
#if defined(WH2_BASIS_ASAN)
#include <sanitizer/asan_interface.h>
#define POISON(p,n) __asan_poison_memory_region(p,n)
#define UNPOISON(p,n) __asan_unpoison_memory_region(p,n)
#else
#define POISON(p,n) ((void)0)
#define UNPOISON(p,n) ((void)0)
#endif

namespace arena {
constexpr unsigned CarrierBytes = 16384;
alignas(4096) std::uint8_t basis[2][CarrierBytes];
alignas(64) std::uint8_t object[512], evaluator[512];
void* live[3] = {};
std::size_t extents[3] = {};
unsigned carrier = 0, offset = 0, calls = 0, work_allocations = 0;
std::size_t message_bytes = 0;
bool active = false, working = false, bad = false;

void* Allocate(std::size_t bytes, bool array)
{
    if (working) ++work_allocations;
    if (!active) {
        void* p = std::malloc(bytes ? bytes : 1);
        if (!p) throw std::bad_alloc();
        return p;
    }
    const unsigned index = calls++;
    if (index >= 3 || (index == 1) != array || (index == 1 ? bytes != message_bytes : bytes > 512)) {
        bad = true; throw std::bad_alloc();
    }
    void* p = index == 0 ? static_cast<void*>(object) : index == 2 ? static_cast<void*>(evaluator) :
        static_cast<void*>(basis[carrier] + 4096 + offset);
    if (live[index]) { bad = true; throw std::bad_alloc(); }
    live[index] = p; extents[index] = bytes;
    UNPOISON(p,bytes);
    return p;
}
bool In(const void* p, const void* start, std::size_t bytes)
{
    const auto address = reinterpret_cast<std::uintptr_t>(p);
    const auto base = reinterpret_cast<std::uintptr_t>(start);
    return address >= base && address - base < bytes;
}
void Release(void* p) noexcept
{
    for (unsigned i = 0; i < 3; ++i) if (p && p == live[i]) { POISON(p,extents[i]); live[i] = nullptr; return; }
    if (In(p,object,sizeof(object)) || In(p,evaluator,sizeof(evaluator)) || In(p,basis,sizeof(basis))) {
        bad = true; return; // Never pass an interior or double-freed arena pointer to free.
    }
    std::free(p);
}
}

void* operator new(std::size_t n) { return arena::Allocate(n,false); }
void* operator new[](std::size_t n) { return arena::Allocate(n,true); }
void operator delete(void* p) noexcept { arena::Release(p); }
void operator delete[](void* p) noexcept { arena::Release(p); }
void* operator new(std::size_t n, const std::nothrow_t&) noexcept
{ try { return arena::Allocate(n,false); } catch (...) { return nullptr; } }
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept
{ try { return arena::Allocate(n,true); } catch (...) { return nullptr; } }
void operator delete(void* p, const std::nothrow_t&) noexcept { arena::Release(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { arena::Release(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { arena::Release(p); }
void operator delete[](void* p, std::size_t) noexcept { arena::Release(p); }
#endif

namespace {
constexpr unsigned Batch = 64, Reps = 12, MaxMessage = 8*1280;
const unsigned Sides[18] = {0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
alignas(64) std::uint8_t source[MaxMessage], expected[2*MaxMessage];
alignas(64) std::uint8_t output[Batch*2*MaxMessage+128];
void Check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> F Symbol(void* lib, const char* name, void* owner)
{
    void* p = dlsym(lib,name); Dl_info info = {};
    Check(p && dladdr(p,&info) && info.dli_fbase == owner,"DSO owner");
    F f; static_assert(sizeof(f) == sizeof(p),"Linux function pointer size");
    std::memcpy(&f,&p,sizeof(f)); return f;
}
struct Library {
    void* lib;
    decltype(&wirehair_v2_encoder_create_profile_id_with_options) create;
    decltype(&wirehair_v2_encode) encode;
    decltype(&wirehair_v2_free) free;
    explicit Library(const char* path): lib(dlopen(path,RTLD_NOW | RTLD_LOCAL))
    {
        Check(lib,"dlopen"); Dl_info owner = {};
        Check(dladdr(dlsym(lib,"wirehair_init_"),&owner),"DSO base");
        Check(Symbol<decltype(&wirehair_init_)>(lib,"wirehair_init_",owner.dli_fbase)(WIREHAIR_VERSION) == Wirehair_Success,"init");
        create = Symbol<decltype(create)>(lib,"wirehair_v2_encoder_create_profile_id_with_options",owner.dli_fbase);
        encode = Symbol<decltype(encode)>(lib,"wirehair_v2_encode",owner.dli_fbase);
        free = Symbol<decltype(free)>(lib,"wirehair_v2_free",owner.dli_fbase);
    }
    ~Library() { dlclose(lib); }
};
struct Shape { unsigned k, bytes, policy; };
Shape ShapeOf(unsigned cell) { return Shape{cell < 4 ? 5u : 8u, (cell/2)%2 ? 1280u : 64u,cell%2}; }
std::uint32_t Packet(unsigned k, unsigned j) { return j < k ? k+j : UINT32_MAX-2*(j-k); }
std::uint64_t Now()
{ return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

WirehairV2Codec Create(Library& lib, Shape s)
{
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = s.policy ? WirehairV2EncoderSource_BorrowedImmutable : WirehairV2EncoderSource_Independent;
    std::uint8_t profile[32]; std::uint32_t written = 0; WirehairV2Codec h = nullptr;
    const auto result = lib.create(s.k == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09,
        source,s.k*s.bytes,s.bytes,&options,profile,32,&written,&h);
    Check(result == WirehairV2_Success && h && written == 32,"create");
    return h;
}
void PrepareExpected(Library& lib, Shape s)
{
    // Ordinary malloc-backed handle; identical public packets are the oracle
    // for this address-only intervention, not an independent field proof.
    WirehairV2Codec h = Create(lib,s);
    for (unsigned j = 0; j < 2*s.k; ++j) {
        std::uint32_t written = 0;
        Check(lib.encode(h,Packet(s.k,j),expected+j*s.bytes,s.bytes,&written) == WirehairV2_Success && written == s.bytes,"oracle repair");
    }
    lib.free(h);
}
void PrepareArena(Shape s, unsigned carrier, unsigned offset)
{
    Check(!arena::active && !arena::working && !arena::bad &&
          !arena::live[0] && !arena::live[1] && !arena::live[2],"arena lifetime");
    Check(carrier < 2 && offset <= 48 && offset%16 == 0,"placement");
    arena::carrier = carrier; arena::offset = offset; arena::message_bytes = s.k*s.bytes;
    Check(4096+offset+arena::message_bytes < arena::CarrierBytes,"carrier extent");
    UNPOISON(arena::basis,sizeof(arena::basis));
    UNPOISON(arena::object,sizeof(arena::object));
    UNPOISON(arena::evaluator,sizeof(arena::evaluator));
    std::memset(arena::basis,173,sizeof(arena::basis));
    std::memset(arena::object,173,sizeof(arena::object));
    std::memset(arena::evaluator,173,sizeof(arena::evaluator));
    std::memset(output,173,sizeof(output));
    arena::calls = 0; arena::work_allocations = 0;
    POISON(arena::basis,sizeof(arena::basis));
    POISON(arena::object,sizeof(arena::object));
    POISON(arena::evaluator,sizeof(arena::evaluator));
}
void CheckGuards()
{
    UNPOISON(arena::basis,sizeof(arena::basis));
    UNPOISON(arena::object,sizeof(arena::object));
    UNPOISON(arena::evaluator,sizeof(arena::evaluator));
    for (unsigned slot = 0; slot < 2; ++slot) {
        const unsigned begin = 4096+arena::offset;
        for (unsigned i = 0; i < arena::CarrierBytes; ++i) {
            const std::uint8_t wanted = slot == arena::carrier && i >= begin && i < begin+arena::message_bytes ?
                source[i-begin] : 173;
            Check(arena::basis[slot][i] == wanted,"basis contents and guards");
        }
    }
    for (unsigned i = 0; i < 512; ++i) {
        if (i >= arena::extents[0]) Check(arena::object[i] == 173,"handle guard");
        if (i >= arena::extents[2]) Check(arena::evaluator[i] == 173,"evaluator guard");
    }
    POISON(arena::basis,sizeof(arena::basis));
    POISON(arena::object,sizeof(arena::object));
    POISON(arena::evaluator,sizeof(arena::evaluator));
    for (unsigned i = 0; i < 3; ++i) if (arena::live[i]) UNPOISON(arena::live[i],arena::extents[i]);
}
__attribute__((noinline)) void Work(Library& lib, WirehairV2Codec h, Shape s)
{
    for (unsigned cycle = 0; cycle < Batch; ++cycle) for (unsigned j = 0; j < 2*s.k; ++j) {
        std::uint32_t written = 0;
        Check(lib.encode(h,Packet(s.k,j),output+64+(cycle*2*s.k+j)*s.bytes,s.bytes,&written) == WirehairV2_Success &&
              written == s.bytes,"repair work");
    }
}
void Verify(Shape s)
{
    const unsigned bytes = 2*s.k*s.bytes;
    for (unsigned cycle = 0; cycle < Batch; ++cycle)
        Check(std::memcmp(output+64+cycle*bytes,expected,bytes) == 0,"repair parity");
    Check(std::all_of(output,output+64,[](std::uint8_t x){return x==173;}) &&
          std::all_of(output+64+Batch*bytes,output+sizeof(output),[](std::uint8_t x){return x==173;}),"output guards");
    CheckGuards();
}
struct Row { unsigned rep,cell,pair,order,pos,side,carrier,offset; std::uint64_t ns; };
std::vector<Row> rows;
bool published = false;
void Publish()
{
    if (published) return;
    published = true;
    std::puts("rep,cell,pair,order,position,side,carrier,offset,ns,handle,evaluator,basis,source,output");
    for (const auto& r : rows)
        std::printf("%u,%u,%u,%u,%u,%u,%u,%u,%llu,%llu,%llu,%llu,%llu,%llu\n",r.rep,r.cell,r.pair,r.order,r.pos,r.side,r.carrier,r.offset,
            static_cast<unsigned long long>(r.ns),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(arena::object)),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(arena::evaluator)),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(arena::basis[r.carrier]+4096+r.offset)),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(source)),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(output+64)));
    Check(!std::ferror(stdout) && std::fflush(stdout) == 0,"output");
}
void Placement(unsigned pair, unsigned side, unsigned& carrier, unsigned& offset)
{
    if (pair < 4) { carrier = side; offset = pair*16; } // Same alignment, different storage slots.
    else { carrier = (pair-4)/3; offset = side ? ((pair-4)%3+1)*16 : 0; }
}
void Run(Library& lib, bool timing)
{
    for (unsigned i = 0; i < MaxMessage; ++i) source[i] = static_cast<std::uint8_t>(i*37+i/11);
    rows.reserve(Reps*8*10*2*18);
    const auto start = Now();
    for (unsigned rep = 0; rep < (timing ? Reps : 1u); ++rep)
        for (unsigned slot = 0; slot < 8; ++slot) {
            const unsigned cell = (slot+rep*3)%8; const Shape s = ShapeOf(cell);
            PrepareExpected(lib,s);
            for (unsigned ps = 0; ps < 10; ++ps) for (unsigned order = 0; order < 2; ++order) {
                const unsigned pair = (ps+rep+cell)%10;
                for (unsigned pos = 0; pos < (timing ? 18u : 2u); ++pos) {
                    const unsigned side = Sides[pos]^order;
                    unsigned carrier,offset; Placement(pair,side,carrier,offset);
                    PrepareArena(s,carrier,offset);
                    arena::active = true;
                    WirehairV2Codec h = Create(lib,s);
                    arena::active = false;
                    Check(!arena::bad && arena::calls == 3 && static_cast<void*>(h) == arena::object &&
                          arena::live[1] == arena::basis[carrier]+4096+offset && arena::live[2] == arena::evaluator,
                          "exact allocation interposition");
                    CheckGuards();
                    arena::working = true;
                    const auto begin = Now(); Work(lib,h,s); const auto end = Now();
                    arena::working = false;
                    Verify(s);
                    Check(arena::work_allocations == 0,"allocation-free work");
                    lib.free(h);
                    Check(!arena::bad && !arena::live[0] && !arena::live[1] && !arena::live[2],"exact deallocation interposition");
                    CheckGuards();
                    if (timing) rows.push_back(Row{rep,cell,pair,order,pos,side,carrier,offset,end-begin});
                    Check(Now()-start < UINT64_C(120000000000),"120-second cap");
                }
            }
        }
}
}

int main(int argc, char** argv)
{
    try {
        Check(argc == 3 && (!std::strcmp(argv[2],"neutral") || !std::strcmp(argv[2],"run")),"arguments");
        const bool timing = !std::strcmp(argv[2],"run");
        cpu_set_t mask; CPU_ZERO(&mask); CPU_SET(50,&mask);
        Check(sched_setaffinity(0,sizeof(mask),&mask) == 0 && sched_getaffinity(0,sizeof(mask),&mask) == 0 &&
              CPU_COUNT(&mask) == 1 && CPU_ISSET(50,&mask) && sched_getcpu() == 50,"CPU50 singleton");
        Library lib(argv[1]);
        Run(lib,timing);
        if (!timing) { std::puts("PASS 320 controlled constructions; repair parity, placements, guards and frees"); return 0; }
        Check(rows.size() == Reps*8*10*2*18,"complete cohort");
        Publish();
    } catch (const std::exception& e) {
        arena::active = false; arena::working = false;
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        try { if (!rows.empty()) Publish(); } catch (...) {}
        return 1;
    }
}
