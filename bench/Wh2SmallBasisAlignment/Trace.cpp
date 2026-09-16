// Non-timing observer of actual, unmodified public-DSO constructor allocations.
#include <wirehair/wirehair.h>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <new>
#include <stdexcept>

namespace {
struct Allocation { void* pointer; std::size_t bytes; bool array; };
Allocation allocations[8];
unsigned allocation_count = 0;
bool tracking = false, overflow = false;

void* Allocate(std::size_t bytes, bool array)
{
    void* p = std::malloc(bytes ? bytes : 1);
    if (!p) throw std::bad_alloc();
    if (tracking) {
        if (allocation_count == 8) overflow = true;
        else allocations[allocation_count++] = Allocation{p, bytes, array};
    }
    return p;
}
}

// Ordinary malloc/free semantics; no address selection or alignment changes.
void* operator new(std::size_t n) { return Allocate(n, false); }
void* operator new[](std::size_t n) { return Allocate(n, true); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void* operator new(std::size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n, false); } catch (...) { return nullptr; } }
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n, true); } catch (...) { return nullptr; } }
void operator delete(void* p, const std::nothrow_t&) noexcept { std::free(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

namespace {
void Check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> F Symbol(void* dso, const char* name, void* owner)
{
    void* p = dlsym(dso, name);
    Dl_info info = {};
    Check(p && dladdr(p, &info) && info.dli_fbase == owner, "symbol owner");
    F f;
    static_assert(sizeof(f) == sizeof(p), "Linux function pointer size");
    std::memcpy(&f, &p, sizeof(f));
    return f;
}
struct Library {
    void* lib;
    decltype(&wirehair_v2_encoder_create_profile_id_with_options) create;
    decltype(&wirehair_v2_encode) encode;
    decltype(&wirehair_v2_free) free;
    explicit Library(const char* path): lib(dlopen(path, RTLD_NOW | RTLD_LOCAL))
    {
        Check(lib, "dlopen");
        Dl_info owner = {};
        Check(dladdr(dlsym(lib,"wirehair_init_"), &owner), "DSO base");
        Check(Symbol<decltype(&wirehair_init_)>(lib,"wirehair_init_",owner.dli_fbase)(WIREHAIR_VERSION) == Wirehair_Success, "init");
        create = Symbol<decltype(create)>(lib,"wirehair_v2_encoder_create_profile_id_with_options",owner.dli_fbase);
        encode = Symbol<decltype(encode)>(lib,"wirehair_v2_encode",owner.dli_fbase);
        free = Symbol<decltype(free)>(lib,"wirehair_v2_free",owner.dli_fbase);
    }
    ~Library() { dlclose(lib); }
    Library(const Library&) = delete;
    Library& operator=(const Library&) = delete;
};
struct Handle {
    Library* library;
    WirehairV2Codec value = nullptr;
    explicit Handle(Library* l): library(l) {}
    ~Handle() { if (value) library->free(value); }
};
struct Row {
    unsigned k, width, policy, creation_order, arm, allocation;
    std::size_t bytes;
    std::uintptr_t pointer;
    bool array;
};
Row rows[144];
unsigned row_count = 0;
std::uint64_t ProfileId(unsigned k)
{
    return k == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
        k == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09;
}
void Trace(Library* baseline, Library* candidate)
{
    alignas(64) std::uint8_t source[8*1280], packets[2][1280];
    for (std::size_t i = 0; i < sizeof(source); ++i) source[i] = static_cast<std::uint8_t>(i*37 + i/11);
    for (unsigned k : {3u,5u,8u}) for (unsigned width : {64u,1280u})
        for (unsigned policy = 0; policy < 2; ++policy)
            for (unsigned order = 0; order < 2; ++order) {
                Handle a(baseline), b(candidate);
                Handle* handles[] = {&a,&b};
                std::array<std::uint8_t,32> profiles[2] = {};
                void* bases[2] = {};
                for (unsigned slot = 0; slot < 2; ++slot) {
                    const unsigned arm = slot ^ order;
                    Handle& h = *handles[arm];
                    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
                    options.source_policy = policy ? WirehairV2EncoderSource_BorrowedImmutable : WirehairV2EncoderSource_Independent;
                    std::uint32_t written = 0;
                    allocation_count = 0; overflow = false; tracking = true;
                    const auto result = h.library->create(ProfileId(k), source, k*width, width, &options,
                        profiles[arm].data(), 32, &written, &h.value);
                    tracking = false;
                    Check(result == WirehairV2_Success && h.value && written == 32, "constructor");
                    Check(!overflow && allocation_count == 3, "exact three intercepted allocations");
                    Check(allocations[0].pointer == h.value && !allocations[0].array &&
                          allocations[1].array && allocations[1].bytes == k*width &&
                          !allocations[2].array, "handle/basis/evaluator allocation order");
                    bases[arm] = allocations[1].pointer;
                    Check(bases[arm] != source && std::memcmp(bases[arm],source,k*width) == 0, "private basis copy");
                    for (unsigned i = 0; i < 3; ++i) {
                        Check(row_count < 144, "row cap");
                        rows[row_count++] = Row{k,width,policy,order,arm,i,allocations[i].bytes,
                            reinterpret_cast<std::uintptr_t>(allocations[i].pointer),allocations[i].array};
                    }
                }
                Check(bases[0] != bases[1] && profiles[0] == profiles[1], "separate private bases, same equations");
                for (unsigned j = 0; j < 2*k; ++j) {
                    const std::uint32_t id = j < k ? k+j : UINT32_MAX-2*(j-k);
                    for (unsigned arm = 0; arm < 2; ++arm) {
                        std::uint32_t written = 0;
                        Handle& h = *handles[arm];
                        Check(h.library->encode(h.value,id,packets[arm],width,&written) == WirehairV2_Success &&
                              written == width, "repair");
                    }
                    Check(std::memcmp(packets[0],packets[1],width) == 0, "repair parity");
                }
            }
    Check(row_count == 144, "complete trace");
}
}

int main(int argc, char** argv)
{
    try {
        Check(argc == 4 && (std::strcmp(argv[3],"normal") == 0 || std::strcmp(argv[3],"reverse") == 0), "arguments");
        const bool reverse = std::strcmp(argv[3],"reverse") == 0;
        Library first(argv[reverse ? 2 : 1]), second(argv[reverse ? 1 : 2]);
        Check(first.lib != second.lib, "distinct DSO files");
        Trace(reverse ? &second : &first, reverse ? &first : &second);
        std::puts("k,width,policy,creation_order,arm,allocation,bytes,array,address,mod64,mod4096");
        for (unsigned i = 0; i < row_count; ++i) {
            const Row& r = rows[i];
            std::printf("%u,%u,%u,%u,%u,%u,%zu,%u,%llu,%llu,%llu\n",r.k,r.width,r.policy,r.creation_order,
                r.arm,r.allocation,r.bytes,unsigned(r.array),static_cast<unsigned long long>(r.pointer),
                static_cast<unsigned long long>(r.pointer%64),static_cast<unsigned long long>(r.pointer%4096));
        }
        Check(!std::ferror(stdout) && std::fflush(stdout) == 0, "output");
    } catch (const std::exception& e) {
        tracking = false;
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        return 1;
    }
}
