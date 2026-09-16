// Non-timing C++ allocation observer. Ordinary malloc/free; no forced placement.
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>

namespace allocation_trace {
struct Event {
    std::uintptr_t pointer, caller, stack;
    std::size_t bytes;
    bool allocate, array;
    unsigned phase;
};
struct Invocation {
    unsigned phase;
    std::uintptr_t handle, first, second;
    std::size_t bytes;
    std::uint32_t packet;
};
const unsigned Capacity = 256; // 32 lifecycles x (3 allocations + 3 frees).
Event events[Capacity];
unsigned count = 0;
Invocation calls[1152]; // 32 x (constructor + up to 32 feeds + recover + free).
unsigned call_count = 0, phase = 4;
bool enabled = false, overflow = false;

void Call(unsigned stage, const void* handle, const void* first, const void* second,
          std::size_t bytes, std::uint32_t packet = 0) noexcept
{
    if (!enabled) return;
    phase = stage;
    if (call_count == 1152) { overflow = true; return; }
    calls[call_count++] = Invocation{stage,reinterpret_cast<std::uintptr_t>(handle),
        reinterpret_cast<std::uintptr_t>(first),reinterpret_cast<std::uintptr_t>(second),bytes,packet};
}

void Record(void* p, std::size_t bytes, bool allocate, bool array,
            void* caller, const void* stack) noexcept
{
    if (!enabled || !p) return;
    if (count == Capacity) { overflow = true; return; }
    events[count++] = Event{reinterpret_cast<std::uintptr_t>(p),
        reinterpret_cast<std::uintptr_t>(caller), reinterpret_cast<std::uintptr_t>(stack),
        bytes, allocate, array, phase};
}
__attribute__((noinline)) void* Allocate(std::size_t bytes, bool array, void* caller)
{
    void* p = std::malloc(bytes ? bytes : 1);
    if (!p) throw std::bad_alloc();
    Record(p, bytes, true, array, caller, &p);
    return p;
}
__attribute__((noinline)) void Release(void* p, bool array, void* caller) noexcept
{
    Record(p, 0, false, array, caller, &p);
    std::free(p);
}
void Start() noexcept { count = 0; call_count = 0; phase = 4; overflow = false; enabled = true; }
void Stop() noexcept { enabled = false; }
}

// Capture the codec caller before entering the common observer helper. These
// markers describe this interposed observer only, not old timed stack addresses.
#define CALLER __builtin_return_address(0)
__attribute__((noinline)) void* operator new(std::size_t n)
{ return allocation_trace::Allocate(n, false, CALLER); }
__attribute__((noinline)) void* operator new[](std::size_t n)
{ return allocation_trace::Allocate(n, true, CALLER); }
__attribute__((noinline)) void operator delete(void* p) noexcept
{ allocation_trace::Release(p, false, CALLER); }
__attribute__((noinline)) void operator delete[](void* p) noexcept
{ allocation_trace::Release(p, true, CALLER); }
__attribute__((noinline)) void* operator new(std::size_t n, const std::nothrow_t&) noexcept
{ try { return allocation_trace::Allocate(n, false, CALLER); } catch (...) { return nullptr; } }
__attribute__((noinline)) void* operator new[](std::size_t n, const std::nothrow_t&) noexcept
{ try { return allocation_trace::Allocate(n, true, CALLER); } catch (...) { return nullptr; } }
__attribute__((noinline)) void operator delete(void* p, const std::nothrow_t&) noexcept
{ allocation_trace::Release(p, false, CALLER); }
__attribute__((noinline)) void operator delete[](void* p, const std::nothrow_t&) noexcept
{ allocation_trace::Release(p, true, CALLER); }
#if defined(__cpp_sized_deallocation)
__attribute__((noinline)) void operator delete(void* p, std::size_t) noexcept
{ allocation_trace::Release(p, false, CALLER); }
__attribute__((noinline)) void operator delete[](void* p, std::size_t) noexcept
{ allocation_trace::Release(p, true, CALLER); }
#endif
#undef CALLER
