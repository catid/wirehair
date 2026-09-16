#include <wirehair/wirehair.h>
#include "WirehairSmallCore.h"
#include "WirehairSmallLookup.h"
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <vector>

namespace {
struct Allocation { void* pointer; size_t bytes; bool array; bool live; };
Allocation records[16];
bool tracking = false;
size_t attempts = 0, used = 0, freed = 0, failure = SIZE_MAX;
void Check(bool value, const char* what)
{
    if (!value) { std::fprintf(stderr, "FAIL: %s\n", what); std::abort(); }
}
void* Allocate(size_t bytes, bool array)
{
    if (tracking && attempts++ == failure) throw std::bad_alloc();
    void* pointer = std::malloc(bytes ? bytes : 1);
    if (!pointer) throw std::bad_alloc();
    if (tracking) {
        Check(used < 16, "bounded allocation ledger");
        records[used++] = {pointer, bytes, array, true};
        std::memset(pointer, 0xa5, bytes);
    }
    return pointer;
}
void Release(void* pointer, bool array, size_t bytes = SIZE_MAX) noexcept
{
    if (pointer && tracking) {
        size_t slot = 0;
        while (slot < used && records[slot].pointer != pointer) ++slot;
        Check(slot < used && records[slot].live, "free exact live allocation base once");
        Check(records[slot].array == array, "scalar/array allocation and delete agree");
        Check(bytes == SIZE_MAX || bytes == records[slot].bytes, "sized delete matches full allocation");
        records[slot].live = false;
        ++freed;
    }
    std::free(pointer);
}
void Begin(size_t fail = SIZE_MAX)
{
    Check(!tracking, "non-nested allocation recording");
    attempts = used = freed = 0;
    failure = fail;
    tracking = true;
}
void End(size_t expected)
{
    tracking = false;
    Check(attempts == expected && freed == used, "complete allocation/failure/free ledger");
}
}

#define NOINLINE __attribute__((noinline))
NOINLINE void* operator new(size_t n) { return Allocate(n, false); }
NOINLINE void* operator new[](size_t n) { return Allocate(n, true); }
NOINLINE void operator delete(void* p) noexcept { Release(p, false); }
NOINLINE void operator delete[](void* p) noexcept { Release(p, true); }
NOINLINE void* operator new(size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n, false); } catch (const std::bad_alloc&) { return nullptr; } }
NOINLINE void* operator new[](size_t n, const std::nothrow_t&) noexcept
{ try { return Allocate(n, true); } catch (const std::bad_alloc&) { return nullptr; } }
NOINLINE void operator delete(void* p, const std::nothrow_t&) noexcept { Release(p, false); }
NOINLINE void operator delete[](void* p, const std::nothrow_t&) noexcept { Release(p, true); }
#if defined(__cpp_sized_deallocation)
NOINLINE void operator delete(void* p, size_t n) noexcept { Release(p, false, n); }
NOINLINE void operator delete[](void* p, size_t n) noexcept { Release(p, true, n); }
#endif

namespace {
namespace S = wirehair_small_core;
#ifdef STORAGE_CANDIDATE
constexpr size_t CoreAllocations = 1;
#else
constexpr size_t CoreAllocations = 2;
#endif
using Byte = uint8_t;
using Descriptor = std::array<Byte, 32>;
size_t core_cases = 0, public_cases = 0, oom_cases = 0;

// Preserve the actual allocating factory boundary. Inlined new-expressions
// may be elided by an optimizing compiler even with replacement allocation
// functions; public-library constructors already have a separate TU boundary.
template<unsigned K> NOINLINE S::Status MakeCore(S::Lookup lookup, uint64_t message,
    uint32_t block, std::unique_ptr<S::Decoder<K>>& output)
{ return S::Decoder<K>::Create(lookup, message, block, output); }

template<unsigned K> void Direct(S::Lookup lookup, uint32_t block, uint32_t tail)
{
    using D = S::Decoder<K>;
    static_assert(alignof(D) <= alignof(std::max_align_t), "Raw allocation provides decoder alignment");
    const size_t message = size_t(K - 1) * block + tail;
    const size_t slab_bytes = size_t(K + 1) * block;
    std::vector<Byte> source(message), output(message + 2, 0xcc), bad(block, 0);
    for (size_t i = 0; i < message; ++i) source[i] = Byte((i * 29 + 17) & 255);
    std::unique_ptr<D> decoder;
    Begin();
    Check(MakeCore<K>(lookup, message, block, decoder) == S::Status::Success, "core create");
    Check(used == CoreAllocations && records[0].pointer == decoder.get() &&
          reinterpret_cast<uintptr_t>(decoder.get()) % alignof(D) == 0, "core allocation count/base/alignment");
    Byte* slab = nullptr;
#ifdef STORAGE_CANDIDATE
    Check(records[0].bytes == sizeof(D) + slab_bytes && !records[0].array, "one exact combined allocation");
    slab = static_cast<Byte*>(records[0].pointer) + sizeof(D);
#else
    Check(records[0].bytes == sizeof(D) && !records[0].array && records[1].bytes == slab_bytes && records[1].array,
          "original object and slab allocations");
    slab = static_cast<Byte*>(records[1].pointer);
#endif
    Check(decoder->Feed(0, decoder.get(), block).status == S::Status::InvalidInput && decoder->Rank() == 0,
          "reject object input alias before access");
    Check(decoder->Feed(0, slab, block).status == S::Status::InvalidInput && decoder->Rank() == 0,
          "reject trailing slab input alias before access");
#ifdef STORAGE_CANDIDATE
    Check(decoder->Feed(0, slab - 1, block).status == S::Status::InvalidInput && decoder->Rank() == 0,
          "reject input spanning adjacent object/slab boundary");
    Check(decoder->Recover(slab - 1, message).status == S::Status::InvalidInput,
          "reject output spanning adjacent object/slab boundary");
#endif
    Check(decoder->Recover(output.data() + 1, message).status == S::Status::NeedMore,
          "incomplete recovery");
    for (Byte b : output) Check(b == 0xcc, "incomplete output unchanged");
    // Reverse pivot order and partial systematic tail must use the same slab.
    for (unsigned i = K; i-- > 0;) {
        const size_t bytes = i == K - 1 ? tail : block;
        const auto result = decoder->Feed(i, source.data() + size_t(i) * block, bytes);
        Check(result.status == (i == 0 ? S::Status::Success : S::Status::NeedMore), "systematic feed/rank");
    }
    Check(decoder->Recover(output.data() + 1, message - 1).status == S::Status::BufferTooSmall,
          "undersized recovery");
    for (Byte b : output) Check(b == 0xcc, "undersized output unchanged");
    Check(decoder->Recover(slab, message).status == S::Status::InvalidInput, "reject output slab alias");
    for (unsigned repeat = 0; repeat < 2; ++repeat) {
        std::memset(output.data(), 0xcc, output.size());
        Check(decoder->Recover(output.data() + 1, message).status == S::Status::Success, "repeat recovery");
        Check(output.front() == 0xcc && output.back() == 0xcc &&
              std::memcmp(source.data(), output.data() + 1, message) == 0, "exact guarded recovered bytes");
    }
    std::memcpy(bad.data(), source.data(), block);
    bad[0] ^= 1;
    Check(decoder->Feed(0, bad.data(), block).status == S::Status::Conflict, "contradiction preserved");
    Check(decoder->Feed(0, source.data(), block).status == S::Status::Success, "valid duplicate preserved");
    std::memset(output.data(), 0xcc, output.size());
    Check(decoder->Recover(output.data() + 1, message).status == S::Status::Success &&
          output.front() == 0xcc && output.back() == 0xcc &&
          std::memcmp(source.data(), output.data() + 1, message) == 0, "conflict preserves solved basis");
    Check(attempts == CoreAllocations, "no feed or recovery allocation");
    decoder.reset();
    End(CoreAllocations);
    ++core_cases;
    for (size_t fail = 0; fail < CoreAllocations; ++fail) {
        Begin(fail);
        Check(MakeCore<K>(lookup, message, block, decoder) == S::Status::OutOfMemory && !decoder,
              "core OOM preserves output");
        End(fail + 1);
        ++oom_cases;
    }
}

void Public(unsigned k, uint64_t profile_id, uint32_t block, uint32_t tail)
{
    const size_t message = size_t(k - 1) * block + tail;
    std::vector<Byte> source(message), output(message + 2, 0xcc), packets(size_t(k) * block);
    for (size_t i = 0; i < message; ++i) source[i] = Byte((i * 31 + 11) & 255);
    for (bool borrowed : {false, true}) for (bool distant : {false, true}) {
        Descriptor profile = {};
        WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
        options.source_policy = borrowed ? WirehairV2EncoderSource_BorrowedImmutable : WirehairV2EncoderSource_Independent;
        WirehairV2Codec encoder = nullptr;
        uint32_t count = 0;
        Check(wirehair_v2_encoder_create_profile_id_with_options(profile_id, source.data(), message, block,
            &options, profile.data(), 32, &count, &encoder) == WirehairV2_Success && count == 32, "public encoder fixture");
        for (unsigned i = 0; i < k; ++i) {
            const uint32_t id = distant ? UINT32_MAX - 2 * i : k + i;
            Check(wirehair_v2_encode(encoder, id, packets.data() + size_t(i) * block, block, &count) ==
                  WirehairV2_Success && count == block, "public repair fixture");
        }
        wirehair_v2_free(encoder);
        WirehairV2Codec decoder = nullptr;
        Begin();
        Check(wirehair_v2_decoder_create(profile.data(), 32, &decoder) == WirehairV2_Success &&
              used == CoreAllocations + 1, "public decoder allocates facade plus core storage");
        WirehairV2Result status = WirehairV2_NeedMore;
        for (unsigned i = 0; i < k && status == WirehairV2_NeedMore; ++i) {
            const uint32_t id = distant ? UINT32_MAX - 2 * i : k + i;
            status = wirehair_v2_decode(decoder, id, packets.data() + size_t(i) * block, block);
        }
        // Never assume arbitrary distant rows are independent: preserve fallback.
        for (unsigned i = 0; i < k && status == WirehairV2_NeedMore; ++i)
            status = wirehair_v2_decode(decoder, i, source.data() + size_t(i) * block, i == k - 1 ? tail : block);
        Check(status == WirehairV2_Success, "public decode own first success");
        for (unsigned repeat = 0; repeat < 2; ++repeat) {
            std::memset(output.data(), 0xcc, output.size());
            uint64_t recovered = 0;
            Check(wirehair_v2_recover(decoder, output.data() + 1, message, &recovered) == WirehairV2_Success &&
                  recovered == message && output.front() == 0xcc && output.back() == 0xcc &&
                  std::memcmp(output.data() + 1, source.data(), message) == 0, "public guarded repeated recovery");
        }
        wirehair_v2_free(decoder);
        End(CoreAllocations + 1);
        ++public_cases;
        for (size_t fail = 0; fail < CoreAllocations + 1; ++fail) {
            decoder = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
            Begin(fail);
            Check(wirehair_v2_decoder_create(profile.data(), 32, &decoder) == WirehairV2_OOM && !decoder,
                  "public constructor OOM unwinds");
            End(fail + 1);
            ++oom_cases;
        }
    }
}
}

int main()
{
    Check(gf256_init() == 0 && wirehair_init() == Wirehair_Success, "GF and library initialization");
    for (uint32_t block : {1u, 2u, 15u, 16u, 17u, 31u, 32u, 33u, 63u, 64u, 65u, 127u,
                           128u, 129u, 255u, 256u, 257u, 1280u, 4096u}) {
        for (uint32_t tail : {1u, block}) {
            Direct<3>(S::K3Lookup(), block, tail);
            Direct<5>(S::K5Lookup(), block, tail);
            Direct<8>(S::K8Lookup(), block, tail);
            Public(3, WIREHAIR_V2_PROFILE_SMALL_K3_2026_09, block, tail);
            Public(5, WIREHAIR_V2_PROFILE_SMALL_K5_2026_09, block, tail);
            Public(8, WIREHAIR_V2_PROFILE_SMALL_K8_2026_09, block, tail);
        }
    }
    std::printf("core_cases=%zu public_cases=%zu oom_cases=%zu core_allocations=%zu\n",
                core_cases, public_cases, oom_cases, CoreAllocations);
    return 0;
}
