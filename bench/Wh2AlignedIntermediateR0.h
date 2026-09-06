#pragma once

// Benchmark-only retained-storage allocator.  This is not a production option.
#include "../WirehairTools.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <type_traits>
#include <vector>

#if defined(WH_ALIGN64) || defined(WH_HUGEPAGE)
#error "The aligned-intermediate prototype forbids allocation overrides"
#endif

static_assert(GF256_ALIGN_BYTES == 32,
    "The aligned-intermediate prototype requires the native 32-byte helper");

namespace wh2_aligned_intermediate_r0 {

template<class T>
class Allocator
{
public:
    static_assert(alignof(T) <= 32, "Unsupported over-aligned allocation");
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::true_type propagate_on_container_move_assignment;
    typedef std::true_type propagate_on_container_swap;
    typedef std::true_type is_always_equal;

    template<class U> struct rebind { typedef Allocator<U> other; };

    Allocator() noexcept = default;
    template<class U> Allocator(const Allocator<U>&) noexcept {}

    size_type max_size() const noexcept
    {
        const size_type byte_limit =
            (std::numeric_limits<size_type>::max() - 32u) / sizeof(T);
        const size_type difference_limit = static_cast<size_type>(
            std::numeric_limits<difference_type>::max()) / sizeof(T);
        return byte_limit < difference_limit ? byte_limit : difference_limit;
    }

    pointer allocate(size_type count, const void* = nullptr)
    {
        // Zero allocations have no ownership and do not call the helper.
        if (count == 0u) {
            return nullptr;
        }
        if (count > max_size()) {
            throw std::bad_alloc();
        }
        // max_size checks both multiplication and the helper's 32-byte pad.
        void* const allocation = wirehair::SIMDSafeAllocate(count * sizeof(T));
        if (!allocation) {
            throw std::bad_alloc();
        }
        return static_cast<pointer>(allocation);
    }

    void deallocate(pointer allocation, size_type) noexcept
    {
        wirehair::SIMDSafeFree(allocation);
    }
};

template<>
class Allocator<void>
{
public:
    typedef void value_type;
    typedef void* pointer;
    typedef const void* const_pointer;
    template<class U> struct rebind { typedef Allocator<U> other; };
    Allocator() noexcept = default;
    template<class U> Allocator(const Allocator<U>&) noexcept {}
};

template<class T, class U>
bool operator==(const Allocator<T>&, const Allocator<U>&) noexcept
{
    return true;
}

template<class T, class U>
bool operator!=(const Allocator<T>&, const Allocator<U>&) noexcept
{
    return false;
}

typedef std::vector<std::uint8_t, Allocator<std::uint8_t> > ByteVector;

} // namespace wh2_aligned_intermediate_r0
