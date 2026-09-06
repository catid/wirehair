#include "Wh2ThueMorseMulRowsR0.h"
#include "gf256.h"

#include <cstring>
#include <limits>

#if defined(_MSC_VER)
#define WH2_MULROWS_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_MULROWS_NOINLINE __attribute__((noinline))
#else
#define WH2_MULROWS_NOINLINE
#endif

namespace wh2_thue_mulrows_r0 {
using wh2_thue_native::Lookup;
using wh2_thue_native::Status;
using wh2_thue_native::kLookupBytes;
using wh2_thue_native::kSourceCount;
namespace {

bool Span(const void* pointer, std::size_t bytes)
{
    return pointer && bytes <= std::numeric_limits<std::uintptr_t>::max() -
        reinterpret_cast<std::uintptr_t>(pointer);
}

// Both spans have already passed Span. Only potentially written bytes need
// alias exclusion; capacity is nevertheless checked separately for wrapping.
bool Overlap(const void* a, std::size_t a_bytes, const void* b, std::size_t b_bytes)
{
    const std::uintptr_t x = reinterpret_cast<std::uintptr_t>(a);
    const std::uintptr_t y = reinterpret_cast<std::uintptr_t>(b);
    return a_bytes && b_bytes && x < y + b_bytes && y < x + a_bytes;
}

bool ValidLookup(Lookup lookup)
{
    return lookup.bytes == kLookupBytes && Span(lookup.data, lookup.bytes);
}

unsigned Parity(unsigned value)
{
    value ^= value >> 4;
    value ^= value >> 2;
    value ^= value >> 1;
    return value & 1;
}

void Apply(const std::uint8_t* matrix, std::uint8_t vector[kSourceCount])
{
    std::uint8_t result[kSourceCount];
    for (unsigned r = 0; r < kSourceCount; ++r) {
        std::uint8_t value = gf256_mul(matrix[r * kSourceCount], vector[0]);
        for (unsigned c = 1; c < kSourceCount; ++c)
            value ^= gf256_mul(matrix[r * kSourceCount + c], vector[c]);
        result[r] = value;
    }
    std::memcpy(vector, result, kSourceCount);
}

// Preserve the original core's out-of-line Row -> Map call boundary.
WH2_MULROWS_NOINLINE void Map(Lookup lookup, std::uint32_t id, std::uint8_t output[kSourceCount])
{
    if (id < 1024) {
        std::memcpy(output, lookup.data + id * kSourceCount, kSourceCount);
        return;
    }
    const unsigned high = id >> 24;
    const unsigned mid17 = (id >> 17) & 127;
    const unsigned mid10 = (id >> 10) & 127;
    const unsigned low = id & 1023;
    const unsigned phase17 = Parity(high);
    const unsigned phase10 = phase17 ^ Parity(mid17);
    const unsigned phase0 = phase10 ^ Parity(mid10);
    std::memcpy(output, lookup.data + phase0 * 6144 + low * kSourceCount, kSourceCount);
    if (mid10) Apply(lookup.data + 12288 + phase10 * 4608 + mid10 * 36, output);
    if (mid17) Apply(lookup.data + 21504 + phase17 * 4608 + mid17 * 36, output);
    if (high) Apply(lookup.data + 30720 + high * 36, output);
}

} // namespace

Status Row(Lookup lookup, std::uint32_t id, std::uint8_t output[kSourceCount])
{
    if (!ValidLookup(lookup) || !Span(output, kSourceCount) ||
        Overlap(output, kSourceCount, lookup.data, lookup.bytes)) return Status::InvalidInput;
    Map(lookup, id, output);
    return Status::Success;
}

} // namespace wh2_thue_mulrows_r0

#undef WH2_MULROWS_NOINLINE
