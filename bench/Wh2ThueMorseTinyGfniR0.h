#ifndef WH2_THUE_MORSE_TINY_GFNI_R0_H
#define WH2_THUE_MORSE_TINY_GFNI_R0_H

#include "gf256.h"
#include <cstdint>

#ifdef WH_COUNT
#error "This benchmark-only payload prototype does not implement GF work counters"
#endif

namespace wh2_thue_tiny_gfni_r0 {

// An isomorphism between two byte bases of the SAME GF(256), not a new field.
// Map x in 0x14d to the smallest root (0x46) of that polynomial in 0x11b.
// GFNI orders output bit i's row in matrix byte 7-i. The neutral test derives
// both constants independently and verifies every possible input product.
static const std::uint64_t kToHardware = UINT64_C(0x0dae86dcbcfca25c);
static const std::uint64_t kFromHardware = UINT64_C(0x43e06c2ef2283088);

// Requires the initialized shared GF runtime. Exactly six readable pointers,
// six readable scale bytes and two readable bytes per source. Readable source
// spans may overlap each other; the two-byte output is disjoint from all input.
bool Available();
void Scalar(void* output, const void* const* sources, const std::uint8_t* scales);
void Apply(void* output, const void* const* sources, const std::uint8_t* scales);

// Preserve the generic function's contract, including nonpositive no-ops.
inline void Payload(void* output, const void* const* sources,
                    const std::uint8_t* scales, int count, int bytes)
{
    if (count == 6 && bytes == 2) {
        Apply(output, sources, scales);
        return;
    }
    gf256_mulset_multi_mem(output, sources, scales, count, bytes);
}

} // namespace wh2_thue_tiny_gfni_r0
#endif
