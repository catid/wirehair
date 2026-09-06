#ifndef WH2_THUE_MORSE_TINY_PAYLOAD_R0_H
#define WH2_THUE_MORSE_TINY_PAYLOAD_R0_H

#include "gf256.h"

#ifdef WH_COUNT
#error "This benchmark-only payload prototype does not implement GF work counters"
#endif

namespace wh2_thue_tiny_payload_r0 {

// Same initialized-runtime, readable-source and disjoint-output contract as
// gf256_mulset_multi_mem. No new storage, ISA requirement or packet equations.
inline void Payload(void* destination, const void* const* sources,
                    const uint8_t* scales, int source_count, int bytes)
{
    if (source_count != 6 || bytes != 2) {
        gf256_mulset_multi_mem(destination, sources, scales, source_count, bytes);
        return;
    }
    uint8_t first = 0, second = 0;
    for (unsigned source = 0; source < 6; ++source) {
        const uint8_t* input = static_cast<const uint8_t*>(sources[source]);
        first ^= gf256_mul(input[0], scales[source]);
        second ^= gf256_mul(input[1], scales[source]);
    }
    uint8_t* output = static_cast<uint8_t*>(destination);
    output[0] = first;
    output[1] = second;
}

} // namespace wh2_thue_tiny_payload_r0
#endif
