#ifndef WH2_SMALL_PAYLOAD2_OUTLINED_H
#define WH2_SMALL_PAYLOAD2_OUTLINED_H

#include "../../codec/WirehairK6Payload.h"

namespace wh2_small_payload2_outlined {
// Private benchmark-only dispatcher. Keeping the existing five-argument call
// boundary lets the compiler retain the surrounding encoder's code shape.
// Same initialized-runtime, readable-input, disjoint-output contract as Payload.
void Apply(void* output, const void* const* sources, const std::uint8_t* scales,
           int count, int bytes);

template<unsigned K>
inline void Payload(void* output, const void* const* sources,
                    const std::uint8_t* scales, int bytes)
{
    Apply(output, sources, scales, K, bytes);
}
} // namespace wh2_small_payload2_outlined
#endif
