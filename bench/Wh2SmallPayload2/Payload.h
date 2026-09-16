#ifndef WH2_SMALL_PAYLOAD2_EXPERIMENT_H
#define WH2_SMALL_PAYLOAD2_EXPERIMENT_H

#include "../../codec/WirehairK6Payload.h"
#include <cstring>

namespace wh2_small_payload2 {
// Benchmark-only candidate. Same disjoint-output/readable-input contract as
// Payload; no equation, allocation, field, lookup or wider-kernel change.
template<unsigned K>
inline void Payload(void* output, const void* const* sources,
                    const std::uint8_t* scales, int bytes)
{
#ifndef WH_COUNT
    if ((K == 3 || K == 5 || K == 8) && bytes == 2) {
        std::uint8_t result[2] = {};
        for (unsigned byte = 0; byte < 2; ++byte)
            for (unsigned source = 0; source < K; ++source)
                result[byte] ^= gf256_mul(
                    static_cast<const std::uint8_t*>(sources[source])[byte], scales[source]);
        std::memcpy(output, result, sizeof(result));
        return;
    }
#endif
    wirehair_k6_payload::Payload(output, sources, scales, K, bytes);
}
} // namespace wh2_small_payload2
#endif
