#include "Payload.h"
#include <cstring>

namespace wh2_small_payload2_outlined {
void Apply(void* output, const void* const* sources, const std::uint8_t* scales,
           int count, int bytes)
{
#ifndef WH_COUNT
    if (bytes == 2 && (count == 3 || count == 5 || count == 8)) {
        std::uint8_t result[2] = {};
        for (int source = 0; source < count; ++source) {
            const auto* input = static_cast<const std::uint8_t*>(sources[source]);
            // Source-major: reuse the coefficient/table base for both bytes.
            result[0] ^= gf256_mul(input[0], scales[source]);
            result[1] ^= gf256_mul(input[1], scales[source]);
        }
        std::memcpy(output, result, sizeof(result));
        return;
    }
#endif
    // Preserve nonpositive no-ops, K6 dispatch, wider kernels and op counts.
    wirehair_k6_payload::Payload(output, sources, scales, count, bytes);
}
} // namespace wh2_small_payload2_outlined
