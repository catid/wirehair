#ifndef WH2_THUE_MORSE_TINY_PAYLOAD_COST_R0_BRIDGE_H
#define WH2_THUE_MORSE_TINY_PAYLOAD_COST_R0_BRIDGE_H

#include <cstddef>
#include <cstdint>
#include <type_traits>

// Benchmark-only POD boundary. Handles belong to exactly the Api that created
// them; callers must never mix the two codec namespaces or free functions.
namespace wh2_tiny_payload_cost_r0 {

enum Status {
    Success = 0, NeedMore = 1, InvalidInput = 2, BufferTooSmall = 3,
    Conflict = 4, OutOfMemory = 5
};

struct Result {
    int status;
    std::size_t required;
    std::size_t written;
};

struct Roundtrip {
    int create_status;
    int feed_status[6];
    unsigned feed_count;
    int recover_status;
    std::size_t required;
    std::size_t written;
    bool recovered;
};

struct Api {
    // Empty output handle required. Failures preserve it. Source and lookup
    // remain borrowed with exactly the original benchmark codec's lifetimes.
    int (*create)(const std::uint8_t*, std::size_t, const void*, std::uint64_t,
                  std::uint32_t, void**);
    Result (*encode)(void*, std::uint32_t, void*, std::size_t);
    void (*encoder_free)(void*);
    int (*row)(const std::uint8_t*, std::size_t, std::uint32_t, std::uint8_t*);

    // Untimed preflight only: six contiguous full-B repairs, IDs6..11, fed
    // until first Success. Creates/frees a real decoder on every call. Output
    // capacity must cover M, and its M-byte range must not overlap source,
    // lookup or the repair corpus. Caller owns all buffers and their guards.
    // Unobserved statuses are -1. Invalid bridge arguments fail before create.
    Roundtrip (*roundtrip)(const std::uint8_t*, std::size_t, const void*,
                          std::uint64_t, std::uint32_t, const void*,
                          std::size_t, void*, std::size_t);
};

static_assert(std::is_standard_layout<Result>::value && std::is_trivial<Result>::value,
              "Result must remain POD");
static_assert(std::is_standard_layout<Roundtrip>::value && std::is_trivial<Roundtrip>::value,
              "Roundtrip must remain POD");
static_assert(std::is_standard_layout<Api>::value && std::is_trivial<Api>::value,
              "Api must remain POD");

} // namespace wh2_tiny_payload_cost_r0

extern "C" const wh2_tiny_payload_cost_r0::Api* wh2_tiny_payload_cost_baseline_api();
extern "C" const wh2_tiny_payload_cost_r0::Api* wh2_tiny_payload_cost_candidate_api();

#endif
