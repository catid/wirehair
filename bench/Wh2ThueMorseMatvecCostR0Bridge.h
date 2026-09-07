#ifndef WH2_THUE_MORSE_MATVEC_COST_R0_BRIDGE_H
#define WH2_THUE_MORSE_MATVEC_COST_R0_BRIDGE_H

#include <cstddef>
#include <cstdint>
#include <type_traits>

// Benchmark-only boundary: no codec class or owning C++ object crosses it.
namespace wh2_matvec_cost_r0 {

enum Status { Success = 0, NeedMore = 1, Failure = 2 };
enum LengthKind { NoLength = 0, ReturnedLength = 1, InferredLength = 2 };

struct Result {
    int status;                 // Shared classification only.
    int code;                   // Unchanged arm-local code; -1 = bridge rejection.
    std::uint64_t required;      // Exposed/inferred size; WH1 short Encode reports 0.
    std::uint64_t written;       // Zero unless the operation succeeds.
    unsigned length_kind;       // WH1 Recover has no returned byte count.
};

// Zero-initialize before first use. A live state is noncopyable by contract;
// its handle belongs to exactly the Api that created it. Do not mutate its
// metadata while live. No wrapper allocation.
// Source and native lookup remain immutable/alive until the relevant free.
// The emitted public profile belongs to this POD, not to encoder-private memory.
struct State {
    void* handle;
    const void* source;
    const std::uint8_t* lookup;
    std::size_t lookup_bytes;
    std::uint64_t message_bytes;
    std::uint32_t block_bytes;
    unsigned arm;
    unsigned kind;              // 0 = empty, 1 = encoder, 2 = decoder.
    std::uint32_t profile_bytes;
    std::uint8_t profile[32];
};

struct Api {
    Result (*create)(const std::uint8_t*, std::size_t, const void*,
                     std::uint64_t, std::uint32_t, State*);
    Result (*encode)(const State*, std::uint32_t, void*, std::size_t);
    void (*encoder_free)(State*);
    Result (*decoder_create)(const State* encoder, State* decoder);
    Result (*feed)(State*, std::uint32_t, const void*, std::size_t);
    Result (*recover)(State*, void*, std::size_t);
    void (*decoder_free)(State*);

    // Untimed metadata validation only. Row is null for WH1/public WH2.
    bool (*valid_profile)(const State*);
    Result (*row)(const std::uint8_t*, std::size_t, std::uint32_t, std::uint8_t*);
};

static_assert(std::is_standard_layout<Result>::value && std::is_trivial<Result>::value,
              "Result must be POD");
static_assert(std::is_standard_layout<State>::value && std::is_trivial<State>::value,
              "State must be POD");
static_assert(std::is_standard_layout<Api>::value && std::is_trivial<Api>::value,
              "Api must be POD");

} // namespace wh2_matvec_cost_r0

extern "C" const wh2_matvec_cost_r0::Api* wh2_matvec_cost_baseline_api();
extern "C" const wh2_matvec_cost_r0::Api* wh2_matvec_cost_candidate_api();
extern "C" const wh2_matvec_cost_r0::Api* wh2_matvec_cost_wh1_api();
extern "C" const wh2_matvec_cost_r0::Api* wh2_matvec_cost_public_api();

#endif
