#include "Wh2ThueMorseTinyPayloadCostR0Bridge.h"

#if !defined(WH2_TINY_COST_CANDIDATE) || \
    (WH2_TINY_COST_CANDIDATE != 0 && WH2_TINY_COST_CANDIDATE != 1)
#error "Compile this bridge with WH2_TINY_COST_CANDIDATE=0 or 1"
#endif
#ifdef WH_COUNT
#error "Tiny payload cost wrappers must not omit GF work counters"
#endif

#if WH2_TINY_COST_CANDIDATE
#define wh2_thue_native wh2_thue_tiny_payload_r0
#include "Wh2ThueMorseNativeCodec.h"
#undef wh2_thue_native
namespace N = wh2_thue_tiny_payload_r0;
#define WH2_TINY_BRIDGE_NAMESPACE wh2_tiny_payload_cost_candidate_bridge
#define WH2_TINY_BRIDGE_GETTER wh2_tiny_payload_cost_candidate_api
#else
#include "Wh2ThueMorseNativeCodec.h"
namespace N = wh2_thue_native;
#define WH2_TINY_BRIDGE_NAMESPACE wh2_tiny_payload_cost_baseline_bridge
#define WH2_TINY_BRIDGE_GETTER wh2_tiny_payload_cost_baseline_api
#endif

#include <cstring>
#include <limits>
#include <memory>

namespace WH2_TINY_BRIDGE_NAMESPACE {
namespace B = wh2_tiny_payload_cost_r0;

static_assert(sizeof(N::Encoder) == 96, "Frozen native handle-span receipt ABI");

static_assert(static_cast<int>(N::Status::Success) == B::Success &&
              static_cast<int>(N::Status::NeedMore) == B::NeedMore &&
              static_cast<int>(N::Status::InvalidInput) == B::InvalidInput &&
              static_cast<int>(N::Status::BufferTooSmall) == B::BufferTooSmall &&
              static_cast<int>(N::Status::Conflict) == B::Conflict &&
              static_cast<int>(N::Status::OutOfMemory) == B::OutOfMemory,
              "Bridge status mapping must match the original codec");

namespace {
bool Span(const void* pointer, std::size_t bytes)
{
    return pointer && bytes <= std::numeric_limits<std::uintptr_t>::max() -
        reinterpret_cast<std::uintptr_t>(pointer);
}

bool Overlap(const void* first, std::size_t first_bytes,
             const void* second, std::size_t second_bytes)
{
    const std::uintptr_t a = reinterpret_cast<std::uintptr_t>(first);
    const std::uintptr_t b = reinterpret_cast<std::uintptr_t>(second);
    return first_bytes && second_bytes && a < b + second_bytes && b < a + first_bytes;
}
} // namespace

int Create(const std::uint8_t* lookup, std::size_t lookup_bytes,
           const void* source, std::uint64_t message_bytes,
           std::uint32_t block_bytes, void** output)
{
    if (message_bytes > std::numeric_limits<std::size_t>::max() ||
        !Span(output, sizeof(*output)) ||
        reinterpret_cast<std::uintptr_t>(output) % alignof(void*) != 0 ||
        !Span(source, static_cast<std::size_t>(message_bytes)) ||
        !Span(lookup, lookup_bytes) ||
        Overlap(output, sizeof(*output), source, static_cast<std::size_t>(message_bytes)) ||
        Overlap(output, sizeof(*output), lookup, lookup_bytes)) return B::InvalidInput;
    if (*output) return B::InvalidInput;
    std::unique_ptr<N::Encoder> encoder;
    const N::Status status = N::Encoder::Create(
        N::Lookup{lookup, lookup_bytes}, source, message_bytes, block_bytes, encoder);
    if (status == N::Status::Success) *output = encoder.release();
    return static_cast<int>(status);
}

B::Result Encode(void* handle, std::uint32_t id, void* output, std::size_t capacity)
{
    if (!handle) return B::Result{B::InvalidInput, 0, 0};
    const N::Result result = static_cast<N::Encoder*>(handle)->Encode(id, output, capacity);
    return B::Result{static_cast<int>(result.status), result.bytes_required, result.bytes_written};
}

void EncoderFree(void* handle)
{
    delete static_cast<N::Encoder*>(handle);
}

int Row(const std::uint8_t* lookup, std::size_t lookup_bytes,
        std::uint32_t id, std::uint8_t* output)
{
    return static_cast<int>(N::Row(N::Lookup{lookup, lookup_bytes}, id, output));
}

B::Roundtrip Roundtrip(const std::uint8_t* lookup, std::size_t lookup_bytes,
                       const void* source, std::uint64_t message_bytes,
                       std::uint32_t block_bytes, const void* repair_packets,
                       std::size_t packet_bytes, void* output, std::size_t capacity)
{
    B::Roundtrip result = {};
    result.create_status = B::InvalidInput;
    result.recover_status = -1;
    for (unsigned i = 0; i < 6; ++i) result.feed_status[i] = -1;
    if (message_bytes > std::numeric_limits<std::size_t>::max() ||
        static_cast<std::uint64_t>(block_bytes) * 6 != packet_bytes ||
        capacity < message_bytes || !Span(output, capacity) ||
        !Span(source, static_cast<std::size_t>(message_bytes)) ||
        !Span(lookup, lookup_bytes) || !Span(repair_packets, packet_bytes) ||
        Overlap(output, static_cast<std::size_t>(message_bytes), source,
                static_cast<std::size_t>(message_bytes)) ||
        Overlap(output, static_cast<std::size_t>(message_bytes), lookup, lookup_bytes) ||
        Overlap(output, static_cast<std::size_t>(message_bytes), repair_packets, packet_bytes))
        return result;

    std::unique_ptr<N::Decoder> decoder;
    result.create_status = static_cast<int>(N::Decoder::Create(
        N::Lookup{lookup, lookup_bytes}, message_bytes, block_bytes, decoder));
    if (result.create_status != B::Success) return result;
    const std::uint8_t* packets = static_cast<const std::uint8_t*>(repair_packets);
    for (unsigned i = 0; i < 6; ++i) {
        const N::Result fed = decoder->Feed(6u + i,
            packets + static_cast<std::size_t>(i) * block_bytes, block_bytes);
        result.feed_status[i] = static_cast<int>(fed.status);
        result.feed_count = i + 1;
        if (fed.status == N::Status::Success) break;
        if (fed.status != N::Status::NeedMore) return result;
    }
    if (result.feed_status[result.feed_count - 1] != B::Success) return result;
    const N::Result recovered = decoder->Recover(output, capacity);
    result.recover_status = static_cast<int>(recovered.status);
    result.required = recovered.bytes_required;
    result.written = recovered.bytes_written;
    result.recovered = recovered.status == N::Status::Success &&
        recovered.bytes_written == message_bytes &&
        std::memcmp(output, source, static_cast<std::size_t>(message_bytes)) == 0;
    return result;
}

const B::Api api = {Create, Encode, EncoderFree, Row, Roundtrip};
} // namespace WH2_TINY_BRIDGE_NAMESPACE

extern "C" const wh2_tiny_payload_cost_r0::Api* WH2_TINY_BRIDGE_GETTER()
{
    return &WH2_TINY_BRIDGE_NAMESPACE::api;
}

#undef WH2_TINY_BRIDGE_NAMESPACE
#undef WH2_TINY_BRIDGE_GETTER
