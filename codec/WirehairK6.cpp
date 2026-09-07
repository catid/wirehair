#include "wirehair/wirehair_k6.h"
#include "WirehairK6Core.h"

#include <cstring>
#include <limits>
#include <new>

namespace N = wirehair_k6_core;
namespace {
alignas(64) const uint8_t kLookup[] = {
#include "WirehairK6Lookup.inc"
};
static_assert(sizeof(kLookup) == N::kLookupBytes, "Exact selected lookup");
const N::Lookup kView = {kLookup, sizeof(kLookup)};
struct Shape { uint64_t message; uint32_t block; };

bool Span(const void* p, size_t n)
{
    return p && n <= std::numeric_limits<uintptr_t>::max() - reinterpret_cast<uintptr_t>(p);
}

bool Overlap(const void* a, size_t an, const void* b, size_t bn)
{
    const uintptr_t x = reinterpret_cast<uintptr_t>(a), y = reinterpret_cast<uintptr_t>(b);
    return an && bn && (x <= y ? y - x < an : x - y < bn);
}

bool Dimensions(Shape s)
{
    return s.block && s.block <= WIREHAIR_K6_MAX_BLOCK_BYTES &&
        uint64_t(s.block) * 7 <= std::numeric_limits<size_t>::max() &&
        s.message > uint64_t(s.block) * 5 && s.message <= uint64_t(s.block) * 6;
}

uint64_t Get(const uint8_t* p, unsigned bytes)
{
    uint64_t value = 0;
    for (unsigned i = 0; i < bytes; ++i) value |= uint64_t(p[i]) << (8 * i);
    return value;
}

void Put(uint8_t* p, uint64_t value, unsigned bytes)
{
    for (unsigned i = 0; i < bytes; ++i) p[i] = static_cast<uint8_t>(value >> (8 * i));
}

WirehairK6Status Parse(const void* profile, size_t bytes, Shape& shape)
{
    if (bytes != WIREHAIR_K6_PROFILE_BYTES || !Span(profile, bytes)) return WirehairK6_InvalidInput;
    const uint8_t* p = static_cast<const uint8_t*>(profile);
    if (std::memcmp(p, "WHK6", 4) || Get(p + 4, 2) != 1 || Get(p + 6, 2) != 32 ||
        Get(p + 8, 8) != WIREHAIR_K6_PROFILE_ID || Get(p + 28, 4)) return WirehairK6_UnsupportedProfile;
    shape = Shape{Get(p + 16, 8), static_cast<uint32_t>(Get(p + 24, 4))};
    return Dimensions(shape) ? WirehairK6_Success : WirehairK6_InvalidDimensions;
}

WirehairK6Status Status(N::Status status)
{
    switch (status) {
    case N::Status::Success: return WirehairK6_Success;
    case N::Status::NeedMore: return WirehairK6_NeedMore;
    case N::Status::InvalidInput: return WirehairK6_InvalidInput;
    case N::Status::BufferTooSmall: return WirehairK6_BufferTooSmall;
    case N::Status::Conflict: return WirehairK6_Conflict;
    case N::Status::OutOfMemory: return WirehairK6_OutOfMemory;
    }
    return WirehairK6_InvalidInput;
}

WirehairK6Result Answer(WirehairK6Status status, uint64_t required = 0, uint64_t written = 0)
{
    return WirehairK6Result{status, required, written};
}
WirehairK6CreateResult Failed(WirehairK6Status status) { return WirehairK6CreateResult{status, nullptr}; }
} // namespace

struct WirehairK6CodecImpl {
    Shape shape;
    const uint8_t* borrowed = nullptr;
    bool poisoned = false;
    // Reverse destruction order keeps source storage alive through Encoder.
    std::unique_ptr<uint8_t[]> source;
    std::unique_ptr<N::Encoder> encoder;
    std::unique_ptr<N::Decoder> decoder;
};

namespace {
WirehairK6CreateResult CreateEncoder(const void* source, Shape shape, uint32_t policy)
{
    if (policy != WirehairK6_Independent && policy != WirehairK6_BorrowedImmutable)
        return Failed(WirehairK6_InvalidInput);
    if (!Dimensions(shape)) return Failed(WirehairK6_InvalidDimensions);
    if (!Span(source, static_cast<size_t>(shape.message))) return Failed(WirehairK6_InvalidInput);
    std::unique_ptr<WirehairK6CodecImpl> result(new (std::nothrow) WirehairK6CodecImpl);
    if (!result) return Failed(WirehairK6_OutOfMemory);
    result->shape = shape;
    if (policy == WirehairK6_Independent) {
        result->source.reset(new (std::nothrow) uint8_t[static_cast<size_t>(shape.message)]);
        if (!result->source) return Failed(WirehairK6_OutOfMemory);
        std::memcpy(result->source.get(), source, static_cast<size_t>(shape.message));
        source = result->source.get();
    } else result->borrowed = static_cast<const uint8_t*>(source);
    const N::Status status = N::Encoder::Create(kView, source, shape.message,
                                               shape.block, result->encoder);
    if (status != N::Status::Success) return Failed(Status(status));
    return WirehairK6CreateResult{WirehairK6_Success, result.release()};
}

bool Writable(WirehairK6Codec codec, void* output, size_t capacity, size_t required)
{
    return Span(output, capacity) && !Overlap(output, required, codec, sizeof(*codec)) &&
        !Overlap(output, required, kLookup, sizeof(kLookup));
}
} // namespace

extern "C" WirehairK6Status wirehair_k6_profile_validate(const void* profile, size_t bytes) noexcept
{
    Shape shape = {};
    return Parse(profile, bytes, shape);
}

extern "C" WirehairK6CreateResult wirehair_k6_encoder_create(const void* source, uint64_t message,
    uint32_t block, uint32_t policy, void* profile, size_t capacity) noexcept
{
    if (policy != WirehairK6_Independent && policy != WirehairK6_BorrowedImmutable)
        return Failed(WirehairK6_InvalidInput);
    const Shape shape = {message, block};
    if (!Dimensions(shape)) return Failed(WirehairK6_InvalidDimensions);
    if (!Span(source, static_cast<size_t>(message))) return Failed(WirehairK6_InvalidInput);
    if (policy == WirehairK6_BorrowedImmutable && profile &&
        Overlap(profile, WIREHAIR_K6_PROFILE_BYTES, source, static_cast<size_t>(message)))
        return Failed(WirehairK6_InvalidInput);
    if (capacity < WIREHAIR_K6_PROFILE_BYTES) return Failed(WirehairK6_BufferTooSmall);
    if (!Span(profile, capacity) || Overlap(profile, WIREHAIR_K6_PROFILE_BYTES, kLookup, sizeof(kLookup)))
        return Failed(WirehairK6_InvalidInput);
    const WirehairK6CreateResult result = CreateEncoder(source, shape, policy);
    if (result.status != WirehairK6_Success) return result;
    uint8_t encoded[WIREHAIR_K6_PROFILE_BYTES] = {};
    std::memcpy(encoded, "WHK6", 4);
    Put(encoded + 4, 1, 2); Put(encoded + 6, sizeof(encoded), 2);
    Put(encoded + 8, WIREHAIR_K6_PROFILE_ID, 8);
    Put(encoded + 16, message, 8); Put(encoded + 24, block, 4);
    std::memcpy(profile, encoded, sizeof(encoded));
    return result;
}

extern "C" WirehairK6CreateResult wirehair_k6_encoder_create_profile(const void* source,
    const void* profile, size_t bytes, uint32_t policy) noexcept
{
    Shape shape = {};
    const WirehairK6Status status = Parse(profile, bytes, shape);
    return status == WirehairK6_Success ? CreateEncoder(source, shape, policy) : Failed(status);
}

extern "C" WirehairK6CreateResult wirehair_k6_decoder_create(const void* profile, size_t bytes) noexcept
{
    Shape shape = {};
    const WirehairK6Status status = Parse(profile, bytes, shape);
    if (status != WirehairK6_Success) return Failed(status);
    std::unique_ptr<WirehairK6CodecImpl> result(new (std::nothrow) WirehairK6CodecImpl);
    if (!result) return Failed(WirehairK6_OutOfMemory);
    result->shape = shape;
    const N::Status created = N::Decoder::Create(kView, shape.message, shape.block, result->decoder);
    if (created != N::Status::Success) return Failed(Status(created));
    return WirehairK6CreateResult{WirehairK6_Success, result.release()};
}

extern "C" WirehairK6Status wirehair_k6_encoder_detach_input(WirehairK6Codec codec) noexcept
{
    if (!codec || !codec->encoder) return WirehairK6_InvalidInput;
    if (!codec->borrowed) return WirehairK6_Success;
    std::unique_ptr<uint8_t[]> source(new (std::nothrow) uint8_t[static_cast<size_t>(codec->shape.message)]);
    if (!source) return WirehairK6_OutOfMemory;
    std::memcpy(source.get(), codec->borrowed, static_cast<size_t>(codec->shape.message));
    std::unique_ptr<N::Encoder> encoder;
    const N::Status status = N::Encoder::Create(kView, source.get(), codec->shape.message,
                                               codec->shape.block, encoder);
    if (status != N::Status::Success) return Status(status);
    codec->encoder.swap(encoder);
    codec->source.swap(source);
    codec->borrowed = nullptr;
    return WirehairK6_Success;
}

extern "C" WirehairK6Result wirehair_k6_encode(WirehairK6Codec codec, uint32_t id, void* output, size_t capacity) noexcept
{
    if (!codec || !codec->encoder) return Answer(WirehairK6_InvalidInput);
    const size_t required = id == 5 ? static_cast<size_t>(codec->shape.message - uint64_t(codec->shape.block) * 5) : codec->shape.block;
    if (!Writable(codec, output, capacity, required)) return Answer(WirehairK6_InvalidInput, required);
    const N::Result result = codec->encoder->Encode(id, output, capacity);
    return Answer(Status(result.status), result.bytes_required, result.bytes_written);
}

extern "C" WirehairK6Status wirehair_k6_decode(WirehairK6Codec codec, uint32_t id, const void* input, size_t bytes) noexcept
{
    if (!codec || !codec->decoder || !Span(input, bytes) ||
        Overlap(input, bytes, codec, sizeof(*codec))) return WirehairK6_InvalidInput;
    const size_t required = id == 5 ? static_cast<size_t>(codec->shape.message - uint64_t(codec->shape.block) * 5) : codec->shape.block;
    if (bytes != required) return WirehairK6_InvalidInput;
    if (codec->poisoned) return WirehairK6_Conflict;
    const N::Status status = codec->decoder->Feed(id, input, bytes).status;
    if (status == N::Status::Conflict) codec->poisoned = true;
    return Status(status);
}

extern "C" WirehairK6Result wirehair_k6_recover(WirehairK6Codec codec, void* output, size_t capacity) noexcept
{
    if (!codec || !codec->decoder) return Answer(WirehairK6_InvalidInput);
    const size_t required = static_cast<size_t>(codec->shape.message);
    if (!Writable(codec, output, capacity, required)) return Answer(WirehairK6_InvalidInput, required);
    if (capacity < required) return Answer(WirehairK6_BufferTooSmall, required);
    if (codec->poisoned) return Answer(WirehairK6_Conflict, required);
    const N::Result result = codec->decoder->Recover(output, capacity);
    return Answer(Status(result.status), result.bytes_required, result.bytes_written);
}

extern "C" void wirehair_k6_free(WirehairK6Codec codec) noexcept { delete codec; }
