#ifndef WIREHAIR_SMALL_FACADE_H
#define WIREHAIR_SMALL_FACADE_H

// Reusable private small-block facade, promoted from the qualified prototype.
// Only K3 is instantiated by the installed small-block API. Existing K6 is
// unchanged. The opaque public pointer is converted back before any access.
#include "wirehair/wirehair_small.h"
#include "WirehairSmallCore.h"

namespace wirehair_small_facade {
namespace N = wirehair_small_core;

template<unsigned K, class Traits> class Facade {
    static_assert(K == 3 || K == 6, "Only screened K3/K6 dimensions");
    typedef N::Encoder<K> Encoder;
    typedef N::Decoder<K> Decoder;
    struct Shape { uint64_t message; uint32_t block; };
    struct Handle {
        Shape shape;
        const uint8_t* borrowed = nullptr;
        bool poisoned = false;
        // Reverse destruction keeps owned source alive through Encoder.
        std::unique_ptr<uint8_t[]> source;
        std::unique_ptr<Encoder> encoder;
        std::unique_ptr<Decoder> decoder;
    };
    static bool Span(const void* p, size_t n) { return N::detail::Span(p, n); }
    static bool Overlap(const void* a, size_t an, const void* b, size_t bn)
    {
        const uintptr_t x = reinterpret_cast<uintptr_t>(a), y = reinterpret_cast<uintptr_t>(b);
        return an && bn && (x <= y ? y - x < an : x - y < bn);
    }
    static bool Dimensions(Shape s)
    {
        return s.block && s.block <= MaxBlockBytes &&
            uint64_t(s.block) * (K + 1) <= std::numeric_limits<size_t>::max() &&
            s.message > uint64_t(s.block) * (K - 1) && s.message <= uint64_t(s.block) * K;
    }
    static uint64_t Get(const uint8_t* p, unsigned bytes)
    {
        uint64_t value = 0;
        for (unsigned i = 0; i < bytes; ++i) value |= uint64_t(p[i]) << (8 * i);
        return value;
    }
    static void Put(uint8_t* p, uint64_t value, unsigned bytes)
    {
        for (unsigned i = 0; i < bytes; ++i) p[i] = static_cast<uint8_t>(value >> (8 * i));
    }
    static WirehairSmallStatus Parse(const void* profile, size_t bytes, Shape& shape)
    {
        if (bytes != WIREHAIR_SMALL_PROFILE_BYTES || !Span(profile, bytes)) return WirehairSmall_InvalidInput;
        const uint8_t* p = static_cast<const uint8_t*>(profile);
        if (std::memcmp(p, "WHK", 3) || p[3] != '0' + K || Get(p + 4, 2) != 1 || Get(p + 6, 2) != 32 ||
            Get(p + 8, 8) != ProfileId || Get(p + 28, 4)) return WirehairSmall_UnsupportedProfile;
        shape = Shape{Get(p + 16, 8), static_cast<uint32_t>(Get(p + 24, 4))};
        return Dimensions(shape) ? WirehairSmall_Success : WirehairSmall_InvalidDimensions;
    }
    static WirehairSmallStatus Status(N::Status status)
    {
        switch (status) {
        case N::Status::Success: return WirehairSmall_Success;
        case N::Status::NeedMore: return WirehairSmall_NeedMore;
        case N::Status::InvalidInput: return WirehairSmall_InvalidInput;
        case N::Status::BufferTooSmall: return WirehairSmall_BufferTooSmall;
        case N::Status::Conflict: return WirehairSmall_Conflict;
        case N::Status::OutOfMemory: return WirehairSmall_OutOfMemory;
        }
        return WirehairSmall_InvalidInput;
    }
    static WirehairSmallResult Answer(WirehairSmallStatus s, uint64_t required = 0, uint64_t written = 0)
    { return WirehairSmallResult{s, required, written}; }
    static WirehairSmallCreateResult Failed(WirehairSmallStatus s) { return WirehairSmallCreateResult{s, nullptr}; }
    static WirehairSmallCreateResult CreateEncoder(const void* source, Shape shape, uint32_t policy)
    {
        if (policy != WirehairSmall_Independent && policy != WirehairSmall_BorrowedImmutable) return Failed(WirehairSmall_InvalidInput);
        if (!Dimensions(shape)) return Failed(WirehairSmall_InvalidDimensions);
        if (!Span(source, static_cast<size_t>(shape.message))) return Failed(WirehairSmall_InvalidInput);
        std::unique_ptr<Handle> result(new (std::nothrow) Handle);
        if (!result) return Failed(WirehairSmall_OutOfMemory);
        result->shape = shape;
        if (policy == WirehairSmall_Independent) {
            result->source.reset(new (std::nothrow) uint8_t[static_cast<size_t>(shape.message)]);
            if (!result->source) return Failed(WirehairSmall_OutOfMemory);
            std::memcpy(result->source.get(), source, static_cast<size_t>(shape.message));
            source = result->source.get();
        } else result->borrowed = static_cast<const uint8_t*>(source);
        const N::Status s = Encoder::Create(Traits::View(), source, shape.message, shape.block, result->encoder);
        return s == N::Status::Success ? WirehairSmallCreateResult{WirehairSmall_Success, reinterpret_cast<WirehairSmallCodec>(result.release())} : Failed(Status(s));
    }
    static bool Writable(Handle* codec, void* output, size_t capacity, size_t required)
    {
        const N::Lookup lookup = Traits::View();
        return Span(output, capacity) && !Overlap(output, required, codec, sizeof(*codec)) &&
            !Overlap(output, required, lookup.data, lookup.bytes);
    }
    static size_t PacketBytes(Handle* codec, uint32_t id)
    {
        return id == K - 1 ? static_cast<size_t>(codec->shape.message - uint64_t(codec->shape.block) * (K - 1)) :
            codec->shape.block;
    }
public:
    static constexpr uint64_t ProfileId = UINT64_C(0x5748324b30544d31) + (uint64_t(K) << 24);
    static constexpr uint32_t MaxBlockBytes = UINT32_C(268435456) / (K + 1);

    static WirehairSmallStatus ProfileValidate(const void* p, size_t n) noexcept
    { Shape shape = {}; return Parse(p, n, shape); }

    static WirehairSmallCreateResult EncoderCreate(const void* source, uint64_t message, uint32_t block,
        uint32_t policy, void* profile, size_t capacity) noexcept
    {
        if (policy != WirehairSmall_Independent && policy != WirehairSmall_BorrowedImmutable) return Failed(WirehairSmall_InvalidInput);
        const Shape shape = {message, block};
        if (!Dimensions(shape)) return Failed(WirehairSmall_InvalidDimensions);
        if (!Span(source, static_cast<size_t>(message))) return Failed(WirehairSmall_InvalidInput);
        if (policy == WirehairSmall_BorrowedImmutable && profile &&
            Overlap(profile, WIREHAIR_SMALL_PROFILE_BYTES, source, static_cast<size_t>(message)))
            return Failed(WirehairSmall_InvalidInput);
        if (capacity < WIREHAIR_SMALL_PROFILE_BYTES) return Failed(WirehairSmall_BufferTooSmall);
        const N::Lookup lookup = Traits::View();
        if (!Span(profile, capacity) || Overlap(profile, WIREHAIR_SMALL_PROFILE_BYTES, lookup.data, lookup.bytes))
            return Failed(WirehairSmall_InvalidInput);
        const WirehairSmallCreateResult result = CreateEncoder(source, shape, policy);
        if (result.status != WirehairSmall_Success) return result;
        uint8_t encoded[WIREHAIR_SMALL_PROFILE_BYTES] = {};
        std::memcpy(encoded, "WHK", 3); encoded[3] = '0' + K;
        Put(encoded + 4, 1, 2); Put(encoded + 6, sizeof(encoded), 2);
        Put(encoded + 8, ProfileId, 8); Put(encoded + 16, message, 8); Put(encoded + 24, block, 4);
        std::memcpy(profile, encoded, sizeof(encoded));
        return result;
    }
    static WirehairSmallCreateResult EncoderCreateProfile(const void* source, const void* profile,
        size_t bytes, uint32_t policy) noexcept
    {
        Shape shape = {};
        const WirehairSmallStatus s = Parse(profile, bytes, shape);
        return s == WirehairSmall_Success ? CreateEncoder(source, shape, policy) : Failed(s);
    }
    static WirehairSmallCreateResult DecoderCreate(const void* profile, size_t bytes) noexcept
    {
        Shape shape = {};
        const WirehairSmallStatus s = Parse(profile, bytes, shape);
        if (s != WirehairSmall_Success) return Failed(s);
        std::unique_ptr<Handle> result(new (std::nothrow) Handle);
        if (!result) return Failed(WirehairSmall_OutOfMemory);
        result->shape = shape;
        const N::Status created = Decoder::Create(Traits::View(), shape.message, shape.block, result->decoder);
        return created == N::Status::Success ? WirehairSmallCreateResult{WirehairSmall_Success, reinterpret_cast<WirehairSmallCodec>(result.release())} : Failed(Status(created));
    }
    static WirehairSmallStatus Detach(WirehairSmallCodec handle) noexcept
    {
        Handle* codec = reinterpret_cast<Handle*>(handle);
        if (!codec || !codec->encoder) return WirehairSmall_InvalidInput;
        if (!codec->borrowed) return WirehairSmall_Success;
        std::unique_ptr<uint8_t[]> source(new (std::nothrow) uint8_t[static_cast<size_t>(codec->shape.message)]);
        if (!source) return WirehairSmall_OutOfMemory;
        std::memcpy(source.get(), codec->borrowed, static_cast<size_t>(codec->shape.message));
        std::unique_ptr<Encoder> encoder;
        const N::Status s = Encoder::Create(Traits::View(), source.get(), codec->shape.message, codec->shape.block, encoder);
        if (s != N::Status::Success) return Status(s);
        codec->encoder.swap(encoder); codec->source.swap(source); codec->borrowed = nullptr;
        return WirehairSmall_Success;
    }
    static WirehairSmallResult Encode(WirehairSmallCodec handle, uint32_t id, void* output, size_t capacity) noexcept
    {
        Handle* codec = reinterpret_cast<Handle*>(handle);
        if (!codec || !codec->encoder) return Answer(WirehairSmall_InvalidInput);
        const size_t required = PacketBytes(codec, id);
        if (!Writable(codec, output, capacity, required)) return Answer(WirehairSmall_InvalidInput, required);
        const N::Result r = codec->encoder->Encode(id, output, capacity);
        return Answer(Status(r.status), r.bytes_required, r.bytes_written);
    }
    static WirehairSmallStatus Decode(WirehairSmallCodec handle, uint32_t id, const void* input, size_t bytes) noexcept
    {
        Handle* codec = reinterpret_cast<Handle*>(handle);
        if (!codec || !codec->decoder || !Span(input, bytes) || Overlap(input, bytes, codec, sizeof(*codec)))
            return WirehairSmall_InvalidInput;
        if (bytes != PacketBytes(codec, id)) return WirehairSmall_InvalidInput;
        if (codec->poisoned) return WirehairSmall_Conflict;
        const N::Status s = codec->decoder->Feed(id, input, bytes).status;
        if (s == N::Status::Conflict) codec->poisoned = true;
        return Status(s);
    }
    static WirehairSmallResult Recover(WirehairSmallCodec handle, void* output, size_t capacity) noexcept
    {
        Handle* codec = reinterpret_cast<Handle*>(handle);
        if (!codec || !codec->decoder) return Answer(WirehairSmall_InvalidInput);
        const size_t required = static_cast<size_t>(codec->shape.message);
        if (!Writable(codec, output, capacity, required)) return Answer(WirehairSmall_InvalidInput, required);
        if (capacity < required) return Answer(WirehairSmall_BufferTooSmall, required);
        if (codec->poisoned) return Answer(WirehairSmall_Conflict, required);
        const N::Result r = codec->decoder->Recover(output, capacity);
        return Answer(Status(r.status), r.bytes_required, r.bytes_written);
    }
    static void Free(WirehairSmallCodec handle) noexcept { delete reinterpret_cast<Handle*>(handle); }
};
} // namespace wirehair_small_facade
#endif
