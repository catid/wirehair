#ifndef WH2_SMALL_SERIALIZED_CORE_H
#define WH2_SMALL_SERIALIZED_CORE_H

// Private compile-time generalization of the qualified K6 C boundary. This
// file is used only by benchmark builds, never by the installed library.
#include "Wh2SmallSerialized.h"
#include "../codec/WirehairSmallCore.h"

namespace wh2_small_serialized {
namespace N = wirehair_small_core;

template<unsigned K, class Traits> class Facade {
    static_assert(K == 2 || K == 3 || K == 4 || K == 5 || K == 6 || K == 8, "Only screened K2/K3/K4/K5/K6/K8 dimensions");
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
    static Wh2SmallStatus Parse(const void* profile, size_t bytes, Shape& shape)
    {
        if (bytes != WH2_SMALL_PROFILE_BYTES || !Span(profile, bytes)) return Wh2Small_InvalidInput;
        const uint8_t* p = static_cast<const uint8_t*>(profile);
        if (std::memcmp(p, "WHK", 3) || p[3] != '0' + K || Get(p + 4, 2) != 1 || Get(p + 6, 2) != 32 ||
            Get(p + 8, 8) != ProfileId || Get(p + 28, 4)) return Wh2Small_UnsupportedProfile;
        shape = Shape{Get(p + 16, 8), static_cast<uint32_t>(Get(p + 24, 4))};
        return Dimensions(shape) ? Wh2Small_Success : Wh2Small_InvalidDimensions;
    }
    static Wh2SmallStatus Status(N::Status status)
    {
        switch (status) {
        case N::Status::Success: return Wh2Small_Success;
        case N::Status::NeedMore: return Wh2Small_NeedMore;
        case N::Status::InvalidInput: return Wh2Small_InvalidInput;
        case N::Status::BufferTooSmall: return Wh2Small_BufferTooSmall;
        case N::Status::Conflict: return Wh2Small_Conflict;
        case N::Status::OutOfMemory: return Wh2Small_OutOfMemory;
        }
        return Wh2Small_InvalidInput;
    }
    static Wh2SmallResult Answer(Wh2SmallStatus s, uint64_t required = 0, uint64_t written = 0)
    { return Wh2SmallResult{s, required, written}; }
    static Wh2SmallCreateResult Failed(Wh2SmallStatus s) { return Wh2SmallCreateResult{s, nullptr}; }
    static Wh2SmallCreateResult CreateEncoder(const void* source, Shape shape, uint32_t policy)
    {
        if (policy != Wh2Small_Independent && policy != Wh2Small_BorrowedImmutable) return Failed(Wh2Small_InvalidInput);
        if (!Dimensions(shape)) return Failed(Wh2Small_InvalidDimensions);
        if (!Span(source, static_cast<size_t>(shape.message))) return Failed(Wh2Small_InvalidInput);
        std::unique_ptr<Handle> result(new (std::nothrow) Handle);
        if (!result) return Failed(Wh2Small_OutOfMemory);
        result->shape = shape;
        if (policy == Wh2Small_Independent) {
            result->source.reset(new (std::nothrow) uint8_t[static_cast<size_t>(shape.message)]);
            if (!result->source) return Failed(Wh2Small_OutOfMemory);
            std::memcpy(result->source.get(), source, static_cast<size_t>(shape.message));
            source = result->source.get();
        } else result->borrowed = static_cast<const uint8_t*>(source);
        const N::Status s = Encoder::Create(Traits::View(), source, shape.message, shape.block, result->encoder);
        return s == N::Status::Success ? Wh2SmallCreateResult{Wh2Small_Success, result.release()} : Failed(Status(s));
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

    static Wh2SmallStatus ProfileValidate(const void* p, size_t n) noexcept
    { Shape shape = {}; return Parse(p, n, shape); }

    static Wh2SmallCreateResult EncoderCreate(const void* source, uint64_t message, uint32_t block,
        uint32_t policy, void* profile, size_t capacity) noexcept
    {
        if (policy != Wh2Small_Independent && policy != Wh2Small_BorrowedImmutable) return Failed(Wh2Small_InvalidInput);
        const Shape shape = {message, block};
        if (!Dimensions(shape)) return Failed(Wh2Small_InvalidDimensions);
        if (!Span(source, static_cast<size_t>(message))) return Failed(Wh2Small_InvalidInput);
        if (policy == Wh2Small_BorrowedImmutable && profile &&
            Overlap(profile, WH2_SMALL_PROFILE_BYTES, source, static_cast<size_t>(message)))
            return Failed(Wh2Small_InvalidInput);
        if (capacity < WH2_SMALL_PROFILE_BYTES) return Failed(Wh2Small_BufferTooSmall);
        const N::Lookup lookup = Traits::View();
        if (!Span(profile, capacity) || Overlap(profile, WH2_SMALL_PROFILE_BYTES, lookup.data, lookup.bytes))
            return Failed(Wh2Small_InvalidInput);
        const Wh2SmallCreateResult result = CreateEncoder(source, shape, policy);
        if (result.status != Wh2Small_Success) return result;
        uint8_t encoded[WH2_SMALL_PROFILE_BYTES] = {};
        std::memcpy(encoded, "WHK", 3); encoded[3] = '0' + K;
        Put(encoded + 4, 1, 2); Put(encoded + 6, sizeof(encoded), 2);
        Put(encoded + 8, ProfileId, 8); Put(encoded + 16, message, 8); Put(encoded + 24, block, 4);
        std::memcpy(profile, encoded, sizeof(encoded));
        return result;
    }
    static Wh2SmallCreateResult EncoderCreateProfile(const void* source, const void* profile,
        size_t bytes, uint32_t policy) noexcept
    {
        Shape shape = {};
        const Wh2SmallStatus s = Parse(profile, bytes, shape);
        return s == Wh2Small_Success ? CreateEncoder(source, shape, policy) : Failed(s);
    }
    static Wh2SmallCreateResult DecoderCreate(const void* profile, size_t bytes) noexcept
    {
        Shape shape = {};
        const Wh2SmallStatus s = Parse(profile, bytes, shape);
        if (s != Wh2Small_Success) return Failed(s);
        std::unique_ptr<Handle> result(new (std::nothrow) Handle);
        if (!result) return Failed(Wh2Small_OutOfMemory);
        result->shape = shape;
        const N::Status created = Decoder::Create(Traits::View(), shape.message, shape.block, result->decoder);
        return created == N::Status::Success ? Wh2SmallCreateResult{Wh2Small_Success, result.release()} : Failed(Status(created));
    }
    static Wh2SmallStatus Detach(Wh2SmallCodec handle) noexcept
    {
        Handle* codec = static_cast<Handle*>(handle);
        if (!codec || !codec->encoder) return Wh2Small_InvalidInput;
        if (!codec->borrowed) return Wh2Small_Success;
        std::unique_ptr<uint8_t[]> source(new (std::nothrow) uint8_t[static_cast<size_t>(codec->shape.message)]);
        if (!source) return Wh2Small_OutOfMemory;
        std::memcpy(source.get(), codec->borrowed, static_cast<size_t>(codec->shape.message));
        std::unique_ptr<Encoder> encoder;
        const N::Status s = Encoder::Create(Traits::View(), source.get(), codec->shape.message, codec->shape.block, encoder);
        if (s != N::Status::Success) return Status(s);
        codec->encoder.swap(encoder); codec->source.swap(source); codec->borrowed = nullptr;
        return Wh2Small_Success;
    }
    static Wh2SmallResult Encode(Wh2SmallCodec handle, uint32_t id, void* output, size_t capacity) noexcept
    {
        Handle* codec = static_cast<Handle*>(handle);
        if (!codec || !codec->encoder) return Answer(Wh2Small_InvalidInput);
        const size_t required = PacketBytes(codec, id);
        if (!Writable(codec, output, capacity, required)) return Answer(Wh2Small_InvalidInput, required);
        const N::Result r = codec->encoder->Encode(id, output, capacity);
        return Answer(Status(r.status), r.bytes_required, r.bytes_written);
    }
    static Wh2SmallStatus Decode(Wh2SmallCodec handle, uint32_t id, const void* input, size_t bytes) noexcept
    {
        Handle* codec = static_cast<Handle*>(handle);
        if (!codec || !codec->decoder || !Span(input, bytes) || Overlap(input, bytes, codec, sizeof(*codec)))
            return Wh2Small_InvalidInput;
        if (bytes != PacketBytes(codec, id)) return Wh2Small_InvalidInput;
        if (codec->poisoned) return Wh2Small_Conflict;
        const N::Status s = codec->decoder->Feed(id, input, bytes).status;
        if (s == N::Status::Conflict) codec->poisoned = true;
        return Status(s);
    }
    static Wh2SmallResult Recover(Wh2SmallCodec handle, void* output, size_t capacity) noexcept
    {
        Handle* codec = static_cast<Handle*>(handle);
        if (!codec || !codec->decoder) return Answer(Wh2Small_InvalidInput);
        const size_t required = static_cast<size_t>(codec->shape.message);
        if (!Writable(codec, output, capacity, required)) return Answer(Wh2Small_InvalidInput, required);
        if (capacity < required) return Answer(Wh2Small_BufferTooSmall, required);
        if (codec->poisoned) return Answer(Wh2Small_Conflict, required);
        const N::Result r = codec->decoder->Recover(output, capacity);
        return Answer(Status(r.status), r.bytes_required, r.bytes_written);
    }
    static void Free(Wh2SmallCodec handle) noexcept { delete static_cast<Handle*>(handle); }
};
} // namespace wh2_small_serialized
#endif
