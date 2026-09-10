#ifndef WIREHAIR_SMALL_CORE_H
#define WIREHAIR_SMALL_CORE_H

// Private compile-time small-block core, shared by the K3/K5 library paths
// and correctness harness (including screened K2/K8). Existing K6 does not include this file.
// Initialize the existing shared GF256 runtime before using these classes.
#include "WirehairK6Payload.h"
#include <climits>
#include <cstring>
#include <limits>
#include <memory>
#include <new>

namespace wirehair_small_core {

struct Lookup { const std::uint8_t* data; std::size_t bytes; };
enum class Status { Success, NeedMore, InvalidInput, BufferTooSmall, Conflict, OutOfMemory };
struct Result { Status status; std::size_t bytes_required, bytes_written; };

namespace detail {
inline bool Span(const void* pointer, std::size_t bytes)
{
    return pointer && bytes <= std::numeric_limits<std::uintptr_t>::max() -
        reinterpret_cast<std::uintptr_t>(pointer);
}
// Callers have checked Span before comparing potentially written ranges.
inline bool Overlap(const void* a, std::size_t an, const void* b, std::size_t bn)
{
    const std::uintptr_t x = reinterpret_cast<std::uintptr_t>(a);
    const std::uintptr_t y = reinterpret_cast<std::uintptr_t>(b);
    return an && bn && x < y + bn && y < x + an;
}
inline unsigned Parity(unsigned value)
{
    value ^= value >> 4; value ^= value >> 2; value ^= value >> 1;
    return value & 1;
}
inline Result Answer(Status status, std::size_t required, std::size_t written = 0)
{
    return Result{status, required, written};
}

template<unsigned K> struct Geometry {
    static_assert(K == 2 || K == 3 || K == 5 || K == 6 || K == 8, "Only the screened K2/K3/K5/K6/K8 dimensions are supported");
    enum : std::size_t {
        LowPhase = 1024 * K, MatrixBytes = K * K, MiddlePhase = 128 * MatrixBytes,
        Middle10 = 2 * LowPhase, Middle17 = Middle10 + 2 * MiddlePhase,
        High24 = Middle17 + 2 * MiddlePhase, LookupBytes = High24 + 256 * MatrixBytes
    };
    static bool ValidLookup(Lookup lookup)
    {
        return lookup.bytes == LookupBytes && Span(lookup.data, lookup.bytes);
    }
    static bool SystematicLookup(Lookup lookup)
    {
        if (!ValidLookup(lookup)) return false;
        for (unsigned r = 0; r < K; ++r) for (unsigned c = 0; c < K; ++c)
            if (lookup.data[r * K + c] != (r == c ? 1 : 0)) return false;
        return true;
    }
    static bool Shape(std::uint64_t message, std::uint32_t block, std::uint32_t& tail)
    {
        if (!block || block > static_cast<std::uint32_t>(INT_MAX) ||
            message > std::numeric_limits<std::size_t>::max() ||
            static_cast<std::uint64_t>(block) * (K + 1) > std::numeric_limits<std::size_t>::max() ||
            message <= static_cast<std::uint64_t>(block) * (K - 1) ||
            message > static_cast<std::uint64_t>(block) * K) return false;
        tail = static_cast<std::uint32_t>(message - static_cast<std::uint64_t>(block) * (K - 1));
        return true;
    }
    static void Apply(const std::uint8_t* matrix, std::uint8_t vector[K])
    {
        std::uint8_t result[K];
        for (unsigned r = 0; r < K; ++r) {
            std::uint8_t value = gf256_mul(vector[0], matrix[r * K]);
            for (unsigned c = 1; c < K; ++c) value ^= gf256_mul(vector[c], matrix[r * K + c]);
            result[r] = value;
        }
        std::memcpy(vector, result, K);
    }
    static void Map(Lookup lookup, std::uint32_t id, std::uint8_t output[K])
    {
        if (id < 1024) { std::memcpy(output, lookup.data + id * K, K); return; }
        const unsigned high = id >> 24, mid17 = (id >> 17) & 127, mid10 = (id >> 10) & 127;
        const unsigned phase17 = Parity(high), phase10 = phase17 ^ Parity(mid17);
        const unsigned phase0 = phase10 ^ Parity(mid10);
        std::memcpy(output, lookup.data + phase0 * LowPhase + (id & 1023) * K, K);
        if (mid10) Apply(lookup.data + Middle10 + phase10 * MiddlePhase + mid10 * MatrixBytes, output);
        if (mid17) Apply(lookup.data + Middle17 + phase17 * MiddlePhase + mid17 * MatrixBytes, output);
        if (high) Apply(lookup.data + High24 + high * MatrixBytes, output);
    }
};
} // namespace detail

template<unsigned K> Status Row(Lookup lookup, std::uint32_t id, std::uint8_t output[K])
{
    if (!detail::Geometry<K>::ValidLookup(lookup) || !detail::Span(output, K) ||
        detail::Overlap(output, K, lookup.data, lookup.bytes)) return Status::InvalidInput;
    detail::Geometry<K>::Map(lookup, id, output);
    return Status::Success;
}

// Source and immutable lookup must outlive the handle. Repairs can read the
// borrowed source; only a partial last block is copied/padded privately.
// Output handle must be empty. Every failed Create preserves it.
template<unsigned K> class Encoder {
    static_assert(K == 2 || K == 3 || K == 5 || K == 6 || K == 8, "Only K2/K3/K5/K6/K8 handles are supported");
    typedef detail::Geometry<K> G;
public:
    static Status Create(Lookup lookup, const void* source, std::uint64_t message,
                         std::uint32_t block, std::unique_ptr<Encoder>& output)
    {
        std::uint32_t tail = 0;
        if (output || !G::Shape(message, block, tail) || !detail::Span(source, static_cast<std::size_t>(message)) ||
            !G::ValidLookup(lookup) ||
            detail::Overlap(&output, sizeof(output), source, static_cast<std::size_t>(message)) ||
            detail::Overlap(&output, sizeof(output), lookup.data, lookup.bytes) || !G::SystematicLookup(lookup))
            return Status::InvalidInput;
        std::unique_ptr<Encoder> encoder(new (std::nothrow) Encoder(lookup, source, static_cast<std::size_t>(message), block, tail));
        if (!encoder) return Status::OutOfMemory;
        if (tail != block) {
            encoder->padding_.reset(new (std::nothrow) std::uint8_t[block]);
            if (!encoder->padding_) return Status::OutOfMemory;
            std::memcpy(encoder->padding_.get(), encoder->sources_[K - 1], tail);
            std::memset(encoder->padding_.get() + tail, 0, block - tail);
            encoder->sources_[K - 1] = encoder->padding_.get();
        }
        output.swap(encoder);
        return Status::Success;
    }
    ~Encoder() = default;
    Encoder(const Encoder&) = delete;
    Encoder& operator=(const Encoder&) = delete;

    // ID K-1 carries the meaningful tail. Other IDs carry one full block.
    // Capacity/alias failures never write; no Encode allocation.
    Result Encode(std::uint32_t id, void* output, std::size_t capacity) const
    {
        const std::size_t required = id == K - 1 ? tail_bytes_ : block_bytes_;
        if (capacity < required) return detail::Answer(Status::BufferTooSmall, required);
        if (!detail::Span(output, capacity) || detail::Overlap(output, required, this, sizeof(*this)) ||
            detail::Overlap(output, required, source_, message_bytes_) ||
            detail::Overlap(output, required, lookup_.data, lookup_.bytes) ||
            (padding_ && detail::Overlap(output, required, padding_.get(), block_bytes_)))
            return detail::Answer(Status::InvalidInput, required);
        if (id < K) std::memcpy(output, source_ + static_cast<std::size_t>(id) * block_bytes_, required);
        else {
            std::uint8_t coefficients[K];
            G::Map(lookup_, id, coefficients);
            wirehair_k6_payload::Payload(output, sources_, coefficients, K, static_cast<int>(block_bytes_));
        }
        return detail::Answer(Status::Success, required, required);
    }
private:
    Encoder(Lookup lookup, const void* source, std::size_t message, std::uint32_t block, std::uint32_t tail)
        : lookup_(lookup), source_(static_cast<const std::uint8_t*>(source)), message_bytes_(message),
          block_bytes_(block), tail_bytes_(tail)
    {
        for (unsigned i = 0; i < K; ++i) sources_[i] = source_ + static_cast<std::size_t>(i) * block;
    }
    Lookup lookup_;
    const std::uint8_t* source_;
    std::size_t message_bytes_;
    std::uint32_t block_bytes_, tail_bytes_;
    const void* sources_[K];
    std::unique_ptr<std::uint8_t[]> padding_;
};

// Same private core contract as the existing K6 implementation: dependent
// contradictions return Conflict without poisoning the retained basis. A
// public facade may add permanent poison; do not substitute this for WH2's
// public semantics. There is no packet-ID ledger or receive allocation.
template<unsigned K> class Decoder {
    static_assert(K == 2 || K == 3 || K == 5 || K == 6 || K == 8, "Only K2/K3/K5/K6/K8 handles are supported");
    typedef detail::Geometry<K> G;
public:
    static Status Create(Lookup lookup, std::uint64_t message, std::uint32_t block,
                         std::unique_ptr<Decoder>& output)
    {
        std::uint32_t tail = 0;
        if (output || !G::Shape(message, block, tail) || !G::ValidLookup(lookup) ||
            detail::Overlap(&output, sizeof(output), lookup.data, lookup.bytes) || !G::SystematicLookup(lookup))
            return Status::InvalidInput;
        std::unique_ptr<Decoder> decoder(new (std::nothrow) Decoder(lookup, static_cast<std::size_t>(message), block, tail));
        if (!decoder) return Status::OutOfMemory;
        decoder->slab_.reset(new (std::nothrow) std::uint8_t[decoder->slab_bytes_]);
        if (!decoder->slab_) return Status::OutOfMemory;
        for (unsigned i = 0; i < K; ++i) decoder->rhs_[i] = decoder->slab_.get() + static_cast<std::size_t>(i) * block;
        decoder->scratch_ = decoder->slab_.get() + static_cast<std::size_t>(K) * block;
        output.swap(decoder);
        return Status::Success;
    }
    ~Decoder() = default;
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;
    unsigned Rank() const { return rank_; }

    Result Feed(std::uint32_t id, const void* input, std::size_t bytes)
    {
        const std::size_t required = id == K - 1 ? tail_bytes_ : block_bytes_;
        if (bytes != required || !detail::Span(input, bytes) ||
            detail::Overlap(input, bytes, this, sizeof(*this)) || detail::Overlap(input, bytes, slab_.get(), slab_bytes_))
            return detail::Answer(Status::InvalidInput, required);
        std::uint8_t row[K]; G::Map(lookup_, id, row);
        std::memcpy(scratch_, input, bytes);
        if (bytes < block_bytes_) std::memset(scratch_ + bytes, 0, block_bytes_ - bytes);
        for (unsigned p = 0; p < K; ++p) {
            if (!(pivot_mask_ & (1u << p)) || !row[p]) continue;
            const std::uint8_t factor = row[p];
            for (unsigned c = p; c < K; ++c) row[c] ^= gf256_mul(coefficients_[p][c], factor);
            gf256_muladd_mem(scratch_, factor, rhs_[p], static_cast<int>(block_bytes_));
        }
        unsigned pivot = 0;
        while (pivot < K && !row[pivot]) ++pivot;
        if (pivot == K) {
            for (std::uint32_t i = 0; i < block_bytes_; ++i)
                if (scratch_[i]) return detail::Answer(Status::Conflict, required);
            return detail::Answer(rank_ == K ? Status::Success : Status::NeedMore, required);
        }
        const std::uint8_t inverse = gf256_inv(row[pivot]);
        for (unsigned c = pivot; c < K; ++c) row[c] = gf256_mul(row[c], inverse);
        if (inverse != 1) gf256_mul_mem(scratch_, scratch_, inverse, static_cast<int>(block_bytes_));
        std::memcpy(coefficients_[pivot], row, K);
        std::uint8_t* previous = rhs_[pivot]; rhs_[pivot] = scratch_; scratch_ = previous;
        pivot_mask_ = static_cast<std::uint8_t>(pivot_mask_ | (1u << pivot));
        ++rank_;
        return detail::Answer(rank_ == K ? Status::Success : Status::NeedMore, required);
    }

    // Solve is deferred to Recover; successful recovery changes the echelon
    // to identity and is repeatable. Failure preserves basis and output.
    Result Recover(void* output, std::size_t capacity)
    {
        if (capacity < message_bytes_) return detail::Answer(Status::BufferTooSmall, message_bytes_);
        if (!detail::Span(output, capacity) || detail::Overlap(output, message_bytes_, this, sizeof(*this)) ||
            detail::Overlap(output, message_bytes_, slab_.get(), slab_bytes_) ||
            detail::Overlap(output, message_bytes_, lookup_.data, lookup_.bytes))
            return detail::Answer(Status::InvalidInput, message_bytes_);
        if (rank_ != K) return detail::Answer(Status::NeedMore, message_bytes_);
        if (!solved_) {
            for (int p = static_cast<int>(K) - 1; p >= 0; --p) for (int r = 0; r < p; ++r) {
                const std::uint8_t factor = coefficients_[r][p];
                if (!factor) continue;
                gf256_muladd_mem(rhs_[r], factor, rhs_[p], static_cast<int>(block_bytes_));
                coefficients_[r][p] = 0;
            }
            solved_ = true;
        }
        std::uint8_t* destination = static_cast<std::uint8_t*>(output);
        for (unsigned i = 0; i < K; ++i)
            std::memcpy(destination + static_cast<std::size_t>(i) * block_bytes_, rhs_[i],
                        i == K - 1 ? tail_bytes_ : block_bytes_);
        return detail::Answer(Status::Success, message_bytes_, message_bytes_);
    }
private:
    Decoder(Lookup lookup, std::size_t message, std::uint32_t block, std::uint32_t tail)
        : lookup_(lookup), message_bytes_(message), block_bytes_(block), tail_bytes_(tail),
          slab_bytes_(static_cast<std::size_t>(block) * (K + 1)), pivot_mask_(0), rank_(0), solved_(false), scratch_(nullptr)
    {
        std::memset(coefficients_, 0, sizeof(coefficients_));
        for (unsigned i = 0; i < K; ++i) rhs_[i] = nullptr;
    }
    Lookup lookup_;
    std::size_t message_bytes_;
    std::uint32_t block_bytes_, tail_bytes_;
    std::size_t slab_bytes_;
    std::uint8_t coefficients_[K][K], pivot_mask_;
    unsigned rank_;
    bool solved_;
    std::unique_ptr<std::uint8_t[]> slab_;
    std::uint8_t* rhs_[K];
    std::uint8_t* scratch_;
};
} // namespace wirehair_small_core
#endif
