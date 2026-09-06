#include "Wh2ThueMorseNativeCodec.h"
#include "gf256.h"

#include <climits>
#include <cstring>
#include <limits>
#include <new>

namespace wh2_thue_native {
namespace {

bool Span(const void* pointer, std::size_t bytes)
{
    return pointer && bytes <= std::numeric_limits<std::uintptr_t>::max() -
        reinterpret_cast<std::uintptr_t>(pointer);
}

// Both spans have already passed Span. Only potentially written bytes need
// alias exclusion; capacity is nevertheless checked separately for wrapping.
bool Overlap(const void* a, std::size_t a_bytes, const void* b, std::size_t b_bytes)
{
    const std::uintptr_t x = reinterpret_cast<std::uintptr_t>(a);
    const std::uintptr_t y = reinterpret_cast<std::uintptr_t>(b);
    return a_bytes && b_bytes && x < y + b_bytes && y < x + a_bytes;
}

bool ValidLookup(Lookup lookup)
{
    return lookup.bytes == kLookupBytes && Span(lookup.data, lookup.bytes);
}

bool SystematicLookup(Lookup lookup)
{
    if (!ValidLookup(lookup)) return false;
    for (unsigned r = 0; r < kSourceCount; ++r)
        for (unsigned c = 0; c < kSourceCount; ++c)
            if (lookup.data[r * kSourceCount + c] != (r == c ? 1 : 0)) return false;
    return true;
}

bool Shape(std::uint64_t message, std::uint32_t block, std::uint32_t& tail)
{
    if (!block || block > static_cast<std::uint32_t>(INT_MAX) ||
        message > std::numeric_limits<std::size_t>::max() ||
        static_cast<std::uint64_t>(block) * 7 > std::numeric_limits<std::size_t>::max() ||
        message <= static_cast<std::uint64_t>(block) * 5 ||
        message > static_cast<std::uint64_t>(block) * 6) return false;
    tail = static_cast<std::uint32_t>(message - static_cast<std::uint64_t>(block) * 5);
    return true;
}

unsigned Parity(unsigned value)
{
    value ^= value >> 4;
    value ^= value >> 2;
    value ^= value >> 1;
    return value & 1;
}

void Apply(const std::uint8_t* matrix, std::uint8_t vector[kSourceCount])
{
    std::uint8_t result[kSourceCount];
    for (unsigned r = 0; r < kSourceCount; ++r) {
        std::uint8_t value = gf256_mul(vector[0], matrix[r * kSourceCount]);
        for (unsigned c = 1; c < kSourceCount; ++c)
            value ^= gf256_mul(vector[c], matrix[r * kSourceCount + c]);
        result[r] = value;
    }
    std::memcpy(vector, result, kSourceCount);
}

void Map(Lookup lookup, std::uint32_t id, std::uint8_t output[kSourceCount])
{
    if (id < 1024) {
        std::memcpy(output, lookup.data + id * kSourceCount, kSourceCount);
        return;
    }
    const unsigned high = id >> 24;
    const unsigned mid17 = (id >> 17) & 127;
    const unsigned mid10 = (id >> 10) & 127;
    const unsigned low = id & 1023;
    const unsigned phase17 = Parity(high);
    const unsigned phase10 = phase17 ^ Parity(mid17);
    const unsigned phase0 = phase10 ^ Parity(mid10);
    std::memcpy(output, lookup.data + phase0 * 6144 + low * kSourceCount, kSourceCount);
    if (mid10) Apply(lookup.data + 12288 + phase10 * 4608 + mid10 * 36, output);
    if (mid17) Apply(lookup.data + 21504 + phase17 * 4608 + mid17 * 36, output);
    if (high) Apply(lookup.data + 30720 + high * 36, output);
}

Result Answer(Status status, std::size_t required, std::size_t written = 0)
{
    return Result{status, required, written};
}

} // namespace

Status Row(Lookup lookup, std::uint32_t id, std::uint8_t output[kSourceCount])
{
    if (!ValidLookup(lookup) || !Span(output, kSourceCount) ||
        Overlap(output, kSourceCount, lookup.data, lookup.bytes)) return Status::InvalidInput;
    Map(lookup, id, output);
    return Status::Success;
}

Encoder::Encoder(Lookup lookup, const void* source, std::size_t message_bytes,
                 std::uint32_t block_bytes, std::uint32_t tail_bytes)
    : lookup_(lookup), source_(static_cast<const std::uint8_t*>(source)),
      message_bytes_(message_bytes), block_bytes_(block_bytes), tail_bytes_(tail_bytes)
{
    for (unsigned i = 0; i < kSourceCount; ++i)
        sources_[i] = source_ + static_cast<std::size_t>(i) * block_bytes;
}

Encoder::~Encoder() = default;

Status Encoder::Create(Lookup lookup, const void* source, std::uint64_t message_bytes,
                       std::uint32_t block_bytes, std::unique_ptr<Encoder>& output)
{
    std::uint32_t tail = 0;
    if (output || !Shape(message_bytes, block_bytes, tail) ||
        !Span(source, static_cast<std::size_t>(message_bytes)) ||
        !ValidLookup(lookup) ||
        Overlap(&output, sizeof(output), source, static_cast<std::size_t>(message_bytes)) ||
        Overlap(&output, sizeof(output), lookup.data, lookup.bytes) ||
        !SystematicLookup(lookup)) return Status::InvalidInput;
    std::unique_ptr<Encoder> encoder(new (std::nothrow) Encoder(
        lookup, source, static_cast<std::size_t>(message_bytes), block_bytes, tail));
    if (!encoder) return Status::OutOfMemory;
    if (tail != block_bytes) {
        encoder->padding_.reset(new (std::nothrow) std::uint8_t[block_bytes]);
        if (!encoder->padding_) return Status::OutOfMemory;
        std::memcpy(encoder->padding_.get(), encoder->sources_[5], tail);
        std::memset(encoder->padding_.get() + tail, 0, block_bytes - tail);
        encoder->sources_[5] = encoder->padding_.get();
    }
    output.swap(encoder);
    return Status::Success;
}

Result Encoder::Encode(std::uint32_t id, void* output, std::size_t capacity) const
{
    const std::size_t required = id == 5 ? tail_bytes_ : block_bytes_;
    if (capacity < required) return Answer(Status::BufferTooSmall, required);
    if (!Span(output, capacity) || Overlap(output, required, this, sizeof(*this)) ||
        Overlap(output, required, source_, message_bytes_) ||
        Overlap(output, required, lookup_.data, lookup_.bytes) ||
        (padding_ && Overlap(output, required, padding_.get(), block_bytes_)))
        return Answer(Status::InvalidInput, required);
    if (id < kSourceCount) {
        std::memcpy(output, source_ + static_cast<std::size_t>(id) * block_bytes_, required);
    } else {
        std::uint8_t coefficients[kSourceCount];
        Map(lookup_, id, coefficients);
        gf256_mulset_multi_mem(output, sources_, coefficients, kSourceCount,
                              static_cast<int>(block_bytes_));
    }
    return Answer(Status::Success, required, required);
}

Decoder::Decoder(Lookup lookup, std::size_t message_bytes,
                 std::uint32_t block_bytes, std::uint32_t tail_bytes)
    : lookup_(lookup), message_bytes_(message_bytes), block_bytes_(block_bytes),
      tail_bytes_(tail_bytes), slab_bytes_(static_cast<std::size_t>(block_bytes) * 7),
      pivot_mask_(0), rank_(0), solved_(false), scratch_(nullptr)
{
    std::memset(coefficients_, 0, sizeof(coefficients_));
    for (unsigned i = 0; i < kSourceCount; ++i) rhs_[i] = nullptr;
}

Decoder::~Decoder() = default;

Status Decoder::Create(Lookup lookup, std::uint64_t message_bytes,
                       std::uint32_t block_bytes, std::unique_ptr<Decoder>& output)
{
    std::uint32_t tail = 0;
    if (output || !Shape(message_bytes, block_bytes, tail) || !ValidLookup(lookup) ||
        Overlap(&output, sizeof(output), lookup.data, lookup.bytes) ||
        !SystematicLookup(lookup)) return Status::InvalidInput;
    std::unique_ptr<Decoder> decoder(new (std::nothrow) Decoder(
        lookup, static_cast<std::size_t>(message_bytes), block_bytes, tail));
    if (!decoder) return Status::OutOfMemory;
    decoder->slab_.reset(new (std::nothrow) std::uint8_t[decoder->slab_bytes_]);
    if (!decoder->slab_) return Status::OutOfMemory;
    for (unsigned i = 0; i < kSourceCount; ++i)
        decoder->rhs_[i] = decoder->slab_.get() + static_cast<std::size_t>(i) * block_bytes;
    decoder->scratch_ = decoder->slab_.get() + static_cast<std::size_t>(kSourceCount) * block_bytes;
    output.swap(decoder);
    return Status::Success;
}

Result Decoder::Feed(std::uint32_t id, const void* input, std::size_t bytes)
{
    const std::size_t required = id == 5 ? tail_bytes_ : block_bytes_;
    if (bytes != required || !Span(input, bytes) ||
        Overlap(input, bytes, this, sizeof(*this)) ||
        Overlap(input, bytes, slab_.get(), slab_bytes_)) return Answer(Status::InvalidInput, required);
    std::uint8_t row[kSourceCount];
    Map(lookup_, id, row);
    std::memcpy(scratch_, input, bytes);
    if (bytes < block_bytes_) std::memset(scratch_ + bytes, 0, block_bytes_ - bytes);
    for (unsigned p = 0; p < kSourceCount; ++p) {
        if (!(pivot_mask_ & (1u << p)) || !row[p]) continue;
        const std::uint8_t factor = row[p];
        for (unsigned c = p; c < kSourceCount; ++c)
            row[c] ^= gf256_mul(coefficients_[p][c], factor);
        gf256_muladd_mem(scratch_, factor, rhs_[p], static_cast<int>(block_bytes_));
    }
    unsigned pivot = 0;
    while (pivot < kSourceCount && !row[pivot]) ++pivot;
    if (pivot == kSourceCount) {
        for (std::uint32_t i = 0; i < block_bytes_; ++i)
            if (scratch_[i]) return Answer(Status::Conflict, required);
        return Answer(rank_ == kSourceCount ? Status::Success : Status::NeedMore, required);
    }
    const std::uint8_t inverse = gf256_inv(row[pivot]);
    for (unsigned c = pivot; c < kSourceCount; ++c) row[c] = gf256_mul(row[c], inverse);
    if (inverse != 1) gf256_mul_mem(scratch_, scratch_, inverse, static_cast<int>(block_bytes_));
    std::memcpy(coefficients_[pivot], row, kSourceCount);
    std::uint8_t* previous = rhs_[pivot];
    rhs_[pivot] = scratch_;
    scratch_ = previous;
    pivot_mask_ = static_cast<std::uint8_t>(pivot_mask_ | (1u << pivot));
    ++rank_;
    return Answer(rank_ == kSourceCount ? Status::Success : Status::NeedMore, required);
}

Result Decoder::Recover(void* output, std::size_t capacity)
{
    if (capacity < message_bytes_) return Answer(Status::BufferTooSmall, message_bytes_);
    if (!Span(output, capacity) || Overlap(output, message_bytes_, this, sizeof(*this)) ||
        Overlap(output, message_bytes_, slab_.get(), slab_bytes_) ||
        Overlap(output, message_bytes_, lookup_.data, lookup_.bytes))
        return Answer(Status::InvalidInput, message_bytes_);
    if (rank_ != kSourceCount) return Answer(Status::NeedMore, message_bytes_);
    if (!solved_) {
        for (int p = static_cast<int>(kSourceCount) - 1; p >= 0; --p) {
            for (int r = 0; r < p; ++r) {
                const std::uint8_t factor = coefficients_[r][p];
                if (!factor) continue;
                gf256_muladd_mem(rhs_[r], factor, rhs_[p], static_cast<int>(block_bytes_));
                coefficients_[r][p] = 0;
            }
        }
        solved_ = true;
    }
    std::uint8_t* destination = static_cast<std::uint8_t*>(output);
    for (unsigned i = 0; i < kSourceCount; ++i)
        std::memcpy(destination + static_cast<std::size_t>(i) * block_bytes_, rhs_[i],
                    i == 5 ? tail_bytes_ : block_bytes_);
    return Answer(Status::Success, message_bytes_, message_bytes_);
}

} // namespace wh2_thue_native
