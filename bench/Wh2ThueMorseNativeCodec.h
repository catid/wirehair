#ifndef WH2_THUE_MORSE_NATIVE_CODEC_H
#define WH2_THUE_MORSE_NATIVE_CODEC_H

#include <cstddef>
#include <cstdint>
#include <memory>

// Benchmark-only K6 codec: not a public wirehair profile or API replacement.
// Call wirehair_init()/gf256_init() once before using this shared GF runtime.
namespace wh2_thue_native {

static const unsigned kSourceCount = 6;
static const std::size_t kLookupBytes = 39936;

// Exactly the packed 10/7/7/8 tables, immutable and alive through all handles.
// The six initial low-table rows must be the systematic identity. Other table
// contents are authenticated by the benchmark, not reconstructed by this core.
struct Lookup {
    const std::uint8_t* data;
    std::size_t bytes;
};

enum class Status {
    Success, NeedMore, InvalidInput, BufferTooSmall, Conflict, OutOfMemory
};

struct Result {
    Status status;
    std::size_t bytes_required;
    std::size_t bytes_written;
};

// Output is exactly six bytes and must not overlap the immutable lookup.
Status Row(Lookup lookup, std::uint32_t id, std::uint8_t output[kSourceCount]);

class Encoder {
public:
    // Output handle must be empty. All failures preserve it. Source is borrowed
    // and must stay readable/immutable until destruction, including for repairs.
    // Only a partial final block is privately copied and zero-padded.
    static Status Create(Lookup lookup, const void* source,
                         std::uint64_t message_bytes, std::uint32_t block_bytes,
                         std::unique_ptr<Encoder>& output);
    ~Encoder();

    // ID5 returns the meaningful tail; every other ID returns block_bytes.
    // Invalid/capacity failures do not write output. The written range must not
    // overlap this handle, lookup, retained source or private padding.
    Result Encode(std::uint32_t id, void* output, std::size_t capacity) const;

    Encoder(const Encoder&) = delete;
    Encoder& operator=(const Encoder&) = delete;

private:
    Encoder(Lookup lookup, const void* source, std::size_t message_bytes,
            std::uint32_t block_bytes, std::uint32_t tail_bytes);
    Lookup lookup_;
    const std::uint8_t* source_;
    std::size_t message_bytes_;
    std::uint32_t block_bytes_;
    std::uint32_t tail_bytes_;
    const void* sources_[kSourceCount];
    std::unique_ptr<std::uint8_t[]> padding_;
};

class Decoder {
public:
    static Status Create(Lookup lookup, std::uint64_t message_bytes,
                         std::uint32_t block_bytes,
                         std::unique_ptr<Decoder>& output);
    ~Decoder();

    // Exact ID5 tail length, full block otherwise. No receive allocation.
    // Rank6 is Success; actual back-substitution is deferred until Recover.
    // Matching dependent packets are idempotent. Any contradictory dependent
    // RHS returns Conflict with the accepted basis unchanged. There is no ID
    // ledger and no claim of public WH2's poisoned-state semantics.
    // Input may not overlap this handle or its owned slab.
    Result Feed(std::uint32_t id, const void* input, std::size_t bytes);

    // NeedMore/invalid/capacity failures preserve basis and output. Success
    // transforms the retained basis to identity, allowing subsequent Feed and
    // repeated Recover. The written range cannot overlap handle, slab or lookup.
    Result Recover(void* output, std::size_t capacity);
    unsigned Rank() const { return rank_; }

    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

private:
    Decoder(Lookup lookup, std::size_t message_bytes,
            std::uint32_t block_bytes, std::uint32_t tail_bytes);
    Lookup lookup_;
    std::size_t message_bytes_;
    std::uint32_t block_bytes_;
    std::uint32_t tail_bytes_;
    std::size_t slab_bytes_;
    std::uint8_t coefficients_[kSourceCount][kSourceCount];
    std::uint8_t pivot_mask_;
    unsigned rank_;
    bool solved_;
    std::unique_ptr<std::uint8_t[]> slab_;
    std::uint8_t* rhs_[kSourceCount];
    std::uint8_t* scratch_;
};

} // namespace wh2_thue_native
#endif
