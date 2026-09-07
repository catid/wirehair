#ifndef WH2_K6_SERIALIZED_H
#define WH2_K6_SERIALIZED_H

#include <stddef.h>
#include <stdint.h>

/* Experimental K6 GF(256) codec, NOT a production Wirehair profile or ABI.
 * Initialize the shared GF runtime once with wirehair_init()/gf256_init().
 * All calls on a handle (including detach/free) are externally serialized.
 * Only live handles returned here may be passed back; null free is accepted.
 * No C++ owning object or caller-supplied equation table crosses this boundary.
 *
 * Descriptor: exactly 32 bytes, little endian. Magic "WHK6" at 0; version
 * u16=1 at 4; size u16=32 at 6; profile ID u64 at 8; message u64 at 16;
 * block u32 at 24; reserved u32=0 at 28. K=ceil(message/block) must equal 6.
 * The ID fixes the selected Thue-Morse pair AND the immutable 39936-byte table
 * (SHA256 27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d).
 * No seed search, external lookup, live encoder, or local profile state is
 * required to reconstruct a decoder. Unknown/current/retired WH2 descriptors
 * are rejected, never reinterpreted. The ID is not an integrity check.
 */
#define WH2_K6_PROFILE_BYTES 32u
#define WH2_K6_PROFILE_ID UINT64_C(0x5748324b36544d31)
#define WH2_K6_MAX_BLOCK_BYTES (UINT32_C(268435456) / 7u)

#ifdef __cplusplus
#define WH2_K6_NOEXCEPT noexcept
extern "C" {
#else
#define WH2_K6_NOEXCEPT
#endif

typedef struct Wh2K6CodecImpl* Wh2K6Codec;
typedef enum Wh2K6Status {
    Wh2K6_Success = 0, Wh2K6_NeedMore, Wh2K6_InvalidInput,
    Wh2K6_BufferTooSmall, Wh2K6_UnsupportedProfile, Wh2K6_InvalidDimensions,
    Wh2K6_Conflict, Wh2K6_OutOfMemory
} Wh2K6Status;
typedef enum Wh2K6SourcePolicy {
    Wh2K6_Independent = 1, Wh2K6_BorrowedImmutable = 2
} Wh2K6SourcePolicy;
typedef struct Wh2K6CreateResult {
    Wh2K6Status status;
    Wh2K6Codec codec; /* Always null on failure. */
} Wh2K6CreateResult;
typedef struct Wh2K6Result {
    Wh2K6Status status;
    uint64_t bytes_required;
    uint64_t bytes_written; /* Zero on failure; buffers remain unchanged. */
} Wh2K6Result;

Wh2K6Status wh2_k6_profile_validate(const void* profile, size_t bytes) WH2_K6_NOEXCEPT;

/* Independent copies the source. BorrowedImmutable requires every message
 * byte to remain readable and immutable through successful detach or free;
 * unlike public WH2, repair packets ALSO read the borrowed source.
 * Profile output is written only on success. It may overlap source only in
 * Independent mode. Short profile capacity performs no allocation.
 */
Wh2K6CreateResult wh2_k6_encoder_create(const void* source, uint64_t message_bytes,
    uint32_t block_bytes, uint32_t source_policy, void* profile,
    size_t profile_capacity) WH2_K6_NOEXCEPT;
Wh2K6CreateResult wh2_k6_encoder_create_profile(const void* source,
    const void* profile, size_t profile_bytes, uint32_t source_policy) WH2_K6_NOEXCEPT;
Wh2K6CreateResult wh2_k6_decoder_create(const void* profile,
    size_t profile_bytes) WH2_K6_NOEXCEPT;

/* Detach is idempotent. Unlike public WH2 detach, the first borrowed detach
 * allocates/copies. Failure is transactional: the old encoder and its source
 * obligation remain intact. Success permits immediate source mutation/free.
 */
Wh2K6Status wh2_k6_encoder_detach_input(Wh2K6Codec codec) WH2_K6_NOEXCEPT;

/* ID5 has the exact meaningful tail length; every other ID has block length.
 * Writable ranges must not overlap the handle, its private storage, immutable
 * table, or retained source. Capacity and pointer-wrap checks precede writing.
 */
Wh2K6Result wh2_k6_encode(Wh2K6Codec codec, uint32_t id, void* output,
    size_t capacity) WH2_K6_NOEXCEPT;

/* Input bytes are consumed during the call; no input pointer is retained.
 * Input cannot overlap the handle or its private storage. Duplicate/dependent
 * matching packets are idempotent. Contradictions permanently poison this
 * decoder: subsequent valid feed/recover calls report Conflict. Invalid input
 * and capacity failures do not poison it. Feed and recover allocate nothing.
 */
Wh2K6Status wh2_k6_decode(Wh2K6Codec codec, uint32_t id, const void* input,
    size_t bytes) WH2_K6_NOEXCEPT;
Wh2K6Result wh2_k6_recover(Wh2K6Codec codec, void* output,
    size_t capacity) WH2_K6_NOEXCEPT;
void wh2_k6_free(Wh2K6Codec codec) WH2_K6_NOEXCEPT;

/* Recovery is NOT authentication: verify trusted framing/digest/MAC. */
#ifdef __cplusplus
}
#endif
#undef WH2_K6_NOEXCEPT
#endif
