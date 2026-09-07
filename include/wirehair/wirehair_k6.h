#ifndef WIREHAIR_K6_SERIALIZED_H
#define WIREHAIR_K6_SERIALIZED_H

#include "wirehair.h"
#include <stddef.h>
#include <stdint.h>

/* Explicit opt-in K6 GF(256) codec. Never selected by existing WH1/WH2 APIs.
 * This has a distinct wire identity and source-ownership contract.
 * Qualification is scoped; see K6_WIRE_PROFILE.md before deployment.
 * Initialize the shared GF runtime once with wirehair_init().
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
#define WIREHAIR_K6_PROFILE_BYTES 32u
#define WIREHAIR_K6_PROFILE_ID UINT64_C(0x5748324b36544d31)
#define WIREHAIR_K6_MAX_BLOCK_BYTES (UINT32_C(268435456) / 7u)

#ifdef __cplusplus
#define WIREHAIR_K6_NOEXCEPT noexcept
extern "C" {
#else
#define WIREHAIR_K6_NOEXCEPT
#endif

typedef struct WirehairK6CodecImpl* WirehairK6Codec;
typedef enum WirehairK6Status {
    WirehairK6_Success = 0, WirehairK6_NeedMore, WirehairK6_InvalidInput,
    WirehairK6_BufferTooSmall, WirehairK6_UnsupportedProfile, WirehairK6_InvalidDimensions,
    WirehairK6_Conflict, WirehairK6_OutOfMemory
} WirehairK6Status;
typedef enum WirehairK6SourcePolicy {
    WirehairK6_Independent = 1, WirehairK6_BorrowedImmutable = 2
} WirehairK6SourcePolicy;
typedef struct WirehairK6CreateResult {
    WirehairK6Status status;
    WirehairK6Codec codec; /* Always null on failure. */
} WirehairK6CreateResult;
typedef struct WirehairK6Result {
    WirehairK6Status status;
    uint64_t bytes_required;
    uint64_t bytes_written; /* Zero on failure; buffers remain unchanged. */
} WirehairK6Result;

WIREHAIR_EXPORT WirehairK6Status wirehair_k6_profile_validate(const void* profile, size_t bytes) WIREHAIR_K6_NOEXCEPT;

/* Independent copies the source. BorrowedImmutable requires every message
 * byte to remain readable and immutable through successful detach or free;
 * unlike public WH2, repair packets ALSO read the borrowed source.
 * Profile output is written only on success. It may overlap source only in
 * Independent mode. Short profile capacity performs no allocation.
 */
WIREHAIR_EXPORT WirehairK6CreateResult wirehair_k6_encoder_create(const void* source, uint64_t message_bytes,
    uint32_t block_bytes, uint32_t source_policy, void* profile,
    size_t profile_capacity) WIREHAIR_K6_NOEXCEPT;
WIREHAIR_EXPORT WirehairK6CreateResult wirehair_k6_encoder_create_profile(const void* source,
    const void* profile, size_t profile_bytes, uint32_t source_policy) WIREHAIR_K6_NOEXCEPT;
WIREHAIR_EXPORT WirehairK6CreateResult wirehair_k6_decoder_create(const void* profile,
    size_t profile_bytes) WIREHAIR_K6_NOEXCEPT;

/* Detach is idempotent. Unlike public WH2 detach, the first borrowed detach
 * allocates/copies. Failure is transactional: the old encoder and its source
 * obligation remain intact. Success permits immediate source mutation/free.
 */
WIREHAIR_EXPORT WirehairK6Status wirehair_k6_encoder_detach_input(WirehairK6Codec codec) WIREHAIR_K6_NOEXCEPT;

/* ID5 has the exact meaningful tail length; every other ID has block length.
 * Writable ranges must not overlap the handle, its private storage, immutable
 * table, or retained source. Capacity and pointer-wrap checks precede writing.
 */
WIREHAIR_EXPORT WirehairK6Result wirehair_k6_encode(WirehairK6Codec codec, uint32_t id, void* output,
    size_t capacity) WIREHAIR_K6_NOEXCEPT;

/* Input bytes are consumed during the call; no input pointer is retained.
 * Input cannot overlap the handle or its private storage. Duplicate/dependent
 * matching packets are idempotent. Contradictions permanently poison this
 * decoder: subsequent valid feed/recover calls report Conflict. Invalid input
 * and capacity failures do not poison it. Feed and recover allocate nothing.
 */
WIREHAIR_EXPORT WirehairK6Status wirehair_k6_decode(WirehairK6Codec codec, uint32_t id, const void* input,
    size_t bytes) WIREHAIR_K6_NOEXCEPT;
WIREHAIR_EXPORT WirehairK6Result wirehair_k6_recover(WirehairK6Codec codec, void* output,
    size_t capacity) WIREHAIR_K6_NOEXCEPT;
WIREHAIR_EXPORT void wirehair_k6_free(WirehairK6Codec codec) WIREHAIR_K6_NOEXCEPT;

/* Recovery is NOT authentication: verify trusted framing/digest/MAC. */
#ifdef __cplusplus
}
#endif
#undef WIREHAIR_K6_NOEXCEPT
#endif
