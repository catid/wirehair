#ifndef WIREHAIR_SMALL_SERIALIZED_H
#define WIREHAIR_SMALL_SERIALIZED_H

#include "wirehair.h"
#include <stddef.h>
#include <stdint.h>

/* Explicit opt-in small-block GF(256) API. Currently supports exactly K=3;
 * other dimensions/profiles are rejected. Never selected by existing
 * WH1/WH2/K6 APIs; their wire identities and source contracts are unchanged.
 * Qualification is scoped: see SMALL_WIRE_PROFILES.md before deployment.
 * Initialize the shared runtime once with wirehair_init(). Calls on a handle
 * and source mutation/free must be externally serialized. Only live handles
 * returned by this API may be passed back; null free is harmless.
 *
 * K3 descriptor: 32 bytes, little endian. Magic WHK3 at 0; u16 version=1 at 4;
 * u16 size=32 at 6; u64 profile ID at 8; u64 message bytes at 16; u32 block
 * bytes at 24; u32 reserved=0 at 28. Require 2*block < message <= 3*block.
 * The ID fixes the pair (8,14,7)/(9,14,7) and the immutable 13056-byte lookup.
 * No caller-supplied table, seed search, live encoder, or process-local profile
 * registry is needed to create a decoder. Unknown and retired descriptors
 * are rejected, never reinterpreted. The ID is not an integrity check.
 */
#define WIREHAIR_SMALL_PROFILE_BYTES 32u
#define WIREHAIR_SMALL_K3_PROFILE_ID UINT64_C(0x5748324b33544d31)
#define WIREHAIR_SMALL_K3_MAX_BLOCK_BYTES (UINT32_C(268435456) / 4u)

#ifdef __cplusplus
#define WIREHAIR_SMALL_NOEXCEPT noexcept
extern "C" {
#else
#define WIREHAIR_SMALL_NOEXCEPT
#endif

typedef struct WirehairSmallCodecImpl* WirehairSmallCodec;
typedef enum WirehairSmallStatus {
    WirehairSmall_Success = 0, WirehairSmall_NeedMore, WirehairSmall_InvalidInput,
    WirehairSmall_BufferTooSmall, WirehairSmall_UnsupportedProfile, WirehairSmall_InvalidDimensions,
    WirehairSmall_Conflict, WirehairSmall_OutOfMemory
} WirehairSmallStatus;
typedef enum WirehairSmallSourcePolicy {
    WirehairSmall_Independent = 1, WirehairSmall_BorrowedImmutable = 2
} WirehairSmallSourcePolicy;
typedef struct WirehairSmallCreateResult {
    WirehairSmallStatus status;
    WirehairSmallCodec codec; /* Null on failure. */
} WirehairSmallCreateResult;
typedef struct WirehairSmallResult {
    WirehairSmallStatus status;
    uint64_t bytes_required;
    uint64_t bytes_written; /* Zero on failure; buffers remain unchanged. */
} WirehairSmallResult;

WIREHAIR_EXPORT WirehairSmallStatus wirehair_small_profile_validate(const void* profile,
    size_t bytes) WIREHAIR_SMALL_NOEXCEPT;

/* Independent copies the source before returning. BorrowedImmutable requires
 * every source byte to remain readable and immutable until successful detach
 * or free: repair packets ALSO read the source, unlike existing WH2.
 * Profile output is written only on success; it may overlap source only in
 * Independent mode. A short profile capacity performs no allocation.
 */
WIREHAIR_EXPORT WirehairSmallCreateResult wirehair_small_encoder_create(const void* source,
    uint64_t message_bytes, uint32_t block_bytes, uint32_t source_policy,
    void* profile, size_t profile_capacity) WIREHAIR_SMALL_NOEXCEPT;
WIREHAIR_EXPORT WirehairSmallCreateResult wirehair_small_encoder_create_profile(const void* source,
    const void* profile, size_t profile_bytes, uint32_t source_policy) WIREHAIR_SMALL_NOEXCEPT;
WIREHAIR_EXPORT WirehairSmallCreateResult wirehair_small_decoder_create(const void* profile,
    size_t profile_bytes) WIREHAIR_SMALL_NOEXCEPT;

/* First borrowed detach allocates/copies transactionally. Failure preserves
 * the old encoder and source obligation. Success permits source mutation/free.
 * Independent and already detached encoders detach without allocation.
 */
WIREHAIR_EXPORT WirehairSmallStatus wirehair_small_encoder_detach_input(WirehairSmallCodec codec)
    WIREHAIR_SMALL_NOEXCEPT;

/* K3 ID2 has the meaningful tail length; every other ID has block length.
 * Writable ranges may not overlap the handle, private storage, immutable
 * lookup or retained source. Capacity and pointer-wrap checks precede writes.
 * Encode, decode and recover allocate no memory.
 */
WIREHAIR_EXPORT WirehairSmallResult wirehair_small_encode(WirehairSmallCodec codec,
    uint32_t id, void* output, size_t capacity) WIREHAIR_SMALL_NOEXCEPT;

/* Feed consumes bytes during the call and retains no input pointer. Input
 * cannot overlap the handle or its private storage. Matching dependent packets
 * are idempotent. Contradictions permanently poison the decoder: subsequent
 * valid feed/recover calls report Conflict. Invalid input and capacity failures
 * do not poison it. NeedMore preserves the accepted equations and recovery can
 * be repeated. Recovery is NOT authentication: use trusted framing/digest/MAC.
 */
WIREHAIR_EXPORT WirehairSmallStatus wirehair_small_decode(WirehairSmallCodec codec,
    uint32_t id, const void* input, size_t bytes) WIREHAIR_SMALL_NOEXCEPT;
WIREHAIR_EXPORT WirehairSmallResult wirehair_small_recover(WirehairSmallCodec codec,
    void* output, size_t capacity) WIREHAIR_SMALL_NOEXCEPT;
WIREHAIR_EXPORT void wirehair_small_free(WirehairSmallCodec codec) WIREHAIR_SMALL_NOEXCEPT;

#ifdef __cplusplus
}
#endif
#undef WIREHAIR_SMALL_NOEXCEPT
#endif
