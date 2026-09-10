#ifndef WH2_SMALL_SERIALIZED_H
#define WH2_SMALL_SERIALIZED_H

#include <stddef.h>
#include <stdint.h>

/* Private benchmark C boundary, built for one of K2, K3 (default), K4, K5 or K8. Not
 * installed, not selected by WH1/WH2/K6, and not a production speed claim.
 * Initialize the shared runtime with wirehair_init() before codec calls.
 * Only live handles returned by this boundary may be passed back. Calls and
 * source mutation/free are externally serialized. A null free is harmless.
 *
 * The 32-byte little-endian descriptor fixes the exact sealed equations:
 * WHK2/WHK3/WHK4/WHK5/WHK8 / u16 version1 / u16 size32 / u64 profile ID / u64 message bytes /
 * u32 block bytes / u32 zero. No external table or live encoder is required.
 * Unknown and retired descriptors are rejected, never reinterpreted.
 */
#define WH2_SMALL_PROFILE_BYTES 32u
#ifndef WH2_SMALL_CODEC_K
#define WH2_SMALL_CODEC_K 3
#endif
#if WH2_SMALL_CODEC_K != 2 && WH2_SMALL_CODEC_K != 3 && WH2_SMALL_CODEC_K != 4 && WH2_SMALL_CODEC_K != 5 && WH2_SMALL_CODEC_K != 8
#error "Unsupported benchmark boundary dimension"
#endif
#define WH2_SMALL_PROFILE_ID (UINT64_C(0x5748324b30544d31) + ((uint64_t)WH2_SMALL_CODEC_K << 24))
#define WH2_SMALL_MAX_BLOCK_BYTES (UINT32_C(268435456) / (WH2_SMALL_CODEC_K + 1u))

#ifdef __cplusplus
#define WH2_SMALL_NOEXCEPT noexcept
extern "C" {
#else
#define WH2_SMALL_NOEXCEPT
#endif

typedef void* Wh2SmallCodec;
typedef enum Wh2SmallStatus {
    Wh2Small_Success = 0, Wh2Small_NeedMore, Wh2Small_InvalidInput,
    Wh2Small_BufferTooSmall, Wh2Small_UnsupportedProfile, Wh2Small_InvalidDimensions,
    Wh2Small_Conflict, Wh2Small_OutOfMemory
} Wh2SmallStatus;
typedef enum Wh2SmallSourcePolicy {
    Wh2Small_Independent = 1, Wh2Small_BorrowedImmutable = 2
} Wh2SmallSourcePolicy;
typedef struct Wh2SmallCreateResult {
    Wh2SmallStatus status;
    Wh2SmallCodec codec; /* Null on failure. */
} Wh2SmallCreateResult;
typedef struct Wh2SmallResult {
    Wh2SmallStatus status;
    uint64_t bytes_required, bytes_written; /* Written is zero on failure. */
} Wh2SmallResult;

Wh2SmallStatus wh2_small_profile_validate(const void*, size_t) WH2_SMALL_NOEXCEPT;
/* Independent copies before descriptor output, which may alias the original
 * source. Borrowed requires every source byte to stay readable and immutable
 * through successful detach/free: repair packets also read source bytes.
 * Failed creation preserves descriptor output and returns no handle. */
Wh2SmallCreateResult wh2_small_encoder_create(const void*, uint64_t, uint32_t,
    uint32_t, void*, size_t) WH2_SMALL_NOEXCEPT;
Wh2SmallCreateResult wh2_small_encoder_create_profile(const void*, const void*,
    size_t, uint32_t) WH2_SMALL_NOEXCEPT;
Wh2SmallCreateResult wh2_small_decoder_create(const void*, size_t) WH2_SMALL_NOEXCEPT;
/* First borrowed detach allocates/copies transactionally; failure preserves
 * the old source obligation. Independent/repeated detach allocates nothing. */
Wh2SmallStatus wh2_small_encoder_detach_input(Wh2SmallCodec) WH2_SMALL_NOEXCEPT;
/* ID K-1 carries its meaningful tail; other packets have block_bytes bytes.
 * Outputs may not overlap handle/private storage/table or retained source.
 * Encode/Decode/Recover allocate nothing. Capacity/alias failures do not write.
 * Decode copies input during the call. Contradictions permanently poison the
 * decoder; matching dependent packets and repeated recovery are idempotent.
 * Recovery is not authentication: use trusted framing and a digest/MAC. */
Wh2SmallResult wh2_small_encode(Wh2SmallCodec, uint32_t, void*, size_t) WH2_SMALL_NOEXCEPT;
Wh2SmallStatus wh2_small_decode(Wh2SmallCodec, uint32_t, const void*, size_t) WH2_SMALL_NOEXCEPT;
Wh2SmallResult wh2_small_recover(Wh2SmallCodec, void*, size_t) WH2_SMALL_NOEXCEPT;
void wh2_small_free(Wh2SmallCodec) WH2_SMALL_NOEXCEPT;

#ifdef __cplusplus
}
#endif
#undef WH2_SMALL_NOEXCEPT
#endif
