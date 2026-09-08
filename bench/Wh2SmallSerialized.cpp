#include "Wh2SmallSerializedCore.h"
#include "Wh2K3NativeData.inc"

namespace {
struct K3 {
    static wirehair_small_core::Lookup View()
    { return {wh2_k3_data::kLookup, sizeof(wh2_k3_data::kLookup)}; }
};
using F = wh2_small_serialized::Facade<3, K3>;
static_assert(F::ProfileId == WH2_SMALL_PROFILE_ID, "K3 descriptor identity");
static_assert(F::MaxBlockBytes == WH2_SMALL_MAX_BLOCK_BYTES, "K3 slab bound");
}

extern "C" Wh2SmallStatus wh2_small_profile_validate(const void* p, size_t n) noexcept
{ return F::ProfileValidate(p, n); }
extern "C" Wh2SmallCreateResult wh2_small_encoder_create(const void* s, uint64_t m, uint32_t b,
    uint32_t policy, void* p, size_t n) noexcept
{ return F::EncoderCreate(s, m, b, policy, p, n); }
extern "C" Wh2SmallCreateResult wh2_small_encoder_create_profile(const void* s, const void* p,
    size_t n, uint32_t policy) noexcept
{ return F::EncoderCreateProfile(s, p, n, policy); }
extern "C" Wh2SmallCreateResult wh2_small_decoder_create(const void* p, size_t n) noexcept
{ return F::DecoderCreate(p, n); }
extern "C" Wh2SmallStatus wh2_small_encoder_detach_input(Wh2SmallCodec h) noexcept
{ return F::Detach(h); }
extern "C" Wh2SmallResult wh2_small_encode(Wh2SmallCodec h, uint32_t id, void* out, size_t n) noexcept
{ return F::Encode(h, id, out, n); }
extern "C" Wh2SmallStatus wh2_small_decode(Wh2SmallCodec h, uint32_t id, const void* in, size_t n) noexcept
{ return F::Decode(h, id, in, n); }
extern "C" Wh2SmallResult wh2_small_recover(Wh2SmallCodec h, void* out, size_t n) noexcept
{ return F::Recover(h, out, n); }
extern "C" void wh2_small_free(Wh2SmallCodec h) noexcept { F::Free(h); }
