#include "WirehairSmallFacade.h"

namespace {
alignas(64) const uint8_t kLookup[] = {
#include "WirehairSmallK3Lookup.inc"
};
static_assert(sizeof(kLookup) == wirehair_small_core::detail::Geometry<3>::LookupBytes,
              "Exact selected K3 lookup");
struct K3 {
    static wirehair_small_core::Lookup View()
    { return {kLookup, sizeof(kLookup)}; }
};
using F = wirehair_small_facade::Facade<3, K3>;
static_assert(F::ProfileId == WIREHAIR_SMALL_K3_PROFILE_ID, "K3 descriptor identity");
static_assert(F::MaxBlockBytes == WIREHAIR_SMALL_K3_MAX_BLOCK_BYTES, "K3 slab bound");
}

extern "C" WirehairSmallStatus wirehair_small_profile_validate(const void* p, size_t n) noexcept
{ return F::ProfileValidate(p, n); }
extern "C" WirehairSmallCreateResult wirehair_small_encoder_create(const void* s, uint64_t m, uint32_t b,
    uint32_t policy, void* p, size_t n) noexcept
{ return F::EncoderCreate(s, m, b, policy, p, n); }
extern "C" WirehairSmallCreateResult wirehair_small_encoder_create_profile(const void* s, const void* p,
    size_t n, uint32_t policy) noexcept
{ return F::EncoderCreateProfile(s, p, n, policy); }
extern "C" WirehairSmallCreateResult wirehair_small_decoder_create(const void* p, size_t n) noexcept
{ return F::DecoderCreate(p, n); }
extern "C" WirehairSmallStatus wirehair_small_encoder_detach_input(WirehairSmallCodec h) noexcept
{ return F::Detach(h); }
extern "C" WirehairSmallResult wirehair_small_encode(WirehairSmallCodec h, uint32_t id, void* out, size_t n) noexcept
{ return F::Encode(h, id, out, n); }
extern "C" WirehairSmallStatus wirehair_small_decode(WirehairSmallCodec h, uint32_t id, const void* in, size_t n) noexcept
{ return F::Decode(h, id, in, n); }
extern "C" WirehairSmallResult wirehair_small_recover(WirehairSmallCodec h, void* out, size_t n) noexcept
{ return F::Recover(h, out, n); }
extern "C" void wirehair_small_free(WirehairSmallCodec h) noexcept { F::Free(h); }
