#include "WirehairSmallLookup.h"

namespace {
alignas(64) const uint8_t kLookup[] = {
#include "WirehairSmallK8Lookup.inc"
};
static_assert(sizeof(kLookup) == wirehair_small_core::detail::Geometry<8>::LookupBytes,
              "Exact selected K8 lookup");
}

wirehair_small_core::Lookup wirehair_small_core::K8Lookup()
{ return {kLookup, sizeof(kLookup)}; }
