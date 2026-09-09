#include "WirehairSmallLookup.h"

namespace {
alignas(64) const uint8_t kLookup[] = {
#include "WirehairSmallK5Lookup.inc"
};
static_assert(sizeof(kLookup) == wirehair_small_core::detail::Geometry<5>::LookupBytes,
              "Exact selected K5 lookup");
}

wirehair_small_core::Lookup wirehair_small_core::K5Lookup()
{ return {kLookup, sizeof(kLookup)}; }
