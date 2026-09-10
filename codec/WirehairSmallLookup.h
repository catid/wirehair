#ifndef WIREHAIR_SMALL_LOOKUP_H
#define WIREHAIR_SMALL_LOOKUP_H

#include "WirehairSmallCore.h"

namespace wirehair_small_core {
// One immutable table shared by the WHK3 and WHV2 facades, never serialized.
Lookup K3Lookup();
// Sealed K5 table, in a separate translation unit from the existing K3 facade.
Lookup K5Lookup();
// Sealed lambda-2 K8 table; no runtime selection or dependency on bench data.
Lookup K8Lookup();
}

#endif
