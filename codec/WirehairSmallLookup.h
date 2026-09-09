#ifndef WIREHAIR_SMALL_LOOKUP_H
#define WIREHAIR_SMALL_LOOKUP_H

#include "WirehairSmallCore.h"

namespace wirehair_small_core {
// One immutable table shared by the WHK3 and WHV2 facades, never serialized.
Lookup K3Lookup();
}

#endif
