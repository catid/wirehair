#pragma once

#include "Wh2RdpruTargetIdentityV2.h"

#include <cstdint>
#include <string>

namespace wirehair_wh2_bench {

// Capture the already-pinned logical CPU using the frozen target-identity v2
// semantic machinery without applying its separate CPU50 equality manifest.
bool CapturePublicBorrowedTargetIdentity(
    int32_t target_cpu,
    TargetIdentityReceiptV2& receipt,
    std::string& diagnostic);

} // namespace wirehair_wh2_bench
