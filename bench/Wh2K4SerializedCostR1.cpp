#error "K4 R1 timing family retired: preserve spent evidence; no claim-path relaunch"
#if 0 // Historical wrapper; original source remains at each producing commit.
// Fixed K4 configuration for the fresh R1 prospective lifecycle gate.
#define WH2_SMALL_COST_K 4
#define WH2_SMALL_COST_REPAIRS 8
#define WH2_SMALL_COST_PROTOCOL "wirehair.wh2.k4-serialized-cost-r0"
#define WH2_SMALL_COST_CLAIM_PATH "/var/tmp/wh2-k4-serialized-cost-r1.R24/science/CLAIM.json"
#include "Wh2SmallLifecycleWorkerR0.h"
#endif
