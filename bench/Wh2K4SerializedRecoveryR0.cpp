// Frozen retained K4 comparison; no new equation/trace selection or timing.
#include "Wh2K4NativeData.inc"
namespace small_recovery_data = wh2_k4_data;
namespace small_recovery_config {
constexpr unsigned K=4, lambda=1, history_count=38, row_count=2252;
constexpr unsigned records=6254, packet_ids=49891;
constexpr unsigned char feedback[4]={64,120,54,15};
constexpr char protocol[]="wirehair.wh2.k4-serialized-recovery-r0";
constexpr char claim_path[]="/var/tmp/wh2-k4-serialized-recovery-r0/CLAIM.json";
}
#include "Wh2SmallRecoveryWorkerR0.h"
