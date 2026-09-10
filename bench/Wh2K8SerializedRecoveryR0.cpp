// Frozen retained K8 comparison; no new equation/trace selection or timing.
#include "Wh2K8NativeData.inc"
namespace small_recovery_data = wh2_k8_data;
namespace small_recovery_config {
constexpr unsigned K=8, lambda=2, history_count=44, row_count=2347;
constexpr unsigned records=6260, packet_ids=74946;
constexpr unsigned char feedback[8]={96,19,186,153,85,252,7,255};
constexpr char protocol[]="wirehair.wh2.k8-serialized-recovery-r0";
constexpr char claim_path[]="/var/tmp/wh2-k8-serialized-recovery-r0/CLAIM.json";
}
#include "Wh2SmallRecoveryWorkerR0.h"
