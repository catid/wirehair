#pragma once

#include "WirehairV2PrecodeEncode.h"

#include <wirehair/wirehair.h>

#include <stdint.h>

namespace wirehair_v2 {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)

/**
    Provisional, test-build-only identities for the two repair-v1 candidates.

    These are deliberately not public profile ids and are not returned by
    FindV2EquationContract().  A later authenticated campaign may publish at
    most one winner under a new public contract; until then dispatch-v1 is
    unchanged.
*/
enum class RepairV1Arm : uint32_t
{
    Pure8S0M1D3 = 0u,
    Pure9S0M1D3 = 1u
};

static const uint32_t kRepairV1AttemptCount = 8u;
static const uint32_t kRepairV1NoAttempt = UINT32_MAX;

struct RepairV1Contract
{
    uint64_t ProvisionalId;
    const char* Name;
    const char* CanonicalName;
    const char* CanonicalSha256;
    const char* RepairPolicySha256;
    uint32_t SeedAttemptCount;
    RepairV1Arm Arm;
    uint32_t GF256Rows;
    uint32_t GF16Rows;
    uint32_t DenseRows;
    uint32_t RecoveryMixCount;
};

/**
    Frozen repair-v1 selector policy identity.

    The digest is SHA-256 over kRepairV1PolicyCanonicalName with no trailing
    NUL or newline.
*/
extern const char kRepairV1PolicyCanonicalName[];
extern const char kRepairV1PolicySha256[];

const RepairV1Contract* RepairV1Contracts(uint32_t& count_out);
const RepairV1Contract* FindRepairV1Contract(RepairV1Arm arm);
const RepairV1Contract* FindRepairV1Contract(uint64_t provisional_id);
const RepairV1Contract* FindRepairV1Contract(const char* name);

struct RepairV1AttemptTelemetry
{
    bool ProbeExecuted = false;
    WirehairResult ProbeResult = Wirehair_InvalidInput;
    bool ProbeStatsAvailable = false;
    PrecodeSolveStats ProbeStats = {};
    bool RealExecuted = false;
    WirehairResult RealResult = Wirehair_InvalidInput;
    bool RealStatsAvailable = false;
    PrecodeSolveStats RealStats = {};
};

/**
    Complete work receipt for one lazy selector invocation.

    AttemptsExecuted is the number of distinct attempt indices considered.
    CallsExecuted counts every real-payload configuration call or structural
    probe call, including a call that exits before algebra on validation/OOM.
    A failed call never commits the destination endpoint; Committed is true
    only after the final successful real-payload solve.
*/
struct RepairV1SelectorTelemetry
{
    RepairV1AttemptTelemetry Attempts[kRepairV1AttemptCount] = {};
    uint32_t AttemptsExecuted = 0u;
    uint32_t CallsExecuted = 0u;
    uint32_t RealConfigurationCalls = 0u;
    uint32_t StructuralProbeCalls = 0u;
    uint32_t SelectedAttempt = kRepairV1NoAttempt;
    WirehairResult Result = Wirehair_InvalidInput;
    bool CapExhausted = false;
    bool FatalAttemptZeroMismatch = false;
    bool Oom = false;
    bool Committed = false;
};

/**
    Bounded per-thread geometry transaction for a provisional contract.

    Construction first verifies that every equation-active control other than
    the arm's row split is canonical, then installs the arm's 8+0 or 9+0
    mixed-completion geometry on this thread.  The complete prior row state is
    restored on destruction.  No process-global hook is changed.
*/
class ScopedRepairV1ContractStateForTesting
{
public:
    explicit ScopedRepairV1ContractStateForTesting(
        const RepairV1Contract& contract);
    explicit ScopedRepairV1ContractStateForTesting(
        uint64_t provisional_id);
    ~ScopedRepairV1ContractStateForTesting() noexcept;

    ScopedRepairV1ContractStateForTesting(
        const ScopedRepairV1ContractStateForTesting&) = delete;
    ScopedRepairV1ContractStateForTesting& operator=(
        const ScopedRepairV1ContractStateForTesting&) = delete;

    bool IsValid() const noexcept { return Valid; }

private:
    void Initialize(const RepairV1Contract* contract);

    uint32_t SavedPeriod = 0u;
    uint32_t SavedGF256Rows = 0u;
    uint32_t SavedGF16Rows = 0u;
    bool HaveSnapshot = false;
    bool GeometryChanged = false;
    bool Valid = false;
    ScopedMixedBandTrackingXForTesting TrackingXSnapshot;
};

/**
    Expand one exact arm/root/attempt descriptor.

    Root is a uint32_t by design: its zero extension is the attempt-zero
    precode seed and its exact value is the attempt-zero packet seed.  Attempts
    outside [0, 8) fail without modifying config_out.  The caller must hold a
    valid ScopedRepairV1ContractStateForTesting for `contract`; the returned
    explicit descriptor is pinned to that bounded state and cannot be used
    after the scope without opening an equivalent contract scope.
*/
bool MakeRepairV1ExplicitConfigForTesting(
    const RepairV1Contract& contract,
    uint32_t block_count,
    uint32_t construction_root,
    uint32_t seed_attempt,
    bool cache_systematic_source,
    bool cache_received_systematic,
    ExplicitMessagePrecodeConfigForTesting& config_out);

/**
    Deterministic policy-machine oracle used to exhaustively test lazy result
    handling without sampling any campaign outcomes.

    ProbeResults are read only until the first success or fatal result.
    SelectedRealResult is consumed only when a probe selects attempt > 0.
*/
WirehairResult EvaluateRepairV1PolicyOracleForTesting(
    WirehairResult attempt_zero_real_result,
    const WirehairResult probe_results[kRepairV1AttemptCount],
    WirehairResult selected_real_result,
    RepairV1SelectorTelemetry& telemetry_out);

enum class RepairV1AccountingOomStageForTesting : uint32_t
{
    AttemptZeroBeforeSystemBuild = 0u,
    ProbeDuringSystemBuild = 1u,
    ProbeDuringSolve = 2u
};

/**
    Exercise begin/finish bookkeeping for pre-algebra and partial-algebra OOM
    without sampling a codec construction seed.
*/
WirehairResult EvaluateRepairV1AccountingOomForTesting(
    RepairV1AccountingOomStageForTesting stage,
    RepairV1SelectorTelemetry& telemetry_out);

#endif // WIREHAIR_V2_ENABLE_TEST_HOOKS

} // namespace wirehair_v2
