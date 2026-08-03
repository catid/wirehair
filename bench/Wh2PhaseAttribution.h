#ifndef WIREHAIR_BENCH_WH2_PHASE_ATTRIBUTION_H
#define WIREHAIR_BENCH_WH2_PHASE_ATTRIBUTION_H

#include "Wh2NativeCodec.h"
#include "Wh2NativePanel.h"

#include <array>
#include <memory>
#include <stdint.h>
#include <string>
#include <vector>

namespace wirehair_wh2_bench {

/** Fixed sides of the diagnostic Two07-versus-head solve panel. */
enum class PhaseSolveArm : uint32_t
{
    Two07 = 0,
    Head
};

/**
    One chronological solve observation captured outside the generic panel.

    ElapsedNanoseconds is the outer cold-solve interval reported by
    NativeSolveFixture.  The five component intervals live in Stats.  Failed
    solves have a null-equivalent zero outer interval and are retained only as
    a weak-cell ledger; they never contribute phase timing evidence.
*/
struct PhaseSolveObservation
{
    PhaseSolveArm Arm = PhaseSolveArm::Two07;
    WirehairResult Result = Wirehair_Error;
    bool BytesVerified = false;
    uint64_t ElapsedNanoseconds = 0u;
    wirehair_v2::PrecodeSolveStats Stats = {};
};

/** Shared append-only chronology for fresh panel invocations. */
class PhaseObservationCollector
{
public:
    void Append(const PhaseSolveObservation& observation);
    const std::vector<PhaseSolveObservation>& Observations() const;

private:
    std::vector<PhaseSolveObservation> Values;
};

/**
    Fresh decoder-solve invocation used only by the phase diagnostic.

    The worker prepares one immutable fixture per arm outside the panel.
    Every generic-panel factory still creates a fresh invocation object, while
    Invoke() calls the shared fixture's const Run(), which creates fresh solve
    state, and appends its result only after the outer timer has stopped.  A
    stable preparation failure may be supplied without a fixture for the weak
    ledger.  The production solver and generic panel remain unchanged.
*/
class PhaseSolveInvocation : public NativePanelInvocation
{
public:
    PhaseSolveInvocation(
        const std::string& identity,
        PhaseSolveArm arm,
        WirehairResult preparation_result,
        const std::shared_ptr<const NativeSolveFixture>& fixture,
        uint32_t fixed_received_overhead,
        const std::shared_ptr<PhaseObservationCollector>& collector);

    std::string Identity() const override;
    NativePanelInvocationResult Invoke() override;

private:
    std::string IdentityValue;
    PhaseSolveArm ArmValue;
    uint32_t FixedReceivedOverhead;
    WirehairResult PreparationResult;
    std::shared_ptr<PhaseObservationCollector> Collector;
    std::shared_ptr<const NativeSolveFixture> Fixture;
};

/** Chronology tag reconstructed independently from the generic panel. */
struct PhaseMeasuredObservation
{
    uint32_t Block = 0u;
    uint32_t Repeat = 0u;
    uint32_t Slot = 0u;
    PhaseSolveObservation Observation;
};

/** Overflow-safe sums for one of the eight counterbalanced slots. */
struct PhaseSlotTotals
{
    bool HasElapsed = false;
    uint64_t OuterNanoseconds = 0u;
    uint64_t BuildNanoseconds = 0u;
    uint64_t PeelNanoseconds = 0u;
    uint64_t ProjectNanoseconds = 0u;
    uint64_t ResidualNanoseconds = 0u;
    uint64_t BackSubNanoseconds = 0u;
};

struct PhasePanelAssembly
{
    bool Comparable = false;
    PhaseSolveObservation LeftWarmup;
    PhaseSolveObservation RightWarmup;
    std::vector<PhaseMeasuredObservation> Measured;
    std::array<PhaseSlotTotals, 8> Slots;
    /** Complete per-arm non-timing counters; timer fields are cleared. */
    wirehair_v2::PrecodeSolveStats LeftCounters = {};
    wirehair_v2::PrecodeSolveStats RightCounters = {};
};

/**
    Authenticate and independently reconstruct a completed phase panel.

    The fixed semantic sides are left=Two07 and right=production head.  The
    function requires exactly two warmups followed by 4*n observations in the
    generic scheduler's documented repeat-major chronology, validates every
    outcome/byte/counter invariant, proves each five-phase sum fits inside its
    outer interval, reconstructs the exact eight slot sums, and checks those
    sums against the generic panel.  Weak cells preserve only stable outcomes
    and null slot timings.  result_out is unchanged on failure.
*/
bool ValidateAndAssemblePhasePanel(
    const NativePanelResult& panel,
    NativePanelOrder expected_order,
    uint32_t invocations_per_slot,
    const std::vector<PhaseSolveObservation>& observations,
    PhasePanelAssembly& result_out,
    std::string& diagnostic);

} // namespace wirehair_wh2_bench

#endif // WIREHAIR_BENCH_WH2_PHASE_ATTRIBUTION_H
