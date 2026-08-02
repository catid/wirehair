#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "Wh2NativePanel.h"

#include <errno.h>
#include <limits>
#include <sstream>
#include <string.h>

#if defined(__linux__)
#include <sched.h>
#endif

namespace wirehair_wh2_bench {
namespace {

static const uint64_t kMaxInt63 =
    (uint64_t)std::numeric_limits<int64_t>::max();

bool IsKnownOrder(NativePanelOrder order)
{
    return order == NativePanelOrder::ABBA ||
        order == NativePanelOrder::BAAB;
}

bool IsKnownDisposition(NativePanelDisposition disposition)
{
    return disposition == NativePanelDisposition::Success ||
        disposition == NativePanelDisposition::PreflightFailure ||
        disposition == NativePanelDisposition::Fatal;
}

bool SameOutcome(
    const NativePanelInvocationResult& expected,
    const NativePanelInvocationResult& actual)
{
    if (expected.Disposition != actual.Disposition ||
        expected.OutcomeCode != actual.OutcomeCode ||
        expected.HasDecodedExtra != actual.HasDecodedExtra)
    {
        return false;
    }
    return !expected.HasDecodedExtra ||
        expected.DecodedExtra == actual.DecodedExtra;
}

void ClearTimings(NativePanelResult& result)
{
    result.PanelComparable = false;
    for (size_t i = 0; i < result.Slots.size(); ++i)
    {
        result.Slots[i].HasElapsedNanoseconds = false;
        result.Slots[i].ElapsedNanoseconds = 0u;
        result.Slots[i].Invocation.ElapsedNanoseconds = 0u;
    }
}

NativePanelResult Fail(
    NativePanelResult result,
    NativePanelStatus status,
    const std::string& diagnostic)
{
    result.Status = status;
    result.Diagnostic = diagnostic;
    ClearTimings(result);
    return result;
}

const NativePanelArm& SelectArm(
    NativePanelSide side,
    const NativePanelArm& left,
    const NativePanelArm& right)
{
    return side == NativePanelSide::Left ? left : right;
}

struct InvocationObservation
{
    NativePanelStatus Status;
    std::string Diagnostic;
    NativePanelInvocationResult Result;

    InvocationObservation()
        : Status(NativePanelStatus::InvocationFailed)
    {
    }
};

InvocationObservation RunFreshInvocation(
    int target_cpu,
    const NativePanelArm& arm,
    NativePanelRuntime& runtime)
{
    InvocationObservation observation;
    std::unique_ptr<NativePanelInvocation> invocation;
    try {
        invocation = arm.MakeInvocation();
    }
    catch (...) {
        observation.Status = NativePanelStatus::FactoryFailed;
        observation.Diagnostic = "invocation factory threw";
        return observation;
    }
    if (!invocation)
    {
        observation.Status = NativePanelStatus::FactoryFailed;
        observation.Diagnostic = "invocation factory returned null";
        return observation;
    }

    std::string identity;
    try {
        identity = invocation->Identity();
    }
    catch (...) {
        observation.Status = NativePanelStatus::IdentityDrift;
        observation.Diagnostic = "identity query threw before invocation";
        return observation;
    }
    if (identity != arm.ExpectedIdentity)
    {
        observation.Status = NativePanelStatus::IdentityDrift;
        observation.Diagnostic = "fresh invocation identity mismatch";
        return observation;
    }

    std::string cpu_diagnostic;
    NativePanelStatus status = runtime.VerifySingletonCpu(
        target_cpu, cpu_diagnostic);
    if (status != NativePanelStatus::Complete)
    {
        observation.Status = status;
        observation.Diagnostic = cpu_diagnostic;
        return observation;
    }

    bool invocation_threw = false;
    try {
        observation.Result = invocation->Invoke();
    }
    catch (...) {
        invocation_threw = true;
    }

    // A throwing callback is still an invocation: Always perform the required
    // post-invocation CPU and identity checks before reporting its failure.
    cpu_diagnostic.clear();
    status = runtime.VerifySingletonCpu(target_cpu, cpu_diagnostic);
    if (status != NativePanelStatus::Complete)
    {
        observation.Status = status;
        observation.Diagnostic = cpu_diagnostic;
        return observation;
    }

    try {
        identity = invocation->Identity();
    }
    catch (...) {
        observation.Status = NativePanelStatus::IdentityDrift;
        observation.Diagnostic = "identity query threw after invocation";
        return observation;
    }
    if (identity != arm.ExpectedIdentity)
    {
        observation.Status = NativePanelStatus::IdentityDrift;
        observation.Diagnostic = "invocation identity changed";
        return observation;
    }

    if (invocation_threw)
    {
        observation.Status = NativePanelStatus::InvocationFailed;
        observation.Diagnostic = "invocation callback threw";
        return observation;
    }

    if (!IsKnownDisposition(observation.Result.Disposition) ||
        observation.Result.Disposition == NativePanelDisposition::Fatal)
    {
        observation.Status = NativePanelStatus::InvocationFailed;
        observation.Diagnostic = "invocation returned a fatal or invalid outcome";
        return observation;
    }
    if (observation.Result.Disposition == NativePanelDisposition::Success &&
        (observation.Result.ElapsedNanoseconds == 0u ||
         observation.Result.ElapsedNanoseconds > kMaxInt63))
    {
        observation.Status = NativePanelStatus::InvalidElapsed;
        observation.Diagnostic =
            "successful invocation elapsed time is not positive int63";
        return observation;
    }

    observation.Status = NativePanelStatus::Complete;
    observation.Diagnostic.clear();
    return observation;
}

class SystemNativePanelRuntime : public NativePanelRuntime
{
public:
    NativePanelStatus PinAndVerifySingletonCpu(
        int target_cpu,
        std::string& diagnostic) override
    {
#if defined(__linux__)
        if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
        {
            diagnostic = "target CPU is outside cpu_set_t range";
            return NativePanelStatus::InvalidArgument;
        }

        cpu_set_t selected;
        CPU_ZERO(&selected);
        CPU_SET(target_cpu, &selected);
        if (sched_setaffinity(0, sizeof(selected), &selected) != 0)
        {
            std::ostringstream message;
            message << "sched_setaffinity failed: " << strerror(errno);
            diagnostic = message.str();
            return NativePanelStatus::PinFailed;
        }
        return VerifySingletonCpu(target_cpu, diagnostic);
#else
        (void)target_cpu;
        diagnostic = "native timing panels require Linux CPU affinity";
        return NativePanelStatus::UnsupportedPlatform;
#endif
    }

    NativePanelStatus VerifySingletonCpu(
        int target_cpu,
        std::string& diagnostic) override
    {
#if defined(__linux__)
        if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
        {
            diagnostic = "target CPU is outside cpu_set_t range";
            return NativePanelStatus::InvalidArgument;
        }

        cpu_set_t observed_set;
        CPU_ZERO(&observed_set);
        if (sched_getaffinity(0, sizeof(observed_set), &observed_set) != 0)
        {
            std::ostringstream message;
            message << "sched_getaffinity failed: " << strerror(errno);
            diagnostic = message.str();
            return NativePanelStatus::AffinityVerificationFailed;
        }

        unsigned selected_count = 0u;
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
            selected_count += CPU_ISSET(cpu, &observed_set) ? 1u : 0u;
        }
        if (selected_count != 1u || !CPU_ISSET(target_cpu, &observed_set))
        {
            diagnostic = "worker affinity is not the requested singleton CPU";
            return NativePanelStatus::AffinityVerificationFailed;
        }

        errno = 0;
        const int observed_cpu = sched_getcpu();
        if (observed_cpu < 0)
        {
            std::ostringstream message;
            message << "sched_getcpu failed: " << strerror(errno);
            diagnostic = message.str();
            return NativePanelStatus::CpuQueryFailed;
        }
        if (observed_cpu != target_cpu)
        {
            std::ostringstream message;
            message << "worker migrated to CPU " << observed_cpu
                    << " from requested CPU " << target_cpu;
            diagnostic = message.str();
            return NativePanelStatus::CpuMigration;
        }

        diagnostic.clear();
        return NativePanelStatus::Complete;
#else
        (void)target_cpu;
        diagnostic = "native timing panels require Linux CPU affinity";
        return NativePanelStatus::UnsupportedPlatform;
#endif
    }
};

} // namespace

NativePanelInvocationResult::NativePanelInvocationResult()
    : Disposition(NativePanelDisposition::Fatal)
    , OutcomeCode(0)
    , HasDecodedExtra(false)
    , DecodedExtra(0u)
    , ElapsedNanoseconds(0u)
{
}

NativePanelInvocationResult::NativePanelInvocationResult(
    NativePanelDisposition disposition,
    int64_t outcome_code,
    bool has_decoded_extra,
    uint32_t decoded_extra,
    uint64_t elapsed_nanoseconds)
    : Disposition(disposition)
    , OutcomeCode(outcome_code)
    , HasDecodedExtra(has_decoded_extra)
    , DecodedExtra(decoded_extra)
    , ElapsedNanoseconds(elapsed_nanoseconds)
{
}

NativePanelInvocation::~NativePanelInvocation()
{
}

NativePanelArm::NativePanelArm()
{
}

NativePanelArm::NativePanelArm(
    const std::string& expected_identity,
    const NativePanelInvocationFactory& make_invocation)
    : ExpectedIdentity(expected_identity)
    , MakeInvocation(make_invocation)
{
}

const char* NativePanelStatusName(NativePanelStatus status)
{
    switch (status)
    {
    case NativePanelStatus::Complete: return "complete";
    case NativePanelStatus::UnsupportedPlatform: return "unsupported_platform";
    case NativePanelStatus::InvalidArgument: return "invalid_argument";
    case NativePanelStatus::PinFailed: return "pin_failed";
    case NativePanelStatus::AffinityVerificationFailed:
        return "affinity_verification_failed";
    case NativePanelStatus::CpuQueryFailed: return "cpu_query_failed";
    case NativePanelStatus::CpuMigration: return "cpu_migration";
    case NativePanelStatus::FactoryFailed: return "factory_failed";
    case NativePanelStatus::IdentityDrift: return "identity_drift";
    case NativePanelStatus::OutcomeDrift: return "outcome_drift";
    case NativePanelStatus::InvocationFailed: return "invocation_failed";
    case NativePanelStatus::InvalidElapsed: return "invalid_elapsed";
    default: return "unknown";
    }
}

NativePanelSlot::NativePanelSlot()
    : Side(NativePanelSide::Left)
    , HasElapsedNanoseconds(false)
    , ElapsedNanoseconds(0u)
{
}

NativePanelResult::NativePanelResult()
    : Status(NativePanelStatus::InvalidArgument)
    , Order(NativePanelOrder::ABBA)
    , TargetCpu(-1)
    , InvocationsPerSlot(0u)
    , HasLeftPreflight(false)
    , HasRightPreflight(false)
    , PanelComparable(false)
{
}

bool NativePanelResult::IsFatal() const
{
    return Status != NativePanelStatus::Complete;
}

bool NativePanelResult::HasEightNullTimings() const
{
    for (size_t i = 0; i < Slots.size(); ++i) {
        if (Slots[i].HasElapsedNanoseconds) {
            return false;
        }
    }
    return true;
}

bool NativePanelPlatformSupported()
{
#if defined(__linux__)
    return true;
#else
    return false;
#endif
}

NativePanelResult ExecuteNativeTimingPanel(
    int target_cpu,
    NativePanelOrder order,
    uint32_t invocations_per_slot,
    const NativePanelArm& left,
    const NativePanelArm& right)
{
    SystemNativePanelRuntime runtime;
    return ExecuteNativeTimingPanelWithRuntime(
        target_cpu, order, invocations_per_slot, left, right, runtime);
}

NativePanelResult ExecuteNativeTimingPanelWithRuntime(
    int target_cpu,
    NativePanelOrder order,
    uint32_t invocations_per_slot,
    const NativePanelArm& left,
    const NativePanelArm& right,
    NativePanelRuntime& runtime)
{
    NativePanelResult result;
    result.Order = order;
    result.TargetCpu = target_cpu;
    result.InvocationsPerSlot = invocations_per_slot;

    if (order == NativePanelOrder::ABBA)
    {
        result.Slots[0].Side = NativePanelSide::Left;
        result.Slots[1].Side = NativePanelSide::Right;
        result.Slots[2].Side = NativePanelSide::Right;
        result.Slots[3].Side = NativePanelSide::Left;
        result.Slots[4].Side = NativePanelSide::Right;
        result.Slots[5].Side = NativePanelSide::Left;
        result.Slots[6].Side = NativePanelSide::Left;
        result.Slots[7].Side = NativePanelSide::Right;
    }
    else if (order == NativePanelOrder::BAAB)
    {
        result.Slots[0].Side = NativePanelSide::Right;
        result.Slots[1].Side = NativePanelSide::Left;
        result.Slots[2].Side = NativePanelSide::Left;
        result.Slots[3].Side = NativePanelSide::Right;
        result.Slots[4].Side = NativePanelSide::Left;
        result.Slots[5].Side = NativePanelSide::Right;
        result.Slots[6].Side = NativePanelSide::Right;
        result.Slots[7].Side = NativePanelSide::Left;
    }

    if (target_cpu < 0 || invocations_per_slot < 2u ||
        !IsKnownOrder(order) ||
        left.ExpectedIdentity.empty() || right.ExpectedIdentity.empty() ||
        !left.MakeInvocation || !right.MakeInvocation)
    {
        return Fail(result, NativePanelStatus::InvalidArgument,
            "invalid CPU, order, arm identity, or invocation factory");
    }

    std::string diagnostic;
    NativePanelStatus status = runtime.PinAndVerifySingletonCpu(
        target_cpu, diagnostic);
    if (status != NativePanelStatus::Complete) {
        return Fail(result, status, diagnostic);
    }

    InvocationObservation left_warm = RunFreshInvocation(
        target_cpu, left, runtime);
    if (left_warm.Status != NativePanelStatus::Complete) {
        return Fail(result, left_warm.Status, left_warm.Diagnostic);
    }
    result.HasLeftPreflight = true;
    result.LeftPreflight = left_warm.Result;

    InvocationObservation right_warm = RunFreshInvocation(
        target_cpu, right, runtime);
    if (right_warm.Status != NativePanelStatus::Complete) {
        return Fail(result, right_warm.Status, right_warm.Diagnostic);
    }
    result.HasRightPreflight = true;
    result.RightPreflight = right_warm.Result;

    const bool comparable =
        result.LeftPreflight.Disposition == NativePanelDisposition::Success &&
        result.RightPreflight.Disposition == NativePanelDisposition::Success;

    std::array<uint64_t, 8> elapsed_sums = {{ 0u }};
    std::array<NativePanelInvocationResult, 8> final_invocations;
    const uint32_t primary_count =
        invocations_per_slot / 2u + invocations_per_slot % 2u;
    const uint32_t secondary_count = invocations_per_slot / 2u;

    for (size_t block = 0u; block < 2u; ++block)
    {
        const size_t first_slot = block * 4u;
        const uint32_t repeat_count = block == 0u ?
            primary_count : secondary_count;
        for (uint32_t repeat = 0u; repeat < repeat_count; ++repeat)
        {
            for (size_t block_slot = 0u; block_slot < 4u; ++block_slot)
            {
                const size_t slot_index = first_slot + block_slot;
                const NativePanelSlot& slot = result.Slots[slot_index];
                const NativePanelArm& arm = SelectArm(slot.Side, left, right);
                const NativePanelInvocationResult& expected =
                    slot.Side == NativePanelSide::Left ?
                        result.LeftPreflight : result.RightPreflight;
                InvocationObservation measured = RunFreshInvocation(
                    target_cpu, arm, runtime);
                if (measured.Status != NativePanelStatus::Complete) {
                    return Fail(result, measured.Status, measured.Diagnostic);
                }
                if (!SameOutcome(expected, measured.Result))
                {
                    return Fail(result, NativePanelStatus::OutcomeDrift,
                        "measured outcome or decoded-extra drifted from warmup");
                }

                final_invocations[slot_index] = measured.Result;
                if (comparable)
                {
                    const uint64_t elapsed =
                        measured.Result.ElapsedNanoseconds;
                    if (elapsed > kMaxInt63 - elapsed_sums[slot_index])
                    {
                        return Fail(result, NativePanelStatus::InvalidElapsed,
                            "measured batch elapsed time exceeds positive int63");
                    }
                    elapsed_sums[slot_index] += elapsed;
                }
            }
        }
    }

    // Publish the eight slot observations only after both counterbalanced
    // blocks have completed.  Thus every earlier failure returns eight null
    // timings without exposing a completed prefix.
    for (size_t slot_index = 0u; slot_index < result.Slots.size(); ++slot_index)
    {
        NativePanelSlot& slot = result.Slots[slot_index];
        slot.Invocation = final_invocations[slot_index];
        if (comparable)
        {
            slot.HasElapsedNanoseconds = true;
            slot.ElapsedNanoseconds = elapsed_sums[slot_index];
            slot.Invocation.ElapsedNanoseconds = elapsed_sums[slot_index];
        }
    }

    result.Status = NativePanelStatus::Complete;
    result.Diagnostic.clear();
    result.PanelComparable = comparable;
    if (!comparable) {
        ClearTimings(result);
    }
    return result;
}

NativePanelRuntime::~NativePanelRuntime()
{
}

} // namespace wirehair_wh2_bench
