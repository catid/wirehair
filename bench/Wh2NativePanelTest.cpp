#include "Wh2NativePanel.h"

#include <cstdio>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

using wirehair_wh2_bench::ExecuteNativeTimingPanelWithRuntime;
using wirehair_wh2_bench::NativePanelArm;
using wirehair_wh2_bench::NativePanelDisposition;
using wirehair_wh2_bench::NativePanelInvocation;
using wirehair_wh2_bench::NativePanelInvocationResult;
using wirehair_wh2_bench::NativePanelOrder;
using wirehair_wh2_bench::NativePanelResult;
using wirehair_wh2_bench::NativePanelRuntime;
using wirehair_wh2_bench::NativePanelSide;
using wirehair_wh2_bench::NativePanelStatus;

int Failures = 0;

bool Check(bool condition, const char* message)
{
    if (!condition)
    {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++Failures;
    }
    return condition;
}

class FakeRuntime : public NativePanelRuntime
{
public:
    NativePanelStatus PinStatus;
    NativePanelStatus VerifyFailureStatus;
    int FailVerifyAt;
    int PinCalls;
    int VerifyCalls;
    int LastTargetCpu;

    FakeRuntime()
        : PinStatus(NativePanelStatus::Complete)
        , VerifyFailureStatus(NativePanelStatus::CpuMigration)
        , FailVerifyAt(-1)
        , PinCalls(0)
        , VerifyCalls(0)
        , LastTargetCpu(-1)
    {
    }

    NativePanelStatus PinAndVerifySingletonCpu(
        int target_cpu,
        std::string& diagnostic) override
    {
        ++PinCalls;
        LastTargetCpu = target_cpu;
        if (PinStatus != NativePanelStatus::Complete) {
            diagnostic = "injected pin failure";
        } else {
            diagnostic.clear();
        }
        return PinStatus;
    }

    NativePanelStatus VerifySingletonCpu(
        int target_cpu,
        std::string& diagnostic) override
    {
        ++VerifyCalls;
        LastTargetCpu = target_cpu;
        if (VerifyCalls == FailVerifyAt)
        {
            diagnostic = "injected CPU migration";
            return VerifyFailureStatus;
        }
        diagnostic.clear();
        return NativePanelStatus::Complete;
    }
};

struct FakeArmState
{
    char Tag;
    std::string Identity;
    NativePanelDisposition Disposition;
    int64_t OutcomeCode;
    bool HasDecodedExtra;
    uint32_t DecodedExtra;
    uint64_t BaseElapsed;
    int Created;
    int Invoked;
    int DriftOutcomeSerial;
    int DriftExtraSerial;
    int DriftIdentityAfterInvokeSerial;
    int InvalidElapsedSerial;
    int ThrowSerial;
    std::vector<char>* Chronology;

    FakeArmState(
        char tag,
        const char* identity,
        std::vector<char>& chronology)
        : Tag(tag)
        , Identity(identity)
        , Disposition(NativePanelDisposition::Success)
        , OutcomeCode(tag == 'L' ? 11 : 22)
        , HasDecodedExtra(true)
        , DecodedExtra(tag == 'L' ? 1u : 2u)
        , BaseElapsed(tag == 'L' ? 100u : 200u)
        , Created(0)
        , Invoked(0)
        , DriftOutcomeSerial(-1)
        , DriftExtraSerial(-1)
        , DriftIdentityAfterInvokeSerial(-1)
        , InvalidElapsedSerial(-1)
        , ThrowSerial(-1)
        , Chronology(&chronology)
    {
    }
};

class FakeInvocation : public NativePanelInvocation
{
public:
    FakeInvocation(FakeArmState& state, int serial)
        : State(state)
        , Serial(serial)
        , WasInvoked(false)
    {
    }

    std::string Identity() const override
    {
        if (WasInvoked && Serial == State.DriftIdentityAfterInvokeSerial) {
            return State.Identity + "-drift";
        }
        return State.Identity;
    }

    NativePanelInvocationResult Invoke() override
    {
        if (WasInvoked)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                -99,
                false,
                0u,
                0u);
        }
        WasInvoked = true;
        ++State.Invoked;
        State.Chronology->push_back(State.Tag);
        if (Serial == State.ThrowSerial) {
            throw 17;
        }

        int64_t outcome_code = State.OutcomeCode;
        uint32_t decoded_extra = State.DecodedExtra;
        if (Serial == State.DriftOutcomeSerial) {
            ++outcome_code;
        }
        if (Serial == State.DriftExtraSerial) {
            ++decoded_extra;
        }
        const uint64_t elapsed = Serial == State.InvalidElapsedSerial ?
            0u : State.BaseElapsed + (uint64_t)Serial;
        return NativePanelInvocationResult(
            State.Disposition,
            outcome_code,
            State.HasDecodedExtra,
            decoded_extra,
            elapsed);
    }

private:
    FakeArmState& State;
    int Serial;
    bool WasInvoked;
};

NativePanelArm MakeFakeArm(FakeArmState& state)
{
    return NativePanelArm(
        state.Identity,
        [&state]() -> std::unique_ptr<NativePanelInvocation> {
            const int serial = state.Created++;
            return std::unique_ptr<NativePanelInvocation>(
                new FakeInvocation(state, serial));
        });
}

bool SameChronology(
    const std::vector<char>& actual,
    const char* expected)
{
    const std::string text(expected);
    return actual == std::vector<char>(text.begin(), text.end());
}

void CheckSuccessfulPanel(
    NativePanelOrder order,
    const char* expected_chronology,
    const NativePanelSide expected_sides[4],
    const uint64_t expected_elapsed[4])
{
    std::vector<char> chronology;
    FakeArmState left('L', "left-id", chronology);
    FakeArmState right('R', "right-id", chronology);
    FakeRuntime runtime;

    const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
        17, order, 2u, MakeFakeArm(left), MakeFakeArm(right), runtime);

    Check(result.Status == NativePanelStatus::Complete,
        "successful panel did not complete");
    Check(!result.IsFatal(), "successful panel reported fatal");
    Check(result.PanelComparable, "successful panel is not comparable");
    Check(!result.HasFourNullTimings(),
        "successful panel emitted four null timings");
    Check(result.HasLeftPreflight && result.HasRightPreflight,
        "successful panel omitted a warmup outcome");
    Check(result.InvocationsPerSlot == 2u,
        "successful panel did not receipt its exact batch count");
    Check(left.Created == 5 && right.Created == 5,
        "panel did not create one warmup plus four fresh callbacks per side");
    Check(left.Invoked == 5 && right.Invoked == 5,
        "panel did not invoke each fresh callback exactly once");
    Check(SameChronology(chronology, expected_chronology),
        "panel callback chronology changed");
    Check(runtime.PinCalls == 1 && runtime.VerifyCalls == 20,
        "panel did not verify CPU around all ten fresh callbacks");
    Check(runtime.LastTargetCpu == 17,
        "panel runtime observed the wrong target CPU");

    for (size_t i = 0; i < result.Slots.size(); ++i)
    {
        Check(result.Slots[i].Side == expected_sides[i],
            "measured side chronology changed");
        Check(result.Slots[i].HasElapsedNanoseconds,
            "successful measured slot has null elapsed time");
        Check(result.Slots[i].ElapsedNanoseconds == expected_elapsed[i],
            "measured elapsed value was reordered");
        Check(result.Slots[i].Invocation.Disposition ==
                  NativePanelDisposition::Success,
            "successful slot lost its invocation outcome");
    }
}

void CheckAbbaAndBaab()
{
    const NativePanelSide abba_sides[4] = {
        NativePanelSide::Left,
        NativePanelSide::Right,
        NativePanelSide::Right,
        NativePanelSide::Left
    };
    const uint64_t abba_elapsed[4] = { 203u, 403u, 407u, 207u };
    CheckSuccessfulPanel(
        NativePanelOrder::ABBA,
        "LRLLRRRRLL",
        abba_sides,
        abba_elapsed);

    const NativePanelSide baab_sides[4] = {
        NativePanelSide::Right,
        NativePanelSide::Left,
        NativePanelSide::Left,
        NativePanelSide::Right
    };
    const uint64_t baab_elapsed[4] = { 403u, 203u, 207u, 407u };
    CheckSuccessfulPanel(
        NativePanelOrder::BAAB,
        "LRRRLLLLRR",
        baab_sides,
        baab_elapsed);
}

void CheckPreflightFailureProducesFourNulls()
{
    std::vector<char> chronology;
    FakeArmState left('L', "left-id", chronology);
    FakeArmState right('R', "right-id", chronology);
    left.Disposition = NativePanelDisposition::PreflightFailure;
    left.OutcomeCode = -7;
    left.HasDecodedExtra = false;
    left.BaseElapsed = 0u;
    FakeRuntime runtime;

    const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
        3, NativePanelOrder::ABBA, 2u,
        MakeFakeArm(left), MakeFakeArm(right), runtime);

    Check(result.Status == NativePanelStatus::Complete,
        "stable preflight failure was treated as fatal");
    Check(!result.PanelComparable,
        "panel with a failed preflight is comparable");
    Check(result.HasFourNullTimings(),
        "preflight failure did not clear all four timings");
    for (size_t i = 0; i < result.Slots.size(); ++i)
    {
        Check(result.Slots[i].ElapsedNanoseconds == 0u &&
              result.Slots[i].Invocation.ElapsedNanoseconds == 0u,
            "null timing retained a hidden elapsed value");
    }
    Check(SameChronology(chronology, "LRLLRRRRLL"),
        "preflight failure skipped or reordered measured callbacks");
    Check(result.InvocationsPerSlot == 2u &&
          left.Created == 5 && right.Created == 5 &&
          left.Invoked == 5 && right.Invoked == 5,
        "preflight failure did not retain fresh callback discipline");
}

void CheckOutcomeAndExtraDriftAreFatal()
{
    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        left.DriftOutcomeSerial = 2;
        FakeRuntime runtime;
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            4, NativePanelOrder::ABBA, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::OutcomeDrift,
            "outcome-code drift was accepted");
        Check(result.HasFourNullTimings(),
            "outcome drift retained partial timings");
        Check(SameChronology(chronology, "LRLL"),
            "later-repeat outcome drift ran an extra callback");
    }

    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        right.DriftExtraSerial = 2;
        FakeRuntime runtime;
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            4, NativePanelOrder::ABBA, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::OutcomeDrift,
            "decoded-extra drift was accepted");
        Check(result.HasFourNullTimings(),
            "decoded-extra drift retained partial timings");
    }
}

void CheckIdentityAndTimingFailuresAreFatal()
{
    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        left.DriftIdentityAfterInvokeSerial = 2;
        FakeRuntime runtime;
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            5, NativePanelOrder::ABBA, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::IdentityDrift,
            "post-callback identity drift was accepted");
        Check(result.HasFourNullTimings(),
            "identity drift retained partial timings");
    }

    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        left.InvalidElapsedSerial = 2;
        FakeRuntime runtime;
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            5, NativePanelOrder::ABBA, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::InvalidElapsed,
            "zero successful elapsed time was accepted");
        Check(result.HasFourNullTimings(),
            "invalid elapsed time retained partial timings");
    }

    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        left.BaseElapsed = UINT64_MAX;
        FakeRuntime runtime;
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            5, NativePanelOrder::ABBA, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::InvalidElapsed,
            "successful elapsed time above int63 was accepted");
        Check(chronology == std::vector<char>(1u, 'L'),
            "invalid warmup elapsed ran another callback");
        Check(result.HasFourNullTimings(),
            "elapsed time above int63 retained timings");
    }

    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        left.BaseElapsed =
            static_cast<uint64_t>(std::numeric_limits<int64_t>::max()) / 2u;
        FakeRuntime runtime;
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            5, NativePanelOrder::ABBA, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::InvalidElapsed,
            "overflowing batch elapsed sum was accepted");
        Check(SameChronology(chronology, "LRLL"),
            "batch overflow did not fail at the exact later repeat");
        Check(result.HasFourNullTimings(),
            "batch elapsed overflow retained partial timings");
    }
}

void CheckThrowStillGetsPostCpuVerification()
{
    std::vector<char> chronology;
    FakeArmState left('L', "left-id", chronology);
    FakeArmState right('R', "right-id", chronology);
    left.ThrowSerial = 2;
    FakeRuntime runtime;
    const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
        6, NativePanelOrder::ABBA, 2u,
        MakeFakeArm(left), MakeFakeArm(right), runtime);
    Check(result.Status == NativePanelStatus::InvocationFailed,
        "throwing callback was not rejected");
    Check(runtime.VerifyCalls == 8,
        "throwing callback skipped its post-invocation CPU check");
    Check(result.HasFourNullTimings(),
        "throwing callback retained timings");
}

void CheckWrongCoreAndMigrationAreFatal()
{
    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        FakeRuntime runtime;
        runtime.PinStatus = NativePanelStatus::CpuMigration;
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            9, NativePanelOrder::ABBA, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::CpuMigration,
            "wrong initial CPU was accepted");
        Check(chronology.empty(),
            "callback ran after initial CPU rejection");
        Check(result.HasFourNullTimings(),
            "initial CPU rejection retained timings");
    }

    {
        std::vector<char> chronology;
        FakeArmState left('L', "left-id", chronology);
        FakeArmState right('R', "right-id", chronology);
        FakeRuntime runtime;
        runtime.FailVerifyAt = 2; // Immediately after warm-left Invoke().
        const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
            9, NativePanelOrder::BAAB, 2u,
            MakeFakeArm(left), MakeFakeArm(right), runtime);
        Check(result.Status == NativePanelStatus::CpuMigration,
            "post-callback migration was accepted");
        Check(chronology == std::vector<char>(1u, 'L'),
            "migration rejection ran an extra callback");
        Check(result.HasFourNullTimings(),
            "migration rejection retained timings");
    }
}

void CheckInvalidInputsFailClosed()
{
    std::vector<char> chronology;
    FakeArmState left('L', "left-id", chronology);
    FakeArmState right('R', "right-id", chronology);
    FakeRuntime runtime;
    const NativePanelResult bad_order = ExecuteNativeTimingPanelWithRuntime(
        1,
        static_cast<NativePanelOrder>(99u),
        1u,
        MakeFakeArm(left),
        MakeFakeArm(right),
        runtime);
    Check(bad_order.Status == NativePanelStatus::InvalidArgument,
        "invalid chronology was accepted");
    Check(bad_order.HasFourNullTimings(),
        "invalid chronology retained timings");
    Check(runtime.PinCalls == 0 && chronology.empty(),
        "invalid chronology reached the runtime or callback");

    const NativePanelResult zero_batch = ExecuteNativeTimingPanelWithRuntime(
        1,
        NativePanelOrder::ABBA,
        0u,
        MakeFakeArm(left),
        MakeFakeArm(right),
        runtime);
    Check(zero_batch.Status == NativePanelStatus::InvalidArgument &&
            zero_batch.InvocationsPerSlot == 0u &&
            zero_batch.HasFourNullTimings(),
        "zero invocation batch did not fail closed");
    Check(runtime.PinCalls == 0 && chronology.empty(),
        "zero invocation batch reached the runtime or callback");
}

void CheckAaBatchUsesFreshInvocations()
{
    std::vector<char> chronology;
    FakeArmState arm('A', "same-id", chronology);
    arm.OutcomeCode = 33;
    arm.DecodedExtra = 4u;
    FakeRuntime runtime;
    const NativePanelResult result = ExecuteNativeTimingPanelWithRuntime(
        11, NativePanelOrder::ABBA, 3u,
        MakeFakeArm(arm), MakeFakeArm(arm), runtime);

    Check(result.Status == NativePanelStatus::Complete &&
            result.PanelComparable && result.InvocationsPerSlot == 3u,
        "A/A batched panel did not complete with its exact count");
    Check(arm.Created == 14 && arm.Invoked == 14 &&
            chronology == std::vector<char>(14u, 'A'),
        "A/A batch reused or skipped a supposedly fresh invocation");
    Check(runtime.VerifyCalls == 28,
        "A/A batch did not verify affinity around every fresh invocation");
    for (std::size_t slot = 0u; slot < result.Slots.size(); ++slot) {
        Check(result.Slots[slot].HasElapsedNanoseconds &&
                result.Slots[slot].ElapsedNanoseconds ==
                    result.Slots[slot].Invocation.ElapsedNanoseconds,
            "A/A batch did not expose the aggregate elapsed sum");
    }
}

} // namespace

int main()
{
    CheckAbbaAndBaab();
    CheckPreflightFailureProducesFourNulls();
    CheckOutcomeAndExtraDriftAreFatal();
    CheckIdentityAndTimingFailuresAreFatal();
    CheckThrowStillGetsPostCpuVerification();
    CheckWrongCoreAndMigrationAreFatal();
    CheckInvalidInputsFailClosed();
    CheckAaBatchUsesFreshInvocations();

    if (Failures != 0)
    {
        std::fprintf(stderr, "Wh2NativePanelTest: %d failure(s)\n", Failures);
        return 1;
    }
    std::printf("Wh2NativePanelTest: pass\n");
    return 0;
}
