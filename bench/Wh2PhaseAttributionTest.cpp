#include "Wh2PhaseAttribution.h"

#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace {

using wirehair_wh2_bench::NativePanelDisposition;
using wirehair_wh2_bench::NativePanelInvocationResult;
using wirehair_wh2_bench::NativePanelOrder;
using wirehair_wh2_bench::NativePanelResult;
using wirehair_wh2_bench::NativePanelSide;
using wirehair_wh2_bench::NativePanelStatus;
using wirehair_wh2_bench::PhasePanelAssembly;
using wirehair_wh2_bench::PhaseObservationCollector;
using wirehair_wh2_bench::PhaseSolveArm;
using wirehair_wh2_bench::PhaseSolveInvocation;
using wirehair_wh2_bench::PhaseSolveObservation;

int Failures = 0;

void Check(bool condition, const char* what)
{
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++Failures;
    }
}

wirehair_v2::PrecodeSolveStats Counters(uint32_t tag)
{
    wirehair_v2::PrecodeSolveStats stats = {};
    stats.PacketRows = 100u + tag;
    stats.PeeledColumns = 200u + tag;
    stats.InactivatedColumns = 30u + tag;
    stats.ResidualRows = 40u + tag;
    stats.ResidualRank = 30u + tag;
    stats.BinaryResidualRank = 20u + tag;
    stats.BinaryRowReferences = 1000u + tag;
    stats.BinaryRowStorageBytes = 2000u + tag;
    stats.BinaryAdjacencyStorageBytes = 3000u + tag;
    stats.BinaryRowStorageAllocations = 3u + tag;
    stats.BinaryAdjacencyStorageAllocations = 2u + tag;
    stats.BlockXors = 4000u + tag;
    stats.BlockMulAdds = 5000u + tag;
    stats.PacketSeedAttempt = tag;
    return stats;
}

PhaseSolveObservation Success(PhaseSolveArm arm, uint64_t outer)
{
    PhaseSolveObservation observation;
    observation.Arm = arm;
    observation.Result = Wirehair_Success;
    observation.BytesVerified = true;
    observation.ElapsedNanoseconds = outer;
    observation.Stats = Counters(
        arm == PhaseSolveArm::Two07 ? 7u : 11u);
    observation.Stats.BuildNanoseconds = 10u;
    observation.Stats.PeelNanoseconds = 20u;
    observation.Stats.ProjectNanoseconds = 15u;
    observation.Stats.ResidualNanoseconds = 25u;
    observation.Stats.BackSubNanoseconds = 20u;
    return observation;
}

PhaseSolveObservation Weak(PhaseSolveArm arm, WirehairResult result)
{
    PhaseSolveObservation observation;
    observation.Arm = arm;
    observation.Result = result;
    observation.Stats = Counters(
        arm == PhaseSolveArm::Two07 ? 7u : 11u);
    // Failed solve phase timers are diagnostic implementation detail.  They
    // are deliberately not promoted into weak-cell timing evidence.
    observation.Stats.BuildNanoseconds = 3u;
    observation.Stats.PeelNanoseconds = 5u;
    return observation;
}

NativePanelSide ExpectedSide(
    NativePanelOrder order,
    uint32_t block,
    uint32_t block_slot)
{
    const bool abba_left = block_slot == 0u || block_slot == 3u;
    const bool primary_left = order == NativePanelOrder::ABBA ?
        abba_left : !abba_left;
    const bool left = block == 0u ? primary_left : !primary_left;
    return left ? NativePanelSide::Left : NativePanelSide::Right;
}

NativePanelInvocationResult PanelObservation(
    const PhaseSolveObservation& observation)
{
    if (observation.Result == Wirehair_Success)
    {
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            static_cast<int64_t>(Wirehair_Success), true, 4u,
            observation.ElapsedNanoseconds);
    }
    return NativePanelInvocationResult(
        NativePanelDisposition::PreflightFailure,
        static_cast<int64_t>(observation.Result), false, 0u, 0u);
}

struct SyntheticPanel
{
    NativePanelResult Panel;
    std::vector<PhaseSolveObservation> Observations;
};

SyntheticPanel BuildPanel(
    NativePanelOrder order,
    uint32_t n,
    bool weak)
{
    SyntheticPanel built;
    NativePanelResult& panel = built.Panel;
    panel.Status = NativePanelStatus::Complete;
    panel.Order = order;
    panel.TargetCpu = 3;
    panel.InvocationsPerSlot = n;
    panel.HasLeftPreflight = true;
    panel.HasRightPreflight = true;
    const PhaseSolveObservation left = weak ?
        Weak(PhaseSolveArm::Two07, Wirehair_NeedMore) :
        Success(PhaseSolveArm::Two07, 101u);
    const PhaseSolveObservation right =
        Success(PhaseSolveArm::Head, 103u);
    built.Observations.push_back(left);
    built.Observations.push_back(right);
    panel.LeftPreflight = PanelObservation(left);
    panel.RightPreflight = PanelObservation(right);
    panel.PanelComparable = !weak;

    uint64_t sums[8] = {};
    const uint32_t repeats[2] = {
        n / 2u + n % 2u,
        n / 2u
    };
    uint64_t sequence = 0u;
    for (uint32_t block = 0u; block < 2u; ++block)
    {
        for (uint32_t repeat = 0u; repeat < repeats[block]; ++repeat)
        {
            for (uint32_t block_slot = 0u; block_slot < 4u; ++block_slot)
            {
                const uint32_t slot = block * 4u + block_slot;
                const NativePanelSide side =
                    ExpectedSide(order, block, block_slot);
                panel.Slots[slot].Side = side;
                PhaseSolveObservation observation;
                if (side == NativePanelSide::Left && weak) {
                    observation = Weak(
                        PhaseSolveArm::Two07, Wirehair_NeedMore);
                }
                else {
                    observation = Success(
                        side == NativePanelSide::Left ?
                            PhaseSolveArm::Two07 : PhaseSolveArm::Head,
                        100u + (++sequence % 37u));
                }
                built.Observations.push_back(observation);
                if (!weak) sums[slot] += observation.ElapsedNanoseconds;
                panel.Slots[slot].Invocation = PanelObservation(observation);
            }
        }
    }
    for (std::size_t slot = 0u; slot < panel.Slots.size(); ++slot)
    {
        if (!weak)
        {
            panel.Slots[slot].HasElapsedNanoseconds = true;
            panel.Slots[slot].ElapsedNanoseconds = sums[slot];
            panel.Slots[slot].Invocation.ElapsedNanoseconds = sums[slot];
        }
        else
        {
            panel.Slots[slot].HasElapsedNanoseconds = false;
            panel.Slots[slot].ElapsedNanoseconds = 0u;
            panel.Slots[slot].Invocation.ElapsedNanoseconds = 0u;
        }
    }
    return built;
}

bool Validate(
    const SyntheticPanel& built,
    NativePanelOrder order,
    uint32_t n,
    PhasePanelAssembly& assembly)
{
    std::string diagnostic;
    return wirehair_wh2_bench::ValidateAndAssemblePhasePanel(
        built.Panel, order, n, built.Observations, assembly, diagnostic);
}

void TestHappyPaths()
{
    const NativePanelOrder orders[] = {
        NativePanelOrder::ABBA, NativePanelOrder::BAAB
    };
    const uint32_t profiles[] = { 16u, 24u };
    for (NativePanelOrder order : orders)
    {
        for (uint32_t n : profiles)
        {
            const SyntheticPanel built = BuildPanel(order, n, false);
            PhasePanelAssembly assembly;
            Check(Validate(built, order, n, assembly),
                "valid phase panel rejected");
            Check(assembly.Comparable &&
                    assembly.Measured.size() == 4u * n,
                "valid phase chronology not assembled");
            for (std::size_t slot = 0u; slot < assembly.Slots.size(); ++slot)
            {
                Check(assembly.Slots[slot].HasElapsed &&
                        assembly.Slots[slot].OuterNanoseconds ==
                            built.Panel.Slots[slot].ElapsedNanoseconds,
                    "valid phase slot sum not reconstructed");
            }
        }
    }
}

void TestChronologyAndTimingRejections()
{
    const uint32_t n = 16u;
    PhasePanelAssembly assembly;

    SyntheticPanel changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations.pop_back();
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "short chronology accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations[2].Arm = PhaseSolveArm::Head;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "arm tag drift accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations[2].Stats.BackSubNanoseconds = 1000u;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "component sum above outer interval accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations[2].ElapsedNanoseconds = 0u;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "zero successful outer interval accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations[2].Stats.BuildNanoseconds = 0u;
    changed.Observations[2].Stats.PeelNanoseconds = 0u;
    changed.Observations[2].Stats.ProjectNanoseconds = 0u;
    changed.Observations[2].Stats.ResidualNanoseconds = 0u;
    changed.Observations[2].Stats.BackSubNanoseconds = 0u;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "all-zero successful phase timing accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations[2].Stats.ProjectNanoseconds = 0u;
    Check(Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "single zero successful phase component rejected");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations[2].BytesVerified = false;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "successful byte-verification failure accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Observations[2].ElapsedNanoseconds =
        static_cast<uint64_t>(std::numeric_limits<int64_t>::max());
    changed.Observations[2].Stats.BuildNanoseconds = 1u;
    changed.Observations[2].Stats.PeelNanoseconds = 0u;
    changed.Observations[2].Stats.ProjectNanoseconds = 0u;
    changed.Observations[2].Stats.ResidualNanoseconds = 0u;
    changed.Observations[2].Stats.BackSubNanoseconds = 0u;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "positive-int63 aggregate overflow accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Panel.Slots[0].ElapsedNanoseconds += 1u;
    changed.Panel.Slots[0].Invocation.ElapsedNanoseconds += 1u;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "generic/diagnostic outer sum mismatch accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Panel.LeftPreflight.ElapsedNanoseconds += 1u;
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "warmup outer interval mismatch accepted");

    changed = BuildPanel(NativePanelOrder::ABBA, n, false);
    changed.Panel.Diagnostic = "forged complete diagnostic";
    Check(!Validate(changed, NativePanelOrder::ABBA, n, assembly),
        "completed panel with failure diagnostic accepted");
}

void TestCounterAndOutcomeRejections()
{
    const uint32_t n = 16u;
    PhasePanelAssembly assembly;
    for (unsigned field = 0u; field < 14u; ++field)
    {
        SyntheticPanel changed = BuildPanel(
            NativePanelOrder::BAAB, n, false);
        wirehair_v2::PrecodeSolveStats& stats = changed.Observations[2].Stats;
        switch (field)
        {
        case 0u: ++stats.PacketRows; break;
        case 1u: ++stats.PeeledColumns; break;
        case 2u: ++stats.InactivatedColumns; break;
        case 3u: ++stats.ResidualRows; break;
        case 4u: ++stats.ResidualRank; break;
        case 5u: ++stats.BinaryResidualRank; break;
        case 6u: ++stats.BinaryRowReferences; break;
        case 7u: ++stats.BinaryRowStorageBytes; break;
        case 8u: ++stats.BinaryAdjacencyStorageBytes; break;
        case 9u: ++stats.BinaryRowStorageAllocations; break;
        case 10u: ++stats.BinaryAdjacencyStorageAllocations; break;
        case 11u: ++stats.BlockXors; break;
        case 12u: ++stats.BlockMulAdds; break;
        case 13u: ++stats.PacketSeedAttempt; break;
        }
        Check(!Validate(changed, NativePanelOrder::BAAB, n, assembly),
            "non-timing counter drift accepted");
    }

    SyntheticPanel changed = BuildPanel(NativePanelOrder::BAAB, n, false);
    changed.Observations[2].Result = Wirehair_NeedMore;
    changed.Observations[2].BytesVerified = false;
    changed.Observations[2].ElapsedNanoseconds = 0u;
    Check(!Validate(changed, NativePanelOrder::BAAB, n, assembly),
        "measured outcome drift accepted");

    changed = BuildPanel(NativePanelOrder::BAAB, n, false);
    changed.Panel.LeftPreflight.DecodedExtra = 3u;
    Check(!Validate(changed, NativePanelOrder::BAAB, n, assembly),
        "preflight decoded-extra drift accepted");
}

void TestWeakPanel()
{
    const uint32_t n = 24u;
    SyntheticPanel weak = BuildPanel(NativePanelOrder::ABBA, n, true);
    PhasePanelAssembly assembly;
    Check(Validate(weak, NativePanelOrder::ABBA, n, assembly),
        "stable NeedMore weak panel rejected");
    Check(!assembly.Comparable && assembly.Measured.size() == 4u * n,
        "weak panel did not preserve complete ledger");
    for (const auto& slot : assembly.Slots) {
        Check(!slot.HasElapsed && slot.OuterNanoseconds == 0u,
            "weak panel promoted elapsed timing");
    }

    weak = BuildPanel(NativePanelOrder::ABBA, n, true);
    weak.Panel.Slots[0].HasElapsedNanoseconds = true;
    weak.Panel.Slots[0].ElapsedNanoseconds = 1u;
    weak.Panel.Slots[0].Invocation.ElapsedNanoseconds = 1u;
    Check(!Validate(weak, NativePanelOrder::ABBA, n, assembly),
        "weak panel with exposed timing accepted");

    weak = BuildPanel(NativePanelOrder::ABBA, n, true);
    weak.Panel.Slots[0].Invocation.OutcomeCode =
        static_cast<int64_t>(Wirehair_BadPeelSeed);
    Check(!Validate(weak, NativePanelOrder::ABBA, n, assembly),
        "weak terminal slot outcome drift accepted");

    weak = BuildPanel(NativePanelOrder::ABBA, n, true);
    weak.Observations[2].Result = Wirehair_BadPeelSeed;
    Check(!Validate(weak, NativePanelOrder::ABBA, n, assembly),
        "weak outcome drift accepted");

    const WirehairResult construction_failures[] = {
        Wirehair_BadDenseSeed, Wirehair_BadPeelSeed
    };
    for (WirehairResult failure : construction_failures)
    {
        weak = BuildPanel(NativePanelOrder::ABBA, n, true);
        weak.Panel.LeftPreflight.OutcomeCode = static_cast<int64_t>(failure);
        for (std::size_t i = 0u; i < weak.Observations.size(); ++i)
        {
            if (weak.Observations[i].Arm != PhaseSolveArm::Two07) {
                continue;
            }
            weak.Observations[i].Result = failure;
            weak.Observations[i].BytesVerified = false;
            weak.Observations[i].ElapsedNanoseconds = 0u;
        }
        for (auto& slot : weak.Panel.Slots)
        {
            if (slot.Side == NativePanelSide::Left) {
                slot.Invocation.OutcomeCode = static_cast<int64_t>(failure);
            }
        }
        Check(Validate(weak, NativePanelOrder::ABBA, n, assembly),
            "stable construction-failure weak panel rejected");
        Check(!assembly.Comparable,
            "construction-failure weak panel became comparable");
    }
}

void TestInvocationFailsSafelyWithoutCollector()
{
    std::shared_ptr<const wirehair_wh2_bench::NativeSolveFixture> no_fixture;
    std::shared_ptr<PhaseObservationCollector> no_collector;
    PhaseSolveInvocation invocation(
        "phase-test", PhaseSolveArm::Two07, Wirehair_NeedMore,
        no_fixture, 4u, no_collector);
    const NativePanelInvocationResult result = invocation.Invoke();
    Check(result.Disposition == NativePanelDisposition::Fatal &&
            result.OutcomeCode ==
                static_cast<int64_t>(Wirehair_InvalidInput),
        "missing phase collector did not fail safely");
}

} // namespace

int main()
{
    TestHappyPaths();
    TestChronologyAndTimingRejections();
    TestCounterAndOutcomeRejections();
    TestWeakPanel();
    TestInvocationFailsSafelyWithoutCollector();
    if (Failures != 0)
    {
        std::fprintf(stderr,
            "Wh2PhaseAttributionTest: %d failure(s)\n", Failures);
        return 1;
    }
    std::printf("Wh2PhaseAttributionTest: pass\n");
    return 0;
}
