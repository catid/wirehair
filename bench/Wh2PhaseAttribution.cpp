#include "Wh2PhaseAttribution.h"

#if !defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
#error "phase timing requires counter-free benchmark equations"
#elif defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "phase timing must not compile with generic WH2 test hooks"
#endif

#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wirehair_wh2_bench {
namespace {

static const uint64_t kMaxInt63 =
    static_cast<uint64_t>(std::numeric_limits<int64_t>::max());

bool PhaseObservationCount(
    uint32_t invocations_per_slot,
    std::size_t& observation_count)
{
    observation_count = 0u;
#if SIZE_MAX < UINT64_MAX
    const uint64_t count_wide =
        UINT64_C(2) + static_cast<uint64_t>(invocations_per_slot) * 4u;
    if (count_wide > SIZE_MAX) {
        return false;
    }
#endif
    observation_count = 2u +
        static_cast<std::size_t>(invocations_per_slot) * 4u;
    return true;
}

bool KnownArm(PhaseSolveArm arm)
{
    return arm == PhaseSolveArm::Two07 || arm == PhaseSolveArm::Head;
}

bool WeakResult(WirehairResult result)
{
    return result == Wirehair_NeedMore ||
        result == Wirehair_BadDenseSeed ||
        result == Wirehair_BadPeelSeed;
}

NativePanelInvocationResult Classify(
    const PhaseSolveObservation& observation,
    uint32_t fixed_received_overhead)
{
    if (observation.Result == Wirehair_Success)
    {
        if (!observation.BytesVerified ||
            observation.ElapsedNanoseconds == 0u ||
            observation.ElapsedNanoseconds > kMaxInt63)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                static_cast<int64_t>(Wirehair_Error), false, 0u, 0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            static_cast<int64_t>(Wirehair_Success), true,
            fixed_received_overhead, observation.ElapsedNanoseconds);
    }
    if (WeakResult(observation.Result) &&
        !observation.BytesVerified &&
        observation.ElapsedNanoseconds == 0u)
    {
        return NativePanelInvocationResult(
            NativePanelDisposition::PreflightFailure,
            static_cast<int64_t>(observation.Result), false, 0u, 0u);
    }
    return NativePanelInvocationResult(
        NativePanelDisposition::Fatal,
        static_cast<int64_t>(observation.Result), false, 0u, 0u);
}

bool SameOutcome(
    const NativePanelInvocationResult& expected,
    const NativePanelInvocationResult& actual)
{
    return expected.Disposition == actual.Disposition &&
        expected.OutcomeCode == actual.OutcomeCode &&
        expected.HasDecodedExtra == actual.HasDecodedExtra &&
        (!expected.HasDecodedExtra ||
         expected.DecodedExtra == actual.DecodedExtra);
}

wirehair_v2::PrecodeSolveStats NonTimingStats(
    wirehair_v2::PrecodeSolveStats stats)
{
    stats.BuildNanoseconds = 0u;
    stats.PeelNanoseconds = 0u;
    stats.ProjectNanoseconds = 0u;
    stats.ResidualNanoseconds = 0u;
    stats.BackSubNanoseconds = 0u;
    return stats;
}

bool SameNonTimingStats(
    const wirehair_v2::PrecodeSolveStats& a,
    const wirehair_v2::PrecodeSolveStats& b)
{
    const wirehair_v2::PrecodeSolveStats x = NonTimingStats(a);
    const wirehair_v2::PrecodeSolveStats y = NonTimingStats(b);
    return x.PacketRows == y.PacketRows &&
        x.PeeledColumns == y.PeeledColumns &&
        x.InactivatedColumns == y.InactivatedColumns &&
        x.ResidualRows == y.ResidualRows &&
        x.ResidualRank == y.ResidualRank &&
        x.BinaryResidualRank == y.BinaryResidualRank &&
        x.BinaryRowReferences == y.BinaryRowReferences &&
        x.BinaryRowStorageBytes == y.BinaryRowStorageBytes &&
        x.BinaryAdjacencyStorageBytes == y.BinaryAdjacencyStorageBytes &&
        x.BinaryRowStorageAllocations == y.BinaryRowStorageAllocations &&
        x.BinaryAdjacencyStorageAllocations ==
            y.BinaryAdjacencyStorageAllocations &&
        x.BlockXors == y.BlockXors &&
        x.BlockMulAdds == y.BlockMulAdds &&
        x.PacketSeedAttempt == y.PacketSeedAttempt;
}

bool AddInt63(uint64_t value, uint64_t& total)
{
    if (value > kMaxInt63 - total) {
        return false;
    }
    total += value;
    return true;
}

bool ValidSuccessfulTiming(const PhaseSolveObservation& observation)
{
    if (observation.Result != Wirehair_Success ||
        !observation.BytesVerified ||
        observation.ElapsedNanoseconds == 0u ||
        observation.ElapsedNanoseconds > kMaxInt63)
    {
        return false;
    }
    uint64_t phases = 0u;
    return AddInt63(observation.Stats.BuildNanoseconds, phases) &&
        AddInt63(observation.Stats.PeelNanoseconds, phases) &&
        AddInt63(observation.Stats.ProjectNanoseconds, phases) &&
        AddInt63(observation.Stats.ResidualNanoseconds, phases) &&
        AddInt63(observation.Stats.BackSubNanoseconds, phases) &&
        phases != 0u &&
        phases <= observation.ElapsedNanoseconds;
}

PhaseSolveArm ExpectedArm(NativePanelSide side)
{
    return side == NativePanelSide::Left ?
        PhaseSolveArm::Two07 : PhaseSolveArm::Head;
}

bool ExactObservationOutcome(
    const PhaseSolveObservation& observation,
    const NativePanelInvocationResult& panel_outcome,
    uint32_t fixed_received_overhead)
{
    if (!KnownArm(observation.Arm)) {
        return false;
    }
    const NativePanelInvocationResult classified =
        Classify(observation, fixed_received_overhead);
    return classified.Disposition != NativePanelDisposition::Fatal &&
        SameOutcome(classified, panel_outcome);
}

bool AddObservationTiming(
    const PhaseSolveObservation& observation,
    PhaseSlotTotals& totals)
{
    return AddInt63(observation.ElapsedNanoseconds, totals.OuterNanoseconds) &&
        AddInt63(observation.Stats.BuildNanoseconds,
            totals.BuildNanoseconds) &&
        AddInt63(observation.Stats.PeelNanoseconds,
            totals.PeelNanoseconds) &&
        AddInt63(observation.Stats.ProjectNanoseconds,
            totals.ProjectNanoseconds) &&
        AddInt63(observation.Stats.ResidualNanoseconds,
            totals.ResidualNanoseconds) &&
        AddInt63(observation.Stats.BackSubNanoseconds,
            totals.BackSubNanoseconds);
}

bool SlotPhasesFitOuter(const PhaseSlotTotals& totals)
{
    uint64_t phases = 0u;
    return AddInt63(totals.BuildNanoseconds, phases) &&
        AddInt63(totals.PeelNanoseconds, phases) &&
        AddInt63(totals.ProjectNanoseconds, phases) &&
        AddInt63(totals.ResidualNanoseconds, phases) &&
        AddInt63(totals.BackSubNanoseconds, phases) &&
        phases != 0u &&
        phases <= totals.OuterNanoseconds;
}

bool ExpectedSide(
    NativePanelOrder order,
    uint32_t block,
    uint32_t slot,
    NativePanelSide& side)
{
    if ((order != NativePanelOrder::ABBA &&
         order != NativePanelOrder::BAAB) ||
        block > 1u || slot > 3u)
    {
        return false;
    }
    const bool abba_left = slot == 0u || slot == 3u;
    const bool primary_left = order == NativePanelOrder::ABBA ?
        abba_left : !abba_left;
    const bool left = block == 0u ? primary_left : !primary_left;
    side = left ? NativePanelSide::Left : NativePanelSide::Right;
    return true;
}

} // namespace

void PhaseObservationCollector::Append(
    const PhaseSolveObservation& observation)
{
    Values.push_back(observation);
}

const std::vector<PhaseSolveObservation>&
PhaseObservationCollector::Observations() const
{
    return Values;
}

PhaseSolveInvocation::PhaseSolveInvocation(
    const std::string& identity,
    PhaseSolveArm arm,
    WirehairResult preparation_result,
    const std::shared_ptr<const NativeSolveFixture>& fixture,
    uint32_t fixed_received_overhead,
    const std::shared_ptr<PhaseObservationCollector>& collector)
    : IdentityValue(identity)
    , ArmValue(arm)
    , FixedReceivedOverhead(fixed_received_overhead)
    , PreparationResult(preparation_result)
    , Collector(collector)
    , Fixture(fixture)
{
    if (identity.empty() || !KnownArm(arm) || !collector ||
        fixed_received_overhead == 0u ||
        (preparation_result == Wirehair_Success) !=
            (fixture && fixture->IsInitialized()))
    {
        PreparationResult = Wirehair_InvalidInput;
        Fixture.reset();
        return;
    }
}

std::string PhaseSolveInvocation::Identity() const
{
    return IdentityValue;
}

NativePanelInvocationResult PhaseSolveInvocation::Invoke()
{
    if (!Collector)
    {
        return NativePanelInvocationResult(
            NativePanelDisposition::Fatal,
            static_cast<int64_t>(Wirehair_InvalidInput), false, 0u, 0u);
    }
    PhaseSolveObservation observation;
    observation.Arm = ArmValue;
    if (PreparationResult == Wirehair_Success && Fixture) {
        const IsolatedSolveResult result = Fixture->Run();
        observation.Result = result.Result;
        observation.BytesVerified = result.BytesVerified;
        observation.ElapsedNanoseconds = result.ElapsedNanoseconds;
        observation.Stats = result.Stats;
    }
    else {
        observation.Result = PreparationResult;
    }
    Collector->Append(observation);
    return Classify(observation, FixedReceivedOverhead);
}

bool ValidateAndAssemblePhasePanel(
    const NativePanelResult& panel,
    NativePanelOrder expected_order,
    uint32_t invocations_per_slot,
    const std::vector<PhaseSolveObservation>& observations,
    PhasePanelAssembly& result_out,
    std::string& diagnostic)
{
    diagnostic.clear();
    std::size_t expected_observations = 0u;
    if ((expected_order != NativePanelOrder::ABBA &&
         expected_order != NativePanelOrder::BAAB) ||
        invocations_per_slot < 2u ||
        !PhaseObservationCount(
            invocations_per_slot, expected_observations) ||
        observations.size() != expected_observations ||
        panel.Status != NativePanelStatus::Complete ||
        !panel.Diagnostic.empty() ||
        !panel.HasLeftPreflight || !panel.HasRightPreflight ||
        panel.Order != expected_order ||
        panel.InvocationsPerSlot != invocations_per_slot)
    {
        diagnostic = "phase panel shape/status differs from the protocol";
        return false;
    }

    try
    {
        PhasePanelAssembly next;
        next.LeftWarmup = observations[0];
        next.RightWarmup = observations[1];
        if (next.LeftWarmup.Arm != PhaseSolveArm::Two07 ||
            next.RightWarmup.Arm != PhaseSolveArm::Head ||
            !ExactObservationOutcome(
                next.LeftWarmup, panel.LeftPreflight, 4u) ||
            !ExactObservationOutcome(
                next.RightWarmup, panel.RightPreflight, 4u) ||
            panel.LeftPreflight.ElapsedNanoseconds !=
                next.LeftWarmup.ElapsedNanoseconds ||
            panel.RightPreflight.ElapsedNanoseconds !=
                next.RightWarmup.ElapsedNanoseconds)
        {
            diagnostic = "phase warmup identity/outcome is inconsistent";
            return false;
        }
        next.LeftCounters = NonTimingStats(next.LeftWarmup.Stats);
        next.RightCounters = NonTimingStats(next.RightWarmup.Stats);

        const bool comparable =
            next.LeftWarmup.Result == Wirehair_Success &&
            next.RightWarmup.Result == Wirehair_Success;
        next.Comparable = comparable;
        if (panel.PanelComparable != comparable ||
            (comparable && panel.HasEightNullTimings()) ||
            (!comparable && !panel.HasEightNullTimings()))
        {
            diagnostic = "phase panel comparability/null timing drift";
            return false;
        }
        if ((next.LeftWarmup.Result == Wirehair_Success &&
             !ValidSuccessfulTiming(next.LeftWarmup)) ||
            (next.RightWarmup.Result == Wirehair_Success &&
             !ValidSuccessfulTiming(next.RightWarmup)))
        {
            diagnostic = "phase warmup component timing is invalid";
            return false;
        }

        next.Measured.reserve(
            static_cast<std::size_t>(invocations_per_slot) * 4u);
        const uint32_t repeats[2] = {
            invocations_per_slot / 2u + invocations_per_slot % 2u,
            invocations_per_slot / 2u
        };
        std::size_t observation_index = 2u;
        for (uint32_t block = 0u; block < 2u; ++block)
        {
            for (uint32_t repeat = 0u; repeat < repeats[block]; ++repeat)
            {
                for (uint32_t block_slot = 0u; block_slot < 4u; ++block_slot)
                {
                    const uint32_t slot = block * 4u + block_slot;
                    NativePanelSide expected_side = NativePanelSide::Left;
                    if (!ExpectedSide(
                            expected_order, block, block_slot,
                            expected_side) ||
                        panel.Slots[slot].Side != expected_side)
                    {
                        diagnostic = "phase panel slot side/order drift";
                        return false;
                    }
                    const PhaseSolveObservation& observation =
                        observations[observation_index++];
                    const NativePanelInvocationResult& expected_outcome =
                        expected_side == NativePanelSide::Left ?
                            panel.LeftPreflight : panel.RightPreflight;
                    if (observation.Arm != ExpectedArm(expected_side) ||
                        !ExactObservationOutcome(
                            observation, expected_outcome, 4u))
                    {
                        diagnostic =
                            "phase measured arm/outcome chronology drift";
                        return false;
                    }
                    const wirehair_v2::PrecodeSolveStats& arm_stats =
                        expected_side == NativePanelSide::Left ?
                            next.LeftWarmup.Stats : next.RightWarmup.Stats;
                    if (!SameNonTimingStats(observation.Stats, arm_stats))
                    {
                        diagnostic = "phase non-timing counters drifted";
                        return false;
                    }
                    if (observation.Result == Wirehair_Success &&
                        !ValidSuccessfulTiming(observation))
                    {
                        diagnostic =
                            "phase measured component timing is invalid";
                        return false;
                    }
                    PhaseMeasuredObservation tagged;
                    tagged.Block = block;
                    tagged.Repeat = repeat;
                    tagged.Slot = slot;
                    tagged.Observation = observation;
                    next.Measured.push_back(tagged);
                    if (comparable &&
                        !AddObservationTiming(observation, next.Slots[slot]))
                    {
                        diagnostic = "phase slot sum overflowed positive int63";
                        return false;
                    }
                }
            }
        }
        if (observation_index != observations.size())
        {
            diagnostic = "phase observation chronology has trailing entries";
            return false;
        }

        for (std::size_t slot = 0u; slot < next.Slots.size(); ++slot)
        {
            PhaseSlotTotals& totals = next.Slots[slot];
            const NativePanelInvocationResult& expected_outcome =
                panel.Slots[slot].Side == NativePanelSide::Left ?
                    panel.LeftPreflight : panel.RightPreflight;
            if (!SameOutcome(
                    panel.Slots[slot].Invocation, expected_outcome))
            {
                diagnostic = "phase terminal slot outcome drift";
                return false;
            }
            if (comparable)
            {
                totals.HasElapsed = true;
                if (!panel.Slots[slot].HasElapsedNanoseconds ||
                    panel.Slots[slot].ElapsedNanoseconds !=
                        totals.OuterNanoseconds ||
                    panel.Slots[slot].Invocation.ElapsedNanoseconds !=
                        totals.OuterNanoseconds ||
                    !SlotPhasesFitOuter(totals))
                {
                    diagnostic = "phase slot reconstruction mismatch";
                    return false;
                }
            }
            else if (panel.Slots[slot].HasElapsedNanoseconds ||
                panel.Slots[slot].ElapsedNanoseconds != 0u ||
                panel.Slots[slot].Invocation.ElapsedNanoseconds != 0u)
            {
                diagnostic = "weak phase panel exposed elapsed timing";
                return false;
            }
        }
        result_out = std::move(next);
        diagnostic.clear();
        return true;
    }
    catch (const std::bad_alloc&)
    {
        diagnostic = "phase assembly allocation failed";
        return false;
    }
    catch (const std::length_error&)
    {
        diagnostic = "phase assembly length is invalid";
        return false;
    }
}

} // namespace wirehair_wh2_bench
