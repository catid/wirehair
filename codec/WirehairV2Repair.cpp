#include "WirehairV2Repair.h"

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)

#include "WirehairV2Contract.h"
#include "WirehairV2PrecodeDecode.h"

#include "../gf256.h"

#include <cstring>
#include <new>
#include <stdexcept>
#include <vector>

namespace wirehair_v2 {

const char kRepairV1PolicyCanonicalName[] =
    "wirehair:v2:repair-v1:cap8:attempts0-7:"
    "precode-root-plus-attempt-times-9e3779b97f4a7c15-mod2^64:"
    "packet-xorfold-root-plus-attempt-times-9e3779b9-mod2^32:"
    "lazy-real-attempt0-then-zero-rhs-first-full-rank-0-through-7:"
    "selected-real-once:decoder-no-search:2026-07";

const char kRepairV1PolicySha256[] =
    "5e67a150d1f909d6ed80468185fa2dd0"
    "e82eb2fc3486c0fa662e213cf3100b42";

namespace {

static const RepairV1Contract kContracts[] = {
    {
        UINT64_C(0x19cccf775ce0bf09),
        "pure8_s0m1_d3_repair_v1",
        "wirehair:v2:provisional-repair-v1:pure8-s0m1-d3:"
            "smallband-minus1-d3:mixed8gf256-0gf65536-p244-frozen-x:"
            "packet-v4-mix3:raw-uniform:"
            "policy-5e67a150d1f909d6ed80468185fa2dd0"
            "e82eb2fc3486c0fa662e213cf3100b42-cap8:2026-07",
        "19cccf775ce0bf098c9a425cb349714c"
            "4c4a880e7cf136c3bc365e13c05089a5",
        kRepairV1PolicySha256,
        kRepairV1AttemptCount,
        RepairV1Arm::Pure8S0M1D3,
        8u,
        0u,
        3u,
        kCertifiedPacketMixCount
    },
    {
        UINT64_C(0xa530f9105beaa450),
        "pure9_s0m1_d3_repair_v1",
        "wirehair:v2:provisional-repair-v1:pure9-s0m1-d3:"
            "smallband-minus1-d3:mixed9gf256-0gf65536-p244-frozen-x:"
            "packet-v4-mix3:raw-uniform:"
            "policy-5e67a150d1f909d6ed80468185fa2dd0"
            "e82eb2fc3486c0fa662e213cf3100b42-cap8:2026-07",
        "a530f9105beaa450dee70ad9b2a5cc54"
            "c944d3cd47f0aa6534630b8971608541",
        kRepairV1PolicySha256,
        kRepairV1AttemptCount,
        RepairV1Arm::Pure9S0M1D3,
        9u,
        0u,
        3u,
        kCertifiedPacketMixCount
    }
};

bool SameText(const char* left, const char* right)
{
    return left && right && std::strcmp(left, right) == 0;
}

bool IsExactContract(const RepairV1Contract& contract)
{
    const RepairV1Contract* canonical =
        FindRepairV1Contract(contract.ProvisionalId);
    return canonical == &contract;
}

bool SetMixedRowsAndPeriod(
    uint32_t period,
    uint32_t gf256_rows,
    uint32_t gf16_rows)
{
    // Normalize through the production geometry with the maximum period.
    // This ordering is valid from every state accepted by the hook setters,
    // including the experimental 12+4 band.
    if ((ActiveMixedGF256Rows() > kMixedGF256Rows &&
         !SetMixedGF256RowsForTesting(kMixedGF256Rows)) ||
        !SetMixedCoefficientPeriodForTesting(kMixedCoefficientPeriod) ||
        !SetMixedGF256RowsForTesting(kMixedGF256Rows) ||
        !SetMixedGF16RowsForTesting(kMixedGF16Rows))
    {
        return false;
    }

    // H12+ extensions require all four extension rows before the extra
    // subfield rows can be selected.  Every other supported state accepts
    // its extension count before its subfield count.
    if (gf256_rows >= kMixedGF256Rows + 2u)
    {
        if (!SetMixedGF16RowsForTesting(gf16_rows) ||
            !SetMixedCoefficientPeriodForTesting(period) ||
            !SetMixedGF256RowsForTesting(gf256_rows))
        {
            return false;
        }
        return true;
    }
    else if (!SetMixedGF16RowsForTesting(gf16_rows) ||
             !SetMixedGF256RowsForTesting(gf256_rows))
    {
        return false;
    }
    return SetMixedCoefficientPeriodForTesting(period);
}

bool IsCanonicalRepairV1NonGeometryState()
{
    return
        ActiveMixedCoefficientGeometry() ==
            MixedCoefficientGeometry::FrozenPowerX &&
        ActiveMixedResidueSchedule() == MixedResidueSchedule::Constant &&
        ActiveMixedResidueSkew() == 0u &&
        !ActiveMixedResiduesRotated() &&
        !ActiveMixedIndependentExtensionResidues() &&
        ActiveMixedGroupedGF256Rows() == 0u &&
        IsCanonicalStableTargetStaircaseState() &&
        IsCanonicalStaircaseDegreeScaleState() &&
        IsCanonicalStableTargetPacketRowState() &&
        IsCanonicalPacketDegreeState();
}

bool MessageGeometry(
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint32_t& block_count_out)
{
    if (message_bytes == 0u ||
        block_bytes == 0u ||
        block_bytes > UINT32_C(0x7fffffff) ||
        (block_bytes & 1u) != 0u)
    {
        return false;
    }
    const uint64_t count = (message_bytes - 1u) / block_bytes + 1u;
    if (count < 2u || count > 100u) {
        return false;
    }
    block_count_out = (uint32_t)count;
    return true;
}

WirehairResult FinishTelemetry(
    RepairV1SelectorTelemetry& telemetry,
    WirehairResult result)
{
    telemetry.Result = result;
    telemetry.Oom = result == Wirehair_OOM;
    return result;
}

/**
    Pure lazy-policy state machine shared by the real selector and the oracle.
*/
class RepairV1PolicyMachine
{
public:
    explicit RepairV1PolicyMachine(RepairV1SelectorTelemetry& telemetry)
        : Telemetry(telemetry)
    {
    }

    void BeginAttemptZeroReal()
    {
        Telemetry.Attempts[0].RealExecuted = true;
        Telemetry.AttemptsExecuted = 1u;
        ++Telemetry.CallsExecuted;
        ++Telemetry.RealConfigurationCalls;
    }

    WirehairResult FinishAttemptZeroReal(
        WirehairResult result,
        const PrecodeSolveStats* stats,
        bool stats_available)
    {
        RepairV1AttemptTelemetry& observation = Telemetry.Attempts[0];
        if (stats_available && stats)
        {
            observation.RealStats = *stats;
            observation.RealStatsAvailable = true;
        }
        observation.RealResult = result;
        if (result == Wirehair_Success)
        {
            Telemetry.SelectedAttempt = 0u;
            Telemetry.Committed = true;
            return FinishTelemetry(Telemetry, result);
        }
        if (result != Wirehair_NeedMore && result != Wirehair_Error) {
            return FinishTelemetry(Telemetry, result);
        }
        return Wirehair_NeedMore;
    }

    void BeginProbe(uint32_t attempt)
    {
        RepairV1AttemptTelemetry& observation =
            Telemetry.Attempts[attempt];
        observation.ProbeExecuted = true;
        ++Telemetry.CallsExecuted;
        ++Telemetry.StructuralProbeCalls;
        if (Telemetry.AttemptsExecuted < attempt + 1u) {
            Telemetry.AttemptsExecuted = attempt + 1u;
        }
    }

    WirehairResult FinishProbe(
        uint32_t attempt,
        WirehairResult result,
        const PrecodeSolveStats* stats,
        bool stats_available)
    {
        RepairV1AttemptTelemetry& observation =
            Telemetry.Attempts[attempt];
        if (stats_available && stats)
        {
            observation.ProbeStats = *stats;
            observation.ProbeStatsAvailable = true;
        }
        observation.ProbeResult = result;
        if (result == Wirehair_NeedMore) {
            return Wirehair_NeedMore;
        }
        if (result != Wirehair_Success) {
            return FinishTelemetry(Telemetry, result);
        }
        if (attempt == 0u)
        {
            // A full-rank zero-RHS solve cannot follow either a NeedMore or
            // Error real solve for the same equations.  The Error case is
            // called out explicitly by repair-v1 because it distinguishes a
            // structural weak seed from a solver invariant failure.
            Telemetry.FatalAttemptZeroMismatch = true;
            return FinishTelemetry(Telemetry, Wirehair_Error);
        }
        Telemetry.SelectedAttempt = attempt;
        return Wirehair_Success;
    }

    void BeginSelectedReal(uint32_t attempt)
    {
        RepairV1AttemptTelemetry& observation =
            Telemetry.Attempts[attempt];
        observation.RealExecuted = true;
        ++Telemetry.CallsExecuted;
        ++Telemetry.RealConfigurationCalls;
    }

    WirehairResult FinishSelectedReal(
        uint32_t attempt,
        WirehairResult result,
        const PrecodeSolveStats* stats,
        bool stats_available)
    {
        RepairV1AttemptTelemetry& observation =
            Telemetry.Attempts[attempt];
        if (stats_available && stats)
        {
            observation.RealStats = *stats;
            observation.RealStatsAvailable = true;
        }
        observation.RealResult = result;
        if (result == Wirehair_Success) {
            Telemetry.Committed = true;
        }
        return FinishTelemetry(Telemetry, result);
    }

    WirehairResult Exhaust()
    {
        Telemetry.CapExhausted = true;
        Telemetry.SelectedAttempt = kRepairV1NoAttempt;
        return FinishTelemetry(Telemetry, Wirehair_BadPeelSeed);
    }

private:
    RepairV1SelectorTelemetry& Telemetry;
};

/**
    Shared active-call receipt used by the selector's exception handlers and
    the deterministic OOM accounting tests.
*/
class RepairV1ActiveCallTracker
{
public:
    explicit RepairV1ActiveCallTracker(RepairV1PolicyMachine& machine)
        : Machine(machine)
    {
    }

    void BeginAttemptZeroReal()
    {
        Machine.BeginAttemptZeroReal();
        Active = Kind::AttemptZeroReal;
        Attempt = 0u;
    }

    WirehairResult FinishAttemptZeroReal(
        WirehairResult result,
        const PrecodeSolveStats* stats,
        bool stats_available)
    {
        Active = Kind::None;
        return Machine.FinishAttemptZeroReal(
            result, stats, stats_available);
    }

    void BeginProbe(uint32_t attempt)
    {
        Machine.BeginProbe(attempt);
        Active = Kind::Probe;
        Attempt = attempt;
    }

    WirehairResult FinishProbe(
        uint32_t attempt,
        WirehairResult result,
        const PrecodeSolveStats* stats,
        bool stats_available)
    {
        Active = Kind::None;
        return Machine.FinishProbe(
            attempt, result, stats, stats_available);
    }

    void BeginSelectedReal(uint32_t attempt)
    {
        Machine.BeginSelectedReal(attempt);
        Active = Kind::SelectedReal;
        Attempt = attempt;
    }

    WirehairResult FinishSelectedReal(
        uint32_t attempt,
        WirehairResult result,
        const PrecodeSolveStats* stats,
        bool stats_available)
    {
        Active = Kind::None;
        return Machine.FinishSelectedReal(
            attempt, result, stats, stats_available);
    }

    WirehairResult FinishOom()
    {
        const Kind active = Active;
        Active = Kind::None;
        switch (active)
        {
        case Kind::AttemptZeroReal:
            return Machine.FinishAttemptZeroReal(
                Wirehair_OOM, nullptr, false);
        case Kind::Probe:
            return Machine.FinishProbe(
                Attempt, Wirehair_OOM, nullptr, false);
        case Kind::SelectedReal:
            return Machine.FinishSelectedReal(
                Attempt, Wirehair_OOM, nullptr, false);
        case Kind::None:
            break;
        }
        return Wirehair_Error;
    }

    bool HasActiveCall() const
    {
        return Active != Kind::None;
    }

private:
    enum class Kind
    {
        None,
        AttemptZeroReal,
        Probe,
        SelectedReal
    };

    RepairV1PolicyMachine& Machine;
    Kind Active = Kind::None;
    uint32_t Attempt = kRepairV1NoAttempt;
};

} // namespace

const RepairV1Contract* RepairV1Contracts(uint32_t& count_out)
{
    count_out = (uint32_t)(sizeof(kContracts) / sizeof(kContracts[0]));
    return kContracts;
}

const RepairV1Contract* FindRepairV1Contract(RepairV1Arm arm)
{
    uint32_t count = 0u;
    const RepairV1Contract* contracts = RepairV1Contracts(count);
    for (uint32_t i = 0u; i < count; ++i) {
        if (contracts[i].Arm == arm) return &contracts[i];
    }
    return nullptr;
}

const RepairV1Contract* FindRepairV1Contract(uint64_t provisional_id)
{
    uint32_t count = 0u;
    const RepairV1Contract* contracts = RepairV1Contracts(count);
    for (uint32_t i = 0u; i < count; ++i) {
        if (contracts[i].ProvisionalId == provisional_id) {
            return &contracts[i];
        }
    }
    return nullptr;
}

const RepairV1Contract* FindRepairV1Contract(const char* name)
{
    if (!name) return nullptr;
    uint32_t count = 0u;
    const RepairV1Contract* contracts = RepairV1Contracts(count);
    for (uint32_t i = 0u; i < count; ++i) {
        if (SameText(contracts[i].Name, name)) return &contracts[i];
    }
    return nullptr;
}

ScopedRepairV1ContractStateForTesting::
    ScopedRepairV1ContractStateForTesting(
        const RepairV1Contract& contract)
    : TrackingXSnapshot(true, false)
{
    Initialize(&contract);
}

ScopedRepairV1ContractStateForTesting::
    ScopedRepairV1ContractStateForTesting(
        uint64_t provisional_id)
    : TrackingXSnapshot(provisional_id != 0u, false)
{
    Initialize(
        provisional_id == 0u ?
            nullptr : FindRepairV1Contract(provisional_id));
}

void ScopedRepairV1ContractStateForTesting::Initialize(
    const RepairV1Contract* contract)
{
    if (!contract) {
        Valid = false;
        return;
    }
    SavedPeriod = ActiveMixedCoefficientPeriod();
    SavedGF256Rows = ActiveMixedGF256Rows();
    SavedGF16Rows = ActiveMixedGF16Rows();
    HaveSnapshot = true;

    // Tracking-X is process-global at the experiment control surface.  Do not
    // inherit or rewrite it: provisional frozen-X contracts reject a set
    // ambient control and use the scoped false snapshot for the operation.
    if (AmbientMixedBandTrackingXForTesting() ||
        !IsCanonicalRepairV1NonGeometryState() ||
        !IsExactContract(*contract) ||
        !SameText(
            contract->RepairPolicySha256,
            kRepairV1PolicySha256) ||
        contract->SeedAttemptCount != kRepairV1AttemptCount ||
        contract->GF16Rows != 0u ||
        (contract->GF256Rows != 8u && contract->GF256Rows != 9u) ||
        contract->DenseRows != 3u ||
        contract->RecoveryMixCount != kCertifiedPacketMixCount)
    {
        Valid = false;
        return;
    }
    GeometryChanged = true;
    if (!SetMixedRowsAndPeriod(
            kMixedCoefficientPeriod,
            contract->GF256Rows,
            contract->GF16Rows))
    {
        Valid = false;
        return;
    }
    Valid = true;
}

ScopedRepairV1ContractStateForTesting::
    ~ScopedRepairV1ContractStateForTesting() noexcept
{
    if (HaveSnapshot && GeometryChanged)
    {
        const bool restored = SetMixedRowsAndPeriod(
            SavedPeriod, SavedGF256Rows, SavedGF16Rows);
        (void)restored;
    }
}

bool MakeRepairV1ExplicitConfigForTesting(
    const RepairV1Contract& contract,
    uint32_t block_count,
    uint32_t construction_root,
    uint32_t seed_attempt,
    bool cache_systematic_source,
    bool cache_received_systematic,
    ExplicitMessagePrecodeConfigForTesting& config_out)
{
    if (!IsExactContract(contract) ||
        block_count < 2u || block_count > 100u ||
        seed_attempt >= kRepairV1AttemptCount ||
        ActiveMixedCoefficientPeriod() != kMixedCoefficientPeriod ||
        ActiveMixedGF256Rows() != contract.GF256Rows ||
        ActiveMixedGF16Rows() != contract.GF16Rows ||
        AmbientMixedBandTrackingXForTesting() ||
        !IsCanonicalRepairV1NonGeometryState())
    {
        return false;
    }

    ExplicitMessagePrecodeConfigForTesting candidate;
    candidate.Params = MakeCertifiedParams(
        block_count, (uint64_t)construction_root);
    const uint32_t small_band = SmallBandStaircaseCount(block_count);
    if (small_band < 2u) {
        return false;
    }
    candidate.Params.Staircase = small_band - 1u;
    candidate.Params.DenseRows = contract.DenseRows;
    candidate.Params.HeavyRows =
        contract.GF256Rows + contract.GF16Rows;
    candidate.Params.Field = CompletionField::MixedGF256GF16;
    candidate.Params.DenseIdentityCorner = false;
    candidate.Params.DenseTwoAnchor = false;
    candidate.Params.DenseTwoAnchorPhase = 0u;
    candidate.Params.SegmentedDenseAnchors =
        DenseAnchorLayout::Disabled;
    candidate.Params.HeavyFamily =
        HeavyCoefficientFamily::PeriodicCauchy;
    candidate.Params = PrecodeParamsForAttempt(
        candidate.Params, seed_attempt);
    candidate.Packet.PeelSeed = construction_root;
    candidate.Packet.MixCount = contract.RecoveryMixCount;
    candidate.Packet = PacketConfigForAttempt(
        candidate.Packet, seed_attempt);
    candidate.CacheSystematicSource = cache_systematic_source;
    candidate.CacheReceivedSystematicPackets =
        cache_received_systematic;
    if (!candidate.PinActiveEquationStateForTesting()) {
        return false;
    }
    // The repair contract fixes tracking-X off and this scope pins the
    // effective equation generator to that value.  The generic descriptor
    // helper records the process-wide ambient bit, which another thread may
    // toggle after our fail-closed guard above.  Bind the contract value
    // directly so that race cannot label frozen-X equations as tracking-X.
    candidate.EquationState.MixedBandTrackingX = false;
    config_out = candidate;
    return true;
}

WirehairResult EvaluateRepairV1PolicyOracleForTesting(
    WirehairResult attempt_zero_real_result,
    const WirehairResult probe_results[kRepairV1AttemptCount],
    WirehairResult selected_real_result,
    RepairV1SelectorTelemetry& telemetry_out)
{
    telemetry_out = RepairV1SelectorTelemetry();
    if (!probe_results) {
        return FinishTelemetry(telemetry_out, Wirehair_InvalidInput);
    }
    RepairV1PolicyMachine machine(telemetry_out);
    PrecodeSolveStats stats;
    stats.PacketSeedAttempt = 0u;
    const bool attempt_zero_stats_available =
        attempt_zero_real_result == Wirehair_Success ||
        attempt_zero_real_result == Wirehair_NeedMore ||
        attempt_zero_real_result == Wirehair_Error;
    machine.BeginAttemptZeroReal();
    WirehairResult result = machine.FinishAttemptZeroReal(
        attempt_zero_real_result,
        attempt_zero_stats_available ? &stats : nullptr,
        attempt_zero_stats_available);
    if (result != Wirehair_NeedMore ||
        (attempt_zero_real_result != Wirehair_NeedMore &&
         attempt_zero_real_result != Wirehair_Error))
    {
        return result;
    }
    for (uint32_t attempt = 0u;
         attempt < kRepairV1AttemptCount;
         ++attempt)
    {
        stats = PrecodeSolveStats();
        stats.PacketSeedAttempt = attempt;
        const bool probe_stats_available =
            probe_results[attempt] == Wirehair_Success ||
            probe_results[attempt] == Wirehair_NeedMore ||
            probe_results[attempt] == Wirehair_Error;
        machine.BeginProbe(attempt);
        result = machine.FinishProbe(
            attempt,
            probe_results[attempt],
            probe_stats_available ? &stats : nullptr,
            probe_stats_available);
        if (result == Wirehair_NeedMore) {
            continue;
        }
        if (result != Wirehair_Success) {
            return result;
        }
        stats = PrecodeSolveStats();
        stats.PacketSeedAttempt = attempt;
        const bool selected_stats_available =
            selected_real_result == Wirehair_Success ||
            selected_real_result == Wirehair_NeedMore ||
            selected_real_result == Wirehair_Error;
        machine.BeginSelectedReal(attempt);
        return machine.FinishSelectedReal(
            attempt,
            selected_real_result,
            selected_stats_available ? &stats : nullptr,
            selected_stats_available);
    }
    return machine.Exhaust();
}

WirehairResult EvaluateRepairV1AccountingOomForTesting(
    RepairV1AccountingOomStageForTesting stage,
    RepairV1SelectorTelemetry& telemetry_out)
{
    telemetry_out = RepairV1SelectorTelemetry();
    if (stage !=
            RepairV1AccountingOomStageForTesting::
                AttemptZeroBeforeSystemBuild &&
        stage !=
            RepairV1AccountingOomStageForTesting::
                ProbeDuringSystemBuild &&
        stage != RepairV1AccountingOomStageForTesting::ProbeDuringSolve)
    {
        return FinishTelemetry(telemetry_out, Wirehair_InvalidInput);
    }

    RepairV1PolicyMachine machine(telemetry_out);
    RepairV1ActiveCallTracker active(machine);
    active.BeginAttemptZeroReal();
    if (stage ==
        RepairV1AccountingOomStageForTesting::AttemptZeroBeforeSystemBuild)
    {
        return active.FinishOom();
    }

    PrecodeSolveStats real_stats;
    real_stats.PacketSeedAttempt = 0u;
    const WirehairResult real_result = active.FinishAttemptZeroReal(
        Wirehair_NeedMore, &real_stats, true);
    if (real_result != Wirehair_NeedMore) {
        return FinishTelemetry(telemetry_out, Wirehair_Error);
    }

    active.BeginProbe(0u);
    if (stage ==
        RepairV1AccountingOomStageForTesting::ProbeDuringSystemBuild)
    {
        return active.FinishOom();
    }

    PrecodeSolveStats partial_probe_stats;
    partial_probe_stats.PacketRows = 1u;
    partial_probe_stats.PacketSeedAttempt = 0u;
    return active.FinishProbe(
        0u, Wirehair_OOM, &partial_probe_stats, true);
}

WirehairResult MessagePrecodeEncoder::InitializeRepairV1ResultForTesting(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    const RepairV1Contract& contract,
    uint32_t construction_root,
    RepairV1SelectorTelemetry* telemetry,
    bool cache_systematic_source)
{
    RepairV1SelectorTelemetry local;
    RepairV1SelectorTelemetry& receipt = telemetry ? *telemetry : local;
    receipt = RepairV1SelectorTelemetry();

    uint32_t block_count = 0u;
    if (!message ||
        !IsExactContract(contract) ||
        !MessageGeometry(message_bytes, block_bytes, block_count))
    {
        return FinishTelemetry(receipt, Wirehair_InvalidInput);
    }
    RepairV1PolicyMachine machine(receipt);
    RepairV1ActiveCallTracker active(machine);

    try
    {
        const ScopedRepairV1ContractStateForTesting scope(contract);
        if (!scope.IsValid()) {
            return FinishTelemetry(receipt, Wirehair_InvalidInput);
        }
        active.BeginAttemptZeroReal();
        ExplicitMessagePrecodeConfigForTesting config;
        if (!MakeRepairV1ExplicitConfigForTesting(
                contract, block_count, construction_root, 0u,
                cache_systematic_source, false, config))
        {
            return active.FinishAttemptZeroReal(
                Wirehair_InvalidInput, nullptr, false);
        }
        SeedProfile profile = MakeExplicitDiagnosticProfile(
            block_count, block_bytes, config);
        MessagePrecodeEncoderOptions options =
            MakeExplicitDiagnosticOptions(config);
        PrecodeSolveStats real_stats;
        bool real_stats_available = false;
        const WirehairResult attempt_zero_result =
            InitializeConfigurationResult(
                message,
                message_bytes,
                block_bytes,
                config.Params,
                config.Packet,
                profile,
                options,
                0u,
                true,
                false,
                ConfigurationFailurePolicy::Raw,
                &real_stats,
                &real_stats_available);

        WirehairResult result = active.FinishAttemptZeroReal(
            attempt_zero_result,
            real_stats_available ? &real_stats : nullptr,
            real_stats_available);
        if (attempt_zero_result == Wirehair_Success)
        {
            ProfileValue.V2SeedAttempt = 0u;
            SolveStatsValue.PacketSeedAttempt = 0u;
            ExplicitEquationStateValue = config.EquationState;
            ProvisionalRepairContractIdValue = contract.ProvisionalId;
            return result;
        }
        if (result != Wirehair_NeedMore ||
            (attempt_zero_result != Wirehair_NeedMore &&
             attempt_zero_result != Wirehair_Error))
        {
            return result;
        }

        const uint8_t zero_rhs[2] = {0u, 0u};
        std::vector<SolvePacket> packets;

        for (uint32_t attempt = 0u;
             attempt < kRepairV1AttemptCount;
             ++attempt)
        {
            active.BeginProbe(attempt);
            if (attempt == 0u)
            {
                packets.resize(block_count);
                for (uint32_t block_id = 0u;
                     block_id < block_count;
                     ++block_id)
                {
                    packets[block_id].BlockId = block_id;
                    packets[block_id].Data = zero_rhs;
                }
            }
            if (!MakeRepairV1ExplicitConfigForTesting(
                    contract, block_count, construction_root, attempt,
                    cache_systematic_source, false, config))
            {
                return active.FinishProbe(
                    attempt, Wirehair_InvalidInput, nullptr, false);
            }
            PrecodeSystem system;
            if (!BuildPrecodeSystem(config.Params, system))
            {
                return active.FinishProbe(
                    attempt, Wirehair_InvalidInput, nullptr, false);
            }
            std::vector<uint8_t> intermediate;
            PrecodeSolveStats probe_stats;
            probe_stats.PacketSeedAttempt = kRepairV1NoAttempt;
            const WirehairResult probe_result = SolvePrecodeSystem(
                system,
                config.Packet,
                packets,
                2u,
                intermediate,
                &probe_stats);
            const bool probe_stats_available =
                probe_stats.PacketSeedAttempt != kRepairV1NoAttempt;
            if (probe_stats_available) {
                probe_stats.PacketSeedAttempt = attempt;
            }
            result = active.FinishProbe(
                attempt,
                probe_result,
                probe_stats_available ? &probe_stats : nullptr,
                probe_stats_available);
            if (result == Wirehair_NeedMore) {
                continue;
            }
            if (result != Wirehair_Success) {
                return result;
            }

            active.BeginSelectedReal(attempt);
            profile = MakeExplicitDiagnosticProfile(
                block_count, block_bytes, config);
            options = MakeExplicitDiagnosticOptions(config);
            real_stats = PrecodeSolveStats();
            real_stats_available = false;
            const WirehairResult selected_real_result =
                InitializeConfigurationResult(
                    message,
                    message_bytes,
                    block_bytes,
                    config.Params,
                    config.Packet,
                    profile,
                    options,
                    attempt,
                    true,
                    false,
                    ConfigurationFailurePolicy::Raw,
                    &real_stats,
                    &real_stats_available);
            result = active.FinishSelectedReal(
                attempt,
                selected_real_result,
                real_stats_available ? &real_stats : nullptr,
                real_stats_available);
            if (selected_real_result == Wirehair_Success)
            {
                ProfileValue.V2SeedAttempt = attempt;
                SolveStatsValue.PacketSeedAttempt = attempt;
                ExplicitEquationStateValue = config.EquationState;
                ProvisionalRepairContractIdValue =
                    contract.ProvisionalId;
            }
            return result;
        }
        return machine.Exhaust();
    }
    catch (const std::bad_alloc&)
    {
        return active.HasActiveCall() ?
            active.FinishOom() :
            FinishTelemetry(receipt, Wirehair_OOM);
    }
    catch (const std::length_error&)
    {
        return active.HasActiveCall() ?
            active.FinishOom() :
            FinishTelemetry(receipt, Wirehair_OOM);
    }
}

WirehairResult
MessagePrecodeEncoder::InitializeRepairV1SelectedResultForTesting(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    const RepairV1Contract& contract,
    uint32_t construction_root,
    uint32_t seed_attempt,
    bool cache_systematic_source)
{
    uint32_t block_count = 0u;
    if (!message ||
        !IsExactContract(contract) ||
        seed_attempt >= kRepairV1AttemptCount ||
        !MessageGeometry(message_bytes, block_bytes, block_count))
    {
        return Wirehair_InvalidInput;
    }
    try
    {
        const ScopedRepairV1ContractStateForTesting scope(contract);
        if (!scope.IsValid()) {
            return Wirehair_InvalidInput;
        }
        ExplicitMessagePrecodeConfigForTesting config;
        if (!MakeRepairV1ExplicitConfigForTesting(
                contract, block_count, construction_root, seed_attempt,
                cache_systematic_source, false, config))
        {
            return Wirehair_InvalidInput;
        }
        SeedProfile profile = MakeExplicitDiagnosticProfile(
            block_count, block_bytes, config);
        profile.V2SeedAttempt = seed_attempt;
        const MessagePrecodeEncoderOptions options =
            MakeExplicitDiagnosticOptions(config);
        const WirehairResult result = InitializeConfigurationResult(
            message,
            message_bytes,
            block_bytes,
            config.Params,
            config.Packet,
            profile,
            options,
            seed_attempt,
            true,
            false,
            ConfigurationFailurePolicy::Raw);
        if (result == Wirehair_Success)
        {
            ProfileValue.V2SeedAttempt = seed_attempt;
            SolveStatsValue.PacketSeedAttempt = seed_attempt;
            ExplicitEquationStateValue = config.EquationState;
            ProvisionalRepairContractIdValue = contract.ProvisionalId;
        }
        return result;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

WirehairResult MessagePrecodeDecoder::InitializeRepairV1ResultForTesting(
    uint64_t message_bytes,
    uint32_t block_bytes,
    const RepairV1Contract& contract,
    uint32_t construction_root,
    uint32_t seed_attempt,
    bool cache_received_systematic)
{
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }
    uint32_t block_count = 0u;
    if (!IsExactContract(contract) ||
        seed_attempt >= kRepairV1AttemptCount ||
        !MessageGeometry(message_bytes, block_bytes, block_count))
    {
        return Wirehair_InvalidInput;
    }
    try
    {
        const ScopedRepairV1ContractStateForTesting scope(contract);
        if (!scope.IsValid()) {
            return Wirehair_InvalidInput;
        }
        ExplicitMessagePrecodeConfigForTesting config;
        if (!MakeRepairV1ExplicitConfigForTesting(
                contract, block_count, construction_root, seed_attempt,
                false, cache_received_systematic, config))
        {
            return Wirehair_InvalidInput;
        }
        PrecodeSystem system;
        if (!BuildPrecodeSystem(config.Params, system) ||
            !MessagePrecodeEncoder::SameExplicitParams(
                config.Params, system.Params))
        {
            return Wirehair_InvalidInput;
        }
        SeedProfile profile =
            MessagePrecodeEncoder::MakeExplicitDiagnosticProfile(
                block_count, block_bytes, config);
        profile.V2SeedAttempt = seed_attempt;
        const MessagePrecodeEncoderOptions options =
            MessagePrecodeEncoder::MakeExplicitDiagnosticOptions(config);
        const WirehairResult result = InitializeConfigurationResult(
            message_bytes,
            block_bytes,
            std::move(system),
            config.Packet,
            profile,
            options,
            seed_attempt,
            false);
        if (result == Wirehair_Success)
        {
            ProfileValue.V2SeedAttempt = seed_attempt;
            ExplicitEquationStateValue = config.EquationState;
            ProvisionalRepairContractIdValue = contract.ProvisionalId;
        }
        return result;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

} // namespace wirehair_v2

#endif // WIREHAIR_V2_ENABLE_TEST_HOOKS
