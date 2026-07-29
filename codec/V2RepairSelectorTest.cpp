#include "WirehairV2Contract.h"
#include "WirehairV2Fingerprint.h"
#include "WirehairV2PrecodeDecode.h"
#include "WirehairV2Repair.h"

#include <atomic>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>

namespace {

class Digest64
{
public:
    void U8(uint8_t value)
    {
        State ^= value;
        State *= UINT64_C(1099511628211);
    }
    void U32(uint32_t value)
    {
        for (unsigned i = 0u; i < 4u; ++i) {
            U8((uint8_t)(value >> (8u * i)));
        }
    }
    void U64(uint64_t value)
    {
        for (unsigned i = 0u; i < 8u; ++i) {
            U8((uint8_t)(value >> (8u * i)));
        }
    }
    uint64_t Get() const { return State; }

private:
    uint64_t State = UINT64_C(1469598103934665603);
};

void FillMessage(
    std::vector<uint8_t>& message,
    uint32_t K,
    uint32_t width,
    uint32_t variant)
{
    message.resize((size_t)K * width);
    uint64_t state =
        UINT64_C(0x243f6a8885a308d3) ^
        ((uint64_t)K << 32) ^ width ^ ((uint64_t)variant << 48);
    for (uint8_t& byte : message)
    {
        state += UINT64_C(0x9e3779b97f4a7c15);
        uint64_t z = state;
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        byte = (uint8_t)(z ^ (z >> 31));
    }
}

bool SameParams(
    const wirehair_v2::PrecodeParams& a,
    const wirehair_v2::PrecodeParams& b)
{
    return a.BlockCount == b.BlockCount &&
        a.Staircase == b.Staircase &&
        a.DenseRows == b.DenseRows &&
        a.HeavyRows == b.HeavyRows &&
        a.SourceHits == b.SourceHits &&
        a.Field == b.Field &&
        a.DegreeBalancedStaircase == b.DegreeBalancedStaircase &&
        a.StaircaseDegreeScale == b.StaircaseDegreeScale &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.DenseTwoAnchor == b.DenseTwoAnchor &&
        a.DenseTwoAnchorPhase == b.DenseTwoAnchorPhase &&
        a.SegmentedDenseAnchors == b.SegmentedDenseAnchors &&
        a.HeavyFamily == b.HeavyFamily &&
        a.Seed == b.Seed;
}

uint64_t EquationDigest(
    const wirehair_v2::RepairV1Contract& contract,
    const wirehair_v2::ExplicitMessagePrecodeConfigForTesting& config,
    bool& valid)
{
    valid = false;
    wirehair_v2::ScopedRepairV1ContractStateForTesting scope(contract);
    wirehair_v2::PrecodeSystem system;
    if (!scope.IsValid() ||
        !wirehair_v2::BuildPrecodeSystem(config.Params, system))
    {
        return 0u;
    }
    Digest64 hash;
    hash.U32(system.Params.BlockCount);
    hash.U32(system.Params.Staircase);
    hash.U32(system.Params.DenseRows);
    hash.U32(system.Params.HeavyRows);
    hash.U64(system.Params.Seed);
    for (const std::vector<uint32_t>& row : system.StaircaseRows)
    {
        hash.U32((uint32_t)row.size());
        for (uint32_t column : row) hash.U32(column);
    }
    for (const std::vector<uint32_t>& row : system.DenseRowColumns)
    {
        hash.U32((uint32_t)row.size());
        for (uint32_t column : row) hash.U32(column);
    }
    const wirehair_v2::MixedCoefficientRows* coefficients =
        wirehair_v2::GetMixedCoefficientRows();
    if (!coefficients) return 0u;
    const uint32_t columns =
        system.Params.BlockCount + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    for (uint32_t row = 0u; row < contract.GF256Rows; ++row)
    {
        for (uint32_t column = 0u; column < columns; ++column)
        {
            hash.U8(coefficients->Subfield[row][
                wirehair_v2::ActiveMixedCoefficientResidue(column)]);
        }
    }
    for (uint32_t block_id = 0u;
         block_id < system.Params.BlockCount + 8u;
         ++block_id)
    {
        const std::vector<uint32_t> row =
            wirehair_v2::GeneratePacketMatrixRow(
                system.Params.BlockCount,
                system.Params.Staircase + system.Params.DenseRows +
                    system.Params.HeavyRows,
                block_id,
                config.Packet);
        if (row.empty()) return 0u;
        hash.U32(block_id);
        hash.U32((uint32_t)row.size());
        for (uint32_t column : row) hash.U32(column);
    }
    valid = true;
    return hash.Get();
}

bool CheckContractsAndVectors()
{
    static const uint64_t expected_precode[8] = {
        UINT64_C(0x0000000089abcdef),
        UINT64_C(0x9e3779ba08f64a04),
        UINT64_C(0x3c6ef3738840c619),
        UINT64_C(0xdaa66d2d078b422e),
        UINT64_C(0x78dde6e686d5be43),
        UINT64_C(0x171560a006203a58),
        UINT64_C(0xb54cda59856ab66d),
        UINT64_C(0x5384541304b53282)
    };
    static const uint32_t expected_packet[8] = {
        UINT32_C(0x89abcdef), UINT32_C(0x27e347a8),
        UINT32_C(0xc61ac161), UINT32_C(0x64523b1a),
        UINT32_C(0x0289b4d3), UINT32_C(0xa0c12e8c),
        UINT32_C(0x3ef8a845), UINT32_C(0xdd3021fe)
    };

    uint32_t count = 0u;
    const wirehair_v2::RepairV1Contract* contracts =
        wirehair_v2::RepairV1Contracts(count);
    wirehair_v2::V2EquationContract digest_contract = {};
    digest_contract.CanonicalName =
        wirehair_v2::kRepairV1PolicyCanonicalName;
    uint8_t digest[wirehair_v2::kEquationFingerprintBytes] = {};
    char digest_hex[wirehair_v2::kEquationFingerprintHexChars + 1u] = {};
    if (!wirehair_v2::ComputeEquationContractNameDigest(
            digest_contract, digest))
    {
        return false;
    }
    wirehair_v2::FormatEquationFingerprintHex(digest, digest_hex);
    uint32_t equation_count = 0u;
    (void)wirehair_v2::V2EquationContracts(equation_count);
    if (count != 2u || equation_count != 5u ||
        std::strcmp(
            wirehair_v2::kRepairV1PolicySha256,
            "5e67a150d1f909d6ed80468185fa2dd0"
            "e82eb2fc3486c0fa662e213cf3100b42") != 0 ||
        std::strcmp(
            wirehair_v2::kRepairV1PolicySha256, digest_hex) != 0 ||
        wirehair_v2::FindV2EquationContract(
            UINT64_C(0x19cccf775ce0bf09)) != nullptr ||
        wirehair_v2::FindV2EquationContract(
            UINT64_C(0xa530f9105beaa450)) != nullptr)
    {
        std::fprintf(stderr, "repair-v1 contract isolation mismatch\n");
        return false;
    }
    static const uint64_t expected_ids[2] = {
        UINT64_C(0x19cccf775ce0bf09),
        UINT64_C(0xa530f9105beaa450)
    };
    static const char* const expected_hashes[2] = {
        "19cccf775ce0bf098c9a425cb349714c"
            "4c4a880e7cf136c3bc365e13c05089a5",
        "a530f9105beaa450dee70ad9b2a5cc54"
            "c944d3cd47f0aa6534630b8971608541"
    };
    for (uint32_t arm = 0u; arm < count; ++arm)
    {
        const wirehair_v2::RepairV1Contract& contract = contracts[arm];
        digest_contract.CanonicalName = contract.CanonicalName;
        if (!wirehair_v2::ComputeEquationContractNameDigest(
                digest_contract, digest))
        {
            return false;
        }
        wirehair_v2::FormatEquationFingerprintHex(digest, digest_hex);
        uint64_t digest_id = 0u;
        for (uint32_t i = 0u; i < 8u; ++i) {
            digest_id = (digest_id << 8) | digest[i];
        }
        if (contract.ProvisionalId != expected_ids[arm] ||
            std::strcmp(contract.CanonicalSha256,
                expected_hashes[arm]) != 0 ||
            std::strcmp(contract.CanonicalSha256, digest_hex) != 0 ||
            contract.ProvisionalId != digest_id ||
            std::strcmp(
                contract.RepairPolicySha256,
                wirehair_v2::kRepairV1PolicySha256) != 0 ||
            contract.SeedAttemptCount != 8u ||
            wirehair_v2::FindRepairV1Contract(contract.Arm) !=
                &contract ||
            wirehair_v2::FindRepairV1Contract(contract.Name) !=
                &contract)
        {
            std::fprintf(stderr,
                "repair-v1 arm identity mismatch arm=%u\n", arm);
            return false;
        }
        wirehair_v2::ScopedRepairV1ContractStateForTesting scope(contract);
        if (!scope.IsValid()) return false;
        wirehair_v2::ExplicitMessagePrecodeConfigForTesting previous;
        previous.Params.Seed = UINT64_C(0xfeedfacecafebeef);
        for (uint32_t attempt = 0u; attempt < 8u; ++attempt)
        {
            wirehair_v2::ExplicitMessagePrecodeConfigForTesting config;
            if (!wirehair_v2::MakeRepairV1ExplicitConfigForTesting(
                    contract, 64u, UINT32_C(0x89abcdef), attempt,
                    false, false, config) ||
                config.Params.Seed != expected_precode[attempt] ||
                config.Packet.PeelSeed != expected_packet[attempt] ||
                config.Params.Staircase !=
                    wirehair_v2::SmallBandStaircaseCount(64u) - 1u ||
                config.Params.DenseRows != 3u ||
                config.Params.HeavyRows != contract.GF256Rows ||
                config.Params.Field !=
                    wirehair_v2::CompletionField::MixedGF256GF16 ||
                config.Packet.MixCount != 3u)
            {
                std::fprintf(stderr,
                    "repair-v1 attempt vector mismatch arm=%u a=%u\n",
                    arm, attempt);
                return false;
            }
            if (attempt == 0u)
            {
                wirehair_v2::PrecodeParams raw =
                    wirehair_v2::MakeMixedParams(
                        64u, UINT64_C(0x89abcdef));
                raw.Staircase =
                    wirehair_v2::SmallBandStaircaseCount(64u) - 1u;
                raw.DenseRows = 3u;
                wirehair_v2::PacketRowConfig raw_packet;
                raw_packet.PeelSeed =
                    wirehair_v2::RawUniformPacketPeelSeed(
                        UINT64_C(0x89abcdef));
                raw_packet.MixCount = 3u;
                if (!SameParams(raw, config.Params) ||
                    raw_packet.PeelSeed != config.Packet.PeelSeed ||
                    raw_packet.MixCount != config.Packet.MixCount)
                {
                    std::fprintf(stderr,
                        "repair-v1 attempt zero changed raw bytes arm=%u\n",
                        arm);
                    return false;
                }
            }
            previous = config;
        }
        const uint64_t previous_seed = previous.Params.Seed;
        if (wirehair_v2::MakeRepairV1ExplicitConfigForTesting(
                contract, 64u, UINT32_C(0x89abcdef), 8u,
                false, false, previous) ||
            previous.Params.Seed != previous_seed)
        {
            std::fprintf(stderr,
                "repair-v1 attempt cap was not transactional arm=%u\n",
                arm);
            return false;
        }
    }
    return true;
}

bool CheckFrozenEquationAggregate()
{
    Digest64 aggregate;
    uint32_t contract_count = 0u;
    const wirehair_v2::RepairV1Contract* contracts =
        wirehair_v2::RepairV1Contracts(contract_count);
    for (uint32_t arm = 0u; arm < contract_count; ++arm)
    {
        const wirehair_v2::RepairV1Contract& contract = contracts[arm];
        wirehair_v2::ScopedRepairV1ContractStateForTesting scope(contract);
        if (!scope.IsValid()) return false;
        aggregate.U64(contract.ProvisionalId);
        for (uint32_t K = 2u; K <= 100u; ++K)
        {
            const uint32_t root =
                UINT32_C(0x510e527f) ^
                (uint32_t)((uint64_t)K * UINT64_C(0x85ebca6b));
            for (uint32_t attempt = 0u;
                 attempt < wirehair_v2::kRepairV1AttemptCount;
                 ++attempt)
            {
                wirehair_v2::ExplicitMessagePrecodeConfigForTesting config;
                bool digest_valid = false;
                if (!wirehair_v2::MakeRepairV1ExplicitConfigForTesting(
                        contract, K, root, attempt,
                        false, false, config))
                {
                    return false;
                }
                const uint64_t equation =
                    EquationDigest(contract, config, digest_valid);
                if (!digest_valid) return false;
                aggregate.U32(K);
                aggregate.U32(attempt);
                aggregate.U64(config.Params.Seed);
                aggregate.U32(config.Packet.PeelSeed);
                aggregate.U32(config.Packet.MixCount);
                aggregate.U64(equation);
            }
        }
    }
    static const uint64_t kExpected =
        UINT64_C(0x0ee90cceafd4944c);
    if (aggregate.Get() != kExpected)
    {
        std::fprintf(stderr,
            "repair-v1 equation aggregate mismatch "
            "expected=%016llx actual=%016llx\n",
            (unsigned long long)kExpected,
            (unsigned long long)aggregate.Get());
        return false;
    }
    return true;
}

bool CheckPolicyOracle()
{
    WirehairResult probes[8];
    for (WirehairResult initial :
        {Wirehair_NeedMore, Wirehair_Error})
    {
        for (uint32_t selected = 1u; selected < 8u; ++selected)
        {
            for (WirehairResult& result : probes) {
                result = Wirehair_NeedMore;
            }
            probes[selected] = Wirehair_Success;
            if (selected + 1u < 8u) {
                probes[selected + 1u] = Wirehair_Success;
            }
            wirehair_v2::RepairV1SelectorTelemetry telemetry;
            if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
                    initial, probes, Wirehair_Success, telemetry) !=
                    Wirehair_Success ||
                telemetry.SelectedAttempt != selected ||
                telemetry.AttemptsExecuted != selected + 1u ||
                telemetry.StructuralProbeCalls != selected + 1u ||
                telemetry.RealConfigurationCalls != 2u ||
                !telemetry.Attempts[0].RealStatsAvailable ||
                !telemetry.Attempts[selected].ProbeStatsAvailable ||
                !telemetry.Attempts[selected].RealStatsAvailable ||
                telemetry.Attempts[selected + (selected + 1u < 8u ? 1u : 0u)]
                    .ProbeExecuted != (selected + 1u >= 8u) ||
                !telemetry.Committed || telemetry.CapExhausted)
            {
                std::fprintf(stderr,
                    "repair-v1 first-success oracle mismatch "
                    "initial=%d selected=%u\n", (int)initial, selected);
                return false;
            }
        }
    }

    for (WirehairResult& result : probes) {
        result = Wirehair_NeedMore;
    }
    wirehair_v2::RepairV1SelectorTelemetry telemetry;
    if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
            Wirehair_Success, probes, Wirehair_Error, telemetry) !=
            Wirehair_Success ||
        telemetry.SelectedAttempt != 0u ||
        telemetry.CallsExecuted != 1u ||
        !telemetry.Attempts[0].RealStatsAvailable ||
        !telemetry.Committed)
    {
        return false;
    }
    if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
            Wirehair_NeedMore, probes, Wirehair_Success, telemetry) !=
            Wirehair_BadPeelSeed ||
        !telemetry.CapExhausted ||
        telemetry.AttemptsExecuted != 8u ||
        telemetry.CallsExecuted != 9u ||
        !telemetry.Attempts[0].RealStatsAvailable ||
        !telemetry.Attempts[7].ProbeStatsAvailable ||
        telemetry.Committed)
    {
        return false;
    }

    probes[0] = Wirehair_Success;
    if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
            Wirehair_Error, probes, Wirehair_Success, telemetry) !=
            Wirehair_Error ||
        !telemetry.FatalAttemptZeroMismatch ||
        telemetry.Committed)
    {
        return false;
    }
    if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
            Wirehair_NeedMore, probes, Wirehair_Success, telemetry) !=
            Wirehair_Error ||
        !telemetry.FatalAttemptZeroMismatch ||
        telemetry.Committed)
    {
        return false;
    }
    probes[0] = Wirehair_OOM;
    if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
            Wirehair_NeedMore, probes, Wirehair_Success, telemetry) !=
            Wirehair_OOM ||
        !telemetry.Oom ||
        !telemetry.Attempts[0].RealStatsAvailable ||
        telemetry.Attempts[0].ProbeStatsAvailable ||
        telemetry.Committed)
    {
        return false;
    }
    probes[0] = Wirehair_NeedMore;
    probes[2] = Wirehair_Success;
    if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
            Wirehair_NeedMore, probes, Wirehair_OOM, telemetry) !=
            Wirehair_OOM ||
        telemetry.SelectedAttempt != 2u ||
        !telemetry.Attempts[2].ProbeStatsAvailable ||
        telemetry.Attempts[2].RealStatsAvailable ||
        !telemetry.Oom || telemetry.Committed)
    {
        return false;
    }

    if (wirehair_v2::EvaluateRepairV1AccountingOomForTesting(
            wirehair_v2::RepairV1AccountingOomStageForTesting::
                AttemptZeroBeforeSystemBuild,
            telemetry) != Wirehair_OOM ||
        telemetry.CallsExecuted != 1u ||
        telemetry.RealConfigurationCalls != 1u ||
        telemetry.StructuralProbeCalls != 0u ||
        telemetry.AttemptsExecuted != 1u ||
        !telemetry.Attempts[0].RealExecuted ||
        telemetry.Attempts[0].RealResult != Wirehair_OOM ||
        telemetry.Attempts[0].RealStatsAvailable ||
        telemetry.Attempts[0].ProbeExecuted ||
        !telemetry.Oom || telemetry.Committed)
    {
        std::fprintf(stderr,
            "repair-v1 pre-build OOM accounting mismatch\n");
        return false;
    }
    if (wirehair_v2::EvaluateRepairV1AccountingOomForTesting(
            wirehair_v2::RepairV1AccountingOomStageForTesting::
                ProbeDuringSystemBuild,
            telemetry) != Wirehair_OOM ||
        telemetry.CallsExecuted != 2u ||
        telemetry.RealConfigurationCalls != 1u ||
        telemetry.StructuralProbeCalls != 1u ||
        telemetry.AttemptsExecuted != 1u ||
        !telemetry.Attempts[0].RealStatsAvailable ||
        !telemetry.Attempts[0].ProbeExecuted ||
        telemetry.Attempts[0].ProbeResult != Wirehair_OOM ||
        telemetry.Attempts[0].ProbeStatsAvailable ||
        !telemetry.Oom || telemetry.Committed)
    {
        std::fprintf(stderr,
            "repair-v1 probe-build OOM accounting mismatch\n");
        return false;
    }
    if (wirehair_v2::EvaluateRepairV1AccountingOomForTesting(
            wirehair_v2::RepairV1AccountingOomStageForTesting::
                ProbeDuringSolve,
            telemetry) != Wirehair_OOM ||
        telemetry.CallsExecuted != 2u ||
        telemetry.RealConfigurationCalls != 1u ||
        telemetry.StructuralProbeCalls != 1u ||
        telemetry.AttemptsExecuted != 1u ||
        !telemetry.Attempts[0].RealStatsAvailable ||
        !telemetry.Attempts[0].ProbeStatsAvailable ||
        telemetry.Attempts[0].ProbeStats.PacketRows != 1u ||
        telemetry.Attempts[0].ProbeResult != Wirehair_OOM ||
        !telemetry.Oom || telemetry.Committed)
    {
        std::fprintf(stderr,
            "repair-v1 probe-solve OOM accounting mismatch\n");
        return false;
    }

    for (int value = 0; value < (int)WirehairResult_Count; ++value)
    {
        const WirehairResult terminal = (WirehairResult)value;
        if (terminal == Wirehair_Success ||
            terminal == Wirehair_NeedMore ||
            terminal == Wirehair_Error)
        {
            continue;
        }
        if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
                terminal, probes, Wirehair_Success, telemetry) != terminal ||
            telemetry.CallsExecuted != 1u ||
            telemetry.RealConfigurationCalls != 1u ||
            telemetry.StructuralProbeCalls != 0u ||
            telemetry.Attempts[0].RealStatsAvailable ||
            telemetry.Committed)
        {
            std::fprintf(stderr,
                "repair-v1 initial terminal result mismatch result=%d\n",
                value);
            return false;
        }
    }

    for (int value = 0; value < (int)WirehairResult_Count; ++value)
    {
        const WirehairResult terminal = (WirehairResult)value;
        if (terminal == Wirehair_Success ||
            terminal == Wirehair_NeedMore)
        {
            continue;
        }
        for (WirehairResult& probe : probes) {
            probe = Wirehair_NeedMore;
        }
        probes[0] = terminal;
        const bool stats_available = terminal == Wirehair_Error;
        if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
                Wirehair_NeedMore, probes, Wirehair_Success, telemetry) !=
                terminal ||
            telemetry.CallsExecuted != 2u ||
            telemetry.RealConfigurationCalls != 1u ||
            telemetry.StructuralProbeCalls != 1u ||
            telemetry.Attempts[0].ProbeStatsAvailable != stats_available ||
            telemetry.Committed)
        {
            std::fprintf(stderr,
                "repair-v1 probe terminal result mismatch result=%d\n",
                value);
            return false;
        }
    }

    for (WirehairResult& probe : probes) {
        probe = Wirehair_NeedMore;
    }
    probes[1] = Wirehair_Success;
    for (int value = 0; value < (int)WirehairResult_Count; ++value)
    {
        const WirehairResult selected_result = (WirehairResult)value;
        const bool stats_available =
            selected_result == Wirehair_Success ||
            selected_result == Wirehair_NeedMore ||
            selected_result == Wirehair_Error;
        if (wirehair_v2::EvaluateRepairV1PolicyOracleForTesting(
                Wirehair_NeedMore, probes, selected_result, telemetry) !=
                selected_result ||
            telemetry.SelectedAttempt != 1u ||
            telemetry.CallsExecuted != 4u ||
            telemetry.RealConfigurationCalls != 2u ||
            telemetry.StructuralProbeCalls != 2u ||
            telemetry.Attempts[1].RealStatsAvailable != stats_available ||
            telemetry.Committed !=
                (selected_result == Wirehair_Success))
        {
            std::fprintf(stderr,
                "repair-v1 selected result mismatch result=%d\n", value);
            return false;
        }
    }
    return true;
}

struct IndependenceBaseline
{
    WirehairResult Result = Wirehair_InvalidInput;
    uint32_t Attempt = wirehair_v2::kRepairV1NoAttempt;
    uint64_t Equation = 0u;
};

bool CheckAllKPayloadWidthIndependence()
{
    uint32_t contract_count = 0u;
    const wirehair_v2::RepairV1Contract* contracts =
        wirehair_v2::RepairV1Contracts(contract_count);
    const uint32_t widths[] = {
        2u, 6u, 32u, 64u, 256u, 1280u, 4096u
    };
    for (uint32_t arm = 0u; arm < contract_count; ++arm)
    {
        const wirehair_v2::RepairV1Contract& contract = contracts[arm];
        for (uint32_t K = 2u; K <= 100u; ++K)
        {
            const uint32_t root =
                UINT32_C(0x6a09e667) ^
                (uint32_t)((uint64_t)K * UINT64_C(0x9e3779b9));
            IndependenceBaseline baseline;
            bool have_baseline = false;
            for (uint32_t width : widths)
            {
                for (uint32_t variant = 0u; variant < 2u; ++variant)
                {
                    for (uint32_t partial = 0u;
                         partial < 3u;
                         ++partial)
                    {
                        std::vector<uint8_t> message;
                        FillMessage(message, K, width, variant);
                        if (partial == 1u) {
                            message.resize(message.size() - 1u);
                        }
                        else if (partial == 2u) {
                            message.resize((size_t)(K - 1u) * width + 1u);
                        }
                        wirehair_v2::MessagePrecodeEncoder encoder;
                        wirehair_v2::RepairV1SelectorTelemetry telemetry;
                        const WirehairResult result =
                            encoder.InitializeRepairV1ResultForTesting(
                                message.data(), message.size(), width,
                                contract, root, &telemetry);
                        if ((result != Wirehair_Success &&
                             result != Wirehair_BadPeelSeed) ||
                            (result == Wirehair_Success &&
                             (!telemetry.Committed ||
                              telemetry.SelectedAttempt >=
                                wirehair_v2::kRepairV1AttemptCount)) ||
                            (result == Wirehair_BadPeelSeed &&
                             (!telemetry.CapExhausted ||
                              telemetry.Committed ||
                              telemetry.SelectedAttempt !=
                                wirehair_v2::kRepairV1NoAttempt)))
                        {
                            std::fprintf(stderr,
                                "repair-v1 all-K terminal mismatch "
                                "arm=%u K=%u result=%d\n",
                                arm, K, (int)result);
                            return false;
                        }
                        uint64_t equation = 0u;
                        if (result == Wirehair_Success)
                        {
                            wirehair_v2::ScopedRepairV1ContractStateForTesting
                                scope(contract);
                            wirehair_v2::
                                ExplicitMessagePrecodeConfigForTesting config;
                            bool digest_valid = false;
                            if (!scope.IsValid() ||
                                !wirehair_v2::
                                    MakeRepairV1ExplicitConfigForTesting(
                                        contract, K, root,
                                        telemetry.SelectedAttempt,
                                        false, false, config))
                            {
                                return false;
                            }
                            equation = EquationDigest(
                                contract, config, digest_valid);
                            if (!digest_valid ||
                                encoder.Profile().V2SeedAttempt !=
                                    telemetry.SelectedAttempt ||
                                !telemetry.Attempts[
                                    telemetry.SelectedAttempt].
                                        RealStatsAvailable ||
                                telemetry.Attempts[
                                    telemetry.SelectedAttempt].
                                        RealStats.PacketSeedAttempt !=
                                    telemetry.SelectedAttempt ||
                                encoder.
                                    ProvisionalRepairContractIdForTesting() !=
                                    contract.ProvisionalId ||
                                encoder.DiagnosticIdentityForTesting() !=
                                    wirehair_v2::
                                    MessagePrecodeDiagnosticIdentityForTesting::
                                        ProvisionalRepairContract)
                            {
                                return false;
                            }
                        }
                        if (!have_baseline)
                        {
                            baseline.Result = result;
                            baseline.Attempt = telemetry.SelectedAttempt;
                            baseline.Equation = equation;
                            have_baseline = true;
                        }
                        else if (result != baseline.Result ||
                                 telemetry.SelectedAttempt != baseline.Attempt ||
                                 equation != baseline.Equation)
                        {
                            std::fprintf(stderr,
                                "repair-v1 payload/width dependence "
                                "arm=%u K=%u width=%u variant=%u partial=%u\n",
                                arm, K, width, variant, partial);
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

bool CheckSemanticBridges()
{
    const uint32_t Ks[] = {2u, 3u, 10u, 64u, 100u};
    uint32_t contract_count = 0u;
    const wirehair_v2::RepairV1Contract* contracts =
        wirehair_v2::RepairV1Contracts(contract_count);
    for (uint32_t arm = 0u; arm < contract_count; ++arm)
    {
        const wirehair_v2::RepairV1Contract& contract = contracts[arm];
        for (uint32_t K : Ks)
        {
            const uint32_t width = 6u;
            const uint32_t root =
                UINT32_C(0xbb67ae85) ^ K * UINT32_C(0x85ebca6b);
            std::vector<uint8_t> message;
            FillMessage(message, K, width, 9u);
            for (uint32_t attempt = 0u; attempt < 8u; ++attempt)
            {
                wirehair_v2::ScopedRepairV1ContractStateForTesting scope(
                    contract);
                wirehair_v2::ExplicitMessagePrecodeConfigForTesting config;
                if (!scope.IsValid() ||
                    !wirehair_v2::MakeRepairV1ExplicitConfigForTesting(
                        contract, K, root, attempt,
                        false, false, config))
                {
                    return false;
                }
                wirehair_v2::MessagePrecodeEncoder explicit_encoder;
                wirehair_v2::MessagePrecodeEncoder contract_encoder;
                const WirehairResult explicit_result =
                    explicit_encoder.InitializeExplicitResultForTesting(
                        message.data(), message.size(), width, config);
                const WirehairResult contract_result =
                    contract_encoder.
                        InitializeRepairV1SelectedResultForTesting(
                            message.data(), message.size(), width,
                            contract, root, attempt);
                if (explicit_result != contract_result) {
                    return false;
                }
                if (contract_result != Wirehair_Success) {
                    continue;
                }
                const size_t intermediate_bytes =
                    (size_t)(K + config.Params.Staircase +
                        config.Params.DenseRows + config.Params.HeavyRows) *
                    width;
                if (!SameParams(
                        explicit_encoder.BlockEncoder().System().Params,
                        contract_encoder.BlockEncoder().System().Params) ||
                    std::memcmp(
                        explicit_encoder.IntermediateBlocks(),
                        contract_encoder.IntermediateBlocks(),
                        intermediate_bytes) != 0)
                {
                    return false;
                }
                for (uint32_t block_id = 0u;
                     block_id < K + 5u;
                     ++block_id)
                {
                    uint8_t explicit_block[6] = {};
                    uint8_t contract_block[6] = {};
                    uint32_t explicit_bytes = 0u;
                    uint32_t contract_bytes = 0u;
                    if (explicit_encoder.EncodeResult(
                            block_id, explicit_block, width,
                            &explicit_bytes) != Wirehair_Success ||
                        contract_encoder.EncodeResult(
                            block_id, contract_block, width,
                            &contract_bytes) != Wirehair_Success ||
                        explicit_bytes != contract_bytes ||
                        std::memcmp(
                            explicit_block, contract_block,
                            contract_bytes) != 0)
                    {
                        return false;
                    }
                }

                wirehair_v2::MessagePrecodeDecoder decoder;
                if (decoder.InitializeRepairV1ResultForTesting(
                        message.size(), width, contract, root, attempt) !=
                        Wirehair_Success ||
                    decoder.PacketSeedAttempt() != attempt ||
                    decoder.ProvisionalRepairContractIdForTesting() !=
                        contract.ProvisionalId)
                {
                    return false;
                }
                WirehairResult decode_result = Wirehair_NeedMore;
                for (uint32_t block_id = 0u; block_id < K; ++block_id)
                {
                    uint8_t block[6] = {};
                    uint32_t bytes = 0u;
                    if (contract_encoder.EncodeResult(
                            block_id, block, width, &bytes) !=
                            Wirehair_Success)
                    {
                        return false;
                    }
                    decode_result =
                        decoder.DecodeResult(block_id, block, bytes);
                }
                std::vector<uint8_t> recovered(message.size());
                if (decode_result != Wirehair_Success ||
                    decoder.RecoverResult(
                        recovered.data(), recovered.size()) !=
                        Wirehair_Success ||
                    recovered != message)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

bool CheckEndpointScopeLifecycle()
{
    if (wirehair_v2::ActiveMixedCoefficientPeriod() !=
            wirehair_v2::kMixedCoefficientPeriod ||
        wirehair_v2::ActiveMixedGF256Rows() !=
            wirehair_v2::kMixedGF256Rows ||
        wirehair_v2::ActiveMixedGF16Rows() !=
            wirehair_v2::kMixedGF16Rows)
    {
        return false;
    }

    uint32_t contract_count = 0u;
    const wirehair_v2::RepairV1Contract* contracts =
        wirehair_v2::RepairV1Contracts(contract_count);
    const uint32_t K = 24u;
    const uint32_t width = 6u;
    for (uint32_t arm = 0u; arm < contract_count; ++arm)
    {
        const wirehair_v2::RepairV1Contract& contract = contracts[arm];
        const uint32_t root =
            UINT32_C(0xbb67ae85) ^ K * UINT32_C(0x85ebca6b);
        std::vector<uint8_t> message;
        FillMessage(message, K, width, arm + 17u);
        wirehair_v2::MessagePrecodeEncoder encoder;
        wirehair_v2::RepairV1SelectorTelemetry telemetry;
        if (encoder.InitializeRepairV1ResultForTesting(
                message.data(), message.size(), width,
                contract, root, &telemetry) != Wirehair_Success)
        {
            std::fprintf(stderr,
                "repair-v1 endpoint lifecycle seed failed arm=%u\n", arm);
            return false;
        }
        wirehair_v2::MessagePrecodeDecoder decoder;
        if (decoder.InitializeRepairV1ResultForTesting(
                message.size(), width, contract, root,
                telemetry.SelectedAttempt) != Wirehair_Success ||
            !encoder.PinnedEquationStateForTesting().IsPinned() ||
            !decoder.PinnedEquationStateForTesting().IsPinned() ||
            decoder.DiagnosticIdentityForTesting() !=
                wirehair_v2::
                MessagePrecodeDiagnosticIdentityForTesting::
                    ProvisionalRepairContract)
        {
            return false;
        }

        uint8_t drift_block[6] = {};
        uint32_t drift_bytes = 0u;
        if (encoder.EncodeResult(
                0u, drift_block, width, &drift_bytes) != Wirehair_Success)
        {
            return false;
        }
        const bool set_drift_pmf =
            wirehair_v2::SetPeelDegreesForTesting(
                std::vector<double>{1.0, 1.0});
        if (!set_drift_pmf)
        {
            wirehair_v2::ClearPeelDegreesForTesting();
            return false;
        }
        uint8_t ignored_block[6] = {};
        uint32_t ignored_bytes = 0u;
        const WirehairResult drift_encode = encoder.EncodeResult(
            0u, ignored_block, width, &ignored_bytes);
        const WirehairResult drift_decode = decoder.DecodeResult(
            0u, drift_block, drift_bytes);
        wirehair_v2::ClearPeelDegreesForTesting();
        if (drift_encode != Wirehair_Error ||
            drift_decode != Wirehair_Error ||
            decoder.ReceivedCount() != 0u)
        {
            std::fprintf(stderr,
                "repair-v1 endpoint accepted equation-state drift arm=%u\n",
                arm);
            return false;
        }

        WirehairResult decode_result = Wirehair_NeedMore;
        for (uint32_t block_id = 0u; block_id < K; ++block_id)
        {
            uint8_t block[6] = {};
            uint32_t bytes = 0u;
            if (encoder.EncodeResult(
                    block_id, block, width, &bytes) != Wirehair_Success ||
                wirehair_v2::ActiveMixedGF256Rows() !=
                    wirehair_v2::kMixedGF256Rows ||
                wirehair_v2::ActiveMixedGF16Rows() !=
                    wirehair_v2::kMixedGF16Rows)
            {
                return false;
            }
            decode_result =
                decoder.DecodeResult(block_id, block, bytes);
            if (wirehair_v2::ActiveMixedGF256Rows() !=
                    wirehair_v2::kMixedGF256Rows ||
                wirehair_v2::ActiveMixedGF16Rows() !=
                    wirehair_v2::kMixedGF16Rows)
            {
                return false;
            }
        }
        std::vector<uint8_t> recovered(message.size(), uint8_t{0xa5});
        if (decode_result != Wirehair_Success) return false;
        const bool set_recover_drift_pmf =
            wirehair_v2::SetPeelDegreesForTesting(
                std::vector<double>{1.0, 1.0});
        if (!set_recover_drift_pmf)
        {
            wirehair_v2::ClearPeelDegreesForTesting();
            return false;
        }
        const WirehairResult drift_recover = decoder.RecoverResult(
            recovered.data(), recovered.size());
        wirehair_v2::ClearPeelDegreesForTesting();
        if (drift_recover != Wirehair_Error ||
            recovered !=
                std::vector<uint8_t>(message.size(), uint8_t{0xa5}) ||
            decoder.RecoverResult(
                recovered.data(), recovered.size()) != Wirehair_Success ||
            recovered != message ||
            wirehair_v2::ActiveMixedGF256Rows() !=
                wirehair_v2::kMixedGF256Rows ||
            wirehair_v2::ActiveMixedGF16Rows() !=
                wirehair_v2::kMixedGF16Rows)
        {
            std::fprintf(stderr,
                "repair-v1 endpoint scope lifecycle failed arm=%u\n", arm);
            return false;
        }
    }
    return true;
}

bool CheckAllRestorableGeometryStates()
{
    const wirehair_v2::RepairV1Contract* contract =
        wirehair_v2::FindRepairV1Contract(
            wirehair_v2::RepairV1Arm::Pure8S0M1D3);
    if (!contract) return false;

    for (uint32_t gf256_rows = 1u;
         gf256_rows <= wirehair_v2::kMixedGF256Rows;
         ++gf256_rows)
    {
        for (uint32_t gf16_rows = 0u; gf16_rows <= 4u; ++gf16_rows)
        {
            if (!wirehair_v2::SetMixedCoefficientPeriodForTesting(
                    wirehair_v2::kMixedCoefficientPeriod) ||
                !wirehair_v2::SetMixedGF256RowsForTesting(
                    wirehair_v2::kMixedGF256Rows) ||
                !wirehair_v2::SetMixedGF16RowsForTesting(gf16_rows) ||
                !wirehair_v2::SetMixedGF256RowsForTesting(gf256_rows) ||
                !wirehair_v2::SetMixedCoefficientPeriodForTesting(
                    gf256_rows + gf16_rows))
            {
                return false;
            }
            {
                wirehair_v2::ScopedRepairV1ContractStateForTesting scope(
                    *contract);
                if (!scope.IsValid() ||
                    wirehair_v2::ActiveMixedCoefficientPeriod() !=
                        wirehair_v2::kMixedCoefficientPeriod ||
                    wirehair_v2::ActiveMixedGF256Rows() !=
                        contract->GF256Rows ||
                    wirehair_v2::ActiveMixedGF16Rows() !=
                        contract->GF16Rows)
                {
                    return false;
                }
            }
            if (wirehair_v2::ActiveMixedCoefficientPeriod() !=
                    gf256_rows + gf16_rows ||
                wirehair_v2::ActiveMixedGF256Rows() != gf256_rows ||
                wirehair_v2::ActiveMixedGF16Rows() != gf16_rows)
            {
                std::fprintf(stderr,
                    "repair-v1 exhaustive geometry restoration failed "
                    "gf256=%u gf16=%u\n", gf256_rows, gf16_rows);
                return false;
            }
        }
    }
    return
        wirehair_v2::SetMixedCoefficientPeriodForTesting(
            wirehair_v2::kMixedCoefficientPeriod) &&
        wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256Rows) &&
        wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16Rows);
}

bool CheckTransactionsAndAmbientIsolation()
{
    const wirehair_v2::RepairV1Contract* contract =
        wirehair_v2::FindRepairV1Contract(
            wirehair_v2::RepairV1Arm::Pure8S0M1D3);
    if (!contract) return false;
    const uint32_t K = 16u;
    const uint32_t width = 6u;
    const uint32_t root = UINT32_C(0x3c6ef372);
    std::vector<uint8_t> message;
    FillMessage(message, K, width, 4u);
    wirehair_v2::MessagePrecodeEncoder encoder;
    wirehair_v2::RepairV1SelectorTelemetry telemetry;
    if (encoder.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry) != Wirehair_Success)
    {
        return false;
    }
    const uint32_t old_attempt = encoder.Profile().V2SeedAttempt;
    const uint64_t old_id =
        encoder.ProvisionalRepairContractIdForTesting();
    const size_t old_bytes =
        (size_t)(K + encoder.Profile().V2StaircaseCount +
            encoder.Profile().V2DenseRowCount +
            encoder.Profile().V2HeavyRowCount) * width;
    std::vector<uint8_t> old_intermediate(
        encoder.IntermediateBlocks(),
        encoder.IntermediateBlocks() + old_bytes);

    if (!wirehair_v2::SetPeelDegreesForTesting(
            std::vector<double>{1.0, 1.0}))
    {
        wirehair_v2::ClearPeelDegreesForTesting();
        return false;
    }
    wirehair_v2::MessagePrecodeEncoder polluted_encoder;
    const WirehairResult pmf_polluted =
        polluted_encoder.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    wirehair_v2::ClearPeelDegreesForTesting();
    if (pmf_polluted != Wirehair_InvalidInput ||
        telemetry.CallsExecuted != 0u)
    {
        std::fprintf(stderr,
            "repair-v1 accepted packet-PMF pollution\n");
        return false;
    }
    if (wirehair_v2::SetPeelDegreesForTesting(
            std::vector<double>{0.0, 0.0}))
    {
        wirehair_v2::ClearPeelDegreesForTesting();
        return false;
    }
    const WirehairResult invalid_pmf =
        polluted_encoder.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    wirehair_v2::ClearPeelDegreesForTesting();
    if (invalid_pmf != Wirehair_InvalidInput ||
        telemetry.CallsExecuted != 0u)
    {
        std::fprintf(stderr,
            "repair-v1 accepted malformed packet-PMF pollution\n");
        return false;
    }
    wirehair_v2::SetStaircaseDegreeScaleForTesting(0.0);
    const WirehairResult scale_polluted =
        polluted_encoder.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    wirehair_v2::ClearStaircaseDegreeScaleForTesting();
    if (scale_polluted != Wirehair_InvalidInput ||
        telemetry.CallsExecuted != 0u)
    {
        std::fprintf(stderr,
            "repair-v1 accepted staircase-scale pollution\n");
        return false;
    }
    wirehair_v2::SetStaircaseDegreeScaleForTesting(-2.0);
    const WirehairResult invalid_scale =
        polluted_encoder.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    wirehair_v2::ClearStaircaseDegreeScaleForTesting();
    if (invalid_scale != Wirehair_InvalidInput ||
        telemetry.CallsExecuted != 0u)
    {
        std::fprintf(stderr,
            "repair-v1 accepted malformed staircase-scale pollution\n");
        return false;
    }

    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult oom = encoder.InitializeRepairV1ResultForTesting(
        message.data(), message.size(), width,
        *contract, root + 1u, &telemetry);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (oom != Wirehair_OOM || !telemetry.Oom ||
        telemetry.CallsExecuted != 1u ||
        telemetry.RealConfigurationCalls != 1u ||
        telemetry.StructuralProbeCalls != 0u ||
        telemetry.AttemptsExecuted != 1u ||
        !telemetry.Attempts[0].RealExecuted ||
        telemetry.Attempts[0].RealResult != Wirehair_OOM ||
        telemetry.Attempts[0].RealStatsAvailable ||
        telemetry.Attempts[0].ProbeExecuted ||
        telemetry.Committed ||
        encoder.Profile().V2SeedAttempt != old_attempt ||
        encoder.ProvisionalRepairContractIdForTesting() != old_id ||
        std::memcmp(
            encoder.IntermediateBlocks(),
            old_intermediate.data(), old_bytes) != 0)
    {
        std::fprintf(stderr, "repair-v1 encoder OOM rollback failed\n");
        return false;
    }
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult selected_oom =
        encoder.InitializeRepairV1SelectedResultForTesting(
            message.data(), message.size(), width,
            *contract, root + 2u, old_attempt);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (selected_oom != Wirehair_OOM ||
        encoder.Profile().V2SeedAttempt != old_attempt ||
        encoder.ProvisionalRepairContractIdForTesting() != old_id ||
        std::memcmp(
            encoder.IntermediateBlocks(),
            old_intermediate.data(), old_bytes) != 0)
    {
        std::fprintf(stderr,
            "repair-v1 forced-selected OOM rollback failed\n");
        return false;
    }
    wirehair_v2::RepairV1Contract forged = *contract;
    forged.GF256Rows = 9u;
    if (encoder.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            forged, root, &telemetry) != Wirehair_InvalidInput ||
        encoder.ProvisionalRepairContractIdForTesting() != old_id ||
        encoder.InitializeRepairV1SelectedResultForTesting(
            message.data(), message.size(), width,
            *contract, root, 8u) != Wirehair_InvalidInput ||
        encoder.ProvisionalRepairContractIdForTesting() != old_id)
    {
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder decoder;
    if (decoder.InitializeRepairV1ResultForTesting(
            message.size(), width, *contract, root,
            old_attempt) != Wirehair_Success)
    {
        return false;
    }
    const uint64_t decoder_id =
        decoder.ProvisionalRepairContractIdForTesting();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
    const WirehairResult decoder_oom =
        decoder.InitializeRepairV1ResultForTesting(
            message.size(), width, *contract, root + 1u,
            old_attempt);
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    if (decoder_oom != Wirehair_OOM ||
        decoder.ProvisionalRepairContractIdForTesting() != decoder_id ||
        decoder.InitializeRepairV1ResultForTesting(
            message.size(), width, *contract, root, 8u) !=
            Wirehair_InvalidInput ||
        decoder.ProvisionalRepairContractIdForTesting() != decoder_id)
    {
        return false;
    }

    if (!wirehair_v2::SetMixedGF16RowsForTesting(0u) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(9u) ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(31u))
    {
        return false;
    }
    wirehair_v2::MessagePrecodeEncoder isolated;
    const WirehairResult isolated_result =
        isolated.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    if (wirehair_v2::ActiveMixedGF256Rows() != 9u ||
        wirehair_v2::ActiveMixedGF16Rows() != 0u ||
        wirehair_v2::ActiveMixedCoefficientPeriod() != 31u ||
        isolated_result != Wirehair_Success ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(
            wirehair_v2::kMixedCoefficientPeriod) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(10u) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(2u))
    {
        std::fprintf(stderr, "repair-v1 row-state restoration failed\n");
        return false;
    }

    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(3u)) {
        return false;
    }
    const WirehairResult polluted =
        isolated.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    (void)wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u);
    if (polluted != Wirehair_InvalidInput) {
        return false;
    }
    if (!wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX))
    {
        return false;
    }
    const WirehairResult geometry_polluted =
        isolated.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    const bool geometry_preserved =
        wirehair_v2::ActiveMixedCoefficientGeometry() ==
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX;
    (void)wirehair_v2::SetMixedCoefficientGeometryForTesting(
        wirehair_v2::MixedCoefficientGeometry::FrozenPowerX);
    if (geometry_polluted != Wirehair_InvalidInput ||
        !geometry_preserved)
    {
        return false;
    }
    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    const WirehairResult tracking =
        isolated.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), width,
            *contract, root, &telemetry);
    wirehair_v2::ClearMixedBandTrackingXForTesting();
    return tracking == Wirehair_InvalidInput;
}

bool CheckConcurrentIsolation()
{
    const wirehair_v2::RepairV1Contract* pure8 =
        wirehair_v2::FindRepairV1Contract(
            wirehair_v2::RepairV1Arm::Pure8S0M1D3);
    const wirehair_v2::RepairV1Contract* pure9 =
        wirehair_v2::FindRepairV1Contract(
            wirehair_v2::RepairV1Arm::Pure9S0M1D3);
    const wirehair_v2::V2EquationContract* dispatch =
        wirehair_v2::FindV2EquationContract("dispatch-v1");
    if (!pure8 || !pure9 || !dispatch) return false;
    std::atomic<bool> ok(true);
    const auto arm_worker = [&ok](
        const wirehair_v2::RepairV1Contract* contract)
    {
        for (uint32_t K = 2u; K <= 100u && ok.load(); K += 7u)
        {
            std::vector<uint8_t> message;
            FillMessage(message, K, 6u, contract->GF256Rows);
            wirehair_v2::MessagePrecodeEncoder encoder;
            wirehair_v2::RepairV1SelectorTelemetry telemetry;
            const WirehairResult result =
                encoder.InitializeRepairV1ResultForTesting(
                    message.data(), message.size(), 6u, *contract,
                    UINT32_C(0xa54ff53a) ^ K, &telemetry);
            if (result != Wirehair_Success &&
                result != Wirehair_BadPeelSeed)
            {
                ok.store(false);
                return;
            }
            if (result == Wirehair_Success &&
                (encoder.ProvisionalRepairContractIdForTesting() !=
                    contract->ProvisionalId ||
                 encoder.Profile().V2HeavyRowCount !=
                    contract->GF256Rows))
            {
                ok.store(false);
                return;
            }
        }
    };
    const auto dispatch_worker = [&ok, dispatch]()
    {
        for (uint32_t K = 2u; K <= 100u && ok.load(); ++K)
        {
            wirehair_v2::SeedProfile profile;
            if (!wirehair_v2::MakeRawContractProfile(
                    *dispatch, K, 6u,
                    UINT32_C(0x510e527f) ^ K, profile) ||
                profile.V2HeavyRowCount !=
                    wirehair_v2::kMixedGF256Rows +
                        wirehair_v2::kMixedGF16Rows ||
                profile.V2SeedAttempt != 0u)
            {
                ok.store(false);
                return;
            }
        }
    };
    std::thread a(arm_worker, pure8);
    std::thread b(arm_worker, pure9);
    std::thread c(dispatch_worker);
    a.join();
    b.join();
    c.join();
    return ok.load() &&
        wirehair_v2::ActiveMixedGF256Rows() ==
            wirehair_v2::kMixedGF256Rows &&
        wirehair_v2::ActiveMixedGF16Rows() ==
        wirehair_v2::kMixedGF16Rows;
}

bool CheckEnvironmentPollutionRejected()
{
    const wirehair_v2::RepairV1Contract* contract =
        wirehair_v2::FindRepairV1Contract(
            wirehair_v2::RepairV1Arm::Pure8S0M1D3);
    if (!contract) return false;
    std::vector<uint8_t> message;
    FillMessage(message, 16u, 6u, 33u);
    wirehair_v2::MessagePrecodeEncoder encoder;
    wirehair_v2::RepairV1SelectorTelemetry telemetry;
    return encoder.InitializeRepairV1ResultForTesting(
            message.data(), message.size(), 6u, *contract,
            UINT32_C(0x3c6ef372), &telemetry) ==
            Wirehair_InvalidInput &&
        telemetry.CallsExecuted == 0u &&
        !telemetry.Committed;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 &&
        std::strcmp(argv[1], "--reject-environment") == 0)
    {
        if (!CheckEnvironmentPollutionRejected()) {
            return 1;
        }
        std::puts("repair-v1 environment rejection: PASS");
        return 0;
    }
    if (argc != 1) return 2;
    if (!CheckContractsAndVectors() ||
        !CheckFrozenEquationAggregate() ||
        !CheckPolicyOracle() ||
        !CheckAllKPayloadWidthIndependence() ||
        !CheckSemanticBridges() ||
        !CheckEndpointScopeLifecycle() ||
        !CheckAllRestorableGeometryStates() ||
        !CheckTransactionsAndAmbientIsolation() ||
        !CheckConcurrentIsolation())
    {
        return 1;
    }
    std::puts(
        "repair-v1 selector contract/oracle/all-K/isolation: PASS");
    return 0;
}
