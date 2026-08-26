#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"
#include "Wh2NativePanel.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#error "systematic emission screen requires an exact source commit receipt"
#endif

namespace {

using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativeEncoderStorageReceipt;
using wirehair_wh2_bench::NativePanelArm;
using wirehair_wh2_bench::NativePanelDisposition;
using wirehair_wh2_bench::NativePanelInvocation;
using wirehair_wh2_bench::NativePanelInvocationResult;
using wirehair_wh2_bench::NativePanelOrder;
using wirehair_wh2_bench::NativePanelResult;
using wirehair_wh2_bench::NativePanelSide;
using wirehair_wh2_bench::NativePanelStatus;
using wirehair_wh2_bench::NativeSystematicEmission;
using wirehair_wh2_bench::ResolvedNativeWh2Configuration;

static const char kConfigSchema[] =
    "wirehair.wh2.systematic-emission-screen.config.v1";
static const char kPanelSchema[] =
    "wirehair.wh2.systematic-emission-screen.panel.v1";
static const char kTerminalSchema[] =
    "wirehair.wh2.systematic-emission-screen.terminal.v1";
static const char kDescriptorSchema[] =
    "wirehair.wh2.systematic-emission-arm.v1";
static const char kSourceStorage[] = "native-arm-owned-kxb-v1";
static const char kSourceAcquisition[] =
    "fixture-copy-before-clock-move-v1";
static const char kRepairSourceAcquisition[] =
    "native-arm-copy-before-repair-clock-v1";
static const uint32_t kPanelReplicates = 12u;
static const uint32_t kInvocationBudget = 65536u;
static const uint32_t kInternalDeadlineSeconds = 115u;

struct CellShape
{
    uint32_t K;
    uint32_t BlockBytes;
};

static const CellShape kCellShapes[] = {
    { 32u, 64u },
    { 1000u, 64u },
    { 2048u, 1280u }
};
static const char* const kExpectedConfigurationSha256[] = {
    "51ddbe10ffbdc0d0dbe534904936c43b0425fe099b94c9ac60a16602a0a29135",
    "f6ae0eeca51c96a4e0fa6ddf01558167fa27f6bf174dfc70ae93d64448905809",
    "546bd612c4fb21e66e4959c8d2c5572cb9dd7c8dfa107e9c1fa60353b0606541"
};
static const char* const kExpectedSourceSha256[] = {
    "c837f2e189651cb420f6d2de1af12d8a5808cc2a9221f7ac1c2a5a12e3864e08",
    "931232d4fbe6172bd35eb9ba8a3090a39403ce98028fc630c3da833c5e82e492",
    "5dddcad4ab034d18b0fa11fd096935dd3f57c2b10d4d0dac1e533c2efdeb3402"
};
static const char* const kExpectedRepairOracleSha256[] = {
    "323a6d8e0e9462ce71f470a0a2e00994b0ea38652912c7440c777dcebed504ba",
    "f7262b5a3c47b416306df10ce56a7933450098587a8980b97a99d79ebc57a6f1",
    "28b358a4dbabea124cd61f239673f078f5cbcdf5c6d92f82b0eaac56f527c154"
};

struct ScreenCell
{
    CellShape Shape;
    std::vector<uint8_t> Source;
    std::vector<uint8_t> ExpectedRepairs;
    std::string EquationConfigurationJson;
    std::string EquationConfigurationSha256;
    std::string RepairOracleSha256;
    std::string SourceSha256;
    bool HighRepairControlsVerified = false;
};

const char* SideName(NativePanelSide side)
{
    return side == NativePanelSide::Left ? "left" : "right";
}

const char* OrderName(NativePanelOrder order)
{
    return order == NativePanelOrder::ABBA ? "ABBA" : "BAAB";
}

bool IsKnownScope(const char* scope)
{
    return scope &&
        (std::strcmp(scope, "encoder-systematic") == 0 ||
         std::strcmp(scope, "repair-evaluate") == 0);
}

uint32_t InvocationsPerSlot(uint32_t K, uint32_t replicate)
{
    const uint32_t total = (kInvocationBudget + K - 1u) / K;
    const uint32_t quotient = total / kPanelReplicates;
    const uint32_t remainder = total % kPanelReplicates;
    return quotient + (replicate < remainder ? 1u : 0u);
}

std::string PanelKeySha256(
    const CellShape& shape,
    const char* comparison,
    const char* scope)
{
    if (!comparison || !*comparison || !IsKnownScope(scope)) {
        return std::string();
    }
    std::ostringstream json;
    json << "{\"K\":" << shape.K
         << ",\"block_bytes\":" << shape.BlockBytes
         << ",\"comparison\":\"" << comparison
         << "\",\"scope\":\"" << scope << "\"}";
    return wirehair::wh2_benchmark::Sha256Hex(json.str());
}

NativePanelOrder PanelOrder(
    const CellShape& shape,
    const char* comparison,
    const char* scope,
    uint32_t replicate)
{
    const std::string digest = PanelKeySha256(shape, comparison, scope);
    if (digest.size() != 64u) {
        return static_cast<NativePanelOrder>(UINT32_MAX);
    }
    const char last = digest.back();
    const uint32_t nibble = last >= '0' && last <= '9' ?
        static_cast<uint32_t>(last - '0') :
        static_cast<uint32_t>(last - 'a' + 10);
    return ((replicate & 1u) == (nibble & 1u)) ?
        NativePanelOrder::ABBA : NativePanelOrder::BAAB;
}

bool ParseCpu(const char* text, int& cpu_out)
{
    if (!text || !*text) return false;
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value < 0 ||
        value > INT_MAX)
    {
        return false;
    }
    cpu_out = static_cast<int>(value);
    return true;
}

bool SamePrecodeParams(
    const wirehair_v2::PrecodeParams& left,
    const wirehair_v2::PrecodeParams& right)
{
    return left.BlockCount == right.BlockCount &&
        left.Staircase == right.Staircase &&
        left.DenseRows == right.DenseRows &&
        left.HeavyRows == right.HeavyRows &&
        left.SourceHits == right.SourceHits &&
        left.DenseIdentityCorner == right.DenseIdentityCorner &&
        left.DenseAnchors == right.DenseAnchors &&
        left.HeavyFamily == right.HeavyFamily &&
        left.Seed == right.Seed;
}

bool SameResolved(
    const ResolvedNativeWh2Configuration& left,
    const ResolvedNativeWh2Configuration& right)
{
    return SamePrecodeParams(left.Params, right.Params) &&
        left.PacketConfig.PeelSeed == right.PacketConfig.PeelSeed &&
        left.PacketConfig.MixCount == right.PacketConfig.MixCount &&
        left.PrecodeAttempt == right.PrecodeAttempt &&
        left.PacketAttempt == right.PacketAttempt;
}

NativeArmSpec ArmSpec(NativeSystematicEmission emission)
{
    NativeArmSpec spec = wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    spec.SystematicEmission = emission;
    return spec;
}

std::string DescriptorJson(
    const char* scope,
    NativeSystematicEmission emission)
{
    if (!IsKnownScope(scope)) return std::string();
    const char* const name =
        wirehair_wh2_bench::NativeSystematicEmissionName(emission);
    if (!name) return std::string();
    std::string json;
    json.reserve(240u);
    const char* const acquisition =
        std::strcmp(scope, "encoder-systematic") == 0 ?
            kSourceAcquisition : kRepairSourceAcquisition;
    json += "{\"arm\":\"wirehair2_head\",\"codec\":\"wirehair2_certified";
    json += "\",\"equation_transform\":\"none\",\"schema\":\"";
    json += kDescriptorSchema;
    json += "\",\"source_acquisition\":\"";
    json += acquisition;
    json += "\",\"source_storage\":\"";
    json += kSourceStorage;
    json += "\",\"systematic_emission\":\"";
    json += name;
    json += "\"}";
    return json;
}

std::string ConfigurationJson(
    const CellShape& shape,
    const ResolvedNativeWh2Configuration& resolved)
{
    const wirehair_v2::PrecodeParams& params = resolved.Params;
    std::ostringstream json;
    json << "{\"K\":" << shape.K
         << ",\"block_bytes\":" << shape.BlockBytes
         << ",\"dense_anchor_layout\":"
         << static_cast<uint32_t>(params.DenseAnchors)
         << ",\"dense_identity_corner\":"
         << (params.DenseIdentityCorner ? "true" : "false")
         << ",\"dense_rows\":" << params.DenseRows
         << ",\"heavy_family\":"
         << static_cast<uint32_t>(params.HeavyFamily)
         << ",\"heavy_rows\":" << params.HeavyRows
         << ",\"mix_count\":" << resolved.PacketConfig.MixCount
         << ",\"packet_attempt\":" << resolved.PacketAttempt
         << ",\"packet_peel_seed\":"
         << resolved.PacketConfig.PeelSeed
         << ",\"precode_attempt\":" << resolved.PrecodeAttempt
         << ",\"precode_seed\":" << params.Seed
         << ",\"source_hits\":" << params.SourceHits
         << ",\"staircase\":" << params.Staircase << "}";
    return json.str();
}

bool ValidReceipt(
    const NativeEncoderStorageReceipt& receipt,
    NativeSystematicEmission expected)
{
    return receipt.SourceStorage && receipt.SourceAcquisition &&
        std::strcmp(receipt.SourceStorage, kSourceStorage) == 0 &&
        std::strcmp(receipt.SourceAcquisition, kSourceAcquisition) == 0 &&
        receipt.SystematicEmission == expected;
}

bool BuildCell(const CellShape& shape, ScreenCell& cell_out)
{
    ScreenCell cell;
    cell.Shape = shape;
    if (!wirehair_wh2_bench::MakeDeterministicSource(
            shape.K,
            shape.BlockBytes,
            UINT64_C(0x26c1f7a98540db3e) ^
                (static_cast<uint64_t>(shape.K) << 17) ^ shape.BlockBytes,
            cell.Source))
    {
        return false;
    }

    const NativeArmSpec equation =
        ArmSpec(NativeSystematicEmission::EquationEvaluation);
    const NativeArmSpec direct =
        ArmSpec(NativeSystematicEmission::RetainedSourceDirect);
    ResolvedNativeWh2Configuration equation_config;
    ResolvedNativeWh2Configuration direct_config;
    if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
            equation, shape.K, shape.BlockBytes, equation_config) ||
        !wirehair_wh2_bench::ResolveNativeWh2Configuration(
            direct, shape.K, shape.BlockBytes, direct_config) ||
        !SameResolved(equation_config, direct_config))
    {
        return false;
    }
    const std::string configuration = ConfigurationJson(shape, equation_config);
    cell.EquationConfigurationJson = configuration;
    cell.EquationConfigurationSha256 =
        wirehair::wh2_benchmark::Sha256Hex(configuration);
    cell.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Source.data(), cell.Source.size());
    if (cell.EquationConfigurationSha256.size() != 64u ||
        cell.SourceSha256.size() != 64u)
    {
        return false;
    }

    NativeArm reference;
    NativeArm direct_state;
    if (reference.Initialize(
            equation, shape.K, shape.BlockBytes, cell.Source) !=
            Wirehair_Success ||
        direct_state.Initialize(
            direct, shape.K, shape.BlockBytes, cell.Source) !=
            Wirehair_Success ||
        !reference.HasSameWh2SolvedState(direct_state))
    {
        return false;
    }
    try {
        cell.ExpectedRepairs.assign(
            static_cast<size_t>(shape.K) * shape.BlockBytes, 0u);
    }
    catch (...) {
        return false;
    }
    for (uint32_t offset = 0u; offset < shape.K; ++offset)
    {
        if (reference.EncodeInto(
                shape.K + offset,
                cell.ExpectedRepairs.data() +
                    static_cast<size_t>(offset) * shape.BlockBytes,
                shape.BlockBytes) != Wirehair_Success)
        {
            return false;
        }
    }
    cell.RepairOracleSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.ExpectedRepairs.data(), cell.ExpectedRepairs.size());
    if (cell.RepairOracleSha256.size() != 64u) return false;
    std::vector<uint8_t> equation_packet(shape.BlockBytes, 0u);
    std::vector<uint8_t> direct_packet(shape.BlockBytes, 0u);
    const uint32_t high_ids[] = {
        shape.K + 7u, UINT32_C(0x80000000), UINT32_MAX
    };
    for (uint32_t id : high_ids)
    {
        if (reference.EncodeInto(
                id, equation_packet.data(), shape.BlockBytes) !=
                Wirehair_Success ||
            direct_state.EncodeInto(
                id, direct_packet.data(), shape.BlockBytes) !=
                Wirehair_Success ||
            equation_packet != direct_packet)
        {
            return false;
        }
    }
    cell.HighRepairControlsVerified = true;
    cell_out = std::move(cell);
    return true;
}

class EncoderInvocation : public NativePanelInvocation
{
public:
    EncoderInvocation(
        const std::shared_ptr<const ScreenCell>& cell,
        NativeSystematicEmission emission,
        const std::string& identity)
        : Cell(cell)
        , Emission(emission)
        , IdentityValue(identity)
    {
        const NativeArmSpec spec = ArmSpec(emission);
        Ready = Cell && Fixture.Initialize(
            spec,
            Cell->Shape.K,
            Cell->Shape.BlockBytes,
            Cell->Source) == Wirehair_Success;
        NativeEncoderStorageReceipt receipt;
        Ready = Ready && Fixture.GetStorageReceipt(receipt) &&
            ValidReceipt(receipt, emission);
    }

    std::string Identity() const override { return IdentityValue; }

    NativePanelInvocationResult Invoke() override
    {
        if (!Ready) {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                Wirehair_Error,
                false,
                0u,
                0u);
        }
        const wirehair_wh2_bench::TimedArmResult timed = Fixture.Run();
        const uint64_t expected_direct =
            Emission == NativeSystematicEmission::RetainedSourceDirect ?
                Cell->Shape.K : 0u;
        if (timed.Result != Wirehair_Success || !timed.BytesVerified ||
            timed.ElapsedNanoseconds == 0u ||
            timed.DirectSystematicPackets != expected_direct)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::PreflightFailure,
                timed.Result,
                true,
                0u,
                0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            timed.Result,
            true,
            static_cast<uint32_t>(timed.DirectSystematicPackets),
            timed.ElapsedNanoseconds);
    }

private:
    std::shared_ptr<const ScreenCell> Cell;
    NativeSystematicEmission Emission;
    std::string IdentityValue;
    NativeEncoderFixture Fixture;
    bool Ready = false;
};

class RepairInvocation : public NativePanelInvocation
{
public:
    RepairInvocation(
        const std::shared_ptr<const ScreenCell>& cell,
        NativeSystematicEmission emission,
        const std::string& identity)
        : Cell(cell)
        , IdentityValue(identity)
    {
        if (!Cell) return;
        try {
            Output.assign(Cell->ExpectedRepairs.size(), 0u);
        }
        catch (...) {
            return;
        }
        Ready = Arm.Initialize(
            ArmSpec(emission),
            Cell->Shape.K,
            Cell->Shape.BlockBytes,
            Cell->Source) == Wirehair_Success &&
            Arm.SystematicEmission() == emission &&
            Arm.UsesRetainedSourceFor(0u) ==
                (emission ==
                    NativeSystematicEmission::RetainedSourceDirect) &&
            Arm.UsesRetainedSourceFor(Cell->Shape.K - 1u) ==
                (emission ==
                    NativeSystematicEmission::RetainedSourceDirect) &&
            !Arm.UsesRetainedSourceFor(Cell->Shape.K) &&
            !Arm.UsesRetainedSourceFor(UINT32_MAX);
    }

    std::string Identity() const override { return IdentityValue; }

    NativePanelInvocationResult Invoke() override
    {
        if (!Ready) {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                Wirehair_Error,
                false,
                0u,
                0u);
        }
        WirehairResult result = Wirehair_Success;
        const std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
        for (uint32_t offset = 0u; offset < Cell->Shape.K; ++offset)
        {
            result = Arm.EncodeInto(
                Cell->Shape.K + offset,
                Output.data() +
                    static_cast<size_t>(offset) * Cell->Shape.BlockBytes,
                Cell->Shape.BlockBytes);
            if (result != Wirehair_Success) break;
        }
        const std::chrono::steady_clock::time_point finish =
            std::chrono::steady_clock::now();
        const std::chrono::nanoseconds elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                finish - start);
        const uint64_t elapsed_ns = elapsed.count() > 0 ?
            static_cast<uint64_t>(elapsed.count()) : 0u;
        const bool verified = result == Wirehair_Success &&
            Output == Cell->ExpectedRepairs;
        if (!verified || elapsed_ns == 0u) {
            return NativePanelInvocationResult(
                NativePanelDisposition::PreflightFailure,
                result,
                true,
                0u,
                0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            result,
            true,
            0u,
            elapsed_ns);
    }

private:
    std::shared_ptr<const ScreenCell> Cell;
    std::string IdentityValue;
    NativeArm Arm;
    std::vector<uint8_t> Output;
    bool Ready = false;
};

NativePanelArm MakePanelArm(
    const std::shared_ptr<const ScreenCell>& cell,
    const char* scope,
    NativeSystematicEmission emission,
    const std::string& descriptor_sha256)
{
    if (!IsKnownScope(scope)) return NativePanelArm();
    std::ostringstream identity;
    identity << scope << ':' << cell->Shape.K << ':' << cell->Shape.BlockBytes
             << ':' << descriptor_sha256;
    const std::string identity_value = identity.str();
    if (std::strcmp(scope, "encoder-systematic") == 0)
    {
        return NativePanelArm(
            identity_value,
            [cell, emission, identity_value]() {
                return std::unique_ptr<NativePanelInvocation>(
                    new EncoderInvocation(cell, emission, identity_value));
            });
    }
    return NativePanelArm(
        identity_value,
        [cell, emission, identity_value]() {
            return std::unique_ptr<NativePanelInvocation>(
                new RepairInvocation(cell, emission, identity_value));
        });
}

void EmitConfig(
    int target_cpu,
    const std::vector<std::shared_ptr<const ScreenCell> >& cells,
    const std::string& encoder_equation_descriptor,
    const std::string& encoder_equation_sha256,
    const std::string& encoder_direct_descriptor,
    const std::string& encoder_direct_sha256,
    const std::string& repair_equation_descriptor,
    const std::string& repair_equation_sha256,
    const std::string& repair_direct_descriptor,
    const std::string& repair_direct_sha256)
{
    std::cout << "{\"cells\":[";
    for (size_t i = 0u; i < cells.size(); ++i)
    {
        if (i != 0u) std::cout << ',';
        std::cout << "{\"K\":" << cells[i]->Shape.K
                  << ",\"block_bytes\":" << cells[i]->Shape.BlockBytes
                  << ",\"construction_equivalent\":true"
                  << ",\"equation_configuration\":"
                  << cells[i]->EquationConfigurationJson
                  << ",\"equation_configuration_sha256\":\""
                  << cells[i]->EquationConfigurationSha256
                  << "\",\"high_repair_controls_verified\":"
                  << (cells[i]->HighRepairControlsVerified ? "true" : "false")
                  << ",\"repair_oracle_sha256\":\""
                  << cells[i]->RepairOracleSha256
                  << "\",\"source_sha256\":\""
                  << cells[i]->SourceSha256 << "\"}";
    }
    std::cout << "],\"encoder_equation_descriptor\":"
              << encoder_equation_descriptor
              << ",\"encoder_equation_descriptor_sha256\":\""
              << encoder_equation_sha256
              << "\",\"encoder_retained_descriptor\":"
              << encoder_direct_descriptor
              << ",\"encoder_retained_descriptor_sha256\":\""
              << encoder_direct_sha256
              << "\",\"internal_deadline_seconds\":"
              << kInternalDeadlineSeconds
              << ",\"invocation_budget\":" << kInvocationBudget
              << ",\"panel_replicates\":" << kPanelReplicates
              << ",\"repair_equation_descriptor\":"
              << repair_equation_descriptor
              << ",\"repair_equation_descriptor_sha256\":\""
              << repair_equation_sha256
              << "\",\"repair_retained_descriptor\":"
              << repair_direct_descriptor
              << ",\"repair_retained_descriptor_sha256\":\""
              << repair_direct_sha256
              << "\",\"schema\":\"" << kConfigSchema
              << "\",\"source_git_commit\":\""
              << WIREHAIR_WH2_SOURCE_GIT_COMMIT
              << "\",\"target_cpu\":" << target_cpu << "}\n";
}

void EmitPanel(
    const ScreenCell& cell,
    const char* comparison,
    uint32_t invocations_per_slot,
    const std::string& left_sha256,
    NativePanelOrder order,
    const std::string& panel_key_sha256,
    uint32_t replicate,
    const std::string& right_sha256,
    const char* scope,
    const NativePanelResult& panel)
{
    std::cout << "{\"K\":" << cell.Shape.K
              << ",\"block_bytes\":" << cell.Shape.BlockBytes
              << ",\"comparison\":\"" << comparison
              << "\",\"invocations_per_slot\":" << invocations_per_slot
              << ",\"left_descriptor_sha256\":\"" << left_sha256
              << "\",\"left_direct_systematic_packets\":"
              << panel.LeftPreflight.DecodedExtra
              << ",\"left_outcome_code\":"
              << panel.LeftPreflight.OutcomeCode
              << ",\"order\":\"" << OrderName(order)
              << "\",\"panel_key_sha256\":\"" << panel_key_sha256
              << "\",\"replicate\":" << replicate
              << ",\"right_descriptor_sha256\":\"" << right_sha256
              << "\",\"right_direct_systematic_packets\":"
              << panel.RightPreflight.DecodedExtra
              << ",\"right_outcome_code\":"
              << panel.RightPreflight.OutcomeCode
              << ",\"schema\":\"" << kPanelSchema
              << "\",\"scope\":\"" << scope << "\",\"slots\":[";
    for (size_t i = 0u; i < panel.Slots.size(); ++i)
    {
        if (i != 0u) std::cout << ',';
        std::cout << "{\"elapsed_ns\":"
                  << panel.Slots[i].ElapsedNanoseconds
                  << ",\"side\":\"" << SideName(panel.Slots[i].Side)
                  << "\"}";
    }
    std::cout << "],\"target_cpu\":" << panel.TargetCpu << "}\n";
}

bool SelfTest()
{
    ScreenCell cell;
    if (!BuildCell(kCellShapes[0], cell)) return false;
    for (size_t i = 0u;
         i < sizeof(kCellShapes) / sizeof(kCellShapes[0]);
         ++i)
    {
        ScreenCell receipt_cell;
        if (!BuildCell(kCellShapes[i], receipt_cell) ||
            receipt_cell.EquationConfigurationSha256 !=
                kExpectedConfigurationSha256[i] ||
            receipt_cell.SourceSha256 != kExpectedSourceSha256[i] ||
            receipt_cell.RepairOracleSha256 !=
                kExpectedRepairOracleSha256[i])
        {
            return false;
        }
    }
    const std::string equation =
        DescriptorJson(
            "encoder-systematic",
            NativeSystematicEmission::EquationEvaluation);
    const std::string direct =
        DescriptorJson(
            "encoder-systematic",
            NativeSystematicEmission::RetainedSourceDirect);
    if (equation.empty() || direct.empty() || equation == direct ||
        !DescriptorJson(
            "unknown", NativeSystematicEmission::EquationEvaluation).empty() ||
        equation.find("\"systematic_emission\":\"equation-eval-v1\"") ==
            std::string::npos ||
        direct.find(
            "\"systematic_emission\":\"retained-source-direct-v1\"") ==
            std::string::npos)
    {
        return false;
    }
    NativeArm equation_arm;
    NativeArm direct_arm;
    if (equation_arm.Initialize(
            ArmSpec(NativeSystematicEmission::EquationEvaluation),
            cell.Shape.K,
            cell.Shape.BlockBytes,
            cell.Source) != Wirehair_Success ||
        direct_arm.Initialize(
            ArmSpec(NativeSystematicEmission::RetainedSourceDirect),
            cell.Shape.K,
            cell.Shape.BlockBytes,
            cell.Source) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> left(cell.Shape.BlockBytes, 0u);
    std::vector<uint8_t> right(cell.Shape.BlockBytes, 0u);
    const uint32_t ids[] = {
        0u, 1u, cell.Shape.K - 1u, cell.Shape.K,
        cell.Shape.K + 7u, UINT32_MAX
    };
    for (uint32_t id : ids)
    {
        if (equation_arm.EncodeInto(
                id, left.data(), cell.Shape.BlockBytes) != Wirehair_Success ||
            direct_arm.EncodeInto(
                id, right.data(), cell.Shape.BlockBytes) != Wirehair_Success ||
            left != right)
        {
            return false;
        }
    }
    return direct_arm.UsesRetainedSourceFor(0u) &&
        direct_arm.UsesRetainedSourceFor(cell.Shape.K - 1u) &&
        !direct_arm.UsesRetainedSourceFor(cell.Shape.K) &&
        !direct_arm.UsesRetainedSourceFor(UINT32_MAX) &&
        !equation_arm.UsesRetainedSourceFor(0u);
}

bool RunScreen(int target_cpu)
{
    if (!wirehair_wh2_bench::NativePanelPlatformSupported()) {
        std::cerr << "systematic emission screen requires Linux affinity\n";
        return false;
    }
    const std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::now() +
        std::chrono::seconds(kInternalDeadlineSeconds);

    std::vector<std::shared_ptr<const ScreenCell> > cells;
    for (const CellShape& shape : kCellShapes)
    {
        std::shared_ptr<ScreenCell> cell(new ScreenCell);
        if (!BuildCell(shape, *cell)) {
            std::cerr << "cannot build exact systematic emission cell\n";
            return false;
        }
        cells.push_back(cell);
    }

    const std::string encoder_equation_descriptor = DescriptorJson(
        "encoder-systematic", NativeSystematicEmission::EquationEvaluation);
    const std::string encoder_direct_descriptor = DescriptorJson(
        "encoder-systematic", NativeSystematicEmission::RetainedSourceDirect);
    const std::string repair_equation_descriptor = DescriptorJson(
        "repair-evaluate", NativeSystematicEmission::EquationEvaluation);
    const std::string repair_direct_descriptor = DescriptorJson(
        "repair-evaluate", NativeSystematicEmission::RetainedSourceDirect);
    const std::string encoder_equation_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(encoder_equation_descriptor);
    const std::string encoder_direct_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(encoder_direct_descriptor);
    const std::string repair_equation_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(repair_equation_descriptor);
    const std::string repair_direct_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(repair_direct_descriptor);
    if (encoder_equation_descriptor.empty() ||
        encoder_direct_descriptor.empty() ||
        repair_equation_descriptor.empty() ||
        repair_direct_descriptor.empty() ||
        encoder_equation_sha256.size() != 64u ||
        encoder_direct_sha256.size() != 64u ||
        repair_equation_sha256.size() != 64u ||
        repair_direct_sha256.size() != 64u ||
        encoder_equation_sha256 == encoder_direct_sha256 ||
        repair_equation_sha256 == repair_direct_sha256)
    {
        std::cerr << "invalid systematic emission descriptor seal\n";
        return false;
    }
    EmitConfig(
        target_cpu,
        cells,
        encoder_equation_descriptor,
        encoder_equation_sha256,
        encoder_direct_descriptor,
        encoder_direct_sha256,
        repair_equation_descriptor,
        repair_equation_sha256,
        repair_direct_descriptor,
        repair_direct_sha256);

    static const char* const scopes[] = {
        "encoder-systematic", "repair-evaluate"
    };
    struct Comparison
    {
        const char* Name;
        NativeSystematicEmission Left;
        NativeSystematicEmission Right;
    };
    static const Comparison comparisons[] = {
        { "baseline-aa",
          NativeSystematicEmission::EquationEvaluation,
          NativeSystematicEmission::EquationEvaluation },
        { "candidate-aa",
          NativeSystematicEmission::RetainedSourceDirect,
          NativeSystematicEmission::RetainedSourceDirect },
        { "candidate-over-baseline",
          NativeSystematicEmission::RetainedSourceDirect,
          NativeSystematicEmission::EquationEvaluation }
    };

    uint32_t panel_count = 0u;
    for (uint32_t replicate = 0u; replicate < kPanelReplicates; ++replicate)
    {
        for (const std::shared_ptr<const ScreenCell>& cell : cells)
        {
            const uint32_t invocations = InvocationsPerSlot(
                cell->Shape.K, replicate);
            for (const char* scope : scopes)
            {
                for (const Comparison& comparison : comparisons)
                {
                    const std::string panel_key_sha256 = PanelKeySha256(
                        cell->Shape, comparison.Name, scope);
                    const NativePanelOrder order = PanelOrder(
                        cell->Shape, comparison.Name, scope, replicate);
                    if (panel_key_sha256.size() != 64u ||
                        (order != NativePanelOrder::ABBA &&
                         order != NativePanelOrder::BAAB))
                    {
                        std::cerr << "invalid systematic emission panel key\n";
                        return false;
                    }
                    if (std::chrono::steady_clock::now() >= deadline) {
                        std::cerr << "systematic emission internal deadline\n";
                        return false;
                    }
                    const bool encoder_scope =
                        std::strcmp(scope, "encoder-systematic") == 0;
                    const std::string& scope_equation_sha = encoder_scope ?
                        encoder_equation_sha256 : repair_equation_sha256;
                    const std::string& scope_direct_sha = encoder_scope ?
                        encoder_direct_sha256 : repair_direct_sha256;
                    const std::string& left_sha = comparison.Left ==
                            NativeSystematicEmission::RetainedSourceDirect ?
                            scope_direct_sha : scope_equation_sha;
                    const std::string& right_sha = comparison.Right ==
                            NativeSystematicEmission::RetainedSourceDirect ?
                            scope_direct_sha : scope_equation_sha;
                    const NativePanelArm left = MakePanelArm(
                        cell, scope, comparison.Left, left_sha);
                    const NativePanelArm right = MakePanelArm(
                        cell, scope, comparison.Right, right_sha);
                    const NativePanelResult panel =
                        wirehair_wh2_bench::ExecuteNativeTimingPanel(
                            target_cpu,
                            order,
                            invocations,
                            left,
                            right);
                    if (panel.Status != NativePanelStatus::Complete ||
                        !panel.PanelComparable || panel.IsFatal())
                    {
                        std::cerr << "systematic emission panel failed: "
                                  << wirehair_wh2_bench::NativePanelStatusName(
                                         panel.Status)
                                  << ' ' << panel.Diagnostic << '\n';
                        return false;
                    }
                    const uint32_t expected_left_direct =
                        encoder_scope && comparison.Left ==
                            NativeSystematicEmission::RetainedSourceDirect ?
                            cell->Shape.K : 0u;
                    const uint32_t expected_right_direct =
                        encoder_scope && comparison.Right ==
                            NativeSystematicEmission::RetainedSourceDirect ?
                            cell->Shape.K : 0u;
                    if (!panel.HasLeftPreflight ||
                        !panel.HasRightPreflight ||
                        panel.LeftPreflight.Disposition !=
                            NativePanelDisposition::Success ||
                        panel.RightPreflight.Disposition !=
                            NativePanelDisposition::Success ||
                        panel.LeftPreflight.OutcomeCode != Wirehair_Success ||
                        panel.RightPreflight.OutcomeCode != Wirehair_Success ||
                        !panel.LeftPreflight.HasDecodedExtra ||
                        !panel.RightPreflight.HasDecodedExtra ||
                        panel.LeftPreflight.DecodedExtra !=
                            expected_left_direct ||
                        panel.RightPreflight.DecodedExtra !=
                            expected_right_direct)
                    {
                        std::cerr <<
                            "systematic emission preflight receipt drifted\n";
                        return false;
                    }
                    for (const wirehair_wh2_bench::NativePanelSlot& slot :
                         panel.Slots)
                    {
                        if (!slot.HasElapsedNanoseconds ||
                            slot.ElapsedNanoseconds == 0u)
                        {
                            std::cerr << "systematic emission panel has null time\n";
                            return false;
                        }
                    }
                    EmitPanel(
                        *cell,
                        comparison.Name,
                        invocations,
                        left_sha,
                        order,
                        panel_key_sha256,
                        replicate,
                        right_sha,
                        scope,
                        panel);
                    ++panel_count;
                }
            }
        }
    }
    if (std::chrono::steady_clock::now() >= deadline) {
        std::cerr << "systematic emission internal deadline\n";
        return false;
    }
    std::cout << "{\"panel_count\":" << panel_count
              << ",\"schema\":\"" << kTerminalSchema
              << "\",\"status\":\"complete\"}\n";
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0)
    {
        if (!SelfTest()) {
            std::cerr << "WH2 systematic emission screen selftest failed\n";
            return 1;
        }
        std::cout << "WH2 systematic emission screen selftest passed\n";
        return 0;
    }
    int target_cpu = -1;
    if (argc != 4 || std::strcmp(argv[1], "--run") != 0 ||
        std::strcmp(argv[2], "--cpu") != 0 ||
        !ParseCpu(argv[3], target_cpu))
    {
        std::cerr << "usage: " << argv[0]
                  << " --selftest | --run --cpu <nonnegative-int>\n";
        return 2;
    }
    return RunScreen(target_cpu) ? 0 : 1;
}
