#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"
#include "Wh2NativePanel.h"

#include <algorithm>
#include <array>
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
#error "broad systematic emission screen requires an exact source commit receipt"
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
    "wirehair.wh2.broad-systematic-emission-screen.config.v1";
static const char kPanelSchema[] =
    "wirehair.wh2.broad-systematic-emission-screen.panel.v1";
static const char kTerminalSchema[] =
    "wirehair.wh2.broad-systematic-emission-screen.terminal.v1";
static const char kDescriptorSchema[] =
    "wirehair.wh2.broad-systematic-emission-arm.v1";
static const char kSourceStorage[] = "native-arm-owned-kxb-v1";
static const uint32_t kPanelReplicates = 12u;
static const uint32_t kInvocationBudget = 65536u;
static const uint32_t kMinimumInvocations = 24u;
static const uint32_t kInternalDeadlineSeconds = 115u;

static const char kFreshSystematicScope[] = "fresh-systematic";
static const char kFreshRepairScope[] = "fresh-repair";
static const char kSteadyRepairScope[] = "steady-repair";
static const char kFreshSystematicMetric[] = "fresh-systematic-total";
static const char kFreshRepairMetric[] = "fresh-repair-total";
static const char kFreshRepairNestedMetric[] =
    "fresh-repair-init-first";
static const char kSteadyRepairMetric[] = "steady-repair-total";

struct CellShape
{
    uint32_t K;
    uint32_t BlockBytes;
};

static const CellShape kCellShapes[] = {
    { 8u, 64u }, { 8u, 1280u },
    { 128u, 64u }, { 128u, 1280u },
    { 512u, 64u }, { 512u, 1280u },
    { 5000u, 64u }, { 5000u, 1280u },
    { 64000u, 64u }, { 64000u, 1280u }
};

// Frozen from the untimed --selftest receipt.  Every array is indexed by
// kCellShapes and intentionally duplicated in the Python adjudicator so
// neither side can silently expand the roster.
static const char* const kExpectedConfigurationSha256[] = {
    "f696d3f825d9d8b56d4f667d423f1fca0d76db5cea156b2bae18ed27e5e3369e",
    "39b3f23adeabdb6de6f63144e74008b526d653ff06fd01748fe0b1f07a7eed24",
    "1cecbe75d7c35c5f4413ba602e8919f6bdc90a6c8b1dfdc467894ba5c9270572",
    "9e1ee0fa094d14e1ced2226af42fa762405df1e4f59a4f212bdd53a522617951",
    "5f201020c4e541b37c31e2cecb20088e1c77794ebb266a5d8f9f42e412bfac71",
    "f3772040f91e7a58eb7fab8b4be44a726afed8b606ab0d21af3b007faeea3af8",
    "96b443b6fe64b4a6d29ca4f307f3372ed13951356e3de9a568952ce7f9d5829b",
    "ad0769817c845f6ed6eda1ada8a5c7d048e7ee31f9303302485831983eea2696",
    "112530ebf8646e224936bfb32552591dbecb0cb133df02e0f841b61eccc267bc",
    "6ef35f5bf51ad0688504b1590c179c476c9d6cd31a73424ea8be364f692d0455"
};
static const char* const kExpectedSourceSha256[] = {
    "da7d0d04aaeb241cc351650a839fd974f99be26766b7697de0af6b23a72c1276",
    "b609fa6e4de08903bd33716155bf3936f53fe9cb14fe4ad9c44a1cbb55604adb",
    "24389e04d91e44f451a9d8f801be054aa6c78485819eb25676838ae2e9765f4d",
    "d8216125397beca9146bb3dddc87a904adca60fcee659506b8546f29a2a1ebb9",
    "5bf428d7f837c9e3983f7116d6dd1516e9abfdc8a0fde2ede78b837653195a26",
    "4fcecbc851e2ab882f98c9c831eae9045cfc481d69232081964f76c704f7a1ac",
    "df192ef37e7f2a187108d7659714570d070b2dcbae74e11f437bf53b64a0c840",
    "a61f8833c1deff31737a411237b47cb2494469dba9ffbe85e3dab7d102beb5f9",
    "a37eee40a067926c7979a159ddaea717cc96a21090d67b712cc8db2d7df16f18",
    "8a0df20ca26df770afc859c858877acf95c435926a7a57d8e238f3ecf076c80d"
};
static const char* const kExpectedRepairOracleSha256[] = {
    "75112775332232f5c4f2e70ef9dcfbabd41afc1aa89217374d2fbf025c42ddbf",
    "893f50a4b50f6136d515ec44229a0808adf006dd41f25db406dd92c0233787e3",
    "04bc95892a4eee3afea70fc33030735a76c5a69af2738dd12c626bb6fe721b47",
    "5353ff791ae97b657d7d0fde48e6ea537bca80097bb6de65253ecb66cf76d6a4",
    "17f5b2ca5321774bb9740b0cfcbb4e06053a1af88800928b067f389d037c0dbb",
    "e75a45a5d8af915c4f24dd61e2348c16474d7452d8e238555b6b83bd287000b8",
    "da5a781cabfcacac7b765c2cd9d87f00e840a33df943619b81a37d6d314a6ee8",
    "cc0ab889b23c9f2a440e9edf46a1f85e4c09e1a1a9620952cdd0791058f50acd",
    "38f7f201fac458b70d3d73e6a8d2c8703119c7bf91f3b8bc949352564c09d18e",
    "a6083d31b6ba6b2dc0e4e5214317b2a16e8b3ce0a7b0643fedd8fc7a5807e4cd"
};
static const char* const kExpectedFirstRepairSha256[] = {
    "422a576a15228c4fe99ae3d72b220c4b9c3eff12981778733ffcda41ff86c4ab",
    "5f5246fb387706ca99faa6e958c916f8882cf3b70339fa249984e5690012b94d",
    "aee24d5a978e4618e4dd9b86f9f8a82273e1200a9ad817fe89c954565c9b2438",
    "e7f4ec656551ae770f10245da4c187b7a16c61d70ffcfc8e7877898c47e404e3",
    "4b78dcec571a957965b27e22fc213b87afb3e5c52059a9fa128fd7946bff67e5",
    "f063214b55a1be96d89bf6f49d7c212d99c30ae6f835d1bd8c2b50cde9aca4b7",
    "061d006b12581421ae9f1cf0f50de6284e138c716b5b01cd29862266f48c05e0",
    "951d348686c4c5b249e009e7e22584757ebfe123b42d0a2260de6bb0453d97b3",
    "f667a2be87515e63e70f9597314fc16c98e0c33e8c0f3870a6a73c62a9134219",
    "fb85f157c691fda71203e13a2ce775a12ef240aadcfd9156efbbe0b115fb4db4"
};
static const char* const kExpectedHighIdOracleSha256[] = {
    "cea33c827c0102fffd20f8d6ddefb39474f8a63ce70a6fae17d5f0eba51eb608",
    "cbd3ea46f46890873ba75bf0cef18b0981ac3f32a526550eff2d311b0ae2ba42",
    "6ab0519fc415a5c6e3832e6e30df31ca7d23b4be3340acd4611ecdaf87aca051",
    "ccb3d70885fca5b305a733c90d9c2b086168b1d689323eb7f10929dc63ce5385",
    "0bd4ee7de5e9d967adf465d0b1795eec26f0dc10a5764fe5f2b6489e3b2c3b5a",
    "93f9a3ebf6c1020620c7384e1cc1db78a81f9ec9353e1240f609489e731638a2",
    "79b7f9d4c4cf2232e750da5661d8cf43d4b38eb9fdb645030d832e77e72dd2ef",
    "6a89ce732b75af1e513c404f6ede3808d020f9fb3fcbabea8b838a48589fb42e",
    "3d2c21f2ddfb32bb088193cf9daac9426ba26bb2d0afbd66f3e41d4c8ffa4e8b",
    "ef170b0bcf6ba809a25a47369018782eac71c2c89d33052268c574d397fd9662"
};
static const char* const kExpectedSolvedStateEquivalenceReceiptSha256[] = {
    "c765f5bfe152e6eb9376bbfcc82eda06164077650e28aa12c92170d57dc3f534",
    "123f643f8c36098aaa2b2a47c4408bf120f0d31a5a66b12bc85defc97a57859b",
    "398551197d004d987ed033e3f2a30e7289de0e6db990128cfdae2eff1519c28a",
    "3e4b62fb4a5b4dfb35045eb0c8c1130c9321fe2bb6fde029ff7da498a02cceb8",
    "8f3dc098fc821babc939505b0d7961dac12f799be6f24604965904055ce92983",
    "e1661d99ce94c5c007fbad170e5c63424c258f647e5318d7ab153a84e94f5165",
    "fce0ca96ed238b58804f7126625150b18cdc9f5c2b6bf1b252ff856b2177c0d2",
    "942cf4574936e99f0af206d3eb0be6e23aefc4995528ff749d40e19c3343635a",
    "e5d820cd7a0fc7ba09dd4d7aa79332f40e5be077dcbf978ae8eae2ec958fc76c",
    "220c9abba2e111e616f9e729eff286aa2d8bdb2a9a8b4c1dcfdc83c321b68b13"
};

struct ScreenCell
{
    CellShape Shape = {};
    std::vector<uint8_t> Source;
    std::vector<uint8_t> ExpectedRepairs;
    std::string EquationConfigurationJson;
    std::string EquationConfigurationSha256;
    std::string FirstRepairSha256;
    std::string HighIdOracleSha256;
    std::string RepairOracleSha256;
    std::string SolvedStateEquivalenceReceiptSha256;
    std::string SourceSha256;
    bool HighRepairControlsVerified = false;
};

struct NestedTimingCollector
{
    std::vector<uint64_t> Elapsed;
};

struct Comparison
{
    const char* Name;
    NativeSystematicEmission Left;
    NativeSystematicEmission Right;
};

static const Comparison kComparisons[] = {
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
        (std::strcmp(scope, kFreshSystematicScope) == 0 ||
         std::strcmp(scope, kFreshRepairScope) == 0 ||
         std::strcmp(scope, kSteadyRepairScope) == 0);
}

const char* PrimaryMetric(const char* scope)
{
    if (!scope) return nullptr;
    if (std::strcmp(scope, kFreshSystematicScope) == 0) {
        return kFreshSystematicMetric;
    }
    if (std::strcmp(scope, kFreshRepairScope) == 0) {
        return kFreshRepairMetric;
    }
    if (std::strcmp(scope, kSteadyRepairScope) == 0) {
        return kSteadyRepairMetric;
    }
    return nullptr;
}

const char* NestedMetric(const char* scope)
{
    return scope && std::strcmp(scope, kFreshRepairScope) == 0 ?
        kFreshRepairNestedMetric : "none";
}

uint32_t TotalInvocations(uint32_t K)
{
    const uint32_t budget = (kInvocationBudget + K - 1u) / K;
    return std::max(kMinimumInvocations, budget);
}

uint32_t InvocationsPerSlot(uint32_t K, uint32_t replicate)
{
    const uint32_t total = TotalInvocations(K);
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

const char* EmissionMode(NativeSystematicEmission emission)
{
    return emission == NativeSystematicEmission::EquationEvaluation ?
        "equation" : "retained";
}

std::string DescriptorJson(
    const char* metric,
    NativeSystematicEmission emission)
{
    if (!metric) return std::string();
    const char* acquisition = nullptr;
    const char* timed_work = nullptr;
    if (std::strcmp(metric, kFreshSystematicMetric) == 0)
    {
        acquisition = "fixture-copy-before-clock-move-v1";
        timed_work = "fresh-eager-init-plus-systematic-ids-0-through-k-minus-1-v1";
    }
    else if (std::strcmp(metric, kFreshRepairMetric) == 0)
    {
        acquisition = "source-copy-in-fresh-init-clock-v1";
        timed_work = "fresh-eager-init-plus-repair-ids-k-through-2k-minus-1-v1";
    }
    else if (std::strcmp(metric, kFreshRepairNestedMetric) == 0)
    {
        acquisition = "source-copy-in-fresh-init-clock-v1";
        timed_work = "nested-fresh-eager-init-plus-first-repair-id-k-v1";
    }
    else if (std::strcmp(metric, kSteadyRepairMetric) == 0)
    {
        acquisition = "native-arm-copy-before-steady-clock-v1";
        timed_work = "prebuilt-repair-ids-k-through-2k-minus-1-v1";
    }
    if (!acquisition || !timed_work) return std::string();
    const char* const emission_name =
        wirehair_wh2_bench::NativeSystematicEmissionName(emission);
    if (!emission_name) return std::string();
    std::string json;
    json.reserve(360u);
    json += "{\"arm\":\"wirehair2_head\",\"codec\":\"wirehair2_certified";
    json += "\",\"equation_transform\":\"none\",\"metric\":\"";
    json += metric;
    json += "\",\"schema\":\"";
    json += kDescriptorSchema;
    json += "\",\"source_acquisition\":\"";
    json += acquisition;
    json += "\",\"source_storage\":\"";
    json += kSourceStorage;
    json += "\",\"systematic_emission\":\"";
    json += emission_name;
    json += "\",\"timed_work\":\"";
    json += timed_work;
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
         << ",\"packet_peel_seed\":" << resolved.PacketConfig.PeelSeed
         << ",\"precode_attempt\":" << resolved.PrecodeAttempt
         << ",\"precode_seed\":" << params.Seed
         << ",\"source_hits\":" << params.SourceHits
         << ",\"staircase\":" << params.Staircase << "}";
    return json.str();
}

std::string SolvedStateEquivalenceReceiptJson(const ScreenCell& cell)
{
    std::string json;
    json.reserve(520u);
    json += "{\"configuration_sha256\":\"";
    json += cell.EquationConfigurationSha256;
    json += "\",\"first_repair_sha256\":\"";
    json += cell.FirstRepairSha256;
    json += "\",\"high_id_oracle_sha256\":\"";
    json += cell.HighIdOracleSha256;
    json += "\",\"repair_oracle_sha256\":\"";
    json += cell.RepairOracleSha256;
    json += "\",\"source_sha256\":\"";
    json += cell.SourceSha256;
    json += "\",\"systematic_oracle_sha256\":\"";
    json += cell.SourceSha256;
    json += "\"}";
    return json;
}

bool BuildCell(const CellShape& shape, ScreenCell& cell_out)
{
    ScreenCell cell;
    cell.Shape = shape;
    if (!wirehair_wh2_bench::MakeDeterministicSource(
            shape.K,
            shape.BlockBytes,
            UINT64_C(0xb6402ee71c8365a9) ^
                (static_cast<uint64_t>(shape.K) << 17) ^ shape.BlockBytes,
            cell.Source))
    {
        return false;
    }

    const NativeArmSpec equation =
        ArmSpec(NativeSystematicEmission::EquationEvaluation);
    const NativeArmSpec retained =
        ArmSpec(NativeSystematicEmission::RetainedSourceDirect);
    ResolvedNativeWh2Configuration equation_config;
    ResolvedNativeWh2Configuration retained_config;
    if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
            equation, shape.K, shape.BlockBytes, equation_config) ||
        !wirehair_wh2_bench::ResolveNativeWh2Configuration(
            retained, shape.K, shape.BlockBytes, retained_config) ||
        !SameResolved(equation_config, retained_config))
    {
        return false;
    }
    cell.EquationConfigurationJson =
        ConfigurationJson(shape, equation_config);
    cell.EquationConfigurationSha256 =
        wirehair::wh2_benchmark::Sha256Hex(
            cell.EquationConfigurationJson);
    cell.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Source.data(), cell.Source.size());
    if (cell.EquationConfigurationSha256.size() != 64u ||
        cell.SourceSha256.size() != 64u)
    {
        return false;
    }

    NativeArm equation_arm;
    NativeArm retained_arm;
    if (equation_arm.Initialize(
            equation, shape.K, shape.BlockBytes, cell.Source) !=
            Wirehair_Success ||
        retained_arm.Initialize(
            retained, shape.K, shape.BlockBytes, cell.Source) !=
            Wirehair_Success ||
        !equation_arm.HasSameWh2SolvedState(retained_arm))
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
    std::vector<uint8_t> equation_packet(shape.BlockBytes, 0u);
    std::vector<uint8_t> retained_packet(shape.BlockBytes, 0u);
    for (uint32_t id = 0u; id < shape.K; ++id)
    {
        if (equation_arm.EncodeInto(
                id, equation_packet.data(), shape.BlockBytes) !=
                Wirehair_Success ||
            retained_arm.EncodeInto(
                id, retained_packet.data(), shape.BlockBytes) !=
                Wirehair_Success ||
            equation_packet != retained_packet ||
            std::memcmp(
                equation_packet.data(),
                cell.Source.data() + static_cast<size_t>(id) *
                    shape.BlockBytes,
                shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    for (uint32_t offset = 0u; offset < shape.K; ++offset)
    {
        uint8_t* const expected = cell.ExpectedRepairs.data() +
            static_cast<size_t>(offset) * shape.BlockBytes;
        if (equation_arm.EncodeInto(
                shape.K + offset, expected, shape.BlockBytes) !=
                Wirehair_Success ||
            retained_arm.EncodeInto(
                shape.K + offset,
                retained_packet.data(),
                shape.BlockBytes) != Wirehair_Success ||
            std::memcmp(expected, retained_packet.data(), shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    cell.RepairOracleSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.ExpectedRepairs.data(), cell.ExpectedRepairs.size());
    cell.FirstRepairSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.ExpectedRepairs.data(), shape.BlockBytes);

    std::vector<uint8_t> high_outputs;
    try {
        high_outputs.assign(static_cast<size_t>(3u) * shape.BlockBytes, 0u);
    }
    catch (...) {
        return false;
    }
    const uint32_t high_ids[] = {
        shape.K + 7u, UINT32_C(0x80000000), UINT32_MAX
    };
    for (size_t index = 0u; index < 3u; ++index)
    {
        uint8_t* const output = high_outputs.data() +
            index * shape.BlockBytes;
        if (equation_arm.EncodeInto(
                high_ids[index], output, shape.BlockBytes) !=
                Wirehair_Success ||
            retained_arm.EncodeInto(
                high_ids[index],
                retained_packet.data(),
                shape.BlockBytes) != Wirehair_Success ||
            std::memcmp(output, retained_packet.data(), shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    cell.HighIdOracleSha256 = wirehair::wh2_benchmark::Sha256Hex(
        high_outputs.data(), high_outputs.size());
    if (cell.RepairOracleSha256.size() != 64u ||
        cell.FirstRepairSha256.size() != 64u ||
        cell.HighIdOracleSha256.size() != 64u)
    {
        return false;
    }
    cell.HighRepairControlsVerified = true;
    cell.SolvedStateEquivalenceReceiptSha256 =
        wirehair::wh2_benchmark::Sha256Hex(
            SolvedStateEquivalenceReceiptJson(cell));
    if (cell.SolvedStateEquivalenceReceiptSha256.size() != 64u) return false;
    cell_out = std::move(cell);
    return true;
}

bool ValidReceipt(
    const NativeEncoderStorageReceipt& receipt,
    NativeSystematicEmission expected)
{
    return receipt.SourceStorage && receipt.SourceAcquisition &&
        std::strcmp(receipt.SourceStorage, kSourceStorage) == 0 &&
        std::strcmp(
            receipt.SourceAcquisition,
            "fixture-copy-before-clock-move-v1") == 0 &&
        receipt.SystematicEmission == expected;
}

class FreshSystematicInvocation : public NativePanelInvocation
{
public:
    FreshSystematicInvocation(
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
        if (!Ready) return Fatal();
        const wirehair_wh2_bench::TimedArmResult timed = Fixture.Run();
        const uint32_t expected_direct =
            Emission == NativeSystematicEmission::RetainedSourceDirect ?
                Cell->Shape.K : 0u;
        if (timed.Result != Wirehair_Success || !timed.BytesVerified ||
            timed.ElapsedNanoseconds == 0u ||
            timed.DirectSystematicPackets != expected_direct)
        {
            return Preflight(timed.Result);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            timed.Result,
            true,
            expected_direct,
            timed.ElapsedNanoseconds);
    }

private:
    static NativePanelInvocationResult Fatal()
    {
        return NativePanelInvocationResult(
            NativePanelDisposition::Fatal, Wirehair_Error, false, 0u, 0u);
    }
    static NativePanelInvocationResult Preflight(WirehairResult result)
    {
        return NativePanelInvocationResult(
            NativePanelDisposition::PreflightFailure,
            result,
            true,
            0u,
            0u);
    }

    std::shared_ptr<const ScreenCell> Cell;
    NativeSystematicEmission Emission;
    std::string IdentityValue;
    NativeEncoderFixture Fixture;
    bool Ready = false;
};

class FreshRepairInvocation : public NativePanelInvocation
{
public:
    FreshRepairInvocation(
        const std::shared_ptr<const ScreenCell>& cell,
        NativeSystematicEmission emission,
        const std::string& identity,
        const std::shared_ptr<NestedTimingCollector>& collector)
        : Cell(cell)
        , Emission(emission)
        , IdentityValue(identity)
        , Collector(collector)
    {
        if (!Cell || !Collector) return;
        try {
            Output.assign(Cell->ExpectedRepairs.size(), 0u);
        }
        catch (...) {
            return;
        }
        Ready = wirehair_init() == Wirehair_Success;
    }

    std::string Identity() const override { return IdentityValue; }

    NativePanelInvocationResult Invoke() override
    {
        if (!Ready) return NativePanelInvocationResult(
            NativePanelDisposition::Fatal,
            Wirehair_Error,
            false,
            0u,
            0u);
        NativeArm arm;
        WirehairResult result = Wirehair_Error;
        const std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
        result = arm.Initialize(
            ArmSpec(Emission),
            Cell->Shape.K,
            Cell->Shape.BlockBytes,
            Cell->Source);
        if (result == Wirehair_Success)
        {
            result = arm.EncodeInto(
                Cell->Shape.K,
                Output.data(),
                Cell->Shape.BlockBytes);
        }
        const std::chrono::steady_clock::time_point nested_finish =
            std::chrono::steady_clock::now();
        if (result == Wirehair_Success)
        {
            for (uint32_t offset = 1u; offset < Cell->Shape.K; ++offset)
            {
                result = arm.EncodeInto(
                    Cell->Shape.K + offset,
                    Output.data() +
                        static_cast<size_t>(offset) * Cell->Shape.BlockBytes,
                    Cell->Shape.BlockBytes);
                if (result != Wirehair_Success) break;
            }
        }
        const std::chrono::steady_clock::time_point finish =
            std::chrono::steady_clock::now();
        const std::chrono::nanoseconds nested_elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                nested_finish - start);
        const std::chrono::nanoseconds total_elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);
        const uint64_t nested_ns = nested_elapsed.count() > 0 ?
            static_cast<uint64_t>(nested_elapsed.count()) : 0u;
        const uint64_t total_ns = total_elapsed.count() > 0 ?
            static_cast<uint64_t>(total_elapsed.count()) : 0u;
        const bool policy_ok = result == Wirehair_Success &&
            arm.SystematicEmission() == Emission &&
            arm.UsesRetainedSourceFor(0u) ==
                (Emission == NativeSystematicEmission::RetainedSourceDirect) &&
            arm.UsesRetainedSourceFor(Cell->Shape.K - 1u) ==
                (Emission == NativeSystematicEmission::RetainedSourceDirect) &&
            !arm.UsesRetainedSourceFor(Cell->Shape.K) &&
            !arm.UsesRetainedSourceFor(UINT32_MAX);
        const bool verified = policy_ok && Output == Cell->ExpectedRepairs;
        if (!verified || nested_ns == 0u || total_ns < nested_ns)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::PreflightFailure,
                result,
                true,
                0u,
                0u);
        }
        try {
            Collector->Elapsed.push_back(nested_ns);
        }
        catch (...) {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                Wirehair_OOM,
                false,
                0u,
                0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            result,
            true,
            0u,
            total_ns);
    }

private:
    std::shared_ptr<const ScreenCell> Cell;
    NativeSystematicEmission Emission;
    std::string IdentityValue;
    std::shared_ptr<NestedTimingCollector> Collector;
    std::vector<uint8_t> Output;
    bool Ready = false;
};

class SteadyRepairInvocation : public NativePanelInvocation
{
public:
    SteadyRepairInvocation(
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
                (emission == NativeSystematicEmission::RetainedSourceDirect) &&
            Arm.UsesRetainedSourceFor(Cell->Shape.K - 1u) ==
                (emission == NativeSystematicEmission::RetainedSourceDirect) &&
            !Arm.UsesRetainedSourceFor(Cell->Shape.K) &&
            !Arm.UsesRetainedSourceFor(UINT32_MAX);
    }

    std::string Identity() const override { return IdentityValue; }

    NativePanelInvocationResult Invoke() override
    {
        if (!Ready) return NativePanelInvocationResult(
            NativePanelDisposition::Fatal,
            Wirehair_Error,
            false,
            0u,
            0u);
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
            std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);
        const uint64_t elapsed_ns = elapsed.count() > 0 ?
            static_cast<uint64_t>(elapsed.count()) : 0u;
        if (result != Wirehair_Success || elapsed_ns == 0u ||
            Output != Cell->ExpectedRepairs)
        {
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
    const std::string& descriptor_sha256,
    const std::string& nested_descriptor_sha256,
    const std::shared_ptr<NestedTimingCollector>& collector)
{
    if (!IsKnownScope(scope)) return NativePanelArm();
    std::ostringstream identity;
    identity << scope << ':' << cell->Shape.K << ':' << cell->Shape.BlockBytes
             << ':' << descriptor_sha256 << ':'
             << (nested_descriptor_sha256.empty() ?
                    "none" : nested_descriptor_sha256);
    const std::string identity_value = identity.str();
    if (std::strcmp(scope, kFreshSystematicScope) == 0)
    {
        return NativePanelArm(
            identity_value,
            [cell, emission, identity_value]() {
                return std::unique_ptr<NativePanelInvocation>(
                    new FreshSystematicInvocation(
                        cell, emission, identity_value));
            });
    }
    if (std::strcmp(scope, kFreshRepairScope) == 0)
    {
        return NativePanelArm(
            identity_value,
            [cell, emission, identity_value, collector]() {
                return std::unique_ptr<NativePanelInvocation>(
                    new FreshRepairInvocation(
                        cell, emission, identity_value, collector));
            });
    }
    return NativePanelArm(
        identity_value,
        [cell, emission, identity_value]() {
            return std::unique_ptr<NativePanelInvocation>(
                new SteadyRepairInvocation(cell, emission, identity_value));
        });
}

bool CollectNestedSlots(
    NativePanelOrder order,
    uint32_t invocations_per_slot,
    const NestedTimingCollector& left,
    const NestedTimingCollector& right,
    std::array<uint64_t, 8>& slots_out)
{
    slots_out.fill(0u);
    if (invocations_per_slot < 2u || left.Elapsed.empty() ||
        right.Elapsed.empty())
    {
        return false;
    }
    static const NativePanelSide abba[] = {
        NativePanelSide::Left, NativePanelSide::Right,
        NativePanelSide::Right, NativePanelSide::Left,
        NativePanelSide::Right, NativePanelSide::Left,
        NativePanelSide::Left, NativePanelSide::Right
    };
    static const NativePanelSide baab[] = {
        NativePanelSide::Right, NativePanelSide::Left,
        NativePanelSide::Left, NativePanelSide::Right,
        NativePanelSide::Left, NativePanelSide::Right,
        NativePanelSide::Right, NativePanelSide::Left
    };
    const NativePanelSide* const sides = order == NativePanelOrder::ABBA ?
        abba : order == NativePanelOrder::BAAB ? baab : nullptr;
    if (!sides) return false;
    size_t left_index = 1u;  // Skip the panel executor's left warmup.
    size_t right_index = 1u; // Skip the panel executor's right warmup.
    const uint32_t primary_count =
        invocations_per_slot / 2u + invocations_per_slot % 2u;
    const uint32_t secondary_count = invocations_per_slot / 2u;
    for (size_t block = 0u; block < 2u; ++block)
    {
        const size_t first_slot = block * 4u;
        const uint32_t repeats = block == 0u ?
            primary_count : secondary_count;
        for (uint32_t repeat = 0u; repeat < repeats; ++repeat)
        {
            for (size_t offset = 0u; offset < 4u; ++offset)
            {
                const size_t slot = first_slot + offset;
                const std::vector<uint64_t>& values =
                    sides[slot] == NativePanelSide::Left ?
                        left.Elapsed : right.Elapsed;
                size_t& index = sides[slot] == NativePanelSide::Left ?
                    left_index : right_index;
                if (index >= values.size() ||
                    values[index] > UINT64_MAX - slots_out[slot])
                {
                    return false;
                }
                slots_out[slot] += values[index++];
            }
        }
    }
    return left_index == left.Elapsed.size() &&
        right_index == right.Elapsed.size() &&
        std::all_of(
            slots_out.begin(), slots_out.end(),
            [](uint64_t value) { return value != 0u; });
}

struct DescriptorSeal
{
    const char* Metric;
    NativeSystematicEmission Emission;
    std::string Json;
    std::string Sha256;
};

const DescriptorSeal* FindDescriptor(
    const std::vector<DescriptorSeal>& seals,
    const char* metric,
    NativeSystematicEmission emission)
{
    if (!metric) return nullptr;
    for (const DescriptorSeal& seal : seals)
    {
        if (seal.Emission == emission &&
            std::strcmp(seal.Metric, metric) == 0)
        {
            return &seal;
        }
    }
    return nullptr;
}

bool BuildDescriptorSeals(std::vector<DescriptorSeal>& seals_out)
{
    static const char* const metrics[] = {
        kFreshSystematicMetric,
        kFreshRepairMetric,
        kFreshRepairNestedMetric,
        kSteadyRepairMetric
    };
    std::vector<DescriptorSeal> seals;
    for (const char* metric : metrics)
    {
        for (NativeSystematicEmission emission : {
                 NativeSystematicEmission::EquationEvaluation,
                 NativeSystematicEmission::RetainedSourceDirect })
        {
            DescriptorSeal seal;
            seal.Metric = metric;
            seal.Emission = emission;
            seal.Json = DescriptorJson(metric, emission);
            seal.Sha256 = wirehair::wh2_benchmark::Sha256Hex(seal.Json);
            if (seal.Json.empty() || seal.Sha256.size() != 64u) return false;
            seals.push_back(std::move(seal));
        }
    }
    seals_out.swap(seals);
    return true;
}

void EmitConfig(
    int target_cpu,
    const std::vector<std::shared_ptr<const ScreenCell> >& cells,
    const std::vector<DescriptorSeal>& descriptors)
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
                  << "\",\"first_repair_sha256\":\""
                  << cells[i]->FirstRepairSha256
                  << "\",\"high_id_oracle_sha256\":\""
                  << cells[i]->HighIdOracleSha256
                  << "\",\"high_repair_controls_verified\":"
                  << (cells[i]->HighRepairControlsVerified ? "true" : "false")
                  << ",\"repair_oracle_sha256\":\""
                  << cells[i]->RepairOracleSha256
                  << "\",\"solved_state_equivalence_receipt_sha256\":\""
                  << cells[i]->SolvedStateEquivalenceReceiptSha256
                  << "\",\"source_sha256\":\""
                  << cells[i]->SourceSha256 << "\"}";
    }
    std::cout << "],\"descriptors\":[";
    for (size_t i = 0u; i < descriptors.size(); ++i)
    {
        if (i != 0u) std::cout << ',';
        std::cout << "{\"descriptor\":" << descriptors[i].Json
                  << ",\"descriptor_sha256\":\""
                  << descriptors[i].Sha256
                  << "\",\"mode\":\""
                  << EmissionMode(descriptors[i].Emission) << "\"}";
    }
    std::cout << "],\"internal_deadline_seconds\":"
              << kInternalDeadlineSeconds
              << ",\"invocation_budget\":" << kInvocationBudget
              << ",\"minimum_invocations\":" << kMinimumInvocations
              << ",\"panel_replicates\":" << kPanelReplicates
              << ",\"schema\":\"" << kConfigSchema
              << "\",\"source_git_commit\":\""
              << WIREHAIR_WH2_SOURCE_GIT_COMMIT
              << "\",\"target_cpu\":" << target_cpu << "}\n";
}

void EmitPanel(
    const ScreenCell& cell,
    const char* comparison,
    uint32_t invocations_per_slot,
    const std::string& left_sha256,
    const std::string& left_nested_sha256,
    NativePanelOrder order,
    const std::string& panel_key_sha256,
    uint32_t replicate,
    const std::string& right_sha256,
    const std::string& right_nested_sha256,
    const char* scope,
    const char* primary_metric,
    const NativePanelResult& panel,
    const std::array<uint64_t, 8>& nested_slots)
{
    std::cout << "{\"K\":" << cell.Shape.K
              << ",\"block_bytes\":" << cell.Shape.BlockBytes
              << ",\"comparison\":\"" << comparison
              << "\",\"invocations_per_slot\":" << invocations_per_slot
              << ",\"left_descriptor_sha256\":\"" << left_sha256
              << "\",\"left_direct_systematic_packets\":"
              << panel.LeftPreflight.DecodedExtra
              << ",\"left_nested_descriptor_sha256\":";
    if (left_nested_sha256.empty()) {
        std::cout << "null";
    }
    else {
        std::cout << '"' << left_nested_sha256 << '"';
    }
    std::cout
              << ",\"left_outcome_code\":"
              << panel.LeftPreflight.OutcomeCode
              << ",\"nested_metric\":\"" << NestedMetric(scope)
              << "\",\"order\":\"" << OrderName(order)
              << "\",\"panel_key_sha256\":\"" << panel_key_sha256
              << "\",\"primary_metric\":\"" << primary_metric
              << "\",\"replicate\":" << replicate
              << ",\"right_descriptor_sha256\":\"" << right_sha256
              << "\",\"right_direct_systematic_packets\":"
              << panel.RightPreflight.DecodedExtra
              << ",\"right_nested_descriptor_sha256\":";
    if (right_nested_sha256.empty()) {
        std::cout << "null";
    }
    else {
        std::cout << '"' << right_nested_sha256 << '"';
    }
    std::cout
              << ",\"right_outcome_code\":"
              << panel.RightPreflight.OutcomeCode
              << ",\"schema\":\"" << kPanelSchema
              << "\",\"scope\":\"" << scope << "\",\"slots\":[";
    for (size_t i = 0u; i < panel.Slots.size(); ++i)
    {
        if (i != 0u) std::cout << ',';
        std::cout << "{\"elapsed_ns\":"
                  << panel.Slots[i].ElapsedNanoseconds
                  << ",\"nested_elapsed_ns\":" << nested_slots[i]
                  << ",\"side\":\"" << SideName(panel.Slots[i].Side)
                  << "\"}";
    }
    std::cout << "],\"target_cpu\":" << panel.TargetCpu << "}\n";
}

void PrintCellSeals(const ScreenCell& cell)
{
    std::cerr << cell.Shape.K << ',' << cell.Shape.BlockBytes
              << " config=" << cell.EquationConfigurationSha256
              << " source=" << cell.SourceSha256
              << " repair=" << cell.RepairOracleSha256
              << " first=" << cell.FirstRepairSha256
              << " high=" << cell.HighIdOracleSha256
              << " equivalence_receipt="
              << cell.SolvedStateEquivalenceReceiptSha256 << '\n';
}

bool CellSealsMatch(size_t index, const ScreenCell& cell)
{
    return cell.EquationConfigurationSha256 ==
            kExpectedConfigurationSha256[index] &&
        cell.SourceSha256 == kExpectedSourceSha256[index] &&
        cell.RepairOracleSha256 == kExpectedRepairOracleSha256[index] &&
        cell.FirstRepairSha256 == kExpectedFirstRepairSha256[index] &&
        cell.HighIdOracleSha256 == kExpectedHighIdOracleSha256[index] &&
        cell.SolvedStateEquivalenceReceiptSha256 ==
            kExpectedSolvedStateEquivalenceReceiptSha256[index];
}

bool TestNestedCollection()
{
    struct Case
    {
        NativePanelOrder Order;
        uint32_t Invocations;
        std::array<uint64_t, 8> Expected;
    };
    static const Case cases[] = {
        { NativePanelOrder::ABBA, 2u,
          {{ 1u, 101u, 102u, 2u, 103u, 3u, 4u, 104u }} },
        { NativePanelOrder::BAAB, 2u,
          {{ 101u, 1u, 2u, 102u, 3u, 103u, 104u, 4u }} },
        { NativePanelOrder::ABBA, 3u,
          {{ 4u, 204u, 206u, 6u, 105u, 5u, 6u, 106u }} },
        { NativePanelOrder::BAAB, 3u,
          {{ 204u, 4u, 6u, 206u, 5u, 105u, 106u, 6u }} }
    };
    for (const Case& test : cases)
    {
        NestedTimingCollector left;
        NestedTimingCollector right;
        left.Elapsed.push_back(UINT64_C(900000001));
        right.Elapsed.push_back(UINT64_C(900000101));
        const uint32_t values_per_side = 2u * test.Invocations;
        for (uint32_t i = 1u; i <= values_per_side; ++i)
        {
            left.Elapsed.push_back(i);
            right.Elapsed.push_back(100u + i);
        }
        std::array<uint64_t, 8> actual = {{ 0u }};
        if (!CollectNestedSlots(
                test.Order,
                test.Invocations,
                left,
                right,
                actual) || actual != test.Expected)
        {
            return false;
        }
        left.Elapsed.push_back(UINT64_C(777));
        if (CollectNestedSlots(
                test.Order,
                test.Invocations,
                left,
                right,
                actual))
        {
            return false;
        }
    }
    return true;
}

bool SelfTest()
{
    if (!TestNestedCollection()) return false;
    bool seals_match = true;
    for (size_t i = 0u;
         i < sizeof(kCellShapes) / sizeof(kCellShapes[0]);
         ++i)
    {
        ScreenCell cell;
        if (!BuildCell(kCellShapes[i], cell)) return false;
        if (!CellSealsMatch(i, cell)) {
            PrintCellSeals(cell);
            seals_match = false;
        }
    }
    if (!seals_match) return false;
    std::vector<DescriptorSeal> descriptors;
    if (!BuildDescriptorSeals(descriptors) || descriptors.size() != 8u) {
        return false;
    }
    for (size_t i = 0u; i < descriptors.size(); i += 2u)
    {
        if (std::strcmp(
                descriptors[i].Metric, descriptors[i + 1u].Metric) != 0 ||
            descriptors[i].Json == descriptors[i + 1u].Json ||
            descriptors[i].Sha256 == descriptors[i + 1u].Sha256)
        {
            return false;
        }
    }
    for (const CellShape& shape : kCellShapes)
    {
        uint32_t sum = 0u;
        for (uint32_t replicate = 0u;
             replicate < kPanelReplicates;
             ++replicate)
        {
            const uint32_t n = InvocationsPerSlot(shape.K, replicate);
            if (n < 2u) return false;
            sum += n;
        }
        if (sum != TotalInvocations(shape.K)) return false;
    }
    return DescriptorJson("unknown", NativeSystematicEmission::EquationEvaluation).empty() &&
        !IsKnownScope("unknown") &&
        std::strcmp(PrimaryMetric(kFreshRepairScope), kFreshRepairMetric) == 0 &&
        std::strcmp(NestedMetric(kFreshRepairScope),
                    kFreshRepairNestedMetric) == 0;
}

bool RunScreen(int target_cpu)
{
    if (!wirehair_wh2_bench::NativePanelPlatformSupported()) {
        std::cerr << "broad systematic emission screen requires Linux affinity\n";
        return false;
    }
    const std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::now() +
        std::chrono::seconds(kInternalDeadlineSeconds);

    std::vector<std::shared_ptr<const ScreenCell> > cells;
    for (size_t index = 0u;
         index < sizeof(kCellShapes) / sizeof(kCellShapes[0]);
         ++index)
    {
        std::shared_ptr<ScreenCell> cell(new ScreenCell);
        if (!BuildCell(kCellShapes[index], *cell) ||
            !CellSealsMatch(index, *cell))
        {
            std::cerr << "cannot build exact broad systematic cell\n";
            return false;
        }
        cells.push_back(cell);
    }
    std::vector<DescriptorSeal> descriptors;
    if (!BuildDescriptorSeals(descriptors)) {
        std::cerr << "invalid broad systematic descriptor seal\n";
        return false;
    }
    EmitConfig(target_cpu, cells, descriptors);

    static const char* const scopes[] = {
        kFreshSystematicScope,
        kFreshRepairScope,
        kSteadyRepairScope
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
                for (const Comparison& comparison : kComparisons)
                {
                    const char* const metric = PrimaryMetric(scope);
                    const DescriptorSeal* const left_descriptor =
                        FindDescriptor(descriptors, metric, comparison.Left);
                    const DescriptorSeal* const right_descriptor =
                        FindDescriptor(descriptors, metric, comparison.Right);
                    const bool has_nested_metric =
                        std::strcmp(scope, kFreshRepairScope) == 0;
                    const DescriptorSeal* const left_nested_descriptor =
                        has_nested_metric ? FindDescriptor(
                            descriptors,
                            kFreshRepairNestedMetric,
                            comparison.Left) : nullptr;
                    const DescriptorSeal* const right_nested_descriptor =
                        has_nested_metric ? FindDescriptor(
                            descriptors,
                            kFreshRepairNestedMetric,
                            comparison.Right) : nullptr;
                    const std::string panel_key_sha256 = PanelKeySha256(
                        cell->Shape, comparison.Name, scope);
                    const NativePanelOrder order = PanelOrder(
                        cell->Shape, comparison.Name, scope, replicate);
                    if (!left_descriptor || !right_descriptor ||
                        (has_nested_metric &&
                         (!left_nested_descriptor ||
                          !right_nested_descriptor)) ||
                        panel_key_sha256.size() != 64u ||
                        (order != NativePanelOrder::ABBA &&
                         order != NativePanelOrder::BAAB))
                    {
                        std::cerr << "invalid broad systematic panel key\n";
                        return false;
                    }
                    if (std::chrono::steady_clock::now() >= deadline) {
                        std::cerr << "broad systematic internal deadline\n";
                        return false;
                    }
                    std::shared_ptr<NestedTimingCollector> left_nested(
                        new NestedTimingCollector);
                    std::shared_ptr<NestedTimingCollector> right_nested(
                        new NestedTimingCollector);
                    const NativePanelArm left = MakePanelArm(
                        cell,
                        scope,
                        comparison.Left,
                        left_descriptor->Sha256,
                        left_nested_descriptor ?
                            left_nested_descriptor->Sha256 : std::string(),
                        left_nested);
                    const NativePanelArm right = MakePanelArm(
                        cell,
                        scope,
                        comparison.Right,
                        right_descriptor->Sha256,
                        right_nested_descriptor ?
                            right_nested_descriptor->Sha256 : std::string(),
                        right_nested);
                    const NativePanelResult panel =
                        wirehair_wh2_bench::ExecuteNativeTimingPanel(
                            target_cpu,
                            order,
                            invocations,
                            left,
                            right);
                    if (panel.Status != NativePanelStatus::Complete ||
                        !panel.PanelComparable || panel.IsFatal() ||
                        !panel.HasLeftPreflight ||
                        !panel.HasRightPreflight ||
                        panel.LeftPreflight.Disposition !=
                            NativePanelDisposition::Success ||
                        panel.RightPreflight.Disposition !=
                            NativePanelDisposition::Success ||
                        panel.LeftPreflight.OutcomeCode != Wirehair_Success ||
                        panel.RightPreflight.OutcomeCode != Wirehair_Success)
                    {
                        std::cerr << "broad systematic panel failed: "
                                  << wirehair_wh2_bench::NativePanelStatusName(
                                         panel.Status)
                                  << ' ' << panel.Diagnostic << '\n';
                        return false;
                    }
                    const bool systematic =
                        std::strcmp(scope, kFreshSystematicScope) == 0;
                    const uint32_t expected_left_direct = systematic &&
                            comparison.Left ==
                                NativeSystematicEmission::RetainedSourceDirect ?
                            cell->Shape.K : 0u;
                    const uint32_t expected_right_direct = systematic &&
                            comparison.Right ==
                                NativeSystematicEmission::RetainedSourceDirect ?
                            cell->Shape.K : 0u;
                    if (!panel.LeftPreflight.HasDecodedExtra ||
                        !panel.RightPreflight.HasDecodedExtra ||
                        panel.LeftPreflight.DecodedExtra !=
                            expected_left_direct ||
                        panel.RightPreflight.DecodedExtra !=
                            expected_right_direct)
                    {
                        std::cerr << "broad systematic preflight receipt drifted\n";
                        return false;
                    }
                    std::array<uint64_t, 8> nested_slots = {{ 0u }};
                    if (std::strcmp(scope, kFreshRepairScope) == 0)
                    {
                        if (!CollectNestedSlots(
                                order,
                                invocations,
                                *left_nested,
                                *right_nested,
                                nested_slots))
                        {
                            std::cerr << "fresh-repair nested receipt drifted\n";
                            return false;
                        }
                    }
                    else if (!left_nested->Elapsed.empty() ||
                             !right_nested->Elapsed.empty())
                    {
                        std::cerr << "unexpected nested timing receipt\n";
                        return false;
                    }
                    for (const wirehair_wh2_bench::NativePanelSlot& slot :
                         panel.Slots)
                    {
                        if (!slot.HasElapsedNanoseconds ||
                            slot.ElapsedNanoseconds == 0u)
                        {
                            std::cerr << "broad systematic panel has null time\n";
                            return false;
                        }
                    }
                    EmitPanel(
                        *cell,
                        comparison.Name,
                        invocations,
                        left_descriptor->Sha256,
                        left_nested_descriptor ?
                            left_nested_descriptor->Sha256 : std::string(),
                        order,
                        panel_key_sha256,
                        replicate,
                        right_descriptor->Sha256,
                        right_nested_descriptor ?
                            right_nested_descriptor->Sha256 : std::string(),
                        scope,
                        metric,
                        panel,
                        nested_slots);
                    ++panel_count;
                }
            }
        }
    }
    if (std::chrono::steady_clock::now() >= deadline) {
        std::cerr << "broad systematic internal deadline\n";
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
            std::cerr << "WH2 broad systematic emission selftest failed\n";
            return 1;
        }
        std::cout << "WH2 broad systematic emission selftest passed\n";
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
