#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#if !defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
#error "contract worker requires counter-free benchmark equations"
#elif defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "contract worker must not compile with generic WH2 test hooks"
#endif

#include "Wh2FrozenTiming.h"
#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"
#include "Wh2NativePanel.h"
#include "Wh2PhaseAttribution.h"

#include "../codec/WirehairV2RawSeed.h"

#include <wirehair/wirehair.h>

#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#endif

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#define WIREHAIR_WH2_SOURCE_GIT_COMMIT \
    "0000000000000000000000000000000000000000"
#endif

namespace {

using wirehair::wh2_benchmark::FrozenPacketTrace;
using wirehair::wh2_benchmark::FrozenPanelKind;
using wirehair::wh2_benchmark::FrozenRecoveryCell;
using wirehair::wh2_benchmark::FrozenScheduleName;
using wirehair::wh2_benchmark::FrozenTimingBaseCell;
using wirehair::wh2_benchmark::FrozenTimingCell;
using wirehair::wh2_benchmark::FrozenTimingOrder;
using wirehair::wh2_benchmark::FrozenTimingPanel;
using wirehair::wh2_benchmark::FrozenTimingScope;
using wirehair::wh2_benchmark::FrozenTraceStatus;
using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativePanelArm;
using wirehair_wh2_bench::NativePanelDisposition;
using wirehair_wh2_bench::NativePanelInvocation;
using wirehair_wh2_bench::NativePanelInvocationResult;
using wirehair_wh2_bench::NativePanelOrder;
using wirehair_wh2_bench::NativePanelResult;
using wirehair_wh2_bench::NativePanelSide;
using wirehair_wh2_bench::NativePanelStatus;
using wirehair_wh2_bench::NativeReceiveFixture;
using wirehair_wh2_bench::NativeSolveFixture;
using wirehair_wh2_bench::NativeTimingControlProbe;
using wirehair_wh2_bench::NativeTimingControlQualification;
using wirehair_wh2_bench::NativeTimingControlQualificationResult;
using wirehair_wh2_bench::PhaseMeasuredObservation;
using wirehair_wh2_bench::PhaseObservationCollector;
using wirehair_wh2_bench::PhasePanelAssembly;
using wirehair_wh2_bench::PhaseSlotTotals;
using wirehair_wh2_bench::PhaseSolveArm;
using wirehair_wh2_bench::PhaseSolveInvocation;
using wirehair_wh2_bench::PhaseSolveObservation;
using wirehair_wh2_bench::RecoveryCellResult;
using wirehair_wh2_bench::RecoveryOutcome;

static const char kTraceSchema[] =
    "wirehair.wh2.native-trace-record.v1";
static const char kRecoverySchema[] =
    "wirehair.wh2.native-recovery-record.v1";
static const char kRawRecoverySchema[] =
    "wirehair.wh2.native-recovery-record.v3";
static const char kTimingSchema[] =
    "wirehair.wh2.native-timing-record.v4";
static const char kTimingQualificationSchema[] =
    "wirehair.wh2.native-timing-qualification-record.v1";
static const char kPhaseAttributionSchema[] =
    "wirehair.wh2.native-phase-attribution-record.v1";
static const char kPhaseTraceSchema[] =
    "wirehair.wh2.native-phase-attribution-trace-record.v1";
static const char kPhaseCellSchema[] =
    "wirehair.wh2.phase-attribution-cell.v1";
static const char kDescriptionSchema[] =
    "wirehair.wh2.native-worker-description.v1";
static const char kRawDescriptionSchema[] =
    "wirehair.wh2.native-worker-description.v3";
static const char kDescriptorSchema[] =
    "wirehair.wh2.native-arm-descriptor.v1";
static const char kWorkSchema[] = "wirehair.wh2.native-work.v1";
static const char kRealizedSchema[] =
    "wirehair.wh2.realized-construction.v1";
static const char kRawRealizedSchema[] =
    "wirehair.wh2.raw-realized-construction.v2";
static const char kTimingProxyWitnessSchema[] =
    "wirehair.wh2.native-timing-proxy-witness.v2";
static const char kTimingProxyDomain[] =
    "wirehair.wh2.timing-proxy-domain.v2\0";
static const char kTimingProxyConfigurationDomain[] =
    "wirehair.wh2.production-timing-proxy-configuration.v2\0";
static const char kTimingProxyEquationDomain[] =
    "wirehair.wh2.production-timing-proxy-equations.v2\0";
static const char kTimingProxyNormalizedStructureDomain[] =
    "wirehair.wh2.timing-proxy-normalized-structure.v2\0";
static const char kRealizedDomain[] =
    "wirehair.wh2.realized-construction.v1\0";
static const char kSourceDomain[] = "wirehair.wh2.source.v1\0";
static const char kTimingQualificationMessageDomain[] =
    "wirehair.wh2.timing-qualification-message.v1\0";
static const char kZeroSha256[] =
    "0000000000000000000000000000000000000000000000000000000000000000";
static const uint64_t kWh2ProfileId = UINT64_C(0x4b295bbb47f4f9c9);
static const uint16_t kWh2ProfileEncodingVersion = 1u;
static const uint32_t kRecoveryOverheadCap = 4u;
static const uint32_t kTimingReceiveOverheadCap = 256u;
static const uint32_t kPhaseCoordinateCount = 24u;
static const uint32_t kPhaseProfileCount = 2u;
static const char kRawControlArm[] =
    "wirehair2_raw_d12_h12_periodic";
static const char kRawControlDescriptorSha256[] =
    "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11";
static const char kTimingCandidateArm[] =
    "wirehair2_dense_two07_basis_v1";
static const char kTimingCandidateTransform[] =
    "dense-anchor-two07+basis-segment-adjacent-symdiff-v1";
static const char kTimingCandidateDescriptorSha256[] =
    "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388";
static const char kProductionSeedBasis[] = "production-profile-v1";
static const char kNotApplicableSeedBasis[] = "not-applicable";
static const char kDenseAnchorDisabled[] = "disabled";
static const char kDenseAnchorTwo07[] = "two07";
static const char kDenseAnchorNotApplicable[] = "not-applicable";

struct ArmDescriptor
{
    const char* Arm;
    const char* Codec;
    std::string CanonicalDescriptor;
    std::string DescriptorSha256;
    const char* ConstructionSeedBasis;
    const char* SeedScheduleSha256;
    const char* DenseAnchorLayout;
};

struct CandidateDefinition
{
    const char* Id;
    const char* Arm;
    const char* Transform;
    const char* DescriptorSha256;
    uint32_t DenseRows;
    uint32_t HeavyRows;
    wirehair_v2::HeavyCoefficientFamily HeavyFamily;
    wirehair_v2::DenseAnchorLayout DenseAnchors;
};

static const CandidateDefinition kRecoveryCandidates[] = {
    {
        "d12-h11-periodic",
        "wirehair2_raw_d12_h11_periodic",
        "d12-h11-periodic",
        "91d7c1a558e1cf93b002fcf2062b7657d301faca03972215495bdf2429499e90",
        12u,
        11u,
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy,
        wirehair_v2::DenseAnchorLayout::Disabled
    },
    {
        "d12-h13-periodic",
        "wirehair2_raw_d12_h13_periodic",
        "d12-h13-periodic",
        "7c7889747a97ac160726b807fb03349344d49d4bec84c9e8220aa4689b00d2ca",
        12u,
        13u,
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy,
        wirehair_v2::DenseAnchorLayout::Disabled
    },
    {
        "d13-h12-periodic",
        "wirehair2_raw_d13_h12_periodic",
        "d13-h12-periodic",
        "c70e0f57bb8d7783fa29b0decbed5da5058a8eb532d57d540f72108e114f091a",
        13u,
        12u,
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy,
        wirehair_v2::DenseAnchorLayout::Disabled
    },
    {
        "two07",
        kTimingCandidateArm,
        kTimingCandidateTransform,
        kTimingCandidateDescriptorSha256,
        12u,
        12u,
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy,
        wirehair_v2::DenseAnchorLayout::Two07
    }
};

struct WorkerContext
{
    int Cpu;
    int Pid;
    uint64_t ProcessStartTicks;
    std::string BinarySha256;
    std::vector<ArmDescriptor> Arms;
    const CandidateDefinition* Candidate;
};

struct TimingQualificationProbeCache
{
    NativeTimingControlProbe Probe;
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    uint32_t ConstructionAttempt = 0u;
    bool Initialized = false;
};

bool EmitLine(const std::string& line);

bool IsLowerSha256(const std::string& value)
{
    if (value.size() != 64u) {
        return false;
    }
    for (std::size_t i = 0u; i < value.size(); ++i)
    {
        const char ch = value[i];
        if (!((ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f'))) {
            return false;
        }
    }
    return true;
}

bool IsLowerGitCommit(const std::string& value)
{
    if (value.size() != 40u) {
        return false;
    }
    for (std::size_t i = 0u; i < value.size(); ++i)
    {
        const char ch = value[i];
        if (!((ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f'))) {
            return false;
        }
    }
    return true;
}

std::string HexSeed(uint64_t value)
{
    static const char hex[] = "0123456789abcdef";
    std::string text("0x0000000000000000");
    for (unsigned i = 0u; i < 16u; ++i)
    {
        text[17u - i] = hex[value & 15u];
        value >>= 4;
    }
    return text;
}

std::string HexPacketSeed(uint32_t value)
{
    static const char hex[] = "0123456789abcdef";
    std::string text("0x00000000");
    for (unsigned i = 0u; i < 8u; ++i)
    {
        text[9u - i] = hex[value & 15u];
        value >>= 4;
    }
    return text;
}

void AppendUint64LittleEndian(uint64_t value, std::string& output)
{
    for (unsigned i = 0u; i < 8u; ++i) {
        output += static_cast<char>(value >> (8u * i));
    }
}

void AppendUint32LittleEndian(uint32_t value, std::string& output)
{
    for (unsigned i = 0u; i < 4u; ++i) {
        output += static_cast<char>(value >> (8u * i));
    }
}

void AppendUint16LittleEndian(uint16_t value, std::string& output)
{
    output += static_cast<char>(value);
    output += static_cast<char>(value >> 8u);
}

bool ParseHexByte(char high, char low, uint8_t& value)
{
    const auto nibble = [](char ch) -> int {
        if (ch >= '0' && ch <= '9') return ch - '0';
        if (ch >= 'a' && ch <= 'f') return ch - 'a' + 10;
        return -1;
    };
    const int h = nibble(high);
    const int l = nibble(low);
    if (h < 0 || l < 0) {
        return false;
    }
    value = static_cast<uint8_t>((h << 4) | l);
    return true;
}

bool Sha256HexBytes(const std::string& hex, std::string& bytes)
{
    if (!IsLowerSha256(hex)) {
        return false;
    }
    std::string decoded;
    decoded.reserve(32u);
    for (std::size_t i = 0u; i < hex.size(); i += 2u)
    {
        uint8_t byte = 0u;
        if (!ParseHexByte(hex[i], hex[i + 1u], byte)) {
            return false;
        }
        decoded += static_cast<char>(byte);
    }
    bytes.swap(decoded);
    return true;
}

bool SourceSeedFromCellJson(const std::string& cell_json, uint64_t& seed)
{
    const std::string input = std::string(
        kSourceDomain, sizeof(kSourceDomain) - 1u) + cell_json;
    const std::string digest = wirehair::wh2_benchmark::Sha256Hex(input);
    if (!IsLowerSha256(digest)) {
        return false;
    }
    uint64_t result = 0u;
    for (unsigned i = 0u; i < 8u; ++i)
    {
        uint8_t byte = 0u;
        if (!ParseHexByte(digest[2u * i], digest[2u * i + 1u], byte)) {
            return false;
        }
        result = (result << 8u) | byte;
    }
    seed = result;
    return true;
}

std::string TimingQualificationMessageSha256(
    const std::string& base_cell_json)
{
    if (base_cell_json.empty()) {
        return std::string();
    }
    const std::string input = std::string(
        kTimingQualificationMessageDomain,
        sizeof(kTimingQualificationMessageDomain) - 1u) + base_cell_json;
    const std::string digest = wirehair::wh2_benchmark::Sha256Hex(input);
    return IsLowerSha256(digest) ? digest : std::string();
}

bool ReadSelfBinary(std::string& bytes)
{
#if defined(__linux__)
    std::ifstream input("/proc/self/exe", std::ios::binary);
    if (!input) {
        return false;
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.eof() && input.fail()) {
        return false;
    }
    bytes = contents.str();
    return !bytes.empty();
#else
    (void)bytes;
    return false;
#endif
}

bool SelfBinarySha256(std::string& digest)
{
    std::string binary;
    if (!ReadSelfBinary(binary)) {
        return false;
    }
    digest = wirehair::wh2_benchmark::Sha256Hex(binary);
    return IsLowerSha256(digest);
}

std::string CanonicalDescriptor(
    const char* arm,
    const char* codec,
    const char* equation_transform)
{
    std::string json;
    json.reserve(192u);
    json += "{\"arm\":\"";
    json += arm;
    json += "\",\"codec\":\"";
    json += codec;
    json += "\",\"equation_transform\":\"";
    json += equation_transform;
    json += "\",\"schema\":\"";
    json += kDescriptorSchema;
    json += "\"}";
    return json;
}

bool BuildArmDescriptors(std::vector<ArmDescriptor>& arms)
{
    static const char* const names[] = {
        "wirehair2_head", "wirehair1", kTimingCandidateArm
    };
    static const char* const codecs[] = {
        "wirehair2_certified", "wirehair1", "wirehair2_experiment"
    };
    static const char* const transforms[] = {
        "none", "none", kTimingCandidateTransform
    };
    std::vector<ArmDescriptor> built;
    built.reserve(3u);
    for (std::size_t i = 0u; i < 3u; ++i)
    {
        ArmDescriptor descriptor;
        descriptor.Arm = names[i];
        descriptor.Codec = codecs[i];
        descriptor.CanonicalDescriptor = CanonicalDescriptor(
            names[i], codecs[i], transforms[i]);
        descriptor.DescriptorSha256 = wirehair::wh2_benchmark::Sha256Hex(
            descriptor.CanonicalDescriptor);
        descriptor.ConstructionSeedBasis = i == 1u ?
            kNotApplicableSeedBasis : kProductionSeedBasis;
        descriptor.SeedScheduleSha256 = kZeroSha256;
        descriptor.DenseAnchorLayout = i == 1u ?
            kDenseAnchorNotApplicable :
            (i == 2u ? kDenseAnchorTwo07 : kDenseAnchorDisabled);
        if (!IsLowerSha256(descriptor.DescriptorSha256) ||
            (i == 2u && descriptor.DescriptorSha256 !=
                kTimingCandidateDescriptorSha256))
        {
            return false;
        }
        built.push_back(descriptor);
    }
    arms.swap(built);
    return true;
}

const CandidateDefinition* FindRecoveryCandidate(const std::string& id)
{
    for (std::size_t i = 0u;
         i < sizeof(kRecoveryCandidates) / sizeof(kRecoveryCandidates[0]);
         ++i)
    {
        if (id == kRecoveryCandidates[i].Id) {
            return &kRecoveryCandidates[i];
        }
    }
    return nullptr;
}

bool ClosedCandidateSemantics(const CandidateDefinition& candidate)
{
    const std::string id(candidate.Id);
    const std::string arm(candidate.Arm);
    const std::string transform(candidate.Transform);
    const bool h11 = id == "d12-h11-periodic" &&
        arm == "wirehair2_raw_d12_h11_periodic" &&
        transform == id && candidate.DenseRows == 12u &&
        candidate.HeavyRows == 11u &&
        candidate.DenseAnchors == wirehair_v2::DenseAnchorLayout::Disabled;
    const bool h13 = id == "d12-h13-periodic" &&
        arm == "wirehair2_raw_d12_h13_periodic" &&
        transform == id && candidate.DenseRows == 12u &&
        candidate.HeavyRows == 13u &&
        candidate.DenseAnchors == wirehair_v2::DenseAnchorLayout::Disabled;
    const bool d13 = id == "d13-h12-periodic" &&
        arm == "wirehair2_raw_d13_h12_periodic" &&
        transform == id && candidate.DenseRows == 13u &&
        candidate.HeavyRows == 12u &&
        candidate.DenseAnchors == wirehair_v2::DenseAnchorLayout::Disabled;
    const bool two07 = id == "two07" &&
        arm == kTimingCandidateArm &&
        transform == kTimingCandidateTransform &&
        candidate.DenseRows == 12u && candidate.HeavyRows == 12u &&
        candidate.DenseAnchors == wirehair_v2::DenseAnchorLayout::Two07;
    return (h11 || h13 || d13 || two07) &&
        candidate.HeavyFamily ==
            wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy;
}

bool ValidateRecoveryCandidateDefinition(
    const CandidateDefinition& candidate);

bool ValidRawSeedSchedule()
{
    static const bool valid = []() -> bool {
        const std::string digest = wirehair::wh2_benchmark::Sha256Hex(
            wirehair_v2::test::kRawArchitectureSeedScheduleCanonical);
        if (!IsLowerSha256(digest) || digest !=
            wirehair_v2::test::kRawArchitectureSeedScheduleSha256)
        {
            return false;
        }
        static const uint32_t attempts[] = { 0u, 1u, 2u, 255u };
        static const uint64_t expected_precode[] = {
            UINT64_C(0x487468302aad7105),
            UINT64_C(0xe6abe1e9a9f7ed1a),
            UINT64_C(0x84e35ba32942692f),
            UINT64_C(0xe1b6a7f5f5df09f0)
        };
        static const uint32_t expected_packet[] = {
            UINT32_C(0x4ec72102),
            UINT32_C(0xecfe9abb),
            UINT32_C(0x8b361474),
            UINT32_C(0xe8096049)
        };
        wirehair_v2::PrecodeParams params;
        wirehair_v2::PacketRowConfig packet;
        wirehair_v2::test::MakeRawArchitectureConfiguration(
            3u, params, packet);
        if (params.BlockCount != 3u || params.DenseRows != 12u ||
            params.HeavyRows != 12u ||
            params.HeavyFamily !=
                wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
            packet.MixCount != wirehair_v2::kCertifiedPacketMixCount)
        {
            return false;
        }
        for (std::size_t i = 0u;
             i < sizeof(attempts) / sizeof(attempts[0]); ++i)
        {
            if (wirehair_v2::PrecodeParamsForAttempt(
                    params, attempts[i]).Seed != expected_precode[i] ||
                wirehair_v2::PacketConfigForAttempt(
                    packet, attempts[i]).PeelSeed != expected_packet[i])
            {
                return false;
            }
        }
        return true;
    }();
    return valid;
}

std::string RawTransformName(const char* structure)
{
    if (!structure || !*structure || !ValidRawSeedSchedule()) {
        return std::string();
    }
    return std::string(structure);
}

const char* DenseAnchorLayoutName(wirehair_v2::DenseAnchorLayout layout)
{
    switch (layout)
    {
    case wirehair_v2::DenseAnchorLayout::Disabled:
        return kDenseAnchorDisabled;
    case wirehair_v2::DenseAnchorLayout::Two07:
        return kDenseAnchorTwo07;
    default:
        return nullptr;
    }
}

bool BuildCandidateArmDescriptors(
    const CandidateDefinition& candidate,
    std::vector<ArmDescriptor>& arms)
{
    if (FindRecoveryCandidate(candidate.Id) != &candidate ||
        !ClosedCandidateSemantics(candidate) ||
        !ValidRawSeedSchedule() ||
        !ValidateRecoveryCandidateDefinition(candidate))
    {
        std::cerr << "wirehair_wh2_contract_worker: invalid raw recovery "
            "candidate definition\n";
        return false;
    }

    std::vector<ArmDescriptor> built;
    if (!BuildArmDescriptors(built) || built.size() != 3u) {
        std::cerr << "wirehair_wh2_contract_worker: cannot build base arm "
            "descriptors\n";
        return false;
    }
    ArmDescriptor raw_control;
    raw_control.Arm = kRawControlArm;
    raw_control.Codec = "wirehair2_experiment";
    const std::string raw_control_transform = RawTransformName(
        "d12-h12-periodic");
    raw_control.CanonicalDescriptor = CanonicalDescriptor(
        raw_control.Arm, raw_control.Codec, raw_control_transform.c_str());
    raw_control.DescriptorSha256 = wirehair::wh2_benchmark::Sha256Hex(
        raw_control.CanonicalDescriptor);
    raw_control.ConstructionSeedBasis =
        wirehair_v2::test::kRawArchitectureSeedBasis;
    raw_control.SeedScheduleSha256 =
        wirehair_v2::test::kRawArchitectureSeedScheduleSha256;
    raw_control.DenseAnchorLayout = kDenseAnchorDisabled;

    ArmDescriptor descriptor;
    descriptor.Arm = candidate.Arm;
    descriptor.Codec = "wirehair2_experiment";
    const std::string candidate_transform = RawTransformName(
        candidate.Transform);
    descriptor.CanonicalDescriptor = CanonicalDescriptor(
        candidate.Arm, descriptor.Codec, candidate_transform.c_str());
    descriptor.DescriptorSha256 = wirehair::wh2_benchmark::Sha256Hex(
        descriptor.CanonicalDescriptor);
    descriptor.ConstructionSeedBasis =
        wirehair_v2::test::kRawArchitectureSeedBasis;
    descriptor.SeedScheduleSha256 =
        wirehair_v2::test::kRawArchitectureSeedScheduleSha256;
    descriptor.DenseAnchorLayout =
        candidate.DenseAnchors == wirehair_v2::DenseAnchorLayout::Two07 ?
            kDenseAnchorTwo07 : kDenseAnchorDisabled;
    if (raw_control_transform.empty() || candidate_transform.empty() ||
        !IsLowerSha256(raw_control.DescriptorSha256) ||
        raw_control.DescriptorSha256 != kRawControlDescriptorSha256 ||
        !IsLowerSha256(descriptor.DescriptorSha256) ||
        descriptor.DescriptorSha256 != candidate.DescriptorSha256 ||
        raw_control.DescriptorSha256 == descriptor.DescriptorSha256 ||
        descriptor.DescriptorSha256 == built[0].DescriptorSha256 ||
        descriptor.DescriptorSha256 == built[1].DescriptorSha256 ||
        raw_control.DescriptorSha256 == built[0].DescriptorSha256 ||
        raw_control.DescriptorSha256 == built[1].DescriptorSha256)
    {
        std::cerr << "wirehair_wh2_contract_worker: raw recovery arm "
            "descriptor mismatch\n";
        return false;
    }
    for (std::size_t i = 0u;
         i < sizeof(kRecoveryCandidates) / sizeof(kRecoveryCandidates[0]);
         ++i)
    {
        const CandidateDefinition& other = kRecoveryCandidates[i];
        if (&other == &candidate) continue;
        if (std::string(other.Id) == candidate.Id ||
            std::string(other.Arm) == candidate.Arm ||
            std::string(other.DescriptorSha256) ==
                candidate.DescriptorSha256)
        {
            std::cerr << "wirehair_wh2_contract_worker: duplicate raw "
                "recovery candidate identity\n";
            return false;
        }
    }
    built[2] = raw_control;
    built.push_back(descriptor);
    arms.swap(built);
    return true;
}

bool TwoAnchorBasisTransform(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig&,
    void*)
{
    if (block_count != params.BlockCount || block_bytes == 0u ||
        params.DenseAnchors != wirehair_v2::DenseAnchorLayout::Disabled)
    {
        return false;
    }
    params.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    return true;
}

bool CanonicalStructureTransformInput(
    uint32_t block_count,
    const wirehair_v2::PrecodeParams& params,
    const wirehair_v2::PacketRowConfig& packet)
{
    const wirehair_v2::PrecodeParams expected =
        wirehair_v2::MakeCertifiedParams(block_count, params.Seed);
    return ValidRawSeedSchedule() &&
        params.BlockCount == expected.BlockCount &&
        params.Staircase == expected.Staircase &&
        params.DenseRows == expected.DenseRows &&
        params.HeavyRows == expected.HeavyRows &&
        params.SourceHits == expected.SourceHits &&
        params.DenseIdentityCorner == expected.DenseIdentityCorner &&
        params.DenseAnchors == expected.DenseAnchors &&
        params.HeavyFamily == expected.HeavyFamily &&
        params.Seed == 0u && packet.PeelSeed == 0u &&
        packet.MixCount == wirehair_v2::kCertifiedPacketMixCount;
}

bool RawBaselineTransform(
    uint32_t block_count,
    uint32_t,
    wirehair_v2::PrecodeParams& base_params,
    wirehair_v2::PacketRowConfig& base_packet_config,
    void*)
{
    if (!CanonicalStructureTransformInput(
            block_count, base_params, base_packet_config))
    {
        return false;
    }
    wirehair_v2::PrecodeParams transformed;
    wirehair_v2::PacketRowConfig transformed_packet;
    wirehair_v2::test::MakeRawArchitectureConfiguration(
        block_count, transformed, transformed_packet);
    base_params = transformed;
    base_packet_config = transformed_packet;
    return true;
}

bool RecoveryCandidateTransform(
    uint32_t block_count,
    uint32_t,
    wirehair_v2::PrecodeParams& base_params,
    wirehair_v2::PacketRowConfig& base_packet_config,
    void* context)
{
    const CandidateDefinition* const candidate =
        static_cast<const CandidateDefinition*>(context);
    if (!candidate || FindRecoveryCandidate(candidate->Id) != candidate ||
        !ClosedCandidateSemantics(*candidate) ||
        !CanonicalStructureTransformInput(
            block_count, base_params, base_packet_config))
    {
        return false;
    }

    const wirehair_v2::PrecodeParams original_params = base_params;
    const wirehair_v2::PacketRowConfig original_packet = base_packet_config;
    wirehair_v2::PrecodeParams transformed;
    wirehair_v2::PacketRowConfig transformed_packet;
    wirehair_v2::test::MakeRawArchitectureConfiguration(
        block_count, transformed, transformed_packet);
    transformed.DenseRows = candidate->DenseRows;
    transformed.HeavyRows = candidate->HeavyRows;
    transformed.HeavyFamily = candidate->HeavyFamily;
    transformed.DenseAnchors = candidate->DenseAnchors;
    if (transformed.BlockCount != original_params.BlockCount ||
        transformed.Staircase != original_params.Staircase ||
        transformed.SourceHits != original_params.SourceHits ||
        transformed.DenseIdentityCorner !=
            original_params.DenseIdentityCorner ||
        transformed.DenseAnchors != candidate->DenseAnchors ||
        transformed.Seed != wirehair_v2::test::kRawArchitecturePrecodeSeed ||
        transformed_packet.PeelSeed !=
            wirehair_v2::test::kRawArchitecturePacketSeed ||
        transformed_packet.MixCount != original_packet.MixCount)
    {
        return false;
    }
    base_params = transformed;
    base_packet_config = transformed_packet;
    return true;
}

bool ValidateRecoveryCandidateDefinition(
    const CandidateDefinition& candidate)
{
    if (!ClosedCandidateSemantics(candidate)) {
        return false;
    }
    static const uint32_t block_counts[] = {
        3u, 4u, 1001u, 3061u, 5550u, 7533u
    };
    static const uint32_t block_byte_values[] = {
        2u, 64u, 1280u, 4096u
    };
    static const uint32_t attempts[] = { 0u, 1u, 2u, 255u };
    static const uint64_t expected_precode[] = {
        UINT64_C(0x487468302aad7105),
        UINT64_C(0xe6abe1e9a9f7ed1a),
        UINT64_C(0x84e35ba32942692f),
        UINT64_C(0xe1b6a7f5f5df09f0)
    };
    static const uint32_t expected_packet[] = {
        UINT32_C(0x4ec72102),
        UINT32_C(0xecfe9abb),
        UINT32_C(0x8b361474),
        UINT32_C(0xe8096049)
    };
    for (std::size_t k_index = 0u;
         k_index < sizeof(block_counts) / sizeof(block_counts[0]); ++k_index)
    {
        const uint32_t K = block_counts[k_index];
        wirehair_v2::PrecodeParams certified;
        wirehair_v2::PacketRowConfig certified_packet;
        certified = wirehair_v2::MakeCertifiedParams(K, 0u);
        certified_packet = wirehair_v2::PacketRowConfig();
        if (certified.BlockCount != K ||
            certified.DenseRows != 12u || certified.HeavyRows != 12u ||
            certified.HeavyFamily !=
                wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
            certified.Seed != 0u || certified_packet.PeelSeed != 0u ||
            certified_packet.MixCount !=
                wirehair_v2::kCertifiedPacketMixCount)
        {
            return false;
        }

        for (std::size_t byte_index = 0u;
             byte_index < sizeof(block_byte_values) /
                 sizeof(block_byte_values[0]); ++byte_index)
        {
            const uint32_t block_bytes = block_byte_values[byte_index];
            for (std::size_t attempt_index = 0u;
                 attempt_index < sizeof(attempts) / sizeof(attempts[0]);
                 ++attempt_index)
            {
                const uint32_t attempt = attempts[attempt_index];
                const NativeArmSpec raw_control_spec =
                    wirehair_wh2_bench::MakeExperimentalWh2Arm(
                        attempt,
                        &RawBaselineTransform,
                        nullptr,
                        nullptr,
                        wirehair_wh2_bench::NativeWh2BaseKind::
                            CanonicalCertifiedStructure);
                const NativeArmSpec candidate_spec =
                    wirehair_wh2_bench::MakeExperimentalWh2Arm(
                        attempt,
                        &RecoveryCandidateTransform,
                        const_cast<CandidateDefinition*>(&candidate),
                        nullptr,
                        wirehair_wh2_bench::NativeWh2BaseKind::
                            CanonicalCertifiedStructure);
                wirehair_wh2_bench::ResolvedNativeWh2Configuration control;
                wirehair_wh2_bench::ResolvedNativeWh2Configuration actual;
                if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
                        raw_control_spec, K, block_bytes, control) ||
                    !wirehair_wh2_bench::ResolveNativeWh2Configuration(
                        candidate_spec, K, block_bytes, actual) ||
                    control.PrecodeAttempt != attempt ||
                    control.PacketAttempt != attempt ||
                    actual.PrecodeAttempt != attempt ||
                    actual.PacketAttempt != attempt ||
                    control.Params.Seed != expected_precode[attempt_index] ||
                    actual.Params.Seed != expected_precode[attempt_index] ||
                    control.PacketConfig.PeelSeed !=
                        expected_packet[attempt_index] ||
                    actual.PacketConfig.PeelSeed !=
                        expected_packet[attempt_index] ||
                    control.Params.BlockCount != K ||
                    actual.Params.BlockCount != K ||
                    control.Params.Staircase != certified.Staircase ||
                    actual.Params.Staircase != certified.Staircase ||
                    control.Params.SourceHits != certified.SourceHits ||
                    actual.Params.SourceHits != certified.SourceHits ||
                    control.Params.DenseIdentityCorner !=
                        certified.DenseIdentityCorner ||
                    actual.Params.DenseIdentityCorner !=
                        certified.DenseIdentityCorner ||
                    control.Params.DenseAnchors !=
                        wirehair_v2::DenseAnchorLayout::Disabled ||
                    actual.Params.DenseAnchors != candidate.DenseAnchors ||
                    control.Params.DenseRows != 12u ||
                    control.Params.HeavyRows != 12u ||
                    control.Params.HeavyFamily !=
                        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
                    actual.Params.DenseRows != candidate.DenseRows ||
                    actual.Params.HeavyRows != candidate.HeavyRows ||
                    actual.Params.HeavyFamily != candidate.HeavyFamily ||
                    control.PacketConfig.MixCount !=
                        wirehair_v2::kCertifiedPacketMixCount ||
                    actual.PacketConfig.MixCount !=
                        wirehair_v2::kCertifiedPacketMixCount)
                {
                    return false;
                }
            }
        }

        // Exercise each fail-closed base guard independently and require a
        // rejected transform to leave both inputs exactly unchanged.
        for (unsigned mutation = 0u; mutation < 10u; ++mutation)
        {
            wirehair_v2::PrecodeParams unexpected = certified;
            if (mutation == 0u) unexpected.DenseRows = 11u;
            else if (mutation == 1u) unexpected.HeavyRows = 11u;
            else if (mutation == 2u) unexpected.HeavyFamily =
                    wirehair_v2::HeavyCoefficientFamily::HashedNonzero;
            else if (mutation == 3u) ++unexpected.Staircase;
            else if (mutation == 4u) unexpected.SourceHits = 1u;
            else if (mutation == 5u) unexpected.DenseIdentityCorner = true;
            else if (mutation == 7u) unexpected.Seed = 1u;
            else if (mutation == 9u) unexpected.DenseAnchors =
                    wirehair_v2::DenseAnchorLayout::Two07;
            wirehair_v2::PacketRowConfig packet = certified_packet;
            if (mutation == 6u) packet.MixCount = 2u;
            else if (mutation == 8u) packet.PeelSeed = 1u;
            const wirehair_v2::PrecodeParams before = unexpected;
            const wirehair_v2::PacketRowConfig before_packet = packet;
            if (RecoveryCandidateTransform(
                    K, 64u, unexpected, packet,
                    const_cast<CandidateDefinition*>(&candidate)) ||
                unexpected.BlockCount != before.BlockCount ||
                unexpected.Staircase != before.Staircase ||
                unexpected.DenseRows != before.DenseRows ||
                unexpected.HeavyRows != before.HeavyRows ||
                unexpected.SourceHits != before.SourceHits ||
                unexpected.DenseIdentityCorner != before.DenseIdentityCorner ||
                unexpected.DenseAnchors != before.DenseAnchors ||
                unexpected.HeavyFamily != before.HeavyFamily ||
                unexpected.Seed != before.Seed ||
                packet.PeelSeed != before_packet.PeelSeed ||
                packet.MixCount != before_packet.MixCount)
            {
                return false;
            }
        }
    }
    return true;
}

bool ArmSpecFor(
    std::size_t arm_index,
    uint32_t construction_attempt,
    const CandidateDefinition* candidate,
    NativeArmSpec& spec)
{
    switch (arm_index)
    {
    case 0u:
        spec = wirehair_wh2_bench::MakeCertifiedWh2Arm(
            construction_attempt);
        return true;
    case 1u:
        if (construction_attempt != 0u) {
            return false;
        }
        spec = wirehair_wh2_bench::MakeWirehair1Arm();
        return true;
    case 2u:
        if (candidate)
        {
            spec = wirehair_wh2_bench::MakeExperimentalWh2Arm(
                construction_attempt,
                &RawBaselineTransform,
                nullptr,
                nullptr,
                wirehair_wh2_bench::NativeWh2BaseKind::
                    CanonicalCertifiedStructure);
        }
        else
        {
            spec = wirehair_wh2_bench::MakeExperimentalWh2Arm(
                construction_attempt, &TwoAnchorBasisTransform);
        }
        return true;
    case 3u:
        if (!candidate || FindRecoveryCandidate(candidate->Id) != candidate ||
            !ClosedCandidateSemantics(*candidate))
        {
            return false;
        }
        spec = wirehair_wh2_bench::MakeExperimentalWh2Arm(
            construction_attempt,
            &RecoveryCandidateTransform,
            const_cast<CandidateDefinition*>(candidate),
            nullptr,
            wirehair_wh2_bench::NativeWh2BaseKind::
                CanonicalCertifiedStructure);
        return true;
    default:
        return false;
    }
}

uint32_t ConstructionAttemptFor(std::size_t arm_index, uint32_t base_attempt)
{
    return arm_index == 1u ? 0u : base_attempt;
}

bool ValidateResolvedProductionArm(
    std::size_t arm_index,
    const NativeArmSpec& spec,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t construction_attempt)
{
    if (arm_index == 1u) {
        return spec.Kind == wirehair_wh2_bench::NativeArmKind::Wirehair1 &&
            construction_attempt == 0u;
    }
    if (arm_index > 2u) {
        return false;
    }

    wirehair_wh2_bench::ResolvedNativeWh2Configuration actual;
    wirehair_wh2_bench::ResolvedNativeWh2Configuration head;
    if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
            spec, K, block_bytes, actual) ||
        !wirehair_wh2_bench::ResolveNativeWh2Configuration(
            wirehair_wh2_bench::MakeCertifiedWh2Arm(
                construction_attempt),
            K, block_bytes, head))
    {
        return false;
    }

    const wirehair_v2::DenseAnchorLayout expected_layout = arm_index == 2u ?
        wirehair_v2::DenseAnchorLayout::Two07 :
        wirehair_v2::DenseAnchorLayout::Disabled;
    const wirehair_v2::PrecodeParams& a = actual.Params;
    const wirehair_v2::PrecodeParams& h = head.Params;
    return actual.PrecodeAttempt == construction_attempt &&
        actual.PacketAttempt == construction_attempt &&
        head.PrecodeAttempt == construction_attempt &&
        head.PacketAttempt == construction_attempt &&
        a.BlockCount == h.BlockCount && a.Staircase == h.Staircase &&
        a.DenseRows == h.DenseRows && a.HeavyRows == h.HeavyRows &&
        a.SourceHits == h.SourceHits &&
        a.DenseIdentityCorner == h.DenseIdentityCorner &&
        a.HeavyFamily == h.HeavyFamily && a.Seed == h.Seed &&
        h.DenseAnchors == wirehair_v2::DenseAnchorLayout::Disabled &&
        a.DenseAnchors == expected_layout &&
        actual.PacketConfig.PeelSeed == head.PacketConfig.PeelSeed &&
        actual.PacketConfig.MixCount == head.PacketConfig.MixCount;
}

bool ValidateProductionTimingArmDefinitions()
{
    static const uint32_t block_counts[] = { 3u, 1001u, 64000u };
    static const uint32_t block_bytes_values[] = { 64u, 1280u };
    static const uint32_t attempts[] = { 0u, 1u, 2u, 255u };
    for (uint32_t K : block_counts)
    {
        for (uint32_t block_bytes : block_bytes_values)
        {
            for (uint32_t attempt : attempts)
            {
                NativeArmSpec head_spec;
                NativeArmSpec candidate_spec;
                if (!ArmSpecFor(0u, attempt, nullptr, head_spec) ||
                    !ArmSpecFor(2u, attempt, nullptr, candidate_spec) ||
                    !ValidateResolvedProductionArm(
                        0u, head_spec, K, block_bytes, attempt) ||
                    !ValidateResolvedProductionArm(
                        2u, candidate_spec, K, block_bytes, attempt))
                {
                    return false;
                }
            }
        }
    }

    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(1001u, UINT64_C(0x5eed));
    wirehair_v2::PacketRowConfig packet;
    const wirehair_v2::PrecodeParams original = params;
    const wirehair_v2::PacketRowConfig original_packet = packet;
    if (!TwoAnchorBasisTransform(
            1001u, 64u, params, packet, nullptr) ||
        params.DenseAnchors != wirehair_v2::DenseAnchorLayout::Two07 ||
        params.BlockCount != original.BlockCount ||
        params.Staircase != original.Staircase ||
        params.DenseRows != original.DenseRows ||
        params.HeavyRows != original.HeavyRows ||
        params.SourceHits != original.SourceHits ||
        params.DenseIdentityCorner != original.DenseIdentityCorner ||
        params.HeavyFamily != original.HeavyFamily ||
        params.Seed != original.Seed ||
        packet.PeelSeed != original_packet.PeelSeed ||
        packet.MixCount != original_packet.MixCount)
    {
        return false;
    }
    wirehair_v2::PrecodeSystem transformed_system;
    if (!wirehair_v2::BuildPrecodeSystem(params, transformed_system) ||
        transformed_system.Params.DenseAnchors !=
            wirehair_v2::DenseAnchorLayout::Two07)
    {
        return false;
    }
    params = original;
    packet = original_packet;
    if (TwoAnchorBasisTransform(1000u, 64u, params, packet, nullptr) ||
        params.DenseAnchors != original.DenseAnchors ||
        packet.PeelSeed != original_packet.PeelSeed ||
        packet.MixCount != original_packet.MixCount)
    {
        return false;
    }
    params = original;
    packet = original_packet;
    if (TwoAnchorBasisTransform(1001u, 0u, params, packet, nullptr) ||
        params.DenseAnchors != original.DenseAnchors ||
        packet.PeelSeed != original_packet.PeelSeed ||
        packet.MixCount != original_packet.MixCount)
    {
        return false;
    }
    params = original;
    packet = original_packet;
    params.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    const wirehair_v2::PrecodeParams invalid_before = params;
    if (TwoAnchorBasisTransform(1001u, 64u, params, packet, nullptr) ||
        params.BlockCount != invalid_before.BlockCount ||
        params.Staircase != invalid_before.Staircase ||
        params.DenseRows != invalid_before.DenseRows ||
        params.HeavyRows != invalid_before.HeavyRows ||
        params.SourceHits != invalid_before.SourceHits ||
        params.DenseIdentityCorner != invalid_before.DenseIdentityCorner ||
        params.DenseAnchors != invalid_before.DenseAnchors ||
        params.HeavyFamily != invalid_before.HeavyFamily ||
        params.Seed != invalid_before.Seed ||
        packet.PeelSeed != original_packet.PeelSeed ||
        packet.MixCount != original_packet.MixCount)
    {
        return false;
    }
    return true;
}

// This is deliberately a structure-only timing transform.  It verifies that
// the D12/H12/periodic/anchors-disabled structure descriptor is a no-op when
// instantiated under the production timing seed policy.  It does not apply
// the uniform-raw recovery seed schedule, whose graph/peel/rank/solve work is
// distinct and is bound separately in the witness.
bool ProductionTimingD12StructureTransform(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig& packet,
    void*)
{
    if (block_count != params.BlockCount || block_bytes == 0u ||
        params.Staircase == 0u || params.DenseRows != 12u ||
        params.HeavyRows != 12u || params.SourceHits == 0u ||
        params.DenseIdentityCorner ||
        params.DenseAnchors != wirehair_v2::DenseAnchorLayout::Disabled ||
        params.HeavyFamily !=
            wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
        packet.MixCount != wirehair_v2::kCertifiedPacketMixCount)
    {
        return false;
    }
    return true;
}

bool SameResolvedConfiguration(
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& left,
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& right)
{
    const wirehair_v2::PrecodeParams& a = left.Params;
    const wirehair_v2::PrecodeParams& b = right.Params;
    return left.PrecodeAttempt == right.PrecodeAttempt &&
        left.PacketAttempt == right.PacketAttempt &&
        a.BlockCount == b.BlockCount && a.Staircase == b.Staircase &&
        a.DenseRows == b.DenseRows && a.HeavyRows == b.HeavyRows &&
        a.SourceHits == b.SourceHits &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.HeavyFamily == b.HeavyFamily &&
        a.DenseAnchors == b.DenseAnchors && a.Seed == b.Seed &&
        left.PacketConfig.PeelSeed == right.PacketConfig.PeelSeed &&
        left.PacketConfig.MixCount == right.PacketConfig.MixCount;
}

bool SameResolvedNonSeedStructure(
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& left,
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& right)
{
    const wirehair_v2::PrecodeParams& a = left.Params;
    const wirehair_v2::PrecodeParams& b = right.Params;
    return left.PrecodeAttempt == right.PrecodeAttempt &&
        left.PacketAttempt == right.PacketAttempt &&
        a.BlockCount == b.BlockCount && a.Staircase == b.Staircase &&
        a.DenseRows == b.DenseRows && a.HeavyRows == b.HeavyRows &&
        a.SourceHits == b.SourceHits &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.HeavyFamily == b.HeavyFamily &&
        a.DenseAnchors == b.DenseAnchors &&
        left.PacketConfig.MixCount == right.PacketConfig.MixCount;
}

std::string TimingProxyNormalizedStructureBytes(
    uint32_t block_bytes,
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& resolved)
{
    const wirehair_v2::PrecodeParams& params = resolved.Params;
    std::string bytes(
        kTimingProxyNormalizedStructureDomain,
        sizeof(kTimingProxyNormalizedStructureDomain) - 1u);
    bytes.reserve(bytes.size() + 56u);
    AppendUint32LittleEndian(params.BlockCount, bytes);
    AppendUint32LittleEndian(block_bytes, bytes);
    AppendUint32LittleEndian(resolved.PrecodeAttempt, bytes);
    AppendUint32LittleEndian(resolved.PacketAttempt, bytes);
    AppendUint32LittleEndian(params.Staircase, bytes);
    AppendUint32LittleEndian(params.DenseRows, bytes);
    AppendUint32LittleEndian(params.HeavyRows, bytes);
    AppendUint32LittleEndian(params.SourceHits, bytes);
    AppendUint32LittleEndian(params.DenseIdentityCorner ? 1u : 0u, bytes);
    AppendUint32LittleEndian(
        static_cast<uint32_t>(params.HeavyFamily), bytes);
    AppendUint32LittleEndian(
        static_cast<uint32_t>(params.DenseAnchors), bytes);
    AppendUint32LittleEndian(resolved.PacketConfig.MixCount, bytes);
    return bytes;
}

std::string TimingProxyConfigurationBytes(
    uint32_t block_bytes,
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& resolved)
{
    const wirehair_v2::PrecodeParams& params = resolved.Params;
    std::string bytes(
        kTimingProxyConfigurationDomain,
        sizeof(kTimingProxyConfigurationDomain) - 1u);
    bytes.reserve(bytes.size() + 64u);
    AppendUint32LittleEndian(params.BlockCount, bytes);
    AppendUint32LittleEndian(block_bytes, bytes);
    AppendUint32LittleEndian(resolved.PrecodeAttempt, bytes);
    AppendUint32LittleEndian(resolved.PacketAttempt, bytes);
    AppendUint32LittleEndian(params.Staircase, bytes);
    AppendUint32LittleEndian(params.DenseRows, bytes);
    AppendUint32LittleEndian(params.HeavyRows, bytes);
    AppendUint32LittleEndian(params.SourceHits, bytes);
    AppendUint32LittleEndian(params.DenseIdentityCorner ? 1u : 0u, bytes);
    AppendUint32LittleEndian(
        static_cast<uint32_t>(params.HeavyFamily), bytes);
    AppendUint32LittleEndian(
        static_cast<uint32_t>(params.DenseAnchors), bytes);
    AppendUint64LittleEndian(params.Seed, bytes);
    AppendUint32LittleEndian(resolved.PacketConfig.PeelSeed, bytes);
    AppendUint32LittleEndian(resolved.PacketConfig.MixCount, bytes);
    return bytes;
}

bool TimingProxyEquationSha256(
    uint32_t block_bytes,
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& resolved,
    std::string& digest)
{
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(resolved.Params, system) ||
        !wirehair_v2::ValidatePrecodeSystem(system))
    {
        return false;
    }
    std::string bytes(
        kTimingProxyEquationDomain,
        sizeof(kTimingProxyEquationDomain) - 1u);
    const std::string configuration =
        TimingProxyConfigurationBytes(block_bytes, resolved);
    bytes.reserve(bytes.size() + configuration.size() +
        system.StaircaseRows.size() * 16u +
        system.DenseBasisRowColumns.size() * 16u);
    AppendUint32LittleEndian(
        static_cast<uint32_t>(configuration.size()), bytes);
    bytes += configuration;
    const auto append_rows = [&bytes](
        const std::vector<std::vector<uint32_t>>& rows)
    {
        AppendUint32LittleEndian(static_cast<uint32_t>(rows.size()), bytes);
        for (const std::vector<uint32_t>& row : rows)
        {
            AppendUint32LittleEndian(static_cast<uint32_t>(row.size()), bytes);
            for (uint32_t column : row) {
                AppendUint32LittleEndian(column, bytes);
            }
        }
    };
    append_rows(system.StaircaseRows);
    append_rows(system.DenseBasisRowColumns);
    digest = wirehair::wh2_benchmark::Sha256Hex(bytes);
    return IsLowerSha256(digest);
}

bool EmitTimingProxyWitness(
    const std::string& binary_sha256,
    const std::string& source_git_commit,
    const std::vector<ArmDescriptor>& arms)
{
    const std::string raw_reference_descriptor = CanonicalDescriptor(
        kRawControlArm, "wirehair2_experiment", "d12-h12-periodic");
    if (arms.size() != 3u || !IsLowerSha256(binary_sha256) ||
        !IsLowerGitCommit(source_git_commit) ||
        arms[0].Arm != std::string("wirehair2_head") ||
        arms[0].DescriptorSha256.empty() ||
        arms[0].ConstructionSeedBasis != std::string(kProductionSeedBasis) ||
        arms[0].SeedScheduleSha256 != std::string(kZeroSha256) ||
        arms[2].Arm != std::string(kTimingCandidateArm) ||
        arms[2].DescriptorSha256 != kTimingCandidateDescriptorSha256 ||
        arms[2].ConstructionSeedBasis != std::string(kProductionSeedBasis) ||
        arms[2].SeedScheduleSha256 != std::string(kZeroSha256) ||
        !ValidRawSeedSchedule() ||
        wirehair::wh2_benchmark::Sha256Hex(raw_reference_descriptor) !=
            kRawControlDescriptorSha256)
    {
        return false;
    }
    const std::vector<FrozenTimingBaseCell> base_cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentTimingBaseCells();
    if (base_cells.size() != 2304u) {
        return false;
    }
    std::vector<FrozenTimingBaseCell> cells;
    for (const FrozenTimingBaseCell& base : base_cells)
    {
        bool duplicate = false;
        for (const FrozenTimingBaseCell& cell : cells) {
            duplicate = duplicate ||
                (cell.K == base.K && cell.block_bytes == base.block_bytes &&
                 cell.base_seed_attempt == base.base_seed_attempt);
        }
        if (!duplicate) cells.push_back(base);
    }
    if (cells.size() != 24u) {
        return false;
    }

    std::string domain_json = "[";
    std::string cell_json = "[";
    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        const FrozenTimingBaseCell& cell = cells[i];
        const NativeArmSpec head =
            wirehair_wh2_bench::MakeCertifiedWh2Arm(
                cell.base_seed_attempt);
        const NativeArmSpec proxy =
            wirehair_wh2_bench::MakeExperimentalWh2Arm(
                cell.base_seed_attempt,
                &ProductionTimingD12StructureTransform);
        const NativeArmSpec raw_reference =
            wirehair_wh2_bench::MakeExperimentalWh2Arm(
                cell.base_seed_attempt, &RawBaselineTransform, nullptr,
                nullptr,
                wirehair_wh2_bench::NativeWh2BaseKind::
                    CanonicalCertifiedStructure);
        wirehair_wh2_bench::ResolvedNativeWh2Configuration head_resolved;
        wirehair_wh2_bench::ResolvedNativeWh2Configuration proxy_resolved;
        wirehair_wh2_bench::ResolvedNativeWh2Configuration raw_resolved;
        if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
                head, cell.K, cell.block_bytes, head_resolved) ||
            !wirehair_wh2_bench::ResolveNativeWh2Configuration(
                proxy, cell.K, cell.block_bytes, proxy_resolved) ||
            !wirehair_wh2_bench::ResolveNativeWh2Configuration(
                raw_reference, cell.K, cell.block_bytes, raw_resolved) ||
            !SameResolvedConfiguration(head_resolved, proxy_resolved) ||
            !SameResolvedNonSeedStructure(head_resolved, raw_resolved) ||
            head_resolved.Params.Seed == raw_resolved.Params.Seed ||
            head_resolved.PacketConfig.PeelSeed ==
                raw_resolved.PacketConfig.PeelSeed)
        {
            return false;
        }
        wirehair_v2::PrecodeSystem head_system;
        wirehair_v2::PrecodeSystem proxy_system;
        if (!wirehair_v2::BuildPrecodeSystem(
                head_resolved.Params, head_system) ||
            !wirehair_v2::BuildPrecodeSystem(
                proxy_resolved.Params, proxy_system) ||
            head_system.StaircaseRows != proxy_system.StaircaseRows ||
            head_system.DenseBasisRowColumns !=
                proxy_system.DenseBasisRowColumns)
        {
            return false;
        }
        const std::string configuration_sha256 =
            wirehair::wh2_benchmark::Sha256Hex(
                TimingProxyConfigurationBytes(
                    cell.block_bytes, head_resolved));
        const std::string normalized_structure_sha256 =
            wirehair::wh2_benchmark::Sha256Hex(
                TimingProxyNormalizedStructureBytes(
                    cell.block_bytes, head_resolved));
        std::string equation_sha256;
        if (!IsLowerSha256(configuration_sha256) ||
            !IsLowerSha256(normalized_structure_sha256) ||
            normalized_structure_sha256 !=
                wirehair::wh2_benchmark::Sha256Hex(
                    TimingProxyNormalizedStructureBytes(
                        cell.block_bytes, raw_resolved)) ||
            !TimingProxyEquationSha256(
                cell.block_bytes, head_resolved, equation_sha256))
        {
            return false;
        }
        if (i != 0u) {
            domain_json += ',';
            cell_json += ',';
        }
        const std::string coordinate =
            std::string("{\"K\":") + std::to_string(cell.K) +
            ",\"block_bytes\":" + std::to_string(cell.block_bytes) +
            ",\"construction_attempt\":" +
            std::to_string(cell.base_seed_attempt) + "}";
        domain_json += coordinate;
        cell_json += std::string("{\"K\":") + std::to_string(cell.K) +
            ",\"block_bytes\":" + std::to_string(cell.block_bytes) +
            ",\"construction_attempt\":" +
            std::to_string(cell.base_seed_attempt) +
            ",\"normalized_structure_sha256\":\"" +
            normalized_structure_sha256 +
            "\",\"production_timing_configuration_sha256\":\"" +
            configuration_sha256 +
            "\",\"production_timing_equation_system_sha256\":\"" +
            equation_sha256 +
            "\",\"production_timing_packet_seed\":\"" +
            HexPacketSeed(head_resolved.PacketConfig.PeelSeed) +
            "\",\"production_timing_precode_seed\":\"" +
            HexSeed(head_resolved.Params.Seed) +
            "\",\"raw_recovery_packet_seed\":\"" +
            HexPacketSeed(raw_resolved.PacketConfig.PeelSeed) +
            "\",\"raw_recovery_precode_seed\":\"" +
            HexSeed(raw_resolved.Params.Seed) +
            "\",\"seeds_differ\":true}";
    }
    domain_json += ']';
    cell_json += ']';
    const std::string domain_input = std::string(
        kTimingProxyDomain, sizeof(kTimingProxyDomain) - 1u) + domain_json;
    const std::string domain_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(domain_input);
    if (!IsLowerSha256(domain_sha256)) {
        return false;
    }
    std::string json;
    json.reserve(cell_json.size() + 800u);
    json += "{\"applicability\":\"development-attempt-0-only-new-witness-required-for-other-attempt-semantics\",\"binary_sha256\":\"";
    json += binary_sha256;
    json += "\",\"cells\":";
    json += cell_json;
    json += ",\"construction_attempts\":[0],\"evidence_phase\":\"development\"";
    json += ",\"production_timing_proxy_arm\":\"wirehair2_head\"";
    json += ",\"production_timing_proxy_arm_descriptor_sha256\":\"";
    json += arms[0].DescriptorSha256;
    json += "\",\"proof_scope\":\"d12-disabled-structure-under-production-timing-seed-policy-v1";
    json += "\",\"raw_recovery_reference_arm\":\"";
    json += kRawControlArm;
    json += "\",\"raw_recovery_reference_arm_descriptor_sha256\":\"";
    json += kRawControlDescriptorSha256;
    json += "\",\"raw_recovery_seed_basis\":\"";
    json += wirehair_v2::test::kRawArchitectureSeedBasis;
    json += "\",\"raw_recovery_seed_schedule_sha256\":\"";
    json += wirehair_v2::test::kRawArchitectureSeedScheduleSha256;
    json += "\",\"schema\":\"";
    json += kTimingProxyWitnessSchema;
    json += "\",\"seed_relationship\":\"raw-recovery-precode-and-packet-seeds-differ-from-production-timing";
    json += "\",\"source_git_commit\":\"";
    json += source_git_commit;
    json += "\",\"timing_candidate_arm\":\"";
    json += kTimingCandidateArm;
    json += "\",\"timing_candidate_arm_descriptor_sha256\":\"";
    json += kTimingCandidateDescriptorSha256;
    json += "\",\"timing_seed_basis\":\"";
    json += kProductionSeedBasis;
    json += "\",\"timing_seed_policy_arms\":[\"wirehair2_head\",\"";
    json += kTimingCandidateArm;
    json += "\"],\"timing_seed_schedule_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"witness_domain_sha256\":\"";
    json += domain_sha256;
    json += "\"}";
    return EmitLine(json);
}

bool RealizedConstructionSha256(
    const ArmDescriptor& arm,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t construction_attempt,
    std::string& digest)
{
    if (arm.Codec == std::string("wirehair2_certified"))
    {
        std::string descriptor_bytes;
        if (!Sha256HexBytes(arm.DescriptorSha256, descriptor_bytes)) {
            return false;
        }
        const uint64_t message_bytes =
            static_cast<uint64_t>(K) * block_bytes;
        std::string profile;
        profile.reserve(32u);
        profile += "WHV2";
        AppendUint16LittleEndian(kWh2ProfileEncodingVersion, profile);
        AppendUint16LittleEndian(32u, profile);
        AppendUint64LittleEndian(kWh2ProfileId, profile);
        AppendUint64LittleEndian(message_bytes, profile);
        AppendUint32LittleEndian(block_bytes, profile);
        profile += static_cast<char>(construction_attempt);
        profile.append(3u, '\0');
        if (profile.size() != 32u) {
            return false;
        }
        const std::string input = std::string(
            kRealizedDomain, sizeof(kRealizedDomain) - 1u) +
            descriptor_bytes + profile;
        digest = wirehair::wh2_benchmark::Sha256Hex(input);
    }
    else
    {
        std::string json;
        json.reserve(288u);
        json += "{\"K\":";
        json += std::to_string(K);
        json += ",\"arm_descriptor_sha256\":\"";
        json += arm.DescriptorSha256;
        json += "\",\"block_bytes\":";
        json += std::to_string(block_bytes);
        json += ",\"codec\":\"";
        json += arm.Codec;
        json += "\",\"construction_attempt\":";
        json += std::to_string(construction_attempt);
        json += ",\"schema\":\"";
        json += kRealizedSchema;
        json += "\"}";
        digest = wirehair::wh2_benchmark::Sha256Hex(json);
    }
    return IsLowerSha256(digest);
}

struct PhaseCoordinate
{
    std::size_t Ordinal;
    FrozenTimingBaseCell Base;
    FrozenTimingCell Qualified;
};

bool SameArmDescriptor(
    const ArmDescriptor& actual,
    const ArmDescriptor& expected)
{
    if (!actual.Arm || !actual.Codec || !actual.ConstructionSeedBasis ||
        !actual.SeedScheduleSha256 || !actual.DenseAnchorLayout ||
        !expected.Arm || !expected.Codec ||
        !expected.ConstructionSeedBasis || !expected.SeedScheduleSha256 ||
        !expected.DenseAnchorLayout)
    {
        return false;
    }
    return std::string(actual.Arm) == expected.Arm &&
        std::string(actual.Codec) == expected.Codec &&
        actual.CanonicalDescriptor == expected.CanonicalDescriptor &&
        actual.DescriptorSha256 == expected.DescriptorSha256 &&
        std::string(actual.ConstructionSeedBasis) ==
            expected.ConstructionSeedBasis &&
        std::string(actual.SeedScheduleSha256) ==
            expected.SeedScheduleSha256 &&
        std::string(actual.DenseAnchorLayout) == expected.DenseAnchorLayout;
}

bool ExactPhaseArmRoster(const std::vector<ArmDescriptor>& arms)
{
    std::vector<ArmDescriptor> expected;
    if (arms.size() != 3u || !BuildArmDescriptors(expected) ||
        expected.size() != arms.size())
    {
        return false;
    }
    for (std::size_t i = 0u; i < arms.size(); ++i) {
        if (!SameArmDescriptor(arms[i], expected[i])) return false;
    }
    return true;
}

bool PhaseProfile(
    std::size_t profile_ordinal,
    const char*& name,
    uint32_t& invocations_per_slot)
{
    if (profile_ordinal == 0u)
    {
        name = "n16";
        invocations_per_slot = 16u;
        return true;
    }
    if (profile_ordinal == 1u)
    {
        name = "n24";
        invocations_per_slot = 24u;
        return true;
    }
    name = nullptr;
    invocations_per_slot = 0u;
    return false;
}

bool BuildPhaseCoordinates(std::vector<PhaseCoordinate>& coordinates)
{
    static const uint32_t expected_k[] = {
        8u, 32u, 100u, 128u, 512u, 1000u,
        2048u, 5000u, 8192u, 20000u, 32768u, 64000u
    };
    static const uint32_t expected_widths[] = { 64u, 1280u };
    const std::vector<FrozenTimingBaseCell> bases =
        wirehair::wh2_benchmark::EnumerateDevelopmentTimingBaseCells();
    if (bases.size() != 2304u) {
        return false;
    }

    std::vector<PhaseCoordinate> built;
    built.reserve(kPhaseCoordinateCount);
    for (std::size_t width = 0u; width < 2u; ++width)
    {
        for (std::size_t k_index = 0u; k_index < 12u; ++k_index)
        {
            // Round zero, then width-major/K-major, lane zero.  This is the
            // strict natural-order projection of the existing timing domain.
            const std::size_t source_ordinal =
                width * 96u + k_index * 8u;
            const FrozenTimingBaseCell& base = bases[source_ordinal];
            FrozenTimingCell qualified;
            if (base.ordinal != source_ordinal ||
                base.K != expected_k[k_index] ||
                base.block_bytes != expected_widths[width] ||
                base.phase != "development" ||
                base.band.empty() ||
                base.loss_ppm != 100000u ||
                base.schedule != wirehair::wh2_benchmark::FrozenSchedule::Iid ||
                base.replicate != 0u || base.base_seed_attempt != 0u ||
                base.fixed_received_overhead != kRecoveryOverheadCap ||
                base.receive_overhead_cap != kTimingReceiveOverheadCap ||
                base.interleave_policy !=
                    "self-counterbalanced-repeat-major-v1" ||
                !wirehair::wh2_benchmark::QualifyDevelopmentTimingCell(
                    base, 0u, qualified) ||
                qualified.ordinal != source_ordinal ||
                qualified.loss_retry_offset != 0u ||
                qualified.loss_seed != qualified.base_loss_seed ||
                qualified.base_cell_sha256 !=
                    wirehair::wh2_benchmark::TimingBaseCellSha256(base) ||
                wirehair::wh2_benchmark::TimingSourceIdentitySha256(
                    qualified) != qualified.base_cell_sha256 ||
                !IsLowerSha256(
                    wirehair::wh2_benchmark::TimingCellSha256(qualified)))
            {
                return false;
            }
            PhaseCoordinate coordinate;
            coordinate.Ordinal = built.size();
            coordinate.Base = base;
            coordinate.Qualified = qualified;
            built.push_back(coordinate);
        }
    }
    if (built.size() != kPhaseCoordinateCount) {
        return false;
    }
    coordinates.swap(built);
    return true;
}

bool GeneratePhaseTrace(
    const PhaseCoordinate& coordinate,
    FrozenPacketTrace& trace,
    wirehair::wh2_benchmark::FrozenTimingTraceReceipt& receipt)
{
    return coordinate.Ordinal < kPhaseCoordinateCount &&
        wirehair::wh2_benchmark::GenerateDevelopmentTimingTrace(
            coordinate.Qualified, trace, receipt) ==
            FrozenTraceStatus::Complete &&
        trace.delivered_ids.size() ==
            static_cast<std::size_t>(coordinate.Qualified.K) +
                kTimingReceiveOverheadCap &&
        receipt.cell_sha256 ==
            wirehair::wh2_benchmark::TimingCellSha256(
                coordinate.Qualified) &&
        receipt.trace_sha256 == trace.trace_sha256 &&
        IsLowerSha256(trace.trace_sha256);
}

std::string CanonicalPhaseCellJson(
    const PhaseCoordinate& coordinate,
    std::size_t profile_ordinal,
    const std::vector<ArmDescriptor>& arms,
    const std::string& trace_sha256,
    const std::string& left_realized_sha256,
    const std::string& right_realized_sha256)
{
    const char* profile = nullptr;
    uint32_t invocations_per_slot = 0u;
    const FrozenTimingCell& cell = coordinate.Qualified;
    const std::string qualified_sha256 =
        wirehair::wh2_benchmark::TimingCellSha256(cell);
    if (!PhaseProfile(
            profile_ordinal, profile, invocations_per_slot) ||
        coordinate.Ordinal >= kPhaseCoordinateCount ||
        !ExactPhaseArmRoster(arms) ||
        !IsLowerSha256(trace_sha256) ||
        !IsLowerSha256(left_realized_sha256) ||
        !IsLowerSha256(right_realized_sha256) ||
        cell.base_seed_attempt != 0u || cell.loss_retry_offset != 0u ||
        cell.fixed_received_overhead != kRecoveryOverheadCap ||
        cell.receive_overhead_cap != kTimingReceiveOverheadCap ||
        !IsLowerSha256(cell.base_cell_sha256) ||
        !IsLowerSha256(qualified_sha256))
    {
        return std::string();
    }

    std::string json;
    json.reserve(1500u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"band\":\"";
    json += cell.band;
    json += "\",\"base_loss_seed\":\"";
    json += HexSeed(cell.base_loss_seed);
    json += "\",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json += ",\"construction_attempt\":0,\"coordinate_ordinal\":";
    json += std::to_string(coordinate.Ordinal);
    json += ",\"diagnostic_phase\":\"development\"";
    json += ",\"fixed_received_overhead\":4";
    json +=
        ",\"interleave_policy\":\"self-counterbalanced-repeat-major-v1\"";
    json += ",\"invocations_per_slot\":";
    json += std::to_string(invocations_per_slot);
    json += ",\"left_arm\":\"";
    json += arms[2].Arm;
    json += "\",\"left_arm_descriptor_sha256\":\"";
    json += arms[2].DescriptorSha256;
    json += "\",\"left_realized_construction_sha256\":\"";
    json += left_realized_sha256;
    json += "\",\"left_repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"loss_retry_offset\":0,\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"order\":\"ABBA\",\"panel_kind\":\"ab\"";
    json += ",\"profile\":\"";
    json += profile;
    json += "\",\"profile_ordinal\":";
    json += std::to_string(profile_ordinal);
    json += ",\"receive_overhead_cap\":256,\"replicate\":0";
    json += ",\"right_arm\":\"";
    json += arms[0].Arm;
    json += "\",\"right_arm_descriptor_sha256\":\"";
    json += arms[0].DescriptorSha256;
    json += "\",\"right_realized_construction_sha256\":\"";
    json += right_realized_sha256;
    json += "\",\"right_repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"schedule\":\"iid\",\"schema\":\"";
    json += kPhaseCellSchema;
    json += "\",\"scope\":\"decoder_solve\"";
    json += ",\"source_base_cell_sha256\":\"";
    json += cell.base_cell_sha256;
    json += "\",\"trace_qualified_timing_cell_sha256\":\"";
    json += qualified_sha256;
    json += "\",\"trace_sha256\":\"";
    json += trace_sha256;
    json += "\"}";
    return json;
}

std::string PhaseCellSha256(
    const PhaseCoordinate& coordinate,
    std::size_t profile_ordinal,
    const std::vector<ArmDescriptor>& arms,
    const std::string& trace_sha256,
    const std::string& left_realized_sha256,
    const std::string& right_realized_sha256)
{
    const std::string json = CanonicalPhaseCellJson(
        coordinate, profile_ordinal, arms, trace_sha256,
        left_realized_sha256, right_realized_sha256);
    return json.empty() ? std::string() :
        wirehair::wh2_benchmark::Sha256Hex(json);
}

bool EmitPhaseAttributionTraces(const std::vector<ArmDescriptor>& arms)
{
    std::vector<PhaseCoordinate> coordinates;
    if (!ExactPhaseArmRoster(arms) ||
        !BuildPhaseCoordinates(coordinates) ||
        coordinates.size() != kPhaseCoordinateCount)
    {
        return false;
    }
    for (const PhaseCoordinate& coordinate : coordinates)
    {
        FrozenPacketTrace trace;
        wirehair::wh2_benchmark::FrozenTimingTraceReceipt receipt;
        std::string left_realized;
        std::string right_realized;
        if (!GeneratePhaseTrace(coordinate, trace, receipt) ||
            !RealizedConstructionSha256(
                arms[2], coordinate.Qualified.K,
                coordinate.Qualified.block_bytes, 0u, left_realized) ||
            !RealizedConstructionSha256(
                arms[0], coordinate.Qualified.K,
                coordinate.Qualified.block_bytes, 0u, right_realized))
        {
            return false;
        }
        std::string phase_cells[2];
        for (std::size_t profile = 0u;
             profile < kPhaseProfileCount;
             ++profile)
        {
            phase_cells[profile] = PhaseCellSha256(
                coordinate, profile, arms, trace.trace_sha256,
                left_realized, right_realized);
            if (!IsLowerSha256(phase_cells[profile])) {
                return false;
            }
        }
        const std::string qualified_sha256 =
            wirehair::wh2_benchmark::TimingCellSha256(
                coordinate.Qualified);
        std::string json;
        json.reserve(800u);
        json += "{\"candidate_count\":";
        json += std::to_string(trace.attempted_candidates);
        json += ",\"coordinate_ordinal\":";
        json += std::to_string(coordinate.Ordinal);
        json += ",\"packet_count\":";
        json += std::to_string(trace.delivered_ids.size());
        json += ",\"phase_cell_sha256_by_profile\":[\"";
        json += phase_cells[0];
        json += "\",\"";
        json += phase_cells[1];
        json += "\"],\"schema\":\"";
        json += kPhaseTraceSchema;
        json += "\",\"source_base_cell_sha256\":\"";
        json += coordinate.Qualified.base_cell_sha256;
        json += "\",\"trace_qualified_timing_cell_sha256\":\"";
        json += qualified_sha256;
        json += "\",\"trace_sha256\":\"";
        json += trace.trace_sha256;
        json += "\"}";
        if (!EmitLine(json)) {
            return false;
        }
    }
    return true;
}

bool RawRealizedConstructionSha256(
    const ArmDescriptor& arm,
    uint32_t K,
    uint32_t block_bytes,
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& resolved,
    std::string& digest)
{
    const wirehair_v2::PrecodeParams& params = resolved.Params;
    const wirehair_v2::PacketRowConfig& packet = resolved.PacketConfig;
    const char* const dense_anchor_layout =
        DenseAnchorLayoutName(params.DenseAnchors);
    if (arm.Codec != std::string("wirehair2_experiment") ||
        !arm.ConstructionSeedBasis || !arm.SeedScheduleSha256 ||
        arm.ConstructionSeedBasis != std::string(
            wirehair_v2::test::kRawArchitectureSeedBasis) ||
        arm.SeedScheduleSha256 != std::string(
            wirehair_v2::test::kRawArchitectureSeedScheduleSha256) ||
        !IsLowerSha256(arm.DescriptorSha256) ||
        !IsLowerSha256(arm.SeedScheduleSha256) ||
        params.BlockCount != K ||
        params.HeavyFamily !=
            wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
        !dense_anchor_layout || !arm.DenseAnchorLayout ||
        arm.DenseAnchorLayout != std::string(dense_anchor_layout) ||
        resolved.PrecodeAttempt > 255u || resolved.PacketAttempt > 255u)
    {
        return false;
    }

    std::string json;
    json.reserve(900u);
    json += "{\"K\":";
    json += std::to_string(K);
    json += ",\"arm\":\"";
    json += arm.Arm;
    json += "\",\"arm_descriptor_sha256\":\"";
    json += arm.DescriptorSha256;
    json += "\",\"binary_dense_rows\":";
    json += std::to_string(params.DenseRows);
    json += ",\"block_bytes\":";
    json += std::to_string(block_bytes);
    json += ",\"codec\":\"";
    json += arm.Codec;
    json += "\",\"construction_seed_basis\":\"";
    json += arm.ConstructionSeedBasis;
    json += "\",\"dense_anchor_layout\":\"";
    json += dense_anchor_layout;
    json += "\",\"dense_identity_corner\":";
    json += params.DenseIdentityCorner ? "true" : "false";
    json += ",\"effective_packet_seed\":\"";
    json += HexPacketSeed(packet.PeelSeed);
    json += "\",\"effective_precode_seed\":\"";
    json += HexSeed(params.Seed);
    json += "\",\"gf256_heavy_rows\":";
    json += std::to_string(params.HeavyRows);
    json += ",\"heavy_family\":\"periodic-cauchy\",\"mix_count\":";
    json += std::to_string(packet.MixCount);
    json += ",\"packet_attempt\":";
    json += std::to_string(resolved.PacketAttempt);
    json += ",\"precode_attempt\":";
    json += std::to_string(resolved.PrecodeAttempt);
    json += ",\"schema\":\"";
    json += kRawRealizedSchema;
    json += "\",\"seed_schedule_sha256\":\"";
    json += arm.SeedScheduleSha256;
    json += "\",\"source_hits\":";
    json += std::to_string(params.SourceHits);
    json += ",\"staircase\":";
    json += std::to_string(params.Staircase);
    json += '}';
    digest = wirehair::wh2_benchmark::Sha256Hex(json);
    return IsLowerSha256(digest);
}

bool MonotonicNanoseconds(uint64_t& value)
{
#if defined(__linux__)
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0 || now.tv_sec < 0 ||
        now.tv_nsec < 0 || now.tv_nsec >= 1000000000L)
    {
        return false;
    }
    const uint64_t seconds = static_cast<uint64_t>(now.tv_sec);
    if (seconds > (std::numeric_limits<uint64_t>::max() -
            static_cast<uint64_t>(now.tv_nsec)) / UINT64_C(1000000000))
    {
        return false;
    }
    value = seconds * UINT64_C(1000000000) +
        static_cast<uint64_t>(now.tv_nsec);
    return value <= static_cast<uint64_t>(
        std::numeric_limits<int64_t>::max());
#else
    (void)value;
    return false;
#endif
}

bool VerifySingletonCpu(int target_cpu, int& actual_cpu, std::string& error)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        error = "target CPU is outside cpu_set_t range";
        return false;
    }
    cpu_set_t observed;
    CPU_ZERO(&observed);
    if (sched_getaffinity(0, sizeof(observed), &observed) != 0)
    {
        error = std::string("sched_getaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    unsigned count = 0u;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        count += CPU_ISSET(static_cast<std::size_t>(cpu), &observed) ?
            1u : 0u;
    }
    if (count != 1u ||
        !CPU_ISSET(static_cast<std::size_t>(target_cpu), &observed))
    {
        error = "worker affinity is not the requested singleton CPU";
        return false;
    }
    errno = 0;
    actual_cpu = sched_getcpu();
    if (actual_cpu < 0)
    {
        error = std::string("sched_getcpu failed: ") + std::strerror(errno);
        return false;
    }
    if (actual_cpu != target_cpu)
    {
        error = "worker migrated from its requested singleton CPU";
        return false;
    }
    error.clear();
    return true;
#else
    (void)target_cpu;
    (void)actual_cpu;
    error = "native contract worker requires Linux CPU affinity";
    return false;
#endif
}

bool PinSingletonCpu(int target_cpu, std::string& error)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        error = "target CPU is outside cpu_set_t range";
        return false;
    }
    cpu_set_t selected;
    CPU_ZERO(&selected);
    CPU_SET(static_cast<std::size_t>(target_cpu), &selected);
    if (sched_setaffinity(0, sizeof(selected), &selected) != 0)
    {
        error = std::string("sched_setaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    int actual_cpu = -1;
    return VerifySingletonCpu(target_cpu, actual_cpu, error);
#else
    (void)target_cpu;
    error = "native contract worker requires Linux CPU affinity";
    return false;
#endif
}

bool ParseCanonicalUnsigned(const std::string& text, std::size_t& value)
{
    if (text.empty() || (text.size() > 1u && text[0] == '0')) {
        return false;
    }
    std::size_t parsed = 0u;
    for (std::size_t i = 0u; i < text.size(); ++i)
    {
        const char ch = text[i];
        if (ch < '0' || ch > '9') {
            return false;
        }
        const unsigned digit = static_cast<unsigned>(ch - '0');
        if (parsed > (std::numeric_limits<std::size_t>::max() - digit) /
                10u)
        {
            return false;
        }
        parsed = parsed * 10u + digit;
    }
    value = parsed;
    return true;
}

bool ReadSelfProcessStartTicks(uint64_t& value)
{
#if defined(__linux__)
    std::ifstream input("/proc/self/stat");
    std::string line;
    if (!input || !std::getline(input, line)) {
        return false;
    }
    const std::size_t close = line.rfind(')');
    if (close == std::string::npos || close + 2u > line.size()) {
        return false;
    }
    // The tokens following "(comm) " start at proc(5) field 3.  Starttime
    // is field 22, hence zero-based token 19 in this suffix.
    std::istringstream suffix(line.substr(close + 2u));
    std::string token;
    for (unsigned index = 0u; index <= 19u; ++index)
    {
        if (!(suffix >> token)) {
            return false;
        }
    }
    if (token.empty() || (token.size() > 1u && token[0] == '0')) {
        return false;
    }
    uint64_t parsed = 0u;
    for (std::size_t i = 0u; i < token.size(); ++i)
    {
        const char ch = token[i];
        if (ch < '0' || ch > '9') {
            return false;
        }
        const unsigned digit = static_cast<unsigned>(ch - '0');
        if (parsed > (std::numeric_limits<uint64_t>::max() - digit) / 10u) {
            return false;
        }
        parsed = parsed * 10u + digit;
    }
    value = parsed;
    return value != 0u && value <= static_cast<uint64_t>(
        std::numeric_limits<int64_t>::max());
#else
    (void)value;
    return false;
#endif
}

bool ParseJobLine(
    const std::string& line,
    char& kind,
    std::size_t& cell_ordinal,
    std::size_t& item_index)
{
    if (line.size() < 5u ||
        (line[0] != 'R' && line[0] != 'L' && line[0] != 'T' &&
         line[0] != 'P') ||
        line[1] != ' ')
    {
        return false;
    }
    const std::size_t separator = line.find(' ', 2u);
    if (separator == std::string::npos || separator == 2u ||
        separator + 1u >= line.size() ||
        line.find(' ', separator + 1u) != std::string::npos)
    {
        return false;
    }
    std::size_t cell = 0u;
    std::size_t item = 0u;
    if (!ParseCanonicalUnsigned(line.substr(2u, separator - 2u), cell) ||
        !ParseCanonicalUnsigned(line.substr(separator + 1u), item))
    {
        return false;
    }
    kind = line[0];
    cell_ordinal = cell;
    item_index = item;
    return true;
}

bool ParseCpu(const std::string& text, int& cpu)
{
    std::size_t parsed = 0u;
    if (!ParseCanonicalUnsigned(text, parsed) ||
        parsed > static_cast<std::size_t>(INT_MAX))
    {
        return false;
    }
    cpu = static_cast<int>(parsed);
    return true;
}

bool EmitLine(const std::string& line)
{
    std::cout << line << '\n';
    std::cout.flush();
    return static_cast<bool>(std::cout);
}

std::string TraceRecord(
    std::size_t ordinal,
    const std::string& cell_sha256,
    const FrozenPacketTrace& trace)
{
    std::string json;
    json.reserve(320u);
    json += "{\"candidate_count\":";
    json += std::to_string(trace.attempted_candidates);
    json += ",\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"packet_count\":";
    json += std::to_string(trace.delivered_ids.size());
    json += ",\"schema\":\"";
    json += kTraceSchema;
    json += "\",\"trace_sha256\":\"";
    json += trace.trace_sha256;
    json += "\"}";
    return json;
}

bool EmitRecoveryTraces()
{
    const std::vector<FrozenRecoveryCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentRecoveryCells();
    if (cells.size() != 360u) {
        return false;
    }
    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        FrozenPacketTrace trace;
        if (cells[i].ordinal != i ||
            wirehair::wh2_benchmark::GenerateFrozenPacketTrace(
                cells[i], trace) != FrozenTraceStatus::Complete ||
            trace.delivered_ids.size() !=
                static_cast<std::size_t>(cells[i].K) + 4u)
        {
            return false;
        }
        const std::string cell_sha256 =
            wirehair::wh2_benchmark::RecoveryCellSha256(cells[i]);
        if (!IsLowerSha256(cell_sha256) ||
            !IsLowerSha256(trace.trace_sha256) ||
            !EmitLine(TraceRecord(i, cell_sha256, trace)))
        {
            return false;
        }
    }
    return true;
}

std::string WorkerDescription(
    const std::string& binary_sha256,
    const std::string& source_git_commit,
    const std::vector<ArmDescriptor>& arms)
{
    if (arms.size() != 3u || !IsLowerGitCommit(source_git_commit))
    {
        return std::string();
    }
    std::string json;
    json.reserve(1400u);
    json += "{\"arms\":[";
    for (std::size_t i = 0u; i < arms.size(); ++i)
    {
        if (i != 0u) json += ',';
        json += "{\"arm\":\"";
        json += arms[i].Arm;
        json += "\",\"arm_descriptor_sha256\":\"";
        json += arms[i].DescriptorSha256;
        json += "\",\"codec\":\"";
        json += arms[i].Codec;
        json += "\"}";
    }
    json += "],\"binary_sha256\":\"";
    json += binary_sha256;
    json += "\",\"schema\":\"";
    json += kDescriptionSchema;
    json += "\",\"source_git_commit\":\"";
    json += source_git_commit;
    json += "\"}";
    return json;
}

std::string RawWorkerDescription(
    const std::string& binary_sha256,
    const std::string& source_git_commit,
    const std::vector<ArmDescriptor>& arms)
{
    if (arms.size() != 4u || !IsLowerSha256(binary_sha256) ||
        !IsLowerGitCommit(source_git_commit))
    {
        return std::string();
    }
    std::string json;
    json.reserve(1940u);
    json += "{\"arms\":[";
    for (std::size_t i = 0u; i < arms.size(); ++i)
    {
        const ArmDescriptor& arm = arms[i];
        if (!arm.ConstructionSeedBasis || !arm.SeedScheduleSha256 ||
            !arm.DenseAnchorLayout ||
            !IsLowerSha256(arm.DescriptorSha256) ||
            !IsLowerSha256(arm.SeedScheduleSha256))
        {
            return std::string();
        }
        if (i != 0u) json += ',';
        json += "{\"arm\":\"";
        json += arm.Arm;
        json += "\",\"arm_descriptor_sha256\":\"";
        json += arm.DescriptorSha256;
        json += "\",\"codec\":\"";
        json += arm.Codec;
        json += "\",\"construction_seed_basis\":\"";
        json += arm.ConstructionSeedBasis;
        json += "\",\"dense_anchor_layout\":\"";
        json += arm.DenseAnchorLayout;
        json += "\",\"seed_schedule_sha256\":\"";
        json += arm.SeedScheduleSha256;
        json += "\"}";
    }
    json += "],\"binary_sha256\":\"";
    json += binary_sha256;
    json += "\",\"schema\":\"";
    json += kRawDescriptionSchema;
    json += "\",\"source_git_commit\":\"";
    json += source_git_commit;
    json += "\"}";
    return json;
}

std::string RecoveryWorkSha256(
    const std::string& cell_sha256,
    std::size_t ordinal)
{
    std::string json;
    json += "{\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"evidence_kind\":\"recovery\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"phase\":\"development\",\"schema\":\"";
    json += kWorkSchema;
    json += "\"}";
    return wirehair::wh2_benchmark::Sha256Hex(json);
}

std::string TimingWorkSha256(
    const std::string& cell_sha256,
    std::size_t ordinal)
{
    std::string json;
    json += "{\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"evidence_kind\":\"timing\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"phase\":\"development\",\"schema\":\"";
    json += kWorkSchema;
    json += "\"}";
    return wirehair::wh2_benchmark::Sha256Hex(json);
}

std::string TimingQualificationWorkSha256(
    const std::string& cell_sha256,
    std::size_t ordinal)
{
    std::string json;
    json += "{\"cell_sha256\":\"";
    json += cell_sha256;
    json +=
        "\",\"evidence_kind\":\"timing_qualification\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"phase\":\"development\",\"schema\":\"";
    json += kWorkSchema;
    json += "\"}";
    return wirehair::wh2_benchmark::Sha256Hex(json);
}

std::string PhaseAttributionWorkSha256(
    const std::string& cell_sha256,
    std::size_t ordinal)
{
    std::string json;
    json += "{\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"evidence_kind\":\"phase_attribution\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"phase\":\"development\",\"schema\":\"";
    json += kWorkSchema;
    json += "\"}";
    return wirehair::wh2_benchmark::Sha256Hex(json);
}

bool RecoveryOutcomeJson(
    const RecoveryCellResult& result,
    const char*& outcome,
    bool& has_decoded_extra,
    uint32_t& decoded_extra)
{
    has_decoded_extra = false;
    decoded_extra = 0u;
    switch (result.Outcome)
    {
    case RecoveryOutcome::Success:
        if (result.Result != Wirehair_Success ||
            result.FirstOverhead > kRecoveryOverheadCap)
        {
            return false;
        }
        outcome = "success";
        has_decoded_extra = true;
        decoded_extra = result.FirstOverhead;
        return true;
    case RecoveryOutcome::NeedMoreAtCap:
        if (result.Result != Wirehair_NeedMore) return false;
        outcome = "need_more_at_cap";
        return true;
    case RecoveryOutcome::ConstructFailed:
        if (result.Result != Wirehair_BadDenseSeed &&
            result.Result != Wirehair_BadPeelSeed)
        {
            return false;
        }
        outcome = "construct_failed";
        return true;
    case RecoveryOutcome::Unsupported:
        if (result.Result != Wirehair_BadInput_SmallN &&
            result.Result != Wirehair_BadInput_LargeN &&
            result.Result != Wirehair_UnsupportedPlatform)
        {
            return false;
        }
        outcome = "unsupported";
        return true;
    case RecoveryOutcome::Fatal:
        return false;
    }
    return false;
}

std::string RecoveryPayload(
    const FrozenRecoveryCell& cell,
    const ArmDescriptor& arm,
    const std::string& binary_sha256,
    const std::string& cell_sha256,
    const std::string& trace_sha256,
    uint32_t construction_attempt,
    const std::string& realized_sha256,
    const char* outcome,
    bool has_decoded_extra,
    uint32_t decoded_extra)
{
    std::string json;
    json.reserve(1000u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"arm\":\"";
    json += arm.Arm;
    json += "\",\"arm_descriptor_sha256\":\"";
    json += arm.DescriptorSha256;
    json += "\",\"band\":\"";
    json += cell.band;
    json += "\",\"base_seed_attempt\":";
    json += std::to_string(cell.base_seed_attempt);
    json += ",\"binary_sha256\":\"";
    json += binary_sha256;
    json += "\",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json += ",\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"construction_attempt\":";
    json += std::to_string(construction_attempt);
    json += ",\"decoded_extra\":";
    if (has_decoded_extra) json += std::to_string(decoded_extra);
    else json += "null";
    json += ",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"outcome\":\"";
    json += outcome;
    json += "\",\"overhead_cap\":";
    json += std::to_string(cell.overhead_cap);
    json += ",\"phase\":\"";
    json += cell.phase;
    json += "\",\"realized_construction_sha256\":\"";
    json += realized_sha256;
    json += "\",\"repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"schedule\":\"";
    json += FrozenScheduleName(cell.schedule);
    json += "\",\"trace_sha256\":\"";
    json += trace_sha256;
    json += "\",\"trial\":";
    json += std::to_string(cell.trial);
    json += '}';
    return json;
}

std::string RawRecoveryPayload(
    const FrozenRecoveryCell& cell,
    const ArmDescriptor& arm,
    const std::string& binary_sha256,
    const std::string& cell_sha256,
    const std::string& trace_sha256,
    uint32_t construction_attempt,
    const wirehair_wh2_bench::ResolvedNativeWh2Configuration& resolved,
    const std::string& realized_sha256,
    const char* outcome,
    bool has_decoded_extra,
    uint32_t decoded_extra)
{
    const wirehair_v2::PrecodeParams& params = resolved.Params;
    const wirehair_v2::PacketRowConfig& packet = resolved.PacketConfig;
    const char* const dense_anchor_layout =
        DenseAnchorLayoutName(params.DenseAnchors);
    if (!arm.ConstructionSeedBasis || !arm.SeedScheduleSha256 ||
        !dense_anchor_layout || !arm.DenseAnchorLayout ||
        arm.DenseAnchorLayout != std::string(dense_anchor_layout) ||
        construction_attempt != resolved.PrecodeAttempt ||
        construction_attempt != resolved.PacketAttempt ||
        params.BlockCount != cell.K ||
        params.HeavyFamily !=
            wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy)
    {
        return std::string();
    }

    std::string json;
    json.reserve(1600u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"arm\":\"";
    json += arm.Arm;
    json += "\",\"arm_descriptor_sha256\":\"";
    json += arm.DescriptorSha256;
    json += "\",\"band\":\"";
    json += cell.band;
    json += "\",\"base_seed_attempt\":";
    json += std::to_string(cell.base_seed_attempt);
    json += ",\"binary_dense_rows\":";
    json += std::to_string(params.DenseRows);
    json += ",\"binary_sha256\":\"";
    json += binary_sha256;
    json += "\",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json += ",\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"construction_attempt\":";
    json += std::to_string(construction_attempt);
    json += ",\"construction_seed_basis\":\"";
    json += arm.ConstructionSeedBasis;
    json += "\",\"decoded_extra\":";
    if (has_decoded_extra) json += std::to_string(decoded_extra);
    else json += "null";
    json += ",\"dense_anchor_layout\":\"";
    json += dense_anchor_layout;
    json += '"';
    json += ",\"dense_identity_corner\":";
    json += params.DenseIdentityCorner ? "true" : "false";
    json += ",\"effective_packet_seed\":\"";
    json += HexPacketSeed(packet.PeelSeed);
    json += "\",\"effective_precode_seed\":\"";
    json += HexSeed(params.Seed);
    json += "\",\"gf256_heavy_rows\":";
    json += std::to_string(params.HeavyRows);
    json += ",\"heavy_family\":\"periodic-cauchy\",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"mix_count\":";
    json += std::to_string(packet.MixCount);
    json += ",\"outcome\":\"";
    json += outcome;
    json += "\",\"overhead_cap\":";
    json += std::to_string(cell.overhead_cap);
    json += ",\"packet_attempt\":";
    json += std::to_string(resolved.PacketAttempt);
    json += ",\"phase\":\"";
    json += cell.phase;
    json += "\",\"precode_attempt\":";
    json += std::to_string(resolved.PrecodeAttempt);
    json += ",\"realized_construction_sha256\":\"";
    json += realized_sha256;
    json += "\",\"repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"schedule\":\"";
    json += FrozenScheduleName(cell.schedule);
    json += "\",\"seed_schedule_sha256\":\"";
    json += arm.SeedScheduleSha256;
    json += "\",\"source_hits\":";
    json += std::to_string(params.SourceHits);
    json += ",\"staircase\":";
    json += std::to_string(params.Staircase);
    json += ",\"trace_sha256\":\"";
    json += trace_sha256;
    json += "\",\"trial\":";
    json += std::to_string(cell.trial);
    json += '}';
    return json;
}

std::string Envelope(
    const char* schema,
    std::size_t ordinal,
    int cpu,
    int pid,
    uint64_t started_ns,
    uint64_t finished_ns,
    uint64_t process_start_ticks,
    const std::string& binary_sha256,
    const std::string& message_sha256,
    const std::string& work_sha256,
    const std::string& payload)
{
    std::string json;
    json.reserve(payload.size() + 640u);
    json += "{\"cpu\":";
    json += std::to_string(cpu);
    json += ",\"finished_monotonic_ns\":";
    json += std::to_string(finished_ns);
    json += ",\"message_sha256\":\"";
    json += message_sha256;
    json += "\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"payload\":";
    json += payload;
    json += ",\"schema\":\"";
    json += schema;
    json += "\",\"started_monotonic_ns\":";
    json += std::to_string(started_ns);
    json += ",\"work_sha256\":\"";
    json += work_sha256;
    json += "\",\"worker_binary_sha256\":\"";
    json += binary_sha256;
    json += "\",\"worker_pid\":";
    json += std::to_string(pid);
    json += ",\"worker_process_start_ticks\":";
    json += std::to_string(process_start_ticks);
    json += '}';
    return json;
}

bool RunRecoveryJob(
    const WorkerContext& context,
    std::size_t cell_ordinal,
    std::size_t arm_index,
    std::string& envelope,
    std::string& error)
{
    const std::vector<FrozenRecoveryCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentRecoveryCells();
    if (cells.size() != 360u || cell_ordinal >= cells.size() ||
        arm_index >= context.Arms.size())
    {
        error = "recovery job index is outside the frozen domain";
        return false;
    }
    const FrozenRecoveryCell& cell = cells[cell_ordinal];
    if (cell.ordinal != cell_ordinal || cell.overhead_cap != 4u) {
        error = "recovery cell differs from the frozen domain";
        return false;
    }

    int actual_cpu = -1;
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t started_ns = 0u;
    if (!MonotonicNanoseconds(started_ns)) {
        error = "cannot read CLOCK_MONOTONIC before recovery job";
        return false;
    }

    FrozenPacketTrace trace;
    if (wirehair::wh2_benchmark::GenerateFrozenPacketTrace(cell, trace) !=
            FrozenTraceStatus::Complete ||
        trace.delivered_ids.size() != static_cast<std::size_t>(cell.K) + 4u)
    {
        error = "cannot regenerate complete frozen recovery trace";
        return false;
    }
    const std::string cell_json =
        wirehair::wh2_benchmark::CanonicalRecoveryCellJson(cell);
    const std::string cell_sha256 =
        wirehair::wh2_benchmark::RecoveryCellSha256(cell);
    uint64_t source_seed = 0u;
    std::vector<uint8_t> source;
    if (cell_json.empty() || !IsLowerSha256(cell_sha256) ||
        !IsLowerSha256(trace.trace_sha256) ||
        !SourceSeedFromCellJson(cell_json, source_seed) ||
        !wirehair_wh2_bench::MakeDeterministicSource(
            cell.K, cell.block_bytes, source_seed, source))
    {
        error = "cannot construct deterministic recovery source";
        return false;
    }
    const std::string message_sha256 = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    if (!IsLowerSha256(message_sha256)) {
        error = "cannot hash deterministic recovery source";
        return false;
    }

    const uint32_t attempt = ConstructionAttemptFor(
        arm_index, cell.base_seed_attempt);
    NativeArmSpec spec;
    if (!ArmSpecFor(arm_index, attempt, context.Candidate, spec)) {
        error = "cannot construct selected native recovery arm";
        return false;
    }
    const bool raw_arm = context.Candidate && arm_index >= 2u;
    wirehair_wh2_bench::ResolvedNativeWh2Configuration resolved;
    if ((raw_arm && !wirehair_wh2_bench::ResolveNativeWh2Configuration(
             spec, cell.K, cell.block_bytes, resolved)) ||
        (!raw_arm && !ValidateResolvedProductionArm(
             arm_index, spec, cell.K, cell.block_bytes, attempt)))
    {
        error = raw_arm ?
            "cannot resolve selected raw recovery construction" :
            "recovery arm differs from its production descriptor";
        return false;
    }
    const RecoveryCellResult result = wirehair_wh2_bench::RunRecoveryCell(
        spec, cell.K, cell.block_bytes, source, trace.delivered_ids);
    const char* outcome = nullptr;
    bool has_decoded_extra = false;
    uint32_t decoded_extra = 0u;
    if (!RecoveryOutcomeJson(
            result, outcome, has_decoded_extra, decoded_extra))
    {
        error = "fatal or internally inconsistent recovery result";
        return false;
    }
    std::string realized_sha256;
    const bool realized_ok = raw_arm ?
        RawRealizedConstructionSha256(
            context.Arms[arm_index], cell.K, cell.block_bytes, resolved,
            realized_sha256) :
        RealizedConstructionSha256(
            context.Arms[arm_index], cell.K, cell.block_bytes, attempt,
            realized_sha256);
    if (!realized_ok)
    {
        error = "cannot hash realized recovery construction";
        return false;
    }
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t finished_ns = 0u;
    if (!MonotonicNanoseconds(finished_ns) || finished_ns < started_ns) {
        error = "cannot read CLOCK_MONOTONIC after recovery job";
        return false;
    }
    const std::string payload = raw_arm ?
        RawRecoveryPayload(
            cell, context.Arms[arm_index], context.BinarySha256,
            cell_sha256, trace.trace_sha256, attempt, resolved,
            realized_sha256, outcome, has_decoded_extra, decoded_extra) :
        RecoveryPayload(
            cell, context.Arms[arm_index], context.BinarySha256,
            cell_sha256, trace.trace_sha256, attempt, realized_sha256,
            outcome, has_decoded_extra, decoded_extra);
    if (payload.empty()) {
        error = "cannot serialize recovery payload";
        return false;
    }
    const std::size_t ordinal = cell_ordinal * context.Arms.size() + arm_index;
    const std::string work_sha256 = RecoveryWorkSha256(
        cell_sha256, ordinal);
    if (!IsLowerSha256(work_sha256)) {
        error = "cannot hash recovery work identity";
        return false;
    }
    envelope = Envelope(
        raw_arm ? kRawRecoverySchema : kRecoverySchema,
        ordinal, actual_cpu, context.Pid, started_ns,
        finished_ns, context.ProcessStartTicks, context.BinarySha256,
        message_sha256, work_sha256, payload);
    return true;
}

NativePanelInvocationResult ClassifyTimedResult(
    WirehairResult result,
    uint64_t elapsed_ns,
    bool has_decoded_extra,
    uint32_t decoded_extra,
    bool allow_need_more)
{
    if (result == Wirehair_Success)
    {
        if (elapsed_ns == 0u || elapsed_ns > static_cast<uint64_t>(
                std::numeric_limits<int64_t>::max()))
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                static_cast<int64_t>(Wirehair_Error), false, 0u, 0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            static_cast<int64_t>(result), has_decoded_extra,
            decoded_extra, elapsed_ns);
    }
    switch (result)
    {
    case Wirehair_NeedMore:
        if (!allow_need_more) {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                static_cast<int64_t>(result), false, 0u, 0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::PreflightFailure,
            static_cast<int64_t>(result), false, 0u, 0u);
    case Wirehair_BadDenseSeed:
    case Wirehair_BadPeelSeed:
    case Wirehair_BadInput_SmallN:
    case Wirehair_BadInput_LargeN:
    case Wirehair_UnsupportedPlatform:
        return NativePanelInvocationResult(
            NativePanelDisposition::PreflightFailure,
            static_cast<int64_t>(result), false, 0u, 0u);
    default:
        return NativePanelInvocationResult(
            NativePanelDisposition::Fatal,
            static_cast<int64_t>(result), false, 0u, 0u);
    }
}

class CodecTimingInvocation : public NativePanelInvocation
{
public:
    CodecTimingInvocation(
        const std::string& identity,
        FrozenTimingScope scope,
        const NativeArmSpec& spec,
        uint32_t K,
        uint32_t block_bytes,
        const std::vector<uint8_t>& source,
        const std::vector<uint32_t>& packet_ids,
        uint32_t fixed_received_overhead,
        uint32_t receive_overhead_cap)
        : IdentityValue(identity)
        , Scope(scope)
        , FixedReceivedOverhead(fixed_received_overhead)
        , ReceiveOverheadCap(receive_overhead_cap)
        , PreparationResult(Wirehair_Error)
    {
        if (scope == FrozenTimingScope::EncoderInitPlusFirstKSymbols)
        {
            Encoder.reset(new NativeEncoderFixture);
            PreparationResult = Encoder->Initialize(
                spec, K, block_bytes, source);
            return;
        }

        NativeArm arm;
        PreparationResult = arm.Initialize(spec, K, block_bytes, source);
        if (PreparationResult != Wirehair_Success) {
            return;
        }
        if (scope == FrozenTimingScope::ReceiveToSuccess)
        {
            Receive.reset(new NativeReceiveFixture);
            PreparationResult = Receive->Initialize(
                arm, packet_ids, receive_overhead_cap);
            return;
        }
        if (scope == FrozenTimingScope::DecoderSolve)
        {
            Solve.reset(new NativeSolveFixture);
            PreparationResult = Solve->Initialize(
                arm, packet_ids, fixed_received_overhead);
            return;
        }
        PreparationResult = Wirehair_InvalidInput;
    }

    std::string Identity() const override
    {
        return IdentityValue;
    }

    NativePanelInvocationResult Invoke() override
    {
        if (PreparationResult != Wirehair_Success) {
            return ClassifyTimedResult(
                PreparationResult, 0u, false, 0u, false);
        }
        if (Scope == FrozenTimingScope::EncoderInitPlusFirstKSymbols)
        {
            const wirehair_wh2_bench::TimedArmResult result = Encoder->Run();
            return ClassifyTimedResult(
                result.Result, result.ElapsedNanoseconds, false, 0u, false);
        }
        if (Scope == FrozenTimingScope::ReceiveToSuccess)
        {
            const wirehair_wh2_bench::TimedArmResult result = Receive->Run();
            const bool has_extra = result.Result == Wirehair_Success &&
                result.BytesVerified &&
                result.DecodedOverhead <= ReceiveOverheadCap;
            if (result.Result == Wirehair_Success && !has_extra) {
                return ClassifyTimedResult(
                    Wirehair_Error, 0u, false, 0u, false);
            }
            return ClassifyTimedResult(
                result.Result, result.ElapsedNanoseconds,
                has_extra, result.DecodedOverhead, true);
        }
        if (Scope == FrozenTimingScope::DecoderSolve)
        {
            const wirehair_wh2_bench::IsolatedSolveResult result =
                Solve->Run();
            if (result.Result == Wirehair_Success && !result.BytesVerified) {
                return ClassifyTimedResult(
                    Wirehair_Error, 0u, false, 0u, false);
            }
            return ClassifyTimedResult(
                result.Result, result.ElapsedNanoseconds,
                result.Result == Wirehair_Success,
                FixedReceivedOverhead, true);
        }
        return ClassifyTimedResult(
            Wirehair_InvalidInput, 0u, false, 0u, false);
    }

private:
    std::string IdentityValue;
    FrozenTimingScope Scope;
    uint32_t FixedReceivedOverhead;
    uint32_t ReceiveOverheadCap;
    WirehairResult PreparationResult;
    std::unique_ptr<NativeEncoderFixture> Encoder;
    std::unique_ptr<NativeReceiveFixture> Receive;
    std::unique_ptr<NativeSolveFixture> Solve;
};

bool ArmIndexForName(
    const std::vector<ArmDescriptor>& arms,
    const std::string& name,
    std::size_t& index)
{
    for (std::size_t i = 0u; i < arms.size(); ++i)
    {
        if (name == arms[i].Arm) {
            index = i;
            return true;
        }
    }
    return false;
}

bool TimingOutcomeJson(
    const NativePanelInvocationResult& invocation,
    const char*& outcome,
    bool& has_decoded_extra,
    uint32_t& decoded_extra)
{
    has_decoded_extra = invocation.HasDecodedExtra;
    decoded_extra = invocation.DecodedExtra;
    if (invocation.Disposition == NativePanelDisposition::Success &&
        invocation.OutcomeCode == static_cast<int64_t>(Wirehair_Success))
    {
        outcome = "success";
        return true;
    }
    if (invocation.Disposition != NativePanelDisposition::PreflightFailure ||
        invocation.HasDecodedExtra)
    {
        return false;
    }
    switch (static_cast<WirehairResult>(invocation.OutcomeCode))
    {
    case Wirehair_NeedMore:
        outcome = "need_more_at_cap";
        return true;
    case Wirehair_BadDenseSeed:
    case Wirehair_BadPeelSeed:
        outcome = "construct_failed";
        return true;
    case Wirehair_BadInput_SmallN:
    case Wirehair_BadInput_LargeN:
    case Wirehair_UnsupportedPlatform:
        outcome = "unsupported";
        return true;
    default:
        return false;
    }
}

std::string TimingInvocationIdentity(
    const ArmDescriptor& arm,
    const std::string& realized_sha256,
    FrozenTimingScope scope)
{
    std::string identity = arm.DescriptorSha256;
    identity += ':';
    identity += realized_sha256;
    identity += ':';
    identity += wirehair::wh2_benchmark::FrozenTimingScopeName(scope);
    return identity;
}

bool TimingQualificationOutcomeJson(
    WirehairResult result,
    bool wirehair1,
    uint32_t decoded_overhead,
    const char*& outcome,
    bool& has_decoded_extra,
    uint32_t& decoded_extra)
{
    has_decoded_extra = false;
    decoded_extra = 0u;
    if (result == Wirehair_NeedMore)
    {
        if (wirehair1 && decoded_overhead != UINT32_MAX) {
            return false;
        }
        outcome = "need_more_at_bound";
        return true;
    }
    if (result != Wirehair_Success) {
        return false;
    }
    if (wirehair1)
    {
        if (decoded_overhead > kTimingReceiveOverheadCap) {
            return false;
        }
        decoded_extra = decoded_overhead;
    }
    else {
        decoded_extra = kRecoveryOverheadCap;
    }
    has_decoded_extra = true;
    outcome = "success";
    return true;
}

std::string TimingQualificationPayload(
    const FrozenTimingCell& cell,
    const FrozenPacketTrace& trace,
    const std::string& cell_sha256,
    const char* wirehair1_outcome,
    bool wirehair1_has_extra,
    uint32_t wirehair1_extra,
    const char* wirehair2_outcome,
    bool wirehair2_has_extra,
    uint32_t wirehair2_extra)
{
    std::string json;
    json.reserve(1024u);
    json += "{\"base_cell_sha256\":\"";
    json += cell.base_cell_sha256;
    json += "\",\"candidate_count\":";
    json += std::to_string(trace.attempted_candidates);
    json += ",\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"loss_retry_offset\":";
    json += std::to_string(cell.loss_retry_offset);
    json += ",\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"ordinal\":";
    json += std::to_string(cell.ordinal);
    json += ",\"packet_count\":";
    json += std::to_string(trace.delivered_ids.size());
    json += ",\"trace_sha256\":\"";
    json += trace.trace_sha256;
    json += "\",\"wirehair1_decoded_extra\":";
    if (wirehair1_has_extra) json += std::to_string(wirehair1_extra);
    else json += "null";
    json += ",\"wirehair1_outcome\":\"";
    json += wirehair1_outcome;
    json += "\",\"wirehair2_head_decoded_extra\":";
    if (wirehair2_has_extra) json += std::to_string(wirehair2_extra);
    else json += "null";
    json += ",\"wirehair2_head_outcome\":\"";
    json += wirehair2_outcome;
    json += "\"}";
    return json;
}

bool RunTimingQualificationJob(
    const WorkerContext& context,
    TimingQualificationProbeCache& cache,
    std::size_t cell_ordinal,
    std::size_t retry_offset,
    std::string& envelope,
    std::string& error)
{
    static const std::vector<FrozenTimingBaseCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentTimingBaseCells();
    if (context.Candidate || cells.size() != 2304u ||
        cell_ordinal >= cells.size() || retry_offset >= 256u)
    {
        error = "timing qualification index is outside the frozen domain";
        return false;
    }
    const FrozenTimingBaseCell& base = cells[cell_ordinal];
    FrozenTimingCell cell;
    if (base.ordinal != cell_ordinal ||
        !wirehair::wh2_benchmark::QualifyDevelopmentTimingCell(
            base, static_cast<uint32_t>(retry_offset), cell) ||
        cell.ordinal != cell_ordinal ||
        cell.fixed_received_overhead != kRecoveryOverheadCap ||
        cell.receive_overhead_cap != kTimingReceiveOverheadCap)
    {
        error = "timing qualification cell differs from the frozen domain";
        return false;
    }

    int actual_cpu = -1;
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t started_ns = 0u;
    if (!MonotonicNanoseconds(started_ns)) {
        error = "cannot read CLOCK_MONOTONIC before timing qualification";
        return false;
    }

    FrozenPacketTrace trace;
    wirehair::wh2_benchmark::FrozenTimingTraceReceipt trace_receipt;
    if (wirehair::wh2_benchmark::GenerateDevelopmentTimingTrace(
            cell, trace, trace_receipt) != FrozenTraceStatus::Complete ||
        trace.delivered_ids.size() !=
            static_cast<std::size_t>(cell.K) + cell.receive_overhead_cap)
    {
        error = "cannot generate complete timing qualification trace";
        return false;
    }
    const std::string base_json =
        wirehair::wh2_benchmark::CanonicalTimingBaseCellJson(base);
    const std::string base_sha256 =
        wirehair::wh2_benchmark::TimingBaseCellSha256(base);
    const std::string cell_sha256 =
        wirehair::wh2_benchmark::TimingCellSha256(cell);
    if (base_json.empty() || !IsLowerSha256(base_sha256) ||
        !IsLowerSha256(cell_sha256) || !IsLowerSha256(trace.trace_sha256) ||
        cell.base_cell_sha256 != base_sha256 ||
        wirehair::wh2_benchmark::CanonicalTimingSourceJson(cell) != base_json ||
        wirehair::wh2_benchmark::TimingSourceIdentitySha256(cell) !=
            base_sha256 ||
        trace_receipt.cell_sha256 != cell_sha256 ||
        trace_receipt.trace_sha256 != trace.trace_sha256)
    {
        error = "timing qualification identity is inconsistent";
        return false;
    }

    std::size_t head_index = 0u;
    NativeArmSpec head_spec;
    const uint32_t construction_attempt =
        ConstructionAttemptFor(0u, cell.base_seed_attempt);
    if (!ArmIndexForName(context.Arms, "wirehair2_head", head_index) ||
        head_index != 0u ||
        !ArmSpecFor(head_index, construction_attempt, nullptr, head_spec))
    {
        error = "cannot construct timing qualification control";
        return false;
    }
    if (!cache.Initialized || cache.K != cell.K ||
        cache.BlockBytes != cell.block_bytes ||
        cache.ConstructionAttempt != construction_attempt)
    {
        NativeTimingControlProbe next;
        const WirehairResult init = next.Initialize(
            head_spec, cell.K, cell.block_bytes);
        if (init != Wirehair_Success || !next.IsInitialized())
        {
            error = "fatal timing qualification control initialization";
            return false;
        }
        cache.Probe = std::move(next);
        cache.K = cell.K;
        cache.BlockBytes = cell.block_bytes;
        cache.ConstructionAttempt = construction_attempt;
        cache.Initialized = true;
    }

    const NativeTimingControlQualificationResult result =
        cache.Probe.Run(trace.delivered_ids);
    if (result.Qualification == NativeTimingControlQualification::Fatal)
    {
        error = "fatal timing qualification control result";
        return false;
    }
    const char* wirehair1_outcome = nullptr;
    const char* wirehair2_outcome = nullptr;
    bool wirehair1_has_extra = false;
    bool wirehair2_has_extra = false;
    uint32_t wirehair1_extra = 0u;
    uint32_t wirehair2_extra = 0u;
    if (!TimingQualificationOutcomeJson(
            result.Wirehair1Result, true,
            result.Wirehair1DecodedOverhead,
            wirehair1_outcome, wirehair1_has_extra, wirehair1_extra) ||
        !TimingQualificationOutcomeJson(
            result.Wirehair2HeadResult, false, UINT32_MAX,
            wirehair2_outcome, wirehair2_has_extra, wirehair2_extra))
    {
        error = "timing qualification result is internally inconsistent";
        return false;
    }
    const bool both_success =
        result.Wirehair1Result == Wirehair_Success &&
        result.Wirehair2HeadResult == Wirehair_Success;
    if ((result.Qualification == NativeTimingControlQualification::Success) !=
            both_success ||
        (result.Qualification == NativeTimingControlQualification::NeedMore) ==
            both_success)
    {
        error = "timing qualification disposition disagrees with controls";
        return false;
    }
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t finished_ns = 0u;
    if (!MonotonicNanoseconds(finished_ns) || finished_ns < started_ns) {
        error = "cannot read CLOCK_MONOTONIC after timing qualification";
        return false;
    }

    const std::size_t ordinal = cell_ordinal * 256u + retry_offset;
    const std::string message_sha256 =
        TimingQualificationMessageSha256(base_json);
    const std::string work_sha256 =
        TimingQualificationWorkSha256(cell_sha256, ordinal);
    const std::string payload = TimingQualificationPayload(
        cell, trace, cell_sha256,
        wirehair1_outcome, wirehair1_has_extra, wirehair1_extra,
        wirehair2_outcome, wirehair2_has_extra, wirehair2_extra);
    if (!IsLowerSha256(message_sha256) || !IsLowerSha256(work_sha256) ||
        payload.empty())
    {
        error = "cannot bind timing qualification evidence";
        return false;
    }
    envelope = Envelope(
        kTimingQualificationSchema, ordinal, actual_cpu, context.Pid,
        started_ns, finished_ns, context.ProcessStartTicks,
        context.BinarySha256, message_sha256, work_sha256, payload);
    return true;
}

std::string TimingPayload(
    const FrozenTimingCell& cell,
    const FrozenTimingPanel& panel,
    const WorkerContext& context,
    std::size_t left_index,
    std::size_t right_index,
    uint32_t left_attempt,
    uint32_t right_attempt,
    const std::string& left_realized,
    const std::string& right_realized,
    const char* order,
    const char* left_outcome,
    const char* right_outcome,
    bool left_has_extra,
    uint32_t left_extra,
    bool right_has_extra,
    uint32_t right_extra,
    const NativePanelResult& result,
    const std::string& cell_sha256,
    const std::string& trace_sha256)
{
    const ArmDescriptor& left = context.Arms[left_index];
    const ArmDescriptor& right = context.Arms[right_index];
    std::string json;
    json.reserve(1900u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"band\":\"";
    json += cell.band;
    json += "\",\"base_cell_sha256\":\"";
    json += cell.base_cell_sha256;
    json += "\",\"base_loss_seed\":\"";
    json += HexSeed(cell.base_loss_seed);
    json += "\",\"base_seed_attempt\":";
    json += std::to_string(cell.base_seed_attempt);
    json += ",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json += ",\"cell_sha256\":\"";
    json += cell_sha256;
    json += "\",\"elapsed_ns\":[";
    for (std::size_t i = 0u; i < result.Slots.size(); ++i)
    {
        if (i != 0u) json += ',';
        if (result.Slots[i].HasElapsedNanoseconds) {
            json += std::to_string(result.Slots[i].ElapsedNanoseconds);
        }
        else {
            json += "null";
        }
    }
    json += "],\"fixed_received_overhead\":";
    json += std::to_string(cell.fixed_received_overhead);
    json += ",\"interleave_policy\":\"";
    json += cell.interleave_policy;
    json += "\"";
    json += ",\"invocations_per_slot\":";
    json += std::to_string(cell.invocations_per_slot);
    json += ",\"left_arm\":\"";
    json += left.Arm;
    json += "\",\"left_arm_descriptor_sha256\":\"";
    json += left.DescriptorSha256;
    json += "\",\"left_binary_sha256\":\"";
    json += context.BinarySha256;
    json += "\",\"left_construction_attempt\":";
    json += std::to_string(left_attempt);
    json += ",\"left_decoded_extra\":";
    if (left_has_extra) json += std::to_string(left_extra);
    else json += "null";
    json += ",\"left_outcome\":\"";
    json += left_outcome;
    json += "\",\"left_realized_construction_sha256\":\"";
    json += left_realized;
    json += "\",\"left_repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"loss_retry_offset\":";
    json += std::to_string(cell.loss_retry_offset);
    json += ",\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"order\":\"";
    json += order;
    json += "\",\"panel_kind\":\"";
    json += wirehair::wh2_benchmark::FrozenPanelKindName(panel.panel_kind);
    json += "\",\"phase\":\"";
    json += cell.phase;
    json += "\",\"receive_overhead_cap\":";
    json += std::to_string(cell.receive_overhead_cap);
    json += ",\"replicate\":";
    json += std::to_string(cell.replicate);
    json += ",\"right_arm\":\"";
    json += right.Arm;
    json += "\",\"right_arm_descriptor_sha256\":\"";
    json += right.DescriptorSha256;
    json += "\",\"right_binary_sha256\":\"";
    json += context.BinarySha256;
    json += "\",\"right_construction_attempt\":";
    json += std::to_string(right_attempt);
    json += ",\"right_decoded_extra\":";
    if (right_has_extra) json += std::to_string(right_extra);
    else json += "null";
    json += ",\"right_outcome\":\"";
    json += right_outcome;
    json += "\",\"right_realized_construction_sha256\":\"";
    json += right_realized;
    json += "\",\"right_repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"schedule\":\"";
    json += FrozenScheduleName(cell.schedule);
    json += "\",\"scope\":\"";
    json += wirehair::wh2_benchmark::FrozenTimingScopeName(panel.scope);
    json += "\",\"trace_sha256\":\"";
    json += trace_sha256;
    json += "\"}";
    return json;
}

const char* PhaseOutcomeName(WirehairResult result)
{
    switch (result)
    {
    case Wirehair_Success: return "success";
    case Wirehair_NeedMore: return "need_more_at_cap";
    case Wirehair_BadDenseSeed:
    case Wirehair_BadPeelSeed: return "construct_failed";
    default: return nullptr;
    }
}

std::string PhaseCountersJson(
    const wirehair_v2::PrecodeSolveStats& stats)
{
    std::string json;
    json.reserve(720u);
    json += "{\"binary_adjacency_storage_allocations\":";
    json += std::to_string(stats.BinaryAdjacencyStorageAllocations);
    json += ",\"binary_adjacency_storage_bytes\":";
    json += std::to_string(stats.BinaryAdjacencyStorageBytes);
    json += ",\"binary_residual_rank\":";
    json += std::to_string(stats.BinaryResidualRank);
    json += ",\"binary_row_references\":";
    json += std::to_string(stats.BinaryRowReferences);
    json += ",\"binary_row_storage_allocations\":";
    json += std::to_string(stats.BinaryRowStorageAllocations);
    json += ",\"binary_row_storage_bytes\":";
    json += std::to_string(stats.BinaryRowStorageBytes);
    json += ",\"block_xors\":";
    json += std::to_string(stats.BlockXors);
    json +=
        ",\"full_block_gf256_multiply_add_divide_normalize_ops\":";
    json += std::to_string(stats.BlockMulAdds);
    json += ",\"inactivated_columns\":";
    json += std::to_string(stats.InactivatedColumns);
    json += ",\"packet_rows\":";
    json += std::to_string(stats.PacketRows);
    json += ",\"packet_seed_attempt\":";
    json += std::to_string(stats.PacketSeedAttempt);
    json += ",\"peeled_columns\":";
    json += std::to_string(stats.PeeledColumns);
    json += ",\"residual_rank\":";
    json += std::to_string(stats.ResidualRank);
    json += ",\"residual_rows\":";
    json += std::to_string(stats.ResidualRows);
    json += "}";
    return json;
}

std::string PhaseCounterSha256(
    const wirehair_v2::PrecodeSolveStats& stats)
{
    return wirehair::wh2_benchmark::Sha256Hex(PhaseCountersJson(stats));
}

bool AppendPhaseObservationJson(
    const PhaseSolveObservation& observation,
    bool timing_visible,
    std::string& json)
{
    const char* const outcome = PhaseOutcomeName(observation.Result);
    const std::string counter_sha256 = PhaseCounterSha256(observation.Stats);
    if ((observation.Arm != PhaseSolveArm::Two07 &&
         observation.Arm != PhaseSolveArm::Head) ||
        !outcome || !IsLowerSha256(counter_sha256))
    {
        return false;
    }
    json += "{\"arm\":\"";
    json += observation.Arm == PhaseSolveArm::Two07 ?
        kTimingCandidateArm : "wirehair2_head";
    json += "\",\"back_sub_ns\":";
    json += timing_visible ?
        std::to_string(observation.Stats.BackSubNanoseconds) : "null";
    json += ",\"build_ns\":";
    json += timing_visible ?
        std::to_string(observation.Stats.BuildNanoseconds) : "null";
    json += ",\"bytes_verified\":";
    json += observation.BytesVerified ? "true" : "false";
    json += ",\"counter_sha256\":\"";
    json += counter_sha256;
    json += "\",\"outcome\":\"";
    json += outcome;
    json += "\",\"outer_ns\":";
    json += timing_visible ?
        std::to_string(observation.ElapsedNanoseconds) : "null";
    json += ",\"peel_ns\":";
    json += timing_visible ?
        std::to_string(observation.Stats.PeelNanoseconds) : "null";
    json += ",\"project_ns\":";
    json += timing_visible ?
        std::to_string(observation.Stats.ProjectNanoseconds) : "null";
    json += ",\"residual_ns\":";
    json += timing_visible ?
        std::to_string(observation.Stats.ResidualNanoseconds) : "null";
    json += ",\"wirehair_result\":";
    json += std::to_string(static_cast<int>(observation.Result));
    json += "}";
    return true;
}

void AppendNullablePhaseTotal(
    bool present,
    uint64_t value,
    std::string& json)
{
    if (present) json += std::to_string(value);
    else json += "null";
}

std::string PhaseAttributionPayload(
    const PhaseCoordinate& coordinate,
    std::size_t profile_ordinal,
    const WorkerContext& context,
    const std::string& phase_cell_sha256,
    const std::string& trace_sha256,
    const std::string& left_realized,
    const std::string& right_realized,
    const PhasePanelAssembly& assembly)
{
    const char* profile = nullptr;
    uint32_t invocations_per_slot = 0u;
    const FrozenTimingCell& cell = coordinate.Qualified;
    const std::string qualified_sha256 =
        wirehair::wh2_benchmark::TimingCellSha256(cell);
    if (!PhaseProfile(
            profile_ordinal, profile, invocations_per_slot) ||
        !ExactPhaseArmRoster(context.Arms) ||
        !IsLowerSha256(phase_cell_sha256) ||
        !IsLowerSha256(trace_sha256) ||
        !IsLowerSha256(left_realized) ||
        !IsLowerSha256(right_realized) ||
        assembly.Measured.size() !=
            static_cast<std::size_t>(invocations_per_slot) * 4u)
    {
        return std::string();
    }

    const std::string left_counters =
        PhaseCountersJson(assembly.LeftCounters);
    const std::string right_counters =
        PhaseCountersJson(assembly.RightCounters);
    std::string json;
    json.reserve(assembly.Measured.size() * 600u + 6000u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"band\":\"";
    json += cell.band;
    json += "\",\"base_loss_seed\":\"";
    json += HexSeed(cell.base_loss_seed);
    json += "\",\"binary_sha256\":\"";
    json += context.BinarySha256;
    json += "\",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json +=
        ",\"block_muladds_semantics\":\"full-block-gf256-multiply-add-divide-normalize-operations\"";
    json += ",\"cell_sha256\":\"";
    json += phase_cell_sha256;
    json += "\",\"construction_attempt\":0,\"coordinate_ordinal\":";
    json += std::to_string(coordinate.Ordinal);
    json += ",\"fixed_received_overhead\":4";
    json +=
        ",\"interleave_policy\":\"self-counterbalanced-repeat-major-v1\"";
    json += ",\"invocations_per_slot\":";
    json += std::to_string(invocations_per_slot);
    json += ",\"left_arm\":\"";
    json += context.Arms[2].Arm;
    json += "\",\"left_arm_descriptor_sha256\":\"";
    json += context.Arms[2].DescriptorSha256;
    json += "\",\"left_non_timing_counters\":";
    json += left_counters;
    json += ",\"left_realized_construction_sha256\":\"";
    json += left_realized;
    json += "\",\"left_repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"loss_retry_offset\":0,\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"measured_observations\":[";
    for (std::size_t i = 0u; i < assembly.Measured.size(); ++i)
    {
        const PhaseMeasuredObservation& measured = assembly.Measured[i];
        if (i != 0u) json += ',';
        json += "{\"block\":";
        json += std::to_string(measured.Block);
        json += ",\"observation\":";
        if (!AppendPhaseObservationJson(
                measured.Observation, assembly.Comparable, json))
        {
            return std::string();
        }
        json += ",\"repeat\":";
        json += std::to_string(measured.Repeat);
        json += ",\"slot\":";
        json += std::to_string(measured.Slot);
        json += "}";
    }
    json += "],\"order\":\"ABBA\",\"panel_comparable\":";
    json += assembly.Comparable ? "true" : "false";
    json += ",\"panel_kind\":\"ab\",\"profile\":\"";
    json += profile;
    json += "\",\"profile_ordinal\":";
    json += std::to_string(profile_ordinal);
    json += ",\"receive_overhead_cap\":256,\"replicate\":0";
    json += ",\"right_arm\":\"";
    json += context.Arms[0].Arm;
    json += "\",\"right_arm_descriptor_sha256\":\"";
    json += context.Arms[0].DescriptorSha256;
    json += "\",\"right_non_timing_counters\":";
    json += right_counters;
    json += ",\"right_realized_construction_sha256\":\"";
    json += right_realized;
    json += "\",\"right_repair_map_sha256\":\"";
    json += kZeroSha256;
    json += "\",\"schedule\":\"iid\",\"scope\":\"decoder_solve\"";
    json += ",\"slot_sums\":[";
    for (std::size_t slot = 0u; slot < assembly.Slots.size(); ++slot)
    {
        if (slot != 0u) json += ',';
        const PhaseSlotTotals& totals = assembly.Slots[slot];
        json += "{\"back_sub_ns\":";
        AppendNullablePhaseTotal(
            totals.HasElapsed, totals.BackSubNanoseconds, json);
        json += ",\"build_ns\":";
        AppendNullablePhaseTotal(
            totals.HasElapsed, totals.BuildNanoseconds, json);
        json += ",\"outer_ns\":";
        AppendNullablePhaseTotal(
            totals.HasElapsed, totals.OuterNanoseconds, json);
        json += ",\"peel_ns\":";
        AppendNullablePhaseTotal(
            totals.HasElapsed, totals.PeelNanoseconds, json);
        json += ",\"project_ns\":";
        AppendNullablePhaseTotal(
            totals.HasElapsed, totals.ProjectNanoseconds, json);
        json += ",\"residual_ns\":";
        AppendNullablePhaseTotal(
            totals.HasElapsed, totals.ResidualNanoseconds, json);
        json += ",\"slot\":";
        json += std::to_string(slot);
        json += "}";
    }
    json += "],\"source_base_cell_sha256\":\"";
    json += cell.base_cell_sha256;
    json += "\",\"trace_qualified_timing_cell_sha256\":\"";
    json += qualified_sha256;
    json += "\",\"trace_sha256\":\"";
    json += trace_sha256;
    json += "\",\"warmups\":{\"left\":";
    if (!AppendPhaseObservationJson(
            assembly.LeftWarmup, assembly.Comparable, json))
    {
        return std::string();
    }
    json += ",\"right\":";
    if (!AppendPhaseObservationJson(
            assembly.RightWarmup, assembly.Comparable, json))
    {
        return std::string();
    }
    json += "},\"weak_ledger\":";
    json += assembly.Comparable ? "false" : "true";
    json += "}";
    return json;
}

WirehairResult PreparePhaseSolveFixture(
    const NativeArmSpec& spec,
    uint32_t K,
    uint32_t block_bytes,
    const std::vector<uint8_t>& source,
    const std::vector<uint32_t>& packet_ids,
    std::shared_ptr<const NativeSolveFixture>& fixture_out)
{
    fixture_out.reset();
    NativeArm arm;
    const WirehairResult arm_result =
        arm.Initialize(spec, K, block_bytes, source);
    if (arm_result != Wirehair_Success) {
        return arm_result;
    }
    try
    {
        const std::shared_ptr<NativeSolveFixture> fixture(
            new NativeSolveFixture);
        const WirehairResult fixture_result = fixture->Initialize(
            arm, packet_ids, kRecoveryOverheadCap);
        if (fixture_result == Wirehair_Success) {
            fixture_out = fixture;
        }
        return fixture_result;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool RunPhaseAttributionJobImpl(
    const WorkerContext& context,
    std::size_t coordinate_ordinal,
    std::size_t profile_ordinal,
    std::string& envelope,
    std::string& error)
{
    std::vector<PhaseCoordinate> coordinates;
    const char* profile = nullptr;
    uint32_t invocations_per_slot = 0u;
    if (context.Candidate ||
        !ExactPhaseArmRoster(context.Arms) ||
        !BuildPhaseCoordinates(coordinates) ||
        coordinates.size() != kPhaseCoordinateCount ||
        coordinate_ordinal >= coordinates.size() ||
        !PhaseProfile(
            profile_ordinal, profile, invocations_per_slot) ||
        !profile)
    {
        error = "phase-attribution coordinate/profile is outside the frozen domain";
        return false;
    }
    const PhaseCoordinate& coordinate = coordinates[coordinate_ordinal];
    const FrozenTimingCell& cell = coordinate.Qualified;
    if (coordinate.Ordinal != coordinate_ordinal ||
        cell.base_seed_attempt != 0u ||
        cell.fixed_received_overhead != kRecoveryOverheadCap ||
        cell.receive_overhead_cap != kTimingReceiveOverheadCap ||
        invocations_per_slot != (profile_ordinal == 0u ? 16u : 24u))
    {
        error = "phase-attribution protocol/profile closure failed";
        return false;
    }

    int actual_cpu = -1;
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t started_ns = 0u;
    if (!MonotonicNanoseconds(started_ns)) {
        error = "cannot read CLOCK_MONOTONIC before phase attribution";
        return false;
    }

    FrozenPacketTrace trace;
    wirehair::wh2_benchmark::FrozenTimingTraceReceipt receipt;
    if (!GeneratePhaseTrace(coordinate, trace, receipt))
    {
        error = "cannot regenerate complete phase-attribution trace";
        return false;
    }
    const std::string source_json =
        wirehair::wh2_benchmark::CanonicalTimingSourceJson(cell);
    const std::string qualified_sha256 =
        wirehair::wh2_benchmark::TimingCellSha256(cell);
    uint64_t source_seed = 0u;
    std::vector<uint8_t> source;
    if (source_json.empty() ||
        source_json !=
            wirehair::wh2_benchmark::CanonicalTimingBaseCellJson(
                coordinate.Base) ||
        !IsLowerSha256(qualified_sha256) ||
        receipt.cell_sha256 != qualified_sha256 ||
        receipt.trace_sha256 != trace.trace_sha256 ||
        !SourceSeedFromCellJson(source_json, source_seed) ||
        !wirehair_wh2_bench::MakeDeterministicSource(
            cell.K, cell.block_bytes, source_seed, source))
    {
        error = "cannot construct phase-attribution source/provenance";
        return false;
    }
    const std::string message_sha256 = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    if (!IsLowerSha256(message_sha256)) {
        error = "cannot hash phase-attribution source";
        return false;
    }

    NativeArmSpec left_spec;
    NativeArmSpec right_spec;
    std::string left_realized;
    std::string right_realized;
    if (!ArmSpecFor(2u, 0u, nullptr, left_spec) ||
        !ArmSpecFor(0u, 0u, nullptr, right_spec) ||
        !ValidateResolvedProductionArm(
            2u, left_spec, cell.K, cell.block_bytes, 0u) ||
        !ValidateResolvedProductionArm(
            0u, right_spec, cell.K, cell.block_bytes, 0u) ||
        !RealizedConstructionSha256(
            context.Arms[2], cell.K, cell.block_bytes,
            0u, left_realized) ||
        !RealizedConstructionSha256(
            context.Arms[0], cell.K, cell.block_bytes,
            0u, right_realized))
    {
        error = "phase-attribution arms differ from production descriptors";
        return false;
    }
    const std::string phase_cell_sha256 = PhaseCellSha256(
        coordinate, profile_ordinal, context.Arms, trace.trace_sha256,
        left_realized, right_realized);
    if (!IsLowerSha256(phase_cell_sha256)) {
        error = "cannot hash phase-attribution cell";
        return false;
    }

    const std::shared_ptr<PhaseObservationCollector> collector(
        new PhaseObservationCollector);
    std::shared_ptr<const NativeSolveFixture> left_fixture;
    std::shared_ptr<const NativeSolveFixture> right_fixture;
    const WirehairResult left_preparation = PreparePhaseSolveFixture(
        left_spec, cell.K, cell.block_bytes, source,
        trace.delivered_ids, left_fixture);
    const WirehairResult right_preparation = PreparePhaseSolveFixture(
        right_spec, cell.K, cell.block_bytes, source,
        trace.delivered_ids, right_fixture);
    const auto valid_preparation = [](
        WirehairResult result,
        const std::shared_ptr<const NativeSolveFixture>& fixture) -> bool
    {
        const bool weak = result == Wirehair_NeedMore ||
            result == Wirehair_BadDenseSeed ||
            result == Wirehair_BadPeelSeed;
        return (result == Wirehair_Success && fixture &&
                fixture->IsInitialized()) ||
            (weak && !fixture);
    };
    if (!valid_preparation(left_preparation, left_fixture) ||
        !valid_preparation(right_preparation, right_fixture))
    {
        error = "fatal phase-attribution fixture preparation";
        return false;
    }
    const std::string left_identity = TimingInvocationIdentity(
        context.Arms[2], left_realized, FrozenTimingScope::DecoderSolve) +
        ":phase-attribution-v1";
    const std::string right_identity = TimingInvocationIdentity(
        context.Arms[0], right_realized, FrozenTimingScope::DecoderSolve) +
        ":phase-attribution-v1";
    const NativePanelArm left_arm(
        left_identity,
        [left_identity, left_preparation, left_fixture, collector]() {
            return std::unique_ptr<NativePanelInvocation>(
                new PhaseSolveInvocation(
                    left_identity, PhaseSolveArm::Two07,
                    left_preparation, left_fixture,
                    kRecoveryOverheadCap, collector));
        });
    const NativePanelArm right_arm(
        right_identity,
        [right_identity, right_preparation, right_fixture, collector]() {
            return std::unique_ptr<NativePanelInvocation>(
                new PhaseSolveInvocation(
                    right_identity, PhaseSolveArm::Head,
                    right_preparation, right_fixture,
                    kRecoveryOverheadCap, collector));
        });
    const NativePanelResult panel =
        wirehair_wh2_bench::ExecuteNativeTimingPanel(
            context.Cpu, NativePanelOrder::ABBA,
            invocations_per_slot, left_arm, right_arm);
    if (panel.Status != NativePanelStatus::Complete ||
        panel.TargetCpu != context.Cpu ||
        panel.Order != NativePanelOrder::ABBA ||
        panel.InvocationsPerSlot != invocations_per_slot ||
        !panel.HasLeftPreflight || !panel.HasRightPreflight)
    {
        error = std::string("fatal phase-attribution panel: ") +
            wirehair_wh2_bench::NativePanelStatusName(panel.Status) +
            (panel.Diagnostic.empty() ? "" : ": " + panel.Diagnostic);
        return false;
    }
    PhasePanelAssembly assembly;
    std::string assembly_error;
    if (!wirehair_wh2_bench::ValidateAndAssemblePhasePanel(
            panel, NativePanelOrder::ABBA, invocations_per_slot,
            collector->Observations(), assembly, assembly_error))
    {
        error = "phase-attribution reconstruction failed: " +
            assembly_error;
        return false;
    }
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t finished_ns = 0u;
    if (!MonotonicNanoseconds(finished_ns) || finished_ns < started_ns) {
        error = "cannot read CLOCK_MONOTONIC after phase attribution";
        return false;
    }

    const std::string payload = PhaseAttributionPayload(
        coordinate, profile_ordinal, context, phase_cell_sha256,
        trace.trace_sha256, left_realized, right_realized, assembly);
    const std::size_t ordinal =
        coordinate_ordinal * kPhaseProfileCount + profile_ordinal;
    const std::string work_sha256 =
        PhaseAttributionWorkSha256(phase_cell_sha256, ordinal);
    if (payload.empty() || !IsLowerSha256(work_sha256))
    {
        error = "cannot bind phase-attribution evidence";
        return false;
    }
    envelope = Envelope(
        kPhaseAttributionSchema, ordinal, actual_cpu, context.Pid,
        started_ns, finished_ns, context.ProcessStartTicks,
        context.BinarySha256, message_sha256, work_sha256, payload);
    return true;
}

bool RunPhaseAttributionJob(
    const WorkerContext& context,
    std::size_t coordinate_ordinal,
    std::size_t profile_ordinal,
    std::string& envelope,
    std::string& error)
{
    try {
        return RunPhaseAttributionJobImpl(
            context, coordinate_ordinal, profile_ordinal, envelope, error);
    }
    catch (const std::bad_alloc&) {
        error = "phase-attribution allocation failed";
    }
    catch (const std::length_error&) {
        error = "phase-attribution allocation length is invalid";
    }
    catch (...) {
        error = "phase-attribution internal failure";
    }
    envelope.clear();
    return false;
}

bool RunTimingJob(
    const WorkerContext& context,
    std::size_t cell_ordinal,
    std::size_t panel_index,
    std::size_t retry_offset,
    std::string& envelope,
    std::string& error)
{
    static const std::vector<FrozenTimingBaseCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentTimingBaseCells();
    static const std::vector<FrozenTimingPanel> panels =
        wirehair::wh2_benchmark::EnumerateOneCandidateTimingPanels(
            kTimingCandidateArm);
    if (cells.size() != 2304u || panels.size() != 11u ||
        cell_ordinal >= cells.size() || panel_index >= panels.size() ||
        retry_offset >= 256u)
    {
        error = "timing job index is outside the frozen domain";
        return false;
    }
    const FrozenTimingBaseCell& base = cells[cell_ordinal];
    FrozenTimingCell cell;
    const FrozenTimingPanel& panel = panels[panel_index];
    if (!wirehair::wh2_benchmark::QualifyDevelopmentTimingCell(
            base, static_cast<uint32_t>(retry_offset), cell) ||
        cell.ordinal != cell_ordinal || panel.ordinal != panel_index ||
        cell.fixed_received_overhead != kRecoveryOverheadCap ||
        cell.receive_overhead_cap != kTimingReceiveOverheadCap ||
        cell.interleave_policy !=
            "self-counterbalanced-repeat-major-v1" ||
        cell.invocations_per_slot < 2u ||
        cell.invocations_per_slot !=
            wirehair::wh2_benchmark::
                DevelopmentTimingInvocationsPerSlot(cell.K))
    {
        error = "timing cell or panel differs from the frozen domain";
        return false;
    }
    int actual_cpu = -1;
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t started_ns = 0u;
    if (!MonotonicNanoseconds(started_ns)) {
        error = "cannot read CLOCK_MONOTONIC before timing job";
        return false;
    }

    FrozenPacketTrace trace;
    wirehair::wh2_benchmark::FrozenTimingTraceReceipt trace_receipt;
    if (wirehair::wh2_benchmark::GenerateDevelopmentTimingTrace(
            cell, trace, trace_receipt) != FrozenTraceStatus::Complete ||
        trace.delivered_ids.size() !=
            static_cast<std::size_t>(cell.K) + cell.receive_overhead_cap)
    {
        error = "cannot regenerate complete frozen timing trace";
        return false;
    }
    const std::string source_json =
        wirehair::wh2_benchmark::CanonicalTimingSourceJson(cell);
    const std::string source_identity_sha256 =
        wirehair::wh2_benchmark::TimingSourceIdentitySha256(cell);
    const std::string cell_sha256 =
        wirehair::wh2_benchmark::TimingCellSha256(cell);
    uint64_t source_seed = 0u;
    std::vector<uint8_t> source;
    if (source_json.empty() || source_json !=
            wirehair::wh2_benchmark::CanonicalTimingBaseCellJson(base) ||
        source_identity_sha256 != cell.base_cell_sha256 ||
        trace_receipt.cell_sha256 != cell_sha256 ||
        trace_receipt.trace_sha256 != trace.trace_sha256 ||
        !IsLowerSha256(cell_sha256) || !IsLowerSha256(trace.trace_sha256) ||
        !IsLowerSha256(cell.base_cell_sha256) ||
        !SourceSeedFromCellJson(source_json, source_seed) ||
        !wirehair_wh2_bench::MakeDeterministicSource(
            cell.K, cell.block_bytes, source_seed, source))
    {
        error = "cannot construct deterministic timing source";
        return false;
    }
    const std::string message_sha256 = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    if (!IsLowerSha256(message_sha256)) {
        error = "cannot hash deterministic timing source";
        return false;
    }

    std::size_t left_index = 0u;
    std::size_t right_index = 0u;
    if (!ArmIndexForName(context.Arms, panel.left_arm, left_index) ||
        !ArmIndexForName(context.Arms, panel.right_arm, right_index))
    {
        error = "timing panel names an unknown arm";
        return false;
    }
    const uint32_t left_attempt = ConstructionAttemptFor(
        left_index, cell.base_seed_attempt);
    const uint32_t right_attempt = ConstructionAttemptFor(
        right_index, cell.base_seed_attempt);
    NativeArmSpec left_spec;
    NativeArmSpec right_spec;
    if (!ArmSpecFor(
            left_index, left_attempt, context.Candidate, left_spec) ||
        !ArmSpecFor(
            right_index, right_attempt, context.Candidate, right_spec) ||
        !ValidateResolvedProductionArm(
            left_index, left_spec, cell.K, cell.block_bytes, left_attempt) ||
        !ValidateResolvedProductionArm(
            right_index, right_spec, cell.K, cell.block_bytes, right_attempt))
    {
        error = "timing arm differs from its production descriptor";
        return false;
    }
    std::string left_realized;
    std::string right_realized;
    if (!RealizedConstructionSha256(
            context.Arms[left_index], cell.K, cell.block_bytes,
            left_attempt, left_realized) ||
        !RealizedConstructionSha256(
            context.Arms[right_index], cell.K, cell.block_bytes,
            right_attempt, right_realized))
    {
        error = "cannot hash realized timing constructions";
        return false;
    }

    const std::string left_identity = TimingInvocationIdentity(
        context.Arms[left_index], left_realized, panel.scope);
    const std::string right_identity = TimingInvocationIdentity(
        context.Arms[right_index], right_realized, panel.scope);
    const NativePanelArm left_arm(
        left_identity,
        [left_identity, panel, left_spec, &cell, &source, &trace]() {
            return std::unique_ptr<NativePanelInvocation>(
                new CodecTimingInvocation(
                    left_identity, panel.scope, left_spec, cell.K,
                    cell.block_bytes, source, trace.delivered_ids,
                    cell.fixed_received_overhead,
                    cell.receive_overhead_cap));
        });
    const NativePanelArm right_arm(
        right_identity,
        [right_identity, panel, right_spec, &cell, &source, &trace]() {
            return std::unique_ptr<NativePanelInvocation>(
                new CodecTimingInvocation(
                    right_identity, panel.scope, right_spec, cell.K,
                    cell.block_bytes, source, trace.delivered_ids,
                    cell.fixed_received_overhead,
                    cell.receive_overhead_cap));
        });
    const FrozenTimingOrder frozen_order =
        wirehair::wh2_benchmark::TimingPanelOrder(panel, cell.replicate);
    NativePanelOrder native_order;
    if (frozen_order == FrozenTimingOrder::ABBA) {
        native_order = NativePanelOrder::ABBA;
    }
    else if (frozen_order == FrozenTimingOrder::BAAB) {
        native_order = NativePanelOrder::BAAB;
    }
    else {
        error = "timing panel has invalid frozen counterbalancing order";
        return false;
    }
    const NativePanelResult result =
        wirehair_wh2_bench::ExecuteNativeTimingPanel(
            context.Cpu, native_order, cell.invocations_per_slot,
            left_arm, right_arm);
    if (result.Status != NativePanelStatus::Complete ||
        !result.HasLeftPreflight || !result.HasRightPreflight)
    {
        error = std::string("fatal native timing panel: ") +
            wirehair_wh2_bench::NativePanelStatusName(result.Status) +
            (result.Diagnostic.empty() ? "" : ": " + result.Diagnostic);
        return false;
    }
    const char* left_outcome = nullptr;
    const char* right_outcome = nullptr;
    bool left_has_extra = false;
    bool right_has_extra = false;
    uint32_t left_extra = 0u;
    uint32_t right_extra = 0u;
    if (!TimingOutcomeJson(
            result.LeftPreflight, left_outcome,
            left_has_extra, left_extra) ||
        !TimingOutcomeJson(
            result.RightPreflight, right_outcome,
            right_has_extra, right_extra))
    {
        error = "fatal or internally inconsistent timing outcome";
        return false;
    }
    const bool left_success = std::string(left_outcome) == "success";
    const bool right_success = std::string(right_outcome) == "success";
    const auto valid_extra = [panel, &cell](
        bool success, bool has_extra, uint32_t extra) -> bool
    {
        if (!success) return !has_extra;
        if (panel.scope ==
            FrozenTimingScope::EncoderInitPlusFirstKSymbols)
        {
            return !has_extra;
        }
        if (panel.scope == FrozenTimingScope::DecoderSolve) {
            return has_extra && extra == cell.fixed_received_overhead;
        }
        return panel.scope == FrozenTimingScope::ReceiveToSuccess &&
            has_extra && extra <= cell.receive_overhead_cap;
    };
    if (!valid_extra(left_success, left_has_extra, left_extra) ||
        !valid_extra(right_success, right_has_extra, right_extra))
    {
        error = "timing decoded-extra differs from its frozen scope";
        return false;
    }
    const bool both_success =
        result.LeftPreflight.Disposition == NativePanelDisposition::Success &&
        result.RightPreflight.Disposition == NativePanelDisposition::Success;
    bool slots_valid = result.Order == native_order &&
        result.TargetCpu == context.Cpu &&
        result.InvocationsPerSlot == cell.invocations_per_slot &&
        result.Slots.size() == 8u;
    for (std::size_t slot = 0u; slot < result.Slots.size(); ++slot)
    {
        const std::size_t local_slot = slot % 4u;
        const bool primary_left = native_order == NativePanelOrder::ABBA ?
            (local_slot == 0u || local_slot == 3u) :
            (local_slot == 1u || local_slot == 2u);
        const bool left_slot = slot < 4u ? primary_left : !primary_left;
        const NativePanelSide expected_side = left_slot ?
            NativePanelSide::Left : NativePanelSide::Right;
        const bool has_elapsed = result.Slots[slot].HasElapsedNanoseconds;
        slots_valid = slots_valid &&
            result.Slots[slot].Side == expected_side &&
            has_elapsed == both_success;
        if (has_elapsed)
        {
            slots_valid = slots_valid &&
                result.Slots[slot].ElapsedNanoseconds != 0u &&
                result.Slots[slot].ElapsedNanoseconds <=
                    static_cast<uint64_t>(
                        std::numeric_limits<int64_t>::max()) &&
                result.Slots[slot].ElapsedNanoseconds ==
                    result.Slots[slot].Invocation.ElapsedNanoseconds;
        }
    }
    if (result.PanelComparable != both_success || !slots_valid ||
        (both_success && result.HasEightNullTimings()) ||
        (!both_success && !result.HasEightNullTimings()))
    {
        error = "timing panel comparability/timing slots are inconsistent";
        return false;
    }
    if (panel.panel_kind == FrozenPanelKind::AA &&
        (std::string(left_outcome) != right_outcome ||
         left_has_extra != right_has_extra ||
         (left_has_extra && left_extra != right_extra)))
    {
        error = "A/A timing sides produced different outcomes";
        return false;
    }
    if (!VerifySingletonCpu(context.Cpu, actual_cpu, error)) {
        return false;
    }
    uint64_t finished_ns = 0u;
    if (!MonotonicNanoseconds(finished_ns) || finished_ns < started_ns) {
        error = "cannot read CLOCK_MONOTONIC after timing job";
        return false;
    }
    const char* order = wirehair::wh2_benchmark::FrozenTimingOrderName(
        frozen_order);
    const std::string payload = TimingPayload(
        cell, panel, context, left_index, right_index,
        left_attempt, right_attempt, left_realized, right_realized,
        order, left_outcome, right_outcome,
        left_has_extra, left_extra, right_has_extra, right_extra,
        result, cell_sha256, trace.trace_sha256);
    const std::size_t ordinal = cell_ordinal * panels.size() + panel_index;
    const std::string work_sha256 = TimingWorkSha256(
        cell_sha256, ordinal);
    if (!IsLowerSha256(work_sha256)) {
        error = "cannot hash timing work identity";
        return false;
    }
    envelope = Envelope(
        kTimingSchema, ordinal, actual_cpu, context.Pid, started_ns,
        finished_ns, context.ProcessStartTicks, context.BinarySha256,
        message_sha256, work_sha256, payload);
    return true;
}

bool RunWorker(
    int cpu,
    const std::string& binary_sha256,
    const std::vector<ArmDescriptor>& arms,
    const CandidateDefinition* candidate)
{
    std::string error;
    if ((!candidate && !ValidateProductionTimingArmDefinitions()) ||
        !PinSingletonCpu(cpu, error))
    {
        std::cerr << "wirehair_wh2_contract_worker: " <<
            (error.empty() ?
                "production timing arm definition is invalid" : error) <<
            '\n';
        return false;
    }
#if defined(__linux__)
    const int pid = static_cast<int>(getpid());
#else
    const int pid = -1;
#endif
    if (pid <= 0) {
        std::cerr << "wirehair_wh2_contract_worker: invalid worker PID\n";
        return false;
    }
    WorkerContext context;
    context.Cpu = cpu;
    context.Pid = pid;
    if (!ReadSelfProcessStartTicks(context.ProcessStartTicks))
    {
        std::cerr << "wirehair_wh2_contract_worker: cannot read worker "
            "process start ticks\n";
        return false;
    }
    context.BinarySha256 = binary_sha256;
    context.Arms = arms;
    context.Candidate = candidate;
    TimingQualificationProbeCache qualification_cache;

    std::string line;
    while (std::getline(std::cin, line))
    {
        if (line == "Q") {
            return true;
        }
        char kind = 0;
        std::size_t cell_ordinal = 0u;
        std::size_t item_index = 0u;
        if (!ParseJobLine(line, kind, cell_ordinal, item_index))
        {
            std::cerr << "wirehair_wh2_contract_worker: malformed command\n";
            return false;
        }
        if (candidate && kind != 'R')
        {
            std::cerr << "wirehair_wh2_contract_worker: recovery candidate "
                "worker rejects qualification/timing jobs\n";
            return false;
        }
        std::string output;
        bool ok = false;
        if (kind == 'R') {
            ok = RunRecoveryJob(
                context, cell_ordinal, item_index, output, error);
        }
        else if (kind == 'L') {
            ok = RunTimingQualificationJob(
                context, qualification_cache,
                cell_ordinal, item_index, output, error);
        }
        else if (kind == 'T') {
            const std::size_t retry_offset = item_index & 255u;
            const std::size_t panel_index = item_index >> 8u;
            ok = RunTimingJob(
                context, cell_ordinal, panel_index, retry_offset,
                output, error);
        }
        else if (kind == 'P') {
            ok = RunPhaseAttributionJob(
                context, cell_ordinal, item_index, output, error);
        }
        else {
            error = "native command kind is unsupported";
            ok = false;
        }
        if (!ok)
        {
            std::cerr << "wirehair_wh2_contract_worker: " << error << '\n';
            return false;
        }
        if (!EmitLine(output))
        {
            std::cerr << "wirehair_wh2_contract_worker: output write failed\n";
            return false;
        }
    }
    if (std::cin.bad()) {
        std::cerr << "wirehair_wh2_contract_worker: input read failed\n";
    }
    else {
        std::cerr << "wirehair_wh2_contract_worker: EOF before Q\n";
    }
    return false;
}

int Usage()
{
    std::cerr << "usage: wirehair_wh2_contract_worker --describe | "
        "--emit-traces recovery|phase-attribution | "
        "--emit-timing-proxy-witness | "
        "--worker CPU | "
        "--describe-recovery-candidate ID | "
        "--recovery-candidate-worker ID CPU\n";
    return 2;
}

} // namespace

int main(int argc, char** argv)
{
    std::string binary_sha256;
    const std::string source_git_commit =
        WIREHAIR_WH2_SOURCE_GIT_COMMIT;
    std::vector<ArmDescriptor> arms;
    if (!SelfBinarySha256(binary_sha256) ||
        !IsLowerGitCommit(source_git_commit) ||
        !BuildArmDescriptors(arms))
    {
        std::cerr << "wirehair_wh2_contract_worker: cannot establish "
            "binary/descriptor identity\n";
        return 1;
    }
    if (argc == 2 && std::string(argv[1]) == "--describe")
    {
        const std::string description = WorkerDescription(
            binary_sha256, source_git_commit, arms);
        return !description.empty() && EmitLine(description) ? 0 : 1;
    }
    if (argc == 3 && std::string(argv[1]) == "--emit-traces")
    {
        const std::string kind(argv[2]);
        if (kind == "recovery") return EmitRecoveryTraces() ? 0 : 1;
        if (kind == "phase-attribution") {
            return EmitPhaseAttributionTraces(arms) ? 0 : 1;
        }
        return Usage();
    }
    if (argc == 2 &&
        std::string(argv[1]) == "--emit-timing-proxy-witness")
    {
        return EmitTimingProxyWitness(
            binary_sha256, source_git_commit, arms) ? 0 : 1;
    }
    if (argc == 3 &&
        std::string(argv[1]) == "--describe-recovery-candidate")
    {
        const CandidateDefinition* const candidate =
            FindRecoveryCandidate(argv[2]);
        std::vector<ArmDescriptor> candidate_arms;
        if (!candidate) {
            return Usage();
        }
        if (!BuildCandidateArmDescriptors(*candidate, candidate_arms)) {
            return 1;
        }
        const std::string description = RawWorkerDescription(
            binary_sha256, source_git_commit, candidate_arms);
        return !description.empty() && EmitLine(description) ? 0 : 1;
    }
    if (argc == 3 && std::string(argv[1]) == "--worker")
    {
        int cpu = -1;
        if (!ParseCpu(argv[2], cpu)) {
            return Usage();
        }
        return RunWorker(cpu, binary_sha256, arms, nullptr) ? 0 : 1;
    }
    if (argc == 4 &&
        std::string(argv[1]) == "--recovery-candidate-worker")
    {
        const CandidateDefinition* const candidate =
            FindRecoveryCandidate(argv[2]);
        int cpu = -1;
        std::vector<ArmDescriptor> candidate_arms;
        if (!candidate || !ParseCpu(argv[3], cpu)) {
            return Usage();
        }
        if (!BuildCandidateArmDescriptors(*candidate, candidate_arms)) {
            return 1;
        }
        return RunWorker(
            cpu, binary_sha256, candidate_arms, candidate) ? 0 : 1;
    }
    return Usage();
}
