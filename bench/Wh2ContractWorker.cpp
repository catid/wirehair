#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "Wh2FrozenTiming.h"
#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"
#include "Wh2NativePanel.h"

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
#include <sstream>
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
using wirehair_wh2_bench::RecoveryCellResult;
using wirehair_wh2_bench::RecoveryOutcome;

static const char kTraceSchema[] =
    "wirehair.wh2.native-trace-record.v1";
static const char kRecoverySchema[] =
    "wirehair.wh2.native-recovery-record.v1";
static const char kTimingSchema[] =
    "wirehair.wh2.native-timing-record.v2";
static const char kDescriptionSchema[] =
    "wirehair.wh2.native-worker-description.v1";
static const char kDescriptorSchema[] =
    "wirehair.wh2.native-arm-descriptor.v1";
static const char kWorkSchema[] = "wirehair.wh2.native-work.v1";
static const char kRealizedSchema[] =
    "wirehair.wh2.realized-construction.v1";
static const char kRealizedDomain[] =
    "wirehair.wh2.realized-construction.v1\0";
static const char kSourceDomain[] = "wirehair.wh2.source.v1\0";
static const char kZeroSha256[] =
    "0000000000000000000000000000000000000000000000000000000000000000";
static const uint64_t kWh2ProfileId = UINT64_C(0x4b295bbb47f4f9c9);
static const uint16_t kWh2ProfileEncodingVersion = 1u;
static const uint32_t kRecoveryOverheadCap = 4u;

struct ArmDescriptor
{
    const char* Arm;
    const char* Codec;
    std::string CanonicalDescriptor;
    std::string DescriptorSha256;
};

struct WorkerContext
{
    int Cpu;
    int Pid;
    uint64_t ProcessStartTicks;
    std::string BinarySha256;
    std::vector<ArmDescriptor> Arms;
};

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
        "wirehair2_head", "wirehair1", "wirehair2_identity"
    };
    static const char* const codecs[] = {
        "wirehair2_certified", "wirehair1", "wirehair2_experiment"
    };
    static const char* const transforms[] = {
        "none", "none", "identity"
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
        if (!IsLowerSha256(descriptor.DescriptorSha256)) {
            return false;
        }
        built.push_back(descriptor);
    }
    arms.swap(built);
    return true;
}

bool IdentityTransform(
    uint32_t,
    uint32_t,
    wirehair_v2::PrecodeParams&,
    wirehair_v2::PacketRowConfig&,
    void*)
{
    return true;
}

bool ArmSpecFor(
    std::size_t arm_index,
    uint32_t construction_attempt,
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
        spec = wirehair_wh2_bench::MakeExperimentalWh2Arm(
            construction_attempt, &IdentityTransform);
        return true;
    default:
        return false;
    }
}

uint32_t ConstructionAttemptFor(std::size_t arm_index, uint32_t base_attempt)
{
    return arm_index == 1u ? 0u : base_attempt;
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
    if (line.size() < 5u || (line[0] != 'R' && line[0] != 'T') ||
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

bool EmitTimingTraces()
{
    const std::vector<FrozenTimingCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentTimingCells();
    if (cells.size() != 192u) {
        return false;
    }
    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        FrozenPacketTrace trace;
        wirehair::wh2_benchmark::FrozenTimingTraceReceipt receipt;
        if (cells[i].ordinal != i ||
            wirehair::wh2_benchmark::GenerateDevelopmentTimingTrace(
                cells[i], trace, receipt) != FrozenTraceStatus::Complete ||
            receipt.ordinal != i ||
            trace.delivered_ids.size() !=
                static_cast<std::size_t>(cells[i].K) + 4u)
        {
            return false;
        }
        if (!IsLowerSha256(receipt.cell_sha256) ||
            !IsLowerSha256(trace.trace_sha256) ||
            receipt.trace_sha256 != trace.trace_sha256 ||
            !EmitLine(TraceRecord(i, receipt.cell_sha256, trace)))
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
    if (arms.size() != 3u || !IsLowerGitCommit(source_git_commit)) {
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
    if (!ArmSpecFor(arm_index, attempt, spec)) {
        error = "cannot construct selected native recovery arm";
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
    if (!RealizedConstructionSha256(
            context.Arms[arm_index], cell.K, cell.block_bytes, attempt,
            realized_sha256))
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
    const std::string payload = RecoveryPayload(
        cell, context.Arms[arm_index], context.BinarySha256, cell_sha256,
        trace.trace_sha256, attempt, realized_sha256, outcome,
        has_decoded_extra, decoded_extra);
    const std::size_t ordinal = cell_ordinal * context.Arms.size() + arm_index;
    const std::string work_sha256 = RecoveryWorkSha256(
        cell_sha256, ordinal);
    if (!IsLowerSha256(work_sha256)) {
        error = "cannot hash recovery work identity";
        return false;
    }
    envelope = Envelope(
        kRecoverySchema, ordinal, actual_cpu, context.Pid, started_ns,
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
        uint32_t fixed_received_overhead)
        : IdentityValue(identity)
        , Scope(scope)
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
            PreparationResult = Receive->Initialize(arm, packet_ids);
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
                result.BytesVerified && result.DecodedOverhead <= 4u;
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
                result.Result == Wirehair_Success, 4u, true);
        }
        return ClassifyTimedResult(
            Wirehair_InvalidInput, 0u, false, 0u, false);
    }

private:
    std::string IdentityValue;
    FrozenTimingScope Scope;
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
    json += ",\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"order\":\"";
    json += order;
    json += "\",\"panel_kind\":\"";
    json += wirehair::wh2_benchmark::FrozenPanelKindName(panel.panel_kind);
    json += "\",\"phase\":\"";
    json += cell.phase;
    json += "\",\"replicate\":";
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

bool RunTimingJob(
    const WorkerContext& context,
    std::size_t cell_ordinal,
    std::size_t panel_index,
    std::string& envelope,
    std::string& error)
{
    const std::vector<FrozenTimingCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentTimingCells();
    const std::vector<FrozenTimingPanel> panels =
        wirehair::wh2_benchmark::EnumerateOneCandidateTimingPanels(
            "wirehair2_identity");
    if (cells.size() != 192u || panels.size() != 11u ||
        cell_ordinal >= cells.size() || panel_index >= panels.size())
    {
        error = "timing job index is outside the frozen domain";
        return false;
    }
    const FrozenTimingCell& cell = cells[cell_ordinal];
    const FrozenTimingPanel& panel = panels[panel_index];
    if (cell.ordinal != cell_ordinal || panel.ordinal != panel_index ||
        cell.fixed_received_overhead != 4u ||
        cell.invocations_per_slot == 0u ||
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
        trace.delivered_ids.size() != static_cast<std::size_t>(cell.K) + 4u)
    {
        error = "cannot regenerate complete frozen timing trace";
        return false;
    }
    const std::string cell_json =
        wirehair::wh2_benchmark::CanonicalTimingCellJson(cell);
    const std::string cell_sha256 =
        wirehair::wh2_benchmark::TimingCellSha256(cell);
    uint64_t source_seed = 0u;
    std::vector<uint8_t> source;
    if (cell_json.empty() || trace_receipt.cell_sha256 != cell_sha256 ||
        trace_receipt.trace_sha256 != trace.trace_sha256 ||
        !IsLowerSha256(cell_sha256) || !IsLowerSha256(trace.trace_sha256) ||
        !SourceSeedFromCellJson(cell_json, source_seed) ||
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
    if (!ArmSpecFor(left_index, left_attempt, left_spec) ||
        !ArmSpecFor(right_index, right_attempt, right_spec))
    {
        error = "cannot construct timing arm specifications";
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
                    cell.fixed_received_overhead));
        });
    const NativePanelArm right_arm(
        right_identity,
        [right_identity, panel, right_spec, &cell, &source, &trace]() {
            return std::unique_ptr<NativePanelInvocation>(
                new CodecTimingInvocation(
                    right_identity, panel.scope, right_spec, cell.K,
                    cell.block_bytes, source, trace.delivered_ids,
                    cell.fixed_received_overhead));
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
    const auto valid_extra = [panel](
        bool success, bool has_extra, uint32_t extra) -> bool
    {
        if (!success) return !has_extra;
        if (panel.scope ==
            FrozenTimingScope::EncoderInitPlusFirstKSymbols)
        {
            return !has_extra;
        }
        if (panel.scope == FrozenTimingScope::DecoderSolve) {
            return has_extra && extra == 4u;
        }
        return panel.scope == FrozenTimingScope::ReceiveToSuccess &&
            has_extra && extra <= 4u;
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
        result.InvocationsPerSlot == cell.invocations_per_slot;
    for (std::size_t slot = 0u; slot < result.Slots.size(); ++slot)
    {
        const bool left_slot = native_order == NativePanelOrder::ABBA ?
            (slot == 0u || slot == 3u) : (slot == 1u || slot == 2u);
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
        (both_success && result.HasFourNullTimings()) ||
        (!both_success && !result.HasFourNullTimings()))
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

bool RunWorker(int cpu, const std::string& binary_sha256,
    const std::vector<ArmDescriptor>& arms)
{
    std::string error;
    if (!PinSingletonCpu(cpu, error))
    {
        std::cerr << "wirehair_wh2_contract_worker: " << error << '\n';
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
        std::string output;
        bool ok = false;
        if (kind == 'R') {
            ok = RunRecoveryJob(
                context, cell_ordinal, item_index, output, error);
        }
        else {
            ok = RunTimingJob(
                context, cell_ordinal, item_index, output, error);
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
        "--emit-traces recovery|timing | --worker CPU\n";
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
        if (kind == "timing") return EmitTimingTraces() ? 0 : 1;
        return Usage();
    }
    if (argc == 3 && std::string(argv[1]) == "--worker")
    {
        int cpu = -1;
        if (!ParseCpu(argv[2], cpu)) {
            return Usage();
        }
        return RunWorker(cpu, binary_sha256, arms) ? 0 : 1;
    }
    return Usage();
}
