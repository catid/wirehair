#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <fcntl.h>
#include <sched.h>
#include <setjmp.h>
#include <signal.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#if defined(__linux__) && defined(__x86_64__) && \
    (defined(__GNUC__) || defined(__clang__))
#include <cpuid.h>
#define WIREHAIR_RDPRU_COMPILED 1
#else
#define WIREHAIR_RDPRU_COMPILED 0
#endif

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#define WIREHAIR_WH2_SOURCE_GIT_COMMIT "unknown"
#endif

#ifndef WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256
#define WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256 "unknown"
#endif

#ifndef WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256
#define WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256 "unknown"
#endif

#ifndef WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256
#define WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256 "unknown"
#endif

#ifndef WIREHAIR_WH2_CMAKE_SHA256
#define WIREHAIR_WH2_CMAKE_SHA256 "unknown"
#endif

namespace {

#if WIREHAIR_RDPRU_COMPILED
extern "C" const unsigned char wh2_rdpru_bracket_opcode_begin[];
extern "C" const unsigned char wh2_rdpru_bracket_opcode_end[];
#endif

using wirehair_wh2_bench::Batch128CounterArmResult;
using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeCounterBracket;
using wirehair_wh2_bench::NativeCounterReadFunction;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativeReceiveFixture;

static const char kProtocolTag[] =
    "wirehair1-native-rdpru-batch128-aa-v1";
#if defined(__linux__)
static const char kLaunchManifestSchema[] =
    "wirehair-rdpru-launch-manifest-v1";
#endif
static const uint32_t kBracketProtocolId = 1u;
static const char kBracketProtocol[] =
    "compiler-barrier+mfence+lfence+rdpru(selector)+lfence+compiler-barrier-v1";
static const char kHarmlessProbeSourceSha256[] =
    "62b72bcd488771d1ed37d86f3222a454ef0b2043c08f1771473268050411a69f";
static const char kHarmlessProbeBinarySha256[] =
    "06119dd319bafd423fa1caa5c99b64dd475d9fcef4cdfa73234eba5271b43fae";
static const uint32_t kK = 2u;
static const uint32_t kLossPpm = 100000u;
static const uint32_t kReceiveOverheadCap = 256u;
static const int kFrozenTargetCpu = 50;
static const uint32_t kFrozenLogicalCpuCount = 128u;
static const uint32_t kFrozenPhysicalCoreCount = 64u;
static const int kFrozenFirstOnlineCpu = 0;
static const int kFrozenLastOnlineCpu = 127;
static const int kFrozenTargetSiblingCpu = 114;
static const char kFrozenOnlineCpuList[] = "0-127";
static const char kFrozenInitialAffinityList[] = "0-127";
static const char kFrozenTargetSiblingList[] = "50,114";
static const uint32_t kFrozenCpuFamily = 26u;
static const uint32_t kFrozenCpuModel = 8u;
static const uint32_t kFrozenCpuStepping = 1u;
static const uint32_t kFrozenLeaf1Eax = UINT32_C(0x00b00f81);
static const uint32_t kFrozenLeaf1Ebx = UINT32_C(0x29800800);
static const uint32_t kFrozenLeaf1Ecx = UINT32_C(0x7efa320b);
static const uint32_t kFrozenLeaf1Edx = UINT32_C(0x178bfbff);
static const uint32_t kFrozenLeaf6Eax = UINT32_C(0x00000004);
static const uint32_t kFrozenLeaf6Ebx = UINT32_C(0x00000000);
static const uint32_t kFrozenLeaf6Ecx = UINT32_C(0x00000001);
static const uint32_t kFrozenLeaf6Edx = UINT32_C(0x00000000);
static const uint32_t kFrozenLeaf80000008Eax = UINT32_C(0x00003934);
static const uint32_t kFrozenLeaf80000008Ebx = UINT32_C(0x79bef25f);
static const uint32_t kFrozenLeaf80000008Ecx = UINT32_C(0x0000707f);
static const uint32_t kFrozenLeaf80000008Edx = UINT32_C(0x00010007);
static const uint32_t kFrozenLeaf80000021Eax = UINT32_C(0xd93fffcf);
static const uint32_t kFrozenLeaf80000021Ebx = UINT32_C(0x00080382);
static const uint32_t kFrozenLeaf80000021Ecx = UINT32_C(0x00000000);
static const uint32_t kFrozenLeaf80000021Edx = UINT32_C(0x00000000);
static const char kCpuSelectionRule[] =
    "sum-numeric-proc-interrupts-rows-per-cpu;minimum-first-smt-cpus-0-63;tie-cpu-ascending";
static const uint32_t kCpuSelectionFirstCpu = 0u;
static const uint32_t kCpuSelectionLastCpu = 63u;
static const uint64_t kCpuSelectionTargetInterruptCount = UINT64_C(289);
static const uint32_t kCpuSelectionRunnerUpCpu = 54u;
static const uint64_t kCpuSelectionRunnerUpInterruptCount = UINT64_C(297);
static const char kCpuSelectionSnapshotPath[] =
    "/dev/shm/wirehair-rdpru-cpu-selection-interrupts-pre-run.txt";
static const char kCpuSelectionSnapshotSha256[] =
    "770e0f948c8ba7d7299c1c994e44e45284062a50f1b962b4104707fab6654df5";
static const uint32_t kRoundCount = 32u;
static const uint32_t kCellCount = 18u;
static const uint32_t kBaseCellCount = 6u;
static const uint32_t kModeCount = 2u;
static const uint32_t kSelectorCount = 2u;
static const uint32_t kSequenceCount = 4u;
static const uint32_t kPeriodsPerSequence = 4u;
static const uint32_t kBatchesPerSequence = 8u;
static const uint32_t kMeasuredBatchesPerSelectorSuperblock = 16u;
static const uint32_t kWashoutBatchesPerSelectorSuperblock = 16u;
static const uint32_t kBatchSize = 128u;
static const uint32_t kEmptySampleCount = 4096u;
static const uint32_t kEmptyP999OneBasedIndex = 4092u;
static const uint32_t kEmptyOverheadMultiplier = 200u;
#if defined(__linux__)
static const size_t kMaximumLaunchManifestBytes = 1024u * 1024u;
#endif
static const size_t kOriginalSuperblockCount =
    static_cast<size_t>(kCellCount) * kRoundCount * kModeCount;
static const size_t kSelectorSuperblockCount =
    kOriginalSuperblockCount * kSelectorCount;
static const size_t kMeasuredBatchCountPerSelector =
    kOriginalSuperblockCount * kMeasuredBatchesPerSelectorSuperblock;
static const size_t kWashoutBatchCountPerSelector =
    kOriginalSuperblockCount * kWashoutBatchesPerSelectorSuperblock;
static const size_t kMeasuredBatchCount =
    kMeasuredBatchCountPerSelector * kSelectorCount;
static const size_t kWashoutBatchCount =
    kWashoutBatchCountPerSelector * kSelectorCount;
static const size_t kMeasuredInnerInvocationCount =
    kMeasuredBatchCount * kBatchSize;
static const size_t kWashoutInnerInvocationCount =
    kWashoutBatchCount * kBatchSize;
static const size_t kPreflightBatchCount =
    static_cast<size_t>(kBaseCellCount) * 3u * 3u;
static const size_t kPreflightInnerInvocationCount =
    kPreflightBatchCount * kBatchSize;
static const size_t kWarmBatchCount =
    static_cast<size_t>(kCellCount) * 3u * kSelectorCount;
static const size_t kWarmInnerInvocationCount =
    kWarmBatchCount * kBatchSize;
static const size_t kPanelActionCount =
    kMeasuredBatchCount + kWashoutBatchCount;
static const size_t kOfficialGateCount = 56u;
static const size_t kExpectedN96GateCount = 48u;
static const size_t kExpectedN576GateCount = 8u;
static const double kT95N96 = 3.430611421941;
static const double kT95N576 = 3.339748611060;
// Inward binary64 margins: each is conservative relative to the exact real
// log(0.98)/log(1.02) boundary.
static const double kLogLowerBound = -0.020202707317519445;
static const double kLogUpperBound = 0.019802627296179713;
static const char kLogLowerBoundDecimal[] = "-0.020202707317519445";
static const char kLogUpperBoundDecimal[] = "0.019802627296179713";
static const char kExpectedScheduleSha256[] =
    "1f115260825a6dfb5e9ce272c6d3cc2eae251586dad74f89f1e8a737b7519c66";
static const char kExpectedConfigSha256[] =
    "53205b3081074aea5407fd500f765f9961015582dfeb1d38228461544bd0af57";
static const char kExpectedFixtureSha256[] =
    "fec01bd85496b25c878519375660a3b9a98f427321108f4a83085e2a5eecae76";
static const uint32_t kBlockBytes[] = { 64u, 1280u };
static const uint64_t kRoots[] = {
    UINT64_C(0xd1b54a32d192ed03),
    UINT64_C(0x94d049bb133111eb),
    UINT64_C(0x8538ecb5bd456ea3)
};
static const uint32_t kExpectedIidOverheads[kBaseCellCount] = {
    0u, 0u, 0u, 0u, 0u, 0u
};

enum class Scope : uint32_t
{
    Encoder = 0u,
    IidReceive = 1u,
    TwoRepairReceive = 2u
};

enum class Mode : uint32_t
{
    Exposed = 0u,
    Isolated = 1u
};

enum class Estimand : uint32_t
{
    Logical = 0u,
    Physical = 1u
};

enum class LogicalLabel : uint32_t
{
    A = 0u,
    B = 1u
};

enum class Direction : uint32_t
{
    Forward = 0u,
    Reverse = 1u
};

enum class PhysicalBank : uint32_t
{
    X = 0u,
    Y = 1u,
    W = 2u
};

enum class PredecessorClass : uint32_t
{
    Washout = 0u,
    Aa = 1u,
    Ab = 2u,
    Ba = 3u,
    Bb = 4u
};

static const LogicalLabel kSequences[kSequenceCount][kPeriodsPerSequence] = {
    { LogicalLabel::A, LogicalLabel::B, LogicalLabel::B, LogicalLabel::A },
    { LogicalLabel::B, LogicalLabel::A, LogicalLabel::A, LogicalLabel::B },
    { LogicalLabel::A, LogicalLabel::A, LogicalLabel::B, LogicalLabel::B },
    { LogicalLabel::B, LogicalLabel::B, LogicalLabel::A, LogicalLabel::A }
};

const char* ScopeName(Scope scope)
{
    switch (scope)
    {
    case Scope::Encoder: return "encoder-init-systematic";
    case Scope::IidReceive: return "iid10-receive-to-success";
    case Scope::TwoRepairReceive: return "two-repair-receive-to-success";
    }
    return "invalid";
}

const char* ModeName(Mode mode)
{
    return mode == Mode::Exposed ? "exposed" : "isolated";
}

const char* EstimandName(Estimand estimand)
{
    return estimand == Estimand::Logical ? "logical-a-b" : "physical-x-y";
}

const char* SelectorName(uint32_t selector)
{
    return selector == 0u ? "mperf" : (selector == 1u ? "aperf" : "invalid");
}

const char* LabelName(LogicalLabel label)
{
    return label == LogicalLabel::A ? "A" : "B";
}

const char* BankName(PhysicalBank bank)
{
    switch (bank)
    {
    case PhysicalBank::X: return "X";
    case PhysicalBank::Y: return "Y";
    case PhysicalBank::W: return "W";
    }
    return "invalid";
}

const char* DirectionName(Direction direction)
{
    return direction == Direction::Forward ? "forward" : "reverse";
}

void AppendLe32(std::string& bytes, uint32_t value)
{
    for (unsigned shift = 0u; shift != 32u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendLe64(std::string& bytes, uint64_t value)
{
    for (unsigned shift = 0u; shift != 64u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendSizedString(std::string& bytes, const std::string& value)
{
    AppendLe64(bytes, static_cast<uint64_t>(value.size()));
    bytes.append(value);
}

bool AddU64(uint64_t value, uint64_t& total)
{
    if (value > std::numeric_limits<uint64_t>::max() - total) return false;
    total += value;
    return true;
}

uint64_t DoubleBits(double value)
{
    uint64_t bits = 0u;
    static_assert(sizeof(bits) == sizeof(value), "binary64 size mismatch");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

bool IsLowerHexSha256(const char* text)
{
    if (!text || std::strlen(text) != 64u) return false;
    for (size_t i = 0u; i < 64u; ++i)
    {
        const char c = text[i];
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) {
            return false;
        }
    }
    return true;
}

bool IsAllZeroHex(const char* text)
{
    if (!text || !*text) return false;
    for (const char* next = text; *next; ++next) {
        if (*next != '0') return false;
    }
    return true;
}

bool IsLowerHex(const std::string& text)
{
    if (text.empty()) return false;
    for (size_t i = 0u; i < text.size(); ++i) {
        if (!((text[i] >= '0' && text[i] <= '9') ||
              (text[i] >= 'a' && text[i] <= 'f'))) return false;
    }
    return true;
}

bool EmbeddedIdentitiesValid()
{
    return std::strcmp(WIREHAIR_WH2_SOURCE_GIT_COMMIT, "unknown") != 0 &&
        std::strlen(WIREHAIR_WH2_SOURCE_GIT_COMMIT) == 40u &&
        IsLowerHex(std::string(WIREHAIR_WH2_SOURCE_GIT_COMMIT)) &&
        !IsAllZeroHex(WIREHAIR_WH2_SOURCE_GIT_COMMIT) &&
        IsLowerHexSha256(WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256) &&
        !IsAllZeroHex(WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256) &&
        IsLowerHexSha256(WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256) &&
        !IsAllZeroHex(WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256) &&
        IsLowerHexSha256(WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256) &&
        !IsAllZeroHex(WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256) &&
        IsLowerHexSha256(WIREHAIR_WH2_CMAKE_SHA256) &&
        !IsAllZeroHex(WIREHAIR_WH2_CMAKE_SHA256);
}

struct SuperblockPlan
{
    uint32_t Cell = 0u;
    uint32_t Round = 0u;
    uint32_t Visit = 0u;
    uint32_t ModePosition = 0u;
    Mode ModeValue = Mode::Exposed;
    uint32_t U = 0u;
    Direction DirectionValue = Direction::Forward;
    uint32_t OrientationFlip = 0u;
    uint32_t Rotation = 0u;
    uint32_t FirstSelector = 0u;
};

typedef std::array<SuperblockPlan, kOriginalSuperblockCount> FrozenSchedule;

uint32_t CellForVisit(uint32_t round, uint32_t visit)
{
    if (((round >> 1u) & 1u) == 0u) {
        return (5u * visit + 7u * round) % kCellCount;
    }
    return (5u * (kCellCount - 1u - visit) + 7u * round) % kCellCount;
}

uint32_t ScheduleU(uint32_t cell, uint32_t round)
{
    return (13u * round + 11u * cell + 5u) % kRoundCount;
}

SuperblockPlan MakePlan(
    uint32_t cell,
    uint32_t round,
    uint32_t visit,
    uint32_t mode_position,
    Mode mode)
{
    SuperblockPlan plan;
    plan.Cell = cell;
    plan.Round = round;
    plan.Visit = visit;
    plan.ModePosition = mode_position;
    plan.ModeValue = mode;
    plan.U = ScheduleU(cell, round);
    plan.DirectionValue = ((round >> 1u) & 1u) == 0u ?
        Direction::Forward : Direction::Reverse;
    plan.OrientationFlip = mode == Mode::Exposed ?
        ((plan.U >> 1u) & 1u) : ((plan.U >> 2u) & 1u);
    const uint32_t q = (plan.U >> 3u) & 3u;
    plan.Rotation = mode == Mode::Exposed ? q :
        (q + 1u + 2u * (plan.U & 1u)) % kSequenceCount;
    // Five-bit parity avoids aliasing selector order with the crossed
    // orientation, mode-position, direction, q/rotation and parity factors.
    // Round/visit long-span residuals remain explicit diagnostics.
    plan.FirstSelector =
        (plan.U ^ (plan.U >> 1u) ^ (plan.U >> 2u) ^
            (plan.U >> 3u) ^ (plan.U >> 4u)) & 1u;
    return plan;
}

bool BuildFrozenSchedule(FrozenSchedule& schedule)
{
    size_t next = 0u;
    for (uint32_t round = 0u; round < kRoundCount; ++round)
    {
        for (uint32_t visit = 0u; visit < kCellCount; ++visit)
        {
            const uint32_t cell = CellForVisit(round, visit);
            const uint32_t u = ScheduleU(cell, round);
            const Mode first = (u & 1u) == 0u ?
                Mode::Exposed : Mode::Isolated;
            const Mode second = first == Mode::Exposed ?
                Mode::Isolated : Mode::Exposed;
            if (next + 2u > schedule.size()) return false;
            schedule[next++] = MakePlan(cell, round, visit, 0u, first);
            schedule[next++] = MakePlan(cell, round, visit, 1u, second);
        }
    }
    return next == schedule.size();
}

PhysicalBank LogicalAPhysicalForSequence(
    const SuperblockPlan& plan,
    uint32_t sequence)
{
    const uint32_t family = sequence < 2u ? 0u : 1u;
    return (family ^ plan.OrientationFlip) == 0u ?
        PhysicalBank::X : PhysicalBank::Y;
}

std::string ScheduleSha256(const FrozenSchedule& schedule)
{
    std::string bytes;
    static const char marker[] = "wirehair-rdpru-schedule-le-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    AppendLe64(bytes, static_cast<uint64_t>(schedule.size()));
    uint64_t global_action = 0u;
    uint64_t selector_actions[2] = { 0u, 0u };
    for (size_t i = 0u; i < schedule.size(); ++i)
    {
        const SuperblockPlan& plan = schedule[i];
        for (uint32_t selector_position = 0u;
             selector_position < kSelectorCount;
             ++selector_position)
        {
            const uint32_t selector = plan.FirstSelector ^ selector_position;
            for (uint32_t sequence_position = 0u;
                 sequence_position < kSequenceCount;
                 ++sequence_position)
            {
                const uint32_t sequence =
                    (plan.Rotation + sequence_position) % kSequenceCount;
                const PhysicalBank logical_a =
                    LogicalAPhysicalForSequence(plan, sequence);
                const PhysicalBank logical_b = logical_a == PhysicalBank::X ?
                    PhysicalBank::Y : PhysicalBank::X;
                for (uint32_t step = 0u; step < kBatchesPerSequence; ++step)
                {
                    const bool measured = plan.ModeValue == Mode::Exposed ?
                        step >= kPeriodsPerSequence : (step & 1u) != 0u;
                    const uint32_t period = plan.ModeValue == Mode::Exposed ?
                        (measured ? step - kPeriodsPerSequence : step) :
                        step / 2u;
                    const LogicalLabel label = measured ?
                        kSequences[sequence][period] : LogicalLabel::A;
                    const PhysicalBank physical = !measured ?
                        PhysicalBank::W :
                        (label == LogicalLabel::A ? logical_a : logical_b);
                    AppendLe64(bytes, global_action++);
                    AppendLe64(bytes, selector_actions[selector]++);
                    AppendLe64(bytes, static_cast<uint64_t>(i));
                    AppendLe32(bytes, plan.Cell);
                    AppendLe32(bytes, plan.Round);
                    AppendLe32(bytes, plan.Visit);
                    AppendLe32(bytes, plan.U);
                    AppendLe32(bytes, static_cast<uint32_t>(plan.DirectionValue));
                    AppendLe32(bytes, static_cast<uint32_t>(plan.ModeValue));
                    AppendLe32(bytes, plan.ModePosition);
                    AppendLe32(bytes, plan.OrientationFlip);
                    AppendLe32(bytes, plan.Rotation);
                    AppendLe32(bytes, plan.FirstSelector);
                    AppendLe32(bytes, selector_position);
                    AppendLe32(bytes, selector);
                    AppendLe32(bytes, sequence_position);
                    AppendLe32(bytes, sequence);
                    AppendLe32(bytes, period);
                    AppendLe32(bytes, measured ? 1u : 0u);
                    AppendLe32(bytes, static_cast<uint32_t>(label));
                    AppendLe32(bytes, static_cast<uint32_t>(physical));
                    AppendLe32(bytes, kBatchSize);
                }
            }
        }
    }
    if (global_action != kPanelActionCount ||
        selector_actions[0] != kMeasuredBatchCountPerSelector +
            kWashoutBatchCountPerSelector ||
        selector_actions[1] != kMeasuredBatchCountPerSelector +
            kWashoutBatchCountPerSelector)
    {
        return std::string();
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

std::string ConfigSha256(const std::string& schedule_sha256)
{
    std::string bytes;
    bytes.append(kProtocolTag, sizeof(kProtocolTag) - 1u);
    AppendLe32(bytes, kK);
    AppendLe32(bytes, kLossPpm);
    AppendLe32(bytes, kReceiveOverheadCap);
    AppendLe32(bytes, static_cast<uint32_t>(kFrozenTargetCpu));
    AppendLe32(bytes, kFrozenLogicalCpuCount);
    AppendLe32(bytes, kFrozenPhysicalCoreCount);
    AppendLe32(bytes, static_cast<uint32_t>(kFrozenFirstOnlineCpu));
    AppendLe32(bytes, static_cast<uint32_t>(kFrozenLastOnlineCpu));
    AppendLe32(bytes, static_cast<uint32_t>(kFrozenTargetSiblingCpu));
    AppendSizedString(bytes, kFrozenOnlineCpuList);
    AppendSizedString(bytes, kFrozenInitialAffinityList);
    AppendSizedString(bytes, kFrozenTargetSiblingList);
    AppendLe32(bytes, kFrozenCpuFamily);
    AppendLe32(bytes, kFrozenCpuModel);
    AppendLe32(bytes, kFrozenCpuStepping);
    AppendLe32(bytes, kFrozenLeaf1Eax);
    AppendLe32(bytes, kFrozenLeaf1Ebx);
    AppendLe32(bytes, kFrozenLeaf1Ecx);
    AppendLe32(bytes, kFrozenLeaf1Edx);
    AppendLe32(bytes, kFrozenLeaf6Eax);
    AppendLe32(bytes, kFrozenLeaf6Ebx);
    AppendLe32(bytes, kFrozenLeaf6Ecx);
    AppendLe32(bytes, kFrozenLeaf6Edx);
    AppendLe32(bytes, kFrozenLeaf80000008Eax);
    AppendLe32(bytes, kFrozenLeaf80000008Ebx);
    AppendLe32(bytes, kFrozenLeaf80000008Ecx);
    AppendLe32(bytes, kFrozenLeaf80000008Edx);
    AppendLe32(bytes, kFrozenLeaf80000021Eax);
    AppendLe32(bytes, kFrozenLeaf80000021Ebx);
    AppendLe32(bytes, kFrozenLeaf80000021Ecx);
    AppendLe32(bytes, kFrozenLeaf80000021Edx);
    AppendSizedString(bytes, kCpuSelectionRule);
    AppendLe32(bytes, kCpuSelectionFirstCpu);
    AppendLe32(bytes, kCpuSelectionLastCpu);
    AppendLe32(bytes, static_cast<uint32_t>(kFrozenTargetCpu));
    AppendLe64(bytes, kCpuSelectionTargetInterruptCount);
    AppendLe32(bytes, kCpuSelectionRunnerUpCpu);
    AppendLe64(bytes, kCpuSelectionRunnerUpInterruptCount);
    AppendSizedString(bytes, kCpuSelectionSnapshotPath);
    AppendSizedString(bytes, kCpuSelectionSnapshotSha256);
    AppendLe32(bytes, kRoundCount);
    AppendLe32(bytes, kCellCount);
    AppendLe32(bytes, kSelectorCount);
    AppendLe32(bytes, kBatchSize);
    AppendLe32(bytes, kEmptySampleCount);
    AppendLe32(bytes, kEmptyP999OneBasedIndex);
    AppendLe32(bytes, kEmptyOverheadMultiplier);
    AppendLe32(bytes, kBracketProtocolId);
    AppendSizedString(bytes, kBracketProtocol);
    AppendLe64(bytes, static_cast<uint64_t>(kOriginalSuperblockCount));
    AppendLe64(bytes, static_cast<uint64_t>(kSelectorSuperblockCount));
    AppendLe64(bytes, static_cast<uint64_t>(kMeasuredBatchCount));
    AppendLe64(bytes, static_cast<uint64_t>(kWashoutBatchCount));
    AppendLe64(bytes, static_cast<uint64_t>(kMeasuredInnerInvocationCount));
    AppendLe64(bytes, static_cast<uint64_t>(kWashoutInnerInvocationCount));
    AppendLe64(bytes, static_cast<uint64_t>(kPreflightBatchCount));
    AppendLe64(bytes, static_cast<uint64_t>(kWarmBatchCount));
    AppendLe64(bytes, static_cast<uint64_t>(kOfficialGateCount));
    AppendLe64(bytes, DoubleBits(kT95N96));
    AppendLe64(bytes, DoubleBits(kT95N576));
    AppendLe64(bytes, DoubleBits(kLogLowerBound));
    AppendLe64(bytes, DoubleBits(kLogUpperBound));
    AppendSizedString(bytes, "tail=1/2240;p=2239/2240;interval=[log(.98),log(1.02)];48xn96;8xn576");
    AppendSizedString(bytes, schedule_sha256);
    for (size_t i = 0u; i < sizeof(kBlockBytes) / sizeof(kBlockBytes[0]); ++i) {
        AppendLe32(bytes, kBlockBytes[i]);
    }
    for (size_t i = 0u; i < sizeof(kRoots) / sizeof(kRoots[0]); ++i) {
        AppendLe64(bytes, kRoots[i]);
    }
    for (size_t i = 0u; i < kBaseCellCount; ++i) {
        AppendLe32(bytes, kExpectedIidOverheads[i]);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

struct CpuFeatures
{
    bool Compiled = false;
    bool VendorAuthenticAmd = false;
    bool Hypervisor = false;
    bool HasLeaf6 = false;
    bool AperfMperf = false;
    bool HasLeaf80000008 = false;
    bool Rdpru = false;
    uint32_t RdpruMax = 0u;
    bool HasLeaf80000021 = false;
    bool LfenceAlwaysSerializing = false;
    uint32_t Family = 0u;
    uint32_t Model = 0u;
    uint32_t Stepping = 0u;
    uint32_t Leaf1Eax = 0u;
    uint32_t Leaf1Ebx = 0u;
    uint32_t Leaf1Ecx = 0u;
    uint32_t Leaf1Edx = 0u;
    uint32_t Leaf6Eax = 0u;
    uint32_t Leaf6Ebx = 0u;
    uint32_t Leaf6Ecx = 0u;
    uint32_t Leaf6Edx = 0u;
    uint32_t Leaf80000008Eax = 0u;
    uint32_t Leaf80000008Ebx = 0u;
    uint32_t Leaf80000008Ecx = 0u;
    uint32_t Leaf80000008Edx = 0u;
    uint32_t Leaf80000021Eax = 0u;
    uint32_t Leaf80000021Ebx = 0u;
    uint32_t Leaf80000021Ecx = 0u;
    uint32_t Leaf80000021Edx = 0u;
};

void DecodeDisplayFamilyModel(
    uint32_t leaf1_eax,
    uint32_t& family_out,
    uint32_t& model_out,
    uint32_t& stepping_out)
{
    stepping_out = leaf1_eax & 0xfu;
    const uint32_t base_family = (leaf1_eax >> 8u) & 0xfu;
    const uint32_t base_model = (leaf1_eax >> 4u) & 0xfu;
    const uint32_t extended_family = (leaf1_eax >> 20u) & 0xffu;
    const uint32_t extended_model = (leaf1_eax >> 16u) & 0xfu;
    family_out = base_family == 0xfu ?
        base_family + extended_family : base_family;
    model_out = base_model;
    if (base_family == 0x6u || base_family == 0xfu) {
        model_out += extended_model << 4u;
    }
}

bool ValidateCpuFeatures(const CpuFeatures& features, std::string& diagnostic)
{
    if (!features.Compiled) {
        diagnostic = "RDPRU controller requires Linux x86_64 GNU-compatible assembly";
        return false;
    }
    if (!features.VendorAuthenticAmd) {
        diagnostic = "CPUID vendor is not AuthenticAMD";
        return false;
    }
    if (features.Hypervisor) {
        diagnostic = "CPUID reports a hypervisor; bare metal is required";
        return false;
    }
    if (!features.HasLeaf6 || !features.AperfMperf) {
        diagnostic = "CPUID does not advertise APERF/MPERF";
        return false;
    }
    if (!features.HasLeaf80000008 || !features.Rdpru) {
        diagnostic = "CPUID does not advertise RDPRU";
        return false;
    }
    if (features.RdpruMax < 1u) {
        diagnostic = "CPUID RdpruMax does not include selector 1";
        return false;
    }
    if (!features.HasLeaf80000021 || !features.LfenceAlwaysSerializing) {
        diagnostic = "CPUID does not guarantee always-serializing LFENCE dispatch";
        return false;
    }
    if (features.Family != kFrozenCpuFamily ||
        features.Model != kFrozenCpuModel ||
        features.Stepping != kFrozenCpuStepping ||
        features.Leaf1Eax != kFrozenLeaf1Eax ||
        features.Leaf1Ebx != kFrozenLeaf1Ebx ||
        features.Leaf1Ecx != kFrozenLeaf1Ecx ||
        features.Leaf1Edx != kFrozenLeaf1Edx ||
        features.Leaf6Eax != kFrozenLeaf6Eax ||
        features.Leaf6Ebx != kFrozenLeaf6Ebx ||
        features.Leaf6Ecx != kFrozenLeaf6Ecx ||
        features.Leaf6Edx != kFrozenLeaf6Edx ||
        features.Leaf80000008Eax != kFrozenLeaf80000008Eax ||
        features.Leaf80000008Ebx != kFrozenLeaf80000008Ebx ||
        features.Leaf80000008Ecx != kFrozenLeaf80000008Ecx ||
        features.Leaf80000008Edx != kFrozenLeaf80000008Edx ||
        features.Leaf80000021Eax != kFrozenLeaf80000021Eax ||
        features.Leaf80000021Ebx != kFrozenLeaf80000021Ebx ||
        features.Leaf80000021Ecx != kFrozenLeaf80000021Ecx ||
        features.Leaf80000021Edx != kFrozenLeaf80000021Edx)
    {
        diagnostic = "CPUID differs from the frozen CPU50 host identity";
        return false;
    }
    diagnostic.clear();
    return true;
}

bool ReadCpuFeaturesRaw(CpuFeatures& features, std::string& diagnostic)
{
    CpuFeatures value;
#if WIREHAIR_RDPRU_COMPILED
    value.Compiled = true;
    unsigned vendor_signature = 0u;
    const unsigned max_basic = __get_cpuid_max(0u, &vendor_signature);
    if (max_basic == 0u) {
        features = value;
        diagnostic = "basic CPUID is unavailable";
        return false;
    }
    unsigned eax = 0u;
    unsigned ebx = 0u;
    unsigned ecx = 0u;
    unsigned edx = 0u;
    __cpuid_count(0u, 0u, eax, ebx, ecx, edx);
    char vendor[13];
    std::memcpy(vendor + 0, &ebx, 4u);
    std::memcpy(vendor + 4, &edx, 4u);
    std::memcpy(vendor + 8, &ecx, 4u);
    vendor[12] = '\0';
    value.VendorAuthenticAmd = std::strcmp(vendor, "AuthenticAMD") == 0;
    if (max_basic >= 1u)
    {
        __cpuid_count(1u, 0u, eax, ebx, ecx, edx);
        value.Leaf1Eax = eax;
        value.Leaf1Ebx = ebx;
        value.Leaf1Ecx = ecx;
        value.Leaf1Edx = edx;
        value.Hypervisor = (ecx & (UINT32_C(1) << 31u)) != 0u;
        DecodeDisplayFamilyModel(
            eax, value.Family, value.Model, value.Stepping);
    }
    value.HasLeaf6 = max_basic >= 6u;
    if (value.HasLeaf6)
    {
        __cpuid_count(6u, 0u, eax, ebx, ecx, edx);
        value.Leaf6Eax = eax;
        value.Leaf6Ebx = ebx;
        value.Leaf6Ecx = ecx;
        value.Leaf6Edx = edx;
        value.AperfMperf = (ecx & 1u) != 0u;
    }
    const unsigned max_extended = __get_cpuid_max(0x80000000u, nullptr);
    value.HasLeaf80000008 = max_extended >= 0x80000008u;
    if (value.HasLeaf80000008)
    {
        __cpuid_count(0x80000008u, 0u, eax, ebx, ecx, edx);
        value.Leaf80000008Eax = eax;
        value.Leaf80000008Ebx = ebx;
        value.Leaf80000008Ecx = ecx;
        value.Leaf80000008Edx = edx;
        value.Rdpru = (ebx & (UINT32_C(1) << 4u)) != 0u;
        value.RdpruMax = edx >> 16u;
    }
    value.HasLeaf80000021 = max_extended >= 0x80000021u;
    if (value.HasLeaf80000021)
    {
        __cpuid_count(0x80000021u, 0u, eax, ebx, ecx, edx);
        value.Leaf80000021Eax = eax;
        value.Leaf80000021Ebx = ebx;
        value.Leaf80000021Ecx = ecx;
        value.Leaf80000021Edx = edx;
        value.LfenceAlwaysSerializing =
            (eax & (UINT32_C(1) << 2u)) != 0u;
    }
#else
    value.Compiled = false;
#endif
    features = value;
    diagnostic.clear();
    return true;
}

#if WIREHAIR_RDPRU_COMPILED
__attribute__((noinline))
#endif
bool ReadRdpru(void* context, uint32_t selector, uint64_t* value_out)
{
    (void)context;
    if (!value_out || selector > 1u) return false;
#if WIREHAIR_RDPRU_COMPILED
    uint64_t value = 0u;
    __asm__ __volatile__("" ::: "memory");
    __asm__ __volatile__(
        ".globl wh2_rdpru_bracket_opcode_begin\n\t"
        "wh2_rdpru_bracket_opcode_begin:\n\t"
        "mfence\n\t"
        "lfence\n\t"
        ".byte 0x0f, 0x01, 0xfd\n\t"
        "shlq $32, %%rdx\n\t"
        "orq %%rdx, %%rax\n\t"
        "lfence\n\t"
        ".globl wh2_rdpru_bracket_opcode_end\n\t"
        "wh2_rdpru_bracket_opcode_end:\n\t"
        : "=&a"(value)
        : "c"(selector)
        : "rdx", "cc", "memory");
    __asm__ __volatile__("" ::: "memory");
    *value_out = value;
    return true;
#else
    return false;
#endif
}

#if defined(__linux__)
static sigjmp_buf gAccessJump;
static volatile sig_atomic_t gAccessProbeActive = 0;
static volatile sig_atomic_t gAccessSignal = 0;

void AccessSignalHandler(int signal_number)
{
    if (gAccessProbeActive != 0)
    {
        gAccessSignal = signal_number;
        siglongjmp(gAccessJump, 1);
    }
    _exit(128 + signal_number);
}
#endif

bool AccessProbeInactive()
{
#if defined(__linux__)
    return gAccessProbeActive == 0;
#else
    return true;
#endif
}

bool GuardedCounterRead(
    wirehair_wh2_bench::NativeCounterReadFunction read,
    void* context,
    uint32_t selector,
    uint64_t& value_out,
    int& signal_out)
{
    signal_out = 0;
    if (!read) return false;
#if defined(__linux__)
    struct sigaction action;
    struct sigaction old_ill;
    struct sigaction old_segv;
    struct sigaction old_bus;
    std::memset(&action, 0, sizeof(action));
    action.sa_handler = AccessSignalHandler;
    sigemptyset(&action.sa_mask);
    if (sigaction(SIGILL, &action, &old_ill) != 0) return false;
    if (sigaction(SIGSEGV, &action, &old_segv) != 0)
    {
        (void)sigaction(SIGILL, &old_ill, nullptr);
        return false;
    }
    if (sigaction(SIGBUS, &action, &old_bus) != 0)
    {
        (void)sigaction(SIGSEGV, &old_segv, nullptr);
        (void)sigaction(SIGILL, &old_ill, nullptr);
        return false;
    }
    bool ok = false;
    gAccessSignal = 0;
    gAccessProbeActive = 1;
    if (sigsetjmp(gAccessJump, 1) == 0)
    {
        try {
            ok = read(context, selector, &value_out);
        }
        catch (...) {
            ok = false;
        }
    }
    gAccessProbeActive = 0;
    signal_out = static_cast<int>(gAccessSignal);
    const bool restored_bus = sigaction(SIGBUS, &old_bus, nullptr) == 0;
    const bool restored_segv = sigaction(SIGSEGV, &old_segv, nullptr) == 0;
    const bool restored_ill = sigaction(SIGILL, &old_ill, nullptr) == 0;
    const bool restored = restored_bus && restored_segv && restored_ill;
    return ok && restored && signal_out == 0;
#else
    (void)context;
    (void)selector;
    (void)value_out;
    return false;
#endif
}

bool VerifyCpu(int target_cpu, std::string& diagnostic)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "target CPU is outside cpu_set_t range";
        return false;
    }
    cpu_set_t observed;
    CPU_ZERO(&observed);
    if (sched_getaffinity(0, sizeof(observed), &observed) != 0)
    {
        diagnostic = std::string("sched_getaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    uint32_t selected = 0u;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        selected += CPU_ISSET(cpu, &observed) ? 1u : 0u;
    }
    if (selected != 1u || !CPU_ISSET(target_cpu, &observed))
    {
        diagnostic = "affinity is not the requested singleton CPU";
        return false;
    }
    const int cpu = sched_getcpu();
    if (cpu < 0)
    {
        diagnostic = std::string("sched_getcpu failed: ") +
            std::strerror(errno);
        return false;
    }
    if (cpu != target_cpu)
    {
        std::ostringstream message;
        message << "CPU migration from " << target_cpu << " to " << cpu;
        diagnostic = message.str();
        return false;
    }
    diagnostic.clear();
    return true;
#else
    (void)target_cpu;
    diagnostic = "singleton CPU affinity requires Linux";
    return false;
#endif
}

bool PinCpuOnce(int target_cpu, std::string& diagnostic)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "target CPU is outside cpu_set_t range";
        return false;
    }
    cpu_set_t selected;
    CPU_ZERO(&selected);
    CPU_SET(target_cpu, &selected);
    if (sched_setaffinity(0, sizeof(selected), &selected) != 0)
    {
        diagnostic = std::string("sched_setaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    return VerifyCpu(target_cpu, diagnostic);
#else
    (void)target_cpu;
    diagnostic = "singleton CPU affinity requires Linux";
    return false;
#endif
}

bool ParseCpu(const char* text, int& cpu_out)
{
    if (!text || !*text) return false;
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value < 0 ||
        value > std::numeric_limits<int>::max())
    {
        return false;
    }
    cpu_out = static_cast<int>(value);
    return true;
}

#if defined(__linux__)
bool ReadSmallTextFile(
    const std::string& path,
    std::string& content,
    std::string& diagnostic)
{
#if defined(__linux__)
    content.clear();
    const int file = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (file < 0)
    {
        diagnostic = std::string("open topology file failed: ") +
            std::strerror(errno);
        return false;
    }
    static const size_t maximum = 4096u;
    char buffer[256];
    bool ok = true;
    while (content.size() <= maximum)
    {
        const ssize_t got = read(file, buffer, sizeof(buffer));
        if (got < 0)
        {
            if (errno == EINTR) continue;
            diagnostic = std::string("read topology file failed: ") +
                std::strerror(errno);
            ok = false;
            break;
        }
        if (got == 0) break;
        if (static_cast<size_t>(got) > maximum -
                std::min(content.size(), maximum))
        {
            diagnostic = "topology file exceeds bound";
            ok = false;
            break;
        }
        content.append(buffer, static_cast<size_t>(got));
    }
    const int close_result = close(file);
    if (ok && close_result != 0)
    {
        diagnostic = std::string("close topology file failed: ") +
            std::strerror(errno);
        ok = false;
    }
    if (ok && (content.empty() || content.size() > maximum))
    {
        diagnostic = "topology file is empty or exceeds bound";
        ok = false;
    }
    return ok;
#else
    (void)path;
    (void)content;
    diagnostic = "host topology requires Linux";
    return false;
#endif
}
#endif

bool ParseUnsignedToken(
    const std::string& text,
    size_t begin,
    size_t end,
    int& value_out)
{
    if (begin >= end) return false;
    if (end - begin > 1u && text[begin] == '0') return false;
    uint64_t value = 0u;
    for (size_t i = begin; i < end; ++i)
    {
        if (text[i] < '0' || text[i] > '9') return false;
        value = value * 10u + static_cast<unsigned>(text[i] - '0');
        if (value > UINT64_C(65535)) return false;
    }
    value_out = static_cast<int>(value);
    return true;
}

bool ParseCpuList(const std::string& input, std::vector<int>& cpus)
{
    cpus.clear();
    std::string text(input);
    if (!text.empty() && text.back() == '\n') text.resize(text.size() - 1u);
    if (text.empty() || text.find('\n') != std::string::npos ||
        text.find('\r') != std::string::npos) return false;
    size_t offset = 0u;
    int previous = -1;
    while (offset < text.size())
    {
        const size_t comma = text.find(',', offset);
        const size_t end = comma == std::string::npos ? text.size() : comma;
        const size_t dash = text.find('-', offset);
        int first = -1;
        int last = -1;
        if (dash != std::string::npos && dash < end)
        {
            if (text.find('-', dash + 1u) < end ||
                !ParseUnsignedToken(text, offset, dash, first) ||
                !ParseUnsignedToken(text, dash + 1u, end, last) ||
                first > last) return false;
        }
        else
        {
            if (!ParseUnsignedToken(text, offset, end, first)) return false;
            last = first;
        }
        if (first <= previous) return false;
        for (int cpu = first; cpu <= last; ++cpu)
        {
            cpus.push_back(cpu);
            previous = cpu;
        }
        if (comma == std::string::npos) break;
        offset = comma + 1u;
        if (offset == text.size()) return false;
    }
    return !cpus.empty();
}

std::string FormatCpuList(const std::vector<int>& cpus)
{
    if (cpus.empty()) return std::string();
    std::ostringstream text;
    size_t first = 0u;
    while (first < cpus.size())
    {
        size_t last = first;
        while (last + 1u < cpus.size() &&
               cpus[last + 1u] == cpus[last] + 1) ++last;
        if (first != 0u) text << ',';
        text << cpus[first];
        if (last != first) text << '-' << cpus[last];
        first = last + 1u;
    }
    return text.str();
}

#if defined(__linux__)
bool ParseTopologyInteger(const std::string& input, int& value_out)
{
    std::string text(input);
    if (!text.empty() && text.back() == '\n') text.resize(text.size() - 1u);
    return ParseUnsignedToken(text, 0u, text.size(), value_out);
}
#endif

struct HostTopology
{
    std::string OnlineRaw;
    std::string SiblingsRaw;
    std::vector<int> Online;
    std::vector<int> InitialAffinity;
    std::vector<int> TargetSiblings;
    uint32_t PhysicalCoreCount = 0u;
};

bool ReadHostTopology(HostTopology& topology, std::string& diagnostic)
{
#if defined(__linux__)
    HostTopology value;
    if (!ReadSmallTextFile(
            "/sys/devices/system/cpu/online", value.OnlineRaw, diagnostic))
    {
        topology = value;
        return false;
    }
    if (!ParseCpuList(value.OnlineRaw, value.Online))
    {
        diagnostic = "invalid online CPU list";
        topology = value;
        return false;
    }
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    if (sched_getaffinity(0, sizeof(affinity), &affinity) != 0)
    {
        diagnostic = std::string("read initial affinity failed: ") +
            std::strerror(errno);
        topology = value;
        return false;
    }
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        if (CPU_ISSET(cpu, &affinity)) value.InitialAffinity.push_back(cpu);
    }
    std::ostringstream siblings_path;
    siblings_path << "/sys/devices/system/cpu/cpu" << kFrozenTargetCpu
                  << "/topology/thread_siblings_list";
    if (!ReadSmallTextFile(
            siblings_path.str(), value.SiblingsRaw, diagnostic))
    {
        topology = value;
        return false;
    }
    if (!ParseCpuList(value.SiblingsRaw, value.TargetSiblings))
    {
        diagnostic = "invalid target sibling CPU list";
        topology = value;
        return false;
    }
    std::set<std::pair<int, int> > physical_cores;
    for (size_t i = 0u; i < value.Online.size(); ++i)
    {
        const int cpu = value.Online[i];
        std::ostringstream package_path;
        std::ostringstream core_path;
        package_path << "/sys/devices/system/cpu/cpu" << cpu
                     << "/topology/physical_package_id";
        core_path << "/sys/devices/system/cpu/cpu" << cpu
                  << "/topology/core_id";
        std::string package_text;
        std::string core_text;
        int package = -1;
        int core = -1;
        if (!ReadSmallTextFile(package_path.str(), package_text, diagnostic) ||
            !ReadSmallTextFile(core_path.str(), core_text, diagnostic) ||
            !ParseTopologyInteger(package_text, package) ||
            !ParseTopologyInteger(core_text, core))
        {
            if (diagnostic.empty()) diagnostic = "invalid topology integer";
            topology = value;
            return false;
        }
        physical_cores.insert(std::make_pair(package, core));
    }
    value.PhysicalCoreCount = static_cast<uint32_t>(physical_cores.size());
    topology = value;
    diagnostic.clear();
    return true;
#else
    (void)topology;
    diagnostic = "host topology requires Linux";
    return false;
#endif
}

bool ValidateHostTopology(
    const HostTopology& topology,
    std::string& diagnostic)
{
    if (topology.Online.size() != kFrozenLogicalCpuCount ||
        topology.InitialAffinity.size() != kFrozenLogicalCpuCount ||
        topology.PhysicalCoreCount != kFrozenPhysicalCoreCount)
    {
        diagnostic = "host logical/core/initial-affinity count mismatch";
        return false;
    }
    for (int cpu = kFrozenFirstOnlineCpu;
         cpu <= kFrozenLastOnlineCpu;
         ++cpu)
    {
        const size_t index = static_cast<size_t>(cpu - kFrozenFirstOnlineCpu);
        if (topology.Online[index] != cpu ||
            topology.InitialAffinity[index] != cpu)
        {
            diagnostic = "host online or initial-affinity CPU set mismatch";
            return false;
        }
    }
    if (topology.TargetSiblings.size() != 2u ||
        topology.TargetSiblings[0] != kFrozenTargetCpu ||
        topology.TargetSiblings[1] != kFrozenTargetSiblingCpu)
    {
        diagnostic = "target CPU thread-sibling set mismatch";
        return false;
    }
    diagnostic.clear();
    return true;
}

struct ContextReceipt
{
    int32_t Cpu = -1;
    uint32_t AffinityCount = 0u;
    uint32_t TargetInAffinity = 0u;
    int64_t Voluntary = -1;
    int64_t Involuntary = -1;
};

typedef bool (*ContextCaptureFunction)(
    int target_cpu,
    ContextReceipt& receipt,
    std::string& diagnostic);

bool CaptureContext(
    int target_cpu,
    ContextReceipt& receipt,
    std::string& diagnostic)
{
#if defined(__linux__)
    cpu_set_t observed;
    CPU_ZERO(&observed);
    if (sched_getaffinity(0, sizeof(observed), &observed) != 0)
    {
        diagnostic = std::string("sched_getaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    uint32_t selected = 0u;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        selected += CPU_ISSET(cpu, &observed) ? 1u : 0u;
    }
    const int current_cpu = sched_getcpu();
    if (current_cpu < 0)
    {
        diagnostic = std::string("sched_getcpu failed: ") +
            std::strerror(errno);
        return false;
    }
    struct rusage usage;
    if (getrusage(RUSAGE_THREAD, &usage) != 0)
    {
        diagnostic = std::string("getrusage(RUSAGE_THREAD) failed: ") +
            std::strerror(errno);
        return false;
    }
    receipt.Cpu = current_cpu;
    receipt.AffinityCount = selected;
    receipt.TargetInAffinity =
        target_cpu >= 0 && target_cpu < CPU_SETSIZE &&
        CPU_ISSET(target_cpu, &observed) ? 1u : 0u;
    receipt.Voluntary = static_cast<int64_t>(usage.ru_nvcsw);
    receipt.Involuntary = static_cast<int64_t>(usage.ru_nivcsw);
    diagnostic.clear();
    return true;
#else
    (void)target_cpu;
    (void)receipt;
    diagnostic = "context receipts require Linux";
    return false;
#endif
}

bool ContextPairValid(
    int target_cpu,
    const ContextReceipt& before,
    const ContextReceipt& after,
    std::string& diagnostic)
{
    if (before.Cpu != target_cpu || after.Cpu != target_cpu ||
        before.AffinityCount != 1u || after.AffinityCount != 1u ||
        before.TargetInAffinity != 1u || after.TargetInAffinity != 1u)
    {
        diagnostic = "CPU migration or affinity change inside bracket";
        return false;
    }
    if (before.Voluntary < 0 || before.Involuntary < 0 ||
        after.Voluntary < before.Voluntary ||
        after.Involuntary < before.Involuntary)
    {
        diagnostic = "invalid context-switch receipt";
        return false;
    }
    if (after.Voluntary != before.Voluntary ||
        after.Involuntary != before.Involuntary)
    {
        diagnostic = "context switch occurred inside bracket";
        return false;
    }
    diagnostic.clear();
    return true;
}

struct EmptyReceipt
{
    uint32_t Selector = 0u;
    uint32_t Sample = 0u;
    uint32_t SelectorPosition = 0u;
    uint64_t Start = 0u;
    uint64_t End = 0u;
    uint64_t Delta = 0u;
    uint32_t BeginReadSucceeded = 0u;
    uint32_t EndReadSucceeded = 0u;
    ContextReceipt Before;
    ContextReceipt After;
};

bool RunEmptyBracket(
    const NativeCounterBracket& bracket,
    int target_cpu,
    EmptyReceipt& receipt,
    std::string& diagnostic,
    ContextCaptureFunction capture_context = CaptureContext)
{
    if (!bracket.Read || !capture_context || bracket.Selector > 1u) {
        diagnostic = "invalid empty-bracket reader";
        return false;
    }
    receipt.Start = 0u;
    receipt.End = 0u;
    receipt.Delta = 0u;
    receipt.BeginReadSucceeded = 0u;
    receipt.EndReadSucceeded = 0u;
    receipt.Before = ContextReceipt();
    receipt.After = ContextReceipt();
    if (!capture_context(target_cpu, receipt.Before, diagnostic)) return false;
    bool begin_ok = false;
    bool end_ok = false;
    bool reader_threw = false;
    try
    {
        begin_ok = bracket.Read(
            bracket.Context, bracket.Selector, &receipt.Start);
        receipt.BeginReadSucceeded = begin_ok ? 1u : 0u;
        end_ok = begin_ok && bracket.Read(
            bracket.Context, bracket.Selector, &receipt.End);
        receipt.EndReadSucceeded = end_ok ? 1u : 0u;
    }
    catch (...) {
        reader_threw = true;
    }
    // The after-context receipt is mandatory even when either counter read
    // throws, so the failed attempt preserves the complete bracket boundary.
    if (!capture_context(target_cpu, receipt.After, diagnostic)) return false;
    if (!ContextPairValid(target_cpu, receipt.Before, receipt.After, diagnostic)) {
        return false;
    }
    if (reader_threw)
    {
        diagnostic = "empty-bracket reader threw";
        return false;
    }
    if (!begin_ok || !end_ok || receipt.End <= receipt.Start)
    {
        diagnostic = "empty bracket read failed, reset, wrapped, or was nonpositive";
        return false;
    }
    receipt.Delta = receipt.End - receipt.Start;
    return true;
}

struct FixtureBank
{
    std::vector<uint8_t> Source;
    std::vector<uint32_t> IidIds;
    std::vector<uint32_t> RepairIds;
    NativeArm Arm;
    NativeEncoderFixture Encoder;
    NativeReceiveFixture IidReceive;
    NativeReceiveFixture RepairReceive;
};

struct CalibrationCell
{
    uint32_t Index = 0u;
    uint32_t WidthIndex = 0u;
    uint32_t RootIndex = 0u;
    uint32_t BlockBytes = 0u;
    uint64_t Root = 0u;
    uint32_t IidOverhead = UINT32_MAX;
    std::array<FixtureBank, 3> Banks;
    std::string SourceSha256;
    std::string IidIdsSha256;
    std::string RepairIdsSha256;
};

typedef std::array<std::unique_ptr<CalibrationCell>, kBaseCellCount> CellArray;

struct LogicalCellInfo
{
    uint32_t Index = 0u;
    uint32_t BaseIndex = 0u;
    uint32_t WidthIndex = 0u;
    uint32_t RootIndex = 0u;
    Scope ScopeValue = Scope::Encoder;
    uint32_t BlockBytes = 0u;
    uint64_t Root = 0u;
};

LogicalCellInfo GetLogicalCellInfo(uint32_t cell)
{
    LogicalCellInfo info;
    info.Index = cell;
    info.WidthIndex = cell / 9u;
    const uint32_t within_width = cell % 9u;
    info.ScopeValue = static_cast<Scope>(within_width / 3u);
    info.RootIndex = within_width % 3u;
    info.BaseIndex = info.WidthIndex * 3u + info.RootIndex;
    info.BlockBytes = kBlockBytes[info.WidthIndex];
    info.Root = kRoots[info.RootIndex];
    return info;
}

const FixtureBank& SelectBank(
    const CalibrationCell& cell,
    PhysicalBank bank)
{
    return cell.Banks[static_cast<uint32_t>(bank)];
}

Batch128CounterArmResult RunFixtureBatch128(
    const FixtureBank& bank,
    Scope scope,
    const NativeCounterBracket* bracket)
{
    if (scope == Scope::Encoder) {
        return bracket ? bank.Encoder.RunWh1Batch128Counter(*bracket) :
            bank.Encoder.PreflightBatch128();
    }
    if (scope == Scope::IidReceive) {
        return bracket ? bank.IidReceive.RunWh1Batch128Counter(*bracket) :
            bank.IidReceive.PreflightBatch128();
    }
    return bracket ? bank.RepairReceive.RunWh1Batch128Counter(*bracket) :
        bank.RepairReceive.PreflightBatch128();
}

bool FullVerificationMask(const Batch128CounterArmResult& result)
{
    return result.VerificationMask[0] == UINT64_MAX &&
        result.VerificationMask[1] == UINT64_MAX;
}

bool ValidSemanticBatch128Result(
    Scope scope,
    uint32_t expected_iid_overhead,
    const Batch128CounterArmResult& result,
    bool measured)
{
    if (wirehair_wh2_bench::kNativeCounterBatchSize != kBatchSize) {
        return false;
    }
    if (measured)
    {
        if (!result.BeginReadSucceeded || !result.EndReadSucceeded ||
            result.CounterEnd <= result.CounterStart ||
            result.CounterDelta == 0u ||
            result.CounterDelta != result.CounterEnd - result.CounterStart)
        {
            return false;
        }
    }
    else if (result.BeginReadSucceeded || result.EndReadSucceeded ||
        result.CounterStart != 0u || result.CounterEnd != 0u ||
        result.CounterDelta != 0u)
    {
        return false;
    }
    if (!FullVerificationMask(result)) return false;
    for (size_t i = 0u; i < kBatchSize; ++i)
    {
        if (result.Results[i] != Wirehair_Success ||
            result.DirectSystematicPackets[i] != 0u)
        {
            return false;
        }
        if (scope == Scope::Encoder)
        {
            if (result.DecodedOverheads[i] != UINT32_MAX) return false;
        }
        else if (scope == Scope::IidReceive)
        {
            if (expected_iid_overhead == UINT32_MAX ||
                expected_iid_overhead > kReceiveOverheadCap ||
                result.DecodedOverheads[i] != expected_iid_overhead)
            {
                return false;
            }
        }
        else if (result.DecodedOverheads[i] != 0u) {
            return false;
        }
    }
    return true;
}

std::string InnerReceiptSha256(const Batch128CounterArmResult& result)
{
    std::string bytes;
    static const char marker[] = "wirehair-rdpru-inner-receipt-le-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    AppendLe64(bytes, result.CounterStart);
    AppendLe64(bytes, result.CounterEnd);
    AppendLe64(bytes, result.CounterDelta);
    AppendLe32(bytes, result.BeginReadSucceeded ? 1u : 0u);
    AppendLe32(bytes, result.EndReadSucceeded ? 1u : 0u);
    for (size_t i = 0u; i < kBatchSize; ++i) {
        AppendLe32(bytes, static_cast<uint32_t>(result.Results[i]));
    }
    AppendLe64(bytes, result.VerificationMask[0]);
    AppendLe64(bytes, result.VerificationMask[1]);
    for (size_t i = 0u; i < kBatchSize; ++i) {
        AppendLe32(bytes, result.DecodedOverheads[i]);
    }
    for (size_t i = 0u; i < kBatchSize; ++i) {
        AppendLe64(bytes, result.DirectSystematicPackets[i]);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

struct SemanticSummary
{
    uint32_t SuccessCount = 0u;
    uint64_t OverheadSum = 0u;
    uint32_t OverheadMinimum = UINT32_MAX;
    uint32_t OverheadMaximum = 0u;
    uint64_t DirectSystematicSum = 0u;
    std::string InnerSha256;
};

bool SummarizeResult(
    const Batch128CounterArmResult& result,
    SemanticSummary& summary)
{
    SemanticSummary value;
    for (size_t i = 0u; i < kBatchSize; ++i)
    {
        value.SuccessCount += result.Results[i] == Wirehair_Success ? 1u : 0u;
        if (!AddU64(result.DecodedOverheads[i], value.OverheadSum) ||
            !AddU64(result.DirectSystematicPackets[i],
                value.DirectSystematicSum))
        {
            return false;
        }
        value.OverheadMinimum = std::min(
            value.OverheadMinimum, result.DecodedOverheads[i]);
        value.OverheadMaximum = std::max(
            value.OverheadMaximum, result.DecodedOverheads[i]);
    }
    value.InnerSha256 = InnerReceiptSha256(result);
    if (!IsLowerHexSha256(value.InnerSha256.c_str())) return false;
    summary = value;
    return true;
}

bool InitializeBank(
    uint32_t block_bytes,
    uint64_t root,
    FixtureBank& bank,
    std::string& diagnostic)
{
    if (!wirehair_wh2_bench::MakeDeterministicSource(
            kK, block_bytes, root, bank.Source))
    {
        diagnostic = "deterministic source generation failed";
        return false;
    }
    wirehair::wh2_benchmark::FrozenPacketTrace iid_trace;
    if (wirehair::wh2_benchmark::GenerateFrozenPacketTrace(
            kK,
            block_bytes,
            kLossPpm,
            root,
            wirehair::wh2_benchmark::FrozenSchedule::Iid,
            kReceiveOverheadCap,
            iid_trace) !=
            wirehair::wh2_benchmark::FrozenTraceStatus::Complete ||
        iid_trace.delivered_ids.size() !=
            static_cast<size_t>(kK) + kReceiveOverheadCap)
    {
        diagnostic = "exact IID10 K+256 trace generation failed";
        return false;
    }
    bank.IidIds = iid_trace.delivered_ids;
    bank.RepairIds.resize(static_cast<size_t>(kK) + kReceiveOverheadCap);
    for (size_t i = 0u; i < bank.RepairIds.size(); ++i) {
        bank.RepairIds[i] = static_cast<uint32_t>(i) + 2u;
    }

    const NativeArmSpec spec = wirehair_wh2_bench::MakeWirehair1Arm();
    WirehairResult result = bank.Arm.Initialize(
        spec, kK, block_bytes, bank.Source);
    if (result != Wirehair_Success ||
        bank.Arm.Kind() != wirehair_wh2_bench::NativeArmKind::Wirehair1 ||
        bank.Arm.BlockCount() != kK || bank.Arm.BlockBytes() != block_bytes)
    {
        diagnostic = "Wirehair1 native arm initialization failed";
        return false;
    }
    result = bank.Encoder.Initialize(spec, kK, block_bytes, bank.Source);
    if (result != Wirehair_Success)
    {
        diagnostic = "Wirehair1 encoder fixture initialization failed";
        return false;
    }
    result = bank.IidReceive.Initialize(
        bank.Arm, bank.IidIds, kReceiveOverheadCap);
    if (result != Wirehair_Success)
    {
        diagnostic = "Wirehair1 IID fixture initialization failed";
        return false;
    }
    result = bank.RepairReceive.Initialize(
        bank.Arm, bank.RepairIds, kReceiveOverheadCap);
    if (result != Wirehair_Success)
    {
        diagnostic = "Wirehair1 repair fixture initialization failed";
        return false;
    }
    return true;
}

bool BuildCells(CellArray& cells, std::string& diagnostic)
{
    for (uint32_t cell_index = 0u; cell_index < kBaseCellCount; ++cell_index)
    {
        std::unique_ptr<CalibrationCell> cell(new CalibrationCell);
        cell->Index = cell_index;
        cell->WidthIndex = cell_index / 3u;
        cell->RootIndex = cell_index % 3u;
        cell->BlockBytes = kBlockBytes[cell->WidthIndex];
        cell->Root = kRoots[cell->RootIndex];
        for (uint32_t bank_index = 0u; bank_index < 3u; ++bank_index)
        {
            if (!InitializeBank(
                    cell->BlockBytes,
                    cell->Root,
                    cell->Banks[bank_index],
                    diagnostic))
            {
                return false;
            }
        }
        const FixtureBank& x = cell->Banks[0];
        const FixtureBank& y = cell->Banks[1];
        const FixtureBank& w = cell->Banks[2];
        if (x.Source != y.Source || x.Source != w.Source ||
            x.IidIds != y.IidIds || x.IidIds != w.IidIds ||
            x.RepairIds != y.RepairIds || x.RepairIds != w.RepairIds ||
            x.Source.empty() || x.IidIds.empty() || x.RepairIds.empty() ||
            x.Source.data() == y.Source.data() ||
            x.Source.data() == w.Source.data() ||
            y.Source.data() == w.Source.data() ||
            x.IidIds.data() == y.IidIds.data() ||
            x.IidIds.data() == w.IidIds.data() ||
            y.IidIds.data() == w.IidIds.data() ||
            x.RepairIds.data() == y.RepairIds.data() ||
            x.RepairIds.data() == w.RepairIds.data() ||
            y.RepairIds.data() == w.RepairIds.data())
        {
            diagnostic = "X/Y/W fixture ownership is not independent and equal";
            return false;
        }
        cell->SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
            x.Source.data(), x.Source.size());
        cell->IidIdsSha256 =
            wirehair::wh2_benchmark::PacketIdsSha256(x.IidIds);
        cell->RepairIdsSha256 =
            wirehair::wh2_benchmark::PacketIdsSha256(x.RepairIds);
        if (!IsLowerHexSha256(cell->SourceSha256.c_str()) ||
            !IsLowerHexSha256(cell->IidIdsSha256.c_str()) ||
            !IsLowerHexSha256(cell->RepairIdsSha256.c_str()))
        {
            diagnostic = "fixture provenance hashing failed";
            return false;
        }
        cells[cell_index] = std::move(cell);
    }
    return true;
}

bool PreflightCells(CellArray& cells, std::string& diagnostic)
{
    size_t count = 0u;
    for (uint32_t cell_index = 0u; cell_index < kBaseCellCount; ++cell_index)
    {
        CalibrationCell& cell = *cells[cell_index];
        const Batch128CounterArmResult x_iid =
            cell.Banks[0].IidReceive.PreflightBatch128();
        ++count;
        if (!FullVerificationMask(x_iid) ||
            x_iid.DecodedOverheads[0] == UINT32_MAX ||
            x_iid.DecodedOverheads[0] > kReceiveOverheadCap)
        {
            diagnostic = "X IID clock-free preflight failed";
            return false;
        }
        cell.IidOverhead = x_iid.DecodedOverheads[0];
        if (cell.IidOverhead != kExpectedIidOverheads[cell_index])
        {
            diagnostic = "IID overhead differs from frozen control";
            return false;
        }
        for (uint32_t bank_index = 0u; bank_index < 3u; ++bank_index)
        {
            const FixtureBank& bank = cell.Banks[bank_index];
            for (uint32_t scope_index = 0u; scope_index < 3u; ++scope_index)
            {
                const Scope scope = static_cast<Scope>(scope_index);
                Batch128CounterArmResult result;
                if (bank_index == 0u && scope == Scope::IidReceive) {
                    result = x_iid;
                }
                else
                {
                    result = RunFixtureBatch128(bank, scope, nullptr);
                    ++count;
                }
                if (!ValidSemanticBatch128Result(
                        scope, cell.IidOverhead, result, false))
                {
                    std::ostringstream message;
                    message << "clock-free preflight failed for base_cell="
                            << cell_index << ",bank=" << bank_index
                            << ",scope=" << ScopeName(scope);
                    diagnostic = message.str();
                    return false;
                }
            }
        }
    }
    if (count != kPreflightBatchCount)
    {
        diagnostic = "clock-free preflight batch count invalid";
        return false;
    }
    return true;
}

std::string FixtureSha256(const CellArray& cells)
{
    std::string bytes;
    static const char marker[] = "wirehair-rdpru-fixtures-le-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    for (size_t i = 0u; i < cells.size(); ++i)
    {
        if (!cells[i]) return std::string();
        const CalibrationCell& cell = *cells[i];
        AppendLe32(bytes, cell.Index);
        AppendLe32(bytes, cell.WidthIndex);
        AppendLe32(bytes, cell.RootIndex);
        AppendLe32(bytes, cell.BlockBytes);
        AppendLe64(bytes, cell.Root);
        AppendLe32(bytes, cell.IidOverhead);
        AppendSizedString(bytes, cell.SourceSha256);
        AppendSizedString(bytes, cell.IidIdsSha256);
        AppendSizedString(bytes, cell.RepairIdsSha256);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

std::string DeterministicFixtureSha256(uint32_t mutated_overhead_cell)
{
    std::string bytes;
    static const char marker[] = "wirehair-rdpru-fixtures-le-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    for (uint32_t cell = 0u; cell < kBaseCellCount; ++cell)
    {
        const uint32_t width_index = cell / 3u;
        const uint32_t root_index = cell % 3u;
        const uint32_t block_bytes = kBlockBytes[width_index];
        const uint64_t root = kRoots[root_index];
        std::vector<uint8_t> source;
        wirehair::wh2_benchmark::FrozenPacketTrace iid_trace;
        if (!wirehair_wh2_bench::MakeDeterministicSource(
                kK, block_bytes, root, source) ||
            wirehair::wh2_benchmark::GenerateFrozenPacketTrace(
                kK,
                block_bytes,
                kLossPpm,
                root,
                wirehair::wh2_benchmark::FrozenSchedule::Iid,
                kReceiveOverheadCap,
                iid_trace) !=
                wirehair::wh2_benchmark::FrozenTraceStatus::Complete ||
            iid_trace.delivered_ids.size() !=
                static_cast<size_t>(kK) + kReceiveOverheadCap)
        {
            return std::string();
        }
        std::vector<uint32_t> repair_ids(
            static_cast<size_t>(kK) + kReceiveOverheadCap);
        for (size_t i = 0u; i < repair_ids.size(); ++i) {
            repair_ids[i] = static_cast<uint32_t>(i) + 2u;
        }
        const std::string source_sha = wirehair::wh2_benchmark::Sha256Hex(
            source.data(), source.size());
        const std::string iid_sha =
            wirehair::wh2_benchmark::PacketIdsSha256(
                iid_trace.delivered_ids);
        const std::string repair_sha =
            wirehair::wh2_benchmark::PacketIdsSha256(repair_ids);
        if (!IsLowerHexSha256(source_sha.c_str()) ||
            !IsLowerHexSha256(iid_sha.c_str()) ||
            !IsLowerHexSha256(repair_sha.c_str())) return std::string();
        uint32_t overhead = kExpectedIidOverheads[cell];
        if (cell == mutated_overhead_cell) ++overhead;
        AppendLe32(bytes, cell);
        AppendLe32(bytes, width_index);
        AppendLe32(bytes, root_index);
        AppendLe32(bytes, block_bytes);
        AppendLe64(bytes, root);
        AppendLe32(bytes, overhead);
        AppendSizedString(bytes, source_sha);
        AppendSizedString(bytes, iid_sha);
        AppendSizedString(bytes, repair_sha);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

struct QualifiedResult
{
    ContextReceipt Before;
    ContextReceipt After;
    Batch128CounterArmResult Batch;
    SemanticSummary Summary;
};

struct FailedAttemptReceipt
{
    bool Present = false;
    bool EmptyAttempt = false;
    std::string Stage;
    std::string Identity;
    uint64_t Ordinal = UINT64_MAX;
    uint32_t Selector = UINT32_MAX;
    uint32_t Cell = UINT32_MAX;
    Scope ScopeValue = Scope::Encoder;
    PhysicalBank Bank = PhysicalBank::W;
    EmptyReceipt Empty;
    QualifiedResult Qualified;
};

void RecordQualifiedFailure(
    FailedAttemptReceipt& failure,
    const char* stage,
    const std::string& identity,
    uint64_t ordinal,
    uint32_t selector,
    uint32_t cell,
    Scope scope,
    PhysicalBank bank,
    const QualifiedResult& result)
{
    failure = FailedAttemptReceipt();
    failure.Present = true;
    failure.Stage = stage;
    failure.Identity = identity;
    failure.Ordinal = ordinal;
    failure.Selector = selector;
    failure.Cell = cell;
    failure.ScopeValue = scope;
    failure.Bank = bank;
    failure.Qualified = result;
}

void RecordEmptyFailure(
    FailedAttemptReceipt& failure,
    const std::string& identity,
    uint64_t ordinal,
    const EmptyReceipt& receipt)
{
    failure = FailedAttemptReceipt();
    failure.Present = true;
    failure.EmptyAttempt = true;
    failure.Stage = "empty";
    failure.Identity = identity;
    failure.Ordinal = ordinal;
    failure.Selector = receipt.Selector;
    failure.Empty = receipt;
}

void RecordStageFailure(
    FailedAttemptReceipt& failure,
    const char* stage,
    const std::string& identity,
    uint64_t ordinal)
{
    failure = FailedAttemptReceipt();
    failure.Present = true;
    failure.Stage = stage;
    failure.Identity = identity;
    failure.Ordinal = ordinal;
}

enum class CompactKind : uint32_t
{
    Warm = 0u,
    Washout = 1u
};

struct CompactHardwareReceipt
{
    CompactKind Kind = CompactKind::Warm;
    uint64_t Ordinal = 0u;
    uint64_t GlobalPanelAction = UINT64_MAX;
    uint64_t SelectorPanelAction = UINT64_MAX;
    uint64_t OriginalSuperblock = UINT64_MAX;
    uint64_t SelectorSuperblock = UINT64_MAX;
    uint32_t Selector = 0u;
    uint32_t Cell = 0u;
    Scope ScopeValue = Scope::Encoder;
    PhysicalBank Bank = PhysicalBank::W;
    ContextReceipt Before;
    ContextReceipt After;
    uint64_t CounterStart = 0u;
    uint64_t CounterEnd = 0u;
    uint64_t CounterDelta = 0u;
    uint32_t BeginReadSucceeded = 0u;
    uint32_t EndReadSucceeded = 0u;
    uint64_t VerificationMask0 = 0u;
    uint64_t VerificationMask1 = 0u;
    SemanticSummary Summary;
};

CompactHardwareReceipt MakeCompactReceipt(
    CompactKind kind,
    uint64_t ordinal,
    uint32_t selector,
    uint32_t cell,
    Scope scope,
    PhysicalBank bank,
    const QualifiedResult& result)
{
    CompactHardwareReceipt receipt;
    receipt.Kind = kind;
    receipt.Ordinal = ordinal;
    receipt.Selector = selector;
    receipt.Cell = cell;
    receipt.ScopeValue = scope;
    receipt.Bank = bank;
    receipt.Before = result.Before;
    receipt.After = result.After;
    receipt.CounterStart = result.Batch.CounterStart;
    receipt.CounterEnd = result.Batch.CounterEnd;
    receipt.CounterDelta = result.Batch.CounterDelta;
    receipt.BeginReadSucceeded = result.Batch.BeginReadSucceeded ? 1u : 0u;
    receipt.EndReadSucceeded = result.Batch.EndReadSucceeded ? 1u : 0u;
    receipt.VerificationMask0 = result.Batch.VerificationMask[0];
    receipt.VerificationMask1 = result.Batch.VerificationMask[1];
    receipt.Summary = result.Summary;
    return receipt;
}

bool ExecuteQualifiedBatch128(
    const CalibrationCell& base_cell,
    const LogicalCellInfo& logical_cell,
    PhysicalBank bank_name,
    const NativeCounterBracket& bracket,
    int target_cpu,
    QualifiedResult& result_out,
    std::string& diagnostic)
{
    QualifiedResult result;
    if (!CaptureContext(target_cpu, result.Before, diagnostic))
    {
        result_out = result;
        return false;
    }
    uint32_t exception_kind = 0u;
    try {
        result.Batch = RunFixtureBatch128(
            SelectBank(base_cell, bank_name),
            logical_cell.ScopeValue,
            &bracket);
    }
    catch (const std::exception&) {
        exception_kind = 1u;
    }
    catch (...) {
        exception_kind = 2u;
    }
    if (!CaptureContext(target_cpu, result.After, diagnostic))
    {
        result_out = result;
        return false;
    }
    if (!ContextPairValid(target_cpu, result.Before, result.After, diagnostic)) {
        result_out = result;
        return false;
    }
    if (exception_kind != 0u)
    {
        diagnostic = exception_kind == 1u ?
            "batch invocation threw a standard exception" :
            "batch invocation threw a non-standard exception";
        result_out = result;
        return false;
    }
    if (!ValidSemanticBatch128Result(
            logical_cell.ScopeValue,
            base_cell.IidOverhead,
            result.Batch,
            true) ||
        !SummarizeResult(result.Batch, result.Summary))
    {
        std::ostringstream message;
        message << "measured batch semantic failure cell="
                << logical_cell.Index << ",bank=" << BankName(bank_name)
                << ",scope=" << ScopeName(logical_cell.ScopeValue)
                << ",selector=" << bracket.Selector;
        diagnostic = message.str();
        result_out = result;
        return false;
    }
    result_out = result;
    return true;
}

typedef bool (*QualifiedBatchRunFunction)(
    const CalibrationCell& base_cell,
    const LogicalCellInfo& logical_cell,
    PhysicalBank bank_name,
    const NativeCounterBracket& bracket,
    int target_cpu,
    QualifiedResult& result_out,
    std::string& diagnostic);

struct CalibrationRunFunctions
{
    void* CounterContext;
    NativeCounterReadFunction CounterRead;
    ContextCaptureFunction Capture;
    QualifiedBatchRunFunction RunQualified;
};

const CalibrationRunFunctions& ProductionRunFunctions()
{
    static const CalibrationRunFunctions functions = {
        nullptr,
        ReadRdpru,
        CaptureContext,
        ExecuteQualifiedBatch128
    };
    return functions;
}

bool ValidRunFunctions(const CalibrationRunFunctions& functions)
{
    return functions.CounterRead && functions.Capture &&
        functions.RunQualified;
}

struct RawTiming
{
    uint64_t RawIndex = 0u;
    uint64_t GlobalPanelAction = 0u;
    uint64_t SelectorPanelAction = 0u;
    uint64_t OriginalSuperblock = 0u;
    uint64_t SelectorSuperblock = 0u;
    uint32_t Cell = 0u;
    uint32_t Round = 0u;
    uint32_t Visit = 0u;
    uint32_t U = 0u;
    Direction DirectionValue = Direction::Forward;
    Mode ModeValue = Mode::Exposed;
    uint32_t ModePosition = 0u;
    uint32_t OrientationFlip = 0u;
    uint32_t Rotation = 0u;
    uint32_t FirstSelector = 0u;
    uint32_t SelectorPosition = 0u;
    uint32_t Selector = 0u;
    uint32_t SequencePosition = 0u;
    uint32_t Sequence = 0u;
    uint32_t SequenceFamily = 0u;
    uint32_t Period = 0u;
    LogicalLabel Label = LogicalLabel::A;
    PhysicalBank LogicalAPhysical = PhysicalBank::X;
    PhysicalBank Physical = PhysicalBank::X;
    PredecessorClass Predecessor = PredecessorClass::Washout;
    QualifiedResult Result;
};

struct SuperblockObservation
{
    SuperblockPlan Plan;
    uint32_t SelectorPosition = 0u;
    uint32_t Selector = 0u;
    std::array<uint64_t, 2> LogicalSums = {{ 0u, 0u }};
    std::array<uint64_t, 2> PhysicalSums = {{ 0u, 0u }};
    std::array<uint64_t, 4> FactorialSums = {{ 0u, 0u, 0u, 0u }};
    std::array<uint32_t, 2> LogicalCounts = {{ 0u, 0u }};
    std::array<uint32_t, 2> PhysicalCounts = {{ 0u, 0u }};
    std::array<uint32_t, 4> FactorialCounts = {{ 0u, 0u, 0u, 0u }};
    double LogLogicalAa = 0.0;
    double LogPhysicalXy = 0.0;
    double LogInteraction = 0.0;
};

bool RatioLog(uint64_t numerator, uint64_t denominator, double& value_out)
{
    if (numerator == 0u || denominator == 0u) return false;
    const long double value = std::log(
        static_cast<long double>(numerator) /
        static_cast<long double>(denominator));
    if (!std::isfinite(value)) return false;
    value_out = static_cast<double>(value);
    return true;
}

bool AccumulateMeasured(
    SuperblockObservation& observation,
    LogicalLabel label,
    PhysicalBank physical,
    uint64_t delta)
{
    const uint32_t logical_index = static_cast<uint32_t>(label);
    const uint32_t physical_index = static_cast<uint32_t>(physical);
    if (logical_index >= 2u || physical_index >= 2u || delta == 0u) {
        return false;
    }
    const uint32_t factorial = logical_index * 2u + physical_index;
    return AddU64(delta, observation.LogicalSums[logical_index]) &&
        AddU64(delta, observation.PhysicalSums[physical_index]) &&
        AddU64(delta, observation.FactorialSums[factorial]) &&
        ++observation.LogicalCounts[logical_index] != 0u &&
        ++observation.PhysicalCounts[physical_index] != 0u &&
        ++observation.FactorialCounts[factorial] != 0u;
}

bool FinalizeObservation(SuperblockObservation& observation)
{
    for (size_t i = 0u; i < 2u; ++i) {
        if (observation.LogicalCounts[i] != 8u ||
            observation.PhysicalCounts[i] != 8u) return false;
    }
    for (size_t i = 0u; i < 4u; ++i) {
        if (observation.FactorialCounts[i] != 4u) return false;
    }
    double left = 0.0;
    double right = 0.0;
    if (!RatioLog(
            observation.LogicalSums[0],
            observation.LogicalSums[1],
            observation.LogLogicalAa) ||
        !RatioLog(
            observation.PhysicalSums[0],
            observation.PhysicalSums[1],
            observation.LogPhysicalXy) ||
        !RatioLog(
            observation.FactorialSums[0],
            observation.FactorialSums[1],
            left) ||
        !RatioLog(
            observation.FactorialSums[2],
            observation.FactorialSums[3],
            right))
    {
        return false;
    }
    observation.LogInteraction = left - right;
    return std::isfinite(observation.LogInteraction);
}

size_t ObservationKey(
    uint32_t selector,
    uint32_t cell,
    uint32_t round,
    Mode mode)
{
    return (((static_cast<size_t>(selector) * kCellCount + cell) *
        kRoundCount + round) * kModeCount +
        static_cast<uint32_t>(mode));
}

bool RunPinnedWarmups(
    const CellArray& cells,
    int target_cpu,
    std::vector<CompactHardwareReceipt>& warm,
    FailedAttemptReceipt& failure,
    std::string& diagnostic,
    const CalibrationRunFunctions& functions = ProductionRunFunctions())
{
    warm.clear();
    if (!ValidRunFunctions(functions))
    {
        diagnostic = "invalid warmup run functions";
        return false;
    }
    warm.reserve(kWarmBatchCount);
    for (uint32_t cell = 0u; cell < kCellCount; ++cell)
    {
        const LogicalCellInfo info = GetLogicalCellInfo(cell);
        const CalibrationCell& base = *cells[info.BaseIndex];
        for (uint32_t bank = 0u; bank < 3u; ++bank)
        {
            const uint32_t first_selector = (cell ^ bank) & 1u;
            for (uint32_t position = 0u; position < kSelectorCount; ++position)
            {
                NativeCounterBracket bracket;
                bracket.Context = functions.CounterContext;
                bracket.Read = functions.CounterRead;
                bracket.Selector = first_selector ^ position;
                QualifiedResult discarded;
                if (!functions.RunQualified(
                        base,
                        info,
                        static_cast<PhysicalBank>(bank),
                        bracket,
                        target_cpu,
                        discarded,
                        diagnostic))
                {
                    std::ostringstream identity;
                    identity << "warm;cell=" << cell
                             << ";bank=" << bank
                             << ";selector=" << bracket.Selector;
                    RecordQualifiedFailure(
                        failure,
                        "warm",
                        identity.str(),
                        static_cast<uint64_t>(warm.size()),
                        bracket.Selector,
                        cell,
                        info.ScopeValue,
                        static_cast<PhysicalBank>(bank),
                        discarded);
                    std::ostringstream message;
                    message << "warm failure ordinal=" << warm.size()
                            << ",cell=" << cell << ",bank=" << bank
                            << ",selector=" << bracket.Selector
                            << ",detail=" << diagnostic;
                    diagnostic = message.str();
                    return false;
                }
                warm.push_back(MakeCompactReceipt(
                    CompactKind::Warm,
                    static_cast<uint64_t>(warm.size()),
                    bracket.Selector,
                    cell,
                    info.ScopeValue,
                    static_cast<PhysicalBank>(bank),
                    discarded));
            }
        }
    }
    if (warm.size() != kWarmBatchCount)
    {
        diagnostic = "pinned warm batch count invalid";
        return false;
    }
    return true;
}

bool RunEmptySamples(
    int target_cpu,
    std::vector<EmptyReceipt>& empty,
    FailedAttemptReceipt& failure,
    std::string& diagnostic,
    const CalibrationRunFunctions& functions = ProductionRunFunctions())
{
    empty.clear();
    if (!ValidRunFunctions(functions))
    {
        diagnostic = "invalid empty-sample run functions";
        return false;
    }
    empty.reserve(static_cast<size_t>(kEmptySampleCount) * kSelectorCount);
    for (uint32_t sample = 0u; sample < kEmptySampleCount; ++sample)
    {
        const uint32_t first_selector = sample & 1u;
        for (uint32_t position = 0u; position < kSelectorCount; ++position)
        {
            EmptyReceipt receipt;
            receipt.Sample = sample;
            receipt.SelectorPosition = position;
            receipt.Selector = first_selector ^ position;
            NativeCounterBracket bracket;
            bracket.Context = functions.CounterContext;
            bracket.Read = functions.CounterRead;
            bracket.Selector = receipt.Selector;
            if (!RunEmptyBracket(
                    bracket,
                    target_cpu,
                    receipt,
                    diagnostic,
                    functions.Capture))
            {
                std::ostringstream identity;
                identity << "empty;sample=" << sample
                         << ";selector_position=" << position
                         << ";selector=" << receipt.Selector;
                RecordEmptyFailure(
                    failure,
                    identity.str(),
                    static_cast<uint64_t>(empty.size()),
                    receipt);
                std::ostringstream message;
                message << "empty failure sample=" << sample
                        << ",selector_position=" << position
                        << ",selector=" << receipt.Selector
                        << ",detail=" << diagnostic;
                diagnostic = message.str();
                return false;
            }
            empty.push_back(receipt);
        }
    }
    if (empty.size() != static_cast<size_t>(kEmptySampleCount) *
            kSelectorCount)
    {
        diagnostic = "empty sample count invalid";
        return false;
    }
    return true;
}

bool RunSelectorSubpanel(
    const SuperblockPlan& plan,
    uint64_t original_superblock,
    uint32_t selector_position,
    uint32_t selector,
    uint64_t selector_superblock,
    const CalibrationCell& base_cell,
    const LogicalCellInfo& logical_cell,
    int target_cpu,
    uint64_t& global_action,
    std::array<uint64_t, 2>& selector_action,
    std::array<size_t, 2>& measured_count,
    std::array<size_t, 2>& washout_count,
    std::vector<CompactHardwareReceipt>& washout,
    std::vector<RawTiming>& raw,
    SuperblockObservation& observation,
    FailedAttemptReceipt& failure,
    std::string& diagnostic,
    const CalibrationRunFunctions& functions)
{
    observation = SuperblockObservation();
    observation.Plan = plan;
    observation.SelectorPosition = selector_position;
    observation.Selector = selector;
    const size_t raw_before = raw.size();
    const size_t measured_before = measured_count[selector];
    const size_t washout_before = washout_count[selector];
    const uint64_t action_before = global_action;

    NativeCounterBracket bracket;
    bracket.Context = functions.CounterContext;
    bracket.Read = functions.CounterRead;
    bracket.Selector = selector;
    for (uint32_t sequence_position = 0u;
         sequence_position < kSequenceCount;
         ++sequence_position)
    {
        const uint32_t sequence =
            (plan.Rotation + sequence_position) % kSequenceCount;
        const PhysicalBank logical_a =
            LogicalAPhysicalForSequence(plan, sequence);
        const PhysicalBank logical_b = logical_a == PhysicalBank::X ?
            PhysicalBank::Y : PhysicalBank::X;
        for (uint32_t step = 0u; step < kBatchesPerSequence; ++step)
        {
            const bool measured = plan.ModeValue == Mode::Exposed ?
                step >= kPeriodsPerSequence : (step & 1u) != 0u;
            const uint32_t period = plan.ModeValue == Mode::Exposed ?
                (measured ? step - kPeriodsPerSequence : step) :
                step / 2u;
            const uint64_t this_global_action = global_action++;
            const uint64_t this_selector_action = selector_action[selector]++;
            if (!measured)
            {
                QualifiedResult discarded;
                if (!functions.RunQualified(
                        base_cell,
                        logical_cell,
                        PhysicalBank::W,
                        bracket,
                        target_cpu,
                        discarded,
                        diagnostic))
                {
                    std::ostringstream identity;
                    identity << "washout;global_action="
                             << this_global_action
                             << ";selector_action=" << this_selector_action
                             << ";original_superblock=" << original_superblock
                             << ";selector_superblock=" << selector_superblock
                             << ";cell=" << plan.Cell
                             << ";selector=" << selector;
                    RecordQualifiedFailure(
                        failure,
                        "washout",
                        identity.str(),
                        static_cast<uint64_t>(washout.size()),
                        selector,
                        plan.Cell,
                        logical_cell.ScopeValue,
                        PhysicalBank::W,
                        discarded);
                    std::ostringstream message;
                    message << "washout failure global_action="
                            << this_global_action
                            << ",selector_action=" << this_selector_action
                            << ",original_superblock=" << original_superblock
                            << ",selector_superblock=" << selector_superblock
                            << ",cell=" << plan.Cell
                            << ",selector=" << selector
                            << ",detail=" << diagnostic;
                    diagnostic = message.str();
                    return false;
                }
                CompactHardwareReceipt compact = MakeCompactReceipt(
                    CompactKind::Washout,
                    static_cast<uint64_t>(washout.size()),
                    selector,
                    plan.Cell,
                    logical_cell.ScopeValue,
                    PhysicalBank::W,
                    discarded);
                compact.GlobalPanelAction = this_global_action;
                compact.SelectorPanelAction = this_selector_action;
                compact.OriginalSuperblock = original_superblock;
                compact.SelectorSuperblock = selector_superblock;
                washout.push_back(compact);
                ++washout_count[selector];
                continue;
            }
            const LogicalLabel label = kSequences[sequence][period];
            const PhysicalBank physical = label == LogicalLabel::A ?
                logical_a : logical_b;
            RawTiming receipt;
            receipt.RawIndex = static_cast<uint64_t>(raw.size());
            receipt.GlobalPanelAction = this_global_action;
            receipt.SelectorPanelAction = this_selector_action;
            receipt.OriginalSuperblock = original_superblock;
            receipt.SelectorSuperblock = selector_superblock;
            receipt.Cell = plan.Cell;
            receipt.Round = plan.Round;
            receipt.Visit = plan.Visit;
            receipt.U = plan.U;
            receipt.DirectionValue = plan.DirectionValue;
            receipt.ModeValue = plan.ModeValue;
            receipt.ModePosition = plan.ModePosition;
            receipt.OrientationFlip = plan.OrientationFlip;
            receipt.Rotation = plan.Rotation;
            receipt.FirstSelector = plan.FirstSelector;
            receipt.SelectorPosition = selector_position;
            receipt.Selector = selector;
            receipt.SequencePosition = sequence_position;
            receipt.Sequence = sequence;
            receipt.SequenceFamily = sequence < 2u ? 0u : 1u;
            receipt.Period = period;
            receipt.Label = label;
            receipt.LogicalAPhysical = logical_a;
            receipt.Physical = physical;
            receipt.Predecessor = PredecessorClass::Washout;
            if (plan.ModeValue == Mode::Exposed && period > 0u)
            {
                const uint32_t previous = static_cast<uint32_t>(
                    kSequences[sequence][period - 1u]);
                const uint32_t current = static_cast<uint32_t>(label);
                receipt.Predecessor = static_cast<PredecessorClass>(
                    1u + previous * 2u + current);
            }
            if (!functions.RunQualified(
                    base_cell,
                    logical_cell,
                    physical,
                    bracket,
                    target_cpu,
                    receipt.Result,
                    diagnostic) ||
                !AccumulateMeasured(
                    observation,
                    label,
                    physical,
                    receipt.Result.Batch.CounterDelta))
            {
                std::ostringstream identity;
                identity << "measured;global_action="
                         << this_global_action
                         << ";selector_action=" << this_selector_action
                         << ";original_superblock=" << original_superblock
                         << ";selector_superblock=" << selector_superblock
                         << ";cell=" << plan.Cell
                         << ";selector=" << selector
                         << ";sequence_position=" << sequence_position
                         << ";period=" << period;
                RecordQualifiedFailure(
                    failure,
                    "measured",
                    identity.str(),
                    static_cast<uint64_t>(raw.size()),
                    selector,
                    plan.Cell,
                    logical_cell.ScopeValue,
                    physical,
                    receipt.Result);
                std::ostringstream message;
                message << "measured failure global_action="
                        << this_global_action
                        << ",selector_action=" << this_selector_action
                        << ",original_superblock=" << original_superblock
                        << ",selector_superblock=" << selector_superblock
                        << ",cell=" << plan.Cell
                        << ",selector=" << selector
                        << ",sequence_position=" << sequence_position
                        << ",period=" << period
                        << ",detail=" << diagnostic;
                diagnostic = message.str();
                return false;
            }
            raw.push_back(receipt);
            ++measured_count[selector];
        }
    }
    if (raw.size() - raw_before != kMeasuredBatchesPerSelectorSuperblock ||
        measured_count[selector] - measured_before !=
            kMeasuredBatchesPerSelectorSuperblock ||
        washout_count[selector] - washout_before !=
            kWashoutBatchesPerSelectorSuperblock ||
        global_action - action_before !=
            kMeasuredBatchesPerSelectorSuperblock +
                kWashoutBatchesPerSelectorSuperblock ||
        !FinalizeObservation(observation))
    {
        RecordStageFailure(
            failure,
            "panel-superblock-validation",
            "selector-superblock-count-or-arithmetic",
            original_superblock * kSelectorCount + selector_position);
        diagnostic = "selector superblock count or arithmetic invalid";
        return false;
    }
    return true;
}

bool RunFrozenPanel(
    const CellArray& cells,
    const FrozenSchedule& schedule,
    int target_cpu,
    std::vector<RawTiming>& raw,
    std::vector<CompactHardwareReceipt>& washout,
    std::vector<SuperblockObservation>& observations,
    FailedAttemptReceipt& failure,
    std::string& diagnostic,
    const CalibrationRunFunctions& functions = ProductionRunFunctions())
{
    raw.clear();
    washout.clear();
    observations.clear();
    if (!ValidRunFunctions(functions))
    {
        diagnostic = "invalid panel run functions";
        return false;
    }
    raw.reserve(kMeasuredBatchCount);
    washout.reserve(kWashoutBatchCount);
    observations.reserve(kSelectorSuperblockCount);
    uint64_t global_action = 0u;
    std::array<uint64_t, 2> selector_action = {{ 0u, 0u }};
    std::array<size_t, 2> measured_count = {{ 0u, 0u }};
    std::array<size_t, 2> washout_count = {{ 0u, 0u }};
    std::array<uint64_t, 2> selector_superblock = {{ 0u, 0u }};
    for (size_t plan_index = 0u; plan_index < schedule.size(); ++plan_index)
    {
        const SuperblockPlan& plan = schedule[plan_index];
        if (plan.Cell >= kCellCount)
        {
            RecordStageFailure(
                failure,
                "panel-schedule",
                "schedule-invalid-cell",
                static_cast<uint64_t>(plan_index));
            diagnostic = "schedule references invalid cell";
            return false;
        }
        const LogicalCellInfo info = GetLogicalCellInfo(plan.Cell);
        const CalibrationCell& base = *cells[info.BaseIndex];
        for (uint32_t selector_position = 0u;
             selector_position < kSelectorCount;
             ++selector_position)
        {
            const uint32_t selector = plan.FirstSelector ^ selector_position;
            SuperblockObservation observation;
            if (!RunSelectorSubpanel(
                    plan,
                    static_cast<uint64_t>(plan_index),
                    selector_position,
                    selector,
                    selector_superblock[selector]++,
                    base,
                    info,
                    target_cpu,
                    global_action,
                    selector_action,
                    measured_count,
                    washout_count,
                    washout,
                    raw,
                    observation,
                    failure,
                    diagnostic,
                    functions))
            {
                return false;
            }
            observations.push_back(observation);
        }
    }
    if (raw.size() != kMeasuredBatchCount ||
        washout.size() != kWashoutBatchCount ||
        observations.size() != kSelectorSuperblockCount ||
        global_action != kPanelActionCount)
    {
        RecordStageFailure(
            failure,
            "panel-cardinality",
            "final-panel-cardinality",
            static_cast<uint64_t>(observations.size()));
        diagnostic = "final panel cardinality invalid";
        return false;
    }
    for (uint32_t selector = 0u; selector < kSelectorCount; ++selector)
    {
        if (selector_action[selector] !=
                kMeasuredBatchCountPerSelector + kWashoutBatchCountPerSelector ||
            measured_count[selector] != kMeasuredBatchCountPerSelector ||
            washout_count[selector] != kWashoutBatchCountPerSelector ||
            selector_superblock[selector] != kOriginalSuperblockCount)
        {
            std::ostringstream identity;
            identity << "selector=" << selector;
            RecordStageFailure(
                failure,
                "panel-selector-cardinality",
                identity.str(),
                selector);
            diagnostic = "per-selector panel cardinality invalid";
            return false;
        }
    }
    return true;
}

bool EqualObservation(
    const SuperblockObservation& left,
    const SuperblockObservation& right)
{
    return left.Plan.Cell == right.Plan.Cell &&
        left.Plan.Round == right.Plan.Round &&
        left.Plan.Visit == right.Plan.Visit &&
        left.Plan.ModePosition == right.Plan.ModePosition &&
        left.Plan.ModeValue == right.Plan.ModeValue &&
        left.Plan.U == right.Plan.U &&
        left.Plan.DirectionValue == right.Plan.DirectionValue &&
        left.Plan.OrientationFlip == right.Plan.OrientationFlip &&
        left.Plan.Rotation == right.Plan.Rotation &&
        left.Plan.FirstSelector == right.Plan.FirstSelector &&
        left.SelectorPosition == right.SelectorPosition &&
        left.Selector == right.Selector &&
        left.LogicalSums == right.LogicalSums &&
        left.PhysicalSums == right.PhysicalSums &&
        left.FactorialSums == right.FactorialSums &&
        left.LogicalCounts == right.LogicalCounts &&
        left.PhysicalCounts == right.PhysicalCounts &&
        left.FactorialCounts == right.FactorialCounts &&
        left.LogLogicalAa == right.LogLogicalAa &&
        left.LogPhysicalXy == right.LogPhysicalXy &&
        left.LogInteraction == right.LogInteraction;
}

bool ValidateRawReceipts(
    const FrozenSchedule& schedule,
    const std::vector<RawTiming>& raw,
    const std::vector<SuperblockObservation>& observations,
    int target_cpu,
    std::string& diagnostic)
{
    if (raw.size() != kMeasuredBatchCount ||
        observations.size() != kSelectorSuperblockCount)
    {
        diagnostic = "post-panel cardinality invalid";
        return false;
    }
    size_t raw_next = 0u;
    size_t observation_next = 0u;
    uint64_t global_action = 0u;
    std::array<uint64_t, 2> selector_action = {{ 0u, 0u }};
    std::array<uint64_t, 2> selector_superblock = {{ 0u, 0u }};
    for (size_t plan_index = 0u; plan_index < schedule.size(); ++plan_index)
    {
        const SuperblockPlan& plan = schedule[plan_index];
        const LogicalCellInfo info = GetLogicalCellInfo(plan.Cell);
        for (uint32_t selector_position = 0u;
             selector_position < kSelectorCount;
             ++selector_position)
        {
            const uint32_t selector = plan.FirstSelector ^ selector_position;
            SuperblockObservation rebuilt;
            rebuilt.Plan = plan;
            rebuilt.SelectorPosition = selector_position;
            rebuilt.Selector = selector;
            const uint64_t expected_selector_superblock =
                selector_superblock[selector]++;
            for (uint32_t sequence_position = 0u;
                 sequence_position < kSequenceCount;
                 ++sequence_position)
            {
                const uint32_t sequence =
                    (plan.Rotation + sequence_position) % kSequenceCount;
                const PhysicalBank logical_a =
                    LogicalAPhysicalForSequence(plan, sequence);
                const PhysicalBank logical_b = logical_a == PhysicalBank::X ?
                    PhysicalBank::Y : PhysicalBank::X;
                for (uint32_t step = 0u; step < kBatchesPerSequence; ++step)
                {
                    const bool measured = plan.ModeValue == Mode::Exposed ?
                        step >= kPeriodsPerSequence : (step & 1u) != 0u;
                    const uint32_t period = plan.ModeValue == Mode::Exposed ?
                        (measured ? step - kPeriodsPerSequence : step) :
                        step / 2u;
                    const uint64_t expected_global_action = global_action++;
                    const uint64_t expected_selector_action =
                        selector_action[selector]++;
                    if (!measured) continue;
                    if (raw_next >= raw.size())
                    {
                        diagnostic = "post-panel raw receipt underflow";
                        return false;
                    }
                    const RawTiming& receipt = raw[raw_next];
                    const LogicalLabel label = kSequences[sequence][period];
                    const PhysicalBank physical = label == LogicalLabel::A ?
                        logical_a : logical_b;
                    PredecessorClass predecessor = PredecessorClass::Washout;
                    if (plan.ModeValue == Mode::Exposed && period > 0u)
                    {
                        const uint32_t previous = static_cast<uint32_t>(
                            kSequences[sequence][period - 1u]);
                        predecessor = static_cast<PredecessorClass>(
                            1u + previous * 2u +
                                static_cast<uint32_t>(label));
                    }
                    SemanticSummary summary;
                    if (receipt.RawIndex != raw_next ||
                        receipt.GlobalPanelAction != expected_global_action ||
                        receipt.SelectorPanelAction !=
                            expected_selector_action ||
                        receipt.OriginalSuperblock != plan_index ||
                        receipt.SelectorSuperblock !=
                            expected_selector_superblock ||
                        receipt.Cell != plan.Cell ||
                        receipt.Round != plan.Round ||
                        receipt.Visit != plan.Visit ||
                        receipt.U != plan.U ||
                        receipt.DirectionValue != plan.DirectionValue ||
                        receipt.ModeValue != plan.ModeValue ||
                        receipt.ModePosition != plan.ModePosition ||
                        receipt.OrientationFlip != plan.OrientationFlip ||
                        receipt.Rotation != plan.Rotation ||
                        receipt.FirstSelector != plan.FirstSelector ||
                        receipt.SelectorPosition != selector_position ||
                        receipt.Selector != selector ||
                        receipt.SequencePosition != sequence_position ||
                        receipt.Sequence != sequence ||
                        receipt.SequenceFamily != (sequence < 2u ? 0u : 1u) ||
                        receipt.Period != period || receipt.Label != label ||
                        receipt.LogicalAPhysical != logical_a ||
                        receipt.Physical != physical ||
                        receipt.Predecessor != predecessor ||
                        !ContextPairValid(
                            target_cpu,
                            receipt.Result.Before,
                            receipt.Result.After,
                            diagnostic) ||
                        !ValidSemanticBatch128Result(
                            info.ScopeValue,
                            kExpectedIidOverheads[info.BaseIndex],
                            receipt.Result.Batch,
                            true) ||
                        !SummarizeResult(receipt.Result.Batch, summary) ||
                        summary.SuccessCount !=
                            receipt.Result.Summary.SuccessCount ||
                        summary.OverheadSum !=
                            receipt.Result.Summary.OverheadSum ||
                        summary.OverheadMinimum !=
                            receipt.Result.Summary.OverheadMinimum ||
                        summary.OverheadMaximum !=
                            receipt.Result.Summary.OverheadMaximum ||
                        summary.DirectSystematicSum !=
                            receipt.Result.Summary.DirectSystematicSum ||
                        summary.InnerSha256 !=
                            receipt.Result.Summary.InnerSha256 ||
                        !AccumulateMeasured(
                            rebuilt,
                            label,
                            physical,
                            receipt.Result.Batch.CounterDelta))
                    {
                        diagnostic = "post-panel raw receipt semantic mismatch";
                        return false;
                    }
                    ++raw_next;
                }
            }
            if (!FinalizeObservation(rebuilt) ||
                observation_next >= observations.size() ||
                !EqualObservation(rebuilt, observations[observation_next]))
            {
                diagnostic = "post-panel observation mismatch";
                return false;
            }
            ++observation_next;
        }
    }
    if (raw_next != raw.size() || observation_next != observations.size() ||
        global_action != kPanelActionCount)
    {
        diagnostic = "post-panel coverage invalid";
        return false;
    }
    return true;
}

struct LogStats
{
    size_t Count = 0u;
    double Mean = 0.0;
    double Lower = 0.0;
    double Upper = 0.0;
    double GeometricMean = 0.0;
    double GeometricLower = 0.0;
    double GeometricUpper = 0.0;
};

double CriticalValue(size_t count)
{
    return count == 96u ? kT95N96 :
        (count == 576u ? kT95N576 : 0.0);
}

bool ComputeLogStats(
    const double* values,
    size_t count,
    LogStats& stats_out)
{
    const double critical = CriticalValue(count);
    if (!values || count < 2u || critical <= 0.0) return false;
    long double sum = 0.0L;
    for (size_t i = 0u; i < count; ++i)
    {
        if (!std::isfinite(values[i])) return false;
        sum += static_cast<long double>(values[i]);
    }
    const long double mean = sum / static_cast<long double>(count);
    long double squares = 0.0L;
    for (size_t i = 0u; i < count; ++i)
    {
        const long double difference =
            static_cast<long double>(values[i]) - mean;
        squares += difference * difference;
    }
    const long double variance = squares /
        static_cast<long double>(count - 1u);
    const long double standard_error = std::sqrt(
        variance / static_cast<long double>(count));
    const long double half_width =
        static_cast<long double>(critical) * standard_error;
    LogStats stats;
    stats.Count = count;
    stats.Mean = static_cast<double>(mean);
    const long double exact_lower = mean - half_width;
    const long double exact_upper = mean + half_width;
    stats.Lower = static_cast<double>(exact_lower);
    stats.Upper = static_cast<double>(exact_upper);
    if (static_cast<long double>(stats.Lower) > exact_lower) {
        stats.Lower = std::nextafter(
            stats.Lower, -std::numeric_limits<double>::infinity());
    }
    if (static_cast<long double>(stats.Upper) < exact_upper) {
        stats.Upper = std::nextafter(
            stats.Upper, std::numeric_limits<double>::infinity());
    }
    stats.GeometricMean = static_cast<double>(std::exp(mean));
    const long double exact_geometric_lower = std::exp(exact_lower);
    const long double exact_geometric_upper = std::exp(exact_upper);
    stats.GeometricLower = static_cast<double>(exact_geometric_lower);
    stats.GeometricUpper = static_cast<double>(exact_geometric_upper);
    if (static_cast<long double>(stats.GeometricLower) >
        exact_geometric_lower)
    {
        stats.GeometricLower = std::nextafter(stats.GeometricLower, 0.0);
    }
    if (static_cast<long double>(stats.GeometricUpper) <
        exact_geometric_upper)
    {
        stats.GeometricUpper = std::nextafter(
            stats.GeometricUpper,
            std::numeric_limits<double>::infinity());
    }
    if (!std::isfinite(stats.Mean) || !std::isfinite(stats.Lower) ||
        !std::isfinite(stats.Upper) ||
        !std::isfinite(stats.GeometricMean) ||
        !std::isfinite(stats.GeometricLower) ||
        !std::isfinite(stats.GeometricUpper) ||
        stats.Lower > stats.Mean || stats.Mean > stats.Upper ||
        stats.GeometricLower <= 0.0 ||
        stats.GeometricLower > stats.GeometricMean ||
        stats.GeometricMean > stats.GeometricUpper)
    {
        return false;
    }
    stats_out = stats;
    return true;
}

bool OfficialGate(const LogStats& stats)
{
    return stats.Lower >= kLogLowerBound && stats.Lower <= 0.0 &&
        stats.Upper >= 0.0 && stats.Upper <= kLogUpperBound;
}

struct OfficialResult
{
    uint32_t Selector = 0u;
    Estimand EstimandValue = Estimand::Logical;
    Mode ModeValue = Mode::Exposed;
    uint32_t Group = 0u;
    bool Pooled = false;
    uint32_t WidthIndex = 0u;
    Scope ScopeValue = Scope::Encoder;
    LogStats Stats;
    bool Pass = false;
    std::string Identity;
};

bool BuildObservationIndex(
    const std::vector<SuperblockObservation>& observations,
    std::vector<size_t>& by_key)
{
    if (observations.size() != kSelectorSuperblockCount) return false;
    by_key.assign(kSelectorSuperblockCount,
        std::numeric_limits<size_t>::max());
    for (size_t i = 0u; i < observations.size(); ++i)
    {
        const SuperblockObservation& observation = observations[i];
        const size_t key = ObservationKey(
            observation.Selector,
            observation.Plan.Cell,
            observation.Plan.Round,
            observation.Plan.ModeValue);
        if (key >= by_key.size() ||
            by_key[key] != std::numeric_limits<size_t>::max())
        {
            return false;
        }
        by_key[key] = i;
    }
    return std::find(
        by_key.begin(),
        by_key.end(),
        std::numeric_limits<size_t>::max()) == by_key.end();
}

const SuperblockObservation* FindObservation(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    uint32_t selector,
    uint32_t cell,
    uint32_t round,
    Mode mode)
{
    const size_t key = ObservationKey(selector, cell, round, mode);
    if (key >= by_key.size() || by_key[key] >= observations.size()) {
        return nullptr;
    }
    return &observations[by_key[key]];
}

double ObservationLog(
    const SuperblockObservation& observation,
    Estimand estimand)
{
    return estimand == Estimand::Logical ? observation.LogLogicalAa :
        observation.LogPhysicalXy;
}

bool BuildOfficialResults(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    std::array<OfficialResult, kOfficialGateCount>& results,
    bool& all_pass)
{
    size_t next = 0u;
    size_t n96 = 0u;
    size_t n576 = 0u;
    all_pass = true;
    std::array<double, 576> logs = {{ 0.0 }};
    for (uint32_t selector = 0u; selector < kSelectorCount; ++selector)
    {
        for (uint32_t estimand_index = 0u; estimand_index < 2u;
             ++estimand_index)
        {
            const Estimand estimand = static_cast<Estimand>(estimand_index);
            for (uint32_t mode_index = 0u; mode_index < kModeCount;
                 ++mode_index)
            {
                const Mode mode = static_cast<Mode>(mode_index);
                for (uint32_t width = 0u; width < 2u; ++width)
                {
                    for (uint32_t scope_index = 0u; scope_index < 3u;
                         ++scope_index)
                    {
                        size_t count = 0u;
                        for (uint32_t root = 0u; root < 3u; ++root)
                        {
                            const uint32_t cell =
                                width * 9u + scope_index * 3u + root;
                            for (uint32_t round = 0u;
                                 round < kRoundCount;
                                 ++round)
                            {
                                const SuperblockObservation* observation =
                                    FindObservation(
                                        observations,
                                        by_key,
                                        selector,
                                        cell,
                                        round,
                                        mode);
                                if (!observation || count >= 96u) return false;
                                logs[count++] = ObservationLog(
                                    *observation, estimand);
                            }
                        }
                        OfficialResult result;
                        result.Selector = selector;
                        result.EstimandValue = estimand;
                        result.ModeValue = mode;
                        result.Group = width * 3u + scope_index;
                        result.WidthIndex = width;
                        result.ScopeValue = static_cast<Scope>(scope_index);
                        std::ostringstream identity;
                        identity << "selector=" << SelectorName(selector)
                                 << ";estimand=" << EstimandName(estimand)
                                 << ";mode=" << ModeName(mode)
                                 << ";width=" << kBlockBytes[width]
                                 << ";scope=" << ScopeName(result.ScopeValue);
                        result.Identity = identity.str();
                        if (count != 96u ||
                            !ComputeLogStats(logs.data(), count, result.Stats))
                        {
                            return false;
                        }
                        result.Pass = OfficialGate(result.Stats);
                        all_pass = all_pass && result.Pass;
                        if (next >= results.size()) return false;
                        results[next++] = result;
                        ++n96;
                    }
                }
                size_t count = 0u;
                for (uint32_t cell = 0u; cell < kCellCount; ++cell)
                {
                    for (uint32_t round = 0u; round < kRoundCount; ++round)
                    {
                        const SuperblockObservation* observation =
                            FindObservation(
                                observations,
                                by_key,
                                selector,
                                cell,
                                round,
                                mode);
                        if (!observation || count >= logs.size()) return false;
                        logs[count++] = ObservationLog(*observation, estimand);
                    }
                }
                OfficialResult pooled;
                pooled.Selector = selector;
                pooled.EstimandValue = estimand;
                pooled.ModeValue = mode;
                pooled.Group = 6u;
                pooled.Pooled = true;
                std::ostringstream identity;
                identity << "selector=" << SelectorName(selector)
                         << ";estimand=" << EstimandName(estimand)
                         << ";mode=" << ModeName(mode)
                         << ";group=pooled";
                pooled.Identity = identity.str();
                if (count != 576u ||
                    !ComputeLogStats(logs.data(), count, pooled.Stats))
                {
                    return false;
                }
                pooled.Pass = OfficialGate(pooled.Stats);
                all_pass = all_pass && pooled.Pass;
                if (next >= results.size()) return false;
                results[next++] = pooled;
                ++n576;
            }
        }
    }
    return next == results.size() && n96 == kExpectedN96GateCount &&
        n576 == kExpectedN576GateCount;
}

struct EmptyGate
{
    uint32_t Selector = 0u;
    uint64_t Minimum = 0u;
    uint64_t Median = 0u;
    uint64_t P999 = 0u;
    uint64_t Maximum = 0u;
    uint64_t MinimumMeasured = 0u;
    bool Pass = false;
};

bool BuildEmptyGates(
    const std::vector<EmptyReceipt>& empty,
    const std::vector<RawTiming>& raw,
    std::array<EmptyGate, 2>& gates)
{
    if (empty.size() != static_cast<size_t>(kEmptySampleCount) * 2u ||
        raw.size() != kMeasuredBatchCount) return false;
    for (uint32_t selector = 0u; selector < 2u; ++selector)
    {
        std::vector<uint64_t> values;
        values.reserve(kEmptySampleCount);
        for (size_t i = 0u; i < empty.size(); ++i) {
            if (empty[i].Selector == selector) values.push_back(empty[i].Delta);
        }
        if (values.size() != kEmptySampleCount) return false;
        std::sort(values.begin(), values.end());
        EmptyGate gate;
        gate.Selector = selector;
        gate.Minimum = values.front();
        gate.Median = values[values.size() / 2u];
        gate.P999 = values[kEmptyP999OneBasedIndex - 1u];
        gate.Maximum = values.back();
        gate.MinimumMeasured = std::numeric_limits<uint64_t>::max();
        for (size_t i = 0u; i < raw.size(); ++i) {
            if (raw[i].Selector == selector) {
                gate.MinimumMeasured = std::min(
                    gate.MinimumMeasured,
                    raw[i].Result.Batch.CounterDelta);
            }
        }
        if (gate.P999 > std::numeric_limits<uint64_t>::max() /
                kEmptyOverheadMultiplier ||
            gate.MinimumMeasured == std::numeric_limits<uint64_t>::max())
        {
            return false;
        }
        gate.Pass = kEmptyOverheadMultiplier * gate.P999 <=
            gate.MinimumMeasured;
        gates[selector] = gate;
    }
    return true;
}

void AppendContext(std::string& bytes, const ContextReceipt& receipt)
{
    AppendLe32(bytes, static_cast<uint32_t>(receipt.Cpu));
    AppendLe32(bytes, receipt.AffinityCount);
    AppendLe32(bytes, receipt.TargetInAffinity);
    AppendLe64(bytes, static_cast<uint64_t>(receipt.Voluntary));
    AppendLe64(bytes, static_cast<uint64_t>(receipt.Involuntary));
}

void AppendBatch(
    std::string& bytes,
    const Batch128CounterArmResult& result)
{
    AppendLe64(bytes, result.CounterStart);
    AppendLe64(bytes, result.CounterEnd);
    AppendLe64(bytes, result.CounterDelta);
    AppendLe32(bytes, result.BeginReadSucceeded ? 1u : 0u);
    AppendLe32(bytes, result.EndReadSucceeded ? 1u : 0u);
    for (size_t i = 0u; i < kBatchSize; ++i) {
        AppendLe32(bytes, static_cast<uint32_t>(result.Results[i]));
    }
    AppendLe64(bytes, result.VerificationMask[0]);
    AppendLe64(bytes, result.VerificationMask[1]);
    for (size_t i = 0u; i < kBatchSize; ++i) {
        AppendLe32(bytes, result.DecodedOverheads[i]);
    }
    for (size_t i = 0u; i < kBatchSize; ++i) {
        AppendLe64(bytes, result.DirectSystematicPackets[i]);
    }
}

std::string RawSha256WithBracketProtocol(
    const std::vector<RawTiming>& raw,
    const std::string& bracket_protocol)
{
    std::string bytes;
    static const char marker[] = "wirehair-rdpru-raw-le-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    AppendLe64(bytes, static_cast<uint64_t>(raw.size()));
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        const RawTiming& receipt = raw[i];
        if (receipt.Cell >= kCellCount) return std::string();
        const LogicalCellInfo info = GetLogicalCellInfo(receipt.Cell);
        AppendLe64(bytes, static_cast<uint64_t>(i));
        AppendLe64(bytes, receipt.RawIndex);
        AppendLe64(bytes, receipt.GlobalPanelAction);
        AppendLe64(bytes, receipt.SelectorPanelAction);
        AppendLe64(bytes, receipt.OriginalSuperblock);
        AppendLe64(bytes, receipt.SelectorSuperblock);
        AppendLe32(bytes, receipt.Cell);
        AppendLe32(bytes, info.WidthIndex);
        AppendLe32(bytes, info.BlockBytes);
        AppendLe32(bytes, static_cast<uint32_t>(info.ScopeValue));
        AppendLe32(bytes, info.RootIndex);
        AppendLe64(bytes, info.Root);
        AppendLe32(bytes, receipt.Round);
        AppendLe32(bytes, receipt.Visit);
        AppendLe32(bytes, receipt.U);
        AppendLe32(bytes, static_cast<uint32_t>(receipt.DirectionValue));
        AppendLe32(bytes, static_cast<uint32_t>(receipt.ModeValue));
        AppendLe32(bytes, receipt.ModePosition);
        AppendLe32(bytes, receipt.OrientationFlip);
        AppendLe32(bytes, receipt.Rotation);
        AppendLe32(bytes, receipt.FirstSelector);
        AppendLe32(bytes, receipt.SelectorPosition);
        AppendLe32(bytes, receipt.Selector);
        AppendLe32(bytes, receipt.SequencePosition);
        AppendLe32(bytes, receipt.Sequence);
        AppendLe32(bytes, receipt.SequenceFamily);
        AppendLe32(bytes, receipt.Period);
        AppendLe32(bytes, static_cast<uint32_t>(receipt.Label));
        AppendLe32(bytes, static_cast<uint32_t>(receipt.LogicalAPhysical));
        AppendLe32(bytes, static_cast<uint32_t>(receipt.Physical));
        AppendLe32(bytes, static_cast<uint32_t>(receipt.Predecessor));
        AppendLe32(bytes, kBatchSize);
        AppendLe32(bytes, kBracketProtocolId);
        AppendSizedString(bytes, bracket_protocol);
        AppendContext(bytes, receipt.Result.Before);
        AppendContext(bytes, receipt.Result.After);
        AppendBatch(bytes, receipt.Result.Batch);
        AppendLe32(bytes, receipt.Result.Summary.SuccessCount);
        AppendLe64(bytes, receipt.Result.Summary.OverheadSum);
        AppendLe32(bytes, receipt.Result.Summary.OverheadMinimum);
        AppendLe32(bytes, receipt.Result.Summary.OverheadMaximum);
        AppendLe64(bytes, receipt.Result.Summary.DirectSystematicSum);
        AppendSizedString(bytes, receipt.Result.Summary.InnerSha256);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

std::string RawSha256(const std::vector<RawTiming>& raw)
{
    return RawSha256WithBracketProtocol(raw, kBracketProtocol);
}

std::string EmptySha256WithBracketProtocol(
    const std::vector<EmptyReceipt>& empty,
    const std::string& bracket_protocol)
{
    std::string bytes;
    static const char marker[] = "wirehair-rdpru-empty-le-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    AppendLe64(bytes, static_cast<uint64_t>(empty.size()));
    for (size_t i = 0u; i < empty.size(); ++i)
    {
        const EmptyReceipt& receipt = empty[i];
        AppendLe64(bytes, static_cast<uint64_t>(i));
        AppendLe32(bytes, receipt.Selector);
        AppendLe32(bytes, receipt.Sample);
        AppendLe32(bytes, receipt.SelectorPosition);
        AppendLe32(bytes, kBracketProtocolId);
        AppendSizedString(bytes, bracket_protocol);
        AppendLe64(bytes, receipt.Start);
        AppendLe64(bytes, receipt.End);
        AppendLe64(bytes, receipt.Delta);
        AppendLe32(bytes, receipt.BeginReadSucceeded);
        AppendLe32(bytes, receipt.EndReadSucceeded);
        AppendContext(bytes, receipt.Before);
        AppendContext(bytes, receipt.After);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

std::string EmptySha256(const std::vector<EmptyReceipt>& empty)
{
    return EmptySha256WithBracketProtocol(empty, kBracketProtocol);
}

bool ValidateEmptyReceipts(
    const std::vector<EmptyReceipt>& empty,
    int target_cpu,
    std::string& diagnostic)
{
    if (empty.size() != static_cast<size_t>(kEmptySampleCount) * 2u)
    {
        diagnostic = "empty receipt cardinality invalid";
        return false;
    }
    for (uint32_t sample = 0u; sample < kEmptySampleCount; ++sample)
    {
        const uint32_t first_selector = sample & 1u;
        for (uint32_t position = 0u; position < 2u; ++position)
        {
            const size_t index = static_cast<size_t>(sample) * 2u + position;
            const EmptyReceipt& receipt = empty[index];
            const uint32_t selector = first_selector ^ position;
            if (receipt.Sample != sample ||
                receipt.SelectorPosition != position ||
                receipt.Selector != selector || receipt.End <= receipt.Start ||
                receipt.Delta == 0u ||
                receipt.Delta != receipt.End - receipt.Start ||
                receipt.BeginReadSucceeded != 1u ||
                receipt.EndReadSucceeded != 1u ||
                !ContextPairValid(
                    target_cpu, receipt.Before, receipt.After, diagnostic))
            {
                diagnostic = "empty receipt order, counter, or context invalid";
                return false;
            }
        }
    }
    return true;
}

bool ValidCompactSemantics(const CompactHardwareReceipt& receipt)
{
    if (receipt.Selector >= 2u || receipt.Cell >= kCellCount ||
        static_cast<uint32_t>(receipt.ScopeValue) >= 3u ||
        static_cast<uint32_t>(receipt.Bank) >= 3u ||
        receipt.CounterEnd <= receipt.CounterStart ||
        receipt.CounterDelta == 0u ||
        receipt.CounterDelta != receipt.CounterEnd - receipt.CounterStart ||
        receipt.BeginReadSucceeded != 1u || receipt.EndReadSucceeded != 1u ||
        receipt.VerificationMask0 != UINT64_MAX ||
        receipt.VerificationMask1 != UINT64_MAX ||
        receipt.Summary.SuccessCount != kBatchSize ||
        receipt.Summary.DirectSystematicSum != 0u ||
        !IsLowerHexSha256(receipt.Summary.InnerSha256.c_str()))
    {
        return false;
    }
    const LogicalCellInfo info = GetLogicalCellInfo(receipt.Cell);
    if (receipt.ScopeValue != info.ScopeValue) return false;
    uint32_t expected_overhead = 0u;
    if (receipt.ScopeValue == Scope::Encoder) expected_overhead = UINT32_MAX;
    else if (receipt.ScopeValue == Scope::IidReceive) {
        expected_overhead = kExpectedIidOverheads[info.BaseIndex];
    }
    const uint64_t expected_sum =
        static_cast<uint64_t>(expected_overhead) * kBatchSize;
    return receipt.Summary.OverheadMinimum == expected_overhead &&
        receipt.Summary.OverheadMaximum == expected_overhead &&
        receipt.Summary.OverheadSum == expected_sum;
}

bool ValidateCompactReceipts(
    const FrozenSchedule& schedule,
    const std::vector<CompactHardwareReceipt>& warm,
    const std::vector<CompactHardwareReceipt>& washout,
    const std::vector<RawTiming>& raw,
    int target_cpu,
    std::string& diagnostic)
{
    if (warm.size() != kWarmBatchCount ||
        washout.size() != kWashoutBatchCount ||
        raw.size() != kMeasuredBatchCount)
    {
        diagnostic = "warm/washout/measured receipt cardinality invalid";
        return false;
    }
    std::array<size_t, 2> warm_by_selector = {{ 0u, 0u }};
    for (size_t i = 0u; i < warm.size(); ++i)
    {
        const CompactHardwareReceipt& receipt = warm[i];
        const uint32_t expected_cell = static_cast<uint32_t>(i / 6u);
        const uint32_t within_cell = static_cast<uint32_t>(i % 6u);
        const uint32_t expected_bank = within_cell / 2u;
        const uint32_t expected_position = within_cell % 2u;
        const uint32_t expected_selector =
            ((expected_cell ^ expected_bank) & 1u) ^ expected_position;
        const LogicalCellInfo expected_info =
            GetLogicalCellInfo(expected_cell);
        if (receipt.Kind != CompactKind::Warm || receipt.Ordinal != i ||
            receipt.GlobalPanelAction != UINT64_MAX ||
            receipt.SelectorPanelAction != UINT64_MAX ||
            receipt.OriginalSuperblock != UINT64_MAX ||
            receipt.SelectorSuperblock != UINT64_MAX ||
            receipt.Cell != expected_cell ||
            receipt.ScopeValue != expected_info.ScopeValue ||
            receipt.Bank != static_cast<PhysicalBank>(expected_bank) ||
            receipt.Selector != expected_selector ||
            !ValidCompactSemantics(receipt) ||
            !ContextPairValid(
                target_cpu, receipt.Before, receipt.After, diagnostic))
        {
            diagnostic = "warm hardware receipt invalid";
            return false;
        }
        ++warm_by_selector[receipt.Selector];
    }
    if (warm_by_selector[0] != kWarmBatchCount / 2u ||
        warm_by_selector[1] != kWarmBatchCount / 2u)
    {
        diagnostic = "warm selector balance invalid";
        return false;
    }

    std::vector<uint8_t> global_seen(kPanelActionCount, 0u);
    std::array<std::vector<uint8_t>, 2> selector_seen;
    selector_seen[0].assign(
        kMeasuredBatchCountPerSelector + kWashoutBatchCountPerSelector, 0u);
    selector_seen[1].assign(
        kMeasuredBatchCountPerSelector + kWashoutBatchCountPerSelector, 0u);
    std::array<size_t, 2> washout_by_selector = {{ 0u, 0u }};
    for (size_t i = 0u; i < washout.size(); ++i)
    {
        const CompactHardwareReceipt& receipt = washout[i];
        if (receipt.Selector >= 2u ||
            receipt.Kind != CompactKind::Washout || receipt.Ordinal != i ||
            receipt.Bank != PhysicalBank::W ||
            receipt.GlobalPanelAction >= global_seen.size() ||
            receipt.SelectorPanelAction >=
                selector_seen[receipt.Selector].size() ||
            receipt.OriginalSuperblock >= kOriginalSuperblockCount ||
            receipt.SelectorSuperblock >= kOriginalSuperblockCount ||
            !ValidCompactSemantics(receipt) ||
            !ContextPairValid(
                target_cpu, receipt.Before, receipt.After, diagnostic) ||
            global_seen[receipt.GlobalPanelAction] != 0u ||
            selector_seen[receipt.Selector][receipt.SelectorPanelAction] != 0u)
        {
            diagnostic = "washout hardware receipt invalid";
            return false;
        }
        global_seen[receipt.GlobalPanelAction] = 1u;
        selector_seen[receipt.Selector][receipt.SelectorPanelAction] = 1u;
        ++washout_by_selector[receipt.Selector];
    }

    size_t expected_washout_next = 0u;
    uint64_t expected_global_action = 0u;
    std::array<uint64_t, 2> expected_selector_action = {{ 0u, 0u }};
    std::array<uint64_t, 2> expected_selector_superblock = {{ 0u, 0u }};
    for (size_t plan_index = 0u; plan_index < schedule.size(); ++plan_index)
    {
        const SuperblockPlan& plan = schedule[plan_index];
        const LogicalCellInfo info = GetLogicalCellInfo(plan.Cell);
        for (uint32_t position = 0u; position < 2u; ++position)
        {
            const uint32_t selector = plan.FirstSelector ^ position;
            const uint64_t selector_superblock =
                expected_selector_superblock[selector]++;
            for (uint32_t sequence_position = 0u;
                 sequence_position < kSequenceCount;
                 ++sequence_position)
            {
                (void)sequence_position;
                for (uint32_t step = 0u; step < kBatchesPerSequence; ++step)
                {
                    const bool measured = plan.ModeValue == Mode::Exposed ?
                        step >= kPeriodsPerSequence : (step & 1u) != 0u;
                    const uint64_t global_action = expected_global_action++;
                    const uint64_t selector_action =
                        expected_selector_action[selector]++;
                    if (measured) continue;
                    if (expected_washout_next >= washout.size())
                    {
                        diagnostic = "washout schedule receipt underflow";
                        return false;
                    }
                    const CompactHardwareReceipt& receipt =
                        washout[expected_washout_next];
                    if (receipt.Ordinal != expected_washout_next ||
                        receipt.GlobalPanelAction != global_action ||
                        receipt.SelectorPanelAction != selector_action ||
                        receipt.OriginalSuperblock != plan_index ||
                        receipt.SelectorSuperblock != selector_superblock ||
                        receipt.Selector != selector ||
                        receipt.Cell != plan.Cell ||
                        receipt.ScopeValue != info.ScopeValue ||
                        receipt.Bank != PhysicalBank::W)
                    {
                        diagnostic = "washout schedule identity mismatch";
                        return false;
                    }
                    ++expected_washout_next;
                }
            }
        }
    }
    if (expected_washout_next != washout.size() ||
        expected_global_action != kPanelActionCount ||
        expected_selector_action[0] !=
            kMeasuredBatchCountPerSelector + kWashoutBatchCountPerSelector ||
        expected_selector_action[1] !=
            kMeasuredBatchCountPerSelector + kWashoutBatchCountPerSelector)
    {
        diagnostic = "washout schedule coverage invalid";
        return false;
    }
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        const RawTiming& receipt = raw[i];
        if (receipt.Selector >= 2u ||
            receipt.GlobalPanelAction >= global_seen.size() ||
            receipt.SelectorPanelAction >=
                selector_seen[receipt.Selector].size() ||
            global_seen[receipt.GlobalPanelAction] != 0u ||
            selector_seen[receipt.Selector][receipt.SelectorPanelAction] != 0u)
        {
            diagnostic = "measured/washout action coverage overlaps";
            return false;
        }
        global_seen[receipt.GlobalPanelAction] = 1u;
        selector_seen[receipt.Selector][receipt.SelectorPanelAction] = 1u;
    }
    if (washout_by_selector[0] != kWashoutBatchCountPerSelector ||
        washout_by_selector[1] != kWashoutBatchCountPerSelector ||
        std::find(global_seen.begin(), global_seen.end(), 0u) !=
            global_seen.end() ||
        std::find(selector_seen[0].begin(), selector_seen[0].end(), 0u) !=
            selector_seen[0].end() ||
        std::find(selector_seen[1].begin(), selector_seen[1].end(), 0u) !=
            selector_seen[1].end())
    {
        diagnostic = "hardware action coverage incomplete";
        return false;
    }
    return true;
}

void AppendCompact(
    std::string& bytes,
    const CompactHardwareReceipt& receipt,
    const std::string& bracket_protocol)
{
    AppendLe32(bytes, static_cast<uint32_t>(receipt.Kind));
    AppendLe64(bytes, receipt.Ordinal);
    AppendLe64(bytes, receipt.GlobalPanelAction);
    AppendLe64(bytes, receipt.SelectorPanelAction);
    AppendLe64(bytes, receipt.OriginalSuperblock);
    AppendLe64(bytes, receipt.SelectorSuperblock);
    AppendLe32(bytes, receipt.Selector);
    AppendLe32(bytes, receipt.Cell);
    AppendLe32(bytes, static_cast<uint32_t>(receipt.ScopeValue));
    AppendLe32(bytes, static_cast<uint32_t>(receipt.Bank));
    AppendLe32(bytes, kBracketProtocolId);
    AppendSizedString(bytes, bracket_protocol);
    AppendContext(bytes, receipt.Before);
    AppendContext(bytes, receipt.After);
    AppendLe64(bytes, receipt.CounterStart);
    AppendLe64(bytes, receipt.CounterEnd);
    AppendLe64(bytes, receipt.CounterDelta);
    AppendLe32(bytes, receipt.BeginReadSucceeded);
    AppendLe32(bytes, receipt.EndReadSucceeded);
    AppendLe64(bytes, receipt.VerificationMask0);
    AppendLe64(bytes, receipt.VerificationMask1);
    AppendLe32(bytes, receipt.Summary.SuccessCount);
    AppendLe64(bytes, receipt.Summary.OverheadSum);
    AppendLe32(bytes, receipt.Summary.OverheadMinimum);
    AppendLe32(bytes, receipt.Summary.OverheadMaximum);
    AppendLe64(bytes, receipt.Summary.DirectSystematicSum);
    AppendSizedString(bytes, receipt.Summary.InnerSha256);
}

std::string CompactSha256WithBracketProtocol(
    const std::vector<CompactHardwareReceipt>& warm,
    const std::vector<CompactHardwareReceipt>& washout,
    const std::string& bracket_protocol)
{
    std::string bytes;
    static const char marker[] = "wirehair-rdpru-compact-hardware-le-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    AppendLe64(bytes, static_cast<uint64_t>(warm.size()));
    for (size_t i = 0u; i < warm.size(); ++i) {
        AppendCompact(bytes, warm[i], bracket_protocol);
    }
    AppendLe64(bytes, static_cast<uint64_t>(washout.size()));
    for (size_t i = 0u; i < washout.size(); ++i) {
        AppendCompact(bytes, washout[i], bracket_protocol);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

std::string CompactSha256(
    const std::vector<CompactHardwareReceipt>& warm,
    const std::vector<CompactHardwareReceipt>& washout)
{
    return CompactSha256WithBracketProtocol(
        warm, washout, kBracketProtocol);
}

bool MeanValues(
    const std::vector<double>& values,
    double& mean_out)
{
    if (values.empty()) return false;
    long double sum = 0.0L;
    for (size_t i = 0u; i < values.size(); ++i)
    {
        if (!std::isfinite(values[i])) return false;
        sum += static_cast<long double>(values[i]);
    }
    const long double mean = sum / static_cast<long double>(values.size());
    if (!std::isfinite(mean)) return false;
    mean_out = static_cast<double>(mean);
    return true;
}

void AppendDiagnosticMean(
    std::ostringstream& output,
    const char* identity,
    const std::vector<double>& values)
{
    double mean = 0.0;
    output << "DIAGNOSTIC,identity=" << identity
           << ",n=" << values.size();
    if (MeanValues(values, mean)) {
        output << ",mean_log=" << std::setprecision(17) << mean
               << ",geometric_mean=" << std::exp(mean);
    }
    else {
        output << ",mean_log=INVALID,geometric_mean=INVALID";
    }
    output << ",rescuing=false\n";
}

bool AppendDiagnostics(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    const std::vector<RawTiming>& raw,
    std::ostringstream& output)
{
    std::vector<double> values;
    for (uint32_t selector = 0u; selector < 2u; ++selector)
    {
        for (uint32_t estimand_index = 0u; estimand_index < 2u;
             ++estimand_index)
        {
            const Estimand estimand = static_cast<Estimand>(estimand_index);
            for (uint32_t mode_index = 0u; mode_index < 2u; ++mode_index)
            {
                const Mode mode = static_cast<Mode>(mode_index);
                for (uint32_t width = 0u; width < 2u; ++width)
                {
                    for (uint32_t scope = 0u; scope < 3u; ++scope)
                    {
                        for (uint32_t root = 0u; root < 3u; ++root)
                        {
                            values.clear();
                            const uint32_t cell = width * 9u + scope * 3u + root;
                            for (uint32_t round = 0u; round < 32u; ++round)
                            {
                                const SuperblockObservation* observation =
                                    FindObservation(
                                        observations, by_key, selector, cell,
                                        round, mode);
                                if (!observation) return false;
                                values.push_back(ObservationLog(
                                    *observation, estimand));
                            }
                            std::ostringstream identity;
                            identity << "root;selector=" << SelectorName(selector)
                                     << ";estimand=" << EstimandName(estimand)
                                     << ";mode=" << ModeName(mode)
                                     << ";width=" << kBlockBytes[width]
                                     << ";scope=" << ScopeName(
                                            static_cast<Scope>(scope))
                                     << ";root=" << root;
                            AppendDiagnosticMean(
                                output, identity.str().c_str(), values);
                        }
                    }
                }
                for (uint32_t half = 0u; half < 2u; ++half)
                {
                    values.clear();
                    for (uint32_t cell = 0u; cell < kCellCount; ++cell)
                    {
                        for (uint32_t round = half * 16u;
                             round < (half + 1u) * 16u;
                             ++round)
                        {
                            const SuperblockObservation* observation =
                                FindObservation(
                                    observations, by_key, selector, cell,
                                    round, mode);
                            if (!observation) return false;
                            values.push_back(ObservationLog(
                                *observation, estimand));
                        }
                    }
                    std::ostringstream identity;
                    identity << "round-half;selector=" << SelectorName(selector)
                             << ";estimand=" << EstimandName(estimand)
                             << ";mode=" << ModeName(mode)
                             << ";half=" << half;
                    AppendDiagnosticMean(output, identity.str().c_str(), values);
                }
                for (uint32_t position = 0u; position < 2u; ++position)
                {
                    values.clear();
                    for (size_t i = 0u; i < observations.size(); ++i)
                    {
                        const SuperblockObservation& observation = observations[i];
                        if (observation.Selector == selector &&
                            observation.SelectorPosition == position &&
                            observation.Plan.ModeValue == mode)
                        {
                            values.push_back(ObservationLog(
                                observation, estimand));
                        }
                    }
                    std::ostringstream identity;
                    identity << "selector-order;selector="
                             << SelectorName(selector)
                             << ";estimand=" << EstimandName(estimand)
                             << ";mode=" << ModeName(mode)
                             << ";position=" << position;
                    AppendDiagnosticMean(output, identity.str().c_str(), values);
                }
            }
            values.clear();
            for (uint32_t cell = 0u; cell < kCellCount; ++cell)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round)
                {
                    const SuperblockObservation* exposed = FindObservation(
                        observations, by_key, selector, cell, round,
                        Mode::Exposed);
                    const SuperblockObservation* isolated = FindObservation(
                        observations, by_key, selector, cell, round,
                        Mode::Isolated);
                    if (!exposed || !isolated) return false;
                    values.push_back(
                        ObservationLog(*exposed, estimand) -
                        ObservationLog(*isolated, estimand));
                }
            }
            std::ostringstream identity;
            identity << "mode-exposed-minus-isolated;selector="
                     << SelectorName(selector)
                     << ";estimand=" << EstimandName(estimand);
            AppendDiagnosticMean(output, identity.str().c_str(), values);
        }
    }
    for (uint32_t estimand_index = 0u; estimand_index < 2u; ++estimand_index)
    {
        const Estimand estimand = static_cast<Estimand>(estimand_index);
        for (uint32_t mode_index = 0u; mode_index < 2u; ++mode_index)
        {
            const Mode mode = static_cast<Mode>(mode_index);
            values.clear();
            for (uint32_t cell = 0u; cell < kCellCount; ++cell)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round)
                {
                    const SuperblockObservation* mperf = FindObservation(
                        observations, by_key, 0u, cell, round, mode);
                    const SuperblockObservation* aperf = FindObservation(
                        observations, by_key, 1u, cell, round, mode);
                    if (!mperf || !aperf) return false;
                    values.push_back(
                        ObservationLog(*mperf, estimand) -
                        ObservationLog(*aperf, estimand));
                }
            }
            std::ostringstream identity;
            identity << "counter-mperf-minus-aperf;estimand="
                     << EstimandName(estimand)
                     << ";mode=" << ModeName(mode);
            AppendDiagnosticMean(output, identity.str().c_str(), values);
        }
    }
    for (uint32_t selector = 0u; selector < 2u; ++selector)
    {
        uint64_t minimum = std::numeric_limits<uint64_t>::max();
        uint64_t maximum = 0u;
        size_t count = 0u;
        for (size_t i = 0u; i < raw.size(); ++i)
        {
            if (raw[i].Selector != selector) continue;
            minimum = std::min(minimum, raw[i].Result.Batch.CounterDelta);
            maximum = std::max(maximum, raw[i].Result.Batch.CounterDelta);
            ++count;
        }
        if (count != kMeasuredBatchCountPerSelector) return false;
        output << "COUNTER_DIAGNOSTIC,selector=" << SelectorName(selector)
               << ",n=" << count << ",minimum=" << minimum
               << ",maximum=" << maximum << ",rescuing=false\n";
    }
    return true;
}

bool IsCanonicalAbsolutePath(const std::string& path)
{
    if (path.size() < 2u || path[0] != '/' || path.back() == '/' ||
        path.find('\0') != std::string::npos) return false;
    size_t begin = 1u;
    while (begin < path.size())
    {
        const size_t slash = path.find('/', begin);
        const size_t end = slash == std::string::npos ? path.size() : slash;
        if (end == begin) return false;
        const std::string component = path.substr(begin, end - begin);
        if (component == "." || component == "..") return false;
        if (slash == std::string::npos) break;
        begin = slash + 1u;
    }
    return true;
}

#if defined(__linux__)
bool SplitOutputPath(
    const std::string& output_path,
    std::string& directory,
    std::string& basename)
{
    if (!IsCanonicalAbsolutePath(output_path)) {
        return false;
    }
    const size_t slash = output_path.find_last_of('/');
    if (slash == std::string::npos)
    {
        directory = ".";
        basename = output_path;
    }
    else
    {
        directory = slash == 0u ? "/" : output_path.substr(0u, slash);
        basename = output_path.substr(slash + 1u);
    }
    return !basename.empty() && basename != "." && basename != "..";
}
#endif

enum class AtomicFaultPoint
{
    None,
    Write,
    FileFsync,
    FileFsyncEintrOnce,
    Rename,
    DirectoryFsync,
    DirectoryFsyncEintrOnce,
    DirectoryFsyncCleanupUnlink,
    DirectoryFsyncCleanupFsync,
    DirectoryFsyncCleanupFsyncEintrOnce
};

#if defined(__linux__)
bool FsyncRetry(int file, bool inject_eintr_once)
{
#if defined(__linux__)
    bool injected = false;
    for (;;)
    {
        int result = 0;
        if (inject_eintr_once && !injected)
        {
            injected = true;
            errno = EINTR;
            result = -1;
        }
        else {
            result = fsync(file);
        }
        if (result == 0) return true;
        if (errno != EINTR) return false;
    }
#else
    (void)file;
    (void)inject_eintr_once;
    return false;
#endif
}

bool WriteAll(
    int file,
    const char* data,
    size_t size,
    AtomicFaultPoint fault,
    std::string& diagnostic)
{
#if defined(__linux__)
    if (fault == AtomicFaultPoint::Write)
    {
        diagnostic = "injected write evidence failure";
        return false;
    }
    size_t offset = 0u;
    while (offset < size)
    {
        const ssize_t written = write(file, data + offset, size - offset);
        if (written < 0)
        {
            if (errno == EINTR) continue;
            diagnostic = std::string("write evidence failed: ") +
                std::strerror(errno);
            return false;
        }
        if (written == 0)
        {
            diagnostic = "write evidence made no progress";
            return false;
        }
        offset += static_cast<size_t>(written);
    }
    return true;
#else
    (void)file;
    (void)data;
    (void)size;
    (void)fault;
    diagnostic = "atomic evidence requires Linux";
    return false;
#endif
}
#endif

bool AtomicPublishNewImpl(
    const std::string& output_path,
    const std::string& evidence,
    AtomicFaultPoint fault,
    std::string& diagnostic)
{
#if defined(__linux__) && defined(SYS_renameat2)
    std::string directory;
    std::string basename;
    if (!SplitOutputPath(output_path, directory, basename))
    {
        diagnostic = "invalid output path";
        return false;
    }
    const int directory_fd = open(
        directory.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (directory_fd < 0)
    {
        diagnostic = std::string("open output directory failed: ") +
            std::strerror(errno);
        return false;
    }
    std::ostringstream temporary_name_builder;
    temporary_name_builder << "." << basename << ".tmp."
                           << static_cast<unsigned long long>(getpid());
    const std::string temporary_name = temporary_name_builder.str();
    const int file = openat(
        directory_fd,
        temporary_name.c_str(),
        O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC,
        S_IRUSR | S_IWUSR);
    if (file < 0)
    {
        diagnostic = std::string("create same-directory temporary failed: ") +
            std::strerror(errno);
        (void)close(directory_fd);
        return false;
    }
    struct stat temporary_status;
    if (fstat(file, &temporary_status) != 0)
    {
        diagnostic = std::string("fstat temporary evidence failed: ") +
            std::strerror(errno);
        (void)close(file);
        (void)unlinkat(directory_fd, temporary_name.c_str(), 0);
        (void)close(directory_fd);
        return false;
    }
    bool ok = WriteAll(
        file, evidence.data(), evidence.size(), fault, diagnostic);
    if (ok && (fault == AtomicFaultPoint::FileFsync ||
            !FsyncRetry(
                file, fault == AtomicFaultPoint::FileFsyncEintrOnce)))
    {
        diagnostic = fault == AtomicFaultPoint::FileFsync ?
            "injected fsync evidence failure" :
            std::string("fsync evidence failed: ") + std::strerror(errno);
        ok = false;
    }
    if (close(file) != 0 && ok)
    {
        diagnostic = std::string("close evidence failed: ") +
            std::strerror(errno);
        ok = false;
    }
    if (ok)
    {
        static const unsigned rename_noreplace = 1u;
        if (fault == AtomicFaultPoint::Rename || syscall(
                SYS_renameat2,
                directory_fd,
                temporary_name.c_str(),
                directory_fd,
                basename.c_str(),
                rename_noreplace) != 0)
        {
            diagnostic = fault == AtomicFaultPoint::Rename ?
                "injected atomic rename failure" :
                std::string("atomic no-replace rename failed: ") +
                    std::strerror(errno);
            ok = false;
        }
    }
    bool renamed = ok;
    const bool inject_directory_fsync =
        fault == AtomicFaultPoint::DirectoryFsync ||
        fault == AtomicFaultPoint::DirectoryFsyncCleanupUnlink ||
        fault == AtomicFaultPoint::DirectoryFsyncCleanupFsync ||
        fault == AtomicFaultPoint::DirectoryFsyncCleanupFsyncEintrOnce;
    if (ok && (inject_directory_fsync ||
            !FsyncRetry(
                directory_fd,
                fault == AtomicFaultPoint::DirectoryFsyncEintrOnce)))
    {
        diagnostic = inject_directory_fsync ?
            "injected fsync output directory failure" :
            std::string("fsync output directory failed: ") +
                std::strerror(errno);
        ok = false;
    }
    if (!ok && renamed)
    {
        struct stat published_status;
        const bool same_identity = fstatat(
                directory_fd,
                basename.c_str(),
                &published_status,
                AT_SYMLINK_NOFOLLOW) == 0 &&
            published_status.st_dev == temporary_status.st_dev &&
            published_status.st_ino == temporary_status.st_ino;
        const bool unlinked = same_identity &&
            fault != AtomicFaultPoint::DirectoryFsyncCleanupUnlink &&
            unlinkat(directory_fd, basename.c_str(), 0) == 0;
        const bool cleanup_synced = unlinked &&
            fault != AtomicFaultPoint::DirectoryFsyncCleanupFsync &&
            FsyncRetry(
                directory_fd,
                fault ==
                    AtomicFaultPoint::DirectoryFsyncCleanupFsyncEintrOnce);
        struct stat remaining;
        errno = 0;
        const bool absent = fstatat(
            directory_fd,
            basename.c_str(),
            &remaining,
            AT_SYMLINK_NOFOLLOW) != 0 && errno == ENOENT;
        if (!same_identity || !unlinked || !cleanup_synced || !absent) {
            diagnostic += "; post-rename cleanup could not be proved";
        }
    }
    if (!ok && !renamed) {
        (void)unlinkat(directory_fd, temporary_name.c_str(), 0);
    }
    // Once rename and directory fsync succeeded, close cannot revoke durable
    // publication and therefore does not downgrade the result.
    (void)close(directory_fd);
    return ok;
#else
    (void)output_path;
    (void)evidence;
    (void)fault;
    diagnostic = "atomic renameat2 evidence publication requires Linux";
    return false;
#endif
}

bool AtomicPublishNew(
    const std::string& output_path,
    const std::string& evidence,
    std::string& diagnostic)
{
    return AtomicPublishNewImpl(
        output_path, evidence, AtomicFaultPoint::None, diagnostic);
}

bool CheckOutputDestinationNew(
    const std::string& output_path,
    std::string& diagnostic)
{
#if defined(__linux__)
    std::string directory;
    std::string basename;
    if (!SplitOutputPath(output_path, directory, basename))
    {
        diagnostic = "invalid output path";
        return false;
    }
    const int directory_fd = open(
        directory.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (directory_fd < 0)
    {
        diagnostic = std::string("open output directory failed: ") +
            std::strerror(errno);
        return false;
    }
    struct stat status;
    errno = 0;
    const int result = fstatat(
        directory_fd, basename.c_str(), &status, AT_SYMLINK_NOFOLLOW);
    const int saved_errno = errno;
    const bool close_ok = close(directory_fd) == 0;
    if (result == 0)
    {
        diagnostic = "output path already exists";
        return false;
    }
    if (saved_errno != ENOENT)
    {
        diagnostic = std::string("inspect output path failed: ") +
            std::strerror(saved_errno);
        return false;
    }
    if (!close_ok)
    {
        diagnostic = std::string("close output directory failed: ") +
            std::strerror(errno);
        return false;
    }
    return true;
#else
    (void)output_path;
    diagnostic = "atomic evidence requires Linux";
    return false;
#endif
}

std::string HexEncode(const std::string& value);

std::string FinalizeCalibrationPayload(
    const std::string& accumulated_evidence,
    bool scientific_valid,
    const std::string& diagnostic)
{
    std::ostringstream payload;
    payload << accumulated_evidence;
    if (!accumulated_evidence.empty() &&
        accumulated_evidence.back() != '\n') payload << '\n';
    payload << "CALIBRATION_DIAGNOSTIC_LENGTH=" << diagnostic.size() << '\n'
            << "CALIBRATION_DIAGNOSTIC_HEX=" << HexEncode(diagnostic) << '\n'
            << "SCIENTIFIC_GATE_VERDICT="
            << (scientific_valid ? "PASS" : "FAIL") << '\n'
            << "EVIDENCE_AUTHORITY="
               "REQUIRES_ATOMIC_PUBLICATION_PASS_RECEIPT\n"
            << "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION="
            << (scientific_valid ? "PENDING_PUBLICATION" : "INVALID")
            << '\n';
    return payload.str();
}

std::string HexEncode(const std::string& value)
{
    static const char digits[] = "0123456789abcdef";
    std::string encoded;
    encoded.reserve(value.size() * 2u);
    for (size_t i = 0u; i < value.size(); ++i)
    {
        const unsigned byte = static_cast<unsigned char>(value[i]);
        encoded.push_back(digits[byte >> 4u]);
        encoded.push_back(digits[byte & 15u]);
    }
    return encoded;
}

#if defined(__linux__)
bool ReadLe16At(const std::string& bytes, size_t offset, uint16_t& value)
{
    if (offset > bytes.size() || bytes.size() - offset < 2u) return false;
    value = static_cast<uint16_t>(
        static_cast<unsigned char>(bytes[offset]) |
        (static_cast<uint16_t>(static_cast<unsigned char>(bytes[offset + 1u]))
            << 8u));
    return true;
}

bool ReadLe32At(const std::string& bytes, size_t offset, uint32_t& value)
{
    if (offset > bytes.size() || bytes.size() - offset < 4u) return false;
    value = 0u;
    for (unsigned i = 0u; i < 4u; ++i) {
        value |= static_cast<uint32_t>(
            static_cast<unsigned char>(bytes[offset + i])) << (8u * i);
    }
    return true;
}

bool ReadLe64At(const std::string& bytes, size_t offset, uint64_t& value)
{
    if (offset > bytes.size() || bytes.size() - offset < 8u) return false;
    value = 0u;
    for (unsigned i = 0u; i < 8u; ++i) {
        value |= static_cast<uint64_t>(
            static_cast<unsigned char>(bytes[offset + i])) << (8u * i);
    }
    return true;
}

bool ExtractElfBuildId(
    const std::string& executable,
    std::string& build_id)
{
    if (executable.size() < 64u ||
        static_cast<unsigned char>(executable[0]) != 0x7fu ||
        executable[1] != 'E' || executable[2] != 'L' || executable[3] != 'F' ||
        static_cast<unsigned char>(executable[4]) != 2u ||
        static_cast<unsigned char>(executable[5]) != 1u)
    {
        return false;
    }
    uint64_t program_offset = 0u;
    uint16_t program_entry_size = 0u;
    uint16_t program_count = 0u;
    if (!ReadLe64At(executable, 32u, program_offset) ||
        !ReadLe16At(executable, 54u, program_entry_size) ||
        !ReadLe16At(executable, 56u, program_count) ||
        program_entry_size < 56u) return false;
    for (uint32_t program = 0u; program < program_count; ++program)
    {
        if (program_offset > std::numeric_limits<uint64_t>::max() -
                static_cast<uint64_t>(program) * program_entry_size)
        {
            return false;
        }
        const uint64_t entry64 = program_offset +
            static_cast<uint64_t>(program) * program_entry_size;
        if (entry64 > executable.size() ||
            executable.size() - static_cast<size_t>(entry64) < 56u)
        {
            return false;
        }
        const size_t entry = static_cast<size_t>(entry64);
        uint32_t type = 0u;
        uint64_t note_offset64 = 0u;
        uint64_t note_size64 = 0u;
        if (!ReadLe32At(executable, entry, type) || type != 4u) continue;
        if (!ReadLe64At(executable, entry + 8u, note_offset64) ||
            !ReadLe64At(executable, entry + 32u, note_size64) ||
            note_offset64 > executable.size() ||
            note_size64 > executable.size() -
                static_cast<size_t>(note_offset64))
        {
            return false;
        }
        size_t next = static_cast<size_t>(note_offset64);
        const size_t end = next + static_cast<size_t>(note_size64);
        while (next < end)
        {
            if (end - next < 12u) return false;
            uint32_t name_size = 0u;
            uint32_t descriptor_size = 0u;
            uint32_t note_type = 0u;
            if (!ReadLe32At(executable, next, name_size) ||
                !ReadLe32At(executable, next + 4u, descriptor_size) ||
                !ReadLe32At(executable, next + 8u, note_type)) return false;
            next += 12u;
            const uint64_t padded_name =
                (static_cast<uint64_t>(name_size) + 3u) & ~UINT64_C(3);
            const uint64_t padded_descriptor =
                (static_cast<uint64_t>(descriptor_size) + 3u) & ~UINT64_C(3);
            if (padded_name > end - next) return false;
            const size_t name_offset = next;
            next += static_cast<size_t>(padded_name);
            if (padded_descriptor > end - next) return false;
            const size_t descriptor_offset = next;
            next += static_cast<size_t>(padded_descriptor);
            if (note_type == 3u && name_size == 4u && descriptor_size > 0u &&
                std::memcmp(executable.data() + name_offset, "GNU\0", 4u) == 0)
            {
                build_id = HexEncode(executable.substr(
                    descriptor_offset, descriptor_size));
                return IsLowerHex(build_id) && !IsAllZeroHex(build_id.c_str());
            }
        }
    }
    return false;
}
#endif

bool ReadSelfExecutableIdentity(
    std::string& sha256,
    std::string& build_id,
    std::string& diagnostic)
{
#if defined(__linux__)
    const int file = open("/proc/self/exe", O_RDONLY | O_CLOEXEC);
    if (file < 0)
    {
        diagnostic = std::string("open /proc/self/exe failed: ") +
            std::strerror(errno);
        return false;
    }
    struct stat status;
    if (fstat(file, &status) != 0 || !S_ISREG(status.st_mode) ||
        status.st_size <= 0 || status.st_size > 256 * 1024 * 1024)
    {
        diagnostic = "self executable is not a bounded regular file";
        (void)close(file);
        return false;
    }
    std::string executable(static_cast<size_t>(status.st_size), '\0');
    size_t offset = 0u;
    while (offset < executable.size())
    {
        const ssize_t got = read(
            file, &executable[offset], executable.size() - offset);
        if (got < 0)
        {
            if (errno == EINTR) continue;
            diagnostic = std::string("read self executable failed: ") +
                std::strerror(errno);
            (void)close(file);
            return false;
        }
        if (got == 0)
        {
            diagnostic = "self executable truncated during read";
            (void)close(file);
            return false;
        }
        offset += static_cast<size_t>(got);
    }
    if (close(file) != 0)
    {
        diagnostic = std::string("close self executable failed: ") +
            std::strerror(errno);
        return false;
    }
    sha256 = wirehair::wh2_benchmark::Sha256Hex(executable);
    if (!IsLowerHexSha256(sha256.c_str()) ||
        !ExtractElfBuildId(executable, build_id))
    {
        diagnostic = "self executable SHA-256/build-ID extraction failed";
        return false;
    }
    return true;
#else
    (void)sha256;
    (void)build_id;
    diagnostic = "self executable identity requires Linux ELF";
    return false;
#endif
}

struct LaunchManifest
{
    std::string Path;
    std::string Sha256;
    std::string Content;
    std::map<std::string, std::string> Fields;
};

std::string DecimalU32(uint32_t value);

#if defined(__linux__)
bool ParseLaunchManifestContent(
    const std::string& content,
    int target_cpu,
    const std::string& output_path,
    const std::string& executable_sha256,
    const std::string& executable_build_id,
    std::map<std::string, std::string>& parsed_fields,
    std::string& diagnostic)
{
    if (content.empty() || content.back() != '\n' ||
        content.find('\0') != std::string::npos)
    {
        diagnostic = "launch manifest is empty, unterminated, or contains NUL";
        return false;
    }
    std::map<std::string, std::string> fields;
    size_t offset = 0u;
    while (offset < content.size())
    {
        const size_t newline = content.find('\n', offset);
        if (newline == std::string::npos) return false;
        const std::string line = content.substr(offset, newline - offset);
        offset = newline + 1u;
        if (line.empty())
        {
            diagnostic = "launch manifest contains a blank line";
            return false;
        }
        const size_t equals = line.find('=');
        if (equals == std::string::npos || equals == 0u ||
            equals + 1u >= line.size() ||
            line.find('=', equals + 1u) != std::string::npos)
        {
            diagnostic = "launch manifest line is malformed";
            return false;
        }
        const std::string key = line.substr(0u, equals);
        const std::string value = line.substr(equals + 1u);
        if (!fields.insert(std::make_pair(key, value)).second)
        {
            diagnostic = "launch manifest contains duplicate key";
            return false;
        }
    }
    // Preserve every syntactically parsed claim even if a later exact-value
    // or schema gate fails, so failed identities survive in INVALID evidence.
    parsed_fields = fields;
    std::ostringstream cpu;
    cpu << target_cpu;
    const struct ExpectedField
    {
        const char* Key;
        std::string Value;
    } expected[] = {
        { "manifest_schema", kLaunchManifestSchema },
        { "protocol_tag", kProtocolTag },
        { "git_commit", WIREHAIR_WH2_SOURCE_GIT_COMMIT },
        { "controller_sha256", WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256 },
        { "codec_header_sha256", WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256 },
        { "codec_source_sha256", WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256 },
        { "schedule_sha256", kExpectedScheduleSha256 },
        { "config_sha256", kExpectedConfigSha256 },
        { "fixture_sha256", kExpectedFixtureSha256 },
        { "bracket_protocol_id", DecimalU32(kBracketProtocolId) },
        { "bracket_protocol", kBracketProtocol },
        { "target_cpu", cpu.str() },
        { "online_cpus", kFrozenOnlineCpuList },
        { "logical_cpu_count", DecimalU32(kFrozenLogicalCpuCount) },
        { "physical_core_count", DecimalU32(kFrozenPhysicalCoreCount) },
        { "initial_allowed_cpus", kFrozenInitialAffinityList },
        { "target_cpu_thread_siblings", kFrozenTargetSiblingList },
        { "cpu_selection_rule", kCpuSelectionRule },
        { "cpu_selection_first_cpu", DecimalU32(kCpuSelectionFirstCpu) },
        { "cpu_selection_last_cpu", DecimalU32(kCpuSelectionLastCpu) },
        { "cpu_selection_target_interrupt_count", "289" },
        { "cpu_selection_runner_up_cpu", DecimalU32(kCpuSelectionRunnerUpCpu) },
        { "cpu_selection_runner_up_interrupt_count", "297" },
        { "cpu_selection_snapshot_path", kCpuSelectionSnapshotPath },
        { "cpu_selection_snapshot_sha256", kCpuSelectionSnapshotSha256 },
        { "output_path_hex", HexEncode(output_path) },
        { "cmake_sha256", WIREHAIR_WH2_CMAKE_SHA256 },
        { "executable_sha256", executable_sha256 },
        { "elf_build_id", executable_build_id },
        { "harmless_probe_source_sha256", kHarmlessProbeSourceSha256 },
        { "harmless_probe_binary_sha256", kHarmlessProbeBinarySha256 }
    };
    for (size_t i = 0u; i < sizeof(expected) / sizeof(expected[0]); ++i)
    {
        const std::map<std::string, std::string>::const_iterator found =
            fields.find(expected[i].Key);
        if (found == fields.end() || found->second != expected[i].Value)
        {
            diagnostic = std::string("launch manifest mismatch for ") +
                expected[i].Key;
            return false;
        }
    }
    const char* const host_fields[] = {
        "cpu_family",
        "cpu_model",
        "cpu_stepping",
        "leaf1_eax_hex",
        "leaf1_ebx_hex",
        "leaf1_ecx_hex",
        "leaf1_edx_hex",
        "leaf6_eax_hex",
        "leaf6_ebx_hex",
        "leaf6_ecx_hex",
        "leaf6_edx_hex",
        "leaf80000008_eax_hex",
        "leaf80000008_ebx_hex",
        "leaf80000008_ecx_hex",
        "leaf80000008_edx_hex",
        "leaf80000021_eax_hex",
        "leaf80000021_ebx_hex",
        "leaf80000021_ecx_hex",
        "leaf80000021_edx_hex"
    };
    for (size_t i = 0u; i < sizeof(host_fields) / sizeof(host_fields[0]); ++i)
    {
        const std::map<std::string, std::string>::const_iterator found =
            fields.find(host_fields[i]);
        if (found == fields.end() || found->second.empty())
        {
            diagnostic = std::string("launch manifest missing host field ") +
                host_fields[i];
            return false;
        }
    }
    const size_t expected_field_count =
        sizeof(expected) / sizeof(expected[0]) +
        sizeof(host_fields) / sizeof(host_fields[0]);
    if (fields.size() != expected_field_count)
    {
        diagnostic = "launch manifest contains an unknown or missing field";
        return false;
    }
    diagnostic.clear();
    return true;
}
#endif

bool ReadLaunchManifest(
    const std::string& path,
    const std::string& expected_sha256,
    int target_cpu,
    const std::string& output_path,
    const std::string& executable_sha256,
    const std::string& executable_build_id,
    LaunchManifest& manifest,
    std::string& diagnostic)
{
#if defined(__linux__)
    manifest = LaunchManifest();
    if (!IsCanonicalAbsolutePath(path) ||
        !IsCanonicalAbsolutePath(output_path) ||
        !IsLowerHexSha256(expected_sha256.c_str()) ||
        IsAllZeroHex(expected_sha256.c_str()))
    {
        diagnostic = "launch manifest path or expected SHA-256 is invalid";
        return false;
    }
    const int file = open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
    if (file < 0)
    {
        diagnostic = std::string("open launch manifest failed: ") +
            std::strerror(errno);
        return false;
    }
    struct stat before;
    if (fstat(file, &before) != 0 || !S_ISREG(before.st_mode) ||
        before.st_size <= 0 ||
        static_cast<uint64_t>(before.st_size) > kMaximumLaunchManifestBytes)
    {
        diagnostic = "launch manifest is not a bounded regular file";
        (void)close(file);
        return false;
    }
    std::string content(static_cast<size_t>(before.st_size), '\0');
    size_t offset = 0u;
    while (offset < content.size())
    {
        const ssize_t got = read(file, &content[offset], content.size() - offset);
        if (got < 0)
        {
            if (errno == EINTR) continue;
            diagnostic = std::string("read launch manifest failed: ") +
                std::strerror(errno);
            (void)close(file);
            return false;
        }
        if (got == 0)
        {
            diagnostic = "launch manifest truncated during read";
            (void)close(file);
            return false;
        }
        offset += static_cast<size_t>(got);
    }
    char extra = 0;
    if (read(file, &extra, 1u) != 0)
    {
        diagnostic = "launch manifest grew during read";
        (void)close(file);
        return false;
    }
    struct stat after;
    const bool stable = fstat(file, &after) == 0 &&
        before.st_dev == after.st_dev && before.st_ino == after.st_ino &&
        before.st_size == after.st_size &&
        before.st_mtim.tv_sec == after.st_mtim.tv_sec &&
        before.st_mtim.tv_nsec == after.st_mtim.tv_nsec;
    const bool close_ok = close(file) == 0;
    if (!stable || !close_ok)
    {
        diagnostic = "launch manifest changed during read or close failed";
        return false;
    }
    const std::string observed_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(content);
    manifest.Path = path;
    manifest.Sha256 = observed_sha256;
    manifest.Content = content;
    if (observed_sha256 != expected_sha256)
    {
        diagnostic = "launch manifest SHA-256 mismatch";
        return false;
    }
    if (!ParseLaunchManifestContent(
            content,
            target_cpu,
            output_path,
            executable_sha256,
            executable_build_id,
            manifest.Fields,
            diagnostic)) return false;
    return true;
#else
    (void)path;
    (void)expected_sha256;
    (void)target_cpu;
    (void)output_path;
    (void)executable_sha256;
    (void)executable_build_id;
    (void)manifest;
    diagnostic = "launch manifests require Linux";
    return false;
#endif
}

std::string DecimalU32(uint32_t value)
{
    std::ostringstream text;
    text << value;
    return text.str();
}

std::string HexU32(uint32_t value)
{
    std::ostringstream text;
    text << std::hex << std::nouppercase << std::setw(8)
         << std::setfill('0') << value;
    return text.str();
}

bool ValidateLaunchManifestHost(
    const LaunchManifest& manifest,
    const CpuFeatures& features,
    std::string& diagnostic)
{
    const struct HostField
    {
        const char* Key;
        std::string Value;
    } expected[] = {
        { "cpu_family", DecimalU32(features.Family) },
        { "cpu_model", DecimalU32(features.Model) },
        { "cpu_stepping", DecimalU32(features.Stepping) },
        { "leaf1_eax_hex", HexU32(features.Leaf1Eax) },
        { "leaf1_ebx_hex", HexU32(features.Leaf1Ebx) },
        { "leaf1_ecx_hex", HexU32(features.Leaf1Ecx) },
        { "leaf1_edx_hex", HexU32(features.Leaf1Edx) },
        { "leaf6_eax_hex", HexU32(features.Leaf6Eax) },
        { "leaf6_ebx_hex", HexU32(features.Leaf6Ebx) },
        { "leaf6_ecx_hex", HexU32(features.Leaf6Ecx) },
        { "leaf6_edx_hex", HexU32(features.Leaf6Edx) },
        { "leaf80000008_eax_hex", HexU32(features.Leaf80000008Eax) },
        { "leaf80000008_ebx_hex", HexU32(features.Leaf80000008Ebx) },
        { "leaf80000008_ecx_hex", HexU32(features.Leaf80000008Ecx) },
        { "leaf80000008_edx_hex", HexU32(features.Leaf80000008Edx) },
        { "leaf80000021_eax_hex", HexU32(features.Leaf80000021Eax) },
        { "leaf80000021_ebx_hex", HexU32(features.Leaf80000021Ebx) },
        { "leaf80000021_ecx_hex", HexU32(features.Leaf80000021Ecx) },
        { "leaf80000021_edx_hex", HexU32(features.Leaf80000021Edx) }
    };
    for (size_t i = 0u; i < sizeof(expected) / sizeof(expected[0]); ++i)
    {
        const std::map<std::string, std::string>::const_iterator found =
            manifest.Fields.find(expected[i].Key);
        if (found == manifest.Fields.end() || found->second != expected[i].Value)
        {
            diagnostic = std::string("launch manifest host mismatch for ") +
                expected[i].Key;
            return false;
        }
    }
    return true;
}

void AppendInitialEvidence(
    std::ostringstream& evidence,
    int target_cpu,
    const std::string& output_path,
    const std::string& launch_manifest_path,
    const std::string& launch_manifest_sha256)
{
    evidence << "PROTOCOL,tag=" << kProtocolTag << '\n'
             << "BRACKET_PROTOCOL,id=" << kBracketProtocolId
             << ",definition=" << kBracketProtocol << '\n'
             << std::setprecision(17)
             << "STATISTICAL_CONTRACT,p=2239/2240,tail=1/2240,"
                "critical_n96=" << kT95N96
             << ",critical_n576=" << kT95N576
             << ",lower_log_margin=" << kLogLowerBoundDecimal
             << ",upper_log_margin=" << kLogUpperBoundDecimal
             << ",critical_n96_bits=0x" << std::hex
             << DoubleBits(kT95N96)
             << ",critical_n576_bits=0x" << DoubleBits(kT95N576)
             << ",lower_log_margin_bits=0x"
             << DoubleBits(kLogLowerBound)
             << ",upper_log_margin_bits=0x"
             << DoubleBits(kLogUpperBound)
             << std::dec
             << ",official_gates=" << kOfficialGateCount << '\n'
             << "PROVENANCE_EXPECTED,git_commit="
             << WIREHAIR_WH2_SOURCE_GIT_COMMIT
             << ",controller_sha256="
             << WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256
             << ",codec_header_sha256="
             << WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256
             << ",codec_source_sha256="
             << WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256
             << ",cmake_sha256=" << WIREHAIR_WH2_CMAKE_SHA256
             << ",schedule_sha256=" << kExpectedScheduleSha256
             << ",config_sha256=" << kExpectedConfigSha256
             << ",fixture_sha256=" << kExpectedFixtureSha256
             << ",harmless_probe_source_sha256="
             << kHarmlessProbeSourceSha256
             << ",harmless_probe_binary_sha256="
             << kHarmlessProbeBinarySha256
             << ",frozen_target_cpu=" << kFrozenTargetCpu << '\n'
             << "TOPOLOGY_EXPECTED,online=" << kFrozenOnlineCpuList
             << ",logical_cpu_count=" << kFrozenLogicalCpuCount
             << ",physical_core_count=" << kFrozenPhysicalCoreCount
             << ",initial_affinity=" << kFrozenInitialAffinityList
             << ",target_siblings=" << kFrozenTargetSiblingList << '\n'
             << "CPU_SELECTION_EXPECTED,rule_length="
             << (sizeof(kCpuSelectionRule) - 1u)
             << ",rule_hex=" << HexEncode(kCpuSelectionRule)
             << ",first_cpu=" << kCpuSelectionFirstCpu
             << ",last_cpu=" << kCpuSelectionLastCpu
             << ",target_cpu=" << kFrozenTargetCpu
             << ",target_interrupt_count="
             << kCpuSelectionTargetInterruptCount
             << ",runner_up_cpu=" << kCpuSelectionRunnerUpCpu
             << ",runner_up_interrupt_count="
             << kCpuSelectionRunnerUpInterruptCount
             << ",snapshot_path_length="
             << (sizeof(kCpuSelectionSnapshotPath) - 1u)
             << ",snapshot_path_hex=" << HexEncode(kCpuSelectionSnapshotPath)
             << ",snapshot_sha256=" << kCpuSelectionSnapshotSha256 << '\n'
             << "HOST_EXPECTED,family=" << kFrozenCpuFamily
             << ",model=" << kFrozenCpuModel
             << ",stepping=" << kFrozenCpuStepping
             << ",leaf1_eax=" << HexU32(kFrozenLeaf1Eax)
             << ",leaf1_ebx=" << HexU32(kFrozenLeaf1Ebx)
             << ",leaf1_ecx=" << HexU32(kFrozenLeaf1Ecx)
             << ",leaf1_edx=" << HexU32(kFrozenLeaf1Edx)
             << ",leaf6_eax=" << HexU32(kFrozenLeaf6Eax)
             << ",leaf6_ebx=" << HexU32(kFrozenLeaf6Ebx)
             << ",leaf6_ecx=" << HexU32(kFrozenLeaf6Ecx)
             << ",leaf6_edx=" << HexU32(kFrozenLeaf6Edx)
             << ",leaf80000008_eax="
             << HexU32(kFrozenLeaf80000008Eax)
             << ",leaf80000008_ebx="
             << HexU32(kFrozenLeaf80000008Ebx)
             << ",leaf80000008_ecx="
             << HexU32(kFrozenLeaf80000008Ecx)
             << ",leaf80000008_edx="
             << HexU32(kFrozenLeaf80000008Edx)
             << ",leaf80000021_eax="
             << HexU32(kFrozenLeaf80000021Eax)
             << ",leaf80000021_ebx="
             << HexU32(kFrozenLeaf80000021Ebx)
             << ",leaf80000021_ecx="
             << HexU32(kFrozenLeaf80000021Ecx)
             << ",leaf80000021_edx="
             << HexU32(kFrozenLeaf80000021Edx) << '\n';
    std::ostringstream canonical_command;
    canonical_command << "--calibrate --cpu " << target_cpu
                      << " --output " << output_path
                      << " --launch-manifest " << launch_manifest_path
                      << " --launch-manifest-sha256 "
                      << launch_manifest_sha256;
    evidence << "COMMAND,target_cpu=" << target_cpu
             << ",output_path_length=" << output_path.size()
             << ",output_path_hex=" << HexEncode(output_path)
             << ",launch_manifest_path_length="
             << launch_manifest_path.size()
             << ",launch_manifest_path_hex="
             << HexEncode(launch_manifest_path)
             << ",launch_manifest_sha256=" << launch_manifest_sha256
             << ",canonical_argv_length=" << canonical_command.str().size()
             << ",canonical_argv_hex=" << HexEncode(canonical_command.str())
             << '\n';
}

void AppendLaunchManifestEvidence(
    std::ostringstream& evidence,
    const LaunchManifest& manifest,
    const char* gate)
{
    evidence << "LAUNCH_MANIFEST,path_length=" << manifest.Path.size()
             << ",path_hex=" << HexEncode(manifest.Path)
             << ",content_length=" << manifest.Content.size()
             << ",sha256=" << manifest.Sha256
             << ",content_hex=" << HexEncode(manifest.Content)
             << ",gate=" << gate << '\n';
}

void AppendLaunchManifestClaimsEvidence(
    std::ostringstream& evidence,
    const LaunchManifest& manifest)
{
    evidence << "LAUNCH_MANIFEST_CLAIMS,count=" << manifest.Fields.size();
    for (std::map<std::string, std::string>::const_iterator next =
             manifest.Fields.begin(); next != manifest.Fields.end(); ++next)
    {
        evidence << ",key_length=" << next->first.size()
                 << ",key_hex=" << HexEncode(next->first)
                 << ",value_length=" << next->second.size()
                 << ",value_hex=" << HexEncode(next->second);
    }
    evidence << '\n';
}

void AppendHostEvidence(
    std::ostringstream& evidence,
    const CpuFeatures& features,
    const char* raw_read_gate)
{
    evidence << "HOST_RAW,compiled=" << (features.Compiled ? 1u : 0u)
             << ",authentic_amd="
             << (features.VendorAuthenticAmd ? 1u : 0u)
             << ",hypervisor=" << (features.Hypervisor ? 1u : 0u)
             << ",has_leaf6=" << (features.HasLeaf6 ? 1u : 0u)
             << ",aperf_mperf=" << (features.AperfMperf ? 1u : 0u)
             << ",has_leaf80000008="
             << (features.HasLeaf80000008 ? 1u : 0u)
             << ",rdpru=" << (features.Rdpru ? 1u : 0u)
             << ",rdpru_max=" << features.RdpruMax
             << ",has_leaf80000021="
             << (features.HasLeaf80000021 ? 1u : 0u)
             << ",lfence_always_serializing="
             << (features.LfenceAlwaysSerializing ? 1u : 0u)
             << ",family=" << features.Family
             << ",model=" << features.Model
             << ",stepping=" << features.Stepping
             << ",leaf1_eax=0x" << std::hex << features.Leaf1Eax
             << ",leaf1_ebx=0x" << features.Leaf1Ebx
             << ",leaf1_ecx=0x" << features.Leaf1Ecx
             << ",leaf1_edx=0x" << features.Leaf1Edx
             << ",leaf6_eax=0x" << features.Leaf6Eax
             << ",leaf6_ebx=0x" << features.Leaf6Ebx
             << ",leaf6_ecx=0x" << features.Leaf6Ecx
             << ",leaf6_edx=0x" << features.Leaf6Edx
             << ",leaf80000008_eax=0x" << features.Leaf80000008Eax
             << ",leaf80000008_ebx=0x" << features.Leaf80000008Ebx
             << ",leaf80000008_ecx=0x" << features.Leaf80000008Ecx
             << ",leaf80000008_edx=0x" << features.Leaf80000008Edx
             << ",leaf80000021_eax=0x" << features.Leaf80000021Eax
             << ",leaf80000021_ebx=0x" << features.Leaf80000021Ebx
             << ",leaf80000021_ecx=0x" << features.Leaf80000021Ecx
             << ",leaf80000021_edx=0x" << features.Leaf80000021Edx
             << std::dec << ",raw_read_gate=" << raw_read_gate << '\n';
}

void AppendTopologyEvidence(
    std::ostringstream& evidence,
    const HostTopology& topology,
    const char* read_gate)
{
    const std::string online = FormatCpuList(topology.Online);
    const std::string affinity = FormatCpuList(topology.InitialAffinity);
    const std::string siblings = FormatCpuList(topology.TargetSiblings);
    evidence << "TOPOLOGY_OBSERVED,online_raw_length="
             << topology.OnlineRaw.size()
             << ",online_raw_hex=" << HexEncode(topology.OnlineRaw)
             << ",online_set_length=" << online.size()
             << ",online_set_hex=" << HexEncode(online)
             << ",logical_cpu_count=" << topology.Online.size()
             << ",physical_core_count=" << topology.PhysicalCoreCount
             << ",initial_affinity_length=" << affinity.size()
             << ",initial_affinity_hex=" << HexEncode(affinity)
             << ",target_siblings_raw_length="
             << topology.SiblingsRaw.size()
             << ",target_siblings_raw_hex="
             << HexEncode(topology.SiblingsRaw)
             << ",target_siblings_length=" << siblings.size()
             << ",target_siblings_hex=" << HexEncode(siblings)
             << ",read_gate=" << read_gate << '\n';
}

struct AccessReceipt
{
    uint32_t Selector = 0u;
    uint64_t First = 0u;
    uint64_t Second = 0u;
    bool FirstOk = false;
    bool SecondOk = false;
    int FirstSignal = 0;
    int SecondSignal = 0;
    bool CpuVerifyOk = false;
    std::string CpuVerifyDetail;
};

bool AccessReceiptPasses(const AccessReceipt& receipt)
{
    return receipt.Selector < 2u && receipt.FirstOk && receipt.SecondOk &&
        receipt.FirstSignal == 0 && receipt.SecondSignal == 0 &&
        receipt.Second > receipt.First && receipt.CpuVerifyOk;
}

void AppendAccessEvidence(
    std::ostringstream& evidence,
    const AccessReceipt& receipt)
{
    evidence << "ACCESS,selector=" << receipt.Selector
             << ",name=" << SelectorName(receipt.Selector)
             << ",first_ok=" << (receipt.FirstOk ? 1u : 0u)
             << ",second_ok=" << (receipt.SecondOk ? 1u : 0u)
             << ",first=" << receipt.First
             << ",second=" << receipt.Second
             << ",signal_first=" << receipt.FirstSignal
             << ",signal_second=" << receipt.SecondSignal
             << ",monotonic="
             << (receipt.Second > receipt.First ? 1u : 0u)
             << ",cpu_verify_ok=" << (receipt.CpuVerifyOk ? 1u : 0u)
             << ",cpu_verify_detail_length="
             << receipt.CpuVerifyDetail.size()
             << ",cpu_verify_detail_hex="
             << HexEncode(receipt.CpuVerifyDetail)
             << ",gate="
             << (AccessReceiptPasses(receipt) ? "PASS" : "FAIL") << '\n';
}

void AppendContextColumns(
    std::ostringstream& output,
    const ContextReceipt& before,
    const ContextReceipt& after)
{
    output << ',' << before.Cpu
           << ',' << after.Cpu
           << ',' << before.AffinityCount
           << ',' << after.AffinityCount
           << ',' << before.TargetInAffinity
           << ',' << after.TargetInAffinity
           << ',' << before.Voluntary
           << ',' << after.Voluntary
           << ',' << before.Involuntary
           << ',' << after.Involuntary;
}

void AppendCompactRows(
    const std::vector<CompactHardwareReceipt>& receipts,
    const char* record,
    std::ostringstream& output)
{
    for (size_t i = 0u; i < receipts.size(); ++i)
    {
        const CompactHardwareReceipt& receipt = receipts[i];
        output << record
               << ',' << receipt.Ordinal
               << ',' << receipt.GlobalPanelAction
               << ',' << receipt.SelectorPanelAction
               << ',' << receipt.OriginalSuperblock
               << ',' << receipt.SelectorSuperblock
               << ',' << receipt.Selector
               << ',' << SelectorName(receipt.Selector)
               << ',' << kBracketProtocolId
               << ',' << kBracketProtocol
               << ',' << receipt.Cell
               << ',' << static_cast<uint32_t>(receipt.ScopeValue)
               << ',' << ScopeName(receipt.ScopeValue)
               << ',' << BankName(receipt.Bank)
               << ',' << receipt.CounterStart
               << ',' << receipt.CounterEnd
               << ',' << receipt.CounterDelta;
        AppendContextColumns(output, receipt.Before, receipt.After);
        output << ',' << receipt.BeginReadSucceeded
               << ',' << receipt.EndReadSucceeded
               << ',' << receipt.VerificationMask0
               << ',' << receipt.VerificationMask1
               << ',' << receipt.Summary.SuccessCount
               << ',' << receipt.Summary.OverheadSum
               << ',' << receipt.Summary.OverheadMinimum
               << ',' << receipt.Summary.OverheadMaximum
               << ',' << receipt.Summary.DirectSystematicSum
               << ',' << receipt.Summary.InnerSha256 << '\n';
    }
}

void AppendEvidenceRows(
    const std::vector<EmptyReceipt>& empty,
    const std::vector<CompactHardwareReceipt>& warm,
    const std::vector<CompactHardwareReceipt>& washout,
    const std::vector<RawTiming>& raw,
    std::ostringstream& output)
{
    output << "EMPTY_SCHEMA,columns=sample;selector_position;selector;name;"
              "bracket_protocol_id;bracket_protocol;"
              "start;end;delta;begin_ok;end_ok;cpu_before;cpu_after;affinity_before;"
              "affinity_after;target_before;target_after;nvcsw_before;"
              "nvcsw_after;nivcsw_before;nivcsw_after\n";
    for (size_t i = 0u; i < empty.size(); ++i)
    {
        const EmptyReceipt& receipt = empty[i];
        output << "EMPTY," << receipt.Sample
               << ',' << receipt.SelectorPosition
               << ',' << receipt.Selector
               << ',' << SelectorName(receipt.Selector)
               << ',' << kBracketProtocolId
               << ',' << kBracketProtocol
               << ',' << receipt.Start
               << ',' << receipt.End
               << ',' << receipt.Delta
               << ',' << receipt.BeginReadSucceeded
               << ',' << receipt.EndReadSucceeded;
        AppendContextColumns(output, receipt.Before, receipt.After);
        output << '\n';
    }
    output << "COMPACT_SCHEMA,columns=ordinal;global_action;selector_action;"
              "original_superblock;selector_superblock;selector;selector_name;"
              "bracket_protocol_id;bracket_protocol;"
              "cell;scope_index;scope;bank;start;end;delta;cpu_before;cpu_after;"
              "affinity_before;affinity_after;target_before;target_after;"
              "nvcsw_before;nvcsw_after;nivcsw_before;nivcsw_after;begin_ok;"
              "end_ok;mask0;mask1;success_count;overhead_sum;overhead_min;"
              "overhead_max;direct_sum;inner_sha256\n";
    AppendCompactRows(warm, "WARM", output);
    AppendCompactRows(washout, "WASHOUT", output);
    output << "RAW_SCHEMA,columns=index;global_action;selector_action;"
              "original_superblock;selector_superblock;cell;width_index;"
              "width_bytes;scope_index;scope;root_index;root;round;visit;u;"
              "direction;mode;mode_position;orientation_flip;rotation;"
              "first_selector;selector_position;selector;selector_name;"
              "bracket_protocol_id;bracket_protocol;"
              "sequence_position;sequence;family;period;label;map_a;physical;"
              "predecessor;batch_size;start;end;delta;cpu_before;cpu_after;"
              "affinity_before;affinity_after;target_before;target_after;"
              "nvcsw_before;nvcsw_after;nivcsw_before;nivcsw_after;begin_ok;"
              "end_ok;mask0;mask1;success_count;overhead_sum;overhead_min;"
              "overhead_max;direct_sum;inner_sha256\n";
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        const RawTiming& receipt = raw[i];
        const LogicalCellInfo info = GetLogicalCellInfo(receipt.Cell);
        output << "RAW," << receipt.RawIndex
               << ',' << receipt.GlobalPanelAction
               << ',' << receipt.SelectorPanelAction
               << ',' << receipt.OriginalSuperblock
               << ',' << receipt.SelectorSuperblock
               << ',' << receipt.Cell
               << ',' << info.WidthIndex
               << ',' << info.BlockBytes
               << ',' << static_cast<uint32_t>(info.ScopeValue)
               << ',' << ScopeName(info.ScopeValue)
               << ',' << info.RootIndex
               << ",0x" << std::hex << info.Root << std::dec
               << ',' << receipt.Round
               << ',' << receipt.Visit
               << ',' << receipt.U
               << ',' << DirectionName(receipt.DirectionValue)
               << ',' << ModeName(receipt.ModeValue)
               << ',' << receipt.ModePosition
               << ',' << receipt.OrientationFlip
               << ',' << receipt.Rotation
               << ',' << receipt.FirstSelector
               << ',' << receipt.SelectorPosition
               << ',' << receipt.Selector
               << ',' << SelectorName(receipt.Selector)
               << ',' << kBracketProtocolId
               << ',' << kBracketProtocol
               << ',' << receipt.SequencePosition
               << ',' << receipt.Sequence
               << ',' << receipt.SequenceFamily
               << ',' << receipt.Period
               << ',' << LabelName(receipt.Label)
               << ',' << BankName(receipt.LogicalAPhysical)
               << ',' << BankName(receipt.Physical)
               << ',' << static_cast<uint32_t>(receipt.Predecessor)
               << ',' << kBatchSize
               << ',' << receipt.Result.Batch.CounterStart
               << ',' << receipt.Result.Batch.CounterEnd
               << ',' << receipt.Result.Batch.CounterDelta;
        AppendContextColumns(
            output, receipt.Result.Before, receipt.Result.After);
        output << ',' << (receipt.Result.Batch.BeginReadSucceeded ? 1u : 0u)
               << ',' << (receipt.Result.Batch.EndReadSucceeded ? 1u : 0u)
               << ',' << receipt.Result.Batch.VerificationMask[0]
               << ',' << receipt.Result.Batch.VerificationMask[1]
               << ',' << receipt.Result.Summary.SuccessCount
               << ',' << receipt.Result.Summary.OverheadSum
               << ',' << receipt.Result.Summary.OverheadMinimum
               << ',' << receipt.Result.Summary.OverheadMaximum
               << ',' << receipt.Result.Summary.DirectSystematicSum
               << ',' << receipt.Result.Summary.InnerSha256 << '\n';
    }
}

void AppendFailedAttempt(
    const FailedAttemptReceipt& failure,
    const std::string& diagnostic,
    std::ostringstream& output)
{
    if (!failure.Present) return;
    output << "FAILED_ATTEMPT,stage_length=" << failure.Stage.size()
           << ",stage_hex=" << HexEncode(failure.Stage)
           << ",identity_length=" << failure.Identity.size()
           << ",identity_hex=" << HexEncode(failure.Identity)
           << ",ordinal=" << failure.Ordinal
           << ",selector=" << failure.Selector
           << ",selector_name=" << SelectorName(failure.Selector)
           << ",cell=" << failure.Cell
           << ",empty_attempt=" << (failure.EmptyAttempt ? 1u : 0u)
           << ",detail_length=" << diagnostic.size()
           << ",detail_hex=" << HexEncode(diagnostic)
           << ",bracket_protocol_id=" << kBracketProtocolId
           << ",bracket_protocol=" << kBracketProtocol;
    if (failure.EmptyAttempt)
    {
        output << ",start=" << failure.Empty.Start
               << ",end=" << failure.Empty.End
               << ",delta=" << failure.Empty.Delta
               << ",begin_ok=" << failure.Empty.BeginReadSucceeded
               << ",end_ok=" << failure.Empty.EndReadSucceeded;
        AppendContextColumns(
            output, failure.Empty.Before, failure.Empty.After);
    }
    else
    {
        output << ",scope=" << ScopeName(failure.ScopeValue)
               << ",bank=" << BankName(failure.Bank)
               << ",start=" << failure.Qualified.Batch.CounterStart
               << ",end=" << failure.Qualified.Batch.CounterEnd
               << ",delta=" << failure.Qualified.Batch.CounterDelta;
        AppendContextColumns(
            output, failure.Qualified.Before, failure.Qualified.After);
        output << ",begin_ok="
               << (failure.Qualified.Batch.BeginReadSucceeded ? 1u : 0u)
               << ",end_ok="
               << (failure.Qualified.Batch.EndReadSucceeded ? 1u : 0u)
               << ",mask0=" << failure.Qualified.Batch.VerificationMask[0]
               << ",mask1=" << failure.Qualified.Batch.VerificationMask[1]
               << ",success_count="
               << failure.Qualified.Summary.SuccessCount
               << ",inner_sha256="
               << failure.Qualified.Summary.InnerSha256;
    }
    output << '\n';
}

struct RunEvidenceState
{
    RunEvidenceState(
        std::ostringstream& evidence,
        std::string& diagnostic)
        : Evidence(evidence), Diagnostic(diagnostic)
    {
    }

    ~RunEvidenceState()
    {
        try {
            Flush();
        }
        catch (...) {
        }
    }

    void Flush()
    {
        if (Flushed) return;
        std::ostringstream staged;
        staged << "RUN_STATE_COUNTS,empty=" << Empty.size()
               << ",warm=" << Warm.size()
               << ",washout=" << Washout.size()
               << ",raw=" << Raw.size()
               << ",failed_attempt=" << (Failure.Present ? 1u : 0u)
               << '\n';
        AppendFailedAttempt(Failure, Diagnostic, staged);
        AppendEvidenceRows(Empty, Warm, Washout, Raw, staged);
        staged << "RUN_STATE_FLUSH_COMPLETE=1\n";
        const std::string materialized = staged.str();
        Evidence << materialized;
        if (!Evidence) throw std::runtime_error("run-state evidence append failed");
        Flushed = true;
    }

    std::ostringstream& Evidence;
    std::string& Diagnostic;
    std::vector<EmptyReceipt> Empty;
    std::vector<CompactHardwareReceipt> Warm;
    std::vector<CompactHardwareReceipt> Washout;
    std::vector<RawTiming> Raw;
    FailedAttemptReceipt Failure;
    bool Flushed = false;
};

void AppendOfficialResults(
    const std::array<OfficialResult, kOfficialGateCount>& official,
    std::ostringstream& output)
{
    for (size_t i = 0u; i < official.size(); ++i)
    {
        const OfficialResult& result = official[i];
        output << std::setprecision(17)
               << "OFFICIAL,index=" << i
               << ",identity=" << result.Identity
               << ",n=" << result.Stats.Count
               << ",critical=" << CriticalValue(result.Stats.Count)
               << ",mean_log=" << result.Stats.Mean
               << ",lower_log=" << result.Stats.Lower
               << ",upper_log=" << result.Stats.Upper
               << ",gm=" << result.Stats.GeometricMean
               << ",lower_ratio=" << result.Stats.GeometricLower
               << ",upper_ratio=" << result.Stats.GeometricUpper
               << ",gate=" << (result.Pass ? "PASS" : "FAIL") << '\n';
        if (!result.Pass) {
            output << "FAILED_OFFICIAL,index=" << i
                   << ",identity=" << result.Identity << '\n';
        }
    }
}

std::string PublicationReceiptLine(
    const char* publication_gate,
    const std::string& output_path,
    const std::string& payload)
{
    std::ostringstream receipt;
    receipt << "EVIDENCE_PUBLICATION=" << publication_gate
            << ",output_path_length=" << output_path.size()
            << ",output_path_hex=" << HexEncode(output_path)
            << ",payload_length=" << payload.size()
            << ",payload_sha256="
            << wirehair::wh2_benchmark::Sha256Hex(payload) << '\n';
    return receipt.str();
}

bool EmitPublicationAuthority(
    std::ostream& output,
    const std::string& publication_receipt,
    bool scientific_valid,
    bool overall_valid)
{
    output << publication_receipt
           << "SCIENTIFIC_GATE_VERDICT="
           << (scientific_valid ? "PASS" : "FAIL") << '\n'
           << "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION="
           << (overall_valid ? "VALID" : "INVALID") << '\n';
    output.flush();
    return output.good();
}

int PublishCalibration(
    const std::string& output_path,
    std::ostringstream& evidence,
    bool valid,
    const std::string& diagnostic)
{
    // Never mutate the shared accumulator with a terminal verdict.  If
    // materialization or publication throws, the boundary can safely retry an
    // INVALID payload without inheriting an earlier VALID/PENDING terminal.
    const std::string payload = FinalizeCalibrationPayload(
        evidence.str(), valid, diagnostic);
    std::string publication_diagnostic;
    if (!AtomicPublishNew(output_path, payload, publication_diagnostic))
    {
        std::cerr << "Wirehair1 native RDPRU batch128 A/A evidence publication "
                  << "failed: " << publication_diagnostic << '\n';
        (void)EmitPublicationAuthority(
            std::cout,
            PublicationReceiptLine("FAIL", output_path, payload),
            valid,
            false);
        return 1;
    }
    if (!EmitPublicationAuthority(
            std::cout,
            PublicationReceiptLine("PASS", output_path, payload),
            valid,
            valid))
    {
        std::cerr << "EVIDENCE_PUBLICATION_AUTHORITY=FAIL\n"
                  << "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION=INVALID\n";
        std::cerr.flush();
        return 1;
    }
    return valid ? 0 : 1;
}

int RunCalibrationImpl(
    int target_cpu,
    const std::string& output_path,
    const std::string& launch_manifest_path,
    const std::string& launch_manifest_sha256,
    std::ostringstream& evidence)
{
    std::string diagnostic;
    if (target_cpu != kFrozenTargetCpu)
    {
        return PublishCalibration(
            output_path,
            evidence,
            false,
            "target CPU differs from frozen CPU 50");
    }
    if (!CheckOutputDestinationNew(output_path, diagnostic)) {
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    if (!EmbeddedIdentitiesValid())
    {
        evidence << "EMBEDDED_IDENTITY_GATE=FAIL\n";
        return PublishCalibration(
            output_path,
            evidence,
            false,
            "embedded commit/source identities are unknown or placeholders");
    }
    evidence << "EMBEDDED_IDENTITY_GATE=PASS\n";
    std::string executable_sha256;
    std::string executable_build_id;
    if (!ReadSelfExecutableIdentity(
            executable_sha256, executable_build_id, diagnostic))
    {
        evidence << "EXECUTABLE_OBSERVED,sha256_length="
                 << executable_sha256.size()
                 << ",sha256_hex=" << HexEncode(executable_sha256)
                 << ",build_id_length=" << executable_build_id.size()
                 << ",build_id_hex=" << HexEncode(executable_build_id)
                 << ",gate=FAIL\n";
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    evidence << "EXECUTABLE_OBSERVED,sha256=" << executable_sha256
             << ",build_id=" << executable_build_id
             << ",gate=pending-manifest-claim\n";
    HostTopology topology;
    if (!ReadHostTopology(topology, diagnostic))
    {
        AppendTopologyEvidence(evidence, topology, "FAIL");
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    AppendTopologyEvidence(evidence, topology, "PASS");
    if (!ValidateHostTopology(topology, diagnostic))
    {
        evidence << "TOPOLOGY_GATE=FAIL\n";
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    evidence << "TOPOLOGY_GATE=PASS\n";
    LaunchManifest launch_manifest;
    if (!ReadLaunchManifest(
            launch_manifest_path,
            launch_manifest_sha256,
            target_cpu,
            output_path,
            executable_sha256,
            executable_build_id,
            launch_manifest,
            diagnostic))
    {
        if (!launch_manifest.Content.empty()) {
            AppendLaunchManifestEvidence(evidence, launch_manifest, "FAIL");
        }
        if (!launch_manifest.Fields.empty()) {
            AppendLaunchManifestClaimsEvidence(evidence, launch_manifest);
        }
        evidence << "EXECUTABLE_MANIFEST_GATE=FAIL\n";
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    AppendLaunchManifestEvidence(evidence, launch_manifest, "PASS");
    AppendLaunchManifestClaimsEvidence(evidence, launch_manifest);
    evidence << "EXECUTABLE_MANIFEST_GATE=PASS\n";
    if (!PinCpuOnce(target_cpu, diagnostic)) {
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    CpuFeatures features;
    if (!ReadCpuFeaturesRaw(features, diagnostic))
    {
        AppendHostEvidence(evidence, features, "FAIL");
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    AppendHostEvidence(evidence, features, "PASS");
    if (!ValidateCpuFeatures(features, diagnostic))
    {
        evidence << "CPU_QUALIFICATION_GATE=FAIL\n";
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    evidence << "CPU_QUALIFICATION_GATE=PASS\n";
    if (!ValidateLaunchManifestHost(launch_manifest, features, diagnostic))
    {
        evidence << "HOST_MANIFEST_GATE=FAIL\n";
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    evidence << "HOST_MANIFEST_GATE=PASS\n";
    for (uint32_t selector = 0u; selector < 2u; ++selector)
    {
        AccessReceipt access;
        access.Selector = selector;
        access.FirstOk = GuardedCounterRead(
            ReadRdpru,
            nullptr,
            selector,
            access.First,
            access.FirstSignal);
        access.SecondOk = GuardedCounterRead(
            ReadRdpru,
            nullptr,
            selector,
            access.Second,
            access.SecondSignal);
        access.CpuVerifyOk = VerifyCpu(
            target_cpu, access.CpuVerifyDetail);
        AppendAccessEvidence(evidence, access);
        if (!AccessReceiptPasses(access))
        {
            std::ostringstream message;
            message << "guarded RDPRU access failed selector=" << selector
                    << ",first_ok=" << (access.FirstOk ? 1u : 0u)
                    << ",second_ok=" << (access.SecondOk ? 1u : 0u)
                    << ",signals=" << access.FirstSignal << ':'
                    << access.SecondSignal
                    << ",monotonic="
                    << (access.Second > access.First ? 1u : 0u)
                    << ",cpu_verify_ok="
                    << (access.CpuVerifyOk ? 1u : 0u);
            if (!access.CpuVerifyDetail.empty()) {
                message << ",detail=" << access.CpuVerifyDetail;
            }
            return PublishCalibration(
                output_path, evidence, false, message.str());
        }
    }

    FrozenSchedule schedule;
    if (!BuildFrozenSchedule(schedule)) {
        return PublishCalibration(
            output_path, evidence, false, "schedule construction failed");
    }
    const std::string schedule_sha256 = ScheduleSha256(schedule);
    const std::string config_sha256 = ConfigSha256(schedule_sha256);
    evidence << "FREEZE_EXPECTED,schedule_sha256="
             << kExpectedScheduleSha256
             << ",config_sha256=" << kExpectedConfigSha256 << '\n'
             << "FREEZE_OBSERVED,schedule_sha256=" << schedule_sha256
             << ",config_sha256=" << config_sha256 << '\n';
    if (schedule_sha256 != kExpectedScheduleSha256 ||
        config_sha256 != kExpectedConfigSha256)
    {
        evidence << "FREEZE_GATE=FAIL\n";
        return PublishCalibration(
            output_path, evidence, false, "frozen schedule/config hash mismatch");
    }
    evidence << "FREEZE_GATE=PASS\n";
    if (wirehair_init() != Wirehair_Success) {
        return PublishCalibration(
            output_path, evidence, false, "wirehair_init failed");
    }
    CellArray cells;
    if (!BuildCells(cells, diagnostic) ||
        !PreflightCells(cells, diagnostic))
    {
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    const std::string fixture_sha256 = FixtureSha256(cells);
    evidence << "FIXTURE,expected_sha256=" << kExpectedFixtureSha256
             << ",observed_sha256=" << fixture_sha256
             << ",preflight_batches=" << kPreflightBatchCount
             << ",preflight_inner=" << kPreflightInnerInvocationCount
             << ",gate="
             << (fixture_sha256 == kExpectedFixtureSha256 ? "PASS" : "FAIL")
             << '\n';
    if (fixture_sha256 != kExpectedFixtureSha256) {
        return PublishCalibration(
            output_path, evidence, false, "frozen fixture hash mismatch");
    }

    RunEvidenceState run_state(evidence, diagnostic);
    if (!RunPinnedWarmups(
            cells,
            target_cpu,
            run_state.Warm,
            run_state.Failure,
            diagnostic))
    {
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    if (!RunEmptySamples(
            target_cpu,
            run_state.Empty,
            run_state.Failure,
            diagnostic))
    {
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    if (!ValidateEmptyReceipts(run_state.Empty, target_cpu, diagnostic))
    {
        RecordStageFailure(
            run_state.Failure,
            "empty-validation",
            "empty-receipt-validation",
            static_cast<uint64_t>(run_state.Empty.size()));
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    std::vector<SuperblockObservation> observations;
    if (!RunFrozenPanel(
            cells,
            schedule,
            target_cpu,
            run_state.Raw,
            run_state.Washout,
            observations,
            run_state.Failure,
            diagnostic))
    {
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    if (!ValidateRawReceipts(
            schedule,
            run_state.Raw,
            observations,
            target_cpu,
            diagnostic))
    {
        RecordStageFailure(
            run_state.Failure,
            "raw-validation",
            "raw-receipt-validation",
            static_cast<uint64_t>(run_state.Raw.size()));
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    if (!ValidateCompactReceipts(
            schedule,
            run_state.Warm,
            run_state.Washout,
            run_state.Raw,
            target_cpu,
            diagnostic))
    {
        RecordStageFailure(
            run_state.Failure,
            "compact-validation",
            "warm-washout-compact-validation",
            static_cast<uint64_t>(
                run_state.Warm.size() + run_state.Washout.size()));
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }

    std::vector<size_t> by_key;
    std::array<OfficialResult, kOfficialGateCount> official;
    bool official_pass = false;
    std::array<EmptyGate, 2> empty_gates;
    if (!BuildObservationIndex(observations, by_key))
    {
        diagnostic = "post-panel observation index failed";
        RecordStageFailure(
            run_state.Failure,
            "statistics-index",
            "observation-index",
            static_cast<uint64_t>(observations.size()));
        run_state.Flush();
        return PublishCalibration(
            output_path, evidence, false, diagnostic);
    }
    if (!BuildOfficialResults(
            observations, by_key, official, official_pass))
    {
        diagnostic = "post-panel official family failed";
        RecordStageFailure(
            run_state.Failure,
            "statistics-official",
            "official-family",
            static_cast<uint64_t>(observations.size()));
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    if (!BuildEmptyGates(run_state.Empty, run_state.Raw, empty_gates))
    {
        diagnostic = "post-panel empty gate failed";
        RecordStageFailure(
            run_state.Failure,
            "statistics-empty",
            "empty-gates",
            static_cast<uint64_t>(run_state.Empty.size()));
        run_state.Flush();
        return PublishCalibration(output_path, evidence, false, diagnostic);
    }
    const std::string empty_sha256 = EmptySha256(run_state.Empty);
    const std::string compact_sha256 = CompactSha256(
        run_state.Warm, run_state.Washout);
    const std::string raw_sha256 = RawSha256(run_state.Raw);
    if (!IsLowerHexSha256(empty_sha256.c_str()) ||
        !IsLowerHexSha256(compact_sha256.c_str()) ||
        !IsLowerHexSha256(raw_sha256.c_str()))
    {
        diagnostic = "raw evidence hashing failed";
        RecordStageFailure(
            run_state.Failure,
            "evidence-hash",
            "empty-compact-raw-hash",
            static_cast<uint64_t>(run_state.Raw.size()));
        run_state.Flush();
        return PublishCalibration(
            output_path, evidence, false, diagnostic);
    }
    evidence << "COUNTS,original_superblocks=" << kOriginalSuperblockCount
             << ",selector_superblocks=" << kSelectorSuperblockCount
             << ",measured_batches=" << run_state.Raw.size()
             << ",measured_inner=" << kMeasuredInnerInvocationCount
             << ",washout_batches=" << run_state.Washout.size()
             << ",washout_inner=" << kWashoutInnerInvocationCount
             << ",warm_batches=" << run_state.Warm.size()
             << ",warm_inner=" << kWarmInnerInvocationCount
             << ",empty_per_selector=" << kEmptySampleCount << '\n'
             << "HASH,empty_sha256=" << empty_sha256
             << ",compact_sha256=" << compact_sha256
             << ",raw_sha256=" << raw_sha256 << '\n';
    bool empty_pass = true;
    for (uint32_t selector = 0u; selector < 2u; ++selector)
    {
        const EmptyGate& gate = empty_gates[selector];
        empty_pass = empty_pass && gate.Pass;
        evidence << "EMPTY_GATE,selector=" << SelectorName(selector)
                 << ",n=" << kEmptySampleCount
                 << ",nearest_rank_one_based_index="
                 << kEmptyP999OneBasedIndex
                 << ",minimum=" << gate.Minimum
                 << ",median=" << gate.Median
                 << ",p999=" << gate.P999
                 << ",maximum=" << gate.Maximum
                 << ",minimum_measured=" << gate.MinimumMeasured
                 << ",multiplier=" << kEmptyOverheadMultiplier << ",gate="
                 << (gate.Pass ? "PASS" : "FAIL") << '\n';
        if (!gate.Pass) {
            evidence << "FAILED_EMPTY_GATE,selector="
                     << SelectorName(selector) << '\n';
        }
    }
    AppendOfficialResults(official, evidence);
    if (!AppendDiagnostics(
            observations, by_key, run_state.Raw, evidence))
    {
        diagnostic = "mandatory diagnostics failed";
        RecordStageFailure(
            run_state.Failure,
            "diagnostics",
            "mandatory-diagnostics",
            static_cast<uint64_t>(observations.size()));
        run_state.Flush();
        return PublishCalibration(
            output_path, evidence, false, diagnostic);
    }
    run_state.Flush();
    const bool valid = official_pass && empty_pass;
    return PublishCalibration(
        output_path,
        evidence,
        valid,
        valid ? "all frozen gates passed" : "one or more frozen gates failed");
}

int RunCalibrationBoundary(
    int target_cpu,
    const std::string& output_path,
    const std::string& launch_manifest_path,
    const std::string& launch_manifest_sha256,
    uint32_t injected_exception)
{
    std::ostringstream evidence;
    try
    {
        AppendInitialEvidence(
            evidence,
            target_cpu,
            output_path,
            launch_manifest_path,
            launch_manifest_sha256);
        if (injected_exception == 1u) {
            throw std::runtime_error("injected standard exception before accumulation");
        }
        if (injected_exception == 2u)
        {
            evidence << "PARTIAL_RECEIPT,identity="
                        "injected-after-accumulation,ordinal=17\n";
            throw 17;
        }
        if (injected_exception == 3u)
        {
            evidence << "PARTIAL_RECEIPT,identity="
                        "injected-during-publication,ordinal=23\n";
            // Materialize a would-be PASS payload without changing the
            // accumulator, then emulate a publication-layer exception.
            const std::string pending = FinalizeCalibrationPayload(
                evidence.str(), true, "injected pending payload");
            if (pending.find("=PENDING_PUBLICATION\n") ==
                    std::string::npos) throw 23;
            throw std::runtime_error("injected publication-time exception");
        }
        return RunCalibrationImpl(
            target_cpu,
            output_path,
            launch_manifest_path,
            launch_manifest_sha256,
            evidence);
    }
    catch (const std::exception& exception)
    {
        try
        {
            return PublishCalibration(
                output_path,
                evidence,
                false,
                std::string("unexpected standard exception: ") +
                    exception.what());
        }
        catch (...) {}
        std::cerr << "native RDPRU batch128 A/A exception publication failed\n";
    }
    catch (...)
    {
        try
        {
            return PublishCalibration(
                output_path,
                evidence,
                false,
                "unexpected non-standard exception");
        }
        catch (...) {}
        std::cerr << "native RDPRU batch128 A/A exception publication failed\n";
    }
    std::cout << "EVIDENCE_PUBLICATION=FAIL\n"
              << "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION=INVALID\n";
    return 1;
}

bool TestScheduleBalance()
{
    FrozenSchedule schedule;
    if (!BuildFrozenSchedule(schedule)) return false;
    uint32_t first_counts[kCellCount][2][2] = {};
    uint32_t orientation_counts[kCellCount][2][2][2] = {};
    uint32_t position_counts[kCellCount][2][2][2] = {};
    uint32_t direction_counts[kCellCount][2][2][2] = {};
    uint32_t rotation_counts[kCellCount][2][4][2] = {};
    uint32_t parity_counts[kCellCount][2][5][2][2] = {};
    std::array<std::array<bool, 32>, kCellCount> seen_u = {{}};
    size_t next = 0u;
    for (uint32_t round = 0u; round < kRoundCount; ++round)
    {
        std::array<bool, kCellCount> seen_cells = {{ false }};
        for (uint32_t visit = 0u; visit < kCellCount; ++visit)
        {
            const uint32_t cell = CellForVisit(round, visit);
            if (cell >= kCellCount || seen_cells[cell]) return false;
            seen_cells[cell] = true;
            const uint32_t u = ScheduleU(cell, round);
            if (u >= 32u || seen_u[cell][u]) return false;
            seen_u[cell][u] = true;
            for (uint32_t mode_position = 0u; mode_position < 2u;
                 ++mode_position)
            {
                if (next >= schedule.size()) return false;
                const SuperblockPlan& plan = schedule[next++];
                const uint32_t mode = static_cast<uint32_t>(plan.ModeValue);
                const uint32_t direction =
                    static_cast<uint32_t>(plan.DirectionValue);
                const uint32_t expected_first =
                    (u ^ (u >> 1u) ^ (u >> 2u) ^ (u >> 3u) ^ (u >> 4u)) & 1u;
                if (plan.Cell != cell || plan.Round != round ||
                    plan.Visit != visit || plan.U != u || mode >= 2u ||
                    direction >= 2u || plan.ModePosition != mode_position ||
                    plan.OrientationFlip >= 2u || plan.Rotation >= 4u ||
                    plan.FirstSelector != expected_first)
                {
                    return false;
                }
                ++first_counts[cell][mode][plan.FirstSelector];
                ++orientation_counts[cell][mode][plan.OrientationFlip]
                    [plan.FirstSelector];
                ++position_counts[cell][mode][mode_position]
                    [plan.FirstSelector];
                ++direction_counts[cell][mode][direction]
                    [plan.FirstSelector];
                ++rotation_counts[cell][mode][plan.Rotation]
                    [plan.FirstSelector];
                for (uint32_t bit = 0u; bit < 5u; ++bit) {
                    ++parity_counts[cell][mode][bit][(u >> bit) & 1u]
                        [plan.FirstSelector];
                }
            }
        }
        if (std::find(seen_cells.begin(), seen_cells.end(), false) !=
            seen_cells.end()) return false;
    }
    if (next != schedule.size()) return false;
    for (uint32_t cell = 0u; cell < kCellCount; ++cell)
    {
        if (std::find(seen_u[cell].begin(), seen_u[cell].end(), false) !=
            seen_u[cell].end()) return false;
        for (uint32_t mode = 0u; mode < 2u; ++mode)
        {
            for (uint32_t selector = 0u; selector < 2u; ++selector)
            {
                if (first_counts[cell][mode][selector] != 16u) return false;
                for (uint32_t factor = 0u; factor < 2u; ++factor)
                {
                    if (orientation_counts[cell][mode][factor][selector] != 8u ||
                        position_counts[cell][mode][factor][selector] != 8u ||
                        direction_counts[cell][mode][factor][selector] != 8u)
                    {
                        return false;
                    }
                }
                for (uint32_t rotation = 0u; rotation < 4u; ++rotation) {
                    if (rotation_counts[cell][mode][rotation][selector] != 4u) {
                        return false;
                    }
                }
                for (uint32_t bit = 0u; bit < 5u; ++bit) {
                    for (uint32_t value = 0u; value < 2u; ++value) {
                        if (parity_counts[cell][mode][bit][value][selector] != 8u) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    const std::string schedule_sha256 = ScheduleSha256(schedule);
    const std::string config_sha256 = ConfigSha256(schedule_sha256);
    std::cout << "SELFTEST_FREEZE,schedule_sha256=" << schedule_sha256
              << ",config_sha256=" << config_sha256 << '\n';
    return schedule_sha256 == kExpectedScheduleSha256 &&
        config_sha256 == kExpectedConfigSha256 &&
        kOriginalSuperblockCount == 1152u &&
        kSelectorSuperblockCount == 2304u &&
        kMeasuredBatchCountPerSelector == 18432u &&
        kMeasuredBatchCount == 36864u &&
        kMeasuredInnerInvocationCount == 4718592u &&
        kWashoutBatchCount == 36864u &&
        kWarmBatchCount == 108u && kPreflightBatchCount == 54u;
}

bool TestSequenceSurface()
{
    std::array<uint32_t, 4> transitions = {{ 0u, 0u, 0u, 0u }};
    for (uint32_t period = 0u; period < 4u; ++period)
    {
        uint32_t a_count = 0u;
        for (uint32_t sequence = 0u; sequence < 4u; ++sequence)
        {
            a_count += kSequences[sequence][period] == LogicalLabel::A ? 1u : 0u;
            if (period > 0u)
            {
                const uint32_t previous = static_cast<uint32_t>(
                    kSequences[sequence][period - 1u]);
                const uint32_t current = static_cast<uint32_t>(
                    kSequences[sequence][period]);
                ++transitions[previous * 2u + current];
            }
        }
        if (a_count != 2u) return false;
    }
    for (size_t i = 0u; i < transitions.size(); ++i) {
        if (transitions[i] != 3u) return false;
    }
    return true;
}

bool TestCanonicalRatioSigns()
{
    SuperblockObservation observation;
    observation.LogicalSums = {{ 200u, 100u }};
    observation.PhysicalSums = {{ 100u, 200u }};
    observation.FactorialSums = {{ 80u, 120u, 20u, 80u }};
    observation.LogicalCounts.fill(8u);
    observation.PhysicalCounts.fill(8u);
    observation.FactorialCounts.fill(4u);
    if (!FinalizeObservation(observation) ||
        std::fabs(observation.LogLogicalAa - std::log(2.0)) > 1e-15 ||
        std::fabs(observation.LogPhysicalXy + std::log(2.0)) > 1e-15 ||
        std::fabs(observation.LogInteraction - std::log(8.0 / 3.0)) > 1e-15 ||
        observation.LogLogicalAa <= 0.0 || observation.LogPhysicalXy >= 0.0 ||
        observation.LogInteraction <= 0.0)
    {
        return false;
    }
    SuperblockObservation permuted = observation;
    std::swap(permuted.LogicalSums[0], permuted.LogicalSums[1]);
    std::swap(permuted.FactorialSums[0], permuted.FactorialSums[2]);
    std::swap(permuted.FactorialSums[1], permuted.FactorialSums[3]);
    if (!FinalizeObservation(permuted) ||
        std::fabs(permuted.LogLogicalAa + std::log(2.0)) > 1e-15 ||
        std::fabs(permuted.LogPhysicalXy + std::log(2.0)) > 1e-15 ||
        std::fabs(permuted.LogInteraction + std::log(8.0 / 3.0)) > 1e-15)
    {
        return false;
    }
    return true;
}

CpuFeatures FrozenCpuFeaturesForTest()
{
    CpuFeatures features;
    features.Compiled = true;
    features.VendorAuthenticAmd = true;
    features.HasLeaf6 = true;
    features.AperfMperf = true;
    features.HasLeaf80000008 = true;
    features.Rdpru = true;
    features.RdpruMax = 1u;
    features.HasLeaf80000021 = true;
    features.LfenceAlwaysSerializing = true;
    features.Family = kFrozenCpuFamily;
    features.Model = kFrozenCpuModel;
    features.Stepping = kFrozenCpuStepping;
    features.Leaf1Eax = kFrozenLeaf1Eax;
    features.Leaf1Ebx = kFrozenLeaf1Ebx;
    features.Leaf1Ecx = kFrozenLeaf1Ecx;
    features.Leaf1Edx = kFrozenLeaf1Edx;
    features.Leaf6Eax = kFrozenLeaf6Eax;
    features.Leaf6Ebx = kFrozenLeaf6Ebx;
    features.Leaf6Ecx = kFrozenLeaf6Ecx;
    features.Leaf6Edx = kFrozenLeaf6Edx;
    features.Leaf80000008Eax = kFrozenLeaf80000008Eax;
    features.Leaf80000008Ebx = kFrozenLeaf80000008Ebx;
    features.Leaf80000008Ecx = kFrozenLeaf80000008Ecx;
    features.Leaf80000008Edx = kFrozenLeaf80000008Edx;
    features.Leaf80000021Eax = kFrozenLeaf80000021Eax;
    features.Leaf80000021Ebx = kFrozenLeaf80000021Ebx;
    features.Leaf80000021Ecx = kFrozenLeaf80000021Ecx;
    features.Leaf80000021Edx = kFrozenLeaf80000021Edx;
    return features;
}

bool CpuRejectionPreservesEvidence(
    const CpuFeatures& features,
    const char* expected_identity)
{
    std::ostringstream evidence;
    AppendHostEvidence(evidence, features, "PASS");
    std::string diagnostic;
    const bool accepted = ValidateCpuFeatures(features, diagnostic);
    evidence << "CPU_QUALIFICATION_GATE="
             << (accepted ? "PASS" : "FAIL") << '\n';
    return !accepted && !diagnostic.empty() &&
        evidence.str().find("HOST_RAW,") != std::string::npos &&
        evidence.str().find(expected_identity) != std::string::npos &&
        evidence.str().find("CPU_QUALIFICATION_GATE=FAIL") !=
            std::string::npos;
}

bool TestCpuQualification()
{
    CpuFeatures valid = FrozenCpuFeaturesForTest();
    std::string diagnostic;
    if (!ValidateCpuFeatures(valid, diagnostic)) return false;
    CpuFeatures changed = valid;
    changed.Compiled = false;
    if (!CpuRejectionPreservesEvidence(changed, "compiled=0")) return false;
    changed = valid; changed.VendorAuthenticAmd = false;
    if (!CpuRejectionPreservesEvidence(changed, "authentic_amd=0")) return false;
    changed = valid; changed.Hypervisor = true;
    if (!CpuRejectionPreservesEvidence(changed, "hypervisor=1")) return false;
    changed = valid; changed.HasLeaf6 = false;
    if (!CpuRejectionPreservesEvidence(changed, "has_leaf6=0")) return false;
    changed = valid; changed.AperfMperf = false;
    if (!CpuRejectionPreservesEvidence(changed, "aperf_mperf=0")) return false;
    changed = valid; changed.HasLeaf80000008 = false;
    if (!CpuRejectionPreservesEvidence(
            changed, "has_leaf80000008=0")) return false;
    changed = valid; changed.Rdpru = false;
    if (!CpuRejectionPreservesEvidence(changed, "rdpru=0")) return false;
    changed = valid; changed.RdpruMax = 0u;
    if (!CpuRejectionPreservesEvidence(changed, "rdpru_max=0")) return false;
    changed = valid; changed.HasLeaf80000021 = false;
    if (!CpuRejectionPreservesEvidence(
            changed, "has_leaf80000021=0")) return false;
    changed = valid; changed.LfenceAlwaysSerializing = false;
    if (!CpuRejectionPreservesEvidence(
            changed, "lfence_always_serializing=0")) return false;
    uint32_t CpuFeatures::* const frozen_fields[] = {
        &CpuFeatures::Family,
        &CpuFeatures::Model,
        &CpuFeatures::Stepping,
        &CpuFeatures::Leaf1Eax,
        &CpuFeatures::Leaf1Ebx,
        &CpuFeatures::Leaf1Ecx,
        &CpuFeatures::Leaf1Edx,
        &CpuFeatures::Leaf6Eax,
        &CpuFeatures::Leaf6Ebx,
        &CpuFeatures::Leaf6Ecx,
        &CpuFeatures::Leaf6Edx,
        &CpuFeatures::Leaf80000008Eax,
        &CpuFeatures::Leaf80000008Ebx,
        &CpuFeatures::Leaf80000008Ecx,
        &CpuFeatures::Leaf80000008Edx,
        &CpuFeatures::Leaf80000021Eax,
        &CpuFeatures::Leaf80000021Ebx,
        &CpuFeatures::Leaf80000021Ecx,
        &CpuFeatures::Leaf80000021Edx
    };
    const char* const frozen_field_receipts[] = {
        "family=", "model=", "stepping=",
        "leaf1_eax=0x", "leaf1_ebx=0x", "leaf1_ecx=0x", "leaf1_edx=0x",
        "leaf6_eax=0x", "leaf6_ebx=0x", "leaf6_ecx=0x", "leaf6_edx=0x",
        "leaf80000008_eax=0x", "leaf80000008_ebx=0x",
        "leaf80000008_ecx=0x", "leaf80000008_edx=0x",
        "leaf80000021_eax=0x", "leaf80000021_ebx=0x",
        "leaf80000021_ecx=0x", "leaf80000021_edx=0x"
    };
    static_assert(
        sizeof(frozen_fields) / sizeof(frozen_fields[0]) ==
            sizeof(frozen_field_receipts) / sizeof(frozen_field_receipts[0]),
        "frozen CPUID mutation receipt mismatch");
    for (size_t i = 0u;
         i < sizeof(frozen_fields) / sizeof(frozen_fields[0]);
         ++i)
    {
        changed = valid;
        changed.*frozen_fields[i] ^= 1u;
        if (!CpuRejectionPreservesEvidence(
                changed, frozen_field_receipts[i])) return false;
    }
    uint32_t family = 0u;
    uint32_t model = 0u;
    uint32_t stepping = 0u;
    const uint32_t host_style =
        1u | (8u << 4u) | (15u << 8u) | (8u << 16u) | (11u << 20u);
    DecodeDisplayFamilyModel(host_style, family, model, stepping);
    return family == 26u && model == 136u && stepping == 1u;
}

bool TestHostTopologyQualification()
{
    HostTopology topology;
    topology.OnlineRaw = "0-127\n";
    topology.SiblingsRaw = "50,114\n";
    for (int cpu = kFrozenFirstOnlineCpu;
         cpu <= kFrozenLastOnlineCpu;
         ++cpu)
    {
        topology.Online.push_back(cpu);
        topology.InitialAffinity.push_back(cpu);
    }
    topology.TargetSiblings.push_back(kFrozenTargetCpu);
    topology.TargetSiblings.push_back(kFrozenTargetSiblingCpu);
    topology.PhysicalCoreCount = kFrozenPhysicalCoreCount;
    std::string diagnostic;
    if (!ValidateHostTopology(topology, diagnostic)) return false;
    std::ostringstream evidence;
    AppendTopologyEvidence(evidence, topology, "PASS");
    if (evidence.str().find("online_set_hex=302d313237") ==
            std::string::npos ||
        evidence.str().find("target_siblings_hex=35302c313134") ==
            std::string::npos) return false;
    HostTopology changed = topology;
    changed.Online.pop_back();
    if (ValidateHostTopology(changed, diagnostic)) return false;
    changed = topology;
    changed.InitialAffinity[0] = 1;
    if (ValidateHostTopology(changed, diagnostic)) return false;
    changed = topology;
    changed.PhysicalCoreCount = kFrozenPhysicalCoreCount - 1u;
    if (ValidateHostTopology(changed, diagnostic)) return false;
    changed = topology;
    changed.TargetSiblings[1] = kFrozenTargetSiblingCpu + 1;
    if (ValidateHostTopology(changed, diagnostic)) return false;
    std::vector<int> parsed;
    return ParseCpuList("0-127\n", parsed) && parsed == topology.Online &&
        ParseCpuList("50,114\n", parsed) &&
        parsed == topology.TargetSiblings &&
        !ParseCpuList("0-50,50-127", parsed) &&
        !ParseCpuList("0-127,", parsed) &&
        !ParseCpuList("+50", parsed) &&
        !ParseCpuList("050", parsed);
}

bool TestAccessReceiptFailures()
{
    AccessReceipt valid;
    valid.Selector = 0u;
    valid.First = 100u;
    valid.Second = 200u;
    valid.FirstOk = true;
    valid.SecondOk = true;
    valid.CpuVerifyOk = true;
    if (!AccessReceiptPasses(valid)) return false;
    std::ostringstream pass_evidence;
    AppendAccessEvidence(pass_evidence, valid);
    if (pass_evidence.str().find("gate=PASS") == std::string::npos) {
        return false;
    }
    AccessReceipt failures[] = {
        valid, valid, valid, valid, valid, valid
    };
    failures[0].FirstOk = false;
    failures[1].SecondOk = false;
    failures[2].FirstSignal = 4;
    failures[3].SecondSignal = 11;
    failures[4].Second = failures[4].First;
    failures[5].CpuVerifyOk = false;
    failures[5].CpuVerifyDetail = "injected verify failure";
    for (size_t i = 0u; i < sizeof(failures) / sizeof(failures[0]); ++i)
    {
        if (AccessReceiptPasses(failures[i])) return false;
        std::ostringstream evidence;
        AppendAccessEvidence(evidence, failures[i]);
        if (evidence.str().find("gate=FAIL") == std::string::npos ||
            evidence.str().find("first=") == std::string::npos ||
            evidence.str().find("second=") == std::string::npos ||
            evidence.str().find("cpu_verify_ok=") == std::string::npos)
        {
            return false;
        }
    }
    return true;
}

#if WIREHAIR_RDPRU_COMPILED
bool ContainsOpcodeSequence(
    const unsigned char* code,
    size_t size,
    const unsigned char* pattern,
    size_t pattern_size,
    size_t& position)
{
    if (!code || !pattern || pattern_size == 0u || pattern_size > size) {
        return false;
    }
    for (size_t i = position; i + pattern_size <= size; ++i)
    {
        if (std::memcmp(code + i, pattern, pattern_size) == 0)
        {
            position = i + pattern_size;
            return true;
        }
    }
    return false;
}
#endif

bool TestReaderOpcodeShape()
{
#if WIREHAIR_RDPRU_COMPILED
    const uintptr_t begin_address = reinterpret_cast<uintptr_t>(
        wh2_rdpru_bracket_opcode_begin);
    const uintptr_t end_address = reinterpret_cast<uintptr_t>(
        wh2_rdpru_bracket_opcode_end);
    if (end_address <= begin_address || end_address - begin_address > 64u) {
        return false;
    }
    const unsigned char* const code =
        reinterpret_cast<const unsigned char*>(begin_address);
    const size_t code_size = static_cast<size_t>(
        end_address - begin_address);
    static const unsigned char mfence[] = { 0x0f, 0xae, 0xf0 };
    static const unsigned char lfence[] = { 0x0f, 0xae, 0xe8 };
    static const unsigned char rdpru[] = { 0x0f, 0x01, 0xfd };
    static const unsigned char cpuid[] = { 0x0f, 0xa2 };
    size_t position = 0u;
    if (!ContainsOpcodeSequence(
            code, code_size, mfence, sizeof(mfence), position) ||
        !ContainsOpcodeSequence(
            code, code_size, lfence, sizeof(lfence), position) ||
        !ContainsOpcodeSequence(
            code, code_size, rdpru, sizeof(rdpru), position) ||
        !ContainsOpcodeSequence(
            code, code_size, lfence, sizeof(lfence), position)) return false;
    size_t mfence_count = 0u;
    size_t lfence_count = 0u;
    size_t rdpru_count = 0u;
    for (size_t i = 0u; i < code_size; ++i)
    {
        if (i + sizeof(mfence) <= code_size &&
            std::memcmp(code + i, mfence, sizeof(mfence)) == 0) ++mfence_count;
        if (i + sizeof(lfence) <= code_size &&
            std::memcmp(code + i, lfence, sizeof(lfence)) == 0) ++lfence_count;
        if (i + sizeof(rdpru) <= code_size &&
            std::memcmp(code + i, rdpru, sizeof(rdpru)) == 0) ++rdpru_count;
        if (i + sizeof(cpuid) <= code_size &&
            std::memcmp(code + i, cpuid, sizeof(cpuid)) == 0) return false;
    }
    return mfence_count == 1u && lfence_count == 2u && rdpru_count == 1u;
#else
    uint64_t value = 0u;
    return !ReadRdpru(nullptr, 0u, &value);
#endif
}

struct FakeReaderState
{
    std::vector<uint64_t> Values;
    size_t Next = 0u;
    size_t FailAt = std::numeric_limits<size_t>::max();
    size_t ThrowAt = std::numeric_limits<size_t>::max();
};

bool FakeReader(void* context, uint32_t selector, uint64_t* value_out)
{
    FakeReaderState* state = static_cast<FakeReaderState*>(context);
    if (!state || !value_out || selector > 1u) return false;
    const size_t index = state->Next++;
    if (index == state->ThrowAt) throw 7;
    if (index == state->FailAt || index >= state->Values.size()) return false;
    *value_out = state->Values[index];
    return true;
}

#if defined(__linux__)
bool SignalReader(void*, uint32_t, uint64_t*)
{
#if defined(__linux__)
    (void)raise(SIGILL);
#endif
    return false;
}
#endif

bool TestGuardedReaderFailureRetry()
{
    FakeReaderState state;
    state.Values.push_back(10u);
    state.Values.push_back(20u);
    uint64_t value = 0u;
    int signal_number = 0;
    state.ThrowAt = 0u;
    if (GuardedCounterRead(
            FakeReader, &state, 0u, value, signal_number) ||
        !AccessProbeInactive()) return false;
    state.Next = 0u;
    state.ThrowAt = std::numeric_limits<size_t>::max();
    if (!GuardedCounterRead(
            FakeReader, &state, 0u, value, signal_number) ||
        value != 10u || signal_number != 0 || !AccessProbeInactive())
    {
        return false;
    }
#if defined(__linux__)
    struct sigaction before;
    struct sigaction after;
    if (sigaction(SIGILL, nullptr, &before) != 0) return false;
    if (GuardedCounterRead(
            SignalReader, nullptr, 0u, value, signal_number) ||
        signal_number != SIGILL || !AccessProbeInactive() ||
        sigaction(SIGILL, nullptr, &after) != 0 ||
        before.sa_handler != after.sa_handler)
    {
        return false;
    }
#endif
    state.Next = 0u;
    return GuardedCounterRead(
        FakeReader, &state, 1u, value, signal_number) && value == 10u;
}

Batch128CounterArmResult SyntheticBatch(
    Scope scope,
    uint32_t iid_overhead,
    uint64_t start,
    uint64_t end,
    bool measured)
{
    Batch128CounterArmResult result;
    result.CounterStart = measured ? start : 0u;
    result.CounterEnd = measured ? end : 0u;
    result.CounterDelta = measured && end > start ? end - start : 0u;
    result.BeginReadSucceeded = measured;
    result.EndReadSucceeded = measured;
    result.VerificationMask[0] = UINT64_MAX;
    result.VerificationMask[1] = UINT64_MAX;
    for (size_t i = 0u; i < kBatchSize; ++i)
    {
        result.Results[i] = Wirehair_Success;
        result.DecodedOverheads[i] = scope == Scope::Encoder ? UINT32_MAX :
            (scope == Scope::IidReceive ? iid_overhead : 0u);
        result.DirectSystematicPackets[i] = 0u;
    }
    return result;
}

bool TestSemanticClassification()
{
    Batch128CounterArmResult result = SyntheticBatch(
        Scope::Encoder, 0u, 100u, 200u, false);
    if (!ValidSemanticBatch128Result(Scope::Encoder, 0u, result, false)) {
        return false;
    }
    result = SyntheticBatch(Scope::IidReceive, 7u, 100u, 200u, true);
    if (!ValidSemanticBatch128Result(Scope::IidReceive, 7u, result, true)) {
        return false;
    }
    Batch128CounterArmResult changed = result;
    changed.EndReadSucceeded = false;
    if (ValidSemanticBatch128Result(Scope::IidReceive, 7u, changed, true)) {
        return false;
    }
    changed = result; changed.CounterDelta = 0u;
    if (ValidSemanticBatch128Result(Scope::IidReceive, 7u, changed, true)) {
        return false;
    }
    changed = result; changed.VerificationMask[1] ^= 1u;
    if (ValidSemanticBatch128Result(Scope::IidReceive, 7u, changed, true)) {
        return false;
    }
    changed = result; changed.Results[127] = Wirehair_Error;
    if (ValidSemanticBatch128Result(Scope::IidReceive, 7u, changed, true)) {
        return false;
    }
    changed = result; ++changed.DecodedOverheads[127];
    if (ValidSemanticBatch128Result(Scope::IidReceive, 7u, changed, true)) {
        return false;
    }
    changed = result; ++changed.DirectSystematicPackets[127];
    return !ValidSemanticBatch128Result(
        Scope::IidReceive, 7u, changed, true);
}

bool TestStatisticsAndOfficialFamily()
{
    if (DoubleBits(kT95N96) != UINT64_C(0x400b71e466b42a4c) ||
        DoubleBits(kT95N576) != UINT64_C(0x400ab7ce1eaae9ba) ||
        DoubleBits(kLogLowerBound) != UINT64_C(0xbf94b004bce0abf1) ||
        DoubleBits(kLogUpperBound) != UINT64_C(0x3f944723d272a7f1))
    {
        return false;
    }
    std::ostringstream contract_evidence;
    AppendInitialEvidence(
        contract_evidence,
        kFrozenTargetCpu,
        "/tmp/wirehair-rdpru-evidence",
        "/tmp/wirehair-rdpru-launch-manifest",
        std::string(64u, '1'));
    static const char expected_margin_evidence[] =
        ",lower_log_margin=-0.020202707317519445"
        ",upper_log_margin=0.019802627296179713"
        ",critical_n96_bits=0x400b71e466b42a4c"
        ",critical_n576_bits=0x400ab7ce1eaae9ba"
        ",lower_log_margin_bits=0xbf94b004bce0abf1"
        ",upper_log_margin_bits=0x3f944723d272a7f1"
        ",official_gates=56\n";
    const std::string contract = contract_evidence.str();
    const size_t contract_marker = contract.find("STATISTICAL_CONTRACT,");
    if (contract.find(expected_margin_evidence) == std::string::npos ||
        contract_marker == std::string::npos ||
        contract.find("STATISTICAL_CONTRACT,", contract_marker + 1u) !=
            std::string::npos)
    {
        return false;
    }
    std::array<double, 96> zeros = {{ 0.0 }};
    LogStats stats;
    if (!ComputeLogStats(zeros.data(), zeros.size(), stats) ||
        stats.Mean != 0.0 || stats.Lower != 0.0 || stats.Upper != 0.0 ||
        stats.GeometricMean != 1.0 || stats.GeometricLower != 1.0 ||
        stats.GeometricUpper != 1.0 || !OfficialGate(stats))
    {
        return false;
    }
    LogStats boundary;
    boundary.Lower = kLogLowerBound;
    boundary.Upper = kLogUpperBound;
    if (!OfficialGate(boundary)) return false;
    boundary.Lower = std::nextafter(
        kLogLowerBound, -std::numeric_limits<double>::infinity());
    if (OfficialGate(boundary)) return false;
    boundary.Lower = kLogLowerBound;
    boundary.Upper = std::nextafter(
        kLogUpperBound, std::numeric_limits<double>::infinity());
    if (OfficialGate(boundary)) return false;
    boundary.Upper = kLogUpperBound;
    boundary.Lower = std::nextafter(kLogLowerBound, 0.0);
    if (!OfficialGate(boundary)) return false;
    boundary.Lower = 0.0;
    boundary.Upper = kLogUpperBound;
    if (!OfficialGate(boundary)) return false;
    boundary.Lower = kLogLowerBound;
    boundary.Upper = 0.0;
    if (!OfficialGate(boundary)) return false;

    FrozenSchedule schedule;
    if (!BuildFrozenSchedule(schedule)) return false;
    std::vector<SuperblockObservation> observations;
    observations.reserve(kSelectorSuperblockCount);
    for (size_t i = 0u; i < schedule.size(); ++i)
    {
        for (uint32_t position = 0u; position < 2u; ++position)
        {
            SuperblockObservation observation;
            observation.Plan = schedule[i];
            observation.SelectorPosition = position;
            observation.Selector = schedule[i].FirstSelector ^ position;
            observation.LogicalSums.fill(800000u);
            observation.PhysicalSums.fill(800000u);
            observation.FactorialSums.fill(400000u);
            observation.LogicalCounts.fill(8u);
            observation.PhysicalCounts.fill(8u);
            observation.FactorialCounts.fill(4u);
            if (!FinalizeObservation(observation)) return false;
            observations.push_back(observation);
        }
    }
    std::vector<size_t> by_key;
    std::array<OfficialResult, kOfficialGateCount> official;
    bool all_pass = false;
    if (!BuildObservationIndex(observations, by_key) ||
        !BuildOfficialResults(observations, by_key, official, all_pass) ||
        !all_pass) return false;
    size_t n96 = 0u;
    size_t n576 = 0u;
    std::set<std::string> identities;
    bool tuple_seen[2][2][2][7] = {};
    for (size_t i = 0u; i < official.size(); ++i)
    {
        n96 += official[i].Stats.Count == 96u ? 1u : 0u;
        n576 += official[i].Stats.Count == 576u ? 1u : 0u;
        if (!official[i].Pass || official[i].Stats.Mean != 0.0 ||
            official[i].Selector >= 2u ||
            static_cast<uint32_t>(official[i].EstimandValue) >= 2u ||
            static_cast<uint32_t>(official[i].ModeValue) >= 2u ||
            official[i].Group >= 7u ||
            official[i].Pooled != (official[i].Group == 6u) ||
            (!official[i].Pooled && official[i].Group !=
                official[i].WidthIndex * 3u +
                    static_cast<uint32_t>(official[i].ScopeValue)) ||
            tuple_seen[official[i].Selector]
                [static_cast<uint32_t>(official[i].EstimandValue)]
                [static_cast<uint32_t>(official[i].ModeValue)]
                [official[i].Group] ||
            !identities.insert(official[i].Identity).second)
        {
            return false;
        }
        tuple_seen[official[i].Selector]
            [static_cast<uint32_t>(official[i].EstimandValue)]
            [static_cast<uint32_t>(official[i].ModeValue)]
            [official[i].Group] = true;
    }
    if (identities.size() != kOfficialGateCount) return false;
    for (uint32_t selector = 0u; selector < 2u; ++selector) {
        for (uint32_t estimand = 0u; estimand < 2u; ++estimand) {
            for (uint32_t mode = 0u; mode < 2u; ++mode) {
                for (uint32_t group = 0u; group < 7u; ++group) {
                    if (!tuple_seen[selector][estimand][mode][group]) return false;
                }
            }
        }
    }

    for (uint32_t selector = 0u; selector < 2u; ++selector)
    {
        for (uint32_t estimand_index = 0u; estimand_index < 2u;
             ++estimand_index)
        {
            for (uint32_t mode_index = 0u; mode_index < 2u; ++mode_index)
            {
                for (uint32_t width = 0u; width < 2u; ++width)
                {
                    for (uint32_t scope = 0u; scope < 3u; ++scope)
                    {
                        const uint32_t cell = width * 9u + scope * 3u;
                        const size_t key = ObservationKey(
                            selector,
                            cell,
                            0u,
                            static_cast<Mode>(mode_index));
                        if (key >= by_key.size() || by_key[key] >= observations.size()) {
                            return false;
                        }
                        SuperblockObservation& changed = observations[by_key[key]];
                        double& log_value = estimand_index == 0u ?
                            changed.LogLogicalAa : changed.LogPhysicalXy;
                        const double saved = log_value;
                        log_value = 10.0;
                        if (!BuildOfficialResults(
                                observations, by_key, official, all_pass) ||
                            all_pass)
                        {
                            return false;
                        }
                        bool stratum_failed = false;
                        bool pooled_failed = false;
                        size_t failed_count = 0u;
                        for (size_t result_index = 0u;
                             result_index < official.size();
                             ++result_index)
                        {
                            const OfficialResult& result = official[result_index];
                            if (!result.Pass) ++failed_count;
                            if (result.Selector != selector ||
                                static_cast<uint32_t>(result.EstimandValue) !=
                                    estimand_index ||
                                static_cast<uint32_t>(result.ModeValue) !=
                                    mode_index) continue;
                            if (!result.Pooled && result.WidthIndex == width &&
                                static_cast<uint32_t>(result.ScopeValue) == scope)
                            {
                                stratum_failed = !result.Pass;
                            }
                            if (result.Pooled) pooled_failed = !result.Pass;
                        }
                        log_value = saved;
                        if (!stratum_failed || !pooled_failed ||
                            failed_count != 2u) return false;

                        std::vector<size_t> stratum_indexes;
                        std::vector<double> stratum_saved;
                        for (size_t observation_index = 0u;
                             observation_index < observations.size();
                             ++observation_index)
                        {
                            SuperblockObservation& candidate =
                                observations[observation_index];
                            const LogicalCellInfo candidate_info =
                                GetLogicalCellInfo(candidate.Plan.Cell);
                            if (candidate.Selector != selector ||
                                static_cast<uint32_t>(
                                    candidate.Plan.ModeValue) != mode_index ||
                                candidate_info.WidthIndex != width ||
                                static_cast<uint32_t>(
                                    candidate_info.ScopeValue) != scope)
                            {
                                continue;
                            }
                            double& candidate_log = estimand_index == 0u ?
                                candidate.LogLogicalAa : candidate.LogPhysicalXy;
                            stratum_indexes.push_back(observation_index);
                            stratum_saved.push_back(candidate_log);
                            candidate_log = (stratum_indexes.size() & 1u) != 0u ?
                                -0.06 : 0.06;
                        }
                        if (stratum_indexes.size() != 96u ||
                            !BuildOfficialResults(
                                observations, by_key, official, all_pass) ||
                            all_pass)
                        {
                            return false;
                        }
                        size_t stratum_only_failed_count = 0u;
                        bool target_stratum_only_failed = false;
                        bool corresponding_pool_passed = false;
                        for (size_t result_index = 0u;
                             result_index < official.size();
                             ++result_index)
                        {
                            const OfficialResult& result = official[result_index];
                            if (!result.Pass) ++stratum_only_failed_count;
                            if (result.Selector != selector ||
                                static_cast<uint32_t>(result.EstimandValue) !=
                                    estimand_index ||
                                static_cast<uint32_t>(result.ModeValue) !=
                                    mode_index) continue;
                            if (!result.Pooled && result.WidthIndex == width &&
                                static_cast<uint32_t>(result.ScopeValue) == scope)
                            {
                                target_stratum_only_failed = !result.Pass;
                            }
                            if (result.Pooled) {
                                corresponding_pool_passed = result.Pass;
                            }
                        }
                        for (size_t restore = 0u;
                             restore < stratum_indexes.size();
                             ++restore)
                        {
                            SuperblockObservation& candidate =
                                observations[stratum_indexes[restore]];
                            double& candidate_log = estimand_index == 0u ?
                                candidate.LogLogicalAa : candidate.LogPhysicalXy;
                            candidate_log = stratum_saved[restore];
                        }
                        if (!target_stratum_only_failed ||
                            !corresponding_pool_passed ||
                            stratum_only_failed_count != 1u) return false;
                    }
                }
            }
        }
    }

    for (size_t i = 0u; i < observations.size(); ++i)
    {
        SuperblockObservation& observation = observations[i];
        if (observation.Selector == 0u &&
            observation.Plan.ModeValue == Mode::Exposed)
        {
            observation.LogLogicalAa = 0.002 +
                ((observation.Plan.Round & 1u) == 0u ? -0.01 : 0.01);
        }
    }
    if (!BuildOfficialResults(
            observations, by_key, official, all_pass) || all_pass) return false;
    size_t target_strata_pass = 0u;
    size_t target_pooled_fail = 0u;
    size_t total_pooled_construction_fail = 0u;
    for (size_t i = 0u; i < official.size(); ++i)
    {
        const OfficialResult& result = official[i];
        total_pooled_construction_fail += !result.Pass ? 1u : 0u;
        if (result.Selector == 0u &&
            result.EstimandValue == Estimand::Logical &&
            result.ModeValue == Mode::Exposed)
        {
            if (result.Pooled) target_pooled_fail += !result.Pass ? 1u : 0u;
            else target_strata_pass += result.Pass ? 1u : 0u;
        }
    }
    return n96 == kExpectedN96GateCount &&
        n576 == kExpectedN576GateCount &&
        target_strata_pass == 6u && target_pooled_fail == 1u &&
        total_pooled_construction_fail == 1u &&
        CriticalValue(95u) == 0.0;
}

QualifiedResult SyntheticQualified(
    Scope scope,
    uint32_t iid_overhead,
    int target_cpu,
    uint64_t ordinal,
    uint64_t delta)
{
    QualifiedResult result;
    result.Before.Cpu = target_cpu;
    result.Before.AffinityCount = 1u;
    result.Before.TargetInAffinity = 1u;
    result.Before.Voluntary = static_cast<int64_t>(ordinal);
    result.Before.Involuntary = static_cast<int64_t>(ordinal / 3u);
    result.After = result.Before;
    const uint64_t start = UINT64_C(1000000) + ordinal * UINT64_C(200000);
    result.Batch = SyntheticBatch(scope, iid_overhead, start, start + delta, true);
    (void)SummarizeResult(result.Batch, result.Summary);
    return result;
}

struct InjectedRunState
{
    uint64_t QualifiedCalls = 0u;
    uint64_t FailQualifiedCall = UINT64_MAX;
    uint64_t CounterReads = 0u;
    uint64_t ThrowCounterRead = UINT64_MAX;
};

bool InjectedCaptureContext(
    int target_cpu,
    ContextReceipt& receipt,
    std::string& diagnostic)
{
    receipt = ContextReceipt();
    receipt.Cpu = target_cpu;
    receipt.AffinityCount = 1u;
    receipt.TargetInAffinity = 1u;
    receipt.Voluntary = 0;
    receipt.Involuntary = 0;
    diagnostic.clear();
    return true;
}

bool InjectedCounterRead(
    void* context,
    uint32_t selector,
    uint64_t* value_out)
{
    InjectedRunState* const state =
        static_cast<InjectedRunState*>(context);
    if (!state || !value_out || selector >= kSelectorCount) return false;
    const uint64_t call = state->CounterReads++;
    *value_out = UINT64_C(1000000) + call * UINT64_C(100) + selector;
    if (call == state->ThrowCounterRead) {
        throw std::runtime_error("injected empty counter-read exception");
    }
    return true;
}

bool InjectedQualifiedBatch128(
    const CalibrationCell& base_cell,
    const LogicalCellInfo& logical_cell,
    PhysicalBank bank_name,
    const NativeCounterBracket& bracket,
    int target_cpu,
    QualifiedResult& result_out,
    std::string& diagnostic)
{
    (void)bank_name;
    InjectedRunState* const state =
        static_cast<InjectedRunState*>(bracket.Context);
    if (!state || bracket.Read != InjectedCounterRead ||
        bracket.Selector >= kSelectorCount)
    {
        diagnostic = "invalid injected qualified runner inputs";
        result_out = QualifiedResult();
        return false;
    }
    const uint64_t call = state->QualifiedCalls++;
    QualifiedResult result = SyntheticQualified(
        logical_cell.ScopeValue,
        base_cell.IidOverhead,
        target_cpu,
        call,
        UINT64_C(100000) + bracket.Selector);
    if (call == state->FailQualifiedCall)
    {
        result.Batch.CounterEnd = 0u;
        result.Batch.CounterDelta = 0u;
        result.Batch.EndReadSucceeded = false;
        (void)SummarizeResult(result.Batch, result.Summary);
        result_out = result;
        std::ostringstream message;
        message << "injected qualified failure call=" << call;
        diagnostic = message.str();
        return false;
    }
    result_out = result;
    diagnostic.clear();
    return true;
}

bool BuildInjectedRunCells(CellArray& cells)
{
    for (uint32_t cell = 0u; cell < kBaseCellCount; ++cell)
    {
        std::unique_ptr<CalibrationCell> value(new CalibrationCell);
        value->Index = cell;
        value->WidthIndex = cell / 3u;
        value->RootIndex = cell % 3u;
        value->BlockBytes = kBlockBytes[value->WidthIndex];
        value->Root = kRoots[value->RootIndex];
        value->IidOverhead = kExpectedIidOverheads[cell];
        cells[cell] = std::move(value);
    }
    return true;
}

bool BuildSyntheticHardwareEvidence(
    const FrozenSchedule& schedule,
    int target_cpu,
    std::vector<CompactHardwareReceipt>& warm,
    std::vector<CompactHardwareReceipt>& washout,
    std::vector<RawTiming>& raw,
    std::vector<SuperblockObservation>& observations)
{
    warm.clear();
    washout.clear();
    raw.clear();
    observations.clear();
    for (uint32_t cell = 0u; cell < kCellCount; ++cell)
    {
        const LogicalCellInfo info = GetLogicalCellInfo(cell);
        for (uint32_t bank = 0u; bank < 3u; ++bank)
        {
            const uint32_t first_selector = (cell ^ bank) & 1u;
            for (uint32_t position = 0u; position < 2u; ++position)
            {
                const uint32_t selector = first_selector ^ position;
                const QualifiedResult qualified = SyntheticQualified(
                    info.ScopeValue,
                    kExpectedIidOverheads[info.BaseIndex],
                    target_cpu,
                    warm.size(),
                    100000u + selector);
                warm.push_back(MakeCompactReceipt(
                    CompactKind::Warm,
                    warm.size(),
                    selector,
                    cell,
                    info.ScopeValue,
                    static_cast<PhysicalBank>(bank),
                    qualified));
            }
        }
    }
    uint64_t global_action = 0u;
    std::array<uint64_t, 2> selector_action = {{ 0u, 0u }};
    std::array<uint64_t, 2> selector_superblock = {{ 0u, 0u }};
    for (size_t plan_index = 0u; plan_index < schedule.size(); ++plan_index)
    {
        const SuperblockPlan& plan = schedule[plan_index];
        const LogicalCellInfo info = GetLogicalCellInfo(plan.Cell);
        for (uint32_t selector_position = 0u;
             selector_position < 2u;
             ++selector_position)
        {
            const uint32_t selector = plan.FirstSelector ^ selector_position;
            const uint64_t selector_sb = selector_superblock[selector]++;
            SuperblockObservation observation;
            observation.Plan = plan;
            observation.SelectorPosition = selector_position;
            observation.Selector = selector;
            for (uint32_t sequence_position = 0u;
                 sequence_position < 4u;
                 ++sequence_position)
            {
                const uint32_t sequence =
                    (plan.Rotation + sequence_position) % 4u;
                const PhysicalBank logical_a =
                    LogicalAPhysicalForSequence(plan, sequence);
                const PhysicalBank logical_b = logical_a == PhysicalBank::X ?
                    PhysicalBank::Y : PhysicalBank::X;
                for (uint32_t step = 0u; step < 8u; ++step)
                {
                    const bool measured = plan.ModeValue == Mode::Exposed ?
                        step >= 4u : (step & 1u) != 0u;
                    const uint32_t period = plan.ModeValue == Mode::Exposed ?
                        (measured ? step - 4u : step) : step / 2u;
                    const uint64_t this_global = global_action++;
                    const uint64_t this_selector = selector_action[selector]++;
                    const QualifiedResult qualified = SyntheticQualified(
                        info.ScopeValue,
                        kExpectedIidOverheads[info.BaseIndex],
                        target_cpu,
                        this_global + kWarmBatchCount,
                        100000u);
                    if (!measured)
                    {
                        CompactHardwareReceipt receipt = MakeCompactReceipt(
                            CompactKind::Washout,
                            washout.size(),
                            selector,
                            plan.Cell,
                            info.ScopeValue,
                            PhysicalBank::W,
                            qualified);
                        receipt.GlobalPanelAction = this_global;
                        receipt.SelectorPanelAction = this_selector;
                        receipt.OriginalSuperblock = plan_index;
                        receipt.SelectorSuperblock = selector_sb;
                        washout.push_back(receipt);
                        continue;
                    }
                    const LogicalLabel label = kSequences[sequence][period];
                    const PhysicalBank physical = label == LogicalLabel::A ?
                        logical_a : logical_b;
                    RawTiming receipt;
                    receipt.RawIndex = raw.size();
                    receipt.GlobalPanelAction = this_global;
                    receipt.SelectorPanelAction = this_selector;
                    receipt.OriginalSuperblock = plan_index;
                    receipt.SelectorSuperblock = selector_sb;
                    receipt.Cell = plan.Cell;
                    receipt.Round = plan.Round;
                    receipt.Visit = plan.Visit;
                    receipt.U = plan.U;
                    receipt.DirectionValue = plan.DirectionValue;
                    receipt.ModeValue = plan.ModeValue;
                    receipt.ModePosition = plan.ModePosition;
                    receipt.OrientationFlip = plan.OrientationFlip;
                    receipt.Rotation = plan.Rotation;
                    receipt.FirstSelector = plan.FirstSelector;
                    receipt.SelectorPosition = selector_position;
                    receipt.Selector = selector;
                    receipt.SequencePosition = sequence_position;
                    receipt.Sequence = sequence;
                    receipt.SequenceFamily = sequence < 2u ? 0u : 1u;
                    receipt.Period = period;
                    receipt.Label = label;
                    receipt.LogicalAPhysical = logical_a;
                    receipt.Physical = physical;
                    receipt.Predecessor = PredecessorClass::Washout;
                    if (plan.ModeValue == Mode::Exposed && period > 0u)
                    {
                        const uint32_t previous = static_cast<uint32_t>(
                            kSequences[sequence][period - 1u]);
                        receipt.Predecessor = static_cast<PredecessorClass>(
                            1u + previous * 2u +
                                static_cast<uint32_t>(label));
                    }
                    receipt.Result = qualified;
                    if (!AccumulateMeasured(
                            observation,
                            label,
                            physical,
                            qualified.Batch.CounterDelta)) return false;
                    raw.push_back(receipt);
                }
            }
            if (!FinalizeObservation(observation)) return false;
            observations.push_back(observation);
        }
    }
    return warm.size() == kWarmBatchCount &&
        washout.size() == kWashoutBatchCount &&
        raw.size() == kMeasuredBatchCount &&
        observations.size() == kSelectorSuperblockCount;
}

bool TestReceiptPipelinesAndMutations()
{
    FrozenSchedule schedule;
    if (!BuildFrozenSchedule(schedule)) return false;
    const int target_cpu = 7;
    std::vector<CompactHardwareReceipt> warm;
    std::vector<CompactHardwareReceipt> washout;
    std::vector<RawTiming> raw;
    std::vector<SuperblockObservation> observations;
    std::string diagnostic;
    if (!BuildSyntheticHardwareEvidence(
            schedule, target_cpu, warm, washout, raw, observations) ||
        !ValidateRawReceipts(
            schedule, raw, observations, target_cpu, diagnostic) ||
        !ValidateCompactReceipts(
            schedule, warm, washout, raw, target_cpu, diagnostic))
    {
        return false;
    }
    const std::string raw_hash = RawSha256(raw);
    const std::string compact_hash = CompactSha256(warm, washout);
    const std::string mutated_protocol =
        std::string(kBracketProtocol) + "-mutation";
    if (!IsLowerHexSha256(raw_hash.c_str()) ||
        !IsLowerHexSha256(compact_hash.c_str()) ||
        RawSha256WithBracketProtocol(raw, mutated_protocol) == raw_hash ||
        CompactSha256WithBracketProtocol(
            warm, washout, mutated_protocol) == compact_hash)
    {
        return false;
    }

    ++raw[0].GlobalPanelAction;
    if (RawSha256(raw) == raw_hash) return false;
    --raw[0].GlobalPanelAction;
    raw[0].Result.Batch.Results[127] = Wirehair_Error;
    if (RawSha256(raw) == raw_hash) return false;
    raw[0].Result.Batch.Results[127] = Wirehair_Success;
    ++raw[0].Result.Before.Voluntary;
    if (RawSha256(raw) == raw_hash) return false;
    --raw[0].Result.Before.Voluntary;
    const std::string original_inner = raw[0].Result.Summary.InnerSha256;
    raw[0].Result.Summary.InnerSha256[0] =
        raw[0].Result.Summary.InnerSha256[0] == '0' ? '1' : '0';
    if (RawSha256(raw) == raw_hash) return false;
    raw[0].Result.Summary.InnerSha256 = original_inner;

    ++warm[0].CounterDelta;
    if (CompactSha256(warm, washout) == compact_hash) return false;
    --warm[0].CounterDelta;
    std::swap(warm[0], warm[1]);
    if (ValidateCompactReceipts(
            schedule, warm, washout, raw, target_cpu, diagnostic)) return false;
    std::swap(warm[0], warm[1]);
    const uint32_t selector = washout[0].Selector;
    washout[0].Selector = 2u;
    if (ValidateCompactReceipts(
            schedule, warm, washout, raw, target_cpu, diagnostic)) return false;
    washout[0].Selector = UINT32_MAX;
    if (ValidateCompactReceipts(
            schedule, warm, washout, raw, target_cpu, diagnostic)) return false;
    washout[0].Selector = selector;
    ++washout[0].GlobalPanelAction;
    if (ValidateCompactReceipts(
            schedule, warm, washout, raw, target_cpu, diagnostic)) return false;
    --washout[0].GlobalPanelAction;
    ++raw[0].Result.Before.Cpu;
    const bool rejected_target = !ValidateRawReceipts(
        schedule, raw, observations, target_cpu, diagnostic);
    --raw[0].Result.Before.Cpu;
    return rejected_target;
}

bool TestFixtureFreeze()
{
    const std::string fixture = DeterministicFixtureSha256(UINT32_MAX);
    const std::string mutation = DeterministicFixtureSha256(5u);
    std::cout << "SELFTEST_FIXTURE,sha256=" << fixture << '\n';
    return fixture == kExpectedFixtureSha256 &&
        IsLowerHexSha256(fixture.c_str()) &&
        IsLowerHexSha256(mutation.c_str()) && mutation != fixture;
}

bool TestEmptyReceiptAndGate()
{
    const int target_cpu = 11;
    std::vector<EmptyReceipt> empty;
    empty.reserve(static_cast<size_t>(kEmptySampleCount) * 2u);
    for (uint32_t sample = 0u; sample < kEmptySampleCount; ++sample)
    {
        const uint32_t first_selector = sample & 1u;
        for (uint32_t position = 0u; position < 2u; ++position)
        {
            EmptyReceipt receipt;
            receipt.Sample = sample;
            receipt.SelectorPosition = position;
            receipt.Selector = first_selector ^ position;
            receipt.Start = UINT64_C(1000000) + sample * 1000u + position * 100u;
            receipt.Delta = 100u;
            receipt.End = receipt.Start + receipt.Delta;
            receipt.BeginReadSucceeded = 1u;
            receipt.EndReadSucceeded = 1u;
            receipt.Before.Cpu = target_cpu;
            receipt.Before.AffinityCount = 1u;
            receipt.Before.TargetInAffinity = 1u;
            receipt.Before.Voluntary = sample;
            receipt.Before.Involuntary = sample / 7u;
            receipt.After = receipt.Before;
            empty.push_back(receipt);
        }
    }
    std::string diagnostic;
    if (!ValidateEmptyReceipts(empty, target_cpu, diagnostic)) return false;
    const std::string hash = EmptySha256(empty);
    if (!IsLowerHexSha256(hash.c_str()) ||
        EmptySha256WithBracketProtocol(
            empty, std::string(kBracketProtocol) + "-mutation") == hash)
    {
        return false;
    }
    ++empty[0].Sample;
    if (ValidateEmptyReceipts(empty, target_cpu, diagnostic) ||
        EmptySha256(empty) == hash) return false;
    --empty[0].Sample;
    ++empty[0].Delta;
    if (ValidateEmptyReceipts(empty, target_cpu, diagnostic)) return false;
    --empty[0].Delta;
    ++empty[0].Before.Cpu;
    if (ValidateEmptyReceipts(empty, target_cpu, diagnostic)) return false;
    --empty[0].Before.Cpu;
    empty[0].BeginReadSucceeded = 0u;
    if (ValidateEmptyReceipts(empty, target_cpu, diagnostic) ||
        EmptySha256(empty) == hash) return false;
    empty[0].BeginReadSucceeded = 1u;
    std::vector<RawTiming> raw(kMeasuredBatchCount);
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        raw[i].Selector = static_cast<uint32_t>(i & 1u);
        raw[i].Result.Batch.CounterDelta = 100000u;
    }
    std::array<EmptyGate, 2> gates;
    if (!BuildEmptyGates(empty, raw, gates) ||
        !gates[0].Pass || !gates[1].Pass ||
        gates[0].P999 != 100u || gates[1].P999 != 100u) return false;
    for (uint32_t sample = 0u; sample < kEmptySampleCount; ++sample)
    {
        const uint64_t delta = sample < 4092u ? 100u : 1000u;
        empty[static_cast<size_t>(sample) * 2u].Delta = delta;
        empty[static_cast<size_t>(sample) * 2u + 1u].Delta = delta;
    }
    if (!BuildEmptyGates(empty, raw, gates) ||
        gates[0].P999 != 100u || gates[1].P999 != 100u ||
        !gates[0].Pass || !gates[1].Pass) return false;
    for (uint32_t sample = 0u; sample < kEmptySampleCount; ++sample)
    {
        const uint64_t delta = sample < 4091u ? 100u : 1000u;
        empty[static_cast<size_t>(sample) * 2u].Delta = delta;
        empty[static_cast<size_t>(sample) * 2u + 1u].Delta = delta;
    }
    if (!BuildEmptyGates(empty, raw, gates) ||
        gates[0].P999 != 1000u || gates[1].P999 != 1000u ||
        gates[0].Pass || gates[1].Pass) return false;
    for (size_t i = 0u; i < empty.size(); ++i) empty[i].Delta = 501u;
    if (!BuildEmptyGates(empty, raw, gates) ||
        gates[0].Pass || gates[1].Pass) return false;
    return kEmptyP999OneBasedIndex == 4092u &&
        kEmptyOverheadMultiplier == 200u;
}

size_t CountTextOccurrences(
    const std::string& text,
    const std::string& needle);
bool HasSingleTerminal(
    const std::string& payload,
    const std::string& value);

bool TestRunStateFailureFinalization()
{
    const CalibrationRunFunctions& production = ProductionRunFunctions();
    if (production.CounterContext != nullptr ||
        production.CounterRead != ReadRdpru ||
        production.Capture != CaptureContext ||
        production.RunQualified != ExecuteQualifiedBatch128)
    {
        return false;
    }

    FrozenSchedule schedule;
    CellArray cells;
    if (!BuildFrozenSchedule(schedule) || !BuildInjectedRunCells(cells) ||
        schedule[0].Cell != 0u || schedule[0].ModeValue != Mode::Isolated ||
        schedule[0].FirstSelector != 0u || schedule[0].Rotation != 3u ||
        schedule[0].OrientationFlip != 1u)
    {
        return false;
    }

    enum class FailureKind : uint32_t
    {
        Warm,
        EmptyFirstRead,
        EmptySecondRead,
        Washout,
        Measured
    };
    struct FailureCase
    {
        FailureKind Kind;
        const char* Stage;
        const char* Identity;
        uint64_t QualifiedFailure;
        uint64_t ReaderThrow;
        size_t EmptyCount;
        size_t WarmCount;
        size_t WashoutCount;
        size_t RawCount;
        uint64_t FailureOrdinal;
        uint32_t Selector;
        uint32_t Cell;
        PhysicalBank Bank;
        uint32_t BeginOk;
        uint32_t EndOk;
    };
    const FailureCase cases[] = {
        {
            FailureKind::Warm,
            "warm",
            "warm;cell=0;bank=1;selector=0",
            3u,
            UINT64_MAX,
            0u, 3u, 0u, 0u,
            3u, 0u, 0u, PhysicalBank::Y, 1u, 0u
        },
        {
            FailureKind::EmptyFirstRead,
            "empty",
            "empty;sample=1;selector_position=0;selector=1",
            UINT64_MAX,
            4u,
            2u, kWarmBatchCount, 0u, 0u,
            2u, 1u, UINT32_MAX, PhysicalBank::W, 0u, 0u
        },
        {
            FailureKind::EmptySecondRead,
            "empty",
            "empty;sample=1;selector_position=0;selector=1",
            UINT64_MAX,
            5u,
            2u, kWarmBatchCount, 0u, 0u,
            2u, 1u, UINT32_MAX, PhysicalBank::W, 1u, 0u
        },
        {
            FailureKind::Washout,
            "washout",
            "washout;global_action=2;selector_action=2;"
            "original_superblock=0;selector_superblock=0;cell=0;selector=0",
            2u,
            UINT64_MAX,
            static_cast<size_t>(kEmptySampleCount) * kSelectorCount,
            kWarmBatchCount, 1u, 1u,
            1u, 0u, 0u, PhysicalBank::W, 1u, 0u
        },
        {
            FailureKind::Measured,
            "measured",
            "measured;global_action=3;selector_action=3;"
            "original_superblock=0;selector_superblock=0;cell=0;selector=0;"
            "sequence_position=0;period=1",
            3u,
            UINT64_MAX,
            static_cast<size_t>(kEmptySampleCount) * kSelectorCount,
            kWarmBatchCount, 2u, 1u,
            1u, 0u, 0u, PhysicalBank::Y, 1u, 0u
        }
    };
    const int target_cpu = 17;
    for (size_t case_index = 0u;
         case_index < sizeof(cases) / sizeof(cases[0]);
         ++case_index)
    {
        const FailureCase& expected = cases[case_index];
        InjectedRunState injection;
        CalibrationRunFunctions functions = {
            &injection,
            InjectedCounterRead,
            InjectedCaptureContext,
            InjectedQualifiedBatch128
        };
        std::ostringstream evidence;
        std::string diagnostic;
        {
            RunEvidenceState state(evidence, diagnostic);
            if (expected.Kind == FailureKind::Warm)
            {
                injection.FailQualifiedCall = expected.QualifiedFailure;
                if (RunPinnedWarmups(
                        cells,
                        target_cpu,
                        state.Warm,
                        state.Failure,
                        diagnostic,
                        functions))
                {
                    return false;
                }
            }
            else
            {
                if (!RunPinnedWarmups(
                        cells,
                        target_cpu,
                        state.Warm,
                        state.Failure,
                        diagnostic,
                        functions) || state.Failure.Present)
                {
                    return false;
                }
                if (expected.Kind == FailureKind::EmptyFirstRead ||
                    expected.Kind == FailureKind::EmptySecondRead)
                {
                    injection.ThrowCounterRead = expected.ReaderThrow;
                    if (RunEmptySamples(
                            target_cpu,
                            state.Empty,
                            state.Failure,
                            diagnostic,
                            functions))
                    {
                        return false;
                    }
                }
                else
                {
                    if (!RunEmptySamples(
                            target_cpu,
                            state.Empty,
                            state.Failure,
                            diagnostic,
                            functions) || state.Failure.Present ||
                        !ValidateEmptyReceipts(
                            state.Empty, target_cpu, diagnostic))
                    {
                        return false;
                    }
                    injection.QualifiedCalls = 0u;
                    injection.FailQualifiedCall = expected.QualifiedFailure;
                    std::vector<SuperblockObservation> observations;
                    if (RunFrozenPanel(
                            cells,
                            schedule,
                            target_cpu,
                            state.Raw,
                            state.Washout,
                            observations,
                            state.Failure,
                            diagnostic,
                            functions))
                    {
                        return false;
                    }
                }
            }

            if (!state.Failure.Present ||
                state.Failure.Stage != expected.Stage ||
                state.Failure.Identity != expected.Identity ||
                state.Failure.Ordinal != expected.FailureOrdinal ||
                state.Failure.Selector != expected.Selector ||
                state.Failure.Cell != expected.Cell ||
                state.Empty.size() != expected.EmptyCount ||
                state.Warm.size() != expected.WarmCount ||
                state.Washout.size() != expected.WashoutCount ||
                state.Raw.size() != expected.RawCount)
            {
                return false;
            }
            if (state.Failure.EmptyAttempt)
            {
                const EmptyReceipt& failed = state.Failure.Empty;
                if (failed.Sample != 1u ||
                    failed.SelectorPosition != 0u ||
                    failed.BeginReadSucceeded != expected.BeginOk ||
                    failed.EndReadSucceeded != expected.EndOk ||
                    failed.Before.Cpu != target_cpu ||
                    failed.After.Cpu != target_cpu ||
                    failed.Before.AffinityCount != 1u ||
                    failed.After.AffinityCount != 1u ||
                    failed.Before.TargetInAffinity != 1u ||
                    failed.After.TargetInAffinity != 1u ||
                    diagnostic.find("empty-bracket reader threw") ==
                        std::string::npos)
                {
                    return false;
                }
                if (expected.Kind == FailureKind::EmptyFirstRead &&
                    (failed.Start == 0u || failed.End != 0u)) return false;
                if (expected.Kind == FailureKind::EmptySecondRead &&
                    (failed.Start == 0u || failed.End <= failed.Start)) {
                    return false;
                }
            }
            else
            {
                const QualifiedResult& failed = state.Failure.Qualified;
                if (state.Failure.Bank != expected.Bank ||
                    failed.Before.Cpu != target_cpu ||
                    failed.After.Cpu != target_cpu ||
                    !failed.Batch.BeginReadSucceeded ||
                    failed.Batch.EndReadSucceeded ||
                    failed.Batch.CounterStart == 0u ||
                    failed.Batch.CounterEnd != 0u ||
                    failed.Batch.CounterDelta != 0u ||
                    diagnostic.find("injected qualified failure") ==
                        std::string::npos)
                {
                    return false;
                }
            }
            if (!state.Warm.empty() &&
                (state.Warm.back().Ordinal + 1u != state.Warm.size() ||
                 !ValidCompactSemantics(state.Warm.back()))) return false;
            if (!state.Washout.empty() &&
                (state.Washout.back().Ordinal + 1u != state.Washout.size() ||
                 !ValidCompactSemantics(state.Washout.back()))) return false;
            if (!state.Raw.empty() &&
                state.Raw.back().RawIndex + 1u != state.Raw.size()) return false;
            state.Flush();
        }
        const std::string payload = FinalizeCalibrationPayload(
            evidence.str(), false, diagnostic);
        std::ostringstream expected_counts;
        expected_counts << "RUN_STATE_COUNTS,empty=" << expected.EmptyCount
                        << ",warm=" << expected.WarmCount
                        << ",washout=" << expected.WashoutCount
                        << ",raw=" << expected.RawCount
                        << ",failed_attempt=1\n";
        if (!HasSingleTerminal(payload, "INVALID") ||
            CountTextOccurrences(payload, "RUN_STATE_FLUSH_COMPLETE=1\n") != 1u ||
            CountTextOccurrences(payload, "FAILED_ATTEMPT,") != 1u ||
            payload.find(expected_counts.str()) == std::string::npos ||
            payload.find(std::string("stage_hex=") +
                HexEncode(expected.Stage)) == std::string::npos ||
            payload.find(std::string("identity_hex=") +
                HexEncode(expected.Identity)) == std::string::npos ||
            CountTextOccurrences(payload, "EMPTY,") != expected.EmptyCount ||
            CountTextOccurrences(payload, "WARM,") != expected.WarmCount ||
            CountTextOccurrences(payload, "WASHOUT,") !=
                expected.WashoutCount ||
            CountTextOccurrences(payload, "RAW,") != expected.RawCount)
        {
            return false;
        }
        std::ostringstream flags;
        flags << "begin_ok=" << expected.BeginOk
              << ",end_ok=" << expected.EndOk;
        if (payload.find(flags.str()) == std::string::npos) return false;
    }
    return true;
}

size_t CountTextOccurrences(
    const std::string& text,
    const std::string& needle)
{
    if (needle.empty()) return 0u;
    size_t count = 0u;
    size_t offset = 0u;
    while ((offset = text.find(needle, offset)) != std::string::npos)
    {
        ++count;
        offset += needle.size();
    }
    return count;
}

bool HasSingleTerminal(
    const std::string& payload,
    const std::string& value)
{
    static const std::string marker =
        "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION=";
    const std::string terminal = marker + value + "\n";
    return CountTextOccurrences(payload, marker) == 1u &&
        payload.size() >= terminal.size() &&
        payload.compare(
            payload.size() - terminal.size(), terminal.size(), terminal) == 0;
}

bool TestAtomicEvidencePublication()
{
#if defined(__linux__) && defined(SYS_renameat2)
    char directory_template[] = "/tmp/wh2-rdpru-selftest-XXXXXX";
    char* const directory = mkdtemp(directory_template);
    if (!directory) return false;
    const std::string output = std::string(directory) + "/evidence.csv";
    const std::string accumulated = "alpha\nline-two\n";
    std::string control_diagnostic;
    control_diagnostic.push_back('\0');
    control_diagnostic.push_back('\n');
    control_diagnostic.push_back('\r');
    control_diagnostic.push_back(',');
    control_diagnostic.push_back('\1');
    const std::string payload = FinalizeCalibrationPayload(
        accumulated, true, control_diagnostic);
    const std::string invalid_payload = FinalizeCalibrationPayload(
        accumulated, false, control_diagnostic);
    const std::string authority_receipt = PublicationReceiptLine(
        "PASS", output, payload);
    std::ostringstream expected_pairing;
    expected_pairing << "payload_length=" << payload.size()
                     << ",payload_sha256="
                     << wirehair::wh2_benchmark::Sha256Hex(payload);
    std::ostringstream authority_output;
    std::ostringstream failed_authority_output;
    failed_authority_output.setstate(std::ios::badbit);
    std::string diagnostic;
    bool ok = HasSingleTerminal(payload, "PENDING_PUBLICATION") &&
        HasSingleTerminal(invalid_payload, "INVALID") &&
        payload.find("=VALID\n") == std::string::npos &&
        invalid_payload.find("=VALID\n") == std::string::npos &&
        payload.find("CALIBRATION_DIAGNOSTIC_HEX=000a0d2c01\n") !=
            std::string::npos &&
        payload.find(control_diagnostic) == std::string::npos &&
        authority_receipt.find(expected_pairing.str()) != std::string::npos &&
        EmitPublicationAuthority(
            authority_output, authority_receipt, true, true) &&
        authority_output.str().find(authority_receipt) == 0u &&
        HasSingleTerminal(authority_output.str(), "VALID") &&
        !EmitPublicationAuthority(
            failed_authority_output, authority_receipt, true, true) &&
        payload.find("EVIDENCE_PUBLICATION=PASS") == std::string::npos &&
        CheckOutputDestinationNew(output, diagnostic) &&
        AtomicPublishNew(output, payload, diagnostic);
    std::ifstream input(output.c_str(), std::ios::binary);
    std::ostringstream captured;
    captured << input.rdbuf();
    ok = ok && input.good() && captured.str() == payload &&
        !CheckOutputDestinationNew(output, diagnostic) &&
        !AtomicPublishNew(output, "replacement", diagnostic);
    const AtomicFaultPoint faults[] = {
        AtomicFaultPoint::Write,
        AtomicFaultPoint::FileFsync,
        AtomicFaultPoint::Rename,
        AtomicFaultPoint::DirectoryFsync
    };
    for (size_t i = 0u; i < sizeof(faults) / sizeof(faults[0]); ++i)
    {
        std::ostringstream path;
        path << directory << "/fault-" << i << ".csv";
        if (AtomicPublishNewImpl(path.str(), payload, faults[i], diagnostic)) {
            ok = false;
        }
        struct stat status;
        errno = 0;
        if (lstat(path.str().c_str(), &status) == 0 || errno != ENOENT) {
            ok = false;
        }
    }
    const AtomicFaultPoint retry_faults[] = {
        AtomicFaultPoint::FileFsyncEintrOnce,
        AtomicFaultPoint::DirectoryFsyncEintrOnce
    };
    std::vector<std::string> retry_paths;
    for (size_t i = 0u;
         i < sizeof(retry_faults) / sizeof(retry_faults[0]);
         ++i)
    {
        std::ostringstream path;
        path << directory << "/retry-" << i << ".csv";
        retry_paths.push_back(path.str());
        if (!AtomicPublishNewImpl(
                path.str(), payload, retry_faults[i], diagnostic)) ok = false;
        std::ifstream retry_input(path.str().c_str(), std::ios::binary);
        std::ostringstream retry_text;
        retry_text << retry_input.rdbuf();
        if (retry_text.str() != payload) ok = false;
    }
    const std::string cleanup_unlink_path =
        std::string(directory) + "/cleanup-unlink.csv";
    if (AtomicPublishNewImpl(
            cleanup_unlink_path,
            payload,
            AtomicFaultPoint::DirectoryFsyncCleanupUnlink,
            diagnostic)) ok = false;
    struct stat cleanup_unlink_status;
    std::ifstream cleanup_unlink_input(
        cleanup_unlink_path.c_str(), std::ios::binary);
    std::ostringstream cleanup_unlink_text;
    cleanup_unlink_text << cleanup_unlink_input.rdbuf();
    if (lstat(cleanup_unlink_path.c_str(), &cleanup_unlink_status) != 0 ||
        diagnostic.find("cleanup could not be proved") == std::string::npos ||
        cleanup_unlink_text.str() != payload ||
        !HasSingleTerminal(
            cleanup_unlink_text.str(), "PENDING_PUBLICATION") ||
        cleanup_unlink_text.str().find("=VALID\n") != std::string::npos ||
        cleanup_unlink_text.str().find("EVIDENCE_PUBLICATION=PASS") !=
            std::string::npos ||
        unlink(cleanup_unlink_path.c_str()) != 0)
    {
        ok = false;
    }
    const std::string cleanup_retry_path =
        std::string(directory) + "/cleanup-retry.csv";
    if (AtomicPublishNewImpl(
            cleanup_retry_path,
            payload,
            AtomicFaultPoint::DirectoryFsyncCleanupFsyncEintrOnce,
            diagnostic)) ok = false;
    struct stat cleanup_retry_status;
    errno = 0;
    if (lstat(cleanup_retry_path.c_str(), &cleanup_retry_status) == 0 ||
        errno != ENOENT ||
        diagnostic.find("cleanup could not be proved") != std::string::npos)
    {
        ok = false;
    }
    const std::string cleanup_fsync_path =
        std::string(directory) + "/cleanup-fsync.csv";
    if (AtomicPublishNewImpl(
            cleanup_fsync_path,
            payload,
            AtomicFaultPoint::DirectoryFsyncCleanupFsync,
            diagnostic)) ok = false;
    struct stat cleanup_fsync_status;
    errno = 0;
    if (lstat(cleanup_fsync_path.c_str(), &cleanup_fsync_status) == 0 ||
        errno != ENOENT ||
        diagnostic.find("cleanup could not be proved") == std::string::npos)
    {
        ok = false;
    }
    const std::string missing = std::string(directory) + "/missing/out.csv";
    ok = ok && !AtomicPublishNew(missing, payload, diagnostic);
    const std::string standard_exception =
        std::string(directory) + "/standard-exception.csv";
    const std::string nonstandard_exception =
        std::string(directory) + "/nonstandard-exception.csv";
    const std::string publication_exception =
        std::string(directory) + "/publication-exception.csv";
    const std::string injected_manifest =
        std::string(directory) + "/not-read-by-injection";
    const std::string injected_manifest_sha(64u, '1');
    ok = ok && RunCalibrationBoundary(
            kFrozenTargetCpu,
            standard_exception,
            injected_manifest,
            injected_manifest_sha,
            1u) == 1 &&
        RunCalibrationBoundary(
            kFrozenTargetCpu,
            nonstandard_exception,
            injected_manifest,
            injected_manifest_sha,
            2u) == 1 &&
        RunCalibrationBoundary(
            kFrozenTargetCpu,
            publication_exception,
            injected_manifest,
            injected_manifest_sha,
            3u) == 1;
    std::ifstream standard_input(standard_exception.c_str(), std::ios::binary);
    std::ifstream nonstandard_input(
        nonstandard_exception.c_str(), std::ios::binary);
    std::ifstream publication_input(
        publication_exception.c_str(), std::ios::binary);
    std::ostringstream standard_text;
    std::ostringstream nonstandard_text;
    std::ostringstream publication_text;
    standard_text << standard_input.rdbuf();
    nonstandard_text << nonstandard_input.rdbuf();
    publication_text << publication_input.rdbuf();
    const std::string exception_payloads[] = {
        standard_text.str(), nonstandard_text.str(), publication_text.str()
    };
    for (size_t i = 0u;
         i < sizeof(exception_payloads) / sizeof(exception_payloads[0]);
         ++i)
    {
        ok = ok && HasSingleTerminal(exception_payloads[i], "INVALID") &&
            exception_payloads[i].find("=VALID\n") == std::string::npos &&
            exception_payloads[i].find("=PENDING_PUBLICATION\n") ==
                std::string::npos;
    }
    ok = ok && nonstandard_text.str().find("PARTIAL_RECEIPT") !=
            std::string::npos &&
        publication_text.str().find("injected-during-publication") !=
            std::string::npos;
    const std::string standard_diagnostic =
        "unexpected standard exception: "
        "injected standard exception before accumulation";
    const std::string nonstandard_diagnostic =
        "unexpected non-standard exception";
    const std::string publication_diagnostic =
        "unexpected standard exception: injected publication-time exception";
    ok = ok && standard_text.str().find(
            "CALIBRATION_DIAGNOSTIC_LENGTH=" +
            DecimalU32(static_cast<uint32_t>(standard_diagnostic.size())) +
            "\nCALIBRATION_DIAGNOSTIC_HEX=" +
            HexEncode(standard_diagnostic) + "\n") != std::string::npos &&
        nonstandard_text.str().find(
            "CALIBRATION_DIAGNOSTIC_LENGTH=" +
            DecimalU32(static_cast<uint32_t>(nonstandard_diagnostic.size())) +
            "\nCALIBRATION_DIAGNOSTIC_HEX=" +
            HexEncode(nonstandard_diagnostic) + "\n") != std::string::npos &&
        publication_text.str().find(
            "CALIBRATION_DIAGNOSTIC_LENGTH=" +
            DecimalU32(static_cast<uint32_t>(publication_diagnostic.size())) +
            "\nCALIBRATION_DIAGNOSTIC_HEX=" +
            HexEncode(publication_diagnostic) + "\n") != std::string::npos;
    const bool removed_file = unlink(output.c_str()) == 0;
    const bool removed_standard = unlink(standard_exception.c_str()) == 0;
    const bool removed_nonstandard = unlink(nonstandard_exception.c_str()) == 0;
    const bool removed_publication = unlink(publication_exception.c_str()) == 0;
    bool removed_retries = true;
    for (size_t i = 0u; i < retry_paths.size(); ++i) {
        removed_retries = unlink(retry_paths[i].c_str()) == 0 && removed_retries;
    }
    const bool removed_directory = rmdir(directory) == 0;
    return ok && removed_file && removed_standard && removed_nonstandard &&
        removed_publication && removed_retries && removed_directory;
#else
    std::string diagnostic;
    return !AtomicPublishNew("unsupported", "x", diagnostic);
#endif
}

enum class CliAction
{
    Invalid,
    SelfTest,
    Calibrate
};

struct CliOptions
{
    CliAction Action = CliAction::Invalid;
    int Cpu = -1;
    std::string Output;
    std::string LaunchManifest;
    std::string LaunchManifestSha256;
};

bool ParseCommandLine(int argc, char** argv, CliOptions& options)
{
    CliOptions value;
    if (argc == 2 && argv && argv[0] && argv[1] &&
        std::strcmp(argv[1], "--selftest") == 0)
    {
        value.Action = CliAction::SelfTest;
        options = value;
        return true;
    }
    if (argc == 10 && argv && argv[0] && argv[1] && argv[2] && argv[3] &&
        argv[4] && argv[5] && argv[6] && argv[7] && argv[8] && argv[9] &&
        std::strcmp(argv[1], "--calibrate") == 0 &&
        std::strcmp(argv[2], "--cpu") == 0 &&
        std::strcmp(argv[4], "--output") == 0 &&
        std::strcmp(argv[6], "--launch-manifest") == 0 &&
        std::strcmp(argv[8], "--launch-manifest-sha256") == 0 &&
        std::strcmp(argv[3], "50") == 0 &&
        ParseCpu(argv[3], value.Cpu) && value.Cpu == kFrozenTargetCpu &&
        IsCanonicalAbsolutePath(argv[5]) &&
        IsCanonicalAbsolutePath(argv[7]) && IsLowerHexSha256(argv[9]) &&
        !IsAllZeroHex(argv[9]))
    {
        value.Action = CliAction::Calibrate;
        value.Output = argv[5];
        value.LaunchManifest = argv[7];
        value.LaunchManifestSha256 = argv[9];
        options = value;
        return true;
    }
    return false;
}

bool TestCliSurface()
{
    CliOptions options;
    char program[] = "program";
    char selftest[] = "--selftest";
    char calibrate[] = "--calibrate";
    char cpu_option[] = "--cpu";
    char frozen_cpu[] = "50";
    char output_option[] = "--output";
    char output[] = "/tmp/new-evidence";
    char manifest_option[] = "--launch-manifest";
    char manifest[] = "/tmp/launch-manifest";
    char manifest_sha_option[] = "--launch-manifest-sha256";
    char manifest_sha[] =
        "123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0";
    char bad_cpu[] = "-1";
    char cpu49[] = "49";
    char cpu51[] = "51";
    char plus50[] = "+50";
    char leading_zero50[] = "050";
    char whitespace50[] = " 50";
    char extra[] = "extra";
    char* self_argv[] = { program, selftest };
    if (!ParseCommandLine(2, self_argv, options) ||
        options.Action != CliAction::SelfTest) return false;
    char* calibration_argv[] = {
        program, calibrate, cpu_option, frozen_cpu, output_option, output,
        manifest_option, manifest, manifest_sha_option, manifest_sha
    };
    if (!ParseCommandLine(10, calibration_argv, options) ||
        options.Action != CliAction::Calibrate ||
        options.Cpu != kFrozenTargetCpu ||
        options.Output != output || options.LaunchManifest != manifest ||
        options.LaunchManifestSha256 != manifest_sha) return false;
    char* bad_cpu_argv[] = {
        program, calibrate, cpu_option, bad_cpu, output_option, output,
        manifest_option, manifest, manifest_sha_option, manifest_sha
    };
    if (ParseCommandLine(10, bad_cpu_argv, options)) return false;
    char* noncanonical_cpus[] = {
        cpu49, cpu51, plus50, leading_zero50, whitespace50
    };
    for (size_t i = 0u;
         i < sizeof(noncanonical_cpus) / sizeof(noncanonical_cpus[0]);
         ++i)
    {
        calibration_argv[3] = noncanonical_cpus[i];
        if (ParseCommandLine(10, calibration_argv, options)) return false;
    }
    calibration_argv[3] = frozen_cpu;
    char* wrong_order[] = {
        program, calibrate, output_option, output, cpu_option, frozen_cpu,
        manifest_option, manifest, manifest_sha_option, manifest_sha
    };
    if (ParseCommandLine(10, wrong_order, options)) return false;
    char* too_many[] = {
        program, calibrate, cpu_option, frozen_cpu, output_option, output,
        manifest_option, manifest, manifest_sha_option, manifest_sha, extra
    };
    char newline_output[] = "/tmp/line\nforged,record";
    char* newline_argv[] = {
        program, calibrate, cpu_option, frozen_cpu, output_option, newline_output,
        manifest_option, manifest, manifest_sha_option, manifest_sha
    };
    if (!ParseCommandLine(10, newline_argv, options) ||
        options.Output != newline_output) return false;
    char relative_output[] = "relative-output";
    char dotted_output[] = "/tmp/./evidence";
    char parent_output[] = "/tmp/../evidence";
    char doubled_output[] = "/tmp//evidence";
    char relative_manifest[] = "relative-manifest";
    char* invalid_paths[] = {
        relative_output, dotted_output, parent_output, doubled_output
    };
    for (size_t i = 0u;
         i < sizeof(invalid_paths) / sizeof(invalid_paths[0]);
         ++i)
    {
        calibration_argv[5] = invalid_paths[i];
        if (ParseCommandLine(10, calibration_argv, options)) return false;
    }
    calibration_argv[5] = output;
    calibration_argv[7] = relative_manifest;
    if (ParseCommandLine(10, calibration_argv, options)) return false;
    calibration_argv[7] = manifest;
    char short_sha[65];
    char long_sha[66];
    char uppercase_sha[65];
    char nonhex_sha[65];
    char zero_sha[65];
    char control_sha[65];
    std::memcpy(short_sha, manifest_sha, sizeof(short_sha));
    short_sha[63] = '\0';
    std::memcpy(long_sha, manifest_sha, 64u);
    long_sha[64] = '0';
    long_sha[65] = '\0';
    std::memcpy(uppercase_sha, manifest_sha, sizeof(uppercase_sha));
    uppercase_sha[9] = 'A';
    std::memcpy(nonhex_sha, manifest_sha, sizeof(nonhex_sha));
    nonhex_sha[0] = 'g';
    std::memset(zero_sha, '0', 64u);
    zero_sha[64] = '\0';
    std::memcpy(control_sha, manifest_sha, sizeof(control_sha));
    control_sha[0] = '\n';
    char* invalid_shas[] = {
        short_sha, long_sha, uppercase_sha, nonhex_sha, zero_sha, control_sha
    };
    for (size_t i = 0u;
         i < sizeof(invalid_shas) / sizeof(invalid_shas[0]);
         ++i)
    {
        calibration_argv[9] = invalid_shas[i];
        if (ParseCommandLine(10, calibration_argv, options)) return false;
    }
    calibration_argv[9] = manifest_sha;
    for (size_t i = 0u; i < 10u; ++i)
    {
        char* saved = calibration_argv[i];
        calibration_argv[i] = nullptr;
        if (ParseCommandLine(10, calibration_argv, options)) return false;
        calibration_argv[i] = saved;
    }
    char empty_path[] = "";
    calibration_argv[5] = empty_path;
    if (ParseCommandLine(10, calibration_argv, options)) return false;
    calibration_argv[5] = output;
    for (int argc = 0; argc < 10; ++argc) {
        if (ParseCommandLine(argc, calibration_argv, options)) return false;
    }
    return !ParseCommandLine(11, too_many, options) &&
        !ParseCommandLine(1, self_argv, options) &&
        !ParseCommandLine(0, nullptr, options) &&
        !ParseCommandLine(2, nullptr, options);
}

#if defined(__linux__)
std::string SyntheticLaunchManifestContent(
    int target_cpu,
    const std::string& output_path)
{
    std::string executable_sha256;
    std::string executable_build_id;
    std::string diagnostic;
    if (!ReadSelfExecutableIdentity(
            executable_sha256, executable_build_id, diagnostic)) {
        return std::string();
    }
    std::ostringstream content;
    content << "manifest_schema=" << kLaunchManifestSchema << '\n'
            << "protocol_tag=" << kProtocolTag << '\n'
            << "git_commit=" << WIREHAIR_WH2_SOURCE_GIT_COMMIT << '\n'
            << "controller_sha256="
            << WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256 << '\n'
            << "codec_header_sha256="
            << WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256 << '\n'
            << "codec_source_sha256="
            << WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256 << '\n'
            << "schedule_sha256=" << kExpectedScheduleSha256 << '\n'
            << "config_sha256=" << kExpectedConfigSha256 << '\n'
            << "fixture_sha256=" << kExpectedFixtureSha256 << '\n'
            << "bracket_protocol_id=" << kBracketProtocolId << '\n'
            << "bracket_protocol=" << kBracketProtocol << '\n'
            << "target_cpu=" << target_cpu << '\n'
            << "online_cpus=" << kFrozenOnlineCpuList << '\n'
            << "logical_cpu_count=" << kFrozenLogicalCpuCount << '\n'
            << "physical_core_count=" << kFrozenPhysicalCoreCount << '\n'
            << "initial_allowed_cpus=" << kFrozenInitialAffinityList << '\n'
            << "target_cpu_thread_siblings="
            << kFrozenTargetSiblingList << '\n'
            << "cpu_selection_rule=" << kCpuSelectionRule << '\n'
            << "cpu_selection_first_cpu=" << kCpuSelectionFirstCpu << '\n'
            << "cpu_selection_last_cpu=" << kCpuSelectionLastCpu << '\n'
            << "cpu_selection_target_interrupt_count="
            << kCpuSelectionTargetInterruptCount << '\n'
            << "cpu_selection_runner_up_cpu="
            << kCpuSelectionRunnerUpCpu << '\n'
            << "cpu_selection_runner_up_interrupt_count="
            << kCpuSelectionRunnerUpInterruptCount << '\n'
            << "cpu_selection_snapshot_path="
            << kCpuSelectionSnapshotPath << '\n'
            << "cpu_selection_snapshot_sha256="
            << kCpuSelectionSnapshotSha256 << '\n'
            << "output_path_hex=" << HexEncode(output_path) << '\n'
            << "cmake_sha256=" << WIREHAIR_WH2_CMAKE_SHA256 << '\n'
            << "executable_sha256=" << executable_sha256 << '\n'
            << "elf_build_id=" << executable_build_id << '\n'
            << "harmless_probe_source_sha256="
            << kHarmlessProbeSourceSha256 << '\n'
            << "harmless_probe_binary_sha256="
            << kHarmlessProbeBinarySha256 << '\n'
            << "cpu_family=" << kFrozenCpuFamily << '\n'
            << "cpu_model=" << kFrozenCpuModel << '\n'
            << "cpu_stepping=" << kFrozenCpuStepping << '\n'
            << "leaf1_eax_hex=" << HexU32(kFrozenLeaf1Eax) << '\n'
            << "leaf1_ebx_hex=" << HexU32(kFrozenLeaf1Ebx) << '\n'
            << "leaf1_ecx_hex=" << HexU32(kFrozenLeaf1Ecx) << '\n'
            << "leaf1_edx_hex=" << HexU32(kFrozenLeaf1Edx) << '\n'
            << "leaf6_eax_hex=" << HexU32(kFrozenLeaf6Eax) << '\n'
            << "leaf6_ebx_hex=" << HexU32(kFrozenLeaf6Ebx) << '\n'
            << "leaf6_ecx_hex=" << HexU32(kFrozenLeaf6Ecx) << '\n'
            << "leaf6_edx_hex=" << HexU32(kFrozenLeaf6Edx) << '\n'
            << "leaf80000008_eax_hex="
            << HexU32(kFrozenLeaf80000008Eax) << '\n'
            << "leaf80000008_ebx_hex="
            << HexU32(kFrozenLeaf80000008Ebx) << '\n'
            << "leaf80000008_ecx_hex="
            << HexU32(kFrozenLeaf80000008Ecx) << '\n'
            << "leaf80000008_edx_hex="
            << HexU32(kFrozenLeaf80000008Edx) << '\n'
            << "leaf80000021_eax_hex="
            << HexU32(kFrozenLeaf80000021Eax) << '\n'
            << "leaf80000021_ebx_hex="
            << HexU32(kFrozenLeaf80000021Ebx) << '\n'
            << "leaf80000021_ecx_hex="
            << HexU32(kFrozenLeaf80000021Ecx) << '\n'
            << "leaf80000021_edx_hex="
            << HexU32(kFrozenLeaf80000021Edx) << '\n';
    return content.str();
}

bool WriteSelfTestFile(
    const std::string& path,
    const std::string& content)
{
#if defined(__linux__)
    const int file = open(
        path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC,
        S_IRUSR | S_IWUSR);
    if (file < 0) return false;
    std::string diagnostic;
    const bool written = WriteAll(
        file,
        content.data(),
        content.size(),
        AtomicFaultPoint::None,
        diagnostic);
    const bool synced = written && FsyncRetry(file, false);
    const bool closed = close(file) == 0;
    return written && synced && closed;
#else
    (void)path;
    (void)content;
    return false;
#endif
}
#endif

bool TestLaunchManifestGate()
{
#if defined(__linux__)
    if (!EmbeddedIdentitiesValid()) return false;
    char directory_template[] = "/tmp/wh2-rdpru-manifest-XXXXXX";
    char* const directory = mkdtemp(directory_template);
    if (!directory) return false;
    const int target_cpu = kFrozenTargetCpu;
    const std::string output = std::string(directory) + "/future-output.csv";
    const std::string manifest_path =
        std::string(directory) + "/launch-manifest.txt";
    const std::string malformed_path =
        std::string(directory) + "/malformed.txt";
    const std::string missing_path =
        std::string(directory) + "/missing.txt";
    const std::string content =
        SyntheticLaunchManifestContent(target_cpu, output);
    const std::string sha = wirehair::wh2_benchmark::Sha256Hex(content);
    const std::string malformed = "protocol_tag=duplicate\nprotocol_tag=again\n";
    const std::string malformed_sha =
        wirehair::wh2_benchmark::Sha256Hex(malformed);
    bool ok = !content.empty() && WriteSelfTestFile(manifest_path, content) &&
        WriteSelfTestFile(malformed_path, malformed);
    LaunchManifest manifest;
    std::string diagnostic;
    std::string executable_sha256;
    std::string executable_build_id;
    ok = ok && ReadSelfExecutableIdentity(
        executable_sha256, executable_build_id, diagnostic);
    ok = ok && ReadLaunchManifest(
        manifest_path,
        sha,
        target_cpu,
        output,
        executable_sha256,
        executable_build_id,
        manifest,
        diagnostic) &&
        manifest.Content == content && manifest.Sha256 == sha;
    CpuFeatures synthetic_host = FrozenCpuFeaturesForTest();
    ok = ok && ValidateLaunchManifestHost(
        manifest, synthetic_host, diagnostic);
    ++synthetic_host.Leaf80000021Eax;
    ok = ok && !ValidateLaunchManifestHost(
        manifest, synthetic_host, diagnostic);
    std::string mismatch = sha;
    mismatch[0] = mismatch[0] == '0' ? '1' : '0';
    ok = ok && !ReadLaunchManifest(
        manifest_path,
        mismatch,
        target_cpu,
        output,
        executable_sha256,
        executable_build_id,
        manifest,
        diagnostic) &&
        manifest.Content == content && manifest.Sha256 == sha;
    std::ostringstream failed_digest_evidence;
    AppendLaunchManifestEvidence(failed_digest_evidence, manifest, "FAIL");
    ok = ok && failed_digest_evidence.str().find(
            "content_length=" + DecimalU32(
                static_cast<uint32_t>(content.size()))) != std::string::npos &&
        failed_digest_evidence.str().find("sha256=" + sha) !=
            std::string::npos;
    ok = ok && !ReadLaunchManifest(
        missing_path,
        sha,
        target_cpu,
        output,
        executable_sha256,
        executable_build_id,
        manifest,
        diagnostic) && manifest.Content.empty();
    ok = ok && !ReadLaunchManifest(
        malformed_path,
        malformed_sha,
        target_cpu,
        output,
        executable_sha256,
        executable_build_id,
        manifest,
        diagnostic) && manifest.Content == malformed;
    ok = ok && !ReadLaunchManifest(
        manifest_path,
        sha,
        target_cpu + 1,
        output,
        executable_sha256,
        executable_build_id,
        manifest,
        diagnostic);
    std::map<std::string, std::string> parsed_fields;
    const std::string unknown_field_content = content + "unknown_field=x\n";
    ok = ok && !ParseLaunchManifestContent(
        unknown_field_content,
        target_cpu,
        output,
        executable_sha256,
        executable_build_id,
        parsed_fields,
        diagnostic);
    std::string blank_line_content = content;
    const size_t first_newline = blank_line_content.find('\n');
    if (first_newline == std::string::npos) ok = false;
    else blank_line_content.insert(first_newline + 1u, "\n");
    ok = ok && !ParseLaunchManifestContent(
        blank_line_content,
        target_cpu,
        output,
        executable_sha256,
        executable_build_id,
        parsed_fields,
        diagnostic);
    std::string forged_executable = content;
    const std::string executable_key = "executable_sha256=";
    const size_t executable_offset = forged_executable.find(executable_key);
    if (executable_offset == std::string::npos) ok = false;
    else
    {
        const size_t value_offset = executable_offset + executable_key.size();
        forged_executable[value_offset] =
            forged_executable[value_offset] == '0' ? '1' : '0';
    }
    ok = ok && !ParseLaunchManifestContent(
        forged_executable,
        target_cpu,
        output,
        executable_sha256,
        executable_build_id,
        parsed_fields,
        diagnostic);
    const std::set<std::string> host_keys = {
        "cpu_family", "cpu_model", "cpu_stepping",
        "leaf1_eax_hex", "leaf1_ebx_hex", "leaf1_ecx_hex", "leaf1_edx_hex",
        "leaf6_eax_hex", "leaf6_ebx_hex", "leaf6_ecx_hex", "leaf6_edx_hex",
        "leaf80000008_eax_hex", "leaf80000008_ebx_hex",
        "leaf80000008_ecx_hex", "leaf80000008_edx_hex",
        "leaf80000021_eax_hex", "leaf80000021_ebx_hex",
        "leaf80000021_ecx_hex", "leaf80000021_edx_hex"
    };
    size_t line_begin = 0u;
    size_t claim_count = 0u;
    while (line_begin < content.size())
    {
        const size_t line_end = content.find('\n', line_begin);
        if (line_end == std::string::npos) {
            ok = false;
            break;
        }
        const size_t equals = content.find('=', line_begin);
        if (equals == std::string::npos || equals >= line_end ||
            equals + 1u >= line_end)
        {
            ok = false;
            break;
        }
        const std::string key = content.substr(
            line_begin, equals - line_begin);
        std::string mutated_claim = content;
        char& first_value = mutated_claim[equals + 1u];
        first_value = first_value == '0' ? '1' : '0';
        parsed_fields.clear();
        const bool parsed_mutation = ParseLaunchManifestContent(
            mutated_claim,
            target_cpu,
            output,
            executable_sha256,
            executable_build_id,
            parsed_fields,
            diagnostic);
        if (host_keys.count(key) != 0u)
        {
            LaunchManifest host_mutation;
            host_mutation.Fields = parsed_fields;
            if (!parsed_mutation ||
                ValidateLaunchManifestHost(
                    host_mutation,
                    FrozenCpuFeaturesForTest(),
                    diagnostic)) ok = false;
        }
        else if (parsed_mutation) {
            ok = false;
        }
        std::string deleted_claim = content;
        deleted_claim.erase(line_begin, line_end + 1u - line_begin);
        parsed_fields.clear();
        if (ParseLaunchManifestContent(
                deleted_claim,
                target_cpu,
                output,
                executable_sha256,
                executable_build_id,
                parsed_fields,
                diagnostic)) ok = false;
        ++claim_count;
        line_begin = line_end + 1u;
    }
    ok = ok && claim_count == manifest.Fields.size();
    std::string nul_content = content;
    nul_content.insert(nul_content.size() / 2u, 1u, '\0');
    std::string unterminated_content = content;
    unterminated_content.resize(unterminated_content.size() - 1u);
    std::string extra_equals_content = content;
    extra_equals_content.insert(extra_equals_content.find('\n'), "=extra");
    const std::string syntax_mutations[] = {
        nul_content, unterminated_content, extra_equals_content,
        blank_line_content, malformed
    };
    for (size_t i = 0u;
         i < sizeof(syntax_mutations) / sizeof(syntax_mutations[0]);
         ++i)
    {
        parsed_fields.clear();
        if (ParseLaunchManifestContent(
                syntax_mutations[i],
                target_cpu,
                output,
                executable_sha256,
                executable_build_id,
                parsed_fields,
                diagnostic)) ok = false;
    }
    const std::string symlink_path =
        std::string(directory) + "/manifest-symlink.txt";
    const std::string oversize_path =
        std::string(directory) + "/manifest-oversize.txt";
    const std::string oversize(kMaximumLaunchManifestBytes + 1u, 'x');
    const std::string oversize_sha =
        wirehair::wh2_benchmark::Sha256Hex(oversize);
    ok = ok && symlink(manifest_path.c_str(), symlink_path.c_str()) == 0 &&
        WriteSelfTestFile(oversize_path, oversize) &&
        !ReadLaunchManifest(
            symlink_path,
            sha,
            target_cpu,
            output,
            executable_sha256,
            executable_build_id,
            manifest,
            diagnostic) &&
        !ReadLaunchManifest(
            oversize_path,
            oversize_sha,
            target_cpu,
            output,
            executable_sha256,
            executable_build_id,
            manifest,
            diagnostic);
    const bool removed_symlink = unlink(symlink_path.c_str()) == 0;
    const bool removed_oversize = unlink(oversize_path.c_str()) == 0;
    const bool removed_manifest = unlink(manifest_path.c_str()) == 0;
    const bool removed_malformed = unlink(malformed_path.c_str()) == 0;
    const bool removed_directory = rmdir(directory) == 0;
    return ok && removed_symlink && removed_oversize && removed_manifest &&
        removed_malformed && removed_directory;
#else
    LaunchManifest manifest;
    std::string diagnostic;
    return !ReadLaunchManifest(
        "/missing",
        std::string(64u, '1'),
        kFrozenTargetCpu,
        "/output",
        std::string(64u, '1'),
        "1",
        manifest,
        diagnostic);
#endif
}

bool TestMainExceptionBoundary();

bool RunSelfTest()
{
    if (!IsLowerHexSha256(WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256) ||
        !IsLowerHexSha256(WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256) ||
        !IsLowerHexSha256(WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256))
    {
        std::cerr << "native RDPRU batch128 A/A selftest failed: embedded hashes\n";
        return false;
    }
    struct Test
    {
        const char* Name;
        bool (*Run)();
    };
    const Test tests[] = {
        { "schedule-exhaustive-balance", TestScheduleBalance },
        { "sequence-transition-surface", TestSequenceSurface },
        { "canonical-ratio-signs", TestCanonicalRatioSigns },
        { "cpu-qualification-fail-closed", TestCpuQualification },
        { "host-topology-frozen", TestHostTopologyQualification },
        { "access-receipt-failures", TestAccessReceiptFailures },
        { "reader-disassembly-opcodes", TestReaderOpcodeShape },
        { "guarded-reader-failure-retry", TestGuardedReaderFailureRetry },
        { "batch128-semantic-classification", TestSemanticClassification },
        { "bonferroni-family-boundaries", TestStatisticsAndOfficialFamily },
        { "deterministic-fixture-freeze", TestFixtureFreeze },
        { "empty-receipts-nearest-rank", TestEmptyReceiptAndGate },
        { "receipt-pipelines-mutations", TestReceiptPipelinesAndMutations },
        { "partial-run-failure-evidence", TestRunStateFailureFinalization },
        { "atomic-output-no-replace", TestAtomicEvidencePublication },
        { "exact-cli-surface", TestCliSurface },
        { "launch-manifest-fail-closed", TestLaunchManifestGate }
        ,{ "outer-main-exception-boundary", TestMainExceptionBoundary }
    };
    for (size_t i = 0u; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        if (!tests[i].Run())
        {
            std::cerr << "native RDPRU batch128 A/A selftest failed: "
                      << tests[i].Name << '\n';
            return false;
        }
        std::cout << "SELFTEST,name=" << tests[i].Name << ",status=PASS\n";
    }
    std::cout << "SELFTEST_SOURCE,sha256="
              << WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256 << '\n'
              << "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_SELFTEST=PASS\n";
    return true;
}

void PrintUsage(std::ostream& output)
{
    output << "Usage: wirehair_wh2_native_rdpru_batch128_aa_calibration "
           << "--selftest\n       "
           << "wirehair_wh2_native_rdpru_batch128_aa_calibration "
           << "--calibrate --cpu 50 --output ABS_PATH "
           << "--launch-manifest ABS_PATH "
           << "--launch-manifest-sha256 HEX\n";
}

int MainBoundary(
    int argc,
    char** argv,
    uint32_t injected_exception,
    std::ostream& terminal_output,
    std::ostream& diagnostic_output)
{
    try
    {
        if (injected_exception == 1u) {
            throw std::runtime_error("injected parse exception");
        }
        if (injected_exception == 2u) throw 2;
        CliOptions options;
        if (!ParseCommandLine(argc, argv, options))
        {
            PrintUsage(diagnostic_output);
            terminal_output <<
                "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION=INVALID\n";
            terminal_output.flush();
            return 2;
        }
        if (options.Action == CliAction::SelfTest) return RunSelfTest() ? 0 : 1;
        return RunCalibrationBoundary(
            options.Cpu,
            options.Output,
            options.LaunchManifest,
            options.LaunchManifestSha256,
            0u);
    }
    catch (const std::exception& exception)
    {
        diagnostic_output << "native RDPRU main exception: "
                          << exception.what() << '\n';
    }
    catch (...)
    {
        diagnostic_output << "native RDPRU main non-standard exception\n";
    }
    terminal_output <<
        "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION=INVALID\n";
    terminal_output.flush();
    return 2;
}

bool TestMainExceptionBoundary()
{
    std::ostringstream output;
    std::ostringstream diagnostic;
    if (MainBoundary(0, nullptr, 1u, output, diagnostic) != 2 ||
        CountTextOccurrences(
            output.str(),
            "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION=INVALID\n") != 1u ||
        output.str().find("=VALID\n") != std::string::npos ||
        diagnostic.str().find("injected parse exception") == std::string::npos)
    {
        return false;
    }
    output.str(std::string());
    output.clear();
    diagnostic.str(std::string());
    diagnostic.clear();
    return MainBoundary(0, nullptr, 2u, output, diagnostic) == 2 &&
        CountTextOccurrences(
            output.str(),
            "WIREHAIR1_NATIVE_RDPRU_BATCH128_AA_CALIBRATION=INVALID\n") == 1u &&
        output.str().find("=VALID\n") == std::string::npos &&
        diagnostic.str().find("non-standard") != std::string::npos;
}

} // namespace

int main(int argc, char** argv)
{
    return MainBoundary(argc, argv, 0u, std::cout, std::cerr);
}
