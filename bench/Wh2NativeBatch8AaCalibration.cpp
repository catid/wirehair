#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <time.h>
#endif

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#define WIREHAIR_WH2_SOURCE_GIT_COMMIT "unknown"
#endif

#ifndef WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256
#define WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256 "unknown"
#endif

namespace {

using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativeReceiveFixture;
using wirehair_wh2_bench::Batch8ArmResult;

static const uint32_t kK = 2u;
static const uint32_t kLossPpm = 100000u;
static const uint32_t kReceiveOverheadCap = 256u;
static const uint32_t kRoundCount = 32u;
static const uint32_t kCellCount = 18u;
static const uint32_t kBaseCellCount = 6u;
static const uint32_t kModeCount = 2u;
static const uint32_t kSequenceCount = 4u;
static const uint32_t kPeriodsPerSequence = 4u;
static const uint32_t kBatchesPerSequence = 8u;
static const uint32_t kMeasuredBatchesPerSuperblock = 16u;
static const uint32_t kWashoutBatchesPerSuperblock = 16u;
static const uint32_t kBatchSize = 8u;
static const size_t kOfficialMarginalGateCount = 28u;
static const size_t kSuperblockCount =
    static_cast<size_t>(kCellCount) * kRoundCount * kModeCount;
static const size_t kMeasuredBatchCount =
    kSuperblockCount * kMeasuredBatchesPerSuperblock;
static const size_t kWashoutBatchCount =
    kSuperblockCount * kWashoutBatchesPerSuperblock;
static const size_t kPanelActionCount =
    kMeasuredBatchCount + kWashoutBatchCount;
static const size_t kMeasuredInnerInvocationCount =
    kMeasuredBatchCount * kBatchSize;
static const size_t kWashoutInnerInvocationCount =
    kWashoutBatchCount * kBatchSize;
static const size_t kPinnedWarmBatchCount =
    static_cast<size_t>(kCellCount) * 3u;
static const size_t kPinnedWarmInnerInvocationCount =
    kPinnedWarmBatchCount * kBatchSize;
static const size_t kUnpinnedPreflightBatchCount =
    static_cast<size_t>(kBaseCellCount) * 3u * 3u;
static const size_t kUnpinnedPreflightInnerInvocationCount =
    kUnpinnedPreflightBatchCount * kBatchSize;
static const size_t kQualifiedBatchCount =
    kPinnedWarmBatchCount + kMeasuredBatchCount + kWashoutBatchCount;
static const uint64_t kExpectedAffinityVerificationCount =
    static_cast<uint64_t>(kQualifiedBatchCount) * 2u;
static const char kExpectedScheduleSha256[] =
    "3867fdb94a3a24c8cd5deb2e5b98eb16261cb928e1661906d142a55c62df0923";
static const char kExpectedConfigSha256[] =
    "20f463a962d0ad8ea075fb5179141a58017c3bb77b6976635e6e4d3e7d12866b";
static const long double kQuantizationBudgetLimit = 0.005L;
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
    None = 0u,
    Washout = 1u,
    Aa = 2u,
    Ab = 3u,
    Ba = 4u,
    Bb = 5u
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

const char* DirectionName(Direction direction)
{
    return direction == Direction::Forward ? "forward" : "reverse";
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

const char* SequenceName(uint32_t sequence)
{
    static const char* const names[] = { "ABBA", "BAAB", "AABB", "BBAA" };
    return sequence < kSequenceCount ? names[sequence] : "invalid";
}

const char* PredecessorName(PredecessorClass predecessor)
{
    switch (predecessor)
    {
    case PredecessorClass::None: return "none";
    case PredecessorClass::Washout: return "W";
    case PredecessorClass::Aa: return "AA";
    case PredecessorClass::Ab: return "AB";
    case PredecessorClass::Ba: return "BA";
    case PredecessorClass::Bb: return "BB";
    }
    return "invalid";
}

bool AddElapsed(uint64_t value, uint64_t& total)
{
    if (value > std::numeric_limits<uint64_t>::max() - total) return false;
    total += value;
    return true;
}

bool IsLowerHexSha256(const char* text)
{
    if (!text || std::strlen(text) != 64u) return false;
    for (size_t i = 0u; i < 64u; ++i)
    {
        const char value = text[i];
        if (!((value >= '0' && value <= '9') ||
              (value >= 'a' && value <= 'f')))
        {
            return false;
        }
    }
    return true;
}

bool FlushEvidenceOutput()
{
    std::cout.flush();
    return static_cast<bool>(std::cout);
}

class FailingStreamBuffer : public std::streambuf
{
protected:
    int_type overflow(int_type character) override
    {
        return traits_type::not_eof(character);
    }

    std::streamsize xsputn(
        const char* characters,
        std::streamsize count) override
    {
        (void)characters;
        return count;
    }

    int sync() override
    {
        return -1;
    }
};

void AppendLe32(std::string& bytes, uint32_t value)
{
    for (unsigned shift = 0u; shift < 32u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendLe64(std::string& bytes, uint64_t value)
{
    for (unsigned shift = 0u; shift < 64u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
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
};

typedef std::array<SuperblockPlan, kSuperblockCount> FrozenSchedule;

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
    return plan;
}

PhysicalBank LogicalAPhysicalForSequence(
    const SuperblockPlan& plan,
    uint32_t sequence)
{
    const uint32_t family = sequence < 2u ? 0u : 1u;
    return (family ^ plan.OrientationFlip) == 0u ?
        PhysicalBank::X : PhysicalBank::Y;
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

std::string ScheduleSha256(const FrozenSchedule& schedule)
{
    std::ostringstream canonical;
    size_t action_index = 0u;
    for (size_t i = 0u; i < schedule.size(); ++i)
    {
        const SuperblockPlan& plan = schedule[i];
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
                const PhysicalBank physical = !measured ? PhysicalBank::W :
                    (label == LogicalLabel::A ? logical_a : logical_b);
                PredecessorClass predecessor = PredecessorClass::Washout;
                if (measured && plan.ModeValue == Mode::Exposed &&
                    period > 0u)
                {
                    const uint32_t previous = static_cast<uint32_t>(
                        kSequences[sequence][period - 1u]);
                    predecessor = static_cast<PredecessorClass>(
                        2u + previous * 2u + static_cast<uint32_t>(label));
                }
                canonical << "action=" << action_index++
                          << ",round=" << plan.Round
                          << ",j=" << plan.Visit
                          << ",c=" << plan.Cell
                          << ",u=" << plan.U
                          << ",direction="
                          << DirectionName(plan.DirectionValue)
                          << ",mode_position=" << plan.ModePosition
                          << ",mode=" << ModeName(plan.ModeValue)
                          << ",orientation_flip=" << plan.OrientationFlip
                          << ",rotation=" << plan.Rotation
                          << ",sequence_position=" << sequence_position
                          << ",sequence=" << SequenceName(sequence)
                          << ",family=" << (sequence < 2u ? 0u : 1u)
                          << ",step=" << step
                          << ",period=" << period
                          << ",kind=" << (measured ? "M" : "W")
                          << ",label="
                          << (measured ? LabelName(label) : "W")
                          << ",physical=" << BankName(physical)
                          << ",predecessor="
                          << (measured ?
                                PredecessorName(predecessor) : "none")
                          << ",batch_size=" << kBatchSize << '\n';
            }
        }
    }
    if (action_index != kPanelActionCount) return std::string();
    return wirehair::wh2_benchmark::Sha256Hex(canonical.str());
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

Batch8ArmResult RunFixtureBatch8(
    const FixtureBank& bank,
    Scope scope,
    bool measured)
{
    if (scope == Scope::Encoder) {
        return measured ?
            bank.Encoder.RunBatch8() : bank.Encoder.PreflightBatch8();
    }
    if (scope == Scope::IidReceive) {
        return measured ?
            bank.IidReceive.RunBatch8() : bank.IidReceive.PreflightBatch8();
    }
    return measured ?
        bank.RepairReceive.RunBatch8() :
        bank.RepairReceive.PreflightBatch8();
}

bool ValidSemanticBatch8Result(
    Scope scope,
    uint32_t expected_iid_overhead,
    const Batch8ArmResult& result,
    bool measured)
{
    if (measured ? result.ElapsedNanoseconds == 0u :
            result.ElapsedNanoseconds != 0u)
    {
        return false;
    }
    if (wirehair_wh2_bench::kNativeBatch8Size != kBatchSize ||
        result.VerificationMask != ((1u << kBatchSize) - 1u))
    {
        return false;
    }
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
        bank.Arm.BlockCount() != kK ||
        bank.Arm.BlockBytes() != block_bytes)
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
    for (uint32_t cell_index = 0u;
         cell_index < kBaseCellCount;
         ++cell_index)
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
        if (cell->SourceSha256.empty() || cell->IidIdsSha256.empty() ||
            cell->RepairIdsSha256.empty())
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
    size_t preflight_batch_count = 0u;
    for (uint32_t cell_index = 0u;
         cell_index < kBaseCellCount;
         ++cell_index)
    {
        CalibrationCell& cell = *cells[cell_index];
        const Batch8ArmResult x_iid =
            cell.Banks[0].IidReceive.PreflightBatch8();
        ++preflight_batch_count;
        if (x_iid.ElapsedNanoseconds != 0u ||
            x_iid.VerificationMask != ((1u << kBatchSize) - 1u) ||
            x_iid.DecodedOverheads[0] == UINT32_MAX ||
            x_iid.DecodedOverheads[0] > kReceiveOverheadCap)
        {
            diagnostic = "X IID clock-free preflight failed";
            return false;
        }
        cell.IidOverhead = x_iid.DecodedOverheads[0];
        if (cell.IidOverhead != kExpectedIidOverheads[cell_index])
        {
            diagnostic = "IID overhead differs from the frozen calibration";
            return false;
        }

        for (uint32_t bank_index = 0u; bank_index < 3u; ++bank_index)
        {
            const FixtureBank& bank = cell.Banks[bank_index];
            for (uint32_t scope_index = 0u; scope_index < 3u; ++scope_index)
            {
                const Scope scope = static_cast<Scope>(scope_index);
                Batch8ArmResult result;
                if (scope == Scope::IidReceive && bank_index == 0u) {
                    result = x_iid;
                }
                else
                {
                    result = RunFixtureBatch8(bank, scope, false);
                    ++preflight_batch_count;
                }
                if (!ValidSemanticBatch8Result(
                        scope, cell.IidOverhead, result, false))
                {
                    std::ostringstream message;
                    message << "clock-free preflight failed for cell="
                            << cell_index << " bank=" << bank_index
                            << " scope=" << ScopeName(scope)
                            << " result0="
                            << static_cast<uint32_t>(result.Results[0])
                            << " overhead0=" << result.DecodedOverheads[0];
                    diagnostic = message.str();
                    return false;
                }
            }
        }
    }
    if (preflight_batch_count != kUnpinnedPreflightBatchCount)
    {
        diagnostic = "clock-free preflight batch count invalid";
        return false;
    }
    return true;
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
    unsigned selected = 0u;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        selected += CPU_ISSET(cpu, &observed) ? 1u : 0u;
    }
    if (selected != 1u || !CPU_ISSET(target_cpu, &observed))
    {
        diagnostic = "affinity is not the requested singleton CPU";
        return false;
    }
    errno = 0;
    const int observed_cpu = sched_getcpu();
    if (observed_cpu < 0)
    {
        diagnostic = std::string("sched_getcpu failed: ") +
            std::strerror(errno);
        return false;
    }
    if (observed_cpu != target_cpu)
    {
        std::ostringstream message;
        message << "CPU migration from " << target_cpu
                << " to " << observed_cpu;
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

struct ClockDefinition
{
    uint64_t SteadyNumerator = 0u;
    uint64_t SteadyDenominator = 0u;
    uint64_t SteadyCeilNanoseconds = 0u;
    uint64_t MonotonicResolutionNanoseconds = 0u;
    uint64_t RhoNanoseconds = 0u;
};

bool ReadClockDefinition(ClockDefinition& clock, std::string& diagnostic)
{
    if (!std::chrono::steady_clock::is_steady ||
        std::chrono::steady_clock::period::num <= 0 ||
        std::chrono::steady_clock::period::den <= 0)
    {
        diagnostic = "steady_clock does not provide a positive steady period";
        return false;
    }
    ClockDefinition value;
    value.SteadyNumerator = static_cast<uint64_t>(
        std::chrono::steady_clock::period::num);
    value.SteadyDenominator = static_cast<uint64_t>(
        std::chrono::steady_clock::period::den);
    const long double steady_ns =
        static_cast<long double>(value.SteadyNumerator) * 1000000000.0L /
        static_cast<long double>(value.SteadyDenominator);
    const long double steady_ceil = std::ceil(steady_ns);
    if (!std::isfinite(steady_ceil) || steady_ceil < 1.0L ||
        steady_ceil >
            static_cast<long double>(std::numeric_limits<uint64_t>::max()))
    {
        diagnostic = "steady_clock period cannot be represented in ns";
        return false;
    }
    value.SteadyCeilNanoseconds = static_cast<uint64_t>(steady_ceil);
#if defined(__linux__)
    struct timespec resolution;
    if (clock_getres(CLOCK_MONOTONIC, &resolution) != 0 ||
        resolution.tv_sec < 0 || resolution.tv_nsec < 0 ||
        resolution.tv_nsec >= 1000000000L)
    {
        diagnostic = std::string("clock_getres(CLOCK_MONOTONIC) failed: ") +
            std::strerror(errno);
        return false;
    }
    const uint64_t seconds = static_cast<uint64_t>(resolution.tv_sec);
    if (seconds >
        (std::numeric_limits<uint64_t>::max() -
            static_cast<uint64_t>(resolution.tv_nsec)) /
            UINT64_C(1000000000))
    {
        diagnostic = "CLOCK_MONOTONIC resolution overflows ns";
        return false;
    }
    value.MonotonicResolutionNanoseconds =
        seconds * UINT64_C(1000000000) +
        static_cast<uint64_t>(resolution.tv_nsec);
#else
    diagnostic = "CLOCK_MONOTONIC resolution requires Linux";
    return false;
#endif
    if (value.MonotonicResolutionNanoseconds == 0u)
    {
        diagnostic = "CLOCK_MONOTONIC reported zero resolution";
        return false;
    }
    value.RhoNanoseconds = std::max(
        value.SteadyCeilNanoseconds,
        value.MonotonicResolutionNanoseconds);
    clock = value;
    return true;
}

bool ExecuteQualifiedBatch8(
    const CalibrationCell& base_cell,
    const LogicalCellInfo& logical_cell,
    PhysicalBank bank_name,
    bool measured,
    int target_cpu,
    uint64_t& affinity_verification_count,
    Batch8ArmResult& result_out,
    std::string& diagnostic)
{
    ++affinity_verification_count;
    if (!VerifyCpu(target_cpu, diagnostic)) return false;

    Batch8ArmResult result;
    uint32_t invocation_exception_kind = 0u;
    try {
        result = RunFixtureBatch8(
            SelectBank(base_cell, bank_name),
            logical_cell.ScopeValue,
            measured);
    }
    catch (const std::exception& exception)
    {
        (void)exception;
        invocation_exception_kind = 1u;
    }
    catch (...)
    {
        invocation_exception_kind = 2u;
    }

    ++affinity_verification_count;
    const bool post_affinity = VerifyCpu(target_cpu, diagnostic);
    if (!post_affinity) return false;
    if (invocation_exception_kind != 0u)
    {
        std::ostringstream message;
        message << (measured ? "measured" : "clock-free")
                << " invocation threw for cell=" << logical_cell.Index
                << " bank=" << BankName(bank_name)
                << " scope=" << ScopeName(logical_cell.ScopeValue)
                << " exception="
                << (invocation_exception_kind == 1u ?
                    "standard" : "non-standard");
        diagnostic = message.str();
        return false;
    }
    if (!ValidSemanticBatch8Result(
            logical_cell.ScopeValue,
            base_cell.IidOverhead,
            result,
            measured))
    {
        std::ostringstream message;
        message << (measured ? "measured" : "clock-free")
                << " invocation failed for cell=" << logical_cell.Index
                << " bank=" << BankName(bank_name)
                << " scope=" << ScopeName(logical_cell.ScopeValue)
                << " result0=" << static_cast<uint32_t>(result.Results[0])
                << " verification_mask=" << result.VerificationMask
                << " overhead0=" << result.DecodedOverheads[0]
                << " elapsed_ns=" << result.ElapsedNanoseconds;
        diagnostic = message.str();
        return false;
    }
    result_out = result;
    return true;
}

bool RunPinnedWarmups(
    const CellArray& cells,
    int target_cpu,
    uint64_t& affinity_verification_count,
    std::string& diagnostic)
{
    const uint64_t initial_verification_count = affinity_verification_count;
    for (uint32_t cell = 0u; cell < kCellCount; ++cell)
    {
        const LogicalCellInfo logical_cell = GetLogicalCellInfo(cell);
        const CalibrationCell& base_cell = *cells[logical_cell.BaseIndex];
        for (uint32_t bank = 0u; bank < 3u; ++bank)
        {
            Batch8ArmResult discarded;
            if (!ExecuteQualifiedBatch8(
                    base_cell,
                    logical_cell,
                    static_cast<PhysicalBank>(bank),
                    false,
                    target_cpu,
                    affinity_verification_count,
                    discarded,
                    diagnostic))
            {
                return false;
            }
        }
    }
    if (affinity_verification_count - initial_verification_count !=
        static_cast<uint64_t>(kPinnedWarmBatchCount) * 2u)
    {
        diagnostic = "pinned warm affinity-verification count invalid";
        return false;
    }
    return true;
}

struct RawTiming
{
    uint64_t PanelActionIndex = 0u;
    uint32_t Cell = 0u;
    uint32_t Round = 0u;
    uint32_t Visit = 0u;
    uint32_t U = 0u;
    Direction DirectionValue = Direction::Forward;
    Mode ModeValue = Mode::Exposed;
    uint32_t ModePosition = 0u;
    uint32_t OrientationFlip = 0u;
    PhysicalBank LogicalAPhysical = PhysicalBank::X;
    uint32_t Rotation = 0u;
    uint32_t SequencePosition = 0u;
    uint32_t Sequence = 0u;
    uint32_t SequenceFamily = 0u;
    uint32_t Period = 0u;
    LogicalLabel Label = LogicalLabel::A;
    PhysicalBank Physical = PhysicalBank::X;
    PredecessorClass Predecessor = PredecessorClass::Washout;
    Batch8ArmResult Batch;
};

struct SuperblockObservation
{
    SuperblockPlan Plan;
    std::array<uint64_t, 2> LogicalSums = {{ 0u, 0u }};
    std::array<uint64_t, 2> PhysicalSums = {{ 0u, 0u }};
    std::array<uint64_t, 4> PeriodSums = {{ 0u, 0u, 0u, 0u }};
    std::array<uint64_t, 4> SequencePositionSums = {{ 0u, 0u, 0u, 0u }};
    std::array<uint64_t, 4> TransitionSums = {{ 0u, 0u, 0u, 0u }};
    std::array<uint64_t, 4> FactorialSums = {{ 0u, 0u, 0u, 0u }};
    std::array<uint32_t, 2> LogicalCounts = {{ 0u, 0u }};
    std::array<uint32_t, 2> PhysicalCounts = {{ 0u, 0u }};
    std::array<uint32_t, 4> FactorialCounts = {{ 0u, 0u, 0u, 0u }};
    uint64_t WashoutBatches = 0u;
    double LogLogicalAa = 0.0;
    double LogPhysicalXy = 0.0;
    double LogFactorialInteraction = 0.0;
    double LogicalQuantizationBudget = 0.0;
    double PhysicalQuantizationBudget = 0.0;
    bool LogicalQuantizationPass = false;
    bool PhysicalQuantizationPass = false;
};

size_t ObservationKey(uint32_t cell, uint32_t round, Mode mode)
{
    return (static_cast<size_t>(cell) * kRoundCount + round) * kModeCount +
        static_cast<uint32_t>(mode);
}

PhysicalBank PhysicalForLabel(
    PhysicalBank logical_a,
    LogicalLabel label)
{
    if (label == LogicalLabel::A) return logical_a;
    return logical_a == PhysicalBank::X ? PhysicalBank::Y : PhysicalBank::X;
}

bool AccumulateMeasured(
    SuperblockObservation& observation,
    uint32_t sequence_position,
    uint32_t sequence,
    uint32_t period,
    LogicalLabel label,
    PhysicalBank physical,
    uint64_t elapsed)
{
    if (elapsed == 0u || sequence_position >= kSequenceCount ||
        sequence >= kSequenceCount || period >= kPeriodsPerSequence ||
        physical == PhysicalBank::W)
    {
        return false;
    }
    const uint32_t label_index = static_cast<uint32_t>(label);
    const uint32_t physical_index = static_cast<uint32_t>(physical);
    const uint32_t factorial_index = label_index * 2u + physical_index;
    if (!AddElapsed(elapsed, observation.LogicalSums[label_index]) ||
        !AddElapsed(elapsed, observation.PhysicalSums[physical_index]) ||
        !AddElapsed(elapsed, observation.PeriodSums[period]) ||
        !AddElapsed(
            elapsed, observation.SequencePositionSums[sequence_position]) ||
        !AddElapsed(elapsed, observation.FactorialSums[factorial_index]))
    {
        return false;
    }
    ++observation.LogicalCounts[label_index];
    ++observation.PhysicalCounts[physical_index];
    ++observation.FactorialCounts[factorial_index];
    if (observation.Plan.ModeValue == Mode::Exposed && period > 0u)
    {
        const uint32_t previous = static_cast<uint32_t>(
            kSequences[sequence][period - 1u]);
        const uint32_t transition = previous * 2u + label_index;
        if (!AddElapsed(elapsed, observation.TransitionSums[transition])) {
            return false;
        }
    }
    return true;
}

bool QuantizationBudgetValue(
    uint64_t rho_nanoseconds,
    uint64_t sum_a,
    uint64_t sum_b,
    double& budget_out,
    bool& pass_out)
{
    if (rho_nanoseconds == 0u || sum_a == 0u || sum_b == 0u) return false;
    const long double budget =
        8.0L * static_cast<long double>(rho_nanoseconds) /
            static_cast<long double>(sum_a) +
        8.0L * static_cast<long double>(rho_nanoseconds) /
            static_cast<long double>(sum_b);
    if (!std::isfinite(budget) || budget < 0.0L)
    {
        return false;
    }
    budget_out = static_cast<double>(budget);
    pass_out = budget <= kQuantizationBudgetLimit;
    return true;
}

bool QuantizationBudget(
    uint64_t rho_nanoseconds,
    uint64_t sum_a,
    uint64_t sum_b,
    double& budget_out)
{
    bool pass = false;
    return QuantizationBudgetValue(
            rho_nanoseconds, sum_a, sum_b, budget_out, pass) && pass;
}

bool FinalizeObservation(
    uint64_t rho_nanoseconds,
    SuperblockObservation& observation)
{
    const uint64_t sum_a = observation.LogicalSums[0];
    const uint64_t sum_b = observation.LogicalSums[1];
    const uint64_t sum_x = observation.PhysicalSums[0];
    const uint64_t sum_y = observation.PhysicalSums[1];
    if (sum_a == 0u || sum_b == 0u || sum_x == 0u || sum_y == 0u ||
        observation.WashoutBatches != kWashoutBatchesPerSuperblock ||
        observation.LogicalCounts[0] != 8u ||
        observation.LogicalCounts[1] != 8u ||
        observation.PhysicalCounts[0] != 8u ||
        observation.PhysicalCounts[1] != 8u)
    {
        return false;
    }
    for (size_t i = 0u; i < observation.FactorialCounts.size(); ++i) {
        if (observation.FactorialCounts[i] != 4u ||
            observation.FactorialSums[i] == 0u) return false;
    }
    observation.LogLogicalAa = static_cast<double>(std::log(
        static_cast<long double>(sum_a) / sum_b));
    observation.LogPhysicalXy = static_cast<double>(std::log(
        static_cast<long double>(sum_x) / sum_y));
    observation.LogFactorialInteraction = static_cast<double>(
        std::log(static_cast<long double>(observation.FactorialSums[0])) +
        std::log(static_cast<long double>(observation.FactorialSums[3])) -
        std::log(static_cast<long double>(observation.FactorialSums[1])) -
        std::log(static_cast<long double>(observation.FactorialSums[2])));
    if (!std::isfinite(observation.LogLogicalAa) ||
        !std::isfinite(observation.LogPhysicalXy) ||
        !std::isfinite(observation.LogFactorialInteraction))
    {
        return false;
    }
    if (!QuantizationBudgetValue(
            rho_nanoseconds,
            sum_a,
            sum_b,
            observation.LogicalQuantizationBudget,
            observation.LogicalQuantizationPass) ||
        !QuantizationBudgetValue(
            rho_nanoseconds,
            sum_x,
            sum_y,
            observation.PhysicalQuantizationBudget,
            observation.PhysicalQuantizationPass))
    {
        return false;
    }
    return true;
}

bool RunOneSuperblock(
    const SuperblockPlan& plan,
    const CalibrationCell& base_cell,
    const LogicalCellInfo& logical_cell,
    int target_cpu,
    uint64_t& affinity_verification_count,
    std::vector<RawTiming>& raw,
    size_t& raw_next,
    size_t& panel_action_next,
    size_t& washout_batch_count,
    std::string& diagnostic)
{
    const size_t initial_raw_next = raw_next;
    const size_t initial_panel_action_next = panel_action_next;
    const size_t initial_washout_batch_count = washout_batch_count;

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
                (measured ? step - kPeriodsPerSequence : step) : step / 2u;
            if (!measured)
            {
                Batch8ArmResult discarded;
                if (!ExecuteQualifiedBatch8(
                        base_cell,
                        logical_cell,
                        PhysicalBank::W,
                        false,
                        target_cpu,
                        affinity_verification_count,
                        discarded,
                        diagnostic))
                {
                    return false;
                }
                ++panel_action_next;
                ++washout_batch_count;
                continue;
            }
            const LogicalLabel label = kSequences[sequence][period];
            const PhysicalBank physical = label == LogicalLabel::A ?
                logical_a : logical_b;
            Batch8ArmResult result;
            if (!ExecuteQualifiedBatch8(
                    base_cell,
                    logical_cell,
                    physical,
                    true,
                    target_cpu,
                    affinity_verification_count,
                    result,
                    diagnostic))
            {
                return false;
            }
            if (raw_next >= raw.size())
            {
                diagnostic = "raw receipt overflow";
                return false;
            }
            RawTiming& receipt = raw[raw_next++];
            receipt.PanelActionIndex = panel_action_next++;
            receipt.Cell = plan.Cell;
            receipt.Round = plan.Round;
            receipt.Visit = plan.Visit;
            receipt.U = plan.U;
            receipt.DirectionValue = plan.DirectionValue;
            receipt.ModeValue = plan.ModeValue;
            receipt.ModePosition = plan.ModePosition;
            receipt.OrientationFlip = plan.OrientationFlip;
            receipt.LogicalAPhysical = logical_a;
            receipt.Rotation = plan.Rotation;
            receipt.SequencePosition = sequence_position;
            receipt.Sequence = sequence;
            receipt.SequenceFamily = sequence < 2u ? 0u : 1u;
            receipt.Period = period;
            receipt.Label = label;
            receipt.Physical = physical;
            if (plan.ModeValue == Mode::Isolated) {
                receipt.Predecessor = PredecessorClass::Washout;
            }
            else if (period == 0u) {
                receipt.Predecessor = PredecessorClass::Washout;
            }
            else
            {
                const uint32_t previous = static_cast<uint32_t>(
                    kSequences[sequence][period - 1u]);
                const uint32_t current = static_cast<uint32_t>(label);
                receipt.Predecessor = static_cast<PredecessorClass>(
                    2u + previous * 2u + current);
            }
            receipt.Batch = result;
        }
    }
    if (raw_next - initial_raw_next != kMeasuredBatchesPerSuperblock ||
        panel_action_next - initial_panel_action_next !=
            kMeasuredBatchesPerSuperblock + kWashoutBatchesPerSuperblock ||
        washout_batch_count - initial_washout_batch_count !=
            kWashoutBatchesPerSuperblock)
    {
        diagnostic = "superblock batch count invalid";
        return false;
    }
    return true;
}

bool RunFrozenPanel(
    const CellArray& cells,
    const FrozenSchedule& schedule,
    int target_cpu,
    uint64_t& affinity_verification_count,
    std::vector<RawTiming>& raw,
    std::string& diagnostic)
{
    const uint64_t initial_verification_count = affinity_verification_count;
    size_t raw_next = 0u;
    size_t panel_action_next = 0u;
    size_t washout_batch_count = 0u;
    if (raw.size() != kMeasuredBatchCount)
    {
        diagnostic = "preallocated receipt cardinality is invalid";
        return false;
    }
    for (size_t i = 0u; i < schedule.size(); ++i)
    {
        const SuperblockPlan& plan = schedule[i];
        if (plan.Cell >= kCellCount)
        {
            diagnostic = "schedule references an invalid cell";
            return false;
        }
        const LogicalCellInfo logical_cell = GetLogicalCellInfo(plan.Cell);
        if (logical_cell.BaseIndex >= cells.size() ||
            !cells[logical_cell.BaseIndex])
        {
            diagnostic = "schedule references an invalid base cell";
            return false;
        }
        if (!RunOneSuperblock(
                plan,
                *cells[logical_cell.BaseIndex],
                logical_cell,
                target_cpu,
                affinity_verification_count,
                raw,
                raw_next,
                panel_action_next,
                washout_batch_count,
                diagnostic))
        {
            return false;
        }
    }
    if (raw_next != raw.size() || panel_action_next != kPanelActionCount ||
        washout_batch_count != kWashoutBatchCount)
    {
        diagnostic = "final measured or washout count invalid";
        return false;
    }
    const uint64_t expected_panel_verifications =
        static_cast<uint64_t>(
            kMeasuredBatchCount + kWashoutBatchCount) * 2u;
    if (affinity_verification_count - initial_verification_count !=
        expected_panel_verifications)
    {
        diagnostic = "panel affinity-verification count invalid";
        return false;
    }
    return true;
}

bool BuildObservations(
    const FrozenSchedule& schedule,
    const std::vector<RawTiming>& raw,
    uint64_t rho_nanoseconds,
    std::vector<SuperblockObservation>& observations,
    std::vector<size_t>& observation_by_key,
    std::string& diagnostic)
{
    if (raw.size() != kMeasuredBatchCount ||
        observations.size() != kSuperblockCount ||
        observation_by_key.size() != kSuperblockCount)
    {
        diagnostic = "post-panel observation cardinality is invalid";
        return false;
    }
    std::fill(
        observation_by_key.begin(),
        observation_by_key.end(),
        std::numeric_limits<size_t>::max());

    size_t raw_next = 0u;
    for (size_t i = 0u; i < schedule.size(); ++i)
    {
        const SuperblockPlan& plan = schedule[i];
        SuperblockObservation& observation = observations[i];
        observation = SuperblockObservation();
        observation.Plan = plan;
        observation.WashoutBatches = kWashoutBatchesPerSuperblock;
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
            for (uint32_t period = 0u;
                 period < kPeriodsPerSequence;
                 ++period)
            {
                if (raw_next >= raw.size())
                {
                    diagnostic = "post-panel raw receipt underflow";
                    return false;
                }
                const RawTiming& receipt = raw[raw_next++];
                const LogicalLabel label = kSequences[sequence][period];
                const PhysicalBank physical = label == LogicalLabel::A ?
                    logical_a : logical_b;
                const uint64_t expected_action = static_cast<uint64_t>(
                    i * (kMeasuredBatchesPerSuperblock +
                        kWashoutBatchesPerSuperblock) +
                    sequence_position * kBatchesPerSequence +
                    (plan.ModeValue == Mode::Exposed ?
                        kPeriodsPerSequence + period : 2u * period + 1u));
                PredecessorClass predecessor = PredecessorClass::Washout;
                if (plan.ModeValue == Mode::Isolated) {
                    predecessor = PredecessorClass::Washout;
                }
                else if (period > 0u)
                {
                    const uint32_t previous = static_cast<uint32_t>(
                        kSequences[sequence][period - 1u]);
                    predecessor = static_cast<PredecessorClass>(
                        2u + previous * 2u + static_cast<uint32_t>(label));
                }
                const LogicalCellInfo info = GetLogicalCellInfo(plan.Cell);
                if (receipt.PanelActionIndex != expected_action ||
                    receipt.Cell != plan.Cell ||
                    receipt.Round != plan.Round ||
                    receipt.Visit != plan.Visit ||
                    receipt.U != plan.U ||
                    receipt.DirectionValue != plan.DirectionValue ||
                    receipt.ModeValue != plan.ModeValue ||
                    receipt.ModePosition != plan.ModePosition ||
                    receipt.OrientationFlip != plan.OrientationFlip ||
                    receipt.LogicalAPhysical != logical_a ||
                    receipt.Rotation != plan.Rotation ||
                    receipt.SequencePosition != sequence_position ||
                    receipt.Sequence != sequence ||
                    receipt.SequenceFamily != (sequence < 2u ? 0u : 1u) ||
                    receipt.Period != period ||
                    receipt.Label != label ||
                    receipt.Physical != physical ||
                    receipt.Predecessor != predecessor ||
                    !ValidSemanticBatch8Result(
                        info.ScopeValue,
                        kExpectedIidOverheads[info.BaseIndex],
                        receipt.Batch,
                        true) ||
                    !AccumulateMeasured(
                        observation,
                        sequence_position,
                        sequence,
                        period,
                        label,
                        physical,
                        receipt.Batch.ElapsedNanoseconds))
                {
                    diagnostic = "post-panel raw receipt semantic mismatch";
                    return false;
                }
            }
        }
        if (!FinalizeObservation(rho_nanoseconds, observation))
        {
            diagnostic = "post-panel observation arithmetic failed";
            return false;
        }
        const size_t key = ObservationKey(
            plan.Cell, plan.Round, plan.ModeValue);
        if (key >= observation_by_key.size() ||
            observation_by_key[key] !=
                std::numeric_limits<size_t>::max())
        {
            diagnostic = "duplicate or invalid observation key";
            return false;
        }
        observation_by_key[key] = i;
    }
    if (raw_next != raw.size() ||
        std::find(
            observation_by_key.begin(),
            observation_by_key.end(),
            std::numeric_limits<size_t>::max()) != observation_by_key.end())
    {
        diagnostic = "post-panel observation coverage is invalid";
        return false;
    }
    return true;
}

struct LogStats
{
    size_t Count = 0u;
    double GeometricMean = 0.0;
    double Lower95 = 0.0;
    double Upper95 = 0.0;
};

double StudentT975(size_t count)
{
    if (count == 32u) return 2.039513446;
    if (count == 96u) return 1.985251004;
    if (count == 576u) return 1.964095;
    return 0.0;
}

bool ComputeLogStats(
    const double* values,
    size_t count,
    LogStats& stats_out)
{
    const double t_value = StudentT975(count);
    if (!values || count < 2u || t_value <= 0.0) return false;
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
        static_cast<long double>(t_value) * standard_error;
    LogStats stats;
    stats.Count = count;
    stats.GeometricMean = static_cast<double>(std::exp(mean));
    stats.Lower95 = static_cast<double>(std::exp(mean - half_width));
    stats.Upper95 = static_cast<double>(std::exp(mean + half_width));
    if (!std::isfinite(stats.GeometricMean) ||
        !std::isfinite(stats.Lower95) ||
        !std::isfinite(stats.Upper95) || stats.Lower95 <= 0.0 ||
        stats.Lower95 > stats.GeometricMean ||
        stats.GeometricMean > stats.Upper95)
    {
        return false;
    }
    stats_out = stats;
    return true;
}

bool OfficialGate(const LogStats& stats)
{
    return stats.Lower95 >= 0.98 && stats.Lower95 <= 1.0 &&
        stats.Upper95 >= 1.0 && stats.Upper95 <= 1.02;
}

const SuperblockObservation* FindObservation(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    uint32_t cell,
    uint32_t round,
    Mode mode)
{
    const size_t key = ObservationKey(cell, round, mode);
    if (key >= by_key.size() || by_key[key] >= observations.size()) {
        return nullptr;
    }
    return &observations[by_key[key]];
}

struct OfficialResult
{
    Estimand EstimandValue = Estimand::Logical;
    Mode ModeValue = Mode::Exposed;
    bool Pooled = false;
    uint32_t WidthIndex = 0u;
    Scope ScopeValue = Scope::Encoder;
    LogStats Stats;
    bool Pass = false;
};

bool BuildOfficialResults(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    std::array<OfficialResult, kOfficialMarginalGateCount>& results,
    bool& all_pass)
{
    size_t next = 0u;
    all_pass = true;
    std::array<double, 576> logs = {{ 0.0 }};
    for (uint32_t estimand_value = 0u; estimand_value < 2u; ++estimand_value)
    {
        const Estimand estimand = static_cast<Estimand>(estimand_value);
        for (uint32_t mode_value = 0u; mode_value < kModeCount; ++mode_value)
        {
            const Mode mode = static_cast<Mode>(mode_value);
            for (uint32_t width = 0u; width < 2u; ++width)
            {
                for (uint32_t scope_value = 0u;
                     scope_value < 3u;
                     ++scope_value)
                {
                    size_t count = 0u;
                    for (uint32_t root = 0u; root < 3u; ++root)
                    {
                        const uint32_t cell =
                            width * 9u + scope_value * 3u + root;
                        for (uint32_t round = 0u;
                             round < kRoundCount;
                             ++round)
                        {
                            const SuperblockObservation* observation =
                                FindObservation(
                                    observations, by_key, cell, round, mode);
                            if (!observation || count >= 96u) return false;
                            logs[count++] = estimand == Estimand::Logical ?
                                observation->LogLogicalAa :
                                observation->LogPhysicalXy;
                        }
                    }
                    OfficialResult result;
                    result.EstimandValue = estimand;
                    result.ModeValue = mode;
                    result.WidthIndex = width;
                    result.ScopeValue = static_cast<Scope>(scope_value);
                    if (count != 96u ||
                        !ComputeLogStats(logs.data(), count, result.Stats))
                    {
                        return false;
                    }
                    result.Pass = OfficialGate(result.Stats);
                    all_pass = all_pass && result.Pass;
                    results[next++] = result;
                }
            }

            size_t pooled_count = 0u;
            for (uint32_t cell = 0u; cell < kCellCount; ++cell)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round)
                {
                    const SuperblockObservation* observation = FindObservation(
                        observations, by_key, cell, round, mode);
                    if (!observation || pooled_count >= logs.size()) {
                        return false;
                    }
                    logs[pooled_count++] = estimand == Estimand::Logical ?
                        observation->LogLogicalAa :
                        observation->LogPhysicalXy;
                }
            }
            OfficialResult pooled;
            pooled.EstimandValue = estimand;
            pooled.ModeValue = mode;
            pooled.Pooled = true;
            if (pooled_count != logs.size() ||
                !ComputeLogStats(logs.data(), pooled_count, pooled.Stats))
            {
                return false;
            }
            pooled.Pass = OfficialGate(pooled.Stats);
            all_pass = all_pass && pooled.Pass;
            results[next++] = pooled;
        }
    }
    return next == results.size();
}

void PrintStatsLine(
    const char* record,
    const char* diagnostic,
    Mode mode,
    const LogStats& stats,
    uint32_t block_bytes = 0u,
    Scope scope = Scope::Encoder,
    uint64_t root = 0u,
    const char* gate = nullptr,
    const char* mode_override = nullptr)
{
    std::cout << std::setprecision(12) << record;
    if (diagnostic) std::cout << ",diagnostic=" << diagnostic;
    std::cout << ",mode=" <<
        (mode_override ? mode_override : ModeName(mode));
    if (block_bytes != 0u) {
        std::cout << ",block_bytes=" << block_bytes
                  << ",scope=" << ScopeName(scope);
    }
    if (root != 0u) {
        std::cout << ",root=0x" << std::hex << root << std::dec;
    }
    std::cout << ",n=" << stats.Count
              << ",gm=" << stats.GeometricMean
              << ",lower95=" << stats.Lower95
              << ",upper95=" << stats.Upper95;
    if (gate) std::cout << ",gate=" << gate;
    std::cout << '\n';
}

bool RatioLog(uint64_t numerator, uint64_t denominator, double& value_out)
{
    if (numerator == 0u || denominator == 0u) return false;
    const long double value = std::log(
        static_cast<long double>(numerator) / denominator);
    if (!std::isfinite(value)) return false;
    value_out = static_cast<double>(value);
    return true;
}

double ObservationLog(
    const SuperblockObservation& observation,
    Estimand estimand)
{
    return estimand == Estimand::Logical ?
        observation.LogLogicalAa : observation.LogPhysicalXy;
}

bool PeriodDiagnosticLog(
    const SuperblockObservation& observation,
    uint32_t metric,
    double& value_out)
{
    uint64_t numerator = 0u;
    uint64_t denominator = 0u;
    if (metric == 0u) {
        numerator = observation.PeriodSums[0];
        denominator = observation.PeriodSums[3];
    }
    else if (metric == 1u) {
        numerator = observation.PeriodSums[1];
        denominator = observation.PeriodSums[2];
    }
    else if (metric == 2u)
    {
        if (!AddElapsed(observation.PeriodSums[0], numerator) ||
            !AddElapsed(observation.PeriodSums[3], numerator) ||
            !AddElapsed(observation.PeriodSums[1], denominator) ||
            !AddElapsed(observation.PeriodSums[2], denominator))
        {
            return false;
        }
    }
    else {
        return false;
    }
    return RatioLog(numerator, denominator, value_out);
}

bool SequencePositionDiagnosticLog(
    const SuperblockObservation& observation,
    uint32_t metric,
    double& value_out)
{
    uint64_t numerator = 0u;
    uint64_t denominator = 0u;
    if (metric == 0u) {
        numerator = observation.SequencePositionSums[0];
        denominator = observation.SequencePositionSums[3];
    }
    else if (metric == 1u) {
        numerator = observation.SequencePositionSums[1];
        denominator = observation.SequencePositionSums[2];
    }
    else if (metric == 2u)
    {
        if (!AddElapsed(observation.SequencePositionSums[0], numerator) ||
            !AddElapsed(observation.SequencePositionSums[3], numerator) ||
            !AddElapsed(observation.SequencePositionSums[1], denominator) ||
            !AddElapsed(observation.SequencePositionSums[2], denominator))
        {
            return false;
        }
    }
    else return false;
    return RatioLog(numerator, denominator, value_out);
}

bool TransitionDiagnosticLog(
    const SuperblockObservation& observation,
    uint32_t metric,
    double& value_out)
{
    const uint64_t aa = observation.TransitionSums[0];
    const uint64_t ab = observation.TransitionSums[1];
    const uint64_t ba = observation.TransitionSums[2];
    const uint64_t bb = observation.TransitionSums[3];
    uint64_t numerator = 0u;
    uint64_t denominator = 0u;
    if (metric == 0u)
    {
        if (!AddElapsed(ab, numerator) || !AddElapsed(ba, numerator) ||
            !AddElapsed(aa, denominator) || !AddElapsed(bb, denominator))
        {
            return false;
        }
    }
    else if (metric == 1u) {
        numerator = ab;
        denominator = ba;
    }
    else if (metric == 2u) {
        numerator = aa;
        denominator = bb;
    }
    else {
        return false;
    }
    return RatioLog(numerator, denominator, value_out);
}

bool PrintHalfDiagnostic(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    Estimand estimand,
    Mode mode,
    int logical_cell)
{
    long double first_sum = 0.0L;
    long double second_sum = 0.0L;
    size_t first_count = 0u;
    size_t second_count = 0u;
    const uint32_t first_cell = logical_cell < 0 ?
        0u : static_cast<uint32_t>(logical_cell);
    const uint32_t end_cell = logical_cell < 0 ?
        kCellCount : first_cell + 1u;
    for (uint32_t cell = first_cell; cell < end_cell; ++cell)
    {
        for (uint32_t round = 0u; round < kRoundCount; ++round)
        {
            const SuperblockObservation* observation = FindObservation(
                observations, by_key, cell, round, mode);
            if (!observation) return false;
            if (round < 16u) {
                first_sum += ObservationLog(*observation, estimand);
                ++first_count;
            }
            else {
                second_sum += ObservationLog(*observation, estimand);
                ++second_count;
            }
        }
    }
    if (first_count == 0u || first_count != second_count) return false;
    const long double interaction = std::exp(
        first_sum / static_cast<long double>(first_count) -
        second_sum / static_cast<long double>(second_count));
    if (!std::isfinite(interaction) || interaction <= 0.0L) return false;

    std::cout << std::setprecision(12)
              << "DIAGNOSTIC,diagnostic=rounds-0-15-vs-16-31"
              << ",estimand=" << EstimandName(estimand)
              << ",mode=" << ModeName(mode);
    if (logical_cell >= 0)
    {
        const LogicalCellInfo info = GetLogicalCellInfo(
            static_cast<uint32_t>(logical_cell));
        std::cout << ",block_bytes=" << info.BlockBytes
                  << ",scope=" << ScopeName(info.ScopeValue)
                  << ",root=0x" << std::hex << info.Root << std::dec;
    }
    else {
        std::cout << ",pool=all-cells";
    }
    std::cout << ",n_first=" << first_count
              << ",n_second=" << second_count
              << ",gm_interaction=" << static_cast<double>(interaction)
              << ",alias=rotation-q-plus-2,gate=NONE\n";
    return true;
}

bool PrintDirectionDiagnostic(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    Estimand estimand,
    Mode mode,
    int logical_cell)
{
    long double forward = 0.0L;
    long double reverse = 0.0L;
    size_t forward_count = 0u;
    size_t reverse_count = 0u;
    const uint32_t first_cell = logical_cell < 0 ?
        0u : static_cast<uint32_t>(logical_cell);
    const uint32_t end_cell = logical_cell < 0 ?
        kCellCount : first_cell + 1u;
    for (uint32_t cell = first_cell; cell < end_cell; ++cell)
    {
        for (uint32_t round = 0u; round < kRoundCount; ++round)
        {
            const SuperblockObservation* observation = FindObservation(
                observations, by_key, cell, round, mode);
            if (!observation) return false;
            if (observation->Plan.DirectionValue == Direction::Forward) {
                forward += ObservationLog(*observation, estimand);
                ++forward_count;
            }
            else {
                reverse += ObservationLog(*observation, estimand);
                ++reverse_count;
            }
        }
    }
    if (forward_count == 0u || forward_count != reverse_count) return false;
    const long double interaction = std::exp(
        forward / static_cast<long double>(forward_count) -
        reverse / static_cast<long double>(reverse_count));
    if (!std::isfinite(interaction) || interaction <= 0.0L) return false;
    std::cout << std::setprecision(12)
              << "DIAGNOSTIC,diagnostic=forward-vs-reverse"
              << ",estimand=" << EstimandName(estimand)
              << ",mode=" << ModeName(mode);
    if (logical_cell >= 0)
    {
        const LogicalCellInfo info = GetLogicalCellInfo(
            static_cast<uint32_t>(logical_cell));
        std::cout << ",block_bytes=" << info.BlockBytes
                  << ",scope=" << ScopeName(info.ScopeValue)
                  << ",root=0x" << std::hex << info.Root << std::dec;
    }
    else std::cout << ",pool=all-cells";
    std::cout << ",n_forward=" << forward_count
              << ",n_reverse=" << reverse_count
              << ",gm_interaction=" << static_cast<double>(interaction)
              << ",gate=NONE\n";
    return true;
}

bool PrintPooledDiagnostics(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key)
{
    std::array<double, 576> logs = {{ 0.0 }};
    const char* period_names[] = {
        "pooled-period-1-vs-4",
        "pooled-period-2-vs-3",
        "pooled-period-edge-vs-middle"
    };
    for (uint32_t mode_value = 0u; mode_value < kModeCount; ++mode_value)
    {
        const Mode mode = static_cast<Mode>(mode_value);
        size_t next = 0u;
        for (uint32_t cell = 0u; cell < kCellCount; ++cell)
        {
            for (uint32_t estimand = 0u; estimand < 2u; ++estimand)
            {
                if (!PrintHalfDiagnostic(
                        observations,
                        by_key,
                        static_cast<Estimand>(estimand),
                        mode,
                        static_cast<int>(cell)) ||
                    !PrintDirectionDiagnostic(
                        observations,
                        by_key,
                        static_cast<Estimand>(estimand),
                        mode,
                        static_cast<int>(cell)))
                {
                    return false;
                }
            }
            for (uint32_t round = 0u; round < kRoundCount; ++round)
            {
                const SuperblockObservation* observation = FindObservation(
                    observations, by_key, cell, round, mode);
                if (!observation || next >= logs.size()) return false;
                logs[next++] = observation->LogLogicalAa;
            }
        }
        LogStats stats;
        if (next != logs.size() ||
            !ComputeLogStats(logs.data(), next, stats))
        {
            return false;
        }
        PrintStatsLine(
            "DIAGNOSTIC", "pooled-logical-a-b", mode, stats);

        next = 0u;
        for (uint32_t cell = 0u; cell < kCellCount; ++cell)
        {
            for (uint32_t round = 0u; round < kRoundCount; ++round)
            {
                const SuperblockObservation* observation = FindObservation(
                    observations, by_key, cell, round, mode);
                if (!observation || next >= logs.size()) return false;
                logs[next++] = observation->LogPhysicalXy;
            }
        }
        if (!ComputeLogStats(logs.data(), next, stats)) return false;
        PrintStatsLine(
            "DIAGNOSTIC", "pooled-physical-x-y", mode, stats);

        next = 0u;
        for (uint32_t cell = 0u; cell < kCellCount; ++cell) {
            for (uint32_t round = 0u; round < kRoundCount; ++round) {
                const SuperblockObservation* observation = FindObservation(
                    observations, by_key, cell, round, mode);
                if (!observation || next >= logs.size()) return false;
                logs[next++] = observation->LogFactorialInteraction;
            }
        }
        if (!ComputeLogStats(logs.data(), next, stats)) return false;
        PrintStatsLine(
            "DIAGNOSTIC", "pooled-label-by-physical-interaction", mode,
            stats);

        for (uint32_t metric = 0u; metric < 3u; ++metric)
        {
            next = 0u;
            for (uint32_t cell = 0u; cell < kCellCount; ++cell)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round)
                {
                    const SuperblockObservation* observation = FindObservation(
                        observations, by_key, cell, round, mode);
                    if (!observation || next >= logs.size() ||
                        !PeriodDiagnosticLog(
                            *observation, metric, logs[next++]))
                    {
                        return false;
                    }
                }
            }
            if (!ComputeLogStats(logs.data(), next, stats)) return false;
            PrintStatsLine(
                "DIAGNOSTIC", period_names[metric], mode, stats);
        }
        const char* sequence_position_names[] = {
            "pooled-sequence-position-1-vs-4",
            "pooled-sequence-position-2-vs-3",
            "pooled-sequence-position-edge-vs-middle"
        };
        for (uint32_t metric = 0u; metric < 3u; ++metric)
        {
            next = 0u;
            for (uint32_t cell = 0u; cell < kCellCount; ++cell) {
                for (uint32_t round = 0u; round < kRoundCount; ++round) {
                    const SuperblockObservation* observation = FindObservation(
                        observations, by_key, cell, round, mode);
                    if (!observation || next >= logs.size() ||
                        !SequencePositionDiagnosticLog(
                            *observation, metric, logs[next++])) return false;
                }
            }
            if (!ComputeLogStats(logs.data(), next, stats)) return false;
            PrintStatsLine(
                "DIAGNOSTIC", sequence_position_names[metric], mode, stats);
        }
        for (uint32_t estimand = 0u; estimand < 2u; ++estimand)
        {
            if (!PrintHalfDiagnostic(
                    observations,
                    by_key,
                    static_cast<Estimand>(estimand),
                    mode,
                    -1) ||
                !PrintDirectionDiagnostic(
                    observations,
                    by_key,
                    static_cast<Estimand>(estimand),
                    mode,
                    -1)) return false;
        }
    }

    const char* paired_names[] = {
        "pooled-combined-mode",
        "pooled-exposed-minus-isolated",
        "pooled-first-minus-second-superblock"
    };
    for (uint32_t estimand_value = 0u; estimand_value < 2u; ++estimand_value)
    {
        const Estimand estimand = static_cast<Estimand>(estimand_value);
        for (uint32_t metric = 0u; metric < 3u; ++metric)
        {
            size_t next = 0u;
            for (uint32_t cell = 0u; cell < kCellCount; ++cell)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round)
                {
                    const SuperblockObservation* exposed = FindObservation(
                        observations, by_key, cell, round, Mode::Exposed);
                    const SuperblockObservation* isolated = FindObservation(
                        observations, by_key, cell, round, Mode::Isolated);
                    if (!exposed || !isolated || next >= logs.size()) {
                        return false;
                    }
                    const double exposed_log =
                        ObservationLog(*exposed, estimand);
                    const double isolated_log =
                        ObservationLog(*isolated, estimand);
                    if (metric == 0u) {
                        logs[next++] = (exposed_log + isolated_log) / 2.0;
                    }
                    else if (metric == 1u) {
                        logs[next++] = exposed_log - isolated_log;
                    }
                    else
                    {
                        const double first =
                            exposed->Plan.ModePosition == 0u ?
                                exposed_log : isolated_log;
                        const double second =
                            exposed->Plan.ModePosition == 1u ?
                                exposed_log : isolated_log;
                        logs[next++] = first - second;
                    }
                }
            }
            LogStats stats;
            if (next != logs.size() ||
                !ComputeLogStats(logs.data(), next, stats)) return false;
            std::string diagnostic = paired_names[metric];
            diagnostic += estimand == Estimand::Logical ?
                "-logical" : "-physical";
            PrintStatsLine(
                "DIAGNOSTIC", diagnostic.c_str(), Mode::Exposed, stats,
                0u, Scope::Encoder, 0u, nullptr, "paired");
        }
    }

    const char* transition_names[] = {
        "pooled-predecessor-cross-vs-same",
        "pooled-transition-ab-vs-ba",
        "pooled-transition-aa-vs-bb"
    };
    for (uint32_t metric = 0u; metric < 3u; ++metric)
    {
        size_t next = 0u;
        for (uint32_t cell = 0u; cell < kCellCount; ++cell)
        {
            for (uint32_t round = 0u; round < kRoundCount; ++round)
            {
                const SuperblockObservation* exposed = FindObservation(
                    observations, by_key, cell, round, Mode::Exposed);
                if (!exposed || !TransitionDiagnosticLog(
                        *exposed, metric, logs[next++]))
                {
                    return false;
                }
            }
        }
        LogStats stats;
        if (next != logs.size() ||
            !ComputeLogStats(logs.data(), next, stats))
        {
            return false;
        }
        PrintStatsLine(
            "DIAGNOSTIC", transition_names[metric], Mode::Exposed, stats);
    }
    return true;
}

bool PrintQuantizationDiagnostics(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key,
    bool& all_pass)
{
    all_pass = true;
    for (uint32_t estimand_value = 0u; estimand_value < 2u; ++estimand_value)
    {
        const Estimand estimand = static_cast<Estimand>(estimand_value);
        for (uint32_t mode_value = 0u;
             mode_value < kModeCount;
             ++mode_value)
        {
            const Mode mode = static_cast<Mode>(mode_value);
            double pooled_maximum = 0.0;
            size_t pooled_count = 0u;
            size_t pooled_failures = 0u;
            for (uint32_t width = 0u; width < 2u; ++width)
            {
                for (uint32_t scope = 0u; scope < 3u; ++scope)
                {
                    double maximum = 0.0;
                    size_t count = 0u;
                    size_t failures = 0u;
                    for (uint32_t root = 0u; root < 3u; ++root)
                    {
                        const uint32_t cell = width * 9u + scope * 3u + root;
                        for (uint32_t round = 0u;
                             round < kRoundCount;
                             ++round)
                        {
                            const SuperblockObservation* observation =
                                FindObservation(
                                    observations, by_key, cell, round, mode);
                            if (!observation) return false;
                            const double budget =
                                estimand == Estimand::Logical ?
                                    observation->LogicalQuantizationBudget :
                                    observation->PhysicalQuantizationBudget;
                            const bool pass =
                                estimand == Estimand::Logical ?
                                    observation->LogicalQuantizationPass :
                                    observation->PhysicalQuantizationPass;
                            if (!std::isfinite(budget) || budget < 0.0) {
                                return false;
                            }
                            maximum = std::max(maximum, budget);
                            pooled_maximum = std::max(pooled_maximum, budget);
                            ++count;
                            ++pooled_count;
                            if (!pass) {
                                ++failures;
                                ++pooled_failures;
                            }
                        }
                    }
                    if (count != 96u) return false;
                    all_pass = all_pass && failures == 0u;
                    std::cout << std::setprecision(12)
                              << "QUANTIZATION,estimand="
                              << EstimandName(estimand)
                              << ",mode=" << ModeName(mode)
                              << ",block_bytes=" << kBlockBytes[width]
                              << ",scope="
                              << ScopeName(static_cast<Scope>(scope))
                              << ",n=" << count
                              << ",maximum_budget=" << maximum
                              << ",limit="
                              << static_cast<double>(
                                    kQuantizationBudgetLimit)
                              << ",failures=" << failures
                              << ",gate=" << (failures == 0u ?
                                    "PASS" : "INVALID") << '\n';
                }
            }
            if (pooled_count != 576u) return false;
            std::cout << std::setprecision(12)
                      << "QUANTIZATION,estimand=" << EstimandName(estimand)
                      << ",mode=" << ModeName(mode)
                      << ",pool=all-cells,n=" << pooled_count
                      << ",maximum_budget=" << pooled_maximum
                      << ",limit="
                      << static_cast<double>(kQuantizationBudgetLimit)
                      << ",failures=" << pooled_failures
                      << ",gate=" << (pooled_failures == 0u ?
                            "PASS" : "INVALID") << '\n';
        }
    }
    return true;
}

bool PrintMandatoryDiagnostics(
    const std::vector<SuperblockObservation>& observations,
    const std::vector<size_t>& by_key)
{
    std::array<double, 32> logs = {{ 0.0 }};
    for (uint32_t cell_index = 0u; cell_index < kCellCount; ++cell_index)
    {
        const LogicalCellInfo info = GetLogicalCellInfo(cell_index);
        for (uint32_t mode_value = 0u; mode_value < kModeCount; ++mode_value)
        {
            const Mode mode = static_cast<Mode>(mode_value);
            for (uint32_t round = 0u; round < kRoundCount; ++round)
            {
                const SuperblockObservation* observation = FindObservation(
                    observations, by_key, cell_index, round, mode);
                if (!observation) return false;
                logs[round] = observation->LogLogicalAa;
            }
            LogStats stats;
            if (!ComputeLogStats(logs.data(), logs.size(), stats)) return false;
            PrintStatsLine(
                "DIAGNOSTIC", "root-logical-a-b", mode, stats,
                info.BlockBytes, info.ScopeValue, info.Root);

            for (uint32_t round = 0u; round < kRoundCount; ++round) {
                logs[round] = FindObservation(
                    observations, by_key, cell_index, round, mode)->
                        LogPhysicalXy;
            }
            if (!ComputeLogStats(logs.data(), logs.size(), stats)) return false;
            PrintStatsLine(
                "DIAGNOSTIC", "physical-x-y", mode, stats,
                info.BlockBytes, info.ScopeValue, info.Root);

            const char* period_names[] = {
                "period-1-vs-4", "period-2-vs-3", "period-edge-vs-middle"
            };
            for (uint32_t period_metric = 0u;
                 period_metric < 3u;
                 ++period_metric)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round)
                {
                    const SuperblockObservation* observation = FindObservation(
                        observations, by_key, cell_index, round, mode);
                    uint64_t numerator = 0u;
                    uint64_t denominator = 0u;
                    if (period_metric == 0u) {
                        numerator = observation->PeriodSums[0];
                        denominator = observation->PeriodSums[3];
                    }
                    else if (period_metric == 1u) {
                        numerator = observation->PeriodSums[1];
                        denominator = observation->PeriodSums[2];
                    }
                    else if (!AddElapsed(
                                observation->PeriodSums[0], numerator) ||
                             !AddElapsed(
                                observation->PeriodSums[3], numerator) ||
                             !AddElapsed(
                                observation->PeriodSums[1], denominator) ||
                             !AddElapsed(
                                observation->PeriodSums[2], denominator))
                    {
                        return false;
                    }
                    if (!RatioLog(numerator, denominator, logs[round])) {
                        return false;
                    }
                }
                if (!ComputeLogStats(logs.data(), logs.size(), stats)) {
                    return false;
                }
                PrintStatsLine(
                    "DIAGNOSTIC", period_names[period_metric], mode, stats,
                    info.BlockBytes, info.ScopeValue, info.Root);
            }

            for (uint32_t round = 0u; round < kRoundCount; ++round) {
                const SuperblockObservation* observation = FindObservation(
                    observations, by_key, cell_index, round, mode);
                if (!observation) return false;
                logs[round] = observation->LogFactorialInteraction;
            }
            if (!ComputeLogStats(logs.data(), logs.size(), stats)) return false;
            PrintStatsLine(
                "DIAGNOSTIC", "label-by-physical-interaction", mode, stats,
                info.BlockBytes, info.ScopeValue, info.Root);

            const char* sequence_position_names[] = {
                "sequence-position-1-vs-4",
                "sequence-position-2-vs-3",
                "sequence-position-edge-vs-middle"
            };
            for (uint32_t metric = 0u; metric < 3u; ++metric)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round) {
                    const SuperblockObservation* observation = FindObservation(
                        observations, by_key, cell_index, round, mode);
                    if (!observation || !SequencePositionDiagnosticLog(
                            *observation, metric, logs[round])) return false;
                }
                if (!ComputeLogStats(logs.data(), logs.size(), stats)) {
                    return false;
                }
                PrintStatsLine(
                    "DIAGNOSTIC", sequence_position_names[metric], mode,
                    stats, info.BlockBytes, info.ScopeValue, info.Root);
            }
        }

        for (uint32_t estimand_value = 0u;
             estimand_value < 2u;
             ++estimand_value)
        {
            const Estimand estimand = static_cast<Estimand>(estimand_value);
            const char* paired_names[] = {
                "combined-mode", "exposed-minus-isolated",
                "first-minus-second-superblock"
            };
            for (uint32_t metric = 0u; metric < 3u; ++metric)
            {
                for (uint32_t round = 0u; round < kRoundCount; ++round)
                {
                    const SuperblockObservation* exposed = FindObservation(
                        observations, by_key, cell_index, round, Mode::Exposed);
                    const SuperblockObservation* isolated = FindObservation(
                        observations, by_key, cell_index, round, Mode::Isolated);
                    if (!exposed || !isolated) return false;
                    const double exposed_log = ObservationLog(*exposed, estimand);
                    const double isolated_log =
                        ObservationLog(*isolated, estimand);
                    if (metric == 0u) {
                        logs[round] = (exposed_log + isolated_log) / 2.0;
                    }
                    else if (metric == 1u) {
                        logs[round] = exposed_log - isolated_log;
                    }
                    else
                    {
                        const double first =
                            exposed->Plan.ModePosition == 0u ?
                                exposed_log : isolated_log;
                        const double second =
                            exposed->Plan.ModePosition == 1u ?
                                exposed_log : isolated_log;
                        logs[round] = first - second;
                    }
                }
                LogStats stats;
                if (!ComputeLogStats(logs.data(), logs.size(), stats)) {
                    return false;
                }
                std::string diagnostic = paired_names[metric];
                diagnostic += estimand == Estimand::Logical ?
                    "-logical" : "-physical";
                PrintStatsLine(
                    "DIAGNOSTIC", diagnostic.c_str(), Mode::Exposed, stats,
                    info.BlockBytes, info.ScopeValue, info.Root,
                    nullptr, "paired");
            }
        }

        const char* transition_names[] = {
            "predecessor-cross-vs-same", "transition-ab-vs-ba",
            "transition-aa-vs-bb"
        };
        LogStats stats;
        for (uint32_t metric = 0u; metric < 3u; ++metric)
        {
            for (uint32_t round = 0u; round < kRoundCount; ++round)
            {
                const SuperblockObservation* exposed = FindObservation(
                    observations, by_key, cell_index, round, Mode::Exposed);
                const uint64_t aa = exposed->TransitionSums[0];
                const uint64_t ab = exposed->TransitionSums[1];
                const uint64_t ba = exposed->TransitionSums[2];
                const uint64_t bb = exposed->TransitionSums[3];
                uint64_t numerator = 0u;
                uint64_t denominator = 0u;
                if (metric == 0u)
                {
                    if (!AddElapsed(ab, numerator) ||
                        !AddElapsed(ba, numerator) ||
                        !AddElapsed(aa, denominator) ||
                        !AddElapsed(bb, denominator))
                    {
                        return false;
                    }
                }
                else if (metric == 1u) {
                    numerator = ab;
                    denominator = ba;
                }
                else {
                    numerator = aa;
                    denominator = bb;
                }
                if (!RatioLog(numerator, denominator, logs[round])) {
                    return false;
                }
            }
            if (!ComputeLogStats(logs.data(), logs.size(), stats)) return false;
            PrintStatsLine(
                "DIAGNOSTIC", transition_names[metric], Mode::Exposed,
                stats, info.BlockBytes, info.ScopeValue, info.Root);
        }
    }
    return PrintPooledDiagnostics(observations, by_key);
}

std::string ConfigSha256(
    const CellArray& cells,
    const std::string& schedule_sha256)
{
    std::string bytes;
    bytes.reserve(8192u);
    static const char marker[] = "wirehair1-native-aa-batch8-config-v1";
    bytes.append(marker, sizeof(marker));
    AppendLe32(bytes, kK);
    AppendLe32(bytes, kLossPpm);
    AppendLe32(bytes, kReceiveOverheadCap);
    AppendLe32(bytes, kRoundCount);
    AppendLe32(bytes, kCellCount);
    AppendLe32(bytes, kBatchSize);
    AppendLe64(bytes, static_cast<uint64_t>(kSuperblockCount));
    AppendLe64(bytes, static_cast<uint64_t>(kPanelActionCount));
    AppendLe64(bytes, static_cast<uint64_t>(kMeasuredBatchCount));
    AppendLe64(bytes, static_cast<uint64_t>(kWashoutBatchCount));
    AppendLe64(bytes, static_cast<uint64_t>(kMeasuredInnerInvocationCount));
    AppendLe64(bytes, static_cast<uint64_t>(kWashoutInnerInvocationCount));
    AppendLe64(bytes, static_cast<uint64_t>(kPinnedWarmBatchCount));
    AppendLe64(bytes, static_cast<uint64_t>(kUnpinnedPreflightBatchCount));
    static const char statistical_contract[] =
        "qmax=0.005;t975_n96=1.985251004;"
        "t975_n576=1.964095;t975_n32=2.039513446;"
        "interval=[0.98,1.02];official_marginal_gates=28";
    bytes.append(statistical_contract, sizeof(statistical_contract));
    bytes.append(schedule_sha256);
    for (uint32_t cell = 0u; cell < kBaseCellCount; ++cell)
    {
        const CalibrationCell& value = *cells[cell];
        AppendLe32(bytes, value.Index);
        AppendLe32(bytes, value.WidthIndex);
        AppendLe32(bytes, value.RootIndex);
        AppendLe32(bytes, value.BlockBytes);
        AppendLe64(bytes, value.Root);
        AppendLe32(bytes, value.IidOverhead);
        bytes.append(value.SourceSha256);
        bytes.append(value.IidIdsSha256);
        bytes.append(value.RepairIdsSha256);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

std::string RawSha256(const std::vector<RawTiming>& raw)
{
    std::string bytes;
    bytes.reserve(raw.size() * 256u);
    static const char marker[] = "wirehair1-native-aa-batch8-raw-v1";
    bytes.append(marker, sizeof(marker) - 1u);
    AppendLe64(bytes, static_cast<uint64_t>(raw.size()));
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        const RawTiming& receipt = raw[i];
        if (receipt.Cell >= kCellCount) return std::string();
        const LogicalCellInfo info = GetLogicalCellInfo(receipt.Cell);
        AppendLe64(bytes, static_cast<uint64_t>(i));
        AppendLe64(bytes, receipt.PanelActionIndex);
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
        AppendLe32(bytes, static_cast<uint32_t>(receipt.LogicalAPhysical));
        AppendLe32(bytes, receipt.Rotation);
        AppendLe32(bytes, receipt.SequencePosition);
        AppendLe32(bytes, receipt.Sequence);
        AppendLe32(bytes, receipt.SequenceFamily);
        AppendLe32(bytes, receipt.Period);
        AppendLe32(bytes, static_cast<uint32_t>(receipt.Label));
        AppendLe32(bytes, static_cast<uint32_t>(receipt.Physical));
        AppendLe32(bytes, static_cast<uint32_t>(receipt.Predecessor));
        AppendLe32(bytes, kBatchSize);
        AppendLe64(bytes, receipt.Batch.ElapsedNanoseconds);
        for (size_t inner = 0u; inner < kBatchSize; ++inner) {
            AppendLe32(bytes, static_cast<uint32_t>(receipt.Batch.Results[inner]));
        }
        AppendLe32(bytes, receipt.Batch.VerificationMask);
        for (size_t inner = 0u; inner < kBatchSize; ++inner) {
            AppendLe32(bytes, receipt.Batch.DecodedOverheads[inner]);
        }
        for (size_t inner = 0u; inner < kBatchSize; ++inner) {
            AppendLe64(bytes, receipt.Batch.DirectSystematicPackets[inner]);
        }
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

void PrintRawReceipts(const std::vector<RawTiming>& raw)
{
    std::cout << "RAW_SCHEMA,columns="
              << "index:panel_action:c:width_index:width_bytes:scope_index:"
              << "scope:root_index:root:round:visit:u:direction:mode_position:mode:"
              << "orientation_flip:map_A:rotation:sequence_position:sequence:"
              << "family:period:label:physical:predecessor:batch_size:"
              << "elapsed_ns:results:verification_mask:overheads:"
              << "direct_systematic_packets\n";
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        const RawTiming& receipt = raw[i];
        const LogicalCellInfo info = GetLogicalCellInfo(receipt.Cell);
        std::cout << "RAW," << i
                  << ',' << receipt.PanelActionIndex
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
                  << ',' << receipt.ModePosition
                  << ',' << ModeName(receipt.ModeValue)
                  << ',' << receipt.OrientationFlip
                  << ',' << BankName(receipt.LogicalAPhysical)
                  << ',' << receipt.Rotation
                  << ',' << receipt.SequencePosition
                  << ',' << SequenceName(receipt.Sequence)
                  << ',' << receipt.SequenceFamily
                  << ',' << (receipt.Period + 1u)
                  << ',' << LabelName(receipt.Label)
                  << ',' << BankName(receipt.Physical)
                  << ',' << PredecessorName(receipt.Predecessor)
                  << ',' << kBatchSize
                  << ',' << receipt.Batch.ElapsedNanoseconds
                  << ',';
        for (size_t inner = 0u; inner < kBatchSize; ++inner) {
            if (inner != 0u) std::cout << ':';
            std::cout << static_cast<uint32_t>(receipt.Batch.Results[inner]);
        }
        std::cout << ",0x" << std::hex << receipt.Batch.VerificationMask
                  << std::dec << ',';
        for (size_t inner = 0u; inner < kBatchSize; ++inner) {
            if (inner != 0u) std::cout << ':';
            std::cout << receipt.Batch.DecodedOverheads[inner];
        }
        std::cout << ',';
        for (size_t inner = 0u; inner < kBatchSize; ++inner) {
            if (inner != 0u) std::cout << ':';
            std::cout << receipt.Batch.DirectSystematicPackets[inner];
        }
        std::cout << '\n';
    }
}

struct TimingSpacing
{
    size_t Count = 0u;
    size_t ZeroCount = 0u;
    size_t DistinctValueCount = 0u;
    size_t DistinctSpacingCount = 0u;
    uint64_t MinimumPositiveSpacing = 0u;
};

bool ComputeTimingSpacing(
    const std::vector<RawTiming>& raw,
    int width_index,
    int scope_index,
    TimingSpacing& spacing_out)
{
    TimingSpacing spacing;
    std::set<uint64_t> distinct_values;
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        if (raw[i].Cell >= kCellCount) return false;
        const LogicalCellInfo info = GetLogicalCellInfo(raw[i].Cell);
        if ((width_index >= 0 &&
                info.WidthIndex != static_cast<uint32_t>(width_index)) ||
            (scope_index >= 0 &&
                static_cast<uint32_t>(info.ScopeValue) !=
                    static_cast<uint32_t>(scope_index)))
        {
            continue;
        }
        ++spacing.Count;
        if (raw[i].Batch.ElapsedNanoseconds == 0u) ++spacing.ZeroCount;
        distinct_values.insert(raw[i].Batch.ElapsedNanoseconds);
    }
    std::set<uint64_t> distinct_spacings;
    uint64_t previous = 0u;
    bool have_previous = false;
    for (std::set<uint64_t>::const_iterator it = distinct_values.begin();
         it != distinct_values.end();
         ++it)
    {
        if (have_previous && *it > previous) {
            distinct_spacings.insert(*it - previous);
        }
        previous = *it;
        have_previous = true;
    }
    spacing.DistinctValueCount = distinct_values.size();
    spacing.DistinctSpacingCount = distinct_spacings.size();
    spacing.MinimumPositiveSpacing = distinct_spacings.empty() ? 0u :
        *distinct_spacings.begin();
    spacing_out = spacing;
    return spacing.Count != 0u && spacing.ZeroCount == 0u &&
        spacing.DistinctValueCount != 0u;
}

bool PrintClockDiagnostics(
    const ClockDefinition& clock,
    const std::vector<RawTiming>& raw)
{
    TimingSpacing global;
    if (!ComputeTimingSpacing(raw, -1, -1, global) ||
        global.Count != kMeasuredBatchCount)
    {
        return false;
    }
    std::cout << "CLOCK,steady_period_numerator=" << clock.SteadyNumerator
              << ",steady_period_denominator=" << clock.SteadyDenominator
              << ",steady_period_ceil_ns="
              << clock.SteadyCeilNanoseconds
              << ",clock_monotonic_resolution_ns="
              << clock.MonotonicResolutionNanoseconds
              << ",rho_ns=" << clock.RhoNanoseconds
              << ",measured_batch_count=" << global.Count
              << ",zero_measured_count=" << global.ZeroCount
              << ",distinct_elapsed_values=" << global.DistinctValueCount
              << ",distinct_positive_spacings="
              << global.DistinctSpacingCount
              << ",minimum_positive_spacing_ns="
              << global.MinimumPositiveSpacing << '\n';
    for (uint32_t width = 0u; width < 2u; ++width)
    {
        for (uint32_t scope = 0u; scope < 3u; ++scope)
        {
            TimingSpacing group;
            if (!ComputeTimingSpacing(
                    raw,
                    static_cast<int>(width),
                    static_cast<int>(scope),
                    group) ||
                group.Count != 3072u)
            {
                return false;
            }
            std::cout << "TICK_SPACING,block_bytes=" << kBlockBytes[width]
                      << ",scope="
                      << ScopeName(static_cast<Scope>(scope))
                      << ",measured_count=" << group.Count
                      << ",zero_measured_count=" << group.ZeroCount
                      << ",distinct_elapsed_values="
                      << group.DistinctValueCount
                      << ",distinct_positive_spacings="
                      << group.DistinctSpacingCount
                      << ",minimum_positive_spacing_ns="
                      << group.MinimumPositiveSpacing << '\n';
        }
    }
    return true;
}

bool TestScheduleBalance()
{
    FrozenSchedule schedule;
    if (!BuildFrozenSchedule(schedule)) return false;
    size_t next = 0u;
    std::array<std::array<bool, kRoundCount>, kCellCount> seen_u = {{}};
    uint32_t flip_counts[kCellCount][kModeCount][2] = {};
    uint32_t rotation_counts[kCellCount][kModeCount][kSequenceCount] = {};
    uint32_t position_counts[kCellCount][kModeCount][kModeCount] = {};
    uint32_t direction_order_counts[kCellCount][2][2] = {};
    for (uint32_t round = 0u; round < kRoundCount; ++round)
    {
        std::array<bool, kCellCount> seen_cells = {{ false }};
        for (uint32_t visit = 0u; visit < kCellCount; ++visit)
        {
            const uint32_t cell = CellForVisit(round, visit);
            if (cell >= kCellCount || seen_cells[cell]) return false;
            seen_cells[cell] = true;
            const uint32_t u = ScheduleU(cell, round);
            if (seen_u[cell][u]) return false;
            seen_u[cell][u] = true;
            const SuperblockPlan& first = schedule[next++];
            const SuperblockPlan& second = schedule[next++];
            if (first.Cell != cell || second.Cell != cell ||
                first.Round != round || second.Round != round ||
                first.U != u || second.U != u ||
                first.ModePosition != 0u || second.ModePosition != 1u ||
                first.ModeValue == second.ModeValue ||
                first.ModeValue != ((u & 1u) == 0u ?
                    Mode::Exposed : Mode::Isolated) ||
                first.DirectionValue != second.DirectionValue ||
                first.DirectionValue !=
                    ((((round >> 1u) & 1u) == 0u) ?
                        Direction::Forward : Direction::Reverse))
            {
                return false;
            }
            ++direction_order_counts[cell]
                [static_cast<uint32_t>(first.DirectionValue)][u & 1u];
            const SuperblockPlan* const plans[] = { &first, &second };
            for (size_t plan_index = 0u; plan_index < 2u; ++plan_index)
            {
                const SuperblockPlan& plan = *plans[plan_index];
                const uint32_t mode = static_cast<uint32_t>(plan.ModeValue);
                if (mode >= kModeCount || plan.OrientationFlip >= 2u ||
                    plan.Rotation >= kSequenceCount ||
                    plan.ModePosition >= kModeCount)
                {
                    return false;
                }
                ++flip_counts[cell][mode][plan.OrientationFlip];
                ++rotation_counts[cell][mode][plan.Rotation];
                ++position_counts[cell][mode][plan.ModePosition];

                for (uint32_t period = 0u;
                     period < kPeriodsPerSequence;
                     ++period)
                {
                    uint32_t factorial[4] = {};
                    for (uint32_t sequence = 0u;
                         sequence < kSequenceCount;
                         ++sequence)
                    {
                        const LogicalLabel label = kSequences[sequence][period];
                        const PhysicalBank logical_a =
                            LogicalAPhysicalForSequence(plan, sequence);
                        const PhysicalBank physical = PhysicalForLabel(
                            logical_a, label);
                        ++factorial[static_cast<uint32_t>(label) * 2u +
                            static_cast<uint32_t>(physical)];
                    }
                    for (size_t factor = 0u; factor < 4u; ++factor) {
                        if (factorial[factor] != 1u) return false;
                    }
                }
            }
        }
        if (std::find(seen_cells.begin(), seen_cells.end(), false) !=
            seen_cells.end())
        {
            return false;
        }
    }
    if (next != schedule.size()) return false;
    for (uint32_t cell = 0u; cell < kCellCount; ++cell)
    {
        if (std::find(seen_u[cell].begin(), seen_u[cell].end(), false) !=
            seen_u[cell].end())
        {
            return false;
        }
        for (uint32_t mode = 0u; mode < kModeCount; ++mode)
        {
            if (flip_counts[cell][mode][0] != 16u ||
                flip_counts[cell][mode][1] != 16u ||
                position_counts[cell][mode][0] != 16u ||
                position_counts[cell][mode][1] != 16u)
            {
                return false;
            }
            for (uint32_t rotation = 0u;
                 rotation < kSequenceCount;
                 ++rotation)
            {
                if (rotation_counts[cell][mode][rotation] != 8u) return false;
            }
        }
        for (uint32_t direction = 0u; direction < 2u; ++direction) {
            for (uint32_t order = 0u; order < 2u; ++order) {
                if (direction_order_counts[cell][direction][order] != 8u) {
                    return false;
                }
            }
        }
    }
    const std::string schedule_sha256 = ScheduleSha256(schedule);
    if (schedule_sha256 != kExpectedScheduleSha256) {
        std::cerr << "batch8 schedule digest mismatch: "
                  << schedule_sha256 << '\n';
        return false;
    }
    return kPanelActionCount == 36864u &&
        kPinnedWarmBatchCount == 54u &&
        kPinnedWarmInnerInvocationCount == 432u &&
        kExpectedAffinityVerificationCount == UINT64_C(73836);
}

bool TestSequenceSurface()
{
    std::array<uint32_t, 4> transitions = {{ 0u, 0u, 0u, 0u }};
    for (uint32_t period = 0u; period < kPeriodsPerSequence; ++period)
    {
        uint32_t a_count = 0u;
        for (uint32_t sequence = 0u; sequence < kSequenceCount; ++sequence)
        {
            if (kSequences[sequence][period] == LogicalLabel::A) ++a_count;
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
    for (uint32_t rotation = 0u; rotation < kSequenceCount; ++rotation)
    {
        std::array<bool, kSequenceCount> seen = {{ false }};
        for (uint32_t position = 0u; position < kSequenceCount; ++position) {
            seen[(rotation + position) % kSequenceCount] = true;
        }
        if (std::find(seen.begin(), seen.end(), false) != seen.end()) {
            return false;
        }
    }
    return kMeasuredBatchCount == 18432u &&
        kMeasuredInnerInvocationCount == 147456u &&
        kWashoutBatchCount == 18432u &&
        kWashoutInnerInvocationCount == 147456u &&
        kUnpinnedPreflightBatchCount == 54u &&
        kUnpinnedPreflightInnerInvocationCount == 432u;
}

bool TestCanonicalAccumulation()
{
    SuperblockObservation observation;
    observation.Plan = MakePlan(0u, 0u, 0u, 0u, Mode::Exposed);
    observation.WashoutBatches = kWashoutBatchesPerSuperblock;
    for (uint32_t sequence_position = 0u;
         sequence_position < kSequenceCount;
         ++sequence_position)
    {
        const uint32_t sequence = sequence_position;
        const PhysicalBank logical_a =
            LogicalAPhysicalForSequence(observation.Plan, sequence);
        for (uint32_t period = 0u; period < kPeriodsPerSequence; ++period)
        {
            const LogicalLabel label = kSequences[sequence][period];
            const PhysicalBank physical = PhysicalForLabel(
                logical_a, label);
            const uint32_t factor = static_cast<uint32_t>(label) * 2u +
                static_cast<uint32_t>(physical);
            static const uint64_t elapsed_by_factor[4] = {
                12000u, 10000u, 9000u, 8000u
            };
            const uint64_t elapsed = elapsed_by_factor[factor];
            if (!AccumulateMeasured(
                    observation,
                    sequence_position,
                    sequence,
                    period,
                    label,
                    physical,
                    elapsed))
            {
                return false;
            }
        }
    }
    if (!FinalizeObservation(1u, observation)) return false;
    return observation.LogicalSums[0] == 88000u &&
        observation.LogicalSums[1] == 68000u &&
        observation.PhysicalSums[0] == 84000u &&
        observation.PhysicalSums[1] == 72000u &&
        observation.LogicalQuantizationPass &&
        observation.PhysicalQuantizationPass &&
        std::fabs(observation.LogLogicalAa - std::log(22.0 / 17.0)) < 1e-12 &&
        std::fabs(observation.LogPhysicalXy - std::log(7.0 / 6.0)) < 1e-12 &&
        std::fabs(observation.LogFactorialInteraction -
            std::log(16.0 / 15.0)) < 1e-12;
}

bool TestStatisticsAndGate()
{
    std::array<double, 96> zeros = {{ 0.0 }};
    LogStats stats;
    if (!ComputeLogStats(zeros.data(), zeros.size(), stats) ||
        stats.GeometricMean != 1.0 || stats.Lower95 != 1.0 ||
        stats.Upper95 != 1.0 || !OfficialGate(stats))
    {
        return false;
    }
    std::array<double, 96> shifted;
    shifted.fill(std::log(1.001));
    if (!ComputeLogStats(shifted.data(), shifted.size(), stats) ||
        OfficialGate(stats))
    {
        return false;
    }
    std::array<double, 576> pooled = {{ 0.0 }};
    for (size_t i = 0u; i < pooled.size(); ++i) {
        pooled[i] = (i & 1u) == 0u ? -0.3 : 0.3;
    }
    if (!ComputeLogStats(pooled.data(), pooled.size(), stats) ||
        OfficialGate(stats) || StudentT975(31u) != 0.0 ||
        StudentT975(96u) != 1.985251004 ||
        StudentT975(576u) != 1.964095)
    {
        return false;
    }
    return !ComputeLogStats(zeros.data(), 95u, stats);
}

bool TestQuantizationGate()
{
    double budget = 0.0;
    if (!QuantizationBudget(1u, 3200u, 3200u, budget) ||
        std::fabs(budget - 0.005) > 1e-15)
    {
        return false;
    }
    SuperblockObservation retained_failure;
    retained_failure.LogicalSums[0] = 1000u;
    retained_failure.LogicalSums[1] = 1000u;
    retained_failure.PhysicalSums[0] = 1000u;
    retained_failure.PhysicalSums[1] = 1000u;
    retained_failure.FactorialSums.fill(250u);
    retained_failure.LogicalCounts.fill(8u);
    retained_failure.PhysicalCounts.fill(8u);
    retained_failure.FactorialCounts.fill(4u);
    retained_failure.WashoutBatches = kWashoutBatchesPerSuperblock;
    return !QuantizationBudget(1u, 3199u, 3199u, budget) &&
        !QuantizationBudget(1u, 0u, 3200u, budget) &&
        QuantizationBudget(100u, 1000000u, 1000000u, budget) &&
        FinalizeObservation(1u, retained_failure) &&
        !retained_failure.LogicalQuantizationPass &&
        !retained_failure.PhysicalQuantizationPass &&
        retained_failure.LogicalQuantizationBudget >
            static_cast<double>(kQuantizationBudgetLimit) &&
        retained_failure.PhysicalQuantizationBudget >
            static_cast<double>(kQuantizationBudgetLimit);
}

bool TestResultClassification()
{
    Batch8ArmResult result;
    for (size_t i = 0u; i < kBatchSize; ++i) {
        result.Results[i] = Wirehair_Success;
        result.DecodedOverheads[i] = UINT32_MAX;
        result.DirectSystematicPackets[i] = 0u;
    }
    result.VerificationMask = (1u << kBatchSize) - 1u;
    result.ElapsedNanoseconds = 0u;
    if (!ValidSemanticBatch8Result(Scope::Encoder, 7u, result, false)) {
        return false;
    }
    result.ElapsedNanoseconds = 1u;
    if (!ValidSemanticBatch8Result(Scope::Encoder, 7u, result, true)) {
        return false;
    }
    result.DecodedOverheads.fill(7u);
    if (!ValidSemanticBatch8Result(Scope::IidReceive, 7u, result, true) ||
        ValidSemanticBatch8Result(Scope::IidReceive, 6u, result, true))
    {
        return false;
    }
    result.DecodedOverheads.fill(0u);
    if (!ValidSemanticBatch8Result(
            Scope::TwoRepairReceive, 7u, result, true))
    {
        return false;
    }
    result.DecodedOverheads[7] = 1u;
    if (ValidSemanticBatch8Result(
            Scope::TwoRepairReceive, 7u, result, true))
    {
        return false;
    }
    result.DecodedOverheads[7] = 0u;
    result.VerificationMask &= ~(1u << 7u);
    if (ValidSemanticBatch8Result(
            Scope::TwoRepairReceive, 7u, result, true)) return false;
    result.VerificationMask |= 1u << 7u;
    result.Results[7] = Wirehair_Error;
    return !ValidSemanticBatch8Result(
        Scope::TwoRepairReceive, 7u, result, true);
}

bool TestOutputFailureClassification()
{
    FailingStreamBuffer failing_buffer;
    std::streambuf* const original_buffer =
        std::cout.rdbuf(&failing_buffer);
    std::cout.clear();
    std::cout << "mandatory-evidence";
    const bool accepted = FlushEvidenceOutput();
    std::cout.rdbuf(original_buffer);
    std::cout.clear();
    return !accepted;
}

Batch8ArmResult SyntheticValidBatch8(
    Scope scope,
    uint32_t iid_overhead,
    uint64_t elapsed)
{
    Batch8ArmResult result;
    result.ElapsedNanoseconds = elapsed;
    result.VerificationMask = (1u << kBatchSize) - 1u;
    for (size_t i = 0u; i < kBatchSize; ++i)
    {
        result.Results[i] = Wirehair_Success;
        result.DecodedOverheads[i] = scope == Scope::Encoder ? UINT32_MAX :
            (scope == Scope::IidReceive ? iid_overhead : 0u);
        result.DirectSystematicPackets[i] = 0u;
    }
    return result;
}

bool TestPostPanelReceiptPipeline()
{
    FrozenSchedule schedule;
    if (!BuildFrozenSchedule(schedule)) return false;
    std::vector<RawTiming> raw(kMeasuredBatchCount);
    size_t next = 0u;
    for (size_t plan_index = 0u; plan_index < schedule.size(); ++plan_index)
    {
        const SuperblockPlan& plan = schedule[plan_index];
        const LogicalCellInfo info = GetLogicalCellInfo(plan.Cell);
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
            for (uint32_t period = 0u;
                 period < kPeriodsPerSequence;
                 ++period)
            {
                if (next >= raw.size()) return false;
                RawTiming& receipt = raw[next++];
                receipt.PanelActionIndex = static_cast<uint64_t>(
                    plan_index * (kMeasuredBatchesPerSuperblock +
                        kWashoutBatchesPerSuperblock) +
                    sequence_position * kBatchesPerSequence +
                    (plan.ModeValue == Mode::Exposed ?
                        kPeriodsPerSequence + period : 2u * period + 1u));
                receipt.Cell = plan.Cell;
                receipt.Round = plan.Round;
                receipt.Visit = plan.Visit;
                receipt.U = plan.U;
                receipt.DirectionValue = plan.DirectionValue;
                receipt.ModeValue = plan.ModeValue;
                receipt.ModePosition = plan.ModePosition;
                receipt.OrientationFlip = plan.OrientationFlip;
                receipt.LogicalAPhysical = logical_a;
                receipt.Rotation = plan.Rotation;
                receipt.SequencePosition = sequence_position;
                receipt.Sequence = sequence;
                receipt.SequenceFamily = sequence < 2u ? 0u : 1u;
                receipt.Period = period;
                receipt.Label = kSequences[sequence][period];
                receipt.Physical = receipt.Label == LogicalLabel::A ?
                    logical_a : logical_b;
                receipt.Predecessor = PredecessorClass::Washout;
                if (plan.ModeValue == Mode::Exposed && period > 0u)
                {
                    const uint32_t previous = static_cast<uint32_t>(
                        kSequences[sequence][period - 1u]);
                    receipt.Predecessor = static_cast<PredecessorClass>(
                        2u + previous * 2u +
                        static_cast<uint32_t>(receipt.Label));
                }
                receipt.Batch = SyntheticValidBatch8(
                    info.ScopeValue,
                    kExpectedIidOverheads[info.BaseIndex],
                    100000u);
            }
        }
    }
    if (next != raw.size()) return false;

    std::array<size_t, 6> predecessor_counts = {{ 0u, 0u, 0u, 0u, 0u, 0u }};
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        const uint32_t predecessor =
            static_cast<uint32_t>(raw[i].Predecessor);
        if (predecessor >= predecessor_counts.size()) return false;
        ++predecessor_counts[predecessor];
    }
    if (predecessor_counts[0] != 0u ||
        predecessor_counts[1] != 11520u ||
        predecessor_counts[2] != 1728u ||
        predecessor_counts[3] != 1728u ||
        predecessor_counts[4] != 1728u ||
        predecessor_counts[5] != 1728u)
    {
        return false;
    }

    std::vector<SuperblockObservation> observations(kSuperblockCount);
    std::vector<size_t> by_key(
        kSuperblockCount, std::numeric_limits<size_t>::max());
    std::string diagnostic;
    if (!BuildObservations(
            schedule, raw, 1u, observations, by_key, diagnostic))
    {
        return false;
    }
    size_t mutation_index = raw.size();
    for (size_t i = 0u; i < raw.size(); ++i)
    {
        const bool logical_a = raw[i].Label == LogicalLabel::A;
        const bool physical_x = raw[i].Physical == PhysicalBank::X;
        if (logical_a != physical_x)
        {
            mutation_index = i;
            break;
        }
    }
    if (mutation_index == raw.size()) return false;
    std::vector<RawTiming> mutated(raw);
    ++mutated[mutation_index].Batch.ElapsedNanoseconds;
    std::vector<SuperblockObservation> mutated_observations(kSuperblockCount);
    std::vector<size_t> mutated_by_key(
        kSuperblockCount, std::numeric_limits<size_t>::max());
    if (!BuildObservations(
            schedule,
            mutated,
            1u,
            mutated_observations,
            mutated_by_key,
            diagnostic))
    {
        return false;
    }
    const RawTiming& changed = raw[mutation_index];
    const SuperblockObservation* before = FindObservation(
        observations,
        by_key,
        changed.Cell,
        changed.Round,
        changed.ModeValue);
    const SuperblockObservation* after = FindObservation(
        mutated_observations,
        mutated_by_key,
        changed.Cell,
        changed.Round,
        changed.ModeValue);
    const long double increased = 800001.0L;
    const long double unchanged = 800000.0L;
    const double positive_delta = static_cast<double>(
        std::log(increased / unchanged));
    const double expected_logical = changed.Label == LogicalLabel::A ?
        positive_delta : -positive_delta;
    const double expected_physical = changed.Physical == PhysicalBank::X ?
        positive_delta : -positive_delta;
    if (!before || !after || before->LogLogicalAa != 0.0 ||
        before->LogPhysicalXy != 0.0 ||
        expected_logical * expected_physical >= 0.0 ||
        std::fabs(after->LogLogicalAa - expected_logical) > 1e-15 ||
        std::fabs(after->LogPhysicalXy - expected_physical) > 1e-15)
    {
        return false;
    }
    std::array<OfficialResult, kOfficialMarginalGateCount> official;
    bool all_pass = false;
    if (!BuildOfficialResults(
            observations, by_key, official, all_pass) || !all_pass)
    {
        return false;
    }
    for (size_t i = 0u; i < official.size(); ++i)
    {
        if (!official[i].Pass || official[i].Stats.GeometricMean != 1.0 ||
            official[i].Stats.Lower95 != 1.0 ||
            official[i].Stats.Upper95 != 1.0)
        {
            return false;
        }
    }
    TimingSpacing global_spacing;
    TimingSpacing group_spacing;
    if (!ComputeTimingSpacing(raw, -1, -1, global_spacing) ||
        global_spacing.Count != kMeasuredBatchCount ||
        global_spacing.ZeroCount != 0u ||
        global_spacing.DistinctValueCount != 1u ||
        global_spacing.MinimumPositiveSpacing != 0u ||
        !ComputeTimingSpacing(raw, 0, 0, group_spacing) ||
        group_spacing.Count != 3072u)
    {
        return false;
    }
    ClockDefinition synthetic_clock;
    synthetic_clock.SteadyNumerator = 1u;
    synthetic_clock.SteadyDenominator = UINT64_C(1000000000);
    synthetic_clock.SteadyCeilNanoseconds = 1u;
    synthetic_clock.MonotonicResolutionNanoseconds = 1u;
    synthetic_clock.RhoNanoseconds = 1u;
    std::ostringstream captured;
    std::streambuf* const previous_buffer = std::cout.rdbuf(captured.rdbuf());
    const bool diagnostic_output_ok = PrintMandatoryDiagnostics(
        observations, by_key);
    bool quantization_pass = false;
    const bool quantization_output_ok = PrintQuantizationDiagnostics(
        observations, by_key, quantization_pass);
    const bool clock_output_ok = PrintClockDiagnostics(synthetic_clock, raw);
    std::cout.rdbuf(previous_buffer);
    const std::string diagnostic_output = captured.str();
    if (!diagnostic_output_ok || !quantization_output_ok ||
        !quantization_pass || !clock_output_ok ||
        diagnostic_output.find("diagnostic=root-logical-a-b") ==
            std::string::npos ||
        diagnostic_output.find("diagnostic=physical-x-y") ==
            std::string::npos ||
        diagnostic_output.find("pooled-logical-a-b") ==
            std::string::npos ||
        diagnostic_output.find("pooled-physical-x-y") == std::string::npos ||
        diagnostic_output.find("pooled-label-by-physical-interaction") ==
            std::string::npos ||
        diagnostic_output.find("pooled-sequence-position-1-vs-4") ==
            std::string::npos ||
        diagnostic_output.find("pooled-period-1-vs-4") ==
            std::string::npos ||
        diagnostic_output.find(
            "DIAGNOSTIC,diagnostic=forward-vs-reverse,estimand=physical-x-y") ==
            std::string::npos ||
        diagnostic_output.find("pooled-combined-mode-logical") ==
            std::string::npos ||
        diagnostic_output.find("pooled-exposed-minus-isolated-logical") ==
            std::string::npos ||
        diagnostic_output.find(
            "pooled-first-minus-second-superblock-logical") ==
            std::string::npos ||
        diagnostic_output.find("pooled-combined-mode-physical") ==
            std::string::npos ||
        diagnostic_output.find("pooled-exposed-minus-isolated-physical") ==
            std::string::npos ||
        diagnostic_output.find(
            "pooled-first-minus-second-superblock-physical") ==
            std::string::npos ||
        diagnostic_output.find("pooled-predecessor-cross-vs-same") ==
            std::string::npos ||
        diagnostic_output.find("pooled-transition-ab-vs-ba") ==
            std::string::npos ||
        diagnostic_output.find("pooled-transition-aa-vs-bb") ==
            std::string::npos ||
        diagnostic_output.find(
            "DIAGNOSTIC,diagnostic=rounds-0-15-vs-16-31,estimand=physical-x-y") ==
            std::string::npos ||
        diagnostic_output.find("alias=rotation-q-plus-2") ==
            std::string::npos ||
        diagnostic_output.find("TICK_SPACING,block_bytes=64") ==
            std::string::npos ||
        diagnostic_output.find(
            "QUANTIZATION,estimand=logical-a-b,mode=exposed") ==
            std::string::npos ||
        diagnostic_output.find(
            "QUANTIZATION,estimand=physical-x-y,mode=exposed") ==
            std::string::npos)
    {
        return false;
    }
    const std::string original_hash = RawSha256(raw);
    if (original_hash.empty()) return false;
    std::vector<RawTiming> hash_mutation(raw);
    ++hash_mutation[0].Batch.ElapsedNanoseconds;
    if (RawSha256(hash_mutation) == original_hash) return false;
    hash_mutation = raw;
    hash_mutation[0].Batch.Results[7] = Wirehair_Error;
    if (RawSha256(hash_mutation) == original_hash) return false;
    hash_mutation = raw;
    hash_mutation[0].Batch.VerificationMask ^= UINT32_C(1) << 7u;
    if (RawSha256(hash_mutation) == original_hash) return false;
    hash_mutation = raw;
    --hash_mutation[0].Batch.DecodedOverheads[7];
    if (RawSha256(hash_mutation) == original_hash) return false;
    hash_mutation = raw;
    ++hash_mutation[0].Batch.DirectSystematicPackets[7];
    return RawSha256(hash_mutation) != original_hash;
}

bool TestFixturePreflights()
{
    CellArray cells;
    std::string diagnostic;
    if (!BuildCells(cells, diagnostic) ||
        !PreflightCells(cells, diagnostic))
    {
        return false;
    }
    for (uint32_t cell_index = 0u;
         cell_index < kBaseCellCount;
         ++cell_index)
    {
        const CalibrationCell& cell = *cells[cell_index];
        if (cell.Index != cell_index ||
            cell.WidthIndex != cell_index / 3u ||
            cell.RootIndex != cell_index % 3u ||
            cell.BlockBytes != kBlockBytes[cell.WidthIndex] ||
            cell.Root != kRoots[cell.RootIndex] ||
            cell.IidOverhead == UINT32_MAX ||
            cell.IidOverhead > kReceiveOverheadCap)
        {
            return false;
        }
        std::cout << "SELFTEST_BASE_CELL,base_cell=" << cell.Index
                  << ",block_bytes=" << cell.BlockBytes
                  << ",root=0x" << std::hex << cell.Root << std::dec
                  << ",iid_overhead=" << cell.IidOverhead
                  << ",source_sha256=" << cell.SourceSha256
                  << ",iid_ids_sha256=" << cell.IidIdsSha256
                  << ",repair_ids_sha256=" << cell.RepairIdsSha256 << '\n';
    }
    FrozenSchedule schedule;
    if (!BuildFrozenSchedule(schedule)) return false;
    const std::string config_sha256 = ConfigSha256(
        cells, ScheduleSha256(schedule));
    if (config_sha256 != kExpectedConfigSha256) {
        std::cerr << "batch8 config digest mismatch: "
                  << config_sha256 << '\n';
        return false;
    }
    std::cout << "SELFTEST_CONFIG,sha256=" << config_sha256 << '\n';
    return true;
}

bool RunSelfTest()
{
    if (!IsLowerHexSha256(WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256))
    {
        std::cerr << "native batch8 A/A selftest failed: embedded source digest\n";
        return false;
    }
    if (wirehair_init() != Wirehair_Success)
    {
        std::cerr << "native batch8 A/A selftest failed: wirehair_init\n";
        return false;
    }
    struct Test
    {
        const char* Name;
        bool (*Run)();
    };
    const Test tests[] = {
        { "schedule-balance-permutations", TestScheduleBalance },
        { "sequence-transition-surface", TestSequenceSurface },
        { "physical-logical-canonicalization", TestCanonicalAccumulation },
        { "statistics-official-gates", TestStatisticsAndGate },
        { "quantization-boundary", TestQuantizationGate },
        { "semantic-classification", TestResultClassification },
        { "output-failure-classification", TestOutputFailureClassification },
        { "post-panel-receipt-pipeline", TestPostPanelReceiptPipeline },
        { "fixture-preflights", TestFixturePreflights }
    };
    for (size_t i = 0u; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        if (!tests[i].Run())
        {
            std::cerr << "native batch8 A/A selftest failed: "
                      << tests[i].Name << '\n';
            return false;
        }
        std::cout << "SELFTEST,name=" << tests[i].Name
                  << ",status=PASS\n";
    }
    std::cout << "SELFTEST_SCHEDULE,sha256="
              << kExpectedScheduleSha256 << '\n';
    std::cout << "SELFTEST_SOURCE,sha256="
              << WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256 << '\n';
    std::cout << "WIREHAIR1_NATIVE_BATCH8_AA_SELFTEST=PASS\n";
    return true;
}

int InvalidCalibration(const std::string& diagnostic)
{
    std::cerr << "Wirehair1 native batch8 A/A calibration invalid: "
              << diagnostic << '\n';
    std::cout << "WIREHAIR1_NATIVE_BATCH8_AA_CALIBRATION=INVALID\n";
    return 1;
}

int RunCalibration(int target_cpu)
{
    try
    {
        std::string diagnostic;
        if (!IsLowerHexSha256(
                WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256)) {
            return InvalidCalibration("embedded source digest is invalid");
        }
        if (wirehair_init() != Wirehair_Success) {
            return InvalidCalibration("wirehair_init failed");
        }
        FrozenSchedule schedule;
        if (!BuildFrozenSchedule(schedule)) {
            return InvalidCalibration("frozen schedule construction failed");
        }
        const std::string schedule_sha256 = ScheduleSha256(schedule);
        if (schedule_sha256 != kExpectedScheduleSha256) {
            return InvalidCalibration("frozen schedule digest mismatch");
        }

        CellArray cells;
        if (!BuildCells(cells, diagnostic)) {
            return InvalidCalibration(diagnostic);
        }
        if (!PreflightCells(cells, diagnostic)) {
            return InvalidCalibration(diagnostic);
        }
        const std::string config_sha256 =
            ConfigSha256(cells, schedule_sha256);
        if (config_sha256 != kExpectedConfigSha256) {
            return InvalidCalibration("frozen configuration digest mismatch");
        }

        ClockDefinition clock;
        if (!ReadClockDefinition(clock, diagnostic)) {
            return InvalidCalibration(diagnostic);
        }
        std::vector<RawTiming> raw(kMeasuredBatchCount);
        std::vector<SuperblockObservation> observations(kSuperblockCount);
        std::vector<size_t> observation_by_key(
            kSuperblockCount, std::numeric_limits<size_t>::max());
        uint64_t affinity_verification_count = 0u;

        if (!PinCpuOnce(target_cpu, diagnostic)) {
            return InvalidCalibration(diagnostic);
        }
        if (!RunPinnedWarmups(
                cells,
                target_cpu,
                affinity_verification_count,
                diagnostic))
        {
            return InvalidCalibration(diagnostic);
        }
        if (!RunFrozenPanel(
                cells,
                schedule,
                target_cpu,
                affinity_verification_count,
                raw,
                diagnostic))
        {
            return InvalidCalibration(diagnostic);
        }
        if (affinity_verification_count !=
            kExpectedAffinityVerificationCount)
        {
            return InvalidCalibration(
                "total affinity-verification count invalid");
        }
        const std::string raw_sha256 = RawSha256(raw);
        const bool raw_hash_valid = !raw_sha256.empty();
        const bool observations_valid = BuildObservations(
                schedule,
                raw,
                clock.RhoNanoseconds,
                observations,
                observation_by_key,
                diagnostic);
        const std::string observation_diagnostic = diagnostic;

        std::array<OfficialResult, kOfficialMarginalGateCount> official = {{}};
        bool official_pass = false;
        const bool official_statistics_valid = observations_valid &&
            BuildOfficialResults(
                observations,
                observation_by_key,
                official,
                official_pass);

        std::cout << "CONFIG,protocol=wirehair1-native-aa-batch8-v1"
                  << ",source_commit=" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
                  << ",source_sha256="
                  << WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256
                  << ",cpu=" << target_cpu
                  << ",K=" << kK
                  << ",widths=64:1280"
                  << ",scopes=encoder:iid10:two-repair"
                  << ",roots=3,cells=" << kCellCount
                  << ",rounds=" << kRoundCount
                  << ",superblocks=" << kSuperblockCount
                  << ",batch_size=" << kBatchSize
                  << ",panel_actions=" << kPanelActionCount
                  << ",measured_batches=" << kMeasuredBatchCount
                  << ",measured_inner_invocations="
                  << kMeasuredInnerInvocationCount
                  << ",clock_free_washout_batches=" << kWashoutBatchCount
                  << ",clock_free_washout_inner_invocations="
                  << kWashoutInnerInvocationCount
                  << ",clock_free_pinned_warm_batches="
                  << kPinnedWarmBatchCount
                  << ",clock_free_pinned_warm_inner_invocations="
                  << kPinnedWarmInnerInvocationCount
                  << ",clock_free_fixture_preflight_batches="
                  << kUnpinnedPreflightBatchCount
                  << ",clock_free_fixture_preflight_inner_invocations="
                  << kUnpinnedPreflightInnerInvocationCount
                  << ",pin_verifications=1"
                  << ",affinity_verifications="
                  << affinity_verification_count
                  << ",quantization_limit="
                  << static_cast<double>(kQuantizationBudgetLimit)
                  << ",official_marginal_gates="
                  << kOfficialMarginalGateCount
                  << ",multiplicity_adjustment=none"
                  << ",official_interval_contract="
                     "each-contains-1-and-within-inclusive-0.98:1.02"
                  << ",schedule_sha256=" << schedule_sha256
                  << ",config_sha256=" << config_sha256
                  << ",raw_sha256="
                  << (raw_hash_valid ? raw_sha256 : "INVALID") << '\n';
        for (uint32_t cell_index = 0u;
             cell_index < kBaseCellCount;
             ++cell_index)
        {
            const CalibrationCell& cell = *cells[cell_index];
            std::cout << "BASE_CELL,base_cell=" << cell.Index
                      << ",width_index=" << cell.WidthIndex
                      << ",block_bytes=" << cell.BlockBytes
                      << ",root_index=" << cell.RootIndex
                      << ",root=0x" << std::hex << cell.Root << std::dec
                      << ",iid_overhead=" << cell.IidOverhead
                      << ",source_sha256=" << cell.SourceSha256
                      << ",iid_ids_sha256=" << cell.IidIdsSha256
                      << ",repair_ids_sha256=" << cell.RepairIdsSha256
                      << ",banks=X:Y:W-independent\n";
        }
        const bool clock_diagnostics_valid =
            PrintClockDiagnostics(clock, raw);
        PrintRawReceipts(raw);
        if (official_statistics_valid)
        {
            for (size_t i = 0u; i < official.size(); ++i)
            {
                const OfficialResult& result = official[i];
                if (result.Pooled)
                {
                    PrintStatsLine(
                        "OFFICIAL", result.EstimandValue == Estimand::Logical ?
                            "pooled-logical-a-b" : "pooled-physical-x-y",
                        result.ModeValue,
                        result.Stats, 0u, Scope::Encoder, 0u,
                        result.Pass ? "PASS" : "INVALID");
                }
                else
                {
                    PrintStatsLine(
                        "OFFICIAL",
                        result.EstimandValue == Estimand::Logical ?
                            "width-scope-logical-a-b" :
                            "width-scope-physical-x-y",
                        result.ModeValue, result.Stats,
                        kBlockBytes[result.WidthIndex], result.ScopeValue, 0u,
                        result.Pass ? "PASS" : "INVALID");
                }
            }
        }
        bool quantization_pass = false;
        const bool quantization_diagnostics_valid = observations_valid &&
            PrintQuantizationDiagnostics(
                observations, observation_by_key, quantization_pass);
        const bool mandatory_diagnostics_valid = observations_valid &&
            PrintMandatoryDiagnostics(observations, observation_by_key);
        std::cout << "INTERPRETATION,label_driver="
                     "follows-logical-A-B-after-X-Y-crossover"
                  << ",preparation_address_allocator=follows-physical-X-Y"
                  << ",carryover_cache=only-exposed-or-predecessor-dependent"
                  << ",order_drift_frequency="
                     "period-sequence-position-direction-or-half-without-predecessor"
                  << ",quantization=shortest-cell-budget-concentration"
                  << ",contradictory=unresolved-cause-non-gating\n";
        const bool calibration_valid = raw_hash_valid &&
            observations_valid && official_statistics_valid &&
            clock_diagnostics_valid && quantization_diagnostics_valid &&
            quantization_pass && mandatory_diagnostics_valid && official_pass;
        std::cout << "POST_PANEL_ANALYSIS,raw_hash="
                  << (raw_hash_valid ? "PASS" : "INVALID")
                  << ",observations="
                  << (observations_valid ? "PASS" : "INVALID")
                  << ",official_statistics="
                  << (official_statistics_valid ? "PASS" : "INVALID")
                  << ",clock_diagnostics="
                  << (clock_diagnostics_valid ? "PASS" : "INVALID")
                  << ",quantization_diagnostics="
                  << (quantization_diagnostics_valid ? "PASS" : "INVALID")
                  << ",quantization_gate="
                  << (quantization_pass ? "PASS" : "INVALID")
                  << ",mandatory_diagnostics="
                  << (mandatory_diagnostics_valid ? "PASS" : "INVALID")
                  << ",official_gates="
                  << (official_pass ? "PASS" : "INVALID");
        if (!observations_valid && !observation_diagnostic.empty()) {
            std::cout << ",observation_failure=" << observation_diagnostic;
        }
        std::cout << '\n';
        if (!FlushEvidenceOutput())
        {
            std::cerr << "Wirehair1 native batch8 A/A calibration invalid: "
                      << "mandatory evidence publication failed\n";
            return 1;
        }
        std::cout << "EVIDENCE_PUBLICATION=PASS\n";
        std::cout << "WIREHAIR1_NATIVE_BATCH8_AA_CALIBRATION="
                  << (calibration_valid ? "VALID" : "INVALID") << '\n';
        if (!FlushEvidenceOutput())
        {
            std::cerr << "Wirehair1 native batch8 A/A calibration invalid: "
                      << "terminal publication failed\n";
            return 1;
        }
        return calibration_valid ? 0 : 1;
    }
    catch (const std::bad_alloc&) {
        return InvalidCalibration("allocation failure");
    }
    catch (const std::length_error&) {
        return InvalidCalibration("length failure");
    }
    catch (const std::exception& exception) {
        return InvalidCalibration(
            std::string("unexpected exception: ") + exception.what());
    }
    catch (...) {
        return InvalidCalibration("unexpected non-standard exception");
    }
}

void PrintUsage(const char* program)
{
    std::cerr << "usage: " << program
              << " --selftest | --calibrate --cpu N\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0) {
        return RunSelfTest() ? 0 : 1;
    }
    if (argc == 4 && std::strcmp(argv[1], "--calibrate") == 0 &&
        std::strcmp(argv[2], "--cpu") == 0)
    {
        int target_cpu = -1;
        if (!ParseCpu(argv[3], target_cpu))
        {
            PrintUsage(argv[0]);
            return 1;
        }
        return RunCalibration(target_cpu);
    }
    PrintUsage(argv[0]);
    return 1;
}
