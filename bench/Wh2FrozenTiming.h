#ifndef WIREHAIR_BENCH_WH2_FROZEN_TIMING_H
#define WIREHAIR_BENCH_WH2_FROZEN_TIMING_H

#include "Wh2FrozenTrace.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace wirehair {
namespace wh2_benchmark {

struct FrozenTimingCell
{
    std::size_t ordinal;
    std::string phase;
    std::string band;
    uint32_t K;
    uint32_t block_bytes;
    uint32_t loss_ppm;
    FrozenSchedule schedule;
    uint32_t replicate;
    uint32_t base_seed_attempt;
    uint64_t loss_seed;
    uint32_t fixed_received_overhead;
    std::string interleave_policy;
    uint32_t invocations_per_slot;

    FrozenTimingCell();
};

// The exact 2 widths x 96 replicates x 12 K values, in round, width, K,
// then eight-lane checker order.
std::vector<FrozenTimingCell> EnumerateDevelopmentTimingCells();

std::string CanonicalTimingCellJson(const FrozenTimingCell& cell);
std::string TimingCellSha256(const FrozenTimingCell& cell);
std::string DevelopmentTimingDomainSha256();

// Deterministic n for the self-counterbalanced measured panel.  Its primary
// four slots receive ceil(n/2) calls each and the opposite four slots receive
// floor(n/2) calls each.  Zero K is invalid and returns zero; all positive K
// values use max(2, ceil(65536 / K)) without overflowing at UINT32_MAX.
uint32_t DevelopmentTimingInvocationsPerSlot(uint32_t K);

// Frozen v4 execution geometry for the one-candidate development timing
// screen: 264 homogeneous (K, width, panel) cohorts, twelve independent
// eight-job rounds, and a barrier after every wave.  Worker slots index the
// sorted eight-CPU roster.
uint32_t DevelopmentTimingWorkerCount();
std::size_t DevelopmentTimingCohortCount();
bool DevelopmentTimingWorkerSlot(
    std::size_t cohort_index,
    uint32_t replicate,
    uint32_t& worker_slot);

// Derive the paired development loss seed and fixed production construction
// attempt for one of the 96 frozen replicates.  Loss roots continue to
// cycle, while timing always uses public construction attempt zero.
bool DevelopmentTimingSeed(
    uint32_t replicate,
    uint32_t& base_seed_attempt,
    uint64_t& loss_seed);

struct FrozenTimingTraceReceipt
{
    std::size_t ordinal;
    uint32_t K;
    std::string cell_sha256;
    std::string trace_sha256;
    uint64_t attempted_candidates;
    uint64_t candidate_limit;

    FrozenTimingTraceReceipt();
};

// Validate the cell, generate its one unsalted-IID K+4 packet prefix, and
// produce the arm-free fields needed by the canonical trace manifest.
FrozenTraceStatus GenerateDevelopmentTimingTrace(
    const FrozenTimingCell& cell,
    FrozenPacketTrace& trace,
    FrozenTimingTraceReceipt& receipt);

// Canonical checker row: cell_sha256, ordinal, trace_sha256.  Diagnostic
// candidate counts remain in FrozenTimingTraceReceipt but outside this row.
std::string CanonicalTimingTraceManifestRow(
    const FrozenTimingTraceReceipt& receipt);

enum class FrozenPanelKind
{
    AA,
    AB
};

enum class FrozenTimingScope
{
    DecoderSolve,
    ReceiveToSuccess,
    EncoderInitPlusFirstKSymbols
};

const char* FrozenPanelKindName(FrozenPanelKind kind);
const char* FrozenTimingScopeName(FrozenTimingScope scope);

struct FrozenTimingPanel
{
    std::size_t ordinal;
    FrozenPanelKind panel_kind;
    FrozenTimingScope scope;
    std::string left_arm;
    std::string right_arm;

    FrozenTimingPanel();
};

// Four control A/A, three candidate A/A, then four candidate/control A/B
// panels, matching timing_panels() for one candidate exactly.
std::vector<FrozenTimingPanel> EnumerateOneCandidateTimingPanels(
    const std::string& candidate_arm);

std::string CanonicalTimingPanelJson(const FrozenTimingPanel& panel);
std::string TimingPanelSha256(const FrozenTimingPanel& panel);

enum class FrozenTimingOrder
{
    ABBA,
    BAAB,
    Invalid
};

const char* FrozenTimingOrderName(FrozenTimingOrder order);

// The primary block is ABBA iff replicate parity equals the low bit of the
// canonical panel SHA.  The panel's second block always uses the opposite.
FrozenTimingOrder TimingPanelOrder(
    const FrozenTimingPanel& panel,
    uint32_t replicate);

} // namespace wh2_benchmark
} // namespace wirehair

#endif // WIREHAIR_BENCH_WH2_FROZEN_TIMING_H
