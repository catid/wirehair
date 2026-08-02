#ifndef WIREHAIR_BENCH_WH2_FROZEN_TRACE_H
#define WIREHAIR_BENCH_WH2_FROZEN_TRACE_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace wirehair {
namespace wh2_benchmark {

// These are the four packet-loss schedules frozen by the pure-GF(256)
// benchmark contract.  This benchmark-only API deliberately excludes the
// codec's other exploratory schedule kinds.
enum class FrozenSchedule
{
    Iid,
    Burst,
    Adversarial,
    RepairOnly
};

const char* FrozenScheduleName(FrozenSchedule schedule);

struct FrozenRecoveryCell
{
    std::size_t ordinal;
    std::string phase;
    std::string band;
    uint32_t K;
    uint32_t block_bytes;
    uint32_t loss_ppm;
    FrozenSchedule schedule;
    uint32_t stratum_index;
    uint32_t trial;
    uint32_t base_seed_attempt;
    uint64_t loss_seed;
    uint32_t overhead_cap;

    FrozenRecoveryCell();
};

// Returns the exact 360 arm-free development cells in checker order.
std::vector<FrozenRecoveryCell> EnumerateDevelopmentRecoveryCells();

// The JSON has exactly the recovery cell-key fields, sorted and compact in the
// same way as Python's json.dumps(sort_keys=True, separators=(",", ":")).
// Enumerator-only metadata such as ordinal and stratum_index is excluded.
std::string CanonicalRecoveryCellJson(const FrozenRecoveryCell& cell);
std::string RecoveryCellSha256(const FrozenRecoveryCell& cell);
std::string DevelopmentRecoveryDomainSha256();

// Small SHA-256 entry points used by the later native emitter.  A null pointer
// is accepted only for an empty byte string.  Invalid non-empty null input
// returns an empty string.
std::string Sha256Hex(const void* data, std::size_t bytes);
std::string Sha256Hex(const std::string& data);

// Hash packet IDs as consecutive explicit little-endian uint32 words.  This
// intentionally does not hash native object representation.
std::string PacketIdsSha256(const std::vector<uint32_t>& packet_ids);

uint64_t FrozenTraceSeed(
    uint32_t K,
    uint32_t block_bytes,
    uint64_t loss_seed);

// Compute 256 * delivered_count + 65536 without unsigned overflow.
bool FrozenCandidateLimit(
    uint64_t delivered_count,
    uint64_t& candidate_limit);

enum class FrozenTraceStatus
{
    Complete,
    CandidateLimitReached,
    InvalidInput
};

struct FrozenPacketTrace
{
    FrozenTraceStatus status;
    uint32_t K;
    uint32_t block_bytes;
    uint32_t loss_ppm;
    FrozenSchedule schedule;
    uint64_t loss_seed;
    uint64_t trace_seed;
    uint64_t attempted_candidates;
    uint64_t candidate_limit;
    std::vector<uint32_t> delivered_ids;
    std::string trace_sha256;

    FrozenPacketTrace();
};

// Generate exactly K+delivered_overhead packet IDs or stop at the derived
// candidate cap.  The explicit overhead is limited to 65536 so malformed
// benchmark input cannot request an effectively unbounded allocation.
// Contract K values are 2..64000.  loss_ppm=1000000 is accepted so the
// cap-failure path can be tested and reported deterministically.
FrozenTraceStatus GenerateFrozenPacketTrace(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t loss_ppm,
    uint64_t loss_seed,
    FrozenSchedule schedule,
    uint32_t delivered_overhead,
    FrozenPacketTrace& trace);

// Recovery-compatible entry point.  This remains exactly K+4 and is kept as
// a separate overload so existing recovery domains, receipts, and hashes do
// not change when timing requests a longer arm-free prefix.
FrozenTraceStatus GenerateFrozenPacketTrace(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t loss_ppm,
    uint64_t loss_seed,
    FrozenSchedule schedule,
    FrozenPacketTrace& trace);

FrozenTraceStatus GenerateFrozenPacketTrace(
    const FrozenRecoveryCell& cell,
    FrozenPacketTrace& trace);

// Copy the first K+overhead IDs from one already generated K+4 trace.  This is
// the only supported way to form the nested K+0..K+4 replay prefixes.
bool CopyNestedPrefix(
    const FrozenPacketTrace& trace,
    uint32_t overhead,
    std::vector<uint32_t>& prefix);

} // namespace wh2_benchmark
} // namespace wirehair

#endif // WIREHAIR_BENCH_WH2_FROZEN_TRACE_H
