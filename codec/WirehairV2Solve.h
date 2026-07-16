#pragma once

#include "WirehairV2Precode.h"
#include "WirehairV2Seeds.h"

#include <wirehair/wirehair.h>

#include <stddef.h>
#include <stdint.h>
#include <vector>

namespace wirehair_v2 {

// Message packets use a distinct final contract number so they cannot be
// confused with the separately versioned experimental recovery-row helper.
static const uint32_t kPacketRowContractVersion = 4u;
static const uint32_t kPrecodeContractVersion = 2u; // existing/default alias
static const uint32_t kMixedPrecodeContractVersion = 3u;

inline uint32_t PrecodeContractVersion(CompletionField field)
{
    if (field == CompletionField::GF256) return kPrecodeContractVersion;
    if (field == CompletionField::MixedGF256GF16)
        return kMixedPrecodeContractVersion;
    return 0u;
}
static const uint32_t kCertifiedPacketMixCount = 3u;
static const uint32_t kMaxPacketSeedAttempts = 256u;
static const uint32_t kMaxInactiveColumns = 4096u;
static const uint32_t kMinPacketPrecodeCount = 2u;
static const uint32_t kMaxPacketPrecodeCount = 65521u;
static const uint32_t kBinaryQuotientMinBlockBytes = 2048u;

struct PacketRowConfig
{
    uint32_t PeelSeed = 0;
    uint32_t MixCount = kCertifiedPacketMixCount;
};

/**
    Validated process-local invariants for one packet-row domain.

    Prime values are deliberately private and are derived only by Initialize;
    matching counts/config therefore prove cached iterators are equation-identical
    to the generic helpers.  This state is never serialized into a V2 profile.
*/
class PacketRowRuntime
{
public:
    bool Initialize(
        uint32_t source_count,
        uint32_t precode_count,
        uint32_t mix_count);
    bool IsValidFor(
        uint32_t source_count,
        uint32_t precode_count,
        uint32_t mix_count) const;

    uint16_t SourcePrime() const { return SourcePrimeValue; }
    uint16_t PrecodePrime() const { return PrecodePrimeValue; }

private:
    uint32_t SourceCount = 0u;
    uint32_t PrecodeCount = 0u;
    uint32_t MixCount = 0u;
    uint16_t SourcePrimeValue = 0u;
    uint16_t PrecodePrimeValue = 0u;
};

struct SolvePacket
{
    uint32_t BlockId = 0;
    const uint8_t* Data = nullptr;
};

struct PrecodeSolveStats
{
    uint32_t PacketRows = 0;
    uint32_t PeeledColumns = 0;
    uint32_t InactivatedColumns = 0;
    uint32_t ResidualRows = 0;
    uint32_t ResidualRank = 0;
    uint32_t BinaryResidualRank = 0;
    uint64_t BinaryRowReferences = 0;
    uint64_t BinaryRowStorageBytes = 0;
    uint64_t BinaryAdjacencyStorageBytes = 0;
    uint32_t BinaryRowStorageAllocations = 0;
    uint32_t BinaryAdjacencyStorageAllocations = 0;
    uint64_t BlockXors = 0;
    uint64_t BlockMulAdds = 0;
    uint64_t BuildNanoseconds = 0;
    uint64_t PeelNanoseconds = 0;
    uint64_t ProjectNanoseconds = 0;
    uint64_t ResidualNanoseconds = 0;
    uint64_t BackSubNanoseconds = 0;
    uint32_t PacketSeedAttempt = 0;
};

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/**
    Multiply packet ids by an odd 32-bit constant before initializing the row
    PRNG.  Odd multiplication is a permutation of the complete id domain, so
    experiments can decorrelate consecutive packet streams without collisions.
    Returns false without changing the active multiplier for zero/even input.
*/
bool SetPacketRowSeedMultiplierForTesting(uint32_t multiplier);

/** Apply a bijective 32-bit avalanche permutation after id multiplication. */
void SetPacketRowSeedAvalancheForTesting(bool enabled);

/**
    XOR the peel seed for odd packet ids in the calling thread.  This is a
    test-only control for measuring whether interleaving two deterministic
    packet-graph phases reduces block-count resonances.  Zero restores the
    production equation mapping.
*/
void SetOddPacketPeelSeedXorForTesting(uint32_t seed_xor);

/**
    Enable an independent dense-expansion oracle for mixed completion
    coefficient projection.  Balanced enable/disable calls are nestable.  Test
    code uses this to compare the optimized residue-bucket projection against
    the original L-by-R expansion exactly.
*/
void SetMixedProjectionOracleForTesting(bool enabled);

/** Reset/read the number of successful optimized-versus-dense comparisons. */
void ResetMixedProjectionOracleComparisonsForTesting();
uint64_t MixedProjectionOracleComparisonsForTesting();

/**
    Enable an exact comparison between low-degree-XOR and row-scan binary
    peeling.  The oracle also rejects duplicate columns within any equation.
*/
void SetBinaryPeelOracleForTesting(bool enabled);

/** Reset/read the number of successful optimized-versus-scan comparisons. */
void ResetBinaryPeelOracleComparisonsForTesting();
uint64_t BinaryPeelOracleComparisonsForTesting();

/**
    Compare packed-GF(2) residual insertion against the byte GF(256) oracle at
    word boundaries, including poisoned tail bits and inconsistent RHS rows.
*/
bool CheckPackedBinaryResidualOracleForTesting();
#endif

/**
    Algebraic checkpoint for appending packet equations after a rank-deficient
    solve.  This is an internal movable value: callers should use
    ResumePrecodeSystem() rather than mutate its fields.
*/
struct PrecodeSolveResumeState
{
    uint32_t SourceCount = 0u;
    uint32_t PrecodeCount = 0u;
    uint32_t ColumnCount = 0u;
    uint32_t BlockBytes = 0u;
    uint32_t InactiveCount = 0u;
    uint32_t ProjectionWords = 0u;
    uint32_t Rank = 0u;
    PacketRowConfig Config = {};
    PacketRowRuntime Runtime = {};
    PrecodeSolveStats Stats = {};
    std::vector<uint32_t> InactiveIndex;
    std::vector<uint32_t> InactiveColumns;
    std::vector<uint64_t> Projection;
    std::vector<uint8_t> Values;
    std::vector<uint8_t> PivotCoefficients;
    std::vector<uint8_t> PivotRhs;
    std::vector<uint8_t> HavePivot;
    std::vector<uint8_t> CoefficientScratch;
    std::vector<uint8_t> RhsScratch;
    bool Active = false;

    void Clear();
    void Swap(PrecodeSolveResumeState& other) noexcept;
    size_t PersistentBytes() const;
};

/**
    Return whether the exact version-4 packet-row iterator supports this
    domain and mix count.

    This is deliberately narrower than BuildPrecodeSystem()'s structural
    domain.  Packet rows require at least two precode columns for parameter
    generation and use NextPrime16(), whose largest supported input is 65521.
    A structurally valid system outside this range may still be inspected or
    used by structure-only tooling, but it cannot be evaluated or solved as a
    version-4 packet profile.
*/
bool IsPacketRowDomainValid(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t mix_count);

/**
    Generate the version-4 packet equation for a public block id.

    The source prefix is contractually bound to production Wirehair's integer
    GeneratePeelRowWeight()/PeelRowIterator rule.  Exactly MixCount distinct
    precode columns are added with RowMixIterator.  Addressing by the public
    block id keeps systematic equations [0,K) disjoint from repair equations.
    Unsupported packet domains return an empty row without invoking either
    iterator.
*/
std::vector<uint32_t> GeneratePacketMatrixRow(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config);

/** Internal row generator using a validated process-local prime cache. */
std::vector<uint32_t> GeneratePacketMatrixRowWithRuntime(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime);

uint32_t PacketPeelSeedFromProfile(
    const SeedProfile& profile,
    uint64_t salt);

PacketRowConfig PacketConfigForAttempt(
    const PacketRowConfig& base,
    uint32_t attempt);

PrecodeParams PrecodeParamsForAttempt(
    const PrecodeParams& base,
    uint32_t attempt);

/**
    Evaluate one packet row over all intermediate blocks.

    `block_out[0, block_bytes)` must not overlap any byte in the complete
    `(K + P) * block_bytes` intermediate-block array.  Overlap is rejected
    before writing either `block_out` or `block_ops_out`.
*/
bool EvaluatePacketBlock(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/**
    Internal fast path for an already validated immutable system/config.

    The same non-overlap and failure no-write contract as
    EvaluatePacketBlock() applies.
*/
bool EvaluatePacketBlockForValidatedSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/** Validated-system evaluator using process-local cached packet primes. */
bool EvaluatePacketBlockForValidatedSystemWithRuntime(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/**
    Solve the complete V2 system over its configured completion field.

    Binary staircase/dense constraints and packet equations are peeled first.
    Unused binary rows are projected onto the inactivated columns.  The
    generic path expands them into a GF(256) residual before inserting its
    Cauchy heavy equations; the mixed path keeps that residual packed in
    GF(2), then solves only its remaining quotient over GF(2^16).  On success
    `intermediate_blocks_out` contains all
    K+S+D2+H block values.  The exact residual solve is bounded to
    kMaxInactiveColumns to contain adversarial memory use.  NeedMore means the
    supplied equations were rank deficient or exceeded that bound; additional
    independent packets can reduce the residual.  Output remains unchanged on
    every failure.  When resume_state is non-null, a rank-deficient residual
    within the cap atomically replaces it with an active affine/pivot
    checkpoint; cap failures leave it unchanged and require a cold retry.
    Mixed GF(256)/GF(2^16) solves do not publish resume checkpoints and leave
    any caller-provided resume state unchanged on NeedMore.
    `stats` is diagnostic rather than transactional output: completed algebraic
    outcomes publish their counters, while validation failures and allocation
    failures before a resumable checkpoint may leave the caller's prior value
    unchanged.  This lets stateful decoders preserve their last committed
    counters when a cold retry cannot be constructed.
*/
WirehairResult SolvePrecodeSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats = nullptr,
    PrecodeSolveResumeState* resume_state = nullptr);

/** Internal cold solve using cached packet primes; still validates system. */
WirehairResult SolvePrecodeSystemWithRuntime(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats = nullptr,
    PrecodeSolveResumeState* resume_state = nullptr);

/**
    Internal cold solve for an immutable system validated by construction.
    The caller must retain exclusive ownership of that trust boundary: unlike
    SolvePrecodeSystemWithRuntime(), this does not inspect the stored row graph.
*/
WirehairResult SolvePrecodeSystemForValidatedSystemWithRuntime(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats = nullptr,
    PrecodeSolveResumeState* resume_state = nullptr);

/**
    Append one packet equation to a rank-deficient solve checkpoint.

    With allow_insert=false this performs a non-mutating duplicate consistency
    check.  With allow_insert=true an independent row is committed and Success
    is returned as soon as the complete intermediate vector is reconstructed.
    Allocations finish before an inserting call changes the algebraic state, so
    OOM is retryable.  On OOM, stats receives the unchanged checkpoint counters
    when non-null.  Output remains unchanged on NeedMore and every failure.
*/
WirehairResult ResumePrecodeSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    uint32_t block_id,
    const uint8_t* block_data,
    uint32_t block_bytes,
    PrecodeSolveResumeState& resume_state,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats = nullptr,
    bool allow_insert = true);

/** Select the first deterministic packet seed whose K systematic rows rank. */
WirehairResult SelectSystematicPacketConfig(
    const PrecodeSystem& system,
    const PacketRowConfig& base_config,
    PacketRowConfig& selected_config,
    uint32_t* attempt_out = nullptr);

/** Select the first full-rank deterministic joint precode/packet seed. */
WirehairResult SelectSystematicConfiguration(
    const PrecodeParams& base_params,
    const PacketRowConfig& base_config,
    PrecodeSystem& selected_system,
    PacketRowConfig& selected_config,
    uint32_t* attempt_out = nullptr);

/** Expensive test/oracle validation of every supplied equation. */
bool VerifyPrecodeSolution(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes);

} // namespace wirehair_v2
