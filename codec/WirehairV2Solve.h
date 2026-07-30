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
static const uint32_t kMixedTwoAnchorPrecodeContractVersion = 4u;
static const uint32_t kMixedSmallBandD4PrecodeContractVersion = 5u;

inline uint32_t PrecodeContractVersion(
    CompletionField field,
    bool adaptive_dense_two_anchor = false,
    V2PrecodeArchitecture architecture =
        V2PrecodeArchitecture::LegacyD12)
{
    if (architecture == V2PrecodeArchitecture::SmallBandD4)
    {
        return field == CompletionField::MixedGF256GF16 &&
                !adaptive_dense_two_anchor ?
            kMixedSmallBandD4PrecodeContractVersion : 0u;
    }
    if (architecture != V2PrecodeArchitecture::LegacyD12) return 0u;
    if (adaptive_dense_two_anchor)
    {
        return field == CompletionField::MixedGF256GF16 ?
            kMixedTwoAnchorPrecodeContractVersion : 0u;
    }
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
// Certified production keeps packed residuals at and above this width even
// when a legacy resume checkpoint may fit.  Below it, a narrow resumable solve
// retains the byte basis to avoid a packed-to-byte conversion on NeedMore.
static const uint32_t kCertifiedPackedResumeMinBlockBytes = 2048u;
// The test-only forced byte reference retains the same historical quotient
// crossover so wide byte-vs-packed differentials compare identical algebra.
static const uint32_t kLegacyByteQuotientMinBlockBytes = 2048u;

/**
    Tiny-K mixed completion fast-path gate.

    At or below these bounds the mixed GF(256)/GF(2^16) completion quotient
    is solved by one direct dense scalar elimination over the inactivated
    system instead of the packed-projection/residue-bucket/factorization
    machinery.  The fast path evaluates the same equations in the same exact
    arithmetic, so every solve outcome (result code, solution bytes, and rank
    receipts) is identical; the bounds are performance tuning only.

    All three bounds are runtime constants with compile-time defaults chosen
    from paired idle-host nonzero-payload medians (forced-on versus
    forced-off, mixed-p244-d4 dispatch config plus mixed-p244-baseline as a
    second probe): the fast path engages when K <= max_source_blocks,
    block_bytes <= max_block_bytes, and K * block_bytes <= max_product_bytes.
    The measured forced-on/forced-off crossover sits at a flat K near 530 for
    2- and 4-byte blocks (a product-only bound would engage ~20%% regressions
    by K=1536 at 2-byte blocks) and at a product near 3072 bytes for 6- and
    8-byte blocks, so the flat K bound stays and the product bound carries
    the narrow-block extension.  Every engaged cell measured win-or-parity
    (worst boundary cell 1.003x on the dispatch config, 1.008x on the
    baseline probe, at the edge of noise); the nearest excluded cells
    measured with forced pairs (K=520 at 6-byte blocks, K=400 at 8-byte
    blocks, K=544 at 4-byte blocks, K=576 at 2-byte blocks) were 1-2%%
    regressions growing steadily with K.  A narrow residual win band
    (K=513..543 at 2-byte blocks, ~2%%) stays excluded: no single product
    bound can engage it without also engaging the 4-byte regression at the
    same K.  Bounds must be set from nonzero-payload timing: zero-payload
    probes let the fast path's zero-block skip flatter it by ~12%% at the
    2-byte seam.
*/
#ifndef WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_SOURCE_BLOCKS
#define WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_SOURCE_BLOCKS 512u
#endif
#ifndef WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_BLOCK_BYTES
#define WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_BLOCK_BYTES 8u
#endif
#ifndef WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_PRODUCT_BYTES
#define WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_PRODUCT_BYTES 3072u
#endif

/** Active tiny-K mixed fast-path bounds (compile-time defaults). */
uint32_t TinyMixedFastPathMaxSourceBlocks();
uint32_t TinyMixedFastPathMaxBlockBytes();
uint32_t TinyMixedFastPathMaxProductBytes();

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/**
    Override the tiny mixed fast-path dispatch on the calling thread:
    -1 disables it, +1 forces it wherever structurally possible, and 0
    restores the automatic K/block-bytes gate.  Because the fast path is
    equation-preserving this changes timing only, never solve results;
    tests use forced on/off pairs to prove exactly that.
*/
bool SetTinyMixedFastPathModeForTesting(int mode);
int TinyMixedFastPathModeForTesting();

/**
    Override the certified residual-basis representation on the calling
    thread: -1 forces the legacy byte matrix, +1 forces the packed GF(2)
    matrix, and 0 restores production's adaptive selection: packed without a
    checkpoint, at wide payloads, or when the decoder proves a legacy
    checkpoint cannot fit; byte at narrow widths where resumability is useful.
    This changes representation and timing only; the GF(256) quotient and
    public equations are identical.
*/
bool SetCertifiedPackedResidualModeForTesting(int mode);
int CertifiedPackedResidualModeForTesting();

/**
    Fail one of the cold certified packed-resume materialization allocations
    on this thread.  Zero fails the next allocation, positive values count
    down, and -1 disables injection.
*/
void SetCertifiedPackedResumeAllocationFailureCountdownForTesting(
    int64_t countdown);

/**
    Exhaustively verify the tiny fast path's scalar GF(2^16) helpers: the
    subfield-norm inverse against the power-chain reference for every field
    element, and the fused row muladd/scale kernels against the scalar
    reference kernels across widths and scale specials.
*/
bool CheckTinyMixedScalarHelpersForTesting();

/**
    Prove one operation's captured packet-row policy owns its hook snapshot
    and one multi-packet solve captures exactly one policy for the row batch.
*/
bool CheckPacketRowPolicySnapshotForTesting();
#endif

struct PacketRowConfig
{
    uint32_t PeelSeed = 0;
    uint32_t MixCount = kCertifiedPacketMixCount;
};

/**
    True when packet ids and packet peel seeds use the frozen production
    mapping.

    The receipted packet-degree PMF overlay is deliberately excluded.  Stable
    target identities use this to reject the unreceipted id multiplier,
    avalanche, and odd-id peel-seed experiments in test-hook builds.
*/
bool IsCanonicalStableTargetPacketRowState();

/**
    True only when no packet-degree PMF overlay (valid or malformed) is active.

    Stable benchmark targets may receipt a PMF separately.  Frozen equation
    contracts that do not name such an overlay combine this with
    IsCanonicalStableTargetPacketRowState().
*/
bool IsCanonicalPacketDegreeState();

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/**
    Compact semantic fingerprint of every ambient packet-row hook.

    The active id permutation, peel-seed phase, and normalized packet-degree
    CDF are included.  Equation-preserving solver/dispatch oracles are not.
*/
uint64_t ActivePacketRowEquationStateFingerprintForTesting();

/**
    Test the automatic joint-bucket decision through the solver's real
    dispatch-width normalization.  In particular, block_bytes == 0 means the
    calibrated counting width rather than a zero-byte production payload.
*/
bool UseAutomaticMixedJointResidueBucketsForSolveTesting(
    uint32_t block_count,
    uint32_t block_bytes,
    uint32_t coefficient_period);
#endif

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
    /**
        Deterministic non-XOR operation counters.

        These complete the block-work receipt so a calibrated cost model can
        predict solve latency without a stopwatch.  Every counter below is a
        pure function of the structural configuration and the delivered packet
        set, so repeated evaluations of one cell reproduce it exactly.

        BlockCopies      whole-payload-block logical copy passes: memcpy,
                         a distinguished first source fused into set-form
                         initialization, and payload-plane conversion passes.
        BlockZeroFills   whole-payload-block zero fills, including the bulk
                         zero-initialization of the value workspace (counted
                         in whole blocks, not bytes).
        PeelAdjacencyVisits  column-to-row adjacency entries visited while
                         resolving peeled and inactivated columns.
        PeelRowScanSteps  row-column steps walked by the degree-two scan and
                         the inactivation fallback cursor.
        PeelHeapOperations  degree-two push_heap calls and invalid-root
                         pop_heap calls.  A compaction rebuild uses make_heap,
                         so its work is reported by the two counters below
                         rather than as fictional logical insertions.
        PeelHeapCompactionRebuildColumnProbes  original binary-system columns
                         inspected while rebuilding the lazy heap.
        PeelHeapCompactionHeapifyInputKeys  retained keys passed to make_heap
                         after the rebuild scan.
        ProjectionWordXors  64-bit words XORed while accumulating the affine
                         projection of a peeled column onto inactive columns.
        ResidualCoeffWordXors  64-bit words XORed by packed GF(2) residual
                         elimination (mixed or certified completion).
        ResidualCoeffByteOps  coefficient bytes touched by GF(256) residual
                         or quotient elimination.  This includes relation
                         updates, fixed-H column-span elimination, and
                         fixed-H syndrome dot products; bit probes used to
                         select a relation are not byte operations.
                         Equivalent coefficient representations may move work
                         between the word and byte counters; terminal
                         classifiers may also eliminate payload-block work.
    */
    uint64_t BlockCopies = 0;
    uint64_t BlockZeroFills = 0;
    // gf256_addset_multi_mem writes dst = XOR(sources) in ONE pass, so the
    // kernel call is neither a BlockXor (which accumulates into an existing
    // dst) nor another BlockCopy.  Its distinguished first source retains the
    // caller-attributed logical BlockCopy documented above.
    // Profiling at a 64 KiB payload put it at 8.0% of runtime with no counter
    // behind it, which is why the cost model could not see it.
    // BlockAddSets counts the calls; BlockAddSetSources counts the source
    // blocks consumed, since cost scales with the latter.
    uint64_t BlockAddSets = 0;
    uint64_t BlockAddSetSources = 0;
    uint64_t PeelAdjacencyVisits = 0;
    uint64_t PeelRowScanSteps = 0;
    uint64_t PeelHeapOperations = 0;
    uint64_t ProjectionWordXors = 0;
    uint64_t ResidualCoeffWordXors = 0;
    uint64_t ResidualCoeffByteOps = 0;
    uint64_t BuildNanoseconds = 0;
    uint64_t PeelNanoseconds = 0;
    uint64_t ProjectNanoseconds = 0;
    uint64_t ResidualNanoseconds = 0;
    uint64_t BackSubNanoseconds = 0;
    uint32_t PacketSeedAttempt = 0;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    (defined(WIREHAIR_V2_PEEL_HEAP_COMPACTION_DENSITY_PERCENT) && \
     WIREHAIR_V2_PEEL_HEAP_COMPACTION_DENSITY_PERCENT != 0)
    // Experiment receipts.  Keep these absent when both test hooks and heap
    // compaction are disabled: the zero experiment must preserve the shipped
    // stats and resume-state layout exactly.  Enabled production experiments
    // still expose the scan and heapify work needed by physical A/B tooling.
    uint64_t PeelHeapCompactionRebuildColumnProbes = 0;
    uint64_t PeelHeapCompactionHeapifyInputKeys = 0;
#endif
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    // Additional experiment-only lazy-heap diagnostics.  The shipped stats
    // layout and ABI do not contain these fields when test hooks are disabled.
    uint64_t PeelHeapResolvedStalePops = 0;
    uint64_t PeelHeapUnresolvedCountMismatchPops = 0;
    // The density gate is checked once at each eligible pre-compaction
    // queue-drained heap-selection seam.  After the configured invalid-root
    // proof pops, the heap must be guaranteed to contain CutoffKeys entries
    // beyond the maximum possible rebuild output before compaction may
    // trigger once.
    uint32_t PeelHeapCompactionDensitySeamChecks = 0;
    uint64_t PeelHeapCompactionDensityMaxInputKeys = 0;
    uint64_t PeelHeapCompactionDensityCutoffKeys = 0;
    uint32_t PeelHeapCompactionStalePopThreshold = 0;
    uint32_t PeelHeapCompactionDensityQualifiedSeams = 0;
    // ProofArmed requires enough density slack to survive T actual pops and
    // still leave at least CutoffKeys entries beyond the maximum possible
    // rebuild output (one per remaining unresolved column).
    // Marginal dense seams defer without entering the proof-aware loop.
    uint32_t PeelHeapCompactionProofArmedSeams = 0;
    uint32_t PeelHeapCompactionInsufficientProofPopSlackDeferrals = 0;
    uint32_t PeelHeapCompactionInsufficientRemovableKeyExcessDeferrals = 0;
    uint64_t PeelHeapCompactionProofStalePops = 0;
    uint32_t PeelHeapCompactionInsufficientProofDeferrals = 0;
    uint32_t PeelHeapCompactionStalePopThresholdCandidates = 0;
    uint32_t PeelHeapCompactionStalePopThresholdTriggers = 0;
    // Trigger-only state: exact same-episode proof prefix, post-pop heap size,
    // and unresolved columns immediately before the exact rebuild.
    uint32_t PeelHeapCompactionTriggerProofStalePops = 0;
    uint64_t PeelHeapCompactionTriggerPostPopKeys = 0;
    uint32_t PeelHeapCompactionTriggerRemainingColumns = 0;
    uint32_t PeelHeapCompactions = 0;
    uint64_t PeelHeapCompactionInputKeys = 0;
    uint64_t PeelHeapCompactionOutputKeys = 0;
    // Per cold solve, incremented exactly once when the tiny mixed completion
    // solver accepts the quotient shape and supplies the terminal result.
    // Forced-path tests use this to distinguish real engagement from a
    // byte-identical fallthrough to the general solver.
    uint32_t TinyMixedFastPathAcceptances = 0;
    // Representation diagnostics, not wire or deterministic-work contracts.
    uint32_t CertifiedPackedResidualUses = 0;
    uint64_t ResidualCoefficientStorageBytes = 0;
    uint32_t CertifiedPackedResumeMaterializations = 0;
    // Certified GF(256) completion diagnostics.  The first two are specific
    // to the q>H classifier; the remaining fields also describe the q<=H and
    // legacy paths.  Together they prove the terminal path retained only
    // fixed H-by-H coefficient state, allocated no square quotient, skipped
    // heavy RHS work when the quotient already spans all H rows, and replayed
    // original heavy rows at most once when publishing a byte checkpoint.
    uint32_t CertifiedFixedHQuotientUses = 0;
    uint64_t CertifiedFixedHCoefficientStorageBytes = 0;
    uint32_t CertifiedSquareQuotientColumns = 0;
    uint32_t CertifiedHeavyRhsRowsBuilt = 0;
    uint32_t CertifiedLegacyHeavyRowsReplayed = 0;
    uint64_t MixedJointSourceXors = 0;
    uint64_t MixedJointMarginalXors = 0;
    uint64_t MixedJointMarginalCopies = 0;
    uint64_t MixedJointScratchBytes = 0;
    uint32_t MixedJointActiveDeltas = 0;
    uint64_t MixedDualSourceColumns = 0;
#endif
};

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/** Largest mixed quotient reconstructed by the null-witness diagnostic. */
static const uint32_t kMaxMixedNullWitnessQuotientColumns = 15u;

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
    Per-thread override of the PEELING (packet) row degree distribution --
    weights for degrees 1, 2, 3, ... -- replacing the frozen table in
    WirehairTools.cpp.  Empty clears the override.

    THIS IS THE LT OUTPUT DEGREE DISTRIBUTION, the one soliton theory is about.
    It is not the staircase row-degree knob in WirehairV2Precode.h: that shapes
    the S-row sparse ladder (S is 5 to 10 at K=128), while this shapes the K
    peeling rows the peeler actually peels.  The two were conflated for a whole
    tuning campaign, so the distinction is spelled out here rather than assumed.

    Thread-local, matching SetStaircaseDegreesForTesting, because a concurrent
    search evaluates a different distribution per worker; a parse-once
    environment global cannot do that.  WIREHAIR_V2_PEEL_DEGREES is honoured as
    a fallback for single-config harnesses, and the thread-local wins.

    SAMPLING IS STREAM-PRESERVING.  The override consumes the same PRNG draw
    the stock weight generator consumes and maps it through the supplied CDF, so
    supplying the stock PMF reproduces the stock construction BIT-IDENTICALLY
    rather than landing on a different random stream.  That identity is the gate
    that distinguishes a real effect from a broken hook.  CmdSelfTest checks the
    recovered native PMF against the generator on uniform probes and every exact
    CDF boundary, including the override sampler itself.

    NOTE that the stock PMF is NOT the textbook ideal soliton, and feeding the
    textbook version does NOT reproduce it.  GeneratePeelRowWeight takes the
    weight-1 mass by `rv -= P1`, a SHIFT rather than a carve-out, so that mass
    comes out of the deep tail; and the table's last entry is 0xffffffff, so
    ALL remaining mass is lumped at d = 64.  Before the later N/2 clamp,
    exactly:

        P1_eff = N <= 2048 ? floor((1/128) * 0xffffffff) : 0
        P(1)  = P1_eff / 2^32
        P(2)  = (T[1] + 1) / 2^32  = 1/2
        P(i)  = (T[i-1] - T[i-2]) / 2^32       for 3 <= i <= 63
        P(64) = (2^32 - P1_eff - 1 - T[62]) / 2^32
              = 0.008061 for N <= 2048, 0.015873 otherwise

    with T = kPeelCountDistribution.

    The d = 64 atom is large -- about thirty-two times the ideal soliton's
    1/(64*63) = 0.000248 through N=2048 and sixty-four times larger above
    N=2048 -- because the table lumps all residual tail mass there.  It is NOT
    load-bearing.  An earlier version of this comment claimed it was, on the
    strength of a full-codec run in
    which replacing that one atom took decoding from 150/150 pass to 150/150
    fail.  That measurement was an artefact of the bug described below, and
    with the bug fixed the same substitution decodes 150/150.  Recorded here
    because a comment asserting an artefact, in the file a reader reaches
    first, is worse than no comment.

    THE BUG, and the reason both consumers now use ONE canonical initializer.
    The override was originally installed only in the DECODER's row builder.
    The ENCODER built its rows with PeelRowParameters::Initialize directly.
    So with any override active the encoder emitted packets built from stock
    degrees while the solver built equations from overridden ones; every block
    id whose degree moved produced a wrong equation, the solve returned
    Wirehair_Success over a system that did not describe the data, and the
    payload came back corrupt.  Deterministic in block id, hence always 0/N or
    N/N and never partial.  Keeping the policy in one initializer prevents a
    future hook from being added to only one side again.

    WHY THE IDENTITY GATE ABOVE CANNOT CATCH IT.  Stock maps to stock degrees
    on both sides by construction, so the two paths agree and the control
    passes.  The identity arm is a faithful test of ONE hook and is
    structurally blind to an asymmetry BETWEEN two, because its only input is
    the fixed point of both.  It passed for an entire campaign while every
    non-stock distribution was being corrupted.

    SO THERE ARE TWO GATES AND THEY TEST DIFFERENT THINGS.  The identity gate
    proves the hook changes nothing when it should change nothing.  A
    FULL-CODEC gate -- `compare --precode --precode-profile mixed` at OH ~ 0,
    run on a NON-STOCK distribution -- proves the hook changes the right things
    everywhere they are built.  Only the second can detect a missing call site,
    and only a non-stock input exercises it.  Run both.

    The N/2 clamp in PeelRowParameters::Initialize is RESPECTED, not superseded:
    a weight above peel_column_count/2 is clamped exactly as the stock path
    clamps it, so an override cannot silently change the row-weight invariant
    the solver relies on.  Degrees are also clamped into [1, kMaxPeelCount].
*/
/**
    Install and validate a thread-local PMF.  Returns false for non-finite,
    negative, all-zero, or overflowing weights and leaves an active-invalid
    override that makes encoder and decoder row construction fail closed.
    Passing an empty vector is the explicit clear operation and returns true.
*/
bool SetPeelDegreesForTesting(const std::vector<double>& weights);
void ClearPeelDegreesForTesting();
/// Validate the active thread-local/environment arm without sampling it.
bool PeelDegreeConfigurationValidForTesting();

/**
    Sample the active peel override from an exact uint32 uniform value.

    Returns zero when no override is active.  This is the boundary oracle used
    by CmdSelfTest; exposing the sampler prevents a test from reimplementing the
    comparison and missing the same lower_bound/upper_bound mistake.
*/
uint32_t SamplePeelDegreeForTesting(
    uint32_t random_value,
    uint16_t peel_column_count);
/// Fixed number of CDF comparisons performed by every override sample.
uint32_t PeelDegreeSampleComparisonCountForTesting();

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
    Override the calling thread's one-shot lazy-heap density percentage.
    190 means 1.90 times the original binary-system column count.  A qualifying
    queue-drained episode counts only roots the baseline loop actually rejects.
    The count resets at each episode.  Compaction requires enough initial
    density slack to survive the proof prefix and then still guarantee that
    many rebuild-removable entries after reserving one possible output key per
    unresolved column.  Zero in either control disables the experiment.
    Values are captured once at binary-peel entry, and the independent oracle
    always remains off.
*/
void SetBinaryPeelHeapCompactionDensityPercentForTesting(
    uint32_t density_percent);
uint32_t BinaryPeelHeapCompactionDensityPercentForTesting();
void SetBinaryPeelHeapCompactionStalePopThresholdForTesting(
    uint32_t stale_pop_threshold);
uint32_t BinaryPeelHeapCompactionStalePopThresholdForTesting();

/**
    Compare packed-GF(2) residual insertion against the byte GF(256) oracle at
    word boundaries, including poisoned tail bits and inconsistent RHS rows.
*/
bool CheckPackedBinaryResidualOracleForTesting();

/**
    Differentially verify the fixed-H certified q>H classifier against
    independently constructed dense quotients.  Covers full-rank and
    chronologically deficient matrices, packed-word boundaries, byte/packed
    binary bases, and every production heavy-row count from zero through 12.
*/
bool CheckCertifiedFixedHQuotientFactorForTesting();

/**
    Compare fused mixed-pivot and unused-row RHS initialization against the
    legacy copy/zero-plus-XOR sequence across every production quotient width,
    extended test-hook widths, packed-word boundaries, SIMD boundary block
    sizes, null packet data, and residual classifications.
*/
bool CheckMixedRhsFusionOracleForTesting();

/**
    Exercise the production-period H=12, q=13 rectangular mixed quotient.
    The fixed free columns span only residues 0..10: a zero RHS is consistent
    NeedMore, while a known residue-12 value must be rejected by its exact
    left-nullspace syndrome without changing the caller's values.
*/
bool CheckMixedQGreaterThanHSyndromeForTesting();

enum class MixedNullWitnessStatus : uint32_t
{
    None = 0,
    Captured = 1,
    AlgebraMismatch = 2,
    AllocationFailure = 3,
    InternalError = 4
};

/**
    Exact test-only description of ker(A) intersect ker(C) for a singular
    mixed quotient.  CanonicalBasis is row-major over all original solver
    columns, with KernelDimension rows and ColumnCount columns.
*/
struct MixedNullWitnessDiagnostic
{
    MixedNullWitnessStatus Status = MixedNullWitnessStatus::None;
    uint32_t ColumnCount = 0u;
    uint32_t InactiveCount = 0u;
    uint32_t BinaryRank = 0u;
    uint32_t QuotientColumns = 0u;
    uint32_t QuotientRank = 0u;
    uint32_t KernelDimension = 0u;
    bool QuotientVerified = false;
    bool BinaryVerified = false;
    bool CompletionVerified = false;
    bool CanonicalVerified = false;
    uint64_t BasisHashLow = 0u;
    uint64_t BasisHashHigh = 0u;
    std::vector<uint16_t> CanonicalBasis;
    /** One byte per original column: one for inactive, zero for peeled. */
    std::vector<uint8_t> InactiveMask;
};

/** Install a non-owning capture sink for mixed solves on the calling thread. */
void SetMixedNullWitnessDiagnosticForTesting(
    MixedNullWitnessDiagnostic* diagnostic);

/** Pure GF(2^16) nullspace/canonicalization invariance oracle. */
bool CheckMixedNullWitnessCanonicalizationForTesting();

/**
    Compare recorded mixed-quotient row plans with direct coefficient/RHS
    Gauss-Jordan elimination, including the derived left transform.
*/
bool CheckMixedQuotientFactorReplayForTesting();

/**
    Exercise the coefficient-rank-deficient mixed-quotient shortcut with a
    consistent zero RHS and an algebraically guaranteed inconsistent RHS.
*/
bool CheckMixedQuotientDeficientSyndromeForTesting();
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

/**
    Generate a packet row without allocating output storage.

    The function fails closed when the runtime is invalid, `columns_out` is
    null, or `capacity` is smaller than the complete row.  In every failure
    case it leaves both caller storage and `count_out` untouched.
    `count_out` must not overlap the output row; overlap is rejected before
    writing.  On success it writes exactly the same ordered columns as
    GeneratePacketMatrixRowWithRuntime() and stores their count in
    `count_out`.  The caller must initialize any enabled test-hook state
    outside an allocation-sensitive interval; production has no such state.
*/
bool GeneratePacketMatrixRowIntoWithRuntime(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    uint32_t* columns_out,
    size_t capacity,
    size_t& count_out);

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
    Mixed completion and the certified fast path keep the binary residual
    packed in GF(2), then solve only its remaining quotient over GF(2^16) or
    GF(256), respectively.  A narrow certified solve that can retain a legacy
    resume checkpoint keeps its byte basis to avoid conversion; test hooks
    force either representation for differential checks.  On success
    `intermediate_blocks_out` contains all
    K+S+D2+H block values.  The exact residual solve is bounded to
    kMaxInactiveColumns to contain adversarial memory use.  NeedMore means the
    supplied equations were rank deficient or exceeded that bound; additional
    independent packets can reduce the residual.  Output remains unchanged on
    every failure.  When resume_state is non-null, a rank-deficient residual
    within the cap and persistent-byte budget atomically replaces it with an
    active affine/pivot checkpoint; cap/budget failures leave it unchanged and
    require a cold retry.
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

    resume_persistent_byte_limit is a conservative checkpoint budget.  A
    deficient solve may return NeedMore without publishing a resume state when
    the minimum legacy checkpoint shape already exceeds it.  The maximum
    size_t value preserves the ordinary unbounded solver contract.
*/
WirehairResult SolvePrecodeSystemForValidatedSystemWithRuntime(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats = nullptr,
    PrecodeSolveResumeState* resume_state = nullptr,
    size_t resume_persistent_byte_limit = (size_t)-1);

/**
    Append one packet equation to a rank-deficient solve checkpoint.

    With allow_insert=false this performs a non-mutating duplicate consistency
    check.  With allow_insert=true an independent row is committed and Success
    is returned as soon as the complete intermediate vector is reconstructed.
    Allocations finish before an inserting call changes the algebraic state, so
    OOM is retryable.  On OOM, stats receives the unchanged checkpoint counters
    when non-null.  Output remains unchanged on NeedMore and every failure.
    A zero block width is the supported payload-free receipt mode; block_data
    must still be non-null, but copied empty checkpoint vectors need no backing
    capacity.
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
