#pragma once

#include "WirehairV2Peel.h"
#include "WirehairV2Precode.h"
#include "WirehairV2Seeds.h"
#include "WirehairV2Solve.h"

#include <wirehair/wirehair.h>

#include <stddef.h>
#include <stdint.h>
#include <vector>

/*
    Encoder value phase for the certified V2 precode (wirehair-axd Phase 3).

    Given the K source block values, computes the S + D2 + H intermediate
    parity block values so that every constraint row of the PrecodeSystem
    sums to zero over the full column value vector
    [source | staircase | dense | heavy]:

    - STAIRCASE rows solve by one forward pass: parity j is the XOR of its
      row's other columns (source hits + the j-1 link), exactly
      min(N1,S)*K + S - 1 block ops.
    - Shuffle-2 DENSE rows couple all D2 dense parity values, because every
      row spans source, staircase AND dense columns.  The rows are
      transformed to consecutive-row differences (an invertible unit
      bidiagonal row transform, so solvability and the solution are
      unchanged): the difference of two consecutive rows is exactly the two
      flipped columns, so each transformed row's known part costs at most 2
      block ops -- the incremental structure behind the certified
      ceil(span/2) + 2*(D2-1) cost model.  What remains is a D2 x D2 GF(2)
      system over the dense columns, solved by Gauss-Jordan with block RHS.
      That corner is NOT guaranteed invertible -- when singular the
      constraints have no (unique) solution for these source values and the
      function returns false.  Use DenseCornerInvertible to measure or gate
      on this without touching block data.
    - Cauchy HEAVY rows: the known part is a GF(256) muladd accumulation
      over the K + S + D2 known columns; the H x H corner over the heavy
      columns (consecutive GE columns, always distinct within one
      256 - H = 244 window) is Cauchy, hence invertible, and is solved by
      GF(256) Gauss-Jordan over blocks.

    Cost accounting (PrecodeEncodeStats) counts every block-sized memory
    operation: a first-term copy into an accumulator counts like an XOR,
    matching the Phase B `precode_gen_xors_mu` convention; zero-filling an
    empty accumulator counts as nothing.
*/

namespace wirehair_v2 {

static const uint32_t kDefaultRecoveryMixCount = kCertifiedPacketMixCount;
static const uint64_t kMessagePrecodeSeedSalt =
    UINT64_C(0x763263707265636f);
static const uint64_t kMessageRecoveryRowSeedSalt =
    UINT64_C(0x76327265636f7665);

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
// Fails the next guarded allocation when countdown is zero; negative disables.
void SetAllocationFailureCountdownForTesting(int64_t countdown);
// Overrides the full residue-bucket scratch cap; production defaults to 64 MiB.
void SetHeavyBucketStorageLimitForTesting(uint64_t bytes);
#endif

struct PrecodeEncodeStats
{
    /// Staircase forward pass block ops
    /// (certified model: min(N1,S)*K + S - 1)
    uint64_t StaircaseBlockOps = 0;

    /// Dense difference-row known-part block ops
    /// (certified model: ceil(span/2) + 2*(D2 - 1), counted even when the
    /// dense corner turns out singular)
    uint64_t DenseKnownBlockOps = 0;

    /// GF(2) Gauss-Jordan block XORs + solution copies for the dense corner
    uint64_t DenseSolveBlockOps = 0;

    /// Plain block XOR-class contributions folding known columns into heavy
    /// residue buckets.  The joint-delta path excludes its write-only
    /// first-marginal copies, reported separately below.
    uint64_t HeavyBucketXors = 0;

    /// GF(256) muladd-class block ops accumulating the heavy known parts
    /// (model: H * min(K + S + D2, 256 - H) bucketed, H * (K + S + D2)
    /// direct)
    uint64_t HeavyMulAdds = 0;

    /// Row-level block ops (muladd/div/copy) solving the H x H completion
    /// corner.  These are GF(256) for the legacy profile and planar
    /// GF(2^16) for the mixed profile.
    uint64_t HeavySolveBlockOps = 0;

    /// Subset of HeavyMulAdds evaluated with extension-field coefficients.
    /// The other mixed-profile rows are ordinary GF(256) muladds.
    uint64_t MixedGF16MulAdds = 0;

    /// Mixed-profile planar GF(2^16) scale/muladd/copy operations used by
    /// the 12 x 12 completion solve.
    uint64_t MixedGF16SolveBlockOps = 0;

    /// Full-block deinterleave/interleave conversions in the mixed path.
    uint64_t MixedPlaneConversions = 0;

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    /// Joint-delta experiment counters.  These stay zero unless that
    /// explicitly selected test path is used.
    uint64_t MixedJointSourceXors = 0;
    uint64_t MixedJointMarginalXors = 0;
    uint64_t MixedJointMarginalCopies = 0;
    uint64_t MixedJointScratchBytes = 0;
    uint32_t MixedJointActiveDeltas = 0;
    uint64_t MixedDualSourceColumns = 0;
#endif
};

/**
    Structure-only feasibility check: GF(2) rank of the D2 x D2 submatrix of
    the dense rows over the dense columns [K + S, K + S + D2).  When this is
    rank-deficient, ComputePrecodeValues must fail for generic source data
    (the staircase part is always solvable and the heavy corner is always
    invertible, so this corner is the ONLY feasibility gate).  Requires
    DenseRows <= 64.  D2 == 0 counts as invertible.
*/
bool DenseCornerInvertible(const PrecodeSystem& system);

/**
    Compute the S + D2 + H intermediate parity block values.

    source_blocks: K contiguous blocks of block_bytes each.
    parity_blocks: out, S + D2 + H contiguous blocks -- staircase parities
    [0, S), dense parities [S, S + D2), heavy parities [S + D2, S + D2 + H),
    i.e. parity i is global column K + i.

    Returns false on invalid arguments (block_bytes 0 or > 2^31 - 1,
    DenseRows > 64, HeavyRows > 128, row vectors inconsistent with Params)
    or when the D2 x D2 dense corner is singular.  On the singular-corner
    path stats (if given) still hold the staircase and dense known-part
    counts, so cost measurements do not require a feasible seed.
*/
bool ComputePrecodeValues(
    const PrecodeSystem& system,
    const uint8_t* source_blocks,
    uint32_t block_bytes,
    uint8_t* parity_blocks,
    PrecodeEncodeStats* stats = nullptr);

/**
    Compute one V2 recovery symbol from a row over the intermediate domain
    [source | staircase | dense | heavy].

    `row_columns` uses the same global column numbering as
    GenerateRecoveryMatrixRows(): columns < K reference source_blocks, and
    columns >= K reference parity_blocks at offset column - K.  The output is
    the GF(2) XOR of all referenced block values.  A first-term copy counts as
    one block op, matching the precode cost convention; an empty row produces
    the zero block and zero ops.

    `block_out` must not overlap the source/parity block arrays.
*/
bool ComputeRecoveryBlock(
    const PrecodeSystem& system,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    const std::vector<uint32_t>& row_columns,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/**
    Compute one encoder output block.

    block_id < K copies the corresponding source block.  block_id >= K treats
    block_id - K as the zero-based recovery row index, generates that row with
    GenerateRecoveryMatrixRow(), then evaluates it with ComputeRecoveryBlock().
    `row_seed` is the seed for the recovery row stream, and `mix_count` is the
    number of precode columns appended to each recovery row.

    `parity_blocks` may be null only for source block ids; recovery ids require
    the S + D2 + H parity blocks already computed by ComputePrecodeValues().
    `block_out` must not overlap the source/parity block arrays.
*/
bool ComputeEncodedBlock(
    const PrecodeSystem& system,
    const PeelingCodec& codec,
    uint64_t row_seed,
    uint32_t mix_count,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/**
    Encoder-side state for the V2 precode block-data path.

    The source block array is borrowed, not copied, and must remain valid while
    the encoder is used.  The computed S + D2 + H parity blocks are owned by
    this object.  Initialize returns false for invalid arguments or
    encoder-infeasible precode systems, e.g. a singular dense corner.
    Accessors return zero/null/default state until initialization succeeds.
    A failed reinitialization preserves the last successfully initialized
    state and its accessors.
*/
class PrecodeEncoder
{
public:
    PrecodeEncoder();

    bool Initialize(
        const PrecodeSystem& system,
        const PeelingCodec& codec,
        uint64_t row_seed,
        uint32_t mix_count,
        const uint8_t* source_blocks,
        uint32_t block_bytes);

    WirehairResult InitializeResult(
        const PrecodeSystem& system,
        const PeelingCodec& codec,
        uint64_t row_seed,
        uint32_t mix_count,
        const uint8_t* source_blocks,
        uint32_t block_bytes);

    bool Encode(
        uint32_t block_id,
        uint8_t* block_out,
        uint64_t* block_ops_out = nullptr) const;

    /// On failure block_out and block_ops_out are left unchanged.
    WirehairResult EncodeResult(
        uint32_t block_id,
        uint8_t* block_out,
        uint64_t* block_ops_out = nullptr) const;

    bool IsInitialized() const;
    uint32_t SourceBlockCount() const;
    uint32_t ParityBlockCount() const;
    uint32_t BlockBytes() const;
    uint64_t RecoveryRowSeed() const;
    uint32_t RecoveryMixCount() const;
    const PrecodeEncodeStats& EncodeStats() const;
    const uint8_t* ParityBlocks() const;
    const uint8_t* IntermediateBlocks() const;

    /**
        True when System() retains the complete validated row graph.

        The older direct InitializeResult path needs those rows for recovery
        encoding and preserves the original complete-System contract.  The
        solved packet-contract path used by MessagePrecodeEncoder needs only
        Params after validation and deliberately releases the row graph.
    */
    bool HasCompleteSystem() const;

    /**
        Returns the retained system descriptor.  Params are complete for every
        initialized encoder; StaircaseRows and DenseRowColumns are complete
        only when HasCompleteSystem() is true.
    */
    const PrecodeSystem& System() const;

private:
    friend class MessagePrecodeEncoder;
    void Swap(PrecodeEncoder& other) noexcept;
    WirehairResult InitializeSolvedSystem(
        const PrecodeSystem& system,
        const PacketRowConfig& packet_config,
        std::vector<uint8_t>& intermediate_blocks,
        uint32_t block_bytes);

    PrecodeSystem SystemValue = {};
    PeelingCodec CodecValue = {};
    uint64_t RowSeed = 0;
    uint32_t MixCount = 0;
    const uint8_t* SourceBlocks = nullptr;
    uint32_t BlockBytesValue = 0;
    std::vector<uint8_t> ParityBlockStorage;
    std::vector<uint8_t> SolvedIntermediateStorage;
    PacketRowConfig PacketConfigValue = {};
    PacketRowRuntime PacketRuntimeValue = {};
    PrecodeEncodeStats StatsValue = {};
    bool UsesPacketContract = false;
    bool Initialized = false;
};

struct MessagePrecodeEncoderOptions
{
    // Named version-4 packet contracts bind this value together with the
    // completion field: three columns for the original profiles, or two for
    // the opt-in mixed/mix2 profile.  Keep it explicit so supplied options can
    // be checked against a selected profile rather than silently ignored.
    uint32_t RecoveryMixCount = kDefaultRecoveryMixCount;
    bool DenseIdentityCorner = false;
    // Versioned equation policy: retain certified D12 below K=4096 and use
    // the second shuffled-deck anchor at and above that boundary.
    bool AdaptiveDenseTwoAnchor = false;
    uint64_t PrecodeSeedSalt = kMessagePrecodeSeedSalt;
    uint64_t RecoveryRowSeedSalt = kMessageRecoveryRowSeedSalt;
    CompletionField Completion = CompletionField::GF256;

    // Local codec policy, deliberately excluded from the serialized packet
    // contract.  The encoder retains exact message bytes for direct packet
    // copies.  The decoder retains accepted systematic payloads so Recover
    // evaluates only missing source ids.  Either cache may be released.
    bool CacheSystematicSource = false;
    bool CacheReceivedSystematicPackets = false;
};

/** True when any selected or mixed V2 precode contract state is present. */
bool HasMessagePrecodeContractState(const SeedProfile& profile);

/**
    Resolve and validate the matrix options bound to a selected V2 profile.

    An unselected profile uses requested_options (or the defaults when null).
    A selected profile carries the complete versioned packet/precode contract;
    null options inherit that contract and explicit options must match it.
*/
bool ResolveMessagePrecodeOptions(
    const SeedProfile& profile,
    const MessagePrecodeEncoderOptions* requested_options,
    MessagePrecodeEncoderOptions& resolved_options);

/** Resolve exact precode dimensions and packet seeds from the profile. */
bool ResolveMessagePrecodeConfiguration(
    const SeedProfile& profile,
    const MessagePrecodeEncoderOptions& options,
    PrecodeParams& params,
    PacketRowConfig& packet_config);

/** Publish the complete matrix contract and selected seed attempt. */
void BindMessagePrecodeProfile(
    SeedProfile& profile,
    const MessagePrecodeEncoderOptions& options,
    const PrecodeSystem& system,
    const PacketRowConfig& packet_config,
    uint32_t packet_seed_attempt);

/**
    Message-level adapter for the V2 precode encoder.

    This borrows complete caller-owned message blocks only for the duration of
    synchronous initialization, owns and zero-pads at most one partial final
    block while solving, and owns the resulting full intermediate vector.
    The initialized encoder is therefore independent of the message buffer
    after Initialize returns.  By default it keeps no source copy; the
    runtime-only CacheSystematicSource option retains exact message bytes for
    direct systematic copies until ReleaseSystematicSourceCache is called.
    It derives the certified precode system and
    packet seed from a SeedProfile and exposes V1-style encode semantics:
    systematic block ids emit only the original byte count for the final
    partial block, while recovery block ids emit a full block.

    The packet degree distribution is fixed by the version-4 contract to the
    production Wirehair integer sampler.  Each named profile also binds its
    precode mix count: three for the original profiles and two for the opt-in
    mixed/mix2 profile, rather than an unbound runtime policy.  The default
    path preserves the certified full-span dense rows.  It solves
    K deterministic systematic packet equations together with all precode
    constraints for the complete intermediate vector, so it does not require
    the dense parity corner to be independently invertible.  The
    DenseIdentityCorner option remains an experimental oracle variant.

    Initialization is transactional: any invalid, singular, or allocation
    failure preserves the last successfully initialized encoder.
*/
class MessagePrecodeEncoder
{
public:
    MessagePrecodeEncoder();
    MessagePrecodeEncoder(const MessagePrecodeEncoder&) = delete;
    MessagePrecodeEncoder& operator=(const MessagePrecodeEncoder&) = delete;

    bool Initialize(
        const void* message,
        uint64_t message_bytes,
        uint32_t block_bytes,
        const SeedProfile* seed_override = nullptr,
        const MessagePrecodeEncoderOptions* options = nullptr);

    WirehairResult InitializeResult(
        const void* message,
        uint64_t message_bytes,
        uint32_t block_bytes,
        const SeedProfile* seed_override = nullptr,
        const MessagePrecodeEncoderOptions* options = nullptr);

    bool Encode(
        uint32_t block_id,
        uint8_t* block_out,
        uint32_t out_bytes,
        uint32_t* data_bytes_out,
        uint64_t* block_ops_out = nullptr) const;

    /// On failure all output buffers and output counters are left unchanged.
    WirehairResult EncodeResult(
        uint32_t block_id,
        uint8_t* block_out,
        uint32_t out_bytes,
        uint32_t* data_bytes_out,
        uint64_t* block_ops_out = nullptr) const;

    bool IsInitialized() const;
    uint64_t MessageBytes() const;
    uint32_t SourceBlockCount() const;
    uint32_t BlockBytes() const;
    const SeedProfile& Profile() const;
    const MessagePrecodeEncoderOptions& Options() const;
    const PrecodeEncodeStats& EncodeStats() const;
    const PrecodeSolveStats& SolveStats() const;
    /** Complete solved [source-domain | precode] intermediate vector. */
    const uint8_t* IntermediateBlocks() const;
    const PrecodeEncoder& BlockEncoder() const;

    /** Release the optional source cache.  Safe to call repeatedly. */
    void ReleaseSystematicSourceCache() noexcept;
    bool HasSystematicSourceCache() const;
    size_t SystematicSourceCacheBytes() const;

private:
    SeedProfile ProfileValue = {};
    MessagePrecodeEncoderOptions OptionsValue = {};
    PrecodeEncoder EncoderValue;
    PrecodeSolveStats SolveStatsValue = {};
    uint64_t MessageBytesValue = 0;
    uint32_t BlockBytesValue = 0;
    std::vector<uint8_t> SystematicSourceCache;
    bool Initialized = false;
};

} // namespace wirehair_v2
