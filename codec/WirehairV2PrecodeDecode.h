#pragma once

#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"

#include <wirehair/wirehair.h>

#include <stddef.h>
#include <stdint.h>
#include <memory>
#include <vector>

namespace wirehair_v2 {

static const uint32_t kDecoderInitialReceiveOverhead = 32u;
static const uint32_t kDecoderMaximumReceiveOverhead = 1024u;

/**
    Requested cold-receive capacities used by MessagePrecodeDecoder.  These
    model the initialized packet reserve and one pending full block.  The block
    count must be in the public codec domain; an invalid domain, overflow, or
    zero block width returns false.  Allocator-exact accounting uses the
    capacities reported by the decoder's vectors after these reserves.
*/
bool DecoderInitialReceiveCapacities(
    uint32_t block_count,
    uint32_t block_bytes,
    size_t& receive_block_capacity,
    size_t& receive_id_capacity,
    size_t& pending_block_capacity);

/** Exact 25%-growth policy for the persistent resume state. */
bool DecoderResumePersistentByteLimit(
    size_t receive_block_capacity,
    size_t receive_id_capacity,
    size_t pending_block_capacity,
    size_t& persistent_byte_limit);

/**
    Bounded flat hash table used by the decoder's accepted-packet set.

    Every uint32_t value is a valid public packet id, so occupancy is stored
    separately rather than stealing a key as an empty sentinel.  Initialize()
    provisions the normal receive window; rare tail growth is transactional
    and bounded by the decoder's maximum accepted-equation count.
*/
class PacketSlotTable
{
public:
    bool Initialize(size_t initial_entries, size_t max_entries);
    bool Find(uint32_t packet_id, uint32_t* slot_out = nullptr) const;
    /** Allocation failure returns false and leaves the table unchanged. */
    bool Insert(uint32_t packet_id, uint32_t slot);
    bool Erase(uint32_t packet_id);
    void ClearAndRelease() noexcept;
    void Swap(PacketSlotTable& other) noexcept;

    size_t Size() const { return EntryCount; }
    size_t Capacity() const { return Keys.size(); }
    /** Saturating sum of all backing-vector capacities in bytes. */
    size_t StorageBytes() const;

private:
    friend class MessagePrecodeDecoder;

    static uint32_t Hash(uint32_t packet_id);
    bool InsertImpl(uint32_t packet_id, uint32_t slot);
    void Grow();

    std::vector<uint32_t> Keys;
    std::vector<uint32_t> Slots;
    std::vector<uint8_t> Occupied;
    size_t EntryCount = 0u;
    size_t EntryLimit = 0u;
};

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void SetDecoderAllocationFailureCountdownForTesting(int64_t countdown);
void SetDecoderIncrementalResumeEnabledForTesting(bool enabled);
void SetDecoderColdSystematicReuseEnabledForTesting(bool enabled);

/**
    Coarse receive-path evidence for tests that must distinguish the real
    message decoder from direct solver adapters.  Counters are thread-local
    and updated only at coarse decoder phase boundaries, never on each
    ordinary cold-receive packet.
*/
struct DecoderReceivePathCounters
{
    uint64_t ValidatedSystemInitializations = 0u;
    uint32_t LastValidatedPacketSeedAttempt = UINT32_MAX;
    uint64_t DirectSystematicCompletions = 0u;
    uint64_t DirectSystematicCanonicalizations = 0u;
    uint64_t DirectSystematicMaterializationAttempts = 0u;
    uint64_t ColdSolveAttempts = 0u;
    uint64_t ColdSolvePacketAssemblies = 0u;
    uint64_t ColdSolveSlotEntries = 0u;
    uint64_t ColdSolvePayloadBytes = 0u;
    uint64_t CheckpointPendingAllocationAttempts = 0u;
    uint64_t CheckpointAdoptions = 0u;
    uint64_t PendingPacketCopies = 0u;
    uint64_t PendingPacketCopyBytes = 0u;
    uint64_t ResumeAttempts = 0u;
    uint64_t RecoveryPreflights = 0u;
    uint64_t RecoveryPacketEvaluations = 0u;
    uint64_t RecoveryColdPacketCopies = 0u;
    uint64_t RecoveryColdPacketCopyBytes = 0u;
    uint64_t RecoveryColdStorageReleases = 0u;
};

void ResetDecoderReceivePathCountersForTesting();
DecoderReceivePathCounters DecoderReceivePathCountersForTesting();
#endif

/**
    Incremental message decoder for the version-4 packet/precode path.

    NeedMore preserves an algebraically exact sparse-projection and reduced
    residual checkpoint when it fits the decoder memory policy.  Later unique
    packets are projected and inserted without rebuilding the original solve.
    Oversized residuals or checkpoints that would exceed the bounded memory
    policy fall back to cold re-solves.  Packet authentication is a caller
    concern:
    exactly K independent equations with altered payloads can define another
    valid message.  Duplicates and packets checked against a completed or
    overdetermined solution are validated for consistency.  An inconsistent
    accepted packet set remains poisoned because an erasure decoder cannot
    identify which independent packet was corrupt; reinitialize and replay a
    caller-authenticated packet set after `Wirehair_Error`.
*/
class MessagePrecodeDecoder
{
public:
    MessagePrecodeDecoder();
    MessagePrecodeDecoder(const MessagePrecodeDecoder&) = delete;
    MessagePrecodeDecoder& operator=(const MessagePrecodeDecoder&) = delete;

    WirehairResult InitializeResult(
        uint64_t message_bytes,
        uint32_t block_bytes,
        const SeedProfile* seed_override = nullptr,
        const MessagePrecodeEncoderOptions* options = nullptr);

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    /**
        Transactionally initialize the normal receive state around one exact,
        already validated native equation system.  This test-only seam skips
        public seed selection but otherwise shares decoder setup and every
        DecodeResult/RecoverResult path with production.  The system is copied
        so the initialized decoder does not borrow benchmark-arm storage.
    */
    WirehairResult InitializeForValidatedSystemForTesting(
        uint64_t message_bytes,
        uint32_t block_bytes,
        const PrecodeSystem& system,
        const PacketRowConfig& packet_config,
        uint32_t packet_seed_attempt);
#endif

#if defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
    /**
        Benchmark-only counterpart of the exact validated-system initializer.
        This narrow capability preserves experimental equations without
        compiling generic test hooks or diagnostic counters into timed code.
    */
    WirehairResult InitializeForValidatedSystemForBenchmark(
        uint64_t message_bytes,
        uint32_t block_bytes,
        const PrecodeSystem& system,
        const PacketRowConfig& packet_config,
        uint32_t packet_seed_attempt);
#endif

    /** Identical duplicate ids are ignored; conflicting duplicates fail. */
    WirehairResult DecodeResult(
        uint32_t block_id,
        const void* block_in,
        uint32_t data_bytes);

    /**
        Recover requires the exact initialized message size.  The output range
        must not overlap decoder-owned payload or intermediate-block storage;
        overlap is rejected before writing.  Allocation failure for a partial
        final block also leaves the complete output range untouched.  An exact
        all-systematic completion retains its source payload for repeatable
        direct recovery.  Its first successful recovery canonicalizes arrival-
        ordered slots; materially sized id/hash metadata is released while a
        bounded tiny map may remain allocated but unused.  Later direct
        recoveries are one bulk copy.  After an ordinary solve, the first
        successful call may consume retained cold receive slots; later calls
        evaluate packets no longer available from that one-shot state.
    */
    WirehairResult RecoverResult(
        void* message_out,
        uint64_t message_bytes) const;

    bool IsInitialized() const;
    bool IsDecoded() const;
    uint32_t ReceivedCount() const;
    uint32_t SolveAttemptCount() const;
    uint32_t PacketSeedAttempt() const;
    uint32_t PacketPeelSeed() const;
    uint64_t MessageBytes() const;
    uint32_t BlockBytes() const;
    const SeedProfile& Profile() const;
    const MessagePrecodeEncoderOptions& Options() const;
    const PrecodeSolveStats& SolveStats() const;
    const PrecodeSystem& System() const;
    /**
        Returns null for a decoded all-systematic message until a later repair
        packet requires lazy materialization of the ordinary solution.
    */
    const uint8_t* IntermediateBlocks() const;

    /** Release the optional received-systematic cache; idempotent. */
    void ReleaseSystematicPacketCache() noexcept;
    bool HasSystematicPacketCache() const;
    size_t SystematicPacketCacheBytes() const;
    uint32_t CachedSystematicPacketCount() const;

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    bool HasIncrementalResumeStateForTesting() const;
    size_t IncrementalResumeBytesForTesting() const;
    /** Payload-plus-id capacity used as the resume-policy release basis. */
    size_t ColdReceiveCapacityBytesForTesting() const;
    /** Saturating total of payload, id, and packet-slot allocations. */
    size_t ColdReceiveAllocationBytesForTesting() const;
    /** Saturating total of id and packet-slot mapping allocations. */
    size_t ColdReceiveMetadataBytesForTesting() const;
    /** Borrowed cold-payload pointer; invalidated by decoder mutation. */
    const uint8_t* ColdReceiveStorageForTesting() const
    {
        return ReceivedBlockStorage.empty() ?
            nullptr : ReceivedBlockStorage.data();
    }
    /** Borrowed optional systematic-cache pointer; null when disabled. */
    const uint8_t* SystematicPacketCacheDataForTesting() const
    {
        return SystematicPacketCache.get();
    }
#endif

private:
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
    WirehairResult InitializeForValidatedSystem(
        uint64_t message_bytes,
        uint32_t block_bytes,
        const PrecodeSystem& system,
        const PacketRowConfig& packet_config,
        uint32_t packet_seed_attempt);
#endif
    WirehairResult InitializeResolvedSystem(
        uint64_t message_bytes,
        uint32_t block_bytes,
        SeedProfile profile,
        MessagePrecodeEncoderOptions options,
        PrecodeSystem&& system,
        const PacketRowConfig& packet_config,
        uint32_t packet_seed_attempt);
    WirehairResult AttemptSolve();
    WirehairResult MaterializeDirectSystematic();
    WirehairResult ValidateDecodedPacketResult(
        uint32_t block_id,
        const void* block_in,
        uint32_t data_bytes);
    WirehairResult RecoverDirectSystematicResult(
        void* message_out,
        uint64_t message_bytes) const;
    void Swap(MessagePrecodeDecoder& other) noexcept;

    SeedProfile ProfileValue = {};
    MessagePrecodeEncoderOptions OptionsValue = {};
    PacketRowConfig PacketConfigValue = {};
    PacketRowRuntime PacketRuntimeValue = {};
    PrecodeSystem SystemValue = {};
    mutable std::vector<uint32_t> ReceivedBlockIds;
    mutable std::vector<uint8_t> ReceivedBlockStorage;
    // Cold receive payloads are stored in fixed-width slots.  Mapping packet
    // id directly to its slot keeps duplicate validation O(1); after checkpoint
    // adoption the keys remain as the bounded received-id set while the slot
    // values are no longer dereferenced.  A successful systematic-rich cold
    // solve may instead retain the complete table through the first successful
    // RecoverResult call.
    mutable PacketSlotTable ReceivedSlots;
    PrecodeSolveResumeState ResumeState;
    std::vector<uint8_t> PendingPacketStorage;
    uint32_t PendingPacketId = 0u;
    std::vector<uint8_t> IntermediateBlockStorage;
    std::unique_ptr<uint8_t[]> SystematicPacketCache;
    size_t SystematicPacketCacheSize = 0u;
    std::vector<uint8_t> HaveSystematicPacket;
    uint32_t CachedSystematicPacketCountValue = 0u;
    PrecodeSolveStats SolveStatsValue = {};
    uint64_t MessageBytesValue = 0;
    uint32_t BlockBytesValue = 0;
    // Before ordinary decode this stores the accepted systematic count.  Two
    // private high bits distinguish direct decoded/source-ordered states;
    // public ReceivedCount() masks them.  Ordinary solve completion replaces
    // the value with the accepted packet count, preserving the legacy layout.
    mutable uint32_t ReceivedCountValue = 0;
    uint32_t SolveAttemptCountValue = 0;
    uint32_t PacketSeedAttemptValue = 0;
    WirehairResult LastSolveResult = Wirehair_NeedMore;
    bool PendingPacket = false;
    bool Initialized = false;
    bool Decoded = false;
};

} // namespace wirehair_v2
