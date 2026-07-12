#pragma once

#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"

#include <wirehair/wirehair.h>

#include <stddef.h>
#include <stdint.h>
#include <memory>
#include <vector>

namespace wirehair_v2 {

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
    bool Insert(uint32_t packet_id, uint32_t slot);
    bool Erase(uint32_t packet_id);
    void ClearAndRelease() noexcept;
    void Swap(PacketSlotTable& other) noexcept;

    size_t Size() const { return EntryCount; }
    size_t Capacity() const { return Keys.size(); }
    size_t StorageBytes() const;

private:
    static uint32_t Hash(uint32_t packet_id);
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

    /** Identical duplicate ids are ignored; conflicting duplicates fail. */
    WirehairResult DecodeResult(
        uint32_t block_id,
        const void* block_in,
        uint32_t data_bytes);

    /** Recover requires the exact initialized message size. */
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
    const uint8_t* IntermediateBlocks() const;

    /** Release the optional received-systematic cache; idempotent. */
    void ReleaseSystematicPacketCache() noexcept;
    bool HasSystematicPacketCache() const;
    size_t SystematicPacketCacheBytes() const;
    uint32_t CachedSystematicPacketCount() const;

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    bool HasIncrementalResumeStateForTesting() const;
    size_t IncrementalResumeBytesForTesting() const;
    size_t ColdReceiveCapacityBytesForTesting() const;
#endif

private:
    WirehairResult AttemptSolve();
    void Swap(MessagePrecodeDecoder& other) noexcept;

    SeedProfile ProfileValue = {};
    MessagePrecodeEncoderOptions OptionsValue = {};
    PacketRowConfig PacketConfigValue = {};
    PacketRowRuntime PacketRuntimeValue = {};
    PrecodeSystem SystemValue = {};
    std::vector<uint32_t> ReceivedBlockIds;
    std::vector<uint8_t> ReceivedBlockStorage;
    // Cold receive payloads are stored in fixed-width slots.  Mapping packet
    // id directly to its slot keeps duplicate validation O(1); after checkpoint
    // adoption the keys remain as the bounded received-id set while the slot
    // values are no longer dereferenced.
    PacketSlotTable ReceivedSlots;
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
    uint32_t ReceivedCountValue = 0;
    uint32_t SolveAttemptCountValue = 0;
    uint32_t PacketSeedAttemptValue = 0;
    WirehairResult LastSolveResult = Wirehair_NeedMore;
    bool PendingPacket = false;
    bool Initialized = false;
    bool Decoded = false;
};

} // namespace wirehair_v2
