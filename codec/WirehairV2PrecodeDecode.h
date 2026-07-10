#pragma once

#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <unordered_set>
#include <vector>

namespace wirehair_v2 {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void SetDecoderAllocationFailureCountdownForTesting(int64_t countdown);
#endif

/**
    Incremental message decoder for the version-4 packet/precode path.

    NeedMore preserves received packets so later calls can resume.  The current
    prototype rebuilds the sparse solve when an extra packet follows a rank
    failure; precodefail/compare report that real cost rather than claiming a
    persistent factorization.  Packet authentication is a caller concern:
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

private:
    WirehairResult AttemptSolve();
    void Swap(MessagePrecodeDecoder& other) noexcept;

    SeedProfile ProfileValue = {};
    MessagePrecodeEncoderOptions OptionsValue = {};
    PacketRowConfig PacketConfigValue = {};
    PrecodeSystem SystemValue = {};
    std::vector<uint32_t> ReceivedBlockIds;
    std::vector<uint8_t> ReceivedBlockStorage;
    std::unordered_set<uint32_t> ReceivedIds;
    std::vector<uint8_t> IntermediateBlockStorage;
    PrecodeSolveStats SolveStatsValue = {};
    uint64_t MessageBytesValue = 0;
    uint32_t BlockBytesValue = 0;
    uint32_t ReceivedCountValue = 0;
    uint32_t SolveAttemptCountValue = 0;
    uint32_t PacketSeedAttemptValue = 0;
    WirehairResult LastSolveResult = Wirehair_NeedMore;
    bool Initialized = false;
    bool Decoded = false;
};

} // namespace wirehair_v2
