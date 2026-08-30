#pragma once

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <vector>

namespace wirehair_wh2_bench {

/**
    One normalized point of P^1(GF(256)).

    Point 0 is [1:0], point 1 is [0:1], and points 2..256 are
    [1:point-1].  Alpha and Beta repeat that normalized coefficient pair so
    tests do not need to duplicate the point convention.
*/
struct K2ProjectiveDirection
{
    uint16_t Point = 0u;
    uint8_t Alpha = 0u;
    uint8_t Beta = 0u;
};

/**
    Frozen all-uint32 packet-id mapping for the bounded K=2 experiment.

    IDs 0..256 enumerate the projective line exactly.  Higher IDs use the
    frozen bijective 32-bit avalanche and balanced reduction documented in
    Wh2K2ProjectiveCodec.cpp.  Higher IDs deliberately select only points
    2..256, so every high-ID repair is independent of either systematic row.
*/
K2ProjectiveDirection K2ProjectiveDirectionForPacketId(
    uint32_t block_id) noexcept;

/** Frozen high-ID salt selected by the bounded K=2 trace scan. */
uint32_t K2ProjectiveHighIdSalt() noexcept;

/**
    Test/benchmark-only K=2 systematic near-MDS codec over GF(256).

    The first 257 packet IDs are exactly MDS.  All uint32 packet IDs are
    accepted, but higher IDs necessarily reuse projective directions because
    P^1(GF(256)) contains only 257 points.  The decoder therefore tracks raw
    packet IDs and projective rank separately.  A matching new-ID alias is a
    redundant equation; an inconsistent alias fails without changing the
    published decoder state.
*/
class K2ProjectiveCodec
{
public:
    K2ProjectiveCodec() noexcept = default;
    K2ProjectiveCodec(const K2ProjectiveCodec&) = delete;
    K2ProjectiveCodec& operator=(const K2ProjectiveCodec&) = delete;

    WirehairResult InitializeEncoder(
        const void* message,
        uint64_t message_bytes,
        uint32_t block_bytes) noexcept;

    /**
        Benchmark timing seam matching NativeEncoderFixture source ownership.
        The caller must complete wirehair_init() before the clock and pass its
        fresh owned source buffer by move.  No global initialization or
        second source allocation/copy occurs for a full K*B source.
    */
    WirehairResult InitializeEncoderOwnedSourceAfterGlobalInit(
        std::vector<uint8_t>&& owned_source,
        uint32_t block_bytes) noexcept;

    WirehairResult InitializeDecoder(
        uint64_t message_bytes,
        uint32_t block_bytes) noexcept;

    /** On failure, block_out and data_bytes_out are unchanged. */
    WirehairResult EncodeResult(
        uint32_t block_id,
        void* block_out,
        uint32_t output_capacity,
        uint32_t* data_bytes_out) const noexcept;

    /**
        Raw-ID duplicates are idempotent and conflicting duplicates are
        InvalidInput.  A matching different-ID projective alias is redundant
        and returns NeedMore; an inconsistent alias returns Error.  Neither
        alias case changes algebraic rank.
    */
    WirehairResult DecodeResult(
        uint32_t block_id,
        const void* block_in,
        uint32_t data_bytes) noexcept;

    WirehairResult RecoverResult(
        void* message_out,
        uint64_t message_bytes) const noexcept;

    bool IsEncoder() const noexcept;
    bool IsDecoder() const noexcept;
    bool IsDecoded() const noexcept;
    uint32_t BlockCount() const noexcept;
    uint32_t BlockBytes() const noexcept;
    uint64_t MessageBytes() const noexcept;
    uint32_t AcceptedIdCount() const noexcept;
    uint32_t Rank() const noexcept;

    static const uint32_t kMaxAcceptedPacketIds = 2u + 1024u;

private:
    enum class Mode : uint8_t
    {
        Uninitialized,
        Encoder,
        Decoder
    };

    void Swap(K2ProjectiveCodec& other) noexcept;

    static bool ValidK2Shape(
        uint64_t message_bytes,
        uint32_t block_bytes) noexcept;

    uint32_t PacketBytes(uint32_t block_id) const noexcept;
    uint8_t* SourceBlock(uint32_t block_id) noexcept;
    const uint8_t* SourceBlock(uint32_t block_id) const noexcept;
    uint8_t* FirstPacketScratch() noexcept;
    const uint8_t* FirstPacketScratch() const noexcept;
    uint8_t* IncomingPacketScratch() noexcept;
    const uint8_t* IncomingPacketScratch() const noexcept;
    bool HasAcceptedId(uint32_t block_id) const noexcept;
    void CanonicalizePacket(
        const void* block_in,
        uint32_t data_bytes,
        uint8_t* canonical_out) const noexcept;
    void EvaluatePacket(
        const K2ProjectiveDirection& direction,
        uint8_t* canonical_out) const noexcept;

    std::vector<uint8_t> Storage;
    std::vector<uint32_t> AcceptedIds;
    uint64_t MessageBytesValue = 0u;
    uint32_t BlockBytesValue = 0u;
    uint16_t FirstPoint = 0u;
    uint8_t FirstAlpha = 0u;
    uint8_t FirstBeta = 0u;
    uint8_t RankValue = 0u;
    Mode ModeValue = Mode::Uninitialized;
    bool Decoded = false;
};

#if defined(WIREHAIR_WH2_K2_PROJECTIVE_TEST_HOOKS)
/** Fail the guarded allocation at countdown zero; negative disables it. */
void SetK2ProjectiveAllocationFailureCountdownForTesting(
    int64_t countdown) noexcept;
#endif

} // namespace wirehair_wh2_bench
