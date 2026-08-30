#pragma once

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <vector>

namespace wirehair_wh2_bench {

/**
    One normalized point of PG(2, GF(256)).

    Point values 0..257 identify the frozen regular hyperoval.  Values
    258..65792 identify the exact 65,535-point complement in bucket order.
    Every direction is normalized so its first nonzero coefficient is one.
*/
struct K3HyperovalDirection
{
    uint32_t Point = 0u;
    uint8_t Alpha = 0u;
    uint8_t Beta = 0u;
    uint8_t Gamma = 0u;
};

/**
    Map one exact complement bucket in [0,65534].  Invalid buckets leave
    direction_out unchanged and return false.
*/
bool K3HyperovalComplementDirectionForBucket(
    uint32_t bucket,
    K3HyperovalDirection& direction_out) noexcept;

/**
    Frozen all-uint32 packet-id mapping for the benchmark-only K=3 lane.

    IDs 0..257 enumerate the regular hyperoval.  IDs 258..UINT32_MAX apply the
    frozen bijective avalanche and balanced projection into its complement.
*/
K3HyperovalDirection K3HyperovalDirectionForPacketId(
    uint32_t block_id) noexcept;

/** Caller must first complete gf256_init() or wirehair_init() successfully. */
uint8_t K3HyperovalDeterminant(
    const K3HyperovalDirection& first,
    const K3HyperovalDirection& second,
    const K3HyperovalDirection& third) noexcept;

/**
    Cross-check the allocation-free local polynomial-0x14d square mapper
    against the initialized repository GF(256) primitive.  The caller must
    first complete gf256_init() or wirehair_init() successfully.
*/
bool K3HyperovalSquareMapperMatchesInitializedGf256() noexcept;

/** Conceptual full-block operations, independent of WH_COUNT special cases. */
struct K3HyperovalBlockOps
{
    uint64_t Copies = 0u;
    uint64_t Xors = 0u;
    uint64_t Scales = 0u;
    uint64_t MulAdds = 0u;
};

struct K3HyperovalWorkCounters
{
    K3HyperovalBlockOps PacketEvaluation;
    K3HyperovalBlockOps DependentValidation;
    K3HyperovalBlockOps Solve;
};

/**
    Benchmark-only systematic K=3 codec over GF(256).

    Raw packet IDs and algebraic rank are tracked independently.  Before
    success, a matching dependent equation is retained as an accepted raw ID
    and returns NeedMore; an inconsistent new dependent equation returns Error.
    A conflicting retained raw ID returns InvalidInput.  After success, every
    equation is checked against the recovered source before its raw ID can be
    published.
*/
class K3HyperovalCodec
{
public:
    K3HyperovalCodec() noexcept = default;
    K3HyperovalCodec(const K3HyperovalCodec&) = delete;
    K3HyperovalCodec& operator=(const K3HyperovalCodec&) = delete;

    WirehairResult InitializeEncoder(
        const void* message,
        uint64_t message_bytes,
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

    K3HyperovalWorkCounters WorkCounters() const noexcept;
    void ResetWorkCounters() noexcept;

    static const uint32_t kMaxAcceptedPacketIds = 3u + 1024u;

private:
    enum class Mode : uint8_t
    {
        Uninitialized,
        Encoder,
        Decoder
    };

    void Swap(K3HyperovalCodec& other) noexcept;

    static bool ValidK3Shape(
        uint64_t message_bytes,
        uint32_t block_bytes) noexcept;

    uint32_t PacketBytes(uint32_t block_id) const noexcept;
    uint8_t* SourceBlock(uint32_t block_id) noexcept;
    const uint8_t* SourceBlock(uint32_t block_id) const noexcept;
    uint8_t* RetainedPacket(uint32_t row) noexcept;
    const uint8_t* RetainedPacket(uint32_t row) const noexcept;
    uint8_t* IncomingPacketScratch() noexcept;
    const uint8_t* IncomingPacketScratch() const noexcept;
    bool HasAcceptedId(uint32_t block_id) const noexcept;

    void CanonicalizePacket(
        const void* block_in,
        uint32_t data_bytes,
        uint8_t* canonical_out) const noexcept;

    bool EvaluatePacket(
        const K3HyperovalDirection& direction,
        uint8_t* canonical_out,
        K3HyperovalBlockOps& counters) const noexcept;

    bool RetainedSpanCoefficients(
        const K3HyperovalDirection& direction,
        uint8_t& first_scale,
        uint8_t& second_scale) const noexcept;

    bool SolveThreeRows(
        const K3HyperovalDirection& third_direction,
        uint8_t determinant) noexcept;

    std::vector<uint8_t> Storage;
    std::vector<uint32_t> AcceptedIds;
    uint64_t MessageBytesValue = 0u;
    uint32_t BlockBytesValue = 0u;
    K3HyperovalDirection FirstDirection;
    K3HyperovalDirection SecondDirection;
    mutable K3HyperovalWorkCounters WorkCountersValue;
    uint8_t RankValue = 0u;
    Mode ModeValue = Mode::Uninitialized;
    bool Decoded = false;
};

#if defined(WIREHAIR_WH2_K3_HYPEROVAL_TEST_HOOKS)
/** Fail the guarded allocation at countdown zero; negative disables it. */
void SetK3HyperovalAllocationFailureCountdownForTesting(
    int64_t countdown) noexcept;
#endif

} // namespace wirehair_wh2_bench
