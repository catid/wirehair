#pragma once

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <vector>

namespace wirehair_v2 {

/**
    Exact K=1/K=2 MDS codec over GF(256).

    For K=2, packet ids 0..256 enumerate P^1(GF(256)):

        0 -> [1:0], 1 -> [0:1], i -> [1:i-1] for i=2..256.

    Any two distinct supported ids are independent.  K=1 is the degenerate
    repetition code and accepts every uint32_t packet id.  Source tails are
    canonically zero padded before repair evaluation.
*/
class TinyMdsCodec
{
public:
    TinyMdsCodec() = default;
    TinyMdsCodec(const TinyMdsCodec&) = delete;
    TinyMdsCodec& operator=(const TinyMdsCodec&) = delete;

    WirehairResult InitializeEncoder(
        const void* message,
        uint64_t message_bytes,
        uint32_t block_bytes);

    WirehairResult InitializeDecoder(
        uint64_t message_bytes,
        uint32_t block_bytes);

    /**
        On failure, block_out and data_bytes_out are unchanged.
    */
    WirehairResult EncodeResult(
        uint32_t block_id,
        void* block_out,
        uint32_t output_capacity,
        uint32_t* data_bytes_out) const;

    /**
        Matching duplicates retain the prior result.  A conflicting duplicate
        before completion is InvalidInput, matching MessagePrecodeDecoder.
    */
    WirehairResult DecodeResult(
        uint32_t block_id,
        const void* block_in,
        uint32_t data_bytes);

    /**
        Requires the exact initialized message size and a completed decoder.
    */
    WirehairResult RecoverResult(
        void* message_out,
        uint64_t message_bytes) const;

    bool IsEncoder() const;
    bool IsDecoder() const;
    bool IsDecoded() const;
    uint32_t BlockCount() const;
    uint32_t BlockBytes() const;
    uint64_t MessageBytes() const;
    uint32_t ReceivedCount() const;

    static const uint32_t kMaxK2PacketId = 256u;

private:
    enum class Mode : uint8_t
    {
        Uninitialized,
        Encoder,
        Decoder
    };

    void Swap(TinyMdsCodec& other) noexcept;

    static bool MessageBlockCount(
        uint64_t message_bytes,
        uint32_t block_bytes,
        uint32_t& block_count_out);
    static bool CoefficientsForK2Packet(
        uint32_t block_id,
        uint8_t& alpha,
        uint8_t& beta);
    uint32_t PacketBytes(uint32_t block_id) const;
    uint8_t* SourceBlock(uint32_t block_id);
    const uint8_t* SourceBlock(uint32_t block_id) const;
    uint8_t* FirstPacketScratch();
    const uint8_t* FirstPacketScratch() const;
    uint8_t* IncomingPacketScratch();

    std::vector<uint8_t> Storage;
    uint64_t MessageBytesValue = 0u;
    uint32_t BlockBytesValue = 0u;
    uint32_t BlockCountValue = 0u;
    uint32_t FirstPacketId = 0u;
    uint32_t FirstPacketBytes = 0u;
    uint32_t ReceivedCountValue = 0u;
    WirehairResult LastDecodeResult = Wirehair_NeedMore;
    Mode ModeValue = Mode::Uninitialized;
    bool HaveFirstPacket = false;
    bool Decoded = false;
};

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/// Fail the next guarded TinyMdsCodec storage allocation at countdown zero.
void SetTinyMdsAllocationFailureCountdownForTesting(int64_t countdown);
#endif

} // namespace wirehair_v2
