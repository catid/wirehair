#include "WirehairV2TinyMDS.h"

#include "../gf256.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wirehair_v2 {
namespace {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local int64_t TinyMdsAllocationFailureCountdown = -1;

void GuardedTinyMdsAllocation()
{
    if (TinyMdsAllocationFailureCountdown == 0) {
        throw std::bad_alloc();
    }
    if (TinyMdsAllocationFailureCountdown > 0) {
        --TinyMdsAllocationFailureCountdown;
    }
}
#else
void GuardedTinyMdsAllocation() {}
#endif

void SetLinearCombination(
    uint8_t* output,
    uint8_t first_scale,
    const uint8_t* first,
    uint8_t second_scale,
    const uint8_t* second,
    uint32_t bytes)
{
    const int kernel_bytes = (int)bytes;
    if (first_scale == 0u) {
        std::memset(output, 0, bytes);
    }
    else if (first_scale == 1u) {
        std::memcpy(output, first, bytes);
    }
    else {
        gf256_mul_mem(output, first, first_scale, kernel_bytes);
    }

    if (second_scale == 0u) {
        return;
    }
    if (second_scale == 1u) {
        gf256_add_mem(output, second, kernel_bytes);
    }
    else {
        gf256_muladd_mem(
            output, second_scale, second, kernel_bytes);
    }
}

} // namespace

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void SetTinyMdsAllocationFailureCountdownForTesting(int64_t countdown)
{
    TinyMdsAllocationFailureCountdown = countdown;
}
#endif

bool TinyMdsCodec::MessageBlockCount(
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint32_t& block_count_out)
{
    block_count_out = 0u;
    if (message_bytes == 0u || block_bytes == 0u ||
        block_bytes > UINT32_C(0x7fffffff))
    {
        return false;
    }
    const uint64_t block_count =
        message_bytes / block_bytes +
        (message_bytes % block_bytes != 0u ? 1u : 0u);
    if (block_count < 1u || block_count > 2u) {
        return false;
    }
    block_count_out = (uint32_t)block_count;
    return true;
}

bool TinyMdsCodec::CoefficientsForK2Packet(
    uint32_t block_id,
    uint8_t& alpha,
    uint8_t& beta)
{
    if (block_id > kMaxK2PacketId) {
        return false;
    }
    if (block_id == 0u)
    {
        alpha = 1u;
        beta = 0u;
        return true;
    }
    if (block_id == 1u)
    {
        alpha = 0u;
        beta = 1u;
        return true;
    }
    alpha = 1u;
    beta = (uint8_t)(block_id - 1u);
    return true;
}

uint32_t TinyMdsCodec::PacketBytes(uint32_t block_id) const
{
    if (block_id == BlockCountValue - 1u)
    {
        return (uint32_t)(MessageBytesValue -
            (uint64_t)(BlockCountValue - 1u) * BlockBytesValue);
    }
    return BlockBytesValue;
}

uint8_t* TinyMdsCodec::SourceBlock(uint32_t block_id)
{
    return Storage.data() + (size_t)block_id * BlockBytesValue;
}

const uint8_t* TinyMdsCodec::SourceBlock(uint32_t block_id) const
{
    return Storage.data() + (size_t)block_id * BlockBytesValue;
}

uint8_t* TinyMdsCodec::FirstPacketScratch()
{
    return Storage.data() + (size_t)BlockCountValue * BlockBytesValue;
}

const uint8_t* TinyMdsCodec::FirstPacketScratch() const
{
    return Storage.data() + (size_t)BlockCountValue * BlockBytesValue;
}

uint8_t* TinyMdsCodec::IncomingPacketScratch()
{
    return FirstPacketScratch() + BlockBytesValue;
}

void TinyMdsCodec::Swap(TinyMdsCodec& other) noexcept
{
    using std::swap;
    Storage.swap(other.Storage);
    swap(MessageBytesValue, other.MessageBytesValue);
    swap(BlockBytesValue, other.BlockBytesValue);
    swap(BlockCountValue, other.BlockCountValue);
    swap(FirstPacketId, other.FirstPacketId);
    swap(FirstPacketBytes, other.FirstPacketBytes);
    swap(ReceivedCountValue, other.ReceivedCountValue);
    swap(LastDecodeResult, other.LastDecodeResult);
    swap(ModeValue, other.ModeValue);
    swap(HaveFirstPacket, other.HaveFirstPacket);
    swap(Decoded, other.Decoded);
}

WirehairResult TinyMdsCodec::InitializeEncoder(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes)
{
    if (!message) {
        return Wirehair_InvalidInput;
    }
    uint32_t block_count = 0u;
    if (!MessageBlockCount(
            message_bytes, block_bytes, block_count))
    {
        return Wirehair_InvalidInput;
    }
    if (message_bytes >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        return Wirehair_UnsupportedPlatform;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }

    try
    {
        if (block_bytes >
            std::numeric_limits<size_t>::max() / block_count)
        {
            return Wirehair_OOM;
        }
        const size_t storage_bytes =
            (size_t)block_count * block_bytes;
        TinyMdsCodec next;
        GuardedTinyMdsAllocation();
        next.Storage.assign(storage_bytes, uint8_t{0});
        std::memcpy(
            next.Storage.data(), message, (size_t)message_bytes);
        next.MessageBytesValue = message_bytes;
        next.BlockBytesValue = block_bytes;
        next.BlockCountValue = block_count;
        next.ModeValue = Mode::Encoder;
        Swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

WirehairResult TinyMdsCodec::InitializeDecoder(
    uint64_t message_bytes,
    uint32_t block_bytes)
{
    uint32_t block_count = 0u;
    if (!MessageBlockCount(
            message_bytes, block_bytes, block_count))
    {
        return Wirehair_InvalidInput;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }

    try
    {
        // Source blocks plus one retained/validation packet.  K=2 needs one
        // additional padded incoming packet while the two sources are solved.
        const size_t storage_blocks =
            (size_t)block_count + 1u + (block_count == 2u ? 1u : 0u);
        if (block_bytes >
            std::numeric_limits<size_t>::max() / storage_blocks)
        {
            return Wirehair_OOM;
        }
        TinyMdsCodec next;
        GuardedTinyMdsAllocation();
        next.Storage.assign(
            storage_blocks * block_bytes, uint8_t{0});
        next.MessageBytesValue = message_bytes;
        next.BlockBytesValue = block_bytes;
        next.BlockCountValue = block_count;
        next.ModeValue = Mode::Decoder;
        Swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

WirehairResult TinyMdsCodec::EncodeResult(
    uint32_t block_id,
    void* block_out,
    uint32_t output_capacity,
    uint32_t* data_bytes_out) const
{
    if (ModeValue != Mode::Encoder || !block_out || !data_bytes_out) {
        return Wirehair_InvalidInput;
    }
    uint8_t alpha = 0u;
    uint8_t beta = 0u;
    if (BlockCountValue == 2u &&
        !CoefficientsForK2Packet(block_id, alpha, beta))
    {
        return Wirehair_InvalidInput;
    }
    const uint32_t data_bytes = PacketBytes(block_id);
    if (output_capacity < data_bytes) {
        return Wirehair_InvalidInput;
    }

    uint8_t* const output = static_cast<uint8_t*>(block_out);
    if (BlockCountValue == 1u) {
        std::memcpy(output, SourceBlock(0u), data_bytes);
    }
    else if (block_id < 2u) {
        std::memcpy(output, SourceBlock(block_id), data_bytes);
    }
    else
    {
        SetLinearCombination(
            output, alpha, SourceBlock(0u),
            beta, SourceBlock(1u), BlockBytesValue);
    }
    *data_bytes_out = data_bytes;
    return Wirehair_Success;
}

WirehairResult TinyMdsCodec::DecodeResult(
    uint32_t block_id,
    const void* block_in,
    uint32_t data_bytes)
{
    if (ModeValue != Mode::Decoder || !block_in) {
        return Wirehair_InvalidInput;
    }
    uint8_t alpha = 0u;
    uint8_t beta = 0u;
    if (BlockCountValue == 2u &&
        !CoefficientsForK2Packet(block_id, alpha, beta))
    {
        return Wirehair_InvalidInput;
    }
    const uint32_t expected_bytes = PacketBytes(block_id);
    if (data_bytes != expected_bytes) {
        return Wirehair_InvalidInput;
    }
    const uint8_t* const input =
        static_cast<const uint8_t*>(block_in);

    if (Decoded)
    {
        uint8_t* const expected = FirstPacketScratch();
        if (BlockCountValue == 1u) {
            std::memcpy(expected, SourceBlock(0u), expected_bytes);
        }
        else if (block_id < 2u) {
            std::memcpy(
                expected, SourceBlock(block_id), expected_bytes);
        }
        else
        {
            SetLinearCombination(
                expected, alpha, SourceBlock(0u),
                beta, SourceBlock(1u), BlockBytesValue);
        }
        return std::memcmp(expected, input, expected_bytes) == 0 ?
            Wirehair_Success : Wirehair_Error;
    }

    if (BlockCountValue == 1u)
    {
        // A repair carries the canonical zero padding too.  Reject a packet
        // that cannot be an encoding of any one-block message.
        if (block_id != 0u)
        {
            for (uint32_t i = (uint32_t)MessageBytesValue;
                 i < BlockBytesValue;
                 ++i)
            {
                if (input[i] != 0u) {
                    return Wirehair_Error;
                }
            }
        }
        std::memcpy(SourceBlock(0u), input, expected_bytes);
        std::fill(
            SourceBlock(0u) + expected_bytes,
            SourceBlock(0u) + BlockBytesValue,
            uint8_t{0});
        ReceivedCountValue = 1u;
        LastDecodeResult = Wirehair_Success;
        Decoded = true;
        return Wirehair_Success;
    }

    if (HaveFirstPacket && block_id == FirstPacketId)
    {
        if (FirstPacketBytes != data_bytes ||
            std::memcmp(
                FirstPacketScratch(), input, data_bytes) != 0)
        {
            return Wirehair_InvalidInput;
        }
        return LastDecodeResult;
    }

    if (!HaveFirstPacket)
    {
        uint8_t* const first = FirstPacketScratch();
        std::memset(first, 0, BlockBytesValue);
        std::memcpy(first, input, data_bytes);
        FirstPacketId = block_id;
        FirstPacketBytes = data_bytes;
        ReceivedCountValue = 1u;
        LastDecodeResult = Wirehair_NeedMore;
        HaveFirstPacket = true;
        return Wirehair_NeedMore;
    }

    uint8_t first_alpha = 0u;
    uint8_t first_beta = 0u;
    if (!CoefficientsForK2Packet(
            FirstPacketId, first_alpha, first_beta))
    {
        return Wirehair_Error;
    }
    uint8_t* const incoming = IncomingPacketScratch();
    std::memset(incoming, 0, BlockBytesValue);
    std::memcpy(incoming, input, data_bytes);

    const uint8_t determinant =
        gf256_mul(first_alpha, beta) ^
        gf256_mul(alpha, first_beta);
    if (determinant == 0u) {
        return Wirehair_Error;
    }
    const uint8_t inverse = gf256_inv(determinant);
    const uint8_t a_first = gf256_mul(beta, inverse);
    const uint8_t a_second = gf256_mul(first_beta, inverse);
    const uint8_t b_first = gf256_mul(alpha, inverse);
    const uint8_t b_second = gf256_mul(first_alpha, inverse);
    SetLinearCombination(
        SourceBlock(0u),
        a_first, FirstPacketScratch(),
        a_second, incoming,
        BlockBytesValue);
    SetLinearCombination(
        SourceBlock(1u),
        b_first, FirstPacketScratch(),
        b_second, incoming,
        BlockBytesValue);

    // Every repair packet is evaluated over the zero-padded final source
    // block.  Two arbitrary projective-line packets always define an
    // algebraic solution, so explicitly reject one that cannot correspond to
    // the canonical encoding of the declared message size.  The retained
    // first packet and all public state remain unchanged, allowing a valid
    // replacement second packet to be retried.
    const uint32_t tail_bytes = PacketBytes(1u);
    for (uint32_t i = tail_bytes; i < BlockBytesValue; ++i)
    {
        if (SourceBlock(1u)[i] != 0u) {
            return Wirehair_Error;
        }
    }

    ReceivedCountValue = 2u;
    LastDecodeResult = Wirehair_Success;
    Decoded = true;
    return Wirehair_Success;
}

WirehairResult TinyMdsCodec::RecoverResult(
    void* message_out,
    uint64_t message_bytes) const
{
    if (ModeValue != Mode::Decoder || !message_out ||
        message_bytes != MessageBytesValue)
    {
        return Wirehair_InvalidInput;
    }
    if (!Decoded) {
        return Wirehair_NeedMore;
    }
    std::memcpy(
        message_out, Storage.data(), (size_t)MessageBytesValue);
    return Wirehair_Success;
}

bool TinyMdsCodec::IsEncoder() const
{
    return ModeValue == Mode::Encoder;
}

bool TinyMdsCodec::IsDecoder() const
{
    return ModeValue == Mode::Decoder;
}

bool TinyMdsCodec::IsDecoded() const
{
    return IsDecoder() && Decoded;
}

uint32_t TinyMdsCodec::BlockCount() const
{
    return ModeValue == Mode::Uninitialized ? 0u : BlockCountValue;
}

uint32_t TinyMdsCodec::BlockBytes() const
{
    return ModeValue == Mode::Uninitialized ? 0u : BlockBytesValue;
}

uint64_t TinyMdsCodec::MessageBytes() const
{
    return ModeValue == Mode::Uninitialized ? 0u : MessageBytesValue;
}

uint32_t TinyMdsCodec::ReceivedCount() const
{
    return IsDecoder() ? ReceivedCountValue : 0u;
}

} // namespace wirehair_v2
