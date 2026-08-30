#include "Wh2K2ProjectiveCodec.h"

#include "../gf256.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wirehair_wh2_bench {
namespace {

static const uint32_t kHighIdSalt = UINT32_C(0x00000000);

#if defined(WIREHAIR_WH2_K2_PROJECTIVE_TEST_HOOKS)
thread_local int64_t K2ProjectiveAllocationFailureCountdown = -1;

void GuardedK2ProjectiveAllocation()
{
    if (K2ProjectiveAllocationFailureCountdown == 0) {
        throw std::bad_alloc();
    }
    if (K2ProjectiveAllocationFailureCountdown > 0) {
        --K2ProjectiveAllocationFailureCountdown;
    }
}
#else
inline void GuardedK2ProjectiveAllocation() noexcept
{
}
#endif

uint32_t Avalanche32(uint32_t value) noexcept
{
    // Each xor-shift is invertible, and both multipliers are odd, so this is
    // a permutation of all 2^32 inputs rather than a truncating hash.
    value = (value ^ (value >> 16)) * UINT32_C(0x7feb352d);
    value = (value ^ (value >> 15)) * UINT32_C(0x846ca68b);
    return value ^ (value >> 16);
}

uint8_t BalancedNonzeroCoefficient(uint32_t value) noexcept
{
    // Multiply-high maps a uniform uint32 into 0..254 with bucket sizes that
    // differ by at most one.  Adding one excludes the two systematic points
    // [1:0] and [0:1] from the high-ID repair range.
    const uint32_t reduced = static_cast<uint32_t>(
        (static_cast<uint64_t>(value) * UINT64_C(255)) >> 32);
    return static_cast<uint8_t>(reduced + 1u);
}

bool CheckedProduct(
    size_t left,
    size_t right,
    size_t& product_out) noexcept
{
    product_out = 0u;
    if (left != 0u &&
        right > std::numeric_limits<size_t>::max() / left)
    {
        return false;
    }
    product_out = left * right;
    return true;
}

bool SameDirection(
    const K2ProjectiveDirection& left,
    uint16_t right_point) noexcept
{
    return left.Point == right_point;
}

void SetLinearCombination(
    uint8_t* output,
    uint8_t first_scale,
    const uint8_t* first,
    uint8_t second_scale,
    const uint8_t* second,
    uint32_t bytes) noexcept
{
    const int kernel_bytes = static_cast<int>(bytes);
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

K2ProjectiveDirection K2ProjectiveDirectionForPacketId(
    uint32_t block_id) noexcept
{
    K2ProjectiveDirection direction;
    if (block_id == 0u)
    {
        direction.Point = 0u;
        direction.Alpha = 1u;
        direction.Beta = 0u;
        return direction;
    }
    if (block_id == 1u)
    {
        direction.Point = 1u;
        direction.Alpha = 0u;
        direction.Beta = 1u;
        return direction;
    }

    uint8_t beta = 0u;
    if (block_id <= 256u) {
        beta = static_cast<uint8_t>(block_id - 1u);
    }
    else {
        beta = BalancedNonzeroCoefficient(
            Avalanche32(block_id ^ kHighIdSalt));
    }
    direction.Point = static_cast<uint16_t>(beta) + 1u;
    direction.Alpha = 1u;
    direction.Beta = beta;
    return direction;
}

uint32_t K2ProjectiveHighIdSalt() noexcept
{
    return kHighIdSalt;
}

#if defined(WIREHAIR_WH2_K2_PROJECTIVE_TEST_HOOKS)
void SetK2ProjectiveAllocationFailureCountdownForTesting(
    int64_t countdown) noexcept
{
    K2ProjectiveAllocationFailureCountdown = countdown;
}
#endif

bool K2ProjectiveCodec::ValidK2Shape(
    uint64_t message_bytes,
    uint32_t block_bytes) noexcept
{
    if (message_bytes == 0u || block_bytes == 0u ||
        block_bytes > UINT32_C(0x7fffffff) ||
        message_bytes >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
    {
        return false;
    }
    const uint64_t block_count =
        message_bytes / block_bytes +
        (message_bytes % block_bytes != 0u ? 1u : 0u);
    return block_count == 2u;
}

void K2ProjectiveCodec::Swap(K2ProjectiveCodec& other) noexcept
{
    using std::swap;
    Storage.swap(other.Storage);
    AcceptedIds.swap(other.AcceptedIds);
    swap(MessageBytesValue, other.MessageBytesValue);
    swap(BlockBytesValue, other.BlockBytesValue);
    swap(FirstPoint, other.FirstPoint);
    swap(FirstAlpha, other.FirstAlpha);
    swap(FirstBeta, other.FirstBeta);
    swap(RankValue, other.RankValue);
    swap(ModeValue, other.ModeValue);
    swap(Decoded, other.Decoded);
}

WirehairResult K2ProjectiveCodec::InitializeEncoder(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes) noexcept
{
    if (!message || !ValidK2Shape(message_bytes, block_bytes)) {
        return Wirehair_InvalidInput;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }

    try
    {
        size_t storage_bytes = 0u;
        if (!CheckedProduct(
                static_cast<size_t>(block_bytes), 2u, storage_bytes)) {
            return Wirehair_OOM;
        }
        K2ProjectiveCodec next;
        GuardedK2ProjectiveAllocation();
        next.Storage.assign(storage_bytes, uint8_t{0});
        std::memcpy(
            next.Storage.data(), message, static_cast<size_t>(message_bytes));
        next.MessageBytesValue = message_bytes;
        next.BlockBytesValue = block_bytes;
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
    catch (...) {
        return Wirehair_Error;
    }
}

WirehairResult K2ProjectiveCodec::
InitializeEncoderOwnedSourceAfterGlobalInit(
    std::vector<uint8_t>&& owned_source,
    uint32_t block_bytes) noexcept
{
    const uint64_t message_bytes = owned_source.size();
    if (!ValidK2Shape(message_bytes, block_bytes)) {
        return Wirehair_InvalidInput;
    }

    try
    {
        size_t storage_bytes = 0u;
        if (!CheckedProduct(
                static_cast<size_t>(block_bytes), 2u, storage_bytes))
        {
            return Wirehair_OOM;
        }
        if (owned_source.size() < storage_bytes)
        {
            GuardedK2ProjectiveAllocation();
            owned_source.resize(storage_bytes, uint8_t{0});
        }
        K2ProjectiveCodec next;
        next.Storage.swap(owned_source);
        next.MessageBytesValue = message_bytes;
        next.BlockBytesValue = block_bytes;
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
    catch (...) {
        return Wirehair_Error;
    }
}

WirehairResult K2ProjectiveCodec::InitializeDecoder(
    uint64_t message_bytes,
    uint32_t block_bytes) noexcept
{
    if (!ValidK2Shape(message_bytes, block_bytes)) {
        return Wirehair_InvalidInput;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }

    try
    {
        size_t storage_bytes = 0u;
        if (!CheckedProduct(
                static_cast<size_t>(block_bytes), 4u, storage_bytes)) {
            return Wirehair_OOM;
        }
        K2ProjectiveCodec next;
        GuardedK2ProjectiveAllocation();
        next.Storage.assign(storage_bytes, uint8_t{0});
        GuardedK2ProjectiveAllocation();
        next.AcceptedIds.reserve(kMaxAcceptedPacketIds);
        next.MessageBytesValue = message_bytes;
        next.BlockBytesValue = block_bytes;
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
    catch (...) {
        return Wirehair_Error;
    }
}

uint32_t K2ProjectiveCodec::PacketBytes(uint32_t block_id) const noexcept
{
    if (block_id == 1u) {
        return static_cast<uint32_t>(
            MessageBytesValue - BlockBytesValue);
    }
    return BlockBytesValue;
}

uint8_t* K2ProjectiveCodec::SourceBlock(uint32_t block_id) noexcept
{
    return Storage.data() + static_cast<size_t>(block_id) * BlockBytesValue;
}

const uint8_t* K2ProjectiveCodec::SourceBlock(
    uint32_t block_id) const noexcept
{
    return Storage.data() + static_cast<size_t>(block_id) * BlockBytesValue;
}

uint8_t* K2ProjectiveCodec::FirstPacketScratch() noexcept
{
    return Storage.data() + static_cast<size_t>(2u) * BlockBytesValue;
}

const uint8_t* K2ProjectiveCodec::FirstPacketScratch() const noexcept
{
    return Storage.data() + static_cast<size_t>(2u) * BlockBytesValue;
}

uint8_t* K2ProjectiveCodec::IncomingPacketScratch() noexcept
{
    return Storage.data() + static_cast<size_t>(3u) * BlockBytesValue;
}

const uint8_t* K2ProjectiveCodec::IncomingPacketScratch() const noexcept
{
    return Storage.data() + static_cast<size_t>(3u) * BlockBytesValue;
}

bool K2ProjectiveCodec::HasAcceptedId(uint32_t block_id) const noexcept
{
    return std::find(
        AcceptedIds.begin(), AcceptedIds.end(), block_id) !=
        AcceptedIds.end();
}

void K2ProjectiveCodec::CanonicalizePacket(
    const void* block_in,
    uint32_t data_bytes,
    uint8_t* canonical_out) const noexcept
{
    std::memset(canonical_out, 0, BlockBytesValue);
    std::memcpy(canonical_out, block_in, data_bytes);
}

void K2ProjectiveCodec::EvaluatePacket(
    const K2ProjectiveDirection& direction,
    uint8_t* canonical_out) const noexcept
{
    SetLinearCombination(
        canonical_out,
        direction.Alpha,
        SourceBlock(0u),
        direction.Beta,
        SourceBlock(1u),
        BlockBytesValue);
}

WirehairResult K2ProjectiveCodec::EncodeResult(
    uint32_t block_id,
    void* block_out,
    uint32_t output_capacity,
    uint32_t* data_bytes_out) const noexcept
{
    if (ModeValue != Mode::Encoder || !block_out || !data_bytes_out) {
        return Wirehair_InvalidInput;
    }
    const uint32_t data_bytes = PacketBytes(block_id);
    if (output_capacity < data_bytes) {
        return Wirehair_InvalidInput;
    }

    uint8_t* const output = static_cast<uint8_t*>(block_out);
    if (block_id < 2u) {
        std::memcpy(output, SourceBlock(block_id), data_bytes);
    }
    else {
        EvaluatePacket(
            K2ProjectiveDirectionForPacketId(block_id), output);
    }
    *data_bytes_out = data_bytes;
    return Wirehair_Success;
}

WirehairResult K2ProjectiveCodec::DecodeResult(
    uint32_t block_id,
    const void* block_in,
    uint32_t data_bytes) noexcept
{
    if (ModeValue != Mode::Decoder || !block_in) {
        return Wirehair_InvalidInput;
    }
    const uint32_t expected_bytes = PacketBytes(block_id);
    if (data_bytes != expected_bytes) {
        return Wirehair_InvalidInput;
    }

    try
    {
        const bool accepted_id = HasAcceptedId(block_id);
        const K2ProjectiveDirection direction =
            K2ProjectiveDirectionForPacketId(block_id);

        uint8_t* const incoming = IncomingPacketScratch();
        CanonicalizePacket(block_in, data_bytes, incoming);

        if (Decoded)
        {
            uint8_t* const expected_packet = FirstPacketScratch();
            EvaluatePacket(direction, expected_packet);
            if (std::memcmp(
                    expected_packet, incoming, BlockBytesValue) == 0)
            {
                if (accepted_id) {
                    return Wirehair_Success;
                }
                if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
                    return Wirehair_ExtraInsufficient;
                }
                // The complete accepted-ID capacity was reserved during
                // initialization, so publishing this newly validated ID is
                // allocation-free and transactional.
                AcceptedIds.push_back(block_id);
                return Wirehair_Success;
            }
            return accepted_id ? Wirehair_InvalidInput : Wirehair_Error;
        }

        if (RankValue == 0u)
        {
            // Decoder initialization reserved the complete accepted-ID cap,
            // so this append cannot allocate or expose a partial state.
            AcceptedIds.push_back(block_id);
            std::memcpy(
                FirstPacketScratch(), incoming, BlockBytesValue);
            FirstPoint = direction.Point;
            FirstAlpha = direction.Alpha;
            FirstBeta = direction.Beta;
            RankValue = 1u;
            return Wirehair_NeedMore;
        }

        if (accepted_id)
        {
            return std::memcmp(
                       FirstPacketScratch(), incoming, BlockBytesValue) == 0 ?
                Wirehair_NeedMore : Wirehair_InvalidInput;
        }

        if (SameDirection(direction, FirstPoint))
        {
            if (std::memcmp(
                    FirstPacketScratch(), incoming, BlockBytesValue) != 0)
            {
                // The new raw ID claims the retained equation direction but
                // carries a different canonical RHS.  Do not publish the ID.
                return Wirehair_Error;
            }
            if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
                return Wirehair_ExtraInsufficient;
            }
            AcceptedIds.push_back(block_id);
            return Wirehair_NeedMore;
        }

        if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
            return Wirehair_ExtraInsufficient;
        }

        const uint8_t determinant =
            gf256_mul(FirstAlpha, direction.Beta) ^
            gf256_mul(direction.Alpha, FirstBeta);
        if (determinant == 0u) {
            return Wirehair_Error;
        }
        const uint8_t inverse = gf256_inv(determinant);
        const uint8_t a_first = gf256_mul(direction.Beta, inverse);
        const uint8_t a_second = gf256_mul(FirstBeta, inverse);
        const uint8_t b_first = gf256_mul(direction.Alpha, inverse);
        const uint8_t b_second = gf256_mul(FirstAlpha, inverse);
        SetLinearCombination(
            SourceBlock(0u),
            a_first,
            FirstPacketScratch(),
            a_second,
            incoming,
            BlockBytesValue);
        SetLinearCombination(
            SourceBlock(1u),
            b_first,
            FirstPacketScratch(),
            b_second,
            incoming,
            BlockBytesValue);

        const uint32_t tail_bytes = static_cast<uint32_t>(
            MessageBytesValue - BlockBytesValue);
        for (uint32_t i = tail_bytes; i < BlockBytesValue; ++i)
        {
            if (SourceBlock(1u)[i] != 0u) {
                return Wirehair_Error;
            }
        }

        AcceptedIds.push_back(block_id);
        RankValue = 2u;
        Decoded = true;
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
    catch (...) {
        return Wirehair_Error;
    }
}

WirehairResult K2ProjectiveCodec::RecoverResult(
    void* message_out,
    uint64_t message_bytes) const noexcept
{
    if (ModeValue != Mode::Decoder || !message_out ||
        message_bytes != MessageBytesValue)
    {
        return Wirehair_InvalidInput;
    }
    if (!Decoded) {
        return Wirehair_NeedMore;
    }
    std::memmove(
        message_out, Storage.data(), static_cast<size_t>(MessageBytesValue));
    return Wirehair_Success;
}

bool K2ProjectiveCodec::IsEncoder() const noexcept
{
    return ModeValue == Mode::Encoder;
}

bool K2ProjectiveCodec::IsDecoder() const noexcept
{
    return ModeValue == Mode::Decoder;
}

bool K2ProjectiveCodec::IsDecoded() const noexcept
{
    return IsDecoder() && Decoded;
}

uint32_t K2ProjectiveCodec::BlockCount() const noexcept
{
    return ModeValue == Mode::Uninitialized ? 0u : 2u;
}

uint32_t K2ProjectiveCodec::BlockBytes() const noexcept
{
    return ModeValue == Mode::Uninitialized ? 0u : BlockBytesValue;
}

uint64_t K2ProjectiveCodec::MessageBytes() const noexcept
{
    return ModeValue == Mode::Uninitialized ? 0u : MessageBytesValue;
}

uint32_t K2ProjectiveCodec::AcceptedIdCount() const noexcept
{
    return IsDecoder() ? static_cast<uint32_t>(AcceptedIds.size()) : 0u;
}

uint32_t K2ProjectiveCodec::Rank() const noexcept
{
    return IsDecoder() ? RankValue : 0u;
}

} // namespace wirehair_wh2_bench
