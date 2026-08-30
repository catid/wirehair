#include "Wh2K3HyperovalCodec.h"

#include "../gf256.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wirehair_wh2_bench {
namespace {

static const uint32_t kLowPointCount = 258u;
static const uint32_t kAffineComplementCount = 65280u;
static const uint32_t kComplementPointCount = 65535u;

#if defined(WIREHAIR_WH2_K3_HYPEROVAL_TEST_HOOKS)
thread_local int64_t K3HyperovalAllocationFailureCountdown = -1;

void GuardedK3HyperovalAllocation()
{
    if (K3HyperovalAllocationFailureCountdown == 0) {
        throw std::bad_alloc();
    }
    if (K3HyperovalAllocationFailureCountdown > 0) {
        --K3HyperovalAllocationFailureCountdown;
    }
}
#else
inline void GuardedK3HyperovalAllocation() noexcept
{
}
#endif

uint8_t Polynomial14dMultiply(uint8_t left, uint8_t right) noexcept
{
    uint8_t product = 0u;
    for (unsigned bit = 0u; bit < 8u; ++bit)
    {
        if ((right & 1u) != 0u) {
            product ^= left;
        }
        const bool high_bit = (left & 0x80u) != 0u;
        left = static_cast<uint8_t>(left << 1u);
        if (high_bit) {
            left ^= UINT8_C(0x4d);
        }
        right = static_cast<uint8_t>(right >> 1u);
    }
    return product;
}

struct Polynomial14dSquareTable
{
    Polynomial14dSquareTable() noexcept
    {
        for (uint32_t value = 0u; value < 256u; ++value)
        {
            const uint8_t byte = static_cast<uint8_t>(value);
            Values[value] = Polynomial14dMultiply(byte, byte);
        }
    }

    uint8_t Values[256];
};

uint8_t Polynomial14dSquare(uint8_t value) noexcept
{
    static const Polynomial14dSquareTable table;
    return table.Values[value];
}

bool EnsureSquareMapperMatchesInitializedGf256() noexcept
{
    static const bool matches =
        K3HyperovalSquareMapperMatchesInitializedGf256();
    return matches;
}

uint32_t Avalanche32(uint32_t value) noexcept
{
    value = (value ^ (value >> 16)) * UINT32_C(0x7feb352d);
    value = (value ^ (value >> 15)) * UINT32_C(0x846ca68b);
    return value ^ (value >> 16);
}

uint32_t ComplementBucketForHighId(uint32_t block_id) noexcept
{
    return static_cast<uint32_t>(
        (static_cast<uint64_t>(Avalanche32(block_id)) *
            UINT64_C(65535)) >> 32);
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
    const K3HyperovalDirection& left,
    const K3HyperovalDirection& right) noexcept
{
    return left.Point == right.Point &&
        left.Alpha == right.Alpha &&
        left.Beta == right.Beta &&
        left.Gamma == right.Gamma;
}

uint8_t Coordinate(
    const K3HyperovalDirection& direction,
    unsigned index) noexcept
{
    if (index == 0u) return direction.Alpha;
    if (index == 1u) return direction.Beta;
    return direction.Gamma;
}

uint8_t TwoByTwo(
    uint8_t a,
    uint8_t b,
    uint8_t c,
    uint8_t d) noexcept
{
    return gf256_mul(a, b) ^ gf256_mul(c, d);
}

bool SetLinearCombination(
    uint8_t* output,
    const uint8_t* const* sources,
    const uint8_t* scales,
    unsigned count,
    uint32_t bytes,
    K3HyperovalBlockOps& counters) noexcept
{
    unsigned first = count;
    for (unsigned i = 0u; i < count; ++i)
    {
        if (scales[i] != 0u) {
            first = i;
            break;
        }
    }
    if (first == count) {
        return false;
    }

    const int kernel_bytes = static_cast<int>(bytes);
    if (scales[first] == 1u)
    {
        std::memcpy(output, sources[first], bytes);
        ++counters.Copies;
    }
    else
    {
        gf256_mul_mem(output, sources[first], scales[first], kernel_bytes);
        ++counters.Scales;
    }

    for (unsigned i = first + 1u; i < count; ++i)
    {
        if (scales[i] == 0u) {
            continue;
        }
        if (scales[i] == 1u)
        {
            gf256_add_mem(output, sources[i], kernel_bytes);
            ++counters.Xors;
        }
        else
        {
            gf256_muladd_mem(
                output, scales[i], sources[i], kernel_bytes);
            ++counters.MulAdds;
        }
    }
    return true;
}

} // namespace

bool K3HyperovalComplementDirectionForBucket(
    uint32_t bucket,
    K3HyperovalDirection& direction_out) noexcept
{
    if (bucket >= kComplementPointCount) {
        return false;
    }

    K3HyperovalDirection direction;
    direction.Point = kLowPointCount + bucket;
    if (bucket < kAffineComplementCount)
    {
        const uint8_t a = static_cast<uint8_t>(bucket / 255u);
        const uint8_t d = static_cast<uint8_t>(1u + bucket % 255u);
        direction.Alpha = 1u;
        direction.Beta = a;
        direction.Gamma = Polynomial14dSquare(a) ^ d;
    }
    else
    {
        direction.Alpha = 0u;
        direction.Beta = 1u;
        direction.Gamma = static_cast<uint8_t>(
            1u + bucket - kAffineComplementCount);
    }
    direction_out = direction;
    return true;
}

K3HyperovalDirection K3HyperovalDirectionForPacketId(
    uint32_t block_id) noexcept
{
    K3HyperovalDirection direction;
    direction.Point = block_id;
    if (block_id == 0u)
    {
        direction.Alpha = 1u;
        return direction;
    }
    if (block_id == 1u)
    {
        direction.Beta = 1u;
        return direction;
    }
    if (block_id == 2u)
    {
        direction.Gamma = 1u;
        return direction;
    }
    if (block_id < kLowPointCount)
    {
        const uint8_t t = static_cast<uint8_t>(block_id - 2u);
        direction.Alpha = 1u;
        direction.Beta = t;
        direction.Gamma = Polynomial14dSquare(t);
        return direction;
    }

    const uint32_t bucket = ComplementBucketForHighId(block_id);
    const bool mapped =
        K3HyperovalComplementDirectionForBucket(bucket, direction);
    (void)mapped;
    return direction;
}

uint8_t K3HyperovalDeterminant(
    const K3HyperovalDirection& first,
    const K3HyperovalDirection& second,
    const K3HyperovalDirection& third) noexcept
{
    const uint8_t minor_alpha = TwoByTwo(
        second.Beta, third.Gamma,
        second.Gamma, third.Beta);
    const uint8_t minor_beta = TwoByTwo(
        second.Alpha, third.Gamma,
        second.Gamma, third.Alpha);
    const uint8_t minor_gamma = TwoByTwo(
        second.Alpha, third.Beta,
        second.Beta, third.Alpha);
    return gf256_mul(first.Alpha, minor_alpha) ^
        gf256_mul(first.Beta, minor_beta) ^
        gf256_mul(first.Gamma, minor_gamma);
}

bool K3HyperovalSquareMapperMatchesInitializedGf256() noexcept
{
    for (uint32_t value = 0u; value < 256u; ++value)
    {
        const uint8_t byte = static_cast<uint8_t>(value);
        if (Polynomial14dSquare(byte) != gf256_sqr(byte)) {
            return false;
        }
    }
    return true;
}

#if defined(WIREHAIR_WH2_K3_HYPEROVAL_TEST_HOOKS)
void SetK3HyperovalAllocationFailureCountdownForTesting(
    int64_t countdown) noexcept
{
    K3HyperovalAllocationFailureCountdown = countdown;
}
#endif

bool K3HyperovalCodec::ValidK3Shape(
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
    return block_count == 3u;
}

void K3HyperovalCodec::Swap(K3HyperovalCodec& other) noexcept
{
    using std::swap;
    Storage.swap(other.Storage);
    AcceptedIds.swap(other.AcceptedIds);
    swap(MessageBytesValue, other.MessageBytesValue);
    swap(BlockBytesValue, other.BlockBytesValue);
    swap(FirstDirection, other.FirstDirection);
    swap(SecondDirection, other.SecondDirection);
    swap(WorkCountersValue, other.WorkCountersValue);
    swap(RankValue, other.RankValue);
    swap(ModeValue, other.ModeValue);
    swap(Decoded, other.Decoded);
}

WirehairResult K3HyperovalCodec::InitializeEncoder(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes) noexcept
{
    if (!message || !ValidK3Shape(message_bytes, block_bytes)) {
        return Wirehair_InvalidInput;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }
    if (!EnsureSquareMapperMatchesInitializedGf256()) {
        return Wirehair_Error;
    }

    try
    {
        size_t storage_bytes = 0u;
        if (!CheckedProduct(
                static_cast<size_t>(block_bytes), 3u, storage_bytes))
        {
            return Wirehair_OOM;
        }
        K3HyperovalCodec next;
        GuardedK3HyperovalAllocation();
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

WirehairResult K3HyperovalCodec::InitializeDecoder(
    uint64_t message_bytes,
    uint32_t block_bytes) noexcept
{
    if (!ValidK3Shape(message_bytes, block_bytes)) {
        return Wirehair_InvalidInput;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }
    if (!EnsureSquareMapperMatchesInitializedGf256()) {
        return Wirehair_Error;
    }

    try
    {
        size_t storage_bytes = 0u;
        if (!CheckedProduct(
                static_cast<size_t>(block_bytes), 6u, storage_bytes))
        {
            return Wirehair_OOM;
        }
        K3HyperovalCodec next;
        GuardedK3HyperovalAllocation();
        next.Storage.assign(storage_bytes, uint8_t{0});
        GuardedK3HyperovalAllocation();
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

uint32_t K3HyperovalCodec::PacketBytes(uint32_t block_id) const noexcept
{
    if (block_id == 2u)
    {
        return static_cast<uint32_t>(
            MessageBytesValue -
                static_cast<uint64_t>(2u) * BlockBytesValue);
    }
    return BlockBytesValue;
}

uint8_t* K3HyperovalCodec::SourceBlock(uint32_t block_id) noexcept
{
    return Storage.data() + static_cast<size_t>(block_id) * BlockBytesValue;
}

const uint8_t* K3HyperovalCodec::SourceBlock(
    uint32_t block_id) const noexcept
{
    return Storage.data() + static_cast<size_t>(block_id) * BlockBytesValue;
}

uint8_t* K3HyperovalCodec::RetainedPacket(uint32_t row) noexcept
{
    return Storage.data() +
        static_cast<size_t>(3u + row) * BlockBytesValue;
}

const uint8_t* K3HyperovalCodec::RetainedPacket(uint32_t row) const noexcept
{
    return Storage.data() +
        static_cast<size_t>(3u + row) * BlockBytesValue;
}

uint8_t* K3HyperovalCodec::IncomingPacketScratch() noexcept
{
    return Storage.data() + static_cast<size_t>(5u) * BlockBytesValue;
}

const uint8_t* K3HyperovalCodec::IncomingPacketScratch() const noexcept
{
    return Storage.data() + static_cast<size_t>(5u) * BlockBytesValue;
}

bool K3HyperovalCodec::HasAcceptedId(uint32_t block_id) const noexcept
{
    return std::find(
        AcceptedIds.begin(), AcceptedIds.end(), block_id) !=
        AcceptedIds.end();
}

void K3HyperovalCodec::CanonicalizePacket(
    const void* block_in,
    uint32_t data_bytes,
    uint8_t* canonical_out) const noexcept
{
    std::memset(canonical_out, 0, BlockBytesValue);
    std::memcpy(canonical_out, block_in, data_bytes);
}

bool K3HyperovalCodec::EvaluatePacket(
    const K3HyperovalDirection& direction,
    uint8_t* canonical_out,
    K3HyperovalBlockOps& counters) const noexcept
{
    const uint8_t* const sources[3] = {
        SourceBlock(0u), SourceBlock(1u), SourceBlock(2u)
    };
    const uint8_t scales[3] = {
        direction.Alpha, direction.Beta, direction.Gamma
    };
    return SetLinearCombination(
        canonical_out,
        sources,
        scales,
        3u,
        BlockBytesValue,
        counters);
}

bool K3HyperovalCodec::RetainedSpanCoefficients(
    const K3HyperovalDirection& direction,
    uint8_t& first_scale,
    uint8_t& second_scale) const noexcept
{
    static const unsigned pairs[3][2] = {
        { 0u, 1u }, { 0u, 2u }, { 1u, 2u }
    };
    for (unsigned pair = 0u; pair < 3u; ++pair)
    {
        const unsigned p = pairs[pair][0];
        const unsigned q = pairs[pair][1];
        const uint8_t first_p = Coordinate(FirstDirection, p);
        const uint8_t first_q = Coordinate(FirstDirection, q);
        const uint8_t second_p = Coordinate(SecondDirection, p);
        const uint8_t second_q = Coordinate(SecondDirection, q);
        const uint8_t determinant = TwoByTwo(
            first_p, second_q, first_q, second_p);
        if (determinant == 0u) {
            continue;
        }

        const uint8_t inverse = gf256_inv(determinant);
        const uint8_t direction_p = Coordinate(direction, p);
        const uint8_t direction_q = Coordinate(direction, q);
        const uint8_t lambda = gf256_mul(
            TwoByTwo(
                direction_p, second_q,
                direction_q, second_p),
            inverse);
        const uint8_t mu = gf256_mul(
            TwoByTwo(
                first_p, direction_q,
                first_q, direction_p),
            inverse);

        for (unsigned coordinate = 0u; coordinate < 3u; ++coordinate)
        {
            const uint8_t expected =
                gf256_mul(lambda, Coordinate(FirstDirection, coordinate)) ^
                gf256_mul(mu, Coordinate(SecondDirection, coordinate));
            if (expected != Coordinate(direction, coordinate)) {
                return false;
            }
        }
        first_scale = lambda;
        second_scale = mu;
        return true;
    }
    return false;
}

bool K3HyperovalCodec::SolveThreeRows(
    const K3HyperovalDirection& third,
    uint8_t determinant) noexcept
{
    if (determinant == 0u) {
        return false;
    }

    const uint8_t a = FirstDirection.Alpha;
    const uint8_t b = FirstDirection.Beta;
    const uint8_t c = FirstDirection.Gamma;
    const uint8_t d = SecondDirection.Alpha;
    const uint8_t e = SecondDirection.Beta;
    const uint8_t f = SecondDirection.Gamma;
    const uint8_t g = third.Alpha;
    const uint8_t h = third.Beta;
    const uint8_t i = third.Gamma;
    const uint8_t inverse = gf256_inv(determinant);

    uint8_t scales[3][3];
    scales[0][0] = gf256_mul(TwoByTwo(e, i, f, h), inverse);
    scales[0][1] = gf256_mul(TwoByTwo(b, i, c, h), inverse);
    scales[0][2] = gf256_mul(TwoByTwo(b, f, c, e), inverse);
    scales[1][0] = gf256_mul(TwoByTwo(d, i, f, g), inverse);
    scales[1][1] = gf256_mul(TwoByTwo(a, i, c, g), inverse);
    scales[1][2] = gf256_mul(TwoByTwo(a, f, c, d), inverse);
    scales[2][0] = gf256_mul(TwoByTwo(d, h, e, g), inverse);
    scales[2][1] = gf256_mul(TwoByTwo(a, h, b, g), inverse);
    scales[2][2] = gf256_mul(TwoByTwo(a, e, b, d), inverse);

    const uint8_t* const packets[3] = {
        RetainedPacket(0u),
        RetainedPacket(1u),
        IncomingPacketScratch()
    };
    for (uint32_t source = 0u; source < 3u; ++source)
    {
        if (!SetLinearCombination(
                SourceBlock(source),
                packets,
                scales[source],
                3u,
                BlockBytesValue,
                WorkCountersValue.Solve))
        {
            return false;
        }
    }
    return true;
}

WirehairResult K3HyperovalCodec::EncodeResult(
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
    if (block_id < 3u) {
        std::memcpy(output, SourceBlock(block_id), data_bytes);
    }
    else if (!EvaluatePacket(
                 K3HyperovalDirectionForPacketId(block_id),
                 output,
                 WorkCountersValue.PacketEvaluation))
    {
        return Wirehair_Error;
    }
    *data_bytes_out = data_bytes;
    return Wirehair_Success;
}

WirehairResult K3HyperovalCodec::DecodeResult(
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
        const K3HyperovalDirection direction =
            K3HyperovalDirectionForPacketId(block_id);
        uint8_t* const incoming = IncomingPacketScratch();
        CanonicalizePacket(block_in, data_bytes, incoming);

        if (Decoded)
        {
            uint8_t* const expected = RetainedPacket(0u);
            if (!EvaluatePacket(
                    direction,
                    expected,
                    WorkCountersValue.PacketEvaluation))
            {
                return Wirehair_Error;
            }
            if (std::memcmp(expected, incoming, BlockBytesValue) != 0) {
                return accepted_id ? Wirehair_InvalidInput : Wirehair_Error;
            }
            if (accepted_id) {
                return Wirehair_Success;
            }
            if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
                return Wirehair_ExtraInsufficient;
            }
            AcceptedIds.push_back(block_id);
            return Wirehair_Success;
        }

        if (RankValue == 0u)
        {
            if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
                return Wirehair_ExtraInsufficient;
            }
            AcceptedIds.push_back(block_id);
            std::memcpy(RetainedPacket(0u), incoming, BlockBytesValue);
            FirstDirection = direction;
            RankValue = 1u;
            return Wirehair_NeedMore;
        }

        if (RankValue == 1u)
        {
            if (SameDirection(direction, FirstDirection))
            {
                if (std::memcmp(
                        RetainedPacket(0u), incoming, BlockBytesValue) != 0)
                {
                    return accepted_id ?
                        Wirehair_InvalidInput : Wirehair_Error;
                }
                if (accepted_id) {
                    return Wirehair_NeedMore;
                }
                if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
                    return Wirehair_ExtraInsufficient;
                }
                AcceptedIds.push_back(block_id);
                return Wirehair_NeedMore;
            }
            if (accepted_id) {
                return Wirehair_Error;
            }
            if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
                return Wirehair_ExtraInsufficient;
            }
            AcceptedIds.push_back(block_id);
            std::memcpy(RetainedPacket(1u), incoming, BlockBytesValue);
            SecondDirection = direction;
            RankValue = 2u;
            return Wirehair_NeedMore;
        }

        if (RankValue != 2u) {
            return Wirehair_Error;
        }

        const uint8_t determinant = K3HyperovalDeterminant(
            FirstDirection, SecondDirection, direction);
        if (determinant == 0u)
        {
            uint8_t first_scale = 0u;
            uint8_t second_scale = 0u;
            if (!RetainedSpanCoefficients(
                    direction, first_scale, second_scale))
            {
                return Wirehair_Error;
            }
            const uint8_t* const packets[2] = {
                RetainedPacket(0u), RetainedPacket(1u)
            };
            const uint8_t scales[2] = { first_scale, second_scale };
            uint8_t* const expected = SourceBlock(0u);
            if (!SetLinearCombination(
                    expected,
                    packets,
                    scales,
                    2u,
                    BlockBytesValue,
                    WorkCountersValue.DependentValidation))
            {
                return Wirehair_Error;
            }
            if (std::memcmp(expected, incoming, BlockBytesValue) != 0) {
                return accepted_id ? Wirehair_InvalidInput : Wirehair_Error;
            }
            if (accepted_id) {
                return Wirehair_NeedMore;
            }
            if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
                return Wirehair_ExtraInsufficient;
            }
            AcceptedIds.push_back(block_id);
            return Wirehair_NeedMore;
        }

        if (accepted_id) {
            return Wirehair_Error;
        }
        if (AcceptedIds.size() >= kMaxAcceptedPacketIds) {
            return Wirehair_ExtraInsufficient;
        }
        if (!SolveThreeRows(direction, determinant)) {
            return Wirehair_Error;
        }

        const uint32_t tail_bytes = static_cast<uint32_t>(
            MessageBytesValue -
                static_cast<uint64_t>(2u) * BlockBytesValue);
        for (uint32_t i = tail_bytes; i < BlockBytesValue; ++i)
        {
            if (SourceBlock(2u)[i] != 0u) {
                return Wirehair_Error;
            }
        }

        AcceptedIds.push_back(block_id);
        RankValue = 3u;
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

WirehairResult K3HyperovalCodec::RecoverResult(
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

bool K3HyperovalCodec::IsEncoder() const noexcept
{
    return ModeValue == Mode::Encoder;
}

bool K3HyperovalCodec::IsDecoder() const noexcept
{
    return ModeValue == Mode::Decoder;
}

bool K3HyperovalCodec::IsDecoded() const noexcept
{
    return IsDecoder() && Decoded;
}

uint32_t K3HyperovalCodec::BlockCount() const noexcept
{
    return ModeValue == Mode::Uninitialized ? 0u : 3u;
}

uint32_t K3HyperovalCodec::BlockBytes() const noexcept
{
    return ModeValue == Mode::Uninitialized ? 0u : BlockBytesValue;
}

uint64_t K3HyperovalCodec::MessageBytes() const noexcept
{
    return ModeValue == Mode::Uninitialized ? 0u : MessageBytesValue;
}

uint32_t K3HyperovalCodec::AcceptedIdCount() const noexcept
{
    return IsDecoder() ? static_cast<uint32_t>(AcceptedIds.size()) : 0u;
}

uint32_t K3HyperovalCodec::Rank() const noexcept
{
    return IsDecoder() ? RankValue : 0u;
}

K3HyperovalWorkCounters K3HyperovalCodec::WorkCounters() const noexcept
{
    return WorkCountersValue;
}

void K3HyperovalCodec::ResetWorkCounters() noexcept
{
    WorkCountersValue = K3HyperovalWorkCounters();
}

} // namespace wirehair_wh2_bench
