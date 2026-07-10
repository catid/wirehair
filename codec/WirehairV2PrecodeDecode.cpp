#include "WirehairV2PrecodeDecode.h"

#include "../WirehairTools.h"
#include "../gf256.h"
#include "WirehairV2Plan.h"

#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wirehair_v2 {
namespace {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local int64_t DecoderAllocationFailureCountdown = -1;

void GuardedDecoderAllocation()
{
    if (DecoderAllocationFailureCountdown == 0) {
        throw std::bad_alloc();
    }
    if (DecoderAllocationFailureCountdown > 0) {
        --DecoderAllocationFailureCountdown;
    }
}
#else
void GuardedDecoderAllocation() {}
#endif

bool MessageBlockCount(
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint32_t& block_count)
{
    block_count = 0u;
    if (message_bytes == 0u || block_bytes == 0u ||
        block_bytes > 0x7fffffffu)
    {
        return false;
    }
    const uint64_t count = message_bytes / block_bytes +
        (message_bytes % block_bytes != 0u ? 1u : 0u);
    if (count < CAT_WIREHAIR_MIN_N || count > CAT_WIREHAIR_MAX_N) {
        return false;
    }
    block_count = (uint32_t)count;
    return true;
}

uint32_t PacketDataBytes(
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint32_t block_count,
    uint32_t block_id)
{
    if (block_id >= block_count || block_id + 1u < block_count) {
        return block_bytes;
    }
    const uint64_t offset = (uint64_t)block_id * block_bytes;
    const uint64_t remaining = message_bytes - offset;
    return remaining < block_bytes ? (uint32_t)remaining : block_bytes;
}

} // namespace

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void SetDecoderAllocationFailureCountdownForTesting(int64_t countdown)
{
    DecoderAllocationFailureCountdown = countdown;
}
#endif

MessagePrecodeDecoder::MessagePrecodeDecoder()
{
}

void MessagePrecodeDecoder::Swap(MessagePrecodeDecoder& other) noexcept
{
    using std::swap;
    swap(ProfileValue, other.ProfileValue);
    swap(OptionsValue, other.OptionsValue);
    swap(PacketConfigValue, other.PacketConfigValue);
    swap(SystemValue.Params, other.SystemValue.Params);
    SystemValue.StaircaseRows.swap(other.SystemValue.StaircaseRows);
    SystemValue.DenseRowColumns.swap(other.SystemValue.DenseRowColumns);
    ReceivedBlockIds.swap(other.ReceivedBlockIds);
    ReceivedBlockStorage.swap(other.ReceivedBlockStorage);
    ReceivedIds.swap(other.ReceivedIds);
    IntermediateBlockStorage.swap(other.IntermediateBlockStorage);
    swap(SolveStatsValue, other.SolveStatsValue);
    swap(MessageBytesValue, other.MessageBytesValue);
    swap(BlockBytesValue, other.BlockBytesValue);
    swap(ReceivedCountValue, other.ReceivedCountValue);
    swap(SolveAttemptCountValue, other.SolveAttemptCountValue);
    swap(PacketSeedAttemptValue, other.PacketSeedAttemptValue);
    swap(LastSolveResult, other.LastSolveResult);
    swap(Initialized, other.Initialized);
    swap(Decoded, other.Decoded);
}

WirehairResult MessagePrecodeDecoder::InitializeResult(
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override,
    const MessagePrecodeEncoderOptions* options)
{
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }
    uint32_t block_count = 0u;
    if (message_bytes >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        !MessageBlockCount(message_bytes, block_bytes, block_count))
    {
        return Wirehair_InvalidInput;
    }
    const SeedProfile profile = seed_override ? *seed_override :
        SelectSeedProfile(block_count, block_bytes);
    if (profile.BlockCount != block_count ||
        profile.BlockBytes != block_bytes)
    {
        return Wirehair_InvalidInput;
    }
    MessagePrecodeEncoderOptions opts;
    if (!ResolveMessagePrecodeOptions(profile, options, opts)) {
        return Wirehair_InvalidInput;
    }

    try
    {
        PrecodeParams params;
        PacketRowConfig base_config;
        if (!ResolveMessagePrecodeConfiguration(
                profile, opts, params, base_config))
        {
            return Wirehair_InvalidInput;
        }
        PacketRowConfig selected_config;
        PrecodeSystem system;
        uint32_t packet_seed_attempt = profile.V2SeedSelected ?
            profile.V2SeedAttempt : 0u;
        if (profile.V2SeedSelected)
        {
            GuardedDecoderAllocation();
            if (!BuildPrecodeSystem(params, system))
            {
                return Wirehair_InvalidInput;
            }
            selected_config = base_config;
        }
        else
        {
            GuardedDecoderAllocation();
            const WirehairResult select_result =
                SelectSystematicConfiguration(
                    params,
                    base_config,
                    system,
                    selected_config,
                    &packet_seed_attempt);
            if (select_result != Wirehair_Success) {
                return select_result;
            }
        }

        MessagePrecodeDecoder next;
        next.ProfileValue = profile;
        BindMessagePrecodeProfile(
            next.ProfileValue,
            opts,
            system,
            selected_config,
            packet_seed_attempt);
        next.OptionsValue = opts;
        next.PacketConfigValue = selected_config;
        next.SystemValue = std::move(system);
        GuardedDecoderAllocation();
        next.ReceivedBlockIds.reserve((size_t)block_count + 32u);
        const uint64_t receive_capacity =
            ((uint64_t)block_count + 32u) * block_bytes;
        if (receive_capacity >
            (uint64_t)std::numeric_limits<size_t>::max())
        {
            return Wirehair_InvalidInput;
        }
        GuardedDecoderAllocation();
        next.ReceivedBlockStorage.reserve((size_t)receive_capacity);
        GuardedDecoderAllocation();
        next.ReceivedIds.reserve((size_t)block_count + 32u);
        next.MessageBytesValue = message_bytes;
        next.BlockBytesValue = block_bytes;
        next.PacketSeedAttemptValue = packet_seed_attempt;
        next.Initialized = true;
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

WirehairResult MessagePrecodeDecoder::AttemptSolve()
{
    const uint32_t K = ProfileValue.BlockCount;
    if (Decoded) {
        return Wirehair_Success;
    }
    if (ReceivedBlockIds.size() < K) {
        LastSolveResult = Wirehair_NeedMore;
        return Wirehair_NeedMore;
    }

    try
    {
        std::vector<SolvePacket> packets;
        packets.reserve(ReceivedBlockIds.size());
        for (size_t i = 0; i < ReceivedBlockIds.size(); ++i)
        {
            SolvePacket packet;
            packet.BlockId = ReceivedBlockIds[i];
            packet.Data = ReceivedBlockStorage.data() +
                i * BlockBytesValue;
            packets.push_back(packet);
        }
        std::vector<uint8_t> intermediate;
        PrecodeSolveStats solve_stats;
        ++SolveAttemptCountValue;
        // Keep the deterministic failure point at the solve boundary so the
        // test hook models a transient solver OOM (and counts as an attempt).
        GuardedDecoderAllocation();
        const WirehairResult result = SolvePrecodeSystem(
            SystemValue,
            PacketConfigValue,
            packets,
            BlockBytesValue,
            intermediate,
            &solve_stats);
        solve_stats.PacketSeedAttempt = PacketSeedAttemptValue;
        SolveStatsValue = solve_stats;
        if (result == Wirehair_Success)
        {
            IntermediateBlockStorage.swap(intermediate);
            Decoded = true;
            ReceivedCountValue = (uint32_t)ReceivedBlockIds.size();
            std::vector<uint32_t>().swap(ReceivedBlockIds);
            std::vector<uint8_t>().swap(ReceivedBlockStorage);
            std::unordered_set<uint32_t>().swap(ReceivedIds);
        }
        LastSolveResult = result;
        return result;
    }
    catch (const std::bad_alloc&) {
        LastSolveResult = Wirehair_OOM;
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        LastSolveResult = Wirehair_OOM;
        return Wirehair_OOM;
    }
}

WirehairResult MessagePrecodeDecoder::DecodeResult(
    uint32_t block_id,
    const void* block_in,
    uint32_t data_bytes)
{
    if (!Initialized || !block_in) {
        return Wirehair_InvalidInput;
    }
    const uint32_t expected = PacketDataBytes(
        MessageBytesValue,
        BlockBytesValue,
        ProfileValue.BlockCount,
        block_id);
    if (data_bytes != expected) {
        return Wirehair_InvalidInput;
    }

    if (Decoded)
    {
        try
        {
            std::vector<uint8_t> expected_data(BlockBytesValue, 0u);
            if (!EvaluatePacketBlockForValidatedSystem(
                    SystemValue,
                    PacketConfigValue,
                    IntermediateBlockStorage.data(),
                    BlockBytesValue,
                    block_id,
                    expected_data.data()))
            {
                return Wirehair_Error;
            }
            return std::memcmp(expected_data.data(), block_in, data_bytes) == 0 ?
                Wirehair_Success : Wirehair_Error;
        }
        catch (const std::bad_alloc&) {
            return Wirehair_OOM;
        }
    }

    if (ReceivedIds.find(block_id) != ReceivedIds.end())
    {
        for (size_t i = 0; i < ReceivedBlockIds.size(); ++i)
        {
            if (ReceivedBlockIds[i] == block_id)
            {
                if (std::memcmp(
                    ReceivedBlockStorage.data() +
                        i * BlockBytesValue,
                    block_in, data_bytes) != 0)
                {
                    return Wirehair_InvalidInput;
                }
                // A solve allocation failure is transient and the packet is
                // already committed, so an identical resubmission is the
                // retry signal.  Rank-deficient NeedMore is intentionally
                // cached: retrying it without another equation cannot help.
                if (LastSolveResult == Wirehair_OOM &&
                    ReceivedBlockIds.size() >= ProfileValue.BlockCount)
                {
                    return AttemptSolve();
                }
                return LastSolveResult;
            }
        }
        return Wirehair_Error;
    }

    if (ReceivedBlockIds.size() >=
        (size_t)ProfileValue.BlockCount + 1024u)
    {
        return Wirehair_ExtraInsufficient;
    }
    if (ReceivedBlockStorage.size() >
        std::numeric_limits<size_t>::max() - BlockBytesValue)
    {
        return Wirehair_OOM;
    }

    try
    {
        const std::pair<std::unordered_set<uint32_t>::iterator, bool> inserted =
            ReceivedIds.insert(block_id);
        if (!inserted.second) {
            return LastSolveResult;
        }
        const size_t old_storage_size = ReceivedBlockStorage.size();
        try
        {
            ReceivedBlockStorage.resize(
                old_storage_size + BlockBytesValue, 0u);
            std::memcpy(
                ReceivedBlockStorage.data() + old_storage_size,
                block_in,
                data_bytes);
            ReceivedBlockIds.push_back(block_id);
        }
        catch (...)
        {
            ReceivedBlockStorage.resize(old_storage_size);
            ReceivedIds.erase(inserted.first);
            throw;
        }
        return AttemptSolve();
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

WirehairResult MessagePrecodeDecoder::RecoverResult(
    void* message_out,
    uint64_t message_bytes) const
{
    if (!Initialized || !message_out || message_bytes != MessageBytesValue) {
        return Wirehair_InvalidInput;
    }
    if (!Decoded) {
        return Wirehair_NeedMore;
    }
    try
    {
        GuardedDecoderAllocation();
        std::vector<uint8_t> block(BlockBytesValue, 0u);
        uint8_t* output = static_cast<uint8_t*>(message_out);
        const uint32_t K = ProfileValue.BlockCount;
        for (uint32_t block_id = 0; block_id < K; ++block_id)
        {
            if (!EvaluatePacketBlockForValidatedSystem(
                    SystemValue,
                    PacketConfigValue,
                    IntermediateBlockStorage.data(),
                    BlockBytesValue,
                    block_id,
                    block.data()))
            {
                return Wirehair_Error;
            }
            const uint32_t bytes = PacketDataBytes(
                MessageBytesValue, BlockBytesValue, K, block_id);
            std::memcpy(
                output + (size_t)block_id * BlockBytesValue,
                block.data(),
                bytes);
        }
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
}

bool MessagePrecodeDecoder::IsInitialized() const { return Initialized; }
bool MessagePrecodeDecoder::IsDecoded() const { return Decoded; }
uint32_t MessagePrecodeDecoder::ReceivedCount() const
{
    return Decoded ? ReceivedCountValue :
        (uint32_t)ReceivedBlockIds.size();
}
uint32_t MessagePrecodeDecoder::SolveAttemptCount() const
{
    return Initialized ? SolveAttemptCountValue : 0u;
}
uint32_t MessagePrecodeDecoder::PacketSeedAttempt() const
{
    return Initialized ? PacketSeedAttemptValue : 0u;
}
uint32_t MessagePrecodeDecoder::PacketPeelSeed() const
{
    return Initialized ? PacketConfigValue.PeelSeed : 0u;
}
uint64_t MessagePrecodeDecoder::MessageBytes() const
{
    return Initialized ? MessageBytesValue : 0u;
}
uint32_t MessagePrecodeDecoder::BlockBytes() const
{
    return Initialized ? BlockBytesValue : 0u;
}
const SeedProfile& MessagePrecodeDecoder::Profile() const
{
    return ProfileValue;
}
const MessagePrecodeEncoderOptions& MessagePrecodeDecoder::Options() const
{
    return OptionsValue;
}
const PrecodeSolveStats& MessagePrecodeDecoder::SolveStats() const
{
    return SolveStatsValue;
}
const PrecodeSystem& MessagePrecodeDecoder::System() const
{
    return SystemValue;
}
const uint8_t* MessagePrecodeDecoder::IntermediateBlocks() const
{
    return Decoded ? IntermediateBlockStorage.data() : nullptr;
}

} // namespace wirehair_v2
