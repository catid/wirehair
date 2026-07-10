#include "WirehairV2PrecodeDecode.h"

#include "../WirehairTools.h"
#include "../gf256.h"
#include "WirehairV2Plan.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wirehair_v2 {
namespace {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local int64_t DecoderAllocationFailureCountdown = -1;
thread_local bool DecoderIncrementalResumeEnabled = true;

void GuardedDecoderAllocation()
{
    if (DecoderAllocationFailureCountdown == 0) {
        throw std::bad_alloc();
    }
    if (DecoderAllocationFailureCountdown > 0) {
        --DecoderAllocationFailureCountdown;
    }
}

bool IncrementalResumeEnabled()
{
    return DecoderIncrementalResumeEnabled;
}
#else
void GuardedDecoderAllocation() {}
bool IncrementalResumeEnabled() { return true; }
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

bool ResumeFitsMemoryPolicy(
    const PrecodeSolveResumeState& state,
    size_t receive_block_capacity,
    size_t receive_id_capacity,
    size_t pending_block_capacity)
{
    if (!state.Active ||
        receive_id_capacity >
            std::numeric_limits<size_t>::max() / sizeof(uint32_t))
    {
        return false;
    }
    const size_t id_bytes = receive_id_capacity * sizeof(uint32_t);
    if (receive_block_capacity >
        std::numeric_limits<size_t>::max() - id_bytes)
    {
        return false;
    }
    const size_t released_bytes = receive_block_capacity + id_bytes;
    const size_t allowed_extra = released_bytes / 4u;
    const size_t checkpoint_bytes = state.PersistentBytes();
    if (checkpoint_bytes >
        std::numeric_limits<size_t>::max() - pending_block_capacity)
    {
        return false;
    }
    const size_t retained_bytes = checkpoint_bytes + pending_block_capacity;
    return retained_bytes <= released_bytes ||
        retained_bytes - released_bytes <= allowed_extra;
}

} // namespace

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void SetDecoderAllocationFailureCountdownForTesting(int64_t countdown)
{
    DecoderAllocationFailureCountdown = countdown;
}

void SetDecoderIncrementalResumeEnabledForTesting(bool enabled)
{
    DecoderIncrementalResumeEnabled = enabled;
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
    ResumeState.Swap(other.ResumeState);
    PendingPacketStorage.swap(other.PendingPacketStorage);
    swap(PendingPacketId, other.PendingPacketId);
    IntermediateBlockStorage.swap(other.IntermediateBlockStorage);
    swap(SolveStatsValue, other.SolveStatsValue);
    swap(MessageBytesValue, other.MessageBytesValue);
    swap(BlockBytesValue, other.BlockBytesValue);
    swap(ReceivedCountValue, other.ReceivedCountValue);
    swap(SolveAttemptCountValue, other.SolveAttemptCountValue);
    swap(PacketSeedAttemptValue, other.PacketSeedAttemptValue);
    swap(LastSolveResult, other.LastSolveResult);
    swap(PendingPacket, other.PendingPacket);
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
    if (ReceivedIds.size() < K) {
        LastSolveResult = Wirehair_NeedMore;
        return Wirehair_NeedMore;
    }

    try
    {
        if (ResumeState.Active)
        {
            if (!PendingPacket ||
                PendingPacketStorage.size() != BlockBytesValue)
            {
                return LastSolveResult;
            }
            std::vector<uint8_t> intermediate;
            // Seed from the last committed attempt as a defense for every
            // early-return path.  ResumePrecodeSystem also republishes the
            // checkpoint counters on allocation failure, so a retryable OOM
            // cannot erase diagnostics while the algebra remains unchanged.
            PrecodeSolveStats solve_stats = SolveStatsValue;
            ++SolveAttemptCountValue;
            GuardedDecoderAllocation();
            const WirehairResult result = ResumePrecodeSystem(
                SystemValue,
                PacketConfigValue,
                PendingPacketId,
                PendingPacketStorage.data(),
                BlockBytesValue,
                ResumeState,
                intermediate,
                &solve_stats,
                true);
            solve_stats.PacketSeedAttempt = PacketSeedAttemptValue;
            SolveStatsValue = solve_stats;
            if (result == Wirehair_Success)
            {
                IntermediateBlockStorage.swap(intermediate);
                Decoded = true;
                ReceivedCountValue = (uint32_t)ReceivedIds.size();
                std::unordered_set<uint32_t>().swap(ReceivedIds);
                std::vector<uint8_t>().swap(PendingPacketStorage);
                PendingPacket = false;
            }
            else if (result == Wirehair_NeedMore)
            {
                PendingPacket = false;
            }
            LastSolveResult = result;
            return result;
        }

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
        // A cold solve has the same transactional stats contract: preserve
        // the last committed counters if the solver cannot allocate.
        PrecodeSolveStats solve_stats = SolveStatsValue;
        PrecodeSolveResumeState resume_state;
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
            &solve_stats,
            IncrementalResumeEnabled() ? &resume_state : nullptr);
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
        else if (result == Wirehair_NeedMore && resume_state.Active)
        {
            // Allocate the sole pending-packet slot before releasing the cold
            // receive buffers.  The allocator-selected capacity, rather than
            // an assumed size, is included in the 25% policy.  Failure leaves
            // the complete cold state available for an identical retry.
            GuardedDecoderAllocation();
            std::vector<uint8_t> pending_storage(BlockBytesValue, 0u);
            if (ResumeFitsMemoryPolicy(
                    resume_state,
                    ReceivedBlockStorage.capacity(),
                    ReceivedBlockIds.capacity(),
                    pending_storage.capacity()))
            {
                ResumeState.Swap(resume_state);
                PendingPacketStorage.swap(pending_storage);
                std::vector<uint32_t>().swap(ReceivedBlockIds);
                std::vector<uint8_t>().swap(ReceivedBlockStorage);
            }
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

    if (ResumeState.Active)
    {
        const bool duplicate =
            ReceivedIds.find(block_id) != ReceivedIds.end();
        if (duplicate)
        {
            if (PendingPacket && block_id == PendingPacketId)
            {
                if (std::memcmp(
                        PendingPacketStorage.data(),
                        block_in,
                        data_bytes) != 0)
                {
                    return Wirehair_InvalidInput;
                }
            }
            else
            {
                try
                {
                    std::vector<uint8_t> duplicate_data(
                        BlockBytesValue, 0u);
                    std::memcpy(
                        duplicate_data.data(), block_in, data_bytes);
                    std::vector<uint8_t> ignored;
                    const WirehairResult consistency = ResumePrecodeSystem(
                        SystemValue,
                        PacketConfigValue,
                        block_id,
                        duplicate_data.data(),
                        BlockBytesValue,
                        ResumeState,
                        ignored,
                        nullptr,
                        false);
                    if (consistency == Wirehair_OOM) {
                        return Wirehair_OOM;
                    }
                    if (consistency != Wirehair_NeedMore) {
                        return Wirehair_InvalidInput;
                    }
                }
                catch (const std::bad_alloc&) {
                    return Wirehair_OOM;
                }
                catch (const std::length_error&) {
                    return Wirehair_OOM;
                }
            }
            if (LastSolveResult == Wirehair_OOM && PendingPacket) {
                return AttemptSolve();
            }
            return LastSolveResult;
        }

        if (LastSolveResult == Wirehair_Error) {
            return Wirehair_Error;
        }
        // A transient solve OOM leaves exactly one accepted equation pending.
        // Before accepting a different id, retry that equation so it cannot be
        // overwritten while its id remains in ReceivedIds.  If it completes
        // the decode, recursively validate the current packet against the
        // completed solution; this recursion is bounded to one level because
        // successful completion clears ResumeState.  Terminal Error keeps the
        // pending bytes only to distinguish an exact retry from a conflict and
        // is handled above without re-running the poisoned solve.
        if (PendingPacket)
        {
            const WirehairResult pending_result = AttemptSolve();
            if (pending_result == Wirehair_Success) {
                return DecodeResult(block_id, block_in, data_bytes);
            }
            if (pending_result != Wirehair_NeedMore) {
                return pending_result;
            }
        }
        if (ReceivedIds.size() >=
            (size_t)ProfileValue.BlockCount + 1024u)
        {
            return Wirehair_ExtraInsufficient;
        }
        try
        {
            if (PendingPacketStorage.size() != BlockBytesValue) {
                return Wirehair_Error;
            }
            std::fill(
                PendingPacketStorage.begin(),
                PendingPacketStorage.end(),
                uint8_t{0});
            std::memcpy(
                PendingPacketStorage.data(), block_in, data_bytes);
            const std::pair<
                std::unordered_set<uint32_t>::iterator, bool> inserted =
                    ReceivedIds.insert(block_id);
            if (!inserted.second) {
                return LastSolveResult;
            }
            PendingPacketId = block_id;
            PendingPacket = true;
            return AttemptSolve();
        }
        catch (const std::bad_alloc&) {
            PendingPacket = false;
            return Wirehair_OOM;
        }
        catch (const std::length_error&) {
            PendingPacket = false;
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

    // An inconsistent accepted equation set cannot be repaired by adding more
    // equations.  Freeze the terminal set on the cold fallback too, matching
    // the checkpoint path while still validating duplicates above.  This also
    // prevents pointless K-sized re-solves after a terminal Error.
    if (LastSolveResult == Wirehair_Error) {
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
        (uint32_t)ReceivedIds.size();
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool MessagePrecodeDecoder::HasIncrementalResumeStateForTesting() const
{
    return ResumeState.Active;
}

size_t MessagePrecodeDecoder::IncrementalResumeBytesForTesting() const
{
    return ResumeState.Active ?
        ResumeState.PersistentBytes() + PendingPacketStorage.capacity() : 0u;
}

size_t MessagePrecodeDecoder::ColdReceiveCapacityBytesForTesting() const
{
    return ReceivedBlockStorage.capacity() +
        ReceivedBlockIds.capacity() * sizeof(uint32_t);
}
#endif

} // namespace wirehair_v2
