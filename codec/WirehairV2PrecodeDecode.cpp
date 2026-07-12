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

uint32_t PacketSlotTable::Hash(uint32_t packet_id)
{
    uint32_t x = packet_id;
    x ^= x >> 16;
    x *= UINT32_C(0x7feb352d);
    x ^= x >> 15;
    x *= UINT32_C(0x846ca68b);
    x ^= x >> 16;
    return x;
}

bool PacketSlotTable::Initialize(
    size_t initial_entries,
    size_t max_entries)
{
    if (initial_entries == 0u || initial_entries > max_entries ||
        max_entries > std::numeric_limits<size_t>::max() / 2u)
    {
        return false;
    }
    const size_t minimum_capacity = initial_entries * 2u;
    size_t capacity = 1u;
    while (capacity < minimum_capacity)
    {
        if (capacity > std::numeric_limits<size_t>::max() / 2u) {
            return false;
        }
        capacity *= 2u;
    }

    try
    {
        PacketSlotTable next;
        next.Keys.resize(capacity);
        next.Slots.resize(capacity);
        next.Occupied.assign(capacity, uint8_t{0});
        next.EntryLimit = max_entries;
        Swap(next);
        return true;
    }
    catch (const std::bad_alloc&) {
        return false;
    }
    catch (const std::length_error&) {
        return false;
    }
}

void PacketSlotTable::Grow()
{
    if (Keys.empty() || Keys.size() >
            std::numeric_limits<size_t>::max() / 2u)
    {
        throw std::length_error("packet slot table capacity overflow");
    }
    PacketSlotTable next;
    const size_t capacity = Keys.size() * 2u;
    next.Keys.resize(capacity);
    next.Slots.resize(capacity);
    next.Occupied.assign(capacity, uint8_t{0});
    next.EntryLimit = EntryLimit;
    const size_t mask = capacity - 1u;
    for (size_t i = 0u; i < Keys.size(); ++i)
    {
        if (!Occupied[i]) {
            continue;
        }
        size_t index = (size_t)Hash(Keys[i]) & mask;
        while (next.Occupied[index]) {
            index = (index + 1u) & mask;
        }
        next.Keys[index] = Keys[i];
        next.Slots[index] = Slots[i];
        next.Occupied[index] = 1u;
        ++next.EntryCount;
    }
    Swap(next);
}

bool PacketSlotTable::Find(uint32_t packet_id, uint32_t* slot_out) const
{
    if (Keys.empty()) {
        return false;
    }
    const size_t mask = Keys.size() - 1u;
    size_t index = (size_t)Hash(packet_id) & mask;
    for (size_t probes = 0u; probes < Keys.size(); ++probes)
    {
        if (!Occupied[index]) {
            return false;
        }
        if (Keys[index] == packet_id)
        {
            if (slot_out) {
                *slot_out = Slots[index];
            }
            return true;
        }
        index = (index + 1u) & mask;
    }
    return false;
}

bool PacketSlotTable::Insert(uint32_t packet_id, uint32_t slot)
{
    if (Keys.empty() || EntryCount >= EntryLimit) {
        return false;
    }
    if (EntryCount >= Keys.size() / 2u)
    {
        if (Find(packet_id)) {
            return false;
        }
        GuardedDecoderAllocation();
        Grow();
    }
    const size_t mask = Keys.size() - 1u;
    size_t index = (size_t)Hash(packet_id) & mask;
    for (size_t probes = 0u; probes < Keys.size(); ++probes)
    {
        if (!Occupied[index])
        {
            Keys[index] = packet_id;
            Slots[index] = slot;
            Occupied[index] = 1u;
            ++EntryCount;
            return true;
        }
        if (Keys[index] == packet_id) {
            return false;
        }
        index = (index + 1u) & mask;
    }
    return false;
}

bool PacketSlotTable::Erase(uint32_t packet_id)
{
    if (Keys.empty()) {
        return false;
    }
    const size_t mask = Keys.size() - 1u;
    size_t hole = (size_t)Hash(packet_id) & mask;
    for (size_t probes = 0u; probes < Keys.size(); ++probes)
    {
        if (!Occupied[hole]) {
            return false;
        }
        if (Keys[hole] == packet_id) {
            break;
        }
        hole = (hole + 1u) & mask;
    }
    if (!Occupied[hole] || Keys[hole] != packet_id) {
        return false;
    }

    // Backward-shift deletion preserves every probe chain without tombstones.
    size_t scan = (hole + 1u) & mask;
    while (Occupied[scan])
    {
        const size_t home = (size_t)Hash(Keys[scan]) & mask;
        const size_t scan_distance = (scan - home) & mask;
        const size_t hole_distance = (scan - hole) & mask;
        if (scan_distance >= hole_distance)
        {
            Keys[hole] = Keys[scan];
            Slots[hole] = Slots[scan];
            Occupied[hole] = 1u;
            hole = scan;
        }
        scan = (scan + 1u) & mask;
    }
    Occupied[hole] = 0u;
    --EntryCount;
    return true;
}

void PacketSlotTable::ClearAndRelease() noexcept
{
    PacketSlotTable empty;
    Swap(empty);
}

void PacketSlotTable::Swap(PacketSlotTable& other) noexcept
{
    Keys.swap(other.Keys);
    Slots.swap(other.Slots);
    Occupied.swap(other.Occupied);
    std::swap(EntryCount, other.EntryCount);
    std::swap(EntryLimit, other.EntryLimit);
}

size_t PacketSlotTable::StorageBytes() const
{
    return Keys.capacity() * sizeof(uint32_t) +
        Slots.capacity() * sizeof(uint32_t) + Occupied.capacity();
}

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
    swap(PacketRuntimeValue, other.PacketRuntimeValue);
    swap(SystemValue.Params, other.SystemValue.Params);
    SystemValue.StaircaseRows.swap(other.SystemValue.StaircaseRows);
    SystemValue.DenseRowColumns.swap(other.SystemValue.DenseRowColumns);
    ReceivedBlockIds.swap(other.ReceivedBlockIds);
    ReceivedBlockStorage.swap(other.ReceivedBlockStorage);
    ReceivedSlots.Swap(other.ReceivedSlots);
    ResumeState.Swap(other.ResumeState);
    PendingPacketStorage.swap(other.PendingPacketStorage);
    swap(PendingPacketId, other.PendingPacketId);
    IntermediateBlockStorage.swap(other.IntermediateBlockStorage);
    SystematicPacketCache.swap(other.SystematicPacketCache);
    swap(SystematicPacketCacheSize, other.SystematicPacketCacheSize);
    HaveSystematicPacket.swap(other.HaveSystematicPacket);
    swap(CachedSystematicPacketCountValue,
        other.CachedSystematicPacketCountValue);
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
        const uint64_t precode_count =
            (uint64_t)next.SystemValue.Params.Staircase +
            next.SystemValue.Params.DenseRows +
            next.SystemValue.Params.HeavyRows;
        if (precode_count > UINT32_MAX ||
            !next.PacketRuntimeValue.Initialize(
                block_count,
                (uint32_t)precode_count,
                selected_config.MixCount))
        {
            return Wirehair_InvalidInput;
        }
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
        if (!next.ReceivedSlots.Initialize(
                (size_t)block_count + 32u,
                (size_t)block_count + 1024u))
        {
            return Wirehair_OOM;
        }
        if (opts.CacheReceivedSystematicPackets)
        {
            GuardedDecoderAllocation();
            next.SystematicPacketCache.reset(
                new uint8_t[(size_t)message_bytes]);
            next.SystematicPacketCacheSize = (size_t)message_bytes;
            GuardedDecoderAllocation();
            next.HaveSystematicPacket.assign(block_count, uint8_t{0});
        }
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
    if (ReceivedSlots.Size() < K) {
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
                ReceivedCountValue = (uint32_t)ReceivedSlots.Size();
                ReceivedSlots.ClearAndRelease();
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
        const WirehairResult result = SolvePrecodeSystemWithRuntime(
            SystemValue,
            PacketConfigValue,
            PacketRuntimeValue,
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
            ReceivedSlots.ClearAndRelease();
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
            if (!EvaluatePacketBlockForValidatedSystemWithRuntime(
                    SystemValue,
                    PacketConfigValue,
                    PacketRuntimeValue,
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
            ReceivedSlots.Find(block_id);
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
        // overwritten while its id remains in ReceivedSlots.  If it completes
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
        if (ReceivedSlots.Size() >=
            (size_t)ProfileValue.BlockCount + 1024u)
        {
            return Wirehair_ExtraInsufficient;
        }
        try
        {
            if (PendingPacketStorage.size() != BlockBytesValue) {
                return Wirehair_Error;
            }
            std::memcpy(
                PendingPacketStorage.data(), block_in, data_bytes);
            // Recovery packets and complete systematic packets overwrite the
            // full slot.  Only the one possible partial systematic tail needs
            // deterministic zero padding before it enters the algebra.
            std::fill(
                PendingPacketStorage.begin() + data_bytes,
                PendingPacketStorage.end(),
                uint8_t{0});
            if (!ReceivedSlots.Insert(block_id, UINT32_MAX)) {
                return LastSolveResult;
            }
            if (SystematicPacketCache && block_id <
                    ProfileValue.BlockCount)
            {
                std::memcpy(
                    SystematicPacketCache.get() +
                        (size_t)block_id * BlockBytesValue,
                    block_in,
                    data_bytes);
                if (!HaveSystematicPacket[block_id])
                {
                    HaveSystematicPacket[block_id] = 1u;
                    ++CachedSystematicPacketCountValue;
                }
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

    uint32_t duplicate_slot = 0u;
    if (ReceivedSlots.Find(block_id, &duplicate_slot))
    {
        const uint32_t slot = duplicate_slot;
        if (slot >= ReceivedBlockIds.size() ||
            slot >= ReceivedBlockStorage.size() / BlockBytesValue ||
            ReceivedBlockIds[slot] != block_id)
        {
            return Wirehair_Error;
        }
        if (std::memcmp(
                ReceivedBlockStorage.data() +
                    (size_t)slot * BlockBytesValue,
                block_in, data_bytes) != 0)
        {
            return Wirehair_InvalidInput;
        }
        // A solve allocation failure is transient and the packet is already
        // committed, so an identical resubmission is the retry signal.
        // Rank-deficient NeedMore is intentionally cached: retrying it without
        // another equation cannot help.
        if (LastSolveResult == Wirehair_OOM &&
            ReceivedBlockIds.size() >= ProfileValue.BlockCount)
        {
            return AttemptSolve();
        }
        return LastSolveResult;
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
        const uint32_t slot = (uint32_t)ReceivedBlockIds.size();
        if (!ReceivedSlots.Insert(block_id, slot)) {
            return LastSolveResult;
        }
        const size_t old_storage_size = ReceivedBlockStorage.size();
        try
        {
            const uint8_t* const bytes =
                static_cast<const uint8_t*>(block_in);
            ReceivedBlockStorage.insert(
                ReceivedBlockStorage.end(), bytes, bytes + data_bytes);
            if (data_bytes < BlockBytesValue) {
                ReceivedBlockStorage.resize(
                    old_storage_size + BlockBytesValue, 0u);
            }
            ReceivedBlockIds.push_back(block_id);
            if (SystematicPacketCache && block_id <
                    ProfileValue.BlockCount)
            {
                std::memcpy(
                    SystematicPacketCache.get() +
                        (size_t)block_id * BlockBytesValue,
                    block_in,
                    data_bytes);
                if (!HaveSystematicPacket[block_id])
                {
                    HaveSystematicPacket[block_id] = 1u;
                    ++CachedSystematicPacketCountValue;
                }
            }
        }
        catch (...)
        {
            ReceivedBlockStorage.resize(old_storage_size);
            ReceivedSlots.Erase(block_id);
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
        const uint32_t K = ProfileValue.BlockCount;
        const bool cache_enabled = SystematicPacketCache != nullptr;
        if (cache_enabled &&
            (SystematicPacketCacheSize != (size_t)MessageBytesValue ||
             HaveSystematicPacket.size() != K))
        {
            return Wirehair_Error;
        }
        std::vector<uint8_t> block(BlockBytesValue, 0u);
        uint8_t* output = static_cast<uint8_t*>(message_out);
        for (uint32_t block_id = 0; block_id < K; ++block_id)
        {
            const uint32_t bytes = PacketDataBytes(
                MessageBytesValue, BlockBytesValue, K, block_id);
            if (cache_enabled && HaveSystematicPacket[block_id])
            {
                std::memcpy(
                    output + (size_t)block_id * BlockBytesValue,
                    SystematicPacketCache.get() +
                        (size_t)block_id * BlockBytesValue,
                    bytes);
                continue;
            }
            if (!EvaluatePacketBlockForValidatedSystemWithRuntime(
                    SystemValue,
                    PacketConfigValue,
                    PacketRuntimeValue,
                    IntermediateBlockStorage.data(),
                    BlockBytesValue,
                    block_id,
                    block.data()))
            {
                return Wirehair_Error;
            }
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

void MessagePrecodeDecoder::ReleaseSystematicPacketCache() noexcept
{
    SystematicPacketCache.reset();
    SystematicPacketCacheSize = 0u;
    std::vector<uint8_t>().swap(HaveSystematicPacket);
    CachedSystematicPacketCountValue = 0u;
    OptionsValue.CacheReceivedSystematicPackets = false;
}

bool MessagePrecodeDecoder::HasSystematicPacketCache() const
{
    return SystematicPacketCache != nullptr;
}

size_t MessagePrecodeDecoder::SystematicPacketCacheBytes() const
{
    return SystematicPacketCacheSize + HaveSystematicPacket.size();
}

uint32_t MessagePrecodeDecoder::CachedSystematicPacketCount() const
{
    return CachedSystematicPacketCountValue;
}

bool MessagePrecodeDecoder::IsInitialized() const { return Initialized; }
bool MessagePrecodeDecoder::IsDecoded() const { return Decoded; }
uint32_t MessagePrecodeDecoder::ReceivedCount() const
{
    return Decoded ? ReceivedCountValue :
        (uint32_t)ReceivedSlots.Size();
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
