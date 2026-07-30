#include "WirehairV2Codec.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool Check(bool condition, const char* message)
{
    if (!condition) {
        std::fprintf(stderr, "%s\n", message);
        return false;
    }
    return true;
}

bool CheckPacketSlotTable()
{
    wirehair_v2::PacketSlotTable table;
    if (!Check(table.Initialize(37u, 37u),
            "packet slot table init failed") ||
        !Check(table.Size() == 0u && table.Capacity() >= 74u,
            "packet slot table initial shape mismatch") ||
        !Check(table.Insert(UINT32_MAX, 7u),
            "packet slot table rejected UINT32_MAX") ||
        !Check(table.Insert(0u, UINT32_MAX),
            "packet slot table rejected UINT32_MAX slot") ||
        !Check(!table.Insert(UINT32_MAX, 8u),
            "packet slot table accepted duplicate"))
    {
        return false;
    }
    uint32_t slot = 0u;
    if (!Check(table.Find(UINT32_MAX, &slot) && slot == 7u,
            "packet slot table lost maximum id") ||
        !Check(table.Find(0u, &slot) && slot == UINT32_MAX,
            "packet slot table lost maximum slot"))
    {
        return false;
    }

    if (!Check(table.Initialize(2u, 8u),
            "packet slot growth fixture init failed") ||
        !Check(table.Insert(11u, 1u) && table.Insert(22u, 2u),
            "packet slot growth fixture fill failed"))
    {
        return false;
    }
    const size_t before_duplicate_capacity = table.Capacity();
    if (!Check(!table.Insert(11u, 3u) &&
            table.Capacity() == before_duplicate_capacity,
            "half-full duplicate grew packet slot table"))
    {
        return false;
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    bool growth_oom = false;
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
    try {
        (void)table.Insert(33u, 3u);
    }
    catch (const std::bad_alloc&) {
        growth_oom = true;
    }
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    if (!Check(growth_oom && table.Size() == 2u &&
            table.Capacity() == before_duplicate_capacity &&
            !table.Find(33u),
            "packet slot growth OOM mutated table"))
    {
        return false;
    }
#endif
    if (!Check(table.Insert(33u, 3u) &&
            table.Capacity() > before_duplicate_capacity,
            "packet slot table did not grow transactionally"))
    {
        return false;
    }

    if (!Check(table.Initialize(37u, 37u),
            "packet slot oracle reinit failed") ||
        !Check(table.Insert(UINT32_MAX, 7u) &&
            table.Insert(0u, UINT32_MAX),
            "packet slot oracle seed failed"))
    {
        return false;
    }

    std::map<uint32_t, uint32_t> oracle;
    oracle[UINT32_MAX] = 7u;
    oracle[0u] = UINT32_MAX;
    uint32_t random = UINT32_C(0x8f31a25c);
    for (uint32_t step = 0u; step < 20000u; ++step)
    {
        random ^= random << 13;
        random ^= random >> 17;
        random ^= random << 5;
        const uint32_t id = random;
        if ((step % 3u) == 0u)
        {
            const bool expected = oracle.erase(id) != 0u;
            if (!Check(table.Erase(id) == expected,
                    "packet slot table erase disagreed with oracle"))
            {
                return false;
            }
        }
        else
        {
            const bool existed = oracle.find(id) != oracle.end();
            const bool room = oracle.size() < 37u;
            const bool inserted = table.Insert(id, step);
            if (!Check(inserted == (!existed && room),
                    "packet slot table insert disagreed with oracle"))
            {
                return false;
            }
            if (inserted) {
                oracle[id] = step;
            }
        }
        if ((step % 29u) == 0u)
        {
            for (const std::pair<const uint32_t, uint32_t>& entry : oracle)
            {
                if (!Check(
                        table.Find(entry.first, &slot) &&
                            slot == entry.second,
                        "packet slot table lookup disagreed with oracle"))
                {
                    return false;
                }
            }
        }
    }
    if (!Check(table.Size() == oracle.size(),
            "packet slot table size disagreed with oracle"))
    {
        return false;
    }
    if (!Check(table.Initialize(4u, 37u),
            "packet slot table reinit failed"))
    {
        return false;
    }
    for (uint32_t id = 0u; id < 37u; ++id) {
        if (!Check(table.Insert(id, id ^ UINT32_C(0xa5a5a5a5)),
                "packet slot table failed bounded fill"))
        {
            return false;
        }
    }
    if (!Check(!table.Insert(1000u, 1u),
            "packet slot table exceeded entry bound"))
    {
        return false;
    }
    for (uint32_t id = 0u; id < 37u; id += 2u) {
        if (!Check(table.Erase(id),
                "packet slot table failed clustered erase"))
        {
            return false;
        }
    }
    for (uint32_t id = 1u; id < 37u; id += 2u) {
        if (!Check(
                table.Find(id, &slot) &&
                    slot == (id ^ UINT32_C(0xa5a5a5a5)),
                "packet slot table erase broke probe chain"))
        {
            return false;
        }
    }
    if (!Check(table.Initialize(4u, 4u),
            "packet slot wraparound init failed"))
    {
        return false;
    }
    std::vector<uint32_t> wrapping_ids;
    for (uint32_t id = 0u; wrapping_ids.size() < 4u; ++id)
    {
        uint32_t hash = id;
        hash ^= hash >> 16;
        hash *= UINT32_C(0x7feb352d);
        hash ^= hash >> 15;
        hash *= UINT32_C(0x846ca68b);
        hash ^= hash >> 16;
        if ((hash & 7u) == 7u) {
            wrapping_ids.push_back(id);
        }
    }
    for (uint32_t i = 0u; i < wrapping_ids.size(); ++i) {
        if (!Check(table.Insert(wrapping_ids[i], i + 100u),
                "packet slot wraparound insert failed"))
        {
            return false;
        }
    }
    if (!Check(table.Erase(wrapping_ids[0]),
            "packet slot wraparound erase failed"))
    {
        return false;
    }
    for (uint32_t i = 1u; i < wrapping_ids.size(); ++i) {
        if (!Check(table.Find(wrapping_ids[i], &slot) && slot == i + 100u,
                "packet slot wraparound erase broke probe chain"))
        {
            return false;
        }
    }
    table.ClearAndRelease();
    return Check(
        table.Size() == 0u && table.Capacity() == 0u &&
            !table.Find(UINT32_MAX),
        "packet slot table clear retained state");
}

bool SameStats(
    const wirehair_v2::PrecodeSolveStats& a,
    const wirehair_v2::PrecodeSolveStats& b)
{
    return a.PacketRows == b.PacketRows &&
        a.PeeledColumns == b.PeeledColumns &&
        a.InactivatedColumns == b.InactivatedColumns &&
        a.ResidualRows == b.ResidualRows &&
        a.ResidualRank == b.ResidualRank &&
        a.BinaryResidualRank == b.BinaryResidualRank &&
        a.BinaryRowReferences == b.BinaryRowReferences &&
        a.BinaryRowStorageBytes == b.BinaryRowStorageBytes &&
        a.BinaryAdjacencyStorageBytes == b.BinaryAdjacencyStorageBytes &&
        a.BinaryRowStorageAllocations == b.BinaryRowStorageAllocations &&
        a.BinaryAdjacencyStorageAllocations ==
            b.BinaryAdjacencyStorageAllocations &&
        a.BlockXors == b.BlockXors &&
        a.BlockMulAdds == b.BlockMulAdds &&
        a.BlockCopies == b.BlockCopies &&
        a.BlockZeroFills == b.BlockZeroFills &&
        a.BlockAddSets == b.BlockAddSets &&
        a.BlockAddSetSources == b.BlockAddSetSources &&
        a.PeelAdjacencyVisits == b.PeelAdjacencyVisits &&
        a.PeelRowScanSteps == b.PeelRowScanSteps &&
        a.PeelHeapOperations == b.PeelHeapOperations &&
        a.PeelHeapCompactionRebuildColumnProbes ==
            b.PeelHeapCompactionRebuildColumnProbes &&
        a.PeelHeapCompactionHeapifyInputKeys ==
            b.PeelHeapCompactionHeapifyInputKeys &&
        a.ProjectionWordXors == b.ProjectionWordXors &&
        a.ResidualCoeffWordXors == b.ResidualCoeffWordXors &&
        a.ResidualCoeffByteOps == b.ResidualCoeffByteOps &&
        a.BuildNanoseconds == b.BuildNanoseconds &&
        a.PeelNanoseconds == b.PeelNanoseconds &&
        a.ProjectNanoseconds == b.ProjectNanoseconds &&
        a.ResidualNanoseconds == b.ResidualNanoseconds &&
        a.BackSubNanoseconds == b.BackSubNanoseconds &&
        a.PacketSeedAttempt == b.PacketSeedAttempt &&
        a.TinyMixedFastPathAcceptances ==
            b.TinyMixedFastPathAcceptances &&
        a.CertifiedPackedResidualUses ==
            b.CertifiedPackedResidualUses &&
        a.ResidualCoefficientStorageBytes ==
            b.ResidualCoefficientStorageBytes &&
        a.CertifiedPackedResumeMaterializations ==
            b.CertifiedPackedResumeMaterializations &&
        a.CertifiedFixedHQuotientUses ==
            b.CertifiedFixedHQuotientUses &&
        a.CertifiedFixedHCoefficientStorageBytes ==
            b.CertifiedFixedHCoefficientStorageBytes &&
        a.CertifiedSquareQuotientColumns ==
            b.CertifiedSquareQuotientColumns &&
        a.CertifiedHeavyRhsRowsBuilt ==
            b.CertifiedHeavyRhsRowsBuilt &&
        a.CertifiedLegacyHeavyRowsReplayed ==
            b.CertifiedLegacyHeavyRowsReplayed &&
        a.MixedJointSourceXors == b.MixedJointSourceXors &&
        a.MixedJointMarginalXors == b.MixedJointMarginalXors &&
        a.MixedJointMarginalCopies == b.MixedJointMarginalCopies &&
        a.MixedJointScratchBytes == b.MixedJointScratchBytes &&
        a.MixedJointActiveDeltas == b.MixedJointActiveDeltas &&
        a.MixedDualSourceColumns == b.MixedDualSourceColumns;
}

std::vector<uint8_t> MakeMessage(size_t bytes)
{
    std::vector<uint8_t> message(bytes);
    for (size_t i = 0; i < bytes; ++i) {
        message[i] = (uint8_t)(i * 37u + i / 11u + 19u);
    }
    return message;
}

struct PacketFixture
{
    uint32_t Id = 0u;
    std::vector<uint8_t> Data;
    uint32_t Bytes = 0u;
};

bool EncodePacket(
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    uint32_t block_bytes,
    uint32_t id,
    PacketFixture& packet)
{
    packet.Id = id;
    packet.Data.assign(block_bytes, 0u);
    packet.Bytes = 0u;
    return encoder.EncodeResult(
        id,
        packet.Data.data(),
        block_bytes,
        &packet.Bytes) == Wirehair_Success;
}

bool FindRepairRowCollisions(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    std::vector<uint32_t>& collision_ids,
    size_t required_pairs = 2u)
{
    struct SeenRow
    {
        uint32_t First = 0u;
        bool Included = false;
    };
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    std::map<std::vector<uint32_t>, SeenRow> seen;
    collision_ids.clear();
    if (required_pairs == 0u) {
        return true;
    }
    const uint32_t search_limit = K + 1000000u;
    for (uint32_t id = K; id < search_limit; ++id)
    {
        std::vector<uint32_t> row =
            wirehair_v2::GeneratePacketMatrixRow(K, P, id, config);
        std::sort(row.begin(), row.end());
        SeenRow candidate;
        candidate.First = id;
        const std::pair<
            std::map<std::vector<uint32_t>, SeenRow>::iterator, bool> found =
                seen.insert(std::make_pair(row, candidate));
        if (!found.second)
        {
            if (found.first->second.Included) {
                continue;
            }
            collision_ids.push_back(found.first->second.First);
            collision_ids.push_back(id);
            found.first->second.Included = true;
            if (collision_ids.size() / 2u >= required_pairs) {
                return true;
            }
        }
    }
    return false;
}

class DecoderPair
{
public:
    bool Initialize(
        uint64_t message_bytes,
        uint32_t block_bytes,
        const wirehair_v2::SeedProfile& profile)
    {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetDecoderIncrementalResumeEnabledForTesting(true);
#endif
        return Warm.InitializeResult(
                   message_bytes, block_bytes, &profile) == Wirehair_Success &&
            Cold.InitializeResult(
                   message_bytes, block_bytes, &profile) == Wirehair_Success;
    }

    bool DecodeSame(
        const PacketFixture& packet,
        WirehairResult expected,
        bool fail_solve_allocation = false)
    {
        WirehairResult actual = Wirehair_Error;
        if (!DecodeMatching(packet, actual, fail_solve_allocation)) {
            return false;
        }
        if (actual != expected)
        {
            std::fprintf(stderr,
                "decode result mismatch id=%u actual=%d expected=%d\n",
                packet.Id,
                (int)actual,
                (int)expected);
            return false;
        }
        return true;
    }

    bool DecodeMatching(
        const PacketFixture& packet,
        WirehairResult& actual,
        bool fail_solve_allocation = false)
    {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetDecoderIncrementalResumeEnabledForTesting(true);
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(
            fail_solve_allocation ? 0 : -1);
#endif
        const WirehairResult warm = Warm.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetDecoderIncrementalResumeEnabledForTesting(false);
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(
            fail_solve_allocation ? 0 : -1);
#endif
        const WirehairResult cold = Cold.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
        wirehair_v2::SetDecoderIncrementalResumeEnabledForTesting(true);
#endif
        if (warm != cold)
        {
            std::fprintf(stderr,
                "warm/cold decode mismatch id=%u warm=%d cold=%d\n",
                packet.Id, (int)warm, (int)cold);
            return false;
        }
        actual = warm;
        return true;
    }

    bool RecoverSame(
        uint64_t message_bytes,
        std::vector<uint8_t>& recovered)
    {
        std::vector<uint8_t> warm((size_t)message_bytes, 0u);
        std::vector<uint8_t> cold((size_t)message_bytes, 0u);
        if (Warm.RecoverResult(warm.data(), message_bytes) !=
                Wirehair_Success ||
            Cold.RecoverResult(cold.data(), message_bytes) !=
                Wirehair_Success ||
            warm != cold)
        {
            return false;
        }
        recovered.swap(warm);
        return true;
    }

    wirehair_v2::MessagePrecodeDecoder Warm;
    wirehair_v2::MessagePrecodeDecoder Cold;
};

bool FeedInitialSystematicHeavy(
    DecoderPair& decoders,
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    uint32_t K,
    uint32_t block_bytes,
    const std::vector<uint32_t>& collision_ids)
{
    PacketFixture packet;
    for (uint32_t id = 0u; id + collision_ids.size() < K; ++id)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet) ||
            !decoders.DecodeSame(packet, Wirehair_NeedMore))
        {
            return false;
        }
    }
    for (uint32_t id : collision_ids)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet) ||
            !decoders.DecodeSame(packet, Wirehair_NeedMore))
        {
            return false;
        }
    }
    return true;
}

bool CheckIncrementalDecoderParity()
{
#if !defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return true;
#else
    const uint32_t K = 1000u;
    const uint32_t block_bytes = 257u;
    // Keep K unchanged while making the final systematic packet partial.  The
    // missing-row sequence below appends K-1 through both checkpoint and cold
    // paths, covering zero-padding and exact recovery at the resume boundary.
    const uint64_t message_bytes =
        (uint64_t)K * block_bytes - block_bytes / 2u;
    const std::vector<uint8_t> message = MakeMessage((size_t)message_bytes);
    wirehair_v2::MessagePrecodeEncoder encoder;
    if (!Check(
            encoder.InitializeResult(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "resume parity encoder initialization failed"))
    {
        return false;
    }
    const wirehair_v2::SeedProfile profile = encoder.Profile();
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = profile.V2PacketPeelSeed;
    config.MixCount = profile.V2RecoveryMixCount;
    std::vector<uint32_t> collision_ids;
    if (!Check(
            FindRepairRowCollisions(
                encoder.BlockEncoder().System(),
                config,
                collision_ids),
            "unable to find deterministic repair-row collisions"))
    {
        return false;
    }
    const uint32_t collision_a = collision_ids[0];
    const uint32_t collision_b = collision_ids[1];
    PacketFixture first_collision;
    PacketFixture second_collision;
    if (!Check(
            EncodePacket(
                encoder, block_bytes, collision_a, first_collision) &&
                EncodePacket(
                    encoder, block_bytes, collision_b, second_collision) &&
                first_collision.Data == second_collision.Data,
            "equal repair rows did not encode equal bytes"))
    {
        return false;
    }

    // Systematic-heavy, rank-deficient, duplicate/conflict, and OOM retry.
    DecoderPair main_pair;
    if (!Check(
            main_pair.Initialize(message_bytes, block_bytes, profile),
            "resume parity decoder initialization failed") ||
        !FeedInitialSystematicHeavy(
            main_pair,
            encoder,
            K,
            block_bytes,
            collision_ids))
    {
        return false;
    }
    const size_t cold_receive_bytes =
        main_pair.Cold.ColdReceiveCapacityBytesForTesting();
    const size_t warm_resume_bytes =
        main_pair.Warm.IncrementalResumeBytesForTesting();
    if (!Check(
            main_pair.Warm.HasIncrementalResumeStateForTesting() &&
                !main_pair.Cold.HasIncrementalResumeStateForTesting(),
            "warm decoder did not retain the incremental checkpoint") ||
        !Check(
            main_pair.Warm.SolveStats().CertifiedPackedResidualUses == 0u &&
                main_pair.Warm.SolveStats().
                    CertifiedPackedResumeMaterializations == 0u,
            "narrow retained checkpoint did not use the byte basis") ||
        !Check(
            cold_receive_bytes != 0u &&
                warm_resume_bytes <=
                    cold_receive_bytes + cold_receive_bytes / 4u,
            "incremental checkpoint exceeded the 25% memory policy"))
    {
        return false;
    }

    // Checkpoint adoption allocates its one-packet retry slot before releasing
    // the cold receive buffer.  A failure at precisely that boundary must keep
    // all K packets and permit an identical retry to build the checkpoint.
    wirehair_v2::MessagePrecodeDecoder adoption_decoder;
    if (adoption_decoder.InitializeResult(
            message_bytes, block_bytes, &profile) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint32_t> initial_ids;
    initial_ids.reserve(K);
    for (uint32_t id = 0u; id + collision_ids.size() < K; ++id) {
        initial_ids.push_back(id);
    }
    initial_ids.insert(
        initial_ids.end(), collision_ids.begin(), collision_ids.end());
    PacketFixture adoption_packet;
    for (size_t i = 0u; i + 1u < initial_ids.size(); ++i)
    {
        if (!EncodePacket(
                encoder,
                block_bytes,
                initial_ids[i],
                adoption_packet) ||
            adoption_decoder.DecodeResult(
                adoption_packet.Id,
                adoption_packet.Data.data(),
                adoption_packet.Bytes) != Wirehair_NeedMore)
        {
            return false;
        }
    }
    if (!EncodePacket(
            encoder,
            block_bytes,
            initial_ids.back(),
            adoption_packet))
    {
        return false;
    }
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(1);
    const WirehairResult adoption_oom = adoption_decoder.DecodeResult(
        adoption_packet.Id,
        adoption_packet.Data.data(),
        adoption_packet.Bytes);
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    if (!Check(
            adoption_oom == Wirehair_OOM &&
                adoption_decoder.ReceivedCount() == K &&
                !adoption_decoder.HasIncrementalResumeStateForTesting() &&
                adoption_decoder.ColdReceiveCapacityBytesForTesting() != 0u &&
                adoption_decoder.DecodeResult(
                    adoption_packet.Id,
                    adoption_packet.Data.data(),
                    adoption_packet.Bytes) == Wirehair_NeedMore &&
                adoption_decoder.HasIncrementalResumeStateForTesting(),
            "checkpoint-slot OOM did not preserve a cold retry"))
    {
        return false;
    }

    if (!main_pair.DecodeSame(first_collision, Wirehair_NeedMore)) {
        return false;
    }
    PacketFixture conflict = first_collision;
    conflict.Data[0] ^= 1u;
    if (!main_pair.DecodeSame(conflict, Wirehair_InvalidInput)) {
        return false;
    }

    std::vector<PacketFixture> missing(collision_ids.size());
    const uint32_t first_missing = K - (uint32_t)missing.size();
    for (uint32_t i = 0u; i < (uint32_t)missing.size(); ++i) {
        if (!EncodePacket(
                encoder, block_bytes, first_missing + i, missing[i]))
        {
            return false;
        }
    }
    const wirehair_v2::PrecodeSolveStats warm_stats_before_oom =
        main_pair.Warm.SolveStats();
    const wirehair_v2::PrecodeSolveStats cold_stats_before_oom =
        main_pair.Cold.SolveStats();
    if (!main_pair.DecodeSame(missing[0], Wirehair_OOM, true) ||
        main_pair.Warm.ReceivedCount() != K + 1u ||
        main_pair.Cold.ReceivedCount() != K + 1u ||
        !SameStats(warm_stats_before_oom, main_pair.Warm.SolveStats()) ||
        !SameStats(cold_stats_before_oom, main_pair.Cold.SolveStats()))
    {
        return false;
    }
    WirehairResult main_result = Wirehair_NeedMore;
    size_t accepted_missing = 1u;
    for (size_t i = 1u;
         i < missing.size() && main_result == Wirehair_NeedMore;
         ++i)
    {
        // Supplying a different packet must first replay the pending OOM row;
        // otherwise the accepted id would silently lose its equation.
        if (!main_pair.DecodeMatching(missing[i], main_result)) {
            return false;
        }
        ++accepted_missing;
    }
    if (!Check(
            main_result == Wirehair_Success &&
                main_pair.Warm.SolveStats().PacketRows ==
                    main_pair.Cold.SolveStats().PacketRows &&
                main_pair.Warm.SolveStats().PacketRows ==
                    K + accepted_missing,
            "pending OOM equation was not replayed exactly once"))
    {
        return false;
    }
    std::vector<uint8_t> recovered;
    if (!Check(
            main_pair.RecoverSame(message_bytes, recovered) &&
                recovered == message,
            "warm/cold OOM retry recovery mismatch"))
    {
        return false;
    }

    PacketFixture overdetermined;
    if (!EncodePacket(
            encoder,
            block_bytes,
            collision_ids.back() + 1u,
            overdetermined) ||
        !main_pair.DecodeSame(overdetermined, Wirehair_Success))
    {
        return false;
    }
    overdetermined.Data[0] ^= 1u;
    if (!main_pair.DecodeSame(overdetermined, Wirehair_Error)) {
        return false;
    }

    // Exercise a terminal inconsistent append while the checkpoint is still
    // deficient.  The fourth repair id has the same coefficient row as the
    // third, so altering its RHS must match a cold solve's Error.  Retaining
    // just this terminal packet distinguishes an exact Error retry from a
    // conflicting duplicate after the bulk receive buffer has been released.
    DecoderPair terminal_pair;
    if (!terminal_pair.Initialize(message_bytes, block_bytes, profile)) {
        return false;
    }
    PacketFixture packet;
    for (uint32_t id = 0u; id + 3u < K; ++id)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet) ||
            !terminal_pair.DecodeSame(packet, Wirehair_NeedMore))
        {
            return false;
        }
    }
    for (size_t i = 0u; i < 3u; ++i)
    {
        if (!EncodePacket(
                encoder, block_bytes, collision_ids[i], packet) ||
            !terminal_pair.DecodeSame(packet, Wirehair_NeedMore))
        {
            return false;
        }
    }
    PacketFixture terminal_authentic;
    if (!EncodePacket(
            encoder,
            block_bytes,
            collision_ids[3],
            terminal_authentic))
    {
        return false;
    }
    PacketFixture terminal_corrupt = terminal_authentic;
    terminal_corrupt.Data[0] ^= 1u;
    if (!terminal_pair.DecodeSame(terminal_corrupt, Wirehair_Error) ||
        terminal_pair.Warm.SolveStats().PacketRows != K + 1u ||
        terminal_pair.Cold.SolveStats().PacketRows != K + 1u)
    {
        return false;
    }
    const wirehair_v2::PrecodeSolveStats warm_terminal_stats =
        terminal_pair.Warm.SolveStats();
    const wirehair_v2::PrecodeSolveStats cold_terminal_stats =
        terminal_pair.Cold.SolveStats();
    if (!terminal_pair.DecodeSame(terminal_corrupt, Wirehair_Error) ||
        !terminal_pair.DecodeSame(
            terminal_authentic, Wirehair_InvalidInput) ||
        !terminal_pair.DecodeSame(overdetermined, Wirehair_Error) ||
        !SameStats(warm_terminal_stats, terminal_pair.Warm.SolveStats()) ||
        !SameStats(cold_terminal_stats, terminal_pair.Cold.SolveStats()))
    {
        return false;
    }

    // An altered independent row can define a different valid message.  Warm
    // and cold solves must still produce identical bytes and later reject the
    // authentic overdetermining equation with the same error.
    DecoderPair corrupt_pair;
    if (!corrupt_pair.Initialize(message_bytes, block_bytes, profile) ||
        !FeedInitialSystematicHeavy(
            corrupt_pair,
            encoder,
            K,
            block_bytes,
            collision_ids))
    {
        return false;
    }
    PacketFixture altered = missing[0];
    altered.Data[0] ^= 1u;
    if (!corrupt_pair.DecodeSame(altered, Wirehair_NeedMore))
    {
        return false;
    }
    WirehairResult corrupt_result = Wirehair_NeedMore;
    for (size_t i = 1u;
         i < missing.size() && corrupt_result == Wirehair_NeedMore;
         ++i)
    {
        if (!corrupt_pair.DecodeMatching(missing[i], corrupt_result)) {
            return false;
        }
    }
    if (corrupt_result != Wirehair_Success) {
        return false;
    }
    recovered.clear();
    if (!Check(
            corrupt_pair.RecoverSame(message_bytes, recovered) &&
                recovered != message,
            "corrupted warm/cold stream did not produce matching altered bytes") ||
        !corrupt_pair.DecodeSame(missing[0], Wirehair_Error))
    {
        return false;
    }

    // Repair-only input uses distinct public ids but contains one exact matrix
    // collision.  The first K equations are therefore provably deficient.
    DecoderPair repair_pair;
    if (!repair_pair.Initialize(message_bytes, block_bytes, profile)) {
        return false;
    }
    std::vector<uint32_t> repair_ids;
    repair_ids.reserve(K);
    for (uint32_t id = K; repair_ids.size() + 2u < K; ++id) {
        if (id != collision_a && id != collision_b) {
            repair_ids.push_back(id);
        }
    }
    repair_ids.push_back(collision_a);
    repair_ids.push_back(collision_b);
    for (uint32_t id : repair_ids)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet) ||
            !repair_pair.DecodeSame(packet, Wirehair_NeedMore))
        {
            return false;
        }
    }
    if (!Check(
            repair_pair.Warm.HasIncrementalResumeStateForTesting(),
            "repair-only stream did not create a resume checkpoint"))
    {
        return false;
    }
    WirehairResult repair_result = Wirehair_NeedMore;
    for (uint32_t id = collision_b + 1000000u;
         repair_result == Wirehair_NeedMore &&
             id < collision_b + 1000064u;
         ++id)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet)) {
            return false;
        }
        wirehair_v2::SetDecoderIncrementalResumeEnabledForTesting(true);
        const WirehairResult warm = repair_pair.Warm.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
        wirehair_v2::SetDecoderIncrementalResumeEnabledForTesting(false);
        const WirehairResult cold = repair_pair.Cold.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
        wirehair_v2::SetDecoderIncrementalResumeEnabledForTesting(true);
        if (warm != cold ||
            (warm != Wirehair_NeedMore && warm != Wirehair_Success))
        {
            return false;
        }
        repair_result = warm;
    }
    recovered.clear();
    if (!Check(
            repair_result == Wirehair_Success &&
                repair_pair.RecoverSame(message_bytes, recovered) &&
                recovered == message,
            "repair-only warm/cold recovery mismatch"))
    {
        return false;
    }

    // With one-byte blocks the projection/pivot metadata is intentionally too
    // large relative to the receive payload buffer.  The decoder must retain
    // the cold packet store and continue to produce the same results.
    const uint32_t small_k = 64u;
    const uint32_t small_block_bytes = 1u;
    const uint64_t small_message_bytes = small_k;
    const std::vector<uint8_t> small_message =
        MakeMessage((size_t)small_message_bytes);
    wirehair_v2::MessagePrecodeEncoder small_encoder;
    if (small_encoder.InitializeResult(
            small_message.data(),
            small_message_bytes,
            small_block_bytes) != Wirehair_Success)
    {
        return false;
    }
    wirehair_v2::PacketRowConfig small_config;
    small_config.PeelSeed = small_encoder.Profile().V2PacketPeelSeed;
    small_config.MixCount = small_encoder.Profile().V2RecoveryMixCount;
    std::vector<uint32_t> small_collisions;
    if (!FindRepairRowCollisions(
            small_encoder.BlockEncoder().System(),
            small_config,
            small_collisions))
    {
        return false;
    }
    DecoderPair fallback_pair;
    if (!fallback_pair.Initialize(
            small_message_bytes,
            small_block_bytes,
            small_encoder.Profile()) ||
        !FeedInitialSystematicHeavy(
            fallback_pair,
            small_encoder,
            small_k,
            small_block_bytes,
            small_collisions) ||
        !Check(
            !fallback_pair.Warm.HasIncrementalResumeStateForTesting() &&
                fallback_pair.Warm.ColdReceiveCapacityBytesForTesting() != 0u &&
                fallback_pair.Warm.SolveStats().
                    CertifiedPackedResidualUses == 1u &&
                fallback_pair.Warm.SolveStats().
                    CertifiedPackedResumeMaterializations == 0u,
            "oversized checkpoint did not skip packed materialization"))
    {
        return false;
    }
    WirehairResult fallback_result = Wirehair_NeedMore;
    const uint32_t small_first_missing =
        small_k - (uint32_t)small_collisions.size();
    for (uint32_t id = small_first_missing;
         id < small_k && fallback_result == Wirehair_NeedMore;
         ++id)
    {
        if (!EncodePacket(
                small_encoder, small_block_bytes, id, packet) ||
            !fallback_pair.DecodeMatching(packet, fallback_result))
        {
            return false;
        }
    }
    recovered.clear();
    if (!Check(
            fallback_result == Wirehair_Success &&
                fallback_pair.RecoverSame(
                    small_message_bytes, recovered) &&
                recovered == small_message,
            "cold memory-policy fallback changed recovery"))
    {
        return false;
    }

    std::printf(
        "incremental decoder parity: collision=%u/%u checkpoint=%zu "
        "cold=%zu\n",
        collision_a,
        collision_b,
        warm_resume_bytes,
        cold_receive_bytes);
    return true;
#endif
}

bool CheckColdDuplicateSlotLookup()
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 37u;
    const uint32_t tail_bytes = 11u;
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    const std::vector<uint8_t> message = MakeMessage((size_t)message_bytes);
    wirehair_v2::MessagePrecodeEncoder encoder;
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!Check(
            encoder.InitializeResult(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success &&
            decoder.InitializeResult(message_bytes, block_bytes) ==
                Wirehair_Success,
            "cold duplicate slot fixture initialization failed"))
    {
        return false;
    }

    PacketFixture packet;
    for (uint32_t id = 0u; id + 3u < K; ++id)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet) ||
            decoder.DecodeResult(
                packet.Id, packet.Data.data(), packet.Bytes) !=
                Wirehair_NeedMore)
        {
            return false;
        }
    }

    PacketFixture high;
    PacketFixture partial;
    if (!EncodePacket(encoder, block_bytes, UINT32_MAX, high) ||
        decoder.DecodeResult(high.Id, high.Data.data(), high.Bytes) !=
            Wirehair_NeedMore ||
        !EncodePacket(encoder, block_bytes, K - 1u, partial) ||
        partial.Bytes != tail_bytes ||
        decoder.DecodeResult(
            partial.Id, partial.Data.data(), partial.Bytes) !=
            Wirehair_NeedMore ||
        decoder.ReceivedCount() != K - 1u)
    {
        std::fprintf(stderr, "cold duplicate slot fixture feed failed\n");
        return false;
    }

    if (!Check(
            decoder.DecodeResult(
                high.Id, high.Data.data(), high.Bytes) == Wirehair_NeedMore &&
            decoder.DecodeResult(
                partial.Id, partial.Data.data(), partial.Bytes) ==
                Wirehair_NeedMore &&
            decoder.ReceivedCount() == K - 1u,
            "cold duplicate slot lookup changed accepted state"))
    {
        return false;
    }

    PacketFixture high_conflict = high;
    PacketFixture partial_conflict = partial;
    high_conflict.Data[high_conflict.Bytes - 1u] ^= 0x80u;
    partial_conflict.Data[partial_conflict.Bytes - 1u] ^= 0x40u;
    if (!Check(
            decoder.DecodeResult(
                high_conflict.Id,
                high_conflict.Data.data(),
                high_conflict.Bytes) == Wirehair_InvalidInput &&
            decoder.DecodeResult(
                partial_conflict.Id,
                partial_conflict.Data.data(),
                partial_conflict.Bytes) == Wirehair_InvalidInput &&
            decoder.ReceivedCount() == K - 1u,
            "cold duplicate slot conflict was not rejected"))
    {
        return false;
    }

    std::printf("cold duplicate slot lookup: PASS\n");
    return true;
}

bool InitializeCompletedRecoveryFixture(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t tail_bytes,
    bool cache_systematic,
    bool omit_final_systematic,
    wirehair_v2::CompletionField completion,
    std::vector<uint8_t>& message,
    wirehair_v2::MessagePrecodeEncoder& encoder,
    wirehair_v2::MessagePrecodeDecoder& decoder)
{
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    message = MakeMessage((size_t)message_bytes);
    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.CacheReceivedSystematicPackets = cache_systematic;
    options.Completion = completion;
    if (encoder.InitializeResult(
            message.data(), message_bytes, block_bytes, nullptr, &options) !=
                Wirehair_Success ||
        decoder.InitializeResult(
            message_bytes, block_bytes, &encoder.Profile(), &options) !=
                Wirehair_Success)
    {
        return false;
    }

    WirehairResult result = Wirehair_NeedMore;
    PacketFixture packet;
    for (uint32_t id = 0u; id < K; ++id)
    {
        if (omit_final_systematic && id + 1u == K) {
            continue;
        }
        if (!EncodePacket(encoder, block_bytes, id, packet)) {
            return false;
        }
        result = decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
        if (result != Wirehair_NeedMore && result != Wirehair_Success) {
            return false;
        }
    }
    for (uint32_t repair = 0u;
         repair < K + 512u && result == Wirehair_NeedMore;
         ++repair)
    {
        if (!EncodePacket(encoder, block_bytes, K + repair, packet)) {
            return false;
        }
        result = decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
    }
    const uint32_t expected_cached = !cache_systematic ? 0u :
        K - (omit_final_systematic ? 1u : 0u);
    return result == Wirehair_Success &&
        decoder.CachedSystematicPacketCount() == expected_cached;
}

bool CheckDirectRecoveryOutput()
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 38u;
    struct RecoveryCase
    {
        uint32_t TailBytes;
        bool CacheSystematic;
        bool OmitFinalSystematic;
        bool ReleaseCacheBeforeRecover;
        wirehair_v2::CompletionField Completion;
        bool ExpectScratchAllocation;
        const char* Name;
    };
    const RecoveryCase cases[] = {
        {block_bytes, false, false, false,
            wirehair_v2::CompletionField::GF256,
            false, "uncached exact"},
        {11u, false, false, false,
            wirehair_v2::CompletionField::GF256,
            true, "uncached partial"},
        {block_bytes, true, false, false,
            wirehair_v2::CompletionField::GF256,
            false, "cached exact"},
        {11u, true, false, false,
            wirehair_v2::CompletionField::GF256,
            false, "cached partial"},
        {11u, true, true, false,
            wirehair_v2::CompletionField::GF256,
            true, "cached missing partial"},
        {11u, true, false, true,
            wirehair_v2::CompletionField::GF256,
            true, "cached partial after release"},
        {11u, true, false, false,
            wirehair_v2::CompletionField::MixedGF256GF16,
            false, "mixed cached partial"},
        {11u, true, true, false,
            wirehair_v2::CompletionField::MixedGF256GF16,
            true, "mixed cached missing partial"}
    };

    for (const RecoveryCase& c : cases)
    {
        std::vector<uint8_t> message;
        wirehair_v2::MessagePrecodeEncoder encoder;
        wirehair_v2::MessagePrecodeDecoder decoder;
        if (!InitializeCompletedRecoveryFixture(
                K, block_bytes, c.TailBytes, c.CacheSystematic,
                c.OmitFinalSystematic, c.Completion,
                message, encoder, decoder))
        {
            std::fprintf(stderr,
                "direct recovery fixture failed: %s\n", c.Name);
            return false;
        }

        const wirehair_v2::PrecodeParams& params = decoder.System().Params;
        const size_t intermediate_bytes = (size_t)(
            (uint64_t)params.BlockCount + params.Staircase +
            params.DenseRows + params.HeavyRows) * block_bytes;
        const uint8_t* const intermediate = decoder.IntermediateBlocks();
        if (!intermediate || intermediate_bytes <= message.size()) {
            std::fprintf(stderr,
                "direct recovery intermediate shape failed: %s\n", c.Name);
            return false;
        }
        const std::vector<uint8_t> intermediate_before(
            intermediate, intermediate + intermediate_bytes);
        const size_t overlap_offsets[] = {
            0u, 1u, intermediate_bytes - 1u
        };
        WirehairResult result = Wirehair_Error;
        for (size_t overlap_offset : overlap_offsets)
        {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
#endif
            result = decoder.RecoverResult(
                const_cast<uint8_t*>(intermediate) + overlap_offset,
                message.size());
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
#endif
            if (!Check(result == Wirehair_InvalidInput &&
                    std::equal(
                        intermediate_before.begin(), intermediate_before.end(),
                        intermediate),
                    "overlapping recovery output was not rejected/no-write"))
            {
                return false;
            }
        }

        std::vector<uint8_t> recovered(message.size(), 0xa7u);
        const std::vector<uint8_t> before = recovered;
        result = decoder.RecoverResult(
            recovered.data(), recovered.size() - 1u);
        if (!Check(result == Wirehair_InvalidInput && recovered == before,
                "wrong-size recovery was not rejected/no-write"))
        {
            return false;
        }
        if (c.ReleaseCacheBeforeRecover) {
            decoder.ReleaseSystematicPacketCache();
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
#endif
        result = decoder.RecoverResult(recovered.data(), recovered.size());
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
        if (c.ExpectScratchAllocation)
        {
            if (!Check(result == Wirehair_OOM && recovered == before,
                    "partial recovery scratch OOM was not no-write"))
            {
                return false;
            }
            wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(1);
            result = decoder.RecoverResult(
                recovered.data(), recovered.size());
            wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
        }
#endif
        if (!Check(result == Wirehair_Success && recovered == message,
                "direct recovery payload/allocation contract mismatch"))
        {
            std::fprintf(stderr, "direct recovery case: %s\n", c.Name);
            return false;
        }
    }
    return true;
}

bool CheckSystematicRecoverCache()
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 37u;
    const uint32_t tail_bytes = 11u;
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    const std::vector<uint8_t> message = MakeMessage((size_t)message_bytes);
    wirehair_v2::MessagePrecodeEncoder encoder;
    if (!Check(
            encoder.InitializeResult(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "recover-cache encoder initialization failed"))
    {
        return false;
    }
    const wirehair_v2::SeedProfile profile = encoder.Profile();
    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.CacheReceivedSystematicPackets = true;

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    for (int64_t countdown = 4; countdown <= 5; ++countdown)
    {
        wirehair_v2::MessagePrecodeDecoder failed;
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(
            countdown);
        const WirehairResult result = failed.InitializeResult(
            message_bytes, block_bytes, &profile, &options);
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
        if (!Check(result == Wirehair_OOM && !failed.IsInitialized() &&
                !failed.HasSystematicPacketCache(),
                "recover-cache initialization OOM was not transactional"))
        {
            return false;
        }
    }
#endif

    wirehair_v2::MessagePrecodeDecoder defaults;
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!Check(
            defaults.InitializeResult(
                message_bytes, block_bytes, &profile) == Wirehair_Success &&
            !defaults.HasSystematicPacketCache() &&
            defaults.SystematicPacketCacheBytes() == 0u,
            "default decoder allocated systematic cache") ||
        !Check(
            decoder.InitializeResult(
                message_bytes, block_bytes, &profile, &options) ==
                Wirehair_Success &&
            decoder.HasSystematicPacketCache() &&
            decoder.SystematicPacketCacheBytes() == message_bytes + K &&
            decoder.CachedSystematicPacketCount() == 0u,
            "enabled decoder cache shape mismatch"))
    {
        return false;
    }

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const size_t cache_bytes_before_reinit =
        decoder.SystematicPacketCacheBytes();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(4);
    const WirehairResult reinit_oom = decoder.InitializeResult(
        message_bytes, block_bytes, &profile, &options);
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    if (!Check(reinit_oom == Wirehair_OOM && decoder.IsInitialized() &&
            decoder.HasSystematicPacketCache() &&
            decoder.SystematicPacketCacheBytes() ==
                cache_bytes_before_reinit &&
            decoder.CachedSystematicPacketCount() == 0u,
            "recover-cache failed reinit changed active decoder"))
    {
        return false;
    }
#endif

    PacketFixture packet;
    WirehairResult result = Wirehair_NeedMore;
    uint32_t expected_cached = 0u;
    for (uint32_t id = 0u; id < K; ++id)
    {
        if ((id % 5u) == 0u) {
            continue;
        }
        if (!EncodePacket(encoder, block_bytes, id, packet)) {
            return false;
        }
        result = decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
        if (!Check(result == Wirehair_NeedMore,
                "recover-cache systematic feed ended early"))
        {
            return false;
        }
        ++expected_cached;
    }
    if (!EncodePacket(encoder, block_bytes, K - 1u, packet) ||
        !Check(packet.Bytes == tail_bytes &&
            decoder.DecodeResult(
                packet.Id, packet.Data.data(), packet.Bytes) ==
                Wirehair_NeedMore &&
            decoder.CachedSystematicPacketCount() == expected_cached,
            "recover-cache duplicate/partial handling changed cache count"))
    {
        return false;
    }
    packet.Data[packet.Bytes - 1u] ^= 0x80u;
    if (!Check(decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes) ==
            Wirehair_InvalidInput &&
            decoder.CachedSystematicPacketCount() == expected_cached,
            "recover-cache conflicting duplicate mutated cache"))
    {
        return false;
    }
    packet.Data[packet.Bytes - 1u] ^= 0x80u;
    for (uint32_t repair = 0u;
         repair < K + 512u && result == Wirehair_NeedMore;
         ++repair)
    {
        if (!EncodePacket(encoder, block_bytes, K + repair, packet)) {
            return false;
        }
        result = decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
    }
    std::vector<uint8_t> recovered((size_t)message_bytes, 0u);
    if (!Check(result == Wirehair_Success &&
            decoder.CachedSystematicPacketCount() == expected_cached,
            "recover-cache decode did not preserve cache count") ||
        !Check(decoder.RecoverResult(
                recovered.data(), message_bytes) == Wirehair_Success &&
            recovered == message,
            "recover-cache partial-tail recovery mismatch"))
    {
        return false;
    }

    decoder.ReleaseSystematicPacketCache();
    decoder.ReleaseSystematicPacketCache();
    std::fill(recovered.begin(), recovered.end(), uint8_t{0});
    if (!Check(!decoder.HasSystematicPacketCache() &&
            decoder.SystematicPacketCacheBytes() == 0u &&
            decoder.CachedSystematicPacketCount() == 0u &&
            !decoder.Options().CacheReceivedSystematicPackets,
            "recover-cache release was not idempotent") ||
        !Check(decoder.RecoverResult(
                recovered.data(), message_bytes) == Wirehair_Success &&
            recovered == message,
            "recovery after cache release changed payload"))
    {
        return false;
    }
    wirehair_v2::Codec facade;
    if (!Check(facade.ReleasePrecodeDecoderSystematicCache() ==
            Wirehair_InvalidInput,
            "recover-cache facade release accepted empty codec") ||
        !Check(facade.InitializePrecodeDecoder(
                message_bytes, block_bytes, &profile, &options) ==
                Wirehair_Success &&
            facade.ReleasePrecodeDecoderSystematicCache() ==
                Wirehair_Success &&
            facade.ReleasePrecodeDecoderSystematicCache() ==
                Wirehair_Success,
            "recover-cache facade release was not idempotent"))
    {
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder released_before_recover;
    if (!Check(released_before_recover.InitializeResult(
            message_bytes, block_bytes, &profile, &options) ==
            Wirehair_Success,
            "release-before-recover decoder init failed"))
    {
        return false;
    }
    result = Wirehair_NeedMore;
    for (uint32_t id = 0u; id < K; ++id)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet)) {
            return false;
        }
        result = released_before_recover.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
    }
    released_before_recover.ReleaseSystematicPacketCache();
    std::fill(recovered.begin(), recovered.end(), uint8_t{0});
    if (!Check(result == Wirehair_Success &&
            released_before_recover.RecoverResult(
                recovered.data(), message_bytes) == Wirehair_Success &&
            recovered == message,
            "release-before-recover fallback changed payload"))
    {
        return false;
    }
    wirehair_v2::MessagePrecodeDecoder repair_only;
    if (!Check(repair_only.InitializeResult(
            message_bytes, block_bytes, &profile, &options) ==
            Wirehair_Success,
            "recover-cache repair-only decoder init failed"))
    {
        return false;
    }
    result = Wirehair_NeedMore;
    for (uint32_t repair = 0u;
         repair < K + 512u && result == Wirehair_NeedMore;
         ++repair)
    {
        if (!EncodePacket(
                encoder, block_bytes, K + repair, packet))
        {
            return false;
        }
        result = repair_only.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes);
    }
    std::fill(recovered.begin(), recovered.end(), uint8_t{0});
    if (!Check(result == Wirehair_Success &&
            repair_only.CachedSystematicPacketCount() == 0u &&
            repair_only.RecoverResult(
                recovered.data(), message_bytes) == Wirehair_Success &&
            recovered == message,
            "recover-cache repair-only fallback changed payload"))
    {
        return false;
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    // Force the K-th accepted equation to duplicate another repair row.  The
    // deficient cold solve must publish a resume checkpoint; only then append
    // a partial systematic packet and verify the resume receive branch caches
    // its exact bytes before completing recovery.
    const uint32_t resume_K = 1000u;
    const uint32_t resume_bb = 257u;
    const uint64_t resume_bytes =
        (uint64_t)resume_K * resume_bb - resume_bb / 2u;
    const std::vector<uint8_t> resume_message =
        MakeMessage((size_t)resume_bytes);
    wirehair_v2::MessagePrecodeEncoder resume_encoder;
    if (!Check(resume_encoder.InitializeResult(
            resume_message.data(), resume_bytes, resume_bb) ==
            Wirehair_Success,
            "resume recover-cache encoder init failed"))
    {
        return false;
    }
    wirehair_v2::PacketRowConfig resume_config;
    resume_config.PeelSeed = resume_encoder.Profile().V2PacketPeelSeed;
    resume_config.MixCount = resume_encoder.Profile().V2RecoveryMixCount;
    std::vector<uint32_t> collision_ids;
    if (!Check(FindRepairRowCollisions(
            resume_encoder.BlockEncoder().System(),
            resume_config,
            collision_ids),
            "resume recover-cache collision search failed"))
    {
        return false;
    }
    wirehair_v2::MessagePrecodeDecoder resume_decoder;
    if (!Check(resume_decoder.InitializeResult(
            resume_bytes,
            resume_bb,
            &resume_encoder.Profile(),
            &options) == Wirehair_Success,
            "resume recover-cache decoder init failed"))
    {
        return false;
    }
    PacketFixture resume_packet;
    for (uint32_t id = 0u; id + 2u < resume_K; ++id)
    {
        if (!EncodePacket(resume_encoder, resume_bb, id, resume_packet) ||
            resume_decoder.DecodeResult(
                resume_packet.Id,
                resume_packet.Data.data(),
                resume_packet.Bytes) != Wirehair_NeedMore)
        {
            return false;
        }
    }
    for (uint32_t i = 0u; i < 2u; ++i)
    {
        if (!EncodePacket(
                resume_encoder, resume_bb, collision_ids[i], resume_packet) ||
            resume_decoder.DecodeResult(
                resume_packet.Id,
                resume_packet.Data.data(),
                resume_packet.Bytes) != Wirehair_NeedMore)
        {
            return false;
        }
    }
    if (!Check(resume_decoder.HasIncrementalResumeStateForTesting() &&
            resume_decoder.CachedSystematicPacketCount() == resume_K - 2u,
            "recover-cache fixture did not enter resume state"))
    {
        return false;
    }
    if (!EncodePacket(
            resume_encoder, resume_bb, resume_K - 1u, resume_packet) ||
        !Check(resume_packet.Bytes == resume_bb - resume_bb / 2u,
            "resume recover-cache partial packet size mismatch"))
    {
        return false;
    }
    WirehairResult resume_result = resume_decoder.DecodeResult(
        resume_packet.Id, resume_packet.Data.data(), resume_packet.Bytes);
    for (uint32_t repair = 0u;
         repair < 512u && resume_result == Wirehair_NeedMore;
         ++repair)
    {
        const uint32_t id = resume_K + 1000000u + repair;
        if (!EncodePacket(resume_encoder, resume_bb, id, resume_packet)) {
            return false;
        }
        resume_result = resume_decoder.DecodeResult(
            resume_packet.Id, resume_packet.Data.data(), resume_packet.Bytes);
    }
    std::vector<uint8_t> resume_recovered((size_t)resume_bytes, 0u);
    if (!Check(resume_result == Wirehair_Success &&
            resume_decoder.CachedSystematicPacketCount() == resume_K - 1u &&
            resume_decoder.RecoverResult(
                resume_recovered.data(), resume_bytes) == Wirehair_Success &&
            resume_recovered == resume_message,
            "systematic append after resume was not cached exactly"))
    {
        return false;
    }
#endif
    return true;
}

bool RunDecodeCase(
    uint32_t block_count,
    uint32_t block_bytes,
    bool recovery_only,
    bool partial_final,
    bool keep_all_sources,
    uint32_t loss_modulus)
{
    const uint64_t full_bytes = (uint64_t)block_count * block_bytes;
    const uint64_t message_bytes = partial_final ?
        full_bytes - (block_bytes > 1u ? block_bytes / 2u : 0u) :
        full_bytes;
    std::vector<uint8_t> message = MakeMessage((size_t)message_bytes);
    std::vector<uint8_t> recovered((size_t)message_bytes, 0u);
    std::vector<uint8_t> block(block_bytes, 0u);

    wirehair_v2::Codec encoder;
    wirehair_v2::Codec decoder;
    if (!Check(
            encoder.InitializePrecodeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "precode encoder initialization failed") ||
        !Check(
            decoder.InitializePrecodeDecoder(
                message_bytes, block_bytes) ==
                Wirehair_Success,
            "precode decoder initialization failed"))
    {
        return false;
    }

    WirehairResult decode_result = Wirehair_NeedMore;
    if (recovery_only)
    {
        const uint32_t shuffled_ids[] = {
            UINT32_MAX, block_count + 17u, block_count + 3u
        };
        for (const uint32_t id : shuffled_ids)
        {
            uint32_t written = 0u;
            if (!Check(
                    encoder.Encode(
                        id, block.data(), block_bytes, &written) ==
                        Wirehair_Success && written == block_bytes,
                    "shuffled recovery encode failed"))
            {
                return false;
            }
            decode_result = decoder.Decode(id, block.data(), written);
            if (!Check(
                    decode_result == Wirehair_NeedMore,
                    "shuffled recovery decode ended early"))
            {
                return false;
            }
        }
        uint32_t written = 0u;
        if (!Check(
                encoder.Encode(
                    shuffled_ids[0], block.data(), block_bytes, &written) ==
                    Wirehair_Success,
                "duplicate high-id recovery encode failed") ||
            !Check(
                decoder.Decode(shuffled_ids[0], block.data(), written) ==
                    Wirehair_NeedMore,
                "duplicate high-id recovery packet was not ignored"))
        {
            return false;
        }
    }
    if (!recovery_only)
    {
        for (uint32_t id = 0; id < block_count; ++id)
        {
            if (!keep_all_sources && id % loss_modulus == 0u) {
                continue;
            }
            uint32_t written = 0u;
            if (!Check(
                    encoder.Encode(
                        id, block.data(), block_bytes, &written) ==
                        Wirehair_Success,
                    "systematic encode failed"))
            {
                return false;
            }
            decode_result = decoder.Decode(id, block.data(), written);
            const WirehairResult expected =
                keep_all_sources && id + 1u == block_count ?
                    Wirehair_Success : Wirehair_NeedMore;
            if (!Check(
                    decode_result == expected,
                    "systematic decode completed at the wrong point"))
            {
                return false;
            }
        }

        uint32_t written = 0u;
        if (!keep_all_sources &&
            (!Check(
                 encoder.Encode(
                     1u, block.data(), block_bytes, &written) ==
                     Wirehair_Success,
                 "duplicate source encode failed") ||
             !Check(
                 decoder.Decode(1u, block.data(), written) ==
                     Wirehair_NeedMore,
                 "duplicate packet was not ignored")))
        {
            return false;
        }
    }

    for (uint32_t recovery = 0;
        recovery < block_count + 512u &&
            decode_result == Wirehair_NeedMore;
        ++recovery)
    {
        const uint32_t id = block_count + recovery;
        uint32_t written = 0u;
        if (!Check(
                encoder.Encode(
                    id, block.data(), block_bytes, &written) ==
                    Wirehair_Success &&
                    written == block_bytes,
                "recovery encode failed"))
        {
            return false;
        }
        decode_result = decoder.Decode(id, block.data(), written);
    }

    if (!Check(
            decode_result == Wirehair_Success,
            "precode decode did not reach full rank") ||
        !Check(
            decoder.Recover(recovered.data(), message_bytes) ==
                Wirehair_Success,
            "precode recovery failed") ||
        !Check(
            recovered == message,
            "precode recovered message mismatch") ||
        !Check(
            decoder.Encode(0u, block.data(), block_bytes, &block_count) ==
                Wirehair_InvalidInput,
            "decoder mode unexpectedly allowed Encode"))
    {
        return false;
    }
    return true;
}

bool CheckInvalidAndOom()
{
    wirehair_v2::MessagePrecodeEncoderOptions options;
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!Check(
            decoder.InitializeResult(1024u, 16u, nullptr, &options) ==
                Wirehair_Success,
            "direct decoder initialization failed") ||
        !Check(
            decoder.InitializeResult(UINT64_MAX, 16u, nullptr, &options) ==
                Wirehair_InvalidInput && decoder.IsInitialized(),
            "decoder overflow input was not rejected transactionally") ||
        !Check(
            decoder.InitializeResult(
                UINT64_C(0x100000000), UINT32_C(0x80000000),
                nullptr, &options) == Wirehair_InvalidInput &&
                decoder.IsInitialized(),
            "decoder oversized block input was not rejected transactionally") ||
        !Check(
            decoder.RecoverResult(nullptr, 1024u) == Wirehair_InvalidInput,
            "null recovery output was accepted"))
    {
        return false;
    }
    uint8_t block[16] = {};
    if (!Check(
            decoder.DecodeResult(0u, nullptr, 16u) == Wirehair_InvalidInput,
            "null decode input was accepted") ||
        !Check(
            decoder.DecodeResult(0u, block, 15u) == Wirehair_InvalidInput,
            "short non-final source packet was accepted"))
    {
        return false;
    }

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t K = 64u;
    const uint32_t solve_block_bytes = 37u;
    std::vector<uint8_t> source((size_t)K * solve_block_bytes, 0u);
    wirehair_v2::MessagePrecodeDecoder solve_oom_decoder;
    if (!Check(
            solve_oom_decoder.InitializeResult(
                (uint64_t)K * solve_block_bytes,
                solve_block_bytes,
                nullptr,
                &options) == Wirehair_Success,
            "solve OOM decoder initialization failed"))
    {
        return false;
    }
    for (uint32_t id = 0; id + 1u < K; ++id) {
        uint8_t* const packet =
            source.data() + (size_t)id * solve_block_bytes;
        for (uint32_t b = 0; b < solve_block_bytes; ++b) {
            packet[b] = (uint8_t)(id * 7u + b * 11u + 3u);
        }
        if (!Check(
                solve_oom_decoder.DecodeResult(
                    id, packet, solve_block_bytes) == Wirehair_NeedMore,
                "solve OOM decoder ended early"))
        {
            return false;
        }
    }
    uint8_t* const final_packet =
        source.data() + (size_t)(K - 1u) * solve_block_bytes;
    for (uint32_t b = 0; b < solve_block_bytes; ++b) {
        final_packet[b] = (uint8_t)((K - 1u) * 7u + b * 11u + 3u);
    }
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
    const WirehairResult solve_oom = solve_oom_decoder.DecodeResult(
        K - 1u, final_packet, solve_block_bytes);
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    const bool retryable_state =
        solve_oom == Wirehair_OOM &&
        solve_oom_decoder.ReceivedCount() == K &&
        !solve_oom_decoder.IsDecoded();
    const WirehairResult retry_result = solve_oom_decoder.DecodeResult(
        K - 1u, final_packet, solve_block_bytes);
    std::vector<uint8_t> recovered(source.size(), 0u);
    if (!Check(
            retryable_state,
            "solve-path OOM did not preserve retryable packet state") ||
        !Check(
            retry_result == Wirehair_Success,
            "duplicate retry after solve-path OOM did not decode") ||
        !Check(
            solve_oom_decoder.RecoverResult(
                recovered.data(), recovered.size()) == Wirehair_Success &&
                recovered == source,
            "duplicate retry after solve-path OOM recovered wrong data"))
    {
        return false;
    }
#endif
    return true;
}

bool CheckFacadeModeTransitions()
{
    const uint32_t K = 16u;
    const uint32_t block_bytes = 11u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message = MakeMessage((size_t)message_bytes);
    std::vector<uint8_t> block(block_bytes, 0u);
    uint32_t written = 0u;

    wirehair_v2::Codec codec;
    if (!Check(
            codec.InitializeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "initial V1 encoder mode failed"))
    {
        return false;
    }
    wirehair_v2::SeedProfile mismatch = codec.Profile();
    ++mismatch.BlockCount;
    if (!Check(
            codec.InitializePrecodeDecoder(
                message_bytes, block_bytes, &mismatch) ==
                Wirehair_InvalidInput,
            "mismatched V2 decoder profile was accepted") ||
        !Check(
            codec.Encode(0u, block.data(), block_bytes, &written) ==
                Wirehair_Success,
            "failed V2 decoder init did not preserve V1 encoder mode"))
    {
        return false;
    }

    if (!Check(
            codec.InitializePrecodeDecoder(message_bytes, block_bytes) ==
                Wirehair_Success,
            "V1 encoder to V2 decoder transition failed") ||
        !Check(
            codec.Encode(0u, block.data(), block_bytes, &written) ==
                Wirehair_InvalidInput,
            "V2 decoder mode allowed Encode") ||
        !Check(
            codec.Recover(message.data(), message_bytes) == Wirehair_NeedMore,
            "fresh V2 decoder did not request packets"))
    {
        return false;
    }

    const wirehair_v2::SeedProfile decoder_profile = codec.Profile();
    if (!Check(
            codec.InitializePrecodeDecoder(0u, block_bytes) ==
                Wirehair_InvalidInput,
            "invalid V2 decoder reinitialize was accepted") ||
        !Check(
            codec.Profile().BlockCount == decoder_profile.BlockCount &&
                codec.Encode(0u, block.data(), block_bytes, &written) ==
                    Wirehair_InvalidInput,
            "failed V2 decoder reinitialize changed active mode"))
    {
        return false;
    }

    if (!Check(
            codec.InitializePrecodeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "V2 decoder to V2 encoder transition failed") ||
        !Check(
            codec.Encode(K, block.data(), block_bytes, &written) ==
                Wirehair_Success,
            "V2 encoder unusable after decoder transition") ||
        !Check(
            codec.Decode(K, block.data(), written) == Wirehair_InvalidInput,
            "V2 encoder mode allowed Decode"))
    {
        return false;
    }

    if (!Check(
            codec.InitializeDecoder(message_bytes, block_bytes) ==
                Wirehair_Success,
            "V2 encoder to V1 decoder transition failed") ||
        !Check(
            codec.InitializeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "V1 decoder to V1 encoder transition failed") ||
        !Check(
            codec.Encode(0u, block.data(), block_bytes, &written) ==
                Wirehair_Success,
            "V1 encoder unusable after mode transitions"))
    {
        return false;
    }
    return true;
}

bool ParseResumeBenchUint(
    const char* text,
    uint32_t minimum,
    uint32_t maximum,
    uint32_t& value)
{
    if (text == nullptr || *text == '\0') {
        return false;
    }
    for (const char* digit = text; *digit != '\0'; ++digit)
    {
        if (*digit < '0' || *digit > '9') {
            return false;
        }
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        parsed < minimum || parsed > maximum)
    {
        return false;
    }
    value = (uint32_t)parsed;
    return true;
}

bool RunCertifiedResumeBenchmark(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t repetitions,
    uint32_t collision_pair_count)
{
#if !defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    (void)K;
    (void)block_bytes;
    (void)repetitions;
    (void)collision_pair_count;
    std::fprintf(stderr,
        "certified resume benchmark requires test hooks\n");
    return false;
#else
    // Distinct duplicate-row pairs retain valid, distinct packet ids and
    // authentic payloads while forcing an exactly reproducible deficient
    // prefix through the public decoder.  One pair models the common
    // one-rank tail; thirteen drive a deliberately extreme q>H tail beyond
    // the certified H=12 completion budget.
    const size_t collision_packet_count =
        (size_t)collision_pair_count * 2u;
    if (K < 32u || collision_pair_count == 0u ||
        collision_pair_count > 64u ||
        collision_packet_count >= K ||
        block_bytes == 0u || repetitions == 0u)
    {
        std::fprintf(stderr,
            "resume benchmark requires K>=32, block_bytes>0, reps>0, "
            "and at most 64 collision pairs\n");
        return false;
    }
    const uint64_t message_bytes_wide = (uint64_t)K * block_bytes;
    static const uint64_t kMaxResumeBenchMessageBytes =
        UINT64_C(64) * 1024u * 1024u;
    if (message_bytes_wide > std::numeric_limits<size_t>::max() ||
        message_bytes_wide > kMaxResumeBenchMessageBytes)
    {
        std::fprintf(stderr,
            "resume benchmark message exceeds the 64 MiB fixture cap\n");
        return false;
    }
    const size_t message_bytes = (size_t)message_bytes_wide;
    const std::vector<uint8_t> message = MakeMessage(message_bytes);
    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.Completion = wirehair_v2::CompletionField::GF256;
    wirehair_v2::MessagePrecodeEncoder encoder;
    if (encoder.InitializeResult(
            message.data(), message_bytes_wide, block_bytes,
            nullptr, &options) != Wirehair_Success)
    {
        std::fprintf(stderr, "resume benchmark encoder init failed\n");
        return false;
    }

    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = encoder.Profile().V2PacketPeelSeed;
    config.MixCount = encoder.Profile().V2RecoveryMixCount;
    std::vector<uint32_t> collision_ids;
    if (!FindRepairRowCollisions(
            encoder.BlockEncoder().System(),
            config,
            collision_ids,
            (size_t)collision_pair_count))
    {
        std::fprintf(stderr,
            "resume benchmark could not find %zu collision pairs\n",
            (size_t)collision_pair_count);
        return false;
    }

    std::vector<PacketFixture> initial_packets;
    initial_packets.reserve(K);
    PacketFixture packet;
    const uint32_t systematic_prefix =
        K - (uint32_t)collision_packet_count;
    for (uint32_t id = 0u; id < systematic_prefix; ++id)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet)) {
            return false;
        }
        initial_packets.push_back(packet);
    }
    for (uint32_t id : collision_ids)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet)) {
            return false;
        }
        initial_packets.push_back(packet);
    }
    if (initial_packets.size() != K) {
        return false;
    }

    std::vector<PacketFixture> missing_packets;
    missing_packets.reserve(collision_packet_count);
    for (uint32_t id = systematic_prefix; id < K; ++id)
    {
        if (!EncodePacket(encoder, block_bytes, id, packet)) {
            return false;
        }
        missing_packets.push_back(packet);
    }

    struct ModeReset
    {
        ~ModeReset()
        {
            (void)wirehair_v2::
                SetCertifiedPackedResidualModeForTesting(0);
            wirehair_v2::
                SetDecoderIncrementalResumeEnabledForTesting(true);
        }
    } reset;
    (void)reset;

    std::printf(
        "# certified_resume_bench,K=%u,bb=%u,pairs=%zu,reps=%u\n"
        "rep,policy_lane,lane,mode,resume_requested,cold_ns,tail_ns,"
        "tail_packets,retained,tail_path,"
        "resume_bytes,cold_bytes,total_solve_attempts,inactivated,"
        "binary_rank,residual_rank,quotient_deficiency,packed_uses,"
        "materializations,result,fixed_h_uses,fixed_h_bytes,"
        "square_q_columns,heavy_rhs_rows,legacy_replay_rows\n",
        K, block_bytes, (size_t)collision_pair_count, repetitions);

    const int modes[4] = {-1, 1, 1, -1};
    const char* mode_names[4] = {"byte", "packed", "packed", "byte"};
    typedef std::chrono::steady_clock BenchClock;
    for (uint32_t rep = 0u; rep < repetitions; ++rep)
    {
        for (uint32_t policy_order = 0u;
             policy_order < 2u;
             ++policy_order)
        {
            const uint32_t policy_lane = (rep & 1u) == 0u ?
                policy_order : 1u - policy_order;
            const bool resume_requested = policy_lane != 0u;
            bool reference_set = false;
            bool reference_retained = false;
            uint32_t reference_tail_packets = 0u;
            uint32_t reference_inactivated = 0u;
            uint32_t reference_binary_rank = 0u;
            uint32_t reference_residual_rank = 0u;
            size_t reference_resume_bytes = 0u;
            size_t reference_cold_bytes = 0u;
            for (uint32_t lane = 0u; lane < 4u; ++lane)
            {
                if (!wirehair_v2::SetCertifiedPackedResidualModeForTesting(
                        modes[lane]))
                {
                    return false;
                }
                wirehair_v2::
                    SetDecoderIncrementalResumeEnabledForTesting(
                        resume_requested);
                wirehair_v2::MessagePrecodeDecoder decoder;
                if (decoder.InitializeResult(
                        message_bytes_wide,
                        block_bytes,
                        &encoder.Profile(),
                        &options) != Wirehair_Success)
                {
                    std::fprintf(stderr,
                        "resume benchmark decoder init failed "
                        "rep=%u policy=%u lane=%u\n",
                        rep, policy_lane, lane);
                    return false;
                }
                for (size_t i = 0u;
                     i + 1u < initial_packets.size();
                     ++i)
                {
                    const PacketFixture& initial = initial_packets[i];
                    if (decoder.DecodeResult(
                            initial.Id,
                            initial.Data.data(),
                            initial.Bytes) != Wirehair_NeedMore)
                    {
                        std::fprintf(stderr,
                            "resume benchmark early packet failed "
                            "rep=%u policy=%u lane=%u index=%zu\n",
                            rep, policy_lane, lane, i);
                        return false;
                    }
                }

                const PacketFixture& last = initial_packets.back();
                const BenchClock::time_point cold_begin =
                    BenchClock::now();
                WirehairResult result = decoder.DecodeResult(
                    last.Id, last.Data.data(), last.Bytes);
                const BenchClock::time_point cold_end =
                    BenchClock::now();
                if (result != Wirehair_NeedMore)
                {
                    std::fprintf(stderr,
                        "resume benchmark hard prefix was not deficient "
                        "rep=%u policy=%u lane=%u result=%d\n",
                        rep, policy_lane, lane, (int)result);
                    return false;
                }
                const wirehair_v2::PrecodeSolveStats cold_stats =
                    decoder.SolveStats();
                const bool retained =
                    decoder.HasIncrementalResumeStateForTesting();
                const size_t resume_bytes =
                    decoder.IncrementalResumeBytesForTesting();
                const size_t cold_bytes =
                    decoder.ColdReceiveCapacityBytesForTesting();

                uint32_t tail_packets = 0u;
                const BenchClock::time_point tail_begin =
                    BenchClock::now();
                for (const PacketFixture& missing : missing_packets)
                {
                    result = decoder.DecodeResult(
                        missing.Id, missing.Data.data(), missing.Bytes);
                    ++tail_packets;
                    if (result != Wirehair_NeedMore) {
                        break;
                    }
                }
                const BenchClock::time_point tail_end =
                    BenchClock::now();
                const uint64_t cold_ns = (uint64_t)
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        cold_end - cold_begin).count();
                const uint64_t tail_ns = (uint64_t)
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        tail_end - tail_begin).count();
                const uint32_t quotient_deficiency =
                    cold_stats.InactivatedColumns >=
                            cold_stats.BinaryResidualRank ?
                        cold_stats.InactivatedColumns -
                        cold_stats.BinaryResidualRank :
                        UINT32_MAX;

                const bool packed_mode = modes[lane] > 0;
                const bool engagement_ok =
                    cold_stats.CertifiedPackedResidualUses ==
                        (packed_mode ? 1u : 0u);
                const bool policy_state_ok = resume_requested ?
                    (retained ?
                        resume_bytes != 0u && cold_bytes == 0u :
                        resume_bytes == 0u && cold_bytes != 0u) :
                    !retained && resume_bytes == 0u && cold_bytes != 0u;
                const uint64_t materializations =
                    cold_stats.CertifiedPackedResumeMaterializations;
                const bool materialization_ok = !packed_mode ?
                    materializations == 0u :
                    (retained ?
                        materializations == 1u :
                        (resume_requested ?
                            materializations <= 1u :
                            materializations == 0u));
                const bool shape_ok =
                    quotient_deficiency >
                        encoder.BlockEncoder().System().Params.HeavyRows &&
                    tail_packets == collision_pair_count &&
                    decoder.SolveAttemptCount() == 1u + tail_packets;
                const uint32_t heavy_rows =
                    encoder.BlockEncoder().System().Params.HeavyRows;
                const bool fixed_h_expected =
                    packed_mode ||
                    block_bytes >=
                        wirehair_v2::kLegacyByteQuotientMinBlockBytes;
                const bool fixed_h_work_pair =
                    (cold_stats.CertifiedHeavyRhsRowsBuilt == 0u &&
                        cold_stats.CertifiedLegacyHeavyRowsReplayed == 0u) ||
                    (cold_stats.CertifiedHeavyRhsRowsBuilt == heavy_rows &&
                        cold_stats.CertifiedLegacyHeavyRowsReplayed ==
                            heavy_rows);
                const bool fixed_h_checkpoint_work_ok =
                    fixed_h_work_pair &&
                    (resume_requested ||
                        cold_stats.CertifiedHeavyRhsRowsBuilt == 0u) &&
                    (!retained ||
                        cold_stats.CertifiedHeavyRhsRowsBuilt == heavy_rows) &&
                    (materializations == 0u ||
                        cold_stats.CertifiedHeavyRhsRowsBuilt == heavy_rows);
                const bool fixed_h_shape_ok = fixed_h_expected ?
                    cold_stats.CertifiedFixedHQuotientUses == 1u &&
                        cold_stats.
                            CertifiedFixedHCoefficientStorageBytes ==
                            (uint64_t)heavy_rows * heavy_rows &&
                        cold_stats.CertifiedSquareQuotientColumns == 0u &&
                        fixed_h_checkpoint_work_ok :
                    cold_stats.CertifiedFixedHQuotientUses == 0u &&
                        cold_stats.
                            CertifiedFixedHCoefficientStorageBytes == 0u &&
                        cold_stats.CertifiedSquareQuotientColumns == 0u &&
                        cold_stats.CertifiedHeavyRhsRowsBuilt == heavy_rows &&
                        cold_stats.CertifiedLegacyHeavyRowsReplayed == 0u;
                const bool fixed_h_rank_ok =
                    cold_stats.ResidualRank ==
                        cold_stats.BinaryResidualRank + heavy_rows;
                if (!engagement_ok || !policy_state_ok ||
                    !materialization_ok || !shape_ok || !fixed_h_shape_ok ||
                    !fixed_h_rank_ok)
                {
                    std::fprintf(stderr,
                        "resume benchmark lane contract failed "
                        "rep=%u policy=%u lane=%u packed=%llu mat=%llu "
                        "fixed=%u fixed_bytes=%llu square_q=%u "
                        "rhs=%u replay=%u "
                        "retained=%u q=%u tail=%u attempts=%u\n",
                        rep,
                        policy_lane,
                        lane,
                        (unsigned long long)
                            cold_stats.CertifiedPackedResidualUses,
                        (unsigned long long)
                            cold_stats.
                                CertifiedPackedResumeMaterializations,
                        cold_stats.CertifiedFixedHQuotientUses,
                        (unsigned long long)
                            cold_stats.
                                CertifiedFixedHCoefficientStorageBytes,
                        cold_stats.CertifiedSquareQuotientColumns,
                        cold_stats.CertifiedHeavyRhsRowsBuilt,
                        cold_stats.CertifiedLegacyHeavyRowsReplayed,
                        retained ? 1u : 0u,
                        quotient_deficiency,
                        tail_packets,
                        decoder.SolveAttemptCount());
                    return false;
                }
                if (!reference_set)
                {
                    reference_set = true;
                    reference_retained = retained;
                    reference_tail_packets = tail_packets;
                    reference_inactivated =
                        cold_stats.InactivatedColumns;
                    reference_binary_rank =
                        cold_stats.BinaryResidualRank;
                    reference_residual_rank =
                        cold_stats.ResidualRank;
                    reference_resume_bytes = resume_bytes;
                    reference_cold_bytes = cold_bytes;
                }
                else if (retained != reference_retained ||
                    tail_packets != reference_tail_packets ||
                    cold_stats.InactivatedColumns !=
                        reference_inactivated ||
                    cold_stats.BinaryResidualRank !=
                        reference_binary_rank ||
                    cold_stats.ResidualRank !=
                        reference_residual_rank ||
                    resume_bytes != reference_resume_bytes ||
                    cold_bytes != reference_cold_bytes)
                {
                    std::fprintf(stderr,
                        "resume benchmark byte/packed lane mismatch "
                        "rep=%u policy=%u lane=%u\n",
                        rep, policy_lane, lane);
                    return false;
                }

                std::printf(
                    "%u,%u,%u,%s,%u,%llu,%llu,%u,%u,%s,%zu,%zu,%u,"
                    "%u,%u,%u,%u,%llu,%llu,%d,%u,%llu,%u,%u,%u\n",
                    rep,
                    policy_lane,
                    lane,
                    mode_names[lane],
                    resume_requested ? 1u : 0u,
                    (unsigned long long)cold_ns,
                    (unsigned long long)tail_ns,
                    tail_packets,
                    retained ? 1u : 0u,
                    retained ? "resume" : "cold",
                    resume_bytes,
                    cold_bytes,
                    decoder.SolveAttemptCount(),
                    cold_stats.InactivatedColumns,
                    cold_stats.BinaryResidualRank,
                    cold_stats.ResidualRank,
                    quotient_deficiency,
                    (unsigned long long)
                        cold_stats.CertifiedPackedResidualUses,
                    (unsigned long long)
                        cold_stats.CertifiedPackedResumeMaterializations,
                    (int)result,
                    cold_stats.CertifiedFixedHQuotientUses,
                    (unsigned long long)
                        cold_stats.CertifiedFixedHCoefficientStorageBytes,
                    cold_stats.CertifiedSquareQuotientColumns,
                    cold_stats.CertifiedHeavyRhsRowsBuilt,
                    cold_stats.CertifiedLegacyHeavyRowsReplayed);

                if (result != Wirehair_Success)
                {
                    std::fprintf(stderr,
                        "resume benchmark tail did not recover "
                        "rep=%u policy=%u lane=%u result=%d\n",
                        rep, policy_lane, lane, (int)result);
                    return false;
                }
                std::vector<uint8_t> recovered(message_bytes, 0u);
                if (decoder.RecoverResult(
                        recovered.data(), message_bytes_wide) !=
                        Wirehair_Success ||
                    recovered != message)
                {
                    std::fprintf(stderr,
                        "resume benchmark recovery mismatch "
                        "rep=%u policy=%u lane=%u\n",
                        rep, policy_lane, lane);
                    return false;
                }
            }
        }
    }
    return std::fflush(stdout) == 0;
#endif
}

} // namespace

int main(int argc, char** argv)
{
    if ((argc == 5 || argc == 6) &&
        std::string(argv[1]) == "--resume-bench")
    {
        uint32_t K = 0u;
        uint32_t block_bytes = 0u;
        uint32_t repetitions = 0u;
        uint32_t collision_pairs = 1u;
        if (!ParseResumeBenchUint(argv[2], 32u, 64000u, K) ||
            !ParseResumeBenchUint(
                argv[3], 1u, UINT32_C(0x7fffffff), block_bytes) ||
            !ParseResumeBenchUint(argv[4], 1u, 999u, repetitions) ||
            (argc == 6 &&
             !ParseResumeBenchUint(
                 argv[5],
                 1u,
                 std::min<uint32_t>(64u, (K - 1u) / 2u),
                 collision_pairs)))
        {
            std::fprintf(stderr,
                "usage: %s --resume-bench K block_bytes repetitions "
                "[collision_pairs]\n",
                argv[0]);
            return 2;
        }
        try {
            return RunCertifiedResumeBenchmark(
                K, block_bytes, repetitions, collision_pairs) ? 0 : 1;
        }
        catch (const std::bad_alloc&) {
            std::fprintf(stderr,
                "resume benchmark fixture allocation failed\n");
            return 1;
        }
        catch (const std::length_error&) {
            std::fprintf(stderr,
                "resume benchmark fixture size overflowed\n");
            return 1;
        }
    }
    const bool large = argc == 2 && std::string(argv[1]) == "--large";
    const bool large_recovery =
        argc == 2 && std::string(argv[1]) == "--large-recovery";
    if (argc > 2 || (argc == 2 && !large && !large_recovery)) {
        std::fprintf(
            stderr,
            "usage: %s [--large|--large-recovery]\n"
            "       %s --resume-bench K block_bytes repetitions "
            "[collision_pairs]\n",
            argv[0], argv[0]);
        return 2;
    }
    if (!CheckPacketSlotTable() ||
        !CheckInvalidAndOom() ||
        !CheckFacadeModeTransitions() ||
        !CheckIncrementalDecoderParity() ||
        !CheckColdDuplicateSlotLookup() ||
        !CheckDirectRecoveryOutput() ||
        !CheckSystematicRecoverCache())
    {
        return 1;
    }
    if (large) {
        return RunDecodeCase(64000u, 1u, false, false, true, 10u) ? 0 : 1;
    }
    if (large_recovery) {
        return RunDecodeCase(64000u, 1u, true, false, false, 10u) ? 0 : 1;
    }
    if (!RunDecodeCase(64u, 37u, false, true, false, 5u) ||
        !RunDecodeCase(1000u, 3u, false, false, false, 5u) ||
        !RunDecodeCase(32u, 19u, true, true, false, 5u))
    {
        return 1;
    }
    return 0;
}
