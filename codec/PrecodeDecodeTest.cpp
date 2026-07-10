#include "WirehairV2Codec.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <map>
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
        a.BlockXors == b.BlockXors &&
        a.BlockMulAdds == b.BlockMulAdds &&
        a.BuildNanoseconds == b.BuildNanoseconds &&
        a.PeelNanoseconds == b.PeelNanoseconds &&
        a.ProjectNanoseconds == b.ProjectNanoseconds &&
        a.ResidualNanoseconds == b.ResidualNanoseconds &&
        a.BackSubNanoseconds == b.BackSubNanoseconds &&
        a.PacketSeedAttempt == b.PacketSeedAttempt;
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
    std::vector<uint32_t>& collision_ids)
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
    std::vector<uint32_t> first_collision_row;
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
            if (collision_ids.empty())
            {
                collision_ids.push_back(found.first->second.First);
                collision_ids.push_back(id);
                found.first->second.Included = true;
                first_collision_row = row;
                continue;
            }
            if (row != first_collision_row)
            {
                collision_ids.push_back(found.first->second.First);
                collision_ids.push_back(id);
                found.first->second.Included = true;
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
                fallback_pair.Warm.ColdReceiveCapacityBytesForTesting() != 0u,
            "oversized checkpoint did not retain the cold fallback"))
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

} // namespace

int main(int argc, char** argv)
{
    const bool large = argc == 2 && std::string(argv[1]) == "--large";
    const bool large_recovery =
        argc == 2 && std::string(argv[1]) == "--large-recovery";
    if (argc > 2 || (argc == 2 && !large && !large_recovery)) {
        std::fprintf(
            stderr, "usage: %s [--large|--large-recovery]\n", argv[0]);
        return 2;
    }
    if (!CheckInvalidAndOom() ||
        !CheckFacadeModeTransitions() ||
        !CheckIncrementalDecoderParity())
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
