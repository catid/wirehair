#include "WirehairV2Plan.h"
#include "WirehairV2Precode.h"
#include "WirehairV2PrecodeDecode.h"
#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"

#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cstring>
#include <vector>

namespace {

bool RunCase(uint32_t K, uint32_t block_bytes, uint32_t loss_stride)
{
    const wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    wirehair_v2::PrecodeParams params = wirehair_v2::MakeCertifiedParams(
        K,
        wirehair_v2::MatrixSeedFromProfile(
            profile, 0u, wirehair_v2::kMessagePrecodeSeedSalt));
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        std::fprintf(stderr, "solve: precode build failed K=%u\n", K);
        return false;
    }

    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
        profile, wirehair_v2::kMessageRecoveryRowSeedSalt);
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;

    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 131u + K * 17u + 29u);
    }
    std::vector<wirehair_v2::SolvePacket> packets;
    for (uint32_t id = 0; id < K; ++id)
    {
        wirehair_v2::SolvePacket packet;
        packet.BlockId = id;
        packet.Data = message.data() + (size_t)id * block_bytes;
        packets.push_back(packet);
    }

    std::vector<uint8_t> intermediate;
    wirehair_v2::PrecodeSolveStats stats;
    const WirehairResult encoded = wirehair_v2::SolvePrecodeSystem(
        system, config, packets, block_bytes, intermediate, &stats);
    if (encoded != Wirehair_Success ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, packets, intermediate.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: systematic solve failed K=%u result=%d R=%u rank=%u\n",
            K, (int)encoded, stats.InactivatedColumns, stats.ResidualRank);
        return false;
    }

    std::vector<uint8_t> block(block_bytes);
    for (uint32_t id = 0; id < K; ++id)
    {
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, intermediate.data(), block_bytes,
                id, block.data()) ||
            std::memcmp(
                block.data(),
                message.data() + (size_t)id * block_bytes,
                block_bytes) != 0)
        {
            std::fprintf(stderr,
                "solve: systematic row mismatch K=%u id=%u\n", K, id);
            return false;
        }
    }

    std::vector<std::vector<uint8_t> > delivered_data;
    std::vector<wirehair_v2::SolvePacket> delivered;
    delivered_data.reserve(K + 32u);
    delivered.reserve(K + 32u);
    for (uint32_t id = K; id-- > 0u;)
    {
        if (id % loss_stride == 0u) {
            continue;
        }
        delivered_data.push_back(std::vector<uint8_t>(block_bytes));
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, intermediate.data(), block_bytes,
                id, delivered_data.back().data()))
        {
            return false;
        }
        wirehair_v2::SolvePacket packet;
        packet.BlockId = id;
        packet.Data = delivered_data.back().data();
        delivered.push_back(packet);
    }
    WirehairResult decoded = Wirehair_NeedMore;
    std::vector<uint8_t> recovered_intermediate;
    uint32_t repair_id = K;
    while (repair_id < K + 32u && decoded != Wirehair_Success)
    {
        delivered_data.push_back(std::vector<uint8_t>(block_bytes));
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, intermediate.data(), block_bytes,
                repair_id, delivered_data.back().data()))
        {
            return false;
        }
        wirehair_v2::SolvePacket packet;
        packet.BlockId = repair_id++;
        packet.Data = delivered_data.back().data();
        delivered.push_back(packet);
        if (delivered.size() >= K) {
            decoded = wirehair_v2::SolvePrecodeSystem(
                system, config, delivered, block_bytes,
                recovered_intermediate, &stats);
        }
    }
    if (decoded != Wirehair_Success ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, delivered,
            recovered_intermediate.data(), block_bytes))
    {
        std::fprintf(stderr,
            "solve: lossy solve failed K=%u result=%d delivered=%zu\n",
            K, (int)decoded, delivered.size());
        return false;
    }
    for (uint32_t id = 0; id < K; ++id)
    {
        if (!wirehair_v2::EvaluatePacketBlock(
                system, config, recovered_intermediate.data(), block_bytes,
                id, block.data()) ||
            std::memcmp(
                block.data(),
                message.data() + (size_t)id * block_bytes,
                block_bytes) != 0)
        {
            std::fprintf(stderr,
                "solve: recovered message mismatch K=%u id=%u\n", K, id);
            return false;
        }
    }

    recovered_intermediate[0] ^= 1u;
    if (wirehair_v2::VerifyPrecodeSolution(
            system, config, delivered,
            recovered_intermediate.data(), block_bytes))
    {
        std::fprintf(stderr, "solve: corrupted solution verified K=%u\n", K);
        return false;
    }

    std::printf(
        "global solve K=%u bb=%u delivered=%zu inact=%u rank=%u: PASS\n",
        K, block_bytes, delivered.size(),
        stats.InactivatedColumns, stats.ResidualRank);
    return true;
}

bool CheckMixDomainValidation()
{
    wirehair_v2::PrecodeParams params;
    params.BlockCount = 2u;
    params.Staircase = 2u;
    params.DenseRows = 0u;
    params.HeavyRows = 0u;
    params.SourceHits = 1u;
    params.Seed = 7u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        return false;
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = 11u;
    config.MixCount = 3u;
    if (!wirehair_v2::GeneratePacketMatrixRow(2u, 2u, 0u, config).empty()) {
        std::fprintf(stderr, "solve: oversized mix generated duplicates\n");
        return false;
    }

    const uint8_t intermediate[4] = {1u, 2u, 3u, 4u};
    uint8_t output = 0xa5u;
    if (wirehair_v2::EvaluatePacketBlock(
            system, config, intermediate, 1u, 0u, &output) ||
        output != 0xa5u)
    {
        std::fprintf(stderr, "solve: invalid mix modified packet output\n");
        return false;
    }
    const uint8_t zero = 0u;
    std::vector<wirehair_v2::SolvePacket> packets(2u);
    packets[0].BlockId = 0u;
    packets[0].Data = &zero;
    packets[1].BlockId = 1u;
    packets[1].Data = &zero;
    std::vector<uint8_t> solved(3u, 0xccu);
    const std::vector<uint8_t> before = solved;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, packets, 1u, solved) != Wirehair_InvalidInput ||
        solved != before)
    {
        std::fprintf(stderr, "solve: invalid mix solve was not no-write\n");
        return false;
    }
    return true;
}

bool CheckLargePacketEvaluationWork()
{
    const uint32_t K = 64000u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x987654321));
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        return false;
    }
    const uint32_t L = K + params.Staircase +
        params.DenseRows + params.HeavyRows;
    std::vector<uint8_t> intermediate(L);
    for (uint32_t i = 0; i < L; ++i) {
        intermediate[i] = (uint8_t)(i * 17u + 3u);
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = 0x5eedu;
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;

    uint64_t total_work = 0u;
    uint64_t digest = 0u;
    uint8_t output = 0u;
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (uint32_t i = 0; i < 4096u; ++i)
    {
        uint64_t work = 0u;
        const uint32_t id = UINT32_MAX - i * 7919u;
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, intermediate.data(), 1u, id, &output, &work) ||
            work < 4u || work > 67u)
        {
            std::fprintf(stderr,
                "solve: K=64000 packet work invalid id=%u work=%llu\n",
                id, (unsigned long long)work);
            return false;
        }
        total_work += work;
        digest = digest * UINT64_C(0x9e3779b97f4a7c15) + output + work;
    }
    const double milliseconds =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin).count();
    if (total_work > UINT64_C(4096) * 67u || digest == 0u) {
        return false;
    }
    std::printf(
        "K=64000 packet evaluation: rows=4096 work=%llu time=%.3fms "
        "digest=%llu: PASS\n",
        (unsigned long long)total_work,
        milliseconds,
        (unsigned long long)digest);
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 3)
    {
        const uint32_t K = (uint32_t)std::strtoul(argv[1], nullptr, 10);
        const uint32_t block_bytes =
            (uint32_t)std::strtoul(argv[2], nullptr, 10);
        if (K < 2u || K > 64000u || block_bytes == 0u) {
            return 2;
        }
        std::vector<uint8_t> message((size_t)K * block_bytes, 0x5au);
        const wirehair_v2::SeedProfile diagnostic_profile =
            wirehair_v2::SelectSeedProfile(K, block_bytes);
        wirehair_v2::PrecodeParams diagnostic_params =
            wirehair_v2::MakeCertifiedParams(
                K,
                wirehair_v2::MatrixSeedFromProfile(
                    diagnostic_profile,
                    0u,
                    wirehair_v2::kMessagePrecodeSeedSalt));
        wirehair_v2::PrecodeSystem diagnostic_system;
        if (!wirehair_v2::BuildPrecodeSystem(
                diagnostic_params, diagnostic_system))
        {
            return 1;
        }
        wirehair_v2::PacketRowConfig diagnostic_config;
        diagnostic_config.PeelSeed =
            wirehair_v2::PacketPeelSeedFromProfile(
                diagnostic_profile,
                wirehair_v2::kMessageRecoveryRowSeedSalt);
        diagnostic_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
        std::printf(
            "seeds profile_peel=%u profile_dense=%u dense_count=%u "
            "precode=0x%llx packet=0x%x\n",
            diagnostic_profile.PeelSeed,
            diagnostic_profile.DenseSeed,
            diagnostic_profile.DenseCount,
            (unsigned long long)diagnostic_params.Seed,
            diagnostic_config.PeelSeed);
        std::vector<wirehair_v2::SolvePacket> diagnostic_packets(K);
        for (uint32_t id = 0; id < K; ++id) {
            diagnostic_packets[id].BlockId = id;
            diagnostic_packets[id].Data =
                message.data() + (size_t)id * block_bytes;
        }
        uint32_t direct_binary_rank = 0u;
        if (K + diagnostic_params.Staircase + diagnostic_params.DenseRows +
                diagnostic_params.HeavyRows <= 64u)
        {
            std::vector<uint64_t> masks;
            const auto add_mask = [&](const std::vector<uint32_t>& row) {
                uint64_t mask = 0u;
                for (uint32_t column : row) {
                    mask ^= UINT64_C(1) << column;
                }
                masks.push_back(mask);
            };
            for (const std::vector<uint32_t>& row :
                    diagnostic_system.StaircaseRows) {
                add_mask(row);
            }
            for (const std::vector<uint32_t>& row :
                    diagnostic_system.DenseRowColumns) {
                add_mask(row);
            }
            for (uint32_t id = 0; id < K; ++id) {
                add_mask(wirehair_v2::GeneratePacketMatrixRow(
                    K,
                    diagnostic_params.Staircase +
                        diagnostic_params.DenseRows +
                        diagnostic_params.HeavyRows,
                    id,
                    diagnostic_config));
            }
            for (uint32_t column = 0; column < 64u; ++column)
            {
                uint32_t pivot = direct_binary_rank;
                while (pivot < masks.size() &&
                       (masks[pivot] & (UINT64_C(1) << column)) == 0u) {
                    ++pivot;
                }
                if (pivot == masks.size()) {
                    continue;
                }
                std::swap(masks[pivot], masks[direct_binary_rank]);
                for (uint32_t r = 0; r < masks.size(); ++r) {
                    if (r != direct_binary_rank &&
                        (masks[r] & (UINT64_C(1) << column)) != 0u) {
                        masks[r] ^= masks[direct_binary_rank];
                    }
                }
                ++direct_binary_rank;
            }
        }
        std::vector<uint8_t> diagnostic_intermediate;
        wirehair_v2::PrecodeSolveStats diagnostic_stats;
        const WirehairResult diagnostic_result =
            wirehair_v2::SolvePrecodeSystem(
                diagnostic_system,
                diagnostic_config,
                diagnostic_packets,
                block_bytes,
                diagnostic_intermediate,
                &diagnostic_stats);
        std::printf(
            "base result=%d peeled=%u inact=%u residual_rows=%u "
            "binary_rank=%u direct_binary=%u rank=%u\n",
            (int)diagnostic_result,
            diagnostic_stats.PeeledColumns,
            diagnostic_stats.InactivatedColumns,
            diagnostic_stats.ResidualRows,
            diagnostic_stats.BinaryResidualRank,
            direct_binary_rank,
            diagnostic_stats.ResidualRank);
        wirehair_v2::MessagePrecodeEncoder encoder;
        const WirehairResult result = encoder.InitializeResult(
            message.data(), message.size(), block_bytes);
        const wirehair_v2::PrecodeSolveStats& stats = encoder.SolveStats();
        std::printf(
            "profile K=%u bb=%u result=%d attempt=%u inact=%u rank=%u "
            "build=%.3fms peel=%.3fms project=%.3fms residual=%.3fms "
            "backsub=%.3fms\n",
            K, block_bytes, (int)result,
            stats.PacketSeedAttempt,
            stats.InactivatedColumns, stats.ResidualRank,
            stats.BuildNanoseconds / 1000000.0,
            stats.PeelNanoseconds / 1000000.0,
            stats.ProjectNanoseconds / 1000000.0,
            stats.ResidualNanoseconds / 1000000.0,
            stats.BackSubNanoseconds / 1000000.0);
        if (result == Wirehair_Success)
        {
            wirehair_v2::MessagePrecodeDecoder decoder;
            if (decoder.InitializeResult(
                    message.size(), block_bytes, &encoder.Profile()) !=
                    Wirehair_Success)
            {
                return 1;
            }
            std::vector<uint8_t> block(block_bytes);
            WirehairResult decode_result = Wirehair_NeedMore;
            uint32_t delivered = 0u;
            for (uint32_t id = 0u;
                 decode_result == Wirehair_NeedMore && id < K * 2u;
                 ++id)
            {
                if (id % 10u == 0u) {
                    continue;
                }
                uint32_t bytes = 0u;
                if (encoder.EncodeResult(
                        id, block.data(), block_bytes, &bytes) !=
                        Wirehair_Success)
                {
                    return 1;
                }
                decode_result = decoder.DecodeResult(
                    id, block.data(), bytes);
                ++delivered;
            }
            const wirehair_v2::PrecodeSolveStats& decode_stats =
                decoder.SolveStats();
            std::printf(
                "decode result=%d delivered=%u inact=%u rank=%u "
                "build=%.3fms peel=%.3fms project=%.3fms residual=%.3fms "
                "backsub=%.3fms\n",
                (int)decode_result, delivered,
                decode_stats.InactivatedColumns,
                decode_stats.ResidualRank,
                decode_stats.BuildNanoseconds / 1000000.0,
                decode_stats.PeelNanoseconds / 1000000.0,
                decode_stats.ProjectNanoseconds / 1000000.0,
                decode_stats.ResidualNanoseconds / 1000000.0,
                decode_stats.BackSubNanoseconds / 1000000.0);
            if (decode_result != Wirehair_Success) {
                return 1;
            }
        }
        return result == Wirehair_Success ? 0 : 1;
    }
    if (argc != 1) {
        return 2;
    }
    bool ok = true;
    ok = CheckMixDomainValidation() && ok;
    ok = CheckLargePacketEvaluationWork() && ok;
    ok = RunCase(64u, 17u, 7u) && ok;
    ok = RunCase(320u, 37u, 11u) && ok;
    return ok ? 0 : 1;
}
