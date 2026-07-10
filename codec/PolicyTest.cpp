#include "WirehairV2Policy.h"
#include "WirehairV2Codec.h"
#include "WirehairV2Peel.h"
#include "WirehairV2Plan.h"
#include "WirehairV2Precode.h"
#include "WirehairV2Seeds.h"
#include "../gf256.h"
#include "../WirehairTools.h"
#include <wirehair/wirehair.h>

#include <chrono>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {

int Failures = 0;

void Check(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << std::endl;
        ++Failures;
    }
}

void CheckPolicy(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::BlockByteClass byte_class,
    wirehair_v2::BlockCountBand count_band,
    wirehair_v2::PeelStructure structure,
    wirehair_v2::DegreeFamily family,
    uint16_t min_degree,
    uint16_t max_degree)
{
    const wirehair_v2::PeelPolicy policy =
        wirehair_v2::SelectPeelPolicy(block_count, block_bytes);
    const wirehair_v2::PeelingCodec codec =
        wirehair_v2::SelectPeelingCodec(block_count, block_bytes);

    Check(policy.Solver == wirehair_v2::PeelSolver::KsBmaxTop16,
        "selected solver should be ks_bmax_top16");
    Check(policy.ByteClass == byte_class, "unexpected block byte class");
    Check(policy.CountBand == count_band, "unexpected block count band");
    Check(policy.Structure == structure, "unexpected peel structure");
    Check(policy.Codec.Structure == structure,
        "policy codec structure should match policy structure");
    Check(codec.Solver == wirehair_v2::PeelSolver::KsBmaxTop16,
        "selected codec solver should be ks_bmax_top16");
    Check(codec.Structure == structure, "unexpected selected codec structure");
    Check(codec.Family == family, "unexpected selected codec family");
    Check(codec.MinDegree == min_degree, "unexpected selected min degree");
    Check(codec.MaxDegree == max_degree, "unexpected selected max degree");
    Check(codec.SolverCandidateLimit == 16u,
        "ks_bmax_top16 should use 16 candidates");
    Check(codec.FullyRandomRows, "selected codec rows should be fully random");
}

void CheckRowHasNoDuplicates(const std::vector<uint16_t>& row)
{
    for (size_t i = 0; i < row.size(); ++i)
    {
        for (size_t j = i + 1u; j < row.size(); ++j) {
            Check(row[i] != row[j], "generated row should not duplicate columns");
        }
    }
}

void CheckRowHasNoDuplicates(const std::vector<uint32_t>& row)
{
    for (size_t i = 0; i < row.size(); ++i)
    {
        for (size_t j = i + 1u; j < row.size(); ++j) {
            Check(row[i] != row[j],
                "generated recovery row should not duplicate columns");
        }
    }
}

void CheckV2RoundTrip()
{
    const uint32_t N = 40u;
    const uint32_t block_bytes = 96u;
    const uint64_t message_bytes = (uint64_t)(N - 1u) * block_bytes + 37u;

    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> decoded((size_t)message_bytes, 0xa5);
    std::vector<uint8_t> block(block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 17u + 3u);
    }

    wirehair_v2::Codec encoder;
    wirehair_v2::Codec decoder;
    Check(encoder.InitializeEncoder(
            &message[0], message_bytes, block_bytes) == Wirehair_Success,
        "v2 encoder initialization should succeed");
    Check(decoder.InitializeDecoder(
            message_bytes, block_bytes) == Wirehair_Success,
        "v2 decoder initialization should succeed");
    Check(encoder.Profile().PeelSeed == decoder.Profile().PeelSeed,
        "encoder and decoder should select the same peel seed");
    Check(encoder.Profile().DenseSeed == decoder.Profile().DenseSeed,
        "encoder and decoder should select the same dense seed");

    bool decoded_ok = false;
    uint32_t delivered = 0;
    for (uint32_t block_id = 0; block_id < N + 16u && !decoded_ok; ++block_id)
    {
        if ((block_id % 7u) == 3u) {
            continue;
        }
        uint32_t write_bytes = 0;
        Check(encoder.Encode(block_id, &block[0], block_bytes, &write_bytes) ==
                Wirehair_Success,
            "v2 encode should succeed");
        ++delivered;
        const WirehairResult result =
            decoder.Decode(block_id, &block[0], write_bytes);
        if (result == Wirehair_Success) {
            decoded_ok = true;
            break;
        }
        Check(result == Wirehair_NeedMore,
            "v2 decode should need more data before success");
    }
    Check(decoded_ok, "v2 decoder should complete under a simple loss pattern");
    Check(delivered >= N, "v2 delivered count should cover at least N blocks");
    Check(decoder.Recover(&decoded[0], message_bytes) == Wirehair_Success,
        "v2 recover should succeed");
    Check(std::memcmp(&decoded[0], &message[0], (size_t)message_bytes) == 0,
        "v2 recovered message should match input");
}

void CheckV2PrecodeEncoderFacade()
{
    using namespace wirehair_v2;

    const uint32_t N = 1000u;
    const uint32_t block_bytes = 37u;
    const uint32_t final_bytes = 13u;
    const uint64_t message_bytes =
        (uint64_t)(N - 1u) * block_bytes + final_bytes;

    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 29u + 11u);
    }

    MessagePrecodeEncoderOptions options;
    options.DenseIdentityCorner = true;
    options.RecoveryMixCount = 5u;

    Codec encoder;
    Check(encoder.InitializePrecodeEncoder(
            &message[0], message_bytes, block_bytes, 0, &options) ==
            Wirehair_Success,
        "v2 precode facade initialization should succeed");
    Check(encoder.Profile().BlockCount == N &&
            encoder.Profile().BlockBytes == block_bytes,
        "v2 precode facade should publish selected profile");

    MessagePrecodeEncoder reference;
    Check(reference.Initialize(
            &message[0], message_bytes, block_bytes,
            &encoder.Profile(), &options),
        "v2 precode reference initialization should succeed");

    std::vector<uint8_t> got(block_bytes, 0xac);
    std::vector<uint8_t> want(block_bytes, 0xbd);
    uint32_t written = 0;
    Check(encoder.Encode(0u, &got[0], block_bytes, &written) ==
            Wirehair_Success,
        "v2 precode facade should encode a source block");
    Check(written == block_bytes &&
            std::memcmp(&got[0], &message[0], block_bytes) == 0,
        "v2 precode facade source block should match input");

    for (size_t i = 0; i < got.size(); ++i) {
        got[i] = 0xac;
    }
    Check(encoder.Encode(N - 1u, &got[0], final_bytes, &written) ==
            Wirehair_Success,
        "v2 precode facade should encode a partial final source block");
    Check(written == final_bytes &&
            std::memcmp(
                &got[0],
                &message[(size_t)(N - 1u) * block_bytes],
                final_bytes) == 0,
        "v2 precode facade final source bytes should match input");
    bool tail_untouched = true;
    for (uint32_t i = final_bytes; i < block_bytes; ++i) {
        tail_untouched = tail_untouched && got[i] == 0xac;
    }
    Check(tail_untouched,
        "v2 precode facade should not write beyond final source bytes");

    const uint32_t recovery_id = N + 5u;
    uint32_t want_written = 0;
    Check(encoder.Encode(recovery_id, &got[0], block_bytes, &written) ==
            Wirehair_Success,
        "v2 precode facade should encode a recovery block");
    Check(reference.Encode(
            recovery_id, &want[0], block_bytes, &want_written),
        "v2 precode reference should encode matching recovery block");
    Check(written == block_bytes &&
            want_written == block_bytes &&
            std::memcmp(&got[0], &want[0], block_bytes) == 0,
        "v2 precode facade recovery block should match reference");
    Check(encoder.Decode(recovery_id, &got[0], block_bytes) ==
            Wirehair_InvalidInput,
        "v2 precode facade should reject decode in encoder-only mode");
    std::vector<uint8_t> recovered((size_t)message_bytes);
    Check(encoder.Recover(&recovered[0], message_bytes) ==
            Wirehair_InvalidInput,
        "v2 precode facade should reject recover in encoder-only mode");

    const SeedProfile good_profile = encoder.Profile();
    SeedProfile mismatch = good_profile;
    ++mismatch.BlockCount;
    Check(encoder.InitializePrecodeEncoder(
            &message[0], message_bytes, block_bytes, &mismatch, &options) ==
            Wirehair_InvalidInput,
        "v2 precode facade should reject mismatched seed overrides");
    Check(encoder.Profile().BlockCount == good_profile.BlockCount &&
            encoder.Profile().BlockBytes == good_profile.BlockBytes,
        "v2 precode failed init should preserve profile");
    Check(encoder.Encode(recovery_id, &got[0], block_bytes, &written) ==
            Wirehair_Success &&
            std::memcmp(&got[0], &want[0], block_bytes) == 0,
        "v2 precode failed init should preserve last good encoder");

    Check(encoder.InitializeEncoder(
            &message[0], message_bytes, block_bytes) == Wirehair_Success,
        "v2 facade should switch back to the existing encoder");
    for (size_t i = 0; i < got.size(); ++i) {
        got[i] = 0xac;
    }
    Check(encoder.Encode(N - 1u, &got[0], final_bytes, &written) ==
            Wirehair_Success &&
            written == final_bytes,
        "v2 existing encoder should remain usable after precode mode");
}

void CheckV2InvalidInitialization()
{
    const uint32_t N = 40u;
    const uint32_t block_bytes = 96u;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes, 0x42);

    wirehair_v2::Codec decoder;
    Check(decoder.InitializeDecoder(message_bytes, block_bytes) ==
            Wirehair_Success,
        "v2 decoder baseline initialization should succeed");
    const wirehair_v2::SeedProfile profile = decoder.Profile();

    Check(decoder.InitializeDecoder(UINT64_MAX, 16u) ==
            Wirehair_BadInput_LargeN,
        "v2 decoder should reject huge message sizes before narrowing");
    Check(decoder.Profile().BlockCount == profile.BlockCount &&
            decoder.Profile().BlockBytes == profile.BlockBytes &&
            decoder.Profile().PeelSeed == profile.PeelSeed &&
            decoder.Profile().DenseSeed == profile.DenseSeed &&
            decoder.Profile().DenseCount == profile.DenseCount,
        "v2 failed input validation should preserve the last good profile");

    wirehair_v2::SeedProfile invalid_profile = profile;
    invalid_profile.DenseCount = 0u;
    Check(decoder.InitializeDecoder(message_bytes, block_bytes,
            &invalid_profile) == Wirehair_InvalidInput,
        "v2 decoder should reject invalid dense-count overrides");
    Check(decoder.Profile().BlockCount == profile.BlockCount &&
            decoder.Profile().BlockBytes == profile.BlockBytes &&
            decoder.Profile().PeelSeed == profile.PeelSeed &&
            decoder.Profile().DenseSeed == profile.DenseSeed &&
            decoder.Profile().DenseCount == profile.DenseCount,
        "v2 failed core initialization should preserve the last good profile");

    wirehair_v2::Codec encoder;
    Check(encoder.InitializeEncoder(&message[0], message_bytes, block_bytes) ==
            Wirehair_Success,
        "v2 encoder baseline initialization should succeed");
    invalid_profile = encoder.Profile();
    invalid_profile.DenseCount = 0u;
    Check(encoder.InitializeEncoder(&message[0], message_bytes, block_bytes,
            &invalid_profile) == Wirehair_InvalidInput,
        "v2 encoder should reject invalid dense-count overrides");

    wirehair_v2::SeedProfile mismatched_profile = encoder.Profile();
    mismatched_profile.BlockCount = N + 1u;
    Check(encoder.InitializeEncoder(&message[0], message_bytes, block_bytes,
            &mismatched_profile) == Wirehair_InvalidInput,
        "v2 encoder should reject seed overrides for a different block count");
    mismatched_profile = encoder.Profile();
    mismatched_profile.BlockBytes = block_bytes + 1u;
    Check(encoder.InitializeEncoder(&message[0], message_bytes, block_bytes,
            &mismatched_profile) == Wirehair_InvalidInput,
        "v2 encoder should reject seed overrides for different block bytes");
    Check(decoder.InitializeDecoder(message_bytes, block_bytes,
            &mismatched_profile) == Wirehair_InvalidInput,
        "v2 decoder should reject seed overrides for different block bytes");
}

void CheckV2EncodeFailureClearsByteCount()
{
    const uint32_t N = 4u;
    const uint32_t block_bytes = 32u;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes, 0x42);
    std::vector<uint8_t> block(block_bytes);

    wirehair_v2::Codec encoder;
    Check(encoder.InitializeEncoder(&message[0], message_bytes, block_bytes) ==
            Wirehair_Success,
        "v2 byte-count test encoder initialization should succeed");

    uint32_t written = 0x12345678u;
    Check(encoder.Encode(0u, &block[0], block_bytes - 1u, &written) ==
            Wirehair_InvalidInput,
        "v2 encode should reject too-small output buffers");
    Check(written == 0u,
        "v2 failed encode should clear the byte-count output");
}

void CheckGeneratedRows(
    const wirehair_v2::PeelingCodec& codec,
    uint32_t block_count,
    uint32_t row_count,
    uint64_t seed)
{
    const std::vector<std::vector<uint16_t> > rows =
        wirehair_v2::GeneratePeelMatrixRows(
            codec, block_count, row_count, seed);
    const std::vector<std::vector<uint16_t> > rows_repeat =
        wirehair_v2::GeneratePeelMatrixRows(
            codec, block_count, row_count, seed);

    Check(rows.size() == row_count,
        "generated matrix should contain requested rows");
    Check(rows == rows_repeat, "generated matrix rows should be deterministic");
    for (size_t i = 0; i < rows.size(); ++i)
    {
        const std::vector<uint16_t> single =
            wirehair_v2::GeneratePeelMatrixRow(
                codec, block_count, (uint32_t)i, seed);
        Check(single == rows[i],
            "single peel row should match generated matrix row");
        Check(!rows[i].empty(), "generated row should not be empty");
        Check(rows[i].size() <= codec.MaxDegree,
            "generated row should obey degree cap");
        Check(rows[i].size() <= block_count,
            "generated row should not exceed block count");
        CheckRowHasNoDuplicates(rows[i]);
    }
}

void CheckGeneratedRecoveryRows()
{
    const uint32_t K = 1000u;
    const uint32_t row_count = 32u;
    const uint32_t mix_count = 3u;
    const uint64_t seed = UINT64_C(0x12345678);
    const wirehair_v2::PeelingCodec codec =
        wirehair_v2::MakePeelingCodec(
            wirehair_v2::PeelStructure::LtM1C32,
            wirehair_v2::PeelSolver::KsBmaxTop16);
    const wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, seed);
    const uint32_t precode_count =
        params.Staircase + params.DenseRows + params.HeavyRows;

    const std::vector<std::vector<uint16_t> > source_rows =
        wirehair_v2::GeneratePeelMatrixRows(codec, K, row_count, seed);
    const std::vector<std::vector<uint32_t> > rows =
        wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, precode_count, row_count, mix_count, seed);
    const std::vector<std::vector<uint32_t> > rows_repeat =
        wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, precode_count, row_count, mix_count, seed);

    Check(rows.size() == row_count,
        "recovery row generator should contain requested rows");
    Check(rows == rows_repeat,
        "recovery row generator should be deterministic");
    for (size_t i = 0; i < rows.size(); ++i)
    {
        const std::vector<uint32_t>& row = rows[i];
        const std::vector<uint32_t> single =
            wirehair_v2::GenerateRecoveryMatrixRow(
                codec, K, precode_count, (uint32_t)i, mix_count, seed);
        Check(single == row,
            "single recovery row should match generated matrix row");
        Check(row.size() == source_rows[i].size() + mix_count,
            "recovery row should append precode mix columns");
        for (size_t j = 0; j < source_rows[i].size(); ++j) {
            Check(row[j] == source_rows[i][j],
                "recovery row source prefix should match peel row");
        }
        for (size_t j = source_rows[i].size(); j < row.size(); ++j)
        {
            Check(row[j] >= K && row[j] < K + precode_count,
                "recovery row mix column should be in precode range");
        }
        CheckRowHasNoDuplicates(row);
    }

    const std::vector<std::vector<uint32_t> > no_precode =
        wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, 0u, row_count, mix_count, seed);
    Check(no_precode.size() == row_count,
        "zero-precode recovery row generator should contain requested rows");
    for (size_t i = 0; i < no_precode.size(); ++i)
    {
        Check(no_precode[i].size() == source_rows[i].size(),
            "zero-precode recovery rows should have only source columns");
        for (size_t j = 0; j < source_rows[i].size(); ++j) {
            Check(no_precode[i][j] == source_rows[i][j],
                "mix stream changes must not retune source rows");
        }
    }

    const std::vector<std::vector<uint32_t> > capped_mix =
        wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, 2u, row_count, 5u, seed);
    Check(capped_mix.size() == row_count,
        "capped-mix recovery row generator should contain requested rows");
    for (size_t i = 0; i < capped_mix.size(); ++i)
    {
        Check(capped_mix[i].size() == source_rows[i].size() + 2u,
            "recovery row mix count should cap at precode column count");
        for (size_t j = 0; j < source_rows[i].size(); ++j) {
            Check(capped_mix[i][j] == source_rows[i][j],
                "separate mix RNG must preserve the source prefix");
        }
    }

    Check(wirehair_v2::GenerateRecoveryMatrixRows(
            codec, 0u, precode_count, row_count, mix_count, seed).empty(),
        "recovery row generator should reject zero source count");
    Check(wirehair_v2::GenerateRecoveryMatrixRows(
            codec, UINT16_MAX + 1u, precode_count, row_count,
            mix_count, seed).empty(),
        "recovery row generator should reject source counts above uint16");
    Check(wirehair_v2::GenerateRecoveryMatrixRows(
            codec, UINT16_MAX, 1u, row_count, mix_count, seed).empty(),
        "recovery row generator should reject total columns above uint16");
    Check(wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, precode_count, wirehair_v2::kMaxPeelMatrixRows + 1u,
            mix_count, seed).empty(),
        "recovery row generator should reject unrepresentable row counts");
    const uint32_t shuffled_ids[] = {
        65536u, UINT32_C(0xf1234567), UINT32_MAX, 7u
    };
    for (uint32_t row_id : shuffled_ids)
    {
        wirehair_v2::RecoveryRowGenerationStats stats;
        const std::vector<uint32_t> high_row =
            wirehair_v2::GenerateRecoveryMatrixRow(
                codec, K, precode_count, row_id, mix_count, seed, &stats);
        const std::vector<uint32_t> repeated =
            wirehair_v2::GenerateRecoveryMatrixRow(
                codec, K, precode_count, row_id, mix_count, seed);
        const std::vector<uint16_t> source =
            wirehair_v2::GeneratePeelMatrixRow(codec, K, row_id, seed);
        Check(!high_row.empty() && high_row == repeated,
            "arbitrary uint32 recovery rows should be deterministic");
        Check(stats.ContractVersion ==
                wirehair_v2::kRecoveryRowContractVersion &&
                stats.SeekWork == 2u &&
                stats.SourceRandomDraws > 0u &&
                stats.MixRandomDraws >= mix_count,
            "recovery row work counters should report constant seek work");
        Check(high_row.size() == source.size() + mix_count,
            "high recovery row should retain the peel source prefix");
        for (size_t i = 0; i < source.size(); ++i) {
            Check(high_row[i] == source[i],
                "high recovery source prefix should match peel row");
        }
    }

    struct GoldenRow
    {
        uint32_t Id;
        std::vector<uint32_t> Columns;
    };
    const GoldenRow golden[] = {
        { 0u, {485u, 591u, 873u, 952u, 972u, 159u, 441u,
            1058u, 1054u, 1015u} },
        { 1u, {131u, 517u, 1020u, 1024u, 1010u} },
        { 65536u, {852u, 457u, 1030u, 1049u, 1045u} },
        { UINT32_C(0xf1234567), {418u, 272u, 1041u, 1022u, 1012u} },
        { UINT32_MAX, {970u, 530u, 388u, 112u, 975u, 603u,
            1035u, 1026u, 1064u} },
    };
    for (const GoldenRow& expected : golden)
    {
        Check(wirehair_v2::GenerateRecoveryMatrixRow(
                codec, K, 81u, expected.Id, 3u, seed) == expected.Columns,
            "recovery row contract-v2 golden vector changed");
    }
}

void CheckRecoveryRowScaling()
{
    using Clock = std::chrono::steady_clock;
    const uint32_t K = 1000u;
    const uint32_t precode_count = 81u;
    const uint32_t count = 4096u;
    const uint64_t seed = UINT64_C(0x5ca1ab1e);
    const wirehair_v2::PeelingCodec codec =
        wirehair_v2::MakePeelingCodec(
            wirehair_v2::PeelStructure::LtM1C32,
            wirehair_v2::PeelSolver::KsBmaxTop16);

    const Clock::time_point batch_start = Clock::now();
    const std::vector<std::vector<uint32_t> > batch =
        wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, precode_count, count, 5u, seed);
    const Clock::time_point batch_end = Clock::now();

    uint64_t sequential_work = 0u;
    uint64_t first_quarter_work = 0u;
    const Clock::time_point sequential_start = Clock::now();
    for (uint32_t i = 0; i < count; ++i)
    {
        wirehair_v2::RecoveryRowGenerationStats stats;
        const std::vector<uint32_t> row =
            wirehair_v2::GenerateRecoveryMatrixRow(
                codec, K, precode_count, i, 5u, seed, &stats);
        Check(row == batch[i],
            "4096-row batch/single recovery identity failed");
        const uint64_t work = stats.SeekWork + stats.SourceRandomDraws +
            stats.MixRandomDraws;
        sequential_work += work;
        if (i < count / 4u) {
            first_quarter_work += work;
        }
    }
    const Clock::time_point sequential_end = Clock::now();

    uint64_t shuffled_work = 0u;
    const Clock::time_point shuffled_start = Clock::now();
    for (uint32_t i = 0; i < count; ++i)
    {
        const uint32_t row_id =
            i * UINT32_C(2654435761) + UINT32_C(0x80000000);
        wirehair_v2::RecoveryRowGenerationStats stats;
        const std::vector<uint32_t> row =
            wirehair_v2::GenerateRecoveryMatrixRow(
                codec, K, precode_count, row_id, 5u, seed, &stats);
        Check(!row.empty() && stats.SeekWork == 2u,
            "shuffled uint32 recovery generation failed");
        shuffled_work += stats.SeekWork + stats.SourceRandomDraws +
            stats.MixRandomDraws;
    }
    const Clock::time_point shuffled_end = Clock::now();

    Check(batch.size() == count,
        "4096-row recovery scaling batch failed");
    Check(sequential_work > first_quarter_work * 3u &&
            sequential_work < first_quarter_work * 5u,
        "recovery row deterministic work did not scale near-linearly");
    Check(shuffled_work * 4u > sequential_work * 3u &&
            shuffled_work * 3u < sequential_work * 4u,
        "shuffled high IDs changed recovery generation work materially");

    const double batch_ms = std::chrono::duration<double, std::milli>(
        batch_end - batch_start).count();
    const double sequential_ms = std::chrono::duration<double, std::milli>(
        sequential_end - sequential_start).count();
    const double shuffled_ms = std::chrono::duration<double, std::milli>(
        shuffled_end - shuffled_start).count();
    std::cout << "recovery-row scaling: batch=" << batch_ms
              << "ms sequential=" << sequential_ms
              << "ms shuffled=" << shuffled_ms
              << "ms deterministic_work=" << sequential_work << '/'
              << shuffled_work << std::endl;
}

void CheckPeelRowCountBounds()
{
    const wirehair_v2::PeelingCodec codec =
        wirehair_v2::MakePeelingCodec(
            wirehair_v2::PeelStructure::LtM1C16,
            wirehair_v2::PeelSolver::KsBmaxTop16);

    const std::vector<std::vector<uint16_t> > too_many_generated =
        wirehair_v2::GeneratePeelMatrixRows(
            codec, 2u, wirehair_v2::kMaxPeelMatrixRows + 1u,
            UINT64_C(0x1234));
    Check(too_many_generated.empty(),
        "peel row generation should reject unrepresentable row indexes");

    const std::vector<uint16_t> high_single =
        wirehair_v2::GeneratePeelMatrixRow(
            codec, 2u, UINT32_MAX,
            UINT64_C(0x1234));
    Check(!high_single.empty(),
        "single peel row generation should support every uint32 row index");

    std::vector<std::vector<uint16_t> > too_many_rows(
        wirehair_v2::kMaxPeelMatrixRows + 1u);
    too_many_rows[wirehair_v2::kMaxPeelMatrixRows].push_back(0u);
    const wirehair_v2::PeelEvaluation eval =
        wirehair_v2::EvaluatePeelingRows(codec, 2u, too_many_rows);
    Check(eval.TotalXorCost == UINT64_MAX &&
            eval.ResidualColumns == 2u,
        "peel evaluator should reject row counts that would wrap uint16_t");

    std::vector<std::vector<uint16_t> > invalid_columns(1u);
    invalid_columns[0].push_back(65535u);
    const wirehair_v2::PeelEvaluation invalid_eval =
        wirehair_v2::EvaluatePeelingRows(codec, 2u, invalid_columns);
    Check(invalid_eval.TotalXorCost == UINT64_MAX &&
            invalid_eval.ResidualColumns == 2u,
        "peel evaluator should reject caller rows with invalid columns");

    const wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(64000u, 1280u);
    const wirehair_v2::PeelSolvePlan plan =
        wirehair_v2::BuildPeelSolvePlan(profile, 1537u, UINT64_C(0x5511));
    Check(plan.RowCount == 0u &&
            plan.Evaluation.TotalXorCost == UINT64_MAX,
        "peel solve plan should reject N + overhead rows above uint16_t range");
}

void CheckCoreSeedHelpers()
{
    wirehair::PCGRandom prng;
    prng.Seed(1u);
    uint16_t sentinel = 0x1234u;
    wirehair::ShuffleDeck16(prng, &sentinel, 0u);
    Check(sentinel == 0x1234u, "zero-count shuffle should not touch output");

    const unsigned fixup_N = 5550u;
    const uint16_t default_dense = wirehair::GetDenseCount(fixup_N);
    Check(wirehair::GetDenseSeed(fixup_N, default_dense) == 23u,
        "default dense count should use dense seed fixup");

    const uint16_t alternate_dense =
        (default_dense >= 6u) ? (uint16_t)(default_dense - 4u) :
                                (uint16_t)(default_dense + 4u);
    Check(wirehair::GetDenseSeed(fixup_N, alternate_dense) ==
            wirehair::kDenseSeeds[alternate_dense / 4u],
        "non-default dense count should bypass per-N dense seed fixup");

    wirehair::Codec invalid_dense;
    invalid_dense.OverrideSeeds(0u, 0u, 0u);
    Check(invalid_dense.InitializeEncoder(40u, 10u) == Wirehair_InvalidInput,
        "zero dense-count override should fail initialization");
    Check(invalid_dense.InitializeEncoder(40u, 10u) == Wirehair_Success,
        "failed dense-count override should not persist");

    invalid_dense.OverrideSeeds((uint16_t)(CAT_MAX_DENSE_ROWS + 1u), 0u, 0u);
    Check(invalid_dense.InitializeEncoder(40u, 10u) == Wirehair_InvalidInput,
        "oversized dense-count override should fail initialization");
    Check(invalid_dense.InitializeEncoder(40u, 10u) == Wirehair_Success,
        "oversized dense-count override should not persist");
}

void CheckDecodeAfterAllOriginalSuccess()
{
    const uint32_t N = 4u;
    const uint32_t block_bytes = 16u;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;

    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> recovered((size_t)message_bytes, 0xa5);
    std::vector<uint8_t> repair(block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 31u + 7u);
    }

    WirehairCodec encoder = wirehair_encoder_create(
        0, &message[0], message_bytes, block_bytes);
    WirehairCodec decoder = wirehair_decoder_create(
        0, message_bytes, block_bytes);
    Check(encoder != 0, "core encoder initialization should succeed");
    Check(decoder != 0, "core decoder initialization should succeed");
    if (!encoder || !decoder) {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return;
    }

    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t block_id = 0; block_id < N; ++block_id)
    {
        result = wirehair_decode(
            decoder,
            block_id,
            &message[(size_t)block_id * block_bytes],
            block_bytes);
        if (block_id + 1u < N) {
            Check(result == Wirehair_NeedMore,
                "core decoder should need more originals before completion");
        }
    }
    Check(result == Wirehair_Success,
        "core decoder should complete on the last original block");

    uint32_t written = 0;
    Check(wirehair_encode(
            encoder, N, &repair[0], block_bytes, &written) ==
            Wirehair_Success,
        "core encoder should generate a repair block");
    Check(written == block_bytes,
        "core repair block should be a full block");
    Check(wirehair_decode(decoder, N, &repair[0], written) ==
            Wirehair_Success,
        "core decoder should ignore valid packets after completion");
    Check(wirehair_recover(decoder, &recovered[0], message_bytes) ==
            Wirehair_Success,
        "core recovery should still succeed after post-completion packet");
    Check(std::memcmp(&recovered[0], &message[0], (size_t)message_bytes) == 0,
        "core recovered message should match input");

    decoder = wirehair_decoder_create(decoder, message_bytes, block_bytes);
    Check(decoder != 0, "core decoder reuse should succeed");
    if (!decoder) {
        wirehair_free(encoder);
        return;
    }
    result = wirehair_decode(decoder, 0u, &message[0], block_bytes);
    Check(result == Wirehair_NeedMore,
        "core decoder reuse should reset completion state");

    wirehair_free(encoder);
    wirehair_free(decoder);
}

void CheckDecodeAfterDuplicateOriginalFailure()
{
    const uint32_t N = 16u;
    const uint32_t block_bytes = 16u;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;

    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> recovered((size_t)message_bytes, 0xa5);
    std::vector<uint8_t> repair(block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 17u + 3u);
    }

    WirehairCodec encoder = wirehair_encoder_create(
        0, &message[0], message_bytes, block_bytes);
    WirehairCodec decoder = wirehair_decoder_create(
        0, message_bytes, block_bytes);
    Check(encoder != 0, "dup-test encoder initialization should succeed");
    Check(decoder != 0, "dup-test decoder initialization should succeed");
    if (!encoder || !decoder) {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return;
    }

    // Feed N original-range ids with one duplicate (0 twice, N-1 missing):
    // the all-original check fails on the Nth feed before GE state exists.
    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t feed = 0; feed < N; ++feed)
    {
        const uint32_t block_id = (feed == 1u) ? 0u : feed;
        result = wirehair_decode(
            decoder,
            block_id,
            &message[(size_t)block_id * block_bytes],
            block_bytes);
    }
    Check(result == Wirehair_InvalidInput,
        "duplicate original set should be rejected as invalid input");

    // Later packets must return an error instead of resuming missing GE state
    uint32_t written = 0;
    Check(wirehair_encode(
            encoder, N + 2u, &repair[0], block_bytes, &written) ==
            Wirehair_Success,
        "dup-test encoder should generate a repair block");
    Check(wirehair_decode(decoder, N + 2u, &repair[0], written) ==
            Wirehair_InvalidInput,
        "feeds after duplicate-original failure should error, not crash");
    Check(wirehair_recover(decoder, &recovered[0], message_bytes) ==
            Wirehair_InvalidInput,
        "recover after duplicate-original failure should report invalid input");
    uint32_t bytes_out = 0x12345678u;
    Check(wirehair_recover_block(decoder, 0u, &repair[0], &bytes_out) ==
            Wirehair_InvalidInput,
        "recover_block after duplicate-original failure should report invalid input");
    Check(bytes_out == 0u,
        "recover_block after duplicate-original failure should clear byte count");

    // Reuse must clear the failed state
    decoder = wirehair_decoder_create(decoder, message_bytes, block_bytes);
    Check(decoder != 0, "dup-test decoder reuse should succeed");
    if (decoder) {
        result = Wirehair_NeedMore;
        for (uint32_t block_id = 0; block_id < N && result == Wirehair_NeedMore;
            ++block_id)
        {
            result = wirehair_decode(
                decoder,
                block_id,
                &message[(size_t)block_id * block_bytes],
                block_bytes);
        }
        Check(result == Wirehair_Success,
            "dup-test decoder reuse should decode a clean original set");
        Check(wirehair_recover(decoder, &recovered[0], message_bytes) ==
                Wirehair_Success,
            "dup-test reuse recovery should succeed");
        Check(std::memcmp(
                &recovered[0], &message[0], (size_t)message_bytes) == 0,
            "dup-test recovered message should match input");
    }

    wirehair_free(encoder);
    wirehair_free(decoder);
}

void CheckBlockBytesUpperBound()
{
    // gf256 bulk routines take int byte counts; 2^31 must be rejected
    std::vector<uint8_t> message(64u);
    WirehairCodec encoder = wirehair_encoder_create(
        0, &message[0], UINT64_C(0x100000000), 0x80000000u);
    WirehairCodec decoder = wirehair_decoder_create(
        0, UINT64_C(0x100000000), 0x80000000u);
    Check(encoder == 0,
        "block_bytes >= 2^31 should be rejected for encoder initialization");
    Check(decoder == 0,
        "block_bytes >= 2^31 should be rejected for decoder initialization");
    wirehair_free(encoder);
    wirehair_free(decoder);
}

void CheckRecoverBeforeDecodeComplete()
{
    const uint32_t N = 4u;
    const uint32_t block_bytes = 16u;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> recovered((size_t)message_bytes, 0xcc);
    std::vector<uint8_t> block(block_bytes, 0xdd);

    WirehairCodec decoder = wirehair_decoder_create(
        0, message_bytes, block_bytes);
    Check(decoder != 0, "fresh-recover decoder initialization should succeed");
    if (!decoder) {
        return;
    }

    Check(wirehair_recover(decoder, &recovered[0], message_bytes) ==
            Wirehair_NeedMore,
        "fresh decoder recover should require a completed decode");
    Check(recovered[0] == 0xcc,
        "fresh decoder recover should not write output");

    uint32_t bytes_out = 0x12345678u;
    Check(wirehair_recover_block(decoder, 0u, &block[0], &bytes_out) ==
            Wirehair_NeedMore,
        "fresh decoder recover_block should require a completed decode");
    Check(bytes_out == 0u,
        "fresh decoder recover_block should clear byte count");
    Check(block[0] == 0xdd,
        "fresh decoder recover_block should not write output");

    wirehair_free(decoder);
}

void CheckGf256RejectsNonPositiveLengths()
{
    uint8_t x[20];
    uint8_t y[20];
    uint8_t z[20];
    uint8_t x0[20];
    uint8_t y0[20];
    uint8_t z0[20];

    for (size_t i = 0; i < sizeof(x); ++i)
    {
        x[i] = (uint8_t)(0x11u + i);
        y[i] = (uint8_t)(0x80u - i);
        z[i] = (uint8_t)(0x40u + i);
    }
    std::memcpy(x0, x, sizeof(x));
    std::memcpy(y0, y, sizeof(y));
    std::memcpy(z0, z, sizeof(z));

    gf256_add_mem(x, y, -1);
    Check(std::memcmp(x, x0, sizeof(x)) == 0,
        "gf256_add_mem should ignore negative byte counts");
    gf256_add2_mem(z, x, y, -1);
    Check(std::memcmp(z, z0, sizeof(z)) == 0,
        "gf256_add2_mem should ignore negative byte counts");
    gf256_addset_mem(z, x, y, -1);
    Check(std::memcmp(z, z0, sizeof(z)) == 0,
        "gf256_addset_mem should ignore negative byte counts");
    gf256_mul_mem(z, x, 0u, -1);
    Check(std::memcmp(z, z0, sizeof(z)) == 0,
        "gf256_mul_mem should ignore negative byte counts");
    gf256_muladd_mem(z, 1u, x, -1);
    Check(std::memcmp(z, z0, sizeof(z)) == 0,
        "gf256_muladd_mem should ignore negative byte counts");
    gf256_memswap(x, y, -1);
    Check(std::memcmp(x, x0, sizeof(x)) == 0 &&
            std::memcmp(y, y0, sizeof(y)) == 0,
        "gf256_memswap should ignore negative byte counts");

    gf256_add_mem(x, y, 0);
    gf256_add2_mem(z, x, y, 0);
    gf256_addset_mem(z, x, y, 0);
    gf256_mul_mem(z, x, 0u, 0);
    gf256_muladd_mem(z, 1u, x, 0);
    gf256_memswap(x, y, 0);
    Check(std::memcmp(x, x0, sizeof(x)) == 0 &&
            std::memcmp(y, y0, sizeof(y)) == 0 &&
            std::memcmp(z, z0, sizeof(z)) == 0,
        "gf256 operations should ignore zero byte counts");
}

void CheckPeelSolverSelectionSemantics()
{
    wirehair_v2::PeelingCodec lowref =
        wirehair_v2::MakePeelingCodec(
            wirehair_v2::PeelStructure::LtM1C32,
            wirehair_v2::PeelSolver::RqccLowref);
    std::vector<std::vector<uint16_t> > lowref_rows;
    lowref_rows.push_back(std::vector<uint16_t>{0u, 2u});
    lowref_rows.push_back(std::vector<uint16_t>{1u, 3u});
    lowref_rows.push_back(std::vector<uint16_t>{1u, 2u, 3u});
    const wirehair_v2::PeelEvaluation lowref_eval =
        wirehair_v2::EvaluatePeelingRows(lowref, 4u, lowref_rows);
    Check(lowref_eval.ResidualColumns == 2u,
        "rqcc_lowref should use the largest D2 component before low-ref scoring");
    Check(lowref_eval.ResidualRows == 1u,
        "rqcc_lowref hand matrix should leave one deferred row");

    wirehair_v2::PeelingCodec ks =
        wirehair_v2::MakePeelingCodec(
            wirehair_v2::PeelStructure::LtM1C32,
            wirehair_v2::PeelSolver::KsBmaxTop16);
    std::vector<std::vector<uint16_t> > ks_rows;
    ks_rows.push_back(std::vector<uint16_t>{0u, 2u});
    ks_rows.push_back(std::vector<uint16_t>{1u, 3u});
    ks_rows.push_back(std::vector<uint16_t>{1u, 2u, 3u});
    const wirehair_v2::PeelEvaluation ks_eval =
        wirehair_v2::EvaluatePeelingRows(ks, 4u, ks_rows);
    Check(ks_eval.ResidualColumns == 2u,
        "ks_bmax_top16 should use the largest D2 component before boundary scoring");
    Check(ks_eval.ResidualRows == 1u,
        "ks_bmax_top16 hand matrix should leave one deferred row");

    std::vector<std::vector<uint16_t> > deferred_rows;
    deferred_rows.push_back(std::vector<uint16_t>{0u, 1u});
    deferred_rows.push_back(std::vector<uint16_t>{1u, 2u});
    deferred_rows.push_back(std::vector<uint16_t>{0u, 1u, 3u});
    deferred_rows.push_back(std::vector<uint16_t>{0u, 1u, 4u});
    const wirehair_v2::PeelEvaluation deferred_eval =
        wirehair_v2::EvaluatePeelingRows(lowref, 5u, deferred_rows);
    Check(deferred_eval.ResidualColumns == 1u,
        "peel evaluator should solve the deferred-row hand matrix with one inactivation");
    Check(deferred_eval.ResidualRows == 0u,
        "peel evaluator should not count rows solved by the peel cascade as residual");
}

} // namespace

int main()
{
    using namespace wirehair_v2;

    Check(wirehair_init() == Wirehair_Success, "wirehair_init should succeed");
    CheckCoreSeedHelpers();
    CheckDecodeAfterAllOriginalSuccess();
    CheckDecodeAfterDuplicateOriginalFailure();
    CheckBlockBytesUpperBound();
    CheckRecoverBeforeDecodeComplete();
    CheckGf256RejectsNonPositiveLengths();
    CheckV2InvalidInitialization();
    CheckV2EncodeFailureClearsByteCount();

    Check(ClassifyBlockBytes(1280u) == BlockByteClass::Small,
        "1280-byte blocks should use small byte class");
    Check(ClassifyBlockBytes(100u * 1024u) == BlockByteClass::Medium,
        "100 KiB blocks should use medium byte class");
    Check(ClassifyBlockBytes(1024u * 1024u) == BlockByteClass::Large,
        "1 MiB blocks should use large byte class");

    Check(ClassifyBlockCount(1000u) == BlockCountBand::UpTo1000,
        "N=1000 should use first count band");
    Check(ClassifyBlockCount(3200u) == BlockCountBand::UpTo3200,
        "N=3200 should use second count band");
    Check(ClassifyBlockCount(12000u) == BlockCountBand::UpTo12000,
        "N=12000 should use third count band");
    Check(ClassifyBlockCount(12001u) == BlockCountBand::Above12000,
        "N=12001 should use final count band");

    CheckPolicy(1000u, 1280u, BlockByteClass::Small,
        BlockCountBand::UpTo1000, PeelStructure::RobustD1_001D2_012,
        DegreeFamily::RobustD1D2, 1u, 64u);
    CheckPolicy(3200u, 1280u, BlockByteClass::Small,
        BlockCountBand::UpTo3200, PeelStructure::LtM1C32,
        DegreeFamily::Lt, 1u, 32u);
    CheckPolicy(6400u, 1280u, BlockByteClass::Small,
        BlockCountBand::UpTo12000, PeelStructure::RobustD1_001D2_003,
        DegreeFamily::RobustD1D2, 1u, 64u);

    const PeelingCodec robust001_012 =
        MakePeelingCodec(PeelStructure::RobustD1_001D2_012,
            PeelSolver::KsBmaxTop16);
    Check(robust001_012.Degree1Mass == 0.01,
        "robust_d1_001 should mean 1 percent degree-1 mass");
    Check(robust001_012.Degree2Mass == 0.12,
        "robust d2_012 should mean 12 percent degree-2 mass");
    CheckPolicy(32000u, 1280u, BlockByteClass::Small,
        BlockCountBand::Above12000, PeelStructure::LtM2C96,
        DegreeFamily::Lt, 2u, 96u);

    CheckPolicy(1000u, 100u * 1024u, BlockByteClass::Medium,
        BlockCountBand::UpTo1000, PeelStructure::LtM1C16,
        DegreeFamily::Lt, 1u, 16u);
    CheckPolicy(3200u, 100u * 1024u, BlockByteClass::Medium,
        BlockCountBand::UpTo3200, PeelStructure::LtM1C32,
        DegreeFamily::Lt, 1u, 32u);
    CheckPolicy(6400u, 100u * 1024u, BlockByteClass::Medium,
        BlockCountBand::UpTo12000, PeelStructure::LtM1C64,
        DegreeFamily::Lt, 1u, 64u);
    CheckPolicy(32000u, 100u * 1024u, BlockByteClass::Medium,
        BlockCountBand::Above12000, PeelStructure::RsC001D50C128,
        DegreeFamily::RobustSoliton, 1u, 128u);

    CheckPolicy(1000u, 1024u * 1024u, BlockByteClass::Large,
        BlockCountBand::UpTo1000, PeelStructure::LtM1C16,
        DegreeFamily::Lt, 1u, 16u);
    CheckPolicy(3200u, 1024u * 1024u, BlockByteClass::Large,
        BlockCountBand::UpTo3200, PeelStructure::LtM1C32,
        DegreeFamily::Lt, 1u, 32u);
    CheckPolicy(6400u, 1024u * 1024u, BlockByteClass::Large,
        BlockCountBand::UpTo12000, PeelStructure::LtM1C64,
        DegreeFamily::Lt, 1u, 64u);
    CheckPolicy(32000u, 1024u * 1024u, BlockByteClass::Large,
        BlockCountBand::Above12000, PeelStructure::RsC003D10C128,
        DegreeFamily::RobustSoliton, 1u, 128u);

    const PeelingCodec lowref =
        MakePeelingCodec(PeelStructure::LtM1C32, PeelSolver::RqccLowref);
    Check(lowref.SolverCandidateLimit == 0u,
        "rqcc_lowref should not use a top-K boundary candidate limit");
    CheckPeelSolverSelectionSemantics();

    Check(ToString(PeelSolver::KsBmaxTop16) == std::string("ks_bmax_top16"),
        "solver string mismatch");
    Check(ToString(PeelStructure::LtM2C96) == std::string("lt_m2_c96"),
        "structure string mismatch");
    Check(ToString(DegreeFamily::RobustSoliton) ==
            std::string("robust_soliton"),
        "degree family string mismatch");

    const SeedProfile profile3061 = SelectSeedProfile(3061u, 1280u);
    Check(profile3061.PeelSeedBucket == 1013u,
        "N=3061 should map to its modulo-2048 bucket");
    Check(profile3061.UsedPeelFixup,
        "N=3061 should use the exact-N peel seed fixup");

    const SeedProfile profile7533 = SelectSeedProfile(7533u, 1280u);
    Check(profile7533.UsedDenseFixup,
        "N=7533 should use the exact-N dense seed fixup");

    Check(CandidatePeelSeed(17u, 5u, 0u) == 5u,
        "candidate peel seed 0 should preserve base seed");
    Check(CandidatePeelSeed(17u, 5u, 1u) ==
            CandidatePeelSeed(17u, 5u, 1u),
        "candidate peel seed should be deterministic");

    const PeelingCodec eval_codec =
        MakePeelingCodec(PeelStructure::LtM1C16, PeelSolver::KsBmaxTop16);
    const PeelStructure all_structures[] = {
        PeelStructure::LtM1C16,
        PeelStructure::LtM1C32,
        PeelStructure::LtM1C64,
        PeelStructure::LtM2C96,
        PeelStructure::RobustD1_001D2_003,
        PeelStructure::RobustD1_001D2_012,
        PeelStructure::RsC001D50C128,
        PeelStructure::RsC003D10C128
    };
    for (size_t i = 0; i < sizeof(all_structures) / sizeof(all_structures[0]); ++i)
    {
        CheckGeneratedRows(
            MakePeelingCodec(all_structures[i], PeelSolver::KsBmaxTop16),
            96u, 97u,
            UINT64_C(0x1234) ^ ((uint64_t)i * UINT64_C(0x9e3779b97f4a7c15)));
    }
    CheckGeneratedRecoveryRows();
    CheckRecoveryRowScaling();
    CheckPeelRowCountBounds();

    const PeelEvaluation eval_a =
        EvaluatePeeling(eval_codec, 64u, UINT64_C(0x998877));
    const PeelEvaluation eval_b =
        EvaluatePeeling(eval_codec, 64u, UINT64_C(0x998877));
    Check(eval_a.ResidualColumns == eval_b.ResidualColumns,
        "peel evaluation should be deterministic");
    Check(eval_a.MatrixRefs >= eval_a.Rows,
        "peel matrix refs should cover every row");
    Check(eval_a.TotalXorCost >= eval_a.MatrixXors,
        "total xor estimate should include matrix xors");

    const PeelSolvePlan plan_a =
        BuildPeelSolvePlan(SelectSeedProfile(96u, 1280u), 0u, UINT64_C(0x5511));
    const PeelSolvePlan plan_b =
        BuildPeelSolvePlan(SelectSeedProfile(96u, 1280u), 0u, UINT64_C(0x5511));
    Check(plan_a.MatrixSeed == plan_b.MatrixSeed,
        "peel solve plan seed should be deterministic");
    Check(plan_a.Rows.size() == 96u, "H=0 plan should generate N rows");
    Check(plan_a.Evaluation.Columns == 96u,
        "peel solve plan should evaluate N columns");

    const PeelSolvePlan plan_h2 =
        BuildPeelSolvePlan(SelectSeedProfile(96u, 1280u), 2u, UINT64_C(0x5511));
    Check(plan_h2.MatrixSeed != plan_a.MatrixSeed,
        "peel solve plan seed should include row count");
    Check(plan_h2.Rows.size() == 98u,
        "overhead plan should generate N + H rows");

    SeedProfile altered_seed = SelectSeedProfile(96u, 1280u);
    altered_seed.PeelSeed =
        (uint16_t)((altered_seed.PeelSeed + 1u) & 0xffu);
    Check(MatrixSeedFromProfile(altered_seed, 96u, UINT64_C(0x5511)) !=
            plan_a.MatrixSeed,
        "peel solve plan seed should include selected peel seed");

    SeedTuningOptions tuning = DefaultSeedTuningOptions();
    tuning.PeelCandidates = 4u;
    tuning.TrialsPerCandidate = 1u;
    const SeedProfile base80 = SelectSeedProfile(80u, 1280u);
    const SeedProfile tuned = TuneSeedProfile(80u, 1280u, tuning);
    Check(tuned.Tuned, "seed tuner should mark tuned profiles");
    Check(tuned.TuningXorCost > 0u, "seed tuner should report xor cost");
    if (tuned.PeelSeed != base80.PeelSeed) {
        Check(tuned.UsedPeelFixup,
            "seed tuner should update peel fixup metadata after tuning");
    }

    CheckV2RoundTrip();
    CheckV2PrecodeEncoderFacade();

    if (Failures != 0) {
        return 1;
    }

    return 0;
}
