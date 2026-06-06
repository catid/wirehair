#include "WirehairV2Policy.h"
#include "WirehairV2Codec.h"
#include "WirehairV2Peel.h"
#include "WirehairV2Seeds.h"

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

} // namespace

int main()
{
    using namespace wirehair_v2;

    Check(wirehair_init() == Wirehair_Success, "wirehair_init should succeed");

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
    const std::vector<std::vector<uint16_t> > rows =
        GeneratePeelMatrixRows(eval_codec, 48u, 48u, UINT64_C(0x1234));
    Check(rows.size() == 48u, "generated matrix should contain requested rows");
    for (size_t i = 0; i < rows.size(); ++i)
    {
        Check(!rows[i].empty(), "generated row should not be empty");
        Check(rows[i].size() <= 16u, "generated row should obey degree cap");
        CheckRowHasNoDuplicates(rows[i]);
    }

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

    SeedTuningOptions tuning = DefaultSeedTuningOptions();
    tuning.PeelCandidates = 4u;
    tuning.TrialsPerCandidate = 1u;
    const SeedProfile tuned = TuneSeedProfile(80u, 1280u, tuning);
    Check(tuned.Tuned, "seed tuner should mark tuned profiles");
    Check(tuned.TuningXorCost > 0u, "seed tuner should report xor cost");

    CheckV2RoundTrip();

    if (Failures != 0) {
        return 1;
    }

    return 0;
}
