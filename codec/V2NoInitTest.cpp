#include "WirehairV2Codec.h"

#include <cstdio>
#include <cstring>
#include <new>
#include <vector>

namespace {

bool IsDefaultProfile(const wirehair_v2::SeedProfile& profile)
{
    return profile.BlockCount == 0u && profile.BlockBytes == 0u &&
        profile.DenseCount == 0u && profile.PeelSeed == 0u &&
        profile.DenseSeed == 0u && profile.PeelSeedBucket == 0u &&
        !profile.UsedPeelFixup && !profile.UsedDenseFixup &&
        !profile.V2SeedSelected && profile.V2SeedAttempt == 0u &&
        profile.V2PrecodeContractVersion == 0u &&
        profile.V2PacketRowContractVersion == 0u &&
        profile.V2StaircaseCount == 0u &&
        profile.V2DenseRowCount == 0u &&
        profile.V2HeavyRowCount == 0u &&
        profile.V2CompletionField ==
            wirehair_v2::CompletionField::GF256 &&
        profile.V2SourceHits == 0u && profile.V2PrecodeSeed == 0u &&
        profile.V2PacketPeelSeed == 0u &&
        profile.V2RecoveryMixCount == 0u &&
        !profile.V2DenseIdentityCorner &&
        !profile.V2DenseTwoAnchor &&
        !profile.V2AdaptiveDenseTwoAnchor &&
        profile.V2PrecodeSeedSalt == 0u &&
        profile.V2RecoveryRowSeedSalt == 0u &&
        !profile.Tuned &&
        profile.TuningResidualMean == 0.0 &&
        profile.TuningResidualColumns == 0u &&
        profile.TuningXorCost == 0u &&
        profile.TuningCandidatesRequested == 0u &&
        profile.TuningCandidatesUnique == 0u &&
        profile.TuningCandidatesCompleted == 0u &&
        profile.TuningTrials == 0u &&
        profile.Policy.Solver == wirehair_v2::PeelSolver::RqccLowref &&
        profile.Policy.Structure == wirehair_v2::PeelStructure::LtM1C16 &&
        profile.Policy.ByteClass == wirehair_v2::BlockByteClass::Small &&
        profile.Policy.CountBand == wirehair_v2::BlockCountBand::UpTo1000 &&
        profile.Policy.Codec.MinDegree == 0u &&
        profile.Policy.Codec.MaxDegree == 0u &&
        profile.Policy.Codec.Degree1Mass == 0.0 &&
        profile.Policy.Codec.Degree2Mass == 0.0 &&
        profile.Policy.Codec.RobustC == 0.0 &&
        profile.Policy.Codec.RobustDelta == 0.0 &&
        !profile.Policy.Codec.FullyRandomRows;
}

bool CheckPatternInitializedAccessors()
{
    using wirehair_v2::Codec;
    using wirehair_v2::MessagePrecodeDecoder;
    using wirehair_v2::MessagePrecodeEncoder;
    using wirehair_v2::PrecodeEncoder;

    alignas(Codec) unsigned char codec_storage[sizeof(Codec)];
    std::memset(codec_storage, 0xfe, sizeof(codec_storage));
    Codec* codec = new (codec_storage) Codec;
    if (!IsDefaultProfile(codec->Profile())) {
        std::fprintf(stderr, "pattern Codec profile was not value-initialized\n");
        codec->~Codec();
        return false;
    }
    codec->~Codec();

    alignas(PrecodeEncoder)
        unsigned char encoder_storage[sizeof(PrecodeEncoder)];
    std::memset(encoder_storage, 0xfe, sizeof(encoder_storage));
    PrecodeEncoder* encoder = new (encoder_storage) PrecodeEncoder;
    const wirehair_v2::PrecodeSystem& system = encoder->System();
    const wirehair_v2::PrecodeEncodeStats& stats = encoder->EncodeStats();
    if (encoder->IsInitialized() || encoder->SourceBlockCount() != 0u ||
        encoder->ParityBlockCount() != 0u || encoder->BlockBytes() != 0u ||
        encoder->RecoveryRowSeed() != 0u ||
        encoder->RecoveryMixCount() != 0u ||
        encoder->ParityBlocks() != nullptr ||
        encoder->IntermediateBlocks() != nullptr ||
        system.Params.BlockCount != 0u || system.Params.Staircase != 0u ||
        system.Params.DenseRows != 0u || system.Params.HeavyRows != 0u ||
        system.Params.SourceHits != 0u || system.Params.DenseIdentityCorner ||
        system.Params.DenseTwoAnchor ||
        system.Params.Seed != 0u || !system.StaircaseRows.empty() ||
        !system.DenseRowColumns.empty() || stats.StaircaseBlockOps != 0u ||
        stats.DenseKnownBlockOps != 0u || stats.DenseSolveBlockOps != 0u ||
        stats.HeavyBucketXors != 0u || stats.HeavyMulAdds != 0u ||
        stats.HeavySolveBlockOps != 0u)
    {
        std::fprintf(stderr,
            "pattern PrecodeEncoder accessors were not value-initialized\n");
        encoder->~PrecodeEncoder();
        return false;
    }
    encoder->~PrecodeEncoder();

    alignas(MessagePrecodeEncoder)
        unsigned char message_storage[sizeof(MessagePrecodeEncoder)];
    std::memset(message_storage, 0xfe, sizeof(message_storage));
    MessagePrecodeEncoder* message_encoder =
        new (message_storage) MessagePrecodeEncoder;
    const wirehair_v2::MessagePrecodeEncoderOptions& options =
        message_encoder->Options();
    if (message_encoder->IsInitialized() ||
        message_encoder->MessageBytes() != 0u ||
        message_encoder->SourceBlockCount() != 0u ||
        message_encoder->BlockBytes() != 0u ||
        message_encoder->IntermediateBlocks() != nullptr ||
        !IsDefaultProfile(message_encoder->Profile()) ||
        options.RecoveryMixCount != wirehair_v2::kDefaultRecoveryMixCount ||
        options.DenseIdentityCorner ||
        options.AdaptiveDenseTwoAnchor ||
        options.PrecodeSeedSalt != wirehair_v2::kMessagePrecodeSeedSalt ||
        options.RecoveryRowSeedSalt !=
            wirehair_v2::kMessageRecoveryRowSeedSalt ||
        message_encoder->BlockEncoder().IsInitialized())
    {
        std::fprintf(stderr,
            "pattern MessagePrecodeEncoder accessors were not initialized\n");
        message_encoder->~MessagePrecodeEncoder();
        return false;
    }
    message_encoder->~MessagePrecodeEncoder();

    alignas(MessagePrecodeDecoder)
        unsigned char decoder_storage[sizeof(MessagePrecodeDecoder)];
    std::memset(decoder_storage, 0xfe, sizeof(decoder_storage));
    MessagePrecodeDecoder* message_decoder =
        new (decoder_storage) MessagePrecodeDecoder;
    const wirehair_v2::MessagePrecodeEncoderOptions& decoder_options =
        message_decoder->Options();
    const wirehair_v2::PrecodeSystem& decoder_system =
        message_decoder->System();
    const wirehair_v2::PrecodeSolveStats& solve_stats =
        message_decoder->SolveStats();
    if (message_decoder->IsInitialized() || message_decoder->IsDecoded() ||
        message_decoder->ReceivedCount() != 0u ||
        message_decoder->SolveAttemptCount() != 0u ||
        message_decoder->MessageBytes() != 0u ||
        message_decoder->BlockBytes() != 0u ||
        message_decoder->IntermediateBlocks() != nullptr ||
        !IsDefaultProfile(message_decoder->Profile()) ||
        decoder_options.RecoveryMixCount !=
            wirehair_v2::kDefaultRecoveryMixCount ||
        decoder_options.DenseIdentityCorner ||
        decoder_options.AdaptiveDenseTwoAnchor ||
        decoder_options.PrecodeSeedSalt !=
            wirehair_v2::kMessagePrecodeSeedSalt ||
        decoder_options.RecoveryRowSeedSalt !=
            wirehair_v2::kMessageRecoveryRowSeedSalt ||
        decoder_system.Params.BlockCount != 0u ||
        !decoder_system.StaircaseRows.empty() ||
        !decoder_system.DenseRowColumns.empty() ||
        solve_stats.PacketRows != 0u ||
        solve_stats.InactivatedColumns != 0u ||
        solve_stats.ResidualRank != 0u)
    {
        std::fprintf(stderr,
            "pattern MessagePrecodeDecoder accessors were not initialized\n");
        message_decoder->~MessagePrecodeDecoder();
        return false;
    }
    message_decoder->~MessagePrecodeDecoder();
    return true;
}

} // namespace

int main()
{
    if (!CheckPatternInitializedAccessors()) {
        return 1;
    }
    const uint32_t N = 4u;
    const uint32_t block_bytes = 16u;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;

    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> recovered((size_t)message_bytes, 0);
    std::vector<uint8_t> block(block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 13u + 7u);
    }

    wirehair_v2::Codec encoder;
    wirehair_v2::Codec decoder;
    WirehairResult result = encoder.InitializeEncoder(
        &message[0], message_bytes, block_bytes);
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "encoder init failed: %d\n", result);
        return 1;
    }
    result = decoder.InitializeDecoder(message_bytes, block_bytes);
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "decoder init failed: %d\n", result);
        return 1;
    }

    for (uint32_t block_id = 0; block_id < N; ++block_id)
    {
        uint32_t written = 0;
        result = encoder.Encode(block_id, &block[0], block_bytes, &written);
        if (result != Wirehair_Success) {
            std::fprintf(stderr, "encode failed: %d\n", result);
            return 1;
        }
        result = decoder.Decode(block_id, &block[0], written);
        if (block_id + 1u < N && result != Wirehair_NeedMore) {
            std::fprintf(stderr, "decode ended early: %d\n", result);
            return 1;
        }
    }
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "decode did not complete: %d\n", result);
        return 1;
    }
    result = decoder.Recover(&recovered[0], message_bytes);
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "recover failed: %d\n", result);
        return 1;
    }
    if (std::memcmp(&recovered[0], &message[0], (size_t)message_bytes) != 0) {
        std::fprintf(stderr, "recovered message mismatch\n");
        return 1;
    }
    return 0;
}
