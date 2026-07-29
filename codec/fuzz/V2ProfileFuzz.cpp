#include "V2FuzzDriver.h"

#include "../V2TinyDenseOracle.h"
#include "../WirehairV2PrecodeEncode.h"
#include "../WirehairV2TinyMDS.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#ifndef WIREHAIR_FUZZ_TARGET_NAME
#define WIREHAIR_FUZZ_TARGET_NAME "profile"
#endif

#ifndef WIREHAIR_FUZZ_CORPUS_MANIFEST
#define WIREHAIR_FUZZ_CORPUS_MANIFEST "codec/fuzz/corpus/profile/manifest.txt"
#endif

namespace {

bool SameParams(
    const wirehair_v2::PrecodeParams& a,
    const wirehair_v2::PrecodeParams& b)
{
    return a.BlockCount == b.BlockCount &&
        a.Staircase == b.Staircase &&
        a.DenseRows == b.DenseRows &&
        a.HeavyRows == b.HeavyRows &&
        a.SourceHits == b.SourceHits &&
        a.Field == b.Field &&
        a.HeavyFamily == b.HeavyFamily &&
        a.DegreeBalancedStaircase == b.DegreeBalancedStaircase &&
        a.StaircaseDegreeScale == b.StaircaseDegreeScale &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.DenseTwoAnchor == b.DenseTwoAnchor &&
        a.DenseTwoAnchorPhase == b.DenseTwoAnchorPhase &&
        a.SegmentedDenseAnchors == b.SegmentedDenseAnchors &&
        a.Seed == b.Seed;
}

bool SameSystem(
    const wirehair_v2::PrecodeSystem& a,
    const wirehair_v2::PrecodeSystem& b)
{
    return SameParams(a.Params, b.Params) &&
        a.StaircaseRows == b.StaircaseRows &&
        a.DenseRowColumns == b.DenseRowColumns;
}

bool Fail(std::string& failure, const char* message)
{
    failure = message;
    return false;
}

bool FuzzParams(wirehair_v2::fuzz::Input& input, std::string& failure)
{
    const uint32_t K = 2u + input.U8() % 127u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, input.U64());
    const unsigned mutation = input.U8() % 22u;
    bool expected_valid = false;
    switch (mutation)
    {
    case 0: params.BlockCount = 0u; break;
    case 1: params.BlockCount = 1u; break;
    case 2: params.BlockCount = 64001u; break;
    case 3: params.BlockCount = UINT32_MAX; break;
    case 4: params.Staircase = 0u; break;
    case 5: params.Staircase = UINT32_MAX; break;
    case 6: params.DenseRows = 65u; break;
    case 7: params.HeavyRows = 129u; break;
    case 8: params.SourceHits = 0u; break;
    case 9: params.SourceHits = 9u; break;
    case 10:
        params.DenseIdentityCorner = true;
        expected_valid = (uint64_t)params.BlockCount + params.Staircase >=
            2u * (uint64_t)(params.DenseRows >> 1);
        break;
    case 11:
        params.Seed ^= input.U64() | 1u;
        expected_valid = true;
        break;
    case 12:
        expected_valid = true;
        break;
    case 13:
        params.DenseRows = 64u;
        expected_valid = true;
        break;
    case 14:
        params.Field = static_cast<wirehair_v2::CompletionField>(UINT32_MAX);
        break;
    case 15:
        params.DegreeBalancedStaircase = true;
        expected_valid = true;
        break;
    case 16:
        params.StaircaseDegreeScale = 0.0;
        expected_valid = true;
        break;
    case 17:
        params.StaircaseDegreeScale =
            wirehair_v2::kStaircaseDegreeScaleMax + 1.0;
        break;
    case 18:
        params.StaircaseDegreeScale =
            std::numeric_limits<double>::quiet_NaN();
        break;
    case 19:
        params.StaircaseDegreeScale = -2.0;
        break;
    case 20:
        params.StaircaseDegreeScale =
            -std::numeric_limits<double>::infinity();
        break;
    default:
        params = wirehair_v2::MakeMixedParams(K, params.Seed);
        expected_valid = true;
        break;
    }

    wirehair_v2::PrecodeSystem sentinel;
    sentinel.Params.BlockCount = 7u;
    sentinel.Params.Staircase = 3u;
    sentinel.StaircaseRows.push_back(std::vector<uint32_t>{1u, 2u});
    sentinel.DenseRowColumns.push_back(std::vector<uint32_t>{3u});
    wirehair_v2::PrecodeSystem output = sentinel;
    const bool built = wirehair_v2::BuildPrecodeSystem(params, output);
    if (built != expected_valid) {
        return Fail(failure, "precode parameter acceptance mismatch");
    }
    if (!built) {
        return SameSystem(output, sentinel) ? true :
            Fail(failure, "invalid parameters modified output system");
    }
    if (!wirehair_v2::ValidatePrecodeSystem(output)) {
        return Fail(failure, "built precode system did not validate");
    }
    wirehair_v2::PrecodeSystem rebuilt;
    if (!wirehair_v2::BuildPrecodeSystem(params, rebuilt) ||
        !SameSystem(output, rebuilt))
    {
        return Fail(failure, "precode construction was not deterministic");
    }
    return true;
}

bool FuzzSystem(wirehair_v2::fuzz::Input& input, std::string& failure)
{
    const uint32_t K = 16u + input.U8() % 17u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(K, input.U64()), system))
    {
        return Fail(failure, "valid system construction failed");
    }
    const unsigned mutation = input.U8() % 9u;
    switch (mutation)
    {
    case 0: break;
    case 1: ++system.Params.BlockCount; break;
    case 2: system.StaircaseRows.pop_back(); break;
    case 3:
        system.StaircaseRows[0].push_back(UINT32_MAX);
        break;
    case 4:
        system.StaircaseRows[0].push_back(system.StaircaseRows[0].back());
        break;
    case 5:
        std::reverse(
            system.DenseRowColumns[0].begin(),
            system.DenseRowColumns[0].end());
        break;
    case 6: ++system.Params.DenseRows; break;
    case 7: system.Params.HeavyRows = 129u; break;
    default: system.Params.SourceHits = 0u; break;
    }
    const bool expected_valid = mutation == 0u;
    if (wirehair_v2::ValidatePrecodeSystem(system) != expected_valid) {
        return Fail(failure, "mutated system validation mismatch");
    }

    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = input.U32();
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    std::vector<wirehair_v2::SolvePacket> packets;
    std::vector<uint8_t> solved(9u, 0xabu);
    const std::vector<uint8_t> before = solved;
    const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
        system, config, packets, 1u, solved);
    const WirehairResult expected_result = expected_valid ?
        Wirehair_NeedMore : Wirehair_InvalidInput;
    if (result != expected_result || solved != before) {
        return Fail(failure, "mutated system solve/no-write mismatch");
    }
    return true;
}

bool FuzzProfileContract(
    wirehair_v2::fuzz::Input& input,
    std::string& failure)
{
    const uint32_t K = (input.U8() & 31u) == 0u ?
        wirehair_v2::kDenseTwoAnchorMinBlockCount + input.U8() % 17u :
        16u + input.U8() % 17u;
    uint32_t block_bytes = 1u + input.U8() % 64u;
    const bool mixed = input.Bool();
    if (mixed) block_bytes += block_bytes & 1u;
    wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.PrecodeSeedSalt = input.U64();
    options.RecoveryRowSeedSalt = input.U64();
    options.Completion = mixed ?
        wirehair_v2::CompletionField::MixedGF256GF16 :
        wirehair_v2::CompletionField::GF256;
    const bool small_band_d4 = mixed && input.Bool();
    options.Architecture = small_band_d4 ?
        wirehair_v2::V2PrecodeArchitecture::SmallBandD4 :
        wirehair_v2::V2PrecodeArchitecture::LegacyD12;
    options.DenseIdentityCorner = !small_band_d4 && input.Bool();
    options.RecoveryMixCount =
        !small_band_d4 && mixed && input.Bool() ? 2u :
            wirehair_v2::kCertifiedPacketMixCount;
    options.AdaptiveDenseTwoAnchor = !small_band_d4 && mixed &&
        options.RecoveryMixCount == 2u &&
        !options.DenseIdentityCorner && input.Bool();

    wirehair_v2::MessagePrecodeEncoderOptions resolved;
    wirehair_v2::PrecodeParams params;
    wirehair_v2::PacketRowConfig config;
    if (!wirehair_v2::ResolveMessagePrecodeOptions(
            profile, &options, resolved) ||
        !wirehair_v2::ResolveMessagePrecodeConfiguration(
            profile, resolved, params, config))
    {
        return Fail(failure, "valid unselected profile did not resolve");
    }
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
        return Fail(failure, "resolved profile system did not build");
    }
    wirehair_v2::BindMessagePrecodeProfile(
        profile, resolved, system, config, 0u);
    wirehair_v2::MessagePrecodeEncoderOptions bound_options;
    wirehair_v2::PrecodeParams bound_params;
    wirehair_v2::PacketRowConfig bound_config;
    if (!wirehair_v2::ResolveMessagePrecodeOptions(
            profile, nullptr, bound_options) ||
        !wirehair_v2::ResolveMessagePrecodeConfiguration(
            profile, bound_options, bound_params, bound_config) ||
        !SameParams(params, bound_params) ||
        config.PeelSeed != bound_config.PeelSeed ||
        config.MixCount != bound_config.MixCount)
    {
        return Fail(failure, "bound profile did not round-trip");
    }

    wirehair_v2::SeedProfile bad = profile;
    wirehair_v2::MessagePrecodeEncoderOptions requested = bound_options;
    const unsigned mutation = input.U8() % 24u;
    switch (mutation)
    {
    case 0: ++bad.V2PrecodeContractVersion; break;
    case 1: ++bad.V2PacketRowContractVersion; break;
    case 2: ++bad.V2StaircaseCount; break;
    case 3: ++bad.V2DenseRowCount; break;
    case 4: ++bad.V2HeavyRowCount; break;
    case 5:
        bad.V2CompletionField =
            static_cast<wirehair_v2::CompletionField>(UINT32_MAX);
        break;
    case 6: ++bad.V2SourceHits; break;
    case 7: bad.V2PrecodeSeed ^= 1u; break;
    case 8: bad.V2PacketPeelSeed ^= 1u; break;
    case 9: bad.V2RecoveryMixCount = 0u; break;
    case 10: bad.V2DenseIdentityCorner = !bad.V2DenseIdentityCorner; break;
    case 11: bad.V2DenseTwoAnchor = !bad.V2DenseTwoAnchor; break;
    case 12:
        bad.V2AdaptiveDenseTwoAnchor = !bad.V2AdaptiveDenseTwoAnchor;
        break;
    case 13: bad.V2PrecodeSeedSalt ^= 1u; break;
    case 14: bad.V2RecoveryRowSeedSalt ^= 1u; break;
    case 15: bad.V2SeedAttempt = wirehair_v2::kMaxPacketSeedAttempts; break;
    case 16: bad.DenseCount ^= 1u; break;
    case 17:
        bad.V2SeedSelected = false;
        break;
    case 18:
        bad.V2Architecture =
            static_cast<wirehair_v2::V2PrecodeArchitecture>(UINT32_MAX);
        break;
    case 19:
        bad.V2SeedPolicy =
            static_cast<wirehair_v2::V2SeedDerivation>(UINT32_MAX);
        break;
    case 20:
        requested.PrecodeSeedSalt ^= 1u;
        break;
    case 21:
        requested.Completion = requested.Completion ==
                wirehair_v2::CompletionField::MixedGF256GF16 ?
            wirehair_v2::CompletionField::GF256 :
            wirehair_v2::CompletionField::MixedGF256GF16;
        break;
    case 22:
        requested.Architecture =
            requested.Architecture ==
                    wirehair_v2::V2PrecodeArchitecture::LegacyD12 ?
                wirehair_v2::V2PrecodeArchitecture::SmallBandD4 :
                wirehair_v2::V2PrecodeArchitecture::LegacyD12;
        break;
    default:
        requested.AdaptiveDenseTwoAnchor =
            !requested.AdaptiveDenseTwoAnchor;
        break;
    }
    wirehair_v2::MessagePrecodeEncoderOptions ignored_options;
    const wirehair_v2::MessagePrecodeEncoderOptions* requested_pointer =
        mutation >= 20u ? &requested : nullptr;
    if (wirehair_v2::ResolveMessagePrecodeOptions(
            bad, requested_pointer, ignored_options))
    {
        wirehair_v2::PrecodeParams ignored_params;
        wirehair_v2::PacketRowConfig ignored_config;
        if (wirehair_v2::ResolveMessagePrecodeConfiguration(
                bad, ignored_options, ignored_params, ignored_config))
        {
            return Fail(failure, "malformed bound profile was accepted");
        }
    }
    return true;
}

bool FuzzPacketDomain(wirehair_v2::fuzz::Input& input, std::string& failure)
{
    const uint32_t K = 2u + input.U8() % 31u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(K, input.U64()), system))
    {
        return Fail(failure, "packet-domain system build failed");
    }
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = input.U32();
    config.MixCount = input.U8() % 6u;
    const uint32_t block_id = input.U32();
    const uint32_t block_bytes = 1u + input.U8() % 64u;
    const bool valid = wirehair_v2::IsPacketRowDomainValid(
        K, P, config.MixCount);
    const std::vector<uint32_t> row =
        wirehair_v2::GeneratePacketMatrixRow(K, P, block_id, config);
    if (row.empty() == valid) {
        return Fail(failure, "packet row domain/generation mismatch");
    }

    std::vector<uint8_t> intermediate((size_t)(K + P) * block_bytes);
    for (size_t i = 0; i < intermediate.size(); ++i) {
        intermediate[i] = (uint8_t)(i * 29u + block_id * 7u + 3u);
    }
    std::vector<uint8_t> output(block_bytes + 32u, 0xa5u);
    const std::vector<uint8_t> before = output;
    uint64_t operations = UINT64_C(0xfeedfacecafebeef);
    const bool evaluated = wirehair_v2::EvaluatePacketBlock(
        system, config, intermediate.data(), block_bytes, block_id,
        output.data() + 16u, &operations);
    if (evaluated != valid) {
        return Fail(failure, "packet evaluation domain mismatch");
    }
    if (!valid) {
        return output == before &&
            operations == UINT64_C(0xfeedfacecafebeef) ? true :
            Fail(failure, "rejected packet evaluation wrote output");
    }
    std::vector<uint8_t> expected(block_bytes, 0u);
    for (uint32_t column : row) {
        for (uint32_t b = 0; b < block_bytes; ++b) {
            expected[b] ^= intermediate[(size_t)column * block_bytes + b];
        }
    }
    if (operations != row.size() ||
        std::memcmp(output.data() + 16u, expected.data(), block_bytes) != 0 ||
        !std::equal(output.begin(), output.begin() + 16u, before.begin()) ||
        !std::equal(
            output.begin() + 16u + block_bytes,
            output.end(),
            before.begin() + 16u + block_bytes))
    {
        return Fail(failure, "packet evaluation oracle/guard mismatch");
    }
    return true;
}

bool FuzzTinyOracle(wirehair_v2::fuzz::Input& input, std::string& failure)
{
    const uint32_t K = 2u + input.U8() % 7u;
    const uint32_t block_bytes = 1u + input.U8() % 8u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(K, input.U64()), system))
    {
        return Fail(failure, "oracle system build failed");
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = input.U32();
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    const unsigned stream = input.U8() % 5u;
    const uint32_t packet_count = stream == 0u ? K - 1u :
        K + (stream >= 3u ? 2u : 0u);
    std::vector<uint8_t> packet_data((size_t)packet_count * block_bytes);
    for (size_t i = 0; i < packet_data.size(); ++i) {
        packet_data[i] = (uint8_t)(input.U8() + i * 17u);
    }
    std::vector<wirehair_v2::SolvePacket> packets(packet_count);
    for (uint32_t i = 0; i < packet_count; ++i)
    {
        packets[i].BlockId = stream == 1u ? 0u :
            (stream == 2u ? i : input.U32());
        packets[i].Data = packet_data.data() + (size_t)i * block_bytes;
    }

    std::vector<uint8_t> production(7u, 0x31u);
    std::vector<uint8_t> oracle = production;
    const WirehairResult production_result =
        wirehair_v2::SolvePrecodeSystem(
            system, config, packets, block_bytes, production);
    const WirehairResult oracle_result =
        wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
            system, config, packets, block_bytes, oracle);
    if (production_result != oracle_result || production != oracle) {
        return Fail(failure, "sparse solve disagreed with tiny dense oracle");
    }
    return true;
}

bool FuzzTinyMds(wirehair_v2::fuzz::Input& input, std::string& failure)
{
    const uint32_t K = input.Bool() ? 1u : 2u;
    const uint32_t block_bytes = 1u + input.U8() % 32u;
    const uint32_t tail_bytes = 1u + input.U8() % block_bytes;
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(input.U8() ^ (uint8_t)(i * 73u + 19u));
    }

    wirehair_v2::TinyMdsCodec encoder;
    if (encoder.InitializeEncoder(
            message.data(), message_bytes, block_bytes) != Wirehair_Success)
    {
        return Fail(failure, "tiny-MDS fuzz encoder initialization failed");
    }

    uint32_t first_id = input.U32();
    uint32_t second_id = 0u;
    if (K == 2u)
    {
        first_id %= wirehair_v2::TinyMdsCodec::kMaxK2PacketId + 1u;
        second_id = input.U32() %
            wirehair_v2::TinyMdsCodec::kMaxK2PacketId;
        if (second_id >= first_id) {
            ++second_id;
        }
    }
    std::vector<uint8_t> first_packet(block_bytes, 0xa5u);
    uint32_t first_bytes = 0u;
    if (encoder.EncodeResult(
            first_id,
            first_packet.data(),
            block_bytes,
            &first_bytes) != Wirehair_Success ||
        first_bytes == 0u || first_bytes > block_bytes)
    {
        return Fail(failure, "tiny-MDS fuzz first encode failed");
    }

    wirehair_v2::TinyMdsCodec decoder;
    if (decoder.InitializeDecoder(
            message_bytes, block_bytes) != Wirehair_Success)
    {
        return Fail(failure, "tiny-MDS fuzz decoder initialization failed");
    }
    if (K == 2u)
    {
        std::vector<uint8_t> invalid(block_bytes, 0xa5u);
        const std::vector<uint8_t> invalid_before = invalid;
        uint32_t invalid_bytes = UINT32_C(0x12345678);
        const uint32_t invalid_id =
            wirehair_v2::TinyMdsCodec::kMaxK2PacketId +
            1u + input.U8();
        if (encoder.EncodeResult(
                invalid_id,
                invalid.data(),
                block_bytes,
                &invalid_bytes) != Wirehair_InvalidInput ||
            invalid != invalid_before ||
            invalid_bytes != UINT32_C(0x12345678) ||
            decoder.DecodeResult(
                invalid_id,
                invalid.data(),
                block_bytes) != Wirehair_InvalidInput ||
            decoder.ReceivedCount() != 0u)
        {
            return Fail(failure,
                "tiny-MDS fuzz unsupported ID changed state/output");
        }
    }

    const WirehairResult first_result = decoder.DecodeResult(
        first_id, first_packet.data(), first_bytes);
    if (first_result !=
            (K == 1u ? Wirehair_Success : Wirehair_NeedMore) ||
        decoder.ReceivedCount() != 1u)
    {
        return Fail(failure, "tiny-MDS fuzz first decode result mismatch");
    }
    if (K == 2u)
    {
        if (decoder.DecodeResult(
                first_id, first_packet.data(), first_bytes) !=
                Wirehair_NeedMore ||
            decoder.ReceivedCount() != 1u)
        {
            return Fail(failure,
                "tiny-MDS fuzz matching duplicate changed state");
        }
        std::vector<uint8_t> conflict = first_packet;
        conflict[0] ^= 1u;
        if (decoder.DecodeResult(
                first_id, conflict.data(), first_bytes) !=
                Wirehair_InvalidInput ||
            decoder.ReceivedCount() != 1u)
        {
            return Fail(failure,
                "tiny-MDS fuzz conflicting duplicate changed state");
        }

        std::vector<uint8_t> second_packet(block_bytes, 0xa5u);
        uint32_t second_bytes = 0u;
        if (encoder.EncodeResult(
                second_id,
                second_packet.data(),
                block_bytes,
                &second_bytes) != Wirehair_Success ||
            decoder.DecodeResult(
                second_id,
                second_packet.data(),
                second_bytes) != Wirehair_Success ||
            decoder.ReceivedCount() != 2u)
        {
            return Fail(failure, "tiny-MDS fuzz second decode failed");
        }
    }
    if (!decoder.IsDecoded()) {
        return Fail(failure, "tiny-MDS fuzz decoder did not complete");
    }

    std::vector<uint8_t> recovered((size_t)message_bytes + 16u, 0xa5u);
    if (decoder.RecoverResult(
            recovered.data() + 8u, message_bytes) != Wirehair_Success ||
        std::memcmp(
            recovered.data() + 8u,
            message.data(),
            (size_t)message_bytes) != 0 ||
        !std::all_of(
            recovered.begin(), recovered.begin() + 8u,
            [](uint8_t byte) { return byte == 0xa5u; }) ||
        !std::all_of(
            recovered.end() - 8u, recovered.end(),
            [](uint8_t byte) { return byte == 0xa5u; }))
    {
        return Fail(failure, "tiny-MDS fuzz recovery/guard mismatch");
    }

    if (decoder.DecodeResult(
            first_id, first_packet.data(), first_bytes) != Wirehair_Success)
    {
        return Fail(failure,
            "tiny-MDS fuzz consistent completed packet failed");
    }
    first_packet[0] ^= 1u;
    if (decoder.DecodeResult(
            first_id, first_packet.data(), first_bytes) != Wirehair_Error)
    {
        return Fail(failure,
            "tiny-MDS fuzz corrupt completed packet was accepted");
    }
    return true;
}

bool FuzzProfileCase(
    const uint8_t* data,
    size_t size,
    std::string& failure)
{
    if (size > wirehair_v2::fuzz::kMaxFuzzInputBytes) {
        return true;
    }
    wirehair_v2::fuzz::Input input(data, size);
    const uint8_t selector = input.U8();
    // Keep the five historical corpus selectors stable.  One eighth of
    // mutated inputs exercises the independently versioned tiny-MDS state
    // machine without remapping those established seeds.
    if ((selector & 7u) == 7u) {
        return FuzzTinyMds(input, failure);
    }
    switch (selector % 5u)
    {
    case 0: return FuzzParams(input, failure);
    case 1: return FuzzSystem(input, failure);
    case 2: return FuzzProfileContract(input, failure);
    case 3: return FuzzPacketDomain(input, failure);
    default: return FuzzTinyOracle(input, failure);
    }
}

} // namespace

#if defined(WIREHAIR_ENABLE_LIBFUZZER)
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    wirehair_v2::fuzz::RunCoverageGuidedCaseOrAbort(
        WIREHAIR_FUZZ_TARGET_NAME, FuzzProfileCase, data, size);
    return 0;
}
#else
int main(int argc, char** argv)
{
    return wirehair_v2::fuzz::RunDeterministicFuzzer(
        argc,
        argv,
        WIREHAIR_FUZZ_TARGET_NAME,
        WIREHAIR_FUZZ_CORPUS_MANIFEST,
        FuzzProfileCase);
}
#endif
