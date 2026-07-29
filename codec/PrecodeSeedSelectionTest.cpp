#include "WirehairV2PrecodeDecode.h"
#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Contract.h"

#include <cstdio>
#include <string>
#include <vector>

namespace {

uint32_t ExpectedAttempt(uint32_t K)
{
    struct Golden { uint32_t K; uint32_t Attempt; };
    const Golden goldens[] = {
        {2u, 2u}, {4u, 1u}, {7u, 1u}, {14u, 2u},
        {120u, 1u}, {196u, 1u}, {336u, 1u}, {598u, 1u},
        {899u, 1u}, {961u, 1u}, {1056u, 1u}, {1133u, 1u},
        {1217u, 1u}
    };
    for (const Golden& golden : goldens) {
        if (K == golden.K) {
            return golden.Attempt;
        }
    }
    return 0u;
}

bool CheckK(uint32_t K, uint32_t expected_attempt, uint32_t& attempt_out)
{
    const uint32_t block_bytes = 1u;
    std::vector<uint8_t> message(K);
    for (uint32_t i = 0; i < K; ++i) {
        message[i] = (uint8_t)(i * 29u + K * 7u + 3u);
    }

    wirehair_v2::MessagePrecodeEncoder encoder;
    const WirehairResult encode_result = encoder.InitializeResult(
        message.data(), message.size(), block_bytes);
    if (encode_result != Wirehair_Success)
    {
        std::fprintf(stderr,
            "seed selection: encoder K=%u failed result=%d\n",
            K, (int)encode_result);
        return false;
    }
    if (encoder.Profile().V2PacketRowContractVersion != 4u ||
        wirehair_v2::kPacketRowContractVersion != 4u)
    {
        std::fprintf(stderr,
            "seed selection: K=%u did not publish packet contract v4\n", K);
        return false;
    }
    const uint32_t attempt = encoder.SolveStats().PacketSeedAttempt;
    attempt_out = attempt;
    if (attempt != expected_attempt)
    {
        std::fprintf(stderr,
            "seed selection: K=%u attempt=%u expected=%u\n",
            K, attempt, expected_attempt);
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder decoder;
    const WirehairResult decode_result = decoder.InitializeResult(
        message.size(), block_bytes, &encoder.Profile());
    if (decode_result != Wirehair_Success ||
        decoder.PacketSeedAttempt() != attempt ||
        decoder.PacketPeelSeed() !=
            (uint32_t)encoder.BlockEncoder().RecoveryRowSeed() ||
        decoder.System().Params.Seed !=
            encoder.BlockEncoder().System().Params.Seed)
    {
        std::fprintf(stderr,
            "seed selection: encoder/decoder contract mismatch K=%u "
            "enc_attempt=%u dec_attempt=%u result=%d\n",
            K, attempt, decoder.PacketSeedAttempt(), (int)decode_result);
        return false;
    }
    if (expected_attempt != 0u)
    {
        wirehair_v2::MessagePrecodeDecoder standalone;
        if (standalone.InitializeResult(message.size(), block_bytes) !=
                Wirehair_Success ||
            standalone.PacketSeedAttempt() != expected_attempt ||
            standalone.PacketPeelSeed() != decoder.PacketPeelSeed() ||
            standalone.System().Params.Seed != decoder.System().Params.Seed)
        {
            std::fprintf(stderr,
                "seed selection: standalone decoder mismatch K=%u\n", K);
            return false;
        }
    }

    uint8_t block = 0u;
    WirehairResult feed_result = Wirehair_NeedMore;
    for (uint32_t block_id = 0; block_id < K; ++block_id)
    {
        uint32_t bytes = 0u;
        if (!encoder.Encode(
                block_id, &block, 1u, &bytes) || bytes != 1u)
        {
            std::fprintf(stderr,
                "seed selection: source encode failed K=%u id=%u\n",
                K, block_id);
            return false;
        }
        feed_result = decoder.DecodeResult(block_id, &block, bytes);
        const WirehairResult expected = block_id + 1u == K ?
            Wirehair_Success : Wirehair_NeedMore;
        if (feed_result != expected)
        {
            std::fprintf(stderr,
                "seed selection: source decode K=%u id=%u result=%d "
                "expected=%d\n",
                K, block_id, (int)feed_result, (int)expected);
            return false;
        }
    }
    std::vector<uint8_t> recovered(K, 0u);
    if (decoder.RecoverResult(recovered.data(), K) != Wirehair_Success ||
        recovered != message)
    {
        std::fprintf(stderr,
            "seed selection: roundtrip mismatch K=%u\n", K);
        return false;
    }
    return true;
}

bool CheckRawPinnedWeakClassification()
{
    static const uint32_t K = 10u;
    static const uint64_t construction_seed =
        UINT64_C(0x78dde6e660de777f);
    const wirehair_v2::V2EquationContract* contract =
        wirehair_v2::FindV2EquationContract("dispatch-v1");
    if (!contract) {
        std::fprintf(stderr,
            "seed selection: could not find raw target contract\n");
        return false;
    }
    const wirehair_v2::MessagePrecodeEncoderOptions options =
        wirehair_v2::MessageOptionsForContract(*contract);
    const uint32_t widths[] = {2u, 1280u, 4096u};
    for (uint32_t block_bytes : widths)
    {
        wirehair_v2::SeedProfile profile;
        if (!wirehair_v2::MakeRawContractProfile(
                *contract, K, block_bytes, construction_seed, profile))
        {
            std::fprintf(stderr,
                "seed selection: could not build raw weak-seed profile "
                "bb=%u\n", block_bytes);
            return false;
        }
        std::vector<uint8_t> zero_message(
            (size_t)K * block_bytes, uint8_t{0});
        std::vector<uint8_t> nonzero_message(zero_message.size());
        for (size_t i = 0u; i < nonzero_message.size(); ++i) {
            nonzero_message[i] = (uint8_t)(i * 29u + 7u);
        }

        wirehair_v2::MessagePrecodeEncoder zero_encoder;
        wirehair_v2::MessagePrecodeEncoder nonzero_encoder;
        const WirehairResult zero_result = zero_encoder.InitializeResult(
            zero_message.data(), zero_message.size(), block_bytes,
            &profile, &options);
        const WirehairResult nonzero_result =
            nonzero_encoder.InitializeResult(
                nonzero_message.data(), nonzero_message.size(), block_bytes,
                &profile, &options);
        if (zero_result != Wirehair_BadPeelSeed ||
            nonzero_result != Wirehair_BadPeelSeed)
        {
            std::fprintf(stderr,
                "seed selection: raw singular construction changed class "
                "with RHS bb=%u zero=%d nonzero=%d\n",
                block_bytes, (int)zero_result, (int)nonzero_result);
            return false;
        }

        wirehair_v2::ExplicitMessagePrecodeConfigForTesting explicit_config;
        if (!wirehair_v2::ResolveMessagePrecodeConfiguration(
                profile,
                options,
                explicit_config.Params,
                explicit_config.Packet) ||
            !explicit_config.PinActiveEquationStateForTesting())
        {
            std::fprintf(stderr,
                "seed selection: weak explicit configuration resolution "
                "failed bb=%u\n", block_bytes);
            return false;
        }
        wirehair_v2::MessagePrecodeEncoder explicit_zero;
        wirehair_v2::MessagePrecodeEncoder explicit_nonzero;
        const WirehairResult explicit_zero_result =
            explicit_zero.InitializeExplicitResultForTesting(
                zero_message.data(), zero_message.size(),
                block_bytes, explicit_config);
        const WirehairResult explicit_nonzero_result =
            explicit_nonzero.InitializeExplicitResultForTesting(
                nonzero_message.data(), nonzero_message.size(),
                block_bytes, explicit_config);
        if (explicit_zero_result != Wirehair_NeedMore ||
            (explicit_nonzero_result != Wirehair_NeedMore &&
             explicit_nonzero_result != Wirehair_Error) ||
            (block_bytes == 2u &&
             explicit_nonzero_result != Wirehair_Error) ||
            explicit_zero.IsInitialized() ||
            explicit_nonzero.IsInitialized())
        {
            std::fprintf(stderr,
                "seed selection: explicit attempt-zero result was not raw "
                "bb=%u zero=%d nonzero=%d\n",
                block_bytes,
                (int)explicit_zero_result,
                (int)explicit_nonzero_result);
            return false;
        }
        if (block_bytes == 2u)
        {
            // The exact-width diagnostic zero block is the fourth guarded
            // allocation on this complete, uncached input.  Classification
            // must not turn failure of that diagnostic into BadPeelSeed.
            wirehair_v2::MessagePrecodeEncoder oom_encoder;
            wirehair_v2::SetAllocationFailureCountdownForTesting(3);
            const WirehairResult oom_result = oom_encoder.InitializeResult(
                nonzero_message.data(), nonzero_message.size(), block_bytes,
                &profile, &options);
            wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
            if (oom_result != Wirehair_OOM ||
                oom_encoder.IsInitialized())
            {
                std::fprintf(stderr,
                    "seed selection: weak-seed diagnostic OOM was not "
                    "transactional result=%d\n", (int)oom_result);
                return false;
            }

            // The explicit path has no fourth diagnostic allocation or
            // zero-RHS rerun.  Countdown three therefore reaches and returns
            // the original Error unchanged.
            wirehair_v2::MessagePrecodeEncoder explicit_no_probe;
            wirehair_v2::SetAllocationFailureCountdownForTesting(3);
            const WirehairResult explicit_no_probe_result =
                explicit_no_probe.InitializeExplicitResultForTesting(
                    nonzero_message.data(), nonzero_message.size(),
                    block_bytes, explicit_config);
            wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
            if (explicit_no_probe_result != Wirehair_Error ||
                explicit_no_probe.IsInitialized())
            {
                std::fprintf(stderr,
                    "seed selection: explicit weak path performed a hidden "
                    "probe result=%d\n", (int)explicit_no_probe_result);
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    const bool sanitizer_mode =
        argc == 2 && std::string(argv[1]) == "--sanitizer";
    if (argc > 2 || (argc == 2 && !sanitizer_mode)) {
        return 2;
    }
    if (!CheckRawPinnedWeakClassification()) {
        return 1;
    }
    uint32_t fixups = 0u;
    uint32_t max_attempt = 0u;
    const uint32_t sanitizer_cases[] = {
        2u, 3u, 4u, 7u, 14u, 64u, 120u, 196u, 320u, 336u,
        511u, 598u, 899u, 961u, 1000u, 1056u, 1133u, 1217u,
        2048u, 3200u, 10000u
    };
    const uint32_t first = sanitizer_mode ? 0u : 2u;
    const uint32_t last = sanitizer_mode ?
        (uint32_t)(sizeof(sanitizer_cases) / sizeof(sanitizer_cases[0])) - 1u :
        2048u;
    for (uint32_t index = first; index <= last; ++index)
    {
        const uint32_t K = sanitizer_mode ?
            sanitizer_cases[index] : index;
        uint32_t attempt = 0u;
        if (!CheckK(K, ExpectedAttempt(K), attempt)) {
            return 1;
        }
        fixups += attempt != 0u ? 1u : 0u;
        if (attempt > max_attempt) {
            max_attempt = attempt;
        }
    }

    const uint32_t representatives[] = {3200u, 10000u, 32000u, 64000u};
    if (!sanitizer_mode)
    {
        for (uint32_t K : representatives)
        {
            uint32_t attempt = 0u;
            if (!CheckK(K, 0u, attempt)) {
                return 1;
            }
        }
    }
    std::printf(
        "default systematic seed selection %s: PASS "
        "(fixups=%u max_attempt=%u)\n",
        sanitizer_mode ? "sanitizer golden set" : "K=2..2048 + large reps",
        fixups, max_attempt);
    return 0;
}
