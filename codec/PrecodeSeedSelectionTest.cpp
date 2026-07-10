#include "WirehairV2PrecodeDecode.h"
#include "WirehairV2PrecodeEncode.h"

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

} // namespace

int main(int argc, char** argv)
{
    const bool sanitizer_mode =
        argc == 2 && std::string(argv[1]) == "--sanitizer";
    if (argc > 2 || (argc == 2 && !sanitizer_mode)) {
        return 2;
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
