#include "V2FuzzDriver.h"

#include "../WirehairV2Codec.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#ifndef WIREHAIR_FUZZ_TARGET_NAME
#define WIREHAIR_FUZZ_TARGET_NAME "stateful"
#endif

#ifndef WIREHAIR_FUZZ_CORPUS_MANIFEST
#define WIREHAIR_FUZZ_CORPUS_MANIFEST "codec/fuzz/corpus/stateful/manifest.txt"
#endif

namespace {

bool Fail(std::string& failure, const char* message)
{
    failure = message;
    return false;
}

struct Fixture
{
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    uint64_t MessageBytes = 0u;
    std::vector<uint8_t> Message;
};

struct Packet
{
    uint32_t Id = 0u;
    uint32_t Bytes = 0u;
    std::vector<uint8_t> Data;
};

Fixture MakeFixture(wirehair_v2::fuzz::Input& input)
{
    Fixture fixture;
    fixture.K = 2u + input.U8() % 7u;
    fixture.BlockBytes = 1u + input.U8() % 32u;
    const uint32_t final_bytes = 1u + input.U8() % fixture.BlockBytes;
    fixture.MessageBytes =
        (uint64_t)(fixture.K - 1u) * fixture.BlockBytes + final_bytes;
    fixture.Message.resize((size_t)fixture.MessageBytes);
    uint64_t state = input.U64() | 1u;
    for (size_t i = 0; i < fixture.Message.size(); ++i)
    {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        fixture.Message[i] = (uint8_t)(
            state * UINT64_C(2685821657736338717));
    }
    return fixture;
}

bool MakeEncoderAndPackets(
    const Fixture& fixture,
    wirehair_v2::MessagePrecodeEncoder& encoder,
    std::vector<Packet>& packets,
    std::string& failure)
{
    if (encoder.InitializeResult(
            fixture.Message.data(), fixture.MessageBytes,
            fixture.BlockBytes) != Wirehair_Success)
    {
        return Fail(failure, "stateful encoder initialization failed");
    }
    packets.resize(fixture.K + 8u);
    for (uint32_t id = 0; id < packets.size(); ++id)
    {
        Packet& packet = packets[id];
        packet.Id = id;
        packet.Data.assign(fixture.BlockBytes, 0u);
        packet.Bytes = UINT32_MAX;
        if (encoder.EncodeResult(
                id,
                packet.Data.data(),
                fixture.BlockBytes,
                &packet.Bytes) != Wirehair_Success ||
            packet.Bytes == 0u || packet.Bytes > fixture.BlockBytes)
        {
            return Fail(failure, "stateful packet generation failed");
        }
    }
    return true;
}

bool InitializeDecoder(
    const Fixture& fixture,
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    wirehair_v2::MessagePrecodeDecoder& decoder,
    std::string& failure)
{
    if (decoder.InitializeResult(
            fixture.MessageBytes,
            fixture.BlockBytes,
            &encoder.Profile()) != Wirehair_Success)
    {
        return Fail(failure, "stateful decoder initialization failed");
    }
    return true;
}

bool RecoverExact(
    const Fixture& fixture,
    const wirehair_v2::MessagePrecodeDecoder& decoder,
    std::string& failure)
{
    std::vector<uint8_t> output(
        (size_t)fixture.MessageBytes + 32u, 0xa5u);
    if (decoder.RecoverResult(
            output.data() + 16u,
            fixture.MessageBytes) != Wirehair_Success ||
        std::memcmp(
            output.data() + 16u,
            fixture.Message.data(),
            (size_t)fixture.MessageBytes) != 0)
    {
        return Fail(failure, "stateful recovered message mismatch");
    }
    for (size_t i = 0; i < 16u; ++i) {
        if (output[i] != 0xa5u ||
            output[16u + fixture.MessageBytes + i] != 0xa5u)
        {
            return Fail(failure, "stateful recovery guard was modified");
        }
    }
    return true;
}

bool FinishSystematic(
    const Fixture& fixture,
    const std::vector<Packet>& packets,
    wirehair_v2::MessagePrecodeDecoder& decoder,
    std::string& failure)
{
    WirehairResult result = decoder.IsDecoded() ?
        Wirehair_Success : Wirehair_NeedMore;
    for (uint32_t id = 0; id < fixture.K; ++id)
    {
        result = decoder.DecodeResult(
            id, packets[id].Data.data(), packets[id].Bytes);
        if (result != Wirehair_NeedMore && result != Wirehair_Success) {
            return Fail(failure, "systematic completion returned error");
        }
    }
    if (!decoder.IsDecoded() || result != Wirehair_Success) {
        return Fail(failure, "systematic completion did not decode");
    }
    return RecoverExact(fixture, decoder, failure);
}

bool GeneralSequence(
    wirehair_v2::fuzz::Input& input,
    const Fixture& fixture,
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    const std::vector<Packet>& packets,
    std::string& failure)
{
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!InitializeDecoder(fixture, encoder, decoder, failure)) {
        return false;
    }
    std::vector<uint8_t> delivered(packets.size(), 0u);
    const unsigned steps = 1u + input.U8() % 32u;
    for (unsigned step = 0; step < steps; ++step)
    {
        const unsigned operation = input.U8() % 5u;
        const uint32_t id = input.U8() % (uint32_t)packets.size();
        const Packet& packet = packets[id];
        const uint32_t before_count = decoder.ReceivedCount();
        if (operation == 0u)
        {
            const WirehairResult result = decoder.DecodeResult(
                id, packet.Data.data(), packet.Bytes);
            if (result != Wirehair_NeedMore && result != Wirehair_Success) {
                return Fail(failure, "valid reordered packet returned error");
            }
            delivered[id] = 1u;
        }
        else if (operation == 1u)
        {
            const uint32_t wrong_bytes = packet.Bytes == fixture.BlockBytes ?
                packet.Bytes - 1u : packet.Bytes + 1u;
            if (decoder.DecodeResult(
                    id, packet.Data.data(), wrong_bytes) !=
                    Wirehair_InvalidInput ||
                decoder.ReceivedCount() != before_count)
            {
                return Fail(failure, "wrong packet size was not no-state");
            }
        }
        else if (operation == 2u)
        {
            WirehairResult first = Wirehair_NeedMore;
            if (!delivered[id]) {
                first = decoder.DecodeResult(
                    id, packet.Data.data(), packet.Bytes);
                if (first != Wirehair_NeedMore && first != Wirehair_Success) {
                    return Fail(failure, "duplicate setup packet failed");
                }
                delivered[id] = 1u;
            }
            const uint32_t committed_count = decoder.ReceivedCount();
            const WirehairResult duplicate = decoder.DecodeResult(
                id, packet.Data.data(), packet.Bytes);
            if (duplicate != Wirehair_NeedMore &&
                duplicate != Wirehair_Success) {
                return Fail(failure, "identical duplicate returned error");
            }
            if (decoder.ReceivedCount() != committed_count) {
                return Fail(failure, "identical duplicate changed count");
            }
        }
        else if (operation == 3u)
        {
            if (!delivered[id]) {
                const WirehairResult first = decoder.DecodeResult(
                    id, packet.Data.data(), packet.Bytes);
                if (first != Wirehair_NeedMore && first != Wirehair_Success) {
                    return Fail(failure, "conflict setup packet failed");
                }
                delivered[id] = 1u;
            }
            std::vector<uint8_t> corrupt = packet.Data;
            corrupt[0] ^= 1u;
            const uint32_t committed_count = decoder.ReceivedCount();
            const WirehairResult conflict = decoder.DecodeResult(
                id, corrupt.data(), packet.Bytes);
            const WirehairResult expected = decoder.IsDecoded() ?
                Wirehair_Error : Wirehair_InvalidInput;
            if (conflict != expected ||
                decoder.ReceivedCount() != committed_count)
            {
                return Fail(failure, "conflicting duplicate contract failed");
            }
        }
        else
        {
            std::vector<uint8_t> output((size_t)fixture.MessageBytes, 0x5au);
            const std::vector<uint8_t> before = output;
            const WirehairResult recovered = decoder.RecoverResult(
                output.data(), fixture.MessageBytes);
            if (decoder.IsDecoded()) {
                if (recovered != Wirehair_Success ||
                    output != fixture.Message) {
                    return Fail(failure, "completed random recovery failed");
                }
            }
            else if (recovered != Wirehair_NeedMore || output != before) {
                return Fail(failure, "early recovery was not no-write");
            }
        }
    }
    return FinishSystematic(fixture, packets, decoder, failure);
}

bool DuplicateAndConflict(
    const Fixture& fixture,
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    const std::vector<Packet>& packets,
    std::string& failure)
{
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!InitializeDecoder(fixture, encoder, decoder, failure)) return false;
    const Packet& packet = packets[0];
    if (decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes) != Wirehair_NeedMore)
    {
        return Fail(failure, "rank-deficient prefix did not need more");
    }
    for (uint32_t i = 0; i < fixture.K + 3u; ++i) {
        if (decoder.DecodeResult(
                packet.Id, packet.Data.data(), packet.Bytes) !=
                Wirehair_NeedMore ||
            decoder.ReceivedCount() != 1u)
        {
            return Fail(failure, "rank-deficient duplicate was not cached");
        }
    }
    std::vector<uint8_t> conflict = packet.Data;
    conflict[0] ^= 0x80u;
    if (decoder.DecodeResult(
            packet.Id, conflict.data(), packet.Bytes) !=
            Wirehair_InvalidInput ||
        decoder.ReceivedCount() != 1u)
    {
        return Fail(failure, "pre-completion conflict was not rejected");
    }
    return FinishSystematic(fixture, packets, decoder, failure);
}

class DecoderAllocationReset
{
public:
    ~DecoderAllocationReset()
    {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
#endif
    }
};

bool OomRetry(
    const Fixture& fixture,
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    const std::vector<Packet>& packets,
    std::string& failure)
{
#if !defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    (void)fixture; (void)encoder; (void)packets;
    return Fail(failure, "OOM fuzz scenario requires test hooks");
#else
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!InitializeDecoder(fixture, encoder, decoder, failure)) return false;
    for (uint32_t id = 0; id + 1u < fixture.K; ++id) {
        if (decoder.DecodeResult(
                id, packets[id].Data.data(), packets[id].Bytes) !=
                Wirehair_NeedMore)
        {
            return Fail(failure, "OOM prefix did not need more");
        }
    }
    DecoderAllocationReset reset;
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
    const Packet& final_packet = packets[fixture.K - 1u];
    if (decoder.DecodeResult(
            final_packet.Id,
            final_packet.Data.data(),
            final_packet.Bytes) != Wirehair_OOM ||
        decoder.IsDecoded() || decoder.ReceivedCount() != fixture.K)
    {
        return Fail(failure, "injected solve OOM state mismatch");
    }
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    if (decoder.DecodeResult(
            final_packet.Id,
            final_packet.Data.data(),
            final_packet.Bytes) != Wirehair_Success ||
        !decoder.IsDecoded())
    {
        return Fail(failure, "identical packet did not retry OOM solve");
    }
    return RecoverExact(fixture, decoder, failure);
#endif
}

bool PostCompletion(
    const Fixture& fixture,
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    const std::vector<Packet>& packets,
    std::string& failure)
{
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!InitializeDecoder(fixture, encoder, decoder, failure) ||
        !FinishSystematic(fixture, packets, decoder, failure))
    {
        return false;
    }
    const uint32_t count = decoder.ReceivedCount();
    for (const Packet& packet : packets)
    {
        if (decoder.DecodeResult(
                packet.Id, packet.Data.data(), packet.Bytes) !=
                Wirehair_Success)
        {
            return Fail(failure, "consistent post-completion packet failed");
        }
        std::vector<uint8_t> corrupt = packet.Data;
        corrupt[0] ^= 1u;
        if (decoder.DecodeResult(
                packet.Id, corrupt.data(), packet.Bytes) != Wirehair_Error ||
            decoder.ReceivedCount() != count)
        {
            return Fail(failure, "corrupt post-completion packet was accepted");
        }
    }
    return RecoverExact(fixture, decoder, failure);
}

class CodecAllocationReset
{
public:
    ~CodecAllocationReset()
    {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::SetCodecAllocationFailureCountdownForTesting(-1);
#endif
    }
};

bool FacadeReuse(
    const Fixture& fixture,
    const std::vector<Packet>& packets,
    std::string& failure)
{
    wirehair_v2::Codec codec;
    if (codec.InitializePrecodeEncoder(
            fixture.Message.data(), fixture.MessageBytes,
            fixture.BlockBytes) != Wirehair_Success)
    {
        return Fail(failure, "facade precode encoder init failed");
    }
    std::vector<uint8_t> expected(fixture.BlockBytes, 0u);
    uint32_t expected_bytes = 0u;
    if (codec.Encode(
            fixture.K + 1u, expected.data(), fixture.BlockBytes,
            &expected_bytes) != Wirehair_Success)
    {
        return Fail(failure, "facade precode encode failed");
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    CodecAllocationReset reset;
    wirehair_v2::SetCodecAllocationFailureCountdownForTesting(0);
    if (codec.InitializeDecoder(
            fixture.MessageBytes, fixture.BlockBytes) != Wirehair_OOM)
    {
        return Fail(failure, "facade injected transition did not OOM");
    }
    wirehair_v2::SetCodecAllocationFailureCountdownForTesting(-1);
#endif
    std::vector<uint8_t> after(fixture.BlockBytes, 0u);
    uint32_t after_bytes = 0u;
    if (codec.Encode(
            fixture.K + 1u, after.data(), fixture.BlockBytes,
            &after_bytes) != Wirehair_Success ||
        after_bytes != expected_bytes || after != expected)
    {
        return Fail(failure, "failed facade transition changed old encoder");
    }
    if (codec.InitializePrecodeDecoder(
            fixture.MessageBytes, fixture.BlockBytes,
            &codec.Profile()) != Wirehair_Success)
    {
        return Fail(failure, "facade encoder-to-decoder transition failed");
    }
    WirehairResult decode_result = Wirehair_NeedMore;
    for (uint32_t id = 0; id < fixture.K; ++id) {
        decode_result = codec.Decode(
            id, packets[id].Data.data(), packets[id].Bytes);
    }
    std::vector<uint8_t> recovered((size_t)fixture.MessageBytes, 0u);
    if (decode_result != Wirehair_Success ||
        codec.Recover(recovered.data(), fixture.MessageBytes) !=
            Wirehair_Success ||
        recovered != fixture.Message)
    {
        return Fail(failure, "facade precode decoder recovery failed");
    }
    if (codec.InitializeEncoder(
            fixture.Message.data(), fixture.MessageBytes,
            fixture.BlockBytes) != Wirehair_Success)
    {
        return Fail(failure, "facade decoder-to-V1 encoder transition failed");
    }
    std::vector<uint8_t> legacy(fixture.BlockBytes, 0u);
    uint32_t legacy_bytes = 0u;
    if (codec.Encode(
            fixture.K + 3u, legacy.data(), fixture.BlockBytes,
            &legacy_bytes) != Wirehair_Success || legacy_bytes == 0u)
    {
        return Fail(failure, "facade V1 encode after reuse failed");
    }
    return true;
}

bool MalformedProfileRetry(
    wirehair_v2::fuzz::Input& input,
    const Fixture& fixture,
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    const std::vector<Packet>& packets,
    std::string& failure)
{
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!InitializeDecoder(fixture, encoder, decoder, failure)) return false;
    wirehair_v2::SeedProfile bad = encoder.Profile();
    switch (input.U8() % 6u)
    {
    case 0: ++bad.V2PrecodeContractVersion; break;
    case 1: ++bad.V2PacketRowContractVersion; break;
    case 2: bad.V2PrecodeSeedSalt ^= 1u; break;
    case 3: bad.V2RecoveryRowSeedSalt ^= 1u; break;
    case 4: bad.V2RecoveryMixCount = 0u; break;
    default: ++bad.V2StaircaseCount; break;
    }
    if (decoder.InitializeResult(
            fixture.MessageBytes, fixture.BlockBytes, &bad) !=
            Wirehair_InvalidInput ||
        !decoder.IsInitialized())
    {
        return Fail(failure, "malformed profile retry was not transactional");
    }
    return FinishSystematic(fixture, packets, decoder, failure);
}

bool FuzzStatefulCase(
    const uint8_t* data,
    size_t size,
    std::string& failure)
{
    if (size > wirehair_v2::fuzz::kMaxFuzzInputBytes) return true;
    wirehair_v2::fuzz::Input input(data, size);
    const unsigned scenario = input.U8() % 6u;
    const Fixture fixture = MakeFixture(input);
    wirehair_v2::MessagePrecodeEncoder encoder;
    std::vector<Packet> packets;
    if (!MakeEncoderAndPackets(
            fixture, encoder, packets, failure))
    {
        return false;
    }
    switch (scenario)
    {
    case 0:
        return GeneralSequence(
            input, fixture, encoder, packets, failure);
    case 1:
        return DuplicateAndConflict(
            fixture, encoder, packets, failure);
    case 2:
        return OomRetry(fixture, encoder, packets, failure);
    case 3:
        return PostCompletion(fixture, encoder, packets, failure);
    case 4:
        return FacadeReuse(fixture, packets, failure);
    default:
        return MalformedProfileRetry(
            input, fixture, encoder, packets, failure);
    }
}

} // namespace

#if defined(WIREHAIR_ENABLE_LIBFUZZER)
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    wirehair_v2::fuzz::RunCoverageGuidedCaseOrAbort(
        WIREHAIR_FUZZ_TARGET_NAME, FuzzStatefulCase, data, size);
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
        FuzzStatefulCase);
}
#endif
