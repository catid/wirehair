#include "WirehairV2TinyMDS.h"
#include "WirehairV2Codec.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

namespace {

bool Check(bool condition, const char* what)
{
    if (condition) {
        return true;
    }
    std::fprintf(stderr, "tiny MDS test failed: %s\n", what);
    return false;
}

void FillMessage(std::vector<uint8_t>& message)
{
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(19u + i * 73u);
    }
}

struct PacketPair
{
    uint32_t FirstId = 0u;
    uint32_t SecondId = 0u;
    std::vector<uint8_t> First;
    std::vector<uint8_t> Second;
};

struct SpeedFixture
{
    uint32_t BlockBytes = 0u;
    uint64_t MessageBytes = 0u;
    std::vector<uint8_t> Message;
    PacketPair TinyPackets;
    PacketPair V2Packets;
    PacketPair LegacyPackets;
    wirehair_v2::SeedProfile V2Profile = {};
};

bool EncodeV2Packet(
    wirehair_v2::Codec& encoder,
    uint32_t id,
    uint32_t block_bytes,
    std::vector<uint8_t>& packet)
{
    packet.assign(block_bytes, uint8_t{0});
    uint32_t packet_bytes = 0u;
    return encoder.Encode(
            id, packet.data(), block_bytes, &packet_bytes) ==
            Wirehair_Success &&
        packet_bytes == block_bytes;
}

bool FindV2RepairPair(
    wirehair_v2::Codec& encoder,
    SpeedFixture& fixture)
{
    // A fountain construction does not promise that any particular two
    // repair rows are independent.  Find and freeze one successful repair
    // pair outside the timed region rather than silently timing a failed
    // decoder or giving the baseline systematic packets.
    for (uint32_t first = 2u; first <= 64u; ++first)
    {
        std::vector<uint8_t> first_packet;
        if (!EncodeV2Packet(
                encoder, first, fixture.BlockBytes, first_packet))
        {
            return false;
        }
        for (uint32_t second = first + 1u; second <= 256u; ++second)
        {
            std::vector<uint8_t> second_packet;
            if (!EncodeV2Packet(
                    encoder, second, fixture.BlockBytes, second_packet))
            {
                return false;
            }
            wirehair_v2::Codec decoder;
            if (decoder.InitializePrecodeDecoder(
                    fixture.MessageBytes,
                    fixture.BlockBytes,
                    &fixture.V2Profile) != Wirehair_Success ||
                decoder.Decode(
                    first, first_packet.data(), fixture.BlockBytes) !=
                    Wirehair_NeedMore ||
                decoder.Decode(
                    second, second_packet.data(), fixture.BlockBytes) !=
                    Wirehair_Success)
            {
                continue;
            }
            std::vector<uint8_t> recovered(fixture.Message.size(), 0u);
            if (decoder.Recover(
                    recovered.data(), recovered.size()) != Wirehair_Success ||
                recovered != fixture.Message)
            {
                return false;
            }
            fixture.V2Packets.FirstId = first;
            fixture.V2Packets.SecondId = second;
            fixture.V2Packets.First.swap(first_packet);
            fixture.V2Packets.Second.swap(second_packet);
            return true;
        }
    }
    return false;
}

bool EncodeLegacyPacket(
    wirehair::Codec& encoder,
    uint32_t id,
    uint32_t block_bytes,
    std::vector<uint8_t>& packet)
{
    packet.assign(block_bytes, uint8_t{0});
    return encoder.Encode(id, packet.data(), block_bytes) == block_bytes;
}

bool FindLegacyRepairPair(
    wirehair::Codec& encoder,
    SpeedFixture& fixture)
{
    for (uint32_t first = 2u; first <= 64u; ++first)
    {
        std::vector<uint8_t> first_packet;
        if (!EncodeLegacyPacket(
                encoder, first, fixture.BlockBytes, first_packet))
        {
            return false;
        }
        for (uint32_t second = first + 1u; second <= 256u; ++second)
        {
            std::vector<uint8_t> second_packet;
            if (!EncodeLegacyPacket(
                    encoder, second, fixture.BlockBytes, second_packet))
            {
                return false;
            }
            wirehair::Codec decoder;
            if (decoder.InitializeDecoder(
                    fixture.MessageBytes,
                    fixture.BlockBytes) != Wirehair_Success ||
                decoder.DecodeFeed(
                    first, first_packet.data(), fixture.BlockBytes) !=
                    Wirehair_NeedMore ||
                decoder.DecodeFeed(
                    second, second_packet.data(), fixture.BlockBytes) !=
                    Wirehair_Success)
            {
                continue;
            }
            std::vector<uint8_t> recovered(fixture.Message.size(), 0u);
            if (decoder.ReconstructOutput(
                    recovered.data(), recovered.size()) != Wirehair_Success ||
                recovered != fixture.Message)
            {
                return false;
            }
            fixture.LegacyPackets.FirstId = first;
            fixture.LegacyPackets.SecondId = second;
            fixture.LegacyPackets.First.swap(first_packet);
            fixture.LegacyPackets.Second.swap(second_packet);
            return true;
        }
    }
    return false;
}

bool PrepareSpeedFixture(uint32_t block_bytes, SpeedFixture& fixture)
{
    fixture = SpeedFixture();
    fixture.BlockBytes = block_bytes;
    fixture.MessageBytes = (uint64_t)block_bytes * 2u;
    fixture.Message.resize((size_t)fixture.MessageBytes);
    FillMessage(fixture.Message);

    wirehair_v2::TinyMdsCodec tiny_encoder;
    fixture.TinyPackets.FirstId = 2u;
    fixture.TinyPackets.SecondId = 3u;
    fixture.TinyPackets.First.assign(block_bytes, uint8_t{0});
    fixture.TinyPackets.Second.assign(block_bytes, uint8_t{0});
    uint32_t first_bytes = 0u;
    uint32_t second_bytes = 0u;
    if (tiny_encoder.InitializeEncoder(
            fixture.Message.data(),
            fixture.MessageBytes,
            block_bytes) != Wirehair_Success ||
        tiny_encoder.EncodeResult(
            fixture.TinyPackets.FirstId,
            fixture.TinyPackets.First.data(),
            block_bytes,
            &first_bytes) != Wirehair_Success ||
        tiny_encoder.EncodeResult(
            fixture.TinyPackets.SecondId,
            fixture.TinyPackets.Second.data(),
            block_bytes,
            &second_bytes) != Wirehair_Success ||
        first_bytes != block_bytes || second_bytes != block_bytes)
    {
        return false;
    }

    wirehair_v2::TinyMdsCodec tiny_decoder;
    std::vector<uint8_t> recovered(fixture.Message.size(), 0u);
    if (tiny_decoder.InitializeDecoder(
            fixture.MessageBytes, block_bytes) != Wirehair_Success ||
        tiny_decoder.DecodeResult(
            fixture.TinyPackets.FirstId,
            fixture.TinyPackets.First.data(),
            block_bytes) != Wirehair_NeedMore ||
        tiny_decoder.DecodeResult(
            fixture.TinyPackets.SecondId,
            fixture.TinyPackets.Second.data(),
            block_bytes) != Wirehair_Success ||
        tiny_decoder.RecoverResult(
            recovered.data(), recovered.size()) != Wirehair_Success ||
        recovered != fixture.Message)
    {
        return false;
    }

    wirehair_v2::Codec v2_encoder;
    if (v2_encoder.InitializePrecodeEncoder(
            fixture.Message.data(),
            fixture.MessageBytes,
            block_bytes) != Wirehair_Success)
    {
        return false;
    }
    fixture.V2Profile = v2_encoder.Profile();
    if (!FindV2RepairPair(v2_encoder, fixture)) {
        return false;
    }

    wirehair::Codec legacy_encoder;
    if (legacy_encoder.InitializeEncoder(
            fixture.MessageBytes, block_bytes) != Wirehair_Success ||
        legacy_encoder.EncodeFeed(fixture.Message.data()) != Wirehair_Success ||
        !FindLegacyRepairPair(legacy_encoder, fixture))
    {
        return false;
    }
    return true;
}

volatile uint64_t SpeedBenchmarkSink = 0u;

bool MeasureTinyBatch(
    const SpeedFixture& fixture,
    uint32_t iterations,
    double& ns_per_operation)
{
    uint64_t completions = 0u;
    uint64_t failures = 0u;
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (uint32_t i = 0u; i < iterations; ++i)
    {
        wirehair_v2::TinyMdsCodec decoder;
        WirehairResult result = decoder.InitializeDecoder(
            fixture.MessageBytes, fixture.BlockBytes);
        if (result == Wirehair_Success) {
            result = decoder.DecodeResult(
                fixture.TinyPackets.FirstId,
                fixture.TinyPackets.First.data(),
                fixture.BlockBytes);
        }
        if (result == Wirehair_NeedMore) {
            result = decoder.DecodeResult(
                fixture.TinyPackets.SecondId,
                fixture.TinyPackets.Second.data(),
                fixture.BlockBytes);
        }
        if (result == Wirehair_Success) {
            ++completions;
        }
        else {
            failures += (uint64_t)(uint32_t)result + 1u;
        }
    }
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    SpeedBenchmarkSink ^= completions + failures;
    if (completions != iterations || failures != 0u) {
        return false;
    }
    const uint64_t elapsed_ns = (uint64_t)
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            end - begin).count();
    if (elapsed_ns == 0u) {
        return false;
    }
    ns_per_operation = (double)elapsed_ns / iterations;
    return true;
}

bool MeasureV2Batch(
    const SpeedFixture& fixture,
    uint32_t iterations,
    double& ns_per_operation)
{
    uint64_t completions = 0u;
    uint64_t failures = 0u;
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (uint32_t i = 0u; i < iterations; ++i)
    {
        wirehair_v2::Codec decoder;
        WirehairResult result = decoder.InitializePrecodeDecoder(
            fixture.MessageBytes,
            fixture.BlockBytes,
            &fixture.V2Profile);
        if (result == Wirehair_Success) {
            result = decoder.Decode(
                fixture.V2Packets.FirstId,
                fixture.V2Packets.First.data(),
                fixture.BlockBytes);
        }
        if (result == Wirehair_NeedMore) {
            result = decoder.Decode(
                fixture.V2Packets.SecondId,
                fixture.V2Packets.Second.data(),
                fixture.BlockBytes);
        }
        if (result == Wirehair_Success) {
            ++completions;
        }
        else {
            failures += (uint64_t)(uint32_t)result + 1u;
        }
    }
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    SpeedBenchmarkSink ^= completions + failures;
    if (completions != iterations || failures != 0u) {
        return false;
    }
    const uint64_t elapsed_ns = (uint64_t)
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            end - begin).count();
    if (elapsed_ns == 0u) {
        return false;
    }
    ns_per_operation = (double)elapsed_ns / iterations;
    return true;
}

bool MeasureLegacyBatch(
    const SpeedFixture& fixture,
    uint32_t iterations,
    double& ns_per_operation)
{
    uint64_t completions = 0u;
    uint64_t failures = 0u;
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (uint32_t i = 0u; i < iterations; ++i)
    {
        wirehair::Codec decoder;
        WirehairResult result = decoder.InitializeDecoder(
            fixture.MessageBytes, fixture.BlockBytes);
        if (result == Wirehair_Success) {
            result = decoder.DecodeFeed(
                fixture.LegacyPackets.FirstId,
                fixture.LegacyPackets.First.data(),
                fixture.BlockBytes);
        }
        if (result == Wirehair_NeedMore) {
            result = decoder.DecodeFeed(
                fixture.LegacyPackets.SecondId,
                fixture.LegacyPackets.Second.data(),
                fixture.BlockBytes);
        }
        if (result == Wirehair_Success) {
            ++completions;
        }
        else {
            failures += (uint64_t)(uint32_t)result + 1u;
        }
    }
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    SpeedBenchmarkSink ^= completions + failures;
    if (completions != iterations || failures != 0u) {
        return false;
    }
    const uint64_t elapsed_ns = (uint64_t)
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            end - begin).count();
    if (elapsed_ns == 0u) {
        return false;
    }
    ns_per_operation = (double)elapsed_ns / iterations;
    return true;
}

enum class BaselineKind : uint8_t
{
    V2Precode,
    Legacy
};

bool MeasureBaselineBatch(
    BaselineKind baseline,
    const SpeedFixture& fixture,
    uint32_t iterations,
    double& ns_per_operation)
{
    return baseline == BaselineKind::V2Precode ?
        MeasureV2Batch(fixture, iterations, ns_per_operation) :
        MeasureLegacyBatch(fixture, iterations, ns_per_operation);
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2u;
    if ((values.size() & 1u) != 0u) {
        return values[middle];
    }
    return (values[middle - 1u] + values[middle]) * 0.5;
}

bool RunPairedSpeedComparison(
    const SpeedFixture& fixture,
    BaselineKind baseline)
{
    static const char kOrder[] = {
        'A', 'B', 'B', 'A', 'B', 'A', 'A', 'B'
    };
    static const uint32_t kCycles = 7u;
    static const uint32_t kIterations = 4096u;
    static const double kMedianRatioLimit = 0.50;
    static const double kWorstRatioLimit = 0.80;

    std::vector<double> ratios;
    double retained_baseline_ns = 0.0;
    double retained_tiny_ns = 0.0;
    uint32_t retained_samples = 0u;
    for (uint32_t cycle = 0u; cycle < kCycles; ++cycle)
    {
        double baseline_ns = 0.0;
        double tiny_ns = 0.0;
        uint32_t baseline_samples = 0u;
        uint32_t tiny_samples = 0u;
        for (char arm : kOrder)
        {
            double sample_ns = 0.0;
            const bool measured =
                arm == 'A' ?
                    MeasureBaselineBatch(
                        baseline, fixture, kIterations, sample_ns) :
                    MeasureTinyBatch(
                        fixture, kIterations, sample_ns);
            if (!measured) {
                return false;
            }
            if (arm == 'A')
            {
                baseline_ns += sample_ns;
                ++baseline_samples;
            }
            else
            {
                tiny_ns += sample_ns;
                ++tiny_samples;
            }
        }
        if (baseline_samples != 4u || tiny_samples != 4u) {
            return false;
        }
        baseline_ns /= baseline_samples;
        tiny_ns /= tiny_samples;
        // Cycle zero warms allocator and instruction/data caches but is never
        // used for acceptance or the reported aggregate.
        if (cycle != 0u)
        {
            ratios.push_back(tiny_ns / baseline_ns);
            retained_baseline_ns += baseline_ns;
            retained_tiny_ns += tiny_ns;
            ++retained_samples;
        }
    }
    if (ratios.size() != kCycles - 1u || retained_samples == 0u) {
        return false;
    }
    const double median_ratio = Median(ratios);
    const double worst_ratio =
        *std::max_element(ratios.begin(), ratios.end());
    retained_baseline_ns /= retained_samples;
    retained_tiny_ns /= retained_samples;
    const PacketPair& baseline_packets =
        baseline == BaselineKind::V2Precode ?
            fixture.V2Packets : fixture.LegacyPackets;
    const char* const baseline_name =
        baseline == BaselineKind::V2Precode ?
            "wh2-precode" : "wirehair-legacy";
    const bool passed =
        median_ratio <= kMedianRatioLimit &&
        worst_ratio <= kWorstRatioLimit;
    std::printf(
        "tiny_mds_speed,scope=fresh-decoder-plus-two-repairs,"
        "order=ABBABAAB,discard_cycles=1,cycles=%u,iterations=%u,"
        "baseline=%s,bb=%u,baseline_ids=%u:%u,tiny_ids=%u:%u,"
        "baseline_mean_ns=%.2f,tiny_mean_ns=%.2f,"
        "median_ratio=%.4f,worst_ratio=%.4f,"
        "median_limit=%.2f,worst_limit=%.2f,result=%s\n",
        kCycles, kIterations, baseline_name, fixture.BlockBytes,
        baseline_packets.FirstId, baseline_packets.SecondId,
        fixture.TinyPackets.FirstId, fixture.TinyPackets.SecondId,
        retained_baseline_ns, retained_tiny_ns,
        median_ratio, worst_ratio,
        kMedianRatioLimit, kWorstRatioLimit,
        passed ? "PASS" : "FAIL");
    return passed;
}

bool RunPairedSpeedGate()
{
    const uint32_t block_bytes_cases[] = {2u, 8u, 64u};
    for (uint32_t block_bytes : block_bytes_cases)
    {
        SpeedFixture fixture;
        if (!PrepareSpeedFixture(block_bytes, fixture))
        {
            std::fprintf(stderr,
                "tiny MDS speed fixture failed for bb=%u\n", block_bytes);
            return false;
        }
        if (!RunPairedSpeedComparison(
                fixture, BaselineKind::V2Precode) ||
            !RunPairedSpeedComparison(
                fixture, BaselineKind::Legacy))
        {
            std::fprintf(stderr,
                "tiny MDS paired speed gate failed for bb=%u\n",
                block_bytes);
            return false;
        }
    }
    std::puts("V2 tiny MDS paired speed gate: PASS");
    return true;
}

bool CheckK1()
{
    const uint32_t block_bytes = 16u;
    std::vector<uint8_t> message(7u);
    FillMessage(message);
    wirehair_v2::TinyMdsCodec encoder;
    if (!Check(encoder.InitializeEncoder(
            message.data(), message.size(), block_bytes) == Wirehair_Success,
            "K1 encoder initialization"))
    {
        return false;
    }

    const uint32_t packet_ids[] = {
        0u, 1u, 2u, 256u, 257u, UINT32_MAX
    };
    for (uint32_t packet_id : packet_ids)
    {
        std::vector<uint8_t> packet(block_bytes, 0xa5u);
        uint32_t packet_bytes = UINT32_MAX;
        const uint32_t expected_bytes =
            packet_id == 0u ? (uint32_t)message.size() : block_bytes;
        if (!Check(encoder.EncodeResult(
                packet_id, packet.data(), packet.size(), &packet_bytes) ==
                    Wirehair_Success &&
                packet_bytes == expected_bytes &&
                std::memcmp(
                    packet.data(), message.data(), message.size()) == 0,
                "K1 packet encoding"))
        {
            return false;
        }
        if (packet_id != 0u &&
            !Check(std::all_of(
                packet.begin() + message.size(), packet.end(),
                [](uint8_t byte) { return byte == 0u; }),
                "K1 canonical repair padding"))
        {
            return false;
        }

        wirehair_v2::TinyMdsCodec decoder;
        std::vector<uint8_t> recovered(message.size(), 0xa5u);
        if (!Check(decoder.InitializeDecoder(
                message.size(), block_bytes) == Wirehair_Success,
                "K1 decoder initialization") ||
            !Check(decoder.DecodeResult(
                packet_id, packet.data(), packet_bytes) == Wirehair_Success &&
                decoder.IsDecoded() && decoder.ReceivedCount() == 1u,
                "K1 decode from arbitrary id") ||
            !Check(decoder.RecoverResult(
                recovered.data(), recovered.size()) == Wirehair_Success &&
                recovered == message,
                "K1 recovery"))
        {
            return false;
        }
    }

    std::vector<uint8_t> corrupt(block_bytes, 0u);
    std::memcpy(corrupt.data(), message.data(), message.size());
    corrupt.back() = 1u;
    wirehair_v2::TinyMdsCodec decoder;
    std::vector<uint8_t> untouched(message.size(), 0xa5u);
    const std::vector<uint8_t> before = untouched;
    return Check(decoder.InitializeDecoder(
            message.size(), block_bytes) == Wirehair_Success,
            "K1 padding decoder initialization") &&
        Check(decoder.DecodeResult(
            1u, corrupt.data(), corrupt.size()) == Wirehair_Error &&
            !decoder.IsDecoded() && decoder.ReceivedCount() == 0u,
            "K1 noncanonical repair padding rejected") &&
        Check(decoder.RecoverResult(
            untouched.data(), untouched.size()) == Wirehair_NeedMore &&
            untouched == before,
            "K1 failed decode recovery is no-write");
}

bool CheckK2BoundaryAndDuplicates()
{
    const uint32_t block_bytes = 16u;
    std::vector<uint8_t> message(31u);
    FillMessage(message);
    wirehair_v2::TinyMdsCodec encoder;
    if (!Check(encoder.InitializeEncoder(
            message.data(), message.size(), block_bytes) == Wirehair_Success,
            "K2 boundary encoder initialization"))
    {
        return false;
    }

    std::vector<uint8_t> packet2(block_bytes);
    std::vector<uint8_t> packet256(block_bytes);
    uint32_t packet2_bytes = 0u;
    uint32_t packet256_bytes = 0u;
    if (!Check(encoder.EncodeResult(
            2u, packet2.data(), packet2.size(), &packet2_bytes) ==
                Wirehair_Success && packet2_bytes == block_bytes,
            "K2 packet 2") ||
        !Check(encoder.EncodeResult(
            256u, packet256.data(), packet256.size(), &packet256_bytes) ==
                Wirehair_Success && packet256_bytes == block_bytes,
            "K2 packet 256"))
    {
        return false;
    }

    std::vector<uint8_t> invalid(block_bytes, 0xa5u);
    const std::vector<uint8_t> invalid_before = invalid;
    uint32_t invalid_bytes = UINT32_C(0x12345678);
    if (!Check(encoder.EncodeResult(
            257u, invalid.data(), invalid.size(), &invalid_bytes) ==
                Wirehair_InvalidInput &&
            invalid == invalid_before &&
            invalid_bytes == UINT32_C(0x12345678),
            "K2 id 257 encode rejection") ||
        !Check(encoder.EncodeResult(
            2u, invalid.data(), block_bytes - 1u, &invalid_bytes) ==
                Wirehair_InvalidInput &&
            invalid == invalid_before &&
            invalid_bytes == UINT32_C(0x12345678),
            "K2 short encode no-write"))
    {
        return false;
    }

    wirehair_v2::TinyMdsCodec decoder;
    if (!Check(decoder.InitializeDecoder(
            message.size(), block_bytes) == Wirehair_Success,
            "K2 boundary decoder initialization") ||
        !Check(decoder.DecodeResult(
            257u, invalid.data(), invalid.size()) == Wirehair_InvalidInput &&
            decoder.ReceivedCount() == 0u,
            "K2 id 257 decode rejection") ||
        !Check(decoder.DecodeResult(
            2u, packet2.data(), packet2_bytes) == Wirehair_NeedMore &&
            decoder.ReceivedCount() == 1u,
            "K2 first packet") ||
        !Check(decoder.DecodeResult(
            2u, packet2.data(), packet2_bytes) == Wirehair_NeedMore &&
            decoder.ReceivedCount() == 1u,
            "K2 identical duplicate"))
    {
        return false;
    }

    std::vector<uint8_t> conflict = packet2;
    conflict[0] ^= 1u;
    if (!Check(decoder.DecodeResult(
            2u, conflict.data(), conflict.size()) ==
                Wirehair_InvalidInput &&
            decoder.ReceivedCount() == 1u && !decoder.IsDecoded(),
            "K2 conflicting duplicate") ||
        !Check(decoder.DecodeResult(
            256u, packet256.data(), packet256_bytes) == Wirehair_Success &&
            decoder.IsDecoded() && decoder.ReceivedCount() == 2u,
            "K2 second independent packet"))
    {
        return false;
    }

    std::vector<uint8_t> recovered(message.size(), 0xa5u);
    if (!Check(decoder.RecoverResult(
            recovered.data(), recovered.size()) == Wirehair_Success &&
            recovered == message,
            "K2 boundary recovery") ||
        !Check(decoder.DecodeResult(
            2u, packet2.data(), packet2_bytes) == Wirehair_Success,
            "K2 completed duplicate validation") ||
        !Check(decoder.DecodeResult(
            2u, conflict.data(), conflict.size()) == Wirehair_Error,
            "K2 completed conflict validation") ||
        !Check(decoder.RecoverResult(
            recovered.data(), recovered.size()) == Wirehair_Success &&
            recovered == message,
            "K2 completed conflict does not corrupt solution"))
    {
        return false;
    }

    if (!Check(decoder.DecodeResult(
            1u, packet2.data(), block_bytes) == Wirehair_InvalidInput,
            "K2 wrong systematic payload length rejected after completion"))
    {
        return false;
    }

    std::vector<uint8_t> noncanonical = packet256;
    noncanonical.back() ^= 1u;
    wirehair_v2::TinyMdsCodec padding_decoder;
    recovered.assign(message.size(), 0xa5u);
    return Check(padding_decoder.InitializeDecoder(
            message.size(), block_bytes) == Wirehair_Success,
            "K2 padding decoder initialization") &&
        Check(padding_decoder.DecodeResult(
            2u, packet2.data(), packet2_bytes) == Wirehair_NeedMore,
            "K2 padding first packet") &&
        Check(padding_decoder.DecodeResult(
            256u, noncanonical.data(), packet256_bytes) == Wirehair_Error &&
            !padding_decoder.IsDecoded() &&
            padding_decoder.ReceivedCount() == 1u,
            "K2 noncanonical solved padding rejected") &&
        Check(padding_decoder.DecodeResult(
            256u, packet256.data(), packet256_bytes) == Wirehair_Success &&
            padding_decoder.IsDecoded() &&
            padding_decoder.ReceivedCount() == 2u,
            "K2 valid replacement after padding error") &&
        Check(padding_decoder.RecoverResult(
            recovered.data(), recovered.size()) == Wirehair_Success &&
            recovered == message,
            "K2 recovery after padding-error retry");
}

bool CheckK2Exhaustive(uint32_t block_bytes, uint64_t message_bytes)
{
    std::vector<uint8_t> message((size_t)message_bytes);
    FillMessage(message);
    wirehair_v2::TinyMdsCodec encoder;
    if (encoder.InitializeEncoder(
            message.data(), message.size(), block_bytes) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "tiny MDS exhaustive encoder initialization failed "
            "(bb=%u message=%llu)\n",
            block_bytes, (unsigned long long)message_bytes);
        return false;
    }

    static const uint32_t symbol_count =
        wirehair_v2::TinyMdsCodec::kMaxK2PacketId + 1u;
    std::vector<uint8_t> symbols(
        (size_t)symbol_count * block_bytes, uint8_t{0});
    std::vector<uint32_t> symbol_bytes(symbol_count, 0u);
    for (uint32_t id = 0u; id < symbol_count; ++id)
    {
        if (encoder.EncodeResult(
                id,
                symbols.data() + (size_t)id * block_bytes,
                block_bytes,
                &symbol_bytes[id]) != Wirehair_Success)
        {
            std::fprintf(stderr,
                "tiny MDS exhaustive encode failed at id=%u\n", id);
            return false;
        }
    }

    std::vector<uint8_t> recovered(message.size());
    uint64_t pair_count = 0u;
    for (uint32_t first = 0u; first < symbol_count; ++first)
    {
        for (uint32_t second = first + 1u;
             second < symbol_count;
             ++second)
        {
            for (unsigned reverse = 0u; reverse < 2u; ++reverse)
            {
                const uint32_t id0 = reverse ? second : first;
                const uint32_t id1 = reverse ? first : second;
                wirehair_v2::TinyMdsCodec decoder;
                if (decoder.InitializeDecoder(
                        message.size(), block_bytes) != Wirehair_Success ||
                    decoder.DecodeResult(
                        id0,
                        symbols.data() + (size_t)id0 * block_bytes,
                        symbol_bytes[id0]) != Wirehair_NeedMore ||
                    decoder.DecodeResult(
                        id1,
                        symbols.data() + (size_t)id1 * block_bytes,
                        symbol_bytes[id1]) != Wirehair_Success ||
                    decoder.RecoverResult(
                        recovered.data(), recovered.size()) !=
                            Wirehair_Success ||
                    recovered != message)
                {
                    std::fprintf(stderr,
                        "tiny MDS exhaustive pair failed "
                        "(bb=%u message=%llu ids=%u,%u)\n",
                        block_bytes,
                        (unsigned long long)message_bytes,
                        id0, id1);
                    return false;
                }
                ++pair_count;
            }
        }
    }
    const uint64_t expected_pairs =
        (uint64_t)symbol_count * (symbol_count - 1u);
    if (pair_count != expected_pairs)
    {
        std::fprintf(stderr,
            "tiny MDS exhaustive pair count mismatch: %llu != %llu\n",
            (unsigned long long)pair_count,
            (unsigned long long)expected_pairs);
        return false;
    }
    return true;
}

bool CheckAllocationFailures()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t block_bytes = 16u;
    std::vector<uint8_t> message(31u);
    FillMessage(message);

    wirehair_v2::TinyMdsCodec empty;
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(0);
    const WirehairResult cold_oom = empty.InitializeEncoder(
        message.data(), message.size(), block_bytes);
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(-1);
    if (!Check(cold_oom == Wirehair_OOM &&
            !empty.IsEncoder() && !empty.IsDecoder(),
            "cold encoder allocation OOM"))
    {
        return false;
    }
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(0);
    const WirehairResult cold_decoder_oom =
        empty.InitializeDecoder(message.size(), block_bytes);
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(-1);
    if (!Check(cold_decoder_oom == Wirehair_OOM &&
            !empty.IsEncoder() && !empty.IsDecoder(),
            "cold decoder allocation OOM"))
    {
        return false;
    }

    wirehair_v2::TinyMdsCodec codec;
    if (!Check(codec.InitializeEncoder(
            message.data(), message.size(), block_bytes) == Wirehair_Success,
            "transactional OOM baseline encoder"))
    {
        return false;
    }
    std::vector<uint8_t> before(block_bytes);
    uint32_t before_bytes = 0u;
    if (!Check(codec.EncodeResult(
            2u, before.data(), before.size(), &before_bytes) ==
                Wirehair_Success,
            "transactional OOM baseline packet"))
    {
        return false;
    }

    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(0);
    const WirehairResult reinit_oom =
        codec.InitializeDecoder(message.size(), block_bytes);
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(-1);
    std::vector<uint8_t> after(block_bytes);
    uint32_t after_bytes = 0u;
    if (!Check(reinit_oom == Wirehair_OOM && codec.IsEncoder(),
            "transactional decoder reinit OOM state") ||
        !Check(codec.EncodeResult(
            2u, after.data(), after.size(), &after_bytes) ==
                Wirehair_Success &&
            after_bytes == before_bytes && after == before,
            "transactional decoder reinit OOM bytes"))
    {
        return false;
    }

    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(0);
    const WirehairResult encoder_reinit_oom = codec.InitializeEncoder(
        message.data(), message.size(), block_bytes);
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(-1);
    if (!Check(encoder_reinit_oom == Wirehair_OOM && codec.IsEncoder(),
            "transactional encoder reinit OOM"))
    {
        return false;
    }

    std::vector<uint8_t> second(block_bytes);
    uint32_t second_bytes = 0u;
    if (!Check(codec.EncodeResult(
            256u, second.data(), second.size(), &second_bytes) ==
                Wirehair_Success,
            "transactional OOM second packet"))
    {
        return false;
    }
    wirehair_v2::TinyMdsCodec decoder;
    if (!Check(decoder.InitializeDecoder(
            message.size(), block_bytes) == Wirehair_Success,
            "transactional OOM baseline decoder") ||
        !Check(decoder.DecodeResult(
            2u, before.data(), before_bytes) == Wirehair_NeedMore,
            "transactional OOM retained first packet"))
    {
        return false;
    }
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(0);
    const WirehairResult decoder_to_encoder_oom = decoder.InitializeEncoder(
        message.data(), message.size(), block_bytes);
    wirehair_v2::SetTinyMdsAllocationFailureCountdownForTesting(-1);
    return Check(decoder_to_encoder_oom == Wirehair_OOM &&
            decoder.IsDecoder() && decoder.ReceivedCount() == 1u &&
            !decoder.IsDecoded(),
            "transactional encoder reinit preserves decoder") &&
        Check(decoder.DecodeResult(
            256u, second.data(), second_bytes) == Wirehair_Success &&
            decoder.IsDecoded(),
            "transactional OOM decoder remains usable");
#else
    return true;
#endif
}

bool CheckInvalidInputsAndModes()
{
    const uint32_t block_bytes = 16u;
    std::vector<uint8_t> message(31u);
    FillMessage(message);
    std::vector<uint8_t> packet(block_bytes, 0xa5u);
    const std::vector<uint8_t> packet_before = packet;
    uint32_t packet_bytes = UINT32_C(0x12345678);

    wirehair_v2::TinyMdsCodec empty;
    if (!Check(empty.InitializeEncoder(
            nullptr, 1u, 1u) == Wirehair_InvalidInput,
            "null encoder message rejected") ||
        !Check(empty.InitializeEncoder(
            message.data(), 0u, block_bytes) == Wirehair_InvalidInput,
            "zero-byte encoder message rejected") ||
        !Check(empty.InitializeEncoder(
            message.data(), 1u, 0u) == Wirehair_InvalidInput,
            "zero encoder block size rejected") ||
        !Check(empty.InitializeEncoder(
            message.data(), 1u, UINT32_C(0x80000000)) ==
                Wirehair_InvalidInput,
            "oversized encoder block rejected") ||
        !Check(empty.InitializeEncoder(
            message.data(), 33u, block_bytes) == Wirehair_InvalidInput,
            "K3 encoder rejected") ||
        !Check(empty.InitializeDecoder(
            0u, block_bytes) == Wirehair_InvalidInput,
            "zero-byte decoder message rejected") ||
        !Check(empty.InitializeDecoder(
            1u, 0u) == Wirehair_InvalidInput,
            "zero decoder block size rejected") ||
        !Check(empty.InitializeDecoder(
            1u, UINT32_C(0x80000000)) == Wirehair_InvalidInput,
            "oversized decoder block rejected") ||
        !Check(empty.InitializeDecoder(
            33u, block_bytes) == Wirehair_InvalidInput,
            "K3 decoder rejected") ||
        !Check(empty.EncodeResult(
            2u, packet.data(), packet.size(), &packet_bytes) ==
                Wirehair_InvalidInput &&
            packet == packet_before &&
            packet_bytes == UINT32_C(0x12345678),
            "uninitialized encode is no-write") ||
        !Check(empty.DecodeResult(
            2u, packet.data(), packet.size()) == Wirehair_InvalidInput,
            "uninitialized decode rejected"))
    {
        return false;
    }
    std::vector<uint8_t> output(message.size(), 0xa5u);
    const std::vector<uint8_t> output_before = output;
    if (!Check(empty.RecoverResult(
            output.data(), output.size()) == Wirehair_InvalidInput &&
            output == output_before,
            "uninitialized recovery is no-write"))
    {
        return false;
    }

    wirehair_v2::TinyMdsCodec encoder;
    if (!Check(encoder.InitializeEncoder(
            message.data(), message.size(), block_bytes) == Wirehair_Success,
            "mode-test encoder initialization") ||
        !Check(encoder.DecodeResult(
            2u, packet.data(), packet.size()) == Wirehair_InvalidInput,
            "encoder rejects decode") ||
        !Check(encoder.RecoverResult(
            output.data(), output.size()) == Wirehair_InvalidInput &&
            output == output_before,
            "encoder rejects recovery without writes"))
    {
        return false;
    }
    std::vector<uint8_t> reference_packet(block_bytes, 0u);
    uint32_t reference_bytes = 0u;
    if (!Check(encoder.EncodeResult(
            2u,
            reference_packet.data(),
            reference_packet.size(),
            &reference_bytes) == Wirehair_Success,
            "mode-test reference packet") ||
        !Check(encoder.InitializeDecoder(
            0u, block_bytes) == Wirehair_InvalidInput &&
            encoder.IsEncoder(),
            "invalid decoder reinit preserves encoder") ||
        !Check(encoder.InitializeEncoder(
            nullptr, message.size(), block_bytes) == Wirehair_InvalidInput &&
            encoder.IsEncoder(),
            "invalid encoder reinit preserves encoder"))
    {
        return false;
    }
    std::vector<uint8_t> after_invalid(block_bytes, 0u);
    uint32_t after_invalid_bytes = 0u;
    if (!Check(encoder.EncodeResult(
            2u,
            after_invalid.data(),
            after_invalid.size(),
            &after_invalid_bytes) == Wirehair_Success &&
            after_invalid_bytes == reference_bytes &&
            after_invalid == reference_packet,
            "invalid reinitialization preserves encoded bytes"))
    {
        return false;
    }

    wirehair_v2::TinyMdsCodec decoder;
    packet.assign(block_bytes, 0xa5u);
    packet_bytes = UINT32_C(0x12345678);
    if (!Check(decoder.InitializeDecoder(
            message.size(), block_bytes) == Wirehair_Success,
            "mode-test decoder initialization") ||
        !Check(decoder.EncodeResult(
            2u, packet.data(), packet.size(), &packet_bytes) ==
                Wirehair_InvalidInput &&
            packet == packet_before &&
            packet_bytes == UINT32_C(0x12345678),
            "decoder rejects encode without writes") ||
        !Check(decoder.DecodeResult(
            2u, nullptr, block_bytes) == Wirehair_InvalidInput &&
            decoder.ReceivedCount() == 0u,
            "decoder rejects null packet") ||
        !Check(decoder.DecodeResult(
            2u, packet.data(), block_bytes - 1u) ==
                Wirehair_InvalidInput &&
            decoder.ReceivedCount() == 0u,
            "decoder rejects wrong packet length"))
    {
        return false;
    }
    output.assign(message.size(), 0xa5u);
    const std::vector<uint8_t> decoder_output_before = output;
    return Check(decoder.RecoverResult(
            output.data(), output.size()) == Wirehair_NeedMore &&
            output == decoder_output_before,
            "incomplete recovery is no-write") &&
        Check(decoder.RecoverResult(
            output.data(), output.size() - 1u) == Wirehair_InvalidInput &&
            output == decoder_output_before,
            "wrong-size recovery is no-write") &&
        Check(decoder.RecoverResult(
            nullptr, output.size()) == Wirehair_InvalidInput,
            "null recovery output rejected");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], "--benchmark") == 0) {
        return RunPairedSpeedGate() ? 0 : 1;
    }
    if (argc != 1)
    {
        std::fprintf(stderr,
            "usage: V2TinyMDSTest [--benchmark]\n");
        return 2;
    }

    if (!CheckK1() ||
        !CheckK2BoundaryAndDuplicates() ||
        !CheckAllocationFailures() ||
        !CheckInvalidInputsAndModes())
    {
        return 1;
    }

    const uint32_t block_bytes_cases[] = {
        1u, 2u, 3u, 8u, 17u, 64u
    };
    for (uint32_t block_bytes : block_bytes_cases)
    {
        // Sweep every possible final-source length, including the full tail,
        // against every ordered pair of distinct projective directions.
        for (uint32_t tail_bytes = 1u;
             tail_bytes <= block_bytes;
             ++tail_bytes)
        {
            if (!CheckK2Exhaustive(
                    block_bytes,
                    (uint64_t)block_bytes + tail_bytes))
            {
                return 1;
            }
        }
    }

    std::puts(
        "V2 tiny MDS: K1/boundaries/faults/OOM and every ordered K2 "
        "P^1(GF256) pair at every tested tail length PASS");
    return 0;
}
