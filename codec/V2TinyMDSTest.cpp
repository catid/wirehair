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

inline uint32_t RotateRight32(uint32_t value, unsigned bits)
{
    return (value >> bits) | (value << (32u - bits));
}

/**
    Test-local SHA-256 for the tiny-MDS wire-mapping golden.

    Keeping this implementation in the test prevents the compatibility witness
    from depending on the production equation-fingerprint code.  The mapping
    assertions below likewise derive their expected rows directly from the
    documented public contract rather than calling TinyMdsCodec's private
    coefficient selector.
*/
class GoldenSha256
{
public:
    GoldenSha256()
    {
        State[0] = 0x6a09e667u;
        State[1] = 0xbb67ae85u;
        State[2] = 0x3c6ef372u;
        State[3] = 0xa54ff53au;
        State[4] = 0x510e527fu;
        State[5] = 0x9b05688cu;
        State[6] = 0x1f83d9abu;
        State[7] = 0x5be0cd19u;
    }

    void Update(const void* data, size_t bytes)
    {
        const uint8_t* input = static_cast<const uint8_t*>(data);
        TotalBytes += bytes;
        if (BufferedBytes != 0u)
        {
            while (bytes != 0u && BufferedBytes < sizeof(Buffer)) {
                Buffer[BufferedBytes++] = *input++;
                --bytes;
            }
            if (BufferedBytes != sizeof(Buffer)) {
                return;
            }
            ProcessBlock(Buffer);
            BufferedBytes = 0u;
        }
        while (bytes >= sizeof(Buffer))
        {
            ProcessBlock(input);
            input += sizeof(Buffer);
            bytes -= sizeof(Buffer);
        }
        while (bytes != 0u) {
            Buffer[BufferedBytes++] = *input++;
            --bytes;
        }
    }

    void Finalize(uint8_t digest[32])
    {
        const uint64_t total_bits = TotalBytes * 8u;
        const uint8_t one = 0x80u;
        Update(&one, 1u);
        const uint8_t zero = 0u;
        while (BufferedBytes != 56u) {
            Update(&zero, 1u);
        }
        uint8_t length[8];
        for (unsigned i = 0u; i < 8u; ++i) {
            length[i] = (uint8_t)(total_bits >> (56u - 8u * i));
        }
        Update(length, sizeof(length));
        for (unsigned i = 0u; i < 8u; ++i)
        {
            digest[4u * i] = (uint8_t)(State[i] >> 24);
            digest[4u * i + 1u] = (uint8_t)(State[i] >> 16);
            digest[4u * i + 2u] = (uint8_t)(State[i] >> 8);
            digest[4u * i + 3u] = (uint8_t)State[i];
        }
    }

private:
    void ProcessBlock(const uint8_t block[64])
    {
        static const uint32_t RoundConstants[64] = {
            0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
            0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
            0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
            0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
            0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
            0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
            0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
            0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
            0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
            0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
            0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
            0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
            0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
            0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
            0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
            0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
        };
        uint32_t words[64];
        for (unsigned i = 0u; i < 16u; ++i)
        {
            words[i] = ((uint32_t)block[4u * i] << 24) |
                ((uint32_t)block[4u * i + 1u] << 16) |
                ((uint32_t)block[4u * i + 2u] << 8) |
                (uint32_t)block[4u * i + 3u];
        }
        for (unsigned i = 16u; i < 64u; ++i)
        {
            const uint32_t s0 =
                RotateRight32(words[i - 15u], 7u) ^
                RotateRight32(words[i - 15u], 18u) ^
                (words[i - 15u] >> 3);
            const uint32_t s1 =
                RotateRight32(words[i - 2u], 17u) ^
                RotateRight32(words[i - 2u], 19u) ^
                (words[i - 2u] >> 10);
            words[i] = words[i - 16u] + s0 + words[i - 7u] + s1;
        }

        uint32_t a = State[0];
        uint32_t b = State[1];
        uint32_t c = State[2];
        uint32_t d = State[3];
        uint32_t e = State[4];
        uint32_t f = State[5];
        uint32_t g = State[6];
        uint32_t h = State[7];
        for (unsigned i = 0u; i < 64u; ++i)
        {
            const uint32_t sum1 =
                RotateRight32(e, 6u) ^ RotateRight32(e, 11u) ^
                RotateRight32(e, 25u);
            const uint32_t choose = (e & f) ^ (~e & g);
            const uint32_t temp1 =
                h + sum1 + choose + RoundConstants[i] + words[i];
            const uint32_t sum0 =
                RotateRight32(a, 2u) ^ RotateRight32(a, 13u) ^
                RotateRight32(a, 22u);
            const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t temp2 = sum0 + majority;
            h = g;
            g = f;
            f = e;
            e = d + temp1;
            d = c;
            c = b;
            b = a;
            a = temp1 + temp2;
        }
        State[0] += a;
        State[1] += b;
        State[2] += c;
        State[3] += d;
        State[4] += e;
        State[5] += f;
        State[6] += g;
        State[7] += h;
    }

    uint32_t State[8] = {};
    uint64_t TotalBytes = 0u;
    uint8_t Buffer[64] = {};
    size_t BufferedBytes = 0u;
};

void AppendU32LE(std::vector<uint8_t>& stream, uint32_t value)
{
    stream.push_back((uint8_t)value);
    stream.push_back((uint8_t)(value >> 8));
    stream.push_back((uint8_t)(value >> 16));
    stream.push_back((uint8_t)(value >> 24));
}

void AppendU64LE(std::vector<uint8_t>& stream, uint64_t value)
{
    AppendU32LE(stream, (uint32_t)value);
    AppendU32LE(stream, (uint32_t)(value >> 32));
}

void FormatDigestHex(const uint8_t digest[32], char output[65])
{
    static const char HexDigits[] = "0123456789abcdef";
    for (unsigned i = 0u; i < 32u; ++i)
    {
        output[2u * i] = HexDigits[digest[i] >> 4];
        output[2u * i + 1u] = HexDigits[digest[i] & 15u];
    }
    output[64] = '\0';
}

bool CheckGoldenSha256KnownAnswers()
{
    uint8_t digest[32];
    char digest_hex[65];
    GoldenSha256 empty;
    empty.Finalize(digest);
    FormatDigestHex(digest, digest_hex);
    if (std::strcmp(
            digest_hex,
            "e3b0c44298fc1c149afbf4c8996fb924"
            "27ae41e4649b934ca495991b7852b855") != 0)
    {
        std::fprintf(stderr,
            "tiny MDS golden SHA-256 empty-input KAT failed: %s\n",
            digest_hex);
        return false;
    }

    static const char Abc[] = "abc";
    GoldenSha256 abc;
    abc.Update(Abc, sizeof(Abc) - 1u);
    abc.Finalize(digest);
    FormatDigestHex(digest, digest_hex);
    if (std::strcmp(
            digest_hex,
            "ba7816bf8f01cfea414140de5dae2223"
            "b00361a396177a9cb410ff61f20015ad") != 0)
    {
        std::fprintf(stderr,
            "tiny MDS golden SHA-256 abc KAT failed: %s\n",
            digest_hex);
        return false;
    }
    return true;
}

void ExpectedK2Coefficients(
    uint32_t id,
    uint8_t& alpha,
    uint8_t& beta)
{
    alpha = 1u;
    beta = 0u;
    if (id == 1u)
    {
        alpha = 0u;
        beta = 1u;
    }
    else if (id >= 2u) {
        beta = (uint8_t)(id - 1u);
    }
}

uint8_t ReferenceGf256Multiply(uint8_t left, uint8_t right)
{
    // Public tiny-MDS arithmetic uses the byte polynomial basis modulo 0x14d.
    uint16_t multiplicand = left;
    uint8_t multiplier = right;
    uint16_t product = 0u;
    while (multiplier != 0u)
    {
        if ((multiplier & 1u) != 0u) {
            product ^= multiplicand;
        }
        multiplier >>= 1;
        multiplicand <<= 1;
        if ((multiplicand & UINT16_C(0x100)) != 0u) {
            multiplicand ^= UINT16_C(0x14d);
        }
    }
    return (uint8_t)product;
}

struct PublicV2CodecGuard
{
    ~PublicV2CodecGuard()
    {
        wirehair_v2_free(Handle);
    }

    WirehairV2Codec Handle = nullptr;
};

bool AppendK2WireGoldenFixture(
    uint32_t fixture_kind,
    const std::vector<uint8_t>& message,
    uint32_t block_bytes,
    std::vector<uint8_t>& stream)
{
    static const uint32_t PacketCount = 257u;
    if (message.size() != (size_t)block_bytes * 2u) {
        return false;
    }

    uint8_t descriptor[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t descriptor_bytes = 0u;
    PublicV2CodecGuard encoder;
    if (wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_TINY_MDS_2026_07,
            message.data(), message.size(), block_bytes,
            descriptor, sizeof(descriptor), &descriptor_bytes,
            &encoder.Handle) != WirehairV2_Success ||
        descriptor_bytes != sizeof(descriptor) ||
        !encoder.Handle)
    {
        std::fprintf(stderr,
            "tiny MDS K2 wire-golden public encoder initialization "
            "failed for fixture=%u\n",
            fixture_kind);
        return false;
    }
    PublicV2CodecGuard decoder;
    if (wirehair_v2_decoder_create(
            descriptor, descriptor_bytes, &decoder.Handle) !=
                WirehairV2_Success ||
        !decoder.Handle)
    {
        std::fprintf(stderr,
            "tiny MDS K2 wire-golden public decoder initialization "
            "failed for fixture=%u\n",
            fixture_kind);
        return false;
    }

    // Fixture record:
    //   kind (u32 LE); message bytes (u64 LE); block bytes (u32 LE);
    //   descriptor bytes (u32 LE) and exact descriptor; exact source message;
    //   packet count (u32 LE);
    //   for ids 0..256: id (u32 LE), expected alpha (u8), expected beta
    //   (u8), encoded byte count (u32 LE), and exact encoded bytes.
    AppendU32LE(stream, fixture_kind);
    AppendU64LE(stream, message.size());
    AppendU32LE(stream, block_bytes);
    AppendU32LE(stream, descriptor_bytes);
    stream.insert(
        stream.end(), descriptor, descriptor + descriptor_bytes);
    stream.insert(stream.end(), message.begin(), message.end());
    AppendU32LE(stream, PacketCount);

    std::vector<uint8_t> packet(block_bytes, 0xa5u);
    for (uint32_t id = 0u; id < PacketCount; ++id)
    {
        uint8_t expected_alpha = 0u;
        uint8_t expected_beta = 0u;
        ExpectedK2Coefficients(
            id, expected_alpha, expected_beta);
        std::fill(packet.begin(), packet.end(), uint8_t{0xa5u});
        uint32_t packet_bytes = 0u;
        if (wirehair_v2_encode(
                encoder.Handle, id, packet.data(), block_bytes,
                &packet_bytes) != WirehairV2_Success ||
            packet_bytes != block_bytes)
        {
            std::fprintf(stderr,
                "tiny MDS K2 wire-golden encode failed for "
                "fixture=%u id=%u bytes=%u\n",
                fixture_kind, id, packet_bytes);
            return false;
        }
        for (uint32_t byte = 0u; byte < block_bytes; ++byte)
        {
            const uint8_t expected =
                ReferenceGf256Multiply(
                    expected_alpha, message[byte]) ^
                ReferenceGf256Multiply(
                    expected_beta, message[block_bytes + byte]);
            if (packet[byte] != expected)
            {
                std::fprintf(stderr,
                    "tiny MDS K2 wire mapping/arithmetic drift at "
                    "fixture=%u id=%u byte=%u coefficients=[%u,%u]: "
                    "expected=%u actual=%u\n",
                    fixture_kind, id, byte,
                    expected_alpha, expected_beta,
                    expected, packet[byte]);
                return false;
            }
        }
        AppendU32LE(stream, id);
        stream.push_back(expected_alpha);
        stream.push_back(expected_beta);
        AppendU32LE(stream, packet_bytes);
        stream.insert(
            stream.end(), packet.begin(), packet.begin() + packet_bytes);
    }

    // The accepted-domain digest cannot by itself pin rejection behavior.
    // Keep both the adjacent and maximum uint32_t out-of-range cases explicit.
    const uint32_t invalid_ids[] = {257u, UINT32_MAX};
    for (uint32_t invalid_id : invalid_ids)
    {
        std::fill(packet.begin(), packet.end(), uint8_t{0xa5u});
        const std::vector<uint8_t> before = packet;
        uint32_t packet_bytes = UINT32_C(0x12345678);
        if (wirehair_v2_encode(
                encoder.Handle, invalid_id,
                packet.data(), block_bytes, &packet_bytes) !=
                    WirehairV2_InvalidInput ||
            packet_bytes != block_bytes ||
            packet != before ||
            wirehair_v2_decode(
                decoder.Handle, invalid_id,
                packet.data(), block_bytes) != WirehairV2_InvalidInput)
        {
            std::fprintf(stderr,
                "tiny MDS K2 wire-golden invalid-id rejection drift "
                "for fixture=%u id=%u\n",
                fixture_kind, invalid_id);
            return false;
        }
    }
    return true;
}

bool CheckK2WireMappingGolden()
{
    if (!CheckGoldenSha256KnownAnswers()) {
        return false;
    }

    // Canonical stream version 2:
    //   ASCII "WH2TMDS2" (no terminator);
    //   profile ID (u64 LE); fixture count (u32 LE);
    //   followed by the fixture records documented above.
    std::vector<uint8_t> stream;
    static const uint8_t Tag[8] = {
        'W', 'H', '2', 'T', 'M', 'D', 'S', '2'
    };
    stream.insert(stream.end(), Tag, Tag + sizeof(Tag));
    AppendU64LE(stream, WIREHAIR_V2_PROFILE_TINY_MDS_2026_07);
    AppendU32LE(stream, 2u);

    // In the basis fixture each two-byte packet is exactly [alpha,beta].
    const std::vector<uint8_t> basis_message = {
        1u, 0u, 0u, 1u
    };
    if (!AppendK2WireGoldenFixture(
            1u, basis_message, 2u, stream))
    {
        return false;
    }

    // Source block zero is all zero.  Source block one enumerates every byte,
    // so IDs 2..256 cover every nonzero coefficient times every field value;
    // ID 0/1 also pin zero and identity.  This is the complete 256x256
    // multiplication table under the independently specified 0x14d field.
    std::vector<uint8_t> arithmetic_message(512u, 0u);
    for (uint32_t value = 0u; value < 256u; ++value) {
        arithmetic_message[256u + value] = (uint8_t)value;
    }
    if (!AppendK2WireGoldenFixture(
            2u, arithmetic_message, 256u, stream))
    {
        return false;
    }
    static const size_t FrozenStreamBytes = 72094u;
    if (stream.size() != FrozenStreamBytes)
    {
        std::fprintf(stderr,
            "tiny MDS K2 wire-mapping stream size drift: "
            "expected=%llu actual=%llu\n",
            (unsigned long long)FrozenStreamBytes,
            (unsigned long long)stream.size());
        return false;
    }

    GoldenSha256 sha256;
    sha256.Update(stream.data(), stream.size());
    uint8_t digest[32];
    sha256.Finalize(digest);
    char digest_hex[65];
    FormatDigestHex(digest, digest_hex);
    static const char FrozenDigest[] =
        "03bced599bfe6bb916cc26d15987ee86ed680e129fed45c3b9b897f767eb67ca";
    if (std::strcmp(digest_hex, FrozenDigest) != 0)
    {
        std::fprintf(stderr,
            "tiny MDS K2 wire-mapping golden drift:\n"
            "  frozen:   %s\n"
            "  computed: %s\n"
            "Changing this stream under the published profile ID is a wire "
            "compatibility bug.\n",
            FrozenDigest, digest_hex);
        return false;
    }
    return true;
}

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

    // These are deliberately separate gates: the wire digest pins every full
    // K2 row and field product, while K1 and every partial-tail length retain
    // their own exact length/padding and all-pairs recovery checks below.
    if (!CheckK2WireMappingGolden() ||
        !CheckK1() ||
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
        "V2 tiny MDS: full wire mapping/arithmetic golden, "
        "K1/boundaries/faults/OOM, and every ordered K2 P^1(GF256) "
        "pair at every tested tail length PASS");
    return 0;
}
