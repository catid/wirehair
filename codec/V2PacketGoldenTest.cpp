#include <wirehair/wirehair.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

/*
    Representative packet goldens for the public serialized V2 profiles.

    For each public profile ID at representative K/BlockBytes this pins the
    exact serialized descriptor and exact encoded repair packets emitted by
    the public API.  Together with the all-K fingerprints in
    V2FingerprintTest.cpp this freezes the July 2026 equation contract: a
    changed descriptor or packet byte under an existing profile ID is a
    compatibility bug and must ship under a new profile ID.

    The constants are filled from a build of the frozen implementation:

        wirehair_v2_packet_golden_test --print-goldens

    prints the complete golden block ready to paste over the marked section.
    Unset goldens FAIL the test; they never skip.
*/

namespace {

struct PacketGoldenCase
{
    const char* Name;
    uint64_t ProfileId;
    uint64_t MessageBytes;
    uint32_t BlockBytes;
};

// Small-K cases at 16-byte blocks are already frozen in V2ProfileTest.cpp.
// These cases cover K=1000 and K=10000 (the certified N1 staircase-hit
// boundary) with a partial final block in every message.
const PacketGoldenCase kPacketCases[] = {
    {"certified_k1000_b32",
        WIREHAIR_V2_PROFILE_CERTIFIED_2026_07, 31985u, 32u},
    {"mixed_k1000_b32",
        WIREHAIR_V2_PROFILE_MIXED_2026_07, 31985u, 32u},
    {"mixed_mix2_k1000_b32",
        WIREHAIR_V2_PROFILE_MIXED_MIX2_2026_07, 31985u, 32u},
    {"certified_k10000_b2",
        WIREHAIR_V2_PROFILE_CERTIFIED_2026_07, 19999u, 2u},
    {"mixed_k10000_b2",
        WIREHAIR_V2_PROFILE_MIXED_2026_07, 19999u, 2u},
    {"mixed_mix2_k10000_b2",
        WIREHAIR_V2_PROFILE_MIXED_MIX2_2026_07, 19999u, 2u}
};
const uint32_t kPacketCaseCount =
    (uint32_t)(sizeof(kPacketCases) / sizeof(kPacketCases[0]));

/// Deep repair id shared by every case.
const uint32_t kDeepRepairId = 54321u;

// Per case: serialized descriptor hex, then repair packets K, K+1, and
// kDeepRepairId as lowercase hex.
// --- BEGIN FROZEN V2 PACKET GOLDENS ---
// Fill by building, then running:
//   wirehair_v2_packet_golden_test --print-goldens
// and pasting the printed block over this one.
static const char* const kPacketGoldens[kPacketCaseCount][4] = {
    { // certified_k1000_b32
        "5748563201002000c9f9f447bb5b294bf17c0000000000002000000000000000",
        "83ceaf2f3d105b82c9ae2c3cbe5a81c8a9e6f249673519f4ce142e81f2ff9328",
        "c9e872d3cb2f8196382d32948bde7767fde270d3f8cd4380db0c5784ad3f04a7",
        "4b7f73abe1d8ea4a7872a8e03aeb49797543325ba31215f06f49d6d82f6e7d14"
    },
    { // mixed_k1000_b32
        "5748563201002000b79b6f455dce61e1f17c0000000000002000000000000000",
        "4c4caf8e758b3f3b93d334dbff0e15ce86c4295563a6dad41a039e33953f25d7",
        "0a47a3044585bb5be10e9bd003c9fd8fd1e54a7a1db81c28afdaeda37dd18290",
        "a730980ab08786e3d86b6308bc1503682030b9103b7f2b65e09281205f5e37af"
    },
    { // mixed_mix2_k1000_b32
        "5748563201002000a21206877af2a420f17c0000000000002000000000000000",
        "459ec5533a5a85eb47e525e6ef85c4cad4c36bf3aa6f7249be0b5c712274ef93",
        "6c8f0e017bf42b11093e057d0ba3b18754cb3b1db58943fbd1deb46a40f3fab5",
        "c9a7ce1fc0daf53c46e4f1ed3a2c54958931765a1874aa8532caf0cfcb5d0d8a"
    },
    { // certified_k10000_b2
        "5748563201002000c9f9f447bb5b294b1f4e0000000000000200000000000000",
        "6279",
        "b185",
        "09e3"
    },
    { // mixed_k10000_b2
        "5748563201002000b79b6f455dce61e11f4e0000000000000200000000000000",
        "54dc",
        "089c",
        "b398"
    },
    { // mixed_mix2_k10000_b2
        "5748563201002000a21206877af2a4201f4e0000000000000200000000000000",
        "c23c",
        "2f5f",
        "fb58"
    }
};
// --- END FROZEN V2 PACKET GOLDENS ---

const char* const kGoldenSlotNames[4] = {
    "descriptor", "repair K", "repair K+1", "deep repair"
};

void FillMessage(std::vector<uint8_t>& message)
{
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 73u + 19u);
    }
}

std::string ToHex(const uint8_t* data, size_t bytes)
{
    static const char kHexDigits[] = "0123456789abcdef";
    std::string hex;
    hex.reserve(bytes * 2u);
    for (size_t i = 0; i < bytes; ++i)
    {
        hex += kHexDigits[data[i] >> 4];
        hex += kHexDigits[data[i] & 0xfu];
    }
    return hex;
}

/// Compute the four golden strings for one case through the public API.
/// Also cross-checks that a descriptor-recreated encoder reproduces the
/// exact packet bytes.
bool ComputeCase(
    const PacketGoldenCase& golden_case,
    std::string computed_out[4])
{
    std::vector<uint8_t> message((size_t)golden_case.MessageBytes);
    FillMessage(message);
    const uint32_t block_count = (uint32_t)(
        (golden_case.MessageBytes + golden_case.BlockBytes - 1u) /
        golden_case.BlockBytes);

    uint8_t descriptor[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t descriptor_bytes = 0u;
    WirehairV2Codec encoder = nullptr;
    if (wirehair_v2_encoder_create_profile_id(
            golden_case.ProfileId,
            message.data(), message.size(), golden_case.BlockBytes,
            descriptor, sizeof(descriptor), &descriptor_bytes, &encoder) !=
                WirehairV2_Success ||
        descriptor_bytes != sizeof(descriptor))
    {
        std::fprintf(stderr,
            "packet golden: %s encoder creation failed\n",
            golden_case.Name);
        wirehair_v2_free(encoder);
        return false;
    }
    computed_out[0] = ToHex(descriptor, sizeof(descriptor));

    WirehairV2Codec recreated = nullptr;
    if (wirehair_v2_encoder_create_profile(
            message.data(), descriptor, sizeof(descriptor), &recreated) !=
        WirehairV2_Success)
    {
        std::fprintf(stderr,
            "packet golden: %s descriptor recreation failed\n",
            golden_case.Name);
        wirehair_v2_free(encoder);
        return false;
    }

    const uint32_t packet_ids[3] = {
        block_count, block_count + 1u, kDeepRepairId
    };
    bool ok = true;
    for (uint32_t i = 0; ok && i < 3u; ++i)
    {
        std::vector<uint8_t> block(golden_case.BlockBytes, 0u);
        std::vector<uint8_t> recreated_block(golden_case.BlockBytes, 0u);
        uint32_t bytes = 0u;
        uint32_t recreated_bytes = 0u;
        if (wirehair_v2_encode(
                encoder, packet_ids[i], block.data(),
                (uint32_t)block.size(), &bytes) != WirehairV2_Success ||
            bytes != golden_case.BlockBytes ||
            wirehair_v2_encode(
                recreated, packet_ids[i], recreated_block.data(),
                (uint32_t)recreated_block.size(), &recreated_bytes) !=
                    WirehairV2_Success ||
            recreated_bytes != bytes ||
            std::memcmp(block.data(), recreated_block.data(), bytes) != 0)
        {
            std::fprintf(stderr,
                "packet golden: %s id=%u encode/reproduction failed\n",
                golden_case.Name, packet_ids[i]);
            ok = false;
            break;
        }
        computed_out[1u + i] = ToHex(block.data(), bytes);
    }

    wirehair_v2_free(encoder);
    wirehair_v2_free(recreated);
    return ok;
}

int PrintGoldens()
{
    std::string lines;
    for (uint32_t case_index = 0; case_index < kPacketCaseCount; ++case_index)
    {
        std::string computed[4];
        if (!ComputeCase(kPacketCases[case_index], computed)) {
            return 1;
        }
        lines += "    { // ";
        lines += kPacketCases[case_index].Name;
        lines += "\n";
        for (uint32_t slot = 0; slot < 4u; ++slot)
        {
            lines += "        \"";
            lines += computed[slot];
            lines += slot + 1u < 4u ? "\",\n" : "\"\n";
        }
        lines += case_index + 1u < kPacketCaseCount ? "    },\n" : "    }\n";
    }
    std::printf(
        "// --- BEGIN FROZEN V2 PACKET GOLDENS ---\n"
        "// Fill by building, then running:\n"
        "//   wirehair_v2_packet_golden_test --print-goldens\n"
        "// and pasting the printed block over this one.\n"
        "static const char* const kPacketGoldens[kPacketCaseCount][4] = {\n"
        "%s"
        "};\n"
        "// --- END FROZEN V2 PACKET GOLDENS ---\n",
        lines.c_str());
    return 0;
}

int CheckGoldens()
{
    bool failed = false;
    for (uint32_t case_index = 0; case_index < kPacketCaseCount; ++case_index)
    {
        const PacketGoldenCase& golden_case = kPacketCases[case_index];
        std::string computed[4];
        if (!ComputeCase(golden_case, computed)) {
            return 1;
        }
        for (uint32_t slot = 0; slot < 4u; ++slot)
        {
            const char* golden = kPacketGoldens[case_index][slot];
            if (std::strcmp(golden, "UNSET") == 0)
            {
                std::fprintf(stderr,
                    "packet golden: %s %s is UNSET -- build, run "
                    "'wirehair_v2_packet_golden_test --print-goldens', and "
                    "paste the printed block into "
                    "codec/V2PacketGoldenTest.cpp\n",
                    golden_case.Name, kGoldenSlotNames[slot]);
                failed = true;
                continue;
            }
            if (computed[slot] != golden)
            {
                std::fprintf(stderr,
                    "packet golden DRIFT for %s %s:\n"
                    "  frozen:   %s\n"
                    "  computed: %s\n"
                    "Equation bytes changed under a frozen public profile "
                    "ID.  Ship the change under a NEW profile ID; do not "
                    "update this golden.\n",
                    golden_case.Name, kGoldenSlotNames[slot],
                    golden, computed[slot].c_str());
                failed = true;
            }
        }
    }
    if (failed) {
        return 1;
    }
    std::puts("V2 representative packet goldens: PASS");
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    const bool print_goldens =
        argc == 2 && std::strcmp(argv[1], "--print-goldens") == 0;
    if (argc > 2 || (argc == 2 && !print_goldens))
    {
        std::fprintf(stderr, "usage: %s [--print-goldens]\n", argv[0]);
        return 2;
    }
    return print_goldens ? PrintGoldens() : CheckGoldens();
}
