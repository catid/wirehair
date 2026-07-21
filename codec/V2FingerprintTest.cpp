#include "WirehairV2Fingerprint.h"

#include <wirehair/wirehair.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

/*
    All-K equation-freeze test for the public serialized V2 profiles.

    For every supported K this digests the complete equation-affecting
    expansion of each public profile ID (see WirehairV2Fingerprint.h for the
    exact stream) and compares one SHA-256 per profile against the frozen
    constants below.  Any mismatch means the equations published under an
    existing profile ID changed, which is a compatibility bug: the change
    must ship under a new profile ID instead of new goldens.

    The constants are filled from a build of the frozen implementation:

        wirehair_v2_fingerprint_test --print-goldens

    prints the complete constant block ready to paste over the marked
    section.  Unset constants FAIL the test; they never skip.
*/

namespace {

// --- BEGIN FROZEN V2 EQUATION FINGERPRINTS ---
// Fill by building, then running:
//   wirehair_v2_fingerprint_test --print-goldens
// and pasting the printed block over this one.
static const char kCertifiedAllKFingerprint[] =
    "e6146e7fea89089689a819c72e7e82f799344b451f5e4653b125f622dee3de0b";
static const char kMixedAllKFingerprint[] =
    "47bca161ff7b51684f39d19db3b5b0d11137f21335ec5c74d818d06756d93627";
static const char kMixedMix2AllKFingerprint[] =
    "858321c2c0a07103b2615bb6586ce310105d37dd9b2120eb5b63527a6fcb5404";
// --- END FROZEN V2 EQUATION FINGERPRINTS ---

struct GoldenBinding
{
    uint64_t ProfileId;
    const char* ConstantName;
    const char* Golden;
};

const GoldenBinding kGoldenBindings[] = {
    {
        WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
        "kCertifiedAllKFingerprint",
        kCertifiedAllKFingerprint
    },
    {
        WIREHAIR_V2_PROFILE_MIXED_2026_07,
        "kMixedAllKFingerprint",
        kMixedAllKFingerprint
    },
    {
        WIREHAIR_V2_PROFILE_MIXED_MIX2_2026_07,
        "kMixedMix2AllKFingerprint",
        kMixedMix2AllKFingerprint
    }
};

const GoldenBinding* BindingForProfileId(uint64_t profile_id)
{
    for (const GoldenBinding& binding : kGoldenBindings) {
        if (binding.ProfileId == profile_id) {
            return &binding;
        }
    }
    return nullptr;
}

void PrintProgress(void* context, uint32_t block_count)
{
    const uint32_t max_block_count = *(const uint32_t*)context;
    if (block_count % 8192u == 0u || block_count == max_block_count)
    {
        std::fprintf(stderr, "  K=%u/%u\n", block_count, max_block_count);
        std::fflush(stderr);
    }
}

bool CheckContractTable()
{
    uint32_t contract_count = 0u;
    const wirehair_v2::EquationFingerprintContract* contracts =
        wirehair_v2::EquationFingerprintContracts(contract_count);
    if (!contracts || contract_count !=
        (uint32_t)(sizeof(kGoldenBindings) / sizeof(kGoldenBindings[0])))
    {
        std::fprintf(stderr,
            "fingerprint: contract table has %u entries, goldens have %u\n",
            contract_count,
            (uint32_t)(sizeof(kGoldenBindings) / sizeof(kGoldenBindings[0])));
        return false;
    }
    for (uint32_t i = 0; i < contract_count; ++i)
    {
        if (!BindingForProfileId(contracts[i].ProfileId))
        {
            std::fprintf(stderr,
                "fingerprint: contract %s has no golden binding\n",
                contracts[i].Name);
            return false;
        }
    }
    if (WIREHAIR_V2_PROFILE_CURRENT != WIREHAIR_V2_PROFILE_CERTIFIED_2026_07)
    {
        std::fprintf(stderr,
            "fingerprint: WIREHAIR_V2_PROFILE_CURRENT no longer aliases the "
            "certified profile\n");
        return false;
    }
    return true;
}

int PrintGoldens(uint32_t max_block_count)
{
    uint32_t contract_count = 0u;
    const wirehair_v2::EquationFingerprintContract* contracts =
        wirehair_v2::EquationFingerprintContracts(contract_count);
    std::string lines;
    for (uint32_t i = 0; i < contract_count; ++i)
    {
        std::fprintf(stderr, "computing %s fingerprint (K=%u..%u)\n",
            contracts[i].Name,
            wirehair_v2::kEquationFingerprintMinBlockCount,
            max_block_count);
        uint8_t digest[wirehair_v2::kEquationFingerprintBytes];
        if (!wirehair_v2::ComputeEquationFingerprint(
                contracts[i],
                wirehair_v2::kEquationFingerprintMinBlockCount,
                max_block_count,
                digest,
                PrintProgress,
                &max_block_count))
        {
            std::fprintf(stderr,
                "fingerprint: computation failed for %s\n",
                contracts[i].Name);
            return 1;
        }
        char hex[wirehair_v2::kEquationFingerprintHexChars + 1u];
        wirehair_v2::FormatEquationFingerprintHex(digest, hex);
        const GoldenBinding* binding =
            BindingForProfileId(contracts[i].ProfileId);
        lines += "static const char ";
        lines += binding->ConstantName;
        lines += "[] =\n    \"";
        lines += hex;
        lines += "\";\n";
    }

    if (max_block_count != wirehair_v2::kEquationFingerprintMaxBlockCount)
    {
        std::printf(
            "// TRUNCATED RANGE (max K %u) -- NOT FROZEN GOLDENS\n%s",
            max_block_count, lines.c_str());
        return 0;
    }
    std::printf(
        "// --- BEGIN FROZEN V2 EQUATION FINGERPRINTS ---\n"
        "// Fill by building, then running:\n"
        "//   wirehair_v2_fingerprint_test --print-goldens\n"
        "// and pasting the printed block over this one.\n"
        "%s"
        "// --- END FROZEN V2 EQUATION FINGERPRINTS ---\n",
        lines.c_str());
    return 0;
}

int CheckGoldens()
{
    bool unset = false;
    for (const GoldenBinding& binding : kGoldenBindings)
    {
        if (std::strcmp(binding.Golden, "UNSET") == 0)
        {
            std::fprintf(stderr,
                "fingerprint: %s is UNSET -- build, run "
                "'wirehair_v2_fingerprint_test --print-goldens', and paste "
                "the printed block into codec/V2FingerprintTest.cpp\n",
                binding.ConstantName);
            unset = true;
        }
    }
    if (unset) {
        return 1;
    }

    uint32_t contract_count = 0u;
    const wirehair_v2::EquationFingerprintContract* contracts =
        wirehair_v2::EquationFingerprintContracts(contract_count);
    uint32_t max_block_count =
        wirehair_v2::kEquationFingerprintMaxBlockCount;
    bool drift = false;
    for (uint32_t i = 0; i < contract_count; ++i)
    {
        std::fprintf(stderr, "checking %s fingerprint (K=%u..%u)\n",
            contracts[i].Name,
            wirehair_v2::kEquationFingerprintMinBlockCount,
            max_block_count);
        uint8_t digest[wirehair_v2::kEquationFingerprintBytes];
        if (!wirehair_v2::ComputeEquationFingerprint(
                contracts[i],
                wirehair_v2::kEquationFingerprintMinBlockCount,
                max_block_count,
                digest,
                PrintProgress,
                &max_block_count))
        {
            std::fprintf(stderr,
                "fingerprint: computation failed for %s\n",
                contracts[i].Name);
            return 1;
        }
        char hex[wirehair_v2::kEquationFingerprintHexChars + 1u];
        wirehair_v2::FormatEquationFingerprintHex(digest, hex);
        const GoldenBinding* binding =
            BindingForProfileId(contracts[i].ProfileId);
        if (std::strcmp(hex, binding->Golden) != 0)
        {
            std::fprintf(stderr,
                "fingerprint DRIFT for %s (profile id %016llx):\n"
                "  frozen:   %s\n"
                "  computed: %s\n"
                "An equation-affecting input changed under a frozen public "
                "profile ID.  Ship the change under a NEW profile ID; do not "
                "update this golden.\n",
                contracts[i].Name,
                (unsigned long long)contracts[i].ProfileId,
                binding->Golden,
                hex);
            drift = true;
        }
    }
    if (drift) {
        return 1;
    }
    std::puts("V2 all-K equation fingerprints: PASS");
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    bool print_goldens = false;
    uint32_t max_block_count =
        wirehair_v2::kEquationFingerprintMaxBlockCount;
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--print-goldens") == 0) {
            print_goldens = true;
        }
        else if (std::strcmp(argv[i], "--max-k") == 0 && i + 1 < argc) {
            const unsigned long parsed = std::strtoul(argv[++i], nullptr, 10);
            if (parsed < wirehair_v2::kEquationFingerprintMinBlockCount ||
                parsed > wirehair_v2::kEquationFingerprintMaxBlockCount)
            {
                std::fprintf(stderr, "invalid --max-k value\n");
                return 2;
            }
            max_block_count = (uint32_t)parsed;
        }
        else {
            std::fprintf(stderr,
                "usage: %s [--print-goldens [--max-k K]]\n", argv[0]);
            return 2;
        }
    }
    if (!print_goldens &&
        max_block_count != wirehair_v2::kEquationFingerprintMaxBlockCount)
    {
        // Frozen comparisons are only defined over the complete K domain.
        std::fprintf(stderr, "--max-k requires --print-goldens\n");
        return 2;
    }

    if (!CheckContractTable()) {
        return 1;
    }
    return print_goldens ? PrintGoldens(max_block_count) : CheckGoldens();
}
