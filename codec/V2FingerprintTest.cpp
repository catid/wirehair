#include "WirehairV2Fingerprint.h"

#include <wirehair/wirehair.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

/*
    All-K equation-freeze test for public serialized V2 profiles and stable
    benchmark-only architecture targets.

    For every supported K this digests the complete equation-affecting
    expansion of each versioned contract (see WirehairV2Fingerprint.h for the
    exact stream) and compares one SHA-256 per contract against the frozen
    constants below.  Public-profile drift is a compatibility bug; stable
    benchmark-target drift invalidates experiment identity.  Either change
    requires a new contract ID instead of silently replacing a golden.

    The constants are filled from a build of the frozen implementation:

        wirehair_v2_fingerprint_test --print-goldens

    prints the complete constant block ready to paste over the marked
    section.  Unset constants FAIL the test; they never skip.
*/

namespace {

bool ParseMaxBlockCount(const char* text, uint32_t& value_out)
{
    if (!text || text[0] < '0' || text[0] > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' ||
        parsed < wirehair_v2::kEquationFingerprintMinBlockCount ||
        parsed > wirehair_v2::kEquationFingerprintMaxBlockCount)
    {
        return false;
    }
    value_out = (uint32_t)parsed;
    return true;
}

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
static const char kMixedMix2TwoAnchorAllKFingerprint[] =
    "272b24d707befdeb08afdba83472ccdaeed9ebfac17b05bc6c251793b6b17a7b";
static const char kDispatchV1AllKFingerprint[] =
    "66902904804e0efb31971855da6d36dc74e06ed1653b28d6d993d271817a69b1";
// --- END FROZEN V2 EQUATION FINGERPRINTS ---

// Fixed experimental seam: mixed K=4096, SharedCauchyX geometry, one
// independently scheduled GF(256) row.  This pins every grouped-stream field
// and its first-heavy-column boundary, rather than merely requiring that the
// digest differ from the ungrouped stream.
static const char kGroupedGF256K4096Fingerprint[] =
    "4c2af25b58d403d6df117a1bb872d25a27b4225ba2c927a57821511b8c42c7be";

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
    },
    {
        WIREHAIR_V2_PROFILE_MIXED_MIX2_TWO_ANCHOR_2026_07,
        "kMixedMix2TwoAnchorAllKFingerprint",
        kMixedMix2TwoAnchorAllKFingerprint
    },
    {
        wirehair_v2::kWh2DispatchV1ContractId,
        "kDispatchV1AllKFingerprint",
        kDispatchV1AllKFingerprint
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
    uint32_t precode_public_count = 0u;
    for (uint32_t i = 0; i < contract_count; ++i)
    {
        const wirehair_v2::EquationFingerprintContract& contract =
            contracts[i];
        if (!contract.Name || !contract.Name[0] ||
            !contract.CanonicalName || !contract.CanonicalName[0] ||
            contract.Kind != wirehair_v2::V2EquationKind::Precode ||
            !BindingForProfileId(contract.ProfileId))
        {
            std::fprintf(stderr,
                "fingerprint: contract %s has no golden binding\n",
                contract.Name ? contract.Name : "(null)");
            return false;
        }
        if (contract.PublicProfile) {
            if (contract.SeedPolicy !=
                    wirehair_v2::V2SeedDerivation::ProfileDerived ||
                contract.SeedAttemptCount !=
                    wirehair_v2::kMaxPacketSeedAttempts)
            {
                std::fprintf(stderr,
                    "fingerprint: public contract %s changed seed policy\n",
                    contract.Name);
                return false;
            }
            ++precode_public_count;
        }
        uint8_t digest[wirehair_v2::kEquationFingerprintBytes];
        if (!wirehair_v2::ComputeEquationContractNameDigest(
                contract, digest))
        {
            std::fprintf(stderr,
                "fingerprint: contract %s canonical digest failed\n",
                contract.Name);
            return false;
        }
        uint64_t digest_id = 0u;
        for (uint32_t byte = 0u; byte < 8u; ++byte) {
            digest_id = (digest_id << 8u) | digest[byte];
        }
        if (digest_id != contract.ProfileId)
        {
            std::fprintf(stderr,
                "fingerprint: contract %s ID %016llx does not match "
                "canonical SHA-256 prefix %016llx\n",
                contract.Name,
                (unsigned long long)contract.ProfileId,
                (unsigned long long)digest_id);
            return false;
        }
        for (uint32_t prior = 0u; prior < i; ++prior)
        {
            if (contracts[prior].ProfileId == contract.ProfileId ||
                std::strcmp(contracts[prior].Name, contract.Name) == 0 ||
                std::strcmp(
                    contracts[prior].CanonicalName,
                    contract.CanonicalName) == 0)
            {
                std::fprintf(stderr,
                    "fingerprint: duplicate contract registry entry %s\n",
                    contract.Name);
                return false;
            }
        }
    }
    if (precode_public_count != 4u)
    {
        std::fprintf(stderr,
            "fingerprint: contract table has %u public precode profiles, "
            "expected 4\n",
            precode_public_count);
        return false;
    }
    const wirehair_v2::EquationFingerprintContract* tiny =
        wirehair_v2::FindV2EquationContract("tiny_mds_2026_07");
    uint8_t tiny_digest[wirehair_v2::kEquationFingerprintBytes] = {};
    if (!tiny ||
        tiny != wirehair_v2::FindV2EquationContract(
            WIREHAIR_V2_PROFILE_TINY_MDS_2026_07, true) ||
        tiny->ProfileId != WIREHAIR_V2_PROFILE_TINY_MDS_2026_07 ||
        tiny->Kind != wirehair_v2::V2EquationKind::TinyMds ||
        !tiny->PublicProfile ||
        tiny->SeedAttemptCount != 1u ||
        tiny->RecoveryMixCount != 0u ||
        !wirehair_v2::ComputeEquationContractNameDigest(
            *tiny, tiny_digest))
    {
        std::fprintf(stderr,
            "fingerprint: tiny-MDS contract registry mismatch\n");
        return false;
    }
    uint64_t tiny_digest_id = 0u;
    for (uint32_t byte = 0u; byte < 8u; ++byte) {
        tiny_digest_id = (tiny_digest_id << 8u) | tiny_digest[byte];
    }
    if (tiny_digest_id != tiny->ProfileId)
    {
        std::fprintf(stderr,
            "fingerprint: tiny-MDS ID %016llx does not match canonical "
            "SHA-256 prefix %016llx\n",
            (unsigned long long)tiny->ProfileId,
            (unsigned long long)tiny_digest_id);
        return false;
    }
    uint8_t rejected_equation_digest[
        wirehair_v2::kEquationFingerprintBytes];
    std::memset(
        rejected_equation_digest, 0xa5,
        sizeof(rejected_equation_digest));
    uint8_t rejected_equation_before[
        wirehair_v2::kEquationFingerprintBytes];
    std::memcpy(
        rejected_equation_before,
        rejected_equation_digest,
        sizeof(rejected_equation_before));
    if (wirehair_v2::ComputeEquationFingerprint(
            *tiny, 2u, 2u, rejected_equation_digest) ||
        std::memcmp(
            rejected_equation_digest,
            rejected_equation_before,
            sizeof(rejected_equation_digest)) != 0)
    {
        std::fprintf(stderr,
            "fingerprint: tiny MDS entered the precode all-K stream\n");
        return false;
    }
    for (uint32_t i = 0u; i < contract_count; ++i)
    {
        if (contracts[i].ProfileId == tiny->ProfileId ||
            std::strcmp(contracts[i].Name, tiny->Name) == 0 ||
            std::strcmp(
                contracts[i].CanonicalName, tiny->CanonicalName) == 0)
        {
            std::fprintf(stderr,
                "fingerprint: tiny-MDS contract collides with %s\n",
                contracts[i].Name);
            return false;
        }
    }
    const wirehair_v2::EquationFingerprintContract* dispatch =
        wirehair_v2::FindV2EquationContract("dispatch-v1");
    if (!dispatch ||
        dispatch->ProfileId != wirehair_v2::kWh2DispatchV1ContractId ||
        dispatch->PublicProfile ||
        dispatch->SeedPolicy != wirehair_v2::V2SeedDerivation::RawUniform ||
        dispatch->SeedAttemptCount != 1u)
    {
        std::fprintf(stderr,
            "fingerprint: dispatch-v1 must be the one internal target\n");
        return false;
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

bool CheckRawDispatchSeams()
{
    struct RawFoldGolden
    {
        uint64_t ConstructionSeed;
        uint32_t PacketSeed;
    };
    const RawFoldGolden fold_goldens[] = {
        { UINT64_C(0), UINT32_C(0) },
        { UINT64_C(1), UINT32_C(1) },
        { UINT64_C(0x0000000100000000), UINT32_C(1) },
        { UINT64_C(0xdeadbeef01234567), UINT32_C(0xdf8efb88) }
    };
    for (const RawFoldGolden& golden : fold_goldens)
    {
        if (wirehair_v2::RawUniformPacketPeelSeed(
                golden.ConstructionSeed) != golden.PacketSeed)
        {
            std::fprintf(stderr,
                "fingerprint: raw packet-seed fold changed at 0x%llx\n",
                (unsigned long long)golden.ConstructionSeed);
            return false;
        }
    }

    const wirehair_v2::EquationFingerprintContract* dispatch =
        wirehair_v2::FindV2EquationContract("dispatch-v1");
    if (!dispatch) return false;

    const uint32_t block_counts[] = {
        2u, 5u, 15u, 23u, 31u, 63u, 64u, 77u, 92u, 100u, 101u,
        4096u, 64000u
    };
    const uint64_t construction_seed =
        wirehair_v2::kEquationFingerprintRawConstructionSeed;
    const uint32_t packet_seed =
        wirehair_v2::RawUniformPacketPeelSeed(construction_seed);
    const wirehair_v2::MessagePrecodeEncoderOptions options =
        wirehair_v2::MessageOptionsForContract(*dispatch);

    for (uint32_t K : block_counts)
    {
        wirehair_v2::SeedProfile raw;
        if (!wirehair_v2::MakeRawContractProfile(
                *dispatch, K,
                wirehair_v2::kEquationFingerprintCanonicalBlockBytes,
                construction_seed, raw))
        {
            std::fprintf(stderr,
                "fingerprint: raw dispatch expansion failed at K=%u\n", K);
            return false;
        }
        wirehair_v2::PrecodeParams params;
        wirehair_v2::PacketRowConfig packet;
        if (raw.BlockCount != K ||
            raw.BlockBytes !=
                wirehair_v2::kEquationFingerprintCanonicalBlockBytes ||
            raw.DenseCount != 0u || raw.PeelSeed != 0u ||
            raw.DenseSeed != 0u || raw.PeelSeedBucket != 0u ||
            raw.UsedPeelFixup || raw.UsedDenseFixup || raw.Tuned ||
            !raw.V2SeedSelected || raw.V2SeedAttempt != 0u ||
            raw.V2SeedPolicy != wirehair_v2::V2SeedDerivation::RawUniform ||
            raw.V2Architecture !=
                wirehair_v2::V2PrecodeArchitecture::SmallBandD4 ||
            raw.V2StaircaseCount !=
                wirehair_v2::SmallBandStaircaseCount(K) ||
            raw.V2DenseRowCount != 4u || raw.V2HeavyRowCount != 12u ||
            raw.V2CompletionField !=
                wirehair_v2::CompletionField::MixedGF256GF16 ||
            raw.V2PrecodeSeed != construction_seed ||
            raw.V2PacketPeelSeed != packet_seed ||
            !wirehair_v2::ResolveMessagePrecodeConfiguration(
                raw, options, params, packet) ||
            params.Staircase != raw.V2StaircaseCount ||
            params.DenseRows != raw.V2DenseRowCount ||
            params.Seed != construction_seed ||
            packet.PeelSeed != packet_seed ||
            packet.MixCount != wirehair_v2::kCertifiedPacketMixCount)
        {
            std::fprintf(stderr,
                "fingerprint: raw dispatch contract mismatch at K=%u\n", K);
            return false;
        }
    }

    wirehair_v2::SeedProfile rejected;
    const wirehair_v2::EquationFingerprintContract* public_contract =
        wirehair_v2::FindV2EquationContract(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07);
    if (!public_contract ||
        wirehair_v2::MakeRawContractProfile(
            *public_contract, 64u, 1280u, construction_seed, rejected) ||
        wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u, 1279u, construction_seed, rejected))
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch accepted a fallback/public profile\n");
        return false;
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const bool configured =
        wirehair_v2::SetMixedGF16RowsForTesting(3u);
    const bool accepted_noncanonical = configured &&
        wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u, 1280u, construction_seed, rejected);
    const bool restored = wirehair_v2::SetMixedGF16RowsForTesting(
        wirehair_v2::kMixedGF16Rows);
    if (!configured || !restored || accepted_noncanonical)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch accepted noncanonical mixed rows\n");
        return false;
    }
    const wirehair_v2::EquationFingerprintContract* mixed =
        wirehair_v2::FindV2EquationContract(
            WIREHAIR_V2_PROFILE_MIXED_2026_07);
    uint8_t canonical[wirehair_v2::kEquationFingerprintBytes] = {};
    uint8_t inert_seed[wirehair_v2::kEquationFingerprintBytes] = {};
    const bool hashed_canonical = mixed &&
        wirehair_v2::ComputeEquationFingerprint(
            *mixed, 2u, 2u, canonical);
    wirehair_v2::SetMixedResidueHashSeedForTesting(UINT32_C(0x12345678));
    const bool accepted_inert_hash_seed =
        wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u, 1280u, construction_seed, rejected);
    const bool hashed_inert = mixed &&
        wirehair_v2::ComputeEquationFingerprint(
            *mixed, 2u, 2u, inert_seed);
    wirehair_v2::SetMixedResidueHashSeedForTesting(0u);
    if (!hashed_canonical || !accepted_inert_hash_seed || !hashed_inert ||
        std::memcmp(canonical, inert_seed, sizeof(canonical)) != 0)
    {
        std::fprintf(stderr,
            "fingerprint: inactive residue hash seed changed production "
            "state\n");
        return false;
    }

    // A raw profile may outlive the thread-local experiment configuration
    // under which it was bound.  Revalidation must reject it if active
    // equation controls later cease to match the frozen dispatch contract.
    wirehair_v2::SeedProfile stored_raw;
    const bool made_stored_raw =
        wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u, 1280u, construction_seed, stored_raw);
    const bool shared_geometry =
        wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
    const bool skewed = shared_geometry &&
        wirehair_v2::SetMixedResidueSkewForTesting(1u);
    wirehair_v2::PrecodeParams changed_params;
    wirehair_v2::PacketRowConfig changed_packet;
    const bool replayed_noncanonical = skewed &&
        wirehair_v2::ResolveMessagePrecodeConfiguration(
            stored_raw, options, changed_params, changed_packet);
    const bool rebound_noncanonical = skewed &&
        wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u, 1280u, construction_seed, rejected);
    const bool restored_skew =
        wirehair_v2::SetMixedResidueSkewForTesting(0u);
    const bool restored_geometry =
        wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::FrozenPowerX);
    const bool replayed_restored = made_stored_raw && restored_skew &&
        restored_geometry &&
        wirehair_v2::ResolveMessagePrecodeConfiguration(
            stored_raw, options, changed_params, changed_packet);
    if (!made_stored_raw || !shared_geometry || !skewed ||
        replayed_noncanonical || rebound_noncanonical ||
        !restored_skew || !restored_geometry || !replayed_restored)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch replay ignored active mixed state\n");
        return false;
    }

    // Exercise every other unreceipted equation control through BOTH stable
    // target entry points.  Always restore a hook before checking its result,
    // so one failed assertion cannot poison the later golden sweep.
    const auto raw_paths_accept = [&]() {
        wirehair_v2::SeedProfile rebound;
        wirehair_v2::PrecodeParams replay_params;
        wirehair_v2::PacketRowConfig replay_packet;
        const bool canonical =
            wirehair_v2::IsCanonicalStableTargetEquationState();
        const bool rebound_ok = wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u, 1280u, construction_seed, rebound);
        const bool replay_ok =
            wirehair_v2::ResolveMessagePrecodeConfiguration(
                stored_raw, options, replay_params, replay_packet);
        return canonical && rebound_ok && replay_ok;
    };
    const auto raw_paths_reject = [&]() {
        wirehair_v2::SeedProfile rebound;
        wirehair_v2::PrecodeParams replay_params;
        wirehair_v2::PacketRowConfig replay_packet;
        const bool canonical =
            wirehair_v2::IsCanonicalStableTargetEquationState();
        const bool rebound_ok = wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u, 1280u, construction_seed, rebound);
        const bool replay_ok =
            wirehair_v2::ResolveMessagePrecodeConfiguration(
                stored_raw, options, replay_params, replay_packet);
        return !canonical && !rebound_ok && !replay_ok;
    };

    const bool set_multiplier =
        wirehair_v2::SetPacketRowSeedMultiplierForTesting(3u);
    const bool rejected_multiplier =
        set_multiplier && raw_paths_reject();
    const bool reset_multiplier =
        wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u);
    const bool restored_multiplier =
        reset_multiplier && raw_paths_accept();
    if (!set_multiplier || !rejected_multiplier ||
        !reset_multiplier || !restored_multiplier)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch ignored packet-id multiplier state\n");
        return false;
    }

    wirehair_v2::SetPacketRowSeedAvalancheForTesting(true);
    const bool rejected_avalanche = raw_paths_reject();
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(false);
    const bool restored_avalanche = raw_paths_accept();
    if (!rejected_avalanche || !restored_avalanche)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch ignored packet-id avalanche state\n");
        return false;
    }

    wirehair_v2::SetOddPacketPeelSeedXorForTesting(UINT32_C(0x5a17c3e9));
    const bool rejected_odd_seed = raw_paths_reject();
    wirehair_v2::SetOddPacketPeelSeedXorForTesting(0u);
    const bool restored_odd_seed = raw_paths_accept();
    if (!rejected_odd_seed || !restored_odd_seed)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch ignored odd packet peel-seed state\n");
        return false;
    }

    wirehair_v2::SetStaircaseDegreesForTesting(
        std::vector<double>{ 1.0 });
    const bool rejected_staircase_shape = raw_paths_reject();
    wirehair_v2::ClearStaircaseDegreesForTesting();
    const bool restored_staircase_shape = raw_paths_accept();
    if (!rejected_staircase_shape || !restored_staircase_shape)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch ignored staircase shape state\n");
        return false;
    }

    // K=64 has S=10 and the stock two-hit budget has 128 edges, so this is a
    // valid realized degree sequence rather than merely an invalid hook value.
    wirehair_v2::SetStaircaseRowDegreesForTesting(
        std::vector<uint32_t>{
            13u, 13u, 13u, 13u, 13u,
            13u, 13u, 13u, 12u, 12u
        });
    const bool rejected_staircase_rows = raw_paths_reject();
    wirehair_v2::ClearStaircaseRowDegreesForTesting();
    const bool restored_staircase_rows = raw_paths_accept();
    if (!rejected_staircase_rows || !restored_staircase_rows)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch ignored pinned staircase rows\n");
        return false;
    }

    // These two overlays are named in the target receipt.  Keep them
    // deliberately outside the canonical-state gate and prove they remain
    // usable through both binding and replay.
    const bool set_receipted_pmf =
        wirehair_v2::SetPeelDegreesForTesting(
            std::vector<double>{ 1.0, 1.0 });
    const bool accepted_receipted_pmf =
        set_receipted_pmf && raw_paths_accept();
    wirehair_v2::ClearPeelDegreesForTesting();
    const bool restored_after_pmf = raw_paths_accept();
    if (!set_receipted_pmf || !accepted_receipted_pmf ||
        !restored_after_pmf)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch rejected receipted peel PMF state\n");
        return false;
    }

    wirehair_v2::SetStaircaseDegreeScaleForTesting(0.0);
    const bool accepted_receipted_scale = raw_paths_accept();
    wirehair_v2::ClearStaircaseDegreeScaleForTesting();
    const bool restored_after_scale = raw_paths_accept();
    if (!accepted_receipted_scale || !restored_after_scale)
    {
        std::fprintf(stderr,
            "fingerprint: raw dispatch rejected receipted staircase scale\n");
        return false;
    }
#endif
    return true;
}

bool CheckGroupedCoefficientFingerprint()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const wirehair_v2::EquationFingerprintContract* mixed =
        wirehair_v2::FindV2EquationContract(
            WIREHAIR_V2_PROFILE_MIXED_2026_07);
    const uint32_t test_block_count = 4096u;
    uint8_t ungrouped[wirehair_v2::kEquationFingerprintBytes] = {};
    uint8_t grouped[wirehair_v2::kEquationFingerprintBytes] = {};

    const bool shared_geometry =
        wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX);
    const bool hashed_ungrouped = shared_geometry &&
        mixed &&
        wirehair_v2::ComputeEquationFingerprint(
            *mixed, test_block_count, test_block_count, ungrouped);
    const bool configured_group =
        hashed_ungrouped &&
        wirehair_v2::SetMixedGroupedGF256RowsForTesting(1u);
    bool boundary_ok = false;
    if (configured_group)
    {
        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile(
                test_block_count,
                wirehair_v2::kEquationFingerprintCanonicalBlockBytes);
        const wirehair_v2::MessagePrecodeEncoderOptions options =
            wirehair_v2::MessageOptionsForContract(*mixed);
        wirehair_v2::PrecodeParams params;
        wirehair_v2::PacketRowConfig packet;
        if (wirehair_v2::ResolveMessagePrecodeConfiguration(
                base, options, params, packet))
        {
            const uint32_t first_heavy_column =
                params.BlockCount + params.Staircase + params.DenseRows;
            bool saw_secondary_schedule = false;
            for (uint32_t column = 0u;
                 column < first_heavy_column; ++column)
            {
                if (wirehair_v2::
                        ActiveMixedGroupedGF256CoefficientResidue(
                            column, first_heavy_column) !=
                    wirehair_v2::ActiveMixedCoefficientResidue(column))
                {
                    saw_secondary_schedule = true;
                    break;
                }
            }
            boundary_ok = first_heavy_column != 0u &&
                saw_secondary_schedule &&
                wirehair_v2::ActiveMixedGroupedGF256CoefficientResidue(
                    first_heavy_column, first_heavy_column) ==
                    wirehair_v2::ActiveMixedCoefficientResidue(
                        first_heavy_column) &&
                wirehair_v2::ActiveMixedGroupedGF256CoefficientResidue(
                    first_heavy_column + 1u, first_heavy_column) ==
                    wirehair_v2::ActiveMixedCoefficientResidue(
                        first_heavy_column + 1u) &&
                wirehair_v2::ActiveMixedGroupedGF256CoefficientResidue(
                    first_heavy_column + params.HeavyRows - 1u,
                    first_heavy_column) ==
                    wirehair_v2::ActiveMixedCoefficientResidue(
                        first_heavy_column + params.HeavyRows - 1u);
        }
    }
    const bool hashed_grouped = configured_group &&
        wirehair_v2::ComputeEquationFingerprint(
            *mixed, test_block_count, test_block_count, grouped);
    char grouped_hex[wirehair_v2::kEquationFingerprintHexChars + 1u] = {};
    if (hashed_grouped) {
        wirehair_v2::FormatEquationFingerprintHex(grouped, grouped_hex);
    }

    // Restore the thread-local hooks before reporting any failure so later
    // golden checks always observe the production coefficient geometry.
    const bool cleared_group =
        wirehair_v2::SetMixedGroupedGF256RowsForTesting(0u);
    const bool restored_geometry =
        wirehair_v2::SetMixedCoefficientGeometryForTesting(
            wirehair_v2::MixedCoefficientGeometry::FrozenPowerX);
    const wirehair_v2::EquationFingerprintContract* dispatch =
        wirehair_v2::FindV2EquationContract("dispatch-v1");
    wirehair_v2::SeedProfile restored_raw;
    const bool raw_construction_restored =
        cleared_group && restored_geometry && dispatch &&
        wirehair_v2::MakeRawContractProfile(
            *dispatch, 64u,
            wirehair_v2::kEquationFingerprintCanonicalBlockBytes,
            wirehair_v2::kEquationFingerprintRawConstructionSeed,
            restored_raw);
    if (!shared_geometry || !mixed || !hashed_ungrouped ||
        !configured_group || !boundary_ok || !hashed_grouped ||
        !cleared_group || !restored_geometry ||
        !raw_construction_restored ||
        std::memcmp(ungrouped, grouped, sizeof(ungrouped)) == 0 ||
        std::strcmp(
            grouped_hex, kGroupedGF256K4096Fingerprint) != 0)
    {
        std::fprintf(stderr,
            "fingerprint: grouped GF(256) schedule was not fingerprinted\n");
        return false;
    }
#endif
    return true;
}

int PrintGoldens(uint32_t max_block_count, const char* profile_name)
{
    uint32_t contract_count = 0u;
    const wirehair_v2::EquationFingerprintContract* contracts =
        wirehair_v2::EquationFingerprintContracts(contract_count);
    std::string lines;
    for (uint32_t i = 0; i < contract_count; ++i)
    {
        if (profile_name &&
            std::strcmp(contracts[i].Name, profile_name) != 0)
        {
            continue;
        }
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
                "fingerprint DRIFT for %s (contract id %016llx):\n"
                "  frozen:   %s\n"
                "  computed: %s\n"
                "An equation-affecting input changed under a frozen "
                "contract ID.  Ship the change under a NEW contract ID; do "
                "not update this golden.\n",
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
    const char* profile_name = nullptr;
    uint32_t max_block_count =
        wirehair_v2::kEquationFingerprintMaxBlockCount;
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--print-goldens") == 0) {
            print_goldens = true;
        }
        else if (std::strcmp(argv[i], "--max-k") == 0 && i + 1 < argc) {
            if (!ParseMaxBlockCount(argv[++i], max_block_count))
            {
                std::fprintf(stderr, "invalid --max-k value\n");
                return 2;
            }
        }
        else if (std::strcmp(argv[i], "--profile") == 0 && i + 1 < argc) {
            profile_name = argv[++i];
        }
        else {
            std::fprintf(stderr,
                "usage: %s [--print-goldens [--max-k K] "
                "[--profile NAME]]\n", argv[0]);
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
    if (profile_name)
    {
        uint32_t contract_count = 0u;
        const wirehair_v2::EquationFingerprintContract* contracts =
            wirehair_v2::EquationFingerprintContracts(contract_count);
        bool found = false;
        for (uint32_t i = 0u; i < contract_count; ++i) {
            if (std::strcmp(contracts[i].Name, profile_name) == 0) {
                found = true;
            }
        }
        if (!print_goldens || !found) {
            std::fprintf(stderr,
                "--profile requires --print-goldens and a known name\n");
            return 2;
        }
    }

    if (!CheckContractTable() || !CheckRawDispatchSeams() ||
        !CheckGroupedCoefficientFingerprint())
    {
        return 1;
    }
    return print_goldens ?
        PrintGoldens(max_block_count, profile_name) : CheckGoldens();
}
