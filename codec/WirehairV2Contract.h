#pragma once

#include "WirehairV2PrecodeEncode.h"

#include <wirehair/wirehair.h>

#include <stdint.h>

namespace wirehair_v2 {

/**
    Complete equation configurations used by serialized public profiles and
    benchmark-only architecture targets.

    The dispatch target has a stable identifier so every experiment can bind
    exact equations without prematurely adding it to the public serialized
    API.  Promotion to a WIREHAIR_V2_PROFILE_* identifier is a separate gate.
*/
struct V2EquationContract
{
    uint64_t ProfileId;
    const char* Name;
    const char* CanonicalName;
    CompletionField Completion;
    uint32_t RecoveryMixCount;
    bool AdaptiveDenseTwoAnchor;
    V2PrecodeArchitecture Architecture;
    V2SeedDerivation SeedPolicy;
    uint32_t SeedAttemptCount;
    bool PublicProfile;
};

static const uint64_t kWh2DispatchV1ContractId =
    UINT64_C(0xa98c37c23ee7feae);
static const char kWh2DispatchV1CanonicalName[] =
    "wirehair:v2:precode-v5-mixed10gf256-2gf65536-even-smallband100-d4:"
    "packet-v4-mix3:seed-raw-uniform-v1-xorfold32-attempt0-no-tables-"
    "no-fixups:dispatch-v1-2026-07";

inline const V2EquationContract* V2EquationContracts(uint32_t& count_out)
{
    static const V2EquationContract kContracts[] = {
        {
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            "certified_2026_07",
            "wirehair:v2:precode-v2:packet-v4:certified-2026-07",
            CompletionField::GF256,
            kCertifiedPacketMixCount,
            false,
            V2PrecodeArchitecture::LegacyD12,
            V2SeedDerivation::ProfileDerived,
            kMaxPacketSeedAttempts,
            true
        },
        {
            WIREHAIR_V2_PROFILE_MIXED_2026_07,
            "mixed_2026_07",
            "wirehair:v2:precode-v3-mixed10gf256-2gf65536-even:"
                "packet-v4:certified-2026-07",
            CompletionField::MixedGF256GF16,
            kCertifiedPacketMixCount,
            false,
            V2PrecodeArchitecture::LegacyD12,
            V2SeedDerivation::ProfileDerived,
            kMaxPacketSeedAttempts,
            true
        },
        {
            WIREHAIR_V2_PROFILE_MIXED_MIX2_2026_07,
            "mixed_mix2_2026_07",
            "wirehair:v2:precode-v3-mixed10gf256-2gf65536-even:"
                "packet-v4-mix2:certified-2026-07",
            CompletionField::MixedGF256GF16,
            2u,
            false,
            V2PrecodeArchitecture::LegacyD12,
            V2SeedDerivation::ProfileDerived,
            kMaxPacketSeedAttempts,
            true
        },
        {
            WIREHAIR_V2_PROFILE_MIXED_MIX2_TWO_ANCHOR_2026_07,
            "mixed_mix2_two_anchor_2026_07",
            "wirehair:v2:precode-v4-mixed10gf256-2gf65536-even-"
                "d12-two-anchor-k4096:packet-v4-mix2:certified-2026-07",
            CompletionField::MixedGF256GF16,
            2u,
            true,
            V2PrecodeArchitecture::LegacyD12,
            V2SeedDerivation::ProfileDerived,
            kMaxPacketSeedAttempts,
            true
        },
        {
            kWh2DispatchV1ContractId,
            "dispatch-v1",
            kWh2DispatchV1CanonicalName,
            CompletionField::MixedGF256GF16,
            kCertifiedPacketMixCount,
            false,
            V2PrecodeArchitecture::SmallBandD4,
            V2SeedDerivation::RawUniform,
            1u,
            false
        }
    };
    count_out = (uint32_t)(sizeof(kContracts) / sizeof(kContracts[0]));
    return kContracts;
}

inline const V2EquationContract* FindV2EquationContract(
    uint64_t id,
    bool public_only = false)
{
    uint32_t count = 0u;
    const V2EquationContract* contracts = V2EquationContracts(count);
    for (uint32_t i = 0u; i < count; ++i) {
        if (contracts[i].ProfileId == id &&
            (!public_only || contracts[i].PublicProfile))
        {
            return &contracts[i];
        }
    }
    return nullptr;
}

inline const V2EquationContract* FindV2EquationContract(
    const char* name)
{
    if (!name) return nullptr;
    uint32_t count = 0u;
    const V2EquationContract* contracts = V2EquationContracts(count);
    for (uint32_t i = 0u; i < count; ++i)
    {
        const char* left = contracts[i].Name;
        const char* right = name;
        while (*left && *left == *right) {
            ++left;
            ++right;
        }
        if (*left == '\0' && *right == '\0') return &contracts[i];
    }
    return nullptr;
}

inline const V2EquationContract* FindPublicV2EquationContract(
    CompletionField completion,
    uint32_t recovery_mix_count,
    bool adaptive_dense_two_anchor,
    V2PrecodeArchitecture architecture)
{
    uint32_t count = 0u;
    const V2EquationContract* contracts = V2EquationContracts(count);
    for (uint32_t i = 0u; i < count; ++i)
    {
        const V2EquationContract& contract = contracts[i];
        if (contract.PublicProfile &&
            contract.SeedPolicy == V2SeedDerivation::ProfileDerived &&
            contract.SeedAttemptCount == kMaxPacketSeedAttempts &&
            contract.Completion == completion &&
            contract.RecoveryMixCount == recovery_mix_count &&
            contract.AdaptiveDenseTwoAnchor ==
                adaptive_dense_two_anchor &&
            contract.Architecture == architecture)
        {
            return &contract;
        }
    }
    return nullptr;
}

inline MessagePrecodeEncoderOptions MessageOptionsForContract(
    const V2EquationContract& contract)
{
    MessagePrecodeEncoderOptions options;
    options.Architecture = contract.Architecture;
    options.Completion = contract.Completion;
    options.RecoveryMixCount = contract.RecoveryMixCount;
    options.AdaptiveDenseTwoAnchor =
        contract.AdaptiveDenseTwoAnchor;
    if (contract.SeedPolicy == V2SeedDerivation::RawUniform)
    {
        // RawUniform binds both construction seeds directly.  Canonicalize
        // the inherited ProfileDerived salt fields instead of carrying inert
        // aliases under the same raw contract identity.
        options.PrecodeSeedSalt = 0u;
        options.RecoveryRowSeedSalt = 0u;
    }
    return options;
}

/**
    Expand one exact table-free raw target profile.

    This is the single construction path shared by the full codec benchmark and
    the equation fingerprint.  It deliberately starts from a zero SeedProfile,
    binds the construction seeds directly, and never invokes SelectSeedProfile
    or packet/precode attempt stepping.
*/
inline bool MakeRawContractProfile(
    const V2EquationContract& contract,
    uint32_t block_count,
    uint32_t block_bytes,
    uint64_t construction_seed,
    SeedProfile& profile_out)
{
    if (contract.ProfileId != kWh2DispatchV1ContractId ||
        contract.PublicProfile ||
        contract.SeedPolicy != V2SeedDerivation::RawUniform ||
        contract.SeedAttemptCount != 1u ||
        contract.Architecture != V2PrecodeArchitecture::SmallBandD4 ||
        contract.Completion != CompletionField::MixedGF256GF16 ||
        contract.RecoveryMixCount != kCertifiedPacketMixCount ||
        contract.AdaptiveDenseTwoAnchor ||
        block_count < 2u || block_count > 64000u ||
        block_bytes == 0u || block_bytes > UINT32_C(0x7fffffff) ||
        (block_bytes & 1u) != 0u)
    {
        return false;
    }

    // Test-hook builds can carry thread-local equation experiments.
    // The stable dispatch-v1 identifier names the frozen production geometry,
    // so fail closed rather than binding that identifier to different heavy
    // rows, coefficient coordinates, residue schedules, staircase shapes, or
    // packet-id mappings.  The separately receipted peel-PMF and
    // staircase-scale hooks remain valid target overlays.
    if (!IsCanonicalStableTargetEquationState())
    {
        return false;
    }

    SeedProfile profile;
    profile.BlockCount = block_count;
    profile.BlockBytes = block_bytes;
    profile.V2SeedPolicy = V2SeedDerivation::RawUniform;

    PrecodeSystem system;
    system.Params = MakeMixedParams(block_count, construction_seed);
    system.Params.Staircase = SmallBandStaircaseCount(block_count);
    system.Params.DenseRows = 4u;
    system.Params.DenseIdentityCorner = false;
    system.Params.DenseTwoAnchor = false;

    PacketRowConfig packet;
    packet.PeelSeed = RawUniformPacketPeelSeed(construction_seed);
    packet.MixCount = contract.RecoveryMixCount;

    const MessagePrecodeEncoderOptions options =
        MessageOptionsForContract(contract);
    BindMessagePrecodeProfile(profile, options, system, packet, 0u);
    profile_out = profile;
    return true;
}

} // namespace wirehair_v2
