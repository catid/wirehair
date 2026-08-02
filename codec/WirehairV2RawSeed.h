#pragma once

#include "WirehairV2Precode.h"
#include "WirehairV2Solve.h"

#include <stdint.h>

namespace wirehair_v2 {
namespace test {

// Closed seed basis for raw architecture comparisons.  These fixed base
// values are used directly instead of deriving precode or packet-row seeds
// from the legacy K-indexed peel/dense tables.  Construction attempts are
// applied only after raw base construction, so the attempt seed schedule is
// identical for every K and every experimental arm.
//
// The bases were fixed without a codec search: they are the little-endian
// first 64 and 32 bits of SHA-256("wirehair.wh2.raw-architecture.precode.v1")
// and SHA-256("wirehair.wh2.raw-architecture.packet.v1"), respectively.
// This is benchmark/test plumbing, not a production profile or wire-format
// choice.  Keep the canonical name and values bound into every experimental
// arm descriptor and result provenance that uses it.
static const uint64_t kRawArchitecturePrecodeSeed =
    UINT64_C(0x487468302aad7105);
static const uint32_t kRawArchitecturePacketSeed = UINT32_C(0x4ec72102);
static const char kRawArchitectureSeedBasis[] =
    "uniform-raw-v1";
static const char kRawArchitectureSeedScheduleCanonical[] =
    "{\"excluded_inputs\":[\"K\",\"block_bytes\","
    "\"production_peel_seed\",\"production_dense_seed\","
    "\"production_seed_fixups\"],\"inputs\":[\"precode_attempt\","
    "\"packet_attempt\"],\"packet_attempt_stride\":\"0x9e3779b9\","
    "\"packet_base_seed\":\"0x4ec72102\",\"precode_attempt_stride\":"
    "\"0x9e3779b97f4a7c15\",\"precode_base_seed\":"
    "\"0x487468302aad7105\",\"schema\":"
    "\"wirehair.wh2.uniform-raw-seed-schedule.v1\"}";
static const char kRawArchitectureSeedScheduleSha256[] =
    "90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7";

// Build the complete raw comparison base without consulting a production
// SeedProfile, a K-indexed seed table, or any production seed fixup.  Only
// MakeCertifiedParams supplies the structural K-dependent dimensions.
inline void MakeRawArchitectureConfiguration(
    uint32_t block_count,
    PrecodeParams& params,
    PacketRowConfig& packet)
{
    params = MakeCertifiedParams(
        block_count, kRawArchitecturePrecodeSeed);
    packet = PacketRowConfig();
    packet.PeelSeed = kRawArchitecturePacketSeed;
    packet.MixCount = kCertifiedPacketMixCount;
}

} // namespace test
} // namespace wirehair_v2
