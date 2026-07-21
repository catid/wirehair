#pragma once

#include "WirehairV2Peel.h"
#include "WirehairV2Seeds.h"

namespace wirehair_v2 {

struct PeelSolvePlan
{
    SeedProfile Profile;
    uint32_t RowCount;
    uint64_t MatrixSeed;
    std::vector<std::vector<uint16_t> > Rows;
    PeelEvaluation Evaluation;
};

/**
    Canonical BlockBytes constant for the experiment-only payload-independent
    equation-seeding profile (wirehair-sxvz.16.1.5.2).

    WH2 equation algebra never consumes the payload width, but
    MatrixSeedFromProfile() hashes SeedProfile::BlockBytes plus the
    block-byte-classified peel policy, so graph-recovery hotspots currently
    depend on the payload width.  The prior cross-width holdout normalized
    graph seeding to bb=2 and made K-indexed fixups portable across
    bb={2,64,256,1280,4096}, which is why the canonical constant is 2.
*/
static const uint32_t kPayloadIndependentEquationSeedBlockBytes = 2u;

/**
    Pure equation-seed normalization.

    Returns a copy of profile whose only changes are the two
    MatrixSeedFromProfile() inputs that depend on the payload width:
    BlockBytes and the block-byte-classified peel policy, both re-derived at
    canonical_block_bytes.  Every other hashed input (BlockCount, PeelSeed,
    DenseSeed, DenseCount) derives from BlockCount alone and is preserved, so
    the result hashes identically to SelectSeedProfile(BlockCount,
    canonical_block_bytes).  A zero canonical width is a documented no-op so
    CLI zero can mean "normalization off".  Memory layout and kernels must
    keep consuming the original profile.
*/
SeedProfile NormalizeProfileForEquationSeeding(
    const SeedProfile& profile,
    uint32_t canonical_block_bytes);

uint64_t MatrixSeedFromProfile(
    const SeedProfile& profile,
    uint32_t row_count,
    uint64_t salt);

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/**
    Enable the experiment-only payload-independent seeding profile on the
    calling thread.

    While enabled, MatrixSeedFromProfile() (and therefore
    PacketPeelSeedFromProfile() and every V2 precode/packet equation seed
    derived from it) hashes the profile normalized to
    kPayloadIndependentEquationSeedBlockBytes instead of the true payload
    width.  Everything outside seed derivation -- memory layout, block
    kernels, solve arenas -- still uses the true BlockBytes.  This is NOT a
    public V2 profile ID: the hook exists only in test builds, is thread
    local, and defaults to off, so all named/public profiles keep their
    frozen equation inputs byte for byte.
*/
void SetPayloadIndependentEquationSeedingForTesting(bool enabled);

/** Active state of the payload-independent seeding hook on this thread. */
bool PayloadIndependentEquationSeedingForTesting();
#endif

PeelSolvePlan BuildPeelSolvePlan(
    const SeedProfile& profile,
    uint32_t overhead_rows,
    uint64_t salt);

} // namespace wirehair_v2
