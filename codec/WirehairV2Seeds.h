#pragma once

#include "WirehairV2Policy.h"
#include "WirehairV2GF16.h"

#include <stdint.h>

namespace wirehair_v2 {

static const uint16_t kSeedTreeSubdivisions = 2048;

struct SeedProfile
{
    uint32_t BlockCount = 0;
    uint32_t BlockBytes = 0;
    PeelPolicy Policy = {};
    uint16_t DenseCount = 0;
    uint16_t PeelSeed = 0;
    uint16_t DenseSeed = 0;
    uint16_t PeelSeedBucket = 0;
    bool UsedPeelFixup = false;
    bool UsedDenseFixup = false;
    bool V2SeedSelected = false;
    uint32_t V2SeedAttempt = 0;
    uint32_t V2PrecodeContractVersion = 0;
    uint32_t V2PacketRowContractVersion = 0;
    uint32_t V2StaircaseCount = 0;
    uint32_t V2DenseRowCount = 0;
    uint32_t V2HeavyRowCount = 0;
    CompletionField V2CompletionField = CompletionField::GF256;
    uint32_t V2SourceHits = 0;
    uint64_t V2PrecodeSeed = 0;
    uint32_t V2PacketPeelSeed = 0;
    uint32_t V2RecoveryMixCount = 0;
    bool V2DenseIdentityCorner = false;
    bool V2DenseTwoAnchor = false;
    // Profile-level policy bit.  It remains true below the K cutoff even
    // though those small systems retain the certified one-anchor equations.
    bool V2AdaptiveDenseTwoAnchor = false;
    uint64_t V2PrecodeSeedSalt = 0;
    uint64_t V2RecoveryRowSeedSalt = 0;
    bool Tuned = false;
    double TuningResidualMean = 0.0;
    uint32_t TuningResidualColumns = 0;
    uint64_t TuningXorCost = 0;
    uint16_t TuningCandidatesRequested = 0;
    uint16_t TuningCandidatesUnique = 0;
    uint16_t TuningCandidatesCompleted = 0;
    uint32_t TuningTrials = 0;
};

struct SeedTuningOptions
{
    uint16_t PeelCandidates = 0;
    uint32_t TrialsPerCandidate = 0;
    uint64_t Seed = 0;
};

struct SeedTuningCandidate
{
    uint16_t PeelSeed = 0;
    uint16_t DenseSeed = 0;
    double ResidualMean = 0.0;
    uint32_t ResidualMax = 0;
    uint64_t XorCostMean = 0;
    double Score = 0.0;
};

SeedTuningOptions DefaultSeedTuningOptions();

SeedProfile SelectSeedProfile(uint32_t block_count, uint32_t block_bytes);
SeedProfile TuneSeedProfile(
    uint32_t block_count,
    uint32_t block_bytes,
    const SeedTuningOptions& options);

uint16_t PeelSeedBucket(uint32_t block_count);

/// Indices 0..255 enumerate the byte seed domain exactly once, starting at
/// base_seed & 0xff.  Callers must keep index within that domain.
uint16_t CandidatePeelSeed(uint16_t bucket, uint16_t base_seed, uint16_t index);

/// Index zero preserves the full base_seed.  Indices 1..255 are unique byte
/// seeds; when base_seed is itself a byte, 0..255 form a full permutation.
uint16_t CandidateDenseSeed(uint16_t base_seed, uint16_t index);

} // namespace wirehair_v2
