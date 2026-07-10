#pragma once

#include "WirehairV2Policy.h"

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
    bool Tuned = false;
    double TuningResidualMean = 0.0;
    uint32_t TuningResidualColumns = 0;
    uint64_t TuningXorCost = 0;
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
uint16_t CandidatePeelSeed(uint16_t bucket, uint16_t base_seed, uint16_t index);
uint16_t CandidateDenseSeed(uint16_t base_seed, uint16_t index);

} // namespace wirehair_v2
