#pragma once

#include "WirehairV2Policy.h"

#include <stdint.h>

namespace wirehair_v2 {

static const uint16_t kSeedTreeSubdivisions = 2048;

struct SeedProfile
{
    uint32_t BlockCount;
    uint32_t BlockBytes;
    PeelPolicy Policy;
    uint16_t DenseCount;
    uint16_t PeelSeed;
    uint16_t DenseSeed;
    uint16_t PeelSeedBucket;
    bool UsedPeelFixup;
    bool UsedDenseFixup;
    bool Tuned;
    double TuningResidualMean;
    uint32_t TuningResidualColumns;
    uint64_t TuningXorCost;
};

struct SeedTuningOptions
{
    uint16_t PeelCandidates;
    uint16_t TrialsPerCandidate;
    uint64_t Seed;
};

struct SeedTuningCandidate
{
    uint16_t PeelSeed;
    uint16_t DenseSeed;
    double ResidualMean;
    uint32_t ResidualMax;
    uint64_t XorCostMean;
    double Score;
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
