#include "WirehairV2Seeds.h"

#include "../WirehairTools.h"

#include "WirehairV2Peel.h"

#include <math.h>

namespace wirehair_v2 {
namespace {

uint32_t Hash32(uint32_t x)
{
    x ^= x >> 16;
    x *= UINT32_C(0x7feb352d);
    x ^= x >> 15;
    x *= UINT32_C(0x846ca68b);
    x ^= x >> 16;
    return x;
}

bool HasPeelFixup(uint32_t block_count, uint16_t selected_seed)
{
    uint16_t table_seed = 0;
    if (block_count < wirehair::kTinyTableCount + wirehair::kSmallTableCount) {
        table_seed = wirehair::kSmallPeelSeeds[block_count];
    }
    else {
        table_seed =
            wirehair::kPeelSeeds[block_count % wirehair::kPeelSeedSubdivisions];
    }
    return selected_seed != table_seed;
}

bool HasDenseFixup(
    uint32_t block_count,
    uint16_t dense_count,
    uint16_t selected_seed)
{
    uint16_t table_seed = 0;
    if (block_count < wirehair::kTinyTableCount) {
        table_seed = wirehair::kTinyDenseSeeds[block_count];
    }
    else if (block_count < wirehair::kTinyTableCount + wirehair::kSmallTableCount) {
        table_seed =
            wirehair::kSmallDenseSeeds[block_count - wirehair::kTinyTableCount];
    }
    else {
        const unsigned table_index = dense_count / 4u;
        if (table_index < wirehair::kDenseSeedCount) {
            table_seed = wirehair::kDenseSeeds[table_index];
        }
    }
    return selected_seed != table_seed;
}

SeedProfile MakeBaseProfile(uint32_t block_count, uint32_t block_bytes)
{
    SeedProfile profile;
    profile.BlockCount = block_count;
    profile.BlockBytes = block_bytes;
    profile.Policy = SelectPeelPolicy(block_count, block_bytes);
    profile.DenseCount = wirehair::GetDenseCount(block_count);
    profile.PeelSeed = wirehair::GetPeelSeed(block_count);
    profile.DenseSeed = wirehair::GetDenseSeed(block_count, profile.DenseCount);
    profile.PeelSeedBucket = PeelSeedBucket(block_count);
    profile.UsedPeelFixup = HasPeelFixup(block_count, profile.PeelSeed);
    profile.UsedDenseFixup =
        HasDenseFixup(block_count, profile.DenseCount, profile.DenseSeed);
    profile.Tuned = false;
    profile.TuningResidualMean = 0.0;
    profile.TuningResidualColumns = 0;
    profile.TuningXorCost = 0;
    return profile;
}

double CandidateScore(
    double residual_mean,
    double residual_population_mean,
    double residual_population_sd,
    uint64_t xor_cost)
{
    const double residual_limit =
        residual_population_mean + 1.5 * residual_population_sd + 2.0;
    const double outlier =
        residual_mean > residual_limit ? residual_mean - residual_limit : 0.0;
    return outlier * outlier * 1000000.0 + (double)xor_cost;
}

} // namespace

SeedTuningOptions DefaultSeedTuningOptions()
{
    SeedTuningOptions options;
    options.PeelCandidates = 16;
    options.TrialsPerCandidate = 3;
    options.Seed = UINT64_C(0x9a7e11a);
    return options;
}

uint16_t PeelSeedBucket(uint32_t block_count)
{
    return (uint16_t)(block_count % kSeedTreeSubdivisions);
}

uint16_t CandidatePeelSeed(uint16_t bucket, uint16_t base_seed, uint16_t index)
{
    if (index == 0) {
        return (uint16_t)(base_seed & 0xffu);
    }
    const uint32_t h = Hash32(
        ((uint32_t)bucket << 16) ^
        ((uint32_t)base_seed << 8) ^
        (uint32_t)index ^
        UINT32_C(0x7065656c));
    return (uint16_t)(h & 0xffu);
}

uint16_t CandidateDenseSeed(uint16_t base_seed, uint16_t index)
{
    if (index == 0) {
        return (uint16_t)(base_seed & 0xffu);
    }
    const uint32_t h = Hash32(
        ((uint32_t)base_seed << 9) ^
        (uint32_t)index ^
        UINT32_C(0x64656e73));
    return (uint16_t)(h & 0xffu);
}

SeedProfile SelectSeedProfile(uint32_t block_count, uint32_t block_bytes)
{
    return MakeBaseProfile(block_count, block_bytes);
}

SeedProfile TuneSeedProfile(
    uint32_t block_count,
    uint32_t block_bytes,
    const SeedTuningOptions& options)
{
    SeedProfile base = MakeBaseProfile(block_count, block_bytes);
    if (block_count < 2u) {
        return base;
    }

    uint16_t peel_candidates =
        options.PeelCandidates > 0u ? options.PeelCandidates : 1u;
    if (peel_candidates > 256u) {
        peel_candidates = 256u;
    }
    const uint16_t trials =
        options.TrialsPerCandidate > 0u ? options.TrialsPerCandidate : 1u;

    SeedTuningCandidate candidates[256];
    double residual_sum = 0.0;
    double residual_sq = 0.0;
    for (uint16_t i = 0; i < peel_candidates; ++i)
    {
        SeedTuningCandidate candidate;
        candidate.PeelSeed =
            CandidatePeelSeed(base.PeelSeedBucket, base.PeelSeed, i);
        candidate.DenseSeed = base.DenseSeed;
        candidate.ResidualMean = 0.0;
        candidate.ResidualMax = 0;
        candidate.XorCostMean = 0;
        candidate.Score = 0.0;

        uint64_t xor_sum = 0;
        for (uint16_t trial = 0; trial < trials; ++trial)
        {
            const uint64_t matrix_seed = options.Seed ^
                ((uint64_t)block_count * UINT64_C(0xbf58476d1ce4e5b9)) ^
                ((uint64_t)candidate.PeelSeed * UINT64_C(0x94d049bb133111eb)) ^
                ((uint64_t)trial * UINT64_C(0xd6e8feb86659fd93));
            const PeelEvaluation eval =
                EvaluatePeeling(base.Policy.Codec, block_count, matrix_seed);
            candidate.ResidualMean += eval.ResidualColumns;
            if (eval.ResidualColumns > candidate.ResidualMax) {
                candidate.ResidualMax = eval.ResidualColumns;
            }
            xor_sum += eval.TotalXorCost;
        }
        candidate.ResidualMean /= (double)trials;
        candidate.XorCostMean = xor_sum / trials;
        residual_sum += candidate.ResidualMean;
        residual_sq += candidate.ResidualMean * candidate.ResidualMean;
        candidates[i] = candidate;
    }

    const double population_mean = residual_sum / (double)peel_candidates;
    double population_var =
        residual_sq / (double)peel_candidates - population_mean * population_mean;
    if (population_var < 0.0) {
        population_var = 0.0;
    }
    const double population_sd = sqrt(population_var);

    uint16_t best_i = 0;
    for (uint16_t i = 0; i < peel_candidates; ++i)
    {
        candidates[i].Score = CandidateScore(
            candidates[i].ResidualMean,
            population_mean,
            population_sd,
            candidates[i].XorCostMean);
        if (candidates[i].Score < candidates[best_i].Score) {
            best_i = i;
        }
    }

    base.PeelSeed = candidates[best_i].PeelSeed;
    base.DenseSeed = candidates[best_i].DenseSeed;
    base.Tuned = true;
    base.TuningResidualMean = candidates[best_i].ResidualMean;
    base.TuningResidualColumns = candidates[best_i].ResidualMax;
    base.TuningXorCost = candidates[best_i].XorCostMean;
    return base;
}

} // namespace wirehair_v2
