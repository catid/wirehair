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

uint64_t MatrixSeedFromProfile(
    const SeedProfile& profile,
    uint32_t row_count,
    uint64_t salt);

PeelSolvePlan BuildPeelSolvePlan(
    const SeedProfile& profile,
    uint32_t overhead_rows,
    uint64_t salt);

} // namespace wirehair_v2
