#pragma once

#include "WirehairV2Policy.h"

#include <stdint.h>
#include <vector>

namespace wirehair_v2 {

struct PeelEvaluation
{
    uint32_t Rows;
    uint32_t Columns;
    uint32_t ResidualRows;
    uint32_t ResidualColumns;
    uint64_t MatrixRefs;
    uint64_t MatrixXors;
    uint64_t SolveDenseXors;
    uint64_t TotalXorCost;
};

std::vector<std::vector<uint16_t> > GeneratePeelMatrixRows(
    const PeelingCodec& codec,
    uint32_t block_count,
    uint32_t row_count,
    uint64_t seed);

PeelEvaluation EvaluatePeeling(
    const PeelingCodec& codec,
    uint32_t block_count,
    uint64_t seed);

PeelEvaluation EvaluatePeelingRows(
    const PeelingCodec& codec,
    uint32_t block_count,
    const std::vector<std::vector<uint16_t> >& rows);

} // namespace wirehair_v2
