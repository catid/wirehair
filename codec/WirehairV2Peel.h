#pragma once

#include "WirehairV2Policy.h"

#include <stdint.h>
#include <vector>

namespace wirehair_v2 {

static const uint32_t kMaxPeelMatrixRows = 65536u;

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

/**
    Generate one peel row by replaying the same deterministic row stream used
    by GeneratePeelMatrixRows().  `row_index` is zero-based and must be below
    kMaxPeelMatrixRows.
*/
std::vector<uint16_t> GeneratePeelMatrixRow(
    const PeelingCodec& codec,
    uint32_t block_count,
    uint32_t row_index,
    uint64_t seed);

/**
    Generate recovery rows over the full V2 intermediate-symbol domain:
    source columns [0, source_count) followed by precode columns
    [source_count, source_count + precode_count).

    The source-column prefix is identical to GeneratePeelMatrixRows() for the
    same codec/source_count/row_count/seed.  Precode mix columns are drawn from
    a separate deterministic stream so changing precode_count does not retune
    the source-row distribution.
*/
std::vector<std::vector<uint32_t> > GenerateRecoveryMatrixRows(
    const PeelingCodec& codec,
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t row_count,
    uint32_t mix_count,
    uint64_t seed);

/**
    Generate one recovery row by replaying the same deterministic row stream
    used by GenerateRecoveryMatrixRows().  This is the encoder-facing mapping
    from recovery block index to the V2 intermediate-symbol row; callers that
    need many adjacent rows should use the batch API above.
*/
std::vector<uint32_t> GenerateRecoveryMatrixRow(
    const PeelingCodec& codec,
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t row_index,
    uint32_t mix_count,
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
