#pragma once

#include "WirehairV2Policy.h"

#include <stdint.h>
#include <vector>

namespace wirehair_v2 {

static const uint32_t kMaxPeelMatrixRows = 65536u;
static const uint32_t kRecoveryRowContractVersion = 2u;

struct RecoveryRowGenerationStats
{
    uint32_t ContractVersion = kRecoveryRowContractVersion;
    uint64_t SeekWork = 0;
    uint64_t SourceRandomDraws = 0;
    uint64_t MixRandomDraws = 0;
};

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
    Generate one row from the version-2 row-index-addressable stream used by
    GeneratePeelMatrixRows().  Every uint32 row index is valid.
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
    Generate one recovery row from the version-2 row-index-addressable stream.
    Every uint32 row index is valid and lookup work is independent of earlier
    row IDs.  `stats` reports deterministic stream-seek and draw counts.
*/
std::vector<uint32_t> GenerateRecoveryMatrixRow(
    const PeelingCodec& codec,
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t row_index,
    uint32_t mix_count,
    uint64_t seed,
    RecoveryRowGenerationStats* stats = nullptr);

PeelEvaluation EvaluatePeeling(
    const PeelingCodec& codec,
    uint32_t block_count,
    uint64_t seed);

PeelEvaluation EvaluatePeelingRows(
    const PeelingCodec& codec,
    uint32_t block_count,
    const std::vector<std::vector<uint16_t> >& rows);

} // namespace wirehair_v2
