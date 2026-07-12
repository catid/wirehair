#pragma once

#include "WirehairV2Solve.h"

#include <cstddef>
#include <vector>

namespace wirehair_v2 {
namespace test {

static const uint32_t kTinyOracleMaxSourceBlocks = 128u;
static const uint32_t kTinyOracleMaxBlockBytes = 4096u;
static const size_t kTinyOracleMaxAllocationBytes =
    (size_t)64u * 1024u * 1024u;

/**
    Independent bounded dense GF(256) oracle for the GF256-only V2 precode
    equations. MixedGF256GF16 systems are rejected as InvalidInput; their
    coefficients cannot be represented by this oracle's byte matrix.

    It accepts arbitrary packet ids, duplicate/conflicting equations, and
    overdetermined inputs.  Success, rank deficiency, and inconsistency map to
    the same WirehairResult values as SolvePrecodeSystem().  On every failure
    intermediate_blocks_out is unchanged.  The source-domain cap is 128 and
    no individual allocation is allowed to exceed 64 MiB, making this safe for
    deterministic and coverage-guided fuzz targets.
*/
WirehairResult SolvePrecodeSystemTinyDenseOracle(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out);

} // namespace test
} // namespace wirehair_v2
