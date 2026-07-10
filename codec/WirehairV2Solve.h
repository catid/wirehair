#pragma once

#include "WirehairV2Precode.h"
#include "WirehairV2Seeds.h"

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <vector>

namespace wirehair_v2 {

// Message packets use a distinct final contract number so they cannot be
// confused with the separately versioned experimental recovery-row helper.
static const uint32_t kPacketRowContractVersion = 4u;
static const uint32_t kPrecodeContractVersion = 2u;
static const uint32_t kCertifiedPacketMixCount = 3u;
static const uint32_t kMaxPacketSeedAttempts = 256u;
static const uint32_t kMaxInactiveColumns = 4096u;
static const uint32_t kMinPacketPrecodeCount = 2u;
static const uint32_t kMaxPacketPrecodeCount = 65521u;

struct PacketRowConfig
{
    uint32_t PeelSeed = 0;
    uint32_t MixCount = kCertifiedPacketMixCount;
};

struct SolvePacket
{
    uint32_t BlockId = 0;
    const uint8_t* Data = nullptr;
};

struct PrecodeSolveStats
{
    uint32_t PacketRows = 0;
    uint32_t PeeledColumns = 0;
    uint32_t InactivatedColumns = 0;
    uint32_t ResidualRows = 0;
    uint32_t ResidualRank = 0;
    uint32_t BinaryResidualRank = 0;
    uint64_t BinaryRowReferences = 0;
    uint64_t BlockXors = 0;
    uint64_t BlockMulAdds = 0;
    uint64_t BuildNanoseconds = 0;
    uint64_t PeelNanoseconds = 0;
    uint64_t ProjectNanoseconds = 0;
    uint64_t ResidualNanoseconds = 0;
    uint64_t BackSubNanoseconds = 0;
    uint32_t PacketSeedAttempt = 0;
};

/**
    Return whether the exact version-4 packet-row iterator supports this
    domain and mix count.

    This is deliberately narrower than BuildPrecodeSystem()'s structural
    domain.  Packet rows require at least two precode columns for parameter
    generation and use NextPrime16(), whose largest supported input is 65521.
    A structurally valid system outside this range may still be inspected or
    used by structure-only tooling, but it cannot be evaluated or solved as a
    version-4 packet profile.
*/
bool IsPacketRowDomainValid(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t mix_count);

/**
    Generate the version-4 packet equation for a public block id.

    The source prefix is contractually bound to production Wirehair's integer
    GeneratePeelRowWeight()/PeelRowIterator rule.  Exactly MixCount distinct
    precode columns are added with RowMixIterator.  Addressing by the public
    block id keeps systematic equations [0,K) disjoint from repair equations.
    Unsupported packet domains return an empty row without invoking either
    iterator.
*/
std::vector<uint32_t> GeneratePacketMatrixRow(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config);

uint32_t PacketPeelSeedFromProfile(
    const SeedProfile& profile,
    uint64_t salt);

PacketRowConfig PacketConfigForAttempt(
    const PacketRowConfig& base,
    uint32_t attempt);

PrecodeParams PrecodeParamsForAttempt(
    const PrecodeParams& base,
    uint32_t attempt);

/** Evaluate one packet row over all intermediate blocks. */
bool EvaluatePacketBlock(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/** Internal fast path for an already validated immutable system/config. */
bool EvaluatePacketBlockForValidatedSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/**
    Solve the complete V2 system over GF(256).

    Binary staircase/dense constraints and packet equations are peeled first.
    Unused binary rows are projected onto the inactivated columns, then the
    actual Cauchy heavy equations are inserted into the same GF(256) residual
    solve.  On success `intermediate_blocks_out` contains all
    K+S+D2+H block values.  The exact residual solve is bounded to
    kMaxInactiveColumns to contain adversarial memory use.  NeedMore means the
    supplied equations were rank deficient or exceeded that bound; additional
    independent packets can reduce the residual.  Output remains unchanged on
    every failure.
*/
WirehairResult SolvePrecodeSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats = nullptr);

/** Select the first deterministic packet seed whose K systematic rows rank. */
WirehairResult SelectSystematicPacketConfig(
    const PrecodeSystem& system,
    const PacketRowConfig& base_config,
    PacketRowConfig& selected_config,
    uint32_t* attempt_out = nullptr);

/** Select the first full-rank deterministic joint precode/packet seed. */
WirehairResult SelectSystematicConfiguration(
    const PrecodeParams& base_params,
    const PacketRowConfig& base_config,
    PrecodeSystem& selected_system,
    PacketRowConfig& selected_config,
    uint32_t* attempt_out = nullptr);

/** Expensive test/oracle validation of every supplied equation. */
bool VerifyPrecodeSolution(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes);

} // namespace wirehair_v2
