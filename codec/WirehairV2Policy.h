#pragma once

#include <stdint.h>

namespace wirehair_v2 {

enum class PeelSolver
{
    RqccLowref,
    KsBmaxTop16
};

enum class PeelStructure
{
    LtM1C16,
    LtM1C32,
    LtM1C64,
    LtM2C96,
    RobustD1_001D2_003,
    RobustD1_001D2_012,
    RsC001D50C128,
    RsC003D10C128
};

enum class BlockByteClass
{
    Small,
    Medium,
    Large
};

enum class BlockCountBand
{
    UpTo1000,
    UpTo3200,
    UpTo12000,
    Above12000
};

enum class DegreeFamily
{
    Lt,
    RobustD1D2,
    RobustSoliton
};

struct PeelingCodec
{
    PeelSolver Solver;
    PeelStructure Structure;
    DegreeFamily Family;
    uint16_t MinDegree;
    uint16_t MaxDegree;
    uint16_t SolverCandidateLimit;
    double Degree1Mass;
    double Degree2Mass;
    double RobustC;
    double RobustDelta;
    bool FullyRandomRows;
};

struct PeelPolicy
{
    PeelSolver Solver;
    PeelStructure Structure;
    BlockByteClass ByteClass;
    BlockCountBand CountBand;
    PeelingCodec Codec;
};

BlockByteClass ClassifyBlockBytes(uint32_t block_bytes);
BlockCountBand ClassifyBlockCount(uint32_t block_count);
PeelingCodec MakePeelingCodec(PeelStructure structure, PeelSolver solver);
PeelingCodec SelectPeelingCodec(uint32_t block_count, uint32_t block_bytes);
PeelPolicy SelectPeelPolicy(uint32_t block_count, uint32_t block_bytes);

const char* ToString(PeelSolver solver);
const char* ToString(PeelStructure structure);
const char* ToString(BlockByteClass byte_class);
const char* ToString(BlockCountBand count_band);
const char* ToString(DegreeFamily family);

} // namespace wirehair_v2
