#include "WirehairV2Policy.h"

namespace wirehair_v2 {
namespace {

const uint32_t kSmallMediumBlockBytes = 12u * 1024u;
const uint32_t kMediumLargeBlockBytes = 320u * 1024u;

PeelingCodec MakeLtCodec(
    PeelStructure structure,
    PeelSolver solver,
    uint16_t min_degree,
    uint16_t max_degree)
{
    PeelingCodec codec;
    codec.Solver = solver;
    codec.Structure = structure;
    codec.Family = DegreeFamily::Lt;
    codec.MinDegree = min_degree;
    codec.MaxDegree = max_degree;
    codec.SolverCandidateLimit =
        solver == PeelSolver::KsBmaxTop16 ? 16u : 0u;
    codec.Degree1Mass = 0.0;
    codec.Degree2Mass = 0.0;
    codec.RobustC = 0.0;
    codec.RobustDelta = 0.0;
    codec.FullyRandomRows = true;
    return codec;
}

PeelingCodec MakeRobustD1D2Codec(
    PeelStructure structure,
    PeelSolver solver,
    uint16_t max_degree,
    double degree1_mass,
    double degree2_mass)
{
    PeelingCodec codec;
    codec.Solver = solver;
    codec.Structure = structure;
    codec.Family = DegreeFamily::RobustD1D2;
    codec.MinDegree = 1u;
    codec.MaxDegree = max_degree;
    codec.SolverCandidateLimit =
        solver == PeelSolver::KsBmaxTop16 ? 16u : 0u;
    codec.Degree1Mass = degree1_mass;
    codec.Degree2Mass = degree2_mass;
    codec.RobustC = 0.0;
    codec.RobustDelta = 0.0;
    codec.FullyRandomRows = true;
    return codec;
}

PeelingCodec MakeRobustSolitonCodec(
    PeelStructure structure,
    PeelSolver solver,
    double c,
    double delta,
    uint16_t max_degree)
{
    PeelingCodec codec;
    codec.Solver = solver;
    codec.Structure = structure;
    codec.Family = DegreeFamily::RobustSoliton;
    codec.MinDegree = 1u;
    codec.MaxDegree = max_degree;
    codec.SolverCandidateLimit =
        solver == PeelSolver::KsBmaxTop16 ? 16u : 0u;
    codec.Degree1Mass = 0.0;
    codec.Degree2Mass = 0.0;
    codec.RobustC = c;
    codec.RobustDelta = delta;
    codec.FullyRandomRows = true;
    return codec;
}

PeelStructure SelectStructure(BlockByteClass byte_class, BlockCountBand count_band)
{
    switch (byte_class) {
    case BlockByteClass::Small:
        switch (count_band) {
        case BlockCountBand::UpTo1000:
            return PeelStructure::RobustD1_001D2_012;
        case BlockCountBand::UpTo3200:
            return PeelStructure::LtM1C32;
        case BlockCountBand::UpTo12000:
            return PeelStructure::RobustD1_001D2_003;
        case BlockCountBand::Above12000:
            return PeelStructure::LtM2C96;
        }
        break;
    case BlockByteClass::Medium:
        switch (count_band) {
        case BlockCountBand::UpTo1000:
            return PeelStructure::LtM1C16;
        case BlockCountBand::UpTo3200:
            return PeelStructure::LtM1C32;
        case BlockCountBand::UpTo12000:
            return PeelStructure::LtM1C64;
        case BlockCountBand::Above12000:
            return PeelStructure::RsC001D50C128;
        }
        break;
    case BlockByteClass::Large:
        switch (count_band) {
        case BlockCountBand::UpTo1000:
            return PeelStructure::LtM1C16;
        case BlockCountBand::UpTo3200:
            return PeelStructure::LtM1C32;
        case BlockCountBand::UpTo12000:
            return PeelStructure::LtM1C64;
        case BlockCountBand::Above12000:
            return PeelStructure::RsC003D10C128;
        }
        break;
    }

    return PeelStructure::LtM1C32;
}

} // namespace

BlockByteClass ClassifyBlockBytes(uint32_t block_bytes)
{
    if (block_bytes < kSmallMediumBlockBytes) {
        return BlockByteClass::Small;
    }
    if (block_bytes < kMediumLargeBlockBytes) {
        return BlockByteClass::Medium;
    }
    return BlockByteClass::Large;
}

BlockCountBand ClassifyBlockCount(uint32_t block_count)
{
    if (block_count <= 1000u) {
        return BlockCountBand::UpTo1000;
    }
    if (block_count <= 3200u) {
        return BlockCountBand::UpTo3200;
    }
    if (block_count <= 12000u) {
        return BlockCountBand::UpTo12000;
    }
    return BlockCountBand::Above12000;
}

PeelingCodec MakePeelingCodec(PeelStructure structure, PeelSolver solver)
{
    switch (structure) {
    case PeelStructure::LtM1C16:
        return MakeLtCodec(structure, solver, 1u, 16u);
    case PeelStructure::LtM1C32:
        return MakeLtCodec(structure, solver, 1u, 32u);
    case PeelStructure::LtM1C64:
        return MakeLtCodec(structure, solver, 1u, 64u);
    case PeelStructure::LtM2C96:
        return MakeLtCodec(structure, solver, 2u, 96u);
    case PeelStructure::RobustD1_001D2_003:
        return MakeRobustD1D2Codec(structure, solver, 64u, 0.01, 0.03);
    case PeelStructure::RobustD1_001D2_012:
        return MakeRobustD1D2Codec(structure, solver, 64u, 0.01, 0.12);
    case PeelStructure::RsC001D50C128:
        return MakeRobustSolitonCodec(structure, solver, 0.01, 0.50, 128u);
    case PeelStructure::RsC003D10C128:
        return MakeRobustSolitonCodec(structure, solver, 0.03, 0.10, 128u);
    }

    return MakeLtCodec(PeelStructure::LtM1C32, solver, 1u, 32u);
}

PeelingCodec SelectPeelingCodec(uint32_t block_count, uint32_t block_bytes)
{
    const BlockByteClass byte_class = ClassifyBlockBytes(block_bytes);
    const BlockCountBand count_band = ClassifyBlockCount(block_count);
    const PeelStructure structure = SelectStructure(byte_class, count_band);
    return MakePeelingCodec(structure, PeelSolver::KsBmaxTop16);
}

PeelPolicy SelectPeelPolicy(uint32_t block_count, uint32_t block_bytes)
{
    const BlockByteClass byte_class = ClassifyBlockBytes(block_bytes);
    const BlockCountBand count_band = ClassifyBlockCount(block_count);
    const PeelingCodec codec = SelectPeelingCodec(block_count, block_bytes);

    PeelPolicy policy;
    policy.Solver = codec.Solver;
    policy.Structure = codec.Structure;
    policy.ByteClass = byte_class;
    policy.CountBand = count_band;
    policy.Codec = codec;
    return policy;
}

const char* ToString(PeelSolver solver)
{
    switch (solver) {
    case PeelSolver::RqccLowref:
        return "rqcc_lowref";
    case PeelSolver::KsBmaxTop16:
        return "ks_bmax_top16";
    }
    return "unknown";
}

const char* ToString(PeelStructure structure)
{
    switch (structure) {
    case PeelStructure::LtM1C16:
        return "lt_m1_c16";
    case PeelStructure::LtM1C32:
        return "lt_m1_c32";
    case PeelStructure::LtM1C64:
        return "lt_m1_c64";
    case PeelStructure::LtM2C96:
        return "lt_m2_c96";
    case PeelStructure::RobustD1_001D2_003:
        return "robust_d1_001_d2_003";
    case PeelStructure::RobustD1_001D2_012:
        return "robust_d1_001_d2_012";
    case PeelStructure::RsC001D50C128:
        return "rs_c001_d50_c128";
    case PeelStructure::RsC003D10C128:
        return "rs_c003_d10_c128";
    }
    return "unknown";
}

const char* ToString(BlockByteClass byte_class)
{
    switch (byte_class) {
    case BlockByteClass::Small:
        return "small";
    case BlockByteClass::Medium:
        return "medium";
    case BlockByteClass::Large:
        return "large";
    }
    return "unknown";
}

const char* ToString(BlockCountBand count_band)
{
    switch (count_band) {
    case BlockCountBand::UpTo1000:
        return "n_le_1000";
    case BlockCountBand::UpTo3200:
        return "n_le_3200";
    case BlockCountBand::UpTo12000:
        return "n_le_12000";
    case BlockCountBand::Above12000:
        return "n_gt_12000";
    }
    return "unknown";
}

const char* ToString(DegreeFamily family)
{
    switch (family) {
    case DegreeFamily::Lt:
        return "lt";
    case DegreeFamily::RobustD1D2:
        return "robust_d1_d2";
    case DegreeFamily::RobustSoliton:
        return "robust_soliton";
    }
    return "unknown";
}

} // namespace wirehair_v2
