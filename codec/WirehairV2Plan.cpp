#include "WirehairV2Plan.h"

namespace wirehair_v2 {
namespace {

uint64_t Mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

} // namespace

uint64_t MatrixSeedFromProfile(
    const SeedProfile& profile,
    uint32_t row_count,
    uint64_t salt)
{
    uint64_t seed = salt;
    seed ^= Mix64((uint64_t)profile.BlockCount * UINT64_C(0x9e3779b97f4a7c15));
    seed ^= Mix64((uint64_t)profile.BlockBytes * UINT64_C(0xbf58476d1ce4e5b9));
    seed ^= Mix64((uint64_t)row_count * UINT64_C(0x94d049bb133111eb));
    seed ^= Mix64((uint64_t)profile.PeelSeed * UINT64_C(0xd6e8feb86659fd93));
    seed ^= Mix64((uint64_t)profile.DenseSeed * UINT64_C(0xa0761d6478bd642f));
    seed ^= Mix64((uint64_t)profile.DenseCount * UINT64_C(0xe7037ed1a0b428db));
    seed ^= Mix64((uint64_t)profile.Policy.Codec.MinDegree << 32 |
        profile.Policy.Codec.MaxDegree);
    seed ^= Mix64((uint64_t)profile.Policy.Codec.SolverCandidateLimit << 48 |
        static_cast<uint64_t>(profile.Policy.Structure) << 16 |
        static_cast<uint64_t>(profile.Policy.Solver));
    return Mix64(seed);
}

PeelSolvePlan BuildPeelSolvePlan(
    const SeedProfile& profile,
    uint32_t overhead_rows,
    uint64_t salt)
{
    PeelSolvePlan plan;
    plan.Profile = profile;
    plan.RowCount = profile.BlockCount + overhead_rows;
    plan.MatrixSeed = MatrixSeedFromProfile(profile, plan.RowCount, salt);
    plan.Rows = GeneratePeelMatrixRows(
        profile.Policy.Codec,
        profile.BlockCount,
        plan.RowCount,
        plan.MatrixSeed);
    plan.Evaluation = EvaluatePeelingRows(
        profile.Policy.Codec,
        profile.BlockCount,
        plan.Rows);
    return plan;
}

} // namespace wirehair_v2
