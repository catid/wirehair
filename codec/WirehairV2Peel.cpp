#include "WirehairV2Peel.h"

#include <algorithm>
#include <limits.h>
#include <math.h>
#include <queue>

namespace wirehair_v2 {
namespace {

struct Rng
{
    uint64_t State;

    explicit Rng(uint64_t seed) : State(seed) {}

    uint64_t Next()
    {
        uint64_t z = (State += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }

    uint32_t U32()
    {
        return (uint32_t)(Next() >> 32);
    }

    double Unit()
    {
        return (Next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

uint32_t ClampDegree(uint32_t degree, uint32_t block_count)
{
    if (degree < 1u) {
        degree = 1u;
    }
    if (degree > block_count) {
        degree = block_count;
    }
    return degree;
}

double LtWeight(uint32_t degree, uint32_t block_count)
{
    if (degree <= 1u) {
        return 1.0 / (double)block_count;
    }
    return 1.0 / ((double)degree * (double)(degree - 1u));
}

double RobustSolitonWeight(
    uint32_t degree,
    uint32_t block_count,
    double c,
    double delta)
{
    if (c <= 0.0) {
        c = 0.01;
    }
    if (delta <= 0.0 || delta >= 1.0) {
        delta = 0.50;
    }

    const double n = (double)block_count;
    double R = c * log(n / delta) * sqrt(n);
    if (R < 1.0) {
        R = 1.0;
    }
    uint32_t spike = (uint32_t)floor(n / R);
    if (spike < 1u) {
        spike = 1u;
    }
    if (spike > block_count) {
        spike = block_count;
    }

    double tau = 0.0;
    if (degree < spike) {
        tau = R / ((double)degree * n);
    }
    else if (degree == spike) {
        tau = R * log(R / delta) / n;
    }
    return LtWeight(degree, block_count) + tau;
}

double DegreeWeight(
    const PeelingCodec& codec,
    uint32_t degree,
    uint32_t block_count)
{
    switch (codec.Family)
    {
    case DegreeFamily::Lt:
        return LtWeight(degree, block_count);
    case DegreeFamily::RobustD1D2:
    {
        double w = LtWeight(degree, block_count);
        if (degree == 1u) {
            w += codec.Degree1Mass;
        }
        else if (degree == 2u) {
            w += codec.Degree2Mass;
        }
        return w;
    }
    case DegreeFamily::RobustSoliton:
        return RobustSolitonWeight(
            degree, block_count, codec.RobustC, codec.RobustDelta);
    }
    return LtWeight(degree, block_count);
}

uint32_t ChooseDegree(
    const PeelingCodec& codec,
    uint32_t block_count,
    Rng& rng)
{
    const uint32_t min_degree = ClampDegree(codec.MinDegree, block_count);
    const uint32_t max_degree = ClampDegree(codec.MaxDegree, block_count);
    double total = 0.0;
    for (uint32_t d = min_degree; d <= max_degree; ++d) {
        total += DegreeWeight(codec, d, block_count);
    }
    if (total <= 0.0) {
        return min_degree;
    }
    double target = rng.Unit() * total;
    for (uint32_t d = min_degree; d <= max_degree; ++d)
    {
        const double w = DegreeWeight(codec, d, block_count);
        if (target <= w) {
            return d;
        }
        target -= w;
    }
    return max_degree;
}

std::vector<uint16_t> RandomRowColumns(
    uint32_t block_count,
    uint32_t degree,
    Rng& rng)
{
    degree = ClampDegree(degree, block_count);
    std::vector<uint16_t> row;
    row.reserve(degree);

    if (degree * 4u >= block_count)
    {
        std::vector<uint16_t> deck;
        deck.reserve(block_count);
        for (uint32_t i = 0; i < block_count; ++i) {
            deck.push_back((uint16_t)i);
        }
        for (uint32_t i = 0; i < degree; ++i)
        {
            const uint32_t j = i + (rng.U32() % (block_count - i));
            const uint16_t value = deck[j];
            deck[j] = deck[i];
            deck[i] = value;
            row.push_back(value);
        }
        return row;
    }

    while (row.size() < degree)
    {
        const uint16_t column = (uint16_t)(rng.U32() % block_count);
        bool duplicate = false;
        for (uint16_t existing : row)
        {
            if (existing == column) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            row.push_back(column);
        }
    }
    return row;
}

struct RowState
{
    std::vector<uint16_t> Columns;
    uint16_t Live;
};

struct ColumnState
{
    std::vector<uint16_t> Rows;
    bool Todo;
};

class PeelSolverState
{
public:
    PeelSolverState(
        const PeelingCodec& codec,
        uint32_t block_count,
        const std::vector<std::vector<uint16_t> >& rows)
        : Codec(codec),
          BlockCount(block_count),
          Rows(rows.size()),
          Columns(block_count),
          DeferredRows(rows.size(), 0)
    {
        for (uint32_t c = 0; c < block_count; ++c) {
            Columns[c].Todo = true;
        }
        for (uint32_t r = 0; r < rows.size(); ++r)
        {
            Rows[r].Columns = rows[r];
            Rows[r].Live = (uint16_t)rows[r].size();
            for (uint16_t c : rows[r]) {
                Columns[c].Rows.push_back((uint16_t)r);
            }
            if (Rows[r].Live == 1u) {
                Queue.push((uint16_t)r);
            }
        }
        TodoColumns = block_count;
    }

    PeelEvaluation Run()
    {
        while (TodoColumns > 0u)
        {
            ProcessSingletons();
            if (TodoColumns == 0u) {
                break;
            }
            const uint16_t column = SelectInactivationColumn();
            if (column >= BlockCount || !Columns[column].Todo) {
                break;
            }
            ++ResidualColumns;
            MarkResidualRows(column);
            RemoveColumn(column);
        }

        PeelEvaluation eval;
        eval.Rows = (uint32_t)Rows.size();
        eval.Columns = BlockCount;
        eval.ResidualRows = ResidualRows;
        eval.ResidualColumns = ResidualColumns;
        eval.MatrixRefs = MatrixRefs;
        eval.MatrixXors = MatrixXors;
        eval.SolveDenseXors =
            (uint64_t)ResidualColumns *
            (uint64_t)std::max(ResidualColumns, ResidualRows);
        eval.TotalXorCost =
            eval.MatrixXors + eval.SolveDenseXors +
            (uint64_t)ResidualColumns * (uint64_t)BlockCount;
        return eval;
    }

    void SetMatrixWork(uint64_t refs, uint64_t xors)
    {
        MatrixRefs = refs;
        MatrixXors = xors;
    }

private:
    void ProcessSingletons()
    {
        while (!Queue.empty())
        {
            const uint16_t row_i = Queue.front();
            Queue.pop();
            if (row_i >= Rows.size() || Rows[row_i].Live != 1u) {
                continue;
            }
            const uint16_t column = OnlyLiveColumn(row_i);
            if (column < BlockCount && Columns[column].Todo) {
                RemoveColumn(column);
            }
        }
    }

    uint16_t OnlyLiveColumn(uint16_t row_i) const
    {
        for (uint16_t c : Rows[row_i].Columns)
        {
            if (Columns[c].Todo) {
                return c;
            }
        }
        return UINT16_MAX;
    }

    void RemoveColumn(uint16_t column)
    {
        Columns[column].Todo = false;
        --TodoColumns;
        for (uint16_t row_i : Columns[column].Rows)
        {
            RowState& row = Rows[row_i];
            if (row.Live == 0u) {
                continue;
            }
            --row.Live;
            if (row.Live == 1u) {
                Queue.push(row_i);
            }
        }
    }

    uint32_t CountLiveRowsForColumn(uint16_t column) const
    {
        uint32_t count = 0;
        for (uint16_t row_i : Columns[column].Rows)
        {
            if (Rows[row_i].Live > 0u) {
                ++count;
            }
        }
        return count;
    }

    void MarkResidualRows(uint16_t column)
    {
        for (uint16_t row_i : Columns[column].Rows)
        {
            if (Rows[row_i].Live > 0u && !DeferredRows[row_i])
            {
                DeferredRows[row_i] = 1;
                ++ResidualRows;
            }
        }
    }

    uint16_t SelectInactivationColumn() const
    {
        if (Codec.Solver == PeelSolver::RqccLowref) {
            return SelectDegree2ComponentColumn(0u);
        }
        return SelectDegree2ComponentColumn(Codec.SolverCandidateLimit);
    }

    uint16_t SelectFallbackColumn() const
    {
        uint16_t best = UINT16_MAX;
        uint32_t best_refs = UINT_MAX;
        for (uint32_t c = 0; c < BlockCount; ++c)
        {
            if (!Columns[c].Todo) {
                continue;
            }
            const uint32_t refs = CountLiveRowsForColumn((uint16_t)c);
            if (best == UINT16_MAX || refs < best_refs)
            {
                best = (uint16_t)c;
                best_refs = refs;
            }
        }
        return best;
    }

    uint16_t SelectDegree2ComponentColumn(uint16_t top_k) const
    {
        std::vector<uint16_t> candidates;
        candidates.reserve(BlockCount);
        for (uint32_t r = 0; r < Rows.size(); ++r)
        {
            if (Rows[r].Live != 2u) {
                continue;
            }
            for (uint16_t c : Rows[r].Columns)
            {
                if (Columns[c].Todo) {
                    candidates.push_back(c);
                }
            }
        }
        if (candidates.empty()) {
            return SelectFallbackColumn();
        }
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()),
            candidates.end());

        std::sort(candidates.begin(), candidates.end(),
            [this](uint16_t a, uint16_t b) {
                const uint32_t ar = CountLiveRowsForColumn(a);
                const uint32_t br = CountLiveRowsForColumn(b);
                if (ar != br) {
                    return ar < br;
                }
                return a < b;
            });

        if (top_k == 0u || top_k > candidates.size()) {
            top_k = (uint16_t)candidates.size();
        }

        uint16_t best = candidates[0];
        uint32_t best_boundary = 0;
        for (uint16_t i = 0; i < top_k; ++i)
        {
            const uint16_t c = candidates[i];
            const uint32_t boundary = CountBoundaryRows(c);
            if (i == 0u || boundary > best_boundary)
            {
                best = c;
                best_boundary = boundary;
            }
        }
        return best;
    }

    uint32_t CountBoundaryRows(uint16_t column) const
    {
        uint32_t boundary = 0;
        for (uint16_t row_i : Columns[column].Rows)
        {
            const RowState& row = Rows[row_i];
            if (row.Live > 2u) {
                ++boundary;
            }
        }
        return boundary;
    }

    const PeelingCodec& Codec;
    uint32_t BlockCount;
    std::vector<RowState> Rows;
    std::vector<ColumnState> Columns;
    std::vector<uint8_t> DeferredRows;
    std::queue<uint16_t> Queue;
    uint32_t TodoColumns = 0;
    uint32_t ResidualRows = 0;
    uint32_t ResidualColumns = 0;
    uint64_t MatrixRefs = 0;
    uint64_t MatrixXors = 0;
};

} // namespace

std::vector<std::vector<uint16_t> > GeneratePeelMatrixRows(
    const PeelingCodec& codec,
    uint32_t block_count,
    uint32_t row_count,
    uint64_t seed)
{
    Rng rng(seed);
    std::vector<std::vector<uint16_t> > rows;
    rows.reserve(row_count);
    for (uint32_t r = 0; r < row_count; ++r)
    {
        const uint32_t degree = ChooseDegree(codec, block_count, rng);
        rows.push_back(RandomRowColumns(block_count, degree, rng));
    }
    return rows;
}

PeelEvaluation EvaluatePeeling(
    const PeelingCodec& codec,
    uint32_t block_count,
    uint64_t seed)
{
    const std::vector<std::vector<uint16_t> > rows =
        GeneratePeelMatrixRows(codec, block_count, block_count, seed);
    return EvaluatePeelingRows(codec, block_count, rows);
}

PeelEvaluation EvaluatePeelingRows(
    const PeelingCodec& codec,
    uint32_t block_count,
    const std::vector<std::vector<uint16_t> >& rows)
{
    PeelSolverState solver(codec, block_count, rows);
    uint64_t refs = 0;
    uint64_t xors = 0;
    for (size_t i = 0; i < rows.size(); ++i)
    {
        refs += rows[i].size();
        if (!rows[i].empty()) {
            xors += rows[i].size() - 1u;
        }
    }
    solver.SetMatrixWork(refs, xors);
    return solver.Run();
}

} // namespace wirehair_v2
