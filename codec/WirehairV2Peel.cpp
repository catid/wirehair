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

class DegreeSampler
{
public:
    DegreeSampler(
        const PeelingCodec& codec,
        uint32_t block_count)
        : MinDegree(ClampDegree(codec.MinDegree, block_count)),
          MaxDegree(ClampDegree(codec.MaxDegree, block_count))
    {
        Cumulative.reserve(MaxDegree - MinDegree + 1u);
        for (uint32_t d = MinDegree; d <= MaxDegree; ++d)
        {
            Total += DegreeWeight(codec, d, block_count);
            Cumulative.push_back(Total);
        }
        if (Total <= 0.0) {
            Cumulative.clear();
        }
    }

    uint32_t Sample(Rng& rng) const
    {
        if (Cumulative.empty()) {
            return MinDegree;
        }

        const double target = rng.Unit() * Total;
        for (size_t i = 0; i < Cumulative.size(); ++i)
        {
            if (target <= Cumulative[i]) {
                return MinDegree + (uint32_t)i;
            }
        }
        return MaxDegree;
    }

private:
    uint32_t MinDegree;
    uint32_t MaxDegree;
    double Total = 0.0;
    std::vector<double> Cumulative;
};

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
    uint32_t Weight2Refs;
    bool Todo;
};

struct ColumnCandidate
{
    uint16_t Column;
    uint32_t Refs;
    uint32_t Degree2Refs;
    uint32_t Weight2Refs;
    uint32_t Lookahead;
    uint32_t Fill;
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
            Columns[c].Weight2Refs = 0u;
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
            else if (Rows[r].Live == 2u) {
                IncrementWeight2Refs(Rows[r]);
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
            else if (row.Live == 2u) {
                IncrementWeight2Refs(row);
            }
        }
    }

    uint32_t CountLiveRowsForColumn(uint16_t column) const
    {
        // For a todo column, every incident row still has at least this
        // column live, so the live-row count is just the stored incidence.
        return (uint32_t)Columns[column].Rows.size();
    }

    unsigned ScanTodoColumns(
        const RowState& row,
        uint16_t* out,
        unsigned out_count) const
    {
        unsigned count = 0u;
        for (uint16_t c : row.Columns)
        {
            if (Columns[c].Todo)
            {
                if (count < out_count) {
                    out[count] = c;
                }
                ++count;
            }
        }
        return count;
    }

    void IncrementWeight2Refs(const RowState& row)
    {
        uint16_t live[2] = {UINT16_MAX, UINT16_MAX};
        const unsigned count = ScanTodoColumns(row, live, 2u);
        if (count != 2u) {
            return;
        }
        ++Columns[live[0]].Weight2Refs;
        ++Columns[live[1]].Weight2Refs;
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
            return SelectRqccLowRefColumn();
        }
        return SelectKsBoundaryTopKColumn(Codec.SolverCandidateLimit);
    }

    bool LowRefBetter(
        const ColumnCandidate& a,
        const ColumnCandidate& b) const
    {
        if (a.Refs != b.Refs) {
            return a.Refs < b.Refs;
        }
        if (a.Weight2Refs != b.Weight2Refs) {
            return a.Weight2Refs > b.Weight2Refs;
        }
        return a.Column > b.Column;
    }

    bool DefaultBetter(
        const ColumnCandidate& a,
        const ColumnCandidate& b) const
    {
        if (a.Weight2Refs != b.Weight2Refs) {
            return a.Weight2Refs > b.Weight2Refs;
        }
        if (a.Refs != b.Refs) {
            return a.Refs > b.Refs;
        }
        return a.Column > b.Column;
    }

    bool BoundaryBetter(
        const ColumnCandidate& a,
        const ColumnCandidate& b) const
    {
        const uint32_t a_boundary =
            a.Refs > a.Degree2Refs ? a.Refs - a.Degree2Refs : 0u;
        const uint32_t b_boundary =
            b.Refs > b.Degree2Refs ? b.Refs - b.Degree2Refs : 0u;
        if (a_boundary != b_boundary) {
            return a_boundary > b_boundary;
        }
        if (a.Degree2Refs != b.Degree2Refs) {
            return a.Degree2Refs > b.Degree2Refs;
        }
        if (a.Lookahead != b.Lookahead) {
            return a.Lookahead > b.Lookahead;
        }
        if (a.Fill != b.Fill) {
            return a.Fill < b.Fill;
        }
        return a.Column > b.Column;
    }

    ColumnCandidate MeasureCandidate(uint16_t column) const
    {
        ColumnCandidate candidate;
        candidate.Column = column;
        candidate.Refs = CountLiveRowsForColumn(column);
        candidate.Degree2Refs = 0u;
        candidate.Weight2Refs = Columns[column].Weight2Refs;
        candidate.Lookahead = 0u;
        candidate.Fill = 0u;

        for (uint16_t row_i : Columns[column].Rows)
        {
            const RowState& row = Rows[row_i];
            if (row.Live == 0u) {
                continue;
            }
            if (row.Live > 1u) {
                candidate.Fill += row.Live - 1u;
            }
            if (row.Live == 2u)
            {
                uint16_t live[2] = {UINT16_MAX, UINT16_MAX};
                const unsigned count = ScanTodoColumns(row, live, 2u);
                if (count != 2u) {
                    continue;
                }
                uint16_t partner = UINT16_MAX;
                if (live[0] == column) {
                    partner = live[1];
                }
                else if (live[1] == column) {
                    partner = live[0];
                }
                else {
                    continue;
                }
                ++candidate.Degree2Refs;
                candidate.Lookahead += Columns[partner].Weight2Refs;
            }
        }
        return candidate;
    }

    uint16_t SelectLowRefFallbackColumn() const
    {
        bool have_best = false;
        ColumnCandidate best = {};
        for (uint32_t c = 0; c < BlockCount; ++c)
        {
            if (!Columns[c].Todo) {
                continue;
            }
            const ColumnCandidate candidate = MeasureCandidate((uint16_t)c);
            if (!have_best || LowRefBetter(candidate, best))
            {
                best = candidate;
                have_best = true;
            }
        }
        return have_best ? best.Column : UINT16_MAX;
    }

    uint16_t SelectMinRowLowRefColumn() const
    {
        uint16_t best_row = UINT16_MAX;
        uint16_t best_live = UINT16_MAX;
        for (uint32_t r = 0; r < Rows.size(); ++r)
        {
            const RowState& row = Rows[r];
            if (row.Live <= 1u || row.Live >= best_live) {
                continue;
            }
            best_row = (uint16_t)r;
            best_live = row.Live;
        }
        if (best_row == UINT16_MAX) {
            return SelectLowRefFallbackColumn();
        }

        bool have_best = false;
        ColumnCandidate best = {};
        for (uint16_t column : Rows[best_row].Columns)
        {
            if (!Columns[column].Todo) {
                continue;
            }
            const ColumnCandidate candidate = MeasureCandidate(column);
            if (!have_best || LowRefBetter(candidate, best))
            {
                best = candidate;
                have_best = true;
            }
        }
        return have_best ? best.Column : SelectLowRefFallbackColumn();
    }

    uint16_t FindRoot(
        std::vector<uint16_t>& parent,
        uint16_t column) const
    {
        while (parent[column] != column)
        {
            parent[column] = parent[parent[column]];
            column = parent[column];
        }
        return column;
    }

    void UnionRoots(
        std::vector<uint16_t>& parent,
        std::vector<uint32_t>& size,
        uint16_t a,
        uint16_t b) const
    {
        uint16_t ra = FindRoot(parent, a);
        uint16_t rb = FindRoot(parent, b);
        if (ra == rb) {
            return;
        }
        if (size[ra] < size[rb])
        {
            const uint16_t temp = ra;
            ra = rb;
            rb = temp;
        }
        parent[rb] = ra;
        size[ra] += size[rb];
    }

    void AddCandidate(
        std::vector<uint16_t>& candidates,
        std::vector<uint8_t>& seen,
        uint16_t column) const
    {
        if (column == UINT16_MAX || !Columns[column].Todo) {
            return;
        }
        if (!seen[column])
        {
            seen[column] = 1u;
            candidates.push_back(column);
        }
    }

    std::vector<uint16_t> CollectLargestDegree2ComponentColumns() const
    {
        std::vector<uint16_t> candidates;
        std::vector<uint16_t> parent(BlockCount);
        std::vector<uint32_t> size(BlockCount, 0u);
        std::vector<uint8_t> active(BlockCount, 0u);
        std::vector<uint16_t> active_columns;
        const auto activate = [&](uint16_t column) {
            if (!active[column])
            {
                active[column] = 1u;
                parent[column] = column;
                size[column] = 1u;
                active_columns.push_back(column);
            }
        };

        for (const RowState& row : Rows)
        {
            if (row.Live != 2u) {
                continue;
            }
            uint16_t live[2] = {UINT16_MAX, UINT16_MAX};
            const unsigned count = ScanTodoColumns(row, live, 2u);
            if (count == 2u)
            {
                activate(live[0]);
                activate(live[1]);
                UnionRoots(parent, size, live[0], live[1]);
            }
        }

        uint16_t best_root = UINT16_MAX;
        uint32_t best_size = 0u;
        for (uint16_t column : active_columns)
        {
            const uint16_t root = FindRoot(parent, column);
            if (size[root] > best_size)
            {
                best_size = size[root];
                best_root = root;
            }
        }
        if (best_root == UINT16_MAX) {
            return candidates;
        }

        std::vector<uint8_t> seen(BlockCount, 0u);
        for (const RowState& row : Rows)
        {
            if (row.Live != 2u) {
                continue;
            }
            uint16_t live[2] = {UINT16_MAX, UINT16_MAX};
            const unsigned count = ScanTodoColumns(row, live, 2u);
            if (count != 2u || FindRoot(parent, live[0]) != best_root) {
                continue;
            }
            AddCandidate(candidates, seen, live[0]);
            AddCandidate(candidates, seen, live[1]);
        }
        return candidates;
    }

    uint16_t SelectRqccLowRefColumn() const
    {
        const std::vector<uint16_t> columns =
            CollectLargestDegree2ComponentColumns();
        if (columns.empty()) {
            return SelectLowRefFallbackColumn();
        }

        bool have_best = false;
        ColumnCandidate best = {};
        for (uint16_t column : columns)
        {
            const ColumnCandidate candidate = MeasureCandidate(column);
            if (!have_best || LowRefBetter(candidate, best))
            {
                best = candidate;
                have_best = true;
            }
        }
        return have_best ? best.Column : SelectLowRefFallbackColumn();
    }

    uint16_t SelectKsBoundaryTopKColumn(uint16_t top_k) const
    {
        const std::vector<uint16_t> columns =
            CollectLargestDegree2ComponentColumns();
        if (columns.empty()) {
            return SelectMinRowLowRefColumn();
        }

        std::vector<ColumnCandidate> candidates;
        candidates.reserve(columns.size());
        for (uint16_t column : columns) {
            candidates.push_back(MeasureCandidate(column));
        }

        if (top_k == 0u || top_k > candidates.size()) {
            top_k = (uint16_t)candidates.size();
        }
        if (top_k < candidates.size())
        {
            std::sort(candidates.begin(), candidates.end(),
                [this](const ColumnCandidate& a, const ColumnCandidate& b) {
                    return DefaultBetter(a, b);
                });
            candidates.resize(top_k);
        }

        bool have_best = false;
        ColumnCandidate best = {};
        for (const ColumnCandidate& candidate : candidates)
        {
            if (!have_best || BoundaryBetter(candidate, best))
            {
                best = candidate;
                have_best = true;
            }
        }
        return have_best ? best.Column : SelectMinRowLowRefColumn();
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
    const DegreeSampler degrees(codec, block_count);
    std::vector<std::vector<uint16_t> > rows;
    rows.reserve(row_count);
    for (uint32_t r = 0; r < row_count; ++r)
    {
        const uint32_t degree = degrees.Sample(rng);
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
