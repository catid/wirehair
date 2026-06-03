// Standalone Wirehair peeling-schedule sweep harness.
//
// This intentionally does not include or modify WirehairCodec.cpp.  It reuses
// WirehairTools row generation and models only the peel graph/schedule.

#include "WirehairTools.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

using namespace wirehair;

namespace {

using Clock = std::chrono::steady_clock;

static const uint16_t kListTerm = 0xffffu;
static const uint8_t kMarkTodo = 0;
static const uint8_t kMarkPeel = 1;
static const uint8_t kMarkDefer = 2;

static double now_sec()
{
    return std::chrono::duration<double>(Clock::now().time_since_epoch()).count();
}

struct Rng
{
    uint64_t State;

    explicit Rng(uint64_t seed) : State(seed) {}

    uint64_t next()
    {
        uint64_t z = (State += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }

    uint32_t u32()
    {
        return (uint32_t)(next() >> 32);
    }

    double unit()
    {
        return (next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

static uint32_t hash32(uint32_t x)
{
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

static int default_threads()
{
    int hw = (int)std::thread::hardware_concurrency();
    if (hw <= 0) {
        hw = 8;
    }
    int t = (hw * 3) / 4;
    if (t < 1) {
        t = 1;
    }
    if (t > 104) {
        t = 104;
    }
    return t;
}

static std::vector<int> parse_int_list(const std::string& text)
{
    std::vector<int> out;
    size_t start = 0;
    while (start < text.size())
    {
        size_t comma = text.find(',', start);
        const std::string part = text.substr(
            start,
            comma == std::string::npos ? std::string::npos : comma - start);
        if (!part.empty()) {
            out.push_back(std::atoi(part.c_str()));
        }
        if (comma == std::string::npos) {
            break;
        }
        start = comma + 1;
    }
    return out;
}

struct Summary
{
    double Mean = 0.0;
    double Sd = 0.0;
    unsigned Min = 0;
    unsigned P50 = 0;
    unsigned P95 = 0;
    unsigned P99 = 0;
    unsigned Max = 0;
};

static unsigned percentile(const std::vector<unsigned>& sorted, double p)
{
    if (sorted.empty()) {
        return 0;
    }
    size_t index = (size_t)(p * (double)(sorted.size() - 1) + 0.5);
    if (index >= sorted.size()) {
        index = sorted.size() - 1;
    }
    return sorted[index];
}

static Summary summarize(std::vector<unsigned> values)
{
    Summary s;
    if (values.empty()) {
        return s;
    }

    std::sort(values.begin(), values.end());
    uint64_t sum = 0;
    for (unsigned v : values) {
        sum += v;
    }
    s.Mean = (double)sum / values.size();

    double var = 0.0;
    for (unsigned v : values)
    {
        const double d = (double)v - s.Mean;
        var += d * d;
    }
    s.Sd = std::sqrt(var / values.size());
    s.Min = values.front();
    s.P50 = percentile(values, 0.50);
    s.P95 = percentile(values, 0.95);
    s.P99 = percentile(values, 0.99);
    s.Max = values.back();
    return s;
}

struct Method
{
    int Id;
    const char* Name;
    const char* Cost;
};

static const Method kMethods[] = {
    {0,  "default",        "O(cols) cached degree-2 refs"},
    {1,  "exact_d2",       "O(cols*rowrefs) exact degree-2"},
    {2,  "lookahead",      "O(cols*rowrefs) degree-2 partner potential"},
    {3,  "weighted",       "O(cols*rowrefs) weighted exact/look/fill"},
    {4,  "min_fill",       "O(cols*rowrefs) Markowitz-style fill"},
    {5,  "min_fillsq",     "O(cols*rowrefs) quadratic fill"},
    {6,  "dup_partner",    "O(cols*rowrefs) repeated degree-2 partners"},
    {7,  "random_tie",     "O(cols) seeded randomized default ties"},
    {8,  "rowrefs",        "O(cols) column reference count"},
    {9,  "min_row_degree", "O(rows*degree) minimum live row"},
    {10, "raptorq_d2cc",   "O(rows*degree + cols) largest degree-2 component"},
    {11, "topk8_lookfill", "O(cols + 8*rowrefs) top-K look/fill"},
    {12, "exact_ratio",    "O(cols*rowrefs) exact degree-2 per fill"},
    {13, "min_maxrow",     "O(cols*rowrefs) minimize worst touched row"},
    {14, "d3_support",     "O(cols*rowrefs) degree-3 support"},
    {15, "rq_minrow",      "O(rows*degree) Raptor-style minimum row degree"},
    {16, "topk4_lookfill", "O(cols + 4*rowrefs) top-K look/fill"},
    {17, "topk16_lookfill","O(cols + 16*rowrefs) top-K look/fill"},
    {18, "rqd2_default",   "O(rows*degree + cols) largest degree-2 component/default tie"},
    {19, "rqd2_minfill",   "O(rows*degree + cols + rowrefs) largest degree-2 component/min-fill tie"},
    {20, "minrow_best",    "O(rows*degree) best column among all minimum rows"},
};

static const Method* find_method(int id)
{
    for (unsigned i = 0; i < sizeof(kMethods) / sizeof(kMethods[0]); ++i) {
        if (kMethods[i].Id == id) {
            return &kMethods[i];
        }
    }
    return nullptr;
}

struct Row
{
    std::vector<uint16_t> Columns;
    uint16_t Live = 0;
    uint16_t PeelColumn = kListTerm;
    bool Solved = false;
    bool Deferred = false;
};

struct Column
{
    std::vector<uint16_t> Rows;
    uint16_t Weight2Refs = 0;
    uint8_t Mark = kMarkTodo;
};

struct Metrics
{
    unsigned ExactD2 = 0;
    unsigned Fill = 0;
    unsigned FillSquare = 0;
    unsigned Lookahead = 0;
    unsigned DistinctPartners = 0;
    unsigned DuplicatePartners = 0;
    unsigned RowRefs = 0;
    unsigned MinLive = 0xffffu;
    unsigned MaxLive = 0;
    unsigned Degree3Rows = 0;
};

struct Key
{
    int64_t A = 0;
    int64_t B = 0;
    int64_t C = 0;
    int64_t D = 0;
    int64_t E = 0;

    Key() {}
    Key(int64_t a, int64_t b, int64_t c, int64_t d, int64_t e)
        : A(a), B(b), C(c), D(d), E(e)
    {
    }
};

static bool better_key(const Key& a, const Key& b)
{
    if (a.A != b.A) return a.A > b.A;
    if (a.B != b.B) return a.B > b.B;
    if (a.C != b.C) return a.C > b.C;
    if (a.D != b.D) return a.D > b.D;
    return a.E > b.E;
}

struct TrialResult
{
    unsigned ResidualCols = 0;
    unsigned ResidualRows = 0;
    unsigned Components = 0;
    unsigned MaxComponent = 0;
    unsigned SumSquares = 0;
    unsigned InitialTodo = 0;
    unsigned InitialD2Rows = 0;
    unsigned Choices = 0;
    unsigned BuildUsec = 0;
    unsigned GreedyUsec = 0;
    unsigned TotalUsec = 0;
};

class PeelSim
{
public:
    explicit PeelSim(unsigned block_count)
        : N(block_count),
          BlockNextPrime(NextPrime16((uint16_t)block_count)),
          PSeed(GetPeelSeed(block_count)),
          MixCount(GetDenseCount(block_count) + kHeavyRows),
          Columns(block_count)
    {
    }

    void feed(uint32_t row_seed)
    {
        Row row;
        PeelRowParameters params;
        params.Initialize(row_seed, PSeed, (uint16_t)N, MixCount);

        PeelRowIterator iter(params, (uint16_t)N, BlockNextPrime);
        uint16_t unmarked[2] = {kListTerm, kListTerm};
        uint16_t unmarked_count = 0;

        const uint16_t row_i = (uint16_t)Rows.size();
        do
        {
            const uint16_t column_i = iter.GetColumn();
            row.Columns.push_back(column_i);
            Columns[column_i].Rows.push_back(row_i);
            if (Columns[column_i].Mark == kMarkTodo) {
                unmarked[unmarked_count & 1] = column_i;
                ++unmarked_count;
            }
        } while (iter.Iterate());

        row.Live = unmarked_count;
        Rows.push_back(row);

        if (unmarked_count == 0) {
            add_deferred_row(row_i);
        }
        else if (unmarked_count == 1) {
            solve_with_peel(row_i, unmarked[0]);
        }
        else if (unmarked_count == 2)
        {
            ++Columns[unmarked[0]].Weight2Refs;
            ++Columns[unmarked[1]].Weight2Refs;
        }
    }

    TrialResult run_greedy(int method_id, uint32_t seed)
    {
        TrialResult result;
        result.InitialTodo = count_todo_columns();
        result.InitialD2Rows = count_live_rows(2);

        const double greedy_start = now_sec();
        for (;;)
        {
            const uint16_t column_i = select_column(method_id, seed);
            if (column_i == kListTerm) {
                break;
            }

            Column& column = Columns[column_i];
            if (column.Mark != kMarkTodo) {
                break;
            }
            column.Mark = kMarkDefer;
            ++DeferCount;
            ++ChoiceCount;
            avalanche(column_i);
        }
        result.GreedyUsec = elapsed_usec(greedy_start);

        result.ResidualCols = DeferCount;
        result.ResidualRows = DeferredRows;
        result.Choices = ChoiceCount;
        measure_components(result.Components, result.MaxComponent, result.SumSquares);
        return result;
    }

private:
    unsigned N;
    uint16_t BlockNextPrime;
    uint16_t PSeed;
    uint16_t MixCount;
    std::vector<Row> Rows;
    std::vector<Column> Columns;
    unsigned DeferCount = 0;
    unsigned DeferredRows = 0;
    unsigned ChoiceCount = 0;

    static unsigned elapsed_usec(double start)
    {
        const double usec = (now_sec() - start) * 1000000.0;
        return usec > 0.0 ? (unsigned)(usec + 0.5) : 0;
    }

    void add_deferred_row(uint16_t row_i)
    {
        Row& row = Rows[row_i];
        if (!row.Solved && !row.Deferred)
        {
            row.Deferred = true;
            ++DeferredRows;
        }
    }

    unsigned scan_todo_columns(const Row& row, uint16_t* out, unsigned out_count) const
    {
        unsigned count = 0;
        for (uint16_t column_i : row.Columns)
        {
            if (Columns[column_i].Mark == kMarkTodo)
            {
                if (count < out_count) {
                    out[count] = column_i;
                }
                ++count;
            }
        }
        return count;
    }

    void solve_with_peel(uint16_t row_i, uint16_t column_i)
    {
        if (column_i == kListTerm || Columns[column_i].Mark != kMarkTodo) {
            add_deferred_row(row_i);
            return;
        }

        Row& row = Rows[row_i];
        row.Solved = true;
        row.PeelColumn = column_i;
        Columns[column_i].Mark = kMarkPeel;
        avalanche(column_i);
    }

    void avalanche(uint16_t first_column)
    {
        std::vector<uint16_t> queue;
        queue.push_back(first_column);

        for (unsigned queue_i = 0; queue_i < queue.size(); ++queue_i)
        {
            const uint16_t column_i = queue[queue_i];
            const std::vector<uint16_t>& refs = Columns[column_i].Rows;

            for (uint16_t row_i : refs)
            {
                Row& row = Rows[row_i];
                if (row.Live == 0) {
                    continue;
                }

                --row.Live;
                if (row.Live == 1)
                {
                    uint16_t found[2] = {kListTerm, kListTerm};
                    const unsigned count = scan_todo_columns(row, found, 2);
                    row.Live = (uint16_t)count;
                    if (count == 1)
                    {
                        if (Columns[found[0]].Mark == kMarkTodo)
                        {
                            row.Solved = true;
                            row.PeelColumn = found[0];
                            Columns[found[0]].Mark = kMarkPeel;
                            queue.push_back(found[0]);
                        }
                    }
                    else if (count == 0) {
                        add_deferred_row(row_i);
                    }
                }
                else if (row.Live == 2)
                {
                    uint16_t found[3] = {kListTerm, kListTerm, kListTerm};
                    const unsigned count = scan_todo_columns(row, found, 3);
                    row.Live = (uint16_t)count;
                    if (count == 2)
                    {
                        ++Columns[found[0]].Weight2Refs;
                        ++Columns[found[1]].Weight2Refs;
                    }
                    else if (count == 1)
                    {
                        if (Columns[found[0]].Mark == kMarkTodo)
                        {
                            row.Solved = true;
                            row.PeelColumn = found[0];
                            Columns[found[0]].Mark = kMarkPeel;
                            queue.push_back(found[0]);
                        }
                    }
                    else if (count == 0) {
                        add_deferred_row(row_i);
                    }
                }
                else if (row.Live == 0) {
                    add_deferred_row(row_i);
                }
            }
        }
    }

    unsigned count_todo_columns() const
    {
        unsigned count = 0;
        for (const Column& column : Columns) {
            if (column.Mark == kMarkTodo) {
                ++count;
            }
        }
        return count;
    }

    unsigned count_live_rows(uint16_t live) const
    {
        unsigned count = 0;
        for (const Row& row : Rows) {
            if (row.Live == live) {
                ++count;
            }
        }
        return count;
    }

    bool row_has_column(const Row& row, uint16_t column_i) const
    {
        for (uint16_t c : row.Columns) {
            if (c == column_i) {
                return true;
            }
        }
        return false;
    }

    bool degree2_partner(const Row& row, uint16_t column_i, uint16_t& partner) const
    {
        if (row.Live != 2) {
            return false;
        }

        uint16_t found[2] = {kListTerm, kListTerm};
        unsigned count = 0;
        for (uint16_t c : row.Columns)
        {
            if (Columns[c].Mark == kMarkTodo)
            {
                if (count < 2) {
                    found[count] = c;
                }
                ++count;
            }
        }

        if (count != 2) {
            return false;
        }
        if (found[0] == column_i) {
            partner = found[1];
            return true;
        }
        if (found[1] == column_i) {
            partner = found[0];
            return true;
        }
        return false;
    }

    Metrics collect_metrics(uint16_t column_i) const
    {
        Metrics m;
        const Column& column = Columns[column_i];
        m.RowRefs = (unsigned)column.Rows.size();
        std::vector<uint16_t> partners;
        partners.reserve(8);

        for (uint16_t row_i : column.Rows)
        {
            const Row& row = Rows[row_i];
            if (row.Live == 0 || row.Deferred) {
                continue;
            }
            if (!row_has_column(row, column_i)) {
                continue;
            }

            if (row.Live > 1)
            {
                const unsigned links = row.Live - 1;
                m.Fill += links;
                m.FillSquare += links * links;
                if (row.Live < m.MinLive) {
                    m.MinLive = row.Live;
                }
                if (row.Live > m.MaxLive) {
                    m.MaxLive = row.Live;
                }
                if (row.Live == 3) {
                    ++m.Degree3Rows;
                }
            }

            uint16_t partner = kListTerm;
            if (degree2_partner(row, column_i, partner))
            {
                ++m.ExactD2;
                m.Lookahead += Columns[partner].Weight2Refs;

                bool seen = false;
                for (uint16_t p : partners)
                {
                    if (p == partner)
                    {
                        seen = true;
                        break;
                    }
                }
                if (seen) {
                    ++m.DuplicatePartners;
                }
                else {
                    partners.push_back(partner);
                    ++m.DistinctPartners;
                }
            }
        }

        if (m.MinLive == 0xffffu) {
            m.MinLive = 0;
        }
        return m;
    }

    Key metric_key(int method_id, uint16_t column_i, uint32_t seed) const
    {
        const Column& column = Columns[column_i];
        const Metrics m = collect_metrics(column_i);

        Key key;
        switch (method_id)
        {
        case 1:
            key = Key((int64_t)m.ExactD2, (int64_t)column.Weight2Refs,
                (int64_t)m.RowRefs, (int64_t)column_i, 0);
            break;
        case 2:
            key = Key((int64_t)m.ExactD2, (int64_t)m.Lookahead,
                (int64_t)column.Weight2Refs, (int64_t)m.RowRefs, (int64_t)column_i);
            break;
        case 3:
            key = Key((int64_t)m.ExactD2 * 1024 + (int64_t)m.Lookahead - (int64_t)m.Fill,
                (int64_t)column.Weight2Refs, (int64_t)m.RowRefs,
                -(int64_t)m.Fill, (int64_t)column_i);
            break;
        case 4:
            key = Key(-(int64_t)m.Fill, (int64_t)m.ExactD2,
                (int64_t)column.Weight2Refs, (int64_t)m.RowRefs, (int64_t)column_i);
            break;
        case 5:
            key = Key(-(int64_t)m.FillSquare, -(int64_t)m.Fill,
                (int64_t)m.ExactD2, (int64_t)m.RowRefs, (int64_t)column_i);
            break;
        case 6:
            key = Key((int64_t)m.DuplicatePartners, (int64_t)m.ExactD2,
                -(int64_t)m.DistinctPartners, (int64_t)m.RowRefs, (int64_t)column_i);
            break;
        case 12:
            key = Key((int64_t)((uint64_t)m.ExactD2 * 1000000u / (m.Fill + 1u)),
                (int64_t)m.ExactD2, (int64_t)m.Lookahead,
                (int64_t)column.Weight2Refs, (int64_t)column_i);
            break;
        case 13:
            key = Key(-(int64_t)m.MaxLive, -(int64_t)m.Fill,
                (int64_t)m.ExactD2, (int64_t)m.RowRefs, (int64_t)column_i);
            break;
        case 14:
            key = Key((int64_t)m.ExactD2, (int64_t)m.Degree3Rows,
                (int64_t)m.Lookahead, (int64_t)m.RowRefs, (int64_t)column_i);
            break;
        default:
            key = default_key(column_i, seed);
            break;
        }
        return key;
    }

    Key default_key(uint16_t column_i, uint32_t seed) const
    {
        const Column& column = Columns[column_i];
        return Key((int64_t)column.Weight2Refs, (int64_t)column.Rows.size(),
            (int64_t)column_i, (int64_t)seed, 0);
    }

    uint16_t select_by_full_scan(int method_id, uint32_t seed) const
    {
        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i = 0; column_i < N; ++column_i)
        {
            if (Columns[column_i].Mark != kMarkTodo) {
                continue;
            }

            Key key;
            if (method_id == 0) {
                key = default_key(column_i, seed);
            }
            else if (method_id == 7)
            {
                key = default_key(column_i, seed);
                key.C = (int64_t)hash32(seed ^ (uint32_t)column_i ^
                    ((uint32_t)ChoiceCount * 0x9e3779b9u));
                key.D = column_i;
            }
            else if (method_id == 8)
            {
                key = Key((int64_t)Columns[column_i].Rows.size(),
                    (int64_t)Columns[column_i].Weight2Refs,
                    (int64_t)column_i, 0, 0);
            }
            else {
                key = metric_key(method_id, column_i, seed);
            }

            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best;
    }

    std::vector<uint16_t> top_default_candidates(unsigned limit, uint32_t seed) const
    {
        std::vector<uint16_t> candidates;
        candidates.reserve(limit);
        for (uint16_t column_i = 0; column_i < N; ++column_i)
        {
            if (Columns[column_i].Mark != kMarkTodo) {
                continue;
            }

            const Key key = default_key(column_i, seed);
            unsigned pos = 0;
            while (pos < candidates.size() &&
                   !better_key(key, default_key(candidates[pos], seed))) {
                ++pos;
            }
            if (pos < limit)
            {
                candidates.insert(candidates.begin() + pos, column_i);
                if (candidates.size() > limit) {
                    candidates.pop_back();
                }
            }
        }
        return candidates;
    }

    uint16_t select_topk_lookfill(uint32_t seed, unsigned top_k) const
    {
        const std::vector<uint16_t> candidates = top_default_candidates(top_k, seed);
        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : candidates)
        {
            const Metrics m = collect_metrics(column_i);
            const Key key((int64_t)m.ExactD2, (int64_t)m.Lookahead,
                -(int64_t)m.Fill, (int64_t)Columns[column_i].Weight2Refs,
                (int64_t)column_i);
            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best;
    }

    uint16_t select_min_row_degree(uint32_t seed, bool raptor_style, bool all_min_rows) const
    {
        uint16_t best_row = kListTerm;
        unsigned best_live = 0xffffu;
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Live <= 1 || row.Deferred) {
                continue;
            }
            if (row.Live < best_live)
            {
                best_live = row.Live;
                best_row = row_i;
            }
        }

        if (best_row == kListTerm) {
            return select_by_full_scan(0, seed);
        }

        uint16_t best_column = kListTerm;
        Key best_key;
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row_i != best_row && (!all_min_rows || row.Live != best_live)) {
                continue;
            }
            for (uint16_t column_i : row.Columns)
            {
                if (Columns[column_i].Mark != kMarkTodo) {
                    continue;
                }

                Key key;
                if (raptor_style) {
                    key = Key(-(int64_t)Columns[column_i].Rows.size(),
                        (int64_t)Columns[column_i].Weight2Refs, (int64_t)column_i, 0, 0);
                }
                else {
                    key = metric_key(2, column_i, seed);
                }

                if (best_column == kListTerm || better_key(key, best_key))
                {
                    best_column = column_i;
                    best_key = key;
                }
            }
            if (!all_min_rows) {
                break;
            }
        }

        return best_column != kListTerm ? best_column : select_by_full_scan(0, seed);
    }

    uint16_t find_root(std::vector<uint16_t>& parent, uint16_t x) const
    {
        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void union_roots(std::vector<uint16_t>& parent, std::vector<uint16_t>& size,
        uint16_t a, uint16_t b) const
    {
        uint16_t ra = find_root(parent, a);
        uint16_t rb = find_root(parent, b);
        if (ra == rb) {
            return;
        }
        if (size[ra] < size[rb]) {
            const uint16_t temp = ra;
            ra = rb;
            rb = temp;
        }
        parent[rb] = ra;
        size[ra] = (uint16_t)(size[ra] + size[rb]);
    }

    uint16_t select_raptorq_d2_component(uint32_t seed, unsigned tie_mode) const
    {
        std::vector<uint16_t> parent(N);
        std::vector<uint16_t> size(N, 1);
        std::vector<uint8_t> active(N, 0);
        for (uint16_t i = 0; i < N; ++i) {
            parent[i] = i;
        }

        for (const Row& row : Rows)
        {
            if (row.Live != 2 || row.Deferred) {
                continue;
            }
            uint16_t found[2] = {kListTerm, kListTerm};
            unsigned count = 0;
            for (uint16_t c : row.Columns)
            {
                if (Columns[c].Mark == kMarkTodo)
                {
                    if (count < 2) {
                        found[count] = c;
                    }
                    ++count;
                }
            }
            if (count == 2)
            {
                active[found[0]] = 1;
                active[found[1]] = 1;
                union_roots(parent, size, found[0], found[1]);
            }
        }

        uint16_t best_root = kListTerm;
        unsigned best_size = 0;
        for (uint16_t column_i = 0; column_i < N; ++column_i)
        {
            if (!active[column_i]) {
                continue;
            }
            const uint16_t root = find_root(parent, column_i);
            if (size[root] > best_size)
            {
                best_size = size[root];
                best_root = root;
            }
        }
        if (best_root == kListTerm) {
            return select_min_row_degree(seed, true, false);
        }

        uint16_t best_column = kListTerm;
        Key best_key;
        for (const Row& row : Rows)
        {
            if (row.Live != 2 || row.Deferred) {
                continue;
            }
            uint16_t found[2] = {kListTerm, kListTerm};
            unsigned count = 0;
            for (uint16_t c : row.Columns)
            {
                if (Columns[c].Mark == kMarkTodo)
                {
                    if (count < 2) {
                        found[count] = c;
                    }
                    ++count;
                }
            }
            if (count != 2 || find_root(parent, found[0]) != best_root) {
                continue;
            }

            for (uint16_t column_i : found)
            {
                Key key;
                if (tie_mode == 1) {
                    key = default_key(column_i, seed);
                }
                else if (tie_mode == 2)
                {
                    const Metrics m = collect_metrics(column_i);
                    key = Key(-(int64_t)m.Fill, (int64_t)m.ExactD2,
                        (int64_t)Columns[column_i].Weight2Refs,
                        (int64_t)m.RowRefs, (int64_t)column_i);
                }
                else {
                    key = Key(-(int64_t)Columns[column_i].Rows.size(),
                        (int64_t)Columns[column_i].Weight2Refs, (int64_t)column_i, 0, 0);
                }
                if (best_column == kListTerm || better_key(key, best_key))
                {
                    best_column = column_i;
                    best_key = key;
                }
            }
        }

        return best_column != kListTerm ? best_column : select_by_full_scan(0, seed);
    }

    uint16_t select_column(int method_id, uint32_t seed) const
    {
        switch (method_id)
        {
        case 9:
            return select_min_row_degree(seed, false, false);
        case 10:
            return select_raptorq_d2_component(seed, 0);
        case 11:
            return select_topk_lookfill(seed, 8);
        case 15:
            return select_min_row_degree(seed, true, false);
        case 16:
            return select_topk_lookfill(seed, 4);
        case 17:
            return select_topk_lookfill(seed, 16);
        case 18:
            return select_raptorq_d2_component(seed, 1);
        case 19:
            return select_raptorq_d2_component(seed, 2);
        case 20:
            return select_min_row_degree(seed, false, true);
        default:
            return select_by_full_scan(method_id, seed);
        }
    }

    void measure_components(unsigned& component_count,
        unsigned& max_component,
        unsigned& sum_squares) const
    {
        component_count = 0;
        max_component = 0;
        sum_squares = 0;
        if (DeferCount == 0) {
            return;
        }

        std::vector<uint16_t> map(N, kListTerm);
        std::vector<uint16_t> parent(DeferCount);
        std::vector<uint16_t> size(DeferCount, 1);
        uint16_t next = 0;
        for (uint16_t column_i = 0; column_i < N; ++column_i)
        {
            if (Columns[column_i].Mark == kMarkDefer)
            {
                map[column_i] = next;
                parent[next] = next;
                ++next;
            }
        }

        for (const Row& row : Rows)
        {
            if (!row.Deferred) {
                continue;
            }

            uint16_t first = kListTerm;
            for (uint16_t column_i : row.Columns)
            {
                const uint16_t mapped = map[column_i];
                if (mapped == kListTerm) {
                    continue;
                }
                if (first == kListTerm) {
                    first = mapped;
                }
                else {
                    union_roots(parent, size, first, mapped);
                }
            }
        }

        for (uint16_t i = 0; i < DeferCount; ++i)
        {
            if (find_root(parent, i) != i) {
                continue;
            }
            ++component_count;
            if (size[i] > max_component) {
                max_component = size[i];
            }
            sum_squares += (unsigned)size[i] * (unsigned)size[i];
        }
    }
};

static std::vector<uint32_t> make_row_ids(
    unsigned N,
    unsigned rows,
    const std::string& source,
    double loss,
    uint64_t seed)
{
    std::vector<uint32_t> ids;
    ids.reserve(rows);
    Rng rng(seed);

    if (source == "systematic")
    {
        for (unsigned i = 0; i < rows; ++i) {
            ids.push_back(i);
        }
        return ids;
    }

    if (source == "repair")
    {
        uint32_t id = (uint32_t)(N + (rng.u32() % 16u));
        for (unsigned i = 0; i < rows; ++i) {
            ids.push_back(id++);
        }
        return ids;
    }

    uint32_t id = 0;
    while (ids.size() < rows)
    {
        const bool drop = rng.unit() < loss;
        if (!drop) {
            ids.push_back(id);
        }
        ++id;
    }
    return ids;
}

static TrialResult run_trial(
    int method_id,
    unsigned N,
    unsigned row_count,
    const std::string& source,
    double loss,
    uint64_t seed)
{
    const double start = now_sec();

    PeelSim sim(N);
    const std::vector<uint32_t> ids = make_row_ids(N, row_count, source, loss, seed);

    const double build_start = now_sec();
    for (uint32_t id : ids) {
        sim.feed(id);
    }
    const unsigned build_usec = (unsigned)((now_sec() - build_start) * 1000000.0 + 0.5);

    TrialResult result = sim.run_greedy(method_id, (uint32_t)seed);
    result.BuildUsec = build_usec;
    result.TotalUsec = (unsigned)((now_sec() - start) * 1000000.0 + 0.5);
    return result;
}

static std::vector<int> expand_methods(const std::string& spec)
{
    if (spec == "all")
    {
        std::vector<int> ids;
        for (const Method& method : kMethods) {
            ids.push_back(method.Id);
        }
        return ids;
    }
    return parse_int_list(spec);
}

static void print_methods()
{
    std::printf("# methods:\n");
    for (const Method& method : kMethods) {
        std::printf("#   %2d %-16s %s\n", method.Id, method.Name, method.Cost);
    }
}

} // namespace

int main(int argc, char** argv)
{
    std::string nlist = "128,2048,32000";
    std::string methods_spec = "all";
    std::string source = "loss";
    int trials = 40;
    int threads = default_threads();
    int overhead = 0;
    int rows_override = 0;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0x9a7e11a);
    bool list_methods = false;

    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--N") && i + 1 < argc) {
            nlist = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--methods") && i + 1 < argc) {
            methods_spec = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--source") && i + 1 < argc) {
            source = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            trials = std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) {
            threads = std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--loss") && i + 1 < argc) {
            loss = std::atof(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed = std::strtoull(argv[++i], nullptr, 0);
        }
        else if (!std::strcmp(argv[i], "--overhead") && i + 1 < argc) {
            overhead = std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--rows") && i + 1 < argc) {
            rows_override = std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--list-methods")) {
            list_methods = true;
        }
        else if (!std::strcmp(argv[i], "--help"))
        {
            std::printf("usage: peel_sweep [--N csv] [--methods all|csv] [--trials n]\n");
            std::printf("                  [--source loss|systematic|repair] [--loss p]\n");
            std::printf("                  [--overhead n|--rows n] [--threads n] [--seed x]\n");
            print_methods();
            return 0;
        }
    }

    if (list_methods)
    {
        print_methods();
        return 0;
    }

    if (trials < 1) {
        trials = 1;
    }
    if (threads < 1) {
        threads = 1;
    }
    if (threads > trials) {
        threads = trials;
    }
    if (source != "loss" && source != "systematic" && source != "repair")
    {
        std::fprintf(stderr, "unknown source: %s\n", source.c_str());
        return 1;
    }
    if (loss < 0.0 || loss >= 1.0)
    {
        std::fprintf(stderr, "--loss must be in [0, 1)\n");
        return 1;
    }

    const std::vector<int> Ns = parse_int_list(nlist);
    const std::vector<int> method_ids = expand_methods(methods_spec);
    if (Ns.empty() || method_ids.empty())
    {
        std::fprintf(stderr, "empty --N or --methods list\n");
        return 1;
    }

    print_methods();
    std::printf("# sweep: trials=%d threads=%d source=%s loss=%.3f overhead=%d rows_override=%d seed=0x%llx\n",
        trials, threads, source.c_str(), loss, overhead, rows_override,
        (unsigned long long)seed);
    std::printf("method_id,method,cost,N,rows,source,trials,"
                "cols_mu,cols_sd,c50,c95,c99,cmax,"
                "rows_mu,rows_sd,r50,r95,r99,rmax,"
                "comp_mu,maxcomp_mu,sumsq_mu,"
                "init_todo_mu,init_d2_mu,choices_mu,"
                "build_us,greedy_us,total_us\n");

    for (int method_id : method_ids)
    {
        const Method* method = find_method(method_id);
        if (!method)
        {
            std::fprintf(stderr, "unknown method id %d\n", method_id);
            return 1;
        }

        for (int n_value : Ns)
        {
            if (n_value < CAT_WIREHAIR_MIN_N || n_value > CAT_WIREHAIR_MAX_N)
            {
                std::fprintf(stderr, "N out of range: %d\n", n_value);
                return 1;
            }
            const unsigned N = (unsigned)n_value;
            const int row_count_signed = rows_override > 0 ?
                rows_override : n_value + overhead;
            if (row_count_signed < 0)
            {
                std::fprintf(stderr, "row count is negative for N=%d\n", n_value);
                return 1;
            }
            const unsigned row_count = (unsigned)row_count_signed;
            if (row_count > 65535u)
            {
                std::fprintf(stderr, "row count %u exceeds 16-bit experiment model limit\n",
                    row_count);
                return 1;
            }

            std::vector<TrialResult> results((size_t)trials);
            std::atomic<int> next{0};
            std::vector<std::thread> workers;
            for (int t = 0; t < threads; ++t)
            {
                workers.emplace_back([&]() {
                    for (;;)
                    {
                        const int trial = next.fetch_add(1);
                        if (trial >= trials) {
                            break;
                        }
                        const uint64_t trial_seed = seed ^
                            ((uint64_t)method_id * UINT64_C(0x9e3779b97f4a7c15)) ^
                            ((uint64_t)N * UINT64_C(0xbf58476d1ce4e5b9)) ^
                            ((uint64_t)trial * UINT64_C(0x94d049bb133111eb));
                        results[(size_t)trial] = run_trial(
                            method_id, N, row_count, source, loss, trial_seed);
                    }
                });
            }
            for (std::thread& worker : workers) {
                worker.join();
            }

            std::vector<unsigned> cols, rows, comps, maxcomps, sumsqs;
            std::vector<unsigned> init_todo, init_d2, choices;
            unsigned long long build_sum = 0;
            unsigned long long greedy_sum = 0;
            unsigned long long total_sum = 0;
            cols.reserve(results.size());
            rows.reserve(results.size());
            comps.reserve(results.size());
            maxcomps.reserve(results.size());
            sumsqs.reserve(results.size());
            init_todo.reserve(results.size());
            init_d2.reserve(results.size());
            choices.reserve(results.size());

            for (const TrialResult& r : results)
            {
                cols.push_back(r.ResidualCols);
                rows.push_back(r.ResidualRows);
                comps.push_back(r.Components);
                maxcomps.push_back(r.MaxComponent);
                sumsqs.push_back(r.SumSquares);
                init_todo.push_back(r.InitialTodo);
                init_d2.push_back(r.InitialD2Rows);
                choices.push_back(r.Choices);
                build_sum += r.BuildUsec;
                greedy_sum += r.GreedyUsec;
                total_sum += r.TotalUsec;
            }

            const Summary col_s = summarize(cols);
            const Summary row_s = summarize(rows);
            const Summary comp_s = summarize(comps);
            const Summary maxcomp_s = summarize(maxcomps);
            const Summary sumsq_s = summarize(sumsqs);
            const Summary init_todo_s = summarize(init_todo);
            const Summary init_d2_s = summarize(init_d2);
            const Summary choice_s = summarize(choices);

            std::printf("%d,%s,\"%s\",%u,%u,%s,%d,"
                        "%.2f,%.2f,%u,%u,%u,%u,"
                        "%.2f,%.2f,%u,%u,%u,%u,"
                        "%.2f,%.2f,%.2f,"
                        "%.2f,%.2f,%.2f,"
                        "%.1f,%.1f,%.1f\n",
                method->Id, method->Name, method->Cost, N, row_count,
                source.c_str(), trials,
                col_s.Mean, col_s.Sd, col_s.P50, col_s.P95, col_s.P99, col_s.Max,
                row_s.Mean, row_s.Sd, row_s.P50, row_s.P95, row_s.P99, row_s.Max,
                comp_s.Mean, maxcomp_s.Mean, sumsq_s.Mean,
                init_todo_s.Mean, init_d2_s.Mean, choice_s.Mean,
                (double)build_sum / trials,
                (double)greedy_sum / trials,
                (double)total_sum / trials);
            std::fflush(stdout);
        }
    }

    return 0;
}
