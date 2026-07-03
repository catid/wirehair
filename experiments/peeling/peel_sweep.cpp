// Standalone Wirehair peeling-schedule sweep harness.
//
// This intentionally does not include or modify WirehairCodec.cpp.  It reuses
// WirehairTools row generation and models only the peel graph/schedule.

#include "WirehairTools.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <utility>
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

static uint64_t hash_string64(const char* text)
{
    uint64_t h = UINT64_C(1469598103934665603);
    while (*text)
    {
        h ^= (uint8_t)*text++;
        h *= UINT64_C(1099511628211);
    }
    return h;
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

static std::vector<int> parse_int_list(const std::string& text, bool* ok = nullptr)
{
    if (ok) {
        *ok = true;
    }
    std::vector<int> out;
    size_t start = 0;
    while (start < text.size())
    {
        size_t comma = text.find(',', start);
        const std::string part = text.substr(
            start,
            comma == std::string::npos ? std::string::npos : comma - start);
        if (part.empty())
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }

        char* end = nullptr;
        const long value = std::strtol(part.c_str(), &end, 10);
        if (end == part.c_str() || *end != '\0' ||
            value < INT_MIN || value > INT_MAX)
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }

        out.push_back((int)value);
        if (comma == std::string::npos) {
            break;
        }
        if (comma + 1 == text.size())
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }
        start = comma + 1;
    }
    return out;
}

static bool parse_int_scalar(const char* text, int& out)
{
    bool ok = false;
    const std::vector<int> values = parse_int_list(text ? text : "", &ok);
    if (!ok || values.size() != 1u) {
        return false;
    }
    out = values[0];
    return true;
}

static bool parse_u64_scalar(const char* text, uint64_t& out)
{
    if (!text || !*text || *text < '0' || *text > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 0);
    if (errno != 0 || end == text || !end || *end != '\0') {
        return false;
    }
    out = (uint64_t)value;
    return true;
}

static std::vector<double> parse_double_list(
    const std::string& text,
    bool* ok = nullptr)
{
    if (ok) {
        *ok = true;
    }
    std::vector<double> out;
    size_t start = 0;
    while (start < text.size())
    {
        size_t comma = text.find(',', start);
        const std::string part = text.substr(
            start,
            comma == std::string::npos ? std::string::npos : comma - start);
        if (part.empty())
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }

        char* end = nullptr;
        const double value = std::strtod(part.c_str(), &end);
        if (end == part.c_str() || *end != '\0' || !std::isfinite(value))
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }

        out.push_back(value);
        if (comma == std::string::npos) {
            break;
        }
        if (comma + 1 == text.size())
        {
            if (ok) {
                *ok = false;
            }
            return out;
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

struct RealSummary
{
    double Mean = 0.0;
    double Sd = 0.0;
    double Min = 0.0;
    double P50 = 0.0;
    double P95 = 0.0;
    double P99 = 0.0;
    double Max = 0.0;
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

static double percentile_real(const std::vector<double>& sorted, double p)
{
    if (sorted.empty()) {
        return 0.0;
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

static RealSummary summarize_real(std::vector<double> values)
{
    RealSummary s;
    if (values.empty()) {
        return s;
    }

    std::sort(values.begin(), values.end());
    double sum = 0.0;
    for (double v : values) {
        sum += v;
    }
    s.Mean = sum / values.size();

    double var = 0.0;
    for (double v : values)
    {
        const double d = v - s.Mean;
        var += d * d;
    }
    s.Sd = std::sqrt(var / values.size());
    s.Min = values.front();
    s.P50 = percentile_real(values, 0.50);
    s.P95 = percentile_real(values, 0.95);
    s.P99 = percentile_real(values, 0.99);
    s.Max = values.back();
    return s;
}

enum PoolKind
{
    kPoolLegacy = 0,
    kPoolAll,
    kPoolTopDefault,
    kPoolD2Rows,
    kPoolD2LargestComponent,
    kPoolMinRow,
    kPoolAllMinRows,
};

enum ScoreKind
{
    kScoreLegacy = 0,
    kScoreDefault,
    kScoreExactD2,
    kScoreLookFill,
    kScoreRatio,
    kScoreMinFill,
    kScoreMinFillSquare,
    kScoreDuplicate,
    kScoreD3Support,
    kScoreLowRef,
    kScoreLiveRefMin,
    kScoreLiveRefMax,
    kScoreMaxRowMin,
    kScoreRandom,
};

struct Method
{
    int Id;
    const char* Name;
    const char* Cost;
    PoolKind Pool;
    ScoreKind Score;
    unsigned Arg;
};

static const Method kMethods[] = {
    {0,  "default",        "O(cols) cached degree-2 refs", kPoolLegacy, kScoreLegacy, 0},
    {1,  "exact_d2",       "O(cols + touched_rows*degree) incremental exact degree-2", kPoolLegacy, kScoreLegacy, 0},
    {2,  "lookahead",      "O(cols + touched_rows*degree) incremental degree-2 partner potential", kPoolLegacy, kScoreLegacy, 0},
    {3,  "weighted",       "O(cols + touched_rows*degree) incremental exact/look/fill", kPoolLegacy, kScoreLegacy, 0},
    {4,  "min_fill",       "O(cols + touched_rows*degree) incremental Markowitz-style fill", kPoolLegacy, kScoreLegacy, 0},
    {5,  "min_fillsq",     "O(cols + touched_rows*degree) incremental quadratic fill", kPoolLegacy, kScoreLegacy, 0},
    {6,  "dup_partner",    "O(cols + touched_rows*degree) incremental repeated degree-2 partners", kPoolLegacy, kScoreLegacy, 0},
    {7,  "random_tie",     "O(cols) seeded randomized default ties", kPoolLegacy, kScoreLegacy, 0},
    {8,  "rowrefs",        "O(cols) column reference count", kPoolLegacy, kScoreLegacy, 0},
    {9,  "min_row_degree", "O(rows*degree) minimum live row", kPoolLegacy, kScoreLegacy, 0},
    {10, "raptorq_d2cc",   "O(rows*degree + cols) largest degree-2 component", kPoolLegacy, kScoreLegacy, 0},
    {11, "topk8_lookfill", "O(cols + 8*rowrefs) top-K look/fill", kPoolLegacy, kScoreLegacy, 0},
    {12, "exact_ratio",    "O(cols + touched_rows*degree) incremental exact degree-2 per fill", kPoolLegacy, kScoreLegacy, 0},
    {13, "min_maxrow",     "O(cols + touched_rows*degree) incremental minimize worst touched row", kPoolLegacy, kScoreLegacy, 0},
    {14, "d3_support",     "O(cols + touched_rows*degree) incremental degree-3 support", kPoolLegacy, kScoreLegacy, 0},
    {15, "rq_minrow",      "O(rows*degree) Raptor-style minimum row degree", kPoolLegacy, kScoreLegacy, 0},
    {16, "topk4_lookfill", "O(cols + 4*rowrefs) top-K look/fill", kPoolLegacy, kScoreLegacy, 0},
    {17, "topk16_lookfill","O(cols + 16*rowrefs) top-K look/fill", kPoolLegacy, kScoreLegacy, 0},
    {18, "rqd2_default",   "O(rows*degree + cols) largest degree-2 component/default tie", kPoolLegacy, kScoreLegacy, 0},
    {19, "rqd2_minfill",   "O(rows*degree + cols + rowrefs) largest degree-2 component/min-fill tie", kPoolLegacy, kScoreLegacy, 0},
    {20, "minrow_best",    "O(rows*degree) best column among all minimum rows", kPoolLegacy, kScoreLegacy, 0},

    {21, "top32_lookfill", "O(cols + 32*rowrefs) top-K/default pool + look/fill score", kPoolTopDefault, kScoreLookFill, 32},
    {22, "top64_lookfill", "O(cols + 64*rowrefs) top-K/default pool + look/fill score", kPoolTopDefault, kScoreLookFill, 64},
    {23, "top8_ratio",     "O(cols + 8*rowrefs) top-K/default pool + d2/fill ratio", kPoolTopDefault, kScoreRatio, 8},
    {24, "top16_ratio",    "O(cols + 16*rowrefs) top-K/default pool + d2/fill ratio", kPoolTopDefault, kScoreRatio, 16},
    {25, "top32_ratio",    "O(cols + 32*rowrefs) top-K/default pool + d2/fill ratio", kPoolTopDefault, kScoreRatio, 32},
    {26, "top16_dup",      "O(cols + 16*rowrefs) top-K/default pool + duplicate partners", kPoolTopDefault, kScoreDuplicate, 16},
    {27, "top16_d3",       "O(cols + 16*rowrefs) top-K/default pool + degree-3 support", kPoolTopDefault, kScoreD3Support, 16},
    {28, "top16_livemin",  "O(cols + 16*rowrefs) top-K/default pool + low live refs", kPoolTopDefault, kScoreLiveRefMin, 16},
    {29, "top16_lowref",   "O(cols) top-K/default pool + low reference count", kPoolTopDefault, kScoreLowRef, 16},
    {30, "top16_livemax",  "O(cols + 16*rowrefs) top-K/default pool + high live refs", kPoolTopDefault, kScoreLiveRefMax, 16},
    {31, "top32_minfill",  "O(cols + 32*rowrefs) top-K/default pool + min-fill", kPoolTopDefault, kScoreMinFill, 32},
    {32, "top32_maxrow",   "O(cols + 32*rowrefs) top-K/default pool + min max-row", kPoolTopDefault, kScoreMaxRowMin, 32},

    {33, "d2pool_default", "O(rows*degree + candidates) all degree-2 endpoints + default score", kPoolD2Rows, kScoreDefault, 0},
    {34, "d2pool_lowref",  "O(rows*degree + candidates) all degree-2 endpoints + low refs", kPoolD2Rows, kScoreLowRef, 0},
    {35, "d2pool_lookfill","O(rows*degree + candidates*rowrefs) all degree-2 endpoints + look/fill", kPoolD2Rows, kScoreLookFill, 0},
    {36, "d2pool_dup",     "O(rows*degree + candidates*rowrefs) all degree-2 endpoints + duplicate partners", kPoolD2Rows, kScoreDuplicate, 0},
    {37, "d2pool_ratio",   "O(rows*degree + candidates*rowrefs) all degree-2 endpoints + d2/fill ratio", kPoolD2Rows, kScoreRatio, 0},
    {38, "d2pool_minfill", "O(rows*degree + candidates*rowrefs) all degree-2 endpoints + min-fill", kPoolD2Rows, kScoreMinFill, 0},
    {39, "d2pool_maxrow",  "O(rows*degree + candidates*rowrefs) all degree-2 endpoints + min max-row", kPoolD2Rows, kScoreMaxRowMin, 0},

    {40, "rqcc_lowref",    "O(rows*degree + cols) largest degree-2 component + low refs", kPoolD2LargestComponent, kScoreLowRef, 0},
    {41, "rqcc_lookfill",  "O(rows*degree + cols + cc*rowrefs) largest degree-2 component + look/fill", kPoolD2LargestComponent, kScoreLookFill, 0},
    {42, "rqcc_dup",       "O(rows*degree + cols + cc*rowrefs) largest degree-2 component + duplicate partners", kPoolD2LargestComponent, kScoreDuplicate, 0},
    {43, "rqcc_ratio",     "O(rows*degree + cols + cc*rowrefs) largest degree-2 component + d2/fill ratio", kPoolD2LargestComponent, kScoreRatio, 0},
    {44, "rqcc_minfill",   "O(rows*degree + cols + cc*rowrefs) largest degree-2 component + min-fill", kPoolD2LargestComponent, kScoreMinFill, 0},
    {45, "rqcc_livemin",   "O(rows*degree + cols + cc*rowrefs) largest degree-2 component + low live refs", kPoolD2LargestComponent, kScoreLiveRefMin, 0},
    {46, "rqcc_maxrow",    "O(rows*degree + cols + cc*rowrefs) largest degree-2 component + min max-row", kPoolD2LargestComponent, kScoreMaxRowMin, 0},

    {47, "minrow_default", "O(rows*degree + row_degree) first minimum row + default score", kPoolMinRow, kScoreDefault, 0},
    {48, "minrow_lowref",  "O(rows*degree + row_degree) first minimum row + low refs", kPoolMinRow, kScoreLowRef, 0},
    {49, "minrow_minfill", "O(rows*degree + row_degree*rowrefs) first minimum row + min-fill", kPoolMinRow, kScoreMinFill, 0},
    {50, "minrow_lookfill","O(rows*degree + row_degree*rowrefs) first minimum row + look/fill", kPoolMinRow, kScoreLookFill, 0},
    {51, "minrow_ratio",   "O(rows*degree + row_degree*rowrefs) first minimum row + d2/fill ratio", kPoolMinRow, kScoreRatio, 0},

    {52, "allmin_default", "O(rows*degree + candidates) all minimum rows + default score", kPoolAllMinRows, kScoreDefault, 0},
    {53, "allmin_lowref",  "O(rows*degree + candidates) all minimum rows + low refs", kPoolAllMinRows, kScoreLowRef, 0},
    {54, "allmin_minfill", "O(rows*degree + candidates*rowrefs) all minimum rows + min-fill", kPoolAllMinRows, kScoreMinFill, 0},
    {55, "allmin_lookfill","O(rows*degree + candidates*rowrefs) all minimum rows + look/fill", kPoolAllMinRows, kScoreLookFill, 0},
    {56, "allmin_ratio",   "O(rows*degree + candidates*rowrefs) all minimum rows + d2/fill ratio", kPoolAllMinRows, kScoreRatio, 0},

    {57, "global_lowref",  "O(cols) all columns + low reference count", kPoolAll, kScoreLowRef, 0},
    {58, "global_livemin", "O(cols + touched_rows*degree) all columns + low live refs", kPoolAll, kScoreLiveRefMin, 0},
    {59, "global_livemax", "O(cols*rowrefs) all columns + high live refs", kPoolAll, kScoreLiveRefMax, 0},
    {60, "global_random",  "O(cols) all columns + seeded randomized tie", kPoolAll, kScoreRandom, 0},
    {61, "top32_lowref",   "O(cols) top-K/default pool + low reference count", kPoolTopDefault, kScoreLowRef, 32},
    {62, "top64_ratio",    "O(cols + 64*rowrefs) top-K/default pool + d2/fill ratio", kPoolTopDefault, kScoreRatio, 64},
    {63, "top64_minfill",  "O(cols + 64*rowrefs) top-K/default pool + min-fill", kPoolTopDefault, kScoreMinFill, 64},
    {64, "top32_d3",       "O(cols + 32*rowrefs) top-K/default pool + degree-3 support", kPoolTopDefault, kScoreD3Support, 32},

    {65, "random_col",     "O(cols) uniformly random inactivation", kPoolLegacy, kScoreLegacy, 0},
    {66, "bm_minrand",     "O(rows*degree) min row, random all-but-one", kPoolLegacy, kScoreLegacy, 0},
    {67, "bm_minbest",     "O(rows*degree) min row, leave best default", kPoolLegacy, kScoreLegacy, 0},
    {68, "bm_minworst",    "O(rows*degree) min row, leave worst default", kPoolLegacy, kScoreLegacy, 0},
    {69, "bm_c_tiny",      "O(rows*degree) Burshtein-Miller omega2-heavy row", kPoolLegacy, kScoreLegacy, 0},
    {70, "bm_c_decay",     "O(rows*degree) Burshtein-Miller decaying row weights", kPoolLegacy, kScoreLegacy, 0},
    {71, "bm_c_uniform",   "O(rows*degree) random residual row all-but-one", kPoolLegacy, kScoreLegacy, 0},
    {72, "bm_local2",      "O(rows + 64*degree) bounded row-action local lookahead", kPoolLegacy, kScoreLegacy, 0},

    {73, "ks_d2_random",   "O(rows*degree) random degree-2 endpoint", kPoolLegacy, kScoreLegacy, 0},
    {74, "ks_smallcc",     "O(rows*degree + cols) smallest degree-2 component", kPoolLegacy, kScoreLegacy, 0},
    {75, "ks_boundary_min","O(rows*degree + cols + cc*rowrefs) largest degree-2 component/min boundary", kPoolLegacy, kScoreLegacy, 0},
    {76, "ks_boundary_max","O(rows*degree + cols + cc*rowrefs) largest degree-2 component/max boundary", kPoolLegacy, kScoreLegacy, 0},
    {77, "ks_fillmin",     "O(rows*degree + cols + cc*rowrefs) largest degree-2 component/min fill", kPoolLegacy, kScoreLegacy, 0},
    {78, "ks_lowref",      "O(rows*degree + cols) largest degree-2 component/low refs", kPoolLegacy, kScoreLegacy, 0},
    {79, "ks_random",      "O(rows*degree + cols) largest degree-2 component/random tie", kPoolLegacy, kScoreLegacy, 0},
    {80, "ks_all_d2_fill", "O(rows*degree + candidates*rowrefs) all degree-2 endpoints/min fill", kPoolLegacy, kScoreLegacy, 0},

    {81, "minrow_rowdeg",  "O(rows*degree + candidates) all min rows/row degree", kPoolLegacy, kScoreLegacy, 0},
    {82, "minrow_ovmin",   "O(rows*degree + candidates) all min rows/min overlap", kPoolLegacy, kScoreLegacy, 0},
    {83, "minrow_ovmax",   "O(rows*degree + candidates) all min rows/max overlap", kPoolLegacy, kScoreLegacy, 0},
    {84, "minrow_comp",    "O(rows*degree + candidates*rowrefs) all min rows/component impact", kPoolLegacy, kScoreLegacy, 0},
    {85, "minrow_fill",    "O(rows*degree + candidates*rowrefs) all min rows/estimated fill", kPoolLegacy, kScoreLegacy, 0},
    {86, "minrow_fillratio","O(rows*degree + candidates*rowrefs) all min rows/d2-fill ratio", kPoolLegacy, kScoreLegacy, 0},
    {87, "minrow_liveref", "O(rows*degree + candidates*rowrefs) all min rows/low live refs", kPoolLegacy, kScoreLegacy, 0},
    {88, "minrow_rand",    "O(rows*degree + candidates) all min rows/random tie", kPoolLegacy, kScoreLegacy, 0},

    {89, "beam_d2w4",      "O(top-K, tail<=8 depth*width*copy) beam depth 2 width 4", kPoolLegacy, kScoreLegacy, 0},
    {90, "beam_d2w8",      "O(top-K, tail<=8 depth*width*copy) beam depth 2 width 8", kPoolLegacy, kScoreLegacy, 0},
    {91, "beam_d3w4",      "O(top-K, tail<=8 depth*width*copy) beam depth 3 width 4", kPoolLegacy, kScoreLegacy, 0},
    {92, "beam_d3w8",      "O(top-K, tail<=8 depth*width*copy) beam depth 3 width 8", kPoolLegacy, kScoreLegacy, 0},

    {93, "hyb_lowref5",    "O(default then rqcc_lowref below 5pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {94, "hyb_lowref10",   "O(default then rqcc_lowref below 10pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {95, "hyb_ratio5",     "O(default then rqcc_ratio below 5pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {96, "hyb_ratio10",    "O(default then rqcc_ratio below 10pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {97, "hyb_beam5",      "O(default then beam depth2/width4 below 5pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {98, "hyb_beam10",     "O(default then beam depth2/width4 below 10pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {99, "hyb_rqbeam",     "O(rqcc_lowref then beam near tail)", kPoolLegacy, kScoreLegacy, 0},

    {100, "ms_rtie2",      "O(2*random_tie schedule) best of two random-tie runs", kPoolLegacy, kScoreLegacy, 0},
    {101, "ms_rcol2",      "O(2*random_col schedule) best of two random-column runs", kPoolLegacy, kScoreLegacy, 0},
    {102, "ms_rtie4",      "O(4*random_tie schedule) best of four random-tie runs", kPoolLegacy, kScoreLegacy, 0},
    {103, "ms_rcol4",      "O(4*random_col schedule) best of four random-column runs", kPoolLegacy, kScoreLegacy, 0},
    {104, "ms_def_rtie",   "O(default + random_tie multistart)", kPoolLegacy, kScoreLegacy, 0},
    {105, "ms_rqcc_rand",  "O(rqcc randomized tie multistart)", kPoolLegacy, kScoreLegacy, 0},

    {106, "rah_le2",       "O(rows*degree) row-action when live degree <= 2", kPoolLegacy, kScoreLegacy, 0},
    {107, "rah_le3",       "O(rows*degree) row-action when live degree <= 3", kPoolLegacy, kScoreLegacy, 0},
    {108, "rah_le4",       "O(rows*degree) row-action when live degree <= 4", kPoolLegacy, kScoreLegacy, 0},
    {109, "rah_bm3",       "O(rows*degree) BM row-action below degree 3 else singleton", kPoolLegacy, kScoreLegacy, 0},
    {110, "rah_local3",    "O(rows + 64*degree) bounded local row-action below degree 3 else rqcc", kPoolLegacy, kScoreLegacy, 0},
    {111, "rah_min4",      "O(rows*degree) min-row action below degree 4 else rqcc", kPoolLegacy, kScoreLegacy, 0},

    {112, "lazy_liveref",  "O(cols + touched_rows*degree) incremental live-ref cache", kPoolLegacy, kScoreLegacy, 0},
    {113, "lazy_fill",     "O(cols + touched_rows*degree) incremental fill cache", kPoolLegacy, kScoreLegacy, 0},
    {114, "lazy_d2",       "O(cols + touched_rows*degree) incremental degree-2 cache", kPoolLegacy, kScoreLegacy, 0},

    {115, "ks_bmax_top16", "O(rows*degree + cols + min(16,cc)*rowrefs) degree-2 CC/default top16 then boundary-max", kPoolLegacy, kScoreLegacy, 0},
    {116, "ks_bmax_top64", "O(rows*degree + cols + min(64,cc)*rowrefs) degree-2 CC/default top64 then boundary-max", kPoolLegacy, kScoreLegacy, 0},
    {117, "hyb_bmax5",     "O(rqd2_default then boundary-max below 5pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {118, "hyb_bmax10",    "O(rqd2_default then boundary-max below 10pct todo)", kPoolLegacy, kScoreLegacy, 0},

    {119, "av_top16",      "O(cols + 16*local-avalanche) top16/default + local cascade estimate", kPoolLegacy, kScoreLegacy, 0},
    {120, "av_top64",      "O(cols + 64*local-avalanche) top64/default + local cascade estimate", kPoolLegacy, kScoreLegacy, 0},
    {121, "av_rqcc16",     "O(rows*degree + cols + min(16,cc)*local-avalanche) degree-2 CC + local cascade estimate", kPoolLegacy, kScoreLegacy, 0},
    {122, "av_rqcc64",     "O(rows*degree + cols + min(64,cc)*local-avalanche) degree-2 CC + local cascade estimate", kPoolLegacy, kScoreLegacy, 0},
    {123, "sim1_top16",    "O(cols + 16*copy) top16/default + exact one-step score", kPoolLegacy, kScoreLegacy, 0},
    {124, "sim1_rqcc16",   "O(rows*degree + cols + min(16,cc)*copy) degree-2 CC + exact one-step score", kPoolLegacy, kScoreLegacy, 0},
    {125, "hyb_rqav5",     "O(rqd2_default then rqcc16 avalanche-est below 5pct todo)", kPoolLegacy, kScoreLegacy, 0},
    {126, "hyb_rqav10",    "O(rqd2_default then rqcc16 avalanche-est below 10pct todo)", kPoolLegacy, kScoreLegacy, 0},

    {127, "av_lazy16",     "O(cols + dirty local-avalanche) av_top16 with dependency-invalidated cache", kPoolLegacy, kScoreLegacy, 0},
    {128, "av_lazy64",     "O(cols + dirty local-avalanche) av_top64 with dependency-invalidated cache", kPoolLegacy, kScoreLegacy, 0},
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

enum StructureKind
{
    kStructureWirehairRandom = 0,
    kStructureUniform,
    kStructureGaussian,
    kStructureLT,
    kStructureHarmonic,
    kStructureOnePlusUniform,
    kStructureOnePlusLT,
    kStructureTwoModeUniform,
    kStructureBalancedUniform,
    kStructureRobustSoliton,
    kStructureLTAlpha,
    kStructureLTExtra34,
    kStructureLTFold,
    kStructureLTFoldScale,
    kStructureLTHighRows,
    kStructureOverheadLT,
    kStructureRaptorQLdpcStruct,
    kStructureRaptorQLdpcRandom,
    kStructureBinaryP50,
    kStructureDeckColumns,
};

struct MatrixStructure
{
    const char* Name;
    const char* Description;
    StructureKind Kind;
    unsigned MinDegree;
    unsigned MaxDegree;
    double A;
    double B;
    double C;
    unsigned AltMinDegree;
    unsigned AltMaxDegree;
    unsigned BalanceMode;
    double BalanceValue;
    // Underlying row-degree law used by wrapper kinds (kStructureDeckColumns);
    // defaults to Kind itself for every ordinary structure.
    StructureKind DegreeKind;

    MatrixStructure(
        const char* name,
        const char* description,
        StructureKind kind,
        unsigned min_degree,
        unsigned max_degree,
        double a = 0.0,
        double b = 0.0,
        double c = 0.0,
        unsigned alt_min_degree = 0,
        unsigned alt_max_degree = 0,
        unsigned balance_mode = 0,
        double balance_value = 0.0,
        int degree_kind = -1)
        : Name(name),
          Description(description),
          Kind(kind),
          MinDegree(min_degree),
          MaxDegree(max_degree),
          A(a),
          B(b),
          C(c),
          AltMinDegree(alt_min_degree),
          AltMaxDegree(alt_max_degree),
          BalanceMode(balance_mode),
          BalanceValue(balance_value),
          DegreeKind(degree_kind >= 0 ? (StructureKind)degree_kind : kind)
    {
    }
};

static const MatrixStructure kStructures[] = {
    {"wirehair_rand", "Wirehair row-weight distribution with random columns", kStructureWirehairRandom, 1, 65535},
    {"binary_p50", "fully random dense binary matrix rows with P(column)=0.5", kStructureBinaryP50, 1, 65535},
    {"uniform1_8", "uniform row degree in [1,8]", kStructureUniform, 1, 8, 0.0, 0.0},
    {"uniform2_8", "uniform row degree in [2,8]", kStructureUniform, 2, 8, 0.0, 0.0},
    {"uniform3_8", "uniform row degree in [3,8]", kStructureUniform, 3, 8, 0.0, 0.0},
    {"uniform1_16", "uniform row degree in [1,16]", kStructureUniform, 1, 16, 0.0, 0.0},
    {"uniform2_16", "uniform row degree in [2,16]", kStructureUniform, 2, 16, 0.0, 0.0},
    {"uniform3_16", "uniform row degree in [3,16]", kStructureUniform, 3, 16, 0.0, 0.0},
    {"gauss4_min1", "rounded gaussian degree around 4, min 1, cap 16", kStructureGaussian, 1, 16, 4.0, 1.5},
    {"gauss4_min2", "rounded gaussian degree around 4, min 2, cap 16", kStructureGaussian, 2, 16, 4.0, 1.5},
    {"gauss6_min2", "rounded gaussian degree around 6, min 2, cap 24", kStructureGaussian, 2, 24, 6.0, 2.0},
    {"lt_trunc64", "truncated ideal-soliton-like degree distribution, cap 64", kStructureLT, 1, 64, 0.0, 0.0},
    {"lt_no1_64", "truncated ideal-soliton-like distribution without degree 1", kStructureLT, 2, 64, 0.0, 0.0},
    {"lt_robust64", "LT-like distribution with extra degree-1/2 mass, cap 64", kStructureLT, 1, 64, 0.03, 0.08},
    {"lt_m1_c16", "LT degree distribution, min 1, cap 16", kStructureLT, 1, 16},
    {"lt_m1_c32", "LT degree distribution, min 1, cap 32", kStructureLT, 1, 32},
    {"lt_m1_c64", "LT degree distribution, min 1, cap 64", kStructureLT, 1, 64},
    {"lt_m1_c128", "LT degree distribution, min 1, cap 128", kStructureLT, 1, 128},
    {"lt_m1_c256", "LT degree distribution, min 1, cap 256", kStructureLT, 1, 256},
    {"lt_m1_c320", "LT degree distribution, min 1, cap 320", kStructureLT, 1, 320},
    {"lt_m1_c3200", "LT degree distribution, min 1, cap 3200", kStructureLT, 1, 3200},
    {"lt_m2_c16", "LT degree distribution, min 2, cap 16", kStructureLT, 2, 16},
    {"lt_m2_c32", "LT degree distribution, min 2, cap 32", kStructureLT, 2, 32},
    {"lt_m2_c64", "LT degree distribution, min 2, cap 64", kStructureLT, 2, 64},
    {"lt_m2_c96", "LT degree distribution, min 2, cap 96", kStructureLT, 2, 96},
    {"lt_m2_c128", "LT degree distribution, min 2, cap 128", kStructureLT, 2, 128},
    {"lt_m2_c192", "LT degree distribution, min 2, cap 192", kStructureLT, 2, 192},
    {"lt_m2_c256", "LT degree distribution, min 2, cap 256", kStructureLT, 2, 256},
    {"lt_m2_c320", "LT degree distribution, min 2, cap 320", kStructureLT, 2, 320},
    {"lt_m2_c512", "LT degree distribution, min 2, cap 512", kStructureLT, 2, 512},
    {"lt_m2_c640", "LT degree distribution, min 2, cap 640", kStructureLT, 2, 640},
    {"lt_m2_c960", "LT degree distribution, min 2, cap 960", kStructureLT, 2, 960},
    {"lt_m2_c1024", "LT degree distribution, min 2, cap 1024", kStructureLT, 2, 1024},
    {"lt_m2_c1280", "LT degree distribution, min 2, cap 1280", kStructureLT, 2, 1280},
    {"lt_m2_c1920", "LT degree distribution, min 2, cap 1920", kStructureLT, 2, 1920},
    {"lt_m2_c2560", "LT degree distribution, min 2, cap 2560", kStructureLT, 2, 2560},
    {"lt_m2_c3200", "LT degree distribution, min 2, cap 3200", kStructureLT, 2, 3200},
    {"lt_m3_c16", "LT degree distribution, min 3, cap 16", kStructureLT, 3, 16},
    {"lt_m3_c32", "LT degree distribution, min 3, cap 32", kStructureLT, 3, 32},
    {"lt_m3_c64", "LT degree distribution, min 3, cap 64", kStructureLT, 3, 64},
    {"lt_m3_c128", "LT degree distribution, min 3, cap 128", kStructureLT, 3, 128},
    {"lt_m3_c256", "LT degree distribution, min 3, cap 256", kStructureLT, 3, 256},
    {"lt_m3_c320", "LT degree distribution, min 3, cap 320", kStructureLT, 3, 320},
    {"lt_m3_c3200", "LT degree distribution, min 3, cap 3200", kStructureLT, 3, 3200},
    {"d1mix_lt_p0", "0 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.00},
    {"d1mix_lt_p005", "0.5 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.005},
    {"d1mix_lt_p01", "1 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.01},
    {"d1mix_lt_p2", "2 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.02},
    {"d1mix_lt_p5", "5 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.05},
    {"d1mix_lt_p10", "10 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.10},
    {"d1mix_lt_p20", "20 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.20},
    {"d1mix_lt_p30", "30 percent degree-1 rows, otherwise LT min 2 cap 64", kStructureOnePlusLT, 2, 64, 0.30},
    {"d1mix_u2_8_p0", "0 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.00},
    {"d1mix_u2_8_p2", "2 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.02},
    {"d1mix_u2_8_p5", "5 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.05},
    {"d1mix_u2_8_p10", "10 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.10},
    {"d1mix_u2_8_p20", "20 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.20},
    {"d1mix_u2_8_p30", "30 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.30},
    {"d1mix_u2_16_p0", "0 percent degree-1 rows, otherwise uniform [2,16]", kStructureOnePlusUniform, 2, 16, 0.00},
    {"d1mix_u2_16_p2", "2 percent degree-1 rows, otherwise uniform [2,16]", kStructureOnePlusUniform, 2, 16, 0.02},
    {"d1mix_u2_16_p5", "5 percent degree-1 rows, otherwise uniform [2,16]", kStructureOnePlusUniform, 2, 16, 0.05},
    {"d1mix_u2_16_p10", "10 percent degree-1 rows, otherwise uniform [2,16]", kStructureOnePlusUniform, 2, 16, 0.10},
    {"d1mix_u2_16_p20", "20 percent degree-1 rows, otherwise uniform [2,16]", kStructureOnePlusUniform, 2, 16, 0.20},
    {"d1mix_u2_16_p30", "30 percent degree-1 rows, otherwise uniform [2,16]", kStructureOnePlusUniform, 2, 16, 0.30},
    {"d1mix_u2_64_p0", "0 percent degree-1 rows, otherwise uniform [2,64]", kStructureOnePlusUniform, 2, 64, 0.00},
    {"d1mix_u2_64_p2", "2 percent degree-1 rows, otherwise uniform [2,64]", kStructureOnePlusUniform, 2, 64, 0.02},
    {"d1mix_u2_64_p5", "5 percent degree-1 rows, otherwise uniform [2,64]", kStructureOnePlusUniform, 2, 64, 0.05},
    {"d1mix_u2_64_p10", "10 percent degree-1 rows, otherwise uniform [2,64]", kStructureOnePlusUniform, 2, 64, 0.10},
    {"d1mix_u2_64_p20", "20 percent degree-1 rows, otherwise uniform [2,64]", kStructureOnePlusUniform, 2, 64, 0.20},
    {"d1mix_u2_64_p30", "30 percent degree-1 rows, otherwise uniform [2,64]", kStructureOnePlusUniform, 2, 64, 0.30},
    {"robust_d1_001_d2_003", "LT cap64 plus degree-1 mass 0.01 and degree-2 mass 0.03", kStructureLT, 1, 64, 0.01, 0.03},
    {"robust_d1_001_d2_008", "LT cap64 plus degree-1 mass 0.01 and degree-2 mass 0.08", kStructureLT, 1, 64, 0.01, 0.08},
    {"robust_d1_001_d2_012", "LT cap64 plus degree-1 mass 0.01 and degree-2 mass 0.12", kStructureLT, 1, 64, 0.01, 0.12},
    {"robust_d1_001_d2_020", "LT cap64 plus degree-1 mass 0.01 and degree-2 mass 0.20", kStructureLT, 1, 64, 0.01, 0.20},
    {"robust_d1_003_d2_003", "LT cap64 plus degree-1 mass 0.03 and degree-2 mass 0.03", kStructureLT, 1, 64, 0.03, 0.03},
    {"robust_d1_003_d2_008", "LT cap64 plus degree-1 mass 0.03 and degree-2 mass 0.08", kStructureLT, 1, 64, 0.03, 0.08},
    {"robust_d1_003_d2_012", "LT cap64 plus degree-1 mass 0.03 and degree-2 mass 0.12", kStructureLT, 1, 64, 0.03, 0.12},
    {"robust_d1_003_d2_020", "LT cap64 plus degree-1 mass 0.03 and degree-2 mass 0.20", kStructureLT, 1, 64, 0.03, 0.20},
    {"robust_d1_005_d2_003", "LT cap64 plus degree-1 mass 0.05 and degree-2 mass 0.03", kStructureLT, 1, 64, 0.05, 0.03},
    {"robust_d1_005_d2_008", "LT cap64 plus degree-1 mass 0.05 and degree-2 mass 0.08", kStructureLT, 1, 64, 0.05, 0.08},
    {"robust_d1_005_d2_012", "LT cap64 plus degree-1 mass 0.05 and degree-2 mass 0.12", kStructureLT, 1, 64, 0.05, 0.12},
    {"robust_d1_005_d2_020", "LT cap64 plus degree-1 mass 0.05 and degree-2 mass 0.20", kStructureLT, 1, 64, 0.05, 0.20},
    {"robust_d1_010_d2_003", "LT cap64 plus degree-1 mass 0.10 and degree-2 mass 0.03", kStructureLT, 1, 64, 0.10, 0.03},
    {"robust_d1_010_d2_008", "LT cap64 plus degree-1 mass 0.10 and degree-2 mass 0.08", kStructureLT, 1, 64, 0.10, 0.08},
    {"robust_d1_010_d2_012", "LT cap64 plus degree-1 mass 0.10 and degree-2 mass 0.12", kStructureLT, 1, 64, 0.10, 0.12},
    {"robust_d1_010_d2_020", "LT cap64 plus degree-1 mass 0.10 and degree-2 mass 0.20", kStructureLT, 1, 64, 0.10, 0.20},
    {"rs_c001_d50_c128", "true robust soliton c=0.01 delta=0.50, cap 128, renormalized", kStructureRobustSoliton, 1, 128, 0.01, 0.50},
    {"rs_c003_d10_c128", "true robust soliton c=0.03 delta=0.10, cap 128, renormalized", kStructureRobustSoliton, 1, 128, 0.03, 0.10},
    {"rs_c005_d05_c128", "true robust soliton c=0.05 delta=0.05, cap 128, renormalized", kStructureRobustSoliton, 1, 128, 0.05, 0.05},
    {"rs_c010_d01_c128", "true robust soliton c=0.10 delta=0.01, cap 128, renormalized", kStructureRobustSoliton, 1, 128, 0.10, 0.01},
    {"sol_alpha075_c128", "ideal soliton weights raised to alpha=0.75, cap 128", kStructureLTAlpha, 1, 128, 0.75},
    {"sol_alpha125_c128", "ideal soliton weights raised to alpha=1.25, cap 128", kStructureLTAlpha, 1, 128, 1.25},
    {"sol_alpha150_c128", "ideal soliton weights raised to alpha=1.50, cap 128", kStructureLTAlpha, 1, 128, 1.50},
    {"sol_m2_alpha075_c320", "LT min 2 cap 320 with weights raised to alpha=0.75", kStructureLTAlpha, 2, 320, 0.75},
    {"sol_m2_alpha125_c320", "LT min 2 cap 320 with weights raised to alpha=1.25", kStructureLTAlpha, 2, 320, 1.25},
    {"sol_m2_alpha150_c320", "LT min 2 cap 320 with weights raised to alpha=1.50", kStructureLTAlpha, 2, 320, 1.50},
    {"sol_m2_alpha075_c3200", "LT min 2 cap 3200 with weights raised to alpha=0.75", kStructureLTAlpha, 2, 3200, 0.75},
    {"sol_m2_alpha125_c3200", "LT min 2 cap 3200 with weights raised to alpha=1.25", kStructureLTAlpha, 2, 3200, 1.25},
    {"sol_m2_alpha150_c3200", "LT min 2 cap 3200 with weights raised to alpha=1.50", kStructureLTAlpha, 2, 3200, 1.50},
    {"lt_m1_c64_fold", "LT min 1 cap 64 with probability mass above cap folded into cap", kStructureLTFold, 1, 64},
    {"lt_m2_c64_fold", "LT min 2 cap 64 with probability mass above cap folded into cap", kStructureLTFold, 2, 64},
    {"lt_m2_c96_fold", "LT min 2 cap 96 with probability mass above cap folded into cap", kStructureLTFold, 2, 96},
    {"lt_m2_c128_fold", "LT min 2 cap 128 with probability mass above cap folded into cap", kStructureLTFold, 2, 128},
    {"lt_m2_c192_fold", "LT min 2 cap 192 with probability mass above cap folded into cap", kStructureLTFold, 2, 192},
    {"lt_m2_c256_fold", "LT min 2 cap 256 with probability mass above cap folded into cap", kStructureLTFold, 2, 256},
    {"lt_m2_c320_fold", "LT min 2 cap 320 with probability mass above cap folded into cap", kStructureLTFold, 2, 320},
    {"lt_m2_c512_fold", "LT min 2 cap 512 with probability mass above cap folded into cap", kStructureLTFold, 2, 512},
    {"lt_m2_c640_fold", "LT min 2 cap 640 with probability mass above cap folded into cap", kStructureLTFold, 2, 640},
    {"lt_m2_c960_fold", "LT min 2 cap 960 with probability mass above cap folded into cap", kStructureLTFold, 2, 960},
    {"lt_m2_c1024_fold", "LT min 2 cap 1024 with probability mass above cap folded into cap", kStructureLTFold, 2, 1024},
    {"lt_m2_c1280_fold", "LT min 2 cap 1280 with probability mass above cap folded into cap", kStructureLTFold, 2, 1280},
    {"lt_m2_c1920_fold", "LT min 2 cap 1920 with probability mass above cap folded into cap", kStructureLTFold, 2, 1920},
    {"lt_m2_c2560_fold", "LT min 2 cap 2560 with probability mass above cap folded into cap", kStructureLTFold, 2, 2560},
    {"lt_m2_c3200_fold", "LT min 2 cap 3200 with probability mass above cap folded into cap", kStructureLTFold, 2, 3200},
    {"lt_m2_c256_fold0", "LT min 2 cap 256 with 0 percent folded tail mass", kStructureLTFoldScale, 2, 256, 0.0},
    {"lt_m2_c256_fold25", "LT min 2 cap 256 with 25 percent folded tail mass", kStructureLTFoldScale, 2, 256, 0.25},
    {"lt_m2_c256_fold50", "LT min 2 cap 256 with 50 percent folded tail mass", kStructureLTFoldScale, 2, 256, 0.50},
    {"lt_m2_c256_fold200", "LT min 2 cap 256 with 200 percent folded tail mass", kStructureLTFoldScale, 2, 256, 2.00},
    {"lt_m2_c1024_fold0", "LT min 2 cap 1024 with 0 percent folded tail mass", kStructureLTFoldScale, 2, 1024, 0.0},
    {"lt_m2_c1024_fold25", "LT min 2 cap 1024 with 25 percent folded tail mass", kStructureLTFoldScale, 2, 1024, 0.25},
    {"lt_m2_c1024_fold50", "LT min 2 cap 1024 with 50 percent folded tail mass", kStructureLTFoldScale, 2, 1024, 0.50},
    {"lt_m2_c1024_fold200", "LT min 2 cap 1024 with 200 percent folded tail mass", kStructureLTFoldScale, 2, 1024, 2.00},
    {"lt_m2_c384_fold10", "LT min 2 cap 384 with 10 percent folded tail mass", kStructureLTFoldScale, 2, 384, 0.10},
    {"lt_m2_c384_fold25", "LT min 2 cap 384 with 25 percent folded tail mass", kStructureLTFoldScale, 2, 384, 0.25},
    {"lt_m2_c384_fold40", "LT min 2 cap 384 with 40 percent folded tail mass", kStructureLTFoldScale, 2, 384, 0.40},
    {"lt_m2_c384_fold60", "LT min 2 cap 384 with 60 percent folded tail mass", kStructureLTFoldScale, 2, 384, 0.60},
    {"lt_m2_c512_fold10", "LT min 2 cap 512 with 10 percent folded tail mass", kStructureLTFoldScale, 2, 512, 0.10},
    {"lt_m2_c512_fold25", "LT min 2 cap 512 with 25 percent folded tail mass", kStructureLTFoldScale, 2, 512, 0.25},
    {"lt_m2_c512_fold40", "LT min 2 cap 512 with 40 percent folded tail mass", kStructureLTFoldScale, 2, 512, 0.40},
    {"lt_m2_c512_fold60", "LT min 2 cap 512 with 60 percent folded tail mass", kStructureLTFoldScale, 2, 512, 0.60},
    {"lt_m2_c768_fold10", "LT min 2 cap 768 with 10 percent folded tail mass", kStructureLTFoldScale, 2, 768, 0.10},
    {"lt_m2_c768_fold25", "LT min 2 cap 768 with 25 percent folded tail mass", kStructureLTFoldScale, 2, 768, 0.25},
    {"lt_m2_c768_fold40", "LT min 2 cap 768 with 40 percent folded tail mass", kStructureLTFoldScale, 2, 768, 0.40},
    {"lt_m2_c768_fold60", "LT min 2 cap 768 with 60 percent folded tail mass", kStructureLTFoldScale, 2, 768, 0.60},
    {"lt_m2_c1024_fold10", "LT min 2 cap 1024 with 10 percent folded tail mass", kStructureLTFoldScale, 2, 1024, 0.10},
    {"lt_m2_c1024_fold40", "LT min 2 cap 1024 with 40 percent folded tail mass", kStructureLTFoldScale, 2, 1024, 0.40},
    {"lt_m2_c1024_fold60", "LT min 2 cap 1024 with 60 percent folded tail mass", kStructureLTFoldScale, 2, 1024, 0.60},
    {"lt_m2_c1536_fold10", "LT min 2 cap 1536 with 10 percent folded tail mass", kStructureLTFoldScale, 2, 1536, 0.10},
    {"lt_m2_c1536_fold25", "LT min 2 cap 1536 with 25 percent folded tail mass", kStructureLTFoldScale, 2, 1536, 0.25},
    {"lt_m2_c1536_fold40", "LT min 2 cap 1536 with 40 percent folded tail mass", kStructureLTFoldScale, 2, 1536, 0.40},
    {"lt_m2_c1536_fold60", "LT min 2 cap 1536 with 60 percent folded tail mass", kStructureLTFoldScale, 2, 1536, 0.60},
    {"lt_m2_c2048_fold10", "LT min 2 cap 2048 with 10 percent folded tail mass", kStructureLTFoldScale, 2, 2048, 0.10},
    {"lt_m2_c2048_fold25", "LT min 2 cap 2048 with 25 percent folded tail mass", kStructureLTFoldScale, 2, 2048, 0.25},
    {"lt_m2_c2048_fold40", "LT min 2 cap 2048 with 40 percent folded tail mass", kStructureLTFoldScale, 2, 2048, 0.40},
    {"lt_m2_c2048_fold60", "LT min 2 cap 2048 with 60 percent folded tail mass", kStructureLTFoldScale, 2, 2048, 0.60},
    {"lt_d3_003_d4_003", "LT cap64 plus independent degree-3 mass 0.03 and degree-4 mass 0.03", kStructureLTExtra34, 1, 64, 0.03, 0.03, 0.0},
    {"lt_d3_008_d4_003", "LT cap64 plus independent degree-3 mass 0.08 and degree-4 mass 0.03", kStructureLTExtra34, 1, 64, 0.08, 0.03, 0.0},
    {"lt_d3_003_d4_008", "LT cap64 plus independent degree-3 mass 0.03 and degree-4 mass 0.08", kStructureLTExtra34, 1, 64, 0.03, 0.08, 0.0},
    {"lt_d3_008_d4_008", "LT cap64 plus independent degree-3 mass 0.08 and degree-4 mass 0.08", kStructureLTExtra34, 1, 64, 0.08, 0.08, 0.0},
    {"lt_d2_008_d3_003_d4_003", "LT cap64 plus degree-2 mass 0.08, degree-3 mass 0.03, degree-4 mass 0.03", kStructureLTExtra34, 1, 64, 0.03, 0.03, 0.08},
    {"lt_c320_d2_003", "LT cap320 plus independent degree-2 mass 0.03", kStructureLTExtra34, 1, 320, 0.0, 0.0, 0.03},
    {"lt_c320_d2_008", "LT cap320 plus independent degree-2 mass 0.08", kStructureLTExtra34, 1, 320, 0.0, 0.0, 0.08},
    {"lt_c320_d3_003", "LT cap320 plus independent degree-3 mass 0.03", kStructureLTExtra34, 1, 320, 0.03, 0.0, 0.0},
    {"lt_c320_d4_003", "LT cap320 plus independent degree-4 mass 0.03", kStructureLTExtra34, 1, 320, 0.0, 0.03, 0.0},
    {"lt_c320_d2_003_d3_003_d4_003", "LT cap320 plus degree-2/3/4 mass 0.03 each", kStructureLTExtra34, 1, 320, 0.03, 0.03, 0.03},
    {"lt_c320_d2_008_d3_003_d4_003", "LT cap320 plus degree-2 mass 0.08, degree-3/4 mass 0.03", kStructureLTExtra34, 1, 320, 0.03, 0.03, 0.08},
    {"lt_m2_c320_d2_003", "LT min 2 cap320 plus independent degree-2 mass 0.03", kStructureLTExtra34, 2, 320, 0.0, 0.0, 0.03},
    {"lt_m2_c320_d2_008", "LT min 2 cap320 plus independent degree-2 mass 0.08", kStructureLTExtra34, 2, 320, 0.0, 0.0, 0.08},
    {"lt_m2_c320_d3_003", "LT min 2 cap320 plus independent degree-3 mass 0.03", kStructureLTExtra34, 2, 320, 0.03, 0.0, 0.0},
    {"lt_m2_c320_d4_003", "LT min 2 cap320 plus independent degree-4 mass 0.03", kStructureLTExtra34, 2, 320, 0.0, 0.03, 0.0},
    {"lt_m2_c320_d2_003_d3_003_d4_003", "LT min 2 cap320 plus degree-2/3/4 mass 0.03 each", kStructureLTExtra34, 2, 320, 0.03, 0.03, 0.03},
    {"lt_m2_c320_d2_008_d3_003_d4_003", "LT min 2 cap320 plus degree-2 mass 0.08, degree-3/4 mass 0.03", kStructureLTExtra34, 2, 320, 0.03, 0.03, 0.08},
    {"lt_c3200_d2_003", "LT cap3200 plus independent degree-2 mass 0.03", kStructureLTExtra34, 1, 3200, 0.0, 0.0, 0.03},
    {"lt_c3200_d2_008", "LT cap3200 plus independent degree-2 mass 0.08", kStructureLTExtra34, 1, 3200, 0.0, 0.0, 0.08},
    {"lt_c3200_d3_003", "LT cap3200 plus independent degree-3 mass 0.03", kStructureLTExtra34, 1, 3200, 0.03, 0.0, 0.0},
    {"lt_c3200_d4_003", "LT cap3200 plus independent degree-4 mass 0.03", kStructureLTExtra34, 1, 3200, 0.0, 0.03, 0.0},
    {"lt_c3200_d2_003_d3_003_d4_003", "LT cap3200 plus degree-2/3/4 mass 0.03 each", kStructureLTExtra34, 1, 3200, 0.03, 0.03, 0.03},
    {"lt_c3200_d2_008_d3_003_d4_003", "LT cap3200 plus degree-2 mass 0.08, degree-3/4 mass 0.03", kStructureLTExtra34, 1, 3200, 0.03, 0.03, 0.08},
    {"lt_m2_c3200_d2_003", "LT min 2 cap3200 plus independent degree-2 mass 0.03", kStructureLTExtra34, 2, 3200, 0.0, 0.0, 0.03},
    {"lt_m2_c3200_d2_008", "LT min 2 cap3200 plus independent degree-2 mass 0.08", kStructureLTExtra34, 2, 3200, 0.0, 0.0, 0.08},
    {"lt_m2_c3200_d3_003", "LT min 2 cap3200 plus independent degree-3 mass 0.03", kStructureLTExtra34, 2, 3200, 0.03, 0.0, 0.0},
    {"lt_m2_c3200_d4_003", "LT min 2 cap3200 plus independent degree-4 mass 0.03", kStructureLTExtra34, 2, 3200, 0.0, 0.03, 0.0},
    {"lt_m2_c3200_d2_003_d3_003_d4_003", "LT min 2 cap3200 plus degree-2/3/4 mass 0.03 each", kStructureLTExtra34, 2, 3200, 0.03, 0.03, 0.03},
    {"lt_m2_c3200_d2_008_d3_003_d4_003", "LT min 2 cap3200 plus degree-2 mass 0.08, degree-3/4 mass 0.03", kStructureLTExtra34, 2, 3200, 0.03, 0.03, 0.08},
    {"lt_m2_c512_d2_003", "LT min 2 cap512 plus independent degree-2 mass 0.03", kStructureLTExtra34, 2, 512, 0.0, 0.0, 0.03},
    {"lt_m2_c512_d3_003", "LT min 2 cap512 plus independent degree-3 mass 0.03", kStructureLTExtra34, 2, 512, 0.03, 0.0, 0.0},
    {"lt_m2_c1024_d2_003", "LT min 2 cap1024 plus independent degree-2 mass 0.03", kStructureLTExtra34, 2, 1024, 0.0, 0.0, 0.03},
    {"lt_m2_c1024_d3_003", "LT min 2 cap1024 plus independent degree-3 mass 0.03", kStructureLTExtra34, 2, 1024, 0.03, 0.0, 0.0},
    {"lt_m2_c2048_d2_003", "LT min 2 cap2048 plus independent degree-2 mass 0.03", kStructureLTExtra34, 2, 2048, 0.0, 0.0, 0.03},
    {"lt_m2_c2048_d3_003", "LT min 2 cap2048 plus independent degree-3 mass 0.03", kStructureLTExtra34, 2, 2048, 0.03, 0.0, 0.0},
    {"lt_m2_c320_hi1_h160", "LT min 2 cap 320 with 1 explicit degree-N/2 row", kStructureLTHighRows, 2, 320, 1.0, 0.5},
    {"lt_m2_c320_hi2_h160", "LT min 2 cap 320 with 2 explicit degree-N/2 rows", kStructureLTHighRows, 2, 320, 2.0, 0.5},
    {"lt_m2_c320_hi4_h160", "LT min 2 cap 320 with 4 explicit degree-N/2 rows", kStructureLTHighRows, 2, 320, 4.0, 0.5},
    {"lt_m2_c320_hi8_h160", "LT min 2 cap 320 with 8 explicit degree-N/2 rows", kStructureLTHighRows, 2, 320, 8.0, 0.5},
    {"lt_m2_c320_hi1_h320", "LT min 2 cap 320 with 1 explicit degree-N row", kStructureLTHighRows, 2, 320, 1.0, 1.0},
    {"lt_m2_c320_hi2_h320", "LT min 2 cap 320 with 2 explicit degree-N rows", kStructureLTHighRows, 2, 320, 2.0, 1.0},
    {"lt_m2_c320_hi4_h320", "LT min 2 cap 320 with 4 explicit degree-N rows", kStructureLTHighRows, 2, 320, 4.0, 1.0},
    {"lt_m2_c320_hi8_h320", "LT min 2 cap 320 with 8 explicit degree-N rows", kStructureLTHighRows, 2, 320, 8.0, 1.0},
    {"lt_m2_c3200_hi1_h1600", "LT min 2 cap 3200 with 1 explicit degree-N/2 row", kStructureLTHighRows, 2, 3200, 1.0, 0.5},
    {"lt_m2_c3200_hi2_h1600", "LT min 2 cap 3200 with 2 explicit degree-N/2 rows", kStructureLTHighRows, 2, 3200, 2.0, 0.5},
    {"lt_m2_c3200_hi4_h1600", "LT min 2 cap 3200 with 4 explicit degree-N/2 rows", kStructureLTHighRows, 2, 3200, 4.0, 0.5},
    {"lt_m2_c3200_hi8_h1600", "LT min 2 cap 3200 with 8 explicit degree-N/2 rows", kStructureLTHighRows, 2, 3200, 8.0, 0.5},
    {"lt_m2_c3200_hi1_h3200", "LT min 2 cap 3200 with 1 explicit degree-N row", kStructureLTHighRows, 2, 3200, 1.0, 1.0},
    {"lt_m2_c3200_hi2_h3200", "LT min 2 cap 3200 with 2 explicit degree-N rows", kStructureLTHighRows, 2, 3200, 2.0, 1.0},
    {"lt_m2_c3200_hi4_h3200", "LT min 2 cap 3200 with 4 explicit degree-N rows", kStructureLTHighRows, 2, 3200, 4.0, 1.0},
    {"lt_m2_c3200_hi8_h3200", "LT min 2 cap 3200 with 8 explicit degree-N rows", kStructureLTHighRows, 2, 3200, 8.0, 1.0},
    {"ohdep_lt_adapt128", "overhead-dependent LT cap128: low-degree mass adapts from +0 to +5 percent overhead", kStructureOverheadLT, 1, 128, 0.03, 0.10, 0.03},
    {"ohdep_lt_no1_128", "overhead-dependent LT cap128: suppresses degree-1 mass at zero/near-zero overhead", kStructureOverheadLT, 1, 128, 0.00, 0.12, 0.04},
    {"ohdep_lt_tail256", "overhead-dependent LT cap256: longer tail at zero overhead with stronger low-degree mass at overhead", kStructureOverheadLT, 1, 256, 0.02, 0.08, 0.02},
    {"raptorq_ldpc_struct", "RaptorQ LDPC check-row pattern over fixed N columns plus random repair rows", kStructureRaptorQLdpcStruct, 1, 65535},
    {"raptorq_ldpc_rand", "randomized control with the same RaptorQ LDPC check-row degrees plus random repair rows", kStructureRaptorQLdpcRandom, 1, 65535},
    {"two_low24_t5_h8_32", "95 percent uniform [2,4], 5 percent uniform [8,32]", kStructureTwoModeUniform, 2, 4, 0.05, 0.0, 0.0, 8, 32},
    {"two_low24_t5_h16_64", "95 percent uniform [2,4], 5 percent uniform [16,64]", kStructureTwoModeUniform, 2, 4, 0.05, 0.0, 0.0, 16, 64},
    {"two_low24_t5_h32_128", "95 percent uniform [2,4], 5 percent uniform [32,128]", kStructureTwoModeUniform, 2, 4, 0.05, 0.0, 0.0, 32, 128},
    {"two_low24_t10_h8_32", "90 percent uniform [2,4], 10 percent uniform [8,32]", kStructureTwoModeUniform, 2, 4, 0.10, 0.0, 0.0, 8, 32},
    {"two_low24_t10_h16_64", "90 percent uniform [2,4], 10 percent uniform [16,64]", kStructureTwoModeUniform, 2, 4, 0.10, 0.0, 0.0, 16, 64},
    {"two_low24_t10_h32_128", "90 percent uniform [2,4], 10 percent uniform [32,128]", kStructureTwoModeUniform, 2, 4, 0.10, 0.0, 0.0, 32, 128},
    {"two_low24_t20_h8_32", "80 percent uniform [2,4], 20 percent uniform [8,32]", kStructureTwoModeUniform, 2, 4, 0.20, 0.0, 0.0, 8, 32},
    {"two_low24_t20_h16_64", "80 percent uniform [2,4], 20 percent uniform [16,64]", kStructureTwoModeUniform, 2, 4, 0.20, 0.0, 0.0, 16, 64},
    {"two_low24_t20_h32_128", "80 percent uniform [2,4], 20 percent uniform [32,128]", kStructureTwoModeUniform, 2, 4, 0.20, 0.0, 0.0, 32, 128},
    {"two_low23_t10_h16_64", "90 percent uniform [2,3], 10 percent uniform [16,64]", kStructureTwoModeUniform, 2, 3, 0.10, 0.0, 0.0, 16, 64},
    {"bal_u2_16_load125", "balanced uniform [2,16], reject above 1.25x mean load", kStructureBalancedUniform, 2, 16, 0.0, 0.0, 0.0, 0, 0, 1, 1.25},
    {"bal_u2_16_load150", "balanced uniform [2,16], reject above 1.50x mean load", kStructureBalancedUniform, 2, 16, 0.0, 0.0, 0.0, 0, 0, 1, 1.50},
    {"bal_u2_16_load200", "balanced uniform [2,16], reject above 2.00x mean load", kStructureBalancedUniform, 2, 16, 0.0, 0.0, 0.0, 0, 0, 1, 2.00},
    {"bal_u2_16_best2", "balanced uniform [2,16], best of 2 row candidates", kStructureBalancedUniform, 2, 16, 0.0, 0.0, 0.0, 0, 0, 2, 2.0},
    {"bal_u2_16_best4", "balanced uniform [2,16], best of 4 row candidates", kStructureBalancedUniform, 2, 16, 0.0, 0.0, 0.0, 0, 0, 2, 4.0},
    {"deck_lt_m2_c256_fold", "deck-dealt exactly balanced columns, LT min 2 cap 256 folded tail degree law", kStructureDeckColumns, 2, 256, 0.0, 0.0, 0.0, 0, 0, 0, 1.0, kStructureLTFold},
    {"deck_lt_m2_c512_fold", "deck-dealt exactly balanced columns, LT min 2 cap 512 folded tail degree law", kStructureDeckColumns, 2, 512, 0.0, 0.0, 0.0, 0, 0, 0, 1.0, kStructureLTFold},
    {"deck_lt_m2_c128", "deck-dealt exactly balanced columns, LT min 2 cap 128 degree law", kStructureDeckColumns, 2, 128, 0.0, 0.0, 0.0, 0, 0, 0, 1.0, kStructureLT},
    {"deck_lt_m2_c16", "deck-dealt exactly balanced columns, LT min 2 cap 16 degree law", kStructureDeckColumns, 2, 16, 0.0, 0.0, 0.0, 0, 0, 0, 1.0, kStructureLT},
    {"deck_lt_m1_c64", "deck-dealt exactly balanced columns, LT min 1 cap 64 degree law", kStructureDeckColumns, 1, 64, 0.0, 0.0, 0.0, 0, 0, 0, 1.0, kStructureLT},
    {"deck_rs_c001_d50_c128", "deck-dealt exactly balanced columns, robust soliton c=0.01 delta=0.50 cap 128 degree law", kStructureDeckColumns, 1, 128, 0.01, 0.50, 0.0, 0, 0, 0, 1.0, kStructureRobustSoliton},
    {"deckloss10_lt_m2_c256_fold", "deck deal over 1.11x virtual rows then keep configured rows, LT min 2 cap 256 folded tail", kStructureDeckColumns, 2, 256, 0.0, 0.0, 0.0, 0, 0, 0, 1.11, kStructureLTFold},
    {"deckloss30_lt_m2_c256_fold", "deck deal over 1.43x virtual rows then keep configured rows, LT min 2 cap 256 folded tail", kStructureDeckColumns, 2, 256, 0.0, 0.0, 0.0, 0, 0, 0, 1.43, kStructureLTFold},
    {"p2c_lt_m2_c256_fold", "power-of-two-choices column balance, LT min 2 cap 256 folded tail degree law", kStructureLTFold, 2, 256, 0.0, 0.0, 0.0, 0, 0, 3, 2.0},
    {"harmonic1_32", "degree probability proportional to 1/d in [1,32]", kStructureHarmonic, 1, 32, 0.0, 0.0},
    {"harmonic2_32", "degree probability proportional to 1/d in [2,32]", kStructureHarmonic, 2, 32, 0.0, 0.0},
    {"one10_u2_8", "10 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.10, 0.0},
    {"one25_u2_8", "25 percent degree-1 rows, otherwise uniform [2,8]", kStructureOnePlusUniform, 2, 8, 0.25, 0.0},
    {"cap4_u2", "uniform [2,4] max-degree cap", kStructureUniform, 2, 4, 0.0, 0.0},
    {"cap8_u2", "uniform [2,8] max-degree cap", kStructureUniform, 2, 8, 0.0, 0.0},
    {"cap16_u2", "uniform [2,16] max-degree cap", kStructureUniform, 2, 16, 0.0, 0.0},
};

static const MatrixStructure* find_structure(const std::string& name)
{
    for (unsigned i = 0; i < sizeof(kStructures) / sizeof(kStructures[0]); ++i) {
        if (name == kStructures[i].Name) {
            return &kStructures[i];
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
    uint16_t TodoPos = kListTerm;
    uint8_t Mark = kMarkTodo;
};

struct Metrics
{
    unsigned ExactD2 = 0;
    unsigned Fill = 0;
    // 64-bit: dense structures at large N overflow a 32-bit sum of links^2
    uint64_t FillSquare = 0;
    unsigned Lookahead = 0;
    unsigned DistinctPartners = 0;
    unsigned DuplicatePartners = 0;
    unsigned RowRefs = 0;
    unsigned LiveRefs = 0;
    unsigned MinLive = 0xffffu;
    unsigned MaxLive = 0;
    unsigned Degree3Rows = 0;
    unsigned Degree4Rows = 0;
};

enum MetricFlags
{
    kMetricExactD2 = 1u << 0,
    kMetricFill = 1u << 1,
    kMetricFillSquare = 1u << 2,
    kMetricLookahead = 1u << 3,
    kMetricPartners = 1u << 4,
    kMetricLiveRefs = 1u << 5,
    kMetricMinMax = 1u << 6,
    kMetricDegree34 = 1u << 7,
};

static const unsigned kMetricAll =
    kMetricExactD2 | kMetricFill | kMetricFillSquare | kMetricLookahead |
    kMetricPartners | kMetricLiveRefs | kMetricMinMax | kMetricDegree34;

static unsigned normalize_metric_flags(unsigned flags)
{
    if (flags & kMetricFillSquare) {
        flags |= kMetricFill;
    }
    if (flags & (kMetricLookahead | kMetricPartners)) {
        flags |= kMetricExactD2;
    }
    return flags;
}

static unsigned legacy_metric_flags(int method_id)
{
    switch (method_id)
    {
    case 1:
        return kMetricExactD2;
    case 2:
        return kMetricExactD2 | kMetricLookahead;
    case 3:
        return kMetricExactD2 | kMetricLookahead | kMetricFill;
    case 4:
        return kMetricExactD2 | kMetricFill;
    case 5:
        return kMetricExactD2 | kMetricFill | kMetricFillSquare;
    case 6:
        return kMetricExactD2 | kMetricPartners;
    case 12:
        return kMetricExactD2 | kMetricLookahead | kMetricFill;
    case 13:
        return kMetricExactD2 | kMetricFill | kMetricMinMax;
    case 14:
        return kMetricExactD2 | kMetricLookahead | kMetricDegree34;
    default:
        return 0;
    }
}

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

static bool worse_key(const Key& a, const Key& b)
{
    return better_key(b, a);
}

struct TrialResult
{
    unsigned ActualN = 0;
    unsigned MatrixRows = 0;
    uint64_t MatrixRefs = 0;
    uint64_t MatrixXors = 0;
    unsigned ResidualCols = 0;
    unsigned ResidualRows = 0;
    unsigned Components = 0;
    unsigned MaxComponent = 0;
    unsigned SumSquares = 0;
    unsigned InitialTodo = 0;
    unsigned InitialD2Rows = 0;
    unsigned Choices = 0;
    uint64_t SparseSolveXors = 0;
    double DenseSolveXorsEst = 0.0;
    double TotalSolveXorsEst = 0.0;
    unsigned BuildUsec = 0;
    unsigned GreedyUsec = 0;
    unsigned TotalUsec = 0;
};

struct MatrixWork
{
    uint64_t Refs = 0;
    uint64_t Xors = 0;
};

struct TrialConfig
{
    unsigned ActualN = 0;
    unsigned MatrixRows = 0;
    uint64_t BaseSeed = 0;
    uint64_t MatrixSeed = 0;
    uint32_t ActionSeed = 0;
};

struct AvalancheEstimate
{
    unsigned Removed = 0;
    unsigned Peeled = 0;
    unsigned SingletonRows = 0;
};

class PeelSim
{
public:
    explicit PeelSim(unsigned block_count)
        : N(block_count),
          Columns(block_count)
    {
        TodoColumns.reserve(block_count);
        for (unsigned column_i = 0; column_i < N; ++column_i)
        {
            Columns[column_i].TodoPos = (uint16_t)column_i;
            TodoColumns.push_back((uint16_t)column_i);
        }
    }

    TrialResult run_greedy(int method_id, uint32_t seed)
    {
        if (is_multistart_method(method_id)) {
            return run_multistart(method_id, seed);
        }

        TrialResult result;
        result.InitialTodo = count_todo_columns();
        result.InitialD2Rows = count_live_rows(2);

        const double greedy_start = now_sec();
        for (;;)
        {
            const std::vector<uint16_t> action = select_inactivation_action(method_id, seed);
            if (action.empty()) {
                break;
            }

            if (!apply_inactivation_action(action)) {
                break;
            }
        }
        result.GreedyUsec = elapsed_usec(greedy_start);

        result.ResidualCols = DeferCount;
        result.ResidualRows = DeferredRows;
        result.Choices = ChoiceCount;
        result.SparseSolveXors = SparseSolveXors;
        measure_components(result.Components, result.MaxComponent, result.SumSquares);
        return result;
    }

    bool run_greedy_checked(
        int method_id,
        uint32_t seed,
        TrialResult* out,
        std::string* why)
    {
        TrialResult result;
        if (is_multistart_method(method_id)) {
            return run_multistart_checked(method_id, seed, out, why);
        }

        result.InitialTodo = count_todo_columns();
        result.InitialD2Rows = count_live_rows(2);

        const double greedy_start = now_sec();
        for (;;)
        {
            std::string detail;
            if (!validate_next_action(method_id, seed, &detail)) {
                return fail(why, detail.c_str());
            }

            const std::vector<uint16_t> action =
                select_inactivation_action(method_id, seed);
            if (action.empty()) {
                break;
            }
            if (!apply_inactivation_action(action)) {
                return fail(why, "validated action did not apply");
            }
            if (!validate_state(false, &detail)) {
                return fail(why, detail.c_str());
            }
        }
        result.GreedyUsec = elapsed_usec(greedy_start);

        result.ResidualCols = DeferCount;
        result.ResidualRows = DeferredRows;
        result.Choices = ChoiceCount;
        result.SparseSolveXors = SparseSolveXors;
        measure_components(result.Components, result.MaxComponent, result.SumSquares);
        if (out) {
            *out = result;
        }
        return validate_result(result, why);
    }

    bool run_greedy_traced(
        int method_id,
        uint32_t seed,
        TrialResult* out,
        std::vector<uint16_t>* trace,
        std::string* why)
    {
        if (is_multistart_method(method_id)) {
            return fail(why, "traced runs do not support multistart methods");
        }

        TrialResult result;
        result.InitialTodo = count_todo_columns();
        result.InitialD2Rows = count_live_rows(2);

        const double greedy_start = now_sec();
        for (;;)
        {
            const std::vector<uint16_t> action =
                select_inactivation_action(method_id, seed);
            if (action.empty()) {
                break;
            }
            if (trace) {
                trace->insert(trace->end(), action.begin(), action.end());
            }
            if (!apply_inactivation_action(action)) {
                return fail(why, "traced action did not apply");
            }
        }
        result.GreedyUsec = elapsed_usec(greedy_start);

        result.ResidualCols = DeferCount;
        result.ResidualRows = DeferredRows;
        result.Choices = ChoiceCount;
        result.SparseSolveXors = SparseSolveXors;
        measure_components(result.Components, result.MaxComponent,
            result.SumSquares);
        if (out) {
            *out = result;
        }
        return true;
    }

    bool validate_lazy_equivalence(uint32_t seed, std::string* why) const
    {
        for (uint16_t column_i : TodoColumns)
        {
            const Metrics exact = collect_metrics(column_i);
            const Metrics cached = cached_metrics(column_i, kMetricAll);
            if (exact.LiveRefs != cached.LiveRefs ||
                exact.Fill != cached.Fill ||
                exact.FillSquare != cached.FillSquare ||
                exact.ExactD2 != cached.ExactD2 ||
                exact.Lookahead != cached.Lookahead ||
                exact.DistinctPartners != cached.DistinctPartners ||
                exact.DuplicatePartners != cached.DuplicatePartners ||
                exact.Degree3Rows != cached.Degree3Rows ||
                exact.Degree4Rows != cached.Degree4Rows ||
                exact.MinLive != cached.MinLive ||
                exact.MaxLive != cached.MaxLive ||
                exact.RowRefs != cached.RowRefs) {
                return fail(why, "incremental metric cache differs from exact scan");
            }
        }

        const Method* liveref = find_method(58);
        if (!liveref) {
            return fail(why, "global live-ref method is missing");
        }
        if (select_column(112, seed) != select_combo(*liveref, seed)) {
            return fail(why, "lazy live-ref selection differs from exact scan");
        }
        if (select_column(113, seed) != select_by_full_scan(4, seed)) {
            return fail(why, "lazy fill selection differs from exact scan");
        }
        if (select_column(114, seed) != select_by_full_scan(1, seed)) {
            return fail(why, "lazy degree-2 selection differs from exact scan");
        }
        return true;
    }

    bool run_lazy_equivalence_walk(uint32_t seed, std::string* why)
    {
        for (;;)
        {
            std::string detail;
            if (!validate_lazy_equivalence(seed, &detail)) {
                return fail(why, detail.c_str());
            }

            const uint16_t column_i = select_by_full_scan(0, seed);
            if (column_i == kListTerm) {
                break;
            }

            std::vector<uint16_t> action(1, column_i);
            if (!apply_inactivation_action(action)) {
                return fail(why, "lazy equivalence walk action did not apply");
            }
            if (!validate_state(false, &detail)) {
                return fail(why, detail.c_str());
            }
        }
        return validate_state(true, why);
    }

    bool validate_avalanche_estimate(uint32_t seed, std::string* why) const
    {
        const std::vector<uint16_t> candidates = top_default_candidates(32, seed);
        for (uint16_t column_i : candidates)
        {
            const AvalancheEstimate estimate = estimate_avalanche(column_i);
            PeelSim copy = *this;
            const unsigned before = copy.count_todo_columns();
            std::vector<uint16_t> action(1, column_i);
            if (!copy.apply_inactivation_action(action)) {
                return fail(why, "avalanche validation action did not apply");
            }
            const unsigned exact_removed = before - copy.count_todo_columns();
            if (estimate.Removed != exact_removed) {
                return fail(why, "local avalanche estimate differs from exact one-step simulation");
            }
        }
        return true;
    }

    bool run_avalanche_estimate_walk(
        int method_id,
        uint32_t seed,
        unsigned steps,
        std::string* why)
    {
        for (unsigned step = 0; step < steps; ++step)
        {
            std::string detail;
            if (!validate_avalanche_estimate(seed + step, &detail)) {
                return fail(why, detail.c_str());
            }
            const uint16_t column_i = select_column(method_id, seed + step);
            if (column_i == kListTerm) {
                break;
            }
            std::vector<uint16_t> action(1, column_i);
            if (!apply_inactivation_action(action)) {
                return fail(why, "avalanche estimate walk action did not apply");
            }
            if (!validate_state(false, &detail)) {
                return fail(why, detail.c_str());
            }
        }
        return true;
    }

    bool add_test_row(const std::vector<uint16_t>& columns, std::string* why)
    {
        if (columns.empty()) {
            return fail(why, "test row has no columns");
        }
        if (Rows.size() >= 65535u) {
            return fail(why, "test row count exceeds 16-bit model limit");
        }

        std::vector<uint8_t> seen(N, 0);
        for (uint16_t column_i : columns)
        {
            if (column_i >= N) {
                return fail(why, "test row column is out of range");
            }
            if (seen[column_i]) {
                return fail(why, "test row contains a duplicate column");
            }
            seen[column_i] = 1;
        }

        Row row;
        uint16_t unmarked[2] = {kListTerm, kListTerm};
        unsigned unmarked_count = 0;
        const uint16_t row_i = (uint16_t)Rows.size();

        for (uint16_t column_i : columns)
        {
            row.Columns.push_back(column_i);
            Columns[column_i].Rows.push_back(row_i);
            if (Columns[column_i].Mark == kMarkTodo)
            {
                unmarked[unmarked_count & 1u] = column_i;
                ++unmarked_count;
            }
            else if (Columns[column_i].Mark == kMarkPeel) {
                ++SparseSolveXors;
            }
        }

        row.Live = (uint16_t)unmarked_count;
        Rows.push_back(row);

        if (unmarked_count == 0) {
            add_deferred_row(row_i);
        }
        else if (unmarked_count == 1) {
            solve_with_peel(row_i, unmarked[0]);
        }
        else if (unmarked_count == 2)
        {
            increment_weight2_ref(unmarked[0]);
            increment_weight2_ref(unmarked[1]);
        }

        invalidate_metric_cache();
        return validate_state(false, why);
    }

    bool feed_columns(const std::vector<uint16_t>& columns, std::string* why)
    {
        if (columns.empty()) {
            return fail(why, "generated row has no columns");
        }
        if (Rows.size() >= 65535u) {
            return fail(why, "generated row count exceeds 16-bit model limit");
        }

        for (unsigned i = 0; i < columns.size(); ++i)
        {
            const uint16_t column_i = columns[i];
            if (column_i >= N) {
                return fail(why, "generated row column is out of range");
            }
            for (unsigned j = 0; j < i; ++j)
            {
                if (columns[j] == column_i) {
                    return fail(why, "generated row contains a duplicate column");
                }
            }
        }

        Row row;
        uint16_t unmarked[2] = {kListTerm, kListTerm};
        unsigned unmarked_count = 0;
        const uint16_t row_i = (uint16_t)Rows.size();

        for (uint16_t column_i : columns)
        {
            row.Columns.push_back(column_i);
            Columns[column_i].Rows.push_back(row_i);
            if (Columns[column_i].Mark == kMarkTodo)
            {
                unmarked[unmarked_count & 1u] = column_i;
                ++unmarked_count;
            }
            else if (Columns[column_i].Mark == kMarkPeel) {
                ++SparseSolveXors;
            }
        }

        row.Live = (uint16_t)unmarked_count;
        Rows.push_back(row);

        if (unmarked_count == 0) {
            add_deferred_row(row_i);
        }
        else if (unmarked_count == 1) {
            solve_with_peel(row_i, unmarked[0]);
        }
        else if (unmarked_count == 2)
        {
            increment_weight2_ref(unmarked[0]);
            increment_weight2_ref(unmarked[1]);
        }
        invalidate_metric_cache();
        return true;
    }

    bool validate_state(bool require_complete, std::string* why) const
    {
        unsigned counted_deferred_cols = 0;
        unsigned counted_deferred_rows = 0;
        unsigned counted_choices = 0;
        unsigned counted_todo_cols = 0;
        std::vector<uint8_t> seen_todo(N, 0);

        for (uint16_t column_i = 0; column_i < N; ++column_i)
        {
            const Column& column = Columns[column_i];
            if (column.Mark > kMarkDefer) {
                return fail(why, "column has an invalid mark");
            }
            if (column.Mark == kMarkDefer) {
                ++counted_deferred_cols;
            }
            if (column.Mark == kMarkTodo) {
                ++counted_todo_cols;
            }

            for (uint16_t row_i : column.Rows)
            {
                if (row_i >= Rows.size()) {
                    return fail(why, "column references an out-of-range row");
                }
                if (!row_has_column(Rows[row_i], column_i)) {
                    return fail(why, "column-row back reference is missing");
                }
            }
        }

        if (TodoColumns.size() != counted_todo_cols) {
            return fail(why, "todo column list size is stale");
        }
        for (unsigned pos = 0; pos < TodoColumns.size(); ++pos)
        {
            const uint16_t column_i = TodoColumns[pos];
            if (column_i >= N) {
                return fail(why, "todo column list contains an out-of-range column");
            }
            const Column& column = Columns[column_i];
            if (column.Mark != kMarkTodo) {
                return fail(why, "todo column list contains a non-todo column");
            }
            if (column.TodoPos != pos) {
                return fail(why, "todo column position is stale");
            }
            if (seen_todo[column_i]) {
                return fail(why, "todo column list contains a duplicate column");
            }
            seen_todo[column_i] = 1;
        }
        for (uint16_t column_i = 0; column_i < N; ++column_i)
        {
            const Column& column = Columns[column_i];
            if (column.Mark == kMarkTodo && !seen_todo[column_i]) {
                return fail(why, "todo column is missing from todo list");
            }
            if (column.Mark != kMarkTodo && column.TodoPos != kListTerm) {
                return fail(why, "non-todo column has a todo position");
            }
        }

        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Solved && row.Deferred) {
                return fail(why, "row is both solved and deferred");
            }
            if (row.Solved)
            {
                if (row.PeelColumn == kListTerm || row.PeelColumn >= N) {
                    return fail(why, "solved row has an invalid peel column");
                }
                if (Columns[row.PeelColumn].Mark != kMarkPeel) {
                    return fail(why, "solved row peel column is not marked peeled");
                }
                if (!row_has_column(row, row.PeelColumn)) {
                    return fail(why, "solved row does not reference its peel column");
                }
            }
            if (row.Deferred) {
                ++counted_deferred_rows;
            }

            unsigned live = 0;
            for (uint16_t column_i : row.Columns)
            {
                if (column_i >= N) {
                    return fail(why, "row references an out-of-range column");
                }
                if (!column_has_row(Columns[column_i], row_i)) {
                    return fail(why, "row-column back reference is missing");
                }
                if (Columns[column_i].Mark == kMarkTodo) {
                    ++live;
                }
            }
            if (row.Live != live) {
                return fail(why, "row live count is stale");
            }
        }

        counted_choices = counted_deferred_cols;
        if (counted_deferred_cols != DeferCount) {
            return fail(why, "deferred column count is stale");
        }
        if (counted_deferred_rows != DeferredRows) {
            return fail(why, "deferred row count is stale");
        }
        if (counted_choices != ChoiceCount) {
            return fail(why, "choice count does not match deferred columns");
        }
        if (require_complete && count_todo_columns() != 0) {
            return fail(why, "greedy run left todo columns");
        }
        return true;
    }

    bool validate_next_action(int method_id, uint32_t seed, std::string* why) const
    {
        if (!find_method(method_id)) {
            return fail(why, "unknown method id");
        }

        const std::vector<uint16_t> action =
            select_inactivation_action(method_id, seed);
        const unsigned todo_count = count_todo_columns();
        if (todo_count == 0) {
            return action.empty() ? true :
                fail(why, "method returned an action with no todo columns");
        }
        if (action.empty()) {
            return fail(why, "method returned no action while todo columns remain");
        }

        std::vector<uint8_t> seen(N, 0);
        for (uint16_t column_i : action)
        {
            if (column_i >= N) {
                return fail(why, "action column is out of range");
            }
            if (Columns[column_i].Mark != kMarkTodo) {
                return fail(why, "action selected a non-todo column");
            }
            if (seen[column_i]) {
                return fail(why, "action contains a duplicate column");
            }
            seen[column_i] = 1;
        }

        const bool row_action =
            (method_id >= 66 && method_id <= 72) ||
            (method_id >= 106 && method_id <= 111);
        if (!row_action && action.size() != 1) {
            return fail(why, "singleton method returned multiple columns");
        }
        if (row_action && action.size() > 1 &&
            !action_matches_all_but_one_row(action)) {
            return fail(why, "row-action method did not choose all-but-one row columns");
        }

        if (action.size() == 1) {
            return validate_singleton_pool(method_id, action[0], seed, why);
        }
        return true;
    }

    bool validate_result(const TrialResult& result, std::string* why) const
    {
        if (!validate_state(true, why)) {
            return false;
        }
        if (result.ResidualCols != DeferCount) {
            return fail(why, "result residual column count is stale");
        }
        if (result.ResidualRows != DeferredRows) {
            return fail(why, "result residual row count is stale");
        }
        if (result.Choices != ChoiceCount) {
            return fail(why, "result choice count is stale");
        }
        if (result.ResidualCols == 0)
        {
            if (result.Components != 0 || result.MaxComponent != 0 ||
                result.SumSquares != 0) {
                return fail(why, "zero residual result has nonzero components");
            }
        }
        else
        {
            if (result.Components == 0) {
                return fail(why, "nonzero residual result has zero components");
            }
            if (result.MaxComponent == 0 ||
                result.MaxComponent > result.ResidualCols) {
                return fail(why, "invalid maximum component size");
            }
            if (result.SumSquares < result.ResidualCols) {
                return fail(why, "invalid component square sum");
            }
        }
        return true;
    }

    static bool validate_multistart_tie_break(std::string* why)
    {
        TrialResult stable;
        stable.ResidualCols = 10;
        stable.ResidualRows = 14;
        stable.Components = 4;
        stable.MaxComponent = 3;
        stable.SumSquares = 28;
        stable.Choices = 10;
        stable.GreedyUsec = 1000;

        TrialResult fast = stable;
        fast.Components = 3;
        fast.MaxComponent = 4;
        fast.Choices = 11;
        fast.GreedyUsec = 1;

        if (!better_trial_result(stable, fast)) {
            return fail(why, "multistart tie-break preferred elapsed time over structure");
        }
        if (better_trial_result(fast, stable)) {
            return fail(why, "multistart tie-break is not antisymmetric");
        }

        TrialResult equal = stable;
        equal.GreedyUsec = 1;
        if (better_trial_result(stable, equal) ||
            better_trial_result(equal, stable)) {
            return fail(why, "multistart tie-break depends on elapsed time");
        }
        return true;
    }

    static bool validate_minrow_random_tie(std::string* why)
    {
        PeelSim sim(4);
        std::string detail;
        if (!sim.add_test_row(std::vector<uint16_t>{0, 1}, &detail) ||
            !sim.add_test_row(std::vector<uint16_t>{0, 2}, &detail)) {
            return fail(why, detail.c_str());
        }

        const std::vector<uint16_t> action =
            sim.select_inactivation_action(88, UINT32_C(6));
        if (action.size() != 1 || action[0] != 2) {
            return fail(why, "minrow_rand weighted duplicate columns");
        }
        return true;
    }

private:
    unsigned N;
    std::vector<Row> Rows;
    std::vector<Column> Columns;
    std::vector<uint16_t> TodoColumns;
    unsigned DeferCount = 0;
    unsigned DeferredRows = 0;
    unsigned ChoiceCount = 0;
    uint64_t SparseSolveXors = 0;
    mutable std::vector<Metrics> MetricCache;
    mutable std::vector<std::vector<std::pair<uint16_t, unsigned> > > PartnerCounts;
    mutable std::vector<uint16_t> LiveDegreeCounts;
    mutable unsigned LiveDegreeStride = 0;
    mutable bool MetricCacheValid = false;
    mutable unsigned MetricCacheFlags = 0;
    mutable std::vector<AvalancheEstimate> AvalancheCacheValue;
    mutable std::vector<uint8_t> AvalancheCacheValid;
    mutable std::vector<std::vector<uint16_t> > AvalancheRowDependents;
    mutable std::vector<uint32_t> AvalancheRemovedEpoch;
    mutable std::vector<uint32_t> AvalancheRowEpoch;
    mutable std::vector<uint16_t> AvalancheQueueScratch;
    mutable std::vector<uint16_t> AvalancheDepScratch;
    mutable uint32_t AvalancheEpoch = 0;
    mutable bool AvalancheCacheActive = false;

    static unsigned elapsed_usec(double start)
    {
        const double usec = (now_sec() - start) * 1000000.0;
        return usec > 0.0 ? (unsigned)(usec + 0.5) : 0;
    }

    static bool is_multistart_method(int method_id)
    {
        return method_id >= 100 && method_id <= 105;
    }

    void invalidate_metric_cache() const
    {
        MetricCacheValid = false;
        MetricCacheFlags = 0;
        LiveDegreeStride = 0;
        invalidate_avalanche_cache();
    }

    void invalidate_avalanche_cache() const
    {
        AvalancheCacheActive = false;
        AvalancheCacheValue.clear();
        AvalancheCacheValid.clear();
        AvalancheRowDependents.clear();
        AvalancheRemovedEpoch.clear();
        AvalancheRowEpoch.clear();
        AvalancheEpoch = 0;
    }

    uint32_t next_avalanche_epoch() const
    {
        // Epoch stamps avoid re-zeroing O(N) scratch on each cache miss.
        if (AvalancheRemovedEpoch.size() != N) {
            AvalancheRemovedEpoch.assign(N, 0);
        }
        if (AvalancheRowEpoch.size() != Rows.size()) {
            AvalancheRowEpoch.assign(Rows.size(), 0);
        }
        if (++AvalancheEpoch == 0)
        {
            std::fill(AvalancheRemovedEpoch.begin(),
                AvalancheRemovedEpoch.end(), 0u);
            std::fill(AvalancheRowEpoch.begin(),
                AvalancheRowEpoch.end(), 0u);
            AvalancheEpoch = 1;
        }
        return AvalancheEpoch;
    }

    void ensure_avalanche_cache() const
    {
        if (AvalancheCacheActive &&
            AvalancheCacheValue.size() == N &&
            AvalancheRowDependents.size() == Rows.size()) {
            return;
        }
        AvalancheCacheValue.assign(N, AvalancheEstimate());
        AvalancheCacheValid.assign(N, 0);
        AvalancheRowDependents.assign(
            Rows.size(), std::vector<uint16_t>());
        AvalancheCacheActive = true;
    }

    void dirty_avalanche_row(uint16_t row_i) const
    {
        if (row_i >= AvalancheRowDependents.size()) {
            return;
        }
        std::vector<uint16_t>& dependents = AvalancheRowDependents[row_i];
        for (uint16_t column_i : dependents)
        {
            if (column_i < AvalancheCacheValid.size()) {
                AvalancheCacheValid[column_i] = 0;
            }
        }
        dependents.clear();
    }

    void dirty_avalanche_column(uint16_t column_i) const
    {
        if (column_i < AvalancheCacheValid.size()) {
            AvalancheCacheValid[column_i] = 0;
        }
    }

    void mark_column(uint16_t column_i, uint8_t mark)
    {
        Column& column = Columns[column_i];
        if (column.Mark == mark) {
            return;
        }

        if (column.Mark == kMarkTodo)
        {
            const uint16_t pos = column.TodoPos;
            if (pos < TodoColumns.size())
            {
                const uint16_t moved = TodoColumns.back();
                TodoColumns[pos] = moved;
                Columns[moved].TodoPos = pos;
                TodoColumns.pop_back();
            }
            column.TodoPos = kListTerm;
        }
        else if (mark == kMarkTodo)
        {
            column.TodoPos = (uint16_t)TodoColumns.size();
            TodoColumns.push_back(column_i);
        }

        column.Mark = mark;
    }

    static bool better_trial_result(
        const TrialResult& a,
        const TrialResult& b)
    {
        if (a.ResidualCols != b.ResidualCols) {
            return a.ResidualCols < b.ResidualCols;
        }
        if (a.ResidualRows != b.ResidualRows) {
            return a.ResidualRows < b.ResidualRows;
        }
        if (a.SumSquares != b.SumSquares) {
            return a.SumSquares < b.SumSquares;
        }
        if (a.MaxComponent != b.MaxComponent) {
            return a.MaxComponent < b.MaxComponent;
        }
        if (a.Components != b.Components) {
            return a.Components > b.Components;
        }
        if (a.Choices != b.Choices) {
            return a.Choices < b.Choices;
        }
        return false;
    }

    static void multistart_plan(
        int method_id,
        int* methods,
        unsigned* method_count)
    {
        *method_count = 0;
        switch (method_id)
        {
        case 100:
            methods[(*method_count)++] = 7;
            methods[(*method_count)++] = 7;
            break;
        case 101:
            methods[(*method_count)++] = 65;
            methods[(*method_count)++] = 65;
            break;
        case 102:
            for (unsigned i = 0; i < 4; ++i) {
                methods[(*method_count)++] = 7;
            }
            break;
        case 103:
            for (unsigned i = 0; i < 4; ++i) {
                methods[(*method_count)++] = 65;
            }
            break;
        case 104:
            methods[(*method_count)++] = 0;
            methods[(*method_count)++] = 7;
            methods[(*method_count)++] = 7;
            methods[(*method_count)++] = 7;
            break;
        case 105:
            for (unsigned i = 0; i < 4; ++i) {
                methods[(*method_count)++] = 79;
            }
            break;
        default:
            break;
        }
    }

    TrialResult run_multistart(int method_id, uint32_t seed) const
    {
        int methods[4] = {0, 0, 0, 0};
        unsigned method_count = 0;
        multistart_plan(method_id, methods, &method_count);
        if (method_count == 0) {
            PeelSim copy = *this;
            return copy.run_greedy(0, seed);
        }

        TrialResult best;
        bool have_best = false;
        unsigned long long greedy_sum = 0;
        for (unsigned i = 0; i < method_count; ++i)
        {
            PeelSim copy = *this;
            const uint32_t run_seed = hash32(seed ^
                ((uint32_t)i * UINT32_C(0x9e3779b9)) ^
                ((uint32_t)method_id * UINT32_C(0x85ebca6b)));
            TrialResult result = copy.run_greedy(methods[i], run_seed);
            greedy_sum += result.GreedyUsec;
            if (!have_best || better_trial_result(result, best))
            {
                best = result;
                have_best = true;
            }
        }
        best.GreedyUsec = (unsigned)std::min<unsigned long long>(
            greedy_sum, (unsigned long long)UINT_MAX);
        return best;
    }

    bool run_multistart_checked(
        int method_id,
        uint32_t seed,
        TrialResult* out,
        std::string* why) const
    {
        int methods[4] = {0, 0, 0, 0};
        unsigned method_count = 0;
        multistart_plan(method_id, methods, &method_count);
        if (method_count == 0) {
            return fail(why, "multistart method has no child methods");
        }

        TrialResult best;
        bool have_best = false;
        unsigned long long greedy_sum = 0;
        for (unsigned i = 0; i < method_count; ++i)
        {
            PeelSim copy = *this;
            TrialResult result;
            std::string detail;
            const uint32_t run_seed = hash32(seed ^
                ((uint32_t)i * UINT32_C(0x9e3779b9)) ^
                ((uint32_t)method_id * UINT32_C(0x85ebca6b)));
            if (!copy.run_greedy_checked(methods[i], run_seed, &result, &detail))
            {
                const std::string message =
                    std::string("multistart child failed: ") + detail;
                return fail(why, message.c_str());
            }
            greedy_sum += result.GreedyUsec;
            if (!have_best || better_trial_result(result, best))
            {
                best = result;
                have_best = true;
            }
        }
        best.GreedyUsec = (unsigned)std::min<unsigned long long>(
            greedy_sum, (unsigned long long)UINT_MAX);
        if (out) {
            *out = best;
        }
        return true;
    }

    static bool fail(std::string* why, const char* message)
    {
        if (why) {
            *why = message;
        }
        return false;
    }

    bool column_has_row(const Column& column, uint16_t row_i) const
    {
        for (uint16_t r : column.Rows) {
            if (r == row_i) {
                return true;
            }
        }
        return false;
    }

    bool apply_inactivation_action(const std::vector<uint16_t>& action)
    {
        std::vector<uint16_t> applied_columns;
        for (uint16_t column_i : action)
        {
            if (column_i == kListTerm) {
                continue;
            }
            Column& column = Columns[column_i];
            if (column.Mark != kMarkTodo) {
                continue;
            }
            mark_column(column_i, kMarkDefer);
            ++DeferCount;
            ++ChoiceCount;
            applied_columns.push_back(column_i);
        }
        if (applied_columns.empty()) {
            return false;
        }

        for (uint16_t column_i : applied_columns)
        {
            avalanche(column_i);
        }
        return true;
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
        mark_column(column_i, kMarkPeel);
        avalanche(column_i);
    }

    void avalanche(uint16_t first_column)
    {
        std::vector<uint16_t> queue;
        queue.push_back(first_column);

        for (unsigned queue_i = 0; queue_i < queue.size(); ++queue_i)
        {
            const uint16_t column_i = queue[queue_i];
            const bool column_is_peel =
                Columns[column_i].Mark == kMarkPeel;
            const std::vector<uint16_t>& refs = Columns[column_i].Rows;
            if (AvalancheCacheActive) {
                dirty_avalanche_column(column_i);
            }

            for (uint16_t row_i : refs)
            {
                Row& row = Rows[row_i];
                if (row.Live == 0) {
                    continue;
                }
                if (column_is_peel &&
                    !(row.Solved && row.PeelColumn == column_i)) {
                    ++SparseSolveXors;
                }

                adjust_row_metric_contribution(row_i, -1, column_i);
                if (AvalancheCacheActive) {
                    dirty_avalanche_row(row_i);
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
                            mark_column(found[0], kMarkPeel);
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
                        increment_weight2_ref(found[0], row_i);
                        increment_weight2_ref(found[1], row_i);
                    }
                    else if (count == 1)
                    {
                        if (Columns[found[0]].Mark == kMarkTodo)
                        {
                            row.Solved = true;
                            row.PeelColumn = found[0];
                            mark_column(found[0], kMarkPeel);
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
                adjust_row_metric_contribution(row_i, 1);
            }
        }
    }

    unsigned count_todo_columns() const
    {
        return (unsigned)TodoColumns.size();
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
            ++m.LiveRefs;

            if (row.Live > 1)
            {
                const unsigned links = row.Live - 1;
                m.Fill += links;
                m.FillSquare += (uint64_t)links * links;
                if (row.Live < m.MinLive) {
                    m.MinLive = row.Live;
                }
                if (row.Live > m.MaxLive) {
                    m.MaxLive = row.Live;
                }
                if (row.Live == 3) {
                    ++m.Degree3Rows;
                }
                else if (row.Live == 4) {
                    ++m.Degree4Rows;
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

    unsigned collect_metric_live_pair(
        const Row& row,
        uint16_t extra_column,
        uint16_t* live) const
    {
        unsigned count = 0;
        for (uint16_t column_i : row.Columns)
        {
            if (Columns[column_i].Mark == kMarkTodo)
            {
                if (count < 2) {
                    live[count] = column_i;
                }
                ++count;
            }
        }
        if (extra_column != kListTerm && row_has_column(row, extra_column))
        {
            bool present = false;
            for (unsigned i = 0; i < count && i < 2; ++i)
            {
                if (live[i] == extra_column) {
                    present = true;
                }
            }
            if (!present)
            {
                if (count < 2) {
                    live[count] = extra_column;
                }
                ++count;
            }
        }
        return count;
    }

    void recompute_metric_minmax(uint16_t column_i) const
    {
        Metrics& m = MetricCache[column_i];
        m.MinLive = 0xffffu;
        m.MaxLive = 0;
        for (uint16_t row_i : Columns[column_i].Rows)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred || row.Live <= 1) {
                continue;
            }
            if (row.Live < m.MinLive) {
                m.MinLive = row.Live;
            }
            if (row.Live > m.MaxLive) {
                m.MaxLive = row.Live;
            }
        }
        if (m.MinLive == 0xffffu) {
            m.MinLive = 0;
        }
    }

    void adjust_partner_count(uint16_t column_i, uint16_t partner, int delta) const
    {
        if (partner == kListTerm || column_i >= PartnerCounts.size()) {
            return;
        }

        Metrics& m = MetricCache[column_i];
        std::vector<std::pair<uint16_t, unsigned> >& counts =
            PartnerCounts[column_i];

        for (unsigned i = 0; i < counts.size(); ++i)
        {
            if (counts[i].first != partner) {
                continue;
            }

            if (delta > 0)
            {
                ++counts[i].second;
                ++m.DuplicatePartners;
            }
            else
            {
                if (counts[i].second <= 1)
                {
                    if (m.DistinctPartners > 0) {
                        --m.DistinctPartners;
                    }
                    counts[i] = counts.back();
                    counts.pop_back();
                }
                else
                {
                    --counts[i].second;
                    if (m.DuplicatePartners > 0) {
                        --m.DuplicatePartners;
                    }
                }
            }
            return;
        }

        if (delta > 0)
        {
            counts.push_back(std::make_pair(partner, 1u));
            ++m.DistinctPartners;
        }
    }

    void adjust_minmax_count(uint16_t column_i, unsigned live, int delta) const
    {
        if (live <= 1 || LiveDegreeStride == 0 ||
            live >= LiveDegreeStride || column_i >= N) {
            return;
        }

        Metrics& m = MetricCache[column_i];
        uint16_t* counts = &LiveDegreeCounts[(size_t)column_i * LiveDegreeStride];
        if (delta > 0)
        {
            if (counts[live] != UINT16_MAX) {
                ++counts[live];
            }
            // MinLive == 0 is the empty sentinel set by the decrement path;
            // a plain numeric compare would leave it stuck at 0 forever
            if (m.MinLive == 0 || live < m.MinLive) {
                m.MinLive = live;
            }
            if (live > m.MaxLive) {
                m.MaxLive = live;
            }
            return;
        }

        if (counts[live] > 0) {
            --counts[live];
        }
        if (counts[live] != 0) {
            return;
        }

        if (m.MinLive == live)
        {
            unsigned next = 0;
            for (unsigned d = live; d < LiveDegreeStride; ++d)
            {
                if (counts[d] != 0)
                {
                    next = d;
                    break;
                }
            }
            m.MinLive = next == 0 ? 0 : next;
        }

        if (m.MaxLive == live)
        {
            unsigned next = 0;
            for (unsigned d = std::min<unsigned>(live, LiveDegreeStride - 1u);;)
            {
                if (counts[d] != 0)
                {
                    next = d;
                    break;
                }
                if (d == 2u) {
                    break;
                }
                --d;
            }
            m.MaxLive = next;
        }
    }

    void adjust_row_metric_contribution(
        uint16_t row_i,
        int delta,
        uint16_t extra_column = kListTerm) const
    {
        if (!MetricCacheValid) {
            return;
        }

        const Row& row = Rows[row_i];
        if (row.Deferred || row.Live == 0) {
            return;
        }

        const unsigned flags = MetricCacheFlags;
        const bool need_live_refs = (flags & kMetricLiveRefs) != 0;
        const bool need_fill = (flags & kMetricFill) != 0;
        const bool need_fill_square = (flags & kMetricFillSquare) != 0;
        const bool need_lookahead = (flags & kMetricLookahead) != 0;
        const bool need_partners = (flags & kMetricPartners) != 0;
        const bool need_exact_d2 = (flags & kMetricExactD2) != 0;
        const bool need_minmax = (flags & kMetricMinMax) != 0;
        const bool need_degree34 = (flags & kMetricDegree34) != 0;

        const unsigned fill = row.Live > 1 ? row.Live - 1u : 0u;
        const uint64_t fill_square = (uint64_t)fill * fill;
        uint16_t live_pair[2] = {kListTerm, kListTerm};
        const bool need_pair = need_lookahead || need_partners;
        const unsigned live_pair_count = row.Live == 2 && need_pair ?
            collect_metric_live_pair(row, extra_column, live_pair) : 0;
        for (uint16_t column_i : row.Columns)
        {
            if (Columns[column_i].Mark != kMarkTodo) {
                continue;
            }

            Metrics& m = MetricCache[column_i];
            if (delta > 0)
            {
                if (need_live_refs) {
                    ++m.LiveRefs;
                }
                if (need_fill) {
                    m.Fill += fill;
                }
                if (need_fill_square) {
                    m.FillSquare += fill_square;
                }
                if (need_minmax) {
                    adjust_minmax_count(column_i, row.Live, delta);
                }
                if (row.Live == 2) {
                    if (need_exact_d2) {
                        ++m.ExactD2;
                    }
                    if (live_pair_count == 2)
                    {
                        const uint16_t partner =
                            live_pair[0] == column_i ? live_pair[1] :
                            live_pair[1] == column_i ? live_pair[0] : kListTerm;
                        if (partner != kListTerm)
                        {
                            if (need_lookahead) {
                                m.Lookahead += Columns[partner].Weight2Refs;
                            }
                            if (need_partners) {
                                adjust_partner_count(column_i, partner, delta);
                            }
                        }
                    }
                }
                else if (row.Live == 3 && need_degree34) {
                    ++m.Degree3Rows;
                }
                else if (row.Live == 4 && need_degree34) {
                    ++m.Degree4Rows;
                }
            }
            else
            {
                if (need_live_refs && m.LiveRefs > 0) {
                    --m.LiveRefs;
                }
                if (need_fill) {
                    m.Fill = m.Fill >= fill ? m.Fill - fill : 0;
                }
                if (need_fill_square) {
                    m.FillSquare = m.FillSquare >= fill_square ?
                        m.FillSquare - fill_square : 0;
                }
                if (row.Live == 2)
                {
                    if (need_exact_d2 && m.ExactD2 > 0) {
                        --m.ExactD2;
                    }
                    if (live_pair_count == 2)
                    {
                        const uint16_t partner =
                            live_pair[0] == column_i ? live_pair[1] :
                            live_pair[1] == column_i ? live_pair[0] : kListTerm;
                        if (partner != kListTerm)
                        {
                            if (need_lookahead)
                            {
                                const unsigned look = Columns[partner].Weight2Refs;
                                m.Lookahead = m.Lookahead >= look ?
                                    m.Lookahead - look : 0;
                            }
                            if (need_partners) {
                                adjust_partner_count(column_i, partner, delta);
                            }
                        }
                    }
                }
                else if (need_degree34 && row.Live == 3 && m.Degree3Rows > 0) {
                    --m.Degree3Rows;
                }
                else if (need_degree34 && row.Live == 4 && m.Degree4Rows > 0) {
                    --m.Degree4Rows;
                }
                if (need_minmax) {
                    adjust_minmax_count(column_i, row.Live, delta);
                }
            }
        }
    }

    void build_metric_cache(unsigned flags) const
    {
        MetricCacheFlags = normalize_metric_flags(flags);
        MetricCache.assign(N, Metrics());
        if (MetricCacheFlags & kMetricPartners) {
            PartnerCounts.assign(
                N, std::vector<std::pair<uint16_t, unsigned> >());
        }
        else {
            PartnerCounts.clear();
        }
        if (MetricCacheFlags & kMetricMinMax)
        {
            unsigned max_live = 1;
            for (const Row& row : Rows)
            {
                if (row.Columns.size() > max_live) {
                    max_live = (unsigned)row.Columns.size();
                }
            }
            LiveDegreeStride = max_live + 1u;
            LiveDegreeCounts.assign((size_t)N * LiveDegreeStride, 0);
        }
        else
        {
            LiveDegreeCounts.clear();
            LiveDegreeStride = 0;
        }
        for (uint16_t column_i = 0; column_i < N; ++column_i) {
            MetricCache[column_i].RowRefs =
                (unsigned)Columns[column_i].Rows.size();
        }
        MetricCacheValid = true;
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i) {
            adjust_row_metric_contribution(row_i, 1);
        }
        if (MetricCacheFlags & kMetricMinMax)
        {
            for (uint16_t column_i = 0; column_i < N; ++column_i)
            {
                if (MetricCache[column_i].MinLive == 0xffffu) {
                    MetricCache[column_i].MinLive = 0;
                }
            }
        }
    }

    Metrics cached_metrics(uint16_t column_i, unsigned flags) const
    {
        const unsigned wanted = normalize_metric_flags(flags);
        if (!MetricCacheValid || MetricCache.size() != N ||
            (wanted & ~MetricCacheFlags) != 0)
        {
            const unsigned preserved = MetricCacheValid ? MetricCacheFlags : 0;
            build_metric_cache(preserved | wanted);
        }
        return MetricCache[column_i];
    }

    void increment_weight2_ref(uint16_t column_i, uint16_t skip_row_i = kListTerm)
    {
        if (MetricCacheValid && (MetricCacheFlags & kMetricLookahead))
        {
            for (uint16_t row_i : Columns[column_i].Rows)
            {
                if (row_i == skip_row_i) {
                    continue;
                }
                const Row& row = Rows[row_i];
                if (row.Deferred || row.Live != 2) {
                    continue;
                }

                uint16_t live[2] = {kListTerm, kListTerm};
                const unsigned count = scan_todo_columns(row, live, 2);
                if (count != 2) {
                    continue;
                }
                if (live[0] == column_i) {
                    ++MetricCache[live[1]].Lookahead;
                }
                else if (live[1] == column_i) {
                    ++MetricCache[live[0]].Lookahead;
                }
            }
        }
        ++Columns[column_i].Weight2Refs;
    }

    Key metric_key(int method_id, uint16_t column_i, uint32_t seed) const
    {
        const Column& column = Columns[column_i];
        const unsigned flags = legacy_metric_flags(method_id);
        const Metrics m = flags != 0 ?
            cached_metrics(column_i, flags) : collect_metrics(column_i);

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
        for (uint16_t column_i : TodoColumns)
        {
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

    std::vector<uint16_t> top_default_candidates_from(
        const std::vector<uint16_t>& source,
        unsigned limit,
        uint32_t seed) const
    {
        std::vector<uint16_t> candidates;
        if (limit == 0) {
            return candidates;
        }
        candidates.reserve(limit);

        struct CandidateEntry
        {
            Key Score;
            uint16_t Column;
        };
        std::vector<CandidateEntry> heap;
        heap.reserve(limit);
        const auto better_entry = [](const CandidateEntry& a,
            const CandidateEntry& b) {
            return better_key(a.Score, b.Score);
        };

        for (uint16_t column_i : source)
        {
            if (column_i >= N || Columns[column_i].Mark != kMarkTodo) {
                continue;
            }
            const Key key = default_key(column_i, seed);
            const CandidateEntry entry = {key, column_i};
            if (heap.size() < limit)
            {
                heap.push_back(entry);
                std::push_heap(heap.begin(), heap.end(), better_entry);
            }
            else if (better_key(key, heap.front().Score))
            {
                std::pop_heap(heap.begin(), heap.end(), better_entry);
                heap.back() = entry;
                std::push_heap(heap.begin(), heap.end(), better_entry);
            }
        }

        std::sort(heap.begin(), heap.end(), better_entry);
        for (const CandidateEntry& entry : heap) {
            candidates.push_back(entry.Column);
        }
        return candidates;
    }

    std::vector<uint16_t> top_default_candidates(unsigned limit, uint32_t seed) const
    {
        return top_default_candidates_from(TodoColumns, limit, seed);
    }

    uint16_t select_topk_lookfill(uint32_t seed, unsigned top_k) const
    {
        const std::vector<uint16_t> candidates = top_default_candidates(top_k, seed);
        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : candidates)
        {
            const Metrics m = cached_metrics(column_i,
                kMetricExactD2 | kMetricLookahead | kMetricFill);
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
            if (row.Live <= 1 || row.Deferred) {
                continue;
            }
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

    void add_candidate(std::vector<uint16_t>& candidates,
        std::vector<uint8_t>& seen,
        uint16_t column_i) const
    {
        if (column_i == kListTerm || Columns[column_i].Mark != kMarkTodo) {
            return;
        }
        if (!seen[column_i])
        {
            seen[column_i] = 1;
            candidates.push_back(column_i);
        }
    }

    std::vector<uint16_t> collect_all_candidates() const
    {
        return TodoColumns;
    }

    std::vector<uint16_t> collect_min_row_candidates(bool all_min_rows) const
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

        std::vector<uint16_t> candidates;
        if (best_row == kListTerm) {
            return candidates;
        }

        std::vector<uint8_t> seen(N, 0);
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Live <= 1 || row.Deferred) {
                continue;
            }
            if (row_i != best_row && (!all_min_rows || row.Live != best_live)) {
                continue;
            }
            for (uint16_t column_i : row.Columns) {
                add_candidate(candidates, seen, column_i);
            }
            if (!all_min_rows) {
                break;
            }
        }
        return candidates;
    }

    std::vector<uint16_t> collect_d2_candidates(bool largest_component) const
    {
        std::vector<uint16_t> candidates;
        std::vector<uint8_t> seen(N, 0);

        if (!largest_component)
        {
            for (const Row& row : Rows)
            {
                if (row.Live != 2 || row.Deferred) {
                    continue;
                }
                uint16_t found[2] = {kListTerm, kListTerm};
                const unsigned count = scan_todo_columns(row, found, 2);
                if (count == 2)
                {
                    add_candidate(candidates, seen, found[0]);
                    add_candidate(candidates, seen, found[1]);
                }
            }
            return candidates;
        }

        std::vector<uint16_t> parent(N);
        std::vector<uint16_t> size(N);
        std::vector<uint8_t> active(N, 0);
        std::vector<uint16_t> active_columns;
        const auto activate = [&](uint16_t column_i) {
            if (!active[column_i])
            {
                active[column_i] = 1;
                parent[column_i] = column_i;
                size[column_i] = 1;
                active_columns.push_back(column_i);
            }
        };

        for (const Row& row : Rows)
        {
            if (row.Live != 2 || row.Deferred) {
                continue;
            }
            uint16_t found[2] = {kListTerm, kListTerm};
            const unsigned count = scan_todo_columns(row, found, 2);
            if (count == 2)
            {
                activate(found[0]);
                activate(found[1]);
                union_roots(parent, size, found[0], found[1]);
            }
        }

        uint16_t best_root = kListTerm;
        unsigned best_size = 0;
        for (uint16_t column_i : active_columns)
        {
            const uint16_t root = find_root(parent, column_i);
            if (size[root] > best_size)
            {
                best_size = size[root];
                best_root = root;
            }
        }
        if (best_root == kListTerm) {
            return candidates;
        }

        for (const Row& row : Rows)
        {
            if (row.Live != 2 || row.Deferred) {
                continue;
            }
            uint16_t found[2] = {kListTerm, kListTerm};
            const unsigned count = scan_todo_columns(row, found, 2);
            if (count != 2 || find_root(parent, found[0]) != best_root) {
                continue;
            }
            add_candidate(candidates, seen, found[0]);
            add_candidate(candidates, seen, found[1]);
        }
        return candidates;
    }

    std::vector<uint16_t> collect_sized_d2_component_candidates(bool largest) const
    {
        std::vector<uint16_t> candidates;
        std::vector<uint8_t> seen(N, 0);
        std::vector<uint16_t> parent(N);
        std::vector<uint16_t> size(N);
        std::vector<uint8_t> active(N, 0);
        std::vector<uint16_t> active_columns;
        const auto activate = [&](uint16_t column_i) {
            if (!active[column_i])
            {
                active[column_i] = 1;
                parent[column_i] = column_i;
                size[column_i] = 1;
                active_columns.push_back(column_i);
            }
        };

        for (const Row& row : Rows)
        {
            if (row.Live != 2 || row.Deferred) {
                continue;
            }
            uint16_t found[2] = {kListTerm, kListTerm};
            const unsigned count = scan_todo_columns(row, found, 2);
            if (count == 2)
            {
                activate(found[0]);
                activate(found[1]);
                union_roots(parent, size, found[0], found[1]);
            }
        }

        uint16_t best_root = kListTerm;
        unsigned best_size = largest ? 0u : UINT_MAX;
        for (uint16_t column_i : active_columns)
        {
            const uint16_t root = find_root(parent, column_i);
            const unsigned root_size = size[root];
            const bool take = largest ?
                root_size > best_size :
                root_size < best_size;
            if (best_root == kListTerm || take)
            {
                best_root = root;
                best_size = root_size;
            }
        }
        if (best_root == kListTerm) {
            return candidates;
        }

        for (const Row& row : Rows)
        {
            if (row.Live != 2 || row.Deferred) {
                continue;
            }
            uint16_t found[2] = {kListTerm, kListTerm};
            const unsigned count = scan_todo_columns(row, found, 2);
            if (count != 2 || find_root(parent, found[0]) != best_root) {
                continue;
            }
            add_candidate(candidates, seen, found[0]);
            add_candidate(candidates, seen, found[1]);
        }
        return candidates;
    }

    std::vector<uint16_t> collect_candidates(const Method& method, uint32_t seed) const
    {
        switch (method.Pool)
        {
        case kPoolAll:
            return collect_all_candidates();
        case kPoolTopDefault:
            return top_default_candidates(method.Arg, seed);
        case kPoolD2Rows:
            return collect_d2_candidates(false);
        case kPoolD2LargestComponent:
            return collect_d2_candidates(true);
        case kPoolMinRow:
            return collect_min_row_candidates(false);
        case kPoolAllMinRows:
            return collect_min_row_candidates(true);
        default:
            return collect_all_candidates();
        }
    }

    Key combo_score_key(const Method& method, uint16_t column_i, uint32_t seed) const
    {
        const Column& column = Columns[column_i];
        switch (method.Score)
        {
        case kScoreDefault:
            return default_key(column_i, seed);
        case kScoreExactD2:
            return metric_key(1, column_i, seed);
        case kScoreLookFill:
        {
            const Metrics m = cached_metrics(column_i,
                kMetricExactD2 | kMetricLookahead | kMetricFill);
            return Key((int64_t)m.ExactD2, (int64_t)m.Lookahead,
                -(int64_t)m.Fill, (int64_t)column.Weight2Refs,
                (int64_t)column_i);
        }
        case kScoreRatio:
            return metric_key(12, column_i, seed);
        case kScoreMinFill:
            return metric_key(4, column_i, seed);
        case kScoreMinFillSquare:
            return metric_key(5, column_i, seed);
        case kScoreDuplicate:
            return metric_key(6, column_i, seed);
        case kScoreD3Support:
        {
            const Metrics m = cached_metrics(column_i,
                kMetricExactD2 | kMetricLookahead | kMetricDegree34);
            return Key((int64_t)m.ExactD2, (int64_t)m.Degree3Rows,
                (int64_t)m.Degree4Rows, (int64_t)m.Lookahead,
                (int64_t)column_i);
        }
        case kScoreLowRef:
            return Key(-(int64_t)column.Rows.size(), (int64_t)column.Weight2Refs,
                (int64_t)column_i, 0, 0);
        case kScoreLiveRefMin:
        {
            const Metrics m = cached_metrics(column_i,
                kMetricLiveRefs | kMetricExactD2 | kMetricFill);
            return Key(-(int64_t)m.LiveRefs, (int64_t)m.ExactD2,
                -(int64_t)m.Fill, (int64_t)column.Weight2Refs,
                (int64_t)column_i);
        }
        case kScoreLiveRefMax:
        {
            const Metrics m = cached_metrics(column_i,
                kMetricLiveRefs | kMetricExactD2 | kMetricLookahead |
                kMetricFill);
            return Key((int64_t)m.LiveRefs, (int64_t)m.ExactD2,
                (int64_t)m.Lookahead, -(int64_t)m.Fill,
                (int64_t)column_i);
        }
        case kScoreMaxRowMin:
            return metric_key(13, column_i, seed);
        case kScoreRandom:
        {
            Key key = default_key(column_i, seed);
            key.C = (int64_t)hash32(seed ^ (uint32_t)column_i ^
                ((uint32_t)ChoiceCount * 0x9e3779b9u));
            key.D = column_i;
            return key;
        }
        default:
            return default_key(column_i, seed);
        }
    }

    uint16_t select_combo(const Method& method, uint32_t seed) const
    {
        std::vector<uint16_t> candidates = collect_candidates(method, seed);
        if (candidates.empty()) {
            candidates = collect_all_candidates();
        }

        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : candidates)
        {
            if (Columns[column_i].Mark != kMarkTodo) {
                continue;
            }
            const Key key = combo_score_key(method, column_i, seed);
            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best != kListTerm ? best : select_by_full_scan(0, seed);
    }

    uint16_t select_raptorq_d2_component(uint32_t seed, unsigned tie_mode) const
    {
        std::vector<uint16_t> parent(N);
        std::vector<uint16_t> size(N);
        std::vector<uint8_t> active(N, 0);
        std::vector<uint16_t> active_columns;
        const auto activate = [&](uint16_t column_i) {
            if (!active[column_i])
            {
                active[column_i] = 1;
                parent[column_i] = column_i;
                size[column_i] = 1;
                active_columns.push_back(column_i);
            }
        };

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
                activate(found[0]);
                activate(found[1]);
                union_roots(parent, size, found[0], found[1]);
            }
        }

        uint16_t best_root = kListTerm;
        unsigned best_size = 0;
        for (uint16_t column_i : active_columns)
        {
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
                    const Metrics m = cached_metrics(column_i,
                        kMetricExactD2 | kMetricFill);
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

    uint16_t select_best_scored_candidate(
        const std::vector<uint16_t>& candidates,
        unsigned score_mode,
        uint32_t seed,
        bool use_cache = true) const
    {
        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : candidates)
        {
            if (Columns[column_i].Mark != kMarkTodo) {
                continue;
            }
            unsigned flags = 0;
            switch (score_mode)
            {
            case 1:
                flags = kMetricLiveRefs | kMetricExactD2 | kMetricFill;
                break;
            case 2:
                flags = kMetricLiveRefs | kMetricExactD2 | kMetricLookahead |
                    kMetricFill;
                break;
            case 3:
                flags = kMetricExactD2 | kMetricLookahead | kMetricFill;
                break;
            case 4:
            case 5:
                flags = kMetricExactD2 | kMetricFill;
                break;
            default:
                break;
            }
            const Metrics m = flags != 0 ?
                (use_cache ? cached_metrics(column_i, flags) :
                    collect_metrics(column_i)) : Metrics();
            const unsigned boundary =
                m.LiveRefs > m.ExactD2 ? m.LiveRefs - m.ExactD2 : 0;
            Key key;
            switch (score_mode)
            {
            case 1:
                key = Key(-(int64_t)boundary, (int64_t)m.ExactD2,
                    -(int64_t)m.Fill, (int64_t)Columns[column_i].Weight2Refs,
                    (int64_t)column_i);
                break;
            case 2:
                key = Key((int64_t)boundary, (int64_t)m.ExactD2,
                    (int64_t)m.Lookahead, -(int64_t)m.Fill,
                    (int64_t)column_i);
                break;
            case 3:
                key = Key(-(int64_t)m.Fill, (int64_t)m.ExactD2,
                    (int64_t)m.Lookahead, (int64_t)Columns[column_i].Weight2Refs,
                    (int64_t)column_i);
                break;
            case 4:
                key = Key(-(int64_t)Columns[column_i].Rows.size(),
                    (int64_t)m.ExactD2, -(int64_t)m.Fill,
                    (int64_t)Columns[column_i].Weight2Refs, (int64_t)column_i);
                break;
            case 5:
                key = Key((int64_t)hash32(seed ^ (uint32_t)column_i ^
                    ((uint32_t)ChoiceCount * UINT32_C(0x9e3779b9))),
                    (int64_t)m.ExactD2, -(int64_t)m.Fill,
                    (int64_t)column_i, 0);
                break;
            default:
                key = default_key(column_i, seed);
                break;
            }
            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best;
    }

    std::vector<uint16_t> top_default_subset(
        const std::vector<uint16_t>& candidates,
        unsigned limit,
        uint32_t seed) const
    {
        if (limit == 0) {
            return std::vector<uint16_t>();
        }
        if (candidates.size() <= limit) {
            return candidates;
        }

        struct CandidateEntry
        {
            Key Score;
            uint16_t Column;
        };
        std::vector<CandidateEntry> heap;
        heap.reserve(limit);
        const auto better_entry = [](const CandidateEntry& a,
            const CandidateEntry& b) {
            return better_key(a.Score, b.Score);
        };

        for (uint16_t column_i : candidates)
        {
            if (Columns[column_i].Mark != kMarkTodo) {
                continue;
            }
            const CandidateEntry entry = {
                default_key(column_i, seed),
                column_i
            };
            if (heap.size() < limit)
            {
                heap.push_back(entry);
                std::push_heap(heap.begin(), heap.end(), better_entry);
            }
            else if (better_key(entry.Score, heap.front().Score))
            {
                std::pop_heap(heap.begin(), heap.end(), better_entry);
                heap.back() = entry;
                std::push_heap(heap.begin(), heap.end(), better_entry);
            }
        }

        std::sort(heap.begin(), heap.end(), better_entry);
        std::vector<uint16_t> subset;
        subset.reserve(heap.size());
        for (const CandidateEntry& entry : heap) {
            subset.push_back(entry.Column);
        }
        return subset;
    }

    uint16_t select_ks_boundary_top(
        unsigned limit,
        uint32_t seed) const
    {
        std::vector<uint16_t> candidates =
            collect_sized_d2_component_candidates(true);
        if (candidates.empty()) {
            return select_min_row_degree(seed, true, false);
        }
        candidates = top_default_subset(candidates, limit, seed);
        if (candidates.empty()) {
            return select_min_row_degree(seed, true, false);
        }
        const uint16_t selected =
            select_best_scored_candidate(candidates, 2, seed, false);
        return selected != kListTerm ? selected : select_by_full_scan(0, seed);
    }

    uint16_t select_ks_variant(unsigned mode, uint32_t seed) const
    {
        std::vector<uint16_t> candidates;
        unsigned score_mode = 0;
        switch (mode)
        {
        case 0:
            candidates = collect_d2_candidates(false);
            score_mode = 5;
            break;
        case 1:
            candidates = collect_sized_d2_component_candidates(false);
            score_mode = 0;
            break;
        case 2:
            candidates = collect_sized_d2_component_candidates(true);
            score_mode = 1;
            break;
        case 3:
            candidates = collect_sized_d2_component_candidates(true);
            score_mode = 2;
            break;
        case 4:
            candidates = collect_sized_d2_component_candidates(true);
            score_mode = 3;
            break;
        case 5:
            candidates = collect_sized_d2_component_candidates(true);
            score_mode = 4;
            break;
        case 6:
            candidates = collect_sized_d2_component_candidates(true);
            score_mode = 5;
            break;
        case 7:
            candidates = collect_d2_candidates(false);
            score_mode = 3;
            break;
        default:
            break;
        }
        if (candidates.empty()) {
            return select_min_row_degree(seed, true, false);
        }
        const uint16_t selected =
            select_best_scored_candidate(candidates, score_mode, seed);
        return selected != kListTerm ? selected : select_by_full_scan(0, seed);
    }

    uint16_t select_minrow_scored(unsigned mode, uint32_t seed) const
    {
        uint16_t best_live = 0xffffu;
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Live > 1 && !row.Deferred && row.Live < best_live) {
                best_live = row.Live;
            }
        }
        if (best_live == 0xffffu) {
            return select_by_full_scan(0, seed);
        }

        uint16_t best_column = kListTerm;
        Key best_key;
        std::vector<uint8_t> seen;
        if (mode == 7) {
            seen.assign(N, 0);
        }
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Live != best_live || row.Deferred) {
                continue;
            }

            unsigned row_overlap = 0;
            for (uint16_t c : row.Columns)
            {
                if (Columns[c].Mark == kMarkTodo) {
                    row_overlap += (unsigned)Columns[c].Rows.size();
                }
            }

            for (uint16_t column_i : row.Columns)
            {
                if (Columns[column_i].Mark != kMarkTodo) {
                    continue;
                }
                if (mode == 7)
                {
                    if (seen[column_i]) {
                        continue;
                    }
                    seen[column_i] = 1;
                }

                unsigned flags = 0;
                switch (mode)
                {
                case 0:
                case 1:
                case 7:
                    flags = kMetricExactD2 | kMetricFill;
                    break;
                case 2:
                case 4:
                case 5:
                    flags = kMetricExactD2 | kMetricLookahead | kMetricFill;
                    break;
                case 3:
                    flags = kMetricExactD2 | kMetricLookahead |
                        kMetricDegree34 | kMetricFill;
                    break;
                case 6:
                    flags = kMetricLiveRefs | kMetricExactD2 | kMetricFill;
                    break;
                default:
                    break;
                }
                const Metrics m = flags != 0 ?
                    cached_metrics(column_i, flags) : Metrics();
                Key key;
                switch (mode)
                {
                case 0:
                    key = Key(-(int64_t)row.Columns.size(), (int64_t)m.ExactD2,
                        -(int64_t)m.Fill, (int64_t)column_i, 0);
                    break;
                case 1:
                    key = Key(-(int64_t)row_overlap, (int64_t)m.ExactD2,
                        -(int64_t)m.Fill, (int64_t)column_i, 0);
                    break;
                case 2:
                    key = Key((int64_t)row_overlap, (int64_t)m.ExactD2,
                        (int64_t)m.Lookahead, -(int64_t)m.Fill, (int64_t)column_i);
                    break;
                case 3:
                    key = Key((int64_t)m.ExactD2, (int64_t)m.Lookahead,
                        (int64_t)m.Degree3Rows, -(int64_t)m.Fill, (int64_t)column_i);
                    break;
                case 4:
                    key = Key(-(int64_t)m.Fill, (int64_t)m.ExactD2,
                        (int64_t)m.Lookahead, (int64_t)column_i, 0);
                    break;
                case 5:
                    key = Key((int64_t)((uint64_t)m.ExactD2 * 1000000u /
                        (m.Fill + 1u)), (int64_t)m.ExactD2,
                        (int64_t)m.Lookahead, -(int64_t)m.Fill, (int64_t)column_i);
                    break;
                case 6:
                    key = Key(-(int64_t)m.LiveRefs, (int64_t)m.ExactD2,
                        -(int64_t)m.Fill, (int64_t)column_i, 0);
                    break;
                case 7:
                    key = Key((int64_t)hash32(seed ^ (uint32_t)column_i ^
                        ((uint32_t)ChoiceCount * UINT32_C(0x165667b1))),
                        (int64_t)m.ExactD2, -(int64_t)m.Fill, (int64_t)column_i, 0);
                    break;
                default:
                    key = default_key(column_i, seed);
                    break;
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

    Key beam_terminal_key() const
    {
        unsigned todo = 0;
        uint64_t fill = 0;
        uint64_t d2 = 0;
        for (uint16_t column_i : TodoColumns)
        {
            ++todo;
            const Metrics m = cached_metrics(column_i,
                kMetricExactD2 | kMetricFill);
            fill += m.Fill;
            d2 += m.ExactD2;
        }
        return Key(-(int64_t)todo, (int64_t)d2, -(int64_t)fill,
            -(int64_t)DeferredRows, -(int64_t)DeferCount);
    }

    Key beam_score(unsigned depth, unsigned width, uint32_t seed) const
    {
        if (depth == 0 || count_todo_columns() == 0) {
            return beam_terminal_key();
        }

        std::vector<uint16_t> candidates = top_default_candidates(width, seed);
        if (candidates.empty()) {
            return beam_terminal_key();
        }

        Key best_key;
        bool have_best = false;
        for (uint16_t column_i : candidates)
        {
            PeelSim copy = *this;
            std::vector<uint16_t> action(1, column_i);
            if (!copy.apply_inactivation_action(action)) {
                continue;
            }
            const Key key = copy.beam_score(depth - 1u, width,
                seed ^ hash32((uint32_t)column_i + depth * 0x9e37u));
            if (!have_best || better_key(key, best_key))
            {
                best_key = key;
                have_best = true;
            }
        }
        return have_best ? best_key : beam_terminal_key();
    }

    uint16_t select_beam(unsigned depth, unsigned width, uint32_t seed) const
    {
        if (count_todo_columns() > 8u) {
            return select_topk_lookfill(seed, width);
        }

        std::vector<uint16_t> candidates = top_default_candidates(width, seed);
        if (candidates.empty()) {
            return select_by_full_scan(0, seed);
        }

        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : candidates)
        {
            PeelSim copy = *this;
            std::vector<uint16_t> action(1, column_i);
            if (!copy.apply_inactivation_action(action)) {
                continue;
            }
            const Key key = copy.beam_score(depth > 0 ? depth - 1u : 0u,
                width, seed ^ hash32((uint32_t)column_i));
            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best != kListTerm ? best : select_by_full_scan(0, seed);
    }

    AvalancheEstimate estimate_avalanche(uint16_t column_i) const
    {
        AvalancheEstimate estimate;
        if (column_i >= N || Columns[column_i].Mark != kMarkTodo) {
            return estimate;
        }

        std::vector<uint8_t> removed(N, 0);
        std::vector<uint16_t> queue;
        queue.reserve(64);

        removed[column_i] = 1;
        queue.push_back(column_i);
        estimate.Removed = 1;

        for (unsigned queue_i = 0; queue_i < queue.size(); ++queue_i)
        {
            const uint16_t removed_column = queue[queue_i];
            const std::vector<uint16_t>& refs = Columns[removed_column].Rows;
            for (uint16_t row_i : refs)
            {
                const Row& row = Rows[row_i];
                if (row.Live == 0) {
                    continue;
                }

                unsigned live_after = 0;
                uint16_t only_live = kListTerm;
                for (uint16_t row_column : row.Columns)
                {
                    if (Columns[row_column].Mark != kMarkTodo ||
                        removed[row_column]) {
                        continue;
                    }
                    only_live = row_column;
                    ++live_after;
                    if (live_after > 1) {
                        break;
                    }
                }

                if (live_after == 1 && only_live != kListTerm &&
                    !removed[only_live])
                {
                    removed[only_live] = 1;
                    queue.push_back(only_live);
                    ++estimate.Removed;
                    ++estimate.Peeled;
                    ++estimate.SingletonRows;
                }
            }
        }
        return estimate;
    }

    AvalancheEstimate estimate_avalanche_collect(
        uint16_t column_i,
        std::vector<uint16_t>& dep_rows) const
    {
        AvalancheEstimate estimate;
        if (column_i >= N || Columns[column_i].Mark != kMarkTodo) {
            return estimate;
        }

        const uint32_t epoch = next_avalanche_epoch();
        std::vector<uint32_t>& removed = AvalancheRemovedEpoch;
        std::vector<uint32_t>& row_seen = AvalancheRowEpoch;
        std::vector<uint16_t>& queue = AvalancheQueueScratch;
        queue.clear();

        removed[column_i] = epoch;
        queue.push_back(column_i);
        estimate.Removed = 1;

        for (unsigned queue_i = 0; queue_i < queue.size(); ++queue_i)
        {
            const uint16_t removed_column = queue[queue_i];
            const std::vector<uint16_t>& refs = Columns[removed_column].Rows;
            for (uint16_t row_i : refs)
            {
                const Row& row = Rows[row_i];
                if (row.Live == 0) {
                    continue;
                }
                if (row_seen[row_i] != epoch)
                {
                    row_seen[row_i] = epoch;
                    dep_rows.push_back(row_i);
                }

                unsigned live_after = 0;
                uint16_t only_live = kListTerm;
                for (uint16_t row_column : row.Columns)
                {
                    if (Columns[row_column].Mark != kMarkTodo ||
                        removed[row_column] == epoch) {
                        continue;
                    }
                    only_live = row_column;
                    ++live_after;
                    if (live_after > 1) {
                        break;
                    }
                }

                if (live_after == 1 && only_live != kListTerm &&
                    removed[only_live] != epoch)
                {
                    removed[only_live] = epoch;
                    queue.push_back(only_live);
                    ++estimate.Removed;
                    ++estimate.Peeled;
                    ++estimate.SingletonRows;
                }
            }
        }
        return estimate;
    }

    AvalancheEstimate lazy_avalanche_estimate(uint16_t column_i) const
    {
        if (column_i >= N) {
            return AvalancheEstimate();
        }
        ensure_avalanche_cache();
        if (column_i < AvalancheCacheValid.size() &&
            AvalancheCacheValid[column_i]) {
            return AvalancheCacheValue[column_i];
        }

        AvalancheDepScratch.clear();
        const AvalancheEstimate estimate =
            estimate_avalanche_collect(column_i, AvalancheDepScratch);
        for (uint16_t row_i : AvalancheDepScratch) {
            AvalancheRowDependents[row_i].push_back(column_i);
        }
        AvalancheCacheValue[column_i] = estimate;
        AvalancheCacheValid[column_i] = 1;
        return estimate;
    }

    Key avalanche_key_from_estimate(
        uint16_t column_i,
        const AvalancheEstimate& estimate,
        uint32_t seed) const
    {
        if (estimate.Removed == 0) {
            return Key(INT64_MIN, INT64_MIN, INT64_MIN, INT64_MIN, INT64_MIN);
        }

        const Metrics m = cached_metrics(column_i,
            kMetricExactD2 | kMetricLookahead | kMetricFill);
        const Key tie = default_key(column_i, seed);
        return Key((int64_t)estimate.Removed, (int64_t)estimate.Peeled,
            (int64_t)estimate.SingletonRows + (int64_t)m.ExactD2,
            (int64_t)m.Lookahead - (int64_t)m.Fill, tie.A);
    }

    Key avalanche_estimate_key(uint16_t column_i, uint32_t seed) const
    {
        return avalanche_key_from_estimate(
            column_i, estimate_avalanche(column_i), seed);
    }

    Key exact_one_step_key(uint16_t column_i, uint32_t seed) const
    {
        if (column_i >= N || Columns[column_i].Mark != kMarkTodo) {
            return Key(INT64_MIN, INT64_MIN, INT64_MIN, INT64_MIN, INT64_MIN);
        }

        PeelSim copy = *this;
        std::vector<uint16_t> action(1, column_i);
        if (!copy.apply_inactivation_action(action)) {
            return Key(INT64_MIN, INT64_MIN, INT64_MIN, INT64_MIN, INT64_MIN);
        }

        Key key = copy.beam_terminal_key();
        const Key tie = default_key(column_i, seed);
        key.E = tie.A;
        return key;
    }

    uint16_t select_by_scored_candidates(
        const std::vector<uint16_t>& candidates,
        uint32_t seed,
        bool exact_one_step) const
    {
        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : candidates)
        {
            const Key key = exact_one_step ?
                exact_one_step_key(column_i, seed) :
                avalanche_estimate_key(column_i, seed);
            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best != kListTerm ? best : select_by_full_scan(0, seed);
    }

    uint16_t select_avalanche_top(unsigned top_k, uint32_t seed) const
    {
        const std::vector<uint16_t> candidates = top_default_candidates(top_k, seed);
        return select_by_scored_candidates(candidates, seed, false);
    }

    uint16_t select_avalanche_lazy_top(unsigned top_k, uint32_t seed) const
    {
        const std::vector<uint16_t> candidates =
            top_default_candidates(top_k, seed);
        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : candidates)
        {
            const Key key = avalanche_key_from_estimate(
                column_i, lazy_avalanche_estimate(column_i), seed);
            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best != kListTerm ? best : select_by_full_scan(0, seed);
    }

    uint16_t select_avalanche_rqcc(unsigned top_k, uint32_t seed) const
    {
        std::vector<uint16_t> candidates = collect_d2_candidates(true);
        if (!candidates.empty()) {
            candidates = top_default_candidates_from(candidates, top_k, seed);
        }
        if (candidates.empty()) {
            candidates = top_default_candidates(top_k, seed);
        }
        return select_by_scored_candidates(candidates, seed, false);
    }

    uint16_t select_exact_one_step_top(unsigned top_k, uint32_t seed) const
    {
        const std::vector<uint16_t> candidates = top_default_candidates(top_k, seed);
        return select_by_scored_candidates(candidates, seed, true);
    }

    uint16_t select_exact_one_step_rqcc(unsigned top_k, uint32_t seed) const
    {
        std::vector<uint16_t> candidates = collect_d2_candidates(true);
        if (!candidates.empty()) {
            candidates = top_default_candidates_from(candidates, top_k, seed);
        }
        if (candidates.empty()) {
            candidates = top_default_candidates(top_k, seed);
        }
        return select_by_scored_candidates(candidates, seed, true);
    }

    uint16_t select_hybrid(unsigned mode, uint32_t seed) const
    {
        const unsigned todo = count_todo_columns();
        const unsigned pct = mode == 1 || mode == 3 || mode == 5 ||
            mode == 8 || mode == 10 ? 10u : 5u;
        const unsigned threshold = std::max(1u, (N * pct + 99u) / 100u);
        if (mode == 6) {
            return todo <= threshold ? select_beam(2, 4, seed) :
                select_raptorq_d2_component(seed, 0);
        }
        if (mode == 7 || mode == 8) {
            return todo <= threshold ? select_ks_variant(3, seed) :
                select_raptorq_d2_component(seed, 1);
        }
        if (mode == 9 || mode == 10) {
            return todo <= threshold ? select_avalanche_rqcc(16, seed) :
                select_raptorq_d2_component(seed, 1);
        }
        if (todo > threshold) {
            return select_by_full_scan(0, seed);
        }
        if (mode == 0 || mode == 1) {
            const Method* lowref = find_method(40);
            return lowref ? select_combo(*lowref, seed) :
                select_by_full_scan(0, seed);
        }
        if (mode == 2 || mode == 3) {
            const Method* ratio = find_method(43);
            return ratio ? select_combo(*ratio, seed) :
                select_by_full_scan(0, seed);
        }
        return select_beam(2, 4, seed);
    }

    uint16_t select_lazy_cached(unsigned mode, uint32_t seed) const
    {
        (void)seed;
        uint16_t best = kListTerm;
        Key best_key;
        for (uint16_t column_i : TodoColumns)
        {
            const Column& column = Columns[column_i];
            unsigned flags = kMetricExactD2;
            if (mode == 0) {
                flags |= kMetricLiveRefs | kMetricFill;
            }
            else if (mode == 1) {
                flags |= kMetricFill;
            }
            const Metrics m = cached_metrics(column_i, flags);
            Key key;
            if (mode == 0)
            {
                key = Key(-(int64_t)m.LiveRefs, (int64_t)m.ExactD2,
                    -(int64_t)m.Fill, (int64_t)column.Weight2Refs,
                    (int64_t)column_i);
            }
            else if (mode == 1)
            {
                key = Key(-(int64_t)m.Fill, (int64_t)m.ExactD2,
                    (int64_t)column.Weight2Refs, (int64_t)m.RowRefs,
                    (int64_t)column_i);
            }
            else
            {
                key = Key((int64_t)m.ExactD2,
                    (int64_t)column.Weight2Refs, (int64_t)m.RowRefs,
                    (int64_t)column_i, 0);
            }

            if (best == kListTerm || better_key(key, best_key))
            {
                best = column_i;
                best_key = key;
            }
        }
        return best;
    }

    uint16_t select_column(int method_id, uint32_t seed) const
    {
        const Method* method = find_method(method_id);
        if (method && method->Pool != kPoolLegacy) {
            return select_combo(*method, seed);
        }

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
        case 73:
            return select_ks_variant(0, seed);
        case 74:
            return select_ks_variant(1, seed);
        case 75:
            return select_ks_variant(2, seed);
        case 76:
            return select_ks_variant(3, seed);
        case 77:
            return select_ks_variant(4, seed);
        case 78:
            return select_ks_variant(5, seed);
        case 79:
            return select_ks_variant(6, seed);
        case 80:
            return select_ks_variant(7, seed);
        case 81:
            return select_minrow_scored(0, seed);
        case 82:
            return select_minrow_scored(1, seed);
        case 83:
            return select_minrow_scored(2, seed);
        case 84:
            return select_minrow_scored(3, seed);
        case 85:
            return select_minrow_scored(4, seed);
        case 86:
            return select_minrow_scored(5, seed);
        case 87:
            return select_minrow_scored(6, seed);
        case 88:
            return select_minrow_scored(7, seed);
        case 89:
            return select_beam(2, 4, seed);
        case 90:
            return select_beam(2, 8, seed);
        case 91:
            return select_beam(3, 4, seed);
        case 92:
            return select_beam(3, 8, seed);
        case 93:
            return select_hybrid(0, seed);
        case 94:
            return select_hybrid(1, seed);
        case 95:
            return select_hybrid(2, seed);
        case 96:
            return select_hybrid(3, seed);
        case 97:
            return select_hybrid(4, seed);
        case 98:
            return select_hybrid(5, seed);
        case 99:
            return select_hybrid(6, seed);
        case 112:
            return select_lazy_cached(0, seed);
        case 113:
            return select_lazy_cached(1, seed);
        case 114:
            return select_lazy_cached(2, seed);
        case 115:
            return select_ks_boundary_top(16, seed);
        case 116:
            return select_ks_boundary_top(64, seed);
        case 117:
            return select_hybrid(7, seed);
        case 118:
            return select_hybrid(8, seed);
        case 119:
            return select_avalanche_top(16, seed);
        case 120:
            return select_avalanche_top(64, seed);
        case 121:
            return select_avalanche_rqcc(16, seed);
        case 122:
            return select_avalanche_rqcc(64, seed);
        case 123:
            return select_exact_one_step_top(16, seed);
        case 124:
            return select_exact_one_step_rqcc(16, seed);
        case 125:
            return select_hybrid(9, seed);
        case 126:
            return select_hybrid(10, seed);
        case 127:
            return select_avalanche_lazy_top(16, seed);
        case 128:
            return select_avalanche_lazy_top(64, seed);
        default:
            return select_by_full_scan(method_id, seed);
        }
    }

    std::vector<uint16_t> live_columns(uint16_t row_i) const
    {
        std::vector<uint16_t> live;
        const Row& row = Rows[row_i];
        live.reserve(row.Columns.size());
        for (uint16_t column_i : row.Columns)
        {
            if (Columns[column_i].Mark == kMarkTodo) {
                live.push_back(column_i);
            }
        }
        return live;
    }

    uint16_t select_random_todo_column(uint32_t seed) const
    {
        const unsigned todo = count_todo_columns();
        if (todo == 0) {
            return kListTerm;
        }

        const uint32_t h = hash32(seed ^
            ((uint32_t)ChoiceCount * UINT32_C(0x9e3779b9)) ^
            UINT32_C(0xa5f1523d));
        return TodoColumns[h % todo];
    }

    unsigned choose_leave_index(
        const std::vector<uint16_t>& live,
        uint16_t row_i,
        unsigned leave_mode,
        uint32_t seed) const
    {
        if (live.size() <= 1) {
            return 0;
        }
        if (leave_mode == 0)
        {
            const uint32_t h = hash32(seed ^
                ((uint32_t)row_i * UINT32_C(0x85ebca6b)) ^
                ((uint32_t)ChoiceCount * UINT32_C(0xc2b2ae35)));
            return h % (unsigned)live.size();
        }

        unsigned best_pos = 0;
        Key best_key = default_key(live[0], seed);
        for (unsigned i = 1; i < live.size(); ++i)
        {
            const Key key = default_key(live[i], seed);
            const bool take = leave_mode == 1 ?
                better_key(key, best_key) : worse_key(key, best_key);
            if (take)
            {
                best_pos = i;
                best_key = key;
            }
        }
        return best_pos;
    }

    std::vector<uint16_t> build_row_action(
        uint16_t row_i,
        unsigned leave_mode,
        uint32_t seed) const
    {
        const std::vector<uint16_t> live = live_columns(row_i);
        std::vector<uint16_t> action;
        if (live.size() <= 1) {
            return action;
        }

        const unsigned leave_pos = choose_leave_index(live, row_i, leave_mode, seed);
        action.reserve(live.size() - 1);
        for (unsigned i = 0; i < live.size(); ++i)
        {
            if (i != leave_pos) {
                action.push_back(live[i]);
            }
        }
        return action;
    }

    std::vector<uint16_t> select_min_row_action(
        unsigned leave_mode,
        uint32_t seed) const
    {
        uint16_t best_row = kListTerm;
        Key best_key;

        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred || row.Live <= 1) {
                continue;
            }

            const int64_t tie = (int64_t)hash32(seed ^
                ((uint32_t)row_i * UINT32_C(0x27d4eb2d)) ^
                ((uint32_t)ChoiceCount * UINT32_C(0x165667b1)));
            const Key key(-(int64_t)row.Live, tie, -(int64_t)row_i, 0, 0);
            if (best_row == kListTerm || better_key(key, best_key))
            {
                best_row = row_i;
                best_key = key;
            }
        }

        if (best_row == kListTerm)
        {
            const uint16_t column_i = select_random_todo_column(seed);
            return column_i == kListTerm ?
                std::vector<uint16_t>() : std::vector<uint16_t>(1, column_i);
        }
        return build_row_action(best_row, leave_mode, seed);
    }

    std::vector<uint16_t> select_bounded_min_row_action(
        unsigned max_live,
        unsigned leave_mode,
        uint32_t seed) const
    {
        uint16_t best_row = kListTerm;
        Key best_key;

        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred || row.Live <= 1 || row.Live > max_live) {
                continue;
            }
            const int64_t tie = (int64_t)hash32(seed ^
                ((uint32_t)row_i * UINT32_C(0x27d4eb2d)) ^
                ((uint32_t)ChoiceCount * UINT32_C(0x165667b1)));
            const Key key(-(int64_t)row.Live, tie, -(int64_t)row_i, 0, 0);
            if (best_row == kListTerm || better_key(key, best_key))
            {
                best_row = row_i;
                best_key = key;
            }
        }
        if (best_row == kListTerm) {
            return std::vector<uint16_t>();
        }
        return build_row_action(best_row, leave_mode, seed);
    }

    unsigned bm_row_weight(int method_id, unsigned live) const
    {
        if (method_id == 71) {
            return 1;
        }
        if (method_id == 70)
        {
            if (live <= 2) return 65536u;
            if (live == 3) return 16384u;
            if (live == 4) return 4096u;
            if (live == 5) return 1024u;
            if (live <= 8) return 256u;
            return 64u;
        }

        if (live <= 2) return 65536u;
        if (live == 3) return 64u;
        return 1u;
    }

    std::vector<uint16_t> select_weighted_row_action(
        int method_id,
        uint32_t seed) const
    {
        uint64_t total = 0;
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred) {
                continue;
            }
            const unsigned live = row.Live;
            if (live > 1) {
                total += bm_row_weight(method_id, live);
            }
        }
        if (total == 0) {
            return select_min_row_action(0, seed);
        }

        const uint64_t r =
            ((uint64_t)hash32(seed ^ UINT32_C(0x7f4a7c15) ^
                ((uint32_t)ChoiceCount * UINT32_C(0x9e3779b9))) << 32) |
            (uint64_t)hash32(seed ^ UINT32_C(0x94d049bb) ^
                ((uint32_t)ChoiceCount * UINT32_C(0x85ebca6b)));
        uint64_t target = r % total;

        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred) {
                continue;
            }
            const unsigned live = row.Live;
            if (live <= 1) {
                continue;
            }
            const unsigned weight = bm_row_weight(method_id, live);
            if (target < weight) {
                return build_row_action(row_i, 0, seed);
            }
            target -= weight;
        }

        return select_min_row_action(0, seed);
    }

    std::vector<uint16_t> select_bounded_weighted_row_action(
        int method_id,
        unsigned max_live,
        uint32_t seed) const
    {
        uint64_t total = 0;
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred) {
                continue;
            }
            const unsigned live = row.Live;
            if (live > 1 && live <= max_live) {
                total += bm_row_weight(method_id, live);
            }
        }
        if (total == 0) {
            return std::vector<uint16_t>();
        }

        const uint64_t r =
            ((uint64_t)hash32(seed ^ UINT32_C(0x4112f52b) ^
                ((uint32_t)ChoiceCount * UINT32_C(0x9e3779b9))) << 32) |
            (uint64_t)hash32(seed ^ UINT32_C(0x27d4eb2d) ^
                ((uint32_t)ChoiceCount * UINT32_C(0x85ebca6b)));
        uint64_t target = r % total;

        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred) {
                continue;
            }
            const unsigned live = row.Live;
            if (live <= 1 || live > max_live) {
                continue;
            }
            const unsigned weight = bm_row_weight(method_id, live);
            if (target < weight) {
                return build_row_action(row_i, 0, seed);
            }
            target -= weight;
        }
        return std::vector<uint16_t>();
    }

    Key local_row_action_key(
        uint16_t row_i,
        const std::vector<uint16_t>& live,
        uint32_t seed) const
    {
        uint64_t exact = 0;
        uint64_t lookahead = 0;
        uint64_t weight2 = 0;
        uint64_t fill = 0;
        for (uint16_t column_i : live)
        {
            const Metrics m = cached_metrics(column_i,
                kMetricExactD2 | kMetricLookahead | kMetricFill);
            exact += m.ExactD2;
            lookahead += m.Lookahead;
            weight2 += Columns[column_i].Weight2Refs;
            fill += m.Fill;
        }

        const unsigned residual_cost = (unsigned)live.size() - 1;
        const int64_t score =
            (int64_t)exact * 1024 +
            (int64_t)lookahead +
            (int64_t)weight2 * 8 -
            (int64_t)fill -
            (int64_t)residual_cost * 256;
        const int64_t tie = (int64_t)hash32(seed ^
            ((uint32_t)row_i * UINT32_C(0xd3a2646c)) ^
            ((uint32_t)ChoiceCount * UINT32_C(0xfd7046c5)));
        return Key(score, -(int64_t)residual_cost,
            (int64_t)exact, (int64_t)weight2, tie);
    }

    std::vector<uint16_t> top_local_action_rows(
        unsigned max_live,
        unsigned limit,
        uint32_t seed) const
    {
        struct RowEntry
        {
            Key Score;
            uint16_t Row;
        };

        std::vector<RowEntry> heap;
        heap.reserve(limit);
        const auto better_entry = [](const RowEntry& a, const RowEntry& b) {
            return better_key(a.Score, b.Score);
        };

        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred || row.Live <= 1 || row.Live > max_live) {
                continue;
            }

            const int64_t tie = (int64_t)hash32(seed ^
                ((uint32_t)row_i * UINT32_C(0xd3a2646c)) ^
                ((uint32_t)ChoiceCount * UINT32_C(0xfd7046c5)));
            const RowEntry entry = {
                Key(-(int64_t)row.Live, tie, -(int64_t)row_i, 0, 0),
                row_i
            };

            if (heap.size() < limit)
            {
                heap.push_back(entry);
                std::push_heap(heap.begin(), heap.end(), better_entry);
            }
            else if (better_key(entry.Score, heap.front().Score))
            {
                std::pop_heap(heap.begin(), heap.end(), better_entry);
                heap.back() = entry;
                std::push_heap(heap.begin(), heap.end(), better_entry);
            }
        }

        std::sort(heap.begin(), heap.end(), better_entry);
        std::vector<uint16_t> rows;
        rows.reserve(heap.size());
        for (const RowEntry& entry : heap) {
            rows.push_back(entry.Row);
        }
        return rows;
    }

    std::vector<uint16_t> select_local2_row_action(uint32_t seed) const
    {
        uint16_t best_row = kListTerm;
        Key best_key;

        const std::vector<uint16_t> candidates =
            top_local_action_rows(8, 64, seed);
        for (uint16_t row_i : candidates)
        {
            const std::vector<uint16_t> live = live_columns(row_i);
            if (live.size() <= 1 || live.size() > 8) {
                continue;
            }

            const Key key = local_row_action_key(row_i, live, seed);
            if (best_row == kListTerm || better_key(key, best_key))
            {
                best_row = row_i;
                best_key = key;
            }
        }

        if (best_row == kListTerm) {
            return select_min_row_action(1, seed);
        }
        return build_row_action(best_row, 1, seed);
    }

    std::vector<uint16_t> select_bounded_local2_row_action(
        unsigned max_live,
        uint32_t seed) const
    {
        uint16_t best_row = kListTerm;
        Key best_key;

        const std::vector<uint16_t> candidates =
            top_local_action_rows(max_live, 64, seed);
        for (uint16_t row_i : candidates)
        {
            const std::vector<uint16_t> live = live_columns(row_i);
            if (live.size() <= 1 || live.size() > max_live) {
                continue;
            }

            const Key key = local_row_action_key(row_i, live, seed);
            if (best_row == kListTerm || better_key(key, best_key))
            {
                best_row = row_i;
                best_key = key;
            }
        }

        if (best_row == kListTerm) {
            return std::vector<uint16_t>();
        }
        return build_row_action(best_row, 1, seed);
    }

    std::vector<uint16_t> singleton_action(uint16_t column_i) const
    {
        return column_i == kListTerm ?
            std::vector<uint16_t>() : std::vector<uint16_t>(1, column_i);
    }

    std::vector<uint16_t> select_row_action_hybrid(
        unsigned mode,
        uint32_t seed) const
    {
        std::vector<uint16_t> action;
        switch (mode)
        {
        case 0:
            action = select_bounded_min_row_action(2, 1, seed);
            return action.empty() ? singleton_action(select_by_full_scan(0, seed)) : action;
        case 1:
            action = select_bounded_min_row_action(3, 1, seed);
            return action.empty() ? singleton_action(select_by_full_scan(0, seed)) : action;
        case 2:
            action = select_bounded_min_row_action(4, 1, seed);
            return action.empty() ? singleton_action(select_by_full_scan(0, seed)) : action;
        case 3:
            action = select_bounded_weighted_row_action(69, 3, seed);
            return action.empty() ? singleton_action(select_by_full_scan(0, seed)) : action;
        case 4:
            action = select_bounded_local2_row_action(3, seed);
            if (!action.empty()) {
                return action;
            }
            {
                const Method* lowref = find_method(40);
                return singleton_action(lowref ? select_combo(*lowref, seed) :
                    select_by_full_scan(0, seed));
            }
        case 5:
            action = select_bounded_min_row_action(4, 1, seed);
            if (!action.empty()) {
                return action;
            }
            {
                const Method* lowref = find_method(40);
                return singleton_action(lowref ? select_combo(*lowref, seed) :
                    select_by_full_scan(0, seed));
            }
        default:
            return singleton_action(select_by_full_scan(0, seed));
        }
    }

    std::vector<uint16_t> select_inactivation_action(int method_id, uint32_t seed) const
    {
        switch (method_id)
        {
        case 65:
        {
            const uint16_t column_i = select_random_todo_column(seed);
            return column_i == kListTerm ?
                std::vector<uint16_t>() : std::vector<uint16_t>(1, column_i);
        }
        case 66:
            return select_min_row_action(0, seed);
        case 67:
            return select_min_row_action(1, seed);
        case 68:
            return select_min_row_action(2, seed);
        case 69:
        case 70:
        case 71:
            return select_weighted_row_action(method_id, seed);
        case 72:
            return select_local2_row_action(seed);
        case 106:
            return select_row_action_hybrid(0, seed);
        case 107:
            return select_row_action_hybrid(1, seed);
        case 108:
            return select_row_action_hybrid(2, seed);
        case 109:
            return select_row_action_hybrid(3, seed);
        case 110:
            return select_row_action_hybrid(4, seed);
        case 111:
            return select_row_action_hybrid(5, seed);
        default:
        {
            const uint16_t column_i = select_column(method_id, seed);
            return column_i == kListTerm ?
                std::vector<uint16_t>() : std::vector<uint16_t>(1, column_i);
        }
        }
    }

    bool contains_candidate(const std::vector<uint16_t>& candidates,
        uint16_t column_i) const
    {
        for (uint16_t c : candidates)
        {
            if (c == column_i) {
                return true;
            }
        }
        return false;
    }

    bool has_live_row_action() const
    {
        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (!row.Deferred && live_columns(row_i).size() > 1) {
                return true;
            }
        }
        return false;
    }

    bool action_matches_all_but_one_row(
        const std::vector<uint16_t>& action) const
    {
        if (action.empty()) {
            return false;
        }

        for (uint16_t row_i = 0; row_i < Rows.size(); ++row_i)
        {
            const Row& row = Rows[row_i];
            if (row.Deferred) {
                continue;
            }

            const std::vector<uint16_t> live = live_columns(row_i);
            if (live.size() != action.size() + 1u) {
                continue;
            }

            bool all_found = true;
            for (uint16_t action_col : action)
            {
                bool found = false;
                for (uint16_t live_col : live)
                {
                    if (live_col == action_col)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    all_found = false;
                    break;
                }
            }
            if (all_found) {
                return true;
            }
        }
        return false;
    }

    bool validate_singleton_pool(
        int method_id,
        uint16_t column_i,
        uint32_t seed,
        std::string* why) const
    {
        const Method* method = find_method(method_id);
        if (method && method->Pool != kPoolLegacy)
        {
            const std::vector<uint16_t> candidates = collect_candidates(*method, seed);
            if (!candidates.empty() && !contains_candidate(candidates, column_i)) {
                return fail(why, "method selected outside its candidate pool");
            }
            return true;
        }

        if (method_id == 10 || method_id == 18 || method_id == 19)
        {
            const std::vector<uint16_t> candidates = collect_d2_candidates(true);
            if (!candidates.empty() && !contains_candidate(candidates, column_i)) {
                return fail(why, "degree-2 component method selected outside component pool");
            }
        }
        else if (method_id == 9 || method_id == 15 || method_id == 20)
        {
            const bool all_min_rows = method_id == 20;
            const std::vector<uint16_t> candidates =
                collect_min_row_candidates(all_min_rows);
            if (!candidates.empty() && !contains_candidate(candidates, column_i)) {
                return fail(why, "minimum-row method selected outside row pool");
            }
        }
        else if (method_id == 11 || method_id == 16 || method_id == 17 ||
            method_id == 119 || method_id == 120 || method_id == 123)
        {
            unsigned top_k = 8;
            if (method_id == 16) {
                top_k = 4;
            }
            else if (method_id == 17 || method_id == 119 ||
                method_id == 123) {
                top_k = 16;
            }
            else if (method_id == 120) {
                top_k = 64;
            }
            const std::vector<uint16_t> candidates =
                top_default_candidates(top_k, seed);
            if (!candidates.empty() && !contains_candidate(candidates, column_i)) {
                return fail(why, "top-K method selected outside top-K pool");
            }
        }
        else if (method_id == 121 || method_id == 122 || method_id == 124)
        {
            const unsigned top_k = method_id == 122 ? 64u : 16u;
            std::vector<uint16_t> candidates = collect_d2_candidates(true);
            if (!candidates.empty()) {
                candidates = top_default_candidates_from(candidates, top_k, seed);
            }
            if (!candidates.empty() && !contains_candidate(candidates, column_i)) {
                return fail(why, "avalanche method selected outside component pool");
            }
        }
        return true;
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

static std::string empirical_pdf(std::vector<unsigned> values)
{
    if (values.empty()) {
        return "";
    }

    std::sort(values.begin(), values.end());
    std::string out;
    char buf[64];
    for (size_t i = 0; i < values.size();)
    {
        size_t j = i + 1;
        while (j < values.size() && values[j] == values[i]) {
            ++j;
        }
        const double p = (double)(j - i) / (double)values.size();
        std::snprintf(buf, sizeof(buf), "%u:%.3f", values[i], p);
        if (!out.empty()) {
            out += ';';
        }
        out += buf;
        i = j;
    }
    return out;
}

static unsigned clamp_degree(unsigned degree, unsigned N)
{
    if (degree < 1) {
        degree = 1;
    }
    if (degree > N) {
        degree = N;
    }
    return degree;
}

static unsigned uniform_degree(unsigned min_degree, unsigned max_degree, Rng& rng)
{
    if (max_degree < min_degree) {
        max_degree = min_degree;
    }
    return min_degree + (rng.u32() % (max_degree - min_degree + 1u));
}

static unsigned gaussian_degree(
    unsigned min_degree,
    unsigned max_degree,
    double mean,
    double sigma,
    Rng& rng)
{
    double sum = 0.0;
    for (unsigned i = 0; i < 6; ++i) {
        sum += rng.unit();
    }
    const double normalish = (sum - 3.0) * std::sqrt(2.0);
    int degree = (int)std::floor(mean + sigma * normalish + 0.5);
    if (degree < (int)min_degree) {
        degree = (int)min_degree;
    }
    if (degree > (int)max_degree) {
        degree = (int)max_degree;
    }
    return (unsigned)degree;
}

static double lt_weight(unsigned degree, unsigned N)
{
    if (degree <= 1) {
        return 1.0 / (double)N;
    }
    return 1.0 / ((double)degree * (double)(degree - 1u));
}

static double lt_folded_tail_weight(unsigned cap, unsigned N)
{
    if (cap >= N) {
        return 0.0;
    }
    if (cap <= 1) {
        return 1.0 - 1.0 / (double)N;
    }
    return 1.0 / (double)cap - 1.0 / (double)N;
}

static double robust_soliton_weight(
    unsigned degree,
    unsigned N,
    double c,
    double delta)
{
    if (c <= 0.0) {
        c = 0.01;
    }
    if (delta <= 0.0 || delta >= 1.0) {
        delta = 0.10;
    }

    const double n = (double)N;
    double R = c * std::log(n / delta) * std::sqrt(n);
    if (R < 1.0) {
        R = 1.0;
    }

    unsigned spike = (unsigned)std::floor(n / R);
    if (spike < 1u) {
        spike = 1u;
    }
    if (spike > N) {
        spike = N;
    }

    double tau = 0.0;
    if (degree < spike) {
        tau = R / ((double)degree * n);
    }
    else if (degree == spike) {
        tau = R * std::log(R / delta) / n;
    }
    if (tau < 0.0) {
        tau = 0.0;
    }
    return lt_weight(degree, N) + tau;
}

static double overhead_lt_weight(
    const MatrixStructure& structure,
    unsigned degree,
    unsigned N,
    unsigned row_count)
{
    const unsigned extra = row_count > N ? row_count - N : 0;
    const double overhead = N > 0 ? (double)extra / (double)N : 0.0;

    double degree1_mass = structure.A;
    double degree2_mass = structure.B;
    double degree34_mass = structure.C;

    if (extra == 0)
    {
        degree1_mass *= 0.0;
        degree2_mass *= 0.50;
        degree34_mass *= 0.50;
    }
    else if (extra <= 2)
    {
        degree1_mass *= 0.25;
        degree2_mass *= 0.75;
        degree34_mass *= 0.75;
    }
    else if (overhead < 0.03)
    {
        degree1_mass *= 0.75;
        degree2_mass *= 1.00;
        degree34_mass *= 1.00;
    }
    else
    {
        degree1_mass *= 1.50;
        degree2_mass *= 1.25;
        degree34_mass *= 1.25;
    }

    double w = lt_weight(degree, N);
    if (degree == 1) {
        w += degree1_mass;
    }
    else if (degree == 2) {
        w += degree2_mass;
    }
    else if (degree == 3 || degree == 4) {
        w += degree34_mass;
    }
    return w;
}

// Wrapper kinds (deck-dealt column balancing) reuse the row-degree law of an
// underlying structure kind; everything degree-related dispatches on this so
// the wrapper changes only column selection, never the degree distribution.
static StructureKind structure_degree_law(const MatrixStructure& structure)
{
    if (structure.Kind == kStructureDeckColumns) {
        return structure.DegreeKind;
    }
    return structure.Kind;
}

static double structure_degree_weight(
    const MatrixStructure& structure,
    unsigned degree,
    unsigned N,
    unsigned row_count,
    unsigned max_degree)
{
    switch (structure_degree_law(structure))
    {
    case kStructureLT:
    case kStructureOnePlusLT:
    {
        double w = lt_weight(degree, N);
        if (degree == 1) {
            w += structure.A;
        }
        else if (degree == 2) {
            w += structure.B;
        }
        return w;
    }
    case kStructureHarmonic:
        return 1.0 / (double)degree;
    case kStructureRobustSoliton:
        return robust_soliton_weight(degree, N, structure.A, structure.B);
    case kStructureLTAlpha:
    {
        double alpha = structure.A;
        if (alpha <= 0.0) {
            alpha = 1.0;
        }
        return std::pow(lt_weight(degree, N), alpha);
    }
    case kStructureLTExtra34:
    {
        double w = lt_weight(degree, N);
        if (degree == 2) {
            w += structure.C;
        }
        else if (degree == 3) {
            w += structure.A;
        }
        else if (degree == 4) {
            w += structure.B;
        }
        return w;
    }
    case kStructureLTHighRows:
        return lt_weight(degree, N);
    case kStructureLTFold:
    {
        double w = lt_weight(degree, N);
        if (degree == max_degree) {
            w += lt_folded_tail_weight(max_degree, N);
        }
        return w;
    }
    case kStructureLTFoldScale:
    {
        double w = lt_weight(degree, N);
        if (degree == max_degree) {
            w += structure.A * lt_folded_tail_weight(max_degree, N);
        }
        return w;
    }
    case kStructureOverheadLT:
        return overhead_lt_weight(structure, degree, N, row_count);
    default:
        return 1.0 / (double)degree;
    }
}

class WeightedDegreeSampler
{
public:
    WeightedDegreeSampler() = default;

    WeightedDegreeSampler(
        const MatrixStructure& structure,
        unsigned N,
        unsigned row_count)
    {
        reset(structure, N, row_count);
    }

    void reset(
        const MatrixStructure& structure,
        unsigned N,
        unsigned row_count)
    {
        MinDegree = clamp_degree(structure.MinDegree, N);
        MaxDegree = clamp_degree(structure.MaxDegree, N);
        Total = 0.0;
        Cumulative.clear();
        Cumulative.reserve(MaxDegree - MinDegree + 1u);
        for (unsigned d = MinDegree; d <= MaxDegree; ++d)
        {
            Total += structure_degree_weight(
                structure, d, N, row_count, MaxDegree);
            Cumulative.push_back(Total);
        }
        if (Total <= 0.0) {
            Cumulative.clear();
        }
    }

    unsigned sample(Rng& rng) const
    {
        if (Cumulative.empty()) {
            return uniform_degree(MinDegree, MaxDegree, rng);
        }

        const double target = rng.unit() * Total;
        for (size_t i = 0; i < Cumulative.size(); ++i)
        {
            if (target <= Cumulative[i]) {
                return MinDegree + (unsigned)i;
            }
        }
        return MaxDegree;
    }

private:
    unsigned MinDegree = 1;
    unsigned MaxDegree = 1;
    double Total = 0.0;
    std::vector<double> Cumulative;
};

static bool uses_weighted_degree_sampler(const MatrixStructure& structure)
{
    switch (structure_degree_law(structure))
    {
    case kStructureLT:
    case kStructureHarmonic:
    case kStructureRobustSoliton:
    case kStructureLTAlpha:
    case kStructureLTExtra34:
    case kStructureLTFold:
    case kStructureLTFoldScale:
    case kStructureLTHighRows:
    case kStructureOverheadLT:
    case kStructureOnePlusLT:
        return true;
    default:
        return false;
    }
}

static unsigned choose_structure_degree(
    const MatrixStructure& structure,
    unsigned N,
    const WeightedDegreeSampler& weighted,
    Rng& rng)
{
    const unsigned min_degree = clamp_degree(structure.MinDegree, N);
    const unsigned max_degree = clamp_degree(structure.MaxDegree, N);

    switch (structure_degree_law(structure))
    {
    case kStructureWirehairRandom:
    {
        unsigned degree = GeneratePeelRowWeight(rng.u32(), (uint16_t)N);
        unsigned max_weight = N / 2u;
        if (max_weight < 1u) {
            max_weight = 1u;
        }
        if (degree > max_weight) {
            degree = max_weight;
        }
        return clamp_degree(degree, N);
    }
    case kStructureUniform:
        return uniform_degree(min_degree, max_degree, rng);
    case kStructureGaussian:
        return gaussian_degree(min_degree, max_degree, structure.A, structure.B, rng);
    case kStructureLT:
    case kStructureHarmonic:
    case kStructureRobustSoliton:
    case kStructureLTAlpha:
    case kStructureLTExtra34:
    case kStructureLTFold:
    case kStructureLTFoldScale:
    case kStructureLTHighRows:
    case kStructureOverheadLT:
        return weighted.sample(rng);
    case kStructureOnePlusUniform:
        if (rng.unit() < structure.A) {
            return 1;
        }
        return uniform_degree(min_degree, max_degree, rng);
    case kStructureOnePlusLT:
        if (rng.unit() < structure.A) {
            return 1;
        }
        return weighted.sample(rng);
    case kStructureTwoModeUniform:
        if (rng.unit() < structure.A)
        {
            return uniform_degree(
                clamp_degree(structure.AltMinDegree, N),
                clamp_degree(structure.AltMaxDegree, N),
                rng);
        }
        return uniform_degree(min_degree, max_degree, rng);
    case kStructureBalancedUniform:
        return uniform_degree(min_degree, max_degree, rng);
    default:
        return 1;
    }
}

static std::vector<uint16_t> random_row_columns(
    unsigned N,
    unsigned degree,
    Rng& rng)
{
    degree = clamp_degree(degree, N);
    std::vector<uint16_t> row;
    row.reserve(degree);
    if (degree == N)
    {
        for (unsigned column_i = 0; column_i < N; ++column_i) {
            row.push_back((uint16_t)column_i);
        }
        return row;
    }
    if (degree * 4u >= N)
    {
        std::vector<uint16_t> deck;
        deck.reserve(N);
        for (unsigned column_i = 0; column_i < N; ++column_i) {
            deck.push_back((uint16_t)column_i);
        }
        for (unsigned i = 0; i < degree; ++i)
        {
            const unsigned j = i + (rng.u32() % (N - i));
            const uint16_t column_i = deck[j];
            deck[j] = deck[i];
            deck[i] = column_i;
            row.push_back(column_i);
        }
        return row;
    }
    while (row.size() < degree)
    {
        const uint16_t column_i = (uint16_t)(rng.u32() % N);
        bool duplicate = false;
        for (uint16_t c : row)
        {
            if (c == column_i)
            {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            row.push_back(column_i);
        }
    }
    return row;
}

static std::vector<uint16_t> binary_p50_row_columns(unsigned N, Rng& rng)
{
    std::vector<uint16_t> row;
    row.reserve((N + 1u) / 2u);
    for (unsigned column_i = 0; column_i < N; ++column_i)
    {
        if (rng.u32() & 1u) {
            row.push_back((uint16_t)column_i);
        }
    }
    if (row.empty()) {
        row.push_back((uint16_t)(rng.u32() % N));
    }
    return row;
}

static unsigned row_projected_max_load(
    const std::vector<uint16_t>& row,
    const std::vector<unsigned>& column_loads)
{
    unsigned max_load = 0;
    for (uint16_t column_i : row)
    {
        const unsigned load = column_loads[column_i] + 1u;
        if (load > max_load) {
            max_load = load;
        }
    }
    return max_load;
}

static unsigned row_projected_sum_load(
    const std::vector<uint16_t>& row,
    const std::vector<unsigned>& column_loads)
{
    unsigned sum = 0;
    for (uint16_t column_i : row) {
        sum += column_loads[column_i] + 1u;
    }
    return sum;
}

static bool row_within_load_limit(
    const std::vector<uint16_t>& row,
    const std::vector<unsigned>& column_loads,
    unsigned total_refs,
    unsigned N,
    double limit)
{
    const double mean = (double)(total_refs + row.size()) / (double)N;
    unsigned max_allowed = (unsigned)std::ceil(mean * limit);
    if (max_allowed < 1u) {
        max_allowed = 1u;
    }

    for (uint16_t column_i : row)
    {
        if (column_loads[column_i] + 1u > max_allowed) {
            return false;
        }
    }
    return true;
}

static uint16_t sample_column_not_in_row(
    unsigned N,
    const std::vector<uint16_t>& row,
    Rng& rng)
{
    for (;;)
    {
        const uint16_t column_i = (uint16_t)(rng.u32() % N);
        bool duplicate = false;
        for (uint16_t c : row)
        {
            if (c == column_i)
            {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            return column_i;
        }
    }
}

static std::vector<uint16_t> balanced_row_columns(
    const MatrixStructure& structure,
    unsigned N,
    unsigned degree,
    const std::vector<unsigned>& column_loads,
    unsigned total_refs,
    Rng& rng)
{
    if (structure.BalanceMode == 1)
    {
        std::vector<uint16_t> fallback;
        for (unsigned attempt = 0; attempt < 16; ++attempt)
        {
            std::vector<uint16_t> row = random_row_columns(N, degree, rng);
            if (fallback.empty()) {
                fallback = row;
            }
            if (row_within_load_limit(row, column_loads, total_refs, N,
                structure.BalanceValue)) {
                return row;
            }
        }
        return fallback;
    }

    if (structure.BalanceMode == 2)
    {
        unsigned trials = (unsigned)(structure.BalanceValue + 0.5);
        if (trials < 1u) {
            trials = 1u;
        }

        std::vector<uint16_t> best;
        unsigned best_max = UINT_MAX;
        unsigned best_sum = UINT_MAX;
        for (unsigned trial = 0; trial < trials; ++trial)
        {
            std::vector<uint16_t> row = random_row_columns(N, degree, rng);
            const unsigned projected_max =
                row_projected_max_load(row, column_loads);
            const unsigned projected_sum =
                row_projected_sum_load(row, column_loads);
            if (best.empty() || projected_max < best_max ||
                (projected_max == best_max && projected_sum < best_sum))
            {
                best = row;
                best_max = projected_max;
                best_sum = projected_sum;
            }
        }
        return best;
    }

    if (structure.BalanceMode == 3)
    {
        // Power-of-two-choices: for each edge sample two iid candidate
        // columns (rejecting duplicates within the row) and keep the one
        // with the lower current column load; ties keep the first.
        degree = clamp_degree(degree, N);
        if (degree >= N) {
            return random_row_columns(N, N, rng);
        }
        std::vector<uint16_t> row;
        row.reserve(degree);
        while (row.size() < degree)
        {
            const uint16_t first = sample_column_not_in_row(N, row, rng);
            const uint16_t second = sample_column_not_in_row(N, row, rng);
            row.push_back(
                column_loads[second] < column_loads[first] ? second : first);
        }
        return row;
    }

    return random_row_columns(N, degree, rng);
}

static bool is_prime(unsigned value)
{
    if (value < 2u) {
        return false;
    }
    if (value == 2u) {
        return true;
    }
    if ((value & 1u) == 0) {
        return false;
    }
    for (unsigned d = 3; d <= value / d; d += 2)
    {
        if (value % d == 0) {
            return false;
        }
    }
    return true;
}

static unsigned next_prime(unsigned value)
{
    if (value <= 2u) {
        return 2u;
    }
    if ((value & 1u) == 0) {
        ++value;
    }
    while (!is_prime(value)) {
        value += 2u;
    }
    return value;
}

static bool is_raptorq_ldpc_structure(const MatrixStructure& structure)
{
    return structure.Kind == kStructureRaptorQLdpcStruct ||
        structure.Kind == kStructureRaptorQLdpcRandom;
}

struct RaptorQLdpcParams
{
    unsigned S;
    unsigned P;
};

static RaptorQLdpcParams raptorq_ldpc_params(unsigned N)
{
    // RFC 6330 Table 2 parameters for the N-jitter bands used here.
    // P is L-W = K' + S + H - W for the selected K' table entry.
    if (N >= 310u && N <= 324u) {
        return {31u, 28u};
    }
    if (N >= 325u && N <= 337u) {
        return {31u, 29u};
    }
    if (N >= 3190u && N <= 3210u) {
        return {113u, 89u};
    }
    if (N >= 31990u && N <= 32010u) {
        return {577u, 284u};
    }

    if (N < 6u) {
        return {1u, 1u};
    }

    unsigned target = (unsigned)std::ceil((double)N * 0.015) + 3u;
    if (target < 5u) {
        target = 5u;
    }
    unsigned s = next_prime(target);

    const unsigned max_s = (N - 1u) / 3u;
    if (s > max_s)
    {
        s = max_s;
        while (s > 1u && !is_prime(s)) {
            --s;
        }
        if (s < 1u) {
            s = 1u;
        }
    }
    return {s, s};
}

static unsigned raptorq_ldpc_symbol_count(unsigned N)
{
    return raptorq_ldpc_params(N).S;
}

static unsigned generated_row_count_bound(
    const MatrixStructure& structure,
    unsigned N,
    unsigned row_count)
{
    if (!is_raptorq_ldpc_structure(structure)) {
        return row_count;
    }
    return row_count + raptorq_ldpc_symbol_count(N);
}

static void add_unique_column(std::vector<uint16_t>& row, unsigned column)
{
    const uint16_t c = (uint16_t)column;
    for (uint16_t existing : row)
    {
        if (existing == c) {
            return;
        }
    }
    row.push_back(c);
}

static void add_raptorq_ldpc_rows(
    const MatrixStructure& structure,
    unsigned N,
    uint64_t seed,
    std::vector<std::vector<uint16_t> >& rows)
{
    Rng rng(seed ^ UINT64_C(0x726170746f72716c));
    const RaptorQLdpcParams params = raptorq_ldpc_params(N);
    const unsigned S = params.S;
    const unsigned P = params.P;
    if (S == 0 || N <= S + P) {
        return;
    }

    const unsigned B = N - S - P;
    const unsigned ldpc_base = B;
    const unsigned pi_base = B + S;
    std::vector<std::vector<uint16_t> > checks(S);

    for (unsigned i = 0; i < S; ++i)
    {
        add_unique_column(checks[i], ldpc_base + i);
        add_unique_column(checks[i], pi_base + (i % P));
        add_unique_column(checks[i], pi_base + ((i + 1u) % P));
    }

    for (unsigned i = 0; i < B; ++i)
    {
        const unsigned a = 1u + i / S;
        unsigned b = i % S;
        add_unique_column(checks[b], i);
        b = (b + a) % S;
        add_unique_column(checks[b], i);
        b = (b + a) % S;
        add_unique_column(checks[b], i);
    }

    for (const std::vector<uint16_t>& check : checks)
    {
        if (structure.Kind == kStructureRaptorQLdpcRandom) {
            rows.push_back(random_row_columns(N, (unsigned)check.size(), rng));
        }
        else {
            rows.push_back(check);
        }
    }
}

// Appends one freshly shuffled permutation of all N column ids to the deck.
static void append_shuffled_deck(
    std::vector<uint16_t>& deck,
    unsigned N,
    Rng& rng)
{
    const size_t base = deck.size();
    for (unsigned column_i = 0; column_i < N; ++column_i) {
        deck.push_back((uint16_t)column_i);
    }
    for (unsigned i = 0; i < N; ++i)
    {
        const unsigned j = i + (rng.u32() % (N - i));
        const uint16_t swapped = deck[base + j];
        deck[base + j] = deck[base + i];
        deck[base + i] = swapped;
    }
}

static unsigned deck_virtual_row_count(
    const MatrixStructure& structure,
    unsigned row_count)
{
    double multiplier = structure.BalanceValue;
    if (multiplier < 1.0) {
        multiplier = 1.0;
    }
    const double virtual_rows = std::ceil((double)row_count * multiplier);
    if (virtual_rows <= (double)row_count) {
        return row_count;
    }
    return (unsigned)virtual_rows;
}

// Rare end-of-deal conflicts (the final cards all duplicating columns of the
// last rows) can leave one column a use short and another a use over.  Move
// single references from globally max-loaded columns to min-loaded columns
// until usage counts differ by at most 1.  Each move strictly reduces the
// usage-count spread potential, so this terminates quickly.
static void rebalance_deck_columns(
    std::vector<std::vector<uint16_t> >& rows,
    unsigned N)
{
    std::vector<unsigned> counts(N, 0);
    for (const std::vector<uint16_t>& row : rows)
    {
        for (uint16_t column_i : row) {
            ++counts[column_i];
        }
    }

    for (;;)
    {
        unsigned min_col = 0;
        unsigned max_col = 0;
        for (unsigned column_i = 1; column_i < N; ++column_i)
        {
            if (counts[column_i] < counts[min_col]) {
                min_col = column_i;
            }
            if (counts[column_i] > counts[max_col]) {
                max_col = column_i;
            }
        }
        if (counts[max_col] - counts[min_col] <= 1u) {
            break;
        }

        bool moved = false;
        for (std::vector<uint16_t>& row : rows)
        {
            bool has_min = false;
            bool has_max = false;
            size_t max_pos = 0;
            for (size_t i = 0; i < row.size(); ++i)
            {
                if (row[i] == (uint16_t)min_col) {
                    has_min = true;
                }
                else if (row[i] == (uint16_t)max_col)
                {
                    has_max = true;
                    max_pos = i;
                }
            }
            if (has_max && !has_min)
            {
                row[max_pos] = (uint16_t)min_col;
                --counts[max_col];
                ++counts[min_col];
                moved = true;
                break;
            }
        }
        if (!moved) {
            break;
        }
    }
}

// Deck-dealt column-degree regularization: row degrees come from the
// underlying degree law exactly as for the iid structure, but columns are
// dealt from repeated shuffled decks of all N column ids so column degrees
// are exactly balanced (max-min <= 1) at the same total edge budget.  Cards
// that would duplicate a column already in the current row are skipped and
// stay queued for later rows.  Loss variants (BalanceValue > 1) deal across
// extra virtual rows sampled from the same law and keep only the configured
// row count, modeling balance surviving random row loss.
static std::vector<std::vector<uint16_t> > generate_deck_structure_rows(
    const MatrixStructure& structure,
    unsigned N,
    unsigned row_count,
    uint64_t seed)
{
    Rng rng(seed);
    WeightedDegreeSampler degree_sampler;
    if (uses_weighted_degree_sampler(structure)) {
        degree_sampler.reset(structure, N, row_count);
    }

    const unsigned virtual_rows = deck_virtual_row_count(structure, row_count);
    std::vector<unsigned> degrees(virtual_rows, 1);
    for (unsigned row_i = 0; row_i < virtual_rows; ++row_i)
    {
        degrees[row_i] = clamp_degree(
            choose_structure_degree(structure, N, degree_sampler, rng), N);
    }

    std::vector<uint16_t> deck;
    size_t deck_pos = 0;
    std::vector<uint8_t> in_row(N, 0);
    std::vector<std::vector<uint16_t> > rows;
    rows.reserve(virtual_rows);
    for (unsigned row_i = 0; row_i < virtual_rows; ++row_i)
    {
        const unsigned degree = degrees[row_i];
        std::vector<uint16_t> row;
        row.reserve(degree);
        while (row.size() < degree)
        {
            size_t scan = deck_pos;
            for (;;)
            {
                if (scan >= deck.size()) {
                    append_shuffled_deck(deck, N, rng);
                }
                if (!in_row[deck[scan]]) {
                    break;
                }
                // Skip-and-requeue: this card duplicates a column already
                // in the row, so leave it queued and look at the next card.
                ++scan;
            }
            const uint16_t column_i = deck[scan];
            deck[scan] = deck[deck_pos];
            ++deck_pos;
            in_row[column_i] = 1;
            row.push_back(column_i);
        }
        for (uint16_t column_i : row) {
            in_row[column_i] = 0;
        }
        rows.push_back(row);
    }

    rebalance_deck_columns(rows, N);

    if (rows.size() > row_count) {
        rows.resize(row_count);
    }
    return rows;
}

static std::vector<std::vector<uint16_t> > generate_random_structure_rows(
    const MatrixStructure& structure,
    unsigned N,
    unsigned row_count,
    uint64_t seed)
{
    if (structure.Kind == kStructureDeckColumns) {
        return generate_deck_structure_rows(structure, N, row_count, seed);
    }

    Rng rng(seed);
    WeightedDegreeSampler degree_sampler;
    if (uses_weighted_degree_sampler(structure)) {
        degree_sampler.reset(structure, N, row_count);
    }
    std::vector<unsigned> column_loads(N, 0);
    unsigned total_refs = 0;
    std::vector<std::vector<uint16_t> > rows;
    rows.reserve(generated_row_count_bound(structure, N, row_count));
    if (is_raptorq_ldpc_structure(structure)) {
        add_raptorq_ldpc_rows(structure, N, seed, rows);
    }
    for (unsigned row_i = 0; row_i < row_count; ++row_i)
    {
        std::vector<uint16_t> row;
        if (structure.Kind == kStructureBinaryP50)
        {
            row = binary_p50_row_columns(N, rng);
        }
        else if (structure.Kind == kStructureLTHighRows &&
            row_i < (unsigned)(structure.A + 0.5))
        {
            unsigned degree = (unsigned)std::floor(
                structure.B * (double)N + 0.5);
            row = random_row_columns(N, clamp_degree(degree, N), rng);
        }
        else
        {
            const unsigned degree = is_raptorq_ldpc_structure(structure) ?
                choose_structure_degree(
                    kStructures[0], N, degree_sampler, rng) :
                choose_structure_degree(
                    structure, N, degree_sampler, rng);
            row = (structure.Kind == kStructureBalancedUniform ||
                structure.BalanceMode != 0) ?
                balanced_row_columns(structure, N, degree, column_loads,
                    total_refs, rng) :
                random_row_columns(N, degree, rng);
        }
        for (uint16_t column_i : row) {
            ++column_loads[column_i];
        }
        total_refs += (unsigned)row.size();
        rows.push_back(row);
    }
    return rows;
}

static MatrixWork measure_matrix_work(
    const std::vector<std::vector<uint16_t> >& rows)
{
    MatrixWork work;
    for (const std::vector<uint16_t>& row : rows)
    {
        work.Refs += (uint64_t)row.size();
        if (!row.empty()) {
            work.Xors += (uint64_t)(row.size() - 1u);
        }
    }
    return work;
}

static double dense_solve_xors_estimate(
    unsigned residual_cols,
    unsigned residual_rows)
{
    if (residual_cols == 0 || residual_rows == 0) {
        return 0.0;
    }
    const unsigned active_rows =
        residual_rows > residual_cols ? residual_rows : residual_cols;
    return (double)residual_cols * (double)active_rows;
}

static bool feed_row_vectors(
    PeelSim& sim,
    const std::vector<std::vector<uint16_t> >& rows,
    std::string* why)
{
    for (const std::vector<uint16_t>& row : rows)
    {
        if (!sim.feed_columns(row, why)) {
            return false;
        }
    }
    return true;
}

static bool feed_random_structure(
    PeelSim& sim,
    const MatrixStructure& structure,
    unsigned N,
    unsigned row_count,
    uint64_t seed,
    std::string* why)
{
    const std::vector<std::vector<uint16_t> > rows =
        generate_random_structure_rows(structure, N, row_count, seed);
    return feed_row_vectors(sim, rows, why);
}

static unsigned sample_actual_n(
    unsigned anchor_n,
    unsigned jitter,
    int trial,
    uint64_t seed)
{
    unsigned low = anchor_n;
    if (jitter > low || low - jitter < CAT_WIREHAIR_MIN_N) {
        low = CAT_WIREHAIR_MIN_N;
    }
    else {
        low -= jitter;
    }

    unsigned high = anchor_n + jitter;
    if (high < anchor_n || high > CAT_WIREHAIR_MAX_N) {
        high = CAT_WIREHAIR_MAX_N;
    }

    if (low > high) {
        low = high = anchor_n;
    }
    if (low == high) {
        return low;
    }

    const uint64_t jitter_seed = seed ^
        ((uint64_t)anchor_n * UINT64_C(0xbf58476d1ce4e5b9)) ^
        ((uint64_t)trial * UINT64_C(0x94d049bb133111eb)) ^
        UINT64_C(0x6a09e667f3bcc909);
    Rng rng(jitter_seed);
    return low + (unsigned)(rng.u32() % (high - low + 1u));
}

static TrialConfig make_trial_config(
    int method_id,
    unsigned anchor_n,
    int trial,
    int overhead,
    int rows_override,
    double overhead_pct,
    unsigned n_jitter,
    uint64_t seed,
    const MatrixStructure& structure,
    const std::string& matrix_seed_mode)
{
    TrialConfig config;
    config.ActualN = sample_actual_n(anchor_n, n_jitter, trial, seed);
    const double pct_rows = std::ceil(
        (double)config.ActualN * overhead_pct / 100.0);
    const double row_count = rows_override > 0 ?
        (double)rows_override :
        (double)config.ActualN + (double)overhead + pct_rows;
    if (!std::isfinite(row_count) || row_count > 65535.0) {
        config.MatrixRows = 65536u;
    }
    else {
        config.MatrixRows = row_count > 0.0 ? (unsigned)row_count : 0;
    }

    config.BaseSeed = seed ^
        ((uint64_t)config.ActualN * UINT64_C(0xbf58476d1ce4e5b9)) ^
        ((uint64_t)config.MatrixRows * UINT64_C(0xd6e8feb86659fd93)) ^
        ((uint64_t)trial * UINT64_C(0x94d049bb133111eb)) ^
        (hash_string64(structure.Name) * UINT64_C(0x9e3779b97f4a7c15));
    config.MatrixSeed = matrix_seed_mode == "paired" ?
        config.BaseSeed :
        (config.BaseSeed ^
         ((uint64_t)method_id * UINT64_C(0x632be59bd9b4e019)));
    config.ActionSeed = hash32(
        (uint32_t)config.BaseSeed ^
        ((uint32_t)method_id * UINT32_C(0x85ebca6b)) ^
        ((uint32_t)trial * UINT32_C(0xc2b2ae35)));
    return config;
}

static bool self_fail(std::string* why, const std::string& message)
{
    if (why) {
        *why = message;
    }
    return false;
}

static bool validate_raw_structure(
    const char* name,
    unsigned N,
    const std::vector<std::vector<uint16_t> >& rows,
    unsigned min_degree,
    unsigned max_degree,
    unsigned min_one_rows,
    std::string* why)
{
    if (rows.empty()) {
        return self_fail(why, std::string(name) + ": empty structure");
    }

    unsigned one_rows = 0;
    for (unsigned row_i = 0; row_i < rows.size(); ++row_i)
    {
        const std::vector<uint16_t>& row = rows[row_i];
        if (row.size() < min_degree || row.size() > max_degree)
        {
            return self_fail(why, std::string(name) +
                ": row degree outside expected bounds");
        }
        if (row.size() == 1) {
            ++one_rows;
        }

        std::vector<uint8_t> seen(N, 0);
        for (uint16_t column_i : row)
        {
            if (column_i >= N) {
                return self_fail(why, std::string(name) +
                    ": row column out of range");
            }
            if (seen[column_i]) {
                return self_fail(why, std::string(name) +
                    ": duplicate column inside row");
            }
            seen[column_i] = 1;
        }
    }

    if (one_rows < min_one_rows) {
        return self_fail(why, std::string(name) +
            ": not enough degree-1 rows");
    }
    return true;
}

static bool run_structure_methods(
    const char* name,
    unsigned N,
    const std::vector<std::vector<uint16_t> >& rows,
    unsigned min_degree,
    unsigned max_degree,
    unsigned min_one_rows,
    std::string* why)
{
    if (!validate_raw_structure(name, N, rows, min_degree, max_degree,
        min_one_rows, why)) {
        return false;
    }

    for (const Method& method : kMethods)
    {
        PeelSim sim(N);
        for (const std::vector<uint16_t>& row : rows)
        {
            std::string detail;
            if (!sim.add_test_row(row, &detail))
            {
                return self_fail(why, std::string(name) + " add row failed for " +
                    method.Name + ": " + detail);
            }
        }

        const uint32_t seed = hash32((uint32_t)method.Id ^
            (uint32_t)(N * 0x9e37u) ^ UINT32_C(0x51f15eed));
        TrialResult result;
        std::string detail;
        if (!sim.run_greedy_checked(method.Id, seed, &result, &detail))
        {
            return self_fail(why, std::string(name) + " checked run failed for " +
                method.Name + ": " + detail);
        }
    }
    return true;
}

static bool run_synthetic_structure_tests(std::string* why)
{
    {
        PeelSim sim(4);
        std::string detail;
        if (sim.add_test_row(std::vector<uint16_t>{0, 0}, &detail)) {
            return self_fail(why, "duplicate-column test row was accepted");
        }
        if (!sim.validate_state(false, &detail)) {
            return self_fail(why, "duplicate-column rejection corrupted state: " + detail);
        }
        if (sim.add_test_row(std::vector<uint16_t>{1, 4}, &detail)) {
            return self_fail(why, "out-of-range test row was accepted");
        }
        if (!sim.validate_state(false, &detail)) {
            return self_fail(why, "out-of-range rejection corrupted state: " + detail);
        }
    }

    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0});
        rows.push_back(std::vector<uint16_t>{1, 2});
        rows.push_back(std::vector<uint16_t>{2, 3, 4});
        rows.push_back(std::vector<uint16_t>{4, 5, 6, 7});
        rows.push_back(std::vector<uint16_t>{7, 8});
        if (!run_structure_methods("min1", 12, rows, 1, 4, 1, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0, 1});
        rows.push_back(std::vector<uint16_t>{1, 2, 3});
        rows.push_back(std::vector<uint16_t>{3, 4});
        rows.push_back(std::vector<uint16_t>{4, 5, 6, 7, 8});
        rows.push_back(std::vector<uint16_t>{8, 9, 10});
        if (!run_structure_methods("min2", 12, rows, 2, 5, 0, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0, 1, 2});
        rows.push_back(std::vector<uint16_t>{2, 3, 4, 5});
        rows.push_back(std::vector<uint16_t>{5, 6, 7});
        rows.push_back(std::vector<uint16_t>{7, 8, 9, 10, 11});
        rows.push_back(std::vector<uint16_t>{1, 6, 10});
        if (!run_structure_methods("min3", 12, rows, 3, 5, 0, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0, 1});
        rows.push_back(std::vector<uint16_t>{1, 2});
        rows.push_back(std::vector<uint16_t>{2, 3});
        rows.push_back(std::vector<uint16_t>{4, 5});
        rows.push_back(std::vector<uint16_t>{5, 6});
        rows.push_back(std::vector<uint16_t>{6, 4});
        rows.push_back(std::vector<uint16_t>{7, 8});
        if (!run_structure_methods("d2_components", 10, rows, 2, 2, 0, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0});
        rows.push_back(std::vector<uint16_t>{1, 2});
        rows.push_back(std::vector<uint16_t>{2, 3});
        rows.push_back(std::vector<uint16_t>{3, 4, 5});
        rows.push_back(std::vector<uint16_t>{5, 6, 7, 8});
        rows.push_back(std::vector<uint16_t>{0, 2, 4, 6, 8, 10, 12, 14});
        if (!run_structure_methods("lt_like", 16, rows, 1, 8, 1, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0});
        rows.push_back(std::vector<uint16_t>{1, 2});
        rows.push_back(std::vector<uint16_t>{2, 3, 4});
        rows.push_back(std::vector<uint16_t>{4, 5, 6, 7});
        rows.push_back(std::vector<uint16_t>{7, 8, 9, 10, 11});
        rows.push_back(std::vector<uint16_t>{0, 3, 6, 9, 12, 15});
        if (!run_structure_methods("uniform_1_6", 16, rows, 1, 6, 1, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0, 1, 2});
        rows.push_back(std::vector<uint16_t>{2, 3, 4, 5});
        rows.push_back(std::vector<uint16_t>{5, 6, 7, 8, 9});
        rows.push_back(std::vector<uint16_t>{1, 4, 7, 10});
        rows.push_back(std::vector<uint16_t>{0, 6, 11});
        rows.push_back(std::vector<uint16_t>{3, 8, 12, 13, 14});
        if (!run_structure_methods("gaussian4", 16, rows, 3, 5, 0, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0});
        rows.push_back(std::vector<uint16_t>{3});
        rows.push_back(std::vector<uint16_t>{6});
        rows.push_back(std::vector<uint16_t>{1, 2, 3});
        rows.push_back(std::vector<uint16_t>{4, 5, 6});
        rows.push_back(std::vector<uint16_t>{7, 8, 9, 10});
        if (!run_structure_methods("one_fraction", 12, rows, 1, 4, 3, why)) {
            return false;
        }
    }
    {
        std::vector<std::vector<uint16_t> > rows;
        rows.push_back(std::vector<uint16_t>{0, 1});
        rows.push_back(std::vector<uint16_t>{1, 2, 3});
        rows.push_back(std::vector<uint16_t>{3, 4, 5, 6});
        rows.push_back(std::vector<uint16_t>{6, 7});
        rows.push_back(std::vector<uint16_t>{8, 9, 10, 11});
        rows.push_back(std::vector<uint16_t>{0, 5, 10});
        if (!run_structure_methods("cap4", 12, rows, 2, 4, 0, why)) {
            return false;
        }
    }
    return true;
}

static bool run_matrix_work_tests(std::string* why)
{
    std::vector<std::vector<uint16_t> > rows;
    rows.push_back(std::vector<uint16_t>{0, 1, 2});
    rows.push_back(std::vector<uint16_t>{3});
    rows.push_back(std::vector<uint16_t>{1, 4});

    const MatrixWork work = measure_matrix_work(rows);
    if (work.Refs != 6u || work.Xors != 3u) {
        return self_fail(why, "matrix row XOR accounting is wrong");
    }
    if (dense_solve_xors_estimate(0, 5) != 0.0 ||
        dense_solve_xors_estimate(4, 6) != 24.0 ||
        dense_solve_xors_estimate(6, 4) != 36.0) {
        return self_fail(why, "dense solve XOR estimate is wrong");
    }

    {
        PeelSim sim(2);
        std::string detail;
        if (!sim.feed_columns(std::vector<uint16_t>{0}, &detail) ||
            !sim.feed_columns(std::vector<uint16_t>{0, 1}, &detail)) {
            return self_fail(why, "sparse solve XOR test feed failed: " + detail);
        }
        const TrialResult result = sim.run_greedy(0, UINT32_C(0x51f15eed));
        if (result.SparseSolveXors != 1u) {
            return self_fail(why, "sparse solve XOR accounting missed known row input");
        }
    }
    return true;
}

static bool run_random_structure_generator_tests(std::string* why)
{
    const unsigned N = 64;
    for (const MatrixStructure& structure : kStructures)
    {
        for (const Method& method : kMethods)
        {
            PeelSim sim(N);
            const uint64_t seed =
                UINT64_C(0x737472756374) ^
                ((uint64_t)method.Id * UINT64_C(0x9e3779b97f4a7c15)) ^
                ((uint64_t)N * UINT64_C(0xbf58476d1ce4e5b9)) ^
                (hash_string64(structure.Name) * UINT64_C(0x94d049bb133111eb));

            std::string detail;
            if (!feed_random_structure(sim, structure, N, N, seed, &detail))
            {
                return self_fail(why, std::string("random structure feed failed for ") +
                    structure.Name + ": " + detail);
            }
            if (!sim.validate_state(false, &detail))
            {
                return self_fail(why, std::string("random structure state failed for ") +
                    structure.Name + ": " + detail);
            }

            TrialResult result;
            if (!sim.run_greedy_checked(method.Id, (uint32_t)seed, &result, &detail))
            {
                return self_fail(why, std::string("random structure checked run failed for ") +
                    structure.Name + "/" + method.Name + ": " + detail);
            }
        }
    }
    return true;
}

static bool run_random_structure_degree_tests(std::string* why)
{
    const unsigned N = 128;
    const unsigned row_count = 512;

    {
        Rng rng(UINT64_C(0x6869676864656772));
        const std::vector<uint16_t> full = random_row_columns(N, N, rng);
        if (!validate_raw_structure(
            "random_row_full_degree",
            N,
            std::vector<std::vector<uint16_t> >(1, full),
            N,
            N,
            0,
            why)) {
            return false;
        }

        const std::vector<uint16_t> high = random_row_columns(N, N / 2u, rng);
        if (!validate_raw_structure(
            "random_row_high_degree",
            N,
            std::vector<std::vector<uint16_t> >(1, high),
            N / 2u,
            N / 2u,
            0,
            why)) {
            return false;
        }
    }

    for (const MatrixStructure& structure : kStructures)
    {
        const uint64_t seed =
            UINT64_C(0x646567726565) ^
            ((uint64_t)N * UINT64_C(0xbf58476d1ce4e5b9)) ^
            (hash_string64(structure.Name) * UINT64_C(0x94d049bb133111eb));
        const std::vector<std::vector<uint16_t> > rows =
            generate_random_structure_rows(structure, N, row_count, seed);

        const bool has_degree_one_mixture =
            (structure.Kind == kStructureOnePlusUniform ||
             structure.Kind == kStructureOnePlusLT) &&
            structure.A > 0.0;
        const unsigned min_degree = has_degree_one_mixture ?
            1 : clamp_degree(structure.MinDegree, N);
        unsigned max_degree = clamp_degree(structure.MaxDegree, N);
        if (structure.Kind == kStructureTwoModeUniform &&
            structure.AltMaxDegree > max_degree) {
            max_degree = clamp_degree(structure.AltMaxDegree, N);
        }
        const unsigned min_one_rows = has_degree_one_mixture ? 1 : 0;
        if (!validate_raw_structure(structure.Name, N, rows, min_degree,
            max_degree, min_one_rows, why)) {
            return false;
        }
    }

    return true;
}

static bool run_paired_generation_tests(std::string* why)
{
    const unsigned N = 96;
    const unsigned row_count = 128;
    for (const MatrixStructure& structure : kStructures)
    {
        const uint64_t seed =
            UINT64_C(0x706169726564) ^
            ((uint64_t)N * UINT64_C(0xbf58476d1ce4e5b9)) ^
            (hash_string64(structure.Name) * UINT64_C(0x94d049bb133111eb));
        const std::vector<std::vector<uint16_t> > a =
            generate_random_structure_rows(structure, N, row_count, seed);
        const std::vector<std::vector<uint16_t> > b =
            generate_random_structure_rows(structure, N, row_count, seed);
        const std::vector<std::vector<uint16_t> > c =
            generate_random_structure_rows(structure, N, row_count,
                seed ^ UINT64_C(0x9e3779b97f4a7c15));

        if (a != b) {
            return self_fail(why, std::string(structure.Name) +
                ": paired generator is not deterministic");
        }
        if (a == c) {
            return self_fail(why, std::string(structure.Name) +
                ": changed seed did not change generated rows");
        }
    }
    return true;
}

static bool run_new_structure_variant_tests(std::string* why)
{
    const MatrixStructure* robust = find_structure("rs_c003_d10_c128");
    const MatrixStructure* alpha = find_structure("sol_alpha075_c128");
    const MatrixStructure* cap256 = find_structure("lt_m2_c256");
    const MatrixStructure* cap512 = find_structure("lt_m2_c512");
    const MatrixStructure* cap1024 = find_structure("lt_m2_c1024");
    const MatrixStructure* cap3200 = find_structure("lt_m2_c3200");
    const MatrixStructure* fold = find_structure("lt_m2_c128_fold");
    const MatrixStructure* fold512 = find_structure("lt_m2_c512_fold");
    const MatrixStructure* fold1024 = find_structure("lt_m2_c1024_fold");
    const MatrixStructure* fold2560 = find_structure("lt_m2_c2560_fold");
    const MatrixStructure* fold3200 = find_structure("lt_m2_c3200_fold");
    const MatrixStructure* fold0 = find_structure("lt_m2_c256_fold0");
    const MatrixStructure* fold50 = find_structure("lt_m2_c256_fold50");
    const MatrixStructure* fold200 = find_structure("lt_m2_c256_fold200");
    const MatrixStructure* fold1024_0 = find_structure("lt_m2_c1024_fold0");
    const MatrixStructure* fold1024_50 = find_structure("lt_m2_c1024_fold50");
    const MatrixStructure* fold1024_200 = find_structure("lt_m2_c1024_fold200");
    const MatrixStructure* alpha3200 = find_structure("sol_m2_alpha125_c3200");
    const MatrixStructure* extra34 = find_structure("lt_d3_008_d4_008");
    const MatrixStructure* extra34_min2 =
        find_structure("lt_m2_c320_d2_003_d3_003_d4_003");
    const MatrixStructure* extra34_min2_3200 =
        find_structure("lt_m2_c3200_d2_003_d3_003_d4_003");
    const MatrixStructure* high_half = find_structure("lt_m2_c320_hi4_h160");
    const MatrixStructure* high_full = find_structure("lt_m2_c320_hi4_h320");
    const MatrixStructure* high_half_3200 =
        find_structure("lt_m2_c3200_hi4_h1600");
    const MatrixStructure* high_full_3200 =
        find_structure("lt_m2_c3200_hi4_h3200");
    const MatrixStructure* overhead = find_structure("ohdep_lt_adapt128");
    const MatrixStructure* ldpc = find_structure("raptorq_ldpc_struct");
    const MatrixStructure* ldpc_random = find_structure("raptorq_ldpc_rand");
    const MatrixStructure* binary = find_structure("binary_p50");
    if (!robust || !alpha || !cap256 || !cap512 || !cap1024 ||
        !cap3200 || !fold || !fold512 || !fold1024 || !fold2560 ||
        !fold3200 || !fold0 || !fold50 || !fold200 || !fold1024_0 ||
        !fold1024_50 || !fold1024_200 || !alpha3200 || !extra34 ||
        !extra34_min2 || !extra34_min2_3200 || !high_half || !high_full ||
        !high_half_3200 || !high_full_3200 || !overhead || !ldpc ||
        !ldpc_random || !binary) {
        return self_fail(why, "new structure variant lookup failed");
    }

    const unsigned N = 320;
    const unsigned row_count = 320;
    const uint64_t seed = UINT64_C(0x6e65777661726961);
    const MatrixStructure* structures[] = {
        robust, alpha, alpha3200, cap256, cap512, cap1024, cap3200,
        fold, fold512, fold1024, fold2560, fold3200, fold0, fold50,
        fold200, fold1024_0, fold1024_50, fold1024_200, extra34,
        extra34_min2, extra34_min2_3200, high_half, high_full,
        high_half_3200, high_full_3200, overhead, ldpc, ldpc_random, binary
    };
    for (const MatrixStructure* structure : structures)
    {
        const std::vector<std::vector<uint16_t> > rows =
            generate_random_structure_rows(*structure, N, row_count,
                seed ^ hash_string64(structure->Name));
        unsigned max_degree = clamp_degree(structure->MaxDegree, N);
        if (is_raptorq_ldpc_structure(*structure)) {
            max_degree = N;
        }
        if (!validate_raw_structure(structure->Name, N, rows, 1,
            max_degree, 0, why)) {
            return false;
        }
    }

    const unsigned cap = clamp_degree(fold->MaxDegree, N);
    const double folded = structure_degree_weight(
        *fold, cap, N, row_count, cap);
    if (folded <= lt_weight(cap, N)) {
        return self_fail(why, "cap-fold structure did not fold tail mass");
    }

    const unsigned larger_n = 3200;
    const unsigned larger_rows = 3200;
    const unsigned cap512_value = clamp_degree(fold512->MaxDegree, larger_n);
    const unsigned cap1024_value = clamp_degree(fold1024->MaxDegree, larger_n);
    const double folded512 = structure_degree_weight(
        *fold512, cap512_value, larger_n, larger_rows, cap512_value);
    const double folded1024 = structure_degree_weight(
        *fold1024, cap1024_value, larger_n, larger_rows, cap1024_value);
    if (folded512 <= lt_weight(cap512_value, larger_n) ||
        folded1024 <= lt_weight(cap1024_value, larger_n)) {
        return self_fail(why, "large cap-fold structures did not fold tail mass");
    }

    const unsigned scale_cap = clamp_degree(fold50->MaxDegree, N);
    const double scale0 = structure_degree_weight(
        *fold0, scale_cap, N, row_count, scale_cap);
    const double scale50 = structure_degree_weight(
        *fold50, scale_cap, N, row_count, scale_cap);
    const double scale200 = structure_degree_weight(
        *fold200, scale_cap, N, row_count, scale_cap);
    if (!(scale0 < scale50 && scale50 < scale200)) {
        return self_fail(why, "scaled fold weights are not monotonic");
    }

    const unsigned scale1024_cap = clamp_degree(fold1024_50->MaxDegree, larger_n);
    const double scale1024_0 = structure_degree_weight(
        *fold1024_0, scale1024_cap, larger_n, larger_rows, scale1024_cap);
    const double scale1024_50 = structure_degree_weight(
        *fold1024_50, scale1024_cap, larger_n, larger_rows, scale1024_cap);
    const double scale1024_200 = structure_degree_weight(
        *fold1024_200, scale1024_cap, larger_n, larger_rows, scale1024_cap);
    if (!(scale1024_0 < scale1024_50 && scale1024_50 < scale1024_200)) {
        return self_fail(why, "scaled 1024 fold weights are not monotonic");
    }

    const MatrixStructure* fold512_10 = find_structure("lt_m2_c512_fold10");
    const MatrixStructure* fold512_40 = find_structure("lt_m2_c512_fold40");
    const MatrixStructure* fold512_60 = find_structure("lt_m2_c512_fold60");
    const MatrixStructure* fold2048_25 = find_structure("lt_m2_c2048_fold25");
    const MatrixStructure* grid_d2 = find_structure("lt_m2_c1024_d2_003");
    const MatrixStructure* grid_d3 = find_structure("lt_m2_c2048_d3_003");
    if (!fold512_10 || !fold512_40 || !fold512_60 || !fold2048_25 ||
        !grid_d2 || !grid_d3) {
        return self_fail(why, "fold/cap fine-grid structure lookup failed");
    }

    const unsigned grid_cap = clamp_degree(fold512_10->MaxDegree, larger_n);
    const double grid512_10 = structure_degree_weight(
        *fold512_10, grid_cap, larger_n, larger_rows, grid_cap);
    const double grid512_40 = structure_degree_weight(
        *fold512_40, grid_cap, larger_n, larger_rows, grid_cap);
    const double grid512_60 = structure_degree_weight(
        *fold512_60, grid_cap, larger_n, larger_rows, grid_cap);
    if (!(grid512_10 < grid512_40 && grid512_40 < grid512_60)) {
        return self_fail(why, "fine-grid 512 fold weights are not monotonic");
    }
    if (!(grid512_10 > lt_weight(grid_cap, larger_n))) {
        return self_fail(why, "fine-grid cap 512 fold10 did not fold tail mass");
    }

    const unsigned grid2048_cap = clamp_degree(fold2048_25->MaxDegree, larger_n);
    const double grid2048_25 = structure_degree_weight(
        *fold2048_25, grid2048_cap, larger_n, larger_rows, grid2048_cap);
    if (!(grid2048_25 > lt_weight(grid2048_cap, larger_n))) {
        return self_fail(why, "fine-grid cap 2048 fold25 did not fold tail mass");
    }

    const double grid_d2_w = structure_degree_weight(
        *grid_d2, 2, larger_n, larger_rows,
        clamp_degree(grid_d2->MaxDegree, larger_n));
    const double grid_d3_w = structure_degree_weight(
        *grid_d3, 3, larger_n, larger_rows,
        clamp_degree(grid_d3->MaxDegree, larger_n));
    if (!(grid_d2->MinDegree == 2 && grid_d2_w > lt_weight(2, larger_n))) {
        return self_fail(why, "fine-grid cap 1024 d2 mass variant is malformed");
    }
    if (!(grid_d3->MinDegree == 2 && grid_d3_w > lt_weight(3, larger_n))) {
        return self_fail(why, "fine-grid cap 2048 d3 mass variant is malformed");
    }

    const double robust_degree2 = structure_degree_weight(
        *robust, 2, N, row_count, clamp_degree(robust->MaxDegree, N));
    if (!(robust_degree2 > lt_weight(2, N))) {
        return self_fail(why, "robust soliton did not add robust tail mass");
    }

    const double alpha_low = structure_degree_weight(
        *alpha, 2, N, row_count, clamp_degree(alpha->MaxDegree, N));
    const double alpha_high = structure_degree_weight(
        *alpha, 32, N, row_count, clamp_degree(alpha->MaxDegree, N));
    if (!(alpha_low > 0.0 && alpha_high > 0.0)) {
        return self_fail(why, "alpha soliton weights are nonpositive");
    }

    const double extra_d3 = structure_degree_weight(
        *extra34, 3, N, row_count, clamp_degree(extra34->MaxDegree, N));
    const double extra_d5 = structure_degree_weight(
        *extra34, 5, N, row_count, clamp_degree(extra34->MaxDegree, N));
    if (!(extra_d3 > extra_d5)) {
        return self_fail(why, "degree-3/4 explicit mass did not affect degree 3");
    }

    const double extra_min2_d2 = structure_degree_weight(
        *extra34_min2, 2, N, row_count,
        clamp_degree(extra34_min2->MaxDegree, N));
    if (!(extra34_min2->MinDegree == 2 && extra_min2_d2 > lt_weight(2, N))) {
        return self_fail(why, "min-2 degree-mass variant is malformed");
    }

    const double extra_min2_3200_d2 = structure_degree_weight(
        *extra34_min2_3200, 2, larger_n, larger_rows,
        clamp_degree(extra34_min2_3200->MaxDegree, larger_n));
    if (!(extra34_min2_3200->MinDegree == 2 &&
        extra_min2_3200_d2 > lt_weight(2, larger_n))) {
        return self_fail(why, "min-2 N=3200 degree-mass variant is malformed");
    }

    const double high_base_d2 = structure_degree_weight(
        *high_half, 2, N, row_count, clamp_degree(high_half->MaxDegree, N));
    if (!(high_base_d2 == lt_weight(2, N))) {
        return self_fail(why, "explicit high-row base distribution is not LT");
    }

    const double oh_zero = structure_degree_weight(
        *overhead, 1, N, N, clamp_degree(overhead->MaxDegree, N));
    const double oh_high = structure_degree_weight(
        *overhead, 1, N, N + 16, clamp_degree(overhead->MaxDegree, N));
    if (!(oh_high > oh_zero)) {
        return self_fail(why, "overhead-dependent PDF did not change degree-1 mass");
    }

    const std::vector<std::vector<uint16_t> > binary_rows =
        generate_random_structure_rows(*binary, 128, 256, seed);
    unsigned binary_refs = 0;
    for (const std::vector<uint16_t>& row : binary_rows) {
        binary_refs += (unsigned)row.size();
    }
    const double binary_density =
        (double)binary_refs / (double)(128u * 256u);
    if (binary_density < 0.40 || binary_density > 0.60) {
        return self_fail(why, "binary_p50 density is outside smoke-test bounds");
    }

    const std::vector<std::vector<uint16_t> > high_half_rows =
        generate_random_structure_rows(*high_half, N, row_count, seed);
    const std::vector<std::vector<uint16_t> > high_full_rows =
        generate_random_structure_rows(*high_full, N, row_count, seed);
    for (unsigned i = 0; i < 4u; ++i)
    {
        if (high_half_rows[i].size() != N / 2u ||
            high_full_rows[i].size() != N) {
            return self_fail(why, "explicit high-degree rows have wrong degree");
        }
    }

    const std::vector<std::vector<uint16_t> > high_half_3200_rows =
        generate_random_structure_rows(
            *high_half_3200, larger_n, larger_rows, seed);
    const std::vector<std::vector<uint16_t> > high_full_3200_rows =
        generate_random_structure_rows(
            *high_full_3200, larger_n, larger_rows, seed);
    for (unsigned i = 0; i < 4u; ++i)
    {
        if (high_half_3200_rows[i].size() != larger_n / 2u ||
            high_full_3200_rows[i].size() != larger_n) {
            return self_fail(why, "explicit high-degree rows do not scale with N");
        }
    }

    const unsigned S = raptorq_ldpc_symbol_count(N);
    const std::vector<std::vector<uint16_t> > ldpc_rows =
        generate_random_structure_rows(*ldpc, N, 8, seed);
    const std::vector<std::vector<uint16_t> > random_rows =
        generate_random_structure_rows(*ldpc_random, N, 8, seed);
    if (ldpc_rows.size() != 8u + S || random_rows.size() != 8u + S) {
        return self_fail(why, "LDPC structures did not add expected check rows");
    }
    if (ldpc_rows == random_rows) {
        return self_fail(why, "LDPC randomized control matched structured rows");
    }
    for (unsigned i = 0; i < S; ++i)
    {
        if (ldpc_rows[i].size() < 3u) {
            return self_fail(why, "structured LDPC check row is too small");
        }
        if (random_rows[i].size() != ldpc_rows[i].size()) {
            return self_fail(why, "LDPC randomized control changed check degree");
        }
    }
    for (unsigned i = S; i < ldpc_rows.size(); ++i)
    {
        if (ldpc_rows[i] != random_rows[i]) {
            return self_fail(why, "LDPC control did not preserve paired packet rows");
        }
    }

    return true;
}

static unsigned column_usage_spread(
    const std::vector<std::vector<uint16_t> >& rows,
    unsigned N)
{
    std::vector<unsigned> counts(N, 0);
    for (const std::vector<uint16_t>& row : rows)
    {
        for (uint16_t column_i : row) {
            ++counts[column_i];
        }
    }
    unsigned min_count = counts[0];
    unsigned max_count = counts[0];
    for (unsigned column_i = 1; column_i < N; ++column_i)
    {
        if (counts[column_i] < min_count) {
            min_count = counts[column_i];
        }
        if (counts[column_i] > max_count) {
            max_count = counts[column_i];
        }
    }
    return max_count - min_count;
}

static uint64_t total_row_refs(const std::vector<std::vector<uint16_t> >& rows)
{
    uint64_t refs = 0;
    for (const std::vector<uint16_t>& row : rows) {
        refs += (uint64_t)row.size();
    }
    return refs;
}

static bool run_deck_balance_tests(std::string* why)
{
    struct DeckCase
    {
        const char* DeckName;
        const char* BaseName;
    };
    const DeckCase deck_cases[] = {
        {"deck_lt_m2_c256_fold", "lt_m2_c256_fold"},
        {"deck_lt_m2_c512_fold", "lt_m2_c512_fold"},
        {"deck_lt_m2_c128", "lt_m2_c128"},
        {"deck_lt_m2_c16", "lt_m2_c16"},
        {"deck_lt_m1_c64", "lt_m1_c64"},
        {"deck_rs_c001_d50_c128", "rs_c001_d50_c128"},
    };

    const unsigned N = 320;
    const unsigned row_count = 352;
    for (const DeckCase& deck_case : deck_cases)
    {
        const MatrixStructure* deck = find_structure(deck_case.DeckName);
        const MatrixStructure* base = find_structure(deck_case.BaseName);
        if (!deck || !base) {
            return self_fail(why, std::string(deck_case.DeckName) +
                ": deck structure lookup failed");
        }
        if (deck->Kind != kStructureDeckColumns ||
            deck->DegreeKind != base->Kind ||
            deck->MinDegree != base->MinDegree ||
            deck->MaxDegree != base->MaxDegree ||
            deck->A != base->A || deck->B != base->B) {
            return self_fail(why, std::string(deck_case.DeckName) +
                ": deck structure does not mirror its base degree law");
        }

        const uint64_t seed = UINT64_C(0x6465636b6465616c) ^
            hash_string64(deck_case.DeckName);
        const std::vector<std::vector<uint16_t> > rows =
            generate_random_structure_rows(*deck, N, row_count, seed);
        const std::vector<std::vector<uint16_t> > again =
            generate_random_structure_rows(*deck, N, row_count, seed);
        if (rows != again) {
            return self_fail(why, std::string(deck_case.DeckName) +
                ": deck generator is not deterministic");
        }
        if (rows.size() != row_count) {
            return self_fail(why, std::string(deck_case.DeckName) +
                ": deck generator changed the row count");
        }
        // Also checks row degree bounds, column range, and that no row
        // contains the same column twice.
        if (!validate_raw_structure(deck_case.DeckName, N, rows,
            clamp_degree(deck->MinDegree, N),
            clamp_degree(deck->MaxDegree, N), 0, why)) {
            return false;
        }
        if (column_usage_spread(rows, N) > 1u) {
            return self_fail(why, std::string(deck_case.DeckName) +
                ": deck column degrees are not balanced within 1");
        }

        // The deck deal must spend exactly the edge budget an iid generator
        // would draw from the same degree-law stream: every row degree has
        // to match the law sample, so total edges match the iid equivalent.
        Rng law_rng(seed);
        WeightedDegreeSampler law_sampler;
        if (uses_weighted_degree_sampler(*deck)) {
            law_sampler.reset(*deck, N, row_count);
        }
        uint64_t expected_refs = 0;
        for (unsigned row_i = 0; row_i < row_count; ++row_i)
        {
            const unsigned degree = clamp_degree(
                choose_structure_degree(*deck, N, law_sampler, law_rng), N);
            if (rows[row_i].size() != degree) {
                return self_fail(why, std::string(deck_case.DeckName) +
                    ": deck row degree does not match the degree-law sample");
            }
            expected_refs += degree;
        }
        if (total_row_refs(rows) != expected_refs) {
            return self_fail(why, std::string(deck_case.DeckName) +
                ": deck edge budget diverged from the iid degree law");
        }
    }

    const char* loss_names[] = {
        "deckloss10_lt_m2_c256_fold",
        "deckloss30_lt_m2_c256_fold",
    };
    for (const char* name : loss_names)
    {
        const MatrixStructure* loss = find_structure(name);
        if (!loss || loss->Kind != kStructureDeckColumns) {
            return self_fail(why, std::string(name) +
                ": deck loss structure lookup failed");
        }
        if (!(loss->BalanceValue > 1.0)) {
            return self_fail(why, std::string(name) +
                ": deck loss virtual row multiplier is missing");
        }
        if (deck_virtual_row_count(*loss, row_count) <= row_count) {
            return self_fail(why, std::string(name) +
                ": deck loss did not add virtual rows");
        }

        const uint64_t seed = UINT64_C(0x6465636b6c6f7373) ^
            hash_string64(name);
        const std::vector<std::vector<uint16_t> > rows =
            generate_random_structure_rows(*loss, N, row_count, seed);
        const std::vector<std::vector<uint16_t> > again =
            generate_random_structure_rows(*loss, N, row_count, seed);
        if (rows != again) {
            return self_fail(why, std::string(name) +
                ": deck loss generator is not deterministic");
        }
        if (rows.size() != row_count) {
            return self_fail(why, std::string(name) +
                ": deck loss generator kept the wrong row count");
        }
        if (!validate_raw_structure(name, N, rows,
            clamp_degree(loss->MinDegree, N),
            clamp_degree(loss->MaxDegree, N), 0, why)) {
            return false;
        }
        // Dropping virtual rows breaks exact balance, but the surviving
        // spread should stay far below an adversarial blowup.
        if (column_usage_spread(rows, N) > 32u) {
            return self_fail(why, std::string(name) +
                ": deck loss column spread is implausibly large");
        }
    }

    const MatrixStructure* p2c = find_structure("p2c_lt_m2_c256_fold");
    const MatrixStructure* p2c_base = find_structure("lt_m2_c256_fold");
    if (!p2c || !p2c_base) {
        return self_fail(why, "p2c structure lookup failed");
    }
    if (p2c->BalanceMode != 3u || p2c->Kind != p2c_base->Kind ||
        p2c->MinDegree != p2c_base->MinDegree ||
        p2c->MaxDegree != p2c_base->MaxDegree) {
        return self_fail(why, "p2c structure does not mirror its base degree law");
    }

    const uint64_t p2c_seed = UINT64_C(0x7032636465616c) ^
        hash_string64(p2c->Name);
    const std::vector<std::vector<uint16_t> > p2c_rows =
        generate_random_structure_rows(*p2c, N, row_count, p2c_seed);
    const std::vector<std::vector<uint16_t> > p2c_again =
        generate_random_structure_rows(*p2c, N, row_count, p2c_seed);
    if (p2c_rows != p2c_again) {
        return self_fail(why, "p2c generator is not deterministic");
    }
    if (p2c_rows.size() != row_count) {
        return self_fail(why, "p2c generator changed the row count");
    }
    if (!validate_raw_structure(p2c->Name, N, p2c_rows,
        clamp_degree(p2c->MinDegree, N),
        clamp_degree(p2c->MaxDegree, N), 0, why)) {
        return false;
    }
    const std::vector<std::vector<uint16_t> > iid_rows =
        generate_random_structure_rows(*p2c_base, N, row_count, p2c_seed);
    if (column_usage_spread(p2c_rows, N) > column_usage_spread(iid_rows, N)) {
        return self_fail(why, "p2c column spread exceeds the iid spread");
    }

    return true;
}

static bool run_n_jitter_protocol_tests(std::string* why)
{
    const MatrixStructure* structure = find_structure("lt_m2_c64");
    if (!structure) {
        return self_fail(why, "N-jitter test structure is missing");
    }

    const unsigned anchor_n = 320;
    const unsigned jitter = 10;
    const uint64_t seed = UINT64_C(0x4e6a6974746572);
    std::vector<unsigned> seen_ns;
    for (int trial = 0; trial < 64; ++trial)
    {
        const TrialConfig a = make_trial_config(
            0, anchor_n, trial, 0, 0, 2.0, jitter, seed,
            *structure, "paired");
        const TrialConfig b = make_trial_config(
            0, anchor_n, trial, 0, 0, 2.0, jitter, seed,
            *structure, "paired");
        const TrialConfig other_method = make_trial_config(
            7, anchor_n, trial, 0, 0, 2.0, jitter, seed,
            *structure, "paired");
        const TrialConfig independent_method = make_trial_config(
            7, anchor_n, trial, 0, 0, 2.0, jitter, seed,
            *structure, "independent");

        if (a.ActualN < anchor_n - jitter || a.ActualN > anchor_n + jitter) {
            return self_fail(why, "N-jitter sampled outside anchor window");
        }
        if (a.ActualN != b.ActualN || a.MatrixRows != b.MatrixRows ||
            a.MatrixSeed != b.MatrixSeed || a.ActionSeed != b.ActionSeed) {
            return self_fail(why, "N-jitter trial config is not deterministic");
        }
        if (a.ActualN != other_method.ActualN ||
            a.MatrixRows != other_method.MatrixRows ||
            a.MatrixSeed != other_method.MatrixSeed) {
            return self_fail(why, "paired N-jitter config differs across methods");
        }
        if (a.ActionSeed == other_method.ActionSeed) {
            return self_fail(why, "paired methods reused the same action seed");
        }
        if (a.MatrixSeed == independent_method.MatrixSeed) {
            return self_fail(why, "independent N-jitter config reused matrix seed");
        }

        const std::vector<std::vector<uint16_t> > rows_a =
            generate_random_structure_rows(
                *structure, a.ActualN, a.MatrixRows, a.MatrixSeed);
        const std::vector<std::vector<uint16_t> > rows_other =
            generate_random_structure_rows(
                *structure, other_method.ActualN, other_method.MatrixRows,
                other_method.MatrixSeed);
        if (rows_a != rows_other) {
            return self_fail(why, "paired N-jitter matrices differ across methods");
        }

        bool already_seen = false;
        for (unsigned n : seen_ns)
        {
            if (n == a.ActualN)
            {
                already_seen = true;
                break;
            }
        }
        if (!already_seen) {
            seen_ns.push_back(a.ActualN);
        }
    }
    if (seen_ns.size() < 2u) {
        return self_fail(why, "N-jitter did not vary across trials");
    }

    const TrialConfig fixed = make_trial_config(
        0, anchor_n, 0, 0, 0, 0.0, 0, seed, *structure, "paired");
    if (fixed.ActualN != anchor_n) {
        return self_fail(why, "zero N-jitter changed anchor N");
    }
    return true;
}

static bool run_overhead_pct_validation_tests(std::string* why)
{
    bool ok = false;
    std::vector<double> parsed = parse_double_list("0,12.5", &ok);
    if (!ok || parsed.size() != 2u || parsed[0] != 0.0 || parsed[1] != 12.5) {
        return self_fail(why, "finite overhead-pct list did not parse");
    }

    parsed = parse_double_list("nan", &ok);
    if (ok) {
        return self_fail(why, "nan overhead-pct parsed as valid");
    }
    parsed = parse_double_list("inf", &ok);
    if (ok) {
        return self_fail(why, "inf overhead-pct parsed as valid");
    }

    const MatrixStructure* structure = find_structure("lt_m2_c64");
    if (!structure) {
        return self_fail(why, "overhead-pct test structure is missing");
    }

    const TrialConfig valid = make_trial_config(
        0, 320, 0, 0, 0, 10.0, 0, UINT64_C(0x6f76657268656164),
        *structure, "paired");
    if (valid.ActualN != 320u || valid.MatrixRows != 352u) {
        return self_fail(why, "finite overhead-pct row count is wrong");
    }

    const TrialConfig huge = make_trial_config(
        0, 64000, 0, 0, 0, 1.0e9, 0, UINT64_C(0x6f76657268656164),
        *structure, "paired");
    if (huge.MatrixRows <= 65535u) {
        return self_fail(why, "huge overhead-pct did not trip row limit");
    }
    return true;
}

static bool run_lazy_metric_tests(std::string* why)
{
    const MatrixStructure* structure = find_structure("lt_m2_c64");
    if (!structure) {
        return self_fail(why, "lazy metric test structure is missing");
    }

    PeelSim sim(96);
    std::string detail;
    if (!feed_random_structure(
        sim, *structure, 96, 112, UINT64_C(0x6c617a796d657472), &detail)) {
        return self_fail(why, "lazy metric structure feed failed: " + detail);
    }
    if (!sim.validate_lazy_equivalence(UINT32_C(0x51f15eed), &detail)) {
        return self_fail(why, "initial lazy metric equivalence failed: " + detail);
    }
    if (!sim.run_lazy_equivalence_walk(UINT32_C(0x51f15eed), &detail)) {
        return self_fail(why, "lazy metric invalidation walk failed: " + detail);
    }
    return true;
}

static bool run_avalanche_estimate_tests(std::string* why)
{
    const MatrixStructure* structure = find_structure("lt_m2_c128");
    if (!structure) {
        return self_fail(why, "avalanche estimate test structure is missing");
    }

    PeelSim sim(160);
    std::string detail;
    if (!feed_random_structure(
        sim, *structure, 160, 160, UINT64_C(0x6176616c616e6368), &detail)) {
        return self_fail(why, "avalanche estimate structure feed failed: " + detail);
    }
    if (!sim.validate_avalanche_estimate(UINT32_C(0x51f15eed), &detail)) {
        return self_fail(why, "initial avalanche estimate validation failed: " + detail);
    }

    TrialResult result;
    if (!sim.run_greedy_checked(18, UINT32_C(0x51f15eed), &result, &detail)) {
        return self_fail(why, "avalanche estimate walk failed: " + detail);
    }

    PeelSim tail(160);
    if (!feed_random_structure(
        tail, *structure, 160, 160, UINT64_C(0x6176616c616e6368), &detail)) {
        return self_fail(why, "avalanche estimate tail feed failed: " + detail);
    }
    if (!tail.run_avalanche_estimate_walk(
        18, UINT32_C(0x51f15eed), 8, &detail)) {
        return self_fail(why, "avalanche estimate step validation failed: " + detail);
    }
    return true;
}

static bool run_lazy_avalanche_pair_tests(std::string* why)
{
    struct MethodPair
    {
        int ExactMethod;
        int LazyMethod;
        const char* Label;
    };
    static const MethodPair kPairs[] = {
        {119, 127, "av_top16/av_lazy16"},
        {120, 128, "av_top64/av_lazy64"},
    };
    static const char* const kStructureNames[] = {
        "wirehair_rand", "lt_m2_c1024", "lt_m2_c256_fold"
    };
    static const unsigned kAnchors[] = {320, 1280, 3200};
    const int pair_trials = 2;
    const uint64_t seed = UINT64_C(0x6c617a7961766131);

    for (const char* structure_name : kStructureNames)
    {
        const MatrixStructure* structure = find_structure(structure_name);
        if (!structure)
        {
            return self_fail(why,
                std::string("lazy avalanche test structure is missing: ") +
                structure_name);
        }

        for (unsigned anchor_n : kAnchors)
        {
            for (int trial = 0; trial < pair_trials; ++trial)
            {
                const TrialConfig config = make_trial_config(
                    0, anchor_n, trial, 0, 0, 0.0, 10, seed,
                    *structure, "paired");
                const std::vector<std::vector<uint16_t> > rows =
                    generate_random_structure_rows(*structure,
                        config.ActualN, config.MatrixRows, config.MatrixSeed);

                for (const MethodPair& pair : kPairs)
                {
                    std::string detail;
                    PeelSim exact_sim(config.ActualN);
                    PeelSim lazy_sim(config.ActualN);
                    if (!feed_row_vectors(exact_sim, rows, &detail) ||
                        !feed_row_vectors(lazy_sim, rows, &detail))
                    {
                        return self_fail(why,
                            "lazy avalanche pair feed failed: " + detail);
                    }

                    TrialResult exact_result;
                    TrialResult lazy_result;
                    std::vector<uint16_t> exact_trace;
                    std::vector<uint16_t> lazy_trace;
                    if (!exact_sim.run_greedy_traced(pair.ExactMethod,
                        config.ActionSeed, &exact_result, &exact_trace,
                        &detail))
                    {
                        return self_fail(why,
                            "exact avalanche traced run failed: " + detail);
                    }
                    if (!lazy_sim.run_greedy_traced(pair.LazyMethod,
                        config.ActionSeed, &lazy_result, &lazy_trace,
                        &detail))
                    {
                        return self_fail(why,
                            "lazy avalanche traced run failed: " + detail);
                    }

                    if (exact_trace != lazy_trace)
                    {
                        return self_fail(why, std::string(pair.Label) +
                            ": selected column sequences differ");
                    }
                    if (exact_result.ResidualCols != lazy_result.ResidualCols ||
                        exact_result.ResidualRows != lazy_result.ResidualRows ||
                        exact_result.Choices != lazy_result.Choices ||
                        exact_result.Components != lazy_result.Components ||
                        exact_result.MaxComponent != lazy_result.MaxComponent ||
                        exact_result.SumSquares != lazy_result.SumSquares ||
                        exact_result.SparseSolveXors != lazy_result.SparseSolveXors ||
                        exact_result.InitialTodo != lazy_result.InitialTodo ||
                        exact_result.InitialD2Rows != lazy_result.InitialD2Rows)
                    {
                        return self_fail(why, std::string(pair.Label) +
                            ": traced results differ");
                    }
                    if (!exact_sim.validate_state(true, &detail) ||
                        !lazy_sim.validate_state(true, &detail))
                    {
                        return self_fail(why,
                            "lazy avalanche pair final state invalid: " + detail);
                    }
                }
            }
        }
    }
    return true;
}

static bool run_self_tests(std::string* why)
{
    if (!run_synthetic_structure_tests(why)) {
        return false;
    }
    if (!run_matrix_work_tests(why)) {
        return false;
    }
    if (!run_random_structure_generator_tests(why)) {
        return false;
    }
    if (!run_random_structure_degree_tests(why)) {
        return false;
    }
    if (!run_paired_generation_tests(why)) {
        return false;
    }
    if (!run_new_structure_variant_tests(why)) {
        return false;
    }
    if (!run_deck_balance_tests(why)) {
        return false;
    }
    if (!run_n_jitter_protocol_tests(why)) {
        return false;
    }
    if (!run_overhead_pct_validation_tests(why)) {
        return false;
    }
    if (!run_lazy_metric_tests(why)) {
        return false;
    }
    if (!run_avalanche_estimate_tests(why)) {
        return false;
    }
    if (!run_lazy_avalanche_pair_tests(why)) {
        return false;
    }
    if (!PeelSim::validate_multistart_tie_break(why)) {
        return false;
    }
    if (!PeelSim::validate_minrow_random_tie(why)) {
        return false;
    }
    return true;
}

static TrialResult run_trial(
    int method_id,
    unsigned N,
    unsigned row_count,
    const MatrixStructure& structure,
    uint64_t matrix_seed,
    uint32_t action_seed)
{
    const double start = now_sec();

    PeelSim sim(N);

    const double build_start = now_sec();
    std::string why;
    const std::vector<std::vector<uint16_t> > rows =
        generate_random_structure_rows(structure, N, row_count, matrix_seed);
    const MatrixWork work = measure_matrix_work(rows);
    if (!feed_row_vectors(sim, rows, &why))
    {
        std::fprintf(stderr, "random structure generation failed: %s\n",
            why.c_str());
        std::abort();
    }
    const unsigned build_usec = (unsigned)((now_sec() - build_start) * 1000000.0 + 0.5);

    TrialResult result = sim.run_greedy(method_id, action_seed);
    result.ActualN = N;
    result.MatrixRows = (unsigned)rows.size();
    result.MatrixRefs = work.Refs;
    result.MatrixXors = work.Xors;
    result.DenseSolveXorsEst = dense_solve_xors_estimate(
        result.ResidualCols, result.ResidualRows);
    result.TotalSolveXorsEst =
        (double)result.SparseSolveXors + result.DenseSolveXorsEst;
    result.BuildUsec = build_usec;
    result.TotalUsec = (unsigned)((now_sec() - start) * 1000000.0 + 0.5);
    return result;
}

static std::vector<int> expand_methods(const std::string& spec, bool* ok = nullptr)
{
    if (ok) {
        *ok = true;
    }
    if (spec == "all")
    {
        std::vector<int> ids;
        for (const Method& method : kMethods) {
            ids.push_back(method.Id);
        }
        return ids;
    }
    if (spec == "legacy")
    {
        std::vector<int> ids;
        for (const Method& method : kMethods)
        {
            if (method.Id <= 20) {
                ids.push_back(method.Id);
            }
        }
        return ids;
    }
    if (spec == "combo")
    {
        std::vector<int> ids;
        for (const Method& method : kMethods)
        {
            if (method.Id > 20) {
                ids.push_back(method.Id);
            }
        }
        return ids;
    }
    if (spec == "fast")
    {
        std::vector<int> ids;
        for (const Method& method : kMethods)
        {
            if (method.Id <= 20)
            {
                if (method.Id == 0 || method.Id == 10 || method.Id == 11 ||
                    method.Id == 16 || method.Id == 17 || method.Id == 18 ||
                    method.Id == 19) {
                    ids.push_back(method.Id);
                }
            }
            else if (method.Id != 58 && method.Id != 59 && method.Id < 66) {
                ids.push_back(method.Id);
            }
        }
        return ids;
    }
    return parse_int_list(spec, ok);
}

static std::vector<const MatrixStructure*> expand_structures(
    const std::string& spec,
    bool* ok = nullptr)
{
    if (ok) {
        *ok = true;
    }

    std::vector<const MatrixStructure*> out;
    if (spec == "all")
    {
        for (const MatrixStructure& structure : kStructures) {
            out.push_back(&structure);
        }
        return out;
    }
    if (spec == "random")
    {
        for (const MatrixStructure& structure : kStructures) {
            out.push_back(&structure);
        }
        return out;
    }

    size_t start = 0;
    while (start < spec.size())
    {
        size_t comma = spec.find(',', start);
        const std::string part = spec.substr(
            start,
            comma == std::string::npos ? std::string::npos : comma - start);
        if (part.empty())
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }

        const MatrixStructure* structure = find_structure(part);
        if (!structure)
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }
        out.push_back(structure);

        if (comma == std::string::npos) {
            break;
        }
        if (comma + 1 == spec.size())
        {
            if (ok) {
                *ok = false;
            }
            return out;
        }
        start = comma + 1;
    }
    return out;
}

static void print_methods()
{
    std::printf("# methods:\n");
    for (const Method& method : kMethods) {
        std::printf("#   %2d %-16s %s\n", method.Id, method.Name, method.Cost);
    }
}

static void print_structures()
{
    std::printf("# structures:\n");
    for (const MatrixStructure& structure : kStructures) {
        std::printf("#   %-14s %s\n", structure.Name, structure.Description);
    }
}

static void print_dimensions()
{
    std::printf("# candidate pools:\n");
    std::printf("#   all columns: every currently undecided column\n");
    std::printf("#   top-K default: best cached degree-2/reference-count columns\n");
    std::printf("#   degree-2 rows: endpoints of all current degree-2 rows\n");
    std::printf("#   largest degree-2 component: RaptorQ/Karp-Sipser-style component pool\n");
    std::printf("#   first minimum row: columns in the first minimum-live-degree row\n");
    std::printf("#   all minimum rows: columns across every minimum-live-degree row\n");
    std::printf("#   row actions: Burshtein-Miller-style all-but-one row inactivation\n");
    std::printf("# score statistics:\n");
    std::printf("#   cached degree-2 references, exact live degree-2 references\n");
    std::printf("#   partner lookahead, fill sum, fill square, max row degree\n");
    std::printf("#   degree-2/fill ratio, duplicate partners, distinct partners\n");
    std::printf("#   total references, live references, degree-3/degree-4 support\n");
    std::printf("#   seeded randomized tie-breaking\n");
    std::printf("#   local row-action lookahead\n");
}

} // namespace

int main(int argc, char** argv)
{
    std::string nlist = "320,3200,32000";
    std::string methods_spec = "all";
    std::string structures_spec = "all";
    std::string overhead_pct_spec = "0";
    std::string matrix_seed_mode = "paired";
    int trials = 40;
    int threads = default_threads();
    int overhead = 0;
    int rows_override = 0;
    int n_jitter = 10;
    uint64_t seed = UINT64_C(0x9a7e11a);
    bool list_methods = false;
    bool list_dimensions = false;
    bool list_structures = false;
    bool self_test = false;
    bool print_pdf = false;

    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--N") && i + 1 < argc) {
            nlist = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--methods") && i + 1 < argc) {
            methods_spec = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--structures") && i + 1 < argc) {
            structures_spec = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            if (!parse_int_scalar(argv[++i], trials)) {
                std::fprintf(stderr, "bad --trials value: %s\n", argv[i]);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) {
            if (!parse_int_scalar(argv[++i], threads)) {
                std::fprintf(stderr, "bad --threads value: %s\n", argv[i]);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            if (!parse_u64_scalar(argv[++i], seed)) {
                std::fprintf(stderr, "bad --seed value: %s\n", argv[i]);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--overhead") && i + 1 < argc) {
            if (!parse_int_scalar(argv[++i], overhead)) {
                std::fprintf(stderr, "bad --overhead value: %s\n", argv[i]);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--overhead-pct") && i + 1 < argc) {
            overhead_pct_spec = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--rows") && i + 1 < argc) {
            if (!parse_int_scalar(argv[++i], rows_override)) {
                std::fprintf(stderr, "bad --rows value: %s\n", argv[i]);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--N-jitter") && i + 1 < argc) {
            if (!parse_int_scalar(argv[++i], n_jitter)) {
                std::fprintf(stderr, "bad --N-jitter value: %s\n", argv[i]);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--n-jitter") && i + 1 < argc) {
            if (!parse_int_scalar(argv[++i], n_jitter)) {
                std::fprintf(stderr, "bad --n-jitter value: %s\n", argv[i]);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--matrix-seeds") && i + 1 < argc) {
            matrix_seed_mode = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--list-methods")) {
            list_methods = true;
        }
        else if (!std::strcmp(argv[i], "--list-dimensions")) {
            list_dimensions = true;
        }
        else if (!std::strcmp(argv[i], "--list-structures")) {
            list_structures = true;
        }
        else if (!std::strcmp(argv[i], "--self-test")) {
            self_test = true;
        }
        else if (!std::strcmp(argv[i], "--pdf")) {
            print_pdf = true;
        }
        else if (!std::strcmp(argv[i], "--help"))
        {
            std::printf("usage: peel_sweep [--N csv] [--methods all|legacy|combo|fast|csv] [--trials n]\n");
            std::printf("                  [--structures all|random|csv]\n");
            std::printf("                  [--overhead n|--overhead-pct csv|--rows n]\n");
            std::printf("                  [--N-jitter n]\n");
            std::printf("                  [--threads n] [--seed x]\n");
            std::printf("                  [--matrix-seeds paired|independent]\n");
            std::printf("                  [--pdf] [--self-test]\n");
            print_methods();
            print_structures();
            print_dimensions();
            return 0;
        }
        else
        {
            std::fprintf(stderr, "unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (list_methods)
    {
        print_methods();
        return 0;
    }
    if (list_dimensions)
    {
        print_dimensions();
        return 0;
    }
    if (list_structures)
    {
        print_structures();
        return 0;
    }
    if (self_test)
    {
        std::string why;
        if (!run_self_tests(&why))
        {
            std::fprintf(stderr, "self-test failed: %s\n", why.c_str());
            return 1;
        }
        std::printf("self-test: ok\n");
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
    if (n_jitter < 0)
    {
        std::fprintf(stderr, "--N-jitter must be nonnegative\n");
        return 1;
    }
    if (matrix_seed_mode != "paired" && matrix_seed_mode != "independent")
    {
        std::fprintf(stderr, "unknown --matrix-seeds mode: %s\n",
            matrix_seed_mode.c_str());
        return 1;
    }

    bool nlist_ok = false;
    bool methods_ok = false;
    bool structures_ok = false;
    bool overhead_pct_ok = false;
    const std::vector<int> Ns = parse_int_list(nlist, &nlist_ok);
    const std::vector<int> method_ids = expand_methods(methods_spec, &methods_ok);
    const std::vector<const MatrixStructure*> structures =
        expand_structures(structures_spec, &structures_ok);
    const std::vector<double> overhead_pcts =
        parse_double_list(overhead_pct_spec, &overhead_pct_ok);
    if (!nlist_ok)
    {
        std::fprintf(stderr, "invalid --N list: %s\n", nlist.c_str());
        return 1;
    }
    if (!methods_ok)
    {
        std::fprintf(stderr, "invalid --methods list: %s\n",
            methods_spec.c_str());
        return 1;
    }
    if (!structures_ok)
    {
        std::fprintf(stderr, "invalid --structures list: %s\n",
            structures_spec.c_str());
        return 1;
    }
    if (!overhead_pct_ok)
    {
        std::fprintf(stderr, "invalid --overhead-pct list: %s\n",
            overhead_pct_spec.c_str());
        return 1;
    }
    for (double pct : overhead_pcts)
    {
        if (pct < 0.0)
        {
            std::fprintf(stderr, "--overhead-pct values must be nonnegative\n");
            return 1;
        }
    }
    if (Ns.empty() || method_ids.empty() || structures.empty() ||
        overhead_pcts.empty())
    {
        std::fprintf(stderr, "empty --N, --methods, --structures, or --overhead-pct list\n");
        return 1;
    }

    print_methods();
    print_structures();
    std::printf("# sweep: trials=%d threads=%d structures=%s matrix_seeds=%s overhead=%d overhead_pct=%s rows_override=%d N_jitter=%d seed=0x%llx\n",
        trials, threads, structures_spec.c_str(), matrix_seed_mode.c_str(),
        overhead, overhead_pct_spec.c_str(), rows_override, n_jitter,
        (unsigned long long)seed);
    std::printf("method_id,method,cost,structure,N,N_jitter,"
                "actual_N_mu,actual_N_min,actual_N_max,"
                "matrix_rows_mu,matrix_rows_min,matrix_rows_max,"
                "matrix_refs_mu,row_refs_mu,row_xors_mu,matrix_xors_mu,"
                "overhead,overhead_pct,source,trials,"
                "cols_mu,cols_sd,c50,c95,c99,cmax,"
                "rows_mu,rows_sd,r50,r95,r99,rmax,"
                "comp_mu,maxcomp_mu,sumsq_mu,"
                "init_todo_mu,init_d2_mu,choices_mu,"
                "solve_sparse_xors_mu,solve_dense_xors_est_mu,solve_total_xors_est_mu,"
                "build_us,greedy_us,total_us");
    if (print_pdf) {
        std::printf(",cols_pdf,rows_pdf");
    }
    std::printf("\n");

    for (const MatrixStructure* structure : structures)
    {
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
                for (double overhead_pct : overhead_pcts)
                {
                    std::vector<TrialConfig> configs((size_t)trials);
                    for (int trial = 0; trial < trials; ++trial)
                    {
                        configs[(size_t)trial] = make_trial_config(
                            method_id, N, trial, overhead, rows_override,
                            overhead_pct, (unsigned)n_jitter, seed,
                            *structure, matrix_seed_mode);
                        if (configs[(size_t)trial].MatrixRows == 0)
                        {
                            std::fprintf(stderr, "row count is nonpositive for N=%u trial=%d\n",
                                N, trial);
                            return 1;
                        }
                        if (configs[(size_t)trial].MatrixRows > 65535u)
                        {
                            std::fprintf(stderr, "row count %u exceeds 16-bit experiment model limit\n",
                                configs[(size_t)trial].MatrixRows);
                            return 1;
                        }
                        const unsigned generated_rows =
                            generated_row_count_bound(
                                *structure,
                                configs[(size_t)trial].ActualN,
                                configs[(size_t)trial].MatrixRows);
                        if (generated_rows > 65535u)
                        {
                            std::fprintf(stderr, "generated row count %u exceeds 16-bit experiment model limit\n",
                                generated_rows);
                            return 1;
                        }
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
                                const TrialConfig config =
                                    configs[(size_t)trial];
                                results[(size_t)trial] = run_trial(
                                    method_id, config.ActualN,
                                    config.MatrixRows, *structure,
                                    config.MatrixSeed, config.ActionSeed);
                            }
                        });
                    }
                    for (std::thread& worker : workers) {
                        worker.join();
                    }

                    std::vector<unsigned> actual_ns, matrix_rows;
                    std::vector<unsigned> cols, rows, comps, maxcomps, sumsqs;
                    std::vector<unsigned> init_todo, init_d2, choices;
                    std::vector<double> matrix_refs, row_refs, row_xors;
                    std::vector<double> matrix_xors, sparse_solve_xors;
                    std::vector<double> dense_solve_xors, total_solve_xors;
                    unsigned long long build_sum = 0;
                    unsigned long long greedy_sum = 0;
                    unsigned long long total_sum = 0;
                    actual_ns.reserve(results.size());
                    matrix_rows.reserve(results.size());
                    matrix_refs.reserve(results.size());
                    row_refs.reserve(results.size());
                    row_xors.reserve(results.size());
                    matrix_xors.reserve(results.size());
                    sparse_solve_xors.reserve(results.size());
                    dense_solve_xors.reserve(results.size());
                    total_solve_xors.reserve(results.size());
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
                        actual_ns.push_back(r.ActualN);
                        matrix_rows.push_back(r.MatrixRows);
                        matrix_refs.push_back((double)r.MatrixRefs);
                        row_refs.push_back(r.MatrixRows == 0 ?
                            0.0 : (double)r.MatrixRefs / (double)r.MatrixRows);
                        row_xors.push_back(r.MatrixRows == 0 ?
                            0.0 : (double)r.MatrixXors / (double)r.MatrixRows);
                        matrix_xors.push_back((double)r.MatrixXors);
                        sparse_solve_xors.push_back((double)r.SparseSolveXors);
                        dense_solve_xors.push_back(r.DenseSolveXorsEst);
                        total_solve_xors.push_back(r.TotalSolveXorsEst);
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

                    const Summary actual_n_s = summarize(actual_ns);
                    const Summary matrix_rows_s = summarize(matrix_rows);
                    const RealSummary matrix_refs_s = summarize_real(matrix_refs);
                    const RealSummary row_refs_s = summarize_real(row_refs);
                    const RealSummary row_xors_s = summarize_real(row_xors);
                    const RealSummary matrix_xors_s = summarize_real(matrix_xors);
                    const Summary col_s = summarize(cols);
                    const Summary row_s = summarize(rows);
                    const Summary comp_s = summarize(comps);
                    const Summary maxcomp_s = summarize(maxcomps);
                    const Summary sumsq_s = summarize(sumsqs);
                    const Summary init_todo_s = summarize(init_todo);
                    const Summary init_d2_s = summarize(init_d2);
                    const Summary choice_s = summarize(choices);
                    const RealSummary sparse_solve_xors_s =
                        summarize_real(sparse_solve_xors);
                    const RealSummary dense_solve_xors_s =
                        summarize_real(dense_solve_xors);
                    const RealSummary total_solve_xors_s =
                        summarize_real(total_solve_xors);

                    const char* source_label = matrix_seed_mode == "paired" ?
                        "random-paired" : "random-independent";
                    std::printf("%d,%s,\"%s\",%s,%u,%d,"
                                "%.2f,%u,%u,"
                                "%.2f,%u,%u,"
                                "%.2f,%.4f,%.4f,%.2f,"
                                "%d,%.6g,%s,%d,"
                                "%.2f,%.2f,%u,%u,%u,%u,"
                                "%.2f,%.2f,%u,%u,%u,%u,"
                                "%.2f,%.2f,%.2f,"
                                "%.2f,%.2f,%.2f,"
                                "%.2f,%.2f,%.2f,"
                                "%.1f,%.1f,%.1f",
                        method->Id, method->Name, method->Cost, structure->Name,
                        N, n_jitter,
                        actual_n_s.Mean, actual_n_s.Min, actual_n_s.Max,
                        matrix_rows_s.Mean, matrix_rows_s.Min, matrix_rows_s.Max,
                        matrix_refs_s.Mean, row_refs_s.Mean,
                        row_xors_s.Mean, matrix_xors_s.Mean,
                        overhead, overhead_pct, source_label, trials,
                        col_s.Mean, col_s.Sd, col_s.P50, col_s.P95, col_s.P99, col_s.Max,
                        row_s.Mean, row_s.Sd, row_s.P50, row_s.P95, row_s.P99, row_s.Max,
                        comp_s.Mean, maxcomp_s.Mean, sumsq_s.Mean,
                        init_todo_s.Mean, init_d2_s.Mean, choice_s.Mean,
                        sparse_solve_xors_s.Mean,
                        dense_solve_xors_s.Mean,
                        total_solve_xors_s.Mean,
                        (double)build_sum / trials,
                        (double)greedy_sum / trials,
                        (double)total_sum / trials);
                    if (print_pdf)
                    {
                        const std::string cols_pdf = empirical_pdf(cols);
                        const std::string rows_pdf = empirical_pdf(rows);
                        std::printf(",\"%s\",\"%s\"", cols_pdf.c_str(), rows_pdf.c_str());
                    }
                    std::printf("\n");
                    std::fflush(stdout);
                }
            }
        }
    }

    return 0;
}
