#include "../test/SiameseTools.h"

#include "../WirehairCodec.h"
#include "../WirehairTools.h"
#include "SeedSelectionPolicy.h"

#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <atomic>
using namespace std;

//#define ENABLE_FULL_SEARCH


#if !defined(ENABLE_FULL_SEARCH)
// This allows me to run the Unit Test to evaluate seeds, and then
// the ones that tend to fail too much can be put in this list and
// refined further.
static const int N_List[] = {
    17,
};
// This only works up to 2047.  After that the GenerateMostDenseSeeds
// and GeneratePeelSeeds programs take over.
#endif // !ENABLE_FULL_SEARCH


/**
    Wirehair seed selection methodology

    Overview:

    Wirehair's solver inverts a hybrid matrix with three types of rows:
    (1) Peel rows - Sparse binary mixing of recovery set.
    (2) Dense rows - Dense binary mixing of recovery set.
    (3) Heavy rows - Small GF(256) Cauchy matrix in the lower right.

    The encoder inverts the matrix to find the recovery set,
    and then it generates new random peel rows for the decoder.

    The decoder inverts the matrix to find the recovery set,
    and then it regenerates the original data (first N peel rows).


    Heavy Rows:

    The fixed 6x18 GF(256) matrix improves rank recovery but has an empirical
    singular floor after binary-row mixing.  Its bytes and limits are
    documented in HEAVY_MATRIX.md and reproduced by HeavyRowGenerator.cpp.


    Dense Rows:

    We have two knobs to tune for these:
    (1) The seed to generate the random rows.
    (2) The number of rows to generate.

    The seed seems somewhat arbitrary.  Some will be better than others,
    but on average I expected them to be pretty similar.  The number of
    rows to generate seemed more important.  I expected that increasing
    the number of dense rows past a certain point would stop helping.
    Also the number of dense rows needed would increase as N increases,
    and should not decrease.  This means that once the number of dense
    rows is determined, then the seed can be optimized further at the
    largest N that uses it, which should also help all the smaller N
    that use the same row count.  So we can use a small list of seeds
    for each of the dense row counts.

    For benchmarking different dense row counts, we would want to pick
    a bunch of random seeds and check invertibility rates.  There should
    be a knee where increasing the row count doesn't help.


    Peel Rows:

    These are for the most part arbitrary.  Optimizing them to be more
    robust to a small number of losses seems worthwhile.
*/

/**
    GenerateSmallDenseSeeds.cpp generates the dense seeds and counts for
    the smallest values of N.

    It uses many random peel seeds, tries all dense seeds, and many dense counts.
*/


//// Entrypoint

uint8_t Message[64000];

static std::atomic<unsigned> FailedTrials(0);

static void RandomTrial(
    unsigned N,
    unsigned count,
    uint64_t seed,
    const uint16_t d_seed,
    unsigned trial)
{
    const uint16_t dense_count = (uint16_t)count;

    siamese::PCGRandom prng;
    prng.Seed(seed, trial);

    const uint16_t p_seed = (uint16_t)(prng.Next() ^ prng.Next());

    wirehair::Codec codec;

    // Override the seeds
    codec.OverrideSeeds(dense_count, p_seed, d_seed);

    // Initialize codec
    WirehairResult result = codec.InitializeEncoder(N, 1);

    // If initialization succeeded:
    if (result == Wirehair_Success) {
        // Feed message to codec
        result = codec.EncodeFeed(&Message[0]);
    }

    if (result != Wirehair_Success) {
        FailedTrials++;
    }
}

static void RandomPeelLoss(
    unsigned N,
    const uint16_t dense_count,
    const uint16_t p_seed,
    const uint16_t d_seed,
    unsigned miss1,
    unsigned miss2)
{
    wirehair::Codec encoder, decoder;

    // Override the seeds
    encoder.OverrideSeeds(dense_count, p_seed, d_seed);
    decoder.OverrideSeeds(dense_count, p_seed, d_seed);

    // Initialize codec
    WirehairResult result = encoder.InitializeEncoder(N, 1);

    // If initialization succeeded:
    if (result == Wirehair_Success) {
        // Feed message to codec
        result = encoder.EncodeFeed(&Message[0]);
    }

    if (result != Wirehair_Success) {
        // Huge penalty
        FailedTrials += 10000;
        return;
    }

    result = decoder.InitializeDecoder(N, 1);

    if (result != Wirehair_Success) {
        cout << "InitializeDecoder failed" << endl;
        // Huge penalty
        FailedTrials += 10000;
        return;
    }

    unsigned added = 0;

    unsigned failures = 0;

    for (unsigned i = 0;; ++i)
    {
        if (i == miss1 || i == miss2) {
            continue;
        }

        ++added;

        uint8_t encodedData[1];

        const uint32_t encodedBytes = encoder.Encode(i, encodedData, 1);

        if (encodedBytes == 0)
        {
            cout << "Encode failed" << endl;
            // Huge penalty
            FailedTrials += 10000;
            return;
        }

        result = decoder.DecodeFeed(i, encodedData, encodedBytes);

        if (result != Wirehair_NeedMore)
        {
            if (result == Wirehair_Success)
            {
                if (added == N) {
                    // No penalty
                }
                else if (added == N + 1) {
                    failures = 1;
                }
                else if (added == N + 2) {
                    failures = 3;
                }
                else if (added >= N + 3) {
                    failures = (added - N) * 3;
                }
                break;
            }

            cout << "DecodeFeed failed" << endl;
            // Huge penalty
            FailedTrials += 10000;
            return;
        }

        if (i > N + 10) {
            // Huge penalty
            FailedTrials += 10000;
            return;
        }
    }

    FailedTrials += failures;

#if 0
    if (failures > 0)
    {
        cout << failures << " failures for " << miss1 << " and " << miss2 << endl;
    }
#endif

    uint8_t DecodeMessage[64000];

    WirehairResult decodeResult = decoder.ReconstructOutput(DecodeMessage, N);

    if (decodeResult != Wirehair_Success)
    {
        cout << "ReconstructOutput failed" << endl;
        // Huge penalty
        FailedTrials += 10000;
        return;
    }

    if (0 != memcmp(DecodeMessage, Message, N))
    {
        cout << "memcmp failed" << endl;
        // Huge penalty
        FailedTrials += 10000;
        return;
    }

    // Success
}

static const unsigned kTinyTableCount = wirehair::kTinyTableCount;
static const unsigned kSmallTableCount = wirehair::kSmallTableCount;

static uint8_t kTinyDenseCounts[kTinyTableCount] = {
};

static uint16_t kTinyDenseSeeds[kTinyTableCount] = {
};

static uint8_t kSmallDenseSeeds[kSmallTableCount] = {
};

static uint8_t kSmallPeelSeeds[kTinyTableCount + kSmallTableCount] = {
};

void FillTables()
{
    memcpy(kTinyDenseCounts, wirehair::kTinyDenseCounts, sizeof(kTinyDenseCounts));
    memcpy(kTinyDenseSeeds, wirehair::kTinyDenseSeeds, sizeof(kTinyDenseSeeds));
    memcpy(kSmallDenseSeeds, wirehair::kSmallDenseSeeds, sizeof(kSmallDenseSeeds));
    memcpy(kSmallPeelSeeds, wirehair::kSmallPeelSeeds, sizeof(kSmallPeelSeeds));
}

static uint16_t LinearInterpolate(
    unsigned N0,
    unsigned N1,
    unsigned Count0,
    unsigned Count1,
    unsigned N)
{
    CAT_DEBUG_ASSERT(N >= N0 && N <= N1);
    CAT_DEBUG_ASSERT(Count1 > Count0);

    const unsigned numerator = (N - N0) * (Count1 - Count0);
    const unsigned denominator = N1 - N0;

    const unsigned count = Count0 + (unsigned)(numerator / denominator);

    return static_cast<uint16_t>(count);
}

static const unsigned kTinyCountsFromGraph[65] = {
    0,
    0,
    2,
    3,
    3,
    5,
    6,
    7,
    7,
    8,
    10,
    11,
    11,
    11,
    14,
    14,
    15,
    13,
    13,
    14,
    17,
    16,
    22,
    14,
    16,
    14,
    20,
    23,
    24,
    22,
    15,
    22,
    15,
    24,
    15,
    19,
    17,
    26,
    31,
    27,
    22,
    18,
    16,
    15,
    23,
    17,
    22,
    19,
    14,
    18,
    22,
    18,
    19,
    20,
    24,
    19,
    18,
    20,
    20,
    19,
    25,
    22,
    21,
    25,
    17
};

static unsigned GetDenseCountGuess(unsigned N)
{
    if (N <= 64) {
        return kTinyCountsFromGraph[N];
    }

    unsigned dense_count = 0;

    if (N <= 500) {
        dense_count = LinearInterpolate(64, 500, 26, 35, N);
    }
    else if (N <= 1000) {
        dense_count = LinearInterpolate(500, 1000, 35, 48, N);
    }
    else if (N <= 2048) {
        dense_count = LinearInterpolate(1000, 2048, 48, 62, N);
    }

    // Round up to the next D s.t. D Mod 4 = 2
    switch (dense_count & 3)
    {
    case 0: dense_count += 2; break;
    case 1: dense_count += 1; break;
    case 2: break;
    case 3: dense_count += 3; break;
    }

    return dense_count;
}

static bool ParseU64(const char* text, uint64_t& out)
{
    if (!text || !*text || *text < '0' || *text > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long value = strtoull(text, &end, 0);
    if (errno != 0 || !end || *end != '\0') {
        return false;
    }
    out = (uint64_t)value;
    return true;
}

static bool ParseUnsigned(const char* text, unsigned& out)
{
    uint64_t value = 0;
    if (!ParseU64(text, value) || value > UINT_MAX) {
        return false;
    }
    out = (unsigned)value;
    return true;
}

static uint64_t DeriveTrialSeed(uint64_t base_seed, unsigned N)
{
    // SplitMix64 finalizer: stable across processes, compilers, and shard plans.
    uint64_t value = (base_seed ^ 0x534d414c4cULL) +
        0x9e3779b97f4a7c15ULL * (uint64_t)N;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

static uint64_t SmallTablesHash()
{
    uint64_t hash = 1469598103934665603ULL;
    for (unsigned N = 2; N < kTinyTableCount + kSmallTableCount; ++N)
    {
        const uint16_t count = wirehair::GetDenseCount(N);
        const uint16_t denseSeed = wirehair::GetBaseDenseSeed(N, count);
        const uint16_t peelSeed = wirehair::GetBasePeelSeed(N);
        hash ^= N;
        hash *= 1099511628211ULL;
        hash ^= count;
        hash *= 1099511628211ULL;
        hash ^= denseSeed;
        hash *= 1099511628211ULL;
        hash ^= peelSeed;
        hash *= 1099511628211ULL;
    }
    return hash;
}

static bool PathExists(const string& path)
{
    ifstream input(path.c_str(), ios::binary);
    return input.good();
}

static bool FinalizeResults(
    ofstream& results,
    const string& temporary_path,
    const string& final_path)
{
    if (final_path.empty()) {
        return true;
    }
    results.flush();
    const bool write_ok = results.good();
    results.close();
    if (!write_ok || !results) {
        cerr << "failed writing --results file: " << temporary_path << endl;
        return false;
    }
    if (rename(temporary_path.c_str(), final_path.c_str()) != 0) {
        cerr << "unable to publish --results file " << final_path << ": "
            << strerror(errno) << endl;
        return false;
    }
    return true;
}

static bool ParseNList(
    const char* text,
    vector<unsigned>& values,
    string& error)
{
    values.clear();
    const string input = text ? text : "";
    if (input.empty() || input.back() == ',')
    {
        error = "--nlist must not be empty or end with a comma";
        return false;
    }
    size_t begin = 0;
    while (begin < input.size())
    {
        const size_t comma = input.find(',', begin);
        const string item = input.substr(
            begin,
            comma == string::npos ? string::npos : comma - begin);
        unsigned value = 0;
        if (!ParseUnsigned(item.c_str(), value))
        {
            error = "--nlist must contain comma-separated unsigned integers";
            return false;
        }
        for (unsigned previous : values)
        {
            if (previous == value)
            {
                error = "--nlist contains duplicate N=" + to_string(value);
                return false;
            }
        }
        values.push_back(value);
        if (comma == string::npos) {
            break;
        }
        begin = comma + 1u;
    }
    return true;
}

static bool PlanSelection(
    unsigned nlo,
    unsigned nhi,
    const vector<unsigned>& candidates,
    vector<unsigned>& selected,
    string& error)
{
    const unsigned maxN = kTinyTableCount + kSmallTableCount - 1u;
    selected.clear();
    if (nlo < 2u || nhi > maxN || nlo > nhi)
    {
        error = "N range must be in [2, " + to_string(maxN) + "]";
        return false;
    }
    for (size_t i = 0; i < candidates.size(); ++i)
    {
        const unsigned value = candidates[i];
        if (value < 2u || value > maxN)
        {
            error = "candidate N=" + to_string(value) + " is outside the table domain";
            return false;
        }
        for (size_t j = 0; j < i; ++j)
        {
            if (candidates[j] == value)
            {
                error = "candidate list contains duplicate N=" + to_string(value);
                return false;
            }
        }
        if (value >= nlo && value <= nhi) {
            selected.push_back(value);
        }
    }
    if (selected.empty())
    {
        error = "requested N range selects no configured candidates";
        return false;
    }
    return true;
}

static bool SelectionSelfTest()
{
    const vector<unsigned> configured = { 17, 19, 21 };
    vector<unsigned> selected;
    string error;
    if (!PlanSelection(17, 17, configured, selected, error) ||
        selected != vector<unsigned>{ 17 }) {
        return false;
    }
    if (!PlanSelection(18, 20, configured, selected, error) ||
        selected != vector<unsigned>{ 19 }) {
        return false;
    }
    if (PlanSelection(18, 18, configured, selected, error) ||
        PlanSelection(20, 19, configured, selected, error) ||
        PlanSelection(1, 17, configured, selected, error)) {
        return false;
    }

    vector<unsigned> explicit_values;
    if (!ParseNList("17,18,20", explicit_values, error) ||
        explicit_values != vector<unsigned>({ 17, 18, 20 }) ||
        !PlanSelection(18, 20, explicit_values, selected, error) ||
        selected != vector<unsigned>({ 18, 20 })) {
        return false;
    }
    return !ParseNList("", explicit_values, error) &&
        !ParseNList("17,17", explicit_values, error) &&
        !ParseNList("17,,18", explicit_values, error) &&
        !ParseNList("17,", explicit_values, error);
}

static void Usage(const char* program)
{
    cerr
        << "usage: " << program << " [--seed U64] [--trials N] "
        << "[--nlo N] [--nhi N] [--nlist N[,N...]] "
        << "[--full-search] [--results FILE] [--selection-self-test]\n";
}

int main(int argc, char** argv)
{
    uint64_t seed = siamese::GetTimeUsec();
    unsigned trials = 3500;
    unsigned nlo = 2;
    unsigned nhi = kTinyTableCount + kSmallTableCount - 1u;
    vector<unsigned> explicit_candidates;
    bool has_explicit_candidates = false;
    bool full_search = false;
    bool selection_self_test = false;
    string results_path;
    string results_temporary_path;

    for (int i = 1; i < argc; ++i)
    {
        const char* arg = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) {
                cerr << "missing value for " << arg << endl;
                Usage(argv[0]);
                exit(1);
            }
            return argv[++i];
        };
        if (!strcmp(arg, "--seed")) {
            if (!ParseU64(next(), seed)) {
                cerr << "bad --seed value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--trials")) {
            if (!ParseUnsigned(next(), trials) || trials == 0 ||
                trials > (unsigned)INT_MAX)
            {
                cerr << "bad --trials value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--nlo")) {
            if (!ParseUnsigned(next(), nlo)) {
                cerr << "bad --nlo value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--nhi")) {
            if (!ParseUnsigned(next(), nhi)) {
                cerr << "bad --nhi value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--nlist")) {
            string error;
            if (!ParseNList(next(), explicit_candidates, error)) {
                cerr << error << endl;
                return 1;
            }
            has_explicit_candidates = true;
        }
        else if (!strcmp(arg, "--full-search")) {
            full_search = true;
        }
        else if (!strcmp(arg, "--results")) {
            results_path = next();
            if (results_path.empty()) {
                cerr << "bad --results value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--selection-self-test")) {
            selection_self_test = true;
        }
        else if (!strcmp(arg, "--help")) {
            Usage(argv[0]);
            return 0;
        }
        else {
            cerr << "unknown argument " << arg << endl;
            Usage(argv[0]);
            return 1;
        }
    }

    if (selection_self_test)
    {
        if (!SelectionSelfTest()) {
            cerr << "selection planner self-test failed" << endl;
            return 1;
        }
        cout << "selection planner self-test passed" << endl;
        return 0;
    }

    if (full_search && has_explicit_candidates) {
        cerr << "--full-search cannot be combined with --nlist" << endl;
        return 1;
    }

    vector<unsigned> candidates;
    const char* selection_source = nullptr;
    if (has_explicit_candidates) {
        candidates = explicit_candidates;
        selection_source = "explicit";
    }
    else if (full_search) {
        const unsigned maxN = kTinyTableCount + kSmallTableCount - 1u;
        for (unsigned N = 2; N <= maxN; ++N) {
            candidates.push_back(N);
        }
        selection_source = "full";
    }
#ifdef ENABLE_FULL_SEARCH
    else {
        const unsigned maxN = kTinyTableCount + kSmallTableCount - 1u;
        for (unsigned N = 2; N <= maxN; ++N) {
            candidates.push_back(N);
        }
        selection_source = "full";
    }
#else
    else {
        const size_t count = sizeof(N_List) / sizeof(N_List[0]);
        for (size_t i = 0; i < count; ++i) {
            candidates.push_back((unsigned)N_List[i]);
        }
        selection_source = "configured";
    }
#endif

    vector<unsigned> selected;
    string selection_error;
    if (!PlanSelection(nlo, nhi, candidates, selected, selection_error)) {
        cerr << selection_error << endl;
        return 1;
    }

    FillTables();
    const uint64_t smallTablesHash = SmallTablesHash();

    ofstream results;
    if (!results_path.empty())
    {
        results_temporary_path = results_path + ".tmp";
        if (PathExists(results_path) || PathExists(results_temporary_path)) {
            cerr << "--results destination already exists: " << results_path
                << endl;
            return 1;
        }
        results.open(results_temporary_path.c_str(), ios::out | ios::trunc);
        if (!results) {
            cerr << "unable to open --results file: " << results_path << endl;
            return 1;
        }
        results << "N\tDenseCount\tDenseSeed\tPeelSeed\tDenseFailures"
            "\tPeelScore\tCurrentDenseCount\tCurrentDenseSeed"
            "\tCurrentPeelSeed\tTrials\tBaseSeed\tTrialSeed"
            "\tSmallTablesHash\tMethodVersion\n";
    }

    const int gfInitResult = gf256_init();

    // If gf256 init failed:
    if (gfInitResult != 0)
    {
        cout << "GF256 init failed" << endl;
        return -1;
    }

    for (unsigned i = 0; i < 64000; ++i) {
        Message[i] = (uint8_t)i;
    }

    cerr
        << "# GenerateSmallDenseSeeds seed=0x" << hex << seed << dec
        << " trials=" << trials
        << " nlo=" << nlo
        << " nhi=" << nhi
        << " selected=" << selected.size()
        << " source=" << selection_source
        << endl;

    unsigned processed = 0;
    for (unsigned selected_n : selected)
    {
        const int N = (int)selected_n;
        const unsigned currentDenseCount = wirehair::GetDenseCount(N);
        const unsigned currentDenseSeed = wirehair::GetBaseDenseSeed(
            N, currentDenseCount);
        const unsigned currentPeelSeed = wirehair::GetBasePeelSeed(N);
        const uint64_t trialSeed = DeriveTrialSeed(seed, (unsigned)N);

        int countGuess = (int)GetDenseCountGuess(N);
        int countGuessMin = countGuess - 2;
        if (countGuessMin < 2) {
            countGuessMin = 2;
        }
        int countGuessMax = countGuess + 2;
        if (countGuessMax > N) {
            countGuessMax = N;
        }

        unsigned d_seed_max = 65536;

        if (N > 64)
        {
            countGuessMin = countGuess;
            countGuessMax = countGuess;
            d_seed_max = 256;
        }

        int best_seed = (int)currentDenseSeed;
        int best_count = (int)currentDenseCount;
        int best_failures = INT_MAX;
        vector<unsigned> countCandidates;
        if (currentDenseCount >= 2u && currentDenseCount <= (unsigned)N) {
            countCandidates.push_back(currentDenseCount);
        }
        for (int count = countGuessMin; count <= countGuessMax; ++count) {
            if ((unsigned)count != currentDenseCount) {
                countCandidates.push_back((unsigned)count);
            }
        }

        bool foundZeroDenseScore = false;
        for (unsigned count : countCandidates)
        {
            cout << "For N = " << N << " trying count = " << count << endl;

            for (unsigned rank = 0; rank < d_seed_max; ++rank)
            {
                const unsigned d_seed = wirehair_table::CurrentFirstCandidate(
                    rank, currentDenseSeed, d_seed_max);
                FailedTrials = 0;

#pragma omp parallel for
                for (int trial = 0; trial < (int)trials; ++trial) {
                    RandomTrial(
                        N, (int)count, trialSeed, (uint16_t)d_seed, trial);
                }

                const int failures = FailedTrials;

                if (failures < best_failures)
                {
                    best_seed = d_seed;
                    best_count = (int)count;
                    best_failures = failures;

                    cout << "*** New best dense: N = " << N << ": Best dense seed = " << best_seed << ", best count = " << best_count << ", best failures = " << best_failures << endl;

                    if (failures <= 0) {
                        foundZeroDenseScore = true;
                        break;
                    }
                }
            }
            if (foundZeroDenseScore) {
                break;
            }
        }

        cout << "N = " << N << ": Best dense seed = " << best_seed << ", best count = " << best_count << ", best failures = " << best_failures << endl;


        int bestPeelSeed = (int)currentPeelSeed;
        int bestPeelFails = INT_MAX;

        // kSmallPeelSeeds stores uint8_t entries, so searching past 255 would
        // pick a winner that truncates to a different (already-rejected) seed
        // at runtime.  Keep the search inside the table's representable range.
        const int peelSeedMax = 256;

        for (unsigned rank = 0; rank < (unsigned)peelSeedMax; ++rank)
        {
            const int peelSeed = (int)wirehair_table::CurrentFirstCandidate(
                rank, currentPeelSeed, (unsigned)peelSeedMax);
            FailedTrials = 0;

            if (N < 300)
            {
#pragma omp parallel for
                for (int miss1 = 0; miss1 < N; ++miss1) {
                    for (int miss2 = miss1; miss2 < N; ++miss2) {
                        RandomPeelLoss(N, (uint16_t)best_count, (uint16_t)peelSeed, (uint16_t)best_seed, miss1, miss2);
                    }
                }
            }
            else
            {
#pragma omp parallel for
                for (int miss1 = 1; miss1 < N; ++miss1) {
                    RandomPeelLoss(N, (uint16_t)best_count, (uint16_t)peelSeed, (uint16_t)best_seed, miss1 - 1, miss1);
                }
#pragma omp parallel for
                for (int miss1 = 0; miss1 < N; ++miss1) {
                    RandomPeelLoss(N, (uint16_t)best_count, (uint16_t)peelSeed, (uint16_t)best_seed, miss1, miss1);
                }
            }

            const int failures = FailedTrials;

            if (failures < bestPeelFails)
            {
                bestPeelFails = failures;
                bestPeelSeed = peelSeed;

                cout << "*** New best peel: N = " << N << ": Best peel seed = " << bestPeelSeed << " fails = " << bestPeelFails << endl;

                if (failures <= 0) {
                    break;
                }
            }
        }

        cout << "N = " << N << ": Best peel seed = " << bestPeelSeed << " fails = " << bestPeelFails << endl;

        if ((unsigned)N < kTinyTableCount) {
            kTinyDenseCounts[N] = (uint8_t)best_count;
            kTinyDenseSeeds[N] = (uint16_t)best_seed;
        }
        else {
            kSmallDenseSeeds[N - kTinyTableCount] = (uint8_t)best_seed;
        }
        kSmallPeelSeeds[N] = (uint8_t)bestPeelSeed;
        if (!results_path.empty())
        {
            results << N << '\t' << best_count << '\t' << best_seed << '\t'
                << bestPeelSeed << '\t' << best_failures << '\t'
                << bestPeelFails << '\t' << currentDenseCount << '\t'
                << currentDenseSeed << '\t' << currentPeelSeed << '\t'
                << trials << '\t' << seed << '\t' << trialSeed << '\t'
                << smallTablesHash << "\tsmall-v3\n";
            if (!results) {
                cerr << "failed writing --results file: " << results_path << endl;
                return 1;
            }
        }
        ++processed;
    }

    cerr << "# processed_candidates=" << processed << endl;

    cout << "static const unsigned kTinyTableCount = " << kTinyTableCount << ";" << endl;
    cout << "static const unsigned kSmallTableCount = " << kSmallTableCount << ";" << endl;

    const unsigned modulus = 32;

    cout << "const uint8_t kTinyDenseCounts[kTinyTableCount] = {" << endl;

    for (unsigned i = 0; i < kTinyTableCount; ++i)
    {
        if (i % modulus == 0) {
            cout << "    ";
        }

        cout << (int)kTinyDenseCounts[i] << ",";

        if ((i + 1) % modulus == 0) {
            cout << endl;
        }
    }

    cout << endl << "};" << endl << endl;

    cout << "const uint16_t kTinyDenseSeeds[kTinyTableCount] = {" << endl;

    for (unsigned i = 0; i < kTinyTableCount; ++i)
    {
        if (i % modulus == 0) {
            cout << "    ";
        }

        cout << (int)kTinyDenseSeeds[i] << ",";

        if ((i + 1) % modulus == 0) {
            cout << endl;
        }
    }

    cout << endl << "};" << endl;

    cout << endl << "// This table skips the first kTinyTableCount elements" << endl;
    cout << "const uint8_t kSmallDenseSeeds[kSmallTableCount] = {" << endl;

    for (unsigned i = 0; i < kSmallTableCount; ++i)
    {
        if (i % modulus == 0) {
            cout << "    ";
        }

        cout << (int)kSmallDenseSeeds[i] << ",";

        if ((i + 1) % modulus == 0) {
            cout << endl;
        }
    }

    cout << endl << "};" << endl << endl;

    cout << "const uint8_t kSmallPeelSeeds[kTinyTableCount + kSmallTableCount] = {" << endl;

    for (unsigned i = 0; i < kTinyTableCount + kSmallTableCount; ++i)
    {
        if (i % modulus == 0) {
            cout << "    ";
        }

        cout << (int)kSmallPeelSeeds[i] << ",";

        if ((i + 1) % modulus == 0) {
            cout << endl;
        }
    }

    cout << "};" << endl;

    return FinalizeResults(
        results, results_temporary_path, results_path) ? 0 : 1;
}
