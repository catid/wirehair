#include "../test/SiameseTools.h"

#include "../WirehairCodec.h"
#include "../WirehairTools.h"

#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <atomic>
using namespace std;


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
    GenerateDenseCount.cpp generates the dense counts.

    It picks random peel and dense seeds.  I found that about 500 trials
    and 10 failures is a good threshold for the "knee" of the curve for
    each N.  To make sure it found the knee I report the fourth dense
    count that is under 10.  Then I ran this a few times to produce
    the dense*.csv files in the repo.

    I approximated the shape of the curve with code in WirehairCodec.cpp
    in the ChooseSeeds() functions.
*/


//// Entrypoint

std::vector<uint8_t> message;

static std::atomic<unsigned> FailedTrials(0);

static void RandomTrial(
    unsigned N,
    unsigned count,
    uint64_t seed,
    unsigned trial)
{
    const uint16_t dense_count = (uint16_t)count;

    siamese::PCGRandom prng;
    prng.Seed(seed, trial);

    const uint16_t p_seed = (uint16_t)prng.Next();
    const uint16_t d_seed = (uint16_t)prng.Next();

    wirehair::Codec codec;

    // Override the seeds
    codec.OverrideSeeds(dense_count, p_seed, d_seed);

    // Initialize codec
    WirehairResult result = codec.InitializeEncoder(N, 1);

    // If initialization succeeded:
    if (result == Wirehair_Success) {
        // Feed message to codec
        result = codec.EncodeFeed(&message[0]);
    }

    if (result != Wirehair_Success) {
        FailedTrials++;
    }
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
    uint64_t value = (base_seed ^ 0x44434f554e54ULL) +
        0x9e3779b97f4a7c15ULL * (uint64_t)N;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

static uint64_t SourceTablesHash()
{
    uint64_t hash = 1469598103934665603ULL;
    for (unsigned N = CAT_WIREHAIR_MIN_N;
         N <= CAT_WIREHAIR_MAX_N;
         ++N)
    {
        hash ^= N;
        hash *= 1099511628211ULL;
        hash ^= wirehair::GetDenseCount(N);
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

static void Usage(const char* program)
{
    cerr
        << "usage: " << program << " [--seed U64] [--trials N] "
        << "[--nlo N] [--nhi N] [--max-failures N] "
        << "[--low-count-run N] [--results FILE]\n";
}

int main(int argc, char** argv)
{
    uint64_t seed = siamese::GetTimeUsec();
    unsigned trials = 500;
    unsigned nlo = 2;
    unsigned nhi = 64000;
    unsigned maxFailures = 10;
    unsigned lowCountRun = 4;
    string resultsPath;
    string resultsTemporaryPath;

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
        else if (!strcmp(arg, "--max-failures")) {
            if (!ParseUnsigned(next(), maxFailures)) {
                cerr << "bad --max-failures value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--low-count-run")) {
            if (!ParseUnsigned(next(), lowCountRun) || lowCountRun == 0) {
                cerr << "bad --low-count-run value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--results")) {
            resultsPath = next();
            if (resultsPath.empty()) {
                cerr << "bad --results value" << endl;
                return 1;
            }
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

    if (nlo < CAT_WIREHAIR_MIN_N || nhi > CAT_WIREHAIR_MAX_N || nlo > nhi) {
        cerr << "N range must be in [" << CAT_WIREHAIR_MIN_N << ", "
            << CAT_WIREHAIR_MAX_N << "]" << endl;
        return 1;
    }

    const int gfInitResult = gf256_init();

    // If gf256 init failed:
    if (gfInitResult != 0)
    {
        cout << "GF256 init failed" << endl;
        return -1;
    }

    message.resize(64000);
    const uint64_t sourceTablesHash = SourceTablesHash();

    ofstream results;
    if (!resultsPath.empty())
    {
        resultsTemporaryPath = resultsPath + ".tmp";
        if (PathExists(resultsPath) || PathExists(resultsTemporaryPath)) {
            cerr << "--results destination already exists: " << resultsPath
                << endl;
            return 1;
        }
        results.open(resultsTemporaryPath.c_str(), ios::out | ios::trunc);
        if (!results) {
            cerr << "unable to open --results file: " << resultsPath << endl;
            return 1;
        }
        results << "N\tDenseCount\tSelectedFailures\tTrials"
            "\tMaxFailures\tLowCountRun\tBaseSeed\tTrialSeed"
            "\tQualifyingCount\tSourceTablesHash\tMethodVersion\tStatus\n";
    }

    cerr
        << "# GenerateDenseCount seed=0x" << hex << seed << dec
        << " trials=" << trials
        << " nlo=" << nlo
        << " nhi=" << nhi
        << " max_failures=" << maxFailures
        << " low_count_run=" << lowCountRun << endl;

    cout << "N\tDenseCount\tSelectedFailureRate" << endl;

    // Always walk the canonical N=2 lattice and filter by the requested
    // range.  Combined with per-N seeds and search starts, concatenated range
    // shards are byte-equivalent to a monolithic full-domain run.
    unsigned emittedRows = 0;
    for (unsigned N = CAT_WIREHAIR_MIN_N; N <= nhi;)
    {
        if (N < nlo)
        {
            unsigned scale = N / 64;
            if (scale < 1) {
                scale = 1;
            }
            const unsigned next = N + scale;
            N = next > CAT_WIREHAIR_MAX_N ? CAT_WIREHAIR_MAX_N : next;
            continue;
        }

        unsigned lowCount = 0;

        // Use a per-N start rather than carrying the prior result, so shard
        // boundaries cannot alter the candidates evaluated for an N.
        int count = (int)wirehair::GetDenseCount(N) - 10;
        if (count < 1) {
            count = 1;
        }

        unsigned lowestFailures = UINT_MAX;
        unsigned lowestFailureCount = (unsigned)count;
        unsigned selectedFailures = UINT_MAX;
        unsigned selectedCount = 0;
        const uint64_t trialSeed = DeriveTrialSeed(seed, N);

        for (; count <= (int)N; ++count)
        {
#if 0
            if (count % 4 != 2 && N > 32) {
                continue;
            }
#endif

            FailedTrials = 0;

#pragma omp parallel for
            for (int trial = 0; trial < (int)trials; ++trial) {
                RandomTrial(N, count, trialSeed, (unsigned)trial);
            }

            const unsigned failures = FailedTrials;
            //cout << " *** " << N << "\t" << count << "\t" << (failures / (float)trials) << endl;

            if (failures < lowestFailures)
            {
                lowestFailures = failures;
                lowestFailureCount = (unsigned)count;
            }

            if (failures <= maxFailures)
            {
                // Select the configured Nth count under the failure ceiling,
                // matching the methodology documented above.  The old code
                // stopped here but reported the unrelated minimum-failure
                // count, which did not implement its stated knee rule.
                if (++lowCount >= lowCountRun) {
                    selectedCount = (unsigned)count;
                    selectedFailures = failures;
                    break;
                }
            }
        }

        // A pathological N with too few qualifying counts still gets the
        // best observed fallback rather than an invalid zero entry.
        if (selectedCount == 0)
        {
            selectedCount = lowestFailureCount;
            selectedFailures = lowestFailures;
        }

        cout << N << "\t" << selectedCount << "\t" <<
            (selectedFailures / (float)trials) << endl;
        if (results)
        {
            results << N << '\t' << selectedCount << '\t'
                << selectedFailures << '\t' << trials << '\t'
                << maxFailures << '\t' << lowCountRun << '\t' << seed
                << '\t' << trialSeed << '\t' << lowCount << '\t'
                << sourceTablesHash << "\tdense-count-v2\t"
                << (lowCount >= lowCountRun ? "selected-knee" : "fallback-minimum")
                << '\n';
            if (!results) {
                cerr << "failed writing --results file: " << resultsPath << endl;
                return 1;
            }
        }
        ++emittedRows;

        if (N == CAT_WIREHAIR_MAX_N) {
            break;
        }
        int scale = N / 64;
        if (scale < 1) {
            scale = 1;
        }
        const unsigned next = N + (unsigned)scale;
        N = next > CAT_WIREHAIR_MAX_N ? CAT_WIREHAIR_MAX_N : next;
    }

    if (emittedRows == 0) {
        cerr << "requested N range selects no canonical dense-count samples"
            << endl;
        return 1;
    }

    return FinalizeResults(
        results, resultsTemporaryPath, resultsPath) ? 0 : 1;
}
