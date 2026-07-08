#include "../test/SiameseTools.h"

#include "../WirehairCodec.h"

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
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

    The heavy rows make up a special 6x18 Cauchy matrix with the property
    that any disturbance from the binary rows above does not affect the
    100% invertibility of the matrix.  Since these rows are not affected
    by any of the rows above, we just pick one of these special matrices
    and use it for all values of N original blocks.

    The code to generate this is in HeavyRowGenerator.cpp


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

static void Usage(const char* program)
{
    cerr
        << "usage: " << program << " [--seed U64] [--trials N] "
        << "[--nlo N] [--nhi N]\n";
}

int main(int argc, char** argv)
{
    uint64_t seed = siamese::GetTimeUsec();
    unsigned trials = 3500;
    unsigned nlo = 2;
    unsigned nhi = kTinyTableCount + kSmallTableCount - 1u;

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

    const unsigned maxN = kTinyTableCount + kSmallTableCount - 1u;
    if (nlo < 2u || nhi > maxN || nlo > nhi) {
        cerr << "N range must be in [2, " << maxN << "]" << endl;
        return 1;
    }

    FillTables();

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
        << " nhi=" << nhi << endl;

#ifdef ENABLE_FULL_SEARCH
    static const int N_Min = 2;
    static const int N_Max = kTinyTableCount + kSmallTableCount - 1;
    for (int N = N_Min; N <= N_Max; ++N)
    {
        if ((unsigned)N < nlo || (unsigned)N > nhi) {
            continue;
        }
#else
    // This allows me to run the Unit Test to evaluate seeds, and then
    // the ones that tend to fail too much can be put in this list and
    // refined further.
    const int nListCount = (int)(sizeof(N_List) / sizeof(N_List[0]));
    for (int N_i = 0; N_i < nListCount; ++N_i)
    {
        int N = N_List[N_i];
        if ((unsigned)N < nlo || (unsigned)N > nhi) {
            continue;
        }
#endif

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

        int best_seed = 0;
        int best_count = 0;
        int best_failures = 10000;

        for (int count = countGuessMin; count <= countGuessMax; ++count)
        {
            cout << "For N = " << N << " trying count = " << count << endl;

            // Advance the trial-stream seed once per count, not once per
            // candidate: candidates must share paired trials or selection by
            // minimum failure count is biased (winner's curse)
            ++seed;

            for (unsigned d_seed = 0; d_seed < d_seed_max; ++d_seed)
            {
                FailedTrials = 0;

#pragma omp parallel for
                for (int trial = 0; trial < (int)trials; ++trial) {
                    RandomTrial(N, count, seed, (uint16_t)d_seed, trial);
                }

                const int failures = FailedTrials;

                if (failures < best_failures)
                {
                    best_seed = d_seed;
                    best_count = count;
                    best_failures = failures;

                    cout << "*** New best dense: N = " << N << ": Best dense seed = " << best_seed << ", best count = " << best_count << ", best failures = " << best_failures << endl;

                    if (failures <= 0) {
                        break;
                    }
                }
            }
        }

        cout << "N = " << N << ": Best dense seed = " << best_seed << ", best count = " << best_count << ", best failures = " << best_failures << endl;


        int bestPeelSeed = 0;
        int bestPeelFails = 100000;

        // kSmallPeelSeeds stores uint8_t entries, so searching past 255 would
        // pick a winner that truncates to a different (already-rejected) seed
        // at runtime.  Keep the search inside the table's representable range.
        const int peelSeedMax = 256;

        for (int peelSeed = 0; peelSeed < peelSeedMax; ++peelSeed)
        {
            ++seed;

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
    }

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

    return 0;
}
