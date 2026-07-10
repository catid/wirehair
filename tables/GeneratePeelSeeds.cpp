#include "../test/SiameseTools.h"

#include "../WirehairCodec.h"
#include "../WirehairTools.h"
#include "PeelSelectionPolicy.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <string>
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
    GeneratePeelSeeds.cpp generates the dense seeds.

    It tries all the peel seeds for subsets of data until it finds one that
    works best.
*/


//// Entrypoint

static bool SkipTuning = false;
static uint64_t CampaignSeed = 123456789;

static const int kMaxTuningTries = 8;

static const unsigned kPeelSeedSubdivisions = wirehair::kPeelSeedSubdivisions;

static uint8_t kPeelSeeds[kPeelSeedSubdivisions] = {
    // TODO
};

void FillSeeds()
{
    memcpy(kPeelSeeds, wirehair::kPeelSeeds, sizeof(kPeelSeeds));
}

static std::atomic<unsigned> FailedTrials(0);
static const unsigned kHardFailurePenalty = 1000000;

uint8_t Message[64000];

static uint64_t DeriveTrialSeed(uint64_t base_seed, unsigned N)
{
    uint64_t value = (base_seed ^ 0x5045454cULL) +
        0x9e3779b97f4a7c15ULL * (uint64_t)N;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

static uint64_t BaseTablesHash()
{
    uint64_t hash = 1469598103934665603ULL;
    for (unsigned N = 2048; N <= CAT_WIREHAIR_MAX_N; ++N)
    {
        const uint16_t count = wirehair::GetDenseCount(N);
        const uint16_t denseSeed = wirehair::GetBaseDenseSeed(N, count);
        hash ^= N;
        hash *= 1099511628211ULL;
        hash ^= count;
        hash *= 1099511628211ULL;
        hash ^= denseSeed;
        hash *= 1099511628211ULL;
    }
    for (unsigned subdivision = 0;
         subdivision < kPeelSeedSubdivisions;
         ++subdivision)
    {
        hash ^= kPeelSeeds[subdivision];
        hash *= 1099511628211ULL;
    }
    return hash;
}

static void QuickReject(
    unsigned N,
    const uint16_t p_seed)
{
    wirehair::Codec encoder;

    const uint16_t dense_count = wirehair::GetDenseCount(N);
    const uint16_t dense_seed = wirehair::GetBaseDenseSeed(N, dense_count);

    // Override the seeds
    encoder.OverrideSeeds(dense_count, p_seed, dense_seed);

    // Initialize codec
    WirehairResult result = encoder.InitializeEncoder(N, 1);

    // If initialization succeeded:
    if (result == Wirehair_Success) {
        // Feed message to codec
        result = encoder.EncodeFeed(&Message[0]);
    }

    if (result != Wirehair_Success) {
        // Huge penalty
        ++FailedTrials;
    }
}

static void QuickPeelTest(
    unsigned N,
    const uint16_t p_seed,
    unsigned trial,
    unsigned loss_percent)
{
    wirehair::Codec encoder, decoder;

    const uint16_t dense_count = wirehair::GetDenseCount(N);
    const uint16_t dense_seed = wirehair::GetBaseDenseSeed(N, dense_count);

    // Override the seeds
    encoder.OverrideSeeds(dense_count, p_seed, dense_seed);
    decoder.OverrideSeeds(dense_count, p_seed, dense_seed);

    // Initialize codec
    WirehairResult result = encoder.InitializeEncoder(N, 1);

    // If initialization succeeded:
    if (result == Wirehair_Success) {
        // Feed message to codec
        result = encoder.EncodeFeed(&Message[0]);
    }

    if (result != Wirehair_Success) {
        // Huge penalty
        cout << "InitializeEncoder/EncodeFeed failed" << endl;
        FailedTrials += kHardFailurePenalty;
        return;
    }

    result = decoder.InitializeDecoder(N, 1);

    if (result != Wirehair_Success) {
        cout << "InitializeDecoder failed" << endl;
        // Huge penalty
        FailedTrials += kHardFailurePenalty;
        return;
    }

    unsigned added = 0;

    // Loss patterns must depend only on (N, trial): folding the candidate
    // p_seed in scores each candidate on different loss sets, biasing the
    // minimum-failure pick toward lucky streams (winner's curse)
    wirehair::PCGRandom prng;
    prng.Seed(DeriveTrialSeed(CampaignSeed, N), trial);

    const unsigned lossCount =
        (N * loss_percent + 99u) / 100u;
    std::vector<uint16_t> losses(N);

    wirehair::ShuffleDeck16(prng, &losses[0], N);
    const std::vector<uint8_t> dropped =
        wirehair_table::BuildPeelDropMask(losses, lossCount);

    const uint64_t packetIdCount = wirehair_table::PeelPacketIdCount(N);
    for (uint64_t packetId = 0; packetId < packetIdCount; ++packetId)
    {
        if (packetId < dropped.size() && dropped[(size_t)packetId]) {
            continue;
        }

        ++added;

        uint8_t encodedData[1];

        const uint32_t encodedBytes = encoder.Encode(
            (uint32_t)packetId, encodedData, 1);

        if (encodedBytes == 0)
        {
            cout << "Encode failed" << endl;
            // Huge penalty
            FailedTrials += kHardFailurePenalty;
            return;
        }

        result = decoder.DecodeFeed(
            (uint32_t)packetId, encodedData, encodedBytes);

        if (result != Wirehair_NeedMore)
        {
            if (result == Wirehair_Success)
            {
                if (added == N) {
                    // No penalty
                }
                else if (added == N + 1) {
                    FailedTrials++;
                }
                else if (added == N + 2) {
                    FailedTrials += 3;
                }
                else if (added >= N + 3) {
                    FailedTrials += (added - N) * 3;
                }
                return;
            }

            cout << "DecodeFeed failed" << endl;
            // Huge penalty
            FailedTrials += kHardFailurePenalty;
            return;
        }
    }
    cout << "packet-id cap reached for N = " << N << endl;
    FailedTrials += kHardFailurePenalty;
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
        << "usage: " << program << " [--seed U64] [--trials N] [--sublo I] "
        << "[--subhi I] [--nlo N] [--nhi N] [--max-tries N] "
        << "[--loss10-percent N] [--loss30-percent N] "
        << "[--skip-tuning] [--results FILE]\n";
}

static unsigned FirstNForSubdivision(unsigned subdivision, unsigned nlo)
{
    const unsigned first = 2048u + subdivision;
    if (first >= nlo) {
        return first;
    }
    const unsigned delta = nlo - first;
    const unsigned steps =
        (delta + kPeelSeedSubdivisions - 1u) / kPeelSeedSubdivisions;
    return first + steps * kPeelSeedSubdivisions;
}

int main(int argc, char** argv)
{
    unsigned trials = 10;
    unsigned sublo = 0;
    unsigned subhi = kPeelSeedSubdivisions - 1u;
    unsigned nlo = 2048;
    unsigned nhi = CAT_WIREHAIR_MAX_N;
    unsigned maxTries = kMaxTuningTries;
    unsigned loss10Percent = 10;
    unsigned loss30Percent = 30;
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
            if (!ParseU64(next(), CampaignSeed)) {
                cerr << "bad --seed value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--trials")) {
            if (!ParseUnsigned(next(), trials) || trials == 0 ||
                trials > wirehair_table::MaxPeelTrials(kHardFailurePenalty))
            {
                cerr << "bad --trials value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--sublo")) {
            if (!ParseUnsigned(next(), sublo)) {
                cerr << "bad --sublo value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--subhi")) {
            if (!ParseUnsigned(next(), subhi)) {
                cerr << "bad --subhi value" << endl;
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
        else if (!strcmp(arg, "--max-tries")) {
            if (!ParseUnsigned(next(), maxTries) || maxTries == 0) {
                cerr << "bad --max-tries value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--loss10-percent")) {
            if (!ParseUnsigned(next(), loss10Percent)) {
                cerr << "bad --loss10-percent value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--loss30-percent")) {
            if (!ParseUnsigned(next(), loss30Percent)) {
                cerr << "bad --loss30-percent value" << endl;
                return 1;
            }
        }
        else if (!strcmp(arg, "--skip-tuning")) {
            SkipTuning = true;
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

    if (sublo >= kPeelSeedSubdivisions ||
        subhi >= kPeelSeedSubdivisions ||
        sublo > subhi)
    {
        cerr << "subdivision range must be in [0, "
            << (kPeelSeedSubdivisions - 1u) << "]" << endl;
        return 1;
    }
    if (nlo < 2048u || nhi > CAT_WIREHAIR_MAX_N || nlo > nhi) {
        cerr << "N range must be in [2048, " << CAT_WIREHAIR_MAX_N
            << "]" << endl;
        return 1;
    }
    if (loss10Percent < 1u || loss30Percent > 99u ||
        loss10Percent >= loss30Percent)
    {
        cerr << "loss percentages must satisfy 1 <= loss10 < loss30 <= 99"
            << endl;
        return 1;
    }

    const int gfInitResult = gf256_init();

    FillSeeds();

    // If gf256 init failed:
    if (gfInitResult != 0)
    {
        cout << "GF256 init failed" << endl;
        return -1;
    }

    for (unsigned i = 0; i < 64000; ++i) {
        Message[i] = (uint8_t)i;
    }

    ofstream results;
    const uint64_t baseTablesHash = BaseTablesHash();
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
        results << "Subdivision\tFirstN\tLastN\tPeelSeed\tWorstScore"
            "\tTotalScore\tTotalScore10\tTotalScore30"
            "\tCurrentPeelSeed\tTrials\tMaxTries"
            "\tEvaluatedCandidates\tEvaluatedScenarios\tPrunedCandidates"
            "\tScenarioCount\tBaseSeed\tRequestedNLo\tRequestedNHi"
            "\tBaseTablesHash\tLoss10Percent\tLoss30Percent"
            "\tMethodVersion\tStatus\n";
    }

    // Subdivisions 0 and 1 map real block counts (N = 2048+i stepping 2048,
    // i.e. 2048, 2049, 4096, 4097, ...) and need tuning like every other
    // residue class; starting at 2 carried forward untested table entries
    cerr
        << "# GeneratePeelSeeds seed=0x" << hex << CampaignSeed << dec
        << " trials=" << trials
        << " sublo=" << sublo
        << " subhi=" << subhi
        << " nlo=" << nlo
        << " nhi=" << nhi
        << " max_tries=" << maxTries
        << " loss10_percent=" << loss10Percent
        << " loss30_percent=" << loss30Percent
        << " skip_tuning=" << (SkipTuning ? 1 : 0) << endl;

    for (unsigned i = sublo; i <= subhi; ++i)
    {
        const unsigned current_peel_seed = kPeelSeeds[i];
        unsigned best_peel_seed = current_peel_seed;
        wirehair_table::PeelCandidateScore bestScore =
            wirehair_table::PeelCandidateScore::WorstPossible();
        bool found_candidate = false;
        bool selectedWithoutTuning = false;

        unsigned evaluatedCandidates = 0;
        const unsigned firstN = FirstNForSubdivision(i, nlo);
        if (firstN > nhi)
        {
            cout << "subdivision = " << i
                << " : Skipped no N in requested range; kept seed = "
                << (int)kPeelSeeds[i] << endl;
            if (!resultsPath.empty())
            {
                results << i << "\t0\t0\t" << current_peel_seed
                    << "\t0\t0\t0\t0\t" << current_peel_seed << '\t' << trials
                    << '\t' << maxTries << "\t0\t0\t0\t0\t" << CampaignSeed
                    << '\t' << nlo << '\t' << nhi << '\t' << baseTablesHash
                    << '\t' << loss10Percent << '\t' << loss30Percent
                    << "\tpeel-v4\tretained-no-N\n";
            }
            continue;
        }
        const unsigned lastN = firstN +
            ((nhi - firstN) / kPeelSeedSubdivisions) *
                kPeelSeedSubdivisions;
        vector<int> scenarioNs;
        for (int N = (int)firstN;
             N <= (int)nhi;
             N += (int)kPeelSeedSubdivisions)
        {
            scenarioNs.push_back(N);
        }
        reverse(scenarioNs.begin(), scenarioNs.end());
        const unsigned scenarioCount = (unsigned)scenarioNs.size() * 2u;
        uint64_t evaluatedScenarios = 0;
        unsigned prunedCandidates = 0;

        // Evaluate the shipped seed first so exact score ties retain it.
        for (unsigned rank = 0; rank < 256; ++rank)
        {
            const unsigned p_seed = rank == 0 ? current_peel_seed :
                (rank <= current_peel_seed ? rank - 1u : rank);
            FailedTrials = 0;

#pragma omp parallel for
            for (int N = (int)firstN;
                 N <= (int)nhi;
                 N += (int)kPeelSeedSubdivisions)
            {
                QuickReject(N, (uint16_t)p_seed);
            }

            if (FailedTrials > 0) {
                continue;
            }

            cout << "For subdivision " << i << " of " << kPeelSeedSubdivisions << " - testing seed " << p_seed << endl;
            ++evaluatedCandidates;

            if (SkipTuning)
            {
                best_peel_seed = p_seed;
                bestScore = wirehair_table::PeelCandidateScore();
                found_candidate = true;
                selectedWithoutTuning = true;
                break;
            }

            wirehair_table::PeelCandidateScore candidateScore;
            bool pruned = false;
            for (int N : scenarioNs)
            {
                // Evaluate the harder loss first to reject weak candidates
                // earlier without changing the complete score.
                FailedTrials = 0;
#pragma omp parallel for
                for (int trial = 0; trial < (int)trials; ++trial) {
                    QuickPeelTest(
                        N, (uint16_t)p_seed, (unsigned)trial, loss30Percent);
                }
                candidateScore.Add(FailedTrials, true);
                ++evaluatedScenarios;
                if (wirehair_table::CannotBeatPeelScore(
                        candidateScore, bestScore))
                {
                    pruned = true;
                    break;
                }

                FailedTrials = 0;
#pragma omp parallel for
                for (int trial = 0; trial < (int)trials; ++trial) {
                    QuickPeelTest(
                        N, (uint16_t)p_seed, (unsigned)trial, loss10Percent);
                }
                candidateScore.Add(FailedTrials, false);
                ++evaluatedScenarios;
                if (wirehair_table::CannotBeatPeelScore(
                        candidateScore, bestScore))
                {
                    pruned = true;
                    break;
                }
            }
            if (pruned) {
                ++prunedCandidates;
                cout << "pruned worst = " << candidateScore.WorstScore
                    << " partial total = " << candidateScore.TotalScore << endl;
                if (evaluatedCandidates >= maxTries) {
                    break;
                }
                continue;
            }

            cout << "total10 = " << candidateScore.TotalScore10
                << " total30 = " << candidateScore.TotalScore30
                << " worst = " << candidateScore.WorstScore
                << " total = " << candidateScore.TotalScore << endl;

            if (wirehair_table::IsBetterPeelScore(candidateScore, bestScore))
            {
                best_peel_seed = p_seed;
                bestScore = candidateScore;
                found_candidate = true;

                cout << "*** subdivision = " << i << " : Picked seed = "
                    << best_peel_seed << " worst = " << bestScore.WorstScore
                    << " total = " << bestScore.TotalScore << endl;

                if (bestScore.WorstScore == 0) {
                    break;
                }
            }

            if (evaluatedCandidates >= maxTries) {
                break;
            }
        }

        if (!found_candidate)
        {
            cout << "subdivision = " << i
                << " : No admissible candidate; retained seed = "
                << current_peel_seed << endl;
        }

        cout << "subdivision = " << i << " : Picked seed = " << best_peel_seed
            << " worst = " << bestScore.WorstScore
            << " total = " << bestScore.TotalScore
            << " total10 = " << bestScore.TotalScore10
            << " total30 = " << bestScore.TotalScore30
            << endl;

        kPeelSeeds[i] = (uint8_t)best_peel_seed;
        if (!resultsPath.empty())
        {
            results << i << '\t' << firstN << '\t' << lastN << '\t'
                << best_peel_seed << '\t' << bestScore.WorstScore << '\t'
                << bestScore.TotalScore << '\t' << bestScore.TotalScore10 << '\t'
                << bestScore.TotalScore30 << '\t' << current_peel_seed << '\t'
                << trials << '\t'
                << maxTries << '\t' << evaluatedCandidates << '\t'
                << evaluatedScenarios << '\t' << prunedCandidates << '\t'
                << scenarioCount << '\t'
                << CampaignSeed << '\t' << nlo << '\t' << nhi << '\t'
                << baseTablesHash << '\t' << loss10Percent << '\t'
                << loss30Percent << "\tpeel-v4\t"
                << (!found_candidate ? "retained-no-candidate" :
                    selectedWithoutTuning ? "selected-skip-tuning" : "searched")
                << '\n';
            if (!results) {
                cerr << "failed writing --results file: " << resultsPath << endl;
                return 1;
            }
        }
    }

    cout << "static const unsigned kPeelSeedSubdivisions = " << kPeelSeedSubdivisions << ";" << endl;

    cout << "const uint8_t kPeelSeeds[kPeelSeedSubdivisions] = {" << endl;

    const unsigned modulus = 32;

    for (unsigned i = 0; i < kPeelSeedSubdivisions; ++i)
    {
        if (i % modulus == 0) {
            cout << "    ";
        }

        cout << (int)kPeelSeeds[i] << ",";

        if ((i + 1) % modulus == 0) {
            cout << endl;
        }
    }

    cout << "};" << endl;


    return FinalizeResults(
        results, resultsTemporaryPath, resultsPath) ? 0 : 1;
}
