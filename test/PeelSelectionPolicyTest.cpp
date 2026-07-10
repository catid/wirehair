#include "../tables/PeelSelectionPolicy.h"
#include "../tables/SeedSelectionPolicy.h"

#include <algorithm>
#include <cstdio>
#include <numeric>
#include <vector>

static int Fail(const char* message)
{
    std::fprintf(stderr, "PeelSelectionPolicyTest: %s\n", message);
    return 1;
}

int main()
{
    for (unsigned current = 0; current < 256; ++current)
    {
        std::vector<unsigned> seen(256, 0);
        unsigned selected = 999;
        unsigned best = 1;
        for (unsigned rank = 0; rank < 256; ++rank)
        {
            const unsigned candidate = wirehair_table::CurrentFirstCandidate(
                rank, current, 256);
            if (candidate >= 256 || seen[candidate]++) {
                return Fail("current-first candidate order is not a permutation");
            }
            const unsigned score = 0;
            if (wirehair_table::IsBetterFailureScore(score, best)) {
                selected = candidate;
                best = score;
            }
        }
        if (selected != current) {
            return Fail("finite-trial score tie did not retain current seed");
        }
    }

    const std::vector<unsigned> counts =
        wirehair_table::CurrentFirstRange(7, 5, 9);
    if (counts.size() != 5 || counts[0] != 7) {
        return Fail("current dense count is not evaluated first");
    }

    wirehair_table::PeelCandidateScore concentrated;
    concentrated.Add(10, false);
    concentrated.Add(0, false);
    wirehair_table::PeelCandidateScore balanced;
    balanced.Add(6, false);
    balanced.Add(6, false);
    if (!wirehair_table::IsBetterPeelScore(balanced, concentrated)) {
        return Fail("per-N worst score did not beat a lower aggregate pathology");
    }

    wirehair_table::PeelCandidateScore tie = balanced;
    if (wirehair_table::IsBetterPeelScore(tie, balanced)) {
        return Fail("exact score tie should retain the first candidate");
    }

    const unsigned synthetic[][4] = {
        { 4, 2, 3, 1 },
        { 5, 0, 0, 0 },
        { 3, 3, 3, 3 },
        { 4, 2, 3, 1 },
    };
    wirehair_table::PeelCandidateScore exhaustive_best =
        wirehair_table::PeelCandidateScore::WorstPossible();
    unsigned exhaustive_index = 999;
    for (unsigned candidate = 0; candidate < 4; ++candidate) {
        wirehair_table::PeelCandidateScore score;
        for (unsigned scenario = 0; scenario < 4; ++scenario) {
            score.Add(synthetic[candidate][scenario], scenario % 2 != 0);
        }
        if (wirehair_table::IsBetterPeelScore(score, exhaustive_best)) {
            exhaustive_best = score;
            exhaustive_index = candidate;
        }
    }
    wirehair_table::PeelCandidateScore bounded_best =
        wirehair_table::PeelCandidateScore::WorstPossible();
    unsigned bounded_index = 999;
    for (unsigned candidate = 0; candidate < 4; ++candidate) {
        wirehair_table::PeelCandidateScore score;
        bool pruned = false;
        for (unsigned scenario = 0; scenario < 4; ++scenario) {
            score.Add(synthetic[candidate][scenario], scenario % 2 != 0);
            if (wirehair_table::CannotBeatPeelScore(score, bounded_best)) {
                pruned = true;
                break;
            }
        }
        if (!pruned && wirehair_table::IsBetterPeelScore(score, bounded_best)) {
            bounded_best = score;
            bounded_index = candidate;
        }
    }
    if (bounded_index != exhaustive_index ||
        bounded_best.WorstScore != exhaustive_best.WorstScore ||
        bounded_best.TotalScore != exhaustive_best.TotalScore)
    {
        return Fail("branch-and-bound changed the exhaustive winner");
    }

    for (unsigned N = 2; N <= 257; ++N)
    {
        std::vector<uint16_t> shuffled(N);
        std::iota(shuffled.begin(), shuffled.end(), 0);
        std::rotate(
            shuffled.begin(), shuffled.begin() + (N * 7u) % N,
            shuffled.end());
        for (unsigned loss_percent : {10u, 30u, 99u})
        {
            const unsigned loss_count =
                (N * loss_percent + 99u) / 100u;
            const std::vector<uint8_t> mask =
                wirehair_table::BuildPeelDropMask(shuffled, loss_count);
            for (unsigned id = 0; id < N + 10u; ++id)
            {
                const bool linear = id < N && std::find(
                    shuffled.begin(), shuffled.begin() + loss_count,
                    (uint16_t)id) != shuffled.begin() + loss_count;
                const bool indexed = id < mask.size() && mask[id] != 0;
                if (linear != indexed) {
                    return Fail("drop mask differs from linear membership oracle");
                }
            }
        }
    }

    if (wirehair_table::PeelPacketIdCount(64000) != 128512) {
        return Fail("packet-id cap does not widen safely at maximum N");
    }
    if (wirehair_table::MaxPeelTrials(1000000) != 4294) {
        return Fail("peel trial bound does not prevent unsigned score overflow");
    }
    return 0;
}
