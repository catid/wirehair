#pragma once

#include <algorithm>
#include <climits>
#include <cstdint>
#include <limits>
#include <vector>

namespace wirehair_table {

struct PeelCandidateScore
{
    unsigned WorstScore;
    uint64_t TotalScore;
    uint64_t TotalScore10;
    uint64_t TotalScore30;

    PeelCandidateScore()
        : WorstScore(0)
        , TotalScore(0)
        , TotalScore10(0)
        , TotalScore30(0)
    {
    }

    static PeelCandidateScore WorstPossible()
    {
        PeelCandidateScore score;
        score.WorstScore = UINT_MAX;
        score.TotalScore = std::numeric_limits<uint64_t>::max();
        score.TotalScore10 = std::numeric_limits<uint64_t>::max();
        score.TotalScore30 = std::numeric_limits<uint64_t>::max();
        return score;
    }

    void Add(unsigned score, bool secondary_loss)
    {
        WorstScore = std::max(WorstScore, score);
        TotalScore += score;
        if (secondary_loss) {
            TotalScore30 += score;
        }
        else {
            TotalScore10 += score;
        }
    }
};

inline bool IsBetterPeelScore(
    const PeelCandidateScore& candidate,
    const PeelCandidateScore& best)
{
    return candidate.WorstScore < best.WorstScore ||
        (candidate.WorstScore == best.WorstScore &&
         candidate.TotalScore < best.TotalScore);
}

inline bool CannotBeatPeelScore(
    const PeelCandidateScore& partial,
    const PeelCandidateScore& best)
{
    return partial.WorstScore > best.WorstScore ||
        (partial.WorstScore == best.WorstScore &&
         partial.TotalScore >= best.TotalScore);
}

inline std::vector<uint8_t> BuildPeelDropMask(
    const std::vector<uint16_t>& shuffled_ids,
    unsigned loss_count)
{
    std::vector<uint8_t> dropped(shuffled_ids.size(), 0);
    const unsigned count = std::min<unsigned>(
        loss_count, (unsigned)shuffled_ids.size());
    for (unsigned i = 0; i < count; ++i) {
        dropped[shuffled_ids[i]] = 1;
    }
    return dropped;
}

inline uint64_t PeelPacketIdCount(unsigned N)
{
    return (uint64_t)N * 2u + 512u;
}

inline unsigned MaxPeelTrials(unsigned hard_failure_penalty)
{
    return hard_failure_penalty == 0 ? 0 : UINT_MAX / hard_failure_penalty;
}

} // namespace wirehair_table
