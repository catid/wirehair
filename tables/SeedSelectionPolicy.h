#pragma once

#include <vector>

namespace wirehair_table {

inline unsigned CurrentFirstCandidate(
    unsigned rank,
    unsigned current,
    unsigned candidate_count)
{
    if (rank >= candidate_count || current >= candidate_count) {
        return candidate_count;
    }
    if (rank == 0) {
        return current;
    }
    return rank <= current ? rank - 1u : rank;
}

inline std::vector<unsigned> CurrentFirstRange(
    unsigned current,
    unsigned first,
    unsigned last)
{
    std::vector<unsigned> values;
    if (first > last) {
        return values;
    }
    if (current >= first && current <= last) {
        values.push_back(current);
    }
    for (unsigned value = first;; ++value) {
        if (value != current) {
            values.push_back(value);
        }
        if (value == last) {
            break;
        }
    }
    return values;
}

inline bool IsBetterFailureScore(unsigned candidate, unsigned best)
{
    return candidate < best;
}

} // namespace wirehair_table
