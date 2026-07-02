#include "WirehairV2Precode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>

namespace wirehair_v2 {
namespace {

// Unbiased uniform in [0, bound) via Lemire multiply-shift rejection,
// matching the simulator's Rng::Below that the certification ran with
uint32_t UniformBelow(wirehair::PCGRandom& prng, uint32_t bound)
{
    if (bound <= 1u) {
        return 0;
    }
    // 2^32 mod bound, computed in 32-bit arithmetic
    const uint32_t threshold = (0u - bound) % bound;
    for (;;)
    {
        const uint32_t x = prng.Next();
        const uint64_t m = (uint64_t)x * bound;
        if ((uint32_t)m >= threshold) {
            return (uint32_t)(m >> 32);
        }
    }
}

// Sattolo-style inside-out deck shuffle with UNBIASED position draws.
// Production ShuffleDeck16 draws 16-bit chunks modulo ii, which is fine at
// its <= ~500-entry dense-count sizes but measurably biased once the deck
// spans the whole precode (set-half membership 0.455-0.577 instead of 0.500
// at span ~64000).  The Phase B certification assumed unbiased decks, so
// this keeps the same deck structure with a rejection-sampled draw.
void UnbiasedShuffleDeck(
    wirehair::PCGRandom& prng,
    uint16_t* deck,
    uint32_t count)
{
    deck[0] = 0;
    for (uint32_t ii = 1; ii < count; ++ii)
    {
        const uint32_t jj = UniformBelow(prng, ii);
        deck[ii] = deck[jj];
        deck[jj] = (uint16_t)ii;
    }
}

} // namespace

PrecodeParams MakeCertifiedParams(uint32_t block_count, uint64_t seed)
{
    PrecodeParams params;
    params.BlockCount = block_count;
    params.Staircase =
        (block_count >= 2u && block_count <= 64000u) ?
        wirehair::GetDenseCount(block_count) : 0u;
    params.DenseRows = 12u;
    params.HeavyRows = 12u;
    params.SourceHits = 2u;
    params.Seed = seed;
    return params;
}

bool BuildPrecodeSystem(const PrecodeParams& params, PrecodeSystem& out)
{
    const uint32_t K = params.BlockCount;
    const uint32_t S = params.Staircase;
    const uint32_t D2 = params.DenseRows;
    const uint32_t N1 = params.SourceHits;
    const uint32_t span = K + S + D2;

    // N1 is capped at 8 like the simulator's n1 clamp (and the picks
    // scratch array below)
    if (K < 2u || K > 64000u || S == 0u || N1 == 0u || N1 > 8u) {
        return false;
    }
    // ShuffleDeck16 decks are uint16_t
    if (span > 65535u) {
        return false;
    }

    out.Params = params;
    out.StaircaseRows.assign(S, std::vector<uint32_t>());
    out.DenseRowColumns.assign(D2, std::vector<uint32_t>());

    wirehair::PCGRandom prng;
    prng.Seed(params.Seed, K);

    // --- Staircase rows ---
    // Reserve for the expected load: K*N1/S source hits + parity + link
    {
        const uint32_t expected = (K * N1) / S + 3u;
        for (uint32_t j = 0; j < S; ++j) {
            out.StaircaseRows[j].reserve(expected);
        }
    }

    // Each source column hits min(N1, S) distinct parities.  Distinctness
    // by retry: the tiny-K dense-count table goes as low as S=2 (K=2), but
    // hits <= S always leaves a free residue, so the loop terminates.
    {
        uint32_t picks[8];
        const uint32_t hits = N1 < S ? N1 : S;
        for (uint32_t c = 0; c < K; ++c)
        {
            for (uint32_t hit = 0; hit < hits; ++hit)
            {
                uint32_t p;
                bool collide;
                do {
                    p = UniformBelow(prng, S);
                    collide = false;
                    for (uint32_t j = 0; j < hit; ++j) {
                        if (picks[j] == p) {
                            collide = true;
                            break;
                        }
                    }
                } while (collide);
                picks[hit] = p;
                out.StaircaseRows[p].push_back(c);
            }
        }
    }

    // Own parity column and staircase link.  Source entries were appended
    // in ascending c, so appending K+j-1 then K+j keeps rows sorted, and
    // no column can appear twice (per-source hits are distinct parities).
    for (uint32_t j = 0; j < S; ++j)
    {
        if (j > 0u) {
            out.StaircaseRows[j].push_back(K + j - 1u);
        }
        out.StaircaseRows[j].push_back(K + j);
    }

    // --- Shuffle-2 dense rows ---
    // Documented rule (experiments/precode/README.md, D2 section):
    //   1. deck = shuffle over span columns; set_count = ceil(span/2);
    //      first row = deck[0 .. set_count).
    //   2. Reshuffle; next floor(D2/2) rows each flip {deck[ii],
    //      deck[set_count + ii]}, ii = 0, 1, ...
    //   3. Reshuffle again; remaining floor(D2/2) - 1 + (D2 & 1) rows by
    //      the same flip rule, ii restarting at 0.
    if (D2 > 0u)
    {
        const uint32_t set_count = (span + 1u) >> 1;
        std::vector<uint16_t> deck(span);
        std::vector<uint8_t> bitmap(span, 0);

        UnbiasedShuffleDeck(prng, deck.data(), span);
        for (uint32_t i = 0; i < set_count; ++i) {
            bitmap[deck[i]] = 1;
        }

        uint32_t row_i = 0;
        const auto emit_row = [&]() {
            std::vector<uint32_t>& columns = out.DenseRowColumns[row_i++];
            columns.reserve(set_count + 8u);
            for (uint32_t col = 0; col < span; ++col) {
                if (bitmap[col]) {
                    columns.push_back(col);
                }
            }
        };
        emit_row();

        const uint32_t halves[2] = {
            D2 >> 1,
            (D2 >> 1) + (D2 & 1u) - 1u
        };
        for (uint32_t half = 0; half < 2u; ++half)
        {
            UnbiasedShuffleDeck(prng, deck.data(), span);
            for (uint32_t ii = 0; ii < halves[half]; ++ii)
            {
                // Deck entries at distinct positions are distinct columns,
                // so each flip pair changes exactly two columns
                bitmap[deck[ii]] ^= 1u;
                bitmap[deck[set_count + ii]] ^= 1u;
                emit_row();
            }
        }
    }

    return true;
}

uint8_t HeavyCoefficient(
    uint32_t heavy_row,
    uint32_t ge_column,
    uint32_t heavy_rows)
{
    // Past 128 rows the X window shrinks below H and the per-window MDS
    // claim silently degrades
    CAT_DEBUG_ASSERT(heavy_rows >= 1u && heavy_rows <= 128u);
    CAT_DEBUG_ASSERT(heavy_row < heavy_rows);

    // Y values occupy [0, H); X values occupy [H, 256).  X ^ Y is never
    // zero because the sets are disjoint, and within one 244-column window
    // all X are distinct, giving the Cauchy MDS property there.
    const uint32_t x = heavy_rows + (ge_column % (256u - heavy_rows));
    const uint32_t y = heavy_row;
    return gf256_inv((uint8_t)(x ^ y));
}

} // namespace wirehair_v2
