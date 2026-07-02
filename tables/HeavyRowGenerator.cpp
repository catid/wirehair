/** \file
    \brief Wirehair : Heavy Row Generator
    \copyright Copyright (c) 2018 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Wirehair nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdint.h>

#include <iostream>
#include <iomanip>
using namespace std;

static const unsigned kRows = 6;
static const unsigned kColumns = 18;

#include "../gf256.h"
#include "../test/SiameseTools.h"
using namespace siamese;

static void EliminateColumn(
    PCGRandom& prng,
    unsigned columnStart,
    uint8_t* matrix)
{
    uint8_t rowMultiplier[kRows];

    // Get value for this column of each row
    for (unsigned row = 0; row < kRows; ++row)
    {
        const uint8_t value = matrix[row * kColumns + columnStart];

        rowMultiplier[row] = value;
    }

    unsigned column = columnStart;

    for (;;)
    {
        uint32_t word = prng.Next();

        for (unsigned i = 0; i < 32; ++i)
        {
            if (0 != (word & 1))
            {
                for (unsigned row = 0; row < kRows; ++row)
                {
                    const uint8_t value = rowMultiplier[row];

                    matrix[row * kColumns + column] ^= value;
                }
            }

            word >>= 1;
            ++column;

            if (column >= kColumns) {
                return;
            }
        }
    }
}

static bool IsFullRank(
    uint8_t * matrix,
    const unsigned offset
)
{
    uint8_t pivot[kRows];

    // Initialize pivot array
    for (unsigned i = 0; i < kRows; ++i) {
        pivot[i] = (uint8_t)i;
    }

    // For each pivot to determine:
    for (unsigned i = 0; i < kRows; ++i)
    {
        bool found = false;

        // Find pivot
        for (unsigned j = i; j < kRows; ++j)
        {
            // Determine if the row contains a nonzero value in this column
            const unsigned ge_row_j = pivot[j];
            const uint8_t * ge_row = matrix + kColumns * ge_row_j;
            const uint8_t row0 = ge_row[offset + i];

            // If the value was nonzero:
            if (0 != row0)
            {
                found = true;

                // Swap out the pivot index for this one
                if (i != j)
                {
                    const uint8_t temp = pivot[i];
                    pivot[i] = pivot[j];
                    pivot[j] = temp;
                }

                // For each remaining unused row:
                for (unsigned k = j + 1; k < kRows; ++k)
                {
                    // Determine if the row contains the bit we want
                    const unsigned ge_row_k = pivot[k];
                    uint8_t * rem_row = matrix + kColumns * ge_row_k;
                    const uint8_t rem_value = rem_row[offset + i];

                    // If the value was nonzero:
                    if (0 != rem_value)
                    {
                        const uint8_t elim = gf256_div(rem_value, row0);

                        // Eliminate within the offset submatrix being
                        // rank-checked: pivots are read at offset + i, so
                        // the row update must start there too.
                        gf256_muladd_mem(
                            rem_row + offset + i,
                            elim,
                            ge_row + offset + i,
                            kRows - i);
                    }
                }

                break;
            }
        }

        // If pivot could not be found,
        if (!found) {
            return false;
        }
    }

    return true;
}

// Historical note: this generator originally searched for a "perfect" seed
// whose right 6x6 submatrix stayed invertible under every random binary
// perturbation of the left columns.  The rank check used to eliminate the
// wrong columns, so it accepted the shipped seed with zero failures.  With
// the corrected elimination, the measured singular rate for the shipped
// matrix is ~0.39% per trial -- the same as a uniformly random 6x6 GF(256)
// matrix (~1/256) -- and no seed can reach zero over 1M trials, so the
// perfect-seed search was removed.  The generator now verifies that it
// reproduces the shipped kHeavyMatrix bytes (the wire-format regression
// gate) and reports the perturbation failure rate as information only.

// Params for a Cauchy matrix that is offset from 0
struct OffsetCauchyMatrixParams
{
    unsigned Rows, Cols, RowWiggle, ColWiggle;

    bool Initialize(const uint64_t seed)
    {
        PCGRandom prng;
        prng.Seed(seed);

        Rows = (prng.Next() % 256) + 1;
        if (Rows < kRows) {
            return false;
        }
        Cols = 256 - Rows;
        if (Cols < kColumns) {
            return false;
        }

        RowWiggle = Rows - kRows;
        ColWiggle = Cols - kColumns;
        // We could add a random offset here too but not needed to find a solution.
        return true;
    }

    void FillMatrix(uint8_t* matrix)
    {
        for (unsigned i = 0; i < kColumns; ++i)
        {
            for (unsigned j = 0; j < kRows; ++j)
            {
                // Cauchy matrix elements:
                // C_ij = 1 / (X_i - Y_j)

                SIAMESE_DEBUG_ASSERT(ColWiggle + i + Rows < 256);
                const uint8_t X_i = (uint8_t)(ColWiggle + i + Rows);

                SIAMESE_DEBUG_ASSERT(RowWiggle + j < Rows);
                const uint8_t Y_j = (uint8_t)(RowWiggle + j);

                const uint8_t C_ij = gf256_inv(X_i ^ Y_j);

                matrix[j * kColumns + i] = C_ij;
            }
        }
    }
};

static void PrintMatrix(const uint8_t* matrix)
{
    cout << endl;
    cout << "static const unsigned kHeavyRows = " << kRows << ";" << endl;
    cout << "static const unsigned kHeavyCols = " << kColumns << ";" << endl;
    cout << "static const uint8_t kHeavyMatrix[kHeavyRows][kHeavyCols] = {" << endl;

    for (unsigned i = 0; i < kRows; ++i)
    {
        cout << "    { ";

        for (unsigned j = 0; j < kColumns; ++j)
        {
            cout << "0x" << hex << setfill('0') << setw(2) << (unsigned)matrix[i * kColumns + j] << ", ";
        }

        cout << " }," << endl;
    }

    cout << "};" << endl;
    cout << dec << endl;
}

// The kHeavyMatrix bytes shipped in WirehairCodec.cpp.  Changing them breaks
// wire-format compatibility, so the generator must reproduce them exactly.
static const uint8_t kShippedHeavyMatrix[kRows][kColumns] = {
    { 0x85, 0xd3, 0x66, 0xf3, 0x38, 0x95, 0x56, 0xad, 0x57, 0xaf, 0x58, 0x48, 0xbc, 0xfa, 0x02, 0xc5, 0x43, 0xe8, },
    { 0xd3, 0x85, 0xf3, 0x66, 0x95, 0x38, 0xad, 0x56, 0xaf, 0x57, 0x48, 0x58, 0xfa, 0xbc, 0xc5, 0x02, 0xe8, 0x43, },
    { 0x82, 0x22, 0x57, 0xaf, 0x56, 0xad, 0x38, 0x95, 0x66, 0xf3, 0x43, 0xe8, 0x02, 0xc5, 0xbc, 0xfa, 0x58, 0x48, },
    { 0x22, 0x82, 0xaf, 0x57, 0xad, 0x56, 0x95, 0x38, 0xf3, 0x66, 0xe8, 0x43, 0xc5, 0x02, 0xfa, 0xbc, 0x48, 0x58, },
    { 0x51, 0x34, 0x56, 0xad, 0x57, 0xaf, 0x66, 0xf3, 0x38, 0x95, 0x02, 0xc5, 0x43, 0xe8, 0x58, 0x48, 0xbc, 0xfa, },
    { 0x34, 0x51, 0xad, 0x56, 0xaf, 0x57, 0xf3, 0x66, 0x95, 0x38, 0xc5, 0x02, 0xe8, 0x43, 0x48, 0x58, 0xfa, 0xbc, },
};

bool Generate_HeavyRows()
{
    // Stored row-first
    uint8_t matrix[kRows * kColumns];

    gf256_init();

    // The seed that produced the shipped kHeavyMatrix
    const uint64_t seed = 2318331135281;

    OffsetCauchyMatrixParams params;
    if (!params.Initialize(seed))
    {
        cout << "FAILURE: Shipped seed no longer initializes" << endl;
        return false;
    }

    // Regression gate: the generator must reproduce the shipped bytes
    params.FillMatrix(matrix);
    for (unsigned i = 0; i < kRows; ++i)
    {
        for (unsigned j = 0; j < kColumns; ++j)
        {
            if (matrix[i * kColumns + j] != kShippedHeavyMatrix[i][j])
            {
                cout << "FAILURE: Generated matrix does not match shipped "
                    "kHeavyMatrix at row " << i << " col " << j << endl;
                return false;
            }
        }
    }

    // Informational: measure how often the right submatrix goes singular
    // under random binary perturbation of the left columns.  Expect ~1/256.
    unsigned failures = 0;
    static const unsigned kTrials = 1000000;

    for (unsigned trial = 0; trial < kTrials; ++trial)
    {
        params.FillMatrix(matrix);

        PCGRandom mix_prng;
        mix_prng.Seed(seed, trial);

        for (unsigned col = 0; col < (kColumns - kRows); ++col) {
            EliminateColumn(mix_prng, col, matrix);
        }

        if (!IsFullRank(matrix, kColumns - kRows)) {
            ++failures;
        }
    }

    cout << "Shipped Cauchy matrix reproduced from seed = " << seed << endl;
    cout << "* Rows = " << params.Rows << endl;
    cout << "* Cols = " << params.Cols << endl;
    cout << "* RowWiggle = " << params.RowWiggle << endl;
    cout << "* ColWiggle = " << params.ColWiggle << endl;
    cout << "* Perturbation singular rate: " << failures << " / " << kTrials
        << " (expected ~1/256 = 3906)" << endl;

    params.FillMatrix(matrix);

    PrintMatrix(matrix);

    return true;
}
