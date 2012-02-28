/*
	Copyright (c) 2012 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of WirehairFEC nor the names of its contributors may be
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

#include "GF2Matrix.hpp"
#include "../codec_source/SmallPRNG.hpp"
using namespace cat;
using namespace wirehair;

#include <iostream>
using namespace std;

GF2Matrix::GF2Matrix()
{
	_matrix = 0;
	_pivot = 0;
}

GF2Matrix::~GF2Matrix()
{
	Cleanup();
}

void GF2Matrix::Cleanup()
{
	if (_matrix)
	{
		delete []_matrix;
		_matrix = 0;
	}
	if (_pivot)
	{
		delete []_pivot;
		_pivot = 0;
	}
}

void GF2Matrix::Zero()
{
	memset(_matrix, 0, _pitch * _n * sizeof(u64));
}

void GF2Matrix::Identity()
{
	Zero();

	// For each row,
	u64 *row = _matrix;
	for (int column_i = 0; column_i < _n; ++column_i, row += _pitch)
	{
		// Flip diagonal bit
		u64 mask = (u64)1 << (column_i & 63);
		row[column_i >> 6] |= mask;
	}
}

void GF2Matrix::Fill()
{
	CatsChoice prng;
	prng.Initialize(_seed);

	u64 *row = _matrix;
	int words = _pitch * _n;
	while (words--)
	{
		u32 rv1 = prng.Next();
		u32 rv2 = prng.Next();
		u64 word = ((u64)rv2 << 32) | rv1;
		*row++ = word;
	}
}

bool GF2Matrix::Triangle()
{
	CAT_IF_DEBUG(cout << endl << "---- Triangle ----" << endl << endl;)

	// Initialize pivot array
	for (u16 pivot_i = 0; pivot_i < _n; ++pivot_i)
		_pivot[pivot_i] = pivot_i;

	// For each pivot to determine,
	u64 ge_mask = 1;
	for (u16 pivot_i = 0; pivot_i < _n; ++pivot_i)
	{
		int word_offset = pivot_i >> 6;
		u64 *ge_matrix_offset = _matrix + word_offset;

		// Find pivot
		bool found = false;
		for (u16 pivot_j = pivot_i; pivot_j < _n; ++pivot_j)
		{
			// Determine if the row contains the bit we want
			u16 ge_row_j = _pivot[pivot_j];
			u64 *ge_row = &ge_matrix_offset[_pitch * ge_row_j];

			// If the bit was found,
			if (*ge_row & ge_mask)
			{
				found = true;
				CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " found on row " << ge_row_j << endl;)

				// Swap out the pivot index for this one
				u16 temp = _pivot[pivot_i];
				_pivot[pivot_i] = _pivot[pivot_j];
				_pivot[pivot_j] = temp;

				// Unroll first word
				u64 row0 = *ge_row;

				// For each remaining unused row,
				for (u16 pivot_k = pivot_j + 1; pivot_k < _n; ++pivot_k)
				{
					// Determine if the row contains the bit we want
					u16 ge_row_k = _pivot[pivot_k];
					u64 *rem_row = &ge_matrix_offset[_pitch * ge_row_k];

					// If the bit was found,
					if (*rem_row & ge_mask)
					{
						// Add the pivot row to eliminate the bit from this row, preserving previous bits
						*rem_row ^= row0;

						for (int ii = 1; ii < _pitch - word_offset; ++ii)
							rem_row[ii] ^= ge_row[ii];
					}
				}

				break;
			}
		}

		// If pivot could not be found,
		if (!found)
		{
			CAT_IF_DEBUG(cout << "Inversion impossible: Pivot " << pivot_i << " not found!" << endl;)
			return false;
		}

		// Generate next mask
		ge_mask = CAT_ROL64(ge_mask, 1);
	}

	return true;
}

void GF2Matrix::Diagonal()
{
	CAT_IF_DEBUG(cout << endl << "---- Diagonal ----" << endl << endl;)

	// For each pivot row from last to first,
	int pivot_i = _n - 1;
	u64 ge_mask = (u64)1 << (pivot_i & 63);
	for (; pivot_i >= 0; --pivot_i)
	{
		u16 pivot_ge_row_i = _pivot[pivot_i];

		CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " solving row " << pivot_ge_row_i << endl;)

		// For each pivot row above it,
		int word_offset = pivot_i >> 6;
		u64 *ge_row = _matrix + word_offset;
		for (int above_i = pivot_i - 1; above_i >= 0; --above_i)
		{
			u16 ge_above_row_i = _pivot[above_i];

			// If bit is set in that row,
			u64 *above_row = &ge_row[_pitch * ge_above_row_i];
			if (*above_row & ge_mask)
				*above_row ^= ge_mask;
		}

		// Generate next mask
		ge_mask = CAT_ROR64(ge_mask, 1);
	}
}

bool GF2Matrix::Initialize(int n)
{
	Cleanup();

	CAT_IF_DEBUG(cout << endl << "GF2Matrix.Initialize: Seed = " << _seed << " n = " << n << endl << endl;)

	_n = n;
	_pitch = (n + 63) / 64;
	int words = _pitch * n;
	_matrix = new u64[words];
	if (!_matrix) return false;
	_pivot = new u16[n];
	if (!_pivot) return false;
	_seed = 0;

	return true;
}

void GF2Matrix::Print()
{
	int rows = _n;
	int cols = _n;

	cout << endl << "GF2Matrix is " << rows << " x " << cols << " (seed " << _seed << "):" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			if (_matrix[_pitch * ii + (jj >> 6)] & ((u64)1 << (jj & 63)))
				cout << '1';
			else
				cout << '0';
		}
		cout << endl;
	}

	cout << endl;
}
