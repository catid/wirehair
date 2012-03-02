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

#include "GF256Matrix.hpp"
#include "../codec_source/SmallPRNG.hpp"
using namespace cat;
using namespace wirehair;

#include <iostream>
#include <iomanip>
using namespace std;


/*
	Branchless multiply and divide construction from
	"Fast Software Implementations of Finite Field Operations (Extended Abstract)"
	by Cheng Huang, Lihao Xu

	Small corrections made to paper (Q = 255):
		+ The ALOG_TABLE needs to have 512*2+1 elements to handle 0*0 = 0 case.
		+ Element 255*2 should be set to 1.

	After these corrections it works properly and reduces the execution time
	to 58% of the usual version that uses branches to handle zero input.

	This table was generated using polynomial 0x15F.  Maybe it's voodoo but
	random GF(256) matrices with this polynomial tended to be more invertible.
*/

static const u16 LOG_TABLE[256] = {
	512, 255, 1, 122, 2, 244, 123, 181, 3, 48, 245, 224, 124, 84, 182, 111,
	4, 233, 49, 19, 246, 107, 225, 206, 125, 56, 85, 170, 183, 91, 112, 250,
	5, 117, 234, 10, 50, 156, 20, 213, 247, 203, 108, 178, 226, 37, 207, 210,
	126, 150, 57, 100, 86, 141, 171, 40, 184, 73, 92, 164, 113, 146, 251, 229,
	6, 96, 118, 15, 235, 193, 11, 13, 51, 68, 157, 195, 21, 31, 214, 237,
	248, 168, 204, 17, 109, 222, 179, 120, 227, 162, 38, 98, 208, 176, 211, 8,
	127, 188, 151, 239, 58, 132, 101, 216, 87, 80, 142, 33, 172, 27, 41, 23,
	185, 77, 74, 197, 93, 65, 165, 159, 114, 200, 147, 70, 252, 45, 230, 53,
	7, 175, 97, 161, 119, 221, 16, 167, 236, 30, 194, 67, 12, 192, 14, 95,
	52, 44, 69, 199, 158, 64, 196, 76, 22, 26, 32, 79, 215, 131, 238, 187,
	249, 90, 169, 55, 205, 106, 18, 232, 110, 83, 223, 47, 180, 243, 121, 254,
	228, 145, 163, 72, 39, 140, 99, 149, 209, 36, 177, 202, 212, 155, 9, 116,
	128, 61, 189, 218, 152, 137, 240, 103, 59, 135, 133, 134, 102, 136, 217, 60,
	88, 104, 81, 241, 143, 138, 34, 153, 173, 219, 28, 190, 42, 62, 24, 129,
	186, 130, 78, 25, 75, 63, 198, 43, 94, 191, 66, 29, 166, 220, 160, 174,
	115, 154, 201, 35, 148, 139, 71, 144, 253, 242, 46, 82, 231, 105, 54, 89,
};

static const u8 ALOG_TABLE[512*2+1] = {
	1, 2, 4, 8, 16, 32, 64, 128, 95, 190, 35, 70, 140, 71, 142, 67,
	134, 83, 166, 19, 38, 76, 152, 111, 222, 227, 153, 109, 218, 235, 137, 77,
	154, 107, 214, 243, 185, 45, 90, 180, 55, 110, 220, 231, 145, 125, 250, 171,
	9, 18, 36, 72, 144, 127, 254, 163, 25, 50, 100, 200, 207, 193, 221, 229,
	149, 117, 234, 139, 73, 146, 123, 246, 179, 57, 114, 228, 151, 113, 226, 155,
	105, 210, 251, 169, 13, 26, 52, 104, 208, 255, 161, 29, 58, 116, 232, 143,
	65, 130, 91, 182, 51, 102, 204, 199, 209, 253, 165, 21, 42, 84, 168, 15,
	30, 60, 120, 240, 191, 33, 66, 132, 87, 174, 3, 6, 12, 24, 48, 96,
	192, 223, 225, 157, 101, 202, 203, 201, 205, 197, 213, 245, 181, 53, 106, 212,
	247, 177, 61, 122, 244, 183, 49, 98, 196, 215, 241, 189, 37, 74, 148, 119,
	238, 131, 89, 178, 59, 118, 236, 135, 81, 162, 27, 54, 108, 216, 239, 129,
	93, 186, 43, 86, 172, 7, 14, 28, 56, 112, 224, 159, 97, 194, 219, 233,
	141, 69, 138, 75, 150, 115, 230, 147, 121, 242, 187, 41, 82, 164, 23, 46,
	92, 184, 47, 94, 188, 39, 78, 156, 103, 206, 195, 217, 237, 133, 85, 170,
	11, 22, 44, 88, 176, 63, 126, 252, 167, 17, 34, 68, 136, 79, 158, 99,
	198, 211, 249, 173, 5, 10, 20, 40, 80, 160, 31, 62, 124, 248, 175, 1,
	2, 4, 8, 16, 32, 64, 128, 95, 190, 35, 70, 140, 71, 142, 67, 134,
	83, 166, 19, 38, 76, 152, 111, 222, 227, 153, 109, 218, 235, 137, 77, 154,
	107, 214, 243, 185, 45, 90, 180, 55, 110, 220, 231, 145, 125, 250, 171, 9,
	18, 36, 72, 144, 127, 254, 163, 25, 50, 100, 200, 207, 193, 221, 229, 149,
	117, 234, 139, 73, 146, 123, 246, 179, 57, 114, 228, 151, 113, 226, 155, 105,
	210, 251, 169, 13, 26, 52, 104, 208, 255, 161, 29, 58, 116, 232, 143, 65,
	130, 91, 182, 51, 102, 204, 199, 209, 253, 165, 21, 42, 84, 168, 15, 30,
	60, 120, 240, 191, 33, 66, 132, 87, 174, 3, 6, 12, 24, 48, 96, 192,
	223, 225, 157, 101, 202, 203, 201, 205, 197, 213, 245, 181, 53, 106, 212, 247,
	177, 61, 122, 244, 183, 49, 98, 196, 215, 241, 189, 37, 74, 148, 119, 238,
	131, 89, 178, 59, 118, 236, 135, 81, 162, 27, 54, 108, 216, 239, 129, 93,
	186, 43, 86, 172, 7, 14, 28, 56, 112, 224, 159, 97, 194, 219, 233, 141,
	69, 138, 75, 150, 115, 230, 147, 121, 242, 187, 41, 82, 164, 23, 46, 92,
	184, 47, 94, 188, 39, 78, 156, 103, 206, 195, 217, 237, 133, 85, 170, 11,
	22, 44, 88, 176, 63, 126, 252, 167, 17, 34, 68, 136, 79, 158, 99, 198,
	211, 249, 173, 5, 10, 20, 40, 80, 160, 31, 62, 124, 248, 175, 1, 0,
};

CAT_INLINE u8 Multiply(u8 x, u8 y)
{
	return ALOG_TABLE[LOG_TABLE[x] + LOG_TABLE[y]];
}

CAT_INLINE u8 Divide(u8 x, u8 y)
{
	// Precondition: y != 0
	return ALOG_TABLE[LOG_TABLE[x] + 255 - LOG_TABLE[y]];
}


GF256Matrix::GF256Matrix()
{
	_matrix = 0;
	_pivot = 0;
}

GF256Matrix::~GF256Matrix()
{
	Cleanup();
}

void GF256Matrix::Cleanup()
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

void GF256Matrix::Zero()
{
	memset(_matrix, 0, _pitch * _n);
}

void GF256Matrix::Identity()
{
	Zero();

	// For each row,
	u8 *row = _matrix;
	for (int column_i = 0; column_i < _n; ++column_i, row += _pitch)
	{
		row[column_i] = 1;
	}
}

void GF256Matrix::Fill()
{
	CatsChoice prng;
	prng.Initialize(_seed);

	u8 *row = _matrix;
	int words = _pitch * _n;
	while (words--)
	{
		*row++ = (u8)prng.Next();
	}
}

bool GF256Matrix::Triangle()
{
	CAT_IF_DEBUG(cout << endl << "---- Triangle ----" << endl << endl;)

	// Initialize pivot array
	for (u16 pivot_i = 0; pivot_i < _n; ++pivot_i)
		_pivot[pivot_i] = pivot_i;

	// For each pivot to determine,
	for (u16 pivot_i = 0; pivot_i < _n; ++pivot_i)
	{
		// Find pivot
		bool found = false;
		for (u16 pivot_j = pivot_i; pivot_j < _n; ++pivot_j)
		{
			// Determine if the row contains the bit we want
			u16 ge_row_j = _pivot[pivot_j];
			u8 *ge_row = _matrix + _pitch * ge_row_j;

			// If the bit was found,
			if (ge_row[pivot_j])
			{
				found = true;
				CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " found on row " << ge_row_j << endl;)

				// Swap out the pivot index for this one
				u16 temp = _pivot[pivot_i];
				_pivot[pivot_i] = _pivot[pivot_j];
				_pivot[pivot_j] = temp;

				// For each remaining unused row,
				for (u16 pivot_k = pivot_j + 1; pivot_k < _n; ++pivot_k)
				{
					// Determine if the row contains the bit we want
					u16 ge_row_k = _pivot[pivot_k];
					u8 *rem_row = _matrix + _pitch * ge_row_k;

					// If the bit was found,
					if (rem_row[pivot_j])
					{
						for (int ii = 0; ii < _pitch; ++ii)
						{
							rem_row[ii] ^= Multiply(Divide(ge_row[ii], ge_row[pivot_j]), rem_row[ii]);
						}
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
	}

	return true;
}

void GF256Matrix::Diagonal()
{
	CAT_IF_DEBUG(cout << endl << "---- Diagonal ----" << endl << endl;)

	// For each pivot row from last to first,
	int pivot_i = _n - 1;
	for (; pivot_i >= 0; --pivot_i)
	{
		u16 pivot_ge_row_i = _pivot[pivot_i];

		CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " solving row " << pivot_ge_row_i << endl;)

		// For each pivot row above it,
		for (int above_i = pivot_i - 1; above_i >= 0; --above_i)
		{
			u16 ge_above_row_i = _pivot[above_i];

			// If bit is set in that row,
			u8 *above_row = _matrix + _pitch * ge_above_row_i;
			if (above_row[pivot_i])
			{
				// Eliminate here...
				above_row[pivot_i] = 0;
			}
		}
	}
}

bool GF256Matrix::Initialize(int n)
{
	Cleanup();

	CAT_IF_DEBUG(cout << endl << "GF256Matrix.Initialize: Seed = " << _seed << " n = " << n << endl << endl;)

	_n = n;
	_pitch = n;
	int words = _pitch * n;
	_matrix = new u8[words];
	if (!_matrix) return false;
	_pivot = new u16[n];
	if (!_pivot) return false;
	_seed = 0;

	return true;
}

void GF256Matrix::Print()
{
	int rows = _n;
	int cols = _n;

	cout << endl << "GF256Matrix is " << rows << " x " << cols << " (seed " << _seed << "):" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			cout << hex << setfill('0') << setw(2) << (int)_matrix[_pitch * ii + jj] << " ";
		}
		cout << endl;
	}

	cout << dec << endl;
}
