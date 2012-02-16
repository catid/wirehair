/*
	Copyright (c) 2012 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of LibCat nor the names of its contributors may be used
	  to endorse or promote products derived from this software without
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

#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "memxor.hpp"
using namespace cat;
using namespace wirehair;

static int seed = 28;

CAT_INLINE static u32 GetGeneratorSeed(int block_count)
{
	// TODO: Needs to be simulated (2)
	return seed++;
}


/*
	TODO:

	1. Implement peeling matrix inverse compression algorithm
	2. Fix bugs in non-windowed multiplication compression algorithm
	3. Fix bugs in windowed multiplication compression algorithm
	4. Implement decoder
	5. Simulate and tune the added check block count for each N
	6. Tune the dense matrix rows for faster multiplication
	7. Implement implicit zeroed data blocks so that table won't be huge
	8. Create a table to store the seed values for each N
	9. Simulate and tune the generator seed for each N
	10. Implement an asymmetric version too

	Benchmark it!
	Document it!
	Ship it!
	???
	Profit!
*/


// Switches:
#define CAT_STEW_HYPERDYNAMIC_PLATTONIC_ITERATOR /* Use Stew's more efficient column iterator */
//#define CAT_DUMP_ROWOP_COUNTERS /* Dump row operations counters to console */
//#define CAT_ENCODER_COPY_FIRST_N /* Copy the first N rows from the input (much faster) */
#define CAT_INVERSE_THRESHOLD 1 /* Block count where peeling inverse version starts getting used */
#define CAT_WINDOW_THRESHOLD 16 /* Compression row count when 4-bit window is employed */


#if defined(CAT_DEBUG) || defined(CAT_DUMP_ROWOP_COUNTERS)
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#endif // CAT_DEBUG


//// Utility: 16-bit Integer Square Root function

/*
	Based on code from http://www.azillionmonkeys.com/qed/sqroot.html

		Contributors include Arne Steinarson for the basic approximation idea, 
		Dann Corbit and Mathew Hendry for the first cut at the algorithm, 
		Lawrence Kirby for the rearrangement, improvments and range optimization
		and Paul Hsieh for the round-then-adjust idea.

	I tried this out, stdlib sqrt() and a few variations on Newton-Raphson
	and determined this one is, by far, the fastest.  I adapted it to 16-bit
	input for additional performance and tweaked the operations to work best
	with the MSVC optimizer, which turned out to be very sensitive to the
	way that the code is written.
*/
static const u8 SQQ_TABLE[] = {
	0,  16,  22,  27,  32,  35,  39,  42,  45,  48,  50,  53,  55,  57,
	59,  61,  64,  65,  67,  69,  71,  73,  75,  76,  78,  80,  81,  83,
	84,  86,  87,  89,  90,  91,  93,  94,  96,  97,  98,  99, 101, 102,
	103, 104, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118,
	119, 120, 121, 122, 123, 124, 125, 126, 128, 128, 129, 130, 131, 132,
	133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 144, 145,
	146, 147, 148, 149, 150, 150, 151, 152, 153, 154, 155, 155, 156, 157,
	158, 159, 160, 160, 161, 162, 163, 163, 164, 165, 166, 167, 167, 168,
	169, 170, 170, 171, 172, 173, 173, 174, 175, 176, 176, 177, 178, 178,
	179, 180, 181, 181, 182, 183, 183, 184, 185, 185, 186, 187, 187, 188,
	189, 189, 190, 191, 192, 192, 193, 193, 194, 195, 195, 196, 197, 197,
	198, 199, 199, 200, 201, 201, 202, 203, 203, 204, 204, 205, 206, 206,
	207, 208, 208, 209, 209, 210, 211, 211, 212, 212, 213, 214, 214, 215,
	215, 216, 217, 217, 218, 218, 219, 219, 220, 221, 221, 222, 222, 223,
	224, 224, 225, 225, 226, 226, 227, 227, 228, 229, 229, 230, 230, 231,
	231, 232, 232, 233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238,
	239, 240, 240, 241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246,
	246, 247, 247, 248, 248, 249, 249, 250, 250, 251, 251, 252, 252, 253,
	253, 254, 254, 255
};

static u16 cat_fred_sqrt16(u16 x)
{
	u16 r;

	if (x >= 0x100)
	{
		if (x >= 0x1000)
		{
			if (x >= 0x4000)
				r = SQQ_TABLE[x >> 8] + 1;
			else
				r = (SQQ_TABLE[x >> 6] >> 1) + 1;
		}
		else
		{
			if (x >= 0x400)
				r = (SQQ_TABLE[x >> 4] >> 2) + 1;
			else
				r = (SQQ_TABLE[x >> 2] >> 3) + 1;
		}
	}
	else
	{
		return SQQ_TABLE[x] >> 4;
	}

	// Correct rounding if necessary (compiler optimizes this form better)
	r -= (r * r > x);

	return r;
}


//// Utility: 16-bit Truncated Sieve of Eratosthenes Next Prime function

/*
	It uses trial division up to the square root of the number to test.
	Uses a truncated sieve table to pick the next number to try, which
	avoids small factors 2, 3, 5, 7.  This can be considered a more
	involved version of incrementing by 2 instead of 1.  It takes about
	25% less time on average than the approach that just increments by 2.

	Because of the tabular increment this is a hybrid approach.  The sieve
	would just use a very large table, but I wanted to limit the size of
	the table to something fairly small and reasonable.  210 bytes for the
	sieve table.  102 bytes for the primes list.

	It also calculates the integer square root faster than cmath sqrt()
	and uses multiplication to update the square root instead for speed.
*/
const int SIEVE_TABLE_SIZE = 2*3*5*7;
static const u8 SIEVE_TABLE[SIEVE_TABLE_SIZE] = {
	1, 0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	7, 6, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
};

static const u16 PRIMES_UNDER_256[] = {
	11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
	67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
	131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
	193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 0x7fff
};

static u16 NextPrime16(u16 n)
{
	// Handle small n
	switch (n)
	{
	case 0:
	case 1:	return 1;
	case 2:	return 2;
	case 3:	return 3;
	case 4:
	case 5:	return 5;
	case 6:
	case 7:	return 7;
	}

	// Choose first n from table
	int offset = n % SIEVE_TABLE_SIZE;
	u32 next = SIEVE_TABLE[offset];
	offset += next + 1;
	n += next;

	// Initialize p_max to sqrt(n)
	int p_max = cat_fred_sqrt16(n);

	// For each number to try,
	for (;;)
	{
		// For each prime to test up to square root,
		const u16 *prime = PRIMES_UNDER_256;
		for (;;)
		{
			// If the next prime is above p_max we are done!
			int p = *prime;
			if (p > p_max)
				return n;

			// If composite, try next n
			if (n % p == 0)
				break;

			// Try next prime
			++prime;
		}

		// Use table to choose next trial number
		if (offset >= SIEVE_TABLE_SIZE) offset -= SIEVE_TABLE_SIZE;
		u32 next = SIEVE_TABLE[offset];
		offset += next + 1;
		n += next + 1;

		// Derivative square root iteration of p_max
		if (p_max * p_max < n)
			++p_max;
	}

	return n;
}


//// MatrixGenerator

static const u32 WEIGHT_DIST[] = {
	0, 5243, 529531, 704294, 791675, 844104, 879057, 904023,
	922747, 937311, 948962, 958494, 966438, 973160, 978921,
	983914, 988283, 992138, 995565, 998631, 1001391, 1003887,
	1006157, 1008229, 1010129, 1011876, 1013490, 1014983,
	1016370, 1017662, 1048576
};

static u16 GenerateWeight(u32 rv, u16 max_weight)
{
	rv &= 0xfffff;

	u16 ii;
	for (ii = 1; rv >= WEIGHT_DIST[ii]; ++ii);

	return ii > max_weight ? max_weight : ii;
}

CAT_INLINE static int GetCheckBlockCount(int block_count)
{
	// TODO: Needs to be simulated (1)
	return cat_fred_sqrt16(block_count) + 1;
}

#if defined(CAT_STEW_HYPERDYNAMIC_PLATTONIC_ITERATOR)
/*
	This is Stewart Platt's excellent loop-less iterator optimization.
	His common cases all require no additional modulus operation, which
	makes it faster than the rare case that I designed.
*/
CAT_INLINE static void IterateNextColumn(u16 &x, u16 b, u16 p, u16 a)
{
	x = (x + a) % p;

	if (x >= b)
	{
		u16 distance = p - x;

		if (a >= distance)
			x = a - distance;
		else // the rare case:
			x = (((u32)a << 16) - distance) % a;
	}
}
#else
CAT_INLINE static void IterateNextColumn(u16 &x, u16 b, u16 p, u16 a)
{
	do x = (x + a) % p;
	while (x >= b);
}
#endif // STEW_HYPERDYNAMIC_PLATTONIC_ITERATOR


#if defined(CAT_DUMP_ROWOP_COUNTERS)
#define CAT_IF_ROWOP(x) x
#else
#define CAT_IF_ROWOP(x)
#endif


//// Encoder

#pragma pack(push)
#pragma pack(1)
struct Encoder::PeelRow
{
	u16 next;				// Linkage in row list

	// Peeling matrix: Column generator
	u16 peel_weight, peel_a, peel_x0;

	// Mixing matrix: Column generator
	u16 mix_weight, mix_a, mix_x0;

	// Peeling state
	u16 unmarked_count;		// Count of columns that have not been marked yet
	union
	{
		// During peeling:
		struct
		{
			u16 unmarked[2];		// Final two unmarked column indices
		};

		// After peeling:
		struct
		{
			u16 peel_column;
			u8 is_copied;
		};
	};
};
#pragma pack(pop)

// During generator matrix precomputation, tune this to be as small as possible and still succeed
static const int ENCODER_REF_LIST_MAX = 64;

enum MarkTypes
{
	MARK_TODO,
	MARK_PEEL,
	MARK_DEFER
};

#pragma pack(push)
#pragma pack(1)
struct Encoder::PeelColumn
{
	u16 next;			// Linkage in column list

	union
	{
		u16 w2_refs;	// Number of weight-2 rows containing this column
		u16 peel_row;	// Row that solves the column
		u16 ge_column;	// Column that a deferred column is mapped to
	};

	u16 row_count;		// Number of rows containing this column
	u16 rows[ENCODER_REF_LIST_MAX];
	u8 mark;			// One of the MarkTypes enumeration
};
#pragma pack(pop)


//// (1) Peeling:

/*
	Peel() and PeelAvalanche() are split up into two functions because I found
	that the PeelAvalanche() function can be reused later during GreedyPeeling().

		Until N rows are received, the peeling algorithm is executed:

		Columns have 3 states:

		(1) Peeled - Solved by a row during peeling process.
		(2) Deferred - Will be solved by a row during Gaussian Elimination.
		(3) Unmarked - Still deciding.

		Initially all columns are unmarked.

		As a row comes in, the count of columns that are Unmarked is calculated.
	If that count is 1, then the column is marked as Peeled and is solved by
	the received row.  Peeling then goes through all other rows that reference
	that peeled column, reducing their count by 1, potentially causing other
	columns to be marked as Peeled.  This "peeling avalanche" is desired.
*/

bool Encoder::PeelSetup()
{
	CAT_IF_DEBUG(cout << endl << "---- PeelSetup ----" << endl << endl;)

	CAT_IF_DEBUG(cout << "Block Count = " << _block_count << endl;)

	// Allocate space for row data
	_peel_rows = new PeelRow[_block_count];
	if (!_peel_rows) return false;

	// Allocate space for column data
	_peel_cols = new PeelColumn[_block_count];
	if (!_peel_cols) return false;

	// Initialize columns
	for (int ii = 0; ii < _block_count; ++ii)
	{
		_peel_cols[ii].row_count = 0;
		_peel_cols[ii].w2_refs = 0;
		_peel_cols[ii].mark = MARK_TODO;
	}

	// Initialize lists
	_peel_head_rows = _peel_tail_rows = LIST_TERM;
	_defer_head_rows = LIST_TERM;

	// Set inverse method flag based on block count
	_use_inverse_method = (_block_count >= CAT_INVERSE_THRESHOLD);

	// Calculate default mix weight
	u16 mix_weight = 3;
	if (mix_weight >= _added_count)
		mix_weight = _added_count - 1;

	// Generate peeling row data
	for (u16 row_i = 0; row_i < _block_count; ++row_i)
	{
		PeelRow *row = &_peel_rows[row_i];

		// Initialize PRNG
		CatsChoice prng;
		prng.Initialize(row_i, _g_seed);

		// Generate peeling matrix row parameters
		row->peel_weight = GenerateWeight(prng.Next(), _block_count - 1);
		u32 rv = prng.Next();
		row->peel_a = ((u16)rv % (_block_count - 1)) + 1;
		row->peel_x0 = (u16)(rv >> 16) % _block_count;

		// Generate mixing matrix row parameters
		row->mix_weight = mix_weight;
		rv = prng.Next();
		row->mix_a = ((u16)rv % (_added_count - 1)) + 1;
		row->mix_x0 = (u16)(rv >> 16) % _added_count;

		CAT_IF_DEBUG(cout << "   Row " << row_i << " of weight " << row->peel_weight << " [a=" << row->peel_a << "] : ";)

		// Iterate columns in peeling matrix
		u16 weight = row->peel_weight;
		u16 column_i = row->peel_x0;
		u16 a = row->peel_a;
		u16 unmarked_count = 0;
		u16 unmarked[2];
		for (;;)
		{
			CAT_IF_DEBUG(cout << column_i << " ";)

			PeelColumn *col = &_peel_cols[column_i];

			// Add row reference to column
			if (col->row_count >= ENCODER_REF_LIST_MAX)
			{
				CAT_IF_DEBUG(cout << "PeelSetup: Failure!  Ran out of space for row references.  ENCODER_REF_LIST_MAX must be increased!" << endl;)
				return false;
			}
			col->rows[col->row_count++] = row_i;

			// If column is unmarked,
			if (col->mark == MARK_TODO)
				unmarked[unmarked_count++ & 1] = column_i;

			if (--weight <= 0) break;

			IterateNextColumn(column_i, _block_count, _block_next_prime, a);
		}
		CAT_IF_DEBUG(cout << endl;)

		// Initialize row state
		row->unmarked_count = unmarked_count;

		switch (unmarked_count)
		{
		case 0:
			// Link at head of defer list
			row->next = _defer_head_rows;
			_defer_head_rows = row_i;
			break;

		case 1:
			// Solve only unmarked column with this row
			Peel(row_i, row, unmarked[0]);
			break;

		case 2:
			// Remember which two columns were unmarked
			row->unmarked[0] = unmarked[0];
			row->unmarked[1] = unmarked[1];

			// Increment weight-2 reference count for unmarked columns
			_peel_cols[unmarked[0]].w2_refs++;
			_peel_cols[unmarked[1]].w2_refs++;
			break;
		}
	}

	return true;
}

void Encoder::PeelAvalanche(u16 column_i, PeelColumn *column)
{
	// Walk list of peeled rows referenced by this newly solved column
	u16 ref_row_count = column->row_count;
	u16 *ref_rows = column->rows;
	while (ref_row_count--)
	{
		// Update unmarked row count for this referenced row
		u16 ref_row_i = *ref_rows++;
		PeelRow *ref_row = &_peel_rows[ref_row_i];
		u16 unmarked_count = --ref_row->unmarked_count;

		// If row may be solving a column now,
		if (unmarked_count == 1)
		{
			// Find other column
			u16 new_column_i = ref_row->unmarked[0];
			if (new_column_i == column_i)
				new_column_i = ref_row->unmarked[1];

			/*
				Rows that are to be deferred will either end up
				here or below where it handles the case of there
				being no columns unmarked in a row.
			*/

			// If column is already solved,
			if (_peel_cols[new_column_i].mark == MARK_TODO)
				Peel(ref_row_i, ref_row, new_column_i);
			else
			{
				CAT_IF_DEBUG(cout << "PeelAvalanche: Deferred(1) with column " << column_i << " at row " << ref_row_i << endl;)

				// Link at head of defer list
				ref_row->next = _defer_head_rows;
				_defer_head_rows = ref_row_i;
			}
		}
		else if (unmarked_count == 2)
		{
			// Regenerate the row columns to discover which are unmarked
			u16 ref_weight = ref_row->peel_weight;
			u16 ref_column_i = ref_row->peel_x0;
			u16 ref_a = ref_row->peel_a;
			u16 unmarked_count = 0;
			for (;;)
			{
				PeelColumn *ref_col = &_peel_cols[ref_column_i];

				// If column is unmarked,
				if (ref_col->mark == MARK_TODO)
				{
					// Store the two unmarked columns in the row
					ref_row->unmarked[unmarked_count++] = ref_column_i;

					// Increment weight-2 reference count (cannot hurt even if not true)
					ref_col->w2_refs++;
				}

				if (--ref_weight <= 0) break;

				IterateNextColumn(ref_column_i, _block_count, _block_next_prime, ref_a);
			}

			/*
				This is a little subtle, but sometimes the avalanche will
				happen here, and sometimes a row will be marked deferred.
			*/

			if (unmarked_count <= 1)
			{
				// Insure that this row won't be processed further during this recursion
				ref_row->unmarked_count = 0;

				// If row is to be deferred,
				if (unmarked_count == 1)
					Peel(ref_row_i, ref_row, ref_row->unmarked[0]);
				else
				{
					CAT_IF_DEBUG(cout << "PeelAvalanche: Deferred(2) with column " << column_i << " at row " << ref_row_i << endl;)

					// Link at head of defer list
					ref_row->next = _defer_head_rows;
					_defer_head_rows = ref_row_i;
				}
			}
		}
	}
}

void Encoder::Peel(u16 row_i, PeelRow *row, u16 column_i)
{
	CAT_IF_DEBUG(cout << "Peel: Solved column " << column_i << " with row " << row_i << endl;)

	PeelColumn *column = &_peel_cols[column_i];

	// Mark this column as solved
	column->mark = MARK_PEEL;

	// Remember which column it solves
	row->peel_column = column_i;

	// If using the inverse compression algorithm,
	if (_use_inverse_method)
	{
		// Link to back of the peeled list
		if (_peel_tail_rows != LIST_TERM)
			_peel_rows[_peel_tail_rows].next = row_i;
		else
			_peel_head_rows = row_i;
		row->next = LIST_TERM;
		_peel_tail_rows = row_i;

		// Indicate that this row hasn't been copied yet
		row->is_copied = 0;
	}
	else
	{
		// Link to front of the peeled list
		row->next = _peel_head_rows;
		_peel_head_rows = row_i;
	}
	// NOTE: In both algorithms, after compression the list will be in forward rather than reverse order.

	// Attempt to avalanche and solve other columns
	PeelAvalanche(column_i, column);

	// Remember which row solves the column, after done with rows list
	column->peel_row = row_i;
}

/*
		After the opportunistic peeling solver has completed, no columns have
	been deferred to Gaussian elimination yet.  Greedy peeling will then take
	over and start choosing columns to defer.  It is greedy in that the selection
	of which columns to defer is based on a greedy initial approximation of the
	best column to choose, rather than the best one that could be chosen.

		In this algorithm, the column that will cause the largest immediate
	avalanche of peeling solutions is the one that is selected.  If there is
	a tie between two or more columns based on just that criterion then, of the
	columns that tied, the one that affects the most rows is selected.

		In practice with a well designed peeling matrix, only about sqrt(N)
	columns must be deferred to Gaussian elimination using this greedy approach.
*/

void Encoder::GreedyPeeling()
{
	CAT_IF_DEBUG(cout << endl << "---- GreedyPeeling ----" << endl << endl;)

	// Initialize list
	_defer_head_columns = LIST_TERM;
	_defer_count = 0;

	// Until all columns are marked,
	for (;;)
	{
		u16 best_column_i = LIST_TERM;
		u16 best_w2_refs = 0, best_row_count = 0;

		// For each column,
		for (u16 column_i = 0; column_i < _block_count; ++column_i)
		{
			PeelColumn *column = &_peel_cols[column_i];
			u16 w2_refs = column->w2_refs;

			// If column is not marked yet,
			if (column->mark == MARK_TODO)
			{
				// And if it may have the most weight-2 references
				if (w2_refs >= best_w2_refs)
				{
					// Or if it has the largest row references overall,
					if (w2_refs > best_w2_refs || column->row_count >= best_row_count)
					{
						// Use that one
						best_column_i = column_i;
						best_w2_refs = w2_refs;
						best_row_count = column->row_count;
					}
				}
			}
		}

		// If done peeling,
		if (best_column_i == LIST_TERM)
			break;

		// Mark column as deferred
		PeelColumn *best_column = &_peel_cols[best_column_i];
		best_column->mark = MARK_DEFER;
		++_defer_count;

		// Add at head of deferred list
		best_column->next = _defer_head_columns;
		_defer_head_columns = best_column_i;

		CAT_IF_DEBUG(cout << "Deferred column " << best_column_i << " for Gaussian elimination, which had " << best_column->w2_refs << " weight-2 row references" << endl;)

		// Peel resuming from where this column left off
		PeelAvalanche(best_column_i, best_column);
	}
}

/*
		After the peeling solver has completed, only Deferred columns remain
	to be solved.  Conceptually the matrix can be re-ordered in the order of
	solution so that the matrix resembles this example:

		+-----+---------+-----+
		|   1 | 1       | 721 | <-- Peeled row 1
		| 1   |   1     | 518 | <-- Peeled row 2
		| 1   | 1   1   | 934 | <-- Peeled row 3
		|   1 | 1   1 1 | 275 | <-- Peeled row 4
		+-----+---------+-----+
		| 1   | 1 1   1 | 123 | <-- Deferred rows
		|   1 |   1 1 1 | 207 | <-- Deferred rows
		+-----+---------+-----+
			^       ^       ^---- Row value (unmodified so far, ie. 1500 bytes)
			|       \------------ Peeled columns
			\-------------------- Deferred columns

		Re-ordering the actual matrix is not required, but the lower-triangular
	form of the peeled matrix is apparent in the diagram above.
*/


//// (2) Compression:

/*
	Surprisingly, this is by far the most complex step of the algorithm.

	At this point the generator matrix has been re-organized
	into peeled and deferred rows and columns:

		+-----------------------+
		| p p p p p | B B | C C |
		+-----------------------+
					X
		+-----------+-----+-----+
		| 5 2 6 1 3 | 4 0 | x y |
		+-----------+-----+-----+---+   +---+
		| 1         | 0 0 | 1 0 | 0 |   | p |
		| 0 1       | 0 0 | 0 1 | 5 |   | p |
		| 1 0 1     | 1 0 | 1 0 | 3 |   | p |
		| 0 1 0 1   | 0 0 | 0 1 | 4 |   | p |
		| 0 1 0 1 1 | 1 1 | 1 0 | 1 |   | p |
		+-----------+-----+-----+---| = |---|
		| 0 1 0 1 1 | 1 1 | 0 1 | 2 |   | A |
		| 0 1 0 1 0 | 0 1 | 1 0 | 6 |   | A |
		+-----------+-----+-----+---|   |---|
		| 1 1 0 1 1 | 1 1 | 1 0 | x |   | 0 |
		| 1 0 1 1 0 | 1 0 | 0 1 | y |   | 0 |
		+-----------+-----+-----+---+   +---+

	P = Peeled rows/columns (re-ordered)
	B = Columns deferred for GE
	C = Mixing columns always deferred for GE
	A = Sparse rows from peeling matrix deferred for GE
	0 = Dense rows from dense matrix deferred for GE
*/

/*
	Compression step objective:

		The goal of compression is to eliminate the peeled columns from
	the deferred rows of the matrix, so that Gaussian elimination can be
	applied to the deferred columns of the matrix.  The GE matrix will
	be a small square matrix with dimensions proportional to SQRT(N),
	which is already roughly upper triangular (see notes below).

		I see two ways to do this step:

	(1) Multiply the peeled matrix into the deferred rows to eliminate them.

		This approach is called Multiplication-Based Matrix Compression here.
		The multiplication step is O(n^^1.5) for all rows, but doesn't
			require the peeled matrix inversion.
		Matrix uses less memory because the longer dimension is bitwise.
		Can use a windowed approach to multiplication for better performance.

	(2) Invert the peeled matrix and then multiply into deferred rows.

		This approach is called Inversion-Based Matrix Compression here.
		The peeled matrix inversion scales O(n).
		The multiplication step is O(sqrt(n)) for deferred rows, and
			O(n^^1.5) for dense rows.  Windowing helps for dense rows.

	Both approaches require repeating matrix multiplication after GE Triangle
	finishes, so that no temporary space is needed for block values.

	Both approaches cost about the same for the dense rows of the matrix, because
	about half of the bits are set in both cases.  The multiplication approach,
	however, destroys structure in the dense rows if it was there, so it may be
	better to use the inversion approach if the dense rows have structure at lower
	n values than if the dense rows have no structure.

	Multiplication is faster than Inversion for small sizes because it uses much
	less memory and because the row op count is comparable or lower.
*/

/*
	Inversion-Based Matrix Compression:

		Since the re-ordered matrix above is in lower triangular form,
	and since the weight of each row is limited to a constant, the cost
	of inverting the peeled matrix is O(n).  Inverting the peeled matrix
	will cause the mix and deferred columns to add up and become dense.
	These columns are stored in the same order as in the GE matrix, in
	a long vertical matrix with N rows.

		After inverting the peeling matrix, causing the mixing columns
	to become dense, the peeling matrix will be the identity matrix.
	Peeled column output blocks can be used to store the temporary blocks
	generated by this inversion.  These temporary blocks will be overwritten
	later during Substitution.  It then becomes easy to eliminate the
	peeled columns of the deferred rows.  These are then eliminated in
	O(sqrt(n)) time since there are k*sqrt(n) of them and each is sparse.

		Initially, only the dense mix and deferred columns are useful for
	forming the GE matrix.  After the GE matrix is inverted, the steps used
	to eliminate the peeled columns of the deferred rows can be followed again
	to generate the columns solved by GE in the correct output block locations.
	Note that it would be possible to generate these columns earlier but the
	output block locations are not known until GE Triangle() completes, so they
	would need to be shuffled into the right locations, which would require
	more memory and/or memcopies.
*/

bool Encoder::InvCompressSetup()
{
	CAT_IF_DEBUG(cout << endl << "---- InvCompressSetup ----" << endl << endl;)

	if (!InvCompressAllocate())
		return false;

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	return true;
}

void Encoder::InvPeelDiagonal()
{
	CAT_IF_DEBUG(cout << endl << "---- InvPeelDiagonal ----" << endl << endl;)

	/*
		To produce the GE matrix and peeled row values:

			(1) Set deferred column bits
			(2) Set mixing columns for each deferred row and mark them as deferred
			(3) Diagonalize the peeled rows

		This function optimizes the block value generation by combining the first
		memcpy and memxor operations together into a three-way memxor if possible,
		using the is_copied row member.
	*/

	// (1) Set deferred column bits:

	CAT_IF_DEBUG(cout << "(1) Set deferred column bits:" << endl;)

	// For each deferred column,
	PeelColumn *column;
	for (u16 ge_column_i = 0, defer_i = _defer_head_columns; defer_i != LIST_TERM; defer_i = column->next, ++ge_column_i)
	{
		column = &_peel_cols[defer_i];

		CAT_IF_DEBUG(cout << "  GE column " << ge_column_i << " mapped to matrix column " << defer_i << " :";)

		// Set bit for each row affected by this deferred column
		u64 *matrix_row_offset = _ge_compress_matrix + (ge_column_i >> 6);
		u64 ge_mask = (u64)1 << (ge_column_i & 63);
		u16 count = column->row_count;
		u16 *ref_row = column->rows;
		while (count--)
		{
			u16 row_i = *ref_row++;

			CAT_IF_DEBUG(cout << " " << row_i;)

			matrix_row_offset[_ge_compress_pitch * row_i] |= ge_mask;
		}

		CAT_IF_DEBUG(cout << endl;)

		// Set column map for this GE column
		_ge_col_map[ge_column_i] = defer_i;

		// Set reverse mapping also
		column->ge_column = ge_column_i;
	}

	// Set column map for each mix column
	for (u16 added_i = 0; added_i < _added_count; ++added_i)
	{
		u16 ge_column_i = _defer_count + added_i;
		u16 column_i = _block_count + added_i;

		CAT_IF_DEBUG(cout << "  GE column(mix) " << ge_column_i << " mapped to matrix column " << column_i << endl;)

		_ge_col_map[ge_column_i] = column_i;
	}

	CAT_IF_DEBUG(cout << "After filling in deferred columns:" << endl;)
	CAT_IF_DEBUG(InvPrintGECompressMatrix();)

	// (2) Set mixing columns for each deferred row and mark them as deferred:

	CAT_IF_DEBUG(cout << "(2) Set mixing columns for each deferred row and mark them as deferred:" << endl;)

	// For each deferred row,
	PeelRow *row;
	for (u16 defer_row_i = _defer_head_rows; defer_row_i != LIST_TERM; defer_row_i = row->next)
	{
		row = &_peel_rows[defer_row_i];

		CAT_IF_DEBUG(cout << "  Deferred row " << defer_row_i << " set mix columns :";)

		// Mark it as deferred for the following loop
		row->peel_column = LIST_TERM;

		// Generate mixing columns for this row
		u16 weight = row->mix_weight;
		u16 a = row->mix_a;
		u16 x = row->mix_x0;
		u64 *ge_row = _ge_compress_matrix + _ge_compress_pitch * defer_row_i;
		for (;;)
		{
			// Flip bit for each mixing column
			u16 ge_column_i = _defer_count + x;
			u64 ge_mask = (u64)1 << (ge_column_i & 63);
			ge_row[ge_column_i >> 6] ^= ge_mask;

			CAT_IF_DEBUG(cout << " " << ge_column_i;)

			if (--weight <= 0) break;

			IterateNextColumn(x, _added_count, _added_next_prime, a);
		}

		CAT_IF_DEBUG(cout << endl;)
	}

	CAT_IF_DEBUG(cout << "After filling in mixing columns for deferred rows:" << endl;)
	CAT_IF_DEBUG(InvPrintGECompressMatrix();)

	// (3) Diagonalize the peeled rows:

	CAT_IF_DEBUG(cout << "(3) Diagonalize the peeled rows:" << endl;)

	// For each peeled row in forward solution order,
	for (u16 peel_row_i = _peel_head_rows; peel_row_i != LIST_TERM; peel_row_i = row->next)
	{
		row = &_peel_rows[peel_row_i];

		// Lookup peeling results
		u16 peel_column_i = row->peel_column;
		u64 *ge_row = _ge_compress_matrix + _ge_compress_pitch * peel_row_i;

		CAT_IF_DEBUG(cout << "  Peeled row " << peel_row_i << " for peeled column " << peel_column_i << " :";)

		// Generate mixing columns for this row
		u16 weight = row->mix_weight;
		u16 a = row->mix_a;
		u16 x = row->mix_x0;
		for (;;)
		{
			// Flip bit for each mixing column
			u16 ge_column_i = _defer_count + x;
			u64 ge_mask = (u64)1 << (ge_column_i & 63);
			ge_row[ge_column_i >> 6] ^= ge_mask;

			CAT_IF_DEBUG(cout << " " << ge_column_i;)

			if (--weight <= 0) break;

			IterateNextColumn(x, _added_count, _added_next_prime, a);
		}

		CAT_IF_DEBUG(cout << endl;)

		// Lookup output block
		u8 *temp_block_src = _check_blocks + _block_bytes * peel_column_i;

		// If row has not been copied yet,
		if (!row->is_copied)
		{
			// Copy it directly to the output symbol
			const u8 *block_src = _message_blocks + _block_bytes * peel_row_i;
			if (peel_row_i != _block_count - 1)
				memcpy(temp_block_src, block_src, _block_bytes);
			else
			{
				memcpy(temp_block_src, block_src, _final_bytes);
				memset(temp_block_src + _final_bytes, 0, _block_bytes - _final_bytes);
			}

			CAT_IF_DEBUG(cout << "  -- Copied from " << peel_row_i << " because has not been copied yet.  Output block = " << (int)temp_block_src[0] << endl;)

			// NOTE: Do not need to set is_copied here because no further rows reference this one
		}

		// For each row that references this one,
		column = &_peel_cols[peel_column_i];
		u16 count = column->row_count;
		u16 *ref_row = column->rows;
		while (count--)
		{
			u16 ref_row_i = *ref_row++;

			// Skip this row
			if (ref_row_i == peel_row_i) continue;

			CAT_IF_DEBUG(cout << "  ++ Adding to referencing row " << ref_row_i << endl;)

			// Add GE row to referencing GE row
			u64 *ge_ref_row = _ge_compress_matrix + _ge_compress_pitch * ref_row_i;
			for (int ii = 0; ii < _ge_compress_pitch; ++ii)
				ge_ref_row[ii] ^= ge_row[ii];

			// If row is peeled,
			PeelRow *ref_row = &_peel_rows[ref_row_i];
			u16 ref_column_i = ref_row->peel_column;
			if (ref_column_i != LIST_TERM)
			{
				// Generate temporary row block value:
				u8 *temp_block_dest = _check_blocks + _block_bytes * ref_column_i;

				// If referencing row is already copied to the check blocks,
				if (ref_row->is_copied)
				{
					// Add this row block value to it
					memxor(temp_block_dest, temp_block_src, _block_bytes);
				}
				else
				{
					// Add this row block value with message block to it (optimization)
					const u8 *block_src = _message_blocks + _block_bytes * ref_row_i;
					if (ref_row_i != _block_count - 1)
						memxor(temp_block_dest, temp_block_src, block_src, _block_bytes);
					else
					{
						memxor(temp_block_dest, temp_block_src, block_src, _final_bytes);
						memcpy(temp_block_dest + _final_bytes, temp_block_src, _block_bytes - _final_bytes);
					}

					ref_row->is_copied = 1;
				}
			} // end if referencing row is peeled
		} // next referencing row
	} // next peeled row

	CAT_IF_DEBUG(cout << "After peeling matrix inversion:" << endl;)
	CAT_IF_DEBUG(InvPrintGECompressMatrix();)
}

void Encoder::InvCopyDeferredRows()
{
	CAT_IF_DEBUG(cout << endl << "---- InvCopyDeferredRows ----" << endl << endl;)

	// For each deferred row,
	u64 *ge_row = _ge_matrix + _ge_pitch * _added_count;
	for (u16 ge_row_i = 0, defer_row_i = _defer_head_rows; defer_row_i != LIST_TERM;
		defer_row_i = _peel_rows[defer_row_i].next, ge_row += _ge_pitch, ++ge_row_i)
	{
		CAT_IF_DEBUG(cout << "Peeled row " << defer_row_i << " for GE row " << ge_row_i << endl;)

		// Copy compress row to GE row
		u64 *compress_row = _ge_compress_matrix + _ge_compress_pitch * defer_row_i;
		memcpy(ge_row, compress_row, _ge_compress_pitch * sizeof(u64));

		// Set row map for this deferred row
		_ge_row_map[ge_row_i] = defer_row_i;
	}

	CAT_IF_DEBUG(cout << "After copying deferred rows:" << endl;)
	CAT_IF_DEBUG(PrintGEMatrix();)
}

void Encoder::InvMultiplyDenseRows()
{
	CAT_IF_DEBUG(cout << endl << "---- InvMultiplyDenseRows ----" << endl << endl;)

	/*
		NOTE: This does not need to generate the same dense rows that
		the multiply-compression algorithm generates because both the
		encoder and decoder will make the same choice on which to use.

		TODO: Use 4-bit window optimization here also
	*/

	// For each dense row of the GE matrix,
	u64 *ge_dest_row = _ge_matrix;
	for (u16 ge_row_i = 0; ge_row_i < _added_count; ++ge_row_i, ge_dest_row += _ge_pitch)
	{
		CAT_IF_DEBUG(cout << "GE row " << ge_row_i << ":";)

		// Initialize the row to zeroes
		for (int ii = 0; ii < _ge_pitch; ++ii)
			ge_dest_row[ii] = 0;

		// Set identity matrix bit for this row
		ge_dest_row[ge_row_i >> 6] ^= (u64)1 << (ge_row_i & 63);

		// Initialize PRNG
		CatsChoice prng;
		prng.Initialize(_g_seed, ~ge_row_i);

		// For each peeling column,
		u32 row_bits;
		PeelColumn *column = _peel_cols;
		for (u16 column_i = 0; column_i < _block_count; ++column_i, row_bits >>= 1, ++column)
		{
			// If more bits are needed from the generator,
			if ((column_i & 31) == 0)
				row_bits = prng.Next();

			// If the column is active,
			if (row_bits & 1)
			{
				// If the column is peeled,
				if (column->mark == MARK_PEEL)
				{
					CAT_IF_DEBUG(cout << " " << column_i;)

					// Add source to destination
					u16 source_row_i = column->peel_row;
					u64 *ge_source_row = _ge_compress_matrix + _ge_compress_pitch * source_row_i;
					for (int ii = 0; ii < _ge_compress_pitch; ++ii)
						ge_dest_row[ii] ^= ge_source_row[ii];
				}
				else
				{
					CAT_IF_DEBUG(cout << " d" << column_i;)

					// Flip deferred bit
					u16 ge_column_i = column->ge_column;
					ge_dest_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
				}
			} // end if column is active
		} // next peeling column

		CAT_IF_DEBUG(cout << endl;)
	} // next dense row

	CAT_IF_DEBUG(cout << "After multiplying dense rows:" << endl;)
	CAT_IF_DEBUG(PrintGEMatrix();)
}

bool Encoder::InvCompressAllocate()
{
	CAT_IF_DEBUG(cout << endl << "---- InvCompressAllocate ----" << endl << endl;)

	// Allocate GE matrix
	int ge_rows = _defer_count + _added_count;
	int ge_pitch = (ge_rows + 63) / 64;
	int ge_matrix_words = ge_rows * ge_pitch;
	_ge_matrix = new u64[ge_matrix_words];
	if (!_ge_matrix) return false;
	_ge_pitch = ge_pitch;

	CAT_IF_DEBUG(cout << "GE matrix is " << ge_rows << " x " << ge_rows << " with pitch " << ge_pitch << " consuming " << ge_matrix_words * sizeof(u64) << " bytes" << endl;)

	// Allocate GE compress matrix
	int ge_compress_columns = ge_rows;
	int ge_compress_pitch = (ge_compress_columns + 63) / 64;
	int ge_compress_rows = _block_count;
	int ge_compress_matrix_words = ge_compress_rows * ge_compress_pitch;
	_ge_compress_matrix = new u64[ge_compress_matrix_words];
	if (!_ge_compress_matrix) return false;
	_ge_compress_pitch = ge_compress_pitch;

	// Clear entire GE compress matrix
	memset(_ge_compress_matrix, 0, ge_compress_matrix_words * sizeof(u64));

	CAT_IF_DEBUG(cout << "GE compress matrix is " << ge_compress_rows << " x " << ge_compress_columns << " with pitch " << ge_compress_pitch << " consuming " << ge_compress_matrix_words * sizeof(u64) << " bytes" << endl;)

	// Allocate the pivots
	int pivot_words = ge_rows * 2 + _defer_count;
	_ge_pivots = new u16[pivot_words];
	if (!_ge_pivots) return false;
	_ge_col_map = _ge_pivots + ge_rows;
	_ge_row_map = _ge_col_map + ge_rows;

	CAT_IF_DEBUG(cout << "Allocated " << ge_rows << " pivots, consuming " << pivot_words*2 << " bytes" << endl;)

	return true;
}

void Encoder::InvCompress()
{
	CAT_IF_DEBUG(cout << endl << "---- InvCompress ----" << endl << endl;)

	/*
		(1) Peeling matrix inversion -> [ Defer | Mix ] rows + row values in peeling columns, check if rows are deferred
		(2) Copy lower deferred rows to GE matrix
		(3) Perform peeling matrix multiplication and generate GE rows, but no row values yet
	*/

	InvPeelDiagonal();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	InvCopyDeferredRows();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	InvMultiplyDenseRows();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )
}

void Encoder::InvSolveTriangleColumns()
{
	CAT_IF_DEBUG(cout << endl << "---- InvSolveTriangleColumns ----" << endl << endl;)

	// For each pivot,
	const u16 ge_rows = _defer_count + _added_count;
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
	{
		u16 pivot_column_i = _ge_col_map[pivot_i];
		u16 ge_row_i = _ge_pivots[pivot_i];
		u8 *buffer_dest = _check_blocks + _block_bytes * pivot_column_i;

		CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " solving column " << pivot_column_i << " with GE row " << ge_row_i << endl;)

		/*
			For each GE row,

			(1) Row value (if deferred), else 0 for check rows
			(2) Regenerate row and add peeled columns
			(3) For each GE matrix bit up to the diagonal, add deferred and mixed
		*/

		// If original row,
		if (ge_row_i >= _added_count)
		{
			// (1) Row value (if deferred), else 0 for check rows:

			u16 pivot_row_i = _ge_row_map[ge_row_i - _added_count];
			const u8 *buffer_src = _message_blocks + _block_bytes * pivot_row_i;

			CAT_IF_DEBUG(cout << " " << pivot_row_i << ":[" << (int)buffer_src[0] << "]";)

			// If copying from final block,
			if (pivot_row_i == _block_count - 1)
			{
				memcpy(buffer_dest, buffer_src, _final_bytes);
				memset(buffer_dest + _final_bytes, 0, _block_bytes - _final_bytes);
			}
			else
			{
				memcpy(buffer_dest, buffer_src, _block_bytes);
			}

			// (2) Regenerate row and add peeled columns:

			// Generate mixing columns for this row
			PeelRow *row = &_peel_rows[pivot_row_i];
			u16 weight = row->peel_weight;
			u16 a = row->peel_a;
			u16 column_i = row->peel_x0;
			for (;;)
			{
				// If column is peeled,
				PeelColumn *column = &_peel_cols[column_i];
				if (column->mark == MARK_PEEL)
				{
					const u8 *peel_src = _check_blocks + _block_bytes * column_i;
					memxor(buffer_dest, peel_src, _block_bytes);

					CAT_IF_DEBUG(cout << " " << column_i << "[" << (int)peel_src[0] << "]";)
				}

				if (--weight <= 0) break;

				IterateNextColumn(column_i, _block_count, _block_next_prime, a);
			}
		}
		else
		{
			// (1) Row value (if deferred), else 0 for check rows:

			CAT_IF_DEBUG(cout << " [0]";)

			memset(buffer_dest, 0, _block_bytes);

			// (2) Regenerate row and add peeled columns:

			// Initialize PRNG
			CatsChoice prng;
			prng.Initialize(_g_seed, ~ge_row_i);

			// For each peeling column,
			u32 row_bits;
			PeelColumn *column = _peel_cols;
			for (u16 column_i = 0; column_i < _block_count; ++column_i, row_bits >>= 1, ++column)
			{
				// If more bits are needed from the generator,
				if ((column_i & 31) == 0)
					row_bits = prng.Next();

				// If the column is active,
				if (row_bits & 1)
				{
					// If the column is peeled,
					if (column->mark == MARK_PEEL)
					{
						const u8 *peel_src = _check_blocks + _block_bytes * column_i;
						memxor(buffer_dest, peel_src, _block_bytes);

						CAT_IF_DEBUG(cout << " " << column_i << "[" << (int)peel_src[0] << "]";)
					}
				} // end if column is active
			} // next peeling column
		}

		// (3) For each GE matrix bit up to the diagonal, add deferred and mixed:

		// For each GE matrix bit in the row,
		u64 *ge_row = _ge_matrix + _ge_pitch * ge_row_i;
		u64 ge_mask = 1;
		for (u16 ge_column_i = 0; ge_column_i < pivot_i; ++ge_column_i)
		{
			// If bit is set,
			if (ge_row[ge_column_i >> 6] & ge_mask)
			{
				u16 column_i = _ge_col_map[ge_column_i];
				const u8 *peel_src = _check_blocks + _block_bytes * column_i;
				memxor(buffer_dest, peel_src, _block_bytes);

				CAT_IF_DEBUG(cout << " " << column_i << "=[" << (int)peel_src[0] << "]";)
			}

			ge_mask = CAT_ROL64(ge_mask, 1);
		}

		CAT_IF_DEBUG(cout << endl;)
	}
}


/*
	Multiplication-Based Matrix Compression:

		The first step of compression is to generate the
	lower matrix in a bitfield in memory.

		+-----------+---------+
		| 0 1 0 1 1 | 1 1 0 1 |
		| 0 1 0 1 0 | 0 1 1 0 |
		| 1 1 0 1 1 | 1 1 1 0 |
		| 1 0 1 1 0 | 1 0 0 1 |
		+-----------+---------+

		This is a conceptual splitting of the lower matrix
	into two parts: A peeled matrix on the left that will
	be zero at the end of compression, and a GE matrix on
	the right that will be used during Triangularization.

	In practice the algorithm is:

	For N = 9 blocks, GE compression matrix:

		+---+---+---+---+---+---+---+---+---+
		| 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
		+-----------------------------------+
		|           Deferred Rows           |
		+-----------------------------------+
		|            Dense Rows             |
		+-----------------------------------+

	The GE compression matrix includes all peeling columns,
	but does not include any of the mixing columns.  It has
	all deferred rows and the dense rows that sum to 0.

	For Additional Count = 4 blocks and Deferred Count = 2,
	GE matrix:

		+---+---+---+---+---+---+
		| 0 | 1 | 2 | 3 | 4 | 5 |
		+---------------+-------+
		| Deferred Rows | DefCol|
		+---------------+-------+
		| Dense Rows    | DenCol|
		+---------------+-------+

	The deferred rows are split into two parts.  The first set
	of columns are initialized to mixing columns for those rows
	during Compression Setup.  Those are from the mix matrix for
	deferred rows, and for dense rows they are the identity matrix.

	The final columns are copied from the GE compression matrix
	after compression completes, so that row operations during
	GE are faster, which is important since GE is O(n^3) bitops.

	So the overall compress process is:

		(1) Allocate matrices
		(2) Generate initial state of compression matrix
			(a) Deferred rows : Set to zero and regenerate peeling rows
			(b) Dense rows : Initialize straight from PRNG
		(3) Generate initial state of GE matrix
			(a) Deferred rows : Set to zero and regenerate mixing columns for first A columns
			(b) Dense rows : Set to identity matrix
			NOTE: Leaves last columns, for deferred columns, uninitialized
		(4) Compress
			(a) Compress using compress matrix, updating mixing rows in GE matrix also
			(b) Copy deferred columns over to the GE matrix

	I split these up into different functions to make it easier to follow.

	Further description:

	The first part stores all peeled rows that were used to
	construct each GE matrix row, so that will include
	enough columns for the original message columns but not
	the additional check blocks.

	It will be generated by regenerating the peeling columns
	for each of the deferred rows, and setting a 1 bit in those
	positions if they are not deferred columns.  If they are
	deferred columns, then they are instead set in the GE matrix.

	Then from the last to first peeled column, the peeled rows
	will be XOR'd in where 1 bits are set in rows of the compression
	matrix for that column, taking care to preserve the leading
	column bit.  Similarly, the deferred columns are XOR'd
	into the GE matrix instead of the GE compression matrix.

	The result is that the GE matrix is generated, and the GE
	compression matrix consists of only the peeled rows that need to
	be summed to produce each row in the GE matrix.  GE is then
	performed to triangularize the GE matrix, which is guaranteed to
	succeed for the encoder (by precomputation).  Once the mapping
	between output columns and input rows is understood, the rows
	are calculated, then the rows are summed using the operation
	order determined during triangularization, and finally the
	diagonalization step performs back-substitution to produce the
	final block values for each of the GE columns.
*/

bool Encoder::MulCompressSetup()
{
	CAT_IF_DEBUG(cout << endl << "---- CompressSetup ----" << endl << endl;)

	if (!MulCompressAllocate())
		return false;

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	MulFillCompressDeferred();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	MulFillCompressDense();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	MulFillGEDeferred();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	MulFillGEDense();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	return true;
}

bool Encoder::MulCompressAllocate()
{
	CAT_IF_DEBUG(cout << endl << "---- MulCompressAllocate ----" << endl << endl;)

	// Allocate GE matrix
	int ge_rows = _defer_count + _added_count;
	int ge_pitch = (ge_rows + 63) / 64;
	int ge_matrix_words = ge_rows * ge_pitch;
	_ge_matrix = new u64[ge_matrix_words];
	if (!_ge_matrix) return false;
	_ge_pitch = ge_pitch;

	// Initialize the entire GE matrix to zero
	memset(_ge_matrix, 0, ge_rows * _ge_pitch * sizeof(u64));

	CAT_IF_DEBUG(cout << "GE matrix is " << ge_rows << " x " << ge_rows << " with pitch " << ge_pitch << " consuming " << ge_matrix_words * sizeof(u64) << " bytes" << endl;)

	// Allocate GE compress matrix
	int ge_compress_columns = _block_count;
	int ge_compress_pitch = (ge_compress_columns + 63) / 64;
	int ge_compress_matrix_words = ge_rows * ge_compress_pitch;
	_ge_compress_matrix = new u64[ge_compress_matrix_words];
	if (!_ge_compress_matrix) return false;
	_ge_compress_pitch = ge_compress_pitch;

	CAT_IF_DEBUG(cout << "GE compress matrix is " << ge_rows << " x " << ge_compress_columns << " with pitch " << ge_compress_pitch << " consuming " << ge_compress_matrix_words * sizeof(u64) << " bytes" << endl;)

	// Allocate the pivots
	int pivot_words = ge_rows * 2 + _defer_count;
	_ge_pivots = new u16[pivot_words];
	if (!_ge_pivots) return false;
	_ge_col_map = _ge_pivots + ge_rows;
	_ge_row_map = _ge_col_map + ge_rows;

	CAT_IF_DEBUG(cout << "Allocated " << ge_rows << " pivots, consuming " << pivot_words*2 << " bytes" << endl;)

	return true;
}

void Encoder::MulFillCompressDeferred()
{
	CAT_IF_DEBUG(cout << endl << "---- MulFillCompressDeferred ----" << endl << endl;)

	// Initialize the deferred rows to zero
	u64 *ge_compress_row = _ge_compress_matrix + _ge_compress_pitch * _added_count;
	memset(ge_compress_row, 0, _defer_count * _ge_compress_pitch * sizeof(u64));

	// Fill the deferred rows from the peel matrix
	PeelRow *row;
	for (u16 row_i = _defer_head_rows; row_i != LIST_TERM; row_i = row->next, ge_compress_row += _ge_compress_pitch)
	{
		row = &_peel_rows[row_i];

		CAT_IF_DEBUG(cout << "  Setting deferred row " << row_i << endl;)

		// For each peeling column in the row,
		u16 a = row->peel_a;
		u16 weight = row->peel_weight;
		u16 column_i = row->peel_x0;
		for (;;)
		{
			// Set bit for that column
			ge_compress_row[column_i >> 6] |= (u64)1 << (column_i & 63);

			if (--weight <= 0) break;

			IterateNextColumn(column_i, _block_count, _block_next_prime, a);
		}
	}
}

void Encoder::MulFillCompressDense()
{
	CAT_IF_DEBUG(cout << endl << "---- MulFillCompressDense ----" << endl << endl;)

	// Initialize PRNG
	CatsChoice prng;
	prng.Initialize(_g_seed, ~_block_count);

	// Fill the dense matrix with random bits
	// TODO: Optimize this for efficient compression
	int fill_count = _ge_compress_pitch * _added_count;
	u64 *ge_compress_row = _ge_compress_matrix;
	while (fill_count--) *ge_compress_row++ = ((u64)prng.Next() << 32) | prng.Next();

	CAT_IF_DEBUG(
		if (_added_count > 64)
		{
			cout << "Cannot have added count > 64 right now" << endl;
		}
	)
}

void Encoder::MulFillGEDeferred()
{
	CAT_IF_DEBUG(cout << endl << "---- MulFillGEDeferred ----" << endl << endl;)

	// Fill the deferred rows from the mix matrix
	u64 *ge_row = _ge_matrix + _ge_pitch * _added_count;
	u16 *ge_row_map = _ge_row_map;
	PeelRow *row;
	for (u16 row_i = _defer_head_rows; row_i != LIST_TERM; row_i = row->next, ge_row += _ge_pitch)
	{
		row = &_peel_rows[row_i];

		CAT_IF_DEBUG(cout << "  Setting deferred row " << row_i << endl;)

		// For each mixing column in the row,
		u16 a = row->mix_a;
		u16 weight = row->mix_weight;
		u16 x = row->mix_x0;
		for (;;)
		{
			// Set bit for that column
			u16 column_i = _defer_count + x;
			ge_row[column_i >> 6] |= (u64)1 << (column_i & 63);

			if (--weight <= 0) break;

			IterateNextColumn(x, _added_count, _added_next_prime, a);
		}

		// Remember which row was deferred into this GE matrix row
		*ge_row_map++ = row_i;
	}
}

void Encoder::MulFillGEDense()
{
	CAT_IF_DEBUG(cout << endl << "---- MulFillGEDense ----" << endl << endl;)

	// Set the dense mixing bits to the identity matrix:
	u64 *ge_row = _ge_matrix;
	for (u16 x = 0; x < _added_count; ++x, ge_row += _ge_pitch)
	{
		u16 column_i = _defer_count + x;
		ge_row[column_i >> 6] |= (u64)1 << (column_i & 63);

		CAT_IF_DEBUG(cout << "  Setting dense row " << column_i << endl;)

		// Initialize column indices
		_ge_col_map[column_i] = _block_count + x;
	}
}

void Encoder::MulCopyDeferredColumns()
{
	CAT_IF_DEBUG(cout << endl << "---- MulCopyDeferredColumns ----" << endl << endl;)

	CAT_IF_DEBUG(cout << "Copy deferred columns list to pivots array:" << endl;)

	// For each deferred column,
	for (u16 ge_column_i = 0, column_i = _defer_head_columns; column_i != LIST_TERM; column_i = _peel_cols[column_i].next, ++ge_column_i)
	{
		CAT_IF_DEBUG(cout << "  Copying deferred column " << column_i << "." << endl;)

		// Copy the deferred columns from the compress matrix to the GE matrix
		u64 *ge_compress_row = _ge_compress_matrix + (column_i >> 6);
		u64 *ge_row = _ge_matrix + (ge_column_i >> 6);
		u64 ge_compress_mask = (u64)1 << (column_i & 63);
		u64 ge_mask = (u64)1 << (ge_column_i & 63);
		for (u16 ii = 0; ii < _defer_count + _added_count; ++ii)
		{
			// If it is set in compress matrix,
			if (*ge_compress_row & ge_compress_mask)
			{
				// Set bit in GE matrix
				*ge_row |= ge_mask;

				// Clear original bit from compress matrix
				*ge_compress_row ^= ge_compress_mask;
			}

			ge_compress_row += _ge_compress_pitch;
			ge_row += _ge_pitch;
		}

		// Set column for this pivot
		_ge_col_map[ge_column_i] = column_i;
	}
}

/*
		The second step of compression is to zero the left
	matrix and generate temporary block values:

		For each column in the left matrix from right to left,
			Regenerate the peeled row that solves that column,
			generating bitfield T.

			For each row with a bit set in that column,
				Add T to that row, including the row values.
			Next
		Next

		It is clear that at the end of this process that the
	left conceptual matrix will be zero, and each row will have
	a block value associated with it.

		+-----------+---------+---+
		| 0 0 0 0 0 | 0 1 1 1 | a |
		| 0 0 0 0 0 | 0 1 1 0 | b |
		| 0 0 0 0 0 | 1 0 0 1 | c |
		| 0 0 0 0 0 | 0 1 0 1 | d |
		+-----------+---------+---+

		In practice, the left matrix has a 1 bit set for each
	peeled column that was added.  And the block values have not
	been calculated yet.  The calculation of the block values are
	deferred so that the output blocks can be used to store them
	without additional temporary storage space.
*/

void Encoder::MulCompress()
{
	CAT_IF_DEBUG(cout << endl << "---- MulCompress ----" << endl << endl;)

	// Relink the peel list in reverse order
	u16 peel_head_rows = LIST_TERM;

	// For each column that has been peeled,
	u16 row_i = _peel_head_rows;
	u16 ge_rows = _defer_count + _added_count;
	while (row_i != LIST_TERM)
	{
		PeelRow *row = &_peel_rows[row_i];
		u16 column_i = row->peel_column;

		CAT_IF_DEBUG(cout << "For peeled row " << row_i << " solved by column " << column_i << ":";)

		// Generate a list of rows that have this bit set
		u16 count = 0;
		u64 *ge_row = _ge_compress_matrix + (column_i >> 6);
		u64 ge_mask = (u64)1 << (column_i & 63);
		for (u16 ge_row_i = 0; ge_row_i < ge_rows; ++ge_row_i)
		{
			// If bit is set in this row,
			if (*ge_row & ge_mask)
			{
				// Clear it
				*ge_row ^= ge_mask;

				CAT_IF_DEBUG(cout << " " << ge_row_i;)

				// And remember this row had it set
				_ge_pivots[count++] = ge_row_i;
			}

			// Iterate next row
			ge_row += _ge_compress_pitch;
		}
		CAT_IF_DEBUG(cout << endl;)

		// If any of the rows contain it,
		if (count > 0)
		{
			// For each peeling column in the row,
			u16 a = row->peel_a;
			u16 weight = row->peel_weight;
			u16 column_i = row->peel_x0;
			for (;;)
			{
				// Set bit for that column in each row
				u64 *ge_row = _ge_compress_matrix + (column_i >> 6);
				u64 ge_mask = (u64)1 << (column_i & 63);
				for (u16 ii = 0; ii < count; ++ii)
					ge_row[_ge_compress_pitch * _ge_pivots[ii]] ^= ge_mask;

				if (--weight <= 0) break;

				IterateNextColumn(column_i, _block_count, _block_next_prime, a);
			}

			// For each mixing column in the row,
			a = row->mix_a;
			weight = row->mix_weight;
			u16 x = row->mix_x0;
			for (;;)
			{
				// Set bit for that column in each row
				u16 column_i = _defer_count + x;
				u64 *ge_row = _ge_matrix + (column_i >> 6);
				u64 ge_mask = (u64)1 << (column_i & 63);
				for (u16 ii = 0; ii < count; ++ii)
					ge_row[_ge_pitch * _ge_pivots[ii]] ^= ge_mask;

				if (--weight <= 0) break;

				IterateNextColumn(column_i, _added_count, _added_next_prime, a);
			}
		} // end if any rows contain it

		// Relink the peel list in reverse order
		u16 prev_row_i = row_i;

		// Iterate row index
		row_i = row->next;

		// Relink the peel list in reverse order
		row->next = peel_head_rows;
		peel_head_rows = prev_row_i;
	} // end for each peeled column

	// Relink the peel list in reverse order
	_peel_head_rows = peel_head_rows;

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	MulCopyDeferredColumns();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )
}

CAT_INLINE void Encoder::GenerateWindowTable16(const u8 **window_table, u16 active, u16 peel_column_i)
{
	// Generate single bit window values
	u16 peel_row_i = _peel_cols[peel_column_i].peel_row;
	if (peel_row_i != _block_count - 1)
		window_table[1] = _message_blocks + _block_bytes * peel_row_i;
	else
	{
		window_table[1] = window_table[0];
		memcpy((u8*)window_table[0], _message_blocks + _block_bytes * peel_row_i, _final_bytes);
		memset((u8*)window_table[0] + _final_bytes, 0, _block_bytes - _final_bytes);
	}

	peel_row_i = _peel_cols[peel_column_i + 1].peel_row;
	if (peel_row_i != _block_count - 1)
		window_table[2] = _message_blocks + _block_bytes * peel_row_i;
	else
	{
		window_table[2] = window_table[0];
		memcpy((u8*)window_table[0], _message_blocks + _block_bytes * peel_row_i, _final_bytes);
		memset((u8*)window_table[0] + _final_bytes, 0, _block_bytes - _final_bytes);
	}

	peel_row_i = _peel_cols[peel_column_i + 2].peel_row;
	if (peel_row_i != _block_count - 1)
		window_table[4] = _message_blocks + _block_bytes * peel_row_i;
	else
	{
		window_table[4] = window_table[0];
		memcpy((u8*)window_table[0], _message_blocks + _block_bytes * peel_row_i, _final_bytes);
		memset((u8*)window_table[0] + _final_bytes, 0, _block_bytes - _final_bytes);
	}

	peel_row_i = _peel_cols[peel_column_i + 3].peel_row;
	if (peel_row_i != _block_count - 1)
		window_table[8] = _message_blocks + _block_bytes * peel_row_i;
	else
	{
		window_table[8] = window_table[0];
		memcpy((u8*)window_table[0], _message_blocks + _block_bytes * peel_row_i, _final_bytes);
		memset((u8*)window_table[0] + _final_bytes, 0, _block_bytes - _final_bytes);
	}

	// Generate table entries for combinations of two bits
	if (active & (1 << (1 + 2))) memxor((u8*)window_table[1+2], window_table[1], window_table[2], _block_bytes);
	if (active & (1 << (1 + 4))) memxor((u8*)window_table[1+4], window_table[1], window_table[4], _block_bytes);
	if (active & (1 << (1 + 8))) memxor((u8*)window_table[1+8], window_table[1], window_table[8], _block_bytes);
	if (active & (1 << (2 + 4))) memxor((u8*)window_table[2+4], window_table[2], window_table[4], _block_bytes);
	if (active & (1 << (2 + 8))) memxor((u8*)window_table[2+8], window_table[2], window_table[8], _block_bytes);
	if (active & (1 << (4 + 8))) memxor((u8*)window_table[4+8], window_table[4], window_table[8], _block_bytes);

	// Generate table entries for combinations of three bits
	if (active & (1 << (1 + 2 + 4)))
	{
		if (active & (1 << (1 + 2)))
			memxor((u8*)window_table[1+2+4], window_table[1+2], window_table[4], _block_bytes);
		else if (active & (1 << (1 + 4)))
			memxor((u8*)window_table[1+2+4], window_table[1+4], window_table[2], _block_bytes);
		else if (active & (1 << (2 + 4)))
			memxor((u8*)window_table[1+2+4], window_table[2+4], window_table[1], _block_bytes);
		else
		{
			memxor((u8*)window_table[1+2+4], window_table[1], window_table[2], _block_bytes);
			memxor((u8*)window_table[1+2+4], window_table[4], _block_bytes);
		}
	}
	if (active & (1 << (1 + 2 + 8)))
	{
		if (active & (1 << (1 + 2)))
			memxor((u8*)window_table[1+2+8], window_table[1+2], window_table[8], _block_bytes);
		else if (active & (1 << (1 + 8)))
			memxor((u8*)window_table[1+2+8], window_table[1+8], window_table[2], _block_bytes);
		else if (active & (1 << (2 + 8)))
			memxor((u8*)window_table[1+2+8], window_table[2+8], window_table[1], _block_bytes);
		else
		{
			memxor((u8*)window_table[1+2+8], window_table[1], window_table[2], _block_bytes);
			memxor((u8*)window_table[1+2+8], window_table[8], _block_bytes);
		}
	}
	if (active & (1 << (1 + 4 + 8)))
	{
		if (active & (1 << (1 + 4)))
			memxor((u8*)window_table[1+4+8], window_table[1+4], window_table[8], _block_bytes);
		else if (active & (1 << (1 + 8)))
			memxor((u8*)window_table[1+4+8], window_table[1+8], window_table[4], _block_bytes);
		else if (active & (1 << (4 + 8)))
			memxor((u8*)window_table[1+4+8], window_table[4+8], window_table[1], _block_bytes);
		else
		{
			memxor((u8*)window_table[1+4+8], window_table[1], window_table[4], _block_bytes);
			memxor((u8*)window_table[1+4+8], window_table[8], _block_bytes);
		}
	}
	if (active & (1 << (2 + 4 + 8)))
	{
		if (active & (1 << (2 + 4)))
			memxor((u8*)window_table[2+4+8], window_table[2+4], window_table[8], _block_bytes);
		else if (active & (1 << (2 + 8)))
			memxor((u8*)window_table[2+4+8], window_table[2+8], window_table[4], _block_bytes);
		else if (active & (1 << (4 + 8)))
			memxor((u8*)window_table[2+4+8], window_table[4+8], window_table[2], _block_bytes);
		else
		{
			memxor((u8*)window_table[2+4+8], window_table[2], window_table[4], _block_bytes);
			memxor((u8*)window_table[2+4+8], window_table[8], _block_bytes);
		}
	}

	// Generate table entries for combinations of four bits
	if (active & (1 << (1 + 2 + 4 + 8)))
	{
		if (active & (1 << (1 + 2 + 4)))
			memxor((u8*)window_table[1+2+4+8], window_table[1+2+4], window_table[8], _block_bytes);
		else if (active & (1 << (1 + 2 + 8)))
			memxor((u8*)window_table[1+2+4+8], window_table[1+2+8], window_table[4], _block_bytes);
		else if (active & (1 << (1 + 4 + 8)))
			memxor((u8*)window_table[1+2+4+8], window_table[1+4+8], window_table[2], _block_bytes);
		else if (active & (1 << (2 + 4 + 8)))
			memxor((u8*)window_table[1+2+4+8], window_table[2+4+8], window_table[1], _block_bytes);
		else
		{
			if (active & (1 << (1 + 2)))
			{
				if (active & (1 << (4 + 8)))
					memxor((u8*)window_table[1+2+4+8], window_table[1+2], window_table[4+8], _block_bytes);
				else
				{
					memxor((u8*)window_table[1+2+4+8], window_table[1+2], window_table[4], _block_bytes);
					memxor((u8*)window_table[1+2+4+8], window_table[8], _block_bytes);
				}
			}
			else if (active & (1 << (1 + 4)))
			{
				if (active & (1 << (2 + 8)))
					memxor((u8*)window_table[1+2+4+8], window_table[1+4], window_table[2+8], _block_bytes);
				else
				{
					memxor((u8*)window_table[1+2+4+8], window_table[1+4], window_table[2], _block_bytes);
					memxor((u8*)window_table[1+2+4+8], window_table[8], _block_bytes);
				}
			}
			else if (active & (1 << (1 + 8)))
			{
				if (active & (1 << (2 + 4)))
					memxor((u8*)window_table[1+2+4+8], window_table[1+8], window_table[2+4], _block_bytes);
				else
				{
					memxor((u8*)window_table[1+2+4+8], window_table[1+8], window_table[2], _block_bytes);
					memxor((u8*)window_table[1+2+4+8], window_table[4], _block_bytes);
				}
			}
			else if (active & (1 << (2 + 4)))
			{
				memxor((u8*)window_table[1+2+4+8], window_table[2+4], window_table[1], _block_bytes);
				memxor((u8*)window_table[1+2+4+8], window_table[8], _block_bytes);
			}
			else if (active & (1 << (2 + 8)))
			{
				memxor((u8*)window_table[1+2+4+8], window_table[2+8], window_table[1], _block_bytes);
				memxor((u8*)window_table[1+2+4+8], window_table[4], _block_bytes);
			}
			else if (active & (1 << (4 + 8)))
			{
				memxor((u8*)window_table[1+2+4+8], window_table[4+8], window_table[1], _block_bytes);
				memxor((u8*)window_table[1+2+4+8], window_table[2], _block_bytes);
			}
			else
			{
				memxor((u8*)window_table[1+2+4+8], window_table[1], window_table[2], _block_bytes);
				memxor((u8*)window_table[1+2+4+8], window_table[4], _block_bytes);
				memxor((u8*)window_table[1+2+4+8], window_table[8], _block_bytes);
			}
		}
	}
}

bool Encoder::MulSolveTriangleColumnsWindowed()
{
	CAT_IF_DEBUG(cout << endl << "---- MulSolveTriangleColumnsWindowed ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	const u8 *window_table[16];

	// Find 15 peeled columns that can be used as scratch space for this function
	u16 column_i = 0;
	for (int window_i = 0; window_i < 16; ++window_i)
	{
		while (_peel_cols[column_i].mark != MARK_PEEL)
		{
			if (++column_i >= _block_count)
				return false;
		}
		window_table[window_i] = _check_blocks + _block_bytes * column_i;
	}

	const u16 ge_rows = _defer_count + _added_count;
	const u8 *buffer_src = _message_blocks;

	// Unroll first loop
	{
		// Check which window values are active
		u64 *ge_compress_row = _ge_compress_matrix;
		u16 active = 0;
		for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i, ge_compress_row += _ge_compress_pitch)
		{
			u32 bits = *ge_compress_row & 15;
			active |= 1 << bits;
		}

		GenerateWindowTable16(window_table, active, 0);

		// Add in window values
		ge_compress_row = _ge_compress_matrix;
		for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i, ge_compress_row += _ge_compress_pitch)
		{
			u16 pivot_ge_row_i = _ge_pivots[pivot_i];
			u16 pivot_column_i = _ge_col_map[pivot_i];
			u8 *buffer_dest = _check_blocks + _block_bytes * pivot_column_i;

			// If bits are zero,
			u32 bits = *ge_compress_row & 15;
			if (bits == 0)
			{
				// If original row,
				if (pivot_ge_row_i >= _added_count)
				{
					u16 pivot_row_i = _ge_row_map[pivot_ge_row_i - _added_count];
					const u8 *buffer_src = _message_blocks + _block_bytes * pivot_row_i;

					// If copying from final block,
					if (pivot_row_i == _block_count - 1)
					{
						memcpy(buffer_dest, buffer_src, _final_bytes);
						memset(buffer_dest + _final_bytes, 0, _block_bytes - _final_bytes);
					}
					else
					{
						memcpy(buffer_dest, buffer_src, _block_bytes);
					}
				}
				else
				{
					memset(buffer_dest, 0, _block_bytes);
				}
				CAT_IF_ROWOP(++rowops;)

				continue;
			}
			const u8 *buffer_src = window_table[bits];

			// If original row,
			if (pivot_ge_row_i >= _added_count)
			{
				u16 pivot_row_i = _ge_row_map[pivot_ge_row_i - _added_count];
				const u8 *message_src = _message_blocks + _block_bytes * pivot_row_i;

				// If copying from final block,
				if (pivot_row_i == _block_count - 1)
				{
					memcpy(buffer_dest, buffer_src, _block_bytes);
					memxor(buffer_dest, message_src, _final_bytes);
				}
				else
				{
					memxor(buffer_dest, buffer_src, message_src, _block_bytes);
				}
			}
			else
			{
				memcpy(buffer_dest, buffer_src, _block_bytes);
			}
			CAT_IF_ROWOP(++rowops;)
		}
	}

	// For each bit of the compress row,
	u64 ge_mask = 15 << 4;
	int ge_shift = 4;
	u16 peel_column_i;
	for (peel_column_i = 4; peel_column_i < _block_count - 4; peel_column_i += 4)
	{
		// Check which window values are active
		u64 *ge_compress_row = _ge_compress_matrix + (peel_column_i >> 6);
		u16 active = 0;
		for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i, ge_compress_row += _ge_compress_pitch)
		{
			u32 bits = (u32)((*ge_compress_row & ge_mask) >> ge_shift);
			active |= 1 << bits;
		}

		GenerateWindowTable16(window_table, active, peel_column_i);

		// Add in window values
		ge_compress_row = _ge_compress_matrix + (peel_column_i >> 6);
		for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i, ge_compress_row += _ge_compress_pitch)
		{
			u32 bits = (u32)((*ge_compress_row & ge_mask) >> ge_shift);
			if (bits == 0) continue;
			const u8 *buffer_src = window_table[bits];

			u16 pivot_column_i = _ge_col_map[pivot_i];
			u16 pivot_ge_row_i = _ge_pivots[pivot_i];
			u8 *buffer_dest = _check_blocks + _block_bytes * pivot_column_i;

			memxor(buffer_dest, buffer_src, _block_bytes);
		}

		// Generate next mask and shift
		ge_mask = CAT_ROL64(ge_mask, 4);
		ge_shift = (ge_shift + 4) & 63;
	}

	// For each remaining row,
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
	{
		u16 pivot_column_i = _ge_col_map[pivot_i];
		u16 pivot_ge_row_i = _ge_pivots[pivot_i];
		u8 *buffer_dest = _check_blocks + _block_bytes * pivot_column_i;

		CAT_IF_DEBUG(cout << "For compress row " << pivot_ge_row_i << " that solves column " << pivot_column_i << ":";)

		// For each bit of the compress row,
		u64 *ge_compress_row = _ge_compress_matrix + _ge_compress_pitch * pivot_ge_row_i;
		u64 ge_mask = (u64)1 << (peel_column_i & 63);
		for (u16 peel_column_j = peel_column_i; peel_column_i < _block_count; ++peel_column_i)
		{
			// If bit is set,
			if (ge_compress_row[peel_column_i >> 6] & ge_mask)
			{
				CAT_IF_DEBUG(cout << " " << peel_column_i;)

				// Look up which peeled row solves this column
				PeelColumn *column = &_peel_cols[peel_column_i];
				u16 peel_row_i = column->peel_row;

				// Look up source data
				const u8 *buffer_src = _message_blocks + _block_bytes * peel_row_i;
				int block_bytes = (peel_row_i == _block_count - 1) ? _final_bytes : _block_bytes;

				memxor(buffer_dest, buffer_src, block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)buffer_src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}

			// Generate next mask
			ge_mask = CAT_ROL64(ge_mask, 1);
		}

		CAT_IF_DEBUG(cout << endl << "  For GE row " << pivot_ge_row_i << " solved by pivot " << pivot_i << ":";)

		// For each bit of the GE matrix up to the diagonal,
		u64 *ge_row = _ge_matrix + _ge_pitch * pivot_ge_row_i;
		ge_mask = 1;
		for (u16 ge_column_i = 0; ge_column_i < pivot_i; ++ge_column_i)
		{
			// If bit is set,
			if (ge_row[ge_column_i >> 6] & ge_mask)
			{
				CAT_IF_DEBUG(cout << " " << ge_column_i;)

				// Look up source data
				u16 ge_pivot_column_i = _ge_col_map[ge_column_i];
				const u8 *buffer_src = _check_blocks + _block_bytes * ge_pivot_column_i;

				memxor(buffer_dest, buffer_src, _block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)buffer_src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}

			// Generate next mask
			ge_mask = CAT_ROL64(ge_mask, 1);
		}

		CAT_IF_DEBUG(cout << endl;)
		CAT_IF_DEBUG(cout << "Result = " << (int)buffer_dest[0] << endl;)
	}

	CAT_IF_ROWOP(cout << "SolveTriangleColumns used " << rowops << " row ops" << endl;)

	return true;
}

void Encoder::MulSolveTriangleColumns()
{
	u16 ge_rows = _defer_count + _added_count;

	// Attempt to solve using window optimization
	if (ge_rows >= CAT_WINDOW_THRESHOLD && MulSolveTriangleColumnsWindowed())
		return;

	CAT_IF_DEBUG(cout << endl << "---- MulSolveTriangleColumns ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// For each pivot column,
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
	{
		u16 pivot_column_i = _ge_col_map[pivot_i];
		u16 pivot_ge_row_i = _ge_pivots[pivot_i];
		u8 *buffer_dest = _check_blocks + _block_bytes * pivot_column_i;

		CAT_IF_DEBUG(cout << "For compress row " << pivot_ge_row_i << " that solves column " << pivot_column_i << ":";)

		// If original row,
		if (pivot_ge_row_i >= _added_count)
		{
			u16 pivot_row_i = _ge_row_map[pivot_ge_row_i - _added_count];
			const u8 *buffer_src = _message_blocks + _block_bytes * pivot_row_i;

			CAT_IF_DEBUG(cout << " " << pivot_row_i << ":[" << (int)buffer_src[0] << "]";)

			// If copying from final block,
			if (pivot_row_i == _block_count - 1)
			{
				memcpy(buffer_dest, buffer_src, _final_bytes);
				memset(buffer_dest + _final_bytes, 0, _block_bytes - _final_bytes);
			}
			else
			{
				memcpy(buffer_dest, buffer_src, _block_bytes);
			}
		}
		else
		{
			CAT_IF_DEBUG(cout << " [0]";)

			memset(buffer_dest, 0, _block_bytes);
		}
		CAT_IF_ROWOP(++rowops;)

		// For each bit of the compress row,
		u64 *ge_compress_row = _ge_compress_matrix + _ge_compress_pitch * pivot_ge_row_i;
		u64 ge_mask = 1;
		for (u16 peel_column_i = 0; peel_column_i < _block_count; ++peel_column_i)
		{
			// If bit is set,
			if (ge_compress_row[peel_column_i >> 6] & ge_mask)
			{
				CAT_IF_DEBUG(cout << " " << peel_column_i;)

				// Look up which peeled row solves this column
				PeelColumn *column = &_peel_cols[peel_column_i];
				u16 peel_row_i = column->peel_row;

				// Look up source data
				const u8 *buffer_src = _message_blocks + _block_bytes * peel_row_i;
				int block_bytes = (peel_row_i == _block_count - 1) ? _final_bytes : _block_bytes;

				memxor(buffer_dest, buffer_src, block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)buffer_src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}

			// Generate next mask
			ge_mask = CAT_ROL64(ge_mask, 1);
		}

		CAT_IF_DEBUG(cout << endl << "  For GE row " << pivot_ge_row_i << " solved by pivot " << pivot_i << ":";)

		// For each bit of the GE matrix up to the diagonal,
		u64 *ge_row = _ge_matrix + _ge_pitch * pivot_ge_row_i;
		ge_mask = 1;
		for (u16 ge_column_i = 0; ge_column_i < pivot_i; ++ge_column_i)
		{
			// If bit is set,
			if (ge_row[ge_column_i >> 6] & ge_mask)
			{
				CAT_IF_DEBUG(cout << " " << ge_column_i;)

				// Look up source data
				u16 ge_pivot_column_i = _ge_col_map[ge_column_i];
				const u8 *buffer_src = _check_blocks + _block_bytes * ge_pivot_column_i;

				memxor(buffer_dest, buffer_src, _block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)buffer_src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}

			// Generate next mask
			ge_mask = CAT_ROL64(ge_mask, 1);
		}

		CAT_IF_DEBUG(cout << endl;)
		CAT_IF_DEBUG(cout << "Result = " << (int)buffer_dest[0] << endl;)
	}

	CAT_IF_ROWOP(cout << "SolveTriangleColumns used " << rowops << " row ops" << endl;)
}


//// Diagnostic functions for compression algorithms

#if defined(CAT_DEBUG)

void Encoder::PrintGEMatrix()
{
	int rows = _defer_count + _added_count;
	int cols = rows;

	cout << endl << "GE matrix is " << rows << " x " << cols << ":" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			if (_ge_matrix[_ge_pitch * ii + (jj >> 6)] & ((u64)1 << (jj & 63)))
				cout << '1';
			else
				cout << '0';
		}
		cout << endl;
	}

	cout << endl;
}

void Encoder::MulPrintGECompressMatrix()
{
	int rows = _defer_count + _added_count;
	int cols = _block_count;

	cout << endl << "GE (Mul) Compress matrix is " << rows << " x " << cols << ":" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			if (_ge_compress_matrix[_ge_compress_pitch * ii + (jj >> 6)] & ((u64)1 << (jj & 63)))
				cout << '1';
			else
				cout << '0';
		}
		cout << endl;
	}

	cout << endl;
}

void Encoder::InvPrintGECompressMatrix()
{
	int rows = _block_count;
	int cols = _defer_count + _added_count;

	cout << endl << "GE (Inv) Compress matrix is " << rows << " x " << cols << ":" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			if (_ge_compress_matrix[_ge_compress_pitch * ii + (jj >> 6)] & ((u64)1 << (jj & 63)))
				cout << '1';
			else
				cout << '0';
		}
		cout << endl;
	}

	cout << endl;
}

#endif // CAT_DEBUG


/*
	One more subtle optimization.  Depending on how the GE matrix
	is constructed, it can be put in a roughly upper-triangular form
	from the start so it looks like this:

		+---------+---------+
		| D D D D | D D D 1 |
		| D D D D | D D 1 M |
		| D D D D | D 1 M M |
		| D D D D | 1 M M M |
		+---------+---------+
		| 0 1 0 1 | M M M M |
		| 0 0 1 0 | M M M M |
		| 0 0 0 1 | M M M M |
		| 0 0 0 0 | M M M M |
		+---------+---------+

	In the example above, the top 4 rows are dense matrix rows.
	The last 4 columns are mixing columns and are also dense.
	The lower left sub-matrix is sparse and roughly upper triangular
	because it is the intersection of sparse rows in the generator
	matrix.  This form is achieved by adding deferred rows starting
	from last deferred to first, and adding deferred columns starting
	from last deferred to first.

	Gaussian elimination will proceed from left to right on the matrix.
	So, having an upper triangular form will prevent left-most zeroes
	from being eaten up when a row is eliminated by one above it.  This
	reduces row operations.  For example, GE on a 64x64 matrix will do
	on average 100 fewer row operations on this form rather than its
	transpose (where the sparse part is put in the upper right).

	This is not a huge improvement but every little bit helps.
*/


//// (3) Gaussian Elimination

bool Encoder::Triangle()
{
	CAT_IF_DEBUG(cout << endl << "---- Triangle ----" << endl << endl;)

	u16 pivot_count = _defer_count + _added_count;

	// Initialize pivot array
	for (u16 pivot_i = 0; pivot_i < pivot_count; ++pivot_i)
		_ge_pivots[pivot_i] = pivot_i;

	// For each pivot to determine,
	u64 ge_mask = 1;
	for (u16 pivot_i = 0; pivot_i < pivot_count; ++pivot_i)
	{
		int word_offset = pivot_i >> 6;
		u64 *ge_matrix_offset = _ge_matrix + word_offset;

		// Find pivot
		bool found = false;
		for (u16 pivot_j = pivot_i; pivot_j < pivot_count; ++pivot_j)
		{
			// Determine if the row contains the bit we want
			u16 ge_row_j = _ge_pivots[pivot_j];
			u64 *ge_row = &ge_matrix_offset[_ge_pitch * ge_row_j];

			// If the bit was found,
			if (*ge_row & ge_mask)
			{
				found = true;
				CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " found on row " << ge_row_j << endl;)

				// Swap out the pivot index for this one
				u16 temp = _ge_pivots[pivot_i];
				_ge_pivots[pivot_i] = _ge_pivots[pivot_j];
				_ge_pivots[pivot_j] = temp;

				// Prepare masked first word
				u64 row0 = (*ge_row & ~(ge_mask - 1)) ^ ge_mask;

				// For each remaining unused row,
				for (u16 pivot_k = pivot_j + 1; pivot_k < pivot_count; ++pivot_k)
				{
					// Determine if the row contains the bit we want
					u16 ge_row_k = _ge_pivots[pivot_k];
					u64 *rem_row = &ge_matrix_offset[_ge_pitch * ge_row_k];

					// If the bit was found,
					if (*rem_row & ge_mask)
					{
						// Add the pivot row to eliminate the bit from this row, preserving previous bits
						*rem_row ^= row0;

						for (int ii = 1; ii < _ge_pitch - word_offset; ++ii)
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

void Encoder::Diagonal()
{
	CAT_IF_DEBUG(cout << endl << "---- Diagonal ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// For each pivot row from last to first,
	int pivot_i = _defer_count + _added_count - 1;
	u64 ge_mask = (u64)1 << (pivot_i & 63);
	for (; pivot_i >= 0; --pivot_i)
	{
		u16 pivot_ge_row_i = _ge_pivots[pivot_i];

		// Calculate source of copy
		u16 pivot_column_i = _ge_col_map[pivot_i];
		const u8 *src = _check_blocks + _block_bytes * pivot_column_i;

		CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " solving column " << pivot_column_i << ":";)

		// For each pivot row above it,
		u64 *ge_row = _ge_matrix + (pivot_i >> 6);
		for (int above_i = pivot_i - 1; above_i >= 0; --above_i)
		{
			u16 ge_above_row_i = _ge_pivots[above_i];

			// If bit is set in that row,
			if (ge_row[_ge_pitch * ge_above_row_i] & ge_mask)
			{
				// Back-substitute
				u16 above_column_i = _ge_col_map[above_i];
				u8 *dest = _check_blocks + _block_bytes * above_column_i;

				CAT_IF_DEBUG(cout << " " << above_column_i;)

				memxor(dest, src, _block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}
		}

		CAT_IF_DEBUG(cout << endl;)

		// Generate next mask
		ge_mask = CAT_ROR64(ge_mask, 1);
	}

	CAT_IF_ROWOP(cout << "Diagonal used " << rowops << " row ops" << endl;)
}


//// (4) Substitute

void Encoder::Substitute()
{
	CAT_IF_DEBUG(cout << endl << "---- Substitute ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// For each column that has been peeled,
	u16 ge_rows = _defer_count + _added_count;
	PeelRow *row;
	for (u16 row_i = _peel_head_rows; row_i != LIST_TERM; row_i = row->next)
	{
		row = &_peel_rows[row_i];
		u16 dest_column_i = row->peel_column;
		u8 *dest = _check_blocks + _block_bytes * dest_column_i;

		CAT_IF_DEBUG(cout << "Generating column " << dest_column_i << ":";)

		const u8 *combo, *src = _message_blocks + _block_bytes * row_i;

		CAT_IF_DEBUG(cout << " " << row_i << ":[" << (int)src[0] << "]";)

		// If copying from final block,
		if (row_i < _block_count - 1)
			combo = src;
		else
		{
			memcpy(dest, src, _final_bytes);
			memset(dest + _final_bytes, 0, _block_bytes - _final_bytes);
			combo = 0;
		}
		CAT_IF_ROWOP(++rowops;)

		// NOTE: Doing mixing columns first because mixing weight >= 1 so
		// the combo is guaranteed to be used up, and the average weight
		// of mixing columns is less than the peeling columns so the inner
		// loop of the peeling columns should be less complex.

		// For each mixing column in the row,
		u16 a = row->mix_a;
		u16 weight = row->mix_weight;
		u16 x = row->mix_x0;
		for (;;)
		{
			u16 column_i = _block_count + x;
			const u8 *src = _check_blocks + _block_bytes * column_i;

			CAT_IF_DEBUG(cout << " " << column_i;)

			if (!combo)
				memxor(dest, src, _block_bytes);
			else
			{
				memxor(dest, src, combo, _block_bytes);
				combo = 0;
			}
			CAT_IF_DEBUG(cout << "[" << (int)src[0] << "]";)
			CAT_IF_ROWOP(++rowops;)

			if (--weight <= 0) break;

			IterateNextColumn(x, _added_count, _added_next_prime, a);
		}

		// For each peeling column in the row,
		a = row->peel_a;
		weight = row->peel_weight;
		u16 column_i = row->peel_x0;
		for (;;)
		{
			const u8 *src = _check_blocks + _block_bytes * column_i;

			CAT_IF_DEBUG(cout << " " << column_i;)

			// If column is not the solved one,
			if (column_i != dest_column_i)
			{
				memxor(dest, src, _block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}
			else
			{
				CAT_IF_DEBUG(cout << "*";)
			}

			if (--weight <= 0) break;

			IterateNextColumn(column_i, _block_count, _block_next_prime, a);
		}

		CAT_IF_DEBUG(cout << endl;)
	}

	CAT_IF_ROWOP(cout << "Substitute used " << rowops << " row ops" << endl;)
}


//// Main Driver

bool Encoder::GenerateCheckBlocks()
{
	CAT_IF_DEBUG(cout << endl << "---- GenerateCheckBlocks ----" << endl << endl;)

	// (1) Peeling

	if (!PeelSetup())
		return false;

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(
	{
		cout << "After PeelSetup, contents of peeled list:" << endl;
		u16 row_i = _peel_head_rows;
		while (row_i != LIST_TERM)
		{
			PeelRow *row = &_peel_rows[row_i];

			cout << "  Row " << row_i << " solves column " << row->peel_column << endl;

			row_i = row->next;
		}
	})

	GreedyPeeling();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(
	{
		cout << "After GreedyPeeling, contents of peeled list (last to first):" << endl;
		u16 row_i = _peel_head_rows;
		while (row_i != LIST_TERM)
		{
			PeelRow *row = &_peel_rows[row_i];

			cout << "  Row " << row_i << " solves column " << row->peel_column << endl;

			row_i = row->next;
		}
	})

	CAT_IF_DEBUG(
	{
		cout << "After GreedyPeeling, contents of deferred columns list:" << endl;
		u16 column_i = _defer_head_columns;
		while (column_i != LIST_TERM)
		{
			PeelColumn *column = &_peel_cols[column_i];

			cout << "  Column " << column_i << " is deferred." << endl;

			column_i = column->next;
		}
	})

	CAT_IF_DEBUG(
	{
		cout << "After GreedyPeeling, contents of deferred rows list:" << endl;
		u16 row_i = _defer_head_rows;
		while (row_i != LIST_TERM)
		{
			PeelRow *row = &_peel_rows[row_i];

			cout << "  Row " << row_i << " is deferred." << endl;

			row_i = row->next;
		}
	})

	// (2) Compression

	if (_use_inverse_method)
	{
		if (!InvCompressSetup())
			return false;

		CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

		CAT_IF_DEBUG(cout << "After InvCompressSetup:" << endl;)
		CAT_IF_DEBUG(PrintGEMatrix();)
		CAT_IF_DEBUG(InvPrintGECompressMatrix();)

		InvCompress();
	}
	else
	{
		if (!MulCompressSetup())
			return false;

		CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

		CAT_IF_DEBUG(cout << "After CompressSetup:" << endl;)
		CAT_IF_DEBUG(PrintGEMatrix();)
		CAT_IF_DEBUG(MulPrintGECompressMatrix();)

		MulCompress();
	}

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(cout << "After Compress:" << endl;)
	CAT_IF_DEBUG(PrintGEMatrix();)

	// (3) Gaussian Elimination

	if (!Triangle())
	{
		CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

		CAT_IF_DEBUG(cout << "After Triangle FAILED:" << endl;)
		CAT_IF_DEBUG(PrintGEMatrix();)

		return false;
	}

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(cout << "After Triangle:" << endl;)
	CAT_IF_DEBUG(PrintGEMatrix();)

	if (_use_inverse_method)
	{
		InvSolveTriangleColumns();
	}
	else
	{
		MulSolveTriangleColumns();
	}

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	Diagonal();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	// (4) Substitution

	Substitute();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	return true;
}

void Encoder::Cleanup()
{
	if (_check_blocks)
	{
		delete []_check_blocks;
		_check_blocks = 0;
	}

	if (_peel_rows)
	{
		delete []_peel_rows;
		_peel_rows = 0;
	}

	if (_peel_cols)
	{
		delete []_peel_cols;
		_peel_cols = 0;
	}

	if (_ge_matrix)
	{
		delete []_ge_matrix;
		_ge_matrix = 0;
	}

	if (_ge_pivots)
	{
		delete []_ge_pivots;
		_ge_pivots = 0;
	}
}

Encoder::Encoder()
{
	_check_blocks = 0;
	_peel_rows = 0;
	_peel_cols = 0;
	_ge_matrix = 0;
	_ge_pivots = 0;
}

Encoder::~Encoder()
{
	Cleanup();
}


bool Encoder::Initialize(const void *message_in, int message_bytes, int block_bytes)
{
	CAT_IF_DEBUG(cout << endl << "---- Initialize ----" << endl << endl;)

	Cleanup();

	// Calculate message block count
	block_bytes -= 3;
	int block_count = (message_bytes + block_bytes - 1) / block_bytes;
	_final_bytes = message_bytes % block_bytes;
	if (_final_bytes <= 0) _final_bytes = block_bytes;

	// Lookup generator matrix parameters
	_added_count = GetCheckBlockCount(block_count);
	_g_seed = GetGeneratorSeed(block_count);

	// Allocate check blocks
	int check_size = (block_count + _added_count) * block_bytes;
	_check_blocks = new u8[check_size];
	if (!_check_blocks) return false;

	// Initialize encoder
	_block_bytes = block_bytes;
	_block_count = block_count;
	_message_blocks = reinterpret_cast<const u8*>( message_in );
	_next_block_id = 0;

	// Calculate next primes after column counts for pseudo-random generation of peeling rows
	_block_next_prime = NextPrime16(_block_count);
	_added_next_prime = NextPrime16(_added_count);

	CAT_IF_DEBUG(cout << "Total message = " << message_bytes << " bytes" << endl;)
	CAT_IF_DEBUG(cout << "Block bytes = " << block_bytes << ".  Final bytes = " << _final_bytes << endl;)
	CAT_IF_DEBUG(cout << "Block count = " << block_count << " +Prime=" << _block_next_prime << endl;)
	CAT_IF_DEBUG(cout << "Added count = " << _added_count << " +Prime=" << _added_next_prime << endl;)
	CAT_IF_DEBUG(cout << "Generator seed = " << _g_seed << endl;)
	CAT_IF_DEBUG(cout << "Memory overhead for check blocks = " << check_size << " bytes" << endl;)

	// Generate check blocks
	return GenerateCheckBlocks();
}

void Encoder::Generate(void *block)
{
	// Insert ID field
	u8 *buffer = reinterpret_cast<u8*>( block );
	u32 id = _next_block_id++;
	buffer[0] = (u8)id;
	buffer[1] = (u8)(id >> 8);
	buffer[2] = (u8)(id >> 16);
	buffer += 3;

#if defined(CAT_ENCODER_COPY_FIRST_N)
	// For the message blocks,
	if (id < _block_count)
	{
		// Until the final block in message blocks,
		if (id < _block_count - 1)
		{
			// Copy from the original file data
			memcpy(buffer, _message_blocks, _block_bytes);
			_message_blocks += _block_bytes;
		}
		else
		{
			// For the final block, copy partial block
			memcpy(buffer, _message_blocks, _final_bytes);

			// Pad with zeroes
			memset(buffer + _final_bytes, 0, _block_bytes - _final_bytes);
		}

		return;
	}
#endif // CAT_ENCODER_COPY_FIRST_N

	CAT_IF_DEBUG(cout << "Generating row " << id << ":";)

	// Initialize PRNG
	CatsChoice prng;
	prng.Initialize(id, _g_seed);

	// Generate peeling matrix row parameters
	u16 peel_weight = GenerateWeight(prng.Next(), _block_count - 1);
	u32 rv = prng.Next();
	u16 peel_a = ((u16)rv % (_block_count - 1)) + 1;
	u16 peel_x = (u16)(rv >> 16) % _block_count;

	// Generate mixing matrix row parameters
	u16 mix_weight = 3;
	if (mix_weight >= _added_count)
		mix_weight = _added_count - 1;
	rv = prng.Next();
	u16 mix_a = ((u16)rv % (_added_count - 1)) + 1;
	u16 mix_x = (u16)(rv >> 16) % _added_count;

	// Remember first column (there is always at least one)
	u8 *first = _check_blocks + _block_bytes * peel_x;

	CAT_IF_DEBUG(cout << " " << peel_x;)

	// If peeler has multiple columns,
	if (peel_weight > 1)
	{
		--peel_weight;

		IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

		CAT_IF_DEBUG(cout << " " << peel_x;)

		// Combine first two columns into output buffer (faster than memcpy + memxor)
		memxor(buffer, first, _check_blocks + _block_bytes * peel_x, _block_bytes);

		// For each remaining peeler column,
		while (--peel_weight > 0)
		{
			IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

			CAT_IF_DEBUG(cout << " " << peel_x;)

			// Mix in each column
			memxor(buffer, _check_blocks + _block_bytes * peel_x, _block_bytes);
		}

		// Mix first mixer block in directly
		memxor(buffer, _check_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	}
	else
	{
		// Mix first with first mixer block (faster than memcpy + memxor)
		memxor(buffer, first, _check_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	}

	CAT_IF_DEBUG(cout << " " << (_block_count + mix_x);)

	// For each remaining mixer column,
	while (--mix_weight > 0)
	{
		IterateNextColumn(mix_x, _added_count, _added_next_prime, mix_a);

		CAT_IF_DEBUG(cout << " " << (_block_count + mix_x);)

		// Mix in each column
		memxor(buffer, _check_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	}

	CAT_IF_DEBUG(cout << endl;)
}


//// Decoder

bool Decoder::Initialize(void *message_out, int message_bytes, int block_bytes)
{
	return true;
}

bool Decoder::Decode(void *block)
{
	return true;
}
