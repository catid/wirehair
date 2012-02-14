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

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;


// Switches:
//#define STEW_HYPERDYNAMIC_PLATTONIC_ITERATOR
//#define CAT_LOOPLESS_ITERATOR
//#define CAT_DUMP_ROWOP_COUNTERS


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

static int seed = 0;

static u32 GetGeneratorSeed(int block_count)
{
	// TODO: Needs to be simulated (2)
	return seed++;
}

static int GetCheckBlockCount(int block_count)
{
	// TODO: Needs to be simulated (1)
	return (int)sqrt((float)block_count) / 2 + 1;
}

#if defined(STEW_HYPERDYNAMIC_PLATTONIC_ITERATOR)
CAT_INLINE static void IterateNextColumn(u16 &x, u16 b, u16 p, u16 a)
{
	x = (x + a) % p;

	if (x >= b)
	{
		u16 distance = p - x;

		if (a < distance)
			x = (((distance + a - 1) / a) * a + x) - p;
		else
			x -= p - a;
	}
}
#elif defined(CAT_LOOPLESS_ITERATOR)
CAT_INLINE static void IterateNextColumn(u16 &x, u16 b, u16 p, u16 a)
{
	x = (x + a) % p;
	if (x >= b)	// Fix roll without a loop
		x = (((u32)a << 16) - p + x) % a;
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
	// Peeling matrix: Column generator
	u16 peel_weight, peel_a, peel_x0;

	// Mixing matrix: Column generator
	u16 mix_weight, mix_a, mix_x0;

	// Peeling state
	u16 unmarked_count;	// Count of columns that have not been marked yet
	u16 unmarked[2];	// Final two unmarked column indices

	u16 next;			// Linkage in row list
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
	u16 w2_refs;	// Number of weight-2 rows containing this column
	u16 row_count;	// Number of rows containing this column
	u16 rows[ENCODER_REF_LIST_MAX];
	u8 mark;		// One of the MarkTypes enumeration

	u16 next;		// Linkage in column list
};
#pragma pack(pop)

/*
	Peel() and PeelAvalanche() are split up into two functions because I found
	that the PeelAvalanche() function can be reused later during GreedyPeeling().

	(1) Peeling:

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
			if (_peel_cols[new_column_i].mark != MARK_TODO)
			{
				CAT_IF_DEBUG(cout << "PeelAvalanche: Deferred(1) with column " << column_i << " at row " << ref_row_i << endl;)

				// Link at head of defer list
				ref_row->next = _defer_head_rows;
				_defer_head_rows = ref_row_i;
			}
			else
			{
				Peel(ref_row_i, ref_row, new_column_i);
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
				if (unmarked_count == 0)
				{
					CAT_IF_DEBUG(cout << "PeelAvalanche: Deferred(2) with column " << column_i << " at row " << ref_row_i << endl;)

					// Link at head of defer list
					ref_row->next = _defer_head_rows;
					_defer_head_rows = ref_row_i;
				}
				else
				{
					Peel(ref_row_i, ref_row, ref_row->unmarked[0]);
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
	row->unmarked[0] = column_i;

	// Link to front of the peeled list
	row->next = _peel_head_rows;
	_peel_head_rows = row_i;

	// Attempt to avalanche and solve other columns
	PeelAvalanche(column_i, column);

	// Remember which row solves the column, after done with rows list
	column->w2_refs = row_i;
}

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
	_peel_head_rows = LIST_TERM;
	_defer_head_rows = LIST_TERM;

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
				cout << "PeelSetup: Failure!  Ran out of space for row references.  ENCODER_REF_LIST_MAX must be increased!" << endl;
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

	(2) Compression:

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

bool Encoder::CompressSetup()
{
	CAT_IF_DEBUG(cout << endl << "---- CompressSetup ----" << endl << endl;)

	if (!CompressAllocate())
		return false;

	FillCompressDeferred();
	FillCompressDense();
	FillGEDeferred();
	FillGEDense();

	return true;
}

bool Encoder::CompressAllocate()
{
	CAT_IF_DEBUG(cout << endl << "---- CompressAllocate ----" << endl << endl;)

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

void Encoder::FillCompressDeferred()
{
	CAT_IF_DEBUG(cout << endl << "---- FillCompressDeferred ----" << endl << endl;)

	// Initialize the deferred rows to zero
	memset(_ge_compress_matrix, 0, _defer_count * _ge_compress_pitch * sizeof(u64));

	// Fill the deferred rows from the peel matrix
	u16 row_i = _defer_head_rows;
	u64 *ge_compress_row = _ge_compress_matrix;
	while (row_i != LIST_TERM)
	{
		PeelRow *row = &_peel_rows[row_i];

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

		// Next row
		row_i = row->next;
		ge_compress_row += _ge_compress_pitch;
	}
}

void Encoder::FillCompressDense()
{
	CAT_IF_DEBUG(cout << endl << "---- FillCompressDense ----" << endl << endl;)

	// Fill the final few rows from the dense matrix:
	u64 *ge_compress_row = _ge_compress_matrix + _ge_compress_pitch * _defer_count;

	// Initialize PRNG
	CatsChoice prng;
	prng.Initialize(_g_seed, ~_block_count);

	// Fill the dense matrix with random bits
	// TODO: Optimize this for efficient compression
	int fill_count = _ge_compress_pitch * _added_count;
	while (fill_count--)
	{
		u32 rv1 = prng.Next();
		u32 rv2 = prng.Next();
		*ge_compress_row++ = ((u64)rv1 << 32) | rv2;
	}

	CAT_IF_DEBUG(
		if (_added_count > 64)
		{
			cout << "Cannot have added count > 64 right now" << endl;
		}
	)
}

void Encoder::FillGEDeferred()
{
	CAT_IF_DEBUG(cout << endl << "---- FillGEDeferred ----" << endl << endl;)

	// Fill the deferred rows from the mix matrix
	u16 row_i = _defer_head_rows;
	u64 *ge_row = _ge_matrix;
	u16 *ge_row_map = _ge_row_map;
	while (row_i != LIST_TERM)
	{
		PeelRow *row = &_peel_rows[row_i];

		CAT_IF_DEBUG(cout << "  Setting deferred row " << row_i << endl;)

		// For each mixing column in the row,
		u16 a = row->mix_a;
		u16 weight = row->mix_weight;
		u16 column_i = row->mix_x0;
		for (;;)
		{
			// Set bit for that column
			ge_row[column_i >> 6] |= (u64)1 << (column_i & 63);

			if (--weight <= 0) break;

			IterateNextColumn(column_i, _added_count, _added_next_prime, a);
		}

		// Remember which row was deferred into this GE matrix row
		*ge_row_map++ = row_i;

		// Next row
		row_i = row->next;
		ge_row += _ge_pitch;
	}
}

void Encoder::FillGEDense()
{
	CAT_IF_DEBUG(cout << endl << "---- FillGEDense ----" << endl << endl;)

	// Set the dense mixing bits to the identity matrix:
	u64 *ge_row = _ge_matrix + _ge_pitch * _defer_count;
	u64 ge_mask = 1;
	for (u16 column_i = 0; column_i < _added_count; ++column_i)
	{
		ge_row[column_i >> 6] |= ge_mask;

		CAT_IF_DEBUG(cout << "  Setting dense row " << column_i << endl;)

		ge_row += _ge_pitch;

		// Initialize column indices
		_ge_col_map[column_i] = _block_count + column_i;

		// NOTE: Not needed because these are known to be zero
		//_ge_row_map[_defer_count + column_i] = _block_count + column_i;

		// Generate next mask
		ge_mask = CAT_ROL64(ge_mask, 1);
	}
}

void Encoder::CopyDeferredColumns()
{
	CAT_IF_DEBUG(cout << endl << "---- CopyDeferredColumns ----" << endl << endl;)

	CAT_IF_DEBUG(cout << "Copy deferred columns list to pivots array:" << endl;)

	// For each deferred column,
	u16 column_i = _defer_head_columns;
	u16 ge_column_i = _added_count;
	while (column_i != LIST_TERM)
	{
		PeelColumn *column = &_peel_cols[column_i];

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

		// Iterate next column
		column_i = column->next;
		++ge_column_i;
	}
}

void Encoder::PrintGEMatrix()
{
	int rows = _defer_count + _added_count;
	int cols = rows;

	cout << endl << "GE matrix is " << rows << " x " << cols << ":" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			if (_ge_matrix[ii * _ge_pitch + (jj >> 6)] & ((u64)1 << (jj & 63)))
				cout << '1';
			else
				cout << '0';
		}
		cout << endl;
	}

	cout << endl;
}

void Encoder::PrintGECompressMatrix()
{
	int rows = _defer_count + _added_count;
	int cols = _block_count;

	cout << endl << "GE Compress matrix is " << rows << " x " << cols << ":" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			if (_ge_compress_matrix[ii * _ge_compress_pitch + (jj >> 6)] & ((u64)1 << (jj & 63)))
				cout << '1';
			else
				cout << '0';
		}
		cout << endl;
	}

	cout << endl;
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
	peeled column that was added.
*/

void Encoder::Compress()
{
	CAT_IF_DEBUG(cout << endl << "---- Compress ----" << endl << endl;)

	// Relink the peel list in reverse order
	u16 peel_head_rows = LIST_TERM;

	// For each column that has been peeled,
	u16 row_i = _peel_head_rows;
	u16 ge_rows = _defer_count + _added_count;
	while (row_i != LIST_TERM)
	{
		PeelRow *row = &_peel_rows[row_i];
		u16 column_i = row->unmarked[0];

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
					ge_row[_ge_pivots[ii] * _ge_compress_pitch] ^= ge_mask;

				if (--weight <= 0) break;

				IterateNextColumn(column_i, _block_count, _block_next_prime, a);
			}

			// For each mixing column in the row,
			a = row->mix_a;
			weight = row->mix_weight;
			column_i = row->mix_x0;
			for (;;)
			{
				// Set bit for that column in each row
				u64 *ge_row = _ge_matrix + (column_i >> 6);
				u64 ge_mask = (u64)1 << (column_i & 63);
				for (u16 ii = 0; ii < count; ++ii)
					ge_row[_ge_pivots[ii] * _ge_pitch] ^= ge_mask;

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

	CopyDeferredColumns();
}

bool Encoder::Triangle()
{
	CAT_IF_DEBUG(cout << endl << "---- Triangle ----" << endl << endl;)

	u16 ge_rows = _defer_count + _added_count;

	// Initialize pivot array
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
		_ge_pivots[pivot_i] = pivot_i;

	// For each pivot to determine,
	u16 ge_column_i = 0;
	u64 ge_mask = (u64)1 << (ge_column_i & 63);
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
	{
		int word_offset = ge_column_i >> 6;
		u64 *ge_row = _ge_matrix + word_offset;

		// For each remaining row,
		bool found = false;
		for (u16 pivot_row_i = pivot_i; pivot_row_i < ge_rows; ++pivot_row_i)
		{
			// Determine if the row contains the bit we want
			u16 ge_row_i = _ge_pivots[pivot_row_i];
			u64 *row = &ge_row[ge_row_i * _ge_pitch];

			// If the bit was found,
			if (*row & ge_mask)
			{
				found = true;
				CAT_IF_DEBUG(cout << "Pivot " << pivot_i << " found on row " << ge_row_i << endl;)

				// Swap out the pivot index for this one
				u16 temp = _ge_pivots[pivot_i];
				_ge_pivots[pivot_i] = _ge_pivots[pivot_row_i];
				_ge_pivots[pivot_row_i] = temp;

				// Prepare masked first word
				u64 row0 = (*row & ~(ge_mask - 1)) ^ ge_mask;

				// For each remaining row,
				++pivot_row_i;
				for (; pivot_row_i < ge_rows; ++pivot_row_i)
				{
					// Determine if the row contains the bit we want
					u16 rem_row_i = _ge_pivots[pivot_row_i];
					u64 *rem_row = &ge_row[rem_row_i * _ge_pitch];

					// If the bit was found,
					if (*rem_row & ge_mask)
					{
						// Add the pivot row to eliminate the bit from this row, preserving previous bits
						*rem_row ^= row0;
						for (int ii = 1; ii < _ge_pitch - word_offset; ++ii)
							rem_row[ii] ^= row[ii];
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

void Encoder::SolveTriangleColumns()
{
	CAT_IF_DEBUG(cout << endl << "---- SolveTriangleColumns ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// For each mix matrix pivot column,
	u16 ge_rows = _defer_count + _added_count;
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
	{
		u16 pivot_column_i = _ge_col_map[pivot_i];
		u16 pivot_ge_row_i = _ge_pivots[pivot_i];
		u8 *dest = _check_blocks + _block_bytes * pivot_column_i;

		CAT_IF_DEBUG(cout << "For compress row " << pivot_ge_row_i << " that solves column " << pivot_column_i << ":";)

		// If original row,
		if (pivot_ge_row_i < _defer_count)
		{
			u16 pivot_row_i = _ge_row_map[pivot_ge_row_i];
			const u8 *src = _message_blocks + _block_bytes * pivot_row_i;
			if (pivot_row_i == 387)
			{
				int x = 0;
			}

			CAT_IF_DEBUG(cout << " " << pivot_row_i << ":[" << (int)src[0] << "]";)

			// If copying from final block,
			if (pivot_row_i == _block_count - 1)
			{
				memcpy(dest, src, _final_bytes);
				memset(dest + _final_bytes, 0, _block_bytes - _final_bytes);
			}
			else
			{
				memcpy(dest, src, _block_bytes);
			}
		}
		else
		{
			CAT_IF_DEBUG(cout << " [0]";)

			memset(dest, 0, _block_bytes);
		}
		CAT_IF_ROWOP(++rowops;)

		// For each bit of the compress row,
		u64 *ge_compress_row = _ge_compress_matrix + _ge_compress_pitch * pivot_ge_row_i;
		for (u16 peel_column_i = 0; peel_column_i < _block_count; ++peel_column_i)
		{
			// If bit is set,
			u64 ge_mask = (u64)1 << (peel_column_i & 63);
			if (ge_compress_row[peel_column_i >> 6] & ge_mask)
			{
				CAT_IF_DEBUG(cout << " " << peel_column_i;)

				// Look up which peeled row solves this column
				PeelColumn *column = &_peel_cols[peel_column_i];
				u16 peel_row_i = column->w2_refs;

				// Look up source data
				const u8 *src = _message_blocks + _block_bytes * peel_row_i;
				int block_bytes = (peel_row_i == _block_count - 1) ? _final_bytes : _block_bytes;

				// TODO: Optimize this
				memxor(dest, src, block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}
		}

		CAT_IF_DEBUG(cout << endl << "  For GE row " << pivot_ge_row_i << " solved by pivot " << pivot_i << ":";)

		// For each bit of the GE matrix up to the diagonal,
		u64 *ge_row = _ge_matrix + _ge_pitch * pivot_ge_row_i;
		for (u16 ge_column_i = 0; ge_column_i < pivot_i; ++ge_column_i)
		{
			// If bit is set,
			u64 ge_mask = (u64)1 << (ge_column_i & 63);
			if (ge_row[ge_column_i >> 6] & ge_mask)
			{
				CAT_IF_DEBUG(cout << " " << ge_column_i;)

				// Look up source data
				u16 ge_pivot_column_i = _ge_col_map[ge_column_i];
				const u8 *src = _check_blocks + _block_bytes * ge_pivot_column_i;

				// TODO: Optimize this
				memxor(dest, src, _block_bytes);
				CAT_IF_DEBUG(cout << "[" << (int)src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}
		}

		CAT_IF_DEBUG(cout << endl;)
		CAT_IF_DEBUG(cout << "Result = " << (int)dest[0] << endl;)
	}

	CAT_IF_ROWOP(cout << "SolveTriangleColumns used " << rowops << " row ops" << endl;)
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
			if (ge_row[ge_above_row_i * _ge_pitch] & ge_mask)
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

void Encoder::Substitute()
{
	CAT_IF_DEBUG(cout << endl << "---- Substitute ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// For each column that has been peeled,
	u16 row_i = _peel_head_rows;
	u16 ge_rows = _defer_count + _added_count;
	while (row_i != LIST_TERM)
	{
		PeelRow *row = &_peel_rows[row_i];
		u16 dest_column_i = row->unmarked[0];
		u8 *dest = _check_blocks + _block_bytes * dest_column_i;

		CAT_IF_DEBUG(cout << "Generating column " << dest_column_i << ":";)

		const u8 *src = _message_blocks + _block_bytes * row_i;

		CAT_IF_DEBUG(cout << " " << row_i << ":[" << (int)src[0] << "]";)

		// If copying from final block,
		if (row_i == _block_count - 1)
		{
			memcpy(dest, src, _final_bytes);
			memset(dest + _final_bytes, 0, _block_bytes - _final_bytes);
		}
		else
		{
			memcpy(dest, src, _block_bytes);
		}
		CAT_IF_ROWOP(++rowops;)

		// For each peeling column in the row,
		u16 a = row->peel_a;
		u16 weight = row->peel_weight;
		u16 column_i = row->peel_x0;
		for (;;)
		{
			const u8 *src = _check_blocks + _block_bytes * column_i;

			CAT_IF_DEBUG(cout << " " << column_i;)

			// If column is not the solved one,
			if (column_i != dest_column_i)
			{
				// TODO: Optimize
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

		// For each mixing column in the row,
		a = row->mix_a;
		weight = row->mix_weight;
		u16 x = row->mix_x0;
		for (;;)
		{
			u16 column_i = _block_count + x;
			const u8 *src = _check_blocks + _block_bytes * column_i;

			CAT_IF_DEBUG(cout << " " << column_i;)

			// TODO: Optimize
			memxor(dest, src, _block_bytes);
			CAT_IF_DEBUG(cout << "[" << (int)src[0] << "]";)
			CAT_IF_ROWOP(++rowops;)

			if (--weight <= 0) break;

			IterateNextColumn(x, _added_count, _added_next_prime, a);
		}

		CAT_IF_DEBUG(cout << endl;)

		row_i = row->next;
	}

	CAT_IF_ROWOP(cout << "Substitute used " << rowops << " row ops" << endl;)
}

bool Encoder::GenerateCheckBlocks()
{
	CAT_IF_DEBUG(cout << endl << "---- GenerateCheckBlocks ----" << endl << endl;)

	// (1) Peeling

	if (!PeelSetup())
		return false;

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(
	{
		cout << "After PeelSetup, contents of peeled list (last to first):" << endl;
		u16 row_i = _peel_head_rows;
		while (row_i != LIST_TERM)
		{
			PeelRow *row = &_peel_rows[row_i];

			cout << "  Row " << row_i << " solves column " << row->unmarked[0] << endl;

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

			cout << "  Row " << row_i << " solves column " << row->unmarked[0] << endl;

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

	if (!CompressSetup())
		return false;

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(cout << "After CompressSetup:" << endl;)
	CAT_IF_DEBUG(PrintGEMatrix();)
	CAT_IF_DEBUG(PrintGECompressMatrix();)

	Compress();

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(cout << "After Compress:" << endl;)
	CAT_IF_DEBUG(PrintGEMatrix();)
	CAT_IF_DEBUG(PrintGECompressMatrix();)

	// (3) Gaussian Elimination

	if (!Triangle())
		return false;

	CAT_IF_DEBUG( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DEBUG(cout << "After Triangle:" << endl;)
	CAT_IF_DEBUG(PrintGEMatrix();)

	SolveTriangleColumns();

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

static const u16 PRIMES[] = {
	3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
	59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
	113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
	239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
	307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367,
	373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433,
	439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
	503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
	587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643,
	647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719,
	727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
	809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
	877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
	953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019,
	1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
	1087, 1091, 1093, 1097,
};

static int NextHighestPrime(int n)
{
	n |= 1;

	int p_max = (int)sqrt((float)n);

	for (;;)
	{
		const u16 *prime = PRIMES;

		for (;;)
		{
			int p = *prime;

			if (p > p_max)
				return n;

			if (n % p == 0)
				break;

			++prime;
		}

		n += 2;

		if (p_max * p_max < n)
			++p_max;
	}

	return n;
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
	_block_next_prime = NextHighestPrime(_block_count);
	_added_next_prime = NextHighestPrime(_added_count);

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
/*
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

			// Pad with zeros
			memset(buffer + _final_bytes, 0, _block_bytes - _final_bytes);
		}

		return;
	}
*/
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
