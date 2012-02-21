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
#include "WirehairUtil.hpp"
#include "SmallPRNG.hpp"
#include "memxor.hpp"
using namespace cat;
using namespace wirehair;

#if defined(CAT_DUMP_ENCODER_DEBUG)
#define CAT_IF_DUMP(x) x
#else
#define CAT_IF_DUMP(x)
#endif // CAT_DUMP_ENCODER_DEBUG

#if defined(CAT_DUMP_ROWOP_COUNTERS)
#define CAT_IF_ROWOP(x) x
#else
#define CAT_IF_ROWOP(x)
#endif

#if defined(CAT_DUMP_ENCODER_DEBUG) || defined(CAT_DUMP_ROWOP_COUNTERS) || defined(CAT_DUMP_GE_MATRIX)
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#endif // CAT_DUMP_ENCODER_DEBUG


//// Encoder

#pragma pack(push)
#pragma pack(1)
struct Encoder::PeelRow
{
	u16 next;					// Linkage in row list

	// Peeling matrix: Column generator
	u16 peel_weight, peel_a, peel_x0;

	// Mixing matrix: Column generator
	u16 mix_weight, mix_a, mix_x0;

	// Peeling state
	u16 unmarked_count;			// Count of columns that have not been marked yet
	union
	{
		// During peeling:
		struct
		{
			u16 unmarked[2];	// Final two unmarked column indices
		};

		// After peeling:
		struct
		{
			u16 peel_column;	// Peeling column that is solved by this row
			u8 is_copied;		// Row value is copied yet?
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
	CAT_IF_DUMP(cout << endl << "---- PeelSetup ----" << endl << endl;)

	CAT_IF_DUMP(cout << "Block Count = " << _block_count << endl;)

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
		row->peel_weight = GeneratePeelRowWeight(prng.Next(), _block_count - 1);
		u32 rv = prng.Next();
		row->peel_a = ((u16)rv % (_block_count - 1)) + 1;
		row->peel_x0 = (u16)(rv >> 16) % _block_count;

		// Generate mixing matrix row parameters
		row->mix_weight = mix_weight;
		rv = prng.Next();
		row->mix_a = ((u16)rv % (_added_count - 1)) + 1;
		row->mix_x0 = (u16)(rv >> 16) % _added_count;

		CAT_IF_DUMP(cout << "   Row " << row_i << " of weight " << row->peel_weight << " [a=" << row->peel_a << "] : ";)

		// Iterate columns in peeling matrix
		u16 weight = row->peel_weight;
		u16 column_i = row->peel_x0;
		u16 a = row->peel_a;
		u16 unmarked_count = 0;
		u16 unmarked[2];
		for (;;)
		{
			CAT_IF_DUMP(cout << column_i << " ";)

			PeelColumn *col = &_peel_cols[column_i];

			// Add row reference to column
			if (col->row_count >= ENCODER_REF_LIST_MAX)
			{
				CAT_IF_DUMP(cout << "PeelSetup: Failure!  Ran out of space for row references.  ENCODER_REF_LIST_MAX must be increased!" << endl;)
				return false;
			}
			col->rows[col->row_count++] = row_i;

			// If column is unmarked,
			if (col->mark == MARK_TODO)
				unmarked[unmarked_count++ & 1] = column_i;

			if (--weight <= 0) break;

			IterateNextColumn(column_i, _block_count, _block_next_prime, a);
		}
		CAT_IF_DUMP(cout << endl;)

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
				CAT_IF_DUMP(cout << "PeelAvalanche: Deferred(1) with column " << column_i << " at row " << ref_row_i << endl;)

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
					CAT_IF_DUMP(cout << "PeelAvalanche: Deferred(2) with column " << column_i << " at row " << ref_row_i << endl;)

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
	CAT_IF_DUMP(cout << "Peel: Solved column " << column_i << " with row " << row_i << endl;)

	PeelColumn *column = &_peel_cols[column_i];

	// Mark this column as solved
	column->mark = MARK_PEEL;

	// Remember which column it solves
	row->peel_column = column_i;

	// Link to back of the peeled list
	if (_peel_tail_rows != LIST_TERM)
		_peel_rows[_peel_tail_rows].next = row_i;
	else
		_peel_head_rows = row_i;
	row->next = LIST_TERM;
	_peel_tail_rows = row_i;

	// Indicate that this row hasn't been copied yet
	row->is_copied = 0;

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
	CAT_IF_DUMP(cout << endl << "---- GreedyPeeling ----" << endl << endl;)

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

		CAT_IF_DUMP(cout << "Deferred column " << best_column_i << " for Gaussian elimination, which had " << best_column->w2_refs << " weight-2 row references" << endl;)

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

bool Encoder::CompressSetup()
{
	CAT_IF_DUMP(cout << endl << "---- CompressSetup ----" << endl << endl;)

	// Allocate GE matrix
	int ge_rows = _defer_count + _added_count;
	int ge_pitch = (ge_rows + 63) / 64;
	int ge_matrix_words = ge_rows * ge_pitch;
	_ge_matrix = new u64[ge_matrix_words];
	if (!_ge_matrix) return false;
	_ge_pitch = ge_pitch;

	// Clear entire GE matrix
	memset(_ge_matrix, 0, ge_matrix_words * sizeof(u64));

	CAT_IF_DUMP(cout << "GE matrix is " << ge_rows << " x " << ge_rows << " with pitch " << ge_pitch << " consuming " << ge_matrix_words * sizeof(u64) << " bytes" << endl;)

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

	CAT_IF_DUMP(cout << "GE compress matrix is " << ge_compress_rows << " x " << ge_compress_columns << " with pitch " << ge_compress_pitch << " consuming " << ge_compress_matrix_words * sizeof(u64) << " bytes" << endl;)

	// Allocate the pivots
	int pivot_words = ge_rows * 3;
	_ge_pivots = new u16[pivot_words];
	if (!_ge_pivots) return false;
	_ge_col_map = _ge_pivots + ge_rows;
	_ge_row_map = _ge_col_map + ge_rows;

	CAT_IF_DUMP(cout << "Allocated " << ge_rows << " pivots, consuming " << pivot_words*2 << " bytes" << endl;)

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	return true;
}

void Encoder::SetDeferredColumns()
{
	CAT_IF_DUMP(cout << endl << "---- SetDeferredColumns ----" << endl << endl;)

	// For each deferred column,
	PeelColumn *column;
	for (u16 ge_column_i = 0, defer_i = _defer_head_columns; defer_i != LIST_TERM; defer_i = column->next, ++ge_column_i)
	{
		column = &_peel_cols[defer_i];

		CAT_IF_DUMP(cout << "GE column " << ge_column_i << " mapped to matrix column " << defer_i << " :";)

		// Set bit for each row affected by this deferred column
		u64 *matrix_row_offset = _ge_compress_matrix + (ge_column_i >> 6);
		u64 ge_mask = (u64)1 << (ge_column_i & 63);
		u16 count = column->row_count;
		u16 *ref_row = column->rows;
		while (count--)
		{
			u16 row_i = *ref_row++;

			CAT_IF_DUMP(cout << " " << row_i;)

			matrix_row_offset[_ge_compress_pitch * row_i] |= ge_mask;
		}

		CAT_IF_DUMP(cout << endl;)

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

		CAT_IF_DUMP(cout << "GE column(mix) " << ge_column_i << " mapped to matrix column " << column_i << endl;)

		_ge_col_map[ge_column_i] = column_i;
	}
}

void Encoder::SetMixingColumnsForDeferredRows()
{
	CAT_IF_DUMP(cout << endl << "---- SetMixingColumnsForDeferredRows ----" << endl << endl;)

	// For each deferred row,
	PeelRow *row;
	for (u16 defer_row_i = _defer_head_rows; defer_row_i != LIST_TERM; defer_row_i = row->next)
	{
		row = &_peel_rows[defer_row_i];

		CAT_IF_DUMP(cout << "Deferred row " << defer_row_i << " set mix columns :";)

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

			CAT_IF_DUMP(cout << " " << ge_column_i;)

			if (--weight <= 0) break;

			IterateNextColumn(x, _added_count, _added_next_prime, a);
		}

		CAT_IF_DUMP(cout << endl;)
	}
}

void Encoder::PeelDiagonal()
{
	/*
		This function optimizes the block value generation by combining the first
		memcpy and memxor operations together into a three-way memxor if possible,
		using the is_copied row member.
	*/

	// For each peeled row in forward solution order,
	PeelRow *row;
	for (u16 peel_row_i = _peel_head_rows; peel_row_i != LIST_TERM; peel_row_i = row->next)
	{
		row = &_peel_rows[peel_row_i];

		// Lookup peeling results
		u16 peel_column_i = row->peel_column;
		u64 *ge_row = _ge_compress_matrix + _ge_compress_pitch * peel_row_i;

		CAT_IF_DUMP(cout << "  Peeled row " << peel_row_i << " for peeled column " << peel_column_i << " :";)

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

			CAT_IF_DUMP(cout << " " << ge_column_i;)

			if (--weight <= 0) break;

			IterateNextColumn(x, _added_count, _added_next_prime, a);
		}

		CAT_IF_DUMP(cout << endl;)

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

			CAT_IF_DUMP(cout << "  -- Copied from " << peel_row_i << " because has not been copied yet.  Output block = " << (int)temp_block_src[0] << endl;)

			// NOTE: Do not need to set is_copied here because no further rows reference this one
		}

		// For each row that references this one,
		PeelColumn *column = &_peel_cols[peel_column_i];
		u16 count = column->row_count;
		u16 *ref_row = column->rows;
		while (count--)
		{
			u16 ref_row_i = *ref_row++;

			// Skip this row
			if (ref_row_i == peel_row_i) continue;

			CAT_IF_DUMP(cout << "  ++ Adding to referencing row " << ref_row_i << endl;)

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
}

void Encoder::CopyDeferredRows()
{
	CAT_IF_DUMP(cout << endl << "---- CopyDeferredRows ----" << endl << endl;)

	// For each deferred row,
	u64 *ge_row = _ge_matrix + _ge_pitch * _added_count;
	for (u16 ge_row_i = _added_count, defer_row_i = _defer_head_rows; defer_row_i != LIST_TERM;
		defer_row_i = _peel_rows[defer_row_i].next, ge_row += _ge_pitch, ++ge_row_i)
	{
		CAT_IF_DUMP(cout << "Peeled row " << defer_row_i << " for GE row " << ge_row_i << endl;)

		// Copy compress row to GE row
		u64 *compress_row = _ge_compress_matrix + _ge_compress_pitch * defer_row_i;
		memcpy(ge_row, compress_row, _ge_compress_pitch * sizeof(u64));

		// Set row map for this deferred row
		_ge_row_map[ge_row_i] = defer_row_i;
	}

	CAT_IF_DUMP(cout << "After copying deferred rows:" << endl;)
	CAT_IF_DUMP(PrintGEMatrix();)
}

void Encoder::MultiplyDenseRows()
{
	CAT_IF_DUMP(cout << endl << "---- MultiplyDenseRows ----" << endl << endl;)

	/*
		NOTE: This does not need to generate the same dense rows that
		the multiply-compression algorithm generates because both the
		encoder and decoder will make the same choice on which to use.

		TODO: Use 2-bit window optimization here also
	*/

	// Initialize PRNG
	CatsChoice prng;
	prng.Initialize(_g_seed, _block_count);

	// For each column,
	PeelColumn *column = _peel_cols;
	for (u16 column_i = 0; column_i < _block_count; ++column_i, ++column)
	{
		// Generate dense column
		u32 dense_rv = prng.Next();

		// Set up for light column generation
		u16 x = column_i % _light_count;
		u16 a = 1 + (column_i / _light_count) % (_light_count - 1);

		// If the column is peeled,
		if (column->mark == MARK_PEEL)
		{
			// Lookup GE compress matrix source row
			u16 source_row_i = column->peel_row;
			u64 *ge_source_row = _ge_compress_matrix + _ge_compress_pitch * source_row_i;
			u64 *ge_dest_row;

			CAT_IF_DUMP(cout << "For peeled column " << column_i << " solved by peel row " << source_row_i << " :";)

			// Light rows:
			ge_dest_row = _ge_matrix + _ge_pitch * x;
			for (int ii = 0; ii < _ge_compress_pitch; ++ii) ge_dest_row[ii] ^= ge_source_row[ii];
			CAT_IF_DUMP(cout << " " << x;)
			IterateNextColumn(x, _light_count, _light_next_prime, a);
			ge_dest_row = _ge_matrix + _ge_pitch * x;
			for (int ii = 0; ii < _ge_compress_pitch; ++ii) ge_dest_row[ii] ^= ge_source_row[ii];
			CAT_IF_DUMP(cout << " " << x;)
			IterateNextColumn(x, _light_count, _light_next_prime, a);
			ge_dest_row = _ge_matrix + _ge_pitch * x;
			for (int ii = 0; ii < _ge_compress_pitch; ++ii) ge_dest_row[ii] ^= ge_source_row[ii];
			CAT_IF_DUMP(cout << " " << x << ",";)

			// Dense rows:
			ge_dest_row = _ge_matrix + _ge_pitch * _light_count;
			for (u16 dense_i = 0; dense_i < _dense_count; ++dense_i, ge_dest_row += _ge_pitch, dense_rv >>= 1)
			{
				if (dense_rv & 1)
				{
					for (int ii = 0; ii < _ge_compress_pitch; ++ii) ge_dest_row[ii] ^= ge_source_row[ii];
					CAT_IF_DUMP(cout << " " << dense_i + _light_count;)
				}
			}
		}
		else
		{
			// Set up bit mask
			u16 ge_column_i = column->ge_column;
			u64 *ge_row = _ge_matrix + (ge_column_i >> 6);
			u64 ge_mask = (u64)1 << (ge_column_i & 63);

			CAT_IF_DUMP(cout << "For deferred column " << column_i << " at GE column " << ge_column_i << " :";)

			// Light rows:
			ge_row[_ge_pitch * x] ^= ge_mask;
			CAT_IF_DUMP(cout << " " << x;)
			IterateNextColumn(x, _light_count , _light_next_prime, a);
			ge_row[_ge_pitch * x] ^= ge_mask;
			CAT_IF_DUMP(cout << " " << x;)
			IterateNextColumn(x, _light_count , _light_next_prime, a);
			ge_row[_ge_pitch * x] ^= ge_mask;
			CAT_IF_DUMP(cout << " " << x << ",";)

			// Dense rows:
			ge_row += _ge_pitch * _light_count;
			for (u16 dense_i = 0; dense_i < _dense_count; ++dense_i, ge_row += _ge_pitch, dense_rv >>= 1)
			{
				if (dense_rv & 1)
				{
					*ge_row ^= ge_mask;
					CAT_IF_DUMP(cout << " " << dense_i + _light_count;)
				}
			}
		}

		CAT_IF_DUMP(cout << endl;)
	}
}

void Encoder::Compress()
{
	CAT_IF_DUMP(cout << endl << "---- Compress ----" << endl << endl;)

	/*
		(1) Peeling matrix inversion -> [ Defer | Mix ] rows + row values in peeling columns, check if rows are deferred
			(1) Set deferred column bits
			(2) Set mixing columns for each deferred row and mark them as deferred
			(3) Diagonalize the peeled rows
		(2) Copy lower deferred rows to GE matrix
		(3) Perform peeling matrix multiplication and generate GE rows, but no row values yet
	*/

	SetDeferredColumns();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	SetMixingColumnsForDeferredRows();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	PeelDiagonal();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	CopyDeferredRows();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	MultiplyDenseRows();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	AddInvertibleGF2Matrix(_ge_matrix, _defer_count, _ge_pitch, _added_count);

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )
}

void Encoder::SolveTriangleColumns()
{
	CAT_IF_DUMP(cout << endl << "---- SolveTriangleColumns ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// (1) Initialize each check block value to the solution row value

	// For each pivot,
	const u16 ge_rows = _defer_count + _added_count;
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
	{
		u16 column_i = _ge_col_map[pivot_i];
		u16 ge_row_i = _ge_pivots[pivot_i];
		u8 *buffer_dest = _check_blocks + _block_bytes * column_i;

		CAT_IF_DUMP(cout << "Pivot " << pivot_i << " solving column " << column_i << " with GE row " << ge_row_i << " : ";)

		if (ge_row_i < _added_count)
		{
			memset(buffer_dest, 0, _block_bytes);

			// Store which column solves the dense row
			_ge_row_map[ge_row_i] = column_i;

			CAT_IF_DUMP(cout << "[0]";)
			CAT_IF_ROWOP(++rowops;)
		}
		else
		{
			u16 pivot_row_i = _ge_row_map[ge_row_i];
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
			CAT_IF_ROWOP(++rowops;)

			CAT_IF_DUMP(cout << "[" << (int)buffer_src[0] << "]";)

			// Eliminate peeled columns:
			PeelRow *row = &_peel_rows[pivot_row_i];
			u16 column_i = row->peel_x0;
			u16 a = row->peel_a;
			u16 weight = row->peel_weight;
			for (;;)
			{
				PeelColumn *column = &_peel_cols[column_i];
				if (column->mark == MARK_PEEL)
				{
					memxor(buffer_dest, _check_blocks + _block_bytes * column_i, _block_bytes);
					CAT_IF_ROWOP(++rowops;)
				}

				if (--weight <= 0) break;

				IterateNextColumn(column_i, _block_count, _block_next_prime, a);
			}
		}
		CAT_IF_DUMP(cout << endl;)
	}

	// (2) Add dense rows
	// TODO: Use 2-bit window optimization here also

	// Initialize PRNG
	CatsChoice prng;
	prng.Initialize(_g_seed, _block_count);

	// For each column,
	const u8 *source_block = _check_blocks;
	PeelColumn *column = _peel_cols;
	for (u16 column_i = 0; column_i < _block_count; ++column_i, ++column, source_block += _block_bytes)
	{
		// Generate dense column
		u32 dense_rv = prng.Next();

		// If the column is peeled,
		if (column->mark == MARK_PEEL)
		{
			CAT_IF_DUMP(cout << "Peeled column " << column_i << "[" << (int)source_block[0] << "] :";)

			// Light rows:
			u16 x = column_i % _light_count;
			u16 a = 1 + (column_i / _light_count) % (_light_count - 1);

			CAT_IF_DUMP(cout << " " << x;)
			memxor(_check_blocks + _block_bytes * _ge_row_map[x], source_block, _block_bytes);
			CAT_IF_ROWOP(++rowops;)
			IterateNextColumn(x, _light_count, _light_next_prime, a);
			CAT_IF_DUMP(cout << " " << x;)
			memxor(_check_blocks + _block_bytes * _ge_row_map[x], source_block, _block_bytes);
			CAT_IF_ROWOP(++rowops;)
			IterateNextColumn(x, _light_count, _light_next_prime, a);
			CAT_IF_DUMP(cout << " " << x;)
			memxor(_check_blocks + _block_bytes * _ge_row_map[x], source_block, _block_bytes);
			CAT_IF_ROWOP(++rowops;)

			// Dense rows:
			for (u16 ii = _light_count; ii < _added_count; ++ii, dense_rv >>= 1)
			{
				if (dense_rv & 1)
				{
					CAT_IF_DUMP(cout << " " << ii;)
					memxor(_check_blocks + _block_bytes * _ge_row_map[ii], source_block, _block_bytes);
					CAT_IF_ROWOP(++rowops;)
				}
			}

			CAT_IF_DUMP(cout << endl;)
		} // end if peeled
	}

	// (3) For each GE matrix bit up to the diagonal, add deferred and mixed:

	CAT_IF_ROWOP(cout << "SolveTriangleColumns add dense rows used " << rowops << " row ops" << endl;)
	CAT_IF_ROWOP(rowops = 0;)

	// For each pivot,
	for (u16 pivot_i = 0; pivot_i < ge_rows; ++pivot_i)
	{
		u16 pivot_column_i = _ge_col_map[pivot_i];
		u16 ge_row_i = _ge_pivots[pivot_i];
		u8 *buffer_dest = _check_blocks + _block_bytes * pivot_column_i;

		CAT_IF_DUMP(cout << "Pivot " << pivot_i << " solving column " << pivot_column_i << "[" << (int)buffer_dest[0] << "] with GE row " << ge_row_i << " :";)

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
				CAT_IF_ROWOP(++rowops;)

				CAT_IF_DUMP(cout << " " << column_i << "=[" << (int)peel_src[0] << "]";)
			}

			ge_mask = CAT_ROL64(ge_mask, 1);
		}

		CAT_IF_DUMP(cout << endl;)
	}

	CAT_IF_ROWOP(cout << "SolveTriangleColumns under diagonal adding used " << rowops << " row ops" << endl;)
}


//// Diagnostic functions for compression algorithms

#if defined(CAT_DUMP_ENCODER_DEBUG) || defined(CAT_DUMP_GE_MATRIX)

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

void Encoder::PrintGECompressMatrix()
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

#endif // CAT_DUMP_ENCODER_DEBUG


/*
	One more subtle optimization.  Depending on how the GE matrix
	is constructed, it can be put in a roughly upper-triangular form
	from the start so it looks like this:

		+---------+---------+
		| D D D D | D D D 1 |
		| D D D D | D D 1 M |
		| d d d d | D 1 M M |
		| d d d d | 1 M M M |
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

	Since some of the added check rows are light rows and others are
	dense rows, which appear first in the matrix must also be decided.
	In the above diagram, the lowercase 'd' rows represent the dense
	rows and the uppercase 'D' rows represent the light rows.  This
	way GE will prefer to select light rows, which will cause less
	mixing down the added rows, leading to fewer row ops.  I haven't
	quantified this optimization but it is clear from looking at the
	GE matrix after Triangle() that this row order is preferable.

	This is not a huge improvement but every little bit helps.
*/


//// (3) Gaussian Elimination

bool Encoder::Triangle()
{
	CAT_IF_DUMP(cout << endl << "---- Triangle ----" << endl << endl;)

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
				CAT_IF_DUMP(cout << "Pivot " << pivot_i << " found on row " << ge_row_j << endl;)

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
			CAT_IF_DUMP(cout << "Inversion impossible: Pivot " << pivot_i << " not found!" << endl;)
			return false;
		}

		// Generate next mask
		ge_mask = CAT_ROL64(ge_mask, 1);
	}

	return true;
}

void Encoder::Diagonal()
{
	CAT_IF_DUMP(cout << endl << "---- Diagonal ----" << endl << endl;)

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

		CAT_IF_DUMP(cout << "Pivot " << pivot_i << " solving column " << pivot_column_i << ":";)

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

				CAT_IF_DUMP(cout << " " << above_column_i;)

				memxor(dest, src, _block_bytes);
				CAT_IF_DUMP(cout << "[" << (int)src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}
		}

		CAT_IF_DUMP(cout << endl;)

		// Generate next mask
		ge_mask = CAT_ROR64(ge_mask, 1);
	}

	CAT_IF_ROWOP(cout << "Diagonal used " << rowops << " row ops" << endl;)
}


//// (4) Substitute

void Encoder::Substitute()
{
	CAT_IF_DUMP(cout << endl << "---- Substitute ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// For each column that has been peeled,
	u16 ge_rows = _defer_count + _added_count;
	PeelRow *row;
	for (u16 row_i = _peel_head_rows; row_i != LIST_TERM; row_i = row->next)
	{
		row = &_peel_rows[row_i];
		u16 dest_column_i = row->peel_column;
		u8 *dest = _check_blocks + _block_bytes * dest_column_i;

		CAT_IF_DUMP(cout << "Generating column " << dest_column_i << ":";)

		const u8 *combo, *src = _message_blocks + _block_bytes * row_i;

		CAT_IF_DUMP(cout << " " << row_i << ":[" << (int)src[0] << "]";)

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

			CAT_IF_DUMP(cout << " " << column_i;)

			if (!combo)
				memxor(dest, src, _block_bytes);
			else
			{
				memxor(dest, src, combo, _block_bytes);
				combo = 0;
			}
			CAT_IF_DUMP(cout << "[" << (int)src[0] << "]";)
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

			CAT_IF_DUMP(cout << " " << column_i;)

			// If column is not the solved one,
			if (column_i != dest_column_i)
			{
				memxor(dest, src, _block_bytes);
				CAT_IF_DUMP(cout << "[" << (int)src[0] << "]";)
				CAT_IF_ROWOP(++rowops;)
			}
			else
			{
				CAT_IF_DUMP(cout << "*";)
			}

			if (--weight <= 0) break;

			IterateNextColumn(column_i, _block_count, _block_next_prime, a);
		}

		CAT_IF_DUMP(cout << endl;)
	}

	CAT_IF_ROWOP(cout << "Substitute used " << rowops << " row ops" << endl;)
}


//// Main Driver

bool Encoder::GenerateCheckBlocks()
{
	CAT_IF_DUMP(cout << endl << "---- GenerateCheckBlocks ----" << endl << endl;)

	// (1) Peeling

	if (!PeelSetup())
		return false;

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DUMP(
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

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DUMP(
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

	CAT_IF_DUMP(
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

	CAT_IF_DUMP(
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

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	CAT_IF_DUMP(cout << "After CompressSetup:" << endl;)
	CAT_IF_DUMP(PrintGEMatrix();)
	CAT_IF_DUMP(PrintGECompressMatrix();)

	Compress();

	CAT_IF_DUMP(cout << "After Compress:" << endl;)
#if defined(CAT_DUMP_ENCODER_DEBUG) || defined(CAT_DUMP_GE_MATRIX)
	cout << "After Compress:" << endl;
	PrintGEMatrix();
#endif
	CAT_IF_DUMP(PrintGECompressMatrix();)

	// (3) Gaussian Elimination

	if (!Triangle())
	{
		CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

		CAT_IF_DUMP(cout << "After Triangle FAILED:" << endl;)
		CAT_IF_DUMP(PrintGEMatrix();)

		return false;
	}

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

#if defined(CAT_DUMP_ENCODER_DEBUG) || defined(CAT_DUMP_GE_MATRIX)
	cout << "After Triangle:" << endl;
	PrintGEMatrix();
#endif

	SolveTriangleColumns();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	Diagonal();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

	// (4) Substitution

	Substitute();

	CAT_IF_DUMP( _ASSERTE( _CrtCheckMemory( ) ); )

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
	CAT_IF_DUMP(cout << endl << "---- Initialize ----" << endl << endl;)

	Cleanup();

	// Calculate message block count
	int block_count = (message_bytes + block_bytes - 1) / block_bytes;
	_final_bytes = message_bytes % block_bytes;
	if (_final_bytes <= 0) _final_bytes = block_bytes;

	// Lookup generator matrix parameters
	if (!GenerateMatrixParameters(block_count, _g_seed, _light_count, _dense_count))
		return false;
	_added_count = _light_count + _dense_count;

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
	_light_next_prime = NextPrime16(_light_count);

	CAT_IF_DUMP(cout << "Total message = " << message_bytes << " bytes" << endl;)
	CAT_IF_DUMP(cout << "Block bytes = " << block_bytes << ".  Final bytes = " << _final_bytes << endl;)
	CAT_IF_DUMP(cout << "Block count = " << block_count << " +Prime=" << _block_next_prime << endl;)
	CAT_IF_DUMP(cout << "Added count = " << _added_count << " +Prime=" << _added_next_prime << endl;)
	CAT_IF_DUMP(cout << "Generator seed = " << _g_seed << endl;)
	CAT_IF_DUMP(cout << "Memory overhead for check blocks = " << check_size << " bytes" << endl;)

	// Generate check blocks
	return GenerateCheckBlocks();
}

void Encoder::Generate(u32 id, void *block)
{
	u8 *buffer = reinterpret_cast<u8*>( block );

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

	CAT_IF_DUMP(cout << "Generating row " << id << ":";)

	// Initialize PRNG
	CatsChoice prng;
	prng.Initialize(id, _g_seed);

	// Generate peeling matrix row parameters
	u16 peel_weight = GeneratePeelRowWeight(prng.Next(), _block_count - 1);
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

	CAT_IF_DUMP(cout << " " << peel_x;)

	// If peeler has multiple columns,
	if (peel_weight > 1)
	{
		--peel_weight;

		IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

		CAT_IF_DUMP(cout << " " << peel_x;)

		// Combine first two columns into output buffer (faster than memcpy + memxor)
		memxor(buffer, first, _check_blocks + _block_bytes * peel_x, _block_bytes);

		// For each remaining peeler column,
		while (--peel_weight > 0)
		{
			IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

			CAT_IF_DUMP(cout << " " << peel_x;)

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

	CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

	// For each remaining mixer column,
	while (--mix_weight > 0)
	{
		IterateNextColumn(mix_x, _added_count, _added_next_prime, mix_a);

		CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

		// Mix in each column
		memxor(buffer, _check_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	}

	CAT_IF_DUMP(cout << endl;)
}
