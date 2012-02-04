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

static u32 GetGeneratorSeed(int block_count)
{
	// TODO: Needs to be simulated (2)
	return 0;
}

static int GetCheckBlockCount(int block_count)
{
	// TODO: Needs to be simulated (1)
	return 8;
}


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

// During generator matrix pre-computation, tune this to be as small as possible and still succeed
static const int ENCODER_REF_LIST_MAX = 8;

#pragma pack(push)
#pragma pack(1)
struct Encoder::PeelColumn
{
	u16 w2_refs;	// Number of weight-2 rows containing this column
	u16 row_count;	// Number of rows containing this column
	u16 rows[ENCODER_REF_LIST_MAX];

	u16 next;		// Linkage in column list
};
#pragma pack(pop)

void Encoder::Peel(u16 row_i, PeelRow *row, u16 column_i)
{
	// Remember which column it solves
	row->unmarked[0] = column_i;

	// Link to front of the peeled list
	row->next = _peel_head;
	_peel_head = row_i;

	// Walk list of peeled rows referenced by this newly solved column
	PeelColumn *col = &_peel_cols[column_i];
	u16 ref_count = col->row_count;
	u16 *ref_rows = col->rows;
	while (ref_count--)
	{
		// Update unmarked row count for this referenced row
		u16 ref_row_i = *ref_rows++;
		PeelRow *ref_row = &_peel_rows[ref_row_i];
		u16 unmarked_count = --row->unmarked_count;

		// If row may be solving a column now,
		if (unmarked_count == 1)
		{
			// Find other column
			u16 new_column_i = ref_row->unmarked[0];
			if (new_column_i == column_i)
				new_column_i = ref_row->unmarked[1];

			// If column is already solved,
			if (_peel_cols[new_column_i].w2_refs == LIST_TERM)
			{
				// TODO: Mark row as deferred
			}
			else
			{
				// TODO: Peel!
			}
		}
		else if (unmarked_count == 2)
		{
			// Regenerate the row columns to discover which are unmarked
			u16 ref_weight = ref_row->peel_weight;
			u16 ref_column_i = ref_row->peel_x0;
			u16 ref_a = ref_row->mix_a;
			u16 unmarked_count = 0;
			while (ref_weight--)
			{
				PeelColumn *ref_col = &_peel_cols[ref_column_i];

				// Generate next column
				do ref_column_i = (ref_column_i + ref_a) % _block_next_prime;
				while (ref_column_i >= _block_count);

				// If column is unmarked,
				if (col->w2_refs != LIST_TERM)
				{
					// Store the two unmarked columns in the row
					ref_row->unmarked[unmarked_count++] = ref_column_i;
				}
			}
		}
	}
}

bool Encoder::PeelSetup()
{
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
	}

	// Declare a list of peeled rows in reverse-solution order
	_peel_head = LIST_TERM;

	// Generate peeling row data
	for (u16 row_i = 0; row_i < _block_count; ++row_i)
	{
		PeelRow *row = &_peel_rows[row_i];

		// Initialize PRNG
		CatsChoice prng;
		prng.Initialize(row_i, _g_seed);

		// Generate peeling matrix row parameters
		row->peel_weight = GenerateWeight(prng.Next());
		u32 rv = prng.Next();
		row->peel_a = ((u16)rv % (_block_count - 1)) + 1;
		row->peel_x0 = (u16)(rv >> 16) % _block_count;

		// Generate mixing matrix row parameters
		row->mix_weight = 3;
		rv = prng.Next();
		row->mix_a = ((u16)rv % (_added_count - 1)) + 1;
		row->mix_x0 = (u16)(rv >> 16) % _added_count;

		// Iterate columns in peeling matrix
		u16 weight = row->peel_weight;
		u16 column_i = row->peel_x0;
		u16 a = row->mix_a;
		u16 unmarked_count = 0;
		u16 unmarked[2];
		while (weight--)
		{
			PeelColumn *col = &_peel_cols[column_i];

			// Add row reference to column
			col->rows[col->row_count++] = row_i;

			// Generate next column
			do column_i = (column_i + a) % _block_next_prime;
			while (column_i >= _block_count);

			// If column is peeled,
			if (col->w2_refs != LIST_TERM)
				unmarked[unmarked_count++ & 1] = column_i;
		}

		// Generate column list
		row->unmarked_count = unmarked_count;

		// If this row solves a column,
		if (unmarked_count == 1)
		{
			Peel(row_i, row, unmarked[0]);
		}
		else if (unmarked_count == 2)
		{
			// Remember which two columns were unmarked
			row->unmarked[0] = unmarked[0];
			row->unmarked[1] = unmarked[1];

			// Increment weight-2 reference count for unmarked columns
			_peel_cols[unmarked[0]].w2_refs++;
			_peel_cols[unmarked[1]].w2_refs++;
		}
	}
}

void Encoder::GreedyPeeling()
{

}

void Encoder::Compress()
{

}

bool Encoder::Triangle()
{

}

void Encoder::Diagonal()
{

}

void Encoder::Substitute()
{

}

bool Encoder::GenerateCheckBlocks()
{
	if (!PeelSetup())
		return false;

	GreedyPeeling();

	Compress();

	if (!Triangle())
		return false;

	Diagonal();

	Substitute();

	return true;
}

bool Encoder::Initialize(const void *message_in, int message_bytes, int block_bytes)
{
	CAT_IF_DEBUG( if (block_bytes < 16) return false; )

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

	// For the message blocks,
	if (id < _block_count)
	{
		// Until the final block in message blocks,
		if (id < _block_bytes - 1)
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

	// Set up generators
	PeelGenerator peeler;
	MixGenerator mixer;
	_generator.Make(_g_seed, id, peeler, mixer);

	// Remember first column (there is always at least one)
	u8 *first = _check_blocks + _block_bytes * peeler.x;

	// If peeler has multiple columns,
	int weight = peeler.weight;
	if (weight > 1)
	{
		--weight;

		// Combine first two columns into output buffer (faster than memcpy + memxor)
		peeler.Next();
		memxor(buffer, first, _check_blocks + _block_bytes * peeler.x, _block_bytes);

		// For each remaining peeler column,
		while (--weight > 0)
		{
			// Mix in each column
			peeler.Next();
			memxor(buffer, _check_blocks + _block_bytes * peeler.x, _block_bytes);
		}

		// Mix first mixer block in directly
		memxor(buffer, _check_blocks + _block_bytes * (_block_count + mixer.x), _block_bytes);
	}
	else
	{
		// Mix first with first mixer block (faster than memcpy + memxor)
		memxor(buffer, first, _check_blocks + _block_bytes * (_block_count + mixer.x), _block_bytes);
	}

	// For each remaining mixer column,
	weight = mixer.weight;
	while (--weight > 0)
	{
		// Mix in each column
		peeler.Next();
		memxor(buffer, _check_blocks + _block_bytes * (_block_count + mixer.x), _block_bytes);
	}
}


//// Decoder

bool Decoder::Initialize(void *message_out, int message_bytes, int block_bytes)
{
}

bool Decoder::Decode(void *block)
{
}
