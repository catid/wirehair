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

#ifndef CAT_WIREHAIR_DETAILS_HPP
#define CAT_WIREHAIR_DETAILS_HPP

// TODO: Remove this
extern int g_p_seed, g_c_seed;

#include "SmallPRNG.hpp"

// Switches:
//#define CAT_DUMP_CODEC_DEBUG /* Turn on debug output for decoder */
//#define CAT_DUMP_ROWOP_COUNTERS /* Dump row operations counters to console */
//#define CAT_DUMP_GE_MATRIX /* Dump GE matrix to console */
//#define CAT_LIGHT_ROWS /* Use light rows for all check columns (slower) */
//#define CAT_REUSE_COMPRESS /* Reuse the compression matrix for back-substitution (slower) */
#define CAT_COPY_FIRST_N /* Copy the first N rows from the input (faster) */
#define CAT_SHUFFLE_HALF /* Reshuffle second half of check rows from a new starting point (slower) */
#define CAT_WINDOWED_BACKSUB /* Use window optimization for back-substitution (faster) */
#define CAT_REF_LIST_MAX 32 /* Tune to be as small as possible and still succeed */
#define CAT_MAX_CHECK_ROWS 1024 /* Maximum check row count */
#define CAT_MAX_EXTRA_ROWS 32 /* Maximum number of extra rows to support before reusing existing rows */

namespace cat {

namespace wirehair {


//// Result object

enum Result
{
	R_WIN,				// Operation: Success!
	R_MORE_BLOCKS,		// Codec wants more blocks

	R_ERROR,			// Return codes higher than this one are errors:
	R_BAD_CHECK_SEED,	// Encoder needs a better check seed
	R_BAD_PEEL_SEED,	// Encoder needs a better peel seed
	R_BAD_INPUT,		// Input parameters were incorrect
	R_TOO_SMALL,		// message_bytes / block_size is too small.  Try reducing block_size or use a larger message
	R_NEED_MORE_EXTRA,	// Not enough extra rows to solve it, must give up
	R_OUT_OF_MEMORY,	// Out of memory, try reducing the message size
};

// Get Result String function
const char *GetResultString(Result r);


//// Utilities

// 16-bit Integer Square Root function
u16 SquareRoot16(u16 x);

// 16-bit Truncated Sieve of Eratosthenes Next Prime function
u16 NextPrime16(u16 n);

// Peeling Row Weight Generator function
u16 GeneratePeelRowWeight(u32 rv, u16 max_weight);

// GF(2) Invertible Matrix Generator function
bool AddInvertibleGF2Matrix(u64 *matrix, int offset, int pitch, int n);

// Matrix Parameter Generator function
bool GenerateMatrixParameters(int block_count, u32 &p_seed, u32 &c_seed, u16 &light_count, u16 &dense_count);

// Deck Shuffling function
void ShuffleDeck16(CatsChoice &prng, u16 *deck, u32 count);

// Peel Matrix Row Generator function
void GeneratePeelRow(u32 id, u32 p_seed, u16 peel_column_count, u16 mix_column_count,
	u16 &peel_weight, u16 &peel_a, u16 &peel_x0, u16 &mix_a, u16 &mix_x0);


//// Utility: Column iterator function

/*
	This implements a very light PRNG (Weyl function) to quickly generate
	a set of random-looking columns without replacement.

	This is Stewart Platt's excellent loop-less iterator optimization.
	His common cases all require no additional modulus operation, which
	makes it faster than the rare case that I designed.
*/

CAT_INLINE void IterateNextColumn(u16 &x, u16 b, u16 p, u16 a)
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


//// Encoder/Decoder Combined Implementation

class Codec
{
	// Parameters
	u32 _block_bytes;			// Number of bytes in a block
	u16 _block_count;			// Number of blocks in the message
	u16 _extra_count;			// Number of extra rows to allocate
	u32 _p_seed;				// Seed for peeled rows of generator matrix
	u32 _c_seed;				// Seed for check rows of generator matrix
	u16 _block_next_prime;		// Next prime number at or above block count
	u16 _added_next_prime;		// Next prime number at or above added count
	u16 _light_next_prime;		// Next prime number at or above light count
	u16 _used_count;			// Number of stored rows
	u16 _light_count;			// Number of code rows that are light
	u16 _dense_count;			// Number of code rows that are dense
	u16 _added_count;			// Sum of light and dense row counts
	u8 *_recovery_blocks;		// Recovery blocks
	u8 *_input_blocks;			// Input message blocks
	u32 _input_final_bytes;		// Number of bytes in final block of input
	u32 _output_final_bytes;	// Number of bytes in final block of output
	bool _input_allocated;		// Flag indicating that the input buffer was allocated

	// Peeling state
	struct PeelRow;
	struct PeelColumn;
	struct PeelRefs;
	PeelRow *_peel_rows;		// Array of N peeling matrix rows
	PeelColumn *_peel_cols;		// Array of N peeling matrix columns
	PeelRefs *_peel_col_refs;	// List of column references
	PeelRow *_peel_tail_rows;	// Tail of peeling solved rows list
	static const u16 LIST_TERM = 0xffff;
	u16 _peel_head_rows;		// Head of peeling solved rows list
	u16 _defer_head_columns;	// Head of peeling deferred columns list
	u16 _defer_head_rows;		// Head of peeling deferred rows list
	u16 _defer_count;			// Count of deferred rows

	// Gaussian elimination state
	u64 *_ge_matrix;			// Gaussian elimination matrix
	u64 *_ge_compress_matrix;	// Gaussian elimination compression matrix
	int _ge_pitch;				// Words per row of GE matrix and compression matrix
	u16 _ge_rows;				// Number of rows in GE matrix, since this grows
	u16 *_ge_pivots;			// Pivots for each column of the GE matrix
	u16 *_ge_col_map;			// Map of GE columns to conceptual matrix columns
	u16 *_ge_row_map;			// Map of GE rows to conceptual matrix rows
	u16 _ge_resume_pivot;		// Pivot to resume Triangle() on after it fails
#if defined(CAT_REUSE_COMPRESS)
	u8 *_win_table_data;		// Values of temporary symbols for substitution table
#endif

#if defined(CAT_DUMP_CODEC_DEBUG) || defined(CAT_DUMP_GE_MATRIX)
	void PrintGEMatrix();
	void PrintCompressMatrix();
	void PrintPeeled();
	void PrintDeferredRows();
	void PrintDeferredColumns();
#endif


	//// (1) Peeling

	// Avalanche peeling from the newly solved column to others
	void PeelAvalanche(u16 column_i);

	// Peel a row using the given column
	void Peel(u16 row_i, PeelRow *row, u16 column_i);

	// Walk forward through rows and solve as many as possible before deferring any
	bool OpportunisticPeeling(u32 row_i, u32 id);

	// Greedy algorithm to select columns to defer and resume peeling until all columns are marked
	void GreedyPeeling();


	//// (2) Compression

	// Set deferred column bits in compression matrix
	void SetDeferredColumns();

	// Set mixing columns for deferred rows
	void SetMixingColumnsForDeferredRows();

	// Diagonalize the peeling matrix, generating compression matrix
	void PeelDiagonal();

	// Copy deferred rows from the compress matrix to the GE matrix
	void CopyDeferredRows();

	// Multiply dense rows by peeling matrix to generate GE rows, but no row values yet
	void MultiplyDenseRows();


	//// (3) Gaussian Elimination

	// Triangularize the GE matrix (may fail if pivot cannot be found)
	bool Triangle();

	// Initialize column values for GE matrix
	void InitializeColumnValues();

	// Add check matrix to triangle columns
	void AddCheckValues();

	// Add values for GE matrix positions under the diagonal
	void AddSubdiagonalValues();


	//// (4) Substitution

#if defined(CAT_REUSE_COMPRESS)
	// Reuse compression matrix to substitute values (best for small N)
	void CompressionBasedSubstitute();
#endif

	// Back-substitute to diagonalize the GE matrix
	void BackSubstituteAboveDiagonal();

	// Regenerate all of the sparse peeled rows to diagonalize them
	void Substitute();


	//// Main Driver

	// Choose matrix to use based on message bytes
	Result ChooseMatrix(int message_bytes, int block_bytes);

	// Solve matrix so that recovery blocks can be generated
	Result SolveMatrix();

	// Resume solver with a new block
	Result ResumeSolveMatrix(u32 id, const void *block);

	// Use matrix solution to generate recovery blocks
	void GenerateRecoveryBlocks();


	//// Memory Management

	void SetInput(const void *message_in);
	bool AllocateInput();
	void FreeInput();

	bool AllocateMatrix();
	void FreeMatrix();

	bool AllocateWorkspace();
	void FreeWorkspace();

public:
	Codec();
	~Codec();


	//// Accessors

	CAT_INLINE u32 PSeed() { return _p_seed; }
	CAT_INLINE u32 CSeed() { return _c_seed; }
	CAT_INLINE u32 BlockCount() { return _block_count; }


	//// Encoder Mode

	// Initialize encoder mode
	Result InitializeEncoder(int message_bytes, int block_bytes);

	// Feed encoder a message
	Result EncodeFeed(const void *message_in);

	// Encode a block
	void Encode(u32 id, void *block_out);


	//// Decoder Mode

	// Initialize decoder mode
	Result InitializeDecoder(int message_bytes, int block_bytes);

	// Feed decoder a block
	Result DecodeFeed(u32 id, const void *block_in);

	// Generate output blocks from the recovered check blocks
	void ReconstructOutput(void *message_out);
};


} // namespace wirehair

} // namespace cat

#endif // CAT_WIREHAIR_DETAILS_HPP
