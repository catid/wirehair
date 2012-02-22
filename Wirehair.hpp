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

#ifndef CAT_WIREHAIR_HPP
#define CAT_WIREHAIR_HPP

/*
	Wirehair Streaming Forward Error Correction

		Wirehair is an FEC codec used to improve reliability of data sent
	over a Binary Erasure Channel (BEC) such as satellite Internet or WiFi.
	The data to transmit over the BEC is encoded into equal-sized blocks.
	When enough blocks are received at the other end of the channel, then
	the original data can be recovered.

	Block format:

		M = Count of bytes per transmission over the unreliable channel.
		ID = 24-bit identifier for each block.

		+----+-----------------------------+
		| ID |   BYTES OF BLOCK DATA (M)   |
		+----+-----------------------------+

		The identifier allows the decoder to synchronize with the encoder
	even when the blocks arrive out of order or packet loss causes a
	block to be dropped.

	Encoding Setup:

		T = Size of original data in bytes.
		N = ceil(T / M) = Count of blocks in the original data.

		(1) Generator Matrix Construction

			A = Original data blocks, N blocks long.
			H = Count of dense matrix rows (see below), chosen based on N.
			E = N + H blocks = Count of encoded data blocks.
			B = Encoded data blocks, E blocks long.
			G = Generator matrix, with E rows and E columns.

			+---------+---+   +---+   +---+
			|         |   |   |   |   |   |
			|    P    | M |   |   |   | A |
			|         |   | x | B | = |   |
			+---------+---+   |   |   +---+
			|    D    | I |   |   |   | 0 |
			+---------+---+   +---+   +---+

			A and B are Ex1 vectors of blocks.
				A has N rows of the original data padded by H zeroes.
				B has E rows of encoded blocks.

			G is the ExE binary matrix on the left.
				P is the NxN peeling matrix
					- Optimized for success of the peeling solver.
				M is the NxH mixing matrix
					- Used to mix the H dense rows into the N peeling rows.
				D is the HxN dense matrix
					- Used to improve recovery properties.
				I is the HxH identity matrix.

			G matrices for each value of N are precomputed offline and used
			based on the length of the input data, which guarantees that G
			is invertible.

		(2) Generating Matrix P

			The Hamming weight of each row of P is a random variable with a
			distribution chosen to optimize the operation of the peeling
			solver (see below).
			For each row of the matrix, this weight is determined and 1 bits
			are then uniformly distributed over the N columns.

		(3) Generating Matrix M

			Rows of M are generated with a constant weight of 2 and 1 bits are
			uniformly distributed over the H columns.

		(4) Generating Matrix D

			Each bit has a 50% chance of being set.

		(5) Generator Matrix Inversion

			An optimized sparse inversion technique is used to solve for B.

	---------------------------------------------------------------------------
	Sparse Matrix Inversion:

		There are 4 phases to this sparse inversion:

		(1) Peeling
			- Opportunistic fast solution for first N rows.
		(2) Compression
			- Setup for Gaussian elimination on a wide rectangular matrix
		(3) Gaussian Elimination
			- Gaussian elimination on a (hopefully) small square matrix
		(4) Substitution
			- Solves for remaining rows from initial peeling

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

		After N rows, an unmarked column is marked Deferred and the count of all
	rows that reference that column are decremented by 1.  Peeling avalanche may
	occur as a result.  This process is repeated until all columns are marked.

	A comment on the matrix state after peeling:

		After the peeling solver has completed, only Deferred columns remain to be
	solved.  Conceptually the matrix can be re-ordered in the order of solution
	so that the matrix resembles this example:

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

		The next step is to compress the GE rows and columns into a smaller GE
	matrix that GE will be applied to.  Note that up to this point, no row ops
	have been performed.  GE row values will be calculated in this step for use
	during GE.

		For each peeled column from the last to the first,
			For each GE row that contains that column,
				Add the corresponding peeled row to it.
			Next
		Next

		The resultant conceptual matrix may resemble this example:

		+-----+---------+-----+
		|   1 | 0 0 0 0 | 285 |
		| 1 1 | 0 0 0 0 | 712 |
		+-----+---------+-----+
		   ^       ^       ^---- Row values
		   |       \------------ Peeled matrix (rectangular)
		   \-------------------- GE square matrix

		The density of the rows is indicated in the example because it is worth
	considering.  Some of the rows may be low density, so mixing them before
	calculating the new row value is undesirable.  By calculating the row values
	at this point, row operations may be saved.

	(3) Gaussian Elimination:

		In this step the GE square matrix is put in upper-triangular form and
	initial GE square matrix row values are calculated.  Back-substitution
	will be performed afterwards to complete solving for the row values.

		For each GE matrix column from left to right,
			For each GE matrix row from top to bottom,
				If GE matrix bit is set,
					Set this row as the pivot for this column.

					For each remaining GE matrix row containing this column,
						Sum other GE square matrix rows with the pivot row,
							including the row values.
					Next
				End
			Next

			If no pivot was found for this column,
				The process fails and needs to resume on the next block.
			End
		Next

		At successful completion of this process, all pivots have been found.
	The GE square matrix is now in upper triangular form.

		To complete solving for the GE row values, the upper triangular part
	of the GE square matrix must be eliminated.  This process will diagonalize
	the GE square matrix, conceptually leaving it as an identity matrix.

		For each pivot column from last to first,
			For each pivot row above it,
				If the pivot row contains the column,
					Add the pivot column value to the pivot row value
				End
			Next
		Next

		Note that the GE square matrix is no longer needed and the pivot row
	values are entirely determined.  Example result:

		+-----+
		| 001 | <- Pivot 1
		| 002 | <- Pivot 2
		+-----+

	(4) Substitution:

		Finally the peeled row values may be determined by regenerating the
	whole matrix and summing up columns in the original order of peeling.

		For each peeled row from first to last,
			For each active column in the row,
				Add corresponding row value to peeled row.
			Next
		Next

	Review:

		After all of these steps, the row values have been determined and the
	matrix inversion is complete.  Let's analyze the complexity of each step:

		(1) Peeling
			- Opportunistic fast solution for first N rows.

			Weight determination : O(k) average
				Column reference update : Amortized O(1) for each column

			If peeling activates,
				Marking as peeled : O(1)

				Reducing weight of rows referencing this column : O(k)
					If other row weight is reduced to 2,
						Regenerate columns and mark potential Deferred : O(k)
					End
			End

			So peeling is O(1) for each row, and O(N) overall.

		(2) Compression
			- Setup for Gaussian elimination on a wide rectangular matrix

			Summing peeled rows into GE rows : O(k*N*N) bit operations

			This algorithm is fairly complex in terms of real-world time
			spent calculating it.  I expect about 10% of the run time to
			be spent executing this code.  However the total execution
			time is dominated by the number of row sum operations (large
			XOR operations) so this non-linear asymptotic performance is
			not a problem.

			Summing to produce the row values : O(S*N) = O(N^1.5) row ops

		(3) Gaussian Elimination
			- Gaussian elimination on a (hopefully) small square matrix

			Assume the GE square matrix is SxS, and S = sqrt(N) on
			average thanks to the peeling solver above.

			Gaussian elimination : O(S^3) = O(N^1.5) bit operations

			This algorithm is not bad because the matrix is small and
			it really doesn't contribute much to the run time.

			- Solves for rows of small square matrix

			Assume the GE square matrix is SxS, and S = sqrt(N) on
			average thanks to the peeling solver above.

			Back-substitute inside the GE matrix : O(S^2) = O(N) row ops

		(4) Substitution
			- Solves for remaining rows from initial peeling

			Regenerate peeled matrix rows and substitute : O(N*k) row ops

		As compared with matrix multiplication, this algorithm has:

		+ Memory overhead for the peeling solver state.
		+ O(N^2) bit operations for Compression.
		+ O(N + N^1.5) heavy row ops for GE.

		It gains back a little because it is not spending additional row ops
		to solve row values involved in GE.  The final substitution step is
		the same as matrix multiplication for those rows (same as the encoder).

	---------------------------------------------------------------------------

	Encoding:

			The first N output blocks of the encoder are the same as the
		original data.  After that the encoder will start producing random-
		looking M-byte blocks by generating new rows for P and M and
		multiplying them by B.

	Decoding:

			Decoding begins by collecting N blocks from the transmitter.  Once
		N blocks are received, the matrix G' (differing in the first N rows
		from the above matrix G) is generated with the rows of P|M that were
		received.  Generator matrix inversion is attempted, failing at the
		Gaussian elimination step if a pivot cannot be found for one of the GE
		matrix columns (see above).

			New rows are received and submitted directly to the GE solver,
		hopefully providing the missing pivot.  Once enough rows have been
		received, back-substitution reconstructs matrix B.
		
			The first N rows of the original matrix G are then used to fill in
		any blocks that were not received from the original N blocks, and the
		original data is received.
*/

#include "Platform.hpp"

// Switches:
//#define CAT_DUMP_ENCODER_DEBUG /* Turn on debug output for encoder */
//#define CAT_DUMP_DECODER_DEBUG /* Turn on debug output for decoder */
//#define CAT_DUMP_ROWOP_COUNTERS /* Dump row operations counters to console */
//#define CAT_DUMP_GE_MATRIX /* Dump GE matrix to console */
//#define CAT_ENCODER_COPY_FIRST_N /* Copy the first N rows from the input (much faster) */
#define CAT_SHUFFLE_DENSE_ROWS /* Use shuffle function to add dense rows (much faster) */

namespace cat {

namespace wirehair {


/*
	Wirehair Encoder

	Encodes message blocks for transmission over the network.
	The initialization function takes a while (say 10 milliseconds), so I
	recommend performing initialization in a separate thread to take advantage
	of modern multi-core processors.

	Example usage pseudocode:

		cat::wirehair::Encoder encoder;

		encoder.Initialize(file_data, file_bytes, 1500);

		while (file not received on other end)
		{
			u8 buffer[1500];

			encoder.Generate(buffer);

			udp_send(buffer);
		}
*/
class Encoder
{
	// Check block state
	u32 _block_bytes;	// Number of bytes in a block
	u32 _final_bytes;	// Number of bytes in final block
	u16 _block_count;	// Number of blocks in the message
	u16 _light_count;	// Number of check rows that are light
	u16 _dense_count;	// Number of check rows that are dense
	u16 _added_count;	// Number of check blocks added overall
	u8 *_check_blocks;	// Pointer to start of check blocks
	u32 _g_seed;		// Seed for nonsingular generator matrix

	// Encoder state
	const u8 *_message_blocks;	// Original message data (final block is partial)
	u32 _next_block_id;			// Next block identifier to transmit
	u16 _block_next_prime;		// Next prime number including or higher than block count
	u16 _added_next_prime;		// Next prime number including or higher than added count
	u16 _light_next_prime;		// Next prime number including or higher than light count

	// Peeling state
	struct PeelRow;
	struct PeelColumn;
	PeelRow *_peel_rows;		// Array of N peeling matrix rows
	PeelColumn *_peel_cols;		// Array of N peeling matrix columns
	PeelRow *_peel_tail_rows;	// Tail of peeling solved rows list

	// Lists
	static const u16 LIST_TERM = 0xffff;
	u16 _peel_head_rows;		// Head of peeling solved rows list
	u16 _defer_head_columns;	// Head of peeling deferred columns list
	u16 _defer_head_rows;		// Head of peeling deferred rows list
	u16 _defer_count;			// Count of deferred rows

	// Gaussian elimination state
	u64 *_ge_matrix;			// Gaussian elimination matrix
	int _ge_pitch;				// Pitch in words of GE matrix and compression matrix
	u16 *_ge_pivots;			// Pivots for each column of the GE matrix
	u16 *_ge_col_map;			// Map of GE columns to check matrix columns
	u64 *_ge_compress_matrix;	// Gaussian elimination compression matrix
	u16 *_ge_row_map;			// Map of GE rows to check matrix rows

#if defined(CAT_DUMP_ENCODER_DEBUG) || defined(CAT_DUMP_GE_MATRIX)
	void PrintGEMatrix();
	void PrintGECompressMatrix();
#endif


	//// (1) Peeling

	// Avalanche peeling from the newly solved column to others
	void PeelAvalanche(u16 column_i, PeelColumn *column);

	// Peel a row using the given column
	void Peel(u16 row_i, PeelRow *row, u16 column_i);

	// Walk forward through rows and solve as many as possible before deferring any
	bool PeelSetup();

	// Greedy algorithm to select columns to defer and resume peeling until all columns are marked
	void GreedyPeeling();


	//// (2) Compression

	// Build GE matrix for compression
	bool CompressSetup();

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

	// Compress rectangular matrix into conceptual square matrix
	void Compress();

	// Solve pivot column values from the row op schedule from Triangle
	void SolveTriangleColumns();


	//// (3) Gaussian Elimination

	// Triangularize the GE matrix (may fail if pivot cannot be found)
	bool Triangle();

	// Diagonalize the GE matrix to complete solving for the GE blocks
	void Diagonal();


	//// (4) Substitution

	// Substitute and solve for all of the peeled columns
	void Substitute();


	//// Misc

	// Main driver: Generate check blocks from message blocks
	bool GenerateCheckBlocks();

	// Free allocated memory
	void Cleanup();

public:
	Encoder();
	~Encoder();

	CAT_INLINE u32 GetSeed() { return _g_seed; }
	CAT_INLINE u32 GetBlockCount() { return _block_count; }

	// Attempt to initialize the encoder for a given message
	// Returns false on initialization failure
	bool Initialize(const void *message_in, int message_bytes, int block_bytes);

	// Generate a block of size block_bytes specified during initialization
	// id = Which block to generate, first N are same as original input
	void Generate(u32 id, void *block);
};


/*
	Wirehair Decoder

	Decodes messages encoded by the Encoder above.
	The Initialize() function does not take much run time.
	The Decode() function will return true when decoding has completed.

	Example usage pseudocode:

		cat::wirehair::Decoder decoder;

		decoder.Initialize(out_file_buffer, file_bytes, 1500);

		do
		{
			u8 buffer[1500];

			udp_recv(buffer);

		} while (!decoder.Decode(buffer));
*/
class Decoder
{
	// Check block state
	u32 _block_bytes;	// Number of bytes in a block
	u32 _final_bytes;	// Number of bytes in final block
	u16 _block_count;	// Number of blocks in the message
	u16 _used_count;	// Number of blocks used
	u16 _alloc_count;	// Number of blocks allocated above block count
	u16 _light_count;	// Number of check rows that are light
	u16 _dense_count;	// Number of check rows that are dense
	u16 _added_count;	// Number of check blocks added overall
	u8 *_check_blocks;	// Pointer to start of check blocks
	u32 _g_seed;		// Seed for nonsingular generator matrix

	// Encoder state
	u8 *_message_blocks;	// Original message data (final block is partial)
	u32 _next_block_id;		// Next block identifier to transmit
	u16 _block_next_prime;	// Next prime number including or higher than block count
	u16 _added_next_prime;	// Next prime number including or higher than added count
	u16 _light_next_prime;	// Next prime number including or higher than light count

	// Peeling state
	struct PeelRow;
	struct PeelColumn;
	PeelRow *_peel_rows;		// Array of N peeling matrix rows
	PeelColumn *_peel_cols;		// Array of N peeling matrix columns
	PeelRow *_peel_tail_rows;	// Tail of peeling solved rows list

	// Lists
	static const u16 LIST_TERM = 0xffff;
	u16 _peel_head_rows;		// Head of peeling solved rows list
	u16 _defer_head_columns;	// Head of peeling deferred columns list
	u16 _defer_head_rows;		// Head of peeling deferred rows list
	u16 _defer_count;			// Count of deferred rows

	// Gaussian elimination state
	u64 *_ge_matrix;			// Gaussian elimination matrix
	int _ge_pitch;				// Pitch in words of GE matrix and compression matrix
	u16 *_ge_pivots;			// Pivots for each column of the GE matrix
	u16 *_ge_col_map;			// Map of GE columns to check matrix columns
	u64 *_ge_compress_matrix;	// Gaussian elimination compression matrix
	u16 *_ge_row_map;			// Map of GE rows to check matrix rows

#if defined(CAT_DUMP_DECODER_DEBUG) || defined(CAT_DUMP_GE_MATRIX)
	void PrintGEMatrix();
	void PrintGECompressMatrix();
#endif


	//// (1) Peeling

	// Avalanche peeling from the newly solved column to others
	void PeelAvalanche(u16 column_i, PeelColumn *column);

	// Peel a row using the given column
	void Peel(u16 row_i, PeelRow *row, u16 column_i);

	// Walk forward through rows and solve as many as possible before deferring any
	bool PeelSetup();

	// Greedy algorithm to select columns to defer and resume peeling until all columns are marked
	void GreedyPeeling();


	//// (2) Compression

	// Build GE matrix for compression
	bool CompressSetup();

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

	// Compress rectangular matrix into conceptual square matrix
	void Compress();

	// Solve pivot column values from the row op schedule from Triangle
	void SolveTriangleColumns();


	//// (3) Gaussian Elimination

	// Triangularize the GE matrix (may fail if pivot cannot be found)
	bool Triangle();

	// Diagonalize the GE matrix to complete solving for the GE blocks
	void Diagonal();


	//// (4) Substitution

	// Substitute and solve for all of the peeled columns
	void Substitute();


	//// Misc

	// Main driver: Generate check blocks from message blocks
	bool GenerateCheckBlocks();

	// Recreate message with the received blocks and the check blocks
	void RecreateMessage();

	// Resume GE with latest information
	bool GEResume(u32 id, const u8 *buffer);

	// Free allocated memory
	void Cleanup();

public:
	Decoder();
	~Decoder();

	CAT_INLINE u32 GetSeed() { return _g_seed; }
	CAT_INLINE u32 GetBlockCount() { return _block_count; }

	// Attempt to initialize the decoder
	// Decoder will write to the given buffer once decoding completes
	// Returns false on initialization failure
	bool Initialize(void *message_out, int message_bytes, int block_bytes);

	// Decode the provided block of size block_bytes specified during initialization
	// id = ID of provided block, first N are same as original message
	bool Decode(u32 id, const void *block);
};


} // namespace wirehair

} // namespace cat

#endif // CAT_WIREHAIR_HPP
