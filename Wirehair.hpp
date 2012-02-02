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

#include "Platform.hpp"

namespace cat {


/*
	Wirehair Streaming Forward Error Correction

		Wirehair is a FEC codec used to improve reliability of data sent
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
				P is the NxN peeling matrix optimized for success of the peeling solver.
				M is the NxH mixing matrix used to mix the H dense rows into the N peeling rows.
				D is the HxN dense matrix used to improve recovery properties.
				I is the HxH identity matrix.

			G matrices for each value of N are precomputed offline and used based on the length
			of the input data, which guarantees that G is invertible.

		(2) Generating Matrix P

			The Hamming weight of each row of P is a random variable with a distribution chosen to
			optimize the operation of the peeling decoder (see below).
			For each row of the matrix, this weight is determined and 1 bits are then uniformly
			distributed over the N columns.

		(3) Generating Matrix M

			Rows of M are generated with a constant weight of 2 and 1 bits are uniformly distributed
			over the H columns.

		(4) Generating Matrix D

			Each bit has a 50% chance of being set.

		(5) Generator Matrix Inversion

			Peeling and Gaussian elimination are used to solve for B given A.

			(a) Iteratively Peeling Matrix P

				Peeling proceeds only on submatrix P of G.
				During peeling, any row with a weight of 1 is added to the solution for that column.
				As a row is added to the solution, the weight of any other row containing that column
				is reduced by 1.

				Eventually peeling will halt because no row has a weight of 1.
				When this happens, a column is selected for Gaussian elimination and is set aside,
				and peeling then resumes with the remaining unsolved rows.

			(b) Gaussian Elimination Setup

				Any row of matrix G that references a peeled column is summed to eliminate peeled
				columns from G.  The remaining rows and columns are used for Gaussian elimination
				and form a GE matrix.

			(c) Gaussian Elimination

				The GE matrix is scanned from top to bottom looking for a 1 bit set in each row.
				The first 1 bit found is set as pivot for the current column, and remaining 1 bits
				found cause row summing.
				If a column is empty then GE fails.
				After all columns are processed the GE matrix is an upper triangular matrix.

			(d) Gaussian Elimination Back-Substitution

				Back-substitution proceeds from final pivot row to first pivot row, summing the
				pivot row with the previous pivot rows as needed to make the GE matrix an identity
				matrix.

			(e) Peeling Back-Substitution

				Peeled rows are summed with GE-solved rows to eliminate remaining active columns.
				The entire G matrix is now the LxL identity matrix and the row blocks correspond
				to the blocks of B.

	Encoding:

			The first N output blocks of the encoder are the same as the original data.
			After that the encoder will start producing random-looking M-byte blocks by
		generating new rows for P and M and multiplying them by B.

	Decoding:

			Decoding begins by collecting N blocks from the transmitter.  Once N blocks are
		received, the matrix G' (differing in the first N rows from the above matrix G) is
		generated with the rows of P|M that were received.  Generator matrix inversion is
		attempted, failing at the Gaussian Elimination step if a pivot cannot be found for
		one of the GE matrix columns (see above).

			New rows are received and submitted directly to the GE solver, hopefully providing
		the missing pivot.  Once enough rows have been received, back-substitution reconstructs
		matrix B.
		
			The first N rows of the original matrix G are then used to fill in any blocks
		that were not received from the original N blocks, and the original data is received.
*/

class WirehairEncoder
{
public:
	// Attempt to initialize the encoder for a given message
	// Returns false on initialization failure
	bool Initialize(const void *message_in, int message_bytes, int block_bytes);

	// Generate a block of size block_bytes specified during initialization
	void Generate(void *block);
};

class WirehairDecoder
{
public:
	// Attempt to initialize the decoder
	// Decoder will write to the given buffer once decoding completes
	// Returns false on initialization failure
	bool Initialize(void *message_out, int message_bytes, int block_bytes);

	// Decode the provided block of size block_bytes specified during initialization
	bool Decode(void *block);
};


} // namespace cat

#endif // CAT_WIREHAIR_HPP
