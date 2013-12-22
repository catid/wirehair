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

/*
	Mathematical Overview:

		S = Size of original data in bytes.
		N = ceil(S / M) = Count of blocks in the original data.

		(1) Check Matrix Construction

			A = Original data blocks, N blocks long.
			D = Count of dense/heavy matrix rows (see below), chosen based on N.
			E = N + D blocks = Count of encoded data blocks.
			R = Recovery blocks, E blocks long.
			C = Check matrix, with E rows and E columns.
			0 = Dense/heavy rows sum to zero.

			+---------+-------+   +---+   +---+
			|         |       |   |   |   |   |
			|    P    |   M   |   |   |   | A |
			|         |       |   |   |   |   |
			+---------+---+---+ x | R | = +---+
			|    D    | J | 0 |   |   |   | 0 |
			+---------+-+-+-+-+   |   |   +---+
			|    0      | H |I|   |   |   | 0 |
			+-----------+---+-+   +---+   +---+

			A and B are Ex1 vectors of blocks.
				A has N rows of the original data padded by H zeroes.
				R has E rows of encoded blocks.

			C is the ExE binary matrix on the left.
				P is the NxN peeling binary matrix
					- Optimized for success of the peeling solver.
				M is the NxD mixing binary matrix
					- Used to mix the D dense/heavy rows into the peeling rows.
				D is the DxN dense binary matrix
					- Used to improve on recovery properties of peeling code.
				J is a DxD random-looking invertible matrix.
				H is the 6x12 heavy byte matrix
					- Used to improve on recovery properties of dense code.
				I is a 6x6 identity byte matrix.

			C matrices for each value of N are precomputed off-line and used
			based on the length of the input data, which guarantees that C
			is invertible.  The I matrix is also selected based on its size.

		(2) Generating Matrix P

			The Hamming weight of each row of P is a random variable with a
			distribution chosen to optimize the operation of the peeling
			solver (see below).
			For each row of the matrix, this weight is determined and 1 bits
			are then uniformly distributed over the N columns.

		(3) Generating Matrix M

			Rows of M are generated with a constant weight of 3 and 1 bits are
			uniformly distributed over the H columns.

		(4) Generating Matrix D

			Matrix D is generated with a Shuffle-2 Code, a novel invention.
			Shuffle-2 codes produce random matrices that offer possibly the
			fastest matrix-matrix multiplication algorithm for this purpose.
			Each bit has a 50% chance of being set.

		(5) Generating Matrix H

			The heavy matrix H is also a novel invention in this context.
			Adding these rows greatly improves invertibility of C while providing
			a constant-time algorithm to solve them.
			H is a 6x12 random byte matrix.
			Each element of H is a byte instead of a bit, representing a number
			in GF(2^^8) with generator polynomial 0x15F.

		(6) Check Matrix Solver

			An optimized sparse technique is used to solve the recovery blocks:

	---------------------------------------------------------------------------
	Sparse Matrix Solver:

	There are 4 phases to this sparse inversion:

	(1) Peeling
		- Opportunistic fast solution for first N rows.
	(2) Compression
		- Setup for Gaussian elimination
	(3) Gaussian Elimination
		- Gaussian elimination on a small square matrix
	(4) Substitution
		- Solves for remaining rows from initial peeling

	See the code comments in Wirehair.cpp for documentation of each step.

		After all of these steps, the row values have been determined and the
	matrix solver is complete.  Let's analyze the complexity of each step:

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

		The dense row multiplication takes O(N / 2 + ceil(N / D) * 2 * (D - 1))
		where D is approximately SQRT(N), so dense row multiplication takes:
		O(N / 2 + ceil(SQRT(N)) * SQRT(N)) = O(1.5*N).

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

		Solving inside the GE matrix : O(S^2) = O(N) row ops

	(4) Substitution
		- Solves for remaining rows from initial peeling

		Regenerate peeled matrix rows and substitute : O(N*k) row ops

		So overall, the codec scales roughly linearly in row operations,
	meaning that the throughput is somewhat stable over a wide number of N.

	---------------------------------------------------------------------------

	Encoding:

		The first N output blocks of the encoder are the same as the
	original data.  After that the encoder will start producing random-
	looking M-byte blocks by generating new rows for P and M and
	multiplying them by B.

	Decoding:

		Decoding begins by collecting N blocks from the transmitter.  Once
	N blocks are received, the matrix C' (differing in the first N rows
	from the above matrix C) is generated with the rows of P|M that were
	received.  Matrix solving is attempted, failing at the Gaussian
	elimination step if a pivot cannot be found for one of the GE matrix
	columns (see above).

		New rows are received and submitted directly to the GE solver,
	hopefully providing the missing pivot.  Once enough rows have been
	received, back-substitution reconstructs matrix B.

		The first N rows of the original matrix G are then used to fill in
	any blocks that were not received from the original N blocks, and the
	original data is recovered.
*/

#include "Wirehair.hpp"
#include "MemXOR.hpp"
#if defined(CAT_HEAVY_WIN_MULT)
#include "EndianNeutral.hpp"
#endif
using namespace cat;
using namespace wirehair;


//// Precompiler-conditional console output

#if defined(CAT_DUMP_PIVOT_FAIL)
#define CAT_IF_PIVOT(x) x
#else
#define CAT_IF_PIVOT(x)
#endif

#if defined(CAT_DUMP_CODEC_DEBUG)
#define CAT_IF_DUMP(x) x
#else
#define CAT_IF_DUMP(x)
#endif

#if defined(CAT_DUMP_ROWOP_COUNTERS)
#define CAT_IF_ROWOP(x) x
#else
#define CAT_IF_ROWOP(x)
#endif

#if defined(CAT_DUMP_CODEC_DEBUG) || defined(CAT_DUMP_PIVOT_FAIL) || \
	defined(CAT_DUMP_ROWOP_COUNTERS) || defined(CAT_DUMP_GE_MATRIX)
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#endif


//// Utility: Get Result String function

const char *cat::wirehair::GetResultString(Result r)
{
	/*
		See comments in Result enumeration in WirehairDetails.hpp
		for suggestions on how to overcome these errors.
	*/
	switch (r)
	{
	case R_WIN:				return "R_WIN";
	case R_MORE_BLOCKS:		return "R_MORE_BLOCKS";
	case R_BAD_DENSE_SEED:	return "R_BAD_DENSE_SEED";
	case R_BAD_PEEL_SEED:	return "R_BAD_PEEL_SEED";
	case R_TOO_SMALL:		return "R_TOO_SMALL";
	case R_TOO_LARGE:		return "R_TOO_LARGE";
	case R_NEED_MORE_EXTRA:	return "R_NEED_MORE_EXTRA";
	case R_BAD_INPUT:		return "R_BAD_INPUT";
	case R_OUT_OF_MEMORY:	return "R_OUT_OF_MEMORY";

	default:				if (r >= R_ERROR) return "R_UNKNOWN_ERROR";
							else return "R_UNKNOWN";
	}
}


//// Utilities

// 16-bit Integer Square Root function
static u16 SquareRoot16(u16 x);

// 16-bit Truncated Sieve of Eratosthenes Next Prime function
static u16 NextPrime16(u16 n);

// Peeling Row Weight Generator function
static u16 GeneratePeelRowWeight(u32 rv, u16 peel_column_count);

// GF(2) Invertible Matrix Generator function
static bool AddInvertibleGF2Matrix(u64 * CAT_RESTRICT matrix, int offset, int pitch, int n);

// Deck Shuffling function
static void ShuffleDeck16(Abyssinian &prng, u16 * CAT_RESTRICT deck, u32 count);

// Peel Matrix Row Generator function
static void GeneratePeelRow(u32 id, u32 p_seed, u16 peel_column_count, u16 mix_column_count,
	u16 & CAT_RESTRICT peel_weight, u16 & CAT_RESTRICT peel_a, u16 & CAT_RESTRICT peel_x0,
	u16 & CAT_RESTRICT mix_a, u16 & CAT_RESTRICT mix_x0);


//// Utility: 16-bit Integer Square Root function

/*
	Based on code from http://www.azillionmonkeys.com/qed/sqroot.html

		"Contributors include Arne Steinarson for the basic approximation idea, 
		Dann Corbit and Mathew Hendry for the first cut at the algorithm, 
		Lawrence Kirby for the rearrangement, improvments and range optimization
		and Paul Hsieh for the round-then-adjust idea."

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

u16 SquareRoot16(u16 x)
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
	int p_max = SquareRoot16(n);

	// For each number to try,
	for (;;)
	{
		// For each prime to test up to square root,
		const u16 * CAT_RESTRICT prime = PRIMES_UNDER_256;
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


//// Utility: GF(2) Invertible Matrix Generator function

/*
	Sometimes it is helpful to be able to quickly generate a GF2 matrix
	that is invertible.  I guess.  Anyway, here's a lookup table of
	seeds that create invertible GF2 matrices and a function that will
	fill a bitfield with the matrix.

	The function generates random-looking invertible matrices for
		0 < N < 512
	And for larger values of N it will just add the identity matrix.

	It will add the generated matrix rather than overwrite what was there.
*/

static const u8 INVERTIBLE_MATRIX_SEEDS[512] = {
	0x0,0,2,2,10,5,6,1,2,0,0,3,5,0,0,1,0,0,0,3,0,1,2,3,0,1,6,6,1,6,0,0,
	0,4,2,7,0,2,4,2,1,1,0,0,2,12,11,3,3,3,2,1,1,4,4,1,13,2,2,1,3,2,1,1,
	3,1,0,0,1,0,0,10,8,6,0,7,3,0,1,1,0,2,6,3,2,2,1,0,5,2,5,1,1,2,4,1,
	2,1,0,0,0,2,0,5,9,17,5,1,2,2,5,4,4,4,4,4,1,2,2,2,1,0,1,0,3,2,2,0,
	1,4,1,3,1,17,3,0,0,0,0,2,2,0,0,0,1,11,4,2,4,2,1,8,2,1,1,2,6,3,0,4,
	3,10,5,3,3,1,0,1,2,6,10,10,6,0,0,0,0,0,0,1,4,2,1,2,2,12,2,2,4,0,0,2,
	0,7,12,1,1,1,0,6,8,0,0,0,0,2,1,8,6,2,0,5,4,2,7,2,10,4,2,6,4,6,6,1,
	0,0,0,0,3,1,0,4,2,6,1,1,4,2,5,1,4,1,0,0,1,8,0,0,6,0,17,4,9,8,4,4,
	3,0,0,3,1,4,3,3,0,0,3,0,0,0,3,4,4,4,3,0,0,12,1,1,2,5,8,4,8,6,2,2,
	0,0,0,13,0,3,4,2,2,1,6,13,3,12,0,0,3,7,8,2,2,2,0,0,4,0,0,0,2,0,3,6,
	7,1,0,2,2,4,4,3,6,3,6,4,4,1,3,7,1,0,0,0,1,3,0,5,4,4,4,3,1,1,7,13,
	4,6,1,1,2,2,2,5,7,1,0,0,2,2,1,2,1,6,6,6,2,2,2,5,3,2,0,0,0,0,0,0,
	0,0,2,3,2,2,0,4,0,0,4,2,0,0,0,2,4,1,2,3,1,1,1,1,1,1,1,1,4,0,0,0,
	1,1,0,0,0,0,0,4,3,0,0,0,0,4,0,0,4,5,2,0,1,0,0,1,7,1,0,0,0,0,1,1,
	1,6,3,0,0,1,3,2,0,3,0,2,1,1,1,0,0,0,0,0,0,8,0,0,6,4,1,3,5,3,0,1,
	1,6,3,3,5,2,2,9,5,1,2,2,1,1,1,1,1,1,2,2,1,3,1,0,0,4,1,7,0,0,0,0,
};

static bool AddInvertibleGF2Matrix(u64 * CAT_RESTRICT matrix, int offset, int pitch, int n)
{
	if (n <= 0) return false;

	// If we have this value of n in the table,
	if (n < 512)
	{
		// Pull a random matrix out of the lookup table
		Abyssinian prng;
		prng.Initialize(INVERTIBLE_MATRIX_SEEDS[n]);

		// If shift is friendly,
		u32 shift = offset & 63;
		u64 * CAT_RESTRICT row = matrix + (offset >> 6);
		if (shift > 0)
		{
			// For each row,
			for (int row_i = 0; row_i < n; ++row_i, row += pitch)
			{
				// For each word in the row,
				int add_pitch = (n + 63) / 64;
				u64 prev = 0;
				for (int ii = 0; ii < add_pitch - 1; ++ii)
				{
					// Generate next word
					u32 rv1 = prng.Next();
					u32 rv2 = prng.Next();
					u64 word = ((u64)rv2 << 32) | rv1;

					// Add it in
					row[ii] ^= (prev >> (64 - shift)) | (word << shift);

					prev = word;
				}

				// Generate next word
				u32 rv1 = prng.Next();
				u32 rv2 = prng.Next();
				u64 word = ((u64)rv2 << 32) | rv1;

				// Add last word if needed
				int last_bit = (shift + n + 63) / 64;
				if (last_bit > add_pitch)
				{
					// Add it in
					row[add_pitch - 1] ^= (prev >> (64 - shift)) | (word << shift);
					prev = word;

					// Add it in preserving trailing bits
					int write_count = (shift + n) & 63;
					word = prev >> (64 - shift);
					row[add_pitch] ^= (write_count == 0) ? word : (word & (((u64)1 << write_count) - 1));
				}
				else
				{
					// Add it in preserving trailing bits
					int write_count = (shift + n) & 63;
					word = (prev >> (64 - shift)) | (word << shift);
					row[add_pitch - 1] ^= (write_count == 0) ? word : (word & (((u64)1 << write_count) - 1));
				}
			}
		}
		else // Rare aligned case:
		{
			// For each row,
			for (int row_i = 0; row_i < n; ++row_i, row += pitch)
			{
				// For each word in the row,
				int add_pitch = (n + 63) / 64;
				for (int ii = 0; ii < add_pitch - 1; ++ii)
				{
					// Generate next word
					u32 rv1 = prng.Next();
					u32 rv2 = prng.Next();
					u64 word = ((u64)rv2 << 32) | rv1;

					// Add it in
					row[ii] ^= word;
				}

				// Generate next word
				u32 rv1 = prng.Next();
				u32 rv2 = prng.Next();
				u64 word = ((u64)rv2 << 32) | rv1;

				// Add it in preserving trailing bits
				int write_count = n & 63;
				row[add_pitch - 1] ^= (write_count == 0) ? word : (word & (((u64)1 << write_count) - 1));
			}
		}
	}
	else
	{
		// For each row,
		u64 * CAT_RESTRICT row = matrix;
		for (int ii = 0; ii < n; ++ii, row += pitch)
		{
			// Flip diagonal bit
			int column_i = offset + ii;
			row[column_i >> 6] ^= (u64)1 << (column_i & 63);
		}
	}

	return true;
}


//// Utility: Deck Shuffling function

/*
	Given a PRNG, generate a deck of cards in a random order.
	The deck will contain elements with values between 0 and count - 1.
*/

static void ShuffleDeck16(Abyssinian &prng, u16 * CAT_RESTRICT deck, u32 count)
{
	deck[0] = 0;

	// If we can unroll 4 times,
	if (count <= 256)
	{
		for (u32 ii = 1;;)
		{
			u32 jj, rv = prng.Next();

			// 8-bit unroll
			switch (count - ii)
			{
			default:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 24) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				break;

			case 3:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 2:
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 1:
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
			case 0:
				return;
			}
		}
	}
	else
	{
		// For each deck entry,
		for (u32 ii = 1;;)
		{
			u32 jj, rv = prng.Next();

			// 16-bit unroll
			switch (count - ii)
			{
			default:
				jj = (u16)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u16)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				break;

			case 1:
				jj = (u16)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
			case 0:
				return;
			}
		}
	}
}


//// Utility: GF(256) Multiply and Divide functions

/*
	Branchless multiply and divide construction from
	"Fast Software Implementations of Finite Field Operations (Extended Abstract)"
	by Cheng Huang and Lihao Xu

	Small corrections made to paper (Q = 255):
		+ The EXP_TABLE needs to have 512*2+1 elements to handle 0*0 = 0 case.
		+ Element 255*2 should be set to 1.

	After these corrections it works properly and reduces the execution time
	to 58% of the usual version that uses branches to handle zero input.

	These tables were generated using polynomial 0x15F.  Maybe it's voodoo but
	random GF(256) matrices with this polynomial tended to be more invertible.
	There are 16 generator polynomials for GF(256), and 0x1F5 was a close second
	in terms of rate of invertibility.

	INV_TABLE[x] was also generated to accelerate GF256Divide(1, x).
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
	115, 154, 201, 35, 148, 139, 71, 144, 253, 242, 46, 82, 231, 105, 54, 89
};

static const u8 EXP_TABLE[512*2+1] = {
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

static const u8 INV_TABLE[256] = {
	0, 1, 175, 202, 248, 70, 101, 114, 124, 46, 35, 77, 157, 54, 57, 247,
	62, 152, 23, 136, 190, 244, 137, 18, 225, 147, 27, 26, 179, 59, 212, 32,
	31, 213, 76, 10, 164, 182, 68, 220, 95, 144, 122, 113, 235, 195, 9, 125,
	223, 253, 230, 189, 162, 120, 13, 156, 246, 14, 178, 29, 106, 84, 16, 153,
	160, 119, 197, 198, 38, 221, 5, 249, 82, 159, 91, 207, 34, 11, 110, 166,
	128, 104, 72, 158, 61, 107, 151, 201, 218, 116, 206, 74, 171, 155, 145, 40,
	192, 139, 209, 134, 115, 6, 241, 180, 81, 129, 60, 85, 169, 176, 78, 167,
	123, 43, 7, 100, 89, 219, 161, 65, 53, 163, 42, 112, 8, 47, 227, 187,
	80, 105, 148, 232, 205, 214, 99, 208, 19, 22, 193, 97, 173, 229, 211, 238,
	41, 94, 224, 25, 130, 233, 200, 86, 17, 63, 170, 93, 55, 12, 83, 73,
	64, 118, 52, 121, 36, 183, 79, 111, 177, 108, 154, 92, 228, 140, 203, 2,
	109, 168, 58, 28, 103, 240, 37, 165, 250, 217, 226, 127, 231, 51, 20, 245,
	96, 138, 234, 45, 199, 66, 67, 196, 150, 87, 3, 174, 215, 132, 90, 75,
	135, 98, 239, 142, 30, 33, 133, 204, 251, 185, 88, 117, 39, 69, 252, 48,
	146, 24, 186, 126, 172, 141, 50, 188, 131, 149, 194, 44, 255, 243, 143, 210,
	181, 102, 254, 237, 21, 191, 56, 15, 4, 71, 184, 216, 222, 49, 242, 236
};

static CAT_INLINE u8 GF256Multiply(u8 x, u8 y)
{
	return EXP_TABLE[LOG_TABLE[x] + LOG_TABLE[y]];
}

static CAT_INLINE u8 GF256Divide(u8 x, u8 y)
{
	// Precondition: y != 0
	return EXP_TABLE[LOG_TABLE[x] + 255 - LOG_TABLE[y]];
}

// Performs "dest[] += src[] * x" operation in GF(256)
static void GF256MemMulAdd(void * CAT_RESTRICT vdest, u8 x, const void * CAT_RESTRICT vsrc, int bytes)
{
	u8 * CAT_RESTRICT dest = reinterpret_cast<u8*>( vdest );
	const u8 * CAT_RESTRICT src = reinterpret_cast<const u8*>( vsrc );
	int log_x = LOG_TABLE[x];

	// For each block of 8 bytes,
	while (bytes >= 8)
	{
		/*
			For smaller messages, this function takes
			up 50% of execution time.  Unfortunately I
			do not see a way to make it run any faster.
		*/
		dest[0] ^= EXP_TABLE[LOG_TABLE[src[0]] + log_x];
		dest[1] ^= EXP_TABLE[LOG_TABLE[src[1]] + log_x];
		dest[2] ^= EXP_TABLE[LOG_TABLE[src[2]] + log_x];
		dest[3] ^= EXP_TABLE[LOG_TABLE[src[3]] + log_x];
		dest[4] ^= EXP_TABLE[LOG_TABLE[src[4]] + log_x];
		dest[5] ^= EXP_TABLE[LOG_TABLE[src[5]] + log_x];
		dest[6] ^= EXP_TABLE[LOG_TABLE[src[6]] + log_x];
		dest[7] ^= EXP_TABLE[LOG_TABLE[src[7]] + log_x];

		src += 8;
		dest += 8;
		bytes -= 8;
	}

	// For each byte,
	while (bytes-- > 0)
	{
		// Multiply source byte by x and add it to destination byte
		*dest++ ^= EXP_TABLE[LOG_TABLE[*src++] + log_x];
	}
}

// Performs "dest[] /= x" operation in GF(256)
static void GF256MemDivide(void * CAT_RESTRICT vdest, u8 x, int bytes)
{
	u8 * CAT_RESTRICT dest = reinterpret_cast<u8*>( vdest );
	int log_x = 255 - LOG_TABLE[x];

	// For each block of 8 bytes,
	while (bytes >= 8)
	{
		dest[0] = EXP_TABLE[LOG_TABLE[dest[0]] + log_x];
		dest[1] = EXP_TABLE[LOG_TABLE[dest[1]] + log_x];
		dest[2] = EXP_TABLE[LOG_TABLE[dest[2]] + log_x];
		dest[3] = EXP_TABLE[LOG_TABLE[dest[3]] + log_x];
		dest[4] = EXP_TABLE[LOG_TABLE[dest[4]] + log_x];
		dest[5] = EXP_TABLE[LOG_TABLE[dest[5]] + log_x];
		dest[6] = EXP_TABLE[LOG_TABLE[dest[6]] + log_x];
		dest[7] = EXP_TABLE[LOG_TABLE[dest[7]] + log_x];

		dest += 8;
		bytes -= 8;
	}

	// For each byte,
	while (bytes-- > 0)
	{
		// Multiply source byte by x and add it to destination byte
		*dest = EXP_TABLE[LOG_TABLE[*dest] + log_x];
		++dest;
	}
}


//// Utility: Column Iterator function

/*
	This implements a very light PRNG (Weyl function) to quickly generate
	a set of random-looking columns without replacement.

	This is Stewart Platt's excellent loop-less iterator optimization.
	His common cases all require no additional modulus operation, which
	makes it faster than the rare case that I designed.
*/

static CAT_INLINE void IterateNextColumn(u16 &x, u16 b, u16 p, u16 a)
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


//// Utility: Peeling Row Weight Generator function

/*
	Ideal Soliton weight distribution from
	"LT Codes" (2002)
	by Michael Luby

		The weight distribution selected for use in this codec is
	the Ideal Soliton distribution.  The PMF for weights 2 and higher
	is 1 / (k*(k - 1)).  Accumulating these yields the WEIGHT_DIST
	table below up to weight 64.  I stuck ~0 at the end of the
	table to make sure the while loop in the function terminates.

		To produce a good code, the probability of weight-1 should be
	added into each element of the WEIGHT_DIST table.  Instead, I add
	it programmatically in the function to allow it to be easily tuned.

		I played around with where to truncate this table, and found
	that for higher block counts, the number of deferred rows after
	greedy peeling is much lower for 64 weights than 32.  And after
	tuning the codec for weight 64, the performance was slightly
	better than with 32.

		I also tried different probabilities for weight-1 rows, and
	settled on 1/128 as having the best performance in a few select
	tests.  Setting it too high or too low (or even to zero) tends
	to reduce the performance of the codec.

		However, once N gets much much larger, it is actually very
	beneficial to switch over to weight-2 as a minimum.
*/

static const u32 MAX_WEIGHT_1 = 4096; // Stop using weight-1 after this

static const u32 WEIGHT_DIST[] = {
	0x00000000, 0x80000000, 0xaaaaaaaa, 0xc0000000, 0xcccccccc, 0xd5555555, 0xdb6db6db, 0xe0000000,
	0xe38e38e3, 0xe6666666, 0xe8ba2e8b, 0xeaaaaaaa, 0xec4ec4ec, 0xedb6db6d, 0xeeeeeeee, 0xefffffff,
	0xf0f0f0f0, 0xf1c71c71, 0xf286bca1, 0xf3333333, 0xf3cf3cf3, 0xf45d1745, 0xf4de9bd3, 0xf5555555,
	0xf5c28f5c, 0xf6276276, 0xf684bda1, 0xf6db6db6, 0xf72c234f, 0xf7777777, 0xf7bdef7b, 0xf7ffffff,
	0xf83e0f83, 0xf8787878, 0xf8af8af8, 0xf8e38e38, 0xf914c1ba, 0xf9435e50, 0xf96f96f9, 0xf9999999,
	0xf9c18f9c, 0xf9e79e79, 0xfa0be82f, 0xfa2e8ba2, 0xfa4fa4fa, 0xfa6f4de9, 0xfa8d9df5, 0xfaaaaaaa,
	0xfac687d6, 0xfae147ae, 0xfafafafa, 0xfb13b13b, 0xfb2b78c1, 0xfb425ed0, 0xfb586fb5, 0xfb6db6db,
	0xfb823ee0, 0xfb9611a7, 0xfba93868, 0xfbbbbbbb, 0xfbcda3ac, 0xfbdef7bd, 0xfbefbefb, 0xffffffff
};

static u16 GeneratePeelRowWeight(u32 rv, u16 peel_column_count)
{
	// Unroll first 3 for speed (common case):

	// If peel column count is small,
	if (peel_column_count <= MAX_WEIGHT_1)
	{
		// Select probability of weight-1 rows here:
		static const u32 P1 = (u32)((1./128) * 0xffffffff);
		if (rv < P1) return 1;

		// Rescale to match table values
		rv -= P1;
	}

	static const u32 P2 = WEIGHT_DIST[1];
	if (rv <= P2) return 2;

	static const u32 P3 = WEIGHT_DIST[2];
	if (rv <= P3) return 3;

	// Find first table entry containing a number smaller than or equal to rv
	u16 weight = 3;
	while (rv > WEIGHT_DIST[weight++]);
	return weight;
}


//// Utility: Peel Matrix Row Generator function

static void GeneratePeelRow(u32 id, u32 p_seed, u16 peel_column_count, u16 mix_column_count,
	u16 & CAT_RESTRICT peel_weight, u16 & CAT_RESTRICT peel_a, u16 & CAT_RESTRICT peel_x0,
	u16 & CAT_RESTRICT mix_a, u16 & CAT_RESTRICT mix_x0)
{
	// Initialize PRNG
	Abyssinian prng;
	prng.Initialize(id, p_seed);

	// Generate peeling matrix row weight
	u16 weight = GeneratePeelRowWeight(prng.Next(), peel_column_count);
	u16 max_weight = peel_column_count / 2; // Do not set more than N/2 at a time
	peel_weight = (weight > max_weight) ? max_weight : weight;

	// Generate peeling matrix column selection parameters for row
	u32 rv = prng.Next();
	peel_a = ((u16)rv % (peel_column_count - 1)) + 1;
	peel_x0 = (u16)(rv >> 16) % peel_column_count;

	// Generate mixing matrix column selection parameters
	rv = prng.Next();
	mix_a = ((u16)rv % (mix_column_count - 1)) + 1;
	mix_x0 = (u16)(rv >> 16) % mix_column_count;
}


//// Data Structures

#pragma pack(push)
#pragma pack(1)
struct Codec::PeelRow
{
	u16 next;					// Linkage in row list
	u32 id;						// Identifier for this row

	// Peeling matrix: Column generator
	u16 peel_weight, peel_a, peel_x0;

	// Mixing matrix: Column generator
	u16 mix_a, mix_x0;

	// Peeling state
	u16 unmarked_count;			// Count of columns that have not been marked yet
	union
	{
		// During peeling:
		u16 unmarked[2];		// Final two unmarked column indices

		// After peeling:
		struct
		{
			u16 peel_column;	// Peeling column that is solved by this row
			u8 is_copied;		// Row value is copied yet?
		};
	};
};
#pragma pack(pop)

// Marks for PeelColumn
enum MarkTypes
{
	MARK_TODO,	// Unmarked
	MARK_PEEL,	// Was solved by peeling
	MARK_DEFER	// Deferred to Gaussian elimination
};

#pragma pack(push)
#pragma pack(1)
struct Codec::PeelColumn
{
	u16 next;			// Linkage in column list

	union
	{
		u16 w2_refs;	// Number of weight-2 rows containing this column
		u16 peel_row;	// Row that solves the column
		u16 ge_column;	// Column that a deferred column is mapped to
	};

	u8 mark;			// One of the MarkTypes enumeration
};
#pragma pack(pop)

#pragma pack(push)
#pragma pack(1)
struct Codec::PeelRefs
{
	u16 row_count;		// Number of rows containing this column
	u16 rows[CAT_REF_LIST_MAX];
};
#pragma pack(pop)


//// (1) Peeling:

/*
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

/*
	OpportunisticPeeling

		This function accepts a new row from the codec input and
	immediately attempts to solve a column opportunistically using
	the graph-based peeling decoding process.

		The row value is assumed to already be present in the
	_input_blocks data, but this function does take care of
	initializing everything else for a new row, including the
	row ID number and the peeling column generator parameters.
*/

bool Codec::OpportunisticPeeling(u32 row_i, u32 id)
{
	PeelRow *row = &_peel_rows[row_i];

	row->id = id;
	GeneratePeelRow(id, _p_seed, _block_count, _mix_count,
		row->peel_weight, row->peel_a, row->peel_x0, row->mix_a, row->mix_x0);

	CAT_IF_DUMP(cout << "Row " << id << " in slot " << row_i << " of weight " << row->peel_weight << " [a=" << row->peel_a << "] : ";)

	// Iterate columns in peeling matrix
	u16 weight = row->peel_weight;
	u16 column_i = row->peel_x0;
	u16 a = row->peel_a;
	u16 unmarked_count = 0;
	u16 unmarked[2];
	for (;;)
	{
		CAT_IF_DUMP(cout << column_i << " ";)

		PeelRefs *refs = &_peel_col_refs[column_i];

		// Add row reference to column
		if (refs->row_count >= CAT_REF_LIST_MAX)
		{
			CAT_IF_DUMP(cout << "OpportunisticPeeling: Failure!  Ran out of space for row references.  CAT_REF_LIST_MAX must be increased!" << endl;)
			FixPeelFailure(row, column_i);
			return false;
		}
		refs->rows[refs->row_count++] = row_i;

		// If column is unmarked,
		if (_peel_cols[column_i].mark == MARK_TODO)
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

	return true;
}

/*
	FixPeelFailure

		This function unreferences previous columns for a row where,
	in one of the columns, the reference list overflowed.  This avoids
	potential data corruption or more severe problems in case of
	unusually distributed peeling matrices.
*/

void Codec::FixPeelFailure(PeelRow * CAT_RESTRICT row, u16 fail_column_i)
{
	CAT_IF_DUMP(cout << "!!Fixing Peel Failure!! Unreferencing columns, ending at " << fail_column_i << " :";)

	// Iterate columns in peeling matrix
	//u16 weight = row->peel_weight;
	u16 column_i = row->peel_x0;
	u16 a = row->peel_a;
	while (column_i != fail_column_i)
	{
		CAT_IF_DUMP(cout << " " << column_i;)

		PeelRefs * CAT_RESTRICT refs = &_peel_col_refs[column_i];

		// Subtract off row count - Invalidates row number that was written earlier
		refs->row_count--;

		// NOTE: Does not need to validate weight here since fail_column_i is guaranteed to come around

		IterateNextColumn(column_i, _block_count, _block_next_prime, a);
	}
	CAT_IF_DUMP(cout << endl;)
}

/*
	PeelAvalanche

		This function is called after a column is solved by a row during
	peeling.  It attempts to find other rows that reference this column
	and resumes opportunistic peeling in an avalanche of solutions.

		OpportunisticPeeling() and PeelAvalanche() are split up into two
	functions because I found that the PeelAvalanche() function can be
	reused later during GreedyPeeling().
*/

void Codec::PeelAvalanche(u16 column_i)
{
	// Walk list of peeled rows referenced by this newly solved column
	PeelRefs * CAT_RESTRICT refs = &_peel_col_refs[column_i];
	u16 ref_row_count = refs->row_count;
	u16 * CAT_RESTRICT ref_rows = refs->rows;
	while (ref_row_count--)
	{
		// Update unmarked row count for this referenced row
		u16 ref_row_i = *ref_rows++;
		PeelRow * CAT_RESTRICT ref_row = &_peel_rows[ref_row_i];
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
				PeelColumn * CAT_RESTRICT ref_col = &_peel_cols[ref_column_i];

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

/*
	Peel

		This function is called exclusively by OpportunisticPeeling()
	to take care of marking columns solved when a row is able to solve
	a column during the peeling process.
*/

void Codec::Peel(u16 row_i, PeelRow * CAT_RESTRICT row, u16 column_i)
{
	CAT_IF_DUMP(cout << "Peel: Solved column " << column_i << " with row " << row_i << endl;)

	PeelColumn * CAT_RESTRICT column = &_peel_cols[column_i];

	// Mark this column as solved
	column->mark = MARK_PEEL;

	// Remember which column it solves
	row->peel_column = column_i;

	// Link to back of the peeled list
	if (_peel_tail_rows)
		_peel_tail_rows->next = row_i;
	else
		_peel_head_rows = row_i;
	row->next = LIST_TERM;
	_peel_tail_rows = row;

	// Indicate that this row hasn't been copied yet
	row->is_copied = 0;

	// Attempt to avalanche and solve other columns
	PeelAvalanche(column_i);

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
	columns that tied, the one that affects the most rows is selected.  This is
	a better way to choose which columns to defer than just selecting the one
	referenced by the most rows.

		In practice with a well designed peeling matrix, about sqrt(N) + N/150
	columns must be deferred to Gaussian elimination using this greedy approach.
*/

void Codec::GreedyPeeling()
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
		PeelColumn *column = _peel_cols;
		for (u16 column_i = 0; column_i < _block_count; ++column_i, ++column)
		{
			// If column is not marked yet,
			if (column->mark == MARK_TODO)
			{
				// And if it may have the most weight-2 references
				u16 w2_refs = column->w2_refs;
				if (w2_refs >= best_w2_refs)
				{
					// Or if it has the largest row references overall,
					u16 row_count = _peel_col_refs[column_i].row_count;
					if (w2_refs > best_w2_refs || row_count >= best_row_count)
					{
						// Use that one
						best_column_i = column_i;
						best_w2_refs = w2_refs;
						best_row_count = row_count;
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
		PeelAvalanche(best_column_i);
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
	If using a naive approach, this is by far the most complex step of
	the matrix solver.  The approach outlined here makes it much
	easier.  At this point the check matrix has been re-organized
	into peeled and deferred rows and columns:

		+-----------------------+
		| P P P P P | D D | M M |
		+-----------------------+
					X
		+-----------+-----+-----+
		| 5 2 6 1 3 | 4 0 | 7 8 |
		+-----------+-----+-----+---+   +---+
		| 1         | 0 0 | 1 0 | 0 |   | P |
		| 0 1       | 0 0 | 0 1 | 5 |   | P |
		| 1 0 1     | 1 0 | 1 0 | 3 |   | P |
		| 0 1 0 1   | 0 0 | 0 1 | 4 |   | P |
		| 0 1 0 1 1 | 1 1 | 1 0 | 1 |   | P |
		+-----------+-----+-----+---| = |---|
		| 1 1 0 1 1 | 1 1 | 1 0 | 7 |   | 0 |
		| 1 0 1 1 0 | 1 0 | 0 1 | 8 |   | 0 |
		+-----------+-----+-----+---+   +---+
		| 0 1 1 0 0 | 1 1 | 0 1 | 2 |   | D |
		| 0 1 0 1 0 | 0 1 | 1 0 | 6 |   | D |
		+-----------+-----+-----+---|   |---|
		      ^          ^     ^- Mixing columns
		      |          \------- Deferred columns
		      \------------------ Peeled columns intersections with deferred rows

	P = Peeled rows/columns (re-ordered)
	D = Deferred rows/columns (order of deferment)
	M = Mixing columns always deferred for GE
	0 = Dense rows always deferred for GE; they sum to 0

		Since the re-ordered matrix above is in lower triangular form,
	and since the weight of each row is limited to a constant, the cost
	of diagonalizing the peeled matrix is O(n).  Diagonalizing the peeled
	matrix will cause the mix and deferred columns to add up and become
	dense.  These columns are stored in the same order as in the GE matrix,
	in a long vertical matrix with N rows, called the Compression matrix.

		After diagonalizing the peeling matrix, the peeling matrix will be
	the identity matrix.  Peeled column output blocks can be used to store
	the temporary block values generated by this process.  These temporary
	blocks will be overwritten and lost later during Substitution.

		To finish compressing the matrix into a form for Gaussian elimination,
	all of the peeled columns of the deferred/dense rows must be zeroed.
	Wherever a column is set to 1 in a deferred/dense row, the Compression
	matrix row that solves that peeled column is added to the deferred row.
	Essentially the peeled column intersection with the deferred rows are
	multiplied by the peeled matrix to produce initial row values and matrix
	rows for the GE matrix in the lower right.

		This process does not actually produce any final column values because
	at this point it is not certain where those values will end up.  Instead,
	a record of what operations were performed needs to be stored and followed
	later after the destination columns are determined by Gaussian elimination.
*/

/*
	SetDeferredColumns

		This function initializes some mappings between GE columns and
	column values, and it sets bits in the Compression matrix in rows
	that reference the deferred columns.  These bits will get mixed
	throughout the Compression matrix and will make it very dense.

	For each deferred column,
		Set bit for each row affected by this column.
		Map GE column to this column.
		Map this column to the GE column.

	For each mixing column,
		Map GE column to this column.
*/

void Codec::SetDeferredColumns()
{
	CAT_IF_DUMP(cout << endl << "---- SetDeferredColumns ----" << endl << endl;)

	// For each deferred column,
	PeelColumn * CAT_RESTRICT column;
	for (u16 ge_column_i = 0, defer_i = _defer_head_columns; defer_i != LIST_TERM; defer_i = column->next, ++ge_column_i)
	{
		column = &_peel_cols[defer_i];

		CAT_IF_DUMP(cout << "GE column " << ge_column_i << " mapped to matrix column " << defer_i << " :";)

		// Set bit for each row affected by this deferred column
		u64 *matrix_row_offset = _compress_matrix + (ge_column_i >> 6);
		u64 ge_mask = (u64)1 << (ge_column_i & 63);
		PeelRefs * CAT_RESTRICT refs = &_peel_col_refs[defer_i];
		u16 count = refs->row_count;
		u16 *ref_row = refs->rows;
		while (count--)
		{
			u16 row_i = *ref_row++;

			CAT_IF_DUMP(cout << " " << row_i;)

			matrix_row_offset[_ge_pitch * row_i] |= ge_mask;
		}

		CAT_IF_DUMP(cout << endl;)

		// Set column map for this GE column
		_ge_col_map[ge_column_i] = defer_i;

		// Set reverse mapping also
		column->ge_column = ge_column_i;
	}

	// Set column map for each mix column
	for (u16 added_i = 0; added_i < _mix_count; ++added_i)
	{
		u16 ge_column_i = _defer_count + added_i;
		u16 column_i = _block_count + added_i;

		CAT_IF_DUMP(cout << "GE column(mix) " << ge_column_i << " mapped to matrix column " << column_i << endl;)

		_ge_col_map[ge_column_i] = column_i;
	}
}

/*
	SetMixingColumnsForDeferredRows

		This function generates the mixing column bits
	for each deferred row in the GE matrix.  It also
	marks the row's peel_column with LIST_TERM so that
	later it will be easy to check if it was deferred.
*/

void Codec::SetMixingColumnsForDeferredRows()
{
	CAT_IF_DUMP(cout << endl << "---- SetMixingColumnsForDeferredRows ----" << endl << endl;)

	// For each deferred row,
	PeelRow * CAT_RESTRICT row;
	for (u16 defer_row_i = _defer_head_rows; defer_row_i != LIST_TERM; defer_row_i = row->next)
	{
		row = &_peel_rows[defer_row_i];

		CAT_IF_DUMP(cout << "Deferred row " << defer_row_i << " set mix columns :";)

		// Mark it as deferred for the following loop
		row->peel_column = LIST_TERM;

		// Set up mixing column generator
		u64 *ge_row = _compress_matrix + _ge_pitch * defer_row_i;
		u16 a = row->mix_a;
		u16 x = row->mix_x0;

		// Generate mixing column 1
		u16 ge_column_i = _defer_count + x;
		ge_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
		CAT_IF_DUMP(cout << " " << ge_column_i;)
		IterateNextColumn(x, _mix_count, _mix_next_prime, a);

		// Generate mixing column 2
		ge_column_i = _defer_count + x;
		ge_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
		CAT_IF_DUMP(cout << " " << ge_column_i;)
		IterateNextColumn(x, _mix_count, _mix_next_prime, a);

		// Generate mixing column 3
		ge_column_i = _defer_count + x;
		ge_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
		CAT_IF_DUMP(cout << " " << ge_column_i;)

		CAT_IF_DUMP(cout << endl;)
	}
}

/*
	PeelDiagonal

		This function diagonalizes the peeled rows and columns of the
	check matrix.  The result is that the peeled submatrix is the
	identity matrix, and that the other columns of the peeled rows are
	very dense and have temporary block values assigned.  These dense
	columns are used to efficiently zero out the peeled columns of the
	other rows.

		This function is one of the most expensive in the whole codec,
	because its memory access patterns are not cache-friendly.

	For each peeled row in forward solution order,
		Set mixing column bits for the row in the Compression matrix.
		Generate row block value.
		For each row that references this row in the peeling matrix,
			Add Compression matrix row to referencing row.
			If row is peeled,
				Add row block value.
*/

void Codec::PeelDiagonal()
{
	CAT_IF_DUMP(cout << endl << "---- PeelDiagonal ----" << endl << endl;)

	/*
		This function optimizes the block value generation by combining the first
		memcpy and memxor operations together into a three-way memxor if possible,
		using the is_copied row member.
	*/

	CAT_IF_ROWOP(int rowops = 0;)

	// For each peeled row in forward solution order,
	PeelRow * CAT_RESTRICT row;
	for (u16 peel_row_i = _peel_head_rows; peel_row_i != LIST_TERM; peel_row_i = row->next)
	{
		row = &_peel_rows[peel_row_i];

		// Lookup peeling results
		u16 peel_column_i = row->peel_column;
		u64 *ge_row = _compress_matrix + _ge_pitch * peel_row_i;

		CAT_IF_DUMP(cout << "Peeled row " << peel_row_i << " for peeled column " << peel_column_i << " :";)

		// Set up mixing column generator
		u16 a = row->mix_a;
		u16 x = row->mix_x0;

		// Generate mixing column 1
		u16 ge_column_i = _defer_count + x;
		ge_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
		CAT_IF_DUMP(cout << " " << ge_column_i;)
		IterateNextColumn(x, _mix_count, _mix_next_prime, a);

		// Generate mixing column 2
		ge_column_i = _defer_count + x;
		ge_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
		CAT_IF_DUMP(cout << " " << ge_column_i;)
		IterateNextColumn(x, _mix_count, _mix_next_prime, a);

		// Generate mixing column 3
		ge_column_i = _defer_count + x;
		ge_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
		CAT_IF_DUMP(cout << " " << ge_column_i << endl;)

		// Lookup output block
		u8 * CAT_RESTRICT temp_block_src = _recovery_blocks + _block_bytes * peel_column_i;

		// If row has not been copied yet,
		if (!row->is_copied)
		{
			// Copy it directly to the output block
			const u8 * CAT_RESTRICT block_src = _input_blocks + _block_bytes * peel_row_i;
			if (peel_row_i != _block_count - 1)
				memcpy(temp_block_src, block_src, _block_bytes);
			else
			{
				memcpy(temp_block_src, block_src, _input_final_bytes);
				memset(temp_block_src + _input_final_bytes, 0, _block_bytes - _input_final_bytes);
			}
			CAT_IF_ROWOP(++rowops;)

			CAT_IF_DUMP(cout << "-- Copied from " << peel_row_i << " because has not been copied yet.  Output block = " << (int)temp_block_src[0] << endl;)

			// NOTE: Do not need to set is_copied here because no further rows reference this one
		}

		CAT_IF_DUMP(cout << "++ Adding to referencing rows:";)

		// For each row that references this one,
		PeelRefs * CAT_RESTRICT refs = &_peel_col_refs[peel_column_i];
		u16 count = refs->row_count;
		u16 * CAT_RESTRICT ref_row = refs->rows;
		while (count--)
		{
			u16 ref_row_i = *ref_row++;

			// Skip this row
			if (ref_row_i == peel_row_i) continue;

			CAT_IF_DUMP(cout << " " << ref_row_i;)

			// Add GE row to referencing GE row
			u64 * CAT_RESTRICT ge_ref_row = _compress_matrix + _ge_pitch * ref_row_i;
			for (int ii = 0; ii < _ge_pitch; ++ii) ge_ref_row[ii] ^= ge_row[ii];

			// If row is peeled,
			PeelRow * CAT_RESTRICT ref_row = &_peel_rows[ref_row_i];
			u16 ref_column_i = ref_row->peel_column;
			if (ref_column_i != LIST_TERM)
			{
				// Generate temporary row block value:
				u8 * CAT_RESTRICT temp_block_dest = _recovery_blocks + _block_bytes * ref_column_i;

				// If referencing row is already copied to the recovery blocks,
				if (ref_row->is_copied)
				{
					// Add this row block value to it
					memxor(temp_block_dest, temp_block_src, _block_bytes);
				}
				else
				{
					// Add this row block value with message block to it (optimization)
					const u8 * CAT_RESTRICT block_src = _input_blocks + _block_bytes * ref_row_i;
					if (ref_row_i != _block_count - 1)
						memxor_set(temp_block_dest, temp_block_src, block_src, _block_bytes);
					else
					{
						memxor_set(temp_block_dest, temp_block_src, block_src, _input_final_bytes);
						memcpy(temp_block_dest + _input_final_bytes, temp_block_src + _input_final_bytes, _block_bytes - _input_final_bytes);
					}

					ref_row->is_copied = 1;
				}
				CAT_IF_ROWOP(++rowops;)
			} // end if referencing row is peeled
		} // next referencing row

		CAT_IF_DUMP(cout << endl;)
	} // next peeled row

	CAT_IF_ROWOP(cout << "PeelDiagonal used " << rowops << " row ops = " << rowops / (double)_block_count << "*N" << endl;)
}

/*
	CopyDeferredRows

		This function copies deferred rows from the Compression matrix
	into their final location in the GE matrix.  It also maps the GE rows
	to the deferred rows.
*/

void Codec::CopyDeferredRows()
{
	CAT_IF_DUMP(cout << endl << "---- CopyDeferredRows ----" << endl << endl;)

	// For each deferred row,
	u64 * CAT_RESTRICT ge_row = _ge_matrix + _ge_pitch * _dense_count;
	for (u16 ge_row_i = _dense_count, defer_row_i = _defer_head_rows; defer_row_i != LIST_TERM;
		defer_row_i = _peel_rows[defer_row_i].next, ge_row += _ge_pitch, ++ge_row_i)
	{
		CAT_IF_DUMP(cout << "Peeled row " << defer_row_i << " for GE row " << ge_row_i << endl;)

		// Copy compress row to GE row
		u64 * CAT_RESTRICT compress_row = _compress_matrix + _ge_pitch * defer_row_i;
		memcpy(ge_row, compress_row, _ge_pitch * sizeof(u64));

		// Set row map for this deferred row
		_ge_row_map[ge_row_i] = defer_row_i;
	}
}

/*
	Important Optimization: Shuffle-2 Codes

		After compression, the GE matrix is a square matrix that
	looks like this:

		+-----------------+----------------+--------------+
		|  Dense Deferred | Dense Deferred | Dense Mixing | <- About half
		|                 |                |              |
		+-----------------+----------------+--------------+
		|                 | Dense Deferred | Dense Mixing | <- About half
		|        0        |                |              |
		|                 |----------------+--------------+
		|                 | Sparse Deferred| Dense Mixing | <- Last few rows
		+-----------------+----------------+--------------+
		         ^                  ^---- Middle third of the columns
				 \------ Left third of the columns

		The dense rows are generated so that they can quickly be
	eliminated with as few row operations as possible.
	This elimination can be visualized as a matrix-matrix multiplication
	between the peeling submatrix and the deferred/dense submatrix
	intersection with the peeled columns.  Using Shuffle-2 Codes, I
	have been able to achieve this matrix-matrix multiplication in
	just 2.5*N row operations, which is better than any other approach
	I have seen so far.

		I needed to find a way to generate a binary matrix that LOOKS
	random but actually only differs by ~2 bits per row.  I looked at
	using normal Gray codes or more Weyl generators but they both are
	restricted to a subset of the total possibilities.  Instead, the
	standard in-place shuffle algorithm is used to shuffle row and
	column orders to make it look random.  This new algorithm is able
	to generate nearly all possible combinations with approximately
	uniform likelihood, and the generated matrix can be addressed by
	a 32-bit seed, so it is easy to regenerate the same matrix again,
	or generate many new random matrices without reseeding.

		Shuffling generates the first row randomly, and each following
	row is XORed by two columns, one with a bit set and one without.
	The order of XOR pairs is decided by the initial shuffling.  The
	order of the generated rows is shuffled separately.

	Example output: A random 17x17 matrix

		10000001111010011
		00111110100101100
		11000001011010011
		00101100111110010
		00111100101101100
		00111100101110100
		11000011010010011
		01111110000101100
		01011111000101100
		00101000111010011
		00101100111110100
		11010011000001011
		00101000111110010
		10100000111010011
		11010111000001101
		11010111000101100
		11000011010001011

		This code I am calling a Perfect Shuffle-2 Code.

		These are "perfect" matrices in that they have the same
	Hamming weight in each row and each column.  The problem is
	that this type of matrix is NEVER invertible, so the perfect
	structure must be destroyed in order to get a good code for
	error correction.

		The resulting code is called a Shuffle-2 code.

		Here is the Shuffle-2 process:

	Split the dense submatrix of the check matrix into DxD squares.
	For each DxD square subsubmatrix,
		Shuffle the destination row order.
		Shuffle the bit flip order.
		Generate a random bit string with weight D/2 for the first output row.
		Reshuffle the bit flip order. <- Helps recovery properties a lot!
		Flip two bits for each row of the first half of the outputs.
		Reshuffle the bit flip order. <- Helps recovery properties a lot!
		Flip two bits for each row of the last half of the outputs.

		This effectively destroys the perfection of the code, and makes
	the square matrices invertible about as often as a random GF(2) code,
	so that using these easily generated square matrices does not hurt
	the error correction properties of the code.  Random GF(2) matrices
	are invertible about 30% of the time, and these are invertible about
	15% of the time when D Mod 4 = 2.  Other choices of D are not so good.

		A Shuffle-3 Code would reshuffle 3 times and flip 3 bits per row,
	and a Shuffle-4 Code would reshuffle 4 times and flip 4 bits per row.
	This is probably not the first time that someone has invented this.
	I believe that Moon Ho Lee has come up with something similar and more
	mathematically rigorous, though I believe he was using Shuffle-4 for
	quantum cryptography.  Shuffle-2 is much faster and works for this
	application because columns from different DxD matrices are randomly
	selected for use in the GE matrix.  And furthermore the rank needed
	from the selected columns is usually much less than the row count.

		MultiplyDenseRows() does not actually use memxor() to generate
	any row block values because it is not certain where the values
	will end up, yet.  So instead this multiplication is done again
	in the MultiplyDenseValues() function after Triangle() succeeds.
*/

void Codec::MultiplyDenseRows()
{
	CAT_IF_DUMP(cout << endl << "---- MultiplyDenseRows ----" << endl << endl;)

	// Initialize PRNG
	Abyssinian prng;
	prng.Initialize(_d_seed);

	// For each block of columns,
	PeelColumn * CAT_RESTRICT column = _peel_cols;
	u64 * CAT_RESTRICT temp_row = _ge_matrix + _ge_pitch * (_dense_count + _defer_count);
	const int dense_count = _dense_count;
	u16 rows[CAT_MAX_DENSE_ROWS], bits[CAT_MAX_DENSE_ROWS];
	for (u16 column_i = 0; column_i < _block_count; column_i += dense_count, column += dense_count)
	{
		CAT_IF_DUMP(cout << "Shuffled dense matrix starting at column " << column_i << ":" << endl;)

		// Handle final columns
		int max_x = dense_count;
		if (column_i + dense_count > _block_count)
			max_x = _block_count - column_i;

		// Shuffle row and bit order
		ShuffleDeck16(prng, rows, dense_count);
		ShuffleDeck16(prng, bits, dense_count);

		// Initialize counters
		const u16 set_count = (dense_count + 1) >> 1;
		const u16 * CAT_RESTRICT set_bits = bits;
		const u16 * CAT_RESTRICT clr_bits = set_bits + set_count;

		CAT_IF_DUMP( u64 disp_row[(CAT_MAX_DENSE_ROWS+63)/64]; CAT_OBJCLR(disp_row); )

		// Generate first row
		memset(temp_row, 0, _ge_pitch * sizeof(u64));
		for (int ii = 0; ii < set_count; ++ii)
		{
			// If bit is peeled,
			int bit_i = set_bits[ii];
			if (bit_i < max_x)
			{
				if (column[bit_i].mark == MARK_PEEL)
				{
					// Add temp row value
					u64 * CAT_RESTRICT ge_source_row = _compress_matrix + _ge_pitch * column[bit_i].peel_row;
					for (int jj = 0; jj < _ge_pitch; ++jj) temp_row[jj] ^= ge_source_row[jj];
				}
				else
				{
					// Set GE bit for deferred column
					u16 ge_column_i = column[bit_i].ge_column;
					temp_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
				}
			}
			CAT_IF_DUMP(disp_row[bit_i >> 6] ^= (u64)1 << (bit_i & 63);)
		} // next bit

		// Set up generator
		const u16 * CAT_RESTRICT row = rows;

		// Store first row
		CAT_IF_DUMP(for (int ii = 0; ii < dense_count; ++ii) cout << ((disp_row[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row << endl;)
		u64 * CAT_RESTRICT ge_dest_row = _ge_matrix + _ge_pitch * *row++;
		for (int jj = 0; jj < _ge_pitch; ++jj) ge_dest_row[jj] ^= temp_row[jj];

		// Reshuffle bit order: Shuffle-2 Code
		ShuffleDeck16(prng, bits, dense_count);

		// Generate first half of rows
		const int loop_count = (dense_count >> 1);
		for (int ii = 0; ii < loop_count; ++ii)
		{
			int bit0 = set_bits[ii], bit1 = clr_bits[ii];

			// Flip bit 1
			if (bit0 < max_x)
			{
				if (column[bit0].mark == MARK_PEEL)
				{
					// Add temp row value
					u64 * CAT_RESTRICT ge_source_row = _compress_matrix + _ge_pitch * column[bit0].peel_row;
					for (int jj = 0; jj < _ge_pitch; ++jj) temp_row[jj] ^= ge_source_row[jj];
				}
				else
				{
					// Set GE bit for deferred column
					u16 ge_column_i = column[bit0].ge_column;
					temp_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
				}
			}
			CAT_IF_DUMP(disp_row[bit0 >> 6] ^= (u64)1 << (bit0 & 63);)

			// Flip bit 2
			if (bit1 < max_x)
			{
				if (column[bit1].mark == MARK_PEEL)
				{
					// Add temp row value
					u64 * CAT_RESTRICT ge_source_row = _compress_matrix + _ge_pitch * column[bit1].peel_row;
					for (int jj = 0; jj < _ge_pitch; ++jj) temp_row[jj] ^= ge_source_row[jj];
				}
				else
				{
					// Set GE bit for deferred column
					u16 ge_column_i = column[bit1].ge_column;
					temp_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
				}
			}
			CAT_IF_DUMP(disp_row[bit1 >> 6] ^= (u64)1 << (bit1 & 63);)

			// Store in row
			CAT_IF_DUMP(for (int ii = 0; ii < dense_count; ++ii) cout << ((disp_row[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row << endl;)
			ge_dest_row = _ge_matrix + _ge_pitch * *row++;
			for (int jj = 0; jj < _ge_pitch; ++jj) ge_dest_row[jj] ^= temp_row[jj];
		} // next row

		// Reshuffle bit order: Shuffle-2 Code
		ShuffleDeck16(prng, bits, dense_count);

		// Generate second half of rows
		const int second_loop_count = loop_count - 1 + (dense_count & 1);
		for (int ii = 0; ii < second_loop_count; ++ii)
		{
			int bit0 = set_bits[ii], bit1 = clr_bits[ii];

			// Flip bit 1
			if (bit0 < max_x)
			{
				if (column[bit0].mark == MARK_PEEL)
				{
					// Add temp row value
					u64 * CAT_RESTRICT ge_source_row = _compress_matrix + _ge_pitch * column[bit0].peel_row;
					for (int jj = 0; jj < _ge_pitch; ++jj) temp_row[jj] ^= ge_source_row[jj];
				}
				else
				{
					// Set GE bit for deferred column
					u16 ge_column_i = column[bit0].ge_column;
					temp_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
				}
			}
			CAT_IF_DUMP(disp_row[bit0 >> 6] ^= (u64)1 << (bit0 & 63);)

			// Flip bit 2
			if (bit1 < max_x)
			{
				if (column[bit1].mark == MARK_PEEL)
				{
					// Add temp row value
					u64 * CAT_RESTRICT ge_source_row = _compress_matrix + _ge_pitch * column[bit1].peel_row;
					for (int jj = 0; jj < _ge_pitch; ++jj) temp_row[jj] ^= ge_source_row[jj];
				}
				else
				{
					// Set GE bit for deferred column
					u16 ge_column_i = column[bit1].ge_column;
					temp_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
				}
			}
			CAT_IF_DUMP(disp_row[bit1 >> 6] ^= (u64)1 << (bit1 & 63);)

			// Store in row
			CAT_IF_DUMP(for (int ii = 0; ii < dense_count; ++ii) cout << ((disp_row[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row << endl;)
			ge_dest_row = _ge_matrix + _ge_pitch * *row++;
			for (int jj = 0; jj < _ge_pitch; ++jj) ge_dest_row[jj] ^= temp_row[jj];
		} // next row

		CAT_IF_DUMP(cout << endl;)
	} // next column
}

/*
	Important Optimization: O(1) Heavy Row Structure

		The heavy rows are made up of bytes instead of bits.  Each byte
	represents a number in the Galois field GF(2^^8) defined by the
	generator polynomial 0x15F.

		The heavy rows are designed to make it easier to find pivots in
	just the last few columns of the GE matrix.  This design choice was
	made because it allows a constant-time algorithm to be employed that
	will reduce the fail rate from >70% to <3%.  It is true that with
	heavy loss rates, the earlier columns can be where the pivot is
	needed.  However in my estimation it would be better to increase the
	number of dense rows instead to handle this problem than to increase
	the number of heavy rows.  The result is that we can assume missing
	pivots occur near the end.

		The number of heavy rows required is at least 5.  This is because
	the heavy rows are used to fill in for missing pivots in the GE matrix
	and the miss rate is about 1/2 because it's random and binary.  The odds
	of the last 5 columns all being zero in the binary rows is 1/32.
	And the odds of a random GF(256) matrix not being invertible is also
	around 1/32, therefore it needs at least 5 heavy rows.  With less than 5
	rows, the binary matrix fail rate would dominate the overall rate of
	invertibility.  After 5 heavy rows, less likely problems can be overcome,
	so 6 heavy rows were chosen for the baseline version.

		An important realization is that almost all of the missing pivots
	occur within the last M columns of the GE matrix, even for large matrices.
	So, the heavy matrix is always 6xM, where M is around 12.  Since the heavy
	matrix never gets any larger, the execution time doesn't vary based on N,
	and for large enough N it only lowers throughput imperceptibly while still
	providing a huge reduction in fail rate.

		The overall check matrix structure can be visualized as:

		+-----------+-----------+---------------------+
		|           |           |                     |
		|  Dense    |  Dense    |    Dense Mixing     |
		|  Deferred |  Mixing   |    Heavy Overlap    |
		|           |           |                     |
		+-----------+-----------+---------------------+
					|           |                     |
			Zero	|  Deferred |    Deferred Mixing  |
			-ish    |  Mixing   |    Heavy Overlap    |
					|           |                     |
		+-----------+-----------+---------------------+
		|                       |                     |
		|   Extra Binary Rows   |   Extra Heavy Rows  | <-- Uninitialized
		|                       |                     |
		+-----------------------+-------------+-------+
								|             |       |
			Implicitly Zero     |      H      |   I   | <-- 6x6 Identity matrix
			(Not Allocated)		|             |       |
								+-------------+-------+

		The heavy matrix H can be chosen so that each element is random.
	To synchronize it between the encoder and decoder, a PRNG is seeded
	based on the block count and H is generated from the PRNG output.
	I toyed with giving these rows more structure but I have not found
	anything worth using yet.
*/

void Codec::SetHeavyRows()
{
	CAT_IF_DUMP(cout << endl << "---- SetHeavyRows ----" << endl << endl;)

	Abyssinian prng;
	prng.Initialize(_p_seed);

	// Skip extra rows
	u8 * CAT_RESTRICT heavy_offset = _heavy_matrix + _heavy_pitch * _extra_count;

	// For each heavy matrix word,
	u8 * CAT_RESTRICT heavy_row = heavy_offset;
	//u32 words = CAT_HEAVY_ROWS * (_heavy_pitch / sizeof(u32));
	for (int row_i = 0; row_i < CAT_HEAVY_ROWS; ++row_i, heavy_row += _heavy_pitch)
	{
		// NOTE: Each heavy row is a multiple of 4 bytes in size
		u32 * CAT_RESTRICT words = reinterpret_cast<u32*>( heavy_row );
		for (int col_i = 0; col_i < _heavy_columns; col_i += 4)
			*words++ = prng.Next();
	}

	// Add identity matrix to tie heavy rows to heavy mixing columns
	u8 * CAT_RESTRICT lower_right = heavy_offset + _heavy_columns - CAT_HEAVY_ROWS;
	for (int ii = 0; ii < CAT_HEAVY_ROWS; ++ii, lower_right += _heavy_pitch)
		for (int jj = 0; jj < CAT_HEAVY_ROWS; ++jj)
			lower_right[jj] = (ii == jj) ? 1 : 0;
}

/*
		One more subtle optimization.  Why not?  Depending on how the GE
	matrix is constructed, it can be put in a roughly upper-triangular
	form from the start so it looks like this:

		+-------+-------+
		| D D D | D D 1 |
		| D D D | D 1 M | <- Dense rows
		| D D D | 1 M M |
		+-------+-------+
		| 0 0 1 | M M M |
		| 0 0 0 | M M M | <- Sparse deferred rows
		| 0 0 0 | M M M |
		+-------+-------+
		    ^       ^------- Mixing columns
		    \--------------- Deferred columns

		In the example above, the top 4 rows are dense matrix rows.
	The last 4 columns are mixing columns and are also dense.
	The lower left sub-matrix is sparse and roughly upper triangular
	because it is the intersection of sparse rows in the generator
	matrix.  This form is achieved by adding deferred rows starting
	from last deferred to first, and adding deferred columns starting
	from last deferred to first.

		Gaussian elimination will proceed from left to right on the
	matrix.  So, having an upper triangular form will prevent left-most
	zeroes from being eaten up when a row is eliminated by one above it.
	This reduces row operations.  For example, GE on a 64x64 matrix will
	do on average 100 fewer row operations on this form rather than its
	transpose (where the sparse part is put in the upper right).
*/


//// (3) Gaussian Elimination

/*
	SetupTriangle

		This function initializes the state variables to begin doing
	triangularization.  It was split from Triangle() at one point during
	development so that Triangle() can be shared with ResumeSolveMatrix(),
	which saves a lot of code duplication.

		The pivot list is longer than the number of columns to determine so
	that it will be useful for keeping track of extra overhead blocks.
	Initially it only contains non-heavy rows.
*/

void Codec::SetupTriangle()
{
	CAT_IF_DUMP(cout << endl << "---- SetupTriangle ----" << endl << endl;)

	// Initialize pivot array to just non-heavy rows
	const u16 pivot_count = _defer_count + _dense_count;
	for (u16 pivot_i = 0; pivot_i < pivot_count; ++pivot_i)
		_pivots[pivot_i] = pivot_i;

	// Set resume point to the first column
	_next_pivot = 0;
	_pivot_count = pivot_count;

	// If heavy rows are used right from the start,
	if (_first_heavy_column <= 0)
		InsertHeavyRows();
}

/*
	InsertHeavyRows

		This function converts remaining extra rows to heavy rows and adds
	heavy rows to the GE matrix.
*/
void Codec::InsertHeavyRows()
{
	CAT_IF_DUMP(cout << endl << "---- InsertHeavyRows ----" << endl << endl;)

	CAT_IF_DUMP(cout << "Converting remaining extra rows to heavy...";)

	// Initialize index of first heavy pivot
	u16 first_heavy_pivot = _pivot_count;

	// For each remaining pivot in the list,
	const u16 column_count = _defer_count + _mix_count;
	const u16 first_heavy_row = _defer_count + _dense_count;
	for (int pivot_j = _pivot_count - 1; pivot_j >= 0; --pivot_j)
	{
		// If row is extra,
		u16 ge_row_j = _pivots[pivot_j];
		if (ge_row_j < first_heavy_row)
			continue;

		// If pivot is still unused,
		if (pivot_j >= _next_pivot)
		{
			// Swap pivot j into last heavy pivot position
			--first_heavy_pivot;
			_pivots[pivot_j] = _pivots[first_heavy_pivot];
			_pivots[first_heavy_pivot] = ge_row_j;
		}

		CAT_IF_DUMP(cout << " row=" << ge_row_j << ", pivot=" << pivot_j;)

		// Copy heavy columns to heavy matrix row
		u8 * CAT_RESTRICT extra_row = _heavy_matrix + _heavy_pitch * (ge_row_j - first_heavy_row);
		u64 * CAT_RESTRICT ge_extra_row = _ge_matrix + _ge_pitch * ge_row_j;
		for (u16 ge_column_j = _first_heavy_column; ge_column_j < column_count; ++ge_column_j)
		{
			extra_row[ge_column_j - _first_heavy_column] = (ge_extra_row[ge_column_j >> 6] >> (ge_column_j & 63)) & 1;
		}
	}

	CAT_IF_DUMP(cout << endl;)

	// Store first heavy pivot index
	_first_heavy_pivot = first_heavy_pivot;

	// Add heavy rows at the end to cause them to be selected last if given a choice
	for (u16 heavy_i = 0; heavy_i < CAT_HEAVY_ROWS; ++heavy_i)
	{
		// Use GE row index after extra count even if not all are used yet
		_pivots[_pivot_count + heavy_i] = first_heavy_row + _extra_count + heavy_i;
	}
	_pivot_count += CAT_HEAVY_ROWS;

	CAT_IF_DUMP(cout << "Added heavy rows to the end of the pivots list." << endl;)
}

/*
	TriangleNonHeavy

		This function performs triangularization for all of the non-heavy
	columns of the GE matrix.  As soon as the first heavy column needs to
	be determined it converts the extra rows to heavy rows and adds in the
	heavy rows of the matrix to the pivot list.
*/

bool Codec::TriangleNonHeavy()
{
	CAT_IF_DUMP(cout << endl << "---- TriangleNonHeavy ----" << endl << endl;)

	const u16 pivot_count = _pivot_count;
	const u16 first_heavy_column = _first_heavy_column;

	// For the columns that are not protected by heavy rows,
	u16 pivot_i = _next_pivot;
	u64 ge_mask = (u64)1 << (pivot_i & 63);
	for (; pivot_i < first_heavy_column; ++pivot_i)
	{
		const int word_offset = pivot_i >> 6;

		bool found = false;

		// For each remaining GE row that might be the pivot,
		u64 * CAT_RESTRICT ge_matrix_offset = _ge_matrix + word_offset;
		for (u16 pivot_j = pivot_i; pivot_j < pivot_count; ++pivot_j)
		{
			// Determine if the row contains the bit we want
			u16 ge_row_j = _pivots[pivot_j];

			// If the bit was not found,
			u64 * CAT_RESTRICT ge_row = &ge_matrix_offset[_ge_pitch * ge_row_j];
			if (!(*ge_row & ge_mask)) continue; // Skip to next

			// Found it!
			found = true;
			CAT_IF_DUMP(cout << "Pivot " << pivot_i << " found on row " << ge_row_j << endl;)

			// Swap out the pivot index for this one
			_pivots[pivot_j] = _pivots[pivot_i];
			_pivots[pivot_i] = ge_row_j;

			// Prepare masked first word
			u64 row0 = (*ge_row & ~(ge_mask - 1)) ^ ge_mask;

			// For each remaining unused row,
			for (u16 pivot_k = pivot_j + 1; pivot_k < pivot_count; ++pivot_k)
			{
				// Determine if the row contains the bit we want
				u16 ge_row_k = _pivots[pivot_k];
				u64 * CAT_RESTRICT rem_row = &ge_matrix_offset[_ge_pitch * ge_row_k];

				// If the bit was found,
				if (*rem_row & ge_mask)
				{
					// Unroll first word to handle masked word and for speed
					*rem_row ^= row0;

					// Add the pivot row to eliminate the bit from this row, preserving previous bits
					for (int ii = 1; ii < _ge_pitch - word_offset; ++ii)
						rem_row[ii] ^= ge_row[ii];
				}
			} // next remaining row

			break;
		}

		// If pivot could not be found,
		if (!found)
		{
			_next_pivot = pivot_i;
			CAT_IF_DUMP(cout << "Singular: Pivot " << pivot_i << " of " << (_defer_count + _mix_count) << " not found!" << endl;)
			CAT_IF_PIVOT(if (pivot_i + 16 < (_defer_count + _mix_count)) cout << ">>>>> Singular: Pivot " << pivot_i << " of " << (_defer_count + _mix_count) << " not found!" << endl << endl;)
			return false;
		}

		// Generate next mask
		ge_mask = CAT_ROL64(ge_mask, 1);
	}

	_next_pivot = pivot_i;

	InsertHeavyRows();

	return true;
}

/*
	Triangle

		This function performs normal Gaussian elimination on the GE
	matrix in order to put it in upper triangular form, which would
	solve the system of linear equations represented by the matrix.
	The only wrinkle here is that instead of zeroing the columns as
	it goes, this algorithm will keep a record of which columns were
	added together so that these steps can be followed later to
	produce the column values.  Note that memxor() is not used here.

		It uses a pivot array that it swaps row pointers around in, as
	opposed to actually swapping rows in the matrix, which would be
	more expensive.  Heavy rows are kept at the end of the pivot array
	so that they are always selected last.
*/

#if defined(CAT_HEAVY_WIN_MULT)

// Flip endianness at compile time if possible
#if defined(CAT_ENDIAN_BIG)
static u32 GF256_MULT_LOOKUP[16] = {
	0x00000000, 0x01000000, 0x00010000, 0x01010000, 
	0x00000100, 0x01000100, 0x00010100, 0x01010100, 
	0x00000001, 0x01000001, 0x00010001, 0x01010001, 
	0x00000101, 0x01000101, 0x00010101, 0x01010101, 
};
#else // Little-endian or unknown bit order:
static u32 GF256_MULT_LOOKUP[16] = {
	0x00000000, 0x00000001, 0x00000100, 0x00000101, 
	0x00010000, 0x00010001, 0x00010100, 0x00010101, 
	0x01000000, 0x01000001, 0x01000100, 0x01000101, 
	0x01010000, 0x01010001, 0x01010100, 0x01010101, 
};
#endif

#endif // CAT_HEAVY_WIN_MULT

bool Codec::Triangle()
{
	CAT_IF_DUMP(cout << endl << "---- Triangle ----" << endl << endl;)

	const u16 first_heavy_column = _first_heavy_column;

	// If next pivot is not heavy,
	if (_next_pivot < first_heavy_column && !TriangleNonHeavy())
		return false;

	const u16 pivot_count = _pivot_count;
	const u16 column_count = _defer_count + _mix_count;
	const u16 first_heavy_row = _defer_count + _dense_count;
	u16 first_heavy_pivot = _first_heavy_pivot;

	// For each heavy pivot to determine,
	u64 ge_mask = (u64)1 << (_next_pivot & 63);
	for (u16 pivot_i = _next_pivot; pivot_i < column_count;
		++pivot_i, ge_mask = CAT_ROL64(ge_mask, 1))
	{
		const u16 heavy_col_i = pivot_i - first_heavy_column;

		// For each remaining GE row that might be the pivot,
		int word_offset = pivot_i >> 6;
		u64 * CAT_RESTRICT ge_matrix_offset = _ge_matrix + word_offset;
		bool found = false;
		u16 pivot_j;
		for (pivot_j = pivot_i; pivot_j < first_heavy_pivot; ++pivot_j)
		{
			// If the bit was not found,
			u16 ge_row_j = _pivots[pivot_j];
			u64 * CAT_RESTRICT ge_row = &ge_matrix_offset[_ge_pitch * ge_row_j];
			if (!(*ge_row & ge_mask)) continue; // Skip to next

			// Found it!
			found = true;
			CAT_IF_DUMP(cout << "Pivot " << pivot_i << " found on row " << ge_row_j << endl;)

			// Swap out the pivot index for this one
			_pivots[pivot_j] = _pivots[pivot_i];
			_pivots[pivot_i] = ge_row_j;

			// Prepare masked first word
			u64 row0 = (*ge_row & ~(ge_mask - 1)) ^ ge_mask;

			// For each remaining light row,
			u16 pivot_k = pivot_j + 1;
			for (; pivot_k < first_heavy_pivot; ++pivot_k)
			{
				// Determine if the row contains the bit we want
				u16 ge_row_k = _pivots[pivot_k];
				u64 * CAT_RESTRICT rem_row = &ge_matrix_offset[_ge_pitch * ge_row_k];

				// If the bit was found,
				if (*rem_row & ge_mask)
				{
					// Unroll first word to handle masked word and for speed
					*rem_row ^= row0;

					// Add the pivot row to eliminate the bit from this row, preserving previous bits
					for (int ii = 1; ii < _ge_pitch - word_offset; ++ii)
						rem_row[ii] ^= ge_row[ii];
				}
			} // next remaining row

			// For each remaining heavy row,
			for (; pivot_k < pivot_count; ++pivot_k)
			{
				// If the column is non-zero,
				u16 heavy_row_k = _pivots[pivot_k] - first_heavy_row;
				u8 * CAT_RESTRICT rem_row = &_heavy_matrix[_heavy_pitch * heavy_row_k];
				u8 code_value = rem_row[heavy_col_i];
				if (!code_value) continue;

				CAT_IF_DUMP(cout << "Eliminating from heavy row " << heavy_row_k << " :";)

				// For each set bit in the binary pivot row, add rem[i] to rem[i+]:
				u64 * CAT_RESTRICT pivot_row = &_ge_matrix[_ge_pitch * ge_row_j];

#if !defined(CAT_HEAVY_WIN_MULT)
				for (int ge_column_i = pivot_i + 1; ge_column_i < column_count; ++ge_column_i)
				{
					if (pivot_row[ge_column_i >> 6] & ((u64)1 << (ge_column_i & 63)))
					{
						rem_row[ge_column_i - first_heavy_column] ^= code_value;
					}
				}
#else // CAT_HEAVY_WIN_MULT
				// Unroll odd columns:
				u16 odd_count = pivot_i & 3, ge_column_i = pivot_i + 1;
				u64 temp_mask = ge_mask;
				switch (odd_count)
				{
				case 0:	temp_mask = CAT_ROL64(temp_mask, 1);
						if (pivot_row[ge_column_i >> 6] & temp_mask)
						{
							rem_row[ge_column_i - _first_heavy_column] ^= code_value;
							CAT_IF_DUMP(cout << " " << ge_column_i;)
						}
						++ge_column_i;

				case 1: temp_mask = CAT_ROL64(temp_mask, 1);
						if (pivot_row[ge_column_i >> 6] & temp_mask)
						{
							rem_row[ge_column_i - _first_heavy_column] ^= code_value;
							CAT_IF_DUMP(cout << " " << ge_column_i;)
						}
						++ge_column_i;

				case 2:	temp_mask = CAT_ROL64(temp_mask, 1);
						if (pivot_row[ge_column_i >> 6] & temp_mask)
						{
							rem_row[ge_column_i - _first_heavy_column] ^= code_value;
							CAT_IF_DUMP(cout << " " << ge_column_i;)
						}

						// Set GE column to next even column
						ge_column_i = pivot_i + (4 - odd_count);
				}

				// For remaining aligned columns,
				u32 * CAT_RESTRICT word = reinterpret_cast<u32*>( rem_row + ge_column_i - first_heavy_column );
				for (; ge_column_i < column_count; ge_column_i += 4, ++word)
				{
					// Look up 4 bit window
					u32 bits = (u32)(pivot_row[ge_column_i >> 6] >> (ge_column_i & 63)) & 15;
#if defined(CAT_ENDIAN_UNKNOWN)
					u32 window = getLE(GF256_MULT_LOOKUP[bits]);
#else
					u32 window = GF256_MULT_LOOKUP[bits];
#endif

					CAT_IF_DUMP(cout << " " << ge_column_i << "x" << hex << setw(8) << setfill('0') << window << dec;)

					*word ^= window * code_value;
				}
#endif // CAT_HEAVY_WIN_MULT
				CAT_IF_DUMP(cout << endl;)
			} // next heavy row

			break;
		} // next row

		// If not found, then for each remaining heavy pivot,
		// NOTE: The pivot array maintains the heavy rows at the end so that they are always tried last
		if (!found) for (; pivot_j < _pivot_count; ++pivot_j)
		{
			// If heavy row doesn't have the pivot,
			u16 ge_row_j = _pivots[pivot_j];
			u16 heavy_row_j = ge_row_j - first_heavy_row;
			u8 * CAT_RESTRICT pivot_row = &_heavy_matrix[_heavy_pitch * heavy_row_j];
			u8 code_value = pivot_row[heavy_col_i];
			if (!code_value) continue; // Skip to next

			// Found it! (common case)
			found = true;
			CAT_IF_DUMP(cout << "Pivot " << pivot_i << " found on heavy row " << ge_row_j << endl;)

			// Swap pivot i and j
			_pivots[pivot_j] = _pivots[pivot_i];
			_pivots[pivot_i] = ge_row_j;

			// If a non-heavy pivot just got moved into heavy pivot list,
			if (pivot_i < first_heavy_pivot)
			{
				// Swap pivot j with first heavy pivot
				u16 temp = _pivots[first_heavy_pivot];
				_pivots[first_heavy_pivot] = _pivots[pivot_j];
				_pivots[pivot_j] = temp;

				// And move the first heavy pivot up one to cover the hole
				++first_heavy_pivot;
			}

			// If there are any remaining rows,
			u16 pivot_k = pivot_j + 1;
			if (pivot_k < pivot_count)
			{
				// Precompute denominator
				int denominator = 255 - LOG_TABLE[code_value];

				// For each remaining unused row,
				// NOTE: All remaining rows are heavy rows by pivot array organization
				for (; pivot_k < pivot_count; ++pivot_k)
				{
					// If the column is zero,
					u16 ge_row_k = _pivots[pivot_k];
					u16 heavy_row_k = ge_row_k - first_heavy_row;
					u8 * CAT_RESTRICT rem_row = &_heavy_matrix[_heavy_pitch * heavy_row_k];
					u8 rem_value = rem_row[heavy_col_i];
					if (!rem_value) continue; // Skip it

					// x = rem_value / code_value
					u8 x = EXP_TABLE[LOG_TABLE[rem_value] + denominator];

					// Store value for later
					rem_row[heavy_col_i] = x;

					// rem[i+] += x * pivot[i+]
					const int offset = heavy_col_i + 1;
					GF256MemMulAdd(rem_row + offset, x, pivot_row + offset, _heavy_columns - offset);
				} // next remaining row
			}

			break;
		} // next heavy row

		// If pivot could not be found,
		if (!found)
		{
			_next_pivot = pivot_i;
			_first_heavy_pivot = first_heavy_pivot;
			CAT_IF_DUMP(cout << "Singular: Pivot " << pivot_i << " of " << (_defer_count + _mix_count) << " not found!" << endl;)
			CAT_IF_PIVOT(if (pivot_i + 16 < (_defer_count + _mix_count)) cout << ">>>>> Singular: Pivot " << pivot_i << " of " << (_defer_count + _mix_count) << " not found!" << endl << endl; )
			return false;
		}

		CAT_IF_DUMP( PrintGEMatrix(); )
		CAT_IF_DUMP( PrintExtraMatrix(); )
	}

	return true;
}


//// (4) Substitute

/*
	InitializeColumnValues

		This function initializes the output block value for each column
	that was solved by Gaussian elimination.  For deferred rows it follows
	the same steps performed earlier in PeelDiagonal() just for those rows.
	The row values were not formed at that point because the destination
	was uncertain.

	For each pivot that solves the GE matrix,
		If GE row is from a dense row,
			Initialize row value to 0.
		Else it is from a deferred row,
			For each peeled column that it references,
				Add in that peeled column's row value from Compression.

	For each remaining unused row, (happens in the decoder for extra rows)
		If the unused row is a dense row,
			Set the GE row map entry to LIST_TERM so it can be ignored later.
*/

void Codec::InitializeColumnValues()
{
	CAT_IF_DUMP(cout << endl << "---- InitializeColumnValues ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	const u16 first_heavy_row = _defer_count + _dense_count;
	const u16 column_count = _defer_count + _mix_count;

	// For each pivot,
	u16 pivot_i;
	for (pivot_i = 0; pivot_i < column_count; ++pivot_i)
	{
		// Lookup pivot column, GE row, and destination buffer
		u16 dest_column_i = _ge_col_map[pivot_i];
		u16 ge_row_i = _pivots[pivot_i];
		u8 * CAT_RESTRICT buffer_dest = _recovery_blocks + _block_bytes * dest_column_i;

		CAT_IF_DUMP(cout << "Pivot " << pivot_i << " solving column " << dest_column_i << " with GE row " << ge_row_i << " : ";)

		// If it is a dense/heavy(non-extra) row,
		if (ge_row_i < _dense_count ||
			ge_row_i >= (first_heavy_row + _extra_count))
		{
			// Dense/heavy rows sum to zero
			memset(buffer_dest, 0, _block_bytes);

			// Store which column solves the dense row
			_ge_row_map[ge_row_i] = dest_column_i;

			CAT_IF_DUMP(cout << "[0]" << endl;)
			CAT_IF_ROWOP(++rowops;)

			continue;
		}

		// Look up row and input value for GE row
		u16 row_i = _ge_row_map[ge_row_i];
		const u8 * CAT_RESTRICT combo = _input_blocks + _block_bytes * row_i;
		PeelRow * CAT_RESTRICT row = &_peel_rows[row_i];

		CAT_IF_DUMP(cout << "[" << (int)combo[0] << "]";)

		// If copying from final input block,
		if (row_i == _block_count - 1)
		{
			memcpy(buffer_dest, combo, _input_final_bytes);
			memset(buffer_dest + _input_final_bytes, 0, _block_bytes - _input_final_bytes);
			CAT_IF_ROWOP(++rowops;)
			combo = 0;
		}

		// Eliminate peeled columns:
		u16 column_i = row->peel_x0;
		u16 a = row->peel_a;
		u16 weight = row->peel_weight;
		for (;;)
		{
			// If column is peeled,
			PeelColumn * CAT_RESTRICT column = &_peel_cols[column_i];
			if (column->mark == MARK_PEEL)
			{
				// If combo unused,
				if (!combo)
					memxor(buffer_dest, _recovery_blocks + _block_bytes * column_i, _block_bytes);
				else
				{
					// Use combo
					memxor_set(buffer_dest, combo, _recovery_blocks + _block_bytes * column_i, _block_bytes);
					combo = 0;
				}
				CAT_IF_ROWOP(++rowops;)
			}

			if (--weight <= 0) break;

			IterateNextColumn(column_i, _block_count, _block_next_prime, a);
		}

		// If combo still unused,
		if (combo) memcpy(buffer_dest, combo, _block_bytes);
		CAT_IF_DUMP(cout << endl;)
	}

	// For each remaining pivot,
	for (; pivot_i < _pivot_count; ++pivot_i)
	{
		u16 ge_row_i = _pivots[pivot_i];

		// If row is a dense row,
		if (ge_row_i < _dense_count ||
			(ge_row_i >= first_heavy_row && ge_row_i < column_count))
		{
			// Mark it for skipping
			_ge_row_map[ge_row_i] = LIST_TERM;

			CAT_IF_DUMP(cout << "Did not use GE row " << ge_row_i << ", which is a dense row." << endl;)
		}
		else
		{
			CAT_IF_DUMP(cout << "Did not use deferred row " << ge_row_i << ", which is not a dense row." << endl;)
		}
	}

	CAT_IF_ROWOP(cout << "InitializeColumnValues used " << rowops << " row ops = " << rowops / (double)_block_count << "*N" << endl;)
}

/*
	MultiplyDenseValues

		This function follows the same order of operations as the
	MultiplyDenseRows() function to add in peeling column values
	to the dense rows.  Deferred rows are handled above in the
	InitializeColumnValues() function.

		See MultiplyDenseRows() comments for justification of the
	design of the dense row structure.
*/

void Codec::MultiplyDenseValues()
{
	CAT_IF_DUMP(cout << endl << "---- MultiplyDenseValues ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// Initialize PRNG
	Abyssinian prng;
	prng.Initialize(_d_seed);

	// For each block of columns,
	const int dense_count = _dense_count;
	u8 * CAT_RESTRICT temp_block = _recovery_blocks + _block_bytes * (_block_count + _mix_count);
	const u8 * CAT_RESTRICT source_block = _recovery_blocks;
	PeelColumn * CAT_RESTRICT column = _peel_cols;
	u16 rows[CAT_MAX_DENSE_ROWS], bits[CAT_MAX_DENSE_ROWS];
	const u16 block_count = _block_count;
	for (u16 column_i = 0; column_i < block_count; column_i += dense_count,
		column += dense_count, source_block += _block_bytes * dense_count)
	{
		// Handle final columns
		int max_x = dense_count;
		if (column_i + dense_count > block_count)
			max_x = _block_count - column_i;

		CAT_IF_DUMP(cout << endl << "For window of columns between " << column_i << " and " << column_i + dense_count - 1 << " (inclusive):" << endl;)

		// Shuffle row and bit order
		ShuffleDeck16(prng, rows, dense_count);
		ShuffleDeck16(prng, bits, dense_count);

		// Initialize counters
		u16 set_count = (dense_count + 1) >> 1;
		u16 * CAT_RESTRICT set_bits = bits;
		u16 * CAT_RESTRICT clr_bits = set_bits + set_count;
		const u16 * CAT_RESTRICT row = rows;

		CAT_IF_DUMP(cout << "Generating first row " << _ge_row_map[*row] << ":";)

		// Generate first row
		const u8 * CAT_RESTRICT combo = 0;
		CAT_IF_ROWOP(++rowops;)
		for (int ii = 0; ii < set_count; ++ii)
		{
			// If bit is peeled,
			int bit_i = set_bits[ii];
			if (bit_i < max_x && column[bit_i].mark == MARK_PEEL)
			{
				CAT_IF_DUMP(cout << " " << column_i + bit_i;)

				// If no combo used yet,
				const u8 * CAT_RESTRICT src = source_block + _block_bytes * bit_i;
				if (!combo)
					combo = src;
				else if (combo == temp_block)
				{
					// Else if combo has been used: XOR it in
					memxor(temp_block, src, _block_bytes);
					CAT_IF_ROWOP(++rowops;)
				}
				else
				{
					// Else if combo needs to be used: Combine into block
					memxor_set(temp_block, combo, src, _block_bytes);
					CAT_IF_ROWOP(++rowops;)
					combo = temp_block;
				}
			}
		}

		CAT_IF_DUMP(cout << endl;)

		// If no combo ever triggered,
		if (!combo)
			memset(temp_block, 0, _block_bytes);
		else
		{
			// Else if never combined two: Just copy it
			if (combo != temp_block)
			{
				memcpy(temp_block, combo, _block_bytes);
				CAT_IF_ROWOP(++rowops;)
			}

			// Store in destination column in recovery blocks
			u16 dest_column_i = _ge_row_map[*row];
			if (dest_column_i != LIST_TERM)
			{
				memxor(_recovery_blocks + _block_bytes * dest_column_i, temp_block, _block_bytes);
				CAT_IF_ROWOP(++rowops;)
			}
		}
		++row;

		// Reshuffle bit order: Shuffle-2 Code
		ShuffleDeck16(prng, bits, dense_count);

		// Generate first half of rows
		const int loop_count = (dense_count >> 1);
		for (int ii = 0; ii < loop_count; ++ii)
		{
			CAT_IF_DUMP(cout << "Flipping bits for derivative row " << _ge_row_map[*row] << ":";)

			// Add in peeled columns
			int bit0 = set_bits[ii], bit1 = clr_bits[ii];
			if (bit0 < max_x && column[bit0].mark == MARK_PEEL)
			{
				if (bit1 < max_x && column[bit1].mark == MARK_PEEL)
				{
					CAT_IF_DUMP(cout << " " << column_i + bit0 << "+" << column_i + bit1;)
					memxor_add(temp_block, source_block + _block_bytes * bit0, source_block + _block_bytes * bit1, _block_bytes);
				}
				else
				{
					CAT_IF_DUMP(cout << " " << column_i + bit0;)
					memxor(temp_block, source_block + _block_bytes * bit0, _block_bytes);
				}
				CAT_IF_ROWOP(++rowops;)
			}
			else if (bit1 < max_x && column[bit1].mark == MARK_PEEL)
			{
				CAT_IF_DUMP(cout << " " << column_i + bit1;)
				memxor(temp_block, source_block + _block_bytes * bit1, _block_bytes);
				CAT_IF_ROWOP(++rowops;)
			}

			CAT_IF_DUMP(cout << endl;)

			// Store in destination column in recovery blocks
			u16 dest_column_i = _ge_row_map[*row++];
			if (dest_column_i != LIST_TERM)
			{
				memxor(_recovery_blocks + _block_bytes * dest_column_i, temp_block, _block_bytes);
				CAT_IF_ROWOP(++rowops;)
			}
		}

		// Reshuffle bit order: Shuffle-2 Code
		ShuffleDeck16(prng, bits, dense_count);

		// Generate second half of rows
		const int second_loop_count = loop_count - 1 + (dense_count & 1);
		for (int ii = 0; ii < second_loop_count; ++ii)
		{
			CAT_IF_DUMP(cout << "Flipping bits for derivative row " << _ge_row_map[*row] << ":";)

			// Add in peeled columns
			int bit0 = set_bits[ii], bit1 = clr_bits[ii];
			if (bit0 < max_x && column[bit0].mark == MARK_PEEL)
			{
				if (bit1 < max_x && column[bit1].mark == MARK_PEEL)
				{
					CAT_IF_DUMP(cout << " " << column_i + bit0 << "+" << column_i + bit1;)
					memxor_add(temp_block, source_block + _block_bytes * bit0, source_block + _block_bytes * bit1, _block_bytes);
				}
				else
				{
					CAT_IF_DUMP(cout << " " << column_i + bit0;)
					memxor(temp_block, source_block + _block_bytes * bit0, _block_bytes);
				}
				CAT_IF_ROWOP(++rowops;)
			}
			else if (bit1 < max_x && column[bit1].mark == MARK_PEEL)
			{
				CAT_IF_DUMP(cout << " " << column_i + bit1;)
				memxor(temp_block, source_block + _block_bytes * bit1, _block_bytes);
				CAT_IF_ROWOP(++rowops;)
			}

			CAT_IF_DUMP(cout << endl;)

			// Store in destination column in recovery blocks
			u16 dest_column_i = _ge_row_map[*row++];
			if (dest_column_i != LIST_TERM)
			{
				memxor(_recovery_blocks + _block_bytes * dest_column_i, temp_block, _block_bytes);
				CAT_IF_ROWOP(++rowops;)
			}
		}
	} // next column

	CAT_IF_ROWOP(cout << "MultiplyDenseValues used " << rowops << " row ops = " << rowops / (double)_block_count << "*N" << endl;)
}

/*
	AddSubdiagonalValues

		This function uses the bits that were left behind by the
	Triangle() function to follow the same order of operations
	to generate the row values for both deferred and dense rows.
	It is aided by the already roughly upper-triangular form
	of the GE matrix, making this function very cheap to execute.
*/

// These are heuristic values.  Choosing better values has little effect on performance.
#define CAT_UNDER_WIN_THRESH_4 (45 + 4) /* Note: Assumes this is higher than CAT_HEAVY_MAX_COLS */
#define CAT_UNDER_WIN_THRESH_5 (65 + 5)
#define CAT_UNDER_WIN_THRESH_6 (85 + 6)
#define CAT_UNDER_WIN_THRESH_7 (138 + 7)

void Codec::AddSubdiagonalValues()
{
	CAT_IF_DUMP(cout << endl << "---- AddSubdiagonalValues ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0; int heavyops = 0;)

	const int column_count = _defer_count + _mix_count;
	int pivot_i = 0;
	const u16 first_heavy_row = _defer_count + _dense_count;

#if defined(CAT_WINDOWED_LOWERTRI)
	const u16 first_non_binary_row = first_heavy_row + _extra_count;

	// Build temporary storage space if windowing is to be used
	if (column_count >= CAT_UNDER_WIN_THRESH_5)
	{
		// Calculate initial window size
		int w, next_check_i;
		if (column_count >= CAT_UNDER_WIN_THRESH_7)
		{
			w = 7;
			next_check_i = column_count - CAT_UNDER_WIN_THRESH_7;
		}
		else if (column_count >= CAT_UNDER_WIN_THRESH_6)
		{
			w = 6;
			next_check_i = column_count - CAT_UNDER_WIN_THRESH_6;
		}
		else if (column_count >= CAT_UNDER_WIN_THRESH_5)
		{
			w = 5;
			next_check_i = column_count - CAT_UNDER_WIN_THRESH_5;
		}
		else
		{
			w = 4;
			next_check_i = CAT_UNDER_WIN_THRESH_4;
		}
		u32 win_lim = 1 << w;

		CAT_IF_DUMP(cout << "Activating windowed lower triangular elimination with initial window " << w << endl;)

		// Use the first few peel column values as window table space
		// NOTE: The peeled column values were previously used up until this point,
		// but now they are unused, and so they can be reused for temporary space.
		u8 * CAT_RESTRICT win_table[128];
		PeelColumn * CAT_RESTRICT column = _peel_cols;
		u8 * CAT_RESTRICT column_src = _recovery_blocks;
		u32 jj = 1;
		for (u32 count = _block_count; count > 0; --count, ++column, column_src += _block_bytes)
		{
			// If column is peeled,
			if (column->mark == MARK_PEEL)
			{
				// Reuse the block value temporarily as window table space
				win_table[jj] = column_src;

				CAT_IF_DUMP(cout << "-- Window table entry " << jj << " set to column " << _block_count - count << endl;)

				// If done,
				if (++jj >= win_lim) break;
			}
		}

		CAT_IF_DUMP(if (jj < win_lim) cout << "!! Not enough space in peeled columns to generate a table.  Going back to normal lower triangular elimination." << endl;)

		// If enough space was found,
		if (jj >= win_lim) for (;;)
		{
			// Calculate first column in window
			u16 final_i = pivot_i + w - 1;

			CAT_IF_DUMP(cout << "-- Windowing from " << pivot_i << " to " << final_i << " (inclusive)" << endl;)

			// Eliminate lower triangular part below windowed bits:

			// For each column,
			u64 ge_mask = (u64)1 << (pivot_i & 63);
			for (int src_pivot_i = pivot_i; src_pivot_i < final_i;
				++src_pivot_i, ge_mask = CAT_ROL64(ge_mask, 1))
			{
				u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * _ge_col_map[src_pivot_i];

				CAT_IF_DUMP(cout << "Back-substituting small triangle from pivot " << src_pivot_i << "[" << (int)src[0] << "] :";)

				// For each row above the diagonal,
				u64 * CAT_RESTRICT ge_row = _ge_matrix + (src_pivot_i >> 6);
				for (int dest_pivot_i = src_pivot_i + 1; dest_pivot_i <= final_i; ++dest_pivot_i)
				{
					// If row is heavy,
					u16 dest_row_i = _pivots[dest_pivot_i];

					// If bit is set in that row,
					if (ge_row[_ge_pitch * dest_row_i] & ge_mask)
					{
						// Back-substitute
						u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[dest_pivot_i];
						memxor(dest, src, _block_bytes);
						CAT_IF_ROWOP(++rowops;)

						CAT_IF_DUMP(cout << " " << dest_pivot_i;)
					}
				} // next pivot above

				CAT_IF_DUMP(cout << endl;)
			} // next pivot

			CAT_IF_DUMP(cout << "-- Generating window table with " << w << " bits" << endl;)

			// Generate window table: 2 bits
			win_table[1] = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i];
			win_table[2] = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i + 1];
			memxor_set(win_table[3], win_table[1], win_table[2], _block_bytes);
			CAT_IF_ROWOP(++rowops;)

			// Generate window table: 3 bits
			win_table[4] = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i + 2];
			memxor_set(win_table[5], win_table[1], win_table[4], _block_bytes);
			memxor_set(win_table[6], win_table[2], win_table[4], _block_bytes);
			memxor_set(win_table[7], win_table[1], win_table[6], _block_bytes);
			CAT_IF_ROWOP(rowops += 3;)

			// Generate window table: 4 bits
			win_table[8] = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i + 3];
			for (int ii = 1; ii < 8; ++ii)
				memxor_set(win_table[8 + ii], win_table[ii], win_table[8], _block_bytes);
			CAT_IF_ROWOP(rowops += 7;)

			// Generate window table: 5+ bits
			if (w >= 5)
			{
				win_table[16] = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i + 4];
				for (int ii = 1; ii < 16; ++ii)
					memxor_set(win_table[16 + ii], win_table[ii], win_table[16], _block_bytes);
				CAT_IF_ROWOP(rowops += 15;)

				if (w >= 6)
				{
					win_table[32] = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i + 5];
					for (int ii = 1; ii < 32; ++ii)
						memxor_set(win_table[32 + ii], win_table[ii], win_table[32], _block_bytes);
					CAT_IF_ROWOP(rowops += 31;)

					if (w >= 7)
					{
						win_table[64] = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i + 6];
						for (int ii = 1; ii < 64; ++ii)
							memxor_set(win_table[64 + ii], win_table[ii], win_table[64], _block_bytes);
						CAT_IF_ROWOP(rowops += 63;)
					}
				}
			}

			// If not straddling words,
			u32 first_word = pivot_i >> 6;
			u32 shift0 = pivot_i & 63;
			u32 last_word = final_i >> 6;
			u16 * CAT_RESTRICT pivot_row = _pivots + final_i + 1;
			if (first_word == last_word)
			{
				// For each pivot row,
				for (u16 ge_below_i = final_i + 1; ge_below_i < column_count; ++ge_below_i)
				{
					// If pivot row is heavy,
					u16 ge_row_i = *pivot_row++;
					if (ge_row_i >= first_non_binary_row) continue;

					// Calculate window bits
					u64 * CAT_RESTRICT ge_row = _ge_matrix + first_word + _ge_pitch * ge_row_i;
					u32 win_bits = (u32)(ge_row[0] >> shift0) & (win_lim - 1);

					// If any XOR needs to be performed,
					if (win_bits != 0)
					{
						CAT_IF_DUMP(cout << "Adding window table " << win_bits << " to pivot " << ge_below_i << endl;)

						// Back-substitute
						u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[ge_below_i];
						memxor(dest, win_table[win_bits], _block_bytes);
						CAT_IF_ROWOP(++rowops;)
					}
				}
			}
			else // Rare: Straddling case
			{
				u32 shift1 = 64 - shift0;

				// For each pivot row,
				for (u16 ge_below_i = final_i + 1; ge_below_i < column_count; ++ge_below_i)
				{
					// If pivot row is heavy,
					u16 ge_row_i = *pivot_row++;
					if (ge_row_i >= first_non_binary_row) continue;

					// Calculate window bits
					u64 * CAT_RESTRICT ge_row = _ge_matrix + first_word + _ge_pitch * ge_row_i;
					u32 win_bits = ( (u32)(ge_row[0] >> shift0) | (u32)(ge_row[1] << shift1) ) & (win_lim - 1);

					// If any XOR needs to be performed,
					if (win_bits != 0)
					{
						CAT_IF_DUMP(cout << "Adding window table " << win_bits << " to pivot " << ge_below_i << endl;)

						// Back-substitute
						u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[ge_below_i];
						memxor(dest, win_table[win_bits], _block_bytes);
						CAT_IF_ROWOP(++rowops;)
					}
				}
			}

			// If column index falls below window size,
			pivot_i += w;
			if (pivot_i >= next_check_i)
			{
				int remaining_columns = column_count - pivot_i;

				if (remaining_columns >= CAT_UNDER_WIN_THRESH_6)
				{
					w = 6;
					next_check_i = remaining_columns - CAT_UNDER_WIN_THRESH_6;
				}
				else if (remaining_columns >= CAT_UNDER_WIN_THRESH_5)
				{
					w = 5;
					next_check_i = remaining_columns - CAT_UNDER_WIN_THRESH_5;
				}
				else if (remaining_columns >= CAT_UNDER_WIN_THRESH_4)
				{
					w = 4;
					next_check_i = remaining_columns - CAT_UNDER_WIN_THRESH_4;
				}
				else break;

				// Update window limit
				win_lim = 1 << w;
			}
		} // next window
	} // end if windowed
#endif // CAT_WINDOWED_LOWERTRI

	// For each row to eliminate,
	for (u16 ge_column_i = pivot_i + 1; ge_column_i < column_count; ++ge_column_i)
	{
		// Lookup pivot column, GE row, and destination buffer
		u16 column_i = _ge_col_map[ge_column_i];
		u16 ge_row_i = _pivots[ge_column_i];
		u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * column_i;

		CAT_IF_DUMP(cout << "Pivot " << ge_column_i << " solving column " << column_i << "[" << (int)dest[0] << "] with GE row " << ge_row_i << " :";)

		u16 ge_limit = ge_column_i;

		// If row is heavy or extra,
		if (ge_row_i >= first_heavy_row)
		{
			u16 heavy_row_i = ge_row_i - first_heavy_row;

			// For each column up to the diagonal,
			u8 * CAT_RESTRICT heavy_row = _heavy_matrix + _heavy_pitch * heavy_row_i;
			for (u16 sub_i = _first_heavy_column; sub_i < ge_limit; ++sub_i)
			{
				// If column is zero,
				u8 code_value = heavy_row[sub_i - _first_heavy_column];
				if (!code_value) continue; // Skip it

				// Look up data source
				const u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * _ge_col_map[sub_i];

				if (code_value == 1)
				{
					memxor(dest, src, _block_bytes);
					CAT_IF_ROWOP(++rowops;)

					CAT_IF_DUMP(cout << " *" << ge_column_i << "=[" << (int)src[0] << "]";)
				}
				else
				{
					GF256MemMulAdd(dest, code_value, src, _block_bytes);
					CAT_IF_ROWOP(++heavyops;)

					CAT_IF_DUMP(cout << " h" << ge_column_i << "=[" << (int)src[0] << "*" << (int)code_value << "]";)
				}
			}

			// If row is not extra,
			if (heavy_row_i >= _extra_count)
			{
				CAT_IF_DUMP(cout << endl;)
				continue; // Skip binary matrix elimination
			}

			// Limit the binary matrix elimination to non-heavy columns
			if (ge_limit > _first_heavy_column)
				ge_limit = _first_heavy_column;

			// fall-thru..
		}

		// For each GE matrix bit in the row,
		u64 * CAT_RESTRICT ge_row = _ge_matrix + _ge_pitch * ge_row_i;
		u64 ge_mask = (u64)1 << (pivot_i & 63);
		for (u16 ge_sub_i = pivot_i; ge_sub_i < ge_limit; ++ge_sub_i, ge_mask = CAT_ROL64(ge_mask, 1))
		{
			// If bit is non-zero,
			if (ge_row[ge_sub_i >> 6] & ge_mask)
			{
				// Add pivot for non-zero bit to destination row value
				u16 column_i = _ge_col_map[ge_sub_i];
				const u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * column_i;
				memxor(dest, src, _block_bytes);
				CAT_IF_ROWOP(++rowops;)

				CAT_IF_DUMP(cout << " " << ge_sub_i << "=[" << (int)src[0] << "]";)
			}
		}

		CAT_IF_DUMP(cout << endl;)
	}

	CAT_IF_ROWOP(cout << "AddSubdiagonalValues used " << rowops << " row ops = " << rowops / (double)_block_count << "*N and " << heavyops << " heavy ops" << endl;)
}

/*
	Windowed Back-Substitution

		The matrix to diagonalize is in upper triangular form now, where each
	row is random and dense.  Each column has a value assigned to it at this
	point but it is necessary to back-substitute up each column to eliminate
	the upper triangular part of the matrix.  Because the matrix is random,
	the optimal window size for a windowed substitution method is approx:

		w = CEIL[0.85 + 0.85*ln(r)]

		But in practice, the window size is selected from simple heuristic
	rules based on testing.  For this function, the window size stays between
	3 and 6 so it is fine to use heuristics.  The matrix may look like:

		+---+---+---+---+
		| A | 1 | 1 | 0 |
		+---+---+---+---+
		| 0 | B | 1 | 1 |
		+---+---+---+---+
		| 0 | 0 | C | 1 |
		+---+---+---+---+
		| 0 | 0 | 0 | D |
		+---+---+---+---+

		To do the back-substitution with a window width of 2 on the above
	matrix, first back-substitute to diagonalize the lower right:

		+---+---+---+---+
		| A | 1 | 1 | 0 |
		+---+---+---+---+
		| 0 | B | 1 | 1 |
		+---+---+---+---+
		| 0 | 0 | C | 0 |
		+---+---+---+---+
		| 0 | 0 | 0 | D |
		+---+---+---+---+

		Then compute the 2-bit window:

		[0 0] = undefined
		[1 0] = C
		[0 1] = D
		[1 1] = C + D

		And substitute up the last two columns to eliminate them:

		B = B + [1 1] = B + C + D
		A = A + [1 0] = A + C

		+---+---+---+---+
		| A | 1 | 0 | 0 |
		+---+---+---+---+
		| 0 | B | 0 | 0 |
		+---+---+---+---+
		| 0 | 0 | C | 0 |
		+---+---+---+---+
		| 0 | 0 | 0 | D |
		+---+---+---+---+

		This operation is performed until the windowing method is no
	longer worthwhile, and the normal back-substitution is used on the
	remaining matrix pivots.
*/

// These are heuristic values.  Choosing better values has little effect on performance.
#define CAT_ABOVE_WIN_THRESH_4 (20 + 4)
#define CAT_ABOVE_WIN_THRESH_5 (40 + 5)
#define CAT_ABOVE_WIN_THRESH_6 (64 + 6)
#define CAT_ABOVE_WIN_THRESH_7 (128 + 7)

/*
	BackSubstituteAboveDiagonal

		This function uses the windowed approach outlined above
	to eliminate all of the bits in the upper triangular half,
	completing solving for these columns.
*/

void Codec::BackSubstituteAboveDiagonal()
{
	CAT_IF_DUMP(cout << endl << "---- BackSubstituteAboveDiagonal ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0; int heavyops = 0;)

	const int pivot_count = _defer_count + _mix_count;
	int pivot_i = pivot_count - 1;
	const u16 first_heavy_row = _defer_count + _dense_count;
	const u16 first_heavy_column = _first_heavy_column;

#if defined(CAT_WINDOWED_BACKSUB)
	// Build temporary storage space if windowing is to be used
	if (pivot_i >= CAT_ABOVE_WIN_THRESH_5)
	{
		// Calculate initial window size
		int w, next_check_i;
		if (pivot_i >= CAT_ABOVE_WIN_THRESH_7)
		{
			w = 7;
			next_check_i = CAT_ABOVE_WIN_THRESH_7;
		}
		else if (pivot_i >= CAT_ABOVE_WIN_THRESH_6)
		{
			w = 6;
			next_check_i = CAT_ABOVE_WIN_THRESH_6;
		}
		else if (pivot_i >= CAT_ABOVE_WIN_THRESH_5)
		{
			w = 5;
			next_check_i = CAT_ABOVE_WIN_THRESH_5;
		}
		else
		{
			w = 4;
			next_check_i = CAT_ABOVE_WIN_THRESH_4;
		}
		u32 win_lim = 1 << w;

		CAT_IF_DUMP(cout << "Activating windowed back-substitution with initial window " << w << endl;)

		// Use the first few peel column values as window table space
		// NOTE: The peeled column values were previously used up until this point,
		// but now they are unused, and so they can be reused for temporary space.
		u8 * CAT_RESTRICT win_table[128];
		PeelColumn * CAT_RESTRICT column = _peel_cols;
		u8 * CAT_RESTRICT column_src = _recovery_blocks;
		u32 jj = 1;
		for (u32 count = _block_count; count > 0; --count, ++column, column_src += _block_bytes)
		{
			// If column is peeled,
			if (column->mark == MARK_PEEL)
			{
				// Reuse the block value temporarily as window table space
				win_table[jj] = column_src;

				CAT_IF_DUMP(cout << "-- Window table entry " << jj << " set to column " << _block_count - count << endl;)

				// If done,
				if (++jj >= win_lim) break;
			}
		}

		CAT_IF_DUMP(if (jj < win_lim) cout << "!! Not enough space in peeled columns to generate a table.  Going back to normal back-substitute." << endl;)

		// If enough space was found,
		if (jj >= win_lim) for (;;)
		{
			// Calculate first column in window
			u16 backsub_i = pivot_i - w + 1;

			CAT_IF_DUMP(cout << "-- Windowing from " << backsub_i << " to " << pivot_i << " (inclusive)" << endl;)

			// Eliminate upper triangular part above windowed bits:

			// For each column,
			u64 ge_mask = (u64)1 << (pivot_i & 63);
			for (int src_pivot_i = pivot_i; src_pivot_i > backsub_i;
				--src_pivot_i, ge_mask = CAT_ROR64(ge_mask, 1))
			{
				u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * _ge_col_map[src_pivot_i];

				// If diagonal element is heavy,
				u16 ge_row_i = _pivots[src_pivot_i];
				if (ge_row_i >= first_heavy_row && src_pivot_i >= first_heavy_column)
				{
					// Look up row value
					u16 heavy_row_i = ge_row_i - first_heavy_row;
					u16 heavy_col_i = src_pivot_i - first_heavy_column;
					u8 code_value = _heavy_matrix[_heavy_pitch * heavy_row_i + heavy_col_i];

					// Normalize code value, setting it to 1 (implicitly nonzero)
					if (code_value != 1)
					{
						GF256MemDivide(src, code_value, _block_bytes);
						CAT_IF_ROWOP(++heavyops;)
					}

					CAT_IF_DUMP(cout << "Normalized diagonal for heavy pivot " << pivot_i << endl;)
				}

				CAT_IF_DUMP(cout << "Back-substituting small triangle from pivot " << src_pivot_i << "[" << (int)src[0] << "] :";)

				// For each row above the diagonal,
				u64 * CAT_RESTRICT ge_row = _ge_matrix + (src_pivot_i >> 6);
				for (int dest_pivot_i = backsub_i; dest_pivot_i < src_pivot_i; ++dest_pivot_i)
				{
					// If row is heavy,
					u16 dest_row_i = _pivots[dest_pivot_i];
					if (dest_row_i >= first_heavy_row && src_pivot_i >= first_heavy_column)
					{
						// If column is zero,
						u16 heavy_row_i = dest_row_i - first_heavy_row;
						u16 heavy_col_i = src_pivot_i - first_heavy_column;
						u8 code_value = _heavy_matrix[_heavy_pitch * heavy_row_i + heavy_col_i];
						if (!code_value) continue; // Skip it

						// Back-substitute
						u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[dest_pivot_i];
						if (code_value != 1)
						{
							GF256MemMulAdd(dest, code_value, src, _block_bytes);
							CAT_IF_ROWOP(++heavyops;)

							CAT_IF_DUMP(cout << " h" << dest_pivot_i;)
						}
						else
						{
							memxor(dest, src, _block_bytes);
							CAT_IF_ROWOP(++rowops;)

							CAT_IF_DUMP(cout << " *" << dest_pivot_i;)
						}
					}
					else
					{
						// If bit is set in that row,
						if (ge_row[_ge_pitch * dest_row_i] & ge_mask)
						{
							// Back-substitute
							u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[dest_pivot_i];
							memxor(dest, src, _block_bytes);
							CAT_IF_ROWOP(++rowops;)

							CAT_IF_DUMP(cout << " " << dest_pivot_i;)
						}
					}
				} // next pivot above

				CAT_IF_DUMP(cout << endl;)
			} // next pivot

			// Normalize the final diagonal element
			u16 ge_row_i = _pivots[backsub_i];
			if (ge_row_i >= first_heavy_row && backsub_i >= first_heavy_column)
			{
				// Look up row value
				u16 heavy_row_i = ge_row_i - first_heavy_row;
				u16 heavy_col_i = backsub_i - first_heavy_column;
				u8 code_value = _heavy_matrix[_heavy_pitch * heavy_row_i + heavy_col_i];

				// Divide by this code value (implicitly nonzero)
				if (code_value != 1)
				{
					u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i];
					GF256MemDivide(src, code_value, _block_bytes);
					CAT_IF_ROWOP(++heavyops;)
				}
			}

			CAT_IF_DUMP(cout << "-- Generating window table with " << w << " bits" << endl;)

			// Generate window table: 2 bits
			win_table[1] = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i];
			win_table[2] = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i + 1];
			memxor_set(win_table[3], win_table[1], win_table[2], _block_bytes);
			CAT_IF_ROWOP(++rowops;)

			// Generate window table: 3 bits
			win_table[4] = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i + 2];
			memxor_set(win_table[5], win_table[1], win_table[4], _block_bytes);
			memxor_set(win_table[6], win_table[2], win_table[4], _block_bytes);
			memxor_set(win_table[7], win_table[1], win_table[6], _block_bytes);
			CAT_IF_ROWOP(rowops += 3;)

			// Generate window table: 4 bits
			win_table[8] = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i + 3];
			for (int ii = 1; ii < 8; ++ii)
				memxor_set(win_table[8 + ii], win_table[ii], win_table[8], _block_bytes);
			CAT_IF_ROWOP(rowops += 7;)

			// Generate window table: 5+ bits
			if (w >= 5)
			{
				win_table[16] = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i + 4];
				for (int ii = 1; ii < 16; ++ii)
					memxor_set(win_table[16 + ii], win_table[ii], win_table[16], _block_bytes);
				CAT_IF_ROWOP(rowops += 15;)

				if (w >= 6)
				{
					win_table[32] = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i + 5];
					for (int ii = 1; ii < 32; ++ii)
						memxor_set(win_table[32 + ii], win_table[ii], win_table[32], _block_bytes);
					CAT_IF_ROWOP(rowops += 31;)

					if (w >= 7)
					{
						win_table[64] = _recovery_blocks + _block_bytes * _ge_col_map[backsub_i + 6];
						for (int ii = 1; ii < 64; ++ii)
							memxor_set(win_table[64 + ii], win_table[ii], win_table[64], _block_bytes);
						CAT_IF_ROWOP(rowops += 63;)
					}
				}
			}

			// If a row above the window may be heavy,
			if (pivot_i >= first_heavy_column)
			{
				// For each pivot in the window,
				u16 * CAT_RESTRICT pivot_row = _pivots;
				for (u16 ge_above_i = 0; ge_above_i < backsub_i; ++ge_above_i)
				{
					// If row is not heavy,
					u16 ge_row_i = *pivot_row++;
					if (ge_row_i < first_heavy_row)
						continue; // Skip it

					u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[ge_above_i];

					// If the first column of window is not heavy,
					u16 ge_column_j = backsub_i;
					if (ge_column_j < first_heavy_column)
					{
						// For each non-heavy column in the extra row,
						u64 ge_mask = (u64)1 << (ge_column_j & 63);
						u64 * CAT_RESTRICT ge_row = _ge_matrix + _ge_pitch * ge_row_i;
						for (; ge_column_j < first_heavy_column && ge_column_j <= pivot_i; ++ge_column_j, ge_mask = CAT_ROL64(ge_mask, 1))
						{
							// If column is non-zero,
							if (ge_row[ge_column_j >> 6] & ge_mask)
							{
								const u8 *src = _recovery_blocks + _block_bytes * _ge_col_map[ge_column_j];
								memxor(dest, src, _block_bytes);
								CAT_IF_ROWOP(++rowops;)
							}
						}
					}

					// For each heavy column,
					u16 heavy_row_i = ge_row_i - first_heavy_row;
					u16 heavy_col_j = ge_column_j - first_heavy_column;
					u8 * CAT_RESTRICT heavy_row = &_heavy_matrix[_heavy_pitch * heavy_row_i + heavy_col_j];
					for (; ge_column_j <= pivot_i; ++ge_column_j)
					{
						// If zero,
						u8 code_value = *heavy_row++;
						if (!code_value) continue; // Skip it

						// Back-substitute
						const u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * _ge_col_map[ge_column_j];
						if (code_value != 1)
						{
							GF256MemMulAdd(dest, code_value, src, _block_bytes);
							CAT_IF_ROWOP(++heavyops;)
						}
						else
						{
							memxor(dest, src, _block_bytes);
							CAT_IF_ROWOP(++rowops;)
						}
					} // next column in row
				} // next pivot in window
			} // end if contains heavy

			// Only add window table entries for rows under this limit
			u16 window_row_limit = (pivot_i >= first_heavy_column) ? first_heavy_row : 0x7fff;

			// If not straddling words,
			u32 first_word = backsub_i >> 6;
			u32 shift0 = backsub_i & 63;
			u32 last_word = pivot_i >> 6;
			u16 * CAT_RESTRICT pivot_row = _pivots;
			if (first_word == last_word)
			{
				// For each pivot row,
				for (u16 above_pivot_i = 0; above_pivot_i < backsub_i; ++above_pivot_i)
				{
					// If pivot row is heavy,
					u16 ge_row_i = *pivot_row++;
					if (ge_row_i >= window_row_limit)
						continue; // Skip it

					// Calculate window bits
					u64 * CAT_RESTRICT ge_row = _ge_matrix + first_word + _ge_pitch * ge_row_i;
					u32 win_bits = (u32)(ge_row[0] >> shift0) & (win_lim - 1);

					// If any XOR needs to be performed,
					if (win_bits != 0)
					{
						CAT_IF_DUMP(cout << "Adding window table " << win_bits << " to pivot " << above_pivot_i << endl;)

						// Back-substitute
						u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[above_pivot_i];
						memxor(dest, win_table[win_bits], _block_bytes);
						CAT_IF_ROWOP(++rowops;)
					}
				}
			}
			else // Rare: Straddling case
			{
				u32 shift1 = 64 - shift0;

				// For each pivot row,
				for (u16 above_pivot_i = 0; above_pivot_i < backsub_i; ++above_pivot_i)
				{
					// If pivot row is heavy,
					u16 ge_row_i = *pivot_row++;
					if (ge_row_i >= window_row_limit)
						continue; // Skip it

					// Calculate window bits
					u64 * CAT_RESTRICT ge_row = _ge_matrix + first_word + _ge_pitch * ge_row_i;
					u32 win_bits = ( (u32)(ge_row[0] >> shift0) | (u32)(ge_row[1] << shift1) ) & (win_lim - 1);

					// If any XOR needs to be performed,
					if (win_bits != 0)
					{
						CAT_IF_DUMP(cout << "Adding window table " << win_bits << " to pivot " << above_pivot_i << endl;)

						// Back-substitute
						u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[above_pivot_i];
						memxor(dest, win_table[win_bits], _block_bytes);
						CAT_IF_ROWOP(++rowops;)
					}
				}
			}

			// If column index falls below window size,
			pivot_i -= w;
			if (pivot_i < next_check_i)
			{
				if (pivot_i >= CAT_ABOVE_WIN_THRESH_6)
				{
					w = 6;
					next_check_i = CAT_ABOVE_WIN_THRESH_6;
				}
				else if (pivot_i >= CAT_ABOVE_WIN_THRESH_5)
				{
					w = 5;
					next_check_i = CAT_ABOVE_WIN_THRESH_5;
				}
				else if (pivot_i >= CAT_ABOVE_WIN_THRESH_4)
				{
					w = 4;
					next_check_i = CAT_ABOVE_WIN_THRESH_4;
				}
				else break;

				// Update window limit
				win_lim = 1 << w;
			}
		} // next window
	} // end if windowed
#endif // CAT_WINDOWED_BACKSUB

	// For each remaining pivot,
	u64 ge_mask = (u64)1 << (pivot_i & 63);
	for (; pivot_i >= 0; --pivot_i, ge_mask = CAT_ROR64(ge_mask, 1))
	{
		// Calculate source
		u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * _ge_col_map[pivot_i];

		// If diagonal element is heavy,
		u16 ge_row_i = _pivots[pivot_i];
		if (ge_row_i >= first_heavy_row && pivot_i >= first_heavy_column)
		{
			// Look up row value
			u16 heavy_row_i = ge_row_i - first_heavy_row;
			u16 heavy_col_i = pivot_i - first_heavy_column;
			u8 code_value = _heavy_matrix[_heavy_pitch * heavy_row_i + heavy_col_i];

			// Normalize code value, setting it to 1 (implicitly nonzero)
			if (code_value != 1)
			{
				GF256MemDivide(src, code_value, _block_bytes);
				CAT_IF_ROWOP(++heavyops;)
			}

			CAT_IF_DUMP(cout << "Normalized diagonal for heavy pivot " << pivot_i << endl;)
		}

		CAT_IF_DUMP(cout << "Pivot " << pivot_i << "[" << (int)src[0] << "]:";)

		// For each pivot row above it,
		u64 * CAT_RESTRICT ge_row = _ge_matrix + (pivot_i >> 6);
		for (int ge_up_i = 0; ge_up_i < pivot_i; ++ge_up_i)
		{
			// If element is heavy,
			u16 up_row_i = _pivots[ge_up_i];
			if (up_row_i >= first_heavy_row && ge_up_i >= first_heavy_column)
			{
				// If column is zero,
				u16 heavy_row_i = up_row_i - first_heavy_row;
				u16 heavy_col_i = pivot_i - first_heavy_column;
				u8 code_value = _heavy_matrix[_heavy_pitch * heavy_row_i + heavy_col_i];
				if (!code_value) continue; // Skip it

				// Back-substitute
				u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * _ge_col_map[ge_up_i];
				if (code_value != 1)
				{
					GF256MemMulAdd(dest, code_value, src, _block_bytes);
					CAT_IF_ROWOP(++heavyops;)

					CAT_IF_DUMP(cout << " h" << up_row_i;)
				}
				else
				{
					memxor(dest, src, _block_bytes);
					CAT_IF_ROWOP(++rowops;)

					CAT_IF_DUMP(cout << " *" << up_row_i;)
				}
			}
			else
			{
				// If bit is set in that row,
				if (ge_row[_ge_pitch * up_row_i] & ge_mask)
				{
					// Back-substitute
					u8 *dest = _recovery_blocks + _block_bytes * _ge_col_map[ge_up_i];
					memxor(dest, src, _block_bytes);
					CAT_IF_ROWOP(++rowops;)

					CAT_IF_DUMP(cout << " " << up_row_i;)
				}
			}
		} // next pivot above

		CAT_IF_DUMP(cout << endl;)
	}

	CAT_IF_ROWOP(cout << "BackSubstituteAboveDiagonal used " << rowops << " row ops = " << rowops / (double)_block_count << "*N and " << heavyops << " heavy ops" << endl;)
}

/*
	Substitute

		This function generates all of the remaining column values.
	At the point this function is called, the GE columns were solved
	by BackSubstituteAboveDiagonal().  The remaining column values
	are solved by regenerating rows in forward order of peeling.

		Note that as each row is generated, that all of the columns
	active in that row have already been solved.  So long as the
	substitution follows in forward solution order, this is guaranteed.

		It is also possible to start from the Compression matrix values
	and substitute into that matrix.  However, because the mixing columns
	are so dense, it is actually faster in every case to just regenerate
	the rows from scratch and throw away those results.
*/

void Codec::Substitute()
{
	CAT_IF_DUMP(cout << endl << "---- Substitute ----" << endl << endl;)

	CAT_IF_ROWOP(u32 rowops = 0;)

	// For each column that has been peeled,
	PeelRow * CAT_RESTRICT row;
	for (u16 row_i = _peel_head_rows; row_i != LIST_TERM; row_i = row->next)
	{
		row = &_peel_rows[row_i];
		u16 dest_column_i = row->peel_column;
		u8 * CAT_RESTRICT dest = _recovery_blocks + _block_bytes * dest_column_i;

		CAT_IF_DUMP(cout << "Generating column " << dest_column_i << ":";)

		const u8 * CAT_RESTRICT input_src = _input_blocks + _block_bytes * row_i;
		CAT_IF_DUMP(cout << " " << row_i << ":[" << (int)input_src[0] << "]";)

		// Set up mixing column generator
		u16 mix_a = row->mix_a;
		u16 mix_x = row->mix_x0;
		const u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * (_block_count + mix_x);

		// If copying from final block,
		if (row_i != _block_count - 1)
			memxor_set(dest, src, input_src, _block_bytes);
		else
		{
			memxor_set(dest, src, input_src, _input_final_bytes);
			memcpy(dest + _input_final_bytes, src + _input_final_bytes, _block_bytes - _input_final_bytes);
		}
		CAT_IF_ROWOP(++rowops;)

		// Add next two mixing columns in
		IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
		const u8 * CAT_RESTRICT src0 = _recovery_blocks + _block_bytes * (_block_count + mix_x);
		IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
		const u8 * CAT_RESTRICT src1 = _recovery_blocks + _block_bytes * (_block_count + mix_x);
		memxor_add(dest, src0, src1, _block_bytes);
		CAT_IF_ROWOP(++rowops;)

		// If at least two peeling columns are set,
		u16 weight = row->peel_weight;
		if (weight >= 2) // common case:
		{
			u16 a = row->peel_a;
			u16 column0 = row->peel_x0;
			--weight;

			u16 column_i = column0;
			IterateNextColumn(column_i, _block_count, _block_next_prime, a);

			// Common case:
			if (column0 != dest_column_i)
			{
				const u8 * CAT_RESTRICT peel0 = _recovery_blocks + _block_bytes * column0;

				// Common case:
				if (column_i != dest_column_i)
					memxor_add(dest, peel0, _recovery_blocks + _block_bytes * column_i, _block_bytes);
				else // rare:
					memxor(dest, peel0, _block_bytes);
			}
			else // rare:
				memxor(dest, _recovery_blocks + _block_bytes * column_i, _block_bytes);
			CAT_IF_ROWOP(++rowops;)

			// For each remaining column,
			while (--weight > 0)
			{
				IterateNextColumn(column_i, _block_count, _block_next_prime, a);
				const u8 * CAT_RESTRICT src = _recovery_blocks + _block_bytes * column_i;

				CAT_IF_DUMP(cout << " " << column_i;)

				// If column is not the solved one,
				if (column_i != dest_column_i)
				{
					memxor(dest, src, _block_bytes);
					CAT_IF_ROWOP(++rowops;)
					CAT_IF_DUMP(cout << "[" << (int)src[0] << "]";)
				}
				else
				{
					CAT_IF_DUMP(cout << "*";)
				}
			}
		} // end if weight 2

		CAT_IF_DUMP(cout << endl;)
	}

	CAT_IF_ROWOP(cout << "Substitute used " << rowops << " row ops = " << rowops / (double)_block_count << "*N" << endl;)
}


//// Main Driver

/*
	Each element of the DENSE_SEEDS table represents the best
	Abyssinian PRNG seed to use to generate a Shuffle-2 Code
	for the dense rows that has "good" properties.

	Element 0 is for D = 14,
	Element 1 is for D = 18,
	Element 2 is for D = 22,
	and so on.

	If the dense matrix size is DxD, then when D Mod 4 = 0,
	a Shuffle-2 Code matrix is not invertible.  And when
	D Mod 4 = 2, then it is better than the random matrix
	for small D.  However for D > 22, even the best sizes
	are not as often invertible as the random matrix.

	However this is not a huge issue.  Skipping over a lot
	of details out of scope of this article, the way that
	the Shuffle-2 Codes are used is as follows:

	How the codes are actually used:

		(1) Several DxD random matrices are produced with
		Shuffle-2 Codes.  Call these matrices {R0, R1, R2, R3...}.

		(2) Only a few columns (M) are selected at random from
		these matrices to form a new matrix, called the GE matrix.
		M is slightly larger than the square root of the total
		number of columns across all random R# matrices.

		(3) The resulting matrix must be rank M or it leads to
		lower error correcting performance in Wirehair.

		(4) Furthermore, due to the other things going on around
		this algorithm, the bits in the GE matrix are somewhat
		randomly flipped after the first third of them.
		This makes it easier for the matrix to be full rank, but
		is also a challenge because it means that if the columns
		have low Hamming weight that they are in danger of all
		being flipped off.

		(5) I want to be able to use the same random-looking R#
		matrices for any given D and I want it to behave well.
		This allows me to use a short table of PRNG seeds for
		each value of D to generate a best-performing set of
		R# matrices.

	From the way the Shuffle-2 Code is used,
	some requirements are apparent:

		(1) Should be able to generate the R# matrices from a seed.

		(2) The average rank of randomly-selected columns should
		be high to satisfy the primary goal of the code.

		(2a) To achieve (2), the average Hamming distance between
		columns should be maximized.

		(2b) The minimum Hamming distance between columns is 2.

		(2c) Based on the empirical data from before, D is chosen
		so that D Mod 4 = 2.
		In practice D will be rounded up to the next "good" one.

		(3) The minimum Hamming weight of each column should be 3
		to avoid being flipped into oblivion.

	To find the best matrices, I tried all seeds from 0..65535
	(time permitting) and generated D of the DxD matrices.
	Then, I verified requirement (2b) and (3).  Of the remaining
	options, the one with the highest average rank was chosen.
	I am only interested in values of D between 14 and 486 for
	practical use.  Smaller D are special cases that do not need
	to be in the table, and larger D are unused.

	I found that (2b) is not possible to satisfy.
	The minimum Hamming distance will always be 1.

	Here's some example data for 22x22:

	Seed 4504 minimum Hamming distance of 1 and
	average = 14.4 and minimum Hamming weight of 13
		Rank 2 at 0.999
		Rank 4 at 0.986
		Rank 6 at 0.96
		Rank 8 at 0.92
		Rank 10 at 0.887
		Rank 12 at 0.84
		Rank 14 at 0.804
		Rank 16 at 0.722
		Rank 18 at 0.648
		Rank 20 at 0.504
		Rank 21 at 0.36
		Rank 22 at 0.174

	I revised the seed search to look at the best two seeds
	in terms of average Hamming distance, and picked the one
	that is more often invertible for rank D-4 through D-2.
	This comes from the fact that often times the average
	Hamming distance being higher doesn't always mean it is
	better.  The real test is how often it is full rank.

	I set a 4 minute timeout for the best seed search and let
	it run overnight.  The result is a small 118-element table
	that determines all of the unchanging dense rows in the
	Wirehair check matrix for N=2 up to N=64000, providing
	best performance for a given number of dense rows.

	Some other random thoughts:

	+ The average rank of randomly selected columns drops off
	pretty sharply near D.  To achieve 90% average invertibility,
	the number of rows needs to roughly double and add one.
	Adding just one row puts it above normal random matrices for
	invertibility rate.

	+ The seeds are not necessarily the best that could be found,
	since for each D, a range of N use that value of D, and this
	is much less than D*D -- it is up to 64,000 tops.
*/

static const u16 DENSE_SEEDS[119] = {
	4181, 26667, 4504, 11009, 3438, 14320, 15822, 50870,
	4234, 1376, 30232, 1177, 8576, 3099, 8178, 52837,
	773, 5032, 10746, 11964, 1005, 1568, 12581, 2820,
	289, 2, 4322, 4097, 481, 1383, 3765, 166,
	3286, 2605, 3101, 851, 465, 1127, 1548, 1771,
	793, 1170, 361, 1151, 27, 159, 460, 14,
	267, 478, 109, 70, 279, 427, 17, 39,
	20, 5, 34, 15, 22, 37, 24, 23,
	18, 0, 30, 25, 4, 19, 9, 13,
	16, 2, 3, 21, 4, 1, 161,
	29, 127, 30, 21, 30, 24, 86, 37,
	6, 43, 0, 48, 35, 12, 16, 1,
	82, 94, 25, 64, 15, 27, 58, 70,
	2, 26, 15, 31, 27, 7, 53, 56,
	30, 54, 18, 79, 31, 5, 41, 12
};

/*
	These tables were lovingly hand-crafted by hard-working indigenous peoples:
*/

const int SMALL_SEED_MAX = 261;
static const u16 SMALL_PEEL_SEEDS[262] = {
	0, 0, 6, 2, 116, 275, 593, 620, 431, 539, 134, 103, 157, 410, 33, 198, 94, 116,
	207, 227, 34, 34, 2, 174, 23, 198, 159, 97, 265, 89, 31, 41, 113, 89, 126, 29,
	70, 33, 56, 140, 163, 109, 124, 161, 135, 163, 19, 6, 158, 27, 107, 22, 122, 129,
	142, 27, 8, 125, 0, 63, 108, 16, 104, 114, 40, 32, 105, 122, 63, 54, 29, 98, 95,
	40, 14, 12, 60, 17, 79, 72, 95, 78, 14, 88, 0, 23, 95, 42, 14, 73, 1, 33, 10, 17,
	80, 26, 8, 16, 2, 66, 17, 80, 30, 69, 4, 5, 29, 12, 71, 38, 14, 55, 22, 72, 2,
	43, 67, 41, 44, 6, 37, 1, 50, 32, 44, 38, 29, 20, 48, 58, 38, 52, 27, 59, 27, 38,
	42, 27, 43, 38, 36, 0, 15, 63, 57, 11, 23, 41, 36, 57, 18, 59, 2, 11, 34, 8, 28,
	0, 9, 42, 26, 3, 55, 6, 55, 22, 18, 17, 8, 29, 31, 43, 29, 20, 25, 15, 23, 31, 0,
	6, 0, 33, 47, 49, 37, 2, 29, 41, 33, 27, 22, 39, 25, 6, 29, 24, 10, 45, 18, 45, 19,
	17, 3, 30, 3, 18, 8, 44, 43, 4, 30, 38, 28, 2, 40, 26, 19, 4, 37, 45, 22, 40, 6,
	1, 24, 7, 24, 38, 20, 38, 1, 17, 22, 38, 5, 6, 30, 32, 0, 2, 39, 32, 18, 38, 3, 4,
	2, 4, 39, 6, 22, 7, 12, 6, 14, 0, 5, 12, 15, 5, 19, 1
};

// 8KB bitfield table for seeds that cause the encoder to choke
static const u64 EXCEPT_SEEDS[1000] = {
0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x400200000000000ULL, 0x200000000000000ULL, 0x9000ULL, 0x440000000000004ULL,
0x100100000000ULL, 0x100000000ULL, 0x2ULL, 0x8020000ULL, 0x404080000200000ULL, 0x200ULL, 0x82000000000000ULL, 0x800800000200ULL,
0x1000000000000ULL, 0x400000ULL, 0x402ULL, 0x0ULL, 0x1000000000000ULL, 0x0ULL, 0x1000000ULL, 0x10000000002ULL,
0x0ULL, 0x208000000200001ULL, 0x20000040000000ULL, 0x10000000000ULL, 0x0ULL, 0x2000000000000ULL, 0x400800000000ULL, 0xa00000000000200ULL,
0x40080000000ULL, 0x800000040000000ULL, 0x104000ULL, 0x1000000020ULL, 0x410000000000000ULL, 0x10000ULL, 0x4000000ULL, 0x1000000000000000ULL,
0x4000ULL, 0x60000008ULL, 0x40ULL, 0x2000040000ULL, 0x0ULL, 0x0ULL, 0x200ULL, 0x100000000000ULL,
0x4000000000000ULL, 0x100500ULL, 0x0ULL, 0x200840000010000ULL, 0x800ULL, 0x400000000ULL, 0x0ULL, 0x800020ULL,
0x1040000000000000ULL, 0x0ULL, 0x0ULL, 0x8000000080000000ULL, 0x4ULL, 0x2000020000000ULL, 0x0ULL, 0x100000ULL,
0x0ULL, 0x80ULL, 0x10040000ULL, 0x80000000000ULL, 0xc0000004000ULL, 0x0ULL, 0x8000000008000ULL, 0x1800000080000000ULL,
0x0ULL, 0x10ULL, 0x10000000000ULL, 0x80ULL, 0x0ULL, 0x0ULL, 0xa0000ULL, 0x40000ULL,
0x200100ULL, 0x2000800000000000ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x5ULL, 0x400ULL, 0x800000080004006ULL,
0x204000120000ULL, 0x400008002080ULL, 0x100ULL, 0x0ULL, 0x101000ULL, 0x80000ULL, 0x1000004000000ULL, 0x2ULL,
0x1000000ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x20000200008ULL, 0x2000000000ULL, 0x20000020000000ULL, 0x0ULL,
0x100ULL, 0x0ULL, 0x8000000000000ULL, 0x4000200000ULL, 0x0ULL, 0x800000080001200ULL, 0x4000000000000ULL, 0x1000000000000ULL,
0x400000000400ULL, 0x8ULL, 0x10ULL, 0x810400000808ULL, 0x0ULL, 0x800000000900ULL, 0x400000000000ULL, 0x8000800000000000ULL,
0x10000000000800ULL, 0x40000ULL, 0x2000000ULL, 0x0ULL, 0x1000000000ULL, 0x0ULL, 0x208000000ULL, 0x0ULL,
0x1ULL, 0x1000000400000000ULL, 0x4006000000000010ULL, 0x200000000000000ULL, 0x8000ULL, 0x0ULL, 0x100000ULL, 0x80820000000000ULL,
0x0ULL, 0x100004000000ULL, 0x10000000000000ULL, 0x1000000ULL, 0x8000000ULL, 0x200ULL, 0x4ULL, 0x80ULL,
0x2000020000002000ULL, 0x0ULL, 0x200000000000000ULL, 0x40000000000004ULL, 0x800000000000ULL, 0x0ULL, 0x0ULL, 0x1000000000000ULL,
0x800000000800000ULL, 0x40020000000ULL, 0xc0000000040ULL, 0x4000002ULL, 0x0ULL, 0x40010000000ULL, 0x0ULL, 0x10000000000000ULL,
0x2000000ULL, 0x0ULL, 0x0ULL, 0x100002010000000ULL, 0x102000ULL, 0x0ULL, 0x0ULL, 0x80000000ULL,
0x60002040800ULL, 0x0ULL, 0x0ULL, 0x4000000000ULL, 0x8000000000000020ULL, 0x100000000000ULL, 0x40ULL, 0x4000000080000000ULL,
0x0ULL, 0x4000000000000ULL, 0x100ULL, 0x200000004000200ULL, 0x200000000900000ULL, 0x0ULL, 0x2100000000000ULL, 0x10000ULL,
0x24010020000000ULL, 0x40000000000002ULL, 0x8000000000000000ULL, 0x0ULL, 0x0ULL, 0x4000ULL, 0x20ULL, 0x420000000000000ULL,
0x0ULL, 0x0ULL, 0x480000000ULL, 0x200000000002000ULL, 0x210000000ULL, 0x80000000001ULL, 0x10000000000000ULL, 0x20000010ULL,
0x0ULL, 0x0ULL, 0x80008000000ULL, 0x400000000200000ULL, 0x100000000ULL, 0x0ULL, 0x101000080ULL, 0x20020000ULL,
0x800000000ULL, 0x0ULL, 0x2ULL, 0x800000ULL, 0x1000000000ULL, 0x0ULL, 0x400000000ULL, 0x0ULL,
0x200000000000000ULL, 0x2000000000000000ULL, 0x0ULL, 0x20000000000000ULL, 0x400004ULL, 0x40000000ULL, 0x2000000004ULL, 0x80000000000ULL,
0x0ULL, 0x2000000000200ULL, 0x401000000000000ULL, 0x2200002800000000ULL, 0x4000000ULL, 0x40000002000000ULL, 0x10020000ULL, 0x40000160ULL,
0x400000000020ULL, 0x80000ULL, 0x80000000000ULL, 0x0ULL, 0x0ULL, 0x1ULL, 0x0ULL, 0x200000000200ULL,
0x400000000000ULL, 0x82000ULL, 0x0ULL, 0x80000000ULL, 0x80010ULL, 0x200000ULL, 0x800000040000000ULL, 0x0ULL,
0x100004000000102ULL, 0x1000100000000ULL, 0x4000000000ULL, 0x0ULL, 0x400ULL, 0x1000000001000ULL, 0x800000000000000ULL, 0x20000000ULL,
0x1000000000000ULL, 0x0ULL, 0x800000000ULL, 0x0ULL, 0x40000040ULL, 0x4000000000000ULL, 0x8000ULL, 0x4000800000000400ULL,
0x1100000000ULL, 0x1000000000ULL, 0x100000000000000ULL, 0x0ULL, 0x2000000008ULL, 0x2002000000ULL, 0x0ULL, 0x8000000ULL,
0x1200ULL, 0x8008000000000000ULL, 0x0ULL, 0x4000100ULL, 0x102000000000ULL, 0x0ULL, 0x1010000000ULL, 0x20000ULL,
0x0ULL, 0x8000000000ULL, 0x0ULL, 0x20000040000000ULL, 0x20000000000001ULL, 0x2000000000020000ULL, 0x200000000000000ULL, 0x100010400100840ULL,
0x2000000000200ULL, 0x1000000000400ULL, 0x200000000000000ULL, 0xa004000000002ULL, 0x4000002004002ULL, 0x0ULL, 0x8000ULL, 0x1000800000808800ULL,
0x0ULL, 0x0ULL, 0x80000000004ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x4000000000000ULL, 0x1000000200ULL,
0x2000000000040ULL, 0x40ULL, 0x4000ULL, 0x100ULL, 0x0ULL, 0x10ULL, 0x0ULL, 0x80000000ULL,
0x40020000000ULL, 0x800000004100000ULL, 0x0ULL, 0x1ULL, 0x40000ULL, 0x400000000ULL, 0x0ULL, 0x0ULL,
0x80ULL, 0x40ULL, 0xa00000ULL, 0x408000ULL, 0x0ULL, 0x0ULL, 0x1004000008004ULL, 0x4000000000ULL,
0x2000000ULL, 0x4000ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x20000ULL, 0x100000010000400ULL,
0x8200000ULL, 0x801000000000ULL, 0x0ULL, 0x0ULL, 0x4000000ULL, 0x4020000ULL, 0x20000000000ULL, 0x1200800200000000ULL,
0x0ULL, 0x200000000080ULL, 0x0ULL, 0x8000088000000000ULL, 0x20000000ULL, 0x0ULL, 0x0ULL, 0x3ULL,
0x8000000000100ULL, 0x80000000100ULL, 0x800000200004ULL, 0x8008000000ULL, 0x0ULL, 0x2000200100000ULL, 0x10000000000ULL, 0x0ULL,
0x802400ULL, 0x400010000000ULL, 0x10000000ULL, 0x1ULL, 0x4011ULL, 0x0ULL, 0x400000080000ULL, 0x8100000ULL,
0x0ULL, 0x2060000000000000ULL, 0x100000000000000ULL, 0x0ULL, 0x200000000ULL, 0x600000080000ULL, 0x400080ULL, 0x10808000ULL,
0x2ULL, 0x10000000000ULL, 0x0ULL, 0x2000ULL, 0x10ULL, 0x80001000000ULL, 0x200000000ULL, 0x100ULL,
0x40020800000ULL, 0x8100000020ULL, 0x200ULL, 0x2001000ULL, 0x0ULL, 0x0ULL, 0x400008ULL, 0xa000040400080000ULL,
0x0ULL, 0x20400040000008ULL, 0x2000000000000000ULL, 0x0ULL, 0x100000000000ULL, 0x400000000000000ULL, 0x1000000000000000ULL, 0x0ULL,
0x0ULL, 0x0ULL, 0x0ULL, 0x20000010000000ULL, 0x0ULL, 0x400000000000ULL, 0x0ULL, 0x0ULL,
0x80000000000000ULL, 0x0ULL, 0x0ULL, 0x40ULL, 0x0ULL, 0x1000000400000000ULL, 0x10800000000008ULL, 0x1000800000000ULL,
0x100000001000000ULL, 0x80ULL, 0x0ULL, 0x4040000020000ULL, 0x10000ULL, 0x0ULL, 0x0ULL, 0x200000000000ULL,
0x108000ULL, 0x800000081000001ULL, 0x8000200000000ULL, 0x4000000000000000ULL, 0x0ULL, 0x4000000000ULL, 0x460000000002000ULL, 0x8000000000000000ULL,
0x0ULL, 0x2000000002ULL, 0x400000000ULL, 0x10001ULL, 0x0ULL, 0x800000200000000ULL, 0x4ULL, 0x80000000010100ULL,
0x4000000000000000ULL, 0x0ULL, 0x1000000ULL, 0x140000000000000ULL, 0x800000ULL, 0x8000000ULL, 0x0ULL, 0x1000000000000000ULL,
0x0ULL, 0x1000000ULL, 0x2000ULL, 0x10000000000000ULL, 0x0ULL, 0x0ULL, 0x40000ULL, 0x80000001ULL,
0xc020000000800ULL, 0x81000ULL, 0x804000000000000ULL, 0x800ULL, 0x4010000000080ULL, 0x18000000000ULL, 0x0ULL, 0x8000000020000000ULL,
0x20000000000ULL, 0x400000810002000ULL, 0x40000000200ULL, 0x0ULL, 0x4000002000000040ULL, 0x0ULL, 0x401004000000002ULL, 0x400000100001ULL,
0x100000000ULL, 0x10000000000ULL, 0x800ULL, 0x10000000ULL, 0x0ULL, 0x8000000001800000ULL, 0x300000108000ULL, 0x2000000000ULL,
0x200000000021ULL, 0x280ULL, 0x80100002000000ULL, 0x8000080000ULL, 0x80000000ULL, 0x0ULL, 0x8000000000000020ULL, 0x0ULL,
0x40000000000ULL, 0x10000000000ULL, 0x20000000ULL, 0x0ULL, 0x0ULL, 0x400000000000ULL, 0x4000000000000ULL, 0x2000000000000ULL,
0x200000ULL, 0x0ULL, 0x4000001000000ULL, 0x20040000ULL, 0x4020000000000ULL, 0x0ULL, 0x8000020000000000ULL, 0x20000000ULL,
0x100000ULL, 0x1ULL, 0x800000000ULL, 0x4008000000ULL, 0x0ULL, 0x0ULL, 0x10000000000000ULL, 0x2800ULL,
0x1000000000000020ULL, 0x0ULL, 0x40000ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x4000ULL, 0x80000009ULL,
0x0ULL, 0x608000100ULL, 0x8000000000000000ULL, 0x1000000002000820ULL, 0x0ULL, 0x0ULL, 0x1000ULL, 0x800ULL,
0x20000020200000ULL, 0x400000000ULL, 0x0ULL, 0x10000000000000ULL, 0x0ULL, 0x20000ULL, 0x0ULL, 0x0ULL,
0x40000002000ULL, 0x2000000000ULL, 0x400000ULL, 0x10000000000000ULL, 0x400000000ULL, 0x400ULL, 0x2000ULL, 0x200000000020ULL,
0x0ULL, 0x0ULL, 0x4000ULL, 0x8000002000ULL, 0x280000000000ULL, 0x2000000000ULL, 0x100ULL, 0x400000000ULL,
0x8000000020000ULL, 0x0ULL, 0x100000ULL, 0x0ULL, 0x2000000000ULL, 0x1210ULL, 0x40000000010000ULL, 0x20ULL,
0x8000000002000ULL, 0x1200ULL, 0x80000000000008ULL, 0xc0000ULL, 0x0ULL, 0x0ULL, 0x2000020000ULL, 0x20000000ULL,
0x0ULL, 0x0ULL, 0x80040ULL, 0x0ULL, 0x8000000280001ULL, 0x0ULL, 0x8000002ULL, 0x0ULL,
0x400000000000ULL, 0x2000000000000000ULL, 0x1000020000ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x1000000000000000ULL, 0x800010000000000ULL,
0x0ULL, 0x0ULL, 0x20ULL, 0x0ULL, 0x0ULL, 0x1000ULL, 0x1000000000020ULL, 0x0ULL,
0x4ULL, 0x8000000000ULL, 0x80081ULL, 0x40000000000020ULL, 0x0ULL, 0x2200000000000000ULL, 0x4040040000040ULL, 0x400000000000002ULL,
0x24ULL, 0x0ULL, 0x8000400000000ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x100000000000000ULL, 0x0ULL,
0x6000400000004ULL, 0x2000000000000000ULL, 0x0ULL, 0x1ULL, 0x440000000c00ULL, 0x8000000ULL, 0x400000ULL, 0x0ULL,
0x400000000ULL, 0x0ULL, 0x0ULL, 0x8000000000ULL, 0x2000000800000000ULL, 0x1040000000000000ULL, 0x100000000000ULL, 0x0ULL,
0x100000000ULL, 0x200ULL, 0x10010000000000ULL, 0x10000000ULL, 0x2000ULL, 0x200000080000ULL, 0x10200002000000ULL, 0x8ULL,
0x800000000ULL, 0x20000000000000ULL, 0x100000000000000ULL, 0x80400000004ULL, 0x0ULL, 0x200000000ULL, 0x1000ULL, 0x40000ULL,
0x140000008000000ULL, 0x0ULL, 0x100000000000ULL, 0x800000000ULL, 0x0ULL, 0x0ULL, 0x800000000000000ULL, 0x4000010000000ULL,
0x4000000000000ULL, 0x8000ULL, 0x8000000001200ULL, 0x100001000000200ULL, 0x1ULL, 0x100000000020ULL, 0x0ULL, 0x0ULL,
0x0ULL, 0x10000000002000ULL, 0x0ULL, 0x0ULL, 0x40001000000000ULL, 0x100000021ULL, 0x0ULL, 0x800100000100100ULL,
0x8000ULL, 0x80000000ULL, 0x40400000002ULL, 0x0ULL, 0x200000000040008ULL, 0x0ULL, 0x0ULL, 0x2000000000ULL,
0x20000000ULL, 0x80ULL, 0x10ULL, 0x0ULL, 0x0ULL, 0x80000000000ULL, 0x2040000000000000ULL, 0x40000000000000ULL,
0x0ULL, 0x20004000000000ULL, 0x0ULL, 0x800000004000000ULL, 0x8000000400000000ULL, 0x80100000000000ULL, 0x8000020000ULL, 0x100000000000020ULL,
0x8800000ULL, 0x1008000000040110ULL, 0x2000000000000000ULL, 0x0ULL, 0x100ULL, 0x10000000ULL, 0x1000004000ULL, 0x200000000800ULL,
0x4000000000000000ULL, 0x0ULL, 0x8000004002000ULL, 0x4080000000000000ULL, 0x1200000800000000ULL, 0x9000000200000000ULL, 0x0ULL, 0x0ULL,
0x0ULL, 0x40000000ULL, 0x0ULL, 0x400000000000ULL, 0x2006002000000ULL, 0x100000040000000ULL, 0x200020040000ULL, 0x2000000000ULL,
0x2000000000008ULL, 0x0ULL, 0x204040040ULL, 0x10080ULL, 0x2000000000ULL, 0x100000ULL, 0x400000000000ULL, 0x100000000000ULL,
0x500ULL, 0x200002ULL, 0x2000000000ULL, 0x0ULL, 0x10000000000000ULL, 0x200000ULL, 0x1000000000000000ULL, 0x100000000000ULL,
0x80000000000000ULL, 0x1000ULL, 0x0ULL, 0x40000000000ULL, 0x2000ULL, 0x200000000000ULL, 0x4000000000ULL, 0x10000000ULL,
0x0ULL, 0x1004000000000040ULL, 0x4000000000ULL, 0x200000aULL, 0x2000000ULL, 0x400040000200000ULL, 0x800001000100ULL, 0x0ULL,
0x20800000000000ULL, 0x4000100000100000ULL, 0x0ULL, 0x400ULL, 0x80000000000ULL, 0x40002000ULL, 0x4000000000082000ULL, 0x2000000000000ULL,
0x0ULL, 0x200ULL, 0x200000000000ULL, 0x1ULL, 0x1000000000400000ULL, 0x20008000ULL, 0x80000000000800ULL, 0x100ULL,
0x20000000000000ULL, 0x10000000000000ULL, 0x22000ULL, 0x40100ULL, 0x200000000000800ULL, 0x84000000000ULL, 0x1000000000000ULL, 0x0ULL,
0x2000000ULL, 0x0ULL, 0x0ULL, 0x100000000ULL, 0x100240000080ULL, 0x200010000ULL, 0x408000000000ULL, 0x1000000000000000ULL,
0x40000000b000000ULL, 0x0ULL, 0x0ULL, 0x20000000040ULL, 0x1ULL, 0x0ULL, 0x20000280000ULL, 0x20000ULL,
0x200000000ULL, 0x1000800000000ULL, 0x0ULL, 0x0ULL, 0x400ULL, 0x8ULL, 0x0ULL, 0x0ULL,
0x40000000ULL, 0x0ULL, 0x8400000000000ULL, 0x40180400ULL, 0x2000000ULL, 0x0ULL, 0x0ULL, 0x10000000000000ULL,
0x8000000ULL, 0x80800ULL, 0x10000000000000ULL, 0x2000000c00ULL, 0x1000000000800ULL, 0x0ULL, 0x10000000ULL, 0x0ULL,
0x400000000000010ULL, 0x20020008030000ULL, 0x100000002008ULL, 0x20000000ULL, 0x0ULL, 0x8000000000010ULL, 0x8000000000ULL, 0x1000ULL,
0x4000000000000ULL, 0x0ULL, 0x4000000ULL, 0x11000000000ULL, 0x100000000ULL, 0x20ULL, 0x0ULL, 0x4ULL,
0x408ULL, 0x0ULL, 0x0ULL, 0x0ULL, 0x800040000100000ULL, 0x402000000000000ULL, 0x20ULL, 0x20000000000000ULL,
0x100008000ULL, 0x400000080ULL, 0x800000100800000ULL, 0x4000000002000000ULL, 0x2000000000100ULL, 0x40000200000ULL, 0x140000000080ULL, 0x4000000000000000ULL,
0x200000000000ULL, 0x0ULL, 0x4008000000ULL, 0x810000000ULL, 0x800400000020002ULL, 0x0ULL, 0x40ULL, 0x8900ULL,
0x20000ULL, 0x200000000ULL, 0x800000200080028ULL, 0x0ULL, 0x20040000000000ULL, 0x20000000000000ULL, 0x8000000ULL, 0x10000010000000ULL,
0x0ULL, 0x100040ULL, 0x2002800000400000ULL, 0x400000000000000ULL, 0x10000000000000ULL, 0x0ULL, 0x4000008ULL, 0x0ULL,
0x40001000ULL, 0x8000000000040000ULL, 0x0ULL, 0x1000000000800ULL, 0x0ULL, 0x400000ULL, 0x20000ULL, 0x0ULL,
0x0ULL, 0x400000080000ULL, 0x40400080000000ULL, 0x4000000000000000ULL, 0x80000000ULL, 0x0ULL, 0x40800000004000ULL, 0x10000000ULL,
0x0ULL, 0x0ULL, 0x0ULL, 0x20000000100000ULL, 0x0ULL, 0x80000000000000ULL, 0x200000000ULL, 0x0ULL,
0x0ULL, 0x8000ULL, 0x200000000000ULL, 0x1000000003000000ULL, 0x30000000ULL, 0x1ULL, 0x0ULL, 0x10000ULL,
0x0ULL, 0x20000000040ULL, 0x1000004000001ULL, 0x800000000ULL, 0x800001000000000ULL, 0x80000000ULL, 0x0ULL, 0x20200ULL,
0x8000000000000ULL, 0x4ULL, 0x0ULL, 0x800800000008000ULL, 0x1000000000ULL, 0x1000000ULL, 0x40000000000008ULL, 0x0ULL,
0x4000ULL, 0x1000008000ULL, 0x400000000ULL, 0x20000000ULL, 0x1090000000ULL, 0x4000000000ULL, 0x400000000ULL, 0x40000001ULL,
0x0ULL, 0x0ULL, 0x2000000000004ULL, 0x4010000000020ULL, 0x2020000000ULL, 0x0ULL, 0x8001000000000ULL, 0x2000000ULL,
0x200000ULL, 0x20ULL, 0x0ULL, 0x80000000000ULL, 0x0ULL, 0x0ULL, 0x1000000000000ULL, 0x0ULL,
0x800000000ULL, 0x4200000000ULL, 0x220ULL, 0x1000000000ULL, 0x80000000ULL, 0x1000ULL, 0x20000000800000ULL, 0x0ULL,
0x1000000ULL, 0x20008ULL, 0x10000000ULL, 0x0ULL, 0x3000000000ULL, 0x4000004000ULL, 0x102000000000000ULL, 0x10000200000ULL,
0x80002000000ULL, 0x0ULL, 0x10000100000a000ULL, 0x1000000000ULL, 0x0ULL, 0x0ULL, 0x400000ULL, 0x8000000ULL,
0x0ULL, 0x0ULL, 0x2ULL, 0x4000000000000ULL, 0x100000000000000ULL, 0x800000000000ULL, 0x40000000000ULL, 0x0ULL,
0x400400000081ULL, 0x2000408000000ULL, 0x80200000000100ULL, 0x40000000ULL, 0x8000000000000ULL, 0x0ULL, 0x0ULL, 0x8100000000000800ULL,
0x2000000ULL, 0x20000000080ULL, 0x0ULL, 0x0ULL, 0x180010000000ULL, 0x10001000000000ULL, 0x4ULL, 0x2080000000ULL,
0x20000000ULL, 0x0ULL, 0x1000000ULL, 0x0ULL, 0x0ULL, 0x80000000000ULL, 0x100000000008ULL, 0x0ULL,
0x0ULL, 0x81000500ULL, 0x100000000ULL, 0x40000001000000ULL, 0x10000000ULL, 0x8ULL, 0x400010001000000ULL, 0x400000000ULL,
};

/*
void GenTable()
{
	const int words = (64000 + 63)/64;
	u64 table[words] = {0};

	for (int ii = 0; ii < EXCEPTION_COUNT; ++ii)
	{
		u32 n = EXCEPTIONS[ii].n;

		table[n >> 6] |= (u64)1 << (n & 63);
	}

	cout << "static const u64 EXCEPT_TABLE[" << words << "] = {";
	for (int ii = 0; ii < words; ++ii)
	{
		if ((ii & 7) == 0) cout << endl;
		cout << hex << "0x" << setfill('0') << table[ii] << dec << "ULL, ";
	}
	cout << endl << "};" << endl;
}
*/

/*
	ChooseMatrix

		This function determines the check matrix to use based on the
	given message bytes and bytes per block.
*/

Result Codec::ChooseMatrix(int message_bytes, int block_bytes)
{
	CAT_IF_DUMP(cout << endl << "---- ChooseMatrix ----" << endl << endl;)

	//GenTable();

	// Validate input
	if CAT_UNLIKELY(message_bytes < 1 || block_bytes < 1)
		return R_BAD_INPUT;

	// Calculate message block count
	_block_bytes = block_bytes;
	_block_count = (message_bytes + _block_bytes - 1) / _block_bytes;
	_block_next_prime = NextPrime16(_block_count);

	// Validate block count
	if CAT_UNLIKELY(_block_count < CAT_WIREHAIR_MIN_N)
		return R_TOO_SMALL;
	if CAT_UNLIKELY(_block_count > CAT_WIREHAIR_MAX_N)
		return R_TOO_LARGE;

	CAT_IF_DUMP(cout << "Total message = " << message_bytes << " bytes.  Block bytes = " << _block_bytes << endl;)
	CAT_IF_DUMP(cout << "Block count = " << _block_count << " +Prime=" << _block_next_prime << endl;)

	/*
		Calculate the number of dense rows

			The choice of dense row count affects the seeding, so
		this had to be determined before picking seeds.  Similarly,
		the peeling matrix parameters greatly affect the best dense
		row count.

			I simulated without heavy rows added for selected
		values of N and determined the best dense row count to
		use to achieve 30% invertibility.

			Then I plugged the data into MATLAB's cftool to do curve
		fitting.  For different ranges of N there are different
		dominant effects and that is reflected below:
	*/

	// If N is small,
	u16 dense_count;
	if (_block_count < 256)
	{
		// Calculate dense count from math expression
		if (_block_count == 2) dense_count = 2;
		else if (_block_count == 3) dense_count = 6;
		else dense_count = 10 + SquareRoot16(_block_count) / 2 + (_block_count / 50);
	}
	else if (_block_count <= 4096) // Medium N:
	{
		// Square root-dominant region
		dense_count = 11 + SquareRoot16(_block_count) + (u16)(_block_count / 300);
	}
	else if (_block_count <= 32768)
	{
		// Linear-dominant region
		dense_count = 22 + (_block_count / 100);
	}
	else if (_block_count <= 44000)
	{
		// Linear-dominant region
		dense_count = 26 + (_block_count / 114);
	}
	else if (_block_count <= 52500)
	{
		// Linear-dominant region
		dense_count = 74 + (_block_count / 128);
	}
	else
	{
		// Avalanche-dominant region
		dense_count = 880 - (_block_count / 128);
	}

	// Round up to the next D s.t. D Mod 4 = 2 (see above)
	switch (dense_count & 3)
	{
	case 0: dense_count += 2; break;
	case 1: dense_count += 1; break;
	case 2: break;
	case 3: dense_count += 3; break;
	}

	if (dense_count < 14)
	{
		switch (dense_count)
		{
		case 2: _d_seed = 0; break; // Seed doesn't matter
		case 6: _d_seed = 67; break;
		default: return R_BAD_DENSE_SEED;
		}
	}
	else
	{
		if (dense_count > 486)
			return R_BAD_DENSE_SEED;

		// Lookup dense seed given D
		_d_seed = DENSE_SEEDS[(dense_count - 14) / 4];
	}

	_dense_count = dense_count;

	/*
		Now we need to select a seed that will cause
		the check matrix to be full rank.

		Because a randomly selected matrix is full rank
		most of the time, we just need to handle the
		exceptional cases, which leads to a small table.

		For small N, use a lookup table for ideal seed
		since tuning is more important for these cases.
	*/

	// If N is small,
	if (_block_count <= SMALL_SEED_MAX)
	{
		// Lookup seeds from table
		_p_seed = SMALL_PEEL_SEEDS[_block_count];
	}
	else
	{
		// If default seed doesn't work,
		if (EXCEPT_SEEDS[_block_count >> 6] & ((u64)1 << (_block_count & 63)))
		{
			switch (_block_count)
			{
			case 1473:
			case 2653:
			case 6421:
			case 6835:
			case 8091:
			case 9462:
			case 9712:
			case 14853:
			case 19428:
			case 22635:
			case 24015:
			case 31954:
			case 32914:
			case 36023:
			case 39741:
			case 41612:
			case 41797:
			case 51792:
			case 54223:
			case 54403:
			case 54959:
				// Use tertiary backup seed
				_p_seed = 1;
				break;
			default:
				// Use secondary backup seed
				_p_seed = 0;
			}
		}
		else
		{
			// Use default seed
			_p_seed = _block_count;
		}
	}

	CAT_IF_DUMP(cout << "Peel seed = " << _p_seed << "  Dense seed = " << _d_seed << endl;)

	_mix_count = _dense_count + CAT_HEAVY_ROWS;
	_mix_next_prime = NextPrime16(_mix_count);

	CAT_IF_DUMP(cout << "Mix count = " << _mix_count << " +Prime=" << _mix_next_prime << endl;)

	// Initialize lists
	_peel_head_rows = LIST_TERM;
	_peel_tail_rows = 0;
	_defer_head_rows = LIST_TERM;

	return R_WIN;
}

/*
	SolveMatrix

		This function attempts to solve the matrix, given that the
	matrix is currently square and may be solvable at this point
	without any additional blocks.

		It performs the peeling, compression, and GE steps:

	(1) Peeling:
		GreedyPeeling()

	(2) Compression:

		Allocate GE and Compression matrix now that the size is known:

			AllocateMatrix()

		Produce the Compression matrix:

			SetDeferredColumns()
			SetMixingColumnsForDeferredRows()
			PeelDiagonal()

		Produce the GE matrix:

			CopyDeferredRows()
			MultiplyDenseRows()
			AddInvertibleGF2Matrix()

	(3) Gaussian Elimination
		Triangle()
*/

Result Codec::SolveMatrix()
{
	// (1) Peeling

	GreedyPeeling();

	CAT_IF_DUMP( PrintPeeled(); )
	CAT_IF_DUMP( PrintDeferredRows(); )
	CAT_IF_DUMP( PrintDeferredColumns(); )

	// (2) Compression

	if (!AllocateMatrix())
		return R_OUT_OF_MEMORY;

	SetDeferredColumns();
	SetMixingColumnsForDeferredRows();
	PeelDiagonal();
	CopyDeferredRows();
	MultiplyDenseRows();
	SetHeavyRows();

	// Add invertible matrix to mathematically tie dense rows to dense mixing columns
	if (!AddInvertibleGF2Matrix(_ge_matrix, _defer_count, _ge_pitch, _dense_count))
		return R_TOO_SMALL;

#if defined(CAT_DUMP_CODEC_DEBUG) || defined(CAT_DUMP_GE_MATRIX)
	cout << "After Compress:" << endl;
	PrintGEMatrix();
#endif
	CAT_IF_DUMP( PrintCompressMatrix(); )

	// (3) Gaussian Elimination

	SetupTriangle();
	if (!Triangle())
	{
		CAT_IF_DUMP( cout << "After Triangle FAILED:" << endl; )
		CAT_IF_DUMP( PrintGEMatrix(); )
		CAT_IF_DUMP( PrintExtraMatrix(); )

		return R_MORE_BLOCKS;
	}

#if defined(CAT_DUMP_CODEC_DEBUG) || defined(CAT_DUMP_GE_MATRIX)
	cout << "After Triangle:" << endl;
	PrintGEMatrix();
	PrintExtraMatrix();
#endif

	return R_WIN;
}

/*
	GenerateRecoveryBlocks

		This function generates the recovery blocks after the
	Triangle() function succeeds in solving the matrix.

		It performs the final Substitution step:

	(4) Substitution:

		Solves across GE matrix rows:

			InitializeColumnValues()
			MultiplyDenseValues()
			AddSubdiagonalValues()

		Solves up GE matrix columns:

			BackSubstituteAboveDiagonal()

		Solves remaining columns:

			Substitute()
*/

void Codec::GenerateRecoveryBlocks()
{
	// (4) Substitution

	InitializeColumnValues();
	MultiplyDenseValues();
	AddSubdiagonalValues();
	BackSubstituteAboveDiagonal();
	Substitute();
}

/*
	ResumeSolveMatrix

		This function resumes solving the matrix if the initial SolveMatrix()
	failed at finding a pivot.  With heavy rows added this is pretty rare, so
	in theory this function could be written inefficiently and it would not
	affect the average run time much.  But that just wouldn't be right.

		Extra rows have been added to both the GE matrix and the heavy matrix.
	Since the heavy matrix does not cover all the columns of the GE matrix,
	the new rows are staged in the GE matrix and then copied into the heavy
	matrix after they get into range of the heavy columns.
*/

Result Codec::ResumeSolveMatrix(u32 id, const void *block)
{
	CAT_IF_DUMP(cout << endl << "---- ResumeSolveMatrix ----" << endl << endl;)

	if (!block) return R_BAD_INPUT;

	// If there is no room for it,
	u16 row_i, ge_row_i, new_pivot_i;
	if (_row_count >= _block_count + _extra_count)
	{
		const u16 first_heavy_row = _defer_count + _dense_count;

		new_pivot_i = 0;

		// For each pivot in the list,
		for (u16 pivot_i = _next_pivot; pivot_i < _pivot_count; ++pivot_i)
		{
			// If unused row is extra,
			u16 ge_row_i = _pivots[pivot_i];
			if (ge_row_i >= first_heavy_row && ge_row_i < (first_heavy_row + _extra_count))
			{
				// Re-use it
				new_pivot_i = pivot_i;
				break;
			}
		}

		// If nothing was found, return error
		if (!new_pivot_i) return R_NEED_MORE_EXTRA;

		// Look up row indices
		ge_row_i = _pivots[new_pivot_i];
		row_i = _ge_row_map[ge_row_i];
	}
	else
	{
		// Add extra rows to the end of the pivot list
		new_pivot_i = _pivot_count++;
		row_i = _row_count++;
		ge_row_i = _defer_count + _dense_count + row_i - _block_count;
		_ge_row_map[ge_row_i] = row_i;
		_pivots[new_pivot_i] = ge_row_i;

		/*
			Before the extra rows are converted to heavy, the new rows
			are added to the end of the pivot list.  And after the extra
			rows are converted to heavy rows, new rows that come in are
			also heavy and should also be at the end of the pivot list.

			So, this doesn't need to change based on what stage of the
			GE solver is running through at this point.
		*/
	}

	CAT_IF_DUMP(cout << "Resuming using row slot " << row_i << " and GE row " << ge_row_i << endl;)

	// Update row data needed at this point
	PeelRow * CAT_RESTRICT row = &_peel_rows[row_i];
	row->id = id;

	// Copy new block to input blocks
	u8 * CAT_RESTRICT block_store_dest = _input_blocks + _block_bytes * row_i;
	if (id != _block_count - 1)
		memcpy(block_store_dest, block, _block_bytes);
	else
	{
		memcpy(block_store_dest, block, _output_final_bytes);
		memset(block_store_dest + _output_final_bytes, 0, _block_bytes - _output_final_bytes);
	}

	// Generate new GE row
	u64 * CAT_RESTRICT ge_new_row = _ge_matrix + _ge_pitch * ge_row_i;
	memset(ge_new_row, 0, _ge_pitch * sizeof(u64));

	u16 peel_weight, peel_a, peel_x, mix_a, mix_x;
	GeneratePeelRow(id, _p_seed, _block_count, _mix_count,
		peel_weight, peel_a, peel_x, mix_a, mix_x);

	// Store row parameters
	row->peel_weight = peel_weight;
	row->peel_a = peel_a;
	row->peel_x0 = peel_x;
	row->mix_a = mix_a;
	row->mix_x0 = mix_x;

	// Generate mixing bits in GE row
	u16 ge_column_i = mix_x + _defer_count;
	ge_new_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
	IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
	ge_column_i = mix_x + _defer_count;
	ge_new_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
	IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
	ge_column_i = mix_x + _defer_count;
	ge_new_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);

	// Generate peeled bits in GE row
	for (;;)
	{
		// If column is peeled,
		PeelColumn * CAT_RESTRICT ref_col = &_peel_cols[peel_x];
		if (ref_col->mark == MARK_PEEL)
		{
			// Add compress row to the new GE row
			u16 row_i = ref_col->peel_row;
			const u64 * CAT_RESTRICT ge_src_row = _compress_matrix + _ge_pitch * row_i;
			for (int ii = 0; ii < _ge_pitch; ++ii) ge_new_row[ii] ^= ge_src_row[ii];
		}
		else
		{
			// Set bit for this deferred column
			u16 ge_column_i = ref_col->ge_column;
			ge_new_row[ge_column_i >> 6] ^= (u64)1 << (ge_column_i & 63);
		}

		if (--peel_weight <= 0) break;

		IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);
	}

	// For each pivot-found column up to the start of the heavy columns,
	u64 ge_mask = 1;
	for (u16 pivot_j = 0; pivot_j < _next_pivot && pivot_j < _first_heavy_column;
		++pivot_j, ge_mask = CAT_ROL64(ge_mask, 1))
	{
		// If bit is set,
		int word_offset = pivot_j >> 6;
		u64 * CAT_RESTRICT rem_row = &ge_new_row[word_offset];
		if (*rem_row & ge_mask)
		{
			u16 ge_row_j = _pivots[pivot_j];
			u64 * CAT_RESTRICT ge_pivot_row = _ge_matrix + word_offset + _ge_pitch * ge_row_j;
			u64 row0 = (*ge_pivot_row & ~(ge_mask - 1)) ^ ge_mask;

			// Add previous pivot row to new row
			*rem_row ^= row0;
			for (int ii = 1; ii < _ge_pitch - word_offset; ++ii)
				rem_row[ii] ^= ge_pivot_row[ii];
		}
	}

	// If next pivot is not heavy,
	if (_next_pivot < _first_heavy_column)
	{
		// If the next pivot was not found on this row,
		if (!(ge_new_row[_next_pivot >> 6] & ((u64)1 << (_next_pivot & 63))))
			return R_MORE_BLOCKS; // Maybe next time...

		// Swap out the pivot index for this one
		_pivots[new_pivot_i] = _pivots[_next_pivot];
		_pivots[_next_pivot] = ge_row_i;
	}
	else
	{
		// For each heavy column,
		const u16 column_count = _defer_count + _mix_count;
		const u16 first_heavy_row = _dense_count + _defer_count;
		u16 heavy_row_i = ge_row_i - first_heavy_row;
		u8 * CAT_RESTRICT heavy_row = _heavy_matrix + _heavy_pitch * heavy_row_i;
		for (u16 ge_column_j = _first_heavy_column; ge_column_j < column_count; ++ge_column_j)
		{
			u16 heavy_col_j = ge_column_j - _first_heavy_column;
			u8 bit_j = (u8)(ge_new_row[ge_column_j >> 6] >> (ge_column_j & 63)) & 1;

			// Copy bit into column byte
			heavy_row[heavy_col_j] = bit_j;
		}

		// For each pivot-found column in the heavy columns,
		for (u16 pivot_j = _first_heavy_column; pivot_j < _next_pivot; ++pivot_j)
		{
			// If column is zero,
			u16 heavy_col_j = pivot_j - _first_heavy_column;
			u8 code_value = heavy_row[heavy_col_j];
			if (!code_value) continue; // Skip it

			// If previous row is heavy,
			u16 ge_row_j = _pivots[pivot_j];
			if (ge_row_j >= first_heavy_row)
			{
				// Calculate coefficient of elimination
				u16 heavy_row_j = ge_row_j - first_heavy_row;
				u8 * CAT_RESTRICT pivot_row = _heavy_matrix + _heavy_pitch * heavy_row_j;
				u8 pivot_code = pivot_row[heavy_col_j];
				const u16 start_column = heavy_col_j + 1;
				if (pivot_code == 1)
				{
					if (code_value != 1)
					{
						// heavy[m+] += exist[m+] * code_value
						GF256MemMulAdd(heavy_row + start_column, code_value, pivot_row + start_column, _heavy_columns - start_column);
					}
					else
					{
						// heavy[m+] += exist[m+]
						memxor(heavy_row + start_column, pivot_row + start_column, _heavy_columns - start_column);
					}
				}
				else
				{
					u8 eliminator = (code_value == 1) ? INV_TABLE[pivot_code] : GF256Divide(code_value, pivot_code);

					// Store eliminator for later
					heavy_row[heavy_col_j] = eliminator;

					// heavy[m+] += exist[m+] * eliminator
					GF256MemMulAdd(heavy_row + start_column, eliminator, pivot_row + start_column, _heavy_columns - start_column);
				}
			}
			else
			{
				// For each remaining column,
				u64 * CAT_RESTRICT other_row = _ge_matrix + _ge_pitch * ge_row_j;
				u16 ge_column_k = pivot_j + 1;
				u64 ge_mask = (u64)1 << (ge_column_k & 63);
				for (; ge_column_k < column_count; ++ge_column_k, ge_mask = CAT_ROL64(ge_mask, 1))
				{
					// If bit is set,
					if (other_row[ge_column_k >> 6] & ge_mask)
					{
						// Add in the code value for this column
						heavy_row[ge_column_k - _first_heavy_column] ^= code_value;
					}
				}
			} // end if row is heavy
		} // next column

		// If the next pivot was not found on this heavy row,
		u16 next_heavy_col = _next_pivot - _first_heavy_column;
		if (!heavy_row[next_heavy_col])
			return R_MORE_BLOCKS; // Maybe next time...

		// If a non-heavy pivot just got moved into heavy pivot list,
		if (_next_pivot < _first_heavy_pivot)
		{
			// Swap out the pivot index for this one
			_pivots[new_pivot_i] = _pivots[_first_heavy_pivot];
			_pivots[_first_heavy_pivot] = _pivots[_next_pivot];

			// And move the first heavy pivot up one to cover the hole
			++_first_heavy_pivot;
		}
		else
		{
			// Swap out the pivot index for this one
			_pivots[new_pivot_i] = _pivots[_next_pivot];
		}

		_pivots[_next_pivot] = ge_row_i;
	}

	// NOTE: Pivot was found and is definitely not set anywhere else
	// so it doesn't need to be cleared from any other GE rows.

	// If just starting heavy columns,
	if (++_next_pivot == _first_heavy_column)
		InsertHeavyRows();

	// Resume Triangle() at next pivot to determine
	return Triangle() ? R_WIN : R_MORE_BLOCKS;
}

#if defined(CAT_ALL_ORIGINAL)

/*
	IsAllOriginalData

		This function verifies that all of the original N data blocks
	have been received.
*/

bool Codec::IsAllOriginalData()
{
	// Re-purpose and initialize an array to store whether or not each row id needs to be regenerated
	u8 * CAT_RESTRICT copied_rows = reinterpret_cast<u8*>( _recovery_blocks );
	memset(copied_rows, 0, _block_count);

	// Copy any original message rows that were received:
	// For each row,
	PeelRow * CAT_RESTRICT row = _peel_rows;
	u32 seen_rows = 0;
	for (u16 row_i = 0; row_i < _row_count; ++row_i, ++row)
	{
		u32 id = row->id;

		// If the row identifier indicates it is part of the original message data,
		if (id < _block_count)
		{
			// If not already marked copied,
			if (!copied_rows[id])
			{
				copied_rows[id] = 1;
				++seen_rows;
			}
		}
	}

	// If all rows were seen,
	return seen_rows >= _block_count;
}

#endif // CAT_ALL_ORIGINAL

/*
	ReconstructBlock

		This function reconstructs an original block from the recovery
	blocks, which is much slower than copying from the input data, so
	should be done selectively.  This is only done during decoding.

	Precondition: DecodeFeed() has returned success
*/

Result Codec::ReconstructBlock(u16 row_i, void * CAT_RESTRICT dest) {
	CAT_IF_DUMP(cout << endl << "---- ReconstructBlock ----" << endl << endl;)

	// Validate input
	if CAT_UNLIKELY(!dest) return R_BAD_INPUT;

	// Regenerate any single row that got lost:

	u32 block_bytes = _block_bytes;

	// For last row, use final byte count
	if (row_i == _block_count - 1)
		block_bytes = _output_final_bytes;

	CAT_IF_DUMP(cout << "Regenerating row " << row_i << ":";)

	u16 peel_weight, peel_a, peel_x, mix_a, mix_x;
	GeneratePeelRow(row_i, _p_seed, _block_count, _mix_count,
		peel_weight, peel_a, peel_x, mix_a, mix_x);

	// Remember first column (there is always at least one)
	u8 * CAT_RESTRICT first = _recovery_blocks + _block_bytes * peel_x;

	CAT_IF_DUMP(cout << " " << peel_x;)

	// If peeler has multiple columns,
	if (peel_weight > 1)
	{
		--peel_weight;

		IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

		CAT_IF_DUMP(cout << " " << peel_x;)

		// Combine first two columns into output buffer (faster than memcpy + memxor)
		memxor_set(dest, first, _recovery_blocks + _block_bytes * peel_x, block_bytes);

		// For each remaining peeler column,
		while (--peel_weight > 0)
		{
			IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

			CAT_IF_DUMP(cout << " " << peel_x;)

			// Mix in each column
			memxor(dest, _recovery_blocks + _block_bytes * peel_x, block_bytes);
		}

		// Mix first mixer block in directly
		memxor(dest, _recovery_blocks + _block_bytes * (_block_count + mix_x), block_bytes);
	}
	else
	{
		// Mix first with first mixer block (faster than memcpy + memxor)
		memxor_set(dest, first, _recovery_blocks + _block_bytes * (_block_count + mix_x), block_bytes);
	}

	CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

	// Combine remaining two mixer columns together:
	IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
	const u8 *mix0_src = _recovery_blocks + _block_bytes * (_block_count + mix_x);
	CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

	IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
	const u8 *mix1_src = _recovery_blocks + _block_bytes * (_block_count + mix_x);
	CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

	memxor_add(dest, mix0_src, mix1_src, block_bytes);

	CAT_IF_DUMP(cout << endl;)

	return R_WIN;
}


/*
	ReconstructOutput

		This function reconstructs the output by copying inputs
	that were from the first N blocks, and regenerating the rest.
	This is only done during decoding.

	Precondition: DecodeFeed() has returned success
*/

Result Codec::ReconstructOutput(void * CAT_RESTRICT message_out)
{
	CAT_IF_DUMP(cout << endl << "---- ReconstructOutput ----" << endl << endl;)

	// Validate input
	if CAT_UNLIKELY(!message_out) return R_BAD_INPUT;
	u8 * CAT_RESTRICT output_blocks = reinterpret_cast<u8 *>( message_out );

#if defined(CAT_COPY_FIRST_N)
	// Re-purpose and initialize an array to store whether or not each row id needs to be regenerated
	u8 * CAT_RESTRICT copied_rows = reinterpret_cast<u8*>( _peel_cols );
	memset(copied_rows, 0, _block_count);

	// Copy any original message rows that were received:
	// For each row,
	PeelRow * CAT_RESTRICT row = _peel_rows;
	const u8 * CAT_RESTRICT src = _input_blocks;
	for (u16 row_i = 0; row_i < _row_count; ++row_i, ++row, src += _block_bytes)
	{
		u32 id = row->id;

		// If the row identifier indicates it is part of the original message data,
		if (id < _block_count)
		{
			CAT_IF_DUMP(cout << "Copying received row " << id << endl;)

			u8 * CAT_RESTRICT dest = output_blocks + _block_bytes * id;
			int bytes = (id != _block_count - 1) ? _block_bytes : _output_final_bytes;
			memcpy(dest, src, bytes);

			copied_rows[id] = 1;
		}
	}
#endif // CAT_COPY_FIRST_N

	// Regenerate any rows that got lost:

	// For each row,
	u8 * CAT_RESTRICT dest = output_blocks;
	u32 block_bytes = _block_bytes;
#if defined(CAT_COPY_FIRST_N)
	u8 * CAT_RESTRICT copied_row = copied_rows;
	for (u16 row_i = 0; row_i < _block_count; ++row_i, dest += _block_bytes, ++copied_row)
	{
		// If already copied, skip it
		if (*copied_row)
			continue;
#else
	for (u16 row_i = 0; row_i < _block_count; ++row_i, dest += _block_bytes)
	{
#endif
		// For last row, use final byte count
		if (row_i == _block_count - 1)
			block_bytes = _output_final_bytes;

		CAT_IF_DUMP(cout << "Regenerating row " << row_i << ":";)

		u16 peel_weight, peel_a, peel_x, mix_a, mix_x;
		GeneratePeelRow(row_i, _p_seed, _block_count, _mix_count,
			peel_weight, peel_a, peel_x, mix_a, mix_x);

		// Remember first column (there is always at least one)
		u8 * CAT_RESTRICT first = _recovery_blocks + _block_bytes * peel_x;

		CAT_IF_DUMP(cout << " " << peel_x;)

		// If peeler has multiple columns,
		if (peel_weight > 1)
		{
			--peel_weight;

			IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

			CAT_IF_DUMP(cout << " " << peel_x;)

			// Combine first two columns into output buffer (faster than memcpy + memxor)
			memxor_set(dest, first, _recovery_blocks + _block_bytes * peel_x, block_bytes);

			// For each remaining peeler column,
			while (--peel_weight > 0)
			{
				IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

				CAT_IF_DUMP(cout << " " << peel_x;)

				// Mix in each column
				memxor(dest, _recovery_blocks + _block_bytes * peel_x, block_bytes);
			}

			// Mix first mixer block in directly
			memxor(dest, _recovery_blocks + _block_bytes * (_block_count + mix_x), block_bytes);
		}
		else
		{
			// Mix first with first mixer block (faster than memcpy + memxor)
			memxor_set(dest, first, _recovery_blocks + _block_bytes * (_block_count + mix_x), block_bytes);
		}

		CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

		// Combine remaining two mixer columns together:
		IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
		const u8 *mix0_src = _recovery_blocks + _block_bytes * (_block_count + mix_x);
		CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

		IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
		const u8 *mix1_src = _recovery_blocks + _block_bytes * (_block_count + mix_x);
		CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

		memxor_add(dest, mix0_src, mix1_src, block_bytes);

		CAT_IF_DUMP(cout << endl;)
	} // next row

	return R_WIN;
}


//// Memory Management

Codec::Codec()
{
	// Workspace
	_recovery_blocks = 0;
	_workspace_allocated = 0;

	// Matrix
	_compress_matrix = 0;
	_ge_allocated = 0;

	// Input
	_input_blocks = 0;
	_input_allocated = 0;
}

Codec::~Codec()
{
	FreeWorkspace();
	FreeMatrix();
	FreeInput();
}

void Codec::SetInput(const void * CAT_RESTRICT message_in)
{
	FreeInput();

	// Set input blocks to the input message
	_input_blocks = (u8*)message_in;
	_input_allocated = 0;
}

bool Codec::AllocateInput()
{
	CAT_IF_DUMP(cout << endl << "---- AllocateInput ----" << endl << endl;)

	// If need to allocate more,
	u32 size = (_block_count + _extra_count) * _block_bytes;
	if (_input_allocated < size)
	{
		FreeInput();

		// Allocate input blocks
		_input_blocks = new u8[size];
		if (!_input_blocks) return false;
		_input_allocated = size;
	}

	return true;
}

void Codec::FreeInput()
{
	if (_input_allocated > 0 && _input_blocks)
	{
		delete []_input_blocks;
		_input_blocks = 0;
	}

	_input_allocated = 0;
}

bool Codec::AllocateMatrix()
{
	CAT_IF_DUMP(cout << endl << "---- AllocateMatrix ----" << endl << endl;)

	// GE matrix
	const int ge_cols = _defer_count + _mix_count;
	const int ge_rows = _defer_count + _dense_count + _extra_count + 1; // One extra for workspace
	const int ge_pitch = (ge_cols + 63) / 64;
	const u32 ge_matrix_words = ge_rows * ge_pitch;

	// Compression matrix
	const int compress_rows = _block_count;
	const u32 compress_matrix_words = compress_rows * ge_pitch;

	// Pivots
	const int pivot_count = ge_cols + _extra_count;
	const int pivot_words = pivot_count * 2 + ge_cols;

	// Heavy
	const int heavy_rows = CAT_HEAVY_ROWS + _extra_count;
	const int heavy_cols = _mix_count < CAT_HEAVY_MAX_COLS ? _mix_count : CAT_HEAVY_MAX_COLS;
	const int heavy_pitch = (heavy_cols + 3 + 3) & ~3; // Round up columns+3 to next multiple of 4
	const int heavy_bytes = heavy_pitch * heavy_rows;

	// Calculate buffer size
	u32 size = ge_matrix_words * sizeof(u64) + compress_matrix_words * sizeof(u64) + pivot_words * sizeof(u16) + heavy_bytes;

	// If need to allocate more,
	if (_ge_allocated < size)
	{
		FreeMatrix();

		u8 * CAT_RESTRICT matrix = new u8[size];
		if (!matrix) return false;
		_ge_allocated = size;
		_compress_matrix = reinterpret_cast<u64 *>( matrix );
	}

	// Store pointers
	_ge_pitch = ge_pitch;
	_ge_matrix = _compress_matrix + compress_matrix_words;
	_heavy_pitch = heavy_pitch;
	_heavy_columns = heavy_cols;
	_first_heavy_column = _defer_count + _mix_count - heavy_cols;
	_heavy_matrix = reinterpret_cast<u8 *>( _ge_matrix + ge_matrix_words );
	_pivots = reinterpret_cast<u16 *>( _heavy_matrix + heavy_bytes );
	_ge_row_map = _pivots + pivot_count;
	_ge_col_map = _ge_row_map + pivot_count;

	CAT_IF_DUMP(cout << "GE matrix is " << ge_rows << " x " << ge_cols << " with pitch " << ge_pitch << " consuming " << ge_matrix_words * sizeof(u64) << " bytes" << endl;)
	CAT_IF_DUMP(cout << "Compress matrix is " << compress_rows << " x " << ge_cols << " with pitch " << ge_pitch << " consuming " << compress_matrix_words * sizeof(u64) << " bytes" << endl;)
	CAT_IF_DUMP(cout << "Allocated " << pivot_count << " pivots, consuming " << pivot_words*2 << " bytes" << endl;)
	CAT_IF_DUMP(cout << "Allocated " << CAT_HEAVY_ROWS << " heavy rows, consuming " << heavy_bytes << " bytes" << endl;)

	// Clear entire Compression matrix
	memset(_compress_matrix, 0, compress_matrix_words * sizeof(u64));

	// Clear entire GE matrix
	memset(_ge_matrix, 0, ge_cols * ge_pitch * sizeof(u64));

	return true;
}

void Codec::FreeMatrix()
{
	if (_compress_matrix)
	{
		u8 * CAT_RESTRICT matrix = reinterpret_cast<u8 *>( _compress_matrix );
		delete []matrix;
		_compress_matrix = 0;
	}

	_ge_allocated = 0;
}

bool Codec::AllocateWorkspace()
{
	CAT_IF_DUMP(cout << endl << "---- AllocateWorkspace ----" << endl << endl;)

	// Count needed rows and columns
	const u32 recovery_size = (_block_count + _mix_count + 1) * _block_bytes; // +1 for temporary space
	const u32 row_count = _block_count + _extra_count;
	const u32 column_count = _block_count;

	// Calculate size
	u32 size = recovery_size + sizeof(PeelRow) * row_count
		+ sizeof(PeelColumn) * column_count + sizeof(PeelRefs) * column_count;
	if (_workspace_allocated < size)
	{
		FreeWorkspace();

		// Allocate workspace
		_recovery_blocks = new u8[size];
		if (!_recovery_blocks) return false;
		_workspace_allocated = size;
	}

	// Set pointers
	_peel_rows = reinterpret_cast<PeelRow *>( _recovery_blocks + recovery_size );
	_peel_cols = reinterpret_cast<PeelColumn *>( _peel_rows + row_count );
	_peel_col_refs = reinterpret_cast<PeelRefs *>( _peel_cols + column_count );

	CAT_IF_DUMP(cout << "Memory overhead for workspace = " << size << " bytes" << endl;)

	// Initialize columns
	for (int ii = 0; ii < _block_count; ++ii)
	{
		_peel_col_refs[ii].row_count = 0;
		_peel_cols[ii].w2_refs = 0;
		_peel_cols[ii].mark = MARK_TODO;
	}

	return true;
}

void Codec::FreeWorkspace()
{
	if (_recovery_blocks)
	{
		delete []_recovery_blocks;
		_recovery_blocks = 0;
	}

	_workspace_allocated = 0;
}


//// Diagnostic

#if defined(CAT_DUMP_CODEC_DEBUG) || defined(CAT_DUMP_GE_MATRIX)

void Codec::PrintGEMatrix()
{
	const int rows = _dense_count + _defer_count;
	const int cols = _defer_count + _mix_count;

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

	const int heavy_rows = CAT_HEAVY_ROWS;
	const int heavy_cols = _heavy_columns;
	cout << "Heavy submatrix is " << heavy_rows << " x " << heavy_cols << ":" << endl;

	for (int ii = 0; ii < heavy_rows; ++ii)
	{
		for (int jj = 0; jj < heavy_cols; ++jj)
		{
			cout << hex << setw(2) << setfill('0') << (int)_heavy_matrix[_heavy_pitch * (ii + _extra_count) + jj] << dec << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void Codec::PrintExtraMatrix()
{
	cout << endl << "Extra rows: " << endl;

	// For each pivot,
	int extra_count = 0;
	const u16 column_count = _defer_count + _mix_count;
	const u16 first_heavy_row = _defer_count + _dense_count;
	for (u16 pivot_i = 0; pivot_i < _pivot_count; ++pivot_i)
	{
		// If row is extra,
		u16 ge_row_i = _pivots[pivot_i];
		if (ge_row_i >= first_heavy_row && ge_row_i < first_heavy_row + _extra_count)
		{
			u64 * CAT_RESTRICT ge_row = _ge_matrix + _ge_pitch * ge_row_i;
			u16 heavy_row_i = ge_row_i - first_heavy_row;
			u8 * CAT_RESTRICT heavy_row = _heavy_matrix + _heavy_pitch * heavy_row_i;

			cout << "row=" << ge_row_i << " : light={ ";

			// For each non-heavy column,
			for (u16 ge_column_i = 0; ge_column_i < _first_heavy_column; ++ge_column_i)
			{
				// If column is non-zero,
				u64 ge_mask = (u64)1 << (ge_column_i & 63);
				if (ge_row[ge_column_i >> 6] & ge_mask)
				{
					cout << '1';
				}
				else
				{
					cout << '0';
				}
			}

			cout << " } heavy=(";

			// For each heavy column,
			for (u16 ge_column_i = _first_heavy_column; ge_column_i < column_count; ++ge_column_i)
			{
				u16 heavy_col_i = ge_column_i - _first_heavy_column;
				u8 code_value = heavy_row[heavy_col_i];

				cout << " " << hex << setfill('0') << setw(2) << (int)code_value << dec;
			}

			cout << " )" << endl;

			++extra_count;
		}
	}

	cout << "Number of extra rows is " << extra_count << endl << endl;
}

void Codec::PrintCompressMatrix()
{
	const int rows = _block_count;
	const int cols = _defer_count + _mix_count;

	cout << endl << "Compress matrix is " << rows << " x " << cols << ":" << endl;

	for (int ii = 0; ii < rows; ++ii)
	{
		for (int jj = 0; jj < cols; ++jj)
		{
			if (_compress_matrix[_ge_pitch * ii + (jj >> 6)] & ((u64)1 << (jj & 63)))
				cout << '1';
			else
				cout << '0';
		}
		cout << endl;
	}

	cout << endl;
}

void Codec::PrintPeeled()
{
	cout << "Peeled elements :";

	u16 row_i = _peel_head_rows;
	while (row_i != LIST_TERM)
	{
		PeelRow *row = &_peel_rows[row_i];

		cout << " " << row_i << "x" << row->peel_column;

		row_i = row->next;
	}

	cout << endl;
}

void Codec::PrintDeferredRows()
{
	cout << "Deferred rows :";

	u16 row_i = _defer_head_rows;
	while (row_i != LIST_TERM)
	{
		PeelRow *row = &_peel_rows[row_i];

		cout << " " << row_i;

		row_i = row->next;
	}

	cout << endl;
}

void Codec::PrintDeferredColumns()
{
	cout << "Deferred columns :";

	u16 column_i = _defer_head_columns;
	while (column_i != LIST_TERM)
	{
		PeelColumn *column = &_peel_cols[column_i];

		cout << " " << column_i;

		column_i = column->next;
	}

	cout << endl;
}

#endif // CAT_DUMP_CODEC_DEBUG


//// Encoder Mode

Result Codec::InitializeEncoder(int message_bytes, int block_bytes)
{
	Result r = ChooseMatrix(message_bytes, block_bytes);
	if (!r)
	{
		// Calculate partial final bytes
		u32 partial_final_bytes = message_bytes % _block_bytes;
		if (partial_final_bytes <= 0) partial_final_bytes = _block_bytes;

		// Encoder-specific
		_input_final_bytes = partial_final_bytes;
		_output_final_bytes = _block_bytes;
		_extra_count = 0;

		if (!AllocateWorkspace())
			r = R_OUT_OF_MEMORY;
	}

	return r;
}

/*
	EncodeFeed

		This function breaks the input message up into blocks
	and opportunistically peels with each one.  After processing
	all of the blocks from the input, it runs the matrix solver
	and if the solver succeeds, it generates the recovery blocks.

		In practice, the solver should always succeed because the
	encoder should be looking up its check matrix parameters from
	a table, which guarantees the matrix is invertible.
*/

Result Codec::EncodeFeed(const void *message_in)
{
	CAT_IF_DUMP(cout << endl << "---- EncodeFeed ----" << endl << endl;)

	// Validate input
	if CAT_UNLIKELY(message_in == 0) return R_BAD_INPUT;

	SetInput(message_in);

	// For each input row,
	for (u16 id = 0; id < _block_count; ++id)
	{
		if (!OpportunisticPeeling(id, id))
			return R_BAD_PEEL_SEED;
	}

	// Solve matrix and generate recovery blocks
	Result r = SolveMatrix();
	if (!r) GenerateRecoveryBlocks();
	else if (r == R_MORE_BLOCKS) r = R_BAD_PEEL_SEED;
	return r;
}

/*
	Encode

		This function encodes a block.  For the first N blocks
	it simply copies the input to the output block.  For other
	block identifiers, it will generate a new random row and
	sum together recovery blocks to produce the new block.
*/

u32 Codec::Encode(u32 id, void *block_out)
{
	if (!block_out) return 0;
	u8 * CAT_RESTRICT block = reinterpret_cast<u8 *>( block_out );

#if defined(CAT_COPY_FIRST_N)
	// For the message blocks,
	if (id < _block_count)
	{
		// Until the final block in message blocks,
		const u8 * CAT_RESTRICT src = _input_blocks + _block_bytes * id;
		if ((int)id == _block_count - 1)
		{
			// For the final block, copy partial block
			memcpy(block, src, _input_final_bytes);
			return _input_final_bytes;
		}

		// Copy from the original file data
		memcpy(block, src, _block_bytes);
		return _block_bytes;
	}
#endif // CAT_COPY_FIRST_N

	CAT_IF_DUMP(cout << "Encode: Generating row " << id << ":";)

	u16 peel_weight, peel_a, peel_x, mix_a, mix_x;
	GeneratePeelRow(id, _p_seed, _block_count, _mix_count,
		peel_weight, peel_a, peel_x, mix_a, mix_x);

	// Remember first column (there is always at least one)
	u8 * CAT_RESTRICT first = _recovery_blocks + _block_bytes * peel_x;

	CAT_IF_DUMP(cout << " " << peel_x;)

	// If peeler has multiple columns,
	if (peel_weight > 1)
	{
		--peel_weight;

		IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

		CAT_IF_DUMP(cout << " " << peel_x;)

		// Combine first two columns into output buffer (faster than memcpy + memxor)
		memxor_set(block, first, _recovery_blocks + _block_bytes * peel_x, _block_bytes);

		// For each remaining peeler column,
		while (--peel_weight > 0)
		{
			IterateNextColumn(peel_x, _block_count, _block_next_prime, peel_a);

			CAT_IF_DUMP(cout << " " << peel_x;)

			// Mix in each column
			memxor(block, _recovery_blocks + _block_bytes * peel_x, _block_bytes);
		}

		// Mix first mixer block in directly
		memxor(block, _recovery_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	}
	else
	{
		// Mix first with first mixer block (faster than memcpy + memxor)
		memxor_set(block, first, _recovery_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	}

	CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

	// For each remaining mixer column,
	IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
	memxor(block, _recovery_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

	IterateNextColumn(mix_x, _mix_count, _mix_next_prime, mix_a);
	memxor(block, _recovery_blocks + _block_bytes * (_block_count + mix_x), _block_bytes);
	CAT_IF_DUMP(cout << " " << (_block_count + mix_x);)

	CAT_IF_DUMP(cout << endl;)

	return _block_bytes;
}


//// Decoder Mode

Result Codec::InitializeDecoder(int message_bytes, int block_bytes)
{
	Result r = ChooseMatrix(message_bytes, block_bytes);
	if (r == R_WIN)
	{
		// Calculate partial final bytes
		u32 partial_final_bytes = message_bytes % _block_bytes;
		if (partial_final_bytes <= 0) partial_final_bytes = _block_bytes;

		// Decoder-specific
		_row_count = 0;
		_output_final_bytes = partial_final_bytes;

		// Hack: Prevents row-based ids from causing partial copies when they happen to be the last block id
		// This is only an issue because the shared codec source code happens to be built with more encoder-like semantics
		_input_final_bytes = _block_bytes;

		_extra_count = CAT_MAX_EXTRA_ROWS;
#if defined(CAT_ALL_ORIGINAL)
		_all_original = true;
#endif

		if (!AllocateInput() || !AllocateWorkspace())
			return R_OUT_OF_MEMORY;
	}

	return r;
}

/*
	DecodeFeed

		This function accumulates the new block in a large input
	buffer.  As soon as N blocks are collected, the matrix solver
	is attempted.  After N blocks, ResumeSolveMatrix() is used.
*/

Result Codec::DecodeFeed(u32 id, const void *block_in)
{
	// Validate input
	if CAT_UNLIKELY(block_in == 0)
		return R_BAD_INPUT;

	// If less than N rows stored,
	u16 row_i = _row_count;
	if (row_i < _block_count)
	{
#if defined(CAT_ALL_ORIGINAL)
		// If original data,
		if (id >= _block_count)
			_all_original = false;
#endif

		// If opportunistic peeling did not fail,
		if (OpportunisticPeeling(row_i, id))
		{
			u8 *block_store = _input_blocks + _block_bytes * row_i;

			// If this is the last block id,
			if (id == _block_count - 1)
			{
				u32 final_bytes = _output_final_bytes;

				// Copy the new row data into the input block area
				memcpy(block_store, block_in, final_bytes);

				// Pad with zeroes
				memset(block_store + final_bytes, 0, _block_bytes - final_bytes);
			}
			else
			{
				// Copy the new row data into the input block area
				memcpy(block_store, block_in, _block_bytes);
			}

			// If just acquired N blocks,
			if (++_row_count == _block_count)
			{
#if defined(CAT_ALL_ORIGINAL)
				// If all original data,
				if (_all_original && IsAllOriginalData())
					return R_WIN;
#endif

				// Attempt to solve the matrix and generate recovery blocks
				Result r = SolveMatrix();
				if (!r) Codec::GenerateRecoveryBlocks();
				return r;
			}
		} // end if opportunistic peeling succeeded

		return R_MORE_BLOCKS;
	}

	// Resume GE from this row
	Result r = ResumeSolveMatrix(id, block_in);
	if (!r) Codec::GenerateRecoveryBlocks();
	return r;
}

