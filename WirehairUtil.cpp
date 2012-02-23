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

#include "WirehairUtil.hpp"
#include "SmallPRNG.hpp"
using namespace cat;
using namespace wirehair;


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

u16 cat::wirehair::SquareRoot16(u16 x)
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

u16 cat::wirehair::NextPrime16(u16 n)
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


//// Utility: Peeling Row Weight Generator function

static const u32 WEIGHT_DIST[] = {
	0, 5243, 529531, 704294, 791675, 844104, 879057, 904023,
	922747, 937311, 948962, 958494, 966438, 973160, 978921,
	983914, 988283, 992138, 995565, 998631, 1001391, 1003887,
	1006157, 1008229, 1010129, 1011876, 1013490, 1014983,
	1016370, 1017662, 1048576
};

u16 cat::wirehair::GeneratePeelRowWeight(u32 rv, u16 max_weight)
{
	rv &= 0xfffff;

	u16 ii;
	for (ii = 1; rv >= WEIGHT_DIST[ii]; ++ii);

	return ii > max_weight ? max_weight : ii;
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

	It's a little messy (will write random bits past the end of the
	matrix padded out to the end of the last word).  But the way I
	am using it, this is acceptable.  If you want to use it too be
	sure to verify this is not a problem for you...
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

bool cat::wirehair::AddInvertibleGF2Matrix(u64 *matrix, int offset, int pitch, int n)
{
	if (n <= 0) return false;

	// If we have this value of n in the table,
	if (n < 512)
	{
		// Pull a random matrix out of the lookup table
		CatsChoice prng;
		prng.Initialize(INVERTIBLE_MATRIX_SEEDS[n]);

		// If shift is friendly,
		u32 shift = offset & 63;
		u64 *row = matrix + (offset >> 6);
		if (shift > 0)
		{
			// For each row,
			for (int row_i = 0; row_i < n; ++row_i, row += pitch)
			{
				// For each word in the row,
				int add_pitch = (n + 63) / 64;
				u64 prev = 0;
				for (int ii = 0; ii < add_pitch; ++ii)
				{
					// Generate next word
					u32 rv1 = prng.Next();
					u32 rv2 = prng.Next();
					u64 word = ((u64)rv2 << 32) | rv1;

					// Add it in
					row[ii] ^= (prev >> (64 - shift)) | (word << shift);

					prev = word;
				}

				// Add last word if needed
				int last_bit = (shift + n + 63) / 64;
				if (last_bit > add_pitch)
					row[add_pitch] ^= prev >> (64 - shift);
			}
		}
		else // Rare aligned case:
		{
			// For each row,
			for (int row_i = 0; row_i < n; ++row_i, row += pitch)
			{
				// For each word in the row,
				for (int add_pitch = (n + 63) / 64, ii = 0; ii < add_pitch; ++ii)
				{
					// Generate next word
					u32 rv1 = prng.Next();
					u32 rv2 = prng.Next();
					u64 word = ((u64)rv2 << 32) | rv1;

					// Add it in
					row[ii] ^= word;
				}
			}
		}
	}
	else
	{
		// For each row,
		u64 *row = matrix;
		for (int ii = 0; ii < n; ++ii, row += pitch)
		{
			int column_i = offset + ii;

			// Flip diagonal bit
			u64 mask = (u64)1 << (column_i & 63);
			row[column_i >> 6] ^= mask;
		}
	}

	return true;
}


//// Utility: Deck Shuffling function

/*
	Given a PRNG, generate a deck of cards in a random order.
	The deck will contain elements with values between 0 and count - 1.
*/

void cat::wirehair::ShuffleDeck16(CatsChoice &prng, u16 *deck, u32 count)
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


//// Utility: Matrix Parameter Generator function

int g_p_seed, g_c_seed; // TODO: Remove these

bool cat::wirehair::GenerateMatrixParameters(int block_count, u32 &p_seed, u32 &c_seed, u16 &light_count, u16 &dense_count)
{
	p_seed = g_p_seed;// TODO: Remove these
	c_seed = g_c_seed;

	switch (block_count)
	{
	case 16:
		light_count = 6; dense_count = 0;
		return true;
	case 256:
		light_count = 16; dense_count = 4;
		return true;
	case 512:
		light_count = 14;
		dense_count = 5;
		return true;
	case 1024:
		light_count = 18;
		dense_count = 12;
		return true;
	case 2048:
		light_count = 45;
		dense_count = 8;
		return true;
	case 4096:
		light_count = 55;
		dense_count = 14;
		return true;
	case 8192:
		light_count = 100;
		dense_count = 16;
		return true;
	case 10000:
		light_count = 120;
		dense_count = 20;
		return true;
	case 16384:
		light_count = 180;
		dense_count = 26;
		return true;
	case 32768:
		light_count = 400;
		dense_count = 30;
		return true;
	case 40000:
		light_count = 460;
		dense_count = 29;
		return true;
	case 50000:
		light_count = 600;
		dense_count = 34;
		return true;
	case 64000:
		light_count = 6;
		dense_count = 750;
		return true;
	}

	return false;
}
