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

#ifndef CAT_WIREHAIR_UTIL_HPP
#define CAT_WIREHAIR_UTIL_HPP

#include "Platform.hpp"

extern int g_seed;

namespace cat {

namespace wirehair {


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
bool GenerateMatrixParameters(int block_count, u32 &seed, u16 &light_count, u16 &dense_count);


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


} // namespace wirehair

} // namespace cat

#endif // CAT_WIREHAIR_UTIL_HPP
