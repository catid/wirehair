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

#ifndef CAT_ABYSSINIAN_PRNG_HPP
#define CAT_ABYSSINIAN_PRNG_HPP

#include "Platform.hpp"

namespace cat {


/*
	This is a unified implementation of my favorite generator
	that is designed to generate up to 2^^32 numbers per seed.

	Its period is about 2^^126 and passes all BigCrush tests.
	It is the fastest generator I could find that passes all tests.

	Furthermore, the input seeds are hashed to avoid linear
	relationships between the input seeds and the low bits of
	the first few outputs.
*/
class CAT_EXPORT Abyssinian
{
	u64 _x, _y;

public:
	CAT_INLINE void Initialize(u32 x, u32 y)
	{
		// Based on the mixing functions of MurmurHash3
		static const u64 C1 = 0xff51afd7ed558ccdULL;
		static const u64 C2 = 0xc4ceb9fe1a85ec53ULL;

		x += y;
		y += x;

		u64 seed_x = 0x9368e53c2f6af274ULL ^ x;
		u64 seed_y = 0x586dcd208f7cd3fdULL ^ y;

		seed_x *= C1;
		seed_x ^= seed_x >> 33;
		seed_x *= C2;
		seed_x ^= seed_x >> 33;

		seed_y *= C1;
		seed_y ^= seed_y >> 33;
		seed_y *= C2;
		seed_y ^= seed_y >> 33;

		_x = seed_x;
		_y = seed_y;

		// Inlined Next(): Discard first output

		_x = (u64)0xfffd21a7 * (u32)_x + (u32)(_x >> 32);
		_y = (u64)0xfffd1361 * (u32)_y + (u32)(_y >> 32);
	}

	CAT_INLINE void Initialize(u32 seed)
	{
		Initialize(seed, seed);
	}

	CAT_INLINE u32 Next()
	{
		_x = (u64)0xfffd21a7 * (u32)_x + (u32)(_x >> 32);
		_y = (u64)0xfffd1361 * (u32)_y + (u32)(_y >> 32);
		return CAT_ROL32((u32)_x, 7) + (u32)_y;
	}
};


} // namespace cat

#endif // CAT_ABYSSINIAN_PRNG_HPP
