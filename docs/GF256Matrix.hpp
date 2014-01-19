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

#ifndef CAT_WIREHAIR_GF256_MATRIX_HPP
#define CAT_WIREHAIR_GF256_MATRIX_HPP

#include "../codec_source/Platform.hpp"

namespace cat {

namespace wirehair {


/*
	This class allows me to generate random square invertible
	matrices separately from the Wirehair codec.  It's useful
	for testing error correcting properties of a code and for
	pregenerating invertible matrices.
*/
class GF256Matrix
{
	int _n;
	u8 *_matrix;
	int _pitch;
	u32 _seed;
	u16 *_pivot;

	void Cleanup();

public:
	GF256Matrix();
	~GF256Matrix();

	CAT_INLINE u8 *GetFront() { return _matrix; }
	CAT_INLINE int GetPitch() { return _pitch; }

	CAT_INLINE void SetSeed(u32 seed) { _seed = seed; }
	CAT_INLINE u32 GetSeed() { return _seed; }

	CAT_INLINE int Size() { return _n; }

	// Initializes an NxN random matrix and tries to invert it, returns false on error
	bool Initialize(int n);

	CAT_INLINE void NextSeed() { _seed++; }

	void Zero();
	void Identity();
	void Fill();
	bool Triangle();
	void Diagonal();

	void Print();
};


} // namespace wirehair

} // namespace cat

#endif // CAT_WIREHAIR_GF256_MATRIX_HPP
