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
using namespace cat;
using namespace wirehair;


//// PeelingMatrixGenerator

static const u32 WEIGHT_DIST[] = {
	0, 5243, 529531, 704294, 791675, 844104, 879057, 904023,
	922747, 937311, 948962, 958494, 966438, 973160, 978921,
	983914, 988283, 992138, 995565, 998631, 1001391, 1003887,
	1006157, 1008229, 1010129, 1011876, 1013490, 1014983,
	1016370, 1017662, 1048576
};

u16 PMMatrixGenerator::StartWeight(u32 seed, u32 row)
{
	CatsChoice prng;

	prng.Initialize(seed, row);

	_x = _x0 = prng.Next() % _kp;
	_a = (prng.Next() % (_kp - 1)) + 1;

	u32 n = prng.Next() & 0xfffff;

	int ii, degree_max = _max_degree < 40 ? _max_degree : 40;

	for (ii = 1; ii < degree_max && n >= WEIGHT_DIST[ii]; ++ii);

	return (u16)ii;
}


//// Encoder

bool Encoder::Initialize(const void *message_in, int message_bytes, int block_bytes)
{
}

void Encoder::Generate(void *block)
{
}


//// Decoder

bool Decoder::Initialize(void *message_out, int message_bytes, int block_bytes)
{
}

bool Decoder::Decode(void *block)
{
}
