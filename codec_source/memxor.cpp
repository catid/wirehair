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

#include "memxor.hpp"
using namespace cat;


void cat::memxor(void *voutput, const void *vinput, int bytes)
{
	/*
		Often times the output is XOR'd in-place so this version is
		faster than the one below with two inputs.

		There would be a decent performance improvement by using MMX
		if the buffers are all aligned.  However, I don't really have
		control over how the buffers are aligned because the block
		sizes are not always multiples of 16 bytes.  So after an hour
		of tuning this is the best version I found.
	*/

	// Primary engine
	u64 *output64 = reinterpret_cast<u64*>( voutput );
	const u64 *input64 = reinterpret_cast<const u64*>( vinput );

	while (bytes >= 128)
	{
		output64[0] ^= input64[0];
		output64[1] ^= input64[1];
		output64[2] ^= input64[2];
		output64[3] ^= input64[3];
		output64[4] ^= input64[4];
		output64[5] ^= input64[5];
		output64[6] ^= input64[6];
		output64[7] ^= input64[7];
		output64[8] ^= input64[8];
		output64[9] ^= input64[9];
		output64[10] ^= input64[10];
		output64[11] ^= input64[11];
		output64[12] ^= input64[12];
		output64[13] ^= input64[13];
		output64[14] ^= input64[14];
		output64[15] ^= input64[15];
		output64 += 16;
		input64 += 16;
		bytes -= 128;
	}

	// Handle remaining multiples of 8 bytes
	while (bytes >= 8)
	{
		*output64++ ^= *input64++;
		bytes -= 8;
	}

	// Handle final <8 bytes
	u8 *output = reinterpret_cast<u8*>( output64 );
	const u8 *input = reinterpret_cast<const u8*>( input64 );

	switch (bytes)
	{
	case 7:	output[6] ^= input[6];
	case 6:	output[5] ^= input[5];
	case 5:	output[4] ^= input[4];
	case 4:	*(u32*)output ^= *(u32*)input;
		break;
	case 3:	output[2] ^= input[2];
	case 2:	output[1] ^= input[1];
	case 1:	output[0] ^= input[0];
	case 0:
	default:
		break;
	}
}

void cat::memxor_set(void *voutput, const void *va, const void *vb, int bytes)
{
	/*
		This version exists to avoid an expensive memory copy operation when
		an input block is being calculated from a row and some other blocks.
	*/

	// Primary engine
	u64 *output64 = reinterpret_cast<u64*>( voutput );
	const u64 *a64 = reinterpret_cast<const u64*>( va );
	const u64 *b64 = reinterpret_cast<const u64*>( vb );

	while (bytes >= 128)
	{
		output64[0] = a64[0] ^ b64[0];
		output64[1] = a64[1] ^ b64[1];
		output64[2] = a64[2] ^ b64[2];
		output64[3] = a64[3] ^ b64[3];
		output64[4] = a64[4] ^ b64[4];
		output64[5] = a64[5] ^ b64[5];
		output64[6] = a64[6] ^ b64[6];
		output64[7] = a64[7] ^ b64[7];
		output64[8] = a64[8] ^ b64[8];
		output64[9] = a64[9] ^ b64[9];
		output64[10] = a64[10] ^ b64[10];
		output64[11] = a64[11] ^ b64[11];
		output64[12] = a64[12] ^ b64[12];
		output64[13] = a64[13] ^ b64[13];
		output64[14] = a64[14] ^ b64[14];
		output64[15] = a64[15] ^ b64[15];
		output64 += 16;
		a64 += 16;
		b64 += 16;
		bytes -= 128;
	}

	// Handle remaining multiples of 8 bytes
	while (bytes >= 8)
	{
		*output64++ = *a64++ ^ *b64++;
		bytes -= 8;
	}

	// Handle final <8 bytes
	u8 *output = reinterpret_cast<u8*>( output64 );
	const u8 *a = reinterpret_cast<const u8*>( a64 );
	const u8 *b = reinterpret_cast<const u8*>( b64 );

	switch (bytes)
	{
	case 7:	output[6] = a[6] ^ b[6];
	case 6:	output[5] = a[5] ^ b[5];
	case 5:	output[4] = a[4] ^ b[4];
	case 4:	*(u32*)output = *(u32*)a ^ *(u32*)b;
		break;
	case 3:	output[2] = a[2] ^ b[2];
	case 2:	output[1] = a[1] ^ b[1];
	case 1:	output[0] = a[0] ^ b[0];
	case 0:
	default:
		break;
	}
}

void cat::memxor_add(void *voutput, const void *va, const void *vb, int bytes)
{
	/*
		This version adds to the output instead of overwriting it.
	*/

	// Primary engine
	u64 *output64 = reinterpret_cast<u64*>( voutput );
	const u64 *a64 = reinterpret_cast<const u64*>( va );
	const u64 *b64 = reinterpret_cast<const u64*>( vb );

	while (bytes >= 128)
	{
		output64[0] ^= a64[0] ^ b64[0];
		output64[1] ^= a64[1] ^ b64[1];
		output64[2] ^= a64[2] ^ b64[2];
		output64[3] ^= a64[3] ^ b64[3];
		output64[4] ^= a64[4] ^ b64[4];
		output64[5] ^= a64[5] ^ b64[5];
		output64[6] ^= a64[6] ^ b64[6];
		output64[7] ^= a64[7] ^ b64[7];
		output64[8] ^= a64[8] ^ b64[8];
		output64[9] ^= a64[9] ^ b64[9];
		output64[10] ^= a64[10] ^ b64[10];
		output64[11] ^= a64[11] ^ b64[11];
		output64[12] ^= a64[12] ^ b64[12];
		output64[13] ^= a64[13] ^ b64[13];
		output64[14] ^= a64[14] ^ b64[14];
		output64[15] ^= a64[15] ^ b64[15];
		output64 += 16;
		a64 += 16;
		b64 += 16;
		bytes -= 128;
	}

	// Handle remaining multiples of 8 bytes
	while (bytes >= 8)
	{
		*output64++ ^= *a64++ ^ *b64++;
		bytes -= 8;
	}

	// Handle final <8 bytes
	u8 *output = reinterpret_cast<u8*>( output64 );
	const u8 *a = reinterpret_cast<const u8*>( a64 );
	const u8 *b = reinterpret_cast<const u8*>( b64 );

	switch (bytes)
	{
	case 7:	output[6] ^= a[6] ^ b[6];
	case 6:	output[5] ^= a[5] ^ b[5];
	case 5:	output[4] ^= a[4] ^ b[4];
	case 4:	*(u32*)output ^= *(u32*)a ^ *(u32*)b;
		break;
	case 3:	output[2] ^= a[2] ^ b[2];
	case 2:	output[1] ^= a[1] ^ b[1];
	case 1:	output[0] ^= a[0] ^ b[0];
	case 0:
	default:
		break;
	}
}
