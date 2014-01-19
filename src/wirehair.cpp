/*
	Copyright (c) 2012-2014 Christopher A. Taylor.  All rights reserved.

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

#include "wirehair.h"
#include "wirehair_codec.hpp"
using namespace cat;
using namespace wirehair;

static bool m_init = false;

int _wirehair_init(int expected_version) {
	// If version mismatch,
	if (expected_version != WIREHAIR_VERSION) {
		return -1;
	}

	m_init = true;

	return 0;
}

int wirehair_encode(wirehair_state *E, const void *message, int bytes, int block_bytes) {
	// If input is invalid,
	if CAT_UNLIKELY(!m_init || !E || !message || bytes < 1 ||
					block_bytes < 1 || block_bytes % 2 == 0) {
		return -1;
	}

	Codec *codec = reinterpret_cast<Codec *>( *E );

	// Allocate a new Codec object
	if (!codec) {
		codec = new Codec;
		*E = codec;
	}

	// Initialize codec
	Result r = codec->InitializeEncoder(bytes, block_bytes);

	if (!r) {
		// Feed message to codec
		r = codec->EncodeFeed(message);
	}

	return r;
}

int wirehair_count(wirehair_state E) {
	// If input is invalid,
	if CAT_UNLIKELY(!E) {
		return 0;
	}

	Codec *codec = reinterpret_cast<Codec *>( E );

	return codec->BlockCount();
}

int wirehair_write(wirehair_state E, unsigned int id, void *block) {
	// If input is invalid,
	if CAT_UNLIKELY(!E || !block) {
		return -1;
	}

	Codec *codec = reinterpret_cast<Codec *>( E );

	return codec->Encode(id, block);
}

int wirehair_decode(wirehair_state *E, int bytes, int block_bytes) {
	// If input is invalid,
	if CAT_UNLIKELY(!E || bytes < 1 || block_bytes < 1 ||
					block_bytes % 2 == 0) {
		return -1;
	}

	Codec *codec = reinterpret_cast<Codec *>( *E );

	// Allocate a new Codec object
	if (!codec) {
		codec = new Codec;
		*E = codec;
	}

	// Allocate memory for decoding
	return codec->InitializeDecoder(bytes, block_bytes);
}

int wirehair_read(wirehair_state E, unsigned int id, const void *block) {
	// If input is invalid,
	if CAT_UNLIKELY(!E || !block) {
		return -1;
	}

	Codec *codec = reinterpret_cast<Codec *>( E );

	return codec->DecodeFeed(id, block);
}

int wirehair_reconstruct(wirehair_state E, void *message) {
	// If input is invalid,
	if CAT_UNLIKELY(!E || !message) {
		return -1;
	}

	Codec *codec = reinterpret_cast<Codec *>( E );

	return codec->ReconstructOutput(message);
}

void wirehair_free(wirehair_state E) {
	Codec *codec = reinterpret_cast<Codec *>( E );

	if (codec) {
		delete codec;
	}
}

