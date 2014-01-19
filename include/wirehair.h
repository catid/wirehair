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

#ifndef CAT_WIREHAIR_HPP
#define CAT_WIREHAIR_HPP

#ifdef __cplusplus
extern "C" {
#endif

#define WIREHAIR_VERSION 2

/*
 * Verify binary compatibility with the Wirehair API on startup.
 *
 * Example:
 * 	if (wirehair_init()) throw "Update wirehair static library";
 *
 * Returns 0 on success.
 * Returns non-zero if the API level does not match.
 */
extern int _wirehair_init(int expected_version);
#define wirehair_init() _wirehair_init(WIREHAIR_VERSION)

typedef void *wirehair_state;


/*
 * Encode the given message into blocks of size block_bytes.
 *
 * The number of blocks in the message N = CEIL(bytes / block_bytes).
 *
 * Preconditions:
 *	block_bytes is a multiple of 2
 *
 * Returns 0 on success.
 * Returns non-zero on invalid input.
 */
extern int wirehair_encode(wirehair_state *E, const void *message, int bytes, int block_bytes);

/*
 * Returns the number of blocks N in the encoded message.
 */
extern int wirehair_count(wirehair_state E);

/*
 * Write an error correction block.
 *
 * The first id < N blocks are the same as the input data.  This can be
 * used to run the encoder in parallel with normal data transmission.
 *
 * Preconditions:
 *	block pointer has block_bytes of space available to store data
 *
 * Returns 0 on success.
 * Returns non-zero on invalid input.
 */
extern int wirehair_write(wirehair_state E, unsigned int id, void *block);

/*
 * Initialize a decoder for a message of size bytes with block_bytes bytes
 * per received block.
 *
 * Returns 0 on success.
 * Returns non-zero on invalid input.
 */
extern int wirehair_decode(wirehair_state *E, int bytes, int block_bytes);

/*
 * Feed a block to the decoder.
 *
 * This function will return 0 when decoding is possible.
 *
 * Preconditions:
 *	block pointer has block_bytes of space available to store data
 *
 * Returns 0 when decoding is likely to be possible.
 * Returns non-zero on invalid input or not enough data received.
 */
extern int wirehair_read(wirehair_state E, unsigned int id, const void *block);

/*
 * Attempt to reconstruct the message.
 *
 * This function will return 0 when decoding succeeds.
 *
 * Preconditions:
 *	message contains enough space to store the entire decoded message
 *
 * Returns 0 when decoding is complete.
 * Returns non-zero on invalid input or not enough data received.
 */
extern int wirehair_reconstruct(wirehair_state E, void *message);

/*
 * Free memory associated with a state object
 */
extern void wirehair_free(wirehair_state E);


#ifdef __cplusplus
}
#endif

#endif // CAT_WIREHAIR_HPP
