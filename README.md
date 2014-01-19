# Wirehair
## Fast and Portable Erasure Codes in C

Wirehair produces a stream of error correction blocks from a data source
using an erasure code.  When enough of these blocks are received,
the original data can be recovered.  As compared to other similar
libraries, an unlimited number of error correction blocks can be produced,
and much larger block counts are supported.

A simple C API is provided to make it easy to incorporate into existing
projects.


##### Building: Quick Setup

The [wirehair-mobile](https://github.com/catid/wirehair/tree/master/wirehair-mobile)
directory contains an easy-to-import set of C code that also
builds properly for mobile devices.


#### Usage: Encoder

On startup, verify that the Wirehair library is linked correctly:

~~~
	if (!wirehair_init()) {
		// Wrong wirehair static library linked
		exit(1);
	}
~~~

Create an encoder object, providing the message to encode, the number of
bytes in the message, and the number of bytes in each block.  The message
will get broken into equal-sized blocks.

~~~
	char *message = ...;
	int bytes = 1000000;
	int block_bytes = 1300;

	wirehair_state encoder = 0;

	encoder = wirehair_encode(0, message, bytes, block_bytes);
	assert(encoder);

	// The encoder object can now be used to write blocks
~~~

To check how many blocks are in the message:

~~~
	int N = wirehair_count(encoder);

	// N = number of blocks in message
~~~

Each time a block will be written it should be assigned a unique ID number,
starting from 0 and incremented by one each time:

~~~
	char block[1300];
	int ID = 0;

	if (!wirehair_write(encoder, ID, block)) {
		exit(1);
	}
~~~

When you are done with the encoder, you can either free the encoder object
to reclaim the memory, or reuse the encoder object again.  To reuse the
object, pass it as the first argument to `wirehair_encode`.  To free the
object:

~~~
	wirehair_free(encoder);
~~~


#### Usage: Decoder

On startup, verify that the Wirehair library is linked correctly as in
the encoder case.

A decoder object should be created by providing the size of the message
and the number of bytes per block:

~~~
	int bytes = 1000000;
	int block_bytes = 1300;

	wirehair_state decoder = 0;

	decoder = wirehair_decode(0, bytes, block_bytes);
	assert(decoder);

	// The decoder object can now be used to decode the message
~~~

To check how many blocks are in the message:

~~~
	int N = wirehair_count(decoder);

	// N = number of blocks in message
~~~

To decode a message one received block at a time:

~~~
	char *message = new char[bytes];
	int ID = 0;

	// If there is a chance of decoding the message now,
	if (wirehair_read(decoder, ID, block)) {

		// If message can be reconstructed,
		if (wirehair_reconstruct(decoder, message)) {

			// Decoding message complete

		}
	}
~~~

Similar to the encoder, you may either free or reuse the decoder object
in the exact same way as the encoder.


#### Benchmarks

##### libwirehair.a on Macbook Air (1.7 GHz Core i5-2557M Sandy Bridge, July 2011):

Turbo Boost is turned on for these computations.  The block size is 1300 bytes
simulating a file transfer protocol over UDP.

~~~
wirehair_encode(N = 12) in 150 usec, 104 MB/s
wirehair_decode(N = 12) average overhead = 0.009 blocks, average reconstruct time = 134.786 usec, 115.739 MB/s

wirehair_encode(N = 32) in 198 usec, 210.101 MB/s
wirehair_decode(N = 32) average overhead = 0.027 blocks, average reconstruct time = 190.22 usec, 218.694 MB/s

wirehair_encode(N = 102) in 310 usec, 427.742 MB/s
wirehair_decode(N = 102) average overhead = 0.016 blocks, average reconstruct time = 378.212 usec, 350.597 MB/s

wirehair_encode(N = 134) in 416 usec, 418.75 MB/s
wirehair_decode(N = 134) average overhead = 0.023 blocks, average reconstruct time = 497.412 usec, 350.213 MB/s

wirehair_encode(N = 169) in 525 usec, 418.476 MB/s
wirehair_decode(N = 169) average overhead = 0.022 blocks, average reconstruct time = 712.778 usec, 308.231 MB/s

wirehair_encode(N = 201) in 654 usec, 399.541 MB/s
wirehair_decode(N = 201) average overhead = 0.019 blocks, average reconstruct time = 798.307 usec, 327.318 MB/s

wirehair_encode(N = 294) in 852 usec, 448.592 MB/s
wirehair_decode(N = 294) average overhead = 0.018 blocks, average reconstruct time = 1048.08 usec, 364.668 MB/s

wirehair_encode(N = 359) in 1176 usec, 396.854 MB/s
wirehair_decode(N = 359) average overhead = 0.021 blocks, average reconstruct time = 1210.73 usec, 385.471 MB/s

wirehair_encode(N = 413) in 1128 usec, 475.975 MB/s
wirehair_decode(N = 413) average overhead = 0.016 blocks, average reconstruct time = 1538.1 usec, 349.067 MB/s
~~~


#### Details

Wirehair is an FEC codec used to improve reliability of data sent
over a Binary Erasure Channel (BEC) such as satellite Internet or WiFi.
The data to transmit over the BEC is encoded into equal-sized blocks.
When enough blocks are received at the other end of the channel, then
the original data can be recovered.

How many additional blocks are needed is random, though the number
of additional blocks is low and does not vary based on the size of the
data.  Typical overhead is 0.03 additional blocks.

Wirehair is designed to be competitive in performance with FEC based
on Vandermonde matrices such as zfec, while allowing far more than
256 file blocks (it works best at N=1024 blocks).  Furthermore, it
can generate a stream of output that goes on forever, which has good
recovery properties for any part of that stream that is received,
rather than a fixed number of output blocks, making it far more
flexible than other libraries.


##### Discussion: Overhead Reductions with GF(2^16)

Wirehair uses a random 6x18 GF(256) matrix to achieve its low 3%
recovery failure rate despite doing most of the calculations on a
large sparse binary matrix.

It is trivially possible to decrease the overhead to 0.1% or less
by switching to a 9x20 GF(2^16) matrix.  I've explored this option
and found that it is roughly twice as slow for small block counts,
but similar in performance for N >= 1000, since the GF(2^16) matrix
does not grow larger.

For maximum MTU UDP file transfer of roughly 1300 bytes, the average
extra bytes needing to be transmitted due to the <3% failure rate is
under 40 bytes.  Paying a steep performance penalty to reduce this
overhead to 4 bytes seems like it is not worth it.  For moderately
sized messages, 40 bytes is dwarfed by other sources of overhead.


#### Credits

This software was written entirely by myself ( Christopher A. Taylor <mrcatid@gmail.com> ).  If you
find it useful and would like to buy me a coffee, consider [tipping](https://www.gittip.com/catid/).

