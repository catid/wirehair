# Wirehair
## Fast and Portable Erasure Code

Wirehair produces a stream of error correction blocks from a data source
using an erasure code.  When enough of these blocks are received,
the original data can be recovered.  As compared to other similar
libraries, an unlimited number of error correction blocks can be produced,
and much larger block counts are supported.

A simple C API is provided to make it easy to incorporate into existing
projects.


#### Usage




#### Benchmarks

##### on iMac (2.7 GHz Core i5-2500S Sandy Bridge, June 2011):

Turbo Boost is turned on for these computations.  The block size is 1300 bytes
simulating a file transfer protocol over UDP.

~~~
~~~


#### Performance Discussion

Note that the performance is best around N=1024 blocks, where it takes roughly
2.5 milliseconds to process about 1.4 MB of data.  The performance for larger N
is cut in half when N is near 10,000, and in half again when N is near 60,000,
because the data no longer fits easily in the processor's cache.


#### Details

Wirehair is an FEC codec used to improve reliability of data sent
over a Binary Erasure Channel (BEC) such as satellite Internet or WiFi.
The data to transmit over the BEC is encoded into equal-sized blocks.
When enough blocks are received at the other end of the channel, then
the original data can be recovered.

How many additional blocks are needed is random, though the number
of additional blocks is low and does not vary based on the size of the
data.  Typical overhead is 0.001 additional blocks, and the overhead
gets massively smaller as the block count increases.

Wirehair is designed to be competitive in performance with FEC based
on Vandermonde matrices such as zfec, while allowing far more than
256 file blocks (it works best at N=1024 blocks).  Furthermore, it
can generate a stream of output that goes on forever, which has good
recovery properties for any part of that stream that is received,
rather than a fixed number of output blocks, making it far more
flexible than other libraries.


#### Credits

This software was written entirely by myself ( Christopher A. Taylor <mrcatid@gmail.com> ).  If you
find it useful and would like to buy me a coffee, consider [tipping](https://www.gittip.com/catid/).

