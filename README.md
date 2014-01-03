# Wirehair FEC
## Fast Forward Error Correction (FEC) Codec

Wirehair FEC produces a stream of error correction blocks from a
data source.  When enough of these blocks are received, the original
data can be recovered.

Wirehair is an FEC codec used to improve reliability of data sent
over a Binary Erasure Channel (BEC) such as satellite Internet or WiFi.
The data to transmit over the BEC is encoded into equal-sized blocks.
When enough blocks are received at the other end of the channel, then
the original data can be recovered.

How many additional blocks are needed is random, though the number
of additional blocks is low and does not vary based on the size of the
data.  Typical overhead is between 0.015 and 0.03 additional blocks.
This may seem appalling to you, but think of it this way: The channel
is already lossy, so this is just like adding < 3% / N loss - a very
very small amount of additional loss (e.g. 0.003%) - to the channel,
down in the noise.

It is designed to be competitive in performance with FEC based on
Vandermonde matrices such as zfec, while allowing far more than 256
file blocks (it works best at N=1024 blocks).  Furthermore, it can
generate a stream of output that goes on forever, which has good
recovery properties for any part of that stream that is received,
rather than a fixed number of output blocks, making it far more
useful than zfec for file transfer or audio/video streaming.

Wirehair is released under the BSD license, which means that I ask
only that if you use my software that in the binary distribution of your
software you include the copyright notice in WIREHAIR.LICENSE and maybe
it would be nice to say thank you in an Email. :}


##### Future improvements:

+ Fix any seeds that have >3% overhead by extending the exception list.
+ C port.
+ Optimizations for very large files > 1GB.


##### Benchmarking on iMac (2.7 GHz Core i5-2500S Sandy Bridge, June 2011):

Turbo Boost is turned on for these computations.  The block size is 1300 bytes
simulating a file transfer protocol over UDP.

~~~
>> OKAY! N=3(0.0039 MB) encoder.BeginEncode in 77 usec, 50.6494 MB/s
N=3 decoder.Decode in 48.343 usec, 80.6735 MB/s.  Average overhead = 0 (0 bytes)

>> OKAY! N=10(0.013 MB) encoder.BeginEncode in 153 usec, 84.9673 MB/s
N=10 decoder.Decode in 112.664 usec, 115.387 MB/s.  Average overhead = 0.008 (10.4 bytes)
>> OKAY! N=11(0.0143 MB) encoder.BeginEncode in 160 usec, 89.375 MB/s
N=11 decoder.Decode in 115.56 usec, 123.745 MB/s.  Average overhead = 0.006 (7.8 bytes)
>> OKAY! N=12(0.0156 MB) encoder.BeginEncode in 163 usec, 95.7055 MB/s
N=12 decoder.Decode in 116.312 usec, 134.122 MB/s.  Average overhead = 0.015 (19.5 bytes)
>> OKAY! N=13(0.0169 MB) encoder.BeginEncode in 159 usec, 106.289 MB/s
N=13 decoder.Decode in 116.456 usec, 145.119 MB/s.  Average overhead = 0.009 (11.7 bytes)

>> OKAY! N=110(0.143 MB) encoder.BeginEncode in 391 usec, 365.729 MB/s
N=110 decoder.Decode in 301.612 usec, 474.119 MB/s.  Average overhead = 0.02 (26 bytes)
>> OKAY! N=111(0.1443 MB) encoder.BeginEncode in 344 usec, 419.477 MB/s
N=111 decoder.Decode in 318.903 usec, 452.489 MB/s.  Average overhead = 0.015 (19.5 bytes)

>> OKAY! N=1110(1.443 MB) encoder.BeginEncode in 2744 usec, 525.875 MB/s
N=1110 decoder.Decode in 2524.68 usec, 571.558 MB/s.  Average overhead = 0.018 (23.4 bytes)
>> OKAY! N=1111(1.4443 MB) encoder.BeginEncode in 2133 usec, 677.121 MB/s
N=1111 decoder.Decode in 2507.08 usec, 576.089 MB/s.  Average overhead = 0.022 (28.6 bytes)

>> OKAY! N=11111(14.4443 MB) encoder.BeginEncode in 44457 usec, 324.905 MB/s
N=11111 decoder.Decode in 46875.5 usec, 308.142 MB/s.  Average overhead = 0.02 (26 bytes)
~~~


##### Performance Discussion

Note that the performance is best around N=1024 blocks, where it takes roughly
2.5 milliseconds to process about 1.4 MB of data.  The performance for larger N
is cut in half when N is near 10,000, and in half again when N is near 60,000,
because the data no longer fits easily in the processor's cache.


#### Credits

This software was written entirely by myself ( Christopher A. Taylor <mrcatid@gmail.com> ).  If you
find it useful and would like to buy me a coffee, consider [tipping](https://www.gittip.com/catid/).

