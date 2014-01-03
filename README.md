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
>> OKAY! N=3(0.0039 MB) encoder.BeginEncode in 60 usec, 65 MB/s
N=3 decoder.Decode in 39.787 usec, 98.022 MB/s.  Average overhead = 0 (0 bytes)

>> OKAY! N=10(0.013 MB) encoder.BeginEncode in 138 usec, 94.2029 MB/s
N=10 decoder.Decode in 100.812 usec, 128.953 MB/s.  Average overhead = 0.008 (10.4 bytes)
>> OKAY! N=11(0.0143 MB) encoder.BeginEncode in 137 usec, 104.38 MB/s
N=11 decoder.Decode in 97.625 usec, 146.479 MB/s.  Average overhead = 0.006 (7.8 bytes)
>> OKAY! N=12(0.0156 MB) encoder.BeginEncode in 130 usec, 120 MB/s
N=12 decoder.Decode in 98.945 usec, 157.663 MB/s.  Average overhead = 0.015 (19.5 bytes)
>> OKAY! N=13(0.0169 MB) encoder.BeginEncode in 124 usec, 136.29 MB/s
N=13 decoder.Decode in 98.696 usec, 171.233 MB/s.  Average overhead = 0.009 (11.7 bytes)

>> OKAY! N=110(0.143 MB) encoder.BeginEncode in 348 usec, 410.92 MB/s
N=110 decoder.Decode in 289.96 usec, 493.171 MB/s.  Average overhead = 0.02 (26 bytes)
>> OKAY! N=111(0.1443 MB) encoder.BeginEncode in 326 usec, 442.638 MB/s
N=111 decoder.Decode in 304.357 usec, 474.114 MB/s.  Average overhead = 0.015 (19.5 bytes)
>> OKAY! N=112(0.1456 MB) encoder.BeginEncode in 419 usec, 347.494 MB/s
N=112 decoder.Decode in 312.172 usec, 466.41 MB/s.  Average overhead = 0.022 (28.6 bytes)
>> OKAY! N=113(0.1469 MB) encoder.BeginEncode in 322 usec, 456.211 MB/s
N=113 decoder.Decode in 302.759 usec, 485.204 MB/s.  Average overhead = 0.016 (20.8 bytes)

>> OKAY! N=1110(1.443 MB) encoder.BeginEncode in 2873 usec, 502.262 MB/s
N=1110 decoder.Decode in 2482.65 usec, 581.233 MB/s.  Average overhead = 0.018 (23.4 bytes)
>> OKAY! N=1111(1.4443 MB) encoder.BeginEncode in 2151 usec, 671.455 MB/s
N=1111 decoder.Decode in 2480.74 usec, 582.206 MB/s.  Average overhead = 0.022 (28.6 bytes)
>> OKAY! N=1112(1.4456 MB) encoder.BeginEncode in 2905 usec, 497.625 MB/s
N=1112 decoder.Decode in 2559.12 usec, 564.881 MB/s.  Average overhead = 0.018 (23.4 bytes)
>> OKAY! N=1113(1.4469 MB) encoder.BeginEncode in 2203 usec, 656.786 MB/s
N=1113 decoder.Decode in 2455.72 usec, 589.195 MB/s.  Average overhead = 0.021 (27.3 bytes)

>> OKAY! N=11111(14.4443 MB) encoder.BeginEncode in 42666 usec, 338.544 MB/s
N=11111 decoder.Decode in 46541.7 usec, 310.352 MB/s.  Average overhead = 0.02 (26 bytes)
~~~


##### Performance Discussion

Note that the performance is best around N=1024 blocks, where it takes roughly
2.5 milliseconds to process about 1.4 MB of data.  The performance for larger N
is cut in half when N is near 10,000, and in half again when N is near 60,000,
because the data no longer fits easily in the processor's cache.


#### Credits

This software was written entirely by myself ( Christopher A. Taylor <mrcatid@gmail.com> ).  If you
find it useful and would like to buy me a coffee, consider [tipping](https://www.gittip.com/catid/).

