# Wirehair
## Fast and Portable Fountain Codes in C

Wirehair produces a stream of error correction blocks from a data source
using an erasure code.  When enough of these blocks are received,
the original data can be recovered.

As compared to other similar libraries, an unlimited number of error
correction blocks can be produced, and much larger block counts are supported.
Furthermore, it gets slower as O(N) in the amount of input data rather
than O(N Log N) like the Leopard block code or O(N^2) like the Fecal fountain code,
so it is well-suited for large data.

This is not an ideal MDS code, so sometimes it will fail to recover N
original data packets from N symbol packets.  It may take N + 1 or N + 2 or more.
On average it takes about N + 0.02 packets to recover.  Overall the overhead
from the code inefficiency is low, compared to LDPC and many other fountain codes.

A simple C API is provided to make it easy to incorporate into existing
projects.  No external dependencies are required.


##### Building: Quick Setup

The source code in this folder (gf256 and wirehair code) can be incorporated
into your project without any other external dependencies.

The proj folder contains Visual Studio project files for Windows builds.
There is a CMakeLists file to build it on other platforms.


#### Example Usage

Here's an example program using Wirehair.  It's included in the UnitTest project and demonstrates both the sender and receiver, which are normally separate programs.  For example the data sender might be a file server and the data receiver might be downloading a file from the sender.

~~~
static bool ReadmeExample()
{
    // Size of packets to produce
    static const int kPacketSize = 1400;

    // Note: Does not need to be an even multiple of packet size or 16 etc
    static const int kMessageBytes = 1000 * 1000 + 333;

    vector<uint8_t> message(kMessageBytes);

    // Fill message contents
    memset(&message[0], 1, message.size());

    // Create encoder
    WirehairCodec encoder = wirehair_encoder_create(nullptr, &message[0], kMessageBytes, kPacketSize);
    if (!encoder)
    {
        cout << "!!! Failed to create encoder" << endl;
        return false;
    }

    // Create decoder
    WirehairCodec decoder = wirehair_decoder_create(nullptr, kMessageBytes, kPacketSize);
    if (!decoder)
    {
        // Free memory for encoder
        wirehair_free(encoder);

        cout << "!!! Failed to create decoder" << endl;
        return false;
    }

    unsigned blockId = 0, needed = 0;

    for (;;)
    {
        // Select which block to encode.
        // Note: First N blocks are the original data, so it's possible to start
        // sending data while wirehair_encoder_create() is getting started.
        blockId++;

        // Simulate 10% packetloss
        if (blockId % 10 == 0) {
            continue;
        }

        // Keep track of how many pieces were needed
        ++needed;

        vector<uint8_t> block(kPacketSize);

        // Encode a packet
        uint32_t writeLen = 0;
        WirehairResult encodeResult = wirehair_encode(
            encoder, // Encoder object
            blockId, // ID of block to generate
            &block[0], // Output buffer
            kPacketSize, // Output buffer size
            &writeLen); // Returned block length

        if (encodeResult != Wirehair_Success)
        {
            cout << "wirehair_encode failed: " << encodeResult << endl;
            return false;
        }

        // Attempt decode
        WirehairResult decodeResult = wirehair_decode(
            decoder, // Decoder object
            blockId, // ID of block that was encoded
            &block[0], // Input block
            writeLen); // Block length

        // If decoder returns success:
        if (decodeResult == Wirehair_Success) {
            // Decoder has enough data to recover now
            break;
        }

        if (decodeResult != Wirehair_NeedMore)
        {
            cout << "wirehair_decode failed: " << decodeResult << endl;
            return false;
        }
    }

    vector<uint8_t> decoded(kMessageBytes);

    // Recover original data on decoder side
    WirehairResult decodeResult = wirehair_recover(
        decoder,
        &decoded[0],
        kMessageBytes);

    if (decodeResult != Wirehair_Success)
    {
        cout << "wirehair_recover failed: " << decodeResult << endl;
        return false;
    }

    // Free memory for encoder and decoder
    wirehair_free(encoder);
    wirehair_free(decoder);

    return true;
}

int main()
{
    const WirehairResult initResult = wirehair_init();

    if (initResult != Wirehair_Success)
    {
        SIAMESE_DEBUG_BREAK();
        cout << "!!! Wirehair initialization failed: " << initResult << endl;
        return -1;
    }

    if (!ReadmeExample())
    {
        SIAMESE_DEBUG_BREAK();
        cout << "!!! Example usage failed" << endl;
        return -2;
    }
...
~~~


#### Benchmarks

Some quick comments:

Benchmarks on my PC do not mean a whole lot.  Right now it's clocked at 3 GHz and has Turbo Boost on, etc.
To run the test yourself just build and run the UnitTest project in Release mode.

For small values of N < 128 or so this is a pretty inefficient codec compared to the Fecal codec.  Fecal is also a fountain code but is limited to repairing a small number of failures or small input block count.

~~~
For N = 2 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 11 usec (236.364 MBPS)
+ Average wirehair_encode() time: 0 usec (8673.84 MBPS)
+ Average wirehair_decode() time: 2 usec (482.724 MBPS)
+ Average overhead piece count beyond N = 0.015
+ Average wirehair_recover() time: 0 usec (9701.49 MBPS)

For N = 4 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 8 usec (650 MBPS)
+ Average wirehair_encode() time: 0 usec (8354.32 MBPS)
+ Average wirehair_decode() time: 1 usec (708.281 MBPS)
+ Average overhead piece count beyond N = 0.0165
+ Average wirehair_recover() time: 0 usec (10547.7 MBPS)

For N = 8 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 13 usec (800 MBPS)
+ Average wirehair_encode() time: 0 usec (8852.28 MBPS)
+ Average wirehair_decode() time: 1 usec (707.544 MBPS)
+ Average overhead piece count beyond N = 0.0045
+ Average wirehair_recover() time: 1 usec (9130.82 MBPS)

For N = 16 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 26 usec (800 MBPS)
+ Average wirehair_encode() time: 0 usec (9029.06 MBPS)
+ Average wirehair_decode() time: 1 usec (715.655 MBPS)
+ Average overhead piece count beyond N = 0.037
+ Average wirehair_recover() time: 2 usec (9363.04 MBPS)

For N = 32 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 40 usec (1040 MBPS)
+ Average wirehair_encode() time: 0 usec (8178.93 MBPS)
+ Average wirehair_decode() time: 1 usec (934.853 MBPS)
+ Average overhead piece count beyond N = 0.0205
+ Average wirehair_recover() time: 5 usec (8192.2 MBPS)

For N = 64 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 83 usec (1002.41 MBPS)
+ Average wirehair_encode() time: 0 usec (7340.91 MBPS)
+ Average wirehair_decode() time: 1 usec (1004.77 MBPS)
+ Average overhead piece count beyond N = 0.024
+ Average wirehair_recover() time: 11 usec (7409.06 MBPS)

For N = 128 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 196 usec (848.98 MBPS)
+ Average wirehair_encode() time: 0 usec (6068.85 MBPS)
+ Average wirehair_decode() time: 1 usec (854.2 MBPS)
+ Average overhead piece count beyond N = 0.02
+ Average wirehair_recover() time: 26 usec (6275.69 MBPS)

For N = 256 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 330 usec (1008.48 MBPS)
+ Average wirehair_encode() time: 0 usec (6709.24 MBPS)
+ Average wirehair_decode() time: 1 usec (994.742 MBPS)
+ Average overhead piece count beyond N = 0.0245
+ Average wirehair_recover() time: 51 usec (6491.95 MBPS)

For N = 512 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 677 usec (983.161 MBPS)
+ Average wirehair_encode() time: 0 usec (6701.32 MBPS)
+ Average wirehair_decode() time: 1 usec (978.883 MBPS)
+ Average overhead piece count beyond N = 0.0245
+ Average wirehair_recover() time: 104 usec (6343.52 MBPS)

For N = 1024 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 1723 usec (772.606 MBPS)
+ Average wirehair_encode() time: 0 usec (5463.55 MBPS)
+ Average wirehair_decode() time: 2 usec (637.944 MBPS)
+ Average overhead piece count beyond N = 0.017
+ Average wirehair_recover() time: 218 usec (6085.25 MBPS)

For N = 2048 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 3307 usec (805.08 MBPS)
+ Average wirehair_encode() time: 0 usec (5199.33 MBPS)
+ Average wirehair_decode() time: 1 usec (650.556 MBPS)
+ Average overhead piece count beyond N = 0.026
+ Average wirehair_recover() time: 448 usec (5929.76 MBPS)

For N = 4096 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 7836 usec (679.53 MBPS)
+ Average wirehair_encode() time: 0 usec (4447.8 MBPS)
+ Average wirehair_decode() time: 2 usec (561.758 MBPS)
+ Average overhead piece count beyond N = 0.013
+ Average wirehair_recover() time: 1239 usec (4295.54 MBPS)

For N = 8192 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 17524 usec (607.715 MBPS)
+ Average wirehair_encode() time: 0 usec (3366.69 MBPS)
+ Average wirehair_decode() time: 2 usec (506.866 MBPS)
+ Average overhead piece count beyond N = 0.0165
+ Average wirehair_recover() time: 2995 usec (3555.25 MBPS)

For N = 16384 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 43180 usec (493.265 MBPS)
+ Average wirehair_encode() time: 0 usec (2629.7 MBPS)
+ Average wirehair_decode() time: 3 usec (421.17 MBPS)
+ Average overhead piece count beyond N = 0.0395
+ Average wirehair_recover() time: 7608 usec (2799.48 MBPS)
~~~


#### Credits

Software by Christopher A. Taylor <mrcatid@gmail.com>

Please reach out if you need support or would like to collaborate on a project.
