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
For N = 12 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 31 usec (503.226 MBPS)
+ Average wirehair_encode() time: 0 usec (7834.57 MBPS)
+ Average wirehair_decode() time: 2 usec (638.988 MBPS)
+ Average overhead piece count beyond N = 0.011
+ Average wirehair_recover() time: 1 usec (8144.09 MBPS)

For N = 32 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 42 usec (990.476 MBPS)
+ Average wirehair_encode() time: 0 usec (7269.13 MBPS)
+ Average wirehair_decode() time: 1 usec (920.302 MBPS)
+ Average overhead piece count beyond N = 0.0205
+ Average wirehair_recover() time: 5 usec (8061.23 MBPS)

For N = 102 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 131 usec (1012.21 MBPS)
+ Average wirehair_encode() time: 0 usec (6587.3 MBPS)
+ Average wirehair_decode() time: 1 usec (999.513 MBPS)
+ Average overhead piece count beyond N = 0.0195
+ Average wirehair_recover() time: 17 usec (7474 MBPS)

For N = 134 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 179 usec (973.184 MBPS)
+ Average wirehair_encode() time: 0 usec (5942.63 MBPS)
+ Average wirehair_decode() time: 1 usec (949.373 MBPS)
+ Average overhead piece count beyond N = 0.02
+ Average wirehair_recover() time: 24 usec (7178.03 MBPS)

For N = 169 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 212 usec (1036.32 MBPS)
+ Average wirehair_encode() time: 0 usec (6016.77 MBPS)
+ Average wirehair_decode() time: 1 usec (975.657 MBPS)
+ Average overhead piece count beyond N = 0.0225
+ Average wirehair_recover() time: 30 usec (7129.99 MBPS)

For N = 201 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 275 usec (950.182 MBPS)
+ Average wirehair_encode() time: 0 usec (6587.98 MBPS)
+ Average wirehair_decode() time: 1 usec (958.909 MBPS)
+ Average overhead piece count beyond N = 0.012
+ Average wirehair_recover() time: 40 usec (6436.36 MBPS)

For N = 294 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 355 usec (1076.62 MBPS)
+ Average wirehair_encode() time: 0 usec (6218.16 MBPS)
+ Average wirehair_decode() time: 1 usec (1049.73 MBPS)
+ Average overhead piece count beyond N = 0.021
+ Average wirehair_recover() time: 55 usec (6894.69 MBPS)

For N = 359 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 457 usec (1021.23 MBPS)
+ Average wirehair_encode() time: 0 usec (7052.25 MBPS)
+ Average wirehair_decode() time: 1 usec (1082.96 MBPS)
+ Average overhead piece count beyond N = 0.0165
+ Average wirehair_recover() time: 66 usec (7004.4 MBPS)

For N = 413 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 524 usec (1024.62 MBPS)
+ Average wirehair_encode() time: 0 usec (6426.07 MBPS)
+ Average wirehair_decode() time: 1 usec (1079.96 MBPS)
+ Average overhead piece count beyond N = 0.021
+ Average wirehair_recover() time: 76 usec (7044.54 MBPS)

For N = 770 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 1201 usec (833.472 MBPS)
+ Average wirehair_encode() time: 0 usec (7265.59 MBPS)
+ Average wirehair_decode() time: 1 usec (865.626 MBPS)
+ Average overhead piece count beyond N = 0.018
+ Average wirehair_recover() time: 144 usec (6927.41 MBPS)

For N = 1000 packets of 1300 bytes:
+ Average wirehair_encoder_create() time: 1564 usec (831.202 MBPS)
+ Average wirehair_encode() time: 0 usec (5314.48 MBPS)
+ Average wirehair_decode() time: 1 usec (696.552 MBPS)
+ Average overhead piece count beyond N = 0.019
+ Average wirehair_recover() time: 195 usec (6653.45 MBPS)
~~~


#### Credits

Software by Christopher A. Taylor <mrcatid@gmail.com>

Please reach out if you need support or would like to collaborate on a project.
