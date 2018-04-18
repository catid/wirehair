# Wirehair
## Fast and Portable Fountain Codes in C

Wirehair produces a stream of error correction blocks from a data source
using an erasure code.  When enough of these blocks are received,
the original data can be recovered.  As compared to other similar
libraries, an unlimited number of error correction blocks can be produced,
and much larger block counts are supported.

A simple C API is provided to make it easy to incorporate into existing
projects.


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
        cout << "!!! Wirehair nitialization failed: " << initResult << endl;
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

TBD: These need to be updated.  It's about 3x faster than it used to be, and the overhead is much lower.

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


#### Credits

Software by Christopher A. Taylor <mrcatid@gmail.com>

Please reach out if you need support or would like to collaborate on a project.
