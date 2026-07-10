/*
    Copyright (c) 2012-2018 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Wirehair nor the names of its contributors may be
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

#ifndef WIREHAIR_H
#define WIREHAIR_H

/** \mainpage
    Wirehair : Fountain Codes in C++

    Wirehair produces a stream of error correction blocks from a data source
    using an erasure code.  When enough of these blocks are received,
    the original data can be recovered.  It supports up to 64000 blocks.

    This codec provides full-message reconstruction using the
    wirehair_recover() function, and can also
    reconstruct single missing blocks using the
    wirehair_recover_block() function.
    It also accelerates retransmission in a store-and-forward system
    with the wirehair_decoder_becomes_encoder() function.
*/

#define WIREHAIR_VERSION 2

// CMake propagates WIREHAIR_DLL or WIREHAIR_STATIC to installed consumers.
#if defined(WIREHAIR_STATIC)
    #define WIREHAIR_EXPORT extern
#elif defined(WIREHAIR_BUILDING)
# if defined(WIREHAIR_DLL) && defined(_WIN32)
    #define WIREHAIR_EXPORT __declspec(dllexport)
# else
    #define WIREHAIR_EXPORT
# endif
#else
# if defined(WIREHAIR_DLL) && defined(_WIN32)
    #define WIREHAIR_EXPORT __declspec(dllimport)
# else
    #define WIREHAIR_EXPORT extern
# endif
#endif

#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif


//------------------------------------------------------------------------------
// Shared Constants/Datatypes

/// These are the result codes that can be returned from the API functions
typedef enum WirehairResult_t
{
    /// Success code
    Wirehair_Success             = 0,

    /// More data is needed to decode.  This is normal and does not indicate a failure
    Wirehair_NeedMore            = 1,

    // Other values are failure codes:

    /// A function parameter was invalid
    Wirehair_InvalidInput        = 2,

    /// Encoder needs a better dense seed
    Wirehair_BadDenseSeed        = 3,

    /// Encoder needs a better peel seed
    Wirehair_BadPeelSeed         = 4,

    /// N = ceil(messageBytes / blockBytes) is too small.
    /// Try reducing block_size or use a larger message
    Wirehair_BadInput_SmallN     = 5,

    /// N = ceil(messageBytes / blockBytes) is too large.
    /// Try increasing block_size or use a smaller message
    Wirehair_BadInput_LargeN     = 6,

    /// Not enough extra rows to solve it, must give up
    Wirehair_ExtraInsufficient   = 7,

    /// An error occurred during the request
    Wirehair_Error               = 8,

    /// Out of memory
    Wirehair_OOM                 = 9,

    /// Platform is not supported yet
    Wirehair_UnsupportedPlatform = 10,

    WirehairResult_Count, /* for asserts */
    WirehairResult_Padding = 0x7fffffff /* int32_t padding */
} WirehairResult;

/// Get WirehairResult string function
WIREHAIR_EXPORT const char *wirehair_result_string(
    WirehairResult result ///< Result code to convert to string
);


//------------------------------------------------------------------------------
// Wirehair API

/**
    wirehair_init()

    Verify binary compatibility with the API on startup.

    Initialization is thread-safe.  Concurrent callers observe the same cached
    success or permanent platform/self-test failure.  A version mismatch is
    rejected for that caller and does not affect later calls using the correct
    version.

    Example:
        if (wirehair_init()) {
            exit(1);
        }

    Returns Wirehair_Success on success.
    Returns other codes on error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_init_(int expected_version);
#define wirehair_init() wirehair_init_(WIREHAIR_VERSION)

/// WirehairCodec: From wirehair_encoder_create() or wirehair_decoder_create()
typedef struct WirehairCodec_t { char impl; }* WirehairCodec;

//------------------------------------------------------------------------------
// Legacy wire-profile contract

/** Version of WirehairWireProfile, independent from the library ABI version. */
#define WIREHAIR_WIRE_PROFILE_VERSION 1u

/**
    Legacy-v2 equation profile immediately before exact-N seed fixups.

    This identifier is the first 64 bits of SHA-256 over the canonical profile
    name recorded in LEGACY_WIRE_PROFILES.md.  It is a compatibility identifier,
    not a security primitive.
*/
#define WIREHAIR_LEGACY_PROFILE_PRE_FIXUP UINT64_C(0xe1b9f77f1c90f680)

/** Legacy-v2 profile containing the frozen 2026-07 seed-fixup set. */
#define WIREHAIR_LEGACY_PROFILE_FIXUPS_2026_07 UINT64_C(0x4d241359db07bb07)

/**
    Profile used by the unframed legacy API.  This alias is frozen for ABI-v2;
    future equation profiles must receive new names and explicit selection.
*/
#define WIREHAIR_LEGACY_PROFILE_CURRENT WIREHAIR_LEGACY_PROFILE_FIXUPS_2026_07

/** Ask wirehair_encoder_create_profile_ex() to copy the input message. */
#define WIREHAIR_ENCODER_OWN_INPUT 1u

/**
    Identifies every equation-affecting choice used by a legacy codec.

    The structure itself is an in-process API object, not a byte serialization.
    Carry profile_id in authenticated or otherwise trusted application framing,
    using an application-defined byte order, then validate it with
    wirehair_wire_profile_init().  A profile identifier does not authenticate
    packet contents; callers still need an application digest or MAC.
*/
typedef struct WirehairWireProfile_t
{
    uint32_t struct_bytes;       ///< Must equal sizeof(WirehairWireProfile)
    uint32_t profile_version;    ///< WIREHAIR_WIRE_PROFILE_VERSION
    uint64_t profile_id;         ///< One of WIREHAIR_LEGACY_PROFILE_*
} WirehairWireProfile;

/**
    Initialize and validate a supported legacy wire profile descriptor.

    profileOut is cleared on failure.  This function is stateless and may be
    called before wirehair_init().  Unknown identifiers return
    Wirehair_InvalidInput.  Builds with private equation experiment knobs
    return Wirehair_UnsupportedPlatform for otherwise supported identifiers.
*/
WIREHAIR_EXPORT WirehairResult wirehair_wire_profile_init(
    uint64_t profileId,
    WirehairWireProfile* profileOut);

/*
    A codec has one checked lifecycle: encoder, active decoder, completed
    decoder, or converted encoder.  Encoder operations are accepted only by an
    encoder/converted encoder.  Decode is accepted by active/completed decoders;
    recovery before completion returns Wirehair_NeedMore.  Conversion is valid
    exactly once after decode completion.  Other mode mismatches return
    Wirehair_InvalidInput without writing caller data.
*/


/**
    wirehair_encoder_create()

    Encode the given message into blocks of size blockBytes.

    This zero-copy entry point borrows message.  The buffer must remain alive
    and immutable until the encoder is freed or successfully reused.  Use
    wirehair_encoder_create_owned() when that lifetime cannot be guaranteed.

    The number of blocks in the message:
        N = (bytes + (blockBytes-1)) / blockBytes

    If N is too high or too low this function will fail.  In particular if N = 1, then
    using this type of error correction does not make sense: Sending the same message
    over and over is just as good.  And if N > 64000 then some internal variables will
    start to overflow, so too many blocks is unsupported.  The most efficient values
    for N are around 1000.

    Preconditions:
        N >= 2
        N <= 64000
        blockBytes <= 2^31 - 1

    Pass 0 for reuseOpt if you do not want to reuse a WirehairCodec object.
    If this function fails, any codec passed via reuseOpt has been freed:
    do not reuse or wirehair_free() it after a failed call.

    Returns a non-zero object pointer on success.
    Returns nullptr(0) on failure.
*/
WIREHAIR_EXPORT WirehairCodec wirehair_encoder_create(
    WirehairCodec reuseOpt, ///< [Optional] Pointer to prior codec object
    const void*    message, ///< Pointer to message
    uint64_t  messageBytes, ///< Bytes in the message
    uint32_t    blockBytes  ///< Bytes in an output block
);

/**
    Result-preserving form of wirehair_encoder_create().

    The input message is borrowed for the lifetime of the encoder and must stay
    allocated and immutable until wirehair_free() or successful codec reuse.
    On success codecOut receives the encoder.  On failure codecOut is set to
    null and reuseOpt is consumed, except that a null codecOut is rejected
    without consuming reuseOpt.  Before initialization it returns Wirehair_Error;
    after a cached initialization failure it returns that platform error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_encoder_create_ex(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    WirehairCodec* codecOut
);

/**
    Owning form of wirehair_encoder_create().

    The message is copied before this function returns.  The caller may modify
    or release its buffer after success.  This costs messageBytes additional
    storage and O(messageBytes) copy time during creation.
*/
WIREHAIR_EXPORT WirehairCodec wirehair_encoder_create_owned(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes
);

/// Result-preserving form of wirehair_encoder_create_owned().
WIREHAIR_EXPORT WirehairResult wirehair_encoder_create_owned_ex(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    WirehairCodec* codecOut
);

/**
    Create an encoder for an explicit, trusted wire profile.

    flags may be zero (borrow message) or WIREHAIR_ENCODER_OWN_INPUT (copy it).
    Unknown/malformed profiles and unknown flags return Wirehair_InvalidInput
    before any packet can be emitted.  Private equation experiment builds
    return Wirehair_UnsupportedPlatform for a valid named profile.  As with the
    other result-preserving create APIs, reuseOpt is consumed on failure unless
    codecOut is null.
*/
WIREHAIR_EXPORT WirehairResult wirehair_encoder_create_profile_ex(
    WirehairCodec reuseOpt,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    const WirehairWireProfile* profile,
    uint32_t flags,
    WirehairCodec* codecOut
);

/**
    wirehair_encode()

    Write an error correction block.

    The first `blockId` < N blocks are the same as the input data.
    This "systematic" property can be used to run the encoder
    in parallel with sending original data.
    The `blockId` >= N blocks are generated on demand.

    Preconditions:
       Block is at least `blockBytes` in size

    Returns Wirehair_Success on success.
    Returns other codes on error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_encode(
    WirehairCodec    codec, ///< Pointer to codec from wirehair_encoder_init()
    unsigned       blockId, ///< Identifier of block to generate
    void*     blockDataOut, ///< Pointer to output block data
    uint32_t      outBytes, ///< Bytes in the output buffer
    uint32_t* dataBytesOut  ///< Number of bytes written <= blockBytes
);

/**
    wirehair_decoder_create()

    Initialize a decoder for a message of size bytes with blockBytes bytes
    per received block.  These parameters should match the ones passed to
    the corresponding wirehair_encoder_create() call.

    Pass 0 for reuseOpt if you do not want to reuse a WirehairCodec object.
    If this function fails, any codec passed via reuseOpt has been freed:
    do not reuse or wirehair_free() it after a failed call.

    Returns a non-zero object pointer on success.
    Returns nullptr(0) on failure.
*/
WIREHAIR_EXPORT WirehairCodec wirehair_decoder_create(
    WirehairCodec reuseOpt, ///< Codec object to reuse
    uint64_t  messageBytes, ///< Bytes in the message to decode
    uint32_t    blockBytes  ///< Bytes in each encoded block
);

/**
    Result-preserving form of wirehair_decoder_create().

    On success codecOut receives the decoder.  On failure codecOut is set to
    null and reuseOpt is consumed, except that a null codecOut is rejected
    without consuming reuseOpt.  Before initialization it returns Wirehair_Error;
    after a cached initialization failure it returns that platform error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_decoder_create_ex(
    WirehairCodec reuseOpt,
    uint64_t messageBytes,
    uint32_t blockBytes,
    WirehairCodec* codecOut
);

/**
    Create a decoder for an explicit, trusted wire profile.

    Descriptor and experiment-build failure behavior matches
    wirehair_encoder_create_profile_ex().
*/
WIREHAIR_EXPORT WirehairResult wirehair_decoder_create_profile_ex(
    WirehairCodec reuseOpt,
    uint64_t messageBytes,
    uint32_t blockBytes,
    const WirehairWireProfile* profile,
    WirehairCodec* codecOut
);

/**
    wirehair_decode()

    Provide the decoder with a block from the wirehair_encode() function.

    Preconditions:
    The application must not provide duplicate data for the same packet.
    In other words, blockId values cannot be used twice.

    Returns Wirehair_Success if data recovery is complete.
    + Use wirehair_recover() or wirehair_recover_block()
      to reconstruct the recovered data.
    Returns Wirehair_NeedMore if more data is needed to decode.
    Returns other codes on error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_decode(
    WirehairCodec   codec, ///< Codec object
    unsigned      blockId, ///< ID number of received block
    const void* blockData, ///< Pointer to block data
    uint32_t    dataBytes  ///< Number of bytes in the data block
);

/**
    wirehair_recover()

    Reconstruct the message after reading is complete.

    Preconditions:
    Message contains enough space to store the entire decoded message (bytes)

    Returns Wirehair_Success if the message was recovered.
    Returns other codes on error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_recover(
    WirehairCodec    codec, ///< Codec object
    void*       messageOut, ///< Buffer where reconstructed message will be written
    uint64_t  messageBytes  ///< Bytes in the message
);

/**
    wirehair_recover_block()

    Reconstruct a single block of the message after reading is complete.
    If the requested original block was received directly, it is copied from
    decoder staging. Otherwise it is regenerated from recovery state, which can
    be much slower than wirehair_recover().

    Preconditions:
    Block ptr buffer contains enough space to hold the block (blockBytes)

    This legacy entry point cannot verify the output capacity.  New code should
    use wirehair_recover_block_ex().

    May return non-zero to indicate a failure.

    Returns Wirehair_Success if the block was recovered.
    Returns other codes on error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_recover_block(
    WirehairCodec codec, ///< Codec object
    unsigned    blockId, ///< ID of the block to reconstruct between 0..N-1
    void*     blockData, ///< Pointer to block data
    uint32_t*  bytesOut  ///< Set to the number of data bytes in the block
);

/**
    Capacity-checked form of wirehair_recover_block().

    outputCapacity is the number of writable bytes at blockData.  On failure,
    blockData is unchanged and bytesOut is set to zero when it is non-null.
*/
WIREHAIR_EXPORT WirehairResult wirehair_recover_block_ex(
    WirehairCodec codec,
    unsigned blockId,
    void* blockData,
    uint32_t outputCapacity,
    uint32_t* bytesOut
);

/**
    wirehair_decoder_becomes_encoder()

    Convert a decoder WirehairCodec into an encoder WirehairCodec
    after recovery is complete.  This enables the application to
    receive a message and then retransmit it without reinitializing
    the encoder from scratch.

    Note that after converting the codec to an encoder, calling the
    wirehair_encode() function to retrieve original data will be
    much slower than accessing the original data directly on the
    application side.

    Preconditions:
    wirehair_decode() returned Wirehair_Success

    Returns Wirehair_Success if the operation was successful.
    Returns other codes on error.
*/
WIREHAIR_EXPORT WirehairResult wirehair_decoder_becomes_encoder(
    WirehairCodec codec ///< Codec to change
);

/**
    wirehair_free()

    Free memory associated with a WirehairCodec object.
*/
WIREHAIR_EXPORT void wirehair_free(
    WirehairCodec codec ///< Codec object to free
);


#ifdef __cplusplus
}
#endif

#endif // WIREHAIR_H
