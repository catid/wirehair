#pragma once

#include "../codec/WirehairV2Solve.h"

#include <cstddef>
#include <cstdint>

namespace wh2_repair_alignment_r0 {

// Benchmark-only borrowed observation.  The original public encoder and its
// immutable source must remain alive and unchanged for every use of this view.
// Capture accepts live handles only; it cannot validate dangling/forged handles.
struct Snapshot
{
    WirehairV2Codec handle = nullptr;
    const uint8_t* source = nullptr;
    const uint8_t* intermediate = nullptr;
    uint64_t message_bytes = 0;
    uint32_t block_bytes = 0;
    uint32_t source_count = 0;
    uint32_t precode_count = 0;
    size_t intermediate_bytes = 0;
    WirehairV2EncoderSourcePolicy source_policy = WirehairV2EncoderSource_Invalid;
    WirehairV2Profile profile = {};
    uint8_t serialized_profile[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    wirehair_v2::SeedProfile expanded_profile = {};
    // Only Params are retained in a solved production packet encoder.  Own an
    // identical immutable descriptor so every shadow treatment shares metadata.
    wirehair_v2::PrecodeSystem system = {};
    wirehair_v2::PacketRowConfig config = {};
    wirehair_v2::PacketRowRuntime runtime = {};
};

struct Evaluation
{
    WirehairV2Result status = WirehairV2_InvalidInput;
    uint32_t bytes = 0;
    uint64_t operations = 0;
};

// Fixed diagnostic domain: borrowed K6, B1280, M7680, certified 2026_07.
// Uses the original facade's full named-profile reconstruction/validation.
WirehairV2Result Capture(WirehairV2Codec live_encoder, Snapshot& output);

// Pure metadata/address checks, with no dereferences, codec/GF work or clocks.
// Call before timing; output must be disjoint from the selected intermediate,
// original intermediate and borrowed source.  Actual capacities are supplied
// by the owner, not discoverable from pointers.
bool ValidateShadowBuffers(
    const Snapshot& snapshot,
    const uint8_t* intermediate,
    size_t intermediate_capacity,
    uint8_t* output,
    size_t output_capacity);

// Evaluate IDs6..11 with the ORIGINAL production WithRuntime entry.  Snapshot
// must be a successful unmodified Capture and buffers must already have passed
// ValidateShadowBuffers.  Same-byte shadow verification is the driver's duty.
// A nonnull disjoint local operation scratch mirrors PrecodeEncoder::EncodeResult.
Evaluation EvaluateShadow(
    const Snapshot& snapshot,
    const uint8_t* intermediate,
    uint32_t block_id,
    uint8_t* output);

// Synthetic dimensions/profile-field consistency and bounds only: no capture,
// profile expansion, GF initialization, encoder construction or packet calls.
bool NeutralSelfTest();

} // namespace wh2_repair_alignment_r0
