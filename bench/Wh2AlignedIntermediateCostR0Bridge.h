#pragma once

// Public C types only.  No baseline or candidate V2 C++ class crosses this ABI.
#include <wirehair/wirehair.h>

#include <type_traits>

struct Wh2AlignedIntermediateCostSnapshot
{
    WirehairV2Codec handle;
    const uint8_t* source;
    const uint8_t* intermediate;
    uint64_t message_bytes;
    uint64_t intermediate_bytes;
    uint64_t intermediate_capacity;
    uint32_t block_bytes;
    uint32_t source_count;
    uint32_t precode_count;
    uint32_t source_policy;
    uint8_t serialized_profile[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
};

static_assert(std::is_standard_layout<Wh2AlignedIntermediateCostSnapshot>::value &&
    std::is_trivial<Wh2AlignedIntermediateCostSnapshot>::value,
    "Inspection must exchange only identical POD state");

// The source and encoder must stay alive and immutable while these observations
// are used.  Only live handles are accepted; forged/dangling handles cannot be
// diagnosed safely.  Inspect leaves output unchanged on failure.  Callers own
// valid, disjoint output objects; none of these helpers belongs inside a timer.
extern "C" WirehairV2Result wh2_aligned_cost_p_inspect(
    WirehairV2Codec encoder, Wh2AlignedIntermediateCostSnapshot* output);
extern "C" WirehairV2Result wh2_aligned_cost_a_inspect(
    WirehairV2Codec encoder, Wh2AlignedIntermediateCostSnapshot* output);

// Canonical original row indices and an independent scalar XOR packet oracle.
// Systematic last-tail output is truncated exactly like the public C API.
// BufferTooSmall publishes only the required count/bytes; other failures leave
// both outputs unchanged.  Null buffer/capacity zero is a supported size query.
extern "C" WirehairV2Result wh2_aligned_cost_p_packet_columns(
    WirehairV2Codec encoder, uint32_t id, uint32_t* columns,
    uint32_t capacity, uint32_t* count_out);
extern "C" WirehairV2Result wh2_aligned_cost_a_packet_columns(
    WirehairV2Codec encoder, uint32_t id, uint32_t* columns,
    uint32_t capacity, uint32_t* count_out);
extern "C" WirehairV2Result wh2_aligned_cost_p_scalar_packet(
    WirehairV2Codec encoder, uint32_t id, void* output,
    uint32_t capacity, uint32_t* bytes_out);
extern "C" WirehairV2Result wh2_aligned_cost_a_scalar_packet(
    WirehairV2Codec encoder, uint32_t id, void* output,
    uint32_t capacity, uint32_t* bytes_out);

// Declaration types come directly from the unchanged C header.  Rename every
// facade entry in the candidate: its wrappers also call these entries internally.
#define WH2_COST_DECLARE_CANDIDATE(name) \
    extern "C" decltype(wirehair_v2_##name) wh2_aligned_r0_##name;
WH2_COST_DECLARE_CANDIDATE(decode)
WH2_COST_DECLARE_CANDIDATE(decoder_create)
WH2_COST_DECLARE_CANDIDATE(encode)
WH2_COST_DECLARE_CANDIDATE(encoder_create)
WH2_COST_DECLARE_CANDIDATE(encoder_create_profile)
WH2_COST_DECLARE_CANDIDATE(encoder_create_profile_id)
WH2_COST_DECLARE_CANDIDATE(encoder_create_profile_id_with_options)
WH2_COST_DECLARE_CANDIDATE(encoder_create_profile_with_options)
WH2_COST_DECLARE_CANDIDATE(encoder_create_with_options)
WH2_COST_DECLARE_CANDIDATE(encoder_detach_input)
WH2_COST_DECLARE_CANDIDATE(free)
WH2_COST_DECLARE_CANDIDATE(profile_deserialize)
WH2_COST_DECLARE_CANDIDATE(profile_serialize)
WH2_COST_DECLARE_CANDIDATE(profile_validate)
WH2_COST_DECLARE_CANDIDATE(recover)
WH2_COST_DECLARE_CANDIDATE(result_string)
#undef WH2_COST_DECLARE_CANDIDATE
