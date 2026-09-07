#include "Wh2ThueMorseMatvecCostR0Bridge.h"

#if !defined(WH2_MATVEC_COST_ARM) || WH2_MATVEC_COST_ARM < 0 || WH2_MATVEC_COST_ARM > 3
#error "Compile this bridge with WH2_MATVEC_COST_ARM=0..3"
#endif
#ifdef WH_COUNT
#error "Cost bridge must retain ordinary GF work"
#endif

#if WH2_MATVEC_COST_ARM == 0
#include "Wh2ThueMorseNativeCodec.h"
namespace N = wh2_thue_native;
#define BRIDGE_NAMESPACE wh2_matvec_cost_baseline_bridge
#define BRIDGE_GETTER wh2_matvec_cost_baseline_api
#elif WH2_MATVEC_COST_ARM == 1
#define wh2_thue_native wh2_thue_matvec_codec_r0
#include "Wh2ThueMorseNativeCodec.h"
#undef wh2_thue_native
namespace N = wh2_thue_matvec_codec_r0;
#define BRIDGE_NAMESPACE wh2_matvec_cost_candidate_bridge
#define BRIDGE_GETTER wh2_matvec_cost_candidate_api
#elif WH2_MATVEC_COST_ARM == 2
#include "wirehair/wirehair.h"
#define BRIDGE_NAMESPACE wh2_matvec_cost_wh1_bridge
#define BRIDGE_GETTER wh2_matvec_cost_wh1_api
#else
#include "wirehair/wirehair.h"
#define BRIDGE_NAMESPACE wh2_matvec_cost_public_bridge
#define BRIDGE_GETTER wh2_matvec_cost_public_api
#endif

#include <algorithm>
#include <cstring>
#include <limits>
#include <memory>

namespace BRIDGE_NAMESPACE {
namespace B = wh2_matvec_cost_r0;
namespace {

bool Span(const void* pointer, std::size_t bytes)
{
    return pointer && bytes <= std::numeric_limits<std::uintptr_t>::max() -
        reinterpret_cast<std::uintptr_t>(pointer);
}

bool Overlap(const void* a, std::size_t an, const void* b, std::size_t bn)
{
    const std::uintptr_t first = reinterpret_cast<std::uintptr_t>(a);
    const std::uintptr_t second = reinterpret_cast<std::uintptr_t>(b);
    return an && bn && first < second + bn && second < first + an;
}

bool StateSpan(const B::State* state)
{
    return Span(state, sizeof(*state)) &&
        reinterpret_cast<std::uintptr_t>(state) % alignof(B::State) == 0;
}

bool Dimensions(std::uint64_t message, std::uint32_t block)
{
    return block && block <= static_cast<std::uint32_t>(std::numeric_limits<int>::max()) &&
        message <= std::numeric_limits<std::size_t>::max() &&
        message > 5u * static_cast<std::uint64_t>(block) &&
        message <= 6u * static_cast<std::uint64_t>(block);
}

bool Live(const B::State* state, unsigned kind)
{
    return StateSpan(state) && state->handle && state->kind == kind &&
        state->arm == WH2_MATVEC_COST_ARM &&
        Dimensions(state->message_bytes, state->block_bytes) &&
        state->lookup_bytes == 39936u && Span(state->lookup, state->lookup_bytes) &&
        (kind != 1u || Span(state->source, static_cast<std::size_t>(state->message_bytes)));
}

B::Result Answer(int code, std::uint64_t required = 0, std::uint64_t written = 0,
                 unsigned length_kind = B::NoLength)
{
    return B::Result{code == 0 ? B::Success : code == 1 ? B::NeedMore : B::Failure,
                     code, required, written, length_kind};
}

std::uint64_t PacketBytes(const B::State& state, std::uint32_t id)
{
    return id == 5u ? state.message_bytes - 5u * static_cast<std::uint64_t>(state.block_bytes) :
        state.block_bytes;
}

bool Writable(const B::State& state, void* output, std::size_t capacity, std::size_t bytes)
{
    return Span(output, capacity) && Span(output, bytes) &&
        !Overlap(output, bytes, &state, sizeof(state)) &&
        !Overlap(output, bytes, state.lookup, state.lookup_bytes) &&
        (!state.source || !Overlap(output, bytes, state.source,
                                   static_cast<std::size_t>(state.message_bytes)));
}

bool ValidProfile(const B::State* state)
{
    if (!StateSpan(state) || (state->kind != 1u && state->kind != 2u) ||
        !Live(state, state->kind)) return false;
#if WH2_MATVEC_COST_ARM == 3
    static_assert(WIREHAIR_V2_PROFILE_SERIALIZED_BYTES == 32u, "Exact public descriptor");
    WirehairV2Profile profile = {};
    return state->profile_bytes == sizeof(state->profile) &&
        wirehair_v2_profile_validate(state->profile, state->profile_bytes) == WirehairV2_Success &&
        wirehair_v2_profile_deserialize(state->profile, state->profile_bytes, &profile) == WirehairV2_Success &&
        profile.profile_id == WIREHAIR_V2_PROFILE_CURRENT &&
        profile.message_bytes == state->message_bytes && profile.block_bytes == state->block_bytes;
#else
    return state->profile_bytes == 0u &&
        std::all_of(state->profile, state->profile + sizeof(state->profile),
                    [](std::uint8_t byte) { return byte == 0; });
#endif
}

B::Result Create(const std::uint8_t* lookup, std::size_t lookup_bytes,
                 const void* source, std::uint64_t message, std::uint32_t block,
                 B::State* output)
{
    if (!StateSpan(output) || output->handle || output->kind ||
        !Dimensions(message, block) || lookup_bytes != 39936u ||
        !Span(lookup, lookup_bytes) || !Span(source, static_cast<std::size_t>(message)) ||
        Overlap(output, sizeof(*output), lookup, lookup_bytes) ||
        Overlap(output, sizeof(*output), source, static_cast<std::size_t>(message)))
        return Answer(-1);
    B::State staged = {};
    staged.source = source;
    staged.lookup = lookup;
    staged.lookup_bytes = lookup_bytes;
    staged.message_bytes = message;
    staged.block_bytes = block;
    staged.arm = WH2_MATVEC_COST_ARM;
    staged.kind = 1u;
    int code;
#if WH2_MATVEC_COST_ARM < 2
    std::unique_ptr<N::Encoder> encoder;
    code = static_cast<int>(N::Encoder::Create(N::Lookup{lookup, lookup_bytes}, source,
                                             message, block, encoder));
    if (code == 0) staged.handle = encoder.release();
#elif WH2_MATVEC_COST_ARM == 2
    WirehairCodec encoder = nullptr;
    code = wirehair_encoder_create_ex(nullptr, source, message, block, &encoder);
    staged.handle = encoder;
#else
    WirehairV2EncoderOptions options = {};
    options.struct_bytes = sizeof(options);
    options.options_version = WIREHAIR_V2_ENCODER_OPTIONS_VERSION;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    WirehairV2Codec encoder = nullptr;
    code = wirehair_v2_encoder_create_with_options(source, message, block, &options,
        staged.profile, sizeof(staged.profile), &staged.profile_bytes, &encoder);
    staged.handle = encoder;
#endif
    if (code == 0) {
        if (!staged.handle) return Answer(-1);
        *output = staged;
    }
    return Answer(code);
}

B::Result Encode(const B::State* state, std::uint32_t id, void* output, std::size_t capacity)
{
    if (!Live(state, 1u) || capacity > UINT32_MAX) return Answer(-1);
    const std::uint64_t required = PacketBytes(*state, id);
    if (!Writable(*state, output, capacity, static_cast<std::size_t>(required))) return Answer(-1);
#if WH2_MATVEC_COST_ARM < 2
    const N::Result result = static_cast<N::Encoder*>(state->handle)->Encode(id, output, capacity);
    return Answer(static_cast<int>(result.status), result.bytes_required, result.bytes_written,
                  B::ReturnedLength);
#else
    std::uint32_t bytes = UINT32_MAX;
#if WH2_MATVEC_COST_ARM == 2
    const int code = wirehair_encode(static_cast<WirehairCodec>(state->handle), id, output,
                                    static_cast<std::uint32_t>(capacity), &bytes);
#else
    const int code = wirehair_v2_encode(static_cast<WirehairV2Codec>(state->handle), id, output,
                                       static_cast<std::uint32_t>(capacity), &bytes);
#endif
    return bytes == UINT32_MAX ? Answer(code) :
        Answer(code, bytes, code == 0 ? bytes : 0u, B::ReturnedLength);
#endif
}

void EncoderFree(B::State* state)
{
    if (!StateSpan(state) || !state->handle || state->arm != WH2_MATVEC_COST_ARM ||
        state->kind != 1u) return;
#if WH2_MATVEC_COST_ARM < 2
    delete static_cast<N::Encoder*>(state->handle);
#elif WH2_MATVEC_COST_ARM == 2
    wirehair_free(static_cast<WirehairCodec>(state->handle));
#else
    wirehair_v2_free(static_cast<WirehairV2Codec>(state->handle));
#endif
    *state = B::State{};
}

B::Result DecoderCreate(const B::State* encoder, B::State* output)
{
    if (!Live(encoder, 1u) || !StateSpan(output) || output->handle || output->kind ||
        Overlap(output, sizeof(*output), encoder, sizeof(*encoder)) ||
        Overlap(output, sizeof(*output), encoder->source, static_cast<std::size_t>(encoder->message_bytes)) ||
        Overlap(output, sizeof(*output), encoder->lookup, encoder->lookup_bytes)) return Answer(-1);
    B::State staged = *encoder;
    staged.handle = nullptr;
    staged.source = nullptr; // Decoders do not borrow the source message.
    staged.kind = 2u;
    int code;
#if WH2_MATVEC_COST_ARM < 2
    std::unique_ptr<N::Decoder> decoder;
    code = static_cast<int>(N::Decoder::Create(N::Lookup{staged.lookup, staged.lookup_bytes},
        staged.message_bytes, staged.block_bytes, decoder));
    if (code == 0) staged.handle = decoder.release();
#elif WH2_MATVEC_COST_ARM == 2
    WirehairCodec decoder = nullptr;
    code = wirehair_decoder_create_ex(nullptr, staged.message_bytes, staged.block_bytes, &decoder);
    staged.handle = decoder;
#else
    if (staged.profile_bytes != sizeof(staged.profile)) return Answer(-1);
    WirehairV2Codec decoder = nullptr;
    code = wirehair_v2_decoder_create(staged.profile, staged.profile_bytes, &decoder);
    staged.handle = decoder;
#endif
    if (code == 0) {
        if (!staged.handle) return Answer(-1);
        *output = staged;
    }
    return Answer(code);
}

B::Result Feed(B::State* state, std::uint32_t id, const void* input, std::size_t bytes)
{
    if (!Live(state, 2u) || bytes > UINT32_MAX || !Span(input, bytes) ||
        Overlap(input, bytes, state, sizeof(*state))) return Answer(-1);
#if WH2_MATVEC_COST_ARM < 2
    const N::Result result = static_cast<N::Decoder*>(state->handle)->Feed(id, input, bytes);
    return Answer(static_cast<int>(result.status), result.bytes_required, result.bytes_written,
                  B::ReturnedLength);
#elif WH2_MATVEC_COST_ARM == 2
    const int code = wirehair_decode(static_cast<WirehairCodec>(state->handle), id, input,
                                    static_cast<std::uint32_t>(bytes));
    return Answer(code, PacketBytes(*state, id), 0, B::InferredLength);
#else
    const int code = wirehair_v2_decode(static_cast<WirehairV2Codec>(state->handle), id, input,
                                       static_cast<std::uint32_t>(bytes));
    return Answer(code, PacketBytes(*state, id), 0, B::InferredLength);
#endif
}

B::Result Recover(B::State* state, void* output, std::size_t capacity)
{
    if (!Live(state, 2u) ||
        !Writable(*state, output, capacity, static_cast<std::size_t>(state->message_bytes)))
        return Answer(-1);
#if WH2_MATVEC_COST_ARM < 2
    const N::Result result = static_cast<N::Decoder*>(state->handle)->Recover(output, capacity);
    return Answer(static_cast<int>(result.status), result.bytes_required, result.bytes_written,
                  B::ReturnedLength);
#elif WH2_MATVEC_COST_ARM == 2
    // WH1's final parameter is exact message length, not output capacity.
    if (capacity < state->message_bytes) return Answer(-1, state->message_bytes, 0, B::InferredLength);
    const int code = wirehair_recover(static_cast<WirehairCodec>(state->handle), output,
                                     state->message_bytes);
    return Answer(code, state->message_bytes, code == 0 ? state->message_bytes : 0u, B::InferredLength);
#else
    std::uint64_t bytes = UINT64_MAX;
    const int code = wirehair_v2_recover(static_cast<WirehairV2Codec>(state->handle), output,
                                       capacity, &bytes);
    return bytes == UINT64_MAX ? Answer(code) :
        Answer(code, bytes, code == 0 ? bytes : 0u, B::ReturnedLength);
#endif
}

void DecoderFree(B::State* state)
{
    if (!StateSpan(state) || !state->handle || state->arm != WH2_MATVEC_COST_ARM ||
        state->kind != 2u) return;
#if WH2_MATVEC_COST_ARM < 2
    delete static_cast<N::Decoder*>(state->handle);
#elif WH2_MATVEC_COST_ARM == 2
    wirehair_free(static_cast<WirehairCodec>(state->handle));
#else
    wirehair_v2_free(static_cast<WirehairV2Codec>(state->handle));
#endif
    *state = B::State{};
}

#if WH2_MATVEC_COST_ARM < 2
B::Result Row(const std::uint8_t* lookup, std::size_t bytes,
              std::uint32_t id, std::uint8_t* output)
{
    const int code = static_cast<int>(N::Row(N::Lookup{lookup, bytes}, id, output));
    return Answer(code, 6u, code == 0 ? 6u : 0u, B::InferredLength);
}
#endif

const B::Api api = {Create, Encode, EncoderFree, DecoderCreate, Feed, Recover, DecoderFree,
    ValidProfile,
#if WH2_MATVEC_COST_ARM < 2
    Row
#else
    nullptr
#endif
};
} // namespace
} // namespace BRIDGE_NAMESPACE

extern "C" const wh2_matvec_cost_r0::Api* BRIDGE_GETTER()
{
    return &BRIDGE_NAMESPACE::api;
}

#undef BRIDGE_NAMESPACE
#undef BRIDGE_GETTER
