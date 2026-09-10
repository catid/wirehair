// Engineering-only retained-corpus replay through the sealed benchmark boundary
// AND the actual WHV2 library. This is not a timing run or a new recovery sample.
// Compile this TU plus Wh2FrozenTrace.cpp with WH2_SMALL_CODEC_K=3,5,8 and the
// corresponding sealed generated-data include directory. Link that dimension's
// qualified prototype archive followed by the current ordinary library; match
// their private GF backend/sanitizer flags. The default installed build does not
// include this harness or depend on retained experiment files.
//
// WHKx and WHV2 are intentionally different descriptors. Translation below is
// explicit test code, never a production fallback. Only valid corpus operations
// (Success/NeedMore) are compared: ownership/detach/conflict contracts differ and
// are independently tested by V2SmallCodecTest, not silently equated here.
#include "Wh2SmallSerialized.h"
#include "wirehair/wirehair.h"
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

static_assert(WH2_SMALL_CODEC_K == 3 || WH2_SMALL_CODEC_K == 5 || WH2_SMALL_CODEC_K == 8,
              "Installed small dimensions only");
namespace wh2_v2_production_parity {
constexpr uint64_t ProfileId = WH2_SMALL_CODEC_K == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
    WH2_SMALL_CODEC_K == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09;
using Profile = std::array<unsigned char, 32>;
void Require(bool value, const char* why)
{
    if (!value) { std::fprintf(stderr, "FAIL WHV2/prototype parity: %s\n", why); std::abort(); }
}
Wh2SmallStatus Status(WirehairV2Result result)
{
    if (result == WirehairV2_Success) return Wh2Small_Success;
    Require(result == WirehairV2_NeedMore, "unexpected installed result on valid corpus operation");
    return Wh2Small_NeedMore;
}
uint64_t Load(const unsigned char* bytes, unsigned count)
{
    uint64_t value = 0;
    for (unsigned i = 0; i < count; ++i) value |= uint64_t(bytes[i]) << (8 * i);
    return value;
}
Profile Translate(const void* input, size_t bytes, uint64_t& message, uint32_t& block)
{
    Require(input && bytes == 32 && wh2_small_profile_validate(input, bytes) == Wh2Small_Success,
            "sealed benchmark descriptor");
    const auto* p = static_cast<const unsigned char*>(input);
    Require(Load(p + 8, 8) == WH2_SMALL_PROFILE_ID, "exact prototype equation identity");
    message = Load(p + 16, 8); block = static_cast<uint32_t>(Load(p + 24, 4));
    Require(block && block <= 4096 && message > uint64_t(WH2_SMALL_CODEC_K - 1) * block &&
            message <= uint64_t(WH2_SMALL_CODEC_K) * block && message <= SIZE_MAX - 2,
            "bounded retained shape");
    WirehairV2Profile host = {};
    host.struct_bytes = sizeof(host); host.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    host.profile_id = ProfileId; host.message_bytes = message; host.block_bytes = block;
    Profile out = {}; uint32_t written = 0;
    Require(wirehair_v2_profile_serialize(&host, out.data(), 32, &written) == WirehairV2_Success &&
            written == 32, "explicit WHV2 descriptor translation");
    return out;
}
struct Pair {
    Wh2SmallCodec prototype = nullptr;
    WirehairV2Codec integrated = nullptr;
    size_t message = 0;
    uint32_t block = 0;
    std::unique_ptr<unsigned char[]> scratch;
    ~Pair() { wh2_small_free(prototype); wirehair_v2_free(integrated); }
    void Prepare(size_t capacity)
    {
        Require(capacity <= message && capacity <= UINT32_MAX, "bounded replay output capacity");
        std::memset(scratch.get(), 0xa5, message + 2);
    }
    void Check(const Wh2SmallResult& old, WirehairV2Result now, uint64_t required, const void* output)
    {
        Require(old.status == Status(now) && old.bytes_required == required && required <= message &&
                old.bytes_written == (now == WirehairV2_Success ? required : 0), "status/length fields");
        const size_t written = static_cast<size_t>(old.bytes_written);
        Require(!std::memcmp(output, scratch.get() + 1, written), "output bytes");
        Require(scratch[0] == 0xa5 &&
                std::all_of(scratch.get() + 1 + written, scratch.get() + message + 2,
                            [](unsigned char b) { return b == 0xa5; }), "installed no-write and guards");
    }
};
Wh2SmallCreateResult Encoder(const void* source, uint64_t message, uint32_t block,
                            uint32_t policy, void* profile, size_t capacity)
{
    Require(capacity == 32 && (policy == Wh2Small_Independent || policy == Wh2Small_BorrowedImmutable),
            "bounded encoder arguments");
    std::unique_ptr<Pair> pair(new Pair);
    const auto old = wh2_small_encoder_create(source, message, block, policy, profile, capacity);
    pair->prototype = old.codec;
    Require(old.status == Wh2Small_Success && old.codec, "prototype encoder");
    uint64_t parsed_message = 0; uint32_t parsed_block = 0;
    const Profile expected = Translate(profile, capacity, parsed_message, parsed_block);
    Require(parsed_message == message && parsed_block == block, "prototype descriptor shape");
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = policy == Wh2Small_Independent ? WirehairV2EncoderSource_Independent :
        WirehairV2EncoderSource_BorrowedImmutable;
    Profile actual = {}; uint32_t written = 0;
    const auto now = wirehair_v2_encoder_create_profile_id_with_options(ProfileId, source, message,
        block, &options, actual.data(), 32, &written, &pair->integrated);
    Require(now == WirehairV2_Success && pair->integrated && written == 32 && actual == expected,
            "installed encoder and exact translated descriptor");
    pair->message = static_cast<size_t>(message); pair->block = block;
    pair->scratch.reset(new unsigned char[pair->message + 2]);
    return {Wh2Small_Success, pair.release()};
}
Wh2SmallCreateResult Decoder(const void* profile, size_t bytes)
{
    uint64_t message = 0; uint32_t block = 0;
    const Profile translated = Translate(profile, bytes, message, block);
    std::unique_ptr<Pair> pair(new Pair);
    const auto old = wh2_small_decoder_create(profile, bytes);
    pair->prototype = old.codec;
    const auto now = wirehair_v2_decoder_create(translated.data(), 32, &pair->integrated);
    Require(old.status == Wh2Small_Success && old.codec && now == WirehairV2_Success && pair->integrated,
            "standalone receivers without surviving encoder");
    pair->message = static_cast<size_t>(message); pair->block = block;
    pair->scratch.reset(new unsigned char[pair->message + 2]);
    return {Wh2Small_Success, pair.release()};
}
Wh2SmallResult Encode(Wh2SmallCodec handle, uint32_t id, void* output, size_t capacity)
{
    auto* pair = static_cast<Pair*>(handle); pair->Prepare(capacity);
    const auto old = wh2_small_encode(pair->prototype, id, output, capacity);
    uint32_t required = 0;
    const auto now = wirehair_v2_encode(pair->integrated, id, pair->scratch.get() + 1,
        static_cast<uint32_t>(capacity), &required);
    const uint64_t expected = id == WH2_SMALL_CODEC_K - 1 ?
        pair->message - uint64_t(WH2_SMALL_CODEC_K - 1) * pair->block : pair->block;
    Require(required == expected, "independent packet length");
    pair->Check(old, now, required, output); return old;
}
Wh2SmallStatus Decode(Wh2SmallCodec handle, uint32_t id, const void* input, size_t bytes)
{
    auto* pair = static_cast<Pair*>(handle);
    Require(bytes <= pair->block, "bounded packet length");
    const auto old = wh2_small_decode(pair->prototype, id, input, bytes);
    const auto now = wirehair_v2_decode(pair->integrated, id, input, static_cast<uint32_t>(bytes));
    Require(old == Status(now), "every prefix rank/status"); return old;
}
Wh2SmallResult Recover(Wh2SmallCodec handle, void* output, size_t capacity)
{
    auto* pair = static_cast<Pair*>(handle); pair->Prepare(capacity);
    const auto old = wh2_small_recover(pair->prototype, output, capacity);
    uint64_t required = 0;
    const auto now = wirehair_v2_recover(pair->integrated, pair->scratch.get() + 1, capacity, &required);
    Require(required == pair->message, "independent recovered length including NeedMore");
    pair->Check(old, now, required, output); return old;
}
void Free(Wh2SmallCodec handle) { delete static_cast<Pair*>(handle); }
} // namespace wh2_v2_production_parity

#define WH2_SMALL_TEST_SERIALIZED 1
#define wh2_small_encoder_create wh2_v2_production_parity::Encoder
#define wh2_small_decoder_create wh2_v2_production_parity::Decoder
#define wh2_small_encode wh2_v2_production_parity::Encode
#define wh2_small_decode wh2_v2_production_parity::Decode
#define wh2_small_recover wh2_v2_production_parity::Recover
#define wh2_small_free wh2_v2_production_parity::Free
#include "Wh2SmallNativeTest.cpp"
