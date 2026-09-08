// Engineering-only replay of the retained native corpus through BOTH the
// qualified serialized prototype archive and the actual integrated library.
// Compile this TU plus Wh2FrozenTrace.cpp with the sealed generated data header
// on the include path, linking the old prototype archive then the new library.
// Neither the ordinary library nor its installed tests depend on this harness.
#include "Wh2SmallSerialized.h"
#include "wirehair/wirehair_small.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

namespace wh2_production_parity {
void Require(bool value, const char* why)
{
    if (!value) { std::fprintf(stderr, "FAIL production/prototype parity: %s\n", why); std::abort(); }
}
struct Pair {
    Wh2SmallCodec prototype = nullptr;
    WirehairSmallCodec integrated = nullptr;
    size_t message = 0;
    std::unique_ptr<unsigned char[]> scratch;
    ~Pair() { wh2_small_free(prototype); wirehair_small_free(integrated); }
    void Prepare(size_t capacity)
    {
        Require(capacity <= message, "bounded replay output capacity");
        std::memset(scratch.get(), 0xa5, message + 2);
    }
    void Check(const Wh2SmallResult& old, const WirehairSmallResult& now, const void* output)
    {
        Require(int(old.status) == int(now.status) && old.bytes_required == now.bytes_required &&
                old.bytes_written == now.bytes_written && now.bytes_written <= message, "result fields");
        const size_t written = static_cast<size_t>(now.bytes_written);
        Require(!std::memcmp(output, scratch.get() + 1, written), "output bytes");
        Require(scratch[0] == 0xa5 &&
                std::all_of(scratch.get() + 1 + written, scratch.get() + message + 2,
                            [](unsigned char b) { return b == 0xa5; }), "integrated no-write and guards");
    }
};
Wh2SmallCreateResult Encoder(const void* source, uint64_t message, uint32_t block,
                            uint32_t policy, void* profile, size_t capacity)
{
    Require(capacity == 32 && message <= SIZE_MAX - 2, "bounded encoder shape");
    std::unique_ptr<Pair> pair(new Pair);
    unsigned char other[32] = {};
    const auto old = wh2_small_encoder_create(source, message, block, policy, profile, capacity);
    pair->prototype = old.codec;
    const auto now = wirehair_small_encoder_create(source, message, block, policy, other, sizeof(other));
    pair->integrated = now.codec;
    Require(old.status == Wh2Small_Success && now.status == WirehairSmall_Success &&
            old.codec && now.codec && !std::memcmp(profile, other, 32), "encoder/descriptor");
    pair->message = static_cast<size_t>(message);
    pair->scratch.reset(new unsigned char[pair->message + 2]);
    return {Wh2Small_Success, pair.release()};
}
Wh2SmallCreateResult Decoder(const void* profile, size_t bytes)
{
    std::unique_ptr<Pair> pair(new Pair);
    const auto old = wh2_small_decoder_create(profile, bytes);
    pair->prototype = old.codec;
    const auto now = wirehair_small_decoder_create(profile, bytes);
    pair->integrated = now.codec;
    Require(old.status == Wh2Small_Success && now.status == WirehairSmall_Success && old.codec && now.codec,
            "standalone decoder");
    const auto* p = static_cast<const unsigned char*>(profile);
    uint64_t message = 0;
    for (unsigned i = 0; i < 8; ++i) message |= uint64_t(p[16 + i]) << (8 * i);
    Require(message <= SIZE_MAX - 2, "bounded decoder shape");
    pair->message = static_cast<size_t>(message);
    pair->scratch.reset(new unsigned char[pair->message + 2]);
    return {Wh2Small_Success, pair.release()};
}
Wh2SmallResult Encode(Wh2SmallCodec handle, uint32_t id, void* output, size_t capacity)
{
    auto* pair = static_cast<Pair*>(handle); pair->Prepare(capacity);
    const auto old = wh2_small_encode(pair->prototype, id, output, capacity);
    const auto now = wirehair_small_encode(pair->integrated, id, pair->scratch.get() + 1, capacity);
    pair->Check(old, now, output); return old;
}
Wh2SmallStatus Decode(Wh2SmallCodec handle, uint32_t id, const void* input, size_t bytes)
{
    auto* pair = static_cast<Pair*>(handle);
    const auto old = wh2_small_decode(pair->prototype, id, input, bytes);
    const auto now = wirehair_small_decode(pair->integrated, id, input, bytes);
    Require(int(old) == int(now), "prefix rank/status"); return old;
}
Wh2SmallResult Recover(Wh2SmallCodec handle, void* output, size_t capacity)
{
    auto* pair = static_cast<Pair*>(handle); pair->Prepare(capacity);
    const auto old = wh2_small_recover(pair->prototype, output, capacity);
    const auto now = wirehair_small_recover(pair->integrated, pair->scratch.get() + 1, capacity);
    pair->Check(old, now, output); return old;
}
void Free(Wh2SmallCodec handle) { delete static_cast<Pair*>(handle); }
} // namespace wh2_production_parity

#define WH2_SMALL_TEST_SERIALIZED 1
#define wh2_small_encoder_create wh2_production_parity::Encoder
#define wh2_small_decoder_create wh2_production_parity::Decoder
#define wh2_small_encode wh2_production_parity::Encode
#define wh2_small_decode wh2_production_parity::Decode
#define wh2_small_recover wh2_production_parity::Recover
#define wh2_small_free wh2_production_parity::Free
#include "Wh2SmallNativeTest.cpp"
