#include "Wh2ThueMorseMatvecCostR0Bridge.h"
#include "Wh2ThueMorseMatvecGfniR0.h"
#include "wirehair/wirehair.h"
#include "gf256.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
namespace B = wh2_matvec_cost_r0;
typedef std::uint8_t Byte;

void Check(bool value, const char* message)
{
    if (!value) throw std::runtime_error(message);
}

void Exact(const B::Result& result, int code, std::uint64_t required,
           std::uint64_t written, unsigned length_kind)
{
    Check(result.code == code && result.status == (code == 0 ? B::Success :
          code == 1 ? B::NeedMore : B::Failure) && result.required == required &&
          result.written == written && result.length_kind == length_kind,
          "exact bridge status, raw code and length provenance");
}

void Rejected(const B::Result& result)
{
    Exact(result, -1, 0, 0, B::NoLength);
}

bool SameState(const B::State& a, const B::State& b)
{
    return a.handle == b.handle && a.source == b.source && a.lookup == b.lookup &&
        a.lookup_bytes == b.lookup_bytes && a.message_bytes == b.message_bytes &&
        a.block_bytes == b.block_bytes && a.arm == b.arm && a.kind == b.kind &&
        a.profile_bytes == b.profile_bytes && std::equal(a.profile, a.profile + 32, b.profile);
}

std::vector<Byte> NeutralLookup()
{
    // Not the selected construction: all matrices identity, direct rows cyclic.
    std::vector<Byte> result(39936, 0);
    for (unsigned phase = 0; phase < 2; ++phase)
        for (unsigned id = 0; id < 1024; ++id)
            result[phase * 6144u + id * 6u + id % 6u] = 1;
    for (std::size_t offset = 12288; offset < result.size(); offset += 36)
        for (unsigned k = 0; k < 6; ++k) result[offset + 7u * k] = 1;
    return result;
}

struct Owned {
    const B::Api* api;
    B::State state;
    explicit Owned(const B::Api* value) : api(value), state{} {}
    ~Owned() {
        if (state.kind == 1) api->encoder_free(&state);
        if (state.kind == 2) api->decoder_free(&state);
    }
    Owned(const Owned&) = delete;
    Owned& operator=(const Owned&) = delete;
};

bool Guards(const std::vector<Byte>& data, std::size_t bytes)
{
    return std::all_of(data.begin(), data.begin() + 17, [](Byte b) { return b == 0xa7; }) &&
        std::all_of(data.begin() + 17 + bytes, data.end(), [](Byte b) { return b == 0xa7; });
}

void Case(const std::array<const B::Api*, 4>& apis, unsigned arm, unsigned block,
          const std::vector<Byte>& lookup)
{
    const B::Api& api = *apis[arm];
    const std::size_t message_bytes = 6u * block;
    std::vector<Byte> source(message_bytes + 34, 0xa7);
    for (std::size_t i = 0; i < message_bytes; ++i) source[17 + i] = static_cast<Byte>(19u * i + 7u);
    const auto saved_source = source;
    const auto saved_lookup = lookup;
    Owned encoder(&api);
    Rejected(api.create(lookup.data(), lookup.size(), source.data() + 17, message_bytes, block, nullptr));
    Rejected(api.create(lookup.data(), lookup.size(), source.data() + 17, message_bytes - block, block,
                        &encoder.state));
    Check(encoder.state.handle == nullptr && encoder.state.kind == 0, "failed create preserves empty state");
    Exact(api.create(lookup.data(), lookup.size(), source.data() + 17, message_bytes, block,
                     &encoder.state), 0, 0, 0, B::NoLength);
    Check(api.valid_profile(&encoder.state) && encoder.state.source == source.data() + 17 &&
          encoder.state.arm == arm && encoder.state.kind == 1, "live encoder metadata and emitted profile");
    const B::State saved_encoder = encoder.state;
    Rejected(api.create(lookup.data(), lookup.size(), source.data() + 17, message_bytes, block,
                        &encoder.state));
    Check(SameState(encoder.state, saved_encoder), "failed reuse preserves entire encoder state");

    std::array<std::uint32_t, 18> ids = {};
    for (unsigned j = 0; j < 12; ++j) ids[j] = j;
    for (unsigned j = 0; j < 6; ++j) ids[12 + j] = UINT32_MAX - 2u * j;
    std::array<std::vector<Byte>, 18> packets;
    for (unsigned j = 0; j < 18; ++j) {
        packets[j].assign(block + 34u, 0xa7);
        Exact(api.encode(&encoder.state, ids[j], packets[j].data() + 17, block),
              0, block, block, B::ReturnedLength);
        Check(Guards(packets[j], block), "encode guards");
        if (arm < 2 || j < 6) {
            const unsigned column = j < 6 ? j : (ids[j] & 1023u) % 6u;
            Check(std::equal(packets[j].begin() + 17, packets[j].begin() + 17 + block,
                             source.begin() + 17 + column * block), "independent neutral packet oracle");
        }
        if (arm < 2) {
            std::array<Byte, 40> row;
            row.fill(0xa7);
            Exact(api.row(lookup.data(), lookup.size(), ids[j], row.data() + 17),
                  0, 6, 6, B::InferredLength);
            for (unsigned k = 0; k < row.size(); ++k) {
                const Byte expected = k >= 17 && k < 23 ?
                    static_cast<Byte>(k - 17 == (ids[j] & 1023u) % 6u) : 0xa7;
                Check(row[k] == expected, "independent neutral row and guards");
            }
        }
    }
    Check((api.row != nullptr) == (arm < 2), "native-only row endpoint");
    const auto saved_packets = packets;
    std::vector<Byte> scratch(message_bytes + 34, 0xa7);
    Rejected(api.encode(nullptr, 6, scratch.data() + 17, block));
    Rejected(api.encode(&encoder.state, 6, source.data() + 17, block));
    Rejected(api.encode(&encoder.state, 6, &encoder.state, block));
    Rejected(apis[(arm + 1) % 4]->encode(&encoder.state, 6, scratch.data() + 17, block));
    apis[(arm + 1) % 4]->encoder_free(&encoder.state);
    api.decoder_free(&encoder.state);
    Check(encoder.state.handle == saved_encoder.handle, "wrong owner/kind free preserves handle");
    const B::Result short_packet = api.encode(&encoder.state, 6, scratch.data() + 17, block - 1u);
    Exact(short_packet, arm == 2 ? 2 : 3, arm == 2 ? 0u : block, 0, B::ReturnedLength);
    Check(short_packet.status == B::Failure && short_packet.code != 0 && short_packet.written == 0 &&
          std::all_of(scratch.begin(), scratch.end(), [](Byte b) { return b == 0xa7; }),
          "short packet capacity leaves entire output untouched");

    for (unsigned distant = 0; distant < 2; ++distant) {
        Owned decoder(&api);
        Rejected(api.decoder_create(&encoder.state, &encoder.state));
        Exact(api.decoder_create(&encoder.state, &decoder.state), 0, 0, 0, B::NoLength);
        Check(api.valid_profile(&decoder.state) && decoder.state.source == nullptr &&
              decoder.state.kind == 2 && decoder.state.profile_bytes == encoder.state.profile_bytes &&
              std::equal(decoder.state.profile, decoder.state.profile + 32, encoder.state.profile),
              "decoder captures exact profile without source borrowing");
        const void* decoder_handle = decoder.state.handle;
        Rejected(api.encode(&decoder.state, 6, scratch.data() + 17, block));
        Rejected(api.feed(&encoder.state, 6, packets[6].data() + 17, block));
        Rejected(apis[(arm + 1) % 4]->feed(&decoder.state, 6, packets[6].data() + 17, block));
        Rejected(api.feed(&decoder.state, 6, nullptr, block));
        Rejected(api.recover(&decoder.state, &decoder.state, message_bytes));
        std::fill(scratch.begin(), scratch.end(), 0xa7);
        Exact(api.recover(&decoder.state, scratch.data() + 17, message_bytes),
              1, message_bytes, 0, arm == 2 ? B::InferredLength : B::ReturnedLength);
        Check(std::all_of(scratch.begin(), scratch.end(), [](Byte b) { return b == 0xa7; }),
              "NeedMore recovery changes no bytes");
        bool complete = false;
        unsigned fed = 0;
        for (unsigned position = 0; position < 12; ++position) {
            const unsigned index = position < 6 ? (distant ? 12u : 6u) + position : position - 6u;
            const B::Result result = api.feed(&decoder.state, ids[index], packets[index].data() + 17, block);
            Check(result.code == 0 || result.code == 1, "fixed neutral receive has no codec error");
            Exact(result, result.code, block, 0, arm < 2 ? B::ReturnedLength : B::InferredLength);
            ++fed;
            if (result.status == B::Success) { complete = true; break; }
        }
        Check(complete && fed >= 6 && fed <= 12, "first-success within fixed 12-packet cap");
        std::fill(scratch.begin(), scratch.end(), 0xa7);
        const B::Result short_recovery = api.recover(&decoder.state, scratch.data() + 17, message_bytes - 1u);
        Exact(short_recovery, arm == 2 ? -1 : 3, message_bytes, 0,
              arm == 2 ? B::InferredLength : B::ReturnedLength);
        Check(short_recovery.status == B::Failure && short_recovery.written == 0 &&
              std::all_of(scratch.begin(), scratch.end(), [](Byte b) { return b == 0xa7; }),
              "short recovery capacity leaves all bytes untouched");
        Exact(api.recover(&decoder.state, scratch.data() + 17, message_bytes),
              0, message_bytes, message_bytes, arm == 2 ? B::InferredLength : B::ReturnedLength);
        Check(scratch == source, "full exact recovered message and both guards");
        Check(decoder.state.handle == decoder_handle, "receiver state handle remains stable");
        api.encoder_free(&decoder.state);
        Check(decoder.state.handle == decoder_handle, "wrong-kind free cannot release decoder");
        api.decoder_free(&decoder.state);
        Check(!decoder.state.handle && decoder.state.kind == 0 && decoder.state.profile_bytes == 0,
              "decoder free clears caller state");
        api.decoder_free(&decoder.state);
    }
    Check(source == saved_source && lookup == saved_lookup && packets == saved_packets,
          "complete source, lookup and packet corpus immutable");
    Check(SameState(encoder.state, saved_encoder) && api.valid_profile(&encoder.state),
          "all negative calls preserve encoder profile");
    if (arm == 3) {
        encoder.state.profile[0] ^= 1u;
        Check(!api.valid_profile(&encoder.state), "corrupt serialized profile fails validation");
        encoder.state.profile[0] ^= 1u;
    }
    api.encoder_free(&encoder.state);
    Check(!encoder.state.handle && encoder.state.kind == 0 && encoder.state.profile_bytes == 0,
          "encoder free clears caller state");
    api.encoder_free(&encoder.state);
    api.encoder_free(nullptr);
    api.decoder_free(nullptr);
    Check(source == saved_source && lookup == saved_lookup && packets == saved_packets,
          "codec destruction preserves all borrowed bytes");
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2 || (std::strcmp(argv[1], "--require-gfni") != 0 &&
                      std::strcmp(argv[1], "--expect-scalar") != 0)) return 2;
    try {
        Check(wirehair_init() == Wirehair_Success && GF256Ctx.Polynomial == 0x14d,
              "one qualified shared GF runtime");
        Check(wh2_matvec_gfni_r0::Available() == (std::strcmp(argv[1], "--require-gfni") == 0),
              "required GFNI or scalar availability");
        const std::array<const B::Api*, 4> apis = {{wh2_matvec_cost_baseline_api(),
            wh2_matvec_cost_candidate_api(), wh2_matvec_cost_wh1_api(), wh2_matvec_cost_public_api()}};
        const std::vector<Byte> lookup = NeutralLookup();
        for (unsigned block : {2u, 65u})
            for (unsigned arm = 0; arm < apis.size(); ++arm) Case(apis, arm, block, lookup);
        std::cout << "PASS eight neutral four-arm B2/B65 cases, 16 bounded receive lifetimes; no selected data/timing\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
