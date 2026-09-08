#include "../codec/WirehairSmallCore.h"
#include "../codec/WirehairK6Core.h"
#include "Wh2FrozenTrace.h"
#include "Wh2K3NativeData.inc"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {
bool tracking = false;
std::size_t allocations = 0, fail_at = SIZE_MAX, allocation_sizes[8] = {};
void* allocation_pointers[8] = {};
}
// Keep replacement allocation boundaries visible to the fault counter. GCC
// otherwise inlines delete's free into callers and diagnoses the intentional
// malloc-backed replacement new/delete pair as mismatched allocation.
#if defined(__GNUC__) || defined(__clang__)
#define WH2_TEST_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define WH2_TEST_NOINLINE __declspec(noinline)
#else
#define WH2_TEST_NOINLINE
#endif
WH2_TEST_NOINLINE void* operator new(std::size_t size)
{
    const std::size_t index = tracking ? allocations++ : SIZE_MAX;
    if (tracking && index == fail_at) throw std::bad_alloc();
    void* pointer = std::malloc(size ? size : 1);
    if (!pointer) throw std::bad_alloc();
    if (index < 8) { allocation_sizes[index] = size; allocation_pointers[index] = pointer; }
    return pointer;
}
WH2_TEST_NOINLINE void* operator new[](std::size_t n) { return ::operator new(n); }
WH2_TEST_NOINLINE void operator delete(void* p) noexcept { std::free(p); }
WH2_TEST_NOINLINE void operator delete[](void* p) noexcept { std::free(p); }
WH2_TEST_NOINLINE void* operator new(std::size_t n, const std::nothrow_t&) noexcept
{ try { return ::operator new(n); } catch (const std::bad_alloc&) { return nullptr; } }
WH2_TEST_NOINLINE void* operator new[](std::size_t n, const std::nothrow_t&) noexcept
{ try { return ::operator new[](n); } catch (const std::bad_alloc&) { return nullptr; } }
WH2_TEST_NOINLINE void operator delete(void* p, const std::nothrow_t&) noexcept { std::free(p); }
WH2_TEST_NOINLINE void operator delete[](void* p, const std::nothrow_t&) noexcept { std::free(p); }
#if defined(__cpp_sized_deallocation)
WH2_TEST_NOINLINE void operator delete(void* p, std::size_t) noexcept { std::free(p); }
WH2_TEST_NOINLINE void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

namespace {
namespace S = wirehair_small_core;
namespace Old = wirehair_k6_core;
using Byte = std::uint8_t;
alignas(64) const Byte kK6Lookup[] = {
#include "../codec/WirehairK6Lookup.inc"
};
std::uint64_t packet_checks = 0, feed_checks = 0, recover_checks = 0, row_checks = 0, cases = 0;

void Check(bool ok, const char* message)
{
    if (!ok) { std::cerr << "FAIL: " << message << '\n'; std::exit(1); }
}
void Start(std::size_t failure = SIZE_MAX)
{
    allocations = 0; fail_at = failure;
    std::fill(allocation_sizes, allocation_sizes + 8, 0);
    std::fill(allocation_pointers, allocation_pointers + 8, nullptr);
    tracking = true;
}
std::size_t Stop() { tracking = false; return allocations; }
template<class Fn> auto NoAlloc(Fn fn) -> decltype(fn())
{
    Start(0); const auto result = fn();
    Check(Stop() == 0, "packet operation allocated");
    return result;
}
Byte Mul(Byte a, Byte b)
{
    unsigned p = 0;
    for (unsigned i = 0; i < 8; ++i) if (b & (1u << i)) p ^= unsigned(a) << i;
    for (int i = 14; i >= 8; --i) if (p & (1u << i)) p ^= 0x14du << (i - 8);
    return static_cast<Byte>(p);
}
Byte Inverse(Byte value)
{
    Byte result = 1;
    for (unsigned i = 0; i < 254; ++i) result = Mul(result, value);
    return result;
}
template<unsigned K> struct Oracle {
    typedef std::array<Byte, K * K> Matrix;
    typedef std::array<Byte, K> Row;
    Matrix powers[2][32];
    static Matrix Product(const Matrix& a, const Matrix& b)
    {
        Matrix p = {};
        for (unsigned r = 0; r < K; ++r) for (unsigned c = 0; c < K; ++c)
            for (unsigned i = 0; i < K; ++i) p[r * K + c] ^= Mul(a[r * K + i], b[i * K + c]);
        return p;
    }
    Oracle()
    {
        const Byte small[6] = {8, 14, 7, 0, 0, 0};
        const Byte six[6] = {124, 127, 152, 84, 241, 63};
        for (unsigned phase = 0; phase < 2; ++phase) {
            powers[phase][0].fill(0);
            for (unsigned i = 0; i + 1 < K; ++i) powers[phase][0][(i + 1) * K + i] = 1;
            for (unsigned i = 0; i < K; ++i)
                powers[phase][0][i * K + K - 1] = (K == 3 ? small[i] : six[i]) ^ (i == 0 ? phase : 0);
        }
        for (unsigned level = 1; level < 32; ++level) {
            powers[0][level] = Product(powers[0][level - 1], powers[1][level - 1]);
            powers[1][level] = Product(powers[1][level - 1], powers[0][level - 1]);
        }
    }
    Row Coefficients(std::uint32_t id) const
    {
        Row row = {}; row[0] = 1;
        for (unsigned bit = 0; bit < 32; ++bit) if (id & (std::uint32_t(1) << bit)) {
            unsigned phase = 0;
            for (unsigned higher = bit + 1; higher < 32; ++higher) phase ^= (id >> higher) & 1u;
            Row next = {};
            for (unsigned r = 0; r < K; ++r) for (unsigned c = 0; c < K; ++c)
                next[r] ^= Mul(powers[phase][bit][r * K + c], row[c]);
            row = next;
        }
        return row;
    }
    std::vector<Byte> Packet(const std::vector<Byte>& source, unsigned B, std::uint32_t id) const
    {
        const Row row = Coefficients(id);
        std::vector<Byte> result(id == K - 1 ? source.size() - std::size_t(K - 1) * B : B, 0);
        for (std::size_t j = 0; j < result.size(); ++j) for (unsigned i = 0; i < K; ++i) {
            const std::size_t offset = std::size_t(i) * B + j;
            if (offset < source.size()) result[j] ^= Mul(row[i], source[offset]);
        }
        return result;
    }
    static unsigned Rank(std::vector<Row> rows)
    {
        unsigned rank = 0;
        for (unsigned col = 0; col < K && rank < rows.size(); ++col) {
            unsigned pivot = rank;
            while (pivot < rows.size() && !rows[pivot][col]) ++pivot;
            if (pivot == rows.size()) continue;
            std::swap(rows[rank], rows[pivot]);
            const Byte inverse = Inverse(rows[rank][col]);
            for (unsigned c = 0; c < K; ++c) rows[rank][c] = Mul(rows[rank][c], inverse);
            for (unsigned r = 0; r < rows.size(); ++r) if (r != rank) {
                const Byte factor = rows[r][col];
                for (unsigned c = 0; c < K; ++c) rows[r][c] ^= Mul(factor, rows[rank][c]);
            }
            ++rank;
        }
        return rank;
    }
};

std::vector<Byte> Message(std::size_t n)
{
    std::vector<Byte> result(n);
    for (std::size_t i = 0; i < n; ++i) result[i] = static_cast<Byte>(37 * i + i / 11);
    return result;
}
S::Lookup Lookup3() { return S::Lookup{wh2_k3_data::kLookup, sizeof(wh2_k3_data::kLookup)}; }
S::Lookup Lookup6() { return S::Lookup{kK6Lookup, sizeof(kK6Lookup)}; }
static_assert(sizeof(wh2_k3_data::kTraces) / sizeof(wh2_k3_data::kTraces[0]) == 6216, "trace roster");
static_assert(sizeof(wh2_k3_data::kHistory) / sizeof(wh2_k3_data::kHistory[0]) == 45, "history roster");
static_assert(sizeof(wh2_k3_data::kWindows) / sizeof(wh2_k3_data::kWindows[0]) == 43, "window roster");
static_assert(sizeof(wh2_k3_data::kRows) / sizeof(wh2_k3_data::kRows[0]) == 2226, "row roster");
void Match(S::Result result, Old::Result old)
{
    Check(static_cast<unsigned>(result.status) == static_cast<unsigned>(old.status) &&
        result.bytes_required == old.bytes_required && result.bytes_written == old.bytes_written, "K6 result parity");
}

template<unsigned K> void Exercise(S::Lookup lookup, const Oracle<K>& oracle, unsigned B, unsigned tail,
                                   const std::vector<std::uint32_t>& ids, const unsigned* expected = nullptr,
                                   bool stress = false)
{
    typedef S::Encoder<K> Encoder;
    typedef S::Decoder<K> Decoder;
    const std::vector<Byte> source = Message(std::size_t(K - 1) * B + tail);
    std::unique_ptr<Encoder> encoder;
    std::unique_ptr<Decoder> decoder;
    Check(Encoder::Create(lookup, source.data(), source.size(), B, encoder) == S::Status::Success, "encoder create");
    Check(Decoder::Create(lookup, source.size(), B, decoder) == S::Status::Success, "decoder create");
    std::unique_ptr<Old::Encoder> old_encoder;
    std::unique_ptr<Old::Decoder> old_decoder;
    if (K == 6) {
        const Old::Lookup old_lookup = {lookup.data, lookup.bytes};
        Check(Old::Encoder::Create(old_lookup, source.data(), source.size(), B, old_encoder) == Old::Status::Success,
              "original K6 encoder");
        Check(Old::Decoder::Create(old_lookup, source.size(), B, old_decoder) == Old::Status::Success, "original K6 decoder");
    }
    std::vector<typename Oracle<K>::Row> rows;
    std::vector<Byte> output(source.size() + 2, 0xa5), old_output(source.size() + 2, 0xa5);
    Check(NoAlloc([&] { return decoder->Recover(output.data() + 1, source.size()); }).status == S::Status::NeedMore,
          "empty recovery");
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const auto packet = oracle.Packet(source, B, ids[i]);
        std::vector<Byte> generated(packet.size() + 2, 0xa5), old_packet(packet.size() + 2, 0xa5);
        if (stress) {
            const auto short_result = NoAlloc([&] { return encoder->Encode(ids[i], generated.data() + 1, packet.size() - 1); });
            Check(short_result.status == S::Status::BufferTooSmall && !short_result.bytes_written &&
                short_result.bytes_required == packet.size() && std::count(generated.begin(), generated.end(), 0xa5) ==
                static_cast<std::ptrdiff_t>(generated.size()), "short Encode preserves output");
        }
        const auto encoded = NoAlloc([&] { return encoder->Encode(ids[i], generated.data() + 1, packet.size()); });
        Check(encoded.status == S::Status::Success && encoded.bytes_required == packet.size() &&
            encoded.bytes_written == packet.size() && std::equal(packet.begin(), packet.end(), generated.begin() + 1) &&
            generated.front() == 0xa5 && generated.back() == 0xa5, "independent packet bytes/guards");
        ++packet_checks;
        if (K == 6) {
            Match(encoded, NoAlloc([&] { return old_encoder->Encode(ids[i], old_packet.data() + 1, packet.size()); }));
            Check(old_packet == generated, "original K6 packet bytes");
        }
        rows.push_back(oracle.Coefficients(ids[i]));
        const unsigned rank = Oracle<K>::Rank(rows);
        if (expected && i + 1 >= K) Check(rank == expected[i + 1 - K], "recorded prefix rank");
        const auto fed = NoAlloc([&] { return decoder->Feed(ids[i], generated.data() + 1, packet.size()); });
        Check(fed.status == (rank == K ? S::Status::Success : S::Status::NeedMore) && decoder->Rank() == rank,
              "native rank/status");
        ++feed_checks;
        if (K == 6) {
            Match(fed, NoAlloc([&] { return old_decoder->Feed(ids[i], generated.data() + 1, packet.size()); }));
            Check(old_decoder->Rank() == rank, "original K6 rank");
        }
        if (stress) {
            Check(NoAlloc([&] { return decoder->Feed(ids[i], packet.data(), packet.size()); }).status == fed.status,
                  "duplicate idempotence");
            std::vector<Byte> bad = packet; bad.back() ^= 1;
            Check(NoAlloc([&] { return decoder->Feed(ids[i], bad.data(), bad.size()); }).status == S::Status::Conflict &&
                  decoder->Rank() == rank, "dependent conflict preserves basis");
            if (K == 6) {
                Check(NoAlloc([&] { return old_decoder->Feed(ids[i], bad.data(), bad.size()); }).status == Old::Status::Conflict,
                      "original K6 conflict parity");
            }
            Check(NoAlloc([&] { return decoder->Recover(output.data() + 1, source.size() - 1); }).status ==
                  S::Status::BufferTooSmall, "short recovery");
        }
        const auto recovered = NoAlloc([&] { return decoder->Recover(output.data() + 1, source.size()); });
        ++recover_checks;
        Check(recovered.status == fed.status, "Recover rank status");
        if (rank == K) Check(std::equal(source.begin(), source.end(), output.begin() + 1), "exact recovered source");
        else Check(std::count(output.begin(), output.end(), 0xa5) == static_cast<std::ptrdiff_t>(output.size()),
                   "NeedMore leaves output unchanged");
        Check(output.front() == 0xa5 && output.back() == 0xa5, "recovery guards");
        if (K == 6) {
            Match(recovered, NoAlloc([&] { return old_decoder->Recover(old_output.data() + 1, source.size()); }));
            Check(old_output == output, "original K6 recovery parity");
        }
    }
    Check(decoder->Rank() == K, "fixture eventually recovers");
    Check(NoAlloc([&] { return decoder->Recover(output.data() + 1, source.size()); }).status == S::Status::Success,
          "repeat recovery");
    ++cases;
}

template<unsigned K> void InvalidAndAllocations(S::Lookup lookup)
{
    typedef S::Encoder<K> E;
    typedef S::Decoder<K> D;
    const std::size_t message = (K - 1) * 64 + 1;
    auto source = Message(K * 64);
    std::unique_ptr<E> encoder;
    std::unique_ptr<D> decoder;
    for (std::uint64_t bad : {std::uint64_t(0), std::uint64_t((K - 1) * 64), std::uint64_t(K * 64 + 1), UINT64_MAX}) {
        Check(E::Create(lookup, source.data(), bad, 64, encoder) == S::Status::InvalidInput && !encoder, "encoder bad shape");
        Check(D::Create(lookup, bad, 64, decoder) == S::Status::InvalidInput && !decoder, "decoder bad shape");
    }
    Check(E::Create(lookup, nullptr, message, 64, encoder) == S::Status::InvalidInput, "null source");
    Check(E::Create(lookup, &encoder, K, 1, encoder) == S::Status::InvalidInput && !encoder, "handle/source alias");
    Check(D::Create(lookup, K, 0, decoder) == S::Status::InvalidInput, "zero block");
    Check(D::Create(lookup, std::uint64_t(K) * 0x80000000u, 0x80000000u, decoder) == S::Status::InvalidInput, "signed byte cap");
    const void* wrap = reinterpret_cast<const void*>(std::numeric_limits<std::uintptr_t>::max() - 1);
    Check(E::Create(lookup, wrap, message, 64, encoder) == S::Status::InvalidInput, "wrapped source");
    S::Lookup bad_lookup = {lookup.data, lookup.bytes - 1};
    Check(D::Create(bad_lookup, message, 64, decoder) == S::Status::InvalidInput, "short lookup");
    Byte row[K]; std::fill(row, row + K, 0xa5);
    Check(S::Row<K>(bad_lookup, 0, row) == S::Status::InvalidInput && std::count(row, row + K, 0xa5) == K,
          "invalid row leaves output unchanged");
    Check(S::Row<K>(lookup, 0, const_cast<Byte*>(lookup.data)) == S::Status::InvalidInput, "row/lookup alias");
    Check(S::Row<K>(lookup, 0, nullptr) == S::Status::InvalidInput, "null row output");
    auto bad_table = std::vector<Byte>(lookup.data, lookup.data + lookup.bytes); bad_table[0] = 0;
    Check(E::Create({bad_table.data(), bad_table.size()}, source.data(), message, 64, encoder) ==
          S::Status::InvalidInput, "nonsystematic lookup");
    for (std::size_t failure = 0; failure < 2; ++failure) {
        Start(failure); auto status = E::Create(lookup, source.data(), message, 64, encoder); const auto ec = Stop();
        Check(status == S::Status::OutOfMemory && !encoder && ec == failure + 1, "encoder OOM transactional");
        Start(failure); status = D::Create(lookup, message, 64, decoder); const auto dc = Stop();
        Check(status == S::Status::OutOfMemory && !decoder && dc == failure + 1, "decoder OOM transactional");
    }
    Start(); auto status = E::Create(lookup, source.data(), K * 64, 64, encoder); auto count = Stop();
    Check(status == S::Status::Success && count == 1 && allocation_sizes[0] == sizeof(E), "full encoder storage");
    encoder.reset();
    Start(); status = E::Create(lookup, source.data(), message, 64, encoder); count = Stop();
    Check(status == S::Status::Success && count == 2 && allocation_sizes[1] == 64, "partial encoder storage");
    Check(encoder->Encode(K, allocation_pointers[1], 64).status == S::Status::InvalidInput, "padding alias");
    E* kept = encoder.get();
    Check(E::Create(lookup, source.data(), message, 64, encoder) == S::Status::InvalidInput && encoder.get() == kept,
          "nonempty handle unchanged");
    for (void* output : {static_cast<void*>(source.data()), static_cast<void*>(encoder.get()),
                        static_cast<void*>(const_cast<Byte*>(lookup.data))})
        Check(encoder->Encode(K, output, 64).status == S::Status::InvalidInput, "encoder output alias");
    Byte output[K * 64] = {};
    Check(encoder->Encode(K, output, SIZE_MAX).status == S::Status::InvalidInput, "wrapped output");
    Check(encoder->Encode(K - 1, nullptr, 0).status == S::Status::BufferTooSmall, "tail size query");
    Check(encoder->Encode(K - 1, nullptr, 1).status == S::Status::InvalidInput, "null output");
    Start(); status = D::Create(lookup, message, 64, decoder); count = Stop();
    Check(status == S::Status::Success && count == 2 && allocation_sizes[0] == sizeof(D) &&
          allocation_sizes[1] == (K + 1) * 64, "decoder storage");
    void* slab = allocation_pointers[1];
    Check(decoder->Feed(0, slab, 64).status == S::Status::InvalidInput, "slab input alias");
    Check(decoder->Recover(slab, message).status == S::Status::InvalidInput, "slab output alias");
    Check(decoder->Feed(0, decoder.get(), 64).status == S::Status::InvalidInput, "metadata input alias");
    Check(decoder->Recover(decoder.get(), message).status == S::Status::InvalidInput, "metadata output alias");
    Check(decoder->Recover(const_cast<Byte*>(lookup.data), message).status == S::Status::InvalidInput, "lookup output alias");
    Check(decoder->Feed(0, wrap, 64).status == S::Status::InvalidInput, "wrapped input");
    Check(decoder->Feed(0, nullptr, 64).status == S::Status::InvalidInput, "null input");
    Check(decoder->Feed(K - 1, output, 64).status == S::Status::InvalidInput, "oversized partial input");
    Check(decoder->Feed(K, output, 63).status == S::Status::InvalidInput && decoder->Rank() == 0, "short repair input");
    Check(source == Message(K * 64), "invalid operations preserve source");
}

template<unsigned K> void DependencyTests(S::Lookup lookup, const Oracle<K>& oracle)
{
    const unsigned B = 64;
    const auto source = Message(K * B);
    auto modified = std::vector<Byte>(lookup.data, lookup.data + lookup.bytes);
    std::copy(modified.begin(), modified.begin() + K, modified.begin() + 10 * K);
    const S::Lookup custom = {modified.data(), modified.size()};
    std::unique_ptr<S::Decoder<K>> decoder;
    Check(S::Decoder<K>::Create(custom, source.size(), B, decoder) == S::Status::Success, "dependent fixture create");
    const auto first = oracle.Packet(source, B, 0);
    Check(NoAlloc([&] { return decoder->Feed(0, first.data(), B); }).status == S::Status::NeedMore, "first unit row");
    Check(NoAlloc([&] { return decoder->Feed(10, first.data(), B); }).status == S::Status::NeedMore && decoder->Rank() == 1,
          "distinct dependent ID before full rank");
    auto bad = first; bad[31] ^= 1;
    Check(NoAlloc([&] { return decoder->Feed(10, bad.data(), B); }).status == S::Status::Conflict && decoder->Rank() == 1,
          "distinct dependent contradiction before full rank");
    for (unsigned id = 1; id < K; ++id) {
        const auto packet = oracle.Packet(source, B, id);
        Check(NoAlloc([&] { return decoder->Feed(id, packet.data(), B); }).status ==
              (id + 1 == K ? S::Status::Success : S::Status::NeedMore), "resume after deficient conflict");
    }
    auto output = std::vector<Byte>(source.size(), 0);
    Check(NoAlloc([&] { return decoder->Recover(output.data(), output.size()); }).status == S::Status::Success &&
          output == source, "dependent fixture recovered");
    decoder.reset();
    // Fresh decoder, no live encoder: all payloads come from the independent
    // polynomial oracle. Leave the nontrivial repair echelon unsolved until
    // after a novel dependent equation and a contradictory copy.
    Check(S::Decoder<K>::Create(lookup, source.size(), B, decoder) == S::Status::Success, "deferred solve create");
    for (unsigned i = 0; i < K; ++i) {
        const auto packet = oracle.Packet(source, B, K + i);
        Check(NoAlloc([&] { return decoder->Feed(K + i, packet.data(), B); }).status ==
              (i + 1 == K ? S::Status::Success : S::Status::NeedMore), "repair echelon rank");
    }
    const auto distant = oracle.Packet(source, B, UINT32_MAX);
    Check(NoAlloc([&] { return decoder->Feed(UINT32_MAX, distant.data(), B); }).status == S::Status::Success,
          "novel dependent equation before solve");
    bad = distant; bad[17] ^= 1;
    Check(NoAlloc([&] { return decoder->Feed(UINT32_MAX, bad.data(), B); }).status == S::Status::Conflict &&
          decoder->Rank() == K, "unsolved contradiction preserves echelon");
    Check(NoAlloc([&] { return decoder->Recover(output.data(), output.size()); }).status == S::Status::Success &&
          output == source, "deferred solve after conflict");
}

template<unsigned K> void Neutral(S::Lookup lookup, const Oracle<K>& oracle)
{
    InvalidAndAllocations<K>(lookup);
    DependencyTests<K>(lookup, oracle);
    // Exercise every packed selector and both lower-table phases, including
    // mixed-bit chunks. Recovery traces alone do not cover this address space.
    std::vector<std::uint32_t> selectors;
    for (unsigned phase = 0; phase < 2; ++phase) {
        for (unsigned value = 0; value < 1024; ++value) selectors.push_back((phase << 10) | value);
        for (unsigned value = 0; value < 128; ++value) {
            selectors.push_back((phase << 17) | (value << 10) | 513u);
            selectors.push_back((phase << 24) | (value << 17) | (85u << 10) | 341u);
        }
    }
    for (unsigned value = 0; value < 256; ++value)
        selectors.push_back((value << 24) | (65u << 17) | (42u << 10) | 777u);
    Check(selectors.size() == 2816, "selector coverage count");
    for (std::uint32_t id : selectors) {
        std::array<Byte, K + 2> actual; actual.fill(0xa5);
        const auto reference = oracle.Coefficients(id);
        Check(S::Row<K>(lookup, id, actual.data() + 1) == S::Status::Success &&
              std::equal(reference.begin(), reference.end(), actual.begin() + 1) &&
              actual.front() == 0xa5 && actual.back() == 0xa5, "complete selector polynomial oracle");
        if (K == 6) {
            Byte old[6] = {};
            Check(Old::Row({lookup.data, lookup.bytes}, id, old) == Old::Status::Success &&
                  std::equal(reference.begin(), reference.end(), old), "K6 selector parity");
        }
        ++row_checks;
    }
    const unsigned widths[] = {1,2,3,4,7,16,31,32,63,64,65,127,128,129,255,256,257,1280,4096};
    for (unsigned B : widths) {
        std::vector<unsigned> tails = {B};
        if (B > 1) tails.push_back(1);
        if (B > 2) tails.push_back(B - 1);
        for (unsigned tail : tails) for (unsigned order = 0; order < 3; ++order) {
            std::vector<std::uint32_t> ids;
            for (unsigned i = 0; i < K; ++i) ids.push_back(order == 0 ? K - 1 - i : order == 1 ? K + i : UINT32_MAX - 2 * i);
            for (unsigned bit : {10u,17u,24u,31u}) {
                ids.push_back((std::uint32_t(1) << bit) - 1); ids.push_back(std::uint32_t(1) << bit);
            }
            ids.push_back(UINT32_MAX);
            Exercise<K>(lookup, oracle, B, tail, ids, nullptr, true);
        }
    }
}

void Corpus(const Oracle<3>& oracle)
{
    for (const auto& row : wh2_k3_data::kRows) {
        Byte mapped[5] = {0xa5,0,0,0,0xa5};
        Check(S::Row<3>(Lookup3(), row.id, mapped + 1) == S::Status::Success, "native recorded row");
        const auto expected = oracle.Coefficients(row.id);
        Check(std::equal(expected.begin(), expected.end(), mapped + 1) &&
              std::equal(expected.begin(), expected.end(), row.values) && mapped[0] == 0xa5 && mapped[4] == 0xa5,
              "three-way row oracle");
        ++row_checks;
    }
    for (const auto& trace : wh2_k3_data::kTraces)
        Exercise<3>(Lookup3(), oracle, trace.B, trace.B, std::vector<std::uint32_t>(trace.ids, trace.ids + 7), trace.ranks);
    const unsigned widths[] = {2,64,1280};
    for (const auto& prefix : wh2_k3_data::kHistory) for (unsigned bit = 0; bit < 3; ++bit) if (prefix.widths & (1u << bit))
        Exercise<3>(Lookup3(), oracle, widths[bit], widths[bit], std::vector<std::uint32_t>(prefix.ids, prefix.ids + prefix.count));
    for (const auto& window : wh2_k3_data::kWindows)
        for (unsigned a = 0; a < 5; ++a) for (unsigned b = a + 1; b < 6; ++b) for (unsigned c = b + 1; c < 7; ++c)
            Exercise<3>(Lookup3(), oracle, 2, 2, {window[a],window[b],window[c]});
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2 || (std::strcmp(argv[1], "--neutral") && std::strcmp(argv[1], "--corpus"))) return 2;
    Check(gf256_init() == 0, "shared GF runtime");
#if defined(WH2_SMALL_EXPECT_PORTABLE)
    gf256_x86_cpu_features features = {};
    gf256_get_active_x86_cpu_features(&features);
    Check(!features.SSSE3 && !features.AVX2 && !features.GFNI && !features.AVX512 &&
          !wirehair_k6_payload::Available(), "portable backend is active");
#endif
    Check(wirehair::wh2_benchmark::Sha256Hex(wh2_k3_data::kLookup, sizeof(wh2_k3_data::kLookup)) ==
          wh2_k3_data::kLookupSha, "native K3 lookup SHA");
    const Oracle<3> three;
    if (!std::strcmp(argv[1], "--neutral")) {
        for (unsigned a = 0; a < 256; ++a) for (unsigned b = 0; b < 256; ++b)
            Check(gf256_mul(static_cast<Byte>(a), static_cast<Byte>(b)) == Mul(static_cast<Byte>(a), static_cast<Byte>(b)),
                  "shared field oracle");
        Neutral<3>(Lookup3(), three);
        const Oracle<6> six;
        Neutral<6>(Lookup6(), six);
        Check(cases == 324 && packet_checks == 4374 && row_checks == 5632, "neutral count accounting");
    } else {
        Corpus(three);
        Check(cases == 7774 && packet_checks == 48187 && row_checks == 2226, "corpus count accounting");
    }
    std::cout << "PASS " << argv[1] << " cases=" << cases << " packet_oracles=" << packet_checks
              << " feeds=" << feed_checks << " recoveries=" << recover_checks << " rows=" << row_checks
              << " raw=" << wh2_k3_data::kRawSha << " GFNI=" << wirehair_k6_payload::Available() << '\n';
    return 0;
}
