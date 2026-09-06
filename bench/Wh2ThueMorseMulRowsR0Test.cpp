// Neutral algebra/lookup tests only: no fixed pair, generated header or clock.
#include "Wh2ThueMorseMulRowsR0.h"
#include "gf256.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
namespace N = wh2_thue_native;
namespace C = wh2_thue_mulrows_r0;
using Row = std::array<std::uint8_t, 6>;
const std::uint8_t kGuard = 0xa5;
std::uint8_t products[256][256];
std::size_t checked_rows = 0;

void Require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

std::uint8_t Product(unsigned a, unsigned b)
{
    unsigned polynomial = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if ((b >> bit) & 1u) polynomial ^= a << bit;
    for (int bit = 14; bit >= 8; --bit)
        if ((polynomial >> bit) & 1u) polynomial ^= 0x14du << (bit - 8);
    return static_cast<std::uint8_t>(polynomial);
}

void Field()
{
    Require(gf256_init() == 0 && GF256Ctx.Polynomial == 0x14d, "GF256 initialization");
    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b) {
            products[a][b] = Product(a, b);
            Require(gf256_mul(static_cast<std::uint8_t>(a), static_cast<std::uint8_t>(b)) == products[a][b] &&
                    gf256_mul(static_cast<std::uint8_t>(b), static_cast<std::uint8_t>(a)) == products[a][b],
                    "exhaustive independent multiplication/commutativity");
        }
}

std::uint64_t Random(std::uint64_t& state)
{
    std::uint64_t value = (state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

struct NeutralLookup {
    std::vector<std::uint8_t> storage;
    NeutralLookup() : storage(N::kLookupBytes + 16, kGuard) {
        std::uint64_t state = UINT64_C(0x4d554c524f575354);
        for (std::size_t i = 0; i < N::kLookupBytes; ++i)
            Data()[i] = static_cast<std::uint8_t>(Random(state));
        // Identity zero-groups make the independent oracle apply all groups,
        // whereas both implementations deliberately skip zero-index groups.
        const unsigned identity_offsets[] = {0, 6144, 12288, 16896, 21504, 26112, 30720};
        for (unsigned offset : identity_offsets)
            for (unsigned r = 0; r < 6; ++r)
                for (unsigned c = 0; c < 6; ++c)
                    Data()[offset + r * 6 + c] = r == c ? 1 : 0;
    }
    std::uint8_t* Data() { return storage.data() + 8; }
    N::Lookup Get() { return N::Lookup{Data(), N::kLookupBytes}; }
};

unsigned Parity(std::uint32_t value)
{
    unsigned ones = 0;
    for (unsigned bit = 0; bit < 32; ++bit) ones += (value >> bit) & 1u;
    return ones % 2;
}

Row Reference(N::Lookup lookup, std::uint32_t id)
{
    Row vector;
    const unsigned low_phase = Parity(id / 1024);
    const unsigned low = id % 1024;
    std::copy(lookup.data + low_phase * 6144 + low * 6,
              lookup.data + low_phase * 6144 + low * 6 + 6, vector.begin());
    const unsigned starts[] = {10, 17, 24};
    const unsigned widths[] = {7, 7, 8};
    const unsigned offsets[] = {12288, 21504, 30720};
    for (unsigned group = 0; group < 3; ++group) {
        const unsigned digit = (id >> starts[group]) % (1u << widths[group]);
        const unsigned phase = group == 2 ? 0 : Parity(id >> (starts[group] + widths[group]));
        const std::uint8_t* matrix = lookup.data + offsets[group] + phase * 4608 + digit * 36;
        Row next = {{0}};
        for (unsigned c = 0; c < 6; ++c)
            for (unsigned r = 0; r < 6; ++r)
                next[r] ^= products[matrix[r * 6 + c]][vector[c]];
        vector = next;
    }
    return vector;
}

void Check(N::Lookup lookup, std::uint32_t id)
{
    std::array<std::uint8_t, 22> old_output, new_output;
    old_output.fill(kGuard); new_output.fill(kGuard);
    const Row expected = Reference(lookup, id);
    Require(N::Row(lookup, id, old_output.data() + 8) == N::Status::Success &&
            C::Row(lookup, id, new_output.data() + 8) == N::Status::Success, "valid row status");
    Require(old_output == new_output && std::equal(expected.begin(), expected.end(), new_output.begin() + 8),
            "independent row bytes");
    for (unsigned i = 0; i < 22; ++i)
        if (i < 8 || i >= 14) Require(new_output[i] == kGuard, "row output guard");
    ++checked_rows;
}

void Ids()
{
    NeutralLookup table;
    const auto original = table.storage;
    for (std::uint32_t id = 0; id < 1024; ++id) Check(table.Get(), id);
    for (std::uint32_t i = 0; i < 65536; ++i)
        // Odd multiplication is a permutation of uint32: all65536 are distinct.
        Check(table.Get(), i * UINT32_C(0x9e3779b9) + UINT32_C(0x73a91e65));
    for (unsigned bit = 0; bit < 32; ++bit) {
        const std::uint32_t seam = std::uint32_t(1) << bit;
        Check(table.Get(), seam - 1); Check(table.Get(), seam); Check(table.Get(), seam + 1);
    }
    Check(table.Get(), UINT32_MAX); Check(table.Get(), UINT32_MAX - 1);
    Require(table.storage == original, "immutable neutral lookup/guards");
}

void Matrices()
{
    NeutralLookup table;
    std::vector<std::uint8_t> expected = table.storage;
    std::uint64_t state = UINT64_C(0x62617369732d3635);
    const unsigned matrix_offset = 12288 + 36;
    const unsigned low_offset = 6144 + 900 * 6;
    const std::uint32_t id = 1024 + 900;
    for (unsigned sample = 0; sample < 8; ++sample) {
        for (unsigned r = 0; r < 6; ++r)
            for (unsigned c = 0; c < 6; ++c) {
                const std::uint8_t value = sample == 0 ? 0 : sample == 1 ?
                    static_cast<std::uint8_t>(r == c) : sample == 2 ? 255 :
                    static_cast<std::uint8_t>(Random(state));
                table.Data()[matrix_offset + r * 6 + c] = value;
                expected[8 + matrix_offset + r * 6 + c] = value;
            }
        for (unsigned c = 0; c < 6; ++c)
            for (unsigned value = 0; value < 256; ++value) {
                std::memset(table.Data() + low_offset, 0, 6);
                std::fill(expected.begin() + 8 + low_offset, expected.begin() + 8 + low_offset + 6, 0);
                table.Data()[low_offset + c] = static_cast<std::uint8_t>(value);
                expected[8 + low_offset + c] = static_cast<std::uint8_t>(value);
                Check(table.Get(), id);
            }
        for (unsigned trial = 0; trial < 256; ++trial) {
            for (unsigned c = 0; c < 6; ++c) {
                const std::uint8_t value = static_cast<std::uint8_t>(Random(state));
                table.Data()[low_offset + c] = value;
                expected[8 + low_offset + c] = value;
            }
            Check(table.Get(), id);
        }
        Require(table.storage == expected, "matrix oracle lookup/guard preservation");
    }
}

void Invalid()
{
    NeutralLookup table;
    const auto original = table.storage;
    std::array<std::uint8_t, 22> output;
    output.fill(kGuard);
    const auto expected = output;
    const std::uintptr_t top = std::numeric_limits<std::uintptr_t>::max();
    const N::Lookup invalid[] = {
        {nullptr, N::kLookupBytes}, {table.Data(), 0}, {table.Data(), N::kLookupBytes - 1},
        {table.Data(), N::kLookupBytes + 1}, {table.Data(), std::numeric_limits<std::size_t>::max()},
        {reinterpret_cast<const std::uint8_t*>(top - 100), N::kLookupBytes}};
    for (N::Lookup lookup : invalid) {
        Require(C::Row(lookup, UINT32_MAX, output.data() + 8) == N::Status::InvalidInput &&
                N::Row(lookup, UINT32_MAX, output.data() + 8) == N::Status::InvalidInput,
                "invalid lookup status");
        Require(output == expected, "invalid lookup wrote output");
    }
    std::uint8_t* const destinations[] = {nullptr, reinterpret_cast<std::uint8_t*>(top - 2),
        table.Data() - 5, table.Data(), table.Data() + N::kLookupBytes - 1};
    for (std::uint8_t* destination : destinations)
        Require(C::Row(table.Get(), UINT32_MAX, destination) == N::Status::InvalidInput &&
                N::Row(table.Get(), UINT32_MAX, destination) == N::Status::InvalidInput,
                "invalid output/span/alias status");
    Require(table.storage == original && output == expected, "invalid call changed bytes");
}
} // namespace

int main()
{
    try {
        Field(); Ids(); Matrices(); Invalid();
        Require(checked_rows == 80994, "complete neutral row roster");
        std::cout << "PASS: 65536 GF products, " << checked_rows << " neutral rows and invalid/guard cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
