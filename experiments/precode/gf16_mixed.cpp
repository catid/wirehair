// Self-contained payload prototype for a future WH2 mixed-field profile.
// Keeps 12 completion variables/constraints: 10 GF(256) rows followed by
// two GF(2^16) rows. Blocks must contain an even number of bytes.

#include "../../gf256.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {

const uint8_t Lambda = 32;
const uint16_t Generator = 266;
const uint32_t Period = 244, Rows = 12, ByteRows = 10;
uint16_t Extension[2][Period] = {};
uint8_t ByteTable[Rows][Period] = {};

uint64_t Mix64(uint64_t x)
{
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

uint16_t Mul(uint16_t x, uint16_t y)
{
    const uint8_t a = (uint8_t)x, b = (uint8_t)(x >> 8);
    const uint8_t c = (uint8_t)y, d = (uint8_t)(y >> 8);
    const uint8_t ac = gf256_mul(a, c), bd = gf256_mul(b, d);
    const uint8_t cross = gf256_mul((uint8_t)(a ^ b), (uint8_t)(c ^ d));
    return (uint16_t)(uint8_t)(ac ^ gf256_mul(Lambda, bd)) |
        (uint16_t)(uint8_t)(cross ^ ac) << 8;
}

uint16_t Pow(uint16_t x, uint32_t exponent)
{
    uint16_t result = 1;
    while (exponent) {
        if (exponent & 1) result = Mul(result, x);
        exponent >>= 1;
        if (exponent) x = Mul(x, x);
    }
    return result;
}

uint16_t Inv(uint16_t x) { return x ? Pow(x, 65534) : 0; }

bool InitializeField()
{
    for (unsigned x = 0; x < 256; ++x)
        if ((uint8_t)(gf256_mul((uint8_t)x, (uint8_t)x) ^ x) == Lambda)
            return false; // u^2 + u + Lambda must be irreducible.
    const uint32_t factors[] = {3, 5, 17, 257};
    for (uint32_t factor : factors)
        if (Pow(Generator, 65535 / factor) == 1) return false;
    for (unsigned row = 0; row < 2; ++row) {
        const uint16_t y = Pow(Generator, 1000 + row);
        for (unsigned column = 0; column < Period; ++column) {
            const uint16_t denominator =
                (uint16_t)(Pow(Generator, column) ^ y);
            if (!denominator) return false;
            Extension[row][column] = Inv(denominator);
        }
    }
    for (uint32_t row = 0; row < Rows; ++row)
        for (uint32_t column = 0; column < Period; ++column) {
            const uint8_t x = (uint8_t)(12 + column);
            ByteTable[row][column] =
                gf256_inv((uint8_t)(x ^ (uint8_t)row));
        }
    // Pin the coefficient family so a future polynomial/generator change is
    // visible rather than silently changing a proposed wire profile.
    return Extension[0][0] == 34916u &&
        Extension[1][0] == 2472u &&
        Extension[0][243] == 59155u;
}

bool MulAdd(uint8_t* z, uint16_t scale, const uint8_t* x, uint32_t bytes)
{
    if (!z || !x || !bytes || (bytes & 1)) return false;
    for (uint32_t i = 0; i < bytes; i += 2) {
        const uint16_t value = (uint16_t)x[i] | (uint16_t)x[i + 1] << 8;
        const uint16_t product = Mul(scale, value);
        z[i] ^= (uint8_t)product;
        z[i + 1] ^= (uint8_t)(product >> 8);
    }
    return true;
}

bool Scale(uint8_t* x, uint16_t scale, uint32_t bytes)
{
    if (!x || !bytes || (bytes & 1)) return false;
    for (uint32_t i = 0; i < bytes; i += 2) {
        const uint16_t value = (uint16_t)x[i] | (uint16_t)x[i + 1] << 8;
        const uint16_t product = Mul(scale, value);
        x[i] = (uint8_t)product;
        x[i + 1] = (uint8_t)(product >> 8);
    }
    return true;
}

void Deinterleave(const uint8_t* x, uint8_t* lo, uint8_t* hi, uint32_t bytes)
{
    for (uint32_t i = 0; i < bytes / 2; ++i) {
        lo[i] = x[2 * i];
        hi[i] = x[2 * i + 1];
    }
}

void Interleave(const uint8_t* lo, const uint8_t* hi, uint8_t* x, uint32_t bytes)
{
    for (uint32_t i = 0; i < bytes / 2; ++i) {
        x[2 * i] = lo[i];
        x[2 * i + 1] = hi[i];
    }
}

void MulAddPlanar(
    uint8_t* zlo, uint8_t* zhi, uint16_t scale,
    const uint8_t* xlo, const uint8_t* xhi, uint32_t elements)
{
    const uint8_t a = (uint8_t)scale, b = (uint8_t)(scale >> 8);
    gf256_muladd_mem(zlo, a, xlo, (int)elements);
    gf256_muladd_mem(zlo, gf256_mul(Lambda, b), xhi, (int)elements);
    gf256_muladd_mem(zhi, b, xlo, (int)elements);
    gf256_muladd_mem(zhi, (uint8_t)(a ^ b), xhi, (int)elements);
}

uint8_t ByteCoefficient(uint32_t row, uint32_t column)
{
    return ByteTable[row][column % Period];
}

uint16_t Coefficient(uint32_t row, uint32_t column)
{
    return row < ByteRows ? ByteCoefficient(row, column) :
        Extension[row - ByteRows][column % Period];
}

uint32_t Rank(std::vector<uint16_t> a, uint32_t rows, uint32_t columns)
{
    uint32_t rank = 0;
    for (uint32_t column = 0; column < columns && rank < rows; ++column) {
        uint32_t pivot = rank;
        while (pivot < rows && !a[(size_t)pivot * columns + column]) ++pivot;
        if (pivot == rows) continue;
        for (uint32_t j = column; j < columns; ++j)
            std::swap(a[(size_t)rank * columns + j],
                      a[(size_t)pivot * columns + j]);
        const uint16_t inverse = Inv(a[(size_t)rank * columns + column]);
        for (uint32_t j = column; j < columns; ++j)
            a[(size_t)rank * columns + j] =
                Mul(a[(size_t)rank * columns + j], inverse);
        for (uint32_t r = 0; r < rows; ++r) {
            if (r == rank) continue;
            const uint16_t scale = a[(size_t)r * columns + column];
            if (!scale) continue;
            for (uint32_t j = column; j < columns; ++j)
                a[(size_t)r * columns + j] ^=
                    Mul(scale, a[(size_t)rank * columns + j]);
        }
        ++rank;
    }
    return rank;
}

std::vector<uint16_t> Project(
    uint32_t columns, uint32_t global_columns, uint64_t seed)
{
    std::vector<uint16_t> matrix((size_t)Rows * columns);
    for (uint32_t global = 0; global < global_columns; ++global) {
        const uint32_t mask = global < columns ? 1u << global :
            (uint32_t)Mix64(seed ^ global) & ((1u << columns) - 1);
        for (uint32_t row = 0; row < Rows; ++row) {
            uint32_t bits = mask;
            const uint16_t coefficient = Coefficient(row, global);
            while (bits) {
                uint32_t bit = 0;
                while ((bits & (1u << bit)) == 0u) {
                    ++bit;
                }
                matrix[(size_t)row * columns + bit] ^= coefficient;
                bits &= bits - 1;
            }
        }
    }
    return matrix;
}

bool Solve(
    std::vector<uint16_t> a, uint32_t columns,
    std::vector<uint8_t>& rhs, uint32_t bytes)
{
    if (!columns || columns > Rows || (bytes & 1) ||
        rhs.size() != (size_t)Rows * bytes) return false;
    for (uint32_t column = 0; column < columns; ++column) {
        uint32_t pivot = column;
        while (pivot < Rows && !a[(size_t)pivot * columns + column]) ++pivot;
        if (pivot == Rows) return false;
        if (pivot != column) {
            for (uint32_t j = 0; j < columns; ++j)
                std::swap(a[(size_t)column * columns + j],
                          a[(size_t)pivot * columns + j]);
            for (uint32_t j = 0; j < bytes; ++j)
                std::swap(rhs[(size_t)column * bytes + j],
                          rhs[(size_t)pivot * bytes + j]);
        }
        const uint16_t inverse = Inv(a[(size_t)column * columns + column]);
        for (uint32_t j = 0; j < columns; ++j)
            a[(size_t)column * columns + j] =
                Mul(a[(size_t)column * columns + j], inverse);
        if (!Scale(rhs.data() + (size_t)column * bytes, inverse, bytes))
            return false;
        for (uint32_t r = 0; r < Rows; ++r) {
            if (r == column) continue;
            const uint16_t scale = a[(size_t)r * columns + column];
            if (!scale) continue;
            for (uint32_t j = 0; j < columns; ++j)
                a[(size_t)r * columns + j] ^=
                    Mul(scale, a[(size_t)column * columns + j]);
            if (!MulAdd(rhs.data() + (size_t)r * bytes, scale,
                        rhs.data() + (size_t)column * bytes, bytes))
                return false;
        }
    }
    rhs.resize((size_t)columns * bytes);
    return true;
}

bool RoundTrip(uint32_t columns, uint32_t bytes, uint64_t seed)
{
    std::vector<uint16_t> matrix;
    unsigned attempt = 0;
    do {
        matrix = Project(columns, 1076, seed ^ ((uint64_t)attempt << 32));
    } while (Rank(matrix, Rows, columns) != columns && ++attempt < 256);
    if (attempt == 256) return false;
    std::vector<uint8_t> unknown((size_t)columns * bytes);
    for (size_t i = 0; i < unknown.size(); ++i)
        unknown[i] = (uint8_t)Mix64(seed ^ i);
    std::vector<uint8_t> rhs((size_t)Rows * bytes);
    for (uint32_t row = 0; row < Rows; ++row)
        for (uint32_t column = 0; column < columns; ++column) {
            const uint16_t scale = matrix[(size_t)row * columns + column];
            if (scale <= 255)
                gf256_muladd_mem(rhs.data() + (size_t)row * bytes,
                    (uint8_t)scale,
                    unknown.data() + (size_t)column * bytes, (int)bytes);
            else if (!MulAdd(rhs.data() + (size_t)row * bytes, scale,
                             unknown.data() + (size_t)column * bytes, bytes))
                return false;
        }
    return Solve(matrix, columns, rhs, bytes) && rhs == unknown;
}

bool FieldAndKernelTests()
{
    for (uint32_t x = 1; x < 65536; ++x)
        if (Mul((uint16_t)x, Inv((uint16_t)x)) != 1) return false;
    for (uint32_t trial = 0; trial < 100000; ++trial) {
        const uint16_t a = (uint16_t)Mix64(3 * trial);
        const uint16_t b = (uint16_t)Mix64(3 * trial + 1);
        const uint16_t c = (uint16_t)Mix64(3 * trial + 2);
        if (Mul(a, b) != Mul(b, a) ||
            Mul(Mul(a, b), c) != Mul(a, Mul(b, c)) ||
            Mul(a, (uint16_t)(b ^ c)) != (uint16_t)(Mul(a, b) ^ Mul(a, c)))
            return false;
    }
    const uint32_t lengths[] = {2, 4, 30, 128, 1280};
    const uint16_t scales[] = {0, 1, 0xd3, 0xc753, 0xffff};
    for (uint32_t bytes : lengths) for (uint16_t scale : scales) {
        std::vector<uint8_t> source(bytes + 2, 0xa5);
        std::vector<uint8_t> direct(bytes + 2, 0x5a), planar = direct;
        for (uint32_t i = 0; i < bytes; ++i) {
            source[i + 1] = (uint8_t)Mix64(i + bytes);
            direct[i + 1] = planar[i + 1] = (uint8_t)Mix64(i + scale);
        }
        std::vector<uint8_t> xlo(bytes / 2), xhi(bytes / 2);
        std::vector<uint8_t> zlo(bytes / 2), zhi(bytes / 2);
        Deinterleave(source.data() + 1, xlo.data(), xhi.data(), bytes);
        Deinterleave(planar.data() + 1, zlo.data(), zhi.data(), bytes);
        if (!MulAdd(direct.data() + 1, scale, source.data() + 1, bytes))
            return false;
        MulAddPlanar(zlo.data(), zhi.data(), scale,
                     xlo.data(), xhi.data(), bytes / 2);
        Interleave(zlo.data(), zhi.data(), planar.data() + 1, bytes);
        if (direct != planar || direct.front() != 0x5a ||
            direct.back() != 0x5a) return false;
    }
    uint8_t source[4] = {1, 2, 3, 4}, destination[4] = {9, 8, 7, 6};
    const uint8_t snapshot[4] = {9, 8, 7, 6};
    return !MulAdd(destination, 7, source, 3) &&
        std::memcmp(destination, snapshot, 4) == 0;
}

bool CompletionCornerTests()
{
    // The encoder solves the 12 completion parity columns before emitting a
    // packet.  Certify every possible 244-period consecutive corner and a
    // deterministic sample of non-consecutive 12-column subsets.
    for (uint32_t start = 0; start < Period; ++start) {
        std::vector<uint16_t> matrix((size_t)Rows * Rows);
        for (uint32_t row = 0; row < Rows; ++row)
            for (uint32_t column = 0; column < Rows; ++column)
                matrix[(size_t)row * Rows + column] =
                    Coefficient(row, (start + column) % Period);
        if (Rank(matrix, Rows, Rows) != Rows) return false;
    }
    for (uint32_t trial = 0; trial < 10000; ++trial) {
        std::vector<uint16_t> matrix((size_t)Rows * Rows);
        uint32_t selected[Rows];
        for (uint32_t column = 0; column < Rows; ++column) {
            uint32_t candidate =
                (uint32_t)Mix64(((uint64_t)trial << 32) ^ column) % Period;
            for (;;) {
                bool duplicate = false;
                for (uint32_t previous = 0; previous < column; ++previous)
                    duplicate |= selected[previous] == candidate;
                if (!duplicate) break;
                candidate = (candidate + 1) % Period;
            }
            selected[column] = candidate;
            for (uint32_t row = 0; row < Rows; ++row)
                matrix[(size_t)row * Rows + column] =
                    Coefficient(row, candidate);
        }
        if (Rank(matrix, Rows, Rows) != Rows) return false;
    }
    return true;
}

bool BoundedRhsTest()
{
    const uint32_t global_columns = 1076, bytes = 128;
    std::vector<uint8_t> values((size_t)global_columns * bytes);
    for (size_t i = 0; i < values.size(); ++i) values[i] = (uint8_t)Mix64(i);
    std::vector<uint8_t> direct((size_t)Rows * bytes);
    for (uint32_t column = 0; column < global_columns; ++column)
        for (uint32_t row = 0; row < Rows; ++row) {
            const uint16_t scale = Coefficient(row, column);
            uint8_t* z = direct.data() + (size_t)row * bytes;
            const uint8_t* x = values.data() + (size_t)column * bytes;
            if (row < ByteRows) gf256_muladd_mem(z, (uint8_t)scale, x, bytes);
            else if (!MulAdd(z, scale, x, bytes)) return false;
        }
    std::vector<uint8_t> bucket(bytes), byte_rhs((size_t)ByteRows * bytes);
    std::vector<uint8_t> xlo(bytes / 2), xhi(bytes / 2);
    std::vector<uint8_t> extlo(bytes), exthi(bytes);
    void* destinations[ByteRows];
    uint8_t scales[ByteRows];
    for (uint32_t row = 0; row < ByteRows; ++row)
        destinations[row] = byte_rhs.data() + (size_t)row * bytes;
    for (uint32_t residue = 0; residue < Period; ++residue) {
        std::fill(bucket.begin(), bucket.end(), uint8_t{0});
        for (uint32_t column = residue; column < global_columns; column += Period)
            gf256_add_mem(bucket.data(),
                values.data() + (size_t)column * bytes, bytes);
        for (uint32_t row = 0; row < ByteRows; ++row)
            scales[row] = ByteCoefficient(row, residue);
        gf256_muladd_multi_mem(
            destinations, scales, ByteRows, bucket.data(), bytes);
        Deinterleave(bucket.data(), xlo.data(), xhi.data(), bytes);
        for (uint32_t row = 0; row < 2; ++row)
            MulAddPlanar(extlo.data() + (size_t)row * bytes / 2,
                         exthi.data() + (size_t)row * bytes / 2,
                         Extension[row][residue], xlo.data(), xhi.data(),
                         bytes / 2);
    }
    std::vector<uint8_t> bucketed((size_t)Rows * bytes);
    std::memcpy(bucketed.data(), byte_rhs.data(), byte_rhs.size());
    for (uint32_t row = 0; row < 2; ++row)
        Interleave(extlo.data() + (size_t)row * bytes / 2,
                   exthi.data() + (size_t)row * bytes / 2,
                   bucketed.data() + (size_t)(ByteRows + row) * bytes, bytes);
    return direct == bucketed;
}

void Benchmark(uint32_t bytes, uint32_t rounds)
{
    std::vector<uint8_t> source(bytes, 7);
    std::vector<uint8_t> baseline((size_t)Rows * bytes);
    std::vector<uint8_t> mixed((size_t)ByteRows * bytes);
    std::vector<uint8_t> xlo(bytes / 2), xhi(bytes / 2);
    std::vector<uint8_t> extlo(bytes), exthi(bytes);
    void* baseline_destinations[Rows];
    void* mixed_destinations[ByteRows];
    uint8_t baseline_scales[Rows], mixed_scales[ByteRows];
    for (uint32_t row = 0; row < Rows; ++row) {
        baseline_destinations[row] = baseline.data() + (size_t)row * bytes;
    }
    for (uint32_t row = 0; row < ByteRows; ++row) {
        mixed_destinations[row] = mixed.data() + (size_t)row * bytes;
    }
    const auto start_baseline = std::chrono::steady_clock::now();
    for (uint32_t round = 0; round < rounds; ++round)
        for (uint32_t residue = 0; residue < Period; ++residue) {
            for (uint32_t row = 0; row < Rows; ++row)
                baseline_scales[row] = ByteCoefficient(row, residue);
            gf256_muladd_multi_mem(baseline_destinations, baseline_scales,
                Rows, source.data(), bytes);
        }
    const auto stop_baseline = std::chrono::steady_clock::now();
    const auto start_mixed = std::chrono::steady_clock::now();
    for (uint32_t round = 0; round < rounds; ++round)
        for (uint32_t residue = 0; residue < Period; ++residue) {
            for (uint32_t row = 0; row < ByteRows; ++row)
                mixed_scales[row] = ByteCoefficient(row, residue);
            gf256_muladd_multi_mem(mixed_destinations, mixed_scales,
                ByteRows, source.data(), bytes);
            Deinterleave(source.data(), xlo.data(), xhi.data(), bytes);
            for (uint32_t row = 0; row < 2; ++row)
                MulAddPlanar(extlo.data() + (size_t)row * bytes / 2,
                             exthi.data() + (size_t)row * bytes / 2,
                             Extension[row][residue], xlo.data(), xhi.data(),
                             bytes / 2);
        }
    const auto stop_mixed = std::chrono::steady_clock::now();
    const double base_seconds =
        std::chrono::duration<double>(stop_baseline - start_baseline).count();
    const double mixed_seconds =
        std::chrono::duration<double>(stop_mixed - start_mixed).count();
    std::printf("kernel,bb=%u,rounds=%u,baseline_s=%.6f,mixed_s=%.6f,"
                "ratio=%.3f,checksum=%u\n", bytes, rounds, base_seconds,
                mixed_seconds, mixed_seconds / base_seconds,
                (unsigned)(baseline[0] ^ mixed[0] ^ extlo[0] ^ exthi[0]));
}

} // namespace

int main(int argc, char** argv)
{
    bool test_only = false;
    if (argc == 2 && std::strcmp(argv[1], "--test-only") == 0) {
        test_only = true;
    }
    else if (argc != 1) {
        std::fprintf(stderr, "usage: gf16_mixed [--test-only]\n");
        return 2;
    }
    if (gf256_init() || !InitializeField() ||
        !FieldAndKernelTests() || !CompletionCornerTests()) {
        std::fprintf(stderr, "field/kernel test failed\n");
        return 1;
    }
    const uint32_t quotient_sizes[] = {1, 2, 10, 12};
    const uint32_t block_sizes[] = {2, 16, 1280};
    for (uint32_t columns : quotient_sizes)
        for (uint32_t bytes : block_sizes)
            if (!RoundTrip(columns, bytes,
                    UINT64_C(0x5eed1610) ^ ((uint64_t)columns << 32) ^ bytes)) {
                std::fprintf(stderr, "roundtrip failed q=%u bb=%u\n",
                    columns, bytes);
                return 1;
            }
    const std::vector<uint16_t> overwide = Project(13, 1076, 0x5eed1613);
    std::vector<uint8_t> overwide_rhs((size_t)Rows * 2);
    if (Rank(overwide, Rows, 13) != Rows ||
        Solve(overwide, 13, overwide_rhs, 2) ||
        !BoundedRhsTest()) {
        std::fprintf(stderr, "rank/RHS test failed\n");
        return 1;
    }
    std::printf("PASS field=GF256[u]/(u^2+u+32) generator=%u period=%u "
                "profile=10gf256+2gf65536 even_blocks_only "
                "coeff_golden=%u/%u/%u residues=244 gf256_ops=2440 "
                "gf16_ops=488\n",
                (unsigned)Generator, Period,
                (unsigned)Extension[0][0], (unsigned)Extension[1][0],
                (unsigned)Extension[0][243]);
    if (test_only) {
        return 0;
    }
    Benchmark(1280, 2001);
    Benchmark(102400, 21);
    return 0;
}
