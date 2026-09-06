#include "Wh2ThueMorseMatvecGfniR0.h"
#include "gf256.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <sys/mman.h>
#include <unistd.h>

namespace {
using Byte = std::uint8_t;
using Matrix = std::array<Byte, 36>;
using Vector = std::array<Byte, 6>;
using ApplyFunction = void (*)(const Byte*, Byte*);
namespace N = wh2_matvec_gfni_r0;
const ApplyFunction kFunctions[] = {N::Apply, N::ApplyScalar};
std::size_t Cases = 0;

void Check(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

// Carryless multiplication and polynomial reduction, never shared GF tables.
Byte Multiply(Byte a, Byte b)
{
    unsigned product = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if (b & (1u << bit)) product ^= static_cast<unsigned>(a) << bit;
    for (int bit = 14; bit >= 8; --bit)
        if (product & (1u << bit)) product ^= 0x14du << (bit - 8);
    return static_cast<Byte>(product);
}

Vector Oracle(const Matrix& matrix, const Vector& vector)
{
    Vector result = {};
    for (unsigned row = 0; row < 6; ++row)
        for (unsigned column = 0; column < 6; ++column)
            result[row] ^= Multiply(matrix[6u * row + column], vector[column]);
    return result;
}

void Separate(const Matrix& matrix, const Vector& vector,
              unsigned matrix_offset = 0, unsigned vector_offset = 0)
{
    Check(matrix_offset < 16 && vector_offset < 16, "alignment fixture bounds");
    const Vector expected = Oracle(matrix, vector);
    for (ApplyFunction apply : kFunctions) {
        alignas(16) std::array<Byte, 96> matrix_storage;
        alignas(16) std::array<Byte, 64> vector_storage;
        matrix_storage.fill(0xa7);
        vector_storage.fill(0x5d);
        const std::size_t mi = 16u + matrix_offset;
        const std::size_t vi = 16u + vector_offset;
        std::copy(matrix.begin(), matrix.end(), matrix_storage.begin() + mi);
        std::copy(vector.begin(), vector.end(), vector_storage.begin() + vi);
        const auto saved_matrix = matrix_storage;
        auto wanted_vector = vector_storage;
        std::copy(expected.begin(), expected.end(), wanted_vector.begin() + vi);
        Check(reinterpret_cast<std::uintptr_t>(matrix_storage.data() + mi) % 16u == matrix_offset &&
              reinterpret_cast<std::uintptr_t>(vector_storage.data() + vi) % 16u == vector_offset,
              "actual independent alignment offsets");
        apply(matrix_storage.data() + mi, vector_storage.data() + vi);
        Check(matrix_storage == saved_matrix, "matrix bytes and guards must remain unchanged");
        Check(vector_storage == wanted_vector, "oracle result and exactly six writable bytes");
    }
    ++Cases;
}

void ExhaustivePairs()
{
    const std::size_t before = Cases;
    Matrix matrix = {};
    Vector vector = {};
    for (unsigned a = 0; a < 256; ++a) {
        matrix[0] = static_cast<Byte>(a);
        for (unsigned b = 0; b < 256; ++b) {
            vector[0] = static_cast<Byte>(b);
            Separate(matrix, vector);
        }
    }
    Check(Cases - before == 65536u, "complete field-pair roster");
}

void BilinearBasis()
{
    const std::size_t before = Cases;
    for (unsigned position = 0; position < 36; ++position)
        for (unsigned coefficient_bit = 0; coefficient_bit < 8; ++coefficient_bit)
            for (unsigned vector_bit = 0; vector_bit < 48; ++vector_bit) {
                Matrix matrix = {};
                Vector vector = {};
                matrix[position] = static_cast<Byte>(1u << coefficient_bit);
                vector[vector_bit / 8u] = static_cast<Byte>(1u << (vector_bit % 8u));
                Separate(matrix, vector);
            }
    Check(Cases - before == 13824u, "complete matrix/vector bilinear bit-basis roster");
}

std::uint64_t Next(std::uint64_t& state)
{
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

void MixedAndSentinels()
{
    const std::size_t before = Cases;
    std::uint64_t state = UINT64_C(0x243f6a8885a308d3);
    for (unsigned trial = 0; trial < 512; ++trial) {
        Matrix matrix;
        Vector vector;
        // Each case consumes 36 row-major matrix bytes, then six vector bytes;
        // each byte is the low byte of the next frozen xorshift64 state.
        for (Byte& value : matrix) value = static_cast<Byte>(Next(state));
        for (Byte& value : vector) value = static_cast<Byte>(Next(state));
        Separate(matrix, vector);
    }
    Check(Cases - before == 512u, "complete mixed-case roster");
    const Vector vector = {{0, 1, 37, 128, 254, 255}};
    Matrix matrix = {};
    Separate(matrix, vector);
    for (unsigned i = 0; i < 6; ++i) matrix[i * 6u + i] = 1;
    Separate(matrix, vector);
    matrix.fill(1);
    Separate(matrix, vector);
    Check(Cases - before == 515u, "three zero/identity/all-ones sentinels");
}

Matrix LiteralMatrix()
{
    Matrix matrix;
    for (unsigned i = 0; i < 36; ++i)
        matrix[i] = static_cast<Byte>(29u * i + 13u * (i / 6u) + 7u);
    return matrix;
}

void AlignmentsAndOverlaps()
{
    const Matrix matrix = LiteralMatrix();
    const Vector vector = {{0, 1, 37, 128, 254, 255}};
    const std::size_t before = Cases;
    for (unsigned matrix_offset = 0; matrix_offset < 16; ++matrix_offset)
        for (unsigned vector_offset = 0; vector_offset < 16; ++vector_offset)
            Separate(matrix, vector, matrix_offset, vector_offset);
    Check(Cases - before == 256u, "all 256 independent offset pairs");

    for (int relative = -5; relative <= 35; ++relative) {
        std::array<Byte, 128> original;
        for (unsigned i = 0; i < original.size(); ++i)
            original[i] = static_cast<Byte>(43u * i + 19u);
        const std::size_t mi = 32u;
        const std::size_t vi = static_cast<std::size_t>(32 + relative);
        Matrix matrix_before;
        Vector vector_before;
        std::copy_n(original.begin() + mi, 36u, matrix_before.begin());
        std::copy_n(original.begin() + vi, 6u, vector_before.begin());
        auto expected = original;
        const Vector result = Oracle(matrix_before, vector_before);
        std::copy(result.begin(), result.end(), expected.begin() + vi);
        for (ApplyFunction apply : kFunctions) {
            auto actual = original;
            apply(actual.data() + mi, actual.data() + vi);
            Check(actual == expected, "overlap must use original inputs and change only vector bytes");
        }
        ++Cases;
    }
    Check(Cases - before == 297u, "all 41 overlapping matrix/vector placements");
}

void ExactHeapSpans()
{
    const Matrix matrix = LiteralMatrix();
    const Vector vector = {{255, 128, 37, 1, 0, 254}};
    const Vector expected = Oracle(matrix, vector);
    // No readable suffix/prefix canary conceals an overread in an ASAN build.
    std::unique_ptr<Byte[]> matrix_heap(new Byte[36]);
    std::unique_ptr<Byte[]> vector_heap(new Byte[6]);
    for (ApplyFunction apply : kFunctions) {
        std::copy(matrix.begin(), matrix.end(), matrix_heap.get());
        std::copy(vector.begin(), vector.end(), vector_heap.get());
        apply(matrix_heap.get(), vector_heap.get());
        Check(std::equal(matrix.begin(), matrix.end(), matrix_heap.get()), "exact heap matrix immutable");
        Check(std::equal(expected.begin(), expected.end(), vector_heap.get()), "exact heap vector oracle");
    }
    ++Cases;
}

class GuardedEnd {
public:
    GuardedEnd()
    {
        const long page = sysconf(_SC_PAGESIZE);
        Check(page > 36 && static_cast<std::uint64_t>(page) <=
              std::numeric_limits<std::size_t>::max() / 2u, "guard-page size");
        page_ = static_cast<std::size_t>(page);
        mapping_ = mmap(nullptr, 2u * page_, PROT_READ | PROT_WRITE,
                        MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        Check(mapping_ != MAP_FAILED, "guard-page mmap");
        if (mprotect(Data() + page_, page_, PROT_NONE) != 0) {
            munmap(mapping_, 2u * page_);
            mapping_ = MAP_FAILED;
            throw std::runtime_error("guard-page mprotect");
        }
    }
    ~GuardedEnd() { if (mapping_ != MAP_FAILED) munmap(mapping_, 2u * page_); }
    GuardedEnd(const GuardedEnd&) = delete;
    GuardedEnd& operator=(const GuardedEnd&) = delete;
    Byte* Data() { return static_cast<Byte*>(mapping_); }
    Byte* EndSpan(std::size_t bytes) { return Data() + page_ - bytes; }
    std::size_t Page() const { return page_; }
private:
    void* mapping_ = MAP_FAILED;
    std::size_t page_ = 0;
};

void GuardPages()
{
    const Matrix matrix = LiteralMatrix();
    const Vector vector = {{1, 255, 0, 128, 254, 37}};
    const Vector expected = Oracle(matrix, vector);
    GuardedEnd matrix_pages, vector_pages;
    for (ApplyFunction apply : kFunctions) {
        std::fill_n(matrix_pages.Data(), matrix_pages.Page(), Byte{0xa7});
        std::fill_n(vector_pages.Data(), vector_pages.Page(), Byte{0x5d});
        std::copy(matrix.begin(), matrix.end(), matrix_pages.EndSpan(36));
        std::copy(vector.begin(), vector.end(), vector_pages.EndSpan(6));
        apply(matrix_pages.EndSpan(36), vector_pages.EndSpan(6));
        Check(std::equal(matrix.begin(), matrix.end(), matrix_pages.EndSpan(36)), "guard-page matrix immutable");
        Check(std::equal(expected.begin(), expected.end(), vector_pages.EndSpan(6)), "guard-page vector oracle");
        for (std::size_t i = 0; i < matrix_pages.Page() - 36u; ++i)
            Check(matrix_pages.Data()[i] == 0xa7, "guard-page matrix prefix preserved");
        for (std::size_t i = 0; i < vector_pages.Page() - 6u; ++i)
            Check(vector_pages.Data()[i] == 0x5d, "guard-page vector prefix preserved");
    }
    ++Cases;
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2 || (std::strcmp(argv[1], "--require-gfni") != 0 &&
                      std::strcmp(argv[1], "--expect-scalar") != 0)) {
        std::cerr << "usage: --require-gfni | --expect-scalar\n";
        return 2;
    }
    try {
        Check(gf256_init() == 0 && GF256Ctx.Polynomial == 0x14d, "shared ordinary GF256 initialization");
        const bool require_gfni = std::strcmp(argv[1], "--require-gfni") == 0;
        Check(N::Available() == require_gfni, "requested native/scalar dispatch availability");
        ExhaustivePairs();
        BilinearBasis();
        MixedAndSentinels();
        AlignmentsAndOverlaps();
        ExactHeapSpans();
        GuardPages();
        Check(Cases == 80174u, "complete frozen neutral roster");
        std::cout << "PASS " << Cases << " neutral cases, auto+scalar oracle, overlaps and exact ends; "
                  << (require_gfni ? "GFNI available" : "scalar fallback") << "; no timing or selected data\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
