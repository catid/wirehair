#include "Wh2ThueMorseTinyGfniR0.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>
#include <stdexcept>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>

namespace {
namespace N = wh2_thue_tiny_gfni_r0;
using Byte = std::uint8_t;
using Data = std::array<std::array<Byte, 2>, 6>;
using Scales = std::array<Byte, 6>;
using Pointers = std::array<const void*, 6>;
using Function = void (*)(void*, const void* const*, const Byte*);

void Dispatch(void* output, const void* const* sources, const Byte* scales)
{
    N::Payload(output, sources, scales, 6, 2);
}
const Function kFunctions[] = {Dispatch, N::Scalar};

void Check(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

// Independent carryless polynomial product/reduction, never shared GF tables.
Byte Multiply(Byte a, Byte b, unsigned polynomial = 0x14d)
{
    unsigned product = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if (b & (1u << bit)) product ^= static_cast<unsigned>(a) << bit;
    for (int bit = 14; bit >= 8; --bit)
        if (product & (1u << bit)) product ^= polynomial << (bit - 8);
    return static_cast<Byte>(product);
}

std::uint64_t AffineMatrix(const std::array<Byte, 256>& mapping)
{
    std::uint64_t matrix = 0;
    for (unsigned output = 0; output < 8; ++output) {
        unsigned row = 0;
        for (unsigned input = 0; input < 8; ++input)
            row |= ((mapping[1u << input] >> output) & 1u) << input;
        matrix |= static_cast<std::uint64_t>(row) << (8u * (7u - output));
    }
    return matrix;
}

void Basis()
{
    unsigned first_root = 256, roots = 0;
    for (unsigned root = 0; root < 256; ++root) {
        Byte value = 1;
        for (int bit = 7; bit >= 0; --bit)
            value = Multiply(value, static_cast<Byte>(root), 0x11b) ^ ((0x14d >> bit) & 1);
        if (!value) { first_root = std::min(first_root, root); ++roots; }
    }
    Check(first_root == 0x46 && roots == 8, "derive smallest hardware-field root");
    std::array<Byte, 8> powers = {{1}};
    for (unsigned i = 1; i < 8; ++i)
        powers[i] = Multiply(powers[i - 1], static_cast<Byte>(first_root), 0x11b);
    std::array<Byte, 256> forward = {}, reverse = {};
    std::array<bool, 256> seen = {};
    for (unsigned x = 0; x < 256; ++x) {
        for (unsigned bit = 0; bit < 8; ++bit)
            if (x & (1u << bit)) forward[x] ^= powers[bit];
        Check(!seen[forward[x]], "basis map is bijective");
        seen[forward[x]] = true;
        reverse[forward[x]] = static_cast<Byte>(x);
    }
    Check(AffineMatrix(forward) == N::kToHardware && AffineMatrix(reverse) == N::kFromHardware,
          "independently derived GFNI bit ordering and both constants");
    for (unsigned a = 0; a < 256; ++a) for (unsigned b = 0; b < 256; ++b)
        Check(reverse[Multiply(forward[a], forward[b], 0x11b)] == Multiply(a, b),
              "all 65536 products preserve original field representation");
}

std::array<Byte, 2> Oracle(const void* const* sources, const Byte* scales)
{
    std::array<Byte, 2> result = {};
    for (unsigned i = 0; i < 6; ++i) for (unsigned b = 0; b < 2; ++b)
        result[b] ^= Multiply(static_cast<const Byte*>(sources[i])[b], scales[i]);
    return result;
}

void Case(const Data& data, const Scales& scales)
{
    Pointers sources;
    for (unsigned i = 0; i < 6; ++i) sources[i] = data[i].data();
    const auto wanted = Oracle(sources.data(), scales.data());
    const auto saved_data = data;
    const auto saved_scales = scales;
    for (Function function : kFunctions) {
        std::array<Byte, 4> output = {{0xa7, 0x35, 0x8c, 0xa7}};
        function(output.data() + 1, sources.data(), scales.data());
        Check(output[1] == wanted[0] && output[2] == wanted[1] &&
              output.front() == 0xa7 && output.back() == 0xa7, "exact two-byte result/guards");
        Check(data == saved_data && scales == saved_scales, "immutable sources and coefficients");
    }
}

void Arithmetic()
{
    Data data = {};
    Scales scales = {};
    for (unsigned position = 0; position < 6; ++position) {
        for (unsigned scale = 0; scale < 256; ++scale) {
            scales[position] = static_cast<Byte>(scale);
            for (unsigned value = 0; value < 256; ++value) {
                data[position] = {{static_cast<Byte>(value), static_cast<Byte>(255 - value)}};
                Case(data, scales);
            }
        }
        scales[position] = 0;
    }
    // Every coefficient bit against every data bit, including cross-source
    // zeros. This checks lane routing and the XOR fold as a bilinear map.
    for (unsigned scale_bit = 0; scale_bit < 48; ++scale_bit)
        for (unsigned data_bit = 0; data_bit < 96; ++data_bit) {
            data = {}; scales = {};
            scales[scale_bit / 8] = static_cast<Byte>(1u << (scale_bit % 8));
            data[data_bit / 16][(data_bit / 8) % 2] = static_cast<Byte>(1u << (data_bit % 8));
            Case(data, scales);
        }
    std::uint64_t state = UINT64_C(0x243f6a8885a308d3);
    const auto next = [&state] {
        state ^= state << 13; state ^= state >> 7; state ^= state << 17;
        return static_cast<Byte>(state);
    };
    for (unsigned trial = 0; trial < 512; ++trial) {
        for (auto& source : data) for (Byte& value : source) value = next();
        for (Byte& value : scales) value = next();
        if (trial < 2) scales.fill(static_cast<Byte>(trial));
        Case(data, scales);
    }
}

void AlignmentsAndAliases()
{
    for (unsigned si = 0; si < 16; ++si) for (unsigned ci = 0; ci < 16; ++ci)
        for (unsigned oi = 0; oi < 16; ++oi) {
            alignas(16) std::array<std::array<Byte, 48>, 6> data;
            alignas(16) std::array<Byte, 48> scales, output;
            Pointers sources;
            scales.fill(0xa7);
            for (unsigned i = 0; i < 6; ++i) {
                data[i].fill(0xa7);
                data[i][16 + si] = static_cast<Byte>(37*i + ci);
                data[i][17 + si] = static_cast<Byte>(91*i + oi);
                scales[16 + ci + i] = static_cast<Byte>(71*i + si);
                sources[i] = data[i].data() + 16 + si;
                Check(reinterpret_cast<std::uintptr_t>(sources[i]) % 16 == si,
                      "actual source alignment");
            }
            Check(reinterpret_cast<std::uintptr_t>(scales.data() + 16 + ci) % 16 == ci &&
                  reinterpret_cast<std::uintptr_t>(output.data() + 16 + oi) % 16 == oi,
                  "actual independent coefficient/output alignment");
            const auto saved_data = data;
            const auto saved_scales = scales;
            const auto wanted = Oracle(sources.data(), scales.data() + 16 + ci);
            for (Function function : kFunctions) {
                output.fill(0xa7);
                auto expected = output;
                expected[16 + oi] = wanted[0]; expected[17 + oi] = wanted[1];
                function(output.data() + 16 + oi, sources.data(), scales.data() + 16 + ci);
                Check(output == expected && data == saved_data && scales == saved_scales,
                      "4096 independent byte alignments and all guards");
            }
        }
    std::array<Byte, 8> data = {{0, 1, 17, 255, 91, 173, 31, 7}};
    const Scales scales = {{0, 1, 7, 29, 191, 255}};
    for (unsigned pattern = 0; pattern < 2; ++pattern) {
        Pointers sources;
        for (unsigned i = 0; i < 6; ++i) sources[i] = data.data() + (pattern ? i : 0);
        const auto expected = Oracle(sources.data(), scales.data());
        const auto saved = data;
        for (Function function : kFunctions) {
            std::array<Byte, 2> output = {};
            function(output.data(), sources.data(), scales.data());
            Check(output == expected && data == saved, "duplicate and partially overlapping sources");
        }
    }
}

void Fallback()
{
    const unsigned widths[] = {1,2,3,15,16,17,31,32,33,63,64,65,1280};
    const unsigned counts[] = {1,5,6,7,16,17};
    for (unsigned count : counts) for (unsigned width : widths) for (unsigned trial = 0; trial < 16; ++trial) {
        std::vector<std::vector<Byte>> data(count, std::vector<Byte>(width + 34, 0xa7));
        std::vector<const void*> sources(count);
        std::vector<Byte> scales(count);
        for (unsigned i = 0; i < count; ++i) {
            for (unsigned b = 0; b < width; ++b)
                data[i][17 + b] = static_cast<Byte>(i*71 + b*29 + trial*13 + (b >> 2));
            sources[i] = data[i].data() + 17;
            scales[i] = trial < 2 ? static_cast<Byte>(trial) :
                (trial == 2 ? static_cast<Byte>(i % 3) : static_cast<Byte>(i*43 + trial*17));
        }
        if (trial == 15) std::fill(sources.begin(), sources.end(), sources.front());
        const auto saved_data = data;
        const auto saved_scales = scales;
        std::vector<Byte> output(width + 34, 0xa7), baseline = output, expected = output;
        for (unsigned b = 0; b < width; ++b) {
            Byte value = 0;
            for (unsigned i = 0; i < count; ++i)
                value ^= Multiply(static_cast<const Byte*>(sources[i])[b], scales[i]);
            expected[17 + b] = value;
        }
        N::Payload(output.data() + 17, sources.data(), scales.data(), count, width);
        gf256_mulset_multi_mem(baseline.data() + 17, sources.data(), scales.data(), count, width);
        Check(output == expected && baseline == expected && data == saved_data && scales == saved_scales,
              "1248 fallback/mixed cases against generic and polynomial oracle");
    }
    Byte output[2] = {0x35, 0x8c};
    N::Payload(output, nullptr, nullptr, -1, -1);
    N::Payload(output, nullptr, nullptr, -1, 0);
    N::Payload(output, nullptr, nullptr, -1, 2);
    N::Payload(output, nullptr, nullptr, 0, -1);
    N::Payload(output, nullptr, nullptr, 0, 0);
    N::Payload(output, nullptr, nullptr, 0, 2);
    N::Payload(output, nullptr, nullptr, 6, -1);
    N::Payload(output, nullptr, nullptr, 6, 0);
    Check(output[0] == 0x35 && output[1] == 0x8c, "eight nonpositive no-ops");
}

class Guarded {
    Byte* base_ = nullptr;
    std::size_t page_ = 0;
public:
    Guarded()
    {
        const long page = sysconf(_SC_PAGESIZE);
        Check(page >= 128, "guard page size");
        page_ = static_cast<std::size_t>(page);
        void* address = mmap(nullptr, 3*page_, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        Check(address != MAP_FAILED, "guard mapping");
        base_ = static_cast<Byte*>(address);
        if (mprotect(base_ + page_, page_, PROT_READ | PROT_WRITE) != 0) {
            munmap(base_, 3*page_); base_ = nullptr;
            throw std::runtime_error("guard data protection");
        }
    }
    ~Guarded() { if (base_) munmap(base_, 3*page_); }
    Guarded(const Guarded&) = delete;
    Guarded& operator=(const Guarded&) = delete;
    Byte* Edge(std::size_t bytes, bool back) { return base_ + page_ + (back ? page_ - bytes : 0); }
    void ReadOnly() { Check(mprotect(base_ + page_, page_, PROT_READ) == 0, "immutable page protection"); }
};

void ExactSpans()
{
    const Scales scales = {{0, 1, 37, 89, 171, 255}};
    for (unsigned offset = 0; offset < 2; ++offset) {
        std::array<std::unique_ptr<Byte[]>, 6> data;
        std::unique_ptr<Byte[]> coefficients(new Byte[6 + offset]), output(new Byte[2 + offset]);
        Pointers sources;
        for (unsigned i = 0; i < 6; ++i) {
            data[i].reset(new Byte[2 + offset]);
            std::fill_n(data[i].get(), 2 + offset, 0xa7);
            data[i][offset] = static_cast<Byte>(i*37 + 13);
            data[i][offset + 1] = static_cast<Byte>(i*71 + 3);
            sources[i] = data[i].get() + offset;
        }
        std::fill_n(coefficients.get(), 6 + offset, 0xa7);
        std::copy(scales.begin(), scales.end(), coefficients.get() + offset);
        const auto expected = Oracle(sources.data(), scales.data());
        for (Function function : kFunctions) {
            std::fill_n(output.get(), 2 + offset, 0xa7);
            function(output.get() + offset, sources.data(), coefficients.get() + offset);
            Check(output[offset] == expected[0] && output[offset + 1] == expected[1] &&
                  (!offset || output[0] == 0xa7), "exact heap spans and unaligned redzones");
        }
    }
    for (bool back : {false, true}) {
        std::array<Guarded, 9> regions; // Six sources, coefficients, output, pointer array.
        Byte* coefficients = regions[6].Edge(6, back);
        Byte* output = regions[7].Edge(2, back);
        using Pointer = const void*;
        auto* sources = reinterpret_cast<Pointer*>(regions[8].Edge(6 * sizeof(Pointer), back));
        for (unsigned i = 0; i < 6; ++i) {
            Byte* data = regions[i].Edge(2, back);
            data[0] = static_cast<Byte>(i*37 + 13); data[1] = static_cast<Byte>(i*71 + 3);
            new (sources + i) Pointer(data);
            coefficients[i] = scales[i];
            regions[i].ReadOnly();
        }
        regions[6].ReadOnly(); regions[8].ReadOnly();
        const auto expected = Oracle(sources, coefficients);
        for (Function function : kFunctions) {
            function(output, sources, coefficients);
            Check(output[0] == expected[0] && output[1] == expected[1], "protected exact front/back spans");
        }
    }
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
        Check(gf256_init() == 0 && GF256Ctx.Polynomial == 0x14d, "shared ordinary GF256 runtime");
        Check(N::Available() == (std::strcmp(argv[1], "--require-gfni") == 0), "requested GFNI/scalar path");
        Basis(); Arithmetic(); AlignmentsAndAliases(); Fallback(); ExactSpans();
        std::cout << "PASS 65536 basis products, 393216 field/position cases, 4608 bilinear bases, "
                     "512 mixed, 4096 alignments, 1248 fallback, aliases/exact spans/no-ops; no timing\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
