#include "Wh2ThueMorseNativeCodec.h"
#include "gf256.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <vector>

// Neutral executable only: no generated data header or actual candidate pair.
namespace {
bool track_allocations = false;
std::size_t allocation_count = 0;
std::size_t fail_allocation = std::numeric_limits<std::size_t>::max();
std::size_t allocation_sizes[8] = {};
void* allocation_pointers[8] = {};
}

void* operator new(std::size_t bytes)
{
    std::size_t index = 8;
    if (track_allocations) {
        index = allocation_count++;
        if (index < 8) allocation_sizes[index] = bytes;
        if (index == fail_allocation) throw std::bad_alloc();
    }
    void* pointer = std::malloc(bytes ? bytes : 1);
    if (!pointer) throw std::bad_alloc();
    if (index < 8) allocation_pointers[index] = pointer;
    return pointer;
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept
{
    try { return ::operator new(bytes); } catch (const std::bad_alloc&) { return nullptr; }
}
void* operator new[](std::size_t bytes, const std::nothrow_t&) noexcept
{
    try { return ::operator new[](bytes); } catch (const std::bad_alloc&) { return nullptr; }
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept { std::free(pointer); }
void operator delete[](void* pointer, const std::nothrow_t&) noexcept { std::free(pointer); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { std::free(pointer); }
#endif

namespace {
using namespace wh2_thue_native;
typedef std::uint8_t Byte;
typedef std::array<Byte, 36> Matrix;
typedef std::array<Byte, 6> Vector;
typedef std::array<std::array<Matrix, 32>, 2> Blocks;

void Check(bool condition, const char* description)
{
    if (!condition) {
        std::cerr << "FAIL: " << description << '\n';
        std::exit(1);
    }
}

// Independent carryless product followed by polynomial long division.
Byte Multiply(Byte a, Byte b)
{
    unsigned product = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if (b & (1u << bit)) product ^= static_cast<unsigned>(a) << bit;
    for (int bit = 14; bit >= 8; --bit)
        if (product & (1u << bit)) product ^= 0x14du << (bit - 8);
    return static_cast<Byte>(product);
}

Matrix Identity()
{
    Matrix result = {};
    for (unsigned i = 0; i < 6; ++i) result[i * 6 + i] = 1;
    return result;
}

Matrix Product(const Matrix& a, const Matrix& b)
{
    Matrix result = {};
    for (unsigned r = 0; r < 6; ++r)
        for (unsigned c = 0; c < 6; ++c)
            for (unsigned k = 0; k < 6; ++k)
                result[r * 6 + c] ^= Multiply(a[r * 6 + k], b[k * 6 + c]);
    return result;
}

Vector ApplyReference(const Matrix& matrix, const Vector& vector)
{
    Vector result = {};
    for (unsigned r = 0; r < 6; ++r)
        for (unsigned c = 0; c < 6; ++c)
            result[r] ^= Multiply(matrix[r * 6 + c], vector[c]);
    return result;
}

unsigned ParityReference(std::uint32_t value)
{
    unsigned parity = 0;
    while (value) { parity ^= value & 1; value >>= 1; }
    return parity;
}

struct Neutral {
    Blocks blocks;
    std::vector<Byte> table;

    Neutral()
    {
        const Byte feedback[2][6] = {{1, 1, 0, 1, 0, 1}, {3, 1, 0, 1, 0, 1}};
        for (unsigned phase = 0; phase < 2; ++phase) {
            blocks[phase][0].fill(0);
            for (unsigned i = 0; i < 5; ++i) blocks[phase][0][(i + 1) * 6 + i] = 1;
            for (unsigned i = 0; i < 6; ++i) blocks[phase][0][i * 6 + 5] = feedback[phase][i];
        }
        for (unsigned level = 1; level < 32; ++level) {
            blocks[0][level] = Product(blocks[0][level - 1], blocks[1][level - 1]);
            blocks[1][level] = Product(blocks[1][level - 1], blocks[0][level - 1]);
        }
        table.reserve(kLookupBytes);
        Append(0, 10, 0, true); Append(0, 10, 1, true);
        Append(10, 7, 0, false); Append(10, 7, 1, false);
        Append(17, 7, 0, false); Append(17, 7, 1, false);
        Append(24, 8, 0, false);
        Check(table.size() == kLookupBytes, "neutral packed size");
    }

    void Append(unsigned start, unsigned width, unsigned phase, bool vectors)
    {
        Matrix product = Identity();
        for (unsigned value = 0; value < (1u << width); ++value) {
            if (vectors) {
                for (unsigned row = 0; row < 6; ++row) table.push_back(product[row * 6]);
            } else table.insert(table.end(), product.begin(), product.end());
            product = Product(product, blocks[phase ^ ParityReference(value)][start]);
        }
    }

    Lookup View() const { return Lookup{table.data(), table.size()}; }

    // Independent 32-bit decomposition instead of the core's grouped lookup.
    Vector Reference(std::uint32_t id) const
    {
        Vector value = {{1, 0, 0, 0, 0, 0}};
        for (unsigned bit = 0; bit < 32; ++bit)
            if (id & (std::uint32_t(1) << bit)) {
                const unsigned phase = bit == 31 ? 0 : ParityReference(id >> (bit + 1));
                value = ApplyReference(blocks[phase][bit], value);
            }
        return value;
    }
};

std::vector<Byte> Message(std::size_t bytes)
{
    std::vector<Byte> result(bytes);
    for (std::size_t i = 0; i < bytes; ++i)
        result[i] = static_cast<Byte>((i * 29 + (i >> 3) * 71 + 13) & 255);
    return result;
}

std::vector<Byte> PacketReference(const Neutral& neutral, const std::vector<Byte>& message,
                                std::uint32_t block, std::uint32_t id)
{
    const std::size_t tail = message.size() - static_cast<std::size_t>(block) * 5;
    std::vector<Byte> output(id == 5 ? tail : block, 0);
    const Vector row = neutral.Reference(id);
    for (std::size_t byte = 0; byte < output.size(); ++byte)
        for (unsigned source = 0; source < 6; ++source) {
            const std::size_t offset = static_cast<std::size_t>(source) * block + byte;
            if (offset < message.size()) output[byte] ^= Multiply(row[source], message[offset]);
        }
    return output;
}

void MapperTests(const Neutral& neutral)
{
    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b)
            Check(gf256_mul(static_cast<Byte>(a), static_cast<Byte>(b)) ==
                  Multiply(static_cast<Byte>(a), static_cast<Byte>(b)), "shared GF field");
    std::vector<std::uint32_t> ids = {0, 1, 5, 6, 7, 1023, 1024, 1025, 0x12345678u,
                                    0xaaaaaaaau, 0x55555555u, 0xffffffffu};
    for (unsigned bit = 1; bit < 32; ++bit) {
        const std::uint32_t seam = std::uint32_t(1) << bit;
        ids.push_back(seam - 1); ids.push_back(seam); ids.push_back(seam + 1);
    }
    for (std::uint32_t id : ids) {
        std::array<Byte, 8> output = {{0xa7, 0, 0, 0, 0, 0, 0, 0xa7}};
        Check(Row(neutral.View(), id, output.data() + 1) == Status::Success, "mapped status");
        const Vector expected = neutral.Reference(id);
        Check(std::equal(expected.begin(), expected.end(), output.begin() + 1), "mapped byte oracle");
        Check(output.front() == 0xa7 && output.back() == 0xa7, "row canaries");
    }
    Matrix sequential = Identity();
    for (std::uint32_t id = 0; id < 40; ++id) {
        Vector row = {};
        Check(Row(neutral.View(), id, row.data()) == Status::Success, "sequential map status");
        for (unsigned r = 0; r < 6; ++r) Check(row[r] == sequential[r * 6], "right product orientation");
        sequential = Product(sequential, neutral.blocks[ParityReference(id)][0]);
    }
}

void PayloadCase(const Neutral& neutral, std::uint32_t block, std::uint32_t tail)
{
    const std::vector<Byte> message = Message(static_cast<std::size_t>(block) * 5 + tail);
    std::unique_ptr<Encoder> encoder;
    Check(Encoder::Create(neutral.View(), message.data(), message.size(), block, encoder) ==
          Status::Success, "encoder create");
    const std::uint32_t ids[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1023, 1024,
        131071, 131072, 16777215, 16777216, 2147483647u, 2147483648u, 4294967295u};
    for (std::uint32_t id : ids) {
        const std::vector<Byte> expected = PacketReference(neutral, message, block, id);
        std::vector<Byte> output(expected.size() + 34, 0xa7);
        const Result result = encoder->Encode(id, output.data() + 17, expected.size());
        Check(result.status == Status::Success && result.bytes_written == expected.size() &&
              result.bytes_required == expected.size(), "packet exact length");
        Check(std::equal(expected.begin(), expected.end(), output.begin() + 17), "packet independent oracle");
        Check(std::count(output.begin(), output.begin() + 17, 0xa7) == 17 &&
              std::count(output.end() - 17, output.end(), 0xa7) == 17, "packet canaries");
        std::fill(output.begin(), output.end(), 0xa7);
        const Result short_result = encoder->Encode(id, output.data() + 17, expected.size() - 1);
        Check(short_result.status == Status::BufferTooSmall && short_result.bytes_written == 0 &&
              std::count(output.begin(), output.end(), 0xa7) == static_cast<std::ptrdiff_t>(output.size()),
              "short packet untouched");
    }
    const std::uint32_t orders[3][6] = {{0, 1, 2, 3, 4, 5}, {5, 3, 1, 4, 2, 0}, {6, 7, 8, 9, 10, 11}};
    for (const auto& order : orders) {
        std::unique_ptr<Decoder> decoder;
        Check(Decoder::Create(neutral.View(), message.size(), block, decoder) == Status::Success,
              "decoder create");
        std::vector<Byte> output(message.size() + 34, 0xa7);
        Check(decoder->Recover(output.data() + 17, message.size()).status == Status::NeedMore,
              "empty recover needs more");
        for (unsigned i = 0; i < 6; ++i) {
            const std::vector<Byte> packet = PacketReference(neutral, message, block, order[i]);
            const Status expected = i == 5 ? Status::Success : Status::NeedMore;
            Check(decoder->Feed(order[i], packet.data(), packet.size()).status == expected &&
                  decoder->Rank() == i + 1, "incremental rank");
            Check(decoder->Feed(order[i], packet.data(), packet.size()).status == expected &&
                  decoder->Rank() == i + 1, "duplicate idempotence before recover");
            std::vector<Byte> conflict = packet;
            conflict.back() ^= 1;
            Check(decoder->Feed(order[i], conflict.data(), conflict.size()).status == Status::Conflict &&
                  decoder->Rank() == i + 1, "conflict preserves rank before recover");
        }
        // New dependent equations must also preserve the unsolved echelon,
        // including a distant-ID row. These IDs already have packet oracles
        // above and are absent from each of the six-packet accepted orders.
        const std::uint32_t novel_ids[] = {ids[12], ids[20]};
        for (std::uint32_t id : novel_ids) {
            Check(std::find(order, order + 6, id) == order + 6, "dependent fixture has novel ID");
            const std::vector<Byte> packet = PacketReference(neutral, message, block, id);
            Check(decoder->Feed(id, packet.data(), packet.size()).status == Status::Success &&
                  decoder->Rank() == 6, "novel dependent packet before recover");
            std::vector<Byte> conflict = packet;
            conflict[conflict.size() / 2] ^= 1;
            Check(decoder->Feed(id, conflict.data(), conflict.size()).status == Status::Conflict &&
                  decoder->Rank() == 6, "novel dependent conflict preserves unsolved rank");
            Check(decoder->Feed(id, packet.data(), packet.size()).status == Status::Success &&
                  decoder->Rank() == 6, "valid continuation after novel dependent conflict");
        }
        Check(decoder->Recover(output.data() + 17, message.size() - 1).status == Status::BufferTooSmall,
              "short recover rejected before solve");
        Check(std::count(output.begin(), output.end(), 0xa7) == static_cast<std::ptrdiff_t>(output.size()),
              "failed recover output untouched");
        Check(decoder->Recover(output.data() + 17, message.size()).status == Status::Success &&
              std::equal(message.begin(), message.end(), output.begin() + 17), "exact recovered message");
        for (std::uint32_t id : ids) {
            const std::vector<Byte> packet = PacketReference(neutral, message, block, id);
            Check(decoder->Feed(id, packet.data(), packet.size()).status == Status::Success,
                  "identity basis accepts subsequent packet");
            std::vector<Byte> conflict = packet;
            conflict.front() ^= 0x80;
            Check(decoder->Feed(id, conflict.data(), conflict.size()).status == Status::Conflict,
                  "identity basis rejects conflict");
        }
        Check(decoder->Recover(output.data() + 17, message.size()).status == Status::Success &&
              std::equal(message.begin(), message.end(), output.begin() + 17), "repeat recovery unchanged");
        Check(std::count(output.begin(), output.begin() + 17, 0xa7) == 17 &&
              std::count(output.end() - 17, output.end(), 0xa7) == 17, "recover canaries");
    }
}

void InvalidTests(const Neutral& neutral)
{
    std::vector<Byte> message = Message(11);
    std::unique_ptr<Encoder> encoder;
    std::unique_ptr<Decoder> decoder;
    const std::uint64_t bad_messages[] = {0, 10, 13, std::numeric_limits<std::uint64_t>::max()};
    for (std::uint64_t bytes : bad_messages) {
        Check(Encoder::Create(neutral.View(), message.data(), bytes, 2, encoder) == Status::InvalidInput && !encoder,
              "bad encoder shape");
        Check(Decoder::Create(neutral.View(), bytes, 2, decoder) == Status::InvalidInput && !decoder,
              "bad decoder shape");
    }
    Check(Encoder::Create(neutral.View(), message.data(), 11, 0, encoder) == Status::InvalidInput, "zero block");
    Check(Decoder::Create(neutral.View(), 6ull * 0x80000000u, 0x80000000u, decoder) == Status::InvalidInput,
          "GF signed byte cap");
    Check(Encoder::Create(neutral.View(), nullptr, 11, 2, encoder) == Status::InvalidInput, "null source");
    Check(Encoder::Create(neutral.View(), &encoder, 6, 1, encoder) == Status::InvalidInput && !encoder,
          "create output cannot overwrite borrowed source");
    const void* wrapped = reinterpret_cast<const void*>(std::numeric_limits<std::uintptr_t>::max() - 1);
    Check(Encoder::Create(neutral.View(), wrapped, 11, 2, encoder) == Status::InvalidInput, "wrapped source");
    const Lookup bad_lookup = {neutral.table.data(), kLookupBytes - 1};
    Check(Decoder::Create(bad_lookup, 11, 2, decoder) == Status::InvalidInput, "short lookup");
    std::vector<Byte> corrupt_table = neutral.table;
    corrupt_table[0] = 0;
    const Lookup corrupt = {corrupt_table.data(), corrupt_table.size()};
    Check(Decoder::Create(corrupt, 11, 2, decoder) == Status::InvalidInput, "nonsystematic lookup");
    Vector row = {{0xa7, 0xa7, 0xa7, 0xa7, 0xa7, 0xa7}};
    Check(Row(bad_lookup, 0, row.data()) == Status::InvalidInput && row[0] == 0xa7, "invalid row untouched");
    Check(Row(neutral.View(), 0, const_cast<Byte*>(neutral.table.data())) == Status::InvalidInput,
          "row lookup overlap");
    Check(Encoder::Create(neutral.View(), message.data(), 11, 2, encoder) == Status::Success, "valid after invalid");
    Encoder* original = encoder.get();
    Check(Encoder::Create(neutral.View(), message.data(), 11, 2, encoder) == Status::InvalidInput &&
          encoder.get() == original, "nonempty create preserves handle");
    Check(encoder->Encode(0, message.data(), 2).status == Status::InvalidInput, "systematic source overlap");
    Check(encoder->Encode(6, message.data() + 1, 2).status == Status::InvalidInput, "repair source overlap");
    Check(encoder->Encode(6, encoder.get(), 2).status == Status::InvalidInput, "encoder metadata overlap");
    Check(encoder->Encode(6, const_cast<Byte*>(neutral.table.data()), 2).status == Status::InvalidInput,
          "encoder lookup overlap");
    Byte output[16] = {};
    Check(encoder->Encode(6, output, std::numeric_limits<std::size_t>::max()).status == Status::InvalidInput,
          "wrapped output capacity");
    Check(encoder->Encode(5, nullptr, 0).status == Status::BufferTooSmall, "zero capacity query");
    Check(encoder->Encode(5, nullptr, 1).status == Status::InvalidInput, "null adequate output");
    Check(Decoder::Create(neutral.View(), 11, 2, decoder) == Status::Success, "invalid test decoder");
    Check(decoder->Feed(5, output, 2).status == Status::InvalidInput && decoder->Rank() == 0,
          "oversized partial input");
    Check(decoder->Feed(6, output, 1).status == Status::InvalidInput && decoder->Rank() == 0,
          "short repair input");
    Check(decoder->Feed(0, nullptr, 2).status == Status::InvalidInput && decoder->Rank() == 0, "null input");
    Check(decoder->Feed(0, wrapped, 2).status == Status::InvalidInput && decoder->Rank() == 0, "wrapped input");
    Check(decoder->Feed(0, decoder.get(), 2).status == Status::InvalidInput, "decoder metadata input");
    Check(decoder->Recover(decoder.get(), 11).status == Status::InvalidInput, "decoder metadata output");
    Check(decoder->Recover(const_cast<Byte*>(neutral.table.data()), 11).status == Status::InvalidInput,
          "decoder lookup output");
    Check(decoder->Recover(output, std::numeric_limits<std::size_t>::max()).status == Status::InvalidInput,
          "recover capacity overflow");
    Check(message == Message(11), "invalid encodes preserve source");
}

void StartTracking(std::size_t fail = std::numeric_limits<std::size_t>::max())
{
    allocation_count = 0;
    fail_allocation = fail;
    std::fill(allocation_sizes, allocation_sizes + 8, 0);
    std::fill(allocation_pointers, allocation_pointers + 8, nullptr);
    track_allocations = true;
}

void AllocationTests(const Neutral& neutral)
{
    const std::vector<Byte> message = Message(384);
    std::unique_ptr<Encoder> encoder;
    std::unique_ptr<Decoder> decoder;
    StartTracking();
    Status status = Encoder::Create(neutral.View(), message.data(), 384, 64, encoder);
    track_allocations = false;
    Check(status == Status::Success && allocation_count == 1 && allocation_sizes[0] == sizeof(Encoder),
          "full borrowed encoder one handle allocation");
    encoder.reset();
    StartTracking();
    status = Encoder::Create(neutral.View(), message.data(), 321, 64, encoder);
    track_allocations = false;
    Check(status == Status::Success && allocation_count == 2 && allocation_sizes[1] == 64,
          "partial encoder only one padding block");
    Check(encoder->Encode(6, allocation_pointers[1], 64).status == Status::InvalidInput,
          "private padding output overlap");
    encoder.reset();
    StartTracking();
    status = Decoder::Create(neutral.View(), 321, 64, decoder);
    track_allocations = false;
    Check(status == Status::Success && allocation_count == 2 && allocation_sizes[0] == sizeof(Decoder) &&
          allocation_sizes[1] == 7 * 64, "decoder handle and seven-block slab");
    Check(decoder->Feed(0, allocation_pointers[1], 64).status == Status::InvalidInput && decoder->Rank() == 0,
          "owned slab input overlap");
    Check(decoder->Recover(allocation_pointers[1], 321).status == Status::InvalidInput,
          "owned slab output overlap");
    decoder.reset();
    for (std::size_t fail = 0; fail < 2; ++fail) {
        StartTracking(fail);
        status = Encoder::Create(neutral.View(), message.data(), 321, 64, encoder);
        track_allocations = false;
        Check(status == Status::OutOfMemory && !encoder, "encoder allocation failure atomic");
        StartTracking(fail);
        status = Decoder::Create(neutral.View(), 321, 64, decoder);
        track_allocations = false;
        Check(status == Status::OutOfMemory && !decoder, "decoder allocation failure atomic");
    }
    Check(Encoder::Create(neutral.View(), message.data(), 384, 64, encoder) == Status::Success, "allocation encoder");
    Check(Decoder::Create(neutral.View(), 384, 64, decoder) == Status::Success, "allocation decoder");
    std::array<std::array<Byte, 64>, 6> packets = {};
    Byte output[384] = {};
    for (unsigned i = 0; i < 6; ++i)
        Check(encoder->Encode(6 + i, packets[i].data(), 64).status == Status::Success, "allocation packet prep");
    StartTracking(0);
    bool good = true;
    for (unsigned i = 0; i < 6; ++i) {
        good = good && encoder->Encode(6 + i, packets[i].data(), 64).status == Status::Success;
        good = good && decoder->Feed(6 + i, packets[i].data(), 64).status ==
            (i == 5 ? Status::Success : Status::NeedMore);
    }
    good = good && decoder->Recover(output, sizeof(output)).status == Status::Success;
    good = good && decoder->Feed(6, packets[0].data(), 64).status == Status::Success;
    track_allocations = false;
    Check(good && allocation_count == 0 && std::equal(message.begin(), message.end(), output),
          "encode feed recover no allocation");
}

} // namespace

int main()
{
    Check(gf256_init() == 0, "shared GF initialization");
    const Neutral neutral;
    MapperTests(neutral);
    InvalidTests(neutral);
    AllocationTests(neutral);
    const std::uint32_t widths[] = {1, 2, 3, 15, 16, 17, 31, 32, 33, 63, 64, 65,
        127, 128, 129, 255, 256, 257, 1023, 1024, 1025, 1279, 1280, 1281};
    unsigned shapes = 0;
    for (std::uint32_t width : widths) {
        PayloadCase(neutral, width, width); ++shapes;
        if (width > 1) { PayloadCase(neutral, width, 1); ++shapes; }
        if (width > 2) { PayloadCase(neutral, width, width - 1); ++shapes; }
    }
    std::cout << "PASS neutral field, mapper, invalid calls, allocations and " << shapes
              << " payload shapes; no actual candidate data\n";
    return 0;
}
