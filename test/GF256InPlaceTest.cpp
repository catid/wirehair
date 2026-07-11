#include "../gf256.h"

#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

static const size_t kAlignment = 64;
static const size_t kGuardBytes = 64;

struct AlignedBuffer
{
    explicit AlignedBuffer(size_t bytes, unsigned offset, uint32_t salt)
        : Storage(bytes + 2 * kGuardBytes + 2 * kAlignment)
        , Data(nullptr)
    {
        for (size_t i = 0; i < Storage.size(); ++i) {
            Storage[i] = static_cast<uint8_t>(
                (i * 131u + salt * 17u + (i >> 3)) & 0xffu);
        }

        const uintptr_t start = reinterpret_cast<uintptr_t>(Storage.data()) +
            kGuardBytes;
        const uintptr_t aligned = (start + kAlignment - 1) &
            ~(static_cast<uintptr_t>(kAlignment) - 1);
        Data = reinterpret_cast<uint8_t*>(aligned) + offset;
    }

    std::vector<uint8_t> Storage;
    uint8_t* Data;
};

enum class Operation
{
    Multiply,
    Divide
};

static const char* OperationName(Operation operation)
{
    return operation == Operation::Multiply ? "multiply" : "divide";
}

static uint8_t Reference(
    Operation operation,
    uint8_t value,
    uint8_t factor)
{
    return operation == Operation::Multiply ?
        gf256_mul(value, factor) : gf256_div(value, factor);
}

static void Apply(
    Operation operation,
    void* destination,
    const void* source,
    uint8_t factor,
    int bytes)
{
    if (operation == Operation::Multiply) {
        gf256_mul_mem(destination, source, factor, bytes);
    } else {
        gf256_div_mem(destination, source, factor, bytes);
    }
}

static bool CheckStorage(
    const std::vector<uint8_t>& actual,
    const std::vector<uint8_t>& expected,
    Operation operation,
    unsigned factor,
    unsigned bytes,
    unsigned source_offset,
    unsigned destination_offset,
    const char* mode)
{
    if (actual == expected) {
        return true;
    }

    size_t mismatch = 0;
    while (mismatch < actual.size() &&
           actual[mismatch] == expected[mismatch])
    {
        ++mismatch;
    }
    std::fprintf(
        stderr,
        "%s %s mismatch: factor=%u bytes=%u source-offset=%u "
        "destination-offset=%u storage-index=%zu actual=%u expected=%u\n",
        mode,
        OperationName(operation),
        factor,
        bytes,
        source_offset,
        destination_offset,
        mismatch,
        mismatch < actual.size() ?
            static_cast<unsigned>(actual[mismatch]) : 0u,
        mismatch < expected.size() ?
            static_cast<unsigned>(expected[mismatch]) : 0u);
    return false;
}

static bool TestInPlaceCase(
    Operation operation,
    unsigned factor,
    unsigned bytes,
    unsigned offset)
{
    AlignedBuffer buffer(bytes, offset, factor + bytes + offset);
    std::vector<uint8_t> expected = buffer.Storage;
    const size_t data_index = static_cast<size_t>(
        buffer.Data - buffer.Storage.data());
    for (unsigned i = 0; i < bytes; ++i) {
        expected[data_index + i] = Reference(
            operation, expected[data_index + i],
            static_cast<uint8_t>(factor));
    }

    Apply(
        operation, buffer.Data, buffer.Data,
        static_cast<uint8_t>(factor), static_cast<int>(bytes));
    return CheckStorage(
        buffer.Storage, expected, operation, factor, bytes,
        offset, offset, "in-place");
}

static bool TestOutOfPlaceCase(
    Operation operation,
    unsigned factor,
    unsigned bytes,
    unsigned source_offset,
    unsigned destination_offset)
{
    AlignedBuffer source(
        bytes, source_offset, factor + bytes + source_offset);
    AlignedBuffer destination(
        bytes, destination_offset,
        factor * 3u + bytes + destination_offset + 1u);
    const std::vector<uint8_t> original_source = source.Storage;
    std::vector<uint8_t> expected_destination = destination.Storage;
    const size_t destination_index = static_cast<size_t>(
        destination.Data - destination.Storage.data());
    for (unsigned i = 0; i < bytes; ++i) {
        expected_destination[destination_index + i] = Reference(
            operation, source.Data[i], static_cast<uint8_t>(factor));
    }

    Apply(
        operation, destination.Data, source.Data,
        static_cast<uint8_t>(factor), static_cast<int>(bytes));
    if (!CheckStorage(
            source.Storage, original_source, operation, factor, bytes,
            source_offset, destination_offset, "out-of-place source"))
    {
        return false;
    }
    return CheckStorage(
        destination.Storage, expected_destination, operation, factor, bytes,
        source_offset, destination_offset, "out-of-place destination");
}

static bool TestOperation(Operation operation)
{
    static const unsigned kLengths[] = {
        0, 1, 2, 3, 4, 5, 7, 8, 9,
        15, 16, 17, 31, 32, 33, 47, 48, 49,
        63, 64, 65, 95, 96, 97, 127, 128, 129,
        255, 256, 257, 511, 512, 513, 1023
    };
    static const unsigned kOffsets[] = {
        0, 1, 2, 3, 7, 15, 16, 31, 32, 63
    };

    const unsigned first_factor =
        operation == Operation::Multiply ? 0u : 1u;
    for (unsigned factor = first_factor; factor < 256; ++factor)
    {
        for (unsigned length_i = 0;
             length_i < sizeof(kLengths) / sizeof(kLengths[0]);
             ++length_i)
        {
            const unsigned bytes = kLengths[length_i];
            for (unsigned offset_i = 0;
                 offset_i < sizeof(kOffsets) / sizeof(kOffsets[0]);
                 ++offset_i)
            {
                const unsigned source_offset = kOffsets[offset_i];
                const unsigned destination_offset =
                    kOffsets[(offset_i * 7u + 3u) %
                        (sizeof(kOffsets) / sizeof(kOffsets[0]))];
                if (!TestInPlaceCase(
                        operation, factor, bytes, source_offset) ||
                    !TestOutOfPlaceCase(
                        operation, factor, bytes, source_offset,
                        destination_offset))
                {
                    return false;
                }
            }
        }
    }
    return true;
}

static bool TestMultiDestination()
{
    static const unsigned kLengths[] = {
        0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33,
        63, 64, 65, 127, 128, 129, 255, 256, 257, 1023
    };
    static const int kCounts[] = { 1, 2, 3, 7, 12, 13, 128 };
    uint8_t scales[128];
    for (unsigned i = 0; i < 128u; ++i) {
        scales[i] = (uint8_t)(i * 73u + 2u);
    }
    scales[1] = 0u;
    scales[2] = 1u;

    for (unsigned length_i = 0;
         length_i < sizeof(kLengths) / sizeof(kLengths[0]);
         ++length_i)
    {
        const unsigned bytes = kLengths[length_i];
        for (unsigned offset = 0; offset < 4u; ++offset)
        {
            std::vector<uint8_t> source(
                bytes + 2u * kGuardBytes + 8u, 0xa5u);
            uint8_t* const source_data =
                source.data() + kGuardBytes + offset;
            for (unsigned i = 0; i < bytes; ++i) {
                source_data[i] = (uint8_t)(i * 131u + bytes + offset);
            }
            const std::vector<uint8_t> source_before = source;

            for (int count : kCounts)
            {
                std::vector<std::vector<uint8_t> > destinations(
                    (size_t)count);
                std::vector<std::vector<uint8_t> > expected(
                    (size_t)count);
                std::vector<void*> pointers((size_t)count);
                std::vector<unsigned> offsets((size_t)count);
                for (int j = 0; j < count; ++j)
                {
                    offsets[(size_t)j] =
                        (offset + (unsigned)j * 7u) & 15u;
                    destinations[(size_t)j].resize(
                        bytes + 2u * kGuardBytes + 16u);
                    for (size_t i = 0;
                         i < destinations[(size_t)j].size();
                         ++i)
                    {
                        destinations[(size_t)j][i] = (uint8_t)(
                            i * 17u + (unsigned)j * 43u + bytes);
                    }
                    expected[(size_t)j] = destinations[(size_t)j];
                    pointers[(size_t)j] =
                        destinations[(size_t)j].data() + kGuardBytes +
                        offsets[(size_t)j];
                    const size_t expected_start = kGuardBytes +
                        offsets[(size_t)j];
                    for (unsigned i = 0; i < bytes; ++i) {
                        expected[(size_t)j][expected_start + i] ^=
                            gf256_mul(source_data[i], scales[j]);
                    }
                }

                gf256_muladd_multi_mem(
                    pointers.data(),
                    scales,
                    count,
                    source_data,
                    (int)bytes);
                if (source != source_before)
                {
                    std::fprintf(stderr,
                        "multi-destination source changed bytes=%u "
                        "offset=%u count=%d\n",
                        bytes, offset, count);
                    return false;
                }
                for (int j = 0; j < count; ++j)
                {
                    if (destinations[(size_t)j] != expected[(size_t)j])
                    {
                        std::fprintf(stderr,
                            "multi-destination mismatch bytes=%u offset=%u "
                            "count=%d destination=%d scale=%u\n",
                            bytes, offset, count, j,
                            (unsigned)scales[j]);
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

static void RunBenchmark()
{
    static const size_t kBytes = 4u * 1024u * 1024u;
    static const unsigned kIterations = 2048;
    AlignedBuffer source(kBytes, 0, 0x51u);
    AlignedBuffer destination(kBytes, 0, 0xa7u);

    const auto start = std::chrono::steady_clock::now();
    for (unsigned i = 0; i < kIterations; ++i) {
        gf256_mul_mem(
            destination.Data, source.Data,
            static_cast<uint8_t>(2u + i % 253u),
            static_cast<int>(kBytes));
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    const double mib_per_second =
        static_cast<double>(kBytes) * kIterations /
        (1024.0 * 1024.0 * seconds);

    uint64_t checksum = 0;
    for (size_t i = 0; i < kBytes; i += 4093) {
        checksum = checksum * UINT64_C(0x100000001b3) ^ destination.Data[i];
    }
    std::printf(
        "GF256 out-of-place benchmark: %.2f MiB/s checksum=%" PRIu64 "\n",
        mib_per_second, checksum);

    static const size_t kMultiBytes = 1024u * 1024u;
    static const unsigned kMultiIterations = 256u;
    static const int kDestinationCount = 12;
    AlignedBuffer multi_source(kMultiBytes, 1u, 0x39u);
    std::vector<std::vector<uint8_t> > multi_storage(
        kDestinationCount,
        std::vector<uint8_t>(kMultiBytes + kAlignment, 0x5au));
    void* destinations[kDestinationCount];
    uint8_t scales[kDestinationCount];
    for (int j = 0; j < kDestinationCount; ++j)
    {
        destinations[j] = multi_storage[(size_t)j].data() + (j & 15);
        scales[j] = (uint8_t)(j * 19u + 2u);
    }

    const auto repeated_start = std::chrono::steady_clock::now();
    for (unsigned iteration = 0; iteration < kMultiIterations; ++iteration) {
        for (int j = 0; j < kDestinationCount; ++j) {
            gf256_muladd_mem(
                destinations[j], scales[j], multi_source.Data,
                (int)kMultiBytes);
        }
    }
    const double repeated_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - repeated_start).count();

    const auto fused_start = std::chrono::steady_clock::now();
    for (unsigned iteration = 0; iteration < kMultiIterations; ++iteration) {
        gf256_muladd_multi_mem(
            destinations, scales, kDestinationCount, multi_source.Data,
            (int)kMultiBytes);
    }
    const double fused_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - fused_start).count();
    uint64_t multi_checksum = 0u;
    for (int j = 0; j < kDestinationCount; ++j) {
        const uint8_t* const data =
            static_cast<const uint8_t*>(destinations[j]);
        for (size_t i = 0; i < kMultiBytes; i += 4093u) {
            multi_checksum =
                multi_checksum * UINT64_C(0x100000001b3) ^ data[i];
        }
    }
    std::printf(
        "GF256 12-destination benchmark: repeated=%.3fs fused=%.3fs "
        "speedup=%.3fx checksum=%" PRIu64 "\n",
        repeated_seconds,
        fused_seconds,
        repeated_seconds / fused_seconds,
        multi_checksum);
}

} // namespace

int main(int argc, char** argv)
{
    if (gf256_init() != 0) {
        std::fprintf(stderr, "gf256_init failed\n");
        return 1;
    }

    gf256_x86_cpu_features active = {};
    gf256_get_active_x86_cpu_features(&active);
    std::printf(
        "Active x86 kernels: SSSE3=%d AVX2=%d GFNI=%d AVX512=%d\n",
        active.SSSE3, active.AVX2, active.GFNI, active.AVX512);

    if (!TestOperation(Operation::Multiply) ||
        !TestOperation(Operation::Divide) ||
        !TestMultiDestination())
    {
        return 1;
    }

    if (argc == 2 && std::string(argv[1]) == "--benchmark") {
        RunBenchmark();
    } else if (argc != 1) {
        std::fprintf(stderr, "Usage: %s [--benchmark]\n", argv[0]);
        return 2;
    }

    std::puts("GF256 in-place differential test passed");
    return 0;
}
