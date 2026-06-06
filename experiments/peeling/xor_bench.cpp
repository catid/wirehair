// Block-XOR calibration benchmark for peeling-solve cost comparisons.
//
// The harness reports microseconds per block XOR for solve-like random block
// dependencies.  A "xor" here means dst_block ^= src_block for the configured
// piece size, matching the block-XOR unit used by peel_sweep's solve estimates.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Rng
{
    uint64_t State;

    explicit Rng(uint64_t seed) : State(seed) {}

    uint64_t next()
    {
        uint64_t z = (State += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
};

struct Op
{
    uint32_t Dst;
    uint32_t Src;
};

static double now_sec()
{
    return std::chrono::duration<double>(Clock::now().time_since_epoch()).count();
}

static uint64_t parse_u64(const char* text, bool* ok)
{
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 0);
    if (ok) {
        *ok = end != text && *end == '\0';
    }
    return (uint64_t)value;
}

static double parse_double(const char* text, bool* ok)
{
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (ok) {
        *ok = end != text && *end == '\0';
    }
    return value;
}

static std::vector<size_t> parse_size_list(const std::string& text, bool* ok)
{
    if (ok) {
        *ok = true;
    }
    std::vector<size_t> values;
    size_t start = 0;
    while (start < text.size())
    {
        const size_t comma = text.find(',', start);
        const std::string part = text.substr(
            start,
            comma == std::string::npos ? std::string::npos : comma - start);
        if (part.empty())
        {
            if (ok) {
                *ok = false;
            }
            return values;
        }

        bool part_ok = false;
        const uint64_t value = parse_u64(part.c_str(), &part_ok);
        if (!part_ok || value == 0)
        {
            if (ok) {
                *ok = false;
            }
            return values;
        }
        values.push_back((size_t)value);

        if (comma == std::string::npos) {
            break;
        }
        start = comma + 1;
    }
    return values;
}

static size_t default_working_blocks(size_t block_bytes)
{
    if (block_bytes <= 4096) {
        return 8192;
    }
    if (block_bytes >= 1024 * 1024) {
        return 128;
    }
    return 1024;
}

static void block_xor(uint8_t* dst, const uint8_t* src, size_t bytes)
{
    if ((((uintptr_t)dst | (uintptr_t)src | bytes) & 7u) != 0)
    {
        for (size_t i = 0; i < bytes; ++i) {
            dst[i] ^= src[i];
        }
        return;
    }

    const size_t words = bytes / sizeof(uint64_t);
    uint64_t* d64 = reinterpret_cast<uint64_t*>(dst);
    const uint64_t* s64 = reinterpret_cast<const uint64_t*>(src);
    for (size_t i = 0; i < words; ++i) {
        d64[i] ^= s64[i];
    }
    for (size_t i = words * sizeof(uint64_t); i < bytes; ++i) {
        dst[i] ^= src[i];
    }
}

static uint64_t checksum_buffer(const uint8_t* data, size_t bytes)
{
    uint64_t sum = UINT64_C(1469598103934665603);
    const size_t step = bytes > 4096 ? 4096 : 1;
    for (size_t i = 0; i < bytes; i += step)
    {
        sum ^= data[i];
        sum *= UINT64_C(1099511628211);
    }
    sum ^= data[bytes - 1];
    return sum;
}

static bool allocate_aligned(size_t bytes, uint8_t** out)
{
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 64, bytes) != 0) {
        return false;
    }
    *out = reinterpret_cast<uint8_t*>(ptr);
    return true;
}

static int run_one_size(
    size_t block_bytes,
    double target_gib,
    unsigned repeats,
    uint64_t seed)
{
    const size_t working_blocks = default_working_blocks(block_bytes);
    const size_t total_bytes = working_blocks * block_bytes;
    uint8_t* data = nullptr;
    if (!allocate_aligned(total_bytes, &data))
    {
        std::fprintf(stderr, "allocation failed for %zu bytes\n", total_bytes);
        return 1;
    }

    Rng init_rng(seed ^ (uint64_t)block_bytes);
    for (size_t i = 0; i < total_bytes; ++i) {
        data[i] = (uint8_t)init_rng.next();
    }

    const double target_bytes = target_gib * 1024.0 * 1024.0 * 1024.0;
    size_t op_count = (size_t)(target_bytes / (double)block_bytes);
    if (op_count < 1024) {
        op_count = 1024;
    }

    std::vector<Op> ops(op_count);
    Rng op_rng(seed ^ (uint64_t)block_bytes ^
        UINT64_C(0x786f725f736f6c76));
    for (size_t i = 0; i < op_count; ++i)
    {
        const uint32_t dst = (uint32_t)(op_rng.next() % working_blocks);
        uint32_t src = (uint32_t)(op_rng.next() % working_blocks);
        if (src == dst) {
            src = (src + 1u) % (uint32_t)working_blocks;
        }
        ops[i].Dst = dst;
        ops[i].Src = src;
    }

    const size_t warmup_ops = std::min<size_t>(ops.size(), 4096);
    for (size_t i = 0; i < warmup_ops; ++i)
    {
        block_xor(
            data + (size_t)ops[i].Dst * block_bytes,
            data + (size_t)ops[i].Src * block_bytes,
            block_bytes);
    }

    std::vector<double> usec_per_xor;
    std::vector<double> total_usec;
    std::vector<double> gib_per_sec;
    usec_per_xor.reserve(repeats);
    total_usec.reserve(repeats);
    gib_per_sec.reserve(repeats);

    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        const double start = now_sec();
        for (const Op& op : ops)
        {
            block_xor(
                data + (size_t)op.Dst * block_bytes,
                data + (size_t)op.Src * block_bytes,
                block_bytes);
        }
        const double elapsed = now_sec() - start;
        const double usec = elapsed * 1000000.0;
        total_usec.push_back(usec);
        usec_per_xor.push_back(usec / (double)ops.size());

        const double traffic_gib =
            (2.0 * (double)block_bytes * (double)ops.size()) /
            (1024.0 * 1024.0 * 1024.0);
        gib_per_sec.push_back(traffic_gib / elapsed);
    }

    std::sort(usec_per_xor.begin(), usec_per_xor.end());
    std::sort(total_usec.begin(), total_usec.end());
    std::sort(gib_per_sec.begin(), gib_per_sec.end());
    const size_t mid = usec_per_xor.size() / 2u;
    const uint64_t sum = checksum_buffer(data, total_bytes);
    std::printf("%zu,%zu,%zu,%u,%.3f,%.9f,%.3f,0x%016llx\n",
        block_bytes,
        working_blocks,
        ops.size(),
        repeats,
        total_usec[mid],
        usec_per_xor[mid],
        gib_per_sec[mid],
        (unsigned long long)sum);

    std::free(data);
    return 0;
}

static void usage(const char* argv0)
{
    std::fprintf(stderr,
        "usage: %s [--sizes 1280,1048576] [--target-gib 16] "
        "[--repeats 3] [--seed N]\n",
        argv0);
}

} // namespace

int main(int argc, char** argv)
{
    std::string sizes_spec = "1280,1048576";
    double target_gib = 16.0;
    unsigned repeats = 3;
    uint64_t seed = UINT64_C(0x51f15eed12345678);

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--sizes" && i + 1 < argc) {
            sizes_spec = argv[++i];
        }
        else if (arg == "--target-gib" && i + 1 < argc)
        {
            bool ok = false;
            target_gib = parse_double(argv[++i], &ok);
            if (!ok || target_gib <= 0.0)
            {
                usage(argv[0]);
                return 1;
            }
        }
        else if (arg == "--repeats" && i + 1 < argc)
        {
            bool ok = false;
            const uint64_t value = parse_u64(argv[++i], &ok);
            if (!ok || value == 0 || value > 1000)
            {
                usage(argv[0]);
                return 1;
            }
            repeats = (unsigned)value;
        }
        else if (arg == "--seed" && i + 1 < argc)
        {
            bool ok = false;
            seed = parse_u64(argv[++i], &ok);
            if (!ok)
            {
                usage(argv[0]);
                return 1;
            }
        }
        else
        {
            usage(argv[0]);
            return 1;
        }
    }

    bool sizes_ok = false;
    const std::vector<size_t> sizes = parse_size_list(sizes_spec, &sizes_ok);
    if (!sizes_ok || sizes.empty())
    {
        usage(argv[0]);
        return 1;
    }

    std::printf("block_bytes,working_blocks,ops_per_repeat,repeats,"
        "median_total_us,median_us_per_xor,median_gib_per_s,checksum\n");
    for (size_t block_bytes : sizes)
    {
        const int rc = run_one_size(block_bytes, target_gib, repeats, seed);
        if (rc != 0) {
            return rc;
        }
    }
    return 0;
}
