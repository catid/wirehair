// Standalone benchmark for replaying block-data row operations by byte tiles.

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

    uint32_t uniform(uint32_t n)
    {
        return n <= 1 ? 0 : (uint32_t)(next() % n);
    }
};

struct Op
{
    uint32_t Dst;
    uint32_t Src;
};

struct Case
{
    const char* Name;
    size_t BlockBytes;
    uint32_t Blocks;
    uint32_t OpsPerSource;
    uint32_t SourcePasses;
};

static double now_sec()
{
    return std::chrono::duration<double>(Clock::now().time_since_epoch()).count();
}

static void xor_words(uint64_t* dst, const uint64_t* src, size_t words)
{
    for (size_t i = 0; i < words; ++i) {
        dst[i] ^= src[i];
    }
}

static uint64_t checksum_words(const std::vector<uint64_t>& data)
{
    uint64_t h = UINT64_C(1469598103934665603);
    for (uint64_t x : data)
    {
        h ^= x;
        h *= UINT64_C(1099511628211);
    }
    return h;
}

static void fill_data(std::vector<uint64_t>& data, uint64_t seed)
{
    Rng rng(seed);
    for (uint64_t& x : data) {
        x = rng.next();
    }
}

static std::vector<Op> make_schedule(
    uint32_t blocks,
    uint32_t ops_per_source,
    uint32_t source_passes,
    uint64_t seed)
{
    Rng rng(seed);
    std::vector<Op> ops;
    ops.reserve((size_t)blocks * source_passes * ops_per_source);

    for (uint32_t pass = 0; pass < source_passes; ++pass)
    {
        for (uint32_t src = 0; src < blocks; ++src)
        {
            for (uint32_t j = 0; j < ops_per_source; ++j)
            {
                uint32_t dst = rng.uniform(blocks);
                if (dst == src) {
                    dst = (dst + 1) % blocks;
                }
                ops.push_back(Op{dst, src});
            }
        }
    }

    return ops;
}

static void replay_untiled(
    std::vector<uint64_t>& data,
    const std::vector<Op>& ops,
    size_t words_per_block)
{
    for (const Op& op : ops)
    {
        uint64_t* dst = &data[(size_t)op.Dst * words_per_block];
        const uint64_t* src = &data[(size_t)op.Src * words_per_block];
        xor_words(dst, src, words_per_block);
    }
}

static void replay_tiled(
    std::vector<uint64_t>& data,
    const std::vector<Op>& ops,
    size_t words_per_block,
    size_t tile_bytes)
{
    size_t tile_words = tile_bytes / sizeof(uint64_t);
    if (tile_words == 0) {
        tile_words = 1;
    }

    for (size_t offset = 0; offset < words_per_block; offset += tile_words)
    {
        const size_t words = std::min(tile_words, words_per_block - offset);
        for (const Op& op : ops)
        {
            uint64_t* dst = &data[(size_t)op.Dst * words_per_block + offset];
            const uint64_t* src = &data[(size_t)op.Src * words_per_block + offset];
            xor_words(dst, src, words);
        }
    }
}

static double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

static bool verify_case(size_t block_bytes, size_t tile_bytes)
{
    const uint32_t blocks = 64;
    const size_t words_per_block = block_bytes / sizeof(uint64_t);
    std::vector<uint64_t> base((size_t)blocks * words_per_block);
    fill_data(base, UINT64_C(0x766572696679));
    const std::vector<Op> ops = make_schedule(
        blocks, 8, 2, UINT64_C(0x7363686564756c65));

    std::vector<uint64_t> untiled = base;
    std::vector<uint64_t> tiled = base;
    replay_untiled(untiled, ops, words_per_block);
    replay_tiled(tiled, ops, words_per_block, tile_bytes);
    return untiled == tiled;
}

static double time_untiled(
    const std::vector<uint64_t>& base,
    const std::vector<Op>& ops,
    size_t words_per_block,
    unsigned repeats,
    uint64_t* checksum)
{
    std::vector<double> samples;
    samples.reserve(repeats);
    std::vector<uint64_t> data;
    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        data = base;
        const double start = now_sec();
        replay_untiled(data, ops, words_per_block);
        samples.push_back(now_sec() - start);
    }
    *checksum = checksum_words(data);
    return median(samples);
}

static double time_tiled(
    const std::vector<uint64_t>& base,
    const std::vector<Op>& ops,
    size_t words_per_block,
    size_t tile_bytes,
    unsigned repeats,
    uint64_t* checksum)
{
    std::vector<double> samples;
    samples.reserve(repeats);
    std::vector<uint64_t> data;
    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        data = base;
        const double start = now_sec();
        replay_tiled(data, ops, words_per_block, tile_bytes);
        samples.push_back(now_sec() - start);
    }
    *checksum = checksum_words(data);
    return median(samples);
}

static void print_result(
    const Case& c,
    const char* mode,
    size_t tile_bytes,
    size_t op_count,
    double seconds,
    uint64_t checksum)
{
    const double logical_gib =
        (double)op_count * (double)c.BlockBytes / 1073741824.0;
    const double memory_gib_3stream = 3.0 * logical_gib;
    const double gib_per_s = memory_gib_3stream / seconds;
    std::printf(
        "%s,%zu,%u,%zu,%s,%zu,%.6f,%.6f,%.6f,%.3f,0x%016llx\n",
        c.Name,
        c.BlockBytes,
        c.Blocks,
        op_count,
        mode,
        tile_bytes,
        seconds * 1000.0,
        logical_gib,
        memory_gib_3stream,
        gib_per_s,
        (unsigned long long)checksum);
}

} // namespace

int main(int argc, char** argv)
{
    const Case cases[] = {
        {"mtu1280", 1280, 32768, 8, 1},
        {"kib100", 102400, 2048, 8, 1},
        {"mib1", 1048576, 256, 8, 1},
    };
    const size_t tiles[] = {
        256,        // sub-block tiles so the small-block cases measure
        512,        // a real tiled replay instead of degenerating to untiled
        16 * 1024,
        32 * 1024,
        64 * 1024,
        128 * 1024,
        256 * 1024,
        512 * 1024,
    };
    const unsigned repeats = 3;
    bool verify_only = false;

    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--verify-only")) {
            verify_only = true;
        }
        else
        {
            std::fprintf(stderr, "usage: rowop_tiling [--verify-only]\n");
            return 1;
        }
    }

    for (const Case& c : cases)
    {
        if (c.BlockBytes % sizeof(uint64_t) != 0)
        {
            std::fprintf(stderr, "block size must be a multiple of 8: %zu\n",
                c.BlockBytes);
            return 1;
        }
        for (size_t tile : tiles)
        {
            if (!verify_case(c.BlockBytes, tile))
            {
                std::fprintf(stderr, "verification failed for %s tile=%zu\n",
                    c.Name, tile);
                return 1;
            }
        }
    }
    if (verify_only)
    {
        std::printf("verify: ok\n");
        return 0;
    }

    std::printf(
        "case,block_bytes,blocks,ops,mode,tile_bytes,median_ms,"
        "logical_gib,memory_gib_3stream,gib_per_s,checksum\n");

    for (const Case& c : cases)
    {
        const size_t words_per_block = c.BlockBytes / sizeof(uint64_t);
        std::vector<uint64_t> base((size_t)c.Blocks * words_per_block);
        fill_data(base, UINT64_C(0x646174615f626173) ^ (uint64_t)c.BlockBytes);
        const std::vector<Op> ops = make_schedule(
            c.Blocks, c.OpsPerSource, c.SourcePasses,
            UINT64_C(0x6f70735f73656564) ^ (uint64_t)c.BlockBytes);

        uint64_t checksum = 0;
        const double untiled = time_untiled(
            base, ops, words_per_block, repeats, &checksum);
        print_result(c, "untiled", 0, ops.size(), untiled, checksum);

        for (size_t tile : tiles)
        {
            // A tile at or above the block size replays the identical loop
            // as untiled; reporting it as "tiled" would present run-to-run
            // noise as a tiling result.
            if (tile >= c.BlockBytes) {
                continue;
            }
            const double tiled = time_tiled(
                base, ops, words_per_block, tile, repeats, &checksum);
            print_result(c, "tiled", tile, ops.size(), tiled, checksum);
        }
    }

    return 0;
}
