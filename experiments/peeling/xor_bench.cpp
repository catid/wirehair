// Block-XOR calibration benchmark for peeling-solve cost comparisons.
//
// The harness reports microseconds per block XOR for solve-like random block
// dependencies.  A "xor" here means dst_block ^= src_block for the configured
// piece size, matching the block-XOR unit used by peel_sweep's solve estimates.
//
// Extensions:
//   --cold    : DRAM-cold variant of the block-XOR calibration.  The default
//               protocol walks a --working-mib (256 MiB) working set, which on
//               large-L3 parts is partially cache-resident; the cold mode
//               walks a random permutation of a >= 4 GiB pool so every op
//               touches lines that were evicted long ago.  Same CSV columns
//               as the default mode.
//   --muladd  : GF(256) dst[] ^= c * src[] calibration using the production
//               gf256_muladd_mem() routine (nonzero random c per op), run
//               against the block-XOR kernel under identical cold-pool
//               conditions so the muladd:xor ratio is meaningful.
//   --fanin   : k-ary gather-XOR microbenchmark, dst ^= s1 ^ ... ^ sk,
//               comparing k chained pairwise passes against one fused pass.
//               Logical traffic is defined as (k+1)*block_bytes per op for
//               BOTH arms so the throughput ratio is comparable.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "gf256.h"

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
        *ok = end != text && *end == '\0' && std::isfinite(value);
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

static size_t working_blocks_for(size_t block_bytes, size_t working_mib)
{
    const size_t mib = 1024u * 1024u;
    if (working_mib > SIZE_MAX / mib) {
        return 0;
    }
    const size_t target_bytes = working_mib * mib;
    if (target_bytes > SIZE_MAX - block_bytes + 1u) {
        return 0;
    }
    size_t blocks = (target_bytes + block_bytes - 1u) / block_bytes;
    if (blocks < 2u) {
        blocks = 2u;
    }
    return blocks;
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

// Byte-wise k-ary gather XOR tail/fallback: dst[i] ^= srcs[0][i]^...^srcs[k-1][i].
static void block_fanin_bytes(
    uint8_t* dst, const uint8_t* const* srcs, size_t k, size_t bytes,
    size_t begin)
{
    for (size_t i = begin; i < bytes; ++i)
    {
        uint8_t acc = dst[i];
        for (size_t j = 0; j < k; ++j) {
            acc ^= srcs[j][i];
        }
        dst[i] = acc;
    }
}

// Runtime-k fused gather XOR over 64-bit words.
static void block_fanin_generic(
    uint8_t* dst, const uint8_t* const* srcs, size_t k, size_t bytes)
{
    uintptr_t mis = (uintptr_t)dst | (uintptr_t)bytes;
    for (size_t j = 0; j < k; ++j) {
        mis |= (uintptr_t)srcs[j];
    }
    if ((mis & 7u) != 0)
    {
        block_fanin_bytes(dst, srcs, k, bytes, 0);
        return;
    }

    const size_t words = bytes / sizeof(uint64_t);
    uint64_t* d64 = reinterpret_cast<uint64_t*>(dst);
    for (size_t i = 0; i < words; ++i)
    {
        uint64_t acc = d64[i];
        for (size_t j = 0; j < k; ++j) {
            acc ^= reinterpret_cast<const uint64_t*>(srcs[j])[i];
        }
        d64[i] = acc;
    }
    block_fanin_bytes(dst, srcs, k, bytes, words * sizeof(uint64_t));
}

// Compile-time-k fused gather XOR so -O3 can unroll/vectorize the inner loop.
template <unsigned kFanIn>
static void block_fanin_fixed(
    uint8_t* dst, const uint8_t* const* srcs, size_t bytes)
{
    uintptr_t mis = (uintptr_t)dst | (uintptr_t)bytes;
    for (unsigned j = 0; j < kFanIn; ++j) {
        mis |= (uintptr_t)srcs[j];
    }
    if ((mis & 7u) != 0)
    {
        block_fanin_bytes(dst, srcs, kFanIn, bytes, 0);
        return;
    }

    const uint64_t* s64[kFanIn];
    for (unsigned j = 0; j < kFanIn; ++j) {
        s64[j] = reinterpret_cast<const uint64_t*>(srcs[j]);
    }
    const size_t words = bytes / sizeof(uint64_t);
    uint64_t* d64 = reinterpret_cast<uint64_t*>(dst);
    for (size_t i = 0; i < words; ++i)
    {
        uint64_t acc = d64[i];
        for (unsigned j = 0; j < kFanIn; ++j) {
            acc ^= s64[j][i];
        }
        d64[i] = acc;
    }
    block_fanin_bytes(dst, srcs, kFanIn, bytes, words * sizeof(uint64_t));
}

// Fused single-pass gather XOR entry point: dst ^= srcs[0] ^ ... ^ srcs[k-1].
static void block_fanin_gather(
    uint8_t* dst, const uint8_t* const* srcs, size_t k, size_t bytes)
{
    switch (k)
    {
    case 2:  block_fanin_fixed<2>(dst, srcs, bytes);  return;
    case 3:  block_fanin_fixed<3>(dst, srcs, bytes);  return;
    case 4:  block_fanin_fixed<4>(dst, srcs, bytes);  return;
    case 5:  block_fanin_fixed<5>(dst, srcs, bytes);  return;
    case 6:  block_fanin_fixed<6>(dst, srcs, bytes);  return;
    case 7:  block_fanin_fixed<7>(dst, srcs, bytes);  return;
    case 8:  block_fanin_fixed<8>(dst, srcs, bytes);  return;
    case 12: block_fanin_fixed<12>(dst, srcs, bytes); return;
    case 16: block_fanin_fixed<16>(dst, srcs, bytes); return;
    default: block_fanin_generic(dst, srcs, k, bytes); return;
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

//------------------------------------------------------------------------------
// Cold-pool infrastructure (shared by --cold / --muladd / --fanin)

static const double kGiBDouble = 1024.0 * 1024.0 * 1024.0;

static uint8_t* allocate_pool(size_t pool_bytes, uint64_t seed)
{
    uint8_t* data = nullptr;
    if (!allocate_aligned(pool_bytes, &data)) {
        return nullptr;
    }

    Rng rng(seed ^ UINT64_C(0xc01db007c01db007));
    uint64_t* words = reinterpret_cast<uint64_t*>(data);
    const size_t word_count = pool_bytes / sizeof(uint64_t);
    for (size_t i = 0; i < word_count; ++i) {
        words[i] = rng.next();
    }
    for (size_t i = word_count * sizeof(uint64_t); i < pool_bytes; ++i) {
        data[i] = (uint8_t)rng.next();
    }
    return data;
}

// Fisher-Yates permutation of block indices: ops walk this so each op touches
// lines that were last referenced a full pool-pass ago (>= pool_bytes of
// intervening traffic), defeating L2/L3 caching.
static std::vector<uint32_t> make_block_permutation(size_t n_blocks, Rng& rng)
{
    std::vector<uint32_t> perm(n_blocks);
    for (size_t i = 0; i < n_blocks; ++i) {
        perm[i] = (uint32_t)i;
    }
    for (size_t i = n_blocks - 1; i > 0; --i)
    {
        const size_t j = (size_t)(rng.next() % (i + 1));
        std::swap(perm[i], perm[j]);
    }
    return perm;
}

static size_t ops_for_target(
    double target_gib, size_t bytes_per_op, size_t min_ops)
{
    const double target_bytes = target_gib * kGiBDouble;
    const double count = target_bytes / (double)bytes_per_op;
    size_t ops = count < 1.0 ? (size_t)1 : (size_t)count;
    if (ops < min_ops) {
        ops = min_ops;
    }
    return ops;
}

static double median_of(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2u];
}

//------------------------------------------------------------------------------
// Functional self-checks (run before any measurement; abort on mismatch)

static bool self_check_xor()
{
    static const size_t kSizes[2] = { 1280, 1283 };
    for (size_t n : kSizes)
    {
        std::vector<uint8_t> dst(n), src(n), ref(n);
        Rng rng(UINT64_C(0x5e1fc4ec0000) ^ (uint64_t)n);
        for (size_t i = 0; i < n; ++i)
        {
            dst[i] = (uint8_t)rng.next();
            src[i] = (uint8_t)rng.next();
            ref[i] = (uint8_t)(dst[i] ^ src[i]);
        }
        block_xor(dst.data(), src.data(), n);
        if (std::memcmp(dst.data(), ref.data(), n) != 0)
        {
            std::fprintf(stderr,
                "self-check FAILED: block_xor mismatch at %zu bytes\n", n);
            return false;
        }
    }
    return true;
}

static bool self_check_muladd()
{
    static const size_t kSizes[2] = { 1280, 1283 };
    static const uint8_t kMultipliers[3] = { 2, 0x8e, 0xff };
    for (size_t n : kSizes)
    {
        for (uint8_t c : kMultipliers)
        {
            std::vector<uint8_t> dst(n), src(n), ref(n);
            Rng rng(UINT64_C(0x5e1fc4ecf256) ^ (uint64_t)n ^ (uint64_t)c);
            for (size_t i = 0; i < n; ++i)
            {
                dst[i] = (uint8_t)rng.next();
                src[i] = (uint8_t)rng.next();
                ref[i] = (uint8_t)(dst[i] ^ gf256_mul(src[i], c));
            }
            gf256_muladd_mem(dst.data(), c, src.data(), (int)n);
            if (std::memcmp(dst.data(), ref.data(), n) != 0)
            {
                std::fprintf(stderr,
                    "self-check FAILED: gf256_muladd_mem mismatch at "
                    "%zu bytes, c=0x%02x (vs scalar GF(256) reference)\n",
                    n, c);
                return false;
            }
        }
    }
    return true;
}

static bool self_check_fanin(size_t k)
{
    static const size_t kSizes[2] = { 1280, 1283 };
    for (size_t n : kSizes)
    {
        Rng rng(UINT64_C(0x5e1fc4ecfa41) ^ (uint64_t)n ^ ((uint64_t)k << 32));
        std::vector<uint8_t> base(n);
        for (size_t i = 0; i < n; ++i) {
            base[i] = (uint8_t)rng.next();
        }

        std::vector<std::vector<uint8_t> > srcs(k, std::vector<uint8_t>(n));
        std::vector<const uint8_t*> src_ptrs(k);
        for (size_t j = 0; j < k; ++j)
        {
            for (size_t i = 0; i < n; ++i) {
                srcs[j][i] = (uint8_t)rng.next();
            }
            src_ptrs[j] = srcs[j].data();
        }

        // Chained arm reference: k pairwise passes.
        std::vector<uint8_t> chained = base;
        for (size_t j = 0; j < k; ++j) {
            block_xor(chained.data(), src_ptrs[j], n);
        }

        // Fused arm under test.
        std::vector<uint8_t> fused = base;
        block_fanin_gather(fused.data(), src_ptrs.data(), k, n);

        if (std::memcmp(chained.data(), fused.data(), n) != 0)
        {
            std::fprintf(stderr,
                "self-check FAILED: fan-in fused vs chained mismatch at "
                "k=%zu, %zu bytes\n", k, n);
            return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------
// --cold: DRAM-cold block-XOR calibration (same CSV columns as default mode)

static int run_cold_xor_size(
    uint8_t* pool,
    size_t pool_bytes,
    size_t block_bytes,
    double target_gib,
    unsigned repeats,
    uint64_t seed)
{
    const size_t n_blocks = pool_bytes / block_bytes;
    if (n_blocks < 2)
    {
        std::fprintf(stderr,
            "pool too small for block size %zu (need >= 2 blocks)\n",
            block_bytes);
        return 1;
    }

    Rng rng(seed ^ (uint64_t)block_bytes ^ UINT64_C(0xc01d5eedc01d5eed));
    const std::vector<uint32_t> perm = make_block_permutation(n_blocks, rng);

    const size_t op_count = ops_for_target(target_gib, block_bytes, 1024);
    std::vector<Op> ops(op_count);
    for (size_t i = 0; i < op_count; ++i)
    {
        ops[i].Dst = perm[(2 * i) % n_blocks];
        ops[i].Src = perm[(2 * i + 1) % n_blocks];
    }

    // Warm up using the TAIL of the op list: the timed loop reaches those
    // blocks last, after >= pool-pass traffic has evicted them again.
    const size_t warmup_ops = std::min<size_t>(op_count, 256);
    for (size_t i = op_count - warmup_ops; i < op_count; ++i)
    {
        block_xor(
            pool + (size_t)ops[i].Dst * block_bytes,
            pool + (size_t)ops[i].Src * block_bytes,
            block_bytes);
    }

    std::vector<double> usec_per_xor, total_usec, gib_per_sec;
    usec_per_xor.reserve(repeats);
    total_usec.reserve(repeats);
    gib_per_sec.reserve(repeats);

    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        const double start = now_sec();
        for (const Op& op : ops)
        {
            block_xor(
                pool + (size_t)op.Dst * block_bytes,
                pool + (size_t)op.Src * block_bytes,
                block_bytes);
        }
        const double elapsed = now_sec() - start;
        const double usec = elapsed * 1000000.0;
        total_usec.push_back(usec);
        usec_per_xor.push_back(usec / (double)ops.size());

        const double traffic_gib =
            (2.0 * (double)block_bytes * (double)ops.size()) / kGiBDouble;
        gib_per_sec.push_back(traffic_gib / elapsed);
    }

    const uint64_t sum = checksum_buffer(pool, pool_bytes);
    std::printf("%zu,%zu,%zu,%u,%.3f,%.9f,%.3f,0x%016llx\n",
        block_bytes,
        n_blocks,
        ops.size(),
        repeats,
        median_of(total_usec),
        median_of(usec_per_xor),
        median_of(gib_per_sec),
        (unsigned long long)sum);
    return 0;
}

//------------------------------------------------------------------------------
// --muladd: GF(256) muladd calibration vs block-XOR under identical conditions

struct MulOp
{
    uint32_t Dst;
    uint32_t Src;
    uint8_t C;
};

static int run_muladd_size(
    uint8_t* pool,
    size_t pool_bytes,
    size_t block_bytes,
    double target_gib,
    unsigned repeats,
    uint64_t seed)
{
    if (block_bytes > (size_t)INT32_MAX)
    {
        std::fprintf(stderr,
            "block size %zu too large for gf256_muladd_mem\n", block_bytes);
        return 1;
    }
    const size_t n_blocks = pool_bytes / block_bytes;
    if (n_blocks < 2)
    {
        std::fprintf(stderr,
            "pool too small for block size %zu (need >= 2 blocks)\n",
            block_bytes);
        return 1;
    }

    Rng rng(seed ^ (uint64_t)block_bytes ^ UINT64_C(0xf256f256f256f256));
    const std::vector<uint32_t> perm = make_block_permutation(n_blocks, rng);

    const size_t op_count = ops_for_target(target_gib, block_bytes, 1024);
    std::vector<MulOp> ops(op_count);
    for (size_t i = 0; i < op_count; ++i)
    {
        ops[i].Dst = perm[(2 * i) % n_blocks];
        ops[i].Src = perm[(2 * i + 1) % n_blocks];
        // Nonzero, non-1 multiplier: gf256_muladd_mem special-cases 0 and 1.
        ops[i].C = (uint8_t)(2 + (rng.next() % 254));
    }

    // Warm up on the tail of the op list (see run_cold_xor_size).
    const size_t warmup_ops = std::min<size_t>(op_count, 256);
    for (size_t i = op_count - warmup_ops; i < op_count; ++i)
    {
        block_xor(
            pool + (size_t)ops[i].Dst * block_bytes,
            pool + (size_t)ops[i].Src * block_bytes,
            block_bytes);
    }

    std::vector<double> xor_usec, muladd_usec;
    xor_usec.reserve(repeats);
    muladd_usec.reserve(repeats);

    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        // XOR arm (same kernel as the default/--cold calibration).
        double start = now_sec();
        for (const MulOp& op : ops)
        {
            block_xor(
                pool + (size_t)op.Dst * block_bytes,
                pool + (size_t)op.Src * block_bytes,
                block_bytes);
        }
        xor_usec.push_back((now_sec() - start) * 1000000.0 /
            (double)ops.size());

        // GF(256) muladd arm: dst[] ^= c * src[].
        start = now_sec();
        for (const MulOp& op : ops)
        {
            gf256_muladd_mem(
                pool + (size_t)op.Dst * block_bytes,
                op.C,
                pool + (size_t)op.Src * block_bytes,
                (int)block_bytes);
        }
        muladd_usec.push_back((now_sec() - start) * 1000000.0 /
            (double)ops.size());
    }

    const double xor_med = median_of(xor_usec);
    const double muladd_med = median_of(muladd_usec);
    const double traffic_gib_per_op =
        (2.0 * (double)block_bytes) / kGiBDouble;
    const uint64_t sum = checksum_buffer(pool, pool_bytes);

    std::printf("%zu,%zu,%zu,%u,%.9f,%.3f,%.9f,%.3f,%.4f,0x%016llx\n",
        block_bytes,
        n_blocks,
        ops.size(),
        repeats,
        xor_med,
        traffic_gib_per_op / (xor_med / 1000000.0),
        muladd_med,
        traffic_gib_per_op / (muladd_med / 1000000.0),
        muladd_med / xor_med,
        (unsigned long long)sum);
    return 0;
}

//------------------------------------------------------------------------------
// --fanin: k-ary gather XOR, chained pairwise passes vs one fused pass

static int run_fanin_size(
    uint8_t* pool,
    size_t pool_bytes,
    size_t k,
    size_t block_bytes,
    double target_gib,
    unsigned repeats,
    uint64_t seed)
{
    const size_t n_blocks = pool_bytes / block_bytes;
    if (n_blocks < k + 2)
    {
        std::fprintf(stderr,
            "pool too small for block size %zu at fan-in %zu "
            "(need >= %zu blocks, have %zu)\n",
            block_bytes, k, k + 2, n_blocks);
        return 1;
    }

    Rng rng(seed ^ (uint64_t)block_bytes ^
        ((uint64_t)k << 40) ^ UINT64_C(0xfa414fa414fa414f));
    const std::vector<uint32_t> perm = make_block_permutation(n_blocks, rng);

    // Logical traffic per op: k source reads + 1 destination write.
    const size_t stride = k + 1;
    const size_t bytes_per_op = stride * block_bytes;
    const size_t op_count = ops_for_target(target_gib, bytes_per_op, 64);

    // idx[i*stride] is the destination block, the next k entries are sources.
    // Consecutive permutation slots are distinct, so dst != any src.
    std::vector<uint32_t> idx(op_count * stride);
    for (size_t i = 0; i < idx.size(); ++i) {
        idx[i] = perm[i % n_blocks];
    }

    std::vector<const uint8_t*> srcs(k);

    // Warm up on the tail of the op list (see run_cold_xor_size).
    const size_t warmup_ops = std::min<size_t>(op_count, 64);
    for (size_t i = op_count - warmup_ops; i < op_count; ++i)
    {
        uint8_t* dst = pool + (size_t)idx[i * stride] * block_bytes;
        for (size_t j = 0; j < k; ++j)
        {
            block_xor(dst,
                pool + (size_t)idx[i * stride + 1 + j] * block_bytes,
                block_bytes);
        }
    }

    std::vector<double> chained_usec, gather_usec;
    chained_usec.reserve(repeats);
    gather_usec.reserve(repeats);

    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        // Chained arm: k separate full passes over the destination block.
        double start = now_sec();
        for (size_t i = 0; i < op_count; ++i)
        {
            uint8_t* dst = pool + (size_t)idx[i * stride] * block_bytes;
            for (size_t j = 0; j < k; ++j)
            {
                block_xor(dst,
                    pool + (size_t)idx[i * stride + 1 + j] * block_bytes,
                    block_bytes);
            }
        }
        chained_usec.push_back((now_sec() - start) * 1000000.0 /
            (double)op_count);

        // Gather arm: one fused pass reading k sources, writing dst once.
        start = now_sec();
        for (size_t i = 0; i < op_count; ++i)
        {
            uint8_t* dst = pool + (size_t)idx[i * stride] * block_bytes;
            for (size_t j = 0; j < k; ++j)
            {
                srcs[j] =
                    pool + (size_t)idx[i * stride + 1 + j] * block_bytes;
            }
            block_fanin_gather(dst, srcs.data(), k, block_bytes);
        }
        gather_usec.push_back((now_sec() - start) * 1000000.0 /
            (double)op_count);
    }

    const double chained_med = median_of(chained_usec);
    const double gather_med = median_of(gather_usec);
    const double logical_gib_per_op = (double)bytes_per_op / kGiBDouble;
    const uint64_t sum = checksum_buffer(pool, pool_bytes);

    std::printf("%zu,%zu,%zu,%zu,%u,%.9f,%.3f,%.9f,%.3f,%.4f,0x%016llx\n",
        k,
        block_bytes,
        n_blocks,
        op_count,
        repeats,
        chained_med,
        logical_gib_per_op / (chained_med / 1000000.0),
        gather_med,
        logical_gib_per_op / (gather_med / 1000000.0),
        chained_med / gather_med,
        (unsigned long long)sum);
    return 0;
}

//------------------------------------------------------------------------------
// Default mode (unchanged): --working-mib working set, random src/dst pairs

static int run_one_size(
    size_t block_bytes,
    double target_gib,
    unsigned repeats,
    size_t working_mib,
    uint64_t seed)
{
    const size_t working_blocks = working_blocks_for(block_bytes, working_mib);
    if (working_blocks == 0 || working_blocks > SIZE_MAX / block_bytes)
    {
        std::fprintf(stderr, "working set is too large for block size %zu\n",
            block_bytes);
        return 1;
    }
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
        "usage: %s [--sizes 1280,102400,1048576] [--target-gib 16] "
        "[--repeats 3] [--working-mib 256] [--seed N]\n"
        "          [--cold] [--muladd] [--fanin 2,4,8] [--pool-gib 4]\n"
        "\n"
        "  default   : block-XOR over a --working-mib working set "
        "(unchanged legacy protocol)\n"
        "  --cold    : DRAM-cold block-XOR over a --pool-gib pool "
        "(permutation walk; same CSV columns)\n"
        "  --muladd  : gf256_muladd_mem vs block-XOR on the cold pool; "
        "reports muladd:xor ratio\n"
        "  --fanin K : k-ary gather XOR (chained vs fused) on the cold pool "
        "for each k in the list\n"
        "  --pool-gib: cold pool size in GiB (default 4; >= 4 recommended "
        "to defeat L3 caching)\n"
        "\n"
        "Any of --cold/--muladd/--fanin suppresses the default mode.\n",
        argv0);
}

} // namespace

int main(int argc, char** argv)
{
    std::string sizes_spec = "1280,102400,1048576";
    double target_gib = 16.0;
    unsigned repeats = 3;
    size_t working_mib = 256;
    uint64_t seed = UINT64_C(0x51f15eed12345678);
    bool cold_mode = false;
    bool muladd_mode = false;
    std::string fanin_spec;
    double pool_gib = 4.0;

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
        else if (arg == "--working-mib" && i + 1 < argc)
        {
            bool ok = false;
            const uint64_t value = parse_u64(argv[++i], &ok);
            if (!ok || value == 0 || value > SIZE_MAX)
            {
                usage(argv[0]);
                return 1;
            }
            working_mib = (size_t)value;
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
        else if (arg == "--cold") {
            cold_mode = true;
        }
        else if (arg == "--muladd") {
            muladd_mode = true;
        }
        else if (arg == "--fanin" && i + 1 < argc) {
            fanin_spec = argv[++i];
        }
        else if (arg == "--pool-gib" && i + 1 < argc)
        {
            bool ok = false;
            pool_gib = parse_double(argv[++i], &ok);
            if (!ok || pool_gib < 0.5)
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

    std::vector<size_t> fanin_ks;
    if (!fanin_spec.empty())
    {
        bool fanin_ok = false;
        fanin_ks = parse_size_list(fanin_spec, &fanin_ok);
        if (!fanin_ok || fanin_ks.empty())
        {
            usage(argv[0]);
            return 1;
        }
        for (size_t k : fanin_ks)
        {
            if (k < 2 || k > 256)
            {
                std::fprintf(stderr, "--fanin k must be in [2, 256]\n");
                return 1;
            }
        }
    }

    const bool pool_modes = cold_mode || muladd_mode || !fanin_ks.empty();
    if (!pool_modes)
    {
        // Legacy default mode, unchanged.
        std::printf("block_bytes,working_blocks,ops_per_repeat,repeats,"
            "median_total_us,median_us_per_xor,median_gib_per_s,checksum\n");
        for (size_t block_bytes : sizes)
        {
            const int rc = run_one_size(
                block_bytes, target_gib, repeats, working_mib, seed);
            if (rc != 0) {
                return rc;
            }
        }
        return 0;
    }

    if (muladd_mode && gf256_init() != 0)
    {
        std::fprintf(stderr, "gf256_init failed\n");
        return 1;
    }

    // Functional self-checks before any measurement.
    if (!self_check_xor()) {
        return 1;
    }
    if (muladd_mode && !self_check_muladd()) {
        return 1;
    }
    for (size_t k : fanin_ks)
    {
        if (!self_check_fanin(k)) {
            return 1;
        }
    }
    std::fprintf(stderr, "self-checks passed\n");

    if (pool_gib < 4.0)
    {
        std::fprintf(stderr,
            "warning: --pool-gib %.2f < 4; results may be cache-inflated\n",
            pool_gib);
    }
    size_t pool_bytes = (size_t)(pool_gib * kGiBDouble);
    pool_bytes -= pool_bytes % 4096u;
    std::fprintf(stderr, "allocating %.2f GiB cold pool...\n",
        (double)pool_bytes / kGiBDouble);
    uint8_t* pool = allocate_pool(pool_bytes, seed);
    if (pool == nullptr)
    {
        std::fprintf(stderr, "pool allocation failed (%zu bytes)\n",
            pool_bytes);
        return 1;
    }

    int rc = 0;
    if (cold_mode && rc == 0)
    {
        std::printf("block_bytes,working_blocks,ops_per_repeat,repeats,"
            "median_total_us,median_us_per_xor,median_gib_per_s,checksum\n");
        for (size_t block_bytes : sizes)
        {
            rc = run_cold_xor_size(
                pool, pool_bytes, block_bytes, target_gib, repeats, seed);
            if (rc != 0) {
                break;
            }
        }
    }

    if (muladd_mode && rc == 0)
    {
        std::printf("block_bytes,pool_blocks,ops_per_repeat,repeats,"
            "xor_median_us_per_op,xor_median_gib_per_s,"
            "muladd_median_us_per_op,muladd_median_gib_per_s,"
            "muladd_xor_ratio,checksum\n");
        for (size_t block_bytes : sizes)
        {
            rc = run_muladd_size(
                pool, pool_bytes, block_bytes, target_gib, repeats, seed);
            if (rc != 0) {
                break;
            }
        }
    }

    if (!fanin_ks.empty() && rc == 0)
    {
        std::printf("k,block_bytes,pool_blocks,ops_per_repeat,repeats,"
            "chained_median_us_per_op,chained_logical_gib_per_s,"
            "gather_median_us_per_op,gather_logical_gib_per_s,"
            "gather_chained_ratio,checksum\n");
        for (size_t k : fanin_ks)
        {
            for (size_t block_bytes : sizes)
            {
                rc = run_fanin_size(
                    pool, pool_bytes, k, block_bytes,
                    target_gib, repeats, seed);
                if (rc != 0) {
                    break;
                }
            }
            if (rc != 0) {
                break;
            }
        }
    }

    std::free(pool);
    return rc;
}
