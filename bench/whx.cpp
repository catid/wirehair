// whx.cpp - Wirehair eXperiment harness
// Modes:
//   micro  : per-core + aggregate throughput of gf256 bulk kernels (+ correctness self-check)
//   bench  : end-to-end encoder_create/encode/decode/recover MBPS + overhead across an N grid
//   fuzz   : parallel correctness fuzzer (random N, blockBytes, loss patterns) verifying exact recovery
//
// Designed to saturate a many-core machine (default ~80% of hw threads).
// Built by bench/build.sh, linking the wirehair sources directly.

#include "wirehair/wirehair.h"
#include "WirehairCodec.h"
#include "WirehairEnvironment.h"
#include "WirehairTools.h"
#include "gf256.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cerrno>
#include <climits>
#include <cmath>
#include <vector>
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
#include <algorithm>
#include <chrono>
#include <limits>

using namespace std;
using Clock = std::chrono::steady_clock;

static const char* seed_fixup_policy() {
    return wirehair::SeedFixupsEnabled() ? "enabled" : "disabled";
}

static double now_sec() {
    return std::chrono::duration<double>(Clock::now().time_since_epoch()).count();
}

static bool parse_u64_strict(const char* s, uint64_t& out) {
    if (!s || !*s || *s < '0' || *s > '9') return false;
    errno = 0;
    char* end = nullptr;
    unsigned long long v = strtoull(s, &end, 0);
    if (errno != 0 || end == s || !end || *end != '\0') return false;
    out = (uint64_t)v;
    return true;
}

static bool parse_int_strict(const char* s, int& out) {
    if (!s || !*s) return false;
    errno = 0;
    char* end = nullptr;
    long v = strtol(s, &end, 0);
    if (errno != 0 || end == s || !end || *end != '\0' ||
        v < INT_MIN || v > INT_MAX) return false;
    out = (int)v;
    return true;
}

static bool parse_long_strict(const char* s, long& out) {
    if (!s || !*s) return false;
    errno = 0;
    char* end = nullptr;
    long v = strtol(s, &end, 0);
    if (errno != 0 || end == s || !end || *end != '\0') return false;
    out = v;
    return true;
}

static bool parse_double_strict(const char* s, double& out) {
    if (!s || !*s) return false;
    errno = 0;
    char* end = nullptr;
    double v = strtod(s, &end);
    if (errno != 0 || end == s || !end || *end != '\0' || !std::isfinite(v)) {
        return false;
    }
    out = v;
    return true;
}

static bool valid_roundtrip_args(int N, int bb, int startMode, double loss) {
    return N >= 2 && N <= 64000 && bb >= 1 &&
        (startMode == 0 || startMode == 1) &&
        loss >= 0.0 && loss <= 0.99 && std::isfinite(loss);
}

static bool valid_seedcheck_args(int N, int bb, int startMode, double loss) {
    return N >= 2 && N <= 64000 && bb >= 0 &&
        (startMode == 0 || startMode == 1) &&
        loss >= 0.0 && loss <= 0.99 && std::isfinite(loss);
}

// ---- splitmix64 PRNG (deterministic, fast, per-thread) ----
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() {
        uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    uint32_t u32() { return (uint32_t)(next() >> 32); }
    double unit() { return (next() >> 11) * (1.0 / 9007199254740992.0); } // [0,1)
    // uniform int in [lo, hi]; degenerate ranges return lo instead of
    // dividing by zero (reachable via e.g. "ohead --bb 0")
    int range(int lo, int hi) { return hi <= lo ? lo : lo + (int)(u32() % (uint32_t)(hi - lo + 1)); }
};

static unsigned repair_start_block_id(Rng& rng, int N, unsigned maxNeeded) {
    const uint32_t minId = (uint32_t)N;
    const uint32_t maxStart = UINT32_MAX - maxNeeded;
    if (maxStart <= minId) return minId;
    const uint32_t span = maxStart - minId + 1u;
    return minId + (rng.u32() % span);
}

static void* xaligned(size_t bytes) {
    size_t a = 64;
    size_t rounded = (bytes + a - 1) & ~(a - 1);
    void* p = nullptr;
    if (posix_memalign(&p, a, rounded) != 0) { fprintf(stderr, "OOM %zu\n", bytes); exit(2); }
    return p;
}

static int g_threads = 0;
// Safety cap: never oversubscribe the machine (running >HW threads, or stacking campaigns,
// previously drove the box into a hardware reset). Leave clear headroom.
static int thread_cap_for_hw(int hw) {
    if (hw <= 0) hw = 8;
    int cap = hw - 16;          // leave 16 HW threads free
    if (cap > 104) cap = 104;   // and never exceed 104 regardless
    if (cap < 1) cap = 1;
    return cap;
}
static int default_threads_for_hw(int hw) {
    if (hw <= 0) hw = 8;
    int t = (hw * 3) / 4; // ~75%
    const int cap = thread_cap_for_hw(hw);
    if (t > cap) t = cap;
    if (t < 1) t = 1;
    return t;
}
static int kThreadCap() {
    return thread_cap_for_hw((int)std::thread::hardware_concurrency());
}
static int default_threads() {
    return default_threads_for_hw((int)std::thread::hardware_concurrency());
}
// resolve requested thread count with the hard cap applied
static int resolve_threads() {
    int t = g_threads ? g_threads : default_threads();
    if (t > kThreadCap()) t = kThreadCap();
    if (t < 1) t = 1;
    return t;
}

// ============================================================================
// Reference scalar GF(256) for cross-checking SIMD kernels
// ============================================================================
static uint8_t ref_mul(uint8_t a, uint8_t b) { return gf256_mul(a, b); } // uses scalar table

// Verify the bulk kernels against scalar reference on random data & sizes.
static bool micro_selfcheck() {
    Rng rng(0xC0FFEEULL);
    bool ok = true;

    for (int iter = 0; iter < 2000 && ok; ++iter) {
        int n = rng.range(1, 600);
        // extra padding so we can detect overruns
        const int pad = 16;
        uint8_t* x = (uint8_t*)xaligned(n + pad);
        uint8_t* y = (uint8_t*)xaligned(n + pad);
        uint8_t* z = (uint8_t*)xaligned(n + pad);
        uint8_t* ref = (uint8_t*)xaligned(n + pad);
        for (int i = 0; i < n + pad; ++i) { x[i] = (uint8_t)rng.u32(); y[i] = (uint8_t)rng.u32(); z[i] = (uint8_t)rng.u32(); }
        uint8_t c = (uint8_t)rng.u32();
        // Snapshot the padding; kernels with tail-handling overruns would
        // scribble here, so verify it after every destination write.
        uint8_t padx[pad], pady[pad], padz[pad];
        memcpy(padx, x + n, pad); memcpy(pady, y + n, pad); memcpy(padz, z + n, pad);
        auto pad_ok = [&](const uint8_t* buf, const uint8_t* snap, const char* what) {
            if (memcmp(buf + n, snap, pad) != 0) {
                fprintf(stderr, "FAIL %s overran n=%d\n", what, n);
                return false;
            }
            return true;
        };

        // add_mem: x ^= y
        memcpy(ref, x, n);
        for (int i = 0; i < n; ++i) ref[i] ^= y[i];
        gf256_add_mem(x, y, n);
        if (memcmp(x, ref, n) != 0) { fprintf(stderr, "FAIL add_mem n=%d\n", n); ok = false; }
        if (!pad_ok(x, padx, "add_mem")) ok = false;

        // addset_mem: z = x ^ y
        for (int i = 0; i < n; ++i) ref[i] = x[i] ^ y[i];
        gf256_addset_mem(z, x, y, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL addset_mem n=%d\n", n); ok = false; }
        if (!pad_ok(z, padz, "addset_mem")) ok = false;

        // add2_mem: z ^= x ^ y
        for (int i = 0; i < n; ++i) z[i] = (uint8_t)rng.u32();
        for (int i = 0; i < n; ++i) ref[i] = (uint8_t)(z[i] ^ x[i] ^ y[i]);
        gf256_add2_mem(z, x, y, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL add2_mem n=%d\n", n); ok = false; }
        if (!pad_ok(z, padz, "add2_mem")) ok = false;

        // mul_mem: z = x * c
        for (int i = 0; i < n; ++i) ref[i] = ref_mul(x[i], c);
        gf256_mul_mem(z, x, c, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL mul_mem n=%d c=%u\n", n, c); ok = false; }
        if (!pad_ok(z, padz, "mul_mem")) ok = false;

        // muladd_mem: z += x * c
        for (int i = 0; i < n; ++i) z[i] = (uint8_t)rng.u32();
        for (int i = 0; i < n; ++i) ref[i] = (uint8_t)(z[i] ^ ref_mul(x[i], c));
        gf256_muladd_mem(z, c, x, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL muladd_mem n=%d c=%u\n", n, c); ok = false; }
        if (!pad_ok(z, padz, "muladd_mem")) ok = false;

        // memswap: x <-> y
        memcpy(ref, x, n);
        memcpy(z, y, n);
        gf256_memswap(x, y, n);
        if (memcmp(x, z, n) != 0 || memcmp(y, ref, n) != 0) { fprintf(stderr, "FAIL memswap n=%d\n", n); ok = false; }
        if (!pad_ok(x, padx, "memswap.x") || !pad_ok(y, pady, "memswap.y")) ok = false;

        free(x); free(y); free(z); free(ref);
    }

    // Exhaustive: every constant c in 0..255, over x = 0..255 (+ unaligned offsets),
    // for mul_mem and muladd_mem. Catches a wrong GFNI matrix for any single constant.
    {
        const int n = 256;
        for (int off = 0; off < 7 && ok; ++off) {
            uint8_t* x = (uint8_t*)xaligned(n + 8);
            uint8_t* y = (uint8_t*)xaligned(n + 8);
            uint8_t* z = (uint8_t*)xaligned(n + 8);
            uint8_t* ref = (uint8_t*)xaligned(n + 8);
            uint8_t* yref = (uint8_t*)xaligned(n + 8);
            uint8_t* xb = x + off;
            uint8_t* yb = y + off;
            for (int i = 0; i < n; ++i) xb[i] = (uint8_t)i;
            for (int c = 0; c < 256 && ok; ++c) {
                uint8_t* zb = z + off;
                for (int i = 0; i < n; ++i) ref[i] = ref_mul(xb[i], (uint8_t)c);
                gf256_mul_mem(zb, xb, (uint8_t)c, n);
                if (memcmp(zb, ref, n) != 0) { fprintf(stderr, "FAIL exhaustive mul_mem c=%d off=%d\n", c, off); ok = false; break; }
                for (int i = 0; i < n; ++i) zb[i] = (uint8_t)(i * 7 + 3);
                for (int i = 0; i < n; ++i) ref[i] = (uint8_t)(zb[i] ^ ref_mul(xb[i], (uint8_t)c));
                gf256_muladd_mem(zb, (uint8_t)c, xb, n);
                if (memcmp(zb, ref, n) != 0) { fprintf(stderr, "FAIL exhaustive muladd_mem c=%d off=%d\n", c, off); ok = false; break; }
            }
            for (int i = 0; i < n; ++i) { xb[i] = (uint8_t)(i * 3 + 1); yb[i] = (uint8_t)(i * 5 + 7); }
            memcpy(ref, xb, n);
            memcpy(yref, yb, n);
            gf256_memswap(xb, yb, n);
            if (memcmp(xb, yref, n) != 0 || memcmp(yb, ref, n) != 0) { fprintf(stderr, "FAIL exhaustive memswap off=%d\n", off); ok = false; }
            free(x); free(y); free(z); free(ref); free(yref);
        }
    }
    return ok;
}

// ============================================================================
// micro: throughput
// ============================================================================
enum Kernel { K_ADD, K_ADDSET, K_ADD2, K_MUL, K_MULADD, K_SWAP };
static const char* kname[] = {"add_mem", "addset_mem", "add2_mem", "mul_mem", "muladd_mem", "memswap"};
static const int kKernelCount = 6;

struct MicroResult { double gbps_per_thread; double gbps_aggregate; };

static MicroResult micro_run(Kernel k, int bytes, int threads, double seconds) {
    std::atomic<uint64_t> total_bytes{0};
    std::atomic<int> ready{0};
    std::atomic<bool> go{false};
    std::atomic<bool> stop{false};
    vector<thread> ts;
    vector<double> perthread(threads, 0.0);
    for (int t = 0; t < threads; ++t) {
        ts.emplace_back([&, t]() {
            uint8_t* x = (uint8_t*)xaligned(bytes);
            uint8_t* y = (uint8_t*)xaligned(bytes);
            uint8_t* z = (uint8_t*)xaligned(bytes);
            Rng rng(0x1234ULL + t * 7919);
            for (int i = 0; i < bytes; ++i) { x[i] = (uint8_t)rng.u32(); y[i] = (uint8_t)rng.u32(); z[i] = (uint8_t)rng.u32(); }
            uint8_t c = (uint8_t)(rng.u32() | 2); // avoid 0/1 special-cases
            ready++;
            while (!go.load()) {}
            uint64_t local = 0;
            double t0 = now_sec();
            while (!stop.load()) {
                switch (k) {
                    case K_ADD:    gf256_add_mem(x, y, bytes); break;
                    case K_ADDSET: gf256_addset_mem(z, x, y, bytes); break;
                    case K_ADD2:   gf256_add2_mem(z, x, y, bytes); break;
                    case K_MUL:    gf256_mul_mem(z, x, c, bytes); break;
                    case K_MULADD: gf256_muladd_mem(z, c, x, bytes); break;
                    case K_SWAP:   gf256_memswap(x, y, bytes); break;
                }
                local += bytes;
                // mutate c occasionally so mul tables aren't trivially constant-folded; cheap
                c = (uint8_t)(c + 2);
                if (c < 2) c = 2;
            }
            double dt = now_sec() - t0;
            perthread[t] = local / 1e9 / dt;
            total_bytes += local;
            // prevent dead-code elimination
            volatile uint8_t sink = z[bytes - 1] ^ x[0] ^ y[0];
            (void)sink;
            free(x); free(y); free(z);
        });
    }
    while (ready.load() < threads) {}
    double start = now_sec();
    go.store(true);
    // run for `seconds`
    while (now_sec() - start < seconds) {}
    stop.store(true);
    for (auto& th : ts) th.join();
    double elapsed = now_sec() - start;
    MicroResult r;
    std::sort(perthread.begin(), perthread.end());
    r.gbps_per_thread = perthread[threads / 2]; // median
    r.gbps_aggregate = total_bytes.load() / 1e9 / elapsed;
    return r;
}

static int cmd_micro(int argc, char** argv) {
    int threads = resolve_threads();
    double seconds = 0.7;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--secs") && i + 1 < argc) {
            if (!parse_double_strict(argv[++i], seconds)) return 2;
        }
    }
    if (seconds < 0.0) {
        fprintf(stderr, "micro requires --secs >= 0\n");
        return 2;
    }
    printf("# micro: threads=%d secs=%.2f  (per-thread median GBPS / aggregate GBPS)\n", threads, seconds);
    if (!micro_selfcheck()) { printf("!!! micro_selfcheck FAILED\n"); return 1; }
    printf("# selfcheck OK\n");
    int sizes[] = {1024, 4096, 16384, 65536, 262144, 1048576};
    printf("%-12s", "bytes");
    for (int k = 0; k < kKernelCount; ++k) printf("%22s", kname[k]);
    printf("\n");
    for (int sz : sizes) {
        printf("%-12d", sz);
        for (int k = 0; k < kKernelCount; ++k) {
            MicroResult r = micro_run((Kernel)k, sz, threads, seconds);
            printf("  %8.1f /%8.1f ", r.gbps_per_thread, r.gbps_aggregate);
        }
        printf("\n");
        fflush(stdout);
    }
    return 0;
}

// ============================================================================
// Round-trip core (shared by fuzz + bench)
// ============================================================================
static uint64_t saturating_add_u64(uint64_t a, uint64_t b) {
    return b > std::numeric_limits<uint64_t>::max() - a ?
        std::numeric_limits<uint64_t>::max() : a + b;
}

static uint64_t nearest_rank_index(uint64_t count, long double p) {
    if (count == 0) return 0;
    uint64_t rank;
    if (!(p > 0)) rank = 1;
    else if (p >= 1) rank = count;
    else {
        const long double scaled = p * (long double)count;
        if (scaled >= (long double)count) rank = count;
        else {
            rank = (uint64_t)scaled;
            if ((long double)rank < scaled) ++rank;
        }
        if (rank < 1) rank = 1;
    }
    return rank - 1;
}

// Nearest-rank quantile for a histogram.  Ranks are one-based (ceil(p * n));
// p <= 0 selects the first observation, p >= 1 selects the last, and an empty
// histogram returns zero.  Saturating accumulation prevents a wrapped rank
// if synthetic/test histograms exceed uint64_t's representable total.
static size_t histogram_quantile(const vector<uint64_t>& histogram, long double p) {
    uint64_t total = 0;
    bool overflow = false;
    for (uint64_t count : histogram) {
        if (count > std::numeric_limits<uint64_t>::max() - total) overflow = true;
        total = saturating_add_u64(total, count);
    }
    if (total == 0 || histogram.empty()) return 0;

    // This path is only reachable for synthetic histograms: Stat::trials and
    // its bins share the same uint64_t sample-count bound.  Long double avoids
    // wrap and preserves coherent ordering if a caller constructs larger data.
    if (overflow) {
        long double wide_total = 0;
        for (uint64_t count : histogram) wide_total += (long double)count;
        long double target;
        if (!(p > 0)) target = 1;
        else if (p >= 1) target = wide_total;
        else target = ceill(p * wide_total);
        long double cumulative = 0;
        for (size_t i = 0; i < histogram.size(); ++i) {
            cumulative += (long double)histogram[i];
            if (cumulative >= target) return i;
        }
        return histogram.size() - 1;
    }

    const uint64_t target = nearest_rank_index(total, p) + 1;

    uint64_t cumulative = 0;
    for (size_t i = 0; i < histogram.size(); ++i) {
        cumulative = saturating_add_u64(cumulative, histogram[i]);
        if (cumulative >= target) return i;
    }
    return histogram.size() - 1;
}

struct Stat {
    uint64_t trials = 0;
    uint64_t overhead_sum = 0;       // sum of (needed - N)
    uint64_t overhead_sq = 0;
    uint32_t overhead_max = 0;
    vector<uint64_t> overhead_hist;  // index = extra blocks, capped
    Stat() : overhead_hist(64, 0) {}
    void add_overhead(uint32_t extra) {
        overhead_sum = saturating_add_u64(overhead_sum, extra);
        overhead_sq = saturating_add_u64(overhead_sq, (uint64_t)extra * extra);
        if (extra > overhead_max) overhead_max = extra;
        const size_t bin = extra < 64 ? extra : 63;
        overhead_hist[bin] = saturating_add_u64(overhead_hist[bin], 1);
        trials = saturating_add_u64(trials, 1);
    }
    void merge(const Stat& o) {
        trials = saturating_add_u64(trials, o.trials);
        overhead_sum = saturating_add_u64(overhead_sum, o.overhead_sum);
        overhead_sq = saturating_add_u64(overhead_sq, o.overhead_sq);
        if (o.overhead_max > overhead_max) overhead_max = o.overhead_max;
        for (size_t i = 0; i < overhead_hist.size(); ++i) {
            overhead_hist[i] = saturating_add_u64(overhead_hist[i], o.overhead_hist[i]);
        }
    }
};

// Single round trip. Returns 0 = success, nonzero = failure code. Fills *extra with overhead.
// failmsg (if non-null) gets a description on failure.
static bool matrix_encoder_ok(int N, char* failmsg) {
    if (N < 2 || N > 64000) {
        fprintf(stderr, "matrix_encoder_ok: invalid N=%d (need 2 <= N <= 64000)\n", N);
        exit(2);
    }

    wirehair::Codec enc;
    WirehairResult er = enc.InitializeEncoder((uint64_t)N, 1);
    if (er == Wirehair_Success) {
        er = enc.EncodeFeedMatrixOnly();
    }
    if (er != Wirehair_Success) {
        if (failmsg) sprintf(failmsg, "matrix encoder rc=%d N=%d", er, N);
        return false;
    }

    return true;
}

static int matrix_decode_roundtrip(uint64_t seed, int N, int startMode, double lossRate,
                                   uint32_t* extra_out, char* failmsg) {
    if (N < 2 || N > 64000) {
        fprintf(stderr, "matrix_decode_roundtrip: invalid N=%d (need 2 <= N <= 64000)\n", N);
        exit(2);
    }

    Rng rng(seed);

    wirehair::Codec dec;
    WirehairResult dr = dec.InitializeDecoder((uint64_t)N, 1);
    if (dr != Wirehair_Success) {
        if (failmsg) sprintf(failmsg, "matrix decoder init rc=%d N=%d", dr, N);
        return 2;
    }

    if (lossRate > 0.99) lossRate = 0.99;
    unsigned needed = 0;
    const unsigned maxNeeded = (unsigned)N * 2 + 512;
    unsigned blockId = (startMode == 1) ?
        repair_start_block_id(rng, N, maxNeeded) : 0;

    for (;;) {
        if (needed >= maxNeeded) {
            if (failmsg) sprintf(failmsg, "no matrix decode after %u delivered (N=%d)", needed, N);
            return 3;
        }

        bool drop = (startMode == 0) && (rng.unit() < lossRate);
        unsigned thisId = blockId++;
        if (drop) continue;

        ++needed;
        dr = dec.DecodeFeedMatrixOnly(thisId);
        if (dr == Wirehair_Success) {
            if (extra_out) *extra_out = needed - (uint32_t)N;
            return 0;
        }
        if (dr != Wirehair_NeedMore) {
            if (failmsg) sprintf(failmsg, "matrix decode rc=%d id=%u N=%d", dr, thisId, N);
            return 4;
        }
    }
}

static int matrix_roundtrip(uint64_t seed, int N, int startMode, double lossRate,
                            uint32_t* extra_out, char* failmsg) {
    if (!matrix_encoder_ok(N, failmsg)) {
        return 1;
    }
    return matrix_decode_roundtrip(seed, N, startMode, lossRate, extra_out, failmsg);
}

static int roundtrip(uint64_t seed, int N, int blockBytes, int startMode, double lossRate,
                     uint32_t* extra_out, char* failmsg) {
    // Fail fast with a usage message on out-of-domain parameters.  Most
    // commands pass --nlo/--nmax/--bb through unvalidated; negative or zero
    // values wrap messageBytes to ~2^64 below and the resulting uncaught
    // vector length_error kills a whole campaign with no explanation.
    if (N < 2 || N > 64000 || blockBytes < 0) {
        fprintf(stderr, "roundtrip: invalid N=%d bb=%d (need 2 <= N <= 64000, bb >= 0)\n",
                N, blockBytes);
        exit(2);
    }
    if (blockBytes == 0) {
        return matrix_roundtrip(seed, N, startMode, lossRate, extra_out, failmsg);
    }
    Rng rng(seed);
    int finalBytes = rng.range(1, blockBytes);
    uint64_t messageBytes = (uint64_t)(N - 1) * blockBytes + finalBytes;

    vector<uint8_t> message(messageBytes);
    for (uint64_t i = 0; i < messageBytes; ++i) message[i] = (uint8_t)rng.u32();
    vector<uint8_t> decoded(messageBytes, 0xAB);

    WirehairCodec enc = wirehair_encoder_create(nullptr, message.data(), messageBytes, blockBytes);
    if (!enc) { if (failmsg) sprintf(failmsg, "encoder_create null N=%d bb=%d", N, blockBytes); return 1; }
    WirehairCodec dec = wirehair_decoder_create(nullptr, messageBytes, blockBytes);
    if (!dec) { wirehair_free(enc); if (failmsg) sprintf(failmsg, "decoder_create null N=%d bb=%d", N, blockBytes); return 2; }

    vector<uint8_t> block(blockBytes);
    // Loop only advances on delivered blocks; loss at or near 1.0 would spin forever.
    if (lossRate > 0.99) lossRate = 0.99;
    // startMode 0: from blockId 0 (systematic + repair). 1: repair-only from a random offset >= N.
    unsigned needed = 0;
    int rc = 0;
    // Cap on *delivered* blocks. Generous so rare-but-legitimate wirehair overhead tails
    // (pathological erasure subsets, ~1e-5 freq) still complete and get counted. A genuinely
    // broken ablation surfaces as decode-error/mismatch (caught immediately) or never-decode.
    const unsigned maxNeeded = (unsigned)N * 2 + 512;
    unsigned blockId = (startMode == 1) ?
        repair_start_block_id(rng, N, maxNeeded) : 0;
    bool decoded_ok = false;
    for (;;) {
        if (needed >= maxNeeded) {
            if (failmsg) sprintf(failmsg, "no decode after %u delivered (N=%d bb=%d)", needed, N, blockBytes);
            rc = 3; break;
        }
        // drop?
        bool drop = (startMode == 0) && (rng.unit() < lossRate);
        unsigned thisId = blockId++;
        if (drop) continue;

        uint32_t writeLen = 0;
        WirehairResult er = wirehair_encode(enc, thisId, block.data(), blockBytes, &writeLen);
        if (er != Wirehair_Success) { if (failmsg) sprintf(failmsg, "encode rc=%d id=%u N=%d bb=%d", er, thisId, N, blockBytes); rc = 4; break; }
        needed++;
        WirehairResult dr = wirehair_decode(dec, thisId, block.data(), writeLen);
        if (dr == Wirehair_Success) { decoded_ok = true; break; }
        if (dr != Wirehair_NeedMore) { if (failmsg) sprintf(failmsg, "decode rc=%d id=%u N=%d bb=%d", dr, thisId, N, blockBytes); rc = 5; break; }
    }

    if (decoded_ok) {
        WirehairResult rr = wirehair_recover(dec, decoded.data(), messageBytes);
        if (rr != Wirehair_Success) { if (failmsg) sprintf(failmsg, "recover rc=%d N=%d bb=%d", rr, N, blockBytes); rc = 6; }
        else if (memcmp(decoded.data(), message.data(), messageBytes) != 0) {
            // find first mismatch
            uint64_t i = 0; while (i < messageBytes && decoded[i] == message[i]) ++i;
            if (failmsg) sprintf(failmsg, "MISMATCH at %llu N=%d bb=%d final=%d needed=%u", (unsigned long long)i, N, blockBytes, finalBytes, needed);
            rc = 7;
        } else {
            if (extra_out) *extra_out = needed - (uint32_t)N;
        }
    }
    wirehair_free(enc);
    wirehair_free(dec);
    return rc;
}

// ============================================================================
// fuzz
// ============================================================================
static std::atomic<uint64_t> g_fuzz_trials{0};
static std::atomic<uint64_t> g_fuzz_fail{0};
static std::mutex g_print_mu;

static int sample_N(Rng& rng, int Nmax) {
    // log-uniform with emphasis on small/medium where structure switches
    double u = rng.unit();
    double lo = log(2.0), hi = log((double)Nmax);
    int N = (int)exp(lo + u * (hi - lo));
    if (N < 2) N = 2;
    if (N > Nmax) N = Nmax;
    return N;
}
static int sample_blockBytes(Rng& rng) {
    // mostly small (exercise partial final block & tails), some larger
    int pick = rng.range(0, 9);
    if (pick == 0) return rng.range(1, 8);
    if (pick <= 3) return rng.range(1, 64);
    if (pick <= 6) return rng.range(1, 600);
    if (pick <= 8) return rng.range(1, 2048);
    return rng.range(1, 8192);
}

static bool parse_positive_int_list(const string& list, vector<int>& values) {
    values.clear();
    if (list.empty() || list[list.size() - 1] == ',') {
        return false;
    }
    size_t p = 0;
    while (p < list.size()) {
        size_t q = list.find(',', p);
        string tok = list.substr(p, q == string::npos ? string::npos : q - p);
        int v = 0;
        if (tok.empty() || tok[0] < '0' || tok[0] > '9' ||
            !parse_int_strict(tok.c_str(), v) || v <= 0)
        {
            return false;
        }
        values.push_back(v);
        if (q == string::npos) break;
        p = q + 1;
    }
    return !values.empty();
}

static bool read_n_file_strict(const string& path, const char* command, vector<int>& values) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) {
        fprintf(stderr, "%s: cannot open --nfile %s\n", command, path.c_str());
        return false;
    }
    char linebuf[512];
    int lineno = 0;
    while (fgets(linebuf, sizeof(linebuf), f)) {
        ++lineno;
        const size_t length = strlen(linebuf);
        if (length == sizeof(linebuf) - 1 && linebuf[length - 1] != '\n' && !feof(f)) {
            fprintf(stderr, "%s: overlong line at %s:%d\n", command, path.c_str(), lineno);
            fclose(f);
            return false;
        }
        const char* s = linebuf;
        while (*s == ' ' || *s == '\t') ++s;
        if (*s == '\0' || *s == '\n' || *s == '#') continue;
        errno = 0;
        char* end = nullptr;
        long n = strtol(s, &end, 10);
        while (*end == ' ' || *end == '\t' || *end == '\r' || *end == '\n') ++end;
        if (errno != 0 || end == s || *end != '\0' || n < 2 || n > 64000) {
            fprintf(stderr, "%s: bad N at %s:%d: %s", command, path.c_str(), lineno, linebuf);
            if (length == 0 || linebuf[length - 1] != '\n') fputc('\n', stderr);
            fclose(f);
            return false;
        }
        values.push_back((int)n);
    }
    if (ferror(f)) {
        fprintf(stderr, "%s: read failed for --nfile %s\n", command, path.c_str());
        fclose(f);
        return false;
    }
    fclose(f);
    if (values.empty()) {
        fprintf(stderr, "%s: --nfile %s has no N values\n", command, path.c_str());
        return false;
    }
    return true;
}

static bool validate_n_bb_lists(const vector<int>& Ns, const vector<int>& BBs, const char* cmd) {
    if (Ns.empty() || BBs.empty()) {
        fprintf(stderr, "%s requires non-empty --N and --bb/--bb-list positive integer lists\n", cmd);
        return false;
    }
    for (int N : Ns) {
        if (N < 2 || N > 64000) {
            fprintf(stderr, "%s requires every N in 2..64000\n", cmd);
            return false;
        }
    }
    for (int bb : BBs) {
        if (bb < 1) {
            fprintf(stderr, "%s requires every block byte count to be positive\n", cmd);
            return false;
        }
    }
    return true;
}

static bool option_in(const char* arg, const char* const* options) {
    for (int i = 0; options[i]; ++i) {
        if (!strcmp(arg, options[i])) return true;
    }
    return false;
}

static bool validate_mode_options(const string& mode, int argc, char** argv) {
    static const char* none[] = { nullptr };
    static const char* micro_values[] = { "--secs", nullptr };
    static const char* bench_values[] = { "--bb", "--bb-list", "--loss", "--N", "--rounds", "--memory-mib", nullptr };
    static const char* fuzz_values[] = { "--secs", "--nmax", "--seed", nullptr };
    static const char* ohead_values[] = { "--nlo", "--nhi", "--nstep", "--nfile", "--trials", "--bb", "--startmode", "--loss", "--seed", "--samples-out", nullptr };
    static const char* scan_values[] = { "--nlo", "--nhi", "--nfile", "--trials", "--bb", "--startmode", "--loss", "--thresh", "--seed", nullptr };
#ifdef WH_SEED_KNOBS
    static const char* seedmean_values[] = { "--N", "--n", "--pseed", "--dseed", "--trials", "--bb", "--loss", "--startmode", "--seed", nullptr };
    static const char* seedsearch_values[] = { "--nlist", "--nfile", "--tsearch", "--tverify", "--nseeds", "--bb", "--loss", "--startmode", "--dseeds", "--goodthr", nullptr };
#endif
#ifdef WH_COUNT
    static const char* peelstat_values[] = { "--N", "--bb", "--bb-list", "--trials", "--loss", "--startmode", "--seed", "--threads", nullptr };
    static const char* count_values[] = { "--N", "--bb", "--bb-list", "--loss", "--startmode", nullptr };
#endif

    const char* const* value_options = none;
    const char* const* flag_options = none;
    if (mode == "micro") value_options = micro_values;
    else if (mode == "bench") value_options = bench_values;
    else if (mode == "fuzz") value_options = fuzz_values;
    else if (mode == "repro") value_options = none;
    else if (mode == "ohead") value_options = ohead_values;
    else if (mode == "scan") value_options = scan_values;
#ifdef WH_SEED_KNOBS
    else if (mode == "seedmean") value_options = seedmean_values;
    else if (mode == "seedsearch") value_options = seedsearch_values;
#endif
#ifdef WH_COUNT
    else if (mode == "peelstat") value_options = peelstat_values;
    else if (mode == "count") value_options = count_values;
#endif
    else {
        return true; // Unknown mode is reported by the dispatcher.
    }

    for (int i = 0; i < argc; ++i) {
        const char* arg = argv[i];
        if (option_in(arg, value_options)) {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s requires a value\n", arg);
                return false;
            }
            ++i;
        }
        else if (option_in(arg, flag_options)) {
            continue;
        }
        else if (arg[0] == '-') {
            fprintf(stderr, "unknown option for %s: %s\n", mode.c_str(), arg);
            return false;
        }
    }
    return true;
}

static int cmd_fuzz(int argc, char** argv) {
    int threads = resolve_threads();
    double seconds = 20.0;
    int Nmax = 4000;
    uint64_t baseSeed = 0xABCDEF12345ULL;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--secs") && i + 1 < argc) { if (!parse_double_strict(argv[++i], seconds)) return 2; }
        else if (!strcmp(argv[i], "--nmax") && i + 1 < argc) { if (!parse_int_strict(argv[++i], Nmax)) return 2; }
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) { if (!parse_u64_strict(argv[++i], baseSeed)) return 2; }
    }
    if (seconds < 0.0 || Nmax < 2 || Nmax > 64000) {
        fprintf(stderr, "fuzz requires --secs >= 0 and 2 <= --nmax <= 64000\n");
        return 2;
    }
    printf("# fuzz: threads=%d secs=%.1f Nmax=%d seed=0x%llx\n", threads, seconds, Nmax, (unsigned long long)baseSeed);
    Stat global; std::mutex gmu;
    double start = now_sec();
    vector<thread> ts;
    for (int t = 0; t < threads; ++t) {
        ts.emplace_back([&, t]() {
            Rng rng(baseSeed + 0x1000003ULL * (t + 1));
            Stat local;
            char failmsg[256];
            while (now_sec() - start < seconds) {
                for (int b = 0; b < 32; ++b) {
                    uint64_t tseed = rng.next();
                    int N = sample_N(rng, Nmax);
                    int bb = sample_blockBytes(rng);
                    int startMode = rng.range(0, 1);
                    double loss = 0.02 + rng.unit() * 0.43; // up to ~45% loss (realistic upper bound)
                    uint32_t extra = 0;
                    failmsg[0] = 0;
                    int rc = roundtrip(tseed, N, bb, startMode, loss, &extra, failmsg);
                    if (rc == 0) {
                        local.add_overhead(extra);
                        if (extra >= 64) {
                            std::lock_guard<std::mutex> lk(g_print_mu);
                            printf("HIGH-OVERHEAD extra=%u seed=0x%llx N=%d bb=%d startMode=%d loss=%.2f\n",
                                   extra, (unsigned long long)tseed, N, bb, startMode, loss);
                        }
                    }
                    else {
                        g_fuzz_fail++;
                        std::lock_guard<std::mutex> lk(g_print_mu);
                        printf("FAIL rc=%d seed=0x%llx N=%d bb=%d startMode=%d loss=%.2f : %s\n",
                               rc, (unsigned long long)tseed, N, bb, startMode, loss, failmsg);
                        fflush(stdout);
                    }
                    g_fuzz_trials++;
                }
            }
            std::lock_guard<std::mutex> lk(gmu);
            global.merge(local);
        });
    }
    for (auto& th : ts) th.join();
    double dt = now_sec() - start;
    uint64_t trials = g_fuzz_trials.load(), fails = g_fuzz_fail.load();
    double mean = global.trials ? (double)global.overhead_sum / global.trials : 0.0;
    printf("# fuzz done: %.1fs trials=%llu fails=%llu rate=%.0f trials/s\n",
           dt, (unsigned long long)trials, (unsigned long long)fails, trials / dt);
    printf("# overhead(success trials=%llu): mean=%.4f p50=%d p99=%d p999=%d max=%u\n",
           (unsigned long long)global.trials, mean,
           (int)histogram_quantile(global.overhead_hist, 0.50L),
           (int)histogram_quantile(global.overhead_hist, 0.99L),
           (int)histogram_quantile(global.overhead_hist, 0.999L),
           global.overhead_max);
    return fails ? 1 : 0;
}

// ============================================================================
// bench (end-to-end), parallel across threads, fixed work per (N)
// ============================================================================
static bool checked_mul_u64(uint64_t a, uint64_t b, uint64_t& out) {
    if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a) return false;
    out = a * b;
    return true;
}

static bool checked_add_u64(uint64_t a, uint64_t b, uint64_t& out) {
    if (b > std::numeric_limits<uint64_t>::max() - a) return false;
    out = a + b;
    return true;
}

static bool checked_add_product_u64(uint64_t& total, uint64_t a, uint64_t b) {
    uint64_t product = 0;
    return checked_mul_u64(a, b, product) &&
        checked_add_u64(total, product, total);
}

struct BenchWorkerPlan {
    int active_workers = 0;
    uint64_t bytes_per_worker = 0;
};

// Bound one encoder, one decoder, and the three explicit benchmark buffers.
// Codec matrices depend on peeling, so use the largest legal deferred/dense
// dimensions.  This intentionally overestimates normal runs but makes the
// configured budget a real upper bound rather than an N*blockBytes proxy that
// severely undercounts small-block workloads.
static bool estimate_bench_worker_bytes(uint64_t N, uint64_t block_bytes,
                                        uint64_t& bytes) {
    const uint64_t max_mix = CAT_MAX_DENSE_ROWS + wirehair::kHeavyRows;

    auto add_codec = [&](uint64_t extra_rows, bool decoder) -> bool {
        uint64_t recovery_rows = 0;
        if (!checked_add_u64(N, max_mix + 1, recovery_rows)) return false;
        uint64_t recovery_bytes = 0;
        if (!checked_mul_u64(recovery_rows, block_bytes, recovery_bytes)) return false;
        uint64_t peel_offset = 0;
        if (!checked_add_u64(recovery_bytes, 7, peel_offset)) return false;
        peel_offset &= ~UINT64_C(7);
        if (!checked_add_u64(bytes, peel_offset, bytes)) return false;

        uint64_t row_count = 0;
        if (!checked_add_u64(N, extra_rows, row_count) ||
            !checked_add_product_u64(bytes, row_count, sizeof(wirehair::PeelRow)) ||
            !checked_add_product_u64(bytes, N, sizeof(wirehair::PeelColumn)) ||
            !checked_add_product_u64(bytes, N, sizeof(wirehair::PeelRefs)) ||
            !checked_add_product_u64(bytes, N, sizeof(uint16_t)) ||
            !checked_add_u64(bytes, row_count, bytes))
        {
            return false;
        }

        uint64_t ge_cols = 0;
        if (!checked_add_u64(N, max_mix, ge_cols)) return false;
        uint64_t rounded_cols = 0;
        if (!checked_add_u64(ge_cols, 63, rounded_cols)) return false;
        uint64_t ge_pitch = rounded_cols / 64;
#if defined(WH_ALIGN64) && (WH_ALIGN64+0)
        uint64_t aligned_pitch = 0;
        if (!checked_add_u64(ge_pitch, 7, aligned_pitch)) return false;
        ge_pitch = aligned_pitch & ~UINT64_C(7);
#endif
        uint64_t ge_rows = 0;
        if (!checked_add_u64(N, CAT_MAX_DENSE_ROWS + extra_rows + 1, ge_rows)) {
            return false;
        }
        uint64_t matrix_rows = 0;
        uint64_t matrix_words = 0;
        if (!checked_add_u64(N, ge_rows, matrix_rows) ||
            !checked_mul_u64(matrix_rows, ge_pitch, matrix_words)) return false;
        if (!checked_add_product_u64(bytes, matrix_words, sizeof(uint64_t))) return false;

        const uint64_t heavy_cols = wirehair::kHeavyCols;
        const uint64_t heavy_pitch = (heavy_cols + 6) & ~UINT64_C(3);
        if (!checked_add_product_u64(
                bytes, heavy_pitch, wirehair::kHeavyRows + extra_rows))
        {
            return false;
        }
        uint64_t pivot_count = 0;
        if (!checked_add_u64(ge_cols, extra_rows, pivot_count)) return false;
        uint64_t pivot_words = 0;
        if (!checked_mul_u64(pivot_count, 2, pivot_words) ||
            !checked_add_u64(pivot_words, ge_cols, pivot_words) ||
            !checked_add_product_u64(bytes, pivot_words, sizeof(uint16_t)))
        {
            return false;
        }

        if (decoder &&
            !checked_add_product_u64(bytes, row_count, block_bytes))
        {
            return false;
        }
        return checked_add_u64(bytes, sizeof(wirehair::Codec), bytes);
    };

    bytes = 0;
    if (!checked_add_product_u64(bytes, N, block_bytes) ||
        !checked_add_product_u64(bytes, N, block_bytes) ||
        !checked_add_u64(bytes, block_bytes, bytes) ||
        !add_codec(0, false) ||
        !add_codec(CAT_MAX_EXTRA_ROWS, true))
    {
        return false;
    }
    return bytes != 0;
}

static bool plan_bench_workers(int requested_threads, uint64_t total_rounds,
                               uint64_t N, uint64_t block_bytes,
                               uint64_t memory_budget, BenchWorkerPlan& plan) {
    uint64_t per_worker = 0;
    if (requested_threads < 1 || total_rounds == 0 || N == 0 || block_bytes == 0 ||
        !estimate_bench_worker_bytes(N, block_bytes, per_worker))
    {
        return false;
    }
    const uint64_t budget_workers = memory_budget / per_worker;
    if (budget_workers == 0) return false;
    uint64_t active = (uint64_t)requested_threads;
    if (active > total_rounds) active = total_rounds;
    if (active > budget_workers) active = budget_workers;
    if (active == 0 || active > (uint64_t)std::numeric_limits<int>::max()) return false;
    plan.active_workers = (int)active;
    plan.bytes_per_worker = per_worker;
    return true;
}

static bool verify_recovered_bytes(const uint8_t* expected, const uint8_t* actual,
                                   uint64_t bytes, uint64_t& mismatch_offset) {
    if (bytes == 0 || memcmp(expected, actual, (size_t)bytes) == 0) {
        mismatch_offset = bytes;
        return true;
    }
    uint64_t i = 0;
    while (i < bytes && expected[i] == actual[i]) ++i;
    mismatch_offset = i;
    return false;
}

enum BenchFaultLocation { BENCH_FAULT_NONE, BENCH_FAULT_FIRST, BENCH_FAULT_MIDDLE,
                          BENCH_FAULT_LAST, BENCH_FAULT_ABSOLUTE };
static BenchFaultLocation g_bench_fault_location = BENCH_FAULT_NONE;
static uint64_t g_bench_fault_offset = 0;

static bool configure_bench_fault() {
    const wirehair::EnvironmentValue environment("WHX_BENCH_CORRUPT");
    const char* value = environment.Get();
    if (!value || !*value) return true;
    if (!strcmp(value, "first")) g_bench_fault_location = BENCH_FAULT_FIRST;
    else if (!strcmp(value, "middle")) g_bench_fault_location = BENCH_FAULT_MIDDLE;
    else if (!strcmp(value, "last")) g_bench_fault_location = BENCH_FAULT_LAST;
    else {
        if (!parse_u64_strict(value, g_bench_fault_offset)) {
            fprintf(stderr, "WHX_BENCH_CORRUPT must be first, middle, last, or a byte offset\n");
            return false;
        }
        g_bench_fault_location = BENCH_FAULT_ABSOLUTE;
    }
    return true;
}

static bool bench_fault_offset(uint64_t bytes, uint64_t& offset) {
    if (g_bench_fault_location == BENCH_FAULT_NONE || bytes == 0) return false;
    if (g_bench_fault_location == BENCH_FAULT_FIRST) offset = 0;
    else if (g_bench_fault_location == BENCH_FAULT_MIDDLE) offset = bytes / 2;
    else if (g_bench_fault_location == BENCH_FAULT_LAST) offset = bytes - 1;
    else offset = g_bench_fault_offset;
    return offset < bytes;
}

struct BenchAccum {
    double enc_create_us = 0; uint64_t enc_create_n = 0;
    double encode_us = 0; uint64_t encode_n = 0; uint64_t encode_bytes = 0;
    double decode_us = 0; uint64_t decode_n = 0; uint64_t decode_bytes = 0;
    double recover_us = 0; uint64_t recover_n = 0; uint64_t recover_bytes = 0;
    uint64_t overhead_sum = 0; uint64_t rounds = 0; uint64_t fails = 0;
    uint64_t mismatches = 0;
    void merge(const BenchAccum& o){
        enc_create_us+=o.enc_create_us; enc_create_n+=o.enc_create_n;
        encode_us+=o.encode_us; encode_n+=o.encode_n; encode_bytes+=o.encode_bytes;
        decode_us+=o.decode_us; decode_n+=o.decode_n; decode_bytes+=o.decode_bytes;
        recover_us+=o.recover_us; recover_n+=o.recover_n; recover_bytes+=o.recover_bytes;
        overhead_sum+=o.overhead_sum; rounds+=o.rounds; fails+=o.fails;
        mismatches+=o.mismatches;
    }
};

// Persistent per-thread state so codec objects + buffers (and their page faults) are
// allocated ONCE and reused across rounds -- otherwise allocation/faulting dominates and
// masks codec compute. wirehair's reuseOpt parameter recycles the prior codec allocation.
struct BenchWorker {
    WirehairCodec enc = nullptr, dec = nullptr;
    vector<uint8_t> message, decoded, block;
    ~BenchWorker() { if (enc) wirehair_free(enc); if (dec) wirehair_free(dec); }
};

static bool bench_one_round(BenchWorker& w, int N, int blockBytes, double lossRate, uint64_t seed,
                            BenchAccum& a, bool timed) {
    Rng rng(seed);
    uint64_t messageBytes = (uint64_t)N * blockBytes; // full blocks for clean throughput accounting
    // light per-round mutation so the solve isn't trivially identical, cheap
    for (uint64_t i = 0; i < messageBytes; i += 64) w.message[i] = (uint8_t)rng.u32();

    double t0 = now_sec();
    w.enc = wirehair_encoder_create(w.enc, w.message.data(), messageBytes, blockBytes); // reuse
    double t1 = now_sec();
    if (!w.enc) { if (timed) a.fails++; return false; }
    if (timed) { a.enc_create_us += (t1 - t0) * 1e6; a.enc_create_n++; }

    w.dec = wirehair_decoder_create(w.dec, messageBytes, blockBytes); // reuse
    if (!w.dec) { if (timed) a.fails++; return false; }

    // measure encode of repair blocks (ids N..N+enc_count-1)
    const int enc_count = 64;
    uint32_t wl = 0;
    double e0 = now_sec();
    for (int i = 0; i < enc_count; ++i) {
        if (wirehair_encode(w.enc, N + i, w.block.data(), blockBytes, &wl) !=
            Wirehair_Success)
        {
            if (timed) a.fails++;
            return false;
        }
    }
    double e1 = now_sec();
    if (timed) { a.encode_us += (e1 - e0) * 1e6; a.encode_n += enc_count; a.encode_bytes += (uint64_t)enc_count * blockBytes; }

    // decode with random losses; time decode feed
    unsigned blockId = 0, needed = 0;
    double d_accum = 0;
    bool ok = false;
    const unsigned maxNeeded = (unsigned)N * 2 + 512;
    for (;;) {
        if (needed >= maxNeeded) break;
        bool drop = rng.unit() < lossRate;
        unsigned thisId = blockId++;
        if (drop) continue;
        if (wirehair_encode(w.enc, thisId, w.block.data(), blockBytes, &wl) !=
            Wirehair_Success)
        {
            break;
        }
        needed++;
        double d0 = now_sec();
        WirehairResult dr = wirehair_decode(w.dec, thisId, w.block.data(), wl);
        d_accum += now_sec() - d0;
        if (dr == Wirehair_Success) { ok = true; break; }
        if (dr != Wirehair_NeedMore) break;
    }
    if (timed) {
        a.decode_us += d_accum * 1e6; a.decode_n += needed; a.decode_bytes += (uint64_t)needed * blockBytes;
    }
    if (ok) {
        double r0 = now_sec();
        WirehairResult rr = wirehair_recover(w.dec, w.decoded.data(), messageBytes);
        double r1 = now_sec();
        if (rr == Wirehair_Success) {
            uint64_t injected_offset = 0;
            if (timed && bench_fault_offset(messageBytes, injected_offset)) {
                w.decoded[(size_t)injected_offset] ^= 1;
            }
            uint64_t mismatch_offset = 0;
            const bool exact = verify_recovered_bytes(
                w.message.data(), w.decoded.data(), messageBytes, mismatch_offset);
            if (!exact) {
                if (timed) {
                    a.fails++;
                    a.mismatches++;
                }
                std::lock_guard<std::mutex> lk(g_print_mu);
                fprintf(stderr,
                        "BENCH_MISMATCH seed=0x%llx N=%d bb=%d bytes=%llu offset=%llu timed=%d\n",
                        (unsigned long long)seed, N, blockBytes,
                        (unsigned long long)messageBytes,
                        (unsigned long long)mismatch_offset, timed ? 1 : 0);
                return false;
            }
            if (timed) {
                // Exact comparison is intentionally after r1, outside recover timing.
                a.recover_us += (r1 - r0) * 1e6;
                a.recover_n++;
                a.recover_bytes += messageBytes;
                a.overhead_sum += (needed >= (unsigned)N) ? (needed - N) : 0;
                a.rounds++;
            }
            return true;
        }
        if (timed) {
            a.fails++;
        }
        return false;
    }
    else if (timed) {
        // Failed rounds are counted separately; folding their capped "needed"
        // into overhead_sum would silently shift the overhead statistic.
        a.fails++;
    }
    return false;
}

static int cmd_bench(int argc, char** argv) {
    int threads = resolve_threads();
    string bblist = "1300";
    double lossRate = 0.10;
    uint64_t memory_mib = 4096;
    int rounds_per_N = 0; // 0 => auto by N
    bool rounds_set = false;
    string nlist = "32,128,512,1024,2048,8192,32000";
    const wirehair::EnvironmentValue memory_environment(
        "WHX_BENCH_MEMORY_MIB");
    const char* memory_env = memory_environment.Get();
    if (memory_env && !parse_u64_strict(memory_env, memory_mib)) {
        fprintf(stderr, "invalid WHX_BENCH_MEMORY_MIB\n");
        return 2;
    }
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--bb") && i + 1 < argc) bblist = argv[++i];
        else if (!strcmp(argv[i], "--bb-list") && i + 1 < argc) bblist = argv[++i];
        else if (!strcmp(argv[i], "--loss") && i + 1 < argc) {
            if (!parse_double_strict(argv[++i], lossRate)) return 2;
        }
        else if (!strcmp(argv[i], "--N") && i + 1 < argc) nlist = argv[++i];
        else if (!strcmp(argv[i], "--rounds") && i + 1 < argc) {
            if (!parse_int_strict(argv[++i], rounds_per_N)) return 2;
            rounds_set = true;
        }
        else if (!strcmp(argv[i], "--memory-mib") && i + 1 < argc) {
            if (!parse_u64_strict(argv[++i], memory_mib)) return 2;
        }
    }
    vector<int> Ns;
    vector<int> BBs;
    if (!parse_positive_int_list(nlist, Ns) ||
        !parse_positive_int_list(bblist, BBs) ||
        !validate_n_bb_lists(Ns, BBs, "bench"))
    {
        return 2;
    }
    uint64_t memory_budget = 0;
    if (!std::isfinite(lossRate) || lossRate < 0.0 || lossRate >= 1.0 ||
        memory_mib == 0 || !checked_mul_u64(memory_mib, 1024ull * 1024ull, memory_budget) ||
        (rounds_set && rounds_per_N < 1))
    {
        fprintf(stderr, "bench requires 0 <= --loss < 1, positive --rounds when specified, and a representable positive --memory-mib\n");
        return 2;
    }
    if (!configure_bench_fault()) return 2;
    if (g_bench_fault_location == BENCH_FAULT_ABSOLUTE) {
        for (int blockBytes : BBs) for (int N : Ns) {
            const uint64_t message_bytes = (uint64_t)N * (uint64_t)blockBytes;
            if (g_bench_fault_offset >= message_bytes) {
                fprintf(stderr,
                        "WHX_BENCH_CORRUPT offset %llu is outside N=%d bb=%d message (%llu bytes)\n",
                        (unsigned long long)g_bench_fault_offset, N, blockBytes,
                        (unsigned long long)message_bytes);
                return 2;
            }
        }
    }
    const bool matrix = BBs.size() > 1;
    printf("# bench: requested_threads=%d bb=%s loss=%.2f memory_mib=%llu verify=exact(outside_timing)\n",
           threads, bblist.c_str(), lossRate, (unsigned long long)memory_mib);
    if (matrix)
        printf("%-8s %-8s %10s %14s %14s %14s %14s %10s\n", "N", "bb", "msg_MiB", "create_MBPS", "encode_MBPS", "decode_MBPS", "recover_MBPS", "overhead");
    else
        printf("%-8s %14s %14s %14s %14s %10s\n", "N", "create_MBPS", "encode_MBPS", "decode_MBPS", "recover_MBPS", "overhead");
    bool failed = false;
    for (int blockBytes : BBs) for (int N : Ns) {
        // total rounds across all threads; scale down for large N to keep runtime bounded
        int total_rounds = rounds_per_N;
        if (total_rounds == 0) {
            if (N <= 128) total_rounds = threads * 40;
            else if (N <= 1024) total_rounds = threads * 16;
            else if (N <= 8192) total_rounds = threads * 4;
            else total_rounds = threads * 2;
            const uint64_t msgBytes = (uint64_t)N * blockBytes;
            const uint64_t roundCapBytes = 64ull * 1024 * 1024;
            if (msgBytes > roundCapBytes && total_rounds > threads) {
                uint64_t scaled = ((uint64_t)total_rounds * roundCapBytes) / msgBytes;
                if (scaled < (uint64_t)threads) scaled = threads;
                total_rounds = (int)scaled;
            }
        }
        BenchWorkerPlan plan;
        if (!plan_bench_workers(threads, (uint64_t)total_rounds, (uint64_t)N,
                                (uint64_t)blockBytes, memory_budget, plan))
        {
            fprintf(stderr,
                    "bench memory budget cannot support one worker: N=%d bb=%d rounds=%d memory_mib=%llu\n",
                    N, blockBytes, total_rounds, (unsigned long long)memory_mib);
            return 2;
        }
        std::atomic<uint64_t> next{0};
        BenchAccum global; std::mutex mu;
        vector<thread> ts;
        for (int t = 0; t < plan.active_workers; ++t) {
            ts.emplace_back([&, t]() {
                BenchWorker w;
                uint64_t mb = (uint64_t)N * blockBytes;
                w.message.assign(mb, (uint8_t)(0x33 + t));
                w.decoded.assign(mb, 0);
                w.block.assign(blockBytes, 0);
                BenchAccum a;
                // warmup (untimed): allocate codecs + pre-fault all pages
                if (!bench_one_round(w, N, blockBytes, lossRate,
                                     0x9999ULL + (uint64_t)N * 131 + t, a, false))
                {
                    a.fails++;
                    std::lock_guard<std::mutex> lk(mu); global.merge(a);
                    return;
                }
                for (;;) {
                    const uint64_t idx = next.fetch_add(1);
                    if (idx >= (uint64_t)total_rounds) break;
                    bench_one_round(w, N, blockBytes, lossRate, 0x9999ULL + (uint64_t)N * 131 + idx * 2654435761ULL, a, true);
                }
                std::lock_guard<std::mutex> lk(mu); global.merge(a);
            });
        }
        for (auto& th : ts) th.join();
        auto mbps = [](uint64_t bytes, double us){ return us > 0 ? bytes / us : 0.0; };
        double cr = global.enc_create_n ? mbps((uint64_t)global.enc_create_n * N * blockBytes / 1, 0) : 0;
        // create MBPS = total message bytes processed / total create time
        double create_MBPS = global.enc_create_us > 0 ? ((double)global.enc_create_n * N * blockBytes) / global.enc_create_us : 0;
        double encode_MBPS = global.encode_us > 0 ? (double)global.encode_bytes / global.encode_us : 0;
        double decode_MBPS = global.decode_us > 0 ? (double)global.decode_bytes / global.decode_us : 0;
        double recover_MBPS = global.recover_us > 0 ? (double)global.recover_bytes / global.recover_us : 0;
        double overhead = global.rounds ? (double)global.overhead_sum / global.rounds : 0;
        (void)cr;
        if (matrix) {
            double msgMiB = (double)((uint64_t)N * blockBytes) / 1048576.0;
            printf("%-8d %-8d %10.2f %14.1f %14.1f %14.1f %14.1f %10.4f",
                   N, blockBytes, msgMiB, create_MBPS, encode_MBPS, decode_MBPS, recover_MBPS, overhead);
        } else {
            printf("%-8d %14.1f %14.1f %14.1f %14.1f %10.4f",
                   N, create_MBPS, encode_MBPS, decode_MBPS, recover_MBPS, overhead);
        }
        if (global.fails) printf("  FAILS=%llu", (unsigned long long)global.fails);
        if (global.mismatches) printf(" MISMATCHES=%llu", (unsigned long long)global.mismatches);
        printf("  workers=%d worker_MiB=%.2f verified=%llu",
               plan.active_workers, (double)plan.bytes_per_worker / 1048576.0,
               (unsigned long long)global.rounds);
        if (global.fails || global.rounds == 0) {
            failed = true;
        }
        printf("\n");
        fflush(stdout);
    }
    return failed ? 1 : 0;
}

// Deterministic single-trial trace: whx repro <seed> <N> <bb> <startMode> <loss>
static int cmd_repro(int argc, char** argv) {
    if (argc < 5) { fprintf(stderr, "repro <seed> <N> <bb> <startMode> <loss>\n"); return 1; }
    uint64_t seed = 0;
    int N = 0, bb = 0, startMode = 0;
    double loss = 0.0;
    if (!parse_u64_strict(argv[0], seed) ||
        !parse_int_strict(argv[1], N) ||
        !parse_int_strict(argv[2], bb) ||
        !parse_int_strict(argv[3], startMode) ||
        !parse_double_strict(argv[4], loss) ||
        !valid_roundtrip_args(N, bb, startMode, loss))
    {
        fprintf(stderr, "repro: invalid arguments (need 2 <= N <= 64000, bb >= 1, startMode 0/1, 0 <= loss <= 0.99)\n");
        return 2;
    }
    Rng rng(seed);
    int finalBytes = rng.range(1, bb);
    uint64_t messageBytes = (uint64_t)(N - 1) * bb + finalBytes;
    vector<uint8_t> message(messageBytes), decoded(messageBytes, 0xAB), block(bb);
    for (uint64_t i = 0; i < messageBytes; ++i) message[i] = (uint8_t)rng.u32();
    WirehairCodec enc = wirehair_encoder_create(nullptr, message.data(), messageBytes, bb);
    WirehairCodec dec = wirehair_decoder_create(nullptr, messageBytes, bb);
    printf("repro N=%d bb=%d final=%d msgBytes=%llu startMode=%d loss=%.2f enc=%p dec=%p\n",
           N, bb, finalBytes, (unsigned long long)messageBytes, startMode, loss, (void*)enc, (void*)dec);
    if (!enc || !dec) {
        wirehair_free(enc);
        wirehair_free(dec);
        return 3;
    }
    const unsigned maxNeeded = (unsigned)N + 512 + (unsigned)N;
    unsigned blockId = (startMode==1)?repair_start_block_id(rng, N, maxNeeded):0, needed=0;
    for (;;) {
        if (needed >= maxNeeded) { printf("GIVING UP needed=%u\n", needed); break; }
        bool drop = (startMode==0) && (rng.unit() < loss);
        unsigned thisId = blockId++;
        if (drop) continue;
        uint32_t wl=0;
        WirehairResult er = wirehair_encode(enc, thisId, block.data(), bb, &wl);
        if (er != Wirehair_Success) {
            printf("ENCODE ERROR er=%d\n", er);
            wirehair_free(enc); wirehair_free(dec);
            return 4;
        }
        needed++;
        WirehairResult dr = wirehair_decode(dec, thisId, block.data(), wl);
        if (needed <= 8 || needed % 256 == 0 || dr != Wirehair_NeedMore)
            printf("  id=%u needed=%u dr=%d (%s)\n", thisId, needed, dr, wirehair_result_string(dr));
        if (dr == Wirehair_Success) {
            WirehairResult rr = wirehair_recover(dec, decoded.data(), messageBytes);
            int cmp = memcmp(decoded.data(), message.data(), messageBytes);
            printf("DECODED at needed=%u (extra=%d) recover=%d cmp=%d\n", needed, (int)needed-N, rr, cmp);
            wirehair_free(enc); wirehair_free(dec);
            return (rr == Wirehair_Success && cmp == 0) ? 0 : 6;
        }
        if (dr != Wirehair_NeedMore) { printf("DECODE ERROR dr=%d\n", dr); wirehair_free(enc); wirehair_free(dec); return 5; }
    }
    wirehair_free(enc); wirehair_free(dec);
    return 7;
}

// Overhead characterization sweep: for each N, run many round-trips and report
// mean + tail overhead. Flags N values with anomalously high overhead (seed/table weak spots).
static int cmd_ohead(int argc, char** argv) {
    int threads = resolve_threads();
    int nlo = 64, nhi = 4000, nstep = 0; // nstep 0 => ~64 log-spaced points
    long trials = 20000;                 // per N
    int bb = 64;                         // small blocks => overhead dominates, fast
    int startMode = 1;                   // repair-only: pure code overhead (most seed-sensitive)
    double loss = 0.10;
    uint64_t baseSeed = 0xACE0;
    string samples_path, nfile;
    bool range_explicit = false;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--nlo") && i+1<argc) { if (!parse_int_strict(argv[++i], nlo)) return 2; range_explicit = true; }
        else if (!strcmp(argv[i], "--nhi") && i+1<argc) { if (!parse_int_strict(argv[++i], nhi)) return 2; range_explicit = true; }
        else if (!strcmp(argv[i], "--nstep") && i+1<argc) { if (!parse_int_strict(argv[++i], nstep)) return 2; range_explicit = true; }
        else if (!strcmp(argv[i], "--nfile") && i+1<argc) nfile = argv[++i];
        else if (!strcmp(argv[i], "--trials") && i+1<argc) { if (!parse_long_strict(argv[++i], trials)) return 2; }
        else if (!strcmp(argv[i], "--bb") && i+1<argc) { if (!parse_int_strict(argv[++i], bb)) return 2; }
        else if (!strcmp(argv[i], "--startmode") && i+1<argc) { if (!parse_int_strict(argv[++i], startMode)) return 2; }
        else if (!strcmp(argv[i], "--loss") && i+1<argc) { if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if (!strcmp(argv[i], "--seed") && i+1<argc) { if (!parse_u64_strict(argv[++i], baseSeed)) return 2; }
        else if (!strcmp(argv[i], "--samples-out") && i+1<argc) samples_path = argv[++i];
    }
    if ((!nfile.empty() && range_explicit) || nstep < 0 || trials < 1 ||
        bb < 0 || (startMode != 0 && startMode != 1) ||
        !std::isfinite(loss) || loss < 0.0 || loss > 0.99)
    {
        fprintf(stderr, "ohead: invalid parameters (--nfile cannot be combined with range options)\n");
        return 2;
    }
    vector<int> Ns;
    if (!nfile.empty()) {
        if (!read_n_file_strict(nfile, "ohead", Ns)) return 2;
        sort(Ns.begin(), Ns.end());
        if (adjacent_find(Ns.begin(), Ns.end()) != Ns.end()) {
            fprintf(stderr, "ohead: duplicate N in --nfile %s\n", nfile.c_str());
            return 2;
        }
        nlo = Ns.front();
        nhi = Ns.back();
    }
    else if (nlo < 2 || nhi > 64000 || nlo > nhi ||
             !valid_seedcheck_args(nlo, bb, startMode, loss))
    {
        fprintf(stderr, "ohead: invalid range or parameters\n");
        return 2;
    }
    else if (nstep > 0) {
        for (int64_t n = nlo; n <= nhi; n += (int64_t)nstep) Ns.push_back((int)n);
    }
    else if (nlo == nhi) { Ns.push_back(nlo); }
    else {
        const int pts = 64;
        const double lo = log((double)nlo), hi = log((double)nhi);
        for (int i = 0; i < pts; ++i) {
            int n = (int)llround(exp(lo + (hi - lo) * i / (pts - 1)));
            if (i == 0) n = nlo;
            if (i == pts - 1) n = nhi;
            if (n < nlo) n = nlo;
            if (n > nhi) n = nhi;
            if (Ns.empty() || n > Ns.back()) Ns.push_back(n);
        }
    }
    FILE* samples_file = nullptr;
    if (!samples_path.empty()) {
        samples_file = fopen(samples_path.c_str(), "w");
        if (!samples_file) {
            fprintf(stderr, "ohead: cannot open --samples-out %s\n", samples_path.c_str());
            return 2;
        }
        fprintf(samples_file, "# whx-ohead-samples-v1 columns=N trial extra\n");
    }
    printf("# ohead: threads=%d trials/N=%ld bb=%d startMode=%d loss=%.2f seed=0x%llx seed_fixups=%s Ns=%zu %s[%d,%d]\n",
           threads, trials, bb, startMode, loss,
           (unsigned long long)baseSeed, seed_fixup_policy(), Ns.size(),
           nfile.empty() ? "range" : "nfile-range", nlo, nhi);
    printf("%-8s %10s %6s %6s %6s %6s %10s\n", "N", "mean", "p50", "p99", "p999", "max", "fail");
    for (int N : Ns) {
        std::atomic<uint64_t> next{0};
        Stat global; std::mutex mu; std::atomic<long> fails{0};
        vector<uint32_t> raw_samples;
        if (samples_file) raw_samples.assign((size_t)trials, UINT32_MAX);
        vector<thread> ts;
        for (int t = 0; t < threads; ++t) {
            ts.emplace_back([&, t]() {
                Stat loc; long lf = 0;
                for (;;) {
                    const uint64_t idx = next.fetch_add(1);
                    if (idx >= (uint64_t)trials) break;
                    uint32_t extra = 0;
                    uint64_t s = baseSeed + (uint64_t)N * 0x9E3779B1u + (uint64_t)idx * 0xD1B54A32D192ED03ull; // distinct strides: equal strides made seed(N, idx) == seed(N+d, idx-d)
                    int rc = roundtrip(s, N, bb, startMode, loss, &extra, nullptr);
                    if (rc == 0) {
                        loc.add_overhead(extra);
                        if (samples_file) raw_samples[(size_t)idx] = extra;
                    }
                    else lf++;
                }
                std::lock_guard<std::mutex> lk(mu); global.merge(loc); fails += lf;
            });
        }
        for (auto& th : ts) th.join();
        double mean = global.trials ? (double)global.overhead_sum/global.trials : 0;
        printf("%-8d %20.17g %6d %6d %6d %6u %10ld\n", N, mean,
               (int)histogram_quantile(global.overhead_hist, 0.50L),
               (int)histogram_quantile(global.overhead_hist, 0.99L),
               (int)histogram_quantile(global.overhead_hist, 0.999L),
               global.overhead_max, fails.load());
        if (samples_file) {
            for (uint64_t idx = 0; idx < (uint64_t)trials; ++idx) {
                const uint32_t extra = raw_samples[(size_t)idx];
                if (extra != UINT32_MAX) {
                    fprintf(samples_file, "%d %llu %u\n", N,
                            (unsigned long long)idx, extra);
                }
            }
        }
        fflush(stdout);
    }
    if (samples_file && fclose(samples_file) != 0) {
        fprintf(stderr, "ohead: failed writing --samples-out %s\n", samples_path.c_str());
        return 1;
    }
    return 0;
}

#ifdef WH_SEED_KNOBS
extern "C" void wh_set_override(int N, int dense, int pseed, int dseed);

// Run `trials` round-trips at fixed N with peel-seed `pseed` and dense-seed `dseed` (<0 = table
// default for each); returns mean overhead (1e9 if any decode failed). Uses thread-local override.
static double seed_mean(int N, int bb, int startMode, double loss, int pseed, int dseed, long trials, uint64_t base) {
    Rng rng(base);
    const int ovrN = (pseed < 0 && dseed < 0) ? -1 : N;
    uint64_t osum = 0; long fail = 0;

    if (bb == 0) {
        wh_set_override(ovrN, -1, pseed, dseed);
        if (!matrix_encoder_ok(N, nullptr)) {
            wh_set_override(-1, -1, -1, -1);
            return 1e9;
        }
        for (long k = 0; k < trials; ++k) {
            uint32_t extra = 0;
            int rc = matrix_decode_roundtrip(rng.next(), N, startMode, loss, &extra, nullptr);
            if (rc == 0) osum += extra;
            else { fail++; break; }
        }
        wh_set_override(-1, -1, -1, -1);
        return fail ? 1e9 : (double)osum / trials;
    }

    for (long k = 0; k < trials; ++k) {
        wh_set_override(ovrN, -1, pseed, dseed);
        uint32_t extra = 0;
        int rc = roundtrip(rng.next(), N, bb, startMode, loss, &extra, nullptr);
        if (rc == 0) osum += extra;
        else { fail++; break; }
    }
    wh_set_override(-1, -1, -1, -1);
    return fail ? 1e9 : (double)osum / trials;
}

// Parallel version: split `trials` across `threads` (each sets its own thread-local override).
static const uint64_t kValidationTrialStride = 2654435761ull;
static uint64_t validation_trial_seed(uint64_t base, uint64_t trial) {
    return base + trial * kValidationTrialStride;
}

static double seed_mean_par(int N, int bb, int startMode, double loss, int pseed, int dseed, long trials,
                            uint64_t base, int threads) {
    std::atomic<uint64_t> idx{0};
    std::atomic<uint64_t> osum{0};
    std::atomic<long> fail{0};
    std::atomic<bool> stop{false};
    const int ovrN = (pseed < 0 && dseed < 0) ? -1 : N;
    vector<thread> ts;
    for (int t = 0; t < threads; ++t) ts.emplace_back([&]() {
        uint64_t lsum = 0; long lf = 0;
        if (bb == 0) {
            wh_set_override(ovrN, -1, pseed, dseed);
            if (!matrix_encoder_ok(N, nullptr)) {
                fail++;
                stop.store(true);
                wh_set_override(-1, -1, -1, -1);
                return;
            }
        }
        for (;;) {
            if (stop.load()) break;
            const uint64_t k = idx.fetch_add(128);
            if (k >= (uint64_t)trials) break;
            uint64_t end = k + 128;
            if (end > (uint64_t)trials) end = (uint64_t)trials;
            for (uint64_t j = k; j < end; ++j) {
                if (stop.load()) break;
                uint32_t extra = 0;
                int rc;
                if (bb == 0) {
                    rc = matrix_decode_roundtrip(validation_trial_seed(base, j), N, startMode, loss, &extra, nullptr);
                }
                else {
                    wh_set_override(ovrN, -1, pseed, dseed);
                    rc = roundtrip(validation_trial_seed(base, j), N, bb, startMode, loss, &extra, nullptr);
                }
                if (rc == 0) lsum += extra;
                else {
                    lf++;
                    stop.store(true);
                    break;
                }
            }
        }
        osum += lsum; fail += lf;
        wh_set_override(-1, -1, -1, -1);
    });
    for (auto& th : ts) th.join();
    return fail.load() ? 1e9 : (double)osum.load() / trials;
}

struct FinalistValidationPlan {
    double primary_loss;
    double secondary_loss;
    uint64_t primary_base;
    uint64_t secondary_base;
};

struct FinalistScores {
    double primary;
    double secondary;
};

template <typename Scorer>
static FinalistScores validate_finalist(const FinalistValidationPlan& plan,
                                        int pseed, int dseed, Scorer scorer) {
    FinalistScores result;
    result.primary = scorer(plan.primary_loss, pseed, dseed, plan.primary_base);
    result.secondary = scorer(plan.secondary_loss, pseed, dseed, plan.secondary_base);
    return result;
}

static bool finalist_is_strictly_better(const FinalistScores& candidate,
                                        const FinalistScores& incumbent) {
    return std::max(candidate.primary, candidate.secondary) <
           std::max(incumbent.primary, incumbent.secondary);
}

static int cmd_seedmean(int argc, char** argv) {
    int threads = resolve_threads();
    int N = 0, pseed = -1, dseed = -1;
    int bb = 64, startMode = 0;
    long trials = 10000;
    double loss = 0.10;
    uint64_t base = 0x51EED00DULL;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--N") && i + 1 < argc) { if (!parse_int_strict(argv[++i], N)) return 2; }
        else if (!strcmp(argv[i], "--n") && i + 1 < argc) { if (!parse_int_strict(argv[++i], N)) return 2; }
        else if (!strcmp(argv[i], "--pseed") && i + 1 < argc) { if (!parse_int_strict(argv[++i], pseed)) return 2; }
        else if (!strcmp(argv[i], "--dseed") && i + 1 < argc) { if (!parse_int_strict(argv[++i], dseed)) return 2; }
        else if (!strcmp(argv[i], "--trials") && i + 1 < argc) { if (!parse_long_strict(argv[++i], trials)) return 2; }
        else if (!strcmp(argv[i], "--bb") && i + 1 < argc) { if (!parse_int_strict(argv[++i], bb)) return 2; }
        else if (!strcmp(argv[i], "--loss") && i + 1 < argc) { if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if (!strcmp(argv[i], "--startmode") && i + 1 < argc) { if (!parse_int_strict(argv[++i], startMode)) return 2; }
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) { if (!parse_u64_strict(argv[++i], base)) return 2; }
    }
    if (N < 2 || N > 64000 || bb < 0 || trials < 1 ||
        pseed < -1 || pseed > 255 || dseed < -1 || dseed > 255 ||
        !std::isfinite(loss) || loss < 0.0 || loss > 0.99 ||
        (startMode != 0 && startMode != 1))
    {
        fprintf(stderr, "seedmean requires 2 <= --N <= 64000, --bb >= 0, positive --trials, seeds in [-1,255], 0 <= --loss <= 0.99, and --startmode 0 or 1\n");
        return 2;
    }
    const double mean = seed_mean_par(N, bb, startMode, loss, pseed, dseed, trials, base, threads);
    printf("# seedmean: threads=%d trials=%ld bb=%d startMode=%d loss=%.2f seed=0x%llx\n",
           threads, trials, bb, startMode, loss, (unsigned long long)base);
    printf("%-8s %8s %8s %8s %12s %8s\n", "N", "pseed", "dseed", "trials", "mean", "failed");
    printf("%-8d %8d %8d %8ld %12.4f %8d\n", N, pseed, dseed, trials, mean, mean >= 1e8 ? 1 : 0);
    return 0;
}

// For each weak N, search peel seeds [0,255] (parallel across seeds), pick the best, and verify
// it on independent loss patterns. Emits a correction table: "N best_pseed default_mean best_mean".
static int cmd_seedsearch(int argc, char** argv) {
    int threads = resolve_threads();
    long tsearch = 800, tverify = 10000;
    int bb = 64, startMode = 0; double loss = 0.10;
    int nseeds = 96; // candidate peel seeds to search [0,nseeds); ~all threads busy, plenty since almost any non-default seed is good
    string nlist, nfile;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i],"--nlist")&&i+1<argc) nlist=argv[++i];
        else if (!strcmp(argv[i],"--nfile")&&i+1<argc) nfile=argv[++i];
        else if (!strcmp(argv[i],"--tsearch")&&i+1<argc) { if (!parse_long_strict(argv[++i], tsearch)) return 2; }
        else if (!strcmp(argv[i],"--tverify")&&i+1<argc) { if (!parse_long_strict(argv[++i], tverify)) return 2; }
        else if (!strcmp(argv[i],"--nseeds")&&i+1<argc) { if (!parse_int_strict(argv[++i], nseeds)) return 2; }
        else if (!strcmp(argv[i],"--bb")&&i+1<argc) { if (!parse_int_strict(argv[++i], bb)) return 2; }
        else if (!strcmp(argv[i],"--loss")&&i+1<argc) { if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if (!strcmp(argv[i],"--startmode")&&i+1<argc) { if (!parse_int_strict(argv[++i], startMode)) return 2; }
    }
    if (nseeds < 1) {
        fprintf(stderr, "seedsearch requires --nseeds >= 1\n");
        return 2;
    }
    if (nseeds > 256) nseeds = 256;
    vector<int> Ns;
    if (!nfile.empty()) {
        if (!read_n_file_strict(nfile, "seedsearch", Ns)) return 2;
    }
    if (!nlist.empty()) {
        vector<int> listed;
        if (!parse_positive_int_list(nlist, listed)) {
            fprintf(stderr, "seedsearch: bad --nlist\n");
            return 2;
        }
        Ns.insert(Ns.end(), listed.begin(), listed.end());
    }
    int dseeds = 0;          // 0 = peel-only (back-compat); >0 = enable joint dense-seed search for hard N
    double goodthr = 0.05;   // peel-only is accepted if max(v10,v30) < goodthr; else search dense
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i],"--dseeds")&&i+1<argc) { if (!parse_int_strict(argv[++i], dseeds)) return 2; }
        else if (!strcmp(argv[i],"--goodthr")&&i+1<argc) { if (!parse_double_strict(argv[++i], goodthr)) return 2; }
    }
    if (Ns.empty() || tsearch < 1 || tverify < 1 || bb < 0 ||
        !std::isfinite(loss) || loss != 0.10 ||
        (startMode != 0 && startMode != 1) ||
        dseeds < 0 || !std::isfinite(goodthr) || goodthr < 0.0)
    {
        fprintf(stderr, "seedsearch: invalid parameters (legacy output requires --loss exactly 0.10)\n");
        return 2;
    }
    for (int N : Ns) {
        if (N < 2 || N > 64000) {
            fprintf(stderr, "seedsearch requires every N in 2..64000\n");
            return 2;
        }
    }
    if (dseeds > 256) dseeds = 256;
    fprintf(stderr,"# seedsearch: threads=%d Ns=%zu bb=%d nseeds=%d dseeds=%d tsearch=%ld tverify=%ld goodthr=%.3f primary_loss=0.10 secondary_loss=0.30 pairing=paired seed_fixups=%s\n",
            threads,Ns.size(),bb,nseeds,dseeds,tsearch,tverify,goodthr,
            seed_fixup_policy());
    printf("# schema=whx-seedsearch-v1 primary_loss=0.10 secondary_loss=0.30 finalist_trials=paired tie_break=peel_then_lowest_seed\n");
    printf("# pairing=paired validation_trials=%ld trial_stride=%llu primary_base=0xBEE5+N*577 secondary_base=0xF00D+N*331\n",
           tverify, (unsigned long long)kValidationTrialStride);
    printf("# N  best_pseed  best_dseed  default_mean  best_mean@0.10  best_mean@0.30\n");
    for (int N : Ns) {
        // Stage 1: search peel seeds (default dense). One thread per candidate seed.
        vector<double> sm(256, 1e9);
        std::atomic<int> si{0};
        vector<thread> ts;
        for (int t=0;t<threads;++t) ts.emplace_back([&](){
            for(;;){ int s=si.fetch_add(1); if(s>=nseeds) break;
                sm[s]=seed_mean(N,bb,startMode,loss,s,-1,tsearch,0x5EED + (uint64_t)N*911); // base independent of s: paired loss patterns across candidates (avoids winner's curse)
            }
        });
        for(auto&th:ts) th.join();
        int bestP=0; for(int s=1;s<nseeds;++s) if(sm[s]<sm[bestP]) bestP=s;
        int bestD=-1;
        const FinalistValidationPlan validation = {
            loss, 0.30,
            0xBEE5 + (uint64_t)N * 577,
            0xF00D + (uint64_t)N * 331
        };
        auto scorer = [&](double validation_loss, int pseed, int dseed, uint64_t base) {
            return seed_mean_par(N, bb, startMode, validation_loss, pseed, dseed,
                                 tverify, base, threads);
        };
        double defMean = scorer(validation.primary_loss, -1, -1, validation.primary_base);
        FinalistScores peel_scores = validate_finalist(validation, bestP, -1, scorer);
        double v10 = peel_scores.primary;
        double v30 = peel_scores.secondary;

        // Stage 2 (Task5): peel-only insufficient -> joint (peel,dense) search over top-K peels x dense seeds.
        if (dseeds > 0 && (v10 >= goodthr || v30 >= goodthr)) {
            // top-K peels by stage-1 @0.10 mean
            vector<int> order(nseeds); for(int s=0;s<nseeds;++s) order[s]=s;
            std::sort(order.begin(), order.end(), [&](int a,int b){
                return sm[a] < sm[b] || (sm[a] == sm[b] && a < b);
            });
            const int K = order.size() < 4 ? (int)order.size() : 4;
            vector<int> peels(order.begin(), order.begin()+K);
            // score each (peel,dense) by the harder loss (0.30); parallelize across combos
            const int ncomb = K * dseeds;
            vector<double> cs(ncomb, 1e9);
            std::atomic<int> ci{0};
            vector<thread> ts2;
            for (int t=0;t<threads;++t) ts2.emplace_back([&](){
                for(;;){ int c=ci.fetch_add(1); if(c>=ncomb) break;
                    int p=peels[c/dseeds], d=c%dseeds;
                    cs[c]=seed_mean(N,bb,startMode,0.30,p,d,tsearch,0xDEED + (uint64_t)N*701); // base independent of c: paired trials across candidates
                }
            });
            for(auto&th:ts2) th.join();
            int bc=0; for(int c=1;c<ncomb;++c) if(cs[c]<cs[bc]) bc=c;
            int jp=peels[bc/dseeds], jd=bc%dseeds;
            const FinalistScores joint_scores = validate_finalist(validation, jp, jd, scorer);
            double jv10 = joint_scores.primary;
            double jv30 = joint_scores.secondary;
            // keep joint result if it beats peel-only on the harder metric
            if (finalist_is_strictly_better(joint_scores, peel_scores)) {
                bestP=jp; bestD=jd; v10=jv10; v30=jv30;
            }
        }
        printf("%-6d %4d %5d %12.4f %14.4f %14.4f\n", N, bestP, bestD, defMean, v10, v30);
        fflush(stdout);
    }
    return 0;
}
#endif

// Scan a range of N for weak-seed spots: parallelize ACROSS N (each thread tests whole N's),
// flag N whose mean overhead exceeds a threshold. Efficient full-range screening.
static int cmd_scan(int argc, char** argv) {
    int threads = resolve_threads();
    int nlo = 2048, nhi = 64000;
    string nfile;
    bool range_explicit = false;
    long trials = 800;
    int bb = 64, startMode = 0;
    double loss = 0.10, thresh = 0.08;
    uint64_t baseSeed = 0x5CA4;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i],"--nlo")&&i+1<argc) { if (!parse_int_strict(argv[++i], nlo)) return 2; range_explicit = true; }
        else if (!strcmp(argv[i],"--nhi")&&i+1<argc) { if (!parse_int_strict(argv[++i], nhi)) return 2; range_explicit = true; }
        else if (!strcmp(argv[i],"--nfile")&&i+1<argc) nfile = argv[++i];
        else if (!strcmp(argv[i],"--trials")&&i+1<argc) { if (!parse_long_strict(argv[++i], trials)) return 2; }
        else if (!strcmp(argv[i],"--bb")&&i+1<argc) { if (!parse_int_strict(argv[++i], bb)) return 2; }
        else if (!strcmp(argv[i],"--startmode")&&i+1<argc) { if (!parse_int_strict(argv[++i], startMode)) return 2; }
        else if (!strcmp(argv[i],"--loss")&&i+1<argc) { if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if (!strcmp(argv[i],"--thresh")&&i+1<argc) { if (!parse_double_strict(argv[++i], thresh)) return 2; }
        else if (!strcmp(argv[i],"--seed")&&i+1<argc) { if (!parse_u64_strict(argv[++i], baseSeed)) return 2; }
    }
    if ((!nfile.empty() && range_explicit) || trials < 1 || bb < 0 ||
        (startMode != 0 && startMode != 1) ||
        loss < 0.0 || loss > 0.99 || !std::isfinite(loss) ||
        thresh < 0.0 || !std::isfinite(thresh))
    {
        fprintf(stderr, "scan: invalid parameters (--nfile cannot be combined with --nlo/--nhi)\n");
        return 2;
    }
    vector<int> Ns;
    if (!nfile.empty()) {
        if (!read_n_file_strict(nfile, "scan", Ns)) return 2;
        sort(Ns.begin(), Ns.end());
        if (adjacent_find(Ns.begin(), Ns.end()) != Ns.end()) {
            fprintf(stderr, "scan: duplicate N in --nfile %s\n", nfile.c_str());
            return 2;
        }
    }
    else {
        if (nlo < 2 || nhi > 64000 || nlo > nhi ||
            !valid_seedcheck_args(nlo, bb, startMode, loss))
        {
            fprintf(stderr, "scan: invalid range or parameters\n");
            return 2;
        }
        Ns.reserve((size_t)(nhi - nlo + 1));
        for (int N = nlo; N <= nhi; ++N) Ns.push_back(N);
    }
    const int active_threads = std::min<int>(threads, (int)Ns.size());
    if (!nfile.empty()) {
        printf("# scan: threads=%d source=%s count=%zu trials/N=%ld bb=%d startMode=%d loss=%.2f thresh=%.3f seed=0x%llx seed_fixups=%s\n",
               active_threads, nfile.c_str(), Ns.size(), trials, bb, startMode, loss, thresh,
               (unsigned long long)baseSeed,
               seed_fixup_policy());
    }
    else {
        printf("# scan: threads=%d N=[%d,%d] trials/N=%ld bb=%d startMode=%d loss=%.2f thresh=%.3f seed=0x%llx seed_fixups=%s\n",
               active_threads, nlo, nhi, trials, bb, startMode, loss, thresh,
               (unsigned long long)baseSeed,
               seed_fixup_policy());
    }
    std::atomic<size_t> next_index{0};
    std::mutex pmu;
    std::atomic<long> scanned{0};
    struct ScanResult {
        uint64_t overhead_sum = 0;
        uint32_t overhead_max = 0;
        long failures = 0;
    };
    vector<ScanResult> scan_results(Ns.size());
    vector<thread> ts;
    for (int t = 0; t < active_threads; ++t) {
        ts.emplace_back([&]() {
            for (;;) {
                const size_t index = next_index.fetch_add(1);
                if (index >= Ns.size()) break;
                const int N = Ns[index];
                uint64_t osum = 0; uint32_t omax = 0; long fail = 0;
                for (long k = 0; k < trials; ++k) {
                    uint32_t extra = 0;
                    uint64_t s = baseSeed + (uint64_t)N * 0x9E3779B1u + (uint64_t)k * 0xD1B54A32D192ED03ull; // distinct strides: equal strides made seed(N, k) == seed(N+d, k-d)
                    int rc = roundtrip(s, N, bb, startMode, loss, &extra, nullptr);
                    if (rc == 0) { osum += extra; if (extra > omax) omax = extra; } else fail++;
                }
                scan_results[index].overhead_sum = osum;
                scan_results[index].overhead_max = omax;
                scan_results[index].failures = fail;
                long sc = ++scanned;
                if (sc % 4000 == 0) { std::lock_guard<std::mutex> lk(pmu); fprintf(stderr, "  ...scanned %ld N\n", sc); }
            }
        });
    }
    for (auto& th : ts) th.join();
    for (size_t index = 0; index < Ns.size(); ++index) {
        const ScanResult& result = scan_results[index];
        const long double mean =
            (long double)result.overhead_sum / (long double)trials;
        if (mean > (long double)thresh || result.failures > 0) {
            printf("WEAK N=%-6d mean=%.17Lg max=%u fail=%ld\n",
                   Ns[index], mean, result.overhead_max, result.failures);
        }
    }
    printf("# scan done: %ld N scanned\n", scanned.load());
    return 0;
}

#ifdef WH_COUNT
extern "C" {
void wh_stage_reset(); int wh_stage_count(); uint64_t wh_stage_bytes(int); const char* wh_stage_label(int);
unsigned wh_graph_defer_count(); unsigned wh_graph_defer_rows(); unsigned wh_graph_component_count();
unsigned wh_graph_max_component(); uint64_t wh_graph_component_sum_squares();
}

static unsigned percentile_value(const vector<unsigned>& sorted, double p) {
    if (sorted.empty()) return 0;
    return sorted[(size_t)nearest_rank_index((uint64_t)sorted.size(), (long double)p)];
}

static void summarize_values(const vector<unsigned>& values, double& mean, double& sd,
                             unsigned& minv, unsigned& p50, unsigned& p95, unsigned& p99,
                             unsigned& maxv) {
    if (values.empty()) {
        mean = sd = 0.0;
        minv = p50 = p95 = p99 = maxv = 0;
        return;
    }

    vector<unsigned> sorted = values;
    sort(sorted.begin(), sorted.end());
    uint64_t sum = 0;
    for (unsigned v : sorted) sum += v;
    mean = (double)sum / sorted.size();
    double var = 0.0;
    for (unsigned v : sorted) {
        const double d = (double)v - mean;
        var += d * d;
    }
    sd = sqrt(var / sorted.size());
    minv = sorted.front();
    p50 = percentile_value(sorted, 0.50);
    p95 = percentile_value(sorted, 0.95);
    p99 = percentile_value(sorted, 0.99);
    maxv = sorted.back();
}

static int cmd_peelstat(int argc, char** argv) {
    int threads = resolve_threads();
    string nlist = "128,2048,32000";
    string bblist = "64";
    int trials = 200;
    double loss = 0.10;
    int startMode = 0;
    uint64_t baseSeed = 0x9A7E11AULL;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--N") && i + 1 < argc) nlist = argv[++i];
        else if (!strcmp(argv[i], "--bb") && i + 1 < argc) bblist = argv[++i];
        else if (!strcmp(argv[i], "--bb-list") && i + 1 < argc) bblist = argv[++i];
        else if (!strcmp(argv[i], "--trials") && i + 1 < argc) { if (!parse_int_strict(argv[++i], trials)) return 2; }
        else if (!strcmp(argv[i], "--loss") && i + 1 < argc) { if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if (!strcmp(argv[i], "--startmode") && i + 1 < argc) { if (!parse_int_strict(argv[++i], startMode)) return 2; }
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) { if (!parse_u64_strict(argv[++i], baseSeed)) return 2; }
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) {
            if (!parse_int_strict(argv[++i], g_threads) || g_threads < 1) return 2;
            threads = resolve_threads();
        }
    }
    vector<int> Ns;
    vector<int> BBs;
    if (!parse_positive_int_list(nlist, Ns) ||
        !parse_positive_int_list(bblist, BBs) ||
        !validate_n_bb_lists(Ns, BBs, "peelstat"))
    {
        return 2;
    }
    if (trials < 1 || !std::isfinite(loss) || loss < 0.0 || loss > 0.99 ||
        (startMode != 0 && startMode != 1))
    {
        fprintf(stderr, "peelstat: invalid parameters\n");
        return 2;
    }

    const wirehair::EnvironmentValue mode_environment("WH_PEEL_MODE");
    const char* mode = mode_environment.Get();
    if (!mode) mode = "0";
    printf("# peelstat: mode=%s threads=%d trials=%d loss=%.2f startMode=%d seed=0x%llx\n",
           mode, threads, trials, loss, startMode, (unsigned long long)baseSeed);
    printf("%-6s %-8s %-6s %-8s %-6s %9s %8s %5s %5s %5s %5s %5s %9s %8s %5s %5s %5s %5s %5s\n",
           "mode", "N", "bb", "trials", "fail",
           "cols_mu", "cols_sd", "cmin", "c50", "c95", "c99", "cmax",
           "rows_mu", "rows_sd", "rmin", "r50", "r95", "r99", "rmax");

    for (int bb : BBs) for (int N : Ns) {
        vector<unsigned> cols((size_t)trials);
        vector<unsigned> rows((size_t)trials);
        vector<uint8_t> ok((size_t)trials, 0);
        std::atomic<int> next{0};
        std::atomic<int> fails{0};
        vector<thread> ts;
        for (int t = 0; t < threads; ++t) {
            ts.emplace_back([&]() {
                char failmsg[256];
                for (;;) {
                    const int idx = next.fetch_add(1);
                    if (idx >= trials) break;
                    uint32_t extra = 0;
                    const uint64_t seed = baseSeed + (uint64_t)N * 0x9e3779b1u +
                        (uint64_t)bb * 0x85ebca6bu + (uint64_t)idx * 0xc2b2ae35u;
                    failmsg[0] = 0;
                    const int rc = roundtrip(seed, N, bb, startMode, loss, &extra, failmsg);
                    if (rc == 0) {
                        cols[(size_t)idx] = wh_graph_defer_count();
                        rows[(size_t)idx] = wh_graph_defer_rows();
                        ok[(size_t)idx] = 1;
                    }
                    else {
                        fails++;
                    }
                }
            });
        }
        for (auto& th : ts) th.join();

        vector<unsigned> ok_cols;
        vector<unsigned> ok_rows;
        ok_cols.reserve((size_t)trials);
        ok_rows.reserve((size_t)trials);
        for (int i = 0; i < trials; ++i) {
            if (ok[(size_t)i]) {
                ok_cols.push_back(cols[(size_t)i]);
                ok_rows.push_back(rows[(size_t)i]);
            }
        }

        double cm = 0, cs = 0, rm = 0, rs = 0;
        unsigned cmin = 0, c50 = 0, c95 = 0, c99 = 0, cmax = 0;
        unsigned rmin = 0, r50 = 0, r95 = 0, r99 = 0, rmax = 0;
        summarize_values(ok_cols, cm, cs, cmin, c50, c95, c99, cmax);
        summarize_values(ok_rows, rm, rs, rmin, r50, r95, r99, rmax);

        printf("%-6s %-8d %-6d %-8zu %-6d %9.2f %8.2f %5u %5u %5u %5u %5u %9.2f %8.2f %5u %5u %5u %5u %5u\n",
               mode, N, bb, ok_cols.size(), fails.load(),
               cm, cs, cmin, c50, c95, c99, cmax,
               rm, rs, rmin, r50, r95, r99, rmax);
        fflush(stdout);
    }

    return 0;
}

// Symbol-XOR traffic breakdown per codec stage (Task 6a). Single-threaded, one round.
static void count_dump(const char* stage, uint64_t msgBytes) {
    static const char* nm[6] = {"add","add2","addset","mul","muladd","memswap"};
    uint64_t xorB = gf256_count_bytes(0)+gf256_count_bytes(1)+gf256_count_bytes(2);
    uint64_t mulB = gf256_count_bytes(3)+gf256_count_bytes(4);
    uint64_t swpB = gf256_count_bytes(5);
    printf("%-16s msg=%lluB | XOR=%.2fMB (%.1fx msg) GFmul=%.2fMB swap=%.2fMB | total=%.2fMB\n",
        stage, (unsigned long long)msgBytes,
        xorB/1e6, msgBytes? (double)xorB/msgBytes:0, mulB/1e6, swpB/1e6, (xorB+mulB+swpB)/1e6);
    printf("    ");
    for (int i=0;i<6;++i) printf("%s=%.2fMB/%lluc  ", nm[i], gf256_count_bytes(i)/1e6, (unsigned long long)gf256_count_calls(i));
    printf("\n");
    if (wh_graph_defer_count() > 0 && strcmp(stage, "recover")) {
        printf("    deferred_graph: cols=%u rows=%u components=%u max=%u max_frac=%.3f sumsq=%llu\n",
               wh_graph_defer_count(), wh_graph_defer_rows(), wh_graph_component_count(),
               wh_graph_max_component(),
               wh_graph_defer_count() ? (double)wh_graph_max_component() / wh_graph_defer_count() : 0.0,
               (unsigned long long)wh_graph_component_sum_squares());
    }
}
static int cmd_count(int argc, char** argv) {
    string nlist = "2048", bblist = "1300";
    double loss=0.10; int startMode=0;
    for (int i=0;i<argc;++i){ if(!strcmp(argv[i],"--N")&&i+1<argc)nlist=argv[++i];
        else if(!strcmp(argv[i],"--bb")&&i+1<argc)bblist=argv[++i];
        else if(!strcmp(argv[i],"--bb-list")&&i+1<argc)bblist=argv[++i];
        else if(!strcmp(argv[i],"--loss")&&i+1<argc){ if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if(!strcmp(argv[i],"--startmode")&&i+1<argc){ if (!parse_int_strict(argv[++i], startMode)) return 2; } }
    vector<int> Ns;
    vector<int> BBs;
    if (!parse_positive_int_list(nlist, Ns) ||
        !parse_positive_int_list(bblist, BBs) ||
        !validate_n_bb_lists(Ns, BBs, "count"))
    {
        return 2;
    }
    if (!std::isfinite(loss) || loss < 0.0 || loss > 0.99 ||
        (startMode != 0 && startMode != 1))
    {
        fprintf(stderr, "count: invalid parameters\n");
        return 2;
    }
    for (int bb : BBs) for (int N : Ns) {
        Rng rng(0xC007 + (uint64_t)N * 131 + (uint64_t)bb * 17);
        uint64_t msgBytes=(uint64_t)N*bb;
        vector<uint8_t> msg(msgBytes), dec(msgBytes); vector<uint8_t> blk(bb);
        for(uint64_t i=0;i<msgBytes;++i) msg[i]=(uint8_t)rng.u32();
        printf("# count: N=%d bb=%d msg_MiB=%.2f loss=%.2f startMode=%d\n",
               N, bb, (double)msgBytes / 1048576.0, loss, startMode);
        gf256_count_reset(); wh_stage_reset();
        WirehairCodec enc=wirehair_encoder_create(nullptr,msg.data(),msgBytes,bb);
        count_dump("encoder_create", msgBytes);
        { uint64_t attributed=0; for(int i=0;i<wh_stage_count();++i){ attributed+=wh_stage_bytes(i);
            printf("      %-28s %.2fMB (%.1fx)\n", wh_stage_label(i), wh_stage_bytes(i)/1e6, msgBytes?(double)wh_stage_bytes(i)/msgBytes:0); }
          uint64_t tot=gf256_count_bytes(0)+gf256_count_bytes(1)+gf256_count_bytes(2)+gf256_count_bytes(3)+gf256_count_bytes(4)+gf256_count_bytes(5);
          printf("      %-28s %.2fMB (peel-during-feed + misc)\n", "[unattributed]", (tot-attributed)/1e6); }
        WirehairCodec dcd=wirehair_decoder_create(nullptr,msgBytes,bb);
        unsigned id=(startMode==1)?(unsigned)N:0, needed=0; uint32_t wl=0; bool ok=false;
        gf256_count_reset();
        // Attribute encoder repair-block generation separately so the
        // decode(feed) stage reports decoder traffic only.
        uint64_t encB=0, encPreB=0;
        for(;;){ if(needed>=(unsigned)N*2+512)break;
            bool drop=(startMode==0)&&(rng.unit()<loss); unsigned t=id++; if(drop)continue;
            encPreB=0; for(int k2=0;k2<6;++k2) encPreB+=gf256_count_bytes(k2);
            wirehair_encode(enc,t,blk.data(),bb,&wl); needed++;
            { uint64_t post=0; for(int k2=0;k2<6;++k2) post+=gf256_count_bytes(k2); encB+=post-encPreB; }
            WirehairResult dr=wirehair_decode(dcd,t,blk.data(),wl);
            if(dr==Wirehair_Success){ok=true;break;} if(dr!=Wirehair_NeedMore)break; }
        printf("%-16s msg=%lluB | total=%.2fMB (excluded from decode stage below)\n",
            "encode(repair)", (unsigned long long)msgBytes, encB/1e6);
        count_dump("decode(feed)+enc", msgBytes);
        if(ok){ gf256_count_reset(); wirehair_recover(dcd,dec.data(),msgBytes); count_dump("recover", msgBytes);
            printf("# recover correct: %s\n", memcmp(dec.data(),msg.data(),msgBytes)==0?"YES":"NO"); }
        wirehair_free(enc); wirehair_free(dcd);
    }
    return 0;
}
#endif

static int cmd_selftest(int argc, char**) {
    if (argc != 0) {
        fprintf(stderr, "selftest takes no options\n");
        return 2;
    }
    int assertions = 0;
    auto expect = [&](bool condition, const char* what) {
        ++assertions;
        if (!condition) fprintf(stderr, "SELFTEST FAIL: %s\n", what);
        return condition;
    };
    bool ok = true;

#if defined(WIREHAIR_DISABLE_SEED_FIXUPS)
    ok &= expect(!wirehair::SeedFixupsEnabled(),
                 "harness and codec both disable seed fixups");
    printf("# seed_fixups=disabled\n");
#else
    ok &= expect(wirehair::SeedFixupsEnabled(),
                 "harness and codec both enable seed fixups");
    printf("# seed_fixups=enabled\n");
#endif

    vector<uint64_t> hist;
    ok &= expect(histogram_quantile(hist, 0.5L) == 0, "empty histogram");
    hist.assign(4, 0); hist[1] = 1;
    ok &= expect(histogram_quantile(hist, 0.0L) == 1, "single sample p0");
    ok &= expect(histogram_quantile(hist, 0.5L) == 1, "single sample p50");
    ok &= expect(histogram_quantile(hist, 0.999L) == 1, "single sample p999");
    hist.assign(5, 0); hist[1] = 1; hist[3] = 1;
    ok &= expect(histogram_quantile(hist, 0.5L) == 1, "two sample p50 boundary");
    ok &= expect(histogram_quantile(hist, 0.5001L) == 3, "two sample above p50");
    ok &= expect(histogram_quantile(hist, 1.0L) == 3, "two sample p100");
    hist.assign(7, 0); hist[2] = 2; hist[5] = 3;
    ok &= expect(histogram_quantile(hist, 0.4L) == 2, "repeated bin boundary");
    ok &= expect(histogram_quantile(hist, 0.4001L) == 5, "repeated bin next rank");
    hist.assign(2, 0); hist[0] = std::numeric_limits<uint64_t>::max(); hist[1] = 2;
    ok &= expect(histogram_quantile(hist, 0.5L) == 0, "overflow-safe histogram accumulation");
    ok &= expect(histogram_quantile(hist, 1.0L) == 1, "overflow-safe histogram final rank");

    Stat capped;
    capped.add_overhead(1); capped.add_overhead(63); capped.add_overhead(64); capped.add_overhead(999);
    ok &= expect(histogram_quantile(capped.overhead_hist, 0.25L) == 1, "capped tail first rank");
    ok &= expect(histogram_quantile(capped.overhead_hist, 0.5L) == 63, "capped tail median");
    ok &= expect(histogram_quantile(capped.overhead_hist, 1.0L) == 63, "capped tail maximum bin");

    vector<unsigned> samples = {9, 1, 9, 4, 4, 9, 2};
    vector<uint64_t> oracle_hist(10, 0);
    for (unsigned value : samples) oracle_hist[value]++;
    sort(samples.begin(), samples.end());
    const long double probabilities[] = {0.0L, 0.01L, 0.5L, 0.99L, 0.999L, 1.0L};
    for (long double p : probabilities) {
        const size_t oracle = samples[(size_t)nearest_rank_index(samples.size(), p)];
        ok &= expect(histogram_quantile(oracle_hist, p) == oracle, "histogram sorted-sample oracle");
    }

    uint8_t expected[10], actual[10];
    for (size_t i = 0; i < sizeof(expected); ++i) expected[i] = actual[i] = (uint8_t)i;
    uint64_t mismatch = 0;
    ok &= expect(verify_recovered_bytes(expected, actual, sizeof(expected), mismatch) && mismatch == sizeof(expected),
                 "exact recovered bytes");
    const size_t corruptions[] = {0, 5, 9};
    for (size_t offset : corruptions) {
        memcpy(actual, expected, sizeof(expected)); actual[offset] ^= 1;
        ok &= expect(!verify_recovered_bytes(expected, actual, sizeof(expected), mismatch) && mismatch == offset,
                     "first/middle/final corruption in partial-final-block message");
    }

    BenchWorkerPlan plan;
    ok &= expect(!plan_bench_workers(4, 0, 10, 10, 10000, plan), "zero rounds rejected");
    ok &= expect(plan_bench_workers(8, 1, 10, 10, std::numeric_limits<uint64_t>::max(), plan) &&
                 plan.active_workers == 1,
                 "one round uses one worker");
    const uint64_t exact_worker_bytes = plan.bytes_per_worker;
    ok &= expect(plan_bench_workers(8, 20, 10, 10, exact_worker_bytes * 3, plan) &&
                 plan.active_workers == 3 && plan.bytes_per_worker == exact_worker_bytes,
                 "budget bounds many rounds");
    ok &= expect(plan_bench_workers(2, 20, 10, 10, exact_worker_bytes * 10, plan) && plan.active_workers == 2,
                 "requested thread boundary");
    ok &= expect(plan_bench_workers(8, 8, 10, 10, exact_worker_bytes, plan) && plan.active_workers == 1,
                 "exact one-worker budget");
    ok &= expect(!plan_bench_workers(8, 8, 10, 10, exact_worker_bytes - 1, plan),
                 "insufficient budget rejected");
    ok &= expect(!plan_bench_workers(1, 1, std::numeric_limits<uint64_t>::max(), 2,
                                     std::numeric_limits<uint64_t>::max(), plan),
                 "worker size overflow rejected");
    ok &= expect(plan_bench_workers(1, 1, 64000, 1,
                                    std::numeric_limits<uint64_t>::max(), plan) &&
                 plan.bytes_per_worker > 4 * UINT64_C(64000) + 1,
                 "small-block plan includes codec metadata and matrices");
    ok &= expect(default_threads_for_hw(0) >= 1 && thread_cap_for_hw(0) >= 1,
                 "hardware-thread fallback");

#ifdef WH_SEED_KNOBS
    const FinalistValidationPlan validation = {0.10, 0.30, 111, 222};
    vector<vector<uint64_t> > seen_trials;
    auto fake_scorer = [&](double, int pseed, int dseed, uint64_t base) {
        vector<uint64_t> ids;
        for (uint64_t trial = 0; trial < 3; ++trial) {
            ids.push_back(validation_trial_seed(base, trial));
        }
        seen_trials.push_back(ids);
        return (double)(pseed + dseed + 10) / 100.0;
    };
    const FinalistScores peel = validate_finalist(validation, 1, -1, fake_scorer);
    const FinalistScores joint = validate_finalist(validation, 1, 0, fake_scorer);
    ok &= expect(seen_trials.size() == 4 && seen_trials[0] == seen_trials[2] &&
                 seen_trials[1] == seen_trials[3], "finalists use identical validation trial IDs");
    const FinalistScores tied = {peel.primary, peel.secondary};
    ok &= expect(!finalist_is_strictly_better(tied, peel), "finalist ties retain peel-only");
    const FinalistScores better = {joint.primary - 1.0, joint.secondary - 1.0};
    ok &= expect(finalist_is_strictly_better(better, peel), "strictly better finalist selected");
#endif

    if (!ok) return 1;
    printf("# selftest passed assertions=%d\n", assertions);
    return 0;
}

int main(int argc, char** argv) {
    if (wirehair_init() != Wirehair_Success) { fprintf(stderr, "wirehair_init failed\n"); return 2; }
    if (argc < 2) {
#ifdef WH_COUNT
#ifdef WH_SEED_KNOBS
        fprintf(stderr, "usage: whx [selftest|micro|bench|fuzz|repro|ohead|scan|seedmean|seedsearch|peelstat] [--threads T] [opts]\n");
#else
        fprintf(stderr, "usage: whx [selftest|micro|bench|fuzz|repro|ohead|scan|peelstat] [--threads T] [opts]\n");
#endif
        fprintf(stderr, "  bench/count/peelstat accept --N csv and --bb/--bb-list csv for block-count x block-size sweeps\n");
        fprintf(stderr, "  bench memory is capped by --memory-mib (default 4096; env WHX_BENCH_MEMORY_MIB)\n");
#else
#ifdef WH_SEED_KNOBS
        fprintf(stderr, "usage: whx [selftest|micro|bench|fuzz|repro|ohead|scan|seedmean|seedsearch] [--threads T] [opts]\n");
#else
        fprintf(stderr, "usage: whx [selftest|micro|bench|fuzz|repro|ohead|scan] [--threads T] [opts]\n");
#endif
        fprintf(stderr, "  bench accepts --N csv and --bb/--bb-list csv for block-count x block-size sweeps\n");
        fprintf(stderr, "  bench memory is capped by --memory-mib (default 4096; env WHX_BENCH_MEMORY_MIB)\n");
#endif
        return 1;
    }
    string mode = argv[1];
    // Parse the global --threads for ALL modes first (previously scan/seedsearch/repro were
    // dispatched before this, silently ignoring --threads and defaulting to ~102 threads).
    vector<char*> rest;
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "--threads")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--threads requires a value\n");
                return 2;
            }
            if (!parse_int_strict(argv[++i], g_threads) || g_threads < 1) {
                fprintf(stderr, "invalid --threads\n");
                return 2;
            }
        }
        else rest.push_back(argv[i]);
    }
    int ac = (int)rest.size(); char** av = rest.data();
    if (!validate_mode_options(mode, ac, av)) return 2;
    if (mode == "repro") return cmd_repro(ac, av);
    if (mode == "scan")  return cmd_scan(ac, av);
#ifdef WH_SEED_KNOBS
    if (mode == "seedmean") return cmd_seedmean(ac, av);
    if (mode == "seedsearch") return cmd_seedsearch(ac, av);
#endif
    if (mode == "micro") return cmd_micro(ac, av);
    if (mode == "selftest") return cmd_selftest(ac, av);
    if (mode == "bench") return cmd_bench(ac, av);
    if (mode == "fuzz")  return cmd_fuzz(ac, av);
    if (mode == "ohead") return cmd_ohead(ac, av);
#ifdef WH_COUNT
    if (mode == "count") return cmd_count(ac, av);
    if (mode == "peelstat") return cmd_peelstat(ac, av);
#endif
    fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 1;
}
