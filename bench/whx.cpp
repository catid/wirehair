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

using namespace std;
using Clock = std::chrono::steady_clock;

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
static int kThreadCap() {
    int hw = (int)std::thread::hardware_concurrency();
    if (hw <= 0) hw = 8;
    int cap = hw - 16;          // leave 16 HW threads free
    if (cap > 104) cap = 104;   // and never exceed 104 regardless
    if (cap < 1) cap = 1;
    return cap;
}
static int default_threads() {
    int hw = (int)std::thread::hardware_concurrency();
    if (hw <= 0) hw = 8;
    int t = (hw * 3) / 4; // ~75%
    if (t > kThreadCap()) t = kThreadCap();
    if (t < 1) t = 1;
    return t;
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
struct Stat {
    uint64_t trials = 0;
    uint64_t overhead_sum = 0;       // sum of (needed - N)
    uint64_t overhead_sq = 0;
    uint32_t overhead_max = 0;
    vector<uint32_t> overhead_hist;  // index = extra blocks, capped
    Stat() : overhead_hist(64, 0) {}
    void add_overhead(uint32_t extra) {
        overhead_sum += extra; overhead_sq += (uint64_t)extra * extra;
        if (extra > overhead_max) overhead_max = extra;
        overhead_hist[extra < 64 ? extra : 63]++;
        trials++;
    }
    void merge(const Stat& o) {
        trials += o.trials; overhead_sum += o.overhead_sum; overhead_sq += o.overhead_sq;
        if (o.overhead_max > overhead_max) overhead_max = o.overhead_max;
        for (size_t i = 0; i < overhead_hist.size(); ++i) overhead_hist[i] += o.overhead_hist[i];
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
    static const char* bench_values[] = { "--bb", "--bb-list", "--loss", "--N", "--rounds", nullptr };
    static const char* fuzz_values[] = { "--secs", "--nmax", "--seed", nullptr };
    static const char* ohead_values[] = { "--nlo", "--nhi", "--nstep", "--trials", "--bb", "--startmode", "--loss", "--seed", nullptr };
    static const char* scan_values[] = { "--nlo", "--nhi", "--trials", "--bb", "--startmode", "--loss", "--thresh", "--seed", nullptr };
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
    // overhead percentiles from histogram
    auto pct = [&](double p) -> int {
        uint64_t target = (uint64_t)(p * global.trials);
        uint64_t cum = 0;
        for (int i = 0; i < (int)global.overhead_hist.size(); ++i) { cum += global.overhead_hist[i]; if (cum >= target) return i; }
        return (int)global.overhead_hist.size() - 1;
    };
    printf("# fuzz done: %.1fs trials=%llu fails=%llu rate=%.0f trials/s\n",
           dt, (unsigned long long)trials, (unsigned long long)fails, trials / dt);
    printf("# overhead(success trials=%llu): mean=%.4f p50=%d p99=%d p999=%d max=%u\n",
           (unsigned long long)global.trials, mean, pct(0.50), pct(0.99), pct(0.999), global.overhead_max);
    return fails ? 1 : 0;
}

// ============================================================================
// bench (end-to-end), parallel across threads, fixed work per (N)
// ============================================================================
struct BenchAccum {
    double enc_create_us = 0; uint64_t enc_create_n = 0;
    double encode_us = 0; uint64_t encode_n = 0; uint64_t encode_bytes = 0;
    double decode_us = 0; uint64_t decode_n = 0; uint64_t decode_bytes = 0;
    double recover_us = 0; uint64_t recover_n = 0; uint64_t recover_bytes = 0;
    uint64_t overhead_sum = 0; uint64_t rounds = 0; uint64_t fails = 0;
    void merge(const BenchAccum& o){
        enc_create_us+=o.enc_create_us; enc_create_n+=o.enc_create_n;
        encode_us+=o.encode_us; encode_n+=o.encode_n; encode_bytes+=o.encode_bytes;
        decode_us+=o.decode_us; decode_n+=o.decode_n; decode_bytes+=o.decode_bytes;
        recover_us+=o.recover_us; recover_n+=o.recover_n; recover_bytes+=o.recover_bytes;
        overhead_sum+=o.overhead_sum; rounds+=o.rounds; fails+=o.fails;
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

static void bench_one_round(BenchWorker& w, int N, int blockBytes, double lossRate, uint64_t seed,
                            BenchAccum& a, bool timed) {
    Rng rng(seed);
    uint64_t messageBytes = (uint64_t)N * blockBytes; // full blocks for clean throughput accounting
    // light per-round mutation so the solve isn't trivially identical, cheap
    for (uint64_t i = 0; i < messageBytes; i += 64) w.message[i] = (uint8_t)rng.u32();

    double t0 = now_sec();
    w.enc = wirehair_encoder_create(w.enc, w.message.data(), messageBytes, blockBytes); // reuse
    double t1 = now_sec();
    if (!w.enc) { if (timed) a.fails++; return; }
    if (timed) { a.enc_create_us += (t1 - t0) * 1e6; a.enc_create_n++; }

    w.dec = wirehair_decoder_create(w.dec, messageBytes, blockBytes); // reuse
    if (!w.dec) { if (timed) a.fails++; return; }

    // measure encode of repair blocks (ids N..N+enc_count-1)
    const int enc_count = 64;
    uint32_t wl = 0;
    double e0 = now_sec();
    for (int i = 0; i < enc_count; ++i) {
        if (wirehair_encode(w.enc, N + i, w.block.data(), blockBytes, &wl) !=
            Wirehair_Success)
        {
            if (timed) a.fails++;
            return;
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
        if (timed && rr == Wirehair_Success) {
            a.recover_us += (r1 - r0) * 1e6;
            a.recover_n++;
            a.recover_bytes += messageBytes;
            a.overhead_sum += (needed >= (unsigned)N) ? (needed - N) : 0;
            a.rounds++;
        }
        else if (timed) {
            a.fails++;
        }
    }
    else if (timed) {
        // Failed rounds are counted separately; folding their capped "needed"
        // into overhead_sum would silently shift the overhead statistic.
        a.fails++;
    }
}

static int cmd_bench(int argc, char** argv) {
    int threads = resolve_threads();
    string bblist = "1300";
    double lossRate = 0.10;
    int rounds_per_N = 0; // 0 => auto by N
    bool rounds_set = false;
    string nlist = "32,128,512,1024,2048,8192,32000";
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
    }
    vector<int> Ns;
    vector<int> BBs;
    if (!parse_positive_int_list(nlist, Ns) ||
        !parse_positive_int_list(bblist, BBs) ||
        !validate_n_bb_lists(Ns, BBs, "bench"))
    {
        return 2;
    }
    if (!std::isfinite(lossRate) || lossRate < 0.0 || lossRate >= 1.0 ||
        (rounds_set && rounds_per_N < 1))
    {
        fprintf(stderr, "bench requires 0 <= --loss < 1 and positive --rounds when specified\n");
        return 2;
    }
    const bool matrix = BBs.size() > 1;
    printf("# bench: threads=%d bb=%s loss=%.2f\n", threads, bblist.c_str(), lossRate);
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
        std::atomic<int> next{0};
        BenchAccum global; std::mutex mu;
        vector<thread> ts;
        for (int t = 0; t < threads; ++t) {
            ts.emplace_back([&, t]() {
                BenchWorker w;
                uint64_t mb = (uint64_t)N * blockBytes;
                w.message.assign(mb, (uint8_t)(0x33 + t));
                w.decoded.assign(mb, 0);
                w.block.assign(blockBytes, 0);
                BenchAccum a;
                // warmup (untimed): allocate codecs + pre-fault all pages
                bench_one_round(w, N, blockBytes, lossRate, 0x9999ULL + (uint64_t)N * 131 + t, a, false);
                for (;;) {
                    int idx = next.fetch_add(1);
                    if (idx >= total_rounds) break;
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
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--nlo") && i+1<argc) { if (!parse_int_strict(argv[++i], nlo)) return 2; }
        else if (!strcmp(argv[i], "--nhi") && i+1<argc) { if (!parse_int_strict(argv[++i], nhi)) return 2; }
        else if (!strcmp(argv[i], "--nstep") && i+1<argc) { if (!parse_int_strict(argv[++i], nstep)) return 2; }
        else if (!strcmp(argv[i], "--trials") && i+1<argc) { if (!parse_long_strict(argv[++i], trials)) return 2; }
        else if (!strcmp(argv[i], "--bb") && i+1<argc) { if (!parse_int_strict(argv[++i], bb)) return 2; }
        else if (!strcmp(argv[i], "--startmode") && i+1<argc) { if (!parse_int_strict(argv[++i], startMode)) return 2; }
        else if (!strcmp(argv[i], "--loss") && i+1<argc) { if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if (!strcmp(argv[i], "--seed") && i+1<argc) { if (!parse_u64_strict(argv[++i], baseSeed)) return 2; }
    }
    if (nlo < 2 || nhi > 64000 || nlo > nhi || nstep < 0 ||
        trials < 1 || !valid_seedcheck_args(nlo, bb, startMode, loss))
    {
        fprintf(stderr, "ohead: invalid range or parameters\n");
        return 2;
    }
    vector<int> Ns;
    if (nstep > 0) { for (int n = nlo; n <= nhi; n += nstep) Ns.push_back(n); }
    else { int pts = 64; double lo=log((double)nlo), hi=log((double)nhi); for (int i=0;i<pts;++i){int n=(int)exp(lo+(hi-lo)*i/(pts-1)); if(Ns.empty()||n>Ns.back()) Ns.push_back(n);} }
    printf("# ohead: threads=%d trials/N=%ld bb=%d startMode=%d loss=%.2f Ns=%zu range[%d,%d]\n",
           threads, trials, bb, startMode, loss, Ns.size(), nlo, nhi);
    printf("%-8s %10s %6s %6s %6s %6s %10s\n", "N", "mean", "p50", "p99", "p999", "max", "fail");
    for (int N : Ns) {
        std::atomic<long> next{0};
        Stat global; std::mutex mu; std::atomic<long> fails{0};
        vector<thread> ts;
        for (int t = 0; t < threads; ++t) {
            ts.emplace_back([&, t]() {
                Stat loc; long lf = 0;
                for (;;) {
                    long idx = next.fetch_add(1);
                    if (idx >= trials) break;
                    uint32_t extra = 0;
                    uint64_t s = baseSeed + (uint64_t)N * 0x9E3779B1u + (uint64_t)idx * 0xD1B54A32D192ED03ull; // distinct strides: equal strides made seed(N, idx) == seed(N+d, idx-d)
                    int rc = roundtrip(s, N, bb, startMode, loss, &extra, nullptr);
                    if (rc == 0) loc.add_overhead(extra); else lf++;
                }
                std::lock_guard<std::mutex> lk(mu); global.merge(loc); fails += lf;
            });
        }
        for (auto& th : ts) th.join();
        auto pct = [&](double p)->int{ uint64_t tgt=(uint64_t)(p*global.trials),c=0; for(int i=0;i<(int)global.overhead_hist.size();++i){c+=global.overhead_hist[i]; if(c>=tgt) return i;} return 63; };
        double mean = global.trials ? (double)global.overhead_sum/global.trials : 0;
        printf("%-8d %10.4f %6d %6d %6d %6u %10ld\n", N, mean, pct(0.5), pct(0.99), pct(0.999), global.overhead_max, fails.load());
        fflush(stdout);
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
static double seed_mean_par(int N, int bb, int startMode, double loss, int pseed, int dseed, long trials,
                            uint64_t base, int threads) {
    std::atomic<long> idx{0};
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
            long k = idx.fetch_add(128);
            if (k >= trials) break;
            long end = k + 128; if (end > trials) end = trials;
            for (long j = k; j < end; ++j) {
                if (stop.load()) break;
                uint32_t extra = 0;
                int rc;
                if (bb == 0) {
                    rc = matrix_decode_roundtrip(base + (uint64_t)j * 2654435761ull, N, startMode, loss, &extra, nullptr);
                }
                else {
                    wh_set_override(ovrN, -1, pseed, dseed);
                    rc = roundtrip(base + (uint64_t)j * 2654435761ull, N, bb, startMode, loss, &extra, nullptr);
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
        // Loud failure on a bad path or malformed file: the old fscanf loop
        // silently produced an empty/truncated N list and the campaign
        // exited 0 looking like a complete search with no weak N.
        FILE* f = fopen(nfile.c_str(), "r");
        if (!f) { fprintf(stderr, "seedsearch: cannot open --nfile %s\n", nfile.c_str()); return 2; }
        char linebuf[512];
        int lineno = 0;
        while (fgets(linebuf, sizeof(linebuf), f)) {
            ++lineno;
            const char* s = linebuf;
            while (*s == ' ' || *s == '\t') ++s;
            if (*s == '\0' || *s == '\n' || *s == '#') continue; // blank/comment
            char* end = nullptr;
            long n = strtol(s, &end, 10);
            while (*end == ' ' || *end == '\t' || *end == '\r' || *end == '\n') ++end;
            if (end == s || *end != '\0' || n < 2 || n > 64000) {
                fprintf(stderr, "seedsearch: bad N at %s:%d: %s", nfile.c_str(), lineno, linebuf);
                fclose(f);
                return 2;
            }
            Ns.push_back((int)n);
        }
        fclose(f);
        if (Ns.empty()) { fprintf(stderr, "seedsearch: --nfile %s has no N values\n", nfile.c_str()); return 2; }
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
        !std::isfinite(loss) || loss < 0.0 || loss > 0.99 ||
        (startMode != 0 && startMode != 1) ||
        dseeds < 0 || !std::isfinite(goodthr) || goodthr < 0.0)
    {
        fprintf(stderr, "seedsearch: invalid parameters\n");
        return 2;
    }
    for (int N : Ns) {
        if (N < 2 || N > 64000) {
            fprintf(stderr, "seedsearch requires every N in 2..64000\n");
            return 2;
        }
    }
    if (dseeds > 256) dseeds = 256;
    fprintf(stderr,"# seedsearch: threads=%d Ns=%zu bb=%d nseeds=%d dseeds=%d tsearch=%ld tverify=%ld goodthr=%.3f\n",
            threads,Ns.size(),bb,nseeds,dseeds,tsearch,tverify,goodthr);
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
        double defMean = seed_mean_par(N, bb, startMode, loss, -1, -1, tverify, 0xD00D + (uint64_t)N*131, threads);
        double v10 = seed_mean_par(N, bb, startMode, loss, bestP, -1, tverify, 0xBEE5 + (uint64_t)N*577, threads);
        double v30 = seed_mean_par(N, bb, startMode, 0.30, bestP, -1, tverify, 0xF00D + (uint64_t)N*331, threads);

        // Stage 2 (Task5): peel-only insufficient -> joint (peel,dense) search over top-K peels x dense seeds.
        if (dseeds > 0 && (v10 >= goodthr || v30 >= goodthr)) {
            // top-K peels by stage-1 @0.10 mean
            vector<int> order(nseeds); for(int s=0;s<nseeds;++s) order[s]=s;
            std::sort(order.begin(), order.end(), [&](int a,int b){return sm[a]<sm[b];});
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
            double jv10 = seed_mean_par(N, bb, startMode, loss, jp, jd, tverify, 0xAB12 + (uint64_t)N*457, threads);
            double jv30 = seed_mean_par(N, bb, startMode, 0.30, jp, jd, tverify, 0xCD34 + (uint64_t)N*149, threads);
            // keep joint result if it beats peel-only on the harder metric
            if (std::max(jv10,jv30) < std::max(v10,v30)) { bestP=jp; bestD=jd; v10=jv10; v30=jv30; }
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
    long trials = 800;
    int bb = 64, startMode = 0;
    double loss = 0.10, thresh = 0.08;
    uint64_t baseSeed = 0x5CA4;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i],"--nlo")&&i+1<argc) { if (!parse_int_strict(argv[++i], nlo)) return 2; }
        else if (!strcmp(argv[i],"--nhi")&&i+1<argc) { if (!parse_int_strict(argv[++i], nhi)) return 2; }
        else if (!strcmp(argv[i],"--trials")&&i+1<argc) { if (!parse_long_strict(argv[++i], trials)) return 2; }
        else if (!strcmp(argv[i],"--bb")&&i+1<argc) { if (!parse_int_strict(argv[++i], bb)) return 2; }
        else if (!strcmp(argv[i],"--startmode")&&i+1<argc) { if (!parse_int_strict(argv[++i], startMode)) return 2; }
        else if (!strcmp(argv[i],"--loss")&&i+1<argc) { if (!parse_double_strict(argv[++i], loss)) return 2; }
        else if (!strcmp(argv[i],"--thresh")&&i+1<argc) { if (!parse_double_strict(argv[++i], thresh)) return 2; }
        else if (!strcmp(argv[i],"--seed")&&i+1<argc) { if (!parse_u64_strict(argv[++i], baseSeed)) return 2; }
    }
    if (nlo < 2 || nhi > 64000 || nlo > nhi || trials < 1 ||
        !valid_seedcheck_args(nlo, bb, startMode, loss) ||
        thresh < 0.0 || !std::isfinite(thresh))
    {
        fprintf(stderr, "scan: invalid range or parameters\n");
        return 2;
    }
    printf("# scan: threads=%d N=[%d,%d] trials/N=%ld bb=%d startMode=%d loss=%.2f thresh=%.3f\n",
           threads, nlo, nhi, trials, bb, startMode, loss, thresh);
    std::atomic<int> nextN{nlo};
    std::mutex pmu;
    std::atomic<long> scanned{0};
    vector<thread> ts;
    for (int t = 0; t < threads; ++t) {
        ts.emplace_back([&]() {
            for (;;) {
                int N = nextN.fetch_add(1);
                if (N > nhi) break;
                uint64_t osum = 0; uint32_t omax = 0; long fail = 0;
                for (long k = 0; k < trials; ++k) {
                    uint32_t extra = 0;
                    uint64_t s = baseSeed + (uint64_t)N * 0x9E3779B1u + (uint64_t)k * 0xD1B54A32D192ED03ull; // distinct strides: equal strides made seed(N, k) == seed(N+d, k-d)
                    int rc = roundtrip(s, N, bb, startMode, loss, &extra, nullptr);
                    if (rc == 0) { osum += extra; if (extra > omax) omax = extra; } else fail++;
                }
                double mean = (double)osum / trials;
                if (mean > thresh || fail > 0) {
                    std::lock_guard<std::mutex> lk(pmu);
                    printf("WEAK N=%-6d mean=%.4f max=%u fail=%ld\n", N, mean, omax, fail);
                    fflush(stdout);
                }
                long sc = ++scanned;
                if (sc % 4000 == 0) { std::lock_guard<std::mutex> lk(pmu); fprintf(stderr, "  ...scanned %ld N\n", sc); }
            }
        });
    }
    for (auto& th : ts) th.join();
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
    size_t idx = (size_t)(p * (double)(sorted.size() - 1) + 0.5);
    if (idx >= sorted.size()) idx = sorted.size() - 1;
    return sorted[idx];
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

    const char* mode = getenv("WH_PEEL_MODE");
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

int main(int argc, char** argv) {
    if (wirehair_init() != Wirehair_Success) { fprintf(stderr, "wirehair_init failed\n"); return 2; }
    if (argc < 2) {
#ifdef WH_COUNT
#ifdef WH_SEED_KNOBS
        fprintf(stderr, "usage: whx [micro|bench|fuzz|repro|ohead|scan|seedmean|seedsearch|peelstat] [--threads T] [opts]\n");
#else
        fprintf(stderr, "usage: whx [micro|bench|fuzz|repro|ohead|scan|peelstat] [--threads T] [opts]\n");
#endif
        fprintf(stderr, "  bench/count/peelstat accept --N csv and --bb/--bb-list csv for block-count x block-size sweeps\n");
#else
#ifdef WH_SEED_KNOBS
        fprintf(stderr, "usage: whx [micro|bench|fuzz|repro|ohead|scan|seedmean|seedsearch] [--threads T] [opts]\n");
#else
        fprintf(stderr, "usage: whx [micro|bench|fuzz|repro|ohead|scan] [--threads T] [opts]\n");
#endif
        fprintf(stderr, "  bench accepts --N csv and --bb/--bb-list csv for block-count x block-size sweeps\n");
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
