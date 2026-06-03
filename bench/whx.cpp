// whx.cpp - Wirehair eXperiment harness
// Modes:
//   micro  : per-core + aggregate throughput of gf256 bulk kernels (+ correctness self-check)
//   bench  : end-to-end encoder_create/encode/decode/recover MBPS + overhead across an N grid
//   fuzz   : parallel correctness fuzzer (random N, blockBytes, loss patterns) verifying exact recovery
//
// Designed to saturate a many-core machine (default ~80% of hw threads).
// Built by bench/build.sh, linking the wirehair sources directly.

#include "wirehair/wirehair.h"
#include "gf256.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
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
    // uniform int in [lo, hi]
    int range(int lo, int hi) { return lo + (int)(u32() % (uint32_t)(hi - lo + 1)); }
};

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

        // add_mem: x ^= y
        memcpy(ref, x, n);
        for (int i = 0; i < n; ++i) ref[i] ^= y[i];
        gf256_add_mem(x, y, n);
        if (memcmp(x, ref, n) != 0) { fprintf(stderr, "FAIL add_mem n=%d\n", n); ok = false; }

        // addset_mem: z = x ^ y
        for (int i = 0; i < n; ++i) ref[i] = x[i] ^ y[i];
        gf256_addset_mem(z, x, y, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL addset_mem n=%d\n", n); ok = false; }

        // add2_mem: z ^= x ^ y
        for (int i = 0; i < n; ++i) z[i] = (uint8_t)rng.u32();
        for (int i = 0; i < n; ++i) ref[i] = (uint8_t)(z[i] ^ x[i] ^ y[i]);
        gf256_add2_mem(z, x, y, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL add2_mem n=%d\n", n); ok = false; }

        // mul_mem: z = x * c
        for (int i = 0; i < n; ++i) ref[i] = ref_mul(x[i], c);
        gf256_mul_mem(z, x, c, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL mul_mem n=%d c=%u\n", n, c); ok = false; }

        // muladd_mem: z += x * c
        for (int i = 0; i < n; ++i) z[i] = (uint8_t)rng.u32();
        for (int i = 0; i < n; ++i) ref[i] = (uint8_t)(z[i] ^ ref_mul(x[i], c));
        gf256_muladd_mem(z, c, x, n);
        if (memcmp(z, ref, n) != 0) { fprintf(stderr, "FAIL muladd_mem n=%d c=%u\n", n, c); ok = false; }

        free(x); free(y); free(z); free(ref);
    }

    // Exhaustive: every constant c in 0..255, over x = 0..255 (+ unaligned offsets),
    // for mul_mem and muladd_mem. Catches a wrong GFNI matrix for any single constant.
    {
        const int n = 256;
        for (int off = 0; off < 7 && ok; ++off) {
            uint8_t* x = (uint8_t*)xaligned(n + 8);
            uint8_t* z = (uint8_t*)xaligned(n + 8);
            uint8_t* ref = (uint8_t*)xaligned(n + 8);
            uint8_t* xb = x + off;
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
            free(x); free(z); free(ref);
        }
    }
    return ok;
}

// ============================================================================
// micro: throughput
// ============================================================================
enum Kernel { K_ADD, K_ADDSET, K_ADD2, K_MUL, K_MULADD };
static const char* kname[] = {"add_mem", "addset_mem", "add2_mem", "mul_mem", "muladd_mem"};
static const int kKernelCount = 5;

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
        if (!strcmp(argv[i], "--secs") && i + 1 < argc) seconds = atof(argv[++i]);
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
static int roundtrip(uint64_t seed, int N, int blockBytes, int startMode, double lossRate,
                     uint32_t* extra_out, char* failmsg) {
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
    // startMode 0: from blockId 0 (systematic + repair). 1: repair-only from a random offset >= N.
    unsigned blockId = (startMode == 1) ? (unsigned)(N + (rng.u32() % 16)) : 0;
    unsigned needed = 0;
    int rc = 0;
    // Cap on *delivered* blocks. Generous so rare-but-legitimate wirehair overhead tails
    // (pathological erasure subsets, ~1e-5 freq) still complete and get counted. A genuinely
    // broken ablation surfaces as decode-error/mismatch (caught immediately) or never-decode.
    const unsigned maxNeeded = (unsigned)N * 2 + 512;
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

static int cmd_fuzz(int argc, char** argv) {
    int threads = resolve_threads();
    double seconds = 20.0;
    int Nmax = 4000;
    uint64_t baseSeed = 0xABCDEF12345ULL;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--secs") && i + 1 < argc) seconds = atof(argv[++i]);
        else if (!strcmp(argv[i], "--nmax") && i + 1 < argc) Nmax = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) baseSeed = strtoull(argv[++i], nullptr, 0);
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
    uint64_t overhead_sum = 0; uint64_t rounds = 0;
    void merge(const BenchAccum& o){
        enc_create_us+=o.enc_create_us; enc_create_n+=o.enc_create_n;
        encode_us+=o.encode_us; encode_n+=o.encode_n; encode_bytes+=o.encode_bytes;
        decode_us+=o.decode_us; decode_n+=o.decode_n; decode_bytes+=o.decode_bytes;
        recover_us+=o.recover_us; recover_n+=o.recover_n; recover_bytes+=o.recover_bytes;
        overhead_sum+=o.overhead_sum; rounds+=o.rounds;
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
    if (!w.enc) return;
    if (timed) { a.enc_create_us += (t1 - t0) * 1e6; a.enc_create_n++; }

    w.dec = wirehair_decoder_create(w.dec, messageBytes, blockBytes); // reuse
    if (!w.dec) return;

    // measure encode of repair blocks (ids N..N+enc_count-1)
    const int enc_count = 64;
    uint32_t wl = 0;
    double e0 = now_sec();
    for (int i = 0; i < enc_count; ++i) wirehair_encode(w.enc, N + i, w.block.data(), blockBytes, &wl);
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
        wirehair_encode(w.enc, thisId, w.block.data(), blockBytes, &wl);
        needed++;
        double d0 = now_sec();
        WirehairResult dr = wirehair_decode(w.dec, thisId, w.block.data(), wl);
        d_accum += now_sec() - d0;
        if (dr == Wirehair_Success) { ok = true; break; }
        if (dr != Wirehair_NeedMore) break;
    }
    if (timed) {
        a.decode_us += d_accum * 1e6; a.decode_n += needed; a.decode_bytes += (uint64_t)needed * blockBytes;
        a.overhead_sum += (needed >= (unsigned)N) ? (needed - N) : 0; a.rounds++;
    }

    if (ok) {
        double r0 = now_sec();
        WirehairResult rr = wirehair_recover(w.dec, w.decoded.data(), messageBytes);
        double r1 = now_sec();
        if (timed && rr == Wirehair_Success) { a.recover_us += (r1 - r0) * 1e6; a.recover_n++; a.recover_bytes += messageBytes; }
    }
}

static int cmd_bench(int argc, char** argv) {
    int threads = resolve_threads();
    int blockBytes = 1300;
    double lossRate = 0.10;
    int rounds_per_N = 0; // 0 => auto by N
    string nlist = "32,128,512,1024,2048,8192,32000";
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--bb") && i + 1 < argc) blockBytes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--loss") && i + 1 < argc) lossRate = atof(argv[++i]);
        else if (!strcmp(argv[i], "--N") && i + 1 < argc) nlist = argv[++i];
        else if (!strcmp(argv[i], "--rounds") && i + 1 < argc) rounds_per_N = atoi(argv[++i]);
    }
    vector<int> Ns;
    { size_t p = 0; while (p < nlist.size()) { size_t q = nlist.find(',', p); string tok = nlist.substr(p, q==string::npos?string::npos:q-p); if(!tok.empty()) Ns.push_back(atoi(tok.c_str())); if (q==string::npos) break; p = q+1; } }
    printf("# bench: threads=%d bb=%d loss=%.2f\n", threads, blockBytes, lossRate);
    printf("%-8s %14s %14s %14s %14s %10s\n", "N", "create_MBPS", "encode_MBPS", "decode_MBPS", "recover_MBPS", "overhead");
    for (int N : Ns) {
        // total rounds across all threads; scale down for large N to keep runtime bounded
        int total_rounds = rounds_per_N;
        if (total_rounds == 0) {
            if (N <= 128) total_rounds = threads * 40;
            else if (N <= 1024) total_rounds = threads * 16;
            else if (N <= 8192) total_rounds = threads * 4;
            else total_rounds = threads * 2;
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
        printf("%-8d %14.1f %14.1f %14.1f %14.1f %10.4f\n", N, create_MBPS, encode_MBPS, decode_MBPS, recover_MBPS, overhead);
        fflush(stdout);
    }
    return 0;
}

// Deterministic single-trial trace: whx repro <seed> <N> <bb> <startMode> <loss>
static int cmd_repro(int argc, char** argv) {
    if (argc < 5) { fprintf(stderr, "repro <seed> <N> <bb> <startMode> <loss>\n"); return 1; }
    uint64_t seed = strtoull(argv[0], nullptr, 0);
    int N = atoi(argv[1]); int bb = atoi(argv[2]); int startMode = atoi(argv[3]); double loss = atof(argv[4]);
    Rng rng(seed);
    int finalBytes = rng.range(1, bb);
    uint64_t messageBytes = (uint64_t)(N - 1) * bb + finalBytes;
    vector<uint8_t> message(messageBytes), decoded(messageBytes, 0xAB), block(bb);
    for (uint64_t i = 0; i < messageBytes; ++i) message[i] = (uint8_t)rng.u32();
    WirehairCodec enc = wirehair_encoder_create(nullptr, message.data(), messageBytes, bb);
    WirehairCodec dec = wirehair_decoder_create(nullptr, messageBytes, bb);
    printf("repro N=%d bb=%d final=%d msgBytes=%llu startMode=%d loss=%.2f enc=%p dec=%p\n",
           N, bb, finalBytes, (unsigned long long)messageBytes, startMode, loss, (void*)enc, (void*)dec);
    unsigned blockId = (startMode==1)?(unsigned)(N+(rng.u32()%16)):0, needed=0;
    for (;;) {
        if (needed >= (unsigned)N + 512 + (unsigned)N) { printf("GIVING UP needed=%u\n", needed); break; }
        bool drop = (startMode==0) && (rng.unit() < loss);
        unsigned thisId = blockId++;
        if (drop) continue;
        uint32_t wl=0;
        wirehair_encode(enc, thisId, block.data(), bb, &wl);
        needed++;
        WirehairResult dr = wirehair_decode(dec, thisId, block.data(), wl);
        if (needed <= 8 || needed % 256 == 0 || dr != Wirehair_NeedMore)
            printf("  id=%u needed=%u dr=%d (%s)\n", thisId, needed, dr, wirehair_result_string(dr));
        if (dr == Wirehair_Success) {
            WirehairResult rr = wirehair_recover(dec, decoded.data(), messageBytes);
            int cmp = memcmp(decoded.data(), message.data(), messageBytes);
            printf("DECODED at needed=%u (extra=%d) recover=%d cmp=%d\n", needed, (int)needed-N, rr, cmp);
            break;
        }
        if (dr != Wirehair_NeedMore) { printf("DECODE ERROR dr=%d\n", dr); break; }
    }
    wirehair_free(enc); wirehair_free(dec);
    return 0;
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
        if (!strcmp(argv[i], "--nlo") && i+1<argc) nlo = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--nhi") && i+1<argc) nhi = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--nstep") && i+1<argc) nstep = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--trials") && i+1<argc) trials = atol(argv[++i]);
        else if (!strcmp(argv[i], "--bb") && i+1<argc) bb = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--startmode") && i+1<argc) startMode = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--loss") && i+1<argc) loss = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i+1<argc) baseSeed = strtoull(argv[++i],nullptr,0);
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
                    uint64_t s = baseSeed + (uint64_t)N * 0x9E3779B1u + (uint64_t)idx * 2654435761ull;
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
    for (long k = 0; k < trials; ++k) {
        wh_set_override(ovrN, -1, pseed, dseed);
        uint32_t extra = 0;
        int rc = roundtrip(rng.next(), N, bb, startMode, loss, &extra, nullptr);
        if (rc == 0) osum += extra; else fail++;
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
    const int ovrN = (pseed < 0 && dseed < 0) ? -1 : N;
    vector<thread> ts;
    for (int t = 0; t < threads; ++t) ts.emplace_back([&]() {
        uint64_t lsum = 0; long lf = 0;
        for (;;) {
            long k = idx.fetch_add(128);
            if (k >= trials) break;
            long end = k + 128; if (end > trials) end = trials;
            for (long j = k; j < end; ++j) {
                wh_set_override(ovrN, -1, pseed, dseed);
                uint32_t extra = 0;
                int rc = roundtrip(base + (uint64_t)j * 2654435761ull, N, bb, startMode, loss, &extra, nullptr);
                if (rc == 0) lsum += extra; else lf++;
            }
        }
        osum += lsum; fail += lf;
        wh_set_override(-1, -1, -1, -1);
    });
    for (auto& th : ts) th.join();
    return fail.load() ? 1e9 : (double)osum.load() / trials;
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
        else if (!strcmp(argv[i],"--tsearch")&&i+1<argc) tsearch=atol(argv[++i]);
        else if (!strcmp(argv[i],"--tverify")&&i+1<argc) tverify=atol(argv[++i]);
        else if (!strcmp(argv[i],"--nseeds")&&i+1<argc) nseeds=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--bb")&&i+1<argc) bb=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--loss")&&i+1<argc) loss=atof(argv[++i]);
        else if (!strcmp(argv[i],"--startmode")&&i+1<argc) startMode=atoi(argv[++i]);
    }
    if (nseeds < 1) nseeds = 1; if (nseeds > 256) nseeds = 256;
    vector<int> Ns;
    if (!nfile.empty()) { FILE* f=fopen(nfile.c_str(),"r"); if(f){ int n; while(fscanf(f,"%d",&n)==1) Ns.push_back(n); fclose(f);} }
    if (!nlist.empty()) { size_t p=0; while(p<nlist.size()){size_t q=nlist.find(',',p); int n=atoi(nlist.substr(p,q==string::npos?string::npos:q-p).c_str()); if(n>0)Ns.push_back(n); if(q==string::npos)break; p=q+1;} }
    int dseeds = 0;          // 0 = peel-only (back-compat); >0 = enable joint dense-seed search for hard N
    double goodthr = 0.05;   // peel-only is accepted if max(v10,v30) < goodthr; else search dense
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i],"--dseeds")&&i+1<argc) dseeds=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--goodthr")&&i+1<argc) goodthr=atof(argv[++i]);
    }
    if (dseeds > 256) dseeds = 256;
    fprintf(stderr,"# seedsearch: threads=%d Ns=%zu nseeds=%d dseeds=%d tsearch=%ld tverify=%ld goodthr=%.3f\n",
            threads,Ns.size(),nseeds,dseeds,tsearch,tverify,goodthr);
    printf("# N  best_pseed  best_dseed  default_mean  best_mean@0.10  best_mean@0.30\n");
    for (int N : Ns) {
        // Stage 1: search peel seeds (default dense). One thread per candidate seed.
        vector<double> sm(256, 1e9);
        std::atomic<int> si{0};
        vector<thread> ts;
        for (int t=0;t<threads;++t) ts.emplace_back([&](){
            for(;;){ int s=si.fetch_add(1); if(s>=nseeds) break;
                sm[s]=seed_mean(N,bb,startMode,loss,s,-1,tsearch,0x5EED + (uint64_t)N*911 + (uint64_t)s*7);
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
                    cs[c]=seed_mean(N,bb,startMode,0.30,p,d,tsearch,0xDEED + (uint64_t)N*701 + (uint64_t)c*13);
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
        if (!strcmp(argv[i],"--nlo")&&i+1<argc) nlo=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--nhi")&&i+1<argc) nhi=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--trials")&&i+1<argc) trials=atol(argv[++i]);
        else if (!strcmp(argv[i],"--bb")&&i+1<argc) bb=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--startmode")&&i+1<argc) startMode=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--loss")&&i+1<argc) loss=atof(argv[++i]);
        else if (!strcmp(argv[i],"--thresh")&&i+1<argc) thresh=atof(argv[++i]);
        else if (!strcmp(argv[i],"--seed")&&i+1<argc) baseSeed=strtoull(argv[++i],nullptr,0);
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
                    uint64_t s = baseSeed + (uint64_t)N * 0x9E3779B1u + (uint64_t)k * 2654435761ull;
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
extern "C" { void wh_stage_reset(); int wh_stage_count(); uint64_t wh_stage_bytes(int); const char* wh_stage_label(int); }
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
}
static int cmd_count(int argc, char** argv) {
    int N=2048, bb=1300; double loss=0.10; int startMode=0;
    for (int i=0;i<argc;++i){ if(!strcmp(argv[i],"--N")&&i+1<argc)N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--bb")&&i+1<argc)bb=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--loss")&&i+1<argc)loss=atof(argv[++i]);
        else if(!strcmp(argv[i],"--startmode")&&i+1<argc)startMode=atoi(argv[++i]); }
    Rng rng(0xC007);
    uint64_t msgBytes=(uint64_t)N*bb;
    vector<uint8_t> msg(msgBytes), dec(msgBytes); vector<uint8_t> blk(bb);
    for(uint64_t i=0;i<msgBytes;++i) msg[i]=(uint8_t)rng.u32();
    printf("# count: N=%d bb=%d loss=%.2f startMode=%d\n", N,bb,loss,startMode);
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
    for(;;){ if(needed>=(unsigned)N*2+512)break;
        bool drop=(startMode==0)&&(rng.unit()<loss); unsigned t=id++; if(drop)continue;
        wirehair_encode(enc,t,blk.data(),bb,&wl); needed++;
        WirehairResult dr=wirehair_decode(dcd,t,blk.data(),wl);
        if(dr==Wirehair_Success){ok=true;break;} if(dr!=Wirehair_NeedMore)break; }
    count_dump("decode(feed)", msgBytes);
    if(ok){ gf256_count_reset(); wirehair_recover(dcd,dec.data(),msgBytes); count_dump("recover", msgBytes);
        printf("# recover correct: %s\n", memcmp(dec.data(),msg.data(),msgBytes)==0?"YES":"NO"); }
    wirehair_free(enc); wirehair_free(dcd);
    return 0;
}
#endif

int main(int argc, char** argv) {
    if (wirehair_init() != Wirehair_Success) { fprintf(stderr, "wirehair_init failed\n"); return 2; }
#ifdef WH_COUNT
    if (argc>=2 && !strcmp(argv[1],"count")) { vector<char*> r(argv+2,argv+argc); return cmd_count((int)r.size(), r.data()); }
#endif
    if (argc < 2) { fprintf(stderr, "usage: whx [micro|bench|fuzz|repro|ohead|scan|seedsearch] [--threads T] [opts]\n"); return 1; }
    string mode = argv[1];
    // Parse the global --threads for ALL modes first (previously scan/seedsearch/repro were
    // dispatched before this, silently ignoring --threads and defaulting to ~102 threads).
    vector<char*> rest;
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "--threads") && i + 1 < argc) { g_threads = atoi(argv[++i]); }
        else rest.push_back(argv[i]);
    }
    int ac = (int)rest.size(); char** av = rest.data();
    if (mode == "repro") return cmd_repro(ac, av);
    if (mode == "scan")  return cmd_scan(ac, av);
#ifdef WH_SEED_KNOBS
    if (mode == "seedsearch") return cmd_seedsearch(ac, av);
#endif
    if (mode == "micro") return cmd_micro(ac, av);
    if (mode == "bench") return cmd_bench(ac, av);
    if (mode == "fuzz")  return cmd_fuzz(ac, av);
    if (mode == "ohead") return cmd_ohead(ac, av);
    fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 1;
}
