#!/usr/bin/env python3
"""Untuned Pareto architecture screen (beads wirehair-g8iv).

Canonical UNTUNED architecture-comparison protocol across WH2 precode
candidates and legacy Wirehair1:

* Construction seeds are UNTUNED and uniform: a fixed seed value is applied
  identically for every K.  There are no per-K seed tables and no
  seed-attempt rescues - exactly one construction attempt per cell.  The WH2
  arm enforces this with `wirehair_v2_bench precodefail --construction-seed`
  (seed_attempt receipts are audited to be 0); the WH1 arm uses
  `whx untunedcells --cseed` (Codec::OverrideSeeds, single ChooseMatrix).
* Shared sealed stratified sample: the K axis 2..64000 is split into six
  operating regions.  Small-K regions are sampled exhaustively with more
  construction seeds; large-K regions are ~10% sampled with fewer seeds.
  The identical (K, construction-seed, loss-stream) cells are used for
  EVERY method and tuning point; everything is drawn from one recorded
  sample seed with an in-file SplitMix64 so the sample is sealed and
  language/version stable.
* EXACT per-cell first-success overhead, capped at
  max(cap_min=8, ceil(50% * K)) delivered extras; beyond the cap the cell
  is censored (censored cells are reported separately and never folded
  into means).  WH1 feeds the delivered iid stream one symbol at a time
  (incremental DecodeFeed) and records extra = received - K exactly.  WH2
  structure-only trials cannot feed incrementally in the mixed completion,
  so precodefail evaluates a per-K overhead ladder: dense levels 0..8 then
  v += ceil(v/2) geometric levels to the cap, with ascending early stop;
  the first-success level is the recorded extra (exact in the dense range,
  a bracketing upper bound in the geometric range, with the last failing
  level recorded as bracket_low).
* Zero-byte/minimal payloads: WH2 runs structure-only solves
  (full_payload_solve=0, probe bytes 1/2); WH1 uses matrix-only feeds.
* Metrics per (method, tuning point, region) are FLOAT aggregates of the
  overhead fraction extra/K: mean over decoded cells, p50 and p99 over the
  full empirical distribution across all cells (censored cells enter as
  +infinity; a percentile landing in censored mass reports null), the
  censored fraction and the encoder-construct-failure fraction, plus CDF
  readouts (success rate at 1%/10%/25%/50% received overhead).  Raw
  per-cell [K, extra, bracket_low] triples are kept in the result JSON so
  any percentile or CDF point is recomputable.  Secondary work counters
  (WH2: solver block XOR/muladd counts; WH1: gf256 byte counters at bb=1
  from a WH_COUNT build, as a documented proxy) and measured encoder
  throughput ride along.  The Pareto plot y-axes are y1 = p99 overhead
  fraction and y2 = mean overhead fraction.
* Encoder throughput = wall time to initialize the encoder and produce the
  FIRST K systematic symbols only (no repair symbols, no decoding), at the
  minimal legal block size (WH2 mixed: 2 bytes, WH2 certified: 1 byte,
  WH1: 1 byte), single-threaded, median of >= 3 reps per K over a
  documented deterministic subsample of each region's sealed K list.
  The timing pass ranks structural per-symbol cost, not large-payload
  memory-bound throughput.  Timing uses the production seed derivation
  (tuning knobs applied); the untuned fixed seed governs only the overhead
  statistics.

Subcommands:
  seal     --out DIR --sample-seed 0xHEX [--smoke]
  run      --out DIR --bench PATH --whx PATH [--whx-count PATH]
           [--workers N] [--methods a,b,...] [--timing-reps N]
           [--max-tasks N] [--skip-timing] [--force-owner]
  collect  --out DIR [--allow-partial]
  methods  (print the tuning-point table)
  selftest

Results: one canonical JSON per method under DIR/results/, rows keyed by
region, carrying hyperparameters, per-overhead success histograms, work
stats, and throughput.  Raw task outputs are kept under DIR/raw/ and the
run is resumable at file granularity.
"""

import argparse
import datetime
import hashlib
import json
import os
import socket
import subprocess
import sys
import tempfile
import threading
import time

MASK64 = (1 << 64) - 1

SCHEMA_MANIFEST = "wirehair.wh2_untuned_pareto.manifest.v2"
SCHEMA_RESULT = "wirehair.wh2_untuned_pareto.result.v2"
SCHEMA_OWNER = "wirehair.environment_owner.v1"

OWNER_MARKER_PATH = "/tmp/wirehair.environment_owner.v1.json"

# (name, lo, hi, sampled_k_count or None for ALL, construction_seed_count)
FULL_REGIONS = [
    ("r1_2_100", 2, 100, None, 10),
    ("r2_100_1000", 101, 1000, None, 5),
    ("r3_1000_5000", 1001, 5000, 400, 3),
    ("r4_5000_10000", 5001, 10000, 500, 3),
    ("r5_10000_20000", 10001, 20000, 1000, 3),
    ("r6_20000_64000", 20001, 64000, 4400, 3),
]

# Same six operating regions with tiny counts: 200 cells per tuning point.
SMOKE_REGIONS = [
    ("r1_2_100", 2, 100, 15, 2),
    ("r2_100_1000", 101, 1000, 15, 2),
    ("r3_1000_5000", 1001, 5000, 10, 2),
    ("r4_5000_10000", 5001, 10000, 10, 2),
    ("r5_10000_20000", 10001, 20000, 15, 2),
    ("r6_20000_64000", 20001, 64000, 35, 2),
]

PROTOCOL = {
    "loss": 0.1,
    "schedule": "iid",
    "pair_bb": 2,
    "overhead_cap_pct": 50,
    "overhead_cap_min": 8,
    "success_rate_readout_pcts": [1, 10, 25, 50],
    "trials_per_cell": 1,
    "trial_index": 0,
    "paired_overhead_stream": True,
    "overhead_early_stop": True,
    "construction_attempts": 1,
    "structure_only": True,
    "wh1_measurement":
        "exact first-success extra via incremental matrix-only DecodeFeed, "
        "censored at max(cap_min, ceil(cap_pct% * K)) extras",
    "wh2_measurement":
        "per-K overhead ladder: dense levels 0..cap_min then v += "
        "ceil(v/2) geometric levels to the cap, ascending early stop; the "
        "first-success level is the recorded extra (exact in the dense "
        "range, a bracketing upper bound above it; the last failing level "
        "is recorded as bracket_low)",
    "percentile_policy":
        "mean over decoded cells only; p50/p99 over all cells with "
        "censored as +infinity (null when the rank lands in censored "
        "mass); success_rate_at_pct counts cells with extra <= "
        "ceil(q% * K) over all cells",
    "wh2_solve_block_bytes": {"certified": 1, "mixed": 2},
    "wh1_block_bytes": 1,
    "timing_block_bytes": {"wh2-certified": 1, "wh2-mixed": 2, "wh1": 1},
    "timing_seed_policy":
        "production seed derivation (tuning knobs applied); the untuned "
        "fixed construction seed governs only the overhead statistics",
    "wh1_work_proxy":
        "gf256 byte counters (add/add2/addset vs mul/muladd) around encoder "
        "init + first K symbols at bb=1 from a WH_COUNT build; bytes at "
        "bb=1 approximate operation counts",
}


def ceil_div(numerator, denominator):
    return -(-numerator // denominator)


def pct_threshold(k, q_pct):
    """Extra symbols for a q% received-overhead threshold: ceil(q% * K)."""
    return ceil_div(q_pct * k, 100)


def overhead_cap(k, cap_pct, cap_min):
    return max(cap_min, pct_threshold(k, cap_pct))


def overhead_ladder(k, cap_pct, cap_min):
    """Per-K WH2 ladder; mirrors wirehair_v2_bench --overhead-ladder-pct."""
    cap = overhead_cap(k, cap_pct, cap_min)
    levels = list(range(0, min(cap, cap_min) + 1))
    level = cap_min
    while level < cap:
        level += (level + 1) // 2
        levels.append(min(level, cap))
    return levels

TIMING_K_PER_REGION_FULL = 24
TIMING_K_PER_REGION_SMOKE = 4

# Chunking targets: keep single tasks in the seconds-to-minutes range while
# leaving many resumable files.
CHUNK_K_BUDGET = 1500000
CHUNK_MAX_VALUES = 400


class SplitMix64(object):
    """Deterministic PRNG matching the C++ splitmix64 used by the benches."""

    def __init__(self, seed):
        self.state = seed & MASK64

    def next_u64(self):
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK64
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
        return (z ^ (z >> 31)) & MASK64

    def below(self, bound):
        """Uniform integer in [0, bound) by rejection (no modulo bias)."""
        if bound <= 0:
            raise ValueError("bound must be positive")
        limit = (MASK64 + 1) - ((MASK64 + 1) % bound)
        while True:
            value = self.next_u64()
            if value < limit:
                return value % bound


def sample_without_replacement(rng, lo, hi, count):
    """Deterministic partial Fisher-Yates sample of [lo, hi], sorted."""
    values = list(range(lo, hi + 1))
    total = len(values)
    if count is None or count >= total:
        return values
    for i in range(count):
        j = i + rng.below(total - i)
        values[i], values[j] = values[j], values[i]
    return sorted(values[:count])


def canonical_json(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n"


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# Tuning-point grid


def _mixed_flags(period=None, geometry=None, gf256_rows=None, gf16_rows=None,
                 schedule=None, hash_seed=None, hash_keyed=False,
                 independent_ext=False, ext_seed_xor=None, skew=None,
                 mix_count=None, dense_rows=None, source_hits=None):
    flags = ["--completion", "mixed"]
    hyper = {"completion": "mixed", "heavy_family": "periodic"}
    if geometry is not None:
        flags += ["--mixed-geometry", geometry]
    hyper["geometry"] = geometry if geometry is not None else "frozen"
    if gf16_rows is not None:
        flags += ["--mixed-gf16-rows", str(gf16_rows)]
    hyper["gf16_rows"] = gf16_rows if gf16_rows is not None else 2
    if period is not None:
        flags += ["--mixed-period", str(period)]
    hyper["period"] = period if period is not None else 244
    if gf256_rows is not None:
        flags += ["--mixed-gf256-rows", str(gf256_rows)]
    hyper["gf256_rows"] = gf256_rows if gf256_rows is not None else 10
    if skew is not None:
        flags += ["--mixed-residue-skew", str(skew)]
    hyper["residue_skew"] = skew if skew is not None else 0
    if schedule is not None:
        flags += ["--mixed-residue-schedule", schedule]
    hyper["residue_schedule"] = schedule if schedule is not None else "constant"
    if hash_seed is not None:
        flags += ["--mixed-residue-hash-seed", str(hash_seed)]
    hyper["residue_hash_seed"] = hash_seed
    if hash_keyed:
        flags += ["--mixed-residue-hash-keyed"]
    hyper["residue_hash_keyed"] = bool(hash_keyed)
    if independent_ext:
        flags += ["--mixed-independent-extension-residues"]
    hyper["independent_extension_residues"] = bool(independent_ext)
    if ext_seed_xor is not None:
        flags += ["--mixed-extension-residue-seed-xor", str(ext_seed_xor)]
    hyper["extension_residue_seed_xor"] = ext_seed_xor
    if dense_rows is not None:
        flags += ["--binary-dense-rows", str(dense_rows)]
    hyper["binary_dense_rows"] = dense_rows if dense_rows is not None else 12
    if source_hits is not None:
        flags += ["--source-hits", str(source_hits)]
    hyper["source_hits"] = source_hits
    if mix_count is not None:
        flags += ["--mix-count", str(mix_count)]
    hyper["mix_count"] = mix_count if mix_count is not None else 3
    return flags, hyper


def _certified_flags(heavy_family=None, heavy_rows=None, dense_rows=None,
                     mix_count=None, source_hits=None):
    flags = ["--completion", "certified"]
    hyper = {"completion": "certified", "geometry": None}
    if heavy_family is not None:
        flags += ["--heavy-family", heavy_family]
    hyper["heavy_family"] = heavy_family if heavy_family else "periodic"
    if heavy_rows is not None:
        flags += ["--gf256-heavy-rows", str(heavy_rows)]
    hyper["gf256_heavy_rows"] = heavy_rows if heavy_rows is not None else 12
    if dense_rows is not None:
        flags += ["--binary-dense-rows", str(dense_rows)]
    hyper["binary_dense_rows"] = dense_rows if dense_rows is not None else 12
    if source_hits is not None:
        flags += ["--source-hits", str(source_hits)]
    hyper["source_hits"] = source_hits
    if mix_count is not None:
        flags += ["--mix-count", str(mix_count)]
    hyper["mix_count"] = mix_count if mix_count is not None else 3
    return flags, hyper


def build_methods():
    """First-screen tuning points (<= 100), biased toward families with
    sealed prior evidence (mixed H15 lineage).  Grouped-GF256 (P x r),
    degree-balanced staircase, and two-anchor variants live on other
    branches and are listed by the report as pending arms."""
    methods = []

    def add(method_id, family, arm, flags, hyper):
        methods.append({
            "id": method_id,
            "family": family,
            "arm": arm,
            "flags": flags,
            "hyperparameters": hyper,
        })

    add("wh1-legacy", "wh1", "wh1", [], {
        "codec": "legacy wirehair.cpp lineage",
        "dense_count": "GetDenseCount(K) table (structural)",
        "construction_seed_mapping":
            "pseed=cseed[15:0], dseed=cseed[31:16]",
    })
    # Tuned control arm: the production per-K seed tables (normal
    # ChooseMatrix path, no override), same cells and same delivered-id
    # streams, quantifying what WH1's per-K seed tuning buys.
    methods[-1]["tuned"] = False
    add("wh1-legacy-tuned", "wh1", "wh1", [], {
        "codec": "legacy wirehair.cpp lineage",
        "dense_count": "GetDenseCount(K) table (structural)",
        "seed_source": "production per-K seed tables (tuned control)",
    })
    methods[-1]["tuned"] = True

    flags, hyper = _certified_flags()
    add("cert-baseline", "wh2-certified", "wh2", flags, hyper)
    for heavy in (6, 18, 24):
        flags, hyper = _certified_flags(heavy_rows=heavy)
        add("cert-h%d" % heavy, "wh2-certified", "wh2", flags, hyper)
    for dense in (8, 16, 24):
        flags, hyper = _certified_flags(dense_rows=dense)
        add("cert-d%d" % dense, "wh2-certified", "wh2", flags, hyper)
    flags, hyper = _certified_flags(dense_rows=16, heavy_rows=18)
    add("cert-d16-h18", "wh2-certified", "wh2", flags, hyper)
    flags, hyper = _certified_flags(heavy_family="hashed")
    add("cert-hashed", "wh2-certified", "wh2", flags, hyper)
    flags, hyper = _certified_flags(mix_count=2)
    add("cert-mix2", "wh2-certified", "wh2", flags, hyper)
    for hits in (2, 4):
        flags, hyper = _certified_flags(source_hits=hits)
        add("cert-sh%d" % hits, "wh2-certified", "wh2", flags, hyper)

    flags, hyper = _mixed_flags()
    add("mixed-p244-baseline", "wh2-mixed-frozen", "wh2", flags, hyper)
    flags, hyper = _mixed_flags(mix_count=2)
    add("mixed-p244-mix2", "wh2-mixed-frozen", "wh2", flags, hyper)
    for period in (16, 24, 32, 48, 64, 96, 128):
        flags, hyper = _mixed_flags(period=period)
        add("mixed-p%d-frozen" % period, "wh2-mixed-frozen", "wh2",
            flags, hyper)

    for period in (24, 32, 48, 64, 96, 128, 244):
        flags, hyper = _mixed_flags(period=period, geometry="shared-x")
        add("mixed-p%d-sharedx" % period, "wh2-mixed-sharedx", "wh2",
            flags, hyper)
    flags, hyper = _mixed_flags(period=32, geometry="shared-x", mix_count=2)
    add("mixed-p32-sharedx-mix2", "wh2-mixed-sharedx", "wh2", flags, hyper)

    for period in (32, 48, 64, 96):
        flags, hyper = _mixed_flags(
            period=period, geometry="shared-x", gf256_rows=11)
        add("mixed-p%d-11x2" % period, "wh2-mixed-11rows", "wh2",
            flags, hyper)

    for period in (31, 32):
        flags, hyper = _mixed_flags(
            period=period, geometry="shared-x", gf256_rows=12, gf16_rows=4)
        add("mixed-p%d-12x4" % period, "wh2-mixed-12rows", "wh2",
            flags, hyper)

    for period in (32, 48, 64):
        flags, hyper = _mixed_flags(
            period=period, geometry="shared-x", gf256_rows=11, gf16_rows=4)
        add("mixed-p%d-11x4" % period, "wh2-mixed-11x4", "wh2", flags, hyper)
    flags, hyper = _mixed_flags(period=32, gf16_rows=4)
    add("mixed-p32-10x4", "wh2-mixed-11x4", "wh2", flags, hyper)

    for period in (32, 48, 64, 96):
        flags, hyper = _mixed_flags(
            period=period, geometry="shared-x", gf256_rows=11, gf16_rows=4,
            schedule="hashed", hash_seed=68, hash_keyed=True,
            independent_ext=True, ext_seed_xor=78, mix_count=2)
        add("h15-p%d" % period, "wh2-mixed-h15", "wh2", flags, hyper)
    flags, hyper = _mixed_flags(
        period=32, geometry="shared-x", gf256_rows=11, gf16_rows=4,
        schedule="hashed", hash_seed=68, hash_keyed=True,
        independent_ext=True, ext_seed_xor=78, mix_count=3)
    add("h15-p32-mix3", "wh2-mixed-h15", "wh2", flags, hyper)
    flags, hyper = _mixed_flags(
        period=32, geometry="shared-x", gf256_rows=11, gf16_rows=4,
        schedule="hashed", hash_seed=68, hash_keyed=True, mix_count=2)
    add("h15-p32-sharedext", "wh2-mixed-h15", "wh2", flags, hyper)

    for skew in (8, 16):
        flags, hyper = _mixed_flags(period=48, geometry="shared-x", skew=skew)
        add("mixed-p48-skew%d" % skew, "wh2-mixed-skew", "wh2", flags, hyper)

    for period in (48, 96):
        flags, hyper = _mixed_flags(
            period=period, geometry="shared-x", schedule="ramp")
        add("mixed-p%d-ramp" % period, "wh2-mixed-ramp", "wh2", flags, hyper)

    for dense in (8, 16):
        flags, hyper = _mixed_flags(dense_rows=dense)
        add("mixed-p244-d%d" % dense, "wh2-mixed-frozen", "wh2", flags, hyper)

    ids = [m["id"] for m in methods]
    if len(ids) != len(set(ids)):
        raise RuntimeError("duplicate method ids in tuning grid")
    if len(ids) > 100:
        raise RuntimeError("tuning grid exceeds the 100-point screen cap")
    return methods


# ---------------------------------------------------------------------------
# seal


def build_sample(sample_seed, smoke):
    regions_spec = SMOKE_REGIONS if smoke else FULL_REGIONS
    timing_per_region = (TIMING_K_PER_REGION_SMOKE if smoke
                         else TIMING_K_PER_REGION_FULL)
    regions = []
    for index, (name, lo, hi, count, seed_count) in enumerate(regions_spec):
        # Independent, recorded stream per region purpose so a change in one
        # region's spec cannot silently reshuffle another region.
        k_rng = SplitMix64(sample_seed ^ (0xA0 + index))
        seed_rng = SplitMix64(sample_seed ^ (0xB0 + index))
        loss_rng = SplitMix64(sample_seed ^ (0xC0 + index))
        timing_rng = SplitMix64(sample_seed ^ (0xD0 + index))
        k_values = sample_without_replacement(k_rng, lo, hi, count)
        construction_seeds = []
        for _ in range(seed_count):
            value = 0
            while value == 0:
                value = seed_rng.next_u64()
            construction_seeds.append(value)
        loss_seeds = [loss_rng.next_u64() for _ in range(seed_count)]
        timing_k = sample_without_replacement(
            SplitMix64(timing_rng.next_u64()), 0, len(k_values) - 1,
            min(timing_per_region, len(k_values)))
        regions.append({
            "name": name,
            "k_lo": lo,
            "k_hi": hi,
            "k_values": k_values,
            "construction_seeds": construction_seeds,
            "loss_seeds": loss_seeds,
            "timing_k_values": [k_values[i] for i in timing_k],
        })
    return regions


def cmd_seal(args):
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_dir, "manifest.json")
    if os.path.exists(manifest_path) and not args.force:
        print("seal: %s already exists (use --force to reseal)" %
              manifest_path, file=sys.stderr)
        return 1
    sample_seed = int(args.sample_seed, 0)
    regions = build_sample(sample_seed, args.smoke)
    methods = build_methods()
    manifest = {
        "schema": SCHEMA_MANIFEST,
        "sample_seed": "0x%x" % sample_seed,
        "smoke": bool(args.smoke),
        "protocol": PROTOCOL,
        "regions": regions,
        "methods": methods,
        "sealed_utc": utc_now(),
    }
    text = canonical_json(manifest)
    with open(manifest_path, "w") as handle:
        handle.write(text)
    with open(os.path.join(out_dir, "manifest.sha256"), "w") as handle:
        handle.write(sha256_text(text) + "\n")
    total = sum(len(r["k_values"]) * len(r["construction_seeds"])
                for r in regions)
    print("sealed %s: %d regions, %d cells/tuning-point, %d methods, "
          "sha256=%s" % (manifest_path, len(regions), total, len(methods),
                         sha256_text(text)))
    return 0


def load_manifest(out_dir):
    manifest_path = os.path.join(out_dir, "manifest.json")
    with open(manifest_path, "r") as handle:
        text = handle.read()
    recorded_path = os.path.join(out_dir, "manifest.sha256")
    if os.path.exists(recorded_path):
        with open(recorded_path, "r") as handle:
            recorded = handle.read().strip()
        if recorded and recorded != sha256_text(text):
            raise RuntimeError("manifest.json does not match manifest.sha256")
    return json.loads(text), sha256_text(text)


# ---------------------------------------------------------------------------
# task construction


def chunk_k_values(k_values):
    chunks = []
    current = []
    budget = 0
    for k in k_values:
        current.append(k)
        budget += k
        if budget >= CHUNK_K_BUDGET or len(current) >= CHUNK_MAX_VALUES:
            chunks.append(current)
            current = []
            budget = 0
    if current:
        chunks.append(current)
    return chunks


def wh2_cell_command(bench, method, k_values, construction_seed, loss_seed,
                     protocol):
    return [
        bench, "precodefail",
        "--N", ",".join(str(k) for k in k_values),
        "--bb-list", str(protocol["pair_bb"]),
        "--overhead-ladder-pct", str(protocol["overhead_cap_pct"]),
        "--overhead-ladder-min", str(protocol["overhead_cap_min"]),
        "--trials", str(protocol["trials_per_cell"]),
        "--threads", "1",
        "--loss", repr(protocol["loss"]),
        "--seed", "0x%x" % loss_seed,
        "--schedule", protocol["schedule"],
        "--paired-overhead-stream",
        "--overhead-early-stop",
        "--construction-seed", "0x%x" % construction_seed,
    ] + list(method["flags"])


def wh2_timing_command(bench, method, k_values, reps, protocol):
    return [
        bench, "precodefail",
        "--N", ",".join(str(k) for k in k_values),
        "--bb-list", str(protocol["pair_bb"]),
        "--overhead", "0",
        "--trials", "1",
        "--threads", "1",
        "--loss", "0",
        "--seed", "1",
        "--encode-timing", str(reps),
    ] + list(method["flags"])


def wh1_cell_command(whx, k_values, construction_seed, loss_seed, protocol,
                     tuned=False):
    command = [
        whx, "untunedcells",
        "--threads", "1",
        "--nlist", ",".join(str(k) for k in k_values),
        "--seed", "0x%x" % loss_seed,
        "--loss", repr(protocol["loss"]),
        "--pairbb", str(protocol["pair_bb"]),
        "--cap-pct", str(protocol["overhead_cap_pct"]),
        "--cap-min", str(protocol["overhead_cap_min"]),
        "--trial", str(protocol["trial_index"]),
    ]
    if tuned:
        command.append("--tuned")
    else:
        command += ["--cseed", "0x%x" % construction_seed]
    return command


def wh1_timing_command(whx, k_values, reps, protocol):
    return [
        whx, "enctime",
        "--nlist", ",".join(str(k) for k in k_values),
        "--reps", str(reps),
        "--bb", str(protocol["wh1_block_bytes"]),
    ]


def build_tasks(manifest, out_dir, bench, whx, whx_count, method_filter,
                timing_reps, skip_timing, region_filter=None):
    protocol = manifest["protocol"]
    tasks = []
    methods = manifest["methods"]
    if method_filter:
        wanted = set(method_filter)
        known = set(m["id"] for m in methods)
        missing = wanted - known
        if missing:
            raise RuntimeError(
                "unknown method ids: %s" % ",".join(sorted(missing)))
        methods = [m for m in methods if m["id"] in wanted]
    regions = manifest["regions"]
    if region_filter:
        wanted_regions = set(region_filter)
        known_regions = set(r["name"] for r in regions)
        missing_regions = wanted_regions - known_regions
        if missing_regions:
            raise RuntimeError(
                "unknown regions: %s" % ",".join(sorted(missing_regions)))
        regions = [r for r in regions if r["name"] in wanted_regions]
    for method in methods:
        for region in regions:
            base = os.path.join(out_dir, "raw", method["id"], region["name"])
            chunks = chunk_k_values(region["k_values"])
            for seed_index, construction_seed in enumerate(
                    region["construction_seeds"]):
                loss_seed = region["loss_seeds"][seed_index]
                for chunk_index, chunk in enumerate(chunks):
                    out_path = os.path.join(
                        base, "cells_cs%d_chunk%04d.csv" %
                        (seed_index, chunk_index))
                    if method["arm"] == "wh2":
                        command = wh2_cell_command(
                            bench, method, chunk, construction_seed,
                            loss_seed, protocol)
                    else:
                        command = wh1_cell_command(
                            whx, chunk, construction_seed, loss_seed,
                            protocol, tuned=method.get("tuned", False))
                    tasks.append({
                        "kind": "cells",
                        "method": method["id"],
                        "region": region["name"],
                        "out": out_path,
                        "command": command,
                        "cost": sum(chunk),
                    })
            if skip_timing:
                continue
            timing_k = region["timing_k_values"]
            if method["arm"] == "wh2":
                command = wh2_timing_command(
                    bench, method, timing_k, timing_reps, protocol)
            else:
                command = wh1_timing_command(
                    whx, timing_k, timing_reps, protocol)
            tasks.append({
                "kind": "timing",
                "method": method["id"],
                "region": region["name"],
                "out": os.path.join(base, "timing.csv"),
                "command": command,
                "cost": sum(timing_k) * timing_reps,
            })
            if method["arm"] == "wh1" and whx_count:
                tasks.append({
                    "kind": "counters",
                    "method": method["id"],
                    "region": region["name"],
                    "out": os.path.join(base, "counters.csv"),
                    "command": wh1_timing_command(
                        whx_count, timing_k, timing_reps, protocol),
                    "cost": sum(timing_k) * timing_reps,
                })
    return tasks


# ---------------------------------------------------------------------------
# environment owner marker


def acquire_owner_marker(out_dir, force):
    marker = {
        "schema": SCHEMA_OWNER,
        "owner": "wh2_untuned_pareto_sweep",
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "started_utc": utc_now(),
        "out_dir": out_dir,
    }
    if os.path.exists(OWNER_MARKER_PATH):
        try:
            with open(OWNER_MARKER_PATH, "r") as handle:
                existing = json.load(handle)
        except (ValueError, OSError):
            existing = {}
        pid = existing.get("pid")
        alive = False
        if isinstance(pid, int) and pid > 0:
            try:
                os.kill(pid, 0)
                alive = True
            except OSError:
                alive = False
        if alive and existing.get("out_dir") != out_dir and not force:
            raise RuntimeError(
                "environment owner marker %s is held by live pid %s "
                "(owner=%s out_dir=%s); pass --force-owner to take over" %
                (OWNER_MARKER_PATH, pid, existing.get("owner"),
                 existing.get("out_dir")))
    with open(OWNER_MARKER_PATH, "w") as handle:
        json.dump(marker, handle, sort_keys=True)
        handle.write("\n")
    return marker


def release_owner_marker(marker):
    try:
        with open(OWNER_MARKER_PATH, "r") as handle:
            existing = json.load(handle)
        if existing.get("pid") == marker["pid"] and \
                existing.get("started_utc") == marker["started_utc"]:
            os.unlink(OWNER_MARKER_PATH)
    except (ValueError, OSError):
        pass


# ---------------------------------------------------------------------------
# run


def run_task(task):
    out_path = task["out"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp_path = out_path + ".tmp"
    err_path = out_path + ".stderr"
    with open(tmp_path, "w") as out_handle, \
            open(err_path, "w") as err_handle:
        result = subprocess.run(
            task["command"], stdout=out_handle, stderr=err_handle)
    if os.path.getsize(err_path) == 0:
        os.unlink(err_path)
    if result.returncode != 0:
        os.rename(tmp_path, out_path + ".failed")
        return False, result.returncode
    os.replace(tmp_path, out_path)
    return True, 0


def cmd_run(args):
    out_dir = os.path.abspath(args.out)
    manifest, manifest_sha = load_manifest(out_dir)
    bench = os.path.abspath(args.bench) if args.bench else None
    whx = os.path.abspath(args.whx) if args.whx else None
    whx_count = os.path.abspath(args.whx_count) if args.whx_count else None
    tasks = build_tasks(
        manifest, out_dir, bench, whx, whx_count,
        args.methods.split(",") if args.methods else None,
        args.timing_reps, args.skip_timing,
        args.regions.split(",") if args.regions else None)
    for task in tasks:
        arm_binary = task["command"][0]
        if arm_binary is None or not os.path.exists(arm_binary):
            print("run: missing binary %s (pass --bench/--whx)" % arm_binary,
                  file=sys.stderr)
            return 1
    pending = [t for t in tasks if not os.path.exists(t["out"])]
    for task in pending:
        failed = task["out"] + ".failed"
        if os.path.exists(failed):
            os.unlink(failed)
    complete = len(tasks) - len(pending)
    if args.max_tasks is not None:
        pending = pending[:args.max_tasks]
    # Large chunks first so the tail of the run is fine-grained.
    pending.sort(key=lambda t: -t["cost"])
    print("run: %d tasks total, %d already complete, %d selected pending, "
          "workers=%d" % (len(tasks), complete, len(pending), args.workers))
    if not pending:
        return 0
    marker = acquire_owner_marker(out_dir, args.force_owner)
    run_manifest = {
        "schema": "wirehair.wh2_untuned_pareto.run.v1",
        "manifest_sha256": manifest_sha,
        "started_utc": utc_now(),
        "host": socket.gethostname(),
        "workers": args.workers,
        "timing_reps": args.timing_reps,
        "binaries": {},
        "git_commit": None,
    }
    for name, path in (("bench", bench), ("whx", whx),
                       ("whx_count", whx_count)):
        if path and os.path.exists(path):
            run_manifest["binaries"][name] = {
                "path": path, "sha256": sha256_file(path)}
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            capture_output=True, text=True)
        if commit.returncode == 0:
            run_manifest["git_commit"] = commit.stdout.strip()
    except OSError:
        pass
    with open(os.path.join(out_dir, "run_manifest.json"), "w") as handle:
        handle.write(canonical_json(run_manifest))

    lock = threading.Lock()
    state = {"next": 0, "done": 0, "failed": 0}
    start = time.time()

    def worker():
        while True:
            with lock:
                index = state["next"]
                if index >= len(pending):
                    return
                state["next"] = index + 1
            task = pending[index]
            ok, code = run_task(task)
            with lock:
                state["done"] += 1
                if not ok:
                    state["failed"] += 1
                    print("run: task failed rc=%d %s" %
                          (code, task["out"]), file=sys.stderr)
                done = state["done"]
            if done % 50 == 0 or done == len(pending):
                elapsed = time.time() - start
                print("run: %d/%d tasks (%.1fs elapsed)" %
                      (done, len(pending), elapsed))

    threads = []
    try:
        for _ in range(max(1, min(args.workers, len(pending)))):
            thread = threading.Thread(target=worker)
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join()
    finally:
        release_owner_marker(marker)
    print("run: complete, %d failed" % state["failed"])
    return 1 if state["failed"] else 0


# ---------------------------------------------------------------------------
# parsing


class ParseError(RuntimeError):
    pass


def parse_wh2_cells(text, expected_cseed, cap_pct, cap_min):
    """Parse a precodefail exact-ladder cell task output.

    Returns {K: {"extra": int or None, "bracket_low": int, "work": ...}}.
    Rows per K must walk the per-K ladder in order; the first success row
    is the recorded extra (upper bound above the dense range) with the
    previous level as bracket_low.  A censored cell must show every ladder
    level failing.  Enforces seed_attempt == 0 on every row and the
    untuned receipt echo.
    """
    cells = {}
    saw_untuned = False
    header_fields = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("# untuned:"):
            saw_untuned = True
            if "construction_seed_explicit=1" not in line:
                raise ParseError("untuned receipt without explicit seed")
            if "construction_seed=0x%x " % expected_cseed not in line + " ":
                raise ParseError(
                    "untuned receipt seed mismatch: %s" % line)
            if "construction_attempts=1" not in line:
                raise ParseError("untuned receipt attempts != 1")
            if "overhead_ladder_pct=%d " % cap_pct not in line + " ":
                raise ParseError("untuned receipt ladder pct mismatch")
            continue
        if line.startswith("#"):
            continue
        if line.startswith("N,bb,"):
            header_fields = line.split(",")
            continue
        if header_fields is None:
            raise ParseError("row before CSV header: %s" % line)
        fields = line.split(",")
        if len(fields) != len(header_fields):
            raise ParseError("row width mismatch: %s" % line)
        row = dict(zip(header_fields, fields))
        k = int(row["N"])
        overhead = int(row["overhead"])
        if int(row["trials"]) != 1:
            raise ParseError("expected trials=1: %s" % line)
        if int(row["seed_attempt"]) != 0:
            raise ParseError(
                "seed_attempt escalation detected (K=%d): %s" % (k, line))
        if int(row["error"]) != 0:
            raise ParseError("solver error at K=%d: %s" % (k, line))
        success = int(row["success"]) == 1
        cell = cells.setdefault(
            k, {"extra": None, "bracket_low": -1, "work": None, "rows": 0})
        ladder = overhead_ladder(k, cap_pct, cap_min)
        if cell["extra"] is not None:
            raise ParseError("row after success at K=%d: %s" % (k, line))
        if cell["rows"] >= len(ladder) or overhead != ladder[cell["rows"]]:
            raise ParseError(
                "K=%d row %d overhead %d does not walk the ladder %r" %
                (k, cell["rows"], overhead, ladder))
        cell["rows"] += 1
        if success:
            cell["extra"] = overhead
            cell["bracket_low"] = (
                ladder[cell["rows"] - 2] if cell["rows"] >= 2 else -1)
            cell["work"] = (float(row["block_xors_mu"]),
                            float(row["block_muladds_mu"]))
    if not saw_untuned:
        raise ParseError("missing untuned receipt")
    for k, cell in cells.items():
        ladder = overhead_ladder(k, cap_pct, cap_min)
        if cell["extra"] is None and cell["rows"] != len(ladder):
            raise ParseError(
                "censored cell K=%d has %d rows, expected %d" %
                (k, cell["rows"], len(ladder)))
        del cell["rows"]
    return cells


def parse_wh1_cells(text, expected_cseed, tuned=False):
    """Parse a whx untunedcells output (exact incremental extras).

    Returns {K: {"extra": int or None, "enc_ok": bool, "bracket_low": ...}}.
    """
    cells = {}
    saw_receipt = False
    saw_header = False
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("# untunedcells:"):
            saw_receipt = True
            expected_arm = "arm=tuned " if tuned else "arm=untuned "
            if expected_arm not in line + " ":
                raise ParseError("untunedcells arm mismatch: %s" % line)
            if not tuned and \
                    "cseed=0x%x " % expected_cseed not in line + " ":
                raise ParseError("untunedcells cseed mismatch: %s" % line)
            if "construction_attempts=1" not in line:
                raise ParseError("untunedcells attempts != 1")
            continue
        if line.startswith("#"):
            continue
        if line.startswith("K,"):
            if line != "K,enc_ok,extra,censored,decode_limit":
                raise ParseError("unexpected untunedcells header: %s" % line)
            saw_header = True
            continue
        fields = line.split(",")
        if len(fields) != 5:
            raise ParseError("bad untunedcells row: %s" % line)
        k, enc_ok, extra, censored, decode_limit = (int(f) for f in fields)
        if (extra < 0) != (censored == 1):
            raise ParseError("inconsistent censoring at K=%d" % k)
        if decode_limit == 1 and censored != 1:
            raise ParseError("decode_limit without censoring at K=%d" % k)
        cells[k] = {
            "extra": None if extra < 0 else extra,
            # Incremental feeding is exact: the previous delivered count
            # failed by construction.
            "bracket_low": -1 if extra <= 0 else extra - 1,
            "enc_ok": enc_ok == 1,
            # WH1's bounded extra-row replacement pool gave up before the
            # cap (Wirehair_ExtraInsufficient): a terminal decode outcome,
            # counted inside the censored fraction and reported separately.
            "decode_limit": decode_limit == 1,
        }
    if not saw_receipt or not saw_header:
        raise ParseError("missing untunedcells receipt/header")
    return cells


def parse_wh2_timing(text):
    """Parse encode_timing receipt lines from precodefail.

    Returns {K: {"median_ns": int, "solve_ok": bool, "block_bytes": int}}.
    """
    timing = {}
    for line in text.splitlines():
        line = line.strip()
        if not line.startswith("# encode_timing,"):
            continue
        entries = {}
        for token in line[2:].split(","):
            if "=" in token:
                key, value = token.split("=", 1)
                entries[key] = value
        k = int(entries["N"])
        timing[k] = {
            "median_ns": int(entries["median_ns"]),
            "solve_ok": entries["solve_ok"] == "1",
            "block_bytes": int(entries["solve_block_bytes"]),
        }
    if not timing:
        raise ParseError("no encode_timing receipts found")
    return timing


def parse_wh1_timing(text):
    """Parse whx enctime output.

    Returns ({K: {...}}, counted) where counted reflects the WH_COUNT build.
    """
    timing = {}
    counted = False
    saw_header = False
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("# enctime:"):
            counted = "counted=1" in line
            continue
        if line.startswith("#"):
            continue
        if line.startswith("K,"):
            saw_header = True
            continue
        fields = line.split(",")
        if len(fields) != 9:
            raise ParseError("bad enctime row: %s" % line)
        k = int(fields[0])
        timing[k] = {
            "median_ns": int(fields[3]),
            "solve_ok": int(fields[2]) == 1,
            "block_bytes": None,
            "xor_bytes": int(fields[6]),
            "muladd_bytes": int(fields[7]),
        }
    if not saw_header or not timing:
        raise ParseError("no enctime rows found")
    return timing, counted


# ---------------------------------------------------------------------------
# aggregation


def percentile_with_censoring(sorted_fractions, censored, p):
    """Nearest-rank percentile over all cells; censored cells sit above
    every finite value.  Returns None when the rank lands in censored
    mass."""
    total = len(sorted_fractions) + censored
    if total == 0:
        return None
    rank = ceil_div(int(round(p * 10000)) * total, 10000)
    if rank < 1:
        rank = 1
    if rank > len(sorted_fractions):
        return None
    return sorted_fractions[rank - 1]


def aggregate_region(cell_outcomes, protocol):
    """cell_outcomes: list of dicts with k, extra (int or None),
    bracket_low, optional work (xors, muladds), optional enc_ok.

    Float schema: mean(extra/K) over decoded cells; p50/p99 over all cells
    with censored as +infinity; success rate at the readout percentages;
    per-cell triples retained for recomputable CDF/percentiles.
    """
    readouts = protocol["success_rate_readout_pcts"]
    construct_failures = 0
    decoder_limit = 0
    censored = 0
    fractions = []
    fraction_sum = 0.0
    success_at = {str(q): 0 for q in readouts}
    cells_detail = []
    work_xors = 0.0
    work_muladds = 0.0
    work_cells = 0
    for outcome in cell_outcomes:
        k = outcome["k"]
        extra = outcome["extra"]
        if outcome.get("enc_ok") is False:
            construct_failures += 1
        if outcome.get("decode_limit"):
            decoder_limit += 1
        cells_detail.append([
            k,
            -1 if extra is None else extra,
            outcome.get("bracket_low", -1),
        ])
        if extra is None:
            censored += 1
            continue
        fraction = extra / k
        fractions.append(fraction)
        fraction_sum += fraction
        for q in readouts:
            if extra <= pct_threshold(k, q):
                success_at[str(q)] += 1
        work = outcome.get("work")
        if work is not None:
            work_xors += work[0]
            work_muladds += work[1]
            work_cells += 1
    cells = len(cell_outcomes)
    decoded = len(fractions)
    fractions.sort()
    return {
        "cells": cells,
        "decoded": decoded,
        "construct_failures": construct_failures,
        "construct_failure_fraction":
            (construct_failures / cells) if cells else 0.0,
        "censored_fraction": (censored / cells) if cells else 0.0,
        "decoder_limit_fraction":
            (decoder_limit / cells) if cells else 0.0,
        "overhead_fraction": {
            "mean_decoded": (fraction_sum / decoded) if decoded else None,
            "p50": percentile_with_censoring(fractions, censored, 0.50),
            "p99": percentile_with_censoring(fractions, censored, 0.99),
        },
        "success_rate_at_pct": {
            q: (success_at[q] / cells) if cells else 0.0
            for q in success_at
        },
        "cells_detail": cells_detail,
        "work": {
            "decoded_cells_with_work": work_cells,
            "block_xors_mu":
                (work_xors / work_cells) if work_cells else None,
            "block_muladds_mu":
                (work_muladds / work_cells) if work_cells else None,
        },
    }


def aggregate_timing(timing_by_k, reps, counted=False):
    ks = sorted(timing_by_k)
    per_k = [[k, timing_by_k[k]["median_ns"]] for k in ks]
    rates = sorted(
        k / (timing_by_k[k]["median_ns"] / 1e9)
        for k in ks if timing_by_k[k]["median_ns"] > 0)
    seconds = [timing_by_k[k]["median_ns"] / 1e9 for k in ks]
    result = {
        "reps": reps,
        "k_count": len(ks),
        "median_symbols_per_sec":
            rates[len(rates) // 2] if rates else None,
        "mean_seconds_per_k_symbols":
            (sum(seconds) / len(seconds)) if seconds else None,
        "all_solves_ok": all(timing_by_k[k]["solve_ok"] for k in ks),
        "per_k_median_ns": per_k,
    }
    if counted:
        result["work_proxy"] = {
            "note": PROTOCOL["wh1_work_proxy"],
            "per_k": [[k, timing_by_k[k]["xor_bytes"],
                       timing_by_k[k]["muladd_bytes"]] for k in ks],
        }
    return result


# ---------------------------------------------------------------------------
# collect


def collect_method(manifest, manifest_sha, out_dir, method, allow_partial):
    protocol = manifest["protocol"]
    cap_pct = protocol["overhead_cap_pct"]
    cap_min = protocol["overhead_cap_min"]
    missing = []
    result_regions = {}
    for region in manifest["regions"]:
        base = os.path.join(out_dir, "raw", method["id"], region["name"])
        chunks = chunk_k_values(region["k_values"])
        outcomes = []
        for seed_index, construction_seed in enumerate(
                region["construction_seeds"]):
            for chunk_index, chunk in enumerate(chunks):
                path = os.path.join(
                    base, "cells_cs%d_chunk%04d.csv" %
                    (seed_index, chunk_index))
                if not os.path.exists(path):
                    missing.append(path)
                    continue
                with open(path, "r") as handle:
                    text = handle.read()
                if method["arm"] == "wh2":
                    cells = parse_wh2_cells(
                        text, construction_seed, cap_pct, cap_min)
                else:
                    cells = parse_wh1_cells(
                        text, construction_seed,
                        tuned=method.get("tuned", False))
                chunk_set = set(chunk)
                if set(cells) != chunk_set:
                    raise ParseError(
                        "%s covers %d K values, expected %d" %
                        (path, len(cells), len(chunk_set)))
                for k in chunk:
                    outcome = cells[k]
                    outcome["k"] = k
                    outcomes.append(outcome)
        entry = aggregate_region(outcomes, protocol)
        timing_path = os.path.join(base, "timing.csv")
        if os.path.exists(timing_path):
            with open(timing_path, "r") as handle:
                text = handle.read()
            if method["arm"] == "wh2":
                timing = parse_wh2_timing(text)
                entry["timing"] = aggregate_timing(
                    timing, reps=None)
                entry["timing"]["block_bytes"] = (
                    2 if "mixed" in method["hyperparameters"].get(
                        "completion", "") else 1)
            else:
                timing, counted = parse_wh1_timing(text)
                counters_path = os.path.join(base, "counters.csv")
                if not counted and os.path.exists(counters_path):
                    with open(counters_path, "r") as handle:
                        counter_timing, counted = parse_wh1_timing(
                            handle.read())
                    if counted:
                        for k, values in counter_timing.items():
                            if k in timing:
                                timing[k]["xor_bytes"] = values["xor_bytes"]
                                timing[k]["muladd_bytes"] = \
                                    values["muladd_bytes"]
                entry["timing"] = aggregate_timing(
                    timing, reps=None, counted=counted)
                entry["timing"]["block_bytes"] = protocol["wh1_block_bytes"]
        else:
            missing.append(timing_path)
        result_regions[region["name"]] = entry
    if missing and not allow_partial:
        raise ParseError(
            "method %s missing %d raw files (first: %s); rerun `run` or "
            "pass --allow-partial" %
            (method["id"], len(missing), missing[0]))
    result = {
        "schema": SCHEMA_RESULT,
        "method": method["id"],
        "family": method["family"],
        "arm": method["arm"],
        "hyperparameters": method["hyperparameters"],
        "flags": method["flags"],
        "manifest_sha256": manifest_sha,
        "protocol": protocol,
        "regions": result_regions,
        "missing_raw_files": len(missing),
        "collected_utc": utc_now(),
    }
    return result


def cmd_collect(args):
    out_dir = os.path.abspath(args.out)
    manifest, manifest_sha = load_manifest(out_dir)
    results_dir = os.path.join(out_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    methods = manifest["methods"]
    if args.methods:
        wanted = set(args.methods.split(","))
        methods = [m for m in methods if m["id"] in wanted]
    index = []
    for method in methods:
        raw_dir = os.path.join(out_dir, "raw", method["id"])
        if not os.path.isdir(raw_dir):
            continue
        result = collect_method(
            manifest, manifest_sha, out_dir, method, args.allow_partial)
        path = os.path.join(results_dir, "%s.json" % method["id"])
        with open(path, "w") as handle:
            handle.write(canonical_json(result))
        index.append(result)
        for region_name in sorted(result["regions"]):
            region = result["regions"][region_name]
            timing = region.get("timing", {})
            rate = timing.get("median_symbols_per_sec")
            fraction = region["overhead_fraction"]

            def fmt(value):
                return "%.6f" % value if value is not None else "cens"

            success = region["success_rate_at_pct"]
            print("%-24s %-16s cells=%-6d mean=%-9s p50=%-9s p99=%-9s "
                  "cens=%-8.5f cfail=%-8.5f sr@1/10/25/50=%s sym/s=%s" % (
                      method["id"], region_name, region["cells"],
                      fmt(fraction["mean_decoded"]),
                      fmt(fraction["p50"]), fmt(fraction["p99"]),
                      region["censored_fraction"],
                      region["construct_failure_fraction"],
                      "/".join("%.4f" % success[str(q)] for q in
                               result["protocol"][
                                   "success_rate_readout_pcts"]),
                      ("%.0f" % rate) if rate else "n/a"))
    with open(os.path.join(results_dir, "index.json"), "w") as handle:
        handle.write(canonical_json({
            "schema": "wirehair.wh2_untuned_pareto.index.v1",
            "manifest_sha256": manifest_sha,
            "methods": [r["method"] for r in index],
            "collected_utc": utc_now(),
        }))
    print("collect: wrote %d method result files to %s" %
          (len(index), results_dir))
    return 0


def cmd_methods(_args):
    methods = build_methods()
    print("%-24s %-20s %-4s flags" % ("id", "family", "arm"))
    for method in methods:
        print("%-24s %-20s %-4s %s" % (
            method["id"], method["family"], method["arm"],
            " ".join(method["flags"])))
    print("total tuning points: %d" % len(methods))
    return 0


# ---------------------------------------------------------------------------
# selftest


def selftest():
    failures = []

    def check(name, condition):
        if condition:
            print("selftest PASS %s" % name)
        else:
            failures.append(name)
            print("selftest FAIL %s" % name)

    rng = SplitMix64(0)
    check("splitmix64 known answer",
          rng.next_u64() == 0xE220A8397B1DCDAF)

    sample_a = build_sample(0x1234, smoke=True)
    sample_b = build_sample(0x1234, smoke=True)
    sample_c = build_sample(0x1235, smoke=True)
    check("sample determinism",
          canonical_json(sample_a) == canonical_json(sample_b))
    check("sample seed sensitivity",
          canonical_json(sample_a) != canonical_json(sample_c))
    ok = True
    seen = set()
    for spec, region in zip(SMOKE_REGIONS, sample_a):
        name, lo, hi, count, seed_count = spec
        ks = region["k_values"]
        ok = ok and region["name"] == name
        ok = ok and len(ks) == (count if count else hi - lo + 1)
        ok = ok and all(lo <= k <= hi for k in ks)
        ok = ok and ks == sorted(ks) and len(set(ks)) == len(ks)
        ok = ok and not (set(ks) & seen)
        seen |= set(ks)
        ok = ok and len(region["construction_seeds"]) == seed_count
        ok = ok and len(region["loss_seeds"]) == seed_count
        ok = ok and all(s != 0 for s in region["construction_seeds"])
        ok = ok and set(region["timing_k_values"]) <= set(ks)
    check("sample structure", ok)
    full = build_sample(0x9, smoke=False)
    check("full sample counts",
          [len(r["k_values"]) for r in full] == [99, 900, 400, 500, 1000,
                                                 4400] and
          [len(r["construction_seeds"]) for r in full] == [10, 5, 3, 3, 3, 3])

    methods = build_methods()
    check("tuning grid within cap", 1 <= len(methods) <= 100)
    check("tuning grid has wh1 arm",
          any(m["arm"] == "wh1" for m in methods))
    check("tuning grid has h15 family",
          any(m["family"] == "wh2-mixed-h15" for m in methods))

    chunks = chunk_k_values(list(range(2, 5000)))
    flattened = [k for chunk in chunks for k in chunk]
    check("chunking preserves K order/multiset",
          flattened == list(range(2, 5000)) and
          all(len(c) <= CHUNK_MAX_VALUES for c in chunks))

    check("ladder shape at K=64",
          overhead_ladder(64, 50, 8) ==
          [0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 18, 27, 32])
    ladder_64000 = overhead_ladder(64000, 50, 8)
    check("ladder reaches the cap exactly",
          ladder_64000[-1] == 32000 and
          ladder_64000[:10] == [0, 1, 2, 3, 4, 5, 6, 7, 8, 12] and
          all(b > a for a, b in zip(ladder_64000, ladder_64000[1:])))
    check("ladder small-K floor",
          overhead_ladder(2, 50, 8) == [0, 1, 2, 3, 4, 5, 6, 7, 8])
    check("pct threshold ceil", pct_threshold(30, 1) == 1 and
          pct_threshold(300, 1) == 3 and pct_threshold(2, 50) == 1)

    header = ("N,bb,heavy_family,mix_count,overhead,trials,success,"
              "rank_fail,error,fail_rate,inact_mu,inact_max,binary_def_mu,"
              "binary_def_max,heavy_gain_mu,heavy_gain_min,heavy_shortfall,"
              "solve_ms_mu,build_ms_mu,peel_ms_mu,project_ms_mu,"
              "residual_ms_mu,backsub_ms_mu,seed_attempt,block_xors_mu,"
              "block_muladds_mu,first_rank_fail,binary_def_hist,"
              "heavy_gain_hist,failure_trials,active_packet_peel_seed_xor")

    def wh2_row(k, overhead, success, seed_attempt=0, xors=100.0,
                muladds=50.0, error=0):
        return ("%d,2,periodic,3,%d,1,%d,%d,%d,0.0,1.0,1,0.0,0,0.0,0,0,"
                "0.1,0.1,0.1,0.1,0.1,0.1,%d,%.3f,%.3f,-1,0:1,0:1,,0x0" %
                (k, overhead, success, 0 if success else 1, error,
                 seed_attempt, xors, muladds))

    def wh2_cell_rows(k, success_extra):
        ladder = overhead_ladder(k, 50, 8)
        rows = []
        for level in ladder:
            if success_extra is not None and level == success_extra:
                rows.append(wh2_row(k, level, 1, xors=30.0, muladds=8.0))
                break
            rows.append(wh2_row(k, level, 0))
        return rows

    receipt = ("# untuned: construction_seed_explicit=1 "
               "construction_seed=0xabc construction_attempts=1 "
               "overhead_early_stop=1 encode_timing_reps=0 "
               "overhead_ladder_pct=50 overhead_ladder_min=8")
    text = "\n".join(
        ["# precodefail: trials=1", receipt, header] +
        wh2_cell_rows(10, 0) +
        wh2_cell_rows(20, 1) +
        wh2_cell_rows(100, 18) +
        wh2_cell_rows(30, None))
    cells = parse_wh2_cells(text, 0xABC, 50, 8)
    check("wh2 parser exact extras",
          cells[10]["extra"] == 0 and cells[10]["bracket_low"] == -1 and
          cells[20]["extra"] == 1 and cells[20]["bracket_low"] == 0 and
          cells[20]["work"] == (30.0, 8.0))
    check("wh2 parser geometric bracketing",
          cells[100]["extra"] == 18 and cells[100]["bracket_low"] == 12)
    check("wh2 parser censoring", cells[30]["extra"] is None)

    for bad_text, name in (
            (text.replace("construction_attempts=1",
                          "construction_attempts=2"),
             "wh2 parser rejects multi-attempt receipt"),
            (text.replace(receipt, "# other"),
             "wh2 parser requires untuned receipt"),
            (text.replace("overhead_ladder_pct=50",
                          "overhead_ladder_pct=0"),
             "wh2 parser rejects ladder mismatch"),
            ("\n".join([receipt, header, wh2_row(9, 0, 1, seed_attempt=3)]),
             "wh2 parser rejects seed escalation"),
            ("\n".join([receipt, header,
                        wh2_row(9, 0, 0), wh2_row(9, 1, 0)]),
             "wh2 parser rejects truncated censored cell"),
            ("\n".join([receipt, header,
                        wh2_row(9, 0, 0), wh2_row(9, 2, 0)]),
             "wh2 parser rejects off-ladder rows"),
            ("\n".join([receipt, header, wh2_row(9, 0, 0, error=1)]),
             "wh2 parser rejects solver errors")):
        try:
            parse_wh2_cells(bad_text, 0xABC, 50, 8)
            check(name, False)
        except ParseError:
            check(name, True)

    wh1_text = "\n".join([
        "# untunedcells: arm=untuned cseed=0x1f construction_attempts=1 "
        "pseed=31 dseed=0 seed=0x1 loss=0.1 pairbb=2 maxoh=4 cap_pct=50 "
        "cap_min=8 trial=0 threads=1 cells=4",
        "K,enc_ok,extra,censored,decode_limit",
        "10,1,0,0,0",
        "20,1,3,0,0",
        "30,0,-1,1,0",
        "40,1,-1,1,1",
    ])
    wh1_cells = parse_wh1_cells(wh1_text, 0x1F)
    check("wh1 parser",
          wh1_cells[10]["extra"] == 0 and
          wh1_cells[20]["extra"] == 3 and
          wh1_cells[20]["bracket_low"] == 2 and
          wh1_cells[30]["extra"] is None and
          wh1_cells[30]["enc_ok"] is False and
          wh1_cells[40]["decode_limit"] is True and
          wh1_cells[40]["extra"] is None)
    try:
        parse_wh1_cells(wh1_text, 0x20)
        check("wh1 parser rejects cseed mismatch", False)
    except ParseError:
        check("wh1 parser rejects cseed mismatch", True)
    try:
        parse_wh1_cells(wh1_text, 0x1F, tuned=True)
        check("wh1 parser rejects arm mismatch", False)
    except ParseError:
        check("wh1 parser rejects arm mismatch", True)
    wh1_tuned_text = wh1_text.replace(
        "arm=untuned cseed=0x1f", "arm=tuned cseed=0x0")
    check("wh1 parser tuned arm",
          parse_wh1_cells(wh1_tuned_text, 0, tuned=True)[10]["extra"] == 0)

    check("percentile with censoring",
          percentile_with_censoring([0.1, 0.2, 0.3, 0.4], 0, 0.50) == 0.2 and
          percentile_with_censoring([0.1, 0.2, 0.3], 1, 0.99) is None and
          percentile_with_censoring([], 4, 0.50) is None and
          percentile_with_censoring([0.1] * 99, 1, 0.99) == 0.1)

    outcomes = ([{"k": 100, "extra": 0}] * 985 +
                [{"k": 100, "extra": 1}] * 10 +
                [{"k": 100, "extra": 30}] * 3 +
                [{"k": 100, "extra": None, "enc_ok": False}] * 2)
    aggregate = aggregate_region(outcomes, PROTOCOL)
    check("aggregate float means over decoded only",
          abs(aggregate["overhead_fraction"]["mean_decoded"] -
              (10 * 0.01 + 3 * 0.30) / 998.0) < 1e-12)
    check("aggregate percentiles",
          aggregate["overhead_fraction"]["p50"] == 0.0 and
          aggregate["overhead_fraction"]["p99"] == 0.01)
    check("aggregate censoring fractions",
          abs(aggregate["censored_fraction"] - 0.002) < 1e-12 and
          abs(aggregate["construct_failure_fraction"] - 0.002) < 1e-12)
    limit_aggregate = aggregate_region(
        [{"k": 10, "extra": 0}] * 3 +
        [{"k": 10, "extra": None, "decode_limit": True}], PROTOCOL)
    check("aggregate decoder-limit fraction",
          abs(limit_aggregate["decoder_limit_fraction"] - 0.25) < 1e-12 and
          abs(limit_aggregate["censored_fraction"] - 0.25) < 1e-12)
    check("aggregate success-rate CDF readouts",
          abs(aggregate["success_rate_at_pct"]["1"] - 995 / 1000.0) < 1e-12
          and abs(aggregate["success_rate_at_pct"]["10"] - 995 / 1000.0)
          < 1e-12 and
          abs(aggregate["success_rate_at_pct"]["50"] - 998 / 1000.0) < 1e-12)
    check("aggregate keeps per-cell detail",
          len(aggregate["cells_detail"]) == 1000 and
          aggregate["cells_detail"][-1][1] == -1)
    all_censored = aggregate_region(
        [{"k": 10, "extra": None}] * 4, PROTOCOL)
    check("aggregate all-censored",
          all_censored["overhead_fraction"]["p50"] is None and
          all_censored["overhead_fraction"]["mean_decoded"] is None)

    wh2_timing_text = "\n".join([
        "# precodefail: ...",
        header,
        "# encode_timing,N=100,bb=2,heavy_family=periodic,mix_count=3,"
        "reps=3,solve_ok=1,solve_block_bytes=2,median_ns=1000000,"
        "seconds_per_k_symbols=0.001,symbols_per_sec=100000,sink=0",
        "100,2,periodic,3,0,1,1,0,0,0.0,1.0,1,0.0,0,0.0,0,0,0.1,0.1,0.1,"
        "0.1,0.1,0.1,0,1.0,1.0,-1,0:1,0:1,,0x0",
    ])
    wh2_timing = parse_wh2_timing(wh2_timing_text)
    check("wh2 timing parser",
          wh2_timing[100]["median_ns"] == 1000000 and
          wh2_timing[100]["solve_ok"])

    wh1_timing_text = "\n".join([
        "# enctime: reps=3 bb=1 cseed_explicit=0 cseed=0x0 counted=1 "
        "single_threaded=1 seed=0x517",
        "K,reps,ok,median_ns,seconds_per_k_symbols,symbols_per_sec,"
        "xor_bytes,muladd_bytes,sink",
        "100,3,1,20000,2e-05,5000000,1600,150,7",
    ])
    wh1_timing, counted = parse_wh1_timing(wh1_timing_text)
    check("wh1 timing parser",
          counted and wh1_timing[100]["xor_bytes"] == 1600)

    timing_aggregate = aggregate_timing(
        {100: {"median_ns": 1000000, "solve_ok": True},
         200: {"median_ns": 1000000, "solve_ok": True},
         300: {"median_ns": 3000000, "solve_ok": True}}, reps=3)
    # Rates are 1e5, 2e5, 1e5 symbols/sec; the median rate is 1e5.
    check("timing aggregation",
          timing_aggregate["median_symbols_per_sec"] == 100000.0 and
          abs(timing_aggregate["mean_seconds_per_k_symbols"] -
              (0.001 + 0.001 + 0.003) / 3) < 1e-12)

    check("canonical json stability",
          canonical_json({"b": 1, "a": [2, 3]}) ==
          '{"a":[2,3],"b":1}\n')

    with tempfile.TemporaryDirectory() as tmp:
        class SealArgs(object):
            out = tmp
            sample_seed = "0x77"
            smoke = True
            force = False
        check("seal writes manifest", cmd_seal(SealArgs()) == 0)
        manifest, sha_a = load_manifest(tmp)
        check("seal refuses overwrite", cmd_seal(SealArgs()) == 1)
        SealArgs.force = True
        check("seal --force reseal", cmd_seal(SealArgs()) == 0)
        _, sha_b = load_manifest(tmp)
        check("seal determinism end-to-end", sha_a == sha_b)
        total = sum(len(r["k_values"]) * len(r["construction_seeds"])
                    for r in manifest["regions"])
        check("smoke sample is 200 cells", total == 200)
        tasks = build_tasks(manifest, tmp, "/bin/bench", "/bin/whx", None,
                            ["wh1-legacy", "wh1-legacy-tuned",
                             "mixed-p244-baseline"],
                            3, False)
        cells_tasks = [t for t in tasks if t["kind"] == "cells"]
        timing_tasks = [t for t in tasks if t["kind"] == "timing"]
        check("task construction",
              len(timing_tasks) == 3 * 6 and
              len(cells_tasks) >= 3 * 6 * 2 and
              all(t["out"].startswith(tmp) for t in tasks))
        wh2_task = next(t for t in cells_tasks
                        if t["method"] == "mixed-p244-baseline")
        check("wh2 task command shape",
              "--construction-seed" in wh2_task["command"] and
              "--overhead-early-stop" in wh2_task["command"] and
              "--paired-overhead-stream" in wh2_task["command"] and
              "--overhead-ladder-pct" in wh2_task["command"] and
              "--overhead" not in wh2_task["command"])
        wh1_task = next(t for t in cells_tasks
                        if t["method"] == "wh1-legacy")
        check("wh1 task command shape",
              "untunedcells" in wh1_task["command"] and
              "--cseed" in wh1_task["command"] and
              "--cap-pct" in wh1_task["command"] and
              "--tuned" not in wh1_task["command"])
        wh1_tuned_task = next(t for t in cells_tasks
                              if t["method"] == "wh1-legacy-tuned")
        check("wh1 tuned task command shape",
              "--tuned" in wh1_tuned_task["command"] and
              "--cseed" not in wh1_tuned_task["command"])
        region_tasks = build_tasks(
            manifest, tmp, "/bin/bench", "/bin/whx", None,
            ["wh1-legacy"], 3, True, ["r6_20000_64000"])
        check("region filter",
              region_tasks and
              all(t["region"] == "r6_20000_64000" for t in region_tasks))

    if failures:
        print("selftest: %d FAILURES: %s" % (len(failures), failures))
        return 1
    print("selftest: all checks passed")
    return 0


# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    seal = sub.add_parser("seal", help="seal the shared stratified sample")
    seal.add_argument("--out", required=True)
    seal.add_argument("--sample-seed", required=True,
                      help="recorded RNG seed, e.g. 0x1234")
    seal.add_argument("--smoke", action="store_true",
                      help="tiny 200-cell sample with the same structure")
    seal.add_argument("--force", action="store_true")

    run = sub.add_parser("run", help="run pending tasks (resumable)")
    run.add_argument("--out", required=True)
    run.add_argument("--bench", help="wirehair_v2_bench binary")
    run.add_argument("--whx", help="whx binary (plain build)")
    run.add_argument("--whx-count", help="whx binary built with -DWH_COUNT")
    run.add_argument("--workers", type=int, default=120)
    run.add_argument("--methods", help="comma-separated method id filter")
    run.add_argument("--regions", help="comma-separated region name filter")
    run.add_argument("--timing-reps", type=int, default=3)
    run.add_argument("--max-tasks", type=int, default=None)
    run.add_argument("--skip-timing", action="store_true")
    run.add_argument("--force-owner", action="store_true")

    collect = sub.add_parser("collect", help="aggregate raw outputs to JSON")
    collect.add_argument("--out", required=True)
    collect.add_argument("--methods")
    collect.add_argument("--allow-partial", action="store_true")

    sub.add_parser("methods", help="print the tuning-point table")
    sub.add_parser("selftest", help="run harness self-tests")

    args = parser.parse_args()
    if args.command == "seal":
        return cmd_seal(args)
    if args.command == "run":
        return cmd_run(args)
    if args.command == "collect":
        return cmd_collect(args)
    if args.command == "methods":
        return cmd_methods(args)
    if args.command == "selftest":
        return selftest()
    return 2


if __name__ == "__main__":
    sys.exit(main())
