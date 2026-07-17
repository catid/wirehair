#!/usr/bin/env python3
"""Audit WH2 normalized-H15-v5 Stage G/T selection campaign outputs."""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import re
import shlex
import sys
import tempfile
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Callable
import wh2_h15_v5_selection_prepare as prepare_contract
ARMS = (('baseline', 0), ('c27', 39), ('c79', 121), ('c6f', 111), ('ca8', 168))
G_SCHEDULES = ('burst', 'adversarial', 'repair-only')
T_SCHEDULES = ('iid', 'burst', 'permutation', 'systematic-first', 'repair-only', 'adversarial')
T_LOSSES = ('0.35', '0.50', '0.65')
JOB_HEADER = ('job_id', 'stem', 'stage', 'schedule', 'seed_index', 'seed', 'loss', 'arm', 'salt', 'trials', 'k_count', 'k_csv')
REVEAL_HEADER = ('partition', 'index', 'domain', 'sha256', 'derived_u64', 'usage')
RAW_HEADER = ('N', 'bb', 'heavy_family', 'mix_count', 'overhead', 'trials', 'success', 'rank_fail', 'error', 'fail_rate', 'inact_mu', 'inact_max', 'binary_def_mu', 'binary_def_max', 'heavy_gain_mu', 'heavy_gain_min', 'heavy_shortfall', 'solve_ms_mu', 'build_ms_mu', 'peel_ms_mu', 'project_ms_mu', 'residual_ms_mu', 'backsub_ms_mu', 'seed_attempt', 'block_xors_mu', 'block_muladds_mu', 'first_rank_fail', 'binary_def_hist', 'heavy_gain_hist', 'failure_trials', 'active_packet_peel_seed_xor')
FAILURE_HEADER = ('witness_id', 'K', 'bb', 'certified_v4_salt', 'schedule', 'loss', 'seed_index', 'seed', 'trial', 'rank_fail', 'error', 'failure_trials', 'heavy_shortfall', 'inact_mu', 'inact_max', 'binary_def_mu', 'binary_def_max', 'seed_attempt', 'block_xors_mu', 'block_muladds_mu', 'source_base', 'source_sha256')
DISCOVERY_K_HEADER = ('K', 'certified_v4_salt', 'witness_count', 'witness_ids')
SCORE_HEADER = ('K', 'certified_v4_salt', 'salt', 'witness_count', 'rank_fail', 'error', 'heavy_shortfall', 'binary_def_mu_sum', 'binary_def_max', 'inact_mu_sum', 'inact_max', 'block_xors_mu_sum', 'block_muladds_mu_sum', 'exact_eligible', 'existing_v4_entry')
THERMAL_HEADER = ('utc', 'monotonic_s', 'cpu_busy_pct', 'cpu_avg_mhz', 'cpu_tctl_c', 'dimm_i2c1_50_c', 'dimm_i2c1_51_c', 'dimm_i2c1_52_c', 'dimm_i2c1_53_c', 'dimm_i2c2_50_c', 'dimm_i2c2_51_c', 'dimm_i2c2_52_c', 'dimm_i2c2_53_c', 'dimm_read_errors', 'load1', 'load5', 'load15', 'edac_ce', 'edac_ue')
PINNED = {'groups.tsv': '940db61d44fc03462f583c7e2e48bea5c93a830918c8f8e5bf3ba502e840cff4', 'discovery_k.csv': 'ecc722c5e44aee2ed468948cc4c7ee903ec4943b49ea4b67778b37a598e0b125', 'discovery_failures.csv': 'c9fe1730f632bc406c229391f137a30db5f933a5d21e405d822658b74e2d5dd7', 'exact_salt_scores.csv': 'ff238d4a21a340c866408ada53d22da474a880f054c1d50ecbba07007bf128e4', 'requested_k_salt.tsv': '93fe6d38b090d9054e0c293276d3caef7e0fd1dc2306105fdfd9351755b2c130'}
BASELINE_SHA256 = '4340fd9444b6dc2adbd94304b55e987d9887efd8f6d86f205511e7fda75cb41c'

class AuditError(RuntimeError):
    pass

def fail(message: str) -> None:
    raise AuditError(message)

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def require_hash(path: Path, expected: str) -> None:
    actual = sha256(path)
    if actual != expected:
        fail(f'{path}: SHA-256 mismatch: expected {expected}, found {actual}')

def uint(text: str, context: str, maximum: int | None=None) -> int:
    if not re.fullmatch('(?:0|[1-9][0-9]*)', text):
        fail(f'{context}: noncanonical unsigned integer {text!r}')
    value = int(text)
    if maximum is not None and value > maximum:
        fail(f'{context}: value exceeds {maximum}')
    return value

def salt_hex(text: str, context: str) -> int:
    if not re.fullmatch('0x(?:0|[1-9a-f][0-9a-f]*)', text):
        fail(f'{context}: noncanonical lowercase hexadecimal salt {text!r}')
    value = int(text, 16)
    if value > 255:
        fail(f'{context}: salt exceeds uint8')
    return value

def decimal(text: str, context: str) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail(f'{context}: invalid decimal {text!r}')
    if not value.is_finite() or value < 0:
        fail(f'{context}: decimal must be finite and nonnegative')
    return value

def parse_failure_ids(text: str, trials: int, context: str) -> tuple[int, ...]:
    if not text:
        return ()
    result = tuple((uint(part, context, trials - 1) for part in text.split('|')))
    if any((a >= b for a, b in zip(result, result[1:]))):
        fail(f'{context}: failure IDs are duplicate or unsorted')
    return result

def histogram(text: str, trials: int, context: str) -> dict[int, int]:
    result: dict[int, int] = {}
    for part in text.split('|'):
        fields = part.split(':')
        if len(fields) != 2:
            fail(f'{context}: malformed histogram')
        key = uint(fields[0], context)
        count = uint(fields[1], context)
        if key in result or count == 0 or (result and key <= max(result)):
            fail(f'{context}: histogram keys/counts are noncanonical')
        result[key] = count
    if sum(result.values()) != trials:
        fail(f'{context}: histogram counts do not sum to trials')
    return result

def read_v4(path: Path, expected_hash: str=PINNED['requested_k_salt.tsv']) -> dict[int, int]:
    require_hash(path, expected_hash)
    result: dict[int, int] = {}
    with path.open(newline='') as stream:
        for line, fields in enumerate(csv.reader(stream, delimiter='\t'), 1):
            if len(fields) != 2:
                fail(f'{path}:{line}: expected two TSV fields')
            k = uint(fields[0], f'{path}:{line}:K', 64000)
            value = salt_hex(fields[1], f'{path}:{line}:salt')
            if k in result:
                fail(f'{path}:{line}: duplicate K={k}')
            result[k] = value
    if set(result) != set(range(2, 64001)):
        fail(f'{path}: expected exact K range 2..64000')
    if sum((value == 0 for value in result.values())) != 63929:
        fail(f'{path}: expected exactly 63929 salt-zero entries')
    return result

def read_discovery(path: Path, v4: dict[int, int], expected_hash: str=PINNED['discovery_failures.csv']) -> tuple[int, ...]:
    require_hash(path, expected_hash)
    ks: set[int] = set()
    seen: set[tuple[int, str, int]] = set()
    seeds: dict[int, str] = {}
    count = 0
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != FAILURE_HEADER:
            fail(f'{path}: discovery-failure header mismatch')
        for line, row in enumerate(reader, 2):
            if None in row or any((value is None for value in row.values())):
                fail(f'{path}:{line}: malformed CSV row')
            k = uint(row['K'], f'{path}:{line}:K', 64000)
            schedule = row['schedule']
            seed_index = uint(row['seed_index'], f'{path}:{line}:seed_index', 2)
            if k < 2 or v4.get(k) != 0 or salt_hex(row['certified_v4_salt'], f'{path}:{line}:salt') != 0:
                fail(f'{path}:{line}: discovery witness is not a salt-zero K')
            if row['bb'] != '64' or schedule not in G_SCHEDULES or row['loss'] != '0.50':
                fail(f'{path}:{line}: discovery witness geometry mismatch')
            if row['trial'] != '0' or row['failure_trials'] != '0':
                fail(f'{path}:{line}: discovery trial identity mismatch')
            rank_fail = uint(row['rank_fail'], f'{path}:{line}:rank_fail', 1)
            errors = uint(row['error'], f'{path}:{line}:error', 1)
            if rank_fail != 1 or errors != 0:
                fail(f'{path}:{line}: discovery baseline must be an error-free rank failure')
            if not re.fullmatch('0x[0-9a-f]{16}', row['seed']):
                fail(f'{path}:{line}: noncanonical discovery seed')
            if seed_index in seeds and seeds[seed_index] != row['seed']:
                fail(f'{path}:{line}: inconsistent discovery seed index')
            seeds[seed_index] = row['seed']
            key = (k, schedule, seed_index)
            if key in seen:
                fail(f'{path}:{line}: duplicate discovery witness stratum')
            seen.add(key)
            ks.add(k)
            count += 1
    if count != 477 or len(ks) != 466 or set(seeds) != {0, 1, 2}:
        fail(f'{path}: expected 477 witnesses / 466 K / seed indices 0..2')
    return tuple(sorted(ks))

def read_selection_seeds(path: Path) -> list[dict[str, str]]:
    try:
        return prepare_contract.read_reveal(path)
    except prepare_contract.AuditError as exc:
        fail(str(exc))

def require_reveal_binding(manifest: dict[str, object]) -> None:
    if manifest.get('selection_reveal_sha256') != manifest.get('selection_seeds_sha256'):
        fail('selection reveal hash differs from canonical copied seed table')

def unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            fail(f'duplicate JSON key {key!r}')
        result[key] = value
    return result

def json_exact(actual: object, expected: object) -> bool:
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return (actual.keys() == expected.keys() and
                all(json_exact(actual[key], value) for key, value in expected.items()))
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(json_exact(a, b) for a, b in zip(actual, expected))
    return actual == expected

def read_jobs(experiment: Path, audit_frozen: bool=True) -> list[dict[str, object]]:
    meta = experiment / 'meta'
    manifest_path = meta / 'manifest.json'
    if not manifest_path.is_file() or manifest_path.is_symlink():
        fail(f'{manifest_path}: manifest must be a regular nonsymlink file')
    manifest_text = manifest_path.read_bytes().decode('utf-8')
    if '\r' in manifest_text or not manifest_text.endswith('\n'):
        fail(f'{manifest_path}: noncanonical line framing')
    manifest = json.loads(manifest_text, object_pairs_hook=unique_json_object)
    if not isinstance(manifest, dict):
        fail(f'{manifest_path}: manifest root must be an object')
    if json.dumps(manifest, indent=2, sort_keys=True) + '\n' != manifest_text:
        fail(f'{manifest_path}: manifest JSON is not canonical')
    require_reveal_binding(manifest)
    expected_keys = {'schema', 'source_commitment_sha256', 'source_commitment_binding', 'groups_sha256', 'discovery_k_sha256', 'selection_reveal_sha256', 'selection_seeds_sha256', 'jobs_sha256', 'immutable_meta_sha256', 'selection_seed_roles', 'arms', 'stage_g', 'stage_t', 'jobs', 'logical_trials'}
    expected_arms = [{'arm': arm, 'salt': salt, 'salt_hex': f'0x{salt:x}'} for arm, salt in ARMS]
    if set(manifest) != expected_keys or manifest.get('schema') != 'wirehair.wh2.normalized_h15_v5.selection.v1' or manifest.get('source_commitment_sha256') != 'f76eedbf4813597fa84af4b21a1236c957d11ff4a745585ea0442d9676645728' or (manifest.get('source_commitment_binding') != 'opaque whole-file procedural provenance only; no membership claim; no source commitment path was accepted or read') or (manifest.get('groups_sha256') != PINNED['groups.tsv']) or (manifest.get('discovery_k_sha256') != PINNED['discovery_k.csv']) or (not re.fullmatch('[0-9a-f]{64}', str(manifest.get('selection_reveal_sha256', '')))) or (not json_exact(manifest.get('jobs'), 1980)) or (not json_exact(manifest.get('logical_trials'), 1126695)) or (manifest.get('immutable_meta_sha256') != PINNED) or (manifest.get('selection_seed_roles') != {'0': 'StageG', '1': 'StageT', '2': 'StageT'}) or (not json_exact(manifest.get('arms'), expected_arms)) or (not json_exact(manifest.get('stage_g'), {'groups': 120, 'k': 63929, 'jobs': 1800, 'logical_trials': 958935})) or (not json_exact(manifest.get('stage_t'), {'seeds': 2, 'k': 466, 'jobs': 180, 'logical_trials': 167760})):
        fail(f'{manifest_path}: selection manifest contract mismatch')
    for name, digest in PINNED.items():
        if not (meta / name).is_file() or (meta / name).is_symlink():
            fail(f'{meta / name}: immutable metadata must not be a symlink')
        if audit_frozen or name in ('groups.tsv', 'discovery_k.csv'):
            require_hash(meta / name, digest)
    jobs_path = meta / 'jobs.tsv'
    seeds_path = meta / 'selection_seeds.tsv'
    if jobs_path.is_symlink() or seeds_path.is_symlink() or manifest_path.is_symlink():
        fail(f'{meta}: manifest/jobs/seeds must not be symlinks')
    if (manifest.get('jobs_sha256') != sha256(jobs_path) or
            manifest.get('selection_seeds_sha256') != sha256(seeds_path)):
        fail(f'{jobs_path}: hash does not match manifest')
    selection_seeds = read_selection_seeds(seeds_path)
    try:
        groups = prepare_contract.read_groups(meta / 'groups.tsv')
        discovery = prepare_contract.read_discovery(meta / 'discovery_k.csv')
        expected_jobs = prepare_contract.make_jobs(groups, discovery, selection_seeds)
    except prepare_contract.AuditError as exc:
        fail(str(exc))
    if audit_frozen:
        v4 = read_v4(meta / 'requested_k_salt.tsv')
        if read_discovery(meta / 'discovery_failures.csv', v4) != discovery:
            fail(f'{meta}: discovery_k.csv differs from discovery-failure support')
    with jobs_path.open(newline='') as stream:
        jobs_text = jobs_path.read_bytes()
        if b'\r' in jobs_text or not jobs_text.endswith(b'\n') or b'\n\n' in jobs_text:
            fail(f'{jobs_path}: noncanonical line framing')
        reader = csv.DictReader(stream, delimiter='\t')
        if tuple(reader.fieldnames or ()) != JOB_HEADER:
            fail(f'{jobs_path}: job header mismatch')
        rows = list(reader)
    if len(rows) != len(expected_jobs):
        fail(f'{jobs_path}: expected exactly 1980 jobs')
    parsed: list[dict[str, object]] = []
    for line, (row, expected) in enumerate(zip(rows, expected_jobs), 2):
        if None in row or any((value is None for value in row.values())) or any((row[key] != str(expected[key]) for key in JOB_HEADER)):
            fail(f'{jobs_path}:{line}: job differs from exact reconstructed grid')
        cooked: dict[str, object] = dict(row)
        cooked.update({'_ks': tuple(map(int, row['k_csv'].split(','))), '_trials': int(row['trials']), '_salt': int(row['salt'])})
        parsed.append(cooked)
    return parsed

def parse_preamble(line: str, path: Path) -> dict[str, str]:
    prefix = '# precodefail: '
    if not line.startswith(prefix):
        fail(f'{path}: missing precodefail preamble')
    result: dict[str, str] = {}
    for token in line[len(prefix):].split():
        if '=' not in token:
            fail(f'{path}: malformed preamble token {token!r}')
        key, value = token.split('=', 1)
        if not key or key in result:
            fail(f'{path}: duplicate or empty preamble key')
        result[key] = value
    return result

def validate_output(path: Path, job: dict[str, object]) -> list[dict[str, object]]:
    if not path.is_file() or path.is_symlink():
        fail(f'{path}: raw output must be a regular nonsymlink file')
    text = path.read_bytes().decode('utf-8')
    if '\r' in text or not text.endswith('\n') or '\n\n' in text:
        fail(f'{path}: noncanonical line endings or blank line')
    stream = io.StringIO(text, newline='')
    preamble = parse_preamble(stream.readline().rstrip('\n'), path)
    expected = {'trials': str(job['_trials']), 'threads': '1', 'completion': 'mixed', 'mixed_period': '32', 'mixed_gf256_rows': '11', 'mixed_gf16_rows': '4', 'mixed_geometry': 'shared-x', 'mixed_residue_skew': '0', 'mixed_residue_schedule': 'hashed', 'mixed_residue_hash_seed': '0x44', 'mixed_residue_hash_keyed': '1', 'mixed_independent_extension_residues': '1', 'mixed_extension_residue_seed_xor': '0x4e', 'source_hits_override': '0', 'packet_peel_seed_table': 'none', 'binary_dense_rows_override': '0', 'gf256_heavy_rows_override': '0', 'odd_packet_peel_seed_xor': '0x0', 'packet_row_seed_multiplier': '0x1', 'packet_row_seed_avalanche': '0', 'seed_block_bytes_override': '1280', 'overhead_stream': 'salted', 'full_payload_solve': '0', 'schedule': str(job['schedule'])}
    if set(preamble) != set(expected) | {'loss', 'seed', 'packet_peel_seed_xor'}:
        fail(f'{path}: unexpected or missing preamble keys')
    for key, value in expected.items():
        if preamble[key] != value:
            fail(f'{path}: {key}={preamble[key]!r}, expected {value!r}')
    expected_loss = format(float(str(job['loss'])), '.17g')
    if preamble['loss'] != expected_loss:
        fail(f'{path}: loss mismatch')
    expected_seed = f"0x{int(str(job['seed']), 16):x}"
    if preamble['seed'] != expected_seed:
        fail(f'{path}: seed mismatch')
    expected_salt = f"0x{int(job['_salt']):x}"
    if preamble['packet_peel_seed_xor'] != expected_salt:
        fail(f'{path}: packet peel salt mismatch')
    reader = csv.DictReader(stream)
    if tuple(reader.fieldnames or ()) != RAW_HEADER:
        fail(f'{path}: raw CSV header mismatch')
    rows: list[dict[str, object]] = []
    expected_ks = job['_ks']
    for line, row in enumerate(reader, 3):
        if None in row or any((value is None for value in row.values())):
            fail(f'{path}:{line}: malformed CSV row')
        k = uint(row['N'], f'{path}:{line}:N', 64000)
        trials = uint(row['trials'], f'{path}:{line}:trials')
        success = uint(row['success'], f'{path}:{line}:success', trials)
        rank_fail = uint(row['rank_fail'], f'{path}:{line}:rank_fail', trials)
        errors = uint(row['error'], f'{path}:{line}:error', trials)
        if (row['bb'], row['heavy_family'], row['mix_count'], row['overhead'], trials) != ('64', 'periodic', '2', '0', job['_trials']):
            fail(f'{path}:{line}: row geometry mismatch')
        if success + rank_fail + errors != trials:
            fail(f'{path}:{line}: outcome ledger mismatch')
        failures = parse_failure_ids(row['failure_trials'], trials, f'{path}:{line}:failure_trials')
        if len(failures) != rank_fail + errors:
            fail(f'{path}:{line}: failure ID/count mismatch')
        first = row['first_rank_fail']
        if rank_fail == 0:
            if first != '-1':
                fail(f'{path}:{line}: first_rank_fail must be -1')
        elif uint(first, f'{path}:{line}:first_rank_fail', trials - 1) not in failures:
            fail(f'{path}:{line}: first rank-failure ID is absent from failure_trials')
        if decimal(row['fail_rate'], f'{path}:{line}:fail_rate') != Decimal(rank_fail + errors) / trials:
            fail(f'{path}:{line}: fail_rate mismatch')
        for key in ('inact_mu', 'binary_def_mu', 'heavy_gain_mu', 'solve_ms_mu', 'build_ms_mu', 'peel_ms_mu', 'project_ms_mu', 'residual_ms_mu', 'backsub_ms_mu', 'block_xors_mu', 'block_muladds_mu'):
            decimal(row[key], f'{path}:{line}:{key}')
        for key in ('inact_max', 'binary_def_max', 'heavy_gain_min', 'seed_attempt'):
            uint(row[key], f'{path}:{line}:{key}')
        uint(row['heavy_shortfall'], f'{path}:{line}:heavy_shortfall', trials)
        binary_hist = histogram(row['binary_def_hist'], trials, f'{path}:{line}:binary histogram')
        heavy_hist = histogram(row['heavy_gain_hist'], trials, f'{path}:{line}:heavy histogram')
        binary_mean = sum((k * count for k, count in binary_hist.items())) / Decimal(trials)
        heavy_mean = sum((k * count for k, count in heavy_hist.items())) / Decimal(trials)
        if max(binary_hist) != int(row['binary_def_max']) or min(heavy_hist) != int(row['heavy_gain_min']) or binary_mean != Decimal(row['binary_def_mu']) or (heavy_mean != Decimal(row['heavy_gain_mu'])):
            fail(f'{path}:{line}: histogram summary mismatch')
        if row['active_packet_peel_seed_xor'] != expected_salt:
            fail(f'{path}:{line}: active salt mismatch')
        xors = decimal(row['block_xors_mu'], f'{path}:{line}:block_xors_mu') * 1000
        muls = decimal(row['block_muladds_mu'], f'{path}:{line}:block_muladds_mu') * 1000
        if xors != xors.to_integral() or muls != muls.to_integral() or int(xors) * trials % 1000 or int(muls) * trials % 1000:
            fail(f'{path}:{line}: work mean cannot represent integral trial totals')
        rows.append({'K': k, 'success': success, 'rank_fail': rank_fail, 'error': errors, 'failure_trials': failures, 'xors': int(xors), 'muls': int(muls)})
    if tuple((row['K'] for row in rows)) != expected_ks:
        fail(f'{path}: exact K row order/content mismatch')
    return rows

def validate_job(args: argparse.Namespace) -> None:
    experiment = args.experiment.resolve()
    jobs = read_jobs(experiment, audit_frozen=False)
    if args.job_id < 0 or args.job_id >= len(jobs):
        fail(f'job ID {args.job_id} is outside [0,{len(jobs)})')
    job = jobs[args.job_id]
    stem = str(job['stem'])
    stderr = experiment / 'raw' / f'{stem}.stderr'
    if not stderr.is_file() or stderr.is_symlink() or stderr.stat().st_size:
        fail(f'{stderr}: missing or nonempty stderr')
    validate_output(experiment / 'raw' / f'{stem}.csv', job)

def validate_temp_output(args: argparse.Namespace) -> None:
    jobs = read_jobs(args.experiment.resolve(), audit_frozen=False)
    if args.job_id < 0 or args.job_id >= len(jobs):
        fail(f'job ID {args.job_id} is outside [0,{len(jobs)})')
    if args.stdout.is_symlink() or args.stderr.is_symlink() or not args.stderr.is_file() or args.stderr.stat().st_size:
        fail(f'{args.stderr}: missing or nonempty stderr')
    validate_output(args.stdout, jobs[args.job_id])

def discovery_evidence(path: Path) -> tuple[dict[int, set[tuple[str, str, str, str]]], dict[int, int]]:
    strata: dict[int, set[tuple[str, str, str, str]]] = defaultdict(set)
    counts: dict[int, int] = defaultdict(int)
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != FAILURE_HEADER:
            fail(f'{path}: discovery-failure header mismatch')
        for row in reader:
            k = int(row['K'])
            strata[k].add(('D', row['seed'], row['schedule'], row['loss']))
            counts[k] += 1
    return (strata, counts)

def read_scores(path: Path, witnesses: dict[int, int]) -> dict[tuple[int, int], bool]:
    require_hash(path, PINNED['exact_salt_scores.csv'])
    masks: dict[int, int] = defaultdict(int)
    result: dict[tuple[int, int], bool] = {}
    count = 0
    candidate_salts = {value for _, value in ARMS[1:]}
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != SCORE_HEADER:
            fail(f'{path}: exact-score header mismatch')
        for line, row in enumerate(reader, 2):
            if None in row or any((value is None for value in row.values())):
                fail(f'{path}:{line}: malformed exact-score row')
            k = uint(row['K'], f'{path}:{line}:K', 64000)
            certified = salt_hex(row['certified_v4_salt'], f'{path}:{line}:certified')
            candidate = salt_hex(row['salt'], f'{path}:{line}:salt')
            rank_fail = uint(row['rank_fail'], f'{path}:{line}:rank_fail')
            errors = uint(row['error'], f'{path}:{line}:error')
            eligible = uint(row['exact_eligible'], f'{path}:{line}:eligible', 1)
            if k not in witnesses or certified != 0 or uint(row['witness_count'], f'{path}:{line}:witness_count') != witnesses[k] or (uint(row['existing_v4_entry'], f'{path}:{line}:existing', 1) != 0) or (eligible != int(rank_fail == 0 and errors == 0)) or masks[k] & 1 << candidate:
                fail(f'{path}:{line}: exact-score geometry/eligibility mismatch')
            for key in ('heavy_shortfall', 'binary_def_max', 'inact_max'):
                uint(row[key], f'{path}:{line}:{key}')
            for key in ('binary_def_mu_sum', 'inact_mu_sum', 'block_xors_mu_sum', 'block_muladds_mu_sum'):
                decimal(row[key], f'{path}:{line}:{key}')
            masks[k] |= 1 << candidate
            if candidate in candidate_salts:
                result[k, candidate] = bool(eligible and rank_fail == 0 and (errors == 0))
            count += 1
    full = (1 << 256) - 1
    if count != 119296 or set(masks) != set(witnesses) or any((mask != full for mask in masks.values())):
        fail(f'{path}: expected exact 466x256 score grid')
    return result

def validate_aux(experiment: Path, job: dict[str, object], binaries: dict[str, str], sealed_baseline: Path) -> None:
    stem = str(job['stem'])
    status = experiment / 'meta/status' / f'{stem}.ok'
    timing = experiment / 'raw' / f'{stem}.time'
    command = experiment / 'meta/commands' / f'{stem}.txt'
    if status.read_bytes() != b'ok\n':
        fail(f'{status}: expected exact ok status')
    time_text = timing.read_bytes().decode('utf-8')
    match = re.fullmatch('elapsed_s=([0-9]+(?:\\.[0-9]+)?) user_s=([0-9]+(?:\\.[0-9]+)?) sys_s=([0-9]+(?:\\.[0-9]+)?) cpu_pct=([0-9]+)% max_rss_kb=([0-9]+) exit=0\\n', time_text)
    if not match:
        fail(f'{timing}: malformed timing record')
    for value in match.groups()[:3]:
        decimal(value, f'{timing}:time')
    uint(match.group(4), f'{timing}:cpu_pct')
    uint(match.group(5), f'{timing}:rss')
    text = command.read_bytes().decode('utf-8')
    prefix = (f"job_id={job['job_id']}\tstage={job['stage']}\tarm={job['arm']}\t"
              f"seed_index={job['seed_index']}\tk_count={job['k_count']}\tcommand=")
    if not text.startswith(prefix) or not text.endswith(' \n') or '\r' in text or text.count('\n') != 1:
        fail(f'{command}: command-record envelope mismatch')
    argv = shlex.split(text[len(prefix):-1])
    if not argv:
        fail(f'{command}: empty command')
    if argv[0] != str(sealed_baseline):
        fail(f'{command}: command executable differs from sealed baseline')
    executable = argv[0]
    if executable not in binaries:
        binaries[executable] = sha256(Path(executable))
    if binaries[executable] != BASELINE_SHA256:
        fail(f'{command}: baseline binary hash mismatch')
    expected = ['precodefail', '--N', str(job['k_csv']), '--bb-list', '64', '--seed-block-bytes', '1280', '--overhead', '0', '--trials', str(job['trials']), '--threads', '1', '--loss', str(job['loss']), '--seed', str(job['seed']), '--schedule', str(job['schedule']), '--completion', 'mixed', '--mix-count', '2', '--packet-peel-seed-xor', str(job['salt']), '--mixed-gf256-rows', '11', '--mixed-gf16-rows', '4', '--mixed-period', '32', '--mixed-geometry', 'shared-x', '--mixed-residue-schedule', 'hashed', '--mixed-residue-hash-seed', '68', '--mixed-residue-hash-keyed', '--mixed-independent-extension-residues', '--mixed-extension-residue-seed-xor', '78']
    if argv[1:] != expected:
        fail(f'{command}: command arguments differ from sealed job geometry')

def artifact_sets(experiment: Path, jobs: list[dict[str, object]]) -> None:
    stems = {str(job['stem']) for job in jobs}
    expected_raw = {f'{stem}{suffix}' for stem in stems for suffix in ('.csv', '.stderr', '.time')}
    expected_status = {f'{stem}.ok' for stem in stems}
    expected_commands = {f'{stem}.txt' for stem in stems}
    for path, expected in ((experiment / 'raw', expected_raw), (experiment / 'meta/status', expected_status), (experiment / 'meta/commands', expected_commands)):
        if not path.is_dir() or path.is_symlink() or {entry.name for entry in path.iterdir()} != expected:
            fail(f'{path}: missing, extra, or non-file campaign artifacts')
        if any((not entry.is_file() or entry.is_symlink() for entry in path.iterdir())):
            fail(f'{path}: artifact directory contains a non-file or symlink')

def verify_artifact_manifest(experiment: Path, jobs: list[dict[str, object]]) -> None:
    expected_list: list[Path] = []
    for job in jobs:
        stem = str(job['stem'])
        expected_list.extend((experiment / 'raw' / f'{stem}.csv',
                              experiment / 'raw' / f'{stem}.stderr',
                              experiment / 'raw' / f'{stem}.time',
                              experiment / 'meta/status' / f'{stem}.ok',
                              experiment / 'meta/commands' / f'{stem}.txt'))
    expected_list.extend((experiment / 'meta/run_start.txt', experiment / 'meta/run_finish.txt',
                          experiment / 'thermal_interval.csv'))
    expected = set(expected_list)
    manifest = experiment / 'meta/artifact_manifest.sha256'
    if not manifest.is_file() or manifest.is_symlink():
        fail(f'{manifest}: artifact manifest must be a regular nonsymlink file')
    found: dict[Path, str] = {}
    manifest_text = manifest.read_bytes().decode('utf-8')
    if '\r' in manifest_text or not manifest_text.endswith('\n') or '\n\n' in manifest_text:
        fail(f'{manifest}: noncanonical line framing')
    for line, text in enumerate(manifest_text.splitlines(), 1):
        match = re.fullmatch('([0-9a-f]{64})  (/.+)', text)
        if not match:
            fail(f'{manifest}:{line}: malformed SHA-256 record')
        path = Path(match.group(2))
        if (str(path) != match.group(2) or path in found or path not in expected or
                not path.is_file() or path.is_symlink()):
            fail(f'{manifest}:{line}: duplicate, unexpected, or symlink path')
        found[path] = match.group(1)
    if list(found) != expected_list or len(found) != 9903:
        fail(f'{manifest}: expected exact ordered 9903 artifact paths')
    for path, digest in found.items():
        require_hash(path, digest)

def expected_selection_seal_paths(experiment: Path) -> set[Path]:
    scripts = Path(__file__).resolve().parent
    return ({experiment / 'meta' / name for name in
             ('manifest.json', 'selection_seeds.tsv', 'jobs.tsv', 'groups.tsv',
              'discovery_k.csv', 'discovery_failures.csv', 'exact_salt_scores.csv',
              'requested_k_salt.tsv')} |
            {scripts / name for name in ('wh2_h15_v5_selection_analyze.py',
                                         'wh2_h15_v5_selection_prepare.py',
                                         'wh2_h15_v5_selection_job.sh',
                                         'wh2_h15_v5_selection_launch.sh')})

def selection_baseline_from_records(experiment: Path, records: dict[Path, str]) -> Path:
    required = expected_selection_seal_paths(experiment)
    remaining = set(records) - required
    if len(records) != 13 or not required.issubset(records) or len(remaining) != 1:
        fail('selection seal path roles mismatch')
    baseline = next(iter(remaining))
    if records[baseline] != BASELINE_SHA256:
        fail('selection seal baseline digest mismatch')
    return baseline.resolve()

def verify_selection_seal(experiment: Path) -> Path:
    seal = experiment / 'meta/selection_seal.sha256'
    if not seal.is_file() or seal.is_symlink():
        fail(f'{seal}: selection seal must be a regular nonsymlink file')
    text = seal.read_bytes().decode('utf-8')
    if '\r' in text or not text.endswith('\n') or '\n\n' in text:
        fail(f'{seal}: noncanonical line framing')
    records: dict[Path, str] = {}
    for line, value in enumerate(text.splitlines(), 1):
        match = re.fullmatch(r'([0-9a-f]{64})  (/.+)', value)
        if not match:
            fail(f'{seal}:{line}: malformed SHA-256 record')
        path = Path(match.group(2))
        if str(path) != match.group(2) or path in records or not path.is_file() or path.is_symlink():
            fail(f'{seal}:{line}: duplicate, missing, or symlink input')
        records[path] = match.group(1)
    baseline = selection_baseline_from_records(experiment, records)
    for path, digest in records.items():
        require_hash(path, digest)
    return baseline

def read_kv(path: Path, keys: set[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    data = path.read_bytes().decode('utf-8')
    if '\r' in data or not data.endswith('\n') or '\n\n' in data:
        fail(f'{path}: noncanonical line framing')
    for line, text in enumerate(data.splitlines(), 1):
        if '=' not in text:
            fail(f'{path}:{line}: malformed key/value')
        key, value = text.split('=', 1)
        if key in result or key not in keys or (not value):
            fail(f'{path}:{line}: duplicate/unexpected/empty field')
        result[key] = value
    if set(result) != keys:
        fail(f'{path}: missing run fields')
    return result

def environment_audit(experiment: Path) -> dict[str, object]:
    start_keys = {'started_utc', 'workers', 'threads_per_invocation', 'launcher_nice', 'job_count', 'thermal_csv', 'thermal_line_start', 'edac_ce_start', 'edac_ue_start'}
    finish_keys = {'finished_utc', 'xargs_rc', 'status_count', 'stdout_count', 'stderr_count', 'stderr_nonempty', 'time_count', 'command_count', 'raw_entry_count', 'status_entry_count', 'command_entry_count', 'artifact_names_ok', 'thermal_line_end', 'edac_ce_end', 'edac_ue_end'}
    start = read_kv(experiment / 'meta/run_start.txt', start_keys)
    finish = read_kv(experiment / 'meta/run_finish.txt', finish_keys)
    utc_pattern = r'[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(?:\.[0-9]+)?Z'
    if ([start[key] for key in ('workers', 'threads_per_invocation', 'launcher_nice', 'job_count')]
        != ['120', '1', '0', '1980'] or start['thermal_csv'] != '/tmp/wirehair-enoq-thermal.csv' or
        [finish[key] for key in ('xargs_rc', 'status_count', 'stdout_count', 'stderr_count',
          'stderr_nonempty', 'time_count', 'command_count', 'raw_entry_count',
          'status_entry_count', 'command_entry_count', 'artifact_names_ok')]
        != ['0', '1980', '1980', '1980', '0', '1980', '1980', '5940', '1980', '1980', '1'] or
        start['edac_ce_start'] != finish['edac_ce_end'] or
        start['edac_ue_start'] != finish['edac_ue_end']):
        fail('run ledger/count/EDAC audit failed')
    if not re.fullmatch(utc_pattern, start['started_utc']) or not re.fullmatch(utc_pattern, finish['finished_utc']):
        fail('run UTC timestamp is noncanonical')
    for value in (start['thermal_line_start'], finish['thermal_line_end'], start['edac_ce_start'], start['edac_ue_start']):
        uint(value, 'run ledger')
    if int(finish['thermal_line_end']) < int(start['thermal_line_start']):
        fail('thermal interval shrank')
    cpu_max = Decimal(0)
    cpu_busy_min: Decimal | None = None
    cpu_busy_sum = Decimal(0)
    dimm_max = Decimal(0)
    dimm_valid = 0
    dimm_gap_rows = 0
    dimm_missing = 0
    gaps = 0
    max_gap = Decimal(0)
    last = None
    rows = 0
    read_errors = 0
    thermal_path = experiment / 'thermal_interval.csv'
    if not thermal_path.is_file() or thermal_path.is_symlink():
        fail(f'{thermal_path}: thermal interval must be a regular nonsymlink file')
    thermal_text = thermal_path.read_bytes()
    if b'\r' in thermal_text or not thermal_text.endswith(b'\n') or b'\n\n' in thermal_text:
        fail(f'{thermal_path}: noncanonical line framing')
    with thermal_path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != THERMAL_HEADER:
            fail('thermal CSV header mismatch')
        for line, row in enumerate(reader, 2):
            if None in row or any((value is None for value in row.values())):
                fail(f'thermal CSV:{line}: malformed row')
            if not re.fullmatch(utc_pattern, row['utc']):
                fail(f'thermal:{line}: noncanonical UTC timestamp')
            now = decimal(row['monotonic_s'], f'thermal:{line}')
            cpu = decimal(row['cpu_tctl_c'], f'thermal:{line}')
            dimms = [decimal(row[key], f'thermal:{line}') for key in THERMAL_HEADER[5:13] if row[key]]
            missing = 8 - len(dimms)
            errors = uint(row['dimm_read_errors'], f'thermal:{line}', 8)
            busy = decimal(row['cpu_busy_pct'], f'thermal:{line}')
            for key in ('cpu_avg_mhz', 'load1', 'load5', 'load15'):
                decimal(row[key], f'thermal:{line}')
            if missing != errors or busy > 100:
                fail(f'thermal:{line}: DIMM gap/error or CPU busy mismatch')
            if uint(row['edac_ce'], f'thermal:{line}') != int(start['edac_ce_start']) or uint(row['edac_ue'], f'thermal:{line}') != int(start['edac_ue_start']):
                fail(f'thermal:{line}: EDAC changed within run')
            read_errors += errors
            cpu_busy_min = busy if cpu_busy_min is None else min(cpu_busy_min, busy)
            cpu_busy_sum += busy
            dimm_missing += missing
            dimm_gap_rows += int(missing > 0)
            if last is not None:
                gap = now - last
                if gap <= 0:
                    fail(f'thermal:{line}: nonincreasing monotonic time')
                if gap > Decimal('1.5'):
                    gaps += 1
                max_gap = max(max_gap, gap)
            last = now
            cpu_max = max(cpu_max, cpu)
            if dimms:
                dimm_max = max(dimm_max, *dimms)
                dimm_valid += len(dimms)
            rows += 1
    if rows < 2 or rows != int(finish['thermal_line_end']) - int(start['thermal_line_start']) or dimm_valid == 0:
        fail('thermal row interval mismatch')
    if max_gap > Decimal('10'):
        fail('thermal sampling has an internal monotonic gap above 10 seconds')
    if cpu_busy_min is None or cpu_busy_min < Decimal('90'):
        fail('thermal sampling CPU busy minimum is below 90 percent')
    cpu_busy_mean = cpu_busy_sum / rows
    if cpu_busy_mean < Decimal('98'):
        fail('thermal sampling CPU busy mean is below 98 percent')
    return {'samples': rows, 'cpu_tctl_max_c': str(cpu_max), 'dimm_max_c': str(dimm_max),
            'cpu_busy_min_pct': str(cpu_busy_min),
            'cpu_busy_mean_pct': str(cpu_busy_mean),
            'dimm_valid_reads': dimm_valid, 'dimm_missing_reads': dimm_missing,
            'dimm_gap_rows': dimm_gap_rows, 'dimm_read_errors': read_errors,
            'sampling_gap_count': gaps, 'max_monotonic_gap_s': str(max_gap),
            'edac_ce_delta': 0, 'edac_ue_delta': 0}

def fresh_evidence(jobs: list[dict[str, object]], experiment: Path, strata: dict[int, set[tuple[str, str, str, str]]], sealed_baseline: Path) -> dict[int, dict[str, object]]:
    evidence: dict[int, dict[str, object]] = {}
    binaries: dict[str, str] = {}
    for offset in range(0, len(jobs), 5):
        block = jobs[offset:offset + 5]
        outputs = []
        for job in block:
            stem = str(job['stem'])
            validate_aux(experiment, job, binaries, sealed_baseline)
            stderr = experiment / 'raw' / f'{stem}.stderr'
            if stderr.read_bytes():
                fail(f'{stderr}: stderr must be empty')
            outputs.append(validate_output(experiment / 'raw' / f'{stem}.csv', job))
        for row_index, baseline in enumerate(outputs[0]):
            k = int(baseline['K'])
            if baseline['error']:
                fail(f'K={k}: baseline selection arm reported an error')
            state = evidence.setdefault(k, {'fresh': 0, 'strata': set(strata.get(k, set())), 'candidate': [[0, 0, 0, 0] for _ in ARMS[1:]]})
            base_failures = set(baseline['failure_trials'])
            if baseline['rank_fail']:
                state['fresh'] = int(state['fresh']) + int(baseline['rank_fail'])
                state['strata'].add((str(block[0]['stage']), str(block[0]['seed']), str(block[0]['schedule']), str(block[0]['loss'])))
            for candidate_index, candidate_rows in enumerate(outputs[1:]):
                candidate = candidate_rows[row_index]
                candidate_failures = set(candidate['failure_trials'])
                values = state['candidate'][candidate_index]
                values[0] += len(base_failures & candidate_failures)
                values[1] += len(candidate_failures - base_failures)
                values[2] += int(candidate['error'])
                values[3] += len(base_failures - candidate_failures)
    return evidence

def candidate_gates(k: int, candidate_index: int, v4: dict[int, int], evidence: dict[int, dict[str, object]], scores: dict[tuple[int, int], bool], discovery_k: set[int]) -> tuple[bool, list[object]]:
    state = evidence.get(k, {'fresh': 0, 'strata': set(), 'candidate': [[0, 0, 0, 0] for _ in ARMS[1:]]})
    values = state['candidate'][candidate_index]
    salt = ARMS[candidate_index + 1][1]
    gates = [v4[k] == 0, int(state['fresh']) > 0, len(state['strata']) >= 2, k not in discovery_k or scores.get((k, salt), False), values[0] == 0, values[2] == 0, values[1] == 0]
    return (all(gates), [*gates, int(state['fresh']), len(state['strata']), *values])

def choose_salts(v4: dict[int, int], evidence: dict[int, dict[str, object]], scores: dict[tuple[int, int], bool], discovery_k: set[int]) -> tuple[dict[int, int], list[list[object]]]:
    selected: dict[int, int] = {}
    ledger: list[list[object]] = []
    for k in range(2, 64001):
        selected[k] = v4[k]
        if v4[k] != 0:
            continue
        for candidate_index, (arm, salt) in enumerate(ARMS[1:]):
            eligible, detail = candidate_gates(k, candidate_index, v4, evidence, scores, discovery_k)
            ledger.append([k, f'0x{salt:x}', arm, candidate_index + 1, *detail, int(eligible), 0])
            if eligible and selected[k] == 0:
                selected[k] = salt
        for row in ledger[-4:]:
            row[-1] = int(row[1] == f'0x{selected[k]:x}' and selected[k] != 0)
    return (selected, ledger)

def add_work(total: list[int], baseline: dict[str, object], selected: dict[str, object], trials: int) -> None:
    if baseline['success'] == trials and selected['success'] == trials:
        total[0] += 1
        total[1] += trials
        total[2] += int(baseline['xors']) * trials
        total[3] += int(selected['xors']) * trials
        total[4] += int(baseline['muls']) * trials
        total[5] += int(selected['muls']) * trials
    else:
        total[6] += 1
        total[7] += trials
        total[8] += int(baseline['success'] != trials)
        total[9] += int(selected['success'] != trials)

def add_recovery(total: list[int], baseline: dict[str, object], selected: dict[str, object], trials: int) -> None:
    base = set(baseline['failure_trials'])
    choice = set(selected['failure_trials'])
    values = [trials, int(baseline['success']), int(baseline['rank_fail']), int(baseline['error']), int(selected['success']), int(selected['rank_fail']), int(selected['error']), len(base - choice), len(choice - base), len(base & choice)]
    for index, value in enumerate(values):
        total[index] += value

def work_metrics(jobs: list[dict[str, object]], experiment: Path, selected_salts: dict[int, int]) -> tuple[dict[str, list[int]], dict[int, list[int]], dict[str, list[int]]]:
    totals: dict[str, list[int]] = defaultdict(lambda: [0] * 10)
    per_k: dict[int, list[int]] = defaultdict(lambda: [0] * 10)
    recovery: dict[str, list[int]] = defaultdict(lambda: [0] * 10)
    arm_by_salt = {salt: arm for arm, salt in ARMS}
    totals['selected-K-only']
    recovery['selected-K-only']
    for arm, _ in ARMS[1:]:
        totals[f'arm:{arm}']
        recovery[f'selected-arm:{arm}']
    for offset in range(0, len(jobs), 5):
        block = jobs[offset:offset + 5]
        outputs = [validate_output(experiment / 'raw' / f"{job['stem']}.csv", job)
                   for job in block]
        for row_index, baseline in enumerate(outputs[0]):
            k = int(baseline['K'])
            salt = selected_salts[k]
            arm_index = next((i for i, (_, value) in enumerate(ARMS) if value == salt)) if salt else 0
            chosen = outputs[arm_index][row_index]
            trials = int(block[0]['_trials'])
            key = (f"{block[0]['stage']}|{block[0]['seed']}|"
                   f"{block[0]['schedule']}|{block[0]['loss']}")
            add_work(totals['final-overlay'], baseline, chosen, trials)
            add_work(totals[f'stratum:{key}'], baseline, chosen, trials)
            add_recovery(recovery['final-overlay'], baseline, chosen, trials)
            add_recovery(recovery[f'final-stratum:{key}'], baseline, chosen, trials)
            add_recovery(recovery[f'final-K:{k}'], baseline, chosen, trials)
            for tested_index, (tested_arm, _) in enumerate(ARMS[1:], 1):
                tested = outputs[tested_index][row_index]
                add_recovery(recovery[f'tested-arm:{tested_arm}'], baseline, tested, trials)
                add_recovery(recovery[f'tested-stratum:{key}|arm:{tested_arm}'], baseline, tested, trials)
                add_recovery(recovery[f'tested-K:{k}|arm:{tested_arm}'], baseline, tested, trials)
            if salt != 0:
                add_work(totals['selected-K-only'], baseline, chosen, trials)
                add_work(totals[f'arm:{arm_by_salt[salt]}'], baseline, chosen, trials)
                add_work(per_k[k], baseline, chosen, trials)
                add_recovery(recovery['selected-K-only'], baseline, chosen, trials)
                add_recovery(recovery[f'selected-arm:{arm_by_salt[salt]}'], baseline, chosen, trials)
    return (totals, per_k, recovery)

def percent(selected: int, baseline: int) -> str:
    return '' if baseline == 0 else f'{(Decimal(selected) / baseline - 1) * 100:.6f}'

def metric_row(scope: str, values: list[int]) -> list[object]:
    return [scope, *values[:2], values[6], values[7], values[8], values[9], values[2], values[3], percent(values[3], values[2]), values[4], values[5], percent(values[5], values[4])]

def caps_pass(values: list[int]) -> bool:
    return values[2] > 0 and values[4] > 0 and (values[3] * 100 <= values[2] * 101) and (values[5] * 10000 <= values[4] * 10025)

def paired_p(repairs: int, introductions: int) -> str:
    total = repairs + introductions
    if total == 0:
        return '1'
    numerator = min(2 ** total, 2 * sum((math.comb(total, index) for index in range(min(repairs, introductions) + 1))))
    return f'{Decimal(numerator) / Decimal(2 ** total):.12g}'

def recovery_row(scope: str, values: list[int]) -> list[object]:
    base_fail = values[2] + values[3]
    discordant = values[7] + values[8]
    return [scope, *values, discordant, '' if base_fail == 0 else f'{Decimal(values[7]) / base_fail:.12g}', '' if values[1] == 0 else f'{Decimal(values[8]) / values[1]:.12g}', paired_p(values[7], values[8])]

def write_analysis(experiment: Path, v4: dict[int, int], selected: dict[int, int], evidence: dict[int, dict[str, object]], ledger: list[list[object]], work: dict[str, list[int]], per_k: dict[int, list[int]], recovery: dict[str, list[int]], environment: dict[str, object]) -> None:
    out = experiment / 'analysis'
    if out.is_symlink() or out.exists() and (not out.is_dir() or any(out.iterdir())):
        fail(f'{out}: analysis output must be absent or empty')
    out.mkdir(exist_ok=True)
    ledger_header = ('K', 'candidate_salt', 'arm', 'priority', 'source_zero', 'fresh_failure', 'recurrence', 'discovery_repair', 'selection_repair', 'zero_errors', 'no_introduction', 'fresh_rank_fail', 'combined_strata', 'missed_repairs', 'introduced_failures', 'candidate_errors', 'repaired_trials', 'eligible', 'selected')
    with (out / 'decision_ledger.csv').open('x', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(ledger_header)
        writer.writerows(ledger)
    with (out / 'selection_table.tsv').open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(('K', 'v4_salt', 'v5_salt', 'changed', 'decision'))
        for k in range(2, 64001):
            changed = selected[k] != v4[k]
            decision = 'immutable-nonzero-v4' if v4[k] else 'selected' if changed else 'retain-zero'
            writer.writerow((k, f'0x{v4[k]:x}', f'0x{selected[k]:x}', int(changed), decision))
    rows = [metric_row(scope, values) for scope, values in sorted(work.items())]
    usable = {k: values for k, values in per_k.items() if values[2] and values[4]}
    worst_x = max(usable, key=lambda k: Decimal(usable[k][3]) / usable[k][2], default=None)
    worst_m = max(usable, key=lambda k: Decimal(usable[k][5]) / usable[k][4], default=None)
    for label, k in (('worst-K-xor', worst_x), ('worst-K-muladds', worst_m)):
        if k is not None:
            rows.append(metric_row(f'{label}:{k}', usable[k]))
    with (out / 'work_metrics.csv').open('x', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(('scope', 'included_rows', 'trial_weight', 'excluded_rows', 'excluded_trial_weight', 'baseline_nonfull', 'selected_nonfull', 'baseline_xor_milli_trials', 'selected_xor_milli_trials', 'xor_delta_pct', 'baseline_muladd_milli_trials', 'selected_muladd_milli_trials', 'muladd_delta_pct'))
        writer.writerows(rows)
    recovery_rows = [recovery_row(scope, values) for scope, values in sorted(recovery.items())]
    worst = max((item for item in recovery.items() if item[0].startswith('final-K:')),
                key=lambda item: (Decimal(item[1][5] + item[1][6]) / item[1][0],
                                  item[1][8], -int(item[0].split(':')[1])), default=None)
    if worst:
        recovery_rows.append(recovery_row('worst-K:' + worst[0].split(':', 1)[1], worst[1]))
    recovery_header = ('scope', 'trials', 'baseline_success', 'baseline_rank_fail',
                       'baseline_error', 'selected_success', 'selected_rank_fail',
                       'selected_error', 'repairs', 'introductions', 'overlaps',
                       'discordants', 'recovery_rate', 'introduction_rate',
                       'paired_binomial_two_sided_p')
    with (out / 'recovery_metrics.csv').open('x', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(recovery_header)
        writer.writerows(recovery_rows)
    final = work['final-overlay']
    caps = caps_pass(final)
    recovery_summary_names = [name for name in sorted(recovery)
                              if name in ('final-overlay', 'selected-K-only') or
                              name.startswith(('tested-arm:', 'selected-arm:'))]
    summary = {'schema': 'wirehair.wh2.normalized_h15_v5.selection_analysis.v1',
               'decision': 'GO' if caps else 'NO-GO', 'xor_cap_pct': '1.00',
               'muladd_cap_pct': '0.25', 'caps_pass': caps,
               'jobs': 1980, 'campaign_logical_trials': 1126695,
               'candidate_priority': ['0x27', '0x79', '0x6f', '0xa8'],
               'changed_k': sum((selected[k] != v4[k] for k in selected)),
               'selected_by_arm': {arm: sum(selected[k] == salt and v4[k] == 0
                                               for k in selected)
                                   for arm, salt in ARMS[1:]},
               'immutable_nonzero_v4_k': sum((v4[k] != 0 for k in v4)),
               'fresh_failure_k': sum((int(state['fresh']) > 0 for state in evidence.values())),
               'environment': environment,
               'recovery': {name: dict(zip(recovery_header[1:], recovery_row(name, recovery[name])[1:]))
                            for name in recovery_summary_names},
               'work': {row[0]: {'included_rows': row[1], 'trial_weight': row[2],
                                  'excluded_rows': row[3], 'xor_delta_pct': row[9],
                                  'muladd_delta_pct': row[12]} for row in rows}}
    with (out / 'summary.json').open('x') as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write('\n')

def analyze(args: argparse.Namespace) -> None:
    experiment = args.experiment.resolve()
    sealed_baseline = verify_selection_seal(experiment)
    jobs = read_jobs(experiment)
    artifact_sets(experiment, jobs)
    verify_artifact_manifest(experiment, jobs)
    environment = environment_audit(experiment)
    v4 = read_v4(experiment / 'meta/requested_k_salt.tsv')
    require_hash(experiment / 'meta/discovery_failures.csv', PINNED['discovery_failures.csv'])
    strata, witnesses = discovery_evidence(experiment / 'meta/discovery_failures.csv')
    scores = read_scores(experiment / 'meta/exact_salt_scores.csv', witnesses)
    evidence = fresh_evidence(jobs, experiment, strata, sealed_baseline)
    selected, ledger = choose_salts(v4, evidence, scores, set(witnesses))
    work, per_k, recovery = work_metrics(jobs, experiment, selected)
    write_analysis(experiment, v4, selected, evidence, ledger, work, per_k, recovery, environment)

def expect_failure(action: Callable[[], object], label: str) -> None:
    try:
        action()
    except AuditError:
        return
    fail(f'selftest mutation unexpectedly passed: {label}')

def selftest(_args: argparse.Namespace) -> None:
    job: dict[str, object] = {'stem': 'fixture', 'schedule': 'burst', 'seed': '0x0123456789abcdef', 'loss': '0.50', '_salt': 39, '_trials': 2, '_ks': (7,)}
    preamble = '# precodefail: trials=2 threads=1 loss=0.5 seed=0x123456789abcdef completion=mixed mixed_period=32 mixed_gf256_rows=11 mixed_gf16_rows=4 mixed_geometry=shared-x mixed_residue_skew=0 mixed_residue_schedule=hashed mixed_residue_hash_seed=0x44 mixed_residue_hash_keyed=1 mixed_independent_extension_residues=1 mixed_extension_residue_seed_xor=0x4e source_hits_override=0 packet_peel_seed_xor=0x27 packet_peel_seed_table=none binary_dense_rows_override=0 gf256_heavy_rows_override=0 odd_packet_peel_seed_xor=0x0 packet_row_seed_multiplier=0x1 packet_row_seed_avalanche=0 seed_block_bytes_override=1280 overhead_stream=salted full_payload_solve=0 schedule=burst\n'
    values = {key: '0' for key in RAW_HEADER}
    values.update({'N': '7', 'bb': '64', 'heavy_family': 'periodic', 'mix_count': '2', 'trials': '2', 'success': '1', 'rank_fail': '1', 'fail_rate': '0.50000000', 'inact_mu': '2.000', 'inact_max': '2', 'binary_def_mu': '1.000', 'binary_def_max': '1', 'heavy_gain_mu': '1.000', 'heavy_gain_min': '1', 'solve_ms_mu': '0.001', 'build_ms_mu': '0.001', 'peel_ms_mu': '0.001', 'project_ms_mu': '0.001', 'residual_ms_mu': '0.001', 'backsub_ms_mu': '0.001', 'block_xors_mu': '10.000', 'block_muladds_mu': '20.000', 'first_rank_fail': '1', 'binary_def_hist': '1:2', 'heavy_gain_hist': '1:2', 'failure_trials': '1', 'active_packet_peel_seed_xor': '0x27'})
    with tempfile.TemporaryDirectory(prefix='wh2-selection-analyze-') as temp_text:
        temp = Path(temp_text)
        valid = temp / 'valid.csv'
        with valid.open('w', newline='') as stream:
            stream.write(preamble)
            writer = csv.DictWriter(stream, fieldnames=RAW_HEADER, lineterminator='\n')
            writer.writeheader()
            writer.writerow(values)
        validate_output(valid, job)
        linked = temp / 'linked.csv'
        linked.symlink_to(valid)
        expect_failure(lambda: validate_output(linked, job), 'symlink raw output')
        mutations = (('seed', 'seed=0x123456789abcdef', 'seed=0x1'), ('loss', 'loss=0.5', 'loss=0.35'), ('salt', 'packet_peel_seed_xor=0x27', 'packet_peel_seed_xor=0x79'), ('K', '\n7,64,', '\n8,64,'), ('failure ID', ',1,0x27\n', ',0|1,0x27\n'), ('ledger', ',1,1,0,0.50000000,', ',1,1,1,0.50000000,'))
        for label, old, new in mutations:
            path = temp / f"bad-{label.replace(' ', '-')}.csv"
            path.write_text(valid.read_text().replace(old, new, 1))
            expect_failure(lambda path=path: validate_output(path, job), label)
        seed_rows = []
        for index, usage in enumerate(('future-StageT-looking', 'misleading-StageG-looking', 'opaque-purpose')):
            digest = hashlib.sha256(f'selection-seed-{index}'.encode()).hexdigest()
            seed_rows.append({'partition': 'selection', 'index': index,
                              'domain': f'wirehair.selection.selftest.{index}', 'sha256': digest,
                              'derived_u64': '0x' + digest[:16], 'usage': usage})
        seeds_path = temp / 'selection_seeds.tsv'
        prepare_contract.write_tsv(seeds_path, prepare_contract.REVEAL_HEADER, seed_rows)
        if len(read_selection_seeds(seeds_path)) != 3:
            fail('selftest arbitrary non-holdout usage was rejected')
        broken = temp / 'broken/meta'
        broken.mkdir(parents=True)
        (broken / 'manifest.json').write_text('{}\n')
        expect_failure(lambda: read_jobs(temp / 'broken'), 'empty manifest')
        (broken / 'manifest.json').write_text('{"schema":"wirehair.wh2.normalized_h15_v5.selection.v1"}\n')
        expect_failure(lambda: read_jobs(temp / 'broken'), 'partial manifest')
        telemetry = temp / 'telemetry'
        (telemetry / 'meta').mkdir(parents=True)
        (telemetry / 'meta/run_start.txt').write_text(
            'started_utc=2026-01-01T00:00:00Z\nworkers=120\nthreads_per_invocation=1\n'
            'launcher_nice=0\njob_count=1980\nthermal_csv=/tmp/wirehair-enoq-thermal.csv\n'
            'thermal_line_start=1\nedac_ce_start=0\nedac_ue_start=0\n')
        (telemetry / 'meta/run_finish.txt').write_text(
            'finished_utc=2026-01-01T00:01:00Z\nxargs_rc=0\nstatus_count=1980\n'
            'stdout_count=1980\nstderr_count=1980\nstderr_nonempty=0\ntime_count=1980\n'
            'command_count=1980\nraw_entry_count=5940\nstatus_entry_count=1980\n'
            'command_entry_count=1980\nartifact_names_ok=1\nthermal_line_end=3\n'
            'edac_ce_end=0\nedac_ue_end=0\n')
        thermal_rows = []
        for utc, monotonic, busy in (('2026-01-01T00:00:01Z', '1', '99'),
                                     ('2026-01-01T00:00:02Z', '2', '97')):
            thermal_rows.append([utc, monotonic, busy, '3000', '60',
                                 *(['40'] * 8), '0', '120', '120', '120', '0', '0'])
        thermal_path = telemetry / 'thermal_interval.csv'
        with thermal_path.open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(THERMAL_HEADER)
            writer.writerows(thermal_rows)
        stats = environment_audit(telemetry)
        if stats['cpu_busy_min_pct'] != '97' or stats['cpu_busy_mean_pct'] != '98':
            fail('selftest CPU-busy telemetry reduction failed')
        thermal_rows[1][1] = '12'
        with thermal_path.open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(THERMAL_HEADER)
            writer.writerows(thermal_rows)
        expect_failure(lambda: environment_audit(telemetry), 'telemetry gap above 10 seconds')
        thermal_rows[1][1] = '2'
        thermal_rows[0][2] = '89'
        with thermal_path.open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(THERMAL_HEADER)
            writer.writerows(thermal_rows)
        expect_failure(lambda: environment_audit(telemetry), 'CPU busy minimum below 90 percent')
        thermal_rows[0][2] = '97'
        thermal_rows[1][2] = '98'
        with thermal_path.open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(THERMAL_HEADER)
            writer.writerows(thermal_rows)
        expect_failure(lambda: environment_audit(telemetry), 'CPU busy mean below 98 percent')
        fake_baseline = temp / 'sealed-baseline'
        seal_records = {path: '0' * 64 for path in expected_selection_seal_paths(temp)}
        seal_records[fake_baseline] = BASELINE_SHA256
        if selection_baseline_from_records(temp, seal_records) != fake_baseline.resolve():
            fail('selftest exact selection seal roles failed')
        wrong_role = dict(seal_records)
        del wrong_role[Path(__file__).resolve().parent / 'wh2_h15_v5_selection_job.sh']
        wrong_role[temp / 'arbitrary-extra'] = '0' * 64
        expect_failure(lambda: selection_baseline_from_records(temp, wrong_role),
                       'arbitrary selection-seal role')
        wrong_baseline = dict(seal_records)
        wrong_baseline[fake_baseline] = '0' * 64
        expect_failure(lambda: selection_baseline_from_records(temp, wrong_baseline),
                       'unpinned selection-seal baseline')
    expect_failure(lambda: parse_failure_ids('1|1', 2, 'selftest'), 'duplicate ID')
    expect_failure(lambda: uint('+1', 'selftest'), 'signed uint')
    if json_exact({'jobs': 1980.0}, {'jobs': 1980}):
        fail('selftest JSON numeric type distinction failed')
    require_reveal_binding({'selection_reveal_sha256': 'a', 'selection_seeds_sha256': 'a'})
    expect_failure(lambda: require_reveal_binding({'selection_reveal_sha256': 'a',
                                                    'selection_seeds_sha256': 'b'}),
                   'reveal/canonical seed hash mutation')
    candidates = [[0, 0, 0, 1] for _ in ARMS[1:]]
    state = {7: {'fresh': 1, 'strata': {('D', 's0', 'burst', '0.50'), ('G', 's1', 'burst', '0.50')}, 'candidate': candidates}}
    v4 = {7: 0, 9: 5}
    scores = {(7, salt): True for _, salt in ARMS[1:]}
    if not candidate_gates(7, 0, v4, state, scores, {7})[0]:
        fail('selftest repair/recurrence gate failed')
    candidates[0][1] = 1
    if candidate_gates(7, 0, v4, state, scores, {7})[0] or not candidate_gates(7, 1, v4, state, scores, {7})[0]:
        fail('selftest introduction/priority fallback failed')
    candidates[0][1] = 0
    candidates[0][2] = 1
    if candidate_gates(7, 0, v4, state, scores, {7})[0]:
        fail('selftest candidate error passed')
    candidates[0][2] = 0
    state[7]['strata'] = {('G', 's1', 'burst', '0.50')}
    if candidate_gates(7, 0, v4, state, scores, {7})[0]:
        fail('selftest one-stratum recurrence passed')
    if candidate_gates(9, 0, v4, {}, {}, set())[0]:
        fail('selftest immutable nonzero source passed')
    passing = [1, 1, 10000, 10100, 20000, 20050, 0, 0, 0, 0]
    failing = passing.copy()
    failing[5] = 20051
    if not caps_pass(passing) or caps_pass(failing):
        fail('selftest cap boundary failed')
    base = {'failure_trials': (0,), 'success': 1, 'rank_fail': 1, 'error': 0,
            'xors': 1000, 'muls': 2000}
    repaired = {'failure_trials': (), 'success': 2, 'rank_fail': 0, 'error': 0,
                'xors': 1010, 'muls': 2005}
    introduced = dict(repaired, failure_trials=(1,), success=1, rank_fail=1)
    recovery = [0] * 10
    add_recovery(recovery, base, repaired, 2)
    add_recovery(recovery, repaired, introduced, 2)
    if recovery[7:10] != [1, 1, 0] or paired_p(2, 0) != '0.5':
        fail('selftest recovery reducer/p-value failed')
    work = [0] * 10
    add_work(work, base, repaired, 2)
    if work[6:10] != [1, 2, 1, 0]:
        fail('selftest mixed-row exclusion failed')
    print('selftest: PASS (repair/introduction/error/priority/recurrence/immutability/caps/reducers + 17 mutations)')

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    validate = commands.add_parser('validate-job')
    validate.add_argument('--experiment', type=Path, required=True)
    validate.add_argument('--job-id', type=int, required=True)
    validate.set_defaults(func=validate_job)
    temporary = commands.add_parser('validate-output')
    temporary.add_argument('--experiment', type=Path, required=True)
    temporary.add_argument('--job-id', type=int, required=True)
    temporary.add_argument('--stdout', type=Path, required=True)
    temporary.add_argument('--stderr', type=Path, required=True)
    temporary.set_defaults(func=validate_temp_output)
    final = commands.add_parser('analyze')
    final.add_argument('--experiment', type=Path, required=True)
    final.set_defaults(func=analyze)
    test = commands.add_parser('selftest')
    test.set_defaults(func=selftest)
    args = parser.parse_args()
    try:
        args.func(args)
    except (AuditError, OSError, UnicodeError, ValueError, csv.Error, json.JSONDecodeError) as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 1
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
