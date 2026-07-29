#!/usr/bin/env python3
"""Pure post-v3 classification for repair-v1 normalized cell evidence."""

import collections
import hashlib
import struct


REPAIR_V1_ATTEMPTS = 8
EXPECTED_ERROR_LOCATIONS = frozenset({
    "roles.raw.encode_result",
    "attempts.0.real_result",
})
NON_RUNTIME_RESULT_CODES = frozenset((-1, 0, 1, 3, 4, 7, 9))


class ClassificationError(ValueError):
    """Raised when normalized evidence cannot support classification."""


def _decode_rows(columns, rows, expected_columns, expected_rows, label):
    if (
        columns != list(expected_columns) or
        not isinstance(rows, list) or
        len(rows) != expected_rows or
        any(
            not isinstance(row, list) or len(row) != len(columns)
            for row in rows
        )
    ):
        raise ClassificationError(f"{label} compact rows are malformed")
    return [dict(zip(columns, row)) for row in rows]


def _decode_attempts(real, parser_module):
    return _decode_rows(
        real.get("attempt_columns"),
        real.get("attempt_rows"),
        parser_module.REPAIRTIMING_ATTEMPT_FIELDS,
        REPAIR_V1_ATTEMPTS,
        "post-v3 attempt",
    )


def _decode_timing(real, parser_module):
    return _decode_rows(
        real.get("timing_columns"),
        real.get("timing_rows"),
        parser_module.REPAIRTIMING_TIMING_FIELDS,
        len(parser_module.REPAIRTIMING_PANELS) * 8,
        "post-v3 timing",
    )


def _descriptor_is_bound(selector, real, selected_attempt):
    descriptor_sha256 = selector.get("descriptor_sha256")
    descriptor_hex = real.get("descriptor_hex")
    if (
        not isinstance(descriptor_sha256, str) or
        len(descriptor_sha256) != 64 or
        any(character not in "0123456789abcdef"
            for character in descriptor_sha256) or
        not isinstance(descriptor_hex, str) or
        len(descriptor_hex) != 34
    ):
        return False
    try:
        descriptor = bytes.fromhex(descriptor_hex)
    except ValueError:
        return False
    if (
        len(descriptor) != 17 or
        descriptor[:4] != b"W2R3" or
        descriptor[-1] != selected_attempt or
        hashlib.sha256(descriptor).hexdigest() != descriptor_sha256
    ):
        return False
    arm_id, construction_root = struct.unpack("<QI", descriptor[4:16])
    return (
        arm_id == selector.get("arm_id") and
        construction_root == selector.get("construction_root")
    )


def _later_success_is_bound(selector, real, roles, attempts):
    selected = selector.get("selected_attempt")
    if type(selected) is not int or not 1 <= selected < REPAIR_V1_ATTEMPTS:
        return False
    selected_attempt = attempts[selected]
    selector_attempts = selector.get("attempts")
    if (
        not isinstance(selector_attempts, list) or
        len(selector_attempts) != REPAIR_V1_ATTEMPTS or
        not isinstance(selector_attempts[selected], dict) or
        selector_attempts[selected].get("attempt") != selected or
        selector_attempts[selected].get("probe_executed") != 1 or
        selector_attempts[selected].get("probe_result") != 0 or
        selected_attempt.get("attempt") != selected or
        selected_attempt.get("probe_executed") != 1 or
        selected_attempt.get("probe_result") != 0 or
        selected_attempt.get("real_executed") != 1 or
        selected_attempt.get("real_result") != 0 or
        selector.get("selector_result") != 0 or
        selector.get("committed") != 1 or
        selector.get("cap_exhausted") != 0 or
        selector.get("fatal_attempt_zero_mismatch") != 0 or
        selector.get("oom") != 0 or
        not _descriptor_is_bound(selector, real, selected)
    ):
        return False
    repaired = roles.get("repaired")
    return (
        isinstance(repaired, dict) and
        repaired.get("encode_result") == 0 and
        repaired.get("decode_construct_result") == 0 and
        repaired.get("feed_result") == 0 and
        repaired.get("recover_result") == 0 and
        repaired.get("recovery_ok") == 1 and
        repaired.get("outcome_class") == "success"
    )


def _attempt_zero_classification(selector, real, roles, attempts):
    raw = roles.get("raw") if isinstance(roles, dict) else None
    if (
        not isinstance(raw, dict) or
        not isinstance(attempts, list) or
        len(attempts) != REPAIR_V1_ATTEMPTS or
        not isinstance(selector, dict)
    ):
        raise ClassificationError(
            "post-v3 attempt-zero evidence is malformed")
    attempt_zero = attempts[0]
    selector_attempts = selector.get("attempts")
    if (
        attempt_zero.get("attempt") != 0 or
        attempt_zero.get("real_executed") != 1 or
        not isinstance(selector_attempts, list) or
        len(selector_attempts) != REPAIR_V1_ATTEMPTS or
        not isinstance(selector_attempts[0], dict)
    ):
        raise ClassificationError(
            "post-v3 attempt-zero evidence is incomplete")
    selector_attempt_zero = selector_attempts[0]
    for field in ("attempt", "probe_executed", "probe_result"):
        if selector_attempt_zero.get(field) != attempt_zero.get(field):
            raise ClassificationError(
                "selector and real attempt-zero probe evidence disagree")

    raw_result = raw.get("encode_result")
    real_result = attempt_zero.get("real_result")
    if (
        type(raw_result) is not int or
        raw_result != real_result or
        raw_result < 0 or raw_result > 10
    ):
        raise ClassificationError(
            "raw role and attempt-zero real result disagree")

    corroborated_error = (
        raw_result == 8 and
        attempt_zero.get("probe_executed") == 1 and
        attempt_zero.get("probe_result") == 1 and
        selector.get("fatal_attempt_zero_mismatch") == 0 and
        _later_success_is_bound(selector, real, roles, attempts)
    )
    if raw_result in (3, 4):
        kind = "explicit_weak"
    elif raw_result == 1:
        kind = "need_more"
    elif corroborated_error:
        kind = "corroborated_error"
    else:
        kind = "none"
    return kind, (
        EXPECTED_ERROR_LOCATIONS if corroborated_error else frozenset()
    )


def _result_observations(selector, roles, controls, attempts, timing_rows):
    yield "selector.selector_result", selector.get("selector_result")
    for role_name in ("raw", "repaired", "dispatch", "wh1"):
        role = roles.get(role_name)
        if not isinstance(role, dict):
            raise ClassificationError(
                "post-v3 recovery role is malformed")
        for field in (
                "encode_result", "decode_construct_result",
                "feed_result", "recover_result"):
            yield f"roles.{role_name}.{field}", role.get(field)
    if not isinstance(controls, dict):
        raise ClassificationError(
            "post-v3 control evidence is malformed")
    for field in (
            "forced_a_result", "forced_b_result",
            "repair_direct_result", "dispatch_direct_result"):
        yield f"controls.{field}", controls.get(field)
    for attempt_index, attempt in enumerate(attempts):
        for kind in ("probe", "real"):
            if attempt.get(f"{kind}_executed") == 1:
                yield (
                    f"attempts.{attempt_index}.{kind}_result",
                    attempt.get(f"{kind}_result"),
                )
    for row in timing_rows:
        for field in (
                "timing_construct_result", "timing_result",
                "timing_recover_result"):
            yield f"timing.{field}", row.get(field)


def classify_repair_v1_cell(selector, real, *, parser_module):
    """Classify one fully replayed normalized repairtiming-v3 cell.

    The two code-8 observations for a corroborated attempt-zero structural
    failure are duplicate views of one expected solver outcome.  No other
    code-8 location is exempted.
    """
    if (
        not isinstance(real, dict) or
        real.get("schema") !=
            "wirehair.wh2.repairtiming.real-witness.v3"
    ):
        raise ClassificationError(
            "post-v3 real witness schema is invalid")
    roles = real.get("roles")
    attempts = _decode_attempts(real, parser_module)
    timing_rows = _decode_timing(real, parser_module)
    kind, exemptions = _attempt_zero_classification(
        selector, real, roles, attempts)
    repaired = roles.get("repaired") if isinstance(roles, dict) else None
    if not isinstance(repaired, dict):
        raise ClassificationError(
            "post-v3 repaired role is malformed")
    repaired_final_weak = int(
        repaired.get("outcome_class") == "weak")

    runtime_observations = collections.Counter()
    for location, code in _result_observations(
            selector, roles, real.get("controls"), attempts, timing_rows):
        if type(code) is not int or not -1 <= code <= 10:
            raise ClassificationError(
                f"post-v3 result code is malformed at {location}")
        if code == 8 and location in exemptions:
            continue
        if code not in NON_RUNTIME_RESULT_CODES:
            runtime_observations[f"{code}:{location}"] += 1
    return {
        "raw_attempt0_kind": kind,
        "raw_attempt0_structural_weak": int(kind != "none"),
        "repaired_final_weak": repaired_final_weak,
        "candidate_runtime_error": int(bool(runtime_observations)),
        "runtime_error_observations": dict(
            sorted(runtime_observations.items())),
    }


def aggregate_selector_classifications(
        records, *, schedules, first_schedule):
    """Verify the schedule panel and count each selector key exactly once.

    Each record contains ``selector_key``, ``schedule``, ``replicate``, and
    ``classification``.  Callers remain responsible for accumulating runtime
    errors over every physical cell; this helper intentionally handles only
    the selector-key-weighted structural metric.
    """
    if (
        not isinstance(records, list) or
        not isinstance(schedules, (list, tuple)) or
        not schedules or len(set(schedules)) != len(schedules) or
        any(not isinstance(schedule, str) or not schedule
            for schedule in schedules) or
        first_schedule != schedules[0]
    ):
        raise ClassificationError(
            "selector classification schedule contract is invalid")
    expected_schedules = set(schedules)
    groups = {}
    for record in records:
        if (
            not isinstance(record, dict) or
            set(record) != {
                "selector_key", "schedule", "replicate", "classification"
            } or
            not isinstance(record["selector_key"], (list, tuple)) or
            len(record["selector_key"]) != 3 or
            record["schedule"] not in expected_schedules or
            type(record["replicate"]) is not int or
            record["replicate"] < 0 or
            not isinstance(record["classification"], dict)
        ):
            raise ClassificationError(
                "selector classification record is malformed")
        key = tuple(record["selector_key"])
        try:
            hash(key)
        except TypeError:
            raise ClassificationError(
                "selector classification key is unhashable")
        schedule_records = groups.setdefault(key, {})
        if record["schedule"] in schedule_records:
            raise ClassificationError(
                "selector classification repeats a schedule")
        classification = record["classification"]
        kind = classification.get("raw_attempt0_kind")
        weak = classification.get("raw_attempt0_structural_weak")
        repaired_final_weak = classification.get(
            "repaired_final_weak")
        if (
            kind not in (
                "none", "explicit_weak", "need_more",
                "corroborated_error",
            ) or
            weak not in (0, 1) or
            weak != int(kind != "none") or
            repaired_final_weak not in (0, 1)
        ):
            raise ClassificationError(
                "selector structural classification is malformed")
        schedule_records[record["schedule"]] = {
            "kind": kind,
            "weak": weak,
            "repaired_final_weak": repaired_final_weak,
            "replicate": record["replicate"],
        }

    by_kind = {
        "explicit_weak": 0,
        "need_more": 0,
        "corroborated_error": 0,
    }
    for schedule_records in groups.values():
        if set(schedule_records) != expected_schedules:
            raise ClassificationError(
                "selector classification lacks a frozen schedule")
        if len({
                record["replicate"]
                for record in schedule_records.values()
        }) != 1:
            raise ClassificationError(
                "selector witness replicate changed by schedule")
        structural_signatures = {
            (
                record["kind"], record["weak"],
                record["repaired_final_weak"],
            )
            for record in schedule_records.values()
        }
        if len(structural_signatures) != 1:
            raise ClassificationError(
                "selector structural classification changed by schedule")
        witness = schedule_records[first_schedule]
        if witness["kind"] != "none":
            by_kind[witness["kind"]] += 1
    return {
        "unique_selectors": len(groups),
        "raw_attempt0_structural_weak": {
            "total": sum(by_kind.values()),
            "by_kind": by_kind,
        },
        "repaired_final_weak": {
            "unique_selectors": sum(
                next(iter(schedule_records.values()))[
                    "repaired_final_weak"]
                for schedule_records in groups.values()
            ),
        },
    }
