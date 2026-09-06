#!/usr/bin/env python3
"""Strict recorded-native-result reduction; never launches or tunes a codec."""
import hashlib
import importlib.util
import math
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
PROTOCOL = "wirehair.wh2.thue-morse-native-r0"
WIDTHS = (2, 64, 1280)
SCOPES = ("encoder_create", "decoder_create", "prebuilt_systematic",
          "prebuilt_sequential_repair", "prebuilt_distant_repair",
          "receive_recover", "decode_end_to_end")
COMPARISONS = ("T/T", "L/L", "P/P", "T/L", "T/P")
CONDITIONS = ((0, 1, 3, 2), (1, 2, 0, 3), (2, 3, 1, 0), (3, 0, 2, 1))
VISIT_SHA = "b9187d801c8c92e31672dc53fe0b9ccc1c186858196c1f9c00a2ad37028a7d61"
LOOKUP_SHA = "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d"


def expected_config():
    return {"K": 6, "widths": list(WIDTHS), "batch_messages": 64,
            "batch_encode_calls": 384, "scopes": list(SCOPES), "replicates": 12,
            "conditions": 4, "comparisons": [name.split("/") for name in COMPARISONS],
            "invocations_per_slot": 2, "condition_orders": [list(row) for row in CONDITIONS],
            "systematic_ids": list(range(6)), "sequential_ids": list(range(6, 12)),
            "distant_ids": [(1 << 32) - 1 - 2 * i for i in range(6)],
            "receive_ids": list(range(6, 16)), "warmup_order": ["left", "right"],
            "lookup_bytes": 39936, "lookup_sha256": LOOKUP_SHA,
            "message_bytes": [12, 384, 7680], "source_seed": "0x6e61746976653633",
            "source_law": "splitmix64_le_bytes_seed_xor_M_xor_B_shift32"}


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Only its pure, reviewed math and identity validator are reused. No old
# controller, campaign namespace, outcome, batch size or control is inherited.
H = sibling("_tm_native_reviewed_reducers", "Wh2PublicBorrowedCurrentVsWh1R1.py")
H.TARGET_CPU = 50


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, low, high, label):
    require(type(value) is int and low <= value <= high, label)
    return value


def keys(value, expected, label):
    require(type(value) is dict and set(value) == set(expected), label + " keys")


def sha(value, label):
    require(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) is not None, label)


def exact(actual, expected, label):
    # Canonical comparison does not mistake bool for an integer.
    require(H.canonical_bytes(actual) == H.canonical_bytes(expected), label)


def expected_cases(data):
    cases = []
    for index, row in enumerate(data["traces"]):
        first = next((i for i, rank in enumerate(row["ranks"]) if rank == 6), None)
        cases.append(([0, index, row["B"], row["B"]], row["ids"],
                      row["ranks"][-1], first, True))
    for kind, key, collection in ((1, "history_cases", "history"),
                                   (2, "fixture_cases", "fixtures"),
                                   (3, "partial_cases", "traces")):
        for index, (origin, B, tail) in enumerate(data[key]):
            row = data[collection][origin]
            rank = row["ranks"][-1] if kind == 3 else row["rank"]
            # Only the traces have a full recorded prefix-rank sequence.
            first = next((i for i, value in enumerate(row["ranks"]) if value == 6), None) if kind == 3 else None
            cases.append(([kind, index, B, tail], row["ids"], rank, first, kind == 3))
    return cases


def validate_arm(value, packet_count):
    require(type(value) is list and len(value) == 8, "correctness arm shape")
    status = integer(value[0], 0, 7, "correctness arm status")
    extra = value[1]
    if extra is not None:
        integer(extra, 0, packet_count - 6, "first-success overhead")
    for code in value[2:7]:
        if code is not None:
            integer(code, 0, (1 << 31) - 1, "codec result code")
    fed = integer(value[7], 0, packet_count, "fed packet count")
    if status == 0:
        require(extra is not None, "success lacks packet count")
        exact(value, [0, extra, 0, 0, 0, 0, 0, 6 + extra], "success API evidence")
    elif status == 1:
        exact(value, [1, None, 0, 0, 0, 1, None, packet_count], "NeedMore API evidence")
    elif status == 2:
        require(value[2] not in (None, 0), "encoder error code")
        exact(value, [2, None, value[2], None, None, None, None, 0], "encoder error evidence")
    elif status == 3:
        require(value[3] not in (None, 0), "decoder error code")
        exact(value, [3, None, 0, value[3], 0, None, None, 0], "decoder error evidence")
    elif status == 4:
        require(value[4] not in (None, 0), "encode error code")
        exact(value, [4, None, 0, None, value[4], None, None, 0], "encode error evidence")
    elif status == 5:
        require(value[5] not in (None, 0, 1) and fed > 0, "feed error code/count")
        exact(value, [5, None, 0, 0, 0, value[5], None, fed], "feed error evidence")
    elif status == 6:
        require(value[6] not in (None, 0) and extra is not None, "recover error code/count")
        exact(value, [6, extra, 0, 0, 0, 0, value[6], 6 + extra], "recover error evidence")
    else:
        raise ValueError("native packet/oracle/recovered-byte mismatch")
    return status


def validate_correctness(value, data):
    keys(value, ("records", "cases", "candidate_packets", "control_packets",
                 "recovered_messages", "status"), "correctness")
    roster = expected_cases(data)
    records = value["records"]
    require(type(records) is list and len(records) == len(roster), "correctness roster count")
    integer(value["cases"], len(roster), len(roster), "correctness declared cases")
    statuses = [[0] * 8 for _ in range(3)]
    first = [[0] * 6 for _ in range(3)]
    control_errors = 0
    for record, (coordinate, ids, rank, expected_first, exact_first) in zip(records, roster):
        require(type(record) is list and len(record) == 7, "correctness record shape")
        exact(record[:4], coordinate, "correctness ordered coordinate")
        for arm, result in enumerate(record[4:]):
            status = validate_arm(result, len(ids))
            statuses[arm][status] += 1
            if arm == 0:
                require(status == (0 if rank == 6 else 1), "candidate rank/native mismatch")
                if exact_first:
                    exact(result[1], expected_first, "candidate first-success rank mismatch")
            elif status >= 2:
                control_errors += 1
            if coordinate[0] == 0 and coordinate[1] < 6144:
                first[arm][result[1] if status == 0 else 5] += 1
    for field in ("candidate_packets", "control_packets", "recovered_messages"):
        integer(value[field], 0, 3 * sum(len(row[1]) for row in roster), field)
    expected_packets = sum(len(row[1]) for row in roster)
    require(value["candidate_packets"] == expected_packets, "candidate packet oracle coverage")
    if control_errors == 0:
        require(value["control_packets"] == 2 * expected_packets, "control encoded packet coverage")
    require(value["recovered_messages"] == sum(row[0] for row in statuses), "recovered byte coverage")
    integer(value["status"], 1 if control_errors else 0, 1 if control_errors else 0,
            "correctness declared status")
    return {"status_counts": statuses, "fresh_first_success": first,
            "control_errors": control_errors, "cases": len(roster)}


def panel_roster():
    return ((rep, B, scope, comparison, condition)
            for rep in range(12) for B in WIDTHS for scope in range(7)
            for comparison in range(5) for condition in CONDITIONS[rep % 4])


def validate_panels(panels):
    require(type(panels) is list and len(panels) == 5040, "panel roster count")
    samples, physical_by_width, noncomparable, receive_outcomes = {}, {}, [], {}
    for panel, coordinate in zip(panels, panel_roster()):
        require(type(panel) is list and len(panel) == 11, "panel record shape")
        exact(panel[:5], list(coordinate), "panel ordered coordinate")
        rep, B, scope, comparison, condition = coordinate
        status = integer(panel[5], 0, 2, "panel status")
        require(status != 2, "invalid/fatal native panel")
        physical = panel[9]
        require(type(physical) is list and len(physical) == 6, "physical context receipt")
        for address in physical:
            integer(address, 4096, (1 << 64) - 1, "physical address")
            require(address % 4096 == 0, "physical buffer page offset")
        require(len(set(physical)) == 6, "physical contexts not distinct")
        sizes = (6 * B + 16, 10 * (B + 16), 64 * (6 * B + 16)) * 2
        spans = sorted((address, address + size) for address, size in zip(physical, sizes))
        require(all(end <= (1 << 64) for _, end in spans) and
                all(left[1] <= right[0] for left, right in zip(spans, spans[1:])),
                "physical buffer spans wrap or overlap")
        if B in physical_by_width:
            exact(physical, physical_by_width[B], "physical contexts changed")
        else:
            physical_by_width[B] = physical
        receipt = panel[10]
        keys(receipt, ("native_status", "left_code", "right_code", "handles_count",
                       "handles_sha256", "handles_first", "handles_last",
                       "handles_min", "handles_max"), "panel receipt")
        integer(receipt["native_status"], 0, 0, "native scheduler status")
        sha(receipt["handles_sha256"], "ordered allocator digest")
        expected_handles = 10 if scope in (2, 3, 4) else 640
        count = integer(receipt["handles_count"], 0, expected_handles, "allocated handle count")
        for name in ("handles_first", "handles_last", "handles_min", "handles_max"):
            if count == 0:
                require(receipt[name] is None, "empty allocator ledger has addresses")
            else:
                integer(receipt[name], 1, (1 << 64) - 1, name)
        if count:
            for name in ("handles_first", "handles_last"):
                require(receipt["handles_min"] <= receipt[name] <= receipt["handles_max"],
                        "allocator address range")
        if status == 0:
            require(count == expected_handles, "scope handle roster")
        for extra in panel[6:8]:
            if scope in (5, 6) and status == 0:
                integer(extra, 0, 4, "timing first-success overhead")
            elif extra is not None:
                require(scope in (5, 6), "nonreceive panel has decoded overhead")
                integer(extra, 0, 4, "noncomparable decoded overhead")
        for code in (receipt["left_code"], receipt["right_code"]):
            integer(code, 0, (1 << 31) - 1, "panel codec outcome")
        for role, code in zip(COMPARISONS[comparison].split("/"),
                              (receipt["left_code"], receipt["right_code"])):
            require(role != "T" or code == 0, "candidate failure cannot be noncomparable")
        if scope in (5, 6):
            roles = COMPARISONS[comparison].split("/")
            for role, extra, code in zip(roles, panel[6:8],
                                         (receipt["left_code"], receipt["right_code"])):
                key = B, role
                observation = (code, extra)
                if key in receive_outcomes:
                    exact(list(observation), list(receive_outcomes[key]),
                          "same receive trace changed first-success/outcome")
                else:
                    receive_outcomes[key] = observation
        if status == 0:
            require(receipt["left_code"] == receipt["right_code"] == 0,
                    "comparable panel has failed codec outcome")
            require(type(panel[8]) is list and len(panel[8]) == 8, "eight native timings")
            for interval in panel[8]:
                integer(interval, 1, 120 * 1000000000, "bounded native interval")
            samples[coordinate] = H.lane_contrast(panel[8], "ABBA" if condition % 2 == 0 else "BAAB")
        else:
            require(panel[8] is None, "failed panel leaked usable timing")
            require(receipt["left_code"] != 0 or receipt["right_code"] != 0,
                    "noncomparable panel has two successful preflights")
            noncomparable.append(list(coordinate))
    return samples, noncomparable


def validate_counts(raw):
    names = ("callbacks", "encoder_creates", "decoder_creates", "encode_calls",
             "feed_calls", "recover_calls")
    keys(raw["timed_scope_counts"], ("measured", "warm"), "timed counts")
    totals = {group: dict.fromkeys(names, 0) for group in ("measured", "warm")}
    lower = {group: dict.fromkeys(names, 0) for group in ("measured", "warm")}
    complete = all(panel[5] == 0 for panel in raw["panels"])
    for panel in raw["panels"]:
        scope = panel[2]
        for group, callbacks_per_arm in (("measured", 4), ("warm", 1)):
            counts = totals[group]
            before = dict(counts)
            calls = 2 * callbacks_per_arm
            counts["callbacks"] += calls
            if scope == 0:
                counts["encoder_creates"] += calls * 64
            if scope in (1, 6):
                counts["decoder_creates"] += calls * 64
            if scope in (2, 3, 4):
                counts["encode_calls"] += calls * 384
            if scope in (5, 6):
                counts["recover_calls"] += calls * 64
                if panel[5] == 0:
                    counts["feed_calls"] += callbacks_per_arm * 64 * (12 + panel[6] + panel[7])
                else:
                    counts["feed_calls"] += calls * 64 * 10
            if panel[5] == 0:
                for name in names:
                    lower[group][name] += counts[name] - before[name]
            else:
                lower[group]["callbacks"] += calls
                # Encoder-create deliberately attempts every one of64 even if
                # an earlier creation in that callback failed.
                if scope == 0:
                    lower[group]["encoder_creates"] += calls * 64
    for group in totals:
        observed = raw["timed_scope_counts"][group]
        keys(observed, names, group + " scope counts")
        for name in names:
            integer(observed[name], lower[group][name], totals[group][name], group + " " + name)
        if complete:
            exact(observed, totals[group], group + " exact scope call counts")
        else:
            exact(observed["callbacks"], totals[group]["callbacks"], "noncomparable callback completeness")
    measured_ns = sum(sum(panel[8]) for panel in raw["panels"] if panel[5] == 0)
    require(measured_ns <= raw["elapsed_seconds"] * 1000000000 + 1,
            "timed intervals exceed whole worker lifetime")


def statistics(samples):
    summaries, failed = [], []
    boundary = math.log1p(0.02)
    for B in WIDTHS:
        for scope in range(7):
            for comparison in range(5):
                conditions = []
                for condition in range(4):
                    coordinates = [(rep, B, scope, comparison, condition) for rep in range(12)]
                    summary = H.sample_summary([samples[key] for key in coordinates]) if all(
                        key in samples for key in coordinates) else None
                    conditions.append(summary)
                    if summary is not None and comparison < 3 and not (-boundary < summary["lower95_log"] and
                                                summary["upper95_log"] < boundary):
                        failed.append([B, scope, comparison, condition])
                matched = H.sample_summary([
                    math.fsum(samples[(rep, B, scope, comparison, condition)]
                              for condition in range(4)) / 4 for rep in range(12)]) if all(
                                  row is not None for row in conditions) else None
                if matched is not None and comparison < 3 and not (-boundary < matched["lower95_log"] and
                                            matched["upper95_log"] < boundary):
                    failed.append([B, scope, comparison, "matched"])
                summaries.append({"B": B, "scope": SCOPES[scope],
                                  "comparison": COMPARISONS[comparison],
                                  "conditions": conditions, "matched": matched,
                                  "all_upper95_below_one": all(row is not None and row["upper95"] < 1 for row in
                                                              conditions + [matched])})
    return summaries, failed


def reduce_raw(raw, data=None, deadline=None):
    keys(raw, ("schema", "protocol", "target_cpu", "gf_runtime", "visit_sha256",
               "data_sha256", "correctness", "panels", "elapsed_seconds", "outcome",
               "diagnostic", "target_identity_before", "target_identity_after", "config",
               "timed_scope_counts"), "native raw")
    exact(raw["schema"], PROTOCOL + ".raw.v1", "raw schema")
    exact(raw["protocol"], PROTOCOL, "raw protocol")
    exact(raw["config"], expected_config(), "native frozen configuration")
    integer(raw["target_cpu"], 50, 50, "worker CPU")
    require(raw["outcome"] in ("COMPLETE", "NONCOMPARABLE") and raw["diagnostic"] == "",
            "native fatal/partial result")
    require(type(raw["elapsed_seconds"]) in (int, float) and
            math.isfinite(raw["elapsed_seconds"]) and 0 < raw["elapsed_seconds"] <= 120,
            "native elapsed bound")
    for key in ("target_identity_before", "target_identity_after"):
        H.validate_target_identity(raw[key], key)
    exact(raw["target_identity_before"], raw["target_identity_after"], "CPU identity drift")
    runtime = raw["gf_runtime"]
    keys(runtime, ("polynomial", "address", "ssse3", "avx2", "gfni", "avx512"), "GF runtime")
    integer(runtime["polynomial"], 0x14d, 0x14d, "GF polynomial")
    integer(runtime["address"], 1, (1 << 64) - 1, "GF runtime address")
    for key in ("ssse3", "avx2", "gfni", "avx512"):
        integer(runtime[key], 0, 1, "GF selected feature")
    exact(raw["visit_sha256"], VISIT_SHA, "native coefficient visit oracle")
    if data is None:
        D = sibling("_tm_native_recorded_data", "Wh2ThueMorseNativeDataR0.py")
        data = D.extract_data(D.authenticated_inputs(deadline=deadline))
    exact(raw["data_sha256"], hashlib.sha256(H.canonical_bytes(data)).hexdigest(),
          "native data input seal")
    correctness = validate_correctness(raw["correctness"], data)
    samples, noncomparable = validate_panels(raw["panels"])
    validate_counts(raw)
    summaries, failed = statistics(samples)
    exact(raw["outcome"], "NONCOMPARABLE" if noncomparable or correctness["control_errors"] else "COMPLETE",
          "raw declared outcome differs from evidence")
    outcome = "NONCOMPARABLE" if noncomparable or correctness["control_errors"] else (
        "CONTROL_FAIL" if failed else "PASS")
    # PASS means a valid short cost measurement, not all scopes faster than WH1.
    return {"protocol": PROTOCOL, "outcome": outcome, "correctness": correctness,
            "statistics": summaries, "failed_controls": failed,
            "noncomparable_panels": noncomparable, "WH1_compared": True,
            "promotion_claimed": False, "public_candidate_api_claimed": False,
            "all_K_claimed": False}
