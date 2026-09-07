#!/usr/bin/env python3
"""Independent synthetic records only; never starts a codec workload."""
import contextlib
import copy
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock
import Wh2OrderedRhsReleaseR0 as M


def fixture(shrinks=3):
    claim = "a"*64
    rows = [dict(type="header", protocol=M.PROTOCOL, claim_sha256=claim, pid=123, speed_claimed=False)]
    old = dict(fixtures=[]); stderr = []; trace = ["123 brk(NULL) = 0x100000"]
    index = 0; address = 0x100000
    for width in (2, 64, 1280):
        profile, packets = "00"*32, "00"*(18*width)
        rows.append(dict(type="fixture", width=width, profile_hex=profile, packets_hex=packets))
        old["fixtures"].append(dict(handles=[dict(arm=3, profile_hex=profile)], packets_hex=[None, None, None, packets]))
        for family in range(2):
            for cycle in range(3):
                ops = [("create", 0, 0)]
                ops += [("feed", 6+i if family == 0 else (1 << 32)-1-2*i, 0 if i == 5 else 1) for i in range(6)]
                ops += [("recover", 0, 0), ("free", 0, 0)]
                for op, packet, result in ops:
                    index += 1
                    rows.append(dict(type="phase", index=index, operation=op, width=width, family=family,
                                     cycle=cycle, packet=packet, result=result, before=[0]*4, after=[0]*4))
                    for side in ("BEGIN", "END"):
                        text = "WH2_RELEASE_PHASE %d %s\n" % (index, side)
                        stderr.append(text)
                        trace.append('123 write(2, "%s", %d) = %d' % (text.replace("\n", "\\n"), len(text), len(text)))
                        if side != "BEGIN" or width != 1280:
                            continue
                        if op == "feed" and result == 0:
                            address += 5*4096
                            trace.append("123 brk(%s) = %s" % (hex(address), hex(address)))
                            for _ in range(shrinks):
                                address -= 4096
                                trace.append("123 brk(%s) = %s" % (hex(address), hex(address)))
                        elif op == "free":
                            address = 0x100000
                            trace.append("123 brk(%s) = %s" % (hex(address), hex(address)))
    rows.append(dict(type="footer", outcome="COMPLETE", phases=index))
    return rows, "".join(stderr).encode(), ("\n".join(trace)+"\n").encode(), claim, old


def decode_fixture(data):
    rows, error, trace, claim, old = data
    return M.verify(b"".join(M.A.canonical(r) for r in rows), error, trace, claim, old)


class ReleaseTest(unittest.TestCase):
    def test_independent_full_ledger_and_decision(self):
        b = decode_fixture(fixture(3)); c = decode_fixture(fixture(1))
        self.assertEqual((b["phases"], len(b["chronology"])), (162, 162))
        self.assertEqual([x["solve_shrinks"] for x in b["cells"]], [3]*6)
        self.assertEqual([x["lifecycle_vm_calls"] for x in b["cells"]], [5]*6)
        self.assertEqual([x["lifecycle_vm_calls"] for x in c["cells"]], [3]*6)
        self.assertEqual(M.decide([b, c, c, b]), "PASS")

    def test_controls_and_no_subset_rescue(self):
        b = decode_fixture(fixture(3)); c = decode_fixture(fixture(1))
        for arm in (0, 1, 2, 3):
            rs = [copy.deepcopy(r) for r in (b, c, c, b)]
            rs[arm]["cells"][0]["solve_shrinks"] += 1
            self.assertEqual(M.decide(rs), "INCONCLUSIVE")
        self.assertEqual(M.decide([c, c, c, c]), "INCONCLUSIVE")
        self.assertEqual(M.decide([b, b, b, b]), "FAIL")
        for field in ("solve_shrinks", "lifecycle_vm_calls"):
            d = copy.deepcopy(c)
            for x in d["cells"]:
                if x["family"] == 1: x[field] = 6
            self.assertEqual(M.decide([b, d, d, b]), "FAIL")

    def test_corrupt_worker_records(self):
        mutations = [(0, "claim_sha256", "b"*64), (0, "pid", True), (-1, "phases", 163),
                     (1, "packets_hex", "ff"), (2, "index", True), (3, "packet", 42),
                     (8, "result", False), (2, "before", [1, 0, 0, 0]), (2, "after", [-1, 0, 0, 0])]
        for index, key, value in mutations:
            with self.subTest(index=index, key=key):
                data = fixture(); data[0][index][key] = value
                with self.assertRaises((ValueError, KeyError, IndexError, TypeError)):
                    decode_fixture(data)

    def test_marker_and_syscall_corruption(self):
        data = fixture(); rows, error, trace, claim, old = data
        for bad in (trace.replace(b"123 brk", b"124 brk", 1),
                    trace.replace(b"0x100000", b"oops", 1),
                    trace.replace(b"PHASE 1 BEGIN", b"PHASE 2 BEGIN", 1),
                    trace.replace(b"PHASE 162 END", b"PHASE 161 END", 1),
                    trace.replace(b"= 26\n", b"= 25\n", 1), b"x"*(4*1024**2+1)):
            with self.subTest(prefix=bad[:50]), self.assertRaises(ValueError):
                decode_fixture((rows, error, bad, claim, old))
        with self.assertRaises(ValueError):
            decode_fixture((rows, error[:-1], trace, claim, old))

    def test_native_pin_policy(self):
        with mock.patch.object(M.A, "pin", return_value={}) as pin:
            for path, installed in ((Path("/usr/bin/bash"), True), (M.ROOT / M.SOURCES[0], False),
                                    (Path("/tmp/worker"), False)):
                M.pin_input(path); pin.assert_called_with(path, installed=installed)

    def test_capture_timeout_and_exit(self):
        for timeout in (False, True):
            child = mock.Mock(pid=123456789, returncode=None if timeout else 1)
            child.communicate.side_effect = [subprocess.TimeoutExpired("neutral", 35), (b"prefix", b"err")] if timeout else [(b"prefix", b"err")]
            with mock.patch.object(M.subprocess, "Popen", return_value=child), mock.patch.object(M.os, "killpg") as killed:
                raw, error, failure = M.capture(["never-executed"])
                self.assertEqual((raw, error), (b"prefix", b"err")); self.assertIsNotNone(failure)
                if timeout: self.assertGreaterEqual(killed.call_count, 1)
                else: killed.assert_not_called()

    def test_failure_artifact_and_spent_namespace(self):
        with tempfile.TemporaryDirectory(prefix="wh2-release-neutral.") as name:
            base = Path(name); receipt = base / "receipt.json"; output = base / "artifacts"
            M.A.publish(receipt, M.A.canonical(dict(build=str(base))))
            with mock.patch.object(M, "OUTPUT", output), mock.patch.object(M, "current"), \
                    mock.patch.object(M, "capture", return_value=(b"prefix", b"err", "injected")), \
                    contextlib.redirect_stdout(io.StringIO()):
                M.run(receipt)
                self.assertEqual((output / "0-raw.jsonl").read_bytes(), b"prefix")
                self.assertEqual(M.A.decode((output / "COMPLETE.json").read_bytes())["outcome"], "INVALID")
                with self.assertRaises(FileExistsError): M.run(receipt)


if __name__ == "__main__":
    unittest.main()
