# K12 GF(256) structural screen

This is a finite, rank-only feasibility screen selected by the audited
K12 recovery gap. It does not load a codec, measure speed, test a candidate
payload, or change production defaults. The sole output namespace is
`/var/tmp/wh2-k12-thue-morse-r1`; failed output is retained permanently.

The earlier r0 namespace is permanently retained as an INVALID harness-only
attempt: its worker expected a newline-terminated claim while the generic
controller correctly wrote compact canonical JSON. It selected no candidate
and is not evidence. r1 fixes only that pre-worker encoding mismatch.

The screen selects the first candidate in ascending parameter order from a
fixed K12 companion-pair family. It checks all ten binary-factor words and all
1,820 twelve-column minors per word, then checks the packed dyadic mapper at
2,049 sequential IDs, 31 power-of-two seams, and `UINT32_MAX`. It rechecks all
63 original-length failed K12 prefixes from the spent baseline inventory.

Only after those gates pass does it score 512 fresh roots at each of three
widths and four frozen loss schedules (6,144 traces total), delivering K+4
IDs with the same 65,536 candidate cap. OH0 must be at most 1% overall and in
each of the twelve cells. Roots exclude the inventory and earlier K2/K3/K5/K8
and recovery protocols. The screen is GF(256) with polynomial `0x14d`, alpha 2.

The worker is bounded to 60 seconds CPU/wall, 512 MiB address space, zero core
dump and a 4 MiB output cap. It is launched once from `Wh2K12ThueMorseRunR0.py`
after exact pushed-source and retained-inventory authentication. No retry,
reselection, native integration, WH1 comparison, timing, holdout, all-K or
promotion claim is permitted by this screen.
