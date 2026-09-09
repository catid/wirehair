# Certified fall-through dispatch: isolated candidate

Tracked in `wirehair-sxvz.16.1.20.84.3.3`. Production is unchanged.

The compile-only first screen replaces four small-profile dispatch conditions
in a hash-pinned copy of `WirehairV2Profile.cpp` with the existing portable
`CAT_UNLIKELY` macro. No equation, allocation, object size, public contract,
per-K implementation, constructor policy or default changes. The baseline
object must reproduce the qualified native object byte-for-byte.

Native assembly shows the intended certified decode fall-through and removal
of its unconditional return-path merge jump. K3/K5 instead take the out-of-line
branch and merge jump: their cost must be measured, not assumed unchanged.
Encode shrinks 602 to 595 bytes, decode stays 118, recover grows 283 to 308,
and free grows 49 to 54. These are code-layout facts, not a speed result.

Configure with `WH2_FALLTHROUGH_NEUTRAL=ON` only after the compile/assembly
screen passes. It builds the existing public K3/K5 correctness tests against
the native, portable-arithmetic and fully instrumented ASan/UBSan archives,
replacing only the profile object. Native baseline/candidate certified-byte
emitters are also built. Check every link map to ensure the archive's original
profile object was not extracted; compare both emitters' complete output and
return codes. Do not run CTest discovery in the preserved archive directories.

No scientific worker, timing run, new recovery sample or promotion is included.
Retaining this candidate requires a separately frozen full-lifecycle speed
screen that includes ordinary K3/K5 retention as well as certified paths and
same-code controls. Existing valid regressions and failed-control results stay
spent; unchanged code placement or faster assembly is not inherited speed.
