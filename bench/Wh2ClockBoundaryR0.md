# Codec-free clock-boundary diagnostic R0

Issue `wirehair-sxvz.16.1.20.82.3.1.4`; protocol
`wirehair.wh2.clock-boundary-r0`, one-shot namespace
`/var/tmp/wh2-clock-boundary-r0`. This is not a codec speed or recovery gate.

The prior complete K5 timing audit found eight large CPU-and-wall excursions
with zero recorded fault/switch deltas across both codecs. Every one started
within the last 50 microseconds of an absolute millisecond and crossed its
boundary. Eighteen codec-free waits also overshot their targets with similar
CPU-charged delays. That association motivates distinguishing clock-read
windows from fixed-count computation; it does not identify an interrupt,
allocator, kernel or clock implementation as the cause.

## Prospectively fixed workload

One pinned CPU50 process, 524,288 consecutive observations. Each observation:

1. Read per-thread rusage and thread CPU time.
2. Take an ordered TSC/AUX stamp, read `CLOCK_MONOTONIC`, stamp again.
3. Execute exactly 65,536 xorshift13/7/17 iterations from `0x9e3779b97f4a7c15`.
4. Stamp, read `CLOCK_MONOTONIC`, stamp again.
5. Read thread CPU time and rusage; verify the result and affinity.

The computation result must be `0xd63028b03dfd593c`, independently reconstructed
with integer arithmetic by the reader. It is an ordinary noinline/noipa
function, not a constant-folded loop. Its native disassembly must contain
register computation only. The diagnostic links no codec; it reuses only the
existing SHA256 helper for authenticating the exact claim bytes.

No phase targeting, waiting, adaptive iteration count, conditioning prelude,
warm exclusion, correction/subtraction or sample filtering. Report all 48
absolute one-millisecond bins, classified by the first returned monotonic
timestamp; all bins must be populated before a complete diagnostic is claimed.
They are not precise pre-call entry times if the clock read itself stalls.
Output uses fixed 256-record chunks outside observation brackets. Startup and
completion CPU/wall anchors and all cumulative fault/switch counters are kept.

The CPU must be AuthenticAMD family 26/model 8/stepping 1/APIC 100, logical CPU50,
without a reported hypervisor. Require CPUID RDTSCP, invariant TSC and
always-serializing LFENCE support before reading TSC; a short actual neutral run
also tests instruction access. Use compiler memory barriers and the bracket
`mfence; lfence; rdtscp; lfence`. Every TSC AUX must be 50 and singleton affinity
must remain CPU50. The ordering rationale follows the
[upstream Linux ordered-TSC implementation](https://github.com/torvalds/linux/blob/v6.8/arch/x86/include/asm/msr.h)
and this host's explicit LFENCE feature. TSC is a different read mechanism,
not an independent physical clock; its ticks are not retired instructions or
core-cycle counts. Clock windows include fences and small capture bookkeeping,
so a delay there does not alone prove a vDSO defect.

Native limits: CPU 100 seconds, wall 120 seconds, address space 128 MiB and
disabled core dumps. Controller: 150-second worker observation limit, raw
192 MiB, stderr 64 KiB. The ASan-only executable may run only 256-observation
neutral modes and omits the virtual-address cap needed by ASan; it can never
enter scientific mode. No host policy or power settings are changed.

## Qualification, retention and interpretation

Before claiming the namespace: strict native and ASan/UBSan builds, complete
256-record neutral checks, independent checksum/clock/AUX/counter/ordinal
validation, negative reader tests, and final-record result/clock failure
injections. Failure records retain completed fields and an explicit capture
stage. Record all compiler commands, dependencies, runtime libraries, producing
sources, executables and neutral outputs; verify pins before and after the one
native launch. Producing sources must be committed before launch.

Raw arrays have 23 unsigned integer fields as documented beside `Record` in
the worker, and a strict header/footer. Four stamps expose three windows:
first clock, computation, second clock. Report medians/maxima and every bin's
count/computation distribution. The full-run TSC/monotonic endpoint ratio gives
descriptive nanosecond equivalents, without modifying raw ticks. Count windows
at least 50 microseconds equivalent and retain full records for either clock
window meeting that label or computation at least 1.5 times its overall median.
Every other observation remains in the immutable raw file and all summaries.
These labels do not exclude records or act as codec acceptance thresholds.

`DIAGNOSTIC_COMPLETE` means the entire frozen sequence, computation, input
bindings and phase coverage passed engineering checks; it does not mean a
cause was established. Delays confined to clock windows and delays seen in
register computation imply different next investigations. Either can contain
unrecorded interference. A null result is inconclusive, not permission to
repeat, resize or avoid boundary phases. K5 cost R0 stays `CONTROL_FAIL`.
No original codec data, equations, source ownership, admission requirement,
recovery sample, all-K objective or speed bound changes.

Build with `python3 bench/Wh2ClockBoundaryR0.py build /absolute/new/build-dir`.
After qualification and commit, use its manifest for the sole `run` command.
`analyze` is read-only; it never invokes the worker.
