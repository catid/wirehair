# V2 fuzz targets

Three bounded entry points exercise independent failure surfaces:

- `wirehair_v2_profile_fuzz` mutates precode dimensions and row storage,
  versioned profile fields and salts, packet domains, and sparse solves checked
  against the extracted dense GF(256) oracle.
- `wirehair_v2_stateful_fuzz` mutates packet order and sizes, identical and
  conflicting duplicates, rank-deficient prefixes, allocation retries,
  post-completion calls, malformed-profile retries, and transactional V1/V2
  facade reuse.
- `wirehair_v2_borrowed_facade_fuzz` drives the installed C and C++ facade
  through malformed storage-policy options, all borrowed constructor routes,
  direct and high-ID packets, retained-source overlap failures, detach, source
  release, and transactional C++ replacement/move operations.

Ordinary and ASan/UBSan CTest execute every fixed corpus entry followed by
10,000 deterministic mutations per target.  Inputs are capped at 64 MiB,
deterministic mutations at 64 KiB, block bytes at 4096, and the dense oracle at
K <= 128.  Corpus manifests reject aggregate artifacts above 5 MiB.  Failing
properties write the exact input plus seed/index metadata beneath
`Testing/Temporary/v2-fuzz-artifacts`; `--replay FILE` and
`--seed S --replay-index I` reproduce them.

Example:

```bash
cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --target wirehair_v2_profile_fuzz wirehair_v2_stateful_fuzz wirehair_v2_borrowed_facade_fuzz
build/codec/wirehair_v2_profile_fuzz --mutations 10000 --seed 0x243f6a8885a308d3
build/codec/wirehair_v2_stateful_fuzz --mutations 10000 --seed 0x13198a2e03707344
build/codec/wirehair_v2_borrowed_facade_fuzz --mutations 10000 --seed 0xa4093822299f31d0
```

Clang coverage-guided binaries use `-DWIREHAIR_ENABLE_LIBFUZZER=ON`.  The
scheduled workflow gives each target a pinned seed and at least 1,800 seconds,
with `-max_len=65536`, a 1,024 MiB RSS ceiling, and a persistent
crash-artifact prefix.
