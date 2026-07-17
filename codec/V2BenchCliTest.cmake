if(NOT DEFINED BENCH)
    message(FATAL_ERROR "BENCH is required")
endif()

function(run_bench result_var out_var err_var)
    execute_process(
        COMMAND "${BENCH}" ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 90)
    set(${result_var} "${result}" PARENT_SCOPE)
    set(${out_var} "${out}" PARENT_SCOPE)
    set(${err_var} "${err}" PARENT_SCOPE)
endfunction()

function(reject_sanitizer output context)
    if("${output}" MATCHES
       "AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer|MemorySanitizer|ThreadSanitizer|runtime error:")
        message(FATAL_ERROR
            "sanitizer failure during ${context}\n${output}")
    endif()
endfunction()

function(expect_failure pattern)
    run_bench(result out err ${ARGN})
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1)
        message(FATAL_ERROR
            "expected exit 1, got '${result}': ${ARGN}\n"
            "stdout=${out}\nstderr=${err}")
    endif()
    if(NOT err MATCHES "${pattern}")
        message(FATAL_ERROR
            "missing diagnostic '${pattern}': ${ARGN}\nstderr=${err}")
    endif()
    reject_sanitizer("${out}${err}" "expected failure: ${ARGN}")
endfunction()

function(expect_success pattern)
    run_bench(result out err ${ARGN})
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0)
        message(FATAL_ERROR
            "expected exit 0, got '${result}': ${ARGN}\n"
            "stdout=${out}\nstderr=${err}")
    endif()
    if(NOT out MATCHES "${pattern}")
        message(FATAL_ERROR
            "missing output '${pattern}': ${ARGN}\nstdout=${out}")
    endif()
    reject_sanitizer("${out}${err}" "expected success: ${ARGN}")
endfunction()

function(expect_output_sha256 expected)
    set(output_file
        "${CMAKE_CURRENT_BINARY_DIR}/wirehair-v2-bench-output-${expected}.tmp")
    file(REMOVE "${output_file}")
    execute_process(
        COMMAND "${BENCH}" ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_FILE "${output_file}"
        ERROR_VARIABLE err
        TIMEOUT 90)
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0)
        file(REMOVE "${output_file}")
        message(FATAL_ERROR
            "expected exit 0, got '${result}': ${ARGN}\nstderr=${err}")
    endif()
    reject_sanitizer("${err}" "SHA-256 output command: ${ARGN}")
    file(SHA256 "${output_file}" actual)
    file(REMOVE "${output_file}")
    if(NOT actual STREQUAL expected)
        message(FATAL_ERROR
            "output SHA-256 mismatch for ${ARGN}: "
            "got ${actual}, expected ${expected}")
    endif()
endfunction()

# The canonical exhaustive dump is byte-identical to the independently
# verified sealed K=2..64000 selection output.  Lookup goldens cover inherited
# v4, retained zero, and every selected v5 salt arm.
expect_output_sha256(
    ef1729d961f71bc604e972194495e9e508459d62db786d0d9df56c3ce4198f3c
    h15table --dump)
expect_success(
    "domain=2..64000 entries=63999 v4_nonzero=70 overrides=98 c27=75 c79=17 c6f=6 ca8=0 logical_table_bytes=756 max_bucket_entries=26.*completion_seal_sha256=4fdc34d4ff18bbbe5e037bd159bbb8c04d0e58e762cc065e5823ce5cf43495e8"
    h15table --summary)
expect_success("4[\t]0x39[\t]0x39[\t]0[\t]immutable-nonzero-v4"
    h15table --lookup 4)
expect_success("42[\t]0x0[\t]0x0[\t]0[\t]retain-zero"
    h15table --lookup 42)
expect_success("955[\t]0x0[\t]0x27[\t]1[\t]selected"
    h15table --lookup 955)
expect_success("39284[\t]0x0[\t]0x79[\t]1[\t]selected"
    h15table --lookup 39284)
expect_success("44131[\t]0x0[\t]0x6f[\t]1[\t]selected"
    h15table --lookup 44131)
expect_failure("requires exactly one" h15table)
expect_failure("requires exactly one" h15table --summary --dump)
expect_failure("--N values must be in" h15table --lookup 1)
expect_failure("--N values must be in" h15table --lookup 64001)
expect_failure("benchmark-lookups must be" h15table --benchmark-lookups 0)
expect_failure("benchmark-lookups must be"
    h15table --benchmark-lookups 1000000001)
expect_failure("--seed requires" h15table --summary --seed 1)
expect_failure("--table requires a lookup benchmark mode"
    h15table --summary --table normalized-h15-v4)
expect_failure("unknown --table" h15table --benchmark-lookups 1
    --table unknown)
expect_failure("--table may be specified only once"
    h15table --benchmark-lookups 1 --table normalized-h15-v4
    --table normalized-h15-v5)
expect_failure("--benchmark-lookups may be specified only once"
    h15table --benchmark-lookups 1 --benchmark-lookups 1)
expect_failure("--seed may be specified only once"
    h15table --benchmark-lookups 1 --seed 1 --seed 1)
expect_failure("--lookup may be specified only once"
    h15table --lookup 64001 --lookup 2)
expect_failure("requires exactly one" h15table --benchmark-lookups 1
    --permutation-step 40501 --domain-repeats 1)
expect_failure("requires both" h15table --permutation-step 40501)
expect_failure("requires both" h15table --domain-repeats 1)
foreach(step IN ITEMS 0 3 63999)
    expect_failure("coprime to 63999" h15table
        --permutation-step ${step} --domain-repeats 1)
endforeach()
expect_failure("--domain-repeats must be nonzero" h15table
    --permutation-step 40501 --domain-repeats 0)
expect_failure("repeat count overflows" h15table
    --permutation-step 40501 --domain-repeats 18446744073709551615)
expect_failure("lookups must be at most" h15table
    --permutation-step 40501 --domain-repeats 15626)
expect_failure("--seed requires --benchmark-lookups" h15table
    --permutation-step 40501 --domain-repeats 1 --seed 1)
expect_failure("--permutation-step may be specified only once" h15table
    --permutation-step 40501 --permutation-step 40501 --domain-repeats 1)
expect_failure("--domain-repeats may be specified only once" h15table
    --permutation-step 40501 --domain-repeats 1 --domain-repeats 1)

# Runtime-generated keys, separate noinline lookup boundaries, and printed
# loop-carried checksums prevent the optimizer from deleting or precomputing
# either lookup-latency arm.  Both arms traverse identical K streams.
expect_success(
    "profile=random-keys arm=normalized-h15-v4 lookups=4096 seed=0x123456789abcdef0 checksum=0x78d1f4d6c778e12d.*selection_table_sha256=ef1729d961f71bc604e972194495e9e508459d62db786d0d9df56c3ce4198f3c"
    h15table --benchmark-lookups 4096 --seed 0x123456789abcdef0
    --table normalized-h15-v4)
expect_success(
    "profile=random-keys arm=normalized-h15-v5 lookups=4096 seed=0x123456789abcdef0 checksum=0x6bcdf5e4c778e135.*selection_table_sha256=ef1729d961f71bc604e972194495e9e508459d62db786d0d9df56c3ce4198f3c"
    h15table --benchmark-lookups 4096 --seed 0x123456789abcdef0
    --table normalized-h15-v5)

# The frozen, seed-free holdout panel precomputes the full-domain permutation
# outside the timed region, then performs division-free nested repeats.  A
# one-repeat golden covers exact order/count/checksum without slowing CI.
expect_success(
    "profile=full-domain-permutation arm=normalized-h15-v4 lookups=63999 permutation_step=40501 domain_repeats=1 checksum=0xdd572d1ed95a951b"
    h15table --permutation-step 40501 --domain-repeats 1
    --table normalized-h15-v4)
expect_success(
    "profile=full-domain-permutation arm=normalized-h15-v5 lookups=63999 permutation_step=40501 domain_repeats=1 checksum=0x3259c1b3d5cd8df4"
    h15table --permutation-step 40501 --domain-repeats 1
    --table normalized-h15-v5)

# Trial count boundaries, including the old uint16 narrowing boundary.
expect_failure("trials must be" seedtable --N 2 --bb-list 1
    --peel-candidates 1 --trials 0)
foreach(trials IN ITEMS 1 65535 65536 1000000)
    expect_success(",${trials},${trials}," seedtable --N 2 --bb-list 1
        --peel-candidates 1 --trials ${trials})
endforeach()
expect_failure("trials must be" seedtable --N 2 --bb-list 1
    --peel-candidates 1 --trials 1000001)
expect_failure("bad --trials value" seedtable --N 2 --bb-list 1
    --peel-candidates 1 --trials 4294967296)

# Candidate work accounting preserves the raw request while reporting the
# bounded, duplicate-free byte-domain work actually completed.
expect_success(",300,256,256,1,1," seedtable --N 2 --bb-list 1
    --peel-candidates 300 --trials 1)
expect_success(",0,1,1,1,1," seedtable --N 2 --bb-list 1
    --peel-candidates 0 --trials 1)
expect_success("2,1,300,256,256," densetune --N 2 --bb-list 1
    --candidates 300 --trials 1 --loss 0)
expect_success("2,1,0,1,1," densetune --N 2 --bb-list 1
    --candidates 0 --trials 1 --loss 0)

# Every loss-accepting mode rejects values above the shared boundary.
expect_failure("loss must be" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 1 --max-message-mib 1 --loss 0.9900001)
expect_failure("loss must be" precodecheck --N 2 --bb-list 1 --trials 1
    --loss 0.9900001)
expect_failure("loss must be" precodefail --N 64 --bb-list 1 --overhead 0
    --trials 1 --threads 1 --loss 0.9900001)
expect_failure("loss must be" densecheck --N 2 --bb 1 --candidates 1
    --trials 1 --loss 0.9900001)
expect_failure("loss must be" densetune --N 2 --bb-list 1 --candidates 1
    --trials 1 --loss 0.9900001)
expect_failure("loss must be" densecount --N 2 --bb-list 1 --deltas 0
    --trials 1 --loss 0.9900001)
expect_failure("loss must be" densegrid --N 2 --bb-list 1 --deltas 0
    --candidates 1 --trials 1 --loss 0.9900001)

# Parser and interval boundaries.
foreach(loss IN ITEMS -0.1 1.0 nan inf malformed)
    expect_failure("loss" densecheck --N 2 --bb 1 --candidates 1
        --trials 1 --loss ${loss})
    expect_failure("loss" compare --nlo 2 --nhi 2 --trials 1
        --bb-list 1 --max-message-mib 1 --loss ${loss})
    expect_failure("loss" precodecheck --N 2 --bb-list 1 --trials 1
        --loss ${loss})
    expect_failure("loss" precodefail --N 64 --bb-list 1 --overhead 0
        --trials 1 --threads 1 --loss ${loss})
    expect_failure("loss" densetune --N 2 --bb-list 1 --candidates 1
        --trials 1 --loss ${loss})
    expect_failure("loss" densecount --N 2 --bb-list 1 --deltas 0
        --trials 1 --loss ${loss})
    expect_failure("loss" densegrid --N 2 --bb-list 1 --deltas 0
        --candidates 1 --trials 1 --loss ${loss})
endforeach()
foreach(loss IN ITEMS 0 0.37)
    expect_success("loss=${loss}" densecheck --N 2 --bb 1 --candidates 1
        --trials 1 --loss ${loss})
    expect_success("loss=${loss}" compare --nlo 2 --nhi 2 --trials 1
        --bb-list 1 --max-message-mib 1 --loss ${loss})
    expect_success("loss=${loss}" precodecheck --N 2 --bb-list 1
        --trials 1 --loss ${loss})
    expect_success("loss=${loss}" precodefail --N 64 --bb-list 1
        --overhead 0 --trials 1 --threads 1 --loss ${loss})
    expect_success("loss=${loss}" densetune --N 2 --bb-list 1
        --candidates 1 --trials 1 --loss ${loss})
    expect_success("loss=${loss}" densecount --N 2 --bb-list 1
        --deltas 0 --trials 1 --loss ${loss})
    expect_success("loss=${loss}" densegrid --N 2 --bb-list 1
        --deltas 0 --candidates 1 --trials 1 --loss ${loss})
endforeach()
expect_success("loss=0.98999999999999999" densecheck --N 2 --bb 1
    --candidates 1 --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" compare --nlo 2 --nhi 2
    --trials 1 --bb-list 1 --max-message-mib 1 --loss 0.99)
expect_success("loss=0.98999999999999999" precodecheck --N 2 --bb-list 1
    --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" precodefail --N 64 --bb-list 1
    --overhead 0 --trials 1 --threads 1 --loss 0.99)
expect_success("loss=0.98999999999999999" densetune --N 2 --bb-list 1
    --candidates 1 --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" densecount --N 2 --bb-list 1
    --deltas 0 --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" densegrid --N 2 --bb-list 1
    --deltas 0 --candidates 1 --trials 1 --loss 0.99)

# Valid one-trial smoke for all remaining modes.
expect_success("loss boundary oracle: PASS" selftest)
expect_success("mixed null-witness exit policy: PASS" selftest)
expect_success("# compare:" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0)
expect_success("64,8,1,1,0,0,0,0,0," precodecheck --N 64 --bb-list 8
    --trials 1 --loss 0)
expect_success("v2_precode[ ]+8[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode)
expect_success("precode_profile=certified" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode)
expect_success("precode_profile_handoff=encoder-selected-v1" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode)
expect_success("v2_mixed[ ]+8[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed)
expect_success("v2_mixed_cached[ ]+8[ ]+1[ ]+0" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode-cache --precode-profile mixed)
expect_success("mixed_period=64" compare --nlo 64 --nhi 64 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 64)
expect_failure("mixed experiment flags require a mixed precode profile" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile certified --mixed-period 64)
expect_failure("--mixed-period must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 11)
expect_failure("--mixed-period must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 245)
expect_success("mixed_geometry=shared-x" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 64 --mixed-geometry shared-x)
expect_success("mixed_residue_skew=14" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 4 --mixed-period 29
    --mixed-geometry shared-x --mixed-residue-skew 14)
expect_success("mixed_residue_schedule=ramp" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule ramp)
expect_success("mixed_residue_schedule=hashed mixed_residue_hash_seed=0x7 mixed_residue_hash_keyed=1"
    compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-residue-hash-keyed)
expect_success("mixed_independent_extension_residues=1" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed --mixed-mix-count 2
    --mixed-gf16-rows 4 --mixed-period 32 --mixed-geometry shared-x
    --mixed-residue-schedule hashed --mixed-residue-hash-seed 7
    --mixed-residue-hash-keyed --mixed-independent-extension-residues)
expect_success("mixed_residue_buckets_requested=joint-delta" compare
    --nlo 3200 --nhi 3200 --trials 1 --bb-list 4096 --max-message-mib 16
    --loss 0 --precode --precode-profile mixed --mixed-mix-count 2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-independent-extension-residues
    --mixed-residue-buckets joint-delta)
expect_failure("encoder_used=0 decoder_used=0" compare
    --nlo 2 --nhi 2 --trials 1 --bb-list 700000 --max-message-mib 2
    --loss 0 --precode --precode-profile mixed --mixed-mix-count 2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-independent-extension-residues
    --mixed-residue-buckets joint-delta)
expect_failure("explicit --mixed-residue-buckets requires" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed
    --mixed-residue-buckets joint-delta)
expect_failure("independent extension residues require" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed --mixed-gf16-rows 4
    --mixed-period 32 --mixed-geometry shared-x
    --mixed-residue-schedule ramp --mixed-independent-extension-residues)
expect_failure("--mixed-residue-skew must be a corner-preserving" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 4
    --mixed-period 29 --mixed-geometry shared-x --mixed-residue-skew 16)
expect_failure("nonconstant --mixed-residue-schedule requires" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 4
    --mixed-period 28 --mixed-geometry shared-x --mixed-residue-skew 1
    --mixed-residue-schedule ramp)
expect_failure("residue hash seed/keying requires hashed" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 4
    --mixed-period 28 --mixed-geometry shared-x
    --mixed-residue-schedule ramp --mixed-residue-hash-seed 7)
expect_success("mixed_gf16_rows=3.*mixed_geometry=shared-x" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 3
    --mixed-period 64 --mixed-geometry shared-x)
expect_failure("--mixed-gf16-rows must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 1)
expect_failure("--mixed-period must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 3 --mixed-period 12)
expect_failure("unknown --mixed-geometry" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-geometry unknown)
expect_failure("mixed experiment flags require a mixed precode profile" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile certified
    --mixed-geometry shared-x)
expect_success("mixed_mix_count=2" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed --mixed-mix-count 2)
expect_failure("--mixed-mix-count must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-mix-count 1)
expect_success("packet_row_seed_multiplier=0x9e3779b1.*packet_row_seed_avalanche=1"
    compare --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed
    --packet-row-seed-multiplier 2654435761 --packet-row-seed-avalanche)
expect_failure("--packet-row-seed-multiplier must be odd and nonzero"
    compare --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed
    --packet-row-seed-multiplier 2)
expect_success("v2_precode[ ]+17[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 17 --max-message-mib 1 --loss 0 --precode
    --precode-profile certified)
expect_success("v2_cached[ ]+8[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode-cache)
expect_success("encoder_cache=0 decoder_cache=1" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode-decoder-cache)
run_bench(result out err compare --nlo 64 --nhi 64 --trials 2
    --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile both --trial-details)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "precode_profile=both" OR
    NOT out MATCHES "v2_precode[ ]+8[ ]+2[ ]+0" OR
    NOT out MATCHES "v2_mixed[ ]+8[ ]+2[ ]+0" OR
    NOT out MATCHES
        "paired_trial:.*precode_profile=certified.*precode_ok=1" OR
    NOT out MATCHES "paired_trial:.*precode_profile=mixed.*precode_ok=1")
    message(FATAL_ERROR
        "paired certified/mixed precode comparison failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "paired certified/mixed precode comparison")
expect_failure("unknown --precode-profile" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile unknown)
expect_failure("--precode-profile requires" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode-profile mixed)
run_bench(result out err compare --nlo 64 --nhi 64 --trials 1
    --bb-list 17 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
    NOT out STREQUAL "" OR
    NOT err MATCHES "mixed precode profile requires even block bytes")
    message(FATAL_ERROR
        "mixed odd-byte rejection emitted partial output or wrong status\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "mixed odd-byte rejection")
expect_failure("mixed precode profile requires even block bytes" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8,17 --max-message-mib 1
    --loss 0 --precode --precode-profile both)
expect_success("# densetune:" densetune --N 2 --bb-list 8 --candidates 1
    --trials 1 --loss 0)
expect_success("# densecount:" densecount --N 2 --bb-list 8 --deltas 0
    --trials 1 --loss 0)
expect_success("# densegrid:" densegrid --N 2 --bb-list 8 --deltas 0
    --candidates 1 --trials 1 --loss 0)
expect_success("# peelcost:" peelcost --N 2 --bb-list 8 --trials 1
    --structures lt_m1_c16 --precode dense --overhead 0)
expect_success("# precodefail:" precodefail --N 64 --bb-list 8
    --overhead 0,1 --trials 4 --threads 2 --loss 0.1)

function(run_mixed_witness result_var out_var err_var seed trial_count)
    run_bench(result out err precodefail
        --N 945 --bb-list 1280 --overhead 0
        --trials ${trial_count} --threads 2 --loss 0.35 --seed ${seed}
        --schedule burst --completion mixed --mix-count 2
        --binary-dense-rows 13 --mixed-gf256-rows 11
        --mixed-gf16-rows 4 --mixed-period 32 --mixed-geometry shared-x
        --mixed-residue-schedule hashed --mixed-residue-hash-seed 68
        --mixed-residue-hash-keyed --mixed-independent-extension-residues
        --mixed-extension-residue-seed-xor 78 --mixed-null-witnesses)
    set(${result_var} "${result}" PARENT_SCOPE)
    set(${out_var} "${out}" PARENT_SCOPE)
    set(${err_var} "${err}" PARENT_SCOPE)
endfunction()

# The primary trials run on two workers.  Only the deterministic post-join
# replay installs the TLS sink; trial 21 is folded into logical trial zero.
run_mixed_witness(result out err 0xa11ce520f84d877e 2)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "trials=2 threads=2.*mixed_null_witnesses=1" OR
    NOT out MATCHES
        "mixed_null_witness,v=2,N=945,bb=1280,trial=0,status=captured,reason=verified,diagnostic_status=1,replay_stats_ok=1,period=32,schedule=hashed,hash_seed=0xf090c9ff,extension_seed_xor=0x4e,independent_extension=1,L=1019,R=90,binary_rank=75,q=15,quotient_rank=14,d=1.*hash=9351dfef8d448cb7540a5901e8c02974" OR
    NOT out MATCHES
        "mixed_null_row,v=1,row=0,nz=2,source=2,precode=0,source_peeled=1,source_inactive=1,precode_peeled=0,precode_inactive=0,gf256=2,gf16=0,sf_occ=1,sf_cancel=1,sf_cancel_terms=2,sf_max=2,sf_hash=4f7f2ed9c09b8250,ex_occ=1,ex_cancel=1,ex_cancel_terms=2,ex_max=2,ex_hash=beb11b1eb16788b2,pair_occ=1,pair_cancel=1,pair_cancel_terms=2,pair_max=2,pair_hash=cf3f4a67f907b9cf,sf_top=11:2:0000,ex_top=13:2:0000")
    message(FATAL_ERROR
        "mixed null-witness post-run replay failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "mixed null-witness post-run replay")

# A second real solve pins canonical full-L ordering and hashes for d=2.
run_bench(result out err precodefail --N 64 --bb-list 8 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --seed 3 --schedule adversarial
    --completion mixed --source-hits 2 --mix-count 2 --mixed-period 12
    --mixed-null-witnesses)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "64,8,periodic,2,0,1,0,1,0," OR
    NOT out MATCHES
        "mixed_null_witness,v=2,N=64,bb=8,trial=0,status=captured,reason=verified,diagnostic_status=1,replay_stats_ok=1.*L=107,R=33,binary_rank=21,q=12,quotient_rank=10,d=2.*exact_words=214,exact_size_ok=1,hash=16daf34d0ed35f45b4ac5bf49004e1a1" OR
    NOT out MATCHES
        "mixed_null_row,v=1,row=0.*sf_hash=cb3ab68d3ac8b9c4.*ex_hash=be52b9be6f9d47a6.*pair_hash=78e77a1db4dacc7b" OR
    NOT out MATCHES
        "mixed_null_row,v=1,row=1.*sf_hash=596023278d8fd107.*ex_hash=9dd32dc060d03dc5.*pair_hash=86fb1d0bc060b438" OR
    out MATCHES "mixed_null_row,v=1,row=2")
    message(FATAL_ERROR
        "mixed null-witness d=2 classification failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "mixed null-witness d=2 classification")

run_mixed_witness(result out err 0x9d2c4a71 2)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "trials=2 threads=2.*mixed_null_witnesses=1" OR
    NOT out MATCHES "945,1280,periodic,2,0,2,2,0,0," OR
    NOT out MATCHES
        "mixed_null_witness,v=2,N=945,bb=1280,trial=-1,status=none,reason=no_need_more,diagnostic_status=0,replay_stats_ok=0.*L=1019,R=0,binary_rank=0,q=0,quotient_rank=0,d=0.*exact_words=0,exact_size_ok=0,hash=00000000000000000000000000000000" OR
    out MATCHES "mixed_null_row")
    message(FATAL_ERROR
        "mixed null-witness benign record failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "mixed null-witness benign record")

run_bench(result out err precodefail --N 64 --bb-list 8 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --seed 1 --schedule adversarial
    --completion mixed --source-hits 2 --mix-count 1
    --mixed-null-witnesses)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "64,8,periodic,1,0,1,0,1,0," OR
    NOT out MATCHES
        "mixed_null_witness,v=2,N=64,bb=8,trial=0,status=skipped,reason=ineligible_need_more,diagnostic_status=0,replay_stats_ok=0.*L=107,R=32,binary_rank=19,q=13,quotient_rank=0,d=13.*exact_words=0,exact_size_ok=0,hash=00000000000000000000000000000000" OR
    out MATCHES "mixed_null_row")
    message(FATAL_ERROR
        "mixed null-witness skipped record failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "mixed null-witness skipped record")
expect_failure("mixed experiment flags require --completion mixed"
    precodefail --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1
    --loss 0.1 --completion certified --mixed-null-witnesses)

string(CONCAT precodefail_paired_pattern
    "precodefail_paired:.*completion=mixed.*mix2_fail=.*mix3_fail=.*"
    "both_fail=.*mix2_only=.*mix3_only=.*both_success=.*mcnemar_p=.*"
    "mix2_seed_attempt=.*mix3_seed_attempt=.*mix2_wilson95=.*"
    "mix3_wilson95=")
expect_success("${precodefail_paired_pattern}" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 4 --threads 2 --loss 0.1
    --mix-count 2,3 --completion mixed --payload-e2e)
expect_success("mixed_period=64 mixed_gf256_rows=10 mixed_gf16_rows=2 mixed_geometry=shared-x.*full_payload_solve=1"
    precodefail
    --N 64 --bb-list 1280 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-period 64 --mixed-geometry shared-x
    --full-payload-solve)
expect_success("mixed_period=64 mixed_gf256_rows=10 mixed_gf16_rows=3 mixed_geometry=shared-x"
    precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 3 --mixed-period 64
    --mixed-geometry shared-x)
expect_success("mixed_period=32 mixed_gf256_rows=10 mixed_gf16_rows=4 mixed_geometry=shared-x"
    precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x)
expect_success("source_hits_override=2" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 2 --threads 2 --loss 0.1 --source-hits 2)
expect_success("packet_peel_seed_xor=0x7" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 2 --threads 2 --loss 0.1
    --packet-peel-seed-xor 7)
expect_success("odd_packet_peel_seed_xor=0x13" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --odd-packet-peel-seed-xor 19 --payload-e2e)
expect_success("packet_row_seed_multiplier=0x9e3779b1 packet_row_seed_avalanche=1"
    precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --packet-row-seed-multiplier 2654435761
    --packet-row-seed-avalanche --payload-e2e)
expect_success("overhead_stream=paired" precodefail
    --N 64 --bb-list 8 --overhead 0,1 --trials 2 --threads 2 --loss 0.1
    --paired-overhead-stream)
expect_failure("--packet-row-seed-multiplier must be odd and nonzero"
    precodefail --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1
    --loss 0.1 --packet-row-seed-multiplier 2)
expect_success("seed_block_bytes_override=2" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --seed-block-bytes 2 --payload-e2e)
expect_failure("--seed-block-bytes must be nonzero" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --seed-block-bytes 0)
expect_success("mixed_residue_skew=14" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 29
    --mixed-geometry shared-x --mixed-residue-skew 14)
expect_success("mixed_residue_schedule=ramp" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule ramp)
expect_success("mixed_residue_schedule=hashed mixed_residue_hash_seed=0x7 mixed_residue_hash_keyed=1"
    precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-residue-hash-keyed)
expect_success("mixed_independent_extension_residues=1" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mix-count 2 --payload-e2e
    --mixed-gf16-rows 4 --mixed-period 32 --mixed-geometry shared-x
    --mixed-residue-schedule hashed --mixed-residue-hash-seed 7
    --mixed-residue-hash-keyed --mixed-independent-extension-residues)
expect_success("mixed_residue_buckets_requested=joint-delta" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mix-count 2 --full-payload-solve
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-independent-extension-residues
    --mixed-residue-buckets joint-delta)
expect_success("2,2,periodic,2,0,2,1,1,0" precodefail
    --N 2 --bb-list 2 --overhead 0 --trials 2 --threads 1 --loss 0.35
    --completion mixed --mix-count 2 --schedule adversarial
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues
    --mixed-residue-buckets joint-delta)
expect_failure("trial 1 fell back" precodefail
    --N 2 --bb-list 700000 --overhead 0 --trials 2 --threads 1 --loss 0.35
    --completion mixed --mix-count 2 --full-payload-solve
    --schedule adversarial --seed-block-bytes 2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues
    --mixed-residue-buckets joint-delta)
expect_success("mixed_gf256_rows=11" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mix-count 2 --mixed-gf256-rows 11
    --mixed-gf16-rows 4 --mixed-period 32 --mixed-geometry shared-x
    --mixed-residue-schedule hashed --mixed-residue-hash-seed 7
    --mixed-residue-hash-keyed --mixed-independent-extension-residues
    --payload-e2e)
expect_success("4,64.*0x39" precodefail
    --N 4 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 2 --threads 2 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v1
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("42,64.*0x0" precodefail
    --N 42 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v1
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("1683,64.*0x0" precodefail
    --N 1683 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v1
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("4,64.*0x39" precodefail
    --N 4 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("1683,64.*0x13" precodefail
    --N 1683 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("15182,64.*0x62" precodefail
    --N 15182 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("21394,64.*0x1a" precodefail
    --N 21394 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("24432,64.*0x4b" precodefail
    --N 24432 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("34207,64.*0xd5" precodefail
    --N 34207 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("62039,64.*0x2" precodefail
    --N 62039 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("4,64.*0x39" precodefail
    --N 4 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("62039,64.*0x2" precodefail
    --N 62039 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("10,64.*0x8b" precodefail
    --N 10 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("20,64.*0x8c" precodefail
    --N 20 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("11414,64.*0x56" precodefail
    --N 11414 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("48567,64.*0xd1" precodefail
    --N 48567 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("49312,64.*0x34" precodefail
    --N 49312 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("49842,64.*0xbc" precodefail
    --N 49842 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("50281,64.*0x79" precodefail
    --N 50281 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("51375,64.*0xc0" precodefail
    --N 51375 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("53503,64.*0xee" precodefail
    --N 53503 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --schedule adversarial
    --completion mixed --mix-count 2 --packet-peel-seed-table normalized-h15-v3
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)

function(expect_normalized_h15_seed table K seed_hex)
    expect_success("${K},64.*0x${seed_hex}[^0-9a-fA-F]" precodefail
        --N ${K} --bb-list 64 --seed-block-bytes 1280 --overhead 0
        --trials 1 --threads 1 --loss 0.35 --schedule adversarial
        --completion mixed --mix-count 2 --packet-peel-seed-table ${table}
        --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
        --mixed-geometry shared-x --mixed-residue-schedule hashed
        --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
        --mixed-independent-extension-residues)
endfunction()

# V4 preserves every earlier mapping and adds exactly the discovery-selected
# rank-one salts that the frozen independent holdout validated.
set(normalized_h15_v4_seed_cases
    "16,6"
    "39559,3c"
    "40831,b3"
    "43742,1b"
    "43751,63"
    "45168,6c"
    "45464,22"
    "45857,31"
    "45903,3a"
    "46296,4"
    "46606,eb"
    "46933,6a"
    "47029,75"
    "47105,51"
    "47307,b2"
    "48231,e1"
    "48311,7a"
    "48466,57"
    "49124,ed"
    "49412,ad"
    "49486,8e"
    "49627,ac"
    "49727,8f"
    "49865,ff"
    "50689,8e"
    "50885,3f"
    "50899,c"
    "51494,d0"
    "52935,8"
    "53613,1e"
    "53697,cc"
    "53804,a9")
foreach(seed_case IN LISTS normalized_h15_v4_seed_cases)
    string(REPLACE "," ";" seed_fields "${seed_case}")
    list(GET seed_fields 0 seed_K)
    list(GET seed_fields 1 seed_hex)
    expect_normalized_h15_seed(normalized-h15-v4 ${seed_K} ${seed_hex})
endforeach()
expect_normalized_h15_seed(normalized-h15-v4 4 39)
expect_normalized_h15_seed(normalized-h15-v4 62039 2)
expect_normalized_h15_seed(normalized-h15-v4 10 8b)
expect_normalized_h15_seed(normalized-h15-v4 42 0)

# V5 is a benchmark-only mapping: one selected salt proves precodefail wiring,
# while inherited nonzero and retained-zero cases prove complete-table parity.
expect_normalized_h15_seed(normalized-h15-v5 955 27)
expect_normalized_h15_seed(normalized-h15-v5 4 39)
expect_normalized_h15_seed(normalized-h15-v5 42 0)

expect_failure("conflicts with --packet-peel-seed-xor" precodefail
    --N 4 --bb-list 64 --seed-block-bytes 1280 --overhead 0
    --trials 1 --threads 1 --loss 0.35 --completion mixed --mix-count 2
    --packet-peel-seed-table normalized-h15-v1 --packet-peel-seed-xor 57
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_failure("requires its normalized H15/mix2 geometry" precodefail
    --N 4 --bb-list 64 --overhead 0 --trials 1 --threads 1 --loss 0.35
    --completion mixed --mix-count 2
    --packet-peel-seed-table normalized-h15-v1
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_failure("unknown --packet-peel-seed-table" precodefail
    --N 4 --bb-list 64 --overhead 0 --trials 1 --threads 1 --loss 0.35
    --packet-peel-seed-table unknown)
expect_success("mixed_gf256_rows=12" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mix-count 2 --mixed-gf256-rows 12
    --mixed-gf16-rows 4 --mixed-period 32 --mixed-geometry shared-x
    --mixed-residue-schedule hashed --mixed-residue-hash-seed 7
    --mixed-residue-hash-keyed --mixed-independent-extension-residues
    --payload-e2e)
expect_success("mixed_period=31 mixed_gf256_rows=12" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mix-count 2 --mixed-gf256-rows 12
    --mixed-gf16-rows 4 --mixed-period 31 --mixed-geometry shared-x
    --mixed-residue-schedule hashed --mixed-residue-hash-seed 7
    --mixed-residue-hash-keyed --mixed-independent-extension-residues
    --payload-e2e)
expect_failure("--mixed-gf256-rows must be in" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf256-rows 9)
expect_failure("--mixed-gf256-rows must be in" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf256-rows 13)
expect_failure("validated 12" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf256-rows 12 --mixed-gf16-rows 4
    --mixed-period 30 --mixed-geometry shared-x)
expect_failure("validated 12" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf256-rows 12 --mixed-gf16-rows 2
    --mixed-period 32 --mixed-geometry shared-x)
expect_failure("use shared-x for an extra row" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf256-rows 11 --mixed-geometry frozen)
expect_success("mixed_gf256_rows=11" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0.1 --precode --precode-profile mixed
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("mixed_gf256_rows=12" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0.1 --precode --precode-profile mixed
    --mixed-gf256-rows 12 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues)
expect_success("mixed_extension_residue_seed_xor=0x17" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mix-count 2 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues
    --mixed-extension-residue-seed-xor 23)
expect_failure("--mixed-extension-residue-seed-xor requires" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-residue-hash-keyed
    --mixed-extension-residue-seed-xor 23)
expect_success("binary_dense_rows_override=16" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-rows 16)
expect_failure("--binary-dense-rows must be in" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --binary-dense-rows 0)
expect_failure("--binary-dense-rows must be in" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --binary-dense-rows 65)
expect_success("gf256_heavy_rows_override=16" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion certified --gf256-heavy-rows 16)
expect_failure("--gf256-heavy-rows must be in" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --gf256-heavy-rows 0)
expect_failure("--gf256-heavy-rows must be in" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --gf256-heavy-rows 129)
expect_failure("--gf256-heavy-rows requires" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --gf256-heavy-rows 16)
expect_failure("independent extension residues require" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule ramp
    --mixed-independent-extension-residues)
expect_failure("--mixed-gf16-rows must be in" precodefail --N 64
    --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf16-rows 5)
expect_failure("--source-hits must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --source-hits 0)
expect_failure("--mixed-period must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-gf16-rows 3 --mixed-period 12)
expect_failure("mixed experiment flags require --completion mixed" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --mixed-period 64)
expect_failure("--mixed-period must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-period 11)
expect_failure("--mixed-period must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-period 245)
expect_failure("unknown --mixed-geometry" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-geometry unknown)
expect_failure("mixed experiment flags require --completion mixed" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --mixed-geometry shared-x)
expect_failure("unknown --completion" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion unknown)
expect_failure("mixed completion requires even block bytes" precodefail --N 64
    --bb-list 7 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed)
expect_failure("mixed completion requires periodic heavy family" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --heavy-family hashed)
expect_failure("working set exceeds|message dimensions are unsupported" precodefail --N 64000
    --bb-list 2147483647 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --payload-e2e)
expect_failure("working set exceeds|message dimensions are unsupported" precodefail --N 64000
    --bb-list 2200000 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --payload-e2e)

# Common packet schedules are accepted by both E2E comparison entry points
# and the rank-only precode failure harness.
foreach(schedule IN ITEMS iid burst permutation systematic-first repair-only
        adversarial)
    expect_success("schedule=${schedule}" compare --nlo 64 --nhi 64
        --trials 1 --bb-list 8 --max-message-mib 1 --loss 0.1
        --schedule ${schedule} --precode)
    expect_success("schedule=${schedule}" precodecheck --N 64 --bb-list 8
        --trials 1 --loss 0.1 --schedule ${schedule})
    expect_success("schedule=${schedule}" precodefail --N 64 --bb-list 8
        --overhead 0 --trials 1 --threads 1 --loss 0.1
        --schedule ${schedule})
endforeach()
expect_success("paired_trial: schedule=burst" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0.1
    --schedule burst --precode --trial-details)
expect_success("cached_ok=1.*cached_oh=[0-9]+.*cached_delta=" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0.1 --schedule iid --precode-cache --trial-details)
expect_success("precode_trial: schedule=permutation" precodecheck --N 64
    --bb-list 8 --trials 1 --loss 0.1 --schedule permutation
    --trial-details)
# High burst loss used to clamp its start probability to one and drop forever.
# Exercise both consumers at the accepted loss boundary.
expect_success("schedule=burst" compare --nlo 8 --nhi 8 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0.99 --schedule burst --precode)
expect_success("schedule=burst" precodecheck --N 64 --bb-list 8
    --trials 1 --loss 0.99 --schedule burst)
expect_success("schedule=burst" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.99 --schedule burst)
expect_failure("unknown compare schedule" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --schedule unknown)
expect_failure("unknown precodecheck schedule" precodecheck --N 64
    --bb-list 8 --trials 1 --schedule unknown)
expect_failure("unknown precodefail schedule" precodefail --N 64
    --bb-list 8 --overhead 0 --trials 1 --threads 1 --schedule unknown)
expect_success("hashed" precodefail --N 64 --bb-list 1
    --overhead 0 --trials 1 --threads 1 --loss 0.1
    --heavy-family periodic,hashed)
expect_failure("unknown --heavy-family" precodefail --N 64 --bb-list 1
    --overhead 0 --trials 1 --threads 1 --loss 0.1
    --heavy-family unknown)

# Thread parsing and the partially-launched worker cleanup path.
expect_failure("bad --threads value" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads malformed --loss 0.1)
expect_failure("threads must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 0 --loss 0.1)
expect_failure("threads must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 257 --loss 0.1)
expect_failure("bad --threads value" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 4294967296 --loss 0.1)
expect_failure("thread launch failed" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 4 --threads 2 --loss 0.1
    --fail-thread-launch-after 1)

# Invalid experiment grids must fail before emitting a result header.
run_bench(result out err peelcost --N 2 --bb-list 8 --trials 1
    --structures lt_m1_c16 --precode dense --overhead -1)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
    NOT out STREQUAL "" OR
    NOT err MATCHES "overhead must be non-negative")
    message(FATAL_ERROR
        "peelcost invalid-grid failure emitted partial output\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "invalid peelcost grid")

# Strict message caps and widened size checks.
expect_failure("cap cannot accommodate" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 524289 --max-message-mib 1 --loss 0)
expect_success("# compare:" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 524288 --max-message-mib 1 --loss 0)
expect_failure("working set" densecheck --N 64000 --bb 2147483647
    --candidates 1 --trials 1 --loss 0)

if(UNIX AND NOT SKIP_CONSTRAINED)
    set(constrained_commands
        "densecheck --N 2 --bb 100000000 --candidates 1 --trials 1 --loss 0"
        "densetune --N 2 --bb-list 100000000 --candidates 1 --trials 1 --loss 0"
        "densecount --N 2 --bb-list 100000000 --deltas 0 --trials 1 --loss 0"
        "densegrid --N 2 --bb-list 100000000 --deltas 0 --candidates 1 --trials 1 --loss 0")
    foreach(command IN LISTS constrained_commands)
        execute_process(
            COMMAND bash -c "ulimit -v 65536; exec '${BENCH}' ${command}"
            RESULT_VARIABLE result
            OUTPUT_VARIABLE out
            ERROR_VARIABLE err
            TIMEOUT 30)
        if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
            NOT out STREQUAL "" OR
            NOT err MATCHES "working set cannot be allocated")
            message(FATAL_ERROR
                "constrained-memory failure was not clean: ${command}\n"
                "result=${result}\nstdout=${out}\nstderr=${err}")
        endif()
        reject_sanitizer("${out}${err}" "constrained-memory failure: ${command}")
    endforeach()

    execute_process(
        COMMAND bash -c
            "ulimit -v 65536; exec '${BENCH}' peelcost --N 2 --bb-list 8 --trials 4294967295 --structures lt_m1_c16 --precode dense --overhead 0"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 30)
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
        NOT out STREQUAL "" OR
        NOT err MATCHES "evaluation storage cannot be allocated")
        message(FATAL_ERROR
            "constrained peelcost failure was not clean\n"
            "result=${result}\nstdout=${out}\nstderr=${err}")
    endif()
    reject_sanitizer("${out}${err}" "constrained peelcost failure")

    execute_process(
        COMMAND bash -c
            "ulimit -v 65536; exec '${BENCH}' compare --nlo 2 --nhi 2 --trials 1 --bb-list 100000000 --max-message-mib 64 --loss 0"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 30)
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
        NOT out STREQUAL "" OR
        NOT err MATCHES "cap cannot accommodate")
        message(FATAL_ERROR
            "constrained compare failure was not clean\n"
            "result=${result}\nstdout=${out}\nstderr=${err}")
    endif()
    reject_sanitizer("${out}${err}" "constrained compare failure")
endif()
