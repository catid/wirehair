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
expect_success(
    "mixed_period=64.*mixed_grouped_gf256_rows=0 mixed_grouped_gf256_hash_seed=0x0 mixed_grouped_final_h_a_columns=0"
    compare --nlo 64 --nhi 64 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 64)
foreach(grouped_period IN ITEMS 32 48 64 96)
    expect_success(
        "mixed_grouped_gf256_rows=3 mixed_grouped_gf256_hash_seed=0x[1-9a-fA-F][0-9a-fA-F]* mixed_grouped_final_h_a_columns=12"
        compare --nlo 64 --nhi 64 --trials 1 --bb-list 8
        --max-message-mib 1 --loss 0 --precode --precode-profile mixed
        --mixed-period ${grouped_period} --mixed-geometry shared-x
        --mixed-grouped-gf256-rows 3)
endforeach()
expect_success(
    "mixed_grouped_gf256_rows=1.*mixed_residue_buckets_requested=joint-delta"
    compare --nlo 64 --nhi 64 --trials 1 --bb-list 8
    --max-message-mib 1 --loss 0 --precode --precode-profile mixed
    --mixed-period 32 --mixed-geometry shared-x
    --mixed-grouped-gf256-rows 1 --mixed-residue-buckets joint-delta)
expect_failure("mixed experiment flags require a mixed precode profile" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile certified
    --mixed-grouped-gf256-rows 0)
expect_failure("--mixed-grouped-gf256-rows must be in" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed --mixed-period 32
    --mixed-geometry shared-x --mixed-grouped-gf256-rows 10)
expect_failure("nonzero grouping requires shared-x" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed --mixed-period 32
    --mixed-grouped-gf256-rows 3)
expect_failure("nonzero grouping requires shared-x constant-A" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-grouped-gf256-rows 3)
expect_failure("--mixed-grouped-gf256-rows requires a value" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile mixed
    --mixed-grouped-gf256-rows)
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

# The same folded loss stream is a true D12 sparse-alias fixture.  Keep this
# separate from the historical D13 diagnostic above: the two-anchor hook must
# repair the default-D12 witness without changing the checked-in D13 golden.
set(d12_witness_args
    precodefail --N 945 --bb-list 1280 --overhead 0 --trials 2 --threads 2
    --loss 0.35 --seed 0xa11ce520f84d877e --schedule burst
    --completion mixed --mix-count 2 --mixed-gf256-rows 11
    --mixed-gf16-rows 4 --mixed-period 32 --mixed-geometry shared-x
    --mixed-residue-schedule hashed --mixed-residue-hash-seed 68
    --mixed-residue-hash-keyed --mixed-independent-extension-residues
    --mixed-extension-residue-seed-xor 78 --mixed-null-witnesses)
run_bench(result out err ${d12_witness_args})
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "945,1280,periodic,2,0,2,1,1,0," OR
    NOT out MATCHES
        "N=945,bb=1280,trial=0,status=captured,reason=verified.*L=1018,R=86,binary_rank=71,q=15,quotient_rank=14,d=1.*hash=fbbd73391d7826a04e44d80115698ea5")
    message(FATAL_ERROR
        "default-D12 sparse witness replay failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "default-D12 sparse witness replay")

run_bench(result out err ${d12_witness_args} --binary-dense-two-anchor)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "binary_dense_two_anchor=1" OR
    NOT out MATCHES "945,1280,periodic,2,0,2,2,0,0," OR
    NOT out MATCHES
        "N=945,bb=1280,trial=-1,status=none,reason=no_need_more")
    message(FATAL_ERROR
        "two-anchor D12 witness repair failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "two-anchor D12 witness repair")

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
run_bench(result out err precodefail --N 64 --bb-list 8 --overhead 12
    --trials 2 --threads 2 --loss 0.1 --completion mixed --mix-count 2
    --payload-e2e --mixed-null-witnesses --mixed-period 32
    --mixed-geometry shared-x --mixed-grouped-gf256-rows 9
    --mixed-residue-buckets dual)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES
        "mixed_grouped_gf256_rows=9 mixed_grouped_gf256_hash_seed=0x[1-9a-fA-F][0-9a-fA-F]* mixed_grouped_final_h_a_columns=12.*mixed_residue_buckets_requested=dual" OR
    NOT out MATCHES
        "mixed_null_witness,v=2.*grouped_gf256_rows=9,grouped_first_row=1,grouped_hash_seed=0x[1-9a-fA-F][0-9a-fA-F]*,grouped_final_h_a_columns=12")
    message(FATAL_ERROR
        "grouped GF256 precodefail/null receipt failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "grouped GF256 precodefail/null receipt")
expect_failure("nonzero grouping requires shared-x constant-A" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-period 32 --mixed-geometry shared-x
    --mixed-residue-schedule hashed --mixed-residue-hash-seed 7
    --mixed-independent-extension-residues --mixed-grouped-gf256-rows 3)
expect_failure("mixed experiment flags require --completion mixed" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion certified --mixed-grouped-gf256-rows 0)
expect_failure("--mixed-grouped-gf256-rows requires a value" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-grouped-gf256-rows)
expect_failure("--mixed-extension-residue-seed-xor requires" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7 --mixed-residue-hash-keyed
    --mixed-extension-residue-seed-xor 23)
expect_success("binary_dense_rows_override=16" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-rows 16)
expect_success(
    "binary_dense_two_anchor=1.*binary_dense_two_anchor_phase=0" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-two-anchor)
expect_success(
    "binary_dense_two_anchor=1.*binary_dense_two_anchor_phase=1" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-two-anchor
    --binary-dense-two-anchor-phase 1)
expect_success(
    "binary_dense_two_anchor=1.*binary_dense_two_anchor_phase=2" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-two-anchor
    --binary-dense-two-anchor-phase 2)
expect_failure("--binary-dense-two-anchor-phase requires" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-two-anchor-phase 1)
expect_failure("--binary-dense-two-anchor-phase must be in" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-two-anchor
    --binary-dense-two-anchor-phase 3)
expect_success(
    "packet_peel_seed_table=normalized-h15-v4.*binary_dense_two_anchor=1"
    precodefail --N 64 --bb-list 64 --overhead 0 --trials 1 --threads 1
    --loss 0.35 --schedule burst --completion mixed --mix-count 2
    --mixed-gf256-rows 11 --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68 --mixed-residue-hash-keyed
    --mixed-independent-extension-residues
    --mixed-extension-residue-seed-xor 78 --seed-block-bytes 1280
    --packet-peel-seed-table normalized-h15-v4
    --binary-dense-two-anchor)

# The post-selection preferred-attempt harness is benchmark-only.  Route mode
# emits canonical bytes; cached candidate/paired modes hash and parse those
# bytes themselves before they skip systematic probes.
set(route_context_sha256
    "89abcdef0123456789abcdef0123456789abcdef0123456789abcdef01234567")
run_bench(route_result route_k3 route_err preferredattempt --mode route
    --N 3,4 --bb-list 64 --preferred-map "3@64=0,1,2|4@64=none"
    --route-context-sha256 "${route_context_sha256}")
if(NOT route_result EQUAL 0 OR NOT "${route_err}" STREQUAL "")
    message(FATAL_ERROR "cannot create K3 route fixture: ${route_err}")
endif()
set(route_k3_path "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k3.csv")
file(WRITE "${route_k3_path}" "${route_k3}")
file(SHA256 "${route_k3_path}" route_k3_sha256)
if(NOT "${route_k3}" MATCHES
   "preferredattempt-route-manifest: schema=v1.*context_sha256=${route_context_sha256}" OR
   NOT "${route_k3}" MATCHES "3,64,preferred,0,1,1,0,1,0,0,2,0" OR
   NOT "${route_k3}" MATCHES "3,64,preferred,1,1,1,1,0,1,0,0,0" OR
   NOT "${route_k3}" MATCHES "3,64,preferred,2,1,1,0,1,0,0,0,1" OR
   NOT "${route_k3}" MATCHES "4,64,control,-1,[0-9]+,[0-9]+,1,0,1,0,[0-9]+,0")
    message(FATAL_ERROR "unexpected K3 route fixture:\n${route_k3}")
endif()

run_bench(route_result route_k3_only route_err preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "3@64=0,1,2"
    --route-context-sha256 "${route_context_sha256}")
if(NOT route_result EQUAL 0 OR NOT "${route_err}" STREQUAL "")
    message(FATAL_ERROR "cannot create K3-only route fixture: ${route_err}")
endif()
set(route_k3_only_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k3-only.csv")
file(WRITE "${route_k3_only_path}" "${route_k3_only}")
file(SHA256 "${route_k3_only_path}" route_k3_only_sha256)

run_bench(route_result route_k4096 route_err preferredattempt --mode route
    --N 4096 --bb-list 64 --preferred-map "4096@64=0,1"
    --route-context-sha256 "${route_context_sha256}")
if(NOT route_result EQUAL 0 OR NOT "${route_err}" STREQUAL "")
    message(FATAL_ERROR "cannot create K4096 route fixture: ${route_err}")
endif()
set(route_k4096_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k4096.csv")
file(WRITE "${route_k4096_path}" "${route_k4096}")
file(SHA256 "${route_k4096_path}" route_k4096_sha256)

# Paired and timing consumers bind canonical selected manifests: one physical
# route row per preferred K/width.  Candidate batches intentionally retain the
# multi-attempt artifacts above.
run_bench(route_k3_p0_result route_k3_p0 route_k3_p0_err
    preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "3@64=0"
    --route-context-sha256 "${route_context_sha256}")
set(route_k3_p0_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k3-p0.csv")
file(WRITE "${route_k3_p0_path}" "${route_k3_p0}")
file(SHA256 "${route_k3_p0_path}" route_k3_p0_sha256)
run_bench(route_k3_p1_result route_k3_p1 route_k3_p1_err
    preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "3@64=1"
    --route-context-sha256 "${route_context_sha256}")
set(route_k3_p1_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k3-p1.csv")
file(WRITE "${route_k3_p1_path}" "${route_k3_p1}")
file(SHA256 "${route_k3_p1_path}" route_k3_p1_sha256)
run_bench(route_k3_p2_result route_k3_p2 route_k3_p2_err
    preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "3@64=2"
    --route-context-sha256 "${route_context_sha256}")
set(route_k3_p2_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k3-p2.csv")
file(WRITE "${route_k3_p2_path}" "${route_k3_p2}")
file(SHA256 "${route_k3_p2_path}" route_k3_p2_sha256)
run_bench(route_k3_k4_p0_result route_k3_k4_p0 route_k3_k4_p0_err
    preferredattempt --mode route
    --N 3,4 --bb-list 64 --preferred-map "3@64=0|4@64=none"
    --route-context-sha256 "${route_context_sha256}")
set(route_k3_k4_p0_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k3-k4-p0.csv")
file(WRITE "${route_k3_k4_p0_path}" "${route_k3_k4_p0}")
file(SHA256 "${route_k3_k4_p0_path}" route_k3_k4_p0_sha256)
run_bench(route_k4096_p0_result route_k4096_p0 route_k4096_p0_err
    preferredattempt --mode route
    --N 4096 --bb-list 64 --preferred-map "4096@64=0"
    --route-context-sha256 "${route_context_sha256}")
set(route_k4096_p0_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k4096-p0.csv")
file(WRITE "${route_k4096_p0_path}" "${route_k4096_p0}")
file(SHA256 "${route_k4096_p0_path}" route_k4096_p0_sha256)
run_bench(route_k4096_p1_result route_k4096_p1 route_k4096_p1_err
    preferredattempt --mode route
    --N 4096 --bb-list 64 --preferred-map "4096@64=1"
    --route-context-sha256 "${route_context_sha256}")
if(NOT route_k3_p0_result EQUAL 0 OR route_k3_p0_err OR
   NOT route_k3_p1_result EQUAL 0 OR route_k3_p1_err OR
   NOT route_k3_p2_result EQUAL 0 OR route_k3_p2_err OR
   NOT route_k3_k4_p0_result EQUAL 0 OR route_k3_k4_p0_err OR
   NOT route_k4096_p0_result EQUAL 0 OR route_k4096_p0_err OR
   NOT route_k4096_p1_result EQUAL 0 OR route_k4096_p1_err)
    message(FATAL_ERROR
        "cannot create selected route fixtures\n"
        "k3p0=${route_k3_p0_err}\nk3p1=${route_k3_p1_err}\n"
        "k3p2=${route_k3_p2_err}\nk3k4=${route_k3_k4_p0_err}\n"
        "k4096p0=${route_k4096_p0_err}\n"
        "k4096p1=${route_k4096_p1_err}")
endif()
set(route_k4096_p1_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k4096-p1.csv")
file(WRITE "${route_k4096_p1_path}" "${route_k4096_p1}")
file(SHA256 "${route_k4096_p1_path}" route_k4096_p1_sha256)

# A selected attempt may equal the canonical attempt at an encountered wider
# width.  The routing protocol calls this a valid neutral alias (while bb64
# must remain direct), so preferredtiming must time two identical actual arms
# instead of rejecting the frozen K or substituting a sample.
run_bench(route_result route_k4096_wide_noop route_err preferredattempt
    --mode route --N 4096 --bb-list 1280
    --preferred-map "4096@1280=0"
    --route-context-sha256 "${route_context_sha256}")
if(NOT route_result EQUAL 0 OR NOT "${route_err}" STREQUAL "" OR
   NOT "${route_k4096_wide_noop}" MATCHES
       "4096,1280,preferred,0,0,0,1,0,1,0,1,0")
    message(FATAL_ERROR
        "cannot create K4096 wide no-op route fixture: ${route_err}\n"
        "${route_k4096_wide_noop}")
endif()
set(route_k4096_wide_noop_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-k4096-wide-noop.csv")
file(WRITE "${route_k4096_wide_noop_path}" "${route_k4096_wide_noop}")
file(SHA256 "${route_k4096_wide_noop_path}"
    route_k4096_wide_noop_sha256)

# Full-width solve timing uses one immutable K+4 packet trace for both arms,
# prebuilds systems outside the timer, and emits four exact CTTCTCCT cycles.
# The first cycle remains in the audit record but is excluded by the external
# exact-rational analyzer.  A small eviction buffer keeps this CLI fixture
# bounded; the frozen campaign enforces max(2*LLC,256MiB).
run_bench(timing_solve_result timing_solve timing_solve_err preferredtiming
    --N 4096 --bb 64 --preferred-attempt 1 --evict-bytes 4096
    --metric solve --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REGEX MATCHALL "\n4096,64,solve," timing_solve_rows
    "${timing_solve}")
list(LENGTH timing_solve_rows timing_solve_row_count)
string(REGEX MATCHALL
    "\n4096,64,solve,[0-3],[0-7],(control|candidate),[0-9]+,0,[0-9]+,0,"
    timing_solve_valid_rows "${timing_solve}")
list(LENGTH timing_solve_valid_rows timing_solve_valid_row_count)
if(NOT timing_solve_result EQUAL 0 OR timing_solve_err OR
   NOT timing_solve_row_count EQUAL 32 OR
   NOT timing_solve_valid_row_count EQUAL 32 OR NOT timing_solve MATCHES
       "metric=solve.*cycles=4 order=CTTCTCCT discard_cycle=0.*cycle_mode=full.*cycle_index=all.*overhead=4.*payload=distinct-zero-v1.*payload_alignment=64.*payload_prefaulted=1" OR
   NOT timing_solve MATCHES ",262144,268160")
    message(FATAL_ERROR
        "preferred solve timing fixture failed\n${timing_solve}\n${timing_solve_err}")
endif()
foreach(cycle RANGE 0 3)
    foreach(slot RANGE 0 7)
        if(slot EQUAL 0 OR slot EQUAL 3 OR slot EQUAL 5 OR slot EQUAL 6)
            set(timing_arm "control")
        else()
            set(timing_arm "candidate")
        endif()
        if(NOT timing_solve MATCHES
           "4096,64,solve,${cycle},${slot},${timing_arm},")
            message(FATAL_ERROR
                "preferred solve timing order mismatch cycle=${cycle} slot=${slot}")
        endif()
    endforeach()
endforeach()

run_bench(timing_retry_result timing_retry timing_retry_err preferredtiming
    --N 4096 --bb 64 --preferred-attempt 1 --evict-bytes 4096
    --metric solve --cycle-index 2
    --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REGEX MATCHALL "\n4096,64,solve,2," timing_retry_rows
    "${timing_retry}")
list(LENGTH timing_retry_rows timing_retry_row_count)
if(NOT timing_retry_result EQUAL 0 OR timing_retry_err OR
   NOT timing_retry_row_count EQUAL 8 OR NOT timing_retry MATCHES
       "cycles=1 order=CTTCTCCT discard_cycle=0 cycle_mode=replacement cycle_index=2")
    message(FATAL_ERROR
        "preferred timing replacement-cycle fixture failed\n"
        "${timing_retry}\n${timing_retry_err}")
endif()

run_bench(timing_setup_result timing_setup timing_setup_err preferredtiming
    --N 4096 --bb 64 --preferred-attempt 1 --evict-bytes 4096
    --metric setup --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REGEX MATCHALL "\n4096,64,setup," timing_setup_rows
    "${timing_setup}")
list(LENGTH timing_setup_rows timing_setup_row_count)
string(REGEX MATCHALL
    "\n4096,64,setup,[0-3],[0-7],(control|candidate),[0-9]+,0,[0-9]+,0,"
    timing_setup_valid_rows "${timing_setup}")
list(LENGTH timing_setup_valid_rows timing_setup_valid_row_count)
if(NOT timing_setup_result EQUAL 0 OR timing_setup_err OR
   NOT timing_setup_row_count EQUAL 32 OR
   NOT timing_setup_valid_row_count EQUAL 32 OR NOT timing_setup MATCHES
       "metric=setup.*cycles=4 order=CTTCTCCT discard_cycle=0.*overhead=none.*payload=none" OR
   NOT timing_setup MATCHES ",0,0")
    message(FATAL_ERROR
        "preferred setup timing fixture failed\n${timing_setup}\n${timing_setup_err}")
endif()

run_bench(timing_noop_result timing_noop timing_noop_err preferredtiming
    --N 4096 --bb 1280 --preferred-attempt 0 --evict-bytes 4096
    --metric setup --route-cache "${route_k4096_wide_noop_path}"
    --route-cache-sha256 "${route_k4096_wide_noop_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REGEX MATCHALL
    "\n4096,1280,setup,[0-3],[0-7],(control|candidate),0,0,"
    timing_noop_rows "${timing_noop}")
list(LENGTH timing_noop_rows timing_noop_row_count)
if(NOT timing_noop_result EQUAL 0 OR timing_noop_err OR
   NOT timing_noop_row_count EQUAL 32)
    message(FATAL_ERROR
        "preferred wide no-op timing fixture failed\n"
        "${timing_noop}\n${timing_noop_err}")
endif()

expect_failure("requires a cached direct route or valid wide no-op alias"
    preferredtiming
    --N 4096 --bb 64 --preferred-attempt 0 --evict-bytes 4096
    --metric solve --route-cache "${route_k4096_p0_path}"
    --route-cache-sha256 "${route_k4096_p0_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("exact selected-attempt route row" preferredtiming
    --N 4096 --bb 64 --preferred-attempt 1 --evict-bytes 4096
    --metric solve --route-cache "${route_k4096_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("argument domain mismatch" preferredtiming
    --N 4096 --bb 256 --preferred-attempt 1 --evict-bytes 4096
    --metric solve --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("argument domain mismatch" preferredtiming
    --N 4096 --bb 64 --preferred-attempt 1 --evict-bytes 4096
    --metric solve --cycle-index 4
    --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("bad --N value" preferredtiming
    --N 04096 --bb 64 --preferred-attempt 1 --evict-bytes 4096
    --metric solve --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("bad --seed value" preferredtiming
    --N 4096 --bb 64 --preferred-attempt 1 --evict-bytes 4096
    --metric solve --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 0x1234 --schedule burst)

# Grouped timing uses one immutable, full-payload packet trace and reapplies
# every TLS experiment setting outside each timed solve.  The cold fixture
# compares the campaign's raw P48/r0 reference with its P48/r3 finalist in
# four exact ABBABAAB cycles.  Outcomes may differ under the hard trace, but
# every physical observation must reproduce its arm's preflight result so the
# common-success classification remains trustworthy.
run_bench(grouped_timing_result grouped_timing grouped_timing_err
    groupedtiming --N 4096 --bb 64 --overhead 4
    --control-geometry shared-x
    --control-period 48 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x
    --candidate-period 48 --candidate-grouped-rows 3
    --candidate-buckets separate --evict-bytes 4096 --cache-state cold
    --loss 0.5 --seed 4660 --schedule adversarial)
string(REGEX MATCHALL
    "\n4096,64,4,adversarial,4660,0.5,cold," grouped_timing_rows
    "${grouped_timing}")
list(LENGTH grouped_timing_rows grouped_timing_row_count)
string(REGEX MATCHALL
    "\n4096,64,4,adversarial,4660,0.5,cold,[0-3],[0-7],(control|candidate),(48),(0|3),(auto|separate),[0-9]+,0x[0-9a-f]+,0x[0-9a-f]+,[01],(common-success|control-only|candidate-only|common-failure),[01],[01],1,[0-9]+,0,"
    grouped_timing_valid_rows "${grouped_timing}")
list(LENGTH grouped_timing_valid_rows grouped_timing_valid_row_count)
if(NOT grouped_timing_result EQUAL 0 OR grouped_timing_err OR
   NOT grouped_timing_row_count EQUAL 32 OR
   NOT grouped_timing_valid_row_count EQUAL 32 OR
   NOT grouped_timing MATCHES
       "schema=v2.*timing_scope=solve.*cycles=4 order=ABBABAAB discard_cycle=0.*cycle_mode=full cycle_index=all.*overhead=4.*overhead_stream=salted.*control_period=48 control_grouped_rows=0 control_buckets=auto control_grouped_hash_seed=0x0 control_final_h_a_columns=0.*candidate_period=48 candidate_grouped_rows=3 candidate_buckets=separate candidate_grouped_hash_seed=0xb7e15162 candidate_final_h_a_columns=12.*dense_two_anchor=1 control_attempt=0 control_matrix_seed=0x136889600063cbf control_peel_seed=0x382fe3a7 candidate_attempt=0 candidate_matrix_seed=0x136889600063cbf candidate_peel_seed=0x382fe3a7.*payload=distinct-packet-zero-v1.*payload_count=4100.*payload_alignment=64 payload_prefaulted=1.*system_build=outside-timer tls_reapply=full-per-slot-outside-timer allocator_tls_state=preflight-warmed solve_value_storage=owned-noinit solve_value_publish=swap" OR
   NOT grouped_timing MATCHES
       "N,bb,overhead,schedule,seed,loss,cache_state,cycle,slot,arm,period,grouped_rows,buckets_requested,seed_attempt,matrix_seed,peel_seed,preflight_result,cell_class,common_success,result,outcome_stable,elapsed_ns,saturated,cpu_before,cpu_after,cpu_migrated,minflt_delta,majflt_delta,fault_contaminated,inactivated,binary_def,heavy_gain,block_xors,block_muladds,build_ns,peel_ns,project_ns,residual_ns,backsub_ns,joint_source_xors,joint_marginal_xors,joint_marginal_copies,joint_active_deltas,joint_scratch_bytes,dual_source_columns,source_bytes,packet_payload_bytes,intermediate_bytes,solve_value_arena_bytes,solve_value_eager_zero_bytes,solve_value_commit_copy_bytes" OR
   NOT grouped_timing MATCHES
       ",control,48,0,auto,0,0x136889600063cbf,0x382fe3a7,0,common-success,1,0,1,[0-9]+,0,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,117,8,8,90775,1973," OR
   NOT grouped_timing MATCHES
       ",candidate,48,3,separate,0,0x136889600063cbf,0x382fe3a7,0,common-success,1,0,1,[0-9]+,0,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,117,8,8,94847,1974,")
    message(FATAL_ERROR
        "grouped cold timing fixture failed\n"
        "${grouped_timing}\n${grouped_timing_err}")
endif()
foreach(cycle RANGE 0 3)
    foreach(slot RANGE 0 7)
        if(slot EQUAL 0 OR slot EQUAL 3 OR slot EQUAL 5 OR slot EQUAL 6)
            set(grouped_timing_arm "control")
        else()
            set(grouped_timing_arm "candidate")
        endif()
        if(NOT grouped_timing MATCHES
           "cold,${cycle},${slot},${grouped_timing_arm},")
            message(FATAL_ERROR
                "grouped timing order mismatch cycle=${cycle} slot=${slot}")
        endif()
    endforeach()
endforeach()

# Replacement cycles support independent warm panels and explicit dispatch
# comparisons without splicing partial output from a full run.
run_bench(grouped_warm_result grouped_warm grouped_warm_err groupedtiming
    --N 3200 --bb 64 --overhead 0
    --control-geometry shared-x
    --control-period 48 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x
    --candidate-period 32 --candidate-grouped-rows 7
    --candidate-buckets joint-delta --evict-bytes 4096 --cache-state warm
    --cycle-index 2 --loss 0.5 --seed 4660 --schedule repair-only)
string(REGEX MATCHALL
    "\n3200,64,0,repair-only,4660,0.5,warm,2," grouped_warm_rows
    "${grouped_warm}")
list(LENGTH grouped_warm_rows grouped_warm_row_count)
if(NOT grouped_warm_result EQUAL 0 OR grouped_warm_err OR
   NOT grouped_warm_row_count EQUAL 8 OR NOT grouped_warm MATCHES
       "cycles=1 order=ABBABAAB discard_cycle=0 cycle_mode=replacement cycle_index=2.*control_period=48 control_grouped_rows=0 control_buckets=auto.*candidate_period=32 candidate_grouped_rows=7 candidate_buckets=joint-delta.*dense_two_anchor=1" OR
   NOT grouped_warm MATCHES
       ",candidate,32,7,joint-delta,0,0x13a1a9dd5eb58b9d,0xf226e3bc,0,common-success,1,0,1,[0-9]+,0,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,-?[0-9]+,115,12,12,73665,1756,[0-9]+,[0-9]+,[0-9]+,[0-9]+,[0-9]+,3172,1984,64,32,6144,0,204800,204800,210304,210304,0,0")
    message(FATAL_ERROR
        "grouped warm replacement fixture failed\n"
        "${grouped_warm}\n${grouped_warm_err}")
endif()

# Dispatch timing needs intermediate practical payload widths to resolve the
# crossover between the original 64-byte and 1280-byte timing points.  Keep
# this bounded to the explicit sealed-campaign width set.
run_bench(grouped_width_result grouped_width_out grouped_width_err groupedtiming
    --N 4096 --bb 512 --overhead 4
    --control-geometry frozen
    --control-period 244 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x
    --candidate-period 32 --candidate-grouped-rows 7
    --candidate-buckets auto --evict-bytes 4096 --cache-state warm
    --cycle-index 2 --loss 0.5 --seed 15111065706836454659
    --schedule burst)
if(NOT grouped_width_result EQUAL 0 OR NOT grouped_width_err STREQUAL "" OR
   NOT grouped_width_out MATCHES
       "schema=v2.*N=4096 bb=512.*control_period=244 control_grouped_rows=0.*candidate_period=32 candidate_grouped_rows=7.*control_geometry=frozen.*control_rhs_route_expected=streamed.*candidate_geometry=shared-x.*candidate_rhs_route_expected=streamed" OR
   NOT grouped_width_out MATCHES
       ",control,244,0,auto,.*frozen,10,2,0,constant,0x0,0,0,0x4e,0x0,0x0,1,0,0x0,0,1,0,.*streamed,streamed" OR
   NOT grouped_width_out MATCHES
       ",candidate,32,7,auto,.*shared-x,10,2,0,constant,0x0,0,0,0x4e,0x0,0xb7e15162,1,0,0x0,0,1,0,.*streamed,streamed")
    message(FATAL_ERROR
        "groupedtiming intermediate-width smoke failed: ${grouped_width_err}\n${grouped_width_out}")
endif()

# The automatic grouped RHS route has two production crossover boundaries.
# Exercise both sides with the exact P32/r7 architecture and require every
# physical solve to receipt the actual route rather than inferring it from
# zero/nonzero work counters.
foreach(route_cell IN ITEMS
        "3199;4096;streamed"
        "3200;4096;joint-delta"
        "9999;1280;streamed"
        "10000;1280;joint-delta")
    list(GET route_cell 0 route_N)
    list(GET route_cell 1 route_bb)
    list(GET route_cell 2 route_expected)
    run_bench(route_result route_out route_err groupedtiming
        --N ${route_N} --bb ${route_bb} --overhead 4
        --control-geometry frozen
        --control-period 244 --control-grouped-rows 0 --control-buckets auto
        --candidate-geometry shared-x
        --candidate-period 32 --candidate-grouped-rows 7
        --candidate-buckets auto --evict-bytes 4096 --cache-state warm
        --cycle-index 2 --loss 0.5 --seed 15111065706836454659
        --schedule burst)
    string(REGEX MATCHALL
        ",candidate,32,7,auto,[^\n]*,${route_expected},${route_expected}"
        route_candidate_rows "${route_out}")
    list(LENGTH route_candidate_rows route_candidate_count)
    if(NOT route_result EQUAL 0 OR NOT route_err STREQUAL "" OR
       NOT route_out MATCHES
           "control_geometry=frozen.*control_rhs_route_expected=streamed.*candidate_geometry=shared-x.*candidate_rhs_route_expected=${route_expected}.*candidate_preflight_rhs_route=${route_expected}" OR
       NOT route_candidate_count EQUAL 4)
        message(FATAL_ERROR
            "groupedtiming automatic route boundary failed "
            "N=${route_N} bb=${route_bb} expected=${route_expected}: "
            "${route_err}\n${route_out}")
    endif()
endforeach()

expect_failure("requires --N" groupedtiming)
expect_failure("requires .*--control-geometry" groupedtiming
    --N 4096 --bb 64 --overhead 4
    --control-period 244 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x --candidate-period 32
    --candidate-grouped-rows 7 --candidate-buckets auto
    --evict-bytes 4096 --cache-state cold --loss 0.5 --seed 4660
    --schedule burst)
expect_failure("requires .*--candidate-geometry" groupedtiming
    --N 4096 --bb 64 --overhead 4
    --control-geometry frozen --control-period 244
    --control-grouped-rows 0 --control-buckets auto
    --candidate-period 32 --candidate-grouped-rows 7
    --candidate-buckets auto --evict-bytes 4096 --cache-state cold
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("argument domain mismatch" groupedtiming
    --N 4096 --bb 64 --overhead 4
    --control-geometry shared-x
    --control-period 11 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x
    --candidate-period 48 --candidate-grouped-rows 3
    --candidate-buckets separate --evict-bytes 4096 --cache-state cold
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("argument domain mismatch" groupedtiming
    --N 4096 --bb 64 --overhead 4
    --control-geometry shared-x
    --control-period 48 --control-grouped-rows 0
    --control-buckets separate
    --candidate-geometry shared-x
    --candidate-period 48 --candidate-grouped-rows 3
    --candidate-buckets separate --evict-bytes 4096 --cache-state cold
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("argument domain mismatch" groupedtiming
    --N 4096 --bb 64 --overhead 4
    --control-geometry shared-x
    --control-period 48 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x
    --candidate-period 48 --candidate-grouped-rows 10
    --candidate-buckets separate --evict-bytes 4096 --cache-state cold
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("bad --seed value" groupedtiming
    --N 4096 --bb 64 --overhead 4
    --control-geometry shared-x
    --control-period 48 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x
    --candidate-period 48 --candidate-grouped-rows 3
    --candidate-buckets separate --evict-bytes 4096 --cache-state cold
    --loss 0.5 --seed 0x1234 --schedule burst)
expect_failure("cache-state must be cold or warm" groupedtiming
    --N 4096 --bb 64 --overhead 4
    --control-geometry shared-x
    --control-period 48 --control-grouped-rows 0 --control-buckets auto
    --candidate-geometry shared-x
    --candidate-period 48 --candidate-grouped-rows 3
    --candidate-buckets separate --evict-bytes 4096 --cache-state tepid
    --loss 0.5 --seed 4660 --schedule burst)

run_bench(route_result route_mixed route_err preferredattempt --mode route
    --N 3,4096 --bb-list 64
    --preferred-map "3@64=none|4096@64=1"
    --route-context-sha256 "${route_context_sha256}")
if(NOT route_result EQUAL 0 OR NOT "${route_err}" STREQUAL "")
    message(FATAL_ERROR "cannot create mixed route fixture: ${route_err}")
endif()
set(route_mixed_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-mixed.csv")
file(WRITE "${route_mixed_path}" "${route_mixed}")
file(SHA256 "${route_mixed_path}" route_mixed_sha256)

expect_success("mode=paired.*systematic_probe_accounting=explicit"
    preferredattempt --mode paired --N 3 --bb-list 64
    --preferred-map "3@64=0" --route-cache "${route_k3_p0_path}"
    --route-cache-sha256 "${route_k3_p0_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.35 --seed 1 --schedule repair-only)
expect_success("3,64,control,-1,1,1,0,1,0,0,0,1,1,1"
    preferredattempt --mode paired --N 3 --bb-list 64
    --preferred-map "3@64=0" --route-cache "${route_k3_p0_path}"
    --route-cache-sha256 "${route_k3_p0_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.35 --seed 1 --schedule repair-only)
expect_success("3,64,candidate,0,1,1,1,0,1,0,0,0,1,1"
    preferredattempt --mode paired --N 3 --bb-list 64
    --preferred-map "3@64=0" --route-cache "${route_k3_p0_path}"
    --route-cache-sha256 "${route_k3_p0_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.35 --seed 1 --schedule repair-only)
expect_success("4,64,candidate,-1,[0-9]+,[0-9]+,0,1,0,1,0,0,[01],[01]"
    preferredattempt --mode paired --N 3,4 --bb-list 64
    --preferred-map "3@64=0" --route-cache "${route_k3_k4_p0_path}"
    --route-cache-sha256 "${route_k3_k4_p0_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.35 --seed 1 --schedule repair-only)
run_bench(mixed_result mixed_out mixed_err preferredattempt --mode paired
    --N 3,4096 --bb-list 64 --preferred-map "4096@64=1"
    --route-cache "${route_mixed_path}"
    --route-cache-sha256 "${route_mixed_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.35 --seed 1 --schedule repair-only)
set(mixed_rows 0)
set(mixed_physical_solves 0)
set(mixed_candidate_rows 0)
set(mixed_candidate_physical_solves 0)
string(REPLACE "\n" ";" mixed_lines "${mixed_out}")
foreach(line IN LISTS mixed_lines)
    if(line MATCHES "^[0-9]+,")
        math(EXPR mixed_rows "${mixed_rows} + 1")
        string(REPLACE "," ";" columns "${line}")
        list(GET columns 2 arm)
        list(GET columns 11 physical_solve)
        math(EXPR mixed_physical_solves
            "${mixed_physical_solves} + ${physical_solve}")
        if(arm STREQUAL "candidate")
            math(EXPR mixed_candidate_rows "${mixed_candidate_rows} + 1")
            math(EXPR mixed_candidate_physical_solves
                "${mixed_candidate_physical_solves} + ${physical_solve}")
        endif()
    endif()
endforeach()
if(NOT mixed_result EQUAL 0 OR mixed_err OR NOT mixed_rows EQUAL 4 OR
   NOT mixed_candidate_rows EQUAL 2 OR
   NOT mixed_physical_solves EQUAL 3 OR
   NOT mixed_candidate_physical_solves EQUAL 1)
    message(FATAL_ERROR
        "mixed preferred/control logical-physical accounting mismatch\n"
        "${mixed_out}\n${mixed_err}")
endif()
expect_success("4096,64,candidate,1,0,1,1,1,0,0,1,1,0,0"
    preferredattempt --mode paired --N 4096 --bb-list 64
    --preferred-map "4096@64=1" --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_success(
    "mode=candidate.*route_cache_sha256=${route_k4096_sha256}"
    preferredattempt --mode candidate --N 4096 --bb-list 64
    --preferred-map "4096@64=0,1" --route-cache "${route_k4096_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_success("4096,64,candidate,1,0,1,1,1,0,0,1,1,0,0"
    preferredattempt --mode candidate --N 4096 --bb-list 64
    --preferred-map "4096@64=0,1" --route-cache "${route_k4096_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("requires explicit --mode" preferredattempt
    --N 3 --bb-list 64)
expect_failure("specified more than once" preferredattempt --mode route
    --mode route --N 3 --bb-list 64 --preferred-map "3@64=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("canonical ascending decimal lists" preferredattempt --mode route
    --N 03 --bb-list 64 --preferred-map "3@64=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("canonical ascending decimal lists" preferredattempt --mode route
    --N 4,3 --bb-list 64
    --preferred-map "4@64=0|3@64=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("canonical K/width record order" preferredattempt --mode route
    --N 3,4 --bb-list 64
    --preferred-map "4@64=0|3@64=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("canonical K/width record order" preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "03@64=00"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("route mode rejects recovery-only" preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "3@64=0"
    --route-context-sha256 "${route_context_sha256}" --loss 0.5)
expect_failure("require explicit --loss" preferredattempt --mode control
    --N 3 --bb-list 64)
expect_failure("requires a nonempty map" preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map none
    --route-context-sha256 "${route_context_sha256}")
expect_failure("malformed --preferred-map" preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "3@64=0,0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("malformed --preferred-map" preferredattempt --mode route
    --N 3 --bb-list 64
    --preferred-map
    "3@64=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("outside --N x --bb-list" preferredattempt --mode route
    --N 3 --bb-list 64 --preferred-map "4@64=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("must cover the exact" preferredattempt --mode route
    --N 3,4 --bb-list 64 --preferred-map "3@64=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("unique even values" preferredattempt --mode control
    --N 3 --bb-list 63)
expect_failure("at most four widths" preferredattempt --mode control
    --N 3 --bb-list 2,4,6,8,10)
expect_failure("exactly one of --probe-route" preferredattempt
    --mode candidate --N 4096 --bb-list 64
    --preferred-map "4096@64=1"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("attempt rows must exactly match" preferredattempt
    --mode candidate --N 3 --bb-list 64 --preferred-map "3@64=0"
    --route-cache "${route_k3_only_path}"
    --route-cache-sha256 "${route_k3_only_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("selected attempt row must exactly match" preferredattempt
    --mode paired --N 3 --bb-list 64 --preferred-map "3@64=0"
    --route-cache "${route_k3_only_path}"
    --route-cache-sha256 "${route_k3_only_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_success("probe_route=1" preferredattempt
    --mode candidate --N 4096 --bb-list 64
    --preferred-map "4096@64=0,1" --probe-route
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_success(
    "preferred_route: N=4096 bb=64 p=0 a0=0 actual=0 valid=1 fallback=0 no_op=1 direct=0 preferred_probe_solves=0"
    preferredattempt --mode candidate --N 4096 --bb-list 64
    --preferred-map "4096@64=0,1" --probe-route
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_success(
    "4096,64,candidate,1,0,1,1,1,0,0,1,1,0,0"
    preferredattempt --mode candidate --N 4096 --bb-list 64
    --preferred-map "4096@64=0,1" --probe-route
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
run_bench(probe_result probe_out probe_err preferredattempt --mode candidate
    --N 4096 --bb-list 64 --preferred-map "4096@64=0,1" --probe-route
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REGEX MATCHALL "canonical_probe_solves=" canonical_probe_events
    "${probe_out}")
string(REGEX MATCHALL "preferred_probe_solves=" preferred_probe_events
    "${probe_out}")
list(LENGTH canonical_probe_events canonical_probe_event_count)
list(LENGTH preferred_probe_events preferred_probe_event_count)
if(NOT probe_result EQUAL 0 OR probe_err OR
   NOT canonical_probe_event_count EQUAL 1 OR
   NOT preferred_probe_event_count EQUAL 2)
    message(FATAL_ERROR
        "preferred probe accounting is not one-event additive\n${probe_out}\n"
        "${probe_err}")
endif()
run_bench(alias_result alias_out alias_err preferredattempt --mode candidate
    --N 3 --bb-list 64
    --preferred-map "3@64=0,1,2" --route-cache "${route_k3_only_path}"
    --route-cache-sha256 "${route_k3_only_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REGEX MATCHALL "preferred_candidate_alias: N=3 bb=64 p=[012][^\r\n]*physical_solve=0"
    alias_events "${alias_out}")
string(REGEX MATCHALL "(^|\n)3,64," alias_numeric_rows "${alias_out}")
string(REGEX MATCHALL "preferred_probe_accounting" alias_probe_events
    "${alias_out}")
list(LENGTH alias_events alias_event_count)
list(LENGTH alias_numeric_rows alias_numeric_row_count)
list(LENGTH alias_probe_events alias_probe_event_count)
if(NOT alias_result EQUAL 0 OR alias_err OR
   NOT alias_event_count EQUAL 3 OR NOT alias_numeric_row_count EQUAL 0 OR
   NOT alias_probe_event_count EQUAL 0)
    message(FATAL_ERROR
        "cached all-alias accounting mismatch\n${alias_out}\n${alias_err}")
endif()

run_bench(candidate_mixed_result candidate_mixed_out candidate_mixed_err
    preferredattempt --mode candidate --N 4096 --bb-list 64
    --preferred-map "4096@64=0,1" --route-cache "${route_k4096_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REGEX MATCHALL
    "preferred_candidate_alias: N=4096 bb=64 p=0[^\r\n]*physical_solve=0"
    candidate_mixed_aliases "${candidate_mixed_out}")
string(REGEX MATCHALL "(^|\n)4096,64,candidate,1,[^\r\n]*"
    candidate_mixed_rows "${candidate_mixed_out}")
string(REGEX MATCHALL "preferred_probe_accounting" candidate_mixed_probes
    "${candidate_mixed_out}")
list(LENGTH candidate_mixed_aliases candidate_mixed_alias_count)
list(LENGTH candidate_mixed_rows candidate_mixed_row_count)
list(LENGTH candidate_mixed_probes candidate_mixed_probe_count)
if(NOT candidate_mixed_result EQUAL 0 OR candidate_mixed_err OR
   NOT candidate_mixed_alias_count EQUAL 1 OR
   NOT candidate_mixed_row_count EQUAL 1 OR
   NOT candidate_mixed_probe_count EQUAL 0 OR
   NOT candidate_mixed_out MATCHES
       "4096,64,candidate,1,0,1,1,1,0,0,1,1,0,0")
    message(FATAL_ERROR
        "cached mixed alias/direct accounting mismatch\n"
        "${candidate_mixed_out}\n${candidate_mixed_err}")
endif()
expect_success("3,64,candidate,1,1,1,1,1,0,1,0,0"
    preferredattempt --mode paired --N 3 --bb-list 64
    --preferred-map "3@64=1" --route-cache "${route_k3_p1_path}"
    --route-cache-sha256 "${route_k3_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_success("3,64,candidate,2,1,1,1,0,1,0,0,0"
    preferredattempt --mode paired --N 3 --bb-list 64
    --preferred-map "3@64=2" --route-cache "${route_k3_p2_path}"
    --route-cache-sha256 "${route_k3_p2_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)

# The native reader independently enforces the same ordered-group and
# first-row-only canonical probe charge accepted by the frozen Python parser.
# Rehashing semantically noncanonical bytes must not turn them into a cache.
string(REPLACE "3,64,preferred,0,1,1,0,1,0,0,2,0"
    "3,64,preferred,0,1,1,0,1,0,0,0,0"
    route_shifted_charge "${route_k3_only}")
string(REPLACE "3,64,preferred,1,1,1,1,0,1,0,0,0"
    "3,64,preferred,1,1,1,1,0,1,0,2,0"
    route_shifted_charge "${route_shifted_charge}")
set(route_shifted_charge_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-shifted-charge.csv")
file(WRITE "${route_shifted_charge_path}" "${route_shifted_charge}")
file(SHA256 "${route_shifted_charge_path}" route_shifted_charge_sha256)
expect_failure("first-row canonical probe accounting" preferredattempt
    --mode candidate --N 3 --bb-list 64 --preferred-map "3@64=0,1,2"
    --route-cache "${route_shifted_charge_path}"
    --route-cache-sha256 "${route_shifted_charge_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)

string(REPLACE "\n" ";" route_k3_lines "${route_k3}")
list(GET route_k3_lines 0 route_k3_preamble)
list(GET route_k3_lines 1 route_k3_header)
list(GET route_k3_lines 2 route_k3_p0)
list(GET route_k3_lines 3 route_k3_p1)
list(GET route_k3_lines 4 route_k3_p2)
list(GET route_k3_lines 5 route_k4_control)
set(route_reordered
    "${route_k3_preamble}\n${route_k3_header}\n${route_k4_control}\n${route_k3_p0}\n${route_k3_p1}\n${route_k3_p2}\n")
set(route_reordered_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-reordered.csv")
file(WRITE "${route_reordered_path}" "${route_reordered}")
file(SHA256 "${route_reordered_path}" route_reordered_sha256)
expect_failure("key groups are reordered" preferredattempt --mode paired
    --N 3,4 --bb-list 64 --preferred-map "3@64=0"
    --route-cache "${route_reordered_path}"
    --route-cache-sha256 "${route_reordered_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)

# The digest is over the actual bytes and the parser independently rejects a
# self-consistently rehashed forged p<a0 direct classification.
set(route_tampered_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-tampered.csv")
file(WRITE "${route_tampered_path}" "${route_k4096}#")
expect_failure("SHA-256 mismatch" preferredattempt --mode candidate
    --N 4096 --bb-list 64 --preferred-map "4096@64=0,1"
    --route-cache "${route_tampered_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
string(REPLACE "3,64,preferred,0,1,1,0,1,0,0,2,0"
    "3,64,preferred,0,1,0,1,0,0,1,2,0"
    route_forged "${route_k3_only}")
set(route_forged_path
    "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-forged.csv")
file(WRITE "${route_forged_path}" "${route_forged}")
file(SHA256 "${route_forged_path}" route_forged_sha256)
expect_failure("inconsistent classification" preferredattempt
    --mode candidate --N 3 --bb-list 64 --preferred-map "3@64=0,1,2"
    --route-cache "${route_forged_path}"
    --route-cache-sha256 "${route_forged_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("selected attempt row must exactly match" preferredattempt
    --mode paired --N 4096 --bb-list 64 --preferred-map "4096@64=2"
    --route-cache "${route_k4096_p1_path}"
    --route-cache-sha256 "${route_k4096_p1_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("route cache keys must exactly match" preferredattempt
    --mode candidate --N 3 --bb-list 64 --preferred-map "3@64=0,1,2"
    --route-cache "${route_k3_path}"
    --route-cache-sha256 "${route_k3_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_failure("requires explicit --preferred-map" preferredattempt
    --mode paired --N 4096 --bb-list 64
    --route-cache "${route_k4096_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("map/cache route-status mismatch" preferredattempt
    --mode paired --N 4096 --bb-list 64 --preferred-map none
    --route-cache "${route_k4096_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)

# Exactly 32 candidates and multiple K/width records are accepted without
# reordering or dropping records; the independent cache parser also rejects a
# validly classified, self-consistently rehashed 33rd attempt for one key.
run_bench(route_32_result route_32 route_32_err preferredattempt --mode route
    --N 3 --bb-list 64
    --preferred-map
    "3@64=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31"
    --route-context-sha256 "${route_context_sha256}")
if(NOT route_32_result EQUAL 0 OR route_32_err OR
   NOT route_32 MATCHES "3,64,preferred,31,")
    message(FATAL_ERROR "32-attempt route fixture failed\n${route_32}\n${route_32_err}")
endif()
run_bench(route_33_tail_result route_33_tail route_33_tail_err
    preferredattempt --mode route --N 3 --bb-list 64
    --preferred-map "3@64=32"
    --route-context-sha256 "${route_context_sha256}")
string(REPLACE "\n" ";" route_33_tail_lines "${route_33_tail}")
list(GET route_33_tail_lines 2 route_33_tail_row)
string(REGEX REPLACE ",[0-9]+,([0-9]+)$" ",0,\\1"
    route_33_tail_row "${route_33_tail_row}")
set(route_33 "${route_32}${route_33_tail_row}\n")
set(route_33_path "${CMAKE_CURRENT_BINARY_DIR}/preferred-route-33.csv")
file(WRITE "${route_33_path}" "${route_33}")
file(SHA256 "${route_33_path}" route_33_sha256)
if(NOT route_33_tail_result EQUAL 0 OR route_33_tail_err)
    message(FATAL_ERROR "cannot create 33rd route row: ${route_33_tail_err}")
endif()
expect_failure("more than 32 attempts" preferredattempt --mode candidate
    --N 3 --bb-list 64 --preferred-map "3@64=0"
    --route-cache "${route_33_path}"
    --route-cache-sha256 "${route_33_sha256}"
    --route-context-sha256 "${route_context_sha256}"
    --loss 0.5 --seed 4660 --schedule burst)
expect_success("4,256,preferred,0," preferredattempt --mode route --N 3,4
    --bb-list 64,256
    --preferred-map "3@64=0|3@256=0|4@64=0|4@256=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("logical p list must be identical" preferredattempt
    --mode route --N 4096 --bb-list 64,256
    --preferred-map "4096@64=0,1|4096@256=0"
    --route-context-sha256 "${route_context_sha256}")
expect_failure("supports only 64,256,1280,4096" preferredattempt
    --mode control --N 3 --bb-list 128)
expect_failure("magic/schema mismatch" preferredattempt --mode candidate
    --N 4096 --bb-list 64 --preferred-map "4096@64=0,1"
    --route-cache "${route_k4096_path}"
    --route-cache-sha256 "${route_k4096_sha256}"
    --route-context-sha256
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
    --loss 0.5 --seed 4660 --schedule burst)

# The dedicated q0 control is equation-identical to the historical
# precodefail path on both sides of the adaptive two-anchor cutoff, at every
# development width.  Timing columns are deliberately ignored.
run_bench(identity_result identity_out identity_err preferredattempt
    --mode control --N 4095,4096 --bb-list 64,256,1280,4096
    --loss 0.5 --seed 4660 --schedule burst)
run_bench(legacy_4095_result legacy_4095 legacy_4095_err precodefail
    --N 4095 --bb-list 64,256,1280,4096 --overhead 0 --trials 1
    --threads 1 --loss 0.5 --seed 4660 --schedule burst
    --completion mixed --mix-count 2)
run_bench(legacy_4096_result legacy_4096 legacy_4096_err precodefail
    --N 4096 --bb-list 64,256,1280,4096 --overhead 0 --trials 1
    --threads 1 --loss 0.5 --seed 4660 --schedule burst
    --completion mixed --mix-count 2 --binary-dense-two-anchor)
if(NOT identity_result EQUAL 0 OR NOT legacy_4095_result EQUAL 0 OR
   NOT legacy_4096_result EQUAL 0 OR identity_err OR legacy_4095_err OR
   legacy_4096_err OR NOT identity_out MATCHES
       "mode=control.*route_context_sha256=none")
    message(FATAL_ERROR
        "q0 boundary identity command failed\n"
        "dedicated=${identity_err}\nlegacy4095=${legacy_4095_err}\n"
        "legacy4096=${legacy_4096_err}")
endif()
foreach(spec IN ITEMS
    "4095|64|129|97043|4465" "4095|256|139|97089|4588"
    "4095|1280|130|95407|4476" "4095|4096|137|95901|4565"
    "4096|64|147|96644|4681" "4096|256|129|98043|4465"
    "4096|1280|148|96289|4694" "4096|4096|144|96038|4651")
    string(REPLACE "|" ";" fields "${spec}")
    list(GET fields 0 identity_K)
    list(GET fields 1 identity_bb)
    list(GET fields 2 identity_inact)
    list(GET fields 3 identity_xors)
    list(GET fields 4 identity_muladds)
    set(dedicated_pattern
        "${identity_K},${identity_bb},control,-1,0,0,0,1,0,0,0,1,0,0,0,0,${identity_inact},12,12,${identity_xors},${identity_muladds}")
    if(identity_K STREQUAL "4095")
        set(legacy_out "${legacy_4095}")
    else()
        set(legacy_out "${legacy_4096}")
    endif()
    set(legacy_pattern
        "${identity_K},${identity_bb},periodic,2,0,1,1,0,0,0\\.00000000,${identity_inact}\\.000,${identity_inact},12\\.000,12,12\\.000,12,0,[^\r\n]*,0,${identity_xors}\\.000,${identity_muladds}\\.000,")
    if(NOT identity_out MATCHES "${dedicated_pattern}" OR
       NOT legacy_out MATCHES "${legacy_pattern}")
        message(FATAL_ERROR
            "q0 identity mismatch for K=${identity_K} bb=${identity_bb}\n"
            "dedicated=${identity_out}\nlegacy=${legacy_out}")
    endif()
endforeach()
expect_failure("schedule must be burst" preferredattempt --mode control
    --N 3 --bb-list 64 --schedule iid)
expect_failure("requires 12 binary dense rows" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --binary-dense-two-anchor --binary-dense-rows 13)
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
