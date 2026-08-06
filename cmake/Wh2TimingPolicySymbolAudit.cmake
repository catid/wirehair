if(NOT DEFINED NM OR NM STREQUAL "")
    message(FATAL_ERROR "timing-policy symbol audit requires CMAKE_NM")
endif()

foreach(_artifact_var IN ITEMS TIMING_POLICY CONTRACT_WORKER PHASE_TIMING)
    if(NOT DEFINED ${_artifact_var} OR
       NOT EXISTS "${${_artifact_var}}")
        message(FATAL_ERROR
            "timing-policy symbol audit is missing ${_artifact_var}")
    endif()
    execute_process(
        COMMAND "${NM}" -C "${${_artifact_var}}"
        RESULT_VARIABLE _nm_result
        OUTPUT_VARIABLE _symbols
        ERROR_VARIABLE _nm_error)
    if(NOT _nm_result EQUAL 0)
        message(FATAL_ERROR
            "nm failed for ${_artifact_var}: ${_nm_error}")
    endif()
    set(${_artifact_var}_SYMBOLS "${_symbols}")
endforeach()

# This list mirrors every timed test-only state family in the codec.  A broad
# ForTesting rejection covers public entry points, while the internal names
# catch retained TLS/atomic state in the defining archive and final binaries.
set(_forbidden_symbols
    "ForTesting"
    "CodecAllocationFailureCountdown"
    "AllocationFailureCountdown"
    "HeavyBucketStorageLimit"
    "ActiveEncodeAllocationFailure"
    "EncodeAllocationFailureHits"
    "DecoderAllocationFailureCountdown"
    "DecoderColdSystematicReuseEnabled"
    "ColdReceiveAllocationBytesForTesting"
    "DecoderIncrementalResumeEnabled"
    "DecoderReceiveCounters"
    "OddPacketPeelSeedXor"
    "PacketRowSeedMultiplier"
    "PacketRowSeedAvalanche"
    "ActiveSolveAllocationFailure"
    "BinaryPeelOracleUsers"
    "BinaryPeelOracleComparisons"
    "HeavyProjectionOracleUsers"
    "HeavyProjectionOracleComparisons"
    "HeavyProjectionLegacyFallbacks"
    "ProjectionAVX2TestMode"
    "ProjectionAVX2BatchUseCount"
    "ProjectionFallbackBatchUseCount"
    "PackedBinaryResidualTestMode"
    "PackedBinaryResidualUseCount"
    "ResumeSystemFingerprintChecks")

foreach(_artifact_var IN ITEMS TIMING_POLICY CONTRACT_WORKER PHASE_TIMING)
    foreach(_forbidden IN LISTS _forbidden_symbols)
        string(FIND "${${_artifact_var}_SYMBOLS}" "${_forbidden}"
            _forbidden_offset)
        if(NOT _forbidden_offset EQUAL -1)
            message(FATAL_ERROR
                "${_artifact_var} contains test-only symbol ${_forbidden}")
        endif()
    endforeach()
endforeach()

# Require the capability in its defining archive.  Optimizing linkers may
# legitimately inline and internalize this wrapper in final executables; the
# consumers also carry compile-time BENCHMARK_EQUATIONS assertions.
string(FIND "${TIMING_POLICY_SYMBOLS}"
    "InitializeForValidatedSystemForBenchmark" _benchmark_offset)
if(_benchmark_offset EQUAL -1)
    message(FATAL_ERROR
        "TIMING_POLICY lacks the exact benchmark-system initializer")
endif()

message(STATUS
    "counter-free WH2 timing policy symbol audit passed")
