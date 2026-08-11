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
    "SingleWordProjectionTestMode"
    "SingleWordProjectionUseCount"
    "GeneralProjectionUseCount"
    "TinyPeriodicHeavyTransposeTestMode"
    "TinyPeriodicHeavyTransposeUseCount"
    "TinyPeriodicHeavyLegacyUseCount"
    "TinyPeriodicHeavyTimingEnabled"
    "TinyPeriodicHeavyTimedCalls"
    "TinyPeriodicHeavyTimedNanoseconds"
    "TinyPeriodicHeavyTimedDataRows"
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

# The small-rank specialization is deliberately isolated from the generic
# solve text.  GCC IPA cloning can silently move templated noinline functions
# back into .text, which measurably perturbs the untouched R >= 65 path.
if(DEFINED SYSTEM_NAME AND SYSTEM_NAME STREQUAL "Linux")
    if(NOT DEFINED OBJDUMP OR OBJDUMP STREQUAL "")
        message(FATAL_ERROR
            "Linux timing-policy layout audit requires CMAKE_OBJDUMP")
    endif()
    string(TOUPPER "${LTO_MODE}" _lto_mode)
    if(_lto_mode STREQUAL "" OR _lto_mode STREQUAL "OFF")
        set(_layout_artifacts TIMING_POLICY)
        set(_require_input_section TRUE)
    else()
        # Slim LTO archives contain compiler IR rather than native symbol-table
        # rows.  Audit both final consumers instead: the output linker merges
        # .text.* sections, but three distinct out-of-line bodies must remain.
        set(_layout_artifacts CONTRACT_WORKER PHASE_TIMING)
        set(_require_input_section FALSE)
    endif()
    foreach(_artifact_var IN LISTS _layout_artifacts)
        execute_process(
            COMMAND "${OBJDUMP}" -tC "${${_artifact_var}}"
            RESULT_VARIABLE _objdump_result
            OUTPUT_VARIABLE _object_symbols
            ERROR_VARIABLE _objdump_error)
        if(NOT _objdump_result EQUAL 0)
            message(FATAL_ERROR
                "objdump failed for ${_artifact_var}: ${_objdump_error}")
        endif()
        # Match function rows only.  Coverage builds also emit gcov object
        # symbols whose names contain the function name; those are metadata,
        # not escaped code bodies.
        string(REGEX MATCHALL
            "[^\n]*[ \t]F[ \t]+[^\n]*AccumulateSingleWord[^\n]*"
            _single_word_symbols "${_object_symbols}")
        list(LENGTH _single_word_symbols _single_word_symbol_count)
        if(NOT _single_word_symbol_count EQUAL 3)
            message(FATAL_ERROR
                "${_artifact_var} expected three single-word accumulator "
                "function symbols, found ${_single_word_symbol_count}")
        endif()
        foreach(_symbol IN LISTS _single_word_symbols)
            if(_require_input_section)
                string(FIND "${_symbol}" ".text.wh2_single_word"
                    _single_word_section_offset)
                if(_single_word_section_offset EQUAL -1)
                    message(FATAL_ERROR
                        "single-word accumulator escaped isolated text: "
                        "${_symbol}")
                endif()
            elseif(_symbol MATCHES
                   "\\.(constprop|isra|part)(\\.|\\]|\\))")
                message(FATAL_ERROR
                    "LTO cloned a single-word accumulator: ${_symbol}")
            endif()
        endforeach()

        string(REGEX MATCHALL
            "[^\n]*[ \t]F[ \t]+[^\n]*PrepareTinyPeriodicHeavy(Transposed|Legacy|Selected)ByPeel[^\n]*"
            _tiny_periodic_symbols "${_object_symbols}")
        list(LENGTH _tiny_periodic_symbols _tiny_periodic_symbol_count)
        if(NOT _tiny_periodic_symbol_count EQUAL 3)
            message(FATAL_ERROR
                "${_artifact_var} expected three tiny-periodic dispatch/helper "
                "function symbols, found ${_tiny_periodic_symbol_count}")
        endif()
        foreach(_symbol IN LISTS _tiny_periodic_symbols)
            if(_require_input_section)
                string(FIND "${_symbol}" ".text.wh2_tiny_periodic"
                    _tiny_periodic_section_offset)
                if(_tiny_periodic_section_offset EQUAL -1)
                    message(FATAL_ERROR
                        "tiny-periodic helper escaped isolated text: "
                        "${_symbol}")
                endif()
            elseif(_symbol MATCHES
                   "\\.(constprop|isra|part)(\\.|\\]|\\))")
                message(FATAL_ERROR
                    "LTO cloned a tiny-periodic helper: ${_symbol}")
            endif()
        endforeach()
    endforeach()
endif()

message(STATUS
    "counter-free WH2 timing policy symbol audit passed")
