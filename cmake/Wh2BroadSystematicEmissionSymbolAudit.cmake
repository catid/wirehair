if(NOT DEFINED NM OR NM STREQUAL "")
    message(FATAL_ERROR
        "broad-systematic-emission symbol audit requires CMAKE_NM")
endif()
if(NOT DEFINED BROAD_SYSTEMATIC_EMISSION OR
   NOT EXISTS "${BROAD_SYSTEMATIC_EMISSION}")
    message(FATAL_ERROR
        "broad-systematic-emission symbol audit is missing its executable")
endif()

execute_process(
    COMMAND "${NM}" -C "${BROAD_SYSTEMATIC_EMISSION}"
    RESULT_VARIABLE _nm_result
    OUTPUT_VARIABLE _symbols
    ERROR_VARIABLE _nm_error)
if(NOT _nm_result EQUAL 0)
    message(FATAL_ERROR
        "nm failed for broad-systematic-emission screen: ${_nm_error}")
endif()
string(STRIP "${_symbols}" _symbols_stripped)
string(TOLOWER "${_symbols_stripped}" _symbols_lower)
if(_symbols_stripped STREQUAL "" OR _symbols_lower MATCHES "no symbols")
    message(FATAL_ERROR
        "nm returned no auditable symbols for broad-systematic-emission screen")
endif()
string(REGEX MATCH "[ \t][Tt][ \t]+_?main([\r\n]|$)"
    _main_symbol "${_symbols}")
if(_main_symbol STREQUAL "")
    message(FATAL_ERROR
        "broad-systematic-emission screen lacks the positive main symbol")
endif()
# `main` is the one final-link positive symbol required in both ordinary and
# LTO layouts.  Implementation-body/layout requirements are intentionally
# delegated to the common LTO-aware timing-policy audit.

# Mirror the timed-code forbidden-state families from the common policy
# audit.  This screen is intentionally current-only and counter-free: neither
# arm may inherit a test-hook switch, TLS counter, allocation injector, or
# timing probe merely because it is linked as a benchmark executable.
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
    "ColdSolveWideXorTestMode"
    "ColdSolveWideXorObservationCount"
    "LastColdSolveWideXorSelection"
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
    "TinyPeriodicHeavyUses"
    "TinyPeriodicHeavyTransposeUseCount"
    "TinyPeriodicHeavyLegacyUseCount"
    "TinyPeriodicHeavyTimingEnabled"
    "TinyPeriodicHeavyTimedCalls"
    "TinyPeriodicHeavyTimedNanoseconds"
    "TinyPeriodicHeavyTimedDataRows"
    "PackedBinaryResidualTestMode"
    "PackedBinaryResidualUseCount"
    "ResumeSystemFingerprintChecks")

foreach(_forbidden IN LISTS _forbidden_symbols)
    string(FIND "${_symbols}" "${_forbidden}" _forbidden_offset)
    if(NOT _forbidden_offset EQUAL -1)
        message(FATAL_ERROR
            "broad-systematic-emission screen contains test-only symbol "
            "${_forbidden}")
    endif()
endforeach()

message(STATUS
    "counter-free broad-systematic-emission symbol audit passed")
