if(NOT DEFINED NM OR NM STREQUAL "")
    message(FATAL_ERROR
        "direct-systematic-complement symbol audit requires CMAKE_NM")
endif()
if(NOT DEFINED DIRECT_SYSTEMATIC_COMPLEMENT OR
   NOT EXISTS "${DIRECT_SYSTEMATIC_COMPLEMENT}")
    message(FATAL_ERROR
        "direct-systematic-complement symbol audit is missing its executable")
endif()

execute_process(
    COMMAND "${NM}" -C "${DIRECT_SYSTEMATIC_COMPLEMENT}"
    RESULT_VARIABLE _nm_result
    OUTPUT_VARIABLE _symbols
    ERROR_VARIABLE _nm_error)
if(NOT _nm_result EQUAL 0)
    message(FATAL_ERROR
        "nm failed for direct-systematic-complement screen: ${_nm_error}")
endif()
string(STRIP "${_symbols}" _symbols_stripped)
string(TOLOWER "${_symbols_stripped}" _symbols_lower)
if(_symbols_stripped STREQUAL "" OR _symbols_lower MATCHES "no symbols")
    message(FATAL_ERROR
        "nm returned no auditable symbols for direct-systematic-complement screen")
endif()
string(REGEX MATCH "[ \t][Tt][ \t]+_?main([\r\n]|$)"
    _main_symbol "${_symbols}")
if(_main_symbol STREQUAL "")
    message(FATAL_ERROR
        "direct-systematic-complement screen lacks the positive main symbol")
endif()
# `main` is the one final-link positive symbol required in both ordinary and
# LTO layouts.  The shared timing-policy audit owns implementation-body and
# isolated-section layout checks.

# Mirror every timed test-only state family from the common policy audit.  The
# complement worker has no runtime arm selector and must stay counter-free.
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
    "ResumeSystemFingerprintChecks"
    "RepairInvocation")

foreach(_forbidden IN LISTS _forbidden_symbols)
    string(FIND "${_symbols}" "${_forbidden}" _forbidden_offset)
    if(NOT _forbidden_offset EQUAL -1)
        message(FATAL_ERROR
            "direct-systematic-complement screen contains test-only symbol "
            "${_forbidden}")
    endif()
endforeach()

message(STATUS
    "counter-free direct-systematic-complement symbol audit passed")
