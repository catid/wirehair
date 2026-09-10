# Compare complete deterministic streams only after both producers succeed.
foreach(_arm IN ITEMS BASELINE CANDIDATE)
    execute_process(COMMAND "${${_arm}}"
        OUTPUT_FILE "${OUTPUT_DIR}/certified-${_arm}.bin"
        ERROR_FILE "${OUTPUT_DIR}/certified-${_arm}.stderr"
        RESULT_VARIABLE _result TIMEOUT 50)
    if(NOT "${_result}" STREQUAL "0")
        message(FATAL_ERROR "${_arm} compatibility producer failed: ${_result}")
    endif()
    file(SIZE "${OUTPUT_DIR}/certified-${_arm}.stderr" _stderr)
    file(SIZE "${OUTPUT_DIR}/certified-${_arm}.bin" _size)
    file(SHA256 "${OUTPUT_DIR}/certified-${_arm}.bin" _hash)
    if(NOT _stderr EQUAL 0 OR NOT _size EQUAL 2180292 OR
       NOT _hash STREQUAL "2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b")
        message(FATAL_ERROR "${_arm} differs from the preserved 48-case certified bytes")
    endif()
endforeach()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${OUTPUT_DIR}/certified-BASELINE.bin" "${OUTPUT_DIR}/certified-CANDIDATE.bin"
    RESULT_VARIABLE _same)
if(NOT "${_same}" STREQUAL "0")
    message(FATAL_ERROR "Certified baseline/candidate streams differ")
endif()
