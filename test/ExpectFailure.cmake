foreach(required IN ITEMS BINARY EXPECT_EXIT EXPECT_PATTERN)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

execute_process(
    COMMAND "${BINARY}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err)
set(combined "${out}${err}")
message(STATUS "${combined}")

if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL EXPECT_EXIT)
    message(FATAL_ERROR
        "expected exit ${EXPECT_EXIT}, got '${result}'\n${combined}")
endif()
if(NOT combined MATCHES "${EXPECT_PATTERN}")
    message(FATAL_ERROR
        "missing expected diagnostic '${EXPECT_PATTERN}'\n${combined}")
endif()
if(combined MATCHES
   "AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer|MemorySanitizer|ThreadSanitizer|runtime error:")
    message(FATAL_ERROR "sanitizer failure in expected-failure path\n${combined}")
endif()
