cmake_minimum_required(VERSION 3.15)

# Parser-only subprocess entry is confined to this selftest. The production
# checker always invokes the actual nm executable on the actual library.
include("${CMAKE_CURRENT_LIST_DIR}/ParseElfExports.cmake")
if(DEFINED PARSER_FIXTURE)
    file(READ "${PARSER_FIXTURE}" fixture)
    wirehair_parse_elf_exports("${fixture}" "TEST_ABI.1" actual)
    if(NOT actual STREQUAL "wirehair_decode@@TEST_ABI.1")
        message(FATAL_ERROR "Synthetic exact export mismatch: '${actual}'")
    endif()
    return()
endif()

foreach(required IN ITEMS CHECKER LIBRARY MANIFEST NM WORK_DIR)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

function(check_parser name text succeeds)
    set(fixture_path "${WORK_DIR}/${name}.txt")
    file(WRITE "${fixture_path}" "${text}")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" "-DPARSER_FIXTURE=${fixture_path}"
            -P "${CMAKE_CURRENT_LIST_FILE}"
        RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(succeeds)
        if(NOT result EQUAL 0)
            message(FATAL_ERROR "Parser rejected ${name}: ${out}${err}")
        endif()
    elseif(result EQUAL 0)
        message(FATAL_ERROR "Parser accepted invalid ${name}")
    elseif(NOT "${out}${err}" MATCHES
           "(Invalid or duplicate ELF version marker|Synthetic exact export mismatch|Unrecognized nm output line)")
        message(FATAL_ERROR "Parser failed ${name} for an unrelated reason: ${out}${err}")
    endif()
endfunction()
set(callable "wirehair_decode@@TEST_ABI.1 T 123 7\n")
check_parser(no_marker "${callable}" TRUE)
check_parser(gnu_marker "TEST_ABI.1 A 0\n${callable}" TRUE)
check_parser(llvm_marker "TEST_ABI.1@@TEST_ABI.1 A 0 0\n${callable}" TRUE)
check_parser(padded_marker "TEST_ABI.1 A 0000 0000\r\n${callable}" TRUE)
check_parser(wrong_type "TEST_ABI.1 T 0\n${callable}" FALSE)
check_parser(wrong_decorated_type "TEST_ABI.1@@TEST_ABI.1 D 0 0\n${callable}" FALSE)
check_parser(nonzero_value "TEST_ABI.1 A 1\n${callable}" FALSE)
check_parser(nonzero_size "TEST_ABI.1@@TEST_ABI.1 A 0 1\n${callable}" FALSE)
check_parser(malformed_marker "TEST_ABI.1 A 0 0 extra\n${callable}" FALSE)
check_parser(duplicate_marker "TEST_ABI.1 A 0\nTEST_ABI.1@@TEST_ABI.1 A 0 0\n${callable}" FALSE)
check_parser(wrong_version "TEST_ABI.1@@TEST_ABI.2 A 0 0\n${callable}" FALSE)
check_parser(single_at "TEST_ABI.1@TEST_ABI.1 A 0 0\n${callable}" FALSE)
check_parser(lookalike "TEST_ABIx1 A 0\n${callable}" FALSE)
check_parser(missing_function "TEST_ABI.1 A 0\n" FALSE)
check_parser(extra_function "${callable}wirehair_encode@@TEST_ABI.1 T 456 7\n" FALSE)
check_parser(unversioned_function "wirehair_decode T 123 7\n" FALSE)
check_parser(wrong_function_version "wirehair_decode@@TEST_ABI.2 T 123 7\n" FALSE)
check_parser(duplicate_function "${callable}${callable}" FALSE)
check_parser(malformed_line "not-an-nm-line\n${callable}" FALSE)

file(READ "${MANIFEST}" original_manifest)
set(reduced_manifest "${original_manifest}")
string(REPLACE "        wirehair_decode;\n" ""
    reduced_manifest "${reduced_manifest}")
if(reduced_manifest STREQUAL original_manifest)
    message(FATAL_ERROR "Could not create reduced export manifest")
endif()
set(reduced_manifest_path "${WORK_DIR}/wirehair-reduced.map")
file(WRITE "${reduced_manifest_path}" "${reduced_manifest}")

execute_process(
    COMMAND "${CMAKE_COMMAND}"
        "-DLIBRARY=${LIBRARY}"
        "-DMANIFEST=${reduced_manifest_path}"
        "-DNM=${NM}"
        -P "${CHECKER}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err)
set(combined "${out}${err}")
if(result EQUAL 0)
    message(FATAL_ERROR
        "Reduced allowlist unexpectedly accepted the full export table\n"
        "${combined}")
endif()
if(NOT combined MATCHES "ELF dynamic export mismatch" OR
   NOT combined MATCHES "Extra: wirehair_decode@@")
    message(FATAL_ERROR
        "Export mismatch failed without the exact-table diagnostic\n"
        "${combined}")
endif()

message(STATUS "Verified 19 parser fixtures and fail-closed allowlist/export mismatch")
