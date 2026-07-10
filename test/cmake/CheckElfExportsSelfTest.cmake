cmake_minimum_required(VERSION 3.15)

foreach(required IN ITEMS CHECKER LIBRARY MANIFEST NM WORK_DIR)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
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

message(STATUS "Verified that an allowlist/export mismatch fails closed")
