cmake_minimum_required(VERSION 3.15)

foreach(required IN ITEMS CHECKER LIBRARY MANIFEST NM WORK_DIR)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(READ "${MANIFEST}" original_manifest)
if(NOT original_manifest MATCHES
   "^([A-Za-z_][A-Za-z0-9_.]*)[ \t]*\\{")
    message(FATAL_ERROR "Could not identify the export-manifest ABI version")
endif()
set(abi_version "${CMAKE_MATCH_1}")

# Capture the real callable export rows once, remove whichever spelling the
# host nm uses for the version namespace, and inject deterministic GNU/LLVM
# marker fixtures through tiny generated nm stand-ins.  This exercises both
# dialects even when CI happens to provide only one implementation.
execute_process(
    COMMAND "${NM}" -D --defined-only --extern-only --format=posix "${LIBRARY}"
    RESULT_VARIABLE nm_result
    OUTPUT_VARIABLE nm_output
    ERROR_VARIABLE nm_error)
if(NOT nm_result EQUAL 0)
    message(FATAL_ERROR
        "Could not capture nm output for export-checker self-test (${nm_result})\n"
        "stdout:\n${nm_output}\nstderr:\n${nm_error}")
endif()
string(REPLACE "\r\n" "\n" nm_output "${nm_output}")
string(REPLACE "\r" "\n" nm_output "${nm_output}")
string(REPLACE "\n" ";" nm_lines "${nm_output}")
set(callable_nm_output "")
foreach(line IN LISTS nm_lines)
    string(STRIP "${line}" line)
    if(line STREQUAL "")
        continue()
    endif()
    if(line MATCHES "^([^ \t]+)[ \t]+")
        set(nm_symbol "${CMAKE_MATCH_1}")
        if(nm_symbol STREQUAL "${abi_version}" OR
           nm_symbol STREQUAL "${abi_version}@@${abi_version}")
            continue()
        endif()
    endif()
    string(APPEND callable_nm_output "${line}\n")
endforeach()
if(callable_nm_output STREQUAL "")
    message(FATAL_ERROR "Captured nm output has no callable exports")
endif()

find_program(CHMOD_EXECUTABLE chmod)
if(NOT CHMOD_EXECUTABLE)
    message(FATAL_ERROR "chmod is required for generated nm self-test fixtures")
endif()

function(write_fake_nm basename contents output_var)
    set(fake_nm "${WORK_DIR}/${basename}")
    file(WRITE "${fake_nm}.out" "${contents}")
    file(WRITE "${fake_nm}"
        "#!/bin/sh\n"
        "exec \"${CMAKE_COMMAND}\" -E cat \"$0.out\"\n")
    execute_process(
        COMMAND "${CHMOD_EXECUTABLE}" 700 "${fake_nm}"
        RESULT_VARIABLE chmod_result
        ERROR_VARIABLE chmod_error)
    if(NOT chmod_result EQUAL 0)
        message(FATAL_ERROR
            "Could not make ${basename} executable (${chmod_result}): "
            "${chmod_error}")
    endif()
    set(${output_var} "${fake_nm}" PARENT_SCOPE)
endfunction()

function(expect_checker_pass label fake_nm)
    execute_process(
        COMMAND "${CMAKE_COMMAND}"
            "-DLIBRARY=${LIBRARY}"
            "-DMANIFEST=${MANIFEST}"
            "-DNM=${fake_nm}"
            -P "${CHECKER}"
        RESULT_VARIABLE fixture_result
        OUTPUT_VARIABLE fixture_out
        ERROR_VARIABLE fixture_err)
    if(NOT fixture_result EQUAL 0)
        message(FATAL_ERROR
            "${label} marker fixture was rejected (${fixture_result})\n"
            "stdout:\n${fixture_out}\nstderr:\n${fixture_err}")
    endif()
endfunction()

function(expect_checker_failure label fake_nm expected_extra)
    execute_process(
        COMMAND "${CMAKE_COMMAND}"
            "-DLIBRARY=${LIBRARY}"
            "-DMANIFEST=${MANIFEST}"
            "-DNM=${fake_nm}"
            -P "${CHECKER}"
        RESULT_VARIABLE fixture_result
        OUTPUT_VARIABLE fixture_out
        ERROR_VARIABLE fixture_err)
    set(fixture_combined "${fixture_out}${fixture_err}")
    if(fixture_result EQUAL 0)
        message(FATAL_ERROR "${label} marker fixture unexpectedly passed")
    endif()
    string(FIND "${fixture_combined}"
        "ELF dynamic export mismatch" mismatch_index)
    string(FIND "${fixture_combined}"
        "Extra: ${expected_extra}" extra_index)
    if(mismatch_index EQUAL -1 OR extra_index EQUAL -1)
        message(FATAL_ERROR
            "${label} marker fixture failed without the exact-table diagnostic\n"
            "${fixture_combined}")
    endif()
endfunction()

write_fake_nm(fake-nm-gnu
    "${abi_version} A 0\n${callable_nm_output}" fake_nm_gnu)
expect_checker_pass("GNU" "${fake_nm_gnu}")
write_fake_nm(fake-nm-llvm
    "${abi_version}@@${abi_version} A 0 0\n${callable_nm_output}"
    fake_nm_llvm)
expect_checker_pass("LLVM" "${fake_nm_llvm}")

write_fake_nm(fake-nm-gnu-wrong-type
    "${abi_version} T 0\n${callable_nm_output}" fake_nm_gnu_wrong_type)
expect_checker_failure(
    "GNU wrong-type" "${fake_nm_gnu_wrong_type}" "${abi_version}")
write_fake_nm(fake-nm-gnu-nonzero
    "${abi_version} A 1\n${callable_nm_output}" fake_nm_gnu_nonzero)
expect_checker_failure(
    "GNU nonzero" "${fake_nm_gnu_nonzero}" "${abi_version}")
write_fake_nm(fake-nm-llvm-wrong-type
    "${abi_version}@@${abi_version} T 0 0\n${callable_nm_output}"
    fake_nm_llvm_wrong_type)
expect_checker_failure(
    "LLVM wrong-type" "${fake_nm_llvm_wrong_type}"
    "${abi_version}@@${abi_version}")
write_fake_nm(fake-nm-llvm-nonzero
    "${abi_version}@@${abi_version} A 1 0\n${callable_nm_output}"
    fake_nm_llvm_nonzero)
expect_checker_failure(
    "LLVM nonzero" "${fake_nm_llvm_nonzero}"
    "${abi_version}@@${abi_version}")

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

message(STATUS
    "Verified GNU/LLVM version markers and fail-closed export mismatches")
