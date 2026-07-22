cmake_minimum_required(VERSION 3.15)

foreach(required IN ITEMS LIBRARY MANIFEST NM)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

if(NOT EXISTS "${LIBRARY}")
    message(FATAL_ERROR "ELF shared library does not exist: ${LIBRARY}")
endif()
if(NOT EXISTS "${MANIFEST}")
    message(FATAL_ERROR "ELF export manifest does not exist: ${MANIFEST}")
endif()

# Parse the deliberately small linker-map grammar strictly.  Keeping this
# parser independent from the linker's parser catches malformed or accidentally
# broadened manifests before they can become ABI policy.
file(STRINGS "${MANIFEST}" manifest_lines)
set(state version)
set(abi_version "")
set(public_symbols "")
foreach(line IN LISTS manifest_lines)
    string(STRIP "${line}" line)
    if(line STREQUAL "")
        continue()
    endif()

    if(state STREQUAL version)
        if(NOT line MATCHES "^([A-Za-z_][A-Za-z0-9_.]*)[ \t]*\\{$")
            message(FATAL_ERROR
                "Malformed export manifest version line: '${line}'")
        endif()
        set(abi_version "${CMAKE_MATCH_1}")
        set(state global)
    elseif(state STREQUAL global)
        if(NOT line STREQUAL "global:")
            message(FATAL_ERROR
                "Malformed export manifest global section: '${line}'")
        endif()
        set(state symbols)
    elseif(state STREQUAL symbols)
        if(line STREQUAL "local:")
            if(NOT public_symbols)
                message(FATAL_ERROR "Export manifest has no public symbols")
            endif()
            set(state wildcard)
        elseif(line MATCHES "^(wirehair_[A-Za-z0-9_]+);$")
            set(symbol "${CMAKE_MATCH_1}")
            list(FIND public_symbols "${symbol}" duplicate_index)
            if(NOT duplicate_index EQUAL -1)
                message(FATAL_ERROR
                    "Duplicate public symbol in export manifest: ${symbol}")
            endif()
            list(APPEND public_symbols "${symbol}")
        else()
            message(FATAL_ERROR
                "Malformed public symbol in export manifest: '${line}'")
        endif()
    elseif(state STREQUAL wildcard)
        if(NOT line STREQUAL "*;")
            message(FATAL_ERROR
                "Export manifest must hide every unlisted symbol; found '${line}'")
        endif()
        set(state close)
    elseif(state STREQUAL close)
        if(NOT line STREQUAL "};")
            message(FATAL_ERROR
                "Malformed export manifest closing line: '${line}'")
        endif()
        set(state done)
    else()
        message(FATAL_ERROR
            "Unexpected content after export manifest: '${line}'")
    endif()
endforeach()
if(NOT state STREQUAL done)
    message(FATAL_ERROR
        "Truncated export manifest; parser stopped in '${state}' state")
endif()

execute_process(
    COMMAND "${NM}" -D --defined-only --extern-only --format=posix "${LIBRARY}"
    RESULT_VARIABLE nm_result
    OUTPUT_VARIABLE nm_output
    ERROR_VARIABLE nm_error)
if(NOT nm_result EQUAL 0)
    message(FATAL_ERROR
        "nm failed for ${LIBRARY} (${nm_result})\n"
        "stdout:\n${nm_output}\nstderr:\n${nm_error}")
endif()

set(actual_exports "")
set(version_marker_count 0)
string(REPLACE "\r\n" "\n" nm_output "${nm_output}")
string(REPLACE "\r" "\n" nm_output "${nm_output}")
string(REPLACE "\n" ";" nm_lines "${nm_output}")
foreach(line IN LISTS nm_lines)
    string(STRIP "${line}" line)
    if(line STREQUAL "")
        continue()
    endif()
    if(NOT line MATCHES
       "^([^ \t]+)[ \t]+([A-Za-z?])[ \t]+([^ \t]+)([ \t]+([^ \t]+))?$")
        message(FATAL_ERROR "Unrecognized nm output line: '${line}'")
    endif()
    set(nm_symbol "${CMAKE_MATCH_1}")
    set(nm_type "${CMAKE_MATCH_2}")
    set(nm_value "${CMAKE_MATCH_3}")
    set(nm_size "${CMAKE_MATCH_5}")

    # GNU nm prints the absolute version namespace without a suffix, while
    # LLVM nm renders that same linker-metadata row with its default-version
    # suffix.  Ignore either spelling only when the rest of the record proves
    # it is the zero-valued absolute marker.  A callable/data symbol or a
    # nonzero lookalike must remain in the exact export comparison.
    if((nm_symbol STREQUAL "${abi_version}" OR
        nm_symbol STREQUAL "${abi_version}@@${abi_version}") AND
       nm_type STREQUAL "A" AND
       nm_value MATCHES "^0+$" AND
       (nm_size STREQUAL "" OR nm_size MATCHES "^0+$"))
        math(EXPR version_marker_count "${version_marker_count} + 1")
        if(version_marker_count GREATER 1)
            message(FATAL_ERROR
                "nm reported duplicate ELF version namespace markers")
        endif()
        continue()
    endif()
    list(APPEND actual_exports "${nm_symbol}")
endforeach()

set(expected_exports "")
foreach(symbol IN LISTS public_symbols)
    list(APPEND expected_exports "${symbol}@@${abi_version}")
endforeach()
list(SORT expected_exports)
list(SORT actual_exports)

if(NOT actual_exports STREQUAL expected_exports)
    set(missing "")
    foreach(symbol IN LISTS expected_exports)
        if(NOT symbol IN_LIST actual_exports)
            list(APPEND missing "${symbol}")
        endif()
    endforeach()
    set(extra "")
    foreach(symbol IN LISTS actual_exports)
        if(NOT symbol IN_LIST expected_exports)
            list(APPEND extra "${symbol}")
        endif()
    endforeach()
    message(FATAL_ERROR
        "ELF dynamic export mismatch for ${LIBRARY}\n"
        "ABI version: ${abi_version}\n"
        "Missing: ${missing}\nExtra: ${extra}\n"
        "Complete nm output:\n${nm_output}")
endif()

list(LENGTH public_symbols public_symbol_count)
message(STATUS
    "Verified ${public_symbol_count} public ELF symbols at ${abi_version}")
