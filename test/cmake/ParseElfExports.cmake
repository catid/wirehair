# Shared by the exact-table checker and its synthetic parser tests. This only
# normalizes linker metadata; it never changes a callable symbol's spelling.
function(wirehair_parse_elf_exports nm_output abi_version output_variable)
    set(actual_exports "")
    set(marker_seen FALSE)
    string(REPLACE "\r\n" "\n" nm_output "${nm_output}")
    string(REPLACE "\r" "\n" nm_output "${nm_output}")
    string(REPLACE "\n" ";" nm_lines "${nm_output}")
    foreach(line IN LISTS nm_lines)
        string(STRIP "${line}" line)
        if(line STREQUAL "")
            continue()
        endif()
        if(NOT line MATCHES "^([^ \t]+)[ \t]+[A-Za-z?][ \t]+")
            message(FATAL_ERROR "Unrecognized nm output line: '${line}'")
        endif()
        set(symbol "${CMAKE_MATCH_1}")
        # GNU nm leaves GNU ld's ABS version node undecorated; LLVM nm prints
        # its default version too. Mold can omit this synthetic row altogether.
        # Require the exact zero-valued/zero-sized absolute marker, not merely
        # its name. A callable lookalike or duplicate must not disappear.
        if(symbol STREQUAL "${abi_version}" OR
           symbol STREQUAL "${abi_version}@@${abi_version}")
            if(marker_seen OR NOT line MATCHES "^[^ \t]+[ \t]+A[ \t]+0+([ \t]+0+)?$")
                message(FATAL_ERROR "Invalid or duplicate ELF version marker: '${line}'")
            endif()
            set(marker_seen TRUE)
            continue()
        endif()
        list(APPEND actual_exports "${symbol}")
    endforeach()
    set(${output_variable} "${actual_exports}" PARENT_SCOPE)
endfunction()
