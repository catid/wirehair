if(NOT DEFINED NM OR NM STREQUAL "")
    message(FATAL_ERROR
        "public-borrowed r1 symbol audit requires CMAKE_NM")
endif()
if(NOT DEFINED PUBLIC_BORROWED_R1 OR
   NOT EXISTS "${PUBLIC_BORROWED_R1}")
    message(FATAL_ERROR
        "public-borrowed r1 symbol audit is missing its executable")
endif()

# This gate is registered only for native Linux shared-library builds.  The
# Required undefined references prove that the timed arms and untimed closure
# oracles cross only the public DSO boundary. A static or private
# implementation-body linkage cannot satisfy this contract.
execute_process(
    COMMAND "${NM}" -D -C "${PUBLIC_BORROWED_R1}"
    RESULT_VARIABLE _nm_result
    OUTPUT_VARIABLE _symbols
    ERROR_VARIABLE _nm_error)
if(NOT _nm_result EQUAL 0)
    message(FATAL_ERROR
        "nm -D failed for public-borrowed r1 executable: ${_nm_error}")
endif()
string(STRIP "${_symbols}" _symbols_stripped)
string(TOLOWER "${_symbols_stripped}" _symbols_lower)
if(_symbols_stripped STREQUAL "" OR _symbols_lower MATCHES "no symbols")
    message(FATAL_ERROR
        "nm -D returned no auditable symbols for public-borrowed r1 executable")
endif()

string(REPLACE "\r\n" "\n" _symbol_lines "${_symbols}")
string(REPLACE "\r" "\n" _symbol_lines "${_symbol_lines}")
string(REPLACE "\n" ";" _symbol_lines "${_symbol_lines}")

set(_required_public_imports
    wirehair_init_
    wirehair_v2_encoder_create_with_options
    wirehair_v2_encode
    wirehair_v2_profile_validate
    wirehair_v2_profile_deserialize
    wirehair_v2_decoder_create
    wirehair_v2_decode
    wirehair_v2_recover
    wirehair_v2_free
    wirehair_encoder_create_ex
    wirehair_encode
    wirehair_decoder_create
    wirehair_decode
    wirehair_recover
    wirehair_free)
foreach(_required IN LISTS _required_public_imports)
    set(_found FALSE)
    foreach(_line IN LISTS _symbol_lines)
        string(STRIP "${_line}" _line)
        if(_line MATCHES
           "^U[ \t]+_?${_required}(@@?[A-Za-z0-9_.-]+)?$")
            set(_found TRUE)
            break()
        endif()
    endforeach()
    if(NOT _found)
        message(FATAL_ERROR
            "public-borrowed r1 executable lacks dynamic undefined public "
            "API reference ${_required}")
    endif()
endforeach()

# The benchmark must not reach around the C ABI into codec implementation
# classes, private benchmark initializers, or test-hook state.  Inspect the
# complete symbol table for this negative gate: private archive linkage would
# not normally be copied into the dynamic table audited above.
execute_process(
    COMMAND "${NM}" -C "${PUBLIC_BORROWED_R1}"
    RESULT_VARIABLE _all_nm_result
    OUTPUT_VARIABLE _all_symbols
    ERROR_VARIABLE _all_nm_error)
if(NOT _all_nm_result EQUAL 0)
    message(FATAL_ERROR
        "nm failed for public-borrowed r1 executable: ${_all_nm_error}")
endif()
string(STRIP "${_all_symbols}" _all_symbols_stripped)
string(TOLOWER "${_all_symbols_stripped}" _all_symbols_lower)
if(_all_symbols_stripped STREQUAL "" OR
   _all_symbols_lower MATCHES "no symbols")
    message(FATAL_ERROR
        "nm returned no full symbol table for public-borrowed r1 executable")
endif()

# Check both demangled and common raw-mangled spellings so a demangler
# mismatch fails closed rather than weakening this audit.
set(_forbidden_linkage_tokens
    "ForTesting"
    "ForBenchmark"
    "TestHook"
    "wirehair_v2::"
    "wirehair::Codec"
    "WirehairCodec::"
    "_ZN11wirehair_v2"
    "_ZN8wirehair"
    "InitializeForValidatedSystem")
foreach(_forbidden IN LISTS _forbidden_linkage_tokens)
    string(FIND "${_all_symbols}" "${_forbidden}" _forbidden_offset)
    if(NOT _forbidden_offset EQUAL -1)
        message(FATAL_ERROR
            "public-borrowed r1 executable contains internal/test-hook "
            "dynamic linkage token ${_forbidden}")
    endif()
endforeach()

message(STATUS
    "public-borrowed r1 dynamic public-API symbol audit passed")
