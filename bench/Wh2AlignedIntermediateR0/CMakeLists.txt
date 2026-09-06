cmake_minimum_required(VERSION 3.16)
project(Wh2AlignedIntermediateR0 LANGUAGES CXX)

# Standalone, non-installed, one-variant-per-build neutral prototype.
option(WH2_ALIGNMENT_CANDIDATE "Use the frozen aligned-storage overlay" ON)
option(WH2_ALIGNMENT_NEUTRAL_TESTS "Build with consistent neutral test hooks" ON)
option(WH2_ALIGNMENT_SANITIZERS "Enable ASan and UBSan for neutral tests" OFF)

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "This frozen prototype requires native Linux/GNU C++")
endif()
if(CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR "Use a single-configuration build")
endif()
if(WH2_ALIGNMENT_SANITIZERS AND NOT WH2_ALIGNMENT_NEUTRAL_TESTS)
    message(FATAL_ERROR "Sanitizers are restricted to neutral tests")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" _config)
if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR
   CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_config})
    message(FATAL_ERROR "LTO is forbidden")
endif()
set(_external_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_config}} ${CMAKE_EXE_LINKER_FLAGS}")
if(_external_flags MATCHES "WH_ALIGN64|WH_HUGEPAGE|-flto|-fprofile|WIREHAIR_TESTING|WIREHAIR_V2_ENABLE_TEST_HOOKS")
    message(FATAL_ERROR "Do not inject allocator, profiling, LTO, or hook flags")
endif()

get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside_source)
if(_inside_source EQUAL 0)
    message(FATAL_ERROR "Use a fresh build directory outside the source tree")
endif()
set(_base "6442bb19ab17a014f6324c650bef89b982e58f63")
set(_mirror "${_binary}/source")
set(_identity "wirehair.wh2.aligned-intermediate-r0.prototype.v1\nbase=${_base}\nroot=${_root}\ncandidate=${WH2_ALIGNMENT_CANDIDATE}\nneutral=${WH2_ALIGNMENT_NEUTRAL_TESTS}\nsanitizers=${WH2_ALIGNMENT_SANITIZERS}\n")
set(_owner "${_binary}/alignment-source-owner.txt")
foreach(_evidence IN ITEMS "${_owner}" "${_binary}/alignment-source-manifest.txt")
    if(IS_SYMLINK "${_evidence}" OR IS_DIRECTORY "${_evidence}")
        message(FATAL_ERROR "Invalid mirror identity file: ${_evidence}")
    endif()
    if(EXISTS "${_evidence}")
        execute_process(COMMAND /usr/bin/test -f "${_evidence}"
            RESULT_VARIABLE _regular TIMEOUT 2)
        if(NOT _regular STREQUAL "0")
            message(FATAL_ERROR "Non-regular mirror identity file: ${_evidence}")
        endif()
    endif()
endforeach()
if(EXISTS "${_mirror}" AND NOT EXISTS "${_owner}")
    message(FATAL_ERROR "Refusing an unowned source mirror")
endif()
if(EXISTS "${_owner}")
    if(IS_SYMLINK "${_owner}" OR IS_DIRECTORY "${_owner}")
        message(FATAL_ERROR "Invalid mirror ownership file")
    endif()
    file(READ "${_owner}" _prior_identity)
    if(NOT _prior_identity STREQUAL _identity)
        message(FATAL_ERROR "Build identity changed: use a fresh build directory")
    endif()
endif()

set(_translation_units
    wirehair.cpp gf256.cpp WirehairCodec.cpp WirehairTools.cpp
    codec/WirehairV2Codec.cpp codec/WirehairV2Peel.cpp
    codec/WirehairV2Plan.cpp codec/WirehairV2Policy.cpp
    codec/WirehairV2Precode.cpp codec/WirehairV2PrecodeDecode.cpp
    codec/WirehairV2PrecodeEncode.cpp codec/WirehairV2Profile.cpp
    codec/WirehairV2Seeds.cpp codec/WirehairV2Solve.cpp)
set(_headers
    gf256.h WirehairCodec.h WirehairCompiler.h WirehairEnvironment.h
    WirehairHeavy.h WirehairTools.h WirehairDenseFixups.inc WirehairPeelFixups.inc
    include/wirehair/wirehair.h include/wirehair/wirehair.hpp
    codec/WirehairV2Codec.h codec/WirehairV2Peel.h codec/WirehairV2Plan.h
    codec/WirehairV2Policy.h codec/WirehairV2Precode.h
    codec/WirehairV2PrecodeDecode.h codec/WirehairV2PrecodeEncode.h
    codec/WirehairV2RawSeed.h codec/WirehairV2Seeds.h codec/WirehairV2Solve.h)
set(_production ${_translation_units} ${_headers})
find_package(Git REQUIRED)

foreach(_directory IN ITEMS "${_mirror}" "${_mirror}/codec" "${_mirror}/include"
        "${_mirror}/include/wirehair" "${_mirror}/bench")
    if(IS_SYMLINK "${_directory}" OR
       (EXISTS "${_directory}" AND NOT IS_DIRECTORY "${_directory}"))
        message(FATAL_ERROR "Invalid mirror directory: ${_directory}")
    endif()
    file(MAKE_DIRECTORY "${_directory}")
endforeach()
file(WRITE "${_owner}" "${_identity}")
set(_inputs ${_production} bench/Wh2AlignedIntermediateR0.h)
if(WH2_ALIGNMENT_NEUTRAL_TESTS)
    list(APPEND _inputs bench/Wh2AlignedIntermediateR0Test.cpp)
endif()
list(SORT _inputs)
set(_manifest "${_identity}")
foreach(_relative IN LISTS _inputs)
    set(_input "${_root}/${_relative}")
    set(_output "${_mirror}/${_relative}")
    if(NOT EXISTS "${_input}" OR IS_SYMLINK "${_input}" OR IS_DIRECTORY "${_input}" OR
       IS_SYMLINK "${_output}" OR IS_DIRECTORY "${_output}")
        message(FATAL_ERROR "Mirror inputs and outputs must be regular files: ${_relative}")
    endif()
    execute_process(COMMAND /usr/bin/test -f "${_input}"
        RESULT_VARIABLE _regular TIMEOUT 2)
    if(NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Non-regular source file: ${_relative}")
    endif()
    if(EXISTS "${_output}")
        execute_process(COMMAND /usr/bin/test -f "${_output}"
            RESULT_VARIABLE _regular TIMEOUT 2)
        if(NOT _regular STREQUAL "0")
            message(FATAL_ERROR "Non-regular mirror file: ${_relative}")
        endif()
    endif()
    configure_file("${_input}" "${_output}" COPYONLY)
    file(SHA256 "${_input}" _input_sha)
    file(SHA256 "${_output}" _output_sha)
    if(NOT _input_sha STREQUAL _output_sha)
        message(FATAL_ERROR "Source changed while copying: ${_relative}")
    endif()
    if(_relative IN_LIST _production)
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${_root}" cat-file blob
            "${_base}:${_relative}" OUTPUT_VARIABLE _committed
            RESULT_VARIABLE _read_commit TIMEOUT 10)
        if(NOT _read_commit STREQUAL "0")
            message(FATAL_ERROR "Cannot read frozen source blob: ${_relative}")
        endif()
        string(SHA256 _committed_sha "${_committed}")
        if(NOT _output_sha STREQUAL _committed_sha)
            message(FATAL_ERROR "Copied source differs from frozen blob: ${_relative}")
        endif()
    endif()
    string(APPEND _manifest "input ${_input_sha} ${_relative}\n")
endforeach()
file(SHA256 "${_mirror}/codec/WirehairV2PrecodeEncode.h" _header_sha)
file(SHA256 "${_mirror}/codec/WirehairV2PrecodeEncode.cpp" _cpp_sha)
if(NOT _header_sha STREQUAL "a48d2e2ed770b6d4b0963361ba11e19ee1e59f5529c599cefbaaf01e6762f403" OR
   NOT _cpp_sha STREQUAL "b12cbd859ab4fa6be515395da265859c62a92325de504599371a21bfaa214218")
    message(FATAL_ERROR "Frozen overlay source hashes do not match")
endif()

function(_replace_once variable old new)
    string(FIND "${${variable}}" "${old}" _first)
    if(_first LESS 0)
        message(FATAL_ERROR "Frozen overlay anchor is absent")
    endif()
    string(LENGTH "${old}" _length)
    math(EXPR _after "${_first}+${_length}")
    string(SUBSTRING "${${variable}}" ${_after} -1 _tail)
    string(FIND "${_tail}" "${old}" _again)
    if(NOT _again EQUAL -1)
        message(FATAL_ERROR "Frozen overlay anchor is not unique")
    endif()
    string(REPLACE "${old}" "${new}" _replaced "${${variable}}")
    set(${variable} "${_replaced}" PARENT_SCOPE)
endfunction()

if(WH2_ALIGNMENT_CANDIDATE)
    file(READ "${_mirror}/codec/WirehairV2PrecodeEncode.h" _header)
    _replace_once(_header "#include <vector>" "#include <vector>\n#include \"../bench/Wh2AlignedIntermediateR0.h\"")
    _replace_once(_header "    std::vector<uint8_t> SolvedIntermediateStorage;"
        "    wh2_aligned_intermediate_r0::ByteVector SolvedIntermediateStorage;")
    file(WRITE "${_mirror}/codec/WirehairV2PrecodeEncode.h" "${_header}")
    file(READ "${_mirror}/codec/WirehairV2PrecodeEncode.cpp" _cpp)
    _replace_once(_cpp "        next.SolvedIntermediateStorage.swap(intermediate_blocks);"
        "        next.SolvedIntermediateStorage.assign(\n            intermediate_blocks.begin(), intermediate_blocks.end());")
    _replace_once(_cpp "        next.UsesPacketContract = true;\n        next.Initialized = true;\n        Swap(next);\n        return Wirehair_Success;"
        "        next.UsesPacketContract = true;\n        next.Initialized = true;\n        Swap(next);\n        std::vector<uint8_t>().swap(intermediate_blocks);\n        return Wirehair_Success;")
    file(WRITE "${_mirror}/codec/WirehairV2PrecodeEncode.cpp" "${_cpp}")
endif()
foreach(_relative IN LISTS _inputs)
    file(SHA256 "${_mirror}/${_relative}" _sha)
    string(APPEND _manifest "mirror ${_sha} ${_relative}\n")
endforeach()
file(SHA256 "${CMAKE_CURRENT_LIST_FILE}" _cmake_sha)
string(APPEND _manifest "cmake ${_cmake_sha}\ncompiler ${CMAKE_CXX_COMPILER}\ncompiler_version ${CMAKE_CXX_COMPILER_VERSION}\n")
file(WRITE "${_binary}/alignment-source-manifest.txt" "${_manifest}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
find_package(Threads REQUIRED)
set(_sources "")
foreach(_relative IN LISTS _translation_units)
    list(APPEND _sources "${_mirror}/${_relative}")
endforeach()
add_library(wh2_aligned_intermediate_r0_library STATIC ${_sources})
target_compile_features(wh2_aligned_intermediate_r0_library PUBLIC cxx_std_11)
set_target_properties(wh2_aligned_intermediate_r0_library PROPERTIES
    POSITION_INDEPENDENT_CODE ON)
target_include_directories(wh2_aligned_intermediate_r0_library PUBLIC
    "${_mirror}" "${_mirror}/codec" "${_mirror}/bench" "${_mirror}/include")
target_compile_definitions(wh2_aligned_intermediate_r0_library
    PRIVATE WIREHAIR_BUILDING=1 PUBLIC WIREHAIR_STATIC=1)
target_compile_options(wh2_aligned_intermediate_r0_library PUBLIC
    -O3 -march=native -Wall -Wextra -Wpedantic -Werror -fno-lto)
target_link_options(wh2_aligned_intermediate_r0_library INTERFACE -fno-lto)
target_link_libraries(wh2_aligned_intermediate_r0_library PUBLIC Threads::Threads m)
if(WH2_ALIGNMENT_NEUTRAL_TESTS)
    target_compile_definitions(wh2_aligned_intermediate_r0_library PUBLIC
        WIREHAIR_TESTING=1 WIREHAIR_V2_ENABLE_TEST_HOOKS=1)
    enable_testing()
    add_executable(wh2_aligned_intermediate_r0_test
        "${_mirror}/bench/Wh2AlignedIntermediateR0Test.cpp")
    target_link_libraries(wh2_aligned_intermediate_r0_test PRIVATE
        wh2_aligned_intermediate_r0_library)
    target_compile_definitions(wh2_aligned_intermediate_r0_test PRIVATE
        WH2_ALIGNMENT_CANDIDATE=$<BOOL:${WH2_ALIGNMENT_CANDIDATE}>)
    add_test(NAME wh2_aligned_intermediate_r0_neutral
        COMMAND wh2_aligned_intermediate_r0_test)
    set_tests_properties(wh2_aligned_intermediate_r0_neutral PROPERTIES TIMEOUT 60)
endif()
if(WH2_ALIGNMENT_SANITIZERS)
    target_compile_options(wh2_aligned_intermediate_r0_library PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(wh2_aligned_intermediate_r0_library INTERFACE
        -fsanitize=address,undefined)
endif()
