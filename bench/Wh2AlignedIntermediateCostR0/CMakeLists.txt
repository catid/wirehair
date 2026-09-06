cmake_minimum_required(VERSION 3.16)
project(Wh2AlignedIntermediateCostR0 LANGUAGES CXX)

option(WH2_ALIGNMENT_COST_SANITIZERS "Sanitize the new driver's neutral mode only" OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Frozen single configuration" FORCE)
endif()
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
   NOT CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR "Use native GNU C++ and a single Release configuration")
endif()
if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE)
    message(FATAL_ERROR "LTO is forbidden")
endif()
set(_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_EXE_LINKER_FLAGS}")
if(_flags MATCHES "WH_ALIGN64|WH_HUGEPAGE|-flto|-fprofile|WIREHAIR_TESTING|WIREHAIR_V2_ENABLE|WH_SEED_KNOBS")
    message(FATAL_ERROR "Do not inject allocator, LTO, profiling, hook, or equation flags")
endif()
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside_source)
if(_inside_source EQUAL 0)
    message(FATAL_ERROR "Use a fresh cost build directory outside the repository")
endif()
set(_base "199fe0b7916869ed28ee1090c36f14b2bcf82b9d")
find_package(Git REQUIRED)
function(_require_regular_if_exists path)
    if(IS_SYMLINK "${path}" OR IS_DIRECTORY "${path}")
        message(FATAL_ERROR "Invalid regular-file path: ${path}")
    endif()
    if(EXISTS "${path}")
        execute_process(COMMAND /usr/bin/test -f "${path}"
            RESULT_VARIABLE _regular TIMEOUT 2)
        if(NOT _regular STREQUAL "0")
            message(FATAL_ERROR "Non-regular file: ${path}")
        endif()
    endif()
endfunction()
set(_prototype_pins "")
foreach(_relative IN ITEMS bench/Wh2AlignedIntermediateR0.h
        bench/Wh2AlignedIntermediateR0/CMakeLists.txt)
    if(IS_SYMLINK "${_root}/${_relative}" OR IS_DIRECTORY "${_root}/${_relative}")
        message(FATAL_ERROR "Original mirror builder/allocator is not regular")
    endif()
    execute_process(COMMAND /usr/bin/test -f "${_root}/${_relative}"
        RESULT_VARIABLE _regular TIMEOUT 2)
    if(NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Missing original mirror builder/allocator")
    endif()
    file(SHA256 "${_root}/${_relative}" _sha)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${_root}" cat-file blob
        "${_base}:${_relative}" OUTPUT_VARIABLE _committed
        RESULT_VARIABLE _read_commit TIMEOUT 10)
    string(SHA256 _committed_sha "${_committed}")
    if(NOT _read_commit STREQUAL "0" OR NOT _sha STREQUAL _committed_sha)
        message(FATAL_ERROR "Original mirror builder/allocator changed")
    endif()
    string(APPEND _prototype_pins "input ${_sha} ${_root}/${_relative}\n")
endforeach()

# Configure, but never build, the original one-variant targets.  This reuses its
# checked full regular-file mirrors and exactly two storage overlays unchanged.
set(_p "${_binary}/ordinary-mirror")
set(_a "${_binary}/aligned-mirror")
foreach(_variant IN ITEMS p a)
    if(_variant STREQUAL "a")
        set(_candidate ON)
    else()
        set(_candidate OFF)
    endif()
    execute_process(COMMAND "${CMAKE_COMMAND}"
        -S "${_root}/bench/Wh2AlignedIntermediateR0" -B "${_${_variant}}"
        -G Ninja -DCMAKE_BUILD_TYPE=Release
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
        -DCMAKE_CXX_FLAGS= "-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG"
        "-DWH2_ALIGNMENT_CANDIDATE=${_candidate}"
        -DWH2_ALIGNMENT_NEUTRAL_TESTS=OFF -DWH2_ALIGNMENT_SANITIZERS=OFF
        RESULT_VARIABLE _configured TIMEOUT 30)
    if(NOT _configured STREQUAL "0")
        message(FATAL_ERROR "Cannot configure the authenticated ${_variant} mirror")
    endif()
endforeach()

set(_support "${_binary}/support")
foreach(_directory IN ITEMS "${_support}" "${_support}/bench")
    if(IS_SYMLINK "${_directory}" OR
       (EXISTS "${_directory}" AND NOT IS_DIRECTORY "${_directory}"))
        message(FATAL_ERROR "Invalid support directory")
    endif()
    file(MAKE_DIRECTORY "${_directory}")
endforeach()
set(_manifest_path "${_binary}/alignment-cost-manifest.txt")
_require_regular_if_exists("${_manifest_path}")
set(_inputs
    bench/Wh2AlignedIntermediateCostR0.cpp
    bench/Wh2AlignedIntermediateCostR0Bridge.h
    bench/Wh2AlignedIntermediateCostR0Inspect.inc
    bench/Wh2FrozenTrace.cpp bench/Wh2FrozenTrace.h
    bench/Wh2RdpruTargetIdentityV2.cpp bench/Wh2RdpruTargetIdentityV2.h
    bench/Wh2PublicBorrowedTargetIdentity.cpp bench/Wh2PublicBorrowedTargetIdentity.h)
set(_manifest "wirehair.wh2.aligned-intermediate-cost-r0.build.v1\nbase=${_base}\nsanitizers=${WH2_ALIGNMENT_COST_SANITIZERS}\n${_prototype_pins}")
foreach(_relative IN LISTS _inputs)
    set(_source "${_root}/${_relative}")
    set(_copy "${_support}/${_relative}")
    if(IS_SYMLINK "${_source}" OR IS_DIRECTORY "${_source}")
        message(FATAL_ERROR "Support files must be regular: ${_relative}")
    endif()
    execute_process(COMMAND /usr/bin/test -f "${_source}"
        RESULT_VARIABLE _regular TIMEOUT 2)
    if(NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Missing/non-regular support source: ${_relative}")
    endif()
    _require_regular_if_exists("${_copy}")
    configure_file("${_source}" "${_copy}" COPYONLY)
    file(SHA256 "${_source}" _sha)
    file(SHA256 "${_copy}" _copy_sha)
    if(NOT _sha STREQUAL _copy_sha)
        message(FATAL_ERROR "Source changed during support copy: ${_relative}")
    endif()
    if(NOT _relative MATCHES "AlignedIntermediateCostR0")
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${_root}" cat-file blob
            "${_base}:${_relative}" OUTPUT_VARIABLE _committed
            RESULT_VARIABLE _read_commit TIMEOUT 10)
        string(SHA256 _committed_sha "${_committed}")
        if(NOT _read_commit STREQUAL "0" OR NOT _sha STREQUAL _committed_sha)
            message(FATAL_ERROR "Helper differs from the frozen source: ${_relative}")
        endif()
    endif()
    string(APPEND _manifest "input ${_sha} ${_source}\nfile ${_copy_sha} ${_copy}\n")
endforeach()

# Save the ordinary facade before appending anything.  Its quotes still resolve
# in the SAME codec directory; compile it for assembly comparison, never link it.
set(_reference "${_p}/source/codec/WirehairV2ProfileReference.cpp")
_require_regular_if_exists("${_reference}")
configure_file("${_p}/source/codec/WirehairV2Profile.cpp" "${_reference}" COPYONLY)
foreach(_variant IN ITEMS p a)
    set(_profile "${_${_variant}}/source/codec/WirehairV2Profile.cpp")
    file(APPEND "${_profile}" "\n#include \"Wh2AlignedIntermediateCostR0Inspect.inc\"\n")
    file(GLOB_RECURSE _files LIST_DIRECTORIES false "${_${_variant}}/source/*")
    list(SORT _files)
    foreach(_file IN LISTS _files)
        _require_regular_if_exists("${_file}")
        file(SHA256 "${_file}" _sha)
        string(APPEND _manifest "file ${_sha} ${_file}\n")
    endforeach()
    file(SHA256 "${_${_variant}}/alignment-source-manifest.txt" _sha)
    string(APPEND _manifest "preappend_manifest ${_sha} ${_${_variant}}/alignment-source-manifest.txt\n")
endforeach()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Threads REQUIRED)
add_library(wh2_alignment_cost_options INTERFACE)
target_compile_features(wh2_alignment_cost_options INTERFACE cxx_std_11)
target_compile_options(wh2_alignment_cost_options INTERFACE
    -O3 -march=native -Wall -Wextra -Wpedantic -Werror -fno-lto)
target_compile_definitions(wh2_alignment_cost_options INTERFACE WIREHAIR_STATIC=1)
target_link_options(wh2_alignment_cost_options INTERFACE -fno-lto)
target_link_libraries(wh2_alignment_cost_options INTERFACE Threads::Threads m)
if(WH2_ALIGNMENT_COST_SANITIZERS)
    target_compile_options(wh2_alignment_cost_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(wh2_alignment_cost_options INTERFACE -fsanitize=address,undefined)
endif()

set(_common_names wirehair.cpp gf256.cpp WirehairCodec.cpp WirehairTools.cpp)
set(_v2_names WirehairV2Codec.cpp WirehairV2Peel.cpp WirehairV2Plan.cpp
    WirehairV2Policy.cpp WirehairV2Precode.cpp WirehairV2PrecodeDecode.cpp
    WirehairV2PrecodeEncode.cpp WirehairV2Profile.cpp WirehairV2Seeds.cpp WirehairV2Solve.cpp)
set(_common_sources "")
foreach(_name IN LISTS _common_names)
    list(APPEND _common_sources "${_p}/source/${_name}")
endforeach()
add_library(wh2_alignment_cost_common STATIC ${_common_sources})
target_include_directories(wh2_alignment_cost_common PRIVATE
    "${_p}/source" "${_p}/source/include")
target_link_libraries(wh2_alignment_cost_common PUBLIC wh2_alignment_cost_options)
target_compile_definitions(wh2_alignment_cost_common PRIVATE WIREHAIR_BUILDING=1)
set_target_properties(wh2_alignment_cost_common PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(_exports decode decoder_create encode encoder_create encoder_create_profile
    encoder_create_profile_id encoder_create_profile_id_with_options
    encoder_create_profile_with_options encoder_create_with_options
    encoder_detach_input free profile_deserialize profile_serialize
    profile_validate recover result_string)
set(_renames wirehair_v2=wirehair_v2_aligned_intermediate_r0)
foreach(_name IN LISTS _exports)
    list(APPEND _renames "wirehair_v2_${_name}=wh2_aligned_r0_${_name}")
endforeach()
foreach(_variant IN ITEMS p a)
    set(_sources "")
    foreach(_name IN LISTS _v2_names)
        list(APPEND _sources "${_${_variant}}/source/codec/${_name}")
    endforeach()
    set(_target "wh2_alignment_cost_${_variant}")
    add_library(${_target} STATIC ${_sources})
    target_include_directories(${_target} PRIVATE "${_${_variant}}/source"
        "${_${_variant}}/source/codec" "${_${_variant}}/source/include"
        "${_${_variant}}/source/bench" "${_support}/bench")
    target_link_libraries(${_target} PUBLIC wh2_alignment_cost_options)
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_definitions(${_target} PRIVATE WIREHAIR_BUILDING=1
        "WH2_ALIGNMENT_INSPECT_FN=wh2_aligned_cost_${_variant}_inspect"
        "WH2_ALIGNMENT_COLUMNS_FN=wh2_aligned_cost_${_variant}_packet_columns"
        "WH2_ALIGNMENT_SCALAR_FN=wh2_aligned_cost_${_variant}_scalar_packet")
    if(_variant STREQUAL "a")
        target_compile_definitions(${_target} PRIVATE ${_renames} WH2_ALIGNMENT_INSPECT_ALIGNED=1)
    else()
        target_compile_definitions(${_target} PRIVATE WH2_ALIGNMENT_INSPECT_ALIGNED=0)
    endif()
endforeach()

add_library(wh2_alignment_cost_profile_reference OBJECT "${_reference}")
target_include_directories(wh2_alignment_cost_profile_reference PRIVATE
    "${_p}/source" "${_p}/source/codec" "${_p}/source/include")
target_compile_definitions(wh2_alignment_cost_profile_reference PRIVATE WIREHAIR_BUILDING=1)
target_link_libraries(wh2_alignment_cost_profile_reference PRIVATE wh2_alignment_cost_options)
set_target_properties(wh2_alignment_cost_profile_reference PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(wh2_aligned_intermediate_cost_r0
    "${_support}/bench/Wh2AlignedIntermediateCostR0.cpp"
    "${_support}/bench/Wh2FrozenTrace.cpp"
    "${_support}/bench/Wh2RdpruTargetIdentityV2.cpp"
    "${_support}/bench/Wh2PublicBorrowedTargetIdentity.cpp")
target_include_directories(wh2_aligned_intermediate_cost_r0 PRIVATE
    "${_support}/bench" "${_p}/source" "${_p}/source/include")
target_compile_definitions(wh2_aligned_intermediate_cost_r0 PRIVATE
    WH2_ALIGNMENT_COST_SANITIZERS=$<BOOL:${WH2_ALIGNMENT_COST_SANITIZERS}>)
target_compile_options(wh2_aligned_intermediate_cost_r0 PRIVATE -fPIC)
target_link_libraries(wh2_aligned_intermediate_cost_r0 PRIVATE
    wh2_alignment_cost_p wh2_alignment_cost_a wh2_alignment_cost_common)
target_link_options(wh2_aligned_intermediate_cost_r0 PRIVATE
    "LINKER:-Map,${_binary}/alignment-cost.map")
add_dependencies(wh2_aligned_intermediate_cost_r0 wh2_alignment_cost_profile_reference)

file(SHA256 "${CMAKE_CURRENT_LIST_FILE}" _sha)
string(APPEND _manifest "input ${_sha} ${CMAKE_CURRENT_LIST_FILE}\nrenames=${_renames}\ncompiler=${CMAKE_CXX_COMPILER}\ncompiler_version=${CMAKE_CXX_COMPILER_VERSION}\n")
file(WRITE "${_manifest_path}" "${_manifest}")
