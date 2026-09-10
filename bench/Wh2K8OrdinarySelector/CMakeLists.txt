cmake_minimum_required(VERSION 3.20)
project(Wh2K8OrdinarySelector LANGUAGES C CXX)
option(WH2_K8_ORDINARY_NEUTRAL "Build bounded correctness checks, never timing" OFF)

# One ordinary-K8 selector insertion, correctness-only first. Never edit/rebuild
# the qualified library or production source; this target does not time code.
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside)
if(_inside EQUAL 0 OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
   NOT CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "13.3.0")
    message(FATAL_ERROR "Use an external native GNU13.3/Linux build")
endif()

set(_source "${_root}/codec/WirehairV2Profile.cpp")
file(SHA256 "${_source}" _sha)
if(NOT _sha STREQUAL "975da8d892363d05de3ad4535ec79b4a386bf3f6fab377b114af393f82d96d15")
    message(FATAL_ERROR "Frozen profile source changed")
endif()
file(READ "${_source}" _candidate)
include("${CMAKE_CURRENT_SOURCE_DIR}/Overlay.cmake")
set(_generated "${_binary}/WirehairV2K8OrdinarySelector.cpp")
file(WRITE "${_generated}" "${_candidate}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_source}")
include("${CMAKE_CURRENT_SOURCE_DIR}/TestOverlay.cmake")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_library(k8_ordinary_baseline OBJECT "${_source}")
add_library(k8_ordinary_candidate OBJECT "${_generated}")
foreach(_target IN ITEMS k8_ordinary_baseline k8_ordinary_candidate)
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_definitions(${_target} PRIVATE WIREHAIR_BUILDING=1 NDEBUG)
    target_compile_options(${_target} PRIVATE -O3 -Wall -Wextra -Wpedantic -Werror)
    # Needed for the external TU's relative includes. Exact baseline object
    # reproduction below verifies that these include paths are neutral.
    target_include_directories(${_target} PRIVATE "${_root}/include" "${_root}/codec" "${_root}")
endforeach()

set(_qualified "/tmp/wh2-v2-k8-admission.YNN2XmAx")
set(_original "${_qualified}/native/CMakeFiles/wirehair_objects.dir/codec/WirehairV2Profile.cpp.o")
file(SHA256 "${_original}" _original_sha)
if(NOT _original_sha STREQUAL "229da2af98169c140774fff3bc970a23b6688bba1586ae49fae862554aa64071")
    message(FATAL_ERROR "Qualified original object changed")
endif()
add_custom_target(check_original ALL
    COMMAND "${CMAKE_COMMAND}" -E compare_files
        "$<TARGET_OBJECTS:k8_ordinary_baseline>" "${_original}"
    DEPENDS k8_ordinary_baseline
    COMMENT "Require byte-identical installed K8 baseline object"
    VERBATIM)

if(WH2_K8_ORDINARY_NEUTRAL)
    # Reuse the prior independent wire oracle/parity harness without editing it.
    file(SHA256 "${_root}/bench/Wh2ProfileClassification/ValidationTest.cpp" _metadata_sha)
    file(SHA256 "${_root}/bench/Wh2ProfileClassification/CheckParity.cmake" _parity_sha)
    if(NOT _metadata_sha STREQUAL "265d65c1f9a638a7251ebf0a1cd71eec2fe7b95bb063b89aecb67277d2c565b7" OR
       NOT _parity_sha STREQUAL "a7fa7dd8cdf57e038c57414cb875a8893c17f5130d7821f4971217481eed1efb")
        message(FATAL_ERROR "Reused neutral oracle changed")
    endif()
    enable_testing()
    find_package(Threads REQUIRED)
    foreach(_mode IN ITEMS native scalar asan)
        set(_base "${_qualified}/${_mode}")
        set(_defines WIREHAIR_STATIC=1)
        set(_flags -O3 -DNDEBUG)
        set(_link_flags "")
        if(_mode STREQUAL "native")
            set(_archive_sha "0d4a0a7ec3cfa4b8e0f08df75e2c12935b51f28653c8c953a32bb49a5658eff5")
            set(_object k8_ordinary_candidate)
        elseif(_mode STREQUAL "scalar")
            set(_archive_sha "bd7eb2eeb0eeaea15e11d131c6fba2921bf44cb1526dea8dff86171b496de867")
            list(APPEND _defines ANDROID=1)
            set(_object k8_ordinary_scalar)
        else()
            set(_archive_sha "a7043eb65eb864aefd210f51b360d4591a365c60710c3aec87baf0abac3d4f1a")
            # Match the archive's sanitizer instrumentation and ISA-dependent
            # private GF layout. This fresh Release target retains -O3; that is
            # not optimization-flag/object identity with its debug archive.
            # No timing uses ASAN.
            set(_flags -g -fsanitize=address,undefined -fno-omit-frame-pointer -march=native)
            set(_link_flags -fsanitize=address,undefined -no-pie)
            set(_object k8_ordinary_asan)
        endif()
        file(SHA256 "${_base}/libwirehair.a" _actual)
        if(NOT _actual STREQUAL _archive_sha)
            message(FATAL_ERROR "Qualified ${_mode} archive changed")
        endif()
        if(NOT _mode STREQUAL "native")
            add_library(${_object} OBJECT "${_generated}")
            set_target_properties(${_object} PROPERTIES POSITION_INDEPENDENT_CODE ON)
            target_include_directories(${_object} PRIVATE "${_root}/include" "${_root}/codec" "${_root}")
            target_compile_options(${_object} PRIVATE ${_flags} -Wall -Wextra -Wpedantic -Werror)
            target_compile_definitions(${_object} PRIVATE ${_defines} WIREHAIR_BUILDING=1)
        endif()
        # This legacy fault suite requires the facade-only allocation hook.
        # Keep its objects separate from production-style candidate evidence.
        foreach(_arm IN ITEMS baseline candidate)
            set(_fault_object k8_ordinary_fault_${_mode}_${_arm})
            if(_arm STREQUAL "baseline")
                set(_fault_source "${_source}")
            else()
                set(_fault_source "${_generated}")
            endif()
            add_library(${_fault_object} OBJECT "${_fault_source}")
            set_target_properties(${_fault_object} PROPERTIES POSITION_INDEPENDENT_CODE ON)
            target_include_directories(${_fault_object} PRIVATE "${_root}/include" "${_root}/codec" "${_root}")
            target_compile_options(${_fault_object} PRIVATE ${_flags} -Wall -Wextra -Wpedantic -Werror)
            target_compile_definitions(${_fault_object} PRIVATE ${_defines} WIREHAIR_BUILDING=1 WIREHAIR_TESTING=1)
        endforeach()
        foreach(_case IN ITEMS k3 k5 k8 contracts c_consumer profile_baseline borrowed borrowed_k8 borrowed_fault_baseline borrowed_fault_candidate metadata_baseline metadata_candidate)
            set(_test k8_ordinary_${_mode}_${_case})
            if(_case MATCHES "^k([358])$")
                set(_test_k "${CMAKE_MATCH_1}")
                add_executable(${_test} "${_small_test_generated}" $<TARGET_OBJECTS:${_object}>)
                target_compile_definitions(${_test} PRIVATE WIREHAIR_V2_SMALL_TEST_K=${_test_k})
            elseif(_case STREQUAL "contracts")
                add_executable(${_test} "${CMAKE_CURRENT_SOURCE_DIR}/ContractTest.cpp" $<TARGET_OBJECTS:${_object}>)
            elseif(_case STREQUAL "c_consumer")
                add_executable(${_test} "${_c_test_generated}" $<TARGET_OBJECTS:${_object}>)
                set_target_properties(${_test} PROPERTIES LINKER_LANGUAGE CXX)
            elseif(_case STREQUAL "profile_baseline")
                # The full legacy suite deliberately freezes ordinary K8=CERTIFIED.
                # Run it on baseline; candidate routing/contracts are tested separately.
                add_executable(${_test} "${_root}/codec/V2ProfileTest.cpp")
            elseif(_case STREQUAL "borrowed")
                add_executable(${_test} "${_root}/codec/V2BorrowedSourceTest.cpp" $<TARGET_OBJECTS:${_object}>)
            elseif(_case STREQUAL "borrowed_k8")
                add_executable(${_test} "${_borrowed_test_generated}" $<TARGET_OBJECTS:${_object}>)
            elseif(_case MATCHES "^borrowed_fault_(baseline|candidate)$")
                add_executable(${_test} "${_root}/test/V2BorrowedFacadeFaultTest.cpp"
                    $<TARGET_OBJECTS:k8_ordinary_fault_${_mode}_${CMAKE_MATCH_1}>)
            else()
                add_executable(${_test} "${_root}/bench/Wh2ProfileClassification/ValidationTest.cpp")
                if(_case STREQUAL "metadata_candidate")
                    target_sources(${_test} PRIVATE $<TARGET_OBJECTS:${_object}>)
                endif()
            endif()
            target_include_directories(${_test} PRIVATE "${_root}/include" "${_root}/test" "${_root}/codec" "${CMAKE_CURRENT_SOURCE_DIR}")
            target_compile_definitions(${_test} PRIVATE ${_defines})
            target_compile_options(${_test} PRIVATE ${_flags} -Wall -Wextra -Wpedantic -Werror)
            target_link_libraries(${_test} PRIVATE "${_base}/libwirehair.a" Threads::Threads)
            target_link_options(${_test} PRIVATE ${_link_flags} "-Wl,-Map,${_binary}/${_test}.map")
            add_test(NAME ${_test} COMMAND ${_test})
            set_tests_properties(${_test} PROPERTIES TIMEOUT 120)
            if(_mode STREQUAL "asan")
                set_tests_properties(${_test} PROPERTIES ENVIRONMENT
                    "ASAN_OPTIONS=detect_leaks=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1")
            endif()
        endforeach()
        # Engineering replay of the original entire K8 corpus through actual
        # ordinary selection. No new roots, rank screen, timing or sample count.
        set(_prototype "/tmp/wh2-k8-serialized.qQHyZqaa/${_mode}-k8")
        if(_mode STREQUAL "native")
            set(_prototype_sha "0ad76d0471b5b3fda81f40bec87e52e328f4cb6c0dca5f6dbdd90f4e84d583c6")
        elseif(_mode STREQUAL "scalar")
            set(_prototype_sha "4313edd59c522fd11b0a24b9b5042fe0df036cf45eaae204ad463b9b80dd1a10")
        else()
            set(_prototype_sha "6f1e29b12999a867fb909c3d017b95c97a3978955bfe5f75e2415bb791341c6b")
        endif()
        file(SHA256 "${_prototype}/libwh2_small_serialized.a" _actual)
        file(SHA256 "${_prototype}/Wh2K8NativeData.inc" _data_sha)
        if(NOT _actual STREQUAL _prototype_sha OR NOT _data_sha STREQUAL
            "99e2b20da1405755136311d58bd7391474866bdee1a10ccf22a2de7482b4aaea")
            message(FATAL_ERROR "Retained ${_mode} K8 prototype/fixture changed")
        endif()
        set(_parity_test k8_ordinary_${_mode}_corpus)
        add_executable(${_parity_test} "${_parity_generated}" "${_root}/bench/Wh2FrozenTrace.cpp"
            $<TARGET_OBJECTS:${_object}>)
        target_include_directories(${_parity_test} PRIVATE "${_root}/include" "${_root}/bench" "${_prototype}")
        target_compile_definitions(${_parity_test} PRIVATE ${_defines} WH2_SMALL_CODEC_K=8)
        target_compile_options(${_parity_test} PRIVATE ${_flags} -Wall -Wextra -Wpedantic -Werror)
        target_link_libraries(${_parity_test} PRIVATE "${_prototype}/libwh2_small_serialized.a"
            "${_base}/libwirehair.a" Threads::Threads)
        target_link_options(${_parity_test} PRIVATE ${_link_flags} "-Wl,-Map,${_binary}/${_parity_test}.map")
        add_test(NAME ${_parity_test} COMMAND ${_parity_test} --corpus)
        set_tests_properties(${_parity_test} PROPERTIES TIMEOUT 120)
        if(_mode STREQUAL "asan")
            set_tests_properties(${_parity_test} PROPERTIES ENVIRONMENT
                "ASAN_OPTIONS=detect_leaks=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1")
        endif()
    endforeach()
    foreach(_arm IN ITEMS baseline candidate)
        set(_test k8_ordinary_compat_${_arm})
        add_executable(${_test} "${_root}/test/V2CertifiedCompatibility.cpp"
            $<TARGET_OBJECTS:k8_ordinary_${_arm}>)
        target_include_directories(${_test} PRIVATE "${_root}/include")
        target_compile_definitions(${_test} PRIVATE WIREHAIR_STATIC=1)
        target_compile_options(${_test} PRIVATE -O3 -Wall -Wextra -Wpedantic -Werror)
        target_link_libraries(${_test} PRIVATE "${_qualified}/native/libwirehair.a" Threads::Threads)
        target_link_options(${_test} PRIVATE "-Wl,-Map,${_binary}/${_test}.map")
    endforeach()
    add_test(NAME k8_ordinary_certified_parity
        COMMAND "${CMAKE_COMMAND}"
            "-DBASELINE=$<TARGET_FILE:k8_ordinary_compat_baseline>"
            "-DCANDIDATE=$<TARGET_FILE:k8_ordinary_compat_candidate>"
            "-DOUTPUT_DIR=${_binary}"
            -P "${_root}/bench/Wh2ProfileClassification/CheckParity.cmake")
    set_tests_properties(k8_ordinary_certified_parity PROPERTIES TIMEOUT 120)
endif()
