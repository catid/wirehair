cmake_minimum_required(VERSION 3.20)
project(Wh2OrderedRhsReleaseR0 LANGUAGES C CXX)

option(WH2_RELEASE_SANITIZERS "Neutral ASan/UBSan build; worker disabled" OFF)
option(WH2_RELEASE_SCALAR "Neutral portable GF backend; worker disabled" OFF)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "This bounded diagnostic requires native Linux/GNU")
endif()
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside)
if(_inside EQUAL 0)
    message(FATAL_ERROR "Use an external build directory")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Threads REQUIRED)
add_compile_options(-Wall -Wextra -Wpedantic -Werror -fno-lto -g1 -fPIC)
add_link_options(-fno-lto)
if(WH2_RELEASE_SCALAR)
    # gf256.h uses ANDROID to select its portable, non-NEON uint64 backend
    # on this x86 host. This is a neutral backend check, not an Android claim.
    add_compile_definitions(ANDROID=1)
endif()
if(WH2_RELEASE_SANITIZERS)
    add_compile_options(-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined -no-pie)
else()
    add_compile_options(-O3)
    add_link_options(-no-pie)
endif()
enable_testing()

# The candidate is a build artifact, not an edit to the production solver.
set(_solve "${_root}/codec/WirehairV2Solve.cpp")
file(SHA256 "${_solve}" _solve_sha)
if(NOT _solve_sha STREQUAL "14e3d5c07c93f3b7ea9a1d6d63d197fe38cddce1d0915a0cdd5c19e5c3b1b058")
    message(FATAL_ERROR "Frozen production solver changed")
endif()
file(READ "${_solve}" _source)
set(_anchor [=[        intermediate_blocks_out.swap(values);
        if (stats) {
            *stats = st;
        }
        return Wirehair_Success;]=])
string(FIND "${_source}" "${_anchor}" _anchor_at)
if(_anchor_at LESS 0)
    message(FATAL_ERROR "Missing successful cold-solve anchor")
endif()
string(REPLACE "${_anchor}" "" _without_anchor "${_source}")
string(LENGTH "${_source}" _source_len)
string(LENGTH "${_without_anchor}" _without_len)
string(LENGTH "${_anchor}" _anchor_len)
math(EXPR _removed "${_source_len} - ${_without_len}")
if(NOT _removed EQUAL _anchor_len)
    message(FATAL_ERROR "Ambiguous successful cold-solve anchor")
endif()
set(_release [=[        // All RHS users are finished. Release in allocation order so
        // contiguous free chunks can coalesce before exposing the heap top.
        std::vector<uint8_t>().swap(pivot_rhs);
        std::vector<uint8_t>().swap(rhs);
        std::vector<uint8_t>().swap(heavy_rhs);
        std::vector<uint8_t>().swap(quotient_pivot_rhs);
        std::vector<uint8_t>().swap(quotient_rhs);

]=])
string(REPLACE "${_anchor}" "${_release}${_anchor}" _candidate "${_source}")
set(_generated "${_binary}/WirehairV2SolveOrdered.cpp")
file(WRITE "${_generated}" "${_candidate}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_solve}")

foreach(_off IN ITEMS BUILD_SHARED_LIBS BUILD_TESTS BUILD_CODEC_V2 MARCH_NATIVE
        WIREHAIR_BUILD_BOTH WIREHAIR_BUILD_TOOLS WIREHAIR_BUILD_BENCHMARKS
        WIREHAIR_ENABLE_SCHEDULED_TESTS WIREHAIR_ENABLE_PUBLIC_BORROWED_R1_NATIVE_CONFIG
        WIREHAIR_ENABLE_THUE_MORSE_NATIVE_R0 WIREHAIR_ENABLE_LIBFUZZER
        WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER)
    set(${_off} OFF CACHE BOOL "Isolated release-order build" FORCE)
endforeach()
set(WIREHAIR_STATIC_PIC ON CACHE BOOL "" FORCE)
set(WIREHAIR_STRICT_WARNINGS ON CACHE BOOL "" FORCE)
set(WH_LTO OFF CACHE STRING "" FORCE)
set(WH_PGO_MODE OFF CACHE STRING "" FORCE)
add_subdirectory("${_root}" "${_binary}/production" EXCLUDE_FROM_ALL)

# Explicit objects precede the shared original archive in both links. Its
# original Solve.cpp object must not be extracted in either executable.
add_library(release_baseline OBJECT "${_solve}")
add_library(release_candidate OBJECT "${_generated}")
add_library(release_candidate_hooks OBJECT "${_generated}")
target_compile_definitions(release_candidate_hooks PRIVATE WIREHAIR_V2_ENABLE_TEST_HOOKS=1)
foreach(_target IN ITEMS release_baseline release_candidate release_candidate_hooks)
    target_include_directories(${_target} PRIVATE "${_root}/codec" "${_root}" "${_root}/include")
    target_compile_definitions(${_target} PRIVATE WIREHAIR_BUILDING=1 WIREHAIR_STATIC=1)
endforeach()
add_library(release_worker OBJECT "${_root}/bench/Wh2OrderedRhsReleaseR0.cpp")
target_include_directories(release_worker PRIVATE "${_root}/include")
target_compile_definitions(release_worker PRIVATE
    WH2_RELEASE_NEUTRAL=$<OR:$<BOOL:${WH2_RELEASE_SANITIZERS}>,$<BOOL:${WH2_RELEASE_SCALAR}>>)
foreach(_arm IN ITEMS baseline candidate)
    add_executable(release_${_arm}_worker $<TARGET_OBJECTS:release_worker> $<TARGET_OBJECTS:release_${_arm}>)
    target_link_libraries(release_${_arm}_worker PRIVATE wirehair)
    target_link_options(release_${_arm}_worker PRIVATE "-Wl,-Map,${_binary}/${_arm}.map")
    add_test(NAME release_${_arm}_selftest COMMAND release_${_arm}_worker --selftest)
endforeach()

# Existing test-hook policy, with only the candidate solve object substituted.
set(_policy_sources "")
foreach(_name IN ITEMS Codec Peel Plan Policy Precode PrecodeDecode PrecodeEncode Seeds)
    list(APPEND _policy_sources "${_root}/codec/WirehairV2${_name}.cpp")
endforeach()
add_library(release_policy STATIC ${_policy_sources})
target_include_directories(release_policy PUBLIC "${_root}/codec" "${_root}" "${_root}/include")
target_compile_definitions(release_policy PUBLIC WIREHAIR_V2_ENABLE_TEST_HOOKS=1)
target_link_libraries(release_policy PUBLIC wirehair)
foreach(_test IN ITEMS PrecodeSolve PrecodeDecode PrecodeRoundTrip PrecodeEncode V2Transactional V2NoInit Policy)
    add_executable(release_${_test}_test "${_root}/codec/${_test}Test.cpp" $<TARGET_OBJECTS:release_candidate_hooks>)
    if(_test STREQUAL "PrecodeSolve")
        target_sources(release_${_test}_test PRIVATE "${_root}/codec/V2TinyDenseOracle.cpp")
    endif()
    target_link_libraries(release_${_test}_test PRIVATE release_policy Threads::Threads)
    add_test(NAME release_${_test}_test COMMAND release_${_test}_test)
endforeach()
foreach(_test IN ITEMS V2Profile V2BorrowedSource)
    add_executable(release_${_test}_test "${_root}/codec/${_test}Test.cpp" $<TARGET_OBJECTS:release_candidate>)
    target_link_libraries(release_${_test}_test PRIVATE wirehair Threads::Threads)
    add_test(NAME release_${_test}_test COMMAND release_${_test}_test)
endforeach()
get_property(_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${_tests} PROPERTIES TIMEOUT 60)
if(WH2_RELEASE_SANITIZERS)
    set_tests_properties(${_tests} PROPERTIES ENVIRONMENT
        "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
endif()
