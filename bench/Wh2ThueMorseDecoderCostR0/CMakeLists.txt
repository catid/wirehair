cmake_minimum_required(VERSION 3.20)
project(Wh2ThueMorseDecoderCostR0 LANGUAGES CXX)

option(WH2_TINY_GFNI_COST_SANITIZERS "Use ASan/UBSan; scientific worker is disabled" OFF)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR "Native single-configuration Linux/GNU build required")
endif()
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside)
if(_inside EQUAL 0)
    message(FATAL_ERROR "Use a fresh external build directory")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" _config)
set(_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_config}} ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_STATIC_LINKER_FLAGS}")
if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_config} OR
   _flags MATCHES "WH_COUNT|WIREHAIR_TESTING|WIREHAIR_V2_ENABLE_TEST_HOOKS|-DGF256|-flto|-fprofile|-m[a-z]")
    message(FATAL_ERROR "Counters, hooks, LTO, PGO and injected ISA flags are forbidden")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Threads REQUIRED)
enable_testing()

# Authenticate the qualified payload-only kernel and unchanged bridge/runtime.
set(_pins
    "CMakeLists.txt|52742b1c4253185688f85340cd2be13eb82df590c7e6cacef520768c0477468a"
    "gf256.cpp|fb9acc85b995714b0484fdfacf73427b3d4e33eef7c395f8745e876d01edfd90"
    "gf256.h|993b5de4a060b4efee4a45ac9ffc298e97b38c5efb00ead79dc21d87e8d470e2"
    "bench/Wh2ThueMorseNativeCodec.cpp|8e8e32e05394c16af0d599d83db4d319aa429d6c6f0377714cdfda540ff60c2d"
    "bench/Wh2ThueMorseNativeCodec.h|0fedb1d611d75a1235abf92b99761c5daee7ecb5b43efe8007c0f5a5113b037b"
    "bench/Wh2ThueMorseNativeCodecTest.cpp|f2b63428d6578b044eb0c50d5428c7a78447f4824ced854f50ee7925d20ada67"
    "bench/Wh2ThueMorseTinyPayloadR0.h|e38c818bfe3be9aa4a7f4bf83bd11c6bcc15f6d5e86a61f1ed9365fbaa0f4cf7"
    "bench/Wh2ThueMorseMatvecCostR0Bridge.cpp|10ac1344dc4632c89c0005be4da06c3b0d05d1c2c631b6b405f4f9af46180471"
    "bench/Wh2ThueMorseMatvecCostR0Bridge.h|1395c903184824c2b056f86e1eaa05792e6dad003633b73c9e50c3d4a04ee941"
    "bench/Wh2ThueMorseMatvecCostR0BridgeTest.cpp|d6faae5a7a580db2ec5ffd22b515fbc53c3b20b5ad893a050e87be357bb642cf"
    "bench/Wh2ThueMorseMatvecGfniR0.h|e515e652234e905b3a32ea403b995f5f173a68b487dcdc6c0eedb12e01a472b7"
    "bench/Wh2ThueMorseTinyGfniR0.h|3c562e03ea7e90061a80c598c22cf8060fa2f32151c7a2fd7120d8764df39068"
    "bench/Wh2ThueMorseTinyGfniR0.cpp|d596555cefc7486f4a101f9bcd10b4ac43f634bdbdb91e3d0f85619f229b54ea"
    "bench/Wh2ThueMorseTinyGfniR0Test.cpp|0b02aa2e831b8a12e56a64151a1715a275a0474d9cfd82cd902d8f6768d65280"
    "bench/Wh2ThueMorseTinyGfniR0Codec.cpp|0609762920b7ddfd108a77473da2538714aecd8e65427bc624a4d3659cf8d71e"
    "bench/Wh2ThueMorseTinyGfniR0CodecTest.cpp|ebcd7924099d8efc020d218ef47112ec752d47624a7f44cf201b9e418dfac45c"
    "bench/Wh2ThueMorseTinyGfniR0/CMakeLists.txt|ba8efd5c1c5f0b1c98bd8f246fc55c893f816aa5be5f4bcba5cc663a721cbe6c")
foreach(_pin IN LISTS _pins)
    string(REPLACE "|" ";" _fields "${_pin}")
    list(GET _fields 0 _relative)
    list(GET _fields 1 _expected)
    execute_process(COMMAND /usr/bin/test -f "${_root}/${_relative}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_root}/${_relative}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Frozen input must be regular: ${_relative}")
    endif()
    file(SHA256 "${_root}/${_relative}" _actual)
    if(NOT _actual STREQUAL _expected)
        message(FATAL_ERROR "Frozen input changed: ${_relative}")
    endif()
endforeach()

# Record the whole native repository dependency closure, plus build/controller
# sources, before generating or compiling anything. No old cost worker is used.
set(_native_inputs
    CMakeLists.txt cmake/wirehair.pc.in cmake/wirehairConfig.cmake.in
    wirehair.cpp include/wirehair/wirehair.h include/wirehair/wirehair.hpp
    gf256.cpp gf256.h WirehairCodec.cpp WirehairCodec.h WirehairCompiler.h
    WirehairEnvironment.h WirehairHeavy.h WirehairTools.cpp WirehairTools.h
    WirehairDenseFixups.inc WirehairPeelFixups.inc
    codec/WirehairV2Codec.cpp codec/WirehairV2Codec.h
    codec/WirehairV2Peel.cpp codec/WirehairV2Peel.h
    codec/WirehairV2Plan.cpp codec/WirehairV2Plan.h
    codec/WirehairV2Policy.cpp codec/WirehairV2Policy.h
    codec/WirehairV2Precode.cpp codec/WirehairV2Precode.h
    codec/WirehairV2PrecodeDecode.cpp codec/WirehairV2PrecodeDecode.h
    codec/WirehairV2PrecodeEncode.cpp codec/WirehairV2PrecodeEncode.h
    codec/WirehairV2Profile.cpp codec/WirehairV2Seeds.cpp codec/WirehairV2Seeds.h
    codec/WirehairV2Solve.cpp codec/WirehairV2Solve.h
    bench/Wh2FrozenTrace.cpp bench/Wh2FrozenTrace.h
    bench/Wh2PublicBorrowedTargetIdentity.cpp bench/Wh2PublicBorrowedTargetIdentity.h
    bench/Wh2RdpruTargetIdentityV2.cpp bench/Wh2RdpruTargetIdentityV2.h
    bench/Wh2ThueMorseNativeCodec.cpp bench/Wh2ThueMorseNativeCodec.h
    bench/Wh2ThueMorseNativeCodecTest.cpp
    bench/Wh2ThueMorseTinyPayloadR0.h bench/Wh2ThueMorseMatvecGfniR0.h
    bench/Wh2ThueMorseTinyGfniR0.cpp bench/Wh2ThueMorseTinyGfniR0.h
    bench/Wh2ThueMorseTinyGfniR0Test.cpp bench/Wh2ThueMorseTinyGfniR0Codec.cpp
    bench/Wh2ThueMorseTinyGfniR0CodecTest.cpp bench/Wh2ThueMorseTinyGfniR0/CMakeLists.txt
    bench/Wh2ThueMorseMatvecCostR0BridgeTest.cpp
    bench/Wh2ThueMorseDecoderCostR0.cpp bench/Wh2ThueMorseMatvecCostR0Bridge.cpp
    bench/Wh2ThueMorseMatvecCostR0Bridge.h bench/Wh2ThueMorseTinyGfniCostR0BridgeTest.cpp
    bench/Wh2ThueMorseDecoderCostR0/CMakeLists.txt
    bench/Wh2ThueMorseDecoderCostR0.py bench/test_Wh2ThueMorseDecoderCostR0.py)
list(SORT _native_inputs)
set(_manifest "{\"schema\":\"wirehair.wh2.thue-decoder-cost-r0.native-inputs.v1\",\"sources\":{")
set(_separator "")
set(_input_hashes "")
foreach(_relative IN LISTS _native_inputs)
    set(_path "${_root}/${_relative}")
    execute_process(COMMAND /usr/bin/test -f "${_path}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_path}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Native inputs must be regular source files: ${_relative}")
    endif()
    file(SHA256 "${_path}" _digest)
    list(APPEND _input_hashes "${_digest}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_path}")
    string(APPEND _manifest "${_separator}\"${_relative}\":\"${_digest}\"")
    set(_separator ",")
endforeach()
string(APPEND _manifest "}}\n")
set(_manifest_path "${_binary}/native-inputs.json")
if(IS_SYMLINK "${_manifest_path}" OR IS_DIRECTORY "${_manifest_path}")
    message(FATAL_ERROR "Invalid native input manifest")
endif()
if(EXISTS "${_manifest_path}")
    execute_process(COMMAND /usr/bin/test -f "${_manifest_path}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Native input manifest must be regular")
    endif()
    file(SIZE "${_manifest_path}" _manifest_size)
    if(_manifest_size GREATER 65536)
        message(FATAL_ERROR "Oversized native input manifest")
    endif()
    file(READ "${_manifest_path}" _old_manifest)
    if(NOT _old_manifest STREQUAL _manifest)
        message(FATAL_ERROR "Native inputs changed: use a fresh build directory")
    endif()
else()
    if(EXISTS "${_binary}/build.ninja" OR EXISTS "${_binary}/Makefile")
        message(FATAL_ERROR "Cannot adopt an unrecorded existing build")
    endif()
    file(WRITE "${_manifest_path}" "${_manifest}")
endif()

# Reuse the exact .74 payload-only overlay and its neutral qualification targets.
set(WH2_TINY_GFNI_SANITIZERS ${WH2_TINY_GFNI_COST_SANITIZERS} CACHE BOOL "" FORCE)
add_subdirectory("${_root}/bench/Wh2ThueMorseTinyGfniR0" "${_binary}/codec-neutral")

# BUILD_CODEC_V2 disables internal experiments, not the public WH2 codec. The
# unchanged archive supplies both real public APIs. None of its old workers run.
foreach(_off IN ITEMS BUILD_SHARED_LIBS BUILD_TESTS BUILD_CODEC_V2 MARCH_NATIVE
        WIREHAIR_BUILD_BOTH WIREHAIR_BUILD_TOOLS WIREHAIR_BUILD_BENCHMARKS
        WIREHAIR_ENABLE_SCHEDULED_TESTS WIREHAIR_ENABLE_PUBLIC_BORROWED_R1_NATIVE_CONFIG
        WIREHAIR_ENABLE_THUE_MORSE_NATIVE_R0 WIREHAIR_ENABLE_LIBFUZZER
        WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER)
    set(${_off} OFF CACHE BOOL "Disabled in isolated tiny GFNI cost build" FORCE)
endforeach()
set(WIREHAIR_STATIC_PIC ON CACHE BOOL "" FORCE)
set(WIREHAIR_STRICT_WARNINGS ON CACHE BOOL "" FORCE)
set(WH_LTO OFF CACHE STRING "" FORCE)
set(WH_PGO_MODE OFF CACHE STRING "" FORCE)
add_subdirectory("${_root}" "${_binary}/production" EXCLUDE_FROM_ALL)
# The unchanged project defines these five offline generators even with tools
# disabled. They are not built here and must not enter the actual-TU receipt.
# This target-scoped property requires CMake 3.20; no compile/link flag changes.
foreach(_generator IN ITEMS gen_small_dseeds gen_peel_seeds gen_most_dseeds gen_dcounts gen_tables)
    if(NOT TARGET ${_generator})
        message(FATAL_ERROR "Expected unbuilt offline generator is absent")
    endif()
    get_target_property(_excluded ${_generator} EXCLUDE_FROM_ALL)
    if(NOT _excluded)
        message(FATAL_ERROR "Offline generator unexpectedly participates in the build")
    endif()
    set_target_properties(${_generator} PROPERTIES EXPORT_COMPILE_COMMANDS OFF)
endforeach()

# Reconstruction is confined to the already selected lookup. It is not a codec
# run, a parameter selection, or a correctness/performance campaign.
execute_process(COMMAND /usr/bin/python3 "${_root}/bench/Wh2ThueMorseDecoderCostR0.py"
    lookup "${_binary}/lookup.bin" RESULT_VARIABLE _lookup_result TIMEOUT 30)
if(NOT _lookup_result STREQUAL "0")
    message(FATAL_ERROR "Fixed lookup reconstruction failed")
endif()
execute_process(COMMAND /usr/bin/test -f "${_binary}/lookup.bin"
    RESULT_VARIABLE _lookup_regular TIMEOUT 1)
if(IS_SYMLINK "${_binary}/lookup.bin" OR NOT _lookup_regular STREQUAL "0")
    message(FATAL_ERROR "Lookup must be a regular non-symlink file")
endif()
file(SIZE "${_binary}/lookup.bin" _lookup_size)
file(SHA256 "${_binary}/lookup.bin" _lookup_sha)
if(NOT _lookup_size EQUAL 39936 OR
   NOT _lookup_sha STREQUAL "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d")
    message(FATAL_ERROR "Selected lookup bytes differ")
endif()

set(_bridge_objects "")
set(_bridges "")
foreach(_arm RANGE 0 3)
    set(_target "tiny_gfni_cost_bridge_${_arm}")
    add_library(${_target} OBJECT "${_root}/bench/Wh2ThueMorseMatvecCostR0Bridge.cpp")
    target_compile_definitions(${_target} PRIVATE WH2_MATVEC_COST_ARM=${_arm})
    if(_arm EQUAL 1)
        target_compile_definitions(${_target} PRIVATE
            wh2_thue_matvec_codec_r0=wh2_thue_tiny_gfni_r0)
    endif()
    target_link_libraries(${_target} PRIVATE wirehair)
    list(APPEND _bridges ${_target})
    list(APPEND _bridge_objects $<TARGET_OBJECTS:${_target}>)
endforeach()
set(_codec_objects $<TARGET_OBJECTS:tiny_gfni_baseline> $<TARGET_OBJECTS:tiny_gfni_candidate>
    $<TARGET_OBJECTS:tiny_gfni_runtime>)
add_executable(wh2_thue_decoder_cost_r0
    "${_root}/bench/Wh2ThueMorseDecoderCostR0.cpp"
    "${_root}/bench/Wh2RdpruTargetIdentityV2.cpp"
    "${_root}/bench/Wh2PublicBorrowedTargetIdentity.cpp"
    "${_root}/bench/Wh2FrozenTrace.cpp" ${_codec_objects} ${_bridge_objects})
target_compile_definitions(wh2_thue_decoder_cost_r0 PRIVATE
    WH2_TINY_GFNI_COST_LOOKUP_PATH="${_binary}/lookup.bin"
    WH2_TINY_GFNI_COST_SANITIZERS=$<BOOL:${WH2_TINY_GFNI_COST_SANITIZERS}>)
target_link_options(wh2_thue_decoder_cost_r0 PRIVATE "-Wl,-Map,${_binary}/tiny-gfni-cost.map")
add_executable(wh2_tiny_gfni_cost_bridge_test
    "${_root}/bench/Wh2ThueMorseTinyGfniCostR0BridgeTest.cpp" ${_codec_objects} ${_bridge_objects})
target_link_options(wh2_tiny_gfni_cost_bridge_test PRIVATE "-Wl,-Map,${_binary}/tiny-gfni-cost-bridge-test.map")
set(_executables wh2_thue_decoder_cost_r0 wh2_tiny_gfni_cost_bridge_test)
foreach(_target IN ITEMS wirehair ${_bridges} ${_executables})
    target_include_directories(${_target} PRIVATE "${_root}" "${_root}/bench" "${_root}/include")
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-lto)
    if(WH2_TINY_GFNI_COST_SANITIZERS)
        target_compile_options(${_target} PRIVATE -O1 -g -fno-omit-frame-pointer
            -fsanitize=address,undefined -fno-pie)
    else()
        target_compile_options(${_target} PRIVATE -O3)
    endif()
endforeach()
foreach(_target IN LISTS _executables)
    # The .74 object supplies all GF symbols. The production archive's
    # gf256.cpp.o must remain unextracted (audit the actual link map).
    target_link_libraries(${_target} PRIVATE wirehair Threads::Threads)
    target_link_options(${_target} PRIVATE -fno-lto)
    if(WH2_TINY_GFNI_COST_SANITIZERS)
        target_link_options(${_target} PRIVATE -fsanitize=address,undefined -no-pie)
    endif()
endforeach()
add_test(NAME tiny_gfni_cost_neutral COMMAND wh2_thue_decoder_cost_r0 --neutral)
add_test(NAME tiny_gfni_cost_bridge_neutral COMMAND wh2_tiny_gfni_cost_bridge_test --require-gfni)
set_tests_properties(tiny_gfni_cost_neutral tiny_gfni_cost_bridge_neutral PROPERTIES TIMEOUT 60)
if(WH2_TINY_GFNI_COST_SANITIZERS)
    set_tests_properties(tiny_gfni_cost_neutral tiny_gfni_cost_bridge_neutral PROPERTIES ENVIRONMENT
        "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
else()
    find_program(_cost_qemu qemu-x86_64)
    if(_cost_qemu)
        add_test(NAME tiny_gfni_cost_bridge_haswell COMMAND "${_cost_qemu}" -cpu Haswell
            $<TARGET_FILE:wh2_tiny_gfni_cost_bridge_test> --expect-scalar)
        set_tests_properties(tiny_gfni_cost_bridge_haswell PROPERTIES TIMEOUT 60)
    endif()
endif()

# A configure-time edit must not produce a manifest for one set of bytes and
# compilation rules for another. The controller additionally checks depfiles.
set(_index 0)
foreach(_relative IN LISTS _native_inputs)
    list(GET _input_hashes ${_index} _before)
    execute_process(COMMAND /usr/bin/test -f "${_root}/${_relative}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_root}/${_relative}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Native input changed type during configuration")
    endif()
    file(SHA256 "${_root}/${_relative}" _after)
    if(NOT _before STREQUAL _after)
        message(FATAL_ERROR "Native input changed during configuration: ${_relative}")
    endif()
    math(EXPR _index "${_index} + 1")
endforeach()
