cmake_minimum_required(VERSION 3.16)
project(Wh2ThueMorseMatvecCodecR0 LANGUAGES CXX)

# Standalone neutral codec tests only. No selected lookup or timing target.
option(WH2_MATVEC_CODEC_SANITIZERS "Enable ASan/UBSan for neutral codec tests" OFF)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR "Native single-configuration Linux/GNU build required")
endif()
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside)
if(_inside EQUAL 0)
    message(FATAL_ERROR "Use an external build directory")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" _config)
set(_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_config}} ${CMAKE_EXE_LINKER_FLAGS}")
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

# Authenticate all eight frozen inputs before reading or transforming the codec.
# The new wrapper is independently reviewed source, not a frozen old input.
set(_pins
    "WirehairEnvironment.h|5e580558a11944fde700c85381fa34fb2bed544a71f9d30ab6b05faf58d5c161"
    "bench/Wh2ThueMorseNativeCodec.cpp|8e8e32e05394c16af0d599d83db4d319aa429d6c6f0377714cdfda540ff60c2d"
    "bench/Wh2ThueMorseNativeCodec.h|0fedb1d611d75a1235abf92b99761c5daee7ecb5b43efe8007c0f5a5113b037b"
    "bench/Wh2ThueMorseNativeCodecTest.cpp|f2b63428d6578b044eb0c50d5428c7a78447f4824ced854f50ee7925d20ada67"
    "bench/Wh2ThueMorseMatvecGfniR0.cpp|977fad90ab6b90d56aab3ee2e1e0b48dd44ab65d4cd8517104f048cebf0e389c"
    "bench/Wh2ThueMorseMatvecGfniR0.h|e515e652234e905b3a32ea403b995f5f173a68b487dcdc6c0eedb12e01a472b7"
    "gf256.cpp|fb9acc85b995714b0484fdfacf73427b3d4e33eef7c395f8745e876d01edfd90"
    "gf256.h|993b5de4a060b4efee4a45ac9ffc298e97b38c5efb00ead79dc21d87e8d470e2")
list(SORT _pins)
set(_sources "{")
set(_separator "")
foreach(_pin IN LISTS _pins)
    string(REPLACE "|" ";" _fields "${_pin}")
    list(GET _fields 0 _relative)
    list(GET _fields 1 _expected)
    set(_path "${_root}/${_relative}")
    execute_process(COMMAND /usr/bin/test -f "${_path}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_path}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Frozen inputs must be regular non-symlink files: ${_relative}")
    endif()
    file(SHA256 "${_path}" _actual)
    if(NOT _actual STREQUAL _expected)
        message(FATAL_ERROR "Frozen source hash differs: ${_relative}")
    endif()
    if(_relative STREQUAL "bench/Wh2ThueMorseNativeCodec.cpp")
        set(_codec_expected "${_expected}")
    endif()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_path}")
    string(APPEND _sources "${_separator}\"${_relative}\":\"${_actual}\"")
    set(_separator ",")
endforeach()
string(APPEND _sources "}")

# Match the complete original function literally, not a broad regular expression.
# Only its body changes. Macro-renaming Apply would also rename its definition.
set(_old_apply [=[void Apply(const std::uint8_t* matrix, std::uint8_t vector[kSourceCount])
{
    std::uint8_t result[kSourceCount];
    for (unsigned r = 0; r < kSourceCount; ++r) {
        std::uint8_t value = gf256_mul(vector[0], matrix[r * kSourceCount]);
        for (unsigned c = 1; c < kSourceCount; ++c)
            value ^= gf256_mul(vector[c], matrix[r * kSourceCount + c]);
        result[r] = value;
    }
    std::memcpy(vector, result, kSourceCount);
}]=])
set(_new_apply [=[void Apply(const std::uint8_t* matrix, std::uint8_t vector[kSourceCount])
{
    ::wh2_matvec_gfni_r0::Apply(matrix, vector);
}]=])
file(READ "${_root}/bench/Wh2ThueMorseNativeCodec.cpp" _original)
string(SHA256 _original_sha "${_original}")
if(NOT _original_sha STREQUAL _codec_expected)
    message(FATAL_ERROR "Original codec changed while being read")
endif()
string(FIND "${_original}" "${_old_apply}" _first)
if(_first LESS 0)
    message(FATAL_ERROR "Exact original Apply function is absent")
endif()
string(LENGTH "${_old_apply}" _old_size)
math(EXPR _after "${_first} + ${_old_size}")
string(SUBSTRING "${_original}" ${_after} -1 _tail)
string(FIND "${_tail}" "${_old_apply}" _second)
if(NOT _second EQUAL -1)
    message(FATAL_ERROR "Exact original Apply function is not unique")
endif()
string(REPLACE "${_old_apply}" "${_new_apply}" _transformed "${_original}")
set(_generated "#include \"Wh2ThueMorseMatvecGfniR0.h\"\n#define wh2_thue_native wh2_thue_matvec_codec_r0\n${_transformed}\n#undef wh2_thue_native\n")
string(SHA256 _generated_sha "${_generated}")
string(LENGTH "${_generated}" _generated_bytes)
set(_generated_name "Wh2ThueMorseMatvecCodecR0.generated.cpp")
set(_generated_path "${_binary}/${_generated_name}")
set(_manifest "{\"generated\":{\"bytes\":${_generated_bytes},\"path\":\"${_generated_name}\",\"sha256\":\"${_generated_sha}\"},\"schema\":\"wirehair.wh2.thue-matvec-codec-r0.neutral-inputs.v1\",\"sources\":${_sources}}\n")

# Refuse to overwrite or follow unexpected existing generated artifacts.
function(_write_identical_or_new _path _contents)
    if(IS_SYMLINK "${_path}" OR IS_DIRECTORY "${_path}")
        message(FATAL_ERROR "Invalid generated artifact: ${_path}")
    endif()
    if(EXISTS "${_path}")
        execute_process(COMMAND /usr/bin/test -f "${_path}"
            RESULT_VARIABLE _regular TIMEOUT 1)
        if(NOT _regular STREQUAL "0")
            message(FATAL_ERROR "Generated artifact must be regular: ${_path}")
        endif()
        file(SIZE "${_path}" _bytes)
        if(_bytes GREATER 65536)
            message(FATAL_ERROR "Oversized generated artifact: ${_path}")
        endif()
        file(READ "${_path}" _existing)
        if(NOT _existing STREQUAL _contents)
            message(FATAL_ERROR "Generated artifact differs; use a fresh build directory")
        endif()
    else()
        file(WRITE "${_path}" "${_contents}")
    endif()
endfunction()
_write_identical_or_new("${_generated_path}" "${_generated}")
_write_identical_or_new("${_binary}/matvec-codec-source-manifest.json" "${_manifest}")
file(SHA256 "${_generated_path}" _written_sha)
if(NOT _written_sha STREQUAL _generated_sha)
    message(FATAL_ERROR "Generated source hash mismatch after publication")
endif()

# The .72 TU includes the unchanged gf256.cpp exactly once. Never add/link a
# separate gf256.cpp object; both codec arms use this same runtime definition.
add_library(matvec_codec_gf STATIC "${_root}/bench/Wh2ThueMorseMatvecGfniR0.cpp")
target_link_libraries(matvec_codec_gf PUBLIC Threads::Threads)
add_library(matvec_codec_baseline OBJECT "${_root}/bench/Wh2ThueMorseNativeCodec.cpp")
add_library(matvec_codec_candidate OBJECT "${_generated_path}")
add_executable(wh2_matvec_codec_baseline_test
    "${_root}/bench/Wh2ThueMorseMatvecCodecR0Test.cpp"
    $<TARGET_OBJECTS:matvec_codec_baseline>)
add_executable(wh2_matvec_codec_candidate_test
    "${_root}/bench/Wh2ThueMorseMatvecCodecR0Test.cpp"
    $<TARGET_OBJECTS:matvec_codec_candidate>)
target_compile_definitions(wh2_matvec_codec_baseline_test PRIVATE WH2_MATVEC_CODEC_CANDIDATE=0)
target_compile_definitions(wh2_matvec_codec_candidate_test PRIVATE WH2_MATVEC_CODEC_CANDIDATE=1)
set(_executables wh2_matvec_codec_baseline_test wh2_matvec_codec_candidate_test)
foreach(_target IN ITEMS matvec_codec_gf matvec_codec_baseline matvec_codec_candidate ${_executables})
    target_include_directories(${_target} PRIVATE "${_root}" "${_root}/bench")
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-lto)
    if(WH2_MATVEC_CODEC_SANITIZERS)
        target_compile_options(${_target} PRIVATE -O1 -g -fno-omit-frame-pointer
            -fsanitize=address,undefined -fno-pie)
    else()
        target_compile_options(${_target} PRIVATE -O3)
    endif()
endforeach()
foreach(_target IN LISTS _executables)
    target_link_libraries(${_target} PRIVATE matvec_codec_gf)
    target_link_options(${_target} PRIVATE -fno-lto "-Wl,-Map,${_binary}/${_target}.map")
    if(WH2_MATVEC_CODEC_SANITIZERS)
        target_link_options(${_target} PRIVATE -fsanitize=address,undefined -no-pie)
    endif()
    add_test(NAME ${_target} COMMAND ${_target} --require-gfni)
    set_tests_properties(${_target} PROPERTIES TIMEOUT 60)
    if(WH2_MATVEC_CODEC_SANITIZERS)
        set_tests_properties(${_target} PROPERTIES ENVIRONMENT
            "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
    endif()
endforeach()
if(NOT WH2_MATVEC_CODEC_SANITIZERS)
    find_program(_qemu qemu-x86_64)
    if(_qemu)
        foreach(_target IN LISTS _executables)
            add_test(NAME ${_target}_haswell COMMAND "${_qemu}" -cpu Haswell
                $<TARGET_FILE:${_target}> --expect-scalar)
            set_tests_properties(${_target}_haswell PROPERTIES TIMEOUT 60)
        endforeach()
    endif()
endif()
