foreach(required IN ITEMS PROJECT_SOURCE_DIR TEST_ROOT TEST_GENERATOR
        TEST_C_COMPILER TEST_CXX_COMPILER TEST_CONFIG TEST_MULTI_CONFIG
        TEST_LINKER_SENTINEL TEST_LINKER_PATTERN TEST_NATIVE_SUPPORTED
        TEST_PIC_FLAG_REQUIRED)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

function(run_checked output_var label)
    execute_process(
        COMMAND ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 120)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR
            "${label} failed (${result})\nstdout:\n${out}\nstderr:\n${err}")
    endif()
    set(${output_var} "${out}\n${err}" PARENT_SCOPE)
endfunction()

function(require_match value pattern label)
    if(NOT "${value}" MATCHES "${pattern}")
        message(FATAL_ERROR "${label}: missing '${pattern}'\n${value}")
    endif()
endfunction()

function(require_cli_rejected_before_output executable label)
    execute_process(
        COMMAND "${executable}" --no-benchmarks --heavy-trials ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 10)
    if(result EQUAL 0)
        message(FATAL_ERROR "${label}: invalid value was accepted")
    endif()
    if(NOT out STREQUAL "")
        message(FATAL_ERROR
            "${label}: generator produced output before rejecting input:\n${out}")
    endif()
    require_match("${err}" "Usage:" "${label} diagnostic")
endfunction()

set(generator_args -G "${TEST_GENERATOR}")
if(DEFINED TEST_GENERATOR_PLATFORM AND NOT TEST_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND generator_args -A "${TEST_GENERATOR_PLATFORM}")
endif()
if(DEFINED TEST_GENERATOR_TOOLSET AND NOT TEST_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND generator_args -T "${TEST_GENERATOR_TOOLSET}")
endif()
set(compiler_args
    "-DCMAKE_C_COMPILER=${TEST_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}")

file(REMOVE_RECURSE "${TEST_ROOT}")
file(MAKE_DIRECTORY "${TEST_ROOT}")

# Standard and custom single-configuration flag variables must survive intact.
foreach(config IN ITEMS Debug Release RelWithDebInfo AuditConfig)
    string(TOUPPER "${config}" config_upper)
    set(build_dir "${TEST_ROOT}/flags-${config}")
    set(sentinel "WIREHAIR_${config_upper}_SENTINEL")
    if(TEST_MULTI_CONFIG)
        set(config_args "-DCMAKE_CONFIGURATION_TYPES=${config}")
    else()
        set(config_args "-DCMAKE_BUILD_TYPE=${config}")
    endif()
    run_checked(configure_output "${config} configure"
        "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${build_dir}"
        ${generator_args} ${compiler_args} ${config_args}
        "-DCMAKE_CXX_FLAGS_${config_upper}=-D${sentinel}"
        "-DCMAKE_SHARED_LINKER_FLAGS_${config_upper}=${TEST_LINKER_SENTINEL}"
        -DBUILD_SHARED_LIBS=ON
        -DBUILD_TESTS=OFF
        -DBUILD_CODEC_V2=OFF
        -DMARCH_NATIVE=OFF
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
    run_checked(build_output "${config} build"
        "${CMAKE_COMMAND}" --build "${build_dir}"
        --config "${config}" --target wirehair --verbose --parallel 2)
    file(READ "${build_dir}/compile_commands.json" compile_commands)
    require_match("${compile_commands}" "${sentinel}"
        "${config} compile flags")
    require_match("${build_output}" "${TEST_LINKER_PATTERN}"
        "${config} shared linker flags")
    if(compile_commands MATCHES "-march=native")
        message(FATAL_ERROR "portable ${config} build used -march=native")
    endif()
endforeach()

# Exercise a true multi-config custom configuration when Ninja is available.
if(TEST_GENERATOR MATCHES "Ninja")
    set(multi_dir "${TEST_ROOT}/flags-multi")
    run_checked(multi_configure_output "multi-config configure"
        "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${multi_dir}"
        -G "Ninja Multi-Config" ${compiler_args}
        "-DCMAKE_CONFIGURATION_TYPES=Debug\;AuditConfig"
        "-DCMAKE_CXX_FLAGS_AUDITCONFIG=-DWIREHAIR_MULTI_SENTINEL"
        "-DCMAKE_SHARED_LINKER_FLAGS_AUDITCONFIG=${TEST_LINKER_SENTINEL}"
        -DBUILD_SHARED_LIBS=ON
        -DBUILD_TESTS=OFF
        -DBUILD_CODEC_V2=OFF
        -DMARCH_NATIVE=OFF
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
    run_checked(multi_build_output "multi-config custom build"
        "${CMAKE_COMMAND}" --build "${multi_dir}"
        --config AuditConfig --target wirehair --verbose --parallel 2)
    file(READ "${multi_dir}/compile_commands.json" multi_commands)
    require_match("${multi_commands}" "WIREHAIR_MULTI_SENTINEL"
        "multi-config compile flags")
    require_match("${multi_build_output}" "${TEST_LINKER_PATTERN}"
        "multi-config linker flags")
endif()

# Native tuning is explicit and visible when the selected compiler supports it.
if(TEST_NATIVE_SUPPORTED)
    set(native_dir "${TEST_ROOT}/native")
    run_checked(native_configure_output "native configure"
        "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${native_dir}"
        ${generator_args} ${compiler_args}
        "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
        -DBUILD_TESTS=OFF
        -DBUILD_CODEC_V2=OFF
        -DMARCH_NATIVE=ON
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
    run_checked(native_build_output "native build"
        "${CMAKE_COMMAND}" --build "${native_dir}"
        --config "${TEST_CONFIG}" --target wirehair --parallel 2)
    file(READ "${native_dir}/compile_commands.json" native_commands)
    require_match("${native_commands}" "-march=native" "native compile options")
endif()

set(pic_on_dir "${TEST_ROOT}/pic-on")
run_checked(pic_on_configure_output "PIC-on configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${pic_on_dir}"
    ${generator_args} ${compiler_args}
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    -DBUILD_SHARED_LIBS=OFF
    -DWIREHAIR_STATIC_PIC=ON
    -DBUILD_TESTS=OFF
    -DBUILD_CODEC_V2=OFF
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
run_checked(pic_on_build_output "PIC-on build"
    "${CMAKE_COMMAND}" --build "${pic_on_dir}"
    --config "${TEST_CONFIG}" --target wirehair --parallel 2)
file(READ "${pic_on_dir}/compile_commands.json" pic_on_commands)
if(TEST_PIC_FLAG_REQUIRED)
    require_match("${pic_on_commands}" "-fPIC" "default static PIC policy")
endif()

set(pic_off_dir "${TEST_ROOT}/pic-off")
run_checked(pic_off_configure_output "PIC-off configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${pic_off_dir}"
    ${generator_args} ${compiler_args}
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    -DBUILD_SHARED_LIBS=OFF
    -DWIREHAIR_STATIC_PIC=OFF
    -DBUILD_TESTS=OFF
    -DBUILD_CODEC_V2=OFF
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
run_checked(pic_off_build_output "PIC-off build"
    "${CMAKE_COMMAND}" --build "${pic_off_dir}"
    --config "${TEST_CONFIG}" --target wirehair --parallel 2)
file(READ "${pic_off_dir}/compile_commands.json" pic_off_commands)
if(TEST_PIC_FLAG_REQUIRED AND pic_off_commands MATCHES "-fPIC")
    message(FATAL_ERROR "WIREHAIR_STATIC_PIC=OFF still compiled with -fPIC")
endif()

# Explicit dual mode shares one PIC core compilation and installs both variants.
set(both_dir "${TEST_ROOT}/both")
set(both_prefix "${TEST_ROOT}/both-prefix")
run_checked(both_configure_output "build-both configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${both_dir}"
    ${generator_args} ${compiler_args}
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    "-DCMAKE_INSTALL_PREFIX=${both_prefix}"
    -DBUILD_SHARED_LIBS=OFF
    -DWIREHAIR_BUILD_BOTH=ON
    -DBUILD_TESTS=OFF
    -DBUILD_CODEC_V2=OFF)
run_checked(both_build_output "build-both build"
    "${CMAKE_COMMAND}" --build "${both_dir}"
    --config "${TEST_CONFIG}" --parallel 2)
file(GLOB_RECURSE codec_objects "${both_dir}/*WirehairCodec.cpp.o")
list(LENGTH codec_objects codec_object_count)
if(NOT codec_object_count EQUAL 1)
    message(FATAL_ERROR
        "build-both compiled WirehairCodec.cpp ${codec_object_count} times: ${codec_objects}")
endif()
run_checked(both_install_output "build-both install"
    "${CMAKE_COMMAND}" --install "${both_dir}" --config "${TEST_CONFIG}")
file(GLOB both_static "${both_prefix}/*/libwirehair.a")
file(GLOB_RECURSE both_shared
    "${both_prefix}/*/libwirehair.so*"
    "${both_prefix}/*/libwirehair.dylib*")
if(NOT both_static OR NOT both_shared)
    message(FATAL_ERROR
        "build-both artifacts missing: static=${both_static}; shared=${both_shared}")
endif()

# Offline tools exist as explicit targets but do not enter the default build.
set(tools_dir "${TEST_ROOT}/tools")
run_checked(tools_configure_output "tools configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${tools_dir}"
    ${generator_args} ${compiler_args}
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    -DBUILD_TESTS=OFF
    -DBUILD_CODEC_V2=ON
    -DWIREHAIR_BUILD_TOOLS=OFF
    -DWIREHAIR_BUILD_BENCHMARKS=OFF)
run_checked(default_build_output "default tools-excluded build"
    "${CMAKE_COMMAND}" --build "${tools_dir}"
    --config "${TEST_CONFIG}" --parallel 2)
if(TEST_MULTI_CONFIG)
    set(exe_dir "${tools_dir}/${TEST_CONFIG}")
    set(codec_exe_dir "${tools_dir}/codec/${TEST_CONFIG}")
else()
    set(exe_dir "${tools_dir}")
    set(codec_exe_dir "${tools_dir}/codec")
endif()
foreach(excluded IN ITEMS gen_small_dseeds gen_peel_seeds
        gen_most_dseeds gen_dcounts gen_tables)
    if(EXISTS "${exe_dir}/${excluded}${TEST_EXE_SUFFIX}")
        message(FATAL_ERROR "default build unexpectedly produced ${excluded}")
    endif()
endforeach()
if(EXISTS "${codec_exe_dir}/wirehair_v2_bench${TEST_EXE_SUFFIX}" OR
   EXISTS "${exe_dir}/unit_test${TEST_EXE_SUFFIX}")
    message(FATAL_ERROR "default BUILD_TESTS=OFF build produced test/benchmark code")
endif()

run_checked(explicit_build_output "explicit tool build"
    "${CMAKE_COMMAND}" --build "${tools_dir}"
    --config "${TEST_CONFIG}" --parallel 2 --target
    gen_small_dseeds gen_peel_seeds gen_most_dseeds gen_dcounts gen_tables
    wirehair_v2_bench)
run_checked(small_output "small dense seed smoke"
    "${exe_dir}/gen_small_dseeds${TEST_EXE_SUFFIX}"
    --selection-self-test)
run_checked(peel_output "peel seed smoke"
    "${exe_dir}/gen_peel_seeds${TEST_EXE_SUFFIX}"
    --trials 1 --sublo 0 --subhi 0 --nlo 2048 --nhi 2048
    --max-tries 1 --skip-tuning)
run_checked(most_output "most dense seed smoke"
    "${exe_dir}/gen_most_dseeds${TEST_EXE_SUFFIX}"
    --seed 1 --trials 1 --dense-index-lo 13 --dense-index-hi 13)
run_checked(dcount_output "dense count smoke"
    "${exe_dir}/gen_dcounts${TEST_EXE_SUFFIX}"
    --seed 1 --trials 1 --nlo 2 --nhi 2
    --max-failures 1 --low-count-run 1)
run_checked(tables_output "table generator smoke"
    "${exe_dir}/gen_tables${TEST_EXE_SUFFIX}"
    --no-benchmarks --heavy-trials 0)
set(gen_tables_exe "${exe_dir}/gen_tables${TEST_EXE_SUFFIX}")
require_cli_rejected_before_output("${gen_tables_exe}"
    "missing heavy-trials value")
foreach(invalid IN ITEMS -1 +1 " 1" "1 " 1x 0x10 10000001
        4294967296 18446744073709551616)
    require_cli_rejected_before_output("${gen_tables_exe}"
        "invalid heavy-trials value '${invalid}'" "${invalid}")
endforeach()
run_checked(tables_maximum_output "maximum heavy-trials parse"
    "${gen_tables_exe}" --heavy-trials 10000000 --help)
require_match("${tables_maximum_output}" "maximum 10000000"
    "maximum heavy-trials help")
run_checked(bench_output "v2 benchmark smoke"
    "${codec_exe_dir}/wirehair_v2_bench${TEST_EXE_SUFFIX}"
    compare --nlo 2 --nhi 2 --trials 1 --bb-list 8
    --max-message-mib 1 --loss 0)

file(READ "${PROJECT_SOURCE_DIR}/CMakeLists.txt" root_cmake)
if(root_cmake MATCHES "add_compile_options[ \t\r\n]*\\([ \t\r\n]*-w" OR
   root_cmake MATCHES "add_compile_options[ \t\r\n]*\\([ \t\r\n]*/W0")
    message(FATAL_ERROR "project contains a blanket warning suppression")
endif()
