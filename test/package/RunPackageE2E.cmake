foreach(required IN ITEMS PROJECT_SOURCE_DIR TEST_ROOT TEST_GENERATOR
        TEST_C_COMPILER TEST_CXX_COMPILER TEST_SHARED TEST_BOTH TEST_CONFIG
        TEST_STRICT TEST_PROJECT_VERSION)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

if(TEST_SHARED AND TEST_BOTH)
    message(FATAL_ERROR "TEST_SHARED and TEST_BOTH cannot both be enabled")
endif()

function(run_checked label)
    execute_process(
        COMMAND ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR
            "${label} failed (${result})\nstdout:\n${out}\nstderr:\n${err}")
    endif()
endfunction()

function(run_capture output_var label)
    execute_process(
        COMMAND ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR
            "${label} failed (${result})\nstdout:\n${out}\nstderr:\n${err}")
    endif()
    string(STRIP "${out}" out)
    set(${output_var} "${out}" PARENT_SCOPE)
endfunction()

if(TEST_SHARED OR TEST_BOTH)
    set(expect_shared ON)
else()
    set(expect_shared OFF)
endif()
if(NOT TEST_SHARED OR TEST_BOTH)
    set(expect_static ON)
else()
    set(expect_static OFF)
endif()

set(generator_args -G "${TEST_GENERATOR}")
if(DEFINED TEST_GENERATOR_PLATFORM AND NOT TEST_GENERATOR_PLATFORM STREQUAL "")
    list(APPEND generator_args -A "${TEST_GENERATOR_PLATFORM}")
endif()
if(DEFINED TEST_GENERATOR_TOOLSET AND NOT TEST_GENERATOR_TOOLSET STREQUAL "")
    list(APPEND generator_args -T "${TEST_GENERATOR_TOOLSET}")
endif()

# Exercise a nested multilib-style directory rather than assuming plain lib.
set(custom_libdir "lib64/wirehair-e2e")
set(producer_build "${TEST_ROOT}/producer-build")
set(original_prefix "${TEST_ROOT}/original-prefix")
set(component_prefix "${TEST_ROOT}/component-prefix")
set(relocated_prefix "${TEST_ROOT}/relocated-prefix")
set(consumer_build "${TEST_ROOT}/consumer-build")
file(REMOVE_RECURSE "${TEST_ROOT}")
file(MAKE_DIRECTORY "${TEST_ROOT}")

function(assert_install_manifest prefix label)
    if(NOT EXISTS
       "${prefix}/${custom_libdir}/cmake/wirehair/wirehairConfig.cmake")
        message(FATAL_ERROR
            "${label}: package config was not installed to the custom libdir")
    endif()
    if(NOT EXISTS "${prefix}/include/wirehair/wirehair.h")
        message(FATAL_ERROR "${label}: public header was not installed")
    endif()
    if(NOT EXISTS "${prefix}/include/wirehair/wirehair.hpp")
        message(FATAL_ERROR "${label}: public C++ header was not installed")
    endif()
    set(pc_file "${prefix}/${custom_libdir}/pkgconfig/wirehair.pc")
    if(NOT EXISTS "${pc_file}")
        message(FATAL_ERROR "${label}: relocatable wirehair.pc was not installed")
    endif()
    file(READ "${pc_file}" pc_contents)
    if(NOT pc_contents MATCHES "(^|\n)Version: ${TEST_PROJECT_VERSION}(\n|$)")
        message(FATAL_ERROR
            "${label}: pkg-config version does not match ${TEST_PROJECT_VERSION}")
    endif()
    set(installed_license "${prefix}/share/licenses/wirehair/LICENSE")
    if(NOT EXISTS "${installed_license}")
        message(FATAL_ERROR "${label}: installed license text is missing")
    endif()
    file(READ "${PROJECT_SOURCE_DIR}/LICENSE" source_license_text)
    file(READ "${installed_license}" installed_license_text)
    if(NOT installed_license_text STREQUAL source_license_text)
        message(FATAL_ERROR "${label}: installed license text was altered")
    endif()
    foreach(contract IN ITEMS LEGACY_WIRE_PROFILES.md V2_WIRE_PROFILE.md)
        set(installed_contract
            "${prefix}/share/doc/wirehair/${contract}")
        if(NOT EXISTS "${installed_contract}")
            message(FATAL_ERROR
                "${label}: installed wire contract is missing: ${contract}")
        endif()
        file(READ "${PROJECT_SOURCE_DIR}/${contract}" source_contract_text)
        file(READ "${installed_contract}" installed_contract_text)
        if(NOT installed_contract_text STREQUAL source_contract_text)
            message(FATAL_ERROR
                "${label}: installed wire contract was altered: ${contract}")
        endif()
    endforeach()

    file(GLOB_RECURSE python_artifacts
        LIST_DIRECTORIES FALSE
        RELATIVE "${prefix}"
        "${prefix}/python/*")
    list(SORT python_artifacts)
    if(expect_shared)
        set(expected_python_artifacts
            "python/whirehair.py;python/wirehair/__init__.py")
        if(NOT "${python_artifacts}" STREQUAL
           "${expected_python_artifacts}")
            message(FATAL_ERROR
                "${label}: shared/dual install must contain both standard and "
                "compatibility Python imports; found '${python_artifacts}'")
        endif()
    elseif(python_artifacts)
        message(FATAL_ERROR
            "${label}: static-only install contains unusable Python artifacts: "
            "${python_artifacts}")
    endif()

    if(WIN32)
        file(GLOB shared_artifacts "${prefix}/bin/*wirehair*.dll")
        file(GLOB archive_artifacts
            "${prefix}/${custom_libdir}/*wirehair*.lib"
            "${prefix}/${custom_libdir}/*wirehair*.a")
        if(NOT archive_artifacts)
            message(FATAL_ERROR
                "${label}: no static archive or DLL import library was installed")
        endif()
        if(TEST_BOTH)
            list(LENGTH archive_artifacts archive_count)
            if(archive_count LESS 2)
                message(FATAL_ERROR
                    "${label}: dual install has fewer than two archive/import "
                    "artifacts: ${archive_artifacts}")
            endif()
        endif()
    elseif(UNIX)
        file(GLOB static_artifacts
            "${prefix}/${custom_libdir}/libwirehair.a")
        file(GLOB shared_artifacts
            "${prefix}/${custom_libdir}/libwirehair.so*"
            "${prefix}/${custom_libdir}/libwirehair.dylib*")
        if(expect_static AND NOT static_artifacts)
            message(FATAL_ERROR
                "${label}: expected static archive was not installed")
        elseif(NOT expect_static AND static_artifacts)
            message(FATAL_ERROR
                "${label}: unexpected static archive: ${static_artifacts}")
        endif()
    else()
        message(FATAL_ERROR
            "${label}: package artifact test does not support this platform")
    endif()

    if(expect_shared AND NOT shared_artifacts)
        message(FATAL_ERROR
            "${label}: expected shared library was not installed")
    elseif(NOT expect_shared AND shared_artifacts)
        message(FATAL_ERROR
            "${label}: unexpected shared library: ${shared_artifacts}")
    endif()
endfunction()

function(assert_metadata_in_install_manifest prefix label)
    file(GLOB manifests "${producer_build}/install_manifest*.txt")
    if(NOT manifests)
        message(FATAL_ERROR "${label}: CMake install manifest is missing")
    endif()
    set(expected
        "${prefix}/${custom_libdir}/pkgconfig/wirehair.pc"
        "${prefix}/share/licenses/wirehair/LICENSE"
        "${prefix}/share/doc/wirehair/LEGACY_WIRE_PROFILES.md"
        "${prefix}/share/doc/wirehair/V2_WIRE_PROFILE.md"
        "${prefix}/include/wirehair/wirehair.h"
        "${prefix}/include/wirehair/wirehair.hpp")
    if(expect_shared)
        list(APPEND expected
            "${prefix}/python/whirehair.py"
            "${prefix}/python/wirehair/__init__.py")
    endif()
    foreach(path IN LISTS expected)
        set(found OFF)
        foreach(manifest IN LISTS manifests)
            file(STRINGS "${manifest}" installed_paths)
            list(FIND installed_paths "${path}" installed_index)
            if(NOT installed_index EQUAL -1)
                set(found ON)
                break()
            endif()
        endforeach()
        if(NOT found)
            message(FATAL_ERROR
                "${label}: uninstall manifest does not contain ${path}")
        endif()
    endforeach()
endfunction()

function(find_installed_shared_library prefix output_var)
    if(WIN32)
        file(GLOB candidates "${prefix}/bin/*wirehair*.dll")
    elseif(APPLE)
        file(GLOB candidates
            "${prefix}/${custom_libdir}/libwirehair.dylib*")
    else()
        file(GLOB candidates
            "${prefix}/${custom_libdir}/libwirehair.so*")
    endif()
    list(SORT candidates)
    if(NOT candidates)
        message(FATAL_ERROR
            "no installed shared library found under ${prefix}")
    endif()
    list(GET candidates 0 library)
    set(${output_var} "${library}" PARENT_SCOPE)
endfunction()

function(assert_installed_elf_exports prefix label)
    if(NOT expect_shared OR
       NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
        return()
    endif()
    if(NOT DEFINED TEST_NM OR TEST_NM STREQUAL "")
        message(FATAL_ERROR
            "${label}: TEST_NM is required for ELF export verification")
    endif()
    find_installed_shared_library("${prefix}" library)
    run_checked("${label} ELF public exports"
        "${CMAKE_COMMAND}"
        "-DLIBRARY=${library}"
        "-DMANIFEST=${PROJECT_SOURCE_DIR}/abi/wirehair.map"
        "-DNM=${TEST_NM}"
        -P "${PROJECT_SOURCE_DIR}/test/cmake/CheckElfExports.cmake")
endfunction()

function(run_installed_python_e2e prefix label)
    if(NOT TEST_PYTHON_EXECUTABLE OR
       TEST_PYTHON_EXECUTABLE MATCHES "-NOTFOUND$")
        message(FATAL_ERROR
            "${label}: Python is required for shared-package validation")
    endif()
    find_installed_shared_library("${prefix}" library)
    run_checked("${label} automatic Python library discovery"
        "${CMAKE_COMMAND}" -E env
        "WIREHAIR_LIBRARY="
        "WIREHAIR_PREFIX="
        "WIREHAIR_EXPECT_PREFIX=${prefix}"
        "PYTHONPATH=${prefix}/python"
        "${TEST_PYTHON_EXECUTABLE}" -c
        "import os, pathlib, wirehair, whirehair; assert wirehair.Encoder is whirehair.Encoder; assert wirehair.__version__ == '${TEST_PROJECT_VERSION}'; native = wirehair.initialize(); pathlib.Path(native._name).resolve().relative_to(pathlib.Path(os.environ['WIREHAIR_EXPECT_PREFIX']).resolve())")
    run_checked("${label} native Python E2E"
        "${TEST_PYTHON_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/ci/python_native_test.py"
        --module-dir "${prefix}/python"
        --library "${library}")
endfunction()

function(run_pkg_config_consumer prefix label mode)
    if(NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
        return()
    endif()
    if(NOT TEST_PKG_CONFIG OR TEST_PKG_CONFIG MATCHES "-NOTFOUND$")
        message(FATAL_ERROR
            "${label}: pkg-config is required for Linux package validation")
    endif()

    set(pc_dir "${prefix}/${custom_libdir}/pkgconfig")
    set(pkg_env
        "PKG_CONFIG_PATH=${pc_dir}"
        "PKG_CONFIG_LIBDIR=${pc_dir}")
    set(query_args "")
    if(mode STREQUAL "static")
        list(APPEND query_args --static)
    elseif(NOT mode STREQUAL "shared")
        message(FATAL_ERROR "Unknown pkg-config consumer mode: ${mode}")
    endif()

    run_capture(pc_version "${label} pkg-config version"
        "${CMAKE_COMMAND}" -E env ${pkg_env}
        "${TEST_PKG_CONFIG}" --modversion wirehair)
    if(NOT "${pc_version}" STREQUAL "${TEST_PROJECT_VERSION}")
        message(FATAL_ERROR
            "${label}: pkg-config reports ${pc_version}, expected "
            "${TEST_PROJECT_VERSION}")
    endif()
    run_capture(pc_prefix "${label} pkg-config relocated prefix"
        "${CMAKE_COMMAND}" -E env ${pkg_env}
        "${TEST_PKG_CONFIG}" --variable=prefix wirehair)
    get_filename_component(pc_prefix_real "${pc_prefix}" REALPATH)
    get_filename_component(expected_prefix_real "${prefix}" REALPATH)
    if(NOT "${pc_prefix_real}" STREQUAL "${expected_prefix_real}")
        message(FATAL_ERROR
            "${label}: pkg-config prefix ${pc_prefix_real} does not match "
            "relocated install ${expected_prefix_real}")
    endif()

    run_capture(cflags "${label} pkg-config cflags"
        "${CMAKE_COMMAND}" -E env ${pkg_env}
        "${TEST_PKG_CONFIG}" ${query_args} --cflags wirehair)
    run_capture(libs "${label} pkg-config libs"
        "${CMAKE_COMMAND}" -E env ${pkg_env}
        "${TEST_PKG_CONFIG}" ${query_args} --libs wirehair)
    if(expect_shared AND mode STREQUAL "shared")
        if(libs MATCHES "(^|[ \\t])-lm($|[ \\t])")
            message(FATAL_ERROR
                "${label}: shared pkg-config query overlinks libm: ${libs}")
        endif()
    elseif(NOT libs MATCHES "(^|[ \\t])-lm($|[ \\t])")
        message(FATAL_ERROR
            "${label}: static pkg-config query omits required libm: ${libs}")
    endif()
    separate_arguments(cflag_list UNIX_COMMAND "${cflags}")
    separate_arguments(lib_list UNIX_COMMAND "${libs}")

    set(executable "${TEST_ROOT}/pkg-config-${mode}")
    set(strict_flags "")
    if(TEST_STRICT)
        list(APPEND strict_flags -Wall -Wextra -Wpedantic -Werror)
    endif()
    run_checked("${label} pkg-config pure-C ${mode} compile"
        "${TEST_C_COMPILER}" ${strict_flags} ${cflag_list}
        "${PROJECT_SOURCE_DIR}/test/package/consumer.c"
        "${PROJECT_SOURCE_DIR}/test/package/roundtrip.c"
        -o "${executable}" ${lib_list})
    run_checked("${label} pkg-config pure-C ${mode} run"
        "${CMAKE_COMMAND}" -E env
        "LD_LIBRARY_PATH=${prefix}/${custom_libdir}"
        "${executable}")
endfunction()

run_checked("producer configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${producer_build}"
    ${generator_args}
    "-DCMAKE_C_COMPILER=${TEST_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}"
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    "-DCMAKE_INSTALL_PREFIX=${original_prefix}"
    "-DCMAKE_INSTALL_LIBDIR=${custom_libdir}"
    "-DBUILD_SHARED_LIBS=${TEST_SHARED}"
    "-DWIREHAIR_BUILD_BOTH=${TEST_BOTH}"
    "-DWIREHAIR_STRICT_WARNINGS=${TEST_STRICT}"
    -DBUILD_TESTS=OFF
    -DBUILD_CODEC_V2=OFF
    -DMARCH_NATIVE=OFF)
run_checked("producer build"
    "${CMAKE_COMMAND}" --build "${producer_build}"
    --config "${TEST_CONFIG}" --parallel 2)
run_checked("producer install"
    "${CMAKE_COMMAND}" --install "${producer_build}"
    --config "${TEST_CONFIG}")
assert_metadata_in_install_manifest("${original_prefix}" "complete install")
run_checked("producer default-component install"
    "${CMAKE_COMMAND}" --install "${producer_build}"
    --config "${TEST_CONFIG}" --prefix "${component_prefix}"
    --component Unspecified)
assert_metadata_in_install_manifest(
    "${component_prefix}" "default-component install")

assert_install_manifest("${original_prefix}" "complete install")
assert_install_manifest("${component_prefix}" "default-component install")
assert_installed_elf_exports("${original_prefix}" "complete install")
assert_installed_elf_exports("${component_prefix}" "default-component install")

file(RENAME "${original_prefix}" "${relocated_prefix}")
assert_install_manifest("${relocated_prefix}" "relocated install")
assert_installed_elf_exports("${relocated_prefix}" "relocated install")
if(expect_shared)
    run_installed_python_e2e("${relocated_prefix}" "relocated install")
    run_installed_python_e2e("${component_prefix}" "default-component install")
endif()
run_pkg_config_consumer("${relocated_prefix}" "relocated install" shared)
run_pkg_config_consumer("${relocated_prefix}" "relocated install" static)
run_checked("consumer configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}/test/package"
    -B "${consumer_build}" ${generator_args}
    "-DCMAKE_C_COMPILER=${TEST_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}"
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    "-DWIREHAIR_STRICT_WARNINGS=${TEST_STRICT}"
    "-DEXPECTED_WIREHAIR_VERSION=${TEST_PROJECT_VERSION}"
    "-Dwirehair_DIR=${relocated_prefix}/${custom_libdir}/cmake/wirehair")
run_checked("consumer build"
    "${CMAKE_COMMAND}" --build "${consumer_build}"
    --config "${TEST_CONFIG}" --parallel 2)
run_checked("consumer tests"
    "${CMAKE_CTEST_COMMAND}" --test-dir "${consumer_build}"
    -C "${TEST_CONFIG}" --output-on-failure)
