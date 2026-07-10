foreach(required IN ITEMS PROJECT_SOURCE_DIR TEST_ROOT TEST_GENERATOR
        TEST_C_COMPILER TEST_CXX_COMPILER TEST_SHARED TEST_BOTH TEST_CONFIG
        TEST_STRICT)
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

set(custom_libdir "lib/wirehair-e2e")
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

    file(GLOB_RECURSE python_artifacts
        LIST_DIRECTORIES FALSE
        RELATIVE "${prefix}"
        "${prefix}/python/*")
    list(SORT python_artifacts)
    if(expect_shared)
        if(NOT "${python_artifacts}" STREQUAL "python/whirehair.py")
            message(FATAL_ERROR
                "${label}: shared/dual install must contain only the runtime "
                "Python wrapper; found '${python_artifacts}'")
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

function(run_installed_python_e2e prefix label)
    if(NOT TEST_PYTHON_EXECUTABLE OR
       TEST_PYTHON_EXECUTABLE MATCHES "-NOTFOUND$")
        message(STATUS
            "${label}: Python interpreter unavailable; native binding E2E skipped")
        return()
    endif()
    find_installed_shared_library("${prefix}" library)
    run_checked("${label} automatic Python library discovery"
        "${CMAKE_COMMAND}" -E env
        "WIREHAIR_LIBRARY="
        "WIREHAIR_EXPECT_PREFIX=${prefix}"
        "PYTHONPATH=${prefix}/python"
        "${TEST_PYTHON_EXECUTABLE}" -c
        "import os, pathlib, whirehair; native = whirehair.initialize(); pathlib.Path(native._name).resolve().relative_to(pathlib.Path(os.environ['WIREHAIR_EXPECT_PREFIX']).resolve())")
    run_checked("${label} native Python E2E"
        "${TEST_PYTHON_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/ci/python_native_test.py"
        --module-dir "${prefix}/python"
        --library "${library}")
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
run_checked("producer default-component install"
    "${CMAKE_COMMAND}" --install "${producer_build}"
    --config "${TEST_CONFIG}" --prefix "${component_prefix}"
    --component Unspecified)

assert_install_manifest("${original_prefix}" "complete install")
assert_install_manifest("${component_prefix}" "default-component install")

file(RENAME "${original_prefix}" "${relocated_prefix}")
assert_install_manifest("${relocated_prefix}" "relocated install")
if(expect_shared)
    run_installed_python_e2e("${relocated_prefix}" "relocated install")
    run_installed_python_e2e("${component_prefix}" "default-component install")
endif()
run_checked("consumer configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}/test/package"
    -B "${consumer_build}" ${generator_args}
    "-DCMAKE_C_COMPILER=${TEST_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}"
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    "-DWIREHAIR_STRICT_WARNINGS=${TEST_STRICT}"
    "-Dwirehair_DIR=${relocated_prefix}/${custom_libdir}/cmake/wirehair")
run_checked("consumer build"
    "${CMAKE_COMMAND}" --build "${consumer_build}"
    --config "${TEST_CONFIG}" --parallel 2)
run_checked("consumer tests"
    "${CMAKE_CTEST_COMMAND}" --test-dir "${consumer_build}"
    -C "${TEST_CONFIG}" --output-on-failure)
