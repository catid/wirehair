foreach(required IN ITEMS PROJECT_SOURCE_DIR TEST_ROOT TEST_GENERATOR
        TEST_C_COMPILER TEST_CXX_COMPILER TEST_SHARED TEST_CONFIG TEST_STRICT)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

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
set(relocated_prefix "${TEST_ROOT}/relocated-prefix")
set(consumer_build "${TEST_ROOT}/consumer-build")
file(REMOVE_RECURSE "${TEST_ROOT}")
file(MAKE_DIRECTORY "${TEST_ROOT}")

run_checked("producer configure"
    "${CMAKE_COMMAND}" -S "${PROJECT_SOURCE_DIR}" -B "${producer_build}"
    ${generator_args}
    "-DCMAKE_C_COMPILER=${TEST_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}"
    "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}"
    "-DCMAKE_INSTALL_PREFIX=${original_prefix}"
    "-DCMAKE_INSTALL_LIBDIR=${custom_libdir}"
    "-DBUILD_SHARED_LIBS=${TEST_SHARED}"
    "-DWIREHAIR_STRICT_WARNINGS=${TEST_STRICT}"
    -DBUILD_TESTS=OFF
    -DBUILD_CODEC_V2=OFF
    -DWIREHAIR_BUILD_BOTH=OFF
    -DMARCH_NATIVE=OFF)
run_checked("producer build"
    "${CMAKE_COMMAND}" --build "${producer_build}"
    --config "${TEST_CONFIG}" --parallel 2)
run_checked("producer install"
    "${CMAKE_COMMAND}" --install "${producer_build}"
    --config "${TEST_CONFIG}")

if(NOT EXISTS
   "${original_prefix}/${custom_libdir}/cmake/wirehair/wirehairConfig.cmake")
    message(FATAL_ERROR "package config was not installed to the custom libdir")
endif()
if(NOT EXISTS "${original_prefix}/python/whirehair.py" OR
   EXISTS "${original_prefix}/python/test_whirehair.py")
    message(FATAL_ERROR "install must contain only the runtime Python wrapper")
endif()
if(TEST_SHARED)
    if(UNIX AND NOT APPLE)
        file(GLOB shared_artifacts
            "${original_prefix}/${custom_libdir}/libwirehair.so*")
        file(GLOB static_artifacts
            "${original_prefix}/${custom_libdir}/libwirehair.a")
        if(NOT shared_artifacts OR static_artifacts)
            message(FATAL_ERROR
                "shared-only artifact mismatch: shared=${shared_artifacts}; static=${static_artifacts}")
        endif()
    endif()
else()
    if(UNIX)
        file(GLOB static_artifacts
            "${original_prefix}/${custom_libdir}/libwirehair.a")
        file(GLOB shared_artifacts
            "${original_prefix}/${custom_libdir}/libwirehair.so*"
            "${original_prefix}/${custom_libdir}/libwirehair.dylib*")
        if(NOT static_artifacts OR shared_artifacts)
            message(FATAL_ERROR
                "static-only artifact mismatch: static=${static_artifacts}; shared=${shared_artifacts}")
        endif()
    endif()
endif()

file(RENAME "${original_prefix}" "${relocated_prefix}")
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
