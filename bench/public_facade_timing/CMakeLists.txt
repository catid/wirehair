cmake_minimum_required(VERSION 3.15)

project(wirehair_wh2_v2_facade_timing_role LANGUAGES CXX)

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
   NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
    message(FATAL_ERROR "sealed facade timing roles require Linux x86-64")
endif()
get_filename_component(role_cxx_compiler
    "${CMAKE_CXX_COMPILER}" REALPATH)
if(NOT role_cxx_compiler STREQUAL
   "/usr/bin/x86_64-linux-gnu-g++-13")
    message(FATAL_ERROR "sealed facade timing roles require exact GCC 13")
endif()
if(NOT CMAKE_GENERATOR STREQUAL "Ninja" OR
   NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR "sealed facade timing roles require Ninja Release")
endif()

set(WIREHAIR_ROLE "" CACHE STRING "Facade timing role: current or parent")
set(WIREHAIR_ROLE_SOURCE_ROOT "" CACHE PATH
    "Exact implementation source root supplying the public header")
set(WIREHAIR_HARNESS_SOURCE_ROOT "" CACHE PATH
    "Exact harness source root supplying the timing worker")
set(WIREHAIR_ROLE_LIBRARY "" CACHE FILEPATH
    "Exact role-specific shared libwirehair image")
set(WIREHAIR_HARNESS_GIT_COMMIT "" CACHE STRING
    "Exact harness Git commit")
set(WIREHAIR_IMPLEMENTATION_GIT_COMMIT "" CACHE STRING
    "Exact role implementation Git commit")

foreach(required
        WIREHAIR_ROLE_SOURCE_ROOT
        WIREHAIR_HARNESS_SOURCE_ROOT
        WIREHAIR_ROLE_LIBRARY)
    if(NOT IS_ABSOLUTE "${${required}}")
        message(FATAL_ERROR "${required} must be an absolute path")
    endif()
endforeach()
if(NOT IS_DIRECTORY "${WIREHAIR_ROLE_SOURCE_ROOT}" OR
   NOT IS_DIRECTORY "${WIREHAIR_HARNESS_SOURCE_ROOT}")
    message(FATAL_ERROR "facade timing source roots must exist")
endif()
if(NOT EXISTS "${WIREHAIR_ROLE_LIBRARY}" OR
   IS_DIRECTORY "${WIREHAIR_ROLE_LIBRARY}")
    message(FATAL_ERROR "WIREHAIR_ROLE_LIBRARY must name a file")
endif()
foreach(commit WIREHAIR_HARNESS_GIT_COMMIT WIREHAIR_IMPLEMENTATION_GIT_COMMIT)
    string(LENGTH "${${commit}}" commit_length)
    if(NOT commit_length EQUAL 40 OR NOT "${${commit}}" MATCHES "^[0-9a-f]+$")
        message(FATAL_ERROR "${commit} must be one lowercase 40-hex commit")
    endif()
endforeach()

if(WIREHAIR_ROLE STREQUAL "current")
    if(NOT WIREHAIR_IMPLEMENTATION_GIT_COMMIT STREQUAL
       "d6ddb35e046956174202584fb5a26c9a79679ea8")
        message(FATAL_ERROR "current role commit differs from frozen feature")
    endif()
    set(role_definition WIREHAIR_WH2_FACADE_ROLE_CURRENT=1)
    set(role_output wirehair_wh2_v2_facade_timing_worker_current)
elseif(WIREHAIR_ROLE STREQUAL "parent")
    if(NOT WIREHAIR_IMPLEMENTATION_GIT_COMMIT STREQUAL
       "101584b7e5c30326b1791429221c331c82a00807")
        message(FATAL_ERROR "parent role commit differs from exact baseline")
    endif()
    set(role_definition WIREHAIR_WH2_FACADE_ROLE_PARENT=1)
    set(role_output wirehair_wh2_v2_facade_timing_worker_parent)
else()
    message(FATAL_ERROR "WIREHAIR_ROLE must be exactly current or parent")
endif()

set(worker_source
    "${WIREHAIR_HARNESS_SOURCE_ROOT}/bench/Wh2V2FacadeTimingWorker.cpp")
set(frozen_trace_source
    "${WIREHAIR_HARNESS_SOURCE_ROOT}/bench/Wh2FrozenTrace.cpp")
foreach(harness_source worker_source frozen_trace_source)
    if(NOT EXISTS "${${harness_source}}")
        message(FATAL_ERROR "facade timing harness source is absent: ${${harness_source}}")
    endif()
endforeach()

add_executable(wirehair_wh2_v2_facade_timing_worker
    "${worker_source}"
    "${frozen_trace_source}")

# Make the image named by the link command identical to the image resolved at
# runtime.  The external build receipt separately retains the original
# implementation-library input and proves byte equality with this regular
# role-local copy.
set(role_runtime_library "${CMAKE_CURRENT_BINARY_DIR}/libwirehair.so.2")
configure_file(
    "${WIREHAIR_ROLE_LIBRARY}"
    "${role_runtime_library}"
    COPYONLY)

set_target_properties(wirehair_wh2_v2_facade_timing_worker PROPERTIES
    OUTPUT_NAME "${role_output}"
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS ON
    BUILD_WITH_INSTALL_RPATH ON
    BUILD_RPATH "$ORIGIN"
    INSTALL_RPATH "$ORIGIN")
target_include_directories(wirehair_wh2_v2_facade_timing_worker PRIVATE
    "${WIREHAIR_ROLE_SOURCE_ROOT}/include")
target_compile_definitions(wirehair_wh2_v2_facade_timing_worker PRIVATE
    "${role_definition}"
    "WIREHAIR_WH2_HARNESS_GIT_COMMIT=\"${WIREHAIR_HARNESS_GIT_COMMIT}\""
    "WIREHAIR_WH2_IMPLEMENTATION_GIT_COMMIT=\"${WIREHAIR_IMPLEMENTATION_GIT_COMMIT}\"")
target_compile_options(wirehair_wh2_v2_facade_timing_worker PRIVATE
    -O3 -DNDEBUG -Wall -Wextra)
target_link_options(wirehair_wh2_v2_facade_timing_worker PRIVATE
    "LINKER:-z,now" "LINKER:-z,relro")
target_link_libraries(wirehair_wh2_v2_facade_timing_worker PRIVATE
    "${role_runtime_library}")
