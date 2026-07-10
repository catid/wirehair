cmake_minimum_required(VERSION 3.15)
project(wirehair_package_consumer LANGUAGES C CXX)

option(WIREHAIR_STRICT_WARNINGS
    "Treat compiler warnings as errors for package consumers" OFF)

find_package(wirehair CONFIG REQUIRED)
enable_testing()

add_executable(package_c_consumer consumer.c roundtrip.c)
set_target_properties(package_c_consumer PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(package_c_consumer PRIVATE wirehair::wirehair)
add_test(NAME package_c_consumer COMMAND package_c_consumer)

add_executable(package_cpp_consumer consumer.cpp roundtrip.c)
target_link_libraries(package_cpp_consumer PRIVATE wirehair::wirehair)
add_test(NAME package_cpp_consumer COMMAND package_cpp_consumer)

add_library(package_plugin MODULE plugin.c roundtrip.c)
set_target_properties(package_plugin PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(package_plugin PRIVATE wirehair::wirehair)

add_executable(package_plugin_host plugin_host.c)
target_link_libraries(package_plugin_host PRIVATE ${CMAKE_DL_LIBS})
add_test(NAME package_plugin_round_trip
    COMMAND package_plugin_host $<TARGET_FILE:package_plugin>)

if(WIREHAIR_STRICT_WARNINGS)
    foreach(target IN ITEMS package_c_consumer package_cpp_consumer
            package_plugin package_plugin_host)
        if(MSVC)
            target_compile_options(${target} PRIVATE /W4 /WX)
        else()
            target_compile_options(${target} PRIVATE
                -Wall -Wextra -Wpedantic -Werror)
        endif()
    endforeach()
endif()

if(WIN32)
    get_target_property(wirehair_imported_type wirehair::wirehair TYPE)
    if(wirehair_imported_type STREQUAL "SHARED_LIBRARY")
        # The relocated install is the source of the DLL.  Stage it beside the
        # test hosts, matching normal application deployment on Windows.
        foreach(host IN ITEMS package_c_consumer package_cpp_consumer
                package_plugin_host)
            add_custom_command(TARGET ${host} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    $<TARGET_FILE:wirehair::wirehair>
                    $<TARGET_FILE_DIR:${host}>)
        endforeach()
    endif()
endif()
