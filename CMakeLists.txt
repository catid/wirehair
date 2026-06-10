cmake_minimum_required(VERSION 3.5)
project(wirehair)


set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(MARCH_NATIVE "Use -march=native option" ON)
option(BUILD_TESTS "build testsuite" ON)
option(BUILD_CODEC_V2 "build isolated v2 codec scaffold" ON)
set(WH_LTO "OFF" CACHE STRING "Link-time optimization mode: OFF, ON, AUTO, or THIN")
set_property(CACHE WH_LTO PROPERTY STRINGS OFF ON AUTO THIN)
set(WH_PGO_MODE "OFF" CACHE STRING "Profile-guided optimization mode: OFF, GENERATE, or USE")
set_property(CACHE WH_PGO_MODE PROPERTY STRINGS OFF GENERATE USE)
set(WH_PGO_DIR "${CMAKE_BINARY_DIR}/pgo" CACHE PATH "Profile data directory for WH_PGO_MODE")

set(LIB_SOURCE_FILES
        wirehair.cpp
        include/wirehair/wirehair.h
        gf256.cpp
        gf256.h
        WirehairCodec.cpp
        WirehairCodec.h
        WirehairTools.cpp
        WirehairTools.h
        )

set(UNIT_TEST_SOURCE_FILES
        test/SiameseTools.cpp
        test/SiameseTools.h
        test/UnitTest.cpp
        )

set(GEN_SMALL_DSEEDS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GenerateSmallDenseSeeds.cpp
        )

set(GEN_PEEL_SEEDS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GeneratePeelSeeds.cpp
        )

set(GEN_MOST_DSEEDS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GenerateMostDenseSeeds.cpp
        )

set(GEN_DCOUNTS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GenerateDenseCount.cpp
        )

set(GEN_TABLES
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/TableGenerator.cpp
        tables/HeavyRowGenerator.cpp
        tables/HeavyRowGenerator.h
        gf256.cpp
        gf256.h
        )

if(MSVC)
else()
    # Append so user-provided flags (e.g. -DCMAKE_CXX_FLAGS="-DWH_GFNI=1
    # -mgfni") survive; overwriting silently disabled opt-in SIMD kernels
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")
endif()

if(MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native")
endif()

if(NOT MSVC)
    string(TOUPPER "${WH_LTO}" WH_LTO_UC)
    set(WH_LTO_FLAG "")
    if(WH_LTO_UC STREQUAL "ON")
        set(WH_LTO_FLAG "-flto")
    elseif(WH_LTO_UC STREQUAL "AUTO")
        set(WH_LTO_FLAG "-flto=auto")
    elseif(WH_LTO_UC STREQUAL "THIN")
        set(WH_LTO_FLAG "-flto=thin")
    elseif(NOT WH_LTO_UC STREQUAL "OFF")
        message(FATAL_ERROR "Unknown WH_LTO='${WH_LTO}' (expected OFF, ON, AUTO, or THIN)")
    endif()
    if(WH_LTO_FLAG)
        add_compile_options(${WH_LTO_FLAG})
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${WH_LTO_FLAG}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${WH_LTO_FLAG}")
    endif()

    string(TOUPPER "${WH_PGO_MODE}" WH_PGO_UC)
    if(WH_PGO_UC STREQUAL "GENERATE")
        file(MAKE_DIRECTORY "${WH_PGO_DIR}")
        add_compile_options("-fprofile-generate=${WH_PGO_DIR}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-generate=${WH_PGO_DIR}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-generate=${WH_PGO_DIR}")
    elseif(WH_PGO_UC STREQUAL "USE")
        add_compile_options("-fprofile-use=${WH_PGO_DIR}" -fprofile-correction -Wno-missing-profile)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-use=${WH_PGO_DIR} -fprofile-correction -Wno-missing-profile")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-use=${WH_PGO_DIR} -fprofile-correction -Wno-missing-profile")
    elseif(NOT WH_PGO_UC STREQUAL "OFF")
        message(FATAL_ERROR "Unknown WH_PGO_MODE='${WH_PGO_MODE}' (expected OFF, GENERATE, or USE)")
    endif()
endif()

add_library(wirehair ${LIB_SOURCE_FILES})
set_target_properties(wirehair PROPERTIES VERSION 2)
set_target_properties(wirehair PROPERTIES SOVERSION 2)
target_include_directories(wirehair PUBLIC ${PROJECT_SOURCE_DIR}/include)

if (BUILD_TESTS)
    add_executable(unit_test ${UNIT_TEST_SOURCE_FILES})
    target_link_libraries(unit_test wirehair)

    add_executable(gen_small_dseeds ${GEN_SMALL_DSEEDS})
    target_link_libraries(gen_small_dseeds wirehair)

    add_executable(gen_peel_seeds ${GEN_PEEL_SEEDS})
    target_link_libraries(gen_peel_seeds wirehair)

    add_executable(gen_most_dseeds ${GEN_MOST_DSEEDS})
    target_link_libraries(gen_most_dseeds wirehair)

    add_executable(gen_dcounts ${GEN_DCOUNTS})
    target_link_libraries(gen_dcounts wirehair)

    add_executable(gen_tables ${GEN_TABLES})
endif()

if (BUILD_CODEC_V2)
    add_subdirectory(codec)
endif()

include(GNUInstallDirs)

install(TARGETS wirehair
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY python DESTINATION ${CMAKE_INSTALL_PREFIX})
