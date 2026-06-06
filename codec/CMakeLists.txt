add_library(wirehair_v2_policy
    WirehairV2Codec.cpp
    WirehairV2Codec.h
    WirehairV2Peel.cpp
    WirehairV2Peel.h
    WirehairV2Policy.cpp
    WirehairV2Policy.h
    WirehairV2Seeds.cpp
    WirehairV2Seeds.h
)

target_include_directories(wirehair_v2_policy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(wirehair_v2_policy PUBLIC wirehair)

if (BUILD_TESTS)
    add_executable(wirehair_v2_policy_test PolicyTest.cpp)
    target_link_libraries(wirehair_v2_policy_test wirehair_v2_policy)

    add_executable(wirehair_v2_bench WirehairV2Bench.cpp)
    target_link_libraries(wirehair_v2_bench wirehair_v2_policy)
endif()
