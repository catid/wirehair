# The complete production source hash is checked by CMakeLists.txt.
set(_old [=[    if (SmallShape<3>(messageBytes, blockBytes)) {
        return wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_SMALL_K3_2026_09, message, messageBytes,
            blockBytes, serializedProfileOut, serializedProfileCapacity,
            serializedProfileBytesOut, codecOut);
    }
]=])
set(_extra [=[    if (SmallShape<8>(messageBytes, blockBytes)) {
        return wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_SMALL_K8_2026_09, message, messageBytes,
            blockBytes, serializedProfileOut, serializedProfileCapacity,
            serializedProfileBytesOut, codecOut);
    }
]=])
string(REPLACE "${_old}" "" _removed "${_candidate}")
string(LENGTH "${_candidate}" _before)
string(LENGTH "${_removed}" _after)
string(LENGTH "${_old}" _expected)
math(EXPR _actual "${_before} - ${_after}")
if(NOT _actual EQUAL _expected)
    message(FATAL_ERROR "Missing or ambiguous ordinary-K3 selector anchor")
endif()
string(REPLACE "${_old}" "${_old}${_extra}" _candidate "${_candidate}")
