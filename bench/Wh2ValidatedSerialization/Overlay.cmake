# Apply only to the complete source hash authenticated by CMakeLists.txt.
# Each textual replacement must occur exactly once; no fallback or retuning.
function(replace_once old new)
    string(REPLACE "${old}" "" _removed "${_candidate}")
    string(LENGTH "${_candidate}" _before)
    string(LENGTH "${_removed}" _after)
    string(LENGTH "${old}" _expected)
    math(EXPR _actual "${_before} - ${_after}")
    if(NOT _actual EQUAL _expected)
        message(FATAL_ERROR "Missing or ambiguous validated-serialization anchor")
    endif()
    string(REPLACE "${old}" "${new}" _updated "${_candidate}")
    set(_candidate "${_updated}" PARENT_SCOPE)
endfunction()

set(_anchor [=[uint64_t BlockCountWide(uint64_t message_bytes, uint32_t block_bytes)]=])
set(_helper [=[// Only for already-validated profiles and a distinct local output array.
// Write all bytes: constructor staging arrays are deliberately uninitialized.
void SerializeValidatedProfile(
    const WirehairV2Profile& profile,
    uint8_t (&encoded)[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES])
{
    std::memset(encoded, 0, sizeof(encoded));
    std::memcpy(encoded, kProfileMagic, sizeof(kProfileMagic));
    Store16LE(encoded + 4,
        (uint16_t)WIREHAIR_V2_PROFILE_ENCODING_VERSION);
    Store16LE(encoded + 6,
        (uint16_t)WIREHAIR_V2_PROFILE_SERIALIZED_BYTES);
    Store64LE(encoded + 8, profile.profile_id);
    Store64LE(encoded + 16, profile.message_bytes);
    Store32LE(encoded + 24, profile.block_bytes);
    encoded[28] = profile.seed_attempt;
}

]=])
replace_once("${_anchor}" "${_helper}${_anchor}")

set(_old [=[    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    std::memcpy(encoded, kProfileMagic, sizeof(kProfileMagic));
    Store16LE(encoded + 4,
        (uint16_t)WIREHAIR_V2_PROFILE_ENCODING_VERSION);
    Store16LE(encoded + 6,
        (uint16_t)WIREHAIR_V2_PROFILE_SERIALIZED_BYTES);
    Store64LE(encoded + 8, profile->profile_id);
    Store64LE(encoded + 16, profile->message_bytes);
    Store32LE(encoded + 24, profile->block_bytes);
    encoded[28] = profile->seed_attempt;]=])
set(_new [=[    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    SerializeValidatedProfile(*profile, encoded);]=])
replace_once("${_old}" "${_new}")

set(_old [=[    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    uint32_t encoded_bytes = 0u;
    if (result == WirehairV2_Success) {
        result = wirehair_v2_profile_serialize(
            &profile, encoded, sizeof(encoded), &encoded_bytes);
    }
    if (result != WirehairV2_Success ||
        encoded_bytes != WIREHAIR_V2_PROFILE_SERIALIZED_BYTES)
    {
        delete codec;
        return result == WirehairV2_Success ? WirehairV2_Error : result;
    }]=])
set(_new [=[    if (result != WirehairV2_Success) {
        delete codec;
        return result;
    }
    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    SerializeValidatedProfile(profile, encoded);]=])
replace_once("${_old}" "${_new}")

set(_old [=[        uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
        uint32_t encoded_bytes = 0u;
        result = wirehair_v2_profile_serialize(
            &requested, encoded, sizeof(encoded), &encoded_bytes);
        if (result != WirehairV2_Success || encoded_bytes != sizeof(encoded)) {
            FreePublicCodec(codec);
            return result == WirehairV2_Success ? WirehairV2_Error : result;
        }]=])
set(_new [=[        uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
        SerializeValidatedProfile(requested, encoded);]=])
replace_once("${_old}" "${_new}")

set(_old [=[    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    uint32_t encoded_bytes = 0u;
    if (result == WirehairV2_Success) {
        result = wirehair_v2_profile_serialize(
            &profile, encoded, sizeof(encoded), &encoded_bytes);
    }
    if (result != WirehairV2_Success ||
        encoded_bytes != WIREHAIR_V2_PROFILE_SERIALIZED_BYTES ||
        profile.profile_id != profileId)
    {
        delete codec;
        return result == WirehairV2_Success ? WirehairV2_Error : result;
    }]=])
set(_new [=[    if (result != WirehairV2_Success || profile.profile_id != profileId) {
        delete codec;
        return result == WirehairV2_Success ? WirehairV2_Error : result;
    }
    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    SerializeValidatedProfile(profile, encoded);]=])
replace_once("${_old}" "${_new}")
