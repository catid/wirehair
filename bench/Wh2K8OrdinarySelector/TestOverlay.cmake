# Keep the existing production test and its historical consumers unchanged.
# This explicit derivation changes only candidate routing expectations and adds
# selector-specific checks; every existing test body remains in the build.
set(_small_test_source "${_root}/test/V2SmallCodecTest.cpp")
file(SHA256 "${_small_test_source}" _test_sha)
if(NOT _test_sha STREQUAL "61041a29d5464e4f824b45cfcd6e26ddb779567a7899a07c300e10fb6159aa8e")
    message(FATAL_ERROR "Frozen small public test changed")
endif()
file(READ "${_small_test_source}" _test_candidate)
function(replace_test_once old new)
    string(REPLACE "${old}" "" _removed "${_test_candidate}")
    string(LENGTH "${_test_candidate}" _before)
    string(LENGTH "${_removed}" _after)
    string(LENGTH "${old}" _expected)
    math(EXPR _actual "${_before} - ${_after}")
    if(NOT _actual EQUAL _expected)
        message(FATAL_ERROR "Missing or ambiguous ordinary-K8 test anchor")
    endif()
    string(REPLACE "${old}" "${new}" _updated "${_test_candidate}")
    set(_test_candidate "${_updated}" PARENT_SCOPE)
endfunction()
replace_test_once("// K5/K8 are explicit-only until ordinary-path admission gates pass.\nconstexpr unsigned FirstRoute = K == 3 ? 0 : 1;"
    "// Isolated candidate selects K3/K8; production defaults are unchanged.\nconstexpr unsigned FirstRoute = (K == 3 || K == 8) ? 0 : 1;")
replace_test_once("const auto alias_result = K == 3 ?" "const auto alias_result = FirstRoute == 0 ?")
replace_test_once("const auto created = K == 3 ?" "const auto created = FirstRoute == 0 ?")
replace_test_once("const auto failed = K == 3 ?" "const auto failed = FirstRoute == 0 ?")
replace_test_once("// Explicit K5/K8 must not silently alter either ordinary selector."
    "// The isolated candidate changes K8 only; K5 remains certified.")
replace_test_once([=[host.profile_id == WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            "K5/K8 defaults stay certified until admission gates pass"]=]
    [=[host.profile_id == (k == 8 ? WIREHAIR_V2_PROFILE_SMALL_K8_2026_09 :
                WIREHAIR_V2_PROFILE_CERTIFIED_2026_07),
            "isolated ordinary K8 selection preserves K5 default"]=])
replace_test_once("int main()\n{" "#include \"SelectorChecks.inc\"\nint main()\n{")
replace_test_once("    Oracle oracle;\n    HandleAllocationIsolation();"
    "    Oracle oracle;\n    if (K == 8) SelectorChecks();\n    HandleAllocationIsolation();")
set(_small_test_generated "${_binary}/WirehairV2K8OrdinarySmallTest.cpp")
file(WRITE "${_small_test_generated}" "${_test_candidate}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_small_test_source}")

# Exercise the complete borrowed-source contract matrix at K8 as well as its
# unchanged K9 baseline. Preserve the separate certified seed/overlap fixture.
set(_borrowed_test_source "${_root}/codec/V2BorrowedSourceTest.cpp")
file(SHA256 "${_borrowed_test_source}" _borrowed_sha)
if(NOT _borrowed_sha STREQUAL "27f6d3d1915766ef8f1df0e9fa6f4628c966294c78ef07b21d32a61319c221e0")
    message(FATAL_ERROR "Frozen borrowed-source test changed")
endif()
file(READ "${_borrowed_test_source}" _test_candidate)
replace_test_once("    BlockCount = 9u," "    BlockCount = 8u,")
replace_test_once("    7u,\n    8u,\n    BlockCount," "    7u,\n    BlockCount,")
# Full input hash closes the set of identifier replacements. These are test
# calls, never changes to the public CURRENT alias or production validation.
string(REPLACE "WIREHAIR_V2_PROFILE_CURRENT" "WIREHAIR_V2_PROFILE_SMALL_K8_2026_09"
    _test_candidate "${_test_candidate}")
replace_test_once("        profile.profile_id = WIREHAIR_V2_PROFILE_SMALL_K8_2026_09;"
    "        profile.profile_id = WIREHAIR_V2_PROFILE_CERTIFIED_2026_07;")
# On narrow platforms this fixture intentionally exceeds the K8 small bound;
# keep its original certified equation contract. No 32-bit execution is claimed.
replace_test_once("    oversized.message_bytes = static_cast<uint64_t>(SIZE_MAX) + 1u;"
    "    oversized.profile_id = WIREHAIR_V2_PROFILE_CERTIFIED_2026_07;\n    oversized.message_bytes = static_cast<uint64_t>(SIZE_MAX) + 1u;")
set(_borrowed_test_generated "${_binary}/WirehairV2K8OrdinaryBorrowedTest.cpp")
file(WRITE "${_borrowed_test_generated}" "${_test_candidate}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_borrowed_test_source}")

file(SHA256 "${_root}/codec/V2ProfileTest.cpp" _profile_test_sha)
if(NOT _profile_test_sha STREQUAL "d8425d201c3b884367750c39c052850519547852a37a5b514c4ec412fb7f188d")
    message(FATAL_ERROR "Reused legacy overlap matrices changed")
endif()
set(_c_test_source "${_root}/test/V2SmallK8CConsumer.c")
file(SHA256 "${_c_test_source}" _c_test_sha)
if(NOT _c_test_sha STREQUAL "93b3475fb609db457c05bf206d07eb224afd4532774607f9b47c31948d788163")
    message(FATAL_ERROR "Frozen C consumer changed")
endif()
file(READ "${_c_test_source}" _test_candidate)
replace_test_once("wirehair_v2_encoder_create_profile_id_with_options(WIREHAIR_V2_PROFILE_SMALL_K8_2026_09,"
    "wirehair_v2_encoder_create_with_options(")
set(_c_test_generated "${_binary}/WirehairV2K8OrdinaryCConsumer.c")
file(WRITE "${_c_test_generated}" "${_test_candidate}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_c_test_source}")

set(_parity_source "${_root}/bench/Wh2V2SmallProductionParity.cpp")
file(SHA256 "${_parity_source}" _parity_source_sha)
if(NOT _parity_source_sha STREQUAL "a76777732b0dd88ed9aed9a5bc5f20f7d76be67aa493237d416ac837a9c2a342")
    message(FATAL_ERROR "Frozen installed/prototype parity adapter changed")
endif()
file(READ "${_parity_source}" _test_candidate)
replace_test_once("wirehair_v2_encoder_create_profile_id_with_options(ProfileId, source, message,"
    "wirehair_v2_encoder_create_with_options(source, message,")
set(_parity_generated "${_binary}/WirehairV2K8OrdinaryParity.cpp")
file(WRITE "${_parity_generated}" "${_test_candidate}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_parity_source}")
