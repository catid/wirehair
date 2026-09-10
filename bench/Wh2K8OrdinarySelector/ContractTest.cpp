// Reuse selected equation-independent legacy overlap matrices, not the old
// suite's hard-coded ordinary K8=CERTIFIED golden descriptor or renamed main.
#define main Wh2LegacyProfileMainNotRun
#include "V2ProfileTest.cpp"
#undef main

int main()
{
    if (!Check(wirehair_init() == Wirehair_Success, "field initialization")) return 1;
    std::vector<uint8_t> message(117);
    FillMessage(message);
    const uint8_t small[32] = {
        0x57,0x48,0x56,0x32,0x01,0x00,0x20,0x00,
        0xe0,0x0a,0x73,0x5c,0xb8,0x76,0x92,0x7a,
        0x75,0,0,0,0,0,0,0,0x10,0,0,0,0,0,0,0
    };
    if (!CheckSelectingConstructorOverlapGuards(message) ||
        !CheckDescriptorConstructorOverlapGuards(message, small)) return 1;
    WirehairV2Codec reference = nullptr;
    if (!Check(wirehair_v2_encoder_create_profile(message.data(),small,32,&reference) ==
        WirehairV2_Success && reference,"unaliasing packet reference")) return 1;
    // Full and partial descriptor/message overlap through the actual ordinary
    // constructor. The original legacy matrix separately retains certified
    // overlap coverage on the baseline; the small suite checks exact overlap.
    for (size_t offset : {size_t(0),size_t(1),message.size()-1}) {
        std::vector<uint8_t> source(message.size()+32,0x5a);
        std::memcpy(source.data(),message.data(),message.size());
        auto expected_storage = source;
        std::memcpy(expected_storage.data()+offset,small,32);
        WirehairV2Codec codec = nullptr; uint32_t bytes = 0;
        if (!Check(wirehair_v2_encoder_create(source.data(),message.size(),16,
            source.data()+offset,32,&bytes,&codec) == WirehairV2_Success && codec &&
            bytes == 32 && source == expected_storage,
            "ordinary K8 staged exact/partial descriptor overlap")) {
            wirehair_v2_free(codec); wirehair_v2_free(reference); return 1;
        }
        std::fill(source.begin(),source.end(),0xcc);
        for (uint32_t id : {0u,1u,2u,3u,4u,5u,6u,7u,8u,UINT32_MAX}) {
            uint8_t expected[16] = {}, actual[16] = {};
            uint32_t expected_bytes = 0, actual_bytes = 0;
            if (!Check(wirehair_v2_encode(reference,id,expected,16,&expected_bytes) == WirehairV2_Success &&
                wirehair_v2_encode(codec,id,actual,16,&actual_bytes) == WirehairV2_Success &&
                expected_bytes == (id == 7 ? 5u : 16u) && actual_bytes == expected_bytes &&
                !std::memcmp(actual,expected,16),"aliased constructor prepared every source block before publication")) {
                wirehair_v2_free(codec); wirehair_v2_free(reference); return 1;
            }
        }
        wirehair_v2_free(codec);
    }
    wirehair_v2_free(reference);
    std::puts("Ordinary K8 selected legacy overlap matrices passed (not the full legacy suite)");
    return 0;
}
