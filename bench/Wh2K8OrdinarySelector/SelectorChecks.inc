// Included by the hash-derived small public test. Reuse its allocation-failure
// hooks, packet helpers and assertions; no production code includes this file.
namespace {
uint64_t SelectedId(unsigned k)
{
    return k == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
        k == 8 ? WIREHAIR_V2_PROFILE_SMALL_K8_2026_09 : WIREHAIR_V2_PROFILE_CURRENT;
}
void NeighborRoutes()
{
    size_t cases = 0;
    for (unsigned k : {2u,3u,4u,5u,6u,7u,8u,9u,16u,64u,128u})
        for (uint32_t block : {1u,2u,64u}) for (uint32_t tail : {1u,block})
            for (unsigned policy : {0u,1u,2u}) {
                const auto source = Message(size_t(k - 1) * block + tail);
                WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
                options.source_policy = policy == 2 ? WirehairV2EncoderSource_BorrowedImmutable :
                    WirehairV2EncoderSource_Independent;
                Profile actual = {}, expected = {};
                WirehairV2Codec ordinary = nullptr, explicit_codec = nullptr;
                uint32_t bytes = 0;
                Check(Create(0, policy, source.data(), source.size(), block, actual, &ordinary) ==
                    WirehairV2_Success && ordinary, "neighbor ordinary constructor");
                Check(wirehair_v2_encoder_create_profile_id_with_options(SelectedId(k),
                    source.data(), source.size(), block, &options, expected.data(), 32, &bytes,
                    &explicit_codec) == WirehairV2_Success && explicit_codec && bytes == 32 &&
                    actual == expected, "neighbor exact selected descriptor");
                for (uint32_t id : {0u,k-1,k,k+1,1023u,1024u,UINT32_MAX}) {
                    std::vector<Byte> packet(block, 0xa5);
                    Check(wirehair_v2_encode(explicit_codec, id, packet.data(), block, &bytes) ==
                        WirehairV2_Success && bytes == (id == k - 1 ? tail : block),
                        "neighbor explicit packet");
                    packet.resize(bytes);
                    // Certified neighbors do not promise the small core's
                    // allocation-free packet path. Check bytes/guards, not a
                    // stronger allocation contract belonging only to K3/K8.
                    std::vector<Byte> ordinary_packet(size_t(block)+2,0x5a);
                    uint32_t ordinary_bytes = 0;
                    Check(wirehair_v2_encode(ordinary,id,ordinary_packet.data()+1,block,
                        &ordinary_bytes) == WirehairV2_Success && ordinary_bytes == bytes &&
                        !std::memcmp(ordinary_packet.data()+1,packet.data(),bytes) &&
                        ordinary_packet.front() == 0x5a &&
                        std::all_of(ordinary_packet.begin()+1+bytes,ordinary_packet.end(),
                            [](Byte b) { return b == 0x5a; }), "neighbor packet bytes and guards");
                }
                wirehair_v2_free(explicit_codec); wirehair_v2_free(ordinary);
                ++cases;
            }
    std::cout << "Ordinary selector neighbors: " << cases << " route/policy cases\n";
}
size_t FirstAllocation(bool ordinary, uint64_t id, unsigned policy,
    const void* source, uint64_t message, uint32_t block)
{
    Profile p; p.fill(0xa5); const auto before = p;
    WirehairV2Codec h = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = policy == 2 ? WirehairV2EncoderSource_BorrowedImmutable :
        WirehairV2EncoderSource_Independent;
    uint32_t bytes = 777;
    Start(0);
    const auto result = ordinary ? (policy ?
        wirehair_v2_encoder_create_with_options(source, message, block, &options,
            p.data(), 32, &bytes, &h) :
        wirehair_v2_encoder_create(source, message, block, p.data(), 32, &bytes, &h)) :
        wirehair_v2_encoder_create_profile_id_with_options(id, source, message, block,
            &options, p.data(), 32, &bytes, &h);
    const size_t count = Stop(), requested = first_allocation_bytes;
    Check(result == WirehairV2_OOM && !h && count == 1 && requested &&
        bytes == 32 && p == before, "first allocation fails before source read/publication");
    return requested;
}
void SelectorBounds()
{
#if defined(__unix__)
    const uint32_t bound = UINT32_C(268435456) / 9;
    const size_t span = size_t(8) * (bound + 1) + 1;
    // Reserve inaccessible virtual address space, not physical source storage.
    // First-allocation failure must prevent every source read on both routes.
    void* source = mmap(nullptr, span, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    Check(source != MAP_FAILED, "selector boundary guard mapping");
    size_t cases = 0;
    for (uint32_t block : {1u,2u,bound,bound+1})
        for (uint64_t message : {uint64_t(7)*block,uint64_t(7)*block+1,
                uint64_t(8)*block,uint64_t(8)*block+1}) for (unsigned policy : {0u,1u,2u}) {
            const bool small = block <= bound && message > uint64_t(7)*block &&
                message <= uint64_t(8)*block;
            const uint64_t id = small ? WIREHAIR_V2_PROFILE_SMALL_K8_2026_09 :
                WIREHAIR_V2_PROFILE_CURRENT;
            const size_t expected = FirstAllocation(false, id, policy, source, message, block);
            const size_t actual = FirstAllocation(true, 0, policy, source, message, block);
            Check(actual == expected, "exact K8 boundary and certified fallback allocation type");
            ++cases;
        }
    Profile p; p.fill(0xa5); const auto before = p;
    WirehairV2Codec h = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    uint32_t bytes = 777;
    Start(0);
    const auto result = wirehair_v2_encoder_create_profile_id(
        WIREHAIR_V2_PROFILE_SMALL_K8_2026_09, source, uint64_t(8)*(bound+1), bound+1,
        p.data(), 32, &bytes, &h);
    const size_t count = Stop();
    Check(result == WirehairV2_InvalidDimensions && count == 0 && !h && bytes == 32 && p == before,
        "explicit K8 rejects just beyond ordinary fallback boundary");
    Check(munmap(source, span) == 0, "release selector boundary mapping");
    std::cout << "Ordinary K8 bound/fallback probes: " << cases << " (no payload reads)\n";
#endif
}
void SelectorCppReplacement()
{
    using namespace wirehair::v2;
    Encoder encoder;
    for (unsigned k : {8u,9u,8u,3u,8u,5u,8u}) {
        const auto source = Message(k * 64 - 1);
        SerializedProfile descriptor;
        Check(encoder.Create(source.data(), source.size(), 64, descriptor) == WirehairV2_Success &&
            encoder, "ordinary independent C++ replacement across codec types");
        WirehairV2Profile parsed = {};
        Check(descriptor.Deserialize(parsed) == WirehairV2_Success && parsed.profile_id == SelectedId(k),
            "C++ replacement chooses expected profile");
        uint8_t saved[64], after[64]; uint32_t saved_bytes = 0, after_bytes = 0;
        Check(encoder.Encode(UINT32_MAX, saved, 64, saved_bytes) == WirehairV2_Success && saved_bytes == 64,
            "C++ packet before failed replacement");
        SerializedProfile failed; std::memset(failed.data(), 0xa5, failed.size());
        Start(0);
        const auto result = encoder.Create(source.data(), source.size(), 64, failed);
        Stop();
        Check(result == WirehairV2_OOM && encoder &&
            std::all_of(failed.data(), failed.data()+failed.size(), [](Byte b) { return b == 0xa5; }),
            "C++ failed replacement preserves descriptor and old handle");
        Check(encoder.Encode(UINT32_MAX, after, 64, after_bytes) == WirehairV2_Success &&
            after_bytes == saved_bytes && !std::memcmp(saved, after, 64),
            "C++ failed replacement preserves packets");
    }
}
void SelectorChecks()
{
    NeighborRoutes(); SelectorBounds(); SelectorCppReplacement();
}
} // namespace
