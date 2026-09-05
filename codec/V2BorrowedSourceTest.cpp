#include <wirehair/wirehair.h>
#include <wirehair/wirehair.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <new>
#include <utility>
#include <vector>

namespace {

static_assert(sizeof(WirehairV2EncoderOptions) == 16,
    "borrowed-source options ABI size changed");
static_assert(offsetof(WirehairV2EncoderOptions, struct_bytes) == 0,
    "borrowed-source options struct_bytes offset changed");
static_assert(offsetof(WirehairV2EncoderOptions, options_version) == 4,
    "borrowed-source options options_version offset changed");
static_assert(offsetof(WirehairV2EncoderOptions, source_policy) == 8,
    "borrowed-source options source_policy offset changed");
static_assert(offsetof(WirehairV2EncoderOptions, reserved) == 12,
    "borrowed-source options reserved offset changed");
static_assert(sizeof(WirehairV2EncoderSourcePolicy) == 4,
    "borrowed-source policy enum ABI size changed");
static_assert(WirehairV2EncoderSource_Invalid == 0,
    "borrowed-source invalid policy value changed");
static_assert(WirehairV2EncoderSource_Independent == 1,
    "borrowed-source independent policy value changed");
static_assert(WirehairV2EncoderSource_BorrowedImmutable == 2,
    "borrowed-source borrowed policy value changed");
static_assert(WirehairV2EncoderSource_Count == 3,
    "borrowed-source policy count value changed");
static_assert(WirehairV2EncoderSource_Padding == 0x7fffffff,
    "borrowed-source policy padding value changed");

enum : uint32_t
{
    BlockBytes = 16u,
    BlockCount = 9u,
    ProfileBytes = WIREHAIR_V2_PROFILE_SERIALIZED_BYTES,
    RepairProbeCount = 6u,
    ProbeCount = BlockCount + RepairProbeCount
};

enum class ConstructorRoute : unsigned
{
    DefaultSelector,
    ExplicitSelector,
    Descriptor,
    Count
};

enum class SourceFlavor : unsigned
{
    Legacy,
    Independent,
    Borrowed,
    Count
};

const uint32_t PacketIds[ProbeCount] = {
    0u,
    1u,
    2u,
    3u,
    4u,
    5u,
    6u,
    7u,
    8u,
    BlockCount,
    BlockCount + 7u,
    2u * BlockCount - 1u,
    65536u + BlockCount - 1u,
    UINT32_C(0x10000000), // B=16: byte offset is 2^32, not source block 0.
    UINT32_MAX
};

const char* RouteName(ConstructorRoute route)
{
    switch (route)
    {
    case ConstructorRoute::DefaultSelector: return "default-selector";
    case ConstructorRoute::ExplicitSelector: return "explicit-selector";
    case ConstructorRoute::Descriptor: return "descriptor";
    default: return "unknown-route";
    }
}

const char* FlavorName(SourceFlavor flavor)
{
    switch (flavor)
    {
    case SourceFlavor::Legacy: return "legacy";
    case SourceFlavor::Independent: return "independent";
    case SourceFlavor::Borrowed: return "borrowed";
    default: return "unknown-flavor";
    }
}

bool Check(bool condition, const char* what)
{
    if (condition) {
        return true;
    }
    std::fprintf(stderr, "V2 borrowed-source test failed: %s\n", what);
    return false;
}

void FillMessage(std::vector<uint8_t>& message, uint8_t salt = 0u)
{
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = static_cast<uint8_t>(i * 73u + 19u + salt);
    }
}

void PoisonAndRelease(
    std::vector<uint8_t>& source,
    const std::vector<uint8_t>& original)
{
    volatile uint8_t* const bytes = source.data();
    for (size_t i = 0; i < source.size(); ++i) {
        bytes[i] = static_cast<uint8_t>(original[i] ^ 0xffu);
    }
    std::vector<uint8_t>().swap(source);
}

WirehairV2EncoderOptions MakeOptions(uint32_t policy)
{
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = policy;
    return options;
}

struct CodecOwner
{
    WirehairV2Codec Handle = nullptr;

    ~CodecOwner()
    {
        wirehair_v2_free(Handle);
    }

    CodecOwner(const CodecOwner&) = delete;
    CodecOwner& operator=(const CodecOwner&) = delete;
    CodecOwner() = default;
};

struct PacketCapture
{
    uint32_t Id = 0u;
    uint32_t Bytes = 0u;
    std::array<uint8_t, BlockBytes> Data{};
};

using PacketSet = std::array<PacketCapture, ProbeCount>;

bool PacketSetsEqual(const PacketSet& a, const PacketSet& b)
{
    for (unsigned i = 0u; i < ProbeCount; ++i) {
        if (a[i].Id != b[i].Id || a[i].Bytes != b[i].Bytes ||
            std::memcmp(a[i].Data.data(), b[i].Data.data(), BlockBytes) != 0)
        {
            return false;
        }
    }
    return true;
}

WirehairV2Result CreateEncoder(
    ConstructorRoute route,
    SourceFlavor flavor,
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint8_t* descriptor,
    uint32_t* descriptor_bytes,
    WirehairV2Codec* codec_out)
{
    if (flavor == SourceFlavor::Legacy)
    {
        switch (route)
        {
        case ConstructorRoute::DefaultSelector:
            return wirehair_v2_encoder_create(
                message, message_bytes, block_bytes,
                descriptor, ProfileBytes, descriptor_bytes, codec_out);
        case ConstructorRoute::ExplicitSelector:
            return wirehair_v2_encoder_create_profile_id(
                WIREHAIR_V2_PROFILE_CURRENT,
                message, message_bytes, block_bytes,
                descriptor, ProfileBytes, descriptor_bytes, codec_out);
        case ConstructorRoute::Descriptor:
            return wirehair_v2_encoder_create_profile(
                message, descriptor, ProfileBytes, codec_out);
        default:
            return WirehairV2_Error;
        }
    }

    const uint32_t policy = flavor == SourceFlavor::Independent ?
        static_cast<uint32_t>(WirehairV2EncoderSource_Independent) :
        static_cast<uint32_t>(WirehairV2EncoderSource_BorrowedImmutable);
    WirehairV2EncoderOptions options = MakeOptions(policy);
    switch (route)
    {
    case ConstructorRoute::DefaultSelector:
        return wirehair_v2_encoder_create_with_options(
            message, message_bytes, block_bytes, &options,
            descriptor, ProfileBytes, descriptor_bytes, codec_out);
    case ConstructorRoute::ExplicitSelector:
        return wirehair_v2_encoder_create_profile_id_with_options(
            WIREHAIR_V2_PROFILE_CURRENT,
            message, message_bytes, block_bytes, &options,
            descriptor, ProfileBytes, descriptor_bytes, codec_out);
    case ConstructorRoute::Descriptor:
        return wirehair_v2_encoder_create_profile_with_options(
            message, descriptor, ProfileBytes, &options, codec_out);
    default:
        return WirehairV2_Error;
    }
}

bool CapturePackets(
    WirehairV2Codec codec,
    uint64_t message_bytes,
    const std::vector<uint8_t>& original,
    PacketSet& packets,
    const char* what)
{
    static const uint8_t Guard = 0x5au;
    for (unsigned probe = 0u; probe < ProbeCount; ++probe)
    {
        PacketCapture& packet = packets[probe];
        packet.Id = PacketIds[probe];
        const uint32_t expected_bytes = packet.Id == BlockCount - 1u ?
            static_cast<uint32_t>(message_bytes -
                static_cast<uint64_t>(BlockCount - 1u) * BlockBytes) :
            BlockBytes;
        uint8_t guarded[BlockBytes + 2u];
        std::memset(guarded, Guard, sizeof(guarded));
        uint32_t bytes = UINT32_C(0xa5a5a5a5);
        const WirehairV2Result result = wirehair_v2_encode(
            codec, packet.Id, guarded + 1u, BlockBytes, &bytes);
        bool bounds_ok = guarded[0] == Guard &&
            guarded[sizeof(guarded) - 1u] == Guard;
        for (uint32_t i = bytes; bounds_ok && i < BlockBytes; ++i) {
            bounds_ok = guarded[1u + i] == Guard;
        }
        const bool systematic_ok = packet.Id >= BlockCount ||
            (static_cast<uint64_t>(packet.Id) * BlockBytes + bytes <=
                original.size() &&
             std::memcmp(
                 guarded + 1u,
                 original.data() + static_cast<size_t>(packet.Id) * BlockBytes,
                 bytes) == 0);
        if (result != WirehairV2_Success || bytes != expected_bytes ||
            !bounds_ok || !systematic_ok)
        {
            std::fprintf(stderr,
                "V2 borrowed-source packet capture failed: %s "
                "probe=%u id=%u result=%d bytes=%u expected=%u\n",
                what, probe, packet.Id, static_cast<int>(result),
                bytes, expected_bytes);
            return false;
        }
        packet.Bytes = bytes;
        std::copy(guarded + 1u, guarded + 1u + bytes, packet.Data.begin());
    }
    return true;
}

bool ReplayPackets(
    WirehairV2Codec codec,
    const PacketSet& expected,
    const char* what)
{
    static const uint8_t Guard = 0x5au;
    for (unsigned probe = 0u; probe < ProbeCount; ++probe)
    {
        uint8_t guarded[BlockBytes + 2u];
        std::memset(guarded, Guard, sizeof(guarded));
        uint32_t bytes = UINT32_C(0xa5a5a5a5);
        // Exercise the public length calculation before delegated encoding,
        // including a one-byte tail (zero capacity) and high repair IDs.
        const WirehairV2Result short_result = wirehair_v2_encode(
            codec, expected[probe].Id,
            guarded + 1u, expected[probe].Bytes - 1u, &bytes);
        const bool untouched = std::all_of(
            guarded, guarded + sizeof(guarded),
            [](uint8_t byte) { return byte == Guard; });
        if (short_result != WirehairV2_BufferTooSmall ||
            bytes != expected[probe].Bytes || !untouched)
        {
            std::fprintf(stderr,
                "V2 borrowed-source short packet failed: %s "
                "probe=%u id=%u result=%d bytes=%u expected=%u\n",
                what, probe, expected[probe].Id,
                static_cast<int>(short_result), bytes, expected[probe].Bytes);
            return false;
        }
        bytes = UINT32_C(0xa5a5a5a5);
        const WirehairV2Result result = wirehair_v2_encode(
            codec, expected[probe].Id,
            guarded + 1u, expected[probe].Bytes, &bytes);
        bool bounds_ok = guarded[0] == Guard &&
            guarded[sizeof(guarded) - 1u] == Guard;
        for (uint32_t i = bytes; bounds_ok && i < BlockBytes; ++i) {
            bounds_ok = guarded[1u + i] == Guard;
        }
        if (result != WirehairV2_Success ||
            bytes != expected[probe].Bytes || !bounds_ok ||
            std::memcmp(
                guarded + 1u, expected[probe].Data.data(), bytes) != 0)
        {
            std::fprintf(stderr,
                "V2 borrowed-source packet replay failed: %s "
                "probe=%u id=%u result=%d bytes=%u expected=%u\n",
                what, probe, expected[probe].Id,
                static_cast<int>(result), bytes, expected[probe].Bytes);
            return false;
        }
    }
    return true;
}

bool BuildReference(
    const std::vector<uint8_t>& original,
    std::array<uint8_t, ProfileBytes>& descriptor,
    PacketSet& packets)
{
    CodecOwner owner;
    uint32_t descriptor_bytes = 0u;
    const WirehairV2Result result = wirehair_v2_encoder_create(
        original.data(), original.size(), BlockBytes,
        descriptor.data(), descriptor.size(),
        &descriptor_bytes, &owner.Handle);
    WirehairV2Profile parsed = {};
    return result == WirehairV2_Success && owner.Handle &&
        descriptor_bytes == ProfileBytes &&
        wirehair_v2_profile_deserialize(
            descriptor.data(), descriptor.size(), &parsed) ==
                WirehairV2_Success &&
        parsed.message_bytes == original.size() &&
        parsed.block_bytes == BlockBytes &&
        CapturePackets(
            owner.Handle, original.size(), original,
            packets, "reference");
}

bool CheckConstructorAndLifetimeMatrix()
{
    static const uint64_t MessageShapes[] = {
        static_cast<uint64_t>(BlockCount - 1u) * BlockBytes + 1u,
        static_cast<uint64_t>(BlockCount - 1u) * BlockBytes + BlockBytes - 1u,
        static_cast<uint64_t>(BlockCount) * BlockBytes
    };
    static const uint8_t Guard = 0x5au;

    for (uint64_t message_bytes : MessageShapes)
    {
        std::vector<uint8_t> original(static_cast<size_t>(message_bytes));
        FillMessage(original);
        std::array<uint8_t, ProfileBytes> expected_descriptor{};
        PacketSet expected_packets{};
        if (!BuildReference(original, expected_descriptor, expected_packets)) {
            return Check(false, "reference construction/capture");
        }

        for (unsigned route_value = 0u;
             route_value < static_cast<unsigned>(ConstructorRoute::Count);
             ++route_value)
        {
            const ConstructorRoute route =
                static_cast<ConstructorRoute>(route_value);
            for (unsigned flavor_value = 0u;
                 flavor_value < static_cast<unsigned>(SourceFlavor::Count);
                 ++flavor_value)
            {
                const SourceFlavor flavor =
                    static_cast<SourceFlavor>(flavor_value);
                std::vector<uint8_t> source = original;
                uint8_t guarded_descriptor[ProfileBytes + 2u];
                std::memset(
                    guarded_descriptor, Guard, sizeof(guarded_descriptor));
                uint8_t* const descriptor = guarded_descriptor + 1u;
                if (route == ConstructorRoute::Descriptor) {
                    std::memcpy(
                        descriptor, expected_descriptor.data(), ProfileBytes);
                }
                uint32_t descriptor_bytes =
                    route == ConstructorRoute::Descriptor ? ProfileBytes : 0u;
                CodecOwner owner;
                const WirehairV2Result result = CreateEncoder(
                    route, flavor,
                    source.data(), source.size(), BlockBytes,
                    descriptor, &descriptor_bytes, &owner.Handle);
                if (result != WirehairV2_Success || !owner.Handle ||
                    descriptor_bytes != ProfileBytes ||
                    guarded_descriptor[0] != Guard ||
                    guarded_descriptor[sizeof(guarded_descriptor) - 1u] !=
                        Guard ||
                    std::memcmp(
                        descriptor, expected_descriptor.data(),
                        ProfileBytes) != 0)
                {
                    std::fprintf(stderr,
                        "V2 borrowed-source matrix create failed: "
                        "message_bytes=%llu route=%s flavor=%s result=%d\n",
                        static_cast<unsigned long long>(message_bytes),
                        RouteName(route), FlavorName(flavor),
                        static_cast<int>(result));
                    return false;
                }

                PacketSet attached_packets{};
                if (!CapturePackets(
                        owner.Handle, message_bytes, original,
                        attached_packets, "attached") ||
                    !PacketSetsEqual(attached_packets, expected_packets) ||
                    !ReplayPackets(
                        owner.Handle, expected_packets,
                        "attached exact/short capacity") ||
                    source != original)
                {
                    std::fprintf(stderr,
                        "V2 borrowed-source attached packet mismatch: "
                        "message_bytes=%llu route=%s flavor=%s\n",
                        static_cast<unsigned long long>(message_bytes),
                        RouteName(route), FlavorName(flavor));
                    return false;
                }

                if (flavor != SourceFlavor::Borrowed)
                {
                    PoisonAndRelease(source, original);
                    if (!ReplayPackets(
                            owner.Handle, expected_packets,
                            "source-independent before detach"))
                    {
                        return false;
                    }
                }

                if (wirehair_v2_encoder_detach_input(owner.Handle) !=
                        WirehairV2_Success ||
                    wirehair_v2_encoder_detach_input(owner.Handle) !=
                        WirehairV2_Success)
                {
                    std::fprintf(stderr,
                        "V2 borrowed-source detach failed: "
                        "message_bytes=%llu route=%s flavor=%s\n",
                        static_cast<unsigned long long>(message_bytes),
                        RouteName(route), FlavorName(flavor));
                    return false;
                }
                if (flavor == SourceFlavor::Borrowed) {
                    PoisonAndRelease(source, original);
                }
                if (!ReplayPackets(
                        owner.Handle, expected_packets,
                        "after detach and source release") ||
                    guarded_descriptor[0] != Guard ||
                    guarded_descriptor[sizeof(guarded_descriptor) - 1u] !=
                        Guard ||
                    std::memcmp(
                        descriptor, expected_descriptor.data(),
                        ProfileBytes) != 0)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

WirehairV2Result CallSelectingWithOptions(
    bool explicit_profile,
    uint64_t profile_id,
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    const WirehairV2EncoderOptions* options,
    void* descriptor,
    uint32_t descriptor_capacity,
    uint32_t* descriptor_bytes,
    WirehairV2Codec* codec_out)
{
    if (explicit_profile) {
        return wirehair_v2_encoder_create_profile_id_with_options(
            profile_id, message, message_bytes, block_bytes, options,
            descriptor, descriptor_capacity, descriptor_bytes, codec_out);
    }
    return wirehair_v2_encoder_create_with_options(
        message, message_bytes, block_bytes, options,
        descriptor, descriptor_capacity, descriptor_bytes, codec_out);
}

bool CheckSelectingOptionPriority()
{
    struct OptionCase
    {
        const char* Name;
        WirehairV2EncoderOptions Options;
        WirehairV2Result Expected;
    };

    const WirehairV2EncoderOptions valid = MakeOptions(
        static_cast<uint32_t>(WirehairV2EncoderSource_Independent));
    OptionCase cases[8];
    for (OptionCase& item : cases) {
        item.Options = valid;
    }
    cases[0].Name = "size";
    cases[0].Options.struct_bytes = sizeof(valid) - 1u;
    cases[0].Expected = WirehairV2_InvalidSize;
    cases[1].Name = "version";
    cases[1].Options.options_version =
        WIREHAIR_V2_ENCODER_OPTIONS_VERSION + 1u;
    cases[1].Expected = WirehairV2_UnsupportedVersion;
    cases[2].Name = "reserved";
    cases[2].Options.reserved = 1u;
    cases[2].Expected = WirehairV2_ReservedNonzero;
    cases[3].Name = "zero-policy";
    cases[3].Options.source_policy = 0u;
    cases[3].Expected = WirehairV2_InvalidInput;
    cases[4].Name = "count-policy";
    cases[4].Options.source_policy =
        static_cast<uint32_t>(WirehairV2EncoderSource_Count);
    cases[4].Expected = WirehairV2_InvalidInput;
    cases[5].Name = "unknown-policy";
    cases[5].Options.source_policy = UINT32_MAX;
    cases[5].Expected = WirehairV2_InvalidInput;
    cases[6].Name = "simultaneous";
    cases[6].Options.struct_bytes = 0u;
    cases[6].Options.options_version = UINT32_MAX;
    cases[6].Options.reserved = UINT32_MAX;
    cases[6].Options.source_policy = UINT32_MAX;
    cases[6].Expected = WirehairV2_InvalidSize;
    cases[7].Name = "all-zero";
    std::memset(&cases[7].Options, 0, sizeof(cases[7].Options));
    cases[7].Expected = WirehairV2_InvalidSize;

    std::vector<uint8_t> message(BlockCount * BlockBytes);
    FillMessage(message);
    for (unsigned explicit_profile = 0u; explicit_profile < 2u;
         ++explicit_profile)
    {
        for (const OptionCase& item : cases)
        {
            uint8_t descriptor[ProfileBytes];
            std::memset(descriptor, 0xa5, sizeof(descriptor));
            const std::array<uint8_t, ProfileBytes> descriptor_before = [&]() {
                std::array<uint8_t, ProfileBytes> copy{};
                std::copy(
                    descriptor, descriptor + ProfileBytes, copy.begin());
                return copy;
            }();
            uint32_t required = UINT32_C(0xa5a5a5a5);
            WirehairV2Codec codec =
                reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
            const WirehairV2EncoderOptions options_before = item.Options;
            const WirehairV2Result result = CallSelectingWithOptions(
                explicit_profile != 0u,
                explicit_profile ? UINT64_MAX : WIREHAIR_V2_PROFILE_CURRENT,
                message.data(), message.size(), 0u,
                &item.Options,
                descriptor, ProfileBytes - 1u, &required, &codec);
            if (result != item.Expected ||
                required != ProfileBytes || codec != nullptr ||
                std::memcmp(
                    descriptor, descriptor_before.data(),
                    sizeof(descriptor)) != 0 ||
                std::memcmp(
                    &item.Options, &options_before,
                    sizeof(item.Options)) != 0)
            {
                std::fprintf(stderr,
                    "V2 borrowed-source selecting option priority failed: "
                    "selector=%u case=%s result=%d expected=%d\n",
                    explicit_profile, item.Name,
                    static_cast<int>(result),
                    static_cast<int>(item.Expected));
                return false;
            }
        }

        uint8_t descriptor[ProfileBytes];
        std::memset(descriptor, 0xa5, sizeof(descriptor));
        uint32_t required = UINT32_C(0xa5a5a5a5);
        WirehairV2Codec codec =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (CallSelectingWithOptions(
                explicit_profile != 0u,
                WIREHAIR_V2_PROFILE_CURRENT,
                message.data(), message.size(), BlockBytes,
                nullptr, descriptor, sizeof(descriptor),
                &required, &codec) != WirehairV2_InvalidInput ||
            required != ProfileBytes || codec != nullptr ||
            !std::all_of(
                descriptor, descriptor + sizeof(descriptor),
                [](uint8_t byte) { return byte == 0xa5u; }))
        {
            return Check(false, "null selecting options priority");
        }

        for (uint32_t policy : {
                 static_cast<uint32_t>(
                     WirehairV2EncoderSource_Independent),
                 static_cast<uint32_t>(
                     WirehairV2EncoderSource_BorrowedImmutable) })
        {
            WirehairV2EncoderOptions options = MakeOptions(policy);
            std::memset(descriptor, 0xa5, sizeof(descriptor));
            required = UINT32_C(0xa5a5a5a5);
            codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
            const WirehairV2Result result = CallSelectingWithOptions(
                explicit_profile != 0u,
                explicit_profile ? UINT64_MAX : WIREHAIR_V2_PROFILE_CURRENT,
                message.data(), message.size(), BlockBytes,
                &options, descriptor, ProfileBytes - 1u,
                &required, &codec);
            if (result != WirehairV2_BufferTooSmall ||
                required != ProfileBytes ||
                codec != nullptr ||
                !std::all_of(
                    descriptor, descriptor + sizeof(descriptor),
                    [](uint8_t byte) { return byte == 0xa5u; }))
            {
                return Check(
                    false,
                    "valid options preserve ordinary selecting priority");
            }
        }

        if (explicit_profile)
        {
            WirehairV2EncoderOptions options = valid;
            std::memset(descriptor, 0xa5, sizeof(descriptor));
            required = UINT32_C(0xa5a5a5a5);
            codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
            if (CallSelectingWithOptions(
                    true, UINT64_MAX,
                    message.data(), message.size(), BlockBytes,
                    &options, descriptor, sizeof(descriptor),
                    &required, &codec) != WirehairV2_UnsupportedProfile ||
                required != ProfileBytes || codec != nullptr ||
                !std::all_of(
                    descriptor, descriptor + sizeof(descriptor),
                    [](uint8_t byte) { return byte == 0xa5u; }))
            {
                return Check(
                    false, "capacity must precede explicit profile selection");
            }
        }

        WirehairV2EncoderOptions simultaneous = valid;
        simultaneous.struct_bytes = 0u;
        simultaneous.options_version = UINT32_MAX;
        simultaneous.source_policy = UINT32_MAX;
        simultaneous.reserved = UINT32_MAX;
        std::memset(descriptor, 0xa5, sizeof(descriptor));
        required = UINT32_C(0xa5a5a5a5);
        if (CallSelectingWithOptions(
                explicit_profile != 0u,
                explicit_profile ? UINT64_MAX :
                    WIREHAIR_V2_PROFILE_CURRENT,
                message.data(), message.size(), 0u,
                &simultaneous, descriptor, ProfileBytes - 1u,
                &required, nullptr) != WirehairV2_InvalidInput ||
            required != ProfileBytes ||
            !std::all_of(
                descriptor, descriptor + sizeof(descriptor),
                [](uint8_t byte) { return byte == 0xa5u; }))
        {
            return Check(false, "null selecting codec output priority");
        }

        codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (CallSelectingWithOptions(
                explicit_profile != 0u,
                explicit_profile ? UINT64_MAX :
                    WIREHAIR_V2_PROFILE_CURRENT,
                message.data(), message.size(), 0u,
                &simultaneous, descriptor, ProfileBytes - 1u,
                nullptr, &codec) != WirehairV2_InvalidInput ||
            codec != nullptr ||
            !std::all_of(
                descriptor, descriptor + sizeof(descriptor),
                [](uint8_t byte) { return byte == 0xa5u; }))
        {
            return Check(false, "null selecting size output priority");
        }
    }
    return true;
}

bool CheckDescriptorOptionPriority(
    const std::vector<uint8_t>& message,
    const std::array<uint8_t, ProfileBytes>& descriptor)
{
    struct OptionCase
    {
        const char* Name;
        WirehairV2EncoderOptions Options;
        WirehairV2Result Expected;
    };
    const WirehairV2EncoderOptions valid = MakeOptions(
        static_cast<uint32_t>(WirehairV2EncoderSource_BorrowedImmutable));
    OptionCase cases[8];
    for (OptionCase& item : cases) {
        item.Options = valid;
    }
    cases[0].Name = "size";
    cases[0].Options.struct_bytes = sizeof(valid) - 1u;
    cases[0].Expected = WirehairV2_InvalidSize;
    cases[1].Name = "version";
    cases[1].Options.options_version =
        WIREHAIR_V2_ENCODER_OPTIONS_VERSION + 1u;
    cases[1].Expected = WirehairV2_UnsupportedVersion;
    cases[2].Name = "reserved";
    cases[2].Options.reserved = 1u;
    cases[2].Expected = WirehairV2_ReservedNonzero;
    cases[3].Name = "zero-policy";
    cases[3].Options.source_policy = 0u;
    cases[3].Expected = WirehairV2_InvalidInput;
    cases[4].Name = "count-policy";
    cases[4].Options.source_policy =
        static_cast<uint32_t>(WirehairV2EncoderSource_Count);
    cases[4].Expected = WirehairV2_InvalidInput;
    cases[5].Name = "unknown-policy";
    cases[5].Options.source_policy = UINT32_MAX;
    cases[5].Expected = WirehairV2_InvalidInput;
    cases[6].Name = "simultaneous";
    cases[6].Options.struct_bytes = 0u;
    cases[6].Options.options_version = UINT32_MAX;
    cases[6].Options.reserved = UINT32_MAX;
    cases[6].Options.source_policy = UINT32_MAX;
    cases[6].Expected = WirehairV2_InvalidSize;
    cases[7].Name = "all-zero";
    std::memset(&cases[7].Options, 0, sizeof(cases[7].Options));
    cases[7].Expected = WirehairV2_InvalidSize;

    std::array<uint8_t, ProfileBytes> malformed = descriptor;
    malformed[0] ^= 0xffu;
    for (const OptionCase& item : cases)
    {
        WirehairV2Codec codec =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (wirehair_v2_encoder_create_profile_with_options(
                message.data(), malformed.data(), malformed.size(),
                &item.Options, &codec) != WirehairV2_InvalidMagic ||
            codec != nullptr)
        {
            std::fprintf(stderr,
                "V2 descriptor/options parse priority failed: case=%s\n",
                item.Name);
            return false;
        }

        codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (wirehair_v2_encoder_create_profile_with_options(
                message.data(), descriptor.data(), descriptor.size(),
                &item.Options, &codec) != item.Expected ||
            codec != nullptr)
        {
            std::fprintf(stderr,
                "V2 descriptor/options validation priority failed: "
                "case=%s\n", item.Name);
            return false;
        }

        codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (wirehair_v2_encoder_create_profile_with_options(
                message.data(), descriptor.data(), ProfileBytes - 1u,
                &item.Options, &codec) != WirehairV2_InvalidSize ||
            codec != nullptr)
        {
            std::fprintf(stderr,
                "V2 short descriptor/options priority failed: case=%s\n",
                item.Name);
            return false;
        }
    }

    WirehairV2Codec codec =
        reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    if (wirehair_v2_encoder_create_profile_with_options(
            message.data(), descriptor.data(), descriptor.size(),
            nullptr, &codec) != WirehairV2_InvalidInput || codec != nullptr)
    {
        return Check(false, "null descriptor options priority");
    }

    codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    if (wirehair_v2_encoder_create_profile_with_options(
            message.data(), malformed.data(), malformed.size(),
            nullptr, &codec) != WirehairV2_InvalidMagic || codec != nullptr)
    {
        return Check(false, "malformed descriptor outranks null options");
    }

    // A null codec output is the first descriptor-route check and must not
    // dereference any other argument.
    const void* const unreadable = reinterpret_cast<const void*>(uintptr_t(1));
    if (wirehair_v2_encoder_create_profile_with_options(
            unreadable, unreadable, ProfileBytes,
            reinterpret_cast<const WirehairV2EncoderOptions*>(unreadable),
            nullptr) != WirehairV2_InvalidInput)
    {
        return Check(false, "null descriptor codec output priority");
    }
    return true;
}

bool CheckSelectingOptionsOverlap()
{
    enum class Overlap : unsigned
    {
        Descriptor,
        Size,
        Codec,
        Count
    };
    std::vector<uint8_t> message(BlockCount * BlockBytes);
    FillMessage(message);
    const std::vector<uint8_t> message_before = message;

    for (unsigned explicit_profile = 0u; explicit_profile < 2u;
         ++explicit_profile)
    {
        for (unsigned overlap_value = 0u;
             overlap_value < static_cast<unsigned>(Overlap::Count);
             ++overlap_value)
        {
            for (unsigned partial = 0u; partial < 2u; ++partial)
            {
                alignas(std::max_align_t) uint8_t storage[512];
                std::memset(storage, 0xa5, sizeof(storage));
                WirehairV2EncoderOptions* const options = new (storage + 64u)
                    WirehairV2EncoderOptions(MakeOptions(UINT32_MAX));
                options->struct_bytes = 0u;
                options->options_version = UINT32_MAX;
                options->reserved = UINT32_MAX;

                void* descriptor = storage + 256u;
                uint32_t* size_out =
                    reinterpret_cast<uint32_t*>(storage + 320u);
                WirehairV2Codec* codec_out =
                    reinterpret_cast<WirehairV2Codec*>(storage + 336u);
                switch (static_cast<Overlap>(overlap_value))
                {
                case Overlap::Descriptor:
                    descriptor = storage + (partial ? 66u : 64u);
                    break;
                case Overlap::Size:
                    size_out = reinterpret_cast<uint32_t*>(
                        storage + (partial ? 76u : 64u));
                    break;
                case Overlap::Codec:
                    codec_out = reinterpret_cast<WirehairV2Codec*>(
                        storage + (partial ? 72u : 64u));
                    break;
                default:
                    return false;
                }
                uint8_t before[sizeof(storage)];
                std::memcpy(before, storage, sizeof(before));
                const WirehairV2Result result = CallSelectingWithOptions(
                    explicit_profile != 0u,
                    explicit_profile ? UINT64_MAX :
                        WIREHAIR_V2_PROFILE_CURRENT,
                    message.data(), message.size(), 0u,
                    options, descriptor, ProfileBytes - 1u,
                    size_out, codec_out);
                if (result != WirehairV2_InvalidInput ||
                    std::memcmp(storage, before, sizeof(storage)) != 0 ||
                    message != message_before)
                {
                    std::fprintf(stderr,
                        "V2 borrowed-source options overlap failed: "
                        "selector=%u overlap=%u partial=%u result=%d\n",
                        explicit_profile, overlap_value, partial,
                        static_cast<int>(result));
                    return false;
                }
            }
        }
    }
    return true;
}

bool CheckSelectingFixedOverlaps()
{
    enum class Overlap : unsigned
    {
        DescriptorSize,
        DescriptorCodec,
        SizeCodec,
        SizeMessage,
        CodecMessage,
        Count
    };
    std::vector<uint8_t> external_message(BlockCount * BlockBytes);
    FillMessage(external_message);
    const std::vector<uint8_t> external_before = external_message;

    for (unsigned explicit_profile = 0u; explicit_profile < 2u;
         ++explicit_profile)
    {
        for (uint32_t policy : {
                 static_cast<uint32_t>(
                     WirehairV2EncoderSource_Independent),
                 static_cast<uint32_t>(
                     WirehairV2EncoderSource_BorrowedImmutable),
                 UINT32_MAX })
        {
            WirehairV2EncoderOptions options = MakeOptions(policy);
            if (policy == UINT32_MAX) {
                options.struct_bytes = 0u;
                options.options_version = UINT32_MAX;
                options.reserved = UINT32_MAX;
            }
            for (unsigned overlap_value = 0u;
                 overlap_value < static_cast<unsigned>(Overlap::Count);
                 ++overlap_value)
            {
                for (unsigned partial = 0u; partial < 2u; ++partial)
                {
                    if (static_cast<Overlap>(overlap_value) ==
                            Overlap::SizeCodec &&
                        partial != 0u &&
                        sizeof(WirehairV2Codec) == sizeof(uint32_t))
                    {
                        continue;
                    }

                    alignas(std::max_align_t) uint8_t storage[512];
                    std::memset(storage, 0xa5, sizeof(storage));
                    const void* message = external_message.data();
                    void* descriptor = storage + 256u;
                    uint32_t* size_out =
                        reinterpret_cast<uint32_t*>(storage + 320u);
                    WirehairV2Codec* codec_out =
                        reinterpret_cast<WirehairV2Codec*>(storage + 336u);
                    switch (static_cast<Overlap>(overlap_value))
                    {
                    case Overlap::DescriptorSize:
                        descriptor = storage + (partial ? 66u : 64u);
                        size_out = reinterpret_cast<uint32_t*>(storage + 64u);
                        break;
                    case Overlap::DescriptorCodec:
                        descriptor = storage + (partial ? 66u : 64u);
                        codec_out = reinterpret_cast<WirehairV2Codec*>(
                            storage + 64u);
                        break;
                    case Overlap::SizeCodec:
                        size_out = reinterpret_cast<uint32_t*>(
                            storage + (partial ? 68u : 64u));
                        codec_out = reinterpret_cast<WirehairV2Codec*>(
                            storage + 64u);
                        break;
                    case Overlap::SizeMessage:
                        message = storage + (partial ? 66u : 64u);
                        size_out = reinterpret_cast<uint32_t*>(storage + 64u);
                        break;
                    case Overlap::CodecMessage:
                        message = storage + (partial ?
                            64u + sizeof(WirehairV2Codec) / 2u : 64u);
                        codec_out = reinterpret_cast<WirehairV2Codec*>(
                            storage + 64u);
                        break;
                    default:
                        return false;
                    }

                    uint8_t before[sizeof(storage)];
                    std::memcpy(before, storage, sizeof(before));
                    const WirehairV2Result result = CallSelectingWithOptions(
                        explicit_profile != 0u,
                        explicit_profile ? UINT64_MAX :
                            WIREHAIR_V2_PROFILE_CURRENT,
                        message, external_message.size(), 0u,
                        &options, descriptor, ProfileBytes - 1u,
                        size_out, codec_out);
                    if (result != WirehairV2_InvalidInput ||
                        std::memcmp(storage, before, sizeof(storage)) != 0 ||
                        external_message != external_before)
                    {
                        std::fprintf(stderr,
                            "V2 borrowed-source fixed overlap failed: "
                            "selector=%u policy=%u overlap=%u partial=%u "
                            "result=%d\n",
                            explicit_profile, policy, overlap_value, partial,
                            static_cast<int>(result));
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool CheckBorrowedDescriptorMessageOverlap(
    const std::array<uint8_t, ProfileBytes>& expected_descriptor,
    const PacketSet& expected_packets)
{
    for (unsigned explicit_profile = 0u; explicit_profile < 2u;
         ++explicit_profile)
    {
        for (unsigned partial = 0u; partial < 2u; ++partial)
        {
            for (SourceFlavor flavor : {
                     SourceFlavor::Legacy,
                     SourceFlavor::Independent,
                     SourceFlavor::Borrowed })
            {
                std::vector<uint8_t> source(BlockCount * BlockBytes);
                FillMessage(source);
                const std::vector<uint8_t> before = source;
                uint8_t* const descriptor =
                    source.data() + (partial ? 1u : 0u);
                uint32_t descriptor_bytes = UINT32_C(0xa5a5a5a5);
                WirehairV2Codec codec =
                    reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
                const WirehairV2Result result = CreateEncoder(
                    explicit_profile ? ConstructorRoute::ExplicitSelector :
                        ConstructorRoute::DefaultSelector,
                    flavor,
                    source.data(), source.size(), BlockBytes,
                    descriptor, &descriptor_bytes, &codec);
                if (flavor == SourceFlavor::Borrowed)
                {
                    if (result != WirehairV2_InvalidInput ||
                        source != before ||
                        descriptor_bytes != UINT32_C(0xa5a5a5a5) ||
                        codec !=
                            reinterpret_cast<WirehairV2Codec>(uintptr_t(1)))
                    {
                        return Check(
                            false,
                            "borrowed descriptor/message overlap no-write");
                    }
                }
                else
                {
                    if (result != WirehairV2_Success || !codec ||
                        codec == reinterpret_cast<WirehairV2Codec>(
                            uintptr_t(1)) ||
                        descriptor_bytes != ProfileBytes ||
                        std::memcmp(
                            descriptor, expected_descriptor.data(),
                            ProfileBytes) != 0)
                    {
                        return Check(
                            false,
                            "legacy/independent descriptor/message alias");
                    }
                    CodecOwner owner;
                    owner.Handle = codec;
                    if (!ReplayPackets(
                            owner.Handle, expected_packets,
                            "allowed descriptor/message alias"))
                    {
                        return false;
                    }
                }
            }

            std::vector<uint8_t> source(BlockCount * BlockBytes);
            FillMessage(source);
            const std::vector<uint8_t> before = source;
            uint8_t* const descriptor =
                source.data() + (partial ? 1u : 0u);
            WirehairV2EncoderOptions options = MakeOptions(
                static_cast<uint32_t>(
                    WirehairV2EncoderSource_BorrowedImmutable));
            const WirehairV2EncoderOptions options_before = options;
            uint32_t descriptor_bytes = UINT32_C(0xa5a5a5a5);
            WirehairV2Codec codec =
                reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
            const WirehairV2Result result = CallSelectingWithOptions(
                explicit_profile != 0u,
                explicit_profile ? UINT64_MAX :
                    WIREHAIR_V2_PROFILE_CURRENT,
                source.data(), source.size(), 0u, &options,
                descriptor, ProfileBytes - 1u,
                &descriptor_bytes, &codec);
            if (result != WirehairV2_InvalidInput || source != before ||
                descriptor_bytes != UINT32_C(0xa5a5a5a5) ||
                codec != reinterpret_cast<WirehairV2Codec>(uintptr_t(1)) ||
                std::memcmp(
                    &options, &options_before, sizeof(options)) != 0)
            {
                return Check(
                    false,
                    "borrowed descriptor/message overlap result priority");
            }
        }
    }
    return true;
}

bool CheckEmptyBorrowedMessageDoesNotOverlapDescriptor()
{
    for (unsigned explicit_profile = 0u; explicit_profile < 2u;
         ++explicit_profile)
    {
        std::array<uint8_t, ProfileBytes> descriptor{};
        descriptor.fill(0xa5u);
        const std::array<uint8_t, ProfileBytes> descriptor_before =
            descriptor;
        WirehairV2EncoderOptions options = MakeOptions(
            static_cast<uint32_t>(
                WirehairV2EncoderSource_BorrowedImmutable));
        const WirehairV2EncoderOptions options_before = options;
        uint32_t descriptor_bytes = UINT32_C(0xa5a5a5a5);
        WirehairV2Codec codec =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));

        const WirehairV2Result result = CallSelectingWithOptions(
            explicit_profile != 0u,
            WIREHAIR_V2_PROFILE_CURRENT,
            descriptor.data() + 1u, 0u, BlockBytes, &options,
            descriptor.data(), descriptor.size(),
            &descriptor_bytes, &codec);
        if (result != WirehairV2_InvalidDimensions ||
            descriptor_bytes != ProfileBytes || codec != nullptr ||
            descriptor != descriptor_before ||
            std::memcmp(
                &options, &options_before, sizeof(options)) != 0)
        {
            std::fprintf(stderr,
                "V2 empty borrowed message overlap failed: "
                "selector=%u result=%d\n",
                explicit_profile, static_cast<int>(result));
            return false;
        }

        if (explicit_profile)
        {
            descriptor.fill(0xa5u);
            descriptor_bytes = UINT32_C(0xa5a5a5a5);
            codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
            const WirehairV2Result unsupported = CallSelectingWithOptions(
                true, UINT64_MAX,
                descriptor.data() + 1u, 0u, BlockBytes, &options,
                descriptor.data(), descriptor.size(),
                &descriptor_bytes, &codec);
            if (unsupported != WirehairV2_UnsupportedProfile ||
                descriptor_bytes != ProfileBytes || codec != nullptr ||
                descriptor != descriptor_before ||
                std::memcmp(
                    &options, &options_before, sizeof(options)) != 0)
            {
                std::fprintf(stderr,
                    "V2 empty borrowed message/profile priority failed: "
                    "result=%d\n", static_cast<int>(unsupported));
                return false;
            }
        }
    }
    return true;
}

bool CheckReadOnlyInputOverlaps(
    const std::array<uint8_t, ProfileBytes>& descriptor)
{
    struct OptionMessageStorage
    {
        WirehairV2EncoderOptions Options;
        uint8_t Tail[BlockCount * BlockBytes];
    };

    for (unsigned explicit_profile = 0u; explicit_profile < 2u;
         ++explicit_profile)
    {
        for (unsigned partial = 0u; partial < 2u; ++partial)
        {
            OptionMessageStorage storage = {
                MakeOptions(static_cast<uint32_t>(
                    WirehairV2EncoderSource_BorrowedImmutable)),
                {}
            };
            for (size_t i = 0u; i < sizeof(storage.Tail); ++i) {
                storage.Tail[i] = static_cast<uint8_t>(i * 73u + 19u);
            }
            const uint8_t* const message =
                reinterpret_cast<const uint8_t*>(&storage) +
                (partial ? sizeof(storage.Options) / 2u : 0u);
            std::vector<uint8_t> expected(
                message, message + BlockCount * BlockBytes);
            std::array<uint8_t, ProfileBytes> serialized{};
            uint32_t serialized_bytes = 0u;
            CodecOwner owner;
            const WirehairV2Result result = CallSelectingWithOptions(
                explicit_profile != 0u, WIREHAIR_V2_PROFILE_CURRENT,
                message, expected.size(), BlockBytes, &storage.Options,
                serialized.data(), serialized.size(),
                &serialized_bytes, &owner.Handle);
            PacketSet packets{};
            if (result != WirehairV2_Success || !owner.Handle ||
                serialized_bytes != ProfileBytes ||
                !CapturePackets(
                    owner.Handle, expected.size(), expected,
                    packets, "message/options read-only overlap"))
            {
                return Check(false, "message/options read-only overlap");
            }
        }
    }

    for (unsigned partial = 0u; partial < 2u; ++partial)
    {
        std::vector<uint8_t> source(BlockCount * BlockBytes + 1u);
        FillMessage(source);
        uint8_t* const serialized = source.data() + (partial ? 1u : 0u);
        std::copy(descriptor.begin(), descriptor.end(), serialized);
        std::vector<uint8_t> expected(
            source.begin(), source.begin() + BlockCount * BlockBytes);
        WirehairV2EncoderOptions options = MakeOptions(
            static_cast<uint32_t>(
                WirehairV2EncoderSource_BorrowedImmutable));
        CodecOwner owner;
        if (wirehair_v2_encoder_create_profile_with_options(
                source.data(), serialized, ProfileBytes,
                &options, &owner.Handle) != WirehairV2_Success ||
            !owner.Handle)
        {
            return Check(false, "descriptor/message read-only overlap");
        }
        PacketSet packets{};
        if (!CapturePackets(
                owner.Handle, expected.size(), expected,
                packets, "descriptor/message read-only overlap"))
        {
            return false;
        }
    }

    const uint16_t endian_probe = 1u;
    const bool little_endian =
        *reinterpret_cast<const uint8_t*>(&endian_probe) == 1u;
    for (uint32_t block_count = 2u; block_count <= 64u; ++block_count)
    {
        alignas(std::max_align_t) uint8_t storage[64] = {};
        uint8_t* const serialized = storage + (little_endian ? 0u : 1u);
        WirehairV2Profile profile = {};
        profile.struct_bytes = sizeof(profile);
        profile.profile_version = WIREHAIR_V2_PROFILE_VERSION;
        profile.profile_id = WIREHAIR_V2_PROFILE_CURRENT;
        profile.message_bytes =
            static_cast<uint64_t>(block_count) * BlockBytes;
        profile.block_bytes = BlockBytes;
        profile.seed_attempt = little_endian ? 1u : 0u;
        uint32_t serialized_bytes = 0u;
        if (wirehair_v2_profile_serialize(
                &profile, serialized, ProfileBytes,
                &serialized_bytes) != WirehairV2_Success ||
            serialized_bytes != ProfileBytes)
        {
            return Check(false, "descriptor/options overlap fixture");
        }

        WirehairV2EncoderOptions* const options = little_endian ?
            new (serialized + 24u) WirehairV2EncoderOptions(MakeOptions(
                static_cast<uint32_t>(
                    WirehairV2EncoderSource_BorrowedImmutable))) :
            new (serialized + 31u) WirehairV2EncoderOptions(MakeOptions(
                static_cast<uint32_t>(
                    WirehairV2EncoderSource_BorrowedImmutable)));
        std::vector<uint8_t> message(
            static_cast<size_t>(profile.message_bytes));
        FillMessage(message);
        WirehairV2Codec codec = nullptr;
        const WirehairV2Result result =
            wirehair_v2_encoder_create_profile_with_options(
                message.data(), serialized, ProfileBytes,
                options, &codec);
        if (result == WirehairV2_BadSeed && codec == nullptr) {
            continue;
        }
        if (result != WirehairV2_Success || !codec) {
            return Check(false, "descriptor/options read-only overlap");
        }
        CodecOwner owner;
        owner.Handle = codec;
        uint8_t output[BlockBytes];
        uint32_t output_bytes = 0u;
        if (wirehair_v2_encode(
                owner.Handle, 0u, output, sizeof(output),
                &output_bytes) != WirehairV2_Success ||
            output_bytes != BlockBytes ||
            std::memcmp(output, message.data(), BlockBytes) != 0)
        {
            return Check(false, "descriptor/options overlap packet");
        }
        return true;
    }
    return Check(false, "descriptor/options overlap found no usable seed");
}

bool CheckDescriptorFixedOverlaps(
    const std::array<uint8_t, ProfileBytes>& descriptor)
{
    enum class Overlap : unsigned
    {
        Descriptor,
        Options,
        Message,
        Count
    };
    for (unsigned overlap_value = 0u;
         overlap_value < static_cast<unsigned>(Overlap::Count);
         ++overlap_value)
    {
        for (unsigned partial = 0u; partial < 2u; ++partial)
        {
            for (unsigned variant = 0u; variant < 3u; ++variant)
            {
                alignas(std::max_align_t) uint8_t storage[512];
                std::memset(storage, 0xa5, sizeof(storage));
                WirehairV2EncoderOptions* const options = new (storage + 64u)
                    WirehairV2EncoderOptions(MakeOptions(
                        static_cast<uint32_t>(
                            variant == 2u ?
                                WirehairV2EncoderSource_Independent :
                                WirehairV2EncoderSource_BorrowedImmutable)));
                uint8_t* const serialized = storage + 128u;
                std::copy(descriptor.begin(), descriptor.end(), serialized);
                uint8_t* message = storage + 256u;
                std::vector<uint8_t> expected_message(
                    BlockCount * BlockBytes);
                FillMessage(expected_message);
                std::copy(
                    expected_message.begin(), expected_message.end(), message);

                if (variant == 1u)
                {
                    options->struct_bytes = 0u;
                    options->options_version = UINT32_MAX;
                    options->source_policy = UINT32_MAX;
                    options->reserved = UINT32_MAX;
                    if (static_cast<Overlap>(overlap_value) !=
                        Overlap::Message)
                    {
                        serialized[0] ^= 0xffu;
                    }
                }

                WirehairV2Codec* codec_out = nullptr;
                switch (static_cast<Overlap>(overlap_value))
                {
                case Overlap::Descriptor:
                    codec_out = reinterpret_cast<WirehairV2Codec*>(
                        storage + (partial ? 152u : 128u));
                    break;
                case Overlap::Options:
                    codec_out = reinterpret_cast<WirehairV2Codec*>(
                        storage + (partial ? 72u : 64u));
                    break;
                case Overlap::Message:
                    codec_out = reinterpret_cast<WirehairV2Codec*>(
                        storage + (partial ? 376u : 256u));
                    break;
                default:
                    return false;
                }

                uint8_t before[sizeof(storage)];
                std::memcpy(before, storage, sizeof(before));
                const WirehairV2Result result =
                    wirehair_v2_encoder_create_profile_with_options(
                        message, serialized, ProfileBytes,
                        options, codec_out);
                if (result != WirehairV2_InvalidInput ||
                    std::memcmp(storage, before, sizeof(storage)) != 0)
                {
                    std::fprintf(stderr,
                        "V2 borrowed-source descriptor overlap failed: "
                        "overlap=%u partial=%u variant=%u result=%d\n",
                        overlap_value, partial, variant,
                        static_cast<int>(result));
                    return false;
                }
            }
        }
    }
    return true;
}

bool CheckWrappedBorrowedMessageRange(
    const std::array<uint8_t, ProfileBytes>& descriptor)
{
    const uintptr_t limit = std::numeric_limits<uintptr_t>::max();
    const void* const wrapped = reinterpret_cast<const void*>(limit - 63u);
    for (unsigned explicit_profile = 0u; explicit_profile < 2u;
         ++explicit_profile)
    {
        WirehairV2EncoderOptions options = MakeOptions(
            static_cast<uint32_t>(
                WirehairV2EncoderSource_BorrowedImmutable));
        uint8_t descriptor[ProfileBytes];
        std::memset(descriptor, 0xa5, sizeof(descriptor));
        uint32_t descriptor_bytes = UINT32_C(0xa5a5a5a5);
        WirehairV2Codec codec =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        const WirehairV2Result result = CallSelectingWithOptions(
            explicit_profile != 0u,
            WIREHAIR_V2_PROFILE_CURRENT,
            wrapped, static_cast<uint64_t>(BlockCount) * BlockBytes,
            BlockBytes, &options,
            descriptor, sizeof(descriptor),
            &descriptor_bytes, &codec);
        if (result != WirehairV2_InvalidInput ||
            descriptor_bytes != UINT32_C(0xa5a5a5a5) ||
            codec != reinterpret_cast<WirehairV2Codec>(uintptr_t(1)) ||
            !std::all_of(
                descriptor, descriptor + sizeof(descriptor),
                [](uint8_t byte) { return byte == 0xa5u; }))
        {
            return Check(
                false, "wrapping borrowed message range is no-write");
        }
    }

    WirehairV2EncoderOptions options = MakeOptions(
        static_cast<uint32_t>(
            WirehairV2EncoderSource_BorrowedImmutable));
    const WirehairV2EncoderOptions options_before = options;
    WirehairV2Codec codec =
        reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    if (wirehair_v2_encoder_create_profile_with_options(
            wrapped, descriptor.data(), descriptor.size(),
            &options, &codec) != WirehairV2_InvalidInput ||
        codec != reinterpret_cast<WirehairV2Codec>(uintptr_t(1)) ||
        std::memcmp(&options, &options_before, sizeof(options)) != 0)
    {
        return Check(
            false, "descriptor-route wrapping message range is no-write");
    }

#if SIZE_MAX < UINT64_MAX
    uint8_t serialized[ProfileBytes];
    std::copy(descriptor.begin(), descriptor.end(), serialized);
    WirehairV2Profile oversized = {};
    if (wirehair_v2_profile_deserialize(
            serialized, sizeof(serialized), &oversized) !=
            WirehairV2_Success)
    {
        return Check(false, "oversized descriptor fixture parse");
    }
    oversized.message_bytes = static_cast<uint64_t>(SIZE_MAX) + 1u;
    oversized.block_bytes = static_cast<uint32_t>(
        oversized.message_bytes / BlockCount +
        (oversized.message_bytes % BlockCount != 0u ? 1u : 0u));
    uint32_t serialized_bytes = 0u;
    if (wirehair_v2_profile_serialize(
            &oversized, serialized, sizeof(serialized),
            &serialized_bytes) != WirehairV2_Success ||
        serialized_bytes != ProfileBytes)
    {
        return Check(false, "oversized descriptor fixture serialize");
    }
    codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    if (wirehair_v2_encoder_create_profile_with_options(
            reinterpret_cast<const void*>(uintptr_t(1)),
            serialized, sizeof(serialized), &options, &codec) !=
                WirehairV2_InvalidInput ||
        codec != reinterpret_cast<WirehairV2Codec>(uintptr_t(1)))
    {
        return Check(false, "narrow-platform oversized borrowed range");
    }
#endif
    return true;
}

bool CheckAttachedEncodeOverlap()
{
    std::vector<uint8_t> source(BlockCount * BlockBytes);
    FillMessage(source);
    const std::vector<uint8_t> original = source;
    WirehairV2EncoderOptions options = MakeOptions(
        static_cast<uint32_t>(
            WirehairV2EncoderSource_BorrowedImmutable));
    uint8_t descriptor[ProfileBytes];
    uint32_t descriptor_bytes = 0u;
    CodecOwner owner;
    if (wirehair_v2_encoder_create_with_options(
            source.data(), source.size(), BlockBytes, &options,
            descriptor, sizeof(descriptor),
            &descriptor_bytes, &owner.Handle) != WirehairV2_Success ||
        !owner.Handle || descriptor_bytes != ProfileBytes)
    {
        return Check(false, "attached overlap fixture creation");
    }

    struct OutputCase
    {
        uint32_t Id;
        size_t Offset;
        uint32_t Capacity;
    };
    const OutputCase cases[] = {
        { 0u, 0u, BlockBytes },
        { 0u, 1u, BlockBytes },
        { 0u, 2u * BlockBytes, BlockBytes },
        { BlockCount, source.size() - 1u, BlockBytes },
        { UINT32_MAX, 3u * BlockBytes, BlockBytes - 1u }
    };
    for (const OutputCase& item : cases)
    {
        uint32_t bytes = UINT32_C(0xa5a5a5a5);
        const WirehairV2Result result = wirehair_v2_encode(
            owner.Handle, item.Id,
            source.data() + item.Offset, item.Capacity, &bytes);
        if (result != WirehairV2_InvalidInput ||
            bytes != UINT32_C(0xa5a5a5a5) || source != original)
        {
            return Check(
                false, "attached packet/source overlap is no-write");
        }
    }

    for (size_t counter_offset : { size_t(32u), size_t(36u) })
    {
        const WirehairV2Result result = wirehair_v2_encode(
            owner.Handle, 0u, nullptr, 0u,
            reinterpret_cast<uint32_t*>(source.data() + counter_offset));
        if (result != WirehairV2_InvalidInput || source != original) {
            return Check(
                false, "attached counter/source overlap precedes null output");
        }
    }

    {
        alignas(std::max_align_t)
            uint8_t backing[BlockCount * BlockBytes + 8u];
        std::memset(backing, 0xa5, sizeof(backing));
        uint8_t* const shifted_source = backing + 2u;
        for (size_t i = 0u; i < BlockCount * BlockBytes; ++i) {
            shifted_source[i] = static_cast<uint8_t>(i * 73u + 19u);
        }
        uint8_t shifted_descriptor[ProfileBytes];
        uint32_t shifted_descriptor_bytes = 0u;
        CodecOwner shifted_owner;
        if (wirehair_v2_encoder_create_with_options(
                shifted_source, BlockCount * BlockBytes, BlockBytes,
                &options, shifted_descriptor, sizeof(shifted_descriptor),
                &shifted_descriptor_bytes, &shifted_owner.Handle) !=
                    WirehairV2_Success ||
            !shifted_owner.Handle || shifted_descriptor_bytes != ProfileBytes)
        {
            return Check(false, "partial counter/source overlap fixture");
        }
        uint8_t before[sizeof(backing)];
        std::memcpy(before, backing, sizeof(before));
        uint32_t* const straddling_counter =
            reinterpret_cast<uint32_t*>(
                backing + BlockCount * BlockBytes);
        if (wirehair_v2_encode(
                shifted_owner.Handle, 0u, nullptr, 0u,
                straddling_counter) != WirehairV2_InvalidInput ||
            std::memcmp(backing, before, sizeof(backing)) != 0)
        {
            return Check(
                false, "partial-boundary counter/source overlap is no-write");
        }
    }

    uint8_t disjoint[BlockBytes];
    std::memset(disjoint, 0xa5, sizeof(disjoint));
    const std::array<uint8_t, BlockBytes> disjoint_before = [&]() {
        std::array<uint8_t, BlockBytes> copy{};
        std::copy(disjoint, disjoint + BlockBytes, copy.begin());
        return copy;
    }();
    uint32_t bytes = UINT32_C(0xa5a5a5a5);
    if (wirehair_v2_encode(
            owner.Handle, BlockCount, disjoint, BlockBytes - 1u, &bytes) !=
                WirehairV2_BufferTooSmall ||
        bytes != BlockBytes ||
        std::memcmp(
            disjoint, disjoint_before.data(), sizeof(disjoint)) != 0)
    {
        return Check(false, "borrowed short packet remains transactional");
    }

    bytes = UINT32_C(0xa5a5a5a5);
    const uintptr_t limit = std::numeric_limits<uintptr_t>::max();
    if (wirehair_v2_encode(
            owner.Handle, 0u,
            reinterpret_cast<void*>(limit - 7u), BlockBytes, &bytes) !=
                WirehairV2_InvalidInput ||
        bytes != UINT32_C(0xa5a5a5a5) || source != original)
    {
        return Check(false, "wrapping packet range is no-write");
    }

    std::memset(disjoint, 0xa5, sizeof(disjoint));
    if (wirehair_v2_encode(
            owner.Handle, BlockCount,
            disjoint, sizeof(disjoint),
            reinterpret_cast<uint32_t*>(limit - 1u)) !=
                WirehairV2_InvalidInput ||
        std::memcmp(
            disjoint, disjoint_before.data(), sizeof(disjoint)) != 0 ||
        source != original)
    {
        return Check(false, "wrapping counter range is no-write");
    }
    return true;
}

bool CheckDetachedFormerSourceOverlap()
{
    std::vector<uint8_t> original(BlockCount * BlockBytes);
    FillMessage(original);
    std::array<uint8_t, ProfileBytes> descriptor{};
    PacketSet expected{};
    if (!BuildReference(original, descriptor, expected)) {
        return Check(false, "detached former-source overlap reference");
    }

    std::vector<uint8_t> source = original;
    WirehairV2EncoderOptions options = MakeOptions(
        static_cast<uint32_t>(
            WirehairV2EncoderSource_BorrowedImmutable));
    uint32_t descriptor_bytes = 0u;
    CodecOwner owner;
    if (wirehair_v2_encoder_create_with_options(
            source.data(), source.size(), BlockBytes, &options,
            descriptor.data(), descriptor.size(),
            &descriptor_bytes, &owner.Handle) != WirehairV2_Success ||
        !owner.Handle || descriptor_bytes != ProfileBytes ||
        wirehair_v2_encoder_detach_input(owner.Handle) != WirehairV2_Success)
    {
        return Check(false, "detached former-source overlap fixture");
    }

    struct Case
    {
        unsigned Probe;
        size_t Offset;
    };
    const Case cases[] = {
        { 0u, 0u },
        { BlockCount - 1u, 1u },
        { BlockCount, 2u * BlockBytes },
        { ProbeCount - 1u, 4u * BlockBytes }
    };
    for (const Case& item : cases)
    {
        uint32_t bytes = UINT32_C(0xa5a5a5a5);
        const PacketCapture& packet = expected[item.Probe];
        const WirehairV2Result result = wirehair_v2_encode(
            owner.Handle, packet.Id,
            source.data() + item.Offset, BlockBytes, &bytes);
        if (result != WirehairV2_Success || bytes != packet.Bytes ||
            std::memcmp(
                source.data() + item.Offset,
                packet.Data.data(), packet.Bytes) != 0)
        {
            std::fprintf(stderr,
                "V2 detached former-source overlap failed: "
                "probe=%u id=%u result=%d\n",
                item.Probe, packet.Id, static_cast<int>(result));
            return false;
        }
    }
    return true;
}

bool CheckDetachHandleKinds(
    const std::array<uint8_t, ProfileBytes>& descriptor)
{
    if (wirehair_v2_encoder_detach_input(nullptr) !=
            WirehairV2_InvalidInput)
    {
        return Check(false, "null detach result");
    }
    CodecOwner decoder;
    if (wirehair_v2_decoder_create(
            descriptor.data(), descriptor.size(), &decoder.Handle) !=
                WirehairV2_Success ||
        wirehair_v2_encoder_detach_input(decoder.Handle) !=
            WirehairV2_InvalidInput)
    {
        return Check(false, "decoder detach result");
    }
    return true;
}

WirehairV2Result CreateCppBorrowed(
    ConstructorRoute route,
    wirehair::v2::Encoder& encoder,
    const std::vector<uint8_t>& source,
    wirehair::v2::SerializedProfile& profile)
{
    switch (route)
    {
    case ConstructorRoute::DefaultSelector:
        return encoder.CreateBorrowed(
            source.data(), source.size(), BlockBytes, profile);
    case ConstructorRoute::ExplicitSelector:
        return encoder.CreateBorrowed(
            WIREHAIR_V2_PROFILE_CURRENT,
            source.data(), source.size(), BlockBytes, profile);
    case ConstructorRoute::Descriptor:
        return encoder.CreateBorrowed(source.data(), profile);
    default:
        return WirehairV2_Error;
    }
}

bool ReplayCppPackets(
    wirehair::v2::Encoder& encoder,
    const PacketSet& expected,
    const char* what)
{
    static const uint8_t Guard = 0x5au;
    for (unsigned probe = 0u; probe < ProbeCount; ++probe)
    {
        uint8_t guarded[BlockBytes + 2u];
        std::memset(guarded, Guard, sizeof(guarded));
        uint32_t bytes = UINT32_C(0xa5a5a5a5);
        const WirehairV2Result result = encoder.Encode(
            expected[probe].Id,
            guarded + 1u, BlockBytes, bytes);
        bool bounds_ok = guarded[0] == Guard &&
            guarded[sizeof(guarded) - 1u] == Guard;
        for (uint32_t i = bytes; bounds_ok && i < BlockBytes; ++i) {
            bounds_ok = guarded[1u + i] == Guard;
        }
        if (result != WirehairV2_Success ||
            bytes != expected[probe].Bytes || !bounds_ok ||
            std::memcmp(
                guarded + 1u, expected[probe].Data.data(), bytes) != 0)
        {
            std::fprintf(stderr,
                "V2 borrowed-source C++ replay failed: %s "
                "probe=%u id=%u result=%d\n",
                what, probe, expected[probe].Id,
                static_cast<int>(result));
            return false;
        }
    }
    return true;
}

bool CheckCppBorrowedOwnership()
{
    std::vector<uint8_t> original(BlockCount * BlockBytes);
    FillMessage(original);
    std::array<uint8_t, ProfileBytes> expected_descriptor{};
    PacketSet expected_packets{};
    if (!BuildReference(
            original, expected_descriptor, expected_packets))
    {
        return Check(false, "C++ ownership reference");
    }

    for (unsigned route_value = 0u;
         route_value < static_cast<unsigned>(ConstructorRoute::Count);
         ++route_value)
    {
        const ConstructorRoute route =
            static_cast<ConstructorRoute>(route_value);
        std::vector<uint8_t> source = original;
        wirehair::v2::SerializedProfile profile;
        if (route == ConstructorRoute::Descriptor) {
            std::copy(
                expected_descriptor.begin(), expected_descriptor.end(),
                profile.data());
        }
        wirehair::v2::Encoder encoder;
        if (CreateCppBorrowed(route, encoder, source, profile) !=
                WirehairV2_Success ||
            !encoder ||
            std::memcmp(
                profile.data(), expected_descriptor.data(),
                ProfileBytes) != 0 ||
            !ReplayCppPackets(
                encoder, expected_packets,
                "initial borrowed C++ encoder"))
        {
            std::fprintf(stderr,
                "V2 borrowed-source C++ create failed: route=%s\n",
                RouteName(route));
            return false;
        }

        std::vector<uint8_t> failed_source = original;
        FillMessage(failed_source, 0x33u);
        wirehair::v2::SerializedProfile failed_profile;
        std::memset(
            failed_profile.data(), 0xa5, failed_profile.size());
        WirehairV2Result failed_result = WirehairV2_Error;
        WirehairV2Result expected_failure = WirehairV2_Error;
        switch (route)
        {
        case ConstructorRoute::DefaultSelector:
            failed_result = encoder.CreateBorrowed(
                failed_source.data(), failed_source.size(), 0u,
                failed_profile);
            expected_failure = WirehairV2_InvalidDimensions;
            break;
        case ConstructorRoute::ExplicitSelector:
            failed_result = encoder.CreateBorrowed(
                UINT64_MAX,
                failed_source.data(), failed_source.size(), BlockBytes,
                failed_profile);
            expected_failure = WirehairV2_UnsupportedProfile;
            break;
        case ConstructorRoute::Descriptor:
            std::copy(
                expected_descriptor.begin(), expected_descriptor.end(),
                failed_profile.data());
            failed_profile.data()[ProfileBytes - 1u] = 1u;
            failed_result = encoder.CreateBorrowed(
                failed_source.data(), failed_profile);
            expected_failure = WirehairV2_ReservedNonzero;
            break;
        default:
            return false;
        }
        if (failed_result != expected_failure || !encoder)
        {
            return Check(
                false, "C++ borrowed failed replacement transactionality");
        }
        const std::vector<uint8_t> failed_original = failed_source;
        PoisonAndRelease(failed_source, failed_original);
        if (!ReplayCppPackets(
                encoder, expected_packets,
                "failed borrowed replacement releases candidate"))
        {
            return false;
        }

        std::vector<uint8_t> replacement_source(original.size());
        FillMessage(replacement_source, 0x71u);
        std::array<uint8_t, ProfileBytes> replacement_descriptor{};
        PacketSet replacement_packets{};
        if (!BuildReference(
                replacement_source,
                replacement_descriptor, replacement_packets))
        {
            return Check(false, "C++ replacement reference");
        }
        wirehair::v2::SerializedProfile replacement_profile;
        if (route == ConstructorRoute::Descriptor) {
            std::copy(
                replacement_descriptor.begin(),
                replacement_descriptor.end(),
                replacement_profile.data());
        }
        if (CreateCppBorrowed(
                route, encoder,
                replacement_source, replacement_profile) !=
                    WirehairV2_Success ||
            std::memcmp(
                replacement_profile.data(),
                replacement_descriptor.data(), ProfileBytes) != 0 ||
            !ReplayCppPackets(
                encoder, replacement_packets,
                "successful borrowed replacement"))
        {
            return Check(false, "C++ borrowed successful replacement");
        }
        PoisonAndRelease(source, original);

        wirehair::v2::Encoder moved(std::move(encoder));
        if (encoder || !moved ||
            !ReplayCppPackets(
                moved, replacement_packets,
                "borrowed move construction"))
        {
            return Check(false, "C++ borrowed move construction ownership");
        }
        uint8_t invalid_output[BlockBytes + 2u];
        std::memset(invalid_output, 0xa5, sizeof(invalid_output));
        const std::array<uint8_t, BlockBytes + 2u> invalid_before = [&]() {
            std::array<uint8_t, BlockBytes + 2u> copy{};
            std::copy(
                invalid_output, invalid_output + sizeof(invalid_output),
                copy.begin());
            return copy;
        }();
        uint32_t invalid_bytes = UINT32_C(0xa5a5a5a5);
        if (encoder.Encode(
                0u, invalid_output + 1u, BlockBytes, invalid_bytes) !=
                    WirehairV2_InvalidInput ||
            invalid_bytes != 0u ||
            std::memcmp(
                invalid_output, invalid_before.data(),
                sizeof(invalid_output)) != 0)
        {
            return Check(false, "C++ moved-from encoder is invalid");
        }

        std::vector<uint8_t> destination_source(original.size());
        FillMessage(destination_source, 0xb3u);
        wirehair::v2::SerializedProfile destination_profile;
        wirehair::v2::Encoder destination;
        if (destination.CreateBorrowed(
                destination_source.data(), destination_source.size(),
                BlockBytes, destination_profile) != WirehairV2_Success ||
            !destination)
        {
            return Check(false, "C++ move-assignment destination fixture");
        }
        const std::vector<uint8_t> destination_original = destination_source;
        destination = std::move(moved);
        PoisonAndRelease(destination_source, destination_original);
        if (!destination || moved ||
            !ReplayCppPackets(
                destination, replacement_packets,
                "borrowed move assignment"))
        {
            return Check(false, "C++ borrowed move assignment ownership");
        }
        if (destination.DetachInput() != WirehairV2_Success ||
            destination.DetachInput() != WirehairV2_Success)
        {
            return Check(false, "C++ borrowed detach idempotence");
        }
        const std::vector<uint8_t> replacement_original = replacement_source;
        PoisonAndRelease(replacement_source, replacement_original);
        if (!ReplayCppPackets(
                destination, replacement_packets,
                "C++ detached source release"))
        {
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    std::vector<uint8_t> message(BlockCount * BlockBytes);
    FillMessage(message);
    std::array<uint8_t, ProfileBytes> descriptor{};
    PacketSet packets{};
    if (!BuildReference(message, descriptor, packets) ||
        !CheckConstructorAndLifetimeMatrix() ||
        !CheckSelectingOptionPriority() ||
        !CheckDescriptorOptionPriority(message, descriptor) ||
        !CheckSelectingOptionsOverlap() ||
        !CheckSelectingFixedOverlaps() ||
        !CheckBorrowedDescriptorMessageOverlap(descriptor, packets) ||
        !CheckEmptyBorrowedMessageDoesNotOverlapDescriptor() ||
        !CheckReadOnlyInputOverlaps(descriptor) ||
        !CheckDescriptorFixedOverlaps(descriptor) ||
        !CheckWrappedBorrowedMessageRange(descriptor) ||
        !CheckAttachedEncodeOverlap() ||
        !CheckDetachedFormerSourceOverlap() ||
        !CheckDetachHandleKinds(descriptor) ||
        !CheckCppBorrowedOwnership())
    {
        return 1;
    }
    std::puts("V2 borrowed-source contract test passed");
    return 0;
}
