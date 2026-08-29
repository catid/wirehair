#include "Wh2FrozenTrace.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <climits>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

#ifndef WIREHAIR_WH2_HARNESS_GIT_COMMIT
#error "facade timing worker requires an exact harness commit receipt"
#endif

#ifndef WIREHAIR_WH2_IMPLEMENTATION_GIT_COMMIT
#error "facade timing worker requires an exact implementation commit receipt"
#endif

#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT) && \
    defined(WIREHAIR_WH2_FACADE_ROLE_PARENT)
#error "facade timing worker role is ambiguous"
#endif

#if !defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT) && \
    !defined(WIREHAIR_WH2_FACADE_ROLE_PARENT)
#error "facade timing worker requires one compile-time role"
#endif

namespace {

using Clock = std::chrono::steady_clock;

static const char kCampaign[] =
    "wh2-v2-facade-default-parent-falsifier-r0";
static const char kConfigSchema[] =
    "wirehair.wh2.v2-facade-timing-worker.config.v1";
static const char kInvocationSchema[] =
    "wirehair.wh2.v2-facade-timing-worker.invocation.v1";
static const char kTerminalSchema[] =
    "wirehair.wh2.v2-facade-timing-worker.terminal.v1";
static const char kArmDescriptorSchema[] =
    "wirehair.wh2.v2-facade-timing-worker.arm.v1";

static const uint32_t kPanelReplicates = 12u;
static const uint32_t kInvocationBudget = 65536u;
static const uint32_t kMinimumInvocations = 24u;
static const uint32_t kInternalDeadlineSeconds = 840u;
static const uint64_t kSourceSeedBase = UINT64_C(0xb6402ee71c8365a9);

struct CellShape
{
    uint32_t K;
    uint32_t BlockBytes;
};

static const CellShape kCellShapes[] = {
    { 8u, 64u }, { 8u, 1280u },
    { 128u, 64u }, { 128u, 1280u },
    { 512u, 64u }, { 512u, 1280u },
    { 5000u, 64u }, { 5000u, 1280u },
    { 64000u, 64u }, { 64000u, 1280u }
};

enum class Arm
{
    C,
    E,
    D,
    P,
    L,
    Invalid
};

enum class Scope
{
    PrebuiltSystematic,
    FreshSystematic,
    FreshRepair,
    PrebuiltRepair,
    Invalid
};

const char* CompiledRoleName()
{
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
    return "current";
#else
    return "parent";
#endif
}

const char* ArmName(Arm arm)
{
    switch (arm)
    {
    case Arm::C: return "C";
    case Arm::E: return "E";
    case Arm::D: return "D";
    case Arm::P: return "P";
    case Arm::L: return "L";
    default: return nullptr;
    }
}

Arm ParseArm(const std::string& text)
{
    if (text == "C") return Arm::C;
    if (text == "E") return Arm::E;
    if (text == "D") return Arm::D;
    if (text == "P") return Arm::P;
    if (text == "L") return Arm::L;
    return Arm::Invalid;
}

bool ArmSupported(Arm arm)
{
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
    return arm == Arm::C || arm == Arm::E || arm == Arm::D;
#else
    return arm == Arm::P || arm == Arm::L;
#endif
}

bool IsV2Arm(Arm arm)
{
    return arm == Arm::C || arm == Arm::E || arm == Arm::D || arm == Arm::P;
}

const char* ScopeName(Scope scope)
{
    switch (scope)
    {
    case Scope::PrebuiltSystematic: return "prebuilt-systematic";
    case Scope::FreshSystematic: return "fresh-systematic";
    case Scope::FreshRepair: return "fresh-repair";
    case Scope::PrebuiltRepair: return "prebuilt-repair";
    default: return nullptr;
    }
}

Scope ParseScope(const std::string& text)
{
    if (text == "prebuilt-systematic") {
        return Scope::PrebuiltSystematic;
    }
    if (text == "fresh-systematic") return Scope::FreshSystematic;
    if (text == "fresh-repair") return Scope::FreshRepair;
    if (text == "prebuilt-repair") return Scope::PrebuiltRepair;
    return Scope::Invalid;
}

bool IsSystematicScope(Scope scope)
{
    return scope == Scope::PrebuiltSystematic ||
        scope == Scope::FreshSystematic;
}

bool IsFreshScope(Scope scope)
{
    return scope == Scope::FreshSystematic || scope == Scope::FreshRepair;
}

bool CheckedProduct(uint64_t left, uint64_t right, size_t& product_out)
{
    product_out = 0u;
    if (left != 0u && right >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()) / left)
    {
        return false;
    }
    product_out = static_cast<size_t>(left * right);
    return true;
}

uint64_t SplitMix64(uint64_t& state)
{
    state += UINT64_C(0x9e3779b97f4a7c15);
    uint64_t value = state;
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

uint64_t SourceSeed(const CellShape& shape)
{
    return kSourceSeedBase ^
        (static_cast<uint64_t>(shape.K) << 17) ^ shape.BlockBytes;
}

std::string Hex64(uint64_t value)
{
    static const char hex[] = "0123456789abcdef";
    std::string result("0x0000000000000000");
    for (unsigned i = 0u; i < 16u; ++i)
    {
        result[17u - i] = hex[value & 15u];
        value >>= 4;
    }
    return result;
}

std::string HexBytes(const void* data, size_t bytes)
{
    if (!data && bytes != 0u) return std::string();
    static const char hex[] = "0123456789abcdef";
    const uint8_t* const source = static_cast<const uint8_t*>(data);
    std::string result;
    try {
        result.resize(bytes * 2u);
    }
    catch (...) {
        return std::string();
    }
    for (size_t i = 0u; i < bytes; ++i)
    {
        result[i * 2u] = hex[source[i] >> 4];
        result[i * 2u + 1u] = hex[source[i] & 15u];
    }
    return result;
}

bool MakeSource(
    const CellShape& shape,
    uint64_t seed,
    std::vector<uint8_t>& source_out)
{
    size_t bytes = 0u;
    if (shape.K < 2u || shape.K > 64000u || shape.BlockBytes == 0u ||
        !CheckedProduct(shape.K, shape.BlockBytes, bytes))
    {
        return false;
    }
    try
    {
        std::vector<uint8_t> source(bytes, 0u);
        uint64_t state = seed ^
            (static_cast<uint64_t>(shape.K) *
                UINT64_C(0xd6e8feb86659fd93)) ^
            (static_cast<uint64_t>(shape.BlockBytes) *
                UINT64_C(0xa0761d6478bd642f));
        size_t offset = 0u;
        while (offset < source.size())
        {
            const uint64_t word = SplitMix64(state);
            for (unsigned byte = 0u; byte < 8u && offset < source.size();
                 ++byte, ++offset)
            {
                source[offset] = static_cast<uint8_t>(word >> (byte * 8u));
            }
        }
        source_out.swap(source);
        return true;
    }
    catch (...) {
        return false;
    }
}

uint32_t TotalInvocations(uint32_t K)
{
    const uint32_t budget = (kInvocationBudget + K - 1u) / K;
    return std::max(kMinimumInvocations, budget);
}

uint32_t InvocationsForReplicate(uint32_t K, uint32_t replicate)
{
    const uint32_t total = TotalInvocations(K);
    return total / kPanelReplicates +
        (replicate < total % kPanelReplicates ? 1u : 0u);
}

uint64_t PositiveNanoseconds(Clock::duration duration)
{
    const std::chrono::nanoseconds elapsed =
        std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
    return elapsed.count() > 0 ? static_cast<uint64_t>(elapsed.count()) : 0u;
}

class EncoderHandle
{
public:
    EncoderHandle() = default;
    ~EncoderHandle() { Reset(); }

    EncoderHandle(const EncoderHandle&) = delete;
    EncoderHandle& operator=(const EncoderHandle&) = delete;

    EncoderHandle(EncoderHandle&& other) noexcept
        : V2(other.V2)
        , Legacy(other.Legacy)
        , IsV2(other.IsV2)
    {
        other.V2 = nullptr;
        other.Legacy = nullptr;
        other.IsV2 = false;
    }

    EncoderHandle& operator=(EncoderHandle&& other) noexcept
    {
        if (this != &other)
        {
            Reset();
            V2 = other.V2;
            Legacy = other.Legacy;
            IsV2 = other.IsV2;
            other.V2 = nullptr;
            other.Legacy = nullptr;
            other.IsV2 = false;
        }
        return *this;
    }

    void AdoptV2(WirehairV2Codec codec)
    {
        V2 = codec;
        IsV2 = true;
    }

    void AdoptLegacy(WirehairCodec codec)
    {
        Legacy = codec;
        IsV2 = false;
    }

    bool Valid() const { return IsV2 ? V2 != nullptr : Legacy != nullptr; }
    bool UsesV2() const { return IsV2; }
    WirehairV2Codec V2Codec() const { return V2; }
    WirehairCodec LegacyCodec() const { return Legacy; }

    void Reset()
    {
        if (V2) wirehair_v2_free(V2);
        if (Legacy) wirehair_free(Legacy);
        V2 = nullptr;
        Legacy = nullptr;
        IsV2 = false;
    }

private:
    WirehairV2Codec V2 = nullptr;
    WirehairCodec Legacy = nullptr;
    bool IsV2 = false;
};

struct ConstructionArtifacts
{
    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> Profile = {{0u}};
    uint32_t ProfileBytes = 0u;
};

class ConstructionRequest
{
public:
    ConstructionRequest(
        Arm arm,
        const void* source,
        uint64_t message_bytes,
        uint32_t block_bytes)
        : ArmValue(arm)
        , Source(source)
        , MessageBytes(message_bytes)
        , BlockBytes(block_bytes)
    {
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
        if (arm == Arm::C || arm == Arm::E)
        {
            Options.struct_bytes = sizeof(Options);
            Options.options_version = WIREHAIR_V2_ENCODER_OPTIONS_VERSION;
            Options.source_policy = arm == Arm::C ?
                WirehairV2EncoderSource_BorrowedImmutable :
                WirehairV2EncoderSource_Independent;
            Options.reserved = 0u;
        }
#endif
    }

    int64_t Execute(
        EncoderHandle& encoder,
        ConstructionArtifacts& artifacts) const
    {
        if (encoder.Valid()) return WirehairV2_InvalidInput;

        WirehairV2Codec codec = nullptr;
        WirehairV2Result result = WirehairV2_InvalidInput;
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
        if (ArmValue == Arm::D)
#else
        if (ArmValue == Arm::P)
#endif
        {
            result = wirehair_v2_encoder_create(
                Source,
                MessageBytes,
                BlockBytes,
                artifacts.Profile.data(),
                static_cast<uint32_t>(artifacts.Profile.size()),
                &artifacts.ProfileBytes,
                &codec);
        }
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
        else if (ArmValue == Arm::C || ArmValue == Arm::E)
        {
            result = wirehair_v2_encoder_create_with_options(
                Source,
                MessageBytes,
                BlockBytes,
                &Options,
                artifacts.Profile.data(),
                static_cast<uint32_t>(artifacts.Profile.size()),
                &artifacts.ProfileBytes,
                &codec);
        }
#else
        else if (ArmValue == Arm::L)
        {
            WirehairCodec legacy = nullptr;
            const WirehairResult legacy_result = wirehair_encoder_create_ex(
                nullptr, Source, MessageBytes, BlockBytes, &legacy);
            if (legacy_result != Wirehair_Success || !legacy) {
                wirehair_free(legacy);
                return legacy_result;
            }
            encoder.AdoptLegacy(legacy);
            return legacy_result;
        }
#endif
        else return WirehairV2_InvalidInput;
        if (result != WirehairV2_Success || !codec)
        {
            wirehair_v2_free(codec);
            return result;
        }
        encoder.AdoptV2(codec);
        return result;
    }

private:
    Arm ArmValue;
    const void* Source;
    uint64_t MessageBytes;
    uint32_t BlockBytes;
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
    WirehairV2EncoderOptions Options;
#endif
};

struct OracleReceipt
{
    std::string ArmDescriptorSha256;
    std::string EquationConfigurationSha256;
    std::string FirstRepairSha256;
    std::string HighIdSha256;
    std::string PublicStateReceiptSha256;
    std::string RepairSha256;
    std::string RoundtripSha256;
    std::string SerializedProfileHex;
    std::string SerializedProfileSha256;
    std::string SystematicSha256;
    uint32_t SerializedProfileBytes = 0u;
    bool RoundtripVerified = false;
};

struct ArmState
{
    Arm ArmValue = Arm::Invalid;
    EncoderHandle Prebuilt;
    ConstructionArtifacts Artifacts;
    OracleReceipt Oracle;
};

struct Cell
{
    CellShape Shape = {};
    uint64_t Seed = 0u;
    std::vector<uint8_t> Source;
    std::vector<uint8_t> Output;
    std::string SourceSha256;
    std::vector<ArmState> Arms;
};

std::string ArmDescriptorJson(Arm arm)
{
    const char* api = nullptr;
    const char* codec = nullptr;
    const char* equation_profile = nullptr;
    const char* source_policy = nullptr;
    switch (arm)
    {
    case Arm::C:
        api = "wirehair_v2_encoder_create_with_options";
        codec = "wirehair2";
        equation_profile = "certified-2026-07";
        source_policy = "borrowed-immutable";
        break;
    case Arm::E:
        api = "wirehair_v2_encoder_create_with_options";
        codec = "wirehair2";
        equation_profile = "certified-2026-07";
        source_policy = "independent";
        break;
    case Arm::D:
        api = "wirehair_v2_encoder_create";
        codec = "wirehair2";
        equation_profile = "certified-2026-07";
        source_policy = "default-independent";
        break;
    case Arm::P:
        api = "wirehair_v2_encoder_create";
        codec = "wirehair2";
        equation_profile = "certified-2026-07";
        source_policy = "pre-feature-independent";
        break;
    case Arm::L:
        api = "wirehair_encoder_create_ex";
        codec = "wirehair1";
        equation_profile = "legacy-current";
        source_policy = "borrowed";
        break;
    default:
        return std::string();
    }
    std::ostringstream json;
    json << "{\"api\":\"" << api
         << "\",\"arm\":\"" << ArmName(arm)
         << "\",\"codec\":\"" << codec
         << "\",\"equation_profile\":\"" << equation_profile
         << "\",\"schema\":\"" << kArmDescriptorSchema
         << "\",\"source_policy\":\"" << source_policy << "\"}";
    return json.str();
}

bool ValidateV2Profile(
    const Cell& cell,
    const ConstructionArtifacts& artifacts)
{
    if (artifacts.ProfileBytes != artifacts.Profile.size() ||
        wirehair_v2_profile_validate(
            artifacts.Profile.data(), artifacts.ProfileBytes) !=
            WirehairV2_Success)
    {
        return false;
    }
    WirehairV2Profile profile = {};
    if (wirehair_v2_profile_deserialize(
            artifacts.Profile.data(), artifacts.ProfileBytes, &profile) !=
            WirehairV2_Success)
    {
        return false;
    }
    return profile.message_bytes == cell.Source.size() &&
        profile.block_bytes == cell.Shape.BlockBytes;
}

int64_t EncodeSweep(
    const EncoderHandle& encoder,
    uint32_t first_id,
    uint32_t count,
    uint32_t block_bytes,
    uint8_t* output)
{
    if (!encoder.Valid() || !output) return Wirehair_Error;
    if (encoder.UsesV2())
    {
        WirehairV2Result result = WirehairV2_Success;
        for (uint32_t offset = 0u; offset < count; ++offset)
        {
            uint32_t bytes = 0u;
            result = wirehair_v2_encode(
                encoder.V2Codec(),
                first_id + offset,
                output + static_cast<size_t>(offset) * block_bytes,
                block_bytes,
                &bytes);
            if (result != WirehairV2_Success || bytes != block_bytes) {
                return result == WirehairV2_Success ?
                    WirehairV2_Error : result;
            }
        }
        return result;
    }

    WirehairResult result = Wirehair_Success;
    for (uint32_t offset = 0u; offset < count; ++offset)
    {
        uint32_t bytes = 0u;
        result = wirehair_encode(
            encoder.LegacyCodec(),
            first_id + offset,
            output + static_cast<size_t>(offset) * block_bytes,
            block_bytes,
            &bytes);
        if (result != Wirehair_Success || bytes != block_bytes) {
            return result == Wirehair_Success ? Wirehair_Error : result;
        }
    }
    return result;
}

int64_t EncodeOne(
    const EncoderHandle& encoder,
    uint32_t id,
    uint32_t block_bytes,
    uint8_t* output)
{
    return EncodeSweep(encoder, id, 1u, block_bytes, output);
}

bool VerifyPublicRoundtrip(
    Cell& cell,
    const ArmState& state,
    OracleReceipt& oracle)
{
    if (!state.Prebuilt.Valid()) return false;
    std::vector<uint8_t> recovered;
    try {
        recovered.assign(cell.Source.size(), 0u);
    }
    catch (...) {
        return false;
    }

    if (state.Prebuilt.UsesV2())
    {
        WirehairV2Codec decoder = nullptr;
        if (wirehair_v2_decoder_create(
                state.Artifacts.Profile.data(),
                state.Artifacts.ProfileBytes,
                &decoder) != WirehairV2_Success || !decoder)
        {
            return false;
        }
        WirehairV2Result result = WirehairV2_NeedMore;
        for (uint32_t id = 0u; id < cell.Shape.K; ++id)
        {
            uint32_t bytes = 0u;
            uint8_t* const packet = cell.Output.data() +
                static_cast<size_t>(id) * cell.Shape.BlockBytes;
            result = wirehair_v2_encode(
                state.Prebuilt.V2Codec(), id, packet,
                cell.Shape.BlockBytes, &bytes);
            if (result != WirehairV2_Success ||
                bytes != cell.Shape.BlockBytes)
            {
                wirehair_v2_free(decoder);
                return false;
            }
            result = wirehair_v2_decode(decoder, id, packet, bytes);
            if (result != WirehairV2_NeedMore &&
                result != WirehairV2_Success)
            {
                wirehair_v2_free(decoder);
                return false;
            }
        }
        uint64_t bytes_out = 0u;
        const bool complete = result == WirehairV2_Success &&
            wirehair_v2_recover(
                decoder, recovered.data(), recovered.size(), &bytes_out) ==
                WirehairV2_Success && bytes_out == recovered.size();
        wirehair_v2_free(decoder);
        if (!complete) return false;
    }
    else
    {
        WirehairCodec decoder = wirehair_decoder_create(
            nullptr, cell.Source.size(), cell.Shape.BlockBytes);
        if (!decoder) return false;
        WirehairResult result = Wirehair_NeedMore;
        for (uint32_t id = 0u; id < cell.Shape.K; ++id)
        {
            uint32_t bytes = 0u;
            uint8_t* const packet = cell.Output.data() +
                static_cast<size_t>(id) * cell.Shape.BlockBytes;
            result = wirehair_encode(
                state.Prebuilt.LegacyCodec(), id, packet,
                cell.Shape.BlockBytes, &bytes);
            if (result != Wirehair_Success ||
                bytes != cell.Shape.BlockBytes)
            {
                wirehair_free(decoder);
                return false;
            }
            result = wirehair_decode(decoder, id, packet, bytes);
            if (result != Wirehair_NeedMore && result != Wirehair_Success)
            {
                wirehair_free(decoder);
                return false;
            }
        }
        const bool complete = result == Wirehair_Success &&
            wirehair_recover(decoder, recovered.data(), recovered.size()) ==
                Wirehair_Success;
        wirehair_free(decoder);
        if (!complete) return false;
    }
    if (recovered != cell.Source) return false;
    oracle.RoundtripSha256 =
        wirehair::wh2_benchmark::Sha256Hex(
            recovered.data(), recovered.size());
    oracle.RoundtripVerified =
        oracle.RoundtripSha256 == cell.SourceSha256;
    return oracle.RoundtripVerified;
}

std::string PublicStateJson(const Cell& cell, const OracleReceipt& oracle)
{
    std::ostringstream json;
    json << "{\"equation_configuration_sha256\":\""
         << oracle.EquationConfigurationSha256
         << "\",\"first_repair_sha256\":\""
         << oracle.FirstRepairSha256
         << "\",\"high_id_sha256\":\"" << oracle.HighIdSha256
         << "\",\"repair_sha256\":\"" << oracle.RepairSha256
         << "\",\"roundtrip_sha256\":\"" << oracle.RoundtripSha256
         << "\",\"source_sha256\":\"" << cell.SourceSha256
         << "\",\"systematic_sha256\":\"" << oracle.SystematicSha256
         << "\"}";
    return json.str();
}

bool BuildOracle(Cell& cell, ArmState& state)
{
    OracleReceipt oracle;
    const std::string descriptor = ArmDescriptorJson(state.ArmValue);
    oracle.ArmDescriptorSha256 =
        wirehair::wh2_benchmark::Sha256Hex(descriptor);
    if (oracle.ArmDescriptorSha256.size() != 64u) return false;

    if (IsV2Arm(state.ArmValue))
    {
        if (!ValidateV2Profile(cell, state.Artifacts)) return false;
        oracle.SerializedProfileBytes = state.Artifacts.ProfileBytes;
        oracle.SerializedProfileSha256 =
            wirehair::wh2_benchmark::Sha256Hex(
                state.Artifacts.Profile.data(), state.Artifacts.ProfileBytes);
        oracle.SerializedProfileHex = HexBytes(
            state.Artifacts.Profile.data(), state.Artifacts.ProfileBytes);
        oracle.EquationConfigurationSha256 =
            oracle.SerializedProfileSha256;
        if (oracle.SerializedProfileHex.size() !=
            static_cast<size_t>(state.Artifacts.ProfileBytes) * 2u)
        {
            return false;
        }
    }
    else
    {
        std::ostringstream configuration;
        configuration << "{\"K\":" << cell.Shape.K
                      << ",\"block_bytes\":" << cell.Shape.BlockBytes
                      << ",\"codec\":\"wirehair1\",\"message_bytes\":"
                      << cell.Source.size() << "}";
        oracle.EquationConfigurationSha256 =
            wirehair::wh2_benchmark::Sha256Hex(configuration.str());
    }
    if (oracle.EquationConfigurationSha256.size() != 64u) return false;

    if (EncodeSweep(
            state.Prebuilt, 0u, cell.Shape.K, cell.Shape.BlockBytes,
            cell.Output.data()) != 0 || cell.Output != cell.Source)
    {
        return false;
    }
    oracle.SystematicSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Output.data(), cell.Output.size());
    if (oracle.SystematicSha256 != cell.SourceSha256) return false;

    if (EncodeSweep(
            state.Prebuilt, cell.Shape.K, cell.Shape.K,
            cell.Shape.BlockBytes, cell.Output.data()) != 0)
    {
        return false;
    }
    oracle.RepairSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Output.data(), cell.Output.size());
    oracle.FirstRepairSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Output.data(), cell.Shape.BlockBytes);

    std::vector<uint8_t> high;
    try {
        high.assign(static_cast<size_t>(3u) * cell.Shape.BlockBytes, 0u);
    }
    catch (...) {
        return false;
    }
    const uint32_t high_ids[] = {
        cell.Shape.K + 7u, UINT32_C(0x80000000), UINT32_MAX
    };
    for (size_t i = 0u; i < 3u; ++i)
    {
        if (EncodeOne(
                state.Prebuilt, high_ids[i], cell.Shape.BlockBytes,
                high.data() + i * cell.Shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    std::vector<uint8_t> high_repeat;
    try {
        high_repeat.resize(high.size());
        std::transform(
            high.begin(), high.end(), high_repeat.begin(),
            [](uint8_t value) { return static_cast<uint8_t>(~value); });
    }
    catch (...) {
        return false;
    }
    for (size_t i = 0u; i < 3u; ++i)
    {
        if (EncodeOne(
                state.Prebuilt, high_ids[i], cell.Shape.BlockBytes,
                high_repeat.data() + i * cell.Shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    if (high_repeat != high) return false;
    oracle.HighIdSha256 = wirehair::wh2_benchmark::Sha256Hex(
        high.data(), high.size());
    if (oracle.RepairSha256.size() != 64u ||
        oracle.FirstRepairSha256.size() != 64u ||
        oracle.HighIdSha256.size() != 64u)
    {
        return false;
    }

    if (!VerifyPublicRoundtrip(cell, state, oracle)) return false;
    oracle.PublicStateReceiptSha256 =
        wirehair::wh2_benchmark::Sha256Hex(PublicStateJson(cell, oracle));
    if (oracle.PublicStateReceiptSha256.size() != 64u) return false;
    state.Oracle = std::move(oracle);
    return true;
}

bool BuildCellSources(std::vector<Cell>& cells_out)
{
    std::vector<Cell> cells;
    try {
        cells.reserve(sizeof(kCellShapes) / sizeof(kCellShapes[0]));
        for (const CellShape& shape : kCellShapes)
        {
            Cell cell;
            cell.Shape = shape;
            cell.Seed = SourceSeed(shape);
            if (!MakeSource(shape, cell.Seed, cell.Source)) return false;
            cell.Output.assign(cell.Source.size(), 0u);
            cell.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
                cell.Source.data(), cell.Source.size());
            if (cell.SourceSha256.size() != 64u) return false;
            cells.push_back(std::move(cell));
        }
    }
    catch (...) {
        return false;
    }
    cells_out.swap(cells);
    return true;
}

bool BuildPrebuiltAndOracles(std::vector<Cell>& cells)
{
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
    static const Arm supported[] = { Arm::C, Arm::E, Arm::D };
#else
    static const Arm supported[] = { Arm::P, Arm::L };
#endif
    for (Cell& cell : cells)
    {
        try {
            cell.Arms.reserve(sizeof(supported) / sizeof(supported[0]));
        }
        catch (...) {
            return false;
        }
        for (Arm arm : supported)
        {
            ArmState state;
            state.ArmValue = arm;
            const ConstructionRequest request(
                arm, cell.Source.data(), cell.Source.size(),
                cell.Shape.BlockBytes);
            if (request.Execute(state.Prebuilt, state.Artifacts) != 0 ||
                !state.Prebuilt.Valid())
            {
                return false;
            }
            try {
                cell.Arms.push_back(std::move(state));
            }
            catch (...) {
                return false;
            }
            if (!BuildOracle(cell, cell.Arms.back())) return false;
        }
    }
    return true;
}

ArmState* FindArmState(Cell& cell, Arm arm)
{
    for (ArmState& state : cell.Arms) {
        if (state.ArmValue == arm) return &state;
    }
    return nullptr;
}

Cell* FindCell(std::vector<Cell>& cells, uint32_t K, uint32_t block_bytes)
{
    for (Cell& cell : cells) {
        if (cell.Shape.K == K && cell.Shape.BlockBytes == block_bytes) {
            return &cell;
        }
    }
    return nullptr;
}

bool CrossbindLocalV2Oracles(const std::vector<Cell>& cells)
{
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
    for (const Cell& cell : cells)
    {
        if (cell.Arms.size() != 3u) return false;
        const OracleReceipt& expected = cell.Arms.front().Oracle;
        for (size_t i = 1u; i < cell.Arms.size(); ++i)
        {
            const OracleReceipt& actual = cell.Arms[i].Oracle;
            if (actual.EquationConfigurationSha256 !=
                    expected.EquationConfigurationSha256 ||
                actual.FirstRepairSha256 != expected.FirstRepairSha256 ||
                actual.HighIdSha256 != expected.HighIdSha256 ||
                actual.PublicStateReceiptSha256 !=
                    expected.PublicStateReceiptSha256 ||
                actual.RepairSha256 != expected.RepairSha256 ||
                actual.RoundtripSha256 != expected.RoundtripSha256 ||
                actual.RoundtripVerified != expected.RoundtripVerified ||
                actual.SerializedProfileBytes !=
                    expected.SerializedProfileBytes ||
                actual.SerializedProfileHex !=
                    expected.SerializedProfileHex ||
                actual.SerializedProfileSha256 !=
                    expected.SerializedProfileSha256 ||
                actual.SystematicSha256 != expected.SystematicSha256)
            {
                return false;
            }
        }
    }
#else
    (void)cells;
#endif
    return true;
}

bool VerifySourceIntegrity(const std::vector<Cell>& cells)
{
    for (const Cell& cell : cells)
    {
        if (cell.Source.empty() ||
            wirehair::wh2_benchmark::Sha256Hex(
                cell.Source.data(), cell.Source.size()) != cell.SourceSha256)
        {
            return false;
        }
    }
    return true;
}

struct InvocationResult
{
    int64_t Result = Wirehair_Error;
    uint64_t ElapsedNanoseconds = 0u;
    uint64_t ConstructorNanoseconds = 0u;
    uint64_t SystematicNanoseconds = 0u;
    uint64_t InitFirstRepairNanoseconds = 0u;
    uint32_t BorrowedSystematicInvocations = 0u;
    std::string CorrectnessSha256;
    std::string DescriptorSha256;
};

bool RunInvocation(
    Cell& cell,
    ArmState& state,
    Scope scope,
    InvocationResult& result_out)
{
    if (!ArmSupported(state.ArmValue) || !ScopeName(scope) ||
        !state.Prebuilt.Valid() || cell.Output.size() != cell.Source.size() ||
        cell.Output.empty())
    {
        return false;
    }
    InvocationResult output;
    const bool systematic = IsSystematicScope(scope);
    const uint32_t first_id = systematic ? 0u : cell.Shape.K;
    const std::string& expected = systematic ?
        state.Oracle.SystematicSha256 : state.Oracle.RepairSha256;
    if (expected.empty()) return false;
    // Startup leaves the systematic oracle in this buffer, and every prior
    // invocation is accepted only after its full-buffer digest matches the
    // selected oracle.  Complementing every byte outside the clock therefore
    // makes a short/no-write success impossible to hide behind stale output.
    // Each panel performs one validated warmup per participating arm before
    // any value contributes to its timed lanes.
    for (uint8_t& byte : cell.Output) byte ^= UINT8_C(0xff);
    EncoderHandle fresh;
    ConstructionArtifacts fresh_artifacts;

    if (!IsFreshScope(scope))
    {
        const Clock::time_point start = Clock::now();
        output.Result = EncodeSweep(
            state.Prebuilt, first_id, cell.Shape.K,
            cell.Shape.BlockBytes, cell.Output.data());
        const Clock::time_point finish = Clock::now();
        output.ElapsedNanoseconds = PositiveNanoseconds(finish - start);
    }
    else if (scope == Scope::FreshSystematic)
    {
        const ConstructionRequest request(
            state.ArmValue, cell.Source.data(), cell.Source.size(),
            cell.Shape.BlockBytes);
        const Clock::time_point start = Clock::now();
        output.Result = request.Execute(fresh, fresh_artifacts);
        const Clock::time_point constructed = Clock::now();
        if (output.Result == 0)
        {
            output.Result = EncodeSweep(
                fresh, 0u, cell.Shape.K, cell.Shape.BlockBytes,
                cell.Output.data());
        }
        const Clock::time_point finish = Clock::now();
        output.ConstructorNanoseconds =
            PositiveNanoseconds(constructed - start);
        output.SystematicNanoseconds =
            PositiveNanoseconds(finish - constructed);
        output.ElapsedNanoseconds = PositiveNanoseconds(finish - start);
        fresh.Reset(); // Destruction is intentionally after the stopped clock.
    }
    else
    {
        const ConstructionRequest request(
            state.ArmValue, cell.Source.data(), cell.Source.size(),
            cell.Shape.BlockBytes);
        const Clock::time_point start = Clock::now();
        output.Result = request.Execute(fresh, fresh_artifacts);
        if (output.Result == 0)
        {
            output.Result = EncodeOne(
                fresh, cell.Shape.K, cell.Shape.BlockBytes,
                cell.Output.data());
        }
        const Clock::time_point first_repair = Clock::now();
        if (output.Result == 0 && cell.Shape.K > 1u)
        {
            output.Result = EncodeSweep(
                fresh, cell.Shape.K + 1u, cell.Shape.K - 1u,
                cell.Shape.BlockBytes,
                cell.Output.data() + cell.Shape.BlockBytes);
        }
        const Clock::time_point finish = Clock::now();
        output.InitFirstRepairNanoseconds =
            PositiveNanoseconds(first_repair - start);
        output.ElapsedNanoseconds = PositiveNanoseconds(finish - start);
        fresh.Reset(); // Destruction is intentionally after the stopped clock.
    }

    if (output.Result != 0 || output.ElapsedNanoseconds == 0u) return false;
    if (scope == Scope::FreshSystematic &&
        (output.ConstructorNanoseconds == 0u ||
         output.SystematicNanoseconds == 0u ||
         output.ConstructorNanoseconds > output.ElapsedNanoseconds ||
         output.SystematicNanoseconds > output.ElapsedNanoseconds ||
         output.ConstructorNanoseconds >
             output.ElapsedNanoseconds - output.SystematicNanoseconds))
    {
        return false;
    }
    if (scope == Scope::FreshRepair &&
        (output.InitFirstRepairNanoseconds == 0u ||
         output.InitFirstRepairNanoseconds > output.ElapsedNanoseconds))
    {
        return false;
    }

    if (IsFreshScope(scope))
    {
        if (IsV2Arm(state.ArmValue))
        {
            if (!ValidateV2Profile(cell, fresh_artifacts) ||
                fresh_artifacts.ProfileBytes !=
                    state.Oracle.SerializedProfileBytes ||
                wirehair::wh2_benchmark::Sha256Hex(
                    fresh_artifacts.Profile.data(),
                    fresh_artifacts.ProfileBytes) !=
                    state.Oracle.SerializedProfileSha256)
            {
                return false;
            }
        }
        else if (fresh_artifacts.ProfileBytes != 0u) {
            return false;
        }
    }

    output.CorrectnessSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Output.data(), cell.Output.size());
    if (output.CorrectnessSha256 != expected) return false;

    // This records calls eligible for the new borrowed-systematic policy.  It
    // is deliberately not an internal path counter or test-hook observation.
    output.BorrowedSystematicInvocations =
        state.ArmValue == Arm::C && systematic ?
        cell.Shape.K : 0u;
    output.DescriptorSha256 = state.Oracle.ArmDescriptorSha256;
    result_out = std::move(output);
    return true;
}

bool ParseUint64(const std::string& text, uint64_t& value_out)
{
    value_out = 0u;
    if (text.empty() || (text.size() > 1u && text.front() == '0')) {
        return false;
    }
    uint64_t value = 0u;
    for (char c : text)
    {
        if (c < '0' || c > '9') return false;
        const uint64_t digit = static_cast<uint64_t>(c - '0');
        if (value > (UINT64_MAX - digit) / 10u) return false;
        value = value * 10u + digit;
    }
    value_out = value;
    return true;
}

bool ParseUint32(const std::string& text, uint32_t& value_out)
{
    uint64_t value = 0u;
    if (!ParseUint64(text, value) || value > UINT32_MAX) return false;
    value_out = static_cast<uint32_t>(value);
    return true;
}

enum class CommandKind
{
    Invoke,
    Quit,
    Invalid
};

struct Command
{
    CommandKind Kind = CommandKind::Invalid;
    uint64_t Sequence = 0u;
    Arm ArmValue = Arm::Invalid;
    Scope ScopeValue = Scope::Invalid;
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
};

Command ParseCommand(const std::string& line)
{
    Command command;
    if (line == "quit")
    {
        command.Kind = CommandKind::Quit;
        return command;
    }
    if (line.empty() || line.front() == ' ' || line.back() == ' ' ||
        line.find("  ") != std::string::npos ||
        line.find_first_of("\t\r\n") != std::string::npos)
    {
        return command;
    }
    for (char raw : line) {
        const unsigned char c = static_cast<unsigned char>(raw);
        if (c < 0x20u || c > 0x7eu) return command;
    }
    std::istringstream input(line);
    std::string verb;
    std::string sequence;
    std::string arm;
    std::string scope;
    std::string K;
    std::string block_bytes;
    std::string extra;
    if (!(input >> verb >> sequence >> arm >> scope >> K >> block_bytes) ||
        input >> extra || verb != "invoke")
    {
        return command;
    }
    command.ArmValue = ParseArm(arm);
    command.ScopeValue = ParseScope(scope);
    if (!ParseUint64(sequence, command.Sequence) ||
        command.ArmValue == Arm::Invalid ||
        command.ScopeValue == Scope::Invalid ||
        !ParseUint32(K, command.K) ||
        !ParseUint32(block_bytes, command.BlockBytes))
    {
        return Command();
    }
    command.Kind = CommandKind::Invoke;
    return command;
}

bool ParseCpu(const char* text, int& cpu_out)
{
    cpu_out = -1;
    if (!text || !*text) return false;
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value < 0 || value > INT_MAX) {
        return false;
    }
    cpu_out = static_cast<int>(value);
    return true;
}

enum class ReadLineStatus
{
    Line,
    EndOfFile,
    Unterminated,
    TooLong,
    Error
};

ReadLineStatus ReadBoundedLine(std::istream& input, std::string& line_out)
{
    static const size_t kMaximumCommandBytes = 256u;
    line_out.clear();
    for (;;)
    {
        const int next = input.get();
        if (next == std::char_traits<char>::eof())
        {
            if (input.bad()) return ReadLineStatus::Error;
            return line_out.empty() ?
                ReadLineStatus::EndOfFile : ReadLineStatus::Unterminated;
        }
        if (next == '\n') return ReadLineStatus::Line;
        if (line_out.size() == kMaximumCommandBytes) {
            return ReadLineStatus::TooLong;
        }
        line_out.push_back(static_cast<char>(next));
    }
}

bool VerifySingletonCpu(int target_cpu, std::string& diagnostic)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "cpu_out_of_range";
        return false;
    }
    cpu_set_t observed;
    CPU_ZERO(&observed);
    if (sched_getaffinity(0, sizeof(observed), &observed) != 0)
    {
        diagnostic = "sched_getaffinity_failed";
        return false;
    }
    unsigned count = 0u;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        count += CPU_ISSET(static_cast<size_t>(cpu), &observed) ? 1u : 0u;
    }
    if (count != 1u ||
        !CPU_ISSET(static_cast<size_t>(target_cpu), &observed))
    {
        diagnostic = "affinity_not_singleton";
        return false;
    }
    if (sched_getcpu() != target_cpu)
    {
        diagnostic = "cpu_migration";
        return false;
    }
    diagnostic.clear();
    return true;
#else
    (void)target_cpu;
    diagnostic = "linux_affinity_required";
    return false;
#endif
}

bool PinSingletonCpu(int target_cpu, std::string& diagnostic)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "cpu_out_of_range";
        return false;
    }
    cpu_set_t selected;
    CPU_ZERO(&selected);
    CPU_SET(static_cast<size_t>(target_cpu), &selected);
    if (sched_setaffinity(0, sizeof(selected), &selected) != 0)
    {
        diagnostic = "sched_setaffinity_failed";
        return false;
    }
    return VerifySingletonCpu(target_cpu, diagnostic);
#else
    (void)target_cpu;
    diagnostic = "linux_affinity_required";
    return false;
#endif
}

void EmitNullableString(const std::string& value)
{
    if (value.empty()) std::cout << "null";
    else std::cout << '"' << value << '"';
}

void EmitConfig(const std::vector<Cell>& cells, int target_cpu)
{
#if defined(WIREHAIR_WH2_FACADE_ROLE_CURRENT)
    static const Arm supported[] = { Arm::C, Arm::E, Arm::D };
#else
    static const Arm supported[] = { Arm::P, Arm::L };
#endif
    std::cout << "{\"arm_descriptors\":[";
    for (size_t i = 0u; i < sizeof(supported) / sizeof(supported[0]); ++i)
    {
        if (i != 0u) std::cout << ',';
        const std::string descriptor = ArmDescriptorJson(supported[i]);
        std::cout << "{\"arm\":\"" << ArmName(supported[i])
                  << "\",\"descriptor\":" << descriptor
                  << ",\"descriptor_sha256\":\""
                  << wirehair::wh2_benchmark::Sha256Hex(descriptor)
                  << "\"}";
    }
    std::cout << "],\"campaign\":\"" << kCampaign << "\",\"cells\":[";
    for (size_t cell_index = 0u; cell_index < cells.size(); ++cell_index)
    {
        if (cell_index != 0u) std::cout << ',';
        const Cell& cell = cells[cell_index];
        std::cout << "{\"K\":" << cell.Shape.K
                  << ",\"block_bytes\":" << cell.Shape.BlockBytes
                  << ",\"invocations_by_replicate\":[";
        for (uint32_t replicate = 0u; replicate < kPanelReplicates;
             ++replicate)
        {
            if (replicate != 0u) std::cout << ',';
            std::cout << InvocationsForReplicate(cell.Shape.K, replicate);
        }
        std::cout << "],\"message_bytes\":" << cell.Source.size()
                  << ",\"oracles\":[";
        for (size_t arm_index = 0u; arm_index < cell.Arms.size(); ++arm_index)
        {
            if (arm_index != 0u) std::cout << ',';
            const ArmState& state = cell.Arms[arm_index];
            const OracleReceipt& oracle = state.Oracle;
            std::cout << "{\"arm\":\"" << ArmName(state.ArmValue)
                      << "\",\"arm_descriptor_sha256\":\""
                      << oracle.ArmDescriptorSha256
                      << "\",\"borrowed_systematic_ids\":"
                      << (state.ArmValue == Arm::C ? cell.Shape.K : 0u)
                      << ",\"equation_configuration_sha256\":\""
                      << oracle.EquationConfigurationSha256
                      << "\",\"first_repair_sha256\":\""
                      << oracle.FirstRepairSha256
                      << "\",\"high_id_sha256\":\""
                      << oracle.HighIdSha256
                      << "\",\"public_state_receipt_sha256\":\""
                      << oracle.PublicStateReceiptSha256
                      << "\",\"repair_sha256\":\""
                      << oracle.RepairSha256
                      << "\",\"roundtrip_sha256\":\""
                      << oracle.RoundtripSha256
                      << "\",\"roundtrip_verified\":"
                      << (oracle.RoundtripVerified ? "true" : "false")
                      << ",\"serialized_profile_bytes\":"
                      << oracle.SerializedProfileBytes
                      << ",\"serialized_profile_hex\":";
            EmitNullableString(oracle.SerializedProfileHex);
            std::cout << ",\"serialized_profile_sha256\":";
            EmitNullableString(oracle.SerializedProfileSha256);
            std::cout << ",\"systematic_sha256\":\""
                      << oracle.SystematicSha256 << "\"}";
        }
        std::cout << "],\"source_seed\":\"" << Hex64(cell.Seed)
                  << "\",\"source_sha256\":\"" << cell.SourceSha256
                  << "\"}";
    }
    std::cout << "],\"implementation_git_commit\":\""
              << WIREHAIR_WH2_IMPLEMENTATION_GIT_COMMIT
              << "\",\"internal_deadline_seconds\":"
              << kInternalDeadlineSeconds
              << ",\"invocation_budget\":" << kInvocationBudget
              << ",\"minimum_invocations\":" << kMinimumInvocations
              << ",\"panel_replicates\":" << kPanelReplicates
              << ",\"role\":\"" << CompiledRoleName()
              << "\",\"schema\":\"" << kConfigSchema
              << "\",\"supported_arms\":[";
    for (size_t i = 0u; i < sizeof(supported) / sizeof(supported[0]); ++i)
    {
        if (i != 0u) std::cout << ',';
        std::cout << '"' << ArmName(supported[i]) << '"';
    }
    std::cout << "],\"target_cpu\":" << target_cpu
              << ",\"worker_git_commit\":\""
              << WIREHAIR_WH2_HARNESS_GIT_COMMIT << "\"}\n" << std::flush;
}

void EmitInvocation(
    uint64_t sequence,
    Arm arm,
    Scope scope,
    const Cell& cell,
    const InvocationResult& result,
    int target_cpu)
{
    std::cout << "{\"K\":" << cell.Shape.K
              << ",\"arm\":\"" << ArmName(arm)
              << "\",\"block_bytes\":" << cell.Shape.BlockBytes
              << ",\"borrowed_systematic_invocations\":"
              << result.BorrowedSystematicInvocations
              << ",\"constructor_ns\":";
    if (scope == Scope::FreshSystematic) {
        std::cout << result.ConstructorNanoseconds;
    }
    else std::cout << "null";
    std::cout << ",\"correctness_sha256\":\""
              << result.CorrectnessSha256
              << "\",\"descriptor_sha256\":\""
              << result.DescriptorSha256
              << "\",\"elapsed_ns\":" << result.ElapsedNanoseconds
              << ",\"init_first_repair_ns\":";
    if (scope == Scope::FreshRepair) {
        std::cout << result.InitFirstRepairNanoseconds;
    }
    else std::cout << "null";
    std::cout << ",\"result\":" << result.Result
              << ",\"role\":\"" << CompiledRoleName()
              << "\",\"schema\":\"" << kInvocationSchema
              << "\",\"scope\":\"" << ScopeName(scope)
              << "\",\"seq\":" << sequence
              << ",\"systematic_ns\":";
    if (scope == Scope::FreshSystematic) {
        std::cout << result.SystematicNanoseconds;
    }
    else std::cout << "null";
    std::cout << ",\"target_cpu\":" << target_cpu << "}\n" << std::flush;
}

void EmitTerminal(
    uint64_t invocation_count,
    const char* status,
    int target_cpu)
{
    std::cout << "{\"invocation_count\":" << invocation_count
              << ",\"role\":\"" << CompiledRoleName()
              << "\",\"schema\":\"" << kTerminalSchema
              << "\",\"status\":\"" << status
              << "\",\"target_cpu\":" << target_cpu << "}\n"
              << std::flush;
}

bool TestJsonEmitters(const std::vector<Cell>& cells)
{
    if (cells.size() != 1u) return false;
    std::cout.flush();
    std::streambuf* const original = std::cout.rdbuf();

    std::ostringstream config;
    std::cout.rdbuf(config.rdbuf());
    EmitConfig(cells, -1);
    std::cout.rdbuf(original);
    const std::string config_text = config.str();
    const std::string config_prefix = "{\"arm_descriptors\":[";
    const std::string config_suffix =
        std::string("\"worker_git_commit\":\"") +
        WIREHAIR_WH2_HARNESS_GIT_COMMIT + "\"}\n";
    if (config_text.size() < config_prefix.size() + config_suffix.size() ||
        config_text.compare(0u, config_prefix.size(), config_prefix) != 0 ||
        config_text.compare(
            config_text.size() - config_suffix.size(),
            config_suffix.size(), config_suffix) != 0 ||
        std::count(config_text.begin(), config_text.end(), '\n') != 1 ||
        config_text.find(std::string("\"campaign\":\"") + kCampaign +
            "\"") == std::string::npos ||
        config_text.find(std::string("\"schema\":\"") + kConfigSchema +
            "\"") == std::string::npos ||
        config_text.find("\"internal_deadline_seconds\":840") ==
            std::string::npos)
    {
        return false;
    }

    InvocationResult sample;
    sample.Result = 0;
    sample.ElapsedNanoseconds = 5u;
    sample.ConstructorNanoseconds = 2u;
    sample.SystematicNanoseconds = 3u;
    sample.BorrowedSystematicInvocations = 8u;
    sample.CorrectnessSha256.assign(64u, 'a');
    sample.DescriptorSha256.assign(64u, 'b');
    std::ostringstream invocation;
    std::cout.rdbuf(invocation.rdbuf());
    EmitInvocation(
        7u, Arm::C, Scope::FreshSystematic, cells.front(), sample, -1);
    std::cout.rdbuf(original);
    const std::string expected_invocation =
        std::string("{\"K\":8,\"arm\":\"C\",\"block_bytes\":64,") +
        "\"borrowed_systematic_invocations\":8,\"constructor_ns\":2," +
        "\"correctness_sha256\":\"" + std::string(64u, 'a') +
        "\",\"descriptor_sha256\":\"" + std::string(64u, 'b') +
        "\",\"elapsed_ns\":5,\"init_first_repair_ns\":null," +
        "\"result\":0,\"role\":\"" + CompiledRoleName() +
        "\",\"schema\":\"" + kInvocationSchema +
        "\",\"scope\":\"fresh-systematic\",\"seq\":7," +
        "\"systematic_ns\":3,\"target_cpu\":-1}\n";
    if (invocation.str() != expected_invocation) return false;

    std::ostringstream terminal;
    std::cout.rdbuf(terminal.rdbuf());
    EmitTerminal(9u, "complete", -1);
    std::cout.rdbuf(original);
    const std::string expected_terminal =
        std::string("{\"invocation_count\":9,\"role\":\"") +
        CompiledRoleName() + "\",\"schema\":\"" + kTerminalSchema +
        "\",\"status\":\"complete\",\"target_cpu\":-1}\n";
    return terminal.str() == expected_terminal && std::cout.good();
}

bool TestParser()
{
    const Command good = ParseCommand(
        "invoke 17 C fresh-systematic 128 64");
    if (good.Kind != CommandKind::Invoke || good.Sequence != 17u ||
        good.ArmValue != Arm::C || good.ScopeValue != Scope::FreshSystematic ||
        good.K != 128u || good.BlockBytes != 64u)
    {
        return false;
    }
    return ParseCommand("quit").Kind == CommandKind::Quit &&
        ParseCommand("").Kind == CommandKind::Invalid &&
        ParseCommand(" invoke 0 C fresh-systematic 8 64").Kind ==
            CommandKind::Invalid &&
        ParseCommand("invoke 00 C fresh-systematic 8 64").Kind ==
            CommandKind::Invalid &&
        ParseCommand("invoke 0 C unknown 8 64").Kind ==
            CommandKind::Invalid &&
        ParseCommand("invoke 0 C fresh-systematic 8 64 extra").Kind ==
            CommandKind::Invalid &&
        ParseCommand("invoke 18446744073709551616 C fresh-systematic 8 64").Kind ==
            CommandKind::Invalid;
}

bool SelfTest()
{
    if (!TestParser() || wirehair_init() != Wirehair_Success) return false;
    for (const CellShape& shape : kCellShapes)
    {
        uint32_t total = 0u;
        for (uint32_t replicate = 0u; replicate < kPanelReplicates;
             ++replicate)
        {
            const uint32_t count =
                InvocationsForReplicate(shape.K, replicate);
            if (count < 2u) return false;
            total += count;
        }
        if (total != TotalInvocations(shape.K)) return false;
    }

    std::vector<Cell> cells;
    Cell cell;
    cell.Shape = kCellShapes[0];
    cell.Seed = SourceSeed(cell.Shape);
    if (!MakeSource(cell.Shape, cell.Seed, cell.Source)) return false;
    cell.Output.assign(cell.Source.size(), 0u);
    cell.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Source.data(), cell.Source.size());
    if (cell.SourceSha256 !=
        "da7d0d04aaeb241cc351650a839fd974f99be26766b7697de0af6b23a72c1276")
    {
        return false;
    }
    cells.push_back(std::move(cell));
    if (!BuildPrebuiltAndOracles(cells) || !CrossbindLocalV2Oracles(cells)) {
        return false;
    }
    for (ArmState& state : cells[0].Arms)
    {
        for (Scope scope : {
                 Scope::PrebuiltSystematic,
                 Scope::FreshSystematic,
                 Scope::FreshRepair,
                 Scope::PrebuiltRepair })
        {
            InvocationResult result;
            if (!RunInvocation(cells[0], state, scope, result) ||
                result.BorrowedSystematicInvocations !=
                    (state.ArmValue == Arm::C && IsSystematicScope(scope) ?
                        cells[0].Shape.K : 0u))
            {
                return false;
            }
        }
    }
    if (!VerifySourceIntegrity(cells)) return false;
    for (ArmState& state : cells[0].Arms) state.Prebuilt.Reset();
    cells[0].Source[0] ^= UINT8_C(0xff);
    if (VerifySourceIntegrity(cells)) return false;
    cells[0].Source[0] ^= UINT8_C(0xff);
    return VerifySourceIntegrity(cells) && TestJsonEmitters(cells) &&
        !ArmSupported(Arm::Invalid) &&
        ArmDescriptorJson(Arm::Invalid).empty();
}

bool Serve(const char* role, int target_cpu)
{
    if (!role || std::strcmp(role, CompiledRoleName()) != 0)
    {
        std::cerr << "compiled facade worker role mismatch\n";
        return false;
    }
    std::string cpu_diagnostic;
    if (!PinSingletonCpu(target_cpu, cpu_diagnostic))
    {
        std::cerr << "facade worker affinity failed: "
                  << cpu_diagnostic << '\n';
        return false;
    }
    if (wirehair_init() != Wirehair_Success)
    {
        std::cerr << "facade worker initialization failed\n";
        return false;
    }
    std::vector<Cell> cells;
    if (!BuildCellSources(cells) || !BuildPrebuiltAndOracles(cells) ||
        !CrossbindLocalV2Oracles(cells) ||
        !VerifySingletonCpu(target_cpu, cpu_diagnostic))
    {
        std::cerr << "facade worker startup closure failed\n";
        return false;
    }
    EmitConfig(cells, target_cpu);
    if (!std::cout.good()) return false;

    uint64_t invocation_count = 0u;
    bool have_sequence = false;
    uint64_t previous_sequence = 0u;
    std::string line;
    for (;;)
    {
        const ReadLineStatus read_status = ReadBoundedLine(std::cin, line);
        if (read_status == ReadLineStatus::EndOfFile)
        {
            EmitTerminal(invocation_count, "unexpected_eof", target_cpu);
            return false;
        }
        if (read_status != ReadLineStatus::Line)
        {
            EmitTerminal(invocation_count,
                read_status == ReadLineStatus::TooLong ?
                    "protocol_error" : "input_error",
                target_cpu);
            return false;
        }
        const Command command = ParseCommand(line);
        if (command.Kind == CommandKind::Quit)
        {
            if (!VerifySourceIntegrity(cells))
            {
                EmitTerminal(invocation_count, "source_error", target_cpu);
                return false;
            }
            if (!VerifySingletonCpu(target_cpu, cpu_diagnostic))
            {
                EmitTerminal(invocation_count, "affinity_error", target_cpu);
                return false;
            }
            EmitTerminal(invocation_count, "complete", target_cpu);
            return std::cout.good();
        }
        if (command.Kind != CommandKind::Invoke ||
            !ArmSupported(command.ArmValue) ||
            (have_sequence && command.Sequence <= previous_sequence))
        {
            EmitTerminal(invocation_count, "protocol_error", target_cpu);
            return false;
        }
        Cell* const cell = FindCell(
            cells, command.K, command.BlockBytes);
        ArmState* const state = cell ?
            FindArmState(*cell, command.ArmValue) : nullptr;
        if (!cell || !state ||
            !VerifySingletonCpu(target_cpu, cpu_diagnostic))
        {
            EmitTerminal(invocation_count,
                !cell || !state ? "protocol_error" : "affinity_error",
                target_cpu);
            return false;
        }
        InvocationResult result;
        if (!RunInvocation(
                *cell, *state, command.ScopeValue, result) ||
            !VerifySingletonCpu(target_cpu, cpu_diagnostic))
        {
            EmitTerminal(invocation_count,
                cpu_diagnostic.empty() ? "invocation_error" :
                    "affinity_error",
                target_cpu);
            return false;
        }
        EmitInvocation(
            command.Sequence, command.ArmValue, command.ScopeValue,
            *cell, result, target_cpu);
        if (!std::cout.good()) return false;
        previous_sequence = command.Sequence;
        have_sequence = true;
        ++invocation_count;
    }
}

} // namespace

int main(int argc, char** argv)
{
#if defined(SIGPIPE)
    if (std::signal(SIGPIPE, SIG_IGN) == SIG_ERR) return 1;
#endif
    if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0)
    {
        if (!SelfTest())
        {
            std::cerr << "WH2 V2 facade timing worker selftest failed\n";
            return 1;
        }
        std::cout << "WH2 V2 facade timing worker selftest passed\n"
                  << std::flush;
        return std::cout.good() ? 0 : 1;
    }
    int target_cpu = -1;
    if (argc != 6 || std::strcmp(argv[1], "--serve") != 0 ||
        std::strcmp(argv[2], "--role") != 0 ||
        std::strcmp(argv[4], "--cpu") != 0 ||
        !ParseCpu(argv[5], target_cpu))
    {
        std::cerr << "usage: " << argv[0]
                  << " --selftest | --serve --role current|parent --cpu N\n";
        return 2;
    }
    return Serve(argv[3], target_cpu) ? 0 : 1;
}
