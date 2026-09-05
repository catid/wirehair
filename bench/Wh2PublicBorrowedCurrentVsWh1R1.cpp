#include "Wh2FrozenTrace.h"
#include "Wh2NativePanel.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "Wh2RdpruTargetIdentityV2.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <climits>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#if defined(__linux__)
#include <dlfcn.h>
#include <fcntl.h>
#include <sched.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#error "public borrowed r1 requires an exact source commit receipt"
#endif

#ifndef WIREHAIR_WH2_CMAKE_CXX_COMPILER_ID
#error "public borrowed r1 requires the CMake compiler identifier"
#endif

#ifndef WIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION
#error "public borrowed r1 requires the CMake compiler version"
#endif

#ifndef WIREHAIR_WH2_CXX_COMPILER_PATH
#error "public borrowed r1 requires the exact compiler path"
#endif

#ifndef WIREHAIR_WH2_CXX_COMPILER_SHA256
#error "public borrowed r1 requires the exact compiler SHA-256"
#endif

namespace {

using wirehair_wh2_bench::NativePanelArm;
using wirehair_wh2_bench::NativePanelDisposition;
using wirehair_wh2_bench::NativePanelInvocation;
using wirehair_wh2_bench::NativePanelInvocationResult;
using wirehair_wh2_bench::NativePanelOrder;
using wirehair_wh2_bench::NativePanelResult;
using wirehair_wh2_bench::NativePanelSide;
using wirehair_wh2_bench::NativePanelStatus;

using Clock = std::chrono::steady_clock;

static const char kCampaign[] =
    "wh2-public-borrowed-current-vs-wh1-r1";
static const char kConfigSchema[] =
    "wirehair.wh2.public-borrowed-current-vs-wh1-r1.config.v1";
static const char kPanelSchema[] =
    "wirehair.wh2.public-borrowed-current-vs-wh1-r1.panel.v1";
static const char kTerminalSchema[] =
    "wirehair.wh2.public-borrowed-current-vs-wh1-r1.terminal.v1";
static const char kPanelKeySchema[] =
    "wirehair.wh2.public-borrowed-current-vs-wh1-r1.panel-key.v1";
static const char kArmSchema[] =
    "wirehair.wh2.public-borrowed-current-vs-wh1-r1.arm.v1";
static const char kWorkerStartedSchema[] =
    "wirehair.wh2.public-borrowed-current-vs-wh1-r1.worker-started.v1";
static const char kFixedEvidenceDir[] =
    "/var/tmp/wh2-public-borrowed-current-vs-wh1-r1";

static const uint32_t kPanelReplicates = 12u;
static const uint32_t kInvocationBudget = 65536u;
static const uint32_t kMinimumInvocations = 24u;
static const uint32_t kInternalDeadlineSeconds = 330u;
static const uint32_t kExpectedPanels = 1440u;
static const uint32_t kExpectedRecords = 1442u;
static const uint64_t kExpectedMeasuredInvocations = UINT64_C(852480);
static const uint64_t kExpectedWarmupInvocations = UINT64_C(2880);
static const uint64_t kExpectedMeasuredEncodeCalls = UINT64_C(185204736);
static const uint64_t kExpectedWarmupEncodeCalls = UINT64_C(41955840);
static const uint64_t kExpectedEncodeCalls = UINT64_C(227160576);
static const uint64_t kSourceSeedBase = UINT64_C(0xc199f24886210f53);
static const int kTargetCpu = 120;
static const int kSiblingCpu = 56;
static const uint8_t kScratchCanary = UINT8_C(0xa5);
static const uint32_t kLengthCanary = UINT32_C(0xa5a5a5a5);

struct CellShape
{
    uint32_t K;
    uint32_t BlockBytes;
};

static const CellShape kCellShapes[] = {
    { 8u, 64u }, { 8u, 1280u },
    { 128u, 64u }, { 128u, 1280u },
    { 512u, 64u }, { 512u, 1280u },
    { 8192u, 64u }, { 8192u, 1280u },
    { 64000u, 64u }, { 64000u, 1280u }
};

enum class Arm : uint32_t
{
    C = 0,
    L,
    Invalid
};

enum class Scope : uint32_t
{
    PrebuiltSystematic = 0,
    FreshSystematic,
    FreshRepair,
    PrebuiltRepair,
    Invalid
};

struct Comparison
{
    const char* Name;
    Arm Left;
    Arm Right;
};

static const Comparison kComparisons[] = {
    { "C/C", Arm::C, Arm::C },
    { "L/L", Arm::L, Arm::L },
    { "C/L", Arm::C, Arm::L }
};

static const Scope kScopes[] = {
    Scope::PrebuiltSystematic,
    Scope::FreshSystematic,
    Scope::FreshRepair,
    Scope::PrebuiltRepair
};

const char* ArmName(Arm arm)
{
    switch (arm)
    {
    case Arm::C: return "C";
    case Arm::L: return "L";
    default: return nullptr;
    }
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

bool IsSystematic(Scope scope)
{
    return scope == Scope::PrebuiltSystematic ||
        scope == Scope::FreshSystematic;
}

bool IsFresh(Scope scope)
{
    return scope == Scope::FreshSystematic ||
        scope == Scope::FreshRepair;
}

const char* SideName(NativePanelSide side)
{
    return side == NativePanelSide::Left ? "left" : "right";
}

const char* OrderName(NativePanelOrder order)
{
    return order == NativePanelOrder::ABBA ? "ABBA" : "BAAB";
}

bool CheckedProduct(uint64_t left, uint64_t right, size_t& output)
{
    output = 0u;
    if (left != 0u && right >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()) / left)
    {
        return false;
    }
    output = static_cast<size_t>(left * right);
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

uint64_t SourceSeed(const CellShape& shape, uint32_t tail_bytes)
{
    return kSourceSeedBase ^
        (static_cast<uint64_t>(shape.K) << 33) ^
        (static_cast<uint64_t>(shape.BlockBytes) << 1) ^ tail_bytes;
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
    const uint8_t* source = static_cast<const uint8_t*>(data);
    std::string output;
    try {
        output.resize(bytes * 2u);
    }
    catch (...) {
        return std::string();
    }
    for (size_t i = 0u; i < bytes; ++i)
    {
        output[i * 2u] = hex[source[i] >> 4];
        output[i * 2u + 1u] = hex[source[i] & 15u];
    }
    return output;
}

bool MakeSource(
    const CellShape& shape,
    uint32_t tail_bytes,
    std::vector<uint8_t>& output)
{
    if (shape.K < 2u || shape.K > 64000u || shape.BlockBytes == 0u ||
        tail_bytes == 0u || tail_bytes > shape.BlockBytes)
    {
        return false;
    }
    size_t prefix = 0u;
    if (!CheckedProduct(shape.K - 1u, shape.BlockBytes, prefix) ||
        prefix > std::numeric_limits<size_t>::max() - tail_bytes)
    {
        return false;
    }
    const size_t bytes = prefix + tail_bytes;
    try {
        output.assign(bytes, 0u);
    }
    catch (...) {
        return false;
    }
    uint64_t state = SourceSeed(shape, tail_bytes);
    size_t offset = 0u;
    while (offset < output.size())
    {
        const uint64_t word = SplitMix64(state);
        for (unsigned byte = 0u; byte < 8u && offset < output.size();
             ++byte, ++offset)
        {
            output[offset] = static_cast<uint8_t>(word >> (byte * 8u));
        }
    }
    return true;
}

uint32_t TotalInvocations(uint32_t K)
{
    const uint32_t quotient = kInvocationBudget / K;
    const uint32_t remainder = kInvocationBudget % K;
    return std::max(kMinimumInvocations,
        quotient + (remainder != 0u ? 1u : 0u));
}

uint32_t InvocationsPerSlot(uint32_t K, uint32_t replicate)
{
    const uint32_t total = TotalInvocations(K);
    return total / kPanelReplicates +
        (replicate < total % kPanelReplicates ? 1u : 0u);
}

uint64_t PositiveNanoseconds(Clock::duration duration)
{
    const std::chrono::nanoseconds ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
    return ns.count() > 0 ? static_cast<uint64_t>(ns.count()) : 0u;
}

bool ParseTargetCpu(const char* text, int& cpu_out)
{
    cpu_out = -1;
    if (!text || std::strcmp(text, "120") != 0) return false;
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value != kTargetCpu) {
        return false;
    }
    cpu_out = static_cast<int>(value);
    return true;
}

#if defined(__linux__)
bool IsExactSingletonCpuSet(const cpu_set_t& set, int target_cpu)
{
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE ||
        !CPU_ISSET(target_cpu, &set))
    {
        return false;
    }
    unsigned count = 0u;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        count += CPU_ISSET(cpu, &set) ? 1u : 0u;
    }
    return count == 1u;
}
#endif

bool PinAndVerifyTargetCpu(int target_cpu, std::string& diagnostic)
{
    diagnostic.clear();
#if defined(__linux__)
    if (target_cpu != kTargetCpu || target_cpu < 0 ||
        target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "target_cpu_outside_frozen_domain";
        return false;
    }
    cpu_set_t selected;
    CPU_ZERO(&selected);
    CPU_SET(target_cpu, &selected);
    if (sched_setaffinity(0, sizeof(selected), &selected) != 0)
    {
        diagnostic = "sched_setaffinity_failed";
        return false;
    }
    cpu_set_t observed;
    CPU_ZERO(&observed);
    if (sched_getaffinity(0, sizeof(observed), &observed) != 0 ||
        !IsExactSingletonCpuSet(observed, target_cpu))
    {
        diagnostic = "affinity_not_frozen_singleton";
        return false;
    }
    if (sched_getcpu() != target_cpu)
    {
        diagnostic = "cpu_migration";
        return false;
    }
    return true;
#else
    (void)target_cpu;
    diagnostic = "linux_affinity_required";
    return false;
#endif
}

bool VerifyTargetCpu(int target_cpu)
{
#if defined(__linux__)
    cpu_set_t observed;
    CPU_ZERO(&observed);
    return target_cpu == kTargetCpu &&
        sched_getaffinity(0, sizeof(observed), &observed) == 0 &&
        IsExactSingletonCpuSet(observed, target_cpu) &&
        sched_getcpu() == target_cpu;
#else
    (void)target_cpu;
    return false;
#endif
}

bool CanonicalExistingPath(const std::string& input, std::string& output)
{
    output.clear();
#if defined(__linux__)
    if (input.empty()) return false;
    std::array<char, PATH_MAX> resolved = {{0}};
    if (!realpath(input.c_str(), resolved.data())) return false;
    output.assign(resolved.data());
    return !output.empty() && output[0] == '/' &&
        output.find_first_of("\"\\\r\n") == std::string::npos;
#else
    (void)input;
    return false;
#endif
}

bool RuntimeWirehairMapsReceipt(
    std::string& receipt,
    std::string& canonical_path)
{
    receipt.clear();
    canonical_path.clear();
#if defined(__linux__)
    std::ifstream maps("/proc/self/maps");
    if (!maps) return false;
    std::string line;
    std::set<std::string> paths;
    while (std::getline(maps, line))
    {
        const size_t marker = line.find("libwirehair.so");
        if (marker == std::string::npos) continue;
        const size_t path_start = line.rfind(' ', marker);
        if (path_start == std::string::npos || path_start + 1u >= line.size()) {
            return false;
        }
        std::string path;
        if (!CanonicalExistingPath(line.substr(path_start + 1u), path)) {
            return false;
        }
        paths.insert(path);
        receipt += line;
        receipt.push_back('\n');
    }
    if (paths.size() != 1u || receipt.empty()) return false;
    canonical_path = *paths.begin();
    return true;
#else
    return false;
#endif
}

bool RuntimeTimedBindingPath(std::string& canonical_path)
{
    canonical_path.clear();
#if defined(__linux__)
    static const char* const kTimedSymbols[] = {
        "wirehair_v2_encoder_create_with_options",
        "wirehair_v2_encode",
        "wirehair_encoder_create_ex",
        "wirehair_encode"
    };
    std::set<std::string> paths;
    for (const char* name : kTimedSymbols)
    {
        (void)dlerror();
        void* symbol = dlsym(RTLD_DEFAULT, name);
        const char* error = dlerror();
        Dl_info info = {};
        std::string path;
        if (error || !symbol || dladdr(symbol, &info) == 0 ||
            !info.dli_fname || !CanonicalExistingPath(info.dli_fname, path))
        {
            return false;
        }
        paths.insert(path);
    }
    if (paths.size() != 1u) return false;
    canonical_path = *paths.begin();
    return true;
#else
    return false;
#endif
}

std::string RuntimeWirehairMapsText()
{
    std::string receipt;
    std::string path;
    return RuntimeWirehairMapsReceipt(receipt, path) ? receipt : std::string();
}

std::string ArmDescriptorJson(Arm arm)
{
    const char* api = nullptr;
    const char* codec = nullptr;
    const char* equation_profile = nullptr;
    const char* source_policy = nullptr;
    if (arm == Arm::C)
    {
        api = "wirehair_v2_encoder_create_with_options";
        codec = "wirehair2";
        equation_profile = "certified-2026-07";
        source_policy = "borrowed-immutable";
    }
    else if (arm == Arm::L)
    {
        api = "wirehair_encoder_create_ex";
        codec = "wirehair1";
        equation_profile = "legacy-current";
        source_policy = "borrowed";
    }
    else return std::string();
    std::ostringstream json;
    json << "{\"api\":\"" << api
         << "\",\"arm\":\"" << ArmName(arm)
         << "\",\"codec\":\"" << codec
         << "\",\"equation_profile\":\"" << equation_profile
         << "\",\"schema\":\"" << kArmSchema
         << "\",\"source_policy\":\"" << source_policy << "\"}";
    return json.str();
}

std::string PanelKeySha256(
    const CellShape& shape,
    const Comparison& comparison,
    Scope scope)
{
    if (!ArmName(comparison.Left) || !ArmName(comparison.Right) ||
        !ScopeName(scope))
    {
        return std::string();
    }
    std::ostringstream json;
    json << "{\"K\":" << shape.K
         << ",\"block_bytes\":" << shape.BlockBytes
         << ",\"campaign\":\"" << kCampaign
         << "\",\"comparison\":\"" << comparison.Name
         << "\",\"schema\":\"" << kPanelKeySchema
         << "\",\"scope\":\"" << ScopeName(scope) << "\"}";
    return wirehair::wh2_benchmark::Sha256Hex(json.str());
}

NativePanelOrder PanelOrder(
    const CellShape& shape,
    const Comparison& comparison,
    Scope scope,
    uint32_t replicate)
{
    const std::string digest = PanelKeySha256(shape, comparison, scope);
    if (digest.size() != 64u) {
        return static_cast<NativePanelOrder>(UINT32_MAX);
    }
    const char last = digest.back();
    const uint32_t nibble = last >= '0' && last <= '9' ?
        static_cast<uint32_t>(last - '0') :
        static_cast<uint32_t>(last - 'a' + 10);
    return ((replicate & 1u) == (nibble & 1u)) ?
        NativePanelOrder::ABBA : NativePanelOrder::BAAB;
}

class EncoderHandle
{
public:
    EncoderHandle() = default;
    ~EncoderHandle() { Reset(); }
    EncoderHandle(const EncoderHandle&) = delete;
    EncoderHandle& operator=(const EncoderHandle&) = delete;

    EncoderHandle(EncoderHandle&& other) noexcept
        : V2(other.V2), Legacy(other.Legacy), V2Arm(other.V2Arm)
    {
        other.V2 = nullptr;
        other.Legacy = nullptr;
        other.V2Arm = false;
    }

    EncoderHandle& operator=(EncoderHandle&& other) noexcept
    {
        if (this != &other)
        {
            Reset();
            V2 = other.V2;
            Legacy = other.Legacy;
            V2Arm = other.V2Arm;
            other.V2 = nullptr;
            other.Legacy = nullptr;
            other.V2Arm = false;
        }
        return *this;
    }

    void Adopt(WirehairV2Codec codec)
    {
        Reset();
        V2 = codec;
        V2Arm = true;
    }

    void Adopt(WirehairCodec codec)
    {
        Reset();
        Legacy = codec;
        V2Arm = false;
    }

    bool Valid() const { return V2Arm ? V2 != nullptr : Legacy != nullptr; }
    bool IsV2() const { return V2Arm; }
    WirehairV2Codec V2Codec() const { return V2; }
    WirehairCodec LegacyCodec() const { return Legacy; }

    void Reset()
    {
        if (V2) wirehair_v2_free(V2);
        if (Legacy) wirehair_free(Legacy);
        V2 = nullptr;
        Legacy = nullptr;
        V2Arm = false;
    }

private:
    WirehairV2Codec V2 = nullptr;
    WirehairCodec Legacy = nullptr;
    bool V2Arm = false;
};

struct ConstructionArtifacts
{
    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> Profile = {{0u}};
    uint32_t ProfileBytes = 0u;
};

int64_t CreateEncoder(
    Arm arm,
    const std::vector<uint8_t>& source,
    uint32_t block_bytes,
    EncoderHandle& encoder,
    ConstructionArtifacts& artifacts)
{
    if (encoder.Valid() || source.empty()) return WirehairV2_InvalidInput;
    if (arm == Arm::C)
    {
        WirehairV2EncoderOptions options = {};
        options.struct_bytes = sizeof(options);
        options.options_version = WIREHAIR_V2_ENCODER_OPTIONS_VERSION;
        options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
        options.reserved = 0u;
        WirehairV2Codec codec = nullptr;
        const WirehairV2Result result =
            wirehair_v2_encoder_create_with_options(
                source.data(), source.size(), block_bytes, &options,
                artifacts.Profile.data(),
                static_cast<uint32_t>(artifacts.Profile.size()),
                &artifacts.ProfileBytes, &codec);
        if (result != WirehairV2_Success || !codec ||
            artifacts.ProfileBytes != WIREHAIR_V2_PROFILE_SERIALIZED_BYTES)
        {
            wirehair_v2_free(codec);
            return result == WirehairV2_Success ? WirehairV2_Error : result;
        }
        encoder.Adopt(codec);
        return result;
    }
    if (arm == Arm::L)
    {
        WirehairCodec codec = nullptr;
        const WirehairResult result = wirehair_encoder_create_ex(
            nullptr, source.data(), source.size(), block_bytes, &codec);
        if (result != Wirehair_Success || !codec)
        {
            wirehair_free(codec);
            return result == Wirehair_Success ? Wirehair_Error : result;
        }
        encoder.Adopt(codec);
        return result;
    }
    return WirehairV2_InvalidInput;
}

bool ValidateProfile(
    const std::vector<uint8_t>& source,
    uint32_t block_bytes,
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
    return profile.struct_bytes == sizeof(profile) &&
        profile.profile_version == WIREHAIR_V2_PROFILE_VERSION &&
        profile.profile_id == WIREHAIR_V2_PROFILE_CERTIFIED_2026_07 &&
        profile.message_bytes == source.size() &&
        profile.block_bytes == block_bytes &&
        profile.reserved[0] == 0u && profile.reserved[1] == 0u &&
        profile.reserved[2] == 0u;
}

uint32_t ExpectedPacketBytes(
    uint32_t id,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t final_bytes)
{
    return id < K && id == K - 1u ? final_bytes : block_bytes;
}

int64_t EncodeOne(
    const EncoderHandle& encoder,
    uint32_t id,
    void* output,
    uint32_t capacity,
    uint32_t& bytes_out)
{
    if (!encoder.Valid()) return WirehairV2_InvalidInput;
    if (encoder.IsV2()) {
        return wirehair_v2_encode(
            encoder.V2Codec(), id, output, capacity, &bytes_out);
    }
    return wirehair_encode(
        encoder.LegacyCodec(), id, output, capacity, &bytes_out);
}

} // namespace

namespace {

bool TestRosterArithmetic()
{
    uint32_t panel_count = 0u;
    uint64_t measured_invocations = 0u;
    uint64_t warmup_invocations = 0u;
    uint64_t measured_encode_calls = 0u;
    uint64_t warmup_encode_calls = 0u;
    std::set<std::string> panel_keys;
    for (const CellShape& shape : kCellShapes)
    {
        uint32_t total = 0u;
        for (uint32_t replicate = 0u;
             replicate < kPanelReplicates; ++replicate)
        {
            const uint32_t n = InvocationsPerSlot(shape.K, replicate);
            if (n < 2u) return false;
            total += n;
            for (Scope scope : kScopes)
            {
                for (const Comparison& comparison : kComparisons)
                {
                    const NativePanelOrder order = PanelOrder(
                        shape, comparison, scope, replicate);
                    if (order != NativePanelOrder::ABBA &&
                        order != NativePanelOrder::BAAB)
                    {
                        return false;
                    }
                    ++panel_count;
                    measured_invocations += UINT64_C(4) * n;
                    warmup_invocations += 2u;
                    measured_encode_calls +=
                        UINT64_C(4) * n * shape.K;
                    warmup_encode_calls += UINT64_C(2) * shape.K;
                }
            }
        }
        if (total != TotalInvocations(shape.K)) return false;
        for (Scope scope : kScopes)
        {
            for (const Comparison& comparison : kComparisons)
            {
                const std::string key = PanelKeySha256(
                    shape, comparison, scope);
                if (key.size() != 64u || !panel_keys.insert(key).second) {
                    return false;
                }
                uint32_t abba = 0u;
                uint32_t baab = 0u;
                for (uint32_t replicate = 0u;
                     replicate < kPanelReplicates; ++replicate)
                {
                    if (PanelOrder(
                            shape, comparison, scope, replicate) ==
                            NativePanelOrder::ABBA)
                    {
                        ++abba;
                    }
                    else ++baab;
                }
                if (abba != 6u || baab != 6u) return false;
            }
        }
    }
    return panel_keys.size() == 120u &&
        panel_count == kExpectedPanels &&
        measured_invocations == kExpectedMeasuredInvocations &&
        warmup_invocations == kExpectedWarmupInvocations &&
        measured_encode_calls == kExpectedMeasuredEncodeCalls &&
        warmup_encode_calls == kExpectedWarmupEncodeCalls &&
        measured_encode_calls + warmup_encode_calls ==
            kExpectedEncodeCalls &&
        kExpectedRecords == kExpectedPanels + 2u;
}

bool TestCpuParserAndMatcher()
{
    int parsed = -1;
    if (!ParseTargetCpu("120", parsed) || parsed != kTargetCpu) return false;
    for (const char* bad : { "", "119", "121", "120x", "0120" })
    {
        parsed = -1;
        if (ParseTargetCpu(bad, parsed) || parsed != -1) return false;
    }
    parsed = -1;
    if (ParseTargetCpu(nullptr, parsed) || parsed != -1) return false;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (IsExactSingletonCpuSet(set, 0)) return false;
    CPU_SET(0, &set);
    if (!IsExactSingletonCpuSet(set, 0)) return false;
    CPU_SET(1, &set);
    if (IsExactSingletonCpuSet(set, 0)) return false;
#endif
    return true;
}

bool IsLowerSha256(const char* text)
{
    if (!text || std::strlen(text) != 64u) return false;
    for (size_t i = 0u; i < 64u; ++i)
    {
        const char c = text[i];
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) {
            return false;
        }
    }
    return true;
}

std::string WorkerStartedJson(
    const char* claim_sha256,
    const char* worker_sha256,
    const char* native_config_identity_sha256)
{
    if (!IsLowerSha256(claim_sha256) || !IsLowerSha256(worker_sha256) ||
        !IsLowerSha256(native_config_identity_sha256))
    {
        return std::string();
    }
    std::ostringstream json;
    json << "{\"campaign\":\"" << kCampaign
         << "\",\"claim_sha256\":\"" << claim_sha256
         << "\",\"native_config_identity_sha256\":\""
         << native_config_identity_sha256
         << "\",\"schema\":\"" << kWorkerStartedSchema
         << "\",\"source_commit\":\""
         << WIREHAIR_WH2_SOURCE_GIT_COMMIT
         << "\",\"status\":\"started\",\"worker_sha256\":\""
         << worker_sha256 << "\"}\n";
    return json.str();
}

#if defined(__linux__)
bool ReadFdBounded(
    int fd,
    size_t maximum,
    std::string& output,
    std::string& diagnostic)
{
    output.clear();
    std::array<char, 65536u> buffer;
    for (;;)
    {
        const ssize_t bytes = ::read(fd, buffer.data(), buffer.size());
        if (bytes < 0)
        {
            if (errno == EINTR) continue;
            diagnostic = std::string("read failed: ") + std::strerror(errno);
            return false;
        }
        if (bytes == 0) return true;
        const size_t count = static_cast<size_t>(bytes);
        if (count > maximum - std::min(maximum, output.size()))
        {
            diagnostic = "file exceeds its authorization byte cap";
            return false;
        }
        output.append(buffer.data(), count);
    }
}

bool WriteFdAll(int fd, const std::string& data, std::string& diagnostic)
{
    size_t offset = 0u;
    while (offset < data.size())
    {
        const ssize_t bytes = ::write(
            fd, data.data() + offset, data.size() - offset);
        if (bytes < 0)
        {
            if (errno == EINTR) continue;
            diagnostic = std::string("write failed: ") + std::strerror(errno);
            return false;
        }
        if (bytes == 0)
        {
            diagnostic = "authorization marker write made no progress";
            return false;
        }
        offset += static_cast<size_t>(bytes);
    }
    return true;
}

bool VerifyRegularAuthorizationFile(
    int fd,
    mode_t expected_mode,
    size_t maximum,
    std::string& data,
    std::string& diagnostic)
{
    struct stat info;
    if (::fstat(fd, &info) != 0)
    {
        diagnostic = std::string("authorization fstat failed: ") +
            std::strerror(errno);
        return false;
    }
    if (!S_ISREG(info.st_mode) ||
        (expected_mode != 0 && (info.st_mode & 07777) != expected_mode) ||
        info.st_nlink != 1 ||
        info.st_size < 0 || static_cast<uint64_t>(info.st_size) > maximum)
    {
        diagnostic = "authorization file metadata differs";
        return false;
    }
    return ReadFdBounded(fd, maximum, data, diagnostic);
}

bool AuthorizeRunInDirectory(
    const char* evidence_dir,
    const char* claim_sha256,
    const char* worker_sha256,
    const char* native_config_identity_sha256,
    std::string& diagnostic)
{
    if (!evidence_dir || !IsLowerSha256(claim_sha256) ||
        !IsLowerSha256(worker_sha256) ||
        !IsLowerSha256(native_config_identity_sha256))
    {
        diagnostic = "launch authorization arguments differ";
        return false;
    }
    const int directory_fd = ::open(
        evidence_dir, O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW);
    if (directory_fd < 0)
    {
        diagnostic = std::string("cannot open evidence directory: ") +
            std::strerror(errno);
        return false;
    }
    bool authorized = false;
    int claim_fd = -1;
    int executable_fd = -1;
    int marker_fd = -1;
    do
    {
        struct stat directory_info;
        if (::fstat(directory_fd, &directory_info) != 0 ||
            !S_ISDIR(directory_info.st_mode) ||
            (directory_info.st_mode & 07777) != 0700)
        {
            diagnostic = "evidence directory metadata differs";
            break;
        }
        claim_fd = ::openat(
            directory_fd, "CLAIM", O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
        if (claim_fd < 0)
        {
            diagnostic = std::string("cannot open CLAIM: ") +
                std::strerror(errno);
            break;
        }
        std::string claim;
        if (!VerifyRegularAuthorizationFile(
                claim_fd, 0400, 1024u * 1024u, claim, diagnostic) ||
            wirehair::wh2_benchmark::Sha256Hex(claim) != claim_sha256)
        {
            if (diagnostic.empty()) diagnostic = "CLAIM SHA-256 differs";
            break;
        }
        executable_fd = ::open("/proc/self/exe", O_RDONLY | O_CLOEXEC);
        if (executable_fd < 0)
        {
            diagnostic = std::string("cannot open /proc/self/exe: ") +
                std::strerror(errno);
            break;
        }
        std::string executable;
        if (!VerifyRegularAuthorizationFile(
                executable_fd, 0, 64u * 1024u * 1024u,
                executable, diagnostic) ||
            wirehair::wh2_benchmark::Sha256Hex(executable) != worker_sha256)
        {
            if (diagnostic.empty()) diagnostic = "worker SHA-256 differs";
            break;
        }
        const std::string marker =
            WorkerStartedJson(
                claim_sha256, worker_sha256,
                native_config_identity_sha256);
        if (marker.empty())
        {
            diagnostic = "cannot construct worker authorization marker";
            break;
        }
        marker_fd = ::openat(
            directory_fd, "WORKER_STARTED",
            O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC | O_NOFOLLOW, 0400);
        if (marker_fd < 0)
        {
            diagnostic = std::string("cannot exclusively create WORKER_STARTED: ") +
                std::strerror(errno);
            break;
        }
        if (::fchmod(marker_fd, 0400) != 0)
        {
            diagnostic = std::string("WORKER_STARTED chmod failed: ") +
                std::strerror(errno);
            break;
        }
        struct stat marker_info;
        if (::fstat(marker_fd, &marker_info) != 0 ||
            !S_ISREG(marker_info.st_mode) ||
            (marker_info.st_mode & 07777) != 0400 ||
            marker_info.st_nlink != 1 ||
            !WriteFdAll(marker_fd, marker, diagnostic) ||
            ::fsync(marker_fd) != 0 || ::fsync(directory_fd) != 0)
        {
            if (diagnostic.empty()) {
                diagnostic = std::string("WORKER_STARTED durability failed: ") +
                    std::strerror(errno);
            }
            break;
        }
        authorized = true;
    } while (false);
    if (marker_fd >= 0) ::close(marker_fd);
    if (executable_fd >= 0) ::close(executable_fd);
    if (claim_fd >= 0) ::close(claim_fd);
    ::close(directory_fd);
    return authorized;
}

bool AuthorizeRun(
    const char* evidence_dir,
    const char* claim_sha256,
    const char* worker_sha256,
    const char* native_config_identity_sha256,
    std::string& diagnostic)
{
    if (!evidence_dir || std::strcmp(evidence_dir, kFixedEvidenceDir) != 0)
    {
        diagnostic = "fixed evidence directory differs";
        return false;
    }
    return AuthorizeRunInDirectory(
        evidence_dir, claim_sha256, worker_sha256,
        native_config_identity_sha256, diagnostic);
}

bool TestLaunchAuthorization()
{
    char directory_template[] = "/tmp/wh2-r1-auth-selftest.XXXXXX";
    char* const directory = ::mkdtemp(directory_template);
    if (!directory) return false;
    const std::string claim_path = std::string(directory) + "/CLAIM";
    const std::string marker_path =
        std::string(directory) + "/WORKER_STARTED";
    const std::string claim = "{\"selftest\":true}\n";
    const std::string claim_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(claim);
    const std::string config_identity_sha256(64u, 'b');
    bool success = false;
    int claim_fd = -1;
    int executable_fd = -1;
    int marker_fd = -1;
    do
    {
        claim_fd = ::open(
            claim_path.c_str(),
            O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC | O_NOFOLLOW, 0400);
        std::string diagnostic;
        if (claim_fd < 0 || ::fchmod(claim_fd, 0400) != 0 ||
            !WriteFdAll(claim_fd, claim, diagnostic) || ::fsync(claim_fd) != 0)
        {
            break;
        }
        ::close(claim_fd);
        claim_fd = -1;
        executable_fd = ::open("/proc/self/exe", O_RDONLY | O_CLOEXEC);
        std::string executable;
        if (executable_fd < 0 || !VerifyRegularAuthorizationFile(
                executable_fd, 0, 64u * 1024u * 1024u,
                executable, diagnostic))
        {
            break;
        }
        ::close(executable_fd);
        executable_fd = -1;
        const std::string worker_sha256 =
            wirehair::wh2_benchmark::Sha256Hex(executable);
        if (!AuthorizeRunInDirectory(
                directory, claim_sha256.c_str(), worker_sha256.c_str(),
                config_identity_sha256.c_str(), diagnostic))
        {
            break;
        }
        diagnostic.clear();
        if (AuthorizeRunInDirectory(
                directory, claim_sha256.c_str(), worker_sha256.c_str(),
                config_identity_sha256.c_str(), diagnostic))
        {
            break;
        }
        marker_fd = ::open(
            marker_path.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
        std::string marker;
        if (marker_fd < 0 || !VerifyRegularAuthorizationFile(
                marker_fd, 0400, 1024u * 1024u, marker, diagnostic) ||
            marker != WorkerStartedJson(
                claim_sha256.c_str(), worker_sha256.c_str(),
                config_identity_sha256.c_str()))
        {
            break;
        }
        success = true;
    } while (false);
    if (marker_fd >= 0) ::close(marker_fd);
    if (executable_fd >= 0) ::close(executable_fd);
    if (claim_fd >= 0) ::close(claim_fd);
    const bool cleanup =
        (::unlink(marker_path.c_str()) == 0 || errno == ENOENT) &&
        (::unlink(claim_path.c_str()) == 0 || errno == ENOENT) &&
        ::rmdir(directory) == 0;
    return success && cleanup;
}
#else
bool AuthorizeRun(
    const char*, const char*, const char*, const char*,
    std::string& diagnostic)
{
    diagnostic = "real public borrowed r1 runs require Linux";
    return false;
}

bool TestLaunchAuthorization()
{
    return true;
}
#endif

bool SelfTest();
bool EmitNativeConfigSelfTest();
bool RunScreen(
    int target_cpu,
    const std::string& expected_native_config_identity_sha256);

} // namespace

int main(int argc, char** argv)
{
#if defined(SIGPIPE)
    if (std::signal(SIGPIPE, SIG_IGN) == SIG_ERR)
    {
        std::cerr << "cannot ignore SIGPIPE for transactional output\n";
        return 1;
    }
#endif
    if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0)
    {
        if (!SelfTest())
        {
            std::cerr << "WH2 public borrowed current-vs-WH1 r1 "
                         "selftest failed\n";
            return 1;
        }
        std::cout << "WH2 public borrowed current-vs-WH1 r1 "
                     "selftest passed\n" << std::flush;
        return std::cout.good() ? 0 : 1;
    }
    if (argc == 2 &&
        std::strcmp(argv[1], "--emit-config-selftest") == 0)
    {
        return EmitNativeConfigSelfTest() ? 0 : 1;
    }
    int target_cpu = -1;
    if (argc != 12 || std::strcmp(argv[1], "--run") != 0 ||
        std::strcmp(argv[2], "--cpu") != 0 ||
        !ParseTargetCpu(argv[3], target_cpu) ||
        std::strcmp(argv[4], "--evidence-dir") != 0 ||
        std::strcmp(argv[6], "--claim-sha256") != 0 ||
        std::strcmp(argv[8], "--worker-sha256") != 0 ||
        std::strcmp(argv[10], "--config-identity-sha256") != 0 ||
        !IsLowerSha256(argv[11]))
    {
        std::cerr << "usage: " << argv[0]
                  << " --selftest | --emit-config-selftest |"
                     " --run --cpu 120 --evidence-dir PATH"
                     " --claim-sha256 HEX --worker-sha256 HEX"
                     " --config-identity-sha256 HEX\n";
        return 2;
    }
    std::string authorization_diagnostic;
    if (!AuthorizeRun(
            argv[5], argv[7], argv[9], argv[11],
            authorization_diagnostic))
    {
        std::cerr << "public borrowed r1 launch authorization failed: "
                  << authorization_diagnostic << '\n';
        return 1;
    }
    return RunScreen(target_cpu, argv[11]) ? 0 : 1;
}

namespace {

struct OracleReceipt
{
    Arm ArmValue = Arm::Invalid;
    bool AttachedOverlapNoWriteVerified = false;
    uint32_t BorrowedEligibleSystematicIds = 0u;
    std::string DescriptorSha256;
    std::string EquationConfigurationSha256;
    std::string FirstRepairSha256;
    std::string FullRepairSha256;
    std::string HighIdSha256;
    uint32_t ProfileBytes = 0u;
    std::string ProfileHex;
    std::string ProfileSha256;
    uint32_t RoundtripFirstId = 0u;
    uint32_t RoundtripPacketCount = 0u;
    bool RoundtripRepairOnlyVerified = false;
    std::string RoundtripSha256;
    bool RoundtripVerified = false;
    std::string SystematicSha256;
};

struct PartialArmOracle
{
    Arm ArmValue = Arm::Invalid;
    bool AttachedOverlapNoWriteVerified = false;
    std::string FirstRepairSha256;
    std::string HighIdSha256;
    std::string ProfileSha256;
    uint32_t RoundtripFirstId = 0u;
    uint32_t RoundtripPacketCount = 0u;
    bool RoundtripRepairOnlyVerified = false;
    std::string RoundtripSha256;
    bool SourceImmutableVerified = false;
    std::string SystematicSha256;
    bool SystematicVerified = false;
};

struct PartialFinalOracle
{
    uint32_t TailBytes = 0u;
    std::string SourceSha256;
    std::array<PartialArmOracle, 2> Arms;
};

struct ScreenCell
{
    CellShape Shape = {};
    uint64_t SourceSeed = 0u;
    std::vector<uint8_t> Source;
    std::string SourceSha256;
    std::array<OracleReceipt, 2> Oracles;
    std::array<PartialFinalOracle, 2> PartialFinalOracles;
};

const OracleReceipt* FindOracle(const ScreenCell& cell, Arm arm);
bool PrepareScratch(
    uint32_t K,
    uint32_t block_bytes,
    std::vector<uint8_t>& scratch,
    std::vector<uint32_t>& lengths);
int64_t EncodeBatch(
    const EncoderHandle& encoder,
    uint32_t first_id,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t final_bytes,
    std::vector<uint8_t>& scratch,
    std::vector<uint32_t>& lengths);
bool VerifyScratch(
    const std::vector<uint8_t>& source,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t final_bytes,
    bool systematic,
    const std::vector<uint8_t>& scratch,
    const std::vector<uint32_t>& lengths,
    std::string& sha256_out);

class PanelLane
{
public:
    PanelLane(
        const ScreenCell& cell,
        Arm arm,
        Scope scope,
        const char* lane_name)
        : Cell(cell)
        , ArmValue(arm)
        , ScopeValue(scope)
        , Oracle(FindOracle(cell, arm))
    {
        std::ostringstream identity;
        identity << "arm=" << (ArmName(arm) ? ArmName(arm) : "invalid")
                 << ";K=" << cell.Shape.K
                 << ";B=" << cell.Shape.BlockBytes
                 << ";scope=" <<
                    (ScopeName(scope) ? ScopeName(scope) : "invalid")
                 << ";lane=" << (lane_name ? lane_name : "invalid")
                 << ";descriptor=" <<
                    (Oracle ? Oracle->DescriptorSha256 : std::string());
        IdentityValue = identity.str();
    }

    bool Initialize()
    {
        if (!Oracle || !ArmName(ArmValue) || !ScopeName(ScopeValue) ||
            !PrepareScratch(
                Cell.Shape.K, Cell.Shape.BlockBytes, Scratch, Lengths))
        {
            return false;
        }
        if (!IsFresh(ScopeValue))
        {
            if (CreateEncoder(
                    ArmValue, Cell.Source, Cell.Shape.BlockBytes,
                    Prebuilt, PrebuiltArtifacts) != 0 ||
                !Prebuilt.Valid())
            {
                return false;
            }
            if (ArmValue == Arm::C &&
                (!ValidateProfile(
                    Cell.Source, Cell.Shape.BlockBytes,
                    PrebuiltArtifacts) ||
                 wirehair::wh2_benchmark::Sha256Hex(
                    PrebuiltArtifacts.Profile.data(),
                    PrebuiltArtifacts.ProfileBytes) !=
                    Oracle->ProfileSha256))
            {
                return false;
            }
        }
        return true;
    }

    bool Prepare()
    {
        if (Scratch.empty() || Lengths.size() != Cell.Shape.K) return false;
        std::fill(Scratch.begin(), Scratch.end(), kScratchCanary);
        std::fill(Lengths.begin(), Lengths.end(), kLengthCanary);
        LastInvocationComplete = false;
        LastOutputSha256.clear();
        return true;
    }

    const std::string& Identity() const { return IdentityValue; }

    NativePanelInvocationResult Invoke()
    {
        EncoderHandle fresh;
        ConstructionArtifacts artifacts;
        EncoderHandle* encoder = &Prebuilt;
        int64_t result = 0;
        const uint32_t first_id = IsSystematic(ScopeValue) ?
            0u : Cell.Shape.K;
        const Clock::time_point start = Clock::now();
        if (IsFresh(ScopeValue))
        {
            result = CreateEncoder(
                ArmValue, Cell.Source, Cell.Shape.BlockBytes,
                fresh, artifacts);
            encoder = &fresh;
        }
        if (result == 0)
        {
            result = EncodeBatch(
                *encoder, first_id, Cell.Shape.K, Cell.Shape.BlockBytes,
                Cell.Shape.BlockBytes, Scratch, Lengths);
        }
        const Clock::time_point finish = Clock::now();
        const uint64_t elapsed = PositiveNanoseconds(finish - start);
        fresh.Reset();

        if (result != 0 || elapsed == 0u)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                result == 0 ?
                    static_cast<int64_t>(WirehairV2_Error) : result,
                true, Cell.Shape.K, elapsed);
        }
        if (IsFresh(ScopeValue))
        {
            if (ArmValue == Arm::C &&
                (!ValidateProfile(
                    Cell.Source, Cell.Shape.BlockBytes, artifacts) ||
                 wirehair::wh2_benchmark::Sha256Hex(
                    artifacts.Profile.data(), artifacts.ProfileBytes) !=
                    Oracle->ProfileSha256))
            {
                return NativePanelInvocationResult(
                    NativePanelDisposition::Fatal,
                    WirehairV2_Error, true, Cell.Shape.K, elapsed);
            }
            if (ArmValue == Arm::L && artifacts.ProfileBytes != 0u)
            {
                return NativePanelInvocationResult(
                    NativePanelDisposition::Fatal,
                    WirehairV2_Error, true, Cell.Shape.K, elapsed);
            }
        }
        LastInvocationComplete = true;
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            0, true, Cell.Shape.K, elapsed);
    }

    bool VerifyAfterPanel()
    {
        if (!LastInvocationComplete || !Oracle) return false;
        const bool systematic = IsSystematic(ScopeValue);
        if (!VerifyScratch(
                Cell.Source, Cell.Shape.K, Cell.Shape.BlockBytes,
                Cell.Shape.BlockBytes, systematic, Scratch, Lengths,
                LastOutputSha256))
        {
            return false;
        }
        const std::string& expected = systematic ?
            Oracle->SystematicSha256 : Oracle->FullRepairSha256;
        return LastOutputSha256 == expected;
    }

    const std::string& OutputSha256() const { return LastOutputSha256; }

private:
    const ScreenCell& Cell;
    Arm ArmValue;
    Scope ScopeValue;
    const OracleReceipt* Oracle;
    std::string IdentityValue;
    EncoderHandle Prebuilt;
    ConstructionArtifacts PrebuiltArtifacts;
    std::vector<uint8_t> Scratch;
    std::vector<uint32_t> Lengths;
    bool LastInvocationComplete = false;
    std::string LastOutputSha256;
};

class PublicInvocation : public NativePanelInvocation
{
public:
    explicit PublicInvocation(PanelLane& lane)
        : Lane(lane)
    {
    }

    std::string Identity() const override { return Lane.Identity(); }

    NativePanelInvocationResult Invoke() override
    {
        if (Invoked)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                WirehairV2_InvalidInput, false, 0u, 0u);
        }
        Invoked = true;
        return Lane.Invoke();
    }

private:
    PanelLane& Lane;
    bool Invoked = false;
};

NativePanelArm MakePanelArm(PanelLane& lane)
{
    return NativePanelArm(
        lane.Identity(),
        [&lane]() -> std::unique_ptr<NativePanelInvocation> {
            if (!lane.Prepare()) {
                return std::unique_ptr<NativePanelInvocation>();
            }
            return std::unique_ptr<NativePanelInvocation>(
                new PublicInvocation(lane));
        });
}

bool VerifyPanelResult(
    const NativePanelResult& panel,
    uint32_t K,
    uint32_t invocations_per_slot,
    int target_cpu)
{
    if (panel.Status != NativePanelStatus::Complete ||
        panel.IsFatal() || !panel.PanelComparable ||
        !panel.HasLeftPreflight || !panel.HasRightPreflight ||
        panel.TargetCpu != target_cpu ||
        panel.InvocationsPerSlot != invocations_per_slot ||
        panel.LeftPreflight.Disposition != NativePanelDisposition::Success ||
        panel.RightPreflight.Disposition != NativePanelDisposition::Success ||
        panel.LeftPreflight.OutcomeCode != 0 ||
        panel.RightPreflight.OutcomeCode != 0 ||
        !panel.LeftPreflight.HasDecodedExtra ||
        !panel.RightPreflight.HasDecodedExtra ||
        panel.LeftPreflight.DecodedExtra != K ||
        panel.RightPreflight.DecodedExtra != K)
    {
        return false;
    }
    for (size_t i = 0u; i < panel.Slots.size(); ++i)
    {
        const uint32_t expected_count = i < 4u ?
            invocations_per_slot / 2u + invocations_per_slot % 2u :
            invocations_per_slot / 2u;
        const auto& slot = panel.Slots[i];
        if (expected_count == 0u || !slot.HasElapsedNanoseconds ||
            slot.ElapsedNanoseconds == 0u ||
            slot.Invocation.Disposition != NativePanelDisposition::Success ||
            slot.Invocation.OutcomeCode != 0 ||
            !slot.Invocation.HasDecodedExtra ||
            slot.Invocation.DecodedExtra != K)
        {
            return false;
        }
    }
    return true;
}

void EmitNullableString(std::ostream& output, const std::string& value)
{
    if (value.empty()) output << "null";
    else output << '"' << value << '"';
}

std::string TargetIdentityJson(
    const wirehair_wh2_bench::TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
    std::string canonical;
    if (!wirehair_wh2_bench::SerializeTargetIdentityV2(
            receipt, canonical, diagnostic) || canonical.empty() ||
        canonical.size() > 4096u ||
        wirehair::wh2_benchmark::Sha256Hex(canonical) !=
            receipt.CanonicalSha256 ||
        receipt.Before.Affinity.size() != 1u ||
        receipt.After.Affinity.size() != 1u ||
        receipt.Before.Affinity[0] != receipt.RequestedCpu ||
        receipt.After.Affinity[0] != receipt.RequestedCpu)
    {
        if (diagnostic.empty()) diagnostic = "target identity receipt differs";
        return std::string();
    }
    const wirehair_wh2_bench::CpuidSnapshotV2& raw = receipt.Raw;
    const wirehair_wh2_bench::DerivedTargetIdentityV2& derived =
        receipt.Derived;
    std::ostringstream json;
    json << "{\"after_cpu\":" << receipt.After.Cpu
         << ",\"before_cpu\":" << receipt.Before.Cpu
         << ",\"canonical_bytes\":" << canonical.size()
         << ",\"canonical_hex\":\""
         << HexBytes(canonical.data(), canonical.size())
         << "\",\"canonical_sha256\":\"" << receipt.CanonicalSha256
         << "\",\"capabilities\":{\"leaf1_ecx\":"
         << raw.Leaf1.Regs.Ecx
         << ",\"leaf1_edx\":" << raw.Leaf1.Regs.Edx
         << ",\"leaf6_eax\":" << raw.Leaf6.Regs.Eax
         << ",\"leaf6_ecx\":" << raw.Leaf6.Regs.Ecx
         << ",\"leaf80000001_ecx\":" << raw.Leaf80000001.Regs.Ecx
         << ",\"leaf80000001_edx\":" << raw.Leaf80000001.Regs.Edx
         << ",\"leaf80000008_ebx\":" << raw.Leaf80000008.Regs.Ebx
         << ",\"leaf80000021_eax\":" << raw.Leaf80000021.Regs.Eax
         << ",\"max_basic_leaf\":" << raw.Leaf0.Regs.Eax
         << ",\"max_extended_leaf\":" << raw.Leaf80000000.Regs.Eax
         << "},\"derived\":{\"ccd_id\":" << derived.CcdId
         << ",\"complex_id\":" << derived.ComplexId
         << ",\"core_id\":" << derived.CoreId
         << ",\"family\":" << derived.Family
         << ",\"full_apic_id\":" << derived.FullApicId
         << ",\"initial_apic_id_8\":" << derived.InitialApicId8
         << ",\"logical_processors_per_package\":"
         << derived.LogicalProcessorsPerPackage
         << ",\"model\":" << derived.Model
         << ",\"package_id\":" << derived.PackageId
         << ",\"stepping\":" << derived.Stepping
         << ",\"thread_id\":" << derived.ThreadId
         << ",\"threads_per_core\":" << derived.ThreadsPerCore
         << "},\"raw_capture_complete\":"
         << (receipt.RawCaptureComplete ? "true" : "false")
         << ",\"requested_cpu\":" << receipt.RequestedCpu
         << ",\"semantic_validation_passed\":"
         << (receipt.SemanticValidationPassed ? "true" : "false")
         << ",\"singleton_affinity_verified\":true}";
    diagnostic.clear();
    return json.str();
}

bool EmitConfig(
    std::ostream& output,
    const std::vector<std::shared_ptr<const ScreenCell> >& cells,
    const std::string& runtime_maps_sha256,
    const std::string& runtime_library_path,
    int target_cpu,
    const std::string& native_config_identity_sha256,
    const std::string& target_identity_json)
{
    if (!IsLowerSha256(runtime_maps_sha256.c_str()) ||
        runtime_library_path.empty() ||
        !IsLowerSha256(native_config_identity_sha256.c_str()) ||
        target_identity_json.empty() || target_identity_json.front() != '{' ||
        target_identity_json.back() != '}')
    {
        return false;
    }
    output << "{\"arm_descriptors\":[";
    const Arm arms[2] = { Arm::C, Arm::L };
    for (size_t i = 0u; i < 2u; ++i)
    {
        if (i != 0u) output << ',';
        const std::string descriptor = ArmDescriptorJson(arms[i]);
        output << "{\"arm\":\"" << ArmName(arms[i])
               << "\",\"descriptor\":" << descriptor
               << ",\"descriptor_sha256\":\""
               << wirehair::wh2_benchmark::Sha256Hex(descriptor)
               << "\"}";
    }
    output << "],\"campaign\":\"" << kCampaign << "\",\"cells\":[";
    for (size_t cell_index = 0u; cell_index < cells.size(); ++cell_index)
    {
        if (cell_index != 0u) output << ',';
        const ScreenCell& cell = *cells[cell_index];
        output << "{\"K\":" << cell.Shape.K
               << ",\"block_bytes\":" << cell.Shape.BlockBytes
               << ",\"final_bytes\":" << cell.Shape.BlockBytes
               << ",\"invocations_by_replicate\":[";
        for (uint32_t replicate = 0u;
             replicate < kPanelReplicates; ++replicate)
        {
            if (replicate != 0u) output << ',';
            output << InvocationsPerSlot(cell.Shape.K, replicate);
        }
        output << "],\"message_bytes\":" << cell.Source.size()
               << ",\"oracles\":[";
        for (size_t arm_index = 0u;
             arm_index < cell.Oracles.size(); ++arm_index)
        {
            if (arm_index != 0u) output << ',';
            const OracleReceipt& oracle = cell.Oracles[arm_index];
            output << "{\"arm\":\"" << ArmName(oracle.ArmValue)
                   << "\",\"attached_overlap_no_write_verified\":";
            if (oracle.ArmValue == Arm::C) output << "true";
            else output << "null";
            output << ",\"borrowed_eligible_systematic_ids\":"
                   << oracle.BorrowedEligibleSystematicIds
                   << ",\"descriptor_sha256\":\""
                   << oracle.DescriptorSha256
                   << "\",\"equation_configuration_sha256\":\""
                   << oracle.EquationConfigurationSha256
                   << "\",\"first_repair_sha256\":\""
                   << oracle.FirstRepairSha256
                   << "\",\"full_repair_sha256\":\""
                   << oracle.FullRepairSha256
                   << "\",\"high_id_sha256\":\""
                   << oracle.HighIdSha256
                   << "\",\"profile_bytes\":"
                   << oracle.ProfileBytes
                   << ",\"profile_hex\":";
            EmitNullableString(output, oracle.ProfileHex);
            output << ",\"profile_sha256\":";
            EmitNullableString(output, oracle.ProfileSha256);
            output << ",\"roundtrip_first_id\":"
                   << oracle.RoundtripFirstId
                   << ",\"roundtrip_packet_count\":"
                   << oracle.RoundtripPacketCount
                   << ",\"roundtrip_repair_only_verified\":"
                   << (oracle.RoundtripRepairOnlyVerified ? "true" : "false")
                   << ",\"roundtrip_sha256\":\""
                   << oracle.RoundtripSha256
                   << "\",\"roundtrip_verified\":"
                   << (oracle.RoundtripVerified ? "true" : "false")
                   << ",\"systematic_sha256\":\""
                   << oracle.SystematicSha256 << "\"}";
        }
        output << "],\"partial_final_oracles\":[";
        for (size_t partial_index = 0u;
             partial_index < cell.PartialFinalOracles.size(); ++partial_index)
        {
            if (partial_index != 0u) output << ',';
            const PartialFinalOracle& partial =
                cell.PartialFinalOracles[partial_index];
            output << "{\"arms\":[";
            for (size_t arm_index = 0u;
                 arm_index < partial.Arms.size(); ++arm_index)
            {
                if (arm_index != 0u) output << ',';
                const PartialArmOracle& arm = partial.Arms[arm_index];
                output << "{\"arm\":\"" << ArmName(arm.ArmValue)
                       << "\",\"attached_overlap_no_write_verified\":";
                if (arm.ArmValue == Arm::C) output << "true";
                else output << "null";
                output << ",\"first_repair_sha256\":\""
                       << arm.FirstRepairSha256
                       << "\",\"high_id_sha256\":\""
                       << arm.HighIdSha256
                       << "\",\"profile_sha256\":";
                EmitNullableString(output, arm.ProfileSha256);
                output << ",\"roundtrip_first_id\":"
                       << arm.RoundtripFirstId
                       << ",\"roundtrip_packet_count\":"
                       << arm.RoundtripPacketCount
                       << ",\"roundtrip_repair_only_verified\":"
                       << (arm.RoundtripRepairOnlyVerified ? "true" : "false")
                       << ",\"roundtrip_sha256\":\""
                       << arm.RoundtripSha256
                       << "\",\"source_immutable_verified\":"
                       << (arm.SourceImmutableVerified ? "true" : "false")
                       << ",\"systematic_sha256\":\""
                       << arm.SystematicSha256
                       << "\",\"systematic_verified\":"
                       << (arm.SystematicVerified ? "true" : "false")
                       << "}";
            }
            output << "],\"source_sha256\":\"" << partial.SourceSha256
                   << "\",\"tail_bytes\":" << partial.TailBytes << "}";
        }
        output << "],\"source_seed\":\"" << Hex64(cell.SourceSeed)
               << "\",\"source_sha256\":\"" << cell.SourceSha256
               << "\"}";
    }
    output << "],\"comparisons\":[";
    for (size_t i = 0u;
         i < sizeof(kComparisons) / sizeof(kComparisons[0]); ++i)
    {
        if (i != 0u) output << ',';
        output << "{\"left_arm\":\"" << ArmName(kComparisons[i].Left)
               << "\",\"name\":\"" << kComparisons[i].Name
               << "\",\"right_arm\":\""
               << ArmName(kComparisons[i].Right) << "\"}";
    }
    output << "],\"compile\":{\"compiler_path\":\""
           << WIREHAIR_WH2_CXX_COMPILER_PATH
           << "\",\"compiler_sha256\":\""
           << WIREHAIR_WH2_CXX_COMPILER_SHA256
           << "\",\"compiler_version\":\""
           << WIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION
           << "\",\"harness_git_commit\":\""
           << WIREHAIR_WH2_SOURCE_GIT_COMMIT
           << "\",\"implementation_git_commit\":\""
           << WIREHAIR_WH2_SOURCE_GIT_COMMIT
           << "\"},\"expected_encode_calls\":" << kExpectedEncodeCalls
           << ",\"expected_measured_encode_calls\":"
           << kExpectedMeasuredEncodeCalls
           << ",\"expected_measured_invocations\":"
           << kExpectedMeasuredInvocations
           << ",\"expected_panels\":" << kExpectedPanels
           << ",\"expected_records\":" << kExpectedRecords
           << ",\"expected_warmup_encode_calls\":"
           << kExpectedWarmupEncodeCalls
           << ",\"expected_warmup_invocations\":"
           << kExpectedWarmupInvocations
           << ",\"internal_deadline_seconds\":"
           << kInternalDeadlineSeconds
           << ",\"invocation_budget\":" << kInvocationBudget
           << ",\"minimum_invocations\":" << kMinimumInvocations
           << ",\"native_config_identity_sha256\":\""
           << native_config_identity_sha256 << "\""
           << ",\"panel_key_schema\":\"" << kPanelKeySchema
           << "\",\"panel_replicates\":" << kPanelReplicates
           << ",\"roster_order\":[\"replicate\",\"cell\",\"scope\","
              "\"comparison\"],\"runtime_library_maps_sha256\":\""
           << runtime_maps_sha256
           << "\",\"runtime_library_path\":\""
           << runtime_library_path
           << "\",\"schema\":\"" << kConfigSchema
           << "\",\"scopes\":[";
    for (size_t i = 0u; i < sizeof(kScopes) / sizeof(kScopes[0]); ++i)
    {
        if (i != 0u) output << ',';
        output << "{\"first_id\":\""
               << (IsSystematic(kScopes[i]) ? "0" : "K")
               << "\",\"name\":\"" << ScopeName(kScopes[i])
               << "\",\"timed_construction\":"
               << (IsFresh(kScopes[i]) ? "true" : "false") << "}";
    }
    output << "],\"sibling_cpu\":" << kSiblingCpu
           << ",\"source_seed_base\":\"" << Hex64(kSourceSeedBase)
           << "\",\"target_cpu\":" << target_cpu
           << ",\"target_identity\":" << target_identity_json << "}\n";
    output.flush();
    return output.good();
}

std::string ConfigIdentitySha256(
    const std::vector<std::shared_ptr<const ScreenCell> >& cells,
    const std::string& runtime_library_path,
    int target_cpu,
    const std::string& target_identity_json)
{
    const std::string zero_hash(64u, '0');
    std::ostringstream projection;
    if (!EmitConfig(
            projection, cells, zero_hash, runtime_library_path, target_cpu,
            zero_hash, target_identity_json))
    {
        return std::string();
    }
    return wirehair::wh2_benchmark::Sha256Hex(projection.str());
}

bool EmitPanel(
    const ScreenCell& cell,
    const Comparison& comparison,
    Scope scope,
    uint32_t replicate,
    uint64_t sequence,
    uint32_t invocations_per_slot,
    NativePanelOrder order,
    const std::string& panel_key_sha256,
    const std::string& runtime_maps_sha256,
    const PanelLane& left_lane,
    const PanelLane& right_lane,
    const NativePanelResult& panel)
{
    const OracleReceipt* left_oracle = FindOracle(cell, comparison.Left);
    const OracleReceipt* right_oracle = FindOracle(cell, comparison.Right);
    if (!left_oracle || !right_oracle) return false;
    const uint64_t calls_per_side =
        (UINT64_C(2) * invocations_per_slot + 1u) * cell.Shape.K;
    std::cout << "{\"K\":" << cell.Shape.K
              << ",\"affinity_verified\":true"
              << ",\"block_bytes\":" << cell.Shape.BlockBytes
              << ",\"comparison\":\"" << comparison.Name
              << "\",\"cpu_migration_verified\":true"
              << ",\"invocations_per_slot\":" << invocations_per_slot
              << ",\"left_arm\":\"" << ArmName(comparison.Left)
              << "\",\"left_descriptor_sha256\":\""
              << left_oracle->DescriptorSha256
              << "\",\"left_outcome_code\":0"
              << ",\"left_output_sha256\":\""
              << left_lane.OutputSha256()
              << "\",\"left_public_encode_calls\":" << calls_per_side
              << ",\"order\":\"" << OrderName(order)
              << "\",\"panel_key_sha256\":\"" << panel_key_sha256
              << "\",\"replicate\":" << replicate
              << ",\"right_arm\":\"" << ArmName(comparison.Right)
              << "\",\"right_descriptor_sha256\":\""
              << right_oracle->DescriptorSha256
              << "\",\"right_outcome_code\":0"
              << ",\"right_output_sha256\":\""
              << right_lane.OutputSha256()
              << "\",\"right_public_encode_calls\":" << calls_per_side
              << ",\"runtime_library_maps_sha256\":\""
              << runtime_maps_sha256
              << "\",\"schema\":\"" << kPanelSchema
              << "\",\"scope\":\"" << ScopeName(scope)
              << "\",\"sequence\":" << sequence
              << ",\"slots\":[";
    for (size_t i = 0u; i < panel.Slots.size(); ++i)
    {
        if (i != 0u) std::cout << ',';
        const uint32_t count = i < 4u ?
            invocations_per_slot / 2u + invocations_per_slot % 2u :
            invocations_per_slot / 2u;
        std::cout << "{\"elapsed_ns\":"
                  << panel.Slots[i].ElapsedNanoseconds
                  << ",\"invocation_count\":" << count
                  << ",\"side\":\"" << SideName(panel.Slots[i].Side)
                  << "\"}";
    }
    std::cout << "],\"target_cpu\":" << panel.TargetCpu << "}\n";
    std::cout.flush();
    return std::cout.good();
}

bool EmitTerminal(
    uint64_t encode_calls,
    uint64_t measured_invocations,
    uint32_t panels,
    uint64_t warmup_invocations)
{
    std::cout << "{\"encode_call_count\":" << encode_calls
              << ",\"measured_invocation_count\":"
              << measured_invocations
              << ",\"panel_count\":" << panels
              << ",\"record_count\":" << kExpectedRecords
              << ",\"schema\":\"" << kTerminalSchema
              << "\",\"status\":\"complete\""
              << ",\"warmup_invocation_count\":"
              << warmup_invocations << "}\n";
    std::cout.flush();
    return std::cout.good();
}

} // namespace

namespace {

const OracleReceipt* FindOracle(const ScreenCell& cell, Arm arm)
{
    for (const OracleReceipt& oracle : cell.Oracles) {
        if (oracle.ArmValue == arm) return &oracle;
    }
    return nullptr;
}

bool PrepareScratch(
    uint32_t K,
    uint32_t block_bytes,
    std::vector<uint8_t>& scratch,
    std::vector<uint32_t>& lengths)
{
    size_t bytes = 0u;
    if (!CheckedProduct(K, block_bytes, bytes)) return false;
    try
    {
        scratch.assign(bytes, kScratchCanary);
        lengths.assign(K, kLengthCanary);
    }
    catch (...) {
        return false;
    }
    return true;
}

int64_t EncodeBatch(
    const EncoderHandle& encoder,
    uint32_t first_id,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t final_bytes,
    std::vector<uint8_t>& scratch,
    std::vector<uint32_t>& lengths)
{
    size_t expected_size = 0u;
    if (!encoder.Valid() ||
        !CheckedProduct(K, block_bytes, expected_size) ||
        scratch.size() != expected_size || lengths.size() != K)
    {
        return WirehairV2_InvalidInput;
    }
    for (uint32_t offset = 0u; offset < K; ++offset)
    {
        const uint32_t id = first_id + offset;
        uint32_t& length = lengths[offset];
        const int64_t result = EncodeOne(
            encoder, id,
            scratch.data() + static_cast<size_t>(offset) * block_bytes,
            block_bytes, length);
        const uint32_t expected = ExpectedPacketBytes(
            id, K, block_bytes, final_bytes);
        if (result != 0 || length != expected) {
            return result == 0 ?
                static_cast<int64_t>(WirehairV2_Error) : result;
        }
    }
    return 0;
}

bool VerifyScratch(
    const std::vector<uint8_t>& source,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t final_bytes,
    bool systematic,
    const std::vector<uint8_t>& scratch,
    const std::vector<uint32_t>& lengths,
    std::string& sha256_out)
{
    size_t expected_size = 0u;
    if (!CheckedProduct(K, block_bytes, expected_size) ||
        scratch.size() != expected_size || lengths.size() != K)
    {
        return false;
    }
    for (uint32_t offset = 0u; offset < K; ++offset)
    {
        const uint32_t id = systematic ? offset : K + offset;
        if (lengths[offset] != ExpectedPacketBytes(
                id, K, block_bytes, final_bytes))
        {
            return false;
        }
    }
    if (systematic)
    {
        if (source.size() !=
                static_cast<size_t>(K - 1u) * block_bytes + final_bytes ||
            std::memcmp(scratch.data(), source.data(), source.size()) != 0)
        {
            return false;
        }
        const size_t suffix_start =
            static_cast<size_t>(K - 1u) * block_bytes + final_bytes;
        if (!std::all_of(
                scratch.begin() + static_cast<std::ptrdiff_t>(suffix_start),
                scratch.end(),
                [](uint8_t byte) { return byte == kScratchCanary; }))
        {
            return false;
        }
    }
    sha256_out = wirehair::wh2_benchmark::Sha256Hex(
        scratch.data(), scratch.size());
    return sha256_out.size() == 64u;
}

bool VerifySystematicBounds(
    const EncoderHandle& encoder,
    const std::vector<uint8_t>& source,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t final_bytes)
{
    static const size_t kGuardBytes = 32u;
    if (!encoder.Valid() || K == 0u || final_bytes == 0u ||
        final_bytes > block_bytes || source.size() !=
            static_cast<size_t>(K - 1u) * block_bytes + final_bytes)
    {
        return false;
    }
    static const size_t kCount = 2u;
    const uint32_t ids[kCount] = { 0u, K - 1u };
    const size_t stride =
        static_cast<size_t>(block_bytes) + 2u * kGuardBytes;
    std::vector<uint8_t> first;
    std::vector<uint8_t> second;
    try
    {
        first.assign(kCount * stride, UINT8_C(0xc3));
        second.assign(kCount * stride, UINT8_C(0x3c));
        for (size_t i = 0u; i < kCount; ++i)
        {
            std::fill(
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                kScratchCanary);
            std::fill(
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                UINT8_C(0x5a));
        }
    }
    catch (...) {
        return false;
    }
    for (size_t i = 0u; i < kCount; ++i)
    {
        uint8_t* const first_payload =
            first.data() + i * stride + kGuardBytes;
        uint8_t* const second_payload =
            second.data() + i * stride + kGuardBytes;
        const uint32_t expected_bytes = ids[i] + 1u == K ?
            final_bytes : block_bytes;
        const uint8_t* const expected = source.data() +
            static_cast<size_t>(ids[i]) * block_bytes;
        uint32_t first_bytes = kLengthCanary;
        uint32_t second_bytes = kLengthCanary;
        if (EncodeOne(
                encoder, ids[i], first_payload,
                block_bytes, first_bytes) != 0 ||
            EncodeOne(
                encoder, ids[i], second_payload,
                block_bytes, second_bytes) != 0 ||
            first_bytes != expected_bytes || second_bytes != expected_bytes ||
            !std::equal(expected, expected + expected_bytes, first_payload) ||
            !std::equal(expected, expected + expected_bytes, second_payload) ||
            !std::all_of(
                first_payload + expected_bytes,
                first_payload + block_bytes,
                [](uint8_t byte) { return byte == kScratchCanary; }) ||
            !std::all_of(
                second_payload + expected_bytes,
                second_payload + block_bytes,
                [](uint8_t byte) { return byte == UINT8_C(0x5a); }) ||
            !std::all_of(
                first.begin() + static_cast<std::ptrdiff_t>(i * stride),
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                [](uint8_t byte) { return byte == UINT8_C(0xc3); }) ||
            !std::all_of(
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                first.begin() + static_cast<std::ptrdiff_t>((i + 1u) * stride),
                [](uint8_t byte) { return byte == UINT8_C(0xc3); }) ||
            !std::all_of(
                second.begin() + static_cast<std::ptrdiff_t>(i * stride),
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                [](uint8_t byte) { return byte == UINT8_C(0x3c); }) ||
            !std::all_of(
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                second.begin() + static_cast<std::ptrdiff_t>((i + 1u) * stride),
                [](uint8_t byte) { return byte == UINT8_C(0x3c); }))
        {
            return false;
        }
    }
    return true;
}

bool VerifyHighIds(
    const EncoderHandle& encoder,
    uint32_t K,
    uint32_t block_bytes,
    std::string& digest_out)
{
    static const size_t kCount = 4u;
    static const size_t kGuardBytes = 32u;
    const uint32_t ids[kCount] = {
        K, K + 7u, UINT32_C(0x80000000), UINT32_MAX
    };
    std::vector<uint8_t> first;
    std::vector<uint8_t> second;
    std::vector<uint8_t> canonical;
    size_t payload_bytes = 0u;
    const size_t stride = static_cast<size_t>(block_bytes) +
        2u * kGuardBytes;
    size_t guarded_bytes = 0u;
    if (!CheckedProduct(kCount, block_bytes, payload_bytes) ||
        !CheckedProduct(kCount, stride, guarded_bytes))
    {
        return false;
    }
    try
    {
        first.assign(guarded_bytes, UINT8_C(0xc3));
        second.assign(guarded_bytes, UINT8_C(0x3c));
        canonical.resize(payload_bytes);
        for (size_t i = 0u; i < kCount; ++i)
        {
            std::fill(
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                kScratchCanary);
            std::fill(
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                UINT8_C(0x5a));
        }
    }
    catch (...) {
        return false;
    }
    for (size_t i = 0u; i < kCount; ++i)
    {
        uint8_t* const first_payload =
            first.data() + i * stride + kGuardBytes;
        uint8_t* const second_payload =
            second.data() + i * stride + kGuardBytes;
        uint32_t first_bytes = kLengthCanary;
        uint32_t second_bytes = kLengthCanary;
        if (EncodeOne(
                encoder, ids[i], first_payload,
                block_bytes, first_bytes) != 0 ||
            EncodeOne(
                encoder, ids[i], second_payload,
                block_bytes, second_bytes) != 0 ||
            first_bytes != block_bytes || second_bytes != block_bytes ||
            !std::equal(
                first_payload, first_payload + block_bytes, second_payload) ||
            !std::all_of(
                first.begin() + static_cast<std::ptrdiff_t>(i * stride),
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                [](uint8_t byte) { return byte == UINT8_C(0xc3); }) ||
            !std::all_of(
                first.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                first.begin() + static_cast<std::ptrdiff_t>((i + 1u) * stride),
                [](uint8_t byte) { return byte == UINT8_C(0xc3); }) ||
            !std::all_of(
                second.begin() + static_cast<std::ptrdiff_t>(i * stride),
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes),
                [](uint8_t byte) { return byte == UINT8_C(0x3c); }) ||
            !std::all_of(
                second.begin() + static_cast<std::ptrdiff_t>(
                    i * stride + kGuardBytes + block_bytes),
                second.begin() + static_cast<std::ptrdiff_t>((i + 1u) * stride),
                [](uint8_t byte) { return byte == UINT8_C(0x3c); }))
        {
            return false;
        }
        std::copy(
            first_payload, first_payload + block_bytes,
            canonical.begin() + static_cast<std::ptrdiff_t>(i * block_bytes));
    }
    digest_out = wirehair::wh2_benchmark::Sha256Hex(
        canonical.data(), canonical.size());
    return digest_out.size() == 64u;
}

bool VerifyRoundtrip(
    Arm arm,
    const EncoderHandle& encoder,
    const ConstructionArtifacts& artifacts,
    const std::vector<uint8_t>& source,
    uint32_t K,
    uint32_t block_bytes,
    uint32_t final_bytes,
    std::string& digest_out,
    uint32_t& packet_count_out)
{
    static const uint32_t kMaximumRepairOverhead = 64u;
    std::vector<uint8_t> packet(block_bytes, kScratchCanary);
    std::vector<uint8_t> recovered(source.size(), 0u);
    packet_count_out = 0u;
    if (arm == Arm::C)
    {
        WirehairV2Codec decoder = nullptr;
        if (wirehair_v2_decoder_create(
                artifacts.Profile.data(), artifacts.ProfileBytes,
                &decoder) != WirehairV2_Success || !decoder)
        {
            wirehair_v2_free(decoder);
            return false;
        }
        WirehairV2Result result = WirehairV2_NeedMore;
        for (uint32_t offset = 0u;
             offset < K + kMaximumRepairOverhead; ++offset)
        {
            const uint32_t id = K + offset;
            std::fill(packet.begin(), packet.end(), kScratchCanary);
            uint32_t bytes = kLengthCanary;
            if (EncodeOne(
                    encoder, id, packet.data(), block_bytes, bytes) != 0 ||
                bytes != ExpectedPacketBytes(
                    id, K, block_bytes, final_bytes))
            {
                wirehair_v2_free(decoder);
                return false;
            }
            result = wirehair_v2_decode(decoder, id, packet.data(), bytes);
            if (result == WirehairV2_Success)
            {
                packet_count_out = offset + 1u;
                break;
            }
            if (result != WirehairV2_NeedMore)
            {
                wirehair_v2_free(decoder);
                return false;
            }
        }
        uint64_t recovered_bytes = 0u;
        const bool success = result == WirehairV2_Success &&
            packet_count_out >= K &&
            wirehair_v2_recover(
                decoder, recovered.data(), recovered.size(),
                &recovered_bytes) == WirehairV2_Success &&
            recovered_bytes == recovered.size();
        wirehair_v2_free(decoder);
        if (!success) return false;
    }
    else if (arm == Arm::L)
    {
        WirehairCodec decoder = wirehair_decoder_create(
            nullptr, source.size(), block_bytes);
        if (!decoder) return false;
        WirehairResult result = Wirehair_NeedMore;
        for (uint32_t offset = 0u;
             offset < K + kMaximumRepairOverhead; ++offset)
        {
            const uint32_t id = K + offset;
            std::fill(packet.begin(), packet.end(), kScratchCanary);
            uint32_t bytes = kLengthCanary;
            if (EncodeOne(
                    encoder, id, packet.data(), block_bytes, bytes) != 0 ||
                bytes != ExpectedPacketBytes(
                    id, K, block_bytes, final_bytes))
            {
                wirehair_free(decoder);
                return false;
            }
            result = wirehair_decode(decoder, id, packet.data(), bytes);
            if (result == Wirehair_Success)
            {
                packet_count_out = offset + 1u;
                break;
            }
            if (result != Wirehair_NeedMore)
            {
                wirehair_free(decoder);
                return false;
            }
        }
        const bool success = result == Wirehair_Success &&
            packet_count_out >= K &&
            wirehair_recover(
                decoder, recovered.data(), recovered.size()) ==
                Wirehair_Success;
        wirehair_free(decoder);
        if (!success) return false;
    }
    else return false;
    if (recovered != source) return false;
    digest_out = wirehair::wh2_benchmark::Sha256Hex(
        recovered.data(), recovered.size());
    return digest_out.size() == 64u;
}

bool VerifyAttachedOverlapNoWrite(
    const EncoderHandle& encoder,
    const std::vector<uint8_t>& source,
    uint32_t block_bytes)
{
    if (!encoder.Valid() || !encoder.IsV2() ||
        source.size() < static_cast<size_t>(block_bytes) + 1u ||
        std::memcmp(source.data(), source.data() + 1u, block_bytes) == 0)
    {
        return false;
    }
    const std::string before = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    uint32_t bytes = kLengthCanary;
    const WirehairV2Result result = wirehair_v2_encode(
        encoder.V2Codec(), 0u,
        const_cast<uint8_t*>(source.data()) + 1u, block_bytes, &bytes);
    const std::string after = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    return result == WirehairV2_InvalidInput &&
        bytes == kLengthCanary && before.size() == 64u && before == after;
}

std::string EquationConfigurationSha256(
    Arm arm,
    const CellShape& shape,
    uint64_t message_bytes,
    const ConstructionArtifacts& artifacts)
{
    if (arm == Arm::C) {
        return wirehair::wh2_benchmark::Sha256Hex(
            artifacts.Profile.data(), artifacts.ProfileBytes);
    }
    if (arm != Arm::L) return std::string();
    std::ostringstream json;
    json << "{\"K\":" << shape.K
         << ",\"block_bytes\":" << shape.BlockBytes
         << ",\"codec\":\"wirehair1\",\"message_bytes\":"
         << message_bytes << "}";
    return wirehair::wh2_benchmark::Sha256Hex(json.str());
}

bool BuildOracle(
    Arm arm,
    const CellShape& shape,
    uint32_t final_bytes,
    const std::vector<uint8_t>& source,
    OracleReceipt& oracle)
{
    const std::string source_before = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    if (source_before.size() != 64u) return false;
    EncoderHandle encoder;
    ConstructionArtifacts artifacts;
    if (CreateEncoder(
            arm, source, shape.BlockBytes, encoder, artifacts) != 0 ||
        !encoder.Valid())
    {
        return false;
    }
    if (arm == Arm::C && !ValidateProfile(
            source, shape.BlockBytes, artifacts))
    {
        return false;
    }
    if (arm == Arm::L && artifacts.ProfileBytes != 0u) return false;

    std::vector<uint8_t> scratch;
    std::vector<uint32_t> lengths;
    if (!PrepareScratch(shape.K, shape.BlockBytes, scratch, lengths) ||
        EncodeBatch(
            encoder, 0u, shape.K, shape.BlockBytes, final_bytes,
            scratch, lengths) != 0 ||
        !VerifyScratch(
            source, shape.K, shape.BlockBytes, final_bytes, true,
            scratch, lengths, oracle.SystematicSha256) ||
        !VerifySystematicBounds(
            encoder, source, shape.K, shape.BlockBytes, final_bytes))
    {
        return false;
    }

    if (!PrepareScratch(shape.K, shape.BlockBytes, scratch, lengths) ||
        EncodeBatch(
            encoder, shape.K, shape.K, shape.BlockBytes, final_bytes,
            scratch, lengths) != 0 ||
        !VerifyScratch(
            source, shape.K, shape.BlockBytes, final_bytes, false,
            scratch, lengths, oracle.FullRepairSha256))
    {
        return false;
    }
    oracle.FirstRepairSha256 = wirehair::wh2_benchmark::Sha256Hex(
        scratch.data(), shape.BlockBytes);
    std::string repeated_repair_sha256;
    std::fill(scratch.begin(), scratch.end(), UINT8_C(0x5a));
    std::fill(lengths.begin(), lengths.end(), kLengthCanary);
    if (oracle.FirstRepairSha256.size() != 64u ||
        EncodeBatch(
            encoder, shape.K, shape.K, shape.BlockBytes, final_bytes,
            scratch, lengths) != 0 ||
        !VerifyScratch(
            source, shape.K, shape.BlockBytes, final_bytes, false,
            scratch, lengths, repeated_repair_sha256) ||
        repeated_repair_sha256 != oracle.FullRepairSha256 ||
        wirehair::wh2_benchmark::Sha256Hex(
            scratch.data(), shape.BlockBytes) != oracle.FirstRepairSha256 ||
        !VerifyHighIds(
            encoder, shape.K, shape.BlockBytes, oracle.HighIdSha256) ||
        !VerifyRoundtrip(
            arm, encoder, artifacts, source, shape.K, shape.BlockBytes,
            final_bytes, oracle.RoundtripSha256,
            oracle.RoundtripPacketCount))
    {
        return false;
    }

    oracle.ArmValue = arm;
    oracle.AttachedOverlapNoWriteVerified = arm == Arm::C ?
        VerifyAttachedOverlapNoWrite(encoder, source, shape.BlockBytes) :
        false;
    if (arm == Arm::C && !oracle.AttachedOverlapNoWriteVerified) return false;
    oracle.BorrowedEligibleSystematicIds = arm == Arm::C ? shape.K : 0u;
    const std::string descriptor = ArmDescriptorJson(arm);
    oracle.DescriptorSha256 =
        wirehair::wh2_benchmark::Sha256Hex(descriptor);
    oracle.EquationConfigurationSha256 = EquationConfigurationSha256(
        arm, shape, source.size(), artifacts);
    if (arm == Arm::C)
    {
        oracle.ProfileBytes = artifacts.ProfileBytes;
        oracle.ProfileHex = HexBytes(
            artifacts.Profile.data(), artifacts.ProfileBytes);
        oracle.ProfileSha256 = wirehair::wh2_benchmark::Sha256Hex(
            artifacts.Profile.data(), artifacts.ProfileBytes);
    }
    oracle.RoundtripVerified = oracle.RoundtripSha256 == source_before;
    oracle.RoundtripFirstId = shape.K;
    oracle.RoundtripRepairOnlyVerified = true;
    const std::string source_after = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    return oracle.DescriptorSha256.size() == 64u &&
        oracle.EquationConfigurationSha256.size() == 64u &&
        oracle.RoundtripVerified && oracle.RoundtripPacketCount >= shape.K &&
        oracle.RoundtripPacketCount <= shape.K + 64u &&
        source_after == source_before &&
        (arm != Arm::C ||
            (oracle.ProfileBytes == WIREHAIR_V2_PROFILE_SERIALIZED_BYTES &&
             oracle.ProfileHex.size() == 2u * oracle.ProfileBytes &&
             oracle.ProfileSha256.size() == 64u));
}

bool BuildPartialFinalOracle(
    const CellShape& shape,
    uint32_t tail_bytes,
    PartialFinalOracle& output)
{
    std::vector<uint8_t> source;
    if (!MakeSource(shape, tail_bytes, source)) return false;
    output.TailBytes = tail_bytes;
    output.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    if (output.SourceSha256.size() != 64u) return false;
    const Arm arms[2] = { Arm::C, Arm::L };
    for (size_t i = 0u; i < 2u; ++i)
    {
        OracleReceipt oracle;
        if (!BuildOracle(
                arms[i], shape, tail_bytes, source, oracle))
        {
            return false;
        }
        PartialArmOracle& partial = output.Arms[i];
        partial.ArmValue = arms[i];
        partial.AttachedOverlapNoWriteVerified =
            oracle.AttachedOverlapNoWriteVerified;
        partial.FirstRepairSha256 = oracle.FirstRepairSha256;
        partial.HighIdSha256 = oracle.HighIdSha256;
        partial.ProfileSha256 = oracle.ProfileSha256;
        partial.RoundtripFirstId = oracle.RoundtripFirstId;
        partial.RoundtripPacketCount = oracle.RoundtripPacketCount;
        partial.RoundtripRepairOnlyVerified =
            oracle.RoundtripRepairOnlyVerified;
        partial.RoundtripSha256 = oracle.RoundtripSha256;
        partial.SourceImmutableVerified = true;
        partial.SystematicSha256 = oracle.SystematicSha256;
        partial.SystematicVerified = true;
    }
    return output.Arms[0].SystematicSha256 ==
            output.Arms[1].SystematicSha256 &&
        output.Arms[0].FirstRepairSha256 !=
            output.Arms[1].FirstRepairSha256 &&
        output.Arms[0].HighIdSha256 != output.Arms[1].HighIdSha256;
}

bool BuildCell(const CellShape& shape, ScreenCell& output)
{
    ScreenCell cell;
    cell.Shape = shape;
    cell.SourceSeed = SourceSeed(shape, shape.BlockBytes);
    if (!MakeSource(shape, shape.BlockBytes, cell.Source)) return false;
    cell.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Source.data(), cell.Source.size());
    if (cell.SourceSha256.size() != 64u ||
        !BuildOracle(
            Arm::C, shape, shape.BlockBytes, cell.Source,
            cell.Oracles[0]) ||
        !BuildOracle(
            Arm::L, shape, shape.BlockBytes, cell.Source,
            cell.Oracles[1]) ||
        cell.Oracles[0].SystematicSha256 !=
            cell.Oracles[1].SystematicSha256 ||
        cell.Oracles[0].FullRepairSha256 ==
            cell.Oracles[1].FullRepairSha256 ||
        !BuildPartialFinalOracle(shape, 1u, cell.PartialFinalOracles[0]) ||
        !BuildPartialFinalOracle(
            shape, shape.BlockBytes - 1u, cell.PartialFinalOracles[1]))
    {
        return false;
    }
    output = std::move(cell);
    return true;
}

bool BuildAllCells(
    std::vector<std::shared_ptr<const ScreenCell> >& cells)
{
    cells.clear();
    try {
        cells.reserve(sizeof(kCellShapes) / sizeof(kCellShapes[0]));
    }
    catch (...) {
        return false;
    }
    for (const CellShape& shape : kCellShapes)
    {
        std::shared_ptr<ScreenCell> cell;
        try {
            cell.reset(new ScreenCell);
        }
        catch (...) {
            return false;
        }
        if (!BuildCell(shape, *cell)) return false;
        cells.push_back(cell);
    }
    return cells.size() == sizeof(kCellShapes) / sizeof(kCellShapes[0]);
}

} // namespace

namespace {

bool SelfTest()
{
    // Deliberately no call to RunScreen(), ExecuteNativeTimingPanel(), or the
    // 1,440-panel roster exists on this path. This mode exercises arithmetic,
    // canonical receipts, and one small untimed public-API oracle only.
    if (wirehair_init() != Wirehair_Success)
    {
        std::cerr << "selftest: wirehair initialization failed\n";
        return false;
    }
    if (!TestRosterArithmetic())
    {
        std::cerr << "selftest: frozen roster arithmetic failed\n";
        return false;
    }
    if (!TestCpuParserAndMatcher())
    {
        std::cerr << "selftest: CPU parser/matcher failed\n";
        return false;
    }
    if (!TestLaunchAuthorization())
    {
        std::cerr << "selftest: one-shot launch authorization failed\n";
        return false;
    }
    if (std::strcmp(WIREHAIR_WH2_SOURCE_GIT_COMMIT, "") == 0 ||
        std::strcmp(WIREHAIR_WH2_CMAKE_CXX_COMPILER_ID, "") == 0 ||
        std::strcmp(WIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION, "") == 0 ||
        !ArmDescriptorJson(Arm::Invalid).empty())
    {
        std::cerr << "selftest: build identity or arm descriptor failed\n";
        return false;
    }
    ScreenCell small;
    if (!BuildCell(kCellShapes[0], small) ||
        small.Source.size() !=
            static_cast<size_t>(small.Shape.K) * small.Shape.BlockBytes ||
        small.PartialFinalOracles[0].TailBytes != 1u ||
        small.PartialFinalOracles[1].TailBytes !=
            small.Shape.BlockBytes - 1u)
    {
        std::cerr << "selftest: public partial/full oracle closure failed\n";
        return false;
    }
    // The selftest must remain valid for ordinary static CI builds. The real
    // run path independently requires a nonempty single-DSO maps receipt.
    const std::string maps_sha256(64u, 'a');
    std::vector<std::shared_ptr<const ScreenCell> > cells;
    cells.push_back(std::shared_ptr<const ScreenCell>(
        new ScreenCell(std::move(small))));
    std::ostringstream config;
    const std::string synthetic_library_path(
        "/synthetic/libwirehair.so.2.0.0");
    const std::string synthetic_target_identity =
        "{\"after_cpu\":120,\"before_cpu\":120,\"canonical_bytes\":1,"
        "\"canonical_hex\":\"00\",\"canonical_sha256\":\""
        + std::string(64u, 'c') +
        "\",\"capabilities\":{},\"derived\":{},"
        "\"raw_capture_complete\":true,\"requested_cpu\":120,"
        "\"semantic_validation_passed\":true,"
        "\"singleton_affinity_verified\":true}";
    const std::string native_config_identity_sha256 =
        ConfigIdentitySha256(
            cells, synthetic_library_path, kTargetCpu,
            synthetic_target_identity);
    if (!EmitConfig(
            config, cells, maps_sha256, synthetic_library_path, kTargetCpu,
            native_config_identity_sha256, synthetic_target_identity))
    {
        std::cerr << "selftest: canonical config emission failed\n";
        return false;
    }
    const std::string json = config.str();
    const std::string prefix = "{\"arm_descriptors\":[";
    const std::string boundary =
        std::string("\"runtime_library_maps_sha256\":\"") +
        maps_sha256 + "\",\"runtime_library_path\":\"" +
        synthetic_library_path + "\",\"schema\":\"" + kConfigSchema +
        "\"";
    const bool valid = json.size() > prefix.size() &&
        json.compare(0u, prefix.size(), prefix) == 0 &&
        json.back() == '\n' &&
        std::count(json.begin(), json.end(), '\n') == 1 &&
        json.find(boundary) != std::string::npos;
    if (!valid) std::cerr << "selftest: canonical config shape failed\n";
    return valid;
}

bool EmitNativeConfigSelfTest()
{
    // This is an untimed, non-authoritative integration oracle. It pins and
    // brackets CPU identity but never executes a timing panel or touches the
    // fixed namespace.
    std::string diagnostic;
    if (!PinAndVerifyTargetCpu(kTargetCpu, diagnostic))
    {
        std::cerr << "native config selftest: affinity failed: "
                  << diagnostic << '\n';
        return false;
    }
    wirehair_wh2_bench::TargetIdentityReceiptV2 target_identity;
    if (!wirehair_wh2_bench::CapturePublicBorrowedTargetIdentity(
            kTargetCpu, target_identity, diagnostic))
    {
        std::cerr << "native config selftest: identity failed: "
                  << diagnostic << '\n';
        return false;
    }
    const std::string target_identity_json =
        TargetIdentityJson(target_identity, diagnostic);
    if (target_identity_json.empty())
    {
        std::cerr << "native config selftest: identity receipt failed: "
                  << diagnostic << '\n';
        return false;
    }
    if (wirehair_init() != Wirehair_Success)
    {
        std::cerr << "native config selftest: wirehair initialization failed\n";
        return false;
    }
    std::vector<std::shared_ptr<const ScreenCell> > cells;
    if (!BuildAllCells(cells))
    {
        std::cerr << "native config selftest: oracle closure failed\n";
        return false;
    }
    wirehair_wh2_bench::TargetIdentityReceiptV2 identity_after_oracles;
    if (!wirehair_wh2_bench::CapturePublicBorrowedTargetIdentity(
            kTargetCpu, identity_after_oracles, diagnostic) ||
        TargetIdentityJson(identity_after_oracles, diagnostic) !=
            target_identity_json)
    {
        std::cerr << "native config selftest: identity drift: "
                  << diagnostic << '\n';
        return false;
    }
    std::string runtime_maps_text;
    std::string runtime_library_path;
    std::string binding_path;
    if (!RuntimeWirehairMapsReceipt(
            runtime_maps_text, runtime_library_path) ||
        !RuntimeTimedBindingPath(binding_path) ||
        binding_path != runtime_library_path)
    {
        std::cerr << "native config selftest: DSO binding differs\n";
        return false;
    }
    const std::string runtime_maps_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(runtime_maps_text);
    const std::string native_config_identity_sha256 =
        ConfigIdentitySha256(
            cells, runtime_library_path, kTargetCpu, target_identity_json);
    return runtime_maps_sha256.size() == 64u &&
        native_config_identity_sha256.size() == 64u && EmitConfig(
            std::cout, cells, runtime_maps_sha256, runtime_library_path,
            kTargetCpu, native_config_identity_sha256,
            target_identity_json);
}

bool RunScreen(
    int target_cpu,
    const std::string& expected_native_config_identity_sha256)
{
    std::string diagnostic;
    if (!PinAndVerifyTargetCpu(target_cpu, diagnostic))
    {
        std::cerr << "public borrowed r1 affinity failed: "
                  << diagnostic << '\n';
        return false;
    }
    wirehair_wh2_bench::TargetIdentityReceiptV2 target_identity;
    if (!wirehair_wh2_bench::CapturePublicBorrowedTargetIdentity(
            target_cpu, target_identity, diagnostic))
    {
        std::cerr << "public borrowed r1 target identity failed: "
                  << diagnostic << '\n';
        return false;
    }
    const std::string target_identity_json =
        TargetIdentityJson(target_identity, diagnostic);
    if (target_identity_json.empty())
    {
        std::cerr << "public borrowed r1 target identity receipt failed: "
                  << diagnostic << '\n';
        return false;
    }
    if (wirehair_init() != Wirehair_Success)
    {
        std::cerr << "public borrowed r1 wirehair initialization failed\n";
        return false;
    }

    std::vector<std::shared_ptr<const ScreenCell> > cells;
    if (!BuildAllCells(cells))
    {
        std::cerr << "public borrowed r1 oracle closure failed\n";
        return false;
    }
    wirehair_wh2_bench::TargetIdentityReceiptV2 identity_after_oracles;
    if (!wirehair_wh2_bench::CapturePublicBorrowedTargetIdentity(
            target_cpu, identity_after_oracles, diagnostic) ||
        TargetIdentityJson(identity_after_oracles, diagnostic) !=
            target_identity_json)
    {
        std::cerr << "public borrowed r1 pre-timing identity drift: "
                  << diagnostic << '\n';
        return false;
    }
    if (!VerifyTargetCpu(target_cpu))
    {
        std::cerr << "public borrowed r1 pre-timing CPU drift\n";
        return false;
    }
    std::string runtime_maps_text;
    std::string runtime_library_path;
    std::string binding_path;
    if (!RuntimeWirehairMapsReceipt(
            runtime_maps_text, runtime_library_path) ||
        !RuntimeTimedBindingPath(binding_path) ||
        binding_path != runtime_library_path)
    {
        std::cerr << "public borrowed r1 shared-library binding invalid\n";
        return false;
    }
    const std::string runtime_maps_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(runtime_maps_text);
    const std::string native_config_identity_sha256 =
        ConfigIdentitySha256(
            cells, runtime_library_path, target_cpu, target_identity_json);
    if (runtime_maps_text.empty() || runtime_maps_sha256.size() != 64u ||
        native_config_identity_sha256.size() != 64u ||
        native_config_identity_sha256 !=
            expected_native_config_identity_sha256)
    {
        std::cerr << "public borrowed r1 config identity invalid\n";
        return false;
    }
    if (!EmitConfig(
            std::cout, cells, runtime_maps_sha256, runtime_library_path,
            target_cpu, native_config_identity_sha256,
            target_identity_json))
    {
        std::cerr << "public borrowed r1 config write failed\n";
        return false;
    }

    const Clock::time_point deadline =
        Clock::now() + std::chrono::seconds(kInternalDeadlineSeconds);
    uint32_t panel_count = 0u;
    uint64_t measured_invocations = 0u;
    uint64_t warmup_invocations = 0u;
    uint64_t encode_calls = 0u;
    uint64_t sequence = 0u;
    for (uint32_t replicate = 0u;
         replicate < kPanelReplicates; ++replicate)
    {
        for (const std::shared_ptr<const ScreenCell>& cell : cells)
        {
            const uint32_t n = InvocationsPerSlot(
                cell->Shape.K, replicate);
            for (Scope scope : kScopes)
            {
                for (const Comparison& comparison : kComparisons)
                {
                    if (Clock::now() >= deadline)
                    {
                        std::cerr <<
                            "public borrowed r1 internal deadline\n";
                        return false;
                    }
                    const std::string panel_key_sha256 = PanelKeySha256(
                        cell->Shape, comparison, scope);
                    const NativePanelOrder order = PanelOrder(
                        cell->Shape, comparison, scope, replicate);
                    PanelLane left_lane(
                        *cell, comparison.Left, scope, "left");
                    PanelLane right_lane(
                        *cell, comparison.Right, scope, "right");
                    if (panel_key_sha256.size() != 64u ||
                        (order != NativePanelOrder::ABBA &&
                         order != NativePanelOrder::BAAB) ||
                        !left_lane.Initialize() ||
                        !right_lane.Initialize())
                    {
                        std::cerr <<
                            "public borrowed r1 panel setup failed\n";
                        return false;
                    }
                    const NativePanelResult panel =
                        wirehair_wh2_bench::ExecuteNativeTimingPanel(
                            target_cpu, order, n,
                            MakePanelArm(left_lane),
                            MakePanelArm(right_lane));
                    if (!VerifyPanelResult(
                            panel, cell->Shape.K, n, target_cpu) ||
                        !left_lane.VerifyAfterPanel() ||
                        !right_lane.VerifyAfterPanel() ||
                        wirehair::wh2_benchmark::Sha256Hex(
                            cell->Source.data(), cell->Source.size()) !=
                            cell->SourceSha256 ||
                        !VerifyTargetCpu(target_cpu) ||
                        RuntimeWirehairMapsText() != runtime_maps_text)
                    {
                        std::cerr << "public borrowed r1 panel failed: "
                                  << wirehair_wh2_bench::NativePanelStatusName(
                                         panel.Status)
                                  << ' ' << panel.Diagnostic << '\n';
                        return false;
                    }
                    if (!EmitPanel(
                            *cell, comparison, scope, replicate, sequence,
                            n, order, panel_key_sha256,
                            runtime_maps_sha256, left_lane, right_lane,
                            panel))
                    {
                        std::cerr << "public borrowed r1 panel write failed\n";
                        return false;
                    }
                    measured_invocations += UINT64_C(4) * n;
                    warmup_invocations += 2u;
                    encode_calls +=
                        (UINT64_C(4) * n + 2u) * cell->Shape.K;
                    ++panel_count;
                    ++sequence;
                }
            }
        }
    }
    if (Clock::now() >= deadline || !VerifyTargetCpu(target_cpu) ||
        RuntimeWirehairMapsText() != runtime_maps_text ||
        panel_count != kExpectedPanels || sequence != kExpectedPanels ||
        measured_invocations != kExpectedMeasuredInvocations ||
        warmup_invocations != kExpectedWarmupInvocations ||
        encode_calls != kExpectedEncodeCalls)
    {
        std::cerr << "public borrowed r1 terminal closure failed\n";
        return false;
    }
    for (const std::shared_ptr<const ScreenCell>& cell : cells)
    {
        if (wirehair::wh2_benchmark::Sha256Hex(
                cell->Source.data(), cell->Source.size()) !=
                cell->SourceSha256)
        {
            std::cerr << "public borrowed r1 terminal source drift\n";
            return false;
        }
    }
    if (Clock::now() >= deadline || !VerifyTargetCpu(target_cpu) ||
        RuntimeWirehairMapsText() != runtime_maps_text)
    {
        std::cerr << "public borrowed r1 final deadline/binding failed\n";
        return false;
    }
    wirehair_wh2_bench::TargetIdentityReceiptV2 identity_after_timing;
    if (!wirehair_wh2_bench::CapturePublicBorrowedTargetIdentity(
            target_cpu, identity_after_timing, diagnostic) ||
        TargetIdentityJson(identity_after_timing, diagnostic) !=
            target_identity_json || Clock::now() >= deadline)
    {
        std::cerr << "public borrowed r1 terminal identity/deadline failed: "
                  << diagnostic << '\n';
        return false;
    }
    return EmitTerminal(
        encode_calls, measured_invocations, panel_count,
        warmup_invocations);
}

} // namespace
