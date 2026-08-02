#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"
#include "Wh2NativePanel.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

using wirehair::wh2_benchmark::FrozenPacketTrace;
using wirehair::wh2_benchmark::FrozenSchedule;
using wirehair::wh2_benchmark::FrozenTraceStatus;
using wirehair_wh2_bench::IsolatedSolveResult;
using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativePanelArm;
using wirehair_wh2_bench::NativePanelDisposition;
using wirehair_wh2_bench::NativePanelInvocation;
using wirehair_wh2_bench::NativePanelInvocationResult;
using wirehair_wh2_bench::NativePanelOrder;
using wirehair_wh2_bench::NativePanelResult;
using wirehair_wh2_bench::NativePanelSide;
using wirehair_wh2_bench::NativePanelStatus;
using wirehair_wh2_bench::NativeReceiveFixture;
using wirehair_wh2_bench::NativeSolveFixture;
using wirehair_wh2_bench::NativeWh2ExecutionMode;
using wirehair_wh2_bench::ResolvedNativeWh2Configuration;
using wirehair_wh2_bench::TimedArmResult;

static const char kSchema[] = "wirehair.wh2.tiny-native-cell.v1";
static const char kEquationSchema[] =
    "wirehair.wh2.tiny-native-equations.v1";
static const char kConstructionSchema[] =
    "wirehair.wh2.tiny-native-construction.v1";
static const uint32_t kLossPpm = 100000u;
static const uint32_t kTraceOverhead = 256u;
static const uint32_t kSolveOverhead = 4u;
static const uint32_t kMaxInvocationsPerSlot = 32768u;

enum class Layout
{
    Disabled,
    Two07
};

enum class Scope
{
    Solve,
    Receive,
    Encoder
};

enum class ArmMode
{
    Wirehair1,
    DirectOff,
    DirectOn
};

struct Request
{
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    Layout LayoutValue = Layout::Disabled;
    Scope ScopeValue = Scope::Solve;
    ArmMode LeftMode = ArmMode::DirectOff;
    ArmMode RightMode = ArmMode::DirectOn;
    uint32_t ConstructionAttempt = 0u;
    uint64_t TraceSeed = 0u;
    uint32_t Replicate = 0u;
    int Cpu = -1;
    uint32_t InvocationsPerSlot = 0u;
};

const char* LayoutName(Layout layout)
{
    return layout == Layout::Disabled ? "disabled" : "two07";
}

const char* ScopeName(Scope scope)
{
    switch (scope)
    {
    case Scope::Solve: return "solve";
    case Scope::Receive: return "receive";
    case Scope::Encoder: return "encoder";
    default: return "unknown";
    }
}

const char* ArmModeName(ArmMode mode)
{
    switch (mode)
    {
    case ArmMode::Wirehair1: return "wh1";
    case ArmMode::DirectOff: return "off";
    case ArmMode::DirectOn: return "on";
    default: return "unknown";
    }
}

const char* OrderName(NativePanelOrder order)
{
    return order == NativePanelOrder::ABBA ? "ABBA" : "BAAB";
}

bool ParseUnsigned(const std::string& text, uint64_t maximum, uint64_t& value)
{
    if (text.empty() || (text.size() > 1u && text[0] == '0')) {
        return false;
    }
    uint64_t parsed = 0u;
    for (char ch : text)
    {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const uint64_t digit = static_cast<uint64_t>(ch - '0');
        if (digit > maximum || parsed > (maximum - digit) / 10u) {
            return false;
        }
        parsed = parsed * 10u + digit;
    }
    value = parsed;
    return true;
}

bool ParseLayout(const std::string& text, Layout& layout)
{
    if (text == "disabled") {
        layout = Layout::Disabled;
        return true;
    }
    if (text == "two07") {
        layout = Layout::Two07;
        return true;
    }
    return false;
}

bool ParseScope(const std::string& text, Scope& scope)
{
    if (text == "solve") {
        scope = Scope::Solve;
        return true;
    }
    if (text == "receive") {
        scope = Scope::Receive;
        return true;
    }
    if (text == "encoder") {
        scope = Scope::Encoder;
        return true;
    }
    return false;
}

bool ParseArmMode(const std::string& text, ArmMode& mode)
{
    if (text == "wh1") {
        mode = ArmMode::Wirehair1;
        return true;
    }
    if (text == "off") {
        mode = ArmMode::DirectOff;
        return true;
    }
    if (text == "on") {
        mode = ArmMode::DirectOn;
        return true;
    }
    return false;
}

void PrintUsage(std::ostream& output)
{
    output << "usage: wirehair_wh2_tiny_native_cell"
        " --k <2..8> --block-bytes <1..4096>"
        " --layout <disabled|two07>"
        " --scope <solve|receive|encoder>"
        " --left <wh1|off|on> --right <wh1|off|on>"
        " --construction-attempt <0..255>"
        " --trace-seed <uint64-decimal> --replicate <uint32>"
        " --cpu <nonnegative> --invocations-per-slot <2..32768>\n";
}

bool ParseRequest(int argc, char** argv, Request& request, std::string& error)
{
    if (argc != 23)
    {
        error = "exactly eleven option/value pairs are required";
        return false;
    }

    bool have_k = false;
    bool have_block_bytes = false;
    bool have_layout = false;
    bool have_scope = false;
    bool have_left = false;
    bool have_right = false;
    bool have_attempt = false;
    bool have_trace_seed = false;
    bool have_replicate = false;
    bool have_cpu = false;
    bool have_invocations = false;
    for (int i = 1; i < argc; i += 2)
    {
        const std::string option(argv[i]);
        const std::string value(argv[i + 1]);
        uint64_t parsed = 0u;
        if (option == "--k" && !have_k &&
            ParseUnsigned(value, 8u, parsed) && parsed >= 2u)
        {
            request.K = static_cast<uint32_t>(parsed);
            have_k = true;
        }
        else if (option == "--block-bytes" && !have_block_bytes &&
            ParseUnsigned(value, 4096u, parsed) && parsed >= 1u)
        {
            request.BlockBytes = static_cast<uint32_t>(parsed);
            have_block_bytes = true;
        }
        else if (option == "--layout" && !have_layout &&
            ParseLayout(value, request.LayoutValue))
        {
            have_layout = true;
        }
        else if (option == "--scope" && !have_scope &&
            ParseScope(value, request.ScopeValue))
        {
            have_scope = true;
        }
        else if (option == "--left" && !have_left &&
            ParseArmMode(value, request.LeftMode))
        {
            have_left = true;
        }
        else if (option == "--right" && !have_right &&
            ParseArmMode(value, request.RightMode))
        {
            have_right = true;
        }
        else if (option == "--construction-attempt" && !have_attempt &&
            ParseUnsigned(value, 255u, parsed))
        {
            request.ConstructionAttempt = static_cast<uint32_t>(parsed);
            have_attempt = true;
        }
        else if (option == "--trace-seed" && !have_trace_seed &&
            ParseUnsigned(value, UINT64_MAX, parsed))
        {
            request.TraceSeed = parsed;
            have_trace_seed = true;
        }
        else if (option == "--replicate" && !have_replicate &&
            ParseUnsigned(value, UINT32_MAX, parsed))
        {
            request.Replicate = static_cast<uint32_t>(parsed);
            have_replicate = true;
        }
        else if (option == "--cpu" && !have_cpu &&
            ParseUnsigned(value, static_cast<uint64_t>(INT_MAX), parsed))
        {
            request.Cpu = static_cast<int>(parsed);
            have_cpu = true;
        }
        else if (option == "--invocations-per-slot" && !have_invocations &&
            ParseUnsigned(value, kMaxInvocationsPerSlot, parsed) &&
            parsed >= 2u)
        {
            request.InvocationsPerSlot = static_cast<uint32_t>(parsed);
            have_invocations = true;
        }
        else
        {
            error = "unknown, duplicate, or invalid option: " + option;
            return false;
        }
    }

    if (!have_k || !have_block_bytes || !have_layout || !have_scope ||
        !have_left || !have_right || !have_attempt || !have_trace_seed ||
        !have_replicate || !have_cpu || !have_invocations)
    {
        error = "all options are required exactly once";
        return false;
    }
    if (request.ScopeValue == Scope::Solve &&
        (request.LeftMode == ArmMode::Wirehair1 ||
         request.RightMode == ArmMode::Wirehair1))
    {
        error = "solve scope does not support Wirehair1";
        return false;
    }
    return true;
}

std::string Hex64(uint64_t value)
{
    static const char digits[] = "0123456789abcdef";
    std::string text(18u, '0');
    text[0] = '0';
    text[1] = 'x';
    for (unsigned i = 0u; i < 16u; ++i) {
        text[17u - i] = digits[value & 15u];
        value >>= 4u;
    }
    return text;
}

std::string JsonString(const std::string& value)
{
    static const char digits[] = "0123456789abcdef";
    std::string escaped;
    escaped.reserve(value.size() + 2u);
    escaped.push_back('"');
    for (unsigned char ch : value)
    {
        switch (ch)
        {
        case '"': escaped += "\\\""; break;
        case '\\': escaped += "\\\\"; break;
        case '\b': escaped += "\\b"; break;
        case '\f': escaped += "\\f"; break;
        case '\n': escaped += "\\n"; break;
        case '\r': escaped += "\\r"; break;
        case '\t': escaped += "\\t"; break;
        default:
            if (ch < 0x20u)
            {
                escaped += "\\u00";
                escaped.push_back(digits[ch >> 4u]);
                escaped.push_back(digits[ch & 15u]);
            }
            else {
                escaped.push_back(static_cast<char>(ch));
            }
        }
    }
    escaped.push_back('"');
    return escaped;
}

std::string CanonicalRequestJson(const Request& request)
{
    std::ostringstream json;
    json << "{\"K\":" << request.K
         << ",\"block_bytes\":" << request.BlockBytes
         << ",\"construction_attempt\":"
         << request.ConstructionAttempt
         << ",\"invocations_per_slot\":"
         << request.InvocationsPerSlot
         << ",\"layout\":\"" << LayoutName(request.LayoutValue)
         << "\",\"left_mode\":\"" << ArmModeName(request.LeftMode)
         << "\",\"loss_ppm\":" << kLossPpm
         << ",\"replicate\":" << request.Replicate
         << ",\"right_mode\":\"" << ArmModeName(request.RightMode)
         << "\",\"schedule\":\"iid\",\"scope\":\""
         << ScopeName(request.ScopeValue)
         << "\",\"target_cpu\":" << request.Cpu
         << ",\"trace_seed\":\"" << Hex64(request.TraceSeed)
         << "\"}";
    return json.str();
}

bool TwoAnchorTransform(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig&,
    void*)
{
    if (block_count != params.BlockCount || block_bytes == 0u ||
        params.DenseAnchors != wirehair_v2::DenseAnchorLayout::Disabled)
    {
        return false;
    }
    params.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    return true;
}

const char* DenseAnchorName(wirehair_v2::DenseAnchorLayout layout)
{
    switch (layout)
    {
    case wirehair_v2::DenseAnchorLayout::Disabled: return "disabled";
    case wirehair_v2::DenseAnchorLayout::Two07: return "two07";
    case wirehair_v2::DenseAnchorLayout::Four0369: return "four0369";
    default: return "unknown";
    }
}

const char* HeavyFamilyName(wirehair_v2::HeavyCoefficientFamily family)
{
    return family == wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ?
        "periodic-cauchy" : "hashed-nonzero";
}

struct CounterSummary
{
    uint64_t Calls = 0u;
    uint64_t AttemptsTotal = 0u;
    uint64_t CompletionsTotal = 0u;
    uint64_t FallbacksTotal = 0u;
    uint64_t AttemptsMin = UINT64_MAX;
    uint64_t AttemptsMax = 0u;
    uint64_t CompletionsMin = UINT64_MAX;
    uint64_t CompletionsMax = 0u;
    uint64_t FallbacksMin = UINT64_MAX;
    uint64_t FallbacksMax = 0u;
    bool IdentityOk = true;

    bool Observe(
        uint64_t attempts,
        uint64_t completions,
        uint64_t fallbacks,
        bool expect_direct,
        bool operation_executed,
        bool allow_fallback,
        uint64_t maximum_direct_attempts)
    {
        ++Calls;
        AttemptsTotal += attempts;
        CompletionsTotal += completions;
        FallbacksTotal += fallbacks;
        AttemptsMin = std::min(AttemptsMin, attempts);
        AttemptsMax = std::max(AttemptsMax, attempts);
        CompletionsMin = std::min(CompletionsMin, completions);
        CompletionsMax = std::max(CompletionsMax, completions);
        FallbacksMin = std::min(FallbacksMin, fallbacks);
        FallbacksMax = std::max(FallbacksMax, fallbacks);

        bool valid = attempts == 0u && completions == 0u && fallbacks == 0u;
        if (expect_direct && operation_executed)
        {
            // Exact-construction error classification may run a second
            // zero-RHS solve.  Receive may start a fresh solve for each
            // prefix through K+4 until it obtains resumable state.  In every
            // case each eligible solver entry must be accounted for as one
            // completed direct result or an explicitly permitted fallback.
            valid = attempts >= 1u &&
                attempts <= maximum_direct_attempts &&
                (allow_fallback ?
                    completions + fallbacks == attempts :
                    completions == attempts && fallbacks == 0u);
        }
        IdentityOk = IdentityOk && valid;
        return valid;
    }
};

struct ArmDefinition
{
    ArmMode Mode = ArmMode::DirectOff;
    NativeArmSpec Spec;
    std::string EquationSha256;
    std::string ConstructionSha256;
    std::string Identity;
    CounterSummary ConstructionCounters;
    CounterSummary OperationCounters;
};

std::string EquationsJson(
    const Request& request,
    ArmMode mode,
    const ResolvedNativeWh2Configuration* configuration)
{
    std::ostringstream json;
    if (mode == ArmMode::Wirehair1)
    {
        json << "{\"K\":" << request.K
             << ",\"block_bytes\":" << request.BlockBytes
             << ",\"codec\":\"wirehair1\",\"schema\":\""
             << kEquationSchema << "\"}";
        return json.str();
    }

    const wirehair_v2::PrecodeParams& params = configuration->Params;
    const wirehair_v2::PacketRowConfig& packet =
        configuration->PacketConfig;
    json << "{\"K\":" << request.K
         << ",\"block_bytes\":" << request.BlockBytes
         << ",\"codec\":\"wirehair2\",\"dense_anchors\":\""
         << DenseAnchorName(params.DenseAnchors)
         << "\",\"dense_identity_corner\":"
         << (params.DenseIdentityCorner ? "true" : "false")
         << ",\"dense_rows\":" << params.DenseRows
         << ",\"heavy_family\":\"" << HeavyFamilyName(params.HeavyFamily)
         << "\",\"heavy_rows\":" << params.HeavyRows
         << ",\"mix_count\":" << packet.MixCount
         << ",\"packet_attempt\":" << configuration->PacketAttempt
         << ",\"peel_seed\":" << packet.PeelSeed
         << ",\"precode_attempt\":" << configuration->PrecodeAttempt
         << ",\"schema\":\"" << kEquationSchema
         << "\",\"seed\":\"" << Hex64(params.Seed)
         << "\",\"source_hits\":" << params.SourceHits
         << ",\"staircase\":" << params.Staircase << "}";
    return json.str();
}

bool BuildArmDefinition(
    const Request& request,
    ArmMode mode,
    ArmDefinition& definition)
{
    definition.Mode = mode;
    ResolvedNativeWh2Configuration configuration;
    const ResolvedNativeWh2Configuration* configuration_ptr = nullptr;
    if (mode == ArmMode::Wirehair1) {
        definition.Spec = wirehair_wh2_bench::MakeWirehair1Arm();
    }
    else
    {
        definition.Spec = request.LayoutValue == Layout::Disabled ?
            wirehair_wh2_bench::MakeCertifiedWh2Arm(
                request.ConstructionAttempt) :
            wirehair_wh2_bench::MakeExperimentalWh2Arm(
                request.ConstructionAttempt, TwoAnchorTransform);
        const NativeWh2ExecutionMode execution_mode =
            mode == ArmMode::DirectOn ?
                NativeWh2ExecutionMode::TinyDirectEnabled :
                NativeWh2ExecutionMode::TinyDirectDisabled;
        if (!wirehair_wh2_bench::ConfigureNativeWh2ExecutionMode(
                definition.Spec,
                execution_mode,
                wirehair_wh2_bench::SetNativeTinyDirectSolveMode) ||
            !wirehair_wh2_bench::ResolveNativeWh2Configuration(
                definition.Spec,
                request.K,
                request.BlockBytes,
                configuration))
        {
            return false;
        }
        configuration_ptr = &configuration;
    }

    const std::string equations = EquationsJson(
        request, mode, configuration_ptr);
    definition.EquationSha256 =
        wirehair::wh2_benchmark::Sha256Hex(equations);
    if (definition.EquationSha256.empty()) {
        return false;
    }
    std::ostringstream construction;
    construction << "{\"K\":" << request.K
                 << ",\"block_bytes\":" << request.BlockBytes
                 << ",\"construction_attempt\":";
    if (mode == ArmMode::Wirehair1) {
        construction << "null";
    }
    else {
        construction << request.ConstructionAttempt;
    }
    construction << ",\"equations_sha256\":\""
                 << definition.EquationSha256 << "\",\"schema\":\""
                 << kConstructionSchema << "\"}";
    definition.ConstructionSha256 =
        wirehair::wh2_benchmark::Sha256Hex(construction.str());
    if (definition.ConstructionSha256.empty()) {
        return false;
    }

    std::string identity = ArmModeName(mode);
    identity += ':';
    identity += definition.ConstructionSha256;
    identity += ':';
    identity += ScopeName(request.ScopeValue);
    definition.Identity = wirehair_wh2_bench::BindNativeWh2ExecutionIdentity(
        identity, definition.Spec.ExecutionMode);
    return !definition.Identity.empty();
}

void ResetDirectCounters()
{
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
}

void ReadDirectCounters(
    uint64_t& attempts,
    uint64_t& completions,
    uint64_t& fallbacks)
{
    attempts = wirehair_v2::TinyDirectSolveAttemptsForTesting();
    completions = wirehair_v2::TinyDirectSolveCompletionsForTesting();
    fallbacks = wirehair_v2::TinyDirectSolveFallbacksForTesting();
}

NativePanelInvocationResult ClassifyTimedResult(
    WirehairResult result,
    uint64_t elapsed_ns,
    bool bytes_verified,
    bool has_decoded_overhead,
    uint32_t decoded_overhead)
{
    if (result == Wirehair_Success)
    {
        if (!bytes_verified || elapsed_ns == 0u ||
            elapsed_ns > static_cast<uint64_t>(
                std::numeric_limits<int64_t>::max()))
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                static_cast<int64_t>(Wirehair_Error), false, 0u, 0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            static_cast<int64_t>(result),
            has_decoded_overhead,
            decoded_overhead,
            elapsed_ns);
    }
    switch (result)
    {
    case Wirehair_NeedMore:
    case Wirehair_BadDenseSeed:
    case Wirehair_BadPeelSeed:
    case Wirehair_BadInput_SmallN:
    case Wirehair_BadInput_LargeN:
    case Wirehair_UnsupportedPlatform:
        return NativePanelInvocationResult(
            NativePanelDisposition::PreflightFailure,
            static_cast<int64_t>(result), false, 0u, 0u);
    default:
        return NativePanelInvocationResult(
            NativePanelDisposition::Fatal,
            static_cast<int64_t>(result), false, 0u, 0u);
    }
}

class CellInvocation : public NativePanelInvocation
{
public:
    CellInvocation(
        const Request& request,
        ArmDefinition& definition,
        const std::vector<uint8_t>& source,
        const std::vector<uint32_t>& packet_ids)
        : RequestValue(request)
        , Definition(definition)
    {
        ResetDirectCounters();
        if (request.ScopeValue == Scope::Encoder)
        {
            Encoder.reset(new NativeEncoderFixture);
            PreparationResult = Encoder->Initialize(
                definition.Spec,
                request.K,
                request.BlockBytes,
                source);
        }
        else
        {
            NativeArm arm;
            PreparationResult = arm.Initialize(
                definition.Spec,
                request.K,
                request.BlockBytes,
                source);
            if (PreparationResult == Wirehair_Success)
            {
                if (request.ScopeValue == Scope::Solve)
                {
                    Solve.reset(new NativeSolveFixture);
                    PreparationResult = Solve->Initialize(
                        arm, packet_ids, kSolveOverhead);
                }
                else
                {
                    Receive.reset(new NativeReceiveFixture);
                    PreparationResult = Receive->Initialize(
                        arm, packet_ids, kTraceOverhead);
                }
            }
        }

        uint64_t attempts = 0u;
        uint64_t completions = 0u;
        uint64_t fallbacks = 0u;
        ReadDirectCounters(attempts, completions, fallbacks);
        // A bad exact construction may be rejected while building the
        // precode system, before either solver is entered.  Preserve that
        // ordinary weak-attempt outcome.  Once direct counters show entry (or
        // construction succeeds), the enabled path must have exact identity.
        const bool codec_construction_executed =
            request.ScopeValue != Scope::Encoder &&
            (definition.Mode != ArmMode::DirectOn ||
             PreparationResult == Wirehair_Success || attempts != 0u);
        const bool identity_ok = Definition.ConstructionCounters.Observe(
            attempts,
            completions,
            fallbacks,
            definition.Mode == ArmMode::DirectOn,
            codec_construction_executed,
            false,
            2u);
        if (!identity_ok) {
            PreparationResult = Wirehair_Error;
        }
    }

    std::string Identity() const override
    {
        return Definition.Identity;
    }

    NativePanelInvocationResult Invoke() override
    {
        ResetDirectCounters();
        WirehairResult result = PreparationResult;
        uint64_t elapsed_ns = 0u;
        bool bytes_verified = false;
        bool has_decoded_overhead = false;
        uint32_t decoded_overhead = 0u;
        const bool operation_executed =
            PreparationResult == Wirehair_Success;
        if (operation_executed && RequestValue.ScopeValue == Scope::Solve)
        {
            const IsolatedSolveResult solve = Solve->Run();
            result = solve.Result;
            elapsed_ns = solve.ElapsedNanoseconds;
            bytes_verified = solve.BytesVerified;
        }
        else if (operation_executed &&
            RequestValue.ScopeValue == Scope::Receive)
        {
            const TimedArmResult receive = Receive->Run();
            result = receive.Result;
            elapsed_ns = receive.ElapsedNanoseconds;
            bytes_verified = receive.BytesVerified;
            has_decoded_overhead = receive.Result == Wirehair_Success;
            decoded_overhead = receive.DecodedOverhead;
        }
        else if (operation_executed)
        {
            const TimedArmResult encoder = Encoder->Run();
            result = encoder.Result;
            elapsed_ns = encoder.ElapsedNanoseconds;
            bytes_verified = encoder.BytesVerified;
        }

        uint64_t attempts = 0u;
        uint64_t completions = 0u;
        uint64_t fallbacks = 0u;
        ReadDirectCounters(attempts, completions, fallbacks);
        const bool identity_ok = Definition.OperationCounters.Observe(
            attempts,
            completions,
            fallbacks,
            Definition.Mode == ArmMode::DirectOn,
            operation_executed,
            RequestValue.ScopeValue == Scope::Receive,
            RequestValue.ScopeValue == Scope::Receive ?
                static_cast<uint64_t>(kSolveOverhead) + 1u :
                (RequestValue.ScopeValue == Scope::Encoder ? 2u : 1u));
        if (!identity_ok)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                static_cast<int64_t>(Wirehair_Error), false, 0u, 0u);
        }
        if (!operation_executed) {
            bytes_verified = false;
        }
        return ClassifyTimedResult(
            result,
            elapsed_ns,
            bytes_verified,
            has_decoded_overhead,
            decoded_overhead);
    }

private:
    Request RequestValue;
    ArmDefinition& Definition;
    WirehairResult PreparationResult = Wirehair_Error;
    std::unique_ptr<NativeEncoderFixture> Encoder;
    std::unique_ptr<NativeReceiveFixture> Receive;
    std::unique_ptr<NativeSolveFixture> Solve;
};

const char* DispositionName(NativePanelDisposition disposition)
{
    switch (disposition)
    {
    case NativePanelDisposition::Success: return "success";
    case NativePanelDisposition::PreflightFailure: return "preflight_failure";
    case NativePanelDisposition::Fatal: return "fatal";
    default: return "unknown";
    }
}

const char* OutcomeName(int64_t result)
{
    switch (static_cast<WirehairResult>(result))
    {
    case Wirehair_Success: return "success";
    case Wirehair_NeedMore: return "need_more";
    case Wirehair_InvalidInput: return "invalid_input";
    case Wirehair_BadDenseSeed: return "bad_dense_seed";
    case Wirehair_BadPeelSeed: return "bad_peel_seed";
    case Wirehair_BadInput_SmallN: return "small_n";
    case Wirehair_BadInput_LargeN: return "large_n";
    case Wirehair_ExtraInsufficient: return "extra_insufficient";
    case Wirehair_Error: return "error";
    case Wirehair_OOM: return "oom";
    case Wirehair_UnsupportedPlatform: return "unsupported_platform";
    default: return "unknown";
    }
}

std::string OutcomeJson(const NativePanelInvocationResult& result)
{
    std::ostringstream json;
    json << "{\"decoded_overhead\":";
    if (result.HasDecodedExtra) {
        json << result.DecodedExtra;
    }
    else {
        json << "null";
    }
    json << ",\"disposition\":\"" << DispositionName(result.Disposition)
         << "\",\"elapsed_ns\":";
    if (result.ElapsedNanoseconds != 0u) {
        json << result.ElapsedNanoseconds;
    }
    else {
        json << "null";
    }
    json << ",\"outcome\":\"" << OutcomeName(result.OutcomeCode)
         << "\",\"result_code\":" << result.OutcomeCode << "}";
    return json.str();
}

std::string CounterSummaryJson(const CounterSummary& summary)
{
    std::ostringstream json;
    json << "{\"attempts_max\":";
    if (summary.Calls) json << summary.AttemptsMax; else json << "null";
    json << ",\"attempts_min\":";
    if (summary.Calls) json << summary.AttemptsMin; else json << "null";
    json << ",\"attempts_total\":" << summary.AttemptsTotal
         << ",\"calls\":" << summary.Calls
         << ",\"completions_max\":";
    if (summary.Calls) json << summary.CompletionsMax; else json << "null";
    json << ",\"completions_min\":";
    if (summary.Calls) json << summary.CompletionsMin; else json << "null";
    json << ",\"completions_total\":" << summary.CompletionsTotal
         << ",\"fallbacks_max\":";
    if (summary.Calls) json << summary.FallbacksMax; else json << "null";
    json << ",\"fallbacks_min\":";
    if (summary.Calls) json << summary.FallbacksMin; else json << "null";
    json << ",\"fallbacks_total\":" << summary.FallbacksTotal
         << ",\"identity_ok\":"
         << (summary.IdentityOk ? "true" : "false") << "}";
    return json.str();
}

std::string ArmJson(const ArmDefinition& arm)
{
    std::ostringstream json;
    json << "{\"construction_sha256\":\"" << arm.ConstructionSha256
         << "\",\"direct\":{\"construction\":"
         << CounterSummaryJson(arm.ConstructionCounters)
         << ",\"operation\":"
         << CounterSummaryJson(arm.OperationCounters)
         << "},\"equations_sha256\":\"" << arm.EquationSha256
         << "\",\"identity\":" << JsonString(arm.Identity)
         << ",\"mode\":\"" << ArmModeName(arm.Mode) << "\"}";
    return json.str();
}

std::string NullableOutcomeJson(
    bool present,
    const NativePanelInvocationResult& result)
{
    return present ? OutcomeJson(result) : "null";
}

std::string PanelJson(const NativePanelResult& panel)
{
    std::ostringstream json;
    json << "{\"comparable\":"
         << (panel.PanelComparable ? "true" : "false")
         << ",\"diagnostic\":" << JsonString(panel.Diagnostic)
         << ",\"left_preflight\":"
         << NullableOutcomeJson(panel.HasLeftPreflight, panel.LeftPreflight)
         << ",\"order\":\"" << OrderName(panel.Order)
         << "\",\"right_preflight\":"
         << NullableOutcomeJson(panel.HasRightPreflight, panel.RightPreflight)
         << ",\"slots\":[";
    const bool slots_present = panel.Status == NativePanelStatus::Complete;
    for (size_t i = 0u; i < panel.Slots.size(); ++i)
    {
        if (i) json << ',';
        const wirehair_wh2_bench::NativePanelSlot& slot = panel.Slots[i];
        json << "{\"decoded_overhead\":";
        if (slots_present && slot.Invocation.HasDecodedExtra) {
            json << slot.Invocation.DecodedExtra;
        }
        else {
            json << "null";
        }
        json << ",\"disposition\":";
        if (slots_present) {
            json << JsonString(DispositionName(slot.Invocation.Disposition));
        }
        else {
            json << "null";
        }
        json << ",\"elapsed_ns\":";
        if (slots_present && slot.HasElapsedNanoseconds) {
            json << slot.ElapsedNanoseconds;
        }
        else {
            json << "null";
        }
        json << ",\"index\":" << i << ",\"outcome\":";
        if (slots_present) {
            json << JsonString(OutcomeName(slot.Invocation.OutcomeCode));
        }
        else {
            json << "null";
        }
        json << ",\"result_code\":";
        if (slots_present) {
            json << slot.Invocation.OutcomeCode;
        }
        else {
            json << "null";
        }
        json << ",\"side\":\""
             << (slot.Side == NativePanelSide::Left ? "left" : "right")
             << "\"}";
    }
    json << "],\"status\":\"" << NativePanelStatusName(panel.Status)
         << "\"}";
    return json.str();
}

std::string RecordJson(
    const Request& request,
    const std::string& request_json,
    const std::string& request_sha256,
    const std::string& source_sha256,
    uint64_t loss_seed,
    const FrozenPacketTrace& trace,
    const ArmDefinition& left,
    const ArmDefinition& right,
    const NativePanelResult& panel)
{
    std::ostringstream json;
    json << "{\"arms\":{\"left\":" << ArmJson(left)
         << ",\"right\":" << ArmJson(right)
         << "},\"panel\":" << PanelJson(panel)
         << ",\"request\":" << request_json
         << ",\"request_sha256\":\"" << request_sha256
         << "\",\"schema\":\"" << kSchema
         << "\",\"source_sha256\":\"" << source_sha256
         << "\",\"trace\":{\"attempted_candidates\":"
         << trace.attempted_candidates
         << ",\"candidate_limit\":" << trace.candidate_limit
         << ",\"derived_seed\":\"" << Hex64(trace.trace_seed)
         << "\",\"loss_seed\":\"" << Hex64(loss_seed)
         << "\",\"packet_count\":" << trace.delivered_ids.size()
         << ",\"requested_seed\":\"" << Hex64(request.TraceSeed)
         << "\",\"sha256\":\"" << trace.trace_sha256 << "\"}}";
    return json.str();
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--help")
    {
        PrintUsage(std::cout);
        return 0;
    }

    Request request;
    std::string error;
    if (!ParseRequest(argc, argv, request, error))
    {
        std::cerr << "wirehair_wh2_tiny_native_cell: " << error << '\n';
        PrintUsage(std::cerr);
        return 2;
    }

    const std::string request_json = CanonicalRequestJson(request);
    const std::string request_sha256 =
        wirehair::wh2_benchmark::Sha256Hex(request_json);
    // The controller freezes one realized loss seed per cell.  Consume that
    // seed exactly: Replicate is bound metadata and controls panel order, but
    // must not qualify an already-qualified seed a second time.
    const uint64_t loss_seed = request.TraceSeed;
    FrozenPacketTrace trace;
    if (request_sha256.empty() ||
        wirehair::wh2_benchmark::GenerateFrozenPacketTrace(
            request.K,
            request.BlockBytes,
            kLossPpm,
            loss_seed,
            FrozenSchedule::Iid,
            kTraceOverhead,
            trace) != FrozenTraceStatus::Complete ||
        trace.delivered_ids.size() !=
            static_cast<size_t>(request.K) + kTraceOverhead)
    {
        std::cerr << "wirehair_wh2_tiny_native_cell: trace generation failed\n";
        return 2;
    }

    std::vector<uint8_t> source;
    const uint64_t source_seed = loss_seed ^
        UINT64_C(0x6a09e667f3bcc909);
    if (!wirehair_wh2_bench::MakeDeterministicSource(
            request.K,
            request.BlockBytes,
            source_seed,
            source))
    {
        std::cerr << "wirehair_wh2_tiny_native_cell: source generation failed\n";
        return 2;
    }
    const std::string source_sha256 = wirehair::wh2_benchmark::Sha256Hex(
        source.empty() ? nullptr : source.data(), source.size());

    ArmDefinition left;
    ArmDefinition right;
    if (source_sha256.empty() ||
        !BuildArmDefinition(request, request.LeftMode, left) ||
        !BuildArmDefinition(request, request.RightMode, right))
    {
        std::cerr << "wirehair_wh2_tiny_native_cell: arm definition failed\n";
        return 2;
    }

    const NativePanelArm left_arm(
        left.Identity,
        [&request, &left, &source, &trace]() {
            return std::unique_ptr<NativePanelInvocation>(
                new CellInvocation(
                    request, left, source, trace.delivered_ids));
        });
    const NativePanelArm right_arm(
        right.Identity,
        [&request, &right, &source, &trace]() {
            return std::unique_ptr<NativePanelInvocation>(
                new CellInvocation(
                    request, right, source, trace.delivered_ids));
        });
    const NativePanelOrder order = (request.Replicate & 1u) ?
        NativePanelOrder::BAAB : NativePanelOrder::ABBA;
    const NativePanelResult panel =
        wirehair_wh2_bench::ExecuteNativeTimingPanel(
            request.Cpu,
            order,
            request.InvocationsPerSlot,
            left_arm,
            right_arm);

    std::cout << RecordJson(
        request,
        request_json,
        request_sha256,
        source_sha256,
        loss_seed,
        trace,
        left,
        right,
        panel) << '\n';
    return panel.Status == NativePanelStatus::Complete ? 0 : 1;
}
