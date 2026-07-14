#include "WirehairV2Codec.h"
#include "WirehairV2Plan.h"
#include "WirehairV2Precode.h"
#include "WirehairV2Seeds.h"
#include "WirehairV2Solve.h"

#include "../WirehairTools.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <new>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
static const uint32_t kMaxSeedTableTrials = 1000000u;

double NowSeconds()
{
    return std::chrono::duration<double>(
        Clock::now().time_since_epoch()).count();
}

bool ParseU32Scalar(const char* text, uint32_t& out)
{
    if (!text || !*text || *text < '0' || *text > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long value = std::strtoul(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value > UINT32_MAX) {
        return false;
    }
    out = (uint32_t)value;
    return true;
}

bool ParseU16Scalar(const char* text, uint16_t& out)
{
    uint32_t value = 0u;
    if (!ParseU32Scalar(text, value) || value > UINT16_MAX) {
        return false;
    }
    out = (uint16_t)value;
    return true;
}

bool ParseU64Scalar(const char* text, uint64_t& out)
{
    if (!text || !*text || *text < '0' || *text > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 0);
    if (errno != 0 || !end || *end != '\0') {
        return false;
    }
    out = (uint64_t)value;
    return true;
}

bool ParseDoubleScalar(const char* text, double& out)
{
    if (!text || !*text || ((*text < '0' || *text > '9') && *text != '.')) {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (errno != 0 || !end || *end != '\0' || !std::isfinite(value)) {
        return false;
    }
    out = value;
    return true;
}

bool BadArg(const char* option, const char* value)
{
    std::fprintf(stderr, "bad %s value: %s\n", option, value ? value : "");
    return false;
}

bool ParseU32Arg(const char* option, const char* value, uint32_t& out)
{
    return ParseU32Scalar(value, out) || BadArg(option, value);
}

bool ParseU16Arg(const char* option, const char* value, uint16_t& out)
{
    return ParseU16Scalar(value, out) || BadArg(option, value);
}

bool ParseU64Arg(const char* option, const char* value, uint64_t& out)
{
    return ParseU64Scalar(value, out) || BadArg(option, value);
}

bool ParseDoubleArg(const char* option, const char* value, double& out)
{
    return ParseDoubleScalar(value, out) || BadArg(option, value);
}

bool TakeArg(
    const char* command,
    const char* option,
    int argc,
    char** argv,
    int& i,
    const char*& value)
{
    if (i + 1 >= argc) {
        std::fprintf(stderr, "%s: %s requires a value\n", command, option);
        return false;
    }
    value = argv[++i];
    return true;
}

bool UnknownArg(const char* command, const char* option)
{
    std::fprintf(stderr, "%s: unknown option %s\n", command, option);
    return false;
}

struct Rng
{
    uint64_t State;

    explicit Rng(uint64_t seed) : State(seed) {}

    uint64_t Next()
    {
        uint64_t z = (State += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }

    uint32_t U32()
    {
        return (uint32_t)(Next() >> 32);
    }

    double Unit()
    {
        return static_cast<double>(Next() >> 11) *
            (1.0 / 9007199254740992.0);
    }
};

bool ShouldDrop(Rng& rng, double loss_rate)
{
    return rng.Unit() < loss_rate;
}

enum class PacketScheduleKind
{
    Iid,
    Burst,
    Permutation,
    SystematicFirst,
    RepairOnly,
    Adversarial
};

const char* PacketScheduleName(PacketScheduleKind kind)
{
    switch (kind)
    {
    case PacketScheduleKind::Iid:             return "iid";
    case PacketScheduleKind::Burst:           return "burst";
    case PacketScheduleKind::Permutation:     return "permutation";
    case PacketScheduleKind::SystematicFirst: return "systematic-first";
    case PacketScheduleKind::RepairOnly:      return "repair-only";
    case PacketScheduleKind::Adversarial:     return "adversarial";
    }
    return "unknown";
}

bool ParsePacketSchedule(const char* text, PacketScheduleKind& kind)
{
    if (!std::strcmp(text, "iid")) {
        kind = PacketScheduleKind::Iid;
    }
    else if (!std::strcmp(text, "burst")) {
        kind = PacketScheduleKind::Burst;
    }
    else if (!std::strcmp(text, "permutation")) {
        kind = PacketScheduleKind::Permutation;
    }
    else if (!std::strcmp(text, "systematic-first")) {
        kind = PacketScheduleKind::SystematicFirst;
    }
    else if (!std::strcmp(text, "repair-only")) {
        kind = PacketScheduleKind::RepairOnly;
    }
    else if (!std::strcmp(text, "adversarial")) {
        kind = PacketScheduleKind::Adversarial;
    }
    else {
        return false;
    }
    return true;
}

Rng MakeCompareLossRng(uint64_t seed);

std::vector<uint32_t> BuildPacketSchedule(
    uint32_t N,
    uint32_t delivered_count,
    double loss_rate,
    uint64_t seed,
    PacketScheduleKind kind)
{
    std::vector<uint32_t> output;
    output.reserve(delivered_count);
    Rng rng = MakeCompareLossRng(seed);
    uint64_t candidate_index = 0u;
    uint32_t burst_remaining = 0u;
    std::vector<uint32_t> permutation;
    size_t permutation_index = 0u;
    uint32_t permutation_base = 0u;

    const auto next_candidate = [&]() -> uint32_t {
        if (kind == PacketScheduleKind::RepairOnly) {
            return N + (uint32_t)candidate_index++;
        }
        if (kind == PacketScheduleKind::Adversarial) {
            return UINT32_MAX - (uint32_t)(candidate_index++ * 2u);
        }
        if (kind == PacketScheduleKind::Permutation)
        {
            if (permutation_index >= permutation.size())
            {
                const uint32_t count = std::min<uint32_t>(N + 512u, 65536u);
                permutation.resize(count);
                for (uint32_t i = 0; i < count; ++i) {
                    permutation[i] = permutation_base + i;
                }
                for (uint32_t i = count; i > 1u; --i) {
                    std::swap(permutation[i - 1u],
                        permutation[rng.U32() % i]);
                }
                permutation_base += count;
                permutation_index = 0u;
            }
            return permutation[permutation_index++];
        }
        if (kind == PacketScheduleKind::SystematicFirst &&
            candidate_index < N)
        {
            if (permutation.empty())
            {
                permutation.resize(N);
                for (uint32_t i = 0; i < N; ++i) {
                    permutation[i] = i;
                }
                for (uint32_t i = N; i > 1u; --i) {
                    std::swap(permutation[i - 1u],
                        permutation[rng.U32() % i]);
                }
            }
            return permutation[(size_t)candidate_index++];
        }
        if (kind == PacketScheduleKind::SystematicFirst) {
            return (uint32_t)candidate_index++;
        }
        return (uint32_t)candidate_index++;
    };

    const uint64_t candidate_limit =
        (uint64_t)delivered_count * 256u + 65536u;
    uint64_t candidates = 0u;
    while (output.size() < delivered_count && candidates++ < candidate_limit)
    {
        const uint32_t id = next_candidate();
        bool drop = false;
        if (kind == PacketScheduleKind::Burst)
        {
            static const uint32_t kBurstLength = 8u;
            if (burst_remaining > 0u) {
                --burst_remaining;
                drop = true;
            }
            else
            {
                // An idle candidate starts an eight-packet drop burst with
                // probability p.  Each renewal cycle then has eight drops
                // and (1-p)/p delivered idle candidates, so its stationary
                // loss fraction is 8p/(1+7p).  Inverting that expression
                // keeps the requested loss rate exact and p < 1 for every
                // accepted loss_rate < 1, including the 0.99 CLI boundary.
                const double start_probability = loss_rate /
                    (kBurstLength -
                        (kBurstLength - 1u) * loss_rate);
                if (rng.Unit() < start_probability) {
                    burst_remaining = kBurstLength - 1u;
                    drop = true;
                }
            }
        }
        else {
            drop = ShouldDrop(rng, loss_rate);
        }
        if (!drop) {
            output.push_back(id);
        }
    }
    return output;
}

struct TrialResult
{
    bool Ok;
    WirehairResult TerminalResult;
    uint32_t Extra;
    double CreateSeconds;
    double EncodeSeconds;
    double DecodeSeconds;
    double RecoverSeconds;
    uint64_t EncodedBytes;
    uint64_t DecodedBytes;
    uint64_t RecoveredBytes;
};

struct Accum
{
    uint64_t Trials = 0;
    uint64_t Failures = 0;
    uint64_t Nsum = 0;
    uint64_t OverheadSum = 0;
    uint64_t OverheadSq = 0;
    uint32_t OverheadMax = 0;
    std::vector<uint32_t> Overheads;
    double CreateSeconds = 0.0;
    double EncodeSeconds = 0.0;
    double DecodeSeconds = 0.0;
    double RecoverSeconds = 0.0;
    uint64_t CreateBytes = 0;
    uint64_t EncodeBytes = 0;
    uint64_t DecodeBytes = 0;
    uint64_t RecoverBytes = 0;
};

enum CompareProfileMode
{
    CompareProfileBase,
    CompareProfileAuto
};

enum PrecodeProfileMode
{
    PrecodeProfileCertified,
    PrecodeProfileMixed,
    PrecodeProfileBoth
};

const char* PrecodeProfileModeName(PrecodeProfileMode mode)
{
    switch (mode)
    {
    case PrecodeProfileCertified: return "certified";
    case PrecodeProfileMixed:     return "mixed";
    case PrecodeProfileBoth:      return "both";
    }
    return "unknown";
}

bool PrecodeProfileIncludes(
    PrecodeProfileMode mode,
    wirehair_v2::CompletionField completion)
{
    if (mode == PrecodeProfileBoth) {
        return true;
    }
    return completion == wirehair_v2::CompletionField::MixedGF256GF16 ?
        mode == PrecodeProfileMixed : mode == PrecodeProfileCertified;
}

const char* MixedCoefficientGeometryName(
    wirehair_v2::MixedCoefficientGeometry geometry)
{
    return geometry == wirehair_v2::MixedCoefficientGeometry::SharedCauchyX ?
        "shared-x" : "frozen";
}

bool ParseMixedCoefficientGeometry(
    const char* text,
    wirehair_v2::MixedCoefficientGeometry& geometry)
{
    if (!std::strcmp(text, "frozen")) {
        geometry = wirehair_v2::MixedCoefficientGeometry::FrozenPowerX;
        return true;
    }
    if (!std::strcmp(text, "shared-x")) {
        geometry = wirehair_v2::MixedCoefficientGeometry::SharedCauchyX;
        return true;
    }
    return false;
}

struct CompareOptions
{
    CompareProfileMode ProfileMode = CompareProfileBase;
    uint16_t PeelCandidates = 16;
    uint16_t PeelTrials = 3;
    uint16_t AutoTrials = 8;
    uint64_t TuneSeed = UINT64_C(0x9a7e11a);
    uint64_t AutoSeed = UINT64_C(0xa570ca1);
    double AutoMinDelta = 0.10;
    bool DenseOverride = false;
    int DenseDelta = 0;
    uint16_t DenseCandidate = 0;
    bool LogAutoChoices = false;
};

struct CachedCompareProfile
{
    bool UseTuned = false;
    wirehair_v2::SeedProfile TunedProfile;
};

struct PeelCostCandidate
{
    std::string Name;
    bool UsePolicy;
    wirehair_v2::PeelingCodec Codec;
};

enum PrecodeModelKind
{
    PrecodeModelDense,
    PrecodeModelLdpc,
    PrecodeModelLdpcDense,
    PrecodeModelCodecPort,
    PrecodeModelCodecPortIdentityCorner
};

struct PrecodeModel
{
    std::string Name;
    PrecodeModelKind Kind;
};

struct PrecodeRecipe
{
    uint32_t Columns;
    uint32_t LdpcColumns;
    uint32_t DenseRows;
    uint32_t HeavyRows;
    double GenerationXors;
};

struct PeelCostAccum
{
    uint64_t Trials = 0;
    double ResidualColumns = 0.0;
    double ResidualRows = 0.0;
    uint32_t ResidualColumnsMax = 0;
    double MatrixRefs = 0.0;
    double MatrixXors = 0.0;
    double LegacyTotalXors = 0.0;
    double SolveWidth = 0.0;
    double PrecodeGenXors = 0.0;
    double BackSubXors = 0.0;
    double GeBlockXors = 0.0;
    double HeavyMuladds = 0.0;
    double TotalBlockXors = 0.0;
};

std::vector<int> ParseIntList(const std::string& text)
{
    std::vector<int> out;
    size_t pos = 0;
    while (pos <= text.size())
    {
        const size_t comma = text.find(',', pos);
        const std::string token = text.substr(
            pos, comma == std::string::npos ? std::string::npos : comma - pos);
        if (token.empty() || token[0] < '0' || token[0] > '9') {
            return std::vector<int>();
        }
        errno = 0;
        char* end = nullptr;
        const long value = std::strtol(token.c_str(), &end, 10);
        if (errno != 0 || !end || *end != '\0' ||
            value < 1 || value > INT_MAX)
        {
            return std::vector<int>();
        }
        out.push_back((int)value);
        if (comma == std::string::npos) {
            break;
        }
        pos = comma + 1u;
    }
    return out;
}

std::vector<int> ParseSignedIntList(const std::string& text)
{
    std::vector<int> out;
    size_t pos = 0;
    while (pos <= text.size())
    {
        const size_t comma = text.find(',', pos);
        const std::string token = text.substr(
            pos, comma == std::string::npos ? std::string::npos : comma - pos);
        if (token.empty()) {
            return std::vector<int>();
        }
        if (token[0] == '-') {
            if (token.size() == 1u || token[1] < '0' || token[1] > '9') {
                return std::vector<int>();
            }
        }
        else if (token[0] < '0' || token[0] > '9') {
            return std::vector<int>();
        }
        errno = 0;
        char* end = nullptr;
        const long value = std::strtol(token.c_str(), &end, 10);
        if (errno != 0 || !end || *end != '\0' ||
            value < INT_MIN || value > INT_MAX)
        {
            return std::vector<int>();
        }
        out.push_back((int)value);
        if (comma == std::string::npos) {
            break;
        }
        pos = comma + 1u;
    }
    return out;
}

// Seed tables cover block counts [2, 64000].  GetDenseCount extrapolates
// past 64000 and GetDenseSeed then indexes the 100-entry kDenseSeeds table
// out of range, so reject out-of-domain N before building any seed profile.
bool ValidateBlockCounts(const std::vector<int>& Ns, const char* command)
{
    for (int n : Ns)
    {
        if (n < 2 || n > 64000) {
            std::fprintf(stderr,
                "%s --N values must be in [2,64000], got %d\n", command, n);
            return false;
        }
    }
    return true;
}

bool ValidateLoss(double loss, const char* command)
{
    if (!std::isfinite(loss) || loss < 0.0 || loss > 0.99)
    {
        std::fprintf(stderr,
            "%s --loss must be finite and in [0,0.99], got %.17g\n",
            command, loss);
        return false;
    }
    return true;
}

bool ValidateMessageDimensions(
    uint32_t N,
    uint32_t block_bytes,
    const char* command,
    uint64_t max_message_bytes = 0u,
    uint64_t working_block_count = 0u)
{
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    if (N < 2u || N > 64000u || block_bytes == 0u ||
        block_bytes > 0x7fffffffu || message_bytes > (uint64_t)SIZE_MAX)
    {
        std::fprintf(stderr,
            "%s message dimensions are unsupported: N=%u bb=%u\n",
            command, N, block_bytes);
        return false;
    }
    if (max_message_bytes > 0u && message_bytes > max_message_bytes)
    {
        std::fprintf(stderr,
            "%s message (%llu bytes) exceeds configured cap (%llu bytes)\n",
            command,
            (unsigned long long)message_bytes,
            (unsigned long long)max_message_bytes);
        return false;
    }
    const uint64_t metadata_bytes = (uint64_t)N * 4096u;
    if (working_block_count == 0u) {
        working_block_count = 2u * (uint64_t)N + 1u;
    }
    if (working_block_count >
        (UINT64_MAX - metadata_bytes) / block_bytes)
    {
        std::fprintf(stderr, "%s working-set size overflows\n", command);
        return false;
    }
    const uint64_t working_bytes =
        working_block_count * block_bytes + metadata_bytes;
    static const uint64_t kMaxBenchmarkWorkingBytes = UINT64_C(1) << 39;
    if (working_bytes > (uint64_t)SIZE_MAX ||
        working_bytes > kMaxBenchmarkWorkingBytes)
    {
        std::fprintf(stderr,
            "%s working set exceeds the supported benchmark limit\n",
            command);
        return false;
    }
#if defined(__unix__) || defined(__APPLE__)
    struct rlimit address_space_limit = {};
    if (getrlimit(RLIMIT_AS, &address_space_limit) == 0 &&
        address_space_limit.rlim_cur != RLIM_INFINITY &&
        working_bytes > (uint64_t)address_space_limit.rlim_cur)
    {
        std::fprintf(stderr,
            "%s working set cannot be allocated for N=%u bb=%u\n",
            command, N, block_bytes);
        return false;
    }
#endif
    try
    {
        // Reserve without touching pages so impossible/RLIMIT-constrained
        // workloads fail before any result header is emitted.
        std::vector<uint8_t> allocation_probe;
        allocation_probe.reserve((size_t)working_bytes);
    }
    catch (const std::bad_alloc&)
    {
        std::fprintf(stderr,
            "%s working set cannot be allocated for N=%u bb=%u\n",
            command, N, block_bytes);
        return false;
    }
    catch (const std::length_error&)
    {
        std::fprintf(stderr,
            "%s working set exceeds container limits for N=%u bb=%u\n",
            command, N, block_bytes);
        return false;
    }
    return true;
}

bool ValidateMessageInputs(
    const std::vector<int>& Ns,
    const std::vector<int>& BBs,
    const char* command)
{
    for (int n : Ns) for (int bb : BBs) {
        if (!ValidateMessageDimensions(
                (uint32_t)n, (uint32_t)bb, command))
        {
            return false;
        }
    }
    return true;
}

bool ValidatePayloadE2EInputs(
    const std::vector<int>& Ns,
    const std::vector<int>& BBs,
    const char* command,
    uint32_t concurrent_copies = 1u)
{
    for (int n : Ns) for (int bb : BBs)
    {
        if (n < 2 || n > 64000 || bb <= 0 ||
            (uint32_t)bb > UINT32_C(0x7fffffff))
        {
            return ValidateMessageDimensions(
                (uint32_t)n, (uint32_t)bb, command);
        }

        // Payload E2E keeps the source message, encoded intermediate values,
        // and delivered packet storage alive while a second solve allocates
        // its output and residual scratch.  Bound all block-sized storage:
        // 2*K covers message+delivery, 2*L covers the retained intermediate
        // and active solve values, 2*R covers the main and deficient-quotient
        // pivot RHS buffers, and the fixed margin covers successful mixed
        // quotient RHS, heavy buckets, and block temporaries.
        const wirehair_v2::SeedProfile profile =
            wirehair_v2::SelectSeedProfile((uint32_t)n, (uint32_t)bb);
        const uint64_t L = (uint64_t)n + profile.DenseCount + 24u;
        const uint64_t max_inactive = std::min<uint64_t>(
            L, wirehair_v2::kMaxInactiveColumns);
        const uint64_t working_blocks =
            2u * (uint64_t)n + 2u * L + 2u * max_inactive + 96u;
        if (working_blocks > UINT64_MAX / concurrent_copies)
        {
            std::fprintf(stderr, "%s working-set size overflows\n", command);
            return false;
        }
        if (!ValidateMessageDimensions(
                (uint32_t)n,
                (uint32_t)bb,
                command,
                0u,
                working_blocks * concurrent_copies))
        {
            return false;
        }
    }
    return true;
}

bool ValidateDenseDeltas(const std::vector<int>& deltas, const char* command)
{
    for (int delta : deltas)
    {
        if ((delta % 4) != 0)
        {
            std::fprintf(stderr,
                "%s --deltas/--dense-delta must preserve D %% 4 == 2; "
                "delta %d is invalid\n",
                command,
                delta);
            return false;
        }
    }
    return true;
}

bool DenseCountForDelta(
    const wirehair_v2::SeedProfile& base,
    int delta,
    const char* command,
    uint32_t N,
    uint32_t block_bytes,
    uint16_t* dense_count_out)
{
    const int dense_count = (int)base.DenseCount + delta;
    // GetDenseSeed reads kDenseSeeds[dense_count / 4] with only 100 entries,
    // and production dense counts are always 2 mod 4.
    if (dense_count < 2 || dense_count > 398 || (dense_count % 4) != 2)
    {
        std::fprintf(stderr,
            "%s dense delta %d gives invalid dense count %d for N=%u bb=%u\n",
            command,
            delta,
            dense_count,
            N,
            block_bytes);
        return false;
    }
    if (dense_count_out) {
        *dense_count_out = (uint16_t)dense_count;
    }
    return true;
}

bool ValidateDenseCountsForInputs(
    const std::vector<int>& Ns,
    const std::vector<int>& BBs,
    const std::vector<int>& Deltas,
    const char* command)
{
    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile(
                (uint32_t)n_value,
                (uint32_t)bb_value);
        for (int delta : Deltas)
        {
            if (!DenseCountForDelta(
                    base,
                    delta,
                    command,
                    (uint32_t)n_value,
                    (uint32_t)bb_value,
                    nullptr))
            {
                return false;
            }
        }
    }
    return true;
}

std::vector<std::string> ParseStringList(const std::string& text)
{
    std::vector<std::string> out;
    size_t pos = 0;
    while (pos < text.size())
    {
        const size_t comma = text.find(',', pos);
        const std::string token = text.substr(
            pos, comma == std::string::npos ? std::string::npos : comma - pos);
        if (!token.empty()) {
            out.push_back(token);
        }
        if (comma == std::string::npos) {
            break;
        }
        pos = comma + 1u;
    }
    return out;
}

void FillMessage(std::vector<uint8_t>& message, uint64_t seed)
{
    Rng rng(seed);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)rng.U32();
    }
}

static const uint64_t kCompareLossTraceSalt = UINT64_C(0x10fade);

Rng MakeCompareLossRng(uint64_t seed)
{
    // BuildPacketSchedule uses this stream to construct one common delivered-
    // id prefix per trial, which every comparison arm then replays unchanged.
    // The trial helpers also retain it for their schedule-free IID fallback.
    return Rng(seed ^ kCompareLossTraceSalt);
}

TrialResult RunBaselineTrial(
    uint32_t N,
    uint32_t block_bytes,
    double loss_rate,
    uint64_t seed,
    const std::vector<uint32_t>* packet_schedule = nullptr)
{
    TrialResult tr = {};
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> decoded((size_t)message_bytes, 0);
    std::vector<uint8_t> block(block_bytes);
    FillMessage(message, seed);
    Rng rng = MakeCompareLossRng(seed);

    const double c0 = NowSeconds();
    WirehairCodec enc =
        wirehair_encoder_create(0, &message[0], message_bytes, block_bytes);
    WirehairCodec dec =
        wirehair_decoder_create(0, message_bytes, block_bytes);
    const double c2 = NowSeconds();
    tr.CreateSeconds = c2 - c0;
    if (!enc || !dec)
    {
        if (enc) {
            wirehair_free(enc);
        }
        if (dec) {
            wirehair_free(dec);
        }
        return tr;
    }

    uint32_t delivered = 0;
    uint32_t block_id = 0;
    size_t schedule_index = 0u;
    uint32_t write_bytes = 0;
    const uint32_t max_delivered = N * 2u + 512u;
    while (delivered < max_delivered)
    {
        uint32_t this_id = 0u;
        if (packet_schedule) {
            if (schedule_index >= packet_schedule->size()) {
                break;
            }
            this_id = (*packet_schedule)[schedule_index++];
        }
        else {
            const bool drop = ShouldDrop(rng, loss_rate);
            this_id = block_id++;
            if (drop) {
                continue;
            }
        }
        const double e0 = NowSeconds();
        const WirehairResult er =
            wirehair_encode(enc, this_id, &block[0], block_bytes, &write_bytes);
        tr.EncodeSeconds += NowSeconds() - e0;
        if (er != Wirehair_Success) {
            break;
        }
        tr.EncodedBytes += write_bytes;
        ++delivered;
        const double d0 = NowSeconds();
        const WirehairResult dr =
            wirehair_decode(dec, this_id, &block[0], write_bytes);
        tr.DecodeSeconds += NowSeconds() - d0;
        tr.DecodedBytes += write_bytes;
        if (dr == Wirehair_Success) {
            tr.Ok = true;
            break;
        }
        if (dr != Wirehair_NeedMore) {
            break;
        }
    }

    if (tr.Ok)
    {
        const double r0 = NowSeconds();
        const WirehairResult rr =
            wirehair_recover(dec, &decoded[0], message_bytes);
        tr.RecoverSeconds = NowSeconds() - r0;
        tr.RecoveredBytes = message_bytes;
        tr.Ok = rr == Wirehair_Success &&
            std::memcmp(&decoded[0], &message[0], (size_t)message_bytes) == 0;
        if (tr.Ok) {
            tr.Extra = delivered - N;
        }
    }

    wirehair_free(enc);
    wirehair_free(dec);
    return tr;
}

TrialResult RunV2Trial(
    uint32_t N,
    uint32_t block_bytes,
    double loss_rate,
    uint64_t seed,
    const wirehair_v2::SeedProfile* profile,
    const std::vector<uint32_t>* packet_schedule = nullptr)
{
    TrialResult tr = {};
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> decoded((size_t)message_bytes, 0);
    std::vector<uint8_t> block(block_bytes);
    FillMessage(message, seed);
    Rng rng = MakeCompareLossRng(seed);

    wirehair_v2::Codec enc;
    wirehair_v2::Codec dec;
    const double c0 = NowSeconds();
    WirehairResult result =
        enc.InitializeEncoder(&message[0], message_bytes, block_bytes, profile);
    if (result == Wirehair_Success) {
        result = dec.InitializeDecoder(message_bytes, block_bytes, profile);
    }
    const double c2 = NowSeconds();
    tr.CreateSeconds = c2 - c0;
    if (result != Wirehair_Success) {
        return tr;
    }

    uint32_t delivered = 0;
    uint32_t block_id = 0;
    size_t schedule_index = 0u;
    uint32_t write_bytes = 0;
    const uint32_t max_delivered = N * 2u + 512u;
    while (delivered < max_delivered)
    {
        uint32_t this_id = 0u;
        if (packet_schedule) {
            if (schedule_index >= packet_schedule->size()) {
                break;
            }
            this_id = (*packet_schedule)[schedule_index++];
        }
        else {
            const bool drop = ShouldDrop(rng, loss_rate);
            this_id = block_id++;
            if (drop) {
                continue;
            }
        }
        const double e0 = NowSeconds();
        const WirehairResult er =
            enc.Encode(this_id, &block[0], block_bytes, &write_bytes);
        tr.EncodeSeconds += NowSeconds() - e0;
        if (er != Wirehair_Success) {
            break;
        }
        tr.EncodedBytes += write_bytes;
        ++delivered;
        const double d0 = NowSeconds();
        const WirehairResult dr =
            dec.Decode(this_id, &block[0], write_bytes);
        tr.DecodeSeconds += NowSeconds() - d0;
        tr.DecodedBytes += write_bytes;
        if (dr == Wirehair_Success) {
            tr.Ok = true;
            break;
        }
        if (dr != Wirehair_NeedMore) {
            break;
        }
    }

    if (tr.Ok)
    {
        const double r0 = NowSeconds();
        const WirehairResult rr = dec.Recover(&decoded[0], message_bytes);
        tr.RecoverSeconds = NowSeconds() - r0;
        tr.RecoveredBytes = message_bytes;
        tr.Ok = rr == Wirehair_Success &&
            std::memcmp(&decoded[0], &message[0], (size_t)message_bytes) == 0;
        if (tr.Ok) {
            tr.Extra = delivered - N;
        }
    }
    return tr;
}

TrialResult RunV2PrecodeTrial(
    uint32_t N,
    uint32_t block_bytes,
    double loss_rate,
    uint64_t seed,
    const std::vector<uint32_t>* packet_schedule = nullptr,
    wirehair_v2::CompletionField completion =
        wirehair_v2::CompletionField::GF256,
    bool cache_encoder_source = false,
    bool cache_decoder_systematic = false)
{
    TrialResult tr = {};
    tr.TerminalResult = Wirehair_Error;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> decoded((size_t)message_bytes, 0u);
    std::vector<uint8_t> block(block_bytes);
    FillMessage(message, seed);
    Rng rng = MakeCompareLossRng(seed);

    wirehair_v2::Codec enc;
    wirehair_v2::Codec dec;
    wirehair_v2::MessagePrecodeEncoderOptions encoder_options;
    encoder_options.Completion = completion;
    encoder_options.CacheSystematicSource = cache_encoder_source;
    encoder_options.CacheReceivedSystematicPackets = cache_decoder_systematic;
    const bool use_encoder_options =
        completion != wirehair_v2::CompletionField::GF256 ||
        cache_encoder_source;
    const double c0 = NowSeconds();
    WirehairResult result = enc.InitializePrecodeEncoder(
        message.data(), message_bytes, block_bytes, nullptr,
        use_encoder_options ? &encoder_options : nullptr);
    if (result == Wirehair_Success)
    {
        // The public wire protocol serializes the encoder-selected seed
        // attempt in its profile descriptor.  Reuse that exact descriptor
        // here instead of charging compare for an artificial second seed
        // search in the decoder.
        result = dec.InitializePrecodeDecoder(
            message_bytes, block_bytes, &enc.Profile(),
            cache_decoder_systematic ? &encoder_options : nullptr);
        if (result == Wirehair_Success &&
            (enc.Profile().V2CompletionField != completion ||
             dec.Profile().V2CompletionField != completion ||
             dec.Profile().V2SeedAttempt != enc.Profile().V2SeedAttempt ||
             dec.Profile().V2PacketPeelSeed !=
                enc.Profile().V2PacketPeelSeed ||
             dec.Profile().V2PrecodeSeed != enc.Profile().V2PrecodeSeed))
        {
            result = Wirehair_Error;
        }
    }
    tr.CreateSeconds = NowSeconds() - c0;
    if (result != Wirehair_Success) {
        tr.TerminalResult = result;
        return tr;
    }

    uint32_t delivered = 0u;
    uint32_t block_id = 0u;
    size_t schedule_index = 0u;
    uint32_t write_bytes = 0u;
    const uint32_t max_delivered = N * 2u + 512u;
    result = Wirehair_NeedMore;
    while (delivered < max_delivered)
    {
        uint32_t this_id = 0u;
        if (packet_schedule) {
            if (schedule_index >= packet_schedule->size()) {
                break;
            }
            this_id = (*packet_schedule)[schedule_index++];
        }
        else {
            const bool drop = ShouldDrop(rng, loss_rate);
            this_id = block_id++;
            if (drop) {
                continue;
            }
        }
        const double e0 = NowSeconds();
        result = enc.Encode(
            this_id, block.data(), block_bytes, &write_bytes);
        tr.EncodeSeconds += NowSeconds() - e0;
        if (result != Wirehair_Success) {
            break;
        }
        tr.EncodedBytes += write_bytes;
        ++delivered;

        const double d0 = NowSeconds();
        result = dec.Decode(this_id, block.data(), write_bytes);
        tr.DecodeSeconds += NowSeconds() - d0;
        tr.DecodedBytes += write_bytes;
        if (result == Wirehair_Success) {
            tr.Ok = true;
            break;
        }
        if (result != Wirehair_NeedMore) {
            break;
        }
    }

    if (tr.Ok)
    {
        const double r0 = NowSeconds();
        result = dec.Recover(decoded.data(), message_bytes);
        tr.RecoverSeconds = NowSeconds() - r0;
        tr.RecoveredBytes = message_bytes;
        tr.Ok = result == Wirehair_Success &&
            std::memcmp(decoded.data(), message.data(),
                (size_t)message_bytes) == 0;
        if (!tr.Ok && result == Wirehair_Success) {
            result = Wirehair_Error;
        }
        if (tr.Ok) {
            tr.Extra = delivered - N;
        }
    }
    tr.TerminalResult = result;
    return tr;
}

void AddTrial(Accum& acc, uint32_t N, const TrialResult& tr)
{
    ++acc.Trials;
    acc.Nsum += N;
    if (!tr.Ok)
    {
        ++acc.Failures;
        return;
    }
    acc.OverheadSum += tr.Extra;
    acc.OverheadSq += (uint64_t)tr.Extra * (uint64_t)tr.Extra;
    acc.Overheads.push_back(tr.Extra);
    if (tr.Extra > acc.OverheadMax) {
        acc.OverheadMax = tr.Extra;
    }
    acc.CreateSeconds += tr.CreateSeconds;
    acc.EncodeSeconds += tr.EncodeSeconds;
    acc.DecodeSeconds += tr.DecodeSeconds;
    acc.RecoverSeconds += tr.RecoverSeconds;
    acc.CreateBytes += tr.RecoveredBytes;
    acc.EncodeBytes += tr.EncodedBytes;
    acc.DecodeBytes += tr.DecodedBytes;
    acc.RecoverBytes += tr.RecoveredBytes;
}

double Mbps(uint64_t bytes, double seconds)
{
    return seconds > 0.0 ? (double)bytes / seconds / 1000000.0 : 0.0;
}

double MeanOverheadOrFailure(const Accum& acc)
{
    const uint64_t ok = acc.Trials - acc.Failures;
    return ok > 0u ? (double)acc.OverheadSum / (double)ok : -1.0;
}

double SelectionMeanOverhead(const Accum& acc)
{
    const uint64_t ok = acc.Trials - acc.Failures;
    return ok > 0u ? (double)acc.OverheadSum / (double)ok : 1e9;
}

bool IsBetterSelection(
    const Accum& candidate,
    const Accum& base,
    double min_mean_delta)
{
    if (candidate.Failures != base.Failures) {
        return candidate.Failures < base.Failures;
    }

    const double candidate_mean = SelectionMeanOverhead(candidate);
    const double base_mean = SelectionMeanOverhead(base);
    if (base_mean - candidate_mean >= min_mean_delta) {
        return true;
    }

    return false;
}

uint64_t ProfileCacheKey(uint32_t N, uint32_t block_bytes)
{
    return ((uint64_t)N << 32) | (uint64_t)block_bytes;
}

const char* CompareProfileModeName(CompareProfileMode mode)
{
    switch (mode) {
    case CompareProfileBase:
        return "base";
    case CompareProfileAuto:
        return "auto";
    }
    return "unknown";
}

uint64_t Mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

wirehair_v2::PeelingCodec MakeBenchLtCodec(
    uint16_t min_degree,
    uint16_t max_degree,
    wirehair_v2::PeelSolver solver)
{
    wirehair_v2::PeelingCodec codec;
    codec.Solver = solver;
    codec.Structure = wirehair_v2::PeelStructure::LtM1C32;
    codec.Family = wirehair_v2::DegreeFamily::Lt;
    codec.MinDegree = min_degree;
    codec.MaxDegree = max_degree;
    codec.SolverCandidateLimit =
        solver == wirehair_v2::PeelSolver::KsBmaxTop16 ? 16u : 0u;
    codec.Degree1Mass = 0.0;
    codec.Degree2Mass = 0.0;
    codec.RobustC = 0.0;
    codec.RobustDelta = 0.0;
    codec.FullyRandomRows = true;
    return codec;
}

wirehair_v2::PeelingCodec MakeBenchRobustD1D2Codec(
    uint16_t max_degree,
    double degree1_mass,
    double degree2_mass,
    wirehair_v2::PeelSolver solver)
{
    wirehair_v2::PeelingCodec codec;
    codec.Solver = solver;
    codec.Structure = wirehair_v2::PeelStructure::RobustD1_001D2_003;
    codec.Family = wirehair_v2::DegreeFamily::RobustD1D2;
    codec.MinDegree = 1u;
    codec.MaxDegree = max_degree;
    codec.SolverCandidateLimit =
        solver == wirehair_v2::PeelSolver::KsBmaxTop16 ? 16u : 0u;
    codec.Degree1Mass = degree1_mass;
    codec.Degree2Mass = degree2_mass;
    codec.RobustC = 0.0;
    codec.RobustDelta = 0.0;
    codec.FullyRandomRows = true;
    return codec;
}

wirehair_v2::PeelingCodec MakeBenchRobustSolitonCodec(
    double c,
    double delta,
    uint16_t max_degree,
    wirehair_v2::PeelSolver solver)
{
    wirehair_v2::PeelingCodec codec;
    codec.Solver = solver;
    codec.Structure = wirehair_v2::PeelStructure::RsC001D50C128;
    codec.Family = wirehair_v2::DegreeFamily::RobustSoliton;
    codec.MinDegree = 1u;
    codec.MaxDegree = max_degree;
    codec.SolverCandidateLimit =
        solver == wirehair_v2::PeelSolver::KsBmaxTop16 ? 16u : 0u;
    codec.Degree1Mass = 0.0;
    codec.Degree2Mass = 0.0;
    codec.RobustC = c;
    codec.RobustDelta = delta;
    codec.FullyRandomRows = true;
    return codec;
}

bool ParsePeelSolver(
    const std::string& name,
    wirehair_v2::PeelSolver& solver)
{
    if (name == "ks_bmax_top16" || name == "ks") {
        solver = wirehair_v2::PeelSolver::KsBmaxTop16;
        return true;
    }
    if (name == "rqcc_lowref" || name == "rqcc") {
        solver = wirehair_v2::PeelSolver::RqccLowref;
        return true;
    }
    return false;
}

bool MakeNamedPeelCostCandidate(
    const std::string& name,
    wirehair_v2::PeelSolver solver,
    PeelCostCandidate& candidate)
{
    candidate.Name = name;
    candidate.UsePolicy = false;
    if (name == "policy") {
        candidate.UsePolicy = true;
        candidate.Codec = MakeBenchLtCodec(1u, 32u, solver);
        return true;
    }
    if (name == "lt_m1_c16") {
        candidate.Codec = MakeBenchLtCodec(1u, 16u, solver);
        return true;
    }
    if (name == "lt_m1_c32") {
        candidate.Codec = MakeBenchLtCodec(1u, 32u, solver);
        return true;
    }
    if (name == "lt_m1_c64") {
        candidate.Codec = MakeBenchLtCodec(1u, 64u, solver);
        return true;
    }
    if (name == "lt_m2_c96") {
        candidate.Codec = MakeBenchLtCodec(2u, 96u, solver);
        return true;
    }
    if (name == "lt_m2_c128") {
        candidate.Codec = MakeBenchLtCodec(2u, 128u, solver);
        return true;
    }
    if (name == "lt_m2_c256") {
        candidate.Codec = MakeBenchLtCodec(2u, 256u, solver);
        return true;
    }
    if (name == "lt_m2_c512") {
        candidate.Codec = MakeBenchLtCodec(2u, 512u, solver);
        return true;
    }
    if (name == "lt_m2_c1024") {
        candidate.Codec = MakeBenchLtCodec(2u, 1024u, solver);
        return true;
    }
    if (name == "robust_d1_001_d2_003") {
        candidate.Codec = MakeBenchRobustD1D2Codec(64u, 0.01, 0.03, solver);
        return true;
    }
    if (name == "robust_d1_001_d2_012") {
        candidate.Codec = MakeBenchRobustD1D2Codec(64u, 0.01, 0.12, solver);
        return true;
    }
    if (name == "rs_c001_d50_c128") {
        candidate.Codec = MakeBenchRobustSolitonCodec(0.01, 0.50, 128u, solver);
        return true;
    }
    if (name == "rs_c003_d10_c128") {
        candidate.Codec = MakeBenchRobustSolitonCodec(0.03, 0.10, 128u, solver);
        return true;
    }
    return false;
}

bool MakePrecodeModel(const std::string& name, PrecodeModel& model)
{
    model.Name = name;
    if (name == "dense") {
        model.Kind = PrecodeModelDense;
        return true;
    }
    if (name == "ldpc") {
        model.Kind = PrecodeModelLdpc;
        return true;
    }
    if (name == "ldpcdense") {
        model.Kind = PrecodeModelLdpcDense;
        return true;
    }
    if (name == "codecport") {
        model.Kind = PrecodeModelCodecPort;
        return true;
    }
    if (name == "codecport_ic") {
        model.Kind = PrecodeModelCodecPortIdentityCorner;
        return true;
    }
    return false;
}

bool IsCodecPortPrecode(const PrecodeModel& model)
{
    return model.Kind == PrecodeModelCodecPort ||
        model.Kind == PrecodeModelCodecPortIdentityCorner;
}

double CertifiedPrecodeGenerationXors(
    const wirehair_v2::PrecodeParams& params)
{
    const uint32_t hits =
        params.SourceHits < params.Staircase ?
        params.SourceHits : params.Staircase;
    const uint64_t staircase_xors =
        (uint64_t)hits * params.BlockCount + params.Staircase - 1u;
    const uint64_t deck_span =
        (uint64_t)params.BlockCount + params.Staircase +
        (params.DenseIdentityCorner ? 0u : params.DenseRows);
    const uint64_t dense_xors =
        params.DenseRows > 0u ?
        ((deck_span + 1u) >> 1) + 2u * (params.DenseRows - 1u) : 0u;
    return (double)(staircase_xors + dense_xors);
}

PrecodeRecipe BuildPrecodeRecipe(
    const PrecodeModel& model,
    uint32_t block_count,
    uint32_t heavy_rows)
{
    const uint32_t dense_count = wirehair::GetDenseCount(block_count);
    PrecodeRecipe recipe;
    recipe.Columns = dense_count;
    recipe.LdpcColumns = 0;
    recipe.DenseRows = dense_count;
    recipe.HeavyRows = heavy_rows;
    recipe.GenerationXors = 2.5 * (double)(block_count + dense_count);

    if (model.Kind == PrecodeModelLdpc)
    {
        recipe.LdpcColumns = dense_count;
        recipe.DenseRows = 0;
        recipe.GenerationXors =
            3.0 * (double)block_count + (double)recipe.LdpcColumns;
    }
    else if (model.Kind == PrecodeModelLdpcDense)
    {
        if (block_count <= 1200u)
        {
            recipe.LdpcColumns = 25u;
            recipe.DenseRows = 16u;
        }
        else if (block_count <= 5000u)
        {
            recipe.LdpcColumns = 43u;
            recipe.DenseRows = 12u;
        }
        else if (block_count <= 20000u)
        {
            recipe.LdpcColumns = 52u;
            recipe.DenseRows = 34u;
        }
        else
        {
            recipe.LdpcColumns = dense_count;
            recipe.DenseRows = 0;
        }
        recipe.Columns = recipe.LdpcColumns + recipe.DenseRows;
        // Only charge Shuffle-2 dense generation when dense rows exist;
        // the pure-LDPC fallback above must cost the same as PrecodeModelLdpc
        recipe.GenerationXors =
            3.0 * (double)block_count +
            (double)recipe.LdpcColumns +
            (recipe.DenseRows > 0u ?
                2.5 * (double)(block_count + recipe.DenseRows) : 0.0);
    }
    else if (IsCodecPortPrecode(model))
    {
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(block_count, 0u);
        params.DenseIdentityCorner =
            model.Kind == PrecodeModelCodecPortIdentityCorner;
        recipe.Columns = params.Staircase + params.DenseRows;
        recipe.LdpcColumns = params.Staircase;
        recipe.DenseRows = params.DenseRows;
        recipe.HeavyRows = params.HeavyRows;
        recipe.GenerationXors = CertifiedPrecodeGenerationXors(params);
    }
    return recipe;
}

uint64_t PeelCostMatrixSeed(
    uint64_t seed,
    uint32_t block_count,
    uint32_t overhead_rows,
    uint32_t trial)
{
    uint64_t x = seed;
    x ^= Mix64((uint64_t)block_count * UINT64_C(0x9e3779b97f4a7c15));
    x ^= Mix64((uint64_t)overhead_rows * UINT64_C(0x94d049bb133111eb));
    x ^= Mix64((uint64_t)trial * UINT64_C(0xd6e8feb86659fd93));
    return Mix64(x);
}

void AddPeelCostTrial(
    PeelCostAccum& acc,
    const wirehair_v2::PeelEvaluation& eval,
    const PrecodeRecipe& recipe)
{
    const double solve_width =
        (double)eval.ResidualColumns +
        (double)recipe.Columns +
        (double)recipe.HeavyRows;
    const double backsub_xors = solve_width * (double)eval.Columns;
    const double ge_block_xors = solve_width * solve_width * 0.5;
    const double heavy_muladds = solve_width * (double)recipe.HeavyRows;
    const double total =
        (double)eval.MatrixXors +
        recipe.GenerationXors +
        backsub_xors +
        ge_block_xors +
        heavy_muladds;

    ++acc.Trials;
    acc.ResidualColumns += (double)eval.ResidualColumns;
    acc.ResidualRows += (double)eval.ResidualRows;
    if (eval.ResidualColumns > acc.ResidualColumnsMax) {
        acc.ResidualColumnsMax = eval.ResidualColumns;
    }
    acc.MatrixRefs += (double)eval.MatrixRefs;
    acc.MatrixXors += (double)eval.MatrixXors;
    acc.LegacyTotalXors += (double)eval.TotalXorCost;
    acc.SolveWidth += solve_width;
    acc.PrecodeGenXors += recipe.GenerationXors;
    acc.BackSubXors += backsub_xors;
    acc.GeBlockXors += ge_block_xors;
    acc.HeavyMuladds += heavy_muladds;
    acc.TotalBlockXors += total;
}

double CostMean(double sum, uint64_t trials)
{
    return trials > 0u ? sum / (double)trials : 0.0;
}

wirehair_v2::SeedProfile BuildTunedProfile(
    uint32_t N,
    uint32_t block_bytes,
    const CompareOptions& options)
{
    wirehair_v2::SeedTuningOptions tuning =
        wirehair_v2::DefaultSeedTuningOptions();
    tuning.PeelCandidates = options.PeelCandidates;
    tuning.TrialsPerCandidate = options.PeelTrials;
    tuning.Seed = options.TuneSeed;
    return wirehair_v2::TuneSeedProfile(N, block_bytes, tuning);
}

void ApplyDenseOverride(
    wirehair_v2::SeedProfile& profile,
    const CompareOptions& options)
{
    if (!options.DenseOverride) {
        return;
    }

    uint16_t dense_count = 0u;
    if (!DenseCountForDelta(
            profile,
            options.DenseDelta,
            "compare",
            profile.BlockCount,
            profile.BlockBytes,
            &dense_count))
    {
        return;
    }

    profile.DenseCount = dense_count;
    const uint16_t table_seed =
        wirehair::GetDenseSeed(profile.BlockCount, profile.DenseCount);
    profile.DenseSeed =
        wirehair_v2::CandidateDenseSeed(table_seed, options.DenseCandidate);
    profile.UsedDenseFixup = profile.DenseSeed != table_seed;
}

void CalibrateAutoProfile(
    uint32_t N,
    uint32_t block_bytes,
    double loss,
    const CompareOptions& options,
    CachedCompareProfile& cached)
{
    // Apply the dense override before calibrating so the A/B decision is
    // made on exactly the profiles that will be measured afterwards.
    cached.TunedProfile = BuildTunedProfile(N, block_bytes, options);
    ApplyDenseOverride(cached.TunedProfile, options);

    wirehair_v2::SeedProfile base_profile;
    const wirehair_v2::SeedProfile* base_arg = 0;
    if (options.DenseOverride)
    {
        base_profile = wirehair_v2::SelectSeedProfile(N, block_bytes);
        ApplyDenseOverride(base_profile, options);
        base_arg = &base_profile;
    }

    Accum base_acc;
    Accum tuned_acc;
    for (uint32_t trial = 0; trial < options.AutoTrials; ++trial)
    {
        const uint64_t trial_seed =
            options.AutoSeed ^
            ((uint64_t)N * UINT64_C(0x9e3779b97f4a7c15)) ^
            ((uint64_t)block_bytes * UINT64_C(0xbf58476d1ce4e5b9)) ^
            ((uint64_t)trial * UINT64_C(0x94d049bb133111eb));
        AddTrial(base_acc, N,
            RunV2Trial(N, block_bytes, loss, trial_seed, base_arg));
        AddTrial(tuned_acc, N,
            RunV2Trial(N, block_bytes, loss, trial_seed,
                &cached.TunedProfile));
    }

    cached.UseTuned =
        IsBetterSelection(tuned_acc, base_acc, options.AutoMinDelta);
    if (options.LogAutoChoices)
    {
        std::printf(
            "# auto_profile N=%u bb=%u choice=%s base_fail=%llu "
            "base_oh=%.4f base_max=%u tuned_fail=%llu tuned_oh=%.4f "
            "tuned_max=%u min_delta=%.4f\n",
            N,
            block_bytes,
            cached.UseTuned ? "tuned" : "base",
            (unsigned long long)base_acc.Failures,
            MeanOverheadOrFailure(base_acc),
            base_acc.OverheadMax,
            (unsigned long long)tuned_acc.Failures,
            MeanOverheadOrFailure(tuned_acc),
            tuned_acc.OverheadMax,
            options.AutoMinDelta);
    }
}

const wirehair_v2::SeedProfile* SelectCompareProfile(
    uint32_t N,
    uint32_t block_bytes,
    double loss,
    const CompareOptions& options,
    std::map<uint64_t, CachedCompareProfile>& cache)
{
    if (options.ProfileMode == CompareProfileBase) {
        if (options.DenseOverride)
        {
            CachedCompareProfile cached;
            cached.UseTuned = true;
            cached.TunedProfile =
                wirehair_v2::SelectSeedProfile(N, block_bytes);
            ApplyDenseOverride(cached.TunedProfile, options);
            const std::map<uint64_t, CachedCompareProfile>::iterator inserted =
                cache.insert(std::make_pair(
                    ProfileCacheKey(N, block_bytes), cached)).first;
            return &inserted->second.TunedProfile;
        }
        return 0;
    }

    const uint64_t key = ProfileCacheKey(N, block_bytes);
    const std::map<uint64_t, CachedCompareProfile>::iterator found =
        cache.find(key);
    if (found != cache.end()) {
        return found->second.UseTuned ? &found->second.TunedProfile : 0;
    }

    CachedCompareProfile cached;
    // CalibrateAutoProfile already applies the dense override to both
    // arms, so the winning profile needs no further mutation here.
    CalibrateAutoProfile(N, block_bytes, loss, options, cached);
    if (!cached.UseTuned && options.DenseOverride)
    {
        cached.UseTuned = true;
        cached.TunedProfile =
            wirehair_v2::SelectSeedProfile(N, block_bytes);
        ApplyDenseOverride(cached.TunedProfile, options);
    }
    const std::map<uint64_t, CachedCompareProfile>::iterator inserted =
        cache.insert(std::make_pair(key, cached)).first;
    return inserted->second.UseTuned ? &inserted->second.TunedProfile : 0;
}

uint32_t Percentile(std::vector<uint32_t> values, double p)
{
    if (values.empty()) {
        return 0;
    }
    std::sort(values.begin(), values.end());
    size_t index = (size_t)(p * (double)(values.size() - 1u) + 0.5);
    if (index >= values.size()) {
        index = values.size() - 1u;
    }
    return values[index];
}

double OverheadStdDev(const Accum& acc)
{
    const uint64_t ok = acc.Trials - acc.Failures;
    const double mean = MeanOverheadOrFailure(acc);
    if (ok == 0u || mean < 0.0) {
        return 0.0;
    }
    double var = (double)acc.OverheadSq / (double)ok - mean * mean;
    if (var < 0.0) {
        var = 0.0;
    }
    return std::sqrt(var);
}

void PrintAccum(const char* name, uint32_t block_bytes, const Accum& acc)
{
    const uint64_t ok = acc.Trials - acc.Failures;
    const double n_mean =
        acc.Trials > 0u ? (double)acc.Nsum / (double)acc.Trials : 0.0;
    const double overhead_mean =
        ok > 0u ? (double)acc.OverheadSum / (double)ok : 0.0;
    double overhead_var = 0.0;
    if (ok > 0u)
    {
        overhead_var = (double)acc.OverheadSq / (double)ok -
            overhead_mean * overhead_mean;
        if (overhead_var < 0.0) {
            overhead_var = 0.0;
        }
    }
    std::printf(
        "%-15s %-8u %-7llu %-7llu %-10.1f %-10.4f %-8.4f "
        "%-6u %-6u %-6u %-8u "
        "%12.1f %12.1f %12.1f %12.1f\n",
        name,
        block_bytes,
        (unsigned long long)acc.Trials,
        (unsigned long long)acc.Failures,
        n_mean,
        overhead_mean,
        std::sqrt(overhead_var),
        Percentile(acc.Overheads, 0.50),
        Percentile(acc.Overheads, 0.95),
        Percentile(acc.Overheads, 0.99),
        acc.OverheadMax,
        Mbps(acc.CreateBytes, acc.CreateSeconds),
        Mbps(acc.EncodeBytes, acc.EncodeSeconds),
        Mbps(acc.DecodeBytes, acc.DecodeSeconds),
        Mbps(acc.RecoverBytes, acc.RecoverSeconds));
}

int CmdCompare(int argc, char** argv)
{
    uint32_t nlo = 64u;
    uint32_t nhi = 2048u;
    uint32_t trials = 6u;
    uint32_t max_message_mib = 128u;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0xc0decafe);
    std::string bb_list = "1280,102400";
    bool include_precode = false;
    bool include_precode_cache = false;
    bool cache_encoder_source = false;
    bool cache_decoder_systematic = false;
    bool trial_details = false;
    PrecodeProfileMode precode_profile = PrecodeProfileCertified;
    bool precode_profile_explicit = false;
    PacketScheduleKind schedule_kind = PacketScheduleKind::Iid;
    CompareOptions compare_options;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    uint32_t mixed_period = wirehair_v2::kMixedCoefficientPeriod;
    bool mixed_period_explicit = false;
    uint32_t mixed_gf16_rows = wirehair_v2::kMixedGF16Rows;
    bool mixed_gf16_rows_explicit = false;
    wirehair_v2::MixedCoefficientGeometry mixed_geometry =
        wirehair_v2::MixedCoefficientGeometry::FrozenPowerX;
    bool mixed_geometry_explicit = false;
#endif

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--nlo")) {
            if (!TakeArg("compare", "--nlo", argc, argv, i, value) ||
                !ParseU32Arg("--nlo", value, nlo)) return 1;
        }
        else if (!std::strcmp(argv[i], "--nhi")) {
            if (!TakeArg("compare", "--nhi", argc, argv, i, value) ||
                !ParseU32Arg("--nhi", value, nhi)) return 1;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg("compare", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg("compare", "--bb-list", argc, argv, i, value)) return 1;
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!TakeArg("compare", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss)) return 1;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg("compare", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed)) return 1;
        }
        else if (!std::strcmp(argv[i], "--max-message-mib")) {
            if (!TakeArg("compare", "--max-message-mib", argc, argv, i, value) ||
                !ParseU32Arg("--max-message-mib", value, max_message_mib)) return 1;
        }
        else if (!std::strcmp(argv[i], "--precode")) {
            include_precode = true;
        }
        else if (!std::strcmp(argv[i], "--precode-cache")) {
            include_precode_cache = true;
            cache_encoder_source = true;
            cache_decoder_systematic = true;
        }
        else if (!std::strcmp(argv[i], "--precode-encoder-cache")) {
            include_precode_cache = true;
            cache_encoder_source = true;
        }
        else if (!std::strcmp(argv[i], "--precode-decoder-cache")) {
            include_precode_cache = true;
            cache_decoder_systematic = true;
        }
        else if (!std::strcmp(argv[i], "--precode-profile")) {
            if (!TakeArg(
                    "compare", "--precode-profile", argc, argv, i, value))
            {
                return 1;
            }
            precode_profile_explicit = true;
            if (!std::strcmp(value, "certified")) {
                precode_profile = PrecodeProfileCertified;
            }
            else if (!std::strcmp(value, "mixed")) {
                precode_profile = PrecodeProfileMixed;
            }
            else if (!std::strcmp(value, "both")) {
                precode_profile = PrecodeProfileBoth;
            }
            else
            {
                std::fprintf(stderr,
                    "unknown --precode-profile '%s' "
                    "(expected certified, mixed, or both)\n",
                    value);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--trial-details")) {
            trial_details = true;
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        else if (!std::strcmp(argv[i], "--mixed-gf16-rows")) {
            if (!TakeArg(
                    "compare", "--mixed-gf16-rows", argc, argv, i, value) ||
                !ParseU32Arg("--mixed-gf16-rows", value, mixed_gf16_rows))
            {
                return 1;
            }
            mixed_gf16_rows_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-period")) {
            if (!TakeArg(
                    "compare", "--mixed-period", argc, argv, i, value) ||
                !ParseU32Arg("--mixed-period", value, mixed_period))
            {
                return 1;
            }
            mixed_period_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-geometry")) {
            if (!TakeArg(
                    "compare", "--mixed-geometry", argc, argv, i, value) ||
                !ParseMixedCoefficientGeometry(value, mixed_geometry))
            {
                std::fprintf(stderr,
                    "compare unknown --mixed-geometry token %s "
                    "(expected frozen or shared-x)\n",
                    value ? value : "");
                return 1;
            }
            mixed_geometry_explicit = true;
        }
#endif
        else if (!std::strcmp(argv[i], "--schedule")) {
            if (!TakeArg("compare", "--schedule", argc, argv, i, value) ||
                !ParsePacketSchedule(value, schedule_kind))
            {
                std::fprintf(stderr,
                    "unknown compare schedule (expected iid, burst, "
                    "permutation, systematic-first, repair-only, or "
                    "adversarial)\n");
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--v2-profile")) {
            if (!TakeArg("compare", "--v2-profile", argc, argv, i, value)) return 1;
            const char* profile = value;
            if (!std::strcmp(profile, "tuned")) {
                std::fprintf(stderr,
                    "--v2-profile tuned is retired: the synthetic peel "
                    "tuner does not transfer to production p_seed values. "
                    "Use --v2-profile auto for real-trial validation, or "
                    "seedtable for offline diagnostics.\n");
                return 1;
            }
            else if (!std::strcmp(profile, "auto")) {
                compare_options.ProfileMode = CompareProfileAuto;
            }
            else if (!std::strcmp(profile, "base")) {
                compare_options.ProfileMode = CompareProfileBase;
            }
            else {
                std::fprintf(stderr,
                    "unknown --v2-profile '%s' (expected base or auto)\n",
                    profile);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--peel-candidates")) {
            if (!TakeArg("compare", "--peel-candidates", argc, argv, i, value) ||
                !ParseU16Arg("--peel-candidates", value,
                    compare_options.PeelCandidates)) return 1;
        }
        else if (!std::strcmp(argv[i], "--peel-trials")) {
            if (!TakeArg("compare", "--peel-trials", argc, argv, i, value) ||
                !ParseU16Arg("--peel-trials", value,
                    compare_options.PeelTrials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--auto-trials")) {
            if (!TakeArg("compare", "--auto-trials", argc, argv, i, value) ||
                !ParseU16Arg("--auto-trials", value,
                    compare_options.AutoTrials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--tune-seed")) {
            if (!TakeArg("compare", "--tune-seed", argc, argv, i, value) ||
                !ParseU64Arg("--tune-seed", value,
                    compare_options.TuneSeed)) return 1;
        }
        else if (!std::strcmp(argv[i], "--auto-seed")) {
            if (!TakeArg("compare", "--auto-seed", argc, argv, i, value) ||
                !ParseU64Arg("--auto-seed", value,
                    compare_options.AutoSeed)) return 1;
        }
        else if (!std::strcmp(argv[i], "--auto-min-delta")) {
            if (!TakeArg("compare", "--auto-min-delta", argc, argv, i, value) ||
                !ParseDoubleArg("--auto-min-delta", value,
                    compare_options.AutoMinDelta)) return 1;
        }
        else if (!std::strcmp(argv[i], "--dense-delta")) {
            if (!TakeArg("compare", "--dense-delta", argc, argv, i, value)) return 1;
            const std::vector<int> deltas = ParseSignedIntList(value);
            if (deltas.size() != 1u) {
                std::fprintf(stderr, "bad --dense-delta value\n");
                return 1;
            }
            compare_options.DenseOverride = true;
            compare_options.DenseDelta = deltas[0];
        }
        else if (!std::strcmp(argv[i], "--dense-candidate")) {
            compare_options.DenseOverride = true;
            if (!TakeArg("compare", "--dense-candidate", argc, argv, i, value) ||
                !ParseU16Arg("--dense-candidate", value,
                    compare_options.DenseCandidate)) return 1;
        }
        else if (!UnknownArg("compare", argv[i])) {
            return 1;
        }
    }

    const std::vector<int> block_bytes_list = ParseIntList(bb_list);
    if (block_bytes_list.empty() || nlo < 2u || nhi < nlo || nhi > 64000u ||
        trials == 0u)
    {
        std::fprintf(stderr, "invalid compare arguments\n");
        return 1;
    }
    if (!ValidateLoss(loss, "compare")) {
        return 1;
    }
    if (precode_profile_explicit &&
        !include_precode && !include_precode_cache)
    {
        std::fprintf(stderr,
            "--precode-profile requires --precode or a precode cache arm\n");
        return 1;
    }
    const bool include_mixed = PrecodeProfileIncludes(
        precode_profile, wirehair_v2::CompletionField::MixedGF256GF16);
    const bool include_certified = PrecodeProfileIncludes(
        precode_profile, wirehair_v2::CompletionField::GF256);
    if (include_mixed)
    {
        for (int bb_value : block_bytes_list)
        {
            if ((bb_value & 1) != 0)
            {
                std::fprintf(stderr,
                    "compare mixed precode profile requires even block bytes "
                    "(odd bb=%d)\n",
                    bb_value);
                return 1;
            }
        }
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    else if (mixed_period_explicit || mixed_geometry_explicit ||
             mixed_gf16_rows_explicit)
    {
        std::fprintf(stderr,
            "compare mixed experiment flags require a mixed precode "
            "profile\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedGF16RowsForTesting(mixed_gf16_rows))
    {
        std::fprintf(stderr,
            "compare --mixed-gf16-rows must be in [%u,%u] and fit the "
            "active period\n",
            wirehair_v2::kMixedGF16Rows,
            wirehair_v2::kMixedGF16RowsMax);
        return 1;
    }
    if (!wirehair_v2::SetMixedCoefficientPeriodForTesting(mixed_period))
    {
        std::fprintf(stderr,
            "compare --mixed-period must be in [%u,%u]\n",
            wirehair_v2::kMixedGF256Rows +
                wirehair_v2::ActiveMixedGF16Rows(),
            wirehair_v2::kMixedCoefficientPeriod);
        return 1;
    }
    if (!wirehair_v2::SetMixedCoefficientGeometryForTesting(mixed_geometry)) {
        return 1;
    }
#endif
    const uint64_t max_message_bytes = max_message_mib > 0u ?
        (uint64_t)max_message_mib * 1024u * 1024u : 0u;
    for (int bb_value : block_bytes_list)
    {
        const uint64_t largest_n = max_message_bytes > 0u ?
            std::min<uint64_t>(
                nhi, max_message_bytes / (uint32_t)bb_value) : nhi;
        if (largest_n < nlo ||
            !ValidateMessageDimensions(
                (uint32_t)largest_n, (uint32_t)bb_value, "compare",
                max_message_bytes))
        {
            std::fprintf(stderr,
                "compare cap cannot accommodate the requested N range for "
                "bb=%d\n", bb_value);
            return 1;
        }
    }
    if (compare_options.PeelCandidates < 1u) {
        compare_options.PeelCandidates = 1u;
    }
    if (compare_options.PeelCandidates > 256u) {
        compare_options.PeelCandidates = 256u;
    }
    if (compare_options.PeelTrials < 1u) {
        compare_options.PeelTrials = 1u;
    }
    if (compare_options.AutoTrials < 1u) {
        compare_options.AutoTrials = 1u;
    }
    if (compare_options.AutoMinDelta < 0.0) {
        compare_options.AutoMinDelta = 0.0;
    }
    if (compare_options.DenseCandidate > 255u) {
        compare_options.DenseCandidate = 255u;
    }
    if (compare_options.DenseOverride &&
        !ValidateDenseDeltas(
            std::vector<int>(1, compare_options.DenseDelta),
            "compare"))
    {
        return 1;
    }
    if (compare_options.DenseOverride)
    {
        std::vector<int> Ns;
        for (uint32_t N = nlo; N <= nhi; ++N) {
            Ns.push_back((int)N);
        }
        if (!ValidateDenseCountsForInputs(
                Ns,
                block_bytes_list,
                std::vector<int>(1, compare_options.DenseDelta),
                "compare"))
        {
            return 1;
        }
    }
    compare_options.LogAutoChoices =
        compare_options.ProfileMode == CompareProfileAuto && nlo == nhi;

    std::printf(
        "# compare: N=[%u,%u] trials/bb=%u loss=%.17g seed=0x%llx "
        "max_message_mib=%u v2_profile=%s peel_candidates=%u peel_trials=%u "
        "auto_trials=%u auto_min_delta=%.4f tune_seed=0x%llx "
        "auto_seed=0x%llx dense_override=%u dense_delta=%d "
        "dense_candidate=%u precode=%u precode_cache=%u "
        "precode_profile=%s encoder_cache=%u decoder_cache=%u schedule=%s "
        "schedule_seed=0x%llx mixed_period=%u mixed_gf16_rows=%u "
        "mixed_geometry=%s "
        "loss_trace=common-id-v2 "
        "precode_profile_handoff=encoder-selected-v1\n",
        nlo,
        nhi,
        trials,
        loss,
        (unsigned long long)seed,
        max_message_mib,
        CompareProfileModeName(compare_options.ProfileMode),
        compare_options.PeelCandidates,
        compare_options.PeelTrials,
        compare_options.AutoTrials,
        compare_options.AutoMinDelta,
        (unsigned long long)compare_options.TuneSeed,
        (unsigned long long)compare_options.AutoSeed,
        compare_options.DenseOverride ? 1u : 0u,
        compare_options.DenseDelta,
        compare_options.DenseCandidate,
        include_precode ? 1u : 0u,
        include_precode_cache ? 1u : 0u,
        PrecodeProfileModeName(precode_profile),
        cache_encoder_source ? 1u : 0u,
        cache_decoder_systematic ? 1u : 0u,
        PacketScheduleName(schedule_kind),
        (unsigned long long)seed,
        wirehair_v2::ActiveMixedCoefficientPeriod(),
        wirehair_v2::ActiveMixedGF16Rows(),
        MixedCoefficientGeometryName(
            wirehair_v2::ActiveMixedCoefficientGeometry()));
    std::printf(
        "%-15s %-8s %-7s %-7s %-10s %-10s %-8s "
        "%-6s %-6s %-6s %-8s "
        "%12s %12s %12s %12s\n",
        "codec",
        "bb",
        "trials",
        "fail",
        "N_mean",
        "OH_mean",
        "OH_sd",
        "OH50",
        "OH95",
        "OH99",
        "OH_max",
        "create_MBps",
        "encode_MBps",
        "decode_MBps",
        "recover_MBps");

    Rng rng(seed);
    for (int bb_value : block_bytes_list)
    {
        const uint32_t block_bytes = (uint32_t)bb_value;
        uint32_t sample_nlo = nlo;
        uint32_t capped_nhi = nhi;
        if (max_message_mib > 0u)
        {
            const uint64_t max_bytes = max_message_bytes;
            const uint32_t max_n = (uint32_t)std::min<uint64_t>(
                64000u, max_bytes / block_bytes);
            if (capped_nhi > max_n) {
                capped_nhi = max_n;
            }
        }

        Accum baseline;
        Accum v2;
        Accum certified_precode;
        Accum mixed_precode;
        Accum certified_precode_cache;
        Accum mixed_precode_cache;
        std::map<uint64_t, CachedCompareProfile> profile_cache;
        for (uint32_t trial = 0; trial < trials; ++trial)
        {
            const uint32_t N =
                sample_nlo + (rng.U32() % (capped_nhi - sample_nlo + 1u));
            const uint64_t trial_seed = rng.Next();
            const wirehair_v2::SeedProfile* profile =
                SelectCompareProfile(
                    N, block_bytes, loss, compare_options, profile_cache);
            const std::vector<uint32_t> packet_schedule = BuildPacketSchedule(
                N, N * 2u + 512u, loss, trial_seed, schedule_kind);
            if (packet_schedule.size() != (size_t)N * 2u + 512u)
            {
                std::fprintf(stderr,
                    "compare schedule generation exceeded its bounded "
                    "candidate budget schedule=%s seed=0x%llx\n",
                    PacketScheduleName(schedule_kind),
                    (unsigned long long)trial_seed);
                return 1;
            }
            const TrialResult baseline_result = RunBaselineTrial(
                N, block_bytes, loss, trial_seed, &packet_schedule);
            const TrialResult v2_result = RunV2Trial(
                N, block_bytes, loss, trial_seed, profile, &packet_schedule);
            TrialResult certified_precode_result = {};
            TrialResult mixed_precode_result = {};
            TrialResult certified_cache_result = {};
            TrialResult mixed_cache_result = {};
            const auto run_precode_profile = [&](
                wirehair_v2::CompletionField completion,
                TrialResult& precode_result,
                TrialResult& cache_result)
            {
                if (include_precode) {
                    precode_result = RunV2PrecodeTrial(
                        N, block_bytes, loss, trial_seed, &packet_schedule,
                        completion);
                }
                if (include_precode_cache) {
                    cache_result = RunV2PrecodeTrial(
                        N, block_bytes, loss, trial_seed, &packet_schedule,
                        completion,
                        cache_encoder_source, cache_decoder_systematic);
                }
            };
            const bool mixed_first =
                precode_profile == PrecodeProfileBoth && (trial & 1u) != 0u;
            if (mixed_first)
            {
                run_precode_profile(
                    wirehair_v2::CompletionField::MixedGF256GF16,
                    mixed_precode_result, mixed_cache_result);
                run_precode_profile(
                    wirehair_v2::CompletionField::GF256,
                    certified_precode_result, certified_cache_result);
            }
            else
            {
                if (include_certified) {
                    run_precode_profile(
                        wirehair_v2::CompletionField::GF256,
                        certified_precode_result, certified_cache_result);
                }
                if (include_mixed) {
                    run_precode_profile(
                        wirehair_v2::CompletionField::MixedGF256GF16,
                        mixed_precode_result, mixed_cache_result);
                }
            }
            AddTrial(baseline, N, baseline_result);
            AddTrial(v2, N, v2_result);
            if (include_precode && include_certified) {
                AddTrial(certified_precode, N, certified_precode_result);
            }
            if (include_precode && include_mixed) {
                AddTrial(mixed_precode, N, mixed_precode_result);
            }
            if (include_precode_cache && include_certified) {
                AddTrial(
                    certified_precode_cache, N, certified_cache_result);
            }
            if (include_precode_cache && include_mixed) {
                AddTrial(mixed_precode_cache, N, mixed_cache_result);
            }
            if (trial_details)
            {
                const int baseline_oh = baseline_result.Ok ?
                    (int)baseline_result.Extra : -1;
                const int v2_oh = v2_result.Ok ?
                    (int)v2_result.Extra : -1;
                const auto print_precode_detail = [&](
                    const char* profile_name,
                    const TrialResult& precode_result,
                    const TrialResult& cache_result)
                {
                    const int precode_oh =
                        include_precode && precode_result.Ok ?
                            (int)precode_result.Extra : -1;
                    const int cached_oh =
                        include_precode_cache && cache_result.Ok ?
                            (int)cache_result.Extra : -1;
                    std::printf(
                        "# paired_trial: schedule=%s seed=0x%llx bb=%u "
                        "trial=%u N=%u precode_profile=%s "
                        "baseline_ok=%u baseline_oh=%d "
                        "v2_ok=%u v2_oh=%d v2_delta=%d precode_ok=%u "
                        "precode_oh=%d precode_delta=%d cached_ok=%u "
                        "cached_oh=%d cached_delta=%d\n",
                        PacketScheduleName(schedule_kind),
                        (unsigned long long)trial_seed,
                        block_bytes, trial, N, profile_name,
                        baseline_result.Ok ? 1u : 0u, baseline_oh,
                        v2_result.Ok ? 1u : 0u, v2_oh,
                        baseline_result.Ok && v2_result.Ok ?
                            v2_oh - baseline_oh : 0,
                        include_precode && precode_result.Ok ? 1u : 0u,
                        precode_oh,
                        include_precode && baseline_result.Ok &&
                            precode_result.Ok ? precode_oh - baseline_oh : 0,
                        include_precode_cache && cache_result.Ok ? 1u : 0u,
                        cached_oh,
                        include_precode_cache && baseline_result.Ok &&
                            cache_result.Ok ? cached_oh - baseline_oh : 0);
                };
                if (include_certified) {
                    print_precode_detail(
                        "certified", certified_precode_result,
                        certified_cache_result);
                }
                if (include_mixed) {
                    print_precode_detail(
                        "mixed", mixed_precode_result, mixed_cache_result);
                }
            }
        }
        PrintAccum("baseline", block_bytes, baseline);
        PrintAccum(
            compare_options.ProfileMode == CompareProfileBase ? "v2" :
            "v2_auto",
            block_bytes,
            v2);
        if (include_precode && include_certified) {
            PrintAccum("v2_precode", block_bytes, certified_precode);
        }
        if (include_precode && include_mixed) {
            PrintAccum("v2_mixed", block_bytes, mixed_precode);
        }
        if (include_precode_cache && include_certified) {
            PrintAccum(
                "v2_cached", block_bytes, certified_precode_cache);
        }
        if (include_precode_cache && include_mixed) {
            PrintAccum(
                "v2_mixed_cached", block_bytes, mixed_precode_cache);
        }
    }
    return 0;
}

int CmdSeedTable(int argc, char** argv)
{
    std::string nlist = "320,1000,3200,6400,12000,32000";
    std::string bb_list = "1280,102400,1048576";
    uint32_t peel_candidates = 32u;
    uint32_t trials = 5u;
    uint64_t seed = UINT64_C(0x5eed7ab1e);

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("seedtable", "--N", argc, argv, i, value)) return 1;
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg("seedtable", "--bb-list", argc, argv, i, value)) return 1;
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--peel-candidates")) {
            if (!TakeArg("seedtable", "--peel-candidates", argc, argv, i, value) ||
                !ParseU32Arg("--peel-candidates", value, peel_candidates)) return 1;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg("seedtable", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg("seedtable", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed)) return 1;
        }
        else if (!UnknownArg("seedtable", argv[i])) {
            return 1;
        }
    }

    const uint32_t requested_peel_candidates = peel_candidates;
    if (peel_candidates < 1u) {
        peel_candidates = 1u;
    }
    if (peel_candidates > 256u) {
        peel_candidates = 256u;
    }
    if (trials < 1u || trials > kMaxSeedTableTrials) {
        std::fprintf(stderr,
            "seedtable --trials must be in [1,%u]\n",
            kMaxSeedTableTrials);
        return 1;
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    if (Ns.empty() || BBs.empty()) {
        std::fprintf(stderr, "seedtable requires non-empty --N and --bb-list\n");
        return 1;
    }
    if (!ValidateBlockCounts(Ns, "seedtable")) {
        return 1;
    }

    std::printf(
        "# seedtable: requested_candidates=%u completed_candidates=%u "
        "trials=%u seed=0x%llx\n",
        requested_peel_candidates,
        peel_candidates,
        trials,
        (unsigned long long)seed);
    std::printf(
        "N,bb,solver,structure,bucket,base_peel,tuned_peel,dense,"
        "requested_candidates,unique_candidates,completed_candidates,"
        "requested_trials,completed_trials,resid_mean,resid_max,xor_cost,"
        "used_peel_fixup,used_dense_fixup\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        wirehair_v2::SeedTuningOptions options =
            wirehair_v2::DefaultSeedTuningOptions();
        options.PeelCandidates = (uint16_t)peel_candidates;
        options.TrialsPerCandidate = trials;
        options.Seed = seed;

        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile((uint32_t)n_value, (uint32_t)bb_value);
        const wirehair_v2::SeedProfile tuned =
            wirehair_v2::TuneSeedProfile((uint32_t)n_value, (uint32_t)bb_value,
                options);
        std::printf(
            "%d,%d,%s,%s,%u,%u,%u,%u,%u,%u,%u,%u,%u,%.4f,%u,%llu,%u,%u\n",
            n_value,
            bb_value,
            wirehair_v2::ToString(tuned.Policy.Solver),
            wirehair_v2::ToString(tuned.Policy.Structure),
            tuned.PeelSeedBucket,
            base.PeelSeed,
            tuned.PeelSeed,
            tuned.DenseSeed,
            requested_peel_candidates,
            (unsigned)tuned.TuningCandidatesUnique,
            (unsigned)tuned.TuningCandidatesCompleted,
            trials,
            tuned.TuningTrials,
            tuned.TuningResidualMean,
            tuned.TuningResidualColumns,
            (unsigned long long)tuned.TuningXorCost,
            tuned.UsedPeelFixup ? 1u : 0u,
            tuned.UsedDenseFixup ? 1u : 0u);
    }
    return 0;
}

int CmdPeelCost(int argc, char** argv)
{
    std::string nlist = "320,3200,32000";
    std::string bb_list = "1280,102400,1048576";
    std::string structure_list =
        "policy,lt_m1_c32,lt_m2_c96,lt_m2_c128,lt_m2_c256,"
        "lt_m2_c512,rs_c001_d50_c128,rs_c003_d10_c128";
    std::string precode_list = "dense,ldpcdense";
    std::string overhead_list = "0,1,2";
    std::string solver_name = "ks_bmax_top16";
    uint32_t trials = 20u;
    uint32_t heavy_rows = 6u;
    uint64_t seed = UINT64_C(0xc057c0de);

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("peelcost", "--N", argc, argv, i, value)) return 1;
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg("peelcost", "--bb-list", argc, argv, i, value)) return 1;
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--structures")) {
            if (!TakeArg("peelcost", "--structures", argc, argv, i, value)) return 1;
            structure_list = value;
        }
        else if (!std::strcmp(argv[i], "--precode")) {
            if (!TakeArg("peelcost", "--precode", argc, argv, i, value)) return 1;
            precode_list = value;
        }
        else if (!std::strcmp(argv[i], "--overhead")) {
            if (!TakeArg("peelcost", "--overhead", argc, argv, i, value)) return 1;
            overhead_list = value;
        }
        else if (!std::strcmp(argv[i], "--solver")) {
            if (!TakeArg("peelcost", "--solver", argc, argv, i, value)) return 1;
            solver_name = value;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg("peelcost", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--heavy")) {
            if (!TakeArg("peelcost", "--heavy", argc, argv, i, value) ||
                !ParseU32Arg("--heavy", value, heavy_rows)) return 1;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg("peelcost", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed)) return 1;
        }
        else if (!UnknownArg("peelcost", argv[i])) {
            return 1;
        }
    }

    if (trials < 1u) {
        trials = 1u;
    }
    if (heavy_rows > 128u) {
        heavy_rows = 128u;
    }

    wirehair_v2::PeelSolver solver = wirehair_v2::PeelSolver::KsBmaxTop16;
    if (!ParsePeelSolver(solver_name, solver)) {
        std::fprintf(stderr, "unknown peel solver: %s\n", solver_name.c_str());
        return 1;
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    const std::vector<int> Overheads = ParseSignedIntList(overhead_list);
    const std::vector<std::string> structure_names =
        ParseStringList(structure_list);
    const std::vector<std::string> precode_names =
        ParseStringList(precode_list);
    if (Ns.empty() || BBs.empty() || Overheads.empty() ||
        structure_names.empty() || precode_names.empty())
    {
        std::fprintf(stderr,
            "peelcost requires non-empty --N, --bb-list, --overhead, "
            "--structures, and --precode\n");
        return 1;
    }
    if (!ValidateBlockCounts(Ns, "peelcost")) {
        return 1;
    }

    std::vector<PeelCostCandidate> candidates;
    for (size_t i = 0; i < structure_names.size(); ++i)
    {
        PeelCostCandidate candidate;
        if (!MakeNamedPeelCostCandidate(structure_names[i], solver, candidate))
        {
            std::fprintf(stderr,
                "unknown peelcost structure: %s\n", structure_names[i].c_str());
            return 1;
        }
        candidates.push_back(candidate);
    }

    std::vector<PrecodeModel> precodes;
    for (size_t i = 0; i < precode_names.size(); ++i)
    {
        PrecodeModel model;
        if (!MakePrecodeModel(precode_names[i], model))
        {
            std::fprintf(stderr,
                "unknown peelcost precode: %s\n", precode_names[i].c_str());
            return 1;
        }
        precodes.push_back(model);
    }

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const uint32_t N = (uint32_t)n_value;
        const uint32_t block_bytes = (uint32_t)bb_value;
        if (N < 2u || N > 65535u || block_bytes < 1u)
        {
            std::fprintf(stderr,
                "peelcost N must be in [2,65535] and bb must be positive\n");
            return 1;
        }
        for (const PrecodeModel& precode : precodes)
        {
            if (IsCodecPortPrecode(precode) && N > CAT_WIREHAIR_MAX_N)
            {
                std::fprintf(stderr,
                    "peelcost precode %s requires N in [%u,%u]\n",
                    precode.Name.c_str(),
                    (unsigned)CAT_WIREHAIR_MIN_N,
                    (unsigned)CAT_WIREHAIR_MAX_N);
                return 1;
            }
        }
        for (int overhead_value : Overheads)
        {
            if (overhead_value < 0 ||
                (uint64_t)N + (uint32_t)overhead_value > UINT16_MAX)
            {
                std::fprintf(stderr,
                    "peelcost overhead must be non-negative and N+overhead "
                    "must be <= 65535\n");
                return 1;
            }
        }
    }

    std::vector<wirehair_v2::PeelEvaluation> evals;
    try {
        evals.reserve(trials);
    }
    catch (const std::bad_alloc&) {
        std::fprintf(stderr,
            "peelcost evaluation storage cannot be allocated for %u trials\n",
            trials);
        return 1;
    }
    catch (const std::length_error&) {
        std::fprintf(stderr,
            "peelcost trial count exceeds evaluation container limits\n");
        return 1;
    }

    std::printf(
        "# peelcost: trials=%u solver=%s heavy=%u seed=0x%llx\n",
        trials,
        solver_name.c_str(),
        heavy_rows,
        (unsigned long long)seed);
    std::printf(
        "# cost model is a block-XOR proxy: matrix_xors + precode_gen + "
        "N*solve_width + solve_width^2/2 + H*solve_width\n");
    std::printf(
        "# precode schemes change width/generation estimates only; "
        "rank and constraint-row effects require precode_sim\n");
    std::printf(
        "# --heavy applies to dense/ldpc/ldpcdense; "
        "codecport schemes use MakeCertifiedParams() H\n");
    std::printf(
        "N,bb,overhead,solver,structure,precode,D,ldpc_cols,dense_rows,H,"
        "trials,resid_cols_mu,resid_cols_max,resid_rows_mu,matrix_refs_mu,"
        "matrix_xors_mu,legacy_total_xors_mu,solve_width_mu,"
        "precode_gen_xors_mu,backsub_xors_mu,ge_block_xors_mu,"
        "heavy_muladds_mu,total_block_xors_mu\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const uint32_t N = (uint32_t)n_value;
        const uint32_t block_bytes = (uint32_t)bb_value;

        for (int overhead_value : Overheads)
        {
            const uint32_t overhead = (uint32_t)overhead_value;
            for (size_t c = 0; c < candidates.size(); ++c)
            {
                wirehair_v2::PeelingCodec codec =
                    candidates[c].UsePolicy ?
                    wirehair_v2::SelectPeelingCodec(N, block_bytes) :
                    candidates[c].Codec;
                codec.Solver = solver;
                codec.SolverCandidateLimit =
                    solver == wirehair_v2::PeelSolver::KsBmaxTop16 ? 16u : 0u;

                evals.clear();
                for (uint32_t trial = 0; trial < trials; ++trial)
                {
                    const uint64_t matrix_seed =
                        PeelCostMatrixSeed(
                            seed, N, overhead, trial);
                    const std::vector<std::vector<uint16_t> > rows =
                        wirehair_v2::GeneratePeelMatrixRows(
                            codec, N, N + overhead, matrix_seed);
                    evals.push_back(
                        wirehair_v2::EvaluatePeelingRows(codec, N, rows));
                }

                for (size_t p = 0; p < precodes.size(); ++p)
                {
                    const PrecodeRecipe recipe =
                        BuildPrecodeRecipe(precodes[p], N, heavy_rows);
                    PeelCostAccum acc;
                    for (size_t trial = 0; trial < evals.size(); ++trial) {
                        AddPeelCostTrial(acc, evals[trial], recipe);
                    }

                    std::printf(
                        "%u,%u,%u,%s,%s,%s,%u,%u,%u,%u,%llu,"
                        "%.4f,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,"
                        "%.4f,%.4f,%.4f,%.4f\n",
                        N,
                        block_bytes,
                        overhead,
                        wirehair_v2::ToString(solver),
                        candidates[c].Name.c_str(),
                        precodes[p].Name.c_str(),
                        recipe.Columns,
                        recipe.LdpcColumns,
                        recipe.DenseRows,
                        recipe.HeavyRows,
                        (unsigned long long)acc.Trials,
                        CostMean(acc.ResidualColumns, acc.Trials),
                        acc.ResidualColumnsMax,
                        CostMean(acc.ResidualRows, acc.Trials),
                        CostMean(acc.MatrixRefs, acc.Trials),
                        CostMean(acc.MatrixXors, acc.Trials),
                        CostMean(acc.LegacyTotalXors, acc.Trials),
                        CostMean(acc.SolveWidth, acc.Trials),
                        CostMean(acc.PrecodeGenXors, acc.Trials),
                        CostMean(acc.BackSubXors, acc.Trials),
                        CostMean(acc.GeBlockXors, acc.Trials),
                        CostMean(acc.HeavyMuladds, acc.Trials),
                        CostMean(acc.TotalBlockXors, acc.Trials));
                }
            }
        }
    }

    return 0;
}

int CmdPrecodeCheck(int argc, char** argv)
{
    std::string nlist = "64,320,1000";
    std::string bb_list = "16,1280";
    uint32_t trials = 10u;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0x7632707265636f64);
    bool trial_details = false;
    PacketScheduleKind schedule_kind = PacketScheduleKind::Iid;

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("precodecheck", "--N", argc, argv, i, value)) {
                return 1;
            }
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg(
                    "precodecheck", "--bb-list", argc, argv, i, value))
            {
                return 1;
            }
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg(
                    "precodecheck", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!TakeArg(
                    "precodecheck", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg(
                    "precodecheck", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--trial-details")) {
            trial_details = true;
        }
        else if (!std::strcmp(argv[i], "--schedule")) {
            if (!TakeArg(
                    "precodecheck", "--schedule", argc, argv, i, value) ||
                !ParsePacketSchedule(value, schedule_kind))
            {
                std::fprintf(stderr,
                    "unknown precodecheck schedule (expected iid, burst, "
                    "permutation, systematic-first, repair-only, or "
                    "adversarial)\n");
                return 1;
            }
        }
        else if (!UnknownArg("precodecheck", argv[i])) {
            return 1;
        }
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    if (Ns.empty() || BBs.empty() || trials == 0u || trials > 10000u) {
        std::fprintf(stderr,
            "precodecheck requires non-empty lists and trials in [1,10000]\n");
        return 1;
    }
    if (!ValidateBlockCounts(Ns, "precodecheck") ||
        !ValidateLoss(loss, "precodecheck") ||
        !ValidateMessageInputs(Ns, BBs, "precodecheck"))
    {
        return 1;
    }

    const wirehair_v2::MessagePrecodeEncoderOptions production_options;
    std::printf(
        "# precodecheck: trials=%u loss=%.17g seed=0x%llx "
        "identity_corner=%u rowdist=wirehair mix=%u schedule=%s "
        "schedule_seed=0x%llx\n",
        trials,
        loss,
        (unsigned long long)seed,
        production_options.DenseIdentityCorner ? 1u : 0u,
        production_options.RecoveryMixCount,
        PacketScheduleName(schedule_kind),
        (unsigned long long)seed);
    std::printf(
        "N,bb,trials,success,fail,terminal_need_more,terminal_invalid,"
        "terminal_extra_insufficient,terminal_other,OH_mean,OH_sd,OH50,"
        "OH95,OH99,OH_max,create_MBps,encode_MBps,decode_MBps,recover_MBps\n");

    for (const int bb_value : BBs) for (const int n_value : Ns)
    {
        const uint32_t N = (uint32_t)n_value;
        const uint32_t block_bytes = (uint32_t)bb_value;
        Accum acc;
        uint64_t terminal_need_more = 0u;
        uint64_t terminal_invalid = 0u;
        uint64_t terminal_extra = 0u;
        uint64_t terminal_other = 0u;
        for (uint32_t trial = 0; trial < trials; ++trial)
        {
            const uint64_t trial_seed = seed ^
                ((uint64_t)N * UINT64_C(0x9e3779b97f4a7c15)) ^
                ((uint64_t)block_bytes * UINT64_C(0xbf58476d1ce4e5b9)) ^
                ((uint64_t)trial * UINT64_C(0xd6e8feb86659fd93));
            const std::vector<uint32_t> packet_schedule = BuildPacketSchedule(
                N, N * 2u + 512u, loss, trial_seed, schedule_kind);
            if (packet_schedule.size() != (size_t)N * 2u + 512u)
            {
                std::fprintf(stderr,
                    "precodecheck schedule generation exceeded its bounded "
                    "candidate budget schedule=%s seed=0x%llx\n",
                    PacketScheduleName(schedule_kind),
                    (unsigned long long)trial_seed);
                return 1;
            }
            const TrialResult result = RunV2PrecodeTrial(
                N, block_bytes, loss, trial_seed, &packet_schedule);
            AddTrial(acc, N, result);
            if (trial_details)
            {
                std::printf(
                    "# precode_trial: schedule=%s seed=0x%llx bb=%u "
                    "trial=%u N=%u ok=%u overhead=%d terminal=%d\n",
                    PacketScheduleName(schedule_kind),
                    (unsigned long long)trial_seed,
                    block_bytes, trial, N, result.Ok ? 1u : 0u,
                    result.Ok ? (int)result.Extra : -1,
                    (int)result.TerminalResult);
            }
            if (result.Ok) {
                continue;
            }
            switch (result.TerminalResult)
            {
            case Wirehair_NeedMore:
                ++terminal_need_more;
                break;
            case Wirehair_InvalidInput:
                ++terminal_invalid;
                break;
            case Wirehair_ExtraInsufficient:
                ++terminal_extra;
                break;
            default:
                ++terminal_other;
                break;
            }
        }
        std::printf(
            "%u,%u,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%.4f,%.4f,%u,%u,"
            "%u,%u,%.1f,%.1f,%.1f,%.1f\n",
            N,
            block_bytes,
            (unsigned long long)acc.Trials,
            (unsigned long long)(acc.Trials - acc.Failures),
            (unsigned long long)acc.Failures,
            (unsigned long long)terminal_need_more,
            (unsigned long long)terminal_invalid,
            (unsigned long long)terminal_extra,
            (unsigned long long)terminal_other,
            MeanOverheadOrFailure(acc),
            OverheadStdDev(acc),
            Percentile(acc.Overheads, 0.50),
            Percentile(acc.Overheads, 0.95),
            Percentile(acc.Overheads, 0.99),
            acc.OverheadMax,
            Mbps(acc.CreateBytes, acc.CreateSeconds),
            Mbps(acc.EncodeBytes, acc.EncodeSeconds),
            Mbps(acc.DecodeBytes, acc.DecodeSeconds),
            Mbps(acc.RecoverBytes, acc.RecoverSeconds));
    }
    return 0;
}

int CmdDenseCheck(int argc, char** argv)
{
    uint32_t N = 7533u;
    uint32_t block_bytes = 1280u;
    uint32_t candidates = 8u;
    uint32_t trials = 2u;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0xded5eed);

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("densecheck", "--N", argc, argv, i, value) ||
                !ParseU32Arg("--N", value, N)) return 1;
        }
        else if (!std::strcmp(argv[i], "--bb")) {
            if (!TakeArg("densecheck", "--bb", argc, argv, i, value) ||
                !ParseU32Arg("--bb", value, block_bytes)) return 1;
        }
        else if (!std::strcmp(argv[i], "--candidates")) {
            if (!TakeArg("densecheck", "--candidates", argc, argv, i, value) ||
                !ParseU32Arg("--candidates", value, candidates)) return 1;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg("densecheck", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!TakeArg("densecheck", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss)) return 1;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg("densecheck", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed)) return 1;
        }
        else if (!UnknownArg("densecheck", argv[i])) {
            return 1;
        }
    }

    if (candidates < 1u) {
        candidates = 1u;
    }
    if (candidates > 256u) {
        candidates = 256u;
    }
    if (trials < 1u) {
        trials = 1u;
    }
    if (!ValidateLoss(loss, "densecheck") ||
        !ValidateMessageDimensions(N, block_bytes, "densecheck"))
    {
        return 1;
    }
    const wirehair_v2::SeedProfile base =
        wirehair_v2::SelectSeedProfile(N, block_bytes);
    std::printf(
        "# densecheck: N=%u bb=%u base_dense=%u candidates=%u trials=%u "
        "loss=%.17g\n",
        N, block_bytes, base.DenseSeed, candidates, trials, loss);
    std::printf("%-8s %-8s %-10s %-8s\n",
        "dense", "fail", "OH_mean", "OH_max");

    for (uint32_t c = 0; c < candidates; ++c)
    {
        wirehair_v2::SeedProfile profile = base;
        profile.DenseSeed =
            wirehair_v2::CandidateDenseSeed(base.DenseSeed, (uint16_t)c);
        Accum acc;
        for (uint32_t trial = 0; trial < trials; ++trial)
        {
            // Candidate index deliberately excluded: paired trials across
            // candidates (common random numbers).
            const uint64_t trial_seed = seed ^
                ((uint64_t)trial * UINT64_C(0xbf58476d1ce4e5b9));
            AddTrial(acc, N,
                RunV2Trial(N, block_bytes, loss, trial_seed, &profile));
        }
        std::printf("%-8u %-8llu %-10.4f %-8u\n",
            profile.DenseSeed,
            (unsigned long long)acc.Failures,
            MeanOverheadOrFailure(acc),
            acc.OverheadMax);
    }
    return 0;
}

int CmdDenseTune(int argc, char** argv)
{
    std::string nlist = "320,1000,3200";
    std::string bb_list = "1280,1048576";
    uint32_t candidates = 32u;
    uint32_t trials = 4u;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0xded5eed);

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("densetune", "--N", argc, argv, i, value)) return 1;
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg("densetune", "--bb-list", argc, argv, i, value)) return 1;
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--candidates")) {
            if (!TakeArg("densetune", "--candidates", argc, argv, i, value) ||
                !ParseU32Arg("--candidates", value, candidates)) return 1;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg("densetune", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!TakeArg("densetune", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss)) return 1;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg("densetune", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed)) return 1;
        }
        else if (!UnknownArg("densetune", argv[i])) {
            return 1;
        }
    }

    const uint32_t requested_candidates = candidates;
    if (candidates < 1u) {
        candidates = 1u;
    }
    if (candidates > 256u) {
        candidates = 256u;
    }
    if (trials < 1u) {
        trials = 1u;
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    if (Ns.empty() || BBs.empty()) {
        std::fprintf(stderr, "densetune requires non-empty --N and --bb-list\n");
        return 1;
    }
    if (!ValidateBlockCounts(Ns, "densetune")) {
        return 1;
    }
    if (!ValidateLoss(loss, "densetune") ||
        !ValidateMessageInputs(Ns, BBs, "densetune"))
    {
        return 1;
    }

    std::printf(
        "# densetune: requested_candidates=%u completed_candidates=%u "
        "trials=%u loss=%.17g seed=0x%llx\n",
        requested_candidates,
        candidates,
        trials,
        loss,
        (unsigned long long)seed);
    std::printf(
        "N,bb,requested_candidates,unique_candidates,completed_candidates,"
        "base_dense,best_dense,best_fail,best_oh_mean,best_oh_max,"
        "base_fail,base_oh_mean,base_oh_max\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile((uint32_t)n_value, (uint32_t)bb_value);
        Accum base_acc;
        uint16_t best_dense = base.DenseSeed;
        Accum best_acc;
        bool have_best = false;

        for (uint32_t c = 0; c < candidates; ++c)
        {
            wirehair_v2::SeedProfile profile = base;
            profile.DenseSeed =
                wirehair_v2::CandidateDenseSeed(base.DenseSeed, (uint16_t)c);
            Accum acc;
            for (uint32_t trial = 0; trial < trials; ++trial)
            {
                // Candidate index is deliberately excluded so every dense
                // seed candidate sees identical message/loss trials, same
                // as the densecount/densegrid paired protocol.
                const uint64_t trial_seed = seed ^
                    ((uint64_t)n_value * UINT64_C(0x9e3779b97f4a7c15)) ^
                    ((uint64_t)bb_value * UINT64_C(0xbf58476d1ce4e5b9)) ^
                    ((uint64_t)trial * UINT64_C(0xd6e8feb86659fd93));
                AddTrial(acc, (uint32_t)n_value,
                    RunV2Trial((uint32_t)n_value, (uint32_t)bb_value,
                        loss, trial_seed, &profile));
            }
            if (c == 0u) {
                base_acc = acc;
            }
            const uint64_t ok = acc.Trials - acc.Failures;
            const uint64_t best_ok = best_acc.Trials - best_acc.Failures;
            const double mean =
                ok > 0u ? (double)acc.OverheadSum / (double)ok : 1e9;
            const double best_mean =
                best_ok > 0u ? (double)best_acc.OverheadSum /
                    (double)best_ok : 1e9;
            if (!have_best ||
                acc.Failures < best_acc.Failures ||
                (acc.Failures == best_acc.Failures && mean < best_mean) ||
                (acc.Failures == best_acc.Failures && mean == best_mean &&
                 acc.OverheadMax < best_acc.OverheadMax))
            {
                best_dense = profile.DenseSeed;
                best_acc = acc;
                have_best = true;
            }
        }

        const double base_mean = MeanOverheadOrFailure(base_acc);
        const double best_mean = MeanOverheadOrFailure(best_acc);
        std::printf(
            "%d,%d,%u,%u,%u,%u,%u,%llu,%.4f,%u,%llu,%.4f,%u\n",
            n_value,
            bb_value,
            requested_candidates,
            candidates,
            candidates,
            base.DenseSeed,
            best_dense,
            (unsigned long long)best_acc.Failures,
            best_mean,
            best_acc.OverheadMax,
            (unsigned long long)base_acc.Failures,
            base_mean,
            base_acc.OverheadMax);
    }
    return 0;
}

int CmdDenseCount(int argc, char** argv)
{
    std::string nlist = "320,1000,3200";
    std::string bb_list = "1280,102400";
    std::string delta_list = "-16,-8,-4,0,4,8,16,32";
    uint32_t trials = 100u;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0xdeca770);

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("densecount", "--N", argc, argv, i, value)) return 1;
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg("densecount", "--bb-list", argc, argv, i, value)) return 1;
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--deltas")) {
            if (!TakeArg("densecount", "--deltas", argc, argv, i, value)) return 1;
            delta_list = value;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg("densecount", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!TakeArg("densecount", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss)) return 1;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg("densecount", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed)) return 1;
        }
        else if (!UnknownArg("densecount", argv[i])) {
            return 1;
        }
    }

    if (trials < 1u) {
        trials = 1u;
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    const std::vector<int> Deltas = ParseSignedIntList(delta_list);
    if (Ns.empty() || BBs.empty() || Deltas.empty()) {
        std::fprintf(stderr,
            "densecount requires non-empty --N, --bb-list, and --deltas\n");
        return 1;
    }
    if (!ValidateBlockCounts(Ns, "densecount")) {
        return 1;
    }
    if (!ValidateDenseDeltas(Deltas, "densecount")) {
        return 1;
    }
    if (!ValidateDenseCountsForInputs(Ns, BBs, Deltas, "densecount")) {
        return 1;
    }
    if (!ValidateLoss(loss, "densecount") ||
        !ValidateMessageInputs(Ns, BBs, "densecount"))
    {
        return 1;
    }

    std::printf(
        "# densecount: trials=%u loss=%.17g seed=0x%llx\n",
        trials,
        loss,
        (unsigned long long)seed);
    std::printf(
        "N,bb,base_dense,dense,delta,dense_seed,fail,OH_mean,OH_sd,"
        "OH50,OH95,OH99,OH_max,create_MBps,decode_MBps,recover_MBps\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const uint32_t N = (uint32_t)n_value;
        const uint32_t block_bytes = (uint32_t)bb_value;
        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile(N, block_bytes);

        for (int delta : Deltas)
        {
            uint16_t dense_count = 0u;
            if (!DenseCountForDelta(
                    base, delta, "densecount", N, block_bytes, &dense_count))
            {
                return 1;
            }

            wirehair_v2::SeedProfile profile = base;
            profile.DenseCount = dense_count;
            profile.DenseSeed = wirehair::GetDenseSeed(N, profile.DenseCount);

            Accum acc;
            for (uint32_t trial = 0; trial < trials; ++trial)
            {
                const uint64_t trial_seed =
                    seed ^
                    ((uint64_t)N * UINT64_C(0x9e3779b97f4a7c15)) ^
                    ((uint64_t)block_bytes * UINT64_C(0xbf58476d1ce4e5b9)) ^
                    ((uint64_t)trial * UINT64_C(0xd6e8feb86659fd93));
                AddTrial(acc, N,
                    RunV2Trial(N, block_bytes, loss, trial_seed, &profile));
            }

            std::printf(
                "%u,%u,%u,%u,%d,%u,%llu,%.4f,%.4f,%u,%u,%u,%u,"
                "%.1f,%.1f,%.1f\n",
                N,
                block_bytes,
                base.DenseCount,
                profile.DenseCount,
                delta,
                profile.DenseSeed,
                (unsigned long long)acc.Failures,
                MeanOverheadOrFailure(acc),
                OverheadStdDev(acc),
                Percentile(acc.Overheads, 0.50),
                Percentile(acc.Overheads, 0.95),
                Percentile(acc.Overheads, 0.99),
                acc.OverheadMax,
                Mbps(acc.CreateBytes, acc.CreateSeconds),
                Mbps(acc.DecodeBytes, acc.DecodeSeconds),
                Mbps(acc.RecoverBytes, acc.RecoverSeconds));
        }
    }

    return 0;
}

int CmdDenseGrid(int argc, char** argv)
{
    std::string nlist = "320,1000";
    std::string bb_list = "1280,102400";
    std::string delta_list = "-8,0,8";
    uint32_t candidates = 8u;
    uint32_t trials = 80u;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0xded6e12d);

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("densegrid", "--N", argc, argv, i, value)) return 1;
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg("densegrid", "--bb-list", argc, argv, i, value)) return 1;
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--deltas")) {
            if (!TakeArg("densegrid", "--deltas", argc, argv, i, value)) return 1;
            delta_list = value;
        }
        else if (!std::strcmp(argv[i], "--candidates")) {
            if (!TakeArg("densegrid", "--candidates", argc, argv, i, value) ||
                !ParseU32Arg("--candidates", value, candidates)) return 1;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg("densegrid", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials)) return 1;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!TakeArg("densegrid", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss)) return 1;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg("densegrid", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed)) return 1;
        }
        else if (!UnknownArg("densegrid", argv[i])) {
            return 1;
        }
    }

    if (candidates < 1u) {
        candidates = 1u;
    }
    if (candidates > 256u) {
        candidates = 256u;
    }
    if (trials < 1u) {
        trials = 1u;
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    const std::vector<int> Deltas = ParseSignedIntList(delta_list);
    if (Ns.empty() || BBs.empty() || Deltas.empty()) {
        std::fprintf(stderr,
            "densegrid requires non-empty --N, --bb-list, and --deltas\n");
        return 1;
    }
    if (!ValidateBlockCounts(Ns, "densegrid")) {
        return 1;
    }
    if (!ValidateDenseDeltas(Deltas, "densegrid")) {
        return 1;
    }
    if (!ValidateDenseCountsForInputs(Ns, BBs, Deltas, "densegrid")) {
        return 1;
    }
    if (!ValidateLoss(loss, "densegrid") ||
        !ValidateMessageInputs(Ns, BBs, "densegrid"))
    {
        return 1;
    }

    std::printf(
        "# densegrid: candidates=%u trials=%u loss=%.17g seed=0x%llx\n",
        candidates,
        trials,
        loss,
        (unsigned long long)seed);
    std::printf(
        "N,bb,base_dense,dense,delta,table_dense_seed,candidate,dense_seed,"
        "fail,OH_mean,OH_sd,OH50,OH95,OH99,OH_max,create_MBps,decode_MBps,"
        "recover_MBps\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const uint32_t N = (uint32_t)n_value;
        const uint32_t block_bytes = (uint32_t)bb_value;
        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile(N, block_bytes);

        for (int delta : Deltas)
        {
            uint16_t dense_count_u = 0u;
            if (!DenseCountForDelta(
                    base, delta, "densegrid", N, block_bytes, &dense_count_u))
            {
                return 1;
            }
            const uint16_t table_seed =
                wirehair::GetDenseSeed(N, dense_count_u);

            for (uint32_t candidate = 0; candidate < candidates; ++candidate)
            {
                wirehair_v2::SeedProfile profile = base;
                profile.DenseCount = dense_count_u;
                profile.DenseSeed = wirehair_v2::CandidateDenseSeed(
                    table_seed, (uint16_t)candidate);

                Accum acc;
                for (uint32_t trial = 0; trial < trials; ++trial)
                {
                    const uint64_t trial_seed =
                        seed ^
                        ((uint64_t)N * UINT64_C(0x9e3779b97f4a7c15)) ^
                        ((uint64_t)block_bytes * UINT64_C(0xbf58476d1ce4e5b9)) ^
                        ((uint64_t)trial * UINT64_C(0xd6e8feb86659fd93));
                    AddTrial(acc, N,
                        RunV2Trial(N, block_bytes, loss, trial_seed, &profile));
                }

                std::printf(
                    "%u,%u,%u,%u,%d,%u,%u,%u,%llu,%.4f,%.4f,%u,%u,%u,%u,"
                    "%.1f,%.1f,%.1f\n",
                    N,
                    block_bytes,
                    base.DenseCount,
                    profile.DenseCount,
                    delta,
                    table_seed,
                    candidate,
                    profile.DenseSeed,
                    (unsigned long long)acc.Failures,
                    MeanOverheadOrFailure(acc),
                    OverheadStdDev(acc),
                    Percentile(acc.Overheads, 0.50),
                    Percentile(acc.Overheads, 0.95),
                    Percentile(acc.Overheads, 0.99),
                    acc.OverheadMax,
                    Mbps(acc.CreateBytes, acc.CreateSeconds),
                    Mbps(acc.DecodeBytes, acc.DecodeSeconds),
                    Mbps(acc.RecoverBytes, acc.RecoverSeconds));
            }
        }
    }

    return 0;
}

struct MatrixFailureTrial
{
    WirehairResult Result = Wirehair_Error;
    uint32_t Inactivated = 0u;
    uint32_t Rank = 0u;
    uint32_t BinaryRank = 0u;
    uint64_t BlockXors = 0u;
    uint64_t BlockMulAdds = 0u;
    uint64_t SolveNanoseconds = 0u;
    uint64_t BuildNanoseconds = 0u;
    uint64_t PeelNanoseconds = 0u;
    uint64_t ProjectNanoseconds = 0u;
    uint64_t ResidualNanoseconds = 0u;
    uint64_t BackSubNanoseconds = 0u;
};

const char* HeavyFamilyName(wirehair_v2::HeavyCoefficientFamily family);

enum class PrecodeFailCompletion
{
    Certified,
    Mixed
};

const char* PrecodeFailCompletionName(PrecodeFailCompletion completion)
{
    return completion == PrecodeFailCompletion::Mixed ?
        "mixed" : "certified";
}

bool ParsePrecodeFailCompletion(
    const char* text,
    PrecodeFailCompletion& completion)
{
    if (!std::strcmp(text, "certified")) {
        completion = PrecodeFailCompletion::Certified;
        return true;
    }
    if (!std::strcmp(text, "mixed")) {
        completion = PrecodeFailCompletion::Mixed;
        return true;
    }
    return false;
}

struct PairedMixOutcomes
{
    std::vector<bool> Mix2Failures;
    std::vector<bool> Mix3Failures;
    uint32_t Mix2SeedAttempt = UINT32_MAX;
    uint32_t Mix3SeedAttempt = UINT32_MAX;
};

void Wilson95(
    uint32_t failures,
    uint32_t trials,
    double& lower,
    double& upper)
{
    if (trials == 0u) {
        lower = upper = 0.0;
        return;
    }
    static const double z = 1.959963984540054;
    const double n = (double)trials;
    const double p = (double)failures / n;
    const double z2_over_n = z * z / n;
    const double center = (p + z2_over_n / 2.0) / (1.0 + z2_over_n);
    const double radius = z * std::sqrt(
        (p * (1.0 - p) / n) + z * z / (4.0 * n * n)) /
        (1.0 + z2_over_n);
    lower = std::max(0.0, center - radius);
    upper = std::min(1.0, center + radius);
}

double ExactMcNemarP(uint32_t first_only, uint32_t second_only)
{
    const uint32_t discordant = first_only + second_only;
    if (discordant == 0u) {
        return 1.0;
    }
    const uint32_t tail = std::min(first_only, second_only);
    const long double n = (long double)discordant;
    const long double k = (long double)tail;
    const long double log_probability_at_tail =
        std::lgamma(n + 1.0L) - std::lgamma(k + 1.0L) -
        std::lgamma(n - k + 1.0L) - n * std::log(2.0L);
    long double relative_sum = 1.0L;
    long double relative_term = 1.0L;
    for (uint32_t j = tail; j > 0u; --j)
    {
        relative_term *= (long double)j /
            (long double)(discordant - j + 1u);
        relative_sum += relative_term;
    }
    const long double doubled_tail = 2.0L * std::exp(
        log_probability_at_tail + std::log(relative_sum));
    return (double)std::min(1.0L, doubled_tail);
}

void PrintPairedMixOutcomes(
    uint32_t K,
    uint32_t block_bytes,
    PrecodeFailCompletion completion,
    wirehair_v2::HeavyCoefficientFamily heavy_family,
    uint32_t overhead,
    const PairedMixOutcomes& outcomes)
{
    if (outcomes.Mix2Failures.empty() ||
        outcomes.Mix2Failures.size() != outcomes.Mix3Failures.size())
    {
        return;
    }
    uint32_t both_success = 0u;
    uint32_t both_failure = 0u;
    uint32_t mix2_only = 0u;
    uint32_t mix3_only = 0u;
    for (size_t i = 0; i < outcomes.Mix2Failures.size(); ++i)
    {
        const bool mix2_failed = outcomes.Mix2Failures[i] != 0u;
        const bool mix3_failed = outcomes.Mix3Failures[i] != 0u;
        if (mix2_failed && mix3_failed) {
            ++both_failure;
        }
        else if (mix2_failed) {
            ++mix2_only;
        }
        else if (mix3_failed) {
            ++mix3_only;
        }
        else {
            ++both_success;
        }
    }
    const uint32_t trials = (uint32_t)outcomes.Mix2Failures.size();
    const uint32_t mix2_failures = both_failure + mix2_only;
    const uint32_t mix3_failures = both_failure + mix3_only;
    double mix2_lower = 0.0, mix2_upper = 0.0;
    double mix3_lower = 0.0, mix3_upper = 0.0;
    Wilson95(mix2_failures, trials, mix2_lower, mix2_upper);
    Wilson95(mix3_failures, trials, mix3_lower, mix3_upper);
    std::printf(
        "# precodefail_paired: N=%u bb=%u completion=%s heavy_family=%s "
        "overhead=%u trials=%u mix2_fail=%u mix3_fail=%u both_fail=%u "
        "mix2_only=%u mix3_only=%u both_success=%u mcnemar_p=%.12g "
        "mix2_seed_attempt=%u mix3_seed_attempt=%u "
        "mix2_wilson95=[%.8f,%.8f] mix3_wilson95=[%.8f,%.8f]\n",
        K, block_bytes, PrecodeFailCompletionName(completion),
        HeavyFamilyName(heavy_family), overhead, trials,
        mix2_failures, mix3_failures, both_failure,
        mix2_only, mix3_only, both_success,
        ExactMcNemarP(mix2_only, mix3_only),
        outcomes.Mix2SeedAttempt, outcomes.Mix3SeedAttempt,
        mix2_lower, mix2_upper, mix3_lower, mix3_upper);
}

const char* HeavyFamilyName(wirehair_v2::HeavyCoefficientFamily family)
{
    return family == wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ?
        "periodic" : "hashed";
}

bool ParseHeavyFamilies(
    const std::string& text,
    std::vector<wirehair_v2::HeavyCoefficientFamily>& families)
{
    families.clear();
    for (const std::string& token : ParseStringList(text))
    {
        wirehair_v2::HeavyCoefficientFamily family;
        if (token == "periodic") {
            family = wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy;
        }
        else if (token == "hashed") {
            family = wirehair_v2::HeavyCoefficientFamily::HashedNonzero;
        }
        else {
            std::fprintf(stderr,
                "precodefail unknown --heavy-family token %s\n",
                token.c_str());
            return false;
        }
        if (std::find(families.begin(), families.end(), family) ==
            families.end())
        {
            families.push_back(family);
        }
    }
    return !families.empty();
}

std::string CountHistogram(const std::map<uint32_t, uint32_t>& histogram)
{
    std::string out;
    for (const std::pair<const uint32_t, uint32_t>& bin : histogram)
    {
        if (!out.empty()) {
            out += '|';
        }
        out += std::to_string(bin.first);
        out += ':';
        out += std::to_string(bin.second);
    }
    return out;
}

bool RunPrecodeFailPayloadE2E(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_bytes,
    double loss,
    uint64_t seed)
{
    const uint32_t K = system.Params.BlockCount;
    const uint64_t message_bytes_wide = (uint64_t)K * block_bytes;
    const uint64_t delivered_bytes_wide =
        ((uint64_t)K + 32u) * block_bytes;
    if (block_bytes == 0u || block_bytes > UINT32_C(0x7fffffff) ||
        message_bytes_wide > (uint64_t)SIZE_MAX ||
        delivered_bytes_wide > (uint64_t)SIZE_MAX)
    {
        return false;
    }

    std::vector<uint8_t> message((size_t)message_bytes_wide);
    Rng data_rng(seed ^ UINT64_C(0x7061796c6f616432));
    for (uint8_t& byte : message) {
        byte = (uint8_t)data_rng.U32();
    }
    std::vector<wirehair_v2::SolvePacket> systematic(K);
    for (uint32_t id = 0; id < K; ++id)
    {
        systematic[id].BlockId = id;
        systematic[id].Data = message.data() + (size_t)id * block_bytes;
    }
    std::vector<uint8_t> intermediate;
    if (wirehair_v2::SolvePrecodeSystem(
            system, config, systematic, block_bytes, intermediate) !=
        Wirehair_Success)
    {
        return false;
    }

    std::vector<uint8_t> delivered_storage((size_t)delivered_bytes_wide);
    std::vector<wirehair_v2::SolvePacket> delivered;
    delivered.reserve((size_t)K + 32u);
    std::vector<uint8_t> recovered;
    Rng loss_rng(seed ^ UINT64_C(0x6c6f737370617932));
    WirehairResult result = Wirehair_NeedMore;
    uint32_t block_id = 0u;
    while (delivered.size() < (size_t)K + 32u &&
           result == Wirehair_NeedMore)
    {
        const uint32_t id = block_id++;
        if (ShouldDrop(loss_rng, loss)) {
            continue;
        }
        uint8_t* block = delivered_storage.data() +
            delivered.size() * block_bytes;
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, intermediate.data(), block_bytes, id, block))
        {
            return false;
        }
        wirehair_v2::SolvePacket packet;
        packet.BlockId = id;
        packet.Data = block;
        delivered.push_back(packet);
        if (delivered.size() >= K) {
            result = wirehair_v2::SolvePrecodeSystem(
                system, config, delivered, block_bytes, recovered);
        }
    }
    if (result != Wirehair_Success ||
        !wirehair_v2::VerifyPrecodeSolution(
            system, config, delivered, recovered.data(), block_bytes))
    {
        return false;
    }

    std::vector<uint8_t> block(block_bytes);
    for (uint32_t id = 0; id < K; ++id)
    {
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystem(
                system, config, recovered.data(), block_bytes,
                id, block.data()) ||
            std::memcmp(
                block.data(), message.data() + (size_t)id * block_bytes,
                block_bytes) != 0)
        {
            return false;
        }
    }
    std::printf(
        "# payload_e2e: N=%u bb=%u mix_count=%u delivered=%zu PASS\n",
        K, block_bytes, config.MixCount, delivered.size());
    return true;
}

class ThreadJoinGuard
{
public:
    explicit ThreadJoinGuard(std::vector<std::thread>& threads)
        : Threads(threads)
    {
    }

    ~ThreadJoinGuard()
    {
        JoinAll();
    }

    void JoinAll() noexcept
    {
        for (std::thread& thread : Threads)
        {
            if (!thread.joinable()) {
                continue;
            }
            try {
                thread.join();
            }
            catch (...) {
                try {
                    thread.detach();
                }
                catch (...) {
                }
            }
        }
    }

private:
    std::vector<std::thread>& Threads;
};

int CmdPrecodeFail(int argc, char** argv)
{
    std::string nlist = "1000,3200,10000,32000,64000";
    std::string bb_list = "1280";
    std::string overhead_list = "0,1";
    std::string heavy_family_list = "periodic";
    std::string mix_count_list = "3";
    PrecodeFailCompletion completion = PrecodeFailCompletion::Certified;
    bool payload_e2e = false;
    bool full_payload_solve = false;
    uint32_t trials = 100u;
    uint32_t threads = 1u;
    double loss = 0.10;
    uint64_t seed = UINT64_C(0x5eedf411);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    uint32_t fail_thread_launch_after = UINT32_MAX;
    uint32_t mixed_period = wirehair_v2::kMixedCoefficientPeriod;
    bool mixed_period_explicit = false;
    uint32_t mixed_gf16_rows = wirehair_v2::kMixedGF16Rows;
    bool mixed_gf16_rows_explicit = false;
    wirehair_v2::MixedCoefficientGeometry mixed_geometry =
        wirehair_v2::MixedCoefficientGeometry::FrozenPowerX;
    bool mixed_geometry_explicit = false;
#endif

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!TakeArg("precodefail", "--N", argc, argv, i, value)) {
                return 1;
            }
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!TakeArg(
                    "precodefail", "--bb-list", argc, argv, i, value))
            {
                return 1;
            }
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--overhead")) {
            if (!TakeArg(
                    "precodefail", "--overhead", argc, argv, i, value))
            {
                return 1;
            }
            overhead_list = value;
        }
        else if (!std::strcmp(argv[i], "--heavy-family")) {
            if (!TakeArg(
                    "precodefail", "--heavy-family", argc, argv, i, value))
            {
                return 1;
            }
            heavy_family_list = value;
        }
        else if (!std::strcmp(argv[i], "--mix-count")) {
            if (!TakeArg(
                    "precodefail", "--mix-count", argc, argv, i, value))
            {
                return 1;
            }
            mix_count_list = value;
        }
        else if (!std::strcmp(argv[i], "--completion")) {
            if (!TakeArg(
                    "precodefail", "--completion", argc, argv, i, value))
            {
                return 1;
            }
            if (!ParsePrecodeFailCompletion(value, completion))
            {
                std::fprintf(stderr,
                    "precodefail unknown --completion token %s "
                    "(expected certified or mixed)\n",
                    value);
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--payload-e2e")) {
            payload_e2e = true;
        }
        else if (!std::strcmp(argv[i], "--trials")) {
            if (!TakeArg(
                    "precodefail", "--trials", argc, argv, i, value) ||
                !ParseU32Arg("--trials", value, trials))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--threads")) {
            if (!TakeArg(
                    "precodefail", "--threads", argc, argv, i, value) ||
                !ParseU32Arg("--threads", value, threads))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!TakeArg(
                    "precodefail", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!TakeArg(
                    "precodefail", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed))
            {
                return 1;
            }
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        else if (!std::strcmp(argv[i], "--full-payload-solve")) {
            full_payload_solve = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-gf16-rows")) {
            if (!TakeArg(
                    "precodefail", "--mixed-gf16-rows",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-gf16-rows", value, mixed_gf16_rows))
            {
                return 1;
            }
            mixed_gf16_rows_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-period")) {
            if (!TakeArg(
                    "precodefail", "--mixed-period", argc, argv, i, value) ||
                !ParseU32Arg("--mixed-period", value, mixed_period))
            {
                return 1;
            }
            mixed_period_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-geometry")) {
            if (!TakeArg(
                    "precodefail", "--mixed-geometry",
                    argc, argv, i, value) ||
                !ParseMixedCoefficientGeometry(value, mixed_geometry))
            {
                std::fprintf(stderr,
                    "precodefail unknown --mixed-geometry token %s "
                    "(expected frozen or shared-x)\n",
                    value ? value : "");
                return 1;
            }
            mixed_geometry_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--fail-thread-launch-after")) {
            if (!TakeArg(
                    "precodefail", "--fail-thread-launch-after",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--fail-thread-launch-after", value,
                    fail_thread_launch_after))
            {
                return 1;
            }
        }
#endif
        else if (!UnknownArg("precodefail", argv[i])) {
            return 1;
        }
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    const std::vector<int> overheads = ParseSignedIntList(overhead_list);
    const std::vector<int> mix_counts = ParseIntList(mix_count_list);
    std::vector<wirehair_v2::HeavyCoefficientFamily> heavy_families;
    if (Ns.empty() || BBs.empty() || overheads.empty() || mix_counts.empty() ||
        !ParseHeavyFamilies(heavy_family_list, heavy_families)) {
        std::fprintf(stderr,
            "precodefail requires non-empty integer lists\n");
        return 1;
    }
    if (!ValidateBlockCounts(Ns, "precodefail") ||
        !ValidateLoss(loss, "precodefail") ||
        ((payload_e2e || full_payload_solve) &&
         !ValidatePayloadE2EInputs(
             Ns, BBs, "precodefail",
             full_payload_solve ?
                 std::max(1u, std::min(threads, trials)) : 1u)))
    {
        return 1;
    }
    if (trials == 0u || trials > kMaxSeedTableTrials) {
        std::fprintf(stderr,
            "precodefail --trials must be in [1,%u]\n",
            kMaxSeedTableTrials);
        return 1;
    }
    if (threads == 0u || threads > 256u) {
        std::fprintf(stderr,
            "precodefail --threads must be in [1,256]\n");
        return 1;
    }
    for (int overhead : overheads) {
        if (overhead < 0 || overhead > 1024) {
            std::fprintf(stderr,
                "precodefail overhead must be in [0,1024]\n");
            return 1;
        }
    }
    for (int mix_count : mix_counts) {
        if (mix_count < 1 ||
            mix_count > (int)wirehair_v2::kCertifiedPacketMixCount)
        {
            std::fprintf(stderr,
                "precodefail mix count must be in [1,%u]\n",
                wirehair_v2::kCertifiedPacketMixCount);
            return 1;
        }
    }
    const bool pair_mix_counts =
        std::find(mix_counts.begin(), mix_counts.end(), 2) !=
            mix_counts.end() &&
        std::find(mix_counts.begin(), mix_counts.end(), 3) !=
            mix_counts.end();
    if (completion == PrecodeFailCompletion::Mixed)
    {
        for (int bb : BBs)
        {
            if ((bb & 1) != 0)
            {
                std::fprintf(stderr,
                    "precodefail mixed completion requires even block bytes, "
                    "got %d\n",
                    bb);
                return 1;
            }
        }
        for (wirehair_v2::HeavyCoefficientFamily family : heavy_families)
        {
            if (family !=
                wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy)
            {
                std::fprintf(stderr,
                    "precodefail mixed completion requires periodic "
                    "heavy family\n");
                return 1;
            }
        }
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    else if (mixed_period_explicit || mixed_geometry_explicit ||
             mixed_gf16_rows_explicit)
    {
        std::fprintf(stderr,
            "precodefail mixed experiment flags require --completion "
            "mixed\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedGF16RowsForTesting(mixed_gf16_rows))
    {
        std::fprintf(stderr,
            "precodefail --mixed-gf16-rows must be in [%u,%u] and fit "
            "the active period\n",
            wirehair_v2::kMixedGF16Rows,
            wirehair_v2::kMixedGF16RowsMax);
        return 1;
    }
    if (!wirehair_v2::SetMixedCoefficientPeriodForTesting(mixed_period))
    {
        std::fprintf(stderr,
            "precodefail --mixed-period must be in [%u,%u]\n",
            wirehair_v2::kMixedGF256Rows +
                wirehair_v2::ActiveMixedGF16Rows(),
            wirehair_v2::kMixedCoefficientPeriod);
        return 1;
    }
    if (!wirehair_v2::SetMixedCoefficientGeometryForTesting(mixed_geometry)) {
        return 1;
    }
#endif

    if (completion == PrecodeFailCompletion::Certified)
    {
        std::printf(
            "# precodefail: trials=%u threads=%u loss=%.17g seed=0x%llx "
            "full_payload_solve=%u\n",
            trials, threads, loss, (unsigned long long)seed,
            full_payload_solve ? 1u : 0u);
    }
    else
    {
        std::printf(
            "# precodefail: trials=%u threads=%u loss=%.17g seed=0x%llx "
            "completion=%s mixed_period=%u mixed_gf16_rows=%u "
            "mixed_geometry=%s "
            "full_payload_solve=%u\n",
            trials, threads, loss, (unsigned long long)seed,
            PrecodeFailCompletionName(completion),
            wirehair_v2::ActiveMixedCoefficientPeriod(),
            wirehair_v2::ActiveMixedGF16Rows(),
            MixedCoefficientGeometryName(
                wirehair_v2::ActiveMixedCoefficientGeometry()),
            full_payload_solve ? 1u : 0u);
    }
    std::printf(
        "N,bb,heavy_family,mix_count,overhead,trials,success,rank_fail,error,"
        "fail_rate,"
        "inact_mu,inact_max,binary_def_mu,binary_def_max,heavy_gain_mu,"
        "heavy_gain_min,heavy_shortfall,solve_ms_mu,build_ms_mu,peel_ms_mu,"
        "project_ms_mu,residual_ms_mu,backsub_ms_mu,seed_attempt,"
        "block_xors_mu,block_muladds_mu,first_rank_fail,binary_def_hist,"
        "heavy_gain_hist,failure_trials\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const uint32_t K = (uint32_t)n_value;
        const uint32_t bb = (uint32_t)bb_value;
        const wirehair_v2::SeedProfile profile =
            wirehair_v2::SelectSeedProfile(K, bb);
        const uint64_t matrix_seed = wirehair_v2::MatrixSeedFromProfile(
            profile, 0u, wirehair_v2::kMessagePrecodeSeedSalt);
        const wirehair_v2::PrecodeParams canonical_params =
            completion == PrecodeFailCompletion::Mixed ?
                wirehair_v2::MakeMixedParams(K, matrix_seed) :
                wirehair_v2::MakeCertifiedParams(K, matrix_seed);
        wirehair_v2::PacketRowConfig base_config;
        base_config.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
            profile, wirehair_v2::kMessageRecoveryRowSeedSalt);
        for (wirehair_v2::HeavyCoefficientFamily heavy_family : heavy_families)
        {
        std::map<uint32_t, PairedMixOutcomes> paired_outcomes;
        for (int mix_count_value : mix_counts)
        {
            base_config.MixCount = (uint32_t)mix_count_value;
            wirehair_v2::PrecodeParams base_params = canonical_params;
            base_params.HeavyFamily = heavy_family;
            wirehair_v2::PrecodeSystem system;
            wirehair_v2::PacketRowConfig config;
            uint32_t seed_attempt = 0u;
            const WirehairResult select_result =
                wirehair_v2::SelectSystematicConfiguration(
                    base_params, base_config, system, config, &seed_attempt);
            if (select_result != Wirehair_Success) {
                std::fprintf(stderr,
                    "precodefail seed selection failed N=%u bb=%u "
                    "heavy_family=%s mix_count=%u result=%d\n",
                    K, bb, HeavyFamilyName(heavy_family),
                    (uint32_t)mix_count_value,
                    (int)select_result);
                return 2;
            }
            // Exclude the mix count so E2E arms share the message and loss
            // stream just like the rank trials below.
            if (payload_e2e && !RunPrecodeFailPayloadE2E(
                    system, config, bb, loss,
                    seed ^ ((uint64_t)K << 32) ^ bb))
            {
                std::fprintf(stderr,
                    "precodefail payload E2E failed N=%u bb=%u "
                    "heavy_family=%s mix_count=%u\n",
                    K, bb, HeavyFamilyName(heavy_family),
                    (uint32_t)mix_count_value);
                return 2;
            }

        for (int overhead_value : overheads)
        {
            const uint32_t overhead = (uint32_t)overhead_value;
            std::vector<MatrixFailureTrial> results(trials);
            std::atomic<uint32_t> next_trial(0u);
            std::atomic<bool> cancel_workers(false);
            std::atomic<bool> worker_failed(false);
            const uint32_t worker_count = std::min(threads, trials);
            std::vector<std::thread> workers;
            ThreadJoinGuard join_guard(workers);
            workers.reserve(worker_count);
            try
            {
                for (uint32_t worker = 0; worker < worker_count; ++worker)
                {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
                    if (worker == fail_thread_launch_after)
                    {
                        throw std::system_error(
                            std::make_error_code(
                                std::errc::resource_unavailable_try_again),
                            "injected precodefail thread launch failure");
                    }
#endif
                    workers.push_back(std::thread([&, overhead]() {
                        try
                        {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
                            if (!wirehair_v2::
                                    SetMixedGF16RowsForTesting(
                                        mixed_gf16_rows))
                            {
                                throw std::runtime_error(
                                    "invalid mixed GF16 row count");
                            }
                            if (!wirehair_v2::
                                    SetMixedCoefficientPeriodForTesting(
                                        mixed_period))
                            {
                                throw std::runtime_error(
                                    "invalid mixed coefficient period");
                            }
                            if (!wirehair_v2::
                                    SetMixedCoefficientGeometryForTesting(
                                        mixed_geometry))
                            {
                                throw std::runtime_error(
                                    "invalid mixed coefficient geometry");
                            }
#endif
                            uint32_t solve_block_bytes =
                                completion == PrecodeFailCompletion::Mixed ?
                                    2u : 1u;
                            if (full_payload_solve) {
                                solve_block_bytes = bb;
                            }
                            const std::vector<uint8_t> zero(
                                solve_block_bytes, uint8_t{0});
                            for (;;)
                            {
                                if (cancel_workers.load()) {
                                    break;
                                }
                                const uint32_t trial = next_trial.fetch_add(1u);
                                if (trial >= trials) {
                                    break;
                                }
                                // Keep the delivered packet IDs paired across
                                // completion, heavy-family, and mix-count arms.
                                Rng rng(
                                    seed ^
                                    ((uint64_t)K *
                                        UINT64_C(0x9e3779b97f4a7c15)) ^
                                    ((uint64_t)bb *
                                        UINT64_C(0xbf58476d1ce4e5b9)) ^
                                    ((uint64_t)overhead *
                                        UINT64_C(0x94d049bb133111eb)) ^
                                    ((uint64_t)trial *
                                        UINT64_C(0xd6e8feb86659fd93)));
                                std::vector<wirehair_v2::SolvePacket> packets;
                                packets.reserve((size_t)K + overhead);
                                uint32_t block_id = 0u;
                                while (packets.size() < (size_t)K + overhead)
                                {
                                    const uint32_t id = block_id++;
                                    if (ShouldDrop(rng, loss)) {
                                        continue;
                                    }
                                    wirehair_v2::SolvePacket packet;
                                    packet.BlockId = id;
                                    packet.Data = zero.data();
                                    packets.push_back(packet);
                                }
                                std::vector<uint8_t> intermediate;
                                wirehair_v2::PrecodeSolveStats solve_stats;
                                MatrixFailureTrial& result = results[trial];
                                result.Result =
                                    wirehair_v2::SolvePrecodeSystem(
                                        system, config, packets,
                                        solve_block_bytes,
                                        intermediate, &solve_stats);
                                result.Inactivated =
                                    solve_stats.InactivatedColumns;
                                result.Rank = solve_stats.ResidualRank;
                                result.BinaryRank =
                                    solve_stats.BinaryResidualRank;
                                result.BlockXors = solve_stats.BlockXors;
                                result.BlockMulAdds = solve_stats.BlockMulAdds;
                                result.SolveNanoseconds =
                                    solve_stats.BuildNanoseconds +
                                    solve_stats.PeelNanoseconds +
                                    solve_stats.ProjectNanoseconds +
                                    solve_stats.ResidualNanoseconds +
                                    solve_stats.BackSubNanoseconds;
                                result.BuildNanoseconds =
                                    solve_stats.BuildNanoseconds;
                                result.PeelNanoseconds =
                                    solve_stats.PeelNanoseconds;
                                result.ProjectNanoseconds =
                                    solve_stats.ProjectNanoseconds;
                                result.ResidualNanoseconds =
                                    solve_stats.ResidualNanoseconds;
                                result.BackSubNanoseconds =
                                    solve_stats.BackSubNanoseconds;
                            }
                        }
                        catch (...) {
                            worker_failed.store(true);
                            cancel_workers.store(true);
                        }
                    }));
                }
            }
            catch (const std::system_error& error)
            {
                cancel_workers.store(true);
                std::fprintf(stderr,
                    "precodefail thread launch failed: %s\n", error.what());
                return 1;
            }
            join_guard.JoinAll();
            if (worker_failed.load()) {
                std::fprintf(stderr, "precodefail worker failed\n");
                return 2;
            }

            uint32_t successes = 0u;
            uint32_t rank_failures = 0u;
            uint32_t errors = 0u;
            uint32_t first_rank_failure = UINT32_MAX;
            uint64_t inact_sum = 0u;
            uint64_t solve_ns_sum = 0u;
            uint64_t build_ns_sum = 0u;
            uint64_t peel_ns_sum = 0u;
            uint64_t project_ns_sum = 0u;
            uint64_t residual_ns_sum = 0u;
            uint64_t backsub_ns_sum = 0u;
            uint64_t block_xors_sum = 0u;
            uint64_t block_muladds_sum = 0u;
            uint32_t inact_max = 0u;
            uint64_t binary_def_sum = 0u;
            uint32_t binary_def_max = 0u;
            uint64_t heavy_gain_sum = 0u;
            uint32_t heavy_gain_min = UINT32_MAX;
            uint32_t heavy_shortfalls = 0u;
            std::map<uint32_t, uint32_t> binary_def_hist;
            std::map<uint32_t, uint32_t> heavy_gain_hist;
            std::string failure_trials;
            for (const MatrixFailureTrial& result : results)
            {
                if (result.Result == Wirehair_Success) {
                    ++successes;
                }
                else if (result.Result == Wirehair_NeedMore) {
                    ++rank_failures;
                    if (first_rank_failure == UINT32_MAX) {
                        first_rank_failure =
                            (uint32_t)(&result - results.data());
                    }
                }
                else {
                    ++errors;
                }
                if (result.Result != Wirehair_Success)
                {
                    if (!failure_trials.empty()) {
                        failure_trials += '|';
                    }
                    failure_trials += std::to_string(
                        (uint32_t)(&result - results.data()));
                }
                inact_sum += result.Inactivated;
                solve_ns_sum += result.SolveNanoseconds;
                build_ns_sum += result.BuildNanoseconds;
                peel_ns_sum += result.PeelNanoseconds;
                project_ns_sum += result.ProjectNanoseconds;
                residual_ns_sum += result.ResidualNanoseconds;
                backsub_ns_sum += result.BackSubNanoseconds;
                block_xors_sum += result.BlockXors;
                block_muladds_sum += result.BlockMulAdds;
                inact_max = std::max(inact_max, result.Inactivated);
                const uint32_t binary_def = result.Inactivated >=
                    result.BinaryRank ?
                    result.Inactivated - result.BinaryRank : 0u;
                const uint32_t heavy_gain = result.Rank >= result.BinaryRank ?
                    result.Rank - result.BinaryRank : 0u;
                binary_def_sum += binary_def;
                binary_def_max = std::max(binary_def_max, binary_def);
                heavy_gain_sum += heavy_gain;
                heavy_gain_min = std::min(heavy_gain_min, heavy_gain);
                ++binary_def_hist[binary_def];
                ++heavy_gain_hist[heavy_gain];
                if (result.Result == Wirehair_NeedMore &&
                    binary_def <= system.Params.HeavyRows &&
                    heavy_gain < binary_def)
                {
                    ++heavy_shortfalls;
                }
            }
            if (pair_mix_counts &&
                (mix_count_value == 2 || mix_count_value == 3))
            {
                PairedMixOutcomes& paired = paired_outcomes[overhead];
                std::vector<bool>& failures = mix_count_value == 2 ?
                    paired.Mix2Failures : paired.Mix3Failures;
                uint32_t& paired_seed_attempt = mix_count_value == 2 ?
                    paired.Mix2SeedAttempt : paired.Mix3SeedAttempt;
                paired_seed_attempt = seed_attempt;
                failures.resize(trials);
                for (uint32_t trial = 0; trial < trials; ++trial) {
                    failures[trial] = results[trial].Result ==
                        Wirehair_Success ? 0u : 1u;
                }
            }
            const std::string binary_hist_text =
                CountHistogram(binary_def_hist);
            const std::string heavy_hist_text = CountHistogram(heavy_gain_hist);
            std::printf(
                "%u,%u,%s,%u,%u,%u,%u,%u,%u,%.8f,%.3f,%u,%.3f,%u,%.3f,"
                "%u,%u,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%u,%.3f,%.3f,%d,"
                "%s,%s,%s\n",
                K, bb, HeavyFamilyName(heavy_family),
                (uint32_t)mix_count_value, overhead, trials,
                successes, rank_failures, errors,
                (double)(rank_failures + errors) / trials,
                (double)inact_sum / trials,
                inact_max,
                (double)binary_def_sum / trials,
                binary_def_max,
                (double)heavy_gain_sum / trials,
                heavy_gain_min == UINT32_MAX ? 0u : heavy_gain_min,
                heavy_shortfalls,
                (double)solve_ns_sum / trials / 1000000.0,
                (double)build_ns_sum / trials / 1000000.0,
                (double)peel_ns_sum / trials / 1000000.0,
                (double)project_ns_sum / trials / 1000000.0,
                (double)residual_ns_sum / trials / 1000000.0,
                (double)backsub_ns_sum / trials / 1000000.0,
                seed_attempt,
                (double)block_xors_sum / trials,
                (double)block_muladds_sum / trials,
                first_rank_failure == UINT32_MAX ?
                    -1 : (int)first_rank_failure,
                binary_hist_text.c_str(), heavy_hist_text.c_str(),
                failure_trials.c_str());
        }
        }
        for (const std::pair<const uint32_t, PairedMixOutcomes>& paired :
            paired_outcomes)
        {
            PrintPairedMixOutcomes(
                K, bb, completion, heavy_family, paired.first, paired.second);
        }
        }
    }
    return 0;
}

int CmdSelfTest()
{
    double wilson_lower = 0.0;
    double wilson_upper = 0.0;
    Wilson95(0u, 4u, wilson_lower, wilson_upper);
    const double mcnemar_02 = ExactMcNemarP(0u, 2u);
    const double mcnemar_15 = ExactMcNemarP(1u, 5u);
    if (!std::isfinite(mcnemar_02) || !std::isfinite(mcnemar_15) ||
        !std::isfinite(wilson_lower) || !std::isfinite(wilson_upper) ||
        std::fabs(mcnemar_02 - 0.5) > 1e-12 ||
        std::fabs(mcnemar_15 - 0.21875) > 1e-12 ||
        std::fabs(wilson_lower) > 1e-12 ||
        std::fabs(wilson_upper - 0.4898908364545973) > 1e-12)
    {
        std::fprintf(stderr, "paired-statistics oracle mismatch\n");
        return 1;
    }

    const double losses[] = {0.0, 0.99};
    const uint32_t expected_drops[] = {0u, 9890u};
    for (size_t case_index = 0; case_index < 2u; ++case_index)
    {
        Rng rng(UINT64_C(0x123456789abcdef0));
        uint32_t drops = 0u;
        for (uint32_t i = 0; i < 10000u; ++i) {
            drops += ShouldDrop(rng, losses[case_index]) ? 1u : 0u;
        }
        if (drops != expected_drops[case_index])
        {
            std::fprintf(stderr,
                "loss oracle mismatch at %.17g: got %u, expected %u\n",
                losses[case_index], drops, expected_drops[case_index]);
            return 1;
        }
    }
    std::printf("loss boundary oracle: PASS\n");
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    if (wirehair_init() != Wirehair_Success) {
        std::fprintf(stderr, "wirehair_init failed\n");
        return 2;
    }

    if (argc < 2) {
        std::fprintf(stderr,
            "usage: wirehair_v2_bench compare|precodecheck|seedtable|"
            "peelcost|densecheck|densetune|densecount|densegrid|precodefail|"
            "selftest [opts]\n");
        return 1;
    }
    try
    {
        if (!std::strcmp(argv[1], "compare")) {
            return CmdCompare(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "precodecheck")) {
            return CmdPrecodeCheck(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "seedtable")) {
            return CmdSeedTable(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "peelcost")) {
            return CmdPeelCost(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "densecheck")) {
            return CmdDenseCheck(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "densetune")) {
            return CmdDenseTune(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "densecount")) {
            return CmdDenseCount(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "densegrid")) {
            return CmdDenseGrid(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "precodefail")) {
            return CmdPrecodeFail(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "selftest")) {
            if (argc != 2) {
                std::fprintf(stderr, "selftest takes no options\n");
                return 1;
            }
            return CmdSelfTest();
        }
    }
    catch (const std::bad_alloc&)
    {
        std::fprintf(stderr,
            "wirehair_v2_bench: allocation failed for requested workload\n");
        return 2;
    }
    catch (const std::length_error&)
    {
        std::fprintf(stderr,
            "wirehair_v2_bench: requested workload exceeds container limits\n");
        return 2;
    }
    std::fprintf(stderr, "unknown mode: %s\n", argv[1]);
    return 1;
}
