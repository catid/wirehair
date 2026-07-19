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
#include <fstream>
#include <iomanip>
#include <map>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <fcntl.h>
#include <io.h>
#endif

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#if defined(__linux__)
#include <sched.h>
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    uint64_t MixedEncoderJointSourceColumns;
    uint64_t MixedDecoderJointSourceColumns;
    uint64_t MixedEncoderDualSourceColumns;
    uint64_t MixedDecoderDualSourceColumns;
#endif
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
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
#endif

const char* MixedResidueScheduleName(
    wirehair_v2::MixedResidueSchedule schedule)
{
    switch (schedule)
    {
    case wirehair_v2::MixedResidueSchedule::Ramp: return "ramp";
    case wirehair_v2::MixedResidueSchedule::Hashed: return "hashed";
    default: return "constant";
    }
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool ParseMixedResidueSchedule(
    const char* text,
    wirehair_v2::MixedResidueSchedule& schedule)
{
    if (!std::strcmp(text, "constant")) {
        schedule = wirehair_v2::MixedResidueSchedule::Constant;
        return true;
    }
    if (!std::strcmp(text, "ramp")) {
        schedule = wirehair_v2::MixedResidueSchedule::Ramp;
        return true;
    }
    if (!std::strcmp(text, "hashed")) {
        schedule = wirehair_v2::MixedResidueSchedule::Hashed;
        return true;
    }
    return false;
}
#endif

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
const char* MixedResidueBucketModeName(
    wirehair_v2::MixedResidueBucketMode mode)
{
    switch (mode)
    {
    case wirehair_v2::MixedResidueBucketMode::Separate: return "separate";
    case wirehair_v2::MixedResidueBucketMode::Dual: return "dual";
    case wirehair_v2::MixedResidueBucketMode::JointDelta:
        return "joint-delta";
    default: return "auto";
    }
}

bool ParseMixedResidueBucketMode(
    const char* text,
    wirehair_v2::MixedResidueBucketMode& mode)
{
    if (!std::strcmp(text, "auto")) {
        mode = wirehair_v2::MixedResidueBucketMode::Automatic;
        return true;
    }
    if (!std::strcmp(text, "separate")) {
        mode = wirehair_v2::MixedResidueBucketMode::Separate;
        return true;
    }
    if (!std::strcmp(text, "dual")) {
        mode = wirehair_v2::MixedResidueBucketMode::Dual;
        return true;
    }
    if (!std::strcmp(text, "joint-delta")) {
        mode = wirehair_v2::MixedResidueBucketMode::JointDelta;
        return true;
    }
    return false;
}
#endif

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
    bool cache_decoder_systematic = false,
    uint32_t recovery_mix_count =
        wirehair_v2::kCertifiedPacketMixCount)
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
    encoder_options.RecoveryMixCount = recovery_mix_count;
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
             enc.Profile().V2RecoveryMixCount != recovery_mix_count ||
             dec.Profile().V2RecoveryMixCount != recovery_mix_count ||
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (const wirehair_v2::PrecodeSolveStats* stats =
            enc.PrecodeEncoderSolveStatsForTesting())
    {
        tr.MixedEncoderJointSourceColumns = stats->MixedJointSourceXors;
        tr.MixedEncoderDualSourceColumns = stats->MixedDualSourceColumns;
    }
    if (const wirehair_v2::PrecodeSolveStats* stats =
            dec.PrecodeSolveStatsForTesting())
    {
        tr.MixedDecoderJointSourceColumns = stats->MixedJointSourceXors;
        tr.MixedDecoderDualSourceColumns = stats->MixedDualSourceColumns;
    }
#endif
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
    bool mixed_residue_hash_keyed = false;
    bool mixed_independent_extension_residues = false;
    uint32_t mixed_mix_count = wirehair_v2::kCertifiedPacketMixCount;
    uint32_t packet_row_seed_multiplier = 1u;
    bool packet_row_seed_avalanche = false;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    bool mixed_mix_count_explicit = false;
    uint32_t mixed_period = wirehair_v2::kMixedCoefficientPeriod;
    bool mixed_period_explicit = false;
    uint32_t mixed_gf256_rows = wirehair_v2::kMixedGF256Rows;
    bool mixed_gf256_rows_explicit = false;
    uint32_t mixed_grouped_gf256_rows = 0u;
    bool mixed_grouped_gf256_rows_explicit = false;
    uint32_t mixed_gf16_rows = wirehair_v2::kMixedGF16Rows;
    bool mixed_gf16_rows_explicit = false;
    wirehair_v2::MixedCoefficientGeometry mixed_geometry =
        wirehair_v2::MixedCoefficientGeometry::FrozenPowerX;
    bool mixed_geometry_explicit = false;
    uint32_t mixed_residue_skew = 0u;
    bool mixed_residue_skew_explicit = false;
    wirehair_v2::MixedResidueSchedule mixed_residue_schedule =
        wirehair_v2::MixedResidueSchedule::Constant;
    bool mixed_residue_schedule_explicit = false;
    uint32_t mixed_residue_hash_seed = 0u;
    bool mixed_residue_hash_seed_explicit = false;
    wirehair_v2::MixedResidueBucketMode mixed_residue_bucket_mode =
        wirehair_v2::MixedResidueBucketMode::Automatic;
    bool mixed_residue_bucket_mode_explicit = false;
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
        else if (!std::strcmp(argv[i], "--mixed-mix-count")) {
            if (!TakeArg(
                    "compare", "--mixed-mix-count", argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-mix-count", value, mixed_mix_count))
            {
                return 1;
            }
            mixed_mix_count_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-gf16-rows")) {
            if (!TakeArg(
                    "compare", "--mixed-gf16-rows", argc, argv, i, value) ||
                !ParseU32Arg("--mixed-gf16-rows", value, mixed_gf16_rows))
            {
                return 1;
            }
            mixed_gf16_rows_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-gf256-rows")) {
            if (!TakeArg(
                    "compare", "--mixed-gf256-rows", argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-gf256-rows", value, mixed_gf256_rows))
            {
                return 1;
            }
            mixed_gf256_rows_explicit = true;
        }
        else if (!std::strcmp(
                     argv[i], "--mixed-grouped-gf256-rows"))
        {
            if (!TakeArg(
                    "compare", "--mixed-grouped-gf256-rows",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-grouped-gf256-rows", value,
                    mixed_grouped_gf256_rows))
            {
                return 1;
            }
            mixed_grouped_gf256_rows_explicit = true;
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
        else if (!std::strcmp(argv[i], "--mixed-residue-skew")) {
            if (!TakeArg(
                    "compare", "--mixed-residue-skew",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-residue-skew", value, mixed_residue_skew))
            {
                return 1;
            }
            mixed_residue_skew_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-schedule")) {
            if (!TakeArg(
                    "compare", "--mixed-residue-schedule",
                    argc, argv, i, value) ||
                !ParseMixedResidueSchedule(value, mixed_residue_schedule))
            {
                std::fprintf(stderr,
                    "compare unknown --mixed-residue-schedule token %s "
                    "(expected constant, ramp, or hashed)\n",
                    value ? value : "");
                return 1;
            }
            mixed_residue_schedule_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-hash-seed")) {
            if (!TakeArg(
                    "compare", "--mixed-residue-hash-seed",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-residue-hash-seed", value,
                    mixed_residue_hash_seed))
            {
                return 1;
            }
            mixed_residue_hash_seed_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-hash-keyed")) {
            mixed_residue_hash_keyed = true;
        }
        else if (!std::strcmp(
                     argv[i],
                     "--mixed-independent-extension-residues"))
        {
            mixed_independent_extension_residues = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-buckets"))
        {
            if (!TakeArg(
                    "compare", "--mixed-residue-buckets",
                    argc, argv, i, value) ||
                !ParseMixedResidueBucketMode(
                    value, mixed_residue_bucket_mode))
            {
                std::fprintf(stderr,
                    "compare unknown --mixed-residue-buckets token %s "
                    "(expected auto, separate, dual, or joint-delta)\n",
                    value ? value : "");
                return 1;
            }
            mixed_residue_bucket_mode_explicit = true;
        }
        else if (!std::strcmp(
                     argv[i], "--packet-row-seed-multiplier"))
        {
            if (!TakeArg(
                    "compare", "--packet-row-seed-multiplier",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--packet-row-seed-multiplier", value,
                    packet_row_seed_multiplier))
            {
                return 1;
            }
        }
        else if (!std::strcmp(
                     argv[i], "--packet-row-seed-avalanche"))
        {
            packet_row_seed_avalanche = true;
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
    else if (mixed_mix_count_explicit || mixed_period_explicit ||
             mixed_geometry_explicit ||
             mixed_gf256_rows_explicit ||
             mixed_grouped_gf256_rows_explicit ||
             mixed_gf16_rows_explicit ||
             mixed_residue_skew_explicit ||
             mixed_residue_schedule_explicit ||
             mixed_residue_hash_seed_explicit ||
             mixed_residue_hash_keyed ||
             mixed_independent_extension_residues ||
             mixed_residue_bucket_mode_explicit)
    {
        std::fprintf(stderr,
            "compare mixed experiment flags require a mixed precode "
            "profile\n");
        return 1;
    }
    if (mixed_mix_count < 2u ||
        mixed_mix_count > wirehair_v2::kCertifiedPacketMixCount)
    {
        std::fprintf(stderr,
            "compare --mixed-mix-count must be in [2,%u]\n",
            wirehair_v2::kCertifiedPacketMixCount);
        return 1;
    }
    if (!wirehair_v2::SetMixedCoefficientGeometryForTesting(mixed_geometry)) {
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
            wirehair_v2::ActiveMixedGF256Rows() +
                wirehair_v2::ActiveMixedGF16Rows(),
            wirehair_v2::kMixedCoefficientPeriod);
        return 1;
    }
    if (!wirehair_v2::SetMixedGF256RowsForTesting(mixed_gf256_rows))
    {
        std::fprintf(stderr,
            "compare --mixed-gf256-rows must be in [%u,%u], fit the "
            "active period, use shared-x for an extra row, and use the "
            "validated 12+4 geometry for twelve rows\n",
            wirehair_v2::kMixedGF256Rows,
            wirehair_v2::kMixedGF256RowsMax);
        return 1;
    }
    if (!wirehair_v2::SetMixedResidueSkewForTesting(mixed_residue_skew)) {
        std::fprintf(stderr,
            "compare --mixed-residue-skew must be a corner-preserving "
            "shared-x skew in [0,P-H]\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedResidueScheduleForTesting(
            mixed_residue_schedule))
    {
        std::fprintf(stderr,
            "compare nonconstant --mixed-residue-schedule requires shared-x, "
            "P>H, and zero constant skew\n");
        return 1;
    }
    if ((mixed_residue_hash_seed_explicit || mixed_residue_hash_keyed) &&
        mixed_residue_schedule != wirehair_v2::MixedResidueSchedule::Hashed)
    {
        std::fprintf(stderr,
            "compare residue hash seed/keying requires hashed "
            "--mixed-residue-schedule\n");
        return 1;
    }
    wirehair_v2::SetMixedResidueHashSeedForTesting(
        mixed_residue_hash_seed);
    if (!wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(
            mixed_independent_extension_residues))
    {
        std::fprintf(stderr,
            "compare independent extension residues require "
            "shared-x hashed scheduling with P>H\n");
        return 1;
    }
    if (mixed_residue_bucket_mode !=
            wirehair_v2::MixedResidueBucketMode::Automatic &&
        !mixed_independent_extension_residues &&
        mixed_grouped_gf256_rows == 0u)
    {
        std::fprintf(stderr,
            "compare explicit --mixed-residue-buckets requires "
            "--mixed-independent-extension-residues or nonzero "
            "--mixed-grouped-gf256-rows\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedResidueBucketModeForTesting(
            mixed_residue_bucket_mode))
    {
        return 1;
    }
    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(
            packet_row_seed_multiplier))
    {
        std::fprintf(stderr,
            "compare --packet-row-seed-multiplier must be odd and "
            "nonzero\n");
        return 1;
    }
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(
        packet_row_seed_avalanche);
    // Every earlier mixed configuration setter may clear schedule experiments.
    // Activate the grouped suffix only after the complete thread state is set.
    if (!wirehair_v2::SetMixedGroupedGF256RowsForTesting(
            mixed_grouped_gf256_rows))
    {
        std::fprintf(stderr,
            "compare --mixed-grouped-gf256-rows must be in [0,9]; "
            "nonzero grouping requires shared-x constant-A 10+2 geometry, "
            "P>H, and no independent extension residues\n");
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
        "schedule_seed=0x%llx mixed_mix_count=%u mixed_period=%u "
        "mixed_gf256_rows=%u mixed_gf16_rows=%u "
        "mixed_geometry=%s mixed_residue_skew=%u "
        "mixed_residue_schedule=%s mixed_residue_hash_seed=0x%x "
        "mixed_residue_hash_keyed=%u "
        "mixed_independent_extension_residues=%u "
        "mixed_grouped_gf256_rows=%u "
        "mixed_grouped_gf256_hash_seed=0x%x "
        "mixed_grouped_final_h_a_columns=%u "
        "mixed_residue_buckets_requested=%s "
        "packet_row_seed_multiplier=0x%x "
        "packet_row_seed_avalanche=%u "
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
        mixed_mix_count,
        wirehair_v2::ActiveMixedCoefficientPeriod(),
        wirehair_v2::ActiveMixedGF256Rows(),
        wirehair_v2::ActiveMixedGF16Rows(),
        MixedCoefficientGeometryName(
            wirehair_v2::ActiveMixedCoefficientGeometry()),
        wirehair_v2::ActiveMixedResidueSkew(),
        MixedResidueScheduleName(
            wirehair_v2::ActiveMixedResidueSchedule()),
        wirehair_v2::ActiveMixedResidueHashSeed(),
        mixed_residue_hash_keyed ? 1u : 0u,
        mixed_independent_extension_residues ? 1u : 0u,
        wirehair_v2::ActiveMixedGroupedGF256Rows(),
        wirehair_v2::ActiveMixedGroupedGF256HashSeed(),
        wirehair_v2::ActiveMixedGroupedGF256Rows() != 0u ?
            wirehair_v2::ActiveMixedGF256Rows() +
                wirehair_v2::ActiveMixedGF16Rows() : 0u,
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        MixedResidueBucketModeName(mixed_residue_bucket_mode),
#else
        "auto",
#endif
        packet_row_seed_multiplier,
        packet_row_seed_avalanche ? 1u : 0u);
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            if (mixed_residue_hash_keyed)
            {
                uint32_t selected_hash_seed = 0u;
                if (!wirehair_v2::
                        SelectFullCycleMixedResidueKeyedSeedForTesting(
                        mixed_residue_hash_seed, N, selected_hash_seed))
                {
                    std::fprintf(stderr,
                        "compare could not select a full-cycle keyed "
                        "residue hash seed for N=%u\n",
                        N);
                    return 1;
                }
            }
            if (!wirehair_v2::
                    SetMixedIndependentExtensionResiduesForTesting(
                        mixed_independent_extension_residues))
            {
                return 1;
            }
            // Keyed seed selection and independent-schedule replay above may
            // clear schedule experiments.  Keep grouped C last for every
            // trial, not only for the command's initial TLS configuration.
            if (!wirehair_v2::SetMixedGroupedGF256RowsForTesting(
                    mixed_grouped_gf256_rows))
            {
                return 1;
            }
#endif
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
                        completion, false, false,
                        completion ==
                                wirehair_v2::CompletionField::MixedGF256GF16 ?
                            mixed_mix_count :
                            wirehair_v2::kCertifiedPacketMixCount);
                }
                if (include_precode_cache) {
                    cache_result = RunV2PrecodeTrial(
                        N, block_bytes, loss, trial_seed, &packet_schedule,
                        completion,
                        cache_encoder_source, cache_decoder_systematic,
                        completion ==
                                wirehair_v2::CompletionField::MixedGF256GF16 ?
                            mixed_mix_count :
                            wirehair_v2::kCertifiedPacketMixCount);
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            const auto validate_effective_bucket_mode = [&](
                const TrialResult& result, const char* arm) -> bool
            {
                if (mixed_residue_bucket_mode ==
                        wirehair_v2::MixedResidueBucketMode::Automatic ||
                    mixed_residue_bucket_mode ==
                        wirehair_v2::MixedResidueBucketMode::Separate)
                {
                    return true;
                }
                const bool joint = mixed_residue_bucket_mode ==
                    wirehair_v2::MixedResidueBucketMode::JointDelta;
                const bool encoder_used = joint ?
                    result.MixedEncoderJointSourceColumns != 0u :
                    result.MixedEncoderDualSourceColumns != 0u;
                const bool decoder_used = joint ?
                    result.MixedDecoderJointSourceColumns != 0u :
                    result.MixedDecoderDualSourceColumns != 0u;
                if (encoder_used && decoder_used) return true;
                std::fprintf(stderr,
                    "compare requested mixed residue bucket mode %s but "
                    "arm %s fell back (N=%u bb=%u encoder_used=%u "
                    "decoder_used=%u)\n",
                    MixedResidueBucketModeName(mixed_residue_bucket_mode),
                    arm, N, block_bytes,
                    encoder_used ? 1u : 0u, decoder_used ? 1u : 0u);
                return false;
            };
            if (include_mixed &&
                ((include_precode && !validate_effective_bucket_mode(
                    mixed_precode_result, "precode")) ||
                 (include_precode_cache && !validate_effective_bucket_mode(
                    mixed_cache_result, "precode-cache"))))
            {
                return 1;
            }
#endif
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
    uint64_t MixedJointSourceXors = 0u;
    uint64_t MixedJointMarginalXors = 0u;
    uint64_t MixedJointMarginalCopies = 0u;
    uint64_t MixedJointScratchBytes = 0u;
    uint32_t MixedJointActiveDeltas = 0u;
    uint64_t MixedDualSourceColumns = 0u;
};

uint64_t TotalSolveNanoseconds(
    const wirehair_v2::PrecodeSolveStats& stats)
{
    return stats.BuildNanoseconds + stats.PeelNanoseconds +
        stats.ProjectNanoseconds + stats.ResidualNanoseconds +
        stats.BackSubNanoseconds;
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool SameNonTimingSolveStats(
    const wirehair_v2::PrecodeSolveStats& a,
    const wirehair_v2::PrecodeSolveStats& b)
{
    return a.PacketRows == b.PacketRows &&
        a.PeeledColumns == b.PeeledColumns &&
        a.InactivatedColumns == b.InactivatedColumns &&
        a.ResidualRows == b.ResidualRows &&
        a.ResidualRank == b.ResidualRank &&
        a.BinaryResidualRank == b.BinaryResidualRank &&
        a.BinaryRowReferences == b.BinaryRowReferences &&
        a.BinaryRowStorageBytes == b.BinaryRowStorageBytes &&
        a.BinaryAdjacencyStorageBytes == b.BinaryAdjacencyStorageBytes &&
        a.BinaryRowStorageAllocations == b.BinaryRowStorageAllocations &&
        a.BinaryAdjacencyStorageAllocations ==
            b.BinaryAdjacencyStorageAllocations &&
        a.BlockXors == b.BlockXors &&
        a.BlockMulAdds == b.BlockMulAdds &&
        a.MixedJointSourceXors == b.MixedJointSourceXors &&
        a.MixedJointMarginalXors == b.MixedJointMarginalXors &&
        a.MixedJointMarginalCopies == b.MixedJointMarginalCopies &&
        a.MixedJointScratchBytes == b.MixedJointScratchBytes &&
        a.MixedJointActiveDeltas == b.MixedJointActiveDeltas &&
        a.MixedDualSourceColumns == b.MixedDualSourceColumns &&
        a.PacketSeedAttempt == b.PacketSeedAttempt;
}

class MixedNullWitnessScope
{
public:
    explicit MixedNullWitnessScope(
        wirehair_v2::MixedNullWitnessDiagnostic* diagnostic)
        : Active(diagnostic != nullptr)
    {
        if (Active) {
            wirehair_v2::SetMixedNullWitnessDiagnosticForTesting(diagnostic);
        }
    }

    ~MixedNullWitnessScope()
    {
        Disable();
    }

    void Disable()
    {
        if (Active) {
            wirehair_v2::SetMixedNullWitnessDiagnosticForTesting(nullptr);
            Active = false;
        }
    }

private:
    bool Active;
};

enum class MixedNullReplayStatus
{
    None,
    Captured,
    Skipped,
    Error
};

static const int kMixedNullReplayInternalErrorExitCode = 3;

int MixedNullReplayExitCode(MixedNullReplayStatus status)
{
    return status == MixedNullReplayStatus::Error ?
        kMixedNullReplayInternalErrorExitCode : 0;
}

const char* MixedNullReplayStatusName(MixedNullReplayStatus status)
{
    switch (status)
    {
    case MixedNullReplayStatus::None: return "none";
    case MixedNullReplayStatus::Captured: return "captured";
    case MixedNullReplayStatus::Skipped: return "skipped";
    default: return "error";
    }
}

struct MixedNullBucket
{
    uint32_t Count = 0u;
    uint16_t Sum = 0u;
};

struct MixedNullBucketSummary
{
    uint32_t Occupied = 0u;
    uint32_t Cancelled = 0u;
    uint32_t CancelledTerms = 0u;
    uint32_t Maximum = 0u;
    uint64_t Hash = 0u;
    std::string Top;
};

void HashMixedNullByte(uint8_t value, uint64_t& hash)
{
    hash ^= value;
    hash *= UINT64_C(0x100000001b3);
}

void HashMixedNullU32(uint32_t value, uint64_t& hash)
{
    for (uint32_t byte = 0u; byte < 4u; ++byte) {
        HashMixedNullByte((uint8_t)(value >> (8u * byte)), hash);
    }
}

MixedNullBucketSummary SummarizeMixedNullBuckets(
    uint8_t domain,
    uint32_t row,
    uint32_t period,
    const std::vector<MixedNullBucket>& buckets,
    bool include_top)
{
    MixedNullBucketSummary summary;
    summary.Hash = UINT64_C(0xcbf29ce484222325);
    HashMixedNullByte(domain, summary.Hash);
    HashMixedNullU32(row, summary.Hash);
    HashMixedNullU32(period, summary.Hash);
    std::vector<uint32_t> occupied;
    occupied.reserve(include_top ? period : 0u);
    for (uint32_t index = 0u; index < (uint32_t)buckets.size(); ++index)
    {
        const MixedNullBucket& bucket = buckets[index];
        HashMixedNullU32(index, summary.Hash);
        HashMixedNullU32(bucket.Count, summary.Hash);
        HashMixedNullByte((uint8_t)bucket.Sum, summary.Hash);
        HashMixedNullByte((uint8_t)(bucket.Sum >> 8), summary.Hash);
        if (bucket.Count == 0u) continue;
        ++summary.Occupied;
        summary.Maximum = std::max(summary.Maximum, bucket.Count);
        if (bucket.Sum == 0u) {
            ++summary.Cancelled;
            summary.CancelledTerms += bucket.Count;
        }
        if (include_top) occupied.push_back(index);
    }
    std::sort(occupied.begin(), occupied.end(), [&](uint32_t a, uint32_t b) {
        if (buckets[a].Count != buckets[b].Count) {
            return buckets[a].Count > buckets[b].Count;
        }
        return a < b;
    });
    const size_t top_count = std::min<size_t>(8u, occupied.size());
    for (size_t i = 0u; i < top_count; ++i)
    {
        const uint32_t residue = occupied[i];
        const MixedNullBucket& bucket = buckets[residue];
        char item[48];
        std::snprintf(
            item, sizeof(item), "%s%u:%u:%04x",
            i == 0u ? "" : "|", residue, bucket.Count, bucket.Sum);
        summary.Top += item;
    }
    if (summary.Top.empty()) summary.Top = "-";
    return summary;
}

bool BuildMixedNullClassification(
    uint32_t source_count,
    const wirehair_v2::MixedNullWitnessDiagnostic& witness,
    std::vector<std::string>& lines)
{
    const uint32_t L = witness.ColumnCount;
    const uint32_t d = witness.KernelDimension;
    const uint32_t period = wirehair_v2::ActiveMixedCoefficientPeriod();
    const uint32_t subfield_rows = wirehair_v2::ActiveMixedGF256Rows();
    const uint32_t extension_rows = wirehair_v2::ActiveMixedGF16Rows();
    const uint32_t heavy_rows = subfield_rows + extension_rows;
    const uint32_t grouped_gf256_rows =
        wirehair_v2::ActiveMixedGroupedGF256Rows();
    const bool basis_size_overflow = d != 0u &&
        (size_t)L > std::numeric_limits<size_t>::max() / d;
    if (witness.Status != wirehair_v2::MixedNullWitnessStatus::Captured ||
        d == 0u ||
        d > wirehair_v2::kMaxMixedNullWitnessQuotientColumns ||
        period == 0u || period > 244u ||
        grouped_gf256_rows > subfield_rows || L < heavy_rows ||
        basis_size_overflow || source_count > L ||
        witness.InactiveMask.size() != L ||
        witness.CanonicalBasis.size() != (size_t)d * L)
    {
        return false;
    }
    uint32_t inactive_count = 0u;
    for (uint8_t inactive : witness.InactiveMask)
    {
        if (inactive > 1u) return false;
        inactive_count += inactive;
    }
    if (inactive_count != witness.InactiveCount) return false;
    const uint32_t first_grouped_gf256_row =
        subfield_rows - grouped_gf256_rows;
    const uint32_t first_heavy_column = L - heavy_rows;
    std::vector<MixedNullBucket> subfield(period), extension(period);
    std::vector<MixedNullBucket> grouped_subfield(period);
    std::vector<MixedNullBucket> joint((size_t)period * period);
    std::vector<MixedNullBucket> grouped_joint((size_t)period * period);
    lines.clear();
    lines.reserve(d);
    for (uint32_t row = 0u; row < d; ++row)
    {
        std::fill(subfield.begin(), subfield.end(), MixedNullBucket{});
        std::fill(extension.begin(), extension.end(), MixedNullBucket{});
        std::fill(
            grouped_subfield.begin(), grouped_subfield.end(),
            MixedNullBucket{});
        std::fill(joint.begin(), joint.end(), MixedNullBucket{});
        std::fill(
            grouped_joint.begin(), grouped_joint.end(), MixedNullBucket{});
        uint32_t parts[4] = {};
        uint32_t gf256_values = 0u;
        uint32_t gf16_values = 0u;
        const uint16_t* vector = witness.CanonicalBasis.data() +
            (size_t)row * L;
        for (uint32_t column = 0u; column < L; ++column)
        {
            const uint16_t value = vector[column];
            if (value == 0u) continue;
            const bool source = column < source_count;
            const bool inactive = witness.InactiveMask[column] != 0u;
            ++parts[(source ? 0u : 2u) + (inactive ? 1u : 0u)];
            if ((value >> 8) == 0u) ++gf256_values;
            else ++gf16_values;
            const uint32_t sf =
                wirehair_v2::ActiveMixedCoefficientResidue(column);
            const uint32_t ex =
                wirehair_v2::ActiveMixedExtensionCoefficientResidue(column);
            // Grouped suffix rows use C only before the final H completion
            // columns.  The accessor intentionally maps that final corner
            // back to canonical A, matching the encoder and direct syndrome
            // verifier rather than treating all L columns as C-scheduled.
            const uint32_t grouped_sf = wirehair_v2::
                ActiveMixedGroupedGF256CoefficientResidue(
                    column, first_heavy_column);
            if (sf >= period || ex >= period || grouped_sf >= period) {
                return false;
            }
            ++subfield[sf].Count;
            subfield[sf].Sum ^= value;
            ++extension[ex].Count;
            extension[ex].Sum ^= value;
            ++grouped_subfield[grouped_sf].Count;
            grouped_subfield[grouped_sf].Sum ^= value;
            MixedNullBucket& pair = joint[(size_t)sf * period + ex];
            ++pair.Count;
            pair.Sum ^= value;
            MixedNullBucket& grouped_pair =
                grouped_joint[(size_t)sf * period + grouped_sf];
            ++grouped_pair.Count;
            grouped_pair.Sum ^= value;
        }
        const uint32_t support =
            parts[0] + parts[1] + parts[2] + parts[3];
        const MixedNullBucketSummary sf = SummarizeMixedNullBuckets(
            UINT8_C(0x53), row, period, subfield, true);
        const MixedNullBucketSummary ex = SummarizeMixedNullBuckets(
            UINT8_C(0x45), row, period, extension, true);
        const MixedNullBucketSummary pair = SummarizeMixedNullBuckets(
            UINT8_C(0x4a), row, period, joint, false);
        const MixedNullBucketSummary grouped_sf = SummarizeMixedNullBuckets(
            UINT8_C(0x43), row, period, grouped_subfield, true);
        const MixedNullBucketSummary grouped_pair =
            SummarizeMixedNullBuckets(
                UINT8_C(0x4b), row, period, grouped_joint, false);
        if (support != gf256_values + gf16_values ||
            pair.Occupied < std::max(sf.Occupied, ex.Occupied) ||
            pair.Occupied > support ||
            grouped_pair.Occupied <
                std::max(sf.Occupied, grouped_sf.Occupied) ||
            grouped_pair.Occupied > support)
        {
            return false;
        }
        char line[2048];
        const int written = std::snprintf(
            line, sizeof(line),
            "# mixed_null_row,v=1,row=%u,nz=%u,source=%u,precode=%u,"
            "source_peeled=%u,source_inactive=%u,"
            "precode_peeled=%u,precode_inactive=%u,"
            "gf256=%u,gf16=%u,sf_occ=%u,sf_cancel=%u,"
            "sf_cancel_terms=%u,sf_max=%u,sf_hash=%016llx,"
            "ex_occ=%u,ex_cancel=%u,ex_cancel_terms=%u,ex_max=%u,"
            "ex_hash=%016llx,pair_occ=%u,pair_cancel=%u,"
            "pair_cancel_terms=%u,pair_max=%u,pair_hash=%016llx,"
            "sf_top=%s,ex_top=%s,grouped_gf256_rows=%u,"
            "grouped_first_row=%u,grouped_final_h_a_columns=%u,"
            "c_occ=%u,c_cancel=%u,c_cancel_terms=%u,c_max=%u,"
            "c_hash=%016llx,ac_pair_occ=%u,ac_pair_cancel=%u,"
            "ac_pair_cancel_terms=%u,ac_pair_max=%u,"
            "ac_pair_hash=%016llx,c_top=%s",
            row, support, parts[0] + parts[1], parts[2] + parts[3],
            parts[0], parts[1], parts[2], parts[3],
            gf256_values, gf16_values,
            sf.Occupied, sf.Cancelled, sf.CancelledTerms, sf.Maximum,
            (unsigned long long)sf.Hash,
            ex.Occupied, ex.Cancelled, ex.CancelledTerms, ex.Maximum,
            (unsigned long long)ex.Hash,
            pair.Occupied, pair.Cancelled, pair.CancelledTerms, pair.Maximum,
            (unsigned long long)pair.Hash,
            sf.Top.c_str(), ex.Top.c_str(),
            grouped_gf256_rows, first_grouped_gf256_row,
            grouped_gf256_rows != 0u ? heavy_rows : 0u,
            grouped_sf.Occupied, grouped_sf.Cancelled,
            grouped_sf.CancelledTerms, grouped_sf.Maximum,
            (unsigned long long)grouped_sf.Hash,
            grouped_pair.Occupied, grouped_pair.Cancelled,
            grouped_pair.CancelledTerms, grouped_pair.Maximum,
            (unsigned long long)grouped_pair.Hash,
            grouped_sf.Top.c_str());
        if (written < 0 || (size_t)written >= sizeof(line)) return false;
        lines.push_back(line);
    }
    return true;
}
#endif

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

uint32_t NormalizedH15V1PacketPeelSeedXor(uint32_t K)
{
    // Offline hard-loss tuning for the normalized 1280-byte seed profile.
    // Unlisted K values deliberately retain the base packet graph.
    switch (K)
    {
    case 4u:  return 57u;
    case 5u:  return 105u;
    case 6u:  return 23u;
    case 7u:  return 83u;
    case 9u:  return 51u;
    case 11u: return 168u;
    case 12u: return 250u;
    case 13u: return 123u;
    case 14u: return 157u;
    case 15u: return 55u;
    case 17u: return 185u;
    case 18u: return 37u;
    case 19u: return 194u;
    case 22u: return 204u;
    case 25u: return 104u;
    case 28u: return 122u;
    case 29u: return 242u;
    case 31u: return 129u;
    case 32u: return 56u;
    case 33u: return 110u;
    case 37u: return 172u;
    case 39u: return 148u;
    case 41u: return 125u;
    default:  return 0u;
    }
}

uint32_t NormalizedH15V2PacketPeelSeedXor(uint32_t K)
{
    // Fresh-seed and cross-payload holdouts for large-K resonances.  Keep v1
    // immutable so previously recorded benchmark commands remain replayable.
    switch (K)
    {
    case 1683u:  return 19u;
    case 15182u: return 98u;
    case 21394u: return 26u;
    case 24432u: return 75u;
    case 34207u: return 213u;
    case 62039u: return 2u;
    default:     return NormalizedH15V1PacketPeelSeedXor(K);
    }
}

uint32_t NormalizedH15V3PacketPeelSeedXor(uint32_t K)
{
    // Independent hard-loss holdouts for residual v2 all-K hotspots.  Keep
    // both earlier tables immutable so recorded experiments remain replayable.
    switch (K)
    {
    case 10u:    return 139u;
    case 20u:    return 140u;
    case 11414u: return 86u;
    case 48567u: return 209u;
    case 49312u: return 52u;
    case 49842u: return 188u;
    case 50281u: return 121u;
    case 51375u: return 192u;
    case 53503u: return 238u;
    default:     return NormalizedH15V2PacketPeelSeedXor(K);
    }
}

uint32_t NormalizedH15V4PacketPeelSeedXor(uint32_t K)
{
    // Discovery-selected rank-one salts, independently validated by a frozen
    // hard-loss holdout over recurrent v3 hotspots.  Preserve every earlier
    // table entry and leave all other packet graphs unchanged.
    switch (K)
    {
    case 16u:    return 6u;
    case 39559u: return 60u;
    case 40831u: return 179u;
    case 43742u: return 27u;
    case 43751u: return 99u;
    case 45168u: return 108u;
    case 45464u: return 34u;
    case 45857u: return 49u;
    case 45903u: return 58u;
    case 46296u: return 4u;
    case 46606u: return 235u;
    case 46933u: return 106u;
    case 47029u: return 117u;
    case 47105u: return 81u;
    case 47307u: return 178u;
    case 48231u: return 225u;
    case 48311u: return 122u;
    case 48466u: return 87u;
    case 49124u: return 237u;
    case 49412u: return 173u;
    case 49486u: return 142u;
    case 49627u: return 172u;
    case 49727u: return 143u;
    case 49865u: return 255u;
    case 50689u: return 142u;
    case 50885u: return 63u;
    case 50899u: return 12u;
    case 51494u: return 208u;
    case 52935u: return 8u;
    case 53613u: return 30u;
    case 53697u: return 204u;
    case 53804u: return 169u;
    default:     return NormalizedH15V3PacketPeelSeedXor(K);
    }
}

enum class PacketPeelSeedTable : uint32_t
{
    None = 0,
    NormalizedH15V1 = 1,
    NormalizedH15V2 = 2,
    NormalizedH15V3 = 3,
    NormalizedH15V4 = 4
};

const char* PacketPeelSeedTableName(PacketPeelSeedTable table)
{
    switch (table)
    {
    case PacketPeelSeedTable::NormalizedH15V1:
        return "normalized-h15-v1";
    case PacketPeelSeedTable::NormalizedH15V2:
        return "normalized-h15-v2";
    case PacketPeelSeedTable::NormalizedH15V3:
        return "normalized-h15-v3";
    case PacketPeelSeedTable::NormalizedH15V4:
        return "normalized-h15-v4";
    default:
        return "none";
    }
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \
    !defined(WIREHAIR_V2_BENCH_DISABLE_PREFERRED_ATTEMPT)

enum class PreferredAttemptMode
{
    Control,
    Route,
    Candidate,
    Paired
};

const char* PreferredAttemptModeName(PreferredAttemptMode mode)
{
    switch (mode)
    {
    case PreferredAttemptMode::Control:   return "control";
    case PreferredAttemptMode::Route:     return "route";
    case PreferredAttemptMode::Candidate: return "candidate";
    case PreferredAttemptMode::Paired:    return "paired";
    }
    return "unknown";
}

bool ParsePreferredAttemptMode(
    const char* text,
    PreferredAttemptMode& mode)
{
    if (!std::strcmp(text, "control")) {
        mode = PreferredAttemptMode::Control;
    }
    else if (!std::strcmp(text, "route")) {
        mode = PreferredAttemptMode::Route;
    }
    else if (!std::strcmp(text, "candidate")) {
        mode = PreferredAttemptMode::Candidate;
    }
    else if (!std::strcmp(text, "paired")) {
        mode = PreferredAttemptMode::Paired;
    }
    else {
        return false;
    }
    return true;
}

bool MarkPreferredOptionOnce(
    const char* option,
    bool& seen)
{
    if (seen) {
        std::fprintf(stderr,
            "preferredattempt %s specified more than once\n", option);
        return false;
    }
    seen = true;
    return true;
}

bool IsCanonicalPreferredIntList(
    const std::string& text,
    const std::vector<int>& values)
{
    std::ostringstream canonical;
    for (size_t i = 0u; i < values.size(); ++i) {
        if (i != 0u) canonical << ',';
        canonical << values[i];
    }
    return text == canonical.str();
}

typedef std::pair<uint32_t, uint32_t> PreferredAttemptKey;
typedef std::map<PreferredAttemptKey, std::vector<uint32_t> >
    PreferredAttemptMap;

// Parse a compact, shell-independent map used only by the frozen experiment
// harness.  Each record is K@block_bytes=p[,p...] and records are separated by
// '|'.  Candidate order is significant to the external survivor ledger.
bool ParsePreferredAttemptMap(
    const std::string& text,
    PreferredAttemptMap& output)
{
    output.clear();
    if (text == "none") {
        return true;
    }
    size_t record_pos = 0u;
    while (record_pos <= text.size())
    {
        const size_t record_end = text.find('|', record_pos);
        const std::string record = text.substr(
            record_pos,
            record_end == std::string::npos ?
                std::string::npos : record_end - record_pos);
        const size_t at = record.find('@');
        const size_t equals = record.find('=');
        if (record.empty() || at == std::string::npos ||
            equals == std::string::npos || at == 0u || equals <= at + 1u ||
            record.find('@', at + 1u) != std::string::npos ||
            record.find('=', equals + 1u) != std::string::npos)
        {
            output.clear();
            return false;
        }
        uint32_t K = 0u;
        uint32_t block_bytes = 0u;
        if (!ParseU32Scalar(record.substr(0u, at).c_str(), K) ||
            !ParseU32Scalar(
                record.substr(at + 1u, equals - at - 1u).c_str(),
                block_bytes))
        {
            output.clear();
            return false;
        }
        std::vector<uint32_t> attempts;
        const bool explicit_control = record.substr(equals + 1u) == "none";
        if (!explicit_control)
        {
            size_t attempt_pos = equals + 1u;
            while (attempt_pos <= record.size())
            {
                const size_t attempt_end = record.find(',', attempt_pos);
                const std::string token = record.substr(
                    attempt_pos,
                    attempt_end == std::string::npos ?
                        std::string::npos : attempt_end - attempt_pos);
                uint32_t attempt = 0u;
                if (!ParseU32Scalar(token.c_str(), attempt) ||
                    attempt >= wirehair_v2::kMaxPacketSeedAttempts ||
                    std::find(attempts.begin(), attempts.end(), attempt) !=
                        attempts.end())
                {
                    output.clear();
                    return false;
                }
                attempts.push_back(attempt);
                if (attempts.size() > 32u) {
                    output.clear();
                    return false;
                }
                if (attempt_end == std::string::npos) {
                    break;
                }
                attempt_pos = attempt_end + 1u;
            }
        }
        if ((!explicit_control && attempts.empty()) ||
            !output.insert(std::make_pair(
                PreferredAttemptKey(K, block_bytes), attempts)).second)
        {
            output.clear();
            return false;
        }
        if (record_end == std::string::npos) {
            break;
        }
        record_pos = record_end + 1u;
    }
    return true;
}

std::string CanonicalPreferredAttemptMapText(
    const std::vector<int>& Ns,
    const std::vector<int>& BBs,
    const PreferredAttemptMap& preferred)
{
    if (preferred.empty()) return "none";
    std::ostringstream canonical;
    bool first_record = true;
    for (int n_value : Ns) for (int bb_value : BBs)
    {
        const PreferredAttemptMap::const_iterator it = preferred.find(
            PreferredAttemptKey(
                (uint32_t)n_value, (uint32_t)bb_value));
        if (it == preferred.end()) continue;
        if (!first_record) canonical << '|';
        first_record = false;
        canonical << n_value << '@' << bb_value << '=';
        if (it->second.empty()) {
            canonical << "none";
            continue;
        }
        for (size_t i = 0u; i < it->second.size(); ++i) {
            if (i != 0u) canonical << ',';
            canonical << it->second[i];
        }
    }
    return canonical.str();
}

bool IsLowerHexSha256(const std::string& text)
{
    if (text.size() != 64u) return false;
    for (char ch : text) {
        if (!((ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f'))) {
            return false;
        }
    }
    return true;
}

uint32_t Sha256RotateRight(uint32_t value, uint32_t shift)
{
    return (value >> shift) | (value << (32u - shift));
}

// Small benchmark-only SHA-256 implementation used to bind a trusted route
// manifest to the exact bytes produced by route mode.  Keeping verification in
// the consumer prevents a stale or forged preferred map from being blessed by
// an unrelated digest supplied on the command line.
std::string Sha256Hex(const std::string& input)
{
    static const uint32_t kRound[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
        0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
        0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
        0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
        0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
        0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
        0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
        0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
        0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
        0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
        0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
        0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u,
    };
    uint32_t state[8] = {
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u,
    };
    std::vector<uint8_t> bytes(input.begin(), input.end());
    const uint64_t bit_count = (uint64_t)bytes.size() * 8u;
    bytes.push_back(0x80u);
    while ((bytes.size() & 63u) != 56u) bytes.push_back(0u);
    for (int shift = 56; shift >= 0; shift -= 8) {
        bytes.push_back((uint8_t)(bit_count >> shift));
    }
    for (size_t offset = 0u; offset < bytes.size(); offset += 64u)
    {
        uint32_t words[64];
        for (uint32_t i = 0u; i < 16u; ++i) {
            const size_t pos = offset + (size_t)i * 4u;
            words[i] = ((uint32_t)bytes[pos] << 24) |
                ((uint32_t)bytes[pos + 1u] << 16) |
                ((uint32_t)bytes[pos + 2u] << 8) |
                (uint32_t)bytes[pos + 3u];
        }
        for (uint32_t i = 16u; i < 64u; ++i) {
            const uint32_t x = words[i - 15u];
            const uint32_t y = words[i - 2u];
            const uint32_t s0 = Sha256RotateRight(x, 7u) ^
                Sha256RotateRight(x, 18u) ^ (x >> 3u);
            const uint32_t s1 = Sha256RotateRight(y, 17u) ^
                Sha256RotateRight(y, 19u) ^ (y >> 10u);
            words[i] = words[i - 16u] + s0 + words[i - 7u] + s1;
        }
        uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
        uint32_t e = state[4], f = state[5], g = state[6], h = state[7];
        for (uint32_t i = 0u; i < 64u; ++i) {
            const uint32_t s1 = Sha256RotateRight(e, 6u) ^
                Sha256RotateRight(e, 11u) ^ Sha256RotateRight(e, 25u);
            const uint32_t choose = (e & f) ^ ((~e) & g);
            const uint32_t temp1 = h + s1 + choose + kRound[i] + words[i];
            const uint32_t s0 = Sha256RotateRight(a, 2u) ^
                Sha256RotateRight(a, 13u) ^ Sha256RotateRight(a, 22u);
            const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t temp2 = s0 + majority;
            h = g; g = f; f = e; e = d + temp1;
            d = c; c = b; b = a; a = temp1 + temp2;
        }
        state[0] += a; state[1] += b; state[2] += c; state[3] += d;
        state[4] += e; state[5] += f; state[6] += g; state[7] += h;
    }
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (uint32_t value : state) output << std::setw(8) << value;
    return output.str();
}

bool ReadBoundedFile(
    const std::string& path,
    size_t max_bytes,
    std::string& bytes)
{
    bytes.clear();
    std::ifstream input(path.c_str(), std::ios::binary | std::ios::ate);
    if (!input) return false;
    const std::ifstream::pos_type end = input.tellg();
    if (end < 0 || (uint64_t)end > max_bytes) return false;
    bytes.resize((size_t)end);
    input.seekg(0, std::ios::beg);
    if (!bytes.empty() &&
        !input.read(&bytes[0], (std::streamsize)bytes.size()))
    {
        bytes.clear();
        return false;
    }
    return true;
}

bool ParseCanonicalU32(const std::string& text, uint32_t& value)
{
    if (!ParseU32Scalar(text.c_str(), value)) return false;
    char canonical[16];
    std::snprintf(canonical, sizeof(canonical), "%u", value);
    return text == canonical;
}

bool ParseCanonicalU64Decimal(const std::string& text, uint64_t& value)
{
    if (text.empty() || text[0] < '0' || text[0] > '9') return false;
    errno = 0;
    char* end = nullptr;
    const unsigned long long parsed =
        std::strtoull(text.c_str(), &end, 10);
    if (errno != 0 || !end || *end != '\0' ||
        parsed > std::numeric_limits<uint64_t>::max())
    {
        return false;
    }
    value = (uint64_t)parsed;
    return text == std::to_string(value);
}

bool ParseCanonicalU32Arg(
    const char* option,
    const char* text,
    uint32_t& value)
{
    return ParseCanonicalU32(text ? text : "", value) || BadArg(option, text);
}

bool ParseCanonicalU64Arg(
    const char* option,
    const char* text,
    uint64_t& value)
{
    return ParseCanonicalU64Decimal(text ? text : "", value) ||
        BadArg(option, text);
}

struct PreferredRouteRecord
{
    uint32_t PreferredAttempt = 0u;
    uint32_t CanonicalAttempt = 0u;
    uint32_t ActualAttempt = 0u;
    bool PreferredValid = false;
    bool Fallback = false;
    bool NoOp = false;
    bool Direct = false;
    uint32_t CanonicalProbeSolves = 0u;
    uint32_t PreferredProbeSolves = 0u;
};

struct PreferredRouteKeyRecords
{
    uint32_t CanonicalAttempt = UINT32_MAX;
    bool ControlOnly = false;
    uint32_t ControlCanonicalProbeSolves = 0u;
    // Preserve manifest order independently of the lookup map.  Candidate
    // jobs must consume the exact ordered attempt batch that produced their
    // bound route artifact; accepting a strict subset silently changes the
    // frozen logical/physical-work accounting.
    std::vector<uint32_t> AttemptOrder;
    std::map<uint32_t, PreferredRouteRecord> Attempts;
};

typedef std::map<PreferredAttemptKey, PreferredRouteKeyRecords>
    PreferredRouteCache;

bool SplitExactCsv(
    const std::string& line,
    size_t field_count,
    std::vector<std::string>& fields)
{
    fields.clear();
    size_t begin = 0u;
    while (begin <= line.size()) {
        const size_t comma = line.find(',', begin);
        fields.push_back(line.substr(
            begin, comma == std::string::npos ? std::string::npos :
                comma - begin));
        if (comma == std::string::npos) break;
        begin = comma + 1u;
    }
    if (fields.size() != field_count) return false;
    for (const std::string& field : fields) {
        if (field.empty()) return false;
    }
    return true;
}

bool ParsePreferredRouteManifest(
    const std::string& bytes,
    const std::string& expected_context_sha256,
    PreferredRouteCache& cache,
    std::string& error)
{
    static const char kMagicPrefix[] =
        "# preferredattempt-route-manifest: schema=v1 "
        "policy=h12-q0-adaptive canonical=ascending-first-valid-v1 "
        "max_attempt=255 context_sha256=";
    static const char kHeader[] =
        "N,bb,route_status,preferred_attempt,canonical_attempt,actual_attempt,"
        "preferred_valid,fallback,no_op,direct,canonical_probe_solves,"
        "preferred_probe_solves";
    cache.clear();
    if (bytes.empty() || bytes.back() != '\n' ||
        bytes.find('\r') != std::string::npos ||
        bytes.find('\0') != std::string::npos)
    {
        error = "route manifest must be nonempty canonical LF text";
        return false;
    }
    std::istringstream input(bytes);
    std::string line;
    if (!std::getline(input, line) ||
        line.size() != std::strlen(kMagicPrefix) + 64u ||
        line.compare(0u, std::strlen(kMagicPrefix), kMagicPrefix) != 0 ||
        line.substr(std::strlen(kMagicPrefix)) != expected_context_sha256 ||
        !std::getline(input, line) || line != kHeader)
    {
        error = "route manifest magic/schema mismatch";
        return false;
    }
    size_t row_count = 0u;
    bool have_previous_key = false;
    PreferredAttemptKey previous_key;
    std::vector<std::string> fields;
    while (std::getline(input, line))
    {
        if (!SplitExactCsv(line, 12u, fields)) {
            error = "route manifest has a malformed row";
            return false;
        }
        uint32_t K = 0u;
        uint32_t block_bytes = 0u;
        uint32_t values[8];
        if (!ParseCanonicalU32(fields[0], K) ||
            !ParseCanonicalU32(fields[1], block_bytes))
        {
            error = "route manifest has a noncanonical integer";
            return false;
        }
        for (size_t i = 0u; i < 8u; ++i) {
            if (!ParseCanonicalU32(fields[i + 4u], values[i])) {
                error = "route manifest has a noncanonical integer";
                return false;
            }
        }
        const bool control_only = fields[2] == "control";
        const bool preferred_route = fields[2] == "preferred";
        uint32_t p = 0u;
        if ((!control_only && !preferred_route) ||
            (control_only && fields[3] != "-1") ||
            (preferred_route && !ParseCanonicalU32(fields[3], p)))
        {
            error = "route manifest has an invalid route status";
            return false;
        }
        const uint32_t a0 = values[0];
        const uint32_t actual = values[1];
        if (K < 2u || K > 64000u ||
            (block_bytes != 64u && block_bytes != 256u &&
             block_bytes != 1280u && block_bytes != 4096u) ||
            (preferred_route &&
             p >= wirehair_v2::kMaxPacketSeedAttempts) ||
            a0 >= wirehair_v2::kMaxPacketSeedAttempts ||
            actual >= wirehair_v2::kMaxPacketSeedAttempts ||
            values[2] > 1u || values[3] > 1u || values[4] > 1u ||
            values[5] > 1u)
        {
            error = "route manifest has an out-of-domain value";
            return false;
        }
        const bool valid = values[2] != 0u;
        const bool fallback = values[3] != 0u;
        const bool no_op = values[4] != 0u;
        const bool direct = values[5] != 0u;
        const uint32_t canonical_probe_solves = values[6];
        const uint32_t preferred_probe_solves = values[7];
        const bool control_consistent = control_only &&
            actual == a0 && valid && !fallback && no_op && !direct &&
            canonical_probe_solves == a0 + 1u &&
            preferred_probe_solves == 0u;
        const bool preferred_consistent = preferred_route &&
            !(p < a0 && valid) && !(p == a0 && !valid) &&
            actual == (valid ? p : a0) && fallback != valid &&
            no_op == (p == a0) && direct == (valid && !no_op) &&
            (canonical_probe_solves == 0u ||
             canonical_probe_solves == a0 + 1u) &&
            preferred_probe_solves == (p > a0 ? 1u : 0u);
        if (!control_consistent && !preferred_consistent)
        {
            error = "route manifest has an inconsistent classification";
            return false;
        }
        const PreferredAttemptKey key(K, block_bytes);
        if (have_previous_key && key != previous_key && key < previous_key) {
            error = "route manifest key groups are reordered";
            return false;
        }
        previous_key = key;
        have_previous_key = true;
        PreferredRouteKeyRecords& key_records = cache[key];
        if (key_records.CanonicalAttempt != UINT32_MAX &&
            key_records.CanonicalAttempt != a0)
        {
            error = "route manifest disagrees on canonical attempt";
            return false;
        }
        key_records.CanonicalAttempt = a0;
        if (control_only) {
            if (key_records.ControlOnly || !key_records.Attempts.empty()) {
                error = "route manifest mixes control and preferred routes";
                return false;
            }
            key_records.ControlOnly = true;
            key_records.ControlCanonicalProbeSolves =
                canonical_probe_solves;
            ++row_count;
            continue;
        }
        if (key_records.ControlOnly) {
            error = "route manifest mixes control and preferred routes";
            return false;
        }
        const uint32_t expected_canonical_probe_solves =
            key_records.AttemptOrder.empty() ? a0 + 1u : 0u;
        if (canonical_probe_solves != expected_canonical_probe_solves) {
            error = "route manifest first-row canonical probe accounting is "
                "noncanonical";
            return false;
        }
        PreferredRouteRecord record;
        record.PreferredAttempt = p;
        record.CanonicalAttempt = a0;
        record.ActualAttempt = actual;
        record.PreferredValid = valid;
        record.Fallback = fallback;
        record.NoOp = no_op;
        record.Direct = direct;
        record.CanonicalProbeSolves = canonical_probe_solves;
        record.PreferredProbeSolves = preferred_probe_solves;
        if (!key_records.Attempts.insert(std::make_pair(p, record)).second) {
            error = "route manifest has a duplicate K@bb,p row";
            return false;
        }
        key_records.AttemptOrder.push_back(p);
        if (key_records.Attempts.size() > 32u) {
            error = "route manifest has more than 32 attempts for one K@bb";
            return false;
        }
        ++row_count;
    }
    if (row_count == 0u) {
        error = "route manifest has no rows";
        return false;
    }
    for (const std::pair<const PreferredAttemptKey,
             PreferredRouteKeyRecords>& key : cache)
    {
        uint64_t canonical_probe_total =
            key.second.ControlCanonicalProbeSolves;
        for (const std::pair<const uint32_t, PreferredRouteRecord>& attempt :
             key.second.Attempts)
        {
            canonical_probe_total += attempt.second.CanonicalProbeSolves;
        }
        if (canonical_probe_total != key.second.CanonicalAttempt + 1u) {
            error = "route manifest canonical probe accounting is not additive";
            return false;
        }
    }
    return true;
}

bool LoadPreferredRouteCache(
    const std::string& path,
    const std::string& expected_sha256,
    const std::string& expected_context_sha256,
    PreferredRouteCache& cache)
{
    std::string bytes;
    if (!ReadBoundedFile(path, 64u * 1024u * 1024u, bytes)) {
        std::fprintf(stderr,
            "preferredattempt cannot read bounded route cache: %s\n",
            path.c_str());
        return false;
    }
    const std::string actual_sha256 = Sha256Hex(bytes);
    if (actual_sha256 != expected_sha256) {
        std::fprintf(stderr,
            "preferredattempt route cache SHA-256 mismatch: got %s\n",
            actual_sha256.c_str());
        return false;
    }
    std::string error;
    if (!ParsePreferredRouteManifest(
            bytes, expected_context_sha256, cache, error))
    {
        std::fprintf(stderr, "preferredattempt %s\n", error.c_str());
        return false;
    }
    return true;
}

struct PreferredRecoveryMetrics
{
    WirehairResult Result = Wirehair_Error;
    uint32_t Inactivated = 0u;
    uint32_t BinaryDeficit = 0u;
    uint32_t HeavyGain = 0u;
    uint64_t BlockXors = 0u;
    uint64_t BlockMulAdds = 0u;
    bool HeavyShortfall = false;
};

bool MakeH12Q0Configuration(
    uint32_t K,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig& config)
{
    const wirehair_v2::SeedProfile profile =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    params = wirehair_v2::MakeMixedParams(
        K,
        wirehair_v2::MatrixSeedFromProfile(
            profile, 0u, wirehair_v2::kMessagePrecodeSeedSalt));
    params.Staircase = profile.DenseCount;
    params.DenseTwoAnchor =
        K >= wirehair_v2::kDenseTwoAnchorMinBlockCount;
    params.DenseTwoAnchorPhase = 0u;
    config.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
        profile, wirehair_v2::kMessageRecoveryRowSeedSalt);
    config.MixCount = 2u;
    return true;
}

WirehairResult ProbePreferredAttempt(
    const wirehair_v2::PrecodeParams& base_params,
    const wirehair_v2::PacketRowConfig& base_config,
    uint32_t attempt,
    const std::vector<wirehair_v2::SolvePacket>& systematic_packets,
    wirehair_v2::PrecodeSystem& system,
    wirehair_v2::PacketRowConfig& config)
{
    if (attempt >= wirehair_v2::kMaxPacketSeedAttempts ||
        systematic_packets.size() != base_params.BlockCount)
    {
        return Wirehair_InvalidInput;
    }
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::PrecodeParamsForAttempt(base_params, attempt),
            system))
    {
        return Wirehair_InvalidInput;
    }
    config = wirehair_v2::PacketConfigForAttempt(base_config, attempt);
    std::vector<uint8_t> intermediate;
    return wirehair_v2::SolvePrecodeSystem(
        system, config, systematic_packets, 2u, intermediate);
}

struct PreferredRecoveryCell
{
    wirehair_v2::PacketRowRuntime Runtime;
    std::vector<wirehair_v2::SolvePacket> Packets;
    uint8_t Zero[2] = {0u, 0u};
};

WirehairResult BuildPreferredRecoveryCell(
    const wirehair_v2::PrecodeParams& params,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_bytes,
    double loss,
    uint64_t external_seed,
    PacketScheduleKind schedule,
    PreferredRecoveryCell& cell)
{
    const uint32_t K = params.BlockCount;
    const uint64_t precode_count_wide =
        (uint64_t)params.Staircase + params.DenseRows + params.HeavyRows;
    if (precode_count_wide > UINT32_MAX) {
        return Wirehair_InvalidInput;
    }
    if (!cell.Runtime.Initialize(
            K, (uint32_t)precode_count_wide, config.MixCount))
    {
        return Wirehair_InvalidInput;
    }
    const uint64_t loss_seed = external_seed ^
        ((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)) ^
        ((uint64_t)block_bytes * UINT64_C(0xbf58476d1ce4e5b9));
    const std::vector<uint32_t> ids = BuildPacketSchedule(
        K, K, loss, loss_seed, schedule);
    if (ids.size() != K) {
        return Wirehair_InvalidInput;
    }
    cell.Packets.resize(K);
    for (uint32_t i = 0u; i < K; ++i) {
        cell.Packets[i].BlockId = ids[i];
        cell.Packets[i].Data = cell.Zero;
    }
    return Wirehair_Success;
}

WirehairResult EvaluatePreferredRecovery(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const PreferredRecoveryCell& cell,
    PreferredRecoveryMetrics& metrics)
{
    if (cell.Packets.size() != system.Params.BlockCount) {
        return Wirehair_InvalidInput;
    }
    std::vector<uint8_t> intermediate;
    wirehair_v2::PrecodeSolveStats stats;
    metrics.Result =
        wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
            system, config, cell.Runtime, cell.Packets, 2u,
            intermediate, &stats);
    metrics.Inactivated = stats.InactivatedColumns;
    metrics.BinaryDeficit = stats.InactivatedColumns >=
            stats.BinaryResidualRank ?
        stats.InactivatedColumns - stats.BinaryResidualRank : 0u;
    metrics.HeavyGain = stats.ResidualRank >= stats.BinaryResidualRank ?
        stats.ResidualRank - stats.BinaryResidualRank : 0u;
    metrics.BlockXors = stats.BlockXors;
    metrics.BlockMulAdds = stats.BlockMulAdds;
    metrics.HeavyShortfall = metrics.Result == Wirehair_NeedMore &&
        metrics.BinaryDeficit <= system.Params.HeavyRows &&
        metrics.HeavyGain < metrics.BinaryDeficit;
    return Wirehair_Success;
}

void PrintPreferredRecoveryRow(
    uint32_t K,
    uint32_t block_bytes,
    const char* arm,
    int preferred_attempt,
    uint32_t canonical_attempt,
    uint32_t actual_attempt,
    bool routed,
    bool preferred_valid,
    bool fallback,
    bool no_op,
    bool direct,
    bool physical_solve,
    const PreferredRecoveryMetrics& metrics)
{
    const uint32_t rank_fail = metrics.Result == Wirehair_NeedMore ? 1u : 0u;
    const uint32_t error =
        metrics.Result == Wirehair_Success ||
        metrics.Result == Wirehair_NeedMore ? 0u : 1u;
    std::printf(
        "%u,%u,%s,%d,%u,%u,%u,%u,%u,%u,%u,%u,%d,%u,%u,%u,%u,%u,%u,%llu,%llu\n",
        K, block_bytes, arm, preferred_attempt, canonical_attempt,
        actual_attempt, routed ? 1u : 0u, preferred_valid ? 1u : 0u,
        fallback ? 1u : 0u, no_op ? 1u : 0u, direct ? 1u : 0u,
        physical_solve ? 1u : 0u, (int)metrics.Result, rank_fail, error,
        metrics.HeavyShortfall ? 1u : 0u, metrics.Inactivated,
        metrics.BinaryDeficit, metrics.HeavyGain,
        (unsigned long long)metrics.BlockXors,
        (unsigned long long)metrics.BlockMulAdds);
}

int CmdPreferredAttempt(int argc, char** argv)
{
    std::string nlist;
    std::string bb_list;
    std::string preferred_map_text;
    std::string route_cache_path;
    std::string route_cache_sha256;
    std::string route_context_sha256;
    PreferredAttemptMode mode = PreferredAttemptMode::Control;
    bool nlist_explicit = false;
    bool bb_list_explicit = false;
    bool mode_explicit = false;
    bool map_explicit = false;
    bool route_cache_path_explicit = false;
    bool route_cache_explicit = false;
    bool route_context_explicit = false;
    bool probe_route = false;
    bool probe_route_explicit = false;
    bool loss_explicit = false;
    bool seed_explicit = false;
    bool schedule_explicit = false;
    double loss = 0.50;
    uint64_t seed = UINT64_C(0x5eedf411);
    PacketScheduleKind schedule = PacketScheduleKind::Burst;

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (!MarkPreferredOptionOnce("--N", nlist_explicit) || !TakeArg(
                    "preferredattempt", "--N", argc, argv, i, value))
            {
                return 1;
            }
            nlist = value;
        }
        else if (!std::strcmp(argv[i], "--bb-list")) {
            if (!MarkPreferredOptionOnce("--bb-list", bb_list_explicit) ||
                !TakeArg(
                    "preferredattempt", "--bb-list", argc, argv, i, value))
            {
                return 1;
            }
            bb_list = value;
        }
        else if (!std::strcmp(argv[i], "--mode")) {
            if (!MarkPreferredOptionOnce("--mode", mode_explicit) || !TakeArg(
                    "preferredattempt", "--mode", argc, argv, i, value) ||
                !ParsePreferredAttemptMode(value, mode))
            {
                std::fprintf(stderr,
                    "preferredattempt --mode must be control, route, "
                    "candidate, or paired\n");
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--preferred-map")) {
            if (!MarkPreferredOptionOnce("--preferred-map", map_explicit) ||
                !TakeArg(
                    "preferredattempt", "--preferred-map",
                    argc, argv, i, value))
            {
                return 1;
            }
            preferred_map_text = value;
        }
        else if (!std::strcmp(argv[i], "--route-cache-sha256")) {
            if (!MarkPreferredOptionOnce(
                    "--route-cache-sha256", route_cache_explicit) ||
                !TakeArg(
                    "preferredattempt", "--route-cache-sha256",
                    argc, argv, i, value))
            {
                return 1;
            }
            route_cache_sha256 = value;
        }
        else if (!std::strcmp(argv[i], "--route-cache")) {
            if (!MarkPreferredOptionOnce(
                    "--route-cache", route_cache_path_explicit) ||
                !TakeArg(
                    "preferredattempt", "--route-cache",
                    argc, argv, i, value))
            {
                return 1;
            }
            route_cache_path = value;
        }
        else if (!std::strcmp(argv[i], "--route-context-sha256")) {
            if (!MarkPreferredOptionOnce(
                    "--route-context-sha256", route_context_explicit) ||
                !TakeArg(
                    "preferredattempt", "--route-context-sha256",
                    argc, argv, i, value))
            {
                return 1;
            }
            route_context_sha256 = value;
        }
        else if (!std::strcmp(argv[i], "--probe-route")) {
            if (!MarkPreferredOptionOnce("--probe-route", probe_route_explicit)) {
                return 1;
            }
            probe_route = true;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (!MarkPreferredOptionOnce("--loss", loss_explicit) || !TakeArg(
                    "preferredattempt", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (!MarkPreferredOptionOnce("--seed", seed_explicit) || !TakeArg(
                    "preferredattempt", "--seed", argc, argv, i, value) ||
                !ParseU64Arg("--seed", value, seed))
            {
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--schedule")) {
            if (!MarkPreferredOptionOnce(
                    "--schedule", schedule_explicit) || !TakeArg(
                    "preferredattempt", "--schedule", argc, argv, i, value) ||
                !ParsePacketSchedule(value, schedule))
            {
                std::fprintf(stderr,
                    "preferredattempt --schedule must be burst, "
                    "adversarial, or repair-only\n");
                return 1;
            }
        }
        else if (!UnknownArg("preferredattempt", argv[i])) {
            return 1;
        }
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    PreferredAttemptMap preferred;
    if (!mode_explicit || Ns.empty() || BBs.empty() ||
        !ValidateBlockCounts(Ns, "preferredattempt") ||
        !ValidateMessageInputs(Ns, BBs, "preferredattempt") ||
        !ValidateLoss(loss, "preferredattempt"))
    {
        if (!mode_explicit || Ns.empty() || BBs.empty()) {
            std::fprintf(stderr,
                "preferredattempt requires explicit --mode, --N, and "
                "--bb-list\n");
        }
        return 1;
    }
    if (!IsCanonicalPreferredIntList(nlist, Ns) ||
        !IsCanonicalPreferredIntList(bb_list, BBs) ||
        !std::is_sorted(Ns.begin(), Ns.end()) ||
        !std::is_sorted(BBs.begin(), BBs.end()))
    {
        std::fprintf(stderr,
            "preferredattempt --N and --bb-list must be canonical ascending "
            "decimal lists\n");
        return 1;
    }
    if (schedule != PacketScheduleKind::Burst &&
        schedule != PacketScheduleKind::Adversarial &&
        schedule != PacketScheduleKind::RepairOnly)
    {
        std::fprintf(stderr,
            "preferredattempt --schedule must be burst, adversarial, or "
            "repair-only\n");
        return 1;
    }
    for (size_t i = 0u; i < Ns.size(); ++i) {
        if (std::find(Ns.begin(), Ns.begin() + (ptrdiff_t)i, Ns[i]) !=
            Ns.begin() + (ptrdiff_t)i)
        {
            std::fprintf(stderr, "preferredattempt --N contains duplicates\n");
            return 1;
        }
    }
    for (size_t i = 0u; i < BBs.size(); ++i) {
        if ((BBs[i] & 1) != 0 ||
            std::find(BBs.begin(), BBs.begin() + (ptrdiff_t)i, BBs[i]) !=
                BBs.begin() + (ptrdiff_t)i)
        {
            std::fprintf(stderr,
                "preferredattempt --bb-list must contain unique even values\n");
            return 1;
        }
    }
    if (BBs.size() > 4u) {
        std::fprintf(stderr,
            "preferredattempt --bb-list supports at most four widths\n");
        return 1;
    }
    for (int block_bytes : BBs) {
        if (block_bytes != 64 && block_bytes != 256 &&
            block_bytes != 1280 && block_bytes != 4096)
        {
            std::fprintf(stderr,
                "preferredattempt --bb-list supports only "
                "64,256,1280,4096\n");
            return 1;
        }
    }
    if (map_explicit && !ParsePreferredAttemptMap(
            preferred_map_text, preferred))
    {
        std::fprintf(stderr,
            "preferredattempt malformed --preferred-map (expected "
            "K@bb=p[,p...] or route-only K@bb=none, max 32 unique p "
            "values)\n");
        return 1;
    }
    for (PreferredAttemptMap::const_iterator it = preferred.begin();
         it != preferred.end(); ++it)
    {
        const bool known_K = std::find_if(
            Ns.begin(), Ns.end(), [&](int value) {
                return (uint32_t)value == it->first.first;
            }) != Ns.end();
        const bool known_width = std::find_if(
            BBs.begin(), BBs.end(), [&](int value) {
                return (uint32_t)value == it->first.second;
            }) != BBs.end();
        if (!known_K || !known_width)
        {
            std::fprintf(stderr,
                "preferredattempt map key is outside --N x --bb-list\n");
            return 1;
        }
        if (mode == PreferredAttemptMode::Paired && it->second.size() != 1u) {
            std::fprintf(stderr,
                "preferredattempt paired mode requires one p per mapped "
                "K@bb\n");
            return 1;
        }
        if (mode == PreferredAttemptMode::Candidate && it->second.empty()) {
            std::fprintf(stderr,
                "preferredattempt candidate logical p list cannot be empty\n");
            return 1;
        }
    }
    if (mode == PreferredAttemptMode::Control && map_explicit) {
        std::fprintf(stderr,
            "preferredattempt control mode rejects --preferred-map\n");
        return 1;
    }
    if (mode == PreferredAttemptMode::Control) {
        if (route_context_explicit) {
            std::fprintf(stderr,
                "preferredattempt control mode rejects route context\n");
            return 1;
        }
    }
    else if (!route_context_explicit ||
             !IsLowerHexSha256(route_context_sha256))
    {
        std::fprintf(stderr,
            "preferredattempt route/candidate/paired modes require a "
            "lowercase --route-context-sha256\n");
        return 1;
    }
    if ((mode == PreferredAttemptMode::Route ||
         mode == PreferredAttemptMode::Candidate) &&
        (!map_explicit || preferred.empty()))
    {
        std::fprintf(stderr,
            "preferredattempt route/candidate mode requires a nonempty map\n");
        return 1;
    }
    if ((mode == PreferredAttemptMode::Route ||
         mode == PreferredAttemptMode::Candidate) &&
        preferred.size() != Ns.size() * BBs.size())
    {
        std::fprintf(stderr,
            "preferredattempt route/candidate map must cover the exact "
            "--N x --bb-list cross-product\n");
        return 1;
    }
    if (map_explicit)
    {
        for (int n_value : Ns)
        {
            const std::vector<uint32_t>* logical_attempts = nullptr;
            size_t present_widths = 0u;
            for (int bb_value : BBs)
            {
                const PreferredAttemptMap::const_iterator it = preferred.find(
                    PreferredAttemptKey(
                        (uint32_t)n_value, (uint32_t)bb_value));
                if (it == preferred.end()) continue;
                ++present_widths;
                if (!logical_attempts) {
                    logical_attempts = &it->second;
                }
                else if (*logical_attempts != it->second) {
                    std::fprintf(stderr,
                        "preferredattempt logical p list must be identical "
                        "across widths for N=%d\n", n_value);
                    return 1;
                }
            }
            if (mode == PreferredAttemptMode::Paired && present_widths != 0u &&
                present_widths != BBs.size())
            {
                std::fprintf(stderr,
                    "preferredattempt paired routing must cover all or no "
                    "widths for N=%d\n", n_value);
                return 1;
            }
        }
    }
    if (map_explicit && preferred_map_text !=
            CanonicalPreferredAttemptMapText(Ns, BBs, preferred))
    {
        std::fprintf(stderr,
            "preferredattempt --preferred-map must use canonical K/width "
            "record order and decimal spelling\n");
        return 1;
    }
    if (mode == PreferredAttemptMode::Paired && !map_explicit) {
        std::fprintf(stderr,
            "preferredattempt paired mode requires explicit "
            "--preferred-map (use none for all aliases)\n");
        return 1;
    }
    if (mode == PreferredAttemptMode::Candidate) {
        const bool have_bound_cache =
            route_cache_path_explicit && route_cache_explicit;
        if (route_cache_path_explicit != route_cache_explicit ||
            probe_route == have_bound_cache ||
            (route_cache_explicit &&
             !IsLowerHexSha256(route_cache_sha256)))
        {
            std::fprintf(stderr,
                "preferredattempt candidate mode requires exactly one of "
                "--probe-route or --route-cache PATH plus a lowercase "
                "--route-cache-sha256\n");
            return 1;
        }
    }
    else if (mode == PreferredAttemptMode::Paired) {
        if (!route_cache_path_explicit || !route_cache_explicit ||
            !IsLowerHexSha256(route_cache_sha256) || probe_route)
        {
            std::fprintf(stderr,
                "preferredattempt paired mode requires --route-cache PATH "
                "and a lowercase --route-cache-sha256\n");
            return 1;
        }
    }
    else if (route_cache_path_explicit || route_cache_explicit || probe_route) {
        std::fprintf(stderr,
            "preferredattempt route-cache options are candidate/paired-only\n");
        return 1;
    }
    if (mode == PreferredAttemptMode::Route) {
        if (loss_explicit || seed_explicit || schedule_explicit) {
            std::fprintf(stderr,
                "preferredattempt route mode rejects recovery-only "
                "--loss/--seed/--schedule options\n");
            return 1;
        }
    }
    else if (!loss_explicit || !seed_explicit || !schedule_explicit) {
        std::fprintf(stderr,
            "preferredattempt control/candidate/paired modes require explicit "
            "--loss, --seed, and --schedule\n");
        return 1;
    }

    PreferredRouteCache route_cache;
    if (route_cache_path_explicit &&
        !LoadPreferredRouteCache(
            route_cache_path, route_cache_sha256,
            route_context_sha256, route_cache))
    {
        return 1;
    }
    if (!route_cache.empty()) {
        if (route_cache.size() != Ns.size() * BBs.size()) {
            std::fprintf(stderr,
                "preferredattempt route cache keys must exactly match "
                "--N x --bb-list\n");
            return 1;
        }
        for (int n_value : Ns) for (int bb_value : BBs) {
            const PreferredAttemptKey key(
                (uint32_t)n_value, (uint32_t)bb_value);
            const PreferredRouteCache::const_iterator cache_key =
                route_cache.find(key);
            if (cache_key == route_cache.end()) {
                std::fprintf(stderr,
                    "preferredattempt route cache lacks N=%u bb=%u\n",
                    key.first, key.second);
                return 1;
            }
            const bool mapped = preferred.find(key) != preferred.end();
            if (mode == PreferredAttemptMode::Candidate &&
                cache_key->second.ControlOnly)
            {
                std::fprintf(stderr,
                    "preferredattempt candidate cache contains a control-only "
                    "key N=%u bb=%u\n", key.first, key.second);
                return 1;
            }
            if (mode == PreferredAttemptMode::Candidate)
            {
                const PreferredAttemptMap::const_iterator requested =
                    preferred.find(key);
                if (requested == preferred.end() ||
                    cache_key->second.AttemptOrder != requested->second)
                {
                    std::fprintf(stderr,
                        "preferredattempt candidate route-cache attempt rows "
                        "must exactly match --preferred-map for N=%u bb=%u\n",
                        key.first, key.second);
                    return 1;
                }
            }
            if (mode == PreferredAttemptMode::Paired &&
                cache_key->second.ControlOnly == mapped)
            {
                std::fprintf(stderr,
                    "preferredattempt paired map/cache route-status mismatch "
                    "N=%u bb=%u\n", key.first, key.second);
                return 1;
            }
            if (mode == PreferredAttemptMode::Paired && mapped)
            {
                const PreferredAttemptMap::const_iterator requested =
                    preferred.find(key);
                if (requested == preferred.end() ||
                    cache_key->second.AttemptOrder != requested->second)
                {
                    std::fprintf(stderr,
                        "preferredattempt paired route-cache selected attempt "
                        "row must exactly match --preferred-map for N=%u "
                        "bb=%u\n", key.first, key.second);
                    return 1;
                }
            }
        }
    }

    if (mode == PreferredAttemptMode::Route) {
#if defined(_WIN32)
        if (_setmode(_fileno(stdout), _O_BINARY) == -1) {
            std::fprintf(stderr,
                "preferredattempt cannot force canonical LF output\n");
            return 1;
        }
#endif
        std::printf(
            "# preferredattempt-route-manifest: schema=v1 "
            "policy=h12-q0-adaptive canonical=ascending-first-valid-v1 "
            "max_attempt=255 context_sha256=%s\n"
            "N,bb,route_status,preferred_attempt,canonical_attempt,actual_attempt,"
            "preferred_valid,fallback,no_op,direct,canonical_probe_solves,"
            "preferred_probe_solves\n",
            route_context_sha256.c_str());
    }
    else {
        std::printf(
            "# preferredattempt: schema=v2 policy=h12-q0-adaptive "
            "canonical=ascending-first-valid-v1 mode=%s loss=%.17g "
            "seed=0x%llx schedule=%s preferred_batch_max=32 "
            "route_cache_sha256=%s route_context_sha256=%s probe_route=%u "
            "physical_solve_accounting=explicit "
            "systematic_probe_accounting=explicit\n",
            PreferredAttemptModeName(mode), loss, (unsigned long long)seed,
            PacketScheduleName(schedule),
            route_cache_explicit ? route_cache_sha256.c_str() : "none",
            mode == PreferredAttemptMode::Control ?
                "none" : route_context_sha256.c_str(),
            probe_route ? 1u : 0u);
        std::printf(
            "N,bb,arm,preferred_attempt,canonical_attempt,actual_attempt,"
            "routed,preferred_valid,fallback,no_op,direct,physical_solve,"
            "result,rank_fail,error,heavy_shortfall,inactivated,binary_def,heavy_gain,"
            "block_xors,block_muladds\n");
    }

    for (int n_value : Ns) for (int bb_value : BBs)
    {
        const uint32_t K = (uint32_t)n_value;
        const uint32_t block_bytes = (uint32_t)bb_value;
        const PreferredAttemptKey key(K, block_bytes);
        const PreferredAttemptMap::const_iterator map_it = preferred.find(key);
        const PreferredRouteCache::const_iterator cache_it =
            route_cache.find(key);
        wirehair_v2::PrecodeParams base_params;
        wirehair_v2::PacketRowConfig base_config;
        MakeH12Q0Configuration(K, block_bytes, base_params, base_config);

        wirehair_v2::PrecodeSystem canonical_system;
        wirehair_v2::PacketRowConfig canonical_config;
        uint32_t canonical_attempt = UINT32_MAX;
        uint32_t canonical_probe_solves = 0u;
        if (cache_it != route_cache.end())
        {
            canonical_attempt = cache_it->second.CanonicalAttempt;
            if (canonical_attempt >= wirehair_v2::kMaxPacketSeedAttempts)
            {
                std::fprintf(stderr,
                    "preferredattempt cached canonical attempt invalid "
                    "N=%u bb=%u a0=%u\n",
                    K, block_bytes, canonical_attempt);
                return 2;
            }
            // Candidate-only cached jobs never consume canonical_system: alias
            // observations are joined to the separately cached control cube.
            if (mode == PreferredAttemptMode::Paired)
            {
                if (!wirehair_v2::BuildPrecodeSystem(
                        wirehair_v2::PrecodeParamsForAttempt(
                            base_params, canonical_attempt), canonical_system))
                {
                    std::fprintf(stderr,
                        "preferredattempt cached canonical build failed "
                        "N=%u bb=%u a0=%u\n",
                        K, block_bytes, canonical_attempt);
                    return 2;
                }
                canonical_config = wirehair_v2::PacketConfigForAttempt(
                    base_config, canonical_attempt);
            }
        }
        else
        {
            const WirehairResult canonical_result =
                wirehair_v2::SelectSystematicConfiguration(
                    base_params, base_config, canonical_system,
                    canonical_config, &canonical_attempt);
            if (canonical_result != Wirehair_Success ||
                canonical_attempt >= wirehair_v2::kMaxPacketSeedAttempts)
            {
                std::fprintf(stderr,
                    "preferredattempt canonical selection failed N=%u bb=%u "
                    "result=%d attempt=%u\n",
                    K, block_bytes, (int)canonical_result,
                    canonical_attempt);
                return 2;
            }
            canonical_probe_solves = canonical_attempt + 1u;
        }
        if (mode == PreferredAttemptMode::Route ||
            (mode == PreferredAttemptMode::Candidate && probe_route))
        {
            // Route-only work consumes the selected attempt number, not the
            // potentially large canonical system.  Drop it before preferred
            // probing so only one full candidate graph is live at a time.
            canonical_system = wirehair_v2::PrecodeSystem{};
        }

        const uint8_t systematic_zero[2] = {0u, 0u};
        std::vector<wirehair_v2::SolvePacket> systematic_packets;
        if (mode == PreferredAttemptMode::Route || probe_route) {
            systematic_packets.resize(K);
            for (uint32_t block_id = 0u; block_id < K; ++block_id) {
                systematic_packets[block_id].BlockId = block_id;
                systematic_packets[block_id].Data = systematic_zero;
            }
        }

        if (mode == PreferredAttemptMode::Route)
        {
            // Exact map coverage was checked before any output was emitted.
            if (map_it->second.empty()) {
                std::printf(
                    "%u,%u,control,-1,%u,%u,1,0,1,0,%u,0\n",
                    K, block_bytes, canonical_attempt, canonical_attempt,
                    canonical_probe_solves);
                continue;
            }
            bool charge_canonical_probe = true;
            for (uint32_t attempt : map_it->second)
            {
                bool valid = attempt == canonical_attempt;
                if (attempt > canonical_attempt)
                {
                    wirehair_v2::PrecodeSystem candidate_system;
                    wirehair_v2::PacketRowConfig candidate_config;
                    const WirehairResult probe_result = ProbePreferredAttempt(
                        base_params, base_config, attempt,
                        systematic_packets, candidate_system,
                        candidate_config);
                    if (probe_result == Wirehair_Success) {
                        valid = true;
                    }
                    else if (probe_result != Wirehair_NeedMore) {
                        std::fprintf(stderr,
                            "preferredattempt route probe failed N=%u bb=%u "
                            "p=%u result=%d\n",
                            K, block_bytes, attempt, (int)probe_result);
                        return 2;
                    }
                }
                const bool no_op = attempt == canonical_attempt;
                const bool fallback = !valid;
                const bool direct = valid && !no_op;
                std::printf(
                    "%u,%u,preferred,%u,%u,%u,%u,%u,%u,%u,%u,%u\n",
                    K, block_bytes, attempt, canonical_attempt,
                    valid ? attempt : canonical_attempt,
                    valid ? 1u : 0u, fallback ? 1u : 0u,
                    no_op ? 1u : 0u, direct ? 1u : 0u,
                    charge_canonical_probe ? canonical_probe_solves : 0u,
                    attempt > canonical_attempt ? 1u : 0u);
                charge_canonical_probe = false;
            }
            continue;
        }

        if (mode == PreferredAttemptMode::Candidate && !probe_route)
        {
            bool any_direct = false;
            for (uint32_t attempt : map_it->second) {
                const std::map<uint32_t, PreferredRouteRecord>::const_iterator
                    record = cache_it->second.Attempts.find(attempt);
                if (record == cache_it->second.Attempts.end()) {
                    std::fprintf(stderr,
                        "preferredattempt route cache lacks N=%u bb=%u p=%u\n",
                        K, block_bytes, attempt);
                    return 1;
                }
                any_direct = any_direct || record->second.Direct;
            }
            if (!any_direct) {
                for (uint32_t attempt : map_it->second) {
                    const PreferredRouteRecord& record =
                        cache_it->second.Attempts.find(attempt)->second;
                    std::printf(
                        "# preferred_candidate_alias: N=%u bb=%u p=%u "
                        "a0=%u actual=%u valid=%u fallback=%u no_op=%u "
                        "direct=0 physical_solve=0\n",
                        K, block_bytes, attempt, canonical_attempt,
                        record.ActualAttempt,
                        record.PreferredValid ? 1u : 0u,
                        record.Fallback ? 1u : 0u,
                        record.NoOp ? 1u : 0u);
                }
                continue;
            }
        }

        PreferredRecoveryCell recovery_cell;
        const WirehairResult cell_result = BuildPreferredRecoveryCell(
            base_params, base_config, block_bytes, loss, seed, schedule,
            recovery_cell);
        if (cell_result != Wirehair_Success) {
            std::fprintf(stderr,
                "preferredattempt recovery cell setup failed N=%u bb=%u "
                "result=%d\n", K, block_bytes, (int)cell_result);
            return 2;
        }

        PreferredRecoveryMetrics control_metrics;
        if (mode == PreferredAttemptMode::Control ||
            mode == PreferredAttemptMode::Paired)
        {
            const WirehairResult evaluate_result = EvaluatePreferredRecovery(
                canonical_system, canonical_config, recovery_cell,
                control_metrics);
            if (evaluate_result != Wirehair_Success) {
                std::fprintf(stderr,
                    "preferredattempt control evaluation failed N=%u bb=%u "
                    "result=%d\n",
                    K, block_bytes, (int)evaluate_result);
                return 2;
            }
            if (canonical_probe_solves != 0u) {
                std::printf(
                    "# preferred_probe_accounting: N=%u bb=%u "
                    "canonical_probe_solves=%u preferred_probe_solves=0\n",
                    K, block_bytes, canonical_probe_solves);
            }
            PrintPreferredRecoveryRow(
                K, block_bytes, "control", -1, canonical_attempt,
                canonical_attempt, false, true, false, false, false, true,
                control_metrics);
        }
        if (mode == PreferredAttemptMode::Paired) {
            // The paired control observation is complete.  Candidate systems
            // are independent, so avoid retaining both full graphs at once.
            canonical_system = wirehair_v2::PrecodeSystem{};
        }
        if (mode == PreferredAttemptMode::Control) continue;

        if (mode == PreferredAttemptMode::Paired && map_it == preferred.end())
        {
            PrintPreferredRecoveryRow(
                K, block_bytes, "candidate", -1, canonical_attempt,
                canonical_attempt, false, true, false, true, false, false,
                control_metrics);
            continue;
        }

        if (probe_route) {
            std::printf(
                "# preferred_probe_accounting: N=%u bb=%u "
                "canonical_probe_solves=%u\n",
                K, block_bytes, canonical_probe_solves);
        }
        for (uint32_t attempt : map_it->second)
        {
            wirehair_v2::PrecodeSystem candidate_system;
            wirehair_v2::PacketRowConfig candidate_config =
                wirehair_v2::PacketConfigForAttempt(base_config, attempt);
            bool valid = false;
            bool fallback = false;
            bool no_op = false;
            bool direct = false;
            uint32_t actual_attempt = canonical_attempt;
            uint32_t preferred_probe_solves = 0u;

            if (probe_route)
            {
                valid = attempt == canonical_attempt;
                if (attempt > canonical_attempt)
                {
                    preferred_probe_solves = 1u;
                    const WirehairResult probe_result = ProbePreferredAttempt(
                        base_params, base_config, attempt,
                        systematic_packets, candidate_system,
                        candidate_config);
                    if (probe_result == Wirehair_Success) {
                        valid = true;
                    }
                    else if (probe_result != Wirehair_NeedMore) {
                        std::fprintf(stderr,
                            "preferredattempt route probe failed N=%u bb=%u "
                            "p=%u result=%d\n",
                            K, block_bytes, attempt, (int)probe_result);
                        return 2;
                    }
                }
                // Attempts below a0 are known invalid from the ascending
                // selector; p==a0 is a neutral alias and is never copied.
                no_op = attempt == canonical_attempt;
                fallback = !valid;
                direct = valid && !no_op;
                actual_attempt = direct ? attempt : canonical_attempt;
                std::printf(
                    "# preferred_route: N=%u bb=%u p=%u a0=%u actual=%u "
                    "valid=%u fallback=%u no_op=%u direct=%u "
                    "preferred_probe_solves=%u\n",
                    K, block_bytes, attempt, canonical_attempt,
                    actual_attempt, valid ? 1u : 0u,
                    fallback ? 1u : 0u, no_op ? 1u : 0u,
                    direct ? 1u : 0u, preferred_probe_solves);
                if (!direct) continue;
            }
            else
            {
                const PreferredRouteKeyRecords& key_records = cache_it->second;
                const std::map<uint32_t, PreferredRouteRecord>::const_iterator
                    record_it = key_records.Attempts.find(attempt);
                if (record_it == key_records.Attempts.end()) {
                    std::fprintf(stderr,
                        "preferredattempt route cache lacks N=%u bb=%u p=%u\n",
                        K, block_bytes, attempt);
                    return 1;
                }
                const PreferredRouteRecord& record = record_it->second;
                valid = record.PreferredValid;
                fallback = record.Fallback;
                no_op = record.NoOp;
                direct = record.Direct;
                actual_attempt = record.ActualAttempt;
                if (!direct)
                {
                    if (mode == PreferredAttemptMode::Candidate) {
                        std::printf(
                            "# preferred_candidate_alias: N=%u bb=%u p=%u "
                            "a0=%u actual=%u valid=%u fallback=%u no_op=%u "
                            "direct=0 physical_solve=0\n",
                            K, block_bytes, attempt, canonical_attempt,
                            actual_attempt, valid ? 1u : 0u,
                            fallback ? 1u : 0u, no_op ? 1u : 0u);
                        continue;
                    }
                    PrintPreferredRecoveryRow(
                        K, block_bytes, "candidate", (int)attempt,
                        canonical_attempt, canonical_attempt, true, valid,
                        fallback, no_op, false, false, control_metrics);
                    continue;
                }
                if (attempt <= canonical_attempt) {
                    std::fprintf(stderr,
                        "preferredattempt direct cache record is not above "
                        "a0 N=%u bb=%u p=%u a0=%u\n",
                        K, block_bytes, attempt, canonical_attempt);
                    return 1;
                }
                if (!wirehair_v2::BuildPrecodeSystem(
                        wirehair_v2::PrecodeParamsForAttempt(
                            base_params, attempt), candidate_system))
                {
                    std::fprintf(stderr,
                        "preferredattempt cached candidate build failed "
                        "N=%u bb=%u p=%u\n", K, block_bytes, attempt);
                    return 2;
                }
            }

            PreferredRecoveryMetrics candidate_metrics;
            const WirehairResult evaluate_result = EvaluatePreferredRecovery(
                candidate_system, candidate_config, recovery_cell,
                candidate_metrics);
            if (evaluate_result != Wirehair_Success) {
                std::fprintf(stderr,
                    "preferredattempt candidate evaluation failed "
                    "N=%u bb=%u p=%u result=%d\n",
                    K, block_bytes, attempt, (int)evaluate_result);
                return 2;
            }
            PrintPreferredRecoveryRow(
                K, block_bytes, "candidate", (int)attempt,
                canonical_attempt, actual_attempt, true, valid,
                fallback, no_op, direct, true, candidate_metrics);
        }
    }
    return 0;
}

enum class PreferredTimingMetric
{
    Solve,
    Setup
};

const char* PreferredTimingMetricName(PreferredTimingMetric metric)
{
    return metric == PreferredTimingMetric::Solve ? "solve" : "setup";
}

struct PreferredTimingCell
{
    wirehair_v2::PacketRowRuntime Runtime;
    std::vector<wirehair_v2::SolvePacket> Packets;
    std::vector<uint8_t> PayloadStorage;
    uint8_t* Payload = nullptr;
    std::string TraceSha256;
};

bool BuildFullPayloadTimingCell(
    const wirehair_v2::PrecodeParams& params,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_bytes,
    uint32_t overhead,
    double loss,
    uint64_t external_seed,
    PacketScheduleKind schedule,
    bool salt_overhead,
    const char* trace_tag,
    PreferredTimingCell& cell)
{
    const uint32_t K = params.BlockCount;
    const uint64_t delivered_wide = (uint64_t)K + overhead;
    const uint64_t precode_count_wide =
        (uint64_t)params.Staircase + params.DenseRows + params.HeavyRows;
    const uint64_t payload_bytes_wide = delivered_wide * block_bytes;
    if (delivered_wide > UINT32_MAX || precode_count_wide > UINT32_MAX ||
        payload_bytes_wide > (uint64_t)std::numeric_limits<size_t>::max() - 63u)
    {
        return false;
    }
    if (!cell.Runtime.Initialize(
            K, (uint32_t)precode_count_wide, config.MixCount))
    {
        return false;
    }
    const uint64_t loss_seed = external_seed ^
        ((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)) ^
        ((uint64_t)block_bytes * UINT64_C(0xbf58476d1ce4e5b9)) ^
        ((uint64_t)(salt_overhead ? overhead : 0u) *
            UINT64_C(0x94d049bb133111eb));
    const std::vector<uint32_t> ids = BuildPacketSchedule(
        K, (uint32_t)delivered_wide, loss, loss_seed, schedule);
    if (ids.size() != (size_t)delivered_wide) return false;

    const size_t payload_bytes = (size_t)payload_bytes_wide;
    cell.PayloadStorage.resize(payload_bytes + 63u);
    const uintptr_t unaligned =
        reinterpret_cast<uintptr_t>(cell.PayloadStorage.data());
    if (unaligned > std::numeric_limits<uintptr_t>::max() - 63u) {
        return false;
    }
    const uintptr_t aligned = (unaligned + 63u) & ~(uintptr_t)63u;
    cell.Payload = reinterpret_cast<uint8_t*>(aligned);
    if ((aligned & 63u) != 0u || aligned < unaligned ||
        aligned - unaligned > 63u)
    {
        return false;
    }
    // Distinct aligned RHS blocks are all zero so every delivered fountain
    // equation is consistent without timing a separate encoder.  resize()
    // zeroes and first-touches the complete allocation before either arm.
    std::memset(cell.Payload, 0, payload_bytes);

    cell.Packets.resize((size_t)delivered_wide);
    std::ostringstream trace;
    trace << trace_tag << "\n"
          << "K=" << K << "\nblock_bytes=" << block_bytes
          << "\nseed=" << external_seed << "\nloss="
          << std::setprecision(17) << loss << "\nschedule="
          << PacketScheduleName(schedule) << "\nids=";
    for (size_t i = 0u; i < ids.size(); ++i) {
        if (i != 0u) trace << ',';
        trace << ids[i];
        cell.Packets[i].BlockId = ids[i];
        cell.Packets[i].Data = cell.Payload + i * (size_t)block_bytes;
    }
    trace << '\n';
    cell.TraceSha256 = Sha256Hex(trace.str());
    return true;
}

bool BuildPreferredTimingCell(
    const wirehair_v2::PrecodeParams& params,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_bytes,
    double loss,
    uint64_t external_seed,
    PacketScheduleKind schedule,
    PreferredTimingCell& cell)
{
    return BuildFullPayloadTimingCell(
        params, config, block_bytes, 4u, loss, external_seed, schedule,
        false, "wirehair-wh2-preferred-timing-trace-v1", cell);
}

static volatile uint8_t PreferredTimingEvictionSink = 0u;

void EvictPreferredTimingCache(std::vector<uint8_t>& eviction)
{
    uint8_t accumulator = PreferredTimingEvictionSink;
    for (size_t offset = 0u; offset < eviction.size(); offset += 64u) {
        eviction[offset] = (uint8_t)(eviction[offset] + 1u);
        accumulator ^= eviction[offset];
    }
    if (!eviction.empty()) {
        eviction.back() = (uint8_t)(eviction.back() + 1u);
        accumulator ^= eviction.back();
    }
    PreferredTimingEvictionSink = accumulator;
}

int PreferredTimingCurrentCpu()
{
#if defined(__linux__)
    return sched_getcpu();
#else
    return -1;
#endif
}

struct PreferredTimingUsage
{
    int64_t MinorFaults = -1;
    int64_t MajorFaults = -1;
};

PreferredTimingUsage ReadPreferredTimingUsage()
{
    PreferredTimingUsage result;
#if defined(__unix__) || defined(__APPLE__)
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        result.MinorFaults = (int64_t)usage.ru_minflt;
        result.MajorFaults = (int64_t)usage.ru_majflt;
    }
#endif
    return result;
}

uint64_t PreferredTimingElapsedNanoseconds(
    Clock::time_point begin,
    Clock::time_point end,
    bool& saturated)
{
    const std::chrono::nanoseconds elapsed =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    if (elapsed.count() < 0) {
        saturated = true;
        return UINT64_MAX;
    }
    saturated = false;
    return (uint64_t)elapsed.count();
}

int CmdPreferredTiming(int argc, char** argv)
{
    uint32_t K = 0u;
    uint32_t block_bytes = 0u;
    uint32_t preferred_attempt = UINT32_MAX;
    uint32_t cycle_index = 0u;
    uint64_t eviction_bytes = 0u;
    double loss = 0.50;
    uint64_t seed = UINT64_C(0x5eedf411);
    PacketScheduleKind schedule = PacketScheduleKind::Burst;
    PreferredTimingMetric metric = PreferredTimingMetric::Solve;
    std::string route_cache_path;
    std::string route_cache_sha256;
    std::string route_context_sha256;
    bool have_K = false;
    bool have_block_bytes = false;
    bool have_preferred = false;
    bool have_eviction = false;
    bool have_metric = false;
    bool have_cache_path = false;
    bool have_cache_sha256 = false;
    bool have_context = false;
    bool have_loss = false;
    bool have_seed = false;
    bool have_schedule = false;
    bool have_cycle_index = false;

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (have_K ||
                !TakeArg("preferredtiming", "--N", argc, argv, i, value) ||
                !ParseCanonicalU32Arg("--N", value, K)) return 1;
            have_K = true;
        }
        else if (!std::strcmp(argv[i], "--bb")) {
            if (have_block_bytes ||
                !TakeArg("preferredtiming", "--bb", argc, argv, i, value) ||
                !ParseCanonicalU32Arg("--bb", value, block_bytes)) return 1;
            have_block_bytes = true;
        }
        else if (!std::strcmp(argv[i], "--preferred-attempt")) {
            if (have_preferred || !TakeArg(
                    "preferredtiming", "--preferred-attempt",
                    argc, argv, i, value) ||
                !ParseCanonicalU32Arg(
                    "--preferred-attempt", value, preferred_attempt))
            {
                return 1;
            }
            have_preferred = true;
        }
        else if (!std::strcmp(argv[i], "--evict-bytes")) {
            if (have_eviction || !TakeArg(
                    "preferredtiming", "--evict-bytes", argc, argv, i,
                    value) || !ParseCanonicalU64Arg(
                        "--evict-bytes", value, eviction_bytes))
            {
                return 1;
            }
            have_eviction = true;
        }
        else if (!std::strcmp(argv[i], "--metric")) {
            if (have_metric || !TakeArg(
                    "preferredtiming", "--metric", argc, argv, i, value))
            {
                return 1;
            }
            if (!std::strcmp(value, "solve")) {
                metric = PreferredTimingMetric::Solve;
            }
            else if (!std::strcmp(value, "setup")) {
                metric = PreferredTimingMetric::Setup;
            }
            else {
                std::fprintf(stderr,
                    "preferredtiming --metric must be solve or setup\n");
                return 1;
            }
            have_metric = true;
        }
        else if (!std::strcmp(argv[i], "--cycle-index")) {
            if (have_cycle_index || !TakeArg(
                    "preferredtiming", "--cycle-index", argc, argv, i,
                    value) || !ParseCanonicalU32Arg(
                        "--cycle-index", value, cycle_index))
            {
                return 1;
            }
            have_cycle_index = true;
        }
        else if (!std::strcmp(argv[i], "--route-cache")) {
            if (have_cache_path || !TakeArg(
                    "preferredtiming", "--route-cache", argc, argv, i,
                    value)) return 1;
            route_cache_path = value;
            have_cache_path = true;
        }
        else if (!std::strcmp(argv[i], "--route-cache-sha256")) {
            if (have_cache_sha256 || !TakeArg(
                    "preferredtiming", "--route-cache-sha256", argc, argv,
                    i, value)) return 1;
            route_cache_sha256 = value;
            have_cache_sha256 = true;
        }
        else if (!std::strcmp(argv[i], "--route-context-sha256")) {
            if (have_context || !TakeArg(
                    "preferredtiming", "--route-context-sha256", argc, argv,
                    i, value)) return 1;
            route_context_sha256 = value;
            have_context = true;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (have_loss || !TakeArg(
                    "preferredtiming", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss)) return 1;
            have_loss = true;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (have_seed || !TakeArg(
                    "preferredtiming", "--seed", argc, argv, i, value) ||
                !ParseCanonicalU64Arg("--seed", value, seed)) return 1;
            have_seed = true;
        }
        else if (!std::strcmp(argv[i], "--schedule")) {
            if (have_schedule || !TakeArg(
                    "preferredtiming", "--schedule", argc, argv, i, value) ||
                !ParsePacketSchedule(value, schedule))
            {
                std::fprintf(stderr,
                    "preferredtiming --schedule must be burst, adversarial, "
                    "or repair-only\n");
                return 1;
            }
            have_schedule = true;
        }
        else if (!UnknownArg("preferredtiming", argv[i])) {
            return 1;
        }
    }
    if (!have_K || !have_block_bytes || !have_preferred || !have_eviction ||
        !have_metric || !have_cache_path || !have_cache_sha256 ||
        !have_context || !have_loss || !have_seed || !have_schedule)
    {
        std::fprintf(stderr,
            "preferredtiming requires --N, --bb, --preferred-attempt, "
            "--evict-bytes, --metric, --route-cache, --route-cache-sha256, "
            "--route-context-sha256, --loss, --seed, and --schedule\n");
        return 1;
    }
    if (K < 2u || K > 64000u ||
        (block_bytes != 64u && block_bytes != 1280u && block_bytes != 4096u) ||
        preferred_attempt >= wirehair_v2::kMaxPacketSeedAttempts ||
        (have_cycle_index && cycle_index > 3u) ||
        eviction_bytes < 4096u ||
        eviction_bytes > (uint64_t)std::numeric_limits<size_t>::max() ||
        !IsLowerHexSha256(route_cache_sha256) ||
        !IsLowerHexSha256(route_context_sha256) ||
        !ValidateLoss(loss, "preferredtiming") ||
        (schedule != PacketScheduleKind::Burst &&
         schedule != PacketScheduleKind::Adversarial &&
         schedule != PacketScheduleKind::RepairOnly))
    {
        std::fprintf(stderr, "preferredtiming argument domain mismatch\n");
        return 1;
    }

    PreferredRouteCache route_cache;
    if (!LoadPreferredRouteCache(
            route_cache_path, route_cache_sha256,
            route_context_sha256, route_cache))
    {
        return 1;
    }
    const PreferredAttemptKey key(K, block_bytes);
    const PreferredRouteCache::const_iterator key_it = route_cache.find(key);
    if (key_it == route_cache.end() || key_it->second.ControlOnly) {
        std::fprintf(stderr,
            "preferredtiming route cache lacks a preferred K@bb key\n");
        return 1;
    }
    const PreferredRouteKeyRecords& route = key_it->second;
    if (route.AttemptOrder.size() != 1u ||
        route.AttemptOrder[0] != preferred_attempt)
    {
        std::fprintf(stderr,
            "preferredtiming requires an exact selected-attempt route row\n");
        return 1;
    }
    const std::map<uint32_t, PreferredRouteRecord>::const_iterator p_it =
        route.Attempts.find(preferred_attempt);
    bool direct_route = false;
    bool wide_no_op_route = false;
    if (p_it != route.Attempts.end())
    {
        const PreferredRouteRecord& record = p_it->second;
        direct_route = record.Direct && !record.Fallback && !record.NoOp &&
            record.PreferredValid &&
            record.ActualAttempt == preferred_attempt &&
            record.CanonicalAttempt == route.CanonicalAttempt &&
            preferred_attempt > route.CanonicalAttempt;
        // Development requires a direct route at bb64, but explicitly treats
        // a valid p==a0 observation at an encountered wider block size as a
        // neutral alias.  Time that alias as two identical actual-attempt
        // arms so the frozen all-width panel can retain it without inventing
        // work or substituting another K.
        wide_no_op_route = block_bytes != 64u && !record.Direct &&
            !record.Fallback && record.NoOp && record.PreferredValid &&
            record.ActualAttempt == preferred_attempt &&
            record.CanonicalAttempt == route.CanonicalAttempt &&
            preferred_attempt == route.CanonicalAttempt;
    }
    if (!direct_route && !wide_no_op_route)
    {
        std::fprintf(stderr,
            "preferredtiming requires a cached direct route or valid "
            "wide no-op alias\n");
        return 1;
    }

    wirehair_v2::PrecodeParams base_params;
    wirehair_v2::PacketRowConfig base_config;
    MakeH12Q0Configuration(K, block_bytes, base_params, base_config);
    const wirehair_v2::PrecodeParams control_params =
        wirehair_v2::PrecodeParamsForAttempt(
            base_params, route.CanonicalAttempt);
    const wirehair_v2::PrecodeParams preferred_params =
        wirehair_v2::PrecodeParamsForAttempt(base_params, preferred_attempt);
    const wirehair_v2::PacketRowConfig control_config =
        wirehair_v2::PacketConfigForAttempt(
            base_config, route.CanonicalAttempt);
    const wirehair_v2::PacketRowConfig preferred_config =
        wirehair_v2::PacketConfigForAttempt(base_config, preferred_attempt);

    wirehair_v2::PrecodeSystem control_system;
    wirehair_v2::PrecodeSystem preferred_system;
    PreferredTimingCell cell;
    if (metric == PreferredTimingMetric::Solve)
    {
        if (!wirehair_v2::BuildPrecodeSystem(
                control_params, control_system) ||
            !wirehair_v2::BuildPrecodeSystem(
                preferred_params, preferred_system) ||
            !BuildPreferredTimingCell(
                base_params, base_config, block_bytes, loss, seed,
                schedule, cell))
        {
            std::fprintf(stderr, "preferredtiming solve setup failed\n");
            return 2;
        }
    }
    else {
        std::ostringstream trace;
        trace << "wirehair-wh2-preferred-timing-setup-v1\nK=" << K
              << "\nblock_bytes=" << block_bytes << "\nseed=" << seed
              << "\nloss=" << std::setprecision(17) << loss
              << "\nschedule=" << PacketScheduleName(schedule) << '\n';
        cell.TraceSha256 = Sha256Hex(trace.str());
    }

    std::vector<uint8_t> eviction((size_t)eviction_bytes, 0u);
    EvictPreferredTimingCache(eviction);
    static const char kOrder[] = {'C', 'T', 'T', 'C', 'T', 'C', 'C', 'T'};
    const std::string cycle_index_text =
        have_cycle_index ? std::to_string(cycle_index) : "all";
    std::printf(
        "# preferredtiming: schema=v1 policy=h12-q0-adaptive metric=%s "
        "cycles=%u order=CTTCTCCT discard_cycle=0 cycle_mode=%s "
        "cycle_index=%s overhead=%s "
        "payload=%s payload_alignment=%s payload_prefaulted=%s "
        "route_cache_sha256=%s route_context_sha256=%s trace_sha256=%s\n",
        PreferredTimingMetricName(metric),
        have_cycle_index ? 1u : 4u,
        have_cycle_index ? "replacement" : "full",
        cycle_index_text.c_str(),
        metric == PreferredTimingMetric::Solve ? "4" : "none",
        metric == PreferredTimingMetric::Solve ? "distinct-zero-v1" : "none",
        metric == PreferredTimingMetric::Solve ? "64" : "none",
        metric == PreferredTimingMetric::Solve ? "1" : "none",
        route_cache_sha256.c_str(),
        route_context_sha256.c_str(), cell.TraceSha256.c_str());
    std::printf(
        "N,bb,metric,cycle,slot,arm,attempt,result,elapsed_ns,saturated,"
        "cpu_before,cpu_after,minflt_delta,majflt_delta,inactivated,"
        "binary_def,heavy_gain,block_xors,block_muladds,source_bytes,"
        "intermediate_bytes\n");

    const uint32_t first_cycle = have_cycle_index ? cycle_index : 0u;
    const uint32_t end_cycle = have_cycle_index ? cycle_index + 1u : 4u;
    for (uint32_t cycle = first_cycle; cycle < end_cycle; ++cycle)
    {
        for (uint32_t slot = 0u; slot < 8u; ++slot)
        {
            const bool control = kOrder[slot] == 'C';
            const uint32_t attempt = control ?
                route.CanonicalAttempt : preferred_attempt;
            EvictPreferredTimingCache(eviction);
            const PreferredTimingUsage usage_before =
                ReadPreferredTimingUsage();
            const int cpu_before = PreferredTimingCurrentCpu();
            const Clock::time_point begin = Clock::now();
            WirehairResult result = Wirehair_Success;
            wirehair_v2::PrecodeSolveStats stats;
            std::vector<uint8_t> intermediate;
            Clock::time_point end;
            if (metric == PreferredTimingMetric::Solve)
            {
                result = wirehair_v2::
                    SolvePrecodeSystemForValidatedSystemWithRuntime(
                        control ? control_system : preferred_system,
                        control ? control_config : preferred_config,
                        cell.Runtime, cell.Packets, block_bytes,
                        intermediate, &stats);
                end = Clock::now();
            }
            else
            {
                wirehair_v2::PrecodeSystem setup_system;
                wirehair_v2::PacketRowRuntime setup_runtime;
                const wirehair_v2::PrecodeParams& setup_params =
                    control ? control_params : preferred_params;
                const wirehair_v2::PacketRowConfig& setup_config =
                    control ? control_config : preferred_config;
                const uint64_t setup_precode_wide =
                    (uint64_t)setup_params.Staircase +
                    setup_params.DenseRows + setup_params.HeavyRows;
                if (setup_precode_wide > UINT32_MAX ||
                    !wirehair_v2::BuildPrecodeSystem(
                        setup_params, setup_system) ||
                    !setup_runtime.Initialize(
                        K, (uint32_t)setup_precode_wide,
                        setup_config.MixCount))
                {
                    result = Wirehair_Error;
                }
                // Stop while both objects are still alive: decoder setup is
                // construction plus runtime initialization, not teardown.
                end = Clock::now();
            }
            const int cpu_after = PreferredTimingCurrentCpu();
            const PreferredTimingUsage usage_after =
                ReadPreferredTimingUsage();
            bool saturated = false;
            const uint64_t elapsed_ns =
                PreferredTimingElapsedNanoseconds(begin, end, saturated);
            const int64_t minflt_delta =
                usage_before.MinorFaults < 0 || usage_after.MinorFaults < 0 ?
                -1 : usage_after.MinorFaults - usage_before.MinorFaults;
            const int64_t majflt_delta =
                usage_before.MajorFaults < 0 || usage_after.MajorFaults < 0 ?
                -1 : usage_after.MajorFaults - usage_before.MajorFaults;
            const uint32_t binary_def =
                stats.InactivatedColumns >= stats.BinaryResidualRank ?
                stats.InactivatedColumns - stats.BinaryResidualRank : 0u;
            const uint32_t heavy_gain =
                stats.ResidualRank >= stats.BinaryResidualRank ?
                stats.ResidualRank - stats.BinaryResidualRank : 0u;
            std::printf(
                "%u,%u,%s,%u,%u,%s,%u,%d,%llu,%u,%d,%d,%lld,%lld,"
                "%u,%u,%u,%llu,%llu,%zu,%zu\n",
                K, block_bytes, PreferredTimingMetricName(metric), cycle,
                slot, control ? "control" : "candidate", attempt,
                (int)result, (unsigned long long)elapsed_ns,
                saturated ? 1u : 0u, cpu_before, cpu_after,
                (long long)minflt_delta, (long long)majflt_delta,
                stats.InactivatedColumns, binary_def, heavy_gain,
                (unsigned long long)stats.BlockXors,
                (unsigned long long)stats.BlockMulAdds,
                metric == PreferredTimingMetric::Solve &&
                    result == Wirehair_Success ?
                    (size_t)K * block_bytes : 0u,
                intermediate.size());
            if (saturated || result != Wirehair_Success)
            {
                std::fprintf(stderr,
                    "preferredtiming run failed cycle=%u slot=%u result=%d\n",
                    cycle, slot, (int)result);
                return 2;
            }
        }
    }
    return 0;
}

enum class GroupedTimingCacheState
{
    Cold,
    Warm
};

const char* GroupedTimingCacheStateName(GroupedTimingCacheState state)
{
    return state == GroupedTimingCacheState::Cold ? "cold" : "warm";
}

struct GroupedTimingArm
{
    wirehair_v2::MixedCoefficientGeometry Geometry =
        wirehair_v2::MixedCoefficientGeometry::SharedCauchyX;
    uint32_t GF256Rows = wirehair_v2::kMixedGF256Rows;
    uint32_t GF16Rows = wirehair_v2::kMixedGF16Rows;
    uint32_t Period = 0u;
    uint32_t ResidueSkew = 0u;
    wirehair_v2::MixedResidueSchedule ResidueSchedule =
        wirehair_v2::MixedResidueSchedule::Constant;
    uint32_t ResidueHashSeed = 0u;
    bool IndependentExtensionResidues = false;
    uint32_t IndependentExtensionSeedXor = 78u;
    uint32_t GroupedRows = 0u;
    wirehair_v2::MixedResidueBucketMode Buckets =
        wirehair_v2::MixedResidueBucketMode::Automatic;
    uint32_t PacketSeedMultiplier = 1u;
    bool PacketSeedAvalanche = false;
    uint32_t OddPacketPeelSeedXor = 0u;
    uint32_t DenseRowsOverride = 0u;
    bool DenseTwoAnchor = true;
    uint32_t DenseTwoAnchorPhase = 0u;
    uint32_t PacketMixCount = 2u;
};

const char* GroupedTimingRhsRouteName(
    wirehair_v2::MixedCompletionRhsRoute route)
{
    switch (route)
    {
    case wirehair_v2::MixedCompletionRhsRoute::NotReached:
        return "not-reached";
    case wirehair_v2::MixedCompletionRhsRoute::Streamed:
        return "streamed";
    case wirehair_v2::MixedCompletionRhsRoute::Dual:
        return "dual";
    case wirehair_v2::MixedCompletionRhsRoute::JointDelta:
        return "joint-delta";
    }
    return "invalid";
}

bool ConfigureGroupedTimingArm(const GroupedTimingArm& arm)
{
    // Reapply the complete thread-local experiment state for every physical
    // solve.  Several setters deliberately clear secondary schedules, so the
    // grouped suffix must remain last.  None of this setup is timed.
    if (!wirehair_v2::SetMixedCoefficientGeometryForTesting(arm.Geometry) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(arm.GF16Rows) ||
        !wirehair_v2::SetMixedCoefficientPeriodForTesting(arm.Period) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(arm.GF256Rows) ||
        !wirehair_v2::SetMixedResidueSkewForTesting(arm.ResidueSkew) ||
        !wirehair_v2::SetMixedResidueScheduleForTesting(arm.ResidueSchedule))
    {
        return false;
    }
    wirehair_v2::SetMixedResidueHashSeedForTesting(arm.ResidueHashSeed);
    wirehair_v2::SetMixedIndependentExtensionSeedXorForTesting(
        arm.IndependentExtensionSeedXor);
    if (!wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(
            arm.IndependentExtensionResidues) ||
        !wirehair_v2::SetMixedResidueBucketModeForTesting(arm.Buckets) ||
        !wirehair_v2::SetPacketRowSeedMultiplierForTesting(
            arm.PacketSeedMultiplier))
    {
        return false;
    }
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(
        arm.PacketSeedAvalanche);
    wirehair_v2::SetOddPacketPeelSeedXorForTesting(
        arm.OddPacketPeelSeedXor);
    if (!wirehair_v2::SetMixedGroupedGF256RowsForTesting(
            arm.GroupedRows))
    {
        return false;
    }
    return wirehair_v2::ActiveMixedCoefficientGeometry() == arm.Geometry &&
        wirehair_v2::ActiveMixedGF256Rows() == arm.GF256Rows &&
        wirehair_v2::ActiveMixedGF16Rows() == arm.GF16Rows &&
        wirehair_v2::ActiveMixedCoefficientPeriod() == arm.Period &&
        wirehair_v2::ActiveMixedResidueSkew() == arm.ResidueSkew &&
        wirehair_v2::ActiveMixedResidueSchedule() == arm.ResidueSchedule &&
        wirehair_v2::ActiveMixedResidueHashSeed() == arm.ResidueHashSeed &&
        wirehair_v2::ActiveMixedIndependentExtensionResidues() ==
            arm.IndependentExtensionResidues &&
        wirehair_v2::ActiveMixedIndependentExtensionSeedXorForTesting() ==
            arm.IndependentExtensionSeedXor &&
        wirehair_v2::ActiveMixedGroupedGF256Rows() == arm.GroupedRows &&
        wirehair_v2::ActiveMixedResidueBucketModeForTesting() == arm.Buckets &&
        wirehair_v2::ActivePacketRowSeedMultiplierForTesting() ==
            arm.PacketSeedMultiplier &&
        wirehair_v2::ActivePacketRowSeedAvalancheForTesting() ==
            arm.PacketSeedAvalanche &&
        wirehair_v2::ActiveOddPacketPeelSeedXorForTesting() ==
            arm.OddPacketPeelSeedXor;
}

wirehair_v2::MixedCompletionRhsRoute GroupedTimingExpectedRhsRoute(
    const GroupedTimingArm& arm,
    const wirehair_v2::PrecodeSystem& system,
    uint32_t block_bytes)
{
    const uint64_t column_count =
        (uint64_t)system.Params.BlockCount + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint64_t first_heavy_column =
        column_count - system.Params.HeavyRows;
    // The solver intentionally collapses grouped C back to A while every
    // non-heavy column still lies in residue block zero.
    const bool effective_grouped =
        arm.GroupedRows != 0u && first_heavy_column > arm.Period;
    const bool secondary =
        arm.IndependentExtensionResidues || effective_grouped;
    if (!secondary) {
        return wirehair_v2::MixedCompletionRhsRoute::Streamed;
    }
    const bool automatic_joint =
        wirehair_v2::UseAutomaticMixedJointResidueBucketsForTesting(
            system.Params.BlockCount, block_bytes, arm.Period);
    const bool request_joint =
        arm.Buckets == wirehair_v2::MixedResidueBucketMode::JointDelta ||
        (arm.Buckets == wirehair_v2::MixedResidueBucketMode::Automatic &&
         automatic_joint);
    if (request_joint && wirehair_v2::MixedJointResidueBucketStorageFits(
            arm.Period, block_bytes,
            wirehair_v2::kMixedJointResidueBucketDataByteCap))
    {
        return wirehair_v2::MixedCompletionRhsRoute::JointDelta;
    }
    const uint64_t one_plane_bytes = (uint64_t)arm.Period * block_bytes;
    const bool dual_requested =
        arm.Buckets == wirehair_v2::MixedResidueBucketMode::Dual ||
        (arm.Buckets == wirehair_v2::MixedResidueBucketMode::Automatic &&
         !automatic_joint && column_count >= 30000u &&
         block_bytes >= 1024u &&
         one_plane_bytes * 2u <= (UINT64_C(128) << 10));
    if (dual_requested &&
        one_plane_bytes <=
            wirehair_v2::kMixedJointResidueBucketDataByteCap / 2u)
    {
        return wirehair_v2::MixedCompletionRhsRoute::Dual;
    }
    return wirehair_v2::MixedCompletionRhsRoute::Streamed;
}

bool SameGroupedTimingBaseGraph(
    const wirehair_v2::PrecodeParams& a_params,
    const wirehair_v2::PacketRowConfig& a_config,
    const wirehair_v2::PrecodeParams& b_params,
    const wirehair_v2::PacketRowConfig& b_config)
{
    return a_params.BlockCount == b_params.BlockCount &&
        a_params.Staircase == b_params.Staircase &&
        a_params.DenseRows == b_params.DenseRows &&
        a_params.HeavyRows == b_params.HeavyRows &&
        a_params.SourceHits == b_params.SourceHits &&
        a_params.Field == b_params.Field &&
        a_params.DenseIdentityCorner == b_params.DenseIdentityCorner &&
        a_params.DenseTwoAnchor == b_params.DenseTwoAnchor &&
        a_params.DenseTwoAnchorPhase == b_params.DenseTwoAnchorPhase &&
        a_params.HeavyFamily == b_params.HeavyFamily &&
        a_params.Seed == b_params.Seed &&
        a_config.PeelSeed == b_config.PeelSeed &&
        a_config.MixCount == b_config.MixCount;
}

struct GroupedTimingObservation
{
    WirehairResult Result = Wirehair_Error;
    wirehair_v2::PrecodeSolveStats Stats;
};

bool RunGroupedTimingObservation(
    const GroupedTimingArm& arm,
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const PreferredTimingCell& cell,
    uint32_t block_bytes,
    GroupedTimingObservation& observation)
{
    if (!ConfigureGroupedTimingArm(arm)) return false;
    wirehair_v2::SolveValueStorage intermediate;
    observation = GroupedTimingObservation{};
    observation.Result =
        wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
            system, config, cell.Runtime, cell.Packets, block_bytes,
            intermediate, &observation.Stats);
    const uint64_t expected_arena_bytes =
        ((uint64_t)system.Params.BlockCount + system.Params.Staircase +
            system.Params.DenseRows + system.Params.HeavyRows) * block_bytes;
    const bool valid = (observation.Result == Wirehair_Success ||
            observation.Result == Wirehair_NeedMore) &&
        observation.Stats.SolveValueArenaBytes == expected_arena_bytes &&
        observation.Stats.SolveValueArenaEagerZeroBytes == 0u &&
        observation.Stats.SolveValueArenaCommitCopyBytes == 0u &&
        (observation.Result != Wirehair_Success ||
         observation.Stats.MixedRhsRoute ==
            GroupedTimingExpectedRhsRoute(arm, system, block_bytes));
    if (!valid) {
        std::fprintf(stderr,
            "groupedtiming solve-arena preflight receipt mismatch "
            "result=%d arena=%llu expected=%llu eager_zero=%llu copy=%llu "
            "rhs_route=%s expected_route=%s\n",
            (int)observation.Result,
            (unsigned long long)observation.Stats.SolveValueArenaBytes,
            (unsigned long long)expected_arena_bytes,
            (unsigned long long)
                observation.Stats.SolveValueArenaEagerZeroBytes,
            (unsigned long long)
                observation.Stats.SolveValueArenaCommitCopyBytes,
            GroupedTimingRhsRouteName(observation.Stats.MixedRhsRoute),
            GroupedTimingRhsRouteName(
                GroupedTimingExpectedRhsRoute(arm, system, block_bytes)));
    }
    return valid;
}

const char* GroupedTimingOutcomeClass(
    WirehairResult control,
    WirehairResult candidate)
{
    const bool control_success = control == Wirehair_Success;
    const bool candidate_success = candidate == Wirehair_Success;
    if (control_success && candidate_success) return "common-success";
    if (control_success) return "control-only";
    if (candidate_success) return "candidate-only";
    return "common-failure";
}

int GroupedTimingCpuMigrated(int before, int after)
{
    return before < 0 || after < 0 ? -1 : (before != after ? 1 : 0);
}

int GroupedTimingFaultContaminated(
    int64_t minor_delta,
    int64_t major_delta)
{
    return minor_delta < 0 || major_delta < 0 ?
        -1 : (minor_delta != 0 || major_delta != 0 ? 1 : 0);
}

int CmdGroupedTiming(int argc, char** argv)
{
    uint32_t K = 0u;
    uint32_t block_bytes = 0u;
    uint32_t overhead = 0u;
    uint32_t cycle_index = 0u;
    uint64_t eviction_bytes = 0u;
    double loss = 0.0;
    uint64_t seed = 0u;
    PacketScheduleKind schedule = PacketScheduleKind::Burst;
    GroupedTimingCacheState cache_state = GroupedTimingCacheState::Cold;
    GroupedTimingArm control_arm;
    GroupedTimingArm candidate_arm;
    bool have_K = false;
    bool have_block_bytes = false;
    bool have_overhead = false;
    bool have_control_geometry = false;
    bool have_control_period = false;
    bool have_control_rows = false;
    bool have_control_buckets = false;
    bool have_candidate_geometry = false;
    bool have_candidate_period = false;
    bool have_candidate_rows = false;
    bool have_candidate_buckets = false;
    bool have_eviction = false;
    bool have_cache_state = false;
    bool have_loss = false;
    bool have_seed = false;
    bool have_schedule = false;
    bool have_cycle_index = false;

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N")) {
            if (have_K ||
                !TakeArg("groupedtiming", "--N", argc, argv, i, value) ||
                !ParseCanonicalU32Arg("--N", value, K)) return 1;
            have_K = true;
        }
        else if (!std::strcmp(argv[i], "--bb")) {
            if (have_block_bytes ||
                !TakeArg("groupedtiming", "--bb", argc, argv, i, value) ||
                !ParseCanonicalU32Arg("--bb", value, block_bytes)) return 1;
            have_block_bytes = true;
        }
        else if (!std::strcmp(argv[i], "--overhead")) {
            if (have_overhead || !TakeArg(
                    "groupedtiming", "--overhead", argc, argv, i, value) ||
                !ParseCanonicalU32Arg("--overhead", value, overhead))
            {
                return 1;
            }
            have_overhead = true;
        }
        else if (!std::strcmp(argv[i], "--control-period")) {
            if (have_control_period || !TakeArg(
                    "groupedtiming", "--control-period",
                    argc, argv, i, value) ||
                !ParseCanonicalU32Arg(
                    "--control-period", value, control_arm.Period))
            {
                return 1;
            }
            have_control_period = true;
        }
        else if (!std::strcmp(argv[i], "--control-geometry")) {
            if (have_control_geometry || !TakeArg(
                    "groupedtiming", "--control-geometry",
                    argc, argv, i, value) ||
                !ParseMixedCoefficientGeometry(value, control_arm.Geometry))
            {
                std::fprintf(stderr,
                    "groupedtiming --control-geometry must be frozen or "
                    "shared-x\n");
                return 1;
            }
            have_control_geometry = true;
        }
        else if (!std::strcmp(argv[i], "--control-grouped-rows")) {
            if (have_control_rows || !TakeArg(
                    "groupedtiming", "--control-grouped-rows",
                    argc, argv, i, value) ||
                !ParseCanonicalU32Arg(
                    "--control-grouped-rows", value,
                    control_arm.GroupedRows))
            {
                return 1;
            }
            have_control_rows = true;
        }
        else if (!std::strcmp(argv[i], "--control-buckets")) {
            if (have_control_buckets || !TakeArg(
                    "groupedtiming", "--control-buckets",
                    argc, argv, i, value) ||
                !ParseMixedResidueBucketMode(value, control_arm.Buckets))
            {
                std::fprintf(stderr,
                    "groupedtiming --control-buckets must be auto, "
                    "separate, dual, or joint-delta\n");
                return 1;
            }
            have_control_buckets = true;
        }
        else if (!std::strcmp(argv[i], "--candidate-period")) {
            if (have_candidate_period || !TakeArg(
                    "groupedtiming", "--candidate-period",
                    argc, argv, i, value) ||
                !ParseCanonicalU32Arg(
                    "--candidate-period", value, candidate_arm.Period))
            {
                return 1;
            }
            have_candidate_period = true;
        }
        else if (!std::strcmp(argv[i], "--candidate-geometry")) {
            if (have_candidate_geometry || !TakeArg(
                    "groupedtiming", "--candidate-geometry",
                    argc, argv, i, value) ||
                !ParseMixedCoefficientGeometry(value, candidate_arm.Geometry))
            {
                std::fprintf(stderr,
                    "groupedtiming --candidate-geometry must be frozen or "
                    "shared-x\n");
                return 1;
            }
            have_candidate_geometry = true;
        }
        else if (!std::strcmp(argv[i], "--candidate-grouped-rows")) {
            if (have_candidate_rows || !TakeArg(
                    "groupedtiming", "--candidate-grouped-rows",
                    argc, argv, i, value) ||
                !ParseCanonicalU32Arg(
                    "--candidate-grouped-rows", value,
                    candidate_arm.GroupedRows))
            {
                return 1;
            }
            have_candidate_rows = true;
        }
        else if (!std::strcmp(argv[i], "--candidate-buckets")) {
            if (have_candidate_buckets || !TakeArg(
                    "groupedtiming", "--candidate-buckets",
                    argc, argv, i, value) ||
                !ParseMixedResidueBucketMode(value, candidate_arm.Buckets))
            {
                std::fprintf(stderr,
                    "groupedtiming --candidate-buckets must be auto, "
                    "separate, dual, or joint-delta\n");
                return 1;
            }
            have_candidate_buckets = true;
        }
        else if (!std::strcmp(argv[i], "--evict-bytes")) {
            if (have_eviction || !TakeArg(
                    "groupedtiming", "--evict-bytes", argc, argv, i,
                    value) || !ParseCanonicalU64Arg(
                        "--evict-bytes", value, eviction_bytes))
            {
                return 1;
            }
            have_eviction = true;
        }
        else if (!std::strcmp(argv[i], "--cache-state")) {
            if (have_cache_state || !TakeArg(
                    "groupedtiming", "--cache-state", argc, argv, i,
                    value))
            {
                return 1;
            }
            if (!std::strcmp(value, "cold")) {
                cache_state = GroupedTimingCacheState::Cold;
            }
            else if (!std::strcmp(value, "warm")) {
                cache_state = GroupedTimingCacheState::Warm;
            }
            else {
                std::fprintf(stderr,
                    "groupedtiming --cache-state must be cold or warm\n");
                return 1;
            }
            have_cache_state = true;
        }
        else if (!std::strcmp(argv[i], "--cycle-index")) {
            if (have_cycle_index || !TakeArg(
                    "groupedtiming", "--cycle-index", argc, argv, i,
                    value) || !ParseCanonicalU32Arg(
                        "--cycle-index", value, cycle_index))
            {
                return 1;
            }
            have_cycle_index = true;
        }
        else if (!std::strcmp(argv[i], "--loss")) {
            if (have_loss || !TakeArg(
                    "groupedtiming", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, loss)) return 1;
            have_loss = true;
        }
        else if (!std::strcmp(argv[i], "--seed")) {
            if (have_seed || !TakeArg(
                    "groupedtiming", "--seed", argc, argv, i, value) ||
                !ParseCanonicalU64Arg("--seed", value, seed)) return 1;
            have_seed = true;
        }
        else if (!std::strcmp(argv[i], "--schedule")) {
            if (have_schedule || !TakeArg(
                    "groupedtiming", "--schedule", argc, argv, i, value) ||
                !ParsePacketSchedule(value, schedule))
            {
                std::fprintf(stderr,
                    "groupedtiming --schedule must be burst, adversarial, "
                    "or repair-only\n");
                return 1;
            }
            have_schedule = true;
        }
        else if (!UnknownArg("groupedtiming", argv[i])) {
            return 1;
        }
    }
    if (!have_K || !have_block_bytes || !have_overhead ||
        !have_control_geometry || !have_control_period || !have_control_rows ||
        !have_control_buckets || !have_candidate_period ||
        !have_candidate_geometry || !have_candidate_rows ||
        !have_candidate_buckets || !have_eviction ||
        !have_cache_state || !have_loss || !have_seed || !have_schedule)
    {
        std::fprintf(stderr,
            "groupedtiming requires --N, --bb, --overhead, "
            "--control-geometry, --control-period, --control-grouped-rows, "
            "--control-buckets, --candidate-geometry, --candidate-period, "
            "--candidate-grouped-rows, "
            "--candidate-buckets, --evict-bytes, --cache-state, --loss, "
            "--seed, and --schedule\n");
        return 1;
    }
    const auto invalid_arm = [](const GroupedTimingArm& arm) {
        return arm.Period <
                arm.GF256Rows + arm.GF16Rows ||
            arm.Period > wirehair_v2::kMixedCoefficientPeriod ||
            (arm.Geometry ==
                    wirehair_v2::MixedCoefficientGeometry::FrozenPowerX &&
             arm.GF256Rows != wirehair_v2::kMixedGF256Rows) ||
            arm.GroupedRows > 9u ||
            (arm.GroupedRows != 0u &&
             arm.Period <= arm.GF256Rows + arm.GF16Rows) ||
            (arm.GroupedRows != 0u &&
             arm.Geometry !=
                wirehair_v2::MixedCoefficientGeometry::SharedCauchyX) ||
            (arm.GroupedRows == 0u &&
             arm.Buckets !=
                wirehair_v2::MixedResidueBucketMode::Automatic) ||
            arm.ResidueSkew != 0u ||
            arm.ResidueSchedule !=
                wirehair_v2::MixedResidueSchedule::Constant ||
            arm.ResidueHashSeed != 0u ||
            arm.IndependentExtensionResidues ||
            arm.IndependentExtensionSeedXor != 78u ||
            arm.PacketSeedMultiplier != 1u ||
            arm.PacketSeedAvalanche ||
            arm.OddPacketPeelSeedXor != 0u ||
            arm.DenseRowsOverride != 0u || !arm.DenseTwoAnchor ||
            arm.DenseTwoAnchorPhase != 0u || arm.PacketMixCount != 2u;
    };
    if (K < 2u || K > 64000u ||
        (block_bytes != 64u && block_bytes != 256u &&
         block_bytes != 512u && block_bytes != 1024u &&
         block_bytes != 1280u && block_bytes != 4096u) ||
        overhead > 1024u || (uint64_t)K + overhead > UINT32_MAX ||
        invalid_arm(control_arm) || invalid_arm(candidate_arm) ||
        (have_cycle_index && cycle_index > 3u) ||
        eviction_bytes < 4096u ||
        eviction_bytes > (uint64_t)std::numeric_limits<size_t>::max() ||
        !ValidateLoss(loss, "groupedtiming") ||
        (schedule != PacketScheduleKind::Burst &&
         schedule != PacketScheduleKind::Adversarial &&
         schedule != PacketScheduleKind::RepairOnly))
    {
        std::fprintf(stderr, "groupedtiming argument domain mismatch\n");
        return 1;
    }
    const std::vector<int> one_K(1u, (int)K);
    const std::vector<int> one_bb(1u, (int)block_bytes);
    if (!ValidatePayloadE2EInputs(
            one_K, one_bb, "groupedtiming"))
    {
        return 1;
    }

    wirehair_v2::PrecodeParams control_base_params;
    wirehair_v2::PacketRowConfig control_base_config;
    wirehair_v2::PrecodeParams candidate_base_params;
    wirehair_v2::PacketRowConfig candidate_base_config;
    wirehair_v2::PrecodeSystem control_system;
    wirehair_v2::PacketRowConfig control_config;
    wirehair_v2::PrecodeSystem candidate_system;
    wirehair_v2::PacketRowConfig candidate_config;
    uint32_t control_attempt = 0u;
    uint32_t candidate_attempt = 0u;
    uint32_t control_grouped_hash_seed = 0u;
    uint32_t candidate_grouped_hash_seed = 0u;
    uint32_t control_extension_hash_seed = 0u;
    uint32_t candidate_extension_hash_seed = 0u;
    if (!ConfigureGroupedTimingArm(control_arm) ||
        !MakeH12Q0Configuration(
            K, block_bytes, control_base_params, control_base_config))
    {
        std::fprintf(stderr,
            "groupedtiming could not construct the control arm\n");
        return 2;
    }
    // The grouped reliability campaign used the raw two-anchor architecture
    // at every K, rather than the later adaptive named-profile cutoff.
    control_base_params.DenseTwoAnchor = control_arm.DenseTwoAnchor;
    control_base_params.DenseTwoAnchorPhase =
        control_arm.DenseTwoAnchorPhase;
    if (control_arm.DenseRowsOverride != 0u) {
        control_base_params.DenseRows = control_arm.DenseRowsOverride;
    }
    control_base_config.MixCount = control_arm.PacketMixCount;
    control_extension_hash_seed =
        wirehair_v2::ActiveMixedExtensionResidueHashSeed();
    if (control_arm.GroupedRows != 0u) {
        control_grouped_hash_seed =
            wirehair_v2::ActiveMixedGroupedGF256HashSeed();
    }
    const WirehairResult control_select =
        wirehair_v2::SelectSystematicConfiguration(
            control_base_params, control_base_config,
            control_system, control_config, &control_attempt);
    if (control_select != Wirehair_Success) {
        std::fprintf(stderr,
            "groupedtiming control seed selection failed result=%d\n",
            (int)control_select);
        return 2;
    }
    if (!ConfigureGroupedTimingArm(candidate_arm) ||
        !MakeH12Q0Configuration(
            K, block_bytes, candidate_base_params, candidate_base_config))
    {
        std::fprintf(stderr,
            "groupedtiming could not construct the candidate arm\n");
        return 2;
    }
    candidate_base_params.DenseTwoAnchor = candidate_arm.DenseTwoAnchor;
    candidate_base_params.DenseTwoAnchorPhase =
        candidate_arm.DenseTwoAnchorPhase;
    if (candidate_arm.DenseRowsOverride != 0u) {
        candidate_base_params.DenseRows = candidate_arm.DenseRowsOverride;
    }
    candidate_base_config.MixCount = candidate_arm.PacketMixCount;
    candidate_extension_hash_seed =
        wirehair_v2::ActiveMixedExtensionResidueHashSeed();
    if (candidate_arm.GroupedRows != 0u) {
        candidate_grouped_hash_seed =
            wirehair_v2::ActiveMixedGroupedGF256HashSeed();
    }
    if (!SameGroupedTimingBaseGraph(
            control_base_params, control_base_config,
            candidate_base_params, candidate_base_config))
    {
        std::fprintf(stderr,
            "groupedtiming arms do not share one base graph/configuration\n");
        return 2;
    }
    const WirehairResult candidate_select =
        wirehair_v2::SelectSystematicConfiguration(
            candidate_base_params, candidate_base_config,
            candidate_system, candidate_config, &candidate_attempt);
    if (candidate_select != Wirehair_Success) {
        std::fprintf(stderr,
            "groupedtiming candidate seed selection failed result=%d\n",
            (int)candidate_select);
        return 2;
    }

    PreferredTimingCell cell;
    if (!BuildFullPayloadTimingCell(
            control_system.Params, control_config, block_bytes, overhead, loss,
            seed, schedule, true,
            "wirehair-wh2-grouped-timing-trace-v1", cell))
    {
        std::fprintf(stderr, "groupedtiming solve setup failed\n");
        return 2;
    }
    const uint64_t candidate_precode_wide =
        (uint64_t)candidate_system.Params.Staircase +
        candidate_system.Params.DenseRows + candidate_system.Params.HeavyRows;
    if (candidate_precode_wide > UINT32_MAX ||
        !cell.Runtime.IsValidFor(
            K, (uint32_t)candidate_precode_wide,
            candidate_config.MixCount) ||
        cell.Packets.size() != (size_t)K + overhead)
    {
        std::fprintf(stderr,
            "groupedtiming selected arms do not share one packet domain\n");
        return 2;
    }

    // Allocate and prefault the eviction working set before either preflight.
    // In warm mode, touching it afterward would make the first recorded solve
    // cold; in cold mode, each timed slot explicitly traverses it again.
    std::vector<uint8_t> eviction((size_t)eviction_bytes, 0u);
    EvictPreferredTimingCache(eviction);

    GroupedTimingObservation control_preflight;
    GroupedTimingObservation candidate_preflight;
    if (!RunGroupedTimingObservation(
            control_arm, control_system, control_config, cell, block_bytes,
            control_preflight) ||
        !RunGroupedTimingObservation(
            candidate_arm, candidate_system, candidate_config, cell,
            block_bytes,
            candidate_preflight))
    {
        std::fprintf(stderr, "groupedtiming preflight failed\n");
        return 2;
    }
    const char* const outcome_class = GroupedTimingOutcomeClass(
        control_preflight.Result, candidate_preflight.Result);
    const bool common_success =
        control_preflight.Result == Wirehair_Success &&
        candidate_preflight.Result == Wirehair_Success;
    const wirehair_v2::MixedCompletionRhsRoute control_expected_route =
        GroupedTimingExpectedRhsRoute(
            control_arm, control_system, block_bytes);
    const wirehair_v2::MixedCompletionRhsRoute candidate_expected_route =
        GroupedTimingExpectedRhsRoute(
            candidate_arm, candidate_system, block_bytes);
    const auto arm_receipt = [](
        const char* prefix,
        const GroupedTimingArm& arm,
        const wirehair_v2::PrecodeSystem& system,
        const wirehair_v2::PacketRowConfig& config,
        uint32_t grouped_hash_seed,
        uint32_t extension_hash_seed,
        wirehair_v2::MixedCompletionRhsRoute expected_route,
        const GroupedTimingObservation& preflight)
    {
        std::ostringstream out;
        out << prefix << "_geometry=" <<
                MixedCoefficientGeometryName(arm.Geometry)
            << ' ' << prefix << "_gf256_rows=" << arm.GF256Rows
            << ' ' << prefix << "_gf16_rows=" << arm.GF16Rows
            << ' ' << prefix << "_residue_skew=" << arm.ResidueSkew
            << ' ' << prefix << "_residue_schedule=" <<
                MixedResidueScheduleName(arm.ResidueSchedule)
            << ' ' << prefix << "_residue_hash_seed=0x" << std::hex <<
                arm.ResidueHashSeed << std::dec
            << ' ' << prefix << "_residues_rotated=" <<
                (arm.ResidueSkew != 0u ||
                 arm.ResidueSchedule !=
                    wirehair_v2::MixedResidueSchedule::Constant ? 1u : 0u)
            << ' ' << prefix << "_independent_extension_residues=" <<
                (arm.IndependentExtensionResidues ? 1u : 0u)
            << ' ' << prefix << "_independent_extension_seed_xor=0x" <<
                std::hex << arm.IndependentExtensionSeedXor << std::dec
            << ' ' << prefix << "_extension_hash_seed=0x" << std::hex <<
                extension_hash_seed << std::dec
            << ' ' << prefix << "_grouped_hash_seed_exact=0x" << std::hex <<
                grouped_hash_seed << std::dec
            << ' ' << prefix << "_packet_seed_multiplier=" <<
                arm.PacketSeedMultiplier
            << ' ' << prefix << "_packet_seed_avalanche=" <<
                (arm.PacketSeedAvalanche ? 1u : 0u)
            << ' ' << prefix << "_odd_packet_peel_seed_xor=0x" << std::hex <<
                arm.OddPacketPeelSeedXor << std::dec
            << ' ' << prefix << "_dense_rows_override=" <<
                arm.DenseRowsOverride
            << ' ' << prefix << "_dense_two_anchor_exact=" <<
                (arm.DenseTwoAnchor ? 1u : 0u)
            << ' ' << prefix << "_dense_two_anchor_phase=" <<
                arm.DenseTwoAnchorPhase
            << ' ' << prefix << "_staircase_rows=" <<
                system.Params.Staircase
            << ' ' << prefix << "_dense_rows=" << system.Params.DenseRows
            << ' ' << prefix << "_heavy_rows=" << system.Params.HeavyRows
            << ' ' << prefix << "_source_hits=" << system.Params.SourceHits
            << ' ' << prefix << "_field=" <<
                static_cast<uint32_t>(system.Params.Field)
            << ' ' << prefix << "_dense_identity_corner=" <<
                (system.Params.DenseIdentityCorner ? 1u : 0u)
            << ' ' << prefix << "_heavy_family=" <<
                static_cast<uint32_t>(system.Params.HeavyFamily)
            << ' ' << prefix << "_mix_count=" << config.MixCount
            << ' ' << prefix << "_rhs_route_expected=" <<
                GroupedTimingRhsRouteName(expected_route)
            << ' ' << prefix << "_preflight_rhs_route=" <<
                GroupedTimingRhsRouteName(preflight.Stats.MixedRhsRoute);
        return out.str();
    };
    const std::string exact_arm_receipts =
        arm_receipt(
            "control", control_arm, control_system, control_config,
            control_grouped_hash_seed, control_extension_hash_seed,
            control_expected_route, control_preflight) + " " +
        arm_receipt(
            "candidate", candidate_arm, candidate_system, candidate_config,
            candidate_grouped_hash_seed, candidate_extension_hash_seed,
            candidate_expected_route, candidate_preflight);

    static const char kOrder[] = {'A', 'B', 'B', 'A', 'B', 'A', 'A', 'B'};
    const std::string cycle_index_text =
        have_cycle_index ? std::to_string(cycle_index) : "all";
    const uint64_t packet_payload_bytes =
        ((uint64_t)K + overhead) * block_bytes;
    std::printf(
        "# groupedtiming: schema=v2 policy=h12-q0-grouped "
        "timing_scope=solve cycles=%u order=ABBABAAB discard_cycle=0 "
        "cycle_mode=%s cycle_index=%s N=%u bb=%u overhead=%u "
        "loss=%.17g seed=%llu schedule=%s cache_state=%s "
        "overhead_stream=salted "
        "evict_bytes=%llu eviction_prefaulted=1 "
        "control_period=%u control_grouped_rows=%u "
        "control_buckets=%s control_grouped_hash_seed=0x%x "
        "control_final_h_a_columns=%u candidate_period=%u "
        "candidate_grouped_rows=%u candidate_buckets=%s "
        "candidate_grouped_hash_seed=0x%x "
        "candidate_final_h_a_columns=%u gf256_rows=10 gf16_rows=2 "
        "dense_two_anchor=1 control_attempt=%u "
        "control_matrix_seed=0x%llx control_peel_seed=0x%x "
        "candidate_attempt=%u candidate_matrix_seed=0x%llx "
        "candidate_peel_seed=0x%x mix=2 "
        "payload=distinct-packet-zero-v1 payload_count=%u "
        "payload_bytes=%llu payload_alignment=64 payload_prefaulted=1 "
        "system_build=outside-timer tls_reapply=full-per-slot-outside-timer "
        "allocator_tls_state=preflight-warmed "
        "solve_value_storage=owned-noinit solve_value_publish=swap "
        "preflight_control_result=%d preflight_candidate_result=%d "
        "cell_class=%s common_success=%u trace_sha256=%s %s\n",
        have_cycle_index ? 1u : 4u,
        have_cycle_index ? "replacement" : "full",
        cycle_index_text.c_str(), K, block_bytes, overhead, loss,
        (unsigned long long)seed, PacketScheduleName(schedule),
        GroupedTimingCacheStateName(cache_state),
        (unsigned long long)eviction_bytes,
        control_arm.Period, control_arm.GroupedRows,
        MixedResidueBucketModeName(control_arm.Buckets),
        control_grouped_hash_seed,
        control_arm.GroupedRows != 0u ?
            control_system.Params.HeavyRows : 0u,
        candidate_arm.Period, candidate_arm.GroupedRows,
        MixedResidueBucketModeName(candidate_arm.Buckets),
        candidate_grouped_hash_seed,
        candidate_arm.GroupedRows != 0u ?
            candidate_system.Params.HeavyRows : 0u,
        control_attempt,
        (unsigned long long)control_system.Params.Seed,
        control_config.PeelSeed, candidate_attempt,
        (unsigned long long)candidate_system.Params.Seed,
        candidate_config.PeelSeed,
        K + overhead,
        (unsigned long long)packet_payload_bytes,
        (int)control_preflight.Result, (int)candidate_preflight.Result,
        outcome_class, common_success ? 1u : 0u,
        cell.TraceSha256.c_str(), exact_arm_receipts.c_str());
    std::printf(
        "N,bb,overhead,schedule,seed,loss,cache_state,cycle,slot,arm,"
        "period,grouped_rows,buckets_requested,seed_attempt,matrix_seed,"
        "peel_seed,preflight_result,cell_class,common_success,result,"
        "outcome_stable,elapsed_ns,saturated,"
        "cpu_before,cpu_after,cpu_migrated,minflt_delta,majflt_delta,"
        "fault_contaminated,inactivated,binary_def,heavy_gain,block_xors,"
        "block_muladds,build_ns,peel_ns,project_ns,residual_ns,backsub_ns,"
        "joint_source_xors,joint_marginal_xors,joint_marginal_copies,"
        "joint_active_deltas,joint_scratch_bytes,dual_source_columns,"
        "source_bytes,packet_payload_bytes,intermediate_bytes,"
        "solve_value_arena_bytes,solve_value_eager_zero_bytes,"
        "solve_value_commit_copy_bytes,geometry,gf256_rows,gf16_rows,"
        "residue_skew,residue_schedule,residue_hash_seed,residues_rotated,"
        "independent_extension_residues,independent_extension_seed_xor,"
        "extension_hash_seed,grouped_hash_seed_exact,packet_seed_multiplier,"
        "packet_seed_avalanche,odd_packet_peel_seed_xor,dense_rows_override,"
        "dense_two_anchor_exact,dense_two_anchor_phase,staircase_rows,"
        "dense_rows,heavy_rows,source_hits,field,dense_identity_corner,"
        "heavy_family,mix_count,rhs_route_expected,rhs_route_actual\n");

    const uint32_t first_cycle = have_cycle_index ? cycle_index : 0u;
    const uint32_t end_cycle = have_cycle_index ? cycle_index + 1u : 4u;
    for (uint32_t cycle = first_cycle; cycle < end_cycle; ++cycle)
    {
        for (uint32_t slot = 0u; slot < 8u; ++slot)
        {
            const bool control = kOrder[slot] == 'A';
            const GroupedTimingArm& arm = control ?
                control_arm : candidate_arm;
            const wirehair_v2::PrecodeSystem& active_system = control ?
                control_system : candidate_system;
            const wirehair_v2::PacketRowConfig& active_config = control ?
                control_config : candidate_config;
            const uint32_t active_attempt = control ?
                control_attempt : candidate_attempt;
            const WirehairResult preflight_result = control ?
                control_preflight.Result : candidate_preflight.Result;
            const wirehair_v2::MixedCompletionRhsRoute preflight_route =
                control ? control_preflight.Stats.MixedRhsRoute :
                    candidate_preflight.Stats.MixedRhsRoute;
            const wirehair_v2::MixedCompletionRhsRoute expected_route =
                control ? control_expected_route : candidate_expected_route;
            const uint32_t active_grouped_hash_seed = control ?
                control_grouped_hash_seed : candidate_grouped_hash_seed;
            const uint32_t active_extension_hash_seed = control ?
                control_extension_hash_seed : candidate_extension_hash_seed;
            if (!ConfigureGroupedTimingArm(arm)) {
                std::fprintf(stderr,
                    "groupedtiming TLS reapplication failed cycle=%u "
                    "slot=%u\n", cycle, slot);
                return 2;
            }
            if (cache_state == GroupedTimingCacheState::Cold) {
                EvictPreferredTimingCache(eviction);
            }
            const PreferredTimingUsage usage_before =
                ReadPreferredTimingUsage();
            const int cpu_before = PreferredTimingCurrentCpu();
            const Clock::time_point begin = Clock::now();
            wirehair_v2::PrecodeSolveStats stats;
            wirehair_v2::SolveValueStorage intermediate;
            const WirehairResult result =
                wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
                    active_system, active_config,
                    cell.Runtime, cell.Packets, block_bytes,
                    intermediate, &stats);
            const Clock::time_point end = Clock::now();
            const int cpu_after = PreferredTimingCurrentCpu();
            const PreferredTimingUsage usage_after =
                ReadPreferredTimingUsage();
            bool saturated = false;
            const uint64_t elapsed_ns =
                PreferredTimingElapsedNanoseconds(begin, end, saturated);
            const int64_t minflt_delta =
                usage_before.MinorFaults < 0 || usage_after.MinorFaults < 0 ?
                -1 : usage_after.MinorFaults - usage_before.MinorFaults;
            const int64_t majflt_delta =
                usage_before.MajorFaults < 0 || usage_after.MajorFaults < 0 ?
                -1 : usage_after.MajorFaults - usage_before.MajorFaults;
            const uint32_t binary_def =
                stats.InactivatedColumns >= stats.BinaryResidualRank ?
                stats.InactivatedColumns - stats.BinaryResidualRank : 0u;
            const uint32_t heavy_gain =
                stats.ResidualRank >= stats.BinaryResidualRank ?
                stats.ResidualRank - stats.BinaryResidualRank : 0u;
            const bool outcome_stable = result == preflight_result;
            std::ostringstream row_arm_receipt;
            row_arm_receipt << MixedCoefficientGeometryName(arm.Geometry)
                << ',' << arm.GF256Rows << ',' << arm.GF16Rows
                << ',' << arm.ResidueSkew
                << ',' << MixedResidueScheduleName(arm.ResidueSchedule)
                << ",0x" << std::hex << arm.ResidueHashSeed << std::dec
                << ',' << (arm.ResidueSkew != 0u ||
                    arm.ResidueSchedule !=
                        wirehair_v2::MixedResidueSchedule::Constant ? 1u : 0u)
                << ',' << (arm.IndependentExtensionResidues ? 1u : 0u)
                << ",0x" << std::hex << arm.IndependentExtensionSeedXor
                << ",0x" << active_extension_hash_seed
                << ",0x" << active_grouped_hash_seed << std::dec
                << ',' << arm.PacketSeedMultiplier
                << ',' << (arm.PacketSeedAvalanche ? 1u : 0u)
                << ",0x" << std::hex << arm.OddPacketPeelSeedXor << std::dec
                << ',' << arm.DenseRowsOverride
                << ',' << (arm.DenseTwoAnchor ? 1u : 0u)
                << ',' << arm.DenseTwoAnchorPhase
                << ',' << active_system.Params.Staircase
                << ',' << active_system.Params.DenseRows
                << ',' << active_system.Params.HeavyRows
                << ',' << active_system.Params.SourceHits
                << ',' << static_cast<uint32_t>(active_system.Params.Field)
                << ',' << (active_system.Params.DenseIdentityCorner ? 1u : 0u)
                << ',' <<
                    static_cast<uint32_t>(active_system.Params.HeavyFamily)
                << ',' << active_config.MixCount
                << ',' << GroupedTimingRhsRouteName(expected_route)
                << ',' << GroupedTimingRhsRouteName(stats.MixedRhsRoute);
            std::printf(
                "%u,%u,%u,%s,%llu,%.17g,%s,%u,%u,%s,%u,%u,%s,%u,"
                "0x%llx,0x%x,%d,%s,%u,%d,%u,%llu,%u,%d,%d,%d,%lld,%lld,"
                "%d,%u,%u,%u,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,"
                "%llu,%llu,%u,%llu,%llu,%llu,%llu,%zu,%llu,%llu,%llu,%s\n",
                K, block_bytes, overhead, PacketScheduleName(schedule),
                (unsigned long long)seed, loss,
                GroupedTimingCacheStateName(cache_state), cycle, slot,
                control ? "control" : "candidate", arm.Period,
                arm.GroupedRows, MixedResidueBucketModeName(arm.Buckets),
                active_attempt,
                (unsigned long long)active_system.Params.Seed,
                active_config.PeelSeed,
                (int)preflight_result, outcome_class,
                common_success ? 1u : 0u, (int)result,
                outcome_stable ? 1u : 0u,
                (unsigned long long)elapsed_ns, saturated ? 1u : 0u,
                cpu_before, cpu_after,
                GroupedTimingCpuMigrated(cpu_before, cpu_after),
                (long long)minflt_delta, (long long)majflt_delta,
                GroupedTimingFaultContaminated(
                    minflt_delta, majflt_delta),
                stats.InactivatedColumns, binary_def, heavy_gain,
                (unsigned long long)stats.BlockXors,
                (unsigned long long)stats.BlockMulAdds,
                (unsigned long long)stats.BuildNanoseconds,
                (unsigned long long)stats.PeelNanoseconds,
                (unsigned long long)stats.ProjectNanoseconds,
                (unsigned long long)stats.ResidualNanoseconds,
                (unsigned long long)stats.BackSubNanoseconds,
                (unsigned long long)stats.MixedJointSourceXors,
                (unsigned long long)stats.MixedJointMarginalXors,
                (unsigned long long)stats.MixedJointMarginalCopies,
                stats.MixedJointActiveDeltas,
                (unsigned long long)stats.MixedJointScratchBytes,
                (unsigned long long)stats.MixedDualSourceColumns,
                (unsigned long long)((uint64_t)K * block_bytes),
                (unsigned long long)packet_payload_bytes,
                intermediate.size(),
                (unsigned long long)stats.SolveValueArenaBytes,
                (unsigned long long)stats.SolveValueArenaEagerZeroBytes,
                (unsigned long long)stats.SolveValueArenaCommitCopyBytes,
                row_arm_receipt.str().c_str());
            const uint64_t expected_arena_bytes =
                ((uint64_t)K + active_system.Params.Staircase +
                    active_system.Params.DenseRows +
                    active_system.Params.HeavyRows) * block_bytes;
            if (saturated || !outcome_stable ||
                stats.SolveValueArenaBytes != expected_arena_bytes ||
                stats.SolveValueArenaEagerZeroBytes != 0u ||
                stats.SolveValueArenaCommitCopyBytes != 0u ||
                stats.MixedRhsRoute != preflight_route ||
                (result == Wirehair_Success &&
                 stats.MixedRhsRoute != expected_route) ||
                (result != Wirehair_Success && result != Wirehair_NeedMore))
            {
                std::fprintf(stderr,
                    "groupedtiming run failed cycle=%u slot=%u result=%d "
                    "preflight=%d\n",
                    cycle, slot, (int)result, (int)preflight_result);
                return 2;
            }
        }
    }
    return 0;
}

#endif // test hooks && !WIREHAIR_V2_BENCH_DISABLE_PREFERRED_ATTEMPT

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
    PacketScheduleKind schedule_kind = PacketScheduleKind::Iid;
    bool mixed_residue_hash_keyed = false;
    bool mixed_independent_extension_residues = false;
    uint32_t mixed_extension_residue_seed_xor = 78u;
    uint32_t source_hits_override = 0u;
    uint32_t binary_dense_rows_override = 0u;
    uint32_t gf256_heavy_rows_override = 0u;
    uint32_t packet_peel_seed_xor = 0u;
    PacketPeelSeedTable packet_peel_seed_table =
        PacketPeelSeedTable::None;
    uint32_t odd_packet_peel_seed_xor = 0u;
    uint32_t packet_row_seed_multiplier = 1u;
    bool packet_row_seed_avalanche = false;
    uint32_t seed_block_bytes_override = 0u;
    bool paired_overhead_stream = false;
    bool mixed_null_witnesses = false;
    bool binary_dense_two_anchor = false;
    uint32_t binary_dense_two_anchor_phase = 0u;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    bool mixed_null_witness_internal_error = false;
    uint32_t fail_thread_launch_after = UINT32_MAX;
    bool source_hits_explicit = false;
    bool binary_dense_rows_explicit = false;
    bool binary_dense_two_anchor_phase_explicit = false;
    bool gf256_heavy_rows_explicit = false;
    bool packet_peel_seed_xor_explicit = false;
    uint32_t mixed_period = wirehair_v2::kMixedCoefficientPeriod;
    bool mixed_period_explicit = false;
    uint32_t mixed_gf256_rows = wirehair_v2::kMixedGF256Rows;
    bool mixed_gf256_rows_explicit = false;
    uint32_t mixed_grouped_gf256_rows = 0u;
    bool mixed_grouped_gf256_rows_explicit = false;
    uint32_t mixed_gf16_rows = wirehair_v2::kMixedGF16Rows;
    bool mixed_gf16_rows_explicit = false;
    wirehair_v2::MixedCoefficientGeometry mixed_geometry =
        wirehair_v2::MixedCoefficientGeometry::FrozenPowerX;
    bool mixed_geometry_explicit = false;
    uint32_t mixed_residue_skew = 0u;
    bool mixed_residue_skew_explicit = false;
    wirehair_v2::MixedResidueSchedule mixed_residue_schedule =
        wirehair_v2::MixedResidueSchedule::Constant;
    bool mixed_residue_schedule_explicit = false;
    uint32_t mixed_residue_hash_seed = 0u;
    bool mixed_residue_hash_seed_explicit = false;
    bool mixed_extension_residue_seed_xor_explicit = false;
    wirehair_v2::MixedResidueBucketMode mixed_residue_bucket_mode =
        wirehair_v2::MixedResidueBucketMode::Automatic;
    bool mixed_residue_bucket_mode_explicit = false;
    bool seed_block_bytes_explicit = false;
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
        else if (!std::strcmp(argv[i], "--schedule")) {
            if (!TakeArg(
                    "precodefail", "--schedule", argc, argv, i, value) ||
                !ParsePacketSchedule(value, schedule_kind))
            {
                std::fprintf(stderr,
                    "unknown precodefail schedule %s "
                    "(expected iid, burst, permutation, systematic-first, "
                    "repair-only, or adversarial)\n",
                    value ? value : "");
                return 1;
            }
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        else if (!std::strcmp(argv[i], "--full-payload-solve")) {
            full_payload_solve = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-null-witnesses")) {
            mixed_null_witnesses = true;
        }
        else if (!std::strcmp(argv[i], "--source-hits")) {
            if (!TakeArg(
                    "precodefail", "--source-hits", argc, argv, i, value) ||
                !ParseU32Arg(
                    "--source-hits", value, source_hits_override))
            {
                return 1;
            }
            source_hits_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--binary-dense-rows")) {
            if (!TakeArg(
                    "precodefail", "--binary-dense-rows",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--binary-dense-rows", value,
                    binary_dense_rows_override))
            {
                return 1;
            }
            binary_dense_rows_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--binary-dense-two-anchor")) {
            binary_dense_two_anchor = true;
        }
        else if (!std::strcmp(
                     argv[i], "--binary-dense-two-anchor-phase"))
        {
            if (!TakeArg(
                    "precodefail", "--binary-dense-two-anchor-phase",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--binary-dense-two-anchor-phase", value,
                    binary_dense_two_anchor_phase))
            {
                return 1;
            }
            binary_dense_two_anchor_phase_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--gf256-heavy-rows")) {
            if (!TakeArg(
                    "precodefail", "--gf256-heavy-rows",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--gf256-heavy-rows", value,
                    gf256_heavy_rows_override))
            {
                return 1;
            }
            gf256_heavy_rows_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--packet-peel-seed-xor")) {
            if (!TakeArg(
                    "precodefail", "--packet-peel-seed-xor",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--packet-peel-seed-xor", value,
                    packet_peel_seed_xor))
            {
                return 1;
            }
            packet_peel_seed_xor_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--packet-peel-seed-table")) {
            if (!TakeArg(
                    "precodefail", "--packet-peel-seed-table",
                    argc, argv, i, value))
            {
                return 1;
            }
            if (std::strcmp(value, "normalized-h15-v1") == 0) {
                packet_peel_seed_table =
                    PacketPeelSeedTable::NormalizedH15V1;
            }
            else if (std::strcmp(value, "normalized-h15-v2") == 0) {
                packet_peel_seed_table =
                    PacketPeelSeedTable::NormalizedH15V2;
            }
            else if (std::strcmp(value, "normalized-h15-v3") == 0) {
                packet_peel_seed_table =
                    PacketPeelSeedTable::NormalizedH15V3;
            }
            else if (std::strcmp(value, "normalized-h15-v4") == 0) {
                packet_peel_seed_table =
                    PacketPeelSeedTable::NormalizedH15V4;
            }
            else
            {
                std::fprintf(stderr,
                    "precodefail unknown --packet-peel-seed-table %s "
                    "(expected normalized-h15-v1, normalized-h15-v2, "
                    "normalized-h15-v3, or normalized-h15-v4)\n",
                    value);
                return 1;
            }
        }
        else if (!std::strcmp(
                     argv[i], "--odd-packet-peel-seed-xor"))
        {
            if (!TakeArg(
                    "precodefail", "--odd-packet-peel-seed-xor",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--odd-packet-peel-seed-xor", value,
                    odd_packet_peel_seed_xor))
            {
                return 1;
            }
        }
        else if (!std::strcmp(
                     argv[i], "--packet-row-seed-multiplier"))
        {
            if (!TakeArg(
                    "precodefail", "--packet-row-seed-multiplier",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--packet-row-seed-multiplier", value,
                    packet_row_seed_multiplier))
            {
                return 1;
            }
        }
        else if (!std::strcmp(
                     argv[i], "--packet-row-seed-avalanche"))
        {
            packet_row_seed_avalanche = true;
        }
        else if (!std::strcmp(argv[i], "--paired-overhead-stream")) {
            paired_overhead_stream = true;
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
        else if (!std::strcmp(argv[i], "--mixed-gf256-rows")) {
            if (!TakeArg(
                    "precodefail", "--mixed-gf256-rows",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-gf256-rows", value, mixed_gf256_rows))
            {
                return 1;
            }
            mixed_gf256_rows_explicit = true;
        }
        else if (!std::strcmp(
                     argv[i], "--mixed-grouped-gf256-rows"))
        {
            if (!TakeArg(
                    "precodefail", "--mixed-grouped-gf256-rows",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-grouped-gf256-rows", value,
                    mixed_grouped_gf256_rows))
            {
                return 1;
            }
            mixed_grouped_gf256_rows_explicit = true;
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
        else if (!std::strcmp(argv[i], "--mixed-residue-skew")) {
            if (!TakeArg(
                    "precodefail", "--mixed-residue-skew",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-residue-skew", value, mixed_residue_skew))
            {
                return 1;
            }
            mixed_residue_skew_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-schedule")) {
            if (!TakeArg(
                    "precodefail", "--mixed-residue-schedule",
                    argc, argv, i, value) ||
                !ParseMixedResidueSchedule(value, mixed_residue_schedule))
            {
                std::fprintf(stderr,
                    "precodefail unknown --mixed-residue-schedule token %s "
                    "(expected constant, ramp, or hashed)\n",
                    value ? value : "");
                return 1;
            }
            mixed_residue_schedule_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-hash-seed")) {
            if (!TakeArg(
                    "precodefail", "--mixed-residue-hash-seed",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-residue-hash-seed", value,
                    mixed_residue_hash_seed))
            {
                return 1;
            }
            mixed_residue_hash_seed_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-hash-keyed")) {
            mixed_residue_hash_keyed = true;
        }
        else if (!std::strcmp(
                     argv[i],
                     "--mixed-independent-extension-residues"))
        {
            mixed_independent_extension_residues = true;
        }
        else if (!std::strcmp(
                     argv[i], "--mixed-extension-residue-seed-xor"))
        {
            if (!TakeArg(
                    "precodefail", "--mixed-extension-residue-seed-xor",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--mixed-extension-residue-seed-xor", value,
                    mixed_extension_residue_seed_xor))
            {
                return 1;
            }
            mixed_extension_residue_seed_xor_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--mixed-residue-buckets"))
        {
            if (!TakeArg(
                    "precodefail", "--mixed-residue-buckets",
                    argc, argv, i, value) ||
                !ParseMixedResidueBucketMode(
                    value, mixed_residue_bucket_mode))
            {
                std::fprintf(stderr,
                    "precodefail unknown --mixed-residue-buckets token %s "
                    "(expected auto, separate, dual, or joint-delta)\n",
                    value ? value : "");
                return 1;
            }
            mixed_residue_bucket_mode_explicit = true;
        }
        else if (!std::strcmp(argv[i], "--seed-block-bytes")) {
            if (!TakeArg(
                    "precodefail", "--seed-block-bytes",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--seed-block-bytes", value,
                    seed_block_bytes_override))
            {
                return 1;
            }
            seed_block_bytes_explicit = true;
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
    else if (mixed_null_witnesses || mixed_period_explicit ||
             mixed_geometry_explicit ||
             mixed_gf256_rows_explicit ||
             mixed_grouped_gf256_rows_explicit ||
             mixed_gf16_rows_explicit ||
             mixed_residue_skew_explicit ||
             mixed_residue_schedule_explicit ||
             mixed_residue_hash_seed_explicit ||
             mixed_residue_hash_keyed ||
             mixed_independent_extension_residues ||
             mixed_residue_bucket_mode_explicit)
    {
        std::fprintf(stderr,
            "precodefail mixed experiment flags require --completion "
            "mixed\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedCoefficientGeometryForTesting(mixed_geometry)) {
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
    if (source_hits_explicit &&
        (source_hits_override == 0u || source_hits_override > 8u))
    {
        std::fprintf(stderr,
            "precodefail --source-hits must be in [1,8]\n");
        return 1;
    }
    if (binary_dense_rows_explicit &&
        (binary_dense_rows_override == 0u ||
         binary_dense_rows_override > 64u))
    {
        std::fprintf(stderr,
            "precodefail --binary-dense-rows must be in [1,64]\n");
        return 1;
    }
    if (binary_dense_two_anchor && binary_dense_rows_override != 0u &&
        binary_dense_rows_override != 12u)
    {
        std::fprintf(stderr,
            "precodefail --binary-dense-two-anchor requires 12 binary "
            "dense rows\n");
        return 1;
    }
    if (binary_dense_two_anchor_phase_explicit &&
        !binary_dense_two_anchor)
    {
        std::fprintf(stderr,
            "precodefail --binary-dense-two-anchor-phase requires "
            "--binary-dense-two-anchor\n");
        return 1;
    }
    if (binary_dense_two_anchor_phase > 2u)
    {
        std::fprintf(stderr,
            "precodefail --binary-dense-two-anchor-phase must be in "
            "[0,2]\n");
        return 1;
    }
    if (gf256_heavy_rows_explicit &&
        (gf256_heavy_rows_override == 0u ||
         gf256_heavy_rows_override > 128u))
    {
        std::fprintf(stderr,
            "precodefail --gf256-heavy-rows must be in [1,128]\n");
        return 1;
    }
    if (gf256_heavy_rows_override != 0u &&
        completion != PrecodeFailCompletion::Certified)
    {
        std::fprintf(stderr,
            "precodefail --gf256-heavy-rows requires --completion "
            "certified\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedCoefficientPeriodForTesting(mixed_period))
    {
        std::fprintf(stderr,
            "precodefail --mixed-period must be in [%u,%u]\n",
            wirehair_v2::ActiveMixedGF256Rows() +
                wirehair_v2::ActiveMixedGF16Rows(),
            wirehair_v2::kMixedCoefficientPeriod);
        return 1;
    }
    if (!wirehair_v2::SetMixedGF256RowsForTesting(mixed_gf256_rows))
    {
        std::fprintf(stderr,
            "precodefail --mixed-gf256-rows must be in [%u,%u], fit the "
            "active period, use shared-x for an extra row, and use the "
            "validated 12+4 geometry for twelve rows\n",
            wirehair_v2::kMixedGF256Rows,
            wirehair_v2::kMixedGF256RowsMax);
        return 1;
    }
    if (!wirehair_v2::SetMixedResidueSkewForTesting(mixed_residue_skew)) {
        std::fprintf(stderr,
            "precodefail --mixed-residue-skew must be a corner-preserving "
            "shared-x skew in [0,P-H]\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedResidueScheduleForTesting(
            mixed_residue_schedule))
    {
        std::fprintf(stderr,
            "precodefail nonconstant --mixed-residue-schedule requires "
            "shared-x, "
            "P>H, and zero constant skew\n");
        return 1;
    }
    if ((mixed_residue_hash_seed_explicit || mixed_residue_hash_keyed) &&
        mixed_residue_schedule != wirehair_v2::MixedResidueSchedule::Hashed)
    {
        std::fprintf(stderr,
            "precodefail residue hash seed/keying requires hashed "
            "--mixed-residue-schedule\n");
        return 1;
    }
    if (mixed_extension_residue_seed_xor_explicit &&
        !mixed_independent_extension_residues)
    {
        std::fprintf(stderr,
            "precodefail --mixed-extension-residue-seed-xor requires "
            "--mixed-independent-extension-residues\n");
        return 1;
    }
    wirehair_v2::SetMixedResidueHashSeedForTesting(
        mixed_residue_hash_seed);
    wirehair_v2::SetMixedIndependentExtensionSeedXorForTesting(
        mixed_extension_residue_seed_xor);
    if (!wirehair_v2::SetMixedIndependentExtensionResiduesForTesting(
            mixed_independent_extension_residues))
    {
        std::fprintf(stderr,
            "precodefail independent extension residues require "
            "shared-x hashed scheduling with P>H\n");
        return 1;
    }
    if (mixed_residue_bucket_mode !=
            wirehair_v2::MixedResidueBucketMode::Automatic &&
        !mixed_independent_extension_residues &&
        mixed_grouped_gf256_rows == 0u)
    {
        std::fprintf(stderr,
            "precodefail explicit --mixed-residue-buckets requires "
            "--mixed-independent-extension-residues or nonzero "
            "--mixed-grouped-gf256-rows\n");
        return 1;
    }
    if (!wirehair_v2::SetMixedResidueBucketModeForTesting(
            mixed_residue_bucket_mode))
    {
        return 1;
    }
    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(
            packet_row_seed_multiplier))
    {
        std::fprintf(stderr,
            "precodefail --packet-row-seed-multiplier must be odd and "
            "nonzero\n");
        return 1;
    }
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(
        packet_row_seed_avalanche);
    wirehair_v2::SetOddPacketPeelSeedXorForTesting(
        odd_packet_peel_seed_xor);
    // All prerequisite setters above intentionally clear schedule experiments.
    // Grouped C must therefore be the final thread-local codec configuration.
    if (!wirehair_v2::SetMixedGroupedGF256RowsForTesting(
            mixed_grouped_gf256_rows))
    {
        std::fprintf(stderr,
            "precodefail --mixed-grouped-gf256-rows must be in [0,9]; "
            "nonzero grouping requires shared-x constant-A 10+2 geometry, "
            "P>H, and no independent extension residues\n");
        return 1;
    }
    if (seed_block_bytes_explicit && seed_block_bytes_override == 0u)
    {
        std::fprintf(stderr,
            "precodefail --seed-block-bytes must be nonzero\n");
        return 1;
    }
    if (packet_peel_seed_table != PacketPeelSeedTable::None &&
        packet_peel_seed_xor_explicit)
    {
        std::fprintf(stderr,
            "precodefail --packet-peel-seed-table conflicts with "
            "--packet-peel-seed-xor\n");
        return 1;
    }
    if (packet_peel_seed_table != PacketPeelSeedTable::None &&
        (completion != PrecodeFailCompletion::Mixed ||
         mix_counts.size() != 1u || mix_counts[0] != 2 ||
         seed_block_bytes_override != 1280u ||
         mixed_gf256_rows != 11u || mixed_gf16_rows != 4u ||
         mixed_period != 32u ||
         mixed_geometry !=
            wirehair_v2::MixedCoefficientGeometry::SharedCauchyX ||
         mixed_residue_schedule !=
            wirehair_v2::MixedResidueSchedule::Hashed ||
         mixed_residue_hash_seed != 68u ||
         !mixed_residue_hash_keyed ||
         !mixed_independent_extension_residues ||
         mixed_extension_residue_seed_xor != 78u ||
         source_hits_override != 0u || binary_dense_rows_override != 0u ||
         odd_packet_peel_seed_xor != 0u ||
         packet_row_seed_multiplier != 1u ||
         packet_row_seed_avalanche))
    {
        std::fprintf(stderr,
            "precodefail normalized H15 packet seed table requires "
            "its normalized H15/mix2 geometry\n");
        return 1;
    }
#endif

    if (completion == PrecodeFailCompletion::Certified)
    {
        std::printf(
            "# precodefail: trials=%u threads=%u loss=%.17g seed=0x%llx "
            "source_hits_override=%u packet_peel_seed_xor=0x%x "
            "packet_peel_seed_table=%s "
            "binary_dense_rows_override=%u binary_dense_two_anchor=%u "
            "binary_dense_two_anchor_phase=%u "
            "gf256_heavy_rows_override=%u "
            "odd_packet_peel_seed_xor=0x%x "
            "packet_row_seed_multiplier=0x%x "
            "packet_row_seed_avalanche=%u seed_block_bytes_override=%u "
            "overhead_stream=%s full_payload_solve=%u schedule=%s%s\n",
            trials, threads, loss, (unsigned long long)seed,
            source_hits_override,
            packet_peel_seed_xor,
            PacketPeelSeedTableName(packet_peel_seed_table),
            binary_dense_rows_override,
            binary_dense_two_anchor ? 1u : 0u,
            binary_dense_two_anchor_phase,
            gf256_heavy_rows_override,
            odd_packet_peel_seed_xor,
            packet_row_seed_multiplier,
            packet_row_seed_avalanche ? 1u : 0u,
            seed_block_bytes_override,
            paired_overhead_stream ? "paired" : "salted",
            full_payload_solve ? 1u : 0u,
            PacketScheduleName(schedule_kind),
            mixed_null_witnesses ? " mixed_null_witnesses=1" : "");
    }
    else
    {
        std::printf(
            "# precodefail: trials=%u threads=%u loss=%.17g seed=0x%llx "
            "completion=%s mixed_period=%u mixed_gf256_rows=%u "
            "mixed_gf16_rows=%u "
            "mixed_geometry=%s mixed_residue_skew=%u "
            "mixed_residue_schedule=%s mixed_residue_hash_seed=0x%x "
            "mixed_residue_hash_keyed=%u "
            "mixed_independent_extension_residues=%u "
            "mixed_grouped_gf256_rows=%u "
            "mixed_grouped_gf256_hash_seed=0x%x "
            "mixed_grouped_final_h_a_columns=%u "
            "mixed_residue_buckets_requested=%s "
            "mixed_extension_residue_seed_xor=0x%x "
            "source_hits_override=%u packet_peel_seed_xor=0x%x "
            "packet_peel_seed_table=%s "
            "binary_dense_rows_override=%u binary_dense_two_anchor=%u "
            "binary_dense_two_anchor_phase=%u "
            "gf256_heavy_rows_override=%u "
            "odd_packet_peel_seed_xor=0x%x "
            "packet_row_seed_multiplier=0x%x "
            "packet_row_seed_avalanche=%u seed_block_bytes_override=%u "
            "overhead_stream=%s full_payload_solve=%u schedule=%s%s\n",
            trials, threads, loss, (unsigned long long)seed,
            PrecodeFailCompletionName(completion),
            wirehair_v2::ActiveMixedCoefficientPeriod(),
            wirehair_v2::ActiveMixedGF256Rows(),
            wirehair_v2::ActiveMixedGF16Rows(),
            MixedCoefficientGeometryName(
                wirehair_v2::ActiveMixedCoefficientGeometry()),
            wirehair_v2::ActiveMixedResidueSkew(),
            MixedResidueScheduleName(
                wirehair_v2::ActiveMixedResidueSchedule()),
            wirehair_v2::ActiveMixedResidueHashSeed(),
            mixed_residue_hash_keyed ? 1u : 0u,
            mixed_independent_extension_residues ? 1u : 0u,
            wirehair_v2::ActiveMixedGroupedGF256Rows(),
            wirehair_v2::ActiveMixedGroupedGF256HashSeed(),
            wirehair_v2::ActiveMixedGroupedGF256Rows() != 0u ?
                wirehair_v2::ActiveMixedGF256Rows() +
                    wirehair_v2::ActiveMixedGF16Rows() : 0u,
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            MixedResidueBucketModeName(mixed_residue_bucket_mode),
#else
            "auto",
#endif
            mixed_extension_residue_seed_xor,
            source_hits_override,
            packet_peel_seed_xor,
            PacketPeelSeedTableName(packet_peel_seed_table),
            binary_dense_rows_override,
            binary_dense_two_anchor ? 1u : 0u,
            binary_dense_two_anchor_phase,
            gf256_heavy_rows_override,
            odd_packet_peel_seed_xor,
            packet_row_seed_multiplier,
            packet_row_seed_avalanche ? 1u : 0u,
            seed_block_bytes_override,
            paired_overhead_stream ? "paired" : "salted",
            full_payload_solve ? 1u : 0u,
            PacketScheduleName(schedule_kind),
            mixed_null_witnesses ? " mixed_null_witnesses=1" : "");
    }
    std::printf(
        "N,bb,heavy_family,mix_count,overhead,trials,success,rank_fail,error,"
        "fail_rate,"
        "inact_mu,inact_max,binary_def_mu,binary_def_max,heavy_gain_mu,"
        "heavy_gain_min,heavy_shortfall,solve_ms_mu,build_ms_mu,peel_ms_mu,"
        "project_ms_mu,residual_ms_mu,backsub_ms_mu,seed_attempt,"
        "block_xors_mu,block_muladds_mu,first_rank_fail,binary_def_hist,"
        "heavy_gain_hist,failure_trials,active_packet_peel_seed_xor,"
        "mixed_joint_source_xors_mu,mixed_joint_marginal_xors_mu,"
        "mixed_joint_marginal_copies_mu,mixed_joint_active_deltas_mu,"
        "mixed_joint_scratch_bytes_mu,mixed_dual_source_columns_mu\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        const uint32_t K = (uint32_t)n_value;
        const uint32_t bb = (uint32_t)bb_value;
        uint32_t active_packet_peel_seed_xor = packet_peel_seed_xor;
        if (packet_peel_seed_table ==
            PacketPeelSeedTable::NormalizedH15V1)
        {
            active_packet_peel_seed_xor =
                NormalizedH15V1PacketPeelSeedXor(K);
        }
        else if (packet_peel_seed_table ==
                 PacketPeelSeedTable::NormalizedH15V2)
        {
            active_packet_peel_seed_xor =
                NormalizedH15V2PacketPeelSeedXor(K);
        }
        else if (packet_peel_seed_table ==
                 PacketPeelSeedTable::NormalizedH15V3)
        {
            active_packet_peel_seed_xor =
                NormalizedH15V3PacketPeelSeedXor(K);
        }
        else if (packet_peel_seed_table ==
                 PacketPeelSeedTable::NormalizedH15V4)
        {
            active_packet_peel_seed_xor =
                NormalizedH15V4PacketPeelSeedXor(K);
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        uint32_t active_hash_seed = mixed_residue_hash_seed;
        if (mixed_residue_hash_keyed &&
            !wirehair_v2::SelectFullCycleMixedResidueKeyedSeedForTesting(
                mixed_residue_hash_seed, K, active_hash_seed))
        {
            std::fprintf(stderr,
                "precodefail could not select a full-cycle keyed residue "
                "hash seed for N=%u\n",
                K);
            return 1;
        }
        const auto configure_test_thread = [&]() -> bool {
            return wirehair_v2::SetMixedCoefficientGeometryForTesting(
                       mixed_geometry) &&
                wirehair_v2::SetMixedGF16RowsForTesting(mixed_gf16_rows) &&
                wirehair_v2::SetMixedCoefficientPeriodForTesting(
                    mixed_period) &&
                wirehair_v2::SetMixedGF256RowsForTesting(mixed_gf256_rows) &&
                wirehair_v2::SetMixedResidueSkewForTesting(
                    mixed_residue_skew) &&
                wirehair_v2::SetMixedResidueScheduleForTesting(
                    mixed_residue_schedule) &&
                (wirehair_v2::SetMixedResidueHashSeedForTesting(
                     active_hash_seed), true) &&
                (wirehair_v2::
                     SetMixedIndependentExtensionSeedXorForTesting(
                         mixed_extension_residue_seed_xor), true) &&
                wirehair_v2::
                    SetMixedIndependentExtensionResiduesForTesting(
                        mixed_independent_extension_residues) &&
                wirehair_v2::SetMixedResidueBucketModeForTesting(
                    mixed_residue_bucket_mode) &&
                wirehair_v2::SetPacketRowSeedMultiplierForTesting(
                    packet_row_seed_multiplier) &&
                (wirehair_v2::SetPacketRowSeedAvalancheForTesting(
                     packet_row_seed_avalanche), true) &&
                (wirehair_v2::SetOddPacketPeelSeedXorForTesting(
                     odd_packet_peel_seed_xor), true) &&
                wirehair_v2::SetMixedGroupedGF256RowsForTesting(
                    mixed_grouped_gf256_rows);
        };
        if (!configure_test_thread())
        {
            std::fprintf(stderr,
                "precodefail could not configure test thread for N=%u\n",
                K);
            return 1;
        }
#endif
        const wirehair_v2::SeedProfile profile =
            wirehair_v2::SelectSeedProfile(K, bb);
        // Research graph portability without changing the payload width used
        // by E2E validation, loss-stream pairing, or full-payload solving.
        const wirehair_v2::SeedProfile seed_profile =
            seed_block_bytes_override != 0u ?
                wirehair_v2::SelectSeedProfile(
                    K, seed_block_bytes_override) :
                profile;
        const uint64_t matrix_seed = wirehair_v2::MatrixSeedFromProfile(
            seed_profile, 0u, wirehair_v2::kMessagePrecodeSeedSalt);
        const wirehair_v2::PrecodeParams canonical_params =
            completion == PrecodeFailCompletion::Mixed ?
                wirehair_v2::MakeMixedParams(K, matrix_seed) :
                wirehair_v2::MakeCertifiedParams(K, matrix_seed);
        wirehair_v2::PacketRowConfig base_config;
        base_config.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
            seed_profile, wirehair_v2::kMessageRecoveryRowSeedSalt) ^
            active_packet_peel_seed_xor;
        for (wirehair_v2::HeavyCoefficientFamily heavy_family : heavy_families)
        {
        std::map<uint32_t, PairedMixOutcomes> paired_outcomes;
        for (int mix_count_value : mix_counts)
        {
            base_config.MixCount = (uint32_t)mix_count_value;
            wirehair_v2::PrecodeParams base_params = canonical_params;
            if (source_hits_override != 0u) {
                base_params.SourceHits = source_hits_override;
            }
            if (binary_dense_rows_override != 0u) {
                base_params.DenseRows = binary_dense_rows_override;
            }
            base_params.DenseTwoAnchor = binary_dense_two_anchor;
            base_params.DenseTwoAnchorPhase =
                binary_dense_two_anchor_phase;
            if (gf256_heavy_rows_override != 0u) {
                base_params.HeavyRows = gf256_heavy_rows_override;
            }
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
            const uint64_t precode_count_wide =
                (uint64_t)system.Params.Staircase +
                system.Params.DenseRows + system.Params.HeavyRows;
            if (precode_count_wide > UINT32_MAX) {
                std::fprintf(stderr,
                    "precodefail precode count overflow N=%u bb=%u\n",
                    K, bb);
                return 2;
            }
            const uint32_t precode_count = (uint32_t)precode_count_wide;
            wirehair_v2::PacketRowRuntime runtime;
            if (!runtime.Initialize(K, precode_count, config.MixCount)) {
                std::fprintf(stderr,
                    "precodefail packet runtime initialization failed "
                    "N=%u bb=%u heavy_family=%s mix_count=%u\n",
                    K, bb, HeavyFamilyName(heavy_family),
                    (uint32_t)mix_count_value);
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
            uint32_t solve_block_bytes =
                completion == PrecodeFailCompletion::Mixed ? 2u : 1u;
            if (full_payload_solve) solve_block_bytes = bb;
            const auto populate_trial_packets = [&, overhead](
                uint32_t trial,
                const uint8_t* packet_data,
                std::vector<wirehair_v2::SolvePacket>& packets) -> bool
            {
                const uint64_t loss_seed =
                    seed ^
                    ((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)) ^
                    ((uint64_t)bb * UINT64_C(0xbf58476d1ce4e5b9)) ^
                    ((uint64_t)(paired_overhead_stream ? 0u : overhead) *
                        UINT64_C(0x94d049bb133111eb)) ^
                    ((uint64_t)trial * UINT64_C(0xd6e8feb86659fd93));
                if (schedule_kind == PacketScheduleKind::Iid)
                {
                    Rng rng(loss_seed);
                    uint32_t block_id = 0u;
                    size_t delivered = 0u;
                    while (delivered < packets.size())
                    {
                        const uint32_t id = block_id++;
                        if (ShouldDrop(rng, loss)) continue;
                        wirehair_v2::SolvePacket& packet =
                            packets[delivered++];
                        packet.BlockId = id;
                        packet.Data = packet_data;
                    }
                    return true;
                }
                const std::vector<uint32_t> ids = BuildPacketSchedule(
                    K, (uint32_t)packets.size(), loss, loss_seed,
                    schedule_kind);
                if (ids.size() != packets.size()) return false;
                for (size_t i = 0u; i < ids.size(); ++i) {
                    packets[i].BlockId = ids[i];
                    packets[i].Data = packet_data;
                }
                return true;
            };
            std::vector<MatrixFailureTrial> results(trials);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            // Preserve the compact always-on benchmark record.  Full solver
            // metadata is only allocated and written for witness runs.
            std::vector<wirehair_v2::PrecodeSolveStats> witness_trial_stats;
            if (mixed_null_witnesses) witness_trial_stats.resize(trials);
            wirehair_v2::MixedNullWitnessDiagnostic captured_witness;
            uint32_t captured_witness_trial = UINT32_MAX;
            MixedNullReplayStatus witness_status =
                MixedNullReplayStatus::None;
            const char* witness_reason = "no_need_more";
            bool witness_replay_stats_ok = false;
            std::vector<std::string> witness_classification_lines;
#endif
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
                            if (!configure_test_thread()) {
                                throw std::runtime_error(
                                    "invalid worker test configuration");
                            }
#endif
                            const std::vector<uint8_t> zero(
                                solve_block_bytes, uint8_t{0});
                            std::vector<wirehair_v2::SolvePacket> packets(
                                (size_t)K + overhead);
                            for (;;)
                            {
                                if (cancel_workers.load()) {
                                    break;
                                }
                                const uint32_t trial = next_trial.fetch_add(1u);
                                if (trial >= trials) {
                                    break;
                                }
                                if (!populate_trial_packets(
                                        trial, zero.data(), packets))
                                {
                                    throw std::runtime_error(
                                        "packet schedule construction failed");
                                }
                                std::vector<uint8_t> intermediate;
                                wirehair_v2::PrecodeSolveStats solve_stats;
                                MatrixFailureTrial& result = results[trial];
                                result.Result =
                                    wirehair_v2::
                                        SolvePrecodeSystemForValidatedSystemWithRuntime(
                                        system, config, runtime, packets,
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
                                    TotalSolveNanoseconds(solve_stats);
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
                                result.MixedJointSourceXors =
                                    solve_stats.MixedJointSourceXors;
                                result.MixedJointMarginalXors =
                                    solve_stats.MixedJointMarginalXors;
                                result.MixedJointMarginalCopies =
                                    solve_stats.MixedJointMarginalCopies;
                                result.MixedJointScratchBytes =
                                    solve_stats.MixedJointScratchBytes;
                                result.MixedJointActiveDeltas =
                                    solve_stats.MixedJointActiveDeltas;
                                result.MixedDualSourceColumns =
                                    solve_stats.MixedDualSourceColumns;
                                if (mixed_null_witnesses) {
                                    witness_trial_stats[trial] = solve_stats;
                                }
#endif
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            if (mixed_residue_bucket_mode ==
                    wirehair_v2::MixedResidueBucketMode::JointDelta ||
                mixed_residue_bucket_mode ==
                    wirehair_v2::MixedResidueBucketMode::Dual)
            {
                const bool joint = mixed_residue_bucket_mode ==
                    wirehair_v2::MixedResidueBucketMode::JointDelta;
                for (uint32_t trial = 0u; trial < trials; ++trial)
                {
                    const MatrixFailureTrial& result = results[trial];
                    // Both q > H and a coefficient-rank-deficient q <= H
                    // return before any full mixed payload buckets are built.
                    const bool stopped_before_mixed_payload_buckets =
                        result.Result == Wirehair_NeedMore &&
                        result.Inactivated >= result.BinaryRank &&
                        result.Rank < result.Inactivated;
                    if (stopped_before_mixed_payload_buckets) {
                        continue;
                    }
                    const bool used = joint ?
                        result.MixedJointSourceXors != 0u &&
                            result.MixedJointMarginalCopies != 0u :
                        result.MixedDualSourceColumns != 0u;
                    if (used) continue;
                    std::fprintf(stderr,
                        "precodefail requested mixed residue bucket mode %s "
                        "but trial %u fell back (N=%u bb=%u result=%d "
                        "inactivated=%u binary_rank=%u heavy_rows=%u)\n",
                        MixedResidueBucketModeName(
                            mixed_residue_bucket_mode),
                        trial, K, bb, (int)result.Result,
                        result.Inactivated, result.BinaryRank,
                        system.Params.HeavyRows);
                    return 1;
                }
            }
            if (mixed_null_witnesses)
            {
                uint32_t first_need_more = UINT32_MAX;
                for (uint32_t trial = 0u; trial < trials; ++trial)
                {
                    const MatrixFailureTrial& candidate = results[trial];
                    if (candidate.Result != Wirehair_NeedMore) continue;
                    if (first_need_more == UINT32_MAX) first_need_more = trial;
                    const wirehair_v2::PrecodeSolveStats& stats =
                        witness_trial_stats[trial];
                    if (stats.InactivatedColumns <
                            stats.BinaryResidualRank ||
                        stats.ResidualRank < stats.BinaryResidualRank)
                    {
                        continue;
                    }
                    const uint32_t q = stats.InactivatedColumns -
                        stats.BinaryResidualRank;
                    const uint32_t quotient_rank = stats.ResidualRank -
                        stats.BinaryResidualRank;
                    if (q != 0u &&
                        q <= wirehair_v2::
                            kMaxMixedNullWitnessQuotientColumns &&
                        q <= system.Params.HeavyRows && quotient_rank < q)
                    {
                        captured_witness_trial = trial;
                        break;
                    }
                }
                if (captured_witness_trial == UINT32_MAX)
                {
                    if (first_need_more != UINT32_MAX) {
                        captured_witness_trial = first_need_more;
                        witness_status = MixedNullReplayStatus::Skipped;
                        witness_reason = "ineligible_need_more";
                    }
                }
                else
                {
                    try
                    {
                        std::vector<wirehair_v2::SolvePacket> packets(
                            (size_t)K + overhead);
                        const std::vector<uint8_t> replay_zero(
                            solve_block_bytes, uint8_t{0});
                        if (!configure_test_thread()) {
                            witness_status = MixedNullReplayStatus::Error;
                            witness_reason = "replay_config";
                        }
                        else if (!populate_trial_packets(
                                     captured_witness_trial,
                                     replay_zero.data(), packets))
                        {
                            witness_status = MixedNullReplayStatus::Error;
                            witness_reason = "replay_schedule";
                        }
                        else
                        {
                            std::vector<uint8_t> intermediate;
                            wirehair_v2::PrecodeSolveStats replay_stats;
                            WirehairResult replay_result;
                            {
                                MixedNullWitnessScope witness_scope(
                                    &captured_witness);
                                replay_result = wirehair_v2::
                                    SolvePrecodeSystemForValidatedSystemWithRuntime(
                                        system, config, runtime, packets,
                                        solve_block_bytes, intermediate,
                                        &replay_stats);
                            }
                            const MatrixFailureTrial& original =
                                results[captured_witness_trial];
                            witness_replay_stats_ok =
                                replay_result == original.Result &&
                                SameNonTimingSolveStats(
                                    replay_stats,
                                    witness_trial_stats[
                                        captured_witness_trial]);
                            if (!witness_replay_stats_ok) {
                                witness_status = MixedNullReplayStatus::Error;
                                witness_reason = "replay_mismatch";
                            }
                            else if (captured_witness.Status !=
                                     wirehair_v2::
                                         MixedNullWitnessStatus::Captured)
                            {
                                witness_status = MixedNullReplayStatus::Error;
                                witness_reason = "diagnostic_failed";
                            }
                            else {
                                witness_status =
                                    MixedNullReplayStatus::Captured;
                                witness_reason = "verified";
                                if (!BuildMixedNullClassification(
                                        K, captured_witness,
                                        witness_classification_lines))
                                {
                                    witness_classification_lines.clear();
                                    witness_status =
                                        MixedNullReplayStatus::Error;
                                    witness_reason = "classification_failed";
                                }
                            }
                        }
                    }
                    catch (...) {
                        witness_classification_lines.clear();
                        witness_status = MixedNullReplayStatus::Error;
                        witness_reason = "replay_exception";
                    }
                }
                if (witness_status == MixedNullReplayStatus::Error) {
                    mixed_null_witness_internal_error = true;
                }
            }
#endif

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
            uint64_t mixed_joint_source_xors_sum = 0u;
            uint64_t mixed_joint_marginal_xors_sum = 0u;
            uint64_t mixed_joint_marginal_copies_sum = 0u;
            uint64_t mixed_joint_active_deltas_sum = 0u;
            uint64_t mixed_joint_scratch_bytes_sum = 0u;
            uint64_t mixed_dual_source_columns_sum = 0u;
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
                mixed_joint_source_xors_sum += result.MixedJointSourceXors;
                mixed_joint_marginal_xors_sum +=
                    result.MixedJointMarginalXors;
                mixed_joint_marginal_copies_sum +=
                    result.MixedJointMarginalCopies;
                mixed_joint_active_deltas_sum +=
                    result.MixedJointActiveDeltas;
                mixed_joint_scratch_bytes_sum +=
                    result.MixedJointScratchBytes;
                mixed_dual_source_columns_sum +=
                    result.MixedDualSourceColumns;
                inact_max = std::max(
                    inact_max, result.Inactivated);
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
                "%s,%s,%s,0x%x,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
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
                failure_trials.c_str(), active_packet_peel_seed_xor,
                (double)mixed_joint_source_xors_sum / trials,
                (double)mixed_joint_marginal_xors_sum / trials,
                (double)mixed_joint_marginal_copies_sum / trials,
                (double)mixed_joint_active_deltas_sum / trials,
                (double)mixed_joint_scratch_bytes_sum / trials,
                (double)mixed_dual_source_columns_sum / trials);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            if (mixed_null_witnesses)
            {
                const wirehair_v2::MixedNullWitnessDiagnostic& witness =
                    captured_witness;
                const bool have_trial = captured_witness_trial < trials;
                const bool have_diagnostic =
                    witness_status == MixedNullReplayStatus::Captured;
                const wirehair_v2::PrecodeSolveStats* selected_stats =
                    have_trial ?
                        &witness_trial_stats[captured_witness_trial] : nullptr;
                const uint32_t selected_R = selected_stats ?
                    selected_stats->InactivatedColumns : 0u;
                const uint32_t selected_binary_rank = selected_stats ?
                    selected_stats->BinaryResidualRank : 0u;
                const uint32_t selected_rank = selected_stats ?
                    selected_stats->ResidualRank : 0u;
                const uint32_t selected_q =
                    selected_R >= selected_binary_rank ?
                    selected_R - selected_binary_rank : 0u;
                const uint32_t selected_qrank =
                    selected_rank >= selected_binary_rank ?
                    selected_rank - selected_binary_rank : 0u;
                const uint32_t L = have_diagnostic ?
                    witness.ColumnCount : K + precode_count;
                const uint32_t R = have_diagnostic ?
                    witness.InactiveCount : selected_R;
                const uint32_t binary_rank = have_diagnostic ?
                    witness.BinaryRank : selected_binary_rank;
                const uint32_t q = have_diagnostic ?
                    witness.QuotientColumns : selected_q;
                const uint32_t qrank = have_diagnostic ?
                    witness.QuotientRank : selected_qrank;
                const uint32_t d = have_diagnostic ?
                    witness.KernelDimension :
                    (q >= qrank ? q - qrank : 0u);
                const uint64_t expected_words = (uint64_t)d * L;
                const size_t exact_words = have_diagnostic ?
                    witness.CanonicalBasis.size() : 0u;
                const uint64_t hash_high = have_diagnostic ?
                    witness.BasisHashHigh : 0u;
                const uint64_t hash_low = have_diagnostic ?
                    witness.BasisHashLow : 0u;
                const uint32_t active_subfield_rows =
                    wirehair_v2::ActiveMixedGF256Rows();
                const uint32_t active_extension_rows =
                    wirehair_v2::ActiveMixedGF16Rows();
                const uint32_t active_grouped_rows =
                    wirehair_v2::ActiveMixedGroupedGF256Rows();
                const uint32_t active_first_grouped_row =
                    active_grouped_rows <= active_subfield_rows ?
                        active_subfield_rows - active_grouped_rows : 0u;
                std::printf(
                    "# mixed_null_witness,v=2,N=%u,bb=%u,trial=%d,"
                    "status=%s,reason=%s,diagnostic_status=%u,"
                    "replay_stats_ok=%u,period=%u,schedule=%s,"
                    "hash_seed=0x%x,extension_seed_xor=0x%x,"
                    "independent_extension=%u,L=%u,R=%u,binary_rank=%u,q=%u,"
                    "quotient_rank=%u,d=%u,q_ok=%u,A_ok=%u,C_ok=%u,"
                    "canonical_ok=%u,exact_words=%zu,exact_size_ok=%u,"
                    "hash=%016llx%016llx,grouped_gf256_rows=%u,"
                    "grouped_first_row=%u,grouped_hash_seed=0x%x,"
                    "grouped_final_h_a_columns=%u\n",
                    K, bb, have_trial ? (int)captured_witness_trial : -1,
                    MixedNullReplayStatusName(witness_status),
                    witness_reason, (uint32_t)witness.Status,
                    witness_replay_stats_ok ? 1u : 0u,
                    wirehair_v2::ActiveMixedCoefficientPeriod(),
                    MixedResidueScheduleName(
                        wirehair_v2::ActiveMixedResidueSchedule()),
                    active_hash_seed, mixed_extension_residue_seed_xor,
                    mixed_independent_extension_residues ? 1u : 0u,
                    L, R, binary_rank, q, qrank, d,
                    have_diagnostic && witness.QuotientVerified ? 1u : 0u,
                    have_diagnostic && witness.BinaryVerified ? 1u : 0u,
                    have_diagnostic && witness.CompletionVerified ? 1u : 0u,
                    have_diagnostic && witness.CanonicalVerified ? 1u : 0u,
                    exact_words,
                    have_diagnostic && expected_words == exact_words ?
                        1u : 0u,
                    (unsigned long long)hash_high,
                    (unsigned long long)hash_low,
                    active_grouped_rows, active_first_grouped_row,
                    wirehair_v2::ActiveMixedGroupedGF256HashSeed(),
                    active_grouped_rows != 0u ?
                        active_subfield_rows + active_extension_rows : 0u);
                if (witness_status == MixedNullReplayStatus::Captured) {
                    for (const std::string& line :
                         witness_classification_lines)
                    {
                        std::printf("%s\n", line.c_str());
                    }
                }
            }
#endif
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedNullReplayExitCode(
        mixed_null_witness_internal_error ? MixedNullReplayStatus::Error :
            MixedNullReplayStatus::None);
#else
    return 0;
#endif
}

int CmdSelfTest()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (kMixedNullReplayInternalErrorExitCode != 3 ||
        MixedNullReplayExitCode(MixedNullReplayStatus::None) != 0 ||
        MixedNullReplayExitCode(MixedNullReplayStatus::Captured) != 0 ||
        MixedNullReplayExitCode(MixedNullReplayStatus::Skipped) != 0 ||
        MixedNullReplayExitCode(MixedNullReplayStatus::Error) !=
            kMixedNullReplayInternalErrorExitCode)
    {
        std::fprintf(stderr, "mixed null-witness exit policy mismatch\n");
        return 1;
    }
    std::printf("mixed null-witness exit policy: PASS\n");
#endif

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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \
    !defined(WIREHAIR_V2_BENCH_DISABLE_PREFERRED_ATTEMPT)
        std::fprintf(stderr,
            "usage: wirehair_v2_bench compare|precodecheck|seedtable|"
            "peelcost|densecheck|densetune|densecount|densegrid|precodefail|"
            "preferredattempt|preferredtiming|groupedtiming|selftest "
            "[opts]\n");
#else
        std::fprintf(stderr,
            "usage: wirehair_v2_bench compare|precodecheck|seedtable|"
            "peelcost|densecheck|densetune|densecount|densegrid|precodefail|"
            "selftest [opts]\n");
#endif
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \
    !defined(WIREHAIR_V2_BENCH_DISABLE_PREFERRED_ATTEMPT)
        if (!std::strcmp(argv[1], "preferredattempt")) {
            return CmdPreferredAttempt(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "preferredtiming")) {
            return CmdPreferredTiming(argc - 2, argv + 2);
        }
        if (!std::strcmp(argv[1], "groupedtiming")) {
            return CmdGroupedTiming(argc - 2, argv + 2);
        }
#endif
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
