#include "WirehairV2Codec.h"
#include "WirehairV2Plan.h"
#include "WirehairV2Precode.h"
#include "WirehairV2Seeds.h"

#include "../WirehairTools.h"

#include <wirehair/wirehair.h>

#include <algorithm>
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

struct TrialResult
{
    bool Ok;
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
    uint64_t max_message_bytes = 0u)
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
    if (message_bytes >
        (UINT64_MAX - metadata_bytes - block_bytes) / 2u)
    {
        std::fprintf(stderr, "%s working-set size overflows\n", command);
        return false;
    }
    const uint64_t working_bytes =
        2u * message_bytes + block_bytes + metadata_bytes;
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

TrialResult RunBaselineTrial(
    uint32_t N,
    uint32_t block_bytes,
    double loss_rate,
    uint64_t seed)
{
    TrialResult tr = {};
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> decoded((size_t)message_bytes, 0);
    std::vector<uint8_t> block(block_bytes);
    FillMessage(message, seed);
    Rng rng(seed ^ UINT64_C(0x10fade));

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
    uint32_t write_bytes = 0;
    const uint32_t max_delivered = N * 2u + 512u;
    while (delivered < max_delivered)
    {
        const bool drop = ShouldDrop(rng, loss_rate);
        const uint32_t this_id = block_id++;
        if (drop) {
            continue;
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
    const wirehair_v2::SeedProfile* profile)
{
    TrialResult tr = {};
    const uint64_t message_bytes = (uint64_t)N * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> decoded((size_t)message_bytes, 0);
    std::vector<uint8_t> block(block_bytes);
    FillMessage(message, seed);
    Rng rng(seed ^ UINT64_C(0x10fade));

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
    uint32_t write_bytes = 0;
    const uint32_t max_delivered = N * 2u + 512u;
    while (delivered < max_delivered)
    {
        const bool drop = ShouldDrop(rng, loss_rate);
        const uint32_t this_id = block_id++;
        if (drop) {
            continue;
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
        "%-9s %-8u %-7llu %-7llu %-10.1f %-10.4f %-8.4f "
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
    CompareOptions compare_options;

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
        "dense_candidate=%u\n",
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
        compare_options.DenseCandidate);
    std::printf(
        "%-9s %-8s %-7s %-7s %-10s %-10s %-8s "
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
        std::map<uint64_t, CachedCompareProfile> profile_cache;
        for (uint32_t trial = 0; trial < trials; ++trial)
        {
            const uint32_t N =
                sample_nlo + (rng.U32() % (capped_nhi - sample_nlo + 1u));
            const uint64_t trial_seed = rng.Next();
            const wirehair_v2::SeedProfile* profile =
                SelectCompareProfile(
                    N, block_bytes, loss, compare_options, profile_cache);
            AddTrial(baseline, N,
                RunBaselineTrial(N, block_bytes, loss, trial_seed));
            AddTrial(v2, N,
                RunV2Trial(N, block_bytes, loss, trial_seed, profile));
        }
        PrintAccum("baseline", block_bytes, baseline);
        PrintAccum(
            compare_options.ProfileMode == CompareProfileBase ? "v2" :
            "v2_auto",
            block_bytes,
            v2);
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
        "# seedtable: peel_candidates=%u trials=%u seed=0x%llx\n",
        peel_candidates,
        trials,
        (unsigned long long)seed);
    std::printf(
        "N,bb,solver,structure,bucket,base_peel,tuned_peel,dense,"
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
            "%d,%d,%s,%s,%u,%u,%u,%u,%u,%u,%.4f,%u,%llu,%u,%u\n",
            n_value,
            bb_value,
            wirehair_v2::ToString(tuned.Policy.Solver),
            wirehair_v2::ToString(tuned.Policy.Structure),
            tuned.PeelSeedBucket,
            base.PeelSeed,
            tuned.PeelSeed,
            tuned.DenseSeed,
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
        "# densetune: candidates=%u trials=%u loss=%.17g seed=0x%llx\n",
        candidates,
        trials,
        loss,
        (unsigned long long)seed);
    std::printf(
        "N,bb,base_dense,best_dense,best_fail,best_oh_mean,best_oh_max,"
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
            "%d,%d,%u,%u,%llu,%.4f,%u,%llu,%.4f,%u\n",
            n_value,
            bb_value,
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

int CmdSelfTest()
{
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
            "usage: wirehair_v2_bench compare|seedtable|peelcost|densecheck|"
            "densetune|densecount|densegrid|selftest [opts]\n");
        return 1;
    }
    try
    {
        if (!std::strcmp(argv[1], "compare")) {
            return CmdCompare(argc - 2, argv + 2);
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
