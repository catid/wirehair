#include "WirehairV2Codec.h"
#include "WirehairV2Plan.h"
#include "WirehairV2Seeds.h"

#include "../WirehairTools.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double NowSeconds()
{
    return std::chrono::duration<double>(
        Clock::now().time_since_epoch()).count();
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
        return (Next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

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
    CompareProfileTuned,
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

std::vector<int> ParseIntList(const std::string& text)
{
    std::vector<int> out;
    size_t pos = 0;
    while (pos < text.size())
    {
        const size_t comma = text.find(',', pos);
        const std::string token = text.substr(
            pos, comma == std::string::npos ? std::string::npos : comma - pos);
        const int value = std::atoi(token.c_str());
        if (value > 0) {
            out.push_back(value);
        }
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
    while (pos < text.size())
    {
        const size_t comma = text.find(',', pos);
        const std::string token = text.substr(
            pos, comma == std::string::npos ? std::string::npos : comma - pos);
        if (!token.empty()) {
            out.push_back(std::atoi(token.c_str()));
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
        const bool drop = rng.Unit() < loss_rate;
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
        const bool drop = rng.Unit() < loss_rate;
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
    case CompareProfileTuned:
        return "tuned";
    case CompareProfileAuto:
        return "auto";
    }
    return "unknown";
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

    int dense_count = (int)profile.DenseCount + options.DenseDelta;
    if (dense_count < 1) {
        dense_count = 1;
    }
    if (dense_count > CAT_MAX_DENSE_ROWS) {
        dense_count = CAT_MAX_DENSE_ROWS;
    }

    profile.DenseCount = (uint16_t)dense_count;
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
    cached.TunedProfile = BuildTunedProfile(N, block_bytes, options);

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
            RunV2Trial(N, block_bytes, loss, trial_seed, 0));
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
    if (options.ProfileMode == CompareProfileTuned)
    {
        cached.UseTuned = true;
        cached.TunedProfile = BuildTunedProfile(N, block_bytes, options);
        ApplyDenseOverride(cached.TunedProfile, options);
    }
    else {
        CalibrateAutoProfile(N, block_bytes, loss, options, cached);
        if (!cached.UseTuned && options.DenseOverride)
        {
            cached.UseTuned = true;
            cached.TunedProfile =
                wirehair_v2::SelectSeedProfile(N, block_bytes);
            ApplyDenseOverride(cached.TunedProfile, options);
        }
        else if (cached.UseTuned) {
            ApplyDenseOverride(cached.TunedProfile, options);
        }
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
        if (!std::strcmp(argv[i], "--nlo") && i + 1 < argc) {
            nlo = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--nhi") && i + 1 < argc) {
            nhi = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            trials = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--bb-list") && i + 1 < argc) {
            bb_list = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--loss") && i + 1 < argc) {
            loss = std::atof(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed = std::strtoull(argv[++i], 0, 0);
        }
        else if (!std::strcmp(argv[i], "--max-message-mib") && i + 1 < argc) {
            max_message_mib = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--v2-profile") && i + 1 < argc) {
            const char* profile = argv[++i];
            if (!std::strcmp(profile, "tuned")) {
                compare_options.ProfileMode = CompareProfileTuned;
            }
            else if (!std::strcmp(profile, "auto")) {
                compare_options.ProfileMode = CompareProfileAuto;
            }
            else {
                compare_options.ProfileMode = CompareProfileBase;
            }
        }
        else if (!std::strcmp(argv[i], "--peel-candidates") && i + 1 < argc) {
            compare_options.PeelCandidates = (uint16_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--peel-trials") && i + 1 < argc) {
            compare_options.PeelTrials = (uint16_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--auto-trials") && i + 1 < argc) {
            compare_options.AutoTrials = (uint16_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--tune-seed") && i + 1 < argc) {
            compare_options.TuneSeed = std::strtoull(argv[++i], 0, 0);
        }
        else if (!std::strcmp(argv[i], "--auto-seed") && i + 1 < argc) {
            compare_options.AutoSeed = std::strtoull(argv[++i], 0, 0);
        }
        else if (!std::strcmp(argv[i], "--auto-min-delta") && i + 1 < argc) {
            compare_options.AutoMinDelta = std::atof(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--dense-delta") && i + 1 < argc) {
            compare_options.DenseOverride = true;
            compare_options.DenseDelta = std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--dense-candidate") && i + 1 < argc) {
            compare_options.DenseOverride = true;
            compare_options.DenseCandidate = (uint16_t)std::atoi(argv[++i]);
        }
    }

    const std::vector<int> block_bytes_list = ParseIntList(bb_list);
    if (block_bytes_list.empty() || nlo < 2u || nhi < nlo || trials == 0u) {
        std::fprintf(stderr, "invalid compare arguments\n");
        return 1;
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
    compare_options.LogAutoChoices =
        compare_options.ProfileMode == CompareProfileAuto && nlo == nhi;

    std::printf(
        "# compare: N=[%u,%u] trials/bb=%u loss=%.2f seed=0x%llx "
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
            const uint64_t max_bytes =
                (uint64_t)max_message_mib * 1024u * 1024u;
            uint32_t max_n = (uint32_t)(max_bytes / block_bytes);
            if (max_n < 2u) {
                max_n = 2u;
            }
            if (capped_nhi > max_n) {
                capped_nhi = max_n;
            }
            if (sample_nlo > max_n) {
                sample_nlo = max_n;
            }
        }
        if (capped_nhi < sample_nlo) {
            capped_nhi = sample_nlo;
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
            compare_options.ProfileMode == CompareProfileTuned ? "v2_tuned" :
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
        if (!std::strcmp(argv[i], "--N") && i + 1 < argc) {
            nlist = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--bb-list") && i + 1 < argc) {
            bb_list = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--peel-candidates") && i + 1 < argc) {
            peel_candidates = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            trials = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed = std::strtoull(argv[++i], 0, 0);
        }
    }

    if (peel_candidates < 1u) {
        peel_candidates = 1u;
    }
    if (peel_candidates > 256u) {
        peel_candidates = 256u;
    }
    if (trials < 1u) {
        trials = 1u;
    }

    const std::vector<int> Ns = ParseIntList(nlist);
    const std::vector<int> BBs = ParseIntList(bb_list);
    if (Ns.empty() || BBs.empty()) {
        std::fprintf(stderr, "seedtable requires non-empty --N and --bb-list\n");
        return 1;
    }

    std::printf(
        "# seedtable: peel_candidates=%u trials=%u seed=0x%llx\n",
        peel_candidates,
        trials,
        (unsigned long long)seed);
    std::printf(
        "N,bb,solver,structure,bucket,base_peel,tuned_peel,dense,"
        "resid_mean,resid_max,xor_cost,used_peel_fixup,used_dense_fixup\n");

    for (int bb_value : BBs) for (int n_value : Ns)
    {
        wirehair_v2::SeedTuningOptions options =
            wirehair_v2::DefaultSeedTuningOptions();
        options.PeelCandidates = (uint16_t)peel_candidates;
        options.TrialsPerCandidate = (uint16_t)trials;
        options.Seed = seed;

        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile((uint32_t)n_value, (uint32_t)bb_value);
        const wirehair_v2::SeedProfile tuned =
            wirehair_v2::TuneSeedProfile((uint32_t)n_value, (uint32_t)bb_value,
                options);
        std::printf(
            "%d,%d,%s,%s,%u,%u,%u,%u,%.4f,%u,%llu,%u,%u\n",
            n_value,
            bb_value,
            wirehair_v2::ToString(tuned.Policy.Solver),
            wirehair_v2::ToString(tuned.Policy.Structure),
            tuned.PeelSeedBucket,
            base.PeelSeed,
            tuned.PeelSeed,
            tuned.DenseSeed,
            tuned.TuningResidualMean,
            tuned.TuningResidualColumns,
            (unsigned long long)tuned.TuningXorCost,
            tuned.UsedPeelFixup ? 1u : 0u,
            tuned.UsedDenseFixup ? 1u : 0u);
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
        if (!std::strcmp(argv[i], "--N") && i + 1 < argc) {
            N = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--bb") && i + 1 < argc) {
            block_bytes = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--candidates") && i + 1 < argc) {
            candidates = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            trials = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--loss") && i + 1 < argc) {
            loss = std::atof(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed = std::strtoull(argv[++i], 0, 0);
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
    if (N < 1u || block_bytes < 1u) {
        std::fprintf(stderr, "densecheck requires positive --N and --bb\n");
        return 1;
    }
    const wirehair_v2::SeedProfile base =
        wirehair_v2::SelectSeedProfile(N, block_bytes);
    std::printf(
        "# densecheck: N=%u bb=%u base_dense=%u candidates=%u trials=%u\n",
        N, block_bytes, base.DenseSeed, candidates, trials);
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
            const uint64_t trial_seed = seed ^
                ((uint64_t)c * UINT64_C(0x9e3779b97f4a7c15)) ^
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
        if (!std::strcmp(argv[i], "--N") && i + 1 < argc) {
            nlist = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--bb-list") && i + 1 < argc) {
            bb_list = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--candidates") && i + 1 < argc) {
            candidates = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            trials = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--loss") && i + 1 < argc) {
            loss = std::atof(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed = std::strtoull(argv[++i], 0, 0);
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

    std::printf(
        "# densetune: candidates=%u trials=%u loss=%.2f seed=0x%llx\n",
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
                const uint64_t trial_seed = seed ^
                    ((uint64_t)n_value * UINT64_C(0x9e3779b97f4a7c15)) ^
                    ((uint64_t)bb_value * UINT64_C(0xbf58476d1ce4e5b9)) ^
                    ((uint64_t)c * UINT64_C(0x94d049bb133111eb)) ^
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
        if (!std::strcmp(argv[i], "--N") && i + 1 < argc) {
            nlist = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--bb-list") && i + 1 < argc) {
            bb_list = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--deltas") && i + 1 < argc) {
            delta_list = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            trials = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--loss") && i + 1 < argc) {
            loss = std::atof(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed = std::strtoull(argv[++i], 0, 0);
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

    std::printf(
        "# densecount: trials=%u loss=%.2f seed=0x%llx\n",
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
            int dense_count = (int)base.DenseCount + delta;
            if (dense_count < 1) {
                dense_count = 1;
            }
            if (dense_count > CAT_MAX_DENSE_ROWS) {
                dense_count = CAT_MAX_DENSE_ROWS;
            }

            wirehair_v2::SeedProfile profile = base;
            profile.DenseCount = (uint16_t)dense_count;
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
        if (!std::strcmp(argv[i], "--N") && i + 1 < argc) {
            nlist = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--bb-list") && i + 1 < argc) {
            bb_list = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--deltas") && i + 1 < argc) {
            delta_list = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--candidates") && i + 1 < argc) {
            candidates = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--trials") && i + 1 < argc) {
            trials = (uint32_t)std::atoi(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--loss") && i + 1 < argc) {
            loss = std::atof(argv[++i]);
        }
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed = std::strtoull(argv[++i], 0, 0);
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

    std::printf(
        "# densegrid: candidates=%u trials=%u loss=%.2f seed=0x%llx\n",
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
            int dense_count = (int)base.DenseCount + delta;
            if (dense_count < 1) {
                dense_count = 1;
            }
            if (dense_count > CAT_MAX_DENSE_ROWS) {
                dense_count = CAT_MAX_DENSE_ROWS;
            }
            const uint16_t dense_count_u = (uint16_t)dense_count;
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

} // namespace

int main(int argc, char** argv)
{
    if (wirehair_init() != Wirehair_Success) {
        std::fprintf(stderr, "wirehair_init failed\n");
        return 2;
    }

    if (argc < 2) {
        std::fprintf(stderr,
            "usage: wirehair_v2_bench compare|seedtable|densecheck|densetune|densecount|densegrid [opts]\n");
        return 1;
    }
    if (!std::strcmp(argv[1], "compare")) {
        return CmdCompare(argc - 2, argv + 2);
    }
    if (!std::strcmp(argv[1], "seedtable")) {
        return CmdSeedTable(argc - 2, argv + 2);
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
    std::fprintf(stderr, "unknown mode: %s\n", argv[1]);
    return 1;
}
