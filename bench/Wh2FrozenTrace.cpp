#include "Wh2FrozenTrace.h"

#include <algorithm>
#include <climits>
#include <cstring>
#include <limits>

namespace wirehair {
namespace wh2_benchmark {
namespace {

static const uint64_t kSplitMixIncrement = UINT64_C(0x9e3779b97f4a7c15);
static const uint64_t kSplitMixMultiplier1 = UINT64_C(0xbf58476d1ce4e5b9);
static const uint64_t kSplitMixMultiplier2 = UINT64_C(0x94d049bb133111eb);
static const uint64_t kHardScheduleSalt = UINT64_C(0x10fade);
static const uint32_t kMaximumFrozenK = 64000u;
static const uint32_t kOverheadCap = 4u;

static_assert(CHAR_BIT == 8, "frozen trace requires eight-bit bytes");
static_assert(std::numeric_limits<double>::radix == 2 &&
        std::numeric_limits<double>::digits == 53,
    "frozen trace requires binary64-compatible double precision");

class Sha256
{
public:
    Sha256()
        : total_bytes_(0u)
        , buffered_(0u)
    {
        state_[0] = UINT32_C(0x6a09e667);
        state_[1] = UINT32_C(0xbb67ae85);
        state_[2] = UINT32_C(0x3c6ef372);
        state_[3] = UINT32_C(0xa54ff53a);
        state_[4] = UINT32_C(0x510e527f);
        state_[5] = UINT32_C(0x9b05688c);
        state_[6] = UINT32_C(0x1f83d9ab);
        state_[7] = UINT32_C(0x5be0cd19);
        std::memset(block_, 0, sizeof(block_));
    }

    void Update(const uint8_t* data, std::size_t bytes)
    {
        if (bytes == 0u) {
            return;
        }

        total_bytes_ += static_cast<uint64_t>(bytes);

        if (buffered_ != 0u)
        {
            const std::size_t take = std::min(bytes, sizeof(block_) - buffered_);
            std::memcpy(block_ + buffered_, data, take);
            buffered_ += take;
            data += take;
            bytes -= take;
            if (buffered_ == sizeof(block_)) {
                Transform(block_);
                buffered_ = 0u;
            }
        }

        while (bytes >= sizeof(block_))
        {
            Transform(data);
            data += sizeof(block_);
            bytes -= sizeof(block_);
        }

        if (bytes != 0u) {
            std::memcpy(block_, data, bytes);
            buffered_ = bytes;
        }
    }

    void Final(uint8_t digest[32])
    {
        const uint64_t message_bits = total_bytes_ * UINT64_C(8);

        block_[buffered_++] = 0x80u;
        if (buffered_ > 56u) {
            std::memset(block_ + buffered_, 0, sizeof(block_) - buffered_);
            Transform(block_);
            buffered_ = 0u;
        }
        std::memset(block_ + buffered_, 0, 56u - buffered_);
        for (unsigned i = 0u; i < 8u; ++i) {
            block_[63u - i] = static_cast<uint8_t>(message_bits >> (i * 8u));
        }
        Transform(block_);

        for (unsigned i = 0u; i < 8u; ++i)
        {
            digest[i * 4u] = static_cast<uint8_t>(state_[i] >> 24);
            digest[i * 4u + 1u] = static_cast<uint8_t>(state_[i] >> 16);
            digest[i * 4u + 2u] = static_cast<uint8_t>(state_[i] >> 8);
            digest[i * 4u + 3u] = static_cast<uint8_t>(state_[i]);
        }
    }

private:
    static uint32_t RotateRight(uint32_t value, unsigned shift)
    {
        return (value >> shift) | (value << (32u - shift));
    }

    void Transform(const uint8_t block[64])
    {
        static const uint32_t round_constants[64] = {
            UINT32_C(0x428a2f98), UINT32_C(0x71374491),
            UINT32_C(0xb5c0fbcf), UINT32_C(0xe9b5dba5),
            UINT32_C(0x3956c25b), UINT32_C(0x59f111f1),
            UINT32_C(0x923f82a4), UINT32_C(0xab1c5ed5),
            UINT32_C(0xd807aa98), UINT32_C(0x12835b01),
            UINT32_C(0x243185be), UINT32_C(0x550c7dc3),
            UINT32_C(0x72be5d74), UINT32_C(0x80deb1fe),
            UINT32_C(0x9bdc06a7), UINT32_C(0xc19bf174),
            UINT32_C(0xe49b69c1), UINT32_C(0xefbe4786),
            UINT32_C(0x0fc19dc6), UINT32_C(0x240ca1cc),
            UINT32_C(0x2de92c6f), UINT32_C(0x4a7484aa),
            UINT32_C(0x5cb0a9dc), UINT32_C(0x76f988da),
            UINT32_C(0x983e5152), UINT32_C(0xa831c66d),
            UINT32_C(0xb00327c8), UINT32_C(0xbf597fc7),
            UINT32_C(0xc6e00bf3), UINT32_C(0xd5a79147),
            UINT32_C(0x06ca6351), UINT32_C(0x14292967),
            UINT32_C(0x27b70a85), UINT32_C(0x2e1b2138),
            UINT32_C(0x4d2c6dfc), UINT32_C(0x53380d13),
            UINT32_C(0x650a7354), UINT32_C(0x766a0abb),
            UINT32_C(0x81c2c92e), UINT32_C(0x92722c85),
            UINT32_C(0xa2bfe8a1), UINT32_C(0xa81a664b),
            UINT32_C(0xc24b8b70), UINT32_C(0xc76c51a3),
            UINT32_C(0xd192e819), UINT32_C(0xd6990624),
            UINT32_C(0xf40e3585), UINT32_C(0x106aa070),
            UINT32_C(0x19a4c116), UINT32_C(0x1e376c08),
            UINT32_C(0x2748774c), UINT32_C(0x34b0bcb5),
            UINT32_C(0x391c0cb3), UINT32_C(0x4ed8aa4a),
            UINT32_C(0x5b9cca4f), UINT32_C(0x682e6ff3),
            UINT32_C(0x748f82ee), UINT32_C(0x78a5636f),
            UINT32_C(0x84c87814), UINT32_C(0x8cc70208),
            UINT32_C(0x90befffa), UINT32_C(0xa4506ceb),
            UINT32_C(0xbef9a3f7), UINT32_C(0xc67178f2)
        };

        uint32_t words[64];
        for (unsigned i = 0u; i < 16u; ++i)
        {
            const unsigned offset = i * 4u;
            words[i] =
                (static_cast<uint32_t>(block[offset]) << 24) |
                (static_cast<uint32_t>(block[offset + 1u]) << 16) |
                (static_cast<uint32_t>(block[offset + 2u]) << 8) |
                static_cast<uint32_t>(block[offset + 3u]);
        }
        for (unsigned i = 16u; i < 64u; ++i)
        {
            const uint32_t s0 =
                RotateRight(words[i - 15u], 7u) ^
                RotateRight(words[i - 15u], 18u) ^
                (words[i - 15u] >> 3);
            const uint32_t s1 =
                RotateRight(words[i - 2u], 17u) ^
                RotateRight(words[i - 2u], 19u) ^
                (words[i - 2u] >> 10);
            words[i] = words[i - 16u] + s0 + words[i - 7u] + s1;
        }

        uint32_t a = state_[0];
        uint32_t b = state_[1];
        uint32_t c = state_[2];
        uint32_t d = state_[3];
        uint32_t e = state_[4];
        uint32_t f = state_[5];
        uint32_t g = state_[6];
        uint32_t h = state_[7];

        for (unsigned i = 0u; i < 64u; ++i)
        {
            const uint32_t sum1 =
                RotateRight(e, 6u) ^ RotateRight(e, 11u) ^
                RotateRight(e, 25u);
            const uint32_t choose = (e & f) ^ ((~e) & g);
            const uint32_t temporary1 =
                h + sum1 + choose + round_constants[i] + words[i];
            const uint32_t sum0 =
                RotateRight(a, 2u) ^ RotateRight(a, 13u) ^
                RotateRight(a, 22u);
            const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t temporary2 = sum0 + majority;

            h = g;
            g = f;
            f = e;
            e = d + temporary1;
            d = c;
            c = b;
            b = a;
            a = temporary1 + temporary2;
        }

        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
        state_[4] += e;
        state_[5] += f;
        state_[6] += g;
        state_[7] += h;
    }

    uint32_t state_[8];
    uint64_t total_bytes_;
    std::size_t buffered_;
    uint8_t block_[64];
};

class SplitMix64
{
public:
    explicit SplitMix64(uint64_t seed)
        : state_(seed)
    {
    }

    uint64_t Next()
    {
        uint64_t value = (state_ += kSplitMixIncrement);
        value = (value ^ (value >> 30)) * kSplitMixMultiplier1;
        value = (value ^ (value >> 27)) * kSplitMixMultiplier2;
        return value ^ (value >> 31);
    }

    double Unit()
    {
        return static_cast<double>(Next() >> 11) *
            (1.0 / 9007199254740992.0);
    }

private:
    uint64_t state_;
};

bool IsFrozenSchedule(FrozenSchedule schedule)
{
    switch (schedule)
    {
    case FrozenSchedule::Iid:
    case FrozenSchedule::Burst:
    case FrozenSchedule::Adversarial:
    case FrozenSchedule::RepairOnly:
        return true;
    }
    return false;
}

std::string DigestHex(const uint8_t digest[32])
{
    static const char hex[] = "0123456789abcdef";
    std::string text(64u, '0');
    for (unsigned i = 0u; i < 32u; ++i) {
        text[i * 2u] = hex[digest[i] >> 4];
        text[i * 2u + 1u] = hex[digest[i] & 15u];
    }
    return text;
}

std::string FinishSha256(Sha256& sha)
{
    uint8_t digest[32];
    sha.Final(digest);
    return DigestHex(digest);
}

std::string HexSeed(uint64_t seed)
{
    static const char hex[] = "0123456789abcdef";
    std::string text("0x0000000000000000");
    for (unsigned i = 0u; i < 16u; ++i) {
        text[17u - i] = hex[seed & 15u];
        seed >>= 4;
    }
    return text;
}

std::string BandForK(uint32_t K)
{
    if (K >= 2u && K <= 100u) {
        return "2-100";
    }
    if (K <= 1000u) {
        return "101-1000";
    }
    if (K <= 5000u) {
        return "1001-5000";
    }
    if (K <= 10000u) {
        return "5001-10000";
    }
    if (K <= 20000u) {
        return "10001-20000";
    }
    if (K <= kMaximumFrozenK) {
        return "20001-64000";
    }
    return std::string();
}

uint32_t CandidatePacketId(
    FrozenSchedule schedule,
    uint32_t K,
    uint64_t candidate_index)
{
    if (schedule == FrozenSchedule::Adversarial) {
        const uint32_t doubled =
            static_cast<uint32_t>(candidate_index * UINT64_C(2));
        return UINT32_MAX - doubled;
    }
    if (schedule == FrozenSchedule::RepairOnly) {
        return static_cast<uint32_t>(
            static_cast<uint64_t>(K) + candidate_index);
    }
    return static_cast<uint32_t>(candidate_index);
}

bool ExactDevelopmentCell(const FrozenRecoveryCell& cell)
{
    static const uint32_t short_K[] = {
        2u, 3u, 4u, 5u, 6u, 8u, 16u, 32u, 64u, 100u,
        101u, 128u, 256u, 512u, 513u, 1000u,
        1001u, 2048u, 4096u, 5000u,
        5001u, 8192u, 10000u,
        10001u, 16384u, 20000u,
        20001u, 32768u, 49152u, 64000u
    };
    static const uint64_t roots[] = {
        UINT64_C(0xd1b54a32d192ed03),
        UINT64_C(0x94d049bb133111eb),
        UINT64_C(0x8538ecb5bd456ea3)
    };
    static const FrozenSchedule schedules[] = {
        FrozenSchedule::Iid,
        FrozenSchedule::Burst,
        FrozenSchedule::Adversarial,
        FrozenSchedule::RepairOnly
    };
    static const uint32_t losses[] = {
        100000u, 500000u, 500000u, 500000u
    };

    const uint32_t* const k_end =
        short_K + sizeof(short_K) / sizeof(short_K[0]);
    const uint32_t* const k_position = std::find(
        short_K,
        k_end,
        cell.K);
    if (k_position == k_end || cell.trial >= 3u ||
        cell.stratum_index >= 4u)
    {
        return false;
    }
    const std::size_t expected_ordinal =
        static_cast<std::size_t>(cell.trial) * 120u +
        static_cast<std::size_t>(cell.stratum_index) * 30u +
        static_cast<std::size_t>(k_position - short_K);
    return cell.ordinal == expected_ordinal &&
        cell.phase == "development" &&
        cell.band == BandForK(cell.K) &&
        cell.block_bytes == 2u &&
        cell.base_seed_attempt == cell.trial &&
        cell.loss_seed == roots[cell.trial] &&
        cell.schedule == schedules[cell.stratum_index] &&
        cell.loss_ppm == losses[cell.stratum_index] &&
        cell.overhead_cap == kOverheadCap;
}

} // namespace

FrozenRecoveryCell::FrozenRecoveryCell()
    : ordinal(0u)
    , K(0u)
    , block_bytes(0u)
    , loss_ppm(0u)
    , schedule(FrozenSchedule::Iid)
    , stratum_index(0u)
    , trial(0u)
    , base_seed_attempt(0u)
    , loss_seed(0u)
    , overhead_cap(0u)
{
}

const char* FrozenScheduleName(FrozenSchedule schedule)
{
    switch (schedule)
    {
    case FrozenSchedule::Iid:         return "iid";
    case FrozenSchedule::Burst:       return "burst";
    case FrozenSchedule::Adversarial: return "adversarial";
    case FrozenSchedule::RepairOnly:  return "repair-only";
    }
    return "";
}

std::string Sha256Hex(const void* data, std::size_t bytes)
{
    if (data == NULL && bytes != 0u) {
        return std::string();
    }
    Sha256 sha;
    if (bytes != 0u) {
        sha.Update(static_cast<const uint8_t*>(data), bytes);
    }
    return FinishSha256(sha);
}

std::string Sha256Hex(const std::string& data)
{
    return Sha256Hex(data.data(), data.size());
}

std::string PacketIdsSha256(const std::vector<uint32_t>& packet_ids)
{
    Sha256 sha;
    for (std::size_t i = 0u; i < packet_ids.size(); ++i)
    {
        const uint32_t id = packet_ids[i];
        const uint8_t little_endian[4] = {
            static_cast<uint8_t>(id),
            static_cast<uint8_t>(id >> 8),
            static_cast<uint8_t>(id >> 16),
            static_cast<uint8_t>(id >> 24)
        };
        sha.Update(little_endian, sizeof(little_endian));
    }
    return FinishSha256(sha);
}

std::vector<FrozenRecoveryCell> EnumerateDevelopmentRecoveryCells()
{
    static const uint32_t short_K[] = {
        2u, 3u, 4u, 5u, 6u, 8u, 16u, 32u, 64u, 100u,
        101u, 128u, 256u, 512u, 513u, 1000u,
        1001u, 2048u, 4096u, 5000u,
        5001u, 8192u, 10000u,
        10001u, 16384u, 20000u,
        20001u, 32768u, 49152u, 64000u
    };
    static const uint64_t roots[] = {
        UINT64_C(0xd1b54a32d192ed03),
        UINT64_C(0x94d049bb133111eb),
        UINT64_C(0x8538ecb5bd456ea3)
    };
    static const FrozenSchedule schedules[] = {
        FrozenSchedule::Iid,
        FrozenSchedule::Burst,
        FrozenSchedule::Adversarial,
        FrozenSchedule::RepairOnly
    };
    static const uint32_t losses[] = {
        100000u, 500000u, 500000u, 500000u
    };

    std::vector<FrozenRecoveryCell> cells;
    cells.reserve(360u);
    for (uint32_t trial = 0u; trial < 3u; ++trial)
    {
        for (uint32_t stratum = 0u; stratum < 4u; ++stratum)
        {
            for (std::size_t k_index = 0u;
                k_index < sizeof(short_K) / sizeof(short_K[0]);
                ++k_index)
            {
                FrozenRecoveryCell cell;
                cell.ordinal = cells.size();
                cell.phase = "development";
                cell.band = BandForK(short_K[k_index]);
                cell.K = short_K[k_index];
                cell.block_bytes = 2u;
                cell.loss_ppm = losses[stratum];
                cell.schedule = schedules[stratum];
                cell.stratum_index = stratum;
                cell.trial = trial;
                cell.base_seed_attempt = trial;
                cell.loss_seed = roots[trial];
                cell.overhead_cap = kOverheadCap;
                cells.push_back(cell);
            }
        }
    }
    return cells;
}

std::string CanonicalRecoveryCellJson(const FrozenRecoveryCell& cell)
{
    if (!ExactDevelopmentCell(cell)) {
        return std::string();
    }

    // Key order is alphabetical and intentionally written out rather than
    // delegated to a locale- or library-dependent JSON serializer.
    std::string json;
    json.reserve(224u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"band\":\"";
    json += cell.band;
    json += "\",\"base_seed_attempt\":";
    json += std::to_string(cell.base_seed_attempt);
    json += ",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json += ",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"overhead_cap\":";
    json += std::to_string(cell.overhead_cap);
    json += ",\"phase\":\"";
    json += cell.phase;
    json += "\",\"schedule\":\"";
    json += FrozenScheduleName(cell.schedule);
    json += "\",\"trial\":";
    json += std::to_string(cell.trial);
    json += '}';
    return json;
}

std::string RecoveryCellSha256(const FrozenRecoveryCell& cell)
{
    const std::string json = CanonicalRecoveryCellJson(cell);
    if (json.empty()) {
        return std::string();
    }
    return Sha256Hex(json);
}

std::string DevelopmentRecoveryDomainSha256()
{
    const std::vector<FrozenRecoveryCell> cells =
        EnumerateDevelopmentRecoveryCells();
    Sha256 sha;
    static const uint8_t newline = '\n';
    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        const std::string json = CanonicalRecoveryCellJson(cells[i]);
        if (json.empty()) {
            return std::string();
        }
        sha.Update(reinterpret_cast<const uint8_t*>(json.data()), json.size());
        sha.Update(&newline, 1u);
    }
    return FinishSha256(sha);
}

uint64_t FrozenTraceSeed(
    uint32_t K,
    uint32_t block_bytes,
    uint64_t loss_seed)
{
    return loss_seed ^
        (static_cast<uint64_t>(K) * kSplitMixIncrement) ^
        (static_cast<uint64_t>(block_bytes) * kSplitMixMultiplier1);
}

bool FrozenCandidateLimit(
    uint64_t delivered_count,
    uint64_t& candidate_limit)
{
    candidate_limit = 0u;
    const uint64_t fixed_allowance = UINT64_C(65536);
    const uint64_t multiplier = UINT64_C(256);
    if (delivered_count >
        (std::numeric_limits<uint64_t>::max() - fixed_allowance) / multiplier)
    {
        return false;
    }
    candidate_limit = delivered_count * multiplier + fixed_allowance;
    return true;
}

FrozenPacketTrace::FrozenPacketTrace()
    : status(FrozenTraceStatus::InvalidInput)
    , K(0u)
    , block_bytes(0u)
    , loss_ppm(0u)
    , schedule(FrozenSchedule::Iid)
    , loss_seed(0u)
    , trace_seed(0u)
    , attempted_candidates(0u)
    , candidate_limit(0u)
{
}

FrozenTraceStatus GenerateFrozenPacketTrace(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t loss_ppm,
    uint64_t loss_seed,
    FrozenSchedule schedule,
    FrozenPacketTrace& trace)
{
    trace = FrozenPacketTrace();
    trace.K = K;
    trace.block_bytes = block_bytes;
    trace.loss_ppm = loss_ppm;
    trace.schedule = schedule;
    trace.loss_seed = loss_seed;

    if (K < 2u || K > kMaximumFrozenK || block_bytes == 0u ||
        loss_ppm > 1000000u || !IsFrozenSchedule(schedule))
    {
        return trace.status;
    }

    const uint64_t delivered_count = static_cast<uint64_t>(K) + kOverheadCap;
    if (!FrozenCandidateLimit(delivered_count, trace.candidate_limit) ||
        delivered_count >
            static_cast<uint64_t>(trace.delivered_ids.max_size()) ||
        delivered_count >
            static_cast<uint64_t>(std::numeric_limits<std::size_t>::max()))
    {
        return trace.status;
    }

    trace.trace_seed = FrozenTraceSeed(K, block_bytes, loss_seed);
    const uint64_t rng_seed = schedule == FrozenSchedule::Iid ?
        trace.trace_seed : trace.trace_seed ^ kHardScheduleSalt;
    SplitMix64 rng(rng_seed);
    trace.delivered_ids.reserve(static_cast<std::size_t>(delivered_count));

    const double loss_rate =
        static_cast<double>(loss_ppm) / 1000000.0;
    uint32_t burst_remaining = 0u;

    while (trace.delivered_ids.size() < delivered_count &&
        trace.attempted_candidates < trace.candidate_limit)
    {
        const uint64_t candidate_index = trace.attempted_candidates++;
        const uint32_t id = CandidatePacketId(schedule, K, candidate_index);

        bool drop = false;
        if (schedule == FrozenSchedule::Burst)
        {
            static const uint32_t burst_length = 8u;
            if (burst_remaining != 0u) {
                --burst_remaining;
                drop = true;
            }
            else
            {
                const double start_probability = loss_rate /
                    (static_cast<double>(burst_length) -
                        static_cast<double>(burst_length - 1u) * loss_rate);
                if (rng.Unit() < start_probability) {
                    burst_remaining = burst_length - 1u;
                    drop = true;
                }
            }
        }
        else {
            drop = rng.Unit() < loss_rate;
        }

        if (!drop) {
            trace.delivered_ids.push_back(id);
        }
    }

    if (trace.delivered_ids.size() == delivered_count) {
        trace.status = FrozenTraceStatus::Complete;
        trace.trace_sha256 = PacketIdsSha256(trace.delivered_ids);
    }
    else {
        trace.status = FrozenTraceStatus::CandidateLimitReached;
    }
    return trace.status;
}

FrozenTraceStatus GenerateFrozenPacketTrace(
    const FrozenRecoveryCell& cell,
    FrozenPacketTrace& trace)
{
    if (!ExactDevelopmentCell(cell)) {
        trace = FrozenPacketTrace();
        return trace.status;
    }
    return GenerateFrozenPacketTrace(
        cell.K,
        cell.block_bytes,
        cell.loss_ppm,
        cell.loss_seed,
        cell.schedule,
        trace);
}

bool CopyNestedPrefix(
    const FrozenPacketTrace& trace,
    uint32_t overhead,
    std::vector<uint32_t>& prefix)
{
    prefix.clear();
    if (trace.status != FrozenTraceStatus::Complete ||
        overhead > kOverheadCap ||
        trace.K < 2u || trace.K > kMaximumFrozenK ||
        trace.delivered_ids.size() !=
            static_cast<std::size_t>(trace.K) + kOverheadCap ||
        trace.trace_sha256.empty() ||
        trace.trace_sha256 != PacketIdsSha256(trace.delivered_ids))
    {
        return false;
    }

    const std::size_t prefix_size =
        static_cast<std::size_t>(trace.K) + overhead;
    const std::vector<uint32_t>::difference_type prefix_distance =
        static_cast<std::vector<uint32_t>::difference_type>(prefix_size);
    prefix.assign(
        trace.delivered_ids.begin(),
        trace.delivered_ids.begin() + prefix_distance);
    return true;
}

} // namespace wh2_benchmark
} // namespace wirehair
