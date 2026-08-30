#include "Wh2FrozenTrace.h"
#include "Wh2K3HyperovalCodec.h"

#include "../gf256.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

// This executable is deliberately clock-free.  It is a bounded structural
// falsifier, not a benchmark or a seed-search tool.  The manifest path is also
// deliberately equation-free: it generates and hashes raw-ID corpora and
// frozen traces, but never constructs a codec or evaluates a determinant.

#ifndef WIREHAIR_WH2_K3_LOW_CORPUS_SHA256
#define WIREHAIR_WH2_K3_LOW_CORPUS_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_HIGH_IID_SHA256
#define WIREHAIR_WH2_K3_HIGH_IID_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_SEQUENTIAL_SHA256
#define WIREHAIR_WH2_K3_SEQUENTIAL_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_DELTA257_SHA256
#define WIREHAIR_WH2_K3_DELTA257_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_DELTA65535_SHA256
#define WIREHAIR_WH2_K3_DELTA65535_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_DESCENDING_SHA256
#define WIREHAIR_WH2_K3_DESCENDING_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_ONE_LOW_SHA256
#define WIREHAIR_WH2_K3_ONE_LOW_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_TWO_LOW_SHA256
#define WIREHAIR_WH2_K3_TWO_LOW_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_CANONICAL_SHA256
#define WIREHAIR_WH2_K3_CANONICAL_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_TRACE_SHA256
#define WIREHAIR_WH2_K3_TRACE_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_K3_COMPLEMENT_SHA256
#define WIREHAIR_WH2_K3_COMPLEMENT_SHA256 ""
#endif
#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#define WIREHAIR_WH2_SOURCE_GIT_COMMIT "unknown"
#endif
#ifndef WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256
#define WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256 "unknown"
#endif
#ifndef WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256
#define WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256 "unknown"
#endif
#ifndef WIREHAIR_WH2_K3_CODEC_HEADER_SHA256
#define WIREHAIR_WH2_K3_CODEC_HEADER_SHA256 "unknown"
#endif

namespace {

using wirehair_wh2_bench::K3HyperovalCodec;
using wirehair_wh2_bench::K3HyperovalDirection;

static const char kProtocolTag[] =
    "wh2-k3-hyperoval-rank-recovery-v1";
static const char kCorpusSchema[] =
    "wirehair.wh2.k3.hyperoval.corpus.v1";
static const char kComplementSchema[] =
    "wirehair.wh2.k3.hyperoval.complement.v1";
static const char kDevelopmentDomainSha256[] =
    "f97f28c211428cd77aed97160073b192d93014cb4a61a844bc7d76375ac61b77";

static const uint32_t kK = 3u;
static const uint32_t kComparatorBlockBytes = 64u;
static const uint64_t kComparatorMessageBytes = 192u;
static const uint64_t kStochasticTriples = UINT64_C(1048576);
static const uint64_t kLowTriples = UINT64_C(2829056);
static const uint32_t kWorkerCount = 32u;
static const uint64_t kHighIidRoot = UINT64_C(0x243f6a8885a308d3);
static const long double kWilsonZ =
    1.959963984540054L;

struct Direction
{
    uint8_t A = 0u;
    uint8_t B = 0u;
    uint8_t C = 0u;
    bool Nonzero = false;
};

struct Triple
{
    uint32_t Id[3] = { 0u, 0u, 0u };
};

enum class Family : uint32_t
{
    HighIid = 0u,
    Sequential,
    Delta257,
    Delta65535,
    Descending,
    OneLowTwoHigh,
    TwoLowOneHigh,
    Count
};

const char* FamilyName(Family family)
{
    switch (family)
    {
    case Family::HighIid:       return "high-iid";
    case Family::Sequential:    return "sequential";
    case Family::Delta257:      return "delta257";
    case Family::Delta65535:    return "delta65535";
    case Family::Descending:    return "descending";
    case Family::OneLowTwoHigh: return "one-low-two-high";
    case Family::TwoLowOneHigh: return "two-low-one-high";
    case Family::Count: break;
    }
    return "invalid";
}

uint8_t ScalarMultiply(uint8_t left, uint8_t right) noexcept
{
    uint8_t result = 0u;
    for (unsigned bit = 0u; bit != 8u; ++bit)
    {
        if ((right & 1u) != 0u) result ^= left;
        const bool carry = (left & 0x80u) != 0u;
        left = static_cast<uint8_t>(left << 1u);
        if (carry) left ^= UINT8_C(0x4d); // x^8 = x^6+x^3+x^2+1
        right = static_cast<uint8_t>(right >> 1u);
    }
    return result;
}

uint8_t ScalarPower(uint8_t value, uint32_t exponent) noexcept
{
    uint8_t result = 1u;
    while (exponent != 0u)
    {
        if ((exponent & 1u) != 0u) {
            result = ScalarMultiply(result, value);
        }
        value = ScalarMultiply(value, value);
        exponent >>= 1u;
    }
    return result;
}

uint8_t ScalarInverse(uint8_t value) noexcept
{
    return value == 0u ? 0u : ScalarPower(value, 254u);
}

Direction NormalizeDirection(uint8_t a, uint8_t b, uint8_t c) noexcept
{
    Direction direction;
    const uint8_t pivot = a != 0u ? a : (b != 0u ? b : c);
    if (pivot == 0u) return direction;
    const uint8_t inverse = ScalarInverse(pivot);
    direction.A = ScalarMultiply(a, inverse);
    direction.B = ScalarMultiply(b, inverse);
    direction.C = ScalarMultiply(c, inverse);
    direction.Nonzero = true;
    return direction;
}

bool SameDirection(const Direction& left, const Direction& right) noexcept
{
    return left.Nonzero == right.Nonzero && left.A == right.A &&
        left.B == right.B && left.C == right.C;
}

uint32_t DirectionKey(const Direction& direction) noexcept
{
    return static_cast<uint32_t>(direction.A) |
        (static_cast<uint32_t>(direction.B) << 8u) |
        (static_cast<uint32_t>(direction.C) << 16u);
}

uint8_t ScalarDeterminant(
    const Direction& first,
    const Direction& second,
    const Direction& third) noexcept
{
    if (!first.Nonzero || !second.Nonzero || !third.Nonzero) return 0u;
    const uint8_t ei_fh = ScalarMultiply(second.B, third.C) ^
        ScalarMultiply(second.C, third.B);
    const uint8_t di_fg = ScalarMultiply(second.A, third.C) ^
        ScalarMultiply(second.C, third.A);
    const uint8_t dh_eg = ScalarMultiply(second.A, third.B) ^
        ScalarMultiply(second.B, third.A);
    return ScalarMultiply(first.A, ei_fh) ^
        ScalarMultiply(first.B, di_fg) ^
        ScalarMultiply(first.C, dh_eg);
}

K3HyperovalDirection ToCodecDirection(const Direction& direction) noexcept
{
    K3HyperovalDirection converted;
    converted.Point = DirectionKey(direction);
    converted.Alpha = direction.A;
    converted.Beta = direction.B;
    converted.Gamma = direction.C;
    return converted;
}

uint32_t Avalanche32(uint32_t value) noexcept
{
    value = (value ^ (value >> 16u)) * UINT32_C(0x7feb352d);
    value = (value ^ (value >> 15u)) * UINT32_C(0x846ca68b);
    return value ^ (value >> 16u);
}

uint32_t ComplementBucketForId(uint32_t id) noexcept
{
    return static_cast<uint32_t>(
        (static_cast<uint64_t>(Avalanche32(id)) * UINT64_C(65535)) >> 32u);
}

Direction ScalarComplementDirection(uint32_t bucket) noexcept
{
    if (bucket < 65280u)
    {
        const uint8_t a = static_cast<uint8_t>(bucket / 255u);
        const uint8_t d = static_cast<uint8_t>(1u + bucket % 255u);
        return NormalizeDirection(
            1u, a, ScalarMultiply(a, a) ^ d);
    }
    if (bucket < 65535u)
    {
        const uint8_t a = static_cast<uint8_t>(
            1u + bucket - 65280u);
        return NormalizeDirection(0u, 1u, a);
    }
    return Direction();
}

Direction ScalarDirectionForId(uint32_t id) noexcept
{
    if (id == 0u) return NormalizeDirection(1u, 0u, 0u);
    if (id == 1u) return NormalizeDirection(0u, 1u, 0u);
    if (id == 2u) return NormalizeDirection(0u, 0u, 1u);
    if (id <= 257u)
    {
        const uint8_t t = static_cast<uint8_t>(id - 2u);
        return NormalizeDirection(1u, t, ScalarMultiply(t, t));
    }
    return ScalarComplementDirection(ComplementBucketForId(id));
}

Direction CodecDirectionForId(uint32_t id) noexcept
{
    const K3HyperovalDirection value =
        wirehair_wh2_bench::K3HyperovalDirectionForPacketId(id);
    return NormalizeDirection(value.Alpha, value.Beta, value.Gamma);
}

bool CrossCheckDeterminant(
    const Direction& first,
    const Direction& second,
    const Direction& third) noexcept
{
    const uint8_t scalar = ScalarDeterminant(first, second, third);
    const uint8_t codec = wirehair_wh2_bench::K3HyperovalDeterminant(
        ToCodecDirection(first),
        ToCodecDirection(second),
        ToCodecDirection(third));
    return scalar == codec;
}

void AppendLe32(std::string& bytes, uint32_t value)
{
    for (unsigned shift = 0u; shift != 32u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendLe64(std::string& bytes, uint64_t value)
{
    for (unsigned shift = 0u; shift != 64u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendTaggedString(std::string& bytes, const char* value)
{
    bytes.append(value);
    bytes.push_back('\0');
}

class SplitMix64
{
public:
    explicit SplitMix64(uint64_t state) : State(state) {}

    uint64_t Next() noexcept
    {
        uint64_t value =
            (State += UINT64_C(0x9e3779b97f4a7c15));
        value = (value ^ (value >> 30u)) *
            UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27u)) *
            UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31u);
    }

private:
    uint64_t State;
};

bool CheckedId(uint64_t value, uint32_t& id_out) noexcept
{
    if (value > UINT32_MAX) return false;
    id_out = static_cast<uint32_t>(value);
    return true;
}

bool CheckedMultiply64(
    uint64_t left,
    uint64_t right,
    uint64_t& product) noexcept
{
    if (left != 0u && right > UINT64_MAX / left) return false;
    product = left * right;
    return true;
}

bool CheckedAdd64(
    uint64_t left,
    uint64_t right,
    uint64_t& sum) noexcept
{
    if (right > UINT64_MAX - left) return false;
    sum = left + right;
    return true;
}

struct AnalyticReceipt
{
    uint64_t SecantLines = 0u;
    uint64_t ExternalLines = 0u;
    uint64_t DuplicateNumerator = 0u;
    uint64_t DuplicateDenominator = 0u;
    uint64_t DistinctCollinearNumerator = 0u;
    uint64_t DistinctCollinearDenominator = 0u;
    uint64_t TotalNumerator = 0u;
    uint64_t TotalDenominator = 0u;
    uint64_t TwoLowNumerator = 0u;
    uint64_t TwoLowDenominator = 0u;
};

bool BuildAnalyticReceipt(AnalyticReceipt& receipt, std::string& error)
{
    receipt = AnalyticReceipt();
    const uint64_t q = 256u;
    uint64_t q_squared = 0u;
    uint64_t projective_lines = 0u;
    uint64_t temporary = 0u;
    if (!CheckedMultiply64(q, q, q_squared) ||
        !CheckedAdd64(q_squared, q, temporary) ||
        !CheckedAdd64(temporary, 1u, projective_lines))
    {
        error = "analytic projective-line arithmetic overflow";
        return false;
    }

    uint64_t secants = 0u;
    if (!CheckedMultiply64(258u / 2u, 257u, secants) ||
        secants > projective_lines)
    {
        error = "analytic secant arithmetic overflow";
        return false;
    }
    receipt.SecantLines = secants;
    receipt.ExternalLines = projective_lines - secants;

    const uint64_t complement_points = 65535u;
    if (!CheckedMultiply64(
            complement_points,
            complement_points,
            receipt.DuplicateDenominator) ||
        !CheckedMultiply64(3u, complement_points, temporary) || temporary < 2u)
    {
        error = "analytic probability arithmetic overflow";
        return false;
    }
    receipt.DuplicateNumerator = temporary - 2u;
    receipt.DistinctCollinearNumerator = 5548546u;
    receipt.DistinctCollinearDenominator =
        receipt.DuplicateDenominator / 3u;
    uint64_t collinear_scaled = 0u;
    if (receipt.DuplicateDenominator % 3u != 0u ||
        !CheckedMultiply64(
            3u,
            receipt.DistinctCollinearNumerator,
            collinear_scaled) ||
        !CheckedAdd64(
            receipt.DuplicateNumerator,
            collinear_scaled,
            receipt.TotalNumerator))
    {
        error = "analytic fraction combination overflow";
        return false;
    }
    receipt.TotalDenominator = receipt.DuplicateDenominator;
    receipt.TwoLowNumerator = 255u;
    receipt.TwoLowDenominator = complement_points;

    if (receipt.SecantLines != 33153u ||
        receipt.ExternalLines != 32640u ||
        receipt.DuplicateNumerator != 196603u ||
        receipt.DuplicateDenominator != UINT64_C(4294836225) ||
        receipt.DistinctCollinearNumerator != 5548546u ||
        receipt.DistinctCollinearDenominator != UINT64_C(1431612075) ||
        receipt.TotalNumerator != 16842241u ||
        receipt.TotalDenominator != UINT64_C(4294836225) ||
        receipt.TwoLowNumerator != 255u ||
        receipt.TwoLowDenominator != 65535u)
    {
        error = "analytic frozen identity mismatch";
        return false;
    }
    return true;
}

std::vector<std::pair<uint16_t, uint16_t> > LowPairs()
{
    std::vector<std::pair<uint16_t, uint16_t> > pairs;
    pairs.reserve(33153u);
    for (uint16_t first = 0u; first < 258u; ++first) {
        for (uint16_t second = static_cast<uint16_t>(first + 1u);
             second < 258u;
             ++second)
        {
            pairs.push_back(std::make_pair(first, second));
        }
    }
    return pairs;
}

bool GenerateFamily(
    Family family,
    uint64_t count,
    std::vector<Triple>& triples_out,
    std::string& error)
{
    triples_out.clear();
    error.clear();
    if (count > static_cast<uint64_t>(
            std::numeric_limits<size_t>::max()))
    {
        error = "corpus count exceeds size_t";
        return false;
    }
    try {
        triples_out.resize(static_cast<size_t>(count));
    }
    catch (const std::bad_alloc&) {
        error = "corpus allocation failed";
        return false;
    }
    catch (const std::length_error&) {
        error = "corpus length failed";
        return false;
    }

    SplitMix64 random(kHighIidRoot);
    const std::vector<std::pair<uint16_t, uint16_t> > pairs =
        family == Family::TwoLowOneHigh ? LowPairs() :
        std::vector<std::pair<uint16_t, uint16_t> >();
    if (family == Family::TwoLowOneHigh && pairs.size() != 33153u)
    {
        error = "low-pair cardinality mismatch";
        return false;
    }

    for (uint64_t i = 0u; i < count; ++i)
    {
        Triple triple;
        if (family == Family::HighIid)
        {
            for (unsigned position = 0u; position != 3u; ++position)
            {
                for (;;)
                {
                    const uint32_t id = static_cast<uint32_t>(random.Next());
                    if (id < 258u) continue;
                    bool duplicate = false;
                    for (unsigned prior = 0u; prior < position; ++prior) {
                        duplicate = duplicate || triple.Id[prior] == id;
                    }
                    if (!duplicate) {
                        triple.Id[position] = id;
                        break;
                    }
                }
            }
        }
        else if (family == Family::Sequential)
        {
            const uint64_t base = UINT64_C(258) + UINT64_C(3) * i;
            if (!CheckedId(base, triple.Id[0]) ||
                !CheckedId(base + 1u, triple.Id[1]) ||
                !CheckedId(base + 2u, triple.Id[2]))
            {
                error = "sequential corpus wrapped";
                return false;
            }
        }
        else if (family == Family::Delta257 ||
                 family == Family::Delta65535)
        {
            const uint64_t delta = family == Family::Delta257 ?
                UINT64_C(257) : UINT64_C(65535);
            const uint64_t quotient = i / delta;
            const uint64_t remainder = i % delta;
            const uint64_t base = UINT64_C(258) + remainder +
                UINT64_C(3) * quotient * delta;
            if (!CheckedId(base, triple.Id[0]) ||
                !CheckedId(base + delta, triple.Id[1]) ||
                !CheckedId(base + 2u * delta, triple.Id[2]))
            {
                error = "delta corpus wrapped";
                return false;
            }
        }
        else if (family == Family::Descending)
        {
            const uint64_t offset = UINT64_C(3) * i;
            if (offset + 2u > UINT32_MAX - UINT64_C(258))
            {
                error = "descending corpus wrapped";
                return false;
            }
            triple.Id[0] = UINT32_MAX - static_cast<uint32_t>(offset);
            triple.Id[1] = UINT32_MAX - static_cast<uint32_t>(offset + 1u);
            triple.Id[2] = UINT32_MAX - static_cast<uint32_t>(offset + 2u);
        }
        else if (family == Family::OneLowTwoHigh)
        {
            const uint64_t first_high = UINT64_C(258) + UINT64_C(2) * i;
            triple.Id[0] = static_cast<uint32_t>(i % UINT64_C(258));
            if (!CheckedId(first_high, triple.Id[1]) ||
                !CheckedId(first_high + 1u, triple.Id[2]))
            {
                error = "one-low corpus wrapped";
                return false;
            }
        }
        else if (family == Family::TwoLowOneHigh)
        {
            const std::pair<uint16_t, uint16_t>& pair =
                pairs[static_cast<size_t>(i % pairs.size())];
            triple.Id[0] = pair.first;
            triple.Id[1] = pair.second;
            if (!CheckedId(UINT64_C(258) + i, triple.Id[2]))
            {
                error = "two-low corpus wrapped";
                return false;
            }
        }
        else
        {
            error = "invalid corpus family";
            return false;
        }

        if (triple.Id[0] == triple.Id[1] ||
            triple.Id[0] == triple.Id[2] ||
            triple.Id[1] == triple.Id[2])
        {
            error = "corpus triple contains duplicate raw IDs";
            return false;
        }
        triples_out[static_cast<size_t>(i)] = triple;
    }
    return true;
}

std::string CorpusHash(
    const char* family_name,
    uint64_t declared_count,
    uint64_t root,
    uint32_t config,
    const std::vector<Triple>& triples)
{
    std::string bytes;
    const uint64_t prefix_bytes = static_cast<uint64_t>(
        std::strlen(kCorpusSchema) + std::strlen(family_name) + 2u +
        8u + 8u + 4u + 4u);
    const uint64_t total = prefix_bytes +
        static_cast<uint64_t>(triples.size()) * 12u;
    if (total <= static_cast<uint64_t>(
            std::numeric_limits<size_t>::max()))
    {
        bytes.reserve(static_cast<size_t>(total));
    }
    AppendTaggedString(bytes, kCorpusSchema);
    AppendTaggedString(bytes, family_name);
    AppendLe64(bytes, declared_count);
    AppendLe64(bytes, root);
    AppendLe32(bytes, config);
    AppendLe32(bytes, 0u);
    for (const Triple& triple : triples) {
        AppendLe32(bytes, triple.Id[0]);
        AppendLe32(bytes, triple.Id[1]);
        AppendLe32(bytes, triple.Id[2]);
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

bool GenerateLowCorpus(
    std::vector<Triple>* triples_out,
    std::string& hash_out,
    std::string& error)
{
    if (triples_out) triples_out->clear();
    std::vector<Triple> local;
    try {
        local.reserve(static_cast<size_t>(kLowTriples));
        for (uint32_t first = 0u; first < 258u; ++first) {
            for (uint32_t second = first + 1u; second < 258u; ++second) {
                for (uint32_t third = second + 1u; third < 258u; ++third) {
                    Triple triple;
                    triple.Id[0] = first;
                    triple.Id[1] = second;
                    triple.Id[2] = third;
                    local.push_back(triple);
                }
            }
        }
    }
    catch (const std::bad_alloc&) {
        error = "low corpus allocation failed";
        return false;
    }
    catch (const std::length_error&) {
        error = "low corpus length failed";
        return false;
    }
    if (local.size() != static_cast<size_t>(kLowTriples)) {
        error = "low corpus cardinality mismatch";
        return false;
    }
    hash_out = CorpusHash("low-only", kLowTriples, 0u, 258u, local);
    if (hash_out.size() != 64u) {
        error = "low corpus hash failed";
        return false;
    }
    if (triples_out) triples_out->swap(local);
    return true;
}

uint32_t FamilyConfig(Family family)
{
    if (family == Family::Delta257) return 257u;
    if (family == Family::Delta65535) return 65535u;
    if (family == Family::OneLowTwoHigh) return 258u;
    if (family == Family::TwoLowOneHigh) return 33153u;
    return 0u;
}

uint64_t FamilyRoot(Family family)
{
    return family == Family::HighIid ? kHighIidRoot : 0u;
}

struct FrozenEntry
{
    wirehair::wh2_benchmark::FrozenRecoveryCell Cell;
    std::vector<uint32_t> Ids;
};

struct FrozenManifest
{
    std::string CanonicalHash;
    std::string TraceHash;
    std::vector<FrozenEntry> Entries;
};

bool BuildFrozenManifest(FrozenManifest& manifest, std::string& error)
{
    manifest = FrozenManifest();
    if (wirehair::wh2_benchmark::DevelopmentRecoveryDomainSha256() !=
            kDevelopmentDomainSha256)
    {
        error = "development recovery domain hash mismatch";
        return false;
    }
    const std::vector<wirehair::wh2_benchmark::FrozenRecoveryCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentRecoveryCells();
    std::string canonical;
    std::string trace_bytes;
    for (const wirehair::wh2_benchmark::FrozenRecoveryCell& cell : cells)
    {
        if (cell.K != kK) continue;
        wirehair::wh2_benchmark::FrozenPacketTrace trace;
        if (wirehair::wh2_benchmark::GenerateFrozenPacketTrace(cell, trace) !=
                wirehair::wh2_benchmark::FrozenTraceStatus::Complete ||
            trace.delivered_ids.size() != 7u ||
            trace.trace_sha256 !=
                wirehair::wh2_benchmark::PacketIdsSha256(trace.delivered_ids))
        {
            error = "frozen K3 trace generation failed";
            return false;
        }
        canonical += wirehair::wh2_benchmark::CanonicalRecoveryCellJson(cell);
        canonical.push_back('\n');
        AppendLe64(trace_bytes, static_cast<uint64_t>(cell.ordinal));
        for (uint32_t id : trace.delivered_ids) AppendLe32(trace_bytes, id);
        FrozenEntry entry;
        entry.Cell = cell;
        entry.Ids.swap(trace.delivered_ids);
        manifest.Entries.push_back(std::move(entry));
    }
    if (manifest.Entries.size() != 12u) {
        error = "frozen K3 cell count is not 12";
        return false;
    }
    manifest.CanonicalHash =
        wirehair::wh2_benchmark::Sha256Hex(canonical);
    manifest.TraceHash =
        wirehair::wh2_benchmark::Sha256Hex(trace_bytes);
    return manifest.CanonicalHash.size() == 64u &&
        manifest.TraceHash.size() == 64u;
}

bool ComplementManifestHash(std::string& hash_out, std::string& error)
{
    std::string bytes;
    AppendTaggedString(bytes, kComplementSchema);
    AppendLe32(bytes, 65535u);
    for (uint32_t bucket = 0u; bucket < 65535u; ++bucket)
    {
        const Direction direction = ScalarComplementDirection(bucket);
        if (!direction.Nonzero) {
            error = "scalar complement produced zero";
            return false;
        }
        AppendLe32(bytes, bucket);
        bytes.push_back(static_cast<char>(direction.A));
        bytes.push_back(static_cast<char>(direction.B));
        bytes.push_back(static_cast<char>(direction.C));
    }
    hash_out = wirehair::wh2_benchmark::Sha256Hex(bytes);
    return hash_out.size() == 64u;
}

struct ManifestData
{
    std::string LowHash;
    std::array<std::string, static_cast<size_t>(Family::Count)> FamilyHashes;
    std::array<std::vector<Triple>, static_cast<size_t>(Family::Count)> Corpora;
    FrozenManifest Frozen;
    std::string ComplementHash;
};

bool BuildManifestData(
    bool retain_corpora,
    ManifestData& manifest,
    std::string& error)
{
    manifest = ManifestData();
    if (!GenerateLowCorpus(nullptr, manifest.LowHash, error)) return false;
    for (uint32_t family_value = 0u;
         family_value < static_cast<uint32_t>(Family::Count);
         ++family_value)
    {
        const Family family = static_cast<Family>(family_value);
        std::vector<Triple> triples;
        if (!GenerateFamily(family, kStochasticTriples, triples, error)) {
            return false;
        }
        manifest.FamilyHashes[family_value] = CorpusHash(
            FamilyName(family),
            kStochasticTriples,
            FamilyRoot(family),
            FamilyConfig(family),
            triples);
        if (manifest.FamilyHashes[family_value].size() != 64u) {
            error = "stochastic corpus hash failed";
            return false;
        }
        if (retain_corpora) {
            manifest.Corpora[family_value].swap(triples);
        }
    }
    return BuildFrozenManifest(manifest.Frozen, error) &&
        ComplementManifestHash(manifest.ComplementHash, error);
}

const char* ExpectedFamilyHash(Family family)
{
    switch (family)
    {
    case Family::HighIid:       return WIREHAIR_WH2_K3_HIGH_IID_SHA256;
    case Family::Sequential:    return WIREHAIR_WH2_K3_SEQUENTIAL_SHA256;
    case Family::Delta257:      return WIREHAIR_WH2_K3_DELTA257_SHA256;
    case Family::Delta65535:    return WIREHAIR_WH2_K3_DELTA65535_SHA256;
    case Family::Descending:    return WIREHAIR_WH2_K3_DESCENDING_SHA256;
    case Family::OneLowTwoHigh: return WIREHAIR_WH2_K3_ONE_LOW_SHA256;
    case Family::TwoLowOneHigh: return WIREHAIR_WH2_K3_TWO_LOW_SHA256;
    case Family::Count: break;
    }
    return "";
}

bool ExpectedManifestConfigured()
{
    if (std::strlen(WIREHAIR_WH2_K3_LOW_CORPUS_SHA256) != 64u ||
        std::strlen(WIREHAIR_WH2_K3_CANONICAL_SHA256) != 64u ||
        std::strlen(WIREHAIR_WH2_K3_TRACE_SHA256) != 64u ||
        std::strlen(WIREHAIR_WH2_K3_COMPLEMENT_SHA256) != 64u)
    {
        return false;
    }
    for (uint32_t value = 0u;
         value < static_cast<uint32_t>(Family::Count);
         ++value)
    {
        if (std::strlen(ExpectedFamilyHash(static_cast<Family>(value))) != 64u)
            return false;
    }
    return true;
}

bool AuthenticateManifest(const ManifestData& manifest, std::string& error)
{
    if (!ExpectedManifestConfigured()) {
        error = "manifest constants are not frozen into this executable";
        return false;
    }
    if (manifest.LowHash != WIREHAIR_WH2_K3_LOW_CORPUS_SHA256) {
        error = "low corpus hash mismatch";
        return false;
    }
    for (uint32_t value = 0u;
         value < static_cast<uint32_t>(Family::Count);
         ++value)
    {
        if (manifest.FamilyHashes[value] !=
            ExpectedFamilyHash(static_cast<Family>(value)))
        {
            error = std::string(FamilyName(static_cast<Family>(value))) +
                " corpus hash mismatch";
            return false;
        }
    }
    if (manifest.Frozen.CanonicalHash != WIREHAIR_WH2_K3_CANONICAL_SHA256 ||
        manifest.Frozen.TraceHash != WIREHAIR_WH2_K3_TRACE_SHA256 ||
        manifest.ComplementHash != WIREHAIR_WH2_K3_COMPLEMENT_SHA256)
    {
        error = "frozen/complement manifest hash mismatch";
        return false;
    }
    return true;
}

int PublishTerminal(const char* line, int status)
{
    std::cout << line << '\n';
    std::cout.flush();
    return std::cout.good() ? status : 1;
}

int RunManifest()
{
    ManifestData manifest;
    std::string error;
    if (!BuildManifestData(false, manifest, error) ||
        (ExpectedManifestConfigured() &&
         !AuthenticateManifest(manifest, error)))
    {
        std::cerr << "K3 manifest invalid: " << error << '\n';
        return PublishTerminal("K3_MANIFEST=INVALID", 1);
    }
    std::cout << "MANIFEST,protocol=" << kProtocolTag
              << ",kind=low-only,n=" << kLowTriples
              << ",sha256=" << manifest.LowHash << '\n';
    for (uint32_t value = 0u;
         value < static_cast<uint32_t>(Family::Count);
         ++value)
    {
        const Family family = static_cast<Family>(value);
        std::cout << "MANIFEST,protocol=" << kProtocolTag
                  << ",kind=" << FamilyName(family)
                  << ",n=" << kStochasticTriples
                  << ",root=0x" << std::hex << FamilyRoot(family)
                  << std::dec << ",config=" << FamilyConfig(family)
                  << ",sha256=" << manifest.FamilyHashes[value] << '\n';
    }
    std::cout << "MANIFEST,protocol=" << kProtocolTag
              << ",kind=frozen-canonical,cells=12,sha256="
              << manifest.Frozen.CanonicalHash << '\n'
              << "MANIFEST,protocol=" << kProtocolTag
              << ",kind=frozen-trace,cells=12,ids_per_cell=7,sha256="
              << manifest.Frozen.TraceHash << '\n'
              << "MANIFEST,protocol=" << kProtocolTag
              << ",kind=complement-points,n=65535,sha256="
              << manifest.ComplementHash << '\n';
    return PublishTerminal("K3_MANIFEST=COMPLETE", 0);
}

struct LegacyOwner
{
    WirehairCodec Codec = nullptr;
    ~LegacyOwner() { wirehair_free(Codec); }
    LegacyOwner() = default;
    LegacyOwner(const LegacyOwner&) = delete;
    LegacyOwner& operator=(const LegacyOwner&) = delete;
};

struct V2Owner
{
    WirehairV2Codec Codec = nullptr;
    ~V2Owner() { wirehair_v2_free(Codec); }
    V2Owner() = default;
    V2Owner(const V2Owner&) = delete;
    V2Owner& operator=(const V2Owner&) = delete;
};

struct PublicProfiles
{
    WirehairWireProfile Legacy = {};
    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> V2 = {};
    std::string V2Hash;
    uint8_t V2Attempt = 0u;
};

std::vector<uint8_t> IdentityLaneMessage(uint32_t block_bytes)
{
    std::vector<uint8_t> source(static_cast<size_t>(kK) * block_bytes, 0u);
    if (block_bytes >= 4u)
    {
        source[0u * block_bytes + 0u] = 1u;
        source[1u * block_bytes + 1u] = 1u;
        source[2u * block_bytes + 2u] = 1u;
        source[0u * block_bytes + 3u] = 1u;
        source[1u * block_bytes + 3u] = 1u;
        source[2u * block_bytes + 3u] = 1u;
    }
    return source;
}

bool InitializeLegacyEncoder(
    const WirehairWireProfile& profile,
    const std::vector<uint8_t>& source,
    uint32_t block_bytes,
    LegacyOwner& owner,
    std::string& error)
{
    if (owner.Codec) {
        error = "legacy owner was not empty";
        return false;
    }
    const WirehairResult result = wirehair_encoder_create_profile_ex(
        nullptr,
        source.data(),
        source.size(),
        block_bytes,
        &profile,
        WIREHAIR_ENCODER_OWN_INPUT,
        &owner.Codec);
    if (result != Wirehair_Success || !owner.Codec) {
        std::ostringstream stream;
        stream << "Wirehair1 encoder construction failed: result="
               << static_cast<uint32_t>(result);
        error = stream.str();
        return false;
    }
    return true;
}

bool InitializeV2EncoderFromProfile(
    const std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES>& profile,
    const std::vector<uint8_t>& source,
    V2Owner& owner,
    std::string& error)
{
    if (owner.Codec) {
        error = "V2 owner was not empty";
        return false;
    }
    const WirehairV2Result result = wirehair_v2_encoder_create_profile(
        source.data(), profile.data(), profile.size(), &owner.Codec);
    if (result != WirehairV2_Success || !owner.Codec) {
        std::ostringstream stream;
        stream << "WH2 exact-profile encoder construction failed: result="
               << static_cast<uint32_t>(result);
        error = stream.str();
        return false;
    }
    return true;
}

bool SelectPublicProfiles(
    const std::vector<uint8_t>& identity_source,
    PublicProfiles& profiles,
    V2Owner& selected_encoder,
    std::string& error)
{
    if (wirehair_wire_profile_init(
            WIREHAIR_LEGACY_PROFILE_CURRENT, &profiles.Legacy) !=
            Wirehair_Success ||
        profiles.Legacy.profile_id != WIREHAIR_LEGACY_PROFILE_CURRENT)
    {
        error = "Wirehair1 current profile initialization failed";
        return false;
    }

    uint32_t profile_bytes = 0u;
    const WirehairV2Result result = wirehair_v2_encoder_create_profile_id(
        WIREHAIR_V2_PROFILE_CURRENT,
        identity_source.data(),
        kComparatorMessageBytes,
        kComparatorBlockBytes,
        profiles.V2.data(),
        profiles.V2.size(),
        &profile_bytes,
        &selected_encoder.Codec);
    if (result != WirehairV2_Success || !selected_encoder.Codec ||
        profile_bytes != profiles.V2.size())
    {
        std::ostringstream stream;
        stream << "WH2 public profile selection failed: result="
               << static_cast<uint32_t>(result)
               << " bytes=" << profile_bytes;
        error = stream.str();
        return false;
    }
    WirehairV2Profile decoded = {};
    if (wirehair_v2_profile_deserialize(
            profiles.V2.data(), profiles.V2.size(), &decoded) !=
            WirehairV2_Success ||
        decoded.profile_id != WIREHAIR_V2_PROFILE_CURRENT ||
        decoded.message_bytes != kComparatorMessageBytes ||
        decoded.block_bytes != kComparatorBlockBytes)
    {
        error = "WH2 selected profile receipt mismatch";
        return false;
    }
    profiles.V2Attempt = decoded.seed_attempt;
    profiles.V2Hash = wirehair::wh2_benchmark::Sha256Hex(
        profiles.V2.data(), profiles.V2.size());
    return profiles.V2Hash.size() == 64u;
}

bool VerifyIdentityPacket(
    const uint8_t* packet,
    uint32_t bytes,
    Direction& direction,
    std::string& error,
    bool normalize = true)
{
    if (!packet || bytes != kComparatorBlockBytes) {
        error = "identity packet length mismatch";
        return false;
    }
    if (packet[3] != static_cast<uint8_t>(
            packet[0] ^ packet[1] ^ packet[2]))
    {
        error = "identity packet checksum lane mismatch";
        return false;
    }
    uint8_t zero_or = 0u;
    for (uint32_t i = 4u; i < bytes; ++i) zero_or |= packet[i];
    if (zero_or != 0u) {
        error = "identity packet zero lane mismatch";
        return false;
    }
    if (normalize) {
        direction = NormalizeDirection(packet[0], packet[1], packet[2]);
    } else {
        direction.A = packet[0];
        direction.B = packet[1];
        direction.C = packet[2];
        direction.Nonzero =
            (packet[0] | packet[1] | packet[2]) != 0u;
    }
    // A public comparator is allowed to expose a zero equation row.  It is a
    // rank failure, not an extraction failure, and is therefore retained for
    // determinant classification and operational concordance.
    return true;
}

bool ExtractLegacyDirection(
    WirehairCodec codec,
    uint32_t id,
    Direction& direction,
    std::string& error,
    bool normalize = true)
{
    std::array<uint8_t, kComparatorBlockBytes> packet = {};
    uint32_t written = 0u;
    const WirehairResult result = wirehair_encode(
        codec, id, packet.data(), packet.size(), &written);
    if (result != Wirehair_Success || written != packet.size()) {
        std::ostringstream stream;
        stream << "Wirehair1 identity encode failed for id=" << id
               << " result=" << static_cast<uint32_t>(result)
               << " bytes=" << written;
        error = stream.str();
        return false;
    }
    return VerifyIdentityPacket(
        packet.data(), written, direction, error, normalize);
}

bool ExtractV2Direction(
    WirehairV2Codec codec,
    uint32_t id,
    Direction& direction,
    std::string& error,
    bool normalize = true)
{
    std::array<uint8_t, kComparatorBlockBytes> packet = {};
    uint32_t written = 0u;
    const WirehairV2Result result = wirehair_v2_encode(
        codec, id, packet.data(), packet.size(), &written);
    if (result != WirehairV2_Success || written != packet.size()) {
        std::ostringstream stream;
        stream << "WH2 identity encode failed for id=" << id
               << " result=" << static_cast<uint32_t>(result)
               << " bytes=" << written;
        error = stream.str();
        return false;
    }
    return VerifyIdentityPacket(
        packet.data(), written, direction, error, normalize);
}

struct WilsonInterval
{
    long double Lower = 0.0L;
    long double Upper = 0.0L;
};

WilsonInterval Wilson(uint64_t failures, uint64_t count)
{
    WilsonInterval interval;
    if (count == 0u || failures > count) return interval;
    const long double n = static_cast<long double>(count);
    const long double p = static_cast<long double>(failures) / n;
    const long double z2 = kWilsonZ * kWilsonZ;
    const long double denominator = 1.0L + z2 / n;
    const long double center = (p + z2 / (2.0L * n)) / denominator;
    const long double radius = kWilsonZ / denominator * std::sqrt(
        p * (1.0L - p) / n + z2 / (4.0L * n * n));
    interval.Lower = std::max(0.0L, center - radius);
    interval.Upper = std::min(1.0L, center + radius);
    return interval;
}

long double ExactMcNemarUpperTail(
    uint64_t candidate_only,
    uint64_t comparator_only)
{
    const uint64_t discordant = candidate_only + comparator_only;
    if (discordant == 0u || comparator_only <= candidate_only) return 1.0L;
    const uint64_t start = comparator_only;
    const long double log_term =
        std::lgammal(static_cast<long double>(discordant) + 1.0L) -
        std::lgammal(static_cast<long double>(start) + 1.0L) -
        std::lgammal(static_cast<long double>(discordant - start) + 1.0L) -
        static_cast<long double>(discordant) * std::log(2.0L);
    long double term = std::exp(log_term);
    if (term == 0.0L) return 0.0L;
    long double sum = term;
    unsigned negligible = 0u;
    for (uint64_t k = start; k < discordant; ++k)
    {
        term *= static_cast<long double>(discordant - k) /
            static_cast<long double>(k + 1u);
        sum += term;
        if (term <= sum * 1e-20L) {
            if (++negligible == 32u) break;
        }
        else {
            negligible = 0u;
        }
    }
    return std::min(1.0L, sum);
}

struct ComparisonCounts
{
    uint64_t Count = 0u;
    uint64_t CandidateFailures = 0u;
    uint64_t LegacyFailures = 0u;
    uint64_t V2Failures = 0u;
    uint64_t CandidateOnlyLegacy = 0u;
    uint64_t LegacyOnlyCandidate = 0u;
    uint64_t CandidateOnlyV2 = 0u;
    uint64_t V2OnlyCandidate = 0u;
    uint64_t LegacyZeroRowOccurrences = 0u;
    uint64_t V2ZeroRowOccurrences = 0u;

    void Add(const ComparisonCounts& other)
    {
        Count += other.Count;
        CandidateFailures += other.CandidateFailures;
        LegacyFailures += other.LegacyFailures;
        V2Failures += other.V2Failures;
        CandidateOnlyLegacy += other.CandidateOnlyLegacy;
        LegacyOnlyCandidate += other.LegacyOnlyCandidate;
        CandidateOnlyV2 += other.CandidateOnlyV2;
        V2OnlyCandidate += other.V2OnlyCandidate;
        LegacyZeroRowOccurrences += other.LegacyZeroRowOccurrences;
        V2ZeroRowOccurrences += other.V2ZeroRowOccurrences;
    }
};

void RecordOutcomes(
    bool candidate_failure,
    bool legacy_failure,
    bool v2_failure,
    ComparisonCounts& counts)
{
    ++counts.Count;
    counts.CandidateFailures += candidate_failure ? 1u : 0u;
    counts.LegacyFailures += legacy_failure ? 1u : 0u;
    counts.V2Failures += v2_failure ? 1u : 0u;
    counts.CandidateOnlyLegacy +=
        candidate_failure && !legacy_failure ? 1u : 0u;
    counts.LegacyOnlyCandidate +=
        !candidate_failure && legacy_failure ? 1u : 0u;
    counts.CandidateOnlyV2 +=
        candidate_failure && !v2_failure ? 1u : 0u;
    counts.V2OnlyCandidate +=
        !candidate_failure && v2_failure ? 1u : 0u;
}

bool DirectionsForTriple(
    const Triple& triple,
    WirehairCodec legacy,
    WirehairV2Codec v2,
    Direction candidate[3],
    Direction legacy_directions[3],
    Direction v2_directions[3],
    std::string& error)
{
    for (unsigned i = 0u; i != 3u; ++i)
    {
        candidate[i] = ScalarDirectionForId(triple.Id[i]);
        if (!SameDirection(candidate[i], CodecDirectionForId(triple.Id[i]))) {
            error = "candidate scalar/codec direction mismatch";
            return false;
        }
        if (!ExtractLegacyDirection(
                legacy, triple.Id[i], legacy_directions[i], error) ||
            !ExtractV2Direction(v2, triple.Id[i], v2_directions[i], error))
        {
            return false;
        }
    }
    if (!CrossCheckDeterminant(candidate[0], candidate[1], candidate[2]) ||
        !CrossCheckDeterminant(
            legacy_directions[0], legacy_directions[1], legacy_directions[2]) ||
        !CrossCheckDeterminant(
            v2_directions[0], v2_directions[1], v2_directions[2]))
    {
        error = "table/scalar determinant mismatch";
        return false;
    }
    return true;
}

struct WorkerResult
{
    ComparisonCounts Counts;
    std::string Error;
};

void ClassifyWorker(
    const std::vector<Triple>* corpus,
    size_t begin,
    size_t end,
    const PublicProfiles* profiles,
    const std::vector<uint8_t>* identity_source,
    std::atomic<bool>* stop,
    WorkerResult* result)
{
    try
    {
        LegacyOwner legacy;
        V2Owner v2;
        if (!InitializeLegacyEncoder(
                profiles->Legacy,
                *identity_source,
                kComparatorBlockBytes,
                legacy,
                result->Error) ||
            !InitializeV2EncoderFromProfile(
                profiles->V2, *identity_source, v2, result->Error))
        {
            stop->store(true, std::memory_order_release);
            return;
        }
        for (size_t index = begin;
             index < end && !stop->load(std::memory_order_acquire);
             ++index)
        {
            Direction candidate[3];
            Direction legacy_directions[3];
            Direction v2_directions[3];
            if (!DirectionsForTriple(
                    (*corpus)[index],
                    legacy.Codec,
                    v2.Codec,
                    candidate,
                    legacy_directions,
                    v2_directions,
                    result->Error))
            {
                stop->store(true, std::memory_order_release);
                return;
            }
            for (unsigned row = 0u; row != 3u; ++row)
            {
                result->Counts.LegacyZeroRowOccurrences +=
                    legacy_directions[row].Nonzero ? 0u : 1u;
                result->Counts.V2ZeroRowOccurrences +=
                    v2_directions[row].Nonzero ? 0u : 1u;
            }
            RecordOutcomes(
                ScalarDeterminant(
                    candidate[0], candidate[1], candidate[2]) == 0u,
                ScalarDeterminant(
                    legacy_directions[0], legacy_directions[1],
                    legacy_directions[2]) == 0u,
                ScalarDeterminant(
                    v2_directions[0], v2_directions[1],
                    v2_directions[2]) == 0u,
                result->Counts);
        }
    }
    catch (...)
    {
        try {
            result->Error = "worker caught an exception";
        }
        catch (...) {
        }
        stop->store(true, std::memory_order_release);
    }
}

bool ClassifyCorpus(
    const std::vector<Triple>& corpus,
    const PublicProfiles& profiles,
    const std::vector<uint8_t>& identity_source,
    ComparisonCounts& counts,
    std::string& error)
{
    counts = ComparisonCounts();
    std::atomic<bool> stop(false);
    std::array<WorkerResult, kWorkerCount> results;
    std::vector<std::thread> threads;
    try {
        threads.reserve(kWorkerCount);
        for (uint32_t worker = 0u; worker < kWorkerCount; ++worker)
        {
            const size_t begin = corpus.size() * worker / kWorkerCount;
            const size_t end = corpus.size() * (worker + 1u) / kWorkerCount;
            threads.push_back(std::thread(
                ClassifyWorker,
                &corpus,
                begin,
                end,
                &profiles,
                &identity_source,
                &stop,
                &results[worker]));
        }
    }
    catch (const std::system_error& exception) {
        stop.store(true, std::memory_order_release);
        for (std::thread& thread : threads) if (thread.joinable()) thread.join();
        error = std::string("worker creation failed: ") + exception.what();
        return false;
    }
    for (std::thread& thread : threads) thread.join();
    for (const WorkerResult& result : results)
    {
        if (!result.Error.empty()) {
            error = result.Error;
            return false;
        }
        counts.Add(result.Counts);
    }
    if (counts.Count != corpus.size()) {
        error = "worker classified count mismatch";
        return false;
    }
    return true;
}

bool ValidateComplementAndMapper(std::string& error)
{
    std::vector<uint8_t> seen(UINT32_C(1) << 24u, 0u);
    for (uint32_t id = 0u; id <= 257u; ++id)
    {
        const Direction scalar = ScalarDirectionForId(id);
        const Direction codec = CodecDirectionForId(id);
        if (!scalar.Nonzero || !SameDirection(scalar, codec)) {
            error = "low candidate direction mismatch";
            return false;
        }
        const uint32_t key = DirectionKey(scalar);
        if (seen[key]) {
            error = "low hyperoval contains a duplicate point";
            return false;
        }
        seen[key] = 1u;
    }
    for (uint32_t bucket = 0u; bucket < 65535u; ++bucket)
    {
        const Direction scalar = ScalarComplementDirection(bucket);
        K3HyperovalDirection codec_value;
        if (!wirehair_wh2_bench::K3HyperovalComplementDirectionForBucket(
                bucket, codec_value))
        {
            error = "codec rejected a valid complement bucket";
            return false;
        }
        const Direction codec = NormalizeDirection(
            codec_value.Alpha, codec_value.Beta, codec_value.Gamma);
        const uint32_t key = DirectionKey(scalar);
        if (!scalar.Nonzero || !SameDirection(scalar, codec) ||
            codec_value.Point != 258u + bucket || seen[key]) {
            error = "complement zero/duplicate/overlap/mapping failure";
            return false;
        }
        seen[key] = 1u;
    }
    K3HyperovalDirection sentinel;
    sentinel.Point = UINT32_MAX;
    sentinel.Alpha = 19u;
    sentinel.Beta = 23u;
    sentinel.Gamma = 29u;
    if (wirehair_wh2_bench::K3HyperovalComplementDirectionForBucket(
            65535u, sentinel) ||
        sentinel.Point != UINT32_MAX || sentinel.Alpha != 19u ||
        sentinel.Beta != 23u || sentinel.Gamma != 29u)
    {
        error = "invalid complement bucket was not transactional";
        return false;
    }
    if (!wirehair_wh2_bench::
            K3HyperovalSquareMapperMatchesInitializedGf256())
    {
        error = "candidate square mapper disagrees with initialized GF256";
        return false;
    }
    return true;
}

bool ClassifyLowExhaustive(
    WirehairCodec legacy,
    WirehairV2Codec v2,
    ComparisonCounts& counts,
    std::string& error)
{
    std::array<Direction, 258u> candidate;
    std::array<Direction, 258u> legacy_directions;
    std::array<Direction, 258u> v2_directions;
    for (uint32_t id = 0u; id < 258u; ++id)
    {
        candidate[id] = ScalarDirectionForId(id);
        if (!SameDirection(candidate[id], CodecDirectionForId(id)) ||
            !ExtractLegacyDirection(
                legacy, id, legacy_directions[id], error) ||
            !ExtractV2Direction(v2, id, v2_directions[id], error))
        {
            return false;
        }
    }
    counts = ComparisonCounts();
    for (uint32_t first = 0u; first < 258u; ++first) {
        for (uint32_t second = first + 1u; second < 258u; ++second) {
            for (uint32_t third = second + 1u; third < 258u; ++third) {
                if (!CrossCheckDeterminant(
                        candidate[first], candidate[second], candidate[third]) ||
                    !CrossCheckDeterminant(
                        legacy_directions[first], legacy_directions[second],
                        legacy_directions[third]) ||
                    !CrossCheckDeterminant(
                        v2_directions[first], v2_directions[second],
                        v2_directions[third]))
                {
                    error = "low exhaustive determinant oracle mismatch";
                    return false;
                }
                counts.LegacyZeroRowOccurrences +=
                    legacy_directions[first].Nonzero ? 0u : 1u;
                counts.LegacyZeroRowOccurrences +=
                    legacy_directions[second].Nonzero ? 0u : 1u;
                counts.LegacyZeroRowOccurrences +=
                    legacy_directions[third].Nonzero ? 0u : 1u;
                counts.V2ZeroRowOccurrences +=
                    v2_directions[first].Nonzero ? 0u : 1u;
                counts.V2ZeroRowOccurrences +=
                    v2_directions[second].Nonzero ? 0u : 1u;
                counts.V2ZeroRowOccurrences +=
                    v2_directions[third].Nonzero ? 0u : 1u;
                RecordOutcomes(
                    ScalarDeterminant(
                        candidate[first], candidate[second], candidate[third]) ==
                        0u,
                    ScalarDeterminant(
                        legacy_directions[first], legacy_directions[second],
                        legacy_directions[third]) == 0u,
                    ScalarDeterminant(
                        v2_directions[first], v2_directions[second],
                        v2_directions[third]) == 0u,
                    counts);
            }
        }
    }
    if (counts.Count != kLowTriples) {
        error = "low exhaustive count mismatch";
        return false;
    }
    return true;
}

std::vector<uint8_t> DeterministicSource(
    uint32_t block_bytes,
    uint64_t seed,
    uint64_t message_bytes = 0u)
{
    if (message_bytes == 0u) {
        message_bytes = static_cast<uint64_t>(kK) * block_bytes;
    }
    std::vector<uint8_t> source(static_cast<size_t>(message_bytes), 0u);
    SplitMix64 random(seed);
    size_t offset = 0u;
    while (offset < source.size())
    {
        uint64_t word = random.Next();
        for (unsigned byte = 0u; byte != 8u && offset < source.size();
             ++byte, ++offset)
        {
            source[offset] = static_cast<uint8_t>(word >> (byte * 8u));
        }
    }
    return source;
}

bool EncodeLegacyPacket(
    WirehairCodec codec,
    uint32_t id,
    uint32_t block_bytes,
    std::vector<uint8_t>& packet,
    uint32_t& written,
    std::string& error)
{
    packet.assign(block_bytes, 0u);
    written = 0u;
    const WirehairResult result = wirehair_encode(
        codec, id, packet.data(), block_bytes, &written);
    if (result != Wirehair_Success || written == 0u || written > block_bytes) {
        error = "Wirehair1 packet encode failed";
        return false;
    }
    return true;
}

bool EncodeV2Packet(
    WirehairV2Codec codec,
    uint32_t id,
    uint32_t block_bytes,
    std::vector<uint8_t>& packet,
    uint32_t& written,
    std::string& error)
{
    packet.assign(block_bytes, 0u);
    written = 0u;
    const WirehairV2Result result = wirehair_v2_encode(
        codec, id, packet.data(), block_bytes, &written);
    if (result != WirehairV2_Success || written == 0u || written > block_bytes) {
        error = "WH2 packet encode failed";
        return false;
    }
    return true;
}

bool CheckSuperposition(
    const PublicProfiles& profiles,
    WirehairCodec legacy_identity,
    WirehairV2Codec v2_identity,
    std::string& error)
{
    const std::vector<uint8_t> source = DeterministicSource(
        kComparatorBlockBytes, UINT64_C(0x13198a2e03707344));
    LegacyOwner legacy;
    V2Owner v2;
    if (!InitializeLegacyEncoder(
            profiles.Legacy, source, kComparatorBlockBytes, legacy, error) ||
        !InitializeV2EncoderFromProfile(profiles.V2, source, v2, error))
    {
        return false;
    }
    std::vector<uint32_t> ids;
    for (uint32_t id = 0u; id <= 511u; ++id) ids.push_back(id);
    for (uint32_t id = UINT32_MAX - 3u;; ++id) {
        ids.push_back(id);
        if (id == UINT32_MAX) break;
    }
    for (uint32_t id : ids)
    {
        Direction legacy_direction;
        Direction v2_direction;
        // Rank comparisons use canonical projective directions, but encoded
        // bytes retain the equation row's scalar.  Preserve the raw identity
        // lanes here so the superposition oracle checks the exact packet, not
        // merely a projectively equivalent row.
        if (!ExtractLegacyDirection(
                legacy_identity, id, legacy_direction, error, false) ||
            !ExtractV2Direction(
                v2_identity, id, v2_direction, error, false))
        {
            return false;
        }
        std::vector<uint8_t> legacy_packet;
        std::vector<uint8_t> v2_packet;
        uint32_t legacy_bytes = 0u;
        uint32_t v2_bytes = 0u;
        if (!EncodeLegacyPacket(
                legacy.Codec, id, kComparatorBlockBytes,
                legacy_packet, legacy_bytes, error) ||
            !EncodeV2Packet(
                v2.Codec, id, kComparatorBlockBytes,
                v2_packet, v2_bytes, error) ||
            legacy_bytes != kComparatorBlockBytes ||
            v2_bytes != kComparatorBlockBytes)
        {
            return false;
        }
        for (uint32_t byte = 0u; byte < kComparatorBlockBytes; ++byte)
        {
            const uint8_t expected_legacy =
                ScalarMultiply(legacy_direction.A, source[byte]) ^
                ScalarMultiply(
                    legacy_direction.B,
                    source[kComparatorBlockBytes + byte]) ^
                ScalarMultiply(
                    legacy_direction.C,
                    source[2u * kComparatorBlockBytes + byte]);
            const uint8_t expected_v2 =
                ScalarMultiply(v2_direction.A, source[byte]) ^
                ScalarMultiply(
                    v2_direction.B,
                    source[kComparatorBlockBytes + byte]) ^
                ScalarMultiply(
                    v2_direction.C,
                    source[2u * kComparatorBlockBytes + byte]);
            if (legacy_packet[byte] != expected_legacy ||
                v2_packet[byte] != expected_v2)
            {
                error = "public coefficient superposition mismatch";
                return false;
            }
        }
    }
    return true;
}

bool RunLegacyDecode(
    const WirehairWireProfile& profile,
    const Triple& triple,
    const std::vector<std::vector<uint8_t> >& packets,
    const std::vector<uint32_t>& packet_bytes,
    const std::vector<uint8_t>& source,
    bool predicted_full_rank,
    const unsigned order[3],
    std::string& error)
{
    LegacyOwner decoder;
    const WirehairResult create = wirehair_decoder_create_profile_ex(
        nullptr,
        source.size(),
        kComparatorBlockBytes,
        &profile,
        &decoder.Codec);
    if (create != Wirehair_Success || !decoder.Codec) {
        error = "Wirehair1 concordance decoder construction failed";
        return false;
    }
    WirehairResult result = Wirehair_NeedMore;
    for (unsigned position = 0u; position != 3u; ++position)
    {
        const unsigned selected = order[position];
        result = wirehair_decode(
            decoder.Codec,
            triple.Id[selected],
            packets[selected].data(),
            packet_bytes[selected]);
        const WirehairResult expected =
            predicted_full_rank && position == 2u ?
                Wirehair_Success : Wirehair_NeedMore;
        if (result != expected) {
            error = "Wirehair1 determinant/decode concordance mismatch";
            return false;
        }
    }
    std::vector<uint8_t> recovered(source.size(), 0u);
    const WirehairResult recovery = wirehair_recover(
        decoder.Codec, recovered.data(), recovered.size());
    if (predicted_full_rank) {
        if (recovery != Wirehair_Success || recovered != source) {
            error = "Wirehair1 successful concordance recovery mismatch";
            return false;
        }
    }
    else if (recovery != Wirehair_NeedMore) {
        error = "Wirehair1 deficient concordance recovery completed";
        return false;
    }
    return true;
}

bool RunV2Decode(
    const std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES>& profile,
    const Triple& triple,
    const std::vector<std::vector<uint8_t> >& packets,
    const std::vector<uint32_t>& packet_bytes,
    const std::vector<uint8_t>& source,
    bool predicted_full_rank,
    const unsigned order[3],
    std::string& error)
{
    V2Owner decoder;
    const WirehairV2Result create = wirehair_v2_decoder_create(
        profile.data(), profile.size(), &decoder.Codec);
    if (create != WirehairV2_Success || !decoder.Codec) {
        error = "WH2 concordance decoder construction failed";
        return false;
    }
    WirehairV2Result result = WirehairV2_NeedMore;
    for (unsigned position = 0u; position != 3u; ++position)
    {
        const unsigned selected = order[position];
        result = wirehair_v2_decode(
            decoder.Codec,
            triple.Id[selected],
            packets[selected].data(),
            packet_bytes[selected]);
        const WirehairV2Result expected =
            predicted_full_rank && position == 2u ?
                WirehairV2_Success : WirehairV2_NeedMore;
        if (result != expected) {
            error = "WH2 determinant/decode concordance mismatch";
            return false;
        }
    }
    std::vector<uint8_t> recovered(source.size(), 0u);
    uint64_t recovered_bytes = 0u;
    const WirehairV2Result recovery = wirehair_v2_recover(
        decoder.Codec,
        recovered.data(),
        recovered.size(),
        &recovered_bytes);
    if (predicted_full_rank) {
        if (recovery != WirehairV2_Success ||
            recovered_bytes != source.size() || recovered != source)
        {
            error = "WH2 successful concordance recovery mismatch";
            return false;
        }
    }
    else if (recovery != WirehairV2_NeedMore ||
             recovered_bytes != source.size())
    {
        error = "WH2 deficient concordance recovery completed";
        return false;
    }
    return true;
}

bool OperationalConcordance(
    const ManifestData& manifest,
    const PublicProfiles& profiles,
    WirehairCodec legacy_identity,
    WirehairV2Codec v2_identity,
    std::string& error)
{
    std::vector<Triple> triples;
    triples.reserve(2048u);
    uint32_t low_added = 0u;
    for (uint32_t first = 0u; first < 258u && low_added < 256u; ++first) {
        for (uint32_t second = first + 1u;
             second < 258u && low_added < 256u;
             ++second)
        {
            for (uint32_t third = second + 1u;
                 third < 258u && low_added < 256u;
                 ++third)
            {
                Triple triple;
                triple.Id[0] = first;
                triple.Id[1] = second;
                triple.Id[2] = third;
                triples.push_back(triple);
                ++low_added;
            }
        }
    }
    for (uint32_t family = 0u;
         family < static_cast<uint32_t>(Family::Count);
         ++family)
    {
        if (manifest.Corpora[family].size() < 256u) {
            error = "operational corpus prefix missing";
            return false;
        }
        triples.insert(
            triples.end(),
            manifest.Corpora[family].begin(),
            manifest.Corpora[family].begin() + 256u);
    }
    if (triples.size() != 2048u) {
        error = "operational concordance count is not 2048";
        return false;
    }

    const std::vector<uint8_t> source = DeterministicSource(
        kComparatorBlockBytes, UINT64_C(0xa4093822299f31d0));
    LegacyOwner legacy_encoder;
    V2Owner v2_encoder;
    if (!InitializeLegacyEncoder(
            profiles.Legacy, source, kComparatorBlockBytes,
            legacy_encoder, error) ||
        !InitializeV2EncoderFromProfile(
            profiles.V2, source, v2_encoder, error))
    {
        return false;
    }

    bool legacy_independent_permuted = false;
    bool legacy_dependent_permuted = false;
    bool v2_independent_permuted = false;
    bool v2_dependent_permuted = false;
    static const unsigned orders[6][3] = {
        {0u,1u,2u}, {0u,2u,1u}, {1u,0u,2u},
        {1u,2u,0u}, {2u,0u,1u}, {2u,1u,0u}
    };
    for (const Triple& triple : triples)
    {
        Direction ignored_candidate[3];
        Direction legacy_directions[3];
        Direction v2_directions[3];
        if (!DirectionsForTriple(
                triple,
                legacy_identity,
                v2_identity,
                ignored_candidate,
                legacy_directions,
                v2_directions,
                error))
        {
            return false;
        }
        const bool legacy_full = ScalarDeterminant(
            legacy_directions[0], legacy_directions[1],
            legacy_directions[2]) != 0u;
        const bool v2_full = ScalarDeterminant(
            v2_directions[0], v2_directions[1], v2_directions[2]) != 0u;
        std::vector<std::vector<uint8_t> > legacy_packets(3u);
        std::vector<std::vector<uint8_t> > v2_packets(3u);
        std::vector<uint32_t> legacy_bytes(3u, 0u);
        std::vector<uint32_t> v2_bytes(3u, 0u);
        for (unsigned i = 0u; i != 3u; ++i) {
            if (!EncodeLegacyPacket(
                    legacy_encoder.Codec, triple.Id[i],
                    kComparatorBlockBytes, legacy_packets[i],
                    legacy_bytes[i], error) ||
                !EncodeV2Packet(
                    v2_encoder.Codec, triple.Id[i],
                    kComparatorBlockBytes, v2_packets[i],
                    v2_bytes[i], error))
            {
                return false;
            }
        }
        if (!RunLegacyDecode(
                profiles.Legacy, triple, legacy_packets, legacy_bytes,
                source, legacy_full, orders[0], error) ||
            !RunV2Decode(
                profiles.V2, triple, v2_packets, v2_bytes,
                source, v2_full, orders[0], error))
        {
            return false;
        }

        if ((legacy_full && !legacy_independent_permuted) ||
            (!legacy_full && !legacy_dependent_permuted))
        {
            for (const auto& order : orders) {
                if (!RunLegacyDecode(
                        profiles.Legacy, triple, legacy_packets, legacy_bytes,
                        source, legacy_full, order, error)) return false;
            }
            if (legacy_full) legacy_independent_permuted = true;
            else legacy_dependent_permuted = true;
        }
        if ((v2_full && !v2_independent_permuted) ||
            (!v2_full && !v2_dependent_permuted))
        {
            for (const auto& order : orders) {
                if (!RunV2Decode(
                        profiles.V2, triple, v2_packets, v2_bytes,
                        source, v2_full, order, error)) return false;
            }
            if (v2_full) v2_independent_permuted = true;
            else v2_dependent_permuted = true;
        }
    }
    if (!legacy_independent_permuted || !legacy_dependent_permuted ||
        !v2_independent_permuted || !v2_dependent_permuted)
    {
        error = "operational corpus lacked comparator permutation representatives";
        return false;
    }
    return true;
}

bool EncodeCandidatePackets(
    K3HyperovalCodec& encoder,
    const std::vector<uint32_t>& ids,
    std::vector<std::vector<uint8_t> >& packets,
    std::vector<uint32_t>& packet_bytes,
    std::string& error)
{
    packets.assign(ids.size(), std::vector<uint8_t>(encoder.BlockBytes(), 0u));
    packet_bytes.assign(ids.size(), 0u);
    for (size_t i = 0u; i < ids.size(); ++i)
    {
        uint32_t written = UINT32_MAX;
        const WirehairResult result = encoder.EncodeResult(
            ids[i], packets[i].data(), packets[i].size(), &written);
        if (result != Wirehair_Success || written == 0u ||
            written > encoder.BlockBytes())
        {
            error = "candidate packet encode failed";
            return false;
        }
        packet_bytes[i] = written;
    }
    return true;
}

bool CandidateDecodeOrder(
    uint64_t message_bytes,
    uint32_t block_bytes,
    const std::vector<uint8_t>& source,
    const Triple& triple,
    const std::vector<std::vector<uint8_t> >& packets,
    const std::vector<uint32_t>& packet_bytes,
    bool full_rank,
    const unsigned order[3],
    std::string& error)
{
    K3HyperovalCodec decoder;
    if (decoder.InitializeDecoder(message_bytes, block_bytes) !=
            Wirehair_Success)
    {
        error = "candidate decoder construction failed";
        return false;
    }
    for (unsigned position = 0u; position != 3u; ++position)
    {
        const unsigned selected = order[position];
        const WirehairResult result = decoder.DecodeResult(
            triple.Id[selected],
            packets[selected].data(),
            packet_bytes[selected]);
        const WirehairResult expected = full_rank && position == 2u ?
            Wirehair_Success : Wirehair_NeedMore;
        if (result != expected) {
            error = "candidate determinant/decode concordance mismatch";
            return false;
        }
    }
    std::vector<uint8_t> recovered(source.size(), 0u);
    const WirehairResult recovery = decoder.RecoverResult(
        recovered.data(), recovered.size());
    if (full_rank) {
        if (recovery != Wirehair_Success || recovered != source) {
            error = "candidate recovered bytes mismatch";
            return false;
        }
    }
    else if (recovery != Wirehair_NeedMore) {
        error = "candidate deficient triple recovered";
        return false;
    }
    return true;
}

bool FrozenCandidateRecovery(
    const FrozenManifest& frozen,
    std::string& error)
{
    if (frozen.Entries.size() != 12u) {
        error = "authenticated frozen corpus does not contain 12 cells";
        return false;
    }
    for (const FrozenEntry& entry : frozen.Entries)
    {
        const std::vector<uint8_t> source = DeterministicSource(
            entry.Cell.block_bytes,
            entry.Cell.loss_seed ^ UINT64_C(0x6a09e667f3bcc909));
        K3HyperovalCodec encoder;
        if (encoder.InitializeEncoder(
                source.data(), source.size(), entry.Cell.block_bytes) !=
                Wirehair_Success)
        {
            error = "frozen candidate encoder construction failed";
            return false;
        }
        std::vector<std::vector<uint8_t> > all_packets;
        std::vector<uint32_t> all_packet_bytes;
        if (!EncodeCandidatePackets(
                encoder, entry.Ids, all_packets, all_packet_bytes, error))
        {
            return false;
        }
        for (unsigned first = 0u; first < 7u; ++first) {
            for (unsigned second = first + 1u; second < 7u; ++second) {
                for (unsigned third = second + 1u; third < 7u; ++third) {
                    Triple triple;
                    triple.Id[0] = entry.Ids[first];
                    triple.Id[1] = entry.Ids[second];
                    triple.Id[2] = entry.Ids[third];
                    const Direction d0 = ScalarDirectionForId(triple.Id[0]);
                    const Direction d1 = ScalarDirectionForId(triple.Id[1]);
                    const Direction d2 = ScalarDirectionForId(triple.Id[2]);
                    if (!CrossCheckDeterminant(d0, d1, d2)) {
                        error = "frozen candidate determinant oracle mismatch";
                        return false;
                    }
                    if (ScalarDeterminant(d0, d1, d2) == 0u) {
                        error = "frozen candidate subset lacks rank three";
                        return false;
                    }
                    std::vector<std::vector<uint8_t> > packets(3u);
                    std::vector<uint32_t> packet_bytes(3u, 0u);
                    const unsigned selected[3] = { first, second, third };
                    for (unsigned i = 0u; i != 3u; ++i) {
                        packets[i] = all_packets[selected[i]];
                        packet_bytes[i] = all_packet_bytes[selected[i]];
                    }
                    static const unsigned order[3] = { 0u, 1u, 2u };
                    if (!CandidateDecodeOrder(
                            source.size(), entry.Cell.block_bytes, source,
                            triple, packets, packet_bytes, true, order, error))
                    {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool FindRawIdForBucket(
    uint32_t bucket,
    uint32_t& id_out,
    uint32_t avoid = UINT32_MAX,
    uint32_t avoid_second = UINT32_MAX)
{
    for (uint64_t id = 258u; id <= UINT32_MAX; ++id) {
        const uint32_t candidate = static_cast<uint32_t>(id);
        if (candidate != avoid && candidate != avoid_second &&
            ComplementBucketForId(candidate) == bucket)
        {
            id_out = candidate;
            return true;
        }
    }
    return false;
}

bool CandidatePermutationRepresentatives(std::string& error)
{
    static const unsigned orders[6][3] = {
        {0u,1u,2u}, {0u,2u,1u}, {1u,0u,2u},
        {1u,2u,0u}, {2u,0u,1u}, {2u,1u,0u}
    };
    const std::vector<uint8_t> source = DeterministicSource(
        64u, UINT64_C(0x082efa98ec4e6c89));
    K3HyperovalCodec encoder;
    if (encoder.InitializeEncoder(source.data(), source.size(), 64u) !=
            Wirehair_Success)
    {
        error = "candidate permutation encoder construction failed";
        return false;
    }
    uint32_t secant_id = 0u;
    if (!FindRawIdForBucket(255u, secant_id)) {
        error = "could not locate frozen secant representative";
        return false;
    }
    Triple representatives[2];
    representatives[0].Id[0] = 0u;
    representatives[0].Id[1] = 1u;
    representatives[0].Id[2] = 2u;
    representatives[1].Id[0] = 0u;
    representatives[1].Id[1] = 1u;
    representatives[1].Id[2] = secant_id;
    for (unsigned kind = 0u; kind != 2u; ++kind)
    {
        std::vector<uint32_t> ids(
            representatives[kind].Id,
            representatives[kind].Id + 3u);
        std::vector<std::vector<uint8_t> > packets;
        std::vector<uint32_t> packet_bytes;
        if (!EncodeCandidatePackets(
                encoder, ids, packets, packet_bytes, error)) return false;
        const bool full_rank = kind == 0u;
        const Direction d0 = ScalarDirectionForId(ids[0]);
        const Direction d1 = ScalarDirectionForId(ids[1]);
        const Direction d2 = ScalarDirectionForId(ids[2]);
        if ((ScalarDeterminant(d0, d1, d2) != 0u) != full_rank) {
            error = "candidate permutation representative rank mismatch";
            return false;
        }
        for (const auto& order : orders) {
            if (!CandidateDecodeOrder(
                    source.size(), 64u, source, representatives[kind],
                    packets, packet_bytes, full_rank, order, error))
            {
                return false;
            }
        }
    }
    return true;
}

bool CandidateOverheadRecovery(uint32_t block_bytes, std::string& error)
{
    const std::vector<uint8_t> source = DeterministicSource(
        block_bytes,
        UINT64_C(0x452821e638d01377) ^ block_bytes);
    K3HyperovalCodec encoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success)
    {
        error = "candidate width encoder construction failed";
        return false;
    }
    std::vector<uint32_t> ids;
    for (uint32_t id = 0u; id != 7u; ++id) ids.push_back(id);
    std::vector<std::vector<uint8_t> > packets;
    std::vector<uint32_t> packet_bytes;
    if (!EncodeCandidatePackets(encoder, ids, packets, packet_bytes, error)) {
        return false;
    }
    for (uint32_t overhead = 0u; overhead <= 4u; ++overhead)
    {
        K3HyperovalCodec decoder;
        if (decoder.InitializeDecoder(source.size(), block_bytes) !=
                Wirehair_Success)
        {
            error = "candidate OH decoder construction failed";
            return false;
        }
        for (uint32_t i = 0u; i < kK + overhead; ++i)
        {
            const WirehairResult result = decoder.DecodeResult(
                ids[i], packets[i].data(), packet_bytes[i]);
            const WirehairResult expected = i < 2u ?
                Wirehair_NeedMore : Wirehair_Success;
            if (result != expected) {
                error = "candidate OH0..4 result mismatch";
                return false;
            }
        }
        std::vector<uint8_t> recovered(source.size(), 0u);
        if (decoder.RecoverResult(recovered.data(), recovered.size()) !=
                Wirehair_Success || recovered != source)
        {
            error = "candidate OH0..4 recovered bytes mismatch";
            return false;
        }
    }

    const uint64_t tail_bytes = block_bytes == 2u ? 1u : block_bytes - 1u;
    const uint64_t partial_bytes = 2u * static_cast<uint64_t>(block_bytes) +
        tail_bytes;
    const std::vector<uint8_t> partial_source = DeterministicSource(
        block_bytes,
        UINT64_C(0xbe5466cf34e90c6c) ^ block_bytes,
        partial_bytes);
    K3HyperovalCodec partial_encoder;
    if (partial_encoder.InitializeEncoder(
            partial_source.data(), partial_source.size(), block_bytes) !=
            Wirehair_Success)
    {
        error = "candidate partial-tail encoder construction failed";
        return false;
    }
    std::vector<uint8_t> final_packet(block_bytes, 0xa5u);
    uint32_t final_written = UINT32_MAX;
    if (partial_encoder.EncodeResult(
            2u, final_packet.data(), final_packet.size(), &final_written) !=
            Wirehair_Success || final_written != tail_bytes ||
        std::memcmp(
            final_packet.data(),
            partial_source.data() + 2u * block_bytes,
            static_cast<size_t>(tail_bytes)) != 0 ||
        !std::all_of(
            final_packet.begin() + static_cast<size_t>(tail_bytes),
            final_packet.end(),
            [](uint8_t value) { return value == UINT8_C(0xa5); }))
    {
        error = "candidate partial systematic packet mismatch";
        return false;
    }
    const uint32_t repair_ids[3] = { 3u, 4u, 5u };
    std::vector<std::vector<uint8_t> > repair_packets(3u);
    std::vector<uint32_t> repair_bytes(3u, 0u);
    for (unsigned i = 0u; i != 3u; ++i) {
        repair_packets[i].assign(block_bytes, 0u);
        if (partial_encoder.EncodeResult(
                repair_ids[i], repair_packets[i].data(), block_bytes,
                &repair_bytes[i]) != Wirehair_Success ||
            repair_bytes[i] != block_bytes)
        {
            error = "candidate partial-tail repair encode mismatch";
            return false;
        }
    }
    const Direction r0 = ScalarDirectionForId(repair_ids[0]);
    const Direction r1 = ScalarDirectionForId(repair_ids[1]);
    const Direction r2 = ScalarDirectionForId(repair_ids[2]);
    if (ScalarDeterminant(r0, r1, r2) == 0u)
    {
        error = "candidate partial-tail repair triple is rank deficient";
        return false;
    }
    {
        K3HyperovalCodec decoder;
        if (decoder.InitializeDecoder(partial_bytes, block_bytes) !=
                Wirehair_Success)
        {
            error = "candidate partial-tail decoder construction failed";
            return false;
        }
        for (unsigned i = 0u; i != 3u; ++i) {
            const WirehairResult result = decoder.DecodeResult(
                repair_ids[i], repair_packets[i].data(), repair_bytes[i]);
            if (result != (i == 2u ? Wirehair_Success : Wirehair_NeedMore)) {
                error = "candidate partial-tail repair decode mismatch";
                return false;
            }
        }
        std::vector<uint8_t> recovered(partial_source.size(), 0u);
        if (decoder.RecoverResult(recovered.data(), recovered.size()) !=
                Wirehair_Success || recovered != partial_source)
        {
            error = "candidate partial-tail padding recovery mismatch";
            return false;
        }
    }

    const uint32_t retry_ids_array[3] = { 0u, 1u, 3u };
    const std::vector<uint32_t> retry_ids(
        retry_ids_array, retry_ids_array + 3u);
    std::vector<std::vector<uint8_t> > retry_packets;
    std::vector<uint32_t> retry_packet_bytes;
    if (!EncodeCandidatePackets(
            partial_encoder,
            retry_ids,
            retry_packets,
            retry_packet_bytes,
            error))
    {
        return false;
    }
    K3HyperovalCodec retry_decoder;
    if (retry_decoder.InitializeDecoder(partial_bytes, block_bytes) !=
            Wirehair_Success ||
        retry_decoder.DecodeResult(
            retry_ids[0], retry_packets[0].data(), retry_packet_bytes[0]) !=
            Wirehair_NeedMore ||
        retry_decoder.DecodeResult(
            retry_ids[1], retry_packets[1].data(), retry_packet_bytes[1]) !=
            Wirehair_NeedMore)
    {
        error = "candidate padding retry rank-two setup failed";
        return false;
    }
    std::vector<uint8_t> corrupt_padding = retry_packets[2];
    corrupt_padding[block_bytes - 1u] ^= 1u;
    if (retry_decoder.DecodeResult(
            retry_ids[2],
            corrupt_padding.data(),
            retry_packet_bytes[2]) != Wirehair_Error ||
        retry_decoder.Rank() != 2u || retry_decoder.AcceptedIdCount() != 2u ||
        retry_decoder.IsDecoded() ||
        retry_decoder.DecodeResult(
            retry_ids[2],
            retry_packets[2].data(),
            retry_packet_bytes[2]) != Wirehair_Success)
    {
        error = "candidate corrupt-padding failure/retry was not transactional";
        return false;
    }
    std::vector<uint8_t> retry_recovered(partial_source.size(), 0u);
    if (retry_decoder.RecoverResult(
            retry_recovered.data(), retry_recovered.size()) !=
            Wirehair_Success || retry_recovered != partial_source)
    {
        error = "candidate corrupt-padding retry recovery mismatch";
        return false;
    }
    return true;
}

bool CandidateSemanticChecks(uint32_t block_bytes, std::string& error)
{
    const std::vector<uint8_t> source = DeterministicSource(
        block_bytes,
        UINT64_C(0xc0ac29b7c97c50dd) ^
            static_cast<uint64_t>(block_bytes));
    K3HyperovalCodec encoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success)
    {
        error = "candidate semantic encoder construction failed";
        return false;
    }

    uint32_t alias_first = 0u;
    uint32_t alias_second = 0u;
    uint32_t alias_third = 0u;
    if (!FindRawIdForBucket(0u, alias_first) ||
        !FindRawIdForBucket(0u, alias_second, alias_first) ||
        !FindRawIdForBucket(
            0u, alias_third, alias_first, alias_second) ||
        alias_first == alias_second || alias_first == alias_third ||
        alias_second == alias_third)
    {
        error = "candidate alias representatives unavailable";
        return false;
    }
    std::vector<uint32_t> ids;
    ids.push_back(alias_first);
    ids.push_back(alias_second);
    ids.push_back(alias_third);
    ids.push_back(0u);
    ids.push_back(1u);
    ids.push_back(2u);
    ids.push_back(UINT32_MAX);
    std::vector<std::vector<uint8_t> > packets;
    std::vector<uint32_t> packet_bytes;
    if (!EncodeCandidatePackets(encoder, ids, packets, packet_bytes, error)) {
        return false;
    }
    if (packets[0] != packets[1] || packets[0] != packets[2]) {
        error = "candidate same-direction aliases encoded differently";
        return false;
    }

    K3HyperovalCodec alias_decoder;
    if (alias_decoder.InitializeDecoder(source.size(), block_bytes) !=
            Wirehair_Success ||
        alias_decoder.DecodeResult(
            ids[0], packets[0].data(), packet_bytes[0]) != Wirehair_NeedMore ||
        alias_decoder.DecodeResult(
            ids[0], packets[0].data(), packet_bytes[0]) != Wirehair_NeedMore ||
        alias_decoder.DecodeResult(
            ids[1], packets[1].data(), packet_bytes[1]) != Wirehair_NeedMore)
    {
        error = "candidate raw duplicate/alias idempotence mismatch";
        return false;
    }
    std::vector<uint8_t> conflict = packets[0];
    conflict[0] ^= 1u;
    if (alias_decoder.DecodeResult(
            ids[0], conflict.data(), packet_bytes[0]) !=
            Wirehair_InvalidInput ||
        alias_decoder.DecodeResult(
            ids[1], conflict.data(), packet_bytes[1]) !=
            Wirehair_InvalidInput ||
        alias_decoder.DecodeResult(
            ids[2], conflict.data(), packet_bytes[2]) != Wirehair_Error ||
        alias_decoder.DecodeResult(
            ids[2], packets[2].data(), packet_bytes[2]) != Wirehair_NeedMore)
    {
        error = "candidate duplicate/alias conflict classification mismatch";
        return false;
    }

    uint32_t secant_id = 0u;
    if (!FindRawIdForBucket(255u, secant_id)) {
        error = "candidate dependent-row representative unavailable";
        return false;
    }
    const uint32_t dependent_ids_array[3] = { 0u, 1u, secant_id };
    const std::vector<uint32_t> dependent_ids(
        dependent_ids_array, dependent_ids_array + 3u);
    std::vector<std::vector<uint8_t> > dependent_packets;
    std::vector<uint32_t> dependent_packet_bytes;
    if (!EncodeCandidatePackets(
            encoder,
            dependent_ids,
            dependent_packets,
            dependent_packet_bytes,
            error))
    {
        return false;
    }
    K3HyperovalCodec dependent_decoder;
    if (dependent_decoder.InitializeDecoder(source.size(), block_bytes) !=
            Wirehair_Success ||
        dependent_decoder.DecodeResult(
            dependent_ids[0], dependent_packets[0].data(),
            dependent_packet_bytes[0]) != Wirehair_NeedMore ||
        dependent_decoder.DecodeResult(
            dependent_ids[1], dependent_packets[1].data(),
            dependent_packet_bytes[1]) != Wirehair_NeedMore)
    {
        error = "candidate dependent-row rank-two setup failed";
        return false;
    }
    std::vector<uint8_t> dependent_conflict = dependent_packets[2];
    dependent_conflict[0] ^= 1u;
    if (dependent_decoder.DecodeResult(
            dependent_ids[2], dependent_conflict.data(),
            dependent_packet_bytes[2]) != Wirehair_Error ||
        dependent_decoder.Rank() != 2u ||
        dependent_decoder.AcceptedIdCount() != 2u ||
        dependent_decoder.DecodeResult(
            dependent_ids[2], dependent_packets[2].data(),
            dependent_packet_bytes[2]) != Wirehair_NeedMore ||
        dependent_decoder.DecodeResult(
            dependent_ids[2], dependent_conflict.data(),
            dependent_packet_bytes[2]) != Wirehair_InvalidInput)
    {
        error = "candidate inconsistent dependent-row semantics mismatch";
        return false;
    }

    K3HyperovalCodec completed;
    if (completed.InitializeDecoder(source.size(), block_bytes) !=
            Wirehair_Success)
    {
        error = "candidate post-success decoder construction failed";
        return false;
    }
    for (unsigned i = 3u; i != 6u; ++i) {
        const WirehairResult result = completed.DecodeResult(
            ids[i], packets[i].data(), packet_bytes[i]);
        if (result != (i == 5u ? Wirehair_Success : Wirehair_NeedMore)) {
            error = "candidate systematic completion mismatch";
            return false;
        }
    }
    std::vector<uint8_t> max_conflict = packets[6];
    max_conflict[0] ^= 1u;
    if (completed.DecodeResult(
            ids[6], max_conflict.data(), packet_bytes[6]) != Wirehair_Error ||
        completed.DecodeResult(
            ids[6], packets[6].data(), packet_bytes[6]) != Wirehair_Success)
    {
        error = "candidate UINT32_MAX post-success validation failed";
        return false;
    }
    if (completed.DecodeResult(
            ids[6], max_conflict.data(), packet_bytes[6]) !=
            Wirehair_InvalidInput)
    {
        error = "candidate post-success conflicting raw ID mismatch";
        return false;
    }

    K3HyperovalCodec cap_decoder;
    if (cap_decoder.InitializeDecoder(source.size(), block_bytes) !=
            Wirehair_Success)
    {
        error = "candidate accepted-ID-cap decoder construction failed";
        return false;
    }
    for (unsigned i = 3u; i != 6u; ++i)
    {
        const WirehairResult result = cap_decoder.DecodeResult(
            ids[i], packets[i].data(), packet_bytes[i]);
        if (result != (i == 5u ? Wirehair_Success : Wirehair_NeedMore)) {
            error = "candidate accepted-ID-cap setup failed";
            return false;
        }
    }
    std::vector<uint8_t> first_cap_packet;
    uint32_t first_cap_bytes = 0u;
    for (uint32_t id = 3u; id <= 1026u; ++id)
    {
        std::vector<uint8_t> packet(block_bytes, 0u);
        uint32_t data_bytes = 0u;
        if (encoder.EncodeResult(
                id, packet.data(), packet.size(), &data_bytes) !=
                Wirehair_Success || data_bytes != block_bytes ||
            cap_decoder.DecodeResult(
                id, packet.data(), data_bytes) != Wirehair_Success)
        {
            error = "candidate accepted-ID-cap fill failed";
            return false;
        }
        if (id == 3u) {
            first_cap_packet = packet;
            first_cap_bytes = data_bytes;
        }
    }
    std::vector<uint8_t> beyond_cap(block_bytes, 0u);
    uint32_t beyond_cap_bytes = 0u;
    std::vector<uint8_t> retained_conflict = first_cap_packet;
    retained_conflict[0] ^= 1u;
    if (cap_decoder.AcceptedIdCount() !=
            K3HyperovalCodec::kMaxAcceptedPacketIds ||
        encoder.EncodeResult(
            1027u,
            beyond_cap.data(),
            beyond_cap.size(),
            &beyond_cap_bytes) != Wirehair_Success ||
        beyond_cap_bytes != block_bytes ||
        cap_decoder.DecodeResult(
            1027u, beyond_cap.data(), beyond_cap_bytes) !=
            Wirehair_ExtraInsufficient ||
        cap_decoder.DecodeResult(
            3u, first_cap_packet.data(), first_cap_bytes) != Wirehair_Success ||
        cap_decoder.DecodeResult(
            3u, retained_conflict.data(), first_cap_bytes) !=
            Wirehair_InvalidInput)
    {
        error = "candidate accepted-ID-cap semantics mismatch";
        return false;
    }

    encoder.ResetWorkCounters();
    std::vector<uint8_t> repair(block_bytes, 0u);
    uint32_t written = 0u;
    if (encoder.EncodeResult(
            258u, repair.data(), repair.size(), &written) != Wirehair_Success ||
        written != block_bytes)
    {
        error = "candidate work-ceiling repair encode failed";
        return false;
    }
    const wirehair_wh2_bench::K3HyperovalWorkCounters encode_work =
        encoder.WorkCounters();
    if (encode_work.PacketEvaluation.Copies != 1u ||
        encode_work.PacketEvaluation.Scales != 0u ||
        encode_work.PacketEvaluation.Xors +
                encode_work.PacketEvaluation.MulAdds > 2u ||
        encode_work.DependentValidation.Copies != 0u ||
        encode_work.DependentValidation.Xors != 0u ||
        encode_work.DependentValidation.Scales != 0u ||
        encode_work.DependentValidation.MulAdds != 0u ||
        encode_work.Solve.Copies != 0u || encode_work.Solve.Xors != 0u ||
        encode_work.Solve.Scales != 0u || encode_work.Solve.MulAdds != 0u)
    {
        error = "candidate repair evaluation exceeded static work ceiling";
        return false;
    }

    const uint32_t dense_ids_array[3] = { 3u, 4u, 5u };
    const std::vector<uint32_t> dense_ids(
        dense_ids_array, dense_ids_array + 3u);
    std::vector<std::vector<uint8_t> > dense_packets;
    std::vector<uint32_t> dense_packet_bytes;
    if (!EncodeCandidatePackets(
            encoder, dense_ids, dense_packets, dense_packet_bytes, error) ||
        ScalarDeterminant(
            ScalarDirectionForId(dense_ids[0]),
            ScalarDirectionForId(dense_ids[1]),
            ScalarDirectionForId(dense_ids[2])) == 0u)
    {
        error = "candidate dense work-ceiling triple invalid";
        return false;
    }

    K3HyperovalCodec solve_decoder;
    if (solve_decoder.InitializeDecoder(source.size(), block_bytes) !=
            Wirehair_Success)
    {
        error = "candidate work-ceiling decoder construction failed";
        return false;
    }
    solve_decoder.ResetWorkCounters();
    for (unsigned i = 0u; i != 3u; ++i) {
        const WirehairResult result = solve_decoder.DecodeResult(
            dense_ids[i], dense_packets[i].data(), dense_packet_bytes[i]);
        if (result != (i == 2u ? Wirehair_Success : Wirehair_NeedMore)) {
            error = "candidate work-ceiling decode failed";
            return false;
        }
    }
    const wirehair_wh2_bench::K3HyperovalWorkCounters solve_work =
        solve_decoder.WorkCounters();
    if (solve_work.PacketEvaluation.Copies != 0u ||
        solve_work.PacketEvaluation.Xors != 0u ||
        solve_work.PacketEvaluation.Scales != 0u ||
        solve_work.PacketEvaluation.MulAdds != 0u ||
        solve_work.DependentValidation.Copies != 0u ||
        solve_work.DependentValidation.Xors != 0u ||
        solve_work.DependentValidation.Scales != 0u ||
        solve_work.DependentValidation.MulAdds != 0u ||
        solve_work.Solve.Copies + solve_work.Solve.Scales != 3u ||
        solve_work.Solve.Xors + solve_work.Solve.MulAdds > 6u)
    {
        error = "candidate 3x3 solve exceeded static work ceiling";
        return false;
    }
    return true;
}

bool SelfTestAllocationHooks(std::string& error)
{
#if defined(WIREHAIR_WH2_K3_HYPEROVAL_TEST_HOOKS)
    const std::vector<uint8_t> source = DeterministicSource(8u, 123u);

    wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(0);
    K3HyperovalCodec encoder;
    const WirehairResult encoder_result = encoder.InitializeEncoder(
        source.data(), source.size(), 8u);
    wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(-1);
    if (encoder_result != Wirehair_OOM || encoder.IsEncoder() ||
        encoder.IsDecoder()) {
        error = "candidate encoder allocation failure was not transactional";
        return false;
    }

    K3HyperovalCodec preserved_encoder;
    if (preserved_encoder.InitializeEncoder(
            source.data(), source.size(), 8u) != Wirehair_Success)
    {
        error = "candidate preserved encoder setup failed";
        return false;
    }
    const std::vector<uint8_t> replacement = DeterministicSource(16u, 456u);
    wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(0);
    const WirehairResult reinit_encoder = preserved_encoder.InitializeEncoder(
        replacement.data(), replacement.size(), 16u);
    wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(-1);
    std::array<uint8_t, 8u> original_packet = {};
    uint32_t original_bytes = 0u;
    if (reinit_encoder != Wirehair_OOM || !preserved_encoder.IsEncoder() ||
        preserved_encoder.BlockBytes() != 8u ||
        preserved_encoder.MessageBytes() != source.size() ||
        preserved_encoder.EncodeResult(
            0u,
            original_packet.data(),
            original_packet.size(),
            &original_bytes) != Wirehair_Success ||
        original_bytes != original_packet.size() ||
        std::memcmp(
            original_packet.data(), source.data(), original_packet.size()) != 0)
    {
        error = "candidate failed encoder reinit did not preserve state";
        return false;
    }

    for (int64_t countdown = 0; countdown <= 1; ++countdown)
    {
        wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(
            countdown);
        K3HyperovalCodec decoder;
        const WirehairResult decoder_result =
            decoder.InitializeDecoder(24u, 8u);
        wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(
            -1);
        if (decoder_result != Wirehair_OOM || decoder.IsEncoder() ||
            decoder.IsDecoder())
        {
            error = "candidate decoder allocation point was not transactional";
            return false;
        }
    }

    K3HyperovalCodec preserved_decoder;
    if (preserved_decoder.InitializeDecoder(24u, 8u) != Wirehair_Success ||
        preserved_decoder.DecodeResult(
            0u, original_packet.data(), original_bytes) != Wirehair_NeedMore)
    {
        error = "candidate preserved decoder setup failed";
        return false;
    }
    for (int64_t countdown = 0; countdown <= 1; ++countdown)
    {
        wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(
            countdown);
        const WirehairResult decoder_result =
            preserved_decoder.InitializeDecoder(48u, 16u);
        wirehair_wh2_bench::SetK3HyperovalAllocationFailureCountdownForTesting(
            -1);
        if (decoder_result != Wirehair_OOM ||
            !preserved_decoder.IsDecoder() || preserved_decoder.IsDecoded() ||
            preserved_decoder.BlockBytes() != 8u ||
            preserved_decoder.MessageBytes() != 24u ||
            preserved_decoder.Rank() != 1u ||
            preserved_decoder.AcceptedIdCount() != 1u ||
            preserved_decoder.DecodeResult(
                0u, original_packet.data(), original_bytes) != Wirehair_NeedMore)
        {
            error = "candidate failed decoder reinit did not preserve state";
            return false;
        }
    }
#else
    (void)error;
#endif
    return true;
}

bool CheckDomainGate(const ComparisonCounts& counts)
{
    const WilsonInterval candidate = Wilson(
        counts.CandidateFailures, counts.Count);
    return counts.Count == kStochasticTriples &&
        candidate.Upper < 0.01L &&
        counts.CandidateFailures <= counts.LegacyFailures &&
        counts.CandidateFailures <= counts.V2Failures;
}

void PrintCounts(const char* family, const ComparisonCounts& counts)
{
    const WilsonInterval candidate = Wilson(
        counts.CandidateFailures, counts.Count);
    const WilsonInterval legacy = Wilson(counts.LegacyFailures, counts.Count);
    const WilsonInterval v2 = Wilson(counts.V2Failures, counts.Count);
    std::cout << std::setprecision(12)
              << "RANK,family=" << family
              << ",n=" << counts.Count
              << ",candidate_fail=" << counts.CandidateFailures
              << ",wirehair1_fail=" << counts.LegacyFailures
              << ",wh2_fail=" << counts.V2Failures
              << ",wirehair1_zero_row_occurrences="
              << counts.LegacyZeroRowOccurrences
              << ",wh2_zero_row_occurrences="
              << counts.V2ZeroRowOccurrences
              << ",candidate_l95=" << static_cast<double>(candidate.Lower)
              << ",candidate_u95=" << static_cast<double>(candidate.Upper)
              << ",wirehair1_l95=" << static_cast<double>(legacy.Lower)
              << ",wirehair1_u95=" << static_cast<double>(legacy.Upper)
              << ",wh2_l95=" << static_cast<double>(v2.Lower)
              << ",wh2_u95=" << static_cast<double>(v2.Upper)
              << '\n';
}

int RunRankRecovery()
{
    ManifestData manifest;
    AnalyticReceipt analytic;
    std::string error;
    if (!BuildManifestData(true, manifest, error) ||
        !AuthenticateManifest(manifest, error) ||
        !BuildAnalyticReceipt(analytic, error))
    {
        std::cerr << "K3 rank/recovery invalid before classification: "
                  << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
    }
    if (wirehair_init() != Wirehair_Success) {
        std::cerr << "K3 rank/recovery global initialization failed\n";
        return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
    }
    if (!ValidateComplementAndMapper(error)) {
        std::cerr << "K3 candidate mapping rejected: " << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=REJECT", 2);
    }

    const std::vector<uint8_t> identity_source =
        IdentityLaneMessage(kComparatorBlockBytes);
    PublicProfiles profiles;
    V2Owner selected_v2;
    if (!SelectPublicProfiles(
            identity_source, profiles, selected_v2, error))
    {
        std::cerr << "K3 public comparator selection invalid: " << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
    }
    LegacyOwner legacy_identity;
    if (!InitializeLegacyEncoder(
            profiles.Legacy, identity_source, kComparatorBlockBytes,
            legacy_identity, error))
    {
        std::cerr << "K3 Wirehair1 comparator invalid: " << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
    }
    std::cout << "CONFIG,protocol=" << kProtocolTag
              << ",git_commit=" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
              << ",source_sha256="
              << WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256
              << ",codec_cpp_sha256="
              << WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256
              << ",codec_header_sha256="
              << WIREHAIR_WH2_K3_CODEC_HEADER_SHA256
              << ",workers=" << kWorkerCount
              << ",K=" << kK
              << ",B=" << kComparatorBlockBytes
              << ",wirehair1_profile=0x" << std::hex
              << profiles.Legacy.profile_id
              << ",wh2_profile=0x" << WIREHAIR_V2_PROFILE_CURRENT
              << std::dec
              << ",wh2_attempt=" << static_cast<uint32_t>(profiles.V2Attempt)
              << ",wh2_descriptor_sha256=" << profiles.V2Hash << '\n';
    std::cout << "AUTH,protocol=" << kProtocolTag
              << ",kind=low-only,sha256=" << manifest.LowHash << '\n';
    for (uint32_t value = 0u;
         value < static_cast<uint32_t>(Family::Count);
         ++value)
    {
        std::cout << "AUTH,protocol=" << kProtocolTag
                  << ",kind=" << FamilyName(static_cast<Family>(value))
                  << ",sha256=" << manifest.FamilyHashes[value] << '\n';
    }
    std::cout << "AUTH,protocol=" << kProtocolTag
              << ",kind=frozen-canonical,sha256="
              << manifest.Frozen.CanonicalHash << '\n'
              << "AUTH,protocol=" << kProtocolTag
              << ",kind=frozen-trace,sha256="
              << manifest.Frozen.TraceHash << '\n'
              << "AUTH,protocol=" << kProtocolTag
              << ",kind=complement-points,sha256="
              << manifest.ComplementHash << '\n'
              << "ANALYTIC,secant_lines=" << analytic.SecantLines
              << ",external_lines=" << analytic.ExternalLines
              << ",duplicate=" << analytic.DuplicateNumerator << '/'
              << analytic.DuplicateDenominator
              << ",distinct_collinear="
              << analytic.DistinctCollinearNumerator << '/'
              << analytic.DistinctCollinearDenominator
              << ",total=" << analytic.TotalNumerator << '/'
              << analytic.TotalDenominator
              << ",two_low=" << analytic.TwoLowNumerator << '/'
              << analytic.TwoLowDenominator << '\n';

    if (!CheckSuperposition(
            profiles,
            legacy_identity.Codec,
            selected_v2.Codec,
            error))
    {
        std::cerr << "K3 public extraction invalid: " << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
    }

    ComparisonCounts low;
    if (!ClassifyLowExhaustive(
            legacy_identity.Codec, selected_v2.Codec, low, error))
    {
        std::cerr << "K3 low classification invalid: " << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
    }
    PrintCounts("low-only", low);
    bool scientific_pass = low.CandidateFailures == 0u &&
        low.CandidateFailures <= low.LegacyFailures &&
        low.CandidateFailures <= low.V2Failures;

    std::array<ComparisonCounts, static_cast<size_t>(Family::Count)> domain;
    ComparisonCounts pooled;
    for (uint32_t value = 0u;
         value < static_cast<uint32_t>(Family::Count);
         ++value)
    {
        if (!ClassifyCorpus(
                manifest.Corpora[value], profiles, identity_source,
                domain[value], error))
        {
            std::cerr << "K3 stochastic classification invalid: "
                      << error << '\n';
            return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
        }
        PrintCounts(FamilyName(static_cast<Family>(value)), domain[value]);
        scientific_pass = scientific_pass && CheckDomainGate(domain[value]);
        pooled.Add(domain[value]);
    }

    const WilsonInterval pooled_candidate = Wilson(
        pooled.CandidateFailures, pooled.Count);
    const WilsonInterval pooled_legacy = Wilson(
        pooled.LegacyFailures, pooled.Count);
    const WilsonInterval pooled_v2 = Wilson(
        pooled.V2Failures, pooled.Count);
    const long double legacy_mcnemar = ExactMcNemarUpperTail(
        pooled.CandidateOnlyLegacy, pooled.LegacyOnlyCandidate);
    const long double v2_mcnemar = ExactMcNemarUpperTail(
        pooled.CandidateOnlyV2, pooled.V2OnlyCandidate);
    const bool pooled_gate =
        pooled_candidate.Upper < pooled_legacy.Lower &&
        pooled_candidate.Upper < pooled_v2.Lower &&
        pooled.LegacyOnlyCandidate > pooled.CandidateOnlyLegacy &&
        pooled.V2OnlyCandidate > pooled.CandidateOnlyV2 &&
        legacy_mcnemar < 0.025L && v2_mcnemar < 0.025L;
    scientific_pass = scientific_pass && pooled_gate;
    std::cout << std::setprecision(12)
              << "POOLED,n=" << pooled.Count
              << ",candidate_fail=" << pooled.CandidateFailures
              << ",wirehair1_fail=" << pooled.LegacyFailures
              << ",wh2_fail=" << pooled.V2Failures
              << ",wirehair1_zero_row_occurrences="
              << pooled.LegacyZeroRowOccurrences
              << ",wh2_zero_row_occurrences="
              << pooled.V2ZeroRowOccurrences
              << ",candidate_u95="
              << static_cast<double>(pooled_candidate.Upper)
              << ",wirehair1_l95=" << static_cast<double>(pooled_legacy.Lower)
              << ",wh2_l95=" << static_cast<double>(pooled_v2.Lower)
              << ",candidate_only_wirehair1="
              << pooled.CandidateOnlyLegacy
              << ",wirehair1_only_candidate="
              << pooled.LegacyOnlyCandidate
              << ",candidate_only_wh2=" << pooled.CandidateOnlyV2
              << ",wh2_only_candidate=" << pooled.V2OnlyCandidate
              << ",wirehair1_mcnemar_p="
              << static_cast<double>(legacy_mcnemar)
              << ",wh2_mcnemar_p=" << static_cast<double>(v2_mcnemar)
              << ",gate=" << (pooled_gate ? "PASS" : "REJECT") << '\n';

    if (!OperationalConcordance(
            manifest, profiles, legacy_identity.Codec,
            selected_v2.Codec, error))
    {
        std::cerr << "K3 operational concordance invalid: " << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=INVALID", 1);
    }
    if (!FrozenCandidateRecovery(manifest.Frozen, error) ||
        !CandidatePermutationRepresentatives(error) ||
        !CandidateOverheadRecovery(2u, error) ||
        !CandidateOverheadRecovery(64u, error) ||
        !CandidateOverheadRecovery(1280u, error) ||
        !CandidateSemanticChecks(2u, error) ||
        !CandidateSemanticChecks(64u, error) ||
        !CandidateSemanticChecks(1280u, error))
    {
        std::cerr << "K3 candidate recovery/semantics rejected: "
                  << error << '\n';
        return PublishTerminal("K3_RANK_RECOVERY=REJECT", 2);
    }
    std::cout << "VALIDATION,operational_concordance=2048"
              << ",frozen_cells=12,frozen_subsets=420"
              << ",widths=2:64:1280,status=PASS\n";
    return PublishTerminal(
        scientific_pass ?
            "K3_RANK_RECOVERY=PASS" : "K3_RANK_RECOVERY=REJECT",
        scientific_pass ? 0 : 2);
}

bool SelfTestScalarAndDirections(std::string& error)
{
    for (uint32_t left = 0u; left < 256u; ++left) {
        for (uint32_t right = 0u; right < 256u; ++right) {
            if (ScalarMultiply(
                    static_cast<uint8_t>(left),
                    static_cast<uint8_t>(right)) !=
                gf256_mul(
                    static_cast<uint8_t>(left),
                    static_cast<uint8_t>(right)))
            {
                error = "scalar GF256 oracle mismatch";
                return false;
            }
        }
    }
    for (uint32_t id = 0u; id <= 511u; ++id) {
        if (!SameDirection(
                ScalarDirectionForId(id), CodecDirectionForId(id)))
        {
            error = "small direction mapper mismatch";
            return false;
        }
    }
    return true;
}

bool SelfTestSmallCorpora(std::string& error)
{
    if (!ExpectedManifestConfigured()) {
        error = "frozen manifest definitions are missing";
        return false;
    }
    for (uint32_t value = 0u;
         value < static_cast<uint32_t>(Family::Count);
         ++value)
    {
        std::vector<Triple> first;
        std::vector<Triple> second;
        const Family family = static_cast<Family>(value);
        if (!GenerateFamily(family, 16u, first, error) ||
            !GenerateFamily(family, 16u, second, error) ||
            first.size() != 16u || second.size() != 16u)
        {
            error = "small corpus determinism mismatch";
            return false;
        }
        for (size_t i = 0u; i < first.size(); ++i) {
            for (unsigned j = 0u; j != 3u; ++j) {
                if (first[i].Id[j] != second[i].Id[j]) {
                    error = "small corpus determinism mismatch";
                    return false;
                }
            }
        }
    }
    return true;
}

bool SelfTestCandidateSmall(std::string& error)
{
    const std::vector<uint8_t> source = DeterministicSource(8u, 999u);
    K3HyperovalCodec encoder;
    K3HyperovalCodec decoder;
    if (encoder.InitializeEncoder(source.data(), source.size(), 8u) !=
            Wirehair_Success ||
        decoder.InitializeDecoder(source.size(), 8u) != Wirehair_Success)
    {
        error = "small candidate construction failed";
        return false;
    }
    const uint32_t ids[3] = { 0u, 1u, 2u };
    for (unsigned i = 0u; i != 3u; ++i)
    {
        std::array<uint8_t, 8u> packet = {};
        uint32_t written = 0u;
        if (encoder.EncodeResult(
                ids[i], packet.data(), packet.size(), &written) !=
                Wirehair_Success || written != packet.size() ||
            decoder.DecodeResult(ids[i], packet.data(), written) !=
                (i == 2u ? Wirehair_Success : Wirehair_NeedMore))
        {
            error = "small candidate encode/decode failed";
            return false;
        }
    }
    std::vector<uint8_t> recovered(source.size(), 0u);
    if (decoder.RecoverResult(recovered.data(), recovered.size()) !=
            Wirehair_Success || recovered != source)
    {
        error = "small candidate recovery failed";
        return false;
    }
    return true;
}

bool SelfTestStatistics(std::string& error)
{
    AnalyticReceipt analytic;
    const WilsonInterval zero = Wilson(0u, kStochasticTriples);
    const WilsonInterval full = Wilson(kStochasticTriples, kStochasticTriples);
    if (!BuildAnalyticReceipt(analytic, error) ||
        !(zero.Lower <= 1e-18L && zero.Upper > 0.0L &&
          zero.Upper < 0.01L) ||
        !(full.Upper >= 1.0L - 1e-18L && full.Lower > 0.99L) ||
        ExactMcNemarUpperTail(0u, 0u) != 1.0L ||
        !(ExactMcNemarUpperTail(0u, 32u) < 0.025L) ||
        ExactMcNemarUpperTail(32u, 0u) != 1.0L)
    {
        error = "statistics selftest mismatch";
        return false;
    }
    return true;
}

int RunSelfTest()
{
    if (wirehair_init() != Wirehair_Success) {
        std::cerr << "K3 selftest global initialization failed\n";
        return PublishTerminal("K3_HYPEROVAL_SELFTEST=FAIL", 1);
    }
    struct Test { const char* Name; bool (*Run)(std::string&); };
    const Test tests[] = {
        { "scalar-directions", SelfTestScalarAndDirections },
        { "small-corpora", SelfTestSmallCorpora },
        { "statistics", SelfTestStatistics },
        { "candidate-small", SelfTestCandidateSmall },
        { "allocation-hooks", SelfTestAllocationHooks }
    };
    for (const Test& test : tests)
    {
        std::string error;
        if (!test.Run(error)) {
            std::cerr << "K3 selftest failed: " << test.Name
                      << ": " << error << '\n';
            return PublishTerminal("K3_HYPEROVAL_SELFTEST=FAIL", 1);
        }
        std::cout << "SELFTEST,name=" << test.Name << ",status=PASS\n";
    }
    return PublishTerminal("K3_HYPEROVAL_SELFTEST=PASS", 0);
}

int MainImpl(int argc, char** argv)
{
    if (argc != 2 || !argv || !argv[1]) {
        std::cerr << "usage: wirehair_wh2_k3_hyperoval_falsifier "
                     "--selftest|--manifest|--rank-recovery\n";
        return 1;
    }
    const std::string mode(argv[1]);
    if (mode == "--selftest") return RunSelfTest();
    if (mode == "--manifest") return RunManifest();
    if (mode == "--rank-recovery") return RunRankRecovery();
    std::cerr << "unsupported K3 falsifier mode\n";
    return 1;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        return MainImpl(argc, argv);
    }
    catch (const std::bad_alloc&) {
        std::cerr << "K3 falsifier allocation failed\n";
        return 1;
    }
    catch (const std::length_error&) {
        std::cerr << "K3 falsifier allocation length failed\n";
        return 1;
    }
    catch (const std::exception& exception) {
        std::cerr << "K3 falsifier exception: " << exception.what() << '\n';
        return 1;
    }
    catch (...) {
        std::cerr << "K3 falsifier failed unexpectedly\n";
        return 1;
    }
}
