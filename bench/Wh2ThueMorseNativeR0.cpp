// Benchmark-only, fixed .63 native experiment.  No production wire profile.
#include "Wh2ThueMorseNativeCodec.h"
#include "Wh2ThueMorseNativeDataR0.h"
#include "Wh2FrozenTrace.h"
#include "Wh2NativePanel.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "gf256.h"
#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>

namespace {
namespace N = wh2_thue_native;
namespace D = wh2_thue_native_data;
namespace P = wirehair_wh2_bench;
using wirehair::wh2_benchmark::Sha256Hex;
using Clock = std::chrono::steady_clock;
static const int kCpu = 50;
static const uint32_t kBatch = 64;
static const uint8_t kCanary = 0xa5;
static const size_t kGuard = 16;
static const char kProtocol[] = "wirehair.wh2.thue-morse-native-r0";
static const uint32_t kWidths[] = {2, 64, 1280};
static const uint32_t kConditions[4][4] = {
    {0, 1, 3, 2}, {1, 2, 0, 3}, {2, 3, 1, 0}, {3, 0, 2, 1}};
static const uint32_t kComparisons[5][2] = {{0,0},{1,1},{2,2},{0,1},{0,2}};
static const uint32_t kDistant[6] = {
    4294967295u,4294967293u,4294967291u,4294967289u,4294967287u,4294967285u};
enum Scope { EncoderCreate, DecoderCreate, Systematic, Sequential, Distant,
             ReceiveRecover, DecodeEndToEnd };
static const char* const kScopes[] = {"encoder_create", "decoder_create",
    "prebuilt_systematic", "prebuilt_sequential_repair", "prebuilt_distant_repair",
    "receive_recover", "decode_end_to_end"};

void Require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

std::string Quote(const std::string& value)
{
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c >= 32 && c < 127) out << c;
        else out << "\\u00" << std::hex << std::setw(2) << std::setfill('0')
                 << static_cast<unsigned>(c) << std::dec;
    }
    out << '"';
    return out.str();
}

void Little(std::vector<uint8_t>& bytes, uint64_t value, unsigned count)
{
    for (unsigned i = 0; i < count; ++i) {
        bytes.push_back(static_cast<uint8_t>(value)); value >>= 8;
    }
}

uint64_t Address(const void* pointer)
{
    return static_cast<uint64_t>(reinterpret_cast<uintptr_t>(pointer));
}

bool OnCpu()
{
    cpu_set_t mask;
    CPU_ZERO(&mask);
    return sched_getaffinity(0, sizeof(mask), &mask) == 0 &&
        CPU_COUNT(&mask) == 1 && CPU_ISSET(kCpu, &mask) && sched_getcpu() == kCpu;
}

void Pin()
{
    cpu_set_t mask;
    CPU_ZERO(&mask); CPU_SET(kCpu, &mask);
    Require(sched_setaffinity(0, sizeof(mask), &mask) == 0 && OnCpu(), "singleton CPU50 pin");
}

std::string Identity()
{
    P::TargetIdentityReceiptV2 identity;
    std::string diagnostic, canonical;
    Require(P::CapturePublicBorrowedTargetIdentity(kCpu, identity, diagnostic),
            "target identity capture");
    Require(P::SerializeTargetIdentityV2(identity, canonical, diagnostic) &&
            canonical.size() <= 4096 && !canonical.empty() &&
            Sha256Hex(canonical) == identity.CanonicalSha256 &&
            identity.Before.Affinity == std::vector<int32_t>(1,kCpu) &&
            identity.After.Affinity == std::vector<int32_t>(1,kCpu), "target identity serialization");
    std::ostringstream hex;
    for (unsigned char byte : canonical)
        hex << std::hex << std::setfill('0') << std::setw(2) << static_cast<unsigned>(byte);
    const auto& raw=identity.Raw;
    const auto& derived=identity.Derived;
    std::ostringstream out;
    out << "{\"after_cpu\":" << identity.After.Cpu << ",\"before_cpu\":" << identity.Before.Cpu
        << ",\"canonical_bytes\":" << canonical.size() << ",\"canonical_hex\":" << Quote(hex.str())
        << ",\"canonical_sha256\":" << Quote(identity.CanonicalSha256)
        << ",\"capabilities\":{\"leaf1_ecx\":" << raw.Leaf1.Regs.Ecx
        << ",\"leaf1_edx\":" << raw.Leaf1.Regs.Edx << ",\"leaf6_eax\":" << raw.Leaf6.Regs.Eax
        << ",\"leaf6_ecx\":" << raw.Leaf6.Regs.Ecx << ",\"leaf80000001_ecx\":" << raw.Leaf80000001.Regs.Ecx
        << ",\"leaf80000001_edx\":" << raw.Leaf80000001.Regs.Edx
        << ",\"leaf80000008_ebx\":" << raw.Leaf80000008.Regs.Ebx
        << ",\"leaf80000021_eax\":" << raw.Leaf80000021.Regs.Eax
        << ",\"max_basic_leaf\":" << raw.Leaf0.Regs.Eax
        << ",\"max_extended_leaf\":" << raw.Leaf80000000.Regs.Eax
        << "},\"derived\":{\"ccd_id\":" << derived.CcdId << ",\"complex_id\":" << derived.ComplexId
        << ",\"core_id\":" << derived.CoreId << ",\"family\":" << derived.Family
        << ",\"full_apic_id\":" << derived.FullApicId << ",\"initial_apic_id_8\":" << derived.InitialApicId8
        << ",\"logical_processors_per_package\":" << derived.LogicalProcessorsPerPackage
        << ",\"model\":" << derived.Model << ",\"package_id\":" << derived.PackageId
        << ",\"stepping\":" << derived.Stepping << ",\"thread_id\":" << derived.ThreadId
        << ",\"threads_per_core\":" << derived.ThreadsPerCore
        << "},\"raw_capture_complete\":" << (identity.RawCaptureComplete ? "true" : "false")
        << ",\"requested_cpu\":" << identity.RequestedCpu
        << ",\"semantic_validation_passed\":" << (identity.SemanticValidationPassed ? "true" : "false")
        << ",\"singleton_affinity_verified\":true}";
    return out.str();
}

struct Budget {
    Clock::time_point started = Clock::now();
    void Check() const {
        Require(Clock::now() - started < std::chrono::seconds(120), "worker deadline");
        Require(OnCpu(), "worker CPU drift");
    }
    double Elapsed() const { return std::chrono::duration<double>(Clock::now() - started).count(); }
};

uint64_t SplitMix(uint64_t& state)
{
    uint64_t x = (state += UINT64_C(0x9e3779b97f4a7c15));
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

std::vector<uint8_t> Source(uint32_t B, uint32_t tail)
{
    const size_t M = 5u * static_cast<size_t>(B) + tail;
    std::vector<uint8_t> bytes(M);
    uint64_t state = UINT64_C(0x6e61746976653633) ^ M ^ (static_cast<uint64_t>(B) << 32);
    for (size_t pos = 0; pos < M;) {
        uint64_t word = SplitMix(state);
        for (unsigned j = 0; j < 8 && pos < M; ++j, ++pos) {
            bytes[pos] = static_cast<uint8_t>(word); word >>= 8;
        }
    }
    return bytes;
}

N::Lookup Lookup() { return N::Lookup{D::kLookup, sizeof(D::kLookup)}; }

// Independent polynomial long division and 8/8/8/8 prefix products.  Neither
// the production GF tables nor the candidate's 10/7/7/8 mapper are used here.
class Oracle {
    using Matrix = std::array<uint8_t, 36>;
    using Row = std::array<uint8_t, 6>;
    uint8_t mul_[256][256];
    Matrix table_[4][2][256];
    static uint8_t Product(unsigned a, unsigned b) {
        unsigned polynomial = 0;
        for (unsigned bit = 0; bit < 8; ++bit)
            if ((b >> bit) & 1u) polynomial ^= a << bit;
        for (int bit = 14; bit >= 8; --bit)
            if ((polynomial >> bit) & 1u) polynomial ^= 0x14du << (bit - 8);
        return static_cast<uint8_t>(polynomial);
    }
    static Matrix Unit() {
        Matrix result = {{0}};
        for (unsigned i = 0; i < 6; ++i) result[7 * i] = 1;
        return result;
    }
    Matrix Multiply(const Matrix& a, const Matrix& b) const {
        Matrix result = {{0}};
        for (unsigned r = 0; r < 6; ++r)
            for (unsigned c = 0; c < 6; ++c)
                for (unsigned k = 0; k < 6; ++k)
                    result[6*r+c] ^= mul_[a[6*r+k]][b[6*k+c]];
        return result;
    }
    Row Apply(const Matrix& matrix, const Row& vector) const {
        Row result = {{0}};
        for (unsigned r = 0; r < 6; ++r)
            for (unsigned k = 0; k < 6; ++k)
                result[r] ^= mul_[matrix[6*r+k]][vector[k]];
        return result;
    }
    static unsigned Parity(uint32_t x) {
        x ^= x >> 16; x ^= x >> 8; x ^= x >> 4;
        return (0x6996u >> (x & 15u)) & 1u;
    }
public:
    Oracle() {
        for (unsigned a = 0; a < 256; ++a)
            for (unsigned b = 0; b < 256; ++b) mul_[a][b] = Product(a,b);
        Matrix dyadic[32][2];
        for (unsigned s = 0; s < 2; ++s)
            for (unsigned r = 0; r < 6; ++r)
                for (unsigned c = 0; c < 6; ++c) dyadic[0][s][6*r+c] = D::kPair[s][r][c];
        for (unsigned level = 1; level < 32; ++level)
            for (unsigned s = 0; s < 2; ++s)
                dyadic[level][s] = Multiply(dyadic[level-1][s], dyadic[level-1][s^1]);
        for (unsigned group = 0; group < 4; ++group)
            for (unsigned s = 0; s < 2; ++s)
                for (unsigned value = 0; value < 256; ++value) {
                    Matrix result = Unit();
                    unsigned phase = s;
                    for (int bit = 7; bit >= 0; --bit) if ((value >> bit) & 1u) {
                        result = Multiply(result, dyadic[8*group+static_cast<unsigned>(bit)][phase]);
                        phase ^= 1;
                    }
                    table_[group][s][value] = result;
                }
    }
    Row Coefficients(uint32_t id) const {
        Row result = {{1,0,0,0,0,0}};
        for (unsigned group = 0; group < 4; ++group) {
            const unsigned phase = group == 3 ? 0 : Parity(id >> (8*(group+1)));
            result = Apply(table_[group][phase][(id >> (8*group)) & 255u], result);
        }
        return result;
    }
    std::vector<uint8_t> Packet(const std::vector<uint8_t>& source, uint32_t B,
                                uint32_t tail, uint32_t id) const {
        const Row row = Coefficients(id);
        std::vector<uint8_t> result(id == 5 ? tail : B, 0);
        for (size_t byte = 0; byte < result.size(); ++byte)
            for (size_t k = 0; k < 6; ++k)
                if (k*B+byte < source.size()) result[byte] ^= mul_[row[k]][source[k*B+byte]];
        return result;
    }
};

struct Encoder {
    unsigned role = 0;
    std::unique_ptr<N::Encoder> t;
    WirehairCodec l = nullptr;
    WirehairV2Codec p = nullptr;
    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> profile = {{0}};
    uint32_t profile_bytes = 0;
    uint64_t M = 0;
    uint32_t B = 0;
    Encoder() = default;
    Encoder(const Encoder&) = delete;
    Encoder& operator=(const Encoder&) = delete;
    ~Encoder() { Reset(); }
    void Reset() {
        t.reset();
        if (l) wirehair_free(l);
        if (p) wirehair_v2_free(p);
        l=nullptr; p=nullptr;
    }
    uint64_t Pointer() const { return role == 0 ? Address(t.get()) : role == 1 ? Address(l) : Address(p); }
    int64_t Create(unsigned which, const void* source, uint64_t bytes, uint32_t width) {
        Require(!Pointer(), "encoder reused without free");
        role = which; M = bytes; B = width;
        int64_t code;
        if (role == 0) code = static_cast<int64_t>(N::Encoder::Create(Lookup(), source, M, B, t));
        else if (role == 1) code = wirehair_encoder_create_ex(nullptr, source, M, B, &l);
        else {
            WirehairV2EncoderOptions options = {};
            options.struct_bytes = sizeof(options);
            options.options_version = WIREHAIR_V2_ENCODER_OPTIONS_VERSION;
            options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
            code = wirehair_v2_encoder_create_with_options(source, M, B, &options,
                profile.data(), static_cast<uint32_t>(profile.size()), &profile_bytes, &p);
        }
        return code == 0 && !Pointer() ? -1001 : code;
    }
    bool ValidProfile() const {
        if (!Pointer()) return false;
        if (role != 2) return true;
        WirehairV2Profile decoded = {};
        return profile_bytes == profile.size() &&
            wirehair_v2_profile_validate(profile.data(), profile_bytes) == WirehairV2_Success &&
            wirehair_v2_profile_deserialize(profile.data(), profile_bytes, &decoded) == WirehairV2_Success &&
            decoded.message_bytes == M && decoded.block_bytes == B &&
            decoded.profile_id == WIREHAIR_V2_PROFILE_CERTIFIED_2026_07;
    }
    int64_t Encode(uint32_t id, void* out, uint32_t capacity, uint32_t& bytes) const {
        if (role == 0) {
            const N::Result result = t->Encode(id, out, capacity);
            bytes = static_cast<uint32_t>(result.bytes_written);
            return static_cast<int64_t>(result.status);
        }
        if (role == 1) return wirehair_encode(l, id, out, capacity, &bytes);
        return wirehair_v2_encode(p, id, out, capacity, &bytes);
    }
};

struct Decoder {
    unsigned role = 0;
    std::unique_ptr<N::Decoder> t;
    WirehairCodec l = nullptr;
    WirehairV2Codec p = nullptr;
    uint64_t M = 0;
    Decoder() = default;
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;
    ~Decoder() { Reset(); }
    void Reset() {
        t.reset();
        if (l) wirehair_free(l);
        if (p) wirehair_v2_free(p);
        l=nullptr; p=nullptr;
    }
    uint64_t Pointer() const { return role == 0 ? Address(t.get()) : role == 1 ? Address(l) : Address(p); }
    int64_t Create(const Encoder& encoder) {
        Require(!Pointer(), "decoder reused without free");
        role = encoder.role; M = encoder.M;
        int64_t code;
        if (role == 0) code = static_cast<int64_t>(N::Decoder::Create(Lookup(), M, encoder.B, t));
        else if (role == 1) code = wirehair_decoder_create_ex(nullptr, M, encoder.B, &l);
        else code = wirehair_v2_decoder_create(encoder.profile.data(), encoder.profile_bytes, &p);
        return code == 0 && !Pointer() ? -1001 : code;
    }
    int64_t Feed(uint32_t id, const void* bytes, uint32_t count) {
        if (role == 0) return static_cast<int64_t>(t->Feed(id, bytes, count).status);
        if (role == 1) return wirehair_decode(l, id, bytes, count);
        return wirehair_v2_decode(p, id, bytes, count);
    }
    int64_t Recover(void* out, uint64_t capacity) {
        if (role == 0) {
            const N::Result result = t->Recover(out, static_cast<size_t>(capacity));
            return result.status == N::Status::Success && result.bytes_written != M ?
                -1002 : static_cast<int64_t>(result.status);
        }
        if (role == 1) return wirehair_recover(l, out, capacity);
        uint64_t bytes = 0;
        const int64_t code = wirehair_v2_recover(p, out, capacity, &bytes);
        return code == 0 && bytes != M ? -1002 : code;
    }
};

bool Canary(const uint8_t* bytes, size_t count)
{
    for (size_t i = 0; i < count; ++i) if (bytes[i] != kCanary) return false;
    return true;
}

struct ArmCheck {
    unsigned status = 0;
    int extra = -1;
    int64_t encoder = INT64_MIN, decoder = INT64_MIN, encode = INT64_MIN;
    int64_t feed = INT64_MIN, recover = INT64_MIN;
    unsigned fed = 0;
    static void Code(std::ostream& out, int64_t code) {
        if (code == INT64_MIN) out << "null"; else out << code;
    }
    void Json(std::ostream& out) const {
        out << '[' << status << ',';
        if (extra < 0) out << "null"; else out << extra;
        out << ','; Code(out, encoder); out << ','; Code(out, decoder);
        out << ','; Code(out, encode); out << ','; Code(out, feed);
        out << ','; Code(out, recover); out << ',' << fed << ']';
    }
};

struct Correctness {
    std::ostringstream records;
    size_t cases = 0, candidate_packets = 0, control_packets = 0, recovered_messages = 0;
    bool noncomparable = false;
    unsigned status = 0;
};

ArmCheck CheckArm(unsigned role, uint32_t B, uint32_t tail, const uint32_t* ids,
                  size_t count, const uint8_t* ranks, unsigned final_rank,
                  const Oracle& oracle, Correctness& totals)
{
    ArmCheck result;
    const std::vector<uint8_t> original = Source(B, tail);
    const size_t M=original.size();
    std::vector<uint8_t> source=original;
    source.resize(M+kGuard,kCanary);
    std::vector<std::vector<uint8_t>> packets(count), original_packets;
    auto finish = [&]() -> ArmCheck {
        if (!std::equal(original.begin(),original.end(),source.begin()) ||
            !Canary(source.data()+M,kGuard) ||
            (!original_packets.empty() && packets != original_packets)) result.status=7;
        for (size_t j=0;j<count;++j) if (!packets[j].empty()) {
            const size_t expected=ids[j] == 5 ? tail : B;
            if (!Canary(packets[j].data()+expected,packets[j].size()-expected)) result.status=7;
        }
        return result;
    };
    Encoder encoder;
    result.encoder = encoder.Create(role, source.data(), M, B);
    if (result.encoder != 0) { result.status = result.encoder < 0 ? 7 : 2; return finish(); }
    if (!encoder.ValidProfile()) { result.status = 7; return finish(); }
    std::vector<uint32_t> lengths(count);
    result.encode = 0;
    for (size_t j = 0; j < count; ++j) {
        const uint32_t expected = ids[j] == 5 ? tail : B;
        packets[j].assign(static_cast<size_t>(B) + kGuard, kCanary);
        uint32_t bytes = UINT32_MAX;
        result.encode = encoder.Encode(ids[j], packets[j].data(), expected, bytes);
        if (role == 0) ++totals.candidate_packets; else ++totals.control_packets;
        if (result.encode != 0) { result.status = result.encode < 0 ? 7 : 4; return finish(); }
        if (bytes != expected || !Canary(packets[j].data()+expected, packets[j].size()-expected) ||
            (role == 0 && !std::equal(packets[j].begin(), packets[j].begin()+bytes,
                                      oracle.Packet(original, B, tail, ids[j]).begin()))) {
            result.status = 7; return finish();
        }
        lengths[j] = bytes;
    }
    original_packets=packets;
    Decoder decoder;
    result.decoder = decoder.Create(encoder);
    if (result.decoder != 0) { result.status = result.decoder < 0 ? 7 : 3; return finish(); }
    result.feed = 1;
    for (size_t j = 0; j < count; ++j) {
        result.feed = decoder.Feed(ids[j], packets[j].data(), lengths[j]);
        ++result.fed;
        if (result.feed != 0 && result.feed != 1) { result.status = result.feed < 0 ? 7 : 5; return finish(); }
        if (role == 0 && (result.feed == 0) != (decoder.t->Rank() == 6)) {
            result.status = 7; return finish();
        }
        if (role == 0 && j >= 5 && ranks && decoder.t->Rank() != ranks[j-5]) {
            result.status = 7; return finish();
        }
        if (result.feed == 0) {
            if (j < 5) { result.status = 7; return finish(); }
            result.extra = static_cast<int>(j)-5;
            if (role == 0 && (final_rank != 6 || (ranks &&
                !std::all_of(ranks+j-5, ranks+5, [](uint8_t rank) { return rank == 6; })))) {
                result.status = 7; return finish();
            }
            std::vector<uint8_t> recovered(M+kGuard, kCanary);
            result.recover = decoder.Recover(recovered.data(), M);
            if (!Canary(recovered.data()+M,kGuard)) { result.status=7; return finish(); }
            if (result.recover != 0) { result.status = result.recover < 0 ? 7 : 6; return finish(); }
            if (!std::equal(original.begin(), original.end(), recovered.begin()) ||
                !Canary(recovered.data()+M, kGuard)) {
                result.status = 7; return finish();
            }
            ++totals.recovered_messages;
            break;
        }
    }
    if (result.feed == 1) {
        result.status = 1;
        if (role == 0 && (decoder.t->Rank() != final_rank || final_rank == 6)) result.status = 7;
    }
    return finish();
}

void CorrectCase(Correctness& result, unsigned kind, size_t index, uint32_t B, uint32_t tail,
                 const uint32_t* ids, size_t count, const uint8_t* ranks, unsigned final_rank,
                 const Oracle& oracle, const Budget& budget)
{
    budget.Check();
    std::ostringstream record;
    record << '[' << kind << ',' << index << ',' << B << ',' << tail;
    bool invalid = false;
    for (unsigned role = 0; role < 3; ++role) {
        const ArmCheck arm = CheckArm(role,B,tail,ids,count,ranks,final_rank,oracle,result);
        record << ','; arm.Json(record);
        invalid |= arm.status == 7 || (role == 0 && arm.status >= 2);
        if (role != 0 && arm.status >= 2) result.noncomparable = true;
    }
    record << ']';
    if (result.cases++) result.records << ',';
    result.records << record.str();
    if (invalid) { result.status = 2; throw std::runtime_error("native correctness disagreement"); }
}

std::string VisitHash(const Oracle& oracle, const Budget& budget)
{
    std::vector<uint8_t> visits;
    visits.reserve(623470);
    auto visit = [&](uint32_t id) {
        uint8_t row[6];
        Require(N::Row(Lookup(), id, row) == N::Status::Success, "native coefficient row");
        const auto expected = oracle.Coefficients(id);
        Require(std::equal(row,row+6,expected.begin()), "independent coefficient oracle");
        Little(visits,id,4); visits.insert(visits.end(),row,row+6);
    };
    for (const auto& trace : D::kTraces) { budget.Check(); for (uint32_t id : trace.ids) visit(id); }
    for (const auto& prefix : D::kHistory) for (unsigned j=0;j<prefix.count;++j) visit(prefix.ids[j]);
    for (const auto& prefix : D::kFixtures) for (unsigned j=0;j<prefix.count;++j) visit(prefix.ids[j]);
    Require(visits.size() == 623470, "coefficient visit count");
    const std::string digest = Sha256Hex(visits.data(), visits.size());
    Require(digest == D::kCoefficientVisitSha256, "coefficient visit SHA256");
    return digest;
}

void RunCorrectness(Correctness& result, const Oracle& oracle, const Budget& budget)
{
    for (size_t i = 0; i < D::kTraceCount; ++i) {
        const auto& trace = D::kTraces[i];
        CorrectCase(result,0,i,trace.B,trace.B,trace.ids,10,trace.ranks,trace.ranks[4],oracle,budget);
    }
    size_t ordinal=0;
    for (const auto& item : D::kHistoryCases) {
        const auto& prefix = D::kHistory[item.index];
        CorrectCase(result,1,ordinal++,item.B,item.tail,prefix.ids,prefix.count,nullptr,prefix.rank,oracle,budget);
    }
    ordinal=0;
    for (const auto& item : D::kFixtureCases) {
        const auto& prefix = D::kFixtures[item.index];
        CorrectCase(result,2,ordinal++,item.B,item.tail,prefix.ids,prefix.count,nullptr,prefix.rank,oracle,budget);
    }
    ordinal=0;
    for (const auto& item : D::kPartialCases) {
        const auto& trace = D::kTraces[item.index];
        CorrectCase(result,3,ordinal++,item.B,item.tail,trace.ids,10,trace.ranks,trace.ranks[4],oracle,budget);
    }
    Require(result.cases == 6456, "complete native correctness roster");
    result.status = result.noncomparable ? 1 : 0;
}

class Aligned {
    uint8_t* bytes_ = nullptr;
    size_t size_;
public:
    explicit Aligned(size_t size) : size_(size) {
        void* allocation = nullptr;
        Require(posix_memalign(&allocation,4096,size) == 0 && allocation, "aligned context allocation");
        bytes_ = static_cast<uint8_t*>(allocation);
        std::memset(bytes_,kCanary,size_);
    }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
    ~Aligned() { std::free(bytes_); }
    uint8_t* Data() { return bytes_; }
    const uint8_t* Data() const { return bytes_; }
    size_t Size() const { return size_; }
    void Reset() { std::memset(bytes_,kCanary,size_); }
};

struct PacketCache {
    int64_t code = 0;
    std::array<int64_t,3> encode_codes = {{0}};
    std::array<std::vector<uint8_t>,10> receive;
    std::array<std::array<std::vector<uint8_t>,6>,3> encode;
    unsigned first_count = 0;
};

struct Context {
    const uint32_t B;
    const size_t M, stride, output_stride;
    const std::vector<uint8_t> original;
    Aligned source, packet, output;
    std::array<Encoder,3> encoders;
    std::array<int64_t,3> create_codes = {{0}};
    std::array<PacketCache,3> caches;
    Context(uint32_t width, const Oracle& oracle) : B(width), M(6u*static_cast<size_t>(width)),
        stride(width+kGuard), output_stride(M+kGuard), original(Source(width,width)),
        source(M+kGuard), packet(10*stride), output(kBatch*output_stride) {
        std::memcpy(source.Data(),original.data(),M);
        for (unsigned role=0;role<3;++role) {
            Encoder& encoder = encoders[role];
            create_codes[role] = encoder.Create(role,source.Data(),M,B);
            PacketCache& cache = caches[role];
            cache.code = create_codes[role];
            cache.encode_codes.fill(cache.code);
            if (cache.code != 0) continue;
            Require(encoder.ValidProfile(), "timing encoder profile");
            for (unsigned scope=0;scope<3;++scope) for (unsigned j=0;j<6;++j) {
                auto& bytes = cache.encode[scope][j];
                bytes.resize(B);
                uint32_t length = UINT32_MAX;
                const uint32_t id = scope == 0 ? j : scope == 1 ? 6+j : kDistant[j];
                const int64_t code = encoder.Encode(id,bytes.data(),B,length);
                if (code != 0 || length != B) cache.encode_codes[scope] = code == 0 ? -1002 : code;
                if (role == 0) Require(code == 0 && length == B &&
                    bytes == oracle.Packet(original,B,B,id), "timing candidate encode oracle");
            }
            for (unsigned j=0;j<10;++j) {
                cache.receive[j].resize(B);
                uint32_t length = UINT32_MAX;
                const int64_t code = encoder.Encode(6+j,cache.receive[j].data(),B,length);
                if (code != 0 || length != B) cache.code = code == 0 ? -1002 : code;
                if (role == 0) Require(code == 0 && length == B &&
                    cache.receive[j] == oracle.Packet(original,B,B,6+j), "timing candidate receive packet oracle");
            }
            if (cache.code != 0) continue;
            Decoder decoder;
            cache.code = decoder.Create(encoder);
            if (cache.code != 0) continue;
            cache.code = 1;
            for (unsigned j=0;j<10;++j) {
                cache.code = decoder.Feed(6+j,cache.receive[j].data(),B);
                if (cache.code == 0) { cache.first_count=j+1; break; }
                if (cache.code != 1) break;
            }
            if (cache.code == 0) {
                std::vector<uint8_t> recovered(M);
                cache.code = decoder.Recover(recovered.data(),M);
                Require(cache.code != 0 || recovered == original, "timing input recovered bytes");
            }
        }
    }
    bool UnchangedSource() const {
        return std::equal(original.begin(),original.end(),source.Data()) &&
            Canary(source.Data()+M,kGuard);
    }
    void Prepare(unsigned role, unsigned scope) {
        Require(UnchangedSource(), "timing source modified");
        packet.Reset(); output.Reset();
        if (scope >= ReceiveRecover && (caches[role].code == 0 || caches[role].code == 1))
            for (unsigned j=0;j<10;++j)
                std::memcpy(packet.Data()+j*stride,caches[role].receive[j].data(),B);
    }
};

struct Counts {
    uint64_t callbacks=0, encoder_creates=0, decoder_creates=0, encode_calls=0,
             feed_calls=0, recover_calls=0;
    void Add(const Counts& other) {
        callbacks+=other.callbacks; encoder_creates+=other.encoder_creates;
        decoder_creates+=other.decoder_creates; encode_calls+=other.encode_calls;
        feed_calls+=other.feed_calls; recover_calls+=other.recover_calls;
    }
    void Json(std::ostream& out) const {
        out << "{\"callbacks\":" << callbacks << ",\"encoder_creates\":" << encoder_creates
            << ",\"decoder_creates\":" << decoder_creates << ",\"encode_calls\":" << encode_calls
            << ",\"feed_calls\":" << feed_calls << ",\"recover_calls\":" << recover_calls << '}';
    }
};

struct Allocations {
    std::vector<uint8_t> bytes;
    uint64_t count=0, first=0, last=0, minimum=UINT64_MAX, maximum=0;
    std::vector<Counts> callbacks;
    Allocations() { bytes.reserve(10*(4+8*kBatch)); callbacks.reserve(10); }
    void Record(const std::array<uint64_t,kBatch>& pointers, unsigned size, const Counts& work) {
        Require(callbacks.size() < 10 && size <= kBatch, "callback receipt bound");
        Little(bytes,size,4);
        for (unsigned i=0;i<size;++i) {
            const uint64_t pointer=pointers[i];
            Require(pointer != 0, "null recorded owned handle");
            Little(bytes,pointer,8);
            if (!count) first=pointer;
            ++count; last=pointer; minimum=std::min(minimum,pointer); maximum=std::max(maximum,pointer);
        }
        callbacks.push_back(work);
    }
    void Json(std::ostream& out, const P::NativePanelResult& panel) const {
        out << "{\"native_status\":" << static_cast<unsigned>(panel.Status)
            << ",\"left_code\":";
        if (panel.HasLeftPreflight) out << panel.LeftPreflight.OutcomeCode; else out << "null";
        out << ",\"right_code\":";
        if (panel.HasRightPreflight) out << panel.RightPreflight.OutcomeCode; else out << "null";
        out << ",\"handles_count\":" << count << ",\"handles_sha256\":"
            << Quote(Sha256Hex(bytes.data(),bytes.size()));
        const char* names[]={"handles_first","handles_last","handles_min","handles_max"};
        const uint64_t values[]={first,last,minimum,maximum};
        for (unsigned j=0;j<4;++j) {
            out << ',' << Quote(names[j]) << ':';
            if (count) out << values[j]; else out << "null";
        }
        out << '}';
    }
};

class Invocation : public P::NativePanelInvocation {
    Context& context_;
    const unsigned role_, scope_;
    Allocations& allocations_;
    const Budget& budget_;
    std::array<Encoder,kBatch> encoders_;
    std::array<Decoder,kBatch> decoders_;
    std::array<uint64_t,kBatch> pointers_ = {{0}};
    unsigned pointer_count_=0;
    int64_t prepare_code_=0;
    bool invoked_=false;
    void Remember(uint64_t pointer) { if (pointer) pointers_[pointer_count_++]=pointer; }
public:
    Invocation(Context& context, unsigned role, unsigned scope, Allocations& allocations, const Budget& budget)
        : context_(context), role_(role), scope_(scope), allocations_(allocations), budget_(budget) {
        budget_.Check(); context_.Prepare(role_,scope_);
        if (scope_ != EncoderCreate) prepare_code_=context_.create_codes[role_];
        if (scope_ >= ReceiveRecover) {
            prepare_code_=context_.caches[role_].code;
            // A fixed-prefix NeedMore is measured honestly for all64 messages,
            // then emitted as noncomparable, never extended to another ID.
            if (prepare_code_ == 1) prepare_code_=0;
        }
        if (scope_ >= Systematic && scope_ <= Distant)
            prepare_code_=context_.caches[role_].encode_codes[scope_-Systematic];
        if (scope_ == ReceiveRecover && prepare_code_ == 0) {
            for (unsigned i=0;i<kBatch;++i) {
                const int64_t code=decoders_[i].Create(context_.encoders[role_]);
                if (prepare_code_ == 0 && code != 0) prepare_code_=code;
            }
        }
    }
    std::string Identity() const override {
        return std::to_string(role_)+":"+std::to_string(scope_)+":"+std::to_string(context_.B)+":"+
            std::to_string(Address(context_.source.Data()));
    }
    P::NativePanelInvocationResult Invoke() override {
        Require(!invoked_, "callback invoked twice"); invoked_=true;
        budget_.Check();
        Counts work; work.callbacks=1;
        int64_t code=prepare_code_;
        unsigned extra=0;
        const bool receiving=scope_ >= ReceiveRecover;
        bool has_extra=false;
        const Clock::time_point start=Clock::now();
        if (scope_ == EncoderCreate) {
            for (unsigned i=0;i<kBatch;++i) {
                const int64_t current=encoders_[i].Create(role_,context_.source.Data(),context_.M,context_.B);
                ++work.encoder_creates;
                if (code == 0 && current != 0) code=current;
            }
        } else if (scope_ == DecoderCreate && code == 0) {
            for (unsigned i=0;i<kBatch;++i) {
                const int64_t current=decoders_[i].Create(context_.encoders[role_]);
                ++work.decoder_creates;
                if (code == 0 && current != 0) code=current;
            }
        } else if (scope_ >= Systematic && scope_ <= Distant && code == 0) {
            for (unsigned cycle=0;cycle<kBatch;++cycle) for (unsigned j=0;j<6;++j) {
                const uint32_t id=scope_ == Systematic ? j : scope_ == Sequential ? 6+j : kDistant[j];
                uint32_t length=UINT32_MAX;
                const int64_t current=context_.encoders[role_].Encode(id,
                    context_.packet.Data()+j*context_.stride,context_.B,length);
                ++work.encode_calls;
                if (code == 0 && (current != 0 || length != context_.B)) code=current == 0 ? -1002 : current;
            }
        } else if (receiving && code == 0) {
            for (unsigned i=0;i<kBatch;++i) {
                Decoder& decoder=decoders_[i];
                int64_t current=0;
                if (scope_ == DecodeEndToEnd) {
                    current=decoder.Create(context_.encoders[role_]); ++work.decoder_creates;
                    Remember(decoder.Pointer());
                }
                unsigned fed=0;
                if (current == 0) {
                    current=1;
                    for (unsigned j=0;j<10;++j) {
                        current=decoder.Feed(6+j,context_.packet.Data()+j*context_.stride,context_.B);
                        ++work.feed_calls; ++fed;
                        if (current != 1) break;
                    }
                    if (current == 0) {
                        if (fed < 6 || fed != context_.caches[role_].first_count) current=-1002;
                        else {
                            has_extra=true; extra=fed-6;
                            current=decoder.Recover(context_.output.Data()+i*context_.output_stride,context_.M);
                            ++work.recover_calls;
                        }
                    }
                }
                if (code == 0 && current != 0) code=current;
                if (scope_ == DecodeEndToEnd) decoder.Reset();
            }
        }
        const Clock::time_point stop=Clock::now();
        const auto signed_ns=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count();
        if (scope_ == EncoderCreate) for (const auto& encoder : encoders_) Remember(encoder.Pointer());
        else if (scope_ == DecoderCreate || scope_ == ReceiveRecover)
            for (const auto& decoder : decoders_) Remember(decoder.Pointer());
        else if (scope_ >= Systematic && scope_ <= Distant) Remember(context_.encoders[role_].Pointer());
        bool valid=context_.UnchangedSource();
        if (code == 0 && scope_ == EncoderCreate)
            for (const auto& encoder : encoders_) valid &= encoder.ValidProfile();
        if (code == 0 && scope_ >= Systematic && scope_ <= Distant)
            for (unsigned j=0;j<6;++j) {
                const auto& expected=context_.caches[role_].encode[scope_-Systematic][j];
                valid &= expected.size() == context_.B &&
                    std::equal(expected.begin(),expected.end(),context_.packet.Data()+j*context_.stride);
            }
        for (unsigned j=0;j<10;++j)
            valid &= Canary(context_.packet.Data()+j*context_.stride+context_.B,kGuard);
        if (receiving && (context_.caches[role_].code == 0 || context_.caches[role_].code == 1))
            for (unsigned j=0;j<10;++j)
                valid &= std::equal(context_.caches[role_].receive[j].begin(),
                    context_.caches[role_].receive[j].end(),context_.packet.Data()+j*context_.stride);
        if (code == 0 && receiving)
            for (unsigned i=0;i<kBatch;++i)
                valid &= std::equal(context_.original.begin(),context_.original.end(),
                                    context_.output.Data()+i*context_.output_stride);
        for (unsigned i=0;i<kBatch;++i)
            valid &= Canary(context_.output.Data()+i*context_.output_stride+context_.M,kGuard);
        allocations_.Record(pointers_,pointer_count_,work);
        budget_.Check();
        if (!valid || signed_ns <= 0 || code < 0 || (role_ == 0 && code != 0))
            return P::NativePanelInvocationResult(P::NativePanelDisposition::Fatal,
                !valid ? -1002 : code,has_extra,extra,0);
        return P::NativePanelInvocationResult(code == 0 ? P::NativePanelDisposition::Success :
            P::NativePanelDisposition::PreflightFailure,code,has_extra,extra,
            code == 0 ? static_cast<uint64_t>(signed_ns) : 0);
    }
};

struct Timing {
    std::ostringstream panels;
    unsigned count=0;
    bool noncomparable=false;
    Counts measured,warm;
};

void RunTiming(Timing& timing, const Oracle& oracle, const Budget& budget)
{
    std::array<std::array<std::unique_ptr<Context>,2>,3> contexts;
    for (unsigned width=0;width<3;++width) for (unsigned physical=0;physical<2;++physical) {
        budget.Check(); contexts[width][physical].reset(new Context(kWidths[width],oracle));
    }
    for (unsigned rep=0;rep<12;++rep) for (unsigned width=0;width<3;++width)
        for (unsigned scope=0;scope<7;++scope) for (unsigned comparison=0;comparison<5;++comparison)
            for (unsigned step=0;step<4;++step) {
                budget.Check();
                const unsigned condition=kConditions[rep%4][step];
                Allocations allocations;
                std::array<P::NativePanelArm,2> arms;
                for (unsigned logical=0;logical<2;++logical) {
                    Context& context=*contexts[width][logical^(condition/2)];
                    const unsigned role=kComparisons[comparison][logical];
                    const std::string identity=std::to_string(role)+":"+std::to_string(scope)+":"+
                        std::to_string(context.B)+":"+std::to_string(Address(context.source.Data()));
                    arms[logical]=P::NativePanelArm(identity,[&context,role,scope,&allocations,&budget]() {
                        return std::unique_ptr<P::NativePanelInvocation>(new Invocation(context,role,scope,allocations,budget));
                    });
                }
                const auto panel=P::ExecuteNativeTimingPanel(kCpu,condition%2 == 0 ?
                    P::NativePanelOrder::ABBA : P::NativePanelOrder::BAAB,2,arms[0],arms[1]);
                const unsigned status=panel.Status != P::NativePanelStatus::Complete ? 2 :
                    panel.PanelComparable ? 0 : 1;
                std::ostringstream out;
                out << '[' << rep << ',' << kWidths[width] << ',' << scope << ',' << comparison << ','
                    << condition << ',' << status << ',';
                if (panel.HasLeftPreflight && panel.LeftPreflight.HasDecodedExtra) out << panel.LeftPreflight.DecodedExtra;
                else out << "null";
                out << ',';
                if (panel.HasRightPreflight && panel.RightPreflight.HasDecodedExtra) out << panel.RightPreflight.DecodedExtra;
                else out << "null";
                out << ',';
                if (status == 0) {
                    out << '[';
                    for (unsigned i=0;i<8;++i) {
                        Require(panel.Slots[i].HasElapsedNanoseconds && panel.Slots[i].ElapsedNanoseconds,
                                "missing comparable panel interval");
                        if (i) out << ',';
                        out << panel.Slots[i].ElapsedNanoseconds;
                    }
                    out << ']';
                } else out << "null";
                out << ",[";
                for (unsigned physical=0;physical<2;++physical) {
                    const Context& context=*contexts[width][physical];
                    if (physical) out << ',';
                    out << Address(context.source.Data()) << ',' << Address(context.packet.Data()) << ','
                        << Address(context.output.Data());
                }
                out << "],"; allocations.Json(out,panel); out << ']';
                if (timing.count++) timing.panels << ',';
                timing.panels << out.str();
                if (status == 2) throw std::runtime_error("native panel invalid: "+panel.Diagnostic);
                Require(allocations.callbacks.size() == 10, "complete callback roster");
                for (size_t i=0;i<10;++i) (i < 2 ? timing.warm : timing.measured).Add(allocations.callbacks[i]);
                timing.noncomparable |= status == 1;
            }
    Require(timing.count == 5040 && timing.measured.callbacks == 40320 && timing.warm.callbacks == 10080,
            "complete native panel accounting");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2 || std::strcmp(argv[1],"--worker") != 0) return 2;
    const Budget budget;
    Correctness correctness;
    Timing timing;
    std::string outcome="INVALID", diagnostic, visits, before="null", after="null";
    gf256_x86_cpu_features features = {};
    bool initialized=false;
    try {
        struct rlimit limit = {512u*1024u*1024u,512u*1024u*1024u};
        Require(setrlimit(RLIMIT_AS,&limit) == 0, "worker address-space limit");
        Pin(); budget.Check(); before=Identity();
        Require(wirehair_init() == Wirehair_Success && gf256_init() == 0, "shared GF initialization");
        initialized=true;
        Require(GF256Ctx.Polynomial == 0x14d, "ordinary GF256 polynomial");
        gf256_get_active_x86_cpu_features(&features);
        Require(sizeof(D::kLookup) == 39936 && Sha256Hex(D::kLookup,sizeof(D::kLookup)) == D::kLookupSha256,
                "immutable lookup hash");
        const Oracle oracle;
        visits=VisitHash(oracle,budget);
        RunCorrectness(correctness,oracle,budget);
        budget.Check(); Require(Identity() == before, "post-correctness target identity");
        RunTiming(timing,oracle,budget);
        budget.Check(); after=Identity(); Require(after == before, "post-timing target identity");
        gf256_x86_cpu_features final_features = {};
        gf256_get_active_x86_cpu_features(&final_features);
        Require(final_features.SSSE3 == features.SSSE3 && final_features.AVX2 == features.AVX2 &&
            final_features.GFNI == features.GFNI && final_features.AVX512 == features.AVX512 &&
            GF256Ctx.Polynomial == 0x14d, "post-timing shared GF runtime");
        Require(Sha256Hex(D::kLookup,sizeof(D::kLookup)) == D::kLookupSha256, "post-timing immutable lookup");
        outcome=correctness.noncomparable || timing.noncomparable ? "NONCOMPARABLE" : "COMPLETE";
    } catch (const std::exception& error) { diagnostic=error.what(); }
      catch (...) { diagnostic="unknown worker failure"; }
    std::ostringstream out;
    out << "{\"schema\":\"wirehair.wh2.thue-morse-native-r0.raw.v1\",\"protocol\":" << Quote(kProtocol)
        << ",\"target_cpu\":" << kCpu << ",\"gf_runtime\":{\"polynomial\":"
        << (initialized ? GF256Ctx.Polynomial : 0) << ",\"address\":" << Address(&GF256Ctx)
        << ",\"ssse3\":" << features.SSSE3 << ",\"avx2\":" << features.AVX2
        << ",\"gfni\":" << features.GFNI << ",\"avx512\":" << features.AVX512 << '}'
        << ",\"data_sha256\":" << Quote(D::kDataSha256) << ",\"visit_sha256\":" << Quote(visits)
        << ",\"target_identity_before\":" << before << ",\"target_identity_after\":" << after
        << ",\"config\":{\"K\":6,\"widths\":[2,64,1280],\"batch_messages\":64,\"batch_encode_calls\":384,\"scopes\":[";
    for (unsigned i=0;i<7;++i) { if (i) out << ','; out << Quote(kScopes[i]); }
    out << "],\"replicates\":12,\"conditions\":4,\"comparisons\":[[\"T\",\"T\"],[\"L\",\"L\"],[\"P\",\"P\"],[\"T\",\"L\"],[\"T\",\"P\"]],\"invocations_per_slot\":2"
        << ",\"condition_orders\":[[0,1,3,2],[1,2,0,3],[2,3,1,0],[3,0,2,1]]"
        << ",\"systematic_ids\":[0,1,2,3,4,5],\"sequential_ids\":[6,7,8,9,10,11]"
        << ",\"distant_ids\":[4294967295,4294967293,4294967291,4294967289,4294967287,4294967285]"
        << ",\"receive_ids\":[6,7,8,9,10,11,12,13,14,15],\"warmup_order\":[\"left\",\"right\"]"
        << ",\"lookup_bytes\":39936,\"lookup_sha256\":" << Quote(D::kLookupSha256)
        << ",\"message_bytes\":[12,384,7680],\"source_seed\":\"0x6e61746976653633\""
        << ",\"source_law\":\"splitmix64_le_bytes_seed_xor_M_xor_B_shift32\"}"
        << ",\"correctness\":{\"status\":" << correctness.status << ",\"cases\":" << correctness.cases
        << ",\"candidate_packets\":" << correctness.candidate_packets << ",\"control_packets\":" << correctness.control_packets
        << ",\"recovered_messages\":" << correctness.recovered_messages << ",\"records\":[" << correctness.records.str() << "]}"
        << ",\"panels\":[" << timing.panels.str() << "],\"timed_scope_counts\":{\"measured\":";
    timing.measured.Json(out); out << ",\"warm\":"; timing.warm.Json(out);
    out << "},\"elapsed_seconds\":" << std::setprecision(17) << budget.Elapsed()
        << ",\"outcome\":" << Quote(outcome) << ",\"diagnostic\":" << Quote(diagnostic) << "}\n";
    const std::string raw=out.str();
    if (raw.size() > 8u*1024u*1024u) return 1;
    std::cout << raw; std::cout.flush();
    return std::cout.good() && outcome != "INVALID" ? 0 : 1;
}
