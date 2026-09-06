// Frozen .65 row-only cost screen. This is not a production codec.
#include "Wh2ThueMorseMulRowsR0.h"
#include "Wh2ThueMorseNativeDataR0.h"
#include "Wh2FrozenTrace.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "gf256.h"

#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
#include <unistd.h>

namespace {
namespace N = wh2_thue_native;
namespace C = wh2_thue_mulrows_r0;
namespace D = wh2_thue_native_data;
namespace P = wirehair_wh2_bench;
using wirehair::wh2_benchmark::Sha256Hex;
using Clock = std::chrono::steady_clock;
using RowFunction = N::Status (*)(N::Lookup, uint32_t, uint8_t*);
using Ids = std::array<uint32_t, 64>;
using Rows = std::array<std::array<uint8_t, 6>, 64>;
static const int kCpu = 50;
static const char kProtocol[] = "wirehair.wh2.thue-morse-mulrows-r0";
static const uint64_t kChecksumSeed = UINT64_C(1469598103934665603);
static const uint64_t kChecksumFactor = UINT64_C(1099511628211);
static const unsigned kComparison[3][2] = {{0,0}, {1,1}, {0,1}};
static const unsigned kSlots[8] = {0,1,1,0,1,0,0,1};

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

std::string Hex(const std::string& bytes)
{
    std::ostringstream out;
    for (unsigned char byte : bytes)
        out << std::hex << std::setfill('0') << std::setw(2) << static_cast<unsigned>(byte);
    return out.str();
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

// Same canonical identity representation as the immutable native R0 worker.
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
    const auto& raw=identity.Raw;
    const auto& derived=identity.Derived;
    Require(derived.Family == 26 && derived.Model == 8 && derived.Stepping == 1 &&
            derived.FullApicId == 100 && derived.CoreId == 50 && derived.PackageId == 0 &&
            derived.ThreadId == 0 && derived.ThreadsPerCore == 2 && derived.CcdId == 6 &&
            derived.ComplexId == 6 && derived.LogicalProcessorsPerPackage == 128,
            "frozen physical target identity");
    std::ostringstream out;
    out << "{\"after_cpu\":" << identity.After.Cpu << ",\"before_cpu\":" << identity.Before.Cpu
        << ",\"canonical_bytes\":" << canonical.size() << ",\"canonical_hex\":" << Quote(Hex(canonical))
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

std::string Runtime()
{
    gf256_x86_cpu_features features = {};
    gf256_get_active_x86_cpu_features(&features);
    Require(GF256Ctx.Polynomial == 0x14d, "ordinary GF256 polynomial");
    Require(features.SSSE3 == 1 && features.AVX2 == 1 && features.GFNI == 1 && features.AVX512 == 1,
            "frozen shared GF dispatch");
    std::ostringstream out;
    out << "{\"polynomial\":" << GF256Ctx.Polynomial << ",\"address\":" << Address(&GF256Ctx)
        << ",\"ssse3\":" << features.SSSE3 << ",\"avx2\":" << features.AVX2
        << ",\"gfni\":" << features.GFNI << ",\"avx512\":" << features.AVX512 << '}';
    return out.str();
}

struct Budget {
    Clock::time_point started = Clock::now();
    void Check() const {
        Require(Clock::now() - started < std::chrono::seconds(30), "worker wall deadline");
        Require(OnCpu(), "worker CPU drift");
    }
    uint64_t Elapsed() const {
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-started).count();
        Require(elapsed > 0 && elapsed < INT64_C(30000000000), "worker elapsed range");
        return static_cast<uint64_t>(elapsed);
    }
};

Ids Family(unsigned family)
{
    Require(family < 4, "invalid family");
    Ids ids;
    for (uint32_t j=0; j<64; ++j) {
        uint32_t id=(13*j)&1023;
        if (family >= 1) id |= (1+(17*j)%127)<<10;
        if (family >= 2) id |= (1+(23*j)%127)<<17;
        if (family >= 3) id |= (1+(37*j)%255)<<24;
        ids[j]=id;
    }
    return ids;
}

struct Slot {
    uint8_t before[16];
    uint8_t row[6];
    uint8_t after[16];
};
static_assert(sizeof(Slot) == 38, "fixed guarded row stride");
struct Output {
    Slot slots[64];
    void Reset() { std::memset(slots,0xa5,sizeof(slots)); }
    bool Matches(const Rows& expected) const {
        for (unsigned j=0; j<64; ++j) {
            if (std::memcmp(slots[j].row,expected[j].data(),6) != 0) return false;
            for (unsigned k=0; k<16; ++k)
                if (slots[j].before[k] != 0xa5 || slots[j].after[k] != 0xa5) return false;
        }
        return true;
    }
};

uint64_t Consume(uint64_t checksum, const Output& output)
{
    for (unsigned j=0; j<64; ++j)
        for (unsigned k=0; k<6; ++k)
            checksum=(checksum*kChecksumFactor)^output.slots[j].row[k];
    return checksum;
}

struct CallResult {
    uint32_t calls;
    bool success;
};

// One out-of-line, non-specializable indirect callsite for both Row functions.
// Status checks are intentionally inside the interval; row/guard checks are not.
__attribute__((noinline,noclone,noipa))
CallResult RunRows(RowFunction function, N::Lookup lookup, const Ids& ids, Output& output)
{
    CallResult result = {0,true};
    for (unsigned cycle=0; cycle<64; ++cycle)
        for (unsigned j=0; j<64; ++j) {
            const N::Status status=function(lookup,ids[j],output.slots[j].row);
            result.success &= status == N::Status::Success;
            ++result.calls;
        }
    return result;
}

template<class Callback>
std::array<uint64_t,8> Panel(RowFunction left, RowFunction right, unsigned order, Callback& callback)
{
    Require(left && right && order < 2, "invalid panel arguments");
    Require(callback(left,true) > 0, "invalid warm-left interval");
    Require(callback(right,true) > 0, "invalid warm-right interval");
    const RowFunction functions[2] = {left,right};
    std::array<uint64_t,8> times;
    for (unsigned slot=0; slot<8; ++slot) {
        times[slot]=callback(functions[kSlots[slot]^order],false);
        Require(times[slot] > 0 && times[slot] <= static_cast<uint64_t>(INT64_MAX),
                "invalid measured interval");
    }
    return times;
}

// Independent carryless product and polynomial long division, without GF tables.
uint8_t Product(unsigned a, unsigned b)
{
    unsigned polynomial=0;
    for (unsigned bit=0; bit<8; ++bit)
        if ((b>>bit)&1u) polynomial ^= a<<bit;
    for (int bit=14; bit>=8; --bit)
        if ((polynomial>>bit)&1u) polynomial ^= 0x14du<<(bit-8);
    return static_cast<uint8_t>(polynomial);
}

unsigned Parity(uint32_t value)
{
    unsigned result=0;
    for (; value; value>>=1) result ^= value&1;
    return result;
}

std::array<uint8_t,6> Oracle(N::Lookup lookup, uint32_t id)
{
    const unsigned low=id&1023, middle10=(id>>10)&127;
    const unsigned middle17=(id>>17)&127, high=id>>24;
    const unsigned phase17=Parity(high), phase10=phase17^Parity(middle17);
    const unsigned phase0=phase10^Parity(middle10);
    std::array<uint8_t,6> row;
    std::memcpy(row.data(),lookup.data+phase0*6144+low*6,6);
    const unsigned values[3]={middle10,middle17,high};
    const size_t offsets[3]={12288+phase10*4608+middle10*36,
        21504+phase17*4608+middle17*36,30720+high*36};
    for (unsigned group=0; group<3; ++group) if (values[group]) {
        std::array<uint8_t,6> next = {{0}};
        for (unsigned r=0; r<6; ++r)
            for (unsigned c=0; c<6; ++c)
                next[r] ^= Product(lookup.data[offsets[group]+r*6+c],row[c]);
        row=next;
    }
    return row;
}

enum class Mode { Invalid, Worker, Selftest };
Mode Arguments(int argc, const char* const* argv)
{
    if (argc != 2 || !argv || !argv[1]) return Mode::Invalid;
    if (std::strcmp(argv[1],"--worker") == 0) return Mode::Worker;
    if (std::strcmp(argv[1],"--selftest") == 0) return Mode::Selftest;
    return Mode::Invalid;
}

// These callbacks do not read a lookup or initialize/use the GF runtime.
uint64_t mock_calls=0;
N::Status MockA(N::Lookup, uint32_t id, uint8_t* output)
{
    ++mock_calls;
    for (unsigned k=0; k<6; ++k) output[k]=static_cast<uint8_t>((id>>(k%4)*8)^k);
    return N::Status::Success;
}
N::Status MockB(N::Lookup lookup, uint32_t id, uint8_t* output)
{
    return MockA(lookup,id,output);
}
N::Status MockFailure(N::Lookup, uint32_t, uint8_t*)
{
    ++mock_calls;
    return N::Status::InvalidInput;
}

int Selftest()
{
    const char* valid[]={"test","--selftest"};
    const char* worker[]={"test","--worker"};
    const char* unknown[]={"test","--help"};
    const char* missing[]={"test",nullptr};
    Require(Arguments(2,valid) == Mode::Selftest && Arguments(2,worker) == Mode::Worker &&
            Arguments(2,unknown) == Mode::Invalid && Arguments(1,valid) == Mode::Invalid &&
            Arguments(3,valid) == Mode::Invalid && Arguments(2,missing) == Mode::Invalid &&
            Arguments(2,nullptr) == Mode::Invalid, "argument selftest");
    for (unsigned family=0; family<4; ++family) {
        const Ids ids=Family(family);
        for (unsigned j=0; j<64; ++j) {
            Require((ids[j]&1023) == 13*j, "low ID roster");
            Require(((ids[j]>>10)&127) == (family>=1 ? 1+(17*j)%127 : 0), "middle10 ID roster");
            Require(((ids[j]>>17)&127) == (family>=2 ? 1+(23*j)%127 : 0), "middle17 ID roster");
            Require((ids[j]>>24) == (family>=3 ? 1+(37*j)%255 : 0), "high ID roster");
            for (unsigned k=0; k<j; ++k) Require(ids[k] != ids[j], "unique family IDs");
        }
    }
    bool rejected=false;
    try { Family(4); } catch (const std::runtime_error&) { rejected=true; }
    Require(rejected,"invalid family selftest");
    const unsigned expected_slots[2][8]={{0,1,1,0,1,0,0,1},{1,0,0,1,0,1,1,0}};
    uint64_t panels=0,warm=0,measured=0;
    for (unsigned rep=0; rep<12; ++rep)
        for (unsigned family=0; family<4; ++family)
            for (unsigned comparison=0; comparison<3; ++comparison)
                for (unsigned step=0; step<2; ++step) {
                    const unsigned order=(rep+family+comparison+step)%2;
                    unsigned index=0;
                    const RowFunction functions[2]={MockA,MockB};
                    auto callback=[&](RowFunction fn,bool is_warm)->uint64_t {
                        Require(index < 10,"excess callback");
                        const unsigned side=index<2 ? index : expected_slots[order][index-2];
                        Require(fn == functions[kComparison[comparison][side]], "callback side order");
                        Require(is_warm == (index<2), "warm callback order");
                        ++index; (is_warm ? warm : measured)++;
                        return index;
                    };
                    const auto times=Panel(functions[kComparison[comparison][0]],
                        functions[kComparison[comparison][1]],order,callback);
                    Require(index == 10 && times.front() == 3 && times.back() == 10, "panel count selftest");
                    ++panels;
                }
    Require(panels == 288 && warm == 576 && measured == 2304, "full neutral panel roster");
    auto unused=[](RowFunction,bool)->uint64_t { throw std::runtime_error("unexpected callback"); };
    rejected=false;
    try { Panel(MockA,MockB,2,unused); } catch (const std::runtime_error& e) {
        rejected=std::string(e.what()) == "invalid panel arguments";
    }
    Require(rejected,"invalid order selftest");
    rejected=false;
    try { Panel(nullptr,MockB,0,unused); } catch (const std::runtime_error& e) {
        rejected=std::string(e.what()) == "invalid panel arguments";
    }
    Require(rejected,"null callback selftest");
    auto zero=[](RowFunction,bool)->uint64_t { return 0; };
    rejected=false;
    try { Panel(MockA,MockB,0,zero); } catch (const std::runtime_error& e) {
        rejected=std::string(e.what()) == "invalid warm-left interval";
    }
    Require(rejected,"zero warm interval selftest");
    auto invalid=[](RowFunction,bool warm)->uint64_t {
        return warm ? 1 : static_cast<uint64_t>(INT64_MAX)+1;
    };
    rejected=false;
    try { Panel(MockA,MockB,0,invalid); } catch (const std::runtime_error& e) {
        rejected=std::string(e.what()) == "invalid measured interval";
    }
    Require(rejected,"overflowed measured interval selftest");
    Output output;
    const Ids ids=Family(3);
    Rows expected;
    for (unsigned j=0; j<64; ++j)
        for (unsigned k=0; k<6; ++k) expected[j][k]=static_cast<uint8_t>((ids[j]>>(k%4)*8)^k);
    const N::Lookup neutral={nullptr,0};
    for (RowFunction function : {MockA,MockB}) {
        output.Reset(); mock_calls=0;
        const CallResult result=RunRows(function,neutral,ids,output);
        Require(result.calls == 4096 && result.success && mock_calls == 4096 && output.Matches(expected),
                "mock call/output accounting");
    }
    const uint64_t checksum=Consume(kChecksumSeed,output);
    uint64_t independent=kChecksumSeed;
    for (const auto& row : expected)
        for (uint8_t byte : row) independent=(independent*kChecksumFactor)^byte;
    Require(checksum == independent,"checksum selftest");
    output.slots[63].after[15]^=1;
    Require(!output.Matches(expected),"trailing guard selftest");
    output.slots[63].after[15]^=1; output.slots[0].before[0]^=1;
    Require(!output.Matches(expected),"leading guard selftest");
    output.slots[0].before[0]^=1; output.slots[32].row[5]^=1;
    Require(!output.Matches(expected),"row corruption selftest");
    output.Reset(); mock_calls=0;
    const CallResult failed=RunRows(MockFailure,neutral,ids,output);
    Require(!failed.success && failed.calls == 4096 && mock_calls == 4096,"failed status complete roster");
    std::cout << "MULROWS R0 neutral selftest PASS\n";
    std::cout.flush();
    return std::cout.good() ? 0 : 1;
}

int Worker()
{
    const Budget budget;
    // No limits, affinity operations, clocks, or fixed-data access in Selftest.
    const struct rlimit memory={512u*1024u*1024u,512u*1024u*1024u};
    const struct rlimit cpu={25,25};
    Require(setrlimit(RLIMIT_AS,&memory) == 0 && setrlimit(RLIMIT_CPU,&cpu) == 0,"worker resource limits");
    alarm(30);
    Pin(); budget.Check();
    const std::string identity_before=Identity();
    Require(gf256_init() == 0,"shared GF initialization");
    const std::string runtime_before=Runtime();
    const N::Lookup lookup={D::kLookup,sizeof(D::kLookup)};
    Require(lookup.bytes == 39936 && Sha256Hex(lookup.data,lookup.bytes) == D::kLookupSha256,
            "immutable lookup hash");
    std::array<Ids,4> families;
    std::array<Rows,4> expected;
    std::string expected_bytes;
    uint64_t preflight=0;
    for (unsigned family=0; family<4; ++family) {
        families[family]=Family(family);
        for (unsigned j=0; j<64; ++j) {
            const auto oracle=Oracle(lookup,families[family][j]);
            Slot baseline;
            std::memset(&baseline,0xa5,sizeof(baseline));
            Require(N::Row(lookup,families[family][j],baseline.row) == N::Status::Success,
                    "baseline oracle preflight status");
            ++preflight;
            Require(std::memcmp(baseline.row,oracle.data(),6) == 0,"baseline scalar oracle mismatch");
            for (unsigned k=0; k<16; ++k)
                Require(baseline.before[k] == 0xa5 && baseline.after[k] == 0xa5,"baseline preflight guard");
            expected[family][j]=oracle;
            expected_bytes.append(reinterpret_cast<const char*>(oracle.data()),6);
        }
        budget.Check();
    }
    Require(preflight == 256 && expected_bytes.size() == 1536,"oracle preflight count");
    const std::string expected_hash=Sha256Hex(expected_bytes);
    Output output;
    const RowFunction functions[2]={N::Row,C::Row};
    uint64_t warm=0,measured=0,row_calls=0,checked_rows=0,checksum=kChecksumSeed;
    unsigned panel_count=0;
    std::ostringstream panels;
    for (unsigned rep=0; rep<12; ++rep)
        for (unsigned family=0; family<4; ++family)
            for (unsigned comparison=0; comparison<3; ++comparison)
                for (unsigned step=0; step<2; ++step) {
                    const unsigned order=(rep+family+comparison+step)%2;
                    auto callback=[&](RowFunction function,bool is_warm)->uint64_t {
                        budget.Check(); output.Reset();
                        std::atomic_signal_fence(std::memory_order_seq_cst);
                        const auto start=Clock::now();
                        const CallResult result=RunRows(function,lookup,families[family],output);
                        const auto stop=Clock::now();
                        std::atomic_signal_fence(std::memory_order_seq_cst);
                        budget.Check();
                        Require(result.success && result.calls == 4096,"Row status/count mismatch");
                        Require(output.Matches(expected[family]),"callback row/guard mismatch");
                        checksum=Consume(checksum,output);
                        row_calls+=result.calls; checked_rows+=64;
                        (is_warm ? warm : measured)++;
                        const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count();
                        Require(elapsed > 0,"nonpositive callback interval");
                        return static_cast<uint64_t>(elapsed);
                    };
                    const auto times=Panel(functions[kComparison[comparison][0]],
                        functions[kComparison[comparison][1]],order,callback);
                    if (panel_count++) panels << ',';
                    panels << '[' << rep << ',' << family << ',' << comparison << ',' << order << ",[";
                    for (unsigned i=0; i<8; ++i) { if (i) panels << ','; panels << times[i]; }
                    panels << "]]";
                }
    Require(panel_count == 288 && warm == 576 && measured == 2304 &&
            row_calls == UINT64_C(11796480) && checked_rows == 184320,"complete worker accounting");
    budget.Check();
    const std::string identity_after=Identity(),runtime_after=Runtime();
    Require(identity_after == identity_before && runtime_after == runtime_before,"post-worker identity/runtime drift");
    Require(Sha256Hex(lookup.data,lookup.bytes) == D::kLookupSha256,"post-worker lookup drift");
    std::ostringstream out;
    out << "{\"schema\":\"wirehair.wh2.thue-morse-mulrows-r0.raw.v1\",\"protocol\":" << Quote(kProtocol)
        << ",\"outcome\":\"COMPLETE\",\"target_cpu\":" << kCpu
        << ",\"lookup_sha256\":" << Quote(D::kLookupSha256)
        << ",\"panels\":[" << panels.str() << "],\"warm_callbacks\":" << warm
        << ",\"measured_callbacks\":" << measured << ",\"row_calls\":" << row_calls
        << ",\"calls_per_callback\":4096,\"preflight_row_calls\":" << preflight
        << ",\"checked_rows\":" << checked_rows << ",\"checksum\":" << checksum
        << ",\"expected_rows_hex\":" << Quote(Hex(expected_bytes))
        << ",\"expected_rows_sha256\":" << Quote(expected_hash)
        << ",\"output_address\":" << Address(output.slots[0].row) << ",\"output_stride\":" << sizeof(Slot)
        << ",\"baseline_function_address\":" << static_cast<uint64_t>(reinterpret_cast<uintptr_t>(functions[0]))
        << ",\"candidate_function_address\":" << static_cast<uint64_t>(reinterpret_cast<uintptr_t>(functions[1]))
        << ",\"identity_before\":" << identity_before << ",\"identity_after\":" << identity_after
        << ",\"runtime_before\":" << runtime_before << ",\"runtime_after\":" << runtime_after
        << ",\"elapsed_ns\":" << budget.Elapsed() << "}\n";
    const std::string raw=out.str();
    Require(raw.size() <= 2u*1024u*1024u,"raw size cap");
    budget.Check();
    std::cout << raw; std::cout.flush();
    return std::cout.good() ? 0 : 1;
}

} // namespace

int main(int argc, char** argv)
{
    const Mode mode=Arguments(argc,argv);
    if (mode == Mode::Invalid) return 2;
    try { return mode == Mode::Selftest ? Selftest() : Worker(); }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
    catch (...) { std::cerr << "unknown worker failure\n"; return 1; }
}
