// Frozen .64 placement diagnostic. No production dispatch or equation changes.
#include "Wh2RepairAlignmentR0Bridge.h"
#include "Wh2FrozenTrace.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "gf256.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

#if defined(__GNUC__) && !defined(__clang__)
#define WH2_ALIGNMENT_NOINLINE __attribute__((noinline, noipa))
#elif defined(__clang__)
#define WH2_ALIGNMENT_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define WH2_ALIGNMENT_NOINLINE __declspec(noinline)
#else
#define WH2_ALIGNMENT_NOINLINE
#endif

namespace alignment {
namespace B = wh2_repair_alignment_r0;
namespace P = wirehair_wh2_bench;
using wirehair::wh2_benchmark::Sha256Hex;
static const char kProtocol[] = "wirehair.wh2.repair-alignment-r0";
static const uint64_t kMissing=UINT64_MAX, kWallLimit=UINT64_C(10000000000);
static const uint64_t kEncodeLimit=UINT64_C(200000000);
static const uint64_t kHashSeed=UINT64_C(1469598103934665603), kHashFactor=UINT64_C(1099511628211);
static const uint64_t kPreludeSeed=UINT64_C(0x9e3779b97f4a7c15);
static const uint64_t kPreludeFinal=UINT64_C(0x43935dad1647741b);
static const unsigned kPreludeCount=1u<<20, kBlock=1280, kMessage=7680, kIntermediate=46080;
static const unsigned kCarrier=65536, kOutput=12288, kStride=1408, kScratch=492544;
static const unsigned kSlots[8]={0,1,1,0,1,0,0,1};
static const unsigned kNatural[8][2]={{0,0},{1,1},{2,2},{3,3},{0,2},{1,3},{0,1},{2,3}};
static const unsigned kShadow[6][2]={{0,0},{1,1},{2,2},{3,3},{0,1},{2,3}};
using Counters=std::array<uint64_t,4>;
using Packets=std::array<uint8_t,kMessage>;

void Require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
uint64_t Address(const void* pointer) { return static_cast<uint64_t>(reinterpret_cast<uintptr_t>(pointer)); }
std::string Quote(const std::string& value) {
    std::ostringstream out; out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c >= 32 && c < 127) out << c;
        else out << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << unsigned(c) << std::dec;
    }
    out << '"'; return out.str();
}
std::string Hex(const void* input, size_t bytes) {
    const uint8_t* data=static_cast<const uint8_t*>(input);
    static const char digits[]="0123456789abcdef";
    std::string result(bytes*2,'0');
    for (size_t i=0;i<bytes;++i) { result[i*2]=digits[data[i]>>4]; result[i*2+1]=digits[data[i]&15]; }
    return result;
}
void Number(std::ostream& out, uint64_t value) { if (value == kMissing) out << "null"; else out << value; }
template<size_t N> void Numbers(std::ostream& out, const std::array<uint64_t,N>& values) {
    out << '['; for (size_t i=0;i<N;++i) { if (i) out << ','; Number(out,values[i]); } out << ']';
}
std::string PointerJson(const void* pointer) {
    const uint64_t value=Address(pointer); std::ostringstream out;
    out << "{\"address\":" << value << ",\"mod32\":" << value%32
        << ",\"mod64\":" << value%64 << ",\"mod4096\":" << value%4096 << '}'; return out.str();
}
bool OnCpu() {
    cpu_set_t mask; CPU_ZERO(&mask);
    return sched_getaffinity(0,sizeof(mask),&mask) == 0 && CPU_COUNT(&mask) == 1 &&
        CPU_ISSET(50,&mask) && sched_getcpu() == 50;
}
void Pin() {
    cpu_set_t mask; CPU_ZERO(&mask); CPU_SET(50,&mask);
    Require(sched_setaffinity(0,sizeof(mask),&mask) == 0 && OnCpu(),"singleton CPU50 pin");
}
std::string Identity() {
    P::TargetIdentityReceiptV2 identity; std::string diagnostic,canonical;
    Require(P::CapturePublicBorrowedTargetIdentity(50,identity,diagnostic),"target identity capture");
    Require(P::SerializeTargetIdentityV2(identity,canonical,diagnostic) && !canonical.empty() &&
        canonical.size() <= 4096 && Sha256Hex(canonical) == identity.CanonicalSha256 &&
        identity.Before.Affinity == std::vector<int32_t>(1,50) &&
        identity.After.Affinity == std::vector<int32_t>(1,50),"target identity serialization");
    const auto& r=identity.Raw; const auto& d=identity.Derived;
    Require(d.Family == 26 && d.Model == 8 && d.Stepping == 1 && d.FullApicId == 100 &&
        d.CoreId == 50 && d.PackageId == 0 && d.ThreadId == 0 && d.ThreadsPerCore == 2 &&
        d.CcdId == 6 && d.ComplexId == 6 && d.LogicalProcessorsPerPackage == 128,"frozen physical target");
    std::ostringstream out;
    out << "{\"after_cpu\":" << identity.After.Cpu << ",\"before_cpu\":" << identity.Before.Cpu
        << ",\"canonical_bytes\":" << canonical.size() << ",\"canonical_hex\":" << Quote(Hex(canonical.data(),canonical.size()))
        << ",\"canonical_sha256\":" << Quote(identity.CanonicalSha256)
        << ",\"capabilities\":{\"leaf1_ecx\":" << r.Leaf1.Regs.Ecx << ",\"leaf1_edx\":" << r.Leaf1.Regs.Edx
        << ",\"leaf6_eax\":" << r.Leaf6.Regs.Eax << ",\"leaf6_ecx\":" << r.Leaf6.Regs.Ecx
        << ",\"leaf80000001_ecx\":" << r.Leaf80000001.Regs.Ecx << ",\"leaf80000001_edx\":" << r.Leaf80000001.Regs.Edx
        << ",\"leaf80000008_ebx\":" << r.Leaf80000008.Regs.Ebx << ",\"leaf80000021_eax\":" << r.Leaf80000021.Regs.Eax
        << ",\"max_basic_leaf\":" << r.Leaf0.Regs.Eax << ",\"max_extended_leaf\":" << r.Leaf80000000.Regs.Eax
        << "},\"derived\":{\"ccd_id\":" << d.CcdId << ",\"complex_id\":" << d.ComplexId << ",\"core_id\":" << d.CoreId
        << ",\"family\":" << d.Family << ",\"full_apic_id\":" << d.FullApicId << ",\"initial_apic_id_8\":" << d.InitialApicId8
        << ",\"logical_processors_per_package\":" << d.LogicalProcessorsPerPackage << ",\"model\":" << d.Model
        << ",\"package_id\":" << d.PackageId << ",\"stepping\":" << d.Stepping << ",\"thread_id\":" << d.ThreadId
        << ",\"threads_per_core\":" << d.ThreadsPerCore << "},\"raw_capture_complete\":" << (identity.RawCaptureComplete?"true":"false")
        << ",\"requested_cpu\":" << identity.RequestedCpu << ",\"semantic_validation_passed\":"
        << (identity.SemanticValidationPassed?"true":"false") << ",\"singleton_affinity_verified\":true}";
    return out.str();
}
std::string Runtime() {
    gf256_x86_cpu_features f={}; gf256_get_active_x86_cpu_features(&f);
    Require(GF256Ctx.Polynomial == 0x14d && f.SSSE3 == 1 && f.AVX2 == 1 && f.GFNI == 1 && f.AVX512 == 1,
        "frozen ordinary GF256 dispatch");
    std::ostringstream out; out << "{\"polynomial\":" << GF256Ctx.Polynomial << ",\"address\":" << Address(&GF256Ctx)
        << ",\"ssse3\":" << f.SSSE3 << ",\"avx2\":" << f.AVX2 << ",\"gfni\":" << f.GFNI << ",\"avx512\":" << f.AVX512 << '}';
    return out.str();
}
void Touch(void* data,size_t bytes) {
    volatile uint8_t* p=static_cast<volatile uint8_t*>(data);
    for (size_t i=0;i<bytes;++i) { const uint8_t value=p[i]; p[i]=value; }
}
bool Nanoseconds(const timespec& t,uint64_t& value) {
    if (t.tv_sec < 0 || t.tv_nsec < 0 || t.tv_nsec >= 1000000000L ||
        uint64_t(t.tv_sec) > (kMissing-1)/UINT64_C(1000000000)) return false;
    const uint64_t base=uint64_t(t.tv_sec)*UINT64_C(1000000000);
    if (uint64_t(t.tv_nsec) >= kMissing-base) return false;
    value=base+uint64_t(t.tv_nsec); return true;
}
struct Reader {
    timespec clock_scratch={}; rusage usage_scratch={};
    void Prepare() { Touch(&clock_scratch,sizeof(clock_scratch)); Touch(&usage_scratch,sizeof(usage_scratch)); }
    bool Clock(clockid_t id,uint64_t& value) { return clock_gettime(id,&clock_scratch) == 0 && Nanoseconds(clock_scratch,value); }
    bool Mono(uint64_t& value) { return Clock(CLOCK_MONOTONIC,value); }
    bool Cpu(uint64_t& value) { return Clock(CLOCK_THREAD_CPUTIME_ID,value); }
    bool Usage(Counters& value) {
        if (getrusage(RUSAGE_THREAD,&usage_scratch) || usage_scratch.ru_minflt < 0 || usage_scratch.ru_majflt < 0 ||
            usage_scratch.ru_nvcsw < 0 || usage_scratch.ru_nivcsw < 0) return false;
        value={{uint64_t(usage_scratch.ru_minflt),uint64_t(usage_scratch.ru_majflt),
                uint64_t(usage_scratch.ru_nvcsw),uint64_t(usage_scratch.ru_nivcsw)}}; return true;
    }
    uint64_t Now() { uint64_t value=kMissing; Require(Mono(value),"worker clock"); return value; }
    uint64_t Resolution(clockid_t id) {
        uint64_t value=kMissing;
        Require(clock_getres(id,&clock_scratch) == 0 && Nanoseconds(clock_scratch,value) && value > 0,"clock resolution");
        return value;
    }
};
struct Observation {
    uint64_t m0=kMissing,c0=kMissing,m1=kMissing,m2=kMissing,c1=kMissing,m3=kMissing;
    Counters ru0={{kMissing,kMissing,kMissing,kMissing}},ru1={{kMissing,kMissing,kMissing,kMissing}};
};
template<class R,class Work> const char* Capture(R& reader,Work work,Observation& o) {
    if (!reader.Mono(o.m0)) return "m0 capture";
    if (!reader.Usage(o.ru0)) return "ru0 capture";
    if (!reader.Cpu(o.c0)) return "c0 capture";
    if (!reader.Mono(o.m1)) return "m1 capture";
    std::atomic_signal_fence(std::memory_order_seq_cst);
    work();
    std::atomic_signal_fence(std::memory_order_seq_cst);
    if (!reader.Mono(o.m2)) return "m2 capture";
    if (!reader.Cpu(o.c1)) return "c1 capture";
    if (!reader.Usage(o.ru1)) return "ru1 capture";
    if (!reader.Mono(o.m3)) return "m3 capture";
    return nullptr;
}
void Validate(const Observation& o,const Observation* previous) {
    for (uint64_t v : {o.m0,o.c0,o.m1,o.m2,o.c1,o.m3}) Require(v != kMissing,"incomplete observation");
    Require(o.m0 <= o.m1 && o.m1 < o.m2 && o.m2 <= o.m3 && o.c0 <= o.c1,"clock ordering");
    for (unsigned k=0;k<4;++k) Require(o.ru0[k] != kMissing && o.ru1[k] != kMissing && o.ru0[k] <= o.ru1[k],"counter ordering");
    if (previous) {
        Require(previous->m3 <= o.m0 && previous->c1 <= o.c0,"cross-record clock ordering");
        for (unsigned k=0;k<4;++k) Require(previous->ru1[k] <= o.ru0[k],"cross-record counter ordering");
    }
}
void Budget(uint64_t start,uint64_t now,uint64_t encode) {
    Require(now >= start && now-start < kWallLimit,"worker wall deadline");
    Require(encode <= kEncodeLimit,"cumulative encode cap");
}
WH2_ALIGNMENT_NOINLINE uint64_t Prelude(uint64_t x,unsigned count) {
    for (unsigned i=0;i<count;++i) { x^=x<<13; x^=x>>7; x^=x<<17; }
    return x;
}
struct Allocation {
    uint8_t* data=nullptr; size_t bytes=0;
    Allocation()=default;
    Allocation(const Allocation&)=delete;
    Allocation& operator=(const Allocation&)=delete;
    ~Allocation() { std::free(data); }
    void Create(size_t length) {
        Require(!data,"allocation reused"); void* p=nullptr;
        Require(posix_memalign(&p,4096,length) == 0 && p,"aligned allocation");
        data=static_cast<uint8_t*>(p); bytes=length; std::memset(data,0xa5,bytes);
    }
};
struct Fixture {
    std::array<uint8_t,kMessage+128> source;
    std::array<uint8_t,kIntermediate> master;
    Packets packets;
    std::array<std::vector<uint32_t>,6> columns;
    std::array<uint64_t,6> operations={{0}};
    std::array<B::Snapshot,2> snapshots;
    std::array<WirehairV2Codec,2> handles={{nullptr,nullptr}};
    std::array<Allocation,2> carriers,outputs;
    Allocation scratch;
    std::array<std::array<uint8_t,32>,2> profiles;
    std::array<uint64_t,4> preflight={{0,0,0,0}};
    volatile uint8_t touch_sink=0;
    Fixture() { source.fill(0xa5); master.fill(0); packets.fill(0); }
    ~Fixture() { for (auto h : handles) if (h) wirehair_v2_free(h); }
    const uint8_t* Source() const { return source.data()+64; }
    uint8_t* Output(unsigned output,unsigned j) { return outputs[output].data+j*kStride+64; }
    uint8_t* View(unsigned cell) { return carriers[cell/2].data+4096+(cell%2)*16; }
    WH2_ALIGNMENT_NOINLINE void Prepare(unsigned phase,unsigned cell) {
        std::memset(carriers[0].data,0xa5,kCarrier); std::memset(carriers[1].data,0xa5,kCarrier);
        std::memcpy(View(phase ? cell : 0),master.data(),kIntermediate);
        for (unsigned c=0;c<2;++c) {
            volatile const uint8_t* p=carriers[c].data;
            for (unsigned i=0;i<kCarrier;i+=64) touch_sink^=p[i];
        }
        std::memset(scratch.data,0,kScratch);
        std::memset(outputs[0].data,0xa5,kOutput); std::memset(outputs[1].data,0xa5,kOutput);
    }
    bool Guards(unsigned phase,unsigned cell) const {
        const unsigned view_cell=phase ? cell : 0, selected_carrier=view_cell/2;
        const unsigned begin=4096+(view_cell%2)*16, selected_output=phase ? 0 : cell%2;
        for (unsigned c=0;c<2;++c) for (unsigned i=0;i<kCarrier;++i) {
            const uint8_t expected=c == selected_carrier && i >= begin && i < begin+kIntermediate ? master[i-begin] : 0xa5;
            if (carriers[c].data[i] != expected) return false;
        }
        for (unsigned o=0;o<2;++o) for (unsigned i=0;i<kOutput;++i) {
            const unsigned j=i/kStride, within=i%kStride;
            if (o == selected_output && j < 6 && within >= 64 && within < 64+kBlock) continue;
            if (outputs[o].data[i] != 0xa5) return false;
        }
        for (unsigned i=0;i<source.size();++i) {
            const uint8_t expected=i >= 64 && i < 64+kMessage ? static_cast<uint8_t>((37*(i-64)+(i-64)/11)&255) : 0xa5;
            if (source[i] != expected) return false;
        }
        // This fixed, untimed read also makes the otherwise unused scratch
        // reset observable. The noipa preparation boundary prevents elision.
        const volatile uint8_t* scratch_bytes=scratch.data;
        for (unsigned i=0;i<kScratch;i+=64) if (scratch_bytes[i] != 0) return false;
        if (scratch_bytes[kScratch-1] != 0) return false;
        return true;
    }
    bool Matches(unsigned phase,unsigned cell) {
        bool valid=Guards(phase,cell); const unsigned out=phase ? 0 : cell%2;
        for (unsigned j=0;j<6;++j) valid &= std::memcmp(Output(out,j),packets.data()+j*kBlock,kBlock) == 0;
        return valid;
    }
};
struct Arm { const B::Snapshot* snapshot; const uint8_t* intermediate; };
using PacketFunction=B::Evaluation (*)(const Arm&,uint32_t,uint8_t*);
B::Evaluation PublicPacket(const Arm& arm,uint32_t id,uint8_t* output) {
    uint32_t bytes=UINT32_MAX;
    const WirehairV2Result code=wirehair_v2_encode(arm.snapshot->handle,id,output,kBlock,&bytes);
    B::Evaluation result; result.status=code; result.bytes=bytes; result.operations=kMissing;
    return result;
}
B::Evaluation ShadowPacket(const Arm& arm,uint32_t id,uint8_t* output) {
    return B::EvaluateShadow(*arm.snapshot,arm.intermediate,id,output);
}
struct CallResult {
    uint64_t calls=0;
    bool completed=false;
    std::array<uint64_t,6> codes={{kMissing,kMissing,kMissing,kMissing,kMissing,kMissing}};
    std::array<uint64_t,6> lengths={{kMissing,kMissing,kMissing,kMissing,kMissing,kMissing}};
    std::array<uint64_t,6> ops={{kMissing,kMissing,kMissing,kMissing,kMissing,kMissing}};
};
WH2_ALIGNMENT_NOINLINE void RunPackets(PacketFunction function,const Arm& arm,uint8_t* output,bool shadow,
                                     const std::array<uint64_t,6>& expected_ops,CallResult& result) {
    for (unsigned cycle=0;cycle<64;++cycle) for (unsigned j=0;j<6;++j) {
        const B::Evaluation value=function(arm,6+j,output+j*kStride+64); ++result.calls;
        result.codes[j]=uint64_t(value.status); result.lengths[j]=value.bytes;
        if (shadow) {
            if (!cycle) result.ops[j]=0;
            if (value.operations != expected_ops[j] || value.operations >= kMissing-result.ops[j]) {
                // INVALID prefixes retain the offending returned value here;
                // only COMPLETE records promise per-slot operation sums.
                result.ops[j]=value.operations; return;
            }
            result.ops[j]+=value.operations;
        }
        if (value.status != WirehairV2_Success || value.bytes != kBlock) return;
    }
    result.completed=true;
}
struct Record : Observation {
    unsigned index=0,rep=0,order=0,phase=0,comparison=0,position=0,cell=0;
    CallResult result;
};
struct State {
    std::array<Record,3360> records;
    size_t started=0;
    uint64_t callbacks=0,encode_calls=0,checked_packets=0,checksum=kHashSeed,sum_encode=0;
    Observation prelude;
    uint64_t prelude_final=kMissing;
    State() {
        unsigned index=0;
        for (unsigned rep=0;rep<12;++rep) for (unsigned step=0;step<2;++step)
            for (unsigned phase_step=0;phase_step<2;++phase_step) {
                const unsigned phase=(rep+step+phase_step)%2,classes=phase ? 6 : 8,order=(rep+step)%2;
                for (unsigned class_step=0;class_step<classes;++class_step) {
                    const unsigned comparison=(2*rep+step+class_step)%classes;
                    for (unsigned position=0;position<10;++position) {
                        Record& r=records[index]; r.index=index++; r.rep=rep; r.order=order;
                        r.phase=phase; r.comparison=comparison; r.position=position;
                        const unsigned side=(position < 2 ? position : kSlots[position-2])^order;
                        r.cell=phase ? kShadow[comparison][side] : kNatural[comparison][side];
                    }
                }
            }
        Require(index == records.size(),"frozen roster count");
    }
};
uint64_t PacketDigest(Fixture& fixture,unsigned output) {
    uint64_t digest=kHashSeed;
    for (unsigned j=0;j<6;++j) for (unsigned i=0;i<kBlock;++i) digest=(digest*kHashFactor)^fixture.Output(output,j)[i];
    return digest;
}
uint64_t FoldDigest(uint64_t value,uint64_t digest) {
    for (unsigned i=0;i<8;++i) value=(value*kHashFactor)^((digest>>(8*i))&255);
    return value;
}
void Account(State& state,const Record& r,Fixture& fixture) {
    if (r.m1 != kMissing && r.m2 != kMissing && r.m2 >= r.m1) {
        Require(r.m2-r.m1 <= UINT64_MAX-state.sum_encode,"encode sum overflow"); state.sum_encode+=r.m2-r.m1;
    }
    if (!r.result.calls) return;
    ++state.callbacks; state.encode_calls+=r.result.calls;
    const bool guards=fixture.Matches(r.phase,r.cell);
    bool results=r.result.completed && r.result.calls == 384;
    for (unsigned j=0;j<6;++j) results &= r.result.codes[j] == 0 && r.result.lengths[j] == kBlock &&
        r.result.ops[j] == (r.phase ? 64*fixture.operations[j] : kMissing);
    Require(results,"packet status/length/operations mismatch");
    Require(guards,"packet/source/carrier/output guard or oracle mismatch");
    state.checked_packets+=6; state.checksum=FoldDigest(state.checksum,PacketDigest(fixture,r.phase ? 0 : r.cell%2));
}
std::string HandleJson(const B::Snapshot& s) {
    const auto& p=s.system.Params; std::ostringstream out;
    out << "{\"profile_hex\":" << Quote(Hex(s.serialized_profile,32))
        << ",\"profile_sha256\":" << Quote(Sha256Hex(s.serialized_profile,32))
        << ",\"source_sha256\":" << Quote(Sha256Hex(s.source,kMessage))
        << ",\"intermediate_sha256\":" << Quote(Sha256Hex(s.intermediate,kIntermediate))
        << ",\"message_bytes\":" << s.message_bytes << ",\"block_bytes\":" << s.block_bytes
        << ",\"source_count\":" << s.source_count << ",\"precode_count\":" << s.precode_count
        << ",\"intermediate_bytes\":" << s.intermediate_bytes << ",\"source_policy\":" << unsigned(s.source_policy)
        << ",\"seed_attempt\":" << unsigned(s.profile.seed_attempt)
        << ",\"params\":{\"block_count\":" << p.BlockCount << ",\"staircase\":" << p.Staircase
        << ",\"dense_rows\":" << p.DenseRows << ",\"heavy_rows\":" << p.HeavyRows << ",\"source_hits\":" << p.SourceHits
        << ",\"dense_identity_corner\":" << (p.DenseIdentityCorner?"true":"false")
        << ",\"heavy_family\":" << unsigned(p.HeavyFamily) << ",\"dense_anchors\":" << unsigned(p.DenseAnchors)
        << ",\"seed\":" << p.Seed << "},\"config\":{\"peel_seed\":" << s.config.PeelSeed
        << ",\"mix_count\":" << s.config.MixCount << "},\"runtime\":{\"source_prime\":" << s.runtime.SourcePrime()
        << ",\"precode_prime\":" << s.runtime.PrecodePrime() << "}}"; return out.str();
}
std::string Metadata(Fixture& f) {
    std::ostringstream out;
    out << "{\"source_hex\":" << Quote(Hex(f.Source(),kMessage)) << ",\"source_sha256\":" << Quote(Sha256Hex(f.Source(),kMessage))
        << ",\"intermediate_hex\":" << Quote(Hex(f.master.data(),kIntermediate))
        << ",\"intermediate_sha256\":" << Quote(Sha256Hex(f.master.data(),kIntermediate))
        << ",\"expected_packets_hex\":" << Quote(Hex(f.packets.data(),kMessage))
        << ",\"expected_packets_sha256\":" << Quote(Sha256Hex(f.packets.data(),kMessage)) << ",\"columns\":[";
    for (unsigned j=0;j<6;++j) {
        if (j) out << ',';
        out << '[';
        for (size_t k=0;k<f.columns[j].size();++k) { if (k) out << ','; out << f.columns[j][k]; } out << ']';
    }
    out << "],\"expected_operations\":"; Numbers(out,f.operations);
    out << ",\"handles\":[" << HandleJson(f.snapshots[0]) << ',' << HandleJson(f.snapshots[1]) << "]}"; return out.str();
}
std::string Addresses(Fixture& f) {
    std::ostringstream out; out << "{\"source\":" << PointerJson(f.Source()) << ",\"master\":" << PointerJson(f.master.data());
    out << ",\"handles\":[" << PointerJson(f.handles[0]) << ',' << PointerJson(f.handles[1]) << ']';
    out << ",\"intermediates\":[" << PointerJson(f.snapshots[0].intermediate) << ',' << PointerJson(f.snapshots[1].intermediate) << ']';
    out << ",\"carriers\":[" << PointerJson(f.carriers[0].data) << ',' << PointerJson(f.carriers[1].data) << ']';
    out << ",\"outputs\":[" << PointerJson(f.outputs[0].data) << ',' << PointerJson(f.outputs[1].data) << ']';
    out << ",\"scratch\":" << PointerJson(f.scratch.data)
        << ",\"public_function\":" << Address(reinterpret_cast<const void*>(&wirehair_v2_encode))
        << ",\"shadow_function\":" << Address(reinterpret_cast<const void*>(&B::EvaluateShadow))
        << ",\"runner_function\":" << Address(reinterpret_cast<const void*>(&RunPackets))
        << ",\"prelude_function\":" << Address(reinterpret_cast<const void*>(&Prelude)) << '}'; return out.str();
}
void Initialize(Fixture& f) {
    for (unsigned i=0;i<kMessage;++i) f.source[i+64]=static_cast<uint8_t>((37*i+i/11)&255);
    WirehairV2EncoderOptions options=WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy=WirehairV2EncoderSource_BorrowedImmutable;
    for (unsigned h=0;h<2;++h) {
        uint32_t bytes=0;
        Require(wirehair_v2_encoder_create_profile_id_with_options(WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            f.Source(),kMessage,kBlock,&options,f.profiles[h].data(),32,&bytes,&f.handles[h]) == WirehairV2_Success &&
            bytes == 32 && f.handles[h],"borrowed public creation");
        Require(B::Capture(f.handles[h],f.snapshots[h]) == WirehairV2_Success,"private snapshot capture");
        const B::Snapshot& s=f.snapshots[h];
        Require(s.handle == f.handles[h] && s.source == f.Source() && s.intermediate &&
            s.message_bytes == kMessage && s.block_bytes == kBlock && s.source_count == 6 &&
            s.precode_count == 30 && s.intermediate_bytes == kIntermediate &&
            s.source_policy == WirehairV2EncoderSource_BorrowedImmutable &&
            std::memcmp(s.serialized_profile,f.profiles[h].data(),32) == 0,"frozen snapshot dimensions/profile");
    }
    Require(f.handles[0] != f.handles[1] && f.snapshots[0].intermediate != f.snapshots[1].intermediate &&
        f.profiles[0] == f.profiles[1] && HandleJson(f.snapshots[0]) == HandleJson(f.snapshots[1]),"natural handles differ semantically");
    std::memcpy(f.master.data(),f.snapshots[0].intermediate,kIntermediate);
    Require(std::memcmp(f.master.data(),f.snapshots[1].intermediate,kIntermediate) == 0,"intermediate equality");
    for (unsigned i=0;i<2;++i) { f.carriers[i].Create(kCarrier); f.outputs[i].Create(kOutput); }
    f.scratch.Create(kScratch);
    for (unsigned j=0;j<6;++j) {
        const auto& s=f.snapshots[0];
        f.columns[j]=wirehair_v2::GeneratePacketMatrixRowWithRuntime(6,30,6+j,s.config,s.runtime);
        Require(!f.columns[j].empty() && f.columns[j].size() <= 36,"packet column roster");
        std::array<bool,36> used={{false}};
        for (uint32_t column : f.columns[j]) {
            Require(column < 36 && !used[column],"duplicate or invalid packet column"); used[column]=true;
            for (unsigned byte=0;byte<kBlock;++byte) f.packets[j*kBlock+byte]^=f.master[column*kBlock+byte];
        }
        // The production counter includes the initial copy/set, not only XORs.
        ++f.preflight[3]; f.operations[j]=f.columns[j].size();
    }
    for (unsigned phase=0;phase<2;++phase) for (unsigned cell=0;cell<(phase ? 4u : 2u);++cell) {
        const unsigned natural_cell=cell*2;
        f.Prepare(phase,phase ? cell : natural_cell);
        const Arm arm={&f.snapshots[phase ? 0 : cell],phase ? f.View(cell) : f.snapshots[cell].intermediate};
        for (unsigned j=0;j<6;++j) {
            if (phase) Require(B::ValidateShadowBuffers(f.snapshots[0],f.View(cell),kIntermediate,f.Output(0,j),kBlock),"preflight shadow extents");
            const B::Evaluation value=(phase ? ShadowPacket : PublicPacket)(arm,6+j,f.Output(0,j)); ++f.preflight[phase];
            Require(value.status == WirehairV2_Success && value.bytes == kBlock &&
                (!phase || value.operations == f.operations[j]),"packet preflight status/operations");
        }
        Require(f.Matches(phase,phase ? cell : natural_cell),"packet preflight guards/oracle");
    }
    f.Prepare(0,0);
    for (unsigned j=0;j<6;++j) {
        Require(B::ValidateShadowBuffers(f.snapshots[0],f.snapshots[0].intermediate,kIntermediate,f.Output(0,j),kBlock),"original-view extents");
        const B::Evaluation v=B::EvaluateShadow(f.snapshots[0],f.snapshots[0].intermediate,6+j,f.Output(0,j)); ++f.preflight[2];
        Require(v.status == WirehairV2_Success && v.bytes == kBlock && v.operations == f.operations[j],"original-view evaluator");
    }
    Require(f.Matches(0,0) && f.preflight == std::array<uint64_t,4>{{12,24,6,6}},"complete preflight");
    for (unsigned cell=0;cell<4;++cell) for (unsigned j=0;j<6;++j)
        Require(B::ValidateShadowBuffers(f.snapshots[0],f.View(cell),kIntermediate,f.Output(0,j),kBlock),"shadow extents");
}
void ObservationJson(std::ostream& out,const Observation& o) {
    for (uint64_t value : {o.m0,o.c0,o.m1,o.m2,o.c1,o.m3}) { out << ','; Number(out,value); }
    out << ','; Numbers(out,o.ru0); out << ','; Numbers(out,o.ru1);
}
std::string RecordsJson(const State& s) {
    std::ostringstream out; out << '[';
    for (size_t i=0;i<s.started;++i) {
        if (i) out << ',';
        const Record& r=s.records[i];
        out << '[' << r.index << ',' << r.rep << ',' << r.order << ',' << r.phase << ',' << r.comparison << ',' << r.position << ',' << r.cell;
        ObservationJson(out,r); out << ',' << r.result.calls << ','; Numbers(out,r.result.codes);
        out << ','; Numbers(out,r.result.lengths); out << ','; Numbers(out,r.result.ops); out << ']';
    }
    out << ']'; return out.str();
}
std::string PreludeJson(const State& s) {
    const Observation& o=s.prelude; std::ostringstream out;
    out << "{\"iterations\":" << kPreludeCount << ",\"seed\":" << kPreludeSeed << ",\"final_state\":"; Number(out,s.prelude_final);
    const char* names[]={"m0","c0","m1","m2","c1","m3"}; unsigned i=0;
    for (uint64_t value : {o.m0,o.c0,o.m1,o.m2,o.c1,o.m3}) { out << ',' << Quote(names[i++]) << ':'; Number(out,value); }
    out << ",\"ru0\":"; Numbers(out,o.ru0); out << ",\"ru1\":"; Numbers(out,o.ru1); out << '}'; return out.str();
}
void FinalClock(bool captured,uint64_t now,uint64_t start,uint64_t latest,uint64_t encode,
                uint64_t& elapsed,std::string& outcome,std::string& failure) {
    const char* error=nullptr;
    if (!captured || now < start || now < latest || (outcome == "COMPLETE" && now-start < elapsed)) error="final clock regression/failure";
    else { elapsed=now-start; if (elapsed >= kWallLimit || encode > kEncodeLimit) error="post-observation deadline"; }
    if (error) { outcome="INVALID"; if (failure.empty()) failure=error; }
}
int Worker() {
    Reader reader; State state; Fixture fixture;
    uint64_t start=kMissing,elapsed=0,mono_resolution=0,cpu_resolution=0;
    std::string outcome="INVALID",failure,identity_before="null",identity_after="null",runtime_before="null",runtime_after="null";
    std::string metadata="null",addresses="null",handles_after="null",addresses_after="null";
    try {
        start=reader.Now();
        const rlimit memory={512u*1024u*1024u,512u*1024u*1024u},cpu={5,5};
        Require(setrlimit(RLIMIT_AS,&memory) == 0 && setrlimit(RLIMIT_CPU,&cpu) == 0,"resource limits"); alarm(10);
        Pin(); identity_before=Identity(); Require(gf256_init() == 0,"shared GF initialization"); runtime_before=Runtime();
        mono_resolution=reader.Resolution(CLOCK_MONOTONIC); cpu_resolution=reader.Resolution(CLOCK_THREAD_CPUTIME_ID);
        Initialize(fixture); metadata=Metadata(fixture); addresses=Addresses(fixture);
        reader.Prepare(); Touch(&state,sizeof(state));
        const char* prelude_error=Capture(reader,[&] { state.prelude_final=Prelude(kPreludeSeed,kPreludeCount); },state.prelude);
        Require(prelude_error == nullptr,prelude_error ? prelude_error : "prelude capture"); Validate(state.prelude,nullptr);
        Require(state.prelude_final == kPreludeFinal,"prelude final state");
        for (size_t i=0;i<state.records.size();++i) {
            Budget(start,reader.Now(),state.sum_encode); Require(OnCpu(),"pre-callback CPU drift");
            Record& r=state.records[i]; fixture.Prepare(r.phase,r.cell);
            const Arm arm={&fixture.snapshots[r.phase ? 0 : r.cell/2],r.phase ? fixture.View(r.cell) : fixture.snapshots[r.cell/2].intermediate};
            uint8_t* output=fixture.outputs[r.phase ? 0 : r.cell%2].data;
            state.started=i+1;
            const char* error=Capture(reader,[&] { RunPackets(r.phase ? ShadowPacket : PublicPacket,arm,output,r.phase != 0,fixture.operations,r.result); },r);
            Account(state,r,fixture);
            Require(error == nullptr,error ? error : "callback capture");
            Validate(r,i ? &state.records[i-1] : &state.prelude);
            Require(OnCpu(),"post-callback CPU drift"); Budget(start,reader.Now(),state.sum_encode);
        }
        Require(state.callbacks == 3360 && state.encode_calls == 1290240 && state.checked_packets == 20160,"complete callback accounting");
        std::array<std::string,2> post_handles;
        for (unsigned h=0;h<2;++h) {
            B::Snapshot current; Require(B::Capture(fixture.handles[h],current) == WirehairV2_Success,"post snapshot capture");
            post_handles[h]=HandleJson(current);
            Require(current.intermediate == fixture.snapshots[h].intermediate && current.source == fixture.Source() &&
                post_handles[h] == HandleJson(fixture.snapshots[h]) && std::memcmp(current.intermediate,fixture.master.data(),kIntermediate) == 0,
                "post private storage/profile drift");
        }
        handles_after="["+post_handles[0]+","+post_handles[1]+"]";
        addresses_after=Addresses(fixture);
        Require(Metadata(fixture) == metadata && Addresses(fixture) == addresses,"post fixture metadata drift");
        identity_after=Identity(); runtime_after=Runtime();
        Require(identity_before == identity_after && runtime_before == runtime_after,"post identity/runtime drift");
        const uint64_t now=reader.Now(); Budget(start,now,state.sum_encode); elapsed=now-start; outcome="COMPLETE";
    } catch (const std::exception& error) { failure=error.what(); }
      catch (...) { failure="unknown alignment failure"; }
    if (start != kMissing) {
        uint64_t latest=start,now=kMissing;
        for (uint64_t value : {state.prelude.m0,state.prelude.m1,state.prelude.m2,state.prelude.m3})
            if (value != kMissing && value > latest) latest=value;
        for (size_t i=0;i<state.started;++i) for (uint64_t value : {state.records[i].m0,state.records[i].m1,state.records[i].m2,state.records[i].m3})
            if (value != kMissing && value > latest) latest=value;
        const bool captured=reader.Mono(now); FinalClock(captured,now,start,latest,state.sum_encode,elapsed,outcome,failure);
    }
    std::ostringstream out;
    out << "{\"protocol\":" << Quote(kProtocol) << ",\"schema\":" << Quote(std::string(kProtocol)+".raw.v1")
        << ",\"outcome\":" << Quote(outcome) << ",\"failure\":" << (failure.empty()?"null":Quote(failure)) << ",\"target_cpu\":50"
        << ",\"identity_before\":" << identity_before << ",\"identity_after\":" << identity_after
        << ",\"runtime_before\":" << runtime_before << ",\"runtime_after\":" << runtime_after
        << ",\"metadata\":" << metadata << ",\"addresses\":" << addresses
        << ",\"handles_after\":" << handles_after << ",\"addresses_after\":" << addresses_after
        << ",\"preflight\":{\"public\":" << fixture.preflight[0]
        << ",\"shadow\":" << fixture.preflight[1] << ",\"original_view\":" << fixture.preflight[2] << ",\"scalar\":" << fixture.preflight[3] << '}'
        << ",\"prelude\":" << PreludeJson(state) << ",\"callbacks\":" << state.callbacks << ",\"encode_calls\":" << state.encode_calls
        << ",\"checked_packets\":" << state.checked_packets << ",\"checksum\":" << state.checksum << ",\"sum_encode_wall_ns\":" << state.sum_encode
        << ",\"worker_start_ns\":"; Number(out,start);
    out << ",\"monotonic_resolution_ns\":" << mono_resolution << ",\"thread_resolution_ns\":" << cpu_resolution
        << ",\"elapsed_ns\":" << elapsed << ",\"records\":" << RecordsJson(state) << "}\n";
    const std::string raw=out.str(); Require(raw.size() <= 2u*1024u*1024u,"raw output cap");
    std::cout << raw; std::cout.flush(); return std::cout.good() && outcome == "COMPLETE" ? 0 : 1;
}
enum class Selection { Invalid, Worker, Neutral };
Selection Select(int argc,const char* const* argv) {
    if (argc != 2 || !argv || !argv[1]) return Selection::Invalid;
    if (std::strcmp(argv[1],"--worker") == 0) return Selection::Worker;
    if (std::strcmp(argv[1],"--selftest") == 0) return Selection::Neutral;
    return Selection::Invalid;
}
template<class F> void Rejected(F function) {
    bool rejected=false; try { function(); } catch (const std::runtime_error&) { rejected=true; }
    Require(rejected,"expected neutral rejection");
}
struct FakeReader {
    unsigned reads=0,fail=0,event_count=0; uint64_t tick=0;
    std::array<char,9> events={{0}};
    bool Step(char event) { Require(event_count < 9,"mock event overflow"); events[event_count++]=event; ++tick; return ++reads != fail; }
    bool Mono(uint64_t& value) { if (!Step('M')) return false; value=tick*100; return true; }
    bool Cpu(uint64_t& value) { if (!Step('C')) return false; value=tick*100; return true; }
    bool Usage(Counters& value) { if (!Step('U')) return false; value={{0,0,0,0}}; return true; }
    void Work() { events[event_count++]='W'; tick+=10; }
};
unsigned fake_packet_calls=0,fake_bad_at=0,fake_bad_kind=0;
B::Evaluation FakePacket(const Arm&,uint32_t id,uint8_t* output) {
    ++fake_packet_calls;
    B::Evaluation value; value.status=WirehairV2_Success; value.bytes=kBlock; value.operations=id;
    std::memset(output,static_cast<int>(id),kBlock);
    if (fake_packet_calls == fake_bad_at) {
        if (fake_bad_kind == 0) value.status=WirehairV2_Error;
        else if (fake_bad_kind == 1) value.bytes=kBlock-1;
        else if (fake_bad_kind == 2) value.operations+=1;
        else value.operations*=64;
    }
    return value;
}
int Neutral() {
    Require(B::NeutralSelfTest(),"neutral bridge tests");
    const char* worker[]={"test","--worker"},*selftest[]={"test","--selftest"},*wrong[]={"test","--old-worker"};
    Require(Select(2,worker) == Selection::Worker && Select(2,selftest) == Selection::Neutral &&
        Select(2,wrong) == Selection::Invalid && Select(1,worker) == Selection::Invalid && Select(2,nullptr) == Selection::Invalid,"main dispatch");
    Require(Prelude(kPreludeSeed,0) == kPreludeSeed && Prelude(kPreludeSeed,1) == UINT64_C(0xdc1b77ae0bf34dad) &&
        Prelude(kPreludeSeed,2) == UINT64_C(0x64f0eeb9026e6076) && Prelude(kPreludeSeed,3) == UINT64_C(0x7b07ce91e5906136) &&
        Prelude(kPreludeSeed,17) == UINT64_C(0x353cfc387dfae6b8),"small neutral recurrence");
    State state; unsigned counts[2][8][2]={};
    static const unsigned expected_sides[10]={0,1,0,1,1,0,1,0,0,1};
    for (size_t i=0;i<state.records.size();++i) {
        const Record& r=state.records[i];
        const unsigned epoch=static_cast<unsigned>(i/140),rep=epoch/2,step=epoch%2,order=(rep+step)%2;
        const unsigned local=static_cast<unsigned>(i%140),first_count=order ? 60 : 80;
        const unsigned phase=local < first_count ? order : 1-order,relative=local < first_count ? local : local-first_count;
        const unsigned comparison=(epoch+relative/10)%(phase ? 6 : 8),side=expected_sides[i%10]^order;
        unsigned cell=comparison;
        if (comparison >= 4) {
            if (!phase && comparison < 6) cell=comparison-4+2*side;
            else cell=2*(comparison-(phase ? 4 : 6))+side;
        }
        Require(r.index == i && r.rep == rep && r.order == order && r.phase == phase &&
            r.comparison == comparison && r.position == i%10 && r.cell == cell,"neutral full roster/rotation");
        if (!r.position) ++counts[r.phase][r.comparison][r.order];
    }
    for (unsigned phase=0;phase<2;++phase) for (unsigned c=0;c<(phase?6u:8u);++c)
        for (unsigned order=0;order<2;++order) Require(counts[phase][c][order] == 12,"neutral exact panel counts");
    // Unrelated byte fixture: no live codec, GF initialization or fixed data.
    Fixture neutral;
    for (unsigned i=0;i<kMessage;++i) neutral.source[i+64]=static_cast<uint8_t>((37*i+i/11)&255);
    for (unsigned i=0;i<kIntermediate;++i) neutral.master[i]=static_cast<uint8_t>(i*11+i/17);
    for (unsigned j=0;j<6;++j) {
        neutral.operations[j]=j+6;
        std::memset(neutral.packets.data()+j*kBlock,static_cast<int>(j+6),kBlock);
    }
    for (unsigned i=0;i<2;++i) { neutral.carriers[i].Create(kCarrier); neutral.outputs[i].Create(kOutput); }
    neutral.scratch.Create(kScratch);
    const Arm mock_arm={nullptr,nullptr};
    for (unsigned phase=0;phase<2;++phase) for (unsigned cell=0;cell<4;++cell) {
        neutral.Prepare(phase,cell); fake_packet_calls=0; fake_bad_at=0; CallResult result;
        RunPackets(FakePacket,mock_arm,neutral.outputs[phase ? 0 : cell%2].data,phase != 0,neutral.operations,result);
        Require(fake_packet_calls == 384 && result.completed && result.calls == 384 && neutral.Matches(phase,cell),"neutral common packet runner");
        for (unsigned j=0;j<6;++j) Require(result.codes[j] == 0 && result.lengths[j] == kBlock &&
            result.ops[j] == (phase ? 64*neutral.operations[j] : kMissing),"neutral runner result fields");
        Record record; record.phase=phase; record.cell=cell; record.result=result; record.m1=100; record.m2=110;
        State accounting; Account(accounting,record,neutral);
        Require(accounting.callbacks == 1 && accounting.encode_calls == 384 && accounting.checked_packets == 6 &&
            accounting.sum_encode == 10,"neutral valid accounting");
        const unsigned chosen_carrier=phase ? cell/2 : 0;
        neutral.carriers[chosen_carrier].data[0]^=1;
        Require(!neutral.Guards(phase,cell),"neutral entire carrier guard"); neutral.carriers[chosen_carrier].data[0]^=1;
        neutral.outputs[phase ? 0 : cell%2].data[64+kBlock]^=1;
        Require(!neutral.Guards(phase,cell),"neutral output tail guard"); neutral.outputs[phase ? 0 : cell%2].data[64+kBlock]^=1;
        neutral.scratch.data[0]=1; Require(!neutral.Guards(phase,cell),"neutral scratch observable reset"); neutral.scratch.data[0]=0;
    }
    for (unsigned kind=0;kind<4;++kind) {
        neutral.Prepare(1,0); fake_packet_calls=0; fake_bad_at=kind == 3 ? 384 : kind ? 101 : 1; fake_bad_kind=kind;
        Record record; record.phase=1; record.m1=100; record.m2=120;
        RunPackets(FakePacket,mock_arm,neutral.outputs[0].data,true,neutral.operations,record.result);
        Require(!record.result.completed && record.result.calls == fake_bad_at && fake_packet_calls == fake_bad_at,"neutral first bad call stops");
        if (kind == 3) Require(record.result.ops[5] == 64*neutral.operations[5] && neutral.Matches(1,0),"last-call sum-shaped bad ops witness");
        State partial; partial.started=1; partial.records[0]=record;
        Rejected([&] { Account(partial,record,neutral); });
        Require(partial.callbacks == 1 && partial.encode_calls == fake_bad_at && partial.checked_packets == 0 &&
            partial.sum_encode == 20 && RecordsJson(partial).size() > 30,"neutral invalid prefix accounting");
    }
    fake_bad_at=0;
    for (unsigned failure=0;failure<=8;++failure) {
        FakeReader reader; reader.fail=failure; Observation o; unsigned work=0;
        const char* error=Capture(reader,[&] { ++work; reader.Work(); },o);
        Require((error != nullptr) == (failure != 0) && work == (failure == 0 || failure >= 5 ? 1u : 0u),"neutral capture failure prefix");
        if (!failure) {
            Validate(o,nullptr);
            Require(reader.events == std::array<char,9>{{'M','U','C','M','W','M','C','U','M'}} &&
                o.m2-o.m1 == 1100 && o.c1-o.c0 == 1300,"neutral bracket order");
            Observation bad=o; bad.m2=bad.m1; Rejected([&] { Validate(bad,nullptr); });
            bad=o; bad.ru0[0]=1; Rejected([&] { Validate(bad,nullptr); });
            bad=o; bad.m0=o.m3-1; Rejected([&] { Validate(bad,&o); });
            bad=o; bad.c0=o.c1-1; Rejected([&] { Validate(bad,&o); });
        }
    }
    uint64_t ns=kMissing; timespec t={1,2}; Require(Nanoseconds(t,ns) && ns == 1000000002,"timespec");
    t.tv_nsec=1000000000L; Require(!Nanoseconds(t,ns),"timespec nanos"); t.tv_nsec=0; t.tv_sec=-1;
    Require(!Nanoseconds(t,ns),"negative timespec"); t.tv_sec=std::numeric_limits<time_t>::max(); Require(!Nanoseconds(t,ns),"timespec overflow");
    Budget(1,2,kEncodeLimit); Rejected([] { Budget(1,0,0); }); Rejected([] { Budget(0,kWallLimit,0); }); Rejected([] { Budget(0,1,kEncodeLimit+1); });
    for (unsigned failure=0;failure<4;++failure) {
        std::string outcome="COMPLETE",error; uint64_t elapsed=5;
        FinalClock(failure != 0,failure == 1 ? 0 : failure == 2 ? 4 : kWallLimit+1,1,3,0,elapsed,outcome,error);
        Require(outcome == "INVALID" && !error.empty(),"neutral final clock invalidation");
    }
    {
        std::string outcome="INVALID",error="first failure"; uint64_t elapsed=0;
        FinalClock(false,0,1,2,0,elapsed,outcome,error); Require(error == "first failure","neutral first failure retained");
        outcome="COMPLETE"; error.clear(); FinalClock(true,6,1,5,0,elapsed,outcome,error);
        Require(outcome == "COMPLETE" && error.empty() && elapsed == 5,"neutral valid final clock");
    }
    std::cout << "REPAIR ALIGNMENT R0 neutral selftest PASS\n"; std::cout.flush(); return std::cout.good() ? 0 : 1;
}
} // namespace alignment

int main(int argc,char** argv) {
    const alignment::Selection choice=alignment::Select(argc,argv);
    if (choice == alignment::Selection::Invalid) return 2;
    try { return choice == alignment::Selection::Neutral ? alignment::Neutral() : alignment::Worker(); }
    catch (const std::exception& error) { std::cerr << std::string(error.what()).substr(0,1000) << '\n'; return 1; }
    catch (...) { std::cerr << "alignment publication failure\n"; return 1; }
}
#undef WH2_ALIGNMENT_NOINLINE
