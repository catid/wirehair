// .73 frozen, one-shot four-arm matvec cost screen. No production source is modified.
#include "Wh2ThueMorseMatvecCostR0Bridge.h"
#include "Wh2FrozenTrace.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "gf256.h"
#include "Wh2ThueMorseMatvecGfniR0.h"
#include "wirehair/wirehair.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
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
#include <fcntl.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#if defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WH_ALIGN64) || defined(WH_HUGEPAGE) || defined(WH_COUNT)
#error "Scientific and neutral cost-driver builds use unmodified runtime policy"
#endif
#ifndef WH2_MATVEC_COST_SANITIZERS
#error "Builder must identify sanitizer versus native execution"
#endif
#ifndef WH2_MATVEC_COST_LOOKUP_PATH
#error "Builder must identify the pinned scientific lookup path"
#endif
#if defined(__GNUC__) && !defined(__clang__)
#define WH2_COST_NOINLINE __attribute__((noinline, noipa))
#else
#define WH2_COST_NOINLINE __attribute__((noinline))
#endif

namespace cost {
namespace P = wirehair_wh2_bench;
using wirehair::wh2_benchmark::Sha256Hex;
using Counters = std::array<uint64_t,4>;
static const char kProtocol[] = "wirehair.wh2.thue-matvec-cost-r0";
static const char kClaim[] = "/var/tmp/wh2-thue-matvec-cost-r0/CLAIM.json";
static const uint64_t kMissing = UINT64_MAX;
static const int64_t kAbsent = INT64_MIN;
static const uint64_t kWall = UINT64_C(125000000000);
static const uint64_t kCpu = UINT64_C(120000000000);
static const uint64_t kWork = UINT64_C(8000000000);
static const uint64_t kPreludeSeed = UINT64_C(0x9e3779b97f4a7c15);
static const uint64_t kPreludeFinal = UINT64_C(0x43935dad1647741b);
static const unsigned kPreludeIterations = 1u << 20;
static const unsigned kBatch = 128u;
static const unsigned kPairs[7][2] = {{0,0},{1,1},{2,2},{3,3},{0,1},{2,1},{3,1}};
static const char kLookupHash[] = "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d";
static const size_t kLookupBytes = 39936u;

void Require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
uint64_t Address(const void* pointer) {
    return static_cast<uint64_t>(reinterpret_cast<uintptr_t>(pointer));
}
std::string Quote(const std::string& value) {
    std::ostringstream out; out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c >= 32 && c < 127) out << c;
        else out << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << unsigned(c) << std::dec;
    }
    out << '"'; return out.str();
}
std::string Hex(const void* bytes, size_t size) {
    const uint8_t* p = static_cast<const uint8_t*>(bytes);
    const char digits[] = "0123456789abcdef";
    std::string result(size * 2u, '0');
    for (size_t i=0; i<size; ++i) { result[2*i]=digits[p[i]>>4]; result[2*i+1]=digits[p[i]&15]; }
    return result;
}
std::string Hash(const uint8_t* data, size_t size) {
    return Sha256Hex(std::string(reinterpret_cast<const char*>(data), size));
}
void Number(std::ostream& out, uint64_t value) {
    if (value == kMissing) out << "null"; else out << value;
}
void Signed(std::ostream& out, int64_t value) {
    if (value == kAbsent) out << "null"; else out << value;
}
void CounterJson(std::ostream& out, const Counters& values) {
    out << '['; for (unsigned i=0; i<4; ++i) { if(i) out << ','; Number(out,values[i]); } out << ']';
}
bool Nanoseconds(const timespec& value, uint64_t& result) {
    if (value.tv_sec < 0 || value.tv_nsec < 0 || value.tv_nsec >= 1000000000L ||
        uint64_t(value.tv_sec) > (kMissing-1)/UINT64_C(1000000000)) return false;
    const uint64_t base=uint64_t(value.tv_sec)*UINT64_C(1000000000);
    if (uint64_t(value.tv_nsec) >= kMissing-base) return false;
    result=base+uint64_t(value.tv_nsec); return true;
}
bool OrderedClock(uint64_t value,uint64_t& previous) {
    if(previous!=kMissing && value<previous) return false;
    previous=value; return true;
}
struct Reader {
    timespec clock_scratch = {}; rusage usage_scratch = {};
    uint64_t last_mono=kMissing,last_cpu=kMissing;
    bool Clock(clockid_t id, uint64_t& value) {
        return clock_gettime(id,&clock_scratch)==0 && Nanoseconds(clock_scratch,value) &&
               OrderedClock(value,id==CLOCK_MONOTONIC?last_mono:last_cpu);
    }
    bool Mono(uint64_t& value) { return Clock(CLOCK_MONOTONIC,value); }
    bool Cpu(uint64_t& value) { return Clock(CLOCK_THREAD_CPUTIME_ID,value); }
    bool Usage(Counters& values) {
        if (getrusage(RUSAGE_THREAD,&usage_scratch) || usage_scratch.ru_minflt<0 ||
            usage_scratch.ru_majflt<0 || usage_scratch.ru_nvcsw<0 || usage_scratch.ru_nivcsw<0) return false;
        values={{uint64_t(usage_scratch.ru_minflt),uint64_t(usage_scratch.ru_majflt),
                 uint64_t(usage_scratch.ru_nvcsw),uint64_t(usage_scratch.ru_nivcsw)}}; return true;
    }
    uint64_t Now() { uint64_t value=kMissing; Require(Mono(value),"monotonic clock"); return value; }
    uint64_t Thread() { uint64_t value=kMissing; Require(Cpu(value),"thread clock"); return value; }
    uint64_t Resolution(clockid_t id) {
        uint64_t value=kMissing;
        Require(clock_getres(id,&clock_scratch)==0 && Nanoseconds(clock_scratch,value) && value>0,"clock resolution");
        return value;
    }
};
struct Observation {
    uint64_t m0=kMissing,c0=kMissing,m1=kMissing,m2=kMissing,c1=kMissing,m3=kMissing;
    Counters ru0={{kMissing,kMissing,kMissing,kMissing}}, ru1={{kMissing,kMissing,kMissing,kMissing}};
};
void ObservationJson(std::ostream& out, const Observation& o) {
    out << '[';
    const uint64_t clocks[]={o.m0,o.c0,o.m1,o.m2,o.c1,o.m3};
    for(unsigned i=0;i<6;++i) { if(i) out << ','; Number(out,clocks[i]); }
    out << ','; CounterJson(out,o.ru0); out << ','; CounterJson(out,o.ru1); out << ']';
}
template<class R,class Work> const char* Capture(R& reader, Work work, Observation& o) {
    if(!reader.Mono(o.m0)) return "m0 capture";
    if(!reader.Usage(o.ru0)) return "ru0 capture";
    if(!reader.Cpu(o.c0)) return "c0 capture";
    if(!reader.Mono(o.m1)) return "m1 capture";
    std::atomic_signal_fence(std::memory_order_seq_cst);
    work();
    std::atomic_signal_fence(std::memory_order_seq_cst);
    if(!reader.Mono(o.m2)) return "m2 capture";
    if(!reader.Cpu(o.c1)) return "c1 capture";
    if(!reader.Usage(o.ru1)) return "ru1 capture";
    if(!reader.Mono(o.m3)) return "m3 capture";
    return nullptr;
}
void Validate(const Observation& o, const Observation* previous) {
    for(uint64_t value : {o.m0,o.c0,o.m1,o.m2,o.c1,o.m3}) Require(value!=kMissing,"incomplete observation");
    Require(o.m0<=o.m1 && o.m1<o.m2 && o.m2<=o.m3 && o.c0<=o.c1,"clock ordering");
    Require(o.c1-o.c0<=o.m3-o.m0,"thread CPU exceeds enclosing wall interval");
    for(unsigned i=0;i<4;++i) Require(o.ru0[i]!=kMissing && o.ru1[i]!=kMissing && o.ru0[i]<=o.ru1[i],"counter ordering");
    if(previous) {
        Require(previous->m3<=o.m0 && previous->c1<=o.c0,"cross-record clocks");
        for(unsigned i=0;i<4;++i) Require(previous->ru1[i]<=o.ru0[i],"cross-record counters");
    }
}
bool OnCpu() {
    cpu_set_t mask; CPU_ZERO(&mask);
    return sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 &&
           CPU_ISSET(50,&mask) && sched_getcpu()==50;
}
void Pin() {
    cpu_set_t mask; CPU_ZERO(&mask); CPU_SET(50,&mask);
    Require(sched_setaffinity(0,sizeof(mask),&mask)==0 && OnCpu(),"singleton CPU50 pin");
}
std::string Identity() {
    P::TargetIdentityReceiptV2 identity; std::string diagnostic,canonical;
    Require(P::CapturePublicBorrowedTargetIdentity(50,identity,diagnostic),"target identity capture");
    Require(P::SerializeTargetIdentityV2(identity,canonical,diagnostic) && !canonical.empty() &&
            canonical.size()<=4096 && Sha256Hex(canonical)==identity.CanonicalSha256 &&
            identity.Before.Affinity==std::vector<int32_t>(1,50) &&
            identity.After.Affinity==std::vector<int32_t>(1,50),"target identity serialization");
    const auto& r=identity.Raw; const auto& d=identity.Derived;
    Require(d.Family==26 && d.Model==8 && d.Stepping==1 && d.FullApicId==100 &&
            d.CoreId==50 && d.PackageId==0 && d.ThreadId==0 && d.ThreadsPerCore==2 &&
            d.CcdId==6 && d.ComplexId==6 && d.LogicalProcessorsPerPackage==128,"frozen physical target");
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
    Require(GF256Ctx.Polynomial==0x14d && f.SSSE3==1 && f.AVX2==1 && f.GFNI==1 && f.AVX512==1,"frozen GF256 dispatch");
    std::ostringstream out;
    out << "{\"polynomial\":" << GF256Ctx.Polynomial << ",\"address\":" << Address(&GF256Ctx)
        << ",\"ssse3\":" << f.SSSE3 << ",\"avx2\":" << f.AVX2 << ",\"gfni\":" << f.GFNI << ",\"avx512\":" << f.AVX512 << '}';
    return out.str();
}
WH2_COST_NOINLINE uint64_t Prelude(uint64_t value, unsigned count) {
    for(unsigned i=0;i<count;++i) { value^=value<<13; value^=value>>7; value^=value<<17; }
    return value;
}

namespace T = wh2_matvec_cost_r0;
const T::Api* kApis[4] = {nullptr,nullptr,nullptr,nullptr};
void InitializeApis() {
    kApis[0]=wh2_matvec_cost_baseline_api();
    kApis[1]=wh2_matvec_cost_candidate_api();
    kApis[2]=wh2_matvec_cost_wh1_api();
    kApis[3]=wh2_matvec_cost_public_api();
    for(unsigned arm=0;arm<4;++arm) {
        const auto* api=kApis[arm];
        Require(api && api->create && api->encode && api->encoder_free &&
            api->decoder_create && api->feed && api->recover && api->decoder_free &&
            api->valid_profile && ((api->row!=nullptr)==(arm<2u)),"complete POD API");
    }
}
bool Same(const T::Result& a,const T::Result& b) {
    return a.status==b.status && a.code==b.code && a.required==b.required &&
        a.written==b.written && a.length_kind==b.length_kind;
}
void ResultJson(std::ostream& out,const T::Result& r) {
    out << "{\"status\":" << r.status << ",\"code\":" << r.code
        << ",\"required\":" << r.required << ",\"written\":" << r.written
        << ",\"length_kind\":" << r.length_kind << '}';
}
bool NoLength(const T::Result& r,int status) {
    return r.status==status && r.code==status && r.required==0 && r.written==0 &&
        r.length_kind==T::NoLength;
}
bool PacketResult(const T::Result& r,unsigned B) {
    return r.status==T::Success && r.code==0 && r.required==B && r.written==B &&
        r.length_kind==T::ReturnedLength;
}
bool RecoveryResult(const T::Result& r,unsigned arm,size_t M) {
    return r.status==T::Success && r.code==0 && r.required==M && r.written==M &&
        r.length_kind==(arm==2u?T::InferredLength:T::ReturnedLength);
}
bool FeedResult(const T::Result& r,unsigned arm,unsigned B) {
    return (r.status==T::Success || r.status==T::NeedMore) && r.code==r.status &&
        r.required==B && r.written==0 && r.length_kind==(arm<2u?T::ReturnedLength:T::InferredLength);
}
struct Counts {
    uint64_t encoder_creates=0,encoder_frees=0,encodes=0,decoder_creates=0,feeds=0,recovers=0,decoder_frees=0;
};
struct Handle {
    unsigned arm=0;
    T::State state={};
    uint64_t* free_count=nullptr;
    Handle()=default; Handle(const Handle&)=delete; Handle& operator=(const Handle&)=delete;
    void Reset() {
        if(state.handle) {
            if(state.kind==2u) kApis[arm]->decoder_free(&state);
            else kApis[arm]->encoder_free(&state);
            if(free_count) ++*free_count;
        }
    }
    ~Handle() { Reset(); }
};
struct Buffer {
    uint8_t* data=nullptr; size_t size=0;
    Buffer()=default; Buffer(const Buffer&)=delete; Buffer& operator=(const Buffer&)=delete;
    ~Buffer() { std::free(data); }
    void Create(size_t bytes) {
        Require(!data,"buffer already allocated"); void* pointer=nullptr;
        Require(posix_memalign(&pointer,4096,bytes)==0 && pointer,"external buffer allocation");
        data=static_cast<uint8_t*>(pointer); size=bytes; std::memset(data,0xa5,size);
    }
};


// Independent carryless multiplication and descending polynomial reduction.
// This is only used outside WORK, never as the timed packet evaluator.
uint8_t Multiply(uint8_t a,uint8_t b) {
    unsigned product=0;
    for(unsigned bit=0;bit<8;++bit) if(b&(1u<<bit)) product^=unsigned(a)<<bit;
    for(int bit=14;bit>=8;--bit) if(product&(1u<<bit)) product^=0x14du<<(bit-8);
    return uint8_t(product);
}
using Matrix = std::array<uint8_t,36>;
Matrix IdentityMatrix() {
    Matrix result={}; for(unsigned i=0;i<6;++i) result[6*i+i]=1; return result;
}
Matrix Product(const Matrix& a,const Matrix& b) {
    Matrix result={};
    for(unsigned r=0;r<6;++r) for(unsigned c=0;c<6;++c) for(unsigned k=0;k<6;++k)
        result[6*r+c]^=Multiply(a[6*r+k],b[6*k+c]);
    return result;
}
unsigned Parity(uint32_t value) {
    unsigned result=0; while(value) { result^=value&1u; value>>=1; } return result;
}
std::vector<uint8_t> NeutralLookup() {
    // The unrelated literal pair from the original neutral codec test.  No
    // selected pair, receipt, namespace or scientific file is read here.
    const uint8_t feedback[2][6]={{1,1,0,1,0,1},{3,1,0,1,0,1}};
    std::array<std::array<Matrix,32>,2> blocks={};
    for(unsigned phase=0;phase<2;++phase) {
        for(unsigned i=0;i<5;++i) blocks[phase][0][(i+1)*6+i]=1;
        for(unsigned i=0;i<6;++i) blocks[phase][0][i*6+5]=feedback[phase][i];
    }
    for(unsigned level=1;level<32;++level) {
        blocks[0][level]=Product(blocks[0][level-1],blocks[1][level-1]);
        blocks[1][level]=Product(blocks[1][level-1],blocks[0][level-1]);
    }
    std::vector<uint8_t> table; table.reserve(kLookupBytes);
    const unsigned specs[7][3]={{0,10,0},{0,10,1},{10,7,0},{10,7,1},{17,7,0},{17,7,1},{24,8,0}};
    for(const auto& spec : specs) {
        Matrix product=IdentityMatrix();
        for(unsigned value=0;value<(1u<<spec[1]);++value) {
            if(spec[0]==0) for(unsigned r=0;r<6;++r) table.push_back(product[r*6]);
            else table.insert(table.end(),product.begin(),product.end());
            product=Product(product,blocks[spec[2]^Parity(value)][spec[0]]);
        }
    }
    Require(table.size()==kLookupBytes,"neutral lookup extent"); return table;
}

// Independent 8/8/8/8 right-product oracle, not the codec's 10/7/7/8 lookup.
// Only Worker constructs the already-frozen pair; neutral execution supplies
// the unrelated original test pair. No selected row is evaluated at load time.
struct Independent {
    std::array<std::array<std::array<Matrix,256>,4>,2> chunks;
    explicit Independent(const uint8_t feedback[2][6]) {
        std::array<std::array<Matrix,32>,2> blocks={};
        for(unsigned phase=0;phase<2;++phase) {
            for(unsigned i=0;i<5;++i) blocks[phase][0][(i+1)*6+i]=1;
            for(unsigned i=0;i<6;++i) blocks[phase][0][i*6+5]=feedback[phase][i];
        }
        for(unsigned level=1;level<32;++level) {
            blocks[0][level]=Product(blocks[0][level-1],blocks[1][level-1]);
            blocks[1][level]=Product(blocks[1][level-1],blocks[0][level-1]);
        }
        for(unsigned phase=0;phase<2;++phase) for(unsigned chunk=0;chunk<4;++chunk) {
            Matrix p=IdentityMatrix();
            for(unsigned byte=0;byte<256;++byte) {
                chunks[phase][chunk][byte]=p;
                p=Product(p,blocks[phase^Parity(byte)][8u*chunk]);
            }
        }
    }
    std::array<uint8_t,6> Row(uint32_t id) const {
        Matrix p=IdentityMatrix();
        for(int chunk=3;chunk>=0;--chunk) {
            const unsigned shift=unsigned(chunk)*8u;
            const unsigned phase=chunk==3?0u:Parity(id>>(shift+8u));
            p=Product(p,chunks[phase][unsigned(chunk)][(id>>shift)&255u]);
        }
        std::array<uint8_t,6> row={};
        for(unsigned i=0;i<6;++i) row[i]=p[6u*i];
        return row;
    }
};
uint32_t PacketId(unsigned slot) { return slot<12u?slot:UINT32_MAX-2u*(slot-12u); }
unsigned DecodeSlot(unsigned family,unsigned step) {
    return step<6u?(family==0u?6u:12u)+step:step-6u;
}
struct DecodeEvidence {
    T::Result create={},recover={};
    std::array<T::Result,12> feeds={};
    unsigned count=0;
};
void DecodeJson(std::ostream& out,const DecodeEvidence& d) {
    out << "{\"create\":"; ResultJson(out,d.create); out << ",\"feeds\":[";
    for(unsigned i=0;i<d.count;++i) { if(i) out << ','; ResultJson(out,d.feeds[i]); }
    out << "],\"recover\":"; ResultJson(out,d.recover); out << '}';
}
struct Snapshot {
    unsigned rep=0,lane=0,arm=0;
    uint64_t address=0;
    std::string profile;
    T::Result create={};
    std::array<DecodeEvidence,2> decode={};
};
void SnapshotJson(std::ostream& out,const Snapshot& s) {
    out << "{\"rep\":" << s.rep << ",\"lane\":" << s.lane << ",\"arm\":" << s.arm
        << ",\"address\":" << s.address << ",\"profile_hex\":" << Quote(s.profile) << ",\"create\":";
    ResultJson(out,s.create); out << ",\"decode\":[";
    DecodeJson(out,s.decode[0]); out << ','; DecodeJson(out,s.decode[1]); out << "]}";
}
struct Fixture {
    unsigned width_index,B,reps; size_t M,stride;
    const std::vector<uint8_t>& lookup;
    Buffer source; std::array<Buffer,2> outputs;
    std::array<Handle,96> handles;
    std::array<std::vector<uint8_t>,4> packets,packet_reference;
    std::array<uint8_t,108> rows={};
    std::vector<Snapshot> snapshots;
    std::string source_hash;
    std::array<std::string,4> packet_hash;
    Fixture(unsigned index,unsigned block,unsigned repeat,const std::vector<uint8_t>& table,
            const Independent& oracle)
        :width_index(index),B(block),reps(repeat),M(size_t(6)*block),stride(block+128u),lookup(table) {
        Require((reps==12 && width_index<3 && B==(width_index==0?2u:width_index==1?64u:1280u)) ||
                (reps==1 && (B==2 || B==65)),"fixed fixture domain");
        Require(lookup.size()==kLookupBytes,"fixture lookup extent");
        source.Create(M+128u);
        for(size_t i=0;i<M;++i) source.data[64u+i]=uint8_t(37u*i+i/11u);
        source_hash=Hash(Source(),M);
        for(auto& output : outputs) output.Create(18u*stride);
        for(unsigned slot=0;slot<18;++slot) {
            const auto row=oracle.Row(PacketId(slot));
            std::copy(row.begin(),row.end(),rows.begin()+6u*slot);
        }
        for(auto& corpus:packets) corpus.assign(18u*B,0u);
        for(unsigned slot=0;slot<18;++slot) for(unsigned byte=0;byte<B;++byte)
            for(unsigned j=0;j<6;++j)
                packets[0][size_t(slot)*B+byte]^=Multiply(rows[6u*slot+j],Source()[size_t(j)*B+byte]);
        packets[1]=packets[0];
        snapshots.reserve(reps*8u);
    }
    const uint8_t* Source() const { return source.data+64u; }
    uint8_t* Output(unsigned lane,unsigned slot=0) { return outputs[lane].data+size_t(slot)*stride+64u; }
    Handle& At(unsigned rep,unsigned lane,unsigned arm) { return handles[(rep*2u+lane)*4u+arm]; }
    const Snapshot& Expected(unsigned rep,unsigned lane,unsigned arm) const {
        for(const Snapshot& s:snapshots) if(s.rep==rep && s.lane==lane && s.arm==arm) return s;
        throw std::runtime_error("missing preflight descriptor");
    }
    void Prepare() { for(auto& output : outputs) std::memset(output.data,0xa5,output.size); }
    void CheckImmutable() const {
        for(size_t i=0;i<source.size;++i) {
            const uint8_t expected=i<64u || i>=M+64u?0xa5u:uint8_t(37u*(i-64u)+(i-64u)/11u);
            Require(source.data[i]==expected,"immutable source/guard changed");
        }
        for(unsigned arm=0;arm<4;++arm)
            if(!packet_reference[arm].empty()) Require(packets[arm]==packet_reference[arm],"immutable packet corpus changed");
    }
    void CheckOutput(unsigned lane,unsigned arm,unsigned count,unsigned first_slot,bool recovered=false) const {
        const Buffer& output=outputs[lane];
        for(size_t i=0;i<output.size;++i) {
            uint8_t expected=0xa5u;
            if(recovered) {
                if(i>=64u && i<64u+M) expected=Source()[i-64u];
            } else {
                const size_t slot=i/stride,offset=i%stride;
                if(slot<count && offset>=64u && offset<64u+B)
                    expected=packets[arm][size_t(first_slot+slot)*B+offset-64u];
            }
            Require(output.data[i]==expected,"packet/recovery output or guard changed");
        }
    }
    void CheckHandles() const {
        for(const Snapshot& s:snapshots) {
            const Handle& h=handles[(s.rep*2u+s.lane)*4u+s.arm];
            Require(Address(h.state.handle)==s.address && kApis[s.arm]->valid_profile(&h.state) &&
                Hex(h.state.profile,h.state.profile_bytes)==s.profile,"persistent profile/handle changed");
        }
    }
    void Preflight(Counts& counts) {
        for(unsigned rep=0;rep<reps;++rep) for(unsigned step=0;step<4;++step) for(unsigned lane=0;lane<2;++lane) {
            const unsigned arm=(rep+width_index+step)%4u;
            const T::Api& api=*kApis[arm]; Handle& h=At(rep,lane,arm); h.arm=arm;
            h.free_count=&counts.encoder_frees;
            Snapshot s; s.rep=rep; s.lane=lane; s.arm=arm;
            s.create=api.create(lookup.data(),lookup.size(),Source(),M,B,&h.state);
            ++counts.encoder_creates;
            Require(NoLength(s.create,T::Success) && h.state.handle && api.valid_profile(&h.state),"preflight encoder create");
            s.address=Address(h.state.handle);
            Require(h.state.profile_bytes==(arm==3u?32u:0u),"public profile extent");
            s.profile=Hex(h.state.profile,h.state.profile_bytes);
            Prepare();
            for(unsigned slot=0;slot<18;++slot) {
                if(arm<2u) {
                    std::array<uint8_t,8> row={{0xa5,0,0,0,0,0,0,0xa5}};
                    const T::Result r=api.row(lookup.data(),lookup.size(),PacketId(slot),row.data()+1u);
                    Require(r.status==T::Success && r.code==0 && r.required==6u && r.written==6u &&
                        r.length_kind==T::InferredLength && row.front()==0xa5 && row.back()==0xa5 &&
                        std::equal(row.begin()+1,row.begin()+7,rows.begin()+6u*slot),"independent coefficient oracle");
                }
                const T::Result value=api.encode(&h.state,PacketId(slot),Output(lane,slot),B);
                ++counts.encodes; Require(PacketResult(value,B),"preflight encode result");
                // Public controls define their own packet corpus. Later live handles
                // must reproduce it; all arms separately recover original source.
                if(arm>=2u && rep==0u && lane==0u)
                    std::copy(Output(lane,slot),Output(lane,slot)+B,packets[arm].begin()+size_t(slot)*B);
                Require(std::memcmp(Output(lane,slot),packets[arm].data()+size_t(slot)*B,B)==0,"preflight packet oracle");
                if(slot<6u) Require(std::memcmp(Output(lane,slot),Source()+size_t(slot)*B,B)==0,"systematic identity");
            }
            if(packet_reference[arm].empty()) packet_reference[arm]=packets[arm];
            CheckOutput(lane,arm,18u,0u); CheckOutput(1u-lane,arm,0u,0u); CheckImmutable();
            for(unsigned family=0;family<2;++family) {
                Prepare(); Handle decoder; decoder.arm=arm; decoder.free_count=&counts.decoder_frees;
                DecodeEvidence& d=s.decode[family];
                d.create=api.decoder_create(&h.state,&decoder.state); ++counts.decoder_creates;
                Require(NoLength(d.create,T::Success) && decoder.state.handle,"preflight decoder create");
                bool success=false;
                for(unsigned step_id=0;step_id<12;++step_id) {
                    const unsigned slot=DecodeSlot(family,step_id);
                    const T::Result value=api.feed(&decoder.state,PacketId(slot),packets[arm].data()+size_t(slot)*B,B);
                    d.feeds[d.count++]=value; ++counts.feeds;
                    Require(FeedResult(value,arm,B),"preflight feed error");
                    if(value.status==T::Success) { success=true; break; }
                }
                Require(success && d.count>=6u,"preflight fixed decoder endpoint");
                d.recover=api.recover(&decoder.state,Output(lane),M); ++counts.recovers;
                Require(RecoveryResult(d.recover,arm,M),"preflight recover result");
                CheckOutput(lane,arm,0u,0u,true); CheckOutput(1u-lane,arm,0u,0u); CheckImmutable();
            }
            for(const Snapshot& earlier:snapshots) if(earlier.arm==arm || (earlier.arm<2u && arm<2u)) {
                if(earlier.arm==arm) Require(earlier.profile==s.profile && Same(earlier.create,s.create),"unstable construction profile");
                for(unsigned f=0;f<2;++f) {
                    Require(earlier.decode[f].count==s.decode[f].count,"unstable first-success prefix");
                    for(unsigned i=0;i<s.decode[f].count;++i)
                        Require(Same(earlier.decode[f].feeds[i],s.decode[f].feeds[i]),"unstable feed statuses");
                    Require(Same(earlier.decode[f].create,s.decode[f].create) &&
                        Same(earlier.decode[f].recover,s.decode[f].recover),"unstable decoder results");
                }
            }
            snapshots.push_back(s);
        }
        Require(snapshots.size()==size_t(reps)*8u,"preflight allocation roster");
        for(unsigned arm=0;arm<4;++arm) packet_hash[arm]=Hash(packets[arm].data(),packets[arm].size());
        CheckHandles();
    }
};

struct Coordinate { unsigned index,rep,order,width,metric,comparison,position,arm; uint64_t q; };
unsigned Side(unsigned position,unsigned order) {
    if(position<2u) return position^order;
    const unsigned pass=(position-2u)/8u,j=((position-2u)%8u)/2u;
    const unsigned first[4]={0,1,1,0};
    return first[j]^pass^((position-2u)%2u)^order;
}
unsigned Phase(unsigned rep,unsigned position) {
    if(position<2u) return rep+12u*(rep%4u);
    const unsigned pass=(position-2u)/8u,j=((position-2u)%8u)/2u;
    return ((rep+6u*pass)%12u)+12u*j;
}
std::vector<Coordinate> Roster() {
    std::vector<Coordinate> result; result.reserve(54432u);
    for(unsigned r=0;r<12;++r) for(unsigned s=0;s<2;++s) for(unsigned ws=0;ws<3;++ws)
        for(unsigned ms=0;ms<6;++ms) for(unsigned cs=0;cs<7;++cs) for(unsigned p=0;p<18;++p) {
            const unsigned order=(r+s)%2u,width=(r+s+ws)%3u,metric=(r+s+ws+ms)%6u;
            const unsigned comparison=(2u*r+s+ws+ms+cs)%7u;
            Coordinate c={unsigned(result.size()),r,order,width,metric,comparison,p,
                kPairs[comparison][Side(p,order)],(uint64_t(2u*Phase(r,p)+1u)*1000000u)/96u};
            result.push_back(c);
        }
    return result;
}
struct WorkResult {
    Counts calls;
    uint64_t handle=kMissing;
    T::Result create={},recover={};
    bool has_create=false,has_recover=false,complete=false;
    std::array<T::Result,18> encode={};
    std::array<T::Result,12> feeds={};
    unsigned encode_count=0,feed_count=0,address_count=0;
    std::array<uint64_t,kBatch> addresses={{0}};
    uint64_t address_min=kMissing,address_max=kMissing;
    std::string address_hash;
};
void FreeInside(Handle& temporary,WorkResult& result) {
    if(temporary.state.handle) {
        const unsigned kind=temporary.state.kind;
        temporary.Reset();
        if(kind==2u) ++result.calls.decoder_frees;
        else ++result.calls.encoder_frees;
    }
}
bool SameProfile(const T::State& a,const T::State& b) {
    return a.profile_bytes==b.profile_bytes && a.profile_bytes<=32u &&
        std::memcmp(a.profile,b.profile,32u)==0;
}
// One out-of-line, label-oblivious implementation for every arm. All scalar
// result checks are nonthrowing: even the first failing call reaches the full
// post-WORK clock/counter bracket. Complete is set only after every cycle.
WH2_COST_NOINLINE void RunWork(Fixture& f,const Coordinate& c,const Snapshot& expected,
                              Handle& temporary,WorkResult& w) noexcept {
    const T::Api& api=*kApis[c.arm]; temporary.arm=c.arm;
    const T::State& retained=f.At(c.rep,c.order,c.arm).state;
    const bool decoder=c.metric>=4u,prebuilt=c.metric==1u || c.metric==2u;
    if(prebuilt || decoder) w.handle=Address(retained.handle);
    try {
        for(unsigned cycle=0;cycle<kBatch;++cycle) {
            const T::State* encoder=&retained;
            if(!prebuilt) {
                w.has_create=true;
                if(decoder) {
                    w.create=api.decoder_create(&retained,&temporary.state); ++w.calls.decoder_creates;
                } else {
                    w.create=api.create(f.lookup.data(),f.lookup.size(),f.Source(),f.M,f.B,&temporary.state);
                    ++w.calls.encoder_creates; encoder=&temporary.state;
                }
                w.addresses[w.address_count++]=Address(temporary.state.handle);
                const T::Result& wanted=decoder?expected.decode[c.metric-4u].create:expected.create;
                if(!Same(w.create,wanted) || !temporary.state.handle || !SameProfile(temporary.state,retained)) {
                    FreeInside(temporary,w); return;
                }
            }
            if(decoder) {
                const unsigned family=c.metric-4u;
                const DecodeEvidence& wanted=expected.decode[family];
                bool success=false; w.feed_count=0;
                for(unsigned step=0;step<12u;++step) {
                    const unsigned slot=DecodeSlot(family,step);
                    const T::Result value=api.feed(&temporary.state,PacketId(slot),
                        f.packets[c.arm].data()+size_t(slot)*f.B,f.B);
                    w.feeds[w.feed_count++]=value; ++w.calls.feeds;
                    if(step>=wanted.count || !Same(value,wanted.feeds[step])) {
                        FreeInside(temporary,w); return;
                    }
                    if(value.status==T::Success) { success=true; break; }
                }
                if(!success || w.feed_count!=wanted.count) { FreeInside(temporary,w); return; }
                w.has_recover=true; w.recover=api.recover(&temporary.state,f.Output(c.order),f.M);
                ++w.calls.recovers;
                if(!Same(w.recover,wanted.recover)) { FreeInside(temporary,w); return; }
            } else {
                const unsigned count=c.metric==0u?0u:(c.metric==3u?18u:6u);
                const unsigned start=c.metric==1u?6u:(c.metric==2u?12u:0u);
                for(unsigned j=0;j<count;++j) {
                    const T::Result value=api.encode(encoder,PacketId(start+j),f.Output(c.order,j),f.B);
                    w.encode[j]=value; w.encode_count=std::max(w.encode_count,j+1u); ++w.calls.encodes;
                    if(!PacketResult(value,f.B)) { FreeInside(temporary,w); return; }
                }
            }
            if(!prebuilt) FreeInside(temporary,w);
        }
        w.complete=true;
    } catch(...) { w.complete=false; }
}
void AddressSummary(WorkResult& result) {
    std::array<uint8_t,8u*kBatch> bytes={{0}};
    for(unsigned i=0;i<result.address_count;++i) {
        const uint64_t address=result.addresses[i];
        result.address_min=result.address_min==kMissing?address:std::min(result.address_min,address);
        result.address_max=result.address_max==kMissing?address:std::max(result.address_max,address);
        for(unsigned j=0;j<8;++j) bytes[8u*i+j]=uint8_t(address>>(8u*j));
    }
    result.address_hash=Hash(bytes.data(),8u*result.address_count);
}
void WorkJson(std::ostream& out,const WorkResult& w) {
    const Counts& c=w.calls;
    out << "{\"create_calls\":" << c.encoder_creates << ",\"encode_calls\":" << c.encodes
        << ",\"free_calls\":" << c.encoder_frees << ",\"decoder_create_calls\":" << c.decoder_creates
        << ",\"feed_calls\":" << c.feeds << ",\"recover_calls\":" << c.recovers
        << ",\"decoder_free_calls\":" << c.decoder_frees << ",\"create\":";
    if(w.has_create) ResultJson(out,w.create); else out << "null";
    out << ",\"recover\":"; if(w.has_recover) ResultJson(out,w.recover); else out << "null";
    out << ",\"handle\":"; Number(out,w.handle); out << ",\"encode\":[";
    for(unsigned i=0;i<w.encode_count;++i) { if(i) out << ','; ResultJson(out,w.encode[i]); }
    out << "],\"feeds\":[";
    for(unsigned i=0;i<w.feed_count;++i) { if(i) out << ','; ResultJson(out,w.feeds[i]); }
    out << "],\"complete\":" << (w.complete?"true":"false") << ",\"addresses\":{\"count\":" << w.address_count
        << ",\"sha256\":" << Quote(w.address_hash) << ",\"min\":"; Number(out,w.address_min);
    out << ",\"max\":"; Number(out,w.address_max); out << "}}";
}
struct Record {
    Coordinate c;
    uint64_t target=kMissing,ready=kMissing;
    std::array<uint64_t,4> wait={{kMissing,kMissing,kMissing,kMissing}};
    Observation observation;
    WorkResult work;
    bool called=false,checked=false;
    explicit Record(const Coordinate& coordinate):c(coordinate) {}
};
void RecordJson(std::ostream& out,const Record& r) {
    const Coordinate& c=r.c;
    out << "{\"type\":\"record\",\"index\":" << c.index << ",\"rep\":" << c.rep << ",\"order\":" << c.order
        << ",\"width\":" << c.width << ",\"metric\":" << c.metric << ",\"class\":" << c.comparison
        << ",\"position\":" << c.position << ",\"arm\":" << c.arm << ",\"q\":" << c.q << ",\"target\":";
    Number(out,r.target); out << ",\"ready\":"; Number(out,r.ready); out << ",\"wait\":[";
    for(unsigned i=0;i<4;++i) { if(i) out << ','; Number(out,r.wait[i]); }
    out << "],\"observation\":"; ObservationJson(out,r.observation);
    out << ",\"called\":" << (r.called?"true":"false") << ",\"checked\":" << (r.checked?"true":"false")
        << ",\"work\":"; WorkJson(out,r.work); out << '}';
}
void CheckWork(Fixture& f,Record& r,const Snapshot& expected,const Handle& temporary) {
    const auto& c=r.c; const auto& w=r.work; const Counts& n=w.calls;
    const bool decoder=c.metric>=4u,prebuilt=c.metric==1u || c.metric==2u;
    const unsigned count=decoder || c.metric==0u?0u:(c.metric==3u?18u:6u);
    Require(w.complete && !temporary.state.handle,"incomplete work or leaked lifecycle handle");
    Require(w.encode_count==count && n.encodes==kBatch*count,"encode accounting");
    Require(n.encoder_creates==(!decoder && !prebuilt?kBatch:0u) && n.encoder_frees==n.encoder_creates &&
            n.decoder_creates==(decoder?kBatch:0u) && n.decoder_frees==n.decoder_creates &&
            w.address_count==(!prebuilt?kBatch:0u),"lifecycle accounting");
    Require(w.has_create==!prebuilt && (!w.has_create ||
            Same(w.create,decoder?expected.decode[c.metric-4u].create:expected.create)),"fresh creation result");
    Require((prebuilt || decoder)?w.handle==Address(f.At(c.rep,c.order,c.arm).state.handle):w.handle==kMissing,
            "persistent handle identity");
    for(unsigned i=0;i<w.encode_count;++i) Require(PacketResult(w.encode[i],f.B),"packet result aggregate");
    if(decoder) {
        const auto& wanted=expected.decode[c.metric-4u];
        Require(w.feed_count==wanted.count && n.feeds==kBatch*wanted.count &&
                n.recovers==kBatch && w.has_recover && Same(w.recover,wanted.recover),"decoder accounting");
        for(unsigned i=0;i<w.feed_count;++i) Require(Same(w.feeds[i],wanted.feeds[i]),"feed prefix aggregate");
    } else Require(w.feed_count==0u && n.feeds==0u && n.recovers==0u && !w.has_recover,"nondecoder accounting");
    for(unsigned i=0;i<w.address_count;++i) Require(w.addresses[i]>0,"null ephemeral handle");
    f.CheckOutput(c.order,c.arm,count,c.metric==1u?6u:(c.metric==2u?12u:0u),decoder);
    f.CheckOutput(1u-c.order,c.arm,0u,0u); f.CheckImmutable(); r.checked=true;
}
template<class R> void Wait(R& reader,uint64_t target,uint64_t start,
                            std::array<uint64_t,4>& times) {
    Require(reader.Mono(times[0]) && reader.Cpu(times[1]),"wait start clocks");
    uint64_t now=times[0];
    while(now<target) {
        Require(reader.Mono(now),"wait monotonic clock");
        Require(now>=start && now-start<kWall,"worker wall deadline during wait");
    }
    Require(reader.Cpu(times[3]) && reader.Mono(times[2]),"wait end clocks");
    Require(times[0]<=times[2] && times[1]<=times[3] && times[2]>=target,"wait ordering");
}
void Budget(Reader& reader,uint64_t start,uint64_t cpu_start,uint64_t work) {
    const uint64_t now=reader.Now(),cpu=reader.Thread();
    Require(now>=start && now-start<kWall,"worker wall cap");
    Require(cpu>=cpu_start && cpu-cpu_start<=kCpu,"worker CPU cap");
    Require(work<=kWork,"cumulative inner work cap");
}
uint64_t Target(uint64_t t0,const Coordinate& c) {
    const uint64_t offset=UINT64_C(2000000)*c.index+c.q;
    Require(t0<kMissing-offset,"schedule overflow"); return t0+offset;
}
void CheckReady(uint64_t ready,uint64_t target,const Observation& previous) {
    Require(ready!=kMissing && target>=100000u && previous.m3<=ready && ready<=target-100000u,
            "preparation/previous cleanup missed slot deadline");
}
void CheckStart(const Record& record) {
    Require(record.observation.m1>=record.target && record.observation.m1-record.target<=5000u,
            "codec start lateness exceeds 5us");
}


void FixtureJson(std::ostream& out,const Fixture& f) {
    out << "{\"width_index\":" << f.width_index << ",\"block_bytes\":" << f.B << ",\"message_bytes\":" << f.M
        << ",\"source_hex\":" << Quote(Hex(f.Source(),f.M)) << ",\"source_sha256\":" << Quote(f.source_hash)
        << ",\"source_address\":" << Address(f.Source()) << ",\"output_addresses\":[" << Address(f.outputs[0].data)
        << ',' << Address(f.outputs[1].data) << "],\"output_stride\":" << f.stride << ",\"packets_hex\":[";
    for(unsigned arm=0;arm<4;++arm) { if(arm) out << ','; out << Quote(Hex(f.packets[arm].data(),f.packets[arm].size())); }
    out << "],\"packets_sha256\":[";
    for(unsigned arm=0;arm<4;++arm) { if(arm) out << ','; out << Quote(f.packet_hash[arm]); }
    out << "],\"rows_hex\":" << Quote(Hex(f.rows.data(),f.rows.size())) << ",\"handles\":[";
    for(size_t i=0;i<f.snapshots.size();++i) { if(i) out << ','; SnapshotJson(out,f.snapshots[i]); }
    out << "]}";
}
void HeaderJson(std::ostream& out,const std::array<Fixture*,3>& fixtures,const std::vector<uint8_t>& lookup,
                const std::string& claim,uint64_t start,uint64_t cpu_start,const std::string& identity,
                const std::string& runtime,const std::array<uint64_t,2>& resolutions,
                const Observation& prelude,uint64_t final_state) {
    out << "{\"type\":\"header\",\"protocol\":" << Quote(kProtocol)
        << ",\"schema\":" << Quote(std::string(kProtocol)+".raw.v1") << ",\"claim_sha256\":" << Quote(claim)
        << ",\"target_cpu\":50,\"worker_start_ns\":" << start << ",\"worker_start_cpu_ns\":" << cpu_start
        << ",\"identity_before\":" << identity << ",\"runtime_before\":" << runtime
        << ",\"clock_resolutions\":[" << resolutions[0] << ',' << resolutions[1] << ']'
        << ",\"lookup_address\":" << Address(lookup.data()) << ",\"lookup_sha256\":" << Quote(Hash(lookup.data(),lookup.size()))
        << ",\"fixtures\":[";
    for(unsigned i=0;i<3;++i) { if(i) out << ','; FixtureJson(out,*fixtures[i]); }
    out << "],\"prelude\":{\"iterations\":" << kPreludeIterations << ",\"seed\":" << kPreludeSeed
        << ",\"final_state\":" << final_state << ",\"observation\":";
    ObservationJson(out,prelude); out << "}}";
}
struct Accounting {
    uint64_t records=0,callbacks=0,checked=0,work=0,wait_wall=0,wait_cpu=0;
    Counts calls;
    void Add(const Record& r) {
        ++records; callbacks+=r.called?1u:0u; checked+=r.checked?1u:0u;
        calls.encoder_creates+=r.work.calls.encoder_creates; calls.encoder_frees+=r.work.calls.encoder_frees;
        calls.encodes+=r.work.calls.encodes; calls.decoder_creates+=r.work.calls.decoder_creates;
        calls.decoder_frees+=r.work.calls.decoder_frees; calls.feeds+=r.work.calls.feeds; calls.recovers+=r.work.calls.recovers;
        const Observation& o=r.observation;
        if(o.m1!=kMissing && o.m2!=kMissing && o.m2>=o.m1) work+=o.m2-o.m1;
        if(r.wait[0]!=kMissing && r.wait[2]!=kMissing && r.wait[2]>=r.wait[0]) wait_wall+=r.wait[2]-r.wait[0];
        if(r.wait[1]!=kMissing && r.wait[3]!=kMissing && r.wait[3]>=r.wait[1]) wait_cpu+=r.wait[3]-r.wait[1];
    }
};
void FooterJson(std::ostream& out,const Accounting& a,const Counts& p,const std::string& failure,
                uint64_t schedule_now,uint64_t t0,uint64_t end,uint64_t cpu_end,
                const std::string& identity,const std::string& runtime) {
    const Counts& n=a.calls;
    out << "{\"type\":\"footer\",\"protocol\":" << Quote(kProtocol)
        << ",\"outcome\":" << Quote(failure.empty()?"COMPLETE":"INVALID")
        << ",\"failure\":" << (failure.empty()?"null":Quote(failure)) << ",\"schedule_now_ns\":";
    Number(out,schedule_now); out << ",\"t0\":"; Number(out,t0);
    out << ",\"records\":" << a.records << ",\"callbacks\":" << a.callbacks << ",\"checked_callbacks\":" << a.checked
        << ",\"create_calls\":" << n.encoder_creates << ",\"encode_calls\":" << n.encodes << ",\"free_calls\":" << n.encoder_frees
        << ",\"decoder_create_calls\":" << n.decoder_creates << ",\"feed_calls\":" << n.feeds
        << ",\"recover_calls\":" << n.recovers << ",\"decoder_free_calls\":" << n.decoder_frees
        << ",\"sum_work_ns\":" << a.work << ",\"sum_wait_wall_ns\":" << a.wait_wall
        << ",\"sum_wait_cpu_ns\":" << a.wait_cpu << ",\"worker_end_ns\":"; Number(out,end);
    out << ",\"worker_end_cpu_ns\":"; Number(out,cpu_end);
    out << ",\"identity_after\":" << (identity.empty()?"null":identity)
        << ",\"runtime_after\":" << (runtime.empty()?"null":runtime)
        << ",\"preflight\":{\"encoder_creates\":" << p.encoder_creates << ",\"encoder_frees\":" << p.encoder_frees
        << ",\"encodes\":" << p.encodes << ",\"decoder_creates\":" << p.decoder_creates << ",\"feeds\":" << p.feeds
        << ",\"recovers\":" << p.recovers << ",\"decoder_frees\":" << p.decoder_frees << "}}";
}
std::string ReadPinned(const char* path,size_t maximum,const std::string& expected) {
    Require(expected.size()==64u && std::all_of(expected.begin(),expected.end(),[](char c) {
        return (c>='0' && c<='9') || (c>='a' && c<='f'); }),"SHA256 argument");
    const int fd=open(path,O_RDONLY|O_NOFOLLOW|O_NONBLOCK|O_CLOEXEC);
    Require(fd>=0,"pinned file absent");
    struct Close { int fd; ~Close(){close(fd);} } cleanup={fd};
    struct stat before={},after={};
    Require(fstat(fd,&before)==0 && S_ISREG(before.st_mode) && before.st_nlink==1 &&
            before.st_size>0 && uint64_t(before.st_size)<=maximum,"pinned regular bounded file");
    std::string data(size_t(before.st_size),'\0'); size_t done=0;
    while(done<data.size()) {
        const ssize_t got=read(fd,&data[done],data.size()-done);
        if(got<0 && errno==EINTR) continue;
        Require(got>0,"pinned file short read"); done+=size_t(got);
    }
    char extra=0; Require(read(fd,&extra,1)==0 && fstat(fd,&after)==0,"pinned file grew");
    Require(before.st_dev==after.st_dev && before.st_ino==after.st_ino && before.st_size==after.st_size &&
            before.st_nlink==after.st_nlink && before.st_mtim.tv_sec==after.st_mtim.tv_sec &&
            before.st_mtim.tv_nsec==after.st_mtim.tv_nsec && before.st_ctim.tv_sec==after.st_ctim.tv_sec &&
            before.st_ctim.tv_nsec==after.st_ctim.tv_nsec && Sha256Hex(data)==expected,"pinned file changed/hash mismatch");
    return data;
}
void Authenticate(const std::string& expected) {
    const std::string data=ReadPinned(kClaim,1024u*1024u,expected);
    Require(data.find(std::string("\"protocol\":\"")+kProtocol+"\"")!=std::string::npos,"claim protocol absent");
    struct stat output={},error={};
    Require(fstat(STDOUT_FILENO,&output)==0 && S_ISFIFO(output.st_mode) &&
            fstat(STDERR_FILENO,&error)==0 && S_ISFIFO(error.st_mode),"worker requires captured output pipes");
}
void Limits() {
    const rlimit cpu={120,120},memory={512u*1024u*1024u,512u*1024u*1024u},core={0,0};
    Require(setrlimit(RLIMIT_CPU,&cpu)==0 && setrlimit(RLIMIT_AS,&memory)==0 &&
            setrlimit(RLIMIT_CORE,&core)==0,"worker resource limits");
}

int Worker(const std::string& claim) {
    Reader reader; Accounting accounting; Counts preflight_counts;
    uint64_t start=kMissing,cpu_start=kMissing,schedule_now=kMissing,t0=kMissing,end=kMissing,cpu_end=kMissing;
    std::string failure,identity_before,runtime_before,identity_after,runtime_after;
    try {
        Require(WH2_MATVEC_COST_SANITIZERS==0,"sanitizer build cannot run scientific worker");
        Authenticate(claim); Limits(); start=reader.Now(); cpu_start=reader.Thread(); Pin();
        Require(wirehair_init()==Wirehair_Success && gf256_init()==0,"shared GF initialization"); InitializeApis();
        identity_before=Identity(); runtime_before=Runtime();
        Require(wh2_matvec_gfni_r0::Available(),"qualified GFNI matvec required");
        const std::array<uint64_t,2> resolutions={{reader.Resolution(CLOCK_MONOTONIC),reader.Resolution(CLOCK_THREAD_CPUTIME_ID)}};
        Require(resolutions[0]==1u && resolutions[1]==1u,"frozen clock resolutions");
        const std::string packed=ReadPinned(WH2_MATVEC_COST_LOOKUP_PATH,kLookupBytes,kLookupHash);
        Require(packed.size()==kLookupBytes,"scientific lookup exact extent");
        const std::vector<uint8_t> lookup(packed.begin(),packed.end()),lookup_reference=lookup;
        const uint8_t fixed_feedback[2][6]={{124,127,152,84,241,63},{125,127,152,84,241,63}};
        const Independent oracle(fixed_feedback);
        {
            Fixture f0(0u,2u,12u,lookup,oracle),f1(1u,64u,12u,lookup,oracle),f2(2u,1280u,12u,lookup,oracle);
            const std::array<Fixture*,3> fixtures={{&f0,&f1,&f2}};
            uint64_t preflight_feeds=0;
            for(Fixture* f : fixtures) {
                f->Preflight(preflight_counts);
                for(const auto& s:f->snapshots) for(const auto& d:s.decode) preflight_feeds+=d.count;
                Budget(reader,start,cpu_start,0u);
            }
            Require(preflight_counts.feeds==preflight_feeds && lookup==lookup_reference,"preflight feeds/lookup");
            const auto roster=Roster(); Require(roster.size()==54432u,"complete roster size");
            uint64_t expected_feeds=0;
            for(const auto& c:roster) if(c.metric>=4u)
                expected_feeds+=uint64_t(kBatch)*fixtures[c.width]->Expected(c.rep,c.order,c.arm).decode[c.metric-4u].count;
            Budget(reader,start,cpu_start,0u);
            Observation prelude; uint64_t final_state=0;
            const char* error=Capture(reader,[&] { final_state=Prelude(kPreludeSeed,kPreludeIterations); },prelude);
            Require(!error,error?error:"prelude capture"); Validate(prelude,nullptr);
            Require(final_state==kPreludeFinal,"prelude checksum");
            HeaderJson(std::cout,fixtures,lookup,claim,start,cpu_start,identity_before,runtime_before,resolutions,prelude,final_state);
            std::cout << '\n' << std::flush; Require(bool(std::cout),"header stream write");
            schedule_now=reader.Now();
            Require(schedule_now/2000000u<(kMissing/2000000u)-2u,"schedule origin overflow");
            t0=(schedule_now/2000000u+2u)*2000000u;
            Observation previous=prelude;
            for(const Coordinate& c : roster) {
                Fixture& f=*fixtures[c.width]; const Snapshot& expected=f.Expected(c.rep,c.order,c.arm);
                Record record(c); Handle temporary;
                try {
                    record.target=Target(t0,c); f.Prepare();
                    Budget(reader,start,cpu_start,accounting.work);
                    record.ready=reader.Now(); CheckReady(record.ready,record.target,previous);
                    Wait(reader,record.target,start,record.wait);
                    error=Capture(reader,[&] { RunWork(f,c,expected,temporary,record.work); record.called=true; },record.observation);
                    Require(!error,error?error:"callback capture"); Validate(record.observation,&previous); CheckStart(record);
                    for(const Fixture* item:fixtures) item->CheckImmutable();
                    Require(lookup==lookup_reference,"immutable global lookup changed");
                    CheckWork(f,record,expected,temporary); AddressSummary(record.work);
                    Budget(reader,start,cpu_start,accounting.work+record.observation.m2-record.observation.m1);
                    Require(OnCpu(),"worker migrated/affinity changed"); previous=record.observation;
                } catch(...) {
                    AddressSummary(record.work); accounting.Add(record); RecordJson(std::cout,record);
                    std::cout << '\n' << std::flush; throw;
                }
                accounting.Add(record); RecordJson(std::cout,record); std::cout << '\n' << std::flush;
                Require(bool(std::cout),"record stream write");
            }
            const Counts& n=accounting.calls;
            Require(accounting.records==54432u && accounting.callbacks==54432u && accounting.checked==54432u &&
                    n.encoder_creates==2322432u && n.encoder_frees==n.encoder_creates && n.encodes==34836480u &&
                    n.decoder_creates==2322432u && n.decoder_frees==n.decoder_creates &&
                    n.recovers==n.decoder_creates && n.feeds==expected_feeds,"complete callback/API accounting");
            for(const Fixture* f:fixtures) { f->CheckImmutable(); f->CheckHandles(); }
            Require(lookup==lookup_reference && OnCpu(),"final lookup/affinity");
            Require(preflight_counts.feeds==preflight_feeds,"preflight feed totals changed");
            // All 288 retained encoders are destroyed here, within worker cost
            // but outside repair-component timing.
        }
        Require(preflight_counts.encoder_creates==288u && preflight_counts.encoder_frees==288u &&
                preflight_counts.encodes==5184u && preflight_counts.decoder_creates==576u &&
                preflight_counts.recovers==576u && preflight_counts.decoder_frees==576u,"complete preflight accounting");
        identity_after=Identity(); runtime_after=Runtime();
        Require(identity_after==identity_before && runtime_after==runtime_before &&
            wh2_matvec_gfni_r0::Available(),"runtime/physical identity changed");
        Budget(reader,start,cpu_start,accounting.work);
    } catch(const std::exception& error) { failure=error.what(); }
      catch(...) { failure="unexpected worker exception"; }
    if(start!=kMissing) {
        if(!reader.Cpu(cpu_end) || !reader.Mono(end)) {
            if(failure.empty()) failure="footer clock capture";
        } else if(end<start || end-start>=kWall || cpu_end<cpu_start || cpu_end-cpu_start>kCpu ||
                  cpu_end-cpu_start>end-start) {
            if(failure.empty()) failure="final cleanup resource cap";
        }
    }
    FooterJson(std::cout,accounting,preflight_counts,failure,schedule_now,t0,end,cpu_end,identity_after,runtime_after);
    std::cout << '\n' << std::flush; return failure.empty() && bool(std::cout)?0:1;
}

// Neutral-only readers and entrypoints below cannot access scientific files,
// affinity, real clocks, the selected companion pair or the one-shot Worker.
struct FakeReader {
    unsigned calls=0,fail=UINT32_MAX; uint64_t mono=1000,cpu=1000000;
    bool Mono(uint64_t& value) { if(calls++==fail) return false; value=mono; mono+=100; return true; }
    bool Cpu(uint64_t& value) { if(calls++==fail) return false; value=cpu; cpu+=200; return true; }
    bool Usage(Counters& value) { if(calls++==fail) return false; value={{0,0,0,0}}; return true; }
};
template<class F> void Rejected(F function,const char* message) {
    bool rejected=false;
    try { function(); } catch(const std::runtime_error&) { rejected=true; }
    Require(rejected,message);
}
void NeutralClocks() {
    uint64_t mono=1000u,cpu=1000000u;
    Require(!OrderedClock(999u,mono) && mono==1000u &&
        !OrderedClock(999999u,cpu) && cpu==1000000u,"neutral independent clock regressions");
    Require(OrderedClock(1000u,mono) && OrderedClock(1000001u,cpu),"neutral independent clock domains");
    FakeReader reader; Observation o; unsigned work=0;
    Require(Capture(reader,[&]{++work;},o)==nullptr && work==1u,"neutral full capture");
    Validate(o,nullptr);
    for(unsigned failure=0;failure<8u;++failure) {
        FakeReader bad; bad.fail=failure; Observation partial; unsigned completed=0;
        Require(Capture(bad,[&]{++completed;},partial)!=nullptr,"neutral failed capture");
        Require(completed==(failure>=4u?1u:0u),"neutral partial callback boundary");
    }
    Observation later=o; later.m0=1u;
    Rejected([&]{Validate(later,&o);},"neutral reverse observations");
    later=o; later.c1=later.c0+(later.m3-later.m0)+1u;
    Rejected([&]{Validate(later,nullptr);},"neutral CPU exceeds enclosing wall interval");
    later=o; later.ru1[0]=kMissing;
    Rejected([&]{Validate(later,nullptr);},"neutral missing counter");
    timespec time={}; uint64_t ns=0;
    time.tv_sec=-1; Rejected([&]{Require(Nanoseconds(time,ns),"negative");},"neutral negative clock");
    time.tv_sec=0; time.tv_nsec=1000000000L;
    Require(!Nanoseconds(time,ns),"neutral nanosecond extent");
    const Coordinate c={0,0,0,0,0,0,0,0,0};
    Record r(c); r.target=10000000u; r.observation.m1=r.target+5000u; CheckStart(r);
    ++r.observation.m1;
    Rejected([&]{CheckStart(r);},"neutral lateness limit");
    Observation previous; previous.m3=1u;
    CheckReady(r.target-100000u,r.target,previous);
    Rejected([&]{CheckReady(r.target-99999u,r.target,previous);},"neutral preparation limit");
    Require(Prelude(kPreludeSeed,1u)==UINT64_C(0xdc1b77ae0bf34dad) &&
        Prelude(kPreludeSeed,2u)==UINT64_C(0x64f0eeb9026e6076),"neutral short prelude");
}
void NeutralRoster() {
    const auto roster=Roster(); Require(roster.size()==54432u,"neutral roster size");
    unsigned seen[3][6][7][2][12]={},phases[3][6][7][2][48]={};
    const unsigned sides[18]={0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
    uint64_t minimum=kMissing,maximum=0;
    for(size_t i=0;i<roster.size();++i) {
        const auto& c=roster[i];
        Require(c.index==i && c.width<3u && c.arm<4u && c.metric<6u &&
            c.comparison<7u && c.order<2u && c.rep<12u && c.position<18u,"neutral coordinate bounds");
        Require(Side(c.position,c.order)==(sides[c.position]^c.order),"neutral literal side order");
        Require(c.arm==kPairs[c.comparison][sides[c.position]^c.order],"neutral actual arm mapping");
        const unsigned bin=c.position<2u?c.rep+12u*(c.rep%4u):
            ((c.rep+(((c.position-2u)/8u)==0u?0u:6u))%12u)+12u*(((c.position-2u)%8u)/2u);
        Require(c.q==((2u*uint64_t(bin)+1u)*1000000u)/96u,"neutral literal phase mapping");
        if(c.position==0u) ++seen[c.width][c.metric][c.comparison][c.order][c.rep];
        if(c.position>=2u && c.position%2u==0u) {
            const auto& next=roster[i+1u];
            Require(next.q==c.q && next.position==c.position+1u &&
                Side(next.position,next.order)!=Side(c.position,c.order),"neutral paired phases");
            ++phases[c.width][c.metric][c.comparison][c.order][bin];
        }
        if(i) {
            const uint64_t gap=Target(10000000u,c)-Target(10000000u,roster[i-1u]);
            minimum=std::min(minimum,gap); maximum=std::max(maximum,gap);
        }
    }
    for(unsigned w=0;w<3u;++w) for(unsigned m=0;m<6u;++m)
        for(unsigned c=0;c<7u;++c) for(unsigned order=0;order<2u;++order) {
            for(unsigned rep=0;rep<12u;++rep)
                Require(seen[w][m][c][order][rep]==1u,"neutral complete panel coverage");
            for(unsigned phase=0;phase<48u;++phase)
                Require(phases[w][m][c][order][phase]==2u,"neutral twice-covered phase catalog");
        }
    Require(minimum>=1125000u && maximum<=2875000u,"neutral fixed slot bounds");
    const uint32_t far[6]={4294967295u,4294967293u,4294967291u,4294967289u,4294967287u,4294967285u};
    for(unsigned j=0;j<6u;++j) {
        Require(PacketId(12u+j)==far[j] && PacketId(j)==j && PacketId(6u+j)==6u+j,
            "neutral packet roster");
        Require(DecodeSlot(0,j)==6u+j && DecodeSlot(1,j)==12u+j &&
            DecodeSlot(0,6u+j)==j && DecodeSlot(1,6u+j)==j,"neutral guaranteed systematic suffix");
    }
}
void NeutralBridges(Fixture& f) {
    for(unsigned arm=0;arm<4u;++arm) {
        const auto& api=*kApis[arm]; Handle& h=f.At(0,0,arm); const void* address=h.state.handle;
        Require(api.create(f.lookup.data(),f.lookup.size(),f.Source(),f.M,f.B,&h.state).status==T::Failure &&
            h.state.handle==address,"neutral occupied create");
        Require(api.create(f.lookup.data(),f.lookup.size(),f.Source(),f.M,f.B,nullptr).status==T::Failure,
            "neutral null create");
        Require(api.encode(nullptr,6u,f.Output(0),f.B).status==T::Failure,"neutral null encoder");
        f.Prepare();
        const auto short_result=api.encode(&h.state,6u,f.Output(0),f.B-1u);
        Require(short_result.status==T::Failure && short_result.written==0u,"neutral short encode");
        f.CheckOutput(0,arm,0,0); f.CheckOutput(1,arm,0,0);
        Require(api.encode(&h.state,6u,const_cast<uint8_t*>(f.Source()),f.B).status==T::Failure,
            "neutral source alias");
        Require(api.decoder_create(&h.state,nullptr).status==T::Failure,"neutral null decoder");
        f.CheckImmutable(); api.encoder_free(nullptr); api.decoder_free(nullptr);
    }
}
void NeutralMetrics(Fixture& f) {
    for(unsigned metric=0;metric<6u;++metric) for(unsigned arm=0;arm<4u;++arm) for(unsigned lane=0;lane<2u;++lane) {
        const Coordinate c={0,0,lane,f.width_index,metric,0,0,arm,0};
        const Snapshot& expected=f.Expected(0,lane,arm);
        Record r(c); Handle temporary; f.Prepare(); FakeReader reader;
        Require(Capture(reader,[&] { RunWork(f,c,expected,temporary,r.work); r.called=true; },r.observation)==nullptr,
            "neutral metric fake capture");
        Validate(r.observation,nullptr); CheckWork(f,r,expected,temporary); AddressSummary(r.work);
        Require(r.checked,"neutral actual metric bytes");
        r.work.complete=false;
        Rejected([&]{CheckWork(f,r,expected,temporary);},"neutral completed flag cannot be masked");
        r.work.complete=true;
        f.outputs[1u-lane].data[0]^=1u;
        Rejected([&]{f.CheckOutput(1u-lane,arm,0,0);},"neutral cross-lane guard");
        f.outputs[1u-lane].data[0]^=1u;
        f.packet_reference[arm][0]^=1u;
        Rejected([&]{f.CheckImmutable();},"neutral packet corpus mutation");
        f.packet_reference[arm][0]^=1u;
        std::ostringstream out; RecordJson(out,r); Require(out.str().size()<4096u,"neutral bounded record");
    }
    f.CheckHandles();
}
namespace neutral_fault {
const T::Api* real=nullptr;
unsigned calls=0,fail=0;
T::Result Encode(const T::State* state,uint32_t id,void* output,size_t capacity) {
    T::Result r=real->encode(state,id,output,capacity);
    if(++calls==fail) ++r.written;
    return r;
}
T::Result Feed(T::State* state,uint32_t id,const void* input,size_t bytes) {
    T::Result r=real->feed(state,id,input,bytes);
    if(++calls==fail) ++r.required;
    return r;
}
}
void NeutralFailedLastCall(Fixture& f) {
    using namespace neutral_fault;
    for(unsigned metric: {3u,5u}) {
        const Coordinate c={0,0,0,f.width_index,metric,0,0,0,0};
        const Snapshot& expected=f.Expected(0,0,0);
        real=kApis[0]; T::Api probe=*real; calls=0;
        if(metric==3u) { probe.encode=Encode; fail=kBatch*18u; }
        else { probe.feed=Feed; fail=kBatch*expected.decode[1].count; }
        struct Restore {
            const T::Api* saved;
            ~Restore() { kApis[0]=saved; }
        } restore={kApis[0]};
        kApis[0]=&probe;
        Record r(c); Handle temporary; f.Prepare(); FakeReader reader;
        Require(Capture(reader,[&] { RunWork(f,c,expected,temporary,r.work); r.called=true; },r.observation)==nullptr,
            "neutral failed-last-call post-clock capture");
        Validate(r.observation,nullptr);
        Require(calls==fail && !r.work.complete && !temporary.state.handle &&
            r.work.address_count==kBatch,"neutral failed-last-call cleanup and retained count");
        Rejected([&]{CheckWork(f,r,expected,temporary);},"neutral failed-last-call aggregate rejection");
    }
}
// A conservative successful-schema bound, including 12-feed prefixes for all
// decoders and INT64_MAX-sized observational fields. No codec work occurs here.
uint64_t NeutralSerialization() {
    const auto roster=Roster(); uint64_t bytes=2u*1024u*1024u;
    const uint64_t large=uint64_t(INT64_MAX);
    for(const Coordinate& c:roster) {
        Record r(c); WorkResult& w=r.work;
        r.target=large; r.ready=large; r.wait.fill(large);
        r.observation.m0=large; r.observation.c0=large; r.observation.m1=large;
        r.observation.m2=large; r.observation.c1=large; r.observation.m3=large;
        r.observation.ru0.fill(large); r.observation.ru1.fill(large);
        r.called=true; r.checked=true; w.complete=true;
        const unsigned B=c.width==0u?2u:c.width==1u?64u:1280u;
        const bool decoder=c.metric>=4u,prebuilt=c.metric==1u || c.metric==2u;
        w.has_create=!prebuilt; w.has_recover=decoder;
        w.handle=(prebuilt || decoder)?large:kMissing;
        w.address_count=prebuilt?0u:kBatch; w.address_hash=std::string(64u,'f');
        if(!prebuilt) { w.address_min=large; w.address_max=large; }
        w.create=T::Result{T::Success,0,0,0,T::NoLength};
        w.recover=T::Result{T::Success,0,6u*B,6u*B,c.arm==2u?T::InferredLength:T::ReturnedLength};
        w.encode_count=c.metric==3u?18u:(prebuilt?6u:0u);
        w.feed_count=decoder?12u:0u;
        for(unsigned i=0;i<w.encode_count;++i)
            w.encode[i]=T::Result{T::Success,0,B,B,T::ReturnedLength};
        for(unsigned i=0;i<w.feed_count;++i)
            w.feeds[i]=T::Result{i==11u?T::Success:T::NeedMore,i==11u?0:1,B,0,c.arm<2u?T::ReturnedLength:T::InferredLength};
        w.calls.encoder_creates=(!prebuilt && !decoder)?kBatch:0u;
        w.calls.encoder_frees=w.calls.encoder_creates;
        w.calls.decoder_creates=decoder?kBatch:0u; w.calls.decoder_frees=w.calls.decoder_creates;
        w.calls.encodes=kBatch*w.encode_count; w.calls.feeds=kBatch*w.feed_count;
        w.calls.recovers=decoder?kBatch:0u;
        std::ostringstream out; RecordJson(out,r); bytes+=out.str().size()+1u;
    }
    Require(bytes<=96u*1024u*1024u,"neutral complete worst successful stream exceeds raw cap");
    return bytes;
}
int Neutral(bool expect_scalar) {
    Require(wirehair_init()==Wirehair_Success && gf256_init()==0,"neutral shared runtime");
    Require(wh2_matvec_gfni_r0::Available()!=expect_scalar,"neutral expected ISA dispatch");
    InitializeApis(); NeutralClocks(); NeutralRoster();
    const std::vector<uint8_t> lookup=NeutralLookup(),saved=lookup;
    const uint8_t feedback[2][6]={{1,1,0,1,0,1},{3,1,0,1,0,1}};
    const Independent oracle(feedback);
    for(unsigned B: {2u,65u}) {
        Counts counts; uint64_t feeds=0;
        {
            Fixture f(0u,B,1u,lookup,oracle); f.Preflight(counts);
            for(const Snapshot& s:f.snapshots) for(const auto& d:s.decode) feeds+=d.count;
            NeutralBridges(f); NeutralMetrics(f); NeutralFailedLastCall(f);
        }
        Require(counts.encoder_creates==8u && counts.encoder_frees==8u && counts.encodes==144u &&
            counts.decoder_creates==16u && counts.decoder_frees==16u && counts.recovers==16u &&
            counts.feeds==feeds,"neutral complete preflight accounting");
        Require(lookup==saved,"neutral lookup unchanged");
    }
    const uint64_t bytes=NeutralSerialization();
    std::cout << "PASS neutral four arms B2/65, six full metric bodies, fake clocks/phase catalogs, "
        << "first-bad result capture; conservative raw bound " << bytes << " bytes\n";
    return 0;
}
} // namespace cost

int main(int argc,char** argv) {
    try {
        if(argc==2 && std::strcmp(argv[1],"--neutral")==0) return cost::Neutral(false);
        if(argc==3 && std::strcmp(argv[1],"--neutral")==0 &&
            std::strcmp(argv[2],"--expect-scalar")==0) return cost::Neutral(true);
        if(argc==3 && std::strcmp(argv[1],"--worker")==0) return cost::Worker(argv[2]);
        std::cerr << "usage: --neutral [--expect-scalar] | --worker CLAIM_SHA256 (no default execution)\n";
        return 2;
    } catch(const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n'; return 1;
    }
}
