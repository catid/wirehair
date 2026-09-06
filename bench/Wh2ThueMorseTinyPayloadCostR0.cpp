// .69 frozen, one-shot tiny-payload cost screen. No production source is modified.
#include "Wh2ThueMorseTinyPayloadCostR0Bridge.h"
#include "Wh2FrozenTrace.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "gf256.h"

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
#ifndef WH2_TINY_COST_SANITIZERS
#error "Builder must identify sanitizer versus native execution"
#endif
#ifndef WH2_TINY_COST_LOOKUP_PATH
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
static const char kProtocol[] = "wirehair.wh2.thue-tiny-payload-cost-r0";
static const char kClaim[] = "/var/tmp/wh2-thue-tiny-payload-cost-r0/CLAIM.json";
static const uint64_t kMissing = UINT64_MAX;
static const int64_t kAbsent = INT64_MIN;
static const uint64_t kWall = UINT64_C(15000000000);
static const uint64_t kCpu = UINT64_C(14000000000);
static const uint64_t kWork = UINT64_C(1000000000);
static const uint64_t kPreludeSeed = UINT64_C(0x9e3779b97f4a7c15);
static const uint64_t kPreludeFinal = UINT64_C(0x43935dad1647741b);
static const unsigned kPreludeIterations = 1u << 20;
static const unsigned kSides[10] = {0,1,0,1,1,0,1,0,0,1};
static const unsigned kPairs[3][2] = {{0,0},{1,1},{0,1}};
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

namespace T = wh2_tiny_payload_cost_r0;
const T::Api* kApis[2] = {nullptr,nullptr};
void InitializeApis() {
    kApis[0]=wh2_tiny_payload_cost_baseline_api();
    kApis[1]=wh2_tiny_payload_cost_candidate_api();
    for(const auto* api : kApis) Require(api && api->create && api->encode &&
        api->encoder_free && api->row && api->roundtrip,"complete POD API");
}
struct Handle {
    unsigned arm=0; void* value=nullptr;
    uint64_t* free_count=nullptr;
    Handle()=default; Handle(const Handle&)=delete; Handle& operator=(const Handle&)=delete;
    ~Handle() { if(value) { kApis[arm]->encoder_free(value); if(free_count) ++*free_count; } }
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
struct Snapshot {
    unsigned rep=0,lane=0,arm=0;
    uint64_t address=0;
    T::Roundtrip roundtrip={};
};
struct PreflightCounts {
    uint64_t encoder_creates=0,encoder_frees=0,encodes=0,decoder_creates=0,feeds=0,recovers=0,decoder_frees=0;
};
void SnapshotJson(std::ostream& out,const Snapshot& s) {
    out << "{\"rep\":" << s.rep << ",\"lane\":" << s.lane << ",\"arm\":" << s.arm
        << ",\"address\":" << s.address << ",\"roundtrip\":{\"create_code\":" << s.roundtrip.create_status
        << ",\"feed_codes\":[";
    for(unsigned i=0;i<s.roundtrip.feed_count;++i) { if(i) out << ','; out << s.roundtrip.feed_status[i]; }
    out << "],\"feed_count\":" << s.roundtrip.feed_count << ",\"recover_code\":" << s.roundtrip.recover_status
        << ",\"recover_bytes\":" << s.roundtrip.written << "}}";
}
struct Fixture {
    unsigned width_index,B,reps; size_t M,stride;
    const std::vector<uint8_t>& lookup;
    Buffer source; std::array<Buffer,2> outputs;
    std::array<Handle,48> handles;
    std::vector<uint8_t> packets,packet_reference;
    std::array<uint8_t,72> rows={};
    std::vector<Snapshot> snapshots;
    std::string source_hash,packet_hash;
    Fixture(unsigned index,unsigned block,unsigned repeat,const std::vector<uint8_t>& table)
        :width_index(index),B(block),reps(repeat),M(size_t(6)*block),stride(block+128u),lookup(table) {
        Require((reps==12 && width_index<3 && B==(width_index==0?2u:width_index==1?64u:1280u)) ||
                (reps==1 && (B==2 || B==65)),"fixed fixture domain");
        Require(lookup.size()==kLookupBytes,"fixture lookup extent");
        source.Create(M+128u);
        for(size_t i=0;i<M;++i) source.data[64u+i]=uint8_t(37u*i+i/11u);
        source_hash=Hash(Source(),M);
        for(auto& output : outputs) output.Create(12u*stride);
        std::copy(lookup.begin(),lookup.begin()+72,rows.begin());
        packets.assign(12u*B,0u);
        for(unsigned id=0;id<12;++id) for(unsigned byte=0;byte<B;++byte)
            for(unsigned j=0;j<6;++j)
                packets[size_t(id)*B+byte]^=Multiply(rows[6u*id+j],Source()[size_t(j)*B+byte]);
        packet_reference=packets; packet_hash=Hash(packets.data(),packets.size());
        snapshots.reserve(reps*4u);
    }
    const uint8_t* Source() const { return source.data+64u; }
    uint8_t* Output(unsigned lane,unsigned slot=0) { return outputs[lane].data+size_t(slot)*stride+64u; }
    Handle& At(unsigned rep,unsigned lane,unsigned arm) { return handles[(rep*2u+lane)*2u+arm]; }
    void Prepare() { for(auto& output : outputs) std::memset(output.data,0xa5,output.size); }
    void CheckImmutable() const {
        for(size_t i=0;i<source.size;++i) {
            const uint8_t expected=i<64u || i>=M+64u?0xa5u:uint8_t(37u*(i-64u)+(i-64u)/11u);
            Require(source.data[i]==expected,"immutable source/guard changed");
        }
        Require(packets==packet_reference,"immutable packet corpus changed");
        Require(std::equal(rows.begin(),rows.end(),lookup.begin()),"coefficient corpus changed");
    }
    void CheckOutput(unsigned lane,unsigned count,unsigned first_id,bool recovered=false) const {
        const Buffer& output=outputs[lane];
        for(size_t i=0;i<output.size;++i) {
            uint8_t expected=0xa5u;
            if(recovered) {
                if(i>=64u && i<64u+M) expected=Source()[i-64u];
            } else {
                const size_t slot=i/stride,offset=i%stride;
                if(slot<count && offset>=64u && offset<64u+B)
                    expected=packets[size_t(first_id+slot)*B+offset-64u];
            }
            Require(output.data[i]==expected,"packet/recovery output or guard changed");
        }
    }
    void Preflight(PreflightCounts& counts) {
        for(unsigned rep=0;rep<reps;++rep) for(unsigned step=0;step<2;++step) for(unsigned lane=0;lane<2;++lane) {
            const unsigned arm=(rep+width_index+step)%2u;
            const T::Api& api=*kApis[arm]; Handle& h=At(rep,lane,arm); h.arm=arm;
            h.free_count=&counts.encoder_frees;
            const int created=api.create(lookup.data(),lookup.size(),Source(),M,B,&h.value);
            ++counts.encoder_creates; Require(created==0 && h.value,"preflight encoder create");
            Snapshot snapshot; snapshot.rep=rep; snapshot.lane=lane; snapshot.arm=arm; snapshot.address=Address(h.value);
            Prepare();
            for(unsigned id=0;id<12;++id) {
                std::array<uint8_t,8> row={{0xa5,0,0,0,0,0,0,0xa5}};
                Require(api.row(lookup.data(),lookup.size(),id,row.data()+1u)==0 &&
                        row.front()==0xa5 && row.back()==0xa5 &&
                        std::equal(row.begin()+1,row.begin()+7,rows.begin()+6u*id),"coefficient byte oracle");
                const T::Result result=api.encode(h.value,id,Output(lane,id),B);
                ++counts.encodes;
                Require(result.status==0 && result.required==B && result.written==B,"preflight encode result");
                Require(std::memcmp(Output(lane,id),packets.data()+size_t(id)*B,B)==0,"independent polynomial payload oracle");
                if(id<6) Require(std::memcmp(Output(lane,id),Source()+size_t(id)*B,B)==0,"systematic identity");
            }
            CheckOutput(lane,12u,0u); CheckOutput(1u-lane,0u,0u); CheckImmutable();
            Prepare();
            snapshot.roundtrip=api.roundtrip(lookup.data(),lookup.size(),Source(),M,B,
                packets.data()+6u*B,6u*B,Output(lane),M);
            const T::Roundtrip& d=snapshot.roundtrip;
            ++counts.decoder_creates; counts.feeds+=d.feed_count;
            counts.recovers+=d.recover_status==-1?0u:1u;
            counts.decoder_frees+=d.create_status==0?1u:0u;
            Require(d.create_status==0 && d.feed_count==6u && d.recover_status==0 &&
                    d.required==M && d.written==M && d.recovered,"preflight fresh decoder result");
            for(unsigned i=0;i<6;++i) Require(d.feed_status[i]==(i==5u?0:1),"preflight exact first-success sequence");
            CheckOutput(lane,0u,0u,true); CheckOutput(1u-lane,0u,0u); CheckImmutable();
            snapshots.push_back(snapshot);
        }
        Require(snapshots.size()==size_t(reps)*4u,"preflight allocation roster");
    }
};

struct Coordinate { unsigned index,rep,order,width,metric,comparison,position,arm; uint64_t q; };
std::vector<Coordinate> Roster() {
    std::vector<Coordinate> result; result.reserve(4320u);
    for(unsigned r=0;r<12;++r) for(unsigned s=0;s<2;++s) for(unsigned ws=0;ws<3;++ws)
        for(unsigned ms=0;ms<2;++ms) for(unsigned cs=0;cs<3;++cs) for(unsigned p=0;p<10;++p) {
            const unsigned order=(r+s)%2u,width=(r+s+ws)%3u,metric=(r+s+ws+ms)%2u;
            const unsigned comparison=(2u*r+s+ws+ms+cs)%3u,j=p<2u?r%4u:(p-2u)/2u;
            Coordinate c={unsigned(result.size()),r,order,width,metric,comparison,p,
                kPairs[comparison][kSides[p]^order],(uint64_t(2u*(r+12u*j)+1u)*1000000u)/96u};
            result.push_back(c);
        }
    return result;
}

struct WorkResult {
    uint64_t creates=0,encodes=0,frees=0,handle=kMissing;
    int64_t create_code=kAbsent;
    std::array<int64_t,12> codes;
    std::array<uint64_t,12> required,written;
    std::array<uint64_t,64> addresses={{0}};
    unsigned count=0,address_count=0;
    uint64_t address_min=kMissing,address_max=kMissing;
    std::string address_hash;
    bool complete=false;
    WorkResult() { codes.fill(kAbsent); required.fill(kMissing); written.fill(kMissing); }
};
void FreeInside(const T::Api& api,Handle& temporary,WorkResult& result) {
    if(temporary.value) {
        void* handle=temporary.value; temporary.value=nullptr;
        api.encoder_free(handle); ++result.frees;
    }
}
WH2_COST_NOINLINE void RunWork(Fixture& fixture,const Coordinate& c,Handle& temporary,WorkResult& result) noexcept {
    const T::Api& api=*kApis[c.arm]; temporary.arm=c.arm;
    result.count=c.metric==0?6u:12u;
    if(c.metric==0) result.handle=Address(fixture.At(c.rep,c.order,c.arm).value);
    try {
        for(unsigned cycle=0;cycle<64;++cycle) {
            void* handle=fixture.At(c.rep,c.order,c.arm).value;
            if(c.metric==1) {
                result.create_code=api.create(fixture.lookup.data(),fixture.lookup.size(),fixture.Source(),fixture.M,fixture.B,&temporary.value);
                ++result.creates; handle=temporary.value;
                result.addresses[result.address_count++]=Address(handle);
                if(result.create_code!=0 || !handle) { FreeInside(api,temporary,result); return; }
            }
            if(!handle) return;
            for(unsigned j=0;j<result.count;++j) {
                const unsigned id=c.metric==0?6u+j:j;
                const T::Result value=api.encode(handle,id,fixture.Output(c.order,j),fixture.B);
                ++result.encodes;
                result.codes[j]=value.status; result.required[j]=value.required; result.written[j]=value.written;
                if(value.status!=0 || value.required!=fixture.B || value.written!=fixture.B) {
                    FreeInside(api,temporary,result); return;
                }
            }
            if(c.metric==1) FreeInside(api,temporary,result);
        }
        result.complete=true;
    } catch(...) {
        // Retain the post-WORK clock edges. Any still-owned temporary is
        // cleaned up outside this failing bracket by Handle's destructor.
        result.complete=false;
    }
}
void AddressSummary(WorkResult& result) {
    std::array<uint8_t,512> bytes={{0}};
    for(unsigned i=0;i<result.address_count;++i) {
        const uint64_t address=result.addresses[i];
        result.address_min=result.address_min==kMissing?address:std::min(result.address_min,address);
        result.address_max=result.address_max==kMissing?address:std::max(result.address_max,address);
        for(unsigned j=0;j<8;++j) bytes[8u*i+j]=uint8_t(address>>(8u*j));
    }
    result.address_hash=Hash(bytes.data(),8u*result.address_count);
}
void WorkJson(std::ostream& out,const WorkResult& w) {
    out << "{\"create_calls\":" << w.creates << ",\"encode_calls\":" << w.encodes << ",\"free_calls\":" << w.frees
        << ",\"create_code\":"; Signed(out,w.create_code);
    out << ",\"handle\":"; Number(out,w.handle);
    out << ",\"codes\":[";
    for(unsigned i=0;i<w.count;++i) { if(i) out << ','; Signed(out,w.codes[i]); }
    out << "],\"required\":[";
    for(unsigned i=0;i<w.count;++i) { if(i) out << ','; Number(out,w.required[i]); }
    out << "],\"written\":[";
    for(unsigned i=0;i<w.count;++i) { if(i) out << ','; Number(out,w.written[i]); }
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
void CheckWork(Fixture& fixture,Record& record,const Handle& temporary) {
    const auto& c=record.c; const auto& w=record.work;
    Require(w.complete && !temporary.value,"incomplete work or leaked lifecycle handle");
    Require(w.count==(c.metric==0?6u:12u) && w.encodes==(c.metric==0?384u:768u),"encode accounting");
    Require(w.creates==(c.metric==0?0u:64u) && w.frees==w.creates && w.address_count==w.creates,"lifecycle accounting");
    if(c.metric==0) Require(w.create_code==kAbsent && w.handle==Address(fixture.At(c.rep,c.order,c.arm).value) &&
                           w.handle>0,"persistent handle identity");
    else Require(w.create_code==0 && w.handle==kMissing,"fresh lifecycle status");
    for(unsigned i=0;i<w.count;++i) Require(w.codes[i]==0 && w.required[i]==fixture.B && w.written[i]==fixture.B,"packet result aggregate");
    for(unsigned i=0;i<w.address_count;++i) Require(w.addresses[i]>0,"null ephemeral handle");
    fixture.CheckOutput(c.order,w.count,c.metric==0?6u:0u);
    fixture.CheckOutput(1u-c.order,0u,0u); fixture.CheckImmutable(); record.checked=true;
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
        << ',' << Address(f.outputs[1].data) << "],\"output_stride\":" << f.stride
        << ",\"packets_hex\":" << Quote(Hex(f.packets.data(),f.packets.size()))
        << ",\"packets_sha256\":" << Quote(f.packet_hash) << ",\"rows_hex\":" << Quote(Hex(f.rows.data(),f.rows.size()))
        << ",\"handles\":[";
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
    uint64_t records=0,callbacks=0,checked=0,creates=0,encodes=0,frees=0,work=0,wait_wall=0,wait_cpu=0;
    void Add(const Record& r) {
        ++records; callbacks+=r.called?1u:0u; checked+=r.checked?1u:0u;
        creates+=r.work.creates; encodes+=r.work.encodes; frees+=r.work.frees;
        const Observation& o=r.observation;
        if(o.m1!=kMissing && o.m2!=kMissing && o.m2>=o.m1) work+=o.m2-o.m1;
        if(r.wait[0]!=kMissing && r.wait[2]!=kMissing && r.wait[2]>=r.wait[0]) wait_wall+=r.wait[2]-r.wait[0];
        if(r.wait[1]!=kMissing && r.wait[3]!=kMissing && r.wait[3]>=r.wait[1]) wait_cpu+=r.wait[3]-r.wait[1];
    }
};
void FooterJson(std::ostream& out,const Accounting& a,const PreflightCounts& p,const std::string& failure,
                uint64_t schedule_now,uint64_t t0,uint64_t end,uint64_t cpu_end,
                const std::string& identity,const std::string& runtime) {
    out << "{\"type\":\"footer\",\"protocol\":" << Quote(kProtocol)
        << ",\"outcome\":" << Quote(failure.empty()?"COMPLETE":"INVALID")
        << ",\"failure\":" << (failure.empty()?"null":Quote(failure)) << ",\"schedule_now_ns\":";
    Number(out,schedule_now); out << ",\"t0\":"; Number(out,t0);
    out << ",\"records\":" << a.records << ",\"callbacks\":" << a.callbacks << ",\"checked_callbacks\":" << a.checked
        << ",\"create_calls\":" << a.creates << ",\"encode_calls\":" << a.encodes << ",\"free_calls\":" << a.frees
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
    const rlimit cpu={14,14},memory={512u*1024u*1024u,512u*1024u*1024u},core={0,0};
    Require(setrlimit(RLIMIT_CPU,&cpu)==0 && setrlimit(RLIMIT_AS,&memory)==0 &&
            setrlimit(RLIMIT_CORE,&core)==0,"worker resource limits");
}
int Worker(const std::string& claim) {
    Reader reader; Accounting accounting; PreflightCounts preflight_counts;
    uint64_t start=kMissing,cpu_start=kMissing,schedule_now=kMissing,t0=kMissing,end=kMissing,cpu_end=kMissing;
    std::string failure,identity_before,runtime_before,identity_after,runtime_after;
    try {
        Require(WH2_TINY_COST_SANITIZERS==0,"sanitizer build cannot run scientific worker");
        Authenticate(claim); Limits(); start=reader.Now(); cpu_start=reader.Thread(); Pin();
        Require(gf256_init()==0,"shared GF initialization"); InitializeApis();
        identity_before=Identity(); runtime_before=Runtime();
        const std::array<uint64_t,2> resolutions={{reader.Resolution(CLOCK_MONOTONIC),reader.Resolution(CLOCK_THREAD_CPUTIME_ID)}};
        Require(resolutions[0]==1u && resolutions[1]==1u,"frozen clock resolutions");
        const std::string packed=ReadPinned(WH2_TINY_COST_LOOKUP_PATH,kLookupBytes,kLookupHash);
        Require(packed.size()==kLookupBytes,"scientific lookup exact extent");
        const std::vector<uint8_t> lookup(packed.begin(),packed.end());
        const std::vector<uint8_t> lookup_reference=lookup;
        {
            Fixture f0(0u,2u,12u,lookup),f1(1u,64u,12u,lookup),f2(2u,1280u,12u,lookup);
            const std::array<Fixture*,3> fixtures={{&f0,&f1,&f2}};
            for(Fixture* f : fixtures) {
                f->Preflight(preflight_counts);
                Budget(reader,start,cpu_start,0u);
            }
            Require(lookup==lookup_reference,"preflight lookup changed");
            const auto roster=Roster(); Require(roster.size()==4320u,"complete roster size");
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
                Fixture& fixture=*fixtures[c.width]; Record record(c); Handle temporary;
                try {
                    record.target=Target(t0,c); fixture.Prepare();
                    Budget(reader,start,cpu_start,accounting.work);
                    record.ready=reader.Now(); CheckReady(record.ready,record.target,previous);
                    Wait(reader,record.target,start,record.wait);
                    error=Capture(reader,[&] { RunWork(fixture,c,temporary,record.work); record.called=true; },record.observation);
                    Require(!error,error?error:"callback capture"); Validate(record.observation,&previous); CheckStart(record);
                    for(const Fixture* f : fixtures) f->CheckImmutable();
                    Require(lookup==lookup_reference,"immutable global lookup changed");
                    CheckWork(fixture,record,temporary); AddressSummary(record.work);
                    Budget(reader,start,cpu_start,accounting.work+record.observation.m2-record.observation.m1);
                    Require(OnCpu(),"worker migrated/affinity changed"); previous=record.observation;
                } catch(...) {
                    AddressSummary(record.work); accounting.Add(record); RecordJson(std::cout,record);
                    std::cout << '\n' << std::flush; throw;
                }
                accounting.Add(record); RecordJson(std::cout,record); std::cout << '\n' << std::flush;
                Require(bool(std::cout),"record stream write");
            }
            Require(accounting.records==4320u && accounting.callbacks==4320u && accounting.checked==4320u &&
                    accounting.creates==138240u && accounting.encodes==2488320u && accounting.frees==138240u,
                    "complete callback/API accounting");
            for(const Fixture* f : fixtures) f->CheckImmutable();
            Require(lookup==lookup_reference && OnCpu(),"final lookup/affinity");
            // Teardown of all144 persistent encoders is in worker elapsed,
            // outside the persistent repair component measurements.
        }
        Require(preflight_counts.encoder_creates==144u && preflight_counts.encoder_frees==144u &&
                preflight_counts.encodes==1728u && preflight_counts.decoder_creates==144u &&
                preflight_counts.feeds==864u && preflight_counts.recovers==144u && preflight_counts.decoder_frees==144u,
                "complete preflight accounting");
        identity_after=Identity(); runtime_after=Runtime();
        Require(identity_after==identity_before && runtime_after==runtime_before,"runtime/physical identity changed");
        Budget(reader,start,cpu_start,accounting.work);
    } catch(const std::exception& error) { failure=error.what(); }
      catch(...) { failure="unexpected worker exception"; }
    if(start!=kMissing) {
        if(!reader.Mono(end) || !reader.Cpu(cpu_end)) {
            if(failure.empty()) failure="footer clock capture";
        } else if(end<start || end-start>=kWall || cpu_end<cpu_start || cpu_end-cpu_start>kCpu) {
            if(failure.empty()) failure="final cleanup resource cap";
        }
    }
    FooterJson(std::cout,accounting,preflight_counts,failure,schedule_now,t0,end,cpu_end,identity_after,runtime_after);
    std::cout << '\n' << std::flush; return failure.empty() && bool(std::cout)?0:1;
}

struct FakeReader {
    unsigned calls=0,fail=UINT32_MAX; uint64_t mono=1000,cpu=1000000;
    bool Mono(uint64_t& value) { if(calls++==fail) return false; value=mono; mono+=100; return true; }
    bool Cpu(uint64_t& value) { if(calls++==fail) return false; value=cpu; cpu+=200; return true; }
    bool Usage(Counters& value) { if(calls++==fail) return false; value={{0,0,0,0}}; return true; }
};
void NeutralClocks() {
    uint64_t last_mono=1000u,last_cpu=1000000u;
    Require(!OrderedClock(999u,last_mono) && last_mono==1000u &&
            !OrderedClock(999999u,last_cpu) && last_cpu==1000000u,"neutral clock regression above worker start");
    Require(OrderedClock(1000u,last_mono) && OrderedClock(1000001u,last_cpu),"neutral independent clock domains");
    FakeReader reader; Observation o; unsigned work=0;
    Require(Capture(reader,[&]{++work;},o)==nullptr && work==1u,"neutral complete capture");
    Validate(o,nullptr);
    for(unsigned failure=0;failure<8;++failure) {
        FakeReader bad; bad.fail=failure; Observation partial; unsigned completed=0;
        Require(Capture(bad,[&]{++completed;},partial)!=nullptr,"neutral failed clock capture");
        Require(completed==(failure>=4u?1u:0u),"neutral partial callback boundary");
    }
    Observation later=o; later.m0=1u;
    bool rejected=false; try { Validate(later,&o); } catch(const std::runtime_error&) { rejected=true; }
    Require(rejected,"neutral reverse clock rejection");
    const auto roster=Roster(); Require(roster.size()==4320u,"neutral roster size");
    unsigned seen[3][2][3][2][12]={},phases[3][2][3][2][48]={};
    uint64_t minimum=kMissing,maximum=0;
    for(size_t i=0;i<roster.size();++i) {
        const auto& c=roster[i];
        Require(c.index==i && c.width<3 && c.arm<2 && c.metric<2 && c.comparison<3 && c.order<2 && c.rep<12,
                "neutral coordinate bounds");
        if(c.position==0) ++seen[c.width][c.metric][c.comparison][c.order][c.rep];
        if(c.position>=2 && c.position%2u==0) {
            const auto& next=roster[i+1];
            Require(next.q==c.q && next.position==c.position+1u &&
                    (kSides[next.position]^next.order)!=(kSides[c.position]^c.order),"neutral paired phases");
            ++phases[c.width][c.metric][c.comparison][c.order][c.rep+12u*((c.position-2u)/2u)];
        }
        if(i) { const uint64_t gap=Target(10000000u,c)-Target(10000000u,roster[i-1]);
            minimum=std::min(minimum,gap); maximum=std::max(maximum,gap); }
    }
    for(unsigned w=0;w<3;++w) for(unsigned m=0;m<2;++m) for(unsigned c=0;c<3;++c) for(unsigned order=0;order<2;++order) {
        for(unsigned rep=0;rep<12;++rep) Require(seen[w][m][c][order][rep]==1u,"neutral exact panel coverage");
        for(unsigned phase=0;phase<48;++phase) Require(phases[w][m][c][order][phase]==1u,"neutral all48 phases");
    }
    Require(minimum==1250000u && maximum==2250000u,"neutral nonuniform slot gaps");
    Record r(roster[0]); r.target=10000000u; r.observation.m1=r.target+5000u; CheckStart(r);
    r.observation.m1++;
    rejected=false; try { CheckStart(r); } catch(const std::runtime_error&) { rejected=true; }
    Require(rejected,"neutral late start rejection");
    Observation previous; previous.m3=1;
    CheckReady(r.target-100000u,r.target,previous);
    rejected=false; try { CheckReady(r.target-99999u,r.target,previous); } catch(const std::runtime_error&) { rejected=true; }
    Require(rejected,"neutral preparation deadline rejection");
    Require(Prelude(kPreludeSeed,1u)==UINT64_C(0xdc1b77ae0bf34dad) &&
            Prelude(kPreludeSeed,2u)==UINT64_C(0x64f0eeb9026e6076),"neutral short prelude literals");
}
void NeutralBridges(Fixture& fixture) {
    for(unsigned arm=0;arm<2;++arm) {
        const auto& api=*kApis[arm]; Handle& h=fixture.At(0,0,arm);
        void* occupied=h.value;
        Require(api.create(fixture.lookup.data(),fixture.lookup.size(),fixture.Source(),fixture.M,fixture.B,&occupied)==2 &&
                occupied==h.value,"neutral nonempty create output");
        Require(api.create(fixture.lookup.data(),fixture.lookup.size(),fixture.Source(),fixture.M,fixture.B,nullptr)==2,
                "neutral null create output");
        Require(api.encode(nullptr,6u,fixture.Output(0),fixture.B).status==2,"neutral null encode handle");
        fixture.Prepare();
        const T::Result short_result=api.encode(h.value,6u,fixture.Output(0),fixture.B-1u);
        Require(short_result.status==3 && short_result.required==fixture.B && short_result.written==0u,"neutral short encode");
        fixture.CheckOutput(0,0,0); fixture.CheckOutput(1,0,0);
        const T::Roundtrip short_decode=api.roundtrip(fixture.lookup.data(),fixture.lookup.size(),fixture.Source(),fixture.M,
            fixture.B,fixture.packets.data()+6u*fixture.B,6u*fixture.B,fixture.Output(0),fixture.M-1u);
        Require(short_decode.create_status==2 && short_decode.feed_count==0 && !short_decode.recovered,"neutral short roundtrip");
        const T::Roundtrip aliased=api.roundtrip(fixture.lookup.data(),fixture.lookup.size(),fixture.Source(),fixture.M,
            fixture.B,fixture.packets.data()+6u*fixture.B,6u*fixture.B,const_cast<uint8_t*>(fixture.Source()),fixture.M);
        Require(aliased.create_status==2 && aliased.feed_count==0 && !aliased.recovered,"neutral aliased roundtrip");
        fixture.CheckImmutable(); api.encoder_free(nullptr);
    }
}
void NeutralMetrics(Fixture& fixture) {
    for(unsigned metric=0;metric<2;++metric) for(unsigned arm=0;arm<2;++arm) {
        const Coordinate c={0,0,arm,fixture.width_index,metric,0,0,arm,0};
        Record record(c); Handle temporary; fixture.Prepare(); FakeReader reader;
        Require(Capture(reader,[&] { RunWork(fixture,c,temporary,record.work); record.called=true; },record.observation)==nullptr,
                "neutral metric synthetic capture");
        Validate(record.observation,nullptr); CheckWork(fixture,record,temporary); AddressSummary(record.work);
        Require(record.checked,"neutral actual API metric bytes");
        fixture.outputs[1u-c.order].data[0]^=1u;
        bool rejected=false; try { fixture.CheckOutput(1u-c.order,0,0); } catch(const std::runtime_error&) { rejected=true; }
        Require(rejected,"neutral cross-lane guard"); fixture.outputs[1u-c.order].data[0]^=1u;
        std::ostringstream out; RecordJson(out,record); Require(out.str().size()<2048u,"neutral compact record");
    }
}
int Neutral() {
    Require(gf256_init()==0,"neutral shared GF initialization"); InitializeApis(); NeutralClocks();
    const std::vector<uint8_t> lookup=NeutralLookup(),saved=lookup;
    for(unsigned B : {2u,65u}) {
        PreflightCounts counts;
        {
            Fixture fixture(0u,B,1u,lookup); fixture.Preflight(counts); NeutralBridges(fixture); NeutralMetrics(fixture);
        }
        Require(counts.encoder_creates==4u && counts.encoder_frees==4u && counts.encodes==48u &&
                counts.decoder_creates==4u && counts.feeds==24u && counts.recovers==4u && counts.decoder_frees==4u,
                "neutral preflight/free accounting");
        Require(lookup==saved,"neutral lookup unchanged");
    }
    std::cout << "PASS neutral unrelated pair B2/65, payload/decoder/bridge and fake-clock/roster checks\n";
    return 0;
}
} // namespace cost

int main(int argc,char** argv) {
    try {
        if(argc==2 && std::strcmp(argv[1],"--neutral")==0) return cost::Neutral();
        if(argc==3 && std::strcmp(argv[1],"--worker")==0) return cost::Worker(argv[2]);
        std::cerr << "usage: --neutral | --worker CLAIM_SHA256 (no default execution)\n"; return 2;
    } catch(const std::exception& error) { std::cerr << "FAIL: " << error.what() << '\n'; return 1; }
}
