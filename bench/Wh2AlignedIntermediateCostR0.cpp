// .68 frozen, one-shot full-cost screen.  No production source is modified.
#include "Wh2AlignedIntermediateCostR0Bridge.h"
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
    defined(WH_ALIGN64) || defined(WH_HUGEPAGE)
#error "Scientific and neutral cost-driver builds use unmodified runtime policy"
#endif
#ifndef WH2_ALIGNMENT_COST_SANITIZERS
#error "Builder must identify sanitizer versus native execution"
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
static const char kProtocol[] = "wirehair.wh2.aligned-intermediate-cost-r0";
static const char kClaim[] = "/var/tmp/wh2-aligned-intermediate-cost-r0/CLAIM.json";
static const uint64_t kMissing = UINT64_MAX;
static const int64_t kAbsent = INT64_MIN;
static const uint64_t kWall = UINT64_C(15000000000);
static const uint64_t kCpu = UINT64_C(14000000000);
static const uint64_t kWork = UINT64_C(1000000000);
static const uint64_t kPreludeSeed = UINT64_C(0x9e3779b97f4a7c15);
static const uint64_t kPreludeFinal = UINT64_C(0x43935dad1647741b);
static const unsigned kPreludeIterations = 1u << 20;
static const unsigned kSides[10] = {0,1,0,1,1,0,1,0,0,1};
static const unsigned kPairs[5][2] = {{0,0},{1,1},{2,2},{1,2},{0,2}};

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

// Common adapters preserve the actual public APIs, including their validation.
// No V2 C++ type or owning object crosses between the isolated namespaces.
struct Api {
    int64_t (*create_encoder)(const void*,uint64_t,uint32_t,uint8_t*,uint32_t*,void**);
    int64_t (*encode)(void*,uint32_t,void*,uint32_t,uint32_t*);
    int64_t (*create_decoder)(uint64_t,uint32_t,const uint8_t*,void**);
    int64_t (*decode)(void*,uint32_t,const void*,uint32_t);
    int64_t (*recover)(void*,void*,uint64_t,uint64_t*);
    void (*free)(void*);
};
int64_t WCreate(const void* source,uint64_t M,uint32_t B,uint8_t*,uint32_t* length,void** out) {
    WirehairCodec handle=nullptr;
    const WirehairResult code=wirehair_encoder_create_ex(nullptr,source,M,B,&handle);
    *out=handle; *length=0u; return code;
}
int64_t WEncode(void* h,uint32_t id,void* output,uint32_t capacity,uint32_t* bytes) {
    return wirehair_encode(static_cast<WirehairCodec>(h),id,output,capacity,bytes);
}
int64_t WDecoder(uint64_t M,uint32_t B,const uint8_t*,void** out) {
    WirehairCodec handle=nullptr;
    const WirehairResult code=wirehair_decoder_create_ex(nullptr,M,B,&handle);
    *out=handle; return code;
}
int64_t WDecode(void* h,uint32_t id,const void* data,uint32_t bytes) {
    return wirehair_decode(static_cast<WirehairCodec>(h),id,data,bytes);
}
int64_t WRecover(void* h,void* data,uint64_t capacity,uint64_t* bytes) {
    const WirehairResult code=wirehair_recover(static_cast<WirehairCodec>(h),data,capacity);
    if(code==Wirehair_Success) *bytes=capacity;
    return code;
}
void WFree(void* h) { wirehair_free(static_cast<WirehairCodec>(h)); }
#define WH2_DEFINE_API(prefix, api) \
int64_t prefix##Create(const void* source,uint64_t M,uint32_t B,uint8_t* profile,uint32_t* length,void** out) { \
    WirehairV2Codec handle=nullptr; WirehairV2EncoderOptions options={}; \
    options.struct_bytes=sizeof(options); options.options_version=WIREHAIR_V2_ENCODER_OPTIONS_VERSION; \
    options.source_policy=WirehairV2EncoderSource_BorrowedImmutable; \
    const WirehairV2Result code=api##_encoder_create_with_options(source,M,B,&options,profile,32u,length,&handle); \
    *out=handle; return code; } \
int64_t prefix##Encode(void* h,uint32_t id,void* output,uint32_t capacity,uint32_t* bytes) { \
    return api##_encode(static_cast<WirehairV2Codec>(h),id,output,capacity,bytes); } \
int64_t prefix##Decoder(uint64_t,uint32_t,const uint8_t* profile,void** out) { \
    WirehairV2Codec handle=nullptr; const WirehairV2Result code=api##_decoder_create(profile,32u,&handle); \
    *out=handle; return code; } \
int64_t prefix##Decode(void* h,uint32_t id,const void* data,uint32_t bytes) { \
    return api##_decode(static_cast<WirehairV2Codec>(h),id,data,bytes); } \
int64_t prefix##Recover(void* h,void* output,uint64_t capacity,uint64_t* bytes) { \
    return api##_recover(static_cast<WirehairV2Codec>(h),output,capacity,bytes); } \
void prefix##Free(void* h) { api##_free(static_cast<WirehairV2Codec>(h)); }
WH2_DEFINE_API(P,wirehair_v2)
WH2_DEFINE_API(A,wh2_aligned_r0)
#undef WH2_DEFINE_API
static const Api kApis[3] = {
    {WCreate,WEncode,WDecoder,WDecode,WRecover,WFree},
    {PCreate,PEncode,PDecoder,PDecode,PRecover,PFree},
    {ACreate,AEncode,ADecoder,ADecode,ARecover,AFree}
};
struct Handle {
    unsigned arm=0; void* value=nullptr; std::array<uint8_t,32> profile={{0}};
    Handle()=default;
    Handle(const Handle&)=delete; Handle& operator=(const Handle&)=delete;
    ~Handle() { if(value) kApis[arm].free(value); }
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
struct Snapshot {
    unsigned arm=0,rep=0,lane=0;
    uint64_t handle=0,source=0,intermediate=kMissing,bytes=kMissing,capacity=kMissing,policy=kMissing;
    std::string profile,source_hash,intermediate_hash;
};
void SnapshotJson(std::ostream& out,const Snapshot& s) {
    out << "{\"arm\":" << s.arm << ",\"rep\":" << s.rep << ",\"lane\":" << s.lane
        << ",\"handle\":" << s.handle << ",\"source\":" << s.source << ",\"intermediate\":";
    Number(out,s.intermediate); out << ",\"intermediate_bytes\":"; Number(out,s.bytes);
    out << ",\"intermediate_capacity\":"; Number(out,s.capacity);
    out << ",\"source_policy\":"; Number(out,s.policy);
    out << ",\"profile_hex\":" << (s.arm?Quote(s.profile):"null")
        << ",\"source_sha256\":" << Quote(s.source_hash)
        << ",\"intermediate_sha256\":" << (s.arm?Quote(s.intermediate_hash):"null") << '}';
}

struct Fixture {
    unsigned K,B,reps; size_t M,stride,packet_count;
    Buffer source; std::array<Buffer,2> outputs;
    std::array<Handle,72> handles;
    std::array<std::vector<uint8_t>,3> packets;
    std::array<std::vector<uint8_t>,3> packet_reference;
    std::array<std::array<uint8_t,32>,3> profiles;
    std::array<std::vector<uint8_t>,2> intermediate;
    std::array<std::vector<std::vector<uint32_t>>,2> columns;
    std::array<std::vector<int64_t>,3> decode_statuses;
    std::array<unsigned,3> first_success={{0,0,0}};
    std::vector<Snapshot> snapshots;
    std::string source_hash;
    std::array<std::string,2> intermediate_hash;
    explicit Fixture(unsigned count,unsigned block,unsigned repeat):K(count),B(block),reps(repeat),
        M(size_t(count)*block),stride(block+128u),packet_count(count+6u) {
        Require((K==6 && B==1280 && reps==12) ||
                ((K==7 || K==9) && (B==17 || B==65) && reps==1),"fixed fixture domain");
        source.Create(M+128u);
        for(size_t i=0;i<M;++i) source.data[64u+i]=uint8_t(37u*i+i/11u);
        source_hash=Hash(Source(),M);
        for(auto& output : outputs) output.Create(packet_count*stride);
        for(auto& profile : profiles) profile.fill(0u);
    }
    const uint8_t* Source() const { return source.data+64u; }
    uint8_t* Output(unsigned lane,unsigned id=0) { return outputs[lane].data+size_t(id)*stride+64u; }
    Handle& At(unsigned rep,unsigned lane,unsigned arm) { return handles[(rep*2u+lane)*3u+arm]; }
    void CheckSource() const {
        for(size_t i=0;i<source.size;++i) {
            const uint8_t expected=i<64u || i>=M+64u ? 0xa5u : uint8_t(37u*(i-64u)+(i-64u)/11u);
            Require(source.data[i]==expected,"immutable source or source guard");
        }
    }
    Snapshot Inspect(unsigned arm,unsigned rep,unsigned lane,void* handle,const uint8_t* profile) {
        Snapshot s; s.arm=arm; s.rep=rep; s.lane=lane; s.handle=Address(handle);
        s.source=Address(Source()); s.source_hash=source_hash;
        Require(handle!=nullptr,"missing inspection handle");
        if(!arm) return s;
        Wh2AlignedIntermediateCostSnapshot observed={};
        const auto inspect=arm==1?wh2_aligned_cost_p_inspect:wh2_aligned_cost_a_inspect;
        Require(inspect(static_cast<WirehairV2Codec>(handle),&observed)==WirehairV2_Success,"public encoder inspection");
        Require(observed.handle==handle && observed.source==Source() && observed.message_bytes==M &&
                observed.block_bytes==B && observed.source_count==K && observed.source_policy==2u &&
                observed.intermediate && observed.intermediate_bytes==uint64_t(K+observed.precode_count)*B &&
                observed.intermediate_capacity>=observed.intermediate_bytes &&
                std::memcmp(observed.serialized_profile,profile,32u)==0,"inspection dimensions/profile/ownership");
        if(K==6u) Require(observed.precode_count==30u && observed.intermediate_bytes==46080u &&
                         observed.intermediate_capacity==46080u,"frozen K6 retained storage extent");
        if(arm==2) Require(Address(observed.intermediate)%32u==0u,"candidate retained alignment");
        s.intermediate=Address(observed.intermediate); s.bytes=observed.intermediate_bytes;
        s.capacity=observed.intermediate_capacity; s.policy=observed.source_policy;
        s.profile=Hex(profile,32u);
        if(intermediate[arm-1].empty()) {
            intermediate[arm-1].assign(observed.intermediate,observed.intermediate+observed.intermediate_bytes);
            intermediate_hash[arm-1]=Hash(observed.intermediate,size_t(observed.intermediate_bytes));
        }
        Require(observed.intermediate_bytes==intermediate[arm-1].size() &&
                std::memcmp(observed.intermediate,intermediate[arm-1].data(),intermediate[arm-1].size())==0,
                "retained intermediate differs");
        s.intermediate_hash=intermediate_hash[arm-1]; return s;
    }
    void CheckOutput(unsigned lane,unsigned arm,unsigned count,unsigned first_id,bool recover=false) const {
        const Buffer& output=outputs[lane];
        for(size_t i=0;i<output.size;++i) {
            uint8_t expected=0xa5u;
            if(recover) {
                if(i>=64u && i<64u+M) expected=Source()[i-64u];
            } else {
                const size_t slot=i/stride,offset=i%stride;
                if(slot<count && offset>=64u && offset<64u+B)
                    expected=packets[arm][size_t(first_id+slot)*B+offset-64u];
            }
            Require(output.data[i]==expected,"packet/recovery bytes or output guards");
        }
    }
    void DecodePreflight(unsigned arm,const uint8_t* profile,std::vector<int64_t>& codes) {
        Handle decoder; decoder.arm=arm;
        Require(kApis[arm].create_decoder(M,B,profile,&decoder.value)==0 && decoder.value,"preflight decoder create");
        std::memset(outputs[0].data,0xa5,outputs[0].size); codes.clear();
        for(unsigned step=0;step<packet_count;++step) {
            const unsigned id=step<6u?K+step:step-6u;
            const int64_t code=kApis[arm].decode(decoder.value,id,packets[arm].data()+size_t(id)*B,B);
            codes.push_back(code);
            Require(code==0 || code==1,"preflight decoder status");
            if(code==0) break;
        }
        Require(!codes.empty() && codes.back()==0,"preflight decoder did not complete fixed sequence");
        uint64_t bytes=kMissing;
        Require(kApis[arm].recover(decoder.value,Output(0),M,&bytes)==0 && bytes==M,"preflight recovery result");
        CheckOutput(0,arm,0,0,true);
    }
    void Preflight() {
        for(unsigned r=0;r<reps;++r) for(unsigned i=0;i<3;++i) for(unsigned lane=0;lane<2;++lane) {
            const unsigned arm=(r+i)%3u;
            Handle& h=At(r,lane,arm); h.arm=arm;
            uint32_t bytes=UINT32_MAX;
            Require(kApis[arm].create_encoder(Source(),M,B,h.profile.data(),&bytes,&h.value)==0 && h.value,
                    "preflight borrowed encoder create");
            Require(bytes==(arm?32u:0u),"preflight profile length");
            if(packets[arm].empty()) profiles[arm]=h.profile;
            Require(h.profile==profiles[arm],"per-handle profile differs");
            snapshots.push_back(Inspect(arm,r,lane,h.value,h.profile.data()));
            std::vector<uint8_t> actual(packet_count*B);
            for(unsigned id=0;id<packet_count;++id) {
                std::memset(outputs[lane].data,0xa5,outputs[lane].size);
                bytes=UINT32_MAX;
                Require(kApis[arm].encode(h.value,id,Output(lane),B,&bytes)==0 && bytes==B,"preflight encode result");
                std::copy(Output(lane),Output(lane)+B,actual.begin()+size_t(id)*B);
                for(size_t j=0;j<outputs[lane].size;++j)
                    if(j<64u || j>=64u+B) Require(outputs[lane].data[j]==0xa5u,"preflight output guard");
                if(id<K) Require(std::memcmp(Output(lane),Source()+size_t(id)*B,B)==0,"systematic source bytes");
                if(arm) {
                    const auto scalar=arm==1?wh2_aligned_cost_p_scalar_packet:wh2_aligned_cost_a_scalar_packet;
                    const auto query=arm==1?wh2_aligned_cost_p_packet_columns:wh2_aligned_cost_a_packet_columns;
                    std::vector<uint8_t> oracle(B+128u,0xa5u); bytes=UINT32_MAX;
                    Require(scalar(static_cast<WirehairV2Codec>(h.value),id,oracle.data()+64u,B,&bytes)==WirehairV2_Success &&
                            bytes==B && std::memcmp(oracle.data()+64u,Output(lane),B)==0,"original scalar packet oracle");
                    for(size_t j=0;j<oracle.size();++j)
                        if(j<64u || j>=64u+B) Require(oracle[j]==0xa5u,"scalar oracle guard");
                    std::vector<uint32_t> row(K+128u,UINT32_MAX); uint32_t count=UINT32_MAX;
                    Require(query(static_cast<WirehairV2Codec>(h.value),id,row.data(),uint32_t(row.size()),&count)==WirehairV2_Success &&
                            count>0u && count<=row.size(),"canonical packet columns"); row.resize(count);
                    if(columns[arm-1].size()<packet_count) columns[arm-1].push_back(row);
                    else Require(columns[arm-1][id]==row,"per-handle columns differ");
                }
            }
            if(packets[arm].empty()) packets[arm]=actual;
            else Require(packets[arm]==actual,"per-handle packet bytes differ");
            std::vector<int64_t> codes; DecodePreflight(arm,h.profile.data(),codes);
            Require(packets[arm]==actual,"preflight decoder modified its packet input");
            if(decode_statuses[arm].empty()) decode_statuses[arm]=codes;
            else Require(decode_statuses[arm]==codes,"per-handle decode statuses differ");
            first_success[arm]=unsigned(codes.size()-1u); CheckSource();
        }
        Require(profiles[1]==profiles[2] && intermediate[0]==intermediate[1] && packets[1]==packets[2] &&
                columns[0]==columns[1] && decode_statuses[1]==decode_statuses[2],"P/A equation and byte identity");
        Require(snapshots.size()==size_t(reps)*6u,"precreated handle roster");
        packet_reference=packets;
    }
};

struct Coordinate { unsigned index,rep,order,metric,comparison,position,arm; uint64_t q; };
std::vector<Coordinate> Roster() {
    std::vector<Coordinate> result; result.reserve(4800u);
    for(unsigned r=0;r<12;++r) for(unsigned s=0;s<2;++s)
        for(unsigned ms=0;ms<4;++ms) for(unsigned cs=0;cs<5;++cs) for(unsigned p=0;p<10;++p) {
            const unsigned order=(r+s)%2u,metric=(r+s+ms)%4u,comparison=(2u*r+s+cs)%5u;
            const unsigned j=p<2u?r%4u:(p-2u)/2u;
            Coordinate c={unsigned(result.size()),r,order,metric,comparison,p,kPairs[comparison][kSides[p]^order],
                          (uint64_t(2u*(r+12u*j)+1u)*1000000u)/96u};
            result.push_back(c);
        }
    return result;
}

struct WorkResult {
    uint64_t calls=0,first_success=kMissing,bytes=kMissing,handle=kMissing;
    int64_t create_code=kAbsent,recover_code=kAbsent;
    std::array<int64_t,16> codes;
    std::array<uint64_t,16> lengths;
    unsigned count=0;
    bool freed=false,complete=false;
    std::array<uint8_t,32> profile={{0}};
    uint32_t profile_bytes=UINT32_MAX;
    WorkResult() { codes.fill(kAbsent); lengths.fill(kMissing); }
};
void FreeInside(const Api& api,Handle& temporary,WorkResult& result) {
    if(temporary.value) {
        void* value=temporary.value; temporary.value=nullptr;
        api.free(value); ++result.calls; result.freed=true;
    }
}
WH2_COST_NOINLINE void RunWork(Fixture& fixture,const Coordinate& c,
                              Handle& temporary,WorkResult& result) noexcept {
    const Api& api=kApis[c.arm]; temporary.arm=c.arm;
    try {
        void* handle=nullptr;
        if(c.metric==1) handle=fixture.At(c.rep,c.order,c.arm).value;
        else if(c.metric==3) {
            result.create_code=api.create_decoder(fixture.M,fixture.B,fixture.profiles[c.arm].data(),&temporary.value);
            ++result.calls; handle=temporary.value;
        } else {
            result.create_code=api.create_encoder(fixture.Source(),fixture.M,fixture.B,
                result.profile.data(),&result.profile_bytes,&temporary.value);
            ++result.calls; handle=temporary.value;
        }
        result.handle=Address(handle);
        if(!handle || (c.metric!=1 && result.create_code!=0)) {
            if(c.metric>=2) FreeInside(api,temporary,result);
            return;
        }
        if(c.metric==0) { result.complete=true; return; }
        if(c.metric==3) {
            for(unsigned step=0;step<fixture.packet_count;++step) {
                const unsigned id=step<6u?fixture.K+step:step-6u;
                result.codes[step]=api.decode(handle,id,fixture.packets[c.arm].data()+size_t(id)*fixture.B,fixture.B);
                result.lengths[step]=fixture.B; ++result.calls; result.count=step+1u;
                if(result.codes[step]==0) { result.first_success=step; break; }
                if(result.codes[step]!=1) { FreeInside(api,temporary,result); return; }
            }
            if(result.first_success==kMissing) { FreeInside(api,temporary,result); return; }
            result.recover_code=api.recover(handle,fixture.Output(c.order),fixture.M,&result.bytes);
            ++result.calls;
            FreeInside(api,temporary,result);
            result.complete=result.recover_code==0 && result.bytes==fixture.M;
            return;
        }
        const unsigned packets=c.metric==1?6u:unsigned(fixture.packet_count);
        const unsigned cycles=c.metric==1?64u:1u;
        result.count=packets;
        for(unsigned cycle=0;cycle<cycles;++cycle) for(unsigned j=0;j<packets;++j) {
            const unsigned id=c.metric==1?fixture.K+j:j;
            uint32_t bytes=UINT32_MAX;
            const int64_t code=api.encode(handle,id,fixture.Output(c.order,j),fixture.B,&bytes);
            ++result.calls;
            if(cycle==0 || code!=0) result.codes[j]=code;
            if(cycle==0 || bytes!=fixture.B) result.lengths[j]=bytes;
            if(code!=0 || bytes!=fixture.B) {
                if(c.metric==2) FreeInside(api,temporary,result);
                return;
            }
        }
        if(c.metric==2) FreeInside(api,temporary,result);
        result.complete=true;
    } catch(...) {
        // Preserve post-clock edges even if an unexpected C++ exception escapes
        // a public C entry.  Any owned temporary is released outside the timer.
        result.complete=false;
    }
}
void WorkJson(std::ostream& out,const WorkResult& w) {
    out << "{\"calls\":" << w.calls << ",\"create_code\":"; Signed(out,w.create_code);
    out << ",\"recover_code\":"; Signed(out,w.recover_code);
    out << ",\"codes\":[";
    for(unsigned i=0;i<w.count;++i) { if(i) out << ','; Signed(out,w.codes[i]); }
    out << "],\"lengths\":[";
    for(unsigned i=0;i<w.count;++i) { if(i) out << ','; Number(out,w.lengths[i]); }
    out << "],\"first_success\":"; Number(out,w.first_success);
    out << ",\"bytes\":"; Number(out,w.bytes);
    out << ",\"handle\":"; Number(out,w.handle);
    out << ",\"freed\":" << (w.freed?"true":"false")
        << ",\"complete\":" << (w.complete?"true":"false") << '}';
}
struct Record {
    Coordinate c;
    uint64_t target=kMissing,ready=kMissing;
    std::array<uint64_t,4> wait={{kMissing,kMissing,kMissing,kMissing}};
    Observation observation;
    WorkResult work;
    bool called=false,checked=false,has_snapshot=false;
    Snapshot snapshot;
    explicit Record(const Coordinate& coordinate):c(coordinate) {}
};
void RecordJson(std::ostream& out,const Record& r) {
    const Coordinate& c=r.c;
    out << "{\"type\":\"record\",\"index\":" << c.index << ",\"rep\":" << c.rep
        << ",\"order\":" << c.order << ",\"metric\":" << c.metric << ",\"class\":" << c.comparison
        << ",\"position\":" << c.position << ",\"arm\":" << c.arm << ",\"q\":" << c.q << ",\"target\":";
    Number(out,r.target); out << ",\"ready\":"; Number(out,r.ready); out << ",\"wait\":[";
    for(unsigned i=0;i<4;++i) { if(i) out << ','; Number(out,r.wait[i]); }
    out << "],\"observation\":"; ObservationJson(out,r.observation);
    out << ",\"work\":"; WorkJson(out,r.work);
    out << ",\"checked\":" << (r.checked?"true":"false") << ",\"snapshot\":";
    if(r.has_snapshot) SnapshotJson(out,r.snapshot); else out << "null";
    out << '}';
}
void CheckWork(Fixture& fixture,Record& record,Handle& temporary) {
    const Coordinate& c=record.c; const WorkResult& w=record.work;
    Require(w.complete && w.handle>0 && w.handle!=kMissing,"callback incomplete public work");
    if(c.metric!=1) Require(w.create_code==0,"fresh constructor failure");
    if(c.metric==0 || c.metric==2) {
        Require(w.profile_bytes==(c.arm?32u:0u) && w.profile==fixture.profiles[c.arm],"fresh selected profile differs");
    }
    if(c.metric==3) {
        Require(w.first_success==fixture.first_success[c.arm] && w.count==fixture.decode_statuses[c.arm].size(),
                "decoder first-success changed");
        for(unsigned i=0;i<w.count;++i)
            Require(w.codes[i]==fixture.decode_statuses[c.arm][i] && w.lengths[i]==fixture.B,"decoder status/length changed");
        Require(w.recover_code==0 && w.bytes==fixture.M && w.calls==w.count+3u && w.freed,"decoder lifecycle accounting");
        fixture.CheckOutput(c.order,c.arm,0,0,true);
    } else {
        const unsigned count=c.metric==0?0u:c.metric==1?6u:unsigned(fixture.packet_count);
        Require(w.count==count,"encoder call-result count");
        for(unsigned i=0;i<count;++i) Require(w.codes[i]==0 && w.lengths[i]==fixture.B,"encode status/length");
        const uint64_t calls=c.metric==0?1u:c.metric==1?384u:fixture.packet_count+2u;
        Require(w.calls==calls && w.freed==(c.metric==2),"encoder lifecycle accounting");
        fixture.CheckOutput(c.order,c.arm,count,c.metric==1?fixture.K:0u);
    }
    const Buffer& other=fixture.outputs[1u-c.order];
    for(size_t i=0;i<other.size;++i) Require(other.data[i]==0xa5u,"unselected output lane changed");
    Require(fixture.packets==fixture.packet_reference,"immutable packet corpus changed");
    fixture.CheckSource();
    if(c.metric==0) {
        Require(temporary.value && Address(temporary.value)==w.handle,"post-create handle lost");
        record.snapshot=fixture.Inspect(c.arm,c.rep,c.order,temporary.value,w.profile.data()); record.has_snapshot=true;
        // Standalone-create cleanup is outside its component bracket.  The
        // lifecycle metric instead frees within RunWork and charges that work.
        void* value=temporary.value; temporary.value=nullptr; kApis[c.arm].free(value);
    } else Require(!temporary.value,"timed lifecycle leaked a handle");
    record.checked=true;
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

void HeaderJson(std::ostream& out,const Fixture& fixture,const std::string& claim,
                uint64_t start,uint64_t cpu_start,const std::string& identity,const std::string& runtime,
                const std::array<uint64_t,2>& resolutions,const Observation& prelude,uint64_t final_state) {
    out << "{\"type\":\"header\",\"protocol\":" << Quote(kProtocol)
        << ",\"schema\":" << Quote(std::string(kProtocol)+".raw.v1") << ",\"claim_sha256\":" << Quote(claim)
        << ",\"target_cpu\":50,\"worker_start_ns\":" << start << ",\"worker_start_cpu_ns\":" << cpu_start
        << ",\"identity_before\":" << identity << ",\"runtime_before\":" << runtime
        << ",\"clock_resolutions\":[" << resolutions[0] << ',' << resolutions[1] << ']'
        << ",\"source_hex\":" << Quote(Hex(fixture.Source(),fixture.M))
        << ",\"source_sha256\":" << Quote(fixture.source_hash)
        << ",\"outputs\":[" << Address(fixture.outputs[0].data) << ',' << Address(fixture.outputs[1].data)
        << "],\"preflight\":{\"packets_hex\":[";
    for(unsigned i=0;i<3;++i) { if(i) out << ','; out << Quote(Hex(fixture.packets[i].data(),fixture.packets[i].size())); }
    out << "],\"packet_sha256\":[";
    for(unsigned i=0;i<3;++i) { if(i) out << ','; out << Quote(Hash(fixture.packets[i].data(),fixture.packets[i].size())); }
    out << "],\"profiles_hex\":[null," << Quote(Hex(fixture.profiles[1].data(),32u)) << ','
        << Quote(Hex(fixture.profiles[2].data(),32u)) << "],\"intermediate_hex\":[";
    for(unsigned i=0;i<2;++i) { if(i) out << ','; out << Quote(Hex(fixture.intermediate[i].data(),fixture.intermediate[i].size())); }
    out << "],\"intermediate_sha256\":[" << Quote(fixture.intermediate_hash[0]) << ','
        << Quote(fixture.intermediate_hash[1]) << "],\"columns\":[";
    for(unsigned arm=0;arm<2;++arm) {
        if(arm) out << ',';
        out << '[';
        for(size_t i=0;i<fixture.columns[arm].size();++i) {
            if(i) out << ',';
            out << '[';
            const auto& row=fixture.columns[arm][i];
            for(size_t j=0;j<row.size();++j) { if(j) out << ','; out << row[j]; }
            out << ']';
        }
        out << ']';
    }
    out << "],\"decode_statuses\":[";
    for(unsigned arm=0;arm<3;++arm) {
        if(arm) out << ',';
        out << '[';
        for(size_t i=0;i<fixture.decode_statuses[arm].size();++i) { if(i) out << ','; out << fixture.decode_statuses[arm][i]; }
        out << ']';
    }
    out << "],\"first_success\":[" << fixture.first_success[0] << ',' << fixture.first_success[1] << ','
        << fixture.first_success[2] << "],\"snapshots\":[";
    for(size_t i=0;i<fixture.snapshots.size();++i) { if(i) out << ','; SnapshotJson(out,fixture.snapshots[i]); }
    out << "]},\"prelude\":{\"iterations\":" << kPreludeIterations << ",\"seed\":" << kPreludeSeed
        << ",\"final_state\":" << final_state << ",\"observation\":";
    ObservationJson(out,prelude); out << "}}";
}
struct Accounting {
    uint64_t records=0,callbacks=0,calls=0,checked=0,work=0,wait_wall=0,wait_cpu=0;
    void Add(const Record& r) {
        ++records; callbacks+=r.called?1u:0u; calls+=r.work.calls; checked+=r.checked?1u:0u;
        const Observation& o=r.observation;
        if(o.m1!=kMissing && o.m2!=kMissing && o.m2>=o.m1) work+=o.m2-o.m1;
        if(r.wait[0]!=kMissing && r.wait[2]!=kMissing && r.wait[2]>=r.wait[0]) wait_wall+=r.wait[2]-r.wait[0];
        if(r.wait[1]!=kMissing && r.wait[3]!=kMissing && r.wait[3]>=r.wait[1]) wait_cpu+=r.wait[3]-r.wait[1];
    }
};
void FooterJson(std::ostream& out,const Accounting& a,const std::string& failure,
                uint64_t schedule_now,uint64_t t0,uint64_t end,uint64_t cpu_end,
                const std::string& identity,const std::string& runtime) {
    out << "{\"type\":\"footer\",\"protocol\":" << Quote(kProtocol)
        << ",\"outcome\":" << Quote(failure.empty()?"COMPLETE":"INVALID")
        << ",\"failure\":" << (failure.empty()?"null":Quote(failure)) << ",\"schedule_now_ns\":";
    Number(out,schedule_now); out << ",\"t0\":"; Number(out,t0);
    out << ",\"records\":" << a.records << ",\"callbacks\":" << a.callbacks << ",\"work_calls\":" << a.calls
        << ",\"checked_callbacks\":" << a.checked << ",\"sum_work_ns\":" << a.work
        << ",\"sum_wait_wall_ns\":" << a.wait_wall << ",\"sum_wait_cpu_ns\":" << a.wait_cpu
        << ",\"worker_end_ns\":"; Number(out,end);
    out << ",\"worker_end_cpu_ns\":"; Number(out,cpu_end);
    out << ",\"identity_after\":" << (identity.empty()?"null":identity)
        << ",\"runtime_after\":" << (runtime.empty()?"null":runtime) << '}';
}

void Authenticate(const std::string& expected) {
    Require(expected.size()==64u && std::all_of(expected.begin(),expected.end(),[](char c) {
        return (c>='0' && c<='9') || (c>='a' && c<='f'); }),"claim SHA argument");
    const int fd=open(kClaim,O_RDONLY|O_NOFOLLOW|O_NONBLOCK|O_CLOEXEC);
    Require(fd>=0,"exclusive claim absent");
    struct Close { int fd; ~Close(){close(fd);} } cleanup={fd};
    struct stat before={},after={};
    Require(fstat(fd,&before)==0 && S_ISREG(before.st_mode) && before.st_nlink==1 &&
            before.st_size>0 && before.st_size<=1024*1024,"claim regular bounded file");
    std::string data(size_t(before.st_size),'\0'); size_t done=0;
    while(done<data.size()) {
        const ssize_t got=read(fd,&data[done],data.size()-done);
        if(got<0 && errno==EINTR) continue;
        Require(got>0,"claim short read"); done+=size_t(got);
    }
    char extra=0; Require(read(fd,&extra,1)==0 && fstat(fd,&after)==0,"claim grew");
    Require(before.st_dev==after.st_dev && before.st_ino==after.st_ino && before.st_size==after.st_size &&
            before.st_nlink==after.st_nlink && before.st_mtim.tv_sec==after.st_mtim.tv_sec &&
            before.st_mtim.tv_nsec==after.st_mtim.tv_nsec && before.st_ctim.tv_sec==after.st_ctim.tv_sec &&
            before.st_ctim.tv_nsec==after.st_ctim.tv_nsec && Sha256Hex(data)==expected,"claim changed/hash mismatch");
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
    Reader reader; Accounting accounting;
    uint64_t start=kMissing,cpu_start=kMissing,schedule_now=kMissing,t0=kMissing,end=kMissing,cpu_end=kMissing;
    std::string failure,identity_before,runtime_before,identity_after,runtime_after;
    try {
        Require(WH2_ALIGNMENT_COST_SANITIZERS==0,"sanitizer build cannot run scientific worker");
        Authenticate(claim); Limits(); start=reader.Now(); cpu_start=reader.Thread(); Pin();
        Require(wirehair_init()==Wirehair_Success,"wirehair initialization");
        identity_before=Identity(); runtime_before=Runtime();
        const std::array<uint64_t,2> resolutions={{reader.Resolution(CLOCK_MONOTONIC),reader.Resolution(CLOCK_THREAD_CPUTIME_ID)}};
        Require(resolutions[0]==1u && resolutions[1]==1u,"frozen clock resolutions");
        Fixture fixture(6,1280,12); fixture.Preflight();
        const auto roster=Roster(); Require(roster.size()==4800u,"complete roster size");
        Budget(reader,start,cpu_start,0u);
        Observation prelude; uint64_t final_state=0;
        const char* error=Capture(reader,[&] { final_state=Prelude(kPreludeSeed,kPreludeIterations); },prelude);
        Require(!error,error?error:"prelude capture"); Validate(prelude,nullptr);
        Require(final_state==kPreludeFinal,"prelude checksum");
        HeaderJson(std::cout,fixture,claim,start,cpu_start,identity_before,runtime_before,resolutions,prelude,final_state);
        std::cout << '\n' << std::flush; Require(bool(std::cout),"header stream write");
        schedule_now=reader.Now();
        Require(schedule_now/2000000u<(kMissing/2000000u)-2u,"schedule origin overflow");
        t0=(schedule_now/2000000u+2u)*2000000u;
        Observation previous=prelude;
        for(const Coordinate& c : roster) {
            Record record(c); Handle temporary;
            try {
                record.target=Target(t0,c);
                for(auto& output : fixture.outputs) std::memset(output.data,0xa5,output.size);
                // All allocation, serialization and checking from the prior
                // callback has finished before this readiness timestamp.
                Budget(reader,start,cpu_start,accounting.work);
                record.ready=reader.Now(); CheckReady(record.ready,record.target,previous);
                Wait(reader,record.target,start,record.wait);
                error=Capture(reader,[&] { RunWork(fixture,c,temporary,record.work); record.called=true; },record.observation);
                Require(!error,error?error:"callback capture");
                Validate(record.observation,&previous); CheckStart(record);
                CheckWork(fixture,record,temporary);
                Budget(reader,start,cpu_start,accounting.work+record.observation.m2-record.observation.m1);
                Require(OnCpu(),"worker migrated/affinity changed");
                previous=record.observation;
            } catch(...) {
                accounting.Add(record); RecordJson(std::cout,record); std::cout << '\n' << std::flush;
                throw;
            }
            accounting.Add(record); RecordJson(std::cout,record); std::cout << '\n' << std::flush;
            Require(bool(std::cout),"record stream write");
        }
        Require(accounting.records==4800u && accounting.callbacks==4800u && accounting.checked==4800u,"complete callback accounting");
        for(unsigned r=0;r<12;++r) for(unsigned lane=0;lane<2;++lane) for(unsigned arm=0;arm<3;++arm) {
            Handle& h=fixture.At(r,lane,arm);
            fixture.Inspect(arm,r,lane,h.value,h.profile.data());
        }
        fixture.CheckSource(); Require(OnCpu(),"final CPU affinity");
        identity_after=Identity(); runtime_after=Runtime();
        Require(runtime_after==runtime_before,"GF256 runtime changed");
        Budget(reader,start,cpu_start,accounting.work);
        // Destruction of all72 precreated handles is part of worker elapsed,
        // but deliberately outside the repair component measurements.
    } catch(const std::exception& error) { failure=error.what(); }
      catch(...) { failure="unexpected worker exception"; }
    if(start!=kMissing) {
        if(!reader.Mono(end) || !reader.Cpu(cpu_end)) {
            if(failure.empty()) failure="footer clock capture";
        } else if(end<start || end-start>=kWall || cpu_end<cpu_start || cpu_end-cpu_start>kCpu) {
            if(failure.empty()) failure="final cleanup resource cap";
        }
    }
    FooterJson(std::cout,accounting,failure,schedule_now,t0,end,cpu_end,identity_after,runtime_after);
    std::cout << '\n' << std::flush;
    return failure.empty() && bool(std::cout)?0:1;
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
            !OrderedClock(999999u,last_cpu) && last_cpu==1000000u,
            "neutral final/intermediate clock regression above worker start");
    Require(OrderedClock(1000u,last_mono) && OrderedClock(1000001u,last_cpu),
            "neutral independent clock domains");
    FakeReader reader; Observation o; unsigned work=0;
    Require(Capture(reader,[&]{++work;},o)==nullptr && work==1u,"neutral complete capture");
    Validate(o,nullptr); // CPU epoch differs and its span exceeds inner wall.
    for(unsigned failure=0;failure<8;++failure) {
        FakeReader bad; bad.fail=failure; Observation partial; unsigned completed=0;
        Require(Capture(bad,[&]{++completed;},partial)!=nullptr,"neutral failed clock capture");
        Require(completed==(failure>=4u?1u:0u),"neutral partial callback boundary");
    }
    Observation later=o; later.m0=1u;
    bool rejected=false; try { Validate(later,&o); } catch(const std::runtime_error&) { rejected=true; }
    Require(rejected,"neutral reverse clock rejection");
    const auto roster=Roster(); Require(roster.size()==4800u,"neutral roster size");
    unsigned seen[4][5][2][12]={};
    uint64_t minimum=kMissing,maximum=0;
    for(size_t i=0;i<roster.size();++i) {
        const auto& c=roster[i];
        Require(c.index==i && c.arm<3 && c.metric<4 && c.comparison<5 && c.order<2 && c.rep<12,"neutral coordinate bounds");
        if(c.position==0) ++seen[c.metric][c.comparison][c.order][c.rep];
        if(c.position>=2 && c.position%2u==0) {
            const auto& next=roster[i+1];
            Require(next.q==c.q && next.position==c.position+1u &&
                    (kSides[next.position]^next.order)!=(kSides[c.position]^c.order),"neutral matched phase pair");
        }
        if(i) { const uint64_t gap=Target(10000000u,c)-Target(10000000u,roster[i-1]);
            minimum=std::min(minimum,gap); maximum=std::max(maximum,gap); }
    }
    for(unsigned m=0;m<4;++m) for(unsigned c=0;c<5;++c) for(unsigned order=0;order<2;++order)
        for(unsigned rep=0;rep<12;++rep) Require(seen[m][c][order][rep]==1u,"neutral exact panel coverage");
    Require(minimum==1250000u && maximum==2250000u,"neutral nonuniform slot gaps");
    Record r(roster[0]); r.target=10000000u; r.observation.m1=r.target+5000u; CheckStart(r);
    r.observation.m1++;
    rejected=false; try { CheckStart(r); } catch(const std::runtime_error&) { rejected=true; }
    Require(rejected,"neutral late start rejection");
    Observation previous; previous.m3=1;
    CheckReady(r.target-100000u,r.target,previous);
    rejected=false; try { CheckReady(r.target-99999u,r.target,previous); } catch(const std::runtime_error&) { rejected=true; }
    Require(rejected,"neutral preparation deadline rejection");
}
void NeutralBridges(Fixture& fixture) {
    for(unsigned arm=1;arm<3;++arm) {
        const auto inspect=arm==1?wh2_aligned_cost_p_inspect:wh2_aligned_cost_a_inspect;
        const auto query=arm==1?wh2_aligned_cost_p_packet_columns:wh2_aligned_cost_a_packet_columns;
        const auto scalar=arm==1?wh2_aligned_cost_p_scalar_packet:wh2_aligned_cost_a_scalar_packet;
        WirehairV2Codec h=static_cast<WirehairV2Codec>(fixture.At(0,0,arm).value);
        Wh2AlignedIntermediateCostSnapshot snapshot={}; snapshot.message_bytes=123u;
        Require(inspect(nullptr,&snapshot)==WirehairV2_InvalidInput && snapshot.message_bytes==123u,"neutral null inspection");
        Require(inspect(h,nullptr)==WirehairV2_InvalidInput,"neutral null inspection output");
        uint32_t bytes=UINT32_MAX;
        Require(scalar(h,fixture.K,nullptr,0,&bytes)==WirehairV2_BufferTooSmall && bytes==fixture.B,"neutral scalar size query");
        std::vector<uint8_t> guarded(fixture.B+128u,0xa5u); bytes=UINT32_MAX;
        Require(scalar(h,fixture.K,guarded.data()+64u,fixture.B-1u,&bytes)==WirehairV2_BufferTooSmall && bytes==fixture.B &&
                std::all_of(guarded.begin(),guarded.end(),[](uint8_t value){return value==0xa5u;}),"neutral scalar short buffer");
        uint32_t count=UINT32_MAX;
        Require(query(h,fixture.K,nullptr,0,&count)==WirehairV2_BufferTooSmall && count>0u,"neutral column size query");
        std::vector<uint32_t> row(count+2u,UINT32_MAX); const uint32_t required=count; count=UINT32_MAX;
        Require(query(h,fixture.K,row.data()+1u,required-1u,&count)==WirehairV2_BufferTooSmall && count==required &&
                std::all_of(row.begin(),row.end(),[](uint32_t value){return value==UINT32_MAX;}),"neutral column short buffer");
        // Aligned separate storage avoids making a misaligned uint32_t pointer.
        std::array<uint32_t,64> alias={{0}}; alias.fill(UINT32_MAX);
        Require(query(h,fixture.K,alias.data(),32u,alias.data())==WirehairV2_InvalidInput &&
                std::all_of(alias.begin(),alias.end(),[](uint32_t value){return value==UINT32_MAX;}),"neutral column output alias");
        Require(scalar(h,fixture.K,alias.data(),fixture.B,alias.data())==WirehairV2_InvalidInput &&
                std::all_of(alias.begin(),alias.end(),[](uint32_t value){return value==UINT32_MAX;}),"neutral scalar output alias");
    }
}
void NeutralMetrics(Fixture& fixture) {
    for(unsigned metric=0;metric<4;++metric) for(unsigned arm=0;arm<3;++arm) {
        const Coordinate coordinate={0,0,arm%2u,metric,0,0,arm,0};
        Record record(coordinate); Handle temporary;
        for(auto& output : fixture.outputs) std::memset(output.data,0xa5,output.size);
        FakeReader reader;
        Require(Capture(reader,[&] { RunWork(fixture,coordinate,temporary,record.work); record.called=true; },
                        record.observation)==nullptr,"neutral metric synthetic capture");
        Validate(record.observation,nullptr); CheckWork(fixture,record,temporary);
        Require(record.checked,"neutral real API metric correctness");
        // A cross-lane write is detected before a later reset can hide it.
        fixture.outputs[1u-coordinate.order].data[0]^=1u;
        bool rejected=false;
        try { fixture.CheckOutput(1u-coordinate.order,arm,0,0); }
        catch(const std::runtime_error&) { rejected=true; }
        Require(rejected,"neutral cross-lane guard rejection");
        fixture.outputs[1u-coordinate.order].data[0]^=1u;
    }
}
int Neutral() {
    Require(wirehair_init()==Wirehair_Success,"neutral global initialization");
    NeutralClocks();
    for(unsigned K : {7u,9u}) for(unsigned B : {17u,65u}) {
        Fixture fixture(K,B,1); fixture.Preflight(); NeutralBridges(fixture); NeutralMetrics(fixture);
    }
    std::cout << "PASS neutral K7/K9 B17/65, byte/profile/decoder/bridge and synthetic-clock/roster checks\n";
    return 0;
}

} // namespace cost

int main(int argc,char** argv) {
    try {
        if(argc==2 && std::strcmp(argv[1],"--neutral")==0) return cost::Neutral();
        if(argc==3 && std::strcmp(argv[1],"--worker")==0) return cost::Worker(argv[2]);
        std::cerr << "usage: --neutral | --worker CLAIM_SHA256 (no default execution)\n";
        return 2;
    } catch(const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n'; return 1;
    }
}
