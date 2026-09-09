// Exact pre/post shared-library admission regression screen.
// No Wirehair library is linked into this worker. Each DSO owns its GF state.
#include "wirehair/wirehair.h"
#include "wirehair/wirehair_small.h"
#include "wirehair/wirehair_k6.h"
#include "Wh2FrozenTrace.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "gf256.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#include <dlfcn.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sched.h>
#include <sys/resource.h>
#include <time.h>
#if !defined(WH2_ADMISSION_REGRESSION_NEUTRAL) || defined(WH_COUNT) || defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "Explicit unmodified native/neutral build required"
#endif
#define NOINLINE __attribute__((noinline, noipa))
#ifndef WH2_ADMISSION_PROTOCOL
#define WH2_ADMISSION_PROTOCOL "wirehair.wh2.admission-regression-cost-r0"
#endif
#ifndef WH2_ADMISSION_CLAIM_PATH
#define WH2_ADMISSION_CLAIM_PATH "/var/tmp/wh2-admission-regression-cost-r0/CLAIM.json"
#endif
namespace {
const char protocol[]=WH2_ADMISSION_PROTOCOL;
const char claim_path[]=WH2_ADMISSION_CLAIM_PATH;
const unsigned case_count=20,max_batch=128,callbacks=51840;
const unsigned pairs[3][2]={{0,0},{1,1},{0,1}};
const unsigned sides[18]={0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
using Profile=std::array<uint8_t,32>;
static_assert(sizeof(void*)==8 && sizeof(size_t)==8 && sizeof(WirehairV2EncoderOptions)==16 &&
              sizeof(WirehairSmallCreateResult)==16 && sizeof(WirehairK6CreateResult)==16 &&
              sizeof(WirehairSmallResult)==24 && sizeof(WirehairK6Result)==24,"native C ABI sizes");
void Check(bool value,const char* why) { if(!value) throw std::runtime_error(why); }
struct SymbolSpec { const char* name; uintptr_t offset; };
struct LibrarySpec {
    const char* path; const char* sha;
    std::array<SymbolSpec,53> exports; std::array<SymbolSpec,6> slots;
    std::array<SymbolSpec,37> runtime_slots;
    uintptr_t context,context_bytes,getter;
};
// Generated only from the hash-authenticated ELF tables by the build tool.
#include "AdmissionLibraryBindings.h"
std::string ReadFile(const char* path,size_t cap) {
    const int fd=open(path,O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    Check(fd>=0,"open regular input");
    struct Close { int fd; ~Close() { close(fd); } } close_fd{fd};
    struct stat st={}; Check(!fstat(fd,&st) && S_ISREG(st.st_mode) && st.st_nlink==1 &&
        st.st_size>=0 && uint64_t(st.st_size)<=cap,"bounded regular input");
    std::string out; out.resize(size_t(st.st_size)); size_t done=0;
    while(done<out.size()) { const ssize_t n=read(fd,&out[done],out.size()-done);
        Check(n>0,"complete input read"); done+=size_t(n); }
    char extra=0; Check(read(fd,&extra,1)==0,"input grew");
    return out;
}
uintptr_t Address(void* p) { return reinterpret_cast<uintptr_t>(p); }
struct Library {
    void* handle=nullptr; uintptr_t base=0;
    std::array<uintptr_t,37> runtime_targets={};
    explicit Library(unsigned index) {
        const auto& spec=library_specs[index];
        Check(wirehair::wh2_benchmark::Sha256Hex(ReadFile(spec.path,4u*1024u*1024u))==spec.sha,"exact DSO hash");
        Check(!dlsym(RTLD_DEFAULT,"wirehair_init_"),"no global Wirehair");
        handle=dlopen(spec.path,RTLD_NOW|RTLD_LOCAL); Check(handle!=nullptr,"load exact DSO");
        Dl_info info={}; void* init=dlsym(handle,"wirehair_init_");
        Check(init && dladdr(init,&info) && info.dli_fname && !strcmp(info.dli_fname,spec.path),"own DSO base");
        base=Address(info.dli_fbase);
        for(size_t i=0;i<runtime_targets.size();++i) {
            memcpy(&runtime_targets[i],reinterpret_cast<const void*>(base+spec.runtime_slots[i].offset),sizeof(uintptr_t));
            Check(runtime_targets[i]!=0,"resolved runtime GOT");
        }
        Validate(index);
        Check(reinterpret_cast<decltype(&wirehair_init_)>(init)(WIREHAIR_VERSION)==0,"own GF init");
        gf256_x86_cpu_features f={};
        reinterpret_cast<void(*)(gf256_x86_cpu_features*)>(base+spec.getter)(&f);
        Check(f.SSSE3==1 && f.AVX2==1 && f.GFNI==1 && f.AVX512==1,"active native GF");
        const auto* ctx=reinterpret_cast<const gf256_ctx*>(base+spec.context);
        Check(ctx->Polynomial==0x14d && spec.context_bytes==sizeof(gf256_ctx),"own GF context");
    }
    void Validate(unsigned index) const {
        const auto& spec=library_specs[index];
        Check(wirehair::wh2_benchmark::Sha256Hex(ReadFile(spec.path,4u*1024u*1024u))==spec.sha,"unchanged DSO hash");
        Check(!dlsym(RTLD_DEFAULT,"wirehair_init_"),"Wirehair remains local");
        for(const auto& s:spec.exports) Check(Address(dlsym(handle,s.name))==base+s.offset,"owned public export");
        for(const auto& s:spec.slots) {
            uintptr_t target=0; memcpy(&target,reinterpret_cast<const void*>(base+s.offset),sizeof(target));
            Check(target==Address(dlsym(handle,s.name)),"owned internal wirehair GOT");
        }
        for(size_t i=0;i<runtime_targets.size();++i) {
            uintptr_t target=0; memcpy(&target,reinterpret_cast<const void*>(base+spec.runtime_slots[i].offset),sizeof(target));
            Check(target==runtime_targets[i],"stable actual runtime GOT");
        }
    }
    Library(const Library&)=delete; Library& operator=(const Library&)=delete;
    // Deliberately keep providers loaded until process exit.
};
Library* libraries[2]={};
void Load(unsigned order) {
    const unsigned first=order?1:0,second=1-first;
    static Library a(first),b(second); libraries[first]=&a; libraries[second]=&b;
    const auto& x=library_specs[0]; const auto& y=library_specs[1];
    const auto xa=libraries[0]->base+x.context,ya=libraries[1]->base+y.context;
    Check(xa+x.context_bytes<=ya || ya+y.context_bytes<=xa,"distinct GF contexts");
    for(size_t i=0;i<x.runtime_slots.size();++i)
        Check(!strcmp(x.runtime_slots[i].name,y.runtime_slots[i].name),"common runtime import roster");
    Check(a.runtime_targets==b.runtime_targets,"same actual C/C++ runtime GOT providers");
    for(const char* name:{"malloc","free","memcpy","_Znwm","_ZdlPv"})
        Check(dlsym(a.handle,name) && dlsym(a.handle,name)==dlsym(b.handle,name),"shared runtime providers");
}
#define API_SYMBOLS(X) \
    X(wirehair_v2_encoder_create_profile_id_with_options) X(wirehair_v2_encode) \
    X(wirehair_v2_decoder_create) X(wirehair_v2_decode) X(wirehair_v2_recover) X(wirehair_v2_free) \
    X(wirehair_encoder_create_ex) X(wirehair_encoder_create_owned_ex) X(wirehair_encode) \
    X(wirehair_decoder_create_ex) X(wirehair_decode) X(wirehair_recover) X(wirehair_free) \
    X(wirehair_small_encoder_create) X(wirehair_small_encode) X(wirehair_small_decoder_create) \
    X(wirehair_small_decode) X(wirehair_small_recover) X(wirehair_small_free) \
    X(wirehair_k6_encoder_create) X(wirehair_k6_encode) X(wirehair_k6_decoder_create) \
    X(wirehair_k6_decode) X(wirehair_k6_recover) X(wirehair_k6_free)
struct Symbols {
#define DECLARE(name) decltype(&::name) name=nullptr;
    API_SYMBOLS(DECLARE)
#undef DECLARE
    void Bind(void* library) {
#define BIND(name) name=reinterpret_cast<decltype(name)>(dlsym(library,#name)); Check(name!=nullptr,"typed API symbol");
        API_SYMBOLS(BIND)
#undef BIND
    }
};
Symbols symbols[2];
enum Family { Certified,Small,Wh1,K6 };
struct Case { unsigned family,k,b,policy; };
const Case cases[case_count]={
    {Certified,2,2,2},{Certified,2,1280,2},{Certified,3,2,2},{Certified,3,1280,2},
    {Certified,4,2,2},{Certified,4,1280,2},{Certified,6,2,2},{Certified,6,1280,2},
    {Certified,128,2,2},{Certified,128,1280,2},{Certified,3,2,1},{Certified,3,1280,1},
    {Small,3,2,1},{Small,3,2,2},{Small,3,1280,1},{Small,3,1280,2},
    {Wh1,3,2,2},{Wh1,3,1280,2},{K6,6,2,2},{K6,6,1280,2}};
unsigned Batch(const Case& c) { return c.k==128?4:128; }
unsigned PacketCount(const Case& c) { return c.k+14; }
uint32_t Packet(const Case& c,unsigned slot) { return slot<c.k+8?slot:UINT32_MAX-2*(slot-c.k-8); }
unsigned Slot(const Case& c,unsigned step) { return step<14?c.k+step:step-14; }
template<class T> struct HandleResult {
    T value=nullptr; void*& output;
    explicit HandleResult(void*& h):output(h) {}
    ~HandleResult() { output=value; }
    HandleResult(const HandleResult&)=delete; HandleResult& operator=(const HandleResult&)=delete;
};
struct Api {
    const Symbols* s=nullptr; Case c={};
    int Create(const void* source,Profile& p,void*& h) const {
        const uint64_t m=uint64_t(c.k)*c.b;
        if(c.family==Certified) {
            WirehairV2EncoderOptions o=WIREHAIR_V2_ENCODER_OPTIONS_INIT; o.source_policy=c.policy;
            HandleResult<WirehairV2Codec> a(h); uint32_t n=0;
            const int r=s->wirehair_v2_encoder_create_profile_id_with_options(
                WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,source,m,c.b,&o,p.data(),32,&n,&a.value);
            return r==0 && n!=32?-1:r;
        }
        if(c.family==Small) { const auto r=s->wirehair_small_encoder_create(source,m,c.b,c.policy,p.data(),32);
            h=r.codec; return r.status; }
        if(c.family==K6) { const auto r=s->wirehair_k6_encoder_create(source,m,c.b,c.policy,p.data(),32);
            h=r.codec; return r.status; }
        HandleResult<WirehairCodec> a(h);
        return (c.policy==1?s->wirehair_encoder_create_owned_ex:s->wirehair_encoder_create_ex)(nullptr,source,m,c.b,&a.value);
    }
    int Encode(void* h,uint32_t id,void* out,uint32_t* n) const {
        if(c.family==Certified) return s->wirehair_v2_encode(static_cast<WirehairV2Codec>(h),id,out,c.b,n);
        if(c.family==Wh1) return s->wirehair_encode(static_cast<WirehairCodec>(h),id,out,c.b,n);
        if(c.family==Small) { const auto r=s->wirehair_small_encode(static_cast<WirehairSmallCodec>(h),id,out,c.b);
            *n=uint32_t(r.bytes_written); return r.status==0 && (r.bytes_written!=c.b || r.bytes_required!=c.b)?-1:r.status; }
        const auto r=s->wirehair_k6_encode(static_cast<WirehairK6Codec>(h),id,out,c.b);
        *n=uint32_t(r.bytes_written); return r.status==0 && (r.bytes_written!=c.b || r.bytes_required!=c.b)?-1:r.status;
    }
    int Decoder(const Profile& p,void*& h) const {
        if(c.family==Certified) { HandleResult<WirehairV2Codec> a(h);
            return s->wirehair_v2_decoder_create(p.data(),32,&a.value); }
        if(c.family==Small) { const auto r=s->wirehair_small_decoder_create(p.data(),32); h=r.codec; return r.status; }
        if(c.family==K6) { const auto r=s->wirehair_k6_decoder_create(p.data(),32); h=r.codec; return r.status; }
        HandleResult<WirehairCodec> a(h);
        return s->wirehair_decoder_create_ex(nullptr,uint64_t(c.k)*c.b,c.b,&a.value);
    }
    int Feed(void* h,uint32_t id,const void* in) const {
        if(c.family==Certified) return s->wirehair_v2_decode(static_cast<WirehairV2Codec>(h),id,in,c.b);
        if(c.family==Small) return s->wirehair_small_decode(static_cast<WirehairSmallCodec>(h),id,in,c.b);
        if(c.family==K6) return s->wirehair_k6_decode(static_cast<WirehairK6Codec>(h),id,in,c.b);
        return s->wirehair_decode(static_cast<WirehairCodec>(h),id,in,c.b);
    }
    int Recover(void* h,void* out,uint64_t* n) const {
        const uint64_t m=uint64_t(c.k)*c.b;
        if(c.family==Certified) return s->wirehair_v2_recover(static_cast<WirehairV2Codec>(h),out,m,n);
        if(c.family==Small) { const auto r=s->wirehair_small_recover(static_cast<WirehairSmallCodec>(h),out,size_t(m));
            *n=r.bytes_written; return r.status==0 && r.bytes_required!=m?-1:r.status; }
        if(c.family==K6) { const auto r=s->wirehair_k6_recover(static_cast<WirehairK6Codec>(h),out,size_t(m));
            *n=r.bytes_written; return r.status==0 && r.bytes_required!=m?-1:r.status; }
        const int r=s->wirehair_recover(static_cast<WirehairCodec>(h),out,m); *n=r==0?m:0; return r;
    }
    void Free(void* h) const {
        if(c.family==Certified) s->wirehair_v2_free(static_cast<WirehairV2Codec>(h));
        else if(c.family==Small) s->wirehair_small_free(static_cast<WirehairSmallCodec>(h));
        else if(c.family==K6) s->wirehair_k6_free(static_cast<WirehairK6Codec>(h));
        else s->wirehair_free(static_cast<WirehairCodec>(h));
    }
};
struct Owner {
    const Api& api; void* h=nullptr;
    explicit Owner(const Api& a):api(a) {}
    ~Owner() { if(h) api.Free(h); }
    Owner(const Owner&)=delete; Owner& operator=(const Owner&)=delete;
};
struct Arm { Profile profile={}; std::vector<uint8_t> packets; unsigned steps=0; };
struct Fixture { Case c={}; std::vector<uint8_t> source; Arm arm[2]; Api api[2]; };
Fixture fixtures[case_count],reference[case_count];
alignas(64) uint8_t outputs[2][max_batch*(20*1280+128)];
size_t Stride(const Case& c) { return size_t(PacketCount(c))*c.b+128; }
void Initialize(unsigned order) {
    Load(order);
    for(unsigned a=0;a<2;++a) symbols[a].Bind(libraries[a]->handle);
    for(unsigned i=0;i<case_count;++i) {
        auto& f=fixtures[i]; f.c=cases[i]; const auto& c=f.c; f.source.resize(size_t(c.k)*c.b);
        for(size_t j=0;j<f.source.size();++j) f.source[j]=uint8_t(37*j+j/11);
        for(unsigned a=0;a<2;++a) {
            f.api[a].s=&symbols[a]; f.api[a].c=c; const Api& api=f.api[a]; Arm& arm=f.arm[a];
            arm.packets.resize(size_t(PacketCount(c))*c.b);
            { auto source=f.source; Owner e(api);
              Check(api.Create(source.data(),arm.profile,e.h)==0 && e.h,"fixture encoder");
              if(c.policy==1) std::fill(source.begin(),source.end(),0xcc);
              for(unsigned j=0;j<PacketCount(c);++j) {
                  std::vector<uint8_t> out(c.b+128,0xa5); uint32_t n=0;
                  Check(api.Encode(e.h,Packet(c,j),out.data()+64,&n)==0 && n==c.b,"fixture packet");
                  Check(std::all_of(out.begin(),out.begin()+64,[](uint8_t v){return v==0xa5;}) &&
                        std::all_of(out.end()-64,out.end(),[](uint8_t v){return v==0xa5;}),"fixture packet guards");
                  memcpy(arm.packets.data()+size_t(j)*c.b,out.data()+64,c.b);
              }
              Check(c.policy==1?std::all_of(source.begin(),source.end(),[](uint8_t v){return v==0xcc;}):source==f.source,
                    "immutable fixture source");
            }
            { Owner d(api); Check(api.Decoder(arm.profile,d.h)==0 && d.h,"standalone fixture decoder");
              for(unsigned j=0;j<PacketCount(c);++j) {
                  const unsigned slot=Slot(c,j); const int r=api.Feed(d.h,Packet(c,slot),arm.packets.data()+size_t(slot)*c.b);
                  Check(r==0 || r==1,"fixture feed");
                  if(r==0) { Check(j+1>=c.k,"premature decode"); arm.steps=j+1; break; }
              }
              Check(arm.steps>0,"fixture endpoint");
              for(unsigned j=0;j<2;++j) { std::vector<uint8_t> out(f.source.size()+128,0xa5); uint64_t n=0;
                  Check(api.Recover(d.h,out.data()+64,&n)==0 && n==f.source.size() &&
                        !memcmp(out.data()+64,f.source.data(),f.source.size()),"fixture recovered source");
                  Check(std::all_of(out.begin(),out.begin()+64,[](uint8_t v){return v==0xa5;}) &&
                        std::all_of(out.end()-64,out.end(),[](uint8_t v){return v==0xa5;}),"fixture recovery guards");
              }
            }
        }
        Check(f.arm[0].profile==f.arm[1].profile && f.arm[0].packets==f.arm[1].packets &&
              f.arm[0].steps==f.arm[1].steps,"pre/post exact fixture");
        Check(Batch(c)*Stride(c)<=sizeof(outputs[0]),"output extent");
        reference[i]=f;
    }
}
struct Coordinate { unsigned index,rep,order,which,metric,comparison,position,arm; uint64_t q; };
Coordinate CoordinateAt(unsigned index) {
    unsigned n=index; const unsigned p=n%18; n/=18; const unsigned cs=n%3; n/=3;
    const unsigned ms=n%2; n/=2; const unsigned ws=n%case_count; n/=case_count; const unsigned s=n%2,r=n/2;
    const unsigned o=(r+s)%2,w=(r+s+ws)%case_count,m=(r+s+ws+ms)%2,c=(2*r+s+ws+m+cs)%3;
    const unsigned bin=p<2?r+12*(r%4):((r+6*((p-2)/8))%12)+12*(((p-2)%8)/2);
    return Coordinate{index,r,o,w,m,c,p,pairs[c][sides[p]^o],uint64_t(2*bin+1)*1000000/96};
}
using Counters=std::array<uint64_t,4>;
struct Observation { uint64_t m0=0,c0=0,m1=0,m2=0,c1=0,m3=0; Counters before={},after={}; };
struct Reader {
    uint64_t Clock(clockid_t clock) {
        timespec t={}; Check(clock_gettime(clock,&t)==0 && t.tv_sec>=0 && t.tv_nsec>=0 && t.tv_nsec<1000000000,"clock");
        Check(uint64_t(t.tv_sec)<UINT64_MAX/1000000000-1,"clock range");
        return uint64_t(t.tv_sec)*1000000000+uint64_t(t.tv_nsec);
    }
    uint64_t Mono() { return Clock(CLOCK_MONOTONIC); }
    uint64_t Cpu() { return Clock(CLOCK_THREAD_CPUTIME_ID); }
    Counters Usage() {
        rusage r={}; Check(getrusage(RUSAGE_THREAD,&r)==0 && r.ru_minflt>=0 && r.ru_majflt>=0 && r.ru_nvcsw>=0 && r.ru_nivcsw>=0,"counters");
        return {{uint64_t(r.ru_minflt),uint64_t(r.ru_majflt),uint64_t(r.ru_nvcsw),uint64_t(r.ru_nivcsw)}};
    }
};
template<class R,class F> void Capture(R& reader,F function,Observation& o) {
    o.m0=reader.Mono(); o.before=reader.Usage(); o.c0=reader.Cpu(); o.m1=reader.Mono();
    std::atomic_signal_fence(std::memory_order_seq_cst); function();
    std::atomic_signal_fence(std::memory_order_seq_cst);
    o.m2=reader.Mono(); o.c1=reader.Cpu(); o.after=reader.Usage(); o.m3=reader.Mono();
}
void Validate(const Observation& o,const Observation& previous) {
    Check(o.m0<=o.m1 && o.m1<o.m2 && o.m2<=o.m3 && o.c0<=o.c1 &&
          o.c1-o.c0<=o.m3-o.m0 && previous.m3<=o.m0 && previous.c1<=o.c0,"clock ordering");
    for(unsigned j=0;j<4;++j) Check(previous.after[j]<=o.before[j] && o.before[j]<=o.after[j],"counter ordering");
}

struct Work {
    std::array<unsigned,6> counts={}; std::array<uint64_t,max_batch> addresses={};
    unsigned address_count=0; bool complete=false;
};
struct Record { Coordinate c={}; uint64_t ready=0,target=0; std::array<uint64_t,4> wait={}; Observation o; Work work; bool checked=false; };
Record records[callbacks]; unsigned retained=0;
NOINLINE void RunWork(const Api& api,const Fixture& f,const Arm& arm,unsigned metric,uint8_t* output,Work& w) noexcept {
    void* handle=nullptr;
    auto release=[&] { if(handle) { void* owned=handle; handle=nullptr; ++w.counts[5]; api.Free(owned); } };
    try {
        for(unsigned cycle=0;cycle<Batch(f.c);++cycle) {
            uint8_t* out=output+cycle*Stride(f.c)+64; Profile p={}; int result;
            if(metric==0) { ++w.counts[0]; result=api.Create(f.source.data(),p,handle); }
            else { ++w.counts[2]; result=api.Decoder(arm.profile,handle); p=arm.profile; }
            w.addresses[w.address_count++]=uint64_t(reinterpret_cast<uintptr_t>(handle));
            if(result!=0 || !handle || p!=arm.profile) { release(); return; }
            if(metric==0) {
                for(unsigned j=0;j<PacketCount(f.c);++j) {
                    uint32_t n=0; ++w.counts[1]; result=api.Encode(handle,Packet(f.c,j),out+size_t(j)*f.c.b,&n);
                    if(result!=0 || n!=f.c.b) { release(); return; }
                }
            } else {
                bool success=false;
                for(unsigned j=0;j<PacketCount(f.c);++j) {
                    const unsigned slot=Slot(f.c,j); ++w.counts[3];
                    result=api.Feed(handle,Packet(f.c,slot),arm.packets.data()+size_t(slot)*f.c.b);
                    if(j>=arm.steps || result!=(j+1==arm.steps?0:1)) { release(); return; }
                    if(result==0) { success=true; break; }
                }
                if(!success) { release(); return; }
                uint64_t n=0; ++w.counts[4]; result=api.Recover(handle,out,&n);
                if(result!=0 || n!=f.source.size()) { release(); return; }
            }
            release();
        }
        w.complete=true;
    } catch(...) { try { release(); } catch(...) {} }
}
void Prepare(const Fixture& f) { for(auto& out:outputs) memset(out,0xa5,Batch(f.c)*Stride(f.c)); }
void CheckWork(const Fixture& f,const Arm& arm,unsigned metric,unsigned lane,const Work& w) {
    const unsigned batch=Batch(f.c);
    const std::array<unsigned,6> expected={{metric?0:batch,metric?0:batch*PacketCount(f.c),
        metric?batch:0,metric?batch*arm.steps:0,metric?batch:0,batch}};
    Check(w.complete && w.counts==expected && w.address_count==batch,"full lifecycle ledger");
    for(unsigned j=0;j<max_batch;++j) Check(j<batch?w.addresses[j]!=0:w.addresses[j]==0,"address roster");
    const auto& bytes=metric?f.source:arm.packets; const size_t stride=Stride(f.c);
    for(unsigned cycle=0;cycle<batch;++cycle) {
        const uint8_t* out=outputs[lane]+cycle*stride;
        Check(!memcmp(out+64,bytes.data(),bytes.size()),"every lifecycle output");
        for(size_t j=0;j<stride;++j) if(j<64 || j>=64+bytes.size()) Check(out[j]==0xa5,"output guards");
    }
    for(size_t j=0;j<batch*stride;++j) Check(outputs[1-lane][j]==0xa5,"unused lane");
    for(unsigned i=0;i<case_count;++i) {
        Check(fixtures[i].source==reference[i].source,"immutable source");
        for(unsigned a=0;a<2;++a) Check(fixtures[i].arm[a].profile==reference[i].arm[a].profile &&
            fixtures[i].arm[a].packets==reference[i].arm[a].packets && fixtures[i].arm[a].steps==reference[i].arm[a].steps,
            "immutable packet corpus");
    }
}
bool OnCpu() {
    cpu_set_t mask; CPU_ZERO(&mask);
    return sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(50,&mask) && sched_getcpu()==50;
}
void Pin() {
    cpu_set_t mask; CPU_ZERO(&mask); CPU_SET(50,&mask);
    Check(sched_setaffinity(0,sizeof(mask),&mask)==0 && OnCpu(),"CPU50 affinity");
}

std::string Identity() {
    namespace P=wirehair_wh2_bench;
    P::TargetIdentityReceiptV2 r; std::string why,canonical;
    Check(P::CapturePublicBorrowedTargetIdentity(50,r,why) && P::SerializeTargetIdentityV2(r,canonical,why),"target identity");
    Check(canonical.size()==617 && wirehair::wh2_benchmark::Sha256Hex(canonical)==
          "3288e0ef61cf3e628dcd827f9cf003c9d6ec6b5a12169e7a8bfc796baacddba7","exact CPU50 identity");
    return canonical;
}
NOINLINE uint64_t Prelude(uint64_t value,unsigned count) {
    for(unsigned i=0;i<count;++i) { value^=value<<13; value^=value>>7; value^=value<<17; } return value;
}
void Hex(const void* data,size_t bytes) {
    const auto* p=static_cast<const uint8_t*>(data); putchar('"');
    for(size_t i=0;i<bytes;++i) { printf("%02x",unsigned(p[i])); }
    putchar('"');
}
template<class T,size_t N> void Array(const std::array<T,N>& a) {
    putchar('['); for(size_t i=0;i<N;++i) { if(i) putchar(','); printf("%llu",static_cast<unsigned long long>(a[i])); } putchar(']');
}
void ObservationJson(const Observation& o) {
    const std::array<uint64_t,6> clocks={{o.m0,o.c0,o.m1,o.m2,o.c1,o.m3}};
    printf("{\"clocks\":"); Array(clocks); printf(",\"before\":"); Array(o.before); printf(",\"after\":"); Array(o.after); putchar('}');
}
void Flush() { Check(fflush(stdout)==0 && !ferror(stdout),"output stream"); }

void BindingsJson() {
    printf("[");
    for(unsigned a=0;a<2;++a) { if(a) putchar(','); const auto& spec=library_specs[a]; const auto& lib=*libraries[a];
        printf("{\"base\":%llu,\"context\":%llu,\"exports\":[",static_cast<unsigned long long>(lib.base),
            static_cast<unsigned long long>(lib.base+spec.context));
        bool comma=false;
        for(const auto& s:spec.exports) { if(comma) putchar(','); comma=true;
            printf("%llu",static_cast<unsigned long long>(Address(dlsym(lib.handle,s.name)))); }
        printf("],\"slots\":["); comma=false;
        for(const auto& s:spec.slots) { if(comma) putchar(','); comma=true; uintptr_t v=0;
            memcpy(&v,reinterpret_cast<const void*>(lib.base+s.offset),sizeof(v));
            printf("%llu",static_cast<unsigned long long>(v)); }
        printf("],\"providers\":["); comma=false;
        for(const char* name:{"malloc","free","memcpy","_Znwm","_ZdlPv"}) { if(comma) putchar(','); comma=true;
            printf("%llu",static_cast<unsigned long long>(Address(dlsym(lib.handle,name)))); }
        printf("],\"runtime_targets\":"); Array(lib.runtime_targets); printf("}");
    }
    putchar(']');
}
void HeaderJson(const std::string& claim,unsigned order,const std::string& identity,const Observation& prelude) {
    printf("{\"type\":\"header\",\"protocol\":\"%s\",\"claim\":\"%s\",\"load_order\":%u,\"identity_hex\":",protocol,claim.c_str(),order);
    Hex(identity.data(),identity.size()); printf(",\"prelude\":"); ObservationJson(prelude);
    printf(",\"bindings\":"); BindingsJson(); printf(",\"fixtures\":[");
    for(unsigned i=0;i<case_count;++i) { if(i) putchar(','); const auto& f=fixtures[i]; const auto& c=f.c;
        printf("{\"case\":[%u,%u,%u,%u],\"batch\":%u,\"source\":",c.family,c.k,c.b,c.policy,Batch(c));
        Hex(f.source.data(),f.source.size()); printf(",\"arms\":[");
        for(unsigned a=0;a<2;++a) { if(a) putchar(','); printf("{\"profile\":"); Hex(f.arm[a].profile.data(),32);
            printf(",\"packets\":"); Hex(f.arm[a].packets.data(),f.arm[a].packets.size());
            printf(",\"steps\":%u}",f.arm[a].steps); }
        printf("]}");
    }
    printf("]}\n"); Flush();
}
void RecordJson(const Record& r) {
        const Coordinate& c=r.c;
        printf("{\"type\":\"record\",\"coordinate\":[%u,%u,%u,%u,%u,%u,%u,%u,%llu],\"ready\":%llu,\"target\":%llu,\"wait\":",
            c.index,c.rep,c.order,c.which,c.metric,c.comparison,c.position,c.arm,static_cast<unsigned long long>(c.q),
            static_cast<unsigned long long>(r.ready),static_cast<unsigned long long>(r.target)); Array(r.wait);
        printf(",\"observation\":"); ObservationJson(r.o); printf(",\"counts\":"); Array(r.work.counts);
        printf(",\"addresses\":"); Array(r.work.addresses);
        printf(",\"address_count\":%u,\"complete\":%s,\"checked\":%s}\n",r.work.address_count,r.work.complete?"true":"false",r.checked?"true":"false");
    Flush();
}
void FooterJson(bool failed,uint64_t work) {
    printf("{\"type\":\"footer\",\"complete\":%s,\"records\":%u,\"work_ns\":%llu}\n",
        failed?"false":"true",retained,static_cast<unsigned long long>(work)); Flush();
}

void ValidateClaim(const std::string& claim,const char* path) {
    Check(claim.size()==64 && claim.find_first_not_of("0123456789abcdef")==std::string::npos,"claim hex");
    Check(wirehair::wh2_benchmark::Sha256Hex(ReadFile(path,1024*1024))==claim,"claim bytes");
}
int Worker(const std::string& claim,unsigned order) {
    Check(!WH2_ADMISSION_REGRESSION_NEUTRAL,"neutral scientific worker disabled");
    ValidateClaim(claim,claim_path);
    const rlimit cpu={180,180},memory={512u*1024u*1024u,512u*1024u*1024u},core={0,0};
    Check(!setrlimit(RLIMIT_CPU,&cpu) && !setrlimit(RLIMIT_AS,&memory) && !setrlimit(RLIMIT_CORE,&core),"worker limits");
    Reader reader; const uint64_t start=reader.Mono(),cpu_start=reader.Cpu();
    std::string identity; Observation prelude,previous; uint64_t work=0; bool failed=false,header=false; unsigned published=0;
    try {
        Pin(); Initialize(order); identity=Identity(); memset(outputs,0xa5,sizeof(outputs));
        for(unsigned i=0;i<callbacks;++i) records[i].c=CoordinateAt(i);
        uint64_t final=0; Capture(reader,[&]{final=Prelude(UINT64_C(0x9e3779b97f4a7c15),1u<<20);},prelude);
        Check(final==UINT64_C(0x43935dad1647741b),"prelude checksum"); Validate(prelude,previous); previous=prelude;
        HeaderJson(claim,order,identity,prelude); header=true;
        for(unsigned i=0;i<callbacks;++i) {
            Record& r=records[i]; retained=i+1; const auto& c=r.c; auto& f=fixtures[c.which]; const auto& arm=f.arm[c.arm]; Prepare(f);
            Check(reader.Mono()-start<UINT64_C(210000000000) && reader.Cpu()-cpu_start<UINT64_C(180000000000),"worker deadline");
            r.ready=reader.Mono(); Check(r.ready<=UINT64_MAX-c.q,"relative target overflow"); r.target=r.ready+c.q;
            r.wait[0]=reader.Mono(); r.wait[1]=reader.Cpu(); uint64_t now=r.wait[0];
            while(now<r.target) { now=reader.Mono(); Check(now-start<UINT64_C(210000000000),"wait deadline"); }
            r.wait[3]=reader.Cpu(); r.wait[2]=reader.Mono();
            Capture(reader,[&]{RunWork(f.api[c.arm],f,arm,c.metric,outputs[c.order],r.work);},r.o);
            Validate(r.o,previous); Check(r.o.m1>=r.target,"early start"); previous=r.o;
            work+=r.o.m2-r.o.m1; Check(work<=UINT64_C(150000000000),"WORK cap");
            CheckWork(f,arm,c.metric,c.order,r.work); r.checked=true; Check(OnCpu(),"affinity changed");
            RecordJson(r); ++published;
        }
        for(unsigned a=0;a<2;++a) libraries[a]->Validate(a);
        Check(Identity()==identity,"target changed");
        Check(reader.Mono()-start<UINT64_C(210000000000) && reader.Cpu()-cpu_start<UINT64_C(180000000000),"final deadline");
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); failed=true; }
    // A failed preflight has no valid fixture/binding header to serialize.
    if(!header) { FooterJson(true,work); return 1; }
    while(published<retained) RecordJson(records[published++]);
    FooterJson(failed,work); return failed?1:0;
}
struct FakeReader {
    uint64_t now=1000;
    uint64_t Mono() { return now+=100; }
    uint64_t Cpu() { return now+=10; }
    Counters Usage() { return {{0,0,0,0}}; }
};

const Symbols* injected=nullptr;
unsigned injected_calls=0,injected_at=0,injected_frees=0;
bool injected_throw=false;
WirehairV2Result InjectResult(WirehairV2Result result) {
    if(++injected_calls!=injected_at) return result;
    if(injected_throw) throw std::runtime_error("neutral last-call exception");
    return WirehairV2_Error;
}
WirehairV2Result InjectCreate(uint64_t id,const void* source,uint64_t m,uint32_t b,
    const WirehairV2EncoderOptions* options,void* p,uint32_t capacity,uint32_t* n,WirehairV2Codec* h) {
    return InjectResult(injected->wirehair_v2_encoder_create_profile_id_with_options(id,source,m,b,options,p,capacity,n,h));
}
WirehairV2Result InjectDecoder(const void* p,uint32_t n,WirehairV2Codec* h) {
    return InjectResult(injected->wirehair_v2_decoder_create(p,n,h));
}
WirehairV2Result InjectEncode(WirehairV2Codec h,uint32_t id,void* out,uint32_t b,uint32_t* n) {
    return InjectResult(injected->wirehair_v2_encode(h,id,out,b,n));
}
WirehairV2Result InjectFeed(WirehairV2Codec h,uint32_t id,const void* in,uint32_t n) {
    return InjectResult(injected->wirehair_v2_decode(h,id,in,n));
}
WirehairV2Result InjectRecover(WirehairV2Codec h,void* out,uint64_t n,uint64_t* written) {
    return InjectResult(injected->wirehair_v2_recover(h,out,n,written));
}
void InjectFree(WirehairV2Codec h) { ++injected_frees; injected->wirehair_v2_free(h); }
void FailureCheck() {
    // Exercise both actual DSOs and both fixed batch sizes through the same
    // WORK body, including a failing create that returned a live handle.
    for(unsigned i:{0u,8u}) for(unsigned a=0;a<2;++a) for(unsigned test=0;test<8;++test) {
        const unsigned stage=test>=6?test-6:test;
        auto& f=fixtures[i]; const auto& arm=f.arm[a]; const unsigned batch=Batch(f.c);
        const unsigned metric=stage==1 || stage==3 || stage==4;
        Symbols probe=symbols[a]; Api api=f.api[a]; api.s=&probe;
        injected=&symbols[a]; injected_calls=0; injected_frees=0; injected_throw=test>=5;
        injected_at=batch*(stage==2 || stage==5?PacketCount(f.c):stage==3?arm.steps:1);
        if(stage==0) probe.wirehair_v2_encoder_create_profile_id_with_options=InjectCreate;
        else if(stage==1) probe.wirehair_v2_decoder_create=InjectDecoder;
        else if(stage==2 || stage==5) probe.wirehair_v2_encode=InjectEncode;
        else if(stage==3) probe.wirehair_v2_decode=InjectFeed;
        else probe.wirehair_v2_recover=InjectRecover;
        probe.wirehair_v2_free=InjectFree;
        Prepare(f); Work w; Observation o; FakeReader reader;
        Capture(reader,[&]{RunWork(api,f,arm,metric,outputs[0],w);},o); Validate(o,Observation{});
        std::array<unsigned,6> counts={{metric?0:batch,metric?0:batch*PacketCount(f.c),
            metric?batch:0,metric?batch*arm.steps:0,metric?batch:0,batch}};
        if(stage==0) counts[1]-=PacketCount(f.c);
        if(stage==1) { counts[3]-=arm.steps; --counts[4]; }
        if(stage==3) --counts[4];
        Check(!w.complete && w.counts==counts && w.address_count==batch-(test>=6?1:0) && injected_calls==injected_at &&
              injected_frees==batch && o.m3>o.m2,"failed last-call capture/cleanup ledger");
    }
    injected=nullptr;
}
void RosterCheck() {
    unsigned seen[case_count][2][3][2][12]={},phases[case_count][2][3][2][48]={};
    for(unsigned i=0;i<callbacks;++i) {
        const auto c=CoordinateAt(i); Check(c.index==i && c.rep<12 && c.arm<2,"roster bounds");
        if(c.position==0) ++seen[c.which][c.metric][c.comparison][c.order][c.rep];
        if(c.position>=2 && c.position%2==0) {
            const auto n=CoordinateAt(i+1); Check(c.q==n.q && (sides[c.position]^sides[n.position])==1,"paired delays");
            ++phases[c.which][c.metric][c.comparison][c.order][unsigned((c.q*96/1000000)/2)];
        }
    }
    for(unsigned w=0;w<case_count;++w) for(unsigned m=0;m<2;++m) for(unsigned c=0;c<3;++c) for(unsigned o=0;o<2;++o) {
        for(unsigned r=0;r<12;++r) Check(seen[w][m][c][o][r]==1,"full roster");
        for(unsigned p=0;p<48;++p) Check(phases[w][m][c][o][p]==2,"all delay phases");
    }
}
int Neutral(unsigned order,bool fixtures_only) {
    RosterCheck(); Pin(); Initialize(order); const std::string identity=Identity();
    unsigned checked=0;
    if(!fixtures_only) for(unsigned i=0;i<case_count;++i) for(unsigned m=0;m<2;++m) for(unsigned a=0;a<2;++a) for(unsigned lane=0;lane<2;++lane) {
        auto& f=fixtures[i]; Prepare(f); Work w; Observation o; FakeReader reader;
        Capture(reader,[&]{RunWork(f.api[a],f,f.arm[a],m,outputs[lane],w);},o);
        Validate(o,Observation{}); CheckWork(f,f.arm[a],m,lane,w); ++checked;
    }
    if(!fixtures_only) FailureCheck();
    for(unsigned a=0;a<2;++a) libraries[a]->Validate(a);
    Check(Identity()==identity,"neutral target stability");
    if(fixtures_only) { Observation o; FakeReader r; Capture(r,[]{},o); HeaderJson(std::string(64,'0'),order,identity,o); }
    else printf("PASS neutral51840-coordinate roster, %u native WORK cases; no timing\n",checked);
    return 0;
}
unsigned Order(const char* s) { Check(!strcmp(s,"old-new") || !strcmp(s,"new-old"),"explicit load order"); return !strcmp(s,"new-old"); }
} // namespace
int main(int argc,char** argv) {
    try {
        for(const char* name:{"MALLOC_TRIM_THRESHOLD_","MALLOC_MMAP_THRESHOLD_","MALLOC_TOP_PAD_","MALLOC_PERTURB_",
                             "GLIBC_TUNABLES","LD_PRELOAD","LD_LIBRARY_PATH","LD_AUDIT","LD_DEBUG"})
            Check(!getenv(name),"clean allocator/loader environment");
        if(argc==2 && !strcmp(argv[1],"--binding")) {
            printf("%s\n%s\n",protocol,claim_path); Flush(); return 0;
        }
        if(argc==4 && !strcmp(argv[1],"--neutral-claim")) { ValidateClaim(argv[2],argv[3]); return 0; }
        if(argc==3 && !strcmp(argv[1],"--neutral")) return Neutral(Order(argv[2]),false);
        if(argc==3 && !strcmp(argv[1],"--neutral-fixtures")) return Neutral(Order(argv[2]),true);
        Check(argc==4 && !strcmp(argv[1],"--worker"),"explicit mode required");
        return Worker(argv[2],Order(argv[3]));
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); return 1; }
}
