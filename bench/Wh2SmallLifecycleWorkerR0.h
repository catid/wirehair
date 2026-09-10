// Future small-codec lifecycle worker. The separately frozen K2 worker is unchanged.
// Equation dimension, repair count and experiment identity come from a static wrapper.
// Both ownership policies, full/partial messages and each decoder's actual
// first-success cost are retained. Output occurs only after measurement ends.
#include "wirehair/wirehair.h"
#include "Wh2FrozenTrace.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "gf256.h"
#include "Wh2SmallSerialized.h"
#include <vector>
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <csignal>
#include <sched.h>
#include <sys/resource.h>
#include <time.h>

#if !defined(WH2_SMALL_COST_NEUTRAL) || defined(WH_COUNT) || defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "Explicit unmodified native/neutral build required"
#endif
#define NOINLINE __attribute__((noinline, noipa))
namespace {
const char protocol[] = WH2_SMALL_COST_PROTOCOL;
static_assert(WH2_SMALL_CODEC_K == WH2_SMALL_COST_K, "matching sealed C boundary required");
static_assert(WH2_SMALL_COST_K >= 2 && WH2_SMALL_COST_K <= 8, "supported fixed dimension");
const unsigned K=WH2_SMALL_COST_K, repair_count=WH2_SMALL_COST_REPAIRS;
const unsigned packet_count=K+2*repair_count, shape_count=6, batch=128, callbacks=77760, widths[3]={2,64,1280};
const unsigned cpu_seconds=240, wall_seconds=300, work_seconds=180, address_space_mib=384;
const char claim_path[]=WH2_SMALL_COST_CLAIM_PATH;
const unsigned pairs[10][2]={{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{0,1},{2,1},{3,4},{5,4}};
const unsigned sides[18]={0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
using Profile = std::array<uint8_t,32>;
void Check(bool value,const char* why) { if(!value) throw std::runtime_error(why); }
uint32_t Packet(unsigned slot) { return slot<K+repair_count?slot:UINT32_MAX-2*(slot-K-repair_count); }
unsigned Slot(unsigned family,unsigned step) { return step<repair_count?(family?K+repair_count:K)+step:step-repair_count; }
struct Api {
    int (*create)(const void*,unsigned,unsigned,Profile&,void*&);
    int (*encode)(void*,uint32_t,void*,unsigned,uint32_t*);
    int (*decoder)(const Profile&,unsigned,unsigned,void*&);
    int (*feed)(void*,uint32_t,const void*,unsigned);
    int (*recover)(void*,void*,unsigned,uint64_t*);
    void (*free)(void*);
};
int PCreate(const void* source,unsigned message,unsigned b,Profile& p,void*& h,uint32_t policy) {
    WirehairV2Codec handle=nullptr; uint32_t n=0;
    WirehairV2EncoderOptions o=WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    o.source_policy=policy;
    const int r=policy==WirehairV2EncoderSource_Independent ?
        wirehair_v2_encoder_create(source,message,b,p.data(),32,&n,&handle) :
        wirehair_v2_encoder_create_with_options(source,message,b,&o,p.data(),32,&n,&handle);
    h=handle; return r==0 && n!=32?-1:r;
}
int PEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    return wirehair_v2_encode(static_cast<WirehairV2Codec>(h),id,out,b,n);
}
int PDecoder(const Profile& p,unsigned,unsigned,void*& h) {
    WirehairV2Codec handle=nullptr;
    const int r=wirehair_v2_decoder_create(p.data(),32,&handle); h=handle; return r;
}
int PFeed(void* h,uint32_t id,const void* in,unsigned n) {
    return wirehair_v2_decode(static_cast<WirehairV2Codec>(h),id,in,n);
}
int PRecover(void* h,void* out,unsigned n,uint64_t* written) {
    return wirehair_v2_recover(static_cast<WirehairV2Codec>(h),out,n,written);
}
void PFree(void* h) { wirehair_v2_free(static_cast<WirehairV2Codec>(h)); }
int PIndependent(const void* s,unsigned message,unsigned b,Profile& p,void*& h) {
    return PCreate(s,message,b,p,h,WirehairV2EncoderSource_Independent);
}
int PBorrowed(const void* s,unsigned message,unsigned b,Profile& p,void*& h) {
    return PCreate(s,message,b,p,h,WirehairV2EncoderSource_BorrowedImmutable);
}
int SCreate(const void* source,unsigned message,unsigned b,Profile& p,void*& h,uint32_t policy) {
    const auto r=wh2_small_encoder_create(source,message,b,policy,p.data(),p.size()); h=r.codec; return r.status;
}
int SIndependent(const void* s,unsigned message,unsigned b,Profile& p,void*& h) { return SCreate(s,message,b,p,h,Wh2Small_Independent); }
int SBorrowed(const void* s,unsigned message,unsigned b,Profile& p,void*& h) { return SCreate(s,message,b,p,h,Wh2Small_BorrowedImmutable); }
int SEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const auto r=wh2_small_encode(h,id,out,b); *n=static_cast<uint32_t>(r.bytes_written);
    return r.status==Wh2Small_Success && (r.bytes_written!=r.bytes_required || !r.bytes_written || r.bytes_written>b)?-1:static_cast<int>(r.status);
}
int SDecoder(const Profile& p,unsigned,unsigned,void*& h) {
    const auto r=wh2_small_decoder_create(p.data(),p.size()); h=r.codec; return r.status;
}
int SFeed(void* h,uint32_t id,const void* in,unsigned n) { return wh2_small_decode(h,id,in,n); }
int SRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const auto r=wh2_small_recover(h,out,n); *written=r.bytes_written;
    return r.status==Wh2Small_Success && (r.bytes_required!=n || r.bytes_written!=n)?-1:static_cast<int>(r.status);
}
void SFree(void* h) { wh2_small_free(h); }
int WCreate(const void* source,unsigned message,unsigned b,Profile&,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_encoder_create_ex(nullptr,source,message,b,&handle); h=handle; return r;
}
int WOwned(const void* source,unsigned message,unsigned b,Profile&,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_encoder_create_owned_ex(nullptr,source,message,b,&handle); h=handle; return r;
}
int WEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    return wirehair_encode(static_cast<WirehairCodec>(h),id,out,b,n);
}
int WDecoder(const Profile&,unsigned message,unsigned b,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_decoder_create_ex(nullptr,message,b,&handle); h=handle; return r;
}
int WFeed(void* h,uint32_t id,const void* in,unsigned n) { return wirehair_decode(static_cast<WirehairCodec>(h),id,in,n); }
int WRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const int r=wirehair_recover(static_cast<WirehairCodec>(h),out,n); *written=r==0?n:0; return r;
}
void WFree(void* h) { wirehair_free(static_cast<WirehairCodec>(h)); }
const Api apis[6]={
    {PIndependent,PEncode,PDecoder,PFeed,PRecover,PFree},
    {SIndependent,SEncode,SDecoder,SFeed,SRecover,SFree},
    {WOwned,WEncode,WDecoder,WFeed,WRecover,WFree},
    {PBorrowed,PEncode,PDecoder,PFeed,PRecover,PFree},
    {SBorrowed,SEncode,SDecoder,SFeed,SRecover,SFree},
    {WCreate,WEncode,WDecoder,WFeed,WRecover,WFree}};
struct Owner {
    const Api& api; void* handle=nullptr;
    explicit Owner(const Api& a):api(a) {}
    ~Owner() { if(handle) api.free(handle); }
    Owner(const Owner&)=delete; Owner& operator=(const Owner&)=delete;
};
struct Arm { Profile profile={}; uint8_t packets[packet_count*1280]={}, rows[packet_count*K]={}; unsigned steps[2]={}; };
struct Fixture { unsigned b=0,message=0; uint8_t source[(K*1280)]={}; Arm arm[6]; };
Fixture fixtures[shape_count], reference[shape_count];
unsigned PacketBytes(const Fixture& f,unsigned slot) { return Packet(slot)==K-1?f.message-(K-1)*f.b:f.b; }
alignas(64) uint8_t outputs[2][batch*(packet_count*1280+128)];
struct Coordinate { unsigned index,rep,order,width,metric,comparison,position,arm; uint64_t q; };
Coordinate CoordinateAt(unsigned index) {
    unsigned n=index; const unsigned p=n%18; n/=18; const unsigned cs=n%10; n/=10;
    const unsigned ms=n%3; n/=3; const unsigned ws=n%shape_count; n/=shape_count; const unsigned s=n%2,r=n/2;
    const unsigned order=(r+s)%2,w=(r+s+ws)%shape_count,m=(r+s+ws+ms)%3,c=(2*r+s+ws+m+cs)%10;
    const unsigned bin=p<2?r+12*(r%4):((r+6*((p-2)/8))%12)+12*(((p-2)%8)/2);
    return Coordinate{index,r,order,w,m,c,p,pairs[c][sides[p]^order],uint64_t(2*bin+1)*1000000/96};
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
    // create, encode, decoder-create, feed, recover, free; every attempted call.
    std::array<unsigned,6> counts={}; std::array<uint64_t,batch> addresses={};
    unsigned address_count=0; bool complete=false;
};
struct Record {
    Coordinate c={}; uint64_t ready=0,target=0; std::array<uint64_t,4> wait={};
    Observation o; Work work; bool checked=false;
};
Record records[callbacks]; unsigned retained=0;
size_t Stride(const Fixture& f) { return size_t(packet_count)*f.b+128; }
NOINLINE void RunWork(const Api& api,const Fixture& f,const Arm& arm,unsigned metric,uint8_t* output,Work& w) noexcept {
    void* handle=nullptr;
    auto release=[&] { if(handle) { void* owned=handle; handle=nullptr; ++w.counts[5]; api.free(owned); } };
    try {
        for(unsigned cycle=0;cycle<batch;++cycle) {
            uint8_t* out=output+cycle*Stride(f)+64;
            Profile p={}; int result;
            if(metric==0) { ++w.counts[0]; result=api.create(f.source,f.message,f.b,p,handle); }
            else { ++w.counts[2]; result=api.decoder(arm.profile,f.message,f.b,handle); p=arm.profile; }
            w.addresses[w.address_count++]=uint64_t(reinterpret_cast<uintptr_t>(handle));
            if(result!=0 || !handle || p!=arm.profile) { release(); return; }
            if(metric==0) {
                for(unsigned j=0;j<packet_count;++j) {
                    uint32_t written=0; ++w.counts[1]; result=api.encode(handle,Packet(j),out+j*f.b,f.b,&written);
                    if(result!=0 || written!=PacketBytes(f,j)) { release(); return; }
                }
            } else {
                const unsigned family=metric-1,steps=arm.steps[family];
                bool success=false;
                for(unsigned j=0;j<repair_count+K;++j) {
                    const unsigned slot=Slot(family,j);
                    ++w.counts[3]; result=api.feed(handle,Packet(slot),arm.packets+slot*f.b,PacketBytes(f,slot));
                    if(j>=steps || result!=(j+1==steps?0:1)) { release(); return; }
                    if(result==0) { success=true; break; }
                }
                if(!success) { release(); return; }
                uint64_t written=0; ++w.counts[4]; result=api.recover(handle,out,f.message,&written);
                if(result!=0 || written!=f.message) { release(); return; }
            }
            release();
        }
        w.complete=true;
    } catch(...) { try { release(); } catch(...) {} }
}
void Initialize() {
    Check(wirehair_init()==Wirehair_Success,"shared GF init");
    for(unsigned wi=0;wi<shape_count;++wi) {
        Fixture& f=fixtures[wi]; f.b=widths[wi/2]; f.message=(K-1)*f.b+(wi%2?1:f.b);
        for(size_t i=0;i<sizeof(f.source);++i) f.source[i]=uint8_t(37*i+i/11);
        for(unsigned a=0;a<6;++a) {
            const Api& api=apis[a]; Arm& arm=f.arm[a];
            std::fill(arm.packets,arm.packets+sizeof(arm.packets),0xa5);
            { std::vector<uint8_t> input(f.source,f.source+f.message); Owner e(api);
              Check(api.create(input.data(),f.message,f.b,arm.profile,e.handle)==0 && e.handle,"fixture create");
              if(a<3) { std::fill(input.begin(),input.end(),0xcc); std::vector<uint8_t>().swap(input); }
              for(unsigned j=0;j<packet_count;++j) { uint32_t n=0;
                Check(api.encode(e.handle,Packet(j),arm.packets+j*f.b,f.b,&n)==0 && n==PacketBytes(f,j),"fixture encode");
                if(a>=3) Check(input.size()==f.message && !memcmp(input.data(),f.source,f.message),"actual borrowed source immutable");
                for(unsigned byte=n;byte<f.b;++byte) Check(arm.packets[j*f.b+byte]==0xa5,"fixture unused packet capacity"); } }
            // Observe native generator rows independently of message payloads.
            // Every basis encoder uses this exact B and is gone before decoding.
            for(unsigned column=0;column<K;++column) {
                std::vector<uint8_t> input(f.message,0); input[column*f.b]=1;
                const std::vector<uint8_t> original=input;
                Owner e(api); Profile p={};
                Check(api.create(input.data(),f.message,f.b,p,e.handle)==0 && e.handle && p==arm.profile,"basis encoder");
                if(a<3) { std::fill(input.begin(),input.end(),0xcc); std::vector<uint8_t>().swap(input); }
                for(unsigned j=0;j<packet_count;++j) {
                    std::array<uint8_t,1282> packet; packet.fill(0xa5); uint32_t n=0;
                    Check(api.encode(e.handle,Packet(j),packet.data()+1,f.b,&n)==0 && n==PacketBytes(f,j),"basis packet");
                    if(a>=3) Check(input==original,"actual borrowed basis immutable");
                    Check(packet[0]==0xa5,"basis prefix guard");
                    for(unsigned byte=1;byte<n;++byte) Check(packet[byte+1]==0,"basis other bytes");
                    for(size_t byte=n+1;byte<packet.size();++byte) Check(packet[byte]==0xa5,"basis tail/capacity guard");
                    arm.rows[j*K+column]=packet[1];
                }
            }
        }
        // No sender remains live when any of the six receivers is created.
        for(unsigned a=0;a<6;++a) {
            const Api& api=apis[a]; Arm& arm=f.arm[a];
            for(unsigned family=0;family<2;++family) {
                Owner d(api); Check(api.decoder(arm.profile,f.message,f.b,d.handle)==0 && d.handle,"fixture decoder");
                for(unsigned j=0;j<repair_count+K;++j) { const unsigned slot=Slot(family,j);
                    const int r=api.feed(d.handle,Packet(slot),arm.packets+slot*f.b,PacketBytes(f,slot));
                    Check(r==0 || r==1,"fixture feed");
                    if(r==0) { Check(j>=K-1,"premature decode"); arm.steps[family]=j+1; break; }
                }
                Check(arm.steps[family]>0,"fixture endpoint");
                for(unsigned repeat=0;repeat<2;++repeat) {
                    std::array<uint8_t,K*1280+2> out; out.fill(uint8_t(0xa5+repeat)); uint64_t n=0;
                    Check(api.recover(d.handle,out.data()+1,f.message,&n)==0 && n==f.message &&
                          !memcmp(out.data()+1,f.source,f.message),"fixture repeated recovery");
                    Check(out[0]==uint8_t(0xa5+repeat),"fixture recovery prefix guard");
                    for(size_t byte=f.message+1;byte<out.size();++byte)
                        Check(out[byte]==uint8_t(0xa5+repeat),"fixture recovery tail guard");
                }
            }
        }
        for(unsigned a : {0u,3u}) {
            WirehairV2Profile profile={};
            Check(wirehair_v2_profile_deserialize(f.arm[a].profile.data(),32,&profile)==WirehairV2_Success &&
                  profile.profile_id==WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,"ordinary control still selects certified equations");
        }
        for(unsigned a : {1u,4u})
            Check(wh2_small_profile_validate(f.arm[a].profile.data(),32)==Wh2Small_Success,"explicit sealed small-codec identity");
        for(unsigned a=0;a<3;++a)
            Check(f.arm[a].profile==f.arm[a+3].profile &&
                  !memcmp(f.arm[a].packets,f.arm[a+3].packets,packet_count*f.b) &&
                  !memcmp(f.arm[a].rows,f.arm[a+3].rows,packet_count*K) &&
                  std::equal(f.arm[a].steps,f.arm[a].steps+2,f.arm[a+3].steps),"policy-independent equations");
        reference[wi]=f;
    }
}
void Prepare(const Fixture& f) { for(auto& out:outputs) memset(out,0xa5,batch*Stride(f)); }
void CheckWork(const Fixture& f,const Arm& arm,unsigned metric,unsigned lane,const Work& w) {
    const unsigned steps=metric?arm.steps[metric-1]:0;
    const std::array<unsigned,6> expected={{metric?0:batch,metric?0:batch*packet_count,metric?batch:0,batch*steps,metric?batch:0,batch}};
    Check(w.complete && w.counts==expected && w.address_count==batch,"full work ledger");
    for(uint64_t address:w.addresses) Check(address!=0,"live addresses");
    const size_t used=metric?f.message:packet_count*f.b,stride=Stride(f);
    for(unsigned cycle=0;cycle<batch;++cycle) {
        const uint8_t* out=outputs[lane]+cycle*stride;
        Check(!memcmp(out+64,metric?f.source:arm.packets,used),"cycle bytes");
        for(size_t j=0;j<stride;++j) if(j<64 || j>=64+used) Check(out[j]==0xa5,"cycle guard");
    }
    for(size_t j=0;j<batch*stride;++j) Check(outputs[1-lane][j]==0xa5,"unused lane");
    for(unsigned i=0;i<shape_count;++i) {
        Check(fixtures[i].message==reference[i].message && fixtures[i].b==reference[i].b && !memcmp(fixtures[i].source,reference[i].source,(K*1280)),"immutable source");
        for(unsigned a=0;a<6;++a) Check(fixtures[i].arm[a].profile==reference[i].arm[a].profile &&
            !memcmp(fixtures[i].arm[a].packets,reference[i].arm[a].packets,packet_count*1280) &&
            !memcmp(fixtures[i].arm[a].rows,reference[i].arm[a].rows,packet_count*K) &&
            std::equal(fixtures[i].arm[a].steps,fixtures[i].arm[a].steps+2,reference[i].arm[a].steps),"immutable corpus");
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
    const auto& d=r.Derived;
    Check(d.Family==26 && d.Model==8 && d.Stepping==1 && d.FullApicId==100 && d.CoreId==50 &&
          d.PackageId==0 && d.ThreadId==0 && d.ThreadsPerCore==2 && d.CcdId==6 && d.ComplexId==6 &&
          d.LogicalProcessorsPerPackage==128,"frozen physical target");
    gf256_x86_cpu_features f={}; gf256_get_active_x86_cpu_features(&f);
    Check(GF256Ctx.Polynomial==0x14d && f.SSSE3==1 && f.AVX2==1 && f.GFNI==1 && f.AVX512==1,"frozen runtime dispatch");
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
void HeaderJson(const std::string& claim,const std::string& identity,const Observation& prelude) {
    printf("{\"type\":\"header\",\"protocol\":\"%s\",\"claim\":\"%s\",\"batch\":%u,\"identity_hex\":",protocol,claim.c_str(),batch);
    Hex(identity.data(),identity.size()); printf(",\"prelude\":"); ObservationJson(prelude); printf(",\"fixtures\":[");
    // Keep the initial corpus even if a later WORK failure corrupts fixtures.
    for(unsigned wi=0;wi<shape_count;++wi) { const Fixture& f=reference[wi]; if(wi) putchar(',');
        printf("{\"width\":%u,\"message\":%u,\"packet_bytes\":[",f.b,f.message);
        for(unsigned slot=0;slot<packet_count;++slot) { if(slot) putchar(','); printf("%u",PacketBytes(f,slot)); }
        printf("],\"source\":"); Hex(f.source,f.message); printf(",\"arms\":[");
        for(unsigned a=0;a<6;++a) { if(a) putchar(','); printf("{\"profile\":"); Hex(f.arm[a].profile.data(),32);
            printf(",\"packets\":"); Hex(f.arm[a].packets,packet_count*f.b);
            printf(",\"rows\":"); Hex(f.arm[a].rows,packet_count*K);
            printf(",\"steps\":[%u,%u],\"fallback\":[%s,%s]}",f.arm[a].steps[0],f.arm[a].steps[1],
                f.arm[a].steps[0]>repair_count?"true":"false",f.arm[a].steps[1]>repair_count?"true":"false"); }
        printf("]}"); }
    printf("]}\n"); Flush();
}
void RecordJson(const Record& r) {
        const Coordinate& c=r.c;
        printf("{\"type\":\"record\",\"coordinate\":[%u,%u,%u,%u,%u,%u,%u,%u,%llu],\"ready\":%llu,\"target\":%llu,\"wait\":",
            c.index,c.rep,c.order,c.width,c.metric,c.comparison,c.position,c.arm,static_cast<unsigned long long>(c.q),
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
struct JsonSink {
    const std::string& claim; const std::string& identity; const Observation& prelude;
    void Header() { HeaderJson(claim,identity,prelude); }
    void Row(const Record& record) { RecordJson(record); }
    void Footer(bool failed,uint64_t work) { FooterJson(failed,work); }
};
template<class Sink> void Publish(Sink& sink,bool measurement_finished,bool failed,uint64_t work) {
    Check(measurement_finished && retained<=callbacks,"publication before measurement end");
    // Single pass: a failing writer must not retry a partially emitted header
    // or record. In-memory capture remains intact; the controller retains the
    // exact emitted prefix and rejects missing records/footer or worker errors.
    sink.Header();
    for(unsigned i=0;i<retained;++i) sink.Row(records[i]);
    sink.Footer(failed,work);
}
void Authenticate(const char* path,const std::string& claim) {
    Check(claim.size()==64 && claim.find_first_not_of("0123456789abcdef")==std::string::npos,"claim hex");
    std::ifstream in(path,std::ios::binary);
    Check(bool(in),"claimed namespace"); std::string receipt; char ch;
    while(in.get(ch)) { Check(receipt.size()<1024*1024,"claim cap"); receipt+=ch; }
    Check(in.eof() && wirehair::wh2_benchmark::Sha256Hex(receipt)==claim,"claim bytes");
}
int Worker(const std::string& claim) {
    Check(!WH2_SMALL_COST_NEUTRAL,"neutral scientific worker disabled");
    Authenticate(claim_path,claim);
    const rlimit cpu={cpu_seconds,cpu_seconds},
        memory={address_space_mib*1024u*1024u,address_space_mib*1024u*1024u},core={0,0};
    Check(!setrlimit(RLIMIT_CPU,&cpu) && !setrlimit(RLIMIT_AS,&memory) && !setrlimit(RLIMIT_CORE,&core),"worker limits");
    Reader reader; const uint64_t start=reader.Mono(),cpu_start=reader.Cpu();
    std::string identity; Observation prelude,previous; uint64_t work=0; const char* failure=nullptr;
    try {
        Pin(); Initialize(); identity=Identity();
        memset(outputs,0xa5,sizeof(outputs));
        for(unsigned i=0;i<callbacks;++i) records[i].c=CoordinateAt(i);
        uint64_t final=0; Capture(reader,[&] { final=Prelude(UINT64_C(0x9e3779b97f4a7c15),1u<<20); },prelude);
        Check(final==UINT64_C(0x43935dad1647741b),"prelude checksum"); Validate(prelude,previous); previous=prelude;
        for(unsigned i=0;i<callbacks;++i) {
            Record& r=records[i]; retained=i+1; const Coordinate& c=r.c;
            Fixture& f=fixtures[c.width]; const Arm& arm=f.arm[c.arm]; Prepare(f);
            Check(reader.Mono()-start<uint64_t(wall_seconds)*1000000000 &&
                  reader.Cpu()-cpu_start<uint64_t(cpu_seconds)*1000000000,"worker deadline");
            r.ready=reader.Mono(); Check(r.ready<=UINT64_MAX-c.q,"relative target overflow"); r.target=r.ready+c.q;
            r.wait[0]=reader.Mono(); r.wait[1]=reader.Cpu();
            uint64_t now=r.wait[0];
            while(now<r.target) { now=reader.Mono(); Check(now-start<uint64_t(wall_seconds)*1000000000,"wait deadline"); }
            r.wait[3]=reader.Cpu(); r.wait[2]=reader.Mono();
            Capture(reader,[&] { RunWork(apis[c.arm],f,arm,c.metric,outputs[c.order],r.work); },r.o);
            Validate(r.o,previous); Check(r.o.m1>=r.target,"early start"); previous=r.o;
            work+=r.o.m2-r.o.m1; Check(work<=uint64_t(work_seconds)*1000000000,"inner work cap");
            CheckWork(f,arm,c.metric,c.order,r.work); r.checked=true; Check(OnCpu(),"affinity changed");
        }
        Check(Identity()==identity,"target changed");
        Check(reader.Mono()-start<uint64_t(wall_seconds)*1000000000 &&
              reader.Cpu()-cpu_start<uint64_t(cpu_seconds)*1000000000,"final worker deadline");
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); failure="invalid"; }
      catch(...) { fprintf(stderr,"INVALID: unknown measurement exception\n"); failure="invalid"; }
    // This is the only publication site; no output-related blocking can occur
    // between measured callbacks. Broken pipes become explicit output errors.
    Check(std::signal(SIGPIPE,SIG_IGN)!=SIG_ERR,"output signal policy");
    JsonSink sink{claim,identity,prelude}; Publish(sink,true,failure!=nullptr,work);
    Check(reader.Mono()-start<uint64_t(wall_seconds)*1000000000 &&
          reader.Cpu()-cpu_start<uint64_t(cpu_seconds)*1000000000,"publication deadline");
    return failure?1:0;
}
struct FakeReader {
    uint64_t now=1000;
    uint64_t Mono() { return now+=100; }
    uint64_t Cpu() { return now+=10; }
    Counters Usage() { return {{0,0,0,0}}; }
};
const Api* injected=nullptr; unsigned injected_calls=0,injected_at=0;
int BadEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const int r=injected->encode(h,id,out,b,n); if(++injected_calls==injected_at) ++*n; return r;
}
int BadFeed(void* h,uint32_t id,const void* in,unsigned n) {
    const int r=injected->feed(h,id,in,n); return ++injected_calls==injected_at?-1:r;
}
int BadRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const int r=injected->recover(h,out,n,written); if(++injected_calls==injected_at) ++*written; return r;
}
int ThrowEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const int r=injected->encode(h,id,out,b,n);
    if(++injected_calls==injected_at) throw std::runtime_error("neutral injected exception");
    return r;
}
int ThrowRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const int r=injected->recover(h,out,n,written);
    if(++injected_calls==injected_at) throw std::runtime_error("neutral injected exception");
    return r;
}
int ThrowFeed(void* h,uint32_t id,const void* in,unsigned n) {
    const int r=injected->feed(h,id,in,n);
    if(++injected_calls==injected_at) throw std::runtime_error("neutral injected feed exception");
    return r;
}
int BadTailCapacity(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const int r=injected->encode(h,id,out,b,n);
    if(++injected_calls==injected_at) {
        Check(id==K-1 && *n<b,"partial systematic fault reached");
        static_cast<uint8_t*>(out)[*n]^=1;
    }
    return r;
}
struct LateClockReader : FakeReader {
    unsigned monos=0;
    uint64_t Mono() {
        Check(++monos!=4,"neutral final clock failure");
        return FakeReader::Mono();
    }
};
int NeutralPublication(const char* mode) {
    const bool success=!strcmp(mode,"success"),bad=!strcmp(mode,"last-recover"),
        throwing=!strcmp(mode,"throw-recover"),clock=!strcmp(mode,"last-clock"),
        corrupt=!strcmp(mode,"last-source");
    Check(success || bad || throwing || clock || corrupt,"neutral publication mode");
    Initialize(); Observation prelude,previous; FakeReader reader;
    Capture(reader,[]{},prelude); previous=prelude; uint64_t work=0; bool failed=false;
    for(unsigned i=0;i<3;++i) {
        Record& r=records[i]; r.c=CoordinateAt(180*i); retained=i+1;
        const Coordinate& c=r.c; Fixture& f=fixtures[c.width]; const Arm& arm=f.arm[c.arm];
        Prepare(f); Api api=apis[c.arm];
        injected=&apis[c.arm]; injected_calls=0; injected_at=batch;
        if(i==2 && bad) api.recover=BadRecover;
        if(i==2 && throwing) api.recover=ThrowRecover;
        try {
            if(i==2 && clock) {
                LateClockReader late; late.now=reader.now;
                Capture(late,[&] { RunWork(api,f,arm,c.metric,outputs[c.order],r.work); },r.o);
            } else Capture(reader,[&] { RunWork(api,f,arm,c.metric,outputs[c.order],r.work); },r.o);
            Validate(r.o,previous); previous=r.o; work+=r.o.m2-r.o.m1;
            if(i==2 && corrupt) fixtures[0].source[0]^=1;
            CheckWork(f,arm,c.metric,c.order,r.work); r.checked=true;
        } catch(const std::exception&) { failed=true; break; }
    }
    Check(failed!=success && retained==3 && records[0].checked && records[1].checked,
          "neutral publication capture");
    const Work& last=records[2].work;
    Check(last.counts[2]==batch && last.counts[4]==batch && last.counts[5]==batch &&
          last.address_count==batch && records[2].o.m2>records[2].o.m1,
          "neutral final lifecycle retention");
    if(clock) Check(last.complete && records[2].o.m3==0 && !records[2].checked,"partial clock retention");
    else Check(last.complete==(success || corrupt) && records[2].checked==success && records[2].o.m3>records[2].o.m2,
               "neutral result retention");
    const std::string claim(64,'0'),identity="neutral-deferred";
    JsonSink sink{claim,identity,prelude}; Publish(sink,true,failed,work);
    return 0; // Neutral injected failures are expected, never science evidence.
}
struct TestSink {
    unsigned calls=0,fail_at=0; bool failed=false; uint64_t work=0;
    void Step() { Check(++calls!=fail_at,"neutral output failure"); }
    void Header() { Step(); }
    void Row(const Record& record) { Check(&record==&records[calls-1],"publication order"); Step(); }
    void Footer(bool f,uint64_t w) { Step(); failed=f; work=w; }
};
void NeutralWriter() {
    retained=3;
    for(unsigned failure=0;failure<=5;++failure) {
        TestSink sink; sink.fail_at=failure; bool caught=false;
        try { Publish(sink,true,true,123); } catch(const std::exception&) { caught=true; }
        Check(caught==(failure!=0) && sink.calls==(failure?failure:5),"no publication retry");
        if(!failure) Check(sink.failed && sink.work==123,"failed footer retained");
    }
    TestSink before; bool caught=false;
    try { Publish(before,false,false,0); } catch(const std::exception&) { caught=true; }
    Check(caught && before.calls==0,"no early publication");
    retained=0; TestSink empty; Publish(empty,true,true,0);
    Check(empty.calls==2 && empty.failed,"failure before first record");
    retained=callbacks+1; TestSink oversized; caught=false;
    try { Publish(oversized,true,false,0); } catch(const std::exception&) { caught=true; }
    Check(caught && oversized.calls==0,"record bound before publication"); retained=0;
}
int Neutral() {
    NeutralWriter();
    unsigned seen[shape_count][3][10][2][12]={},phases[shape_count][3][10][2][48]={};
    for(unsigned i=0;i<callbacks;++i) {
        const auto c=CoordinateAt(i); Check(c.index==i && c.rep<12 && c.arm<6,"roster bounds");
        if(c.position==0) ++seen[c.width][c.metric][c.comparison][c.order][c.rep];
        if(c.position>=2 && c.position%2==0) {
            const auto next=CoordinateAt(i+1); Check(c.q==next.q && (sides[c.position]^sides[next.position])==1,"adjacent pair");
            const unsigned bin=unsigned((c.q*96/1000000)/2);
            ++phases[c.width][c.metric][c.comparison][c.order][bin];
        }
    }
    for(unsigned w=0;w<shape_count;++w) for(unsigned m=0;m<3;++m) for(unsigned c=0;c<10;++c) for(unsigned o=0;o<2;++o) {
        for(unsigned r=0;r<12;++r) Check(seen[w][m][c][o][r]==1,"complete roster");
        for(unsigned p=0;p<48;++p) Check(phases[w][m][c][o][p]==2,"phase coverage"); }
    Initialize(); unsigned checked=0;
    for(unsigned wi=0;wi<shape_count;++wi) for(unsigned m=0;m<3;++m) for(unsigned a=0;a<6;++a) for(unsigned lane=0;lane<2;++lane) {
        auto& f=fixtures[wi]; Prepare(f); Work w; Observation o; FakeReader reader;
        Capture(reader,[&] { RunWork(apis[a],f,f.arm[a],m,outputs[lane],w); },o);
        Validate(o,Observation{}); CheckWork(f,f.arm[a],m,lane,w); ++checked;
    }
    auto& f=fixtures[0]; Prepare(f); Work w; Observation o; FakeReader reader; Api bad=apis[0];
    injected=&apis[0]; injected_calls=0; injected_at=batch*packet_count; bad.encode=BadEncode;
    Capture(reader,[&] { RunWork(bad,f,f.arm[0],0,outputs[0],w); },o); Validate(o,Observation{});
    Check(!w.complete && w.counts[0]==batch && w.counts[1]==batch*packet_count && w.counts[5]==batch && o.m3>o.m2,"failed last-call retention");
    { Prepare(f); Work failed; Observation observation; FakeReader fake; Api probe=apis[0];
      injected_calls=0; injected_at=batch*packet_count; probe.encode=ThrowEncode;
      Capture(fake,[&] { RunWork(probe,f,f.arm[0],0,outputs[0],failed); },observation);
      Validate(observation,Observation{});
      Check(!failed.complete && failed.counts[0]==batch && failed.counts[1]==batch*packet_count &&
            failed.counts[5]==batch && observation.m3>observation.m2,"throwing last-call retention"); }
    for(unsigned failure=0;failure<4;++failure) {
        Prepare(f); Work failed; Observation observation; FakeReader fake; Api probe=apis[0];
        const bool recovery=failure%2!=0,throwing=failure>=2;
        injected_calls=0; injected_at=recovery?batch:batch*f.arm[0].steps[1];
        if(recovery) probe.recover=throwing?ThrowRecover:BadRecover;
        else probe.feed=throwing?ThrowFeed:BadFeed;
        Capture(fake,[&] { RunWork(probe,f,f.arm[0],2,outputs[0],failed); },observation);
        Validate(observation,Observation{});
        Check(!failed.complete && failed.counts[2]==batch && failed.counts[5]==batch &&
              failed.counts[4]==batch-(recovery?0:1) && injected_calls==injected_at &&
              observation.m3>observation.m2,"decoder failed-last-call retention");
    }
    for(unsigned a=0;a<6;++a) {
        auto& partial=fixtures[1]; Prepare(partial); Work wtail; Api probe=apis[a];
        injected=&apis[a]; injected_calls=0; injected_at=(batch-1)*packet_count+K;
        probe.encode=BadTailCapacity;
        RunWork(probe,partial,partial.arm[a],0,outputs[0],wtail);
        Check(wtail.complete && injected_calls==batch*packet_count && wtail.counts[5]==batch,
              "tail fault final encoder completely freed");
        bool rejected=false;
        try { CheckWork(partial,partial.arm[a],0,0,wtail); }
        catch(const std::exception&) { rejected=true; }
        Check(rejected,"unused partial packet capacity corruption detected");
    }
    printf("PASS neutral77760-coordinate roster,216 actual WORK cases, last-call cleanup/capture (%u checked)\n",checked); return 0;
}
} // namespace
int main(int argc,char** argv) {
    try {
        if(argc==2 && !strcmp(argv[1],"--contract")) {
            printf("{\"K\":%u,\"batch\":%u,\"callbacks\":%u,\"cpu_seconds\":%u,\"wall_seconds\":%u,\"work_seconds\":%u,\"address_space_mib\":%u,\"claim_path\":\"%s\",\"neutral_only\":%s,\"gf_context_bytes\":%zu}\n",
                K,batch,callbacks,cpu_seconds,wall_seconds,work_seconds,address_space_mib,claim_path,
                WH2_SMALL_COST_NEUTRAL?"true":"false",sizeof(GF256Ctx)); return 0;
        }
        if(argc==2 && !strcmp(argv[1],"--claim-path")) { puts(claim_path); return 0; }
        if(argc==4 && !strcmp(argv[1],"--neutral-claim")) {
            Authenticate(argv[2],argv[3]); puts("PASS claim authentication"); return 0;
        }
        if(argc==2 && !strcmp(argv[1],"--neutral")) return Neutral();
        if(argc==3 && !strcmp(argv[1],"--neutral-publication")) {
            Check(std::signal(SIGPIPE,SIG_IGN)!=SIG_ERR,"output signal policy");
            return NeutralPublication(argv[2]);
        }
        if(argc==2 && !strcmp(argv[1],"--neutral-target")) {
            Check(!WH2_SMALL_COST_NEUTRAL,"native target capture only");
            Pin(); Initialize(); Observation o; FakeReader r; Capture(r,[]{},o);
            // Identity is binary (embedded NULs), not a C string. Reuse the
            // exact length-aware header transport used by the scientific run.
            HeaderJson(std::string(64,'0'),Identity(),o); return 0;
        }
        if(argc==2 && !strcmp(argv[1],"--neutral-fixtures")) {
            Initialize(); Observation o; FakeReader r; Capture(r,[]{},o);
            HeaderJson(std::string(64,'0'),"neutral",o); return 0;
        }
        Check(argc==3 && !strcmp(argv[1],"--worker"),"explicit mode required"); return Worker(argv[2]);
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); return 1; }
}
