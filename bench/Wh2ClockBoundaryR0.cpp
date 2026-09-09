// Codec-free, one-shot timing diagnostic. Never a codec speed gate.
// Ordered TSC rationale: Linux v6.8 arch/x86/include/asm/msr.h,
// https://github.com/torvalds/linux/blob/v6.8/arch/x86/include/asm/msr.h
#include "Wh2FrozenTrace.h"
#include <array>
#include <cpuid.h>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sched.h>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <time.h>

#if !defined(WH2_CLOCK_NEUTRAL_ONLY) || !defined(__x86_64__) || !defined(__linux__)
#error "Explicit Linux x86-64 diagnostic build required"
#endif
#define NOINLINE __attribute__((noinline, noipa))
namespace {
constexpr unsigned Samples=524288, NeutralSamples=256, Iterations=65536, Chunk=256;
constexpr uint64_t Seed=UINT64_C(0x9e3779b97f4a7c15), Answer=UINT64_C(0xd63028b03dfd593c);
const char Protocol[]="wirehair.wh2.clock-boundary-r0";
const char ClaimPath[]="/var/tmp/wh2-clock-boundary-r0/CLAIM.json";
void Check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
uint64_t Clock(clockid_t clock) {
    timespec t={}; Check(!clock_gettime(clock,&t) && t.tv_sec>=0 && t.tv_nsec>=0 && t.tv_nsec<1000000000,"clock");
    Check(uint64_t(t.tv_sec)<UINT64_MAX/1000000000-1,"clock overflow");
    return uint64_t(t.tv_sec)*1000000000+uint64_t(t.tv_nsec);
}
using Usage=std::array<uint64_t,4>;
Usage Counters() {
    rusage r={}; Check(!getrusage(RUSAGE_THREAD,&r) && r.ru_minflt>=0 && r.ru_majflt>=0 &&
                      r.ru_nvcsw>=0 && r.ru_nivcsw>=0,"rusage");
    return {{uint64_t(r.ru_minflt),uint64_t(r.ru_majflt),uint64_t(r.ru_nvcsw),uint64_t(r.ru_nivcsw)}};
}
bool OnCpu() {
    cpu_set_t mask; CPU_ZERO(&mask);
    return !sched_getaffinity(0,sizeof(mask),&mask) && CPU_COUNT(&mask)==1 &&
        CPU_ISSET(50,&mask) && sched_getcpu()==50;
}
void Pin() {
    cpu_set_t mask; CPU_ZERO(&mask); CPU_SET(50,&mask);
    Check(!sched_setaffinity(0,sizeof(mask),&mask) && OnCpu(),"CPU50 affinity");
}
using Leaf=std::array<unsigned,5>;
using Identity=std::array<Leaf,7>;
Identity Features() {
    Check(__get_cpuid_max(0x80000000u,nullptr)>=0x80000021u,"extended CPUID");
    Identity result={}; const unsigned ids[7]={0,1,0x80000000u,0x80000001u,0x80000007u,0x80000008u,0x80000021u};
    for(unsigned i=0;i<7;++i) {
        auto& l=result[i]; l[0]=ids[i]; __cpuid_count(ids[i],0,l[1],l[2],l[3],l[4]);
    }
    char vendor[13]={}; std::memcpy(vendor,&result[0][2],4); std::memcpy(vendor+4,&result[0][4],4);
    std::memcpy(vendor+8,&result[0][3],4);
    Check(!std::strcmp(vendor,"AuthenticAMD") && result[1][1]==0x00b00f81u &&
          (result[1][2]>>24)==100 && !(result[1][3]&(1u<<31)),"frozen physical CPU");
    Check((result[3][4]&(1u<<27)) && (result[4][4]&(1u<<8)) && (result[6][1]&(1u<<2)),
          "RDTSCP, invariant TSC and always-serializing LFENCE required");
    return result;
}
struct Stamp { uint64_t tsc; unsigned aux; };
inline Stamp Tsc() {
    unsigned lo,hi,aux;
    __asm__ __volatile__("mfence\n\tlfence\n\trdtscp\n\tlfence"
        : "=a"(lo),"=d"(hi),"=c"(aux) : : "memory");
    return Stamp{uint64_t(lo)|(uint64_t(hi)<<32),aux};
}
NOINLINE uint64_t Compute(uint64_t value,unsigned count) {
    for(unsigned i=0;i<count;++i) { value^=value<<13; value^=value>>7; value^=value<<17; }
    return value;
}
// [index,cpu0,cpu1,mono0,mono1,t0,t1,t2,t3,result,
//  minor0,minor1,major0,major1,voluntary0,voluntary1,involuntary0,involuntary1,
//  aux0,aux1,aux2,aux3,stage]. Stages expose incomplete exception captures.
using Record=std::array<uint64_t,23>;
Record records[Chunk];
NOINLINE void Observe(Record& r,bool fail_clock) {
    const Usage before=Counters(); for(unsigned j=0;j<4;++j) r[10+2*j]=before[j]; r[22]=1;
    r[1]=Clock(CLOCK_THREAD_CPUTIME_ID); r[22]=2;
    Stamp t=Tsc(); r[5]=t.tsc; r[18]=t.aux; r[22]=3;
    r[3]=Clock(CLOCK_MONOTONIC); r[22]=4;
    t=Tsc(); r[6]=t.tsc; r[19]=t.aux; r[22]=5;
    r[9]=Compute(Seed,Iterations); r[22]=6;
    t=Tsc(); r[7]=t.tsc; r[20]=t.aux; r[22]=7;
    if(fail_clock) throw std::runtime_error("neutral clock failure");
    r[4]=Clock(CLOCK_MONOTONIC); r[22]=8;
    t=Tsc(); r[8]=t.tsc; r[21]=t.aux; r[22]=9;
    r[2]=Clock(CLOCK_THREAD_CPUTIME_ID); r[22]=10;
    const Usage after=Counters(); for(unsigned j=0;j<4;++j) r[11+2*j]=after[j]; r[22]=11;
}
void Validate(const Record& r,unsigned aux) {
    Check(r[22]==11 && r[9]==Answer,"complete computation");
    Check(r[1]<=r[2] && r[3]<r[4],"clock order");
    for(unsigned j=5;j<8;++j) Check(r[j]<r[j+1],"TSC order");
    for(unsigned j=18;j<22;++j) Check(r[j]==aux,"TSC AUX migration");
    for(unsigned j=0;j<4;++j) Check(r[10+2*j]<=r[11+2*j],"counter order");
    Check(OnCpu(),"affinity changed");
}
template<class T,size_t N> void Array(const std::array<T,N>& a) {
    putchar('['); for(size_t i=0;i<N;++i) { if(i) putchar(','); printf("%llu",static_cast<unsigned long long>(a[i])); } putchar(']');
}
void Flush() { Check(fflush(stdout)==0 && !ferror(stdout),"output stream"); }
void Publish(unsigned count) { for(unsigned i=0;i<count;++i) { Array(records[i]); putchar('\n'); } Flush(); }
void Authenticate(const char* path,const std::string& claim) {
    Check(claim.size()==64 && claim.find_first_not_of("0123456789abcdef")==std::string::npos,"claim hex");
    std::ifstream input(path,std::ios::binary); Check(bool(input),"claimed namespace");
    std::string bytes; char c;
    while(input.get(c)) { Check(bytes.size()<1024*1024,"claim cap"); bytes+=c; }
    Check(input.eof() && wirehair::wh2_benchmark::Sha256Hex(bytes)==claim,"claim binding");
}
int Run(bool neutral,bool bad_result,bool fail_clock,const std::string& claim) {
    if(!neutral) { Check(!WH2_CLOCK_NEUTRAL_ONLY,"neutral worker cannot measure"); Authenticate(ClaimPath,claim); }
    const rlimit cpu={100,100},memory={128*1024*1024,128*1024*1024},core={0,0};
    Check(!setrlimit(RLIMIT_CPU,&cpu) && !setrlimit(RLIMIT_CORE,&core),"CPU/core limits");
    // ASan reserves a large virtual address space; it can never enter science mode.
    if(!WH2_CLOCK_NEUTRAL_ONLY) Check(!setrlimit(RLIMIT_AS,&memory),"memory limit");
    Pin(); const Identity identity=Features(); const Stamp initial=Tsc();
    Check(initial.aux==50,"frozen TSC AUX");
    const uint64_t start=Clock(CLOCK_MONOTONIC),cpu_start=Clock(CLOCK_THREAD_CPUTIME_ID);
    const unsigned count=neutral?NeutralSamples:Samples;
    printf("{\"type\":\"header\",\"protocol\":\"%s\",\"claim\":\"%s\",\"neutral\":%s,\"samples\":%u,\"iterations\":%u,\"seed\":%llu,\"answer\":%llu,\"aux\":%u,\"start\":[%llu,%llu],\"cpuid\":[",
        Protocol,claim.c_str(),neutral?"true":"false",count,Iterations,
        static_cast<unsigned long long>(Seed),static_cast<unsigned long long>(Answer),initial.aux,
        static_cast<unsigned long long>(start),static_cast<unsigned long long>(cpu_start));
    for(unsigned j=0;j<7;++j) { if(j) putchar(','); Array(identity[j]); } printf("]}\n"); Flush();
    unsigned pending=0,retained=0; bool complete=false; const char* why=nullptr;
    try {
        for(unsigned i=0;i<count;++i) {
            Record& r=records[pending++]; r.fill(0); r[0]=i; ++retained;
            Observe(r,fail_clock && i+1==count);
            if(bad_result && i+1==count) r[9]^=1;
            Validate(r,initial.aux);
            Check(r[4]-start<UINT64_C(120000000000),"wall limit");
            if(pending==Chunk) { Publish(pending); pending=0; }
        }
        Check(Features()==identity,"target changed"); complete=true;
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); why="diagnostic failure"; }
    if(pending) Publish(pending);
    printf("{\"type\":\"footer\",\"complete\":%s,\"records\":%u,\"codec_calls\":0,\"end\":[%llu,%llu]}\n",
        complete?"true":"false",retained,static_cast<unsigned long long>(Clock(CLOCK_MONOTONIC)),
        static_cast<unsigned long long>(Clock(CLOCK_THREAD_CPUTIME_ID))); Flush();
    return why?1:0;
}
} // namespace
int main(int argc,char** argv) {
    try {
        if(argc==2 && !std::strcmp(argv[1],"--contract")) {
            printf("{\"samples\":%u,\"neutral_samples\":%u,\"iterations\":%u,\"cpu_seconds\":100,\"wall_seconds\":120,\"address_space_mib\":128,\"claim_path\":\"%s\",\"neutral_only\":%s}\n",
                Samples,NeutralSamples,Iterations,ClaimPath,WH2_CLOCK_NEUTRAL_ONLY?"true":"false"); return 0;
        }
        if(argc==2 && !std::strcmp(argv[1],"--neutral")) return Run(true,false,false,std::string(64,'0'));
        if(argc==2 && !std::strcmp(argv[1],"--neutral-bad-result")) return Run(true,true,false,std::string(64,'0'));
        if(argc==2 && !std::strcmp(argv[1],"--neutral-fail-clock")) return Run(true,false,true,std::string(64,'0'));
        if(argc==4 && !std::strcmp(argv[1],"--neutral-claim")) {
            Authenticate(argv[2],argv[3]); puts("PASS claim authentication"); return 0;
        }
        Check(argc==3 && !std::strcmp(argv[1],"--worker"),"explicit mode required");
        return Run(false,false,false,argv[2]);
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); return 1; }
}
