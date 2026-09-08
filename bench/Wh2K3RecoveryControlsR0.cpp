// Retained K3 recovery comparison through actual public library APIs.
// No new trace generation, candidate selection, timing or default change.
#include "wirehair/wirehair_small.h"
#include "Wh2FrozenTrace.h"
#include "Wh2K3NativeData.inc"
#include "gf256.h"
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>
#include <time.h>

#if !defined(WH2_K3_RECOVERY_BACKEND) || defined(WH_COUNT) || defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "Explicit unmodified recovery build required"
#endif
static_assert(WH2_K3_RECOVERY_BACKEND >= 0 && WH2_K3_RECOVERY_BACKEND <= 2, "backend");
namespace {
using Byte = uint8_t;
using Profile = std::array<Byte,32>;
const char protocol[] = "wirehair.wh2.k3-recovery-controls-r0";
const char* const backend_names[] = {"native","scalar","asan"};
const unsigned widths[] = {2,64,1280};
void Check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
std::string Sha(const void* data,size_t size) { return wirehair::wh2_benchmark::Sha256Hex(data,size); }
struct Api {
    int (*create)(const void*,unsigned,Profile&,void*&);
    int (*encode)(void*,uint32_t,void*,unsigned,uint32_t*);
    int (*decoder)(const Profile&,unsigned,void*&);
    int (*feed)(void*,uint32_t,const void*,unsigned);
    int (*recover)(void*,void*,unsigned,uint64_t*);
    void (*free)(void*);
};
int PCreate(const void* source,unsigned b,Profile& p,void*& h) {
    WirehairV2Codec handle=nullptr; uint32_t n=0;
    WirehairV2EncoderOptions o=WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    o.source_policy=WirehairV2EncoderSource_BorrowedImmutable;
    const int r=wirehair_v2_encoder_create_with_options(source,3u*b,b,&o,p.data(),32,&n,&handle);
    h=handle; return r==0 && n!=32?-1:r;
}
int PEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    return wirehair_v2_encode(static_cast<WirehairV2Codec>(h),id,out,b,n);
}
int PDecoder(const Profile& p,unsigned,void*& h) {
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
int CCreate(const void* source,unsigned b,Profile& p,void*& h) {
    const auto r=wirehair_small_encoder_create(source,3u*b,b,WirehairSmall_BorrowedImmutable,p.data(),32);
    h=r.codec; return r.status;
}
int CEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const auto r=wirehair_small_encode(static_cast<WirehairSmallCodec>(h),id,out,b);
    *n=static_cast<uint32_t>(r.bytes_written);
    return r.status==WirehairSmall_Success && (r.bytes_required!=b || r.bytes_written!=b)?-1:int(r.status);
}
int CDecoder(const Profile& p,unsigned,void*& h) {
    const auto r=wirehair_small_decoder_create(p.data(),32); h=r.codec; return r.status;
}
int CFeed(void* h,uint32_t id,const void* in,unsigned n) {
    return wirehair_small_decode(static_cast<WirehairSmallCodec>(h),id,in,n);
}
int CRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const auto r=wirehair_small_recover(static_cast<WirehairSmallCodec>(h),out,n); *written=r.bytes_written;
    return r.status==WirehairSmall_Success && r.bytes_required!=n?-1:int(r.status);
}
void CFree(void* h) { wirehair_small_free(static_cast<WirehairSmallCodec>(h)); }
int WCreate(const void* source,unsigned b,Profile&,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_encoder_create_ex(nullptr,source,3u*b,b,&handle); h=handle; return r;
}
int WEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    return wirehair_encode(static_cast<WirehairCodec>(h),id,out,b,n);
}
int WDecoder(const Profile&,unsigned b,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_decoder_create_ex(nullptr,3u*b,b,&handle); h=handle; return r;
}
int WFeed(void* h,uint32_t id,const void* in,unsigned n) { return wirehair_decode(static_cast<WirehairCodec>(h),id,in,n); }
int WRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const int r=wirehair_recover(static_cast<WirehairCodec>(h),out,n); *written=r==0?n:0; return r;
}
void WFree(void* h) { wirehair_free(static_cast<WirehairCodec>(h)); }
const Api apis[3]={
    {PCreate,PEncode,PDecoder,PFeed,PRecover,PFree},
    {CCreate,CEncode,CDecoder,CFeed,CRecover,CFree},
    {WCreate,WEncode,WDecoder,WFeed,WRecover,WFree}};
struct Owner {
    const Api& api; void* handle=nullptr;
    explicit Owner(const Api& a):api(a) {}
    ~Owner() { if(handle) api.free(handle); }
    Owner(const Owner&)=delete; Owner& operator=(const Owner&)=delete;
};

struct Oracle {
    using Row = std::array<Byte,3>;
    using Matrix = std::array<Byte,9>;
    Byte products[256][256] = {};
    Matrix powers[2][32];
    static Byte Multiply(Byte a,Byte b) {
        unsigned p=0;
        for(unsigned i=0;i<8;++i) if(b & (1u<<i)) p^=unsigned(a)<<i;
        for(int i=14;i>=8;--i) if(p & (1u<<i)) p^=0x14du<<(i-8);
        return Byte(p);
    }
    Matrix Product(const Matrix& a,const Matrix& b) const {
        Matrix p={};
        for(unsigned r=0;r<3;++r) for(unsigned c=0;c<3;++c)
            for(unsigned k=0;k<3;++k) p[3*r+c]^=products[a[3*r+k]][b[3*k+c]];
        return p;
    }
    Oracle() {
        for(unsigned a=0;a<256;++a) for(unsigned b=0;b<256;++b) products[a][b]=Multiply(Byte(a),Byte(b));
        for(unsigned p=0;p<2;++p) powers[p][0]={{0,0,Byte(8+p),1,0,14,0,1,7}};
        for(unsigned i=1;i<32;++i) {
            powers[0][i]=Product(powers[0][i-1],powers[1][i-1]);
            powers[1][i]=Product(powers[1][i-1],powers[0][i-1]);
        }
    }
    Row Coefficients(uint32_t id) const {
        Row row={{1,0,0}};
        for(unsigned bit=0;bit<32;++bit) if(id & (uint32_t(1)<<bit)) {
            unsigned phase=0;
            for(unsigned higher=bit+1;higher<32;++higher) phase^=(id>>higher)&1u;
            Row next={};
            for(unsigned r=0;r<3;++r) for(unsigned c=0;c<3;++c)
                next[r]^=products[powers[phase][bit][3*r+c]][row[c]];
            row=next;
        }
        return row;
    }
    std::vector<Byte> Packet(const std::vector<Byte>& source,unsigned B,uint32_t id) const {
        const Row row=Coefficients(id); std::vector<Byte> result(B,0);
        for(unsigned j=0;j<B;++j) for(unsigned k=0;k<3;++k)
            result[j]^=products[row[k]][source[k*B+j]];
        return result;
    }
    unsigned Rank(std::vector<Row> rows) const {
        unsigned rank=0;
        for(unsigned c=0;c<3 && rank<rows.size();++c) {
            unsigned p=rank; while(p<rows.size() && !rows[p][c]) ++p;
            if(p==rows.size()) continue;
            std::swap(rows[p],rows[rank]);
            unsigned inverse=1; while(inverse<256 && products[rows[rank][c]][inverse]!=1) ++inverse;
            Check(inverse<256,"polynomial inverse");
            for(Byte& b:rows[rank]) b=products[b][inverse];
            for(unsigned r=0;r<rows.size();++r) if(r!=rank) {
                const Byte factor=rows[r][c];
                for(unsigned k=0;k<3;++k) rows[r][k]^=products[factor][rows[rank][k]];
            }
            ++rank;
        }
        return rank;
    }
};
std::vector<Byte> Message(unsigned B) {
    std::vector<Byte> source(3u*B);
    for(size_t i=0;i<source.size();++i) source[i]=Byte(37*i+i/11);
    return source;
}
struct ArmResult {
    Profile profile={};
    std::vector<std::string> packets;
    std::vector<int> feed;
    unsigned first=0,recoveries=0;
};
using Result = std::array<ArmResult,3>;
Result Exercise(const Oracle& oracle,unsigned B,const uint32_t* ids,unsigned count) {
    Check((B==2 || B==64 || B==1280) && count>=3 && count<=7,"bounded retained shape");
    const auto source=Message(B);
    const auto source_copy=source;
    std::vector<Oracle::Row> rows; unsigned candidate_first=0;
    std::vector<std::vector<Byte>> independent;
    for(unsigned i=0;i<count;++i) {
        rows.push_back(oracle.Coefficients(ids[i]));
        if(!candidate_first && oracle.Rank(rows)==3) candidate_first=i+1;
        independent.push_back(oracle.Packet(source,B,ids[i]));
    }
    Result result;
    for(unsigned a=0;a<3;++a) {
        auto& r=result[a]; const auto& api=apis[a];
        std::vector<std::vector<Byte>> packets(count,std::vector<Byte>(B+2,0xa5));
        {
            Owner encoder(api);
            Check(api.create(source.data(),B,r.profile,encoder.handle)==0 && encoder.handle,"encoder create");
            for(unsigned i=0;i<count;++i) {
                uint32_t written=0;
                Check(api.encode(encoder.handle,ids[i],packets[i].data()+1,B,&written)==0 && written==B,
                      "encode status/length");
                Check(packets[i].front()==0xa5 && packets[i].back()==0xa5,"packet guards");
                if(a==1) Check(!memcmp(packets[i].data()+1,independent[i].data(),B),"independent polynomial payload");
                r.packets.push_back(Sha(packets[i].data()+1,B));
            }
        } // No encoder or source-derived private state survives into the receiver.
        if(a==0) Check(wirehair_v2_profile_validate(r.profile.data(),32)==WirehairV2_Success,"WH2 profile");
        if(a==1) Check(wirehair_small_profile_validate(r.profile.data(),32)==WirehairSmall_Success,"K3 profile");
        Owner decoder(api);
        Check(api.decoder(r.profile,B,decoder.handle)==0 && decoder.handle,"standalone decoder create");
        for(unsigned i=0;i<count;++i) {
            const int status=api.feed(decoder.handle,ids[i],packets[i].data()+1,B);
            Check(status==0 || status==1,"feed status");
            r.feed.push_back(status);
            if(status==0) { Check(i>=2,"premature success"); r.first=i+1; break; }
        }
        if(r.first) for(unsigned repeat=0;repeat<2;++repeat) {
            std::vector<Byte> output(source.size()+2,0xa5); uint64_t written=0;
            Check(api.recover(decoder.handle,output.data()+1,unsigned(source.size()),&written)==0 &&
                  written==source.size(),"recover status/length");
            Check(output.front()==0xa5 && output.back()==0xa5 &&
                  !memcmp(output.data()+1,source.data(),source.size()),"recovered payload/guards");
            ++r.recoveries;
        }
        if(a==1) Check(r.first==candidate_first,"independent first-success rank");
        Check(source==source_copy,"immutable source");
        for(unsigned i=0;i<count;++i)
            Check(packets[i].front()==0xa5 && packets[i].back()==0xa5 &&
                  Sha(packets[i].data()+1,B)==r.packets[i],"immutable packet");
    }
    return result;
}
void Flush() { Check(fflush(stdout)==0 && !ferror(stdout),"output stream"); }
void Hex(const void* data,size_t n) {
    const auto* p=static_cast<const Byte*>(data); putchar('"');
    for(size_t i=0;i<n;++i) printf("%02x",unsigned(p[i]));
    putchar('"');
}
void Record(unsigned group,unsigned index,unsigned B,const uint32_t* ids,unsigned count,const Result& result) {
    printf("{\"type\":\"record\",\"group\":%u,\"index\":%u,\"B\":%u,\"ids\":[",group,index,B);
    for(unsigned i=0;i<count;++i) { if(i) putchar(','); printf("%u",ids[i]); }
    printf("],\"arms\":[");
    for(unsigned a=0;a<3;++a) {
        if(a) putchar(',');
        const auto& r=result[a];
        printf("{\"profile\":"); Hex(r.profile.data(),32); printf(",\"packets\":[");
        for(size_t i=0;i<r.packets.size();++i) { if(i) putchar(','); printf("\"%s\"",r.packets[i].c_str()); }
        printf("],\"feed\":[");
        for(size_t i=0;i<r.feed.size();++i) { if(i) putchar(','); printf("%d",r.feed[i]); }
        printf("],\"first\":%u,\"recoveries\":%u,\"checked\":true}",r.first,r.recoveries);
    }
    printf("]}\n"); Flush();
}
gf256_x86_cpu_features Initialize(const Oracle& oracle) {
    Check(wirehair_init()==Wirehair_Success && GF256Ctx.Polynomial==0x14d,"shared GF initialization");
    Check(Sha(wh2_k3_data::kLookup,sizeof(wh2_k3_data::kLookup))==wh2_k3_data::kLookupSha,"sealed lookup");
    for(const auto& row:wh2_k3_data::kRows) {
        const auto expected=oracle.Coefficients(row.id);
        Check(std::equal(expected.begin(),expected.end(),row.values),"sealed polynomial coefficients");
    }
    gf256_x86_cpu_features f={}; gf256_get_active_x86_cpu_features(&f);
    if(WH2_K3_RECOVERY_BACKEND==1) Check(!f.SSSE3 && !f.AVX2 && !f.GFNI && !f.AVX512,"portable backend");
    return f;
}
uint64_t Now() {
    timespec t={};
    Check(!clock_gettime(CLOCK_MONOTONIC,&t) && t.tv_sec>=0,"monotonic clock");
    return uint64_t(t.tv_sec)*1000000000+uint64_t(t.tv_nsec);
}
int Worker(const std::string& claim) {
    Check(claim.size()==64 && claim.find_first_not_of("0123456789abcdef")==std::string::npos,"claim hex");
    std::ifstream in("/var/tmp/wh2-k3-recovery-controls-r0/CLAIM.json",std::ios::binary);
    Check(bool(in),"claimed namespace"); std::string receipt; char ch;
    while(in.get(ch)) { Check(receipt.size()<1024*1024,"claim cap"); receipt+=ch; }
    Check(in.eof() && Sha(receipt.data(),receipt.size())==claim,"claim bytes");
    const rlimit cpu={60,60},memory={256u*1024u*1024u,256u*1024u*1024u},core={0,0};
    Check(!setrlimit(RLIMIT_CPU,&cpu) && !setrlimit(RLIMIT_CORE,&core),"worker limits");
    if(WH2_K3_RECOVERY_BACKEND!=2) Check(!setrlimit(RLIMIT_AS,&memory),"worker memory cap");
    const uint64_t start=Now(); const Oracle oracle; const auto f=Initialize(oracle);
    printf("{\"type\":\"header\",\"protocol\":\"%s\",\"claim\":\"%s\",\"backend\":\"%s\","
           "\"retained_raw_sha256\":\"%s\",\"features\":[%u,%u,%u,%u],\"sources\":[",
           protocol,claim.c_str(),backend_names[WH2_K3_RECOVERY_BACKEND],wh2_k3_data::kRawSha,
           unsigned(f.SSSE3),unsigned(f.AVX2),unsigned(f.GFNI),unsigned(f.AVX512));
    for(unsigned B:widths) { if(B!=2) putchar(','); const auto m=Message(B); printf("\"%s\"",Sha(m.data(),m.size()).c_str()); }
    printf("]}\n"); Flush();
    unsigned records=0;
    for(unsigned i=0;i<6216;++i) {
        const auto& t=wh2_k3_data::kTraces[i];
        Check(Now()-start<UINT64_C(90000000000),"worker deadline");
        const auto result=Exercise(oracle,t.B,t.ids,7);
        Record(i<6144?0:1,i<6144?i:i-6144,t.B,t.ids,7,result); ++records;
    }
    for(unsigned i=0;i<45;++i) for(unsigned w=0;w<3;++w) {
        const auto& p=wh2_k3_data::kHistory[i]; if(!(p.widths & (1u<<w))) continue;
        Check(Now()-start<UINT64_C(90000000000),"worker deadline");
        const auto result=Exercise(oracle,widths[w],p.ids,p.count);
        Record(2,i,widths[w],p.ids,p.count,result); ++records;
    }
    Check(records==6269 && Now()-start<UINT64_C(90000000000),"complete roster/deadline");
    printf("{\"type\":\"footer\",\"records\":%u,\"checked\":true}\n",records); Flush(); return 0;
}
int Neutral() {
    const Oracle oracle; Initialize(oracle);
    for(unsigned a=0;a<256;++a) for(unsigned b=0;b<256;++b)
        Check(gf256_mul(Byte(a),Byte(b))==oracle.products[a][b],"exhaustive independent field");
    const uint32_t streams[4][7]={{0,1,2,3,4,5,6},{3,4,5,6,7,8,9},
        {UINT32_MAX,UINT32_MAX-2,UINT32_MAX-4,UINT32_MAX-6,UINT32_MAX-8,UINT32_MAX-10,UINT32_MAX-12},
        {0,0,0,0,0,0,0}};
    unsigned cases=0;
    for(unsigned B:widths) for(unsigned s=0;s<4;++s) {
        const auto result=Exercise(oracle,B,streams[s],7);
        Check(result[1].first==(s==3?0u:3u),"neutral endpoint"); ++cases;
    }
    printf("PASS 12 neutral cases, 3 actual APIs, independent field/payload/rank and recovery guards (%u checked)\n",cases);
    return 0;
}
static_assert(sizeof(wh2_k3_data::kTraces)/sizeof(wh2_k3_data::kTraces[0])==6216,"trace roster");
static_assert(sizeof(wh2_k3_data::kHistory)/sizeof(wh2_k3_data::kHistory[0])==45,"history roster");
static_assert(sizeof(wh2_k3_data::kWindows)/sizeof(wh2_k3_data::kWindows[0])==43,"window provenance");
static_assert(sizeof(wh2_k3_data::kRows)/sizeof(wh2_k3_data::kRows[0])==2226,"coefficient provenance");
} // namespace
int main(int argc,char** argv) {
    try {
        if(argc==2 && !strcmp(argv[1],"--neutral")) return Neutral();
        Check(argc==3 && !strcmp(argv[1],"--worker"),"explicit mode required");
        return Worker(argv[2]);
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); return 1; }
}
