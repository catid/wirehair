// Frozen K5 retained recovery comparison; no timing or new trace selection.
#include "wirehair/wirehair.h"
#include "Wh2SmallSerialized.h"
#include "Wh2FrozenTrace.h"
#include "Wh2K5NativeData.inc"
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

#if !defined(WH2_K5_SERIALIZED_RECOVERY_BACKEND) || defined(WH_COUNT) || defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "Explicit unmodified recovery backend required"
#endif
static_assert(WH2_SMALL_CODEC_K==5, "sealed K5 boundary");
static_assert(WH2_K5_SERIALIZED_RECOVERY_BACKEND>=0 && WH2_K5_SERIALIZED_RECOVERY_BACKEND<=2, "backend");
namespace {
using Byte=uint8_t;
using Profile=std::array<Byte,32>;
const unsigned K=5, arms=6, cpu_seconds=120, wall_seconds=150, address_space_mib=384;
const char protocol[]="wirehair.wh2.k5-serialized-recovery-r0";
const char claim_path[]="/var/tmp/wh2-k5-serialized-recovery-r0/CLAIM.json";
const char* const backends[]={"native","scalar","asan"};
const unsigned widths[]={2,64,1280};
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
int PCreate(const void* source,unsigned b,Profile& p,void*& h,uint32_t policy) {
    WirehairV2Codec handle=nullptr; uint32_t n=0;
    WirehairV2EncoderOptions o=WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    o.source_policy=policy;
    const int r=wirehair_v2_encoder_create_profile_id_with_options(
        WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,source,K*b,b,&o,p.data(),32,&n,&handle);
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
int PIndependent(const void* s,unsigned b,Profile& p,void*& h) {
    return PCreate(s,b,p,h,WirehairV2EncoderSource_Independent);
}
int PBorrowed(const void* s,unsigned b,Profile& p,void*& h) {
    return PCreate(s,b,p,h,WirehairV2EncoderSource_BorrowedImmutable);
}
int SCreate(const void* source,unsigned b,Profile& p,void*& h,uint32_t policy) {
    const auto r=wh2_small_encoder_create(source,K*b,b,policy,p.data(),p.size()); h=r.codec; return r.status;
}
int SIndependent(const void* s,unsigned b,Profile& p,void*& h) { return SCreate(s,b,p,h,Wh2Small_Independent); }
int SBorrowed(const void* s,unsigned b,Profile& p,void*& h) { return SCreate(s,b,p,h,Wh2Small_BorrowedImmutable); }
int SEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const auto r=wh2_small_encode(h,id,out,b); *n=static_cast<uint32_t>(r.bytes_written);
    return r.status==Wh2Small_Success && (r.bytes_written!=b || r.bytes_required!=b)?-1:static_cast<int>(r.status);
}
int SDecoder(const Profile& p,unsigned,void*& h) {
    const auto r=wh2_small_decoder_create(p.data(),p.size()); h=r.codec; return r.status;
}
int SFeed(void* h,uint32_t id,const void* in,unsigned n) { return wh2_small_decode(h,id,in,n); }
int SRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const auto r=wh2_small_recover(h,out,n); *written=r.bytes_written;
    return r.status==Wh2Small_Success && (r.bytes_required!=n || r.bytes_written!=n)?-1:static_cast<int>(r.status);
}
void SFree(void* h) { wh2_small_free(h); }
int WCreate(const void* source,unsigned b,Profile&,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_encoder_create_ex(nullptr,source,K*b,b,&handle); h=handle; return r;
}
int WOwned(const void* source,unsigned b,Profile&,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_encoder_create_owned_ex(nullptr,source,K*b,b,&handle); h=handle; return r;
}
int WEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    return wirehair_encode(static_cast<WirehairCodec>(h),id,out,b,n);
}
int WDecoder(const Profile&,unsigned b,void*& h) {
    WirehairCodec handle=nullptr;
    const int r=wirehair_decoder_create_ex(nullptr,K*b,b,&handle); h=handle; return r;
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

using Ledger=std::array<unsigned,6>;
struct Owner {
    const Api& api; Ledger& ledger; void* handle=nullptr;
    Owner(const Api& a,Ledger& l):api(a),ledger(l) {}
    ~Owner() { if(handle) { ++ledger[5]; api.free(handle); } }
    Owner(const Owner&)=delete; Owner& operator=(const Owner&)=delete;
};
struct Oracle {
    using Row=std::array<Byte,K>;
    using Matrix=std::array<Byte,K*K>;
    Byte products[256][256]={};
    Matrix powers[2][32];
    static Byte Multiply(Byte a,Byte b) {
        unsigned p=0;
        for(unsigned bit=0;bit<8;++bit) if(b & (1u<<bit)) p^=unsigned(a)<<bit;
        for(int bit=14;bit>=8;--bit) if(p & (1u<<bit)) p^=0x14du<<(bit-8);
        return Byte(p);
    }
    Matrix Product(const Matrix& a,const Matrix& b) const {
        Matrix out={};
        for(unsigned i=0;i<K;++i) for(unsigned j=0;j<K;++j)
            for(unsigned k=0;k<K;++k) out[K*i+j]^=products[a[K*i+k]][b[K*k+j]];
        return out;
    }
    Oracle() {
        for(unsigned a=0;a<256;++a) for(unsigned b=0;b<256;++b) products[a][b]=Multiply(Byte(a),Byte(b));
        const Byte feedback[K]={121,110,207,198,31};
        for(unsigned phase=0;phase<2;++phase) {
            Matrix m={};
            for(unsigned i=0;i<K-1;++i) m[K*(i+1)+i]=1;
            for(unsigned i=0;i<K;++i) m[K*i+K-1]=Byte(feedback[i]^(i==0?phase:0));
            powers[phase][0]=m;
        }
        for(unsigned bit=1;bit<32;++bit) {
            powers[0][bit]=Product(powers[0][bit-1],powers[1][bit-1]);
            powers[1][bit]=Product(powers[1][bit-1],powers[0][bit-1]);
        }
    }
    Row Coefficients(uint32_t id) const {
        Row row={{1,0,0,0,0}};
        for(unsigned bit=0;bit<32;++bit) if(id & (uint32_t(1)<<bit)) {
            unsigned phase=0;
            for(unsigned higher=bit+1;higher<32;++higher) phase^=(id>>higher)&1u;
            Row next={};
            for(unsigned r=0;r<K;++r) for(unsigned c=0;c<K;++c)
                next[r]^=products[powers[phase][bit][K*r+c]][row[c]];
            row=next;
        }
        return row;
    }
    std::vector<Byte> Packet(const std::vector<Byte>& source,unsigned B,const Row& row) const {
        std::vector<Byte> packet(B,0);
        for(unsigned j=0;j<B;++j) for(unsigned k=0;k<K;++k)
            packet[j]^=products[row[k]][source[k*B+j]];
        return packet;
    }
    unsigned Rank(std::vector<Row> rows) const {
        unsigned rank=0;
        for(unsigned c=0;c<K && rank<rows.size();++c) {
            unsigned p=rank; while(p<rows.size() && !rows[p][c]) ++p;
            if(p==rows.size()) continue;
            std::swap(rows[p],rows[rank]);
            unsigned inverse=1; while(inverse<256 && products[rows[rank][c]][inverse]!=1) ++inverse;
            Check(inverse<256,"polynomial inverse");
            for(Byte& b:rows[rank]) b=products[b][inverse];
            for(unsigned r=0;r<rows.size();++r) if(r!=rank) {
                const Byte factor=rows[r][c];
                for(unsigned k=0;k<K;++k) rows[r][k]^=products[factor][rows[rank][k]];
            }
            ++rank;
        }
        return rank;
    }
};
std::vector<Byte> Message(unsigned B) {
    std::vector<Byte> source(K*B);
    for(size_t i=0;i<source.size();++i) source[i]=Byte(37*i+i/11);
    return source;
}
struct ArmResult {
    Profile profile={};
    std::vector<std::string> packets;
    std::vector<Oracle::Row> rows;
    std::vector<int> feed;
    unsigned first=0,recoveries=0;
    Ledger counts={};
};
using Result=std::array<ArmResult,arms>;
Result Exercise(const Oracle& oracle,unsigned B,const uint32_t* ids,unsigned count,const Api* api_set=apis) {
    Check((B==2 || B==64 || B==1280) && count>=K && count<=K+4,"bounded retained shape");
    const auto source=Message(B), original=source;
    Result result;
    std::array<std::vector<std::vector<Byte>>,arms> packets;
    for(unsigned a=0;a<arms;++a) {
        auto& r=result[a]; const auto& api=api_set[a];
        packets[a].assign(count,std::vector<Byte>(B+2,0xa5)); r.rows.resize(count);
        for(unsigned basis=0;basis<=K;++basis) {
            std::vector<Byte> input=basis?std::vector<Byte>(K*B,0):source;
            if(basis) input[(basis-1)*B]=1;
            Owner encoder(api,r.counts); Profile profile={};
            ++r.counts[0];
            Check(api.create(input.data(),B,profile,encoder.handle)==0 && encoder.handle,"encoder create");
            if(!basis) r.profile=profile; else Check(profile==r.profile,"basis profile identity");
            if(a<3) { std::fill(input.begin(),input.end(),0xcc); std::vector<Byte>().swap(input); }
            for(unsigned i=0;i<count;++i) {
                std::vector<Byte> packet(B+2,0xa5); uint32_t written=0; ++r.counts[1];
                Check(api.encode(encoder.handle,ids[i],packet.data()+1,B,&written)==0 && written==B,"encode status/length");
                Check(packet.front()==0xa5 && packet.back()==0xa5,"packet guards");
                if(basis) {
                    r.rows[i][basis-1]=packet[1];
                    for(unsigned j=1;j<B;++j) Check(packet[j+1]==0,"basis other bytes");
                } else {
                    r.packets.push_back(Sha(packet.data()+1,B));
                    packets[a][i]=std::move(packet);
                }
            }
        }
        if(a==1 || a==4) Check(wh2_small_profile_validate(r.profile.data(),32)==Wh2Small_Success,"sealed K5 descriptor");
        if(a==0 || a==3) {
            WirehairV2Profile p={};
            Check(wirehair_v2_profile_deserialize(r.profile.data(),32,&p)==WirehairV2_Success &&
                  p.profile_id==WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,"explicit certified profile");
        }
        for(unsigned i=0;i<count;++i) {
            if(a==1 || a==4) Check(r.rows[i]==oracle.Coefficients(ids[i]),"independent sealed coefficients");
            const auto expected=oracle.Packet(source,B,r.rows[i]);
            Check(!memcmp(packets[a][i].data()+1,expected.data(),B),"independent every-arm payload");
        }
    } // Every real and basis encoder across all six arms is now gone.
    for(unsigned a=0;a<arms;++a) {
        auto& r=result[a]; const auto& api=api_set[a];
        unsigned expected_first=0; std::vector<Oracle::Row> prefix;
        for(unsigned i=0;i<count;++i) {
            prefix.push_back(r.rows[i]);
            if(!expected_first && oracle.Rank(prefix)==K) expected_first=i+1;
        }
        {
            Owner decoder(api,r.counts); ++r.counts[2];
            Check(api.decoder(r.profile,B,decoder.handle)==0 && decoder.handle,"standalone decoder create");
            for(unsigned i=0;i<count;++i) {
                ++r.counts[3];
                const int status=api.feed(decoder.handle,ids[i],packets[a][i].data()+1,B);
                Check(status==0 || status==1,"feed status"); r.feed.push_back(status);
                if(status==0) { Check(i>=K-1,"premature success"); r.first=i+1; break; }
            }
            Check(r.first==expected_first,"every-arm first-success rank");
            if(r.first) for(unsigned repeat=0;repeat<2;++repeat) {
                std::vector<Byte> out(source.size()+2,0xa5); uint64_t written=0; ++r.counts[4];
                Check(api.recover(decoder.handle,out.data()+1,unsigned(source.size()),&written)==0 &&
                      written==source.size(),"recover status/length");
                Check(out.front()==0xa5 && out.back()==0xa5 &&
                      !memcmp(out.data()+1,source.data(),source.size()),"recovered payload/guards");
                ++r.recoveries;
            }
        }
        const Ledger expected={{K+1,(K+1)*count,1,r.first?r.first:count,r.first?2u:0u,K+2}};
        Check(r.counts==expected,"complete attempted API ledger");
        Check(source==original,"immutable source");
        for(unsigned i=0;i<count;++i)
            Check(packets[a][i].front()==0xa5 && packets[a][i].back()==0xa5 &&
                  Sha(packets[a][i].data()+1,B)==r.packets[i],"immutable packet");
    }
    for(unsigned a=0;a<3;++a) {
        const auto& x=result[a]; const auto& y=result[a+3];
        Check(x.profile==y.profile && x.packets==y.packets && x.rows==y.rows && x.feed==y.feed &&
              x.first==y.first && x.recoveries==y.recoveries && x.counts==y.counts,"policy-independent complete result");
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
    for(unsigned a=0;a<arms;++a) {
        if(a) putchar(',');
        const auto& r=result[a]; printf("{\"profile\":"); Hex(r.profile.data(),32); printf(",\"packets\":[");
        for(size_t i=0;i<r.packets.size();++i) { if(i) putchar(','); printf("\"%s\"",r.packets[i].c_str()); }
        printf("],\"rows\":\"");
        for(const auto& row:r.rows) for(Byte value:row) printf("%02x",unsigned(value));
        printf("\",\"feed\":[");
        for(size_t i=0;i<r.feed.size();++i) { if(i) putchar(','); printf("%d",r.feed[i]); }
        printf("],\"first\":%u,\"recoveries\":%u,\"counts\":[",r.first,r.recoveries);
        for(unsigned i=0;i<6;++i) { if(i) putchar(','); printf("%u",r.counts[i]); }
        printf("],\"checked\":true}");
    }
    printf("]}\n"); Flush();
}
gf256_x86_cpu_features Initialize(const Oracle& oracle) {
    Check(wirehair_init()==Wirehair_Success && GF256Ctx.Polynomial==0x14d,"shared GF initialization");
    Check(Sha(wh2_k5_data::kLookup,sizeof(wh2_k5_data::kLookup))==wh2_k5_data::kLookupSha,"sealed lookup");
    for(const auto& row:wh2_k5_data::kRows) {
        const auto expected=oracle.Coefficients(row.id);
        Check(std::equal(expected.begin(),expected.end(),row.values),"sealed polynomial coefficients");
    }
    gf256_x86_cpu_features f={}; gf256_get_active_x86_cpu_features(&f);
    if(WH2_K5_SERIALIZED_RECOVERY_BACKEND==1) Check(!f.SSSE3 && !f.AVX2 && !f.GFNI && !f.AVX512,"portable backend");
    else Check(f.SSSE3 && f.AVX2 && f.GFNI && f.AVX512,"frozen native dispatch");
    return f;
}
void Header(const std::string& claim,const char* scope,const gf256_x86_cpu_features& f) {
    printf("{\"type\":\"header\",\"protocol\":\"%s\",\"claim\":\"%s\",\"backend\":\"%s\",\"scope\":\"%s\","
           "\"retained_raw_sha256\":\"%s\",\"features\":[%u,%u,%u,%u],\"sources\":[",
           protocol,claim.c_str(),backends[WH2_K5_SERIALIZED_RECOVERY_BACKEND],scope,wh2_k5_data::kRawSha,
           unsigned(f.SSSE3),unsigned(f.AVX2),unsigned(f.GFNI),unsigned(f.AVX512));
    for(unsigned B:widths) { if(B!=2) putchar(','); const auto m=Message(B); printf("\"%s\"",Sha(m.data(),m.size()).c_str()); }
    printf("]}\n"); Flush();
}
void Footer(unsigned records) { printf("{\"type\":\"footer\",\"records\":%u,\"checked\":true}\n",records); Flush(); }
uint64_t Now() {
    timespec t={};
    Check(!clock_gettime(CLOCK_MONOTONIC,&t) && t.tv_sec>=0 && t.tv_nsec>=0 &&
          t.tv_nsec<1000000000 && uint64_t(t.tv_sec)<UINT64_MAX/1000000000-1,"monotonic clock");
    return uint64_t(t.tv_sec)*1000000000+uint64_t(t.tv_nsec);
}
void Authenticate(const char* path,const std::string& claim) {
    Check(claim.size()==64 && claim.find_first_not_of("0123456789abcdef")==std::string::npos,"claim hex");
    std::ifstream in(path,std::ios::binary); Check(bool(in),"claimed namespace");
    std::string receipt; char ch;
    while(in.get(ch)) { Check(receipt.size()<1024*1024,"claim cap"); receipt+=ch; }
    Check(in.eof() && Sha(receipt.data(),receipt.size())==claim,"claim bytes");
}
int Worker(const std::string& claim) {
    Authenticate(claim_path,claim);
    const rlimit cpu={cpu_seconds,cpu_seconds},memory={address_space_mib*1024u*1024u,address_space_mib*1024u*1024u},core={0,0};
    Check(!setrlimit(RLIMIT_CPU,&cpu) && !setrlimit(RLIMIT_CORE,&core),"worker limits");
    if(WH2_K5_SERIALIZED_RECOVERY_BACKEND!=2) Check(!setrlimit(RLIMIT_AS,&memory),"worker memory cap");
    const uint64_t start=Now(); const Oracle oracle; const auto f=Initialize(oracle); Header(claim,"retained",f);
    unsigned records=0;
    for(unsigned i=0;i<6216;++i) {
        const auto& t=wh2_k5_data::kTraces[i];
        Check(Now()-start<uint64_t(wall_seconds)*1000000000,"worker deadline");
        Record(i<6144?0:1,i<6144?i:i-6144,t.B,t.ids,9,Exercise(oracle,t.B,t.ids,9)); ++records;
    }
    for(unsigned i=0;i<54;++i) for(unsigned w=0;w<3;++w) {
        const auto& p=wh2_k5_data::kHistory[i]; if(!(p.widths & (1u<<w))) continue;
        Check(Now()-start<uint64_t(wall_seconds)*1000000000,"worker deadline");
        Record(2,i,widths[w],p.ids,p.count,Exercise(oracle,widths[w],p.ids,p.count)); ++records;
    }
    Check(records==6273 && Now()-start<uint64_t(wall_seconds)*1000000000,"complete roster/deadline");
    Footer(records); return 0;
}
unsigned fault_kind=0,fault_calls=0,fault_created=0,fault_freed=0;
int TrackedCreate(const void* s,unsigned b,Profile& p,void*& h) {
    const int r=SBorrowed(s,b,p,h); if(h) ++fault_created; return r;
}
int TrackedDecoder(const Profile& p,unsigned b,void*& h) {
    const int r=SDecoder(p,b,h); if(h) ++fault_created; return r;
}
void TrackedFree(void* h) { if(h) ++fault_freed; SFree(h); }
int FaultEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const int r=SEncode(h,id,out,b,n);
    if(++fault_calls==54) { if(fault_kind==1) throw std::runtime_error("neutral injected exception"); ++*n; }
    return r;
}
int FaultFeed(void* h,uint32_t id,const void* in,unsigned n) {
    const int r=SFeed(h,id,in,n); return ++fault_calls==5?-1:r;
}
int FaultRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const int r=SRecover(h,out,n,written); if(++fault_calls==2) ++*written; return r;
}
int Neutral(bool emit) {
    const Oracle oracle; const auto f=Initialize(oracle);
    for(unsigned a=0;a<256;++a) for(unsigned b=0;b<256;++b)
        Check(gf256_mul(Byte(a),Byte(b))==oracle.products[a][b],"exhaustive independent field");
    if(emit) Header(std::string(64,'0'),"neutral",f);
    const uint32_t streams[4][9]={{0,1,2,3,4,5,6,7,8},{5,6,7,8,9,10,11,12,13},
        {UINT32_MAX,UINT32_MAX-2,UINT32_MAX-4,UINT32_MAX-6,UINT32_MAX-8,UINT32_MAX-10,UINT32_MAX-12,UINT32_MAX-14,UINT32_MAX-16},
        {0,0,0,0,0,0,0,0,0}};
    unsigned cases=0;
    for(unsigned B:widths) for(unsigned s=0;s<4;++s) {
        const auto result=Exercise(oracle,B,streams[s],9);
        Check(result[1].first==(s==3?0u:K),"neutral endpoint");
        if(emit) Record(3,s,B,streams[s],9,result);
        ++cases;
    }
    for(unsigned block=0;block<2;++block) for(unsigned j=0;j<12;++j) {
        const unsigned i=block*6132+j; const auto& t=wh2_k5_data::kTraces[i];
        const auto result=Exercise(oracle,t.B,t.ids,9);
        if(emit) Record(4,i,t.B,t.ids,9,result);
        ++cases;
    }
    for(fault_kind=0;fault_kind<4;++fault_kind) {
        Api probe[arms]; std::copy(apis,apis+arms,probe);
        probe[4].create=TrackedCreate; probe[4].decoder=TrackedDecoder; probe[4].free=TrackedFree;
        if(fault_kind<2) probe[4].encode=FaultEncode;
        else if(fault_kind==2) probe[4].feed=FaultFeed;
        else probe[4].recover=FaultRecover;
        fault_calls=fault_created=fault_freed=0; bool failed=false;
        try { Exercise(oracle,2,streams[0],9,probe); } catch(const std::exception&) { failed=true; }
        Check(failed && fault_created==fault_freed && fault_created==(fault_kind<2?6u:7u),"late-call failure cleanup");
    }
    Check(cases==36,"neutral case count");
    if(emit) Footer(cases);
    else printf("PASS 36 neutral cases, six ownership-matched APIs, all-arm packet/rank oracle, four late-call cleanup checks\n");
    return 0;
}
static_assert(sizeof(wh2_k5_data::kTraces)/sizeof(wh2_k5_data::kTraces[0])==6216,"trace roster");
static_assert(sizeof(wh2_k5_data::kHistory)/sizeof(wh2_k5_data::kHistory[0])==54,"history roster");
static_assert(sizeof(wh2_k5_data::kWindows)/sizeof(wh2_k5_data::kWindows[0])==30,"window provenance");
static_assert(sizeof(wh2_k5_data::kRows)/sizeof(wh2_k5_data::kRows[0])==2270,"coefficient provenance");
} // namespace
int main(int argc,char** argv) {
    try {
        if(argc==2 && !strcmp(argv[1],"--contract")) {
            printf("{\"K\":%u,\"arms\":%u,\"records\":6273,\"cpu_seconds\":%u,\"wall_seconds\":%u,\"address_space_mib\":%u,\"asan_shadow_exempt\":%s,\"backend\":\"%s\",\"claim_path\":\"%s\"}\n",
                   K,arms,cpu_seconds,wall_seconds,address_space_mib,WH2_K5_SERIALIZED_RECOVERY_BACKEND==2?"true":"false",
                   backends[WH2_K5_SERIALIZED_RECOVERY_BACKEND],claim_path); return 0;
        }
        if(argc==4 && !strcmp(argv[1],"--neutral-claim")) {
            Authenticate(argv[2],argv[3]); puts("PASS claim authentication"); return 0;
        }
        if(argc==2 && !strcmp(argv[1],"--neutral")) return Neutral(false);
        if(argc==2 && !strcmp(argv[1],"--neutral-fixtures")) return Neutral(true);
        Check(argc==3 && !strcmp(argv[1],"--worker"),"explicit mode required"); return Worker(argv[2]);
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); return 1; }
}
