// Common future small-codec retained recovery worker. Configuration is compile-time.
// Adapted from the immutable K5 worker; old experiment sources remain unchanged.
#include "wirehair/wirehair.h"
#include "Wh2SmallSerialized.h"
#include "Wh2FrozenTrace.h"
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

#if !defined(WH2_SMALL_RECOVERY_BACKEND) || defined(WH_COUNT) || defined(WIREHAIR_TESTING) || defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "Explicit unmodified recovery backend required"
#endif
static_assert(WH2_SMALL_CODEC_K==small_recovery_config::K, "matching sealed boundary");
static_assert(WH2_SMALL_RECOVERY_BACKEND>=0 && WH2_SMALL_RECOVERY_BACKEND<=2, "backend");
namespace {
using Byte=uint8_t;
using Profile=std::array<Byte,32>;
const unsigned K=small_recovery_config::K, arms=6, cpu_seconds=180, wall_seconds=210, address_space_mib=384;
const char* const protocol=small_recovery_config::protocol;
const char* const claim_path=small_recovery_config::claim_path;
const char* const backends[]={"native","scalar","asan"};
const unsigned widths[]={2,64,1280};
void Check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
std::string Sha(const void* data,size_t size) { return wirehair::wh2_benchmark::Sha256Hex(data,size); }
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
    const auto r=wh2_small_encode(h,id,out,b);
    if(r.bytes_written>b || (r.status==Wh2Small_Success &&
        (!r.bytes_written || r.bytes_written!=r.bytes_required))) return -1;
    *n=static_cast<uint32_t>(r.bytes_written); return static_cast<int>(r.status);
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
        const auto& feedback=small_recovery_config::feedback;
        for(unsigned phase=0;phase<2;++phase) {
            Matrix m={};
            for(unsigned i=0;i<K-1;++i) m[K*(i+1)+i]=1;
            for(unsigned i=0;i<K;++i) m[K*i+K-1]=Byte(feedback[i]^(i==0?small_recovery_config::lambda*phase:0));
            powers[phase][0]=m;
        }
        for(unsigned bit=1;bit<32;++bit) {
            powers[0][bit]=Product(powers[0][bit-1],powers[1][bit-1]);
            powers[1][bit]=Product(powers[1][bit-1],powers[0][bit-1]);
        }
    }
    Row Coefficients(uint32_t id) const {
        Row row={}; row[0]=1;
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
            if(k*B+j<source.size()) packet[j]^=products[row[k]][source[k*B+j]];
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
std::vector<Byte> Message(unsigned B,unsigned tail) {
    Check((B==2 || B==64 || B==1280) && tail>=1 && tail<=B,"bounded message");
    std::vector<Byte> source((K-1)*B+tail);
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
unsigned PacketBytes(uint32_t id,unsigned B,unsigned tail) { return id==K-1?tail:B; }
void PacketGuards(const std::vector<Byte>& packet,unsigned bytes) {
    Check(packet.front()==0xa5 && bytes+1<packet.size(),"packet leading guard");
    Check(std::all_of(packet.begin()+bytes+1,packet.end(),[](Byte b){return b==0xa5;}),"packet unused capacity guards");
}
Result Exercise(const Oracle& oracle,unsigned B,unsigned tail,const uint32_t* ids,unsigned count,const Api* api_set=apis) {
    Check((B==2 || B==64 || B==1280) && count>=K && count<=K+4,"bounded retained shape");
    const auto source=Message(B,tail);
    Result result;
    std::array<std::vector<std::vector<Byte>>,arms> packets;
    for(unsigned a=0;a<arms;++a) {
        auto& r=result[a]; const auto& api=api_set[a];
        packets[a].assign(count,std::vector<Byte>(B+2,0xa5)); r.rows.resize(count);
        for(unsigned basis=0;basis<=K;++basis) {
            std::vector<Byte> input=basis?std::vector<Byte>(source.size(),0):source;
            if(basis) input[(basis-1)*B]=1;
            const auto input_copy=input;
            Owner encoder(api,r.counts); Profile profile={};
            ++r.counts[0];
            Check(api.create(input.data(),unsigned(source.size()),B,profile,encoder.handle)==0 && encoder.handle,"encoder create");
            if(!basis) r.profile=profile; else Check(profile==r.profile,"basis profile identity");
            if(a<3) { std::fill(input.begin(),input.end(),0xcc); std::vector<Byte>().swap(input); }
            for(unsigned i=0;i<count;++i) {
                std::vector<Byte> packet(B+2,0xa5); uint32_t written=0; ++r.counts[1];
                const unsigned bytes=PacketBytes(ids[i],B,tail);
                Check(api.encode(encoder.handle,ids[i],packet.data()+1,B,&written)==0 && written==bytes,"encode status/length");
                PacketGuards(packet,bytes);
                if(a>=3) Check(input==input_copy,"actual borrowed input immutable");
                if(basis) {
                    r.rows[i][basis-1]=packet[1];
                    for(unsigned j=1;j<bytes;++j) Check(packet[j+1]==0,"basis other bytes");
                } else {
                    r.packets.push_back(Sha(packet.data()+1,bytes));
                    packets[a][i]=std::move(packet);
                }
            }
        }
        if(a==1 || a==4) Check(wh2_small_profile_validate(r.profile.data(),32)==Wh2Small_Success,"sealed candidate descriptor");
        if(a==0 || a==3) {
            WirehairV2Profile p={};
            Check(wirehair_v2_profile_deserialize(r.profile.data(),32,&p)==WirehairV2_Success &&
                  p.profile_id==WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,"ordinary selects certified profile");
        }
        for(unsigned i=0;i<count;++i) {
            if(a==1 || a==4) Check(r.rows[i]==oracle.Coefficients(ids[i]),"independent sealed coefficients");
            const auto expected=oracle.Packet(source,B,r.rows[i]);
            Check(!memcmp(packets[a][i].data()+1,expected.data(),PacketBytes(ids[i],B,tail)),"independent every-arm payload");
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
            Check(api.decoder(r.profile,unsigned(source.size()),B,decoder.handle)==0 && decoder.handle,"standalone decoder create");
            for(unsigned i=0;i<count;++i) {
                ++r.counts[3];
                const int status=api.feed(decoder.handle,ids[i],packets[a][i].data()+1,PacketBytes(ids[i],B,tail));
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
        for(unsigned i=0;i<count;++i) {
            const unsigned bytes=PacketBytes(ids[i],B,tail);
            PacketGuards(packets[a][i],bytes);
            Check(Sha(packets[a][i].data()+1,bytes)==r.packets[i],"immutable packet");
        }
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
void Record(unsigned group,unsigned index,unsigned B,unsigned tail,const uint32_t* ids,unsigned count,const Result& result) {
    printf("{\"type\":\"record\",\"group\":%u,\"index\":%u,\"B\":%u,\"tail\":%u,\"ids\":[",group,index,B,tail);
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
    Check(Sha(small_recovery_data::kLookup,sizeof(small_recovery_data::kLookup))==small_recovery_data::kLookupSha,"sealed lookup");
    for(const auto& row:small_recovery_data::kRows) {
        const auto expected=oracle.Coefficients(row.id);
        Check(std::equal(expected.begin(),expected.end(),row.values),"sealed polynomial coefficients");
    }
    gf256_x86_cpu_features f={}; gf256_get_active_x86_cpu_features(&f);
    if(WH2_SMALL_RECOVERY_BACKEND==1) Check(!f.SSSE3 && !f.AVX2 && !f.GFNI && !f.AVX512,"portable backend");
    else Check(f.SSSE3 && f.AVX2 && f.GFNI && f.AVX512,"frozen native dispatch");
    return f;
}
void Header(const std::string& claim,const char* scope,const gf256_x86_cpu_features& f) {
    printf("{\"type\":\"header\",\"protocol\":\"%s\",\"claim\":\"%s\",\"backend\":\"%s\",\"scope\":\"%s\","
           "\"retained_raw_sha256\":\"%s\",\"features\":[%u,%u,%u,%u],\"sources\":[",
           protocol,claim.c_str(),backends[WH2_SMALL_RECOVERY_BACKEND],scope,small_recovery_data::kRawSha,
           unsigned(f.SSSE3),unsigned(f.AVX2),unsigned(f.GFNI),unsigned(f.AVX512));
    for(unsigned B:widths) { if(B!=2) putchar(','); const auto m=Message(B,B); printf("\"%s\"",Sha(m.data(),m.size()).c_str()); }
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
    if(WH2_SMALL_RECOVERY_BACKEND!=2) Check(!setrlimit(RLIMIT_AS,&memory),"worker memory cap");
    const uint64_t start=Now(); const Oracle oracle; const auto f=Initialize(oracle); Header(claim,"retained",f);
    unsigned records=0;
    unsigned packet_ids=0;
    for(unsigned i=0;i<6216;++i) {
        const auto& t=small_recovery_data::kTraces[i];
        Check(Now()-start<uint64_t(wall_seconds)*1000000000,"worker deadline");
        const unsigned count=K+4;
        Record(i<6144?0:1,i<6144?i:i-6144,t.B,t.B,t.ids,count,Exercise(oracle,t.B,t.B,t.ids,count));
        ++records; packet_ids+=count;
    }
    for(unsigned i=0;i<small_recovery_config::history_count;++i) {
        const auto& p=small_recovery_data::kHistory[i];
        Check(Now()-start<uint64_t(wall_seconds)*1000000000,"worker deadline");
        Record(2,i,p.B,p.tail,p.ids,p.count,Exercise(oracle,p.B,p.tail,p.ids,p.count));
        ++records; packet_ids+=p.count;
    }
    Check(records==small_recovery_config::records && packet_ids==small_recovery_config::packet_ids &&
          Now()-start<uint64_t(wall_seconds)*1000000000,"complete roster/deadline");
    Footer(records); return 0;
}
unsigned fault_kind=0,fault_calls=0,fault_created=0,fault_freed=0;
int TrackedCreate(const void* s,unsigned message,unsigned b,Profile& p,void*& h) {
    const int r=SBorrowed(s,message,b,p,h); if(h) ++fault_created; return r;
}
int TrackedDecoder(const Profile& p,unsigned message,unsigned b,void*& h) {
    const int r=SDecoder(p,message,b,h); if(h) ++fault_created; return r;
}
void TrackedFree(void* h) { if(h) ++fault_freed; SFree(h); }
int FaultEncode(void* h,uint32_t id,void* out,unsigned b,uint32_t* n) {
    const int r=SEncode(h,id,out,b,n);
    if(++fault_calls==(K+1)*(K+4)) {
        if(fault_kind==1) throw std::runtime_error("neutral injected exception");
        if(fault_kind==2) static_cast<Byte*>(out)[*n]=0;
        else ++*n;
    }
    return r;
}
int FaultFeed(void* h,uint32_t id,const void* in,unsigned n) {
    const int r=SFeed(h,id,in,n); return ++fault_calls==K?-1:r;
}
int FaultRecover(void* h,void* out,unsigned n,uint64_t* written) {
    const int r=SRecover(h,out,n,written);
    if(++fault_calls==2) { if(fault_kind==5) static_cast<Byte*>(out)[n]=0; else ++*written; }
    return r;
}
int Neutral(bool emit) {
    const Oracle oracle; const auto f=Initialize(oracle);
    for(unsigned a=0;a<256;++a) for(unsigned b=0;b<256;++b)
        Check(gf256_mul(Byte(a),Byte(b))==oracle.products[a][b],"exhaustive independent field");
    if(emit) Header(std::string(64,'0'),"neutral",f);
    uint32_t streams[4][K+4]={};
    for(unsigned j=0;j<K+4;++j) {
        streams[0][j]=j; streams[1][j]=K+j; streams[2][j]=UINT32_MAX-2*j;
    }
    unsigned cases=0;
    for(unsigned B:widths) for(unsigned shape=0;shape<2;++shape) for(unsigned s=0;s<4;++s) {
        const unsigned tail=shape?1:B;
        const auto result=Exercise(oracle,B,tail,streams[s],K+4);
        if(s==0) Check(result[1].first==K,"systematic neutral endpoint");
        if(s==3) Check(result[1].first==0,"deficient neutral endpoint");
        if(emit) Record(3,shape*4+s,B,tail,streams[s],K+4,result);
        ++cases;
    }
    for(unsigned block=0;block<2;++block) for(unsigned j=0;j<12;++j) {
        const unsigned i=block*6132+j; const auto& t=small_recovery_data::kTraces[i];
        const auto result=Exercise(oracle,t.B,t.B,t.ids,K+4);
        if(emit) Record(4,i,t.B,t.B,t.ids,K+4,result);
        ++cases;
    }
    for(unsigned shape=0;shape<2;++shape) for(fault_kind=0;fault_kind<6;++fault_kind) {
        Api probe[arms]; std::copy(apis,apis+arms,probe);
        probe[4].create=TrackedCreate; probe[4].decoder=TrackedDecoder; probe[4].free=TrackedFree;
        if(fault_kind<3) probe[4].encode=FaultEncode;
        else if(fault_kind==3) probe[4].feed=FaultFeed;
        else probe[4].recover=FaultRecover;
        fault_calls=fault_created=fault_freed=0; bool failed=false;
        // Last encode is ID K-1, so tail1 checks the un-emitted capacity too.
        uint32_t fault_ids[K+4];
        for(unsigned i=0;i<K+4;++i) fault_ids[i]=streams[0][i];
        std::swap(fault_ids[K-1],fault_ids[K+3]);
        try { Exercise(oracle,2,shape?1:2,fault_kind<3?fault_ids:streams[0],K+4,probe); }
        catch(const std::exception&) { failed=true; }
        Check(failed && fault_created==fault_freed && fault_created==(fault_kind<3?K+1:K+2) &&
              fault_calls==(fault_kind<3?(K+1)*(K+4):(fault_kind==3?K:2)),"late-call failure cleanup");
    }
    Check(cases==48,"neutral case count");
    if(emit) Footer(cases);
    else printf("PASS 48 neutral cases, six ownership-matched APIs, all-arm packet/rank oracle, twelve late-call cleanup checks\n");
    return 0;
}
static_assert(sizeof(small_recovery_data::kTraces)/sizeof(small_recovery_data::kTraces[0])==6216,"trace roster");
static_assert(sizeof(small_recovery_data::kHistory)/sizeof(small_recovery_data::kHistory[0])==small_recovery_config::history_count,"history roster");
// K4 additionally retains the exponent-two window; the existing K8 roster
// stays at 30. These are provenance checks, not extra recovery records.
static_assert(sizeof(small_recovery_data::kWindows)/sizeof(small_recovery_data::kWindows[0])==
              (small_recovery_config::K==4?31u:30u),"window provenance");
static_assert(sizeof(small_recovery_data::kRows)/sizeof(small_recovery_data::kRows[0])==small_recovery_config::row_count,"coefficient provenance");
} // namespace
int main(int argc,char** argv) {
    try {
        if(argc==2 && !strcmp(argv[1],"--contract")) {
            printf("{\"K\":%u,\"arms\":%u,\"records\":%u,\"cpu_seconds\":%u,\"wall_seconds\":%u,\"address_space_mib\":%u,\"asan_shadow_exempt\":%s,\"backend\":\"%s\",\"claim_path\":\"%s\"}\n",
                   K,arms,small_recovery_config::records,cpu_seconds,wall_seconds,address_space_mib,WH2_SMALL_RECOVERY_BACKEND==2?"true":"false",
                   backends[WH2_SMALL_RECOVERY_BACKEND],claim_path); return 0;
        }
        if(argc==4 && !strcmp(argv[1],"--neutral-claim")) {
            Authenticate(argv[2],argv[3]); puts("PASS claim authentication"); return 0;
        }
        if(argc==2 && !strcmp(argv[1],"--neutral")) return Neutral(false);
        if(argc==2 && !strcmp(argv[1],"--neutral-fixtures")) return Neutral(true);
        Check(argc==3 && !strcmp(argv[1],"--worker"),"explicit mode required"); return Worker(argv[2]);
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); return 1; }
}
