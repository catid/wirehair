// Current installed preserved paths, with output only after measurement.
// The historical worker and all of its defaults remain intact.
#include <csignal>
#define main wh2_historical_admission_main
#include "Wh2AdmissionRegressionCostR0.cpp"
#undef main

namespace {
void ReferenceHeaderJson(const std::string& claim,unsigned order,const std::string& identity,const Observation& prelude) {
    printf("{\"type\":\"header\",\"protocol\":\"%s\",\"claim\":\"%s\",\"load_order\":%u,\"identity_hex\":",protocol,claim.c_str(),order);
    Hex(identity.data(),identity.size()); printf(",\"prelude\":"); ObservationJson(prelude);
    printf(",\"bindings\":"); BindingsJson(); printf(",\"fixtures\":[");
    for(unsigned i=0;i<case_count;++i) { if(i) putchar(','); const auto& f=reference[i]; const auto& c=f.c;
        printf("{\"case\":[%u,%u,%u,%u],\"batch\":%u,\"source\":",c.family,c.k,c.b,c.policy,Batch(c));
        Hex(f.source.data(),f.source.size()); printf(",\"arms\":[");
        for(unsigned a=0;a<2;++a) { if(a) putchar(','); printf("{\"profile\":"); Hex(f.arm[a].profile.data(),32);
            printf(",\"packets\":"); Hex(f.arm[a].packets.data(),f.arm[a].packets.size());
            printf(",\"steps\":%u}",f.arm[a].steps); }
        printf("]}");
    }
    printf("]}\n"); Flush();
}

struct DeferredSink {
    const std::string& claim; unsigned order; const std::string& identity; const Observation& prelude;
    void Header() { ReferenceHeaderJson(claim,order,identity,prelude); }
    void Row(const Record& r) { RecordJson(r); }
    void Footer(bool failed,uint64_t work) { FooterJson(failed,work); }
};
template<class Sink> void PublishDeferred(Sink& sink,bool ended,bool failed,uint64_t work) {
    Check(ended && retained<=callbacks,"publication before measurement end");
    sink.Header();
    for(unsigned i=0;i<retained;++i) sink.Row(records[i]);
    sink.Footer(failed,work);
}
template<class R> void CaptureAndCheck(R& reader,Record& r,Observation& previous,uint64_t& work,
                                      const Api& api,bool corrupt_source=false) {
    const auto& c=r.c; auto& f=fixtures[c.which]; const auto& arm=f.arm[c.arm];
    Capture(reader,[&]{RunWork(api,f,arm,c.metric,outputs[c.order],r.work);},r.o);
    Validate(r.o,previous); Check(r.o.m1>=r.target,"early start"); previous=r.o;
    work+=r.o.m2-r.o.m1; Check(work<=UINT64_C(150000000000),"WORK cap");
    if(corrupt_source) fixtures[0].source[0]^=1;
    CheckWork(f,arm,c.metric,c.order,r.work); r.checked=true;
}
int DeferredWorker(const std::string& claim,unsigned order) {
    Check(!WH2_ADMISSION_REGRESSION_NEUTRAL,"neutral scientific worker disabled");
    ValidateClaim(claim,claim_path);
    const rlimit cpu={180,180},memory={512u*1024u*1024u,512u*1024u*1024u},core={0,0};
    Check(!setrlimit(RLIMIT_CPU,&cpu) && !setrlimit(RLIMIT_AS,&memory) && !setrlimit(RLIMIT_CORE,&core),"worker limits");
    Reader reader; const uint64_t start=reader.Mono(),cpu_start=reader.Cpu();
    std::string identity; Observation prelude,previous; uint64_t work=0; bool failed=false,initialized=false;
    try {
        Pin(); Initialize(order); identity=Identity(); memset(outputs,0xa5,sizeof(outputs));
        for(unsigned i=0;i<callbacks;++i) records[i].c=CoordinateAt(i);
        uint64_t final=0; Capture(reader,[&]{final=Prelude(UINT64_C(0x9e3779b97f4a7c15),1u<<20);},prelude);
        Check(final==UINT64_C(0x43935dad1647741b),"prelude checksum"); Validate(prelude,previous); previous=prelude;
        initialized=true;
        for(unsigned i=0;i<callbacks;++i) {
            Record& r=records[i]; retained=i+1; const auto& c=r.c; auto& f=fixtures[c.which]; Prepare(f);
            Check(reader.Mono()-start<UINT64_C(210000000000) && reader.Cpu()-cpu_start<UINT64_C(180000000000),"worker deadline");
            r.ready=reader.Mono(); Check(r.ready<=UINT64_MAX-c.q,"relative target overflow"); r.target=r.ready+c.q;
            r.wait[0]=reader.Mono(); r.wait[1]=reader.Cpu(); uint64_t now=r.wait[0];
            while(now<r.target) { now=reader.Mono(); Check(now-start<UINT64_C(210000000000),"wait deadline"); }
            r.wait[3]=reader.Cpu(); r.wait[2]=reader.Mono();
            CaptureAndCheck(reader,r,previous,work,f.api[c.arm]); Check(OnCpu(),"affinity changed");
        }
        for(unsigned a=0;a<2;++a) libraries[a]->Validate(a);
        Check(Identity()==identity,"target changed");
        Check(reader.Mono()-start<UINT64_C(210000000000) && reader.Cpu()-cpu_start<UINT64_C(180000000000),"final deadline");
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); failed=true; }
      catch(...) { fprintf(stderr,"INVALID: unknown measurement exception\n"); failed=true; }
    Check(std::signal(SIGPIPE,SIG_IGN)!=SIG_ERR,"output signal policy");
    if(!initialized) { FooterJson(true,work); return 1; }
    DeferredSink sink{claim,order,identity,prelude};
    PublishDeferred(sink,true,failed,work);
    Check(reader.Mono()-start<UINT64_C(210000000000) && reader.Cpu()-cpu_start<UINT64_C(180000000000),"publication deadline");
    return failed?1:0;
}
struct PublicationReader : FakeReader {
    bool fail=false; unsigned mono_calls=0;
    uint64_t Mono() {
        if(fail && ++mono_calls==3) throw std::runtime_error("neutral final clock");
        return FakeReader::Mono();
    }
};
struct CountingSink {
    unsigned calls=0,fail_at=0;
    void Step() { Check(++calls!=fail_at,"neutral sink failure"); }
    void Header() { Step(); } void Row(const Record&) { Step(); } void Footer(bool,uint64_t) { Step(); }
};
int PublicationNeutral(unsigned order,const char* mode) {
    const bool success=!strcmp(mode,"success"),bad=!strcmp(mode,"last-recover"),
        throws=!strcmp(mode,"throw-recover"),clock=!strcmp(mode,"last-clock"),source=!strcmp(mode,"last-source");
    Check(success || bad || throws || clock || source,"neutral publication mode");
    Pin(); Initialize(order); const std::string identity=Identity();
    PublicationReader reader; Observation prelude,previous;
    Capture(reader,[]{},prelude); Validate(prelude,previous); previous=prelude;
    uint64_t work=0; bool failed=false;
    for(unsigned i=0;i<3;++i) {
        Record& r=records[i]; r.c=CoordinateAt(54+i); retained=i+1;
        auto& f=fixtures[r.c.which]; const unsigned a=r.c.arm;
        Check(r.c.which==0 && r.c.metric==1 && a==1,"neutral uses exact original decoder coordinates");
        Prepare(f); r.ready=previous.m3; r.target=r.ready; r.wait={{r.ready,previous.c1,r.ready,previous.c1}};
        Symbols probe=symbols[a]; Api api=f.api[a]; api.s=&probe;
        if(i==2) {
            injected=&symbols[a]; injected_calls=injected_frees=0; injected_at=Batch(f.c); injected_throw=throws;
            if(bad || throws) probe.wirehair_v2_recover=InjectRecover;
            probe.wirehair_v2_free=InjectFree; reader.fail=clock; reader.mono_calls=0;
        }
        try { CaptureAndCheck(reader,r,previous,work,api,i==2 && source); }
        catch(const std::exception&) { failed=true; }
        if(failed) break;
    }
    Check(retained==3 && failed==!success,"neutral terminal state");
    Check(injected_frees==128 && ((!bad && !throws) || injected_calls==128),"neutral final handle cleanup");
    Check(records[0].checked && records[1].checked && records[2].checked==success,"neutral full prefix retention");
    Check(records[2].work.complete==!(bad || throws),"neutral completed work retained");
    Check((records[2].o.m2==0)==clock,"neutral partial clock retained");
    CountingSink early; bool rejected=false;
    try { PublishDeferred(early,false,failed,work); } catch(const std::exception&) { rejected=true; }
    Check(rejected && early.calls==0,"early publication forbidden");
    for(unsigned at=1;at<=5;++at) {
        CountingSink broken; broken.fail_at=at; rejected=false;
        try { PublishDeferred(broken,true,failed,work); } catch(const std::exception&) { rejected=true; }
        Check(rejected && broken.calls==at,"output failure never retried");
    }
    Check(std::signal(SIGPIPE,SIG_IGN)!=SIG_ERR,"output signal policy");
    const std::string neutral_claim(64,'0');
    DeferredSink stable{neutral_claim,order,identity,prelude};
    PublishDeferred(stable,true,failed,work);
    return 0;
}
} // namespace
int main(int argc,char** argv) {
    try {
        for(const char* name:{"MALLOC_TRIM_THRESHOLD_","MALLOC_MMAP_THRESHOLD_","MALLOC_TOP_PAD_","MALLOC_PERTURB_",
                             "GLIBC_TUNABLES","LD_PRELOAD","LD_LIBRARY_PATH","LD_AUDIT","LD_DEBUG"})
            Check(!getenv(name),"clean allocator/loader environment");
        if(argc==4 && !strcmp(argv[1],"--neutral-publication")) return PublicationNeutral(Order(argv[2]),argv[3]);
        if(argc==4 && !strcmp(argv[1],"--worker")) return DeferredWorker(argv[2],Order(argv[3]));
        return wh2_historical_admission_main(argc,argv);
    } catch(const std::exception& e) { fprintf(stderr,"INVALID: %s\n",e.what()); return 1; }
}
