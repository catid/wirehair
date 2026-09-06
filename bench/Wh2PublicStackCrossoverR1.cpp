// Additive artifact-policy successor; the frozen R0 source and run remain immutable.
#include "Wh2InstalledCodeArtifact.h"
// One unchanged RTLD_LOCAL baseline; never invoke the included retired launcher.
#define main Wh2OffsetR0UnusedMain
#define offset_r0 stack_offset_private
#include "Wh2PublicOffsetR0.cpp"
#undef offset_r0
#undef main

namespace stack_crossover_r1 {
namespace old = stack_offset_private;
using old::A; using old::Address; using old::Api; using old::Arena;
using old::Cell; using old::Clock; using old::Encoder; using old::Fd;
using old::N; using old::O; using old::Q; using old::Require;
using old::Sha256Hex; using old::Shape;
static const char Name[] = "wh2-public-stack-crossover-r1";
static const char Prefix[] = "wirehair.wh2.public-stack-crossover-r1";
static const char Directory[] = "/var/tmp/wh2-public-stack-crossover-r1";
static const char Baseline[] = "/tmp/wh2-page-aa.uWaCGu/libwirehair.so.2.0.0";
static const char BaselineHash[] = "c0bf0666e3cc51523e18847b8a07384f2eac312518d0b4f9ac48cd14e63fa038";
static const char Libc[] = "/usr/lib/x86_64-linux-gnu/libc.so.6";
static const char LibcHash[] = "8db37cf3f2169f59a0f07ef1fea308c35656668c64c8ff294e1860f4121eb161";
static const uint32_t TargetOrders[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};

std::string CellKey(uint32_t rep, uint32_t cell)
{ return Sha256Hex("stack-crossover-r0:rep=" + N(rep) + ":cell=" + N(cell)); }
std::array<uint32_t,6> CellOrder(uint32_t rep)
{
    std::array<std::pair<std::string,uint32_t>,6> keyed;
    for (uint32_t cell=0;cell<6;++cell) keyed[cell]={CellKey(rep,cell),cell};
    std::sort(keyed.begin(),keyed.end());
    std::array<uint32_t,6> result;
    for (size_t i=0;i<6;++i) result[i]=keyed[i].second;
    return result;
}
std::string Describe()
{
    std::vector<std::string> sequences, conditions, cells, targets;
    for (const auto& row:old::Conditions) sequences.push_back(A({N(row[0]),N(row[1]),N(row[2]),N(row[3])}));
    for (uint32_t c=0;c<4;++c)
        conditions.push_back(O({{"condition",N(c)},{"mapping",Q(old::Mapping(c))},{"order",Q(old::Order(c))},
            {"warmup",A({Q(old::Side(c%2)),Q(old::Side((c%2)^1))})}}));
    for (const Shape& shape:old::Shapes) cells.push_back(O({{"K",N(shape.K)},{"block_bytes",N(shape.B)},
        {"scope",Q("prebuilt-systematic")},{"scope_invocations_per_batch",N(8192/shape.K)}}));
    for (const auto& row:TargetOrders) targets.push_back(A({N(row[0]),N(row[1]),N(row[2])}));
    return O({{"campaign",Q(Name)},{"cell_order",Q("sha256:stack-crossover-r0:rep=R:cell=I")},
        {"comparison_order",Q("(replicate+cell+index)%3")},
        {"comparisons",A({Q("h0/h0"),Q("h1/h1"),Q("h0/h1")})},
        {"condition_sequences",A(sequences)},{"conditions",A(conditions)},
        {"constructor_order",A({N(0),N(1),N(2),N(3)})},{"counter_phase",N(1024)},
        {"expected_encode_calls",N(212336640)},{"expected_measured_batches",N(20736)},
        {"expected_measured_encode_calls",N(169869312)},{"expected_measured_invocations",N(7527168)},
        {"expected_panels",N(2592)},{"expected_records",N(2594)},
        {"expected_warmup_batches",N(5184)},{"expected_warmup_encode_calls",N(42467328)},
        {"expected_warmup_invocations",N(1881792)},{"fixed_envelope_bytes",N(256)},
        {"internal_deadline_seconds",N(60)},{"moving_store_bytes",N(40)},
        {"output_phase",N(2336)},{"page_bytes",N(4096)},
        {"panel_replicates",N(12)},{"schema",Q(std::string(Prefix)+".describe.v1")},
        {"sibling_cpu",N(56)},{"source_phase",N(2048)},{"source_seed_base",Q("0xc199f24886210f53")},
        {"stack_limit_bytes",N(16384)},{"target_cpu",N(120)},
        {"target_order",Q("permutations[(replicate+cell)%6]")},{"target_orders",A(targets)},
        {"timing_granularity",Q("whole-batch")},{"tuples",A(cells)}})+'\n';
}

#if defined(__linux__) && defined(__x86_64__) && !defined(__ILP32__) && \
    defined(__GNUC__) && __GNUC__ == 13 && !defined(__clang__)
#define STACK_NOIPA __attribute__((noipa,noinline,noclone))
#define STACK_READ_RSP(value) __asm__ __volatile__("mov %%rsp, %0" : "=r"(value) : : "memory")
struct StackFacts {
    uintptr_t Pre=0, Hot=0, Final=0, Frame=0;
    size_t Correction=0;
    std::string Json() const
    {
        return O({{"correction",N(Correction)},{"final_rsp",Q(Address(Final))},
            {"frame_address",Q(Address(Frame))},{"hot_rsp",Q(Address(Hot))},{"pre_rsp",Q(Address(Pre))}});
    }
};
struct BatchResult {
    StackFacts Stack;
    uint64_t Nanoseconds=0;
    unsigned Error=0;
};
struct BatchRequest {
    decltype(&wirehair_v2_encode) Encode=nullptr;
    WirehairV2Codec Handle=nullptr;
    uint8_t* Output=nullptr;
    uint32_t* Counters=nullptr;
    uint32_t K=0, B=0, Traversals=0, Mode=0, Target=0;
    const StackFacts* Expected=nullptr;
};

// Mode 0 probes this exact frame without alloca, calls, or clocks.  Modes 1/2
// share a single public callsite, with mode 2 alone enabling whole-batch clocks.
// Frozen GCC rounding is D=(requested+23)&~15; request aligned D-16.
extern "C" STACK_NOIPA void StackCrossoverBatch(const BatchRequest* request, BatchResult* result) noexcept
{
    if (!result) return;
    result->Error=1;
    if (!request || request->Mode>2 || request->Target>4095 || (request->Target&15u)) return;
    uintptr_t pre;
    STACK_READ_RSP(pre);
    result->Stack.Pre=pre;
    result->Stack.Frame=reinterpret_cast<uintptr_t>(__builtin_frame_address(0));
    result->Stack.Hot=pre; result->Stack.Correction=0; result->Nanoseconds=0;
    if ((pre&15u) || pre>UINTPTR_MAX-255 || (result->Stack.Frame&15u) ||
        result->Stack.Frame<pre || result->Stack.Frame-pre>240) return;
    if (request->Mode==0) {
        uintptr_t after; STACK_READ_RSP(after);
        result->Stack.Final=after; result->Error=after!=pre;
        return;
    }
    if (!request->Expected || request->Expected->Pre!=pre || request->Expected->Frame!=result->Stack.Frame ||
        !request->Encode || !request->Handle || !request->Output || !request->Counters ||
        !request->K || request->K>8192 || !request->B || !request->Traversals ||
        request->Traversals>8192/request->K ||
        (request->Mode==2 && request->K*request->Traversals!=8192)) return;
    const size_t correction=4096+((pre-request->Target)&4095u);
    result->Stack.Correction=correction;
    if (result->Stack.Frame+16-pre+correction+40>16384) return;
    const size_t requested=correction-16;
    volatile uint8_t* padding=static_cast<volatile uint8_t*>(__builtin_alloca_with_align(requested,128));
    __asm__ __volatile__("" : : "r"(padding) : "memory");
    uintptr_t hot; STACK_READ_RSP(hot);
    result->Stack.Hot=hot;
    if (hot>pre || pre-hot!=correction || (hot&4095u)!=request->Target) return;
    for (size_t i=0;i<requested;i+=64) padding[i]=static_cast<uint8_t>(i);
    padding[requested-1]=0;
    const auto encode=request->Encode;
    const WirehairV2Codec handle=request->Handle;
    const uint32_t k=request->K, block=request->B, traversals=request->Traversals;
    uint8_t* const output=request->Output;
    uint32_t* const counters=request->Counters;
    const bool timed=request->Mode==2;
    Clock::time_point start;
    if (timed) start=Clock::now();
    uintptr_t before; STACK_READ_RSP(before);
    unsigned error=before!=hot;
    for (uint32_t traversal=0;traversal<traversals && !error;++traversal)
        for (uint32_t id=0;id<k;++id)
            if (encode(handle,id,output+static_cast<size_t>(id)*block,block,counters+id)!=WirehairV2_Success ||
                counters[id]!=block) { error=2; break; }
    uintptr_t after; STACK_READ_RSP(after);
    if (after!=hot) error=3;
    if (timed) {
        const int64_t elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-start).count();
        if (elapsed<=0) error=4;
        else result->Nanoseconds=static_cast<uint64_t>(elapsed);
    }
    // Receipt storage belongs to the caller, outside the fixed-frame envelope.
    // Publish it only after the stop clock, never inside the measured interval.
    result->Stack.Final=after;
    __asm__ __volatile__("" : : "r"(padding) : "memory");
    result->Error=error;
}

struct Invocation {
    const Cell* Data=nullptr;
    Arena* Memory=nullptr;
    const Encoder* Codec=nullptr;
    const StackFacts* Probe=nullptr;
    uint32_t Mode=0, Target=0;
};
// One callback callsite for probe/warm/timed modes; no extra stack arguments.
STACK_NOIPA BatchResult Invoke(const Invocation& call)
{
    Require(call.Data && call.Memory,"invocation lacks cell/arena");
    BatchRequest request;
    request.Mode=call.Mode; request.Target=call.Target; request.Expected=call.Probe;
    if (call.Mode) {
        Require(call.Codec && call.Codec->Owner && call.Codec->Handle,"invocation lacks encoder");
        call.Memory->Prepare(call.Data->Dimensions.K);
        request.Encode=call.Codec->Owner->Encode; request.Handle=call.Codec->Handle;
        request.Output=call.Memory->Output; request.Counters=call.Memory->Counters->Values;
        request.K=call.Data->Dimensions.K; request.B=call.Data->Dimensions.B;
        request.Traversals=8192/request.K;
    }
    BatchResult result;
    StackCrossoverBatch(&request,&result);
    Require(!result.Error,"batch stack geometry or public encode failed");
    if (call.Mode) call.Memory->Verify(*call.Data,false);
    return result;
}

bool Intersects(uintptr_t first, size_t bytes, uintptr_t other, size_t otherBytes)
{
    for (size_t i=0;i<bytes;++i) if (((first+i-other)&4095u)<otherBytes) return true;
    return false;
}
struct Geometry {
    std::array<uint32_t,2> Pair={{0,0}};
    std::array<uint32_t,3> Targets={{0,0,0}};
    StackFacts Probe;
    std::string Json() const
    {
        return O({{"fixed_envelope_bytes",N(256)},{"probe",Probe.Json()},
            {"selected_pair",A({N(Pair[0]),N(Pair[1])})},
            {"target_phases",A({N(Targets[0]),N(Targets[1]),N(Targets[2])})}});
    }
};
Geometry SelectGeometry(const std::array<uintptr_t,4>& handles, uintptr_t arena, size_t arenaBytes,
    const StackFacts& probe)
{
    Require(probe.Pre && !(probe.Pre&15u) && probe.Pre<=UINTPTR_MAX-255 && !(probe.Frame&15u) &&
        probe.Frame>=probe.Pre && probe.Frame-probe.Pre<=240 &&
        probe.Hot==probe.Pre && probe.Final==probe.Pre && !probe.Correction,"probe geometry differs");
    Require(arena && arenaBytes && arena<=UINTPTR_MAX-arenaBytes,"arena address overflows");
    std::set<uintptr_t> unique;
    for (uintptr_t h:handles) Require(h && !(h&15u) && h<=UINTPTR_MAX-0x110 &&
        (h<arena || h>=arena+arenaBytes) && unique.insert(h).second,"natural handles invalid");
    Geometry geometry; geometry.Probe=probe;
    bool found=false;
    for (uint32_t i=0;i<4 && !found;++i) for (uint32_t j=i+1;j<4 && !found;++j) {
        const uintptr_t wi=handles[i]+0xe8, wj=handles[j]+0xe8;
        if (Intersects(wi,40,wj,40) || Intersects(wi,40,probe.Pre,256) ||
            Intersects(wj,40,probe.Pre,256)) continue;
        geometry.Pair={{i,j}};
        geometry.Targets[0]=static_cast<uint32_t>((handles[i]+0x110)&4095u);
        geometry.Targets[1]=static_cast<uint32_t>((handles[j]+0x110)&4095u);
        found=true;
    }
    Require(found,"no qualifying natural handle pair");
    for (uint32_t target=0;target<4096;target+=16) {
        const uintptr_t start=uintptr_t(target)-40;
        if (Intersects(start,40,uintptr_t(geometry.Targets[0])-40,40) ||
            Intersects(start,40,uintptr_t(geometry.Targets[1])-40,40) ||
            Intersects(start,40,probe.Pre,256)) continue;
        geometry.Targets[2]=target;
        return geometry;
    }
    throw std::runtime_error("no qualifying null target");
}
uint32_t HandleIndex(const Geometry& geometry, uint32_t comparison, uint32_t physical)
{ return geometry.Pair[comparison==0?0:comparison==1?1:physical]; }
std::string Slot(const BatchResult& batch, uint32_t logical, uint32_t physical, uint32_t handle)
{
    return O({{"elapsed_ns",N(batch.Nanoseconds)},{"handle_index",N(handle)},
        {"logical_lane",Q(old::Side(logical))},{"physical_lane",N(physical)},{"stack",batch.Stack.Json()}});
}
std::string Terminal()
{
    return O({{"campaign",Q(Name)},{"encode_call_count",N(212336640)},{"measured_batch_count",N(20736)},
        {"measured_invocation_count",N(7527168)},{"panel_count",N(2592)},{"record_count",N(2594)},
        {"schema",Q(std::string(Prefix)+".terminal.v1")},{"status",Q("complete")},
        {"warmup_batch_count",N(5184)},{"warmup_invocation_count",N(1881792)}});
}
void Run(Api& baseline, const Clock::time_point& deadline)
{
    cpu_set_t affinity; CPU_ZERO(&affinity); CPU_SET(120,&affinity);
    Require(::sched_setaffinity(0,sizeof(affinity),&affinity)==0,"cannot pin worker");
    old::Guard(deadline); const std::string identity=old::Identity(), descriptionHash=Sha256Hex(Describe());
    Require(old::Field(baseline.Copy,"path")==Libc && old::Field(baseline.Copy,"elf_offset")=="0x1a14c0",
        "runtime memcpy implementation differs");
    std::vector<Cell> cells; std::vector<std::string> receipts;
    for (Shape shape:old::Shapes) {
        old::Guard(deadline); cells.push_back(old::BuildCell(baseline,baseline,shape)); receipts.push_back(cells.back().Json);
    }
    Require(old::Identity()==identity,"identity changed after oracles"); baseline.Stable(); old::Guard(deadline);
    old::Emit(O({{"bindings",O({{"baseline",baseline.Receipt}})},{"campaign",Q(Name)},
        {"cells",A(receipts)},{"description_sha256",Q(descriptionHash)},
        {"schema",Q(std::string(Prefix)+".config.v1")},{"target_identity",identity}}));
    uint32_t sequence=0;
    for (uint32_t rep=0;rep<12;++rep) for (uint32_t cellIndex:CellOrder(rep)) {
        old::Guard(deadline); const Cell& cell=cells[cellIndex]; Arena arena; arena.Initialize(cell);
        Encoder encoders[4]; std::array<uintptr_t,4> handles;
        for (uint32_t index=0;index<4;++index) {
            encoders[index].Create(baseline,arena.Input,arena.Bytes,cell.Dimensions.B);
            Require(Sha256Hex(encoders[index].Profile.data(),32)==cell.ProfileHash,"timed handle profile differs");
            handles[index]=reinterpret_cast<uintptr_t>(encoders[index].Handle);
        }
        Invocation call; call.Data=&cell; call.Memory=&arena;
        const BatchResult probe=Invoke(call);
        const Geometry geometry=SelectGeometry(handles,reinterpret_cast<uintptr_t>(arena.Base),arena.ArenaBytes,probe.Stack);
        call.Probe=&geometry.Probe;
        const std::string addresses=arena.Json(encoders), geometryJson=geometry.Json();
        for (uint32_t targetIndex:TargetOrders[(rep+cellIndex)%6]) {
            call.Target=geometry.Targets[targetIndex];
            for (uint32_t position=0;position<3;++position) {
                const uint32_t comparison=(rep+cellIndex+position)%3;
                for (uint32_t condition:old::Conditions[rep%4]) {
                    old::Guard(deadline); std::vector<std::string> warmup, slots;
                    call.Mode=1;
                    for (uint32_t warm=0;warm<2;++warm) {
                        const uint32_t logical=warm^(condition%2), physical=logical^(condition/2);
                        const uint32_t handle=HandleIndex(geometry,comparison,physical); call.Codec=&encoders[handle];
                        warmup.push_back(Slot(Invoke(call),logical,physical,handle));
                    }
                    call.Mode=2;
                    for (uint32_t logical:old::Slots(condition)) {
                        old::Guard(deadline); const uint32_t physical=logical^(condition/2);
                        const uint32_t handle=HandleIndex(geometry,comparison,physical); call.Codec=&encoders[handle];
                        slots.push_back(Slot(Invoke(call),logical,physical,handle));
                    }
                    arena.Verify(cell,true);
                    for (const Encoder& encoder:encoders)
                        Require(Sha256Hex(encoder.Profile.data(),32)==cell.ProfileHash,"panel profile changed");
                    old::Guard(deadline);
                    old::Emit(O({{"K",N(cell.Dimensions.K)},{"addresses",addresses},{"block_bytes",N(cell.Dimensions.B)},
                        {"campaign",Q(Name)},{"cell_index",N(cellIndex)},{"cell_key_sha256",Q(CellKey(rep,cellIndex))},
                        {"comparison_index",N(comparison)},{"condition",N(condition)},{"description_sha256",Q(descriptionHash)},
                        {"geometry",geometryJson},{"mapping",Q(old::Mapping(condition))},{"order",Q(old::Order(condition))},
                        {"output_sha256",Q(cell.SourceHash)},{"profile_sha256",Q(cell.ProfileHash)},{"replicate",N(rep)},
                        {"schema",Q(std::string(Prefix)+".panel.v1")},{"scope_invocations_per_batch",N(8192/cell.Dimensions.K)},
                        {"sequence",N(sequence++)},{"slots",A(slots)},{"source_immutable_verified","true"},
                        {"source_sha256",Q(cell.SourceHash)},{"target_index",N(targetIndex)},
                        {"target_phase",N(call.Target)},{"warmup",A(warmup)}}));
                }
            }
        }
    }
    Require(sequence==2592 && old::Identity()==identity,"terminal count or identity differs");
    baseline.Stable(); wirehair::wh2_benchmark::VerifyInstalledCodeArtifact(Libc,LibcHash);
    old::Guard(deadline); old::Emit(Terminal());
}
std::string Authorize(const std::string& directory, const std::string& claimHash,
    const std::string& workerHash, const std::string& descriptionHash)
{
    Require(directory==Directory && old::Hash(claimHash) && old::Hash(workerHash) &&
        descriptionHash==Sha256Hex(Describe()),"launch authorization differs");
    Fd dir(::open(directory.c_str(),O_RDONLY|O_DIRECTORY|O_NOFOLLOW|O_CLOEXEC));
    struct stat statbuf;
    Require(::fstat(dir.Value,&statbuf)==0 && (statbuf.st_mode&07777)==0700,"directory mode differs");
    Fd claimFd(::openat(dir.Value,"CLAIM",O_RDONLY|O_NOFOLLOW|O_CLOEXEC));
    const std::string claim=old::ReadFd(claimFd.Value,1024*1024,0400);
    Require(Sha256Hex(claim)==claimHash && old::Field(claim,"campaign")==Name &&
        old::Field(claim,"worker_sha256")==workerHash && old::Field(claim,"description_sha256")==descriptionHash &&
        old::FileHash("/proc/self/exe")==workerHash,"claim/executable hash differs");
    Require(old::Field(claim,"baseline_path")==Baseline && old::Field(claim,"baseline_sha256")==BaselineHash,
        "baseline authorization differs");
    wirehair::wh2_benchmark::VerifyInstalledCodeArtifact(Libc,LibcHash);
    Fd marker(::openat(dir.Value,"WORKER_STARTED",O_WRONLY|O_CREAT|O_EXCL|O_NOFOLLOW|O_CLOEXEC,0400));
    Require(::fchmod(marker.Value,0400)==0,"marker mode failed");
    old::WriteAll(marker.Value,O({{"campaign",Q(Name)},{"claim_sha256",Q(claimHash)},
        {"description_sha256",Q(descriptionHash)},{"schema",Q(std::string(Prefix)+".worker-started.v1")},
        {"worker_sha256",Q(workerHash)}})+'\n');
    Require(::fsync(marker.Value)==0 && ::fsync(dir.Value)==0,"marker sync failed");
    return claim;
}
struct MockContext { uintptr_t Entry=0; unsigned Calls=0, Error=0; };
extern "C" STACK_NOIPA WirehairV2Result StackCrossoverMock(
    WirehairV2Codec handle, uint32_t id, void* output, uint32_t capacity, uint32_t* written) noexcept
{
    uintptr_t entry; STACK_READ_RSP(entry);
    MockContext* context=reinterpret_cast<MockContext*>(handle);
    if (entry!=context->Entry || id!=context->Calls%8 || capacity!=64 || !output || !written) {
        context->Error=1; return WirehairV2_Error;
    }
    static_cast<uint8_t*>(output)[0]=static_cast<uint8_t>(id^0xa5u);
    *written=capacity; ++context->Calls;
    return WirehairV2_Success;
}
void Selftest()
{
    uint64_t panels=0, measured=0, warm=0;
    for (uint32_t rep=0;rep<12;++rep) {
        const auto cells=CellOrder(rep); Require(std::set<uint32_t>(cells.begin(),cells.end()).size()==6,"cell order differs");
        for (uint32_t cell:cells) for (uint32_t target:TargetOrders[(rep+cell)%6]) {
            Require(target<3,"target order differs");
            for (uint32_t comparison=0;comparison<3;++comparison) for (uint32_t condition:old::Conditions[rep%4]) {
                const auto slots=old::Slots(condition);
                Require(std::count(slots.begin(),slots.end(),0u)==4,"slot balance differs");
                ++panels; measured+=8*(8192/old::Shapes[cell].K); warm+=2*(8192/old::Shapes[cell].K);
            }
        }
    }
    Require(panels==2592 && measured==7527168 && warm==1881792,"schedule counts differ");
    StackFacts probe; probe.Pre=0x100800; probe.Frame=probe.Pre+128; probe.Hot=probe.Final=probe.Pre;
    const Geometry geometry=SelectGeometry({{0x200000,0x200100,0x200300,0x200800}},0x400000,8192,probe);
    Require(geometry.Pair[0]==0 && geometry.Pair[1]==1 && geometry.Targets[0]==0x110 &&
        geometry.Targets[1]==0x210 && geometry.Targets[2]==0,"first pair/null differs");
    Require(Intersects(4090,40,0,16) && !Intersects(4090,6,0,16),"modular byte overlap differs");
    Require(HandleIndex(geometry,0,0)==HandleIndex(geometry,0,1) &&
        HandleIndex(geometry,1,0)==HandleIndex(geometry,1,1) &&
        HandleIndex(geometry,2,0)!=HandleIndex(geometry,2,1),"comparison mapping differs");
    BatchRequest request;
    BatchResult captured;
    StackCrossoverBatch(&request,&captured);
    Require(!captured.Error && captured.Stack.Pre==captured.Stack.Hot &&
        captured.Stack.Pre==captured.Stack.Final && !captured.Stack.Correction,"callback probe failed");
    std::array<uint8_t,512> output={{0}};
    std::array<uint32_t,8> counters={{0}};
    MockContext mock;
    request.Encode=&StackCrossoverMock; request.Handle=reinterpret_cast<WirehairV2Codec>(&mock);
    request.Output=output.data(); request.Counters=counters.data();
    request.K=8; request.B=64; request.Traversals=4; request.Mode=1; request.Expected=&captured.Stack;
    for (unsigned sweep=0;sweep<2;++sweep) for (uint32_t target=0;target<4096;target+=16) {
        request.Target=target; mock.Calls=mock.Error=0;
        mock.Entry=captured.Stack.Pre-(4096+((captured.Stack.Pre-target)&4095u))-8;
        BatchResult batch; StackCrossoverBatch(&request,&batch);
        Require(!batch.Error && !mock.Error && mock.Calls==32 && !batch.Nanoseconds &&
            batch.Stack.Pre==captured.Stack.Pre && batch.Stack.Frame==captured.Stack.Frame &&
            batch.Stack.Final==batch.Stack.Hot && batch.Stack.Hot==mock.Entry+8,
            "mock callback geometry failed");
        for (uint32_t id=0;id<8;++id)
            Require(counters[id]==64 && output[id*64]==static_cast<uint8_t>(id^0xa5u),"mock ABI output differs");
    }
}
#endif
} // namespace stack_crossover_r1

int main(int argc, char** argv)
{
    using namespace stack_crossover_r1;
    try {
#if defined(__linux__) && defined(__x86_64__) && !defined(__ILP32__) && \
    defined(__GNUC__) && __GNUC__ == 13 && !defined(__clang__)
        if (argc==2 && std::string(argv[1])=="--describe") std::cout<<Describe();
        else if (argc==2 && std::string(argv[1])=="--artifact-preflight") {
            wirehair::wh2_benchmark::VerifyInstalledCodeArtifact(Libc,LibcHash);
            std::cout<<"stack-crossover-r1 installed libc verified (no codec or clocks)\n";
        }
        else if (argc==2 && std::string(argv[1])=="--selftest") {
            Selftest(); std::cout<<"stack-crossover-r1 pure selftest passed (no codec or clocks)\n";
        } else {
            Require(argc==12 && std::string(argv[1])=="--run" && std::string(argv[2])=="--cpu" &&
                std::string(argv[3])=="120" && std::string(argv[4])=="--evidence-dir" &&
                std::string(argv[6])=="--claim-sha256" && std::string(argv[8])=="--worker-sha256" &&
                std::string(argv[10])=="--config-identity-sha256","invalid arguments");
            const Clock::time_point deadline=Clock::now()+std::chrono::seconds(60);
            const std::string claim=Authorize(argv[5],argv[7],argv[9],argv[11]);
            Api baseline; baseline.Open(old::Field(claim,"baseline_path"),old::Field(claim,"baseline_sha256"));
            Run(baseline,deadline);
        }
        std::cout.flush(); Require(bool(std::cout),"cannot publish result"); return 0;
#else
        (void)argc; (void)argv;
        throw std::runtime_error("requires GCC 13/Linux x86-64 LP64");
#endif
    } catch (const std::exception& error) { std::cerr<<"stack-crossover-r1: "<<error.what()<<'\n'; return 1; }
}
