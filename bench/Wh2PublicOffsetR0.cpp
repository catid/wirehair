// Standalone dual-RTLD_LOCAL public facade screen: NEVER link libwirehair here.
#include "Wh2FrozenTrace.h"
#include "Wh2PublicBorrowedTargetIdentity.h"
#include "Wh2RdpruTargetIdentityV2.h"
#include <wirehair/wirehair.h>
#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <new>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__linux__)
#include <dlfcn.h>
#include <fcntl.h>
#include <sched.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace offset_r0 {
using wirehair::wh2_benchmark::Sha256Hex;
using Clock = std::chrono::steady_clock;
static const char Name[] = "wh2-public-offset-r0";
static const char Prefix[] = "wirehair.wh2.public-offset-r0";
static const char Directory[] = "/var/tmp/wh2-public-offset-r0";
static const uint64_t Seed = UINT64_C(0xc199f24886210f53);
struct Shape { uint32_t K, B; };
static const Shape Shapes[] = {{8,64},{8,1280},{128,64},{128,1280},{8192,64},{8192,1280}};
static const uint32_t Conditions[4][4] = {{0,1,3,2},{1,2,0,3},{2,3,1,0},{3,0,2,1}};
static const uint32_t Pairs[4][2] = {{0,1},{2,3},{2,0},{3,1}};
void Require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
std::string Q(const std::string& text)
{
    std::string out = "\"";
    for (unsigned char c : text) {
        Require(c >= 32 && c < 127, "non-ASCII JSON text");
        if (c == '"' || c == '\\') out += '\\';
        out += static_cast<char>(c);
    }
    return out + '"';
}
std::string N(uint64_t number) { return std::to_string(number); }
std::string A(const std::vector<std::string>& items)
{
    std::string out = "[";
    for (size_t i = 0; i < items.size(); ++i) { if (i) out += ','; out += items[i]; }
    return out + ']';
}
std::string O(const std::map<std::string, std::string>& items)
{
    std::string out = "{";
    for (const auto& item : items) { if (out.size() > 1) out += ','; out += Q(item.first) + ':' + item.second; }
    return out + '}';
}
std::string Address(uintptr_t address) { std::ostringstream s; s << "0x" << std::hex << address; return s.str(); }
std::string Hex(const void* bytes, size_t size)
{
    static const char digits[] = "0123456789abcdef";
    const uint8_t* data = static_cast<const uint8_t*>(bytes);
    std::string out;
    out.reserve(size * 2);
    for (size_t i = 0; i < size; ++i) { out += digits[data[i] >> 4]; out += digits[data[i] & 15]; }
    return out;
}
bool Hash(const std::string& hash)
{ return hash.size() == 64 && hash.find_first_not_of("0123456789abcdef") == std::string::npos; }
const char* Order(uint32_t condition) { return condition % 2 ? "BAAB" : "ABBA"; }
const char* Mapping(uint32_t condition) { return condition < 2 ? "direct" : "swapped"; }
const char* Side(uint32_t side) { return side ? "right" : "left"; }
std::array<uint32_t,8> Slots(uint32_t condition)
{
    std::array<uint32_t,8> slots = {{0,1,1,0,1,0,0,1}};
    for (auto& side : slots) side ^= condition % 2;
    return slots;
}
std::string CellKey(uint32_t rep, uint32_t cell)
{ return Sha256Hex("offset-r0:rep=" + N(rep) + ":cell=" + N(cell)); }
std::array<uint32_t,6> CellOrder(uint32_t rep)
{
    std::array<std::pair<std::string,uint32_t>,6> keyed;
    for (uint32_t cell = 0; cell < 6; ++cell) keyed[cell] = {CellKey(rep,cell),cell};
    std::sort(keyed.begin(),keyed.end());
    std::array<uint32_t,6> result;
    for (size_t i = 0; i < 6; ++i) result[i] = keyed[i].second;
    return result;
}
std::string Describe()
{
    std::vector<std::string> sequences, conditions, cells;
    for (const auto& row : Conditions) sequences.push_back(A({N(row[0]),N(row[1]),N(row[2]),N(row[3])}));
    for (uint32_t c = 0; c < 4; ++c)
        conditions.push_back(O({{"condition",N(c)},{"mapping",Q(Mapping(c))},{"order",Q(Order(c))},
            {"warmup",A({Q(Side(c%2)),Q(Side((c%2)^1))})}}));
    for (const auto& shape : Shapes) cells.push_back(O({{"K",N(shape.K)},{"block_bytes",N(shape.B)},
        {"scope",Q("prebuilt-systematic")},{"scope_invocations_per_batch",N(8192/shape.K)}}));
    return O({{"campaign",Q(Name)},{"cell_order",Q("sha256:offset-r0:rep=R:cell=I")},
        {"comparison_order",Q("(replicate+cell+index)%4")},
        {"comparisons",A({Q("B0/B1"),Q("C0/C1"),Q("C0/B0"),Q("C1/B1")})},
        {"condition_sequences",A(sequences)},{"conditions",A(conditions)},
        {"constructor_orders",A({A({Q("B0"),Q("C0"),Q("C1"),Q("B1")}),A({Q("C0"),Q("B0"),Q("B1"),Q("C1")})})},
        {"counter_phase",N(1024)},{"expected_encode_calls",N(94371840)},
        {"expected_measured_batches",N(9216)},{"expected_measured_encode_calls",N(75497472)},
        {"expected_measured_invocations",N(3345408)},{"expected_panels",N(1152)},{"expected_records",N(1154)},
        {"expected_warmup_batches",N(2304)},{"expected_warmup_encode_calls",N(18874368)},
        {"expected_warmup_invocations",N(836352)},{"internal_deadline_seconds",N(60)},
        {"output_phase",N(2336)},{"page_bytes",N(4096)},{"panel_replicates",N(12)},
        {"schema",Q(std::string(Prefix)+".describe.v1")},{"sibling_cpu",N(56)},
        {"source_phase",N(2048)},{"source_seed_base",Q("0xc199f24886210f53")},
        {"target_cpu",N(120)},{"timing_granularity",Q("whole-batch")},{"tuples",A(cells)}}) + '\n';
}
uint64_t SplitMix(uint64_t& state)
{
    uint64_t word = (state += UINT64_C(0x9e3779b97f4a7c15));
    word = (word ^ (word >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    word = (word ^ (word >> 27)) * UINT64_C(0x94d049bb133111eb);
    return word ^ (word >> 31);
}
std::vector<uint8_t> Source(Shape shape, uint32_t tail)
{
    std::vector<uint8_t> out(static_cast<size_t>(shape.K-1)*shape.B+tail);
    uint64_t state = Seed ^ (uint64_t(shape.K)<<33) ^ (uint64_t(shape.B)<<1) ^ tail;
    for (size_t i = 0; i < out.size();) {
        uint64_t word = SplitMix(state);
        for (unsigned byte = 0; byte < 8 && i < out.size(); ++byte,++i)
            out[i] = static_cast<uint8_t>(word >> (byte*8));
    }
    return out;
}

#if defined(__linux__)
std::string RealPath(const std::string& path)
{
    char* resolved = ::realpath(path.c_str(),nullptr);
    Require(resolved != nullptr,"cannot resolve path");
    std::string result(resolved); std::free(resolved); return result;
}
struct Fd
{
    int Value;
    explicit Fd(int fd) : Value(fd) { Require(fd >= 0,"cannot open file"); }
    ~Fd() { ::close(Value); }
    Fd(const Fd&) = delete; Fd& operator=(const Fd&) = delete;
};
std::string ReadFd(int fd, size_t limit, unsigned mode = 0)
{
    struct stat statbuf;
    Require(::fstat(fd,&statbuf)==0 && S_ISREG(statbuf.st_mode) && statbuf.st_nlink==1 &&
        (!mode || (statbuf.st_mode & 07777)==mode) && statbuf.st_size>=0 &&
        static_cast<uint64_t>(statbuf.st_size)<=limit,"file metadata differs");
    std::string result; char bytes[16384];
    for (;;) {
        const ssize_t n = ::read(fd,bytes,sizeof(bytes));
        if (n < 0 && errno == EINTR) continue;
        Require(n >= 0,"file read failed");
        if (!n) break;
        Require(result.size()+static_cast<size_t>(n)<=limit,"file exceeds cap");
        result.append(bytes,static_cast<size_t>(n));
    }
    Require(result.size()==static_cast<size_t>(statbuf.st_size),"file changed while reading");
    return result;
}
std::string FileHash(const std::string& path)
{
    Fd fd(::open(path.c_str(),O_RDONLY|O_CLOEXEC));
    return Sha256Hex(ReadFd(fd.Value,64u*1024u*1024u));
}
std::string Binding(void* symbol)
{
    Dl_info info = {};
    Require(symbol && ::dladdr(symbol,&info) && info.dli_fname && info.dli_fbase,"symbol binding absent");
    uintptr_t address = reinterpret_cast<uintptr_t>(symbol), base = reinterpret_cast<uintptr_t>(info.dli_fbase);
    Require(address>=base,"symbol offset underflows");
    return O({{"address",Q(Address(address))},{"elf_offset",Q(Address(address-base))},{"path",Q(RealPath(info.dli_fname))}});
}
struct Api
{
    void* Dso = nullptr;
    std::string Path, Sha, Receipt, Copy;
    decltype(&wirehair_init_) Init = nullptr;
    decltype(&wirehair_v2_encoder_create_with_options) Create = nullptr;
    decltype(&wirehair_v2_encode) Encode = nullptr;
    decltype(&wirehair_v2_free) Free = nullptr;
    std::map<std::string,std::string> Symbols;
    ~Api() { if (Dso) ::dlclose(Dso); }
    Api() = default; Api(const Api&) = delete; Api& operator=(const Api&) = delete;
    template<class Function> void Load(Function& function, const char* name)
    {
        ::dlerror(); void* symbol = ::dlsym(Dso,name);
        Require(!::dlerror() && symbol,"missing public symbol");
        Dl_info info = {};
        Require(::dladdr(symbol,&info) && info.dli_fname && RealPath(info.dli_fname)==Path,"public symbol preempted");
        static_assert(sizeof(function)==sizeof(symbol),"POSIX function pointer size");
        std::memcpy(&function,&symbol,sizeof(function)); Symbols[name]=Binding(symbol);
    }
    void Open(const std::string& path, const std::string& expected = "")
    {
        Path = RealPath(path); Require(Path==path,"DSO path is not canonical"); Sha=FileHash(Path);
        Require(expected.empty() || Sha==expected,"DSO hash differs");
        Require(!::dlsym(RTLD_DEFAULT,"wirehair_v2_encode"),"wirehair is globally linked");
        Dso = ::dlopen(Path.c_str(),RTLD_NOW|RTLD_LOCAL); Require(Dso,"dlopen failed");
        Load(Init,"wirehair_init_"); Load(Create,"wirehair_v2_encoder_create_with_options");
        Load(Encode,"wirehair_v2_encode"); Load(Free,"wirehair_v2_free");
        void* copy = ::dlsym(Dso,"memcpy");
        Require(copy && copy==::dlsym(RTLD_DEFAULT,"memcpy"),"memcpy providers differ");
        Copy=Binding(copy);
        Require(Init(WIREHAIR_VERSION)==Wirehair_Success,"public init failed");
        Receipt=O({{"memcpy",Copy},{"path",Q(Path)},{"sha256",Q(Sha)},{"symbols",O(Symbols)}});
    }
    void Stable() const
    {
        Require(!::dlsym(RTLD_DEFAULT,"wirehair_v2_encode"),"wirehair became globally linked");
        for (const auto& symbol : Symbols)
            Require(Binding(::dlsym(Dso,symbol.first.c_str()))==symbol.second,"public binding changed");
        Require(Binding(::dlsym(Dso,"memcpy"))==Copy,"memcpy binding changed");
    }
};
struct Encoder
{
    Api* Owner = nullptr;
    WirehairV2Codec Handle = nullptr;
    std::array<uint8_t,32> Profile = {{0}};
    ~Encoder() { if (Handle) Owner->Free(Handle); }
    Encoder() = default; Encoder(const Encoder&) = delete; Encoder& operator=(const Encoder&) = delete;
    void Create(Api& api, const uint8_t* source, size_t bytes, uint32_t block)
    {
        Require(!Handle,"encoder already created"); Owner=&api;
        WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
        options.source_policy=WirehairV2EncoderSource_BorrowedImmutable;
        uint32_t written=0;
        Require(api.Create(source,bytes,block,&options,Profile.data(),32,&written,&Handle)==WirehairV2_Success &&
            Handle && written==32,"public create failed");
    }
};

struct Cell
{
    Shape Dimensions;
    std::vector<uint8_t> Bytes;
    std::string SourceHash, ProfileHash, Json;
};
std::string Oracle(Api& baseline, Api& candidate, Shape shape, uint32_t tail, bool fullRepair)
{
    const std::vector<uint8_t> source=Source(shape,tail);
    const std::vector<uint8_t> original=source;
    Encoder encoders[2]; encoders[0].Create(baseline,source.data(),source.size(),shape.B);
    encoders[1].Create(candidate,source.data(),source.size(),shape.B);
    Require(encoders[0].Profile==encoders[1].Profile,"oracle profiles differ");
    std::vector<uint8_t> systematic(source.size()), repairs(static_cast<size_t>(fullRepair?shape.K:1)*shape.B);
    std::vector<uint8_t> packet[2] = {std::vector<uint8_t>(shape.B+32),std::vector<uint8_t>(shape.B+32)};
    std::string high;
    const uint32_t repairCount=fullRepair?shape.K:1;
    for (uint32_t index=0; index<shape.K+repairCount+1; ++index) {
        const uint32_t id=index<shape.K+repairCount?index:UINT32_MAX;
        const uint32_t required=id==shape.K-1?tail:shape.B;
        for (unsigned arm=0; arm<2; ++arm) {
            std::fill(packet[arm].begin(),packet[arm].end(),0xa5);
            uint32_t written=0;
            Require(encoders[arm].Owner->Encode(encoders[arm].Handle,id,packet[arm].data()+16,
                required,&written)==WirehairV2_Success && written==required,"oracle encode failed");
            Require(std::all_of(packet[arm].begin(),packet[arm].begin()+16,[](uint8_t x){return x==0xa5;}) &&
                std::all_of(packet[arm].begin()+16+required,packet[arm].end(),[](uint8_t x){return x==0xa5;}),"oracle canary changed");
        }
        Require(packet[0]==packet[1],"oracle payload differs");
        if (id<shape.K) {
            Require(std::memcmp(packet[0].data()+16,source.data()+static_cast<size_t>(id)*shape.B,required)==0,
                "systematic oracle mismatch");
            std::memcpy(systematic.data()+static_cast<size_t>(id)*shape.B,packet[0].data()+16,required);
        } else if (id==UINT32_MAX) high=Sha256Hex(packet[0].data()+16,required);
        else std::memcpy(repairs.data()+static_cast<size_t>(id-shape.K)*shape.B,packet[0].data()+16,required);
    }
    Require(source==original,"oracle source changed");
    return O({{"high_id_sha256",Q(high)},{"profile_hex",Q(Hex(encoders[0].Profile.data(),32))},
        {"profile_sha256",Q(Sha256Hex(encoders[0].Profile.data(),32))},
        {"repair_sha256",Q(Sha256Hex(repairs.data(),repairs.size()))},
        {"source_sha256",Q(Sha256Hex(source.data(),source.size()))},
        {"systematic_sha256",Q(Sha256Hex(systematic.data(),systematic.size()))},{"tail_bytes",N(tail)}});
}
// Deliberately tiny extractor: only exact unique ASCII string fields produced by
// the canonical controller are accepted, and the complete CLAIM is hash-bound.
std::string Field(const std::string& text, const std::string& key)
{
    const std::string token=Q(key)+":\"";
    size_t start=text.find(token);
    Require(start!=std::string::npos && text.find(token,start+1)==std::string::npos,"missing/duplicate launch field");
    start+=token.size(); size_t end=text.find('"',start);
    Require(end!=std::string::npos,"unterminated launch field");
    std::string value=text.substr(start,end-start);
    Require(!value.empty() && value.find('\\')==std::string::npos &&
        value.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789/._-:")==std::string::npos,
        "unsupported launch field characters");
    return value;
}
Cell BuildCell(Api& baseline, Api& candidate, Shape shape)
{
    Cell cell; cell.Dimensions=shape; cell.Bytes=Source(shape,shape.B);
    cell.SourceHash=Sha256Hex(cell.Bytes.data(),cell.Bytes.size());
    const std::string full=Oracle(baseline,candidate,shape,shape.B,true);
    cell.ProfileHash=Field(full,"profile_sha256");
    cell.Json=O({{"K",N(shape.K)},{"block_bytes",N(shape.B)},
        {"high_id_sha256",Q(Field(full,"high_id_sha256"))},
        {"partial",A({Oracle(baseline,candidate,shape,1,false),Oracle(baseline,candidate,shape,shape.B-1,false)})},
        {"profile_hex",Q(Field(full,"profile_hex"))},{"profile_sha256",Q(cell.ProfileHash)},
        {"repair_sha256",Q(Field(full,"repair_sha256"))},{"source_sha256",Q(cell.SourceHash)}});
    return cell;
}

struct CounterBlock { uint32_t Values[8192]; };
struct Arena
{
    std::unique_ptr<uint8_t[]> Storage;
    uint8_t *Base=nullptr, *Input=nullptr, *Output=nullptr;
    CounterBlock* Counters=nullptr;
    size_t Span=0, Bytes=0, ArenaBytes=0;
    ~Arena() { if (Counters) Counters->~CounterBlock(); }
    void Initialize(const Cell& cell)
    {
        Bytes=cell.Bytes.size(); Span=((Bytes+8191)/4096)*4096;
        ArenaBytes=2*Span+1024+sizeof(CounterBlock)+64;
        Storage.reset(new uint8_t[ArenaBytes+4095]);
        Base=reinterpret_cast<uint8_t*>((reinterpret_cast<uintptr_t>(Storage.get())+4095)&~uintptr_t(4095));
        std::fill(Base,Base+ArenaBytes,0xa5);
        Input=Base+2048; Output=Base+Span+2336;
        Counters=new (Base+2*Span+1024) CounterBlock;
        Require(Input+Bytes<=Output && Output+Bytes<=reinterpret_cast<uint8_t*>(Counters),"arena overlaps");
        std::copy(cell.Bytes.begin(),cell.Bytes.end(),Input);
        std::fill(Counters->Values,Counters->Values+8192,0xa5a5a5a5u);
    }
    void Prepare(uint32_t k)
    { std::fill(Output,Output+Bytes,0xa5); std::fill(Counters->Values,Counters->Values+k,0xa5a5a5a5u); }
    bool Gap(const uint8_t* first, const uint8_t* last) const
    { return first<=last && std::all_of(first,last,[](uint8_t x){return x==0xa5;}); }
    void Verify(const Cell& cell, bool guards) const
    {
        Require(std::memcmp(Input,cell.Bytes.data(),Bytes)==0 && std::memcmp(Output,Input,Bytes)==0,"batch bytes changed");
        for (size_t i=0; i<8192; ++i)
            Require(Counters->Values[i]==(i<cell.Dimensions.K?cell.Dimensions.B:0xa5a5a5a5u),"counter changed");
        if (guards) Require(Gap(Base,Input) && Gap(Input+Bytes,Output) &&
            Gap(Output+Bytes,reinterpret_cast<uint8_t*>(Counters)) &&
            Gap(reinterpret_cast<uint8_t*>(Counters)+sizeof(CounterBlock),Base+ArenaBytes),"arena canary changed");
    }
    std::string Json(const Encoder (&encoders)[4]) const
    {
        std::vector<std::string> handles;
        for (const Encoder& encoder : encoders) handles.push_back(Q(Address(reinterpret_cast<uintptr_t>(encoder.Handle))));
        return O({{"arena",Q(Address(reinterpret_cast<uintptr_t>(Base)))},{"arena_bytes",N(ArenaBytes)},
            {"counters",Q(Address(reinterpret_cast<uintptr_t>(Counters)))},{"handles",A(handles)},
            {"output",Q(Address(reinterpret_cast<uintptr_t>(Output)))},{"source",Q(Address(reinterpret_cast<uintptr_t>(Input)))},
            {"span",N(Span)}});
    }
};
uint64_t Batch(const Cell& cell, Arena& arena, const Encoder& encoder, bool timed, uint32_t traversals)
{
    arena.Prepare(cell.Dimensions.K);
    const auto encode=encoder.Owner->Encode;
    const WirehairV2Codec handle=encoder.Handle;
    const uint32_t k=cell.Dimensions.K, block=cell.Dimensions.B;
    uint8_t* const output=arena.Output;
    uint32_t* const counters=arena.Counters->Values;
    Clock::time_point start;
    if (timed) start=Clock::now();
    for (uint32_t traversal=0; traversal<traversals; ++traversal)
        for (uint32_t id=0; id<k; ++id)
            if (encode(handle,id,output+static_cast<size_t>(id)*block,block,counters+id)!=WirehairV2_Success ||
                counters[id]!=block) throw std::runtime_error("batch encode failed");
    const int64_t nanoseconds=timed?
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-start).count():0;
    Require(!timed || nanoseconds>0,"nonpositive batch time");
    arena.Verify(cell,false);
    return static_cast<uint64_t>(nanoseconds);
}
std::string Identity()
{
    using namespace wirehair_wh2_bench;
    TargetIdentityReceiptV2 receipt; std::string diagnostic, bytes;
    Require(CapturePublicBorrowedTargetIdentity(120,receipt,diagnostic) &&
        SerializeTargetIdentityV2(receipt,bytes,diagnostic) && !bytes.empty() && bytes.size()<=4096 &&
        Sha256Hex(bytes)==receipt.CanonicalSha256,"target identity capture failed");
    const auto& raw=receipt.Raw; const auto& d=receipt.Derived;
    return O({{"after_cpu",N(static_cast<uint32_t>(receipt.After.Cpu))},{"before_cpu",N(static_cast<uint32_t>(receipt.Before.Cpu))},
        {"canonical_bytes",N(bytes.size())},{"canonical_hex",Q(Hex(bytes.data(),bytes.size()))},
        {"canonical_sha256",Q(receipt.CanonicalSha256)},
        {"capabilities",O({{"leaf1_ecx",N(raw.Leaf1.Regs.Ecx)},{"leaf1_edx",N(raw.Leaf1.Regs.Edx)},
            {"leaf6_eax",N(raw.Leaf6.Regs.Eax)},{"leaf6_ecx",N(raw.Leaf6.Regs.Ecx)},
            {"leaf80000001_ecx",N(raw.Leaf80000001.Regs.Ecx)},{"leaf80000001_edx",N(raw.Leaf80000001.Regs.Edx)},
            {"leaf80000008_ebx",N(raw.Leaf80000008.Regs.Ebx)},{"leaf80000021_eax",N(raw.Leaf80000021.Regs.Eax)},
            {"max_basic_leaf",N(raw.Leaf0.Regs.Eax)},{"max_extended_leaf",N(raw.Leaf80000000.Regs.Eax)}})},
        {"derived",O({{"ccd_id",N(d.CcdId)},{"complex_id",N(d.ComplexId)},{"core_id",N(d.CoreId)},
            {"family",N(d.Family)},{"full_apic_id",N(d.FullApicId)},{"initial_apic_id_8",N(d.InitialApicId8)},
            {"logical_processors_per_package",N(d.LogicalProcessorsPerPackage)},{"model",N(d.Model)},
            {"package_id",N(d.PackageId)},{"stepping",N(d.Stepping)},{"thread_id",N(d.ThreadId)},{"threads_per_core",N(d.ThreadsPerCore)}})},
        {"raw_capture_complete",receipt.RawCaptureComplete?"true":"false"},{"requested_cpu",N(120)},
        {"semantic_validation_passed",receipt.SemanticValidationPassed?"true":"false"},{"singleton_affinity_verified","true"}});
}
void Guard(const Clock::time_point& deadline)
{
    cpu_set_t affinity; CPU_ZERO(&affinity);
    Require(Clock::now()<deadline && ::sched_getcpu()==120 && ::sched_getaffinity(0,sizeof(affinity),&affinity)==0 &&
        CPU_COUNT(&affinity)==1 && CPU_ISSET(120,&affinity),"worker deadline or affinity differs");
}
void Emit(const std::string& line)
{ std::cout << line << '\n'; Require(bool(std::cout),"raw write failed"); }
std::string Terminal()
{
    return O({{"campaign",Q(Name)},{"encode_call_count",N(94371840)},{"measured_batch_count",N(9216)},
        {"measured_invocation_count",N(3345408)},{"panel_count",N(1152)},{"record_count",N(1154)},
        {"schema",Q(std::string(Prefix)+".terminal.v1")},{"status",Q("complete")},
        {"warmup_batch_count",N(2304)},{"warmup_invocation_count",N(836352)}});
}
void Run(Api& baseline, Api& candidate, const Clock::time_point& deadline)
{
    cpu_set_t affinity; CPU_ZERO(&affinity); CPU_SET(120,&affinity);
    Require(::sched_setaffinity(0,sizeof(affinity),&affinity)==0,"cannot pin worker");
    Guard(deadline); const std::string identity=Identity(), descriptionHash=Sha256Hex(Describe());
    std::vector<Cell> cells; std::vector<std::string> receipts;
    for (Shape shape : Shapes) { Guard(deadline); cells.push_back(BuildCell(baseline,candidate,shape)); receipts.push_back(cells.back().Json); }
    Require(Identity()==identity,"identity changed after oracles"); baseline.Stable(); candidate.Stable();
    Require(baseline.Copy==candidate.Copy,"DSO memcpy providers differ"); Guard(deadline);
    Emit(O({{"bindings",O({{"B",baseline.Receipt},{"C",candidate.Receipt}})},{"campaign",Q(Name)},
        {"cells",A(receipts)},{"description_sha256",Q(descriptionHash)},
        {"schema",Q(std::string(Prefix)+".config.v1")},{"target_identity",identity}}));
    uint32_t sequence=0;
    for (uint32_t rep=0; rep<12; ++rep) for (uint32_t cellIndex : CellOrder(rep)) {
        Guard(deadline); const Cell& cell=cells[cellIndex]; Arena arena; arena.Initialize(cell);
        Encoder encoders[4]; const uint32_t even[4]={0,2,3,1}, odd[4]={2,0,1,3};
        for (uint32_t position=0;position<4;++position) {
            const uint32_t index=rep%2?odd[position]:even[position];
            encoders[index].Create(index<2?baseline:candidate,arena.Input,arena.Bytes,cell.Dimensions.B);
            Require(Sha256Hex(encoders[index].Profile.data(),32)==cell.ProfileHash,"timed handle profile differs");
        }
        const std::string addresses=arena.Json(encoders);
        for (uint32_t comparisonPosition=0; comparisonPosition<4; ++comparisonPosition) {
            const uint32_t comparison=(rep+cellIndex+comparisonPosition)%4;
            for (uint32_t condition : Conditions[rep%4]) {
                Guard(deadline);
                for (uint32_t warm=0; warm<2; ++warm) {
                    const uint32_t physical=(warm^(condition%2))^(condition/2);
                    Batch(cell,arena,encoders[Pairs[comparison][physical]],false,8192/cell.Dimensions.K);
                }
                std::vector<std::string> slots;
                for (uint32_t logical : Slots(condition)) {
                    Guard(deadline); const uint32_t physical=logical^(condition/2);
                    const uint64_t elapsed=Batch(cell,arena,encoders[Pairs[comparison][physical]],true,8192/cell.Dimensions.K);
                    slots.push_back(O({{"elapsed_ns",N(elapsed)},{"logical_lane",Q(Side(logical))},{"physical_lane",N(physical)}}));
                }
                arena.Verify(cell,true);
                for (const Encoder& encoder : encoders)
                    Require(Sha256Hex(encoder.Profile.data(),32)==cell.ProfileHash,"panel profile changed");
                Guard(deadline);
                Emit(O({{"K",N(cell.Dimensions.K)},{"addresses",addresses},{"block_bytes",N(cell.Dimensions.B)},
                    {"campaign",Q(Name)},{"cell_index",N(cellIndex)},{"cell_key_sha256",Q(CellKey(rep,cellIndex))},
                    {"comparison_index",N(comparison)},{"condition",N(condition)},{"description_sha256",Q(descriptionHash)},
                    {"mapping",Q(Mapping(condition))},{"order",Q(Order(condition))},{"output_sha256",Q(cell.SourceHash)},
                    {"profile_sha256",Q(cell.ProfileHash)},{"replicate",N(rep)},{"schema",Q(std::string(Prefix)+".panel.v1")},
                    {"scope_invocations_per_batch",N(8192/cell.Dimensions.K)},{"sequence",N(sequence++)},{"slots",A(slots)},
                    {"source_immutable_verified","true"},{"source_sha256",Q(cell.SourceHash)}}));
            }
        }
    }
    Require(sequence==1152 && Identity()==identity,"terminal count or identity differs");
    baseline.Stable(); candidate.Stable(); Guard(deadline); Emit(Terminal());
}

void WriteAll(int fd, const std::string& bytes)
{
    size_t offset=0;
    while (offset<bytes.size()) {
        ssize_t written=::write(fd,bytes.data()+offset,bytes.size()-offset);
        if (written<0 && errno==EINTR) continue;
        Require(written>0,"marker write failed"); offset+=static_cast<size_t>(written);
    }
}
std::string Authorize(const std::string& directory, const std::string& claimHash,
    const std::string& workerHash, const std::string& descriptionHash)
{
    Require(directory==Directory && Hash(claimHash) && Hash(workerHash) && descriptionHash==Sha256Hex(Describe()),
        "launch authorization differs");
    Fd dir(::open(directory.c_str(),O_RDONLY|O_DIRECTORY|O_NOFOLLOW|O_CLOEXEC));
    struct stat statbuf;
    Require(::fstat(dir.Value,&statbuf)==0 && (statbuf.st_mode&07777)==0700,"directory mode differs");
    Fd claimFd(::openat(dir.Value,"CLAIM",O_RDONLY|O_NOFOLLOW|O_CLOEXEC));
    const std::string claim=ReadFd(claimFd.Value,1024*1024,0400);
    Require(Sha256Hex(claim)==claimHash && Field(claim,"campaign")==Name &&
        Field(claim,"worker_sha256")==workerHash && Field(claim,"description_sha256")==descriptionHash &&
        FileHash("/proc/self/exe")==workerHash,"claim/executable hash differs");
    Require(Hash(Field(claim,"baseline_sha256")) && Hash(Field(claim,"candidate_sha256")),"library hashes malformed");
    Fd marker(::openat(dir.Value,"WORKER_STARTED",O_WRONLY|O_CREAT|O_EXCL|O_NOFOLLOW|O_CLOEXEC,0400));
    Require(::fchmod(marker.Value,0400)==0,"marker mode failed");
    WriteAll(marker.Value,O({{"campaign",Q(Name)},{"claim_sha256",Q(claimHash)},{"description_sha256",Q(descriptionHash)},
        {"schema",Q(std::string(Prefix)+".worker-started.v1")},{"worker_sha256",Q(workerHash)}})+'\n');
    Require(::fsync(marker.Value)==0 && ::fsync(dir.Value)==0,"marker sync failed");
    return claim;
}
void PublicSelftest(const std::string& baselinePath, const std::string& candidatePath)
{
    Api baseline,candidate; baseline.Open(baselinePath); candidate.Open(candidatePath);
    Require(baseline.Copy==candidate.Copy,"selftest memcpy differs");
    for (Shape shape : Shapes) {
        Cell cell=BuildCell(baseline,candidate,shape); Arena arena; arena.Initialize(cell);
        Encoder encoders[4];
        for (uint32_t i=0;i<4;++i) encoders[i].Create(i<2?baseline:candidate,arena.Input,arena.Bytes,shape.B);
        for (uint32_t pair=0;pair<4;++pair) for (uint32_t cond=0;cond<4;++cond)
            for (uint32_t logical : Slots(cond)) Batch(cell,arena,encoders[Pairs[pair][logical^(cond/2)]],false,1);
        arena.Verify(cell,true);
    }
    baseline.Stable(); candidate.Stable();
}
#endif

void Selftest()
{
    Require(Hash(Sha256Hex(Describe())),"description hash malformed");
    uint64_t panels=0, measured=0, warm=0;
    for (uint32_t rep=0;rep<12;++rep) {
        const auto order=CellOrder(rep); std::set<uint32_t> unique(order.begin(),order.end());
        Require(unique.size()==6,"cell order not permutation");
        for (uint32_t cell : order) for (uint32_t comparison=0;comparison<4;++comparison)
            for (uint32_t condition : Conditions[rep%4]) {
                auto slots=Slots(condition); Require(std::count(slots.begin(),slots.end(),0u)==4,"slot balance differs");
                ++panels; measured+=8*(8192/Shapes[cell].K); warm+=2*(8192/Shapes[cell].K);
            }
    }
    Require(panels==1152 && measured==3345408 && warm==836352,"schedule counts differ");
    Require(Sha256Hex(Source(Shapes[0],64).data(),512)==
        "cba6894442c1997e328ad878e0e770b03a68b8eab70f2b03c630ae3c1c5eb5b4","source generator differs");
}
} // namespace offset_r0

int main(int argc, char** argv)
{
    using namespace offset_r0;
    try {
        if (argc==2 && std::string(argv[1])=="--describe") { std::cout<<Describe(); return std::cout?0:1; }
        if (argc>=2 && std::string(argv[1])=="--selftest") {
            Selftest();
            if (argc==6 && std::string(argv[2])=="--baseline" && std::string(argv[4])=="--candidate") {
#if defined(__linux__)
                PublicSelftest(argv[3],argv[5]);
#else
                throw std::runtime_error("dual DSO selftest requires Linux");
#endif
            } else Require(argc==2,"invalid selftest arguments");
            std::cout<<"offset-r0 selftest passed (untimed)\n"; return 0;
        }
#if defined(__linux__)
        const Clock::time_point deadline=Clock::now()+std::chrono::seconds(60);
        Require(argc==12 && std::string(argv[1])=="--run" && std::string(argv[2])=="--cpu" &&
            std::string(argv[3])=="120" && std::string(argv[4])=="--evidence-dir" &&
            std::string(argv[6])=="--claim-sha256" && std::string(argv[8])=="--worker-sha256" &&
            std::string(argv[10])=="--config-identity-sha256","invalid run arguments");
        const std::string claim=Authorize(argv[5],argv[7],argv[9],argv[11]);
        Api baseline,candidate;
        baseline.Open(Field(claim,"baseline_path"),Field(claim,"baseline_sha256"));
        candidate.Open(Field(claim,"candidate_path"),Field(claim,"candidate_sha256"));
        Require(baseline.Path!=candidate.Path,"two library paths required");
        Run(baseline,candidate,deadline); return 0;
#else
        throw std::runtime_error("authorized run requires Linux");
#endif
    } catch (const std::exception& error) { std::cerr<<"offset-r0: "<<error.what()<<'\n'; return 1; }
}
