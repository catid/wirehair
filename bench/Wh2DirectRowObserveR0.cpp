// Additive .66 observer; the immutable .65 entry point is never dispatched.
#define main Wh2MulRowsR0DormantMain
#include "Wh2ThueMorseMulRowsR0Bench.cpp"
#undef main

#include <limits>
#include <time.h>

// The old included worker has a dormant reference to this symbol. Do not link
// the actual candidate object: any attempt to enter that path is fatal.
namespace wh2_thue_mulrows_r0 {
wh2_thue_native::Status Row(wh2_thue_native::Lookup, std::uint32_t, std::uint8_t*)
{
    throw std::runtime_error("candidate path forbidden in direct-row observation");
}
}

namespace observe {
static const char kProtocol[] = "wirehair.wh2.direct-row-observe-r0";
static const uint64_t kMissing = UINT64_MAX;
static const uint64_t kWallLimit = UINT64_C(10000000000);
static const uint64_t kRowLimit = UINT64_C(100000000);
using Counters = std::array<uint64_t,4>;

struct Record {
    unsigned index=0,rep=0,order=0,position=0;
    uint64_t m0=kMissing,c0=kMissing,m1=kMissing,m2=kMissing,c1=kMissing,m3=kMissing;
    Counters ru0={{kMissing,kMissing,kMissing,kMissing}};
    Counters ru1={{kMissing,kMissing,kMissing,kMissing}};
    CallResult result={0,false};
};

struct State {
    std::array<Record,240> records;
    size_t started=0;
    uint64_t callbacks=0,row_calls=0,checked_rows=0,checksum=kChecksumSeed,sum_row_wall=0;
    State() {
        unsigned index=0;
        for (unsigned rep=0; rep<12; ++rep)
            for (unsigned step=0; step<2; ++step)
                for (unsigned position=0; position<10; ++position) {
                    Record& record=records[index];
                    record.index=index++; record.rep=rep;
                    record.order=(rep+step)%2; record.position=position;
                }
    }
};

// Retain physical page touches even under the frozen optimizing compiler.
void Touch(void* storage, size_t bytes)
{
    volatile uint8_t* pointer=static_cast<volatile uint8_t*>(storage);
    for (size_t i=0; i<bytes; ++i) {
        const uint8_t value=pointer[i];
        pointer[i]=value;
    }
}

bool Nanoseconds(const struct timespec& value, uint64_t& output)
{
    if (value.tv_sec < 0 || value.tv_nsec < 0 || value.tv_nsec >= 1000000000L ||
        static_cast<uint64_t>(value.tv_sec) > (kMissing-1)/UINT64_C(1000000000)) return false;
    const uint64_t base=static_cast<uint64_t>(value.tv_sec)*UINT64_C(1000000000);
    if (static_cast<uint64_t>(value.tv_nsec) >= kMissing-base) return false;
    output=base+static_cast<uint64_t>(value.tv_nsec);
    return true;
}

struct SystemReader {
    struct timespec time_storage = {};
    struct rusage usage_storage = {};
    void Prepare() {
        Touch(&time_storage,sizeof(time_storage));
        Touch(&usage_storage,sizeof(usage_storage));
    }
    bool Clock(clockid_t id, uint64_t& value) {
        return clock_gettime(id,&time_storage) == 0 && Nanoseconds(time_storage,value);
    }
    bool Mono(uint64_t& value) { return Clock(CLOCK_MONOTONIC,value); }
    bool Cpu(uint64_t& value) { return Clock(CLOCK_THREAD_CPUTIME_ID,value); }
    bool Usage(Counters& value) {
        if (getrusage(RUSAGE_THREAD,&usage_storage) != 0 || usage_storage.ru_minflt < 0 ||
            usage_storage.ru_majflt < 0 || usage_storage.ru_nvcsw < 0 || usage_storage.ru_nivcsw < 0)
            return false;
        value={{static_cast<uint64_t>(usage_storage.ru_minflt),
                static_cast<uint64_t>(usage_storage.ru_majflt),
                static_cast<uint64_t>(usage_storage.ru_nvcsw),
                static_cast<uint64_t>(usage_storage.ru_nivcsw)}};
        return true;
    }
    uint64_t Resolution(clockid_t id) {
        uint64_t value=kMissing;
        Require(clock_getres(id,&time_storage) == 0 && Nanoseconds(time_storage,value) &&
                value > 0,"clock resolution");
        return value;
    }
    uint64_t Now() {
        uint64_t value=kMissing;
        Require(Mono(value),"worker monotonic clock");
        return value;
    }
};

void CheckBudget(uint64_t started, uint64_t now, uint64_t row_wall)
{
    Require(now >= started && now-started < kWallLimit,"observer wall deadline");
    Require(row_wall <= kRowLimit,"cumulative row wall cap");
}

void FinalClock(bool captured, uint64_t now, uint64_t started, uint64_t latest,
                uint64_t row_wall, uint64_t& elapsed, std::string& outcome, std::string& failure)
{
    const char* error=nullptr;
    if (!captured || now < started || now < latest ||
        (outcome == "DIAGNOSTIC_COMPLETE" && now-started < elapsed)) error="final observer clock";
    else {
        elapsed=now-started;
        if (elapsed >= kWallLimit || row_wall > kRowLimit) error="post-observation deadline";
    }
    if (error) {
        outcome="INVALID";
        if (failure.empty()) failure=error;
    }
}

// No allocation, CPU checks, row verification, or checksum in this bracket.
// Failure returns a literal diagnostic, leaving every available field intact.
template<class Reader>
const char* Capture(Reader& reader, RowFunction function, N::Lookup lookup,
                    const Ids& ids, Output& output, Record& record)
{
    if (!reader.Mono(record.m0)) return "m0 capture";
    if (!reader.Usage(record.ru0)) return "ru0 capture";
    if (!reader.Cpu(record.c0)) return "c0 capture";
    if (!reader.Mono(record.m1)) return "m1 capture";
    record.result=RunRows(function,lookup,ids,output);
    if (!reader.Mono(record.m2)) return "m2 capture";
    if (!reader.Cpu(record.c1)) return "c1 capture";
    if (!reader.Usage(record.ru1)) return "ru1 capture";
    if (!reader.Mono(record.m3)) return "m3 capture";
    return nullptr;
}

void Validate(const Record& record, const Record* previous)
{
    Require(record.m0 != kMissing && record.c0 != kMissing && record.m1 != kMissing &&
            record.m2 != kMissing && record.c1 != kMissing && record.m3 != kMissing,
            "incomplete clock observation");
    Require(record.m0 <= record.m1 && record.m1 < record.m2 && record.m2 <= record.m3 &&
            record.c0 <= record.c1,"clock ordering");
    // CPU and monotonic timestamps belong to different clock domains. In
    // particular, thread span may legitimately exceed the inner row wall.
    for (unsigned k=0; k<4; ++k)
        Require(record.ru0[k] != kMissing && record.ru1[k] != kMissing &&
                record.ru0[k] <= record.ru1[k],"counter ordering");
    if (previous) {
        Require(previous->m3 <= record.m0 && previous->c1 <= record.c0,
                "cross-record clock ordering");
        for (unsigned k=0; k<4; ++k)
            Require(previous->ru1[k] <= record.ru0[k],"cross-record counter ordering");
    }
}

// This also runs after a failed capture, so a completed Row invocation and its
// available output/timestamps are not silently erased from the invalid prefix.
void Account(State& state, const Record& record, const Output& output, const Rows& expected)
{
    if (record.m1 != kMissing && record.m2 != kMissing && record.m2 >= record.m1) {
        const uint64_t elapsed=record.m2-record.m1;
        Require(elapsed <= UINT64_MAX-state.sum_row_wall,"row duration overflow");
        state.sum_row_wall+=elapsed;
    }
    if (record.result.calls) {
        ++state.callbacks;
        state.row_calls+=record.result.calls;
        Require(record.result.calls == 4096 && record.result.success,"Row status/count mismatch");
        Require(output.Matches(expected),"observation row/guard mismatch");
        state.checked_rows+=64;
        state.checksum=Consume(state.checksum,output);
    }
}

void Number(std::ostream& out, uint64_t value)
{
    if (value == kMissing) out << "null";
    else out << value;
}

void UsageJson(std::ostream& out, const Counters& value)
{
    out << '[';
    for (unsigned k=0; k<4; ++k) { if (k) out << ','; Number(out,value[k]); }
    out << ']';
}

std::string RecordsJson(const State& state)
{
    std::ostringstream out;
    out << '[';
    for (size_t i=0; i<state.started; ++i) {
        if (i) out << ',';
        const Record& r=state.records[i];
        out << '[' << r.index << ',' << r.rep << ',' << r.order << ',' << r.position;
        for (uint64_t value : {r.m0,r.c0,r.m1,r.m2,r.c1,r.m3}) {
            out << ','; Number(out,value);
        }
        out << ','; UsageJson(out,r.ru0); out << ','; UsageJson(out,r.ru1); out << ']';
    }
    out << ']';
    return out.str();
}

// Only the new observer/selftest can be selected. The included dormant main
// and Worker functions are never used by either path below.
enum class Selection { Invalid, Diagnostic, Neutral };
Selection Select(int argc, const char* const* argv)
{
    if (argc != 2 || !argv || !argv[1]) return Selection::Invalid;
    if (std::strcmp(argv[1],"--worker") == 0) return Selection::Diagnostic;
    if (std::strcmp(argv[1],"--selftest") == 0) return Selection::Neutral;
    return Selection::Invalid;
}

// Fixed-size fake instrumentation: no OS calls, GF initialization, or lookup.
struct FakeReader {
    unsigned reads=0,fail_at=0,event_count=0;
    uint64_t tick=0,usage_count=0;
    std::array<char,9> events = {{0}};
    bool Step(char event) {
        Require(event_count < events.size(),"mock event overflow");
        events[event_count++]=event; ++reads; ++tick;
        return reads != fail_at;
    }
    bool Mono(uint64_t& value) {
        if (!Step('M')) return false;
        value=tick*100; return true;
    }
    bool Cpu(uint64_t& value) {
        if (!Step('C')) return false;
        value=tick*100; return true;
    }
    bool Usage(Counters& value) {
        if (!Step('U')) return false;
        value={{usage_count,usage_count,usage_count,usage_count}};
        ++usage_count; return true;
    }
    void Row() {
        Require(event_count < events.size(),"mock row event overflow");
        events[event_count++]='R'; tick+=10;
    }
};
FakeReader* fake_reader=nullptr;
uint64_t fake_calls=0;
N::Status FakeRow(N::Lookup, uint32_t id, uint8_t* output)
{
    if (fake_calls++ == 0) fake_reader->Row();
    for (unsigned k=0; k<6; ++k) output[k]=static_cast<uint8_t>((id>>((k%4)*8))^k);
    return N::Status::Success;
}

template<class Function>
void Rejected(Function function, const char* expected)
{
    bool rejected=false;
    try { function(); } catch (const std::runtime_error& error) {
        rejected=std::strcmp(error.what(),expected) == 0;
    }
    Require(rejected,"expected neutral rejection");
}

int Neutral()
{
    const char* worker[]={"test","--worker"};
    const char* selftest[]={"test","--selftest"};
    const char* wrong[]={"test","--old-worker"};
    const char* missing[]={"test",nullptr};
    Require(Select(2,worker) == Selection::Diagnostic && Select(2,selftest) == Selection::Neutral &&
            Select(2,wrong) == Selection::Invalid && Select(1,selftest) == Selection::Invalid &&
            Select(3,selftest) == Selection::Invalid && Select(2,missing) == Selection::Invalid &&
            Select(2,nullptr) == Selection::Invalid,"new main selection");
    const N::Lookup neutral={nullptr,0};
    Rejected([&] { C::Row(neutral,0,nullptr); },
             "candidate path forbidden in direct-row observation");
    uint64_t ns=kMissing;
    struct timespec time_value={1,2};
    Require(Nanoseconds(time_value,ns) && ns == 1000000002,"timespec conversion");
    time_value.tv_nsec=1000000000L; Require(!Nanoseconds(time_value,ns),"bad timespec nanoseconds");
    time_value.tv_nsec=0; time_value.tv_sec=-1; Require(!Nanoseconds(time_value,ns),"negative timespec");
    time_value.tv_sec=std::numeric_limits<time_t>::max();
    Require(!Nanoseconds(time_value,ns),"timespec overflow");
    CheckBudget(1,kWallLimit,kRowLimit);
    Rejected([] { CheckBudget(1,0,0); },"observer wall deadline");
    Rejected([] { CheckBudget(1,kWallLimit+1,0); },"observer wall deadline");
    Rejected([] { CheckBudget(1,2,kRowLimit+1); },"cumulative row wall cap");
    for (unsigned fault=0; fault<4; ++fault) {
        std::string outcome="DIAGNOSTIC_COMPLETE",failure;
        uint64_t elapsed=123;
        FinalClock(fault != 0,fault == 1 ? 0 : fault == 2 ? 2 : kWallLimit+1,
                   1,3,0,elapsed,outcome,failure);
        Require(outcome == "INVALID" && !failure.empty(),"final clock failure invalidates");
    }
    {
        std::string outcome="INVALID",failure="original failure";
        uint64_t elapsed=0;
        FinalClock(false,0,1,2,0,elapsed,outcome,failure);
        Require(failure == "original failure","first failure preserved");
        outcome="DIAGNOSTIC_COMPLETE"; failure.clear();
        FinalClock(true,5,1,4,0,elapsed,outcome,failure);
        Require(outcome == "DIAGNOSTIC_COMPLETE" && failure.empty() && elapsed == 4,
                "valid final clock");
        elapsed=5;
        FinalClock(true,5,1,4,0,elapsed,outcome,failure);
        Require(outcome == "INVALID" && failure == "final observer clock",
                "regression below prior successful elapsed");
    }
    const Ids ids=Family(0);
    Rows expected;
    for (unsigned j=0; j<64; ++j) {
        Require(ids[j] == 13*j,"neutral low roster");
        for (unsigned k=0; k<6; ++k) expected[j][k]=static_cast<uint8_t>((ids[j]>>((k%4)*8))^k);
    }
    State state;
    Output output;
    output.Reset(); Touch(&output,sizeof(output)); Touch(&state.records,sizeof(state.records));
    for (size_t i=0; i<state.records.size(); ++i) {
        Record& r=state.records[i];
        Require(r.index == i && r.rep == i/20 && r.order == (i/20+(i/10)%2)%2 &&
                r.position == i%10,"neutral complete observation roster");
        FakeReader reader;
        reader.tick=static_cast<uint64_t>(i)*100;
        reader.usage_count=static_cast<uint64_t>(i)*2;
        fake_reader=&reader; fake_calls=0; output.Reset();
        state.started=i+1;
        Require(Capture(reader,FakeRow,neutral,ids,output,r) == nullptr,"neutral capture");
        Require(reader.events == std::array<char,9>{{'M','U','C','M','R','M','C','U','M'}} &&
                fake_calls == 4096,"exact observation ordering");
        Validate(r,i ? &state.records[i-1] : nullptr);
        Require(r.m2-r.m1 == 1100 && r.c1-r.c0 == 1300 && r.m3-r.m0 == 1700,
                "enclosing clock spans");
        Account(state,r,output,expected);
    }
    Require(state.started == 240 && state.callbacks == 240 && state.row_calls == 983040 &&
            state.checked_rows == 15360 && state.sum_row_wall == 264000,"neutral exact counts");
    uint64_t expected_checksum=kChecksumSeed;
    for (unsigned i=0; i<240; ++i)
        for (const auto& row : expected)
            for (uint8_t byte : row) expected_checksum=(expected_checksum*kChecksumFactor)^byte;
    Require(state.checksum == expected_checksum,"neutral checksum");
    Record bad=state.records[0];
    bad.m2=bad.m1;
    Rejected([&] { Validate(bad,nullptr); },"clock ordering");
    bad=state.records[0]; bad.ru1[0]=0; bad.ru0[0]=1;
    Rejected([&] { Validate(bad,nullptr); },"counter ordering");
    bad=state.records[1]; bad.c0=state.records[0].c1-1;
    Rejected([&] { Validate(bad,&state.records[0]); },"cross-record clock ordering");
    bad=state.records[1]; bad.ru0[2]=0;
    Rejected([&] { Validate(bad,&state.records[0]); },"cross-record counter ordering");
    for (unsigned fail=1; fail<=8; ++fail) {
        State partial;
        FakeReader reader;
        reader.fail_at=fail; fake_reader=&reader; fake_calls=0; output.Reset();
        partial.started=1;
        const char* error=Capture(reader,FakeRow,neutral,ids,output,partial.records[0]);
        Require(error != nullptr && partial.records[0].m3 == kMissing,"capture error retained");
        Account(partial,partial.records[0],output,expected);
        Require(partial.callbacks == (fail>=5 ? 1u : 0u) &&
                partial.row_calls == (fail>=5 ? 4096u : 0u),"invalid prefix call accounting");
        const std::string json=RecordsJson(partial);
        Require(json.find("[[0,0,0,0,") == 0 && json.find("null") != std::string::npos,
                "partial observation publication");
    }
    State corrupt;
    corrupt.started=1; corrupt.records[0]=state.records[0];
    output.slots[63].after[15]^=1;
    Rejected([&] { Account(corrupt,corrupt.records[0],output,expected); },"observation row/guard mismatch");
    Require(corrupt.callbacks == 1 && corrupt.row_calls == 4096 && corrupt.checked_rows == 0 &&
            RecordsJson(corrupt).find("null") == std::string::npos,"guard failure complete observation retained");
    fake_reader=nullptr;
    std::cout << "DIRECT ROW OBSERVE R0 neutral selftest PASS\n";
    std::cout.flush();
    return std::cout.good() ? 0 : 1;
}

int Diagnostic()
{
    SystemReader reader;
    State state;
    Output output;
    std::string outcome="INVALID",failure,identity_before="null",identity_after="null";
    std::string runtime_before="null",runtime_after="null",expected_bytes,expected_hash;
    uint64_t started=kMissing,elapsed=0,monotonic_resolution=0,thread_resolution=0,preflight=0;
    const N::Lookup lookup={D::kLookup,sizeof(D::kLookup)};
    try {
        started=reader.Now();
        const struct rlimit memory={512u*1024u*1024u,512u*1024u*1024u},cpu={5,5};
        Require(setrlimit(RLIMIT_AS,&memory) == 0 && setrlimit(RLIMIT_CPU,&cpu) == 0,
                "observer resource limits");
        alarm(10);
        Pin(); identity_before=Identity();
        Require(gf256_init() == 0,"shared GF initialization");
        runtime_before=Runtime();
        Require(lookup.bytes == 39936 && Sha256Hex(lookup.data,lookup.bytes) == D::kLookupSha256,
                "immutable lookup hash");
        monotonic_resolution=reader.Resolution(CLOCK_MONOTONIC);
        thread_resolution=reader.Resolution(CLOCK_THREAD_CPUTIME_ID);
        const Ids ids=Family(0);
        Rows expected;
        for (unsigned j=0; j<64; ++j) {
            const auto oracle=Oracle(lookup,ids[j]);
            Slot baseline;
            std::memset(&baseline,0xa5,sizeof(baseline));
            ++preflight;
            Require(N::Row(lookup,ids[j],baseline.row) == N::Status::Success,"low preflight status");
            Require(std::memcmp(baseline.row,oracle.data(),6) == 0,"low preflight oracle");
            for (unsigned k=0; k<16; ++k)
                Require(baseline.before[k] == 0xa5 && baseline.after[k] == 0xa5,"low preflight guard");
            expected[j]=oracle;
            expected_bytes.append(reinterpret_cast<const char*>(oracle.data()),6);
        }
        expected_hash=Sha256Hex(expected_bytes);
        Require(preflight == 64 && expected_bytes.size() == 384,"low preflight counts");
        output.Reset(); reader.Prepare();
        Touch(&state,sizeof(state)); Touch(&output,sizeof(output));
        // Read-only input pages and stack temporaries are also touched before
        // the first recorded callback; no additional Row warmup is introduced.
        volatile uint8_t input_touch=0;
        for (unsigned j=0; j<64; ++j) input_touch^=expected[j][0];
        (void)input_touch;
        for (size_t i=0; i<state.records.size(); ++i) {
            CheckBudget(started,reader.Now(),state.sum_row_wall);
            Require(OnCpu(),"observer CPU drift");
            output.Reset();
            Record& record=state.records[i];
            state.started=i+1;
            const char* capture_failure=Capture(reader,N::Row,lookup,ids,output,record);
            Account(state,record,output,expected);
            Require(capture_failure == nullptr,capture_failure ? capture_failure : "capture failure");
            Validate(record,i ? &state.records[i-1] : nullptr);
            Require(OnCpu(),"observer CPU drift");
            CheckBudget(started,reader.Now(),state.sum_row_wall);
        }
        Require(state.started == 240 && state.callbacks == 240 && state.row_calls == 983040 &&
                state.checked_rows == 15360,"complete observer accounting");
        identity_after=Identity(); runtime_after=Runtime();
        Require(identity_after == identity_before && runtime_after == runtime_before,
                "post-observation identity/runtime drift");
        Require(Sha256Hex(lookup.data,lookup.bytes) == D::kLookupSha256,"post-observation lookup drift");
        elapsed=reader.Now()-started;
        CheckBudget(0,elapsed,state.sum_row_wall);
        outcome="DIAGNOSTIC_COMPLETE";
    } catch (const std::exception& error) { failure=error.what(); }
      catch (...) { failure="unknown observer failure"; }
    if (started != kMissing) {
        uint64_t now=kMissing,latest=started;
        for (size_t i=0; i<state.started; ++i) {
            const Record& r=state.records[i];
            for (uint64_t value : {r.m0,r.m1,r.m2,r.m3})
                if (value != kMissing && value > latest) latest=value;
        }
        const bool captured=reader.Mono(now);
        FinalClock(captured,now,started,latest,state.sum_row_wall,elapsed,outcome,failure);
    }
    std::ostringstream out;
    out << "{\"protocol\":" << Quote(kProtocol) << ",\"schema\":" << Quote(std::string(kProtocol)+".raw.v1")
        << ",\"outcome\":" << Quote(outcome) << ",\"target_cpu\":50"
        << ",\"lookup_sha256\":" << Quote(D::kLookupSha256)
        << ",\"expected_rows_hex\":" << Quote(Hex(expected_bytes))
        << ",\"expected_rows_sha256\":" << Quote(expected_hash)
        << ",\"output_address\":" << Address(output.slots[0].row) << ",\"output_stride\":" << sizeof(Slot)
        << ",\"baseline_function_address\":" << static_cast<uint64_t>(reinterpret_cast<uintptr_t>(&N::Row))
        << ",\"identity_before\":" << identity_before << ",\"identity_after\":" << identity_after
        << ",\"runtime_before\":" << runtime_before << ",\"runtime_after\":" << runtime_after
        << ",\"monotonic_resolution_ns\":" << monotonic_resolution << ",\"thread_resolution_ns\":" << thread_resolution
        << ",\"preflight_row_calls\":" << preflight << ",\"callbacks\":" << state.callbacks
        << ",\"row_calls\":" << state.row_calls << ",\"checked_rows\":" << state.checked_rows
        << ",\"checksum\":" << state.checksum << ",\"sum_row_wall_ns\":" << state.sum_row_wall
        << ",\"elapsed_ns\":" << elapsed << ",\"records\":" << RecordsJson(state)
        << ",\"failure\":" << (failure.empty() ? "null" : Quote(failure)) << "}\n";
    const std::string raw=out.str();
    Require(raw.size() <= 1024u*1024u,"observer raw cap");
    std::cout << raw; std::cout.flush();
    return std::cout.good() && outcome == "DIAGNOSTIC_COMPLETE" ? 0 : 1;
}

} // namespace observe

int main(int argc, char** argv)
{
    const observe::Selection selection=observe::Select(argc,argv);
    if (selection == observe::Selection::Invalid) return 2;
    try {
        return selection == observe::Selection::Neutral ? observe::Neutral() : observe::Diagnostic();
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
      catch (...) { std::cerr << "observer publication failure\n"; return 1; }
}
