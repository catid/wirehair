// .77: bounded allocation attribution only. Never a speed/recovery-rate gate.
#include "wirehair/wirehair.h"

#include <cstdarg>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <initializer_list>
#include <stdexcept>
#include <sys/resource.h>
#include <unistd.h>

namespace diagnostic {
static const unsigned kEventCap = 8192;
struct Event { unsigned phase, kind; size_t bytes; uintptr_t pointer, prior, caller; };
static Event events[kEventCap];
static unsigned active = 0, event_count = 0, phase_count = 0;
static bool overflow = false;
static char output[2u * 1024u * 1024u];
static size_t output_size = 0;

void Require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
void Record(unsigned kind, size_t bytes, const void* pointer, uintptr_t prior,
            const void* caller) noexcept {
    if (!active) return;
    if (event_count == kEventCap) { overflow = true; return; }
    events[event_count++] = Event{active, kind, bytes,
        reinterpret_cast<uintptr_t>(pointer), prior, reinterpret_cast<uintptr_t>(caller)};
}
void Append(const char* format, ...) {
    va_list args; va_start(args, format);
    const int n = vsnprintf(output + output_size, sizeof(output) - output_size, format, args);
    va_end(args);
    Require(n >= 0 && size_t(n) < sizeof(output) - output_size, "output cap");
    output_size += size_t(n);
}
void Write(int fd, const char* data, size_t bytes) {
    while (bytes) {
        const ssize_t n = write(fd, data, bytes);
        if (n < 0 && errno == EINTR) continue;
        Require(n > 0, "output write"); data += n; bytes -= size_t(n);
    }
}
void Marker(unsigned phase, bool begin) {
    char line[64];
    const int n = snprintf(line, sizeof(line), "WH2_ALLOC_PHASE %u %s\n", phase, begin ? "BEGIN" : "END");
    Require(n > 0 && size_t(n) < sizeof(line), "marker extent");
    Write(STDERR_FILENO, line, size_t(n));
}
void Touch() {
    // Logger storage, not codec storage. Prevent logger first-touch faults from
    // masquerading as codec allocations. No allocator setting is changed.
    volatile unsigned char* p = reinterpret_cast<volatile unsigned char*>(events);
    for (size_t i = 0; i < sizeof(events); ++i) p[i] = 0;
    p = reinterpret_cast<volatile unsigned char*>(output);
    for (size_t i = 0; i < sizeof(output); ++i) p[i] = 0;
}
template<class F> WirehairV2Result Phase(const char* operation, unsigned width,
    unsigned family, unsigned cycle, uint32_t packet, F function) {
    Require(active == 0 && !overflow && phase_count < 270, "phase entry/cap");
    const unsigned phase = ++phase_count, first = event_count;
    Marker(phase, true);
    rusage before = {}, after = {};
    Require(getrusage(RUSAGE_THREAD, &before) == 0, "rusage before");
    active = phase;
    WirehairV2Result result;
    try { result = function(); }
    catch (...) { active = 0; throw; }
    active = 0;
    Require(getrusage(RUSAGE_THREAD, &after) == 0, "rusage after");
    Marker(phase, false);
    Append("{\"type\":\"phase\",\"index\":%u,\"operation\":\"%s\",\"width\":%u,"
        "\"family\":%u,\"cycle\":%u,\"packet\":%u,\"result\":%d,\"first_event\":%u,"
        "\"end_event\":%u,\"counters_before\":[%ld,%ld,%ld,%ld],"
        "\"counters_after\":[%ld,%ld,%ld,%ld]}\n", phase, operation, width,
        family, cycle, packet, int(result), first, event_count,
        before.ru_minflt, before.ru_majflt, before.ru_nvcsw, before.ru_nivcsw,
        after.ru_minflt, after.ru_majflt, after.ru_nvcsw, after.ru_nivcsw);
    Require(!overflow, "event cap");
    return result;
}
void Hex(const uint8_t* data, size_t bytes) {
    static const char hex[] = "0123456789abcdef";
    Require(bytes <= (sizeof(output) - output_size - 1u) / 2u, "hex cap");
    for (size_t i = 0; i < bytes; ++i) {
        output[output_size++] = hex[data[i] >> 4];
        output[output_size++] = hex[data[i] & 15];
    }
}
uint32_t Packet(unsigned slot) { return slot < 12 ? slot : UINT32_MAX - 2u * (slot - 12); }
struct Handle {
    WirehairV2Codec value = nullptr;
    ~Handle() { wirehair_v2_free(value); }
    void Free() { wirehair_v2_free(value); value = nullptr; }
};
void Width(unsigned width) {
    uint8_t source[6u * 1280u], packets[18u * 1280u], recovered[6u * 1280u + 128u];
    const size_t bytes = 6u * width;
    for (size_t i = 0; i < sizeof(source); ++i) source[i] = uint8_t(37u * i + i / 11u);
    memset(packets, 0, sizeof(packets)); memset(recovered, 0xa5, sizeof(recovered));
    uint8_t profile[32] = {}; uint32_t profile_bytes = 0;
    WirehairV2EncoderOptions options = {};
    options.struct_bytes = sizeof(options); options.options_version = WIREHAIR_V2_ENCODER_OPTIONS_VERSION;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    Handle encoder;
    Require(wirehair_v2_encoder_create_with_options(source, bytes, width, &options,
        profile, sizeof(profile), &profile_bytes, &encoder.value) == WirehairV2_Success &&
        encoder.value && profile_bytes == sizeof(profile), "encoder fixture");
    for (unsigned slot = 0; slot < 18; ++slot) {
        uint32_t written = 0;
        Require(wirehair_v2_encode(encoder.value, Packet(slot), packets + size_t(slot) * width,
            width, &written) == WirehairV2_Success && written == width, "packet fixture");
    }
    encoder.Free();
    Append("{\"type\":\"fixture\",\"width\":%u,\"profile_hex\":\"", width); Hex(profile, sizeof(profile));
    Append("\",\"packets_hex\":\""); Hex(packets, 18u * width); Append("\"}\n");
    for (unsigned family = 0; family < 2; ++family) for (unsigned cycle = 0; cycle < 3; ++cycle) {
        Handle decoder;
        Require(Phase("create", width, family, cycle, 0, [&] {
            return wirehair_v2_decoder_create(profile, sizeof(profile), &decoder.value);
        }) == WirehairV2_Success && decoder.value, "decoder create");
        bool success = false;
        for (unsigned step = 0; step < 12; ++step) {
            const unsigned slot = step < 6 ? (family == 0 ? 6u : 12u) + step : step - 6u;
            const WirehairV2Result result = Phase("feed", width, family, cycle, Packet(slot), [&] {
                return wirehair_v2_decode(decoder.value, Packet(slot), packets + size_t(slot) * width, width);
            });
            Require(result == WirehairV2_Success || result == WirehairV2_NeedMore, "decoder feed");
            if (result == WirehairV2_Success) { Require(step >= 5, "early success"); success = true; break; }
        }
        Require(success, "decoder endpoint");
        memset(recovered, 0xa5, sizeof(recovered));
        uint64_t written = 0;
        Require(Phase("recover", width, family, cycle, 0, [&] {
            return wirehair_v2_recover(decoder.value, recovered + 64u, bytes, &written);
        }) == WirehairV2_Success && written == bytes &&
            memcmp(recovered + 64u, source, bytes) == 0, "recovered bytes");
        for (size_t i = 0; i < sizeof(recovered); ++i)
            if (i < 64u || i >= 64u + bytes) Require(recovered[i] == 0xa5, "recovery guard");
        Require(Phase("free", width, family, cycle, 0, [&] {
            decoder.Free(); return WirehairV2_Success;
        }) == WirehairV2_Success, "decoder free");
    }
}
void EmitEvents(bool success = true) {
    for (unsigned i = 0; i < event_count; ++i) {
        const Event& e = events[i];
        Append("{\"type\":\"allocation\",\"index\":%u,\"phase\":%u,\"kind\":%u,\"bytes\":%zu,"
            "\"pointer\":%llu,\"prior\":%llu,\"caller\":%llu}\n", i, e.phase, e.kind, e.bytes,
            (unsigned long long)e.pointer, (unsigned long long)e.prior, (unsigned long long)e.caller);
    }
    Append("{\"type\":\"footer\",\"outcome\":\"%s\",\"phases\":%u,\"events\":%u}\n",
        success ? "COMPLETE" : "INVALID", phase_count, event_count);
    Write(STDOUT_FILENO, output, output_size);
}
} // namespace diagnostic

// Link wrapping observes references in the linked objects; it does not claim
// visibility into allocator-internal/shared-library allocations. Forward the
// original operators (including their new-handler semantics) unchanged.
#define WRAP __attribute__((noinline, noipa))
#define RECORD(k,n,p,old) diagnostic::Record(k,n,p,old,__builtin_return_address(0))
extern "C" {
void* __real__Znwm(size_t); void* __real__Znam(size_t);
void __real__ZdlPv(void*); void __real__ZdaPv(void*);
void __real__ZdlPvm(void*,size_t); void __real__ZdaPvm(void*,size_t);
void* __real_malloc(size_t); void* __real_calloc(size_t,size_t);
void* __real_realloc(void*,size_t); void __real_free(void*);
WRAP void* __wrap__Znwm(size_t n) { void* p=__real__Znwm(n); RECORD(0,n,p,0); return p; }
WRAP void* __wrap__Znam(size_t n) { void* p=__real__Znam(n); RECORD(1,n,p,0); return p; }
WRAP void __wrap__ZdlPv(void* p) { RECORD(2,0,p,0); __real__ZdlPv(p); }
WRAP void __wrap__ZdaPv(void* p) { RECORD(3,0,p,0); __real__ZdaPv(p); }
WRAP void __wrap__ZdlPvm(void* p,size_t n) { RECORD(4,n,p,0); __real__ZdlPvm(p,n); }
WRAP void __wrap__ZdaPvm(void* p,size_t n) { RECORD(5,n,p,0); __real__ZdaPvm(p,n); }
WRAP void* __wrap_malloc(size_t n) { void* p=__real_malloc(n); RECORD(6,n,p,0); return p; }
WRAP void* __wrap_calloc(size_t n,size_t size) {
    void* p=__real_calloc(n,size);
    RECORD(7,size && n>SIZE_MAX/size ? SIZE_MAX : n*size,p,0); return p;
}
WRAP void* __wrap_realloc(void* p,size_t n) {
    const uintptr_t prior=reinterpret_cast<uintptr_t>(p);
    void* next=__real_realloc(p,n); RECORD(8,n,next,prior); return next;
}
WRAP void __wrap_free(void* p) { RECORD(9,0,p,0); __real_free(p); }
}

int main(int argc,char** argv) {
    using namespace diagnostic;
    try {
        if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
            Touch(); active = 1;
            void* p=__wrap__Znwm(17); Require(p!=nullptr,"neutral new"); __wrap__ZdlPv(p);
            p=__wrap__Znam(19); Require(p!=nullptr,"neutral array new"); __wrap__ZdaPv(p);
            p=__wrap__Znwm(23); Require(p!=nullptr,"neutral sized new"); __wrap__ZdlPvm(p,23);
            p=__wrap__Znam(29); Require(p!=nullptr,"neutral sized array new"); __wrap__ZdaPvm(p,29);
            p=__wrap_calloc(3,7); Require(p!=nullptr,"neutral calloc"); __wrap_free(p);
            p=__wrap_malloc(31); Require(p!=nullptr,"neutral malloc");
            memset(p,0x5a,31);
            void* next=__wrap_realloc(p,67); Require(next!=nullptr,"neutral realloc");
            for (unsigned i=0;i<31;++i) Require(static_cast<unsigned char*>(next)[i]==0x5a,"realloc bytes");
            __wrap_free(next);
            active=0;
            const unsigned kinds[]={0,2,1,3,0,4,1,5,7,9,6,8,9};
            const size_t sizes[]={17,0,19,0,23,23,29,29,21,0,31,67,0};
            Require(event_count==13 && !overflow,"neutral event count");
            for (unsigned i=0;i<13;++i) Require(events[i].kind==kinds[i] &&
                events[i].bytes==sizes[i] && events[i].phase==1 && events[i].caller!=0,"neutral wrapper ledger");
            for (unsigned i=0;i<10;i+=2) Require(events[i].pointer==events[i+1].pointer,"neutral free identity");
            Require(events[11].prior==events[10].pointer && events[12].pointer==events[11].pointer,"neutral realloc identity");
            Record(0,0,nullptr,0,nullptr); Require(event_count==13,"inactive log");
            event_count=kEventCap; active=1; Record(0,0,nullptr,0,nullptr); active=0;
            Require(overflow && event_count==kEventCap,"neutral bounded log");
            output_size=sizeof(output)-1;
            bool capped=false;
            try { Append("xx"); } catch (const std::exception&) { capped=true; }
            Require(capped && output_size==sizeof(output)-1,"neutral output cap");
            output_size=0;
            puts("PASS neutral wrapper forwarding, ledger and overflow; no codec work"); return 0;
        }
        Require(argc==3 && strcmp(argv[1],"--worker")==0 && strlen(argv[2])==64,"explicit worker claim required");
        for (unsigned i=0;i<64;++i) Require((argv[2][i]>='0' && argv[2][i]<='9') ||
            (argv[2][i]>='a' && argv[2][i]<='f'),"claim hex");
        const rlimit cpu={10,10},memory={256u*1024u*1024u,256u*1024u*1024u},core={0,0};
        Require(setrlimit(RLIMIT_CPU,&cpu)==0 && setrlimit(RLIMIT_AS,&memory)==0 &&
            setrlimit(RLIMIT_CORE,&core)==0,"worker limits");
        Touch(); Require(wirehair_init()==Wirehair_Success,"GF initialization");
        Append("{\"type\":\"header\",\"protocol\":\"wirehair.wh2.decoder-allocation-r0\","
            "\"claim_sha256\":\"%s\",\"pid\":%ld,\"speed_claimed\":false}\n",argv[2],long(getpid()));
        for (unsigned width : {2u,64u,1280u}) Width(width);
        EmitEvents(); return 0;
    } catch (const std::exception& error) {
        diagnostic::active=0;
        // Preserve the diagnostic prefix and recorded events on ordinary API
        // failures too. A killed process may retain only its syscall markers.
        if (diagnostic::output_size) {
            try { diagnostic::EmitEvents(false); } catch (...) {}
        }
        fprintf(stderr,"INVALID: %s\n",error.what()); return 1;
    }
}
