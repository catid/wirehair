// Mechanism diagnostic only: no timing or allocation-footprint inference.
#include "wirehair/wirehair.h"
#include <cerrno>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <sys/resource.h>
#include <unistd.h>

namespace {
char output[256u * 1024u];
size_t used = 0;
unsigned phases = 0;

void Check(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}
void Append(const char* format, ...) {
    va_list args; va_start(args, format);
    const int size = vsnprintf(output + used, sizeof(output) - used, format, args);
    va_end(args);
    Check(size >= 0 && size_t(size) < sizeof(output) - used, "output cap");
    used += size_t(size);
}
void Write(int fd, const char* data, size_t bytes) {
    while (bytes) {
        const ssize_t size = write(fd, data, bytes);
        if (size < 0 && errno == EINTR) continue;
        Check(size > 0, "write"); data += size; bytes -= size_t(size);
    }
}
void Hex(const uint8_t* data, size_t bytes) {
    const char* hex = "0123456789abcdef";
    Check(bytes <= (sizeof(output) - used - 1u) / 2u, "hex cap");
    for (size_t i = 0; i < bytes; ++i) {
        output[used++] = hex[data[i] >> 4]; output[used++] = hex[data[i] & 15];
    }
}
void Marker(unsigned index, const char* side) {
    char text[64];
    const int bytes = snprintf(text, sizeof(text), "WH2_RELEASE_PHASE %u %s\n", index, side);
    Check(bytes > 0 && size_t(bytes) < sizeof(text), "marker cap");
    Write(STDERR_FILENO, text, size_t(bytes));
}
template<class F> WirehairV2Result Phase(const char* operation, unsigned width,
    unsigned family, unsigned cycle, uint32_t packet, F function) {
    Check(phases < 270, "phase cap");
    const unsigned index = ++phases;
    Marker(index, "BEGIN");
    rusage before = {}, after = {};
    Check(getrusage(RUSAGE_THREAD, &before) == 0, "before counters");
    const WirehairV2Result result = function();
    Check(getrusage(RUSAGE_THREAD, &after) == 0, "after counters");
    Marker(index, "END");
    Append("{\"type\":\"phase\",\"index\":%u,\"operation\":\"%s\",\"width\":%u,"
        "\"family\":%u,\"cycle\":%u,\"packet\":%u,\"result\":%d,"
        "\"before\":[%ld,%ld,%ld,%ld],\"after\":[%ld,%ld,%ld,%ld]}\n",
        index, operation, width, family, cycle, packet, int(result),
        before.ru_minflt, before.ru_majflt, before.ru_nvcsw, before.ru_nivcsw,
        after.ru_minflt, after.ru_majflt, after.ru_nvcsw, after.ru_nivcsw);
    return result;
}
uint32_t Packet(unsigned slot) {
    return slot < 12u ? slot : UINT32_MAX - 2u * (slot - 12u);
}
struct Handle {
    WirehairV2Codec value = nullptr;
    ~Handle() { wirehair_v2_free(value); }
    void Free() { wirehair_v2_free(value); value = nullptr; }
};
void Width(unsigned width) {
    uint8_t source[7680], packets[23040], recovered[7808];
    const size_t bytes = 6u * width;
    for (size_t i = 0; i < sizeof(source); ++i) source[i] = uint8_t(37u * i + i / 11u);
    memset(packets, 0, sizeof(packets)); memset(recovered, 0xa5, sizeof(recovered));
    uint8_t profile[32] = {}; uint32_t profile_bytes = 0;
    WirehairV2EncoderOptions options = {};
    options.struct_bytes = sizeof(options); options.options_version = WIREHAIR_V2_ENCODER_OPTIONS_VERSION;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    Handle encoder;
    Check(wirehair_v2_encoder_create_with_options(source, bytes, width, &options,
        profile, sizeof(profile), &profile_bytes, &encoder.value) == WirehairV2_Success &&
        encoder.value && profile_bytes == sizeof(profile), "fixture encoder");
    for (unsigned slot = 0; slot < 18; ++slot) {
        uint32_t written = 0;
        Check(wirehair_v2_encode(encoder.value, Packet(slot), packets + size_t(slot)*width,
            width, &written) == WirehairV2_Success && written == width, "fixture packet");
    }
    encoder.Free();
    Append("{\"type\":\"fixture\",\"width\":%u,\"profile_hex\":\"", width); Hex(profile, sizeof(profile));
    Append("\",\"packets_hex\":\""); Hex(packets, 18u*width); Append("\"}\n");
    for (unsigned family = 0; family < 2; ++family) for (unsigned cycle = 0; cycle < 3; ++cycle) {
        Handle decoder;
        Check(Phase("create", width, family, cycle, 0, [&] {
            return wirehair_v2_decoder_create(profile, sizeof(profile), &decoder.value);
        }) == WirehairV2_Success && decoder.value, "create");
        bool success = false;
        for (unsigned step = 0; step < 12; ++step) {
            const unsigned slot = step < 6 ? (family == 0 ? 6u : 12u) + step : step - 6u;
            const WirehairV2Result result = Phase("feed", width, family, cycle, Packet(slot), [&] {
                return wirehair_v2_decode(decoder.value, Packet(slot), packets + size_t(slot)*width, width);
            });
            Check(result == WirehairV2_Success || result == WirehairV2_NeedMore, "feed status");
            if (result == WirehairV2_Success) { Check(step >= 5, "premature success"); success = true; break; }
        }
        Check(success, "endpoint");
        memset(recovered, 0xa5, sizeof(recovered)); uint64_t written = 0;
        Check(Phase("recover", width, family, cycle, 0, [&] {
            return wirehair_v2_recover(decoder.value, recovered + 64u, bytes, &written);
        }) == WirehairV2_Success && written == bytes &&
            memcmp(recovered + 64u, source, bytes) == 0, "recovered bytes");
        for (size_t i = 0; i < sizeof(recovered); ++i)
            if (i < 64u || i >= 64u + bytes) Check(recovered[i] == 0xa5, "recovery guard");
        Phase("free", width, family, cycle, 0, [&] { decoder.Free(); return WirehairV2_Success; });
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
            Check(Packet(0) == 0 && Packet(11) == 11 && Packet(12) == UINT32_MAX &&
                Packet(17) == UINT32_MAX-10u, "ID endpoints");
            const uint8_t data[] = {0, 15, 128, 255}; Hex(data, sizeof(data));
            Check(used == 8 && memcmp(output, "000f80ff", 8) == 0, "hex bytes");
            used = sizeof(output)-1; bool capped = false;
            try { Append("xx"); } catch (const std::exception&) { capped = true; }
            Check(capped && used == sizeof(output)-1, "bounded formatting"); used = 0;
            puts("PASS neutral ID/serialization/cap checks; no codec work"); return 0;
        }
        Check(!WH2_RELEASE_NEUTRAL, "neutral scientific worker disabled");
        Check(argc == 3 && strcmp(argv[1], "--worker") == 0 && strlen(argv[2]) == 64, "explicit worker claim");
        for (unsigned i = 0; i < 64; ++i) Check((argv[2][i] >= '0' && argv[2][i] <= '9') ||
            (argv[2][i] >= 'a' && argv[2][i] <= 'f'), "claim hex");
        const rlimit cpu = {10, 10}, memory = {256u*1024u*1024u, 256u*1024u*1024u}, core = {0, 0};
        Check(setrlimit(RLIMIT_CPU, &cpu) == 0 && setrlimit(RLIMIT_AS, &memory) == 0 &&
            setrlimit(RLIMIT_CORE, &core) == 0, "worker limits");
        volatile char* touch = output;
        for (size_t i = 0; i < sizeof(output); ++i) touch[i] = 0;
        Check(wirehair_init() == Wirehair_Success, "GF init");
        Append("{\"type\":\"header\",\"protocol\":\"wirehair.wh2.ordered-rhs-release-r0\","
            "\"claim_sha256\":\"%s\",\"pid\":%ld,\"speed_claimed\":false}\n", argv[2], long(getpid()));
        for (unsigned width : {2u, 64u, 1280u}) Width(width);
        Append("{\"type\":\"footer\",\"outcome\":\"COMPLETE\",\"phases\":%u}\n", phases);
        Write(STDOUT_FILENO, output, used); return 0;
    } catch (const std::exception& error) {
        if (used) { try { Write(STDOUT_FILENO, output, used); } catch (...) {} }
        fprintf(stderr, "INVALID: %s\n", error.what()); return 1;
    }
}
