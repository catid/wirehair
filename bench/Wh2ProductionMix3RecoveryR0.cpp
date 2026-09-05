// Development-only finite seed-map falsifier.  No production equation changes.
#include "Wh2FrozenTrace.h"
#include <wirehair/wirehair.h>

#include <array>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#error "Compile with the exact committed source identity"
#endif

namespace {
using namespace wirehair::wh2_benchmark;
using Clock = std::chrono::steady_clock;
const char kProtocol[] = "wirehair.wh2.production-mix3-k3k6-r0";
const std::array<uint32_t, 2> kKs = {{3u, 6u}};
const std::array<uint32_t, 4> kOverheads = {{0u, 1u, 2u, 4u}};
const std::vector<uint64_t> kTraining = {
    UINT64_C(0x7ccd510f122fc160), UINT64_C(0xb889883a79549774),
    UINT64_C(0xb5666de0987896af), UINT64_C(0x30a0ac4ab2e861cc),
    UINT64_C(0x20e6b10b1cc3838e), UINT64_C(0x894e0a8fedcd6cb5)};
const std::vector<uint64_t> kHoldout = {
    UINT64_C(0x688e4a7ca826b448), UINT64_C(0xa2c6fb6d887f8efe),
    UINT64_C(0xe07d30ba3a9d921f), UINT64_C(0x958b941967c7fbba),
    UINT64_C(0x6f00ea1b271de4eb), UINT64_C(0x754a9eb98cd69323),
    UINT64_C(0x3f3f217a15779be6), UINT64_C(0x8307ffc9acf675e0),
    UINT64_C(0xb7af2e5f6ccb44fa)};
const std::array<FrozenSchedule, 3> kSchedules = {{
    FrozenSchedule::Burst, FrozenSchedule::Adversarial, FrozenSchedule::RepairOnly}};

void Require(bool value, const char* message)
{
    if (!value) throw std::runtime_error(message);
}
std::string Quoted(const std::string& value)
{
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32u) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c);
        else out << c;
    }
    out << '"';
    return out.str();
}
std::string Hex(const uint8_t* data, size_t bytes)
{
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (size_t i = 0; i < bytes; ++i) out << std::setw(2) << unsigned(data[i]);
    return out.str();
}
std::string RootHex(uint64_t root)
{
    std::ostringstream out;
    out << "0x" << std::hex << std::setfill('0') << std::setw(16) << root;
    return out.str();
}
std::string LibraryPath()
{
    const std::array<void*, 10> functions = {{
        reinterpret_cast<void*>(&wirehair_v2_encode),
        reinterpret_cast<void*>(&wirehair_v2_profile_serialize),
        reinterpret_cast<void*>(&wirehair_v2_profile_deserialize),
        reinterpret_cast<void*>(&wirehair_v2_encoder_create),
        reinterpret_cast<void*>(&wirehair_v2_encoder_create_profile),
        reinterpret_cast<void*>(&wirehair_v2_decoder_create),
        reinterpret_cast<void*>(&wirehair_v2_decode),
        reinterpret_cast<void*>(&wirehair_v2_recover),
        reinterpret_cast<void*>(&wirehair_v2_free),
        reinterpret_cast<void*>(&wirehair_init_)}};
    std::string path;
    for (void* function : functions) {
        Dl_info info = {};
        char resolved[PATH_MAX];
        Require(dladdr(function, &info) != 0 && info.dli_fname &&
            realpath(info.dli_fname, resolved), "library identity");
        if (path.empty()) path = resolved;
        Require(path == resolved, "split public API binding");
    }
    return path;
}
struct Handle
{
    WirehairV2Codec Value = nullptr;
    ~Handle() { wirehair_v2_free(Value); }
    Handle() = default;
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
};
std::vector<uint8_t> Source(uint32_t K)
{
    std::vector<uint8_t> source(K * 2u);
    for (size_t i = 0; i < source.size(); ++i)
        source[i] = static_cast<uint8_t>(73u * i + 19u * K + 11u);
    return source;
}
using ProfileBytes = std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES>;
ProfileBytes ExactProfile(uint32_t K, uint32_t attempt)
{
    Require((K == 3u || K == 6u) && attempt <= 255u, "exact profile range");
    WirehairV2Profile profile = {};
    profile.struct_bytes = sizeof(profile);
    profile.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    profile.profile_id = WIREHAIR_V2_PROFILE_CERTIFIED_2026_07;
    profile.message_bytes = K * 2u;
    profile.block_bytes = 2u;
    profile.seed_attempt = static_cast<uint8_t>(attempt);
    ProfileBytes bytes = {{0u}};
    uint32_t written = 0u;
    Require(wirehair_v2_profile_serialize(&profile, bytes.data(), bytes.size(), &written) ==
        WirehairV2_Success && written == bytes.size(), "serialize exact profile");
    return bytes;
}
uint32_t ProfileAttempt(const ProfileBytes& bytes, uint32_t K)
{
    WirehairV2Profile profile = {};
    Require(wirehair_v2_profile_deserialize(bytes.data(), bytes.size(), &profile) ==
        WirehairV2_Success && profile.profile_id == WIREHAIR_V2_PROFILE_CERTIFIED_2026_07 &&
        profile.message_bytes == K * 2u && profile.block_bytes == 2u &&
        bytes == ExactProfile(K, profile.seed_attempt), "profile identity");
    return profile.seed_attempt;
}
void VerifySystematic(Handle& encoder, const std::vector<uint8_t>& source)
{
    for (uint32_t id = 0; id < source.size() / 2u; ++id) {
        uint8_t packet[2] = {};
        uint32_t written = 0u;
        Require(wirehair_v2_encode(encoder.Value, id, packet, sizeof(packet), &written) ==
            WirehairV2_Success && written == 2u &&
            std::memcmp(packet, source.data() + id * 2u, 2u) == 0, "systematic bytes");
    }
}
void Guard(const Clock::time_point& deadline)
{
    Require(Clock::now() < deadline, "worker deadline");
}

// A row contains every nested OH outcome; a new decoder is used for each prefix.
bool EmitCell(const std::string& phase, const std::string& arm, uint32_t K,
    const ProfileBytes& profile, Handle& encoder, bool bad_seed,
    const std::vector<uint8_t>& source, uint32_t root_index, uint64_t root,
    FrozenSchedule schedule, const Clock::time_point& deadline, uint32_t& rows)
{
    Guard(deadline);
    FrozenPacketTrace trace;
    Require(GenerateFrozenPacketTrace(K, 2u, 500000u, root, schedule, trace) ==
        FrozenTraceStatus::Complete && trace.delivered_ids.size() == K + 4u,
        "frozen trace generation");
    std::vector<uint8_t> packets;
    if (!bad_seed) {
        packets.resize(trace.delivered_ids.size() * 2u);
        for (size_t i = 0; i < trace.delivered_ids.size(); ++i) {
            uint32_t written = 0u;
            Require(wirehair_v2_encode(encoder.Value, trace.delivered_ids[i],
                packets.data() + i * 2u, 2u, &written) == WirehairV2_Success &&
                written == 2u, "packet encode");
        }
    }
    std::array<int, 4> outcomes = {{WirehairV2_BadSeed, WirehairV2_BadSeed,
        WirehairV2_BadSeed, WirehairV2_BadSeed}};
    for (size_t h = 0; !bad_seed && h < kOverheads.size(); ++h) {
        Guard(deadline);
        Handle decoder;
        Require(wirehair_v2_decoder_create(profile.data(), profile.size(), &decoder.Value) ==
            WirehairV2_Success && decoder.Value, "decoder create");
        WirehairV2Result result = WirehairV2_NeedMore;
        for (uint32_t i = 0; i < K + kOverheads[h] && result == WirehairV2_NeedMore; ++i) {
            result = wirehair_v2_decode(decoder.Value, trace.delivered_ids[i],
                packets.data() + i * 2u, 2u);
            Require(result == WirehairV2_Success || result == WirehairV2_NeedMore,
                "unexpected decoder result");
        }
        outcomes[h] = result;
        if (result == WirehairV2_Success) {
            std::vector<uint8_t> recovered(source.size(), 0xa5u);
            uint64_t bytes = 0u;
            Require(wirehair_v2_recover(decoder.Value, recovered.data(), recovered.size(), &bytes) ==
                WirehairV2_Success && bytes == source.size() && recovered == source,
                "recovered bytes mismatch");
        }
    }
    Require(source == Source(K), "source mutated during packet/decode work");
    std::cout << "{\"type\":\"cell\",\"phase\":" << Quoted(phase)
        << ",\"arm\":" << Quoted(arm) << ",\"K\":" << K
        << ",\"attempt\":" << ProfileAttempt(profile, K)
        << ",\"root_index\":" << root_index << ",\"root\":" << Quoted(RootHex(root))
        << ",\"schedule\":" << Quoted(FrozenScheduleName(schedule))
        << ",\"profile_hex\":" << Quoted(Hex(profile.data(), profile.size()))
        << ",\"profile_sha256\":" << Quoted(Sha256Hex(profile.data(), profile.size()))
        << ",\"source_sha256\":" << Quoted(Sha256Hex(source.data(), source.size()))
        << ",\"attempted_candidates\":" << trace.attempted_candidates
        << ",\"trace_sha256\":" << Quoted(trace.trace_sha256) << ",\"ids\":[";
    for (size_t i = 0; i < trace.delivered_ids.size(); ++i)
        std::cout << (i ? "," : "") << trace.delivered_ids[i];
    std::cout << "],\"packets_hex\":" << (bad_seed ? "null" : Quoted(Hex(packets.data(), packets.size())))
        << ",\"outcomes\":[";
    for (size_t i = 0; i < outcomes.size(); ++i) std::cout << (i ? "," : "") << outcomes[i];
    std::cout << "],\"recovered_sha256\":[";
    for (size_t i = 0; i < outcomes.size(); ++i)
        std::cout << (i ? "," : "") << (outcomes[i] == 0 ?
            Quoted(Sha256Hex(source.data(), source.size())) : "null");
    std::cout << "]}\n";
    Require(std::cout.good(), "output stream");
    ++rows;
    return outcomes[0] == WirehairV2_Success;
}

bool EmitArm(const std::string& phase, const std::string& arm, uint32_t K,
    int attempt, const std::vector<uint64_t>& roots, const Clock::time_point& deadline,
    uint32_t& rows, uint32_t& actual_attempt)
{
    Guard(deadline);
    const std::vector<uint8_t> source = Source(K);
    ProfileBytes profile = {{0u}};
    Handle encoder;
    WirehairV2Result result;
    if (attempt < 0) {
        uint32_t bytes = 0u;
        result = wirehair_v2_encoder_create(source.data(), source.size(), 2u,
            profile.data(), profile.size(), &bytes, &encoder.Value);
        Require(result == WirehairV2_Success && bytes == profile.size(), "baseline create");
    } else {
        profile = ExactProfile(K, static_cast<uint32_t>(attempt));
        result = wirehair_v2_encoder_create_profile(source.data(), profile.data(),
            profile.size(), &encoder.Value);
        Require(result == WirehairV2_Success || result == WirehairV2_BadSeed,
            "exact constructor result");
    }
    const bool bad_seed = result == WirehairV2_BadSeed;
    Require(bad_seed ? encoder.Value == nullptr : encoder.Value != nullptr, "constructor handle");
    Require(source == Source(K), "source mutated during construction");
    actual_attempt = ProfileAttempt(profile, K);
    if (!bad_seed) VerifySystematic(encoder, source);
    bool all = true;
    for (uint32_t r = 0; r < roots.size(); ++r)
        for (FrozenSchedule schedule : kSchedules)
            if (!EmitCell(phase, arm, K, profile, encoder, bad_seed, source,
                    r, roots[r], schedule, deadline, rows)) all = false;
    return all;
}
void Begin(const std::string& phase)
{
    std::cout << "{\"type\":\"begin\",\"protocol\":" << Quoted(kProtocol)
        << ",\"phase\":" << Quoted(phase) << ",\"source_commit\":"
        << Quoted(WIREHAIR_WH2_SOURCE_GIT_COMMIT) << ",\"library_path\":"
        << Quoted(LibraryPath()) << "}\n";
}
uint32_t Number(const char* value, uint32_t maximum)
{
    Require(value && *value, "empty integer");
    uint32_t result = 0u;
    for (const char* p = value; *p; ++p) {
        Require(*p >= '0' && *p <= '9', "invalid integer");
        const uint32_t digit = static_cast<uint32_t>(*p - '0');
        Require(result <= maximum / 10u &&
            (result < maximum / 10u || digit <= maximum % 10u), "integer range");
        result = result * 10u + digit;
    }
    return result;
}
void Run(const std::string& phase, uint32_t seconds, const std::array<uint32_t, 2>& map)
{
    const Clock::time_point deadline = Clock::now() + std::chrono::seconds(seconds);
    Require(wirehair_init() == Wirehair_Success, "initialization");
    Begin(phase);
    uint32_t rows = 0u;
    std::array<int, 2> selected = {{-1, -1}};
    for (size_t k = 0; k < kKs.size(); ++k) {
        uint32_t actual = 0u;
        EmitArm(phase, "baseline", kKs[k], -1,
            phase == "train" ? kTraining : kHoldout, deadline, rows, actual);
        if (phase == "train") {
            for (uint32_t attempt = 0u; attempt < 256u; ++attempt) {
                if (EmitArm(phase, "candidate", kKs[k], static_cast<int>(attempt),
                        kTraining, deadline, rows, actual)) {
                    selected[k] = static_cast<int>(attempt);
                    break;
                }
            }
        } else {
            selected[k] = static_cast<int>(map[k]);
            EmitArm(phase, "candidate", kKs[k], selected[k], kHoldout, deadline, rows, actual);
        }
    }
    Guard(deadline);
    std::cout << "{\"type\":\"terminal\",\"phase\":" << Quoted(phase)
        << ",\"rows\":" << rows << ",\"selected_attempts\":[" << selected[0] << ','
        << selected[1] << "]}\n";
    std::cout.flush();
    Require(std::cout.good(), "terminal stream");
}
void Selftest()
{
    Require(Sha256Hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "SHA selftest");
    Require(wirehair_init() == Wirehair_Success, "selftest init");
    uint32_t rows = 0u, actual = 0u;
    const std::vector<uint64_t> roots = {UINT64_C(0x1234567890abcdef)};
    for (uint32_t K : kKs) {
        EmitArm("selftest", "baseline", K, -1, roots, Clock::time_point::max(), rows, actual);
        EmitArm("selftest", "candidate", K, static_cast<int>(actual), roots,
            Clock::time_point::max(), rows, actual);
        Require(ProfileAttempt(ExactProfile(K, 255u), K) == 255u, "profile255 selftest");
    }
    for (const char* invalid : {"256", "-1", "1x", "", "999999999999999999999"}) {
        bool rejected = false;
        try { (void)Number(invalid, 255u); } catch (const std::runtime_error&) { rejected = true; }
        Require(rejected, "integer rejection selftest");
    }
    bool profile_rejected = false;
    try { (void)ExactProfile(3u, 256u); } catch (const std::runtime_error&) { profile_rejected = true; }
    Require(profile_rejected && rows == 12u && Number("255", 255u) == 255u, "selftest rows/ranges");
    std::cout << "{\"type\":\"selftest\",\"pass\":true}\n";
}
} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc == 2 && std::string(argv[1]) == "--describe") { Begin("describe"); return 0; }
        if (argc == 2 && std::string(argv[1]) == "--selftest") { Selftest(); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--train") {
            const uint32_t seconds = Number(argv[2], 60u);
            Require(seconds != 0u, "zero deadline");
            Run("train", seconds, {{0u, 0u}}); return 0;
        }
        if (argc == 5 && std::string(argv[1]) == "--holdout") {
            const uint32_t seconds = Number(argv[2], 60u);
            Require(seconds != 0u, "zero deadline");
            Run("holdout", seconds, {{Number(argv[3], 255u), Number(argv[4], 255u)}}); return 0;
        }
        throw std::runtime_error("usage: --describe | --selftest | --train seconds | --holdout seconds a3 a6");
    } catch (const std::exception& error) {
        std::cerr << "production MIX3 recovery r0: " << error.what() << '\n';
        return 2;
    }
}
