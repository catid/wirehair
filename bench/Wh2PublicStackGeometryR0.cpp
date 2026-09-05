// Untimed GCC 13/Linux x86-64 compiler/ABI preflight; never a codec benchmark.
// Build: /usr/bin/g++-13 -std=c++11 -O3 -DNDEBUG -Wall -Wextra -Wpedantic
//        -Werror bench/Wh2PublicStackGeometryR0.cpp -o <temporary>/stack_geometry
// Inspect the exact executable's disassembly before invoking --selftest.
#include "../include/wirehair/wirehair.h"
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>

#if defined(__linux__) && defined(__x86_64__) && !defined(__ILP32__) && \
    defined(__GNUC__) && __GNUC__ == 13 && !defined(__clang__)

#define NO_IPA __attribute__((noipa, noinline, noclone))
#define READ_RSP(value) __asm__ __volatile__("mov %%rsp, %0" : "=r"(value) : : "memory")

namespace {
constexpr std::uintptr_t kPageMask = 4095;
constexpr std::uintptr_t kStackLimit = 16384;
constexpr unsigned kCalls = 32;
using Encode = decltype(&wirehair_v2_encode);

struct MockState {
    std::uintptr_t expected_entry;
    std::uintptr_t observed_entry;
    uint32_t* counter;
    unsigned calls;
    unsigned error;
    uint8_t output;
};

struct Evidence {
    std::uintptr_t pre_rsp;
    std::uintptr_t call_rsp;
    std::uintptr_t final_rsp;
    std::uintptr_t counter_address;
    std::size_t correction;
    unsigned error;
};
}

// The reviewed GCC 13 body must have no stack adjustment before this RSP read:
// observed_entry = ordinary public CALL-site RSP - sizeof(return address).
extern "C" NO_IPA WirehairV2Result MockPublicEncode(
    WirehairV2Codec codec, uint32_t id, void* output, uint32_t capacity,
    uint32_t* written) noexcept
{
    std::uintptr_t entry;
    READ_RSP(entry);
    MockState* state = reinterpret_cast<MockState*>(codec);
    state->observed_entry = entry;
    if (entry != state->expected_entry || id != state->calls ||
        id >= kCalls || output != &state->output || capacity != 1 ||
        written != state->counter) {
        state->error = 1;
        return WirehairV2_Error;
    }
    state->output = static_cast<uint8_t>(id ^ 0xa5u);
    *written = 1;
    ++state->calls;
    return WirehairV2_Success;
}

// One ordinary indirect public-ABI CALL site.  RSP is only read by assembly;
// its movement and restoration belong entirely to the compiler.
extern "C" NO_IPA bool StackBatch(
    unsigned target, Encode encode, MockState* state, Evidence* evidence) noexcept
{
    if (target > kPageMask || (target & 15u) != 0 || !encode || !state || !evidence)
        return false;
    uint32_t written = 0;
    evidence->counter_address = reinterpret_cast<std::uintptr_t>(&written);
    std::uintptr_t pre;
    READ_RSP(pre);
    evidence->pre_rsp = pre;
    // Reviewed GCC 13.3: actual decrement D = (requested_bytes + 23) & ~15.
    // Thus request D-16 for aligned D.  A whole extra page makes this one
    // strictly downward correction; stack-clash page probes remain enabled.
    const std::size_t correction = 4096u + ((pre - target) & kPageMask);
    evidence->correction = correction;
    if ((pre & 15u) != 0 || correction > kStackLimit) {
        evidence->error = 1;
        return false;
    }
    // Alignment is in BITS: 128 bits = 16 bytes.  Same lexical scope until
    // final_rsp is captured; the pointer remains live through the final barrier.
    const std::size_t allocation_bytes = correction - 16;
    volatile uint8_t* padding = static_cast<volatile uint8_t*>(
        __builtin_alloca_with_align(allocation_bytes, 128));
    __asm__ __volatile__("" : : "r"(padding) : "memory");
    std::uintptr_t actual;
    READ_RSP(actual);
    evidence->call_rsp = actual;
    if (actual > pre || pre - actual != correction ||
        (actual & kPageMask) != target || pre - actual > kStackLimit) {
        evidence->error = 2;
        return false;
    }
    // Prefault bounded padding outside the prospective batch loop.
    for (std::size_t i = 0; i < allocation_bytes; i += 64)
        padding[i] = static_cast<uint8_t>(i);
    padding[allocation_bytes - 1] = 0;
    state->expected_entry = actual - sizeof(std::uintptr_t);
    state->counter = &written;
    state->calls = 0;
    state->error = 0;
    for (unsigned i = 0; i < kCalls; ++i) {
        std::uintptr_t before, after;
        READ_RSP(before);
        const WirehairV2Result result = encode(
            reinterpret_cast<WirehairV2Codec>(state), i, &state->output, 1, &written);
        READ_RSP(after);
        if (before != actual || after != actual || result != WirehairV2_Success ||
            state->error || state->calls != i + 1 || written != 1 ||
            state->output != static_cast<uint8_t>(i ^ 0xa5u)) {
            evidence->error = 3;
            return false;
        }
    }
    __asm__ __volatile__("" : : "r"(padding) : "memory");
    std::uintptr_t final_rsp;
    READ_RSP(final_rsp);
    evidence->final_rsp = final_rsp;
    if (final_rsp != actual) {
        evidence->error = 4;
        return false;
    }
    return true;
}

int main(int argc, char** argv)
{
    if (argc != 2 || std::strcmp(argv[1], "--selftest") != 0) {
        std::fputs("INVALID: expected --selftest only\n", stderr);
        return 2;
    }
    MockState state = {};
    Evidence guard = {};
    if (StackBatch(1, &MockPublicEncode, &state, &guard) ||
        StackBatch(4096, &MockPublicEncode, &state, &guard) ||
        StackBatch(0, nullptr, &state, &guard) ||
        StackBatch(0, &MockPublicEncode, nullptr, &guard) ||
        StackBatch(0, &MockPublicEncode, &state, nullptr) || state.calls != 0 ||
        state.observed_entry != 0 || state.error != 0) {
        std::fputs("INVALID: argument guard accepted input or invoked mock\n", stderr);
        return 1;
    }
    Evidence reference = {};
    std::uintptr_t caller_reference = 0;
    std::size_t minimum = kStackLimit, maximum = 0;
    unsigned total_calls = 0;
    for (unsigned sweep = 0; sweep < 2; ++sweep) {
        for (unsigned target = 0; target < 4096; target += 16) {
            Evidence evidence = {};
            std::uintptr_t caller_before, caller_after;
            READ_RSP(caller_before);
            const bool okay = StackBatch(target, &MockPublicEncode, &state, &evidence);
            READ_RSP(caller_after);
            if (sweep == 0 && target == 0) {
                reference = evidence;
                caller_reference = caller_before;
            }
            if (!okay || evidence.error || state.error || state.calls != kCalls ||
                caller_before != caller_after || caller_before != caller_reference ||
                evidence.pre_rsp != reference.pre_rsp ||
                evidence.counter_address != reference.counter_address ||
                evidence.call_rsp != evidence.final_rsp ||
                state.observed_entry + 8 != evidence.call_rsp ||
                caller_before < state.observed_entry ||
                caller_before - state.observed_entry > kStackLimit) {
                std::fprintf(stderr,
                    "INVALID: sweep=%u target=%u error=%u mock_error=%u calls=%u "
                    "actual_phase=%" PRIuPTR " correction=%zu\n",
                    sweep, target, evidence.error, state.error, state.calls,
                    evidence.call_rsp & kPageMask, evidence.correction);
                return 1;
            }
            if (evidence.correction < minimum) minimum = evidence.correction;
            if (evidence.correction > maximum) maximum = evidence.correction;
            total_calls += state.calls;
        }
    }
    const int printed = std::printf(
        "{\"outcome\":\"PASS\",\"sweeps\":2,\"targets\":256,\"calls_per_target\":32,"
        "\"mock_calls\":%u,\"minimum_correction\":%zu,\"maximum_correction\":%zu,"
        "\"fixed_frame_bytes\":%" PRIuPTR ",\"counter_from_pre_rsp\":%" PRIuPTR ","
        "\"argument_guards\":5,\"mock_prologue_bytes\":0,\"codec_calls\":0,"
        "\"timing_claimed\":false}\n",
        total_calls, minimum, maximum, caller_reference - reference.pre_rsp,
        reference.counter_address - reference.pre_rsp);
    const int flushed = std::fflush(stdout);
    if (printed < 0 || flushed != 0) {
        std::fputs("INVALID: cannot publish result\n", stderr);
        return 1;
    }
    return 0;
}

#else
int main()
{
    std::fputs("INVALID: requires GCC 13 on Linux x86-64 LP64\n", stderr);
    return 2;
}
#endif
