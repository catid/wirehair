#include "WirehairV2Codec.h"

#include "../gf256.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <vector>

namespace {

enum class State
{
    LegacyEncoder,
    LegacyDecoder,
    PrecodeEncoder,
    PrecodeDecoder
};

const State kStates[] = {
    State::LegacyEncoder,
    State::LegacyDecoder,
    State::PrecodeEncoder,
    State::PrecodeDecoder
};

const char* StateName(State state)
{
    switch (state)
    {
    case State::LegacyEncoder: return "legacy-encoder";
    case State::LegacyDecoder: return "legacy-decoder";
    case State::PrecodeEncoder: return "precode-encoder";
    case State::PrecodeDecoder: return "precode-decoder";
    }
    return "unknown";
}

bool IsEncoder(State state)
{
    return state == State::LegacyEncoder || state == State::PrecodeEncoder;
}

bool IsPrecode(State state)
{
    return state == State::PrecodeEncoder || state == State::PrecodeDecoder;
}

struct Fixture
{
    uint32_t K;
    uint32_t BlockBytes;
    uint64_t MessageBytes;
    std::vector<uint8_t> Message;

    Fixture(uint32_t k, uint32_t block_bytes, uint32_t short_bytes,
            uint64_t seed)
        : K(k)
        , BlockBytes(block_bytes)
        , MessageBytes((uint64_t)k * block_bytes - short_bytes)
        , Message((size_t)MessageBytes)
    {
        uint64_t x = seed;
        for (size_t i = 0; i < Message.size(); ++i)
        {
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            Message[i] = (uint8_t)(x * UINT64_C(2685821657736338717));
        }
    }

    uint32_t DataBytes(uint32_t block_id) const
    {
        if (block_id + 1u < K) {
            return BlockBytes;
        }
        return (uint32_t)(MessageBytes - (uint64_t)(K - 1u) * BlockBytes);
    }

    const uint8_t* Block(uint32_t block_id) const
    {
        return Message.data() + (size_t)block_id * BlockBytes;
    }
};

bool SameCodecPolicy(
    const wirehair_v2::PeelingCodec& a,
    const wirehair_v2::PeelingCodec& b)
{
    return a.Solver == b.Solver && a.Structure == b.Structure &&
        a.Family == b.Family && a.MinDegree == b.MinDegree &&
        a.MaxDegree == b.MaxDegree &&
        a.SolverCandidateLimit == b.SolverCandidateLimit &&
        a.Degree1Mass == b.Degree1Mass &&
        a.Degree2Mass == b.Degree2Mass && a.RobustC == b.RobustC &&
        a.RobustDelta == b.RobustDelta &&
        a.FullyRandomRows == b.FullyRandomRows &&
        a.UseWirehairRowDistribution == b.UseWirehairRowDistribution;
}

bool SameProfile(
    const wirehair_v2::SeedProfile& a,
    const wirehair_v2::SeedProfile& b)
{
    return a.BlockCount == b.BlockCount && a.BlockBytes == b.BlockBytes &&
        a.Policy.Solver == b.Policy.Solver &&
        a.Policy.Structure == b.Policy.Structure &&
        a.Policy.ByteClass == b.Policy.ByteClass &&
        a.Policy.CountBand == b.Policy.CountBand &&
        SameCodecPolicy(a.Policy.Codec, b.Policy.Codec) &&
        a.DenseCount == b.DenseCount && a.PeelSeed == b.PeelSeed &&
        a.DenseSeed == b.DenseSeed &&
        a.PeelSeedBucket == b.PeelSeedBucket &&
        a.UsedPeelFixup == b.UsedPeelFixup &&
        a.UsedDenseFixup == b.UsedDenseFixup &&
        a.V2SeedSelected == b.V2SeedSelected &&
        a.V2SeedAttempt == b.V2SeedAttempt &&
        a.V2PrecodeContractVersion == b.V2PrecodeContractVersion &&
        a.V2PacketRowContractVersion == b.V2PacketRowContractVersion &&
        a.V2StaircaseCount == b.V2StaircaseCount &&
        a.V2DenseRowCount == b.V2DenseRowCount &&
        a.V2HeavyRowCount == b.V2HeavyRowCount &&
        a.V2CompletionField == b.V2CompletionField &&
        a.V2SourceHits == b.V2SourceHits &&
        a.V2PrecodeSeed == b.V2PrecodeSeed &&
        a.V2PacketPeelSeed == b.V2PacketPeelSeed &&
        a.V2RecoveryMixCount == b.V2RecoveryMixCount &&
        a.V2DenseIdentityCorner == b.V2DenseIdentityCorner &&
        a.V2PrecodeSeedSalt == b.V2PrecodeSeedSalt &&
        a.V2RecoveryRowSeedSalt == b.V2RecoveryRowSeedSalt &&
        a.Tuned == b.Tuned &&
        a.TuningResidualMean == b.TuningResidualMean &&
        a.TuningResidualColumns == b.TuningResidualColumns &&
        a.TuningXorCost == b.TuningXorCost &&
        a.TuningCandidatesRequested == b.TuningCandidatesRequested &&
        a.TuningCandidatesUnique == b.TuningCandidatesUnique &&
        a.TuningCandidatesCompleted == b.TuningCandidatesCompleted &&
        a.TuningTrials == b.TuningTrials;
}

WirehairResult Initialize(
    wirehair_v2::Codec& codec,
    State state,
    const Fixture& fixture,
    const wirehair_v2::SeedProfile* profile = nullptr)
{
    switch (state)
    {
    case State::LegacyEncoder:
        return codec.InitializeEncoder(
            fixture.Message.data(), fixture.MessageBytes,
            fixture.BlockBytes, profile);
    case State::LegacyDecoder:
        return codec.InitializeDecoder(
            fixture.MessageBytes, fixture.BlockBytes, profile);
    case State::PrecodeEncoder:
        return codec.InitializePrecodeEncoder(
            fixture.Message.data(), fixture.MessageBytes,
            fixture.BlockBytes, profile);
    case State::PrecodeDecoder:
        return codec.InitializePrecodeDecoder(
            fixture.MessageBytes, fixture.BlockBytes, profile);
    }
    return Wirehair_Error;
}

bool FeedPrefix(
    wirehair_v2::Codec& codec,
    State state,
    const Fixture& fixture,
    uint32_t prefix)
{
    if (IsEncoder(state)) {
        return prefix == 0u;
    }
    for (uint32_t id = 0; id < prefix; ++id)
    {
        const WirehairResult result = codec.Decode(
            id, fixture.Block(id), fixture.DataBytes(id));
        if (result != Wirehair_NeedMore) {
            std::fprintf(stderr,
                "%s prefix packet %u returned %d\n",
                StateName(state), id, result);
            return false;
        }
    }
    return true;
}

bool InitializeActive(
    wirehair_v2::Codec& codec,
    State state,
    const Fixture& fixture,
    uint32_t prefix)
{
    const WirehairResult result = Initialize(codec, state, fixture);
    if (result != Wirehair_Success)
    {
        std::fprintf(stderr, "%s initialization returned %d\n",
            StateName(state), result);
        return false;
    }
    return FeedPrefix(codec, state, fixture, prefix);
}

bool VerifyEncoder(
    wirehair_v2::Codec& candidate,
    State state,
    const Fixture& fixture)
{
    wirehair_v2::Codec oracle;
    if (Initialize(oracle, state, fixture) != Wirehair_Success ||
        !SameProfile(candidate.Profile(), oracle.Profile()))
    {
        return false;
    }

    for (uint32_t id = 0; id < fixture.K + 9u; ++id)
    {
        std::vector<uint8_t> actual(fixture.BlockBytes + 16u, 0xa5u);
        std::vector<uint8_t> expected(fixture.BlockBytes + 16u, 0xa5u);
        uint32_t actual_bytes = UINT32_MAX;
        uint32_t expected_bytes = UINT32_MAX;
        const WirehairResult actual_result = candidate.Encode(
            id, actual.data(), fixture.BlockBytes, &actual_bytes);
        const WirehairResult expected_result = oracle.Encode(
            id, expected.data(), fixture.BlockBytes, &expected_bytes);
        if (actual_result != expected_result ||
            actual_bytes != expected_bytes || actual != expected)
        {
            std::fprintf(stderr,
                "%s packet %u changed after failed transition\n",
                StateName(state), id);
            return false;
        }
    }
    return true;
}

bool VerifyDecoder(
    wirehair_v2::Codec& codec,
    State state,
    const Fixture& fixture,
    uint32_t prefix)
{
    const State encoder_state = state == State::LegacyDecoder ?
        State::LegacyEncoder : State::PrecodeEncoder;
    wirehair_v2::Codec encoder;
    wirehair_v2::Codec oracle;
    if (Initialize(encoder, encoder_state, fixture) != Wirehair_Success ||
        Initialize(oracle, state, fixture) != Wirehair_Success ||
        !SameProfile(codec.Profile(), oracle.Profile()) ||
        !FeedPrefix(oracle, state, fixture, prefix))
    {
        return false;
    }

    bool decoded = false;
    for (uint32_t packet_i = 0; packet_i < fixture.K + 64u; ++packet_i)
    {
        const uint32_t block_id = fixture.K + packet_i;
        std::vector<uint8_t> block(fixture.BlockBytes, 0u);
        uint32_t data_bytes = 0u;
        if (encoder.Encode(
                block_id, block.data(), fixture.BlockBytes, &data_bytes) !=
                Wirehair_Success || data_bytes != fixture.BlockBytes)
        {
            std::fprintf(stderr,
                "%s source packet %u failed\n", StateName(state), block_id);
            return false;
        }
        const WirehairResult actual =
            codec.Decode(block_id, block.data(), data_bytes);
        const WirehairResult expected =
            oracle.Decode(block_id, block.data(), data_bytes);
        if (actual != expected ||
            (actual != Wirehair_NeedMore && actual != Wirehair_Success))
        {
            std::fprintf(stderr,
                "%s packet %u returned %d, oracle returned %d\n",
                StateName(state), block_id, actual, expected);
            return false;
        }
        if (actual == Wirehair_Success) {
            decoded = true;
            break;
        }
    }
    if (!decoded)
    {
        std::fprintf(stderr, "%s did not decode bounded repair stream\n",
            StateName(state));
        return false;
    }

    std::vector<uint8_t> recovered(fixture.Message.size(), 0x5au);
    std::vector<uint8_t> expected_message(fixture.Message.size(), 0xa5u);
    const WirehairResult result =
        codec.Recover(recovered.data(), fixture.MessageBytes);
    const WirehairResult expected_result =
        oracle.Recover(expected_message.data(), fixture.MessageBytes);
    if (result != expected_result || result != Wirehair_Success ||
        recovered != expected_message || recovered != fixture.Message)
    {
        std::fprintf(stderr, "%s recovery changed after transition\n",
            StateName(state));
        return false;
    }
    return true;
}

bool VerifyOperational(
    wirehair_v2::Codec& codec,
    State state,
    const Fixture& fixture,
    uint32_t prefix)
{
    return IsEncoder(state) ? VerifyEncoder(codec, state, fixture) :
        VerifyDecoder(codec, state, fixture, prefix);
}

bool VerifyPreserved(
    wirehair_v2::Codec& codec,
    State state,
    const Fixture& fixture,
    uint32_t prefix,
    const wirehair_v2::SeedProfile& before)
{
    if (!SameProfile(codec.Profile(), before))
    {
        std::fprintf(stderr, "%s profile changed on failed transition\n",
            StateName(state));
        return false;
    }
    return VerifyOperational(codec, state, fixture, prefix);
}

bool CheckSuccessfulTransitions(const Fixture& from, const Fixture& to)
{
    for (State old_state : kStates)
    {
        for (State next_state : kStates)
        {
            wirehair_v2::Codec codec;
            const uint32_t old_prefix = IsEncoder(old_state) ? 0u : 2u;
            if (!InitializeActive(codec, old_state, from, old_prefix) ||
                Initialize(codec, next_state, to) != Wirehair_Success ||
                codec.Profile().BlockCount != to.K ||
                codec.Profile().BlockBytes != to.BlockBytes ||
                !VerifyOperational(codec, next_state, to, 0u))
            {
                std::fprintf(stderr, "successful transition %s -> %s failed\n",
                    StateName(old_state), StateName(next_state));
                return false;
            }
        }
    }
    return true;
}

bool CheckValidationRollback(const Fixture& active, const Fixture& target)
{
    wirehair_v2::SeedProfile mismatch =
        wirehair_v2::SelectSeedProfile(target.K + 1u, target.BlockBytes);
    for (State old_state : kStates)
    {
        for (State target_state : kStates)
        {
            wirehair_v2::Codec codec;
            const uint32_t prefix = IsEncoder(old_state) ? 0u : 2u;
            if (!InitializeActive(codec, old_state, active, prefix)) {
                return false;
            }
            const wirehair_v2::SeedProfile before = codec.Profile();
            const WirehairResult result =
                Initialize(codec, target_state, target, &mismatch);
            if (result != Wirehair_InvalidInput ||
                !VerifyPreserved(codec, old_state, active, prefix, before))
            {
                std::fprintf(stderr, "validation rollback %s -> %s failed\n",
                    StateName(old_state), StateName(target_state));
                return false;
            }
        }
    }
    return true;
}

bool CheckCoreSetupRollback(const Fixture& active, const Fixture& target)
{
    wirehair_v2::SeedProfile invalid =
        wirehair_v2::SelectSeedProfile(target.K, target.BlockBytes);
    invalid.DenseCount = 0u;
    for (State old_state : kStates)
    {
        for (State target_state :
            { State::LegacyEncoder, State::LegacyDecoder })
        {
            wirehair_v2::Codec codec;
            const uint32_t prefix = IsEncoder(old_state) ? 0u : 2u;
            if (!InitializeActive(codec, old_state, active, prefix)) {
                return false;
            }
            const wirehair_v2::SeedProfile before = codec.Profile();
            const WirehairResult result =
                Initialize(codec, target_state, target, &invalid);
            if (result != Wirehair_InvalidInput ||
                !VerifyPreserved(codec, old_state, active, prefix, before))
            {
                std::fprintf(stderr, "core rollback %s -> %s failed\n",
                    StateName(old_state), StateName(target_state));
                return false;
            }
        }
    }
    return true;
}

bool CheckPrecodeSetupRollback(const Fixture& active, const Fixture& target)
{
    wirehair_v2::SeedProfile invalid =
        wirehair_v2::SelectSeedProfile(target.K, target.BlockBytes);
    // Contract metadata is all-or-nothing.  This passes facade dimensions,
    // then fails the precode seed/configuration validation in the candidate.
    invalid.V2PrecodeContractVersion = 1u;
    for (State old_state : kStates)
    {
        for (State target_state :
            { State::PrecodeEncoder, State::PrecodeDecoder })
        {
            wirehair_v2::Codec codec;
            const uint32_t prefix = IsEncoder(old_state) ? 0u : 2u;
            if (!InitializeActive(codec, old_state, active, prefix)) {
                return false;
            }
            const wirehair_v2::SeedProfile before = codec.Profile();
            const WirehairResult result =
                Initialize(codec, target_state, target, &invalid);
            if (result != Wirehair_InvalidInput ||
                !VerifyPreserved(codec, old_state, active, prefix, before))
            {
                std::fprintf(stderr, "precode rollback %s -> %s failed\n",
                    StateName(old_state), StateName(target_state));
                return false;
            }
        }
    }
    return true;
}

enum class AllocationHook
{
    Facade,
    PrecodeEncoder,
    PrecodeDecoder
};

void SetHook(AllocationHook hook, int64_t countdown)
{
    switch (hook)
    {
    case AllocationHook::Facade:
        wirehair_v2::SetCodecAllocationFailureCountdownForTesting(countdown);
        break;
    case AllocationHook::PrecodeEncoder:
        wirehair_v2::SetAllocationFailureCountdownForTesting(countdown);
        break;
    case AllocationHook::PrecodeDecoder:
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(countdown);
        break;
    }
}

unsigned ExpectedAllocationFailures(AllocationHook hook, State target_state)
{
    if (hook == AllocationHook::Facade) {
        // Legacy: facade object + two core allocations.  Precode: facade
        // object only; the implementation allocations have their own hook.
        return IsPrecode(target_state) ? 1u : 3u;
    }
    // The solved precode encoder validates but no longer copies the complete
    // nested row graph during ownership transfer, removing one real guarded
    // allocation point from this transactional sweep.
    return 4u;
}

bool CheckAllocationHook(
    AllocationHook hook,
    State old_state,
    State target_state,
    const Fixture& active,
    const Fixture& target)
{
    unsigned failures = 0u;
    for (int64_t countdown = 0; countdown < 32; ++countdown)
    {
        wirehair_v2::Codec codec;
        const uint32_t prefix = IsEncoder(old_state) ? 0u : 2u;
        if (!InitializeActive(codec, old_state, active, prefix)) {
            return false;
        }
        const wirehair_v2::SeedProfile before = codec.Profile();
        SetHook(hook, countdown);
        const WirehairResult result = Initialize(codec, target_state, target);
        SetHook(hook, -1);

        if (result == Wirehair_OOM)
        {
            ++failures;
            if (!VerifyPreserved(codec, old_state, active, prefix, before)) {
                return false;
            }
            continue;
        }
        const unsigned expected_failures =
            ExpectedAllocationFailures(hook, target_state);
        if (result != Wirehair_Success || failures != expected_failures ||
            !VerifyOperational(codec, target_state, target, 0u))
        {
            std::fprintf(stderr,
                "%s allocation sweep %s -> %s ended with %d after %u/%u failures\n",
                hook == AllocationHook::Facade ? "facade" :
                    (hook == AllocationHook::PrecodeEncoder ?
                        "precode encoder" : "precode decoder"),
                StateName(old_state), StateName(target_state),
                result, failures, expected_failures);
            return false;
        }
        return true;
    }
    std::fprintf(stderr, "allocation sweep did not reach success\n");
    SetHook(hook, -1);
    return false;
}

bool CheckAllocationRollback(const Fixture& active, const Fixture& target)
{
    for (State old_state : kStates)
    {
        for (State target_state : kStates)
        {
            if (!CheckAllocationHook(
                    AllocationHook::Facade, old_state, target_state,
                    active, target))
            {
                return false;
            }
            if (target_state == State::PrecodeEncoder &&
                !CheckAllocationHook(
                    AllocationHook::PrecodeEncoder, old_state, target_state,
                    active, target))
            {
                return false;
            }
            if (target_state == State::PrecodeDecoder &&
                !CheckAllocationHook(
                    AllocationHook::PrecodeDecoder, old_state, target_state,
                    active, target))
            {
                return false;
            }
        }
    }
    return true;
}

uint64_t NextRandom(uint64_t& x)
{
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    return x * UINT64_C(2685821657736338717);
}

bool StressIteration(uint64_t& random, unsigned iteration)
{
    const uint32_t old_k = 4u + (uint32_t)(NextRandom(random) % 64u);
    const uint32_t old_bb = 8u + (uint32_t)(NextRandom(random) % 128u);
    const uint32_t new_k = 4u + (uint32_t)(NextRandom(random) % 64u);
    const uint32_t new_bb = 8u + (uint32_t)(NextRandom(random) % 128u);
    Fixture active(old_k, old_bb,
        1u + (uint32_t)(NextRandom(random) % (old_bb - 1u)),
        random ^ iteration);
    Fixture target(new_k, new_bb,
        1u + (uint32_t)(NextRandom(random) % (new_bb - 1u)),
        random ^ UINT64_C(0x6a09e667f3bcc909));
    const State old_state = kStates[NextRandom(random) % 4u];
    const State target_state = kStates[NextRandom(random) % 4u];
    const uint32_t prefix = IsEncoder(old_state) ? 0u : 2u;

    wirehair_v2::Codec codec;
    if (!InitializeActive(codec, old_state, active, prefix)) {
        return false;
    }
    const wirehair_v2::SeedProfile before = codec.Profile();

    const unsigned action = (unsigned)(NextRandom(random) % 4u);
    if (action == 0u)
    {
        if (Initialize(codec, target_state, target) != Wirehair_Success) {
            return false;
        }
        return VerifyOperational(codec, target_state, target, 0u);
    }

    WirehairResult result = Wirehair_Error;
    if (action == 1u)
    {
        wirehair_v2::SeedProfile mismatch =
            wirehair_v2::SelectSeedProfile(target.K + 1u, target.BlockBytes);
        result = Initialize(codec, target_state, target, &mismatch);
        if (result != Wirehair_InvalidInput) {
            return false;
        }
    }
    else if (action == 2u || !IsPrecode(target_state))
    {
        const int64_t maximum = IsPrecode(target_state) ? 0 : 2;
        const int64_t countdown =
            maximum == 0 ? 0 : (int64_t)(NextRandom(random) % 3u);
        SetHook(AllocationHook::Facade, countdown);
        result = Initialize(codec, target_state, target);
        SetHook(AllocationHook::Facade, -1);
        if (result != Wirehair_OOM) {
            return false;
        }
    }
    else
    {
        const AllocationHook hook = target_state == State::PrecodeEncoder ?
            AllocationHook::PrecodeEncoder : AllocationHook::PrecodeDecoder;
        SetHook(hook, (int64_t)(NextRandom(random) %
            ExpectedAllocationFailures(hook, target_state)));
        result = Initialize(codec, target_state, target);
        SetHook(hook, -1);
        if (result != Wirehair_OOM) {
            return false;
        }
    }
    return VerifyPreserved(codec, old_state, active, prefix, before);
}

bool RunStress(unsigned thread_count, unsigned iterations)
{
    std::atomic<bool> failed(false);
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (unsigned thread_i = 0; thread_i < thread_count; ++thread_i)
    {
        threads.emplace_back([&, thread_i]() {
            uint64_t random = UINT64_C(0x243f6a8885a308d3) ^
                ((uint64_t)thread_i * UINT64_C(0x9e3779b97f4a7c15));
            for (unsigned i = 0; i < iterations && !failed.load(); ++i)
            {
                if (!StressIteration(random, i))
                {
                    std::fprintf(stderr,
                        "state-machine stress failed in thread %u iteration %u\n",
                        thread_i, i);
                    failed.store(true);
                    return;
                }
            }
        });
    }
    for (std::thread& thread : threads) {
        thread.join();
    }
    return !failed.load();
}

bool ParseUnsigned(const char* text, unsigned limit, unsigned& value)
{
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(text, &end, 10);
    if (!text[0] || !end || *end || parsed == 0u || parsed > limit) {
        return false;
    }
    value = (unsigned)parsed;
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    unsigned threads = 8u;
    unsigned iterations = 20u;
    if (argc == 4 && std::strcmp(argv[1], "--stress") == 0)
    {
        if (!ParseUnsigned(argv[2], 1024u, threads) ||
            !ParseUnsigned(argv[3], 1000000u, iterations))
        {
            std::fprintf(stderr, "usage: %s [--stress THREADS ITERATIONS]\n",
                argv[0]);
            return 2;
        }
    }
    else if (argc != 1)
    {
        std::fprintf(stderr, "usage: %s [--stress THREADS ITERATIONS]\n",
            argv[0]);
        return 2;
    }

    if (gf256_init() != 0) {
        return 1;
    }
    const Fixture active(8u, 31u, 7u, UINT64_C(0x123456789abcdef0));
    const Fixture target(11u, 43u, 13u, UINT64_C(0xfedcba9876543210));

    const bool ok =
        CheckSuccessfulTransitions(active, target) &&
        CheckValidationRollback(active, target) &&
        CheckCoreSetupRollback(active, target) &&
        CheckPrecodeSetupRollback(active, target) &&
        CheckAllocationRollback(active, target) &&
        RunStress(threads, iterations);
    if (!ok) {
        return 1;
    }
    std::printf(
        "transactional facade: PASS (%u threads x %u iterations)\n",
        threads, iterations);
    return 0;
}
