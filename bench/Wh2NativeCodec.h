#pragma once

#include "../codec/WirehairV2PrecodeEncode.h"

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <memory>
#include <vector>

namespace wirehair_wh2_bench {

/**
    Codec identities understood by the frozen WH2 benchmark worker.

    Wirehair2Experiment deliberately uses the same pure-GF(256) execution
    machinery as the certified control.  Its equation transform is applied to
    the resolved base construction before the public attempt offset, allowing
    one process to execute runtime-selected architecture candidates without
    claiming the certified profile identifier.
*/
enum class NativeArmKind : uint32_t
{
    Invalid = 0,
    Wirehair1,
    Wirehair2Certified,
    Wirehair2Experiment
};

/**
    Base equation source for a native WH2 arm.

    ProductionProfile preserves the public codec path.  The canonical
    certified structure is an experiment-only seam: it starts from the
    certified K-derived geometry with zero precode and packet seeds, so a
    transform can install a closed diagnostic seed schedule without first
    traversing production seed selection.
*/
enum class NativeWh2BaseKind : uint32_t
{
    ProductionProfile = 0,
    CanonicalCertifiedStructure
};

typedef bool (*Wh2EquationTransform)(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& base_params,
    wirehair_v2::PacketRowConfig& base_packet_config,
    void* context);

struct NativeArmSpec
{
    NativeArmKind Kind = NativeArmKind::Invalid;

    /** Controls how the pre-transform WH2 equations are constructed. */
    NativeWh2BaseKind BaseKind = NativeWh2BaseKind::ProductionProfile;

    /**
        Exact WH2 public construction attempt in [0,255].  Wirehair1 has no
        construction-attempt axis and requires zero.
    */
    uint32_t ConstructionAttempt = 0u;

    /**
        Pure-GF(256) message-precode options.  Certified arms require the
        equation-affecting defaults.  Runtime-only facade cache flags must
        remain false because this direct adapter rejects, rather than silently
        ignores, axes that do not affect its storage path.
    */
    wirehair_v2::MessagePrecodeEncoderOptions Wh2Options = {};

    /**
        Required to be null for both controls.  Experimental arms may use a
        deterministic transform to alter the resolved base equations.  One
        NativeArm construction invokes the callback synchronously and does
        not retain it; a timing fixture retains its complete NativeArmSpec as
        documented below.
    */
    Wh2EquationTransform Transform = nullptr;
    void* TransformContext = nullptr;
};

NativeArmSpec MakeWirehair1Arm();
NativeArmSpec MakeCertifiedWh2Arm(uint32_t construction_attempt);
NativeArmSpec MakeExperimentalWh2Arm(
    uint32_t construction_attempt,
    Wh2EquationTransform transform,
    void* transform_context = nullptr,
    const wirehair_v2::MessagePrecodeEncoderOptions* options = nullptr,
    NativeWh2BaseKind base_kind = NativeWh2BaseKind::ProductionProfile);

/**
    Generate an arm-independent, byte-stable message of exactly
    block_count * block_bytes bytes.  SplitMix64 output words are serialized
    little-endian so the same inputs produce identical bytes on every host.
*/
bool MakeDeterministicSource(
    uint32_t block_count,
    uint32_t block_bytes,
    uint64_t source_seed,
    std::vector<uint8_t>& source_out);

class NativeArm
{
public:
    NativeArm();
    ~NativeArm();
    NativeArm(NativeArm&& other) noexcept;
    NativeArm& operator=(NativeArm&& other) noexcept;

    NativeArm(const NativeArm&) = delete;
    NativeArm& operator=(const NativeArm&) = delete;

    /**
        Construct exactly the requested arm once.  WH2 never calls a seed
        selector or retries another attempt.  A rank-deficient exact WH2
        attempt is returned as Wirehair_BadPeelSeed.
    */
    WirehairResult Initialize(
        const NativeArmSpec& spec,
        uint32_t block_count,
        uint32_t block_bytes,
        const std::vector<uint8_t>& source);

    bool IsInitialized() const;
    NativeArmKind Kind() const;
    uint32_t BlockCount() const;
    uint32_t BlockBytes() const;
    uint32_t ConstructionAttempt() const;

    /** Encode one full-size packet into caller-owned storage. */
    WirehairResult EncodeInto(
        uint32_t block_id,
        uint8_t* packet_out,
        uint32_t packet_capacity) const;

    /** Transactional allocating convenience wrapper around EncodeInto(). */
    WirehairResult Encode(
        uint32_t block_id,
        std::vector<uint8_t>& packet_out) const;

private:
    WirehairResult InitializeOwnedSourceAfterGlobalInit(
        const NativeArmSpec& spec,
        uint32_t block_count,
        uint32_t block_bytes,
        std::vector<uint8_t>&& source);

    struct Impl;
    std::shared_ptr<Impl> ImplValue;

    friend class NativeEncoderFixture;
    friend class NativeRecoveryFixture;
    friend class NativeReceiveFixture;
    friend class NativeSolveFixture;
};

enum class RecoveryOutcome : uint32_t
{
    Success = 0,
    NeedMoreAtCap,
    ConstructFailed,
    Unsupported,
    Fatal
};

struct RecoveryCellResult
{
    RecoveryOutcome Outcome = RecoveryOutcome::Fatal;
    /** Valid only for Success; UINT32_MAX otherwise. */
    uint32_t FirstOverhead = UINT32_MAX;
    /** Raw codec result for diagnostics. */
    WirehairResult Result = Wirehair_Error;
};

/**
    Owns the arm-free K+4 packet prefix and all encoded payloads for one cell.
    Distinct packet IDs are required.  WH2 solves fresh nested prefixes K..K+4;
    Wirehair1 likewise creates a fresh decoder for each prefix.
*/
class NativeRecoveryFixture
{
public:
    NativeRecoveryFixture();
    ~NativeRecoveryFixture();
    NativeRecoveryFixture(NativeRecoveryFixture&& other) noexcept;
    NativeRecoveryFixture& operator=(NativeRecoveryFixture&& other) noexcept;

    NativeRecoveryFixture(const NativeRecoveryFixture&) = delete;
    NativeRecoveryFixture& operator=(const NativeRecoveryFixture&) = delete;

    WirehairResult Initialize(
        const NativeArm& arm,
        const std::vector<uint32_t>& packet_ids);

    bool IsInitialized() const;
    RecoveryCellResult RunNested() const;

private:
    struct Impl;
    std::unique_ptr<Impl> ImplValue;
};

/** Build and execute a complete frozen recovery cell, including construction. */
RecoveryCellResult RunRecoveryCell(
    const NativeArmSpec& spec,
    uint32_t block_count,
    uint32_t block_bytes,
    const std::vector<uint8_t>& source,
    const std::vector<uint32_t>& packet_ids);

struct IsolatedSolveResult
{
    WirehairResult Result = Wirehair_Error;
    uint64_t ElapsedNanoseconds = 0u;
    bool BytesVerified = false;
    wirehair_v2::PrecodeSolveStats Stats = {};
};

struct TimedArmResult
{
    WirehairResult Result = Wirehair_Error;
    uint64_t ElapsedNanoseconds = 0u;
    bool BytesVerified = false;
    /** Valid for a successful receive-to-success invocation only. */
    uint32_t DecodedOverhead = UINT32_MAX;
};

/**
    Timing-ready fresh encoder fixture.  Global library initialization,
    fixture ownership of the source bytes, and output-buffer allocation are
    outside the measured interval.  Each Run() constructs a new exact arm and
    encodes exactly IDs 0 through K-1; comparison with the source occurs after
    the clock stops.  A transform context stored in the arm spec must remain
    valid for the lifetime of this fixture.
*/
class NativeEncoderFixture
{
public:
    NativeEncoderFixture();
    ~NativeEncoderFixture();
    NativeEncoderFixture(NativeEncoderFixture&& other) noexcept;
    NativeEncoderFixture& operator=(NativeEncoderFixture&& other) noexcept;

    NativeEncoderFixture(const NativeEncoderFixture&) = delete;
    NativeEncoderFixture& operator=(const NativeEncoderFixture&) = delete;

    WirehairResult Initialize(
        const NativeArmSpec& spec,
        uint32_t block_count,
        uint32_t block_bytes,
        const std::vector<uint8_t>& source);

    bool IsInitialized() const;
    TimedArmResult Run() const;

private:
    struct Impl;
    std::unique_ptr<Impl> ImplValue;
};

/**
    Timing-ready receive-to-success fixture over the exact prepared K+4 trace.
    Each Run() creates fresh decoder bookkeeping before the clock, then times
    all feeds through first success plus full Recover.  WH2 uses the same
    pure-GF(256) cold/resume solver path for certified and runtime candidates;
    Wirehair1 uses its public fresh decoder.  Byte comparison occurs only after
    the clock stops.
*/
class NativeReceiveFixture
{
public:
    NativeReceiveFixture();
    ~NativeReceiveFixture();
    NativeReceiveFixture(NativeReceiveFixture&& other) noexcept;
    NativeReceiveFixture& operator=(NativeReceiveFixture&& other) noexcept;

    NativeReceiveFixture(const NativeReceiveFixture&) = delete;
    NativeReceiveFixture& operator=(const NativeReceiveFixture&) = delete;

    WirehairResult Initialize(
        const NativeArm& arm,
        const std::vector<uint32_t>& packet_ids);

    bool IsInitialized() const;
    TimedArmResult Run() const;

private:
    struct Impl;
    std::unique_ptr<Impl> ImplValue;
};

/**
    Timing-ready WH2 cold-solve fixture.  Construction, packet evaluation,
    packet-pointer assembly, and expected bytes are prepared once outside the
    measured interval.  Each Run() creates fresh decoder output and times only
    SolvePrecodeSystemForValidatedSystemWithRuntime(); byte recovery is checked
    after the clock stops.  Wirehair1 is intentionally rejected because its
    fair cross-codec metric is receive-to-success, not isolated solve.
*/
class NativeSolveFixture
{
public:
    NativeSolveFixture();
    ~NativeSolveFixture();
    NativeSolveFixture(NativeSolveFixture&& other) noexcept;
    NativeSolveFixture& operator=(NativeSolveFixture&& other) noexcept;

    NativeSolveFixture(const NativeSolveFixture&) = delete;
    NativeSolveFixture& operator=(const NativeSolveFixture&) = delete;

    WirehairResult Initialize(
        const NativeArm& arm,
        const std::vector<uint32_t>& packet_ids,
        uint32_t fixed_received_overhead);

    bool IsInitialized() const;
    IsolatedSolveResult Run() const;

private:
    struct Impl;
    std::unique_ptr<Impl> ImplValue;
};

} // namespace wirehair_wh2_bench
