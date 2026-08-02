#ifndef WIREHAIR_BENCH_WH2_NATIVE_PANEL_H
#define WIREHAIR_BENCH_WH2_NATIVE_PANEL_H

#include <array>
#include <functional>
#include <memory>
#include <stdint.h>
#include <string>

namespace wirehair_wh2_bench {

enum class NativePanelOrder : uint32_t
{
    ABBA = 0,
    BAAB
};

enum class NativePanelSide : uint32_t
{
    Left = 0,
    Right
};

/**
    Success and PreflightFailure are reproducible codec outcomes.  Fatal is
    reserved for an invocation which cannot provide benchmark evidence.
*/
enum class NativePanelDisposition : uint32_t
{
    Success = 0,
    PreflightFailure,
    Fatal
};

/**
    OutcomeCode is an adapter-defined stable codec result.  HasDecodedExtra
    distinguishes JSON null from a numeric decoded-extra value.  The callback
    owns the scope timer and reports its elapsed interval here; this keeps the
    panel scheduler independent of codec-specific setup and verification.
*/
struct NativePanelInvocationResult
{
    NativePanelDisposition Disposition;
    int64_t OutcomeCode;
    bool HasDecodedExtra;
    uint32_t DecodedExtra;
    uint64_t ElapsedNanoseconds;

    NativePanelInvocationResult();
    NativePanelInvocationResult(
        NativePanelDisposition disposition,
        int64_t outcome_code,
        bool has_decoded_extra,
        uint32_t decoded_extra,
        uint64_t elapsed_nanoseconds);
};

/**
    One freshly-created invocation.  Invoke() deliberately receives no side,
    panel label, order, or slot, so measured code cannot specialize by its
    position in ABBA/BAAB.  Identity() is checked before and after Invoke().
*/
class NativePanelInvocation
{
public:
    virtual ~NativePanelInvocation();
    virtual std::string Identity() const = 0;
    virtual NativePanelInvocationResult Invoke() = 0;
};

typedef std::function<std::unique_ptr<NativePanelInvocation>()>
    NativePanelInvocationFactory;

struct NativePanelArm
{
    std::string ExpectedIdentity;
    NativePanelInvocationFactory MakeInvocation;

    NativePanelArm();
    NativePanelArm(
        const std::string& expected_identity,
        const NativePanelInvocationFactory& make_invocation);
};

enum class NativePanelStatus : uint32_t
{
    Complete = 0,
    UnsupportedPlatform,
    InvalidArgument,
    PinFailed,
    AffinityVerificationFailed,
    CpuQueryFailed,
    CpuMigration,
    FactoryFailed,
    IdentityDrift,
    OutcomeDrift,
    InvocationFailed,
    InvalidElapsed
};

const char* NativePanelStatusName(NativePanelStatus status);

struct NativePanelSlot
{
    NativePanelSide Side;
    NativePanelInvocationResult Invocation;
    /** False is the native equivalent of a null timing receipt value. */
    bool HasElapsedNanoseconds;
    uint64_t ElapsedNanoseconds;

    NativePanelSlot();
};

struct NativePanelResult
{
    NativePanelStatus Status;
    std::string Diagnostic;
    NativePanelOrder Order;
    int TargetCpu;
    bool HasLeftPreflight;
    bool HasRightPreflight;
    NativePanelInvocationResult LeftPreflight;
    NativePanelInvocationResult RightPreflight;
    bool PanelComparable;
    std::array<NativePanelSlot, 4> Slots;

    NativePanelResult();
    bool IsFatal() const;
    bool HasFourNullTimings() const;
};

/**
    OS seam used by the deterministic unit test.  Production callers should
    use ExecuteNativeTimingPanel(), which supplies the Linux implementation.
*/
class NativePanelRuntime
{
public:
    virtual ~NativePanelRuntime();

    virtual NativePanelStatus PinAndVerifySingletonCpu(
        int target_cpu,
        std::string& diagnostic) = 0;

    virtual NativePanelStatus VerifySingletonCpu(
        int target_cpu,
        std::string& diagnostic) = 0;
};

/** True only where the production singleton-affinity path is implemented. */
bool NativePanelPlatformSupported();

/**
    Pin the calling worker to target_cpu, execute warm-left, warm-right, then
    exactly four fresh invocations in the declared chronology, and verify the
    singleton CPU immediately before and after every Invoke().

    If either warmup returns PreflightFailure, all four measured invocations
    are still executed and checked for stable identity/outcome, but all four
    timing flags are cleared.  Migration, affinity drift, identity drift,
    outcome drift, fatal callbacks, and invalid successful timings fail the
    panel and clear every timing flag.
*/
NativePanelResult ExecuteNativeTimingPanel(
    int target_cpu,
    NativePanelOrder order,
    const NativePanelArm& left,
    const NativePanelArm& right);

/** Identical executor with an injected runtime for deterministic tests. */
NativePanelResult ExecuteNativeTimingPanelWithRuntime(
    int target_cpu,
    NativePanelOrder order,
    const NativePanelArm& left,
    const NativePanelArm& right,
    NativePanelRuntime& runtime);

} // namespace wirehair_wh2_bench

#endif // WIREHAIR_BENCH_WH2_NATIVE_PANEL_H
