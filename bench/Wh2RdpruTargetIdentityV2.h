#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace wirehair_wh2_bench {

// Raw CPUID register state.  The identity module always serializes these
// fields explicitly and never hashes native object representation.
struct CpuidRegsV2
{
    CpuidRegsV2() = default;
    CpuidRegsV2(uint32_t eax, uint32_t ebx, uint32_t ecx, uint32_t edx)
        : Eax(eax), Ebx(ebx), Ecx(ecx), Edx(edx) {}

    uint32_t Eax = 0u;
    uint32_t Ebx = 0u;
    uint32_t Ecx = 0u;
    uint32_t Edx = 0u;
};

struct CpuidRecordV2
{
    uint32_t Leaf = 0u;
    uint32_t Subleaf = 0u;
    CpuidRegsV2 Regs;
    bool Captured = false;
};

// ValidLevels excludes the invalid terminator.  Sentinel is exactly the first
// invalid subleaf and no later reserved subleaf is queried or retained.
struct CpuidTopologyEnumerationV2
{
    uint32_t Leaf = 0u;
    std::vector<CpuidRecordV2> ValidLevels;
    CpuidRecordV2 Sentinel;
    bool SentinelCaptured = false;
};

struct CpuidSnapshotV2
{
    CpuidRecordV2 Leaf0;
    CpuidRecordV2 Leaf1;
    CpuidRecordV2 Leaf6;
    CpuidRecordV2 Leaf80000000;
    CpuidRecordV2 Leaf80000001;
    CpuidRecordV2 Leaf80000008;
    CpuidRecordV2 Leaf8000001e;
    CpuidRecordV2 Leaf80000021;
    CpuidTopologyEnumerationV2 LeafB;
    CpuidTopologyEnumerationV2 Leaf80000026;
};

struct TargetContextReceiptV2
{
    int32_t Cpu = -1;
    std::vector<int32_t> Affinity;
    int64_t VoluntaryContextSwitches = -1;
    int64_t InvoluntaryContextSwitches = -1;
    bool Captured = false;
};

struct DerivedTargetIdentityV2
{
    uint32_t Family = 0u;
    uint32_t Model = 0u;
    uint32_t Stepping = 0u;
    uint32_t InitialApicId8 = 0u;
    uint32_t FullApicId = 0u;
    uint32_t ThreadId = 0u;
    uint32_t CoreId = 0u;
    uint32_t ComplexId = 0u;
    uint32_t CcdId = 0u;
    uint32_t PackageId = 0u;
    uint32_t ThreadsPerCore = 0u;
    uint32_t LogicalProcessorsPerPackage = 0u;
};

struct TargetIdentityReceiptV2
{
    int32_t RequestedCpu = -1;
    TargetContextReceiptV2 Before;
    TargetContextReceiptV2 After;
    CpuidSnapshotV2 Raw;
    DerivedTargetIdentityV2 Derived;
    std::string CanonicalSha256;
    bool RawCaptureComplete = false;
    bool SemanticValidationPassed = false;
};

typedef bool (*CpuidReadFunctionV2)(
    void* context,
    uint32_t leaf,
    uint32_t subleaf,
    CpuidRegsV2& registers,
    std::string& diagnostic);

typedef bool (*TargetContextCaptureFunctionV2)(
    void* context,
    int32_t target_cpu,
    TargetContextReceiptV2& receipt,
    std::string& diagnostic);

// Capture helpers are public so isolated synthetic tests can prove the exact
// enumeration surface without executing CPUID or touching affinity.
bool EnumerateCpuidLeafBV2(
    CpuidReadFunctionV2 reader,
    void* reader_context,
    CpuidTopologyEnumerationV2& enumeration,
    std::string& diagnostic);

bool EnumerateCpuidLeaf80000026V2(
    CpuidReadFunctionV2 reader,
    void* reader_context,
    CpuidTopologyEnumerationV2& enumeration,
    std::string& diagnostic);

bool CaptureCpuidSnapshotV2(
    CpuidReadFunctionV2 reader,
    void* reader_context,
    CpuidSnapshotV2& snapshot,
    std::string& diagnostic);

bool ValidateTargetContextPairV2(
    int32_t target_cpu,
    const TargetContextReceiptV2& before,
    const TargetContextReceiptV2& after,
    std::string& diagnostic);

bool ValidateTargetIdentitySemanticsV2(
    const TargetIdentityReceiptV2& receipt,
    DerivedTargetIdentityV2& derived,
    std::string& diagnostic);

bool ValidateFrozenCpu50IdentityV2(
    const TargetIdentityReceiptV2& receipt,
    std::string& diagnostic);

bool CaptureTargetIdentityV2(
    int32_t target_cpu,
    CpuidReadFunctionV2 cpuid_reader,
    void* cpuid_context,
    TargetContextCaptureFunctionV2 context_reader,
    void* context,
    TargetIdentityReceiptV2& receipt,
    std::string& diagnostic);

// This does not set affinity.  The caller must pin CPU50 once; every capture
// then verifies the effective singleton mask, sched_getcpu(), and unchanged
// thread context-switch counters.  A true return proves a complete semantic
// receipt, not equality to the frozen host.  Every production caller MUST next
// call ValidateFrozenCpu50IdentityV2 before accepting the checkpoint.
// Unsupported platforms fail deterministically.
bool CaptureNativeTargetIdentityV2(
    TargetIdentityReceiptV2& receipt,
    std::string& diagnostic);

CpuidSnapshotV2 FrozenCpu50CpuidSnapshotV2();
DerivedTargetIdentityV2 FrozenCpu50DerivedIdentityV2();
TargetIdentityReceiptV2 FrozenCpu50TargetIdentityReceiptV2();

bool SerializeTargetIdentityV2(
    const TargetIdentityReceiptV2& receipt,
    std::string& canonical_bytes,
    std::string& diagnostic);

std::string TargetIdentitySha256V2(
    const TargetIdentityReceiptV2& receipt,
    std::string& diagnostic);

std::string FrozenCpu50TargetIdentitySha256V2();
uint64_t FrozenCpu50TargetIdentityCanonicalBytesV2();

std::string FormatTargetIdentityEvidenceV2(
    const std::string& checkpoint,
    const TargetIdentityReceiptV2& receipt,
    bool gate_passed,
    const std::string& diagnostic);

struct CpuInfoMapRecordV2
{
    uint32_t Cpu = 0u;
    uint32_t Package = 0u;
    uint32_t Core = 0u;
    uint32_t ApicId = 0u;
    uint32_t InitialApicId = 0u;
};

bool ParseProcCpuInfoV2(
    const std::string& content,
    std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic);

bool CanonicalizeCpuInfoMapV2(
    const std::vector<CpuInfoMapRecordV2>& records,
    std::string& canonical,
    std::string& diagnostic);

std::string CpuInfoMapSha256V2(
    const std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic);

bool ReadProcCpuInfoV2(
    // Production identity capture accepts only this exact path:
    // /proc/cpuinfo.  Tests inject text through ParseProcCpuInfoV2 instead.
    const std::string& path,
    std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic);

std::vector<CpuInfoMapRecordV2> FrozenCpuInfoMapV2();

bool ValidateFrozenCpuInfoMapV2(
    const std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic);

const char* FrozenCpuInfoMapSha256V2();
uint64_t FrozenCpuInfoMapCanonicalBytesV2();

// Runs only deterministic, injected logic.  It performs no native CPUID,
// affinity, codec, counter, timing, or calibration operation.
bool RunTargetIdentityV2SelfTests(std::string& diagnostic);

} // namespace wirehair_wh2_bench
