#include "Wh2RdpruTargetIdentityV2.h"

#include "Wh2FrozenTrace.h"

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#if defined(__linux__)
#include <fcntl.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/vfs.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#if defined(__linux__) && defined(__x86_64__) && \
    (defined(__GNUC__) || defined(__clang__))
#include <cpuid.h>
#define WIREHAIR_RDPRU_TARGET_IDENTITY_V2_NATIVE 1
#else
#define WIREHAIR_RDPRU_TARGET_IDENTITY_V2_NATIVE 0
#endif

namespace wirehair_wh2_bench {
namespace {

static const char kIdentityDomain[] =
    "wirehair-rdpru-target-identity-v2";
static const uint32_t kIdentityVersion = 2u;
static const int32_t kFrozenTargetCpu = 50;
static const uint32_t kFrozenFullApicId = UINT32_C(0x64);
static const uint32_t kMaximumTopologyLevels = 32u;
static const uint64_t kFrozenIdentityCanonicalBytes = UINT64_C(617);
static const char kFrozenIdentitySha256[] =
    "3288e0ef61cf3e628dcd827f9cf003c9d6ec6b5a12169e7a8bfc796baacddba7";
static const std::size_t kMaximumCpuInfoBytes = 1024u * 1024u;
static const std::size_t kMaximumCpuInfoLineBytes = 4096u;
static const std::size_t kFrozenCpuInfoRecords = 128u;
static const uint64_t kFrozenCpuInfoCanonicalBytes = UINT64_C(6818);
static const char kFrozenCpuInfoSha256[] =
    "d51bf8bf7c1b249d081d8447fd35f118f2d5bfeee49d2a8f3bc2aba507744d8b";

bool RegsEqual(const CpuidRegsV2& a, const CpuidRegsV2& b)
{
    return a.Eax == b.Eax && a.Ebx == b.Ebx &&
        a.Ecx == b.Ecx && a.Edx == b.Edx;
}

bool RecordEqual(const CpuidRecordV2& a, const CpuidRecordV2& b)
{
    return a.Leaf == b.Leaf && a.Subleaf == b.Subleaf &&
        a.Captured == b.Captured && RegsEqual(a.Regs, b.Regs);
}

bool TopologyEqual(
    const CpuidTopologyEnumerationV2& a,
    const CpuidTopologyEnumerationV2& b)
{
    if (a.Leaf != b.Leaf ||
        a.SentinelCaptured != b.SentinelCaptured ||
        a.ValidLevels.size() != b.ValidLevels.size() ||
        !RecordEqual(a.Sentinel, b.Sentinel)) return false;
    for (std::size_t i = 0u; i < a.ValidLevels.size(); ++i) {
        if (!RecordEqual(a.ValidLevels[i], b.ValidLevels[i])) return false;
    }
    return true;
}

bool SnapshotEqual(const CpuidSnapshotV2& a, const CpuidSnapshotV2& b)
{
    return RecordEqual(a.Leaf0, b.Leaf0) &&
        RecordEqual(a.Leaf1, b.Leaf1) &&
        RecordEqual(a.Leaf6, b.Leaf6) &&
        RecordEqual(a.Leaf80000000, b.Leaf80000000) &&
        RecordEqual(a.Leaf80000001, b.Leaf80000001) &&
        RecordEqual(a.Leaf80000008, b.Leaf80000008) &&
        RecordEqual(a.Leaf8000001e, b.Leaf8000001e) &&
        RecordEqual(a.Leaf80000021, b.Leaf80000021) &&
        TopologyEqual(a.LeafB, b.LeafB) &&
        TopologyEqual(a.Leaf80000026, b.Leaf80000026);
}

bool DerivedEqual(
    const DerivedTargetIdentityV2& a,
    const DerivedTargetIdentityV2& b)
{
    return a.Family == b.Family && a.Model == b.Model &&
        a.Stepping == b.Stepping &&
        a.InitialApicId8 == b.InitialApicId8 &&
        a.FullApicId == b.FullApicId &&
        a.ThreadId == b.ThreadId && a.CoreId == b.CoreId &&
        a.ComplexId == b.ComplexId && a.CcdId == b.CcdId &&
        a.PackageId == b.PackageId &&
        a.ThreadsPerCore == b.ThreadsPerCore &&
        a.LogicalProcessorsPerPackage == b.LogicalProcessorsPerPackage;
}

void AppendLe32(std::string& bytes, uint32_t value)
{
    for (unsigned shift = 0u; shift != 32u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendLe64(std::string& bytes, uint64_t value)
{
    for (unsigned shift = 0u; shift != 64u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendSizedString(std::string& bytes, const std::string& value)
{
    AppendLe64(bytes, static_cast<uint64_t>(value.size()));
    bytes.append(value);
}

void AppendRecord(std::string& bytes, const CpuidRecordV2& record)
{
    AppendLe32(bytes, record.Captured ? 1u : 0u);
    AppendLe32(bytes, record.Leaf);
    AppendLe32(bytes, record.Subleaf);
    AppendLe32(bytes, record.Regs.Eax);
    AppendLe32(bytes, record.Regs.Ebx);
    AppendLe32(bytes, record.Regs.Ecx);
    AppendLe32(bytes, record.Regs.Edx);
}

void AppendTopology(
    std::string& bytes,
    const CpuidTopologyEnumerationV2& topology)
{
    AppendLe32(bytes, topology.Leaf);
    AppendLe32(bytes, static_cast<uint32_t>(topology.ValidLevels.size()));
    for (std::size_t i = 0u; i < topology.ValidLevels.size(); ++i) {
        AppendRecord(bytes, topology.ValidLevels[i]);
    }
    AppendLe32(bytes, topology.SentinelCaptured ? 1u : 0u);
    AppendRecord(bytes, topology.Sentinel);
}

std::string HexU32(uint32_t value)
{
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::hex << std::nouppercase << std::setw(8)
         << std::setfill('0') << value;
    return text.str();
}

std::string HexEncode(const std::string& value)
{
    static const char alphabet[] = "0123456789abcdef";
    std::string encoded;
    encoded.reserve(value.size() * 2u);
    for (std::size_t i = 0u; i < value.size(); ++i)
    {
        const unsigned byte = static_cast<unsigned char>(value[i]);
        encoded.push_back(alphabet[byte >> 4u]);
        encoded.push_back(alphabet[byte & 15u]);
    }
    return encoded;
}

bool IsLowerHexSha256Value(const std::string& value)
{
    if (value.size() != 64u) return false;
    for (std::size_t i = 0u; i < value.size(); ++i) {
        if (!((value[i] >= '0' && value[i] <= '9') ||
              (value[i] >= 'a' && value[i] <= 'f'))) return false;
    }
    return true;
}

std::string FormatAffinity(const std::vector<int32_t>& affinity)
{
    std::ostringstream text;
    text.imbue(std::locale::classic());
    for (std::size_t i = 0u; i < affinity.size(); ++i)
    {
        if (i != 0u) text << ',';
        text << affinity[i];
    }
    return text.str();
}

void DecodeDisplayFamilyModel(
    uint32_t eax,
    uint32_t& family,
    uint32_t& model,
    uint32_t& stepping)
{
    stepping = eax & 0xfu;
    const uint32_t base_family = (eax >> 8u) & 0xfu;
    const uint32_t base_model = (eax >> 4u) & 0xfu;
    family = base_family == 0xfu ?
        base_family + ((eax >> 20u) & 0xffu) : base_family;
    model = base_model;
    if (base_family == 0x6u || base_family == 0xfu) {
        model += ((eax >> 16u) & 0xfu) << 4u;
    }
}

bool IsAuthenticAmd(const CpuidRegsV2& regs)
{
    return regs.Ebx == UINT32_C(0x68747541) &&
        regs.Edx == UINT32_C(0x69746e65) &&
        regs.Ecx == UINT32_C(0x444d4163);
}

uint32_t TopologyCount(const CpuidRegsV2& regs)
{
    return regs.Ebx & 0xffffu;
}

uint32_t TopologyType(const CpuidRegsV2& regs)
{
    return (regs.Ecx >> 8u) & 0xffu;
}

uint32_t TopologyLevel(const CpuidRegsV2& regs)
{
    return regs.Ecx & 0xffu;
}

uint32_t TopologyShift(const CpuidRegsV2& regs)
{
    return regs.Eax & 0x1fu;
}

bool CaptureRecord(
    CpuidReadFunctionV2 reader,
    void* reader_context,
    uint32_t leaf,
    uint32_t subleaf,
    CpuidRecordV2& record,
    std::string& diagnostic)
{
    if (!reader)
    {
        diagnostic = "CPUID reader is null";
        return false;
    }
    record = CpuidRecordV2();
    record.Leaf = leaf;
    record.Subleaf = subleaf;
    CpuidRegsV2 regs;
    if (!reader(reader_context, leaf, subleaf, regs, diagnostic)) return false;
    record.Regs = regs;
    record.Captured = true;
    return true;
}

bool EnumerateTopology(
    uint32_t leaf,
    CpuidReadFunctionV2 reader,
    void* reader_context,
    CpuidTopologyEnumerationV2& enumeration,
    std::string& diagnostic)
{
    enumeration = CpuidTopologyEnumerationV2();
    enumeration.Leaf = leaf;
    for (uint32_t subleaf = 0u; subleaf < kMaximumTopologyLevels; ++subleaf)
    {
        CpuidRecordV2 record;
        if (!CaptureRecord(
                reader, reader_context, leaf, subleaf, record, diagnostic))
        {
            return false;
        }
        if (TopologyCount(record.Regs) == 0u ||
            TopologyType(record.Regs) == 0u)
        {
            enumeration.Sentinel = record;
            enumeration.SentinelCaptured = true;
            return true;
        }
        enumeration.ValidLevels.push_back(record);
    }
    diagnostic = "CPUID topology enumeration has no bounded sentinel";
    return false;
}

bool ValidateRecordTag(
    const CpuidRecordV2& record,
    uint32_t leaf,
    uint32_t subleaf,
    std::string& diagnostic)
{
    if (!record.Captured || record.Leaf != leaf || record.Subleaf != subleaf)
    {
        diagnostic = "CPUID record is missing or has the wrong leaf tag";
        return false;
    }
    return true;
}

bool ValidateTopologyShape(
    const CpuidTopologyEnumerationV2& topology,
    const uint32_t* expected_types,
    std::size_t expected_levels,
    uint32_t& full_apic_id,
    std::string& diagnostic)
{
    if (topology.ValidLevels.size() != expected_levels ||
        !topology.SentinelCaptured ||
        topology.Sentinel.Leaf != topology.Leaf ||
        topology.Sentinel.Subleaf != expected_levels ||
        !topology.Sentinel.Captured)
    {
        diagnostic = "CPUID topology level count or sentinel mismatch";
        return false;
    }
    const CpuidRegsV2& sentinel = topology.Sentinel.Regs;
    if (TopologyCount(sentinel) != 0u ||
        TopologyType(sentinel) != 0u ||
        TopologyLevel(sentinel) != expected_levels)
    {
        diagnostic = "CPUID topology sentinel is invalid";
        return false;
    }
    if (topology.Leaf == 0xbu && sentinel.Ebx != 0u)
    {
        diagnostic = "CPUID leaf B sentinel has nonzero reserved EBX bits";
        return false;
    }
    uint32_t previous_shift = 0u;
    uint32_t previous_count = 0u;
    for (std::size_t i = 0u; i < topology.ValidLevels.size(); ++i)
    {
        const CpuidRecordV2& record = topology.ValidLevels[i];
        const CpuidRegsV2& regs = record.Regs;
        if (!record.Captured || record.Leaf != topology.Leaf ||
            record.Subleaf != i || TopologyLevel(regs) != i)
        {
            diagnostic = "CPUID topology subleaf numbering mismatch";
            return false;
        }
        const uint32_t count = TopologyCount(regs);
        const uint32_t type = TopologyType(regs);
        const uint32_t shift = TopologyShift(regs);
        if (count == 0u || type != expected_types[i] ||
            (topology.Leaf == 0xbu && (regs.Ebx >> 16u) != 0u) ||
            (regs.Eax & ~UINT32_C(0x1f)) != 0u ||
            (regs.Ecx & UINT32_C(0xffff0000)) != 0u ||
            (i != 0u && (shift < previous_shift || count < previous_count)) ||
            static_cast<uint64_t>(count) > (UINT64_C(1) << shift))
        {
            diagnostic = "CPUID topology type/count/shift fields are invalid";
            return false;
        }
        if (i == 0u) full_apic_id = regs.Edx;
        else if (regs.Edx != full_apic_id)
        {
            diagnostic = "CPUID topology APIC ID changes across valid levels";
            return false;
        }
        previous_shift = shift;
        previous_count = count;
    }
    if (topology.Sentinel.Regs.Edx != full_apic_id)
    {
        diagnostic = "CPUID topology APIC ID changes at the sentinel";
        return false;
    }
    return true;
}

bool ContextEqualNormalized(
    const TargetContextReceiptV2& context,
    int32_t target_cpu)
{
    return context.Captured && context.Cpu == target_cpu &&
        context.Affinity.size() == 1u &&
        context.Affinity[0] == target_cpu;
}

void AppendContextNormalized(
    std::string& bytes,
    const TargetContextReceiptV2& context)
{
    AppendLe32(bytes, context.Captured ? 1u : 0u);
    AppendLe32(bytes, static_cast<uint32_t>(context.Cpu));
    AppendLe32(bytes, static_cast<uint32_t>(context.Affinity.size()));
    for (std::size_t i = 0u; i < context.Affinity.size(); ++i) {
        AppendLe32(bytes, static_cast<uint32_t>(context.Affinity[i]));
    }
}

void AppendEvidenceRecord(
    std::ostringstream& evidence,
    const std::string& checkpoint,
    const CpuidRecordV2& record,
    const char* role)
{
    if (!record.Captured) return;
    evidence << "TARGET_CPUID_V2,checkpoint_length=" << checkpoint.size()
             << ",checkpoint_hex=" << HexEncode(checkpoint)
             << ",role=" << role
             << ",leaf=0x" << HexU32(record.Leaf)
             << ",subleaf=0x" << HexU32(record.Subleaf)
             << ",eax=0x" << HexU32(record.Regs.Eax)
             << ",ebx=0x" << HexU32(record.Regs.Ebx)
             << ",ecx=0x" << HexU32(record.Regs.Ecx)
             << ",edx=0x" << HexU32(record.Regs.Edx) << '\n';
}

std::string Trim(const std::string& text)
{
    std::size_t first = 0u;
    while (first < text.size() &&
           (text[first] == ' ' || text[first] == '\t' || text[first] == '\r')) {
        ++first;
    }
    std::size_t last = text.size();
    while (last > first &&
           (text[last - 1u] == ' ' || text[last - 1u] == '\t' ||
            text[last - 1u] == '\r')) {
        --last;
    }
    return text.substr(first, last - first);
}

bool ParseDecimalU32(const std::string& input, uint32_t& value)
{
    const std::string text = Trim(input);
    if (text.empty()) return false;
    uint64_t parsed = 0u;
    for (std::size_t i = 0u; i < text.size(); ++i)
    {
        const char c = text[i];
        if (c < '0' || c > '9') return false;
        parsed = parsed * 10u + static_cast<unsigned>(c - '0');
        if (parsed > std::numeric_limits<uint32_t>::max()) return false;
    }
    value = static_cast<uint32_t>(parsed);
    return true;
}

bool CpuInfoRecordEqual(
    const CpuInfoMapRecordV2& a,
    const CpuInfoMapRecordV2& b)
{
    return a.Cpu == b.Cpu && a.Package == b.Package &&
        a.Core == b.Core && a.ApicId == b.ApicId &&
        a.InitialApicId == b.InitialApicId;
}

#if WIREHAIR_RDPRU_TARGET_IDENTITY_V2_NATIVE
void NativeCompilerBarrier()
{
    __asm__ __volatile__("" ::: "memory");
}

bool NativeCpuidRead(
    void*,
    uint32_t leaf,
    uint32_t subleaf,
    CpuidRegsV2& registers,
    std::string& diagnostic)
{
    unsigned eax = 0u;
    unsigned ebx = 0u;
    unsigned ecx = 0u;
    unsigned edx = 0u;
    NativeCompilerBarrier();
    __cpuid_count(leaf, subleaf, eax, ebx, ecx, edx);
    NativeCompilerBarrier();
    registers.Eax = eax;
    registers.Ebx = ebx;
    registers.Ecx = ecx;
    registers.Edx = edx;
    diagnostic.clear();
    return true;
}

bool NativeContextCapture(
    void*,
    int32_t target_cpu,
    TargetContextReceiptV2& receipt,
    std::string& diagnostic)
{
    receipt = TargetContextReceiptV2();
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "target CPU is outside cpu_set_t range";
        return false;
    }
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    if (sched_getaffinity(0, sizeof(affinity), &affinity) != 0)
    {
        diagnostic = std::string("sched_getaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        if (CPU_ISSET(cpu, &affinity)) receipt.Affinity.push_back(cpu);
    }
    const int cpu = sched_getcpu();
    if (cpu < 0)
    {
        diagnostic = std::string("sched_getcpu failed: ") +
            std::strerror(errno);
        return false;
    }
    struct rusage usage;
    if (getrusage(RUSAGE_THREAD, &usage) != 0)
    {
        diagnostic = std::string("getrusage(RUSAGE_THREAD) failed: ") +
            std::strerror(errno);
        return false;
    }
    receipt.Cpu = cpu;
    receipt.VoluntaryContextSwitches =
        static_cast<int64_t>(usage.ru_nvcsw);
    receipt.InvoluntaryContextSwitches =
        static_cast<int64_t>(usage.ru_nivcsw);
    receipt.Captured = true;
    diagnostic.clear();
    return true;
}
#endif

} // namespace

bool EnumerateCpuidLeafBV2(
    CpuidReadFunctionV2 reader,
    void* reader_context,
    CpuidTopologyEnumerationV2& enumeration,
    std::string& diagnostic)
{
    return EnumerateTopology(
        UINT32_C(0x0000000b), reader, reader_context, enumeration, diagnostic);
}

bool EnumerateCpuidLeaf80000026V2(
    CpuidReadFunctionV2 reader,
    void* reader_context,
    CpuidTopologyEnumerationV2& enumeration,
    std::string& diagnostic)
{
    return EnumerateTopology(
        UINT32_C(0x80000026), reader, reader_context, enumeration, diagnostic);
}

bool CaptureCpuidSnapshotV2(
    CpuidReadFunctionV2 reader,
    void* reader_context,
    CpuidSnapshotV2& snapshot,
    std::string& diagnostic)
{
    snapshot = CpuidSnapshotV2();
    if (!CaptureRecord(reader, reader_context, 0u, 0u, snapshot.Leaf0, diagnostic) ||
        !CaptureRecord(
            reader,
            reader_context,
            UINT32_C(0x80000000),
            0u,
            snapshot.Leaf80000000,
            diagnostic)) return false;
    if (!IsAuthenticAmd(snapshot.Leaf0.Regs) ||
        !IsAuthenticAmd(snapshot.Leaf80000000.Regs))
    {
        diagnostic = "CPUID vendor is not consistently AuthenticAMD";
        return false;
    }
    const uint32_t maximum_basic = snapshot.Leaf0.Regs.Eax;
    const uint32_t maximum_extended = snapshot.Leaf80000000.Regs.Eax;
    if (maximum_basic < 1u || maximum_extended < UINT32_C(0x80000001))
    {
        diagnostic = "required CPUID identity leaves are unavailable";
        return false;
    }
    if (!CaptureRecord(reader, reader_context, 1u, 0u, snapshot.Leaf1, diagnostic)) {
        return false;
    }
    if (maximum_basic < 6u)
    {
        diagnostic = "CPUID leaf 6 is unavailable";
        return false;
    }
    if (!CaptureRecord(reader, reader_context, 6u, 0u, snapshot.Leaf6, diagnostic) ||
        !CaptureRecord(
            reader,
            reader_context,
            UINT32_C(0x80000001),
            0u,
            snapshot.Leaf80000001,
            diagnostic)) return false;
    if (maximum_extended < UINT32_C(0x80000008))
    {
        diagnostic = "CPUID leaf 80000008 is unavailable";
        return false;
    }
    if (!CaptureRecord(
            reader,
            reader_context,
            UINT32_C(0x80000008),
            0u,
            snapshot.Leaf80000008,
            diagnostic)) return false;
    if (maximum_extended < UINT32_C(0x8000001e) ||
        (snapshot.Leaf80000001.Regs.Ecx & (UINT32_C(1) << 22u)) == 0u)
    {
        diagnostic = "CPUID 8000001E is not enabled by TopologyExtensions";
        return false;
    }
    if (!CaptureRecord(
            reader,
            reader_context,
            UINT32_C(0x8000001e),
            0u,
            snapshot.Leaf8000001e,
            diagnostic)) return false;
    if (maximum_extended < UINT32_C(0x80000021))
    {
        diagnostic = "CPUID leaf 80000021 is unavailable";
        return false;
    }
    if (!CaptureRecord(
            reader,
            reader_context,
            UINT32_C(0x80000021),
            0u,
            snapshot.Leaf80000021,
            diagnostic)) return false;
    if (maximum_basic < UINT32_C(0x0b))
    {
        diagnostic = "CPUID leaf B is unavailable";
        return false;
    }
    if (!EnumerateCpuidLeafBV2(
            reader, reader_context, snapshot.LeafB, diagnostic)) return false;
    if (maximum_extended < UINT32_C(0x80000026))
    {
        diagnostic = "CPUID leaf 80000026 is unavailable";
        return false;
    }
    if (!EnumerateCpuidLeaf80000026V2(
            reader, reader_context, snapshot.Leaf80000026, diagnostic)) {
        return false;
    }
    diagnostic.clear();
    return true;
}

bool ValidateTargetContextPairV2(
    int32_t target_cpu,
    const TargetContextReceiptV2& before,
    const TargetContextReceiptV2& after,
    std::string& diagnostic)
{
    if (target_cpu < 0 || !ContextEqualNormalized(before, target_cpu) ||
        !ContextEqualNormalized(after, target_cpu))
    {
        diagnostic = "CPU migration or effective affinity change in CPUID bracket";
        return false;
    }
    if (before.VoluntaryContextSwitches < 0 ||
        before.InvoluntaryContextSwitches < 0 ||
        after.VoluntaryContextSwitches < before.VoluntaryContextSwitches ||
        after.InvoluntaryContextSwitches < before.InvoluntaryContextSwitches)
    {
        diagnostic = "invalid target context-switch receipt";
        return false;
    }
    if (after.VoluntaryContextSwitches != before.VoluntaryContextSwitches ||
        after.InvoluntaryContextSwitches != before.InvoluntaryContextSwitches)
    {
        diagnostic = "context switch occurred in CPUID bracket";
        return false;
    }
    diagnostic.clear();
    return true;
}

bool ValidateTargetIdentitySemanticsV2(
    const TargetIdentityReceiptV2& receipt,
    DerivedTargetIdentityV2& derived,
    std::string& diagnostic)
{
    derived = DerivedTargetIdentityV2();
    if (!receipt.RawCaptureComplete)
    {
        diagnostic = "CPUID snapshot is incomplete";
        return false;
    }
    if (!ValidateTargetContextPairV2(
            receipt.RequestedCpu, receipt.Before, receipt.After, diagnostic)) {
        return false;
    }
    const CpuidSnapshotV2& raw = receipt.Raw;
    const struct RequiredRecord
    {
        const CpuidRecordV2* Record;
        uint32_t Leaf;
    } required[] = {
        { &raw.Leaf0, 0u }, { &raw.Leaf1, 1u }, { &raw.Leaf6, 6u },
        { &raw.Leaf80000000, UINT32_C(0x80000000) },
        { &raw.Leaf80000001, UINT32_C(0x80000001) },
        { &raw.Leaf80000008, UINT32_C(0x80000008) },
        { &raw.Leaf8000001e, UINT32_C(0x8000001e) },
        { &raw.Leaf80000021, UINT32_C(0x80000021) }
    };
    for (std::size_t i = 0u; i < sizeof(required) / sizeof(required[0]); ++i) {
        if (!ValidateRecordTag(
                *required[i].Record, required[i].Leaf, 0u, diagnostic)) {
            return false;
        }
    }
    if (!IsAuthenticAmd(raw.Leaf0.Regs) ||
        !IsAuthenticAmd(raw.Leaf80000000.Regs) ||
        raw.Leaf0.Regs.Eax < UINT32_C(0x0b) ||
        raw.Leaf80000000.Regs.Eax < UINT32_C(0x80000026))
    {
        diagnostic = "CPUID vendor or maximum-leaf gate failed";
        return false;
    }
    if ((raw.Leaf1.Regs.Ecx & (UINT32_C(1) << 31u)) != 0u)
    {
        diagnostic = "CPUID reports a hypervisor";
        return false;
    }
    if ((raw.Leaf6.Regs.Ecx & 1u) == 0u ||
        (raw.Leaf80000001.Regs.Ecx & (UINT32_C(1) << 22u)) == 0u ||
        (raw.Leaf80000008.Regs.Ebx & (UINT32_C(1) << 4u)) == 0u ||
        (raw.Leaf80000008.Regs.Edx >> 16u) < 1u ||
        (raw.Leaf80000021.Regs.Eax & (UINT32_C(1) << 2u)) == 0u)
    {
        diagnostic = "required APERF/MPERF, topology, RDPRU, or LFENCE feature is missing";
        return false;
    }
    static const uint32_t b_types[] = { 1u, 2u };
    static const uint32_t x26_types[] = { 1u, 2u, 3u, 4u };
    uint32_t b_apic = 0u;
    uint32_t x26_apic = 0u;
    if (!ValidateTopologyShape(
            raw.LeafB, b_types, 2u, b_apic, diagnostic) ||
        !ValidateTopologyShape(
            raw.Leaf80000026, x26_types, 4u, x26_apic, diagnostic)) {
        return false;
    }
    if (raw.LeafB.Leaf != UINT32_C(0x0b) ||
        raw.Leaf80000026.Leaf != UINT32_C(0x80000026) ||
        b_apic != x26_apic || raw.Leaf8000001e.Regs.Eax != x26_apic ||
        ((raw.Leaf1.Regs.Ebx >> 24u) & 0xffu) != (x26_apic & 0xffu))
    {
        diagnostic = "CPUID APIC identities disagree across available leaves";
        return false;
    }
    const CpuidRegsV2& b_core = raw.LeafB.ValidLevels[0].Regs;
    const CpuidRegsV2& b_package = raw.LeafB.ValidLevels[1].Regs;
    const CpuidRegsV2& x26_core = raw.Leaf80000026.ValidLevels[0].Regs;
    const CpuidRegsV2& x26_complex = raw.Leaf80000026.ValidLevels[1].Regs;
    const CpuidRegsV2& x26_ccd = raw.Leaf80000026.ValidLevels[2].Regs;
    const CpuidRegsV2& x26_package = raw.Leaf80000026.ValidLevels[3].Regs;
    if (TopologyShift(b_core) != TopologyShift(x26_core) ||
        TopologyCount(b_core) != TopologyCount(x26_core) ||
        TopologyShift(b_package) != TopologyShift(x26_package) ||
        TopologyCount(b_package) != TopologyCount(x26_package))
    {
        diagnostic = "CPUID B and 80000026 topology shapes disagree";
        return false;
    }
    const uint32_t package_shift = TopologyShift(x26_package);
    const uint32_t core_shift = TopologyShift(x26_core);
    if (package_shift < core_shift || package_shift > 31u)
    {
        diagnostic = "CPUID package/core shifts are invalid";
        return false;
    }
    const uint64_t package_mask = package_shift == 0u ? 0u :
        ((UINT64_C(1) << package_shift) - 1u);
    const uint32_t package_local_apic = static_cast<uint32_t>(
        static_cast<uint64_t>(x26_apic) & package_mask);
    derived.FullApicId = x26_apic;
    derived.InitialApicId8 = (raw.Leaf1.Regs.Ebx >> 24u) & 0xffu;
    derived.ThreadId = core_shift == 0u ? 0u :
        x26_apic & static_cast<uint32_t>((UINT64_C(1) << core_shift) - 1u);
    derived.CoreId = package_local_apic >> core_shift;
    derived.ComplexId = package_local_apic >> TopologyShift(x26_complex);
    derived.CcdId = package_local_apic >> TopologyShift(x26_ccd);
    derived.PackageId = x26_apic >> package_shift;
    derived.ThreadsPerCore =
        ((raw.Leaf8000001e.Regs.Ebx >> 8u) & 0xffu) + 1u;
    derived.LogicalProcessorsPerPackage =
        TopologyCount(x26_package);
    DecodeDisplayFamilyModel(
        raw.Leaf1.Regs.Eax, derived.Family, derived.Model, derived.Stepping);
    if ((raw.Leaf8000001e.Regs.Ebx & 0xffu) != derived.CoreId ||
        derived.ThreadsPerCore != TopologyCount(x26_core) ||
        ((raw.Leaf1.Regs.Ebx >> 16u) & 0xffu) !=
            derived.LogicalProcessorsPerPackage ||
        ((raw.Leaf80000008.Regs.Ecx & 0xffu) + 1u) !=
            derived.LogicalProcessorsPerPackage ||
        ((raw.Leaf80000008.Regs.Ecx >> 12u) & 0xfu) != package_shift)
    {
        diagnostic = "CPUID derived core/thread/package topology disagrees";
        return false;
    }
    diagnostic.clear();
    return true;
}

bool ValidateFrozenCpu50IdentityV2(
    const TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
    DerivedTargetIdentityV2 derived;
    if (!ValidateTargetIdentitySemanticsV2(receipt, derived, diagnostic)) {
        return false;
    }
    if (receipt.RequestedCpu != kFrozenTargetCpu ||
        !receipt.SemanticValidationPassed ||
        !SnapshotEqual(receipt.Raw, FrozenCpu50CpuidSnapshotV2()) ||
        !DerivedEqual(receipt.Derived, derived) ||
        !DerivedEqual(derived, FrozenCpu50DerivedIdentityV2()))
    {
        diagnostic = "target identity differs from frozen Linux CPU50/APIC64";
        return false;
    }
    std::string hash_diagnostic;
    const std::string recomputed =
        TargetIdentitySha256V2(receipt, hash_diagnostic);
    if (recomputed.empty() || receipt.CanonicalSha256 != recomputed ||
        recomputed != kFrozenIdentitySha256)
    {
        diagnostic = "stored target identity SHA-256 disagrees with raw receipt";
        return false;
    }
    diagnostic.clear();
    return true;
}

bool CaptureTargetIdentityV2(
    int32_t target_cpu,
    CpuidReadFunctionV2 cpuid_reader,
    void* cpuid_context,
    TargetContextCaptureFunctionV2 context_reader,
    void* context,
    TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
    receipt = TargetIdentityReceiptV2();
    receipt.RequestedCpu = target_cpu;
    if (!cpuid_reader || !context_reader)
    {
        diagnostic = "target identity capture callback is null";
        return false;
    }
    try {
        if (!context_reader(context, target_cpu, receipt.Before, diagnostic)) {
            if (diagnostic.empty()) diagnostic = "before-context capture failed";
            return false;
        }
    }
    catch (const std::exception& error) {
        diagnostic = std::string("before-context reader threw: ") + error.what();
        return false;
    }
    catch (...) {
        diagnostic = "before-context reader threw a non-standard exception";
        return false;
    }
    if (!ContextEqualNormalized(receipt.Before, target_cpu) ||
        receipt.Before.VoluntaryContextSwitches < 0 ||
        receipt.Before.InvoluntaryContextSwitches < 0)
    {
        diagnostic = "before-context is not the requested singleton CPU";
        return false;
    }
    bool raw_ok = false;
    std::string raw_diagnostic;
    try {
        raw_ok = CaptureCpuidSnapshotV2(
            cpuid_reader, cpuid_context, receipt.Raw, raw_diagnostic);
    }
    catch (const std::exception& error) {
        raw_diagnostic = std::string("CPUID reader threw: ") + error.what();
    }
    catch (...) {
        raw_diagnostic = "CPUID reader threw a non-standard exception";
    }
    receipt.RawCaptureComplete = raw_ok;
    std::string after_diagnostic;
    bool after_ok = false;
    try {
        after_ok = context_reader(
            context, target_cpu, receipt.After, after_diagnostic);
    }
    catch (const std::exception& error) {
        after_diagnostic = std::string("after-context reader threw: ") +
            error.what();
    }
    catch (...) {
        after_diagnostic =
            "after-context reader threw a non-standard exception";
    }
    if (!after_ok)
    {
        diagnostic = after_diagnostic.empty() ?
            "after-context capture failed" : after_diagnostic;
        return false;
    }
    if (!ValidateTargetContextPairV2(
            target_cpu, receipt.Before, receipt.After, diagnostic)) {
        return false;
    }
    if (!raw_ok)
    {
        diagnostic = raw_diagnostic.empty() ?
            "CPUID snapshot capture failed" : raw_diagnostic;
        return false;
    }
    DerivedTargetIdentityV2 derived;
    if (!ValidateTargetIdentitySemanticsV2(receipt, derived, diagnostic)) {
        return false;
    }
    receipt.Derived = derived;
    receipt.SemanticValidationPassed = true;
    receipt.CanonicalSha256 = TargetIdentitySha256V2(receipt, diagnostic);
    return !receipt.CanonicalSha256.empty();
}

bool CaptureNativeTargetIdentityV2(
    TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
#if WIREHAIR_RDPRU_TARGET_IDENTITY_V2_NATIVE
    NativeCompilerBarrier();
    const bool captured = CaptureTargetIdentityV2(
        kFrozenTargetCpu,
        NativeCpuidRead,
        nullptr,
        NativeContextCapture,
        nullptr,
        receipt,
        diagnostic);
    NativeCompilerBarrier();
    return captured;
#else
    receipt = TargetIdentityReceiptV2();
    receipt.RequestedCpu = kFrozenTargetCpu;
    diagnostic =
        "native target identity requires Linux x86_64 GNU-compatible CPUID";
    return false;
#endif
}

CpuidSnapshotV2 FrozenCpu50CpuidSnapshotV2()
{
    CpuidSnapshotV2 snapshot;
    const struct Fixed
    {
        uint32_t Leaf;
        CpuidRegsV2 Regs;
        CpuidRecordV2* Record;
    } fixed[] = {
        { 0u, { UINT32_C(0x00000010), UINT32_C(0x68747541),
                UINT32_C(0x444d4163), UINT32_C(0x69746e65) }, &snapshot.Leaf0 },
        { 1u, { UINT32_C(0x00b00f81), UINT32_C(0x64800800),
                UINT32_C(0x7efa320b), UINT32_C(0x178bfbff) }, &snapshot.Leaf1 },
        { 6u, { UINT32_C(0x00000004), UINT32_C(0x00000000),
                UINT32_C(0x00000001), UINT32_C(0x00000000) }, &snapshot.Leaf6 },
        { UINT32_C(0x80000000),
            { UINT32_C(0x80000028), UINT32_C(0x68747541),
              UINT32_C(0x444d4163), UINT32_C(0x69746e65) },
            &snapshot.Leaf80000000 },
        { UINT32_C(0x80000001),
            { UINT32_C(0x00b00f81), UINT32_C(0x70000000),
              UINT32_C(0x75c237ff), UINT32_C(0x2fd3fbff) },
            &snapshot.Leaf80000001 },
        { UINT32_C(0x80000008),
            { UINT32_C(0x00003934), UINT32_C(0x79bef25f),
              UINT32_C(0x0000707f), UINT32_C(0x00010007) },
            &snapshot.Leaf80000008 },
        { UINT32_C(0x8000001e),
            { UINT32_C(0x00000064), UINT32_C(0x00000132),
              UINT32_C(0x00000000), UINT32_C(0x00000000) },
            &snapshot.Leaf8000001e },
        { UINT32_C(0x80000021),
            { UINT32_C(0xd93fffcf), UINT32_C(0x00080382),
              UINT32_C(0x00000000), UINT32_C(0x00000000) },
            &snapshot.Leaf80000021 }
    };
    for (std::size_t i = 0u; i < sizeof(fixed) / sizeof(fixed[0]); ++i)
    {
        fixed[i].Record->Leaf = fixed[i].Leaf;
        fixed[i].Record->Subleaf = 0u;
        fixed[i].Record->Regs = fixed[i].Regs;
        fixed[i].Record->Captured = true;
    }
    snapshot.LeafB.Leaf = UINT32_C(0x0000000b);
    const CpuidRegsV2 b[] = {
        { 1u, 2u, UINT32_C(0x00000100), UINT32_C(0x64) },
        { 7u, UINT32_C(0x80), UINT32_C(0x00000201), UINT32_C(0x64) },
        { 0u, 0u, UINT32_C(0x00000002), UINT32_C(0x64) }
    };
    for (uint32_t i = 0u; i < 2u; ++i)
    {
        CpuidRecordV2 record;
        record.Leaf = UINT32_C(0x0000000b);
        record.Subleaf = i;
        record.Regs = b[i];
        record.Captured = true;
        snapshot.LeafB.ValidLevels.push_back(record);
    }
    snapshot.LeafB.Sentinel.Leaf = UINT32_C(0x0000000b);
    snapshot.LeafB.Sentinel.Subleaf = 2u;
    snapshot.LeafB.Sentinel.Regs = b[2];
    snapshot.LeafB.Sentinel.Captured = true;
    snapshot.LeafB.SentinelCaptured = true;

    snapshot.Leaf80000026.Leaf = UINT32_C(0x80000026);
    const CpuidRegsV2 x26[] = {
        { 1u, 2u, UINT32_C(0x00000100), UINT32_C(0x64) },
        { 4u, UINT32_C(0x10), UINT32_C(0x00000201), UINT32_C(0x64) },
        { 4u, UINT32_C(0x10), UINT32_C(0x00000302), UINT32_C(0x64) },
        { 7u, UINT32_C(0x80), UINT32_C(0x00000403), UINT32_C(0x64) },
        { 0u, 0u, UINT32_C(0x00000004), UINT32_C(0x64) }
    };
    for (uint32_t i = 0u; i < 4u; ++i)
    {
        CpuidRecordV2 record;
        record.Leaf = UINT32_C(0x80000026);
        record.Subleaf = i;
        record.Regs = x26[i];
        record.Captured = true;
        snapshot.Leaf80000026.ValidLevels.push_back(record);
    }
    snapshot.Leaf80000026.Sentinel.Leaf = UINT32_C(0x80000026);
    snapshot.Leaf80000026.Sentinel.Subleaf = 4u;
    snapshot.Leaf80000026.Sentinel.Regs = x26[4];
    snapshot.Leaf80000026.Sentinel.Captured = true;
    snapshot.Leaf80000026.SentinelCaptured = true;
    return snapshot;
}

DerivedTargetIdentityV2 FrozenCpu50DerivedIdentityV2()
{
    DerivedTargetIdentityV2 derived;
    derived.Family = 26u;
    derived.Model = 8u;
    derived.Stepping = 1u;
    derived.InitialApicId8 = UINT32_C(0x64);
    derived.FullApicId = UINT32_C(0x64);
    derived.ThreadId = 0u;
    derived.CoreId = 50u;
    derived.ComplexId = 6u;
    derived.CcdId = 6u;
    derived.PackageId = 0u;
    derived.ThreadsPerCore = 2u;
    derived.LogicalProcessorsPerPackage = 128u;
    return derived;
}

TargetIdentityReceiptV2 FrozenCpu50TargetIdentityReceiptV2()
{
    TargetIdentityReceiptV2 receipt;
    receipt.RequestedCpu = kFrozenTargetCpu;
    receipt.Before.Cpu = kFrozenTargetCpu;
    receipt.Before.Affinity.push_back(kFrozenTargetCpu);
    receipt.Before.VoluntaryContextSwitches = 0;
    receipt.Before.InvoluntaryContextSwitches = 0;
    receipt.Before.Captured = true;
    receipt.After = receipt.Before;
    receipt.Raw = FrozenCpu50CpuidSnapshotV2();
    receipt.Derived = FrozenCpu50DerivedIdentityV2();
    receipt.RawCaptureComplete = true;
    receipt.SemanticValidationPassed = true;
    std::string diagnostic;
    receipt.CanonicalSha256 = TargetIdentitySha256V2(receipt, diagnostic);
    return receipt;
}

bool SerializeTargetIdentityV2(
    const TargetIdentityReceiptV2& receipt,
    std::string& canonical_bytes,
    std::string& diagnostic)
{
    canonical_bytes.clear();
    DerivedTargetIdentityV2 derived;
    if (!ValidateTargetIdentitySemanticsV2(receipt, derived, diagnostic)) {
        return false;
    }
    if (!receipt.SemanticValidationPassed ||
        !DerivedEqual(receipt.Derived, derived))
    {
        diagnostic = "stored derived identity disagrees with raw CPUID data";
        return false;
    }
    AppendSizedString(canonical_bytes, kIdentityDomain);
    AppendLe32(canonical_bytes, kIdentityVersion);
    AppendLe32(canonical_bytes, static_cast<uint32_t>(receipt.RequestedCpu));
    AppendContextNormalized(canonical_bytes, receipt.Before);
    AppendContextNormalized(canonical_bytes, receipt.After);
    AppendLe64(canonical_bytes, 0u); // normalized voluntary switch delta
    AppendLe64(canonical_bytes, 0u); // normalized involuntary switch delta
    AppendRecord(canonical_bytes, receipt.Raw.Leaf0);
    AppendRecord(canonical_bytes, receipt.Raw.Leaf1);
    AppendRecord(canonical_bytes, receipt.Raw.Leaf6);
    AppendRecord(canonical_bytes, receipt.Raw.Leaf80000000);
    AppendRecord(canonical_bytes, receipt.Raw.Leaf80000001);
    AppendRecord(canonical_bytes, receipt.Raw.Leaf80000008);
    AppendRecord(canonical_bytes, receipt.Raw.Leaf8000001e);
    AppendRecord(canonical_bytes, receipt.Raw.Leaf80000021);
    AppendTopology(canonical_bytes, receipt.Raw.LeafB);
    AppendTopology(canonical_bytes, receipt.Raw.Leaf80000026);
    AppendLe32(canonical_bytes, derived.Family);
    AppendLe32(canonical_bytes, derived.Model);
    AppendLe32(canonical_bytes, derived.Stepping);
    AppendLe32(canonical_bytes, derived.InitialApicId8);
    AppendLe32(canonical_bytes, derived.FullApicId);
    AppendLe32(canonical_bytes, derived.ThreadId);
    AppendLe32(canonical_bytes, derived.CoreId);
    AppendLe32(canonical_bytes, derived.ComplexId);
    AppendLe32(canonical_bytes, derived.CcdId);
    AppendLe32(canonical_bytes, derived.PackageId);
    AppendLe32(canonical_bytes, derived.ThreadsPerCore);
    AppendLe32(canonical_bytes, derived.LogicalProcessorsPerPackage);
    diagnostic.clear();
    return true;
}

std::string TargetIdentitySha256V2(
    const TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
    std::string bytes;
    if (!SerializeTargetIdentityV2(receipt, bytes, diagnostic)) {
        return std::string();
    }
    return wirehair::wh2_benchmark::Sha256Hex(bytes);
}

std::string FrozenCpu50TargetIdentitySha256V2()
{
    return kFrozenIdentitySha256;
}

uint64_t FrozenCpu50TargetIdentityCanonicalBytesV2()
{
    return kFrozenIdentityCanonicalBytes;
}

std::string FormatTargetIdentityEvidenceV2(
    const std::string& checkpoint,
    const TargetIdentityReceiptV2& receipt,
    bool gate_passed,
    const std::string& diagnostic)
{
    std::ostringstream evidence;
    evidence.imbue(std::locale::classic());
    const std::string before_affinity = FormatAffinity(receipt.Before.Affinity);
    const std::string after_affinity = FormatAffinity(receipt.After.Affinity);
    const bool canonical_sha_safe =
        IsLowerHexSha256Value(receipt.CanonicalSha256);
    evidence << "TARGET_CONTEXT_V2,checkpoint_length=" << checkpoint.size()
             << ",checkpoint_hex=" << HexEncode(checkpoint)
             << ",target_cpu=" << receipt.RequestedCpu
             << ",before_captured=" << (receipt.Before.Captured ? 1u : 0u)
             << ",before_cpu=" << receipt.Before.Cpu
             << ",before_affinity_count=" << receipt.Before.Affinity.size()
             << ",before_affinity_length=" << before_affinity.size()
             << ",before_affinity_hex=" << HexEncode(before_affinity)
             << ",after_captured=" << (receipt.After.Captured ? 1u : 0u)
             << ",after_cpu=" << receipt.After.Cpu
             << ",after_affinity_count=" << receipt.After.Affinity.size()
             << ",after_affinity_length=" << after_affinity.size()
             << ",after_affinity_hex=" << HexEncode(after_affinity)
             << ",before_voluntary=" << receipt.Before.VoluntaryContextSwitches
             << ",after_voluntary=" << receipt.After.VoluntaryContextSwitches
             << ",before_involuntary="
             << receipt.Before.InvoluntaryContextSwitches
             << ",after_involuntary="
             << receipt.After.InvoluntaryContextSwitches << '\n';
    const CpuidRecordV2* const fixed[] = {
        &receipt.Raw.Leaf0, &receipt.Raw.Leaf1, &receipt.Raw.Leaf6,
        &receipt.Raw.Leaf80000000, &receipt.Raw.Leaf80000001,
        &receipt.Raw.Leaf80000008, &receipt.Raw.Leaf8000001e,
        &receipt.Raw.Leaf80000021
    };
    for (std::size_t i = 0u; i < sizeof(fixed) / sizeof(fixed[0]); ++i) {
        AppendEvidenceRecord(evidence, checkpoint, *fixed[i], "fixed");
    }
    for (std::size_t i = 0u; i < receipt.Raw.LeafB.ValidLevels.size(); ++i) {
        AppendEvidenceRecord(
            evidence, checkpoint, receipt.Raw.LeafB.ValidLevels[i], "valid-b");
    }
    AppendEvidenceRecord(
        evidence, checkpoint, receipt.Raw.LeafB.Sentinel, "sentinel-b");
    for (std::size_t i = 0u;
         i < receipt.Raw.Leaf80000026.ValidLevels.size();
         ++i)
    {
        AppendEvidenceRecord(
            evidence,
            checkpoint,
            receipt.Raw.Leaf80000026.ValidLevels[i],
            "valid-80000026");
    }
    AppendEvidenceRecord(
        evidence,
        checkpoint,
        receipt.Raw.Leaf80000026.Sentinel,
        "sentinel-80000026");
    evidence << "TARGET_IDENTITY_V2,checkpoint_length=" << checkpoint.size()
             << ",checkpoint_hex=" << HexEncode(checkpoint)
             << ",raw_complete=" << (receipt.RawCaptureComplete ? 1u : 0u)
             << ",semantic_pass="
             << (receipt.SemanticValidationPassed ? 1u : 0u)
             << ",full_apic_id=0x" << HexU32(receipt.Derived.FullApicId)
             << ",initial_apic_id=0x"
             << HexU32(receipt.Derived.InitialApicId8)
             << ",thread=" << receipt.Derived.ThreadId
             << ",core=" << receipt.Derived.CoreId
             << ",complex=" << receipt.Derived.ComplexId
             << ",ccd=" << receipt.Derived.CcdId
             << ",package=" << receipt.Derived.PackageId
             << ",canonical_sha256_valid=" << (canonical_sha_safe ? 1u : 0u)
             << ",canonical_sha256="
             << (canonical_sha_safe ? receipt.CanonicalSha256 : std::string())
             << ",canonical_sha256_length=" << receipt.CanonicalSha256.size()
             << ",canonical_sha256_hex="
             << HexEncode(receipt.CanonicalSha256)
             << ",diagnostic_length=" << diagnostic.size()
             << ",diagnostic_hex=" << HexEncode(diagnostic)
             << ",gate=" << (gate_passed ? "PASS" : "FAIL") << '\n';
    return evidence.str();
}

bool ParseProcCpuInfoV2(
    const std::string& content,
    std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic)
{
    records.clear();
    std::vector<CpuInfoMapRecordV2> parsed;
    if (content.empty() || content.size() > kMaximumCpuInfoBytes ||
        content.find('\0') != std::string::npos ||
        content.find('\r') != std::string::npos)
    {
        diagnostic =
            "cpuinfo content is empty, oversized, or contains NUL/CR";
        return false;
    }
    struct Pending
    {
        CpuInfoMapRecordV2 Record;
        bool Cpu = false;
        bool Package = false;
        bool Core = false;
        bool Apic = false;
        bool Initial = false;
        bool Any = false;
    } pending;
    const auto finish_record = [&parsed, &pending, &diagnostic]() -> bool {
        if (!pending.Any) return true;
        if (!pending.Cpu || !pending.Package || !pending.Core ||
            !pending.Apic || !pending.Initial)
        {
            diagnostic = "cpuinfo record is missing a required topology field";
            return false;
        }
        if (parsed.size() >= kFrozenCpuInfoRecords)
        {
            diagnostic = "cpuinfo contains more than 128 records";
            return false;
        }
        parsed.push_back(pending.Record);
        pending = Pending();
        return true;
    };
    std::size_t offset = 0u;
    while (offset < content.size())
    {
        const std::size_t newline = content.find('\n', offset);
        const std::size_t end = newline == std::string::npos ?
            content.size() : newline;
        if (end - offset > kMaximumCpuInfoLineBytes)
        {
            diagnostic = "cpuinfo line is oversized";
            return false;
        }
        const std::string line = content.substr(offset, end - offset);
        offset = newline == std::string::npos ? content.size() : newline + 1u;
        if (Trim(line).empty())
        {
            if (!finish_record()) return false;
            continue;
        }
        pending.Any = true;
        const std::size_t colon = line.find(':');
        if (colon == std::string::npos)
        {
            diagnostic = "cpuinfo line is missing a colon";
            return false;
        }
        const std::string key = Trim(line.substr(0u, colon));
        const std::string value = line.substr(colon + 1u);
        bool* seen = nullptr;
        uint32_t* destination = nullptr;
        if (key == "processor") {
            seen = &pending.Cpu; destination = &pending.Record.Cpu;
        }
        else if (key == "physical id") {
            seen = &pending.Package; destination = &pending.Record.Package;
        }
        else if (key == "core id") {
            seen = &pending.Core; destination = &pending.Record.Core;
        }
        else if (key == "apicid") {
            seen = &pending.Apic; destination = &pending.Record.ApicId;
        }
        else if (key == "initial apicid") {
            seen = &pending.Initial;
            destination = &pending.Record.InitialApicId;
        }
        if (!seen) continue;
        if (*seen || !ParseDecimalU32(value, *destination))
        {
            diagnostic = "cpuinfo topology field is duplicate or malformed";
            return false;
        }
        *seen = true;
    }
    if (!finish_record()) return false;
    if (parsed.empty())
    {
        diagnostic = "cpuinfo contains no topology records";
        return false;
    }
    std::set<uint32_t> cpus;
    std::set<uint32_t> apic_ids;
    std::set<uint32_t> initial_apic_ids;
    for (std::size_t i = 0u; i < parsed.size(); ++i)
    {
        if (!cpus.insert(parsed[i].Cpu).second ||
            !apic_ids.insert(parsed[i].ApicId).second ||
            !initial_apic_ids.insert(parsed[i].InitialApicId).second)
        {
            diagnostic = "cpuinfo contains a duplicate CPU or APIC identity";
            return false;
        }
    }
    records.swap(parsed);
    diagnostic.clear();
    return true;
}

bool CanonicalizeCpuInfoMapV2(
    const std::vector<CpuInfoMapRecordV2>& records,
    std::string& canonical,
    std::string& diagnostic)
{
    canonical.clear();
    if (records.empty() || records.size() > kFrozenCpuInfoRecords)
    {
        diagnostic = "CPU/APIC map record count is outside the bounded domain";
        return false;
    }
    std::vector<CpuInfoMapRecordV2> sorted(records);
    std::sort(
        sorted.begin(), sorted.end(),
        [](const CpuInfoMapRecordV2& a, const CpuInfoMapRecordV2& b) {
            return a.Cpu < b.Cpu;
        });
    std::set<uint32_t> apic_ids;
    std::set<uint32_t> initial_apic_ids;
    for (std::size_t i = 0u; i < sorted.size(); ++i)
    {
        if ((i != 0u && sorted[i - 1u].Cpu == sorted[i].Cpu) ||
            !apic_ids.insert(sorted[i].ApicId).second ||
            !initial_apic_ids.insert(sorted[i].InitialApicId).second)
        {
            diagnostic = "CPU/APIC map contains a duplicate CPU or APIC ID";
            canonical.clear();
            return false;
        }
        std::ostringstream line;
        line.imbue(std::locale::classic());
        line << "cpu=" << sorted[i].Cpu
             << ",package=" << sorted[i].Package
             << ",core=" << sorted[i].Core
             << ",apicid=" << sorted[i].ApicId
             // This spelling is part of the preregistered 6,818-byte/d51bf8
             // canonical domain.  Changing it to initial_apic_id produces a
             // different 6,946-byte stream and must fail the golden test.
             << ",initial_apicid=" << sorted[i].InitialApicId << '\n';
        canonical += line.str();
        if (canonical.size() > kMaximumCpuInfoBytes)
        {
            diagnostic = "canonical CPU/APIC map is oversized";
            canonical.clear();
            return false;
        }
    }
    diagnostic.clear();
    return true;
}

std::string CpuInfoMapSha256V2(
    const std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic)
{
    std::string canonical;
    if (!CanonicalizeCpuInfoMapV2(records, canonical, diagnostic)) {
        return std::string();
    }
    return wirehair::wh2_benchmark::Sha256Hex(canonical);
}

bool ReadProcCpuInfoV2(
    const std::string& path,
    std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic)
{
#if defined(__linux__)
    records.clear();
    if (path != "/proc/cpuinfo")
    {
        diagnostic = "cpuinfo path is not the frozen /proc/cpuinfo path";
        return false;
    }
    const int file = open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
    if (file < 0)
    {
        diagnostic = std::string("open cpuinfo failed: ") + std::strerror(errno);
        return false;
    }
    struct stat before;
    if (fstat(file, &before) != 0 || !S_ISREG(before.st_mode))
    {
        diagnostic = "cpuinfo descriptor is not a regular procfs file";
        (void)close(file);
        return false;
    }
    struct statfs filesystem;
    if (fstatfs(file, &filesystem) != 0 ||
        static_cast<uint64_t>(filesystem.f_type) != UINT64_C(0x9fa0))
    {
        diagnostic = "cpuinfo descriptor is not on procfs";
        (void)close(file);
        return false;
    }
    std::string content;
    char buffer[4096];
    bool ok = true;
    while (true)
    {
        const ssize_t got = read(file, buffer, sizeof(buffer));
        if (got < 0)
        {
            if (errno == EINTR) continue;
            diagnostic = std::string("read cpuinfo failed: ") +
                std::strerror(errno);
            ok = false;
            break;
        }
        if (got == 0) break;
        if (static_cast<std::size_t>(got) >
            kMaximumCpuInfoBytes - content.size())
        {
            diagnostic = "cpuinfo exceeds the bounded input size";
            ok = false;
            break;
        }
        content.append(buffer, static_cast<std::size_t>(got));
    }
    struct stat after;
    const bool same_descriptor = fstat(file, &after) == 0 &&
        before.st_dev == after.st_dev && before.st_ino == after.st_ino;
    const bool close_ok = close(file) == 0;
    if (!ok || !same_descriptor || !close_ok)
    {
        if (ok) diagnostic = same_descriptor ?
            "close cpuinfo failed" : "cpuinfo descriptor identity changed";
        return false;
    }
    return ParseProcCpuInfoV2(content, records, diagnostic);
#else
    (void)path;
    records.clear();
    diagnostic = "cpuinfo map reading requires Linux";
    return false;
#endif
}

std::vector<CpuInfoMapRecordV2> FrozenCpuInfoMapV2()
{
    std::vector<CpuInfoMapRecordV2> records;
    records.reserve(kFrozenCpuInfoRecords);
    for (uint32_t cpu = 0u; cpu < 128u; ++cpu)
    {
        CpuInfoMapRecordV2 record;
        record.Cpu = cpu;
        record.Package = 0u;
        record.Core = cpu & 63u;
        record.ApicId = (record.Core << 1u) | (cpu >= 64u ? 1u : 0u);
        record.InitialApicId = record.ApicId;
        records.push_back(record);
    }
    return records;
}

bool ValidateFrozenCpuInfoMapV2(
    const std::vector<CpuInfoMapRecordV2>& records,
    std::string& diagnostic)
{
    const std::vector<CpuInfoMapRecordV2> expected = FrozenCpuInfoMapV2();
    if (records.size() != expected.size())
    {
        diagnostic = "CPU/APIC map does not contain exactly 128 records";
        return false;
    }
    std::vector<CpuInfoMapRecordV2> sorted(records);
    std::sort(
        sorted.begin(), sorted.end(),
        [](const CpuInfoMapRecordV2& a, const CpuInfoMapRecordV2& b) {
            return a.Cpu < b.Cpu;
        });
    for (std::size_t i = 0u; i < expected.size(); ++i) {
        if (!CpuInfoRecordEqual(sorted[i], expected[i]))
        {
            diagnostic = "CPU/APIC map differs from the frozen 128-entry map";
            return false;
        }
    }
    std::string canonical;
    if (!CanonicalizeCpuInfoMapV2(sorted, canonical, diagnostic)) return false;
    if (canonical.size() != kFrozenCpuInfoCanonicalBytes ||
        wirehair::wh2_benchmark::Sha256Hex(canonical) != kFrozenCpuInfoSha256)
    {
        diagnostic = "CPU/APIC map canonical length or SHA-256 mismatch";
        return false;
    }
    diagnostic.clear();
    return true;
}

const char* FrozenCpuInfoMapSha256V2()
{
    return kFrozenCpuInfoSha256;
}

uint64_t FrozenCpuInfoMapCanonicalBytesV2()
{
    return kFrozenCpuInfoCanonicalBytes;
}

namespace {

struct FakeCpuidState
{
    std::map<std::pair<uint32_t, uint32_t>, CpuidRegsV2> Values;
    std::vector<std::pair<uint32_t, uint32_t> > Calls;
    std::size_t FailAt = std::numeric_limits<std::size_t>::max();
    std::size_t ThrowAt = std::numeric_limits<std::size_t>::max();
};

bool FakeCpuidRead(
    void* context,
    uint32_t leaf,
    uint32_t subleaf,
    CpuidRegsV2& registers,
    std::string& diagnostic)
{
    FakeCpuidState* const state = static_cast<FakeCpuidState*>(context);
    if (!state)
    {
        diagnostic = "fake CPUID state is null";
        return false;
    }
    const std::size_t call = state->Calls.size();
    state->Calls.push_back(std::make_pair(leaf, subleaf));
    if (call == state->ThrowAt) {
        throw std::runtime_error("injected CPUID exception");
    }
    if (call == state->FailAt)
    {
        diagnostic = "injected CPUID failure";
        return false;
    }
    const std::map<std::pair<uint32_t, uint32_t>, CpuidRegsV2>::const_iterator
        found = state->Values.find(std::make_pair(leaf, subleaf));
    if (found == state->Values.end())
    {
        diagnostic = "unexpected CPUID leaf/subleaf";
        return false;
    }
    registers = found->second;
    diagnostic.clear();
    return true;
}

void AddFakeRecord(FakeCpuidState& state, const CpuidRecordV2& record)
{
    state.Values[std::make_pair(record.Leaf, record.Subleaf)] = record.Regs;
}

FakeCpuidState FrozenFakeCpuidState()
{
    FakeCpuidState state;
    const CpuidSnapshotV2 snapshot = FrozenCpu50CpuidSnapshotV2();
    const CpuidRecordV2* const fixed[] = {
        &snapshot.Leaf0, &snapshot.Leaf80000000, &snapshot.Leaf1,
        &snapshot.Leaf6, &snapshot.Leaf80000001, &snapshot.Leaf80000008,
        &snapshot.Leaf8000001e, &snapshot.Leaf80000021
    };
    for (std::size_t i = 0u; i < sizeof(fixed) / sizeof(fixed[0]); ++i) {
        AddFakeRecord(state, *fixed[i]);
    }
    for (std::size_t i = 0u; i < snapshot.LeafB.ValidLevels.size(); ++i) {
        AddFakeRecord(state, snapshot.LeafB.ValidLevels[i]);
    }
    AddFakeRecord(state, snapshot.LeafB.Sentinel);
    for (std::size_t i = 0u;
         i < snapshot.Leaf80000026.ValidLevels.size();
         ++i) {
        AddFakeRecord(state, snapshot.Leaf80000026.ValidLevels[i]);
    }
    AddFakeRecord(state, snapshot.Leaf80000026.Sentinel);
    return state;
}

std::vector<std::pair<uint32_t, uint32_t> > FrozenCpuidCallOrder()
{
    const uint32_t pairs[][2] = {
        { 0u, 0u }, { UINT32_C(0x80000000), 0u },
        { 1u, 0u }, { 6u, 0u },
        { UINT32_C(0x80000001), 0u },
        { UINT32_C(0x80000008), 0u },
        { UINT32_C(0x8000001e), 0u },
        { UINT32_C(0x80000021), 0u },
        { UINT32_C(0x0000000b), 0u },
        { UINT32_C(0x0000000b), 1u },
        { UINT32_C(0x0000000b), 2u },
        { UINT32_C(0x80000026), 0u },
        { UINT32_C(0x80000026), 1u },
        { UINT32_C(0x80000026), 2u },
        { UINT32_C(0x80000026), 3u },
        { UINT32_C(0x80000026), 4u }
    };
    std::vector<std::pair<uint32_t, uint32_t> > order;
    for (std::size_t i = 0u; i < sizeof(pairs) / sizeof(pairs[0]); ++i) {
        order.push_back(std::make_pair(pairs[i][0], pairs[i][1]));
    }
    return order;
}

TargetContextReceiptV2 ValidFakeContext(int32_t cpu, int64_t switches)
{
    TargetContextReceiptV2 receipt;
    receipt.Cpu = cpu;
    receipt.Affinity.push_back(cpu);
    receipt.VoluntaryContextSwitches = switches;
    receipt.InvoluntaryContextSwitches = switches / 2;
    receipt.Captured = true;
    return receipt;
}

struct FakeContextState
{
    TargetContextReceiptV2 Before;
    TargetContextReceiptV2 After;
    std::size_t Calls = 0u;
    std::size_t FailAt = std::numeric_limits<std::size_t>::max();
    std::size_t ThrowAt = std::numeric_limits<std::size_t>::max();
};

bool FakeContextCapture(
    void* context,
    int32_t,
    TargetContextReceiptV2& receipt,
    std::string& diagnostic)
{
    FakeContextState* const state = static_cast<FakeContextState*>(context);
    if (!state)
    {
        diagnostic = "fake context state is null";
        return false;
    }
    const std::size_t call = state->Calls++;
    if (call == state->ThrowAt) {
        throw std::runtime_error("injected context exception");
    }
    if (call == state->FailAt)
    {
        diagnostic = "injected context failure";
        return false;
    }
    receipt = call == 0u ? state->Before : state->After;
    diagnostic.clear();
    return true;
}

FakeContextState FrozenFakeContextState()
{
    FakeContextState state;
    state.Before = ValidFakeContext(kFrozenTargetCpu, 20);
    state.After = state.Before;
    return state;
}

bool CaptureFrozenFake(
    TargetIdentityReceiptV2& receipt,
    FakeCpuidState& cpuid,
    FakeContextState& context,
    std::string& diagnostic)
{
    return CaptureTargetIdentityV2(
        kFrozenTargetCpu,
        FakeCpuidRead,
        &cpuid,
        FakeContextCapture,
        &context,
        receipt,
        diagnostic);
}

CpuidRecordV2* MutableRecordAt(CpuidSnapshotV2& snapshot, std::size_t index)
{
    CpuidRecordV2* const fixed[] = {
        &snapshot.Leaf0, &snapshot.Leaf1, &snapshot.Leaf6,
        &snapshot.Leaf80000000, &snapshot.Leaf80000001,
        &snapshot.Leaf80000008, &snapshot.Leaf8000001e,
        &snapshot.Leaf80000021
    };
    if (index < sizeof(fixed) / sizeof(fixed[0])) return fixed[index];
    index -= sizeof(fixed) / sizeof(fixed[0]);
    if (index < snapshot.LeafB.ValidLevels.size()) {
        return &snapshot.LeafB.ValidLevels[index];
    }
    index -= snapshot.LeafB.ValidLevels.size();
    if (index == 0u) return &snapshot.LeafB.Sentinel;
    --index;
    if (index < snapshot.Leaf80000026.ValidLevels.size()) {
        return &snapshot.Leaf80000026.ValidLevels[index];
    }
    index -= snapshot.Leaf80000026.ValidLevels.size();
    return index == 0u ? &snapshot.Leaf80000026.Sentinel : nullptr;
}

std::size_t CapturedRecordCount(const CpuidSnapshotV2& snapshot)
{
    const CpuidRecordV2* const fixed[] = {
        &snapshot.Leaf0, &snapshot.Leaf80000000, &snapshot.Leaf1,
        &snapshot.Leaf6, &snapshot.Leaf80000001, &snapshot.Leaf80000008,
        &snapshot.Leaf8000001e, &snapshot.Leaf80000021
    };
    std::size_t count = 0u;
    for (std::size_t i = 0u; i < sizeof(fixed) / sizeof(fixed[0]); ++i) {
        if (fixed[i]->Captured) ++count;
    }
    for (std::size_t i = 0u; i < snapshot.LeafB.ValidLevels.size(); ++i) {
        if (snapshot.LeafB.ValidLevels[i].Captured) ++count;
    }
    if (snapshot.LeafB.SentinelCaptured && snapshot.LeafB.Sentinel.Captured) {
        ++count;
    }
    for (std::size_t i = 0u;
         i < snapshot.Leaf80000026.ValidLevels.size();
         ++i) {
        if (snapshot.Leaf80000026.ValidLevels[i].Captured) ++count;
    }
    if (snapshot.Leaf80000026.SentinelCaptured &&
        snapshot.Leaf80000026.Sentinel.Captured) {
        ++count;
    }
    return count;
}

uint32_t* MutableDerivedWordAt(
    DerivedTargetIdentityV2& derived,
    std::size_t index)
{
    uint32_t* const words[] = {
        &derived.Family, &derived.Model, &derived.Stepping,
        &derived.InitialApicId8, &derived.FullApicId, &derived.ThreadId,
        &derived.CoreId, &derived.ComplexId, &derived.CcdId,
        &derived.PackageId, &derived.ThreadsPerCore,
        &derived.LogicalProcessorsPerPackage
    };
    return index < sizeof(words) / sizeof(words[0]) ? words[index] : nullptr;
}

const std::size_t kFrozenRecordCount = 8u + 2u + 1u + 4u + 1u;

bool RefreshSyntheticReceipt(
    TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
    DerivedTargetIdentityV2 derived;
    receipt.RawCaptureComplete = true;
    if (!ValidateTargetIdentitySemanticsV2(receipt, derived, diagnostic)) {
        return false;
    }
    receipt.Derived = derived;
    receipt.SemanticValidationPassed = true;
    receipt.CanonicalSha256.clear();
    receipt.CanonicalSha256 = TargetIdentitySha256V2(receipt, diagnostic);
    return !receipt.CanonicalSha256.empty();
}

bool TestFrozenIdentityBaseline(std::string& diagnostic)
{
    FakeCpuidState cpuid = FrozenFakeCpuidState();
    FakeContextState context = FrozenFakeContextState();
    TargetIdentityReceiptV2 receipt;
    if (!CaptureFrozenFake(receipt, cpuid, context, diagnostic))
    {
        diagnostic = "fake capture: " + diagnostic;
        return false;
    }
    if (cpuid.Calls != FrozenCpuidCallOrder() || context.Calls != 2u)
    {
        diagnostic = "frozen identity capture call order mismatch";
        return false;
    }
    if (!ValidateFrozenCpu50IdentityV2(receipt, diagnostic))
    {
        diagnostic = "frozen exact validation: " + diagnostic;
        return false;
    }
    if (!IsLowerHexSha256Value(receipt.CanonicalSha256) ||
        receipt.CanonicalSha256 != FrozenCpu50TargetIdentitySha256V2())
    {
        diagnostic = "frozen identity SHA-256 mismatch";
        return false;
    }
    std::string canonical;
    if (!SerializeTargetIdentityV2(receipt, canonical, diagnostic) ||
        canonical.size() != FrozenCpu50TargetIdentityCanonicalBytesV2()) {
        if (diagnostic.empty()) diagnostic = "identity canonical length mismatch";
        return false;
    }
    TargetIdentityReceiptV2 shifted = receipt;
    shifted.Before.VoluntaryContextSwitches += 100;
    shifted.After.VoluntaryContextSwitches += 100;
    shifted.Before.InvoluntaryContextSwitches += 50;
    shifted.After.InvoluntaryContextSwitches += 50;
    shifted.CanonicalSha256.clear();
    const std::string shifted_hash = TargetIdentitySha256V2(shifted, diagnostic);
    if (shifted_hash != receipt.CanonicalSha256)
    {
        diagnostic = "absolute context counters changed normalized identity hash";
        return false;
    }
    return true;
}

bool TestFullWidthApicSemantics(std::string& diagnostic)
{
    TargetIdentityReceiptV2 receipt = FrozenCpu50TargetIdentityReceiptV2();
    receipt.Raw.Leaf1.Regs.Ebx =
        (receipt.Raw.Leaf1.Regs.Ebx & UINT32_C(0x00ffffff)) |
        UINT32_C(0x64000000);
    receipt.Raw.Leaf8000001e.Regs.Eax = UINT32_C(0x164);
    for (std::size_t i = 0u; i < receipt.Raw.LeafB.ValidLevels.size(); ++i) {
        receipt.Raw.LeafB.ValidLevels[i].Regs.Edx = UINT32_C(0x164);
    }
    receipt.Raw.LeafB.Sentinel.Regs.Edx = UINT32_C(0x164);
    for (std::size_t i = 0u;
         i < receipt.Raw.Leaf80000026.ValidLevels.size();
         ++i) {
        receipt.Raw.Leaf80000026.ValidLevels[i].Regs.Edx = UINT32_C(0x164);
    }
    receipt.Raw.Leaf80000026.Sentinel.Regs.Edx = UINT32_C(0x164);
    if (!RefreshSyntheticReceipt(receipt, diagnostic) ||
        receipt.Derived.FullApicId != UINT32_C(0x164) ||
        receipt.Derived.InitialApicId8 != UINT32_C(0x64) ||
        receipt.Derived.ThreadId != 0u || receipt.Derived.CoreId != 50u ||
        receipt.Derived.PackageId != 2u)
    {
        if (diagnostic.empty()) diagnostic = "full-width APIC derivation failed";
        return false;
    }
    std::string exact_diagnostic;
    if (ValidateFrozenCpu50IdentityV2(receipt, exact_diagnostic))
    {
        diagnostic = "non-frozen full-width APIC identity was accepted as CPU50";
        return false;
    }
    return true;
}

bool TestFrozenRawMutationSurface(std::string& diagnostic)
{
    const TargetIdentityReceiptV2 baseline =
        FrozenCpu50TargetIdentityReceiptV2();
    for (std::size_t record_index = 0u;
         record_index < kFrozenRecordCount;
         ++record_index)
    {
        for (unsigned field = 0u; field < 4u; ++field)
        {
            for (unsigned bit = 0u; bit < 32u; ++bit)
            {
                TargetIdentityReceiptV2 changed = baseline;
                CpuidRecordV2* const record =
                    MutableRecordAt(changed.Raw, record_index);
                if (!record)
                {
                    diagnostic = "raw mutation record lookup failed";
                    return false;
                }
                uint32_t* words[] = {
                    &record->Regs.Eax, &record->Regs.Ebx,
                    &record->Regs.Ecx, &record->Regs.Edx
                };
                *words[field] ^= UINT32_C(1) << bit;
                std::string rejected;
                if (ValidateFrozenCpu50IdentityV2(changed, rejected))
                {
                    diagnostic = "frozen CPUID register bit mutation was accepted";
                    return false;
                }
            }
        }
        TargetIdentityReceiptV2 changed = baseline;
        MutableRecordAt(changed.Raw, record_index)->Captured = false;
        std::string rejected;
        if (ValidateFrozenCpu50IdentityV2(changed, rejected))
        {
            diagnostic = "missing CPUID record was accepted";
            return false;
        }
        changed = baseline;
        MutableRecordAt(changed.Raw, record_index)->Leaf ^= 1u;
        if (ValidateFrozenCpu50IdentityV2(changed, rejected))
        {
            diagnostic = "CPUID leaf-tag mutation was accepted";
            return false;
        }
        changed = baseline;
        MutableRecordAt(changed.Raw, record_index)->Subleaf ^= 1u;
        if (ValidateFrozenCpu50IdentityV2(changed, rejected))
        {
            diagnostic = "CPUID subleaf-tag mutation was accepted";
            return false;
        }
    }
    TargetIdentityReceiptV2 changed = baseline;
    changed.RawCaptureComplete = false;
    if (ValidateFrozenCpu50IdentityV2(changed, diagnostic)) return false;
    changed = baseline;
    changed.SemanticValidationPassed = false;
    if (ValidateFrozenCpu50IdentityV2(changed, diagnostic)) return false;
    for (std::size_t field = 0u; field < 12u; ++field)
    {
        for (unsigned bit = 0u; bit < 32u; ++bit)
        {
            changed = baseline;
            uint32_t* const word = MutableDerivedWordAt(changed.Derived, field);
            if (!word)
            {
                diagnostic = "derived mutation field lookup failed";
                return false;
            }
            *word ^= UINT32_C(1) << bit;
            std::string rejected;
            if (ValidateFrozenCpu50IdentityV2(changed, rejected))
            {
                diagnostic = "stored derived-identity bit mutation was accepted";
                return false;
            }
        }
    }
    changed = baseline;
    changed.CanonicalSha256[0] =
        changed.CanonicalSha256[0] == '0' ? '1' : '0';
    if (ValidateFrozenCpu50IdentityV2(changed, diagnostic)) return false;
    changed = baseline;
    changed.RequestedCpu = 51;
    if (ValidateFrozenCpu50IdentityV2(changed, diagnostic)) return false;
    diagnostic.clear();
    return true;
}

bool SemanticsRejects(
    const TargetIdentityReceiptV2& receipt)
{
    DerivedTargetIdentityV2 derived;
    std::string diagnostic;
    return !ValidateTargetIdentitySemanticsV2(receipt, derived, diagnostic) &&
        !diagnostic.empty();
}

bool TestSemanticMutationSurface(std::string& diagnostic)
{
    const TargetIdentityReceiptV2 valid = FrozenCpu50TargetIdentityReceiptV2();
    TargetIdentityReceiptV2 changed = valid;
    changed.Raw.Leaf0.Regs.Ebx ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf0.Regs.Eax = UINT32_C(0x0a);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid;
    changed.Raw.Leaf80000000.Regs.Eax = UINT32_C(0x80000025);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf1.Regs.Ecx |= UINT32_C(1) << 31u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf6.Regs.Ecx &= ~UINT32_C(1);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid;
    changed.Raw.Leaf80000001.Regs.Ecx &= ~(UINT32_C(1) << 22u);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid;
    changed.Raw.Leaf80000008.Regs.Ebx &= ~(UINT32_C(1) << 4u);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000008.Regs.Edx &= UINT32_C(0x0000ffff);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid;
    changed.Raw.Leaf80000021.Regs.Eax &= ~(UINT32_C(1) << 2u);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.LeafB.ValidLevels[1].Regs.Edx ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.LeafB.Sentinel.Regs.Edx ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000026.ValidLevels[2].Regs.Edx ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000026.Sentinel.Regs.Edx ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf8000001e.Regs.Eax ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf1.Regs.Ebx ^= UINT32_C(0x01000000);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000026.ValidLevels[2].Regs.Eax = 3u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000026.ValidLevels[2].Regs.Ebx = 8u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000026.ValidLevels[1].Regs.Ebx = 17u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.LeafB.ValidLevels[1].Regs.Eax = 6u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.LeafB.Sentinel.Regs.Ecx = UINT32_C(0x00000102);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000026.SentinelCaptured = false;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.LeafB.ValidLevels.pop_back();
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid;
    std::swap(changed.Raw.LeafB.ValidLevels[0],
              changed.Raw.LeafB.ValidLevels[1]);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid;
    changed.Raw.Leaf80000026.ValidLevels[2] =
        changed.Raw.Leaf80000026.ValidLevels[1];
    changed.Raw.Leaf80000026.ValidLevels[2].Subleaf = 2u;
    changed.Raw.Leaf80000026.ValidLevels[2].Regs.Ecx =
        UINT32_C(0x00000202);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid;
    changed.Raw.Leaf80000026.ValidLevels.push_back(
        changed.Raw.Leaf80000026.Sentinel);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.LeafB.Leaf ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000026.Leaf ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf8000001e.Regs.Ebx ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf1.Regs.Ebx ^= UINT32_C(0x00010000);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000008.Regs.Ecx ^= 1u;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Raw.Leaf80000008.Regs.Ecx ^= UINT32_C(0x1000);
    if (!SemanticsRejects(changed)) goto failure;
    // The frozen leaf intentionally has equal level-1/level-2 widths/counts.
    {
        DerivedTargetIdentityV2 derived;
        if (!ValidateTargetIdentitySemanticsV2(valid, derived, diagnostic)) {
            return false;
        }
    }
    return true;

failure:
    diagnostic = "semantic CPUID/topology mutation was accepted";
    return false;
}

bool TestContextMutationSurface(std::string& diagnostic)
{
    const TargetIdentityReceiptV2 valid = FrozenCpu50TargetIdentityReceiptV2();
    TargetIdentityReceiptV2 changed = valid;
    changed.Before.Cpu = 51;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.After.Cpu = 51;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Before.Affinity.push_back(51);
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.After.Affinity.clear();
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Before.Captured = false;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; ++changed.After.VoluntaryContextSwitches;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; ++changed.After.InvoluntaryContextSwitches;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.After.VoluntaryContextSwitches = -1;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.After.InvoluntaryContextSwitches = -1;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; --changed.After.VoluntaryContextSwitches;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; --changed.After.InvoluntaryContextSwitches;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Before.VoluntaryContextSwitches = -1;
    if (!SemanticsRejects(changed)) goto failure;
    changed = valid; changed.Before.InvoluntaryContextSwitches = -1;
    if (!SemanticsRejects(changed)) goto failure;
    return true;

failure:
    diagnostic = "context mutation was accepted";
    return false;
}

bool TestInjectedCaptureFailures(std::string& diagnostic)
{
    const std::vector<std::pair<uint32_t, uint32_t> > order =
        FrozenCpuidCallOrder();
    for (std::size_t fail = 0u; fail < order.size(); ++fail)
    {
        for (unsigned throws = 0u; throws < 2u; ++throws)
        {
            FakeCpuidState cpuid = FrozenFakeCpuidState();
            if (throws == 0u) cpuid.FailAt = fail;
            else cpuid.ThrowAt = fail;
            FakeContextState context = FrozenFakeContextState();
            TargetIdentityReceiptV2 receipt;
            std::string rejected;
            if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
                context.Calls != 2u || rejected.empty() ||
                cpuid.Calls.size() != fail + 1u ||
                !receipt.After.Captured || receipt.RawCaptureComplete ||
                CapturedRecordCount(receipt.Raw) != fail)
            {
                diagnostic = "injected CPUID failure did not preserve bracket";
                return false;
            }
        }
    }
    for (std::size_t context_call = 0u; context_call < 2u; ++context_call)
    {
        for (unsigned throws = 0u; throws < 2u; ++throws)
        {
            FakeCpuidState cpuid = FrozenFakeCpuidState();
            FakeContextState context = FrozenFakeContextState();
            if (throws == 0u) context.FailAt = context_call;
            else context.ThrowAt = context_call;
            TargetIdentityReceiptV2 receipt;
            std::string rejected;
            if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
                context.Calls != context_call + 1u || rejected.empty() ||
                (context_call == 0u && !cpuid.Calls.empty()) ||
                (context_call == 1u && cpuid.Calls != order))
            {
                diagnostic = "injected context failure call order mismatch";
                return false;
            }
        }
    }
    {
        FakeCpuidState cpuid = FrozenFakeCpuidState();
        FakeContextState context = FrozenFakeContextState();
        TargetIdentityReceiptV2 receipt;
        std::string rejected;
        if (CaptureTargetIdentityV2(
                kFrozenTargetCpu,
                nullptr,
                &cpuid,
                FakeContextCapture,
                &context,
                receipt,
                rejected) ||
            context.Calls != 0u || !cpuid.Calls.empty() || rejected.empty())
        {
            diagnostic = "null CPUID callback did not fail before capture";
            return false;
        }
        if (CaptureTargetIdentityV2(
                kFrozenTargetCpu,
                FakeCpuidRead,
                &cpuid,
                nullptr,
                &context,
                receipt,
                rejected) ||
            context.Calls != 0u || !cpuid.Calls.empty() || rejected.empty())
        {
            diagnostic = "null context callback did not fail before capture";
            return false;
        }
    }
    {
        FakeCpuidState cpuid = FrozenFakeCpuidState();
        FakeContextState context = FrozenFakeContextState();
        context.Before.Affinity.push_back(51);
        TargetIdentityReceiptV2 receipt;
        std::string rejected;
        if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
            context.Calls != 1u || !cpuid.Calls.empty() || rejected.empty())
        {
            diagnostic = "invalid before-context did not prevent CPUID capture";
            return false;
        }
    }
    {
        FakeCpuidState cpuid = FrozenFakeCpuidState();
        FakeContextState context = FrozenFakeContextState();
        context.After.Cpu = 51;
        TargetIdentityReceiptV2 receipt;
        std::string rejected;
        if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
            context.Calls != 2u || cpuid.Calls != order || rejected.empty())
        {
            diagnostic = "invalid after-context did not fail closed";
            return false;
        }
    }
    return true;
}

bool TestCaptureLeafGates(std::string& diagnostic)
{
    const std::vector<std::pair<uint32_t, uint32_t> > full_order =
        FrozenCpuidCallOrder();
    struct GateCase
    {
        enum Mutation
        {
            BasicLeafB,
            ExtendedLeaf1e,
            TopologyExtensions,
            ExtendedLeaf26
        } Change;
        std::size_t ExpectedCalls;
    };
    const GateCase cases[] = {
        { GateCase::BasicLeafB, 8u },
        { GateCase::ExtendedLeaf1e, 6u },
        { GateCase::TopologyExtensions, 6u },
        { GateCase::ExtendedLeaf26, 11u }
    };
    for (std::size_t index = 0u;
         index < sizeof(cases) / sizeof(cases[0]);
         ++index)
    {
        FakeCpuidState cpuid = FrozenFakeCpuidState();
        if (cases[index].Change == GateCase::BasicLeafB) {
            cpuid.Values[std::make_pair(0u, 0u)].Eax = UINT32_C(0x0a);
        }
        else if (cases[index].Change == GateCase::ExtendedLeaf1e) {
            cpuid.Values[std::make_pair(UINT32_C(0x80000000), 0u)].Eax =
                UINT32_C(0x8000001d);
        }
        else if (cases[index].Change == GateCase::TopologyExtensions) {
            cpuid.Values[std::make_pair(UINT32_C(0x80000001), 0u)].Ecx &=
                ~(UINT32_C(1) << 22u);
        }
        else {
            cpuid.Values[std::make_pair(UINT32_C(0x80000000), 0u)].Eax =
                UINT32_C(0x80000025);
        }
        FakeContextState context = FrozenFakeContextState();
        TargetIdentityReceiptV2 receipt;
        std::string rejected;
        if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
            rejected.empty() || context.Calls != 2u ||
            cpuid.Calls.size() != cases[index].ExpectedCalls ||
            !std::equal(
                cpuid.Calls.begin(), cpuid.Calls.end(), full_order.begin()))
        {
            diagnostic = "CPUID maximum-leaf/feature gate queried reserved leaves";
            return false;
        }
    }
    return true;
}

bool TestTopologyEnumerationBounds(std::string& diagnostic)
{
    FakeCpuidState cpuid = FrozenFakeCpuidState();
    for (uint32_t subleaf = 0u; subleaf < kMaximumTopologyLevels; ++subleaf)
    {
        CpuidRegsV2 regs;
        regs.Eax = 5u;
        regs.Ebx = 2u;
        regs.Ecx = UINT32_C(0x00000100) | subleaf;
        regs.Edx = kFrozenFullApicId;
        cpuid.Values[std::make_pair(UINT32_C(0x0b), subleaf)] = regs;
    }
    FakeContextState context = FrozenFakeContextState();
    TargetIdentityReceiptV2 receipt;
    std::string rejected;
    if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
        rejected.find("bounded sentinel") == std::string::npos ||
        cpuid.Calls.empty() ||
        cpuid.Calls.back() !=
            std::make_pair(UINT32_C(0x0b), kMaximumTopologyLevels - 1u))
    {
        diagnostic = "topology sentinel bound was not enforced";
        return false;
    }
    cpuid = FrozenFakeCpuidState();
    cpuid.Values[std::make_pair(UINT32_C(0x0b), 2u)].Ecx =
        UINT32_C(0x00000102);
    context = FrozenFakeContextState();
    if (CaptureFrozenFake(receipt, cpuid, context, rejected))
    {
        diagnostic = "half-valid topology sentinel was accepted";
        return false;
    }
    cpuid = FrozenFakeCpuidState();
    cpuid.Values[std::make_pair(UINT32_C(0x0b), 2u)].Ebx = 1u;
    cpuid.Values[std::make_pair(UINT32_C(0x0b), 2u)].Ecx =
        UINT32_C(0x00000002);
    context = FrozenFakeContextState();
    if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
        cpuid.Calls.empty() ||
        std::find(
            cpuid.Calls.begin(),
            cpuid.Calls.end(),
            std::make_pair(UINT32_C(0x0b), 3u)) != cpuid.Calls.end())
    {
        diagnostic =
            "zero-type/nonzero-count sentinel queried reserved subleaves";
        return false;
    }
    cpuid = FrozenFakeCpuidState();
    cpuid.Values[std::make_pair(UINT32_C(0x80000026), 4u)].Ebx = 1u;
    cpuid.Values[std::make_pair(UINT32_C(0x80000026), 4u)].Ecx =
        UINT32_C(0x00000004);
    context = FrozenFakeContextState();
    if (CaptureFrozenFake(receipt, cpuid, context, rejected) ||
        std::find(
            cpuid.Calls.begin(),
            cpuid.Calls.end(),
            std::make_pair(UINT32_C(0x80000026), 5u)) != cpuid.Calls.end())
    {
        diagnostic =
            "80000026 zero-type sentinel queried a reserved subleaf";
        return false;
    }
    return true;
}

class GroupEveryDigitNumpunct : public std::numpunct<char>
{
protected:
    char do_thousands_sep() const override { return '_'; }
    std::string do_grouping() const override { return std::string(1u, '\1'); }
};

std::string ProcCpuInfoText(
    const std::vector<CpuInfoMapRecordV2>& records,
    bool reverse)
{
    std::ostringstream text;
    text.imbue(std::locale::classic());
    for (std::size_t offset = 0u; offset < records.size(); ++offset)
    {
        const std::size_t index = reverse ? records.size() - 1u - offset : offset;
        const CpuInfoMapRecordV2& record = records[index];
        text << "processor\t: " << record.Cpu << '\n'
             << "vendor_id\t: AuthenticAMD\n"
             << "physical id\t: " << record.Package << '\n'
             << "core id\t\t: " << record.Core << '\n'
             << "apicid\t\t: " << record.ApicId << '\n'
             << "initial apicid\t: " << record.InitialApicId << "\n\n";
    }
    return text.str();
}

bool TestCpuInfoMapAndParser(std::string& diagnostic)
{
    const std::vector<CpuInfoMapRecordV2> frozen = FrozenCpuInfoMapV2();
    if (!ValidateFrozenCpuInfoMapV2(frozen, diagnostic)) return false;
    std::string canonical;
    if (!CanonicalizeCpuInfoMapV2(frozen, canonical, diagnostic) ||
        canonical.size() != FrozenCpuInfoMapCanonicalBytesV2() ||
        wirehair::wh2_benchmark::Sha256Hex(canonical) !=
            FrozenCpuInfoMapSha256V2()) return false;
    {
        const TargetIdentityReceiptV2 identity =
            FrozenCpu50TargetIdentityReceiptV2();
        const std::string baseline_evidence = FormatTargetIdentityEvidenceV2(
            "locale-checkpoint", identity, true, "locale-diagnostic");
        const std::locale previous = std::locale();
        const std::locale grouped(
            std::locale::classic(), new GroupEveryDigitNumpunct());
        std::locale::global(grouped);
        std::string localized;
        std::string localized_diagnostic;
        const bool locale_stable =
            CanonicalizeCpuInfoMapV2(
                frozen, localized, localized_diagnostic) &&
            localized == canonical &&
            CpuInfoMapSha256V2(frozen, localized_diagnostic) ==
                FrozenCpuInfoMapSha256V2() &&
            FormatTargetIdentityEvidenceV2(
                "locale-checkpoint",
                identity,
                true,
                "locale-diagnostic") == baseline_evidence;
        std::locale::global(previous);
        if (!locale_stable)
        {
            diagnostic = "canonical/evidence formatting depends on global locale";
            return false;
        }
    }
    std::vector<CpuInfoMapRecordV2> parsed;
    const std::string reversed_text = ProcCpuInfoText(frozen, true);
    if (!ParseProcCpuInfoV2(reversed_text, parsed, diagnostic) ||
        !ValidateFrozenCpuInfoMapV2(parsed, diagnostic) ||
        CpuInfoMapSha256V2(parsed, diagnostic) != FrozenCpuInfoMapSha256V2()) {
        return false;
    }
    std::vector<CpuInfoMapRecordV2> reversed(frozen.rbegin(), frozen.rend());
    std::string reversed_canonical;
    if (!CanonicalizeCpuInfoMapV2(reversed, reversed_canonical, diagnostic) ||
        reversed_canonical != canonical) return false;

    std::string malformed = reversed_text;
    const std::size_t physical = malformed.find("physical id");
    const std::size_t physical_end = malformed.find('\n', physical);
    malformed.erase(physical, physical_end + 1u - physical);
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed = reversed_text;
    malformed.replace(malformed.find("processor\t: 127"), 15u,
        "processor\t: +127");
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed = reversed_text;
    malformed.replace(malformed.find("processor\t: 127"), 15u,
        "processor\t: 4294967296");
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed = reversed_text;
    malformed.insert(malformed.find('\n') + 1u, "processor: 127\n");
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed = reversed_text;
    malformed.insert(malformed.size() / 2u, 1u, '\0');
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed = reversed_text;
    malformed.insert(malformed.size() / 2u, 1u, '\r');
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed = "processor: 0\nphysical id: 0\ncore id: 0\n"
        "apicid: 1\ninitial apicid: 1\n\n"
        "processor: 1\nphysical id: 0\ncore id: 1\n"
        "apicid: 1\ninitial apicid: 1\n\n";
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed.assign(kMaximumCpuInfoLineBytes + 1u, 'x');
    malformed += ": 0\n";
    if (ParseProcCpuInfoV2(malformed, parsed, diagnostic)) goto failure;
    malformed.assign(kMaximumCpuInfoLineBytes - 1u, 'x');
    malformed += ":\nprocessor: 0\nphysical id: 0\ncore id: 0\n"
        "apicid: 0\ninitial apicid: 0\n";
    if (!ParseProcCpuInfoV2(malformed, parsed, diagnostic) ||
        parsed.size() != 1u) {
        goto failure;
    }
    {
        std::vector<CpuInfoMapRecordV2> too_many = frozen;
        CpuInfoMapRecordV2 extra;
        extra.Cpu = 128u;
        extra.Package = 1u;
        extra.Core = 0u;
        extra.ApicId = 128u;
        extra.InitialApicId = 128u;
        too_many.push_back(extra);
        if (ParseProcCpuInfoV2(
                ProcCpuInfoText(too_many, false), parsed, diagnostic)) {
            goto failure;
        }
    }
    {
        std::vector<CpuInfoMapRecordV2> duplicate = frozen;
        duplicate[1].Cpu = duplicate[0].Cpu;
        if (CanonicalizeCpuInfoMapV2(duplicate, canonical, diagnostic)) {
            goto failure;
        }
    }
    {
        std::vector<CpuInfoMapRecordV2> duplicate = frozen;
        duplicate[1].ApicId = duplicate[0].ApicId;
        if (CanonicalizeCpuInfoMapV2(duplicate, canonical, diagnostic)) {
            goto failure;
        }
        duplicate = frozen;
        duplicate[1].InitialApicId = duplicate[0].InitialApicId;
        if (CanonicalizeCpuInfoMapV2(duplicate, canonical, diagnostic)) {
            goto failure;
        }
    }
    {
        malformed = ProcCpuInfoText(frozen, false);
        const std::size_t second = malformed.find("processor\t: 1");
        malformed.replace(second, 9u, "brokenxx");
        parsed.push_back(frozen[0]);
        if (ParseProcCpuInfoV2(malformed, parsed, diagnostic) || !parsed.empty()) {
            goto failure;
        }
    }
    for (unsigned field = 0u; field < 5u; ++field)
    {
        std::vector<CpuInfoMapRecordV2> changed = frozen;
        uint32_t* const fields[] = {
            &changed[50].Cpu,
            &changed[50].Package,
            &changed[50].Core,
            &changed[50].ApicId,
            &changed[50].InitialApicId
        };
        *fields[field] ^= 1u;
        if (ValidateFrozenCpuInfoMapV2(changed, diagnostic)) goto failure;
    }
    diagnostic.clear();
    return true;

failure:
    diagnostic = "malformed CPU/APIC map input was accepted";
    return false;
}

bool TestEvidenceEncoding(std::string& diagnostic)
{
    TargetIdentityReceiptV2 receipt =
        FrozenCpu50TargetIdentityReceiptV2();
    receipt.CanonicalSha256 = std::string("sha,\n\0", 6u);
    const std::string checkpoint("checkpoint,\n\0", 13u);
    const std::string detail("detail,\n\0", 9u);
    const std::string evidence = FormatTargetIdentityEvidenceV2(
        checkpoint, receipt, false, detail);
    if (evidence.find("checkpoint_hex=") == std::string::npos ||
        evidence.find("diagnostic_hex=") == std::string::npos ||
        evidence.find("canonical_sha256_hex=") == std::string::npos ||
        evidence.find("gate=FAIL") == std::string::npos ||
        evidence.find("leaf=0x0000000b") == std::string::npos ||
        evidence.find("subleaf=0x00000004") == std::string::npos ||
        evidence.find(checkpoint) != std::string::npos ||
        evidence.find(detail) != std::string::npos ||
        evidence.find(receipt.CanonicalSha256) != std::string::npos)
    {
        diagnostic = "target identity evidence is incomplete or unsafely encoded";
        return false;
    }
    return true;
}

} // namespace

bool RunTargetIdentityV2SelfTests(std::string& diagnostic)
{
    struct Test
    {
        const char* Name;
        bool (*Run)(std::string&);
    };
    const Test tests[] = {
        { "frozen-baseline-and-normalized-hash", TestFrozenIdentityBaseline },
        { "full-width-apic-semantics", TestFullWidthApicSemantics },
        { "frozen-raw-bit-mutations", TestFrozenRawMutationSurface },
        { "semantic-topology-mutations", TestSemanticMutationSurface },
        { "context-mutations", TestContextMutationSurface },
        { "injected-capture-failures", TestInjectedCaptureFailures },
        { "capture-leaf-feature-gates", TestCaptureLeafGates },
        { "topology-first-sentinel-bound", TestTopologyEnumerationBounds },
        { "cpuinfo-map-parser-and-hash", TestCpuInfoMapAndParser },
        { "evidence-binary-encoding", TestEvidenceEncoding }
    };
    for (std::size_t i = 0u; i < sizeof(tests) / sizeof(tests[0]); ++i)
    {
        std::string detail;
        if (!tests[i].Run(detail))
        {
            diagnostic = std::string(tests[i].Name) + ": " +
                (detail.empty() ? "failed without diagnostic" : detail);
            return false;
        }
    }
    diagnostic.clear();
    return true;
}

} // namespace wirehair_wh2_bench
