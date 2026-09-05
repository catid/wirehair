#include "Wh2PublicBorrowedTargetIdentity.h"

#include <cerrno>
#include <cstring>

#if defined(__linux__)
#include <sched.h>
#include <sys/resource.h>
#endif

#if defined(__linux__) && defined(__x86_64__) && \
    (defined(__GNUC__) || defined(__clang__))
#include <cpuid.h>
#define WIREHAIR_PUBLIC_BORROWED_TARGET_IDENTITY_NATIVE 1
#else
#define WIREHAIR_PUBLIC_BORROWED_TARGET_IDENTITY_NATIVE 0
#endif

namespace wirehair_wh2_bench {
namespace {

#if WIREHAIR_PUBLIC_BORROWED_TARGET_IDENTITY_NATIVE
void CompilerBarrier()
{
    __asm__ __volatile__("" ::: "memory");
}

bool ReadCpuid(
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
    CompilerBarrier();
    __cpuid_count(leaf, subleaf, eax, ebx, ecx, edx);
    CompilerBarrier();
    registers.Eax = eax;
    registers.Ebx = ebx;
    registers.Ecx = ecx;
    registers.Edx = edx;
    diagnostic.clear();
    return true;
}

bool CaptureContext(
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

bool CapturePublicBorrowedTargetIdentity(
    int32_t target_cpu,
    TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
#if WIREHAIR_PUBLIC_BORROWED_TARGET_IDENTITY_NATIVE
    CompilerBarrier();
    const bool captured = CaptureTargetIdentityV2(
        target_cpu,
        ReadCpuid,
        nullptr,
        CaptureContext,
        nullptr,
        receipt,
        diagnostic);
    CompilerBarrier();
    return captured;
#else
    receipt = TargetIdentityReceiptV2();
    receipt.RequestedCpu = target_cpu;
    diagnostic =
        "target identity requires Linux x86_64 GNU-compatible CPUID";
    return false;
#endif
}

} // namespace wirehair_wh2_bench
