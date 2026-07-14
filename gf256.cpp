/** \file
    \brief GF(256) Main C API Source
    \copyright Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of GF256 nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include "gf256.h"
#include "WirehairEnvironment.h"

#include <cstdlib>
#include <cstring>
#include <mutex>

#ifdef WH_COUNT
// Task 6a: thread-local per-op byte/call counters.
static thread_local uint64_t t_gf_bytes[6] = {0,0,0,0,0,0};
static thread_local uint64_t t_gf_calls[6] = {0,0,0,0,0,0};
extern "C" uint64_t gf256_count_bytes(int op) { return (op>=0&&op<6) ? t_gf_bytes[op] : 0; }
extern "C" uint64_t gf256_count_calls(int op) { return (op>=0&&op<6) ? t_gf_calls[op] : 0; }
extern "C" void gf256_count_reset(void) { for (int i=0;i<6;++i){t_gf_bytes[i]=0;t_gf_calls[i]=0;} }
#define WH_BUMP(op, n) do { t_gf_bytes[op] += (uint64_t)(n); t_gf_calls[op]++; } while(0)
#else
#define WH_BUMP(op, n) do {} while(0)
#endif

static GF256_FORCE_INLINE uint64_t gf256_loadu64(const void* p)
{
    uint64_t v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

static GF256_FORCE_INLINE uint32_t gf256_loadu32(const void* p)
{
    uint32_t v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

static GF256_FORCE_INLINE void gf256_storeu64(void* p, uint64_t v)
{
    std::memcpy(p, &v, sizeof(v));
}

static GF256_FORCE_INLINE void gf256_storeu32(void* p, uint32_t v)
{
    std::memcpy(p, &v, sizeof(v));
}

#if __cplusplus >= 201703L && defined(__has_cpp_attribute)
# if __has_cpp_attribute(fallthrough)
#  define GF256_FALLTHROUGH [[fallthrough]]
# endif
#endif
#ifndef GF256_FALLTHROUGH
# if defined(__GNUC__) && __GNUC__ >= 7
#  define GF256_FALLTHROUGH __attribute__((fallthrough))
# else
#  define GF256_FALLTHROUGH ((void)0)
# endif
#endif

#ifdef LINUX_ARM
#include <unistd.h>
#include <fcntl.h>
#include <elf.h>
#include <linux/auxvec.h>
#endif

//------------------------------------------------------------------------------
// Detect host byte order.
// This check works with GCC and LLVM; assume little-endian byte order when
// using any other compiler.
// The result is verified during initialization.
//
#if defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) \
    && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define GF256_IS_BIG_ENDIAN
#endif

//------------------------------------------------------------------------------
// Workaround for ARMv7 that doesn't provide vqtbl1_*
// This comes from linux-raid (https://www.spinics.net/lists/raid/msg58403.html)
//
#ifdef GF256_TRY_NEON
#if __ARM_ARCH <= 7 && !defined(__aarch64__)
static GF256_FORCE_INLINE uint8x16_t vqtbl1q_u8(uint8x16_t a, uint8x16_t b)
{
    union {
        uint8x16_t    val;
        uint8x8x2_t    pair;
    } __a = { a };

    return vcombine_u8(vtbl2_u8(__a.pair, vget_low_u8(b)),
                       vtbl2_u8(__a.pair, vget_high_u8(b)));
}
#endif
#endif

//------------------------------------------------------------------------------
// Self-Test
//
// This is executed during initialization to make sure the library is working

// Cover at least one 64-byte iteration so the GFNI and AVX-512 kernels
// (which only engage for bytes >= 64) are exercised, plus an odd tail.
static const unsigned kTestBufferBytes = 64 + 32 + 16 + 8 + 4 + 2 + 1;
static const unsigned kTestBufferAllocated = 128;
struct SelfTestBuffersT
{
    GF256_ALIGNED uint8_t A[kTestBufferAllocated];
    GF256_ALIGNED uint8_t B[kTestBufferAllocated];
    GF256_ALIGNED uint8_t C[kTestBufferAllocated];
    GF256_ALIGNED uint8_t D[kTestBufferAllocated];
    GF256_ALIGNED uint8_t E[kTestBufferAllocated];
    GF256_ALIGNED uint8_t F[kTestBufferAllocated];
};
static GF256_ALIGNED SelfTestBuffersT m_SelfTestBuffers;

static bool gf256_self_test()
{
    if ((uintptr_t)m_SelfTestBuffers.A % GF256_ALIGN_BYTES != 0)
        return false;
    if ((uintptr_t)m_SelfTestBuffers.A % GF256_ALIGN_BYTES != 0)
        return false;
    if ((uintptr_t)m_SelfTestBuffers.B % GF256_ALIGN_BYTES != 0)
        return false;
    if ((uintptr_t)m_SelfTestBuffers.C % GF256_ALIGN_BYTES != 0)
        return false;
    if ((uintptr_t)m_SelfTestBuffers.D % GF256_ALIGN_BYTES != 0)
        return false;
    if ((uintptr_t)m_SelfTestBuffers.E % GF256_ALIGN_BYTES != 0)
        return false;
    if ((uintptr_t)m_SelfTestBuffers.F % GF256_ALIGN_BYTES != 0)
        return false;

    // Check multiplication/division
    for (unsigned i = 0; i < 256; ++i)
    {
        for (unsigned j = 0; j < 256; ++j)
        {
            uint8_t prod = gf256_mul((uint8_t)i, (uint8_t)j);
            if (i != 0 && j != 0)
            {
                uint8_t div1 = gf256_div(prod, (uint8_t)i);
                if (div1 != j)
                    return false;
                uint8_t div2 = gf256_div(prod, (uint8_t)j);
                if (div2 != i)
                    return false;
            }
            else if (prod != 0)
                return false;
            if (j == 1 && prod != i)
                return false;
        }
    }

    // Check for overruns
    m_SelfTestBuffers.A[kTestBufferBytes] = 0x5a;
    m_SelfTestBuffers.B[kTestBufferBytes] = 0x5a;
    m_SelfTestBuffers.C[kTestBufferBytes] = 0x5a;
    m_SelfTestBuffers.D[kTestBufferBytes] = 0x5a;
    m_SelfTestBuffers.E[kTestBufferBytes] = 0x5a;
    m_SelfTestBuffers.F[kTestBufferBytes] = 0x5a;

    // Test gf256_add_mem()
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        m_SelfTestBuffers.A[i] = 0x1f;
        m_SelfTestBuffers.B[i] = 0xf7;
    }
    gf256_add_mem(m_SelfTestBuffers.A, m_SelfTestBuffers.B, kTestBufferBytes);
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
        if (m_SelfTestBuffers.A[i] != (0x1f ^ 0xf7))
            return false;

    // Test gf256_add2_mem()
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        m_SelfTestBuffers.A[i] = 0x1f;
        m_SelfTestBuffers.B[i] = 0xf7;
        m_SelfTestBuffers.C[i] = 0x71;
    }
    gf256_add2_mem(m_SelfTestBuffers.A, m_SelfTestBuffers.B, m_SelfTestBuffers.C, kTestBufferBytes);
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
        if (m_SelfTestBuffers.A[i] != (0x1f ^ 0xf7 ^ 0x71))
            return false;

    // Test gf256_addset_mem()
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        m_SelfTestBuffers.A[i] = 0x55;
        m_SelfTestBuffers.B[i] = 0xaa;
        m_SelfTestBuffers.C[i] = 0x6c;
    }
    gf256_addset_mem(m_SelfTestBuffers.A, m_SelfTestBuffers.B, m_SelfTestBuffers.C, kTestBufferBytes);
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
        if (m_SelfTestBuffers.A[i] != (0xaa ^ 0x6c))
            return false;

    // Test gf256_add_multi_mem()
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        m_SelfTestBuffers.A[i] = (uint8_t)(i * 3u + 0x11u);
        m_SelfTestBuffers.B[i] = (uint8_t)(i * 5u + 0x23u);
        m_SelfTestBuffers.C[i] = (uint8_t)(i * 7u + 0x35u);
        m_SelfTestBuffers.D[i] = (uint8_t)(i * 11u + 0x47u);
        m_SelfTestBuffers.E[i] = (uint8_t)(i * 13u + 0x59u);
        m_SelfTestBuffers.F[i] = (uint8_t)(i * 17u + 0x6bu);
    }
    {
        const void* srcs[] = {
            m_SelfTestBuffers.B,
            m_SelfTestBuffers.C,
            m_SelfTestBuffers.D,
            m_SelfTestBuffers.E,
            m_SelfTestBuffers.F
        };
        gf256_add_multi_mem(m_SelfTestBuffers.A, srcs, 5, kTestBufferBytes);
    }
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        const uint8_t expected =
            (uint8_t)(i * 3u + 0x11u) ^
            (uint8_t)(i * 5u + 0x23u) ^
            (uint8_t)(i * 7u + 0x35u) ^
            (uint8_t)(i * 11u + 0x47u) ^
            (uint8_t)(i * 13u + 0x59u) ^
            (uint8_t)(i * 17u + 0x6bu);
        if (m_SelfTestBuffers.A[i] != expected)
            return false;
    }

    // Test gf256_addset_multi_mem()
    for (unsigned i = 0; i < kTestBufferBytes; ++i) {
        m_SelfTestBuffers.A[i] = (uint8_t)(i * 19u + 0x7du);
    }
    {
        const void* srcs[] = {
            m_SelfTestBuffers.B,
            m_SelfTestBuffers.C,
            m_SelfTestBuffers.D,
            m_SelfTestBuffers.E,
            m_SelfTestBuffers.F
        };
        gf256_addset_multi_mem(
            m_SelfTestBuffers.A, srcs, 5, kTestBufferBytes);
    }
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        const uint8_t expected =
            (uint8_t)(i * 5u + 0x23u) ^
            (uint8_t)(i * 7u + 0x35u) ^
            (uint8_t)(i * 11u + 0x47u) ^
            (uint8_t)(i * 13u + 0x59u) ^
            (uint8_t)(i * 17u + 0x6bu);
        if (m_SelfTestBuffers.A[i] != expected)
            return false;
    }

    // Test gf256_muladd_mem()
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        m_SelfTestBuffers.A[i] = 0xff;
        m_SelfTestBuffers.B[i] = 0xaa;
    }
    const uint8_t expectedMulAdd = gf256_mul(0xaa, 0x6c);
    gf256_muladd_mem(m_SelfTestBuffers.A, 0x6c, m_SelfTestBuffers.B, kTestBufferBytes);
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
        if (m_SelfTestBuffers.A[i] != (expectedMulAdd ^ 0xff))
            return false;

    // Test gf256_muladd_multi_mem(), including the zero/identity scales and
    // the odd-sized SIMD tail used by this self-test.
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        m_SelfTestBuffers.A[i] = (uint8_t)(i * 3u + 7u);
        m_SelfTestBuffers.B[i] = (uint8_t)(i * 5u + 11u);
        m_SelfTestBuffers.C[i] = (uint8_t)(i * 7u + 13u);
        m_SelfTestBuffers.D[i] = (uint8_t)(i * 11u + 17u);
        m_SelfTestBuffers.F[i] = (uint8_t)(i * 13u + 19u);
    }
    {
        void* destinations[] = {
            m_SelfTestBuffers.A,
            m_SelfTestBuffers.B,
            m_SelfTestBuffers.C,
            m_SelfTestBuffers.D
        };
        const uint8_t scales[] = { 0u, 1u, 2u, 0x6cu };
        gf256_muladd_multi_mem(
            destinations, scales, 4, m_SelfTestBuffers.F, kTestBufferBytes);
    }
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        const uint8_t source = (uint8_t)(i * 13u + 19u);
        if (m_SelfTestBuffers.A[i] != (uint8_t)(i * 3u + 7u) ||
            m_SelfTestBuffers.B[i] !=
                ((uint8_t)(i * 5u + 11u) ^ source) ||
            m_SelfTestBuffers.C[i] !=
                ((uint8_t)(i * 7u + 13u) ^ gf256_mul(source, 2u)) ||
            m_SelfTestBuffers.D[i] !=
                ((uint8_t)(i * 11u + 17u) ^ gf256_mul(source, 0x6cu)))
        {
            return false;
        }
    }

    // Test gf256_mul_mem()
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
    {
        m_SelfTestBuffers.A[i] = 0xff;
        m_SelfTestBuffers.B[i] = 0x55;
    }
    const uint8_t expectedMul = gf256_mul(0xa2, 0x55);
    gf256_mul_mem(m_SelfTestBuffers.A, m_SelfTestBuffers.B, 0xa2, kTestBufferBytes);
    for (unsigned i = 0; i < kTestBufferBytes; ++i)
        if (m_SelfTestBuffers.A[i] != expectedMul)
            return false;

    if (m_SelfTestBuffers.A[kTestBufferBytes] != 0x5a)
        return false;
    if (m_SelfTestBuffers.B[kTestBufferBytes] != 0x5a)
        return false;
    if (m_SelfTestBuffers.C[kTestBufferBytes] != 0x5a)
        return false;
    if (m_SelfTestBuffers.D[kTestBufferBytes] != 0x5a)
        return false;
    if (m_SelfTestBuffers.E[kTestBufferBytes] != 0x5a)
        return false;
    if (m_SelfTestBuffers.F[kTestBufferBytes] != 0x5a)
        return false;

    return true;
}


//------------------------------------------------------------------------------
// Runtime CPU Architecture Check
//
// Feature checks stolen shamelessly from
// https://github.com/jedisct1/libsodium/blob/master/src/libsodium/sodium/runtime.c

#if defined(HAVE_ANDROID_GETCPUFEATURES)
#include <cpu-features.h>
#endif

#if defined(GF256_TRY_NEON)
# if defined(IOS) && (defined(__ARM_NEON) || defined(__ARM_NEON__))
// Requires iPhone 5S or newer
static const bool CpuHasNeon = true;
# else // ANDROID or LINUX_ARM
#  if defined(__aarch64__)
static bool CpuHasNeon = true;      // if AARCH64, then we have NEON for sure...
#  else
static bool CpuHasNeon = false;     // if not, then we have to check at runtime.
#  endif
# endif
#endif

#if !defined(GF256_TARGET_MOBILE)

#ifdef _MSC_VER
    #include <intrin.h> // __cpuid
    #pragma warning(disable: 4752) // found Intel(R) Advanced Vector Extensions; consider using /arch:AVX
#endif

#if defined(GF256_TRY_AVX2) || defined(GF256_TRY_TARGET_AVX2)
static bool CpuHasAVX2 = false;
#endif
#ifdef GF256_TRY_GFNI
static bool CpuHasGFNI = false;
#endif
#ifdef GF256_TRY_AVX512
static bool CpuHasAVX512 = false;
#endif
#if defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3)
static bool CpuHasSSSE3 = false;
#endif

#define CPUID_EBX_AVX2    0x00000020
#if defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3)
#define CPUID_ECX_SSSE3   0x00000200
#endif

static void _cpuid(unsigned int cpu_info[4U], const unsigned int cpu_info_type)
{
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    __cpuidex((int *) cpu_info, cpu_info_type, 0);
#else //if defined(HAVE_CPUID)
    cpu_info[0] = cpu_info[1] = cpu_info[2] = cpu_info[3] = 0;
# ifdef __i386__
    __asm__ __volatile__ ("pushfl; pushfl; "
                          "popl %0; "
                          "movl %0, %1; xorl %2, %0; "
                          "pushl %0; "
                          "popfl; pushfl; popl %0; popfl" :
                          "=&r" (cpu_info[0]), "=&r" (cpu_info[1]) :
                          "i" (0x200000));
    if (((cpu_info[0] ^ cpu_info[1]) & 0x200000) == 0) {
        return; /* LCOV_EXCL_LINE */
    }
# endif
# ifdef __i386__
    __asm__ __volatile__ ("xchgl %%ebx, %k1; cpuid; xchgl %%ebx, %k1" :
                          "=a" (cpu_info[0]), "=&r" (cpu_info[1]),
                          "=c" (cpu_info[2]), "=d" (cpu_info[3]) :
                          "0" (cpu_info_type), "2" (0U));
# elif defined(__x86_64__)
    __asm__ __volatile__ ("xchgq %%rbx, %q1; cpuid; xchgq %%rbx, %q1" :
                          "=a" (cpu_info[0]), "=&r" (cpu_info[1]),
                          "=c" (cpu_info[2]), "=d" (cpu_info[3]) :
                          "0" (cpu_info_type), "2" (0U));
# else
    __asm__ __volatile__ ("cpuid" :
                          "=a" (cpu_info[0]), "=b" (cpu_info[1]),
                          "=c" (cpu_info[2]), "=d" (cpu_info[3]) :
                          "0" (cpu_info_type), "2" (0U));
# endif
#endif
}

static uint64_t _xgetbv0()
{
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    return _xgetbv(0);
#elif defined(__i386__) || defined(__x86_64__)
    uint32_t eax, edx;
    __asm__ __volatile__ (".byte 0x0f, 0x01, 0xd0" :
                          "=a" (eax), "=d" (edx) : "c" (0));
    return (static_cast<uint64_t>(edx) << 32) | eax;
#else
    return 0;
#endif
}

#else
#if defined(LINUX_ARM)
static void checkLinuxARMNeonCapabilities( bool& cpuHasNeon )
{
#if defined(__aarch64__)
    // AArch64 mandates Advanced SIMD (NEON).  Parsing 64-bit auxv records
    // with the 32-bit Elf32_auxv_t layout would misread AT_HWCAP and
    // wrongly disable NEON here.
    cpuHasNeon = true;
#else
    auto cpufile = open("/proc/self/auxv", O_RDONLY);
    Elf32_auxv_t auxv;
    if (cpufile >= 0)
    {
        const auto size_auxv_t = sizeof(Elf32_auxv_t);
        while (read(cpufile, &auxv, size_auxv_t) == size_auxv_t)
        {
            if (auxv.a_type == AT_HWCAP)
            {
                cpuHasNeon = (auxv.a_un.a_val & 4096) != 0;
                break;
            }
        }
        close(cpufile);
    }
    else
    {
        cpuHasNeon = false;
    }
#endif
}
#endif
#endif // defined(GF256_TARGET_MOBILE)

extern "C" void gf256_select_x86_cpu_features(
    const gf256_x86_cpu_snapshot* snapshot,
    gf256_x86_cpu_features* features)
{
    if (!features) {
        return;
    }
    std::memset(features, 0, sizeof(*features));
    if (!snapshot) {
        return;
    }

    const bool has_leaf1 = snapshot->MaxBasicLeaf >= 1;
    const bool has_leaf7 = snapshot->MaxBasicLeaf >= 7;
    const bool osxsave = has_leaf1 && (snapshot->Leaf1ECX & (1u << 27)) != 0;
    const bool cpu_avx = has_leaf1 && (snapshot->Leaf1ECX & (1u << 28)) != 0;
    const bool ymm_state = osxsave && cpu_avx &&
        (snapshot->XCR0 & UINT64_C(0x6)) == UINT64_C(0x6);
    const bool zmm_state = ymm_state &&
        (snapshot->XCR0 & UINT64_C(0xe6)) == UINT64_C(0xe6);

    features->SSSE3 = has_leaf1 &&
        (snapshot->Leaf1ECX & (1u << 9)) != 0;
    features->AVX2 = has_leaf7 && ymm_state &&
        (snapshot->Leaf7EBX & (1u << 5)) != 0;
    features->AVX512 = has_leaf7 && zmm_state &&
        (snapshot->Leaf7EBX & (1u << 16)) != 0;
    features->GFNI = has_leaf7 && zmm_state &&
        (snapshot->Leaf7EBX & (1u << 16)) != 0 &&
        (snapshot->Leaf7EBX & (1u << 30)) != 0 &&
        (snapshot->Leaf7ECX & (1u << 8)) != 0;
}

extern "C" void gf256_get_active_x86_cpu_features(
    gf256_x86_cpu_features* features)
{
    if (!features) {
        return;
    }
    std::memset(features, 0, sizeof(*features));
    if (gf256_init() != 0) {
        return;
    }
#if !defined(GF256_TARGET_MOBILE)
# if defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3)
    features->SSSE3 = CpuHasSSSE3;
# endif
# if defined(GF256_TRY_AVX2) || defined(GF256_TRY_TARGET_AVX2)
    features->AVX2 = CpuHasAVX2;
# endif
# if defined(GF256_TRY_GFNI)
    features->GFNI = CpuHasGFNI;
# endif
# if defined(GF256_TRY_AVX512)
    features->AVX512 = CpuHasAVX512;
# endif
#endif
}

static void gf256_architecture_init()
{
#if defined(GF256_TRY_NEON)

    // Check for NEON support on Android platform
#if defined(HAVE_ANDROID_GETCPUFEATURES)
    AndroidCpuFamily family = android_getCpuFamily();
    if (family == ANDROID_CPU_FAMILY_ARM)
    {
        if (android_getCpuFeatures() & ANDROID_CPU_ARM_FEATURE_NEON)
            CpuHasNeon = true;
    }
    else if (family == ANDROID_CPU_FAMILY_ARM64)
    {
        CpuHasNeon = true;
    }
#endif

#if defined(LINUX_ARM)
    // Check for NEON support on other ARM/Linux platforms
    checkLinuxARMNeonCapabilities(CpuHasNeon);
#endif

#endif //GF256_TRY_NEON

#if !defined(GF256_TARGET_MOBILE)
    unsigned int cpu_info[4] = { 0, 0, 0, 0 };
    gf256_x86_cpu_snapshot snapshot = {};
    _cpuid(cpu_info, 0);
    snapshot.MaxBasicLeaf = cpu_info[0];
    if (snapshot.MaxBasicLeaf >= 1)
    {
        _cpuid(cpu_info, 1);
        snapshot.Leaf1ECX = cpu_info[2];
        const uint32_t xsave_mask = (1u << 27) | (1u << 28);
        if ((snapshot.Leaf1ECX & xsave_mask) == xsave_mask) {
            snapshot.XCR0 = _xgetbv0();
        }
    }
    if (snapshot.MaxBasicLeaf >= 7)
    {
        _cpuid(cpu_info, 7);
        snapshot.Leaf7EBX = cpu_info[1];
        snapshot.Leaf7ECX = cpu_info[2];
    }

    gf256_x86_cpu_features features;
    gf256_select_x86_cpu_features(&snapshot, &features);
#if defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3)
    CpuHasSSSE3 = features.SSSE3 != 0;
#endif

#if defined(GF256_TRY_AVX2) || defined(GF256_TRY_TARGET_AVX2)
    CpuHasAVX2 = features.AVX2 != 0;
#endif // GF256_TRY_AVX2 || GF256_TRY_TARGET_AVX2

#if defined(GF256_TRY_GFNI)
    CpuHasGFNI = features.GFNI != 0;
#endif // GF256_TRY_GFNI

#if defined(GF256_TRY_AVX512)
    CpuHasAVX512 = features.AVX512 != 0;
#endif // GF256_TRY_AVX512

    // When AVX2 and SSSE3 are unavailable at compile or run time, Siamese
    // takes 4x longer to decode and 2.6x longer to encode.  Encoding requires
    // a lot more simple XOR ops so it is still pretty fast.  Decoding is
    // usually really quick because average loss rates are low, but when needed
    // it requires a lot more GF multiplies requiring table lookups.

#endif // GF256_TARGET_MOBILE
}


//------------------------------------------------------------------------------
// Context Object

// Context object for GF(2^^8) math
GF256_ALIGNED gf256_ctx GF256Ctx;
static std::once_flag InitOnce;
static int InitResult = -3;


//------------------------------------------------------------------------------
// Generator Polynomial

// There are only 16 irreducible polynomials for GF(2^^8)
static const int GF256_GEN_POLY_COUNT = 16;
static const uint8_t GF256_GEN_POLY[GF256_GEN_POLY_COUNT] = {
    0x8e, 0x95, 0x96, 0xa6, 0xaf, 0xb1, 0xb2, 0xb4,
    0xb8, 0xc3, 0xc6, 0xd4, 0xe1, 0xe7, 0xf3, 0xfa
};

static const int kDefaultPolynomialIndex = 3;

// Select which polynomial to use
static void gf256_poly_init(int polynomialIndex)
{
    if (polynomialIndex < 0 || polynomialIndex >= GF256_GEN_POLY_COUNT)
        polynomialIndex = kDefaultPolynomialIndex;

    GF256Ctx.Polynomial = (GF256_GEN_POLY[polynomialIndex] << 1) | 1;
}


//------------------------------------------------------------------------------
// Exponential and Log Tables

// Construct EXP and LOG tables from polynomial
static void gf256_explog_init()
{
    unsigned poly = GF256Ctx.Polynomial;
    uint8_t* exptab = GF256Ctx.GF256_EXP_TABLE;
    uint16_t* logtab = GF256Ctx.GF256_LOG_TABLE;

    logtab[0] = 512;
    exptab[0] = 1;
    for (unsigned jj = 1; jj < 255; ++jj)
    {
        unsigned next = (unsigned)exptab[jj - 1] * 2;
        if (next >= 256)
            next ^= poly;

        exptab[jj] = static_cast<uint8_t>( next );
        logtab[exptab[jj]] = static_cast<uint16_t>( jj );
    }
    exptab[255] = exptab[0];
    logtab[exptab[255]] = 255;
    for (unsigned jj = 256; jj < 2 * 255; ++jj)
        exptab[jj] = exptab[jj % 255];
    exptab[2 * 255] = 1;
    for (unsigned jj = 2 * 255 + 1; jj < 4 * 255; ++jj)
        exptab[jj] = 0;
}


//------------------------------------------------------------------------------
// Multiply and Divide Tables

// Initialize MUL and DIV tables using LOG and EXP tables
static void gf256_muldiv_init()
{
    // Allocate table memory 65KB x 2
    uint8_t* m = GF256Ctx.GF256_MUL_TABLE;
    uint8_t* d = GF256Ctx.GF256_DIV_TABLE;

    // Unroll y = 0 subtable
    for (int x = 0; x < 256; ++x)
        m[x] = d[x] = 0;

    // For each other y value:
    for (int y = 1; y < 256; ++y)
    {
        // Calculate log(y) for mult and 255 - log(y) for div
        const uint8_t log_y = static_cast<uint8_t>(GF256Ctx.GF256_LOG_TABLE[y]);
        const uint8_t log_yn = 255 - log_y;

        // Next subtable
        m += 256, d += 256;

        // Unroll x = 0
        m[0] = 0, d[0] = 0;

        // Calculate x * y, x / y
        for (int x = 1; x < 256; ++x)
        {
            uint16_t log_x = GF256Ctx.GF256_LOG_TABLE[x];

            m[x] = GF256Ctx.GF256_EXP_TABLE[log_x + log_y];
            d[x] = GF256Ctx.GF256_EXP_TABLE[log_x + log_yn];
        }
    }
}


//------------------------------------------------------------------------------
// Inverse Table

// Initialize INV table using DIV table
static void gf256_inv_init()
{
    for (int x = 0; x < 256; ++x)
        GF256Ctx.GF256_INV_TABLE[x] = gf256_div(1, static_cast<uint8_t>(x));
}


//------------------------------------------------------------------------------
// Square Table

// Initialize SQR table using MUL table
static void gf256_sqr_init()
{
    for (int x = 0; x < 256; ++x)
        GF256Ctx.GF256_SQR_TABLE[x] = gf256_mul(static_cast<uint8_t>(x), static_cast<uint8_t>(x));
}


//------------------------------------------------------------------------------
// Multiply and Add Memory Tables

/*
    Fast algorithm to compute m[1..8] = a[1..8] * b in GF(256)
    using SSE3 SIMD instruction set:

    Consider z = x * y in GF(256).
    This operation can be performed bit-by-bit.  Usefully, the partial product
    of each bit is combined linearly with the rest.  This means that the 8-bit
    number x can be split into its high and low 4 bits, and partial products
    can be formed from each half.  Then the halves can be linearly combined:

        z = x[0..3] * y + x[4..7] * y

    The multiplication of each half can be done efficiently via table lookups,
    and the addition in GF(256) is XOR.  There must be two tables that map 16
    input elements for the low or high 4 bits of x to the two partial products.
    Each value for y has a different set of two tables:

        z = TABLE_LO_y(x[0..3]) xor TABLE_HI_y(x[4..7])

    This means that we need 16 * 2 * 256 = 8192 bytes for precomputed tables.

    Computing z[] = x[] * y can be performed 16 bytes at a time by using the
    128-bit register operations supported by modern processors.

    This is efficiently realized in SSE3 using the _mm_shuffle_epi8() function
    provided by Visual Studio 2010 or newer in <tmmintrin.h>.  This function
    uses the low bits to do a table lookup on each byte.  Unfortunately the
    high bit of each mask byte has the special feature that it clears the
    output byte when it is set, so we need to make sure it's cleared by masking
    off the high bit of each byte before using it:

        clr_mask = _mm_set1_epi8(0x0f) = 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f

    For the low half of the partial product, clear the high bit of each byte
    and perform the table lookup:

        p_lo = _mm_and_si128(x, clr_mask)
        p_lo = _mm_shuffle_epi8(p_lo, TABLE_LO_y)

    For the high half of the partial product, shift the high 4 bits of each
    byte into the low 4 bits and clear the high bit of each byte, and then
    perform the table lookup:

        p_hi = _mm_srli_epi64(x, 4)
        p_hi = _mm_and_si128(p_hi, clr_mask)
        p_hi = _mm_shuffle_epi8(p_hi, TABLE_HI_y)

    Finally add the two partial products to form the product, recalling that
    addition is XOR in a Galois field:

        result = _mm_xor_si128(p_lo, p_hi)

    This crunches 16 bytes of x at a time, and the result can be stored in z.
*/

/*
    Intrinsic reference:

    SSE3, VS2010+, tmmintrin.h:

    GF256_M128 _mm_shuffle_epi8(GF256_M128 a, GF256_M128 mask);
        Emits the Supplemental Streaming SIMD Extensions 3 (SSSE3) instruction pshufb. This instruction shuffles 16-byte parameters from a 128-bit parameter.

        Pseudo-code for PSHUFB (with 128 bit operands):

            for i = 0 to 15 {
                 if (SRC[(i * 8)+7] = 1 ) then
                      DEST[(i*8)+7..(i*8)+0] <- 0;
                  else
                      index[3..0] <- SRC[(i*8)+3 .. (i*8)+0];
                      DEST[(i*8)+7..(i*8)+0] <- DEST[(index*8+7)..(index*8+0)];
                 endif
            }

    SSE2, VS2008+, emmintrin.h:

    GF256_M128 _mm_slli_epi64 (GF256_M128 a, int count);
        Shifts the 2 signed or unsigned 64-bit integers in a left by count bits while shifting in zeros.
    GF256_M128 _mm_srli_epi64 (GF256_M128 a, int count);
        Shifts the 2 signed or unsigned 64-bit integers in a right by count bits while shifting in zeros.
    GF256_M128 _mm_set1_epi8 (char b);
        Sets the 16 signed 8-bit integer values to b.
    GF256_M128 _mm_and_si128 (GF256_M128 a, GF256_M128 b);
        Computes the bitwise AND of the 128-bit value in a and the 128-bit value in b.
    GF256_M128 _mm_xor_si128 ( GF256_M128 a, GF256_M128 b);
        Computes the bitwise XOR of the 128-bit value in a and the 128-bit value in b.
*/

// Initialize the multiplication tables using gf256_mul()
static void gf256_mul_mem_init()
{
    // Reuse aligned self test buffers to load table data
    uint8_t* lo = m_SelfTestBuffers.A;
    uint8_t* hi = m_SelfTestBuffers.B;

    for (int y = 0; y < 256; ++y)
    {
        // TABLE_LO_Y maps 0..15 to 8-bit partial product based on y.
        for (unsigned char x = 0; x < 16; ++x)
        {
            lo[x] = gf256_mul(x, static_cast<uint8_t>( y ));
            hi[x] = gf256_mul(x << 4, static_cast<uint8_t>( y ));
        }

#if defined(GF256_TRY_NEON)
        if (CpuHasNeon)
        {
            GF256Ctx.MM128.TABLE_LO_Y[y] = vld1q_u8(lo);
            GF256Ctx.MM128.TABLE_HI_Y[y] = vld1q_u8(hi);
        }
#elif !defined(GF256_TARGET_MOBILE)
        const GF256_M128 table_lo = _mm_loadu_si128((GF256_M128*)lo);
        const GF256_M128 table_hi = _mm_loadu_si128((GF256_M128*)hi);
        _mm_storeu_si128(GF256Ctx.MM128.TABLE_LO_Y + y, table_lo);
        _mm_storeu_si128(GF256Ctx.MM128.TABLE_HI_Y + y, table_hi);
# ifdef GF256_TRY_AVX2
        if (CpuHasAVX2)
        {
            const GF256_M256 table_lo2 = _mm256_broadcastsi128_si256(table_lo);
            const GF256_M256 table_hi2 = _mm256_broadcastsi128_si256(table_hi);
            _mm256_storeu_si256(GF256Ctx.MM256.TABLE_LO_Y + y, table_lo2);
            _mm256_storeu_si256(GF256Ctx.MM256.TABLE_HI_Y + y, table_hi2);
        }
# endif // GF256_TRY_AVX2
#endif // GF256_TARGET_MOBILE
    }
}


//------------------------------------------------------------------------------
// GFNI affine matrices (Ablation S1)

#ifdef GF256_TRY_GFNI
// Build the 8x8 GF(2)-affine matrix (packed into a qword) such that, under the
// vgf2p8affineqb convention [ out_bit j = parity(x & ((M>>(8*j))&0xff)) ], applying
// M to byte x yields gf256_mul(x, c) in Wirehair's field (poly 0x14D).
//
// Intel SDM convention: out.bit[i] = parity(x AND M.byte[7-i]) XOR imm.bit[i].
// gf256_mul(x,c) is GF(2)-linear in x: column k (image of basis bit k) is
// p_k = gf256_mul(c, 1<<k). We need parity(x AND W_i)==out.bit[i] with W_i.bit[k]=p_k.bit[i],
// and M.byte[7-i]==W_i, i.e. output-bit j lands in matrix byte (7-j).
// Verified: M(c=0)=0, M(c=1)=0x0102040810204080 (the GFNI identity matrix).
static uint64_t gf256_gfni_matrix(uint8_t c)
{
    uint64_t m = 0;
    for (unsigned i = 0; i < 8; ++i) // i = input bit index k
    {
        const uint8_t p = gf256_mul(c, (uint8_t)(1u << i));
        for (unsigned j = 0; j < 8; ++j) // j = output bit index
            m |= (uint64_t)((p >> j) & 1) << (8 * (7 - j) + i);
    }
    return m;
}

static void gf256_gfni_init()
{
    for (int c = 0; c < 256; ++c)
        GF256Ctx.GFNI_MUL_MATRIX[c] = gf256_gfni_matrix((uint8_t)c);
}
#endif // GF256_TRY_GFNI


//------------------------------------------------------------------------------
// Initialization

#ifdef GF256_IS_BIG_ENDIAN
static unsigned char kEndianTestData[4] = { 1, 2, 3, 4 };
#else
static unsigned char kEndianTestData[4] = { 4, 3, 2, 1 };
#endif

union UnionType
{
    uint32_t IntValue;
    char CharArray[4];
};

static bool IsExpectedEndian()
{
    UnionType type;
    for (unsigned i = 0; i < 4; ++i)
        type.CharArray[i] = kEndianTestData[i];
    return 0x01020304 == type.IntValue;
}

extern "C" int gf256_init_(int version)
{
    if (version != GF256_VERSION)
        return -1; // User's header does not match library version.

    std::call_once(InitOnce, []() {
#if defined(WIREHAIR_TESTING)
        const wirehair::EnvironmentValue environment(
            "WIREHAIR_GF256_TEST_INIT_RESULT");
        const char* injected = environment.Get();
        if (injected && *injected) {
            const int result = std::atoi(injected);
            if (result != 0) {
                InitResult = result;
                return;
            }
        }
#endif
        if (!IsExpectedEndian()) {
            InitResult = -2;
            return;
        }

        gf256_architecture_init();
        gf256_poly_init(kDefaultPolynomialIndex);
        gf256_explog_init();
        gf256_muldiv_init();
        gf256_inv_init();
        gf256_sqr_init();
        gf256_mul_mem_init();
#ifdef GF256_TRY_GFNI
        gf256_gfni_init();
#endif

        if (!gf256_self_test()) {
            InitResult = -3;
            return;
        }

        InitResult = 0;
    });

    return InitResult;
}


//------------------------------------------------------------------------------
// Operations

// When the GFNI path handles the 64-byte bulk, the AVX2 mul/muladd blocks (which in
// gf256_mul_mem re-derive their pointers from the original args) must be skipped; the
// <64-byte remainder is finished by the SSSE3 + scalar tails using the advanced pointers.
#ifdef GF256_TRY_GFNI
# define GF256_GFNI_ACTIVE CpuHasGFNI
#else
# define GF256_GFNI_ACTIVE false
#endif

// When the AVX-512 XOR path consumes the bulk, the AVX2 XOR blocks must be skipped;
// the <64-byte remainder is finished by the SSSE3 + scalar tails (advanced pointers).
#ifdef GF256_TRY_AVX512
# define GF256_AVX512_ACTIVE CpuHasAVX512
#else
# define GF256_AVX512_ACTIVE false
#endif

#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_TARGET_AVX2)

static GF256_AVX2_TARGET int gf256_add_mem_avx2_target(
    void* GF256_RESTRICT vx,
    const void* GF256_RESTRICT vy,
    int bytes)
{
    uint8_t* const x = reinterpret_cast<uint8_t*>(vx);
    const uint8_t* const y = reinterpret_cast<const uint8_t*>(vy);
    const int vector_bytes = bytes & ~31;
    int offset = 0;
    for (; vector_bytes - offset >= 128; offset += 128)
    {
        for (int lane = 0; lane < 4; ++lane)
        {
            const int lane_offset = offset + lane * 32;
            const __m256i result = _mm256_xor_si256(
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(x + lane_offset)),
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(y + lane_offset)));
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(x + lane_offset), result);
        }
    }
    for (; offset < vector_bytes; offset += 32)
    {
        const __m256i result = _mm256_xor_si256(
            _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x + offset)),
            _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(y + offset)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x + offset), result);
    }
    return vector_bytes;
}

static GF256_AVX2_TARGET int gf256_add2_mem_avx2_target(
    void* GF256_RESTRICT vz,
    const void* GF256_RESTRICT vx,
    const void* GF256_RESTRICT vy,
    int bytes)
{
    uint8_t* const z = reinterpret_cast<uint8_t*>(vz);
    const uint8_t* const x = reinterpret_cast<const uint8_t*>(vx);
    const uint8_t* const y = reinterpret_cast<const uint8_t*>(vy);
    const int vector_bytes = bytes & ~31;
    int offset = 0;
    for (; vector_bytes - offset >= 128; offset += 128)
    {
        for (int lane = 0; lane < 4; ++lane)
        {
            const int lane_offset = offset + lane * 32;
            const __m256i result = _mm256_xor_si256(
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(z + lane_offset)),
                _mm256_xor_si256(
                    _mm256_loadu_si256(
                        reinterpret_cast<const __m256i*>(x + lane_offset)),
                    _mm256_loadu_si256(
                        reinterpret_cast<const __m256i*>(y + lane_offset))));
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(z + lane_offset), result);
        }
    }
    for (; offset < vector_bytes; offset += 32)
    {
        const __m256i result = _mm256_xor_si256(
            _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(z + offset)),
            _mm256_xor_si256(
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(x + offset)),
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(y + offset))));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(z + offset), result);
    }
    return vector_bytes;
}

static GF256_AVX2_TARGET int gf256_addset_mem_avx2_target(
    void* GF256_RESTRICT vz,
    const void* GF256_RESTRICT vx,
    const void* GF256_RESTRICT vy,
    int bytes)
{
    uint8_t* const z = reinterpret_cast<uint8_t*>(vz);
    const uint8_t* const x = reinterpret_cast<const uint8_t*>(vx);
    const uint8_t* const y = reinterpret_cast<const uint8_t*>(vy);
    const int vector_bytes = bytes & ~31;
    int offset = 0;
    for (; vector_bytes - offset >= 128; offset += 128)
    {
        for (int lane = 0; lane < 4; ++lane)
        {
            const int lane_offset = offset + lane * 32;
            const __m256i result = _mm256_xor_si256(
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(x + lane_offset)),
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(y + lane_offset)));
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(z + lane_offset), result);
        }
    }
    for (; offset < vector_bytes; offset += 32)
    {
        const __m256i result = _mm256_xor_si256(
            _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x + offset)),
            _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(y + offset)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(z + offset), result);
    }
    return vector_bytes;
}

#endif

extern "C" void gf256_add_mem(void * GF256_RESTRICT vx,
                              const void * GF256_RESTRICT vy, int bytes)
{
    if (bytes <= 0) {
        return;
    }

    WH_BUMP(0, bytes);
    GF256_M128 * GF256_RESTRICT x16 = reinterpret_cast<GF256_M128 *>(vx);
    const GF256_M128 * GF256_RESTRICT y16 = reinterpret_cast<const GF256_M128 *>(vy);

#if defined(GF256_TARGET_MOBILE)
# if defined(GF256_TRY_NEON)
    // Handle multiples of 64 bytes
    if (CpuHasNeon)
    {
        while (bytes >= 64)
        {
            GF256_M128 x0 = vld1q_u8((uint8_t*) x16);
            GF256_M128 x1 = vld1q_u8((uint8_t*)(x16 + 1) );
            GF256_M128 x2 = vld1q_u8((uint8_t*)(x16 + 2) );
            GF256_M128 x3 = vld1q_u8((uint8_t*)(x16 + 3) );
            GF256_M128 y0 = vld1q_u8((uint8_t*)y16);
            GF256_M128 y1 = vld1q_u8((uint8_t*)(y16 + 1));
            GF256_M128 y2 = vld1q_u8((uint8_t*)(y16 + 2));
            GF256_M128 y3 = vld1q_u8((uint8_t*)(y16 + 3));

            vst1q_u8((uint8_t*)x16,     veorq_u8(x0, y0));
            vst1q_u8((uint8_t*)(x16 + 1), veorq_u8(x1, y1));
            vst1q_u8((uint8_t*)(x16 + 2), veorq_u8(x2, y2));
            vst1q_u8((uint8_t*)(x16 + 3), veorq_u8(x3, y3));

            bytes -= 64, x16 += 4, y16 += 4;
        }

        // Handle multiples of 16 bytes
        while (bytes >= 16)
        {
            GF256_M128 x0 = vld1q_u8((uint8_t*)x16);
            GF256_M128 y0 = vld1q_u8((uint8_t*)y16);

            vst1q_u8((uint8_t*)x16, veorq_u8(x0, y0));

            bytes -= 16, ++x16, ++y16;
        }
    }
    else
# endif // GF256_TRY_NEON
    {
        uint8_t * GF256_RESTRICT x8 = reinterpret_cast<uint8_t *>(x16);
        const uint8_t * GF256_RESTRICT y8 = reinterpret_cast<const uint8_t *>(y16);

        const unsigned count = (unsigned)bytes / 8;
        for (unsigned ii = 0; ii < count; ++ii)
        {
            const unsigned offset = ii * 8u;
            gf256_storeu64(
                x8 + offset,
                gf256_loadu64(x8 + offset) ^ gf256_loadu64(y8 + offset));
        }

        x16 = reinterpret_cast<GF256_M128 *>(x8 + count * 8u);
        y16 = reinterpret_cast<const GF256_M128 *>(y8 + count * 8u);

        bytes -= (count * 8);
    }
#else // GF256_TARGET_MOBILE
# if defined(GF256_TRY_AVX512)
    // Ablation S2: 64-byte ZMM XOR, 4x unrolled (256 bytes/iter).
    if (CpuHasAVX512)
    {
        while (bytes >= 256)
        {
            __m512i a0 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 0)),  _mm512_loadu_si512((const void*)(y16 + 0)));
            __m512i a1 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 4)),  _mm512_loadu_si512((const void*)(y16 + 4)));
            __m512i a2 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 8)),  _mm512_loadu_si512((const void*)(y16 + 8)));
            __m512i a3 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 12)), _mm512_loadu_si512((const void*)(y16 + 12)));
            _mm512_storeu_si512((void*)(x16 + 0),  a0);
            _mm512_storeu_si512((void*)(x16 + 4),  a1);
            _mm512_storeu_si512((void*)(x16 + 8),  a2);
            _mm512_storeu_si512((void*)(x16 + 12), a3);
            bytes -= 256, x16 += 16, y16 += 16; // 256 = 16 * M128
        }
        while (bytes >= 64)
        {
            _mm512_storeu_si512((void*)x16, _mm512_xor_si512(_mm512_loadu_si512((const void*)x16), _mm512_loadu_si512((const void*)y16)));
            bytes -= 64, x16 += 4, y16 += 4;
        }
    }
# endif // GF256_TRY_AVX512
# if defined(GF256_TRY_TARGET_AVX2)
    if (bytes >= 32 && CpuHasAVX2 && !GF256_AVX512_ACTIVE)
    {
        const int vector_bytes =
            gf256_add_mem_avx2_target(x16, y16, bytes);
        bytes -= vector_bytes;
        x16 += vector_bytes / 16;
        y16 += vector_bytes / 16;
    }
# endif // GF256_TRY_TARGET_AVX2
# if defined(GF256_TRY_AVX2)
    if (CpuHasAVX2 && !GF256_AVX512_ACTIVE)
    {
        GF256_M256 * GF256_RESTRICT x32 = reinterpret_cast<GF256_M256 *>(x16);
        const GF256_M256 * GF256_RESTRICT y32 = reinterpret_cast<const GF256_M256 *>(y16);

        while (bytes >= 128)
        {
            GF256_M256 x0 = _mm256_loadu_si256(x32);
            GF256_M256 y0 = _mm256_loadu_si256(y32);
            x0 = _mm256_xor_si256(x0, y0);
            GF256_M256 x1 = _mm256_loadu_si256(x32 + 1);
            GF256_M256 y1 = _mm256_loadu_si256(y32 + 1);
            x1 = _mm256_xor_si256(x1, y1);
            GF256_M256 x2 = _mm256_loadu_si256(x32 + 2);
            GF256_M256 y2 = _mm256_loadu_si256(y32 + 2);
            x2 = _mm256_xor_si256(x2, y2);
            GF256_M256 x3 = _mm256_loadu_si256(x32 + 3);
            GF256_M256 y3 = _mm256_loadu_si256(y32 + 3);
            x3 = _mm256_xor_si256(x3, y3);

            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
            _mm256_storeu_si256(x32 + 2, x2);
            _mm256_storeu_si256(x32 + 3, x3);

            bytes -= 128, x32 += 4, y32 += 4;
        }

        // Handle multiples of 32 bytes
        while (bytes >= 32)
        {
            // x[i] = x[i] xor y[i]
            _mm256_storeu_si256(x32,
                _mm256_xor_si256(
                    _mm256_loadu_si256(x32),
                    _mm256_loadu_si256(y32)));

            bytes -= 32, ++x32, ++y32;
        }

        x16 = reinterpret_cast<GF256_M128 *>(x32);
        y16 = reinterpret_cast<const GF256_M128 *>(y32);
    }
    else
# endif // GF256_TRY_AVX2
    {
        while (bytes >= 64)
        {
            GF256_M128 x0 = _mm_loadu_si128(x16);
            GF256_M128 y0 = _mm_loadu_si128(y16);
            x0 = _mm_xor_si128(x0, y0);
            GF256_M128 x1 = _mm_loadu_si128(x16 + 1);
            GF256_M128 y1 = _mm_loadu_si128(y16 + 1);
            x1 = _mm_xor_si128(x1, y1);
            GF256_M128 x2 = _mm_loadu_si128(x16 + 2);
            GF256_M128 y2 = _mm_loadu_si128(y16 + 2);
            x2 = _mm_xor_si128(x2, y2);
            GF256_M128 x3 = _mm_loadu_si128(x16 + 3);
            GF256_M128 y3 = _mm_loadu_si128(y16 + 3);
            x3 = _mm_xor_si128(x3, y3);

            _mm_storeu_si128(x16, x0);
            _mm_storeu_si128(x16 + 1, x1);
            _mm_storeu_si128(x16 + 2, x2);
            _mm_storeu_si128(x16 + 3, x3);

            bytes -= 64, x16 += 4, y16 += 4;
        }
    }
#endif // GF256_TARGET_MOBILE

#if !defined(GF256_TARGET_MOBILE)
    // Handle multiples of 16 bytes
    while (bytes >= 16)
    {
        // x[i] = x[i] xor y[i]
        _mm_storeu_si128(x16,
            _mm_xor_si128(
                _mm_loadu_si128(x16),
                _mm_loadu_si128(y16)));

        bytes -= 16, ++x16, ++y16;
    }
#endif

    uint8_t * GF256_RESTRICT x1 = reinterpret_cast<uint8_t *>(x16);
    const uint8_t * GF256_RESTRICT y1 = reinterpret_cast<const uint8_t *>(y16);

    // Handle a block of 8 bytes
    const int eight = bytes & 8;
    if (eight)
    {
        gf256_storeu64(x1, gf256_loadu64(x1) ^ gf256_loadu64(y1));
    }

    // Handle a block of 4 bytes
    const int four = bytes & 4;
    if (four)
    {
        gf256_storeu32(
            x1 + eight,
            gf256_loadu32(x1 + eight) ^ gf256_loadu32(y1 + eight));
    }

    // Handle final bytes
    const int offset = eight + four;
    switch (bytes & 3)
    {
    case 3: x1[offset + 2] ^= y1[offset + 2]; GF256_FALLTHROUGH;
    case 2: x1[offset + 1] ^= y1[offset + 1]; GF256_FALLTHROUGH;
    case 1: x1[offset] ^= y1[offset];
    default:
        break;
    }
}

extern "C" void gf256_add2_mem(void * GF256_RESTRICT vz, const void * GF256_RESTRICT vx,
                               const void * GF256_RESTRICT vy, int bytes)
{
    if (bytes <= 0) {
        return;
    }

    WH_BUMP(1, bytes);
    GF256_M128 * GF256_RESTRICT z16 = reinterpret_cast<GF256_M128*>(vz);
    const GF256_M128 * GF256_RESTRICT x16 = reinterpret_cast<const GF256_M128*>(vx);
    const GF256_M128 * GF256_RESTRICT y16 = reinterpret_cast<const GF256_M128*>(vy);

#if defined(GF256_TARGET_MOBILE)
# if defined(GF256_TRY_NEON)
    // Handle multiples of 64 bytes
    if (CpuHasNeon)
    {
        // Handle multiples of 16 bytes
        while (bytes >= 16)
        {
            // z[i] = z[i] xor x[i] xor y[i]
            vst1q_u8((uint8_t*)z16,
                veorq_u8(
                    vld1q_u8((uint8_t*)z16),
                    veorq_u8(
                        vld1q_u8((uint8_t*)x16),
                        vld1q_u8((uint8_t*)y16))));

            bytes -= 16, ++x16, ++y16, ++z16;
        }
    }
    else
# endif // GF256_TRY_NEON
    {
        uint8_t * GF256_RESTRICT z8 = reinterpret_cast<uint8_t *>(z16);
        const uint8_t * GF256_RESTRICT x8 = reinterpret_cast<const uint8_t *>(x16);
        const uint8_t * GF256_RESTRICT y8 = reinterpret_cast<const uint8_t *>(y16);

        const unsigned count = (unsigned)bytes / 8;
        for (unsigned ii = 0; ii < count; ++ii)
        {
            const unsigned offset = ii * 8u;
            gf256_storeu64(
                z8 + offset,
                gf256_loadu64(z8 + offset) ^
                gf256_loadu64(x8 + offset) ^
                gf256_loadu64(y8 + offset));
        }

        z16 = reinterpret_cast<GF256_M128 *>(z8 + count * 8u);
        x16 = reinterpret_cast<const GF256_M128 *>(x8 + count * 8u);
        y16 = reinterpret_cast<const GF256_M128 *>(y8 + count * 8u);

        bytes -= (count * 8);
    }
#else // GF256_TARGET_MOBILE
# if defined(GF256_TRY_AVX512)
    // Ablation S2: z[] ^= x[] ^ y[], 64-byte ZMM, 4x unrolled (256 bytes/iter).
    if (CpuHasAVX512)
    {
        while (bytes >= 256)
        {
            __m512i z0 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(z16 + 0)),  _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 0)),  _mm512_loadu_si512((const void*)(y16 + 0))));
            __m512i z1 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(z16 + 4)),  _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 4)),  _mm512_loadu_si512((const void*)(y16 + 4))));
            __m512i z2 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(z16 + 8)),  _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 8)),  _mm512_loadu_si512((const void*)(y16 + 8))));
            __m512i z3 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(z16 + 12)), _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 12)), _mm512_loadu_si512((const void*)(y16 + 12))));
            _mm512_storeu_si512((void*)(z16 + 0),  z0);
            _mm512_storeu_si512((void*)(z16 + 4),  z1);
            _mm512_storeu_si512((void*)(z16 + 8),  z2);
            _mm512_storeu_si512((void*)(z16 + 12), z3);
            bytes -= 256, z16 += 16, x16 += 16, y16 += 16;
        }
        while (bytes >= 64)
        {
            __m512i z0 = _mm512_xor_si512(_mm512_loadu_si512((const void*)z16), _mm512_xor_si512(_mm512_loadu_si512((const void*)x16), _mm512_loadu_si512((const void*)y16)));
            _mm512_storeu_si512((void*)z16, z0);
            bytes -= 64, z16 += 4, x16 += 4, y16 += 4;
        }
    }
# endif // GF256_TRY_AVX512
# if defined(GF256_TRY_TARGET_AVX2)
    if (bytes >= 32 && CpuHasAVX2 && !GF256_AVX512_ACTIVE)
    {
        const int vector_bytes =
            gf256_add2_mem_avx2_target(z16, x16, y16, bytes);
        bytes -= vector_bytes;
        z16 += vector_bytes / 16;
        x16 += vector_bytes / 16;
        y16 += vector_bytes / 16;
    }
# endif // GF256_TRY_TARGET_AVX2
# if defined(GF256_TRY_AVX2)
    if (CpuHasAVX2 && !GF256_AVX512_ACTIVE)
    {
        GF256_M256 * GF256_RESTRICT z32 = reinterpret_cast<GF256_M256 *>(z16);
        const GF256_M256 * GF256_RESTRICT x32 = reinterpret_cast<const GF256_M256 *>(x16);
        const GF256_M256 * GF256_RESTRICT y32 = reinterpret_cast<const GF256_M256 *>(y16);

        const unsigned count = bytes / 32;
        for (unsigned i = 0; i < count; ++i)
        {
            _mm256_storeu_si256(z32 + i,
                _mm256_xor_si256(
                    _mm256_loadu_si256(z32 + i),
                    _mm256_xor_si256(
                        _mm256_loadu_si256(x32 + i),
                        _mm256_loadu_si256(y32 + i))));
        }

        bytes -= count * 32;
        z16 = reinterpret_cast<GF256_M128 *>(z32 + count);
        x16 = reinterpret_cast<const GF256_M128 *>(x32 + count);
        y16 = reinterpret_cast<const GF256_M128 *>(y32 + count);
    }
# endif // GF256_TRY_AVX2

    // Handle multiples of 16 bytes
    while (bytes >= 16)
    {
        // z[i] = z[i] xor x[i] xor y[i]
        _mm_storeu_si128(z16,
            _mm_xor_si128(
                _mm_loadu_si128(z16),
                _mm_xor_si128(
                    _mm_loadu_si128(x16),
                    _mm_loadu_si128(y16))));

        bytes -= 16, ++x16, ++y16, ++z16;
    }
#endif // GF256_TARGET_MOBILE

    uint8_t * GF256_RESTRICT z1 = reinterpret_cast<uint8_t *>(z16);
    const uint8_t * GF256_RESTRICT x1 = reinterpret_cast<const uint8_t *>(x16);
    const uint8_t * GF256_RESTRICT y1 = reinterpret_cast<const uint8_t *>(y16);

    // Handle a block of 8 bytes
    const int eight = bytes & 8;
    if (eight)
    {
        gf256_storeu64(
            z1,
            gf256_loadu64(z1) ^ gf256_loadu64(x1) ^ gf256_loadu64(y1));
    }

    // Handle a block of 4 bytes
    const int four = bytes & 4;
    if (four)
    {
        gf256_storeu32(
            z1 + eight,
            gf256_loadu32(z1 + eight) ^
            gf256_loadu32(x1 + eight) ^
            gf256_loadu32(y1 + eight));
    }

    // Handle final bytes
    const int offset = eight + four;
    switch (bytes & 3)
    {
    case 3: z1[offset + 2] ^= x1[offset + 2] ^ y1[offset + 2]; GF256_FALLTHROUGH;
    case 2: z1[offset + 1] ^= x1[offset + 1] ^ y1[offset + 1]; GF256_FALLTHROUGH;
    case 1: z1[offset] ^= x1[offset] ^ y1[offset];
    default:
        break;
    }
}

template<unsigned SrcCount, bool SetDestination>
static GF256_FORCE_INLINE void gf256_xor_multi_fixed(
    uint8_t * GF256_RESTRICT z,
    const uint8_t * const * GF256_RESTRICT srcs,
    int bytes)
{
    unsigned offset = 0;

#if defined(GF256_TARGET_MOBILE)
# if defined(GF256_TRY_NEON)
    if (CpuHasNeon)
    {
        while (bytes >= 16)
        {
            uint8x16_t acc = SetDestination ?
                vld1q_u8(srcs[0] + offset) : vld1q_u8(z + offset);
            const unsigned first_source = SetDestination ? 1u : 0u;
            for (unsigned j = first_source; j < SrcCount; ++j) {
                acc = veorq_u8(acc, vld1q_u8(srcs[j] + offset));
            }
            vst1q_u8(z + offset, acc);
            offset += 16u;
            bytes -= 16;
        }
    }
# endif
#else
# if defined(GF256_TRY_AVX512)
    if (CpuHasAVX512)
    {
        while (bytes >= 64)
        {
            __m512i acc = _mm512_loadu_si512((const void*)(
                SetDestination ? srcs[0] + offset : z + offset));
            const unsigned first_source = SetDestination ? 1u : 0u;
            for (unsigned j = first_source; j < SrcCount; ++j) {
                acc = _mm512_xor_si512(
                    acc,
                    _mm512_loadu_si512((const void*)(srcs[j] + offset)));
            }
            _mm512_storeu_si512((void*)(z + offset), acc);
            offset += 64u;
            bytes -= 64;
        }
    }
# endif
# if defined(GF256_TRY_AVX2)
    if (CpuHasAVX2)
    {
        // Keep two independent dependency chains in flight.  The fixed
        // source count lets the compiler interleave their loads and XORs,
        // while one loop branch now covers a full cache line.
        while (bytes >= 64)
        {
            __m256i acc0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    SetDestination ? srcs[0] + offset : z + offset));
            __m256i acc1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    SetDestination ?
                        srcs[0] + offset + 32u : z + offset + 32u));
            const unsigned first_source = SetDestination ? 1u : 0u;
            for (unsigned j = first_source; j < SrcCount; ++j)
            {
                acc0 = _mm256_xor_si256(
                    acc0,
                    _mm256_loadu_si256(
                        reinterpret_cast<const __m256i*>(
                            srcs[j] + offset)));
                acc1 = _mm256_xor_si256(
                    acc1,
                    _mm256_loadu_si256(
                        reinterpret_cast<const __m256i*>(
                            srcs[j] + offset + 32u)));
            }
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(z + offset), acc0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(z + offset + 32u), acc1);
            offset += 64u;
            bytes -= 64;
        }
        while (bytes >= 32)
        {
            __m256i acc = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    SetDestination ? srcs[0] + offset : z + offset));
            const unsigned first_source = SetDestination ? 1u : 0u;
            for (unsigned j = first_source; j < SrcCount; ++j) {
                acc = _mm256_xor_si256(
                    acc,
                    _mm256_loadu_si256(
                        reinterpret_cast<const __m256i*>(
                            srcs[j] + offset)));
            }
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(z + offset), acc);
            offset += 32u;
            bytes -= 32;
        }
    }
# endif
    while (bytes >= 16)
    {
        __m128i acc = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(
                SetDestination ? srcs[0] + offset : z + offset));
        const unsigned first_source = SetDestination ? 1u : 0u;
        for (unsigned j = first_source; j < SrcCount; ++j) {
            acc = _mm_xor_si128(
                acc,
                _mm_loadu_si128(
                    reinterpret_cast<const __m128i*>(srcs[j] + offset)));
        }
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(z + offset), acc);
        offset += 16u;
        bytes -= 16;
    }
#endif

    while (bytes >= 8)
    {
        uint64_t acc = gf256_loadu64(
            SetDestination ? srcs[0] + offset : z + offset);
        const unsigned first_source = SetDestination ? 1u : 0u;
        for (unsigned j = first_source; j < SrcCount; ++j) {
            acc ^= gf256_loadu64(srcs[j] + offset);
        }
        gf256_storeu64(z + offset, acc);
        offset += 8;
        bytes -= 8;
    }

    while (bytes-- > 0)
    {
        uint8_t acc = SetDestination ? srcs[0][offset] : z[offset];
        const unsigned first_source = SetDestination ? 1u : 0u;
        for (unsigned j = first_source; j < SrcCount; ++j) {
            acc ^= srcs[j][offset];
        }
        z[offset++] = acc;
    }
}

template<bool SetDestination>
static GF256_FORCE_INLINE void gf256_xor_multi_generic(
    uint8_t * GF256_RESTRICT z,
    const uint8_t * const * GF256_RESTRICT srcs,
    int src_count,
    int bytes)
{
    unsigned offset = 0;

    while (bytes >= 8)
    {
        uint64_t acc = gf256_loadu64(
            SetDestination ? srcs[0] + offset : z + offset);
        const int first_source = SetDestination ? 1 : 0;
        for (int j = first_source; j < src_count; ++j) {
            acc ^= gf256_loadu64(srcs[j] + offset);
        }
        gf256_storeu64(z + offset, acc);
        offset += 8;
        bytes -= 8;
    }

    while (bytes-- > 0)
    {
        uint8_t acc = SetDestination ? srcs[0][offset] : z[offset];
        const int first_source = SetDestination ? 1 : 0;
        for (int j = first_source; j < src_count; ++j) {
            acc ^= srcs[j][offset];
        }
        z[offset++] = acc;
    }
}

extern "C" void gf256_add_multi_mem(
    void * GF256_RESTRICT vz,
    const void * const * GF256_RESTRICT vsrcs,
    int src_count,
    int bytes)
{
    if (bytes <= 0 || src_count <= 0) {
        return;
    }
    if (src_count == 1)
    {
        gf256_add_mem(vz, vsrcs[0], bytes);
        return;
    }
    if (src_count == 2)
    {
        gf256_add2_mem(vz, vsrcs[0], vsrcs[1], bytes);
        return;
    }

    WH_BUMP(0, bytes);

    uint8_t * GF256_RESTRICT z = reinterpret_cast<uint8_t *>(vz);
    const uint8_t * srcs[16];
    if (src_count <= (int)(sizeof(srcs) / sizeof(srcs[0])))
    {
        for (int j = 0; j < src_count; ++j) {
            srcs[j] = reinterpret_cast<const uint8_t *>(vsrcs[j]);
        }

        switch (src_count)
        {
        case 3:  gf256_xor_multi_fixed<3, false>(z, srcs, bytes);  return;
        case 4:  gf256_xor_multi_fixed<4, false>(z, srcs, bytes);  return;
        case 5:  gf256_xor_multi_fixed<5, false>(z, srcs, bytes);  return;
        case 6:  gf256_xor_multi_fixed<6, false>(z, srcs, bytes);  return;
        case 7:  gf256_xor_multi_fixed<7, false>(z, srcs, bytes);  return;
        case 8:  gf256_xor_multi_fixed<8, false>(z, srcs, bytes);  return;
        case 12: gf256_xor_multi_fixed<12, false>(z, srcs, bytes); return;
        case 16: gf256_xor_multi_fixed<16, false>(z, srcs, bytes); return;
        default:
            gf256_xor_multi_generic<false>(z, srcs, src_count, bytes);
            return;
        }
    }

    const uint8_t * src_window[16];
    int remaining = src_count;
    const void * const * next = vsrcs;
    while (remaining > 0)
    {
        const int count = remaining > 16 ? 16 : remaining;
        for (int j = 0; j < count; ++j) {
            src_window[j] = reinterpret_cast<const uint8_t *>(next[j]);
        }
        gf256_xor_multi_generic<false>(z, src_window, count, bytes);
        next += count;
        remaining -= count;
    }
}

extern "C" void gf256_addset_multi_mem(
    void * GF256_RESTRICT vz,
    const void * const * GF256_RESTRICT vsrcs,
    int src_count,
    int bytes)
{
    if (bytes <= 0 || src_count <= 0) {
        return;
    }
    if (src_count == 1)
    {
        std::memcpy(vz, vsrcs[0], (size_t)bytes);
        return;
    }
    if (src_count == 2)
    {
        gf256_addset_mem(vz, vsrcs[0], vsrcs[1], bytes);
        return;
    }

    WH_BUMP(2, bytes);

    uint8_t * GF256_RESTRICT z = reinterpret_cast<uint8_t *>(vz);
    const uint8_t * srcs[16];
    if (src_count <= (int)(sizeof(srcs) / sizeof(srcs[0])))
    {
        for (int j = 0; j < src_count; ++j) {
            srcs[j] = reinterpret_cast<const uint8_t *>(vsrcs[j]);
        }

        switch (src_count)
        {
        case 3:  gf256_xor_multi_fixed<3, true>(z, srcs, bytes);  return;
        case 4:  gf256_xor_multi_fixed<4, true>(z, srcs, bytes);  return;
        case 5:  gf256_xor_multi_fixed<5, true>(z, srcs, bytes);  return;
        case 6:  gf256_xor_multi_fixed<6, true>(z, srcs, bytes);  return;
        case 7:  gf256_xor_multi_fixed<7, true>(z, srcs, bytes);  return;
        case 8:  gf256_xor_multi_fixed<8, true>(z, srcs, bytes);  return;
        case 9:  gf256_xor_multi_fixed<9, true>(z, srcs, bytes);  return;
        case 10: gf256_xor_multi_fixed<10, true>(z, srcs, bytes); return;
        case 11: gf256_xor_multi_fixed<11, true>(z, srcs, bytes); return;
        case 12: gf256_xor_multi_fixed<12, true>(z, srcs, bytes); return;
        case 13: gf256_xor_multi_fixed<13, true>(z, srcs, bytes); return;
        case 14: gf256_xor_multi_fixed<14, true>(z, srcs, bytes); return;
        case 15: gf256_xor_multi_fixed<15, true>(z, srcs, bytes); return;
        case 16: gf256_xor_multi_fixed<16, true>(z, srcs, bytes); return;
        default:
            gf256_xor_multi_generic<true>(z, srcs, src_count, bytes);
            return;
        }
    }

    const uint8_t * src_window[16];
    int remaining = src_count;
    const void * const * next = vsrcs;
    bool initialized = false;
    while (remaining > 0)
    {
        const int count = remaining > 16 ? 16 : remaining;
        for (int j = 0; j < count; ++j) {
            src_window[j] = reinterpret_cast<const uint8_t *>(next[j]);
        }
        if (initialized) {
            gf256_xor_multi_generic<false>(z, src_window, count, bytes);
        }
        else {
            gf256_xor_multi_generic<true>(z, src_window, count, bytes);
            initialized = true;
        }
        next += count;
        remaining -= count;
    }
}

extern "C" void gf256_addset_mem(void * GF256_RESTRICT vz, const void * GF256_RESTRICT vx,
                                 const void * GF256_RESTRICT vy, int bytes)
{
    if (bytes <= 0) {
        return;
    }

    WH_BUMP(2, bytes);
    GF256_M128 * GF256_RESTRICT z16 = reinterpret_cast<GF256_M128*>(vz);
    const GF256_M128 * GF256_RESTRICT x16 = reinterpret_cast<const GF256_M128*>(vx);
    const GF256_M128 * GF256_RESTRICT y16 = reinterpret_cast<const GF256_M128*>(vy);

#if defined(GF256_TARGET_MOBILE)
# if defined(GF256_TRY_NEON)
    // Handle multiples of 64 bytes
    if (CpuHasNeon)
    {
        while (bytes >= 64)
        {
            GF256_M128 x0 = vld1q_u8((uint8_t*)x16);
            GF256_M128 x1 = vld1q_u8((uint8_t*)(x16 + 1));
            GF256_M128 x2 = vld1q_u8((uint8_t*)(x16 + 2));
            GF256_M128 x3 = vld1q_u8((uint8_t*)(x16 + 3));
            GF256_M128 y0 = vld1q_u8((uint8_t*)(y16));
            GF256_M128 y1 = vld1q_u8((uint8_t*)(y16 + 1));
            GF256_M128 y2 = vld1q_u8((uint8_t*)(y16 + 2));
            GF256_M128 y3 = vld1q_u8((uint8_t*)(y16 + 3));

            vst1q_u8((uint8_t*)z16,     veorq_u8(x0, y0));
            vst1q_u8((uint8_t*)(z16 + 1), veorq_u8(x1, y1));
            vst1q_u8((uint8_t*)(z16 + 2), veorq_u8(x2, y2));
            vst1q_u8((uint8_t*)(z16 + 3), veorq_u8(x3, y3));

            bytes -= 64, x16 += 4, y16 += 4, z16 += 4;
        }

        // Handle multiples of 16 bytes
        while (bytes >= 16)
        {
            // z[i] = x[i] xor y[i]
            vst1q_u8((uint8_t*)z16,
                     veorq_u8(
                         vld1q_u8((uint8_t*)x16),
                         vld1q_u8((uint8_t*)y16)));

            bytes -= 16, ++x16, ++y16, ++z16;
        }
    }
    else
# endif // GF256_TRY_NEON
    {
        uint8_t * GF256_RESTRICT z8 = reinterpret_cast<uint8_t *>(z16);
        const uint8_t * GF256_RESTRICT x8 = reinterpret_cast<const uint8_t *>(x16);
        const uint8_t * GF256_RESTRICT y8 = reinterpret_cast<const uint8_t *>(y16);

        const unsigned count = (unsigned)bytes / 8;
        for (unsigned ii = 0; ii < count; ++ii)
        {
            const unsigned offset = ii * 8u;
            gf256_storeu64(
                z8 + offset,
                gf256_loadu64(x8 + offset) ^ gf256_loadu64(y8 + offset));
        }

        x16 = reinterpret_cast<const GF256_M128 *>(x8 + count * 8u);
        y16 = reinterpret_cast<const GF256_M128 *>(y8 + count * 8u);
        z16 = reinterpret_cast<GF256_M128 *>(z8 + count * 8u);

        bytes -= (count * 8);
    }
#else // GF256_TARGET_MOBILE
# if defined(GF256_TRY_AVX512)
    // Ablation S2: z[] = x[] ^ y[], 64-byte ZMM, 4x unrolled (256 bytes/iter).
    if (CpuHasAVX512)
    {
        while (bytes >= 256)
        {
            __m512i z0 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 0)),  _mm512_loadu_si512((const void*)(y16 + 0)));
            __m512i z1 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 4)),  _mm512_loadu_si512((const void*)(y16 + 4)));
            __m512i z2 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 8)),  _mm512_loadu_si512((const void*)(y16 + 8)));
            __m512i z3 = _mm512_xor_si512(_mm512_loadu_si512((const void*)(x16 + 12)), _mm512_loadu_si512((const void*)(y16 + 12)));
            _mm512_storeu_si512((void*)(z16 + 0),  z0);
            _mm512_storeu_si512((void*)(z16 + 4),  z1);
            _mm512_storeu_si512((void*)(z16 + 8),  z2);
            _mm512_storeu_si512((void*)(z16 + 12), z3);
            bytes -= 256, z16 += 16, x16 += 16, y16 += 16;
        }
        while (bytes >= 64)
        {
            _mm512_storeu_si512((void*)z16, _mm512_xor_si512(_mm512_loadu_si512((const void*)x16), _mm512_loadu_si512((const void*)y16)));
            bytes -= 64, z16 += 4, x16 += 4, y16 += 4;
        }
    }
# endif // GF256_TRY_AVX512
# if defined(GF256_TRY_TARGET_AVX2)
    if (bytes >= 32 && CpuHasAVX2 && !GF256_AVX512_ACTIVE)
    {
        const int vector_bytes =
            gf256_addset_mem_avx2_target(z16, x16, y16, bytes);
        bytes -= vector_bytes;
        z16 += vector_bytes / 16;
        x16 += vector_bytes / 16;
        y16 += vector_bytes / 16;
    }
# endif // GF256_TRY_TARGET_AVX2
# if defined(GF256_TRY_AVX2)
    if (CpuHasAVX2 && !GF256_AVX512_ACTIVE)
    {
        GF256_M256 * GF256_RESTRICT z32 = reinterpret_cast<GF256_M256 *>(z16);
        const GF256_M256 * GF256_RESTRICT x32 = reinterpret_cast<const GF256_M256 *>(x16);
        const GF256_M256 * GF256_RESTRICT y32 = reinterpret_cast<const GF256_M256 *>(y16);

        const unsigned count = bytes / 32;
        for (unsigned i = 0; i < count; ++i)
        {
            _mm256_storeu_si256(z32 + i,
                _mm256_xor_si256(
                    _mm256_loadu_si256(x32 + i),
                    _mm256_loadu_si256(y32 + i)));
        }

        bytes -= count * 32;
        z16 = reinterpret_cast<GF256_M128 *>(z32 + count);
        x16 = reinterpret_cast<const GF256_M128 *>(x32 + count);
        y16 = reinterpret_cast<const GF256_M128 *>(y32 + count);
    }
    else
# endif // GF256_TRY_AVX2
    {
        // Handle multiples of 64 bytes
        while (bytes >= 64)
        {
            GF256_M128 x0 = _mm_loadu_si128(x16);
            GF256_M128 x1 = _mm_loadu_si128(x16 + 1);
            GF256_M128 x2 = _mm_loadu_si128(x16 + 2);
            GF256_M128 x3 = _mm_loadu_si128(x16 + 3);
            GF256_M128 y0 = _mm_loadu_si128(y16);
            GF256_M128 y1 = _mm_loadu_si128(y16 + 1);
            GF256_M128 y2 = _mm_loadu_si128(y16 + 2);
            GF256_M128 y3 = _mm_loadu_si128(y16 + 3);

            _mm_storeu_si128(z16,     _mm_xor_si128(x0, y0));
            _mm_storeu_si128(z16 + 1, _mm_xor_si128(x1, y1));
            _mm_storeu_si128(z16 + 2, _mm_xor_si128(x2, y2));
            _mm_storeu_si128(z16 + 3, _mm_xor_si128(x3, y3));

            bytes -= 64, x16 += 4, y16 += 4, z16 += 4;
        }
    }

    // Handle multiples of 16 bytes
    while (bytes >= 16)
    {
        // z[i] = x[i] xor y[i]
        _mm_storeu_si128(z16,
            _mm_xor_si128(
                _mm_loadu_si128(x16),
                _mm_loadu_si128(y16)));

        bytes -= 16, ++x16, ++y16, ++z16;
    }
#endif // GF256_TARGET_MOBILE

    uint8_t * GF256_RESTRICT z1 = reinterpret_cast<uint8_t *>(z16);
    const uint8_t * GF256_RESTRICT x1 = reinterpret_cast<const uint8_t *>(x16);
    const uint8_t * GF256_RESTRICT y1 = reinterpret_cast<const uint8_t *>(y16);

    // Handle a block of 8 bytes
    const int eight = bytes & 8;
    if (eight)
    {
        gf256_storeu64(z1, gf256_loadu64(x1) ^ gf256_loadu64(y1));
    }

    // Handle a block of 4 bytes
    const int four = bytes & 4;
    if (four)
    {
        gf256_storeu32(
            z1 + eight,
            gf256_loadu32(x1 + eight) ^ gf256_loadu32(y1 + eight));
    }

    // Handle final bytes
    const int offset = eight + four;
    switch (bytes & 3)
    {
    case 3: z1[offset + 2] = x1[offset + 2] ^ y1[offset + 2]; GF256_FALLTHROUGH;
    case 2: z1[offset + 1] = x1[offset + 1] ^ y1[offset + 1]; GF256_FALLTHROUGH;
    case 1: z1[offset] = x1[offset] ^ y1[offset];
    default:
        break;
    }
}

#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_TARGET_AVX2)

// Runtime AVX2 helpers deliberately broadcast the portable 128-bit tables.
// This avoids changing gf256_ctx layout in baseline translation units while
// retaining the 32-byte nibble-shuffle kernels on capable x86 CPUs.
static GF256_AVX2_TARGET int gf256_mul_mem_avx2_target(
    void* vz, const void* vx, uint8_t y, int bytes)
{
    const __m256i table_low = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_LO_Y + y));
    const __m256i table_high = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_HI_Y + y));
    const __m256i clear_low = _mm256_set1_epi8(0x0f);
    const int vector_bytes = bytes & ~31;
    const uint8_t* const source = reinterpret_cast<const uint8_t*>(vx);
    uint8_t* const destination = reinterpret_cast<uint8_t*>(vz);
    for (int offset = 0; offset < vector_bytes; offset += 32)
    {
        __m256i input = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
        __m256i low = _mm256_and_si256(input, clear_low);
        input = _mm256_srli_epi64(input, 4);
        const __m256i high = _mm256_and_si256(input, clear_low);
        low = _mm256_shuffle_epi8(table_low, low);
        const __m256i product = _mm256_xor_si256(
            low, _mm256_shuffle_epi8(table_high, high));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(destination + offset), product);
    }
    return vector_bytes;
}

static GF256_AVX2_TARGET int gf256_muladd_mem_avx2_target(
    void* GF256_RESTRICT vz,
    const void* GF256_RESTRICT vx,
    uint8_t y,
    int bytes)
{
    const __m256i table_low = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_LO_Y + y));
    const __m256i table_high = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_HI_Y + y));
    const __m256i clear_low = _mm256_set1_epi8(0x0f);
    const int vector_bytes = bytes & ~31;
    const uint8_t* const source = reinterpret_cast<const uint8_t*>(vx);
    uint8_t* const destination = reinterpret_cast<uint8_t*>(vz);
    int offset = 0;
    for (; vector_bytes - offset >= 64; offset += 64)
    {
        __m256i input0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
        __m256i input1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset + 32));
        __m256i low0 = _mm256_and_si256(input0, clear_low);
        __m256i low1 = _mm256_and_si256(input1, clear_low);
        input0 = _mm256_srli_epi64(input0, 4);
        input1 = _mm256_srli_epi64(input1, 4);
        const __m256i high0 = _mm256_and_si256(input0, clear_low);
        const __m256i high1 = _mm256_and_si256(input1, clear_low);
        low0 = _mm256_shuffle_epi8(table_low, low0);
        low1 = _mm256_shuffle_epi8(table_low, low1);
        const __m256i product0 = _mm256_xor_si256(
            low0, _mm256_shuffle_epi8(table_high, high0));
        const __m256i product1 = _mm256_xor_si256(
            low1, _mm256_shuffle_epi8(table_high, high1));
        __m256i* const destination0 =
            reinterpret_cast<__m256i*>(destination + offset);
        __m256i* const destination1 =
            reinterpret_cast<__m256i*>(destination + offset + 32);
        _mm256_storeu_si256(
            destination0,
            _mm256_xor_si256(
                _mm256_loadu_si256(destination0), product0));
        _mm256_storeu_si256(
            destination1,
            _mm256_xor_si256(
                _mm256_loadu_si256(destination1), product1));
    }
    if (offset < vector_bytes)
    {
        __m256i input = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
        __m256i low = _mm256_and_si256(input, clear_low);
        input = _mm256_srli_epi64(input, 4);
        const __m256i high = _mm256_and_si256(input, clear_low);
        low = _mm256_shuffle_epi8(table_low, low);
        const __m256i product = _mm256_xor_si256(
            low, _mm256_shuffle_epi8(table_high, high));
        __m256i* const output =
            reinterpret_cast<__m256i*>(destination + offset);
        _mm256_storeu_si256(
            output,
            _mm256_xor_si256(_mm256_loadu_si256(output), product));
    }
    return vector_bytes;
}

static GF256_AVX2_TARGET int gf256_muladd_multi_mem_avx2_target(
    void* const* GF256_RESTRICT destinations,
    const uint8_t* GF256_RESTRICT scales,
    int destination_count,
    const uint8_t* GF256_RESTRICT source,
    int begin,
    int bytes)
{
    const __m256i clear_low = _mm256_set1_epi8(0x0f);
    const int vector_end = begin + ((bytes - begin) & ~31);
    for (int base = 0; base < destination_count; base += 4)
    {
        const int group_count = destination_count - base < 4 ?
            destination_count - base : 4;
        __m256i table_low[4];
        __m256i table_high[4];
        for (int local = 0; local < group_count; ++local)
        {
            const uint8_t scale = scales[base + local];
            table_low[local] = _mm256_broadcastsi128_si256(
                _mm_loadu_si128(GF256Ctx.MM128.TABLE_LO_Y + scale));
            table_high[local] = _mm256_broadcastsi128_si256(
                _mm_loadu_si128(GF256Ctx.MM128.TABLE_HI_Y + scale));
        }
        for (int offset = begin; offset < vector_end; offset += 32)
        {
            __m256i input = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(source + offset));
            const __m256i low = _mm256_and_si256(input, clear_low);
            input = _mm256_srli_epi64(input, 4);
            const __m256i high = _mm256_and_si256(input, clear_low);
            for (int local = 0; local < group_count; ++local)
            {
                if (scales[base + local] == 0u) continue;
                const __m256i product = _mm256_xor_si256(
                    _mm256_shuffle_epi8(table_low[local], low),
                    _mm256_shuffle_epi8(table_high[local], high));
                __m256i* const destination =
                    reinterpret_cast<__m256i*>(
                        reinterpret_cast<uint8_t*>(
                            destinations[base + local]) + offset);
                _mm256_storeu_si256(
                    destination,
                    _mm256_xor_si256(
                        _mm256_loadu_si256(destination), product));
            }
        }
    }
    return vector_end;
}

#endif

#if !defined(GF256_TARGET_MOBILE) && \
    (defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3))

// These helpers are individually target-qualified on baseline GCC/Clang x86
// builds.  Keeping every SSSE3 instruction in a function reached only after
// the CPUID check preserves the portable SSE2 fallback.
static GF256_SSSE3_TARGET int gf256_mul_mem_ssse3(
    void* vz, const void* vx, uint8_t y, int bytes)
{
    GF256_M128* z16 = reinterpret_cast<GF256_M128*>(vz);
    const GF256_M128* x16 = reinterpret_cast<const GF256_M128*>(vx);
    const GF256_M128 table_lo_y =
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_LO_Y + y);
    const GF256_M128 table_hi_y =
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_HI_Y + y);
    const GF256_M128 clear_low = _mm_set1_epi8(0x0f);
    const int vector_bytes = bytes & ~15;

    for (int offset = 0; offset < vector_bytes; offset += 16)
    {
        GF256_M128 input = _mm_loadu_si128(
            reinterpret_cast<const GF256_M128*>(
                reinterpret_cast<const uint8_t*>(x16) + offset));
        GF256_M128 low = _mm_and_si128(input, clear_low);
        input = _mm_srli_epi64(input, 4);
        const GF256_M128 high = _mm_and_si128(input, clear_low);
        low = _mm_shuffle_epi8(table_lo_y, low);
        const GF256_M128 product = _mm_xor_si128(
            low, _mm_shuffle_epi8(table_hi_y, high));
        _mm_storeu_si128(
            reinterpret_cast<GF256_M128*>(
                reinterpret_cast<uint8_t*>(z16) + offset),
            product);
    }
    return vector_bytes;
}

static GF256_SSSE3_TARGET int gf256_muladd_mem_ssse3(
    void* GF256_RESTRICT vz,
    const void* GF256_RESTRICT vx,
    uint8_t y,
    int bytes)
{
    GF256_M128* GF256_RESTRICT z16 =
        reinterpret_cast<GF256_M128*>(vz);
    const GF256_M128* GF256_RESTRICT x16 =
        reinterpret_cast<const GF256_M128*>(vx);
    const GF256_M128 table_lo_y =
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_LO_Y + y);
    const GF256_M128 table_hi_y =
        _mm_loadu_si128(GF256Ctx.MM128.TABLE_HI_Y + y);
    const GF256_M128 clear_low = _mm_set1_epi8(0x0f);
    const int vector_bytes = bytes & ~15;
    int offset = 0;

    // Two independent dependency chains improve throughput on the SSE path.
    for (; offset + 32 <= vector_bytes; offset += 32)
    {
        GF256_M128 input0 = _mm_loadu_si128(
            reinterpret_cast<const GF256_M128*>(
                reinterpret_cast<const uint8_t*>(x16) + offset));
        GF256_M128 input1 = _mm_loadu_si128(
            reinterpret_cast<const GF256_M128*>(
                reinterpret_cast<const uint8_t*>(x16) + offset + 16));
        GF256_M128 low0 = _mm_and_si128(input0, clear_low);
        GF256_M128 low1 = _mm_and_si128(input1, clear_low);
        input0 = _mm_srli_epi64(input0, 4);
        input1 = _mm_srli_epi64(input1, 4);
        const GF256_M128 high0 = _mm_and_si128(input0, clear_low);
        const GF256_M128 high1 = _mm_and_si128(input1, clear_low);
        low0 = _mm_shuffle_epi8(table_lo_y, low0);
        low1 = _mm_shuffle_epi8(table_lo_y, low1);
        const GF256_M128 product0 = _mm_xor_si128(
            low0, _mm_shuffle_epi8(table_hi_y, high0));
        const GF256_M128 product1 = _mm_xor_si128(
            low1, _mm_shuffle_epi8(table_hi_y, high1));
        GF256_M128* const destination0 =
            reinterpret_cast<GF256_M128*>(
                reinterpret_cast<uint8_t*>(z16) + offset);
        GF256_M128* const destination1 =
            reinterpret_cast<GF256_M128*>(
                reinterpret_cast<uint8_t*>(z16) + offset + 16);
        _mm_storeu_si128(
            destination0,
            _mm_xor_si128(_mm_loadu_si128(destination0), product0));
        _mm_storeu_si128(
            destination1,
            _mm_xor_si128(_mm_loadu_si128(destination1), product1));
    }
    if (offset < vector_bytes)
    {
        GF256_M128 input = _mm_loadu_si128(
            reinterpret_cast<const GF256_M128*>(
                reinterpret_cast<const uint8_t*>(x16) + offset));
        GF256_M128 low = _mm_and_si128(input, clear_low);
        input = _mm_srli_epi64(input, 4);
        const GF256_M128 high = _mm_and_si128(input, clear_low);
        low = _mm_shuffle_epi8(table_lo_y, low);
        const GF256_M128 product = _mm_xor_si128(
            low, _mm_shuffle_epi8(table_hi_y, high));
        GF256_M128* const destination =
            reinterpret_cast<GF256_M128*>(
                reinterpret_cast<uint8_t*>(z16) + offset);
        _mm_storeu_si128(
            destination,
            _mm_xor_si128(_mm_loadu_si128(destination), product));
    }
    return vector_bytes;
}

static GF256_SSSE3_TARGET int gf256_muladd_multi_mem_ssse3(
    void* const* GF256_RESTRICT destinations,
    const uint8_t* GF256_RESTRICT scales,
    int destination_count,
    const uint8_t* GF256_RESTRICT source,
    int begin,
    int bytes)
{
    const GF256_M128 clear_low = _mm_set1_epi8(0x0f);
    const int vector_end = begin + ((bytes - begin) & ~15);
    for (int base = 0; base < destination_count; base += 4)
    {
        const int group_count = destination_count - base < 4 ?
            destination_count - base : 4;
        GF256_M128 table_low[4];
        GF256_M128 table_high[4];
        for (int local = 0; local < group_count; ++local)
        {
            const uint8_t scale = scales[base + local];
            table_low[local] = _mm_loadu_si128(
                GF256Ctx.MM128.TABLE_LO_Y + scale);
            table_high[local] = _mm_loadu_si128(
                GF256Ctx.MM128.TABLE_HI_Y + scale);
        }
        for (int offset = begin; offset < vector_end; offset += 16)
        {
            GF256_M128 input = _mm_loadu_si128(
                reinterpret_cast<const GF256_M128*>(source + offset));
            const GF256_M128 low = _mm_and_si128(input, clear_low);
            input = _mm_srli_epi64(input, 4);
            const GF256_M128 high = _mm_and_si128(input, clear_low);
            for (int local = 0; local < group_count; ++local)
            {
                if (scales[base + local] == 0u) continue;
                const GF256_M128 product = _mm_xor_si128(
                    _mm_shuffle_epi8(table_low[local], low),
                    _mm_shuffle_epi8(table_high[local], high));
                GF256_M128* const destination =
                    reinterpret_cast<GF256_M128*>(
                        reinterpret_cast<uint8_t*>(
                            destinations[base + local]) + offset);
                _mm_storeu_si128(
                    destination,
                    _mm_xor_si128(
                        _mm_loadu_si128(destination), product));
            }
        }
    }
    return vector_end;
}

#endif

// vz == vx (in-place) is supported: no restrict qualifiers here.
extern "C" void gf256_mul_mem(void * vz, const void * vx, uint8_t y, int bytes)
{
    if (bytes <= 0) {
        return;
    }

    WH_BUMP(3, bytes);
    // Special cases are rare in codec hot paths.
    if (GF256_UNLIKELY(y <= 1))
    {
        if (y == 0)
            memset(vz, 0, bytes);
        else if (vz != vx)
            memcpy(vz, vx, bytes);
        return;
    }

    // Every pointer view in this function must preserve the documented
    // vz == vx contract.  In particular, do not add restrict to these
    // derived source/destination pointers: they are simultaneously live in
    // the supported in-place case.
    GF256_M128 * z16 = reinterpret_cast<GF256_M128 *>(vz);
    const GF256_M128 * x16 = reinterpret_cast<const GF256_M128 *>(vx);

#if defined(GF256_TARGET_MOBILE)
#if defined(GF256_TRY_NEON)
    if (bytes >= 16 && CpuHasNeon)
    {
        // Partial product tables; see above
        const GF256_M128 table_lo_y = vld1q_u8((uint8_t*)(GF256Ctx.MM128.TABLE_LO_Y + y));
        const GF256_M128 table_hi_y = vld1q_u8((uint8_t*)(GF256Ctx.MM128.TABLE_HI_Y + y));

        // clr_mask = 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
        const GF256_M128 clr_mask = vdupq_n_u8(0x0f);

        // Handle multiples of 16 bytes
        do
        {
            // See above comments for details
            GF256_M128 x0 = vld1q_u8((uint8_t*)x16);
            GF256_M128 l0 = vandq_u8(x0, clr_mask);
            x0 = vshrq_n_u8(x0, 4);
            GF256_M128 h0 = vandq_u8(x0, clr_mask);
            l0 = vqtbl1q_u8(table_lo_y, l0);
            h0 = vqtbl1q_u8(table_hi_y, h0);
            vst1q_u8((uint8_t*)z16, veorq_u8(l0, h0));

            bytes -= 16, ++x16, ++z16;
        } while (bytes >= 16);
    }
#endif
#else
# if defined(GF256_TRY_GFNI)
    // Ablation S1: one vgf2p8affineqb multiplies 64 bytes by the constant y under poly 0x14D.
    if (CpuHasGFNI && bytes >= 64)
    {
        const __m512i A = _mm512_set1_epi64((long long)GF256Ctx.GFNI_MUL_MATRIX[y]);
#  if defined(WH_GFNI_UNROLL4) && (WH_GFNI_UNROLL4+0)
        while (bytes >= 256)
        {
            const __m512i v0 = _mm512_loadu_si512((const void*)(x16 + 0));
            const __m512i v1 = _mm512_loadu_si512((const void*)(x16 + 4));
            const __m512i v2 = _mm512_loadu_si512((const void*)(x16 + 8));
            const __m512i v3 = _mm512_loadu_si512((const void*)(x16 + 12));
            _mm512_storeu_si512((void*)(z16 + 0), _mm512_gf2p8affine_epi64_epi8(v0, A, 0));
            _mm512_storeu_si512((void*)(z16 + 4), _mm512_gf2p8affine_epi64_epi8(v1, A, 0));
            _mm512_storeu_si512((void*)(z16 + 8), _mm512_gf2p8affine_epi64_epi8(v2, A, 0));
            _mm512_storeu_si512((void*)(z16 + 12), _mm512_gf2p8affine_epi64_epi8(v3, A, 0));
            bytes -= 256, x16 += 16, z16 += 16; // 256 bytes = 16 * M128
        }
#  endif
        while (bytes >= 64)
        {
            const __m512i v = _mm512_loadu_si512((const void*)x16);
            _mm512_storeu_si512((void*)z16, _mm512_gf2p8affine_epi64_epi8(v, A, 0));
            bytes -= 64, x16 += 4, z16 += 4; // 64 bytes = 4 * M128
        }
    }
# endif // GF256_TRY_GFNI
# if defined(GF256_TRY_AVX2)
    if (bytes >= 32 && CpuHasAVX2 && !GF256_GFNI_ACTIVE)
    {
        // Partial product tables; see above
        const GF256_M256 table_lo_y = _mm256_loadu_si256(GF256Ctx.MM256.TABLE_LO_Y + y);
        const GF256_M256 table_hi_y = _mm256_loadu_si256(GF256Ctx.MM256.TABLE_HI_Y + y);

        // clr_mask = 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
        const GF256_M256 clr_mask = _mm256_set1_epi8(0x0f);

        GF256_M256 * z32 = reinterpret_cast<GF256_M256 *>(vz);
        const GF256_M256 * x32 = reinterpret_cast<const GF256_M256 *>(vx);

        // Handle multiples of 32 bytes
        do
        {
            // See above comments for details
            GF256_M256 x0 = _mm256_loadu_si256(x32);
            GF256_M256 l0 = _mm256_and_si256(x0, clr_mask);
            x0 = _mm256_srli_epi64(x0, 4);
            GF256_M256 h0 = _mm256_and_si256(x0, clr_mask);
            l0 = _mm256_shuffle_epi8(table_lo_y, l0);
            h0 = _mm256_shuffle_epi8(table_hi_y, h0);
            _mm256_storeu_si256(z32, _mm256_xor_si256(l0, h0));

            bytes -= 32, ++x32, ++z32;
        } while (bytes >= 32);

        z16 = reinterpret_cast<GF256_M128 *>(z32);
        x16 = reinterpret_cast<const GF256_M128 *>(x32);
    }
# endif // GF256_TRY_AVX2
# if defined(GF256_TRY_TARGET_AVX2)
    if (bytes >= 32 && CpuHasAVX2 && !GF256_GFNI_ACTIVE)
    {
        const int vector_bytes =
            gf256_mul_mem_avx2_target(z16, x16, y, bytes);
        bytes -= vector_bytes;
        z16 += vector_bytes / 16;
        x16 += vector_bytes / 16;
    }
# endif // GF256_TRY_TARGET_AVX2
# if defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3)
    if (bytes >= 16 && CpuHasSSSE3)
    {
        const int vector_bytes = gf256_mul_mem_ssse3(z16, x16, y, bytes);
        bytes -= vector_bytes;
        z16 += vector_bytes / 16;
        x16 += vector_bytes / 16;
    }
# endif // GF256_TRY_SSSE3 || GF256_TRY_TARGET_SSSE3
#endif

    uint8_t * z1 = reinterpret_cast<uint8_t*>(z16);
    const uint8_t * x1 = reinterpret_cast<const uint8_t*>(x16);
    const uint8_t * GF256_RESTRICT table = GF256Ctx.GF256_MUL_TABLE + ((unsigned)y << 8);

    // Handle blocks of 8 bytes
    while (bytes >= 8)
    {
#ifdef GF256_IS_BIG_ENDIAN
        uint64_t word = (uint64_t)table[x1[0]] << 56;
        word |= (uint64_t)table[x1[1]] << 48;
        word |= (uint64_t)table[x1[2]] << 40;
        word |= (uint64_t)table[x1[3]] << 32;
        word |= (uint64_t)table[x1[4]] << 24;
        word |= (uint64_t)table[x1[5]] << 16;
        word |= (uint64_t)table[x1[6]] << 8;
        word |= (uint64_t)table[x1[7]];
#else
        uint64_t word = table[x1[0]];
        word |= (uint64_t)table[x1[1]] << 8;
        word |= (uint64_t)table[x1[2]] << 16;
        word |= (uint64_t)table[x1[3]] << 24;
        word |= (uint64_t)table[x1[4]] << 32;
        word |= (uint64_t)table[x1[5]] << 40;
        word |= (uint64_t)table[x1[6]] << 48;
        word |= (uint64_t)table[x1[7]] << 56;
#endif
        gf256_storeu64(z1, word);

        bytes -= 8, x1 += 8, z1 += 8;
    }

    // Handle a block of 4 bytes
    const int four = bytes & 4;
    if (four)
    {
#ifdef GF256_IS_BIG_ENDIAN
        uint32_t word = (uint32_t)table[x1[0]] << 24;
        word |= (uint32_t)table[x1[1]] << 16;
        word |= (uint32_t)table[x1[2]] << 8;
        word |= (uint32_t)table[x1[3]];
#else
        uint32_t word = table[x1[0]];
        word |= (uint32_t)table[x1[1]] << 8;
        word |= (uint32_t)table[x1[2]] << 16;
        word |= (uint32_t)table[x1[3]] << 24;
#endif
        gf256_storeu32(z1, word);
    }

    // Handle single bytes
    const int offset = four;
    switch (bytes & 3)
    {
    case 3: z1[offset + 2] = table[x1[offset + 2]]; GF256_FALLTHROUGH;
    case 2: z1[offset + 1] = table[x1[offset + 1]]; GF256_FALLTHROUGH;
    case 1: z1[offset] = table[x1[offset]];
    default:
        break;
    }
}

extern "C" void gf256_muladd_mem(void * GF256_RESTRICT vz, uint8_t y,
                                 const void * GF256_RESTRICT vx, int bytes)
{
    if (bytes <= 0) {
        return;
    }

    // Special cases are rare in codec hot paths.
    if (GF256_UNLIKELY(y <= 1))
    {
        // Count AFTER the special cases: y == 1 delegates to
        // gf256_add_mem, which does its own op-0 accounting (bumping
        // op-4 here too double-counted those bytes under WH_COUNT), and
        // y == 0 moves no bytes at all.
        if (y == 1)
            gf256_add_mem(vz, vx, bytes);
        return;
    }
    WH_BUMP(4, bytes);

    GF256_M128 * GF256_RESTRICT z16 = reinterpret_cast<GF256_M128 *>(vz);
    const GF256_M128 * GF256_RESTRICT x16 = reinterpret_cast<const GF256_M128 *>(vx);

#if defined(GF256_TARGET_MOBILE)
#if defined(GF256_TRY_NEON)
    if (bytes >= 16 && CpuHasNeon)
    {
        // Partial product tables; see above
        const GF256_M128 table_lo_y = vld1q_u8((uint8_t*)(GF256Ctx.MM128.TABLE_LO_Y + y));
        const GF256_M128 table_hi_y = vld1q_u8((uint8_t*)(GF256Ctx.MM128.TABLE_HI_Y + y));

        // clr_mask = 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
        const GF256_M128 clr_mask = vdupq_n_u8(0x0f);

        // Handle multiples of 16 bytes
        do
        {
            // See above comments for details
            GF256_M128 x0 = vld1q_u8((uint8_t*)x16);
            GF256_M128 l0 = vandq_u8(x0, clr_mask);

            // x0 = vshrq_n_u8(x0, 4);
            x0 = (GF256_M128)vshrq_n_u64( (uint64x2_t)x0, 4);
            GF256_M128 h0 = vandq_u8(x0, clr_mask);
            l0 = vqtbl1q_u8(table_lo_y, l0);
            h0 = vqtbl1q_u8(table_hi_y, h0);
            const GF256_M128 p0 = veorq_u8(l0, h0);
            const GF256_M128 z0 = vld1q_u8((uint8_t*)z16);
            vst1q_u8((uint8_t*)z16, veorq_u8(p0, z0));
            bytes -= 16, ++x16, ++z16;
        } while (bytes >= 16);
    }
#endif
#else // GF256_TARGET_MOBILE
# if defined(GF256_TRY_GFNI)
    // Ablation S1: z[] += x[] * y, 64 bytes per vgf2p8affineqb under poly 0x14D.
    if (CpuHasGFNI && bytes >= 64)
    {
        const __m512i A = _mm512_set1_epi64((long long)GF256Ctx.GFNI_MUL_MATRIX[y]);
#  if defined(WH_GFNI_UNROLL4) && (WH_GFNI_UNROLL4+0)
        while (bytes >= 256)
        {
            const __m512i v0 = _mm512_loadu_si512((const void*)(x16 + 0));
            const __m512i v1 = _mm512_loadu_si512((const void*)(x16 + 4));
            const __m512i v2 = _mm512_loadu_si512((const void*)(x16 + 8));
            const __m512i v3 = _mm512_loadu_si512((const void*)(x16 + 12));
            const __m512i z0 = _mm512_loadu_si512((const void*)(z16 + 0));
            const __m512i z1 = _mm512_loadu_si512((const void*)(z16 + 4));
            const __m512i z2 = _mm512_loadu_si512((const void*)(z16 + 8));
            const __m512i z3 = _mm512_loadu_si512((const void*)(z16 + 12));
            const __m512i p0 = _mm512_gf2p8affine_epi64_epi8(v0, A, 0);
            const __m512i p1 = _mm512_gf2p8affine_epi64_epi8(v1, A, 0);
            const __m512i p2 = _mm512_gf2p8affine_epi64_epi8(v2, A, 0);
            const __m512i p3 = _mm512_gf2p8affine_epi64_epi8(v3, A, 0);
            _mm512_storeu_si512((void*)(z16 + 0), _mm512_xor_si512(p0, z0));
            _mm512_storeu_si512((void*)(z16 + 4), _mm512_xor_si512(p1, z1));
            _mm512_storeu_si512((void*)(z16 + 8), _mm512_xor_si512(p2, z2));
            _mm512_storeu_si512((void*)(z16 + 12), _mm512_xor_si512(p3, z3));
            bytes -= 256, x16 += 16, z16 += 16; // 256 bytes = 16 * M128
        }
#  endif
        while (bytes >= 64)
        {
            const __m512i v = _mm512_loadu_si512((const void*)x16);
            const __m512i p = _mm512_gf2p8affine_epi64_epi8(v, A, 0);
            const __m512i zz = _mm512_loadu_si512((const void*)z16);
            _mm512_storeu_si512((void*)z16, _mm512_xor_si512(p, zz));
            bytes -= 64, x16 += 4, z16 += 4; // 64 bytes = 4 * M128
        }
    }
# endif // GF256_TRY_GFNI
# if defined(GF256_TRY_AVX2)
    if (bytes >= 32 && CpuHasAVX2 && !GF256_GFNI_ACTIVE)
    {
        // Partial product tables; see above
        const GF256_M256 table_lo_y = _mm256_loadu_si256(GF256Ctx.MM256.TABLE_LO_Y + y);
        const GF256_M256 table_hi_y = _mm256_loadu_si256(GF256Ctx.MM256.TABLE_HI_Y + y);

        // clr_mask = 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
        const GF256_M256 clr_mask = _mm256_set1_epi8(0x0f);

        GF256_M256 * GF256_RESTRICT z32 = reinterpret_cast<GF256_M256 *>(z16);
        const GF256_M256 * GF256_RESTRICT x32 = reinterpret_cast<const GF256_M256 *>(x16);

        // On my Reed Solomon codec, the encoder unit test runs in 640 usec without and 550 usec with the optimization (86% of the original time)
        const unsigned count = bytes / 64;
        for (unsigned i = 0; i < count; ++i)
        {
            // See above comments for details
            GF256_M256 x0 = _mm256_loadu_si256(x32 + i * 2);
            GF256_M256 l0 = _mm256_and_si256(x0, clr_mask);
            x0 = _mm256_srli_epi64(x0, 4);
            const GF256_M256 z0 = _mm256_loadu_si256(z32 + i * 2);
            GF256_M256 h0 = _mm256_and_si256(x0, clr_mask);
            l0 = _mm256_shuffle_epi8(table_lo_y, l0);
            h0 = _mm256_shuffle_epi8(table_hi_y, h0);
            const GF256_M256 p0 = _mm256_xor_si256(l0, h0);
            _mm256_storeu_si256(z32 + i * 2, _mm256_xor_si256(p0, z0));

            GF256_M256 x1 = _mm256_loadu_si256(x32 + i * 2 + 1);
            GF256_M256 l1 = _mm256_and_si256(x1, clr_mask);
            x1 = _mm256_srli_epi64(x1, 4);
            const GF256_M256 z1 = _mm256_loadu_si256(z32 + i * 2 + 1);
            GF256_M256 h1 = _mm256_and_si256(x1, clr_mask);
            l1 = _mm256_shuffle_epi8(table_lo_y, l1);
            h1 = _mm256_shuffle_epi8(table_hi_y, h1);
            const GF256_M256 p1 = _mm256_xor_si256(l1, h1);
            _mm256_storeu_si256(z32 + i * 2 + 1, _mm256_xor_si256(p1, z1));
        }
        bytes -= count * 64;
        z32 += count * 2;
        x32 += count * 2;

        if (bytes >= 32)
        {
            GF256_M256 x0 = _mm256_loadu_si256(x32);
            GF256_M256 l0 = _mm256_and_si256(x0, clr_mask);
            x0 = _mm256_srli_epi64(x0, 4);
            GF256_M256 h0 = _mm256_and_si256(x0, clr_mask);
            l0 = _mm256_shuffle_epi8(table_lo_y, l0);
            h0 = _mm256_shuffle_epi8(table_hi_y, h0);
            const GF256_M256 p0 = _mm256_xor_si256(l0, h0);
            const GF256_M256 z0 = _mm256_loadu_si256(z32);
            _mm256_storeu_si256(z32, _mm256_xor_si256(p0, z0));

            bytes -= 32;
            z32++;
            x32++;
        }

        z16 = reinterpret_cast<GF256_M128 *>(z32);
        x16 = reinterpret_cast<const GF256_M128 *>(x32);
    }
# endif // GF256_TRY_AVX2
# if defined(GF256_TRY_TARGET_AVX2)
    if (bytes >= 32 && CpuHasAVX2 && !GF256_GFNI_ACTIVE)
    {
        const int vector_bytes =
            gf256_muladd_mem_avx2_target(z16, x16, y, bytes);
        bytes -= vector_bytes;
        z16 += vector_bytes / 16;
        x16 += vector_bytes / 16;
    }
# endif // GF256_TRY_TARGET_AVX2
# if defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3)
    if (bytes >= 16 && CpuHasSSSE3)
    {
        const int vector_bytes =
            gf256_muladd_mem_ssse3(z16, x16, y, bytes);
        bytes -= vector_bytes;
        z16 += vector_bytes / 16;
        x16 += vector_bytes / 16;
    }
# endif // GF256_TRY_SSSE3 || GF256_TRY_TARGET_SSSE3
#endif // GF256_TARGET_MOBILE

    uint8_t * GF256_RESTRICT z1 = reinterpret_cast<uint8_t*>(z16);
    const uint8_t * GF256_RESTRICT x1 = reinterpret_cast<const uint8_t*>(x16);
    const uint8_t * GF256_RESTRICT table = GF256Ctx.GF256_MUL_TABLE + ((unsigned)y << 8);

    // Handle blocks of 8 bytes
    while (bytes >= 8)
    {
#ifdef GF256_IS_BIG_ENDIAN
        uint64_t word = (uint64_t)table[x1[0]] << 56;
        word |= (uint64_t)table[x1[1]] << 48;
        word |= (uint64_t)table[x1[2]] << 40;
        word |= (uint64_t)table[x1[3]] << 32;
        word |= (uint64_t)table[x1[4]] << 24;
        word |= (uint64_t)table[x1[5]] << 16;
        word |= (uint64_t)table[x1[6]] << 8;
        word |= (uint64_t)table[x1[7]];
#else
        uint64_t word = table[x1[0]];
        word |= (uint64_t)table[x1[1]] << 8;
        word |= (uint64_t)table[x1[2]] << 16;
        word |= (uint64_t)table[x1[3]] << 24;
        word |= (uint64_t)table[x1[4]] << 32;
        word |= (uint64_t)table[x1[5]] << 40;
        word |= (uint64_t)table[x1[6]] << 48;
        word |= (uint64_t)table[x1[7]] << 56;
#endif
        gf256_storeu64(z1, gf256_loadu64(z1) ^ word);

        bytes -= 8, x1 += 8, z1 += 8;
    }

    // Handle a block of 4 bytes
    const int four = bytes & 4;
    if (four)
    {
#ifdef GF256_IS_BIG_ENDIAN
        uint32_t word = (uint32_t)table[x1[0]] << 24;
        word |= (uint32_t)table[x1[1]] << 16;
        word |= (uint32_t)table[x1[2]] << 8;
        word |= (uint32_t)table[x1[3]];
#else
        uint32_t word = table[x1[0]];
        word |= (uint32_t)table[x1[1]] << 8;
        word |= (uint32_t)table[x1[2]] << 16;
        word |= (uint32_t)table[x1[3]] << 24;
#endif
        gf256_storeu32(z1, gf256_loadu32(z1) ^ word);
    }

    // Handle single bytes
    const int offset = four;
    switch (bytes & 3)
    {
    case 3: z1[offset + 2] ^= table[x1[offset + 2]]; GF256_FALLTHROUGH;
    case 2: z1[offset + 1] ^= table[x1[offset + 1]]; GF256_FALLTHROUGH;
    case 1: z1[offset] ^= table[x1[offset]];
    default:
        break;
    }
}

extern "C" void gf256_muladd_multi_mem(
    void * const * GF256_RESTRICT destinations,
    const uint8_t * GF256_RESTRICT scales,
    int destination_count,
    const void * GF256_RESTRICT source,
    int bytes)
{
    if (bytes <= 0 || destination_count <= 0) {
        return;
    }

#if defined(GF256_TRY_GFNI)
    // GFNI already reduces each constant multiply to one instruction.  On
    // current x86 hardware, keeping its tuned per-destination unrolling is
    // faster than a source-fused loop with many simultaneously live affine
    // matrices, so dispatch without adding duplicate instrumentation here.
    if (CpuHasGFNI)
    {
        for (int j = 0; j < destination_count; ++j) {
            gf256_muladd_mem(
                destinations[j], scales[j], source, bytes);
        }
        return;
    }
#endif

    for (int j = 0; j < destination_count; ++j)
    {
        if (scales[j] == 1u) {
            WH_BUMP(0, bytes);
        }
        else if (scales[j] > 1u) {
            WH_BUMP(4, bytes);
        }
    }

    const uint8_t* const x = reinterpret_cast<const uint8_t*>(source);
    int offset = 0;

#if defined(GF256_TARGET_MOBILE)
# if defined(GF256_TRY_NEON)
    if (CpuHasNeon)
    {
        const GF256_M128 clear_low = vdupq_n_u8(0x0f);
        const int vector_end = bytes & ~15;
        for (int base = 0; base < destination_count; base += 4)
        {
            const int group_count =
                destination_count - base < 4 ? destination_count - base : 4;
            GF256_M128 table_low[4];
            GF256_M128 table_high[4];
            for (int local = 0; local < group_count; ++local)
            {
                const uint8_t scale = scales[base + local];
                table_low[local] = vld1q_u8(
                    (const uint8_t*)(GF256Ctx.MM128.TABLE_LO_Y + scale));
                table_high[local] = vld1q_u8(
                    (const uint8_t*)(GF256Ctx.MM128.TABLE_HI_Y + scale));
            }
            for (int vector_offset = 0;
                 vector_offset < vector_end;
                 vector_offset += 16)
            {
                GF256_M128 input = vld1q_u8(x + vector_offset);
                const GF256_M128 low = vandq_u8(input, clear_low);
                input = (GF256_M128)vshrq_n_u64((uint64x2_t)input, 4);
                const GF256_M128 high = vandq_u8(input, clear_low);
                for (int local = 0; local < group_count; ++local)
                {
                    if (scales[base + local] == 0u) {
                        continue;
                    }
                    const GF256_M128 product = veorq_u8(
                        vqtbl1q_u8(table_low[local], low),
                        vqtbl1q_u8(table_high[local], high));
                    uint8_t* const destination =
                        reinterpret_cast<uint8_t*>(
                            destinations[base + local]) + vector_offset;
                    vst1q_u8(
                        destination,
                        veorq_u8(vld1q_u8(destination), product));
                }
            }
        }
        offset = vector_end;
    }
# endif
#else
# if defined(GF256_TRY_AVX2)
    if (CpuHasAVX2 && !GF256_GFNI_ACTIVE)
    {
        const GF256_M256 clear_low = _mm256_set1_epi8(0x0f);
        const int vector_end = bytes & ~31;
        for (int base = 0; base < destination_count; base += 4)
        {
            const int group_count =
                destination_count - base < 4 ? destination_count - base : 4;
            GF256_M256 table_low[4];
            GF256_M256 table_high[4];
            for (int local = 0; local < group_count; ++local)
            {
                const uint8_t scale = scales[base + local];
                table_low[local] = _mm256_loadu_si256(
                    GF256Ctx.MM256.TABLE_LO_Y + scale);
                table_high[local] = _mm256_loadu_si256(
                    GF256Ctx.MM256.TABLE_HI_Y + scale);
            }
            for (int vector_offset = 0;
                 vector_offset < vector_end;
                 vector_offset += 32)
            {
                GF256_M256 input = _mm256_loadu_si256(
                    reinterpret_cast<const GF256_M256*>(x + vector_offset));
                const GF256_M256 low = _mm256_and_si256(input, clear_low);
                input = _mm256_srli_epi64(input, 4);
                const GF256_M256 high = _mm256_and_si256(input, clear_low);
                for (int local = 0; local < group_count; ++local)
                {
                    if (scales[base + local] == 0u) {
                        continue;
                    }
                    const GF256_M256 product = _mm256_xor_si256(
                        _mm256_shuffle_epi8(table_low[local], low),
                        _mm256_shuffle_epi8(table_high[local], high));
                    GF256_M256* const destination =
                        reinterpret_cast<GF256_M256*>(
                            reinterpret_cast<uint8_t*>(
                                destinations[base + local]) + vector_offset);
                    _mm256_storeu_si256(
                        destination,
                        _mm256_xor_si256(
                            _mm256_loadu_si256(destination), product));
                }
            }
        }
        offset = vector_end;
    }
# endif
# if defined(GF256_TRY_TARGET_AVX2)
    if (CpuHasAVX2 && !GF256_GFNI_ACTIVE)
    {
        offset = gf256_muladd_multi_mem_avx2_target(
            destinations, scales, destination_count, x, offset, bytes);
    }
# endif // GF256_TRY_TARGET_AVX2
# if defined(GF256_TRY_SSSE3) || defined(GF256_TRY_TARGET_SSSE3)
    if (CpuHasSSSE3)
    {
        offset = gf256_muladd_multi_mem_ssse3(
            destinations, scales, destination_count, x, offset, bytes);
    }
# endif // GF256_TRY_SSSE3 || GF256_TRY_TARGET_SSSE3
#endif

    for (int base = 0; base < destination_count; base += 16)
    {
        const int group_count = destination_count - base < 16 ?
            destination_count - base : 16;
        const uint8_t* tables[16];
        for (int local = 0; local < group_count; ++local) {
            tables[local] = GF256Ctx.GF256_MUL_TABLE +
                ((unsigned)scales[base + local] << 8);
        }

        int scalar_offset = offset;
        while (bytes - scalar_offset >= 8)
        {
            const uint8_t* const input = x + scalar_offset;
            for (int local = 0; local < group_count; ++local)
            {
                if (scales[base + local] == 0u) {
                    continue;
                }
                const uint8_t* const table = tables[local];
#ifdef GF256_IS_BIG_ENDIAN
                uint64_t product = (uint64_t)table[input[0]] << 56;
                product |= (uint64_t)table[input[1]] << 48;
                product |= (uint64_t)table[input[2]] << 40;
                product |= (uint64_t)table[input[3]] << 32;
                product |= (uint64_t)table[input[4]] << 24;
                product |= (uint64_t)table[input[5]] << 16;
                product |= (uint64_t)table[input[6]] << 8;
                product |= (uint64_t)table[input[7]];
#else
                uint64_t product = table[input[0]];
                product |= (uint64_t)table[input[1]] << 8;
                product |= (uint64_t)table[input[2]] << 16;
                product |= (uint64_t)table[input[3]] << 24;
                product |= (uint64_t)table[input[4]] << 32;
                product |= (uint64_t)table[input[5]] << 40;
                product |= (uint64_t)table[input[6]] << 48;
                product |= (uint64_t)table[input[7]] << 56;
#endif
                uint8_t* const destination =
                    reinterpret_cast<uint8_t*>(
                        destinations[base + local]) + scalar_offset;
                gf256_storeu64(
                    destination,
                    gf256_loadu64(destination) ^ product);
            }
            scalar_offset += 8;
        }
        for (; scalar_offset < bytes; ++scalar_offset)
        {
            const uint8_t value = x[scalar_offset];
            for (int local = 0; local < group_count; ++local) {
                if (scales[base + local] != 0u) {
                    reinterpret_cast<uint8_t*>(
                        destinations[base + local])[scalar_offset] ^=
                        tables[local][value];
                }
            }
        }
    }
}

extern "C" void gf256_memswap(void * GF256_RESTRICT vx, void * GF256_RESTRICT vy, int bytes)
{
    if (bytes <= 0) {
        return;
    }

    WH_BUMP(5, bytes);
#if defined(GF256_TARGET_MOBILE)
    uint8_t * GF256_RESTRICT x16 = reinterpret_cast<uint8_t *>(vx);
    uint8_t * GF256_RESTRICT y16 = reinterpret_cast<uint8_t *>(vy);

    const unsigned count = (unsigned)bytes / 8;
    for (unsigned ii = 0; ii < count; ++ii)
    {
        const unsigned offset = ii * 8u;
        const uint64_t temp = gf256_loadu64(x16 + offset);
        gf256_storeu64(x16 + offset, gf256_loadu64(y16 + offset));
        gf256_storeu64(y16 + offset, temp);
    }

    x16 += count * 8u;
    y16 += count * 8u;
    bytes -= count * 8;
#else
    GF256_M128 * GF256_RESTRICT x16 = reinterpret_cast<GF256_M128 *>(vx);
    GF256_M128 * GF256_RESTRICT y16 = reinterpret_cast<GF256_M128 *>(vy);

# if defined(GF256_TRY_AVX512)
    if (CpuHasAVX512)
    {
        while (bytes >= 256)
        {
            const __m512i x0 = _mm512_loadu_si512((const void*)(x16 + 0));
            const __m512i x1 = _mm512_loadu_si512((const void*)(x16 + 4));
            const __m512i x2 = _mm512_loadu_si512((const void*)(x16 + 8));
            const __m512i x3 = _mm512_loadu_si512((const void*)(x16 + 12));
            const __m512i y0 = _mm512_loadu_si512((const void*)(y16 + 0));
            const __m512i y1 = _mm512_loadu_si512((const void*)(y16 + 4));
            const __m512i y2 = _mm512_loadu_si512((const void*)(y16 + 8));
            const __m512i y3 = _mm512_loadu_si512((const void*)(y16 + 12));
            _mm512_storeu_si512((void*)(x16 + 0), y0);
            _mm512_storeu_si512((void*)(x16 + 4), y1);
            _mm512_storeu_si512((void*)(x16 + 8), y2);
            _mm512_storeu_si512((void*)(x16 + 12), y3);
            _mm512_storeu_si512((void*)(y16 + 0), x0);
            _mm512_storeu_si512((void*)(y16 + 4), x1);
            _mm512_storeu_si512((void*)(y16 + 8), x2);
            _mm512_storeu_si512((void*)(y16 + 12), x3);
            bytes -= 256, x16 += 16, y16 += 16;
        }
        while (bytes >= 64)
        {
            const __m512i x0 = _mm512_loadu_si512((const void*)x16);
            const __m512i y0 = _mm512_loadu_si512((const void*)y16);
            _mm512_storeu_si512((void*)x16, y0);
            _mm512_storeu_si512((void*)y16, x0);
            bytes -= 64, x16 += 4, y16 += 4;
        }
    }
# endif // GF256_TRY_AVX512

# if defined(GF256_TRY_AVX2)
    if (CpuHasAVX2 && !GF256_AVX512_ACTIVE)
    {
        GF256_M256 * GF256_RESTRICT x32 = reinterpret_cast<GF256_M256 *>(x16);
        GF256_M256 * GF256_RESTRICT y32 = reinterpret_cast<GF256_M256 *>(y16);

        while (bytes >= 128)
        {
            const GF256_M256 x0 = _mm256_loadu_si256(x32 + 0);
            const GF256_M256 x1 = _mm256_loadu_si256(x32 + 1);
            const GF256_M256 x2 = _mm256_loadu_si256(x32 + 2);
            const GF256_M256 x3 = _mm256_loadu_si256(x32 + 3);
            const GF256_M256 y0 = _mm256_loadu_si256(y32 + 0);
            const GF256_M256 y1 = _mm256_loadu_si256(y32 + 1);
            const GF256_M256 y2 = _mm256_loadu_si256(y32 + 2);
            const GF256_M256 y3 = _mm256_loadu_si256(y32 + 3);
            _mm256_storeu_si256(x32 + 0, y0);
            _mm256_storeu_si256(x32 + 1, y1);
            _mm256_storeu_si256(x32 + 2, y2);
            _mm256_storeu_si256(x32 + 3, y3);
            _mm256_storeu_si256(y32 + 0, x0);
            _mm256_storeu_si256(y32 + 1, x1);
            _mm256_storeu_si256(y32 + 2, x2);
            _mm256_storeu_si256(y32 + 3, x3);
            bytes -= 128, x32 += 4, y32 += 4;
        }
        while (bytes >= 32)
        {
            const GF256_M256 x0 = _mm256_loadu_si256(x32);
            const GF256_M256 y0 = _mm256_loadu_si256(y32);
            _mm256_storeu_si256(x32, y0);
            _mm256_storeu_si256(y32, x0);
            bytes -= 32, ++x32, ++y32;
        }

        x16 = reinterpret_cast<GF256_M128 *>(x32);
        y16 = reinterpret_cast<GF256_M128 *>(y32);
    }
# endif // GF256_TRY_AVX2

    while (bytes >= 16)
    {
        GF256_M128 x0 = _mm_loadu_si128(x16);
        GF256_M128 y0 = _mm_loadu_si128(y16);
        _mm_storeu_si128(x16, y0);
        _mm_storeu_si128(y16, x0);

        bytes -= 16, ++x16, ++y16;
    }
#endif

    uint8_t * GF256_RESTRICT x1 = reinterpret_cast<uint8_t *>(x16);
    uint8_t * GF256_RESTRICT y1 = reinterpret_cast<uint8_t *>(y16);

    // Handle a block of 8 bytes
    const int eight = bytes & 8;
    if (eight)
    {
        const uint64_t temp = gf256_loadu64(x1);
        gf256_storeu64(x1, gf256_loadu64(y1));
        gf256_storeu64(y1, temp);
    }

    // Handle a block of 4 bytes
    const int four = bytes & 4;
    if (four)
    {
        const uint32_t temp = gf256_loadu32(x1 + eight);
        gf256_storeu32(x1 + eight, gf256_loadu32(y1 + eight));
        gf256_storeu32(y1 + eight, temp);
    }

    // Handle final bytes
    const int offset = eight + four;
    uint8_t temp;
    switch (bytes & 3)
    {
    case 3: temp = x1[offset + 2]; x1[offset + 2] = y1[offset + 2]; y1[offset + 2] = temp; GF256_FALLTHROUGH;
    case 2: temp = x1[offset + 1]; x1[offset + 1] = y1[offset + 1]; y1[offset + 1] = temp; GF256_FALLTHROUGH;
    case 1: temp = x1[offset]; x1[offset] = y1[offset]; y1[offset] = temp;
    default:
        break;
    }
}
