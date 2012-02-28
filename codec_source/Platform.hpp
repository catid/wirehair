/*
	Copyright (c) 2009-2012 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of WirehairFEC nor the names of its contributors may be
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

#ifndef CAT_PLATFORM_HPP
#define CAT_PLATFORM_HPP

// This definition overrides CAT_BUILD_DLL below.  Neuters CAT_EXPORT macro so symbols are
// neither exported or imported.
#define CAT_NEUTER_EXPORT

// This definition changes the meaning of the CAT_EXPORT macro on Windows.  When defined,
// the CAT_EXPORT macro will export the associated symbol.  When undefined, it will import it.
#define CAT_BUILD_DLL

namespace cat {


// Used to detect whether or not WirehairFEC is included
#define CAT_HAS_TOYS


//// Compiler ////

// Mac OS X additional compilation flags
#if defined(__APPLE__)
# include <TargetConditionals.h>
#endif

#if defined(__MINGW32__)
# define CAT_COMPILER_MINGW
#endif

//-----------------------------------------------------------------------------
// Intel C++ Compiler : Interoperates with MSVC and GCC
#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
# define CAT_COMPILER_ICC
# define CAT_FENCE_COMPILER __memory_barrier();
#endif

//-----------------------------------------------------------------------------
// Borland C++ Compiler : Compatible with MSVC syntax
#if defined(__BORLANDC__)
# define CAT_COMPILER_BORLAND
# define CAT_COMPILER_COMPAT_MSVC
# define CAT_INLINE __inline
# define CAT_ASM_EMIT __emit__

//-----------------------------------------------------------------------------
// Digital Mars C++ Compiler (previously known as Symantec C++)
#elif defined(__DMC__) || defined(__SC__) || defined(__SYMANTECC__)
# define CAT_COMPILER_DMARS
# define CAT_COMPILER_COMPAT_MSVC
# define CAT_INLINE __inline
# define CAT_ASM_EMIT __emit__

//-----------------------------------------------------------------------------
// Codeplay VectorC C++ Compiler : Compatible with GCC and MSVC syntax, prefer GCC
#elif defined(__VECTORC__)
# define CAT_COMPILER_CODEPLAY
# define CAT_COMPILER_COMPAT_GCC

//-----------------------------------------------------------------------------
// Pathscale C++ Compiler : Compatible with GCC syntax
#elif defined(__PATHSCALE__)
# define CAT_COMPILER_PATHSCALE
# define CAT_COMPILER_COMPAT_GCC

//-----------------------------------------------------------------------------
// Watcom C++ Compiler : Compatible with GCC and MSVC syntax, prefer GCC
#elif defined(__WATCOMC__)
# define CAT_COMPILER_WATCOM
# define CAT_COMPILER_COMPAT_GCC

//-----------------------------------------------------------------------------
// SUN C++ Compiler : Compatible with GCC syntax
#elif defined(__SUNPRO_CC)
# define CAT_COMPILER_SUN
# define CAT_COMPILER_COMPAT_GCC

//-----------------------------------------------------------------------------
// Metrowerks C++ Compiler : Compatible with MSVC syntax
#elif defined(__MWERKS__)
# define CAT_COMPILER_MWERKS
# define CAT_COMPILER_COMPAT_MSVC
# define CAT_INLINE inline
# define CAT_ASM_BEGIN _asm {
# define CAT_ASM_EMIT __emit__

//-----------------------------------------------------------------------------
// GNU C++ Compiler
// SN Systems ProDG C++ Compiler : Compatible with GCC
#elif defined(__GNUC__) || defined(__APPLE_CC__) || defined(__SNC__)
# define CAT_COMPILER_GCC
# define CAT_COMPILER_COMPAT_GCC
# define CAT_FASTCALL __attribute__ ((fastcall))

//-----------------------------------------------------------------------------
// Microsoft Visual Studio C++ Compiler
#elif defined(_MSC_VER)
# define CAT_COMPILER_MSVC
# define CAT_COMPILER_COMPAT_MSVC
# define CAT_FASTCALL __fastcall

# if _MSC_VER >= 1300
#  define CAT_FUNCTION __FUNCTION__
# endif

} // namespace cat
# include <cstdlib> // Intrinsics
# include <intrin.h> // Intrinsics
namespace cat {

//-----------------------------------------------------------------------------
// Otherwise unknown compiler
#else
# define CAT_COMPILER_UNKNOWN
# define CAT_ALIGNED(n) /* no way to detect alignment syntax */
# define CAT_PACKED /* no way to detect packing syntax */
# define CAT_INLINE inline
// No way to support inline assembly code here
# define CAT_RESTRICT

#endif

/*
	A lot of compilers have similar syntax to MSVC or GCC,
	so for simplicity I have those two defined below, and
	any deviations are implemented with overrides above.
*/

// MSVC-compatible compilers
#if defined(CAT_COMPILER_COMPAT_MSVC)

#if !defined(CAT_ALIGNED)
# define CAT_ALIGNED(n) __declspec(align(n))
#endif
#if !defined(CAT_PACKED)
# define CAT_PACKED
# define CAT_PRAGMA_PACK
#endif
#if !defined(CAT_INLINE)
# define CAT_INLINE __forceinline
#endif
#if !defined(CAT_ASM_INTEL)
# define CAT_ASM_INTEL
#endif
#if !defined(CAT_ASM_BEGIN)
# define CAT_ASM_BEGIN __asm {
#endif
#if !defined(CAT_ASM_EMIT)
# define CAT_ASM_EMIT _emit
#endif
#if !defined(CAT_ASM_END)
# define CAT_ASM_END }
#endif
#if !defined(CAT_TLS)
# define CAT_TLS __declspec( thread )
#endif
#if !defined(CAT_RESTRICT)
# define CAT_RESTRICT __restrict
#endif
#if !defined(CAT_FENCE_COMPILER)
# if defined(CAT_COMPILER_MSVC)
#  pragma intrinsic(_ReadWriteBarrier)
# endif
# define CAT_FENCE_COMPILER _ReadWriteBarrier();
#endif
#if !defined(CAT_DLL_EXPORT)
# define CAT_DLL_EXPORT __declspec(dllexport)
#endif
#if !defined(CAT_DLL_IMPORT)
# define CAT_DLL_IMPORT __declspec(dllimport)
#endif
#if !defined(CAT_FUNCTION)
# define CAT_FUNCTION "???"
#endif

// GCC-compatible compilers
#elif defined(CAT_COMPILER_COMPAT_GCC)

#if !defined(CAT_ALIGNED)
# define CAT_ALIGNED(n) __attribute__ ((aligned (n)))
#endif
#if !defined(CAT_PACKED)
# define CAT_PACKED __attribute__ ((packed))
#endif
#if !defined(CAT_INLINE)
# define CAT_INLINE inline /* __inline__ __attribute__((always_inline)) */
#endif
#if !defined(CAT_ASM_ATT)
# define CAT_ASM_ATT
#endif
#if !defined(CAT_ASM_BEGIN)
# define CAT_ASM_BEGIN __asm__ __volatile__ (
#endif
#if !defined(CAT_ASM_EMIT)
# define CAT_ASM_EMIT .byte
#endif
#if !defined(CAT_ASM_END)
# define CAT_ASM_END );
#endif
#if !defined(CAT_TLS)
# define CAT_TLS __thread
#endif
#if !defined(CAT_RESTRICT)
# define CAT_RESTRICT __restrict__
#endif
#if !defined(CAT_FENCE_COMPILER)
# define CAT_FENCE_COMPILER CAT_ASM_BEGIN "" ::: "memory" CAT_ASM_END
#endif
#if !defined(CAT_DLL_EXPORT)
# define CAT_DLL_EXPORT __attribute__((dllexport))
#endif
#if !defined(CAT_DLL_IMPORT)
# define CAT_DLL_IMPORT __attribute__((dllimport))
#endif
#if !defined(CAT_FUNCTION)
# define CAT_FUNCTION "???"
#endif

#endif // CAT_COMPILER_COMPAT_*


//// Debug Flag ////

#if defined(CAT_COMPILER_MSVC)

# if defined(_DEBUG)
#  define CAT_DEBUG
# endif

#else

# if !defined(NDEBUG)
#  define CAT_DEBUG
# endif

#endif


#if defined(CAT_DEBUG)
#define CAT_IF_DEBUG(x) x
#else
#define CAT_IF_DEBUG(x)
#endif


//// Instruction Set Architecture ////

#if defined(__powerpc__) || defined(__ppc__) || defined(_POWER) || defined(_M_PPC) || \
	defined(_M_MPPC) || defined(__POWERPC) || defined(powerpc) || defined(__ppc64__) || \
	defined(_PS3) || defined(__PS3__) || defined(SN_TARGET_PS3) || defined(__POWERPC__)
# define CAT_ISA_PPC

#elif defined(__i386__) || defined(i386) || defined(intel) || defined(_M_IX86) || \
	  defined(__ia64) || defined(__ia64__) || defined(__x86_64) || defined(_M_IA64) || \
	  defined(_M_X64) || defined(_M_I86) || defined(sun386) || defined(__OS2__)
# define CAT_ISA_X86

#elif defined(TARGET_CPU_ARM)
# define CAT_ISA_ARM

#elif defined(__mips__)
# define CAT_ISA_MIPS

#elif defined(__alpha__) || defined(__alpha)
# define CAT_ISA_ALPHA

#elif defined(vax) || defined(vms) || defined(VMS) || defined(__VMS)
# define CAT_ISA_VMS

#elif defined(AVR)
# define CAT_ISA_ATMEGA

#elif defined(AMIGA)
# define CAT_ISA_AMIGA

#else
# define CAT_ISA_UNKNOWN
#endif


//// Memory Access Alignment ////

#if !defined(CAT_ISA_PPC) && !defined(CAT_ISA_X86)
# define CAT_NEED_ALIGNED_READS
#endif


//// Endianness ////

// If no override for endianness from Config.hpp,
#if !defined(CAT_ENDIAN_BIG) && !defined(CAT_ENDIAN_LITTLE)

#if defined(BIG_ENDIAN) && defined(LITTLE_ENDIAN)
# if defined(BYTE_ORDER) && BYTE_ORDER == BIG_ENDIAN
#  define CAT_ENDIAN_BIG
# elif defined(BYTE_ORDER) && BYTE_ORDER == LITTLE_ENDIAN
#  define CAT_ENDIAN_LITTLE
# endif
#elif defined(BIG_ENDIAN)
# define CAT_ENDIAN_BIG
#elif defined(LITTLE_ENDIAN)
# define CAT_ENDIAN_LITTLE
#endif

#if defined(_BIG_ENDIAN) && defined(_LITTLE_ENDIAN)
# if defined(_BYTE_ORDER) && _BYTE_ORDER == _BIG_ENDIAN
#  define CAT_ENDIAN_BIG
# elif defined(_BYTE_ORDER) && _BYTE_ORDER == _LITTLE_ENDIAN
#  define CAT_ENDIAN_LITTLE
# endif
#elif defined(_BIG_ENDIAN)
# define CAT_ENDIAN_BIG
#elif defined(_LITTLE_ENDIAN)
# define CAT_ENDIAN_LITTLE
#endif

#if defined(__BIG_ENDIAN) && defined(__LITTLE_ENDIAN)
# if defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN
#  define CAT_ENDIAN_BIG
# elif defined(__BYTE_ORDER) && __BYTE_ORDER == __LITTLE_ENDIAN
#  define CAT_ENDIAN_LITTLE
# endif
#elif defined(__BIG_ENDIAN)
# define CAT_ENDIAN_BIG
#elif defined(__LITTLE_ENDIAN)
# define CAT_ENDIAN_LITTLE
#endif

#if defined(__BIG_ENDIAN__) && defined(__LITTLE_ENDIAN__)
# if defined(__BYTE_ORDER__) && __BYTE_ORDER__== __BIG_ENDIAN__
#  define CAT_ENDIAN_BIG
# elif defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __LITTLE_ENDIAN__
#  define CAT_ENDIAN_LITTLE
# endif
#elif defined(__BIG_ENDIAN__)
# define CAT_ENDIAN_BIG
#elif defined(__LITTLE_ENDIAN__)
# define CAT_ENDIAN_LITTLE
#endif

// If no overrides defined, try to base it on ISA
// Okay -- Technically IA64 and PPC can switch endianness with an MSR bit
// But come on no one does that!  ...Right?
// If it's not right, make sure that one of the override flags in Config.hpp are defined.

#if !defined(CAT_ENDIAN_BIG) && !defined(CAT_ENDIAN_LITTLE)
# if defined(CAT_ISA_X86) || defined(CAT_ISA_VMS) || defined(CAT_ISA_ATMEGA)
#  define CAT_ENDIAN_LITTLE
# elif defined(CAT_ISA_PPC) || defined(CAT_ISA_AMIGA)
#  define CAT_ENDIAN_BIG
// Here are some other ISA that are not explicitly defined above:
# elif defined(applec) || defined(__AS400__) || defined(_CRAY) || defined(__hppa) || \
	   defined(__hp9000) || defined(ibm370) || defined(mc68000) || defined(m68k) || \
	   defined(__MRC__) || defined(__MVS__) || defined(sparc) || defined(__sparc) || \
	   defined(__VOS__) || defined(__TANDEM) || defined(__VMCMS__)
#  define CAT_ENDIAN_BIG
# else
#  define CAT_ENDIAN_UNKNOWN /* Must be detected at runtime */
# endif
#endif

#endif // #if !defined(CAT_ENDIAN_BIG) && !defined(CAT_ENDIAN_LITTLE)


//// Word Size ////

#if defined(_LP64) || defined(__LP64__) || defined(__arch64__) || \
	defined(_WIN64) || defined(_M_X64) || defined(__ia64) || \
	defined(__ia64__) || defined(__x86_64) || defined(_M_IA64) || \
	defined(__mips64)

# define CAT_WORD_64

// 64-bit MSVC does not support inline assembly
# if defined(CAT_COMPILER_MSVC)
#  undef CAT_ASM_INTEL
# endif

#else // Assuming 32-bit otherwise!

# define CAT_WORD_32

#endif

// __fastcall calling convention is rarely supported, and doesn't make sense for 64-bit targets
#if !defined(CAT_FASTCALL)
# define CAT_FASTCALL
#elif !defined(CAT_ISA_X86) || defined(CAT_WORD_64)
# undef CAT_FASTCALL
# define CAT_FASTCALL
#endif


//// Operating System ////

#if defined(__SVR4) && defined(__sun)
# define CAT_OS_SOLARIS

#elif defined(__APPLE__) && defined(TARGET_OS_IPHONE)
# define CAT_OS_IPHONE
# define CAT_OS_APPLE

#elif defined(__APPLE__) && (defined(__MACH__) || defined(__DARWIN__))
# define CAT_OS_OSX
# define CAT_OS_APPLE

#elif defined(__OpenBSD__) || defined(__NetBSD__) || defined(__FreeBSD__)
# define CAT_OS_BSD

#elif defined(__linux__) || defined(__unix__)
# define CAT_OS_LINUX

#elif defined(_WIN32_WCE)
# define CAT_OS_WINDOWS_CE
# define CAT_OS_WINDOWS /* Also defined */

#elif defined(_WIN32)
# define CAT_OS_WINDOWS

#elif defined(_XBOX) || defined(_X360)
# define CAT_OS_XBOX

#elif defined(_PS3) || defined(__PS3__) || defined(SN_TARGET_PS3)
# define CAT_OS_PS3

#elif defined(__OS2__)
# define CAT_OS_OS2

#elif defined(__APPLE__)
# define CAT_OS_APPLE

#elif defined(_aix) || defined(aix)
# define CAT_OS_AIX

#elif defined(HPUX)
# define CAT_OS_HPUX

#elif defined(IRIX)
# define CAT_OS_IRIX

#else
# define CAT_OS_UNKNOWN
#endif

// Detect CYGWIN environment
#if defined(__CYGWIN__) || defined(__CYGWIN32__)
# define CAT_CYGWIN
#endif


//// Dynamic Linkage ////

#if defined(CAT_OS_WINDOWS) || defined(CAT_CYGWIN)

# if defined(CAT_NEUTER_EXPORT)
#  define CAT_EXPORT /* Do not import or export any symbols */
# elif defined(CAT_BUILD_DLL)
#  define CAT_EXPORT CAT_DLL_EXPORT /* Implementing a DLL so export this symbol */
# else
#  define CAT_EXPORT CAT_DLL_IMPORT /* Using a DLL so import this symbol, faster on Windows */
# endif

#else

# undef CAT_DLL_EXPORT
# undef CAT_DLL_IMPORT
# define CAT_DLL_EXPORT
# define CAT_DLL_IMPORT
# define CAT_EXPORT

#endif

// Extern tag for dynamic linkage:

#if defined(CAT_NEUTER_EXPORT)
# define CAT_EXPORT_EXTERN
#elif defined(CAT_BUILD_DLL)
# define CAT_EXPORT_EXTERN
#else
# define CAT_EXPORT_EXTERN extern
#endif

// Export statically-linked templates: CAT_EXPORT_TEMPLATE(cat::Dictionary<int, int>);

#define CAT_EXPORT_TEMPLATE(T) CAT_EXPORT_EXTERN template class CAT_EXPORT T;


//// Basic types ////

#if defined(CAT_COMPILER_MSVC)

	// MSVC does not ship with stdint.h (C99 standard...)
	typedef unsigned __int8  u8;
	typedef signed __int8	 s8;
	typedef unsigned __int16 u16;
	typedef signed __int16   s16;
	typedef unsigned __int32 u32;
	typedef signed __int32   s32;
	typedef unsigned __int64 u64;
	typedef signed __int64   s64;

#else

} // namespace cat
#include <stdint.h>
namespace cat {

	// All other compilers use this
	typedef uint8_t  u8;
	typedef int8_t   s8;
	typedef uint16_t u16;
	typedef int16_t  s16;
	typedef uint32_t u32;
	typedef int32_t  s32;
	typedef uint64_t u64;
	typedef int64_t  s64;

#endif

#if defined(CAT_COMPILER_GCC) && defined(CAT_WORD_64)

	// GCC also adds 128-bit types :D
	typedef __uint128_t u128;
	typedef __int128_t  s128;

#endif

typedef float f32;
typedef double f64;

union Float32 {
	float f;
	u32 i;

	Float32(float n) { f = n; }
	Float32(u32 n) { i = n; }
};


//// String and buffer macros ////

} // namespace cat
#include <string.h>
namespace cat {

#if defined(CAT_OS_WINDOWS)
# define CAT_NEWLINE "\r\n"
#else
# define CAT_NEWLINE "\n"
#endif

#define CAT_FOREVER for (;;)

// Same as strncpy() in all ways except that the result is guaranteed to
// be a nul-terminated C string
#if defined(CAT_COMPILER_MSVC)
# define CAT_STRNCPY(dest, src, size) { strncpy_s(dest, size, src, size); (dest)[(size)-1] = '\0'; }
#else
# define CAT_STRNCPY(dest, src, size) { strncpy(dest, src, size); (dest)[(size)-1] = '\0'; }
#endif

// Because memory clearing is a frequent operation
#define CAT_CLR(dest, size) memset(dest, 0, size)

// Works for arrays, also
#define CAT_OBJCLR(obj) CAT_CLR((void*)&(obj), sizeof(obj))

// More secure memory clearing
#if defined(CAT_OS_WINDOWS) && !defined(CAT_OS_WINDOWS_CE) && !defined(CAT_COMPILER_MINGW)
# define CAT_SECURE_CLR(dest, size) SecureZeroMemory(dest, size)
#else

	// From https://www.securecoding.cert.org/confluence/display/cplusplus/MSC06-CPP.+Be+aware+of+compiler+optimization+when+dealing+with+sensitive+data
	CAT_INLINE void *cat_memset_s(void *v, int c, unsigned int n)
	{
		volatile unsigned char *p = (volatile unsigned char *)v;

		while (n--)
			*p++ = c;

		return v;
	}

# define CAT_SECURE_CLR(dest, size) cat_memset_s(dest, 0, size)
#endif

#define CAT_SECURE_OBJCLR(obj) CAT_SECURE_CLR((void*)&(obj), sizeof(obj))

// Stringize
#define CAT_STRINGIZE(X) DO_CAT_STRINGIZE(X)
#define DO_CAT_STRINGIZE(X) #X

// Fine and line string macro
#define CAT_FILE_LINE_STRING __FILE__ ":" CAT_STRINGIZE(__LINE__)

// Variable-length data trailing a struct
template<typename T> CAT_INLINE u8 *GetTrailingBytes(T *t) { return reinterpret_cast<u8*>( t ) + sizeof(T); }

// Bounds
template<typename T> CAT_INLINE T BoundMin(const T &minimum, const T &x)
{
	if (x < minimum) return minimum;
	return x;
}

template<typename T> CAT_INLINE T BoundMax(const T &maximum, const T &x)
{
	if (x > maximum) return maximum;
	return x;
}

template<typename T> CAT_INLINE T Bound(const T &minimum, const T &maximum, const T &x)
{
	if (x < minimum) return minimum;
	if (x > maximum) return maximum;
	return x;
}


//// Miscellaneous bitwise macros ////

#define CAT_BITCLRHI8(reg, count) ((u8)((u8)(reg) << (count)) >> (count)) /* sets to zero a number of high bits in a byte */
#define CAT_BITCLRLO8(reg, count) ((u8)((u8)(reg) >> (count)) << (count)) /* sets to zero a number of low bits in a byte */
#define CAT_BITCLRHI16(reg, count) ((u16)((u16)(reg) << (count)) >> (count)) /* sets to zero a number of high bits in a 16-bit word */
#define CAT_BITCLRLO16(reg, count) ((u16)((u16)(reg) >> (count)) << (count)) /* sets to zero a number of low bits in a 16-bit word */
#define CAT_BITCLRHI32(reg, count) ((u32)((u32)(reg) << (count)) >> (count)) /* sets to zero a number of high bits in a 32-bit word */
#define CAT_BITCLRLO32(reg, count) ((u32)((u32)(reg) >> (count)) << (count)) /* sets to zero a number of low bits in a 32-bit word */


//// Integer macros ////

#define CAT_AT_LEAST_2_BITS(n) ( (n) & ((n) - 1) )
#define CAT_LSB32(n) ( (n) & (u32)(-(s32)(n)) ) /* 0 -> 0 */
#define CAT_LSB64(n) ( (n) & (u64)(-(s64)(n)) ) /* 0 -> 0 */
#define CAT_IS_POWER_OF_2(n) ( n && !CAT_AT_LEAST_2_BITS(n) )

// Safely take the average of two numbers without possibility of overflow
#define CAT_SAFE_AVERAGE(A, B) (((A) & (B)) + (((A) ^ (B)) >> 1))

// Bump 'n' to the next unit of 'width'
// 0=CAT_CEIL_UNIT(0, 16), 1=CAT_CEIL_UNIT(1, 16), 1=CAT_CEIL_UNIT(16, 16), 2=CAT_CEIL_UNIT(17, 16)
#define CAT_CEIL_UNIT(n, width) ( ( (n) + (width) - 1 ) / (width) )
// 0=CAT_CEIL(0, 16), 16=CAT_CEIL(1, 16), 16=CAT_CEIL(16, 16), 32=CAT_CEIL(17, 16)
#define CAT_CEIL(n, width) ( CAT_CEIL_UNIT(n, width) * (width) )


//// Rotation macros ////

#define CAT_ROL8(n, r)  ( ((u8)(n) << (r)) | ((u8)(n) >> ( 8 - (r))) ) /* only works for u8 */
#define CAT_ROR8(n, r)  ( ((u8)(n) >> (r)) | ((u8)(n) << ( 8 - (r))) ) /* only works for u8 */
#define CAT_ROL16(n, r) ( ((u16)(n) << (r)) | ((u16)(n) >> (16 - (r))) ) /* only works for u16 */
#define CAT_ROR16(n, r) ( ((u16)(n) >> (r)) | ((u16)(n) << (16 - (r))) ) /* only works for u16 */
#define CAT_ROL32(n, r) ( ((u32)(n) << (r)) | ((u32)(n) >> (32 - (r))) ) /* only works for u32 */
#define CAT_ROR32(n, r) ( ((u32)(n) >> (r)) | ((u32)(n) << (32 - (r))) ) /* only works for u32 */
#define CAT_ROL64(n, r) ( ((u64)(n) << (r)) | ((u64)(n) >> (64 - (r))) ) /* only works for u64 */
#define CAT_ROR64(n, r) ( ((u64)(n) >> (r)) | ((u64)(n) << (64 - (r))) ) /* only works for u64 */


//// Byte-order swapping ////

#define CAT_BOSWAP16(n) CAT_ROL16(n, 8)
#define CAT_BOSWAP32(n) ( (CAT_ROL32(n, 8) & 0x00ff00ff) | (CAT_ROL32(n, 24) & 0xff00ff00) )
#define CAT_BOSWAP64(n) ( ((u64)CAT_BOSWAP32((u32)n) << 32) | CAT_BOSWAP32((u32)(n >> 32)) )


//// Intrinsics ////

#if defined(CAT_OS_WINDOWS_CE)

#pragma intrinsic(_lrotl)
#pragma intrinsic(_lrotr)

#undef CAT_ROL32
#undef CAT_ROR32

#define CAT_ROL32(n, r) _lrotl(n, r)
#define CAT_ROR32(n, r) _lrotr(n, r)

#elif defined(CAT_COMPILER_MSVC)

#pragma intrinsic(_rotl, _rotr, _rotl64, _rotr64)
#pragma intrinsic(_byteswap_ushort, _byteswap_ulong, _byteswap_uint64)
#pragma intrinsic(_BitScanForward, _BitScanReverse, _bittestandset)
#pragma intrinsic(__emulu)
#pragma intrinsic(_InterlockedExchange)
#pragma intrinsic(_interlockedbittestandset)
#pragma intrinsic(_interlockedbittestandreset)

#if defined(CAT_WORD_64)
#pragma intrinsic(__rdtsc)
#pragma intrinsic(_umul128)
#pragma intrinsic(_BitScanForward64, _BitScanReverse64, _bittestandset64)
#pragma intrinsic(_InterlockedCompareExchange128)
#else
#pragma intrinsic(_InterlockedCompareExchange64)
#endif

#undef CAT_ROL32
#undef CAT_ROR32
#undef CAT_ROL64
#undef CAT_ROR64
#undef CAT_BOSWAP16
#undef CAT_BOSWAP32
#undef CAT_BOSWAP64

#define CAT_ROL32(n, r) _rotl(n, r)
#define CAT_ROR32(n, r) _rotr(n, r)
#define CAT_ROL64(n, r) _rotl64(n, r)
#define CAT_ROR64(n, r) _rotr64(n, r)
#define CAT_BOSWAP16(n) _byteswap_ushort(n)
#define CAT_BOSWAP32(n) _byteswap_ulong(n)
#define CAT_BOSWAP64(n) _byteswap_uint64(n)

#endif


//// No Copy ////

#define CAT_NO_COPY(T) \
private: \
	CAT_INLINE T(const T&) {} \
	CAT_INLINE T& operator=(const T&) { return *this; }


} // namespace cat


//// Memory Leaks ////

#define CAT_DEBUG_LEAKS_DUMP()
#if defined(CAT_DEBUG_LEAKS)
# if defined(CAT_COMPILER_MSVC)
#  define _CRTDBG_MAP_ALLOC
#  include <stdlib.h>
#  include <crtdbg.h>
#  undef CAT_DEBUG_LEAKS_DUMP
#  define CAT_DEBUG_LEAKS_DUMP() _CrtDumpMemoryLeaks()
# endif
#endif

#endif // CAT_PLATFORM_HPP
