/*
	Copyright (c) 2009-2012 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of LibCat nor the names of its contributors may be used
	  to endorse or promote products derived from this software without
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

#ifndef CAT_CONFIG_HPP
#define CAT_CONFIG_HPP

namespace cat {


// This definition overrides CAT_BUILD_DLL below.  Neuters CAT_EXPORT macro so symbols are
// neither exported or imported.
#define CAT_NEUTER_EXPORT

// Define this to greatly weaken the Fortuna implementation, but consume much less CPU.
// Only disable this if you are serious about security of the implementation at expense of performance.
#define CAT_NO_ENTROPY_THREAD

// This definition changes the meaning of the CAT_EXPORT macro on Windows.  When defined,
// the CAT_EXPORT macro will export the associated symbol.  When undefined, it will import it.
//#define CAT_BUILD_DLL

// If you want to remove server-side code from a binary distribution of a client program:
//#define CAT_OMIT_SERVER_CODE

// If you know the endianness of your target, uncomment one of these for better performance.
//#define CAT_ENDIAN_BIG
//#define CAT_ENDIAN_LITTLE

// If you want to use faster 384-bit or 512-bit math, define this:
//#define CAT_UNROLL_OVER_256_BITS

// Adjust if your architecture uses larger than 128-byte cache line
#define CAT_DEFAULT_CACHE_LINE_SIZE 128
#define CAT_DEFAULT_CPU_COUNT 1
#define CAT_DEFAULT_PAGE_SIZE 65536
#define CAT_DEFAULT_ALLOCATION_GRANULARITY CAT_DEFAULT_PAGE_SIZE
#define CAT_DEFAULT_SECTOR_SIZE 512

// Enable leak debug mode of the common runtime heap allocator
#define CAT_DEBUG_LEAKS

// Enable RefObject debug trace mode
//#define CAT_TRACE_REFOBJECT

// Enable event re-ordering for better batching in WorkerThreads
#define CAT_WORKER_THREADS_REORDER_EVENTS

// Dump extra settings information to the console for debugging
#define CAT_SETTINGS_VERBOSE

// Specify which file to use for persisting settings between sessions
#if !defined(CAT_SETTINGS_FILE)
#define CAT_SETTINGS_FILE "Settings.cfg"
#endif
#if !defined(CAT_SETTINGS_OVERRIDE_FILE)
#define CAT_SETTINGS_OVERRIDE_FILE "Override.cfg"
#endif

// Enable Ragdoll-based files to store empty keys
#define CAT_RAGDOLL_STORE_EMPTY

// Enable multi-threaded version of the logger, which reduces latency
#define CAT_THREADED_LOGGER

// Enable Re-use UDP Send Allocator
#define CAT_UDP_SEND_ALLOCATOR

// When not in debug mode, enable/disable levels of logging
#define CAT_RELEASE_DISABLE_INANE
//#define CAT_RELEASE_DISABLE_INFO
//#define CAT_RELEASE_DISABLE_WARN
//#define CAT_RELEASE_DISABLE_OOPS
//#define CAT_RELEASE_DISABLE_FATAL

// Define this to enable auditing of data security code
//#define CAT_AUDIT


} // namespace cat

#endif // CAT_CONFIG_HPP
