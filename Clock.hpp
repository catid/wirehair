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

#ifndef CAT_CLOCK_HPP
#define CAT_CLOCK_HPP

#include "codec_source/Platform.hpp"
#include <string>

namespace cat {


class CAT_EXPORT Clock
{
#ifdef CAT_OS_WINDOWS
	// Windows version requires some initialization
	static const int LOWEST_ACCEPTABLE_PERIOD = 10;

    u32 _period;		// timegettime() and Windows scheduler period
	double _inv_freq;	// Performance counter frequency (does not change, so cache it)
#endif

public:
	bool OnInitialize();
	void OnFinalize();

	static u32 sec();			// Timestamp in seconds
    u32 msec_fast();			// Timestamp in milliseconds, less accurate than msec() but faster
    u32 msec();					// Timestamp in milliseconds
	double usec();				// Timestamp in microseconds
	static u32 cycles();		// Timestamp in cycles

    static std::string format(const char *format_string);

    static void sleep(u32 milliseconds);

    static u32 MeasureClocks(int iterations, void (*FunctionPtr)());
};


} // namespace cat

#endif // CAT_CLOCK_HPP
