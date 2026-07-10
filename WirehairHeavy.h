/** \file
    \brief Portable helpers for heavy-row elimination
    \copyright Copyright (c) 2012-2018 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Wirehair nor the names of its contributors may be
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

#ifndef WIREHAIR_HEAVY_H
#define WIREHAIR_HEAVY_H

#include "gf256.h"

#include <cstring>

namespace wirehair {

/**
    Return a native uint32_t whose four in-memory bytes are the low four bits
    of `bits`, least-significant bit first.

    Loading a byte-defined table with memcpy keeps the memory layout identical
    on little- and big-endian hosts.  Multiplication by an 8-bit coefficient
    then broadcasts that coefficient to each selected byte without carries.
*/
static GF256_FORCE_INLINE uint32_t HeavyWindowWord(unsigned bits)
{
    static const uint8_t kByteLookup[16][4] = {
        { 0, 0, 0, 0 }, { 1, 0, 0, 0 },
        { 0, 1, 0, 0 }, { 1, 1, 0, 0 },
        { 0, 0, 1, 0 }, { 1, 0, 1, 0 },
        { 0, 1, 1, 0 }, { 1, 1, 1, 0 },
        { 0, 0, 0, 1 }, { 1, 0, 0, 1 },
        { 0, 1, 0, 1 }, { 1, 1, 0, 1 },
        { 0, 0, 1, 1 }, { 1, 0, 1, 1 },
        { 0, 1, 1, 1 }, { 1, 1, 1, 1 }
    };

    uint32_t word;
    std::memcpy(&word, kByteLookup[bits & 15u], sizeof(word));
    return word;
}

/** XOR one four-column binary window, scaled by coefficient, into bytes. */
static GF256_FORCE_INLINE void AddHeavyWindow(
    uint8_t* GF256_RESTRICT destination,
    unsigned bits,
    uint8_t coefficient)
{
    uint32_t word;
    std::memcpy(&word, destination, sizeof(word));
    word ^= HeavyWindowWord(bits) * coefficient;
    std::memcpy(destination, &word, sizeof(word));
}

} // namespace wirehair

#endif // WIREHAIR_HEAVY_H
