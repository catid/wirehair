#ifndef WH2_THUE_MORSE_MULROWS_R0_H
#define WH2_THUE_MORSE_MULROWS_R0_H

#include "Wh2ThueMorseNativeCodec.h"

// Benchmark-only .65 candidate. Initialize the shared GF256 runtime first.
// The lookup, validation, equations and output contract are unchanged.
namespace wh2_thue_mulrows_r0 {

wh2_thue_native::Status Row(
    wh2_thue_native::Lookup lookup, std::uint32_t id,
    std::uint8_t output[wh2_thue_native::kSourceCount]);

} // namespace wh2_thue_mulrows_r0
#endif
