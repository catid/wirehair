// Compile the unchanged benchmark codec with precisely one payload-call seam.
// Include the real GF declarations before substitution; the shared runtime is
// neither renamed nor recompiled with a candidate-only ISA/behavior setting.
#include "Wh2ThueMorseTinyPayloadR0.h"

#define wh2_thue_native wh2_thue_tiny_payload_r0
#define gf256_mulset_multi_mem wh2_thue_tiny_payload_r0::Payload
#include "Wh2ThueMorseNativeCodec.cpp"
#undef gf256_mulset_multi_mem
#undef wh2_thue_native
