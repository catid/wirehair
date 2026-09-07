// The only equation-neutral codec seam is the existing six-source payload call.
#include "Wh2ThueMorseTinyGfniR0.h"
#define wh2_thue_native wh2_thue_tiny_gfni_r0
#define gf256_mulset_multi_mem wh2_thue_tiny_gfni_r0::Payload
#include "Wh2ThueMorseNativeCodec.cpp"
#undef gf256_mulset_multi_mem
#undef wh2_thue_native
