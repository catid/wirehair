// Reuse the exact four-arm bridge contract suite with only its ISA query redirected.
// Process the old declaration guard first; no old GFNI implementation is linked.
#include "Wh2ThueMorseMatvecGfniR0.h"
#include "Wh2ThueMorseTinyGfniR0.h"
#define wh2_matvec_gfni_r0 wh2_thue_tiny_gfni_r0
#include "Wh2ThueMorseMatvecCostR0BridgeTest.cpp"
#undef wh2_matvec_gfni_r0
