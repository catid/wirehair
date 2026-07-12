#pragma once
#include <stdint.h>

namespace wirehair_v2 {

// Quadratic extension GF(256)[u]/(u^2+u+32).  A uint16 stores a+b*u as
// low byte a, high byte b.  Payload symbols are adjacent byte pairs.
static const uint8_t kGF16Lambda = 32u;
static const uint16_t kGF16Generator = 266u;
static const uint32_t kMixedCoefficientPeriod = 244u;
static const uint32_t kMixedGF256Rows = 10u;
static const uint32_t kMixedGF16Rows = 2u;

enum class CompletionField : uint32_t
{
    GF256 = 0,
    MixedGF256GF16 = 1
};

bool InitializeGF16();
uint16_t GF16Multiply(uint16_t x, uint16_t y);
uint16_t GF16Inverse(uint16_t x);
// Scalar fast paths for hot elimination loops after one successful
// InitializeGF16() call.  Calling them before initialization is invalid.
uint16_t GF16MultiplyInitialized(uint16_t x, uint16_t y);
uint16_t GF16InverseInitialized(uint16_t x);

// extension_row is 0 or 1. X=g^(column mod 244), Y=g^(1000+row).
uint16_t MixedGF16Coefficient(uint32_t extension_row, uint32_t column);

// Block operations reject odd/zero sizes before writing.  MulAdd supports
// exact in-place destination==source, but rejects partial overlap.
bool GF16MulAddMem(
    void* destination, uint16_t scale, const void* source, uint32_t bytes);
bool GF16ScaleMem(void* block, uint16_t scale, uint32_t bytes);
bool GF16Deinterleave(
    const void* source, void* low, void* high, uint32_t bytes);
bool GF16Interleave(
    const void* low, const void* high, void* destination, uint32_t bytes);

// Each plane contains elements bytes.  This uses four optimized GF(256)
// bulk operations for one extension-field constant muladd.  All plane
// buffers and all conversion inputs/outputs must be pairwise non-overlapping;
// overlap is rejected before writing.
bool GF16MulAddPlanar(
    void* destination_low,
    void* destination_high,
    uint16_t scale,
    const void* source_low,
    const void* source_high,
    uint32_t elements);
// In-place planar scale. scratch contains elements bytes and must not overlap
// either plane. Rejected arguments leave both planes and scratch untouched.
bool GF16ScalePlanar(
    void* low,
    void* high,
    uint16_t scale,
    void* scratch,
    uint32_t elements);

} // namespace wirehair_v2
