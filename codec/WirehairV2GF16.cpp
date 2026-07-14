#include "WirehairV2GF16.h"
#include "../gf256.h"

#include <limits>
#include <cstring>
#include <mutex>
#include <stddef.h>
#include <stdint.h>

namespace wirehair_v2 {
namespace {

std::once_flag InitOnce;
bool InitResult = false;
uint16_t Coefficients[kMixedGF16RowsMax][kMixedCoefficientPeriod] = {};
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
uint16_t SharedXCoefficients[
    kMixedGF16RowsMax][kMixedCoefficientPeriod] = {};
#endif

uint16_t MultiplyUnchecked(uint16_t x, uint16_t y)
{
    const uint8_t a = (uint8_t)x, b = (uint8_t)(x >> 8);
    const uint8_t c = (uint8_t)y, d = (uint8_t)(y >> 8);
    const uint8_t ac = gf256_mul(a, c);
    const uint8_t bd = gf256_mul(b, d);
    const uint8_t cross = gf256_mul((uint8_t)(a ^ b), (uint8_t)(c ^ d));
    const uint8_t low = (uint8_t)(ac ^ gf256_mul(kGF16Lambda, bd));
    const uint8_t high = (uint8_t)(cross ^ ac); // ad+bc+bd
    return (uint16_t)low | (uint16_t)high << 8;
}

uint16_t PowerUnchecked(uint16_t x, uint32_t exponent)
{
    uint16_t result = 1u;
    while (exponent != 0u) {
        if (exponent & 1u) result = MultiplyUnchecked(result, x);
        exponent >>= 1;
        if (exponent != 0u) x = MultiplyUnchecked(x, x);
    }
    return result;
}

uint16_t InverseUnchecked(uint16_t x)
{
    return x == 0u ? 0u : PowerUnchecked(x, 65534u);
}

bool InitializeUnchecked()
{
    if (gf256_init() != 0) return false;
    for (uint32_t x = 0; x < 256u; ++x) {
        if ((uint8_t)(gf256_mul((uint8_t)x, (uint8_t)x) ^ x) ==
            kGF16Lambda) return false;
    }
    const uint32_t factors[] = {3u, 5u, 17u, 257u};
    for (uint32_t factor : factors) {
        if (PowerUnchecked(kGF16Generator, 65535u / factor) == 1u)
            return false;
    }
    for (uint32_t row = 0; row < kMixedGF16RowsMax; ++row) {
        const uint16_t y = PowerUnchecked(kGF16Generator, 1000u + row);
        if ((y >> 8) == 0u) return false;
        for (uint32_t column = 0; column < kMixedCoefficientPeriod; ++column) {
            const uint16_t denominator =
                (uint16_t)(PowerUnchecked(kGF16Generator, column) ^ y);
            if (denominator == 0u) return false;
            Coefficients[row][column] = InverseUnchecked(denominator);

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            // The first ten mixed rows use X=12+column and Y=row in the
            // GF(256) subfield.  Reusing that X here with Y coordinates
            // outside the subfield makes every active row one Cauchy matrix.
            const uint16_t shared_x = (uint16_t)(
                kMixedGF256Rows + kMixedGF16Rows + column);
            SharedXCoefficients[row][column] =
                InverseUnchecked((uint16_t)(shared_x ^ y));
#endif
        }
    }
    return Coefficients[0][0] == 34916u &&
        Coefficients[1][0] == 2472u &&
        Coefficients[0][243] == 59155u;
}

bool ValidBlock(const void* destination, const void* source, uint32_t bytes)
{
    return destination && source && bytes != 0u && (bytes & 1u) == 0u;
}

bool RangesOverlap(const void* a, size_t a_bytes, const void* b, size_t b_bytes)
{
    const uintptr_t a0 = (uintptr_t)a;
    const uintptr_t b0 = (uintptr_t)b;
    if (a_bytes > UINTPTR_MAX - a0 || b_bytes > UINTPTR_MAX - b0)
        return true;
    return a0 < b0 + b_bytes && b0 < a0 + a_bytes;
}

} // namespace

bool InitializeGF16()
{
    std::call_once(InitOnce, []() { InitResult = InitializeUnchecked(); });
    return InitResult;
}

uint16_t GF16Multiply(uint16_t x, uint16_t y)
{
    return InitializeGF16() ? MultiplyUnchecked(x, y) : 0u;
}

uint16_t GF16Inverse(uint16_t x)
{
    return InitializeGF16() ? InverseUnchecked(x) : 0u;
}

uint16_t GF16MultiplyInitialized(uint16_t x, uint16_t y)
{
    return MultiplyUnchecked(x, y);
}

uint16_t GF16InverseInitialized(uint16_t x)
{
    return InverseUnchecked(x);
}

uint16_t MixedGF16Coefficient(uint32_t extension_row, uint32_t column)
{
    if (!InitializeGF16() || extension_row >= kMixedGF16RowsMax) return 0u;
    return Coefficients[extension_row][column % kMixedCoefficientPeriod];
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
uint16_t MixedGF16SharedXCoefficient(
    uint32_t extension_row,
    uint32_t column)
{
    if (!InitializeGF16() || extension_row >= kMixedGF16RowsMax) return 0u;
    return SharedXCoefficients[
        extension_row][column % kMixedCoefficientPeriod];
}
#endif

bool GF16MulAddMem(
    void* destination, uint16_t scale, const void* source, uint32_t bytes)
{
    if (!InitializeGF16() || !ValidBlock(destination, source, bytes) ||
        (destination != source &&
         RangesOverlap(destination, bytes, source, bytes)))
        return false;
    uint8_t* z = (uint8_t*)destination;
    const uint8_t* x = (const uint8_t*)source;
    for (uint32_t offset = 0; offset < bytes; offset += 2u) {
        const uint16_t value =
            (uint16_t)x[offset] | (uint16_t)x[offset + 1u] << 8;
        const uint16_t product = MultiplyUnchecked(scale, value);
        z[offset] ^= (uint8_t)product;
        z[offset + 1u] ^= (uint8_t)(product >> 8);
    }
    return true;
}

bool GF16ScaleMem(void* block, uint16_t scale, uint32_t bytes)
{
    if (!InitializeGF16() || !ValidBlock(block, block, bytes)) return false;
    uint8_t* x = (uint8_t*)block;
    for (uint32_t offset = 0; offset < bytes; offset += 2u) {
        const uint16_t value =
            (uint16_t)x[offset] | (uint16_t)x[offset + 1u] << 8;
        const uint16_t product = MultiplyUnchecked(scale, value);
        x[offset] = (uint8_t)product;
        x[offset + 1u] = (uint8_t)(product >> 8);
    }
    return true;
}

bool GF16Deinterleave(
    const void* source, void* low, void* high, uint32_t bytes)
{
    const size_t plane_bytes = bytes / 2u;
    if (!InitializeGF16() || !source || !low || !high ||
        bytes == 0u || (bytes & 1u) != 0u ||
        RangesOverlap(source, bytes, low, plane_bytes) ||
        RangesOverlap(source, bytes, high, plane_bytes) ||
        RangesOverlap(low, plane_bytes, high, plane_bytes)) return false;
    const uint8_t* x = (const uint8_t*)source;
    uint8_t* lo = (uint8_t*)low;
    uint8_t* hi = (uint8_t*)high;
    for (uint32_t i = 0; i < bytes / 2u; ++i) {
        lo[i] = x[2u * i];
        hi[i] = x[2u * i + 1u];
    }
    return true;
}

bool GF16Interleave(
    const void* low, const void* high, void* destination, uint32_t bytes)
{
    const size_t plane_bytes = bytes / 2u;
    if (!InitializeGF16() || !low || !high || !destination ||
        bytes == 0u || (bytes & 1u) != 0u ||
        RangesOverlap(low, plane_bytes, high, plane_bytes) ||
        RangesOverlap(low, plane_bytes, destination, bytes) ||
        RangesOverlap(high, plane_bytes, destination, bytes)) return false;
    const uint8_t* lo = (const uint8_t*)low;
    const uint8_t* hi = (const uint8_t*)high;
    uint8_t* x = (uint8_t*)destination;
    for (uint32_t i = 0; i < bytes / 2u; ++i) {
        x[2u * i] = lo[i];
        x[2u * i + 1u] = hi[i];
    }
    return true;
}

bool GF16MulAddPlanar(
    void* destination_low,
    void* destination_high,
    uint16_t scale,
    const void* source_low,
    const void* source_high,
    uint32_t elements)
{
    if (!InitializeGF16() || !destination_low || !destination_high ||
        !source_low || !source_high || elements == 0u ||
        elements > (uint32_t)std::numeric_limits<int>::max() ||
        RangesOverlap(destination_low, elements, destination_high, elements) ||
        RangesOverlap(destination_low, elements, source_low, elements) ||
        RangesOverlap(destination_low, elements, source_high, elements) ||
        RangesOverlap(destination_high, elements, source_low, elements) ||
        RangesOverlap(destination_high, elements, source_high, elements) ||
        RangesOverlap(source_low, elements, source_high, elements)) return false;
    if (scale == 0u) return true;
    if (scale == 1u) {
        gf256_add_mem(destination_low, source_low, (int)elements);
        gf256_add_mem(destination_high, source_high, (int)elements);
        return true;
    }
    const uint8_t a = (uint8_t)scale, b = (uint8_t)(scale >> 8);
    void* destinations[2] = {destination_low, destination_high};
    const uint8_t low_scales[2] = {a, b};
    const uint8_t high_scales[2] = {
        gf256_mul(kGF16Lambda, b), (uint8_t)(a ^ b)
    };
    gf256_muladd_multi_mem(
        destinations, low_scales, 2, source_low, (int)elements);
    gf256_muladd_multi_mem(
        destinations, high_scales, 2, source_high, (int)elements);
    return true;
}

bool GF16MulAddPlanar2(
    void* destination0_low,
    void* destination0_high,
    uint16_t scale0,
    void* destination1_low,
    void* destination1_high,
    uint16_t scale1,
    const void* source_low,
    const void* source_high,
    uint32_t elements)
{
    void* destinations[4] = {
        destination0_low, destination0_high,
        destination1_low, destination1_high
    };
    if (!InitializeGF16() || !source_low || !source_high || elements == 0u ||
        elements > (uint32_t)std::numeric_limits<int>::max() ||
        RangesOverlap(source_low, elements, source_high, elements))
    {
        return false;
    }
    for (uint32_t i = 0; i < 4u; ++i)
    {
        if (!destinations[i] ||
            RangesOverlap(destinations[i], elements, source_low, elements) ||
            RangesOverlap(destinations[i], elements, source_high, elements))
        {
            return false;
        }
        for (uint32_t j = 0; j < i; ++j) {
            if (RangesOverlap(
                    destinations[i], elements, destinations[j], elements))
            {
                return false;
            }
        }
    }
    if (scale0 == 0u && scale1 == 0u) return true;

    const uint8_t a0 = (uint8_t)scale0;
    const uint8_t b0 = (uint8_t)(scale0 >> 8);
    const uint8_t a1 = (uint8_t)scale1;
    const uint8_t b1 = (uint8_t)(scale1 >> 8);
    const uint8_t low_scales[4] = {a0, b0, a1, b1};
    const uint8_t high_scales[4] = {
        gf256_mul(kGF16Lambda, b0), (uint8_t)(a0 ^ b0),
        gf256_mul(kGF16Lambda, b1), (uint8_t)(a1 ^ b1)
    };
    gf256_muladd_multi_mem(
        destinations, low_scales, 4, source_low, (int)elements);
    gf256_muladd_multi_mem(
        destinations, high_scales, 4, source_high, (int)elements);
    return true;
}

bool GF16ScalePlanar(
    void* low,
    void* high,
    uint16_t scale,
    void* scratch,
    uint32_t elements)
{
    if (!InitializeGF16() || !low || !high || !scratch || elements == 0u ||
        elements > (uint32_t)std::numeric_limits<int>::max() ||
        RangesOverlap(low, elements, high, elements) ||
        RangesOverlap(low, elements, scratch, elements) ||
        RangesOverlap(high, elements, scratch, elements)) return false;
    if (scale == 1u) return true;
    if (scale == 0u) {
        std::memset(low, 0, elements);
        std::memset(high, 0, elements);
        return true;
    }
    const int bytes = (int)elements;
    const uint8_t a = (uint8_t)scale;
    const uint8_t b = (uint8_t)(scale >> 8);
    std::memcpy(scratch, low, elements);
    gf256_mul_mem(low, scratch, a, bytes);
    gf256_muladd_mem(low, gf256_mul(kGF16Lambda, b), high, bytes);
    gf256_mul_mem(high, high, (uint8_t)(a ^ b), bytes);
    gf256_muladd_mem(high, b, scratch, bytes);
    return true;
}

} // namespace wirehair_v2
