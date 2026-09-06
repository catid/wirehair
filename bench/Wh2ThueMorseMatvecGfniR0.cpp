#include "Wh2ThueMorseMatvecGfniR0.h"

#ifdef WH_COUNT
#error "This benchmark-only matvec does not implement GF work counters"
#endif

// This translation unit supplies the ONE shared GF runtime in the experiment.
// Do not link a second gf256.cpp object alongside it. Including the unchanged
// source keeps its private dispatch flags and multiplication matrices private.
#include "../gf256.cpp"

namespace wh2_matvec_gfni_r0 {

void ApplyScalar(const std::uint8_t* matrix, std::uint8_t* vector)
{
    std::uint8_t result[6];
    for (unsigned r = 0; r < 6; ++r) {
        std::uint8_t value = gf256_mul(vector[0], matrix[r * 6]);
        for (unsigned c = 1; c < 6; ++c)
            value ^= gf256_mul(vector[c], matrix[r * 6 + c]);
        result[r] = value;
    }
    std::memcpy(vector, result, sizeof(result));
}

namespace {

#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_TARGET_GFNI)

static GF256_GFNI_TARGET void ApplyTarget(
    const std::uint8_t* matrix, std::uint8_t* vector)
{
    // Lane c contains the six entries of column c, in byte positions 0..5.
    // Zero padding makes both the unused row bytes and lanes 6/7 inert.
    // Make vector-store alignment explicit, including ASAN's fake-stack layout.
    // GCC may vectorize these initializers even though the loads below are unaligned.
    alignas(64) std::uint64_t columns[8] = {};
    alignas(64) std::uint64_t multipliers[8] = {};
    for (unsigned c = 0; c < 6; ++c) {
        for (unsigned r = 0; r < 6; ++r)
            columns[c] |= static_cast<std::uint64_t>(matrix[r * 6 + c]) << (8 * r);
        multipliers[c] = TargetGFNIMulMatrix[vector[c]];
    }

    // T(vector[c]) is the original field-correct affine multiplication map,
    // not the hardware 0x11b polynomial multiply. All 36 products are formed
    // by this one affine instruction without changing the byte representation.
    const __m512i products = _mm512_gf2p8affine_epi64_epi8(
        _mm512_loadu_si512(static_cast<const void*>(columns)),
        _mm512_loadu_si512(static_cast<const void*>(multipliers)), 0);
    const __m256i fold256 = _mm256_xor_si256(
        _mm512_castsi512_si256(products), _mm512_extracti64x4_epi64(products, 1));
    const __m128i fold128 = _mm_xor_si128(
        _mm256_castsi256_si128(fold256), _mm256_extracti128_si256(fold256, 1));
    const __m128i fold64 = _mm_xor_si128(fold128, _mm_srli_si128(fold128, 8));
    std::uint64_t result = 0;
    _mm_storel_epi64(reinterpret_cast<__m128i*>(&result), fold64);

    // matrix and vector may overlap. No caller storage is written before all
    // matrix and vector bytes above have been consumed; never store eight bytes.
    std::memcpy(vector, &result, 6);
}

#endif

} // namespace

bool Available()
{
#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_TARGET_GFNI)
    return CpuHasTargetGFNI && CpuHasAVX2;
#else
    return false;
#endif
}

void Apply(const std::uint8_t* matrix, std::uint8_t* vector)
{
#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_TARGET_GFNI)
    if (CpuHasTargetGFNI && CpuHasAVX2) {
        ApplyTarget(matrix, vector);
        return;
    }
#endif
    ApplyScalar(matrix, vector);
}

} // namespace wh2_matvec_gfni_r0
