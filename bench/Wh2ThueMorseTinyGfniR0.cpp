#include "Wh2ThueMorseTinyGfniR0.h"
#include "Wh2ThueMorseTinyPayloadR0.h"

// One shared runtime only. The experiment adds no initialization, field table,
// production edit or public API. Do not also link a separate gf256.cpp.
#include "../gf256.cpp"

namespace wh2_thue_tiny_gfni_r0 {

void Scalar(void* output, const void* const* sources, const std::uint8_t* scales)
{
    wh2_thue_tiny_payload_r0::Payload(output, sources, scales, 6, 2);
}

namespace {
#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_TARGET_GFNI)

inline std::uint16_t ReadWord(const void* source)
{
    std::uint16_t word;
    std::memcpy(&word, source, sizeof(word));
    return word;
}

static GF256_GFNI_TARGET void ApplyTarget(
    void* output, const void* const* sources, const std::uint8_t* scales)
{
    // Each source occupies one word. No wide read crosses its two-byte span.
    __m128i data = _mm_cvtsi32_si128(ReadWord(sources[0]));
    data = _mm_insert_epi16(data, ReadWord(sources[1]), 1);
    data = _mm_insert_epi16(data, ReadWord(sources[2]), 2);
    data = _mm_insert_epi16(data, ReadWord(sources[3]), 3);
    data = _mm_insert_epi16(data, ReadWord(sources[4]), 4);
    data = _mm_insert_epi16(data, ReadWord(sources[5]), 5);

    std::int32_t first_scales;
    std::memcpy(&first_scales, scales, sizeof(first_scales));
    __m128i coefficients = _mm_cvtsi32_si128(first_scales);
    coefficients = _mm_insert_epi16(coefficients, ReadWord(scales + 4), 2);
    coefficients = _mm_unpacklo_epi8(coefficients, coefficients);

    const __m128i forward = _mm_set1_epi64x(kToHardware);
    data = _mm_gf2p8affine_epi64_epi8(data, forward, 0);
    coefficients = _mm_gf2p8affine_epi64_epi8(coefficients, forward, 0);
    __m128i product = _mm_gf2p8mul_epi8(data, coefficients);

    // Add all eight word lanes (last two are zero). The inverse basis map is
    // linear, so convert only after reducing the products by XOR.
    product = _mm_xor_si128(product, _mm_srli_si128(product, 8));
    product = _mm_xor_si128(product, _mm_srli_si128(product, 4));
    product = _mm_xor_si128(product, _mm_srli_si128(product, 2));
    product = _mm_gf2p8affine_epi64_epi8(
        product, _mm_set1_epi64x(kFromHardware), 0);
    const std::uint16_t result = static_cast<std::uint16_t>(_mm_cvtsi128_si32(product));
    std::memcpy(output, &result, sizeof(result));
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

void Apply(void* output, const void* const* sources, const std::uint8_t* scales)
{
#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_TARGET_GFNI)
    if (CpuHasTargetGFNI && CpuHasAVX2) {
        ApplyTarget(output, sources, scales);
        return;
    }
#endif
    Scalar(output, sources, scales);
}

} // namespace wh2_thue_tiny_gfni_r0
