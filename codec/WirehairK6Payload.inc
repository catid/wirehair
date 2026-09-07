// Included only by gf256.cpp: use its initialized capability gates and tables.
#include "WirehairK6Payload.h"

namespace wirehair_k6_payload {

void Scalar(void* output, const void* const* sources, const std::uint8_t* scales)
{
    uint8_t first = 0, second = 0;
    for (unsigned source = 0; source < 6; ++source) {
        const uint8_t* input = static_cast<const uint8_t*>(sources[source]);
        first ^= gf256_mul(input[0], scales[source]);
        second ^= gf256_mul(input[1], scales[source]);
    }
    uint8_t* destination = static_cast<uint8_t*>(output);
    destination[0] = first;
    destination[1] = second;
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

} // namespace wirehair_k6_payload
