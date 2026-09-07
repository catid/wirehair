#include "Wh2ThueMorseTinyGfniR0.h"

// Preserve the established test allocation boundary when main is wrapped.
__attribute__((noinline)) void operator delete(void*) noexcept;

#if !defined(WH2_TINY_GFNI_CANDIDATE) || \
    (WH2_TINY_GFNI_CANDIDATE != 0 && WH2_TINY_GFNI_CANDIDATE != 1)
#error "Compile with WH2_TINY_GFNI_CANDIDATE=0 or 1"
#endif
#if WH2_TINY_GFNI_CANDIDATE
#define wh2_thue_native wh2_thue_tiny_gfni_r0
#endif
#define main Wh2OriginalNeutralMain
#include "Wh2ThueMorseNativeCodecTest.cpp"
#undef main
#if WH2_TINY_GFNI_CANDIDATE
#undef wh2_thue_native
#endif

int main(int argc, char** argv)
{
    if (argc != 2 || (std::strcmp(argv[1], "--require-gfni") != 0 &&
                      std::strcmp(argv[1], "--expect-scalar") != 0)) {
        std::cerr << "usage: --require-gfni | --expect-scalar\n";
        return 2;
    }
    Check(gf256_init() == 0 && GF256Ctx.Polynomial == 0x14d, "shared ordinary GF256");
    Check(wh2_thue_tiny_gfni_r0::Available() == (std::strcmp(argv[1], "--require-gfni") == 0),
          "requested native/scalar codec path");
    return Wh2OriginalNeutralMain();
}
