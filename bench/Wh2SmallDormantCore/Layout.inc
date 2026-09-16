#include <cstdio>
#include <cstdlib>
#include <type_traits>

namespace {
// Model the original definitions, not a typed view of any live public handle.
struct PriorPublic {
    explicit PriorPublic(uint8_t k = 0) : SmallK(k) {}
    wirehair_v2::Codec Impl;
    uint64_t MessageBytes = 0u;
    uint32_t BlockBytes = 0u;
    CodecMode Mode = CodecMode::Encoder;
    bool Decoded = false;
    const uint8_t SmallK;
    EncoderSourceState SourceState = EncoderSourceState::Invalid;
    const uint8_t* BorrowedSource = nullptr;
};
template<unsigned K> struct PriorSmall : PriorPublic {
    PriorSmall() : PriorPublic(K) {}
    std::unique_ptr<uint8_t[]> SmallBasis;
    std::unique_ptr<wirehair_small_core::Encoder<K>> SmallEncoder;
    std::unique_ptr<wirehair_small_core::Decoder<K>> SmallDecoder;
};
void LayoutCheck(bool ok)
{ if (!ok) { std::fputs("FAIL layout\n",stderr); std::abort(); } }
std::size_t Offset(const void* object, const void* member)
{ return reinterpret_cast<std::uintptr_t>(member)-reinterpret_cast<std::uintptr_t>(object); }
template<class A, class B> void HeaderLayout(const A& a, const B& b)
{
    LayoutCheck(sizeof(A) == sizeof(B) && alignof(A) == alignof(B));
#define SLOT(M) LayoutCheck(Offset(&a,&a.M) == Offset(&b,&b.M))
    SLOT(Impl); SLOT(MessageBytes); SLOT(BlockBytes); SLOT(Mode); SLOT(Decoded);
    SLOT(SmallK); SLOT(SourceState); SLOT(BorrowedSource);
#undef SLOT
    LayoutCheck(Offset(&a,&a.Impl) == 0);
    static_assert(!std::is_copy_constructible<A>::value && !std::is_copy_assignable<A>::value &&
                  !std::is_move_constructible<A>::value && !std::is_move_assignable<A>::value,
                  "Opaque codec lifetime must not be copied or moved");
}
template<unsigned K> void SmallLayout()
{
    SmallPublicCodec<K> a;
    PriorSmall<K> b;
    HeaderLayout(a,b);
    LayoutCheck(Offset(&a,&a.SmallBasis) == Offset(&b,&b.SmallBasis) &&
                Offset(&a,&a.SmallEncoder) == Offset(&b,&b.SmallEncoder) &&
                Offset(&a,&a.SmallDecoder) == Offset(&b,&b.SmallDecoder));
}
}
int main()
{
    PublicCodec a; PriorPublic b;
    HeaderLayout(a,b);
    SmallLayout<3>(); SmallLayout<5>(); SmallLayout<8>();
    LayoutCheck(sizeof(wirehair_v2::Codec) == 240 && sizeof(PublicCodec) == 272 &&
                sizeof(SmallPublicCodec<3>) == 296);
    std::puts("PASS core240 certified272 small296; all sizes/alignments/member offsets and noncopyability match");
}
