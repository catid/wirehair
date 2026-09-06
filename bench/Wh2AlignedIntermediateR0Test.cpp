// Neutral ownership/correctness only.  The same source runs in separate
// control and candidate mirrors; it contains no clock or benchmark campaign.
#include "Wh2AlignedIntermediateR0.h"
#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2PrecodeDecode.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

#if !defined(WIREHAIR_TESTING) || !defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#error "Every test and library translation unit must enable both test hooks"
#endif
#if !defined(WH2_ALIGNMENT_CANDIDATE) || \
    (WH2_ALIGNMENT_CANDIDATE != 0 && WH2_ALIGNMENT_CANDIDATE != 1)
#error "Select exactly one control or candidate source mirror"
#endif

namespace {

using wirehair_v2::PrecodeEncoder;
using wirehair_v2::PrecodeSystem;
using wirehair_v2::PacketRowConfig;
using wh2_aligned_intermediate_r0::Allocator;
using wh2_aligned_intermediate_r0::ByteVector;

void Check(bool condition, const char* expression, unsigned line)
{
    if (!condition) {
        std::fprintf(stderr, "line %u: %s\n", line, expression);
        throw std::runtime_error("neutral assertion failed");
    }
}
#define CHECK(expression) Check((expression), #expression, __LINE__)

struct Digest
{
    uint64_t State = UINT64_C(14695981039346656037);
    void Bytes(const uint8_t* data, size_t size)
    {
        for (size_t i = 0; i < size; ++i) {
            State = (State ^ data[i]) * UINT64_C(1099511628211);
        }
    }
    void Integer(uint64_t value)
    {
        for (unsigned i = 0; i < 8u; ++i) {
            const uint8_t byte = static_cast<uint8_t>(value >> (8u * i));
            Bytes(&byte, 1u);
        }
    }
};

size_t CallocCalls = 0u;
size_t CallocCount = 0u;
size_t CallocSize = 0u;
bool RejectCalloc = false;

void* CheckedCalloc(size_t count, size_t size)
{
    ++CallocCalls;
    CallocCount = count;
    CallocSize = size;
    // Boundary tests must never attempt a huge real allocation.
    if (RejectCalloc || count != 1u || size > 1024u * 1024u) {
        return nullptr;
    }
    return std::calloc(count, size);
}

struct CallocScope
{
    explicit CallocScope(bool reject)
    {
        CallocCalls = CallocCount = CallocSize = 0u;
        RejectCalloc = reject;
        wirehair::detail::SetSIMDSafeCallocForTesting(CheckedCalloc);
    }
    ~CallocScope()
    {
        wirehair::detail::SetSIMDSafeCallocForTesting(nullptr);
    }
    CallocScope(const CallocScope&) = delete;
    CallocScope& operator=(const CallocScope&) = delete;
};

struct ValidationFailureScope
{
    explicit ValidationFailureScope(
        wirehair_v2::test::EncodeAllocationFailureException exception)
    {
        wirehair_v2::test::SetEncodeAllocationFailurePointForTesting(
            wirehair_v2::test::EncodeAllocationFailurePoint::InitializeSolvedValidation,
            exception);
    }
    ~ValidationFailureScope()
    {
        wirehair_v2::test::SetEncodeAllocationFailurePointForTesting(
            wirehair_v2::test::EncodeAllocationFailurePoint::None,
            wirehair_v2::test::EncodeAllocationFailureException::BadAlloc);
    }
    ValidationFailureScope(const ValidationFailureScope&) = delete;
    ValidationFailureScope& operator=(const ValidationFailureScope&) = delete;
};

bool Aligned(const void* address)
{
    return address && reinterpret_cast<uintptr_t>(address) % 32u == 0u;
}

template<class T> void RejectAllocation(Allocator<T>& allocator, size_t count)
{
    bool failed = false;
    try {
        T* pointer = allocator.allocate(count);
        allocator.deallocate(pointer, count);
    }
    catch (const std::bad_alloc&) {
        failed = true;
    }
    CHECK(failed);
}

void AllocatorCases()
{
    Allocator<uint8_t> allocator;
    {
        CallocScope calls(false);
        CHECK(allocator.allocate(0u) == nullptr);
        allocator.deallocate(nullptr, 0u);
        CHECK(CallocCalls == 0u);
    }
    const size_t sizes[] = {1u, 15u, 16u, 17u, 31u, 32u, 33u,
                            63u, 64u, 65u, 127u, 128u, 129u};
    for (size_t size : sizes)
    {
        CallocScope calls(false);
        uint8_t* pointer = allocator.allocate(size);
        CHECK(Aligned(pointer));
        CHECK(CallocCalls == 1u && CallocCount == 1u);
        CHECK(CallocSize == size + 32u);
        for (size_t j = 0; j < size; ++j) {
            CHECK(pointer[j] == 0u);
            pointer[j] = static_cast<uint8_t>(17u * j + size);
        }
        allocator.deallocate(pointer, size);
    }
    {
        CallocScope calls(true);
        RejectAllocation(allocator, 33u);
        CHECK(CallocCalls == 1u && CallocSize == 65u);
    }
    {
        CallocScope calls(true);
        RejectAllocation(allocator, allocator.max_size() + 1u);
        RejectAllocation(allocator, std::numeric_limits<size_t>::max());
        CHECK(CallocCalls == 0u);
        RejectAllocation(allocator, allocator.max_size());
        CHECK(CallocCalls == 1u && CallocCount == 1u);
        CHECK(CallocSize == allocator.max_size() + 32u);
    }
    Allocator<uint64_t> wide;
    CHECK(wide == allocator && !(wide != allocator));
    {
        CallocScope calls(true);
        RejectAllocation(wide, wide.max_size() + 1u);
        CHECK(CallocCalls == 0u);
        RejectAllocation(wide, wide.max_size());
        CHECK(CallocCalls == 1u);
        CHECK(CallocSize == wide.max_size() * sizeof(uint64_t) + 32u);
    }
    {
        CallocScope calls(false);
        uint64_t* pointer = wide.allocate(3u);
        CHECK(Aligned(pointer) && CallocSize == 3u * sizeof(uint64_t) + 32u);
        typedef std::allocator_traits<Allocator<uint64_t> > Traits;
        Traits::construct(wide, pointer, UINT64_C(1));
        Traits::construct(wide, pointer + 1u, UINT64_MAX);
        Traits::construct(wide, pointer + 2u, UINT64_C(3));
        CHECK(pointer[0] == 1u && pointer[1] == UINT64_MAX && pointer[2] == 3u);
        for (size_t i = 0; i < 3u; ++i) {
            Traits::destroy(wide, pointer + i);
        }
        wide.deallocate(pointer, 3u);
    }
    std::vector<uint8_t> expected(129u);
    for (size_t i = 0; i < expected.size(); ++i) {
        expected[i] = static_cast<uint8_t>(31u * i + 7u);
    }
    ByteVector surviving;
    {
        ByteVector first;
        {
            CallocScope calls(false);
            first.assign(expected.begin(), expected.end());
            CHECK(CallocCalls == 1u && CallocSize == expected.size() + 32u);
        }
        CHECK(Aligned(first.data()));
        ByteVector copied(first);
        CHECK(copied == first && copied.data() != first.data());
        CHECK(Aligned(copied.data()));
        ByteVector assigned(5u, 0x55u);
        assigned = first;
        CHECK(assigned == first && assigned.data() != first.data());
        CHECK(Aligned(assigned.data()));
        ByteVector moved(std::move(copied));
        CHECK(moved == first && Aligned(moved.data()));
        ByteVector moved_to(3u, 0xaau);
        moved_to = std::move(assigned);
        CHECK(moved_to == first && Aligned(moved_to.data()));
        ByteVector different(17u, 0x5au);
        moved.swap(different);
        CHECK(moved == ByteVector(17u, 0x5au));
        CHECK(different == first && Aligned(different.data()));
        surviving = different;
        first[0] ^= 0xffu;
        CHECK(surviving[0] == expected[0]);
        const ByteVector before(surviving);
        const uint8_t* old_data = surviving.data();
        {
            CallocScope calls(true);
            bool failed = false;
            try { surviving.reserve(surviving.capacity() + 1u); }
            catch (const std::bad_alloc&) { failed = true; }
            CHECK(failed && CallocCalls == 1u);
        }
        CHECK(surviving == before && surviving.data() == old_data);
    }
    CHECK(std::equal(surviving.begin(), surviving.end(), expected.begin()));
    CHECK(Aligned(surviving.data()));
}

// Only this neutral test obtains the exact private member pointer.  No layout
// offset, facade replacement, or production class definition is changed.
struct SolvedTag
{
    typedef WirehairResult (PrecodeEncoder::*Type)(
        const PrecodeSystem&, const PacketRowConfig&,
        std::vector<uint8_t>&, uint32_t);
    friend Type SolvedMember(SolvedTag);
};
template<class Tag, typename Tag::Type Member> struct MemberAccess
{
    friend typename Tag::Type SolvedMember(Tag) { return Member; }
};
template struct MemberAccess<SolvedTag, &PrecodeEncoder::InitializeSolvedSystem>;

struct Fixture
{
    uint32_t K;
    uint32_t B;
    std::vector<uint8_t> Message;
    std::vector<uint8_t> Intermediate;
    wirehair_v2::SeedProfile Profile;
    PacketRowConfig Config;
};

void AddProfile(Digest& digest, const wirehair_v2::SeedProfile& profile)
{
    // Explicit contract integers avoid struct padding and address-dependent data.
    const uint64_t fields[] = {
        profile.BlockCount, profile.BlockBytes, profile.DenseCount,
        profile.PeelSeed, profile.DenseSeed, profile.PeelSeedBucket,
        profile.UsedPeelFixup, profile.UsedDenseFixup, profile.V2SeedSelected,
        profile.V2SeedAttempt, profile.V2PrecodeContractVersion,
        profile.V2PacketRowContractVersion, profile.V2StaircaseCount,
        profile.V2DenseRowCount, profile.V2HeavyRowCount, profile.V2SourceHits,
        profile.V2PrecodeSeed, profile.V2PacketPeelSeed,
        profile.V2RecoveryMixCount, profile.V2DenseAnchorLayout,
        profile.V2DenseIdentityCorner, profile.V2PrecodeSeedSalt,
        profile.V2RecoveryRowSeedSalt
    };
    for (uint64_t value : fields) { digest.Integer(value); }
}

std::vector<uint8_t> ScalarPacket(const Fixture& fixture, uint32_t id,
                                uint64_t& operations)
{
    const uint32_t columns = static_cast<uint32_t>(fixture.Intermediate.size() /
                                                    fixture.B);
    const std::vector<uint32_t> row = wirehair_v2::GeneratePacketMatrixRow(
        fixture.K, columns - fixture.K, id, fixture.Config);
    CHECK(!row.empty());
    std::vector<uint8_t> output(fixture.B, 0u);
    for (uint32_t column : row) {
        CHECK(column < columns);
        for (uint32_t j = 0; j < fixture.B; ++j) {
            output[j] ^= fixture.Intermediate[static_cast<size_t>(column) *
                                               fixture.B + j];
        }
    }
    operations = row.size();
    return output;
}

std::vector<uint8_t> CheckBlock(const PrecodeEncoder& encoder,
                               const Fixture& fixture, uint32_t id,
                               Digest& digest)
{
    uint64_t expected_ops = 0u;
    const std::vector<uint8_t> expected = ScalarPacket(fixture, id, expected_ops);
    std::vector<uint8_t> guarded(fixture.B + 64u, 0xa5u);
    uint64_t operations = UINT64_MAX;
    CHECK(encoder.EncodeResult(id, guarded.data() + 32u, &operations) ==
          Wirehair_Success);
    CHECK(operations == expected_ops);
    CHECK(std::equal(expected.begin(), expected.end(), guarded.begin() + 32u));
    CHECK(std::all_of(guarded.begin(), guarded.begin() + 32u,
                     [](uint8_t value) { return value == 0xa5u; }));
    CHECK(std::all_of(guarded.begin() + 32u + fixture.B, guarded.end(),
                     [](uint8_t value) { return value == 0xa5u; }));
    if (id < fixture.K) {
        for (uint32_t j = 0; j < fixture.B; ++j) {
            const size_t offset = static_cast<size_t>(id) * fixture.B + j;
            CHECK(expected[j] == (offset < fixture.Message.size() ?
                                  fixture.Message[offset] : 0u));
        }
    }
    digest.Integer(id);
    digest.Integer(operations);
    digest.Bytes(expected.data(), expected.size());
    return expected;
}

void CheckEncoder(const PrecodeEncoder& encoder, const Fixture& fixture,
                  Digest& digest)
{
    CHECK(encoder.IsInitialized() && !encoder.HasCompleteSystem());
    CHECK(encoder.SourceBlockCount() == fixture.K && encoder.BlockBytes() == fixture.B);
    CHECK(encoder.RecoveryRowSeed() == fixture.Config.PeelSeed);
    CHECK(encoder.RecoveryMixCount() == fixture.Config.MixCount);
    CHECK(encoder.IntermediateBlocks() != nullptr);
    CHECK(std::memcmp(encoder.IntermediateBlocks(), fixture.Intermediate.data(),
                      fixture.Intermediate.size()) == 0);
    CHECK(encoder.ParityBlocks() == encoder.IntermediateBlocks() +
          static_cast<size_t>(fixture.K) * fixture.B);
#if WH2_ALIGNMENT_CANDIDATE
    CHECK(Aligned(encoder.IntermediateBlocks()));
#endif
    for (uint32_t id = 0; id < fixture.K; ++id) {
        CheckBlock(encoder, fixture, id, digest);
    }
    const uint32_t repairs[] = {fixture.K, fixture.K + 1u, fixture.K + 9u,
                               UINT32_C(0x7fffffff), UINT32_MAX};
    for (uint32_t id : repairs) { CheckBlock(encoder, fixture, id, digest); }
}

void TransactionCase(PrecodeEncoder& encoder, const Fixture& fixture,
                     Digest& digest)
{
    PrecodeSystem system;
    CHECK(wirehair_v2::BuildPrecodeSystem(encoder.System().Params, system));
    std::vector<uint8_t> input = fixture.Intermediate;
    const uint8_t* input_address = input.data();
    const size_t input_capacity = input.capacity();
    const uint8_t* old_address = encoder.IntermediateBlocks();
    const SolvedTag::Type initialize = SolvedMember(SolvedTag());
    const wirehair_v2::test::EncodeAllocationFailureException exceptions[] = {
        wirehair_v2::test::EncodeAllocationFailureException::BadAlloc,
        wirehair_v2::test::EncodeAllocationFailureException::LengthError
    };
    for (auto exception : exceptions) {
        ValidationFailureScope failure(exception);
        CHECK((encoder.*initialize)(system, fixture.Config, input, fixture.B) ==
              Wirehair_OOM);
        CHECK(wirehair_v2::test::EncodeAllocationFailureHitsForTesting() == 1u);
        CHECK(input == fixture.Intermediate && input.data() == input_address);
        CHECK(input.capacity() == input_capacity);
        CHECK(encoder.IntermediateBlocks() == old_address);
    }
#if WH2_ALIGNMENT_CANDIDATE
    {
        CallocScope calls(true);
        CHECK((encoder.*initialize)(system, fixture.Config, input, fixture.B) ==
              Wirehair_OOM);
        CHECK(CallocCalls == 1u && CallocCount == 1u);
        CHECK(CallocSize == input.size() + 32u);
    }
    CHECK(input == fixture.Intermediate && input.data() == input_address);
    CHECK(input.capacity() == input_capacity);
    CHECK(encoder.IntermediateBlocks() == old_address);
    CheckEncoder(encoder, fixture, digest);
#else
    // Control consumes the input by ordinary vector swap and never calls calloc.
    // Equal digest work on both variants; only candidate injects its new failure.
    CheckEncoder(encoder, fixture, digest);
#endif
    {
        CallocScope calls(WH2_ALIGNMENT_CANDIDATE == 0);
        CHECK((encoder.*initialize)(system, fixture.Config, input, fixture.B) ==
              Wirehair_Success);
        CHECK(CallocCalls == static_cast<size_t>(WH2_ALIGNMENT_CANDIDATE));
#if WH2_ALIGNMENT_CANDIDATE
        CHECK(CallocCount == 1u && CallocSize == fixture.Intermediate.size() + 32u);
#endif
    }
    CHECK(input.empty() && input.capacity() == 0u);
    CheckEncoder(encoder, fixture, digest);
}

void RoundTrip(const PrecodeEncoder& encoder, const Fixture& fixture,
               Digest& digest)
{
    wirehair_v2::MessagePrecodeDecoder decoder;
    CHECK(decoder.InitializeResult(fixture.Message.size(), fixture.B,
                                   &fixture.Profile) == Wirehair_Success);
    std::vector<uint8_t> recovered(fixture.Message.size() + 64u, 0xa5u);
    CHECK(decoder.RecoverResult(recovered.data() + 32u, fixture.Message.size()) ==
          Wirehair_NeedMore);
    CHECK(std::all_of(recovered.begin(), recovered.end(),
                     [](uint8_t value) { return value == 0xa5u; }));
    // Fixed repair-first schedule exercises the unchanged solver, then includes
    // every systematic id if needed.  This is not a recovery-frequency screen.
    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t step = 0; step < fixture.K + 6u && result != Wirehair_Success;
         ++step)
    {
        const uint32_t id = step < 6u ? fixture.K + step : step - 6u;
        const std::vector<uint8_t> packet = CheckBlock(encoder, fixture, id, digest);
        const uint32_t bytes = id + 1u == fixture.K ?
            static_cast<uint32_t>(fixture.Message.size() -
                static_cast<size_t>(fixture.K - 1u) * fixture.B) : fixture.B;
        result = decoder.DecodeResult(id, packet.data(), bytes);
        CHECK(result == Wirehair_NeedMore || result == Wirehair_Success);
        digest.Integer(static_cast<uint64_t>(result));
    }
    CHECK(result == Wirehair_Success);
    CHECK(decoder.RecoverResult(recovered.data() + 32u, fixture.Message.size()) ==
          Wirehair_Success);
    CHECK(std::equal(fixture.Message.begin(), fixture.Message.end(),
                     recovered.begin() + 32u));
    CHECK(std::all_of(recovered.begin(), recovered.begin() + 32u,
                     [](uint8_t value) { return value == 0xa5u; }));
    CHECK(std::all_of(recovered.end() - 32u, recovered.end(),
                     [](uint8_t value) { return value == 0xa5u; }));
    digest.Bytes(recovered.data() + 32u, fixture.Message.size());
}

void ShapeCase(uint32_t K, uint32_t B, bool partial, Digest& digest)
{
    Fixture fixture;
    fixture.K = K;
    fixture.B = B;
    fixture.Message.resize(static_cast<size_t>(K) * B - (partial ? 1u : 0u));
    for (size_t i = 0; i < fixture.Message.size(); ++i) {
        fixture.Message[i] = static_cast<uint8_t>(i * 67u + K * 13u + B + 5u);
    }
    std::unique_ptr<PrecodeEncoder> copied;
    PrecodeEncoder assigned;
    {
        std::vector<uint8_t> source(fixture.Message.size() + 64u, 0xa5u);
        std::copy(fixture.Message.begin(), fixture.Message.end(), source.begin() + 32u);
        const std::vector<uint8_t> before = source;
        wirehair_v2::MessagePrecodeEncoder owner;
        CHECK(owner.InitializeResult(source.data() + 32u, fixture.Message.size(), B) ==
              Wirehair_Success);
        CHECK(source == before && !owner.HasSystematicSourceCache());
        fixture.Profile = owner.Profile();
        CHECK(fixture.Profile.V2SeedSelected && fixture.Profile.BlockCount == K);
        fixture.Config.PeelSeed = static_cast<uint32_t>(owner.BlockEncoder().RecoveryRowSeed());
        fixture.Config.MixCount = owner.BlockEncoder().RecoveryMixCount();
        const size_t size = static_cast<size_t>(K + owner.BlockEncoder().ParityBlockCount()) * B;
        fixture.Intermediate.assign(owner.IntermediateBlocks(), owner.IntermediateBlocks() + size);
        AddProfile(digest, fixture.Profile);
        digest.Integer(fixture.Message.size());
        digest.Bytes(fixture.Intermediate.data(), fixture.Intermediate.size());
        // Destroy the caller's data content before exercising message encoding.
        std::fill(source.begin() + 32u, source.end() - 32u, 0x3cu);
        for (uint32_t id = 0; id <= K; ++id)
        {
            uint64_t expected_ops = 0u;
            const std::vector<uint8_t> expected = ScalarPacket(fixture, id, expected_ops);
            const uint32_t length = id + 1u == K && partial ? B - 1u : B;
            std::vector<uint8_t> output(B + 64u, 0xa5u);
            uint32_t bytes = UINT32_MAX;
            uint64_t operations = UINT64_MAX;
            CHECK(owner.EncodeResult(id, output.data() + 32u, length,
                                     &bytes, &operations) == Wirehair_Success);
            CHECK(bytes == length && operations == expected_ops);
            CHECK(std::equal(expected.begin(), expected.begin() + length,
                             output.begin() + 32u));
            CHECK(std::all_of(output.begin(), output.begin() + 32u,
                             [](uint8_t value) { return value == 0xa5u; }));
            CHECK(std::all_of(output.begin() + 32u + length, output.end(),
                             [](uint8_t value) { return value == 0xa5u; }));
            const std::vector<uint8_t> successful = output;
            bytes = UINT32_MAX;
            operations = UINT64_MAX;
            CHECK(owner.EncodeResult(id, output.data() + 32u, length - 1u,
                                     &bytes, &operations) == Wirehair_InvalidInput);
            CHECK(output == successful && bytes == UINT32_MAX && operations == UINT64_MAX);
            digest.Bytes(successful.data() + 32u, length);
        }
        copied.reset(new PrecodeEncoder(owner.BlockEncoder()));
        assigned = owner.BlockEncoder();
        CHECK(copied->IntermediateBlocks() != owner.IntermediateBlocks());
        CHECK(assigned.IntermediateBlocks() != owner.IntermediateBlocks());
    }
    // The owning MessagePrecodeEncoder and original caller buffer are now gone.
    CheckEncoder(*copied, fixture, digest);
    CheckEncoder(assigned, fixture, digest);
    PrecodeEncoder moved(std::move(*copied));
    copied.reset(); // No use of moved-from initialized encoder is assumed.
    CheckEncoder(moved, fixture, digest);
    PrecodeEncoder moved_to;
    moved_to = std::move(assigned);
    CheckEncoder(moved_to, fixture, digest);
    PrecodeEncoder overwritten;
    overwritten = moved;
    overwritten = moved_to;
    CheckEncoder(overwritten, fixture, digest);
    std::swap(moved, moved_to);
    CheckEncoder(moved, fixture, digest);
    TransactionCase(overwritten, fixture, digest);
    RoundTrip(moved, fixture, digest);
}

} // namespace

int main()
{
    try {
        CHECK(gf256_init() == 0);
        AllocatorCases();
        Digest digest;
        unsigned shapes = 0u;
        const uint32_t widths[] = {1u, 15u, 16u, 17u, 31u, 32u, 33u,
                                   63u, 64u, 65u, 127u, 128u, 129u};
        const uint32_t counts[] = {7u, 9u};
        for (uint32_t K : counts) {
            for (uint32_t B : widths) {
                ShapeCase(K, B, false, digest);
                ++shapes;
                if (B > 1u) {
                    ShapeCase(K, B, true, digest);
                    ++shapes;
                }
            }
        }
        CHECK(shapes == 50u);
        // Identical stdout is the cross-executable semantic comparison.  No
        // address/alignment/hook counts enter this digest.
        std::printf("PASS shapes=%u semantic_fnv1a64=%016llx\n", shapes,
                    static_cast<unsigned long long>(digest.State));
        return 0;
    }
    catch (const std::exception& error) {
        std::fprintf(stderr, "FAIL: %s\n", error.what());
        return 1;
    }
}
