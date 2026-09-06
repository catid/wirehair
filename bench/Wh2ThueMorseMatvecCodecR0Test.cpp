#include "Wh2ThueMorseMatvecGfniR0.h"

// Keep the test's malloc-backed new/delete pair at a delete call boundary.
// Otherwise GCC can inline free while treating the matching replacement new
// as an opaque builtin allocation and report a mismatched-new-delete warning.
__attribute__((noinline)) void operator delete(void*) noexcept;

#if !defined(WH2_MATVEC_CODEC_CANDIDATE) || \
    (WH2_MATVEC_CODEC_CANDIDATE != 0 && WH2_MATVEC_CODEC_CANDIDATE != 1)
#error "Compile this wrapper with WH2_MATVEC_CODEC_CANDIDATE=0 or 1"
#endif

#if WH2_MATVEC_CODEC_CANDIDATE
#define wh2_thue_native wh2_thue_matvec_codec_r0
#endif
#define main Wh2OriginalNeutralMain
#include "Wh2ThueMorseNativeCodecTest.cpp"
#undef main
#if WH2_MATVEC_CODEC_CANDIDATE
#undef wh2_thue_native
#endif

namespace {
// All oracle buffers and codec allocations precede these guarded calls.
// Even an unexpected exception disables allocation tracking before unwinding.
template<class Operation>
Result DistantNoAllocation(Operation operation)
{
    StartTracking(0);
    Result result;
    try {
        result = operation();
    } catch (...) {
        track_allocations = false;
        throw;
    }
    track_allocations = false;
    Check(allocation_count == 0, "distant-first Encode/Feed/Recover must not allocate");
    return result;
}

void DistantResult(const Result& result, Status status,
                   std::size_t required, std::size_t written)
{
    Check(result.status == status && result.bytes_required == required &&
          result.bytes_written == written, "distant-first exact result/lengths");
}

void DistantFirstCase()
{
    const std::uint32_t block = 64u;
    const std::size_t message_bytes = 383u; // K6, tail63.
    const std::uint32_t distant_id = UINT32_C(0xffffffff);
    const Neutral neutral;
    const std::vector<Byte> message = Message(message_bytes);
    const auto saved_source = message;
    const auto saved_lookup = neutral.table;
    const Vector far_row = neutral.Reference(distant_id);
    Check(std::any_of(far_row.begin(), far_row.end(), [](Byte value) { return value != 0; }),
          "independent distant-first row is nonzero");

    std::unique_ptr<Encoder> encoder;
    std::unique_ptr<Decoder> decoder;
    Check(Encoder::Create(neutral.View(), message.data(), message.size(), block, encoder) == Status::Success,
          "distant-first partial encoder precreated");
    Check(Decoder::Create(neutral.View(), message.size(), block, decoder) == Status::Success &&
          decoder->Rank() == 0, "distant-first decoder precreated empty");

    typedef std::array<Byte, 64u + 34u> GuardedPacket;
    GuardedPacket distant, expected_distant;
    distant.fill(0xa7);
    expected_distant.fill(0xa7);
    const std::vector<Byte> far_reference = PacketReference(neutral, message, block, distant_id);
    Check(far_reference.size() == block, "distant packet never uses truncated systematic tail");
    std::copy(far_reference.begin(), far_reference.end(), expected_distant.begin() + 17);

    std::array<GuardedPacket, 6> packets, expected_packets;
    for (unsigned id = 0; id < 6; ++id) {
        packets[id].fill(0xa7);
        expected_packets[id].fill(0xa7);
        const std::vector<Byte> reference = PacketReference(neutral, message, block, id);
        std::copy(reference.begin(), reference.end(), expected_packets[id].begin() + 17);
        DistantResult(encoder->Encode(id, packets[id].data() + 17, reference.size()),
                      Status::Success, reference.size(), reference.size());
        Check(packets[id] == expected_packets[id], "systematic preflight oracle and packet guards");
    }
    GuardedPacket conflicting = expected_distant;
    conflicting[17u + block / 2u] ^= 1u;
    const GuardedPacket expected_conflict = conflicting;
    std::vector<Byte> recovered(message_bytes + 34u, 0xa7);
    std::vector<Byte> expected_recovered = recovered;
    std::copy(message.begin(), message.end(), expected_recovered.begin() + 17);

    DistantResult(DistantNoAllocation([&] {
        return encoder->Encode(distant_id, distant.data() + 17, block);
    }), Status::Success, block, block);
    Check(distant == expected_distant, "distant Encode independent payload and exact guards");
    const auto immutable = [&] {
        Check(message == saved_source && neutral.table == saved_lookup &&
              packets == expected_packets && distant == expected_distant &&
              conflicting == expected_conflict, "distant-first source/lookup/packet immutability");
    };
    immutable();

    DistantResult(DistantNoAllocation([&] {
        return decoder->Feed(distant_id, distant.data() + 17, block);
    }), Status::NeedMore, block, 0);
    Check(decoder->Rank() == 1, "distant-first equation creates rank one");
    immutable();
    DistantResult(DistantNoAllocation([&] {
        return decoder->Feed(distant_id, distant.data() + 17, block);
    }), Status::NeedMore, block, 0);
    Check(decoder->Rank() == 1, "distant duplicate preserves rank one");
    immutable();
    DistantResult(DistantNoAllocation([&] {
        return decoder->Feed(distant_id, conflicting.data() + 17, block);
    }), Status::Conflict, block, 0);
    Check(decoder->Rank() == 1, "distant conflicting duplicate preserves rank one");
    immutable();
    DistantResult(DistantNoAllocation([&] {
        return decoder->Recover(recovered.data() + 17, message_bytes);
    }), Status::NeedMore, message_bytes, 0);
    Check(std::all_of(recovered.begin(), recovered.end(), [](Byte value) { return value == 0xa7; }),
          "rank-one Recover leaves output and guards untouched");
    immutable();

    for (unsigned id = 0; id < 6; ++id) {
        const unsigned prefix = id + 1u;
        const bool distant_adds_rank = std::any_of(far_row.begin() + prefix, far_row.end(),
            [](Byte value) { return value != 0; });
        const unsigned expected_rank = prefix + (distant_adds_rank ? 1u : 0u);
        const Status expected_status = expected_rank == 6u ? Status::Success : Status::NeedMore;
        const std::size_t bytes = id == 5u ? 63u : block;
        DistantResult(DistantNoAllocation([&] {
            return decoder->Feed(id, packets[id].data() + 17, bytes);
        }), expected_status, bytes, 0);
        Check(decoder->Rank() == expected_rank, "systematic-prefix rank matches independent far-row support");
        immutable();
        std::fill(recovered.begin(), recovered.end(), 0xa7);
        DistantResult(DistantNoAllocation([&] {
            return decoder->Recover(recovered.data() + 17, message_bytes);
        }), expected_status, message_bytes, expected_rank == 6u ? message_bytes : 0u);
        if (expected_rank == 6u) {
            Check(recovered == expected_recovered, "distant-first exact eventual recovery and guards");
        } else {
            Check(std::all_of(recovered.begin(), recovered.end(), [](Byte value) { return value == 0xa7; }),
                  "incomplete systematic-prefix Recover leaves all bytes untouched");
        }
        immutable();
    }
    Check(decoder->Rank() == 6u, "complete distant-first/systematic neutral case");
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2 || (std::strcmp(argv[1], "--require-gfni") != 0 &&
                      std::strcmp(argv[1], "--expect-scalar") != 0)) {
        std::cerr << "usage: --require-gfni | --expect-scalar\n";
        return 2;
    }
    Check(gf256_init() == 0 && GF256Ctx.Polynomial == 0x14d, "qualified shared GF256 initialization");
    Check(wh2_matvec_gfni_r0::Available() == (std::strcmp(argv[1], "--require-gfni") == 0),
          "requested GFNI/scalar codec availability");
    const int result = Wh2OriginalNeutralMain();
    if (result != 0) return result;
    DistantFirstCase();
    std::cout << "PASS one added B64/tail63 distant-first case; 70 total neutral shapes; no selected data/timing\n";
    return 0;
}
