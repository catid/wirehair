#include <wirehair/wirehair.h>

#include "../WirehairCodec.h"
#include "../WirehairTools.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
#endif

namespace {

static const uint16_t kBlockCount = 2048;
static const uint32_t kBlockBytes = 17;
static const uint16_t kPathHead = 0;
static const size_t kSmallStackBytes = 64u * 1024u;

struct AdversarialPath
{
    std::vector<uint32_t> EdgeIds;
    std::vector<uint16_t> Columns;
    uint32_t RootId = 0;
    uint16_t Tail = 0;
};

uint64_t FnvByte(uint64_t hash, uint8_t value)
{
    return (hash ^ value) * UINT64_C(1099511628211);
}

uint64_t FnvU16(uint64_t hash, uint16_t value)
{
    hash = FnvByte(hash, static_cast<uint8_t>(value));
    return FnvByte(hash, static_cast<uint8_t>(value >> 8));
}

uint64_t FnvU32(uint64_t hash, uint32_t value)
{
    hash = FnvU16(hash, static_cast<uint16_t>(value));
    return FnvU16(hash, static_cast<uint16_t>(value >> 16));
}

bool BuildAdversarialPath(AdversarialPath& path)
{
    const uint16_t peel_seed = wirehair::GetPeelSeed(kBlockCount);
    const uint16_t mix_count = static_cast<uint16_t>(
        wirehair::GetDenseCount(kBlockCount) + wirehair::kHeavyRows);
    const uint16_t next_prime = wirehair::NextPrime16(kBlockCount);

    std::vector<uint8_t> seen(kBlockCount, 0);
    path.EdgeIds.reserve(kBlockCount - 1);
    path.Columns.reserve(kBlockCount);
    path.Columns.push_back(kPathHead);
    seen[kPathHead] = 1;
    uint16_t tail = kPathHead;

    uint32_t id = kBlockCount;
    for (; path.EdgeIds.size() + 1 < kBlockCount &&
           id != std::numeric_limits<uint32_t>::max(); ++id)
    {
        wirehair::PeelRowParameters params;
        params.Initialize(id, peel_seed, kBlockCount, mix_count);
        if (params.PeelCount != 2) {
            continue;
        }

        wirehair::PeelRowIterator iter(params, kBlockCount, next_prime);
        const uint16_t a = iter.GetColumn();
        if (!iter.Iterate()) {
            return false;
        }
        const uint16_t b = iter.GetColumn();
        if (iter.Iterate()) {
            return false;
        }

        uint16_t next = std::numeric_limits<uint16_t>::max();
        if (a == tail && !seen[b]) {
            next = b;
        }
        else if (b == tail && !seen[a]) {
            next = a;
        }
        if (next == std::numeric_limits<uint16_t>::max()) {
            continue;
        }

        path.EdgeIds.push_back(id);
        path.Columns.push_back(next);
        seen[next] = 1;
        tail = next;
    }

    if (path.EdgeIds.size() + 1 != kBlockCount) {
        std::fprintf(stderr, "failed to construct path: edges=%zu next_id=%u\n",
            path.EdgeIds.size(), id);
        return false;
    }

    for (; id != std::numeric_limits<uint32_t>::max(); ++id)
    {
        wirehair::PeelRowParameters params;
        params.Initialize(id, peel_seed, kBlockCount, mix_count);
        if (params.PeelCount == 1 && params.PeelFirst == kPathHead) {
            path.RootId = id;
            break;
        }
    }
    if (path.RootId == 0) {
        std::fprintf(stderr, "failed to find path root row\n");
        return false;
    }

    path.Tail = tail;
    for (uint8_t value : seen) {
        if (value != 1) {
            std::fprintf(stderr, "path does not visit every column once\n");
            return false;
        }
    }
    return true;
}

uint64_t PathFingerprint(const AdversarialPath& path)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (uint32_t id : path.EdgeIds) {
        hash = FnvU32(hash, id);
    }
    hash = FnvU32(hash, path.RootId);
    for (uint16_t column : path.Columns) {
        hash = FnvU16(hash, column);
    }
    return hash;
}

uint64_t ExpectedPeelOrder(const AdversarialPath& path)
{
    uint64_t hash = UINT64_C(14695981039346656037);

    // The root packet occupies decoder row N-1 and solves the path head.
    hash = FnvU16(hash, kBlockCount - 1);
    hash = FnvU16(hash, path.Columns[0]);

    // Recursive DFS then pauses each parent and walks the path in feed order.
    for (uint16_t row = 0; row + 1 < kBlockCount; ++row)
    {
        hash = FnvU16(hash, row);
        hash = FnvU16(hash, path.Columns[row + 1]);
    }
    return hash;
}

std::vector<uint8_t> MakeMessage()
{
    const size_t message_bytes =
        static_cast<size_t>(kBlockCount - 1) * kBlockBytes + 12;
    std::vector<uint8_t> message(message_bytes);
    uint32_t state = UINT32_C(0x6f2c91d5);
    for (size_t i = 0; i < message.size(); ++i)
    {
        state = state * UINT32_C(1664525) + UINT32_C(1013904223);
        message[i] = static_cast<uint8_t>(state >> 24);
    }
    return message;
}

struct TestContext
{
    const AdversarialPath* Path = nullptr;
    const std::vector<uint8_t>* Message = nullptr;
    const std::vector<uint8_t>* Packets = nullptr;
    const std::vector<uint8_t>* ExtraPacket = nullptr;
    uint64_t ExpectedOrderHash = 0;
    size_t ObservedStackBytes = 0;
    int Result = 1;
};

int RunAvalanche(TestContext& context)
{
    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, context.Message->size(), kBlockBytes);
    if (!decoder) {
        std::fprintf(stderr, "decoder creation failed\n");
        return 1;
    }

    WirehairResult result = Wirehair_NeedMore;
    for (size_t i = 0; i < context.Path->EdgeIds.size(); ++i)
    {
        result = wirehair_decode(
            decoder,
            context.Path->EdgeIds[i],
            context.Packets->data() + i * kBlockBytes,
            kBlockBytes);
        if (result != Wirehair_NeedMore)
        {
            std::fprintf(stderr, "edge %zu returned %d\n",
                i, static_cast<int>(result));
            wirehair_free(decoder);
            return 1;
        }
    }

    result = wirehair_decode(
        decoder,
        context.Path->RootId,
        context.Packets->data() +
            context.Path->EdgeIds.size() * kBlockBytes,
        kBlockBytes);
    if (result != Wirehair_Success)
    {
        std::fprintf(stderr, "root returned %d\n", static_cast<int>(result));
        wirehair_free(decoder);
        return 1;
    }

    const wirehair::Codec* internal =
        reinterpret_cast<const wirehair::Codec*>(decoder);
    if (internal->TestingPeelMaxDepth() != kBlockCount ||
        internal->TestingPeelOrderHash() != context.ExpectedOrderHash)
    {
        std::fprintf(stderr,
            "peel mismatch: depth=%u expected=%u order=%016llx expected=%016llx\n",
            internal->TestingPeelMaxDepth(), kBlockCount,
            static_cast<unsigned long long>(internal->TestingPeelOrderHash()),
            static_cast<unsigned long long>(context.ExpectedOrderHash));
        wirehair_free(decoder);
        return 1;
    }

    std::vector<uint8_t> recovered(context.Message->size());
    if (wirehair_recover(decoder, recovered.data(), recovered.size()) !=
            Wirehair_Success ||
        recovered != *context.Message)
    {
        std::fprintf(stderr, "recovered message mismatch\n");
        wirehair_free(decoder);
        return 1;
    }

    // A completed decoder must remain stable under a duplicate and an extra
    // valid repair row; this also catches accidental worklist state leakage.
    if (wirehair_decode(
            decoder,
            context.Path->RootId,
            context.Packets->data() +
                context.Path->EdgeIds.size() * kBlockBytes,
            kBlockBytes) != Wirehair_Success ||
        wirehair_decode(
            decoder,
            context.Path->RootId + 1,
            context.ExtraPacket->data(),
            kBlockBytes) != Wirehair_Success)
    {
        std::fprintf(stderr, "completed decoder rejected duplicate/extra row\n");
        wirehair_free(decoder);
        return 1;
    }

    std::fill(recovered.begin(), recovered.end(), uint8_t{0});
    if (wirehair_recover(decoder, recovered.data(), recovered.size()) !=
            Wirehair_Success ||
        recovered != *context.Message)
    {
        std::fprintf(stderr, "post-completion recovery changed\n");
        wirehair_free(decoder);
        return 1;
    }

    wirehair_free(decoder);
    return 0;
}

#if defined(__unix__) || defined(__APPLE__)
bool RecordCurrentStackSize(TestContext& context)
{
#if defined(__linux__)
    pthread_attr_t attributes;
    if (pthread_getattr_np(pthread_self(), &attributes) != 0) {
        std::fprintf(stderr, "pthread_getattr_np failed\n");
        return false;
    }
    size_t stack_bytes = 0;
    const int result = pthread_attr_getstacksize(&attributes, &stack_bytes);
    pthread_attr_destroy(&attributes);
    if (result != 0) {
        std::fprintf(stderr, "pthread_attr_getstacksize failed\n");
        return false;
    }
    context.ObservedStackBytes = stack_bytes;
#elif defined(__APPLE__)
    context.ObservedStackBytes = pthread_get_stacksize_np(pthread_self());
#else
    // pthread_attr_setstacksize() still enforces the requested bound on the
    // platforms covered above; some other pthread implementations do not
    // expose a portable current-thread stack query.
    context.ObservedStackBytes = kSmallStackBytes;
#endif
    if (context.ObservedStackBytes > kSmallStackBytes)
    {
        std::fprintf(stderr, "pthread stack too large: %zu > %zu\n",
            context.ObservedStackBytes, kSmallStackBytes);
        return false;
    }
    return true;
}

void* RunAvalancheThread(void* opaque)
{
    TestContext& context = *static_cast<TestContext*>(opaque);
    if (!RecordCurrentStackSize(context)) {
        context.Result = 1;
        return nullptr;
    }
    context.Result = RunAvalanche(context);
    return nullptr;
}

int RunOnSmallStack(TestContext& context)
{
    pthread_attr_t attributes;
    if (pthread_attr_init(&attributes) != 0) {
        std::fprintf(stderr, "pthread_attr_init failed\n");
        return 1;
    }
    if (pthread_attr_setstacksize(&attributes, kSmallStackBytes) != 0) {
        std::fprintf(stderr, "pthread_attr_setstacksize(%zu) failed\n",
            kSmallStackBytes);
        pthread_attr_destroy(&attributes);
        return 1;
    }

    pthread_t thread;
    const int create_result = pthread_create(
        &thread, &attributes, RunAvalancheThread, &context);
    pthread_attr_destroy(&attributes);
    if (create_result != 0) {
        std::fprintf(stderr, "pthread_create failed: %d\n", create_result);
        return 1;
    }
    if (pthread_join(thread, nullptr) != 0) {
        std::fprintf(stderr, "pthread_join failed\n");
        return 1;
    }
    return context.Result;
}
#endif

} // namespace

int main(int argc, char** argv)
{
    bool small_stack = false;
    if (argc == 2 && std::strcmp(argv[1], "--small-stack") == 0) {
        small_stack = true;
    }
    else if (argc == 2 && std::strcmp(argv[1], "--normal-stack") == 0) {
        small_stack = false;
    }
    else if (argc != 1) {
        std::fprintf(stderr, "Usage: %s [--small-stack|--normal-stack]\n", argv[0]);
        return 2;
    }

    if (wirehair_init() != Wirehair_Success) {
        std::fprintf(stderr, "wirehair_init failed\n");
        return 1;
    }

    AdversarialPath path;
    if (!BuildAdversarialPath(path)) {
        return 1;
    }

    const uint64_t path_hash = PathFingerprint(path);
    const uint64_t expected_order_hash = ExpectedPeelOrder(path);
    static const uint64_t kExpectedPathHash =
        UINT64_C(0x2c92c543d159be0e);
    static const uint64_t kExpectedOrderHash =
        UINT64_C(0xff84f157d47b6aad);
    if (path_hash != kExpectedPathHash ||
        expected_order_hash != kExpectedOrderHash ||
        path.RootId != UINT32_C(38719404) || path.Tail != 973)
    {
        std::fprintf(stderr,
            "path golden mismatch: path=%016llx order=%016llx root=%u tail=%u\n",
            static_cast<unsigned long long>(path_hash),
            static_cast<unsigned long long>(expected_order_hash),
            path.RootId, path.Tail);
        return 1;
    }

    const std::vector<uint8_t> message = MakeMessage();
    WirehairCodec encoder = wirehair_encoder_create_owned(
        nullptr, message.data(), message.size(), kBlockBytes);
    if (!encoder) {
        std::fprintf(stderr, "encoder creation failed\n");
        return 1;
    }

    std::vector<uint8_t> packets(
        static_cast<size_t>(kBlockCount) * kBlockBytes);
    uint64_t packet_hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < path.EdgeIds.size() + 1; ++i)
    {
        const uint32_t id = i < path.EdgeIds.size()
            ? path.EdgeIds[i]
            : path.RootId;
        uint32_t written = 0;
        uint8_t* packet = packets.data() + i * kBlockBytes;
        if (wirehair_encode(
                encoder, id, packet, kBlockBytes, &written) !=
                Wirehair_Success ||
            written != kBlockBytes)
        {
            std::fprintf(stderr, "encode failed at packet %zu\n", i);
            wirehair_free(encoder);
            return 1;
        }
        packet_hash = FnvU32(packet_hash, id);
        for (uint32_t j = 0; j < written; ++j) {
            packet_hash = FnvByte(packet_hash, packet[j]);
        }
    }

    std::vector<uint8_t> extra_packet(kBlockBytes);
    uint32_t extra_written = 0;
    if (wirehair_encode(
            encoder, path.RootId + 1, extra_packet.data(), kBlockBytes,
            &extra_written) != Wirehair_Success ||
        extra_written != kBlockBytes)
    {
        std::fprintf(stderr, "extra packet encode failed\n");
        wirehair_free(encoder);
        return 1;
    }
    wirehair_free(encoder);

    static const uint64_t kExpectedPacketHash =
        UINT64_C(0x697d28ea99bfc679);
    if (packet_hash != kExpectedPacketHash)
    {
        std::fprintf(stderr, "packet golden mismatch: %016llx\n",
            static_cast<unsigned long long>(packet_hash));
        return 1;
    }

    TestContext context;
    context.Path = &path;
    context.Message = &message;
    context.Packets = &packets;
    context.ExtraPacket = &extra_packet;
    context.ExpectedOrderHash = expected_order_hash;

    int result = 0;
#if defined(__unix__) || defined(__APPLE__)
    result = small_stack ? RunOnSmallStack(context) : RunAvalanche(context);
#else
    if (small_stack) {
        std::printf("small-stack pthread check unavailable on this platform\n");
    }
    result = RunAvalanche(context);
#endif

    if (result == 0)
    {
        std::printf(
            "peel_avalanche_ok: depth=%u stack_bytes=%zu path=%016llx "
            "order=%016llx packets=%016llx\n",
            kBlockCount,
            small_stack ? context.ObservedStackBytes : size_t(0),
            static_cast<unsigned long long>(path_hash),
            static_cast<unsigned long long>(expected_order_hash),
            static_cast<unsigned long long>(packet_hash));
    }
    return result;
}
