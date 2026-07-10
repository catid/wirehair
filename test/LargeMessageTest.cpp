#include <wirehair/wirehair.h>

#include "../WirehairEnvironment.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace {

const uint64_t kMiB = UINT64_C(1024) * 1024;

uint64_t Mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

bool ReadEnvU64(const char* name, uint64_t& out)
{
    const wirehair::EnvironmentValue environment(name);
    const char* value = environment.Get();
    if (!value) {
        return true;
    }
    if (!value[0] || value[0] < '0' || value[0] > '9') {
        std::cerr << "Invalid " << name << " = " << value << std::endl;
        return false;
    }

    errno = 0;
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(value, &end, 0);
    if (errno != 0 || end == value || !end || *end != '\0' ||
        parsed > std::numeric_limits<uint64_t>::max())
    {
        std::cerr << "Invalid " << name << " = " << value << std::endl;
        return false;
    }

    out = static_cast<uint64_t>(parsed);
    return true;
}

struct Config
{
    uint32_t block_count;
    uint32_t block_bytes;
    uint32_t final_bytes;
    uint32_t loss_count;
    uint32_t max_repair;
    uint64_t seed;
    uint64_t max_payload_bytes;
    uint64_t max_milliseconds;
    bool inject_terminal;
    bool inject_timeout;
};

bool LoadConfig(Config& config)
{
    uint64_t block_count = 256;
    uint64_t block_bytes = 257;
    if (!ReadEnvU64("WIREHAIR_LARGE_N", block_count) ||
        !ReadEnvU64("WIREHAIR_LARGE_BLOCK_BYTES", block_bytes))
    {
        return false;
    }

    uint64_t final_bytes = std::min<uint64_t>(block_bytes, 113);
    uint64_t loss_count = std::min<uint64_t>(5, block_count > 0 ? block_count - 1 : 0);
    uint64_t max_repair = 64;
    uint64_t seed = UINT64_C(0x8af09d31e7c4526b);
    uint64_t max_payload_bytes = 64 * kMiB;
    uint64_t max_milliseconds = 60000;
    uint64_t inject_terminal = 0;
    uint64_t inject_timeout = 0;

    if (!ReadEnvU64("WIREHAIR_LARGE_FINAL_BYTES", final_bytes) ||
        !ReadEnvU64("WIREHAIR_LARGE_LOSS_COUNT", loss_count) ||
        !ReadEnvU64("WIREHAIR_LARGE_MAX_REPAIR", max_repair) ||
        !ReadEnvU64("WIREHAIR_LARGE_SEED", seed) ||
        !ReadEnvU64("WIREHAIR_LARGE_MAX_PAYLOAD_BYTES", max_payload_bytes) ||
        !ReadEnvU64("WIREHAIR_LARGE_MAX_MILLISECONDS", max_milliseconds) ||
        !ReadEnvU64("WIREHAIR_LARGE_INJECT_TERMINAL", inject_terminal) ||
        !ReadEnvU64("WIREHAIR_LARGE_INJECT_TIMEOUT", inject_timeout))
    {
        return false;
    }

    if (block_count < 2 || block_count > 64000 || block_bytes == 0 ||
        block_bytes > 0x7fffffffu || final_bytes == 0 ||
        final_bytes > block_bytes || loss_count == 0 ||
        loss_count >= block_count || max_repair == 0 ||
        max_repair > std::numeric_limits<uint32_t>::max() ||
        max_repair > std::numeric_limits<uint32_t>::max() - block_count ||
        max_payload_bytes == 0 || max_milliseconds == 0 ||
        inject_terminal > 1 || inject_timeout > 1)
    {
        std::cerr << "Invalid large-message test parameters" << std::endl;
        return false;
    }

    if (block_count > std::numeric_limits<uint64_t>::max() / block_bytes) {
        std::cerr << "N * block_bytes overflows" << std::endl;
        return false;
    }
    const uint64_t payload_span = block_count * block_bytes;
    if (payload_span > max_payload_bytes) {
        std::cerr << "Resource policy rejected N * block_bytes = "
                  << payload_span << " > " << max_payload_bytes << std::endl;
        return false;
    }

    config.block_count = static_cast<uint32_t>(block_count);
    config.block_bytes = static_cast<uint32_t>(block_bytes);
    config.final_bytes = static_cast<uint32_t>(final_bytes);
    config.loss_count = static_cast<uint32_t>(loss_count);
    config.max_repair = static_cast<uint32_t>(max_repair);
    config.seed = seed;
    config.max_payload_bytes = max_payload_bytes;
    config.max_milliseconds = max_milliseconds;
    config.inject_terminal = inject_terminal != 0;
    config.inject_timeout = inject_timeout != 0;
    return true;
}

typedef std::chrono::steady_clock Clock;

uint64_t ElapsedMilliseconds(const Clock::time_point& start)
{
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - start).count());
}

uint64_t PeakRssKiB()
{
#if defined(__unix__) || defined(__APPLE__)
    struct rusage usage = {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0;
    }
#if defined(__APPLE__)
    return static_cast<uint64_t>(usage.ru_maxrss) / 1024;
#else
    return static_cast<uint64_t>(usage.ru_maxrss);
#endif
#else
    return 0;
#endif
}

class MetricsReporter
{
public:
    explicit MetricsReporter(const Clock::time_point& start)
        : Start(start)
    {
    }

    ~MetricsReporter()
    {
        std::cout << "large_message_metrics: elapsed_ms="
                  << ElapsedMilliseconds(Start)
                  << " peak_rss_kib=" << PeakRssKiB() << std::endl;
    }

private:
    Clock::time_point Start;
};

bool CheckDeadline(
    const Clock::time_point& start,
    uint64_t max_milliseconds,
    const char* phase)
{
    const uint64_t elapsed = ElapsedMilliseconds(start);
    if (elapsed <= max_milliseconds) {
        return true;
    }
    std::cerr << "Large-message timeout in " << phase
              << " after " << elapsed << " ms (cap "
              << max_milliseconds << " ms)" << std::endl;
    return false;
}

uint32_t OriginalBytes(const Config& config, uint32_t block_id)
{
    return block_id + 1 == config.block_count
        ? config.final_bytes
        : config.block_bytes;
}

void FillBlock(uint32_t block_id, uint8_t* block, uint32_t bytes)
{
    uint64_t state = Mix64(UINT64_C(0x99f0b41d5b6a71e3) ^ block_id);

    uint32_t offset = 0;
    while (offset + sizeof(uint64_t) <= bytes) {
        state = Mix64(state + UINT64_C(0x9e3779b97f4a7c15));
        std::memcpy(block + offset, &state, sizeof(state));
        offset += sizeof(state);
    }

    if (offset < bytes) {
        state = Mix64(state + UINT64_C(0x9e3779b97f4a7c15));
        std::memcpy(block + offset, &state, bytes - offset);
    }
}

bool FillMessage(
    uint8_t* message,
    const Config& config,
    const Clock::time_point& start)
{
    size_t offset = 0;
    for (uint32_t block_id = 0; block_id < config.block_count; ++block_id) {
        const uint32_t bytes = OriginalBytes(config, block_id);
        FillBlock(block_id, message + offset, bytes);
        offset += bytes;

        if ((block_id & 1023u) == 1023u &&
            !CheckDeadline(start, config.max_milliseconds, "message fill"))
        {
            return false;
        }
    }
    return CheckDeadline(start, config.max_milliseconds, "message fill");
}

bool VerifyMessage(
    const uint8_t* recovered,
    const Config& config,
    const Clock::time_point& start)
{
    std::unique_ptr<uint8_t[]> expected(
        new (std::nothrow) uint8_t[config.block_bytes]);
    if (!expected) {
        std::cerr << "Failed to allocate verification block" << std::endl;
        return false;
    }

    size_t offset = 0;
    for (uint32_t block_id = 0; block_id < config.block_count; ++block_id) {
        const uint32_t bytes = OriginalBytes(config, block_id);
        FillBlock(block_id, expected.get(), bytes);
        if (std::memcmp(recovered + offset, expected.get(), bytes) != 0) {
            std::cerr << "Recovered block mismatch at block " << block_id
                      << std::endl;
            return false;
        }
        offset += bytes;

        if ((block_id & 1023u) == 1023u &&
            !CheckDeadline(start, config.max_milliseconds, "verification"))
        {
            return false;
        }
    }
    return CheckDeadline(start, config.max_milliseconds, "verification");
}

std::vector<uint32_t> BuildPacketOrder(const Config& config, uint64_t& digest)
{
    std::vector<uint32_t> order(config.block_count);
    for (uint32_t i = 0; i < config.block_count; ++i) {
        order[i] = i;
    }

    uint64_t state = config.seed;
    for (uint32_t i = config.block_count - 1; i > 0; --i) {
        state = Mix64(state + UINT64_C(0x9e3779b97f4a7c15));
        const uint32_t j = static_cast<uint32_t>(state % (i + 1u));
        std::swap(order[i], order[j]);
    }

    digest = UINT64_C(0x6035f12b3a9874cd);
    for (size_t i = 0; i < order.size(); ++i) {
        digest = Mix64(digest ^ (static_cast<uint64_t>(order[i]) << 1) ^ i);
    }
    return order;
}

class CodecGuard
{
public:
    explicit CodecGuard(WirehairCodec codec = nullptr)
        : Codec(codec)
    {
    }

    ~CodecGuard()
    {
        if (Codec) {
            wirehair_free(Codec);
        }
    }

    WirehairCodec get() const
    {
        return Codec;
    }

private:
    CodecGuard(const CodecGuard&);
    CodecGuard& operator=(const CodecGuard&);
    WirehairCodec Codec;
};

bool IsDecodeResult(WirehairResult result)
{
    return result == Wirehair_NeedMore || result == Wirehair_Success;
}

} // namespace

int main()
{
    Config config;
    if (!LoadConfig(config)) {
        return 2;
    }

    const uint64_t payload_span =
        static_cast<uint64_t>(config.block_count) * config.block_bytes;
    const uint64_t message_bytes =
        static_cast<uint64_t>(config.block_count - 1) * config.block_bytes +
        config.final_bytes;
    if (message_bytes >
        static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
    {
        std::cerr << "Message size overflows this platform" << std::endl;
        return 2;
    }

    const Clock::time_point start = Clock::now();
    MetricsReporter metrics(start);

    std::cout << "large_message_test: N=" << config.block_count
              << " block_bytes=" << config.block_bytes
              << " final_bytes=" << config.final_bytes
              << " message_bytes=" << message_bytes
              << " payload_span=" << payload_span
              << " payload_cap=" << config.max_payload_bytes
              << " loss_count=" << config.loss_count
              << " max_repair=" << config.max_repair
              << " seed=" << config.seed
              << " max_ms=" << config.max_milliseconds << std::endl;

    if (config.inject_timeout) {
        std::cerr << "Injected timeout failure" << std::endl;
        return 3;
    }

    const WirehairResult init_result = wirehair_init();
    if (init_result != Wirehair_Success) {
        std::cerr << "wirehair_init failed: "
                  << wirehair_result_string(init_result) << std::endl;
        return 1;
    }

    std::unique_ptr<uint8_t[]> message(
        new (std::nothrow) uint8_t[static_cast<size_t>(message_bytes)]);
    std::unique_ptr<uint8_t[]> recovered(
        new (std::nothrow) uint8_t[static_cast<size_t>(message_bytes)]);
    std::unique_ptr<uint8_t[]> packet(
        new (std::nothrow) uint8_t[config.block_bytes]);
    if (!message || !recovered || !packet) {
        std::cerr << "Failed to allocate large-message buffers" << std::endl;
        return 1;
    }
    std::memset(recovered.get(), 0xa5, static_cast<size_t>(message_bytes));

    if (!FillMessage(message.get(), config, start)) {
        return 3;
    }
    const uint64_t fill_ms = ElapsedMilliseconds(start);

    CodecGuard encoder(wirehair_encoder_create(
        nullptr, message.get(), message_bytes, config.block_bytes));
    if (!encoder.get()) {
        std::cerr << "wirehair_encoder_create failed" << std::endl;
        return 1;
    }

    CodecGuard decoder(wirehair_decoder_create(
        nullptr, message_bytes, config.block_bytes));
    if (!decoder.get()) {
        std::cerr << "wirehair_decoder_create failed" << std::endl;
        return 1;
    }

    uint64_t order_digest = 0;
    const std::vector<uint32_t> order = BuildPacketOrder(config, order_digest);
    std::cout << "packet_order_digest=" << order_digest << " dropped_ids=";
    const uint32_t printed_losses = std::min<uint32_t>(config.loss_count, 32);
    for (uint32_t i = 0; i < printed_losses; ++i) {
        if (i != 0) {
            std::cout << ',';
        }
        std::cout << order[i];
    }
    if (printed_losses != config.loss_count) {
        std::cout << ",...(+" << (config.loss_count - printed_losses) << ')';
    }
    std::cout << std::endl;

    WirehairResult decode_result = Wirehair_NeedMore;
    uint32_t received_originals = 0;
    for (uint32_t position = 0; position < config.block_count; ++position) {
        const uint32_t block_id = order[position];
        uint32_t bytes_out = 0;
        const WirehairResult encode_result = wirehair_encode(
            encoder.get(), block_id, packet.get(), config.block_bytes, &bytes_out);
        if (encode_result != Wirehair_Success) {
            std::cerr << "wirehair_encode original failed at id " << block_id
                      << ": " << wirehair_result_string(encode_result)
                      << std::endl;
            return 1;
        }
        const uint32_t expected_bytes = OriginalBytes(config, block_id);
        if (bytes_out != expected_bytes) {
            std::cerr << "Original packet length mismatch at id " << block_id
                      << ": " << bytes_out << " != " << expected_bytes
                      << std::endl;
            return 1;
        }

        if (position < config.loss_count) {
            continue;
        }

        const uint32_t decode_bytes =
            config.inject_terminal && received_originals == 0
                ? 0
                : bytes_out;
        decode_result = wirehair_decode(
            decoder.get(), block_id, packet.get(), decode_bytes);
        if (config.inject_terminal && received_originals == 0) {
            if (IsDecodeResult(decode_result)) {
                std::cerr << "Terminal decode injection was not rejected"
                          << std::endl;
                return 2;
            }
            std::cerr << "Injected terminal decode result: "
                      << wirehair_result_string(decode_result) << std::endl;
            return 1;
        }
        if (!IsDecodeResult(decode_result)) {
            std::cerr << "wirehair_decode original failed at id " << block_id
                      << ": " << wirehair_result_string(decode_result)
                      << std::endl;
            return 1;
        }
        ++received_originals;

        if ((position & 1023u) == 1023u &&
            !CheckDeadline(start, config.max_milliseconds, "original decode"))
        {
            return 3;
        }
    }

    uint32_t repair_count = 0;
    while (decode_result != Wirehair_Success &&
           repair_count < config.max_repair)
    {
        const uint32_t repair_id = config.block_count + repair_count;
        uint32_t bytes_out = 0;
        const WirehairResult encode_result = wirehair_encode(
            encoder.get(), repair_id, packet.get(), config.block_bytes, &bytes_out);
        if (encode_result != Wirehair_Success || bytes_out == 0 ||
            bytes_out > config.block_bytes)
        {
            std::cerr << "wirehair_encode repair failed at id " << repair_id
                      << ": " << wirehair_result_string(encode_result)
                      << " bytes=" << bytes_out << std::endl;
            return 1;
        }

        decode_result = wirehair_decode(
            decoder.get(), repair_id, packet.get(), bytes_out);
        ++repair_count;
        if (!IsDecodeResult(decode_result)) {
            std::cerr << "wirehair_decode repair failed at id " << repair_id
                      << ": " << wirehair_result_string(decode_result)
                      << std::endl;
            return 1;
        }
        if (!CheckDeadline(start, config.max_milliseconds, "repair decode")) {
            return 3;
        }
    }

    if (decode_result != Wirehair_Success) {
        std::cerr << "Decoder still needs more after " << repair_count
                  << " repair blocks (cap reached)" << std::endl;
        return 1;
    }
    const uint64_t decode_ms = ElapsedMilliseconds(start);

    const WirehairResult recover_result =
        wirehair_recover(decoder.get(), recovered.get(), message_bytes);
    if (recover_result != Wirehair_Success) {
        std::cerr << "wirehair_recover failed: "
                  << wirehair_result_string(recover_result) << std::endl;
        return 1;
    }
    const uint64_t recover_ms = ElapsedMilliseconds(start);

    if (!VerifyMessage(recovered.get(), config, start)) {
        return 1;
    }
    const uint64_t verify_ms = ElapsedMilliseconds(start);

    std::cout << "large_message_phases_ms: fill=" << fill_ms
              << " decode=" << (decode_ms - fill_ms)
              << " recover=" << (recover_ms - decode_ms)
              << " verify=" << (verify_ms - recover_ms) << std::endl;
    std::cout << "large_message_test: PASS received_originals="
              << received_originals << " repair_blocks=" << repair_count
              << std::endl;
    return 0;
}
