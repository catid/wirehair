#include "WirehairV2Contract.h"
#include "WirehairV2PrecodeDecode.h"
#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <signal.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#endif

#if defined(__linux__)
#include <sched.h>
#endif

namespace {

static const char kSchema[] = "wirehair.wh2.hook-timing.v1";
static const char* const kEquationEnvironment[] = {
    "WIREHAIR_V2_PEEL_DEGREES",
    "WIREHAIR_V2_STAIRCASE_DEGREES",
    "WIREHAIR_V2_STAIRCASE_ROW_DEGREES",
    "WIREHAIR_V2_STAIRCASE_DEGREE_SCALE",
    "WIREHAIR_V2_BAND_TRACKING_X"
};

const uint32_t kSha256K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
    0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
    0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
    0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
    0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
    0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

inline uint32_t Rotr32(uint32_t value, unsigned bits)
{
    return (value >> bits) | (value << (32u - bits));
}

class Sha256
{
public:
    Sha256()
        : TotalBytes(0u)
        , BufferBytes(0u)
    {
        State[0] = 0x6a09e667u;
        State[1] = 0xbb67ae85u;
        State[2] = 0x3c6ef372u;
        State[3] = 0xa54ff53au;
        State[4] = 0x510e527fu;
        State[5] = 0x9b05688cu;
        State[6] = 0x1f83d9abu;
        State[7] = 0x5be0cd19u;
    }

    void Update(const void* data, size_t bytes)
    {
        const uint8_t* input = static_cast<const uint8_t*>(data);
        TotalBytes += bytes;
        if (BufferBytes != 0u)
        {
            while (bytes != 0u && BufferBytes < sizeof(Buffer)) {
                Buffer[BufferBytes++] = *input++;
                --bytes;
            }
            if (BufferBytes != sizeof(Buffer)) return;
            ProcessBlock(Buffer);
            BufferBytes = 0u;
        }
        while (bytes >= sizeof(Buffer))
        {
            ProcessBlock(input);
            input += sizeof(Buffer);
            bytes -= sizeof(Buffer);
        }
        while (bytes != 0u) {
            Buffer[BufferBytes++] = *input++;
            --bytes;
        }
    }

    void Finalize(uint8_t* digest_out)
    {
        const uint64_t total_bits = TotalBytes * 8u;
        const uint8_t first = 0x80u;
        const uint8_t zero = 0u;
        Update(&first, 1u);
        while (BufferBytes != 56u) Update(&zero, 1u);
        uint8_t length[8];
        for (unsigned i = 0u; i < 8u; ++i) {
            length[i] = (uint8_t)(total_bits >> (56u - 8u * i));
        }
        Update(length, sizeof(length));
        for (unsigned i = 0u; i < 8u; ++i)
        {
            digest_out[4u * i] = (uint8_t)(State[i] >> 24);
            digest_out[4u * i + 1u] = (uint8_t)(State[i] >> 16);
            digest_out[4u * i + 2u] = (uint8_t)(State[i] >> 8);
            digest_out[4u * i + 3u] = (uint8_t)State[i];
        }
    }

private:
    void ProcessBlock(const uint8_t* block)
    {
        uint32_t words[64];
        for (unsigned i = 0u; i < 16u; ++i)
        {
            words[i] = ((uint32_t)block[4u * i] << 24) |
                ((uint32_t)block[4u * i + 1u] << 16) |
                ((uint32_t)block[4u * i + 2u] << 8) |
                (uint32_t)block[4u * i + 3u];
        }
        for (unsigned i = 16u; i < 64u; ++i)
        {
            const uint32_t s0 = Rotr32(words[i - 15u], 7u) ^
                Rotr32(words[i - 15u], 18u) ^ (words[i - 15u] >> 3);
            const uint32_t s1 = Rotr32(words[i - 2u], 17u) ^
                Rotr32(words[i - 2u], 19u) ^ (words[i - 2u] >> 10);
            words[i] =
                words[i - 16u] + s0 + words[i - 7u] + s1;
        }
        uint32_t a = State[0], b = State[1], c = State[2], d = State[3];
        uint32_t e = State[4], f = State[5], g = State[6], h = State[7];
        for (unsigned i = 0u; i < 64u; ++i)
        {
            const uint32_t s1 =
                Rotr32(e, 6u) ^ Rotr32(e, 11u) ^ Rotr32(e, 25u);
            const uint32_t ch = (e & f) ^ (~e & g);
            const uint32_t temp1 =
                h + s1 + ch + kSha256K[i] + words[i];
            const uint32_t s0 =
                Rotr32(a, 2u) ^ Rotr32(a, 13u) ^ Rotr32(a, 22u);
            const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t temp2 = s0 + maj;
            h = g; g = f; f = e; e = d + temp1;
            d = c; c = b; b = a; a = temp1 + temp2;
        }
        State[0] += a; State[1] += b; State[2] += c; State[3] += d;
        State[4] += e; State[5] += f; State[6] += g; State[7] += h;
    }

    uint32_t State[8];
    uint64_t TotalBytes;
    uint8_t Buffer[64];
    size_t BufferBytes;
};

class Hasher
{
public:
    explicit Hasher(const char* tag)
    {
        String(tag);
    }

    void U8(uint8_t value) { Digest.Update(&value, sizeof(value)); }
    void U16(uint16_t value)
    {
        const uint8_t bytes[2] = {
            (uint8_t)value, (uint8_t)(value >> 8)
        };
        Digest.Update(bytes, sizeof(bytes));
    }
    void U32(uint32_t value)
    {
        const uint8_t bytes[4] = {
            (uint8_t)value, (uint8_t)(value >> 8),
            (uint8_t)(value >> 16), (uint8_t)(value >> 24)
        };
        Digest.Update(bytes, sizeof(bytes));
    }
    void U64(uint64_t value)
    {
        U32((uint32_t)value);
        U32((uint32_t)(value >> 32));
    }
    void Double(double value)
    {
        uint64_t bits = 0u;
        static_assert(sizeof(bits) == sizeof(value), "double bits");
        std::memcpy(&bits, &value, sizeof(bits));
        U64(bits);
    }
    void Bytes(const void* data, size_t bytes)
    {
        U64((uint64_t)bytes);
        if (bytes != 0u) Digest.Update(data, bytes);
    }
    void String(const std::string& text)
    {
        Bytes(text.data(), text.size());
    }
    void String(const char* text)
    {
        String(std::string(text ? text : ""));
    }
    std::string Finish()
    {
        uint8_t digest[32];
        Digest.Finalize(digest);
        static const char hex[] = "0123456789abcdef";
        std::string output(64u, '0');
        for (size_t i = 0u; i < sizeof(digest); ++i) {
            output[2u * i] = hex[digest[i] >> 4];
            output[2u * i + 1u] = hex[digest[i] & 15u];
        }
        return output;
    }

private:
    Sha256 Digest;
};

std::string Sha256Bytes(const void* data, size_t bytes, const char* tag)
{
    Hasher hasher(tag);
    hasher.Bytes(data, bytes);
    return hasher.Finish();
}

bool CheckSha256()
{
    Sha256 sha;
    static const char abc[] = "abc";
    sha.Update(abc, 3u);
    uint8_t digest[32];
    sha.Finalize(digest);
    static const uint8_t expected[32] = {
        0xba, 0x78, 0x16, 0xbf, 0x8f, 0x01, 0xcf, 0xea,
        0x41, 0x41, 0x40, 0xde, 0x5d, 0xae, 0x22, 0x23,
        0xb0, 0x03, 0x61, 0xa3, 0x96, 0x17, 0x7a, 0x9c,
        0xb4, 0x10, 0xff, 0x61, 0xf2, 0x00, 0x15, 0xad
    };
    return std::memcmp(digest, expected, sizeof(expected)) == 0;
}

std::string RawSha256(const std::string& input)
{
    Sha256 sha;
    sha.Update(input.data(), input.size());
    uint8_t digest[32];
    sha.Finalize(digest);
    static const char hex[] = "0123456789abcdef";
    std::string output(64u, '0');
    for (size_t i = 0u; i < sizeof(digest); ++i) {
        output[2u * i] = hex[digest[i] >> 4];
        output[2u * i + 1u] = hex[digest[i] & 15u];
    }
    return output;
}

enum class Schedule
{
    Iid,
    Burst,
    Permutation,
    SystematicFirst,
    RepairOnly,
    Adversarial
};

enum ScopeBits
{
    ScopeSemantic = 0,
    ScopeRow = 1,
    ScopeEncoder = 2,
    ScopeDecoderFeed = 4,
    ScopeDecoderFull = 8,
    ScopeDirect = 16,
    ScopeAll = 31
};

struct Options
{
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    uint64_t ConstructionSeed = 0u;
    uint64_t LossSeed = 0u;
    uint32_t LossPpm = 0u;
    Schedule PacketSchedule = Schedule::Iid;
    uint32_t Scopes = 0u;
    uint32_t WarmupReps = 0u;
    uint32_t InnerReps = 0u;
    uint64_t MaxWorkingMiB = 0u;
    std::string ContextSha256;
    int ReadyFd = -1;
    int GoFd = -1;
};

const char* ScheduleName(Schedule schedule)
{
    switch (schedule)
    {
    case Schedule::Iid: return "iid";
    case Schedule::Burst: return "burst";
    case Schedule::Permutation: return "permutation";
    case Schedule::SystematicFirst: return "systematic-first";
    case Schedule::RepairOnly: return "repair-only";
    case Schedule::Adversarial: return "adversarial";
    }
    return "invalid";
}

bool ParseCanonicalU64(const char* text, uint64_t& value)
{
    if (!text || !*text || (text[0] == '0' && text[1] != '\0')) {
        return false;
    }
    uint64_t parsed = 0u;
    for (const char* p = text; *p; ++p)
    {
        if (*p < '0' || *p > '9') return false;
        const uint32_t digit = (uint32_t)(*p - '0');
        if (parsed > (UINT64_MAX - digit) / 10u) return false;
        parsed = parsed * 10u + digit;
    }
    value = parsed;
    return true;
}

bool ParseCanonicalU32(const char* text, uint32_t& value)
{
    uint64_t parsed = 0u;
    if (!ParseCanonicalU64(text, parsed) || parsed > UINT32_MAX) return false;
    value = (uint32_t)parsed;
    return true;
}

bool ParseCanonicalFd(const char* text, int& value)
{
    uint64_t parsed = 0u;
    if (!ParseCanonicalU64(text, parsed) ||
        parsed < 3u || parsed > (uint64_t)INT_MAX)
    {
        return false;
    }
    value = (int)parsed;
    return true;
}

bool IsLowerHexSha256(const char* text)
{
    if (!text || std::strlen(text) != 64u) return false;
    for (size_t i = 0u; i < 64u; ++i) {
        if (!((text[i] >= '0' && text[i] <= '9') ||
              (text[i] >= 'a' && text[i] <= 'f'))) return false;
    }
    return true;
}

bool ParseSchedule(const char* text, Schedule& schedule)
{
    if (!std::strcmp(text, "iid")) schedule = Schedule::Iid;
    else if (!std::strcmp(text, "burst")) schedule = Schedule::Burst;
    else if (!std::strcmp(text, "permutation")) schedule = Schedule::Permutation;
    else if (!std::strcmp(text, "systematic-first"))
        schedule = Schedule::SystematicFirst;
    else if (!std::strcmp(text, "repair-only"))
        schedule = Schedule::RepairOnly;
    else if (!std::strcmp(text, "adversarial"))
        schedule = Schedule::Adversarial;
    else return false;
    return true;
}

bool ParseScope(const char* text, uint32_t& scopes)
{
    if (!std::strcmp(text, "semantic")) scopes = ScopeSemantic;
    else if (!std::strcmp(text, "row")) scopes = ScopeRow;
    else if (!std::strcmp(text, "encoder")) scopes = ScopeEncoder;
    else if (!std::strcmp(text, "decoder-feed")) scopes = ScopeDecoderFeed;
    else if (!std::strcmp(text, "decoder-full")) scopes = ScopeDecoderFull;
    else if (!std::strcmp(text, "direct")) scopes = ScopeDirect;
    else if (!std::strcmp(text, "all")) scopes = ScopeAll;
    else return false;
    return true;
}

const char* ScopeRequestName(uint32_t scopes)
{
    switch (scopes)
    {
    case ScopeSemantic: return "semantic";
    case ScopeRow: return "row";
    case ScopeEncoder: return "encoder";
    case ScopeDecoderFeed: return "decoder-feed";
    case ScopeDecoderFull: return "decoder-full";
    case ScopeDirect: return "direct";
    case ScopeAll: return "all";
    }
    return "invalid";
}

const char* TimingLifecycle(const char* scope)
{
    if (!std::strcmp(scope, "row")) {
        return "caller-row-buffer-reuse-v1";
    }
    if (!std::strcmp(scope, "encoder")) {
        return "fresh-first-then-transactional-reinitialize-"
            "including-prior-release-v1";
    }
    if (!std::strcmp(scope, "decoder-feed")) {
        return "distinct-preinitialized-endpoints-v1";
    }
    if (!std::strcmp(scope, "decoder-full")) {
        return "fresh-first-then-transactional-reinitialize-"
            "including-prior-release-v1";
    }
    if (!std::strcmp(scope, "direct")) {
        return "fresh-first-then-transactional-output-reuse-v1";
    }
    return "invalid";
}

void Usage(const char* program)
{
    std::fprintf(stderr,
        "usage: %s --K N --bb BYTES --construction-seed N --loss-seed N "
        "--loss-ppm N --schedule iid|burst|permutation|systematic-first|"
        "repair-only|adversarial --scope semantic|row|encoder|decoder-feed|"
        "decoder-full|direct|all --warmup-reps N --inner-reps N "
        "--max-working-mib N --context-sha256 HEX "
        "[--ready-fd N --go-fd N]\n", program);
}

bool ParseOptions(int argc, char** argv, Options& options)
{
    enum SeenBits
    {
        SeenK = 1 << 0, SeenBb = 1 << 1, SeenConstruction = 1 << 2,
        SeenLossSeed = 1 << 3, SeenLossPpm = 1 << 4,
        SeenSchedule = 1 << 5, SeenScope = 1 << 6,
        SeenWarmup = 1 << 7, SeenInner = 1 << 8,
        SeenMemory = 1 << 9, SeenContext = 1 << 10,
        SeenReadyFd = 1 << 11, SeenGoFd = 1 << 12,
        SeenRequired = (1 << 11) - 1
    };
    uint32_t seen = 0u;
    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--help")) return false;
        if (i + 1 >= argc) return false;
        const char* option = argv[i];
        const char* value = argv[++i];
        uint32_t bit = 0u;
        bool ok = false;
        if (!std::strcmp(option, "--K")) {
            bit = SeenK; ok = ParseCanonicalU32(value, options.K);
        } else if (!std::strcmp(option, "--bb")) {
            bit = SeenBb; ok = ParseCanonicalU32(value, options.BlockBytes);
        } else if (!std::strcmp(option, "--construction-seed")) {
            bit = SeenConstruction;
            ok = ParseCanonicalU64(value, options.ConstructionSeed);
        } else if (!std::strcmp(option, "--loss-seed")) {
            bit = SeenLossSeed; ok = ParseCanonicalU64(value, options.LossSeed);
        } else if (!std::strcmp(option, "--loss-ppm")) {
            bit = SeenLossPpm; ok = ParseCanonicalU32(value, options.LossPpm);
        } else if (!std::strcmp(option, "--schedule")) {
            bit = SeenSchedule; ok = ParseSchedule(value, options.PacketSchedule);
        } else if (!std::strcmp(option, "--scope")) {
            bit = SeenScope; ok = ParseScope(value, options.Scopes);
        } else if (!std::strcmp(option, "--warmup-reps")) {
            bit = SeenWarmup; ok = ParseCanonicalU32(value, options.WarmupReps);
        } else if (!std::strcmp(option, "--inner-reps")) {
            bit = SeenInner; ok = ParseCanonicalU32(value, options.InnerReps);
        } else if (!std::strcmp(option, "--max-working-mib")) {
            bit = SeenMemory; ok = ParseCanonicalU64(value, options.MaxWorkingMiB);
        } else if (!std::strcmp(option, "--context-sha256")) {
            bit = SeenContext; ok = IsLowerHexSha256(value);
            if (ok) options.ContextSha256 = value;
        } else if (!std::strcmp(option, "--ready-fd")) {
            bit = SeenReadyFd;
            ok = ParseCanonicalFd(value, options.ReadyFd);
        } else if (!std::strcmp(option, "--go-fd")) {
            bit = SeenGoFd;
            ok = ParseCanonicalFd(value, options.GoFd);
        } else {
            return false;
        }
        if (!ok || (seen & bit) != 0u) return false;
        seen |= bit;
    }
    const bool ready_seen = (seen & SeenReadyFd) != 0u;
    const bool go_seen = (seen & SeenGoFd) != 0u;
    const bool barrier_enabled = ready_seen && go_seen;
    const bool one_timed_scope =
        options.Scopes != ScopeSemantic &&
        options.Scopes != ScopeAll &&
        (options.Scopes & (options.Scopes - 1u)) == 0u;
    return (seen & SeenRequired) == SeenRequired &&
        ready_seen == go_seen &&
        (!barrier_enabled ||
            (options.ReadyFd != options.GoFd &&
             options.WarmupReps == 0u &&
             one_timed_scope)) &&
        options.K >= 2u && options.K <= 64000u &&
        options.BlockBytes > 0u &&
        options.BlockBytes <= UINT32_C(0x7fffffff) &&
        (options.BlockBytes & 1u) == 0u &&
        options.LossPpm < 1000000u &&
        options.WarmupReps <= 1000u &&
        options.InnerReps > 0u && options.InnerReps <= 1000000u &&
        options.MaxWorkingMiB > 0u &&
        options.MaxWorkingMiB <= UINT64_C(1048576);
}

bool RejectEquationEnvironment()
{
    for (size_t i = 0u;
         i < sizeof(kEquationEnvironment) / sizeof(kEquationEnvironment[0]);
         ++i)
    {
        if (std::getenv(kEquationEnvironment[i]) != nullptr) {
            std::fprintf(stderr, "refusing ambient equation variable: %s\n",
                kEquationEnvironment[i]);
            return false;
        }
    }
    return true;
}

class StartBarrier
{
public:
    StartBarrier(int ready_fd, int go_fd)
        : ReadyFd(ready_fd)
        , GoFd(go_fd)
    {
    }

    ~StartBarrier()
    {
        CloseIgnoringErrors(ReadyFd);
        CloseIgnoringErrors(GoFd);
    }

    StartBarrier(const StartBarrier&) = delete;
    StartBarrier& operator=(const StartBarrier&) = delete;

    bool Enabled() const
    {
        return ReadyFd >= 0;
    }

    bool Consumed() const
    {
        return WasConsumed;
    }

    bool Validate() const
    {
        if (!Enabled()) return GoFd < 0;
#if defined(__unix__) || defined(__APPLE__)
        return ReadyFd >= 3 && GoFd >= 3 && ReadyFd != GoFd &&
            ValidatePipeEnd(ReadyFd, O_WRONLY) &&
            ValidatePipeEnd(GoFd, O_RDONLY);
#else
        return false;
#endif
    }

    bool Rendezvous()
    {
        if (!Enabled() || WasConsumed || !Validate()) return false;
        WasConsumed = true;
#if defined(__unix__) || defined(__APPLE__)
        const bool ready_ok = WriteReadyAndClose();
        if (!ready_ok)
        {
            CloseIgnoringErrors(GoFd);
            return false;
        }
        return ReadGoAndClose();
#else
        return false;
#endif
    }

private:
#if defined(__unix__) || defined(__APPLE__)
    static bool ValidatePipeEnd(int fd, int access_mode)
    {
        struct stat status = {};
        const int descriptor_flags = fcntl(fd, F_GETFD);
        const int file_flags = fcntl(fd, F_GETFL);
        if (descriptor_flags < 0 || file_flags < 0 ||
            fstat(fd, &status) != 0 ||
            !S_ISFIFO(status.st_mode) ||
            (file_flags & O_ACCMODE) != access_mode ||
            (file_flags & O_NONBLOCK) != 0)
        {
            return false;
        }
#if defined(O_PATH)
        if ((file_flags & O_PATH) != 0) return false;
#endif
        return true;
    }

    static ssize_t ReadRetry(int fd, void* data, size_t bytes)
    {
        for (;;)
        {
            const ssize_t result = read(fd, data, bytes);
            if (result >= 0 || errno != EINTR) return result;
        }
    }

    static ssize_t WriteRetry(int fd, const void* data, size_t bytes)
    {
        for (;;)
        {
            const ssize_t result = write(fd, data, bytes);
            if (result >= 0 || errno != EINTR) return result;
        }
    }

    bool WriteReadyAndClose()
    {
        struct sigaction ignored = {};
        struct sigaction previous = {};
        ignored.sa_handler = SIG_IGN;
        if (sigemptyset(&ignored.sa_mask) != 0 ||
            sigaction(SIGPIPE, &ignored, &previous) != 0)
        {
            CloseIgnoringErrors(ReadyFd);
            return false;
        }
        const char ready = 'R';
        const ssize_t written = WriteRetry(ReadyFd, &ready, 1u);
        const bool restored =
            sigaction(SIGPIPE, &previous, nullptr) == 0;
        const bool closed = CloseChecked(ReadyFd);
        return written == 1 && restored && closed;
    }

    bool ReadGoAndClose()
    {
        char go = '\0';
        const ssize_t first = ReadRetry(GoFd, &go, 1u);
        char extra = '\0';
        const ssize_t second =
            first == 1 && go == 'G' ?
            ReadRetry(GoFd, &extra, 1u) : -1;
        const bool closed = CloseChecked(GoFd);
        return first == 1 && go == 'G' && second == 0 && closed;
    }

    static bool CloseChecked(int& fd)
    {
        if (fd < 0) return false;
        const int closing = fd;
        fd = -1;
        return close(closing) == 0;
    }

    static void CloseIgnoringErrors(int& fd)
    {
        if (fd < 0) return;
        const int closing = fd;
        fd = -1;
        (void)close(closing);
    }
#else
    static void CloseIgnoringErrors(int& fd)
    {
        fd = -1;
    }
#endif

    int ReadyFd = -1;
    int GoFd = -1;
    bool WasConsumed = false;
};

const char* StartBarrierName(const Options& options)
{
    return options.ReadyFd >= 0 ? "ready-go-pipe-v1" : "none";
}

bool PinToFirstAllowedCpu(int& cpu_out)
{
#if defined(__linux__)
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0) return false;
    int selected = -1;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        if (CPU_ISSET(cpu, &allowed)) { selected = cpu; break; }
    }
    if (selected < 0) return false;
    cpu_set_t single;
    CPU_ZERO(&single);
    CPU_SET(selected, &single);
    if (sched_setaffinity(0, sizeof(single), &single) != 0) return false;
    if (sched_getcpu() != selected) return false;
    cpu_out = selected;
    return true;
#else
    (void)cpu_out;
    return false;
#endif
}

struct Rng
{
    uint64_t State;

    explicit Rng(uint64_t seed) : State(seed) {}

    uint64_t Next()
    {
        uint64_t value =
            (State += UINT64_C(0x9e3779b97f4a7c15));
        value = (value ^ (value >> 30)) *
            UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27)) *
            UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31);
    }

    uint32_t U32() { return (uint32_t)(Next() >> 32); }
};

uint64_t DropThreshold(uint32_t loss_ppm)
{
    return ((uint64_t)loss_ppm << 32) / UINT64_C(1000000);
}

bool ShouldDrop(Rng& rng, uint32_t loss_ppm)
{
    return (uint64_t)rng.U32() < DropThreshold(loss_ppm);
}

bool BuildSchedule(
    uint32_t source_count,
    uint32_t delivered_count,
    uint32_t loss_ppm,
    uint64_t seed,
    Schedule kind,
    std::vector<uint32_t>& output)
{
    output.clear();
    output.reserve(delivered_count);
    Rng rng(seed ^ UINT64_C(0x5748324c4f535331)); // "WH2LOSS1"
    uint64_t candidate_index = 0u;
    uint32_t burst_remaining = 0u;
    std::vector<uint32_t> permutation;
    size_t permutation_index = 0u;
    uint32_t permutation_base = 0u;

    const auto next_candidate = [&]() -> uint32_t {
        if (kind == Schedule::RepairOnly) {
            return source_count + (uint32_t)candidate_index++;
        }
        if (kind == Schedule::Adversarial) {
            return UINT32_MAX - (uint32_t)(candidate_index++ * 2u);
        }
        if (kind == Schedule::Permutation)
        {
            if (permutation_index >= permutation.size())
            {
                const uint32_t count =
                    std::min<uint32_t>(source_count + 512u, 65536u);
                permutation.resize(count);
                for (uint32_t i = 0u; i < count; ++i) {
                    permutation[i] = permutation_base + i;
                }
                for (uint32_t i = count; i > 1u; --i) {
                    std::swap(
                        permutation[i - 1u],
                        permutation[rng.U32() % i]);
                }
                permutation_base += count;
                permutation_index = 0u;
            }
            return permutation[permutation_index++];
        }
        if (kind == Schedule::SystematicFirst &&
            candidate_index < source_count)
        {
            if (permutation.empty())
            {
                permutation.resize(source_count);
                for (uint32_t i = 0u; i < source_count; ++i) {
                    permutation[i] = i;
                }
                for (uint32_t i = source_count; i > 1u; --i) {
                    std::swap(
                        permutation[i - 1u],
                        permutation[rng.U32() % i]);
                }
            }
            return permutation[(size_t)candidate_index++];
        }
        return (uint32_t)candidate_index++;
    };

    const uint64_t candidate_limit =
        (uint64_t)delivered_count * 256u + 65536u;
    uint64_t candidates = 0u;
    while (output.size() < delivered_count &&
           candidates++ < candidate_limit)
    {
        const uint32_t id = next_candidate();
        bool drop = false;
        if (kind == Schedule::Burst)
        {
            static const uint32_t kBurstLength = 8u;
            if (burst_remaining != 0u) {
                --burst_remaining;
                drop = true;
            }
            else
            {
                const uint64_t denominator =
                    UINT64_C(8000000) - 7u * loss_ppm;
                const uint32_t start_ppm = (uint32_t)(
                    (uint64_t)loss_ppm * UINT64_C(1000000) /
                    denominator);
                if (ShouldDrop(rng, start_ppm)) {
                    burst_remaining = kBurstLength - 1u;
                    drop = true;
                }
            }
        }
        else {
            drop = ShouldDrop(rng, loss_ppm);
        }
        if (!drop) output.push_back(id);
    }
    return output.size() == delivered_count;
}

bool CheckMemoryBudget(const Options& options)
{
    const uint64_t block_factor = (uint64_t)options.K * 8u + 4096u;
    if (block_factor >
        (UINT64_MAX - UINT64_C(16) * 1024u * 1024u) /
            options.BlockBytes)
    {
        return false;
    }
    const uint64_t estimated =
        block_factor * options.BlockBytes +
        (uint64_t)options.K * 64u +
        UINT64_C(16) * 1024u * 1024u;
    if (options.MaxWorkingMiB > (UINT64_MAX >> 20)) return false;
    const uint64_t allowed = options.MaxWorkingMiB << 20;
    if (estimated > allowed)
    {
        std::fprintf(stderr,
            "working-set estimate exceeds --max-working-mib "
            "(estimate=%llu limit=%llu)\n",
            (unsigned long long)estimated,
            (unsigned long long)allowed);
        return false;
    }
    return true;
}

void HashPeelingCodec(
    Hasher& hash,
    const wirehair_v2::PeelingCodec& codec)
{
    hash.U32((uint32_t)codec.Solver);
    hash.U32((uint32_t)codec.Structure);
    hash.U32((uint32_t)codec.Family);
    hash.U16(codec.MinDegree);
    hash.U16(codec.MaxDegree);
    hash.U16(codec.SolverCandidateLimit);
    hash.Double(codec.Degree1Mass);
    hash.Double(codec.Degree2Mass);
    hash.Double(codec.RobustC);
    hash.Double(codec.RobustDelta);
    hash.U8(codec.FullyRandomRows ? 1u : 0u);
    hash.U8(codec.UseWirehairRowDistribution ? 1u : 0u);
}

std::string HashProfile(const wirehair_v2::SeedProfile& profile)
{
    Hasher hash("WH2HTPROFILE1");
    hash.U32(profile.BlockCount);
    hash.U32(profile.BlockBytes);
    hash.U32((uint32_t)profile.Policy.Solver);
    hash.U32((uint32_t)profile.Policy.Structure);
    hash.U32((uint32_t)profile.Policy.ByteClass);
    hash.U32((uint32_t)profile.Policy.CountBand);
    HashPeelingCodec(hash, profile.Policy.Codec);
    hash.U16(profile.DenseCount);
    hash.U16(profile.PeelSeed);
    hash.U16(profile.DenseSeed);
    hash.U16(profile.PeelSeedBucket);
    hash.U8(profile.UsedPeelFixup ? 1u : 0u);
    hash.U8(profile.UsedDenseFixup ? 1u : 0u);
    hash.U8(profile.V2SeedSelected ? 1u : 0u);
    hash.U32(profile.V2SeedAttempt);
    hash.U32(profile.V2PrecodeContractVersion);
    hash.U32(profile.V2PacketRowContractVersion);
    hash.U32((uint32_t)profile.V2Architecture);
    hash.U32((uint32_t)profile.V2SeedPolicy);
    hash.U32(profile.V2StaircaseCount);
    hash.U32(profile.V2DenseRowCount);
    hash.U32(profile.V2HeavyRowCount);
    hash.U32((uint32_t)profile.V2CompletionField);
    hash.U32(profile.V2SourceHits);
    hash.U64(profile.V2PrecodeSeed);
    hash.U32(profile.V2PacketPeelSeed);
    hash.U32(profile.V2RecoveryMixCount);
    hash.U8(profile.V2DenseIdentityCorner ? 1u : 0u);
    hash.U8(profile.V2DenseTwoAnchor ? 1u : 0u);
    hash.U8(profile.V2AdaptiveDenseTwoAnchor ? 1u : 0u);
    hash.U64(profile.V2PrecodeSeedSalt);
    hash.U64(profile.V2RecoveryRowSeedSalt);
    hash.U8(profile.Tuned ? 1u : 0u);
    hash.Double(profile.TuningResidualMean);
    hash.U32(profile.TuningResidualColumns);
    hash.U64(profile.TuningXorCost);
    hash.U16(profile.TuningCandidatesRequested);
    hash.U16(profile.TuningCandidatesUnique);
    hash.U16(profile.TuningCandidatesCompleted);
    hash.U32(profile.TuningTrials);
    return hash.Finish();
}

void HashParams(
    Hasher& hash,
    const wirehair_v2::PrecodeParams& params)
{
    hash.U32(params.BlockCount);
    hash.U32(params.Staircase);
    hash.U32(params.DenseRows);
    hash.U32(params.HeavyRows);
    hash.U32(params.SourceHits);
    hash.U32((uint32_t)params.Field);
    hash.U8(params.DegreeBalancedStaircase ? 1u : 0u);
    hash.Double(params.StaircaseDegreeScale);
    hash.U8(params.DenseIdentityCorner ? 1u : 0u);
    hash.U8(params.DenseTwoAnchor ? 1u : 0u);
    hash.U32(params.DenseTwoAnchorPhase);
    hash.U32((uint32_t)params.SegmentedDenseAnchors);
    hash.U32((uint32_t)params.HeavyFamily);
    hash.U64(params.Seed);
}

std::string HashSystem(const wirehair_v2::PrecodeSystem& system)
{
    Hasher hash("WH2HTSYSTEM1");
    HashParams(hash, system.Params);
    hash.U32((uint32_t)system.StaircaseRows.size());
    for (const std::vector<uint32_t>& row : system.StaircaseRows)
    {
        hash.U32((uint32_t)row.size());
        for (uint32_t column : row) hash.U32(column);
    }
    hash.U32((uint32_t)system.DenseRowColumns.size());
    for (const std::vector<uint32_t>& row : system.DenseRowColumns)
    {
        hash.U32((uint32_t)row.size());
        for (uint32_t column : row) hash.U32(column);
    }
    return hash.Finish();
}

std::string HashParamsOnly(const wirehair_v2::PrecodeParams& params)
{
    Hasher hash("WH2HTPARAMS1");
    HashParams(hash, params);
    return hash.Finish();
}

std::string HashCoefficients()
{
    const wirehair_v2::MixedCoefficientRows* rows =
        wirehair_v2::GetMixedCoefficientRows();
    if (!rows) return std::string();
    const uint32_t gf256_rows = wirehair_v2::ActiveMixedGF256Rows();
    const uint32_t gf16_rows = wirehair_v2::ActiveMixedGF16Rows();
    const uint32_t period = wirehair_v2::ActiveMixedCoefficientPeriod();
    Hasher hash("WH2HTCOEFFICIENT1");
    hash.U32(gf256_rows);
    hash.U32(gf16_rows);
    hash.U32(period);
    hash.U32((uint32_t)wirehair_v2::ActiveMixedCoefficientGeometry());
    hash.U32((uint32_t)wirehair_v2::ActiveMixedResidueSchedule());
    hash.U32(wirehair_v2::ActiveMixedResidueSkew());
    hash.U32(wirehair_v2::ActiveMixedResidueHashSeed());
    hash.U8(wirehair_v2::ActiveMixedResiduesRotated() ? 1u : 0u);
    hash.U8(
        wirehair_v2::ActiveMixedIndependentExtensionResidues() ? 1u : 0u);
    hash.U32(wirehair_v2::ActiveMixedGroupedGF256Rows());
    hash.U32(wirehair_v2::ActiveMixedGroupedGF256RowMask());
    hash.U32(wirehair_v2::ActiveMixedGroupedGF256HashSeed());
    for (uint32_t row = 0u; row < gf256_rows; ++row) {
        for (uint32_t column = 0u; column < period; ++column) {
            hash.U8(rows->Subfield[row][column]);
        }
    }
    for (uint32_t row = 0u; row < gf16_rows; ++row) {
        for (uint32_t column = 0u; column < period; ++column) {
            hash.U16(rows->Extension[row][column]);
        }
    }
    hash.U32(2u * period);
    for (uint32_t column = 0u; column < 2u * period; ++column)
    {
        hash.U32(wirehair_v2::ActiveMixedCoefficientResidue(column));
        hash.U32(
            wirehair_v2::ActiveMixedExtensionCoefficientResidue(column));
    }
    return hash.Finish();
}

std::string HashSolveStats(const wirehair_v2::PrecodeSolveStats& stats)
{
    Hasher hash("WH2HTSOLVESTATS1");
    hash.U32(stats.PacketRows);
    hash.U32(stats.PeeledColumns);
    hash.U32(stats.InactivatedColumns);
    hash.U32(stats.ResidualRows);
    hash.U32(stats.ResidualRank);
    hash.U32(stats.BinaryResidualRank);
    hash.U64(stats.BinaryRowReferences);
    hash.U64(stats.BinaryRowStorageBytes);
    hash.U64(stats.BinaryAdjacencyStorageBytes);
    hash.U32(stats.BinaryRowStorageAllocations);
    hash.U32(stats.BinaryAdjacencyStorageAllocations);
    hash.U64(stats.BlockXors);
    hash.U64(stats.BlockMulAdds);
    hash.U64(stats.BlockCopies);
    hash.U64(stats.BlockZeroFills);
    hash.U64(stats.BlockAddSets);
    hash.U64(stats.BlockAddSetSources);
    hash.U64(stats.PeelAdjacencyVisits);
    hash.U64(stats.PeelRowScanSteps);
    hash.U64(stats.PeelHeapOperations);
    hash.U64(stats.ProjectionWordXors);
    hash.U64(stats.ResidualCoeffWordXors);
    hash.U64(stats.ResidualCoeffByteOps);
    hash.U32(stats.PacketSeedAttempt);
    return hash.Finish();
}

std::string HashPacketRows(
    uint32_t source_count,
    uint32_t precode_count,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    const std::vector<uint32_t>& ids)
{
    static const size_t kCapacity =
        64u + wirehair_v2::kCertifiedPacketMixCount;
    Hasher hash("WH2HTROWS1");
    hash.U32(source_count);
    hash.U32(precode_count);
    hash.U32(config.PeelSeed);
    hash.U32(config.MixCount);
    hash.U64((uint64_t)ids.size());
    for (uint32_t id : ids)
    {
        const std::vector<uint32_t> expected =
            wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                source_count, precode_count, id, config, runtime);
        uint32_t columns[kCapacity];
        size_t count = 0u;
        if (expected.empty() || expected.size() > kCapacity ||
            !wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
                source_count, precode_count, id, config, runtime,
                columns, kCapacity, count) ||
            count != expected.size() ||
            !std::equal(expected.begin(), expected.end(), columns))
        {
            return std::string();
        }
        hash.U32(id);
        hash.U32((uint32_t)count);
        for (size_t i = 0u; i < count; ++i) hash.U32(columns[i]);
    }
    return hash.Finish();
}

std::string HashTrace(const Options& options, const std::vector<uint32_t>& ids)
{
    Hasher hash("WH2HTTRACE1");
    hash.String(kSchema);
    hash.String("dispatch-v1");
    hash.U32(options.K);
    hash.U32(options.BlockBytes);
    hash.U64(options.ConstructionSeed);
    hash.U64(options.LossSeed);
    hash.U32(options.LossPpm);
    hash.String(ScheduleName(options.PacketSchedule));
    hash.U64((uint64_t)ids.size());
    for (uint32_t id : ids) hash.U32(id);
    return hash.Finish();
}

std::string HashPayload(
    const std::vector<uint32_t>& ids,
    const std::vector<uint8_t>& payload,
    uint32_t block_bytes)
{
    if ((uint64_t)ids.size() * block_bytes != payload.size()) {
        return std::string();
    }
    Hasher hash("WH2HTPAYLOAD1");
    hash.U32(block_bytes);
    hash.U64((uint64_t)ids.size());
    for (size_t i = 0u; i < ids.size(); ++i)
    {
        hash.U32(ids[i]);
        hash.U32(block_bytes);
        hash.Bytes(payload.data() + i * block_bytes, block_bytes);
    }
    return hash.Finish();
}

struct Fixture
{
    std::string Status = "error";
    WirehairResult EncoderResult = WirehairResult_Count;
    WirehairResult DecoderResult = WirehairResult_Count;
    WirehairResult RecoverResult = WirehairResult_Count;
    WirehairResult DirectResult = WirehairResult_Count;
    wirehair_v2::SeedProfile Profile;
    wirehair_v2::MessagePrecodeEncoderOptions CodecOptions;
    wirehair_v2::PrecodeParams Params;
    wirehair_v2::PacketRowConfig PacketConfig;
    wirehair_v2::PacketRowRuntime PacketRuntime;
    wirehair_v2::PrecodeSystem System;
    std::vector<uint8_t> Message;
    std::vector<uint32_t> CandidateIds;
    std::vector<uint32_t> DeliveredIds;
    std::vector<uint8_t> Payload;
    std::vector<wirehair_v2::SolvePacket> Packets;
    std::vector<uint8_t> Intermediate;
    wirehair_v2::PrecodeSolveStats EncoderStats = {};
    wirehair_v2::PrecodeSolveStats DecoderStats = {};
    wirehair_v2::PrecodeSolveStats DirectStats = {};
    uint32_t DecoderSolveAttempts = 0u;
    uint32_t TracePacketCount = 0u;
    std::string ProfileSha256;
    std::string SystemSha256;
    std::string CoefficientsSha256;
    std::string TraceSha256;
    std::string RowSha256;
    std::string MessageSha256;
    std::string PayloadSha256;
    std::string IntermediateSha256;
    std::string RecoveredSha256;
    std::string EncoderStatsSha256;
    std::string DecoderStatsSha256;
    std::string DirectStatsSha256;
    std::string SemanticSha256;
};

enum class PrepareResult
{
    Success,
    WeakRoot,
    Failure
};

struct UsageSnapshot
{
    uint64_t MinorFaults = 0u;
    uint64_t MajorFaults = 0u;
    uint64_t VoluntaryContextSwitches = 0u;
    uint64_t InvoluntaryContextSwitches = 0u;
};

struct TimingSample
{
    uint64_t ElapsedNanoseconds = 0u;
    UsageSnapshot Usage;
};

struct TimingResult
{
    const char* Scope = "";
    uint64_t WorkItemsPerRep = 0u;
    uint32_t MeasuredReps = 0u;
    uint64_t ElapsedNanoseconds = 0u;
    uint64_t MinNanoseconds = UINT64_MAX;
    uint64_t MaxNanoseconds = 0u;
    UsageSnapshot Usage;
    uint64_t Sink = 0u;
};

volatile uint64_t TimingSink = 0u;

bool MonotonicNanoseconds(uint64_t& value)
{
#if defined(__unix__) || defined(__APPLE__)
    struct timespec now = {};
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0 ||
        now.tv_sec < 0 || now.tv_nsec < 0 ||
        now.tv_nsec >= 1000000000L)
    {
        return false;
    }
    const uint64_t seconds = (uint64_t)now.tv_sec;
    if (seconds > (UINT64_MAX - (uint64_t)now.tv_nsec) /
            UINT64_C(1000000000))
    {
        return false;
    }
    value = seconds * UINT64_C(1000000000) + (uint64_t)now.tv_nsec;
    return true;
#else
    (void)value;
    return false;
#endif
}

bool CaptureUsage(UsageSnapshot& usage)
{
#if defined(__unix__) || defined(__APPLE__)
    struct rusage value = {};
    if (getrusage(RUSAGE_SELF, &value) != 0 ||
        value.ru_minflt < 0 || value.ru_majflt < 0 ||
        value.ru_nvcsw < 0 || value.ru_nivcsw < 0)
    {
        return false;
    }
    usage.MinorFaults = (uint64_t)value.ru_minflt;
    usage.MajorFaults = (uint64_t)value.ru_majflt;
    usage.VoluntaryContextSwitches = (uint64_t)value.ru_nvcsw;
    usage.InvoluntaryContextSwitches = (uint64_t)value.ru_nivcsw;
    return true;
#else
    (void)usage;
    return false;
#endif
}

bool UsageDelta(
    const UsageSnapshot& before,
    const UsageSnapshot& after,
    UsageSnapshot& delta)
{
    if (after.MinorFaults < before.MinorFaults ||
        after.MajorFaults < before.MajorFaults ||
        after.VoluntaryContextSwitches <
            before.VoluntaryContextSwitches ||
        after.InvoluntaryContextSwitches <
            before.InvoluntaryContextSwitches)
    {
        return false;
    }
    delta.MinorFaults = after.MinorFaults - before.MinorFaults;
    delta.MajorFaults = after.MajorFaults - before.MajorFaults;
    delta.VoluntaryContextSwitches =
        after.VoluntaryContextSwitches -
        before.VoluntaryContextSwitches;
    delta.InvoluntaryContextSwitches =
        after.InvoluntaryContextSwitches -
        before.InvoluntaryContextSwitches;
    return true;
}

template<class TimedOperation>
bool MeasureInvocation(
    int pinned_cpu,
    StartBarrier* start_barrier,
    const TimedOperation& operation,
    TimingSample& sample)
{
#if defined(__linux__)
    if (sched_getcpu() != pinned_cpu) return false;
#else
    (void)pinned_cpu;
    return false;
#endif
    if (start_barrier && !start_barrier->Rendezvous()) return false;
    UsageSnapshot before_usage;
    UsageSnapshot after_usage;
    uint64_t before_ns = 0u;
    uint64_t after_ns = 0u;
    if (!CaptureUsage(before_usage) ||
        !MonotonicNanoseconds(before_ns))
    {
        return false;
    }
    const bool operation_ok = operation();
    if (!MonotonicNanoseconds(after_ns) ||
        !CaptureUsage(after_usage))
    {
        return false;
    }
#if defined(__linux__)
    if (sched_getcpu() != pinned_cpu) return false;
#endif
    if (!operation_ok || after_ns < before_ns ||
        !UsageDelta(before_usage, after_usage, sample.Usage))
    {
        return false;
    }
    sample.ElapsedNanoseconds = after_ns - before_ns;
    return true;
}

void SetTimingSample(
    const TimingSample& sample,
    uint64_t sink,
    TimingResult& result)
{
    result.MeasuredReps = 1u;
    result.ElapsedNanoseconds = sample.ElapsedNanoseconds;
    result.MinNanoseconds = sample.ElapsedNanoseconds;
    result.MaxNanoseconds = sample.ElapsedNanoseconds;
    result.Usage = sample.Usage;
    result.Sink = sink;
    TimingSink ^= sink;
}

void FillMessage(std::vector<uint8_t>& message, const Options& options)
{
    Rng rng(
        options.ConstructionSeed ^
        (options.LossSeed + UINT64_C(0x5748324d53473131)));
    uint64_t value = 0u;
    for (size_t i = 0u; i < message.size(); ++i)
    {
        if ((i & 7u) == 0u) value = rng.Next();
        message[i] = (uint8_t)(value >> (8u * (i & 7u)));
    }
}

bool CheckRowInto(
    uint32_t source_count,
    uint32_t precode_count,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime)
{
    static const uint32_t kCanary = UINT32_C(0xdecafbad);
    static const size_t kCapacity = 64u + wirehair_v2::kCertifiedPacketMixCount;
    const uint32_t ids[] = {
        0u, 1u, source_count - 1u, source_count,
        source_count + 7u, UINT32_MAX
    };
    bool checked_short = false;
    for (size_t i = 0u; i < sizeof(ids) / sizeof(ids[0]); ++i)
    {
        const std::vector<uint32_t> expected =
            wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                source_count, precode_count, ids[i], config, runtime);
        if (expected.empty() || expected.size() > kCapacity) return false;
        uint32_t storage[kCapacity + 2u];
        std::fill(storage, storage + kCapacity + 2u, kCanary);
        size_t count = std::numeric_limits<size_t>::max();
        if (!wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
                source_count, precode_count, ids[i], config, runtime,
                storage + 1u, kCapacity, count) ||
            count != expected.size() ||
            !std::equal(expected.begin(), expected.end(), storage + 1u) ||
            storage[0] != kCanary ||
            !std::all_of(
                storage + 1u + count,
                storage + kCapacity + 2u,
                [](uint32_t value) { return value == kCanary; }))
        {
            return false;
        }
        if (!checked_short)
        {
            std::fill(storage, storage + kCapacity + 2u, kCanary);
            const size_t count_canary = std::numeric_limits<size_t>::max();
            count = count_canary;
            if (wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
                    source_count, precode_count, ids[i], config, runtime,
                    storage + 1u, expected.size() - 1u, count) ||
                count != count_canary ||
                !std::all_of(storage, storage + kCapacity + 2u,
                    [](uint32_t value) { return value == kCanary; }))
            {
                return false;
            }
            count = count_canary;
            if (wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
                    source_count, precode_count, ids[i], config, runtime,
                    nullptr, expected.size(), count) ||
                count != count_canary)
            {
                return false;
            }
            std::fill(storage, storage + kCapacity + 2u, kCanary);
            count = count_canary;
            if (wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
                    source_count, precode_count, ids[i], config, runtime,
                    storage + 1u, 0u, count) ||
                count != count_canary ||
                !std::all_of(storage, storage + kCapacity + 2u,
                    [](uint32_t value) { return value == kCanary; }))
            {
                return false;
            }
            count = count_canary;
            uint32_t* const aliased_columns =
                reinterpret_cast<uint32_t*>(&count);
            if (wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
                    source_count, precode_count, ids[i], config, runtime,
                    aliased_columns, kCapacity, count) ||
                count != count_canary)
            {
                return false;
            }
            checked_short = true;
        }
    }
    wirehair_v2::PacketRowRuntime invalid;
    uint32_t storage[kCapacity + 2u];
    std::fill(storage, storage + kCapacity + 2u, kCanary);
    size_t count = std::numeric_limits<size_t>::max();
    if (wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
            source_count, precode_count, 0u, config, invalid,
            storage + 1u, kCapacity, count) ||
        count != std::numeric_limits<size_t>::max() ||
        !std::all_of(storage, storage + kCapacity + 2u,
            [](uint32_t value) { return value == kCanary; }))
    {
        return false;
    }
    return checked_short;
}

void FinalizeSemantic(
    const wirehair_v2::V2EquationContract& contract,
    Fixture& fixture)
{
    Hasher semantic("WH2HTSEMANTIC1");
    semantic.String(kSchema);
    semantic.String(contract.CanonicalName);
    semantic.String(fixture.Status);
    semantic.U32((uint32_t)fixture.EncoderResult);
    semantic.U32((uint32_t)fixture.DecoderResult);
    semantic.U32((uint32_t)fixture.RecoverResult);
    semantic.U32((uint32_t)fixture.DirectResult);
    semantic.String(fixture.ProfileSha256);
    semantic.String(fixture.SystemSha256);
    semantic.String(fixture.CoefficientsSha256);
    semantic.String(fixture.TraceSha256);
    semantic.String(fixture.RowSha256);
    semantic.String(fixture.MessageSha256);
    semantic.String(fixture.PayloadSha256);
    semantic.String(fixture.IntermediateSha256);
    semantic.String(fixture.RecoveredSha256);
    semantic.String(fixture.EncoderStatsSha256);
    semantic.String(fixture.DecoderStatsSha256);
    semantic.String(fixture.DirectStatsSha256);
    semantic.U32(fixture.TracePacketCount);
    semantic.U32((uint32_t)fixture.DeliveredIds.size());
    semantic.U32(fixture.DecoderSolveAttempts);
    fixture.SemanticSha256 = semantic.Finish();
}

bool StructurallyWeakEncoderError(
    const Options& options,
    const Fixture& fixture)
{
    std::vector<uint8_t> zero(options.BlockBytes, 0u);
    std::vector<wirehair_v2::SolvePacket> packets(options.K);
    for (uint32_t id = 0u; id < options.K; ++id)
    {
        packets[id].BlockId = id;
        packets[id].Data = zero.data();
    }
    std::vector<uint8_t> intermediate;
    return wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
        fixture.System,
        fixture.PacketConfig,
        fixture.PacketRuntime,
        packets,
        options.BlockBytes,
        intermediate) == Wirehair_NeedMore;
}

PrepareResult PrepareFixture(
    const Options& options,
    const wirehair_v2::V2EquationContract& contract,
    Fixture& fixture)
{
    if (!CheckMemoryBudget(options) ||
        !wirehair_v2::IsCanonicalMixedCompletionState() ||
        !wirehair_v2::MakeRawContractProfile(
            contract, options.K, options.BlockBytes,
            options.ConstructionSeed, fixture.Profile))
    {
        return PrepareResult::Failure;
    }
    fixture.CodecOptions =
        wirehair_v2::MessageOptionsForContract(contract);
    if (fixture.CodecOptions.CacheSystematicSource ||
        fixture.CodecOptions.CacheReceivedSystematicPackets ||
        !wirehair_v2::ResolveMessagePrecodeConfiguration(
            fixture.Profile,
            fixture.CodecOptions,
            fixture.Params,
            fixture.PacketConfig) ||
        !wirehair_v2::BuildPrecodeSystem(
            fixture.Params, fixture.System) ||
        !wirehair_v2::ValidatePrecodeSystem(fixture.System))
    {
        return PrepareResult::Failure;
    }
    const uint32_t precode_count =
        fixture.Params.Staircase +
        fixture.Params.DenseRows +
        fixture.Params.HeavyRows;
    if (!fixture.PacketRuntime.Initialize(
            options.K, precode_count, fixture.PacketConfig.MixCount) ||
        !CheckRowInto(
            options.K, precode_count,
            fixture.PacketConfig, fixture.PacketRuntime))
    {
        return PrepareResult::Failure;
    }

    const uint64_t message_bytes_u64 =
        (uint64_t)options.K * options.BlockBytes;
    const uint64_t intermediate_bytes_u64 =
        (uint64_t)(options.K + precode_count) * options.BlockBytes;
    if (message_bytes_u64 > std::numeric_limits<size_t>::max() ||
        intermediate_bytes_u64 > std::numeric_limits<size_t>::max())
    {
        return PrepareResult::Failure;
    }
    const size_t message_bytes = (size_t)message_bytes_u64;
    const size_t intermediate_bytes = (size_t)intermediate_bytes_u64;
    fixture.Message.resize(message_bytes);
    FillMessage(fixture.Message, options);

    const uint32_t max_delivered = options.K * 2u + 512u;
    if (!BuildSchedule(
            options.K, max_delivered, options.LossPpm,
            options.LossSeed, options.PacketSchedule,
            fixture.CandidateIds))
    {
        return PrepareResult::Failure;
    }
    fixture.ProfileSha256 = HashProfile(fixture.Profile);
    fixture.SystemSha256 = HashSystem(fixture.System);
    fixture.CoefficientsSha256 = HashCoefficients();
    fixture.MessageSha256 = Sha256Bytes(
        fixture.Message.data(), fixture.Message.size(), "WH2HTBYTESEQ1");
    if (fixture.CoefficientsSha256.empty() ||
        fixture.MessageSha256.empty())
    {
        return PrepareResult::Failure;
    }

    wirehair_v2::MessagePrecodeEncoder encoder;
    fixture.EncoderResult = encoder.InitializeResult(
        fixture.Message.data(), message_bytes_u64, options.BlockBytes,
        &fixture.Profile, &fixture.CodecOptions);
    if (fixture.EncoderResult != Wirehair_Success)
    {
        const bool weak =
            fixture.EncoderResult == Wirehair_BadPeelSeed ||
            fixture.EncoderResult == Wirehair_NeedMore ||
            (fixture.EncoderResult == Wirehair_Error &&
             StructurallyWeakEncoderError(options, fixture));
        if (!weak)
        {
            std::fprintf(stderr, "raw encoder initialization failed: %s\n",
                wirehair_result_string(fixture.EncoderResult));
            return PrepareResult::Failure;
        }
        fixture.Status = "weak-root";
        fixture.TracePacketCount =
            (uint32_t)fixture.CandidateIds.size();
        fixture.TraceSha256 =
            HashTrace(options, fixture.CandidateIds);
        fixture.RowSha256 = HashPacketRows(
            options.K, precode_count,
            fixture.PacketConfig, fixture.PacketRuntime,
            fixture.CandidateIds);
        if (fixture.TraceSha256.empty() || fixture.RowSha256.empty()) {
            return PrepareResult::Failure;
        }
        fixture.PayloadSha256 = "not_applicable";
        fixture.IntermediateSha256 = "not_applicable";
        fixture.RecoveredSha256 = "not_applicable";
        fixture.EncoderStatsSha256 = "not_applicable";
        fixture.DecoderStatsSha256 = "not_applicable";
        fixture.DirectStatsSha256 = "not_applicable";
        FinalizeSemantic(contract, fixture);
        return PrepareResult::WeakRoot;
    }
    if (
        !encoder.IsInitialized() ||
        encoder.SourceBlockCount() != options.K ||
        encoder.BlockBytes() != options.BlockBytes ||
        encoder.IntermediateBlocks() == nullptr)
    {
        return PrepareResult::Failure;
    }
    const uint64_t payload_capacity_u64 =
        (uint64_t)max_delivered * options.BlockBytes;
    if (payload_capacity_u64 > std::numeric_limits<size_t>::max()) {
        return PrepareResult::Failure;
    }
    fixture.Payload.resize((size_t)payload_capacity_u64);
    fixture.DeliveredIds.clear();
    fixture.DeliveredIds.reserve(max_delivered);

    wirehair_v2::MessagePrecodeDecoder decoder;
    fixture.DecoderResult = decoder.InitializeResult(
            message_bytes_u64, options.BlockBytes,
            &fixture.Profile, &fixture.CodecOptions);
    if (fixture.DecoderResult != Wirehair_Success)
    {
        return PrepareResult::Failure;
    }
    WirehairResult decode_result = Wirehair_NeedMore;
    for (uint32_t id : fixture.CandidateIds)
    {
        const size_t offset =
            fixture.DeliveredIds.size() * options.BlockBytes;
        uint32_t data_bytes = UINT32_MAX;
        if (encoder.EncodeResult(
                id, fixture.Payload.data() + offset, options.BlockBytes,
                &data_bytes, nullptr) != Wirehair_Success ||
            data_bytes != options.BlockBytes)
        {
            return PrepareResult::Failure;
        }
        fixture.DeliveredIds.push_back(id);
        decode_result = decoder.DecodeResult(
            id, fixture.Payload.data() + offset, data_bytes);
        if (decode_result == Wirehair_Success) break;
        if (decode_result != Wirehair_NeedMore) {
            return PrepareResult::Failure;
        }
    }
    if (decode_result != Wirehair_Success ||
        fixture.DeliveredIds.size() < options.K)
    {
        std::fprintf(stderr,
            "decode trace did not recover within %u delivered packets\n",
            max_delivered);
        return PrepareResult::Failure;
    }
    fixture.DecoderResult = decode_result;
    fixture.Payload.resize(
        fixture.DeliveredIds.size() * options.BlockBytes);

    std::vector<uint8_t> recovered(message_bytes, 0u);
    fixture.RecoverResult = decoder.RecoverResult(
        recovered.data(), message_bytes_u64);
    if (fixture.RecoverResult != Wirehair_Success ||
        recovered != fixture.Message ||
        decoder.IntermediateBlocks() == nullptr)
    {
        return PrepareResult::Failure;
    }

    fixture.Packets.resize(fixture.DeliveredIds.size());
    for (size_t i = 0u; i < fixture.Packets.size(); ++i)
    {
        fixture.Packets[i].BlockId = fixture.DeliveredIds[i];
        fixture.Packets[i].Data =
            fixture.Payload.data() + i * options.BlockBytes;
    }
    fixture.DirectResult =
        wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
            fixture.System,
            fixture.PacketConfig,
            fixture.PacketRuntime,
            fixture.Packets,
            options.BlockBytes,
            fixture.Intermediate,
            &fixture.DirectStats);
    if (fixture.DirectResult != Wirehair_Success ||
        fixture.Intermediate.size() != intermediate_bytes ||
        std::memcmp(
            fixture.Intermediate.data(),
            encoder.IntermediateBlocks(),
            intermediate_bytes) != 0 ||
        std::memcmp(
            fixture.Intermediate.data(),
            decoder.IntermediateBlocks(),
            intermediate_bytes) != 0)
    {
        return PrepareResult::Failure;
    }

    fixture.EncoderStats = encoder.SolveStats();
    fixture.DecoderStats = decoder.SolveStats();
    fixture.DecoderSolveAttempts = decoder.SolveAttemptCount();
    if (encoder.BlockEncoder().HasCompleteSystem() ||
        HashParamsOnly(fixture.System.Params) !=
            HashParamsOnly(encoder.BlockEncoder().System().Params) ||
        fixture.SystemSha256 != HashSystem(decoder.System()))
    {
        return PrepareResult::Failure;
    }
    fixture.TraceSha256 = HashTrace(options, fixture.DeliveredIds);
    fixture.TracePacketCount =
        (uint32_t)fixture.DeliveredIds.size();
    fixture.RowSha256 = HashPacketRows(
        options.K, precode_count,
        fixture.PacketConfig, fixture.PacketRuntime,
        fixture.DeliveredIds);
    fixture.MessageSha256 = Sha256Bytes(
        fixture.Message.data(), fixture.Message.size(), "WH2HTBYTESEQ1");
    fixture.PayloadSha256 = HashPayload(
        fixture.DeliveredIds, fixture.Payload, options.BlockBytes);
    fixture.IntermediateSha256 = Sha256Bytes(
        fixture.Intermediate.data(), fixture.Intermediate.size(),
        "WH2HTBYTESEQ1");
    fixture.RecoveredSha256 = Sha256Bytes(
        recovered.data(), recovered.size(), "WH2HTBYTESEQ1");
    fixture.EncoderStatsSha256 = HashSolveStats(fixture.EncoderStats);
    fixture.DecoderStatsSha256 = HashSolveStats(fixture.DecoderStats);
    fixture.DirectStatsSha256 = HashSolveStats(fixture.DirectStats);
    if (fixture.CoefficientsSha256.empty() ||
        fixture.RowSha256.empty() ||
        fixture.PayloadSha256.empty() ||
        fixture.MessageSha256 != fixture.RecoveredSha256)
    {
        return PrepareResult::Failure;
    }
    fixture.Status = "success";
    FinalizeSemantic(contract, fixture);
    return PrepareResult::Success;
}

bool RunRowBatch(
    const Options& options,
    const Fixture& fixture,
    int pinned_cpu,
    uint32_t repetitions,
    StartBarrier* start_barrier,
    TimingSample& sample,
    uint64_t& sink)
{
    static const size_t kCapacity =
        64u + wirehair_v2::kCertifiedPacketMixCount;
    uint32_t columns[kCapacity];
    size_t final_count = 0u;
    const uint32_t precode_count =
        fixture.Params.Staircase +
        fixture.Params.DenseRows +
        fixture.Params.HeavyRows;
    const bool measured = MeasureInvocation(
        pinned_cpu,
        start_barrier,
        [&]() {
            for (uint32_t repetition = 0u;
                 repetition < repetitions;
                 ++repetition)
            {
                for (uint32_t id : fixture.DeliveredIds)
                {
                    if (!wirehair_v2::GeneratePacketMatrixRowIntoWithRuntime(
                            options.K,
                            precode_count,
                            id,
                            fixture.PacketConfig,
                            fixture.PacketRuntime,
                            columns,
                            kCapacity,
                            final_count))
                    {
                        return false;
                    }
                }
            }
            return true;
        },
        sample);
    if (!measured) return false;
    sink = UINT64_C(0x574832524f573131) ^
        ((uint64_t)repetitions << 48) ^
        ((uint64_t)fixture.DeliveredIds.back() << 16) ^
        (uint64_t)columns[0] ^
        ((uint64_t)columns[final_count - 1u] << 32) ^
        final_count;
    return measured;
}

bool RunEncoderBatch(
    const Options& options,
    const Fixture& fixture,
    int pinned_cpu,
    uint32_t repetitions,
    StartBarrier* start_barrier,
    TimingSample& sample,
    uint64_t& sink)
{
    wirehair_v2::MessagePrecodeEncoder encoder;
    std::vector<uint8_t> output(fixture.Message.size(), 0u);
    const bool measured = MeasureInvocation(
        pinned_cpu,
        start_barrier,
        [&]() {
            for (uint32_t repetition = 0u;
                 repetition < repetitions;
                 ++repetition)
            {
                if (encoder.InitializeResult(
                        fixture.Message.data(),
                        fixture.Message.size(),
                        options.BlockBytes,
                        &fixture.Profile,
                        &fixture.CodecOptions) != Wirehair_Success)
                {
                    return false;
                }
                for (uint32_t id = 0u; id < options.K; ++id)
                {
                    uint32_t data_bytes = UINT32_MAX;
                    if (encoder.EncodeResult(
                            id,
                            output.data() +
                                (size_t)id * options.BlockBytes,
                            options.BlockBytes,
                            &data_bytes,
                            nullptr) != Wirehair_Success ||
                        data_bytes != options.BlockBytes)
                    {
                        return false;
                    }
                }
            }
            return true;
        },
        sample);
    if (!measured ||
        output != fixture.Message ||
        encoder.IntermediateBlocks() == nullptr ||
        std::memcmp(
            encoder.IntermediateBlocks(),
            fixture.Intermediate.data(),
            fixture.Intermediate.size()) != 0 ||
        HashSolveStats(encoder.SolveStats()) !=
            fixture.EncoderStatsSha256)
    {
        return false;
    }
    sink = UINT64_C(0x574832454e433131) ^
        ((uint64_t)repetitions << 48) ^
        encoder.SolveStats().BlockXors ^
        ((uint64_t)output.front() << 8) ^
        output.back();
    return true;
}

bool FeedBatchFitsBudget(
    const Options& options,
    const Fixture& fixture,
    uint32_t repetitions)
{
    const uint64_t total_blocks =
        (uint64_t)options.K +
        fixture.Params.Staircase +
        fixture.Params.DenseRows +
        fixture.Params.HeavyRows;
    const uint64_t bytes_per_block =
        (uint64_t)options.BlockBytes * 3u + 64u;
    if (total_blocks >
        (UINT64_MAX - UINT64_C(1048576)) /
            bytes_per_block)
    {
        return false;
    }
    const uint64_t per_decoder =
        total_blocks * bytes_per_block +
        UINT64_C(1048576);
    if (repetitions != 0u &&
        per_decoder > UINT64_MAX / repetitions)
    {
        return false;
    }
    const uint64_t decoder_bytes = per_decoder * repetitions;
    const uint64_t base_factor = (uint64_t)options.K * 8u + 4096u;
    if (base_factor >
        (UINT64_MAX - (uint64_t)options.K * 64u -
            UINT64_C(16777216)) / options.BlockBytes)
    {
        return false;
    }
    const uint64_t base_bytes =
        base_factor * options.BlockBytes +
        (uint64_t)options.K * 64u +
        UINT64_C(16777216);
    return decoder_bytes <= UINT64_MAX - base_bytes &&
        options.MaxWorkingMiB <= (UINT64_MAX >> 20) &&
        decoder_bytes + base_bytes <=
            (options.MaxWorkingMiB << 20);
}

bool RunDecoderFeedBatch(
    const Options& options,
    const Fixture& fixture,
    int pinned_cpu,
    uint32_t repetitions,
    StartBarrier* start_barrier,
    TimingSample& sample,
    uint64_t& sink)
{
    if (!FeedBatchFitsBudget(options, fixture, repetitions)) {
        std::fprintf(stderr,
            "decoder-feed batch exceeds --max-working-mib\n");
        return false;
    }
    std::vector<
        std::unique_ptr<wirehair_v2::MessagePrecodeDecoder> > decoders;
    decoders.reserve(repetitions);
    for (uint32_t repetition = 0u; repetition < repetitions; ++repetition)
    {
        std::unique_ptr<wirehair_v2::MessagePrecodeDecoder> decoder(
            new wirehair_v2::MessagePrecodeDecoder);
        if (decoder->InitializeResult(
                fixture.Message.size(),
                options.BlockBytes,
                &fixture.Profile,
                &fixture.CodecOptions) != Wirehair_Success)
        {
            return false;
        }
        decoders.push_back(std::move(decoder));
    }
    const bool measured = MeasureInvocation(
        pinned_cpu,
        start_barrier,
        [&]() {
            for (uint32_t repetition = 0u;
                 repetition < repetitions;
                 ++repetition)
            {
                wirehair_v2::MessagePrecodeDecoder& decoder =
                    *decoders[repetition];
                for (size_t i = 0u;
                     i < fixture.DeliveredIds.size();
                     ++i)
                {
                    const WirehairResult result = decoder.DecodeResult(
                        fixture.DeliveredIds[i],
                        fixture.Payload.data() +
                            i * options.BlockBytes,
                        options.BlockBytes);
                    const WirehairResult expected =
                        i + 1u == fixture.DeliveredIds.size() ?
                        Wirehair_Success : Wirehair_NeedMore;
                    if (result != expected) return false;
                }
            }
            return true;
        },
        sample);
    if (!measured) return false;
    std::vector<uint8_t> recovered(fixture.Message.size(), 0u);
    for (const std::unique_ptr<
             wirehair_v2::MessagePrecodeDecoder>& decoder : decoders)
    {
        if (decoder->RecoverResult(
                recovered.data(), recovered.size()) != Wirehair_Success ||
            recovered != fixture.Message ||
            decoder->IntermediateBlocks() == nullptr ||
            std::memcmp(
                decoder->IntermediateBlocks(),
                fixture.Intermediate.data(),
                fixture.Intermediate.size()) != 0 ||
            HashSolveStats(decoder->SolveStats()) !=
                fixture.DecoderStatsSha256)
        {
            return false;
        }
    }
    sink = UINT64_C(0x5748324645454431) ^
        ((uint64_t)repetitions << 48) ^
        ((uint64_t)decoders.back()->SolveAttemptCount() << 32) ^
        decoders.back()->ReceivedCount();
    return true;
}

bool RunDecoderFullBatch(
    const Options& options,
    const Fixture& fixture,
    int pinned_cpu,
    uint32_t repetitions,
    StartBarrier* start_barrier,
    TimingSample& sample,
    uint64_t& sink)
{
    wirehair_v2::MessagePrecodeDecoder decoder;
    std::vector<uint8_t> recovered(fixture.Message.size(), 0u);
    const bool measured = MeasureInvocation(
        pinned_cpu,
        start_barrier,
        [&]() {
            for (uint32_t repetition = 0u;
                 repetition < repetitions;
                 ++repetition)
            {
                if (decoder.InitializeResult(
                        fixture.Message.size(),
                        options.BlockBytes,
                        &fixture.Profile,
                        &fixture.CodecOptions) != Wirehair_Success)
                {
                    return false;
                }
                for (size_t i = 0u;
                     i < fixture.DeliveredIds.size();
                     ++i)
                {
                    const WirehairResult result = decoder.DecodeResult(
                        fixture.DeliveredIds[i],
                        fixture.Payload.data() +
                            i * options.BlockBytes,
                        options.BlockBytes);
                    const WirehairResult expected =
                        i + 1u == fixture.DeliveredIds.size() ?
                        Wirehair_Success : Wirehair_NeedMore;
                    if (result != expected) return false;
                }
                if (decoder.RecoverResult(
                        recovered.data(), recovered.size()) !=
                    Wirehair_Success)
                {
                    return false;
                }
            }
            return true;
        },
        sample);
    if (!measured ||
        recovered != fixture.Message ||
        decoder.IntermediateBlocks() == nullptr ||
        std::memcmp(
            decoder.IntermediateBlocks(),
            fixture.Intermediate.data(),
            fixture.Intermediate.size()) != 0 ||
        HashSolveStats(decoder.SolveStats()) !=
            fixture.DecoderStatsSha256)
    {
        return false;
    }
    sink = UINT64_C(0x57483246554c4c31) ^
        ((uint64_t)repetitions << 48) ^
        ((uint64_t)decoder.SolveAttemptCount() << 32) ^
        decoder.ReceivedCount() ^
        ((uint64_t)recovered.front() << 8) ^
        recovered.back();
    return true;
}

bool RunDirectBatch(
    const Options& options,
    const Fixture& fixture,
    int pinned_cpu,
    uint32_t repetitions,
    StartBarrier* start_barrier,
    TimingSample& sample,
    uint64_t& sink)
{
    std::vector<uint8_t> intermediate;
    intermediate.reserve(fixture.Intermediate.size());
    wirehair_v2::PrecodeSolveStats stats;
    const bool measured = MeasureInvocation(
        pinned_cpu,
        start_barrier,
        [&]() {
            for (uint32_t repetition = 0u;
                 repetition < repetitions;
                 ++repetition)
            {
                stats = wirehair_v2::PrecodeSolveStats();
                if (wirehair_v2::
                        SolvePrecodeSystemForValidatedSystemWithRuntime(
                            fixture.System,
                            fixture.PacketConfig,
                            fixture.PacketRuntime,
                            fixture.Packets,
                            options.BlockBytes,
                            intermediate,
                            &stats) != Wirehair_Success)
                {
                    return false;
                }
            }
            return true;
        },
        sample);
    if (!measured ||
        intermediate != fixture.Intermediate ||
        HashSolveStats(stats) != fixture.DirectStatsSha256)
    {
        return false;
    }
    sink = UINT64_C(0x5748324449523131) ^
        ((uint64_t)repetitions << 48) ^
        stats.BlockXors ^
        ((uint64_t)stats.PacketRows << 32) ^
        ((uint64_t)intermediate.front() << 8) ^
        intermediate.back();
    return true;
}

template<class RunBatch>
bool RunScope(
    const Options& options,
    const char* name,
    uint64_t work_items_per_rep,
    int pinned_cpu,
    StartBarrier* start_barrier,
    const RunBatch& run_batch,
    TimingResult& result)
{
    result.Scope = name;
    result.WorkItemsPerRep = work_items_per_rep;
    for (uint32_t warmup = 0u;
         warmup < options.WarmupReps;
         ++warmup)
    {
        TimingSample ignored_sample;
        uint64_t ignored_sink = 0u;
        if (!run_batch(
                options, pinned_cpu, 1u,
                nullptr,
                ignored_sample, ignored_sink))
        {
            return false;
        }
        TimingSink ^= ignored_sink;
    }
    TimingSample sample;
    uint64_t sink = 0u;
    if (!run_batch(
            options, pinned_cpu, options.InnerReps,
            start_barrier,
            sample, sink))
    {
        return false;
    }
    SetTimingSample(sample, sink, result);
    return true;
}

bool RunSelectedTimings(
    const Options& options,
    const Fixture& fixture,
    int pinned_cpu,
    StartBarrier* start_barrier,
    std::vector<TimingResult>& results)
{
    results.clear();
    const auto row = [&](const Options& local, int cpu, uint32_t reps,
                         StartBarrier* barrier,
                         TimingSample& sample, uint64_t& sink) {
        return RunRowBatch(
            local, fixture, cpu, reps, barrier, sample, sink);
    };
    const auto encoder = [&](const Options& local, int cpu, uint32_t reps,
                             StartBarrier* barrier,
                             TimingSample& sample, uint64_t& sink) {
        return RunEncoderBatch(
            local, fixture, cpu, reps, barrier, sample, sink);
    };
    const auto feed = [&](const Options& local, int cpu, uint32_t reps,
                          StartBarrier* barrier,
                          TimingSample& sample, uint64_t& sink) {
        return RunDecoderFeedBatch(
            local, fixture, cpu, reps, barrier, sample, sink);
    };
    const auto full = [&](const Options& local, int cpu, uint32_t reps,
                          StartBarrier* barrier,
                          TimingSample& sample, uint64_t& sink) {
        return RunDecoderFullBatch(
            local, fixture, cpu, reps, barrier, sample, sink);
    };
    const auto direct = [&](const Options& local, int cpu, uint32_t reps,
                            StartBarrier* barrier,
                            TimingSample& sample, uint64_t& sink) {
        return RunDirectBatch(
            local, fixture, cpu, reps, barrier, sample, sink);
    };
    if ((options.Scopes & ScopeRow) != 0u)
    {
        TimingResult result;
        if (!RunScope(
                options, "row", fixture.DeliveredIds.size(),
                pinned_cpu, start_barrier, row, result))
        {
            return false;
        }
        results.push_back(result);
    }
    if ((options.Scopes & ScopeEncoder) != 0u)
    {
        TimingResult result;
        if (!RunScope(
                options, "encoder", options.K,
                pinned_cpu, start_barrier, encoder, result))
        {
            return false;
        }
        results.push_back(result);
    }
    if ((options.Scopes & ScopeDecoderFeed) != 0u)
    {
        TimingResult result;
        if (!RunScope(
                options, "decoder-feed", fixture.DeliveredIds.size(),
                pinned_cpu, start_barrier, feed, result))
        {
            return false;
        }
        results.push_back(result);
    }
    if ((options.Scopes & ScopeDecoderFull) != 0u)
    {
        TimingResult result;
        if (!RunScope(
                options, "decoder-full", fixture.DeliveredIds.size(),
                pinned_cpu, start_barrier, full, result))
        {
            return false;
        }
        results.push_back(result);
    }
    if ((options.Scopes & ScopeDirect) != 0u)
    {
        TimingResult result;
        if (!RunScope(
                options, "direct", fixture.DeliveredIds.size(),
                pinned_cpu, start_barrier, direct, result))
        {
            return false;
        }
        results.push_back(result);
    }
    return true;
}

bool EmitReceipt(
    const Options& options,
    const Fixture& fixture,
    int cpu,
    const std::vector<TimingResult>& timings)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    static const uint32_t kHooksCompiled = 1u;
#else
    static const uint32_t kHooksCompiled = 0u;
#endif
    std::ostringstream body;
    body
        << "{\"record\":\"header\",\"schema\":\"" << kSchema
        << "\",\"profile\":\"dispatch-v1\",\"hooks_compiled\":"
        << kHooksCompiled
        << ",\"K\":" << options.K
        << ",\"block_bytes\":" << options.BlockBytes
        << ",\"construction_seed\":" << options.ConstructionSeed
        << ",\"loss_seed\":" << options.LossSeed
        << ",\"loss_ppm\":" << options.LossPpm
        << ",\"schedule\":\"" << ScheduleName(options.PacketSchedule)
        << "\",\"scope_request\":\"" << ScopeRequestName(options.Scopes)
        << "\",\"warmup_reps\":" << options.WarmupReps
        << ",\"inner_reps\":" << options.InnerReps
        << ",\"max_working_mib\":" << options.MaxWorkingMiB
        << ",\"context_sha256\":\"" << options.ContextSha256
        << "\",\"cpu\":" << cpu
        << ",\"clock\":\"CLOCK_MONOTONIC\""
        << ",\"start_barrier\":\"" << StartBarrierName(options)
        << "\"}\n";
    body
        << "{\"record\":\"semantic\",\"schema\":\"" << kSchema
        << "\",\"profile\":\"dispatch-v1\",\"status\":\""
        << fixture.Status
        << "\",\"encoder_result\":" << (uint32_t)fixture.EncoderResult
        << ",\"decoder_result\":";
    if (fixture.DecoderResult == WirehairResult_Count) {
        body << "null";
    } else {
        body << (uint32_t)fixture.DecoderResult;
    }
    body << ",\"recover_result\":";
    if (fixture.RecoverResult == WirehairResult_Count) {
        body << "null";
    } else {
        body << (uint32_t)fixture.RecoverResult;
    }
    body << ",\"direct_result\":";
    if (fixture.DirectResult == WirehairResult_Count) {
        body << "null";
    } else {
        body << (uint32_t)fixture.DirectResult;
    }
    body
        << ",\"trace_packets\":" << fixture.TracePacketCount
        << ",\"delivered_packets\":" << fixture.DeliveredIds.size()
        << ",\"overhead_packets\":";
    if (fixture.Status == "success") {
        body << fixture.DeliveredIds.size() - options.K;
    } else {
        body << "null";
    }
    body
        << ",\"decoder_solve_attempts\":" << fixture.DecoderSolveAttempts
        << ",\"profile_sha256\":\"" << fixture.ProfileSha256
        << "\",\"system_sha256\":\"" << fixture.SystemSha256
        << "\",\"coefficients_sha256\":\""
        << fixture.CoefficientsSha256
        << "\",\"trace_sha256\":\"" << fixture.TraceSha256
        << "\",\"row_sha256\":\"" << fixture.RowSha256
        << "\",\"message_sha256\":\"" << fixture.MessageSha256
        << "\",\"payload_sha256\":\"" << fixture.PayloadSha256
        << "\",\"intermediate_sha256\":\""
        << fixture.IntermediateSha256
        << "\",\"recovered_sha256\":\"" << fixture.RecoveredSha256
        << "\",\"encoder_stats_sha256\":\""
        << fixture.EncoderStatsSha256
        << "\",\"decoder_stats_sha256\":\""
        << fixture.DecoderStatsSha256
        << "\",\"direct_stats_sha256\":\""
        << fixture.DirectStatsSha256
        << "\",\"semantic_sha256\":\"" << fixture.SemanticSha256
        << "\"}\n";
    for (const TimingResult& timing : timings)
    {
        body
            << "{\"record\":\"timing\",\"schema\":\"" << kSchema
            << "\",\"profile\":\"dispatch-v1\",\"hooks_compiled\":"
            << kHooksCompiled
            << ",\"scope\":\"" << timing.Scope
            << "\",\"lifecycle\":\"" << TimingLifecycle(timing.Scope)
            << "\",\"semantic_sha256\":\"" << fixture.SemanticSha256
            << "\",\"unit\":\"ns\",\"clock\":\"CLOCK_MONOTONIC\""
            << ",\"warmup_reps\":" << options.WarmupReps
            << ",\"measured_reps\":" << timing.MeasuredReps
            << ",\"inner_reps\":" << options.InnerReps
            << ",\"work_items_per_rep\":" << timing.WorkItemsPerRep
            << ",\"elapsed_ns\":" << timing.ElapsedNanoseconds
            << ",\"min_ns\":" << timing.MinNanoseconds
            << ",\"max_ns\":" << timing.MaxNanoseconds
            << ",\"minor_faults\":" << timing.Usage.MinorFaults
            << ",\"major_faults\":" << timing.Usage.MajorFaults
            << ",\"voluntary_context_switches\":"
            << timing.Usage.VoluntaryContextSwitches
            << ",\"involuntary_context_switches\":"
            << timing.Usage.InvoluntaryContextSwitches
            << ",\"sink\":" << timing.Sink
            << ",\"result\":" << (uint32_t)Wirehair_Success
            << "}\n";
    }
    const std::string prefix = body.str();
    body
        << "{\"record\":\"done\",\"schema\":\"" << kSchema
        << "\",\"status\":\"" << fixture.Status
        << "\",\"records_before_done\":" << 2u + timings.size()
        << ",\"stream_sha256\":\"" << RawSha256(prefix)
        << "\"}\n";
    const std::string output = body.str();
    return std::fwrite(
            output.data(), 1u, output.size(), stdout) == output.size() &&
        std::fflush(stdout) == 0;
}

int Main(int argc, char** argv)
{
    if (argc == 2 && !std::strcmp(argv[1], "--help")) {
        Usage(argv[0]);
        return 0;
    }
    Options options;
    if (!ParseOptions(argc, argv, options)) {
        Usage(argv[0]);
        return 2;
    }
    StartBarrier start_barrier(options.ReadyFd, options.GoFd);
    if (!start_barrier.Validate()) {
        std::fprintf(stderr, "start barrier descriptor validation failed\n");
        return 2;
    }
    if (!RejectEquationEnvironment()) return 2;
    if (!CheckSha256()) {
        std::fprintf(stderr, "SHA-256 self-check failed\n");
        return 1;
    }
    if (wirehair_init() != Wirehair_Success) {
        std::fprintf(stderr, "wirehair_init failed\n");
        return 1;
    }
    int cpu = -1;
    if (!PinToFirstAllowedCpu(cpu)) {
        std::fprintf(stderr, "CPU affinity pinning failed\n");
        return 1;
    }
    const wirehair_v2::V2EquationContract* contract =
        wirehair_v2::FindV2EquationContract("dispatch-v1");
    if (!contract ||
        contract->ProfileId != wirehair_v2::kWh2DispatchV1ContractId ||
        contract->PublicProfile ||
        contract->SeedPolicy != wirehair_v2::V2SeedDerivation::RawUniform ||
        contract->SeedAttemptCount != 1u)
    {
        std::fprintf(stderr, "dispatch-v1 contract lookup failed\n");
        return 1;
    }
    Fixture fixture;
    const PrepareResult prepared =
        PrepareFixture(options, *contract, fixture);
    if (prepared == PrepareResult::Failure)
    {
        std::fprintf(stderr, "semantic preflight failed\n");
        return 1;
    }
    std::vector<TimingResult> timings;
    StartBarrier* measured_barrier =
        start_barrier.Enabled() ? &start_barrier : nullptr;
    if (prepared == PrepareResult::Success &&
        !RunSelectedTimings(
            options, fixture, cpu, measured_barrier, timings))
    {
        std::fprintf(stderr, "timing scope failed\n");
        return 1;
    }
    if (prepared == PrepareResult::WeakRoot &&
        start_barrier.Enabled() &&
        !start_barrier.Rendezvous())
    {
        std::fprintf(stderr, "weak-root start barrier failed\n");
        return 1;
    }
    if (start_barrier.Enabled() && !start_barrier.Consumed()) {
        std::fprintf(stderr, "start barrier was not consumed\n");
        return 1;
    }
    if (prepared == PrepareResult::WeakRoot && !timings.empty()) {
        return 1;
    }
    if (!EmitReceipt(options, fixture, cpu, timings)) {
        std::fprintf(stderr, "receipt write failed\n");
        return 1;
    }
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        return Main(argc, argv);
    }
    catch (const std::bad_alloc&) {
        std::fprintf(stderr, "out of memory\n");
    }
    catch (const std::length_error&) {
        std::fprintf(stderr, "allocation length overflow\n");
    }
    catch (const std::exception& exception) {
        std::fprintf(stderr, "exception: %s\n", exception.what());
    }
    return 1;
}
