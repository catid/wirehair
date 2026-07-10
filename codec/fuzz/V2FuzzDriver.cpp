#include "V2FuzzDriver.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <new>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

namespace wirehair_v2 {
namespace fuzz {
namespace {

using Clock = std::chrono::steady_clock;

struct Options
{
    std::string CorpusManifest;
    std::string ArtifactDirectory = "v2-fuzz-artifacts";
    std::string ReplayPath;
    uint64_t Seed = UINT64_C(0x6a09e667f3bcc909);
    uint64_t Mutations = 10000u;
    uint64_t ReplayIndex = 0u;
    uint32_t MaximumSeconds = 55u;
    uint32_t DurationSeconds = 0u;
    bool HaveReplayIndex = false;
    bool TraceCurrent = false;
};

uint64_t NextRandom(uint64_t& state)
{
    uint64_t z = (state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

bool ParseU64(const char* text, uint64_t& value)
{
    if (!text || !*text || *text == '-' || *text == '+') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(text, &end, 0);
    if (errno != 0 || !end || end == text || *end != '\0') {
        return false;
    }
    value = (uint64_t)parsed;
    return true;
}

bool ParseU32(const char* text, uint32_t& value)
{
    uint64_t parsed = 0u;
    if (!ParseU64(text, parsed) || parsed > UINT32_MAX) {
        return false;
    }
    value = (uint32_t)parsed;
    return true;
}

bool ReadFile(
    const std::string& path,
    size_t maximum_bytes,
    std::vector<uint8_t>& output,
    std::string& error)
{
    std::ifstream input(path.c_str(), std::ios::binary | std::ios::ate);
    if (!input) {
        error = "cannot open " + path;
        return false;
    }
    const std::streamoff end = input.tellg();
    if (end < 0 || (uint64_t)end > maximum_bytes) {
        error = "oversized input " + path;
        return false;
    }
    output.assign((size_t)end, 0u);
    input.seekg(0, std::ios::beg);
    if (!output.empty()) {
        input.read(
            reinterpret_cast<char*>(output.data()),
            (std::streamsize)output.size());
    }
    if (!input) {
        error = "cannot read " + path;
        return false;
    }
    return true;
}

std::string ParentDirectory(const std::string& path)
{
    const size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? std::string(".") :
        path.substr(0u, slash);
}

bool LoadCorpus(
    const std::string& manifest_path,
    std::vector<std::vector<uint8_t> >& corpus,
    std::string& error)
{
    std::vector<uint8_t> manifest_bytes;
    if (!ReadFile(
            manifest_path, kMaxCorpusArtifactBytes,
            manifest_bytes, error))
    {
        return false;
    }
    const std::string text = manifest_bytes.empty() ? std::string() :
        std::string(
            reinterpret_cast<const char*>(manifest_bytes.data()),
            manifest_bytes.size());
    const std::string directory = ParentDirectory(manifest_path);
    size_t total_bytes = manifest_bytes.size();
    size_t offset = 0u;
    while (offset <= text.size())
    {
        const size_t end = text.find('\n', offset);
        const size_t length = end == std::string::npos ?
            text.size() - offset : end - offset;
        std::string line = text.substr(offset, length);
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && line[0] != '#')
        {
            std::vector<uint8_t> entry;
            if (!ReadFile(
                    directory + "/" + line,
                    kMaxFuzzInputBytes,
                    entry,
                    error))
            {
                return false;
            }
            if (total_bytes > kMaxCorpusArtifactBytes ||
                entry.size() > kMaxCorpusArtifactBytes - total_bytes)
            {
                error = "corpus artifacts exceed 5 MiB";
                return false;
            }
            total_bytes += entry.size();
            corpus.push_back(std::move(entry));
        }
        if (end == std::string::npos) {
            break;
        }
        offset = end + 1u;
    }
    if (corpus.empty()) {
        error = "corpus manifest is empty";
        return false;
    }
    return true;
}

std::vector<uint8_t> Mutate(
    const std::vector<std::vector<uint8_t> >& corpus,
    uint64_t seed,
    uint64_t mutation_index)
{
    uint64_t state = seed ^
        (mutation_index + 1u) * UINT64_C(0xd6e8feb86659fd93);
    std::vector<uint8_t> value =
        corpus[(size_t)(NextRandom(state) % corpus.size())];
    if (value.size() > kMaxDeterministicMutationBytes) {
        value.resize(kMaxDeterministicMutationBytes);
    }
    const unsigned rounds = 1u + (unsigned)(NextRandom(state) % 8u);
    for (unsigned round = 0; round < rounds; ++round)
    {
        const unsigned operation = (unsigned)(NextRandom(state) % 8u);
        if (value.empty()) {
            value.push_back((uint8_t)NextRandom(state));
            continue;
        }
        const size_t position =
            (size_t)(NextRandom(state) % value.size());
        switch (operation)
        {
        case 0:
            value[position] ^= (uint8_t)(1u << (NextRandom(state) & 7u));
            break;
        case 1:
            value[position] = (uint8_t)NextRandom(state);
            break;
        case 2:
            if (value.size() < kMaxDeterministicMutationBytes) {
                value.insert(
                    value.begin() + (ptrdiff_t)position,
                    (uint8_t)NextRandom(state));
            }
            break;
        case 3:
            value.erase(value.begin() + (ptrdiff_t)position);
            break;
        case 4:
            value.resize(position + 1u);
            break;
        case 5:
            if (value.size() < kMaxDeterministicMutationBytes) {
                const size_t room =
                    kMaxDeterministicMutationBytes - value.size();
                const size_t count = std::min<size_t>(
                    room, 1u + (size_t)(NextRandom(state) % 32u));
                for (size_t i = 0; i < count; ++i) {
                    value.push_back((uint8_t)NextRandom(state));
                }
            }
            break;
        case 6:
            std::reverse(value.begin(), value.end());
            break;
        default:
            if (value.size() > 1u) {
                const size_t other =
                    (size_t)(NextRandom(state) % value.size());
                std::swap(value[position], value[other]);
            }
            break;
        }
    }
    return value;
}

bool MakeDirectory(const std::string& path)
{
    if (path.empty() || path == ".") {
        return true;
    }
#if defined(_WIN32)
    const int result = _mkdir(path.c_str());
#else
    const int result = mkdir(path.c_str(), 0777);
#endif
    return result == 0 || errno == EEXIST;
}

std::string SafeTargetName(const char* target_name)
{
    std::string safe = target_name ? target_name : "unknown";
    for (char& c : safe) {
        if (!((c >= 'a' && c <= 'z') ||
              (c >= 'A' && c <= 'Z') ||
              (c >= '0' && c <= '9') || c == '-' || c == '_'))
        {
            c = '_';
        }
    }
    return safe;
}

void PreserveFailure(
    const Options& options,
    const char* target_name,
    uint64_t mutation_index,
    const std::vector<uint8_t>& input,
    const std::string& failure)
{
    (void)MakeDirectory(options.ArtifactDirectory);
    char suffix[160];
    std::snprintf(
        suffix, sizeof(suffix),
        "%s-seed-%016" PRIx64 "-case-%" PRIu64,
        SafeTargetName(target_name).c_str(),
        options.Seed,
        mutation_index);
    const std::string base =
        options.ArtifactDirectory + "/" + suffix;
    std::ofstream binary((base + ".bin").c_str(), std::ios::binary);
    if (!input.empty()) {
        binary.write(
            reinterpret_cast<const char*>(input.data()),
            (std::streamsize)input.size());
    }
    std::ofstream metadata((base + ".txt").c_str());
    metadata << "target=" << target_name << "\n"
             << "seed=0x" << std::hex << options.Seed << std::dec << "\n"
             << "mutation=" << mutation_index << "\n"
             << "bytes=" << input.size() << "\n"
             << "failure=" << failure << "\n";
    std::fprintf(
        stderr,
        "FUZZ FAILURE target=%s seed=0x%016" PRIx64
        " mutation=%" PRIu64 " bytes=%zu artifact=%s.bin: %s\n",
        target_name,
        options.Seed,
        mutation_index,
        input.size(),
        base.c_str(),
        failure.c_str());
}

bool ExecuteCase(
    FuzzCaseFunction fuzz_case,
    const std::vector<uint8_t>& input,
    std::string& failure)
{
    if (input.size() > kMaxFuzzInputBytes) {
        failure = "input exceeds 64 MiB";
        return false;
    }
    try {
        return fuzz_case(input.data(), input.size(), failure);
    }
    catch (const std::bad_alloc&) {
        failure = "unexpected std::bad_alloc";
        return false;
    }
    catch (const std::exception& error) {
        failure = std::string("unexpected exception: ") + error.what();
        return false;
    }
    catch (...) {
        failure = "unexpected non-standard exception";
        return false;
    }
}

bool ParseOptions(
    int argc,
    char** argv,
    const char* default_corpus_manifest,
    Options& options)
{
    options.CorpusManifest = default_corpus_manifest ?
        default_corpus_manifest : "";
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        const auto require_value = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "%s requires a value\n", name);
                return nullptr;
            }
            return argv[++i];
        };
        if (argument == "--trace-current") {
            options.TraceCurrent = true;
            continue;
        }
        const char* value = nullptr;
        if (argument == "--corpus-manifest") {
            value = require_value("--corpus-manifest");
            if (!value) return false;
            options.CorpusManifest = value;
        }
        else if (argument == "--artifact-dir") {
            value = require_value("--artifact-dir");
            if (!value) return false;
            options.ArtifactDirectory = value;
        }
        else if (argument == "--replay") {
            value = require_value("--replay");
            if (!value) return false;
            options.ReplayPath = value;
        }
        else if (argument == "--seed") {
            value = require_value("--seed");
            if (!value || !ParseU64(value, options.Seed)) return false;
        }
        else if (argument == "--mutations") {
            value = require_value("--mutations");
            if (!value || !ParseU64(value, options.Mutations) ||
                options.Mutations == 0u)
            {
                return false;
            }
        }
        else if (argument == "--replay-index") {
            value = require_value("--replay-index");
            if (!value || !ParseU64(value, options.ReplayIndex)) return false;
            options.HaveReplayIndex = true;
        }
        else if (argument == "--max-seconds") {
            value = require_value("--max-seconds");
            if (!value || !ParseU32(value, options.MaximumSeconds) ||
                options.MaximumSeconds == 0u)
            {
                return false;
            }
        }
        else if (argument == "--duration-seconds") {
            value = require_value("--duration-seconds");
            if (!value || !ParseU32(value, options.DurationSeconds) ||
                options.DurationSeconds == 0u)
            {
                return false;
            }
        }
        else {
            std::fprintf(stderr, "unknown argument: %s\n", argument.c_str());
            return false;
        }
    }
    if (options.CorpusManifest.empty() ||
        (!options.ReplayPath.empty() && options.HaveReplayIndex))
    {
        return false;
    }
    return true;
}

} // namespace

Input::Input(const uint8_t* data, size_t size)
    : Data(data)
    , Size(size)
{
}

uint8_t Input::Byte()
{
    if (!Data || Offset >= Size) {
        ++Offset;
        return 0u;
    }
    return Data[Offset++];
}

uint8_t Input::U8() { return Byte(); }

uint16_t Input::U16()
{
    uint16_t value = Byte();
    value |= (uint16_t)Byte() << 8;
    return value;
}

uint32_t Input::U32()
{
    uint32_t value = U16();
    value |= (uint32_t)U16() << 16;
    return value;
}

uint64_t Input::U64()
{
    uint64_t value = U32();
    value |= (uint64_t)U32() << 32;
    return value;
}

bool Input::Bool() { return (Byte() & 1u) != 0u; }

size_t Input::Remaining() const
{
    return Offset < Size ? Size - Offset : 0u;
}

int RunDeterministicFuzzer(
    int argc,
    char** argv,
    const char* target_name,
    const char* default_corpus_manifest,
    FuzzCaseFunction fuzz_case)
{
    Options options;
    if (!ParseOptions(argc, argv, default_corpus_manifest, options)) {
        std::fprintf(stderr,
            "usage: %s [--corpus-manifest path] [--artifact-dir path] "
            "[--seed n] [--mutations n] [--max-seconds n] "
            "[--duration-seconds n] [--replay path|--replay-index n] "
            "[--trace-current]\n",
            argc > 0 ? argv[0] : "v2-fuzz");
        return 2;
    }

    std::vector<std::vector<uint8_t> > corpus;
    std::string failure;
    if (!LoadCorpus(options.CorpusManifest, corpus, failure)) {
        std::fprintf(stderr, "corpus error: %s\n", failure.c_str());
        return 2;
    }

    if (!options.ReplayPath.empty())
    {
        std::vector<uint8_t> replay;
        if (!ReadFile(
                options.ReplayPath, kMaxFuzzInputBytes, replay, failure))
        {
            std::fprintf(stderr, "replay error: %s\n", failure.c_str());
            return 2;
        }
        if (!ExecuteCase(fuzz_case, replay, failure)) {
            PreserveFailure(options, target_name, 0u, replay, failure);
            return 1;
        }
        return 0;
    }

    if (options.HaveReplayIndex)
    {
        const std::vector<uint8_t> replay = Mutate(
            corpus, options.Seed, options.ReplayIndex);
        if (!ExecuteCase(fuzz_case, replay, failure)) {
            PreserveFailure(
                options, target_name, options.ReplayIndex, replay, failure);
            return 1;
        }
        return 0;
    }

    const Clock::time_point start = Clock::now();
    for (size_t i = 0; i < corpus.size(); ++i)
    {
        if (options.TraceCurrent)
        {
            std::fprintf(
                stderr,
                "FUZZ CURRENT target=%s corpus=%zu\n",
                target_name,
                i);
            std::fflush(stderr);
        }
        if (!ExecuteCase(fuzz_case, corpus[i], failure)) {
            PreserveFailure(options, target_name, (uint64_t)i, corpus[i], failure);
            return 1;
        }
    }

    uint64_t completed = 0u;
    for (;;)
    {
        if (options.DurationSeconds == 0u) {
            if (completed >= options.Mutations) break;
        }
        else
        {
            const double elapsed = std::chrono::duration<double>(
                Clock::now() - start).count();
            if (elapsed >= options.DurationSeconds) break;
        }
        if (options.TraceCurrent) {
            std::fprintf(
                stderr,
                "FUZZ CURRENT target=%s seed=0x%016" PRIx64
                " mutation=%" PRIu64 "\n",
                target_name, options.Seed, completed);
            std::fflush(stderr);
        }
        const std::vector<uint8_t> mutation = Mutate(
            corpus, options.Seed, completed);
        if (!ExecuteCase(fuzz_case, mutation, failure)) {
            PreserveFailure(
                options, target_name, completed, mutation, failure);
            return 1;
        }
        ++completed;
    }

    const double elapsed = std::chrono::duration<double>(
        Clock::now() - start).count();
    if (options.DurationSeconds == 0u && elapsed > options.MaximumSeconds)
    {
        std::fprintf(
            stderr,
            "fuzz target %s exceeded %u seconds: %.3f seconds for %" PRIu64
            " mutations\n",
            target_name, options.MaximumSeconds, elapsed, completed);
        return 1;
    }
    std::printf(
        "v2 fuzz target=%s corpus=%zu mutations=%" PRIu64
        " seed=0x%016" PRIx64 " seconds=%.3f: PASS\n",
        target_name, corpus.size(), completed, options.Seed, elapsed);
    return 0;
}

void RunCoverageGuidedCaseOrAbort(
    const char* target_name,
    FuzzCaseFunction fuzz_case,
    const uint8_t* data,
    size_t size)
{
    if (size > kMaxFuzzInputBytes) {
        return;
    }
    const uint8_t empty = 0u;
    if (!data) {
        if (size != 0u) {
            std::abort();
        }
        data = &empty;
    }
    std::string failure;
    if (!ExecuteCase(
            fuzz_case,
            std::vector<uint8_t>(data, data + size),
            failure))
    {
        std::fprintf(
            stderr, "coverage fuzz failure target=%s bytes=%zu: %s\n",
            target_name, size, failure.c_str());
        std::abort();
    }
}

} // namespace fuzz
} // namespace wirehair_v2
