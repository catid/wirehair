#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace wirehair_v2 {
namespace fuzz {

static const size_t kMaxFuzzInputBytes =
    (size_t)64u * 1024u * 1024u;
static const size_t kMaxDeterministicMutationBytes = 64u * 1024u;
static const size_t kMaxCorpusArtifactBytes =
    (size_t)5u * 1024u * 1024u;

typedef bool (*FuzzCaseFunction)(
    const uint8_t* data,
    size_t size,
    std::string& failure);

class Input
{
public:
    Input(const uint8_t* data, size_t size);

    uint8_t U8();
    uint16_t U16();
    uint32_t U32();
    uint64_t U64();
    bool Bool();
    size_t Remaining() const;

private:
    uint8_t Byte();

    const uint8_t* Data;
    size_t Size;
    size_t Offset = 0u;
};

/** Deterministic corpus/mutation CLI used by ordinary and sanitizer CTest. */
int RunDeterministicFuzzer(
    int argc,
    char** argv,
    const char* target_name,
    const char* default_corpus_manifest,
    FuzzCaseFunction fuzz_case);

/** Coverage-guided entry points call this and abort on false. */
void RunCoverageGuidedCaseOrAbort(
    const char* target_name,
    FuzzCaseFunction fuzz_case,
    const uint8_t* data,
    size_t size);

} // namespace fuzz
} // namespace wirehair_v2
