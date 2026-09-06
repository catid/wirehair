#ifndef WIREHAIR_BENCH_WH2_INSTALLED_CODE_ARTIFACT_H
#define WIREHAIR_BENCH_WH2_INSTALLED_CODE_ARTIFACT_H

#include <cstddef>
#include <string>

namespace wirehair {
namespace wh2_benchmark {

// Benchmark-only Linux installed-code check, not an evidence-file reader.
// Requires a canonical absolute nonsymlink path and lowercase SHA-256. Stable
// hard links are allowed; every other file property except atime must remain
// unchanged while reading. maxBytes can lower, never raise, the 8 MiB cap.
// Throws on invalid input, unsupported platforms, or any verification failure.
void VerifyInstalledCodeArtifact(
    const std::string& canonicalPath,
    const std::string& expectedSha256,
    std::size_t maxBytes = 8u * 1024u * 1024u);

} // namespace wh2_benchmark
} // namespace wirehair

#endif // WIREHAIR_BENCH_WH2_INSTALLED_CODE_ARTIFACT_H
