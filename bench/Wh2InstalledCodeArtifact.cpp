#include "Wh2InstalledCodeArtifact.h"
#include "Wh2FrozenTrace.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <memory>
#include <stdexcept>

#if defined(__linux__)
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace wirehair {
namespace wh2_benchmark {
namespace {

void Require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

#if defined(__linux__)
void RequireCanonical(const std::string& path)
{
    char* resolved;
    do { resolved = ::realpath(path.c_str(), nullptr); }
    while (!resolved && errno == EINTR);
    const std::unique_ptr<char, decltype(&std::free)> owned(resolved, &std::free);
    Require(resolved && path == resolved, "installed code path is not canonical");
}

struct Descriptor
{
    int value;
    explicit Descriptor(const std::string& path)
    {
        do { value = ::open(path.c_str(), O_RDONLY | O_NONBLOCK | O_NOFOLLOW | O_CLOEXEC); }
        while (value < 0 && errno == EINTR);
        Require(value >= 0, "cannot open installed code");
    }
    // On Linux an interrupted close has already released the descriptor.
    ~Descriptor() { ::close(value); }
    Descriptor(const Descriptor&) = delete;
    Descriptor& operator=(const Descriptor&) = delete;
};

struct stat NamedStat(const std::string& path)
{
    struct stat result = {};
    int status;
    do { status = ::lstat(path.c_str(), &result); }
    while (status < 0 && errno == EINTR);
    Require(status == 0, "cannot stat installed code path");
    return result;
}

struct stat DescriptorStat(int fd)
{
    struct stat result = {};
    int status;
    do { status = ::fstat(fd, &result); }
    while (status < 0 && errno == EINTR);
    Require(status == 0, "cannot stat installed code descriptor");
    return result;
}

bool SameMetadata(const struct stat& a, const struct stat& b)
{
    // Compare fields, not padding; atime alone can change because of this read.
    return a.st_dev == b.st_dev && a.st_ino == b.st_ino && a.st_mode == b.st_mode &&
        a.st_nlink == b.st_nlink && a.st_uid == b.st_uid && a.st_gid == b.st_gid &&
        a.st_rdev == b.st_rdev && a.st_size == b.st_size &&
        a.st_blksize == b.st_blksize && a.st_blocks == b.st_blocks &&
        a.st_mtim.tv_sec == b.st_mtim.tv_sec && a.st_mtim.tv_nsec == b.st_mtim.tv_nsec &&
        a.st_ctim.tv_sec == b.st_ctim.tv_sec && a.st_ctim.tv_nsec == b.st_ctim.tv_nsec;
}
#endif

} // namespace

void VerifyInstalledCodeArtifact(const std::string& canonicalPath,
    const std::string& expectedSha256, std::size_t maxBytes)
{
    Require(!canonicalPath.empty() && canonicalPath[0] == '/' &&
        canonicalPath.find('\0') == std::string::npos, "invalid installed code path");
    Require(expectedSha256.size() == 64 &&
        expectedSha256.find_first_not_of("0123456789abcdef") == std::string::npos,
        "invalid installed code SHA-256");
#if defined(__linux__)
    const std::size_t limit = std::min(maxBytes, std::size_t(8u * 1024u * 1024u));
    RequireCanonical(canonicalPath);
    const struct stat namedBefore = NamedStat(canonicalPath);
    Require(S_ISREG(namedBefore.st_mode) && namedBefore.st_nlink >= 1 &&
        namedBefore.st_size >= 0 && static_cast<uint64_t>(namedBefore.st_size) <= limit,
        "installed code type, link count, or size differs");
    const Descriptor file(canonicalPath);
    const struct stat before = DescriptorStat(file.value);
    Require(SameMetadata(namedBefore, before), "installed code path/descriptor differs");

    const std::size_t expectedBytes = static_cast<std::size_t>(before.st_size);
    std::string bytes;
    bytes.reserve(expectedBytes);
    char buffer[16384];
    for (;;) {
        const std::size_t remaining = expectedBytes - bytes.size();
        // One extra byte at the frozen length distinguishes EOF from growth.
        const std::size_t request = remaining ? std::min(remaining, sizeof(buffer)) : 1;
        ssize_t count;
        do { count = ::read(file.value, buffer, request); }
        while (count < 0 && errno == EINTR);
        Require(count >= 0, "cannot read installed code");
        if (count == 0) break;
        Require(static_cast<std::size_t>(count) <= request &&
            static_cast<std::size_t>(count) <= remaining, "installed code grew while reading");
        bytes.append(buffer, static_cast<std::size_t>(count));
    }
    Require(bytes.size() == expectedBytes, "installed code length changed while reading");
    const std::string actualSha256 = Sha256Hex(bytes);
    const struct stat after = DescriptorStat(file.value);
    RequireCanonical(canonicalPath);
    const struct stat namedAfter = NamedStat(canonicalPath);
    Require(SameMetadata(before, after) && SameMetadata(before, namedAfter),
        "installed code metadata or named identity changed while reading");
    Require(actualSha256 == expectedSha256, "installed code SHA-256 differs");
#else
    (void)maxBytes;
    throw std::runtime_error("installed code verification requires Linux");
#endif
}

} // namespace wh2_benchmark
} // namespace wirehair
