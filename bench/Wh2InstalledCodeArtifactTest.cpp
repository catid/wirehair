// Unit-only filesystem tests. Never invoke the included launcher or any codec.
#define main Wh2InstalledCodeUnusedOffsetMain
#define offset_r0 installed_code_evidence_private
#include "Wh2PublicOffsetR0.cpp"
#undef offset_r0
#undef main
#include "Wh2InstalledCodeArtifact.h"
#include <cstdarg>
#include <functional>

extern "C" {
int __real_open(const char*, int, ...);
int __real_close(int);
int __real_fstat(int, struct stat*);
int __real_lstat(const char*, struct stat*);
ssize_t __real_read(int, void*, size_t);
ssize_t __real___read_chk(int, void*, size_t, size_t);
char* __real_realpath(const char*, char*);
}

namespace {
using wirehair::wh2_benchmark::Sha256Hex;
using wirehair::wh2_benchmark::VerifyInstalledCodeArtifact;
void Check(bool ok, const char* text)
{ if (!ok) throw std::runtime_error(text); }

enum class Fault {
    None, OpenEintr, StatEintr, PathEintr, RealpathEintr, ReadEintr,
    ShortRead, ReadError, EarlyEof, Growth, CorruptRead, InitialIdentity, Metadata, NamedIdentity,
    CanonicalChange, AtimeOnly
};
struct Observation {
    std::string path;
    Fault fault = Fault::None;
    int field = 0, fd = -1, flags = 0;
    unsigned opens = 0, closes = 0, stats = 0, paths = 0, resolutions = 0;
    unsigned reads = 0, hits = 0;
} observed;

bool Target(const char* path)
{ return !observed.path.empty() && observed.path == path; }
bool Interrupt(Fault fault, unsigned call)
{
    if (observed.fault != fault || call != 1) return false;
    ++observed.hits; errno = EINTR; return true;
}
void ChangeMetadata(struct stat& value, int field)
{
    switch (field) {
    case 0: ++value.st_dev; break;
    case 1: ++value.st_ino; break;
    case 2: value.st_mode ^= S_IXUSR; break;
    case 3: ++value.st_nlink; break;
    case 4: ++value.st_uid; break;
    case 5: ++value.st_gid; break;
    case 6: ++value.st_rdev; break;
    case 7: ++value.st_size; break;
    case 8: ++value.st_blksize; break;
    case 9: ++value.st_blocks; break;
    case 10: ++value.st_mtim.tv_sec; break;
    case 11: value.st_mtim.tv_nsec ^= 1; break;
    case 12: ++value.st_ctim.tv_sec; break;
    case 13: value.st_ctim.tv_nsec ^= 1; break;
    default: Check(false, "unknown metadata field");
    }
}

struct Watch {
    Watch(const std::string& path, Fault fault = Fault::None, int field = 0)
    {
        Check(observed.path.empty(), "nested syscall observation");
        observed = Observation(); observed.path = path;
        observed.fault = fault; observed.field = field;
    }
    ~Watch()
    {
        if (observed.fd >= 0) __real_close(observed.fd);
        observed = Observation();
    }
    void Finished(bool opened, bool injected = false) const
    {
        Check(observed.fd == -1 && observed.opens == observed.closes,
            "descriptor not closed exactly once");
        Check(!opened || observed.opens == 1, "expected one successful open");
        Check(!injected || observed.hits > 0, "injection was not exercised");
    }
};

struct Temp {
    std::string directory;
    std::vector<std::string> entries;
    Temp()
    {
        char name[] = "/tmp/wh2-installed-code-test.XXXXXX";
        Check(::mkdtemp(name) != nullptr, "mkdtemp failed");
        directory = name;
    }
    ~Temp()
    {
        for (auto it = entries.rbegin(); it != entries.rend(); ++it) ::unlink(it->c_str());
        ::rmdir(directory.c_str());
    }
    std::string Name(const char* name)
    { entries.push_back(directory + '/' + name); return entries.back(); }
    std::string File(const char* name, const std::string& bytes)
    {
        const std::string path = Name(name);
        const int fd = __real_open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
        Check(fd >= 0, "fixture open failed");
        const ssize_t size = ::write(fd, bytes.data(), bytes.size());
        const int closed = __real_close(fd);
        Check(size >= 0 && static_cast<size_t>(size) == bytes.size() && closed == 0,
            "fixture write failed");
        return path;
    }
};
unsigned tests = 0;
void Case(const char* name, const std::function<void()>& body)
{
    try { body(); ++tests; }
    catch (const std::exception& error) {
        throw std::runtime_error(std::string(name) + ": " + error.what());
    }
}
void Reject(const std::function<void()>& body)
{
    bool rejected = false;
    try { body(); } catch (const std::runtime_error&) { rejected = true; }
    Check(rejected, "invalid artifact was accepted");
}

void Units()
{
    Temp temp;
    const std::string bytes = "installed code\n";
    const std::string path = temp.File("code", bytes), hash = Sha256Hex(bytes);
    const std::string hard = temp.Name("hardlink");
    Check(::link(path.c_str(), hard.c_str()) == 0, "fixture hardlink failed");
    Case("installed hardlink and exact cap", [&] {
        struct stat value = {};
        Check(__real_lstat(path.c_str(), &value) == 0 && value.st_nlink == 2,
            "fixture is not hardlinked");
        Watch watch(path);
        VerifyInstalledCodeArtifact(path, hash, bytes.size());
        const int needed = O_NONBLOCK | O_NOFOLLOW | O_CLOEXEC;
        Check((observed.flags & O_ACCMODE) == O_RDONLY &&
            (observed.flags & needed) == needed, "unsafe installed-code open flags");
        watch.Finished(true);
        VerifyInstalledCodeArtifact(hard, hash);
    });
    Case("old evidence helper stays strict", [&] {
        Reject([&] { installed_code_evidence_private::FileHash(path); });
        const std::string single = temp.File("single", bytes);
        Check(installed_code_evidence_private::FileHash(single) == hash,
            "strict evidence positive control failed");
    });
    Case("wrong content hash", [&] {
        Watch watch(path);
        Reject([&] { VerifyInstalledCodeArtifact(path, std::string(64, '0')); });
        watch.Finished(true);
    });
    Case("hash and path syntax", [&] {
        for (const std::string& bad : {std::string(), std::string(63, 'a'),
            std::string(64, 'G'), hash + std::string(1, '\0')})
            Reject([&] { VerifyInstalledCodeArtifact(path, bad); });
        for (const std::string& bad : {std::string(), std::string("relative-code"),
            temp.directory + "/./code", temp.directory + "/../" +
                temp.directory.substr(temp.directory.rfind('/') + 1) + "/code",
            path + std::string("\0suffix", 7)})
            Reject([&] { VerifyInstalledCodeArtifact(bad, hash); });
        Reject([&] { VerifyInstalledCodeArtifact(temp.directory + "/absent", hash); });
    });
    Case("symlink, directory, and FIFO", [&] {
        const std::string link = temp.Name("symlink"), fifo = temp.Name("fifo");
        Check(::symlink(path.c_str(), link.c_str()) == 0 &&
            ::mkfifo(fifo.c_str(), 0600) == 0, "special fixture creation failed");
        Reject([&] { VerifyInstalledCodeArtifact(link, hash); });
        Reject([&] { VerifyInstalledCodeArtifact(temp.directory, hash); });
        // Reject FIFO before opening; the regular-file case checks O_NONBLOCK.
        Reject([&] { VerifyInstalledCodeArtifact(fifo, hash); });
    });
    Case("byte caps and caller cannot widen hard cap", [&] {
        Reject([&] { VerifyInstalledCodeArtifact(path, hash, bytes.size() - 1); });
        Reject([&] { VerifyInstalledCodeArtifact(path, hash, 0); });
        const std::string large = temp.File("oversize", "");
        const int fd = __real_open(large.c_str(), O_WRONLY | O_CLOEXEC);
        Check(fd >= 0, "oversize fixture open failed");
        const int resized = ::ftruncate(fd, 8 * 1024 * 1024 + 1);
        const int closed = __real_close(fd);
        Check(resized == 0 && closed == 0, "oversize fixture resize failed");
        Watch watch(large);
        Reject([&] { VerifyInstalledCodeArtifact(large, hash, 16u * 1024u * 1024u); });
        Check(observed.reads == 0, "oversize artifact read before rejection");
        watch.Finished(false);
    });
    Case("empty regular artifact", [&] {
        VerifyInstalledCodeArtifact(temp.File("empty", ""), Sha256Hex(""));
    });
    Case("multi-chunk regular artifact", [&] {
        const std::string content(32769, 'x');
        const std::string multiple = temp.File("multiple", content);
        Watch watch(multiple);
        VerifyInstalledCodeArtifact(multiple, Sha256Hex(content));
        if (observed.reads < 4)
            throw std::runtime_error("too few chunk/EOF observations: " + std::to_string(observed.reads));
        watch.Finished(true);
    });
    for (Fault fault : {Fault::OpenEintr, Fault::StatEintr, Fault::PathEintr,
        Fault::RealpathEintr, Fault::ReadEintr, Fault::ShortRead, Fault::AtimeOnly})
        Case("benign syscall interruption/short read/atime", [&] {
            Watch watch(path, fault);
            VerifyInstalledCodeArtifact(path, hash);
            watch.Finished(true, true);
        });
    for (Fault fault : {Fault::ReadError, Fault::EarlyEof, Fault::Growth, Fault::CorruptRead, Fault::InitialIdentity,
        Fault::NamedIdentity, Fault::CanonicalChange})
        Case("read or named identity changes", [&] {
            Watch watch(path, fault);
            Reject([&] { VerifyInstalledCodeArtifact(path, hash); });
            watch.Finished(true, true);
        });
    for (int field = 0; field < 14; ++field)
        Case("stable descriptor metadata", [&] {
            Watch watch(path, Fault::Metadata, field);
            Reject([&] { VerifyInstalledCodeArtifact(path, hash); });
            watch.Finished(true, true);
        });
}
} // namespace

extern "C" int __wrap_open(const char* path, int flags, ...)
{
    mode_t mode = 0;
    if ((flags & O_CREAT) || (flags & O_TMPFILE) == O_TMPFILE) {
        va_list args; va_start(args, flags); mode = va_arg(args, mode_t); va_end(args);
    }
    if (Target(path)) {
        observed.flags = flags;
        if (observed.fault == Fault::OpenEintr && !observed.hits) {
            ++observed.hits; errno = EINTR; return -1;
        }
    }
    const int fd = __real_open(path, flags, mode);
    if (Target(path) && fd >= 0) { observed.fd = fd; ++observed.opens; }
    return fd;
}
extern "C" int __wrap_close(int fd)
{
    if (fd == observed.fd) { ++observed.closes; observed.fd = -1; }
    return __real_close(fd);
}
extern "C" int __wrap_fstat(int fd, struct stat* value)
{
    const bool target = fd == observed.fd;
    if (target && Interrupt(Fault::StatEintr, ++observed.stats)) return -1;
    const int result = __real_fstat(fd, value);
    if (target && result == 0 && observed.stats == 1 && observed.fault == Fault::InitialIdentity) {
        ++value->st_ino; ++observed.hits;
    }
    if (target && result == 0 && observed.stats >= 2) {
        if (observed.fault == Fault::Metadata) {
            ChangeMetadata(*value, observed.field); ++observed.hits;
        } else if (observed.fault == Fault::AtimeOnly) {
            ++value->st_atim.tv_sec; ++observed.hits;
        }
    }
    return result;
}
extern "C" int __wrap_lstat(const char* path, struct stat* value)
{
    const bool target = Target(path);
    if (target && Interrupt(Fault::PathEintr, ++observed.paths)) return -1;
    const int result = __real_lstat(path, value);
    if (target && result == 0 && observed.paths >= 2 && observed.fault == Fault::NamedIdentity) {
        ++value->st_ino; ++observed.hits;
    }
    return result;
}
extern "C" char* __wrap_realpath(const char* path, char* result)
{
    const bool target = Target(path);
    if (target && Interrupt(Fault::RealpathEintr, ++observed.resolutions)) return nullptr;
    char* resolved = __real_realpath(path, result);
    if (target && resolved && observed.resolutions >= 2 && observed.fault == Fault::CanonicalChange) {
        resolved[1] = resolved[1] == 'x' ? 'y' : 'x'; ++observed.hits;
    }
    return resolved;
}
extern "C" ssize_t __wrap_read(int fd, void* bytes, size_t count)
{
    if (fd != observed.fd) return __real_read(fd, bytes, count);
    if (Interrupt(Fault::ReadEintr, ++observed.reads)) return -1;
    if (observed.fault == Fault::ReadError) { ++observed.hits; errno = EIO; return -1; }
    if (observed.fault == Fault::EarlyEof) { ++observed.hits; return 0; }
    if (observed.fault == Fault::ShortRead && count > 1) { count = 1; ++observed.hits; }
    const ssize_t size = __real_read(fd, bytes, count);
    if (size == 0 && observed.fault == Fault::Growth) {
        Check(count == 1, "growth did not target EOF probe");
        static_cast<unsigned char*>(bytes)[0] = 0; ++observed.hits; return 1;
    }
    if (size > 0 && observed.fault == Fault::CorruptRead && !observed.hits) {
        static_cast<unsigned char*>(bytes)[0] ^= 1; ++observed.hits;
    }
    return size;
}

// Fortified O1/sanitizer builds can call __read_chk instead of read. Preserve
// its overflow guard and observe the same valid reads in both compiler forms.
extern "C" ssize_t __wrap___read_chk(int fd, void* bytes, size_t count, size_t capacity)
{
    if (count > capacity) return __real___read_chk(fd, bytes, count, capacity);
    return __wrap_read(fd, bytes, count);
}

int main(int argc, char** argv)
{
    try {
        if (argc == 2 && std::string(argv[1]) == "--libc-preflight") {
            VerifyInstalledCodeArtifact("/usr/lib/x86_64-linux-gnu/libc.so.6",
                "8db37cf3f2169f59a0f07ef1fea308c35656668c64c8ff294e1860f4121eb161");
            std::cout << "installed-code pinned-libc preflight PASS (no codec calls)\n";
        } else {
            Check(argc == 1, "usage: installed-code-test [--libc-preflight]");
            Units(); std::cout << "installed-code units PASS " << tests << " cases\n";
        }
        std::cout.flush(); Check(bool(std::cout), "cannot publish test result");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "installed-code test: " << error.what() << '\n'; return 1;
    }
}
