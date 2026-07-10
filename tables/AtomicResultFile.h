#pragma once

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace wirehair_table {

class AtomicResultFile
{
public:
    explicit AtomicResultFile(const std::string& final_path)
        : FinalPath(final_path)
        , TemporaryPath(final_path.empty() ? "" : final_path + ".tmp")
    {
    }

    ~AtomicResultFile()
    {
        CloseDescriptor();
        RemoveOwnedTemporary();
    }

    bool Open(std::ostream& errors)
    {
        if (FinalPath.empty()) {
            return true;
        }
        if (PathExists(FinalPath)) {
            errors << "--results destination already exists: "
                << FinalPath << std::endl;
            return false;
        }
#ifdef _WIN32
        Descriptor = _open(
            TemporaryPath.c_str(),
            _O_WRONLY | _O_CREAT | _O_EXCL | _O_BINARY,
            _S_IREAD | _S_IWRITE);
#else
        int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_CLOEXEC
        flags |= O_CLOEXEC;
#endif
        Descriptor = open(TemporaryPath.c_str(), flags, 0666);
#endif
        if (Descriptor < 0) {
            errors << "unable to reserve --results temporary file "
                << TemporaryPath << ": " << LastError() << std::endl;
            return false;
        }
        if (!IdentityFromDescriptor(Descriptor, ReservedIdentity)) {
            errors << "unable to identify --results temporary file "
                << TemporaryPath << std::endl;
            CloseDescriptor();
            return false;
        }
        Opened = true;
        return true;
    }

    bool Enabled() const
    {
        return !FinalPath.empty();
    }

    std::ostream& Stream()
    {
        return Output;
    }

    bool Commit(std::ostream& errors)
    {
        if (!Enabled()) {
            return true;
        }
        if (!Opened || Descriptor < 0) {
            errors << "--results file was not opened" << std::endl;
            return false;
        }
        FileIdentity current;
        if (!IdentityFromPath(TemporaryPath, current) ||
            !SameIdentity(current, ReservedIdentity))
        {
            errors << "--results temporary file identity changed: "
                << TemporaryPath << std::endl;
            CloseDescriptor();
            return false;
        }
        if (!Output.good()) {
            errors << "failed buffering --results file: "
                << FinalPath << std::endl;
            CloseDescriptor();
            return false;
        }
        const std::string data = Output.str();
        if (!WriteAll(Descriptor, data.data(), data.size()) ||
            !FlushDescriptor(Descriptor))
        {
            errors << "failed writing --results file " << TemporaryPath
                << ": " << LastError() << std::endl;
            CloseDescriptor();
            return false;
        }
        if (!CloseDescriptor()) {
            errors << "failed closing --results file " << TemporaryPath
                << ": " << LastError() << std::endl;
            return false;
        }
        if (!IdentityFromPath(TemporaryPath, current) ||
            !SameIdentity(current, ReservedIdentity))
        {
            errors << "--results temporary file changed before publication: "
                << TemporaryPath << std::endl;
            return false;
        }

#ifdef _WIN32
        if (!MoveFileExA(
                TemporaryPath.c_str(), FinalPath.c_str(),
                MOVEFILE_WRITE_THROUGH))
        {
            errors << "unable to publish --results file " << FinalPath
                << ": Windows error " << GetLastError() << std::endl;
            return false;
        }
        Published = true;
#else
        if (link(TemporaryPath.c_str(), FinalPath.c_str()) != 0) {
            errors << "unable to publish --results file " << FinalPath
                << ": " << LastError() << std::endl;
            return false;
        }
        Published = true;
        FileIdentity published_identity;
        if (!IdentityFromPath(FinalPath, published_identity) ||
            !SameIdentity(published_identity, ReservedIdentity))
        {
            errors << "published --results identity mismatch: "
                << FinalPath << std::endl;
            RemoveOwnedFinal();
            Published = false;
            return false;
        }
        if (unlink(TemporaryPath.c_str()) != 0) {
            errors << "unable to remove published temporary file "
                << TemporaryPath << ": " << LastError() << std::endl;
            RemoveOwnedFinal();
            Published = false;
            return false;
        }
#endif
        Opened = false;
        return true;
    }

private:
    struct FileIdentity
    {
        uint64_t Device = 0;
        uint64_t Inode = 0;
        bool Regular = false;
        bool Valid = false;
    };

    static bool SameIdentity(const FileIdentity& a, const FileIdentity& b)
    {
        return a.Valid && b.Valid && a.Regular && b.Regular &&
            a.Device == b.Device && a.Inode == b.Inode;
    }

    static bool PathExists(const std::string& path)
    {
#ifdef _WIN32
        const DWORD attributes = GetFileAttributesA(path.c_str());
        if (attributes != INVALID_FILE_ATTRIBUTES) {
            return true;
        }
        const DWORD error = GetLastError();
        return error != ERROR_FILE_NOT_FOUND && error != ERROR_PATH_NOT_FOUND;
#else
        struct stat info;
        if (lstat(path.c_str(), &info) == 0) {
            return true;
        }
        return errno != ENOENT;
#endif
    }

    static bool IdentityFromDescriptor(int descriptor, FileIdentity& identity)
    {
#ifdef _WIN32
        struct _stat64 info;
        if (_fstat64(descriptor, &info) != 0) {
            return false;
        }
        identity.Device = (uint64_t)info.st_dev;
        identity.Inode = (uint64_t)info.st_ino;
        identity.Regular = (info.st_mode & _S_IFMT) == _S_IFREG;
#else
        struct stat info;
        if (fstat(descriptor, &info) != 0) {
            return false;
        }
        identity.Device = (uint64_t)info.st_dev;
        identity.Inode = (uint64_t)info.st_ino;
        identity.Regular = S_ISREG(info.st_mode);
#endif
        identity.Valid = true;
        return identity.Regular;
    }

    static bool IdentityFromPath(
        const std::string& path, FileIdentity& identity)
    {
#ifdef _WIN32
        struct _stat64 info;
        if (_stat64(path.c_str(), &info) != 0) {
            return false;
        }
        identity.Device = (uint64_t)info.st_dev;
        identity.Inode = (uint64_t)info.st_ino;
        identity.Regular = (info.st_mode & _S_IFMT) == _S_IFREG;
#else
        struct stat info;
        if (lstat(path.c_str(), &info) != 0) {
            return false;
        }
        identity.Device = (uint64_t)info.st_dev;
        identity.Inode = (uint64_t)info.st_ino;
        identity.Regular = S_ISREG(info.st_mode);
#endif
        identity.Valid = true;
        return identity.Regular;
    }

    static bool WriteAll(int descriptor, const char* data, size_t bytes)
    {
        size_t offset = 0;
        while (offset < bytes)
        {
#ifdef _WIN32
            const unsigned chunk = (unsigned)std::min<size_t>(
                bytes - offset, 1u << 30);
            const int written = _write(descriptor, data + offset, chunk);
#else
            const ssize_t written = write(descriptor, data + offset, bytes - offset);
#endif
            if (written < 0) {
#ifndef _WIN32
                if (errno == EINTR) {
                    continue;
                }
#endif
                return false;
            }
            if (written == 0) {
                return false;
            }
            offset += (size_t)written;
        }
        return true;
    }

    static bool FlushDescriptor(int descriptor)
    {
#ifdef _WIN32
        return _commit(descriptor) == 0;
#else
        return fsync(descriptor) == 0;
#endif
    }

    bool CloseDescriptor()
    {
        if (Descriptor < 0) {
            return true;
        }
#ifdef _WIN32
        const bool ok = _close(Descriptor) == 0;
#else
        const bool ok = close(Descriptor) == 0;
#endif
        Descriptor = -1;
        return ok;
    }

    void RemoveOwnedTemporary()
    {
        if (!Opened || TemporaryPath.empty()) {
            return;
        }
        FileIdentity current;
        if (IdentityFromPath(TemporaryPath, current) &&
            SameIdentity(current, ReservedIdentity))
        {
            std::remove(TemporaryPath.c_str());
        }
        Opened = false;
    }

    void RemoveOwnedFinal()
    {
        if (!Published) {
            return;
        }
        FileIdentity current;
        if (IdentityFromPath(FinalPath, current) &&
            SameIdentity(current, ReservedIdentity))
        {
            std::remove(FinalPath.c_str());
        }
    }

    static std::string LastError()
    {
#ifdef _WIN32
        return "CRT error " + std::to_string(errno);
#else
        return std::strerror(errno);
#endif
    }

    std::string FinalPath;
    std::string TemporaryPath;
    std::ostringstream Output;
    FileIdentity ReservedIdentity;
    int Descriptor = -1;
    bool Opened = false;
    bool Published = false;
};

} // namespace wirehair_table
