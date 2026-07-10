#include "../tables/AtomicResultFile.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifndef _WIN32
#include <sys/stat.h>
#include <unistd.h>
#endif

static int Fail(const char* message)
{
    std::cerr << "AtomicResultFileTest: " << message << std::endl;
    return 1;
}

static std::string Read(const std::string& path)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    return std::string(
        (std::istreambuf_iterator<char>(input)),
        std::istreambuf_iterator<char>());
}

int main()
{
#ifdef _WIN32
    char temp_directory[MAX_PATH];
    char temporary_path[MAX_PATH];
    if (!GetTempPathA(MAX_PATH, temp_directory) ||
        !GetTempFileNameA(
            temp_directory, "whr", 0, temporary_path))
    {
        return Fail("temporary path creation failed");
    }
    DeleteFileA(temporary_path);
    const std::string root = std::string(temporary_path) + ".dir";
    if (!CreateDirectoryA(root.c_str(), nullptr)) {
        return Fail("CreateDirectory failed");
    }
#else
    char directory_template[] = "/tmp/wirehair-atomic-result.XXXXXX";
    const char* created = mkdtemp(directory_template);
    if (!created) {
        return Fail("mkdtemp failed");
    }
    const std::string root = created;
#endif
    const std::string final_path = root + "/results.tsv";
    const std::string temporary = final_path + ".tmp";
    std::ostringstream errors;

    {
        wirehair_table::AtomicResultFile output(final_path);
        if (!output.Open(errors)) return Fail("initial open failed");
        output.Stream() << "header\nvalue\n";
        if (!output.Commit(errors)) return Fail("initial commit failed");
        if (Read(final_path) != "header\nvalue\n") {
            return Fail("published content mismatch");
        }
    }
#ifndef _WIN32
    chmod(final_path.c_str(), 0);
#endif
    {
        wirehair_table::AtomicResultFile output(final_path);
        if (output.Open(errors)) return Fail("existing final was accepted");
#ifndef _WIN32
        chmod(final_path.c_str(), 0600);
#endif
        if (Read(final_path) != "header\nvalue\n") {
            return Fail("existing final was changed");
        }
    }

    std::remove(final_path.c_str());
    {
        std::ofstream sentinel(temporary.c_str(), std::ios::binary);
        sentinel << "temporary sentinel";
    }
#ifndef _WIN32
    chmod(temporary.c_str(), 0);
#endif
    {
        wirehair_table::AtomicResultFile output(final_path);
        if (output.Open(errors)) return Fail("existing temporary was accepted");
    }
#ifndef _WIN32
    chmod(temporary.c_str(), 0600);
#endif
    if (Read(temporary) != "temporary sentinel") {
        return Fail("existing temporary was changed");
    }
    std::remove(temporary.c_str());

    {
        wirehair_table::AtomicResultFile output(final_path);
        if (!output.Open(errors)) return Fail("concurrent-final open failed");
        output.Stream() << "candidate";
        std::ofstream concurrent(final_path.c_str(), std::ios::binary);
        concurrent << "concurrent";
        concurrent.close();
        if (output.Commit(errors)) return Fail("concurrent final was replaced");
        if (Read(final_path) != "concurrent") {
            return Fail("concurrent final was changed");
        }
    }
    if (Read(temporary) != "") {
        return Fail("owned temporary was not cleaned");
    }
    std::remove(final_path.c_str());

#ifndef _WIN32
    {
        wirehair_table::AtomicResultFile output(final_path);
        if (!output.Open(errors)) return Fail("temp-swap open failed");
        output.Stream() << "candidate";
        if (unlink(temporary.c_str()) != 0) return Fail("temp unlink failed");
        std::ofstream attacker(temporary.c_str(), std::ios::binary);
        attacker << "attacker";
        attacker.close();
        if (output.Commit(errors)) return Fail("swapped temporary was accepted");
        if (Read(temporary) != "attacker") {
            return Fail("swapped temporary was removed");
        }
        std::remove(temporary.c_str());
    }
#endif

#ifdef _WIN32
    RemoveDirectoryA(root.c_str());
#else
    rmdir(root.c_str());
#endif
    return 0;
}
