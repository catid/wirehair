#include <stdio.h>
#include <string.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

typedef int (*round_trip_fn)(void);

int main(int argc, char** argv)
{
    round_trip_fn round_trip;
    int result;
    if (argc != 2) {
        fprintf(stderr, "usage: package_plugin_host PLUGIN\n");
        return 2;
    }

#if defined(_WIN32)
    {
        HMODULE plugin = LoadLibraryA(argv[1]);
        FARPROC symbol;
        if (!plugin) {
            fprintf(stderr, "LoadLibrary failed: %lu\n",
                (unsigned long)GetLastError());
            return 3;
        }
        symbol = GetProcAddress(plugin, "wirehair_plugin_round_trip");
        if (!symbol) {
            fprintf(stderr, "GetProcAddress failed: %lu\n",
                (unsigned long)GetLastError());
            FreeLibrary(plugin);
            return 4;
        }
        if (sizeof(symbol) != sizeof(round_trip)) {
            fprintf(stderr, "function pointer size mismatch\n");
            FreeLibrary(plugin);
            return 5;
        }
        memcpy(&round_trip, &symbol, sizeof(round_trip));
        result = round_trip();
        FreeLibrary(plugin);
    }
#else
    {
        void* plugin = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
        const char* error;
        if (!plugin) {
            fprintf(stderr, "dlopen failed: %s\n", dlerror());
            return 3;
        }
        dlerror();
        *(void**)(&round_trip) = dlsym(plugin, "wirehair_plugin_round_trip");
        error = dlerror();
        if (error) {
            fprintf(stderr, "dlsym failed: %s\n", error);
            dlclose(plugin);
            return 4;
        }
        result = round_trip();
        dlclose(plugin);
    }
#endif
    return result;
}
