// Every direct shared consumer checks its actual runtime provider before main.
// No static codec implementation or private GF arithmetic is linked here.
#include <wirehair/wirehair.h>
#include <dlfcn.h>
#include <limits.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <initializer_list>

#ifndef WH2_SHARED_EXPECTED_PROVIDER
#error "An exact shared-library provider is required"
#endif

namespace {
struct ProviderCheck {
    ProviderCheck()
    {
        char expected[PATH_MAX];
        if (!realpath(WH2_SHARED_EXPECTED_PROVIDER, expected)) Fail();
        for (auto pointer : {reinterpret_cast<void*>(wirehair_init_),
                             reinterpret_cast<void*>(wirehair_v2_encoder_create),
                             reinterpret_cast<void*>(wirehair_v2_encoder_create_with_options),
                             reinterpret_cast<void*>(wirehair_v2_decoder_create),
                             reinterpret_cast<void*>(wirehair_v2_free)}) {
            Dl_info info = {};
            char observed[PATH_MAX];
            if (!dladdr(pointer, &info) || !info.dli_fname ||
                !realpath(info.dli_fname, observed) || std::strcmp(expected, observed)) Fail();
        }
    }
    static void Fail()
    {
        std::fputs("Shared test resolved the wrong Wirehair provider\n", stderr);
        std::exit(1);
    }
};
ProviderCheck provider_check;
}
