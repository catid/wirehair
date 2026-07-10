#include "roundtrip.h"

#if defined(_WIN32)
#define WIREHAIR_PLUGIN_EXPORT __declspec(dllexport)
#else
#define WIREHAIR_PLUGIN_EXPORT
#endif

WIREHAIR_PLUGIN_EXPORT int wirehair_plugin_round_trip(void)
{
    return wirehair_package_round_trip();
}
