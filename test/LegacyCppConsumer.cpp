#define main legacy_c_consumer_main
#include "LegacyCConsumer.c"
#undef main

int main()
{
    return legacy_c_consumer_main();
}
