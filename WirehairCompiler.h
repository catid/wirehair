#pragma once

#if __cplusplus >= 201703L && defined(__has_cpp_attribute)
# if __has_cpp_attribute(fallthrough)
#  define WIREHAIR_FALLTHROUGH [[fallthrough]]
# endif
#endif
#ifndef WIREHAIR_FALLTHROUGH
# if defined(__GNUC__) && __GNUC__ >= 7
#  define WIREHAIR_FALLTHROUGH __attribute__((fallthrough))
# else
#  define WIREHAIR_FALLTHROUGH ((void)0)
# endif
#endif
