/*
 * strcasecmp, strncasecmp
 */

#pragma once

#if defined(_WIN32) && defined(_MSC_VER)
  #include <string.h>
  #define strcasecmp(s1, s2) _stricmp(s1, s2)
  #define strncasecmp(s1, s2, n) _strnicmp(s1, s2, n)
#else
  #include <strings.h>
#endif
