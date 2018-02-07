/*
 * strcasecmp, strncasecmp
 */

#pragma once

#ifdef _WIN32
  #include <string.h>
  #define strcasecmp(s1, s2) _stricmp(s1, s2)
  #define strncasecmp(s1, s2, n) _strnicmp(s1, s2, n)
#else
  #include <strings.h>
#endif
