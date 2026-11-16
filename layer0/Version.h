#ifndef _PyMOL_VERSION
#define _PyMOL_VERSION "4.12.9"
#endif

/* for session file compatibility */

#ifndef _PyMOL_VERSION_int
// X.Y.Z   -> XYYYZZZ
#define _PyMOL_VERSION_int 4012009
// Note: There should have never been a "double" version, it's
// the least useful variant to specify a version.
#define _PyMOL_VERSION_double (_PyMOL_VERSION_int / 1000000.)
#endif
