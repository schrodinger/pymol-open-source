

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_os_predef
#define _H_os_predef

/* Macros used by Fortify source in GCC 4.1.x are incompatible with
   PyMOL's Feedback system... */

#ifdef _FORTIFY_SOURCE
#undef _FORTIFY_SOURCE
#endif


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */

#ifdef WIN32
#define PATH_SEP "\\"
#else
#define PATH_SEP "/"
#endif


/* commercial product */

#ifdef PYMOL_COMM
#ifndef _PYMOL_IP_SPLASH
#define _PYMOL_IP_SPLASH
#endif
#ifndef _PYMOL_IP_EXTRAS
#define _PYMOL_IP_EXTRAS
#endif
#endif


/* collaboration product (placarded) */

#ifdef PYMOL_COLL
#ifndef _PYMOL_IP_SPLASH
#define _PYMOL_IP_SPLASH
#endif
#ifndef _PYMOL_IP_EXTRAS
#define _PYMOL_IP_EXTRAS
#endif
#endif


/* educational product (placarded) */

#ifdef PYMOL_EDU
#ifndef _PYMOL_IP_SPLASH
#define _PYMOL_IP_SPLASH
#endif
#ifndef _PYMOL_IP_EXTRAS
#define _PYMOL_IP_EXTRAS
#endif
#endif


/* evaluation product (placarded) */

#ifdef PYMOL_EVAL
#ifndef _PYMOL_IP_SPLASH
#define _PYMOL_IP_SPLASH
#endif
#endif


/* END PROPRIETARY CODE SEGMENT */

#include <stddef.h>

#if defined(_MSC_VER)
// conversion from '...' to '...', possible loss of data
#pragma warning (disable:4244)

// conversion from 'size_t' to '...', possible loss of data
#pragma warning (disable:4267)

// truncation from 'double' to 'float'
#pragma warning (disable:4305)

// forcing value to bool (performance warning)
#pragma warning (disable:4800)

// '_snprintf', 'sscanf', 'sprintf', 'strcpy', 'strncpy': This function or
// variable may be unsafe. To disable deprecation, use _CRT_SECURE_NO_WARNINGS.
#pragma warning (disable:4996)

// TODO: remove, use std::snprintf instead (C++11)
#if !defined(snprintf) && (_MSC_VER < 1900)
#define snprintf sprintf_s
#endif
#endif

#include "ov_types.h"

// alternative to std::swap if references are not allowed (e.g. bit fields)
#define SWAP_NOREF(a, b) {auto _t=(a);(a)=(b);(b)=_t;}

#endif
