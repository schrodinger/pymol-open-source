

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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


/* Macros used by Fortify source in GCC 4.1.x are incompatible with
   PyMOL's Feedback system... */

#ifdef _FORTIFY_SOURCE
#undef _FORTIFY_SOURCE
#endif


/* Alias-able typedefs */

#ifdef __GNUC__
#if __GNUC__ > 3
typedef int aliased_int __attribute__ ((may_alias));
typedef float aliased_float __attribute__ ((may_alias));
#else
typedef int aliased_int;
typedef float aliased_float;
#endif
#else
typedef int aliased_int;
typedef float aliased_float;
#endif


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */

#ifdef WIN32
#define __inline__ __inline
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

#ifdef __linux__
#include <malloc.h>
#else
#include <stddef.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define fmax max
#define fmin min
#pragma warning (disable:4996)
#define snprintf sprintf_s
#endif

#include "ov_types.h"

#endif
