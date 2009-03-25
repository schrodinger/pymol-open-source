/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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
typedef int aliased_int __attribute__((may_alias));
typedef float aliased_float __attribute__((may_alias));
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

#ifndef _PYMOL_MINGW
#ifdef _WIN32
#ifndef _AFXDLL
	#ifndef _UNICODE
		#ifdef _DEBUG
			#pragma comment(lib, "nafxcwd.lib")
		#else
			#pragma comment(lib, "nafxcw.lib")
		#endif
	#else
		#ifdef _DEBUG
			#pragma comment(lib, "uafxcwd.lib")
		#else
			#pragma comment(lib, "uafxcw.lib")
		#endif
	#endif
#else
	#ifndef _UNICODE
		#ifdef _DEBUG
			#pragma comment(lib, "mfc42d.lib")
			#pragma comment(lib, "mfcs42d.lib")
		#else
			#pragma comment(lib, "mfc42.lib")
			#pragma comment(lib, "mfcs42.lib")
		#endif
	#else
		#ifdef _DEBUG
			#pragma comment(lib, "mfc42ud.lib")
			#pragma comment(lib, "mfcs42ud.lib")
		#else
			#pragma comment(lib, "mfc42u.lib")
			#pragma comment(lib, "mfcs42u.lib")
		#endif
	#endif
#endif

#ifdef _DLL
	#if !defined(_AFX_NO_DEBUG_CRT) && defined(_DEBUG)
		#pragma comment(lib, "msvcrtd.lib")
	#else
		#pragma comment(lib, "msvcrt.lib")
	#endif
#else
#ifdef _MT
	#if !defined(_AFX_NO_DEBUG_CRT) && defined(_DEBUG)
		#pragma comment(lib, "libcmtd.lib")
	#else
		#pragma comment(lib, "libcmt.lib")
	#endif
#else
	#if !defined(_AFX_NO_DEBUG_CRT) && defined(_DEBUG)
		#pragma comment(lib, "libcd.lib")
	#else
		#pragma comment(lib, "libc.lib")
	#endif
#endif
#endif

#pragma comment(lib, "kernel32.lib")
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "comdlg32.lib")
#pragma comment(lib, "winspool.lib")
#pragma comment(lib, "advapi32.lib")
#pragma comment(lib, "shell32.lib")
#pragma comment(lib, "comctl32.lib")

/* force inclusion of NOLIB.OBJ for /disallowlib directives */
#ifndef _PYMOL_PYOMM
#pragma comment(linker, "/include:__afxForceEXCLUDE")
#endif

/* force inclusion of DLLMODUL.OBJ for _USRDLL */
#ifdef _USRDLL
#pragma comment(linker, "/include:__afxForceUSRDLL")
#endif

/* force inclusion of STDAFX.OBJ for precompiled types */
#ifdef _AFXDLL
#pragma comment(linker, "/include:__afxForceSTDAFX")
#endif

#endif
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

#include "ov_types.h"

#endif


