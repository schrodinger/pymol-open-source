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
#ifndef _H_Base
#define _H_Base

#include"os_limits.h"
#include"os_types.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef NULL
#define NULL ((void*)0)
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif


#ifndef uchar
#define uchar unsigned char
#endif

#define MAX_VDW 2.5  /* this has to go */


#ifndef MAXFLOAT
#define MAXFLOAT FLT_MAX
#endif

#ifndef R_SMALL4
#define R_SMALL4 0.0001
#endif

typedef struct { 
  void *ptr;
  int index; /* NOTE: that first record contains the list count...not pick info */
} Pickable;

#define MAXLINELEN 1024

#ifndef _PYMOL_NO_XRAY
#define _PYMOL_XRAY
#endif

#endif
