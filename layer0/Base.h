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

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define MAX_VDW 2.5

typedef struct { 
  void *ptr;
  int index; /* NOTE: that first record contains the list count...not pick info */
} Pickable;


#endif
