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
#ifndef _H_Isosurf
#define _H_Isosurf

#include"Map.h"
#include"MemoryDebug.h"
#include"Crystal.h"

typedef struct {
  int dimensions[3];
  float *points;
  float *data;
} Isofield;

#define F3(Name,P1,P2,P3,dims) \
(*((Name)+((P1)+((dims[0])*((P2)+((dims[1])*(P3)))))))

#define F3Ptr(Name,P1,P2,P3,dims) \
((Name)+((P1)+((dims[0])*((P2)+((dims[1])*(P3))))))

#define F4(Name,P1,P2,P3,P4,dims) \
(*((Name)+((P1)+((dims[0])*((P2)+((dims[1])*((P3)+((dims[2])*(P4)))))))))

#define F4Ptr(Name,P1,P2,P3,P4,dims) \
((Name)+((P1)+((dims[0])*((P2)+((dims[1])*((P3)+((dims[2])*(P4))))))))

Isofield *IsosurfFieldAlloc(int *dims);
void IsosurfFieldFree(Isofield *field);

int	IsosurfVolume(Isofield *field,float level,int **num,float **vert,int *range,int mode);
void IsosurfGetRange(Isofield *field,CCrystal *cryst,float *mn,float *mx,int *range);

int	IsosurfInit(void);

#endif
