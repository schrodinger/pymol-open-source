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
#ifndef _H_Executive
#define _H_Executive

#include"Object.h"
#include"Ortho.h"
#include"Word.h"

void ExecutiveDistance(char *sele1,char *sele2);
void ExecutiveAlter(char *s1,char *expr);
void ExecutiveColor(char *name,char *color,int flags);
void ExecutiveInit(void);
void ExecutiveFree(void);
void ExecutiveManageObject(struct Object *obj);
void ExecutiveUpdateObjectSelection(struct Object *obj);
void ExecutiveManageSelection(char *name);
Block *ExecutiveGetBlock(void);
Object *ExecutiveFindObjectByName(char *name);
int ExecutiveIterateObject(Object **obj,void **hidden);
void ExecutiveDelete(char *name);
void ExecutiveSetControlsOff(char *name);
void ExecutiveSetSetting(char *sname,char *value);
void ExecutiveRay(void);
void ExecutiveFit(char *sele1,char *sele2);
void ExecutiveReset(int cmd);
void ExecutiveDrawNow(void);
void ExecutiveSetAllVisib(int state);
void ExecutiveSetRepVisib(char *name,int rep,int state);
void ExecutiveSetObjVisib(char *name,int state);
void ExecutiveCenter(char *name,int preserve);
void ExecutiveWindowZoom(char *name);
int ExecutiveGetMoment(char *name,Matrix33d mi);
void ExecutiveOrient(char *sele,Matrix33d mi);
char *ExecutiveSeleToPDBStr(char *s1,int state);
void ExecutiveStereo(int flag);
void ExecutiveCopy(char *src,char *dst);

#endif



