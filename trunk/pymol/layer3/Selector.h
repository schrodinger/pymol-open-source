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
#ifndef _H_Selector
#define _H_Selector

#include"ObjectMolecule.h"

void SelectorInit(void);
int *SelectorSelect(char *sele);
void SelectorCreate(char *name,char *sele,ObjectMolecule *obj);
void SelectorToggle(int rep,char *name);
void SelectorCylinder(char *sele,char *onoff);
int SelectorUpdateTable(void);
int SelectorIndexByName(char *sele);
int SelectorMatch(int ref,int sele);
int SelectorNext(int ref);
void SelectorFree(void);
void SelectorDelete(char *sele);
void SelectorFreeTmp(char *name);
void SelectorGetTmp(char *input,char *store);
int SelectorGetPDB(char **charVLA,int sele,int state,int conectFlag);
float SelectorSumVDWOverlap(int sele1,int state1,int sele2,int state2);
int SelectorGetSeleNCSet(int sele);

#endif
