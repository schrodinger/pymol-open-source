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

/* Module for internal C-level PyMOL tests...*/

#include"os_predef.h"
#include"os_std.h"
#include"Base.h"

#include"MemoryDebug.h"
#include"Feedback.h"
#include"TestPyMOL.h"

#include"ObjectMap.h"
#include"Executive.h"

int TestPyMOL_00_00(CTestPyMOL *I);

int TestPyMOL_00_00(CTestPyMOL *I) 
{
  ObjectMap *obj;
  ObjectMapDesc _md,*md;
  int a;

  md=&_md;

  md->mode = cObjectMap_OrthoMinMaxGrid;

  for(a=0;a<3;a++) {
    md->Grid[a] = 0.1F;
    md->MinCorner[a] = 0.0F;
    md->MaxCorner[a] = a+1.0F;
  }
  md->init_mode = -2;
  
  obj = ObjectMapNewFromDesc(md);
  if(obj) {
    ObjectSetName((CObject*)obj,"00_00");
    ExecutiveManageObject((CObject*)obj);
  }
  return (obj!=NULL);

}

int TestPyMOLRun(CTestPyMOL *I,int group,int test)
{
  switch(group) {
  case 0: /* development tests */
    switch(test) {
    case 0: TestPyMOL_00_00(I); break;
    }
  }
  return 1;
}
