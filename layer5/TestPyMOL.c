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

#include"ObjectCGO.h"
#include"VFont.h"
#include"ObjectGadget.h"
#include"P.h"

#include"ObjectMap.h"
#include"Executive.h"

int TestPyMOL_00_00(CTestPyMOL *I);

int TestPyMOL_00_00(CTestPyMOL *I) 
{
  ObjectMap *obj;
  ObjectMapDesc _md,*md;
  ObjectMapState *ms =NULL;

  int a;

  md=&_md;

  md->mode = cObjectMap_OrthoMinMaxGrid;

  for(a=0;a<3;a++) {
    md->Grid[a] = 0.1F;
    md->MinCorner[a] = 0.0F;
    md->MaxCorner[a] = a+1.0F;
  }
  md->init_mode = -2;
  
  obj = ObjectMapNew();
  if(obj) {
    ms = ObjectMapNewStateFromDesc(obj,md,0);    
  }
  if(obj) {
    ObjectSetName((CObject*)obj,"00_00");
    ExecutiveManageObject((CObject*)obj,true);
  }
  return (obj!=NULL);

}

int TestPyMOLRun(CTestPyMOL *I,int group,int test)
{
  switch(group) {
  case 0: /* development tests */
    switch(test) {
    case 0:
      TestPyMOL_00_00(I); break;
      
    case 1: 
      PBlock();
      VFontLoad(1,0,0,true); 
      PUnblock();
      break;
    case 2: {
      CObject *obj = NULL;
      float pos[3] = {0.0,0.0,0.0};
      PBlock();
      obj = (CObject*)ObjectCGONewVFontTest("hello",pos);
      PUnblock();
      if(obj) {
        ObjectSetName(obj,"hello");
        ExecutiveManageObject(obj,true);
      }
    }
    break;
    case 3: {
      CObject *obj = NULL;
      obj = (CObject*)ObjectGadgetTest();
      if(obj)  {
        ObjectSetName(obj,"gadget");
        ExecutiveManageObject(obj,true);
      }
    }
    }
  }
  return 1;
}






