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
/* Example of embedding Python in another program */

#include<Python.h>

#include"os_std.h"

#include"MemoryDebug.h"
#include"Base.h"
#include"Err.h"
#include"Export.h"
#include"RepDot.h"
#include"Executive.h"
#include"ObjectMolecule.h"

void ExportDotsObjFree(ExportDotsObj *rec);

/*--------------------------------------------------------------------- */

ExportDotsObj *ExportDots(char *name,int csIndex)
{
  Object *obj;
  ObjectMolecule *objMol;
  RepDot *rep;
  CoordSet *cs;
  ExportDotsObj *result = NULL;
  int ok = true;

  obj=ExecutiveFindObjectByName(name);
  if(!obj) 
	 ok=ErrMessage("ExportDots","Not a valid object name.");
  else if(obj->type!=cObjectMolecule)
	 ok=ErrMessage("ExportDots","Not molecule object.");

  if(ok) {
	 ExecutiveSetRepVisib(name,cRepDot,1); 
	 objMol = (ObjectMolecule*)obj;
	 cs = ObjectMoleculeGetCoordSet(objMol,csIndex);
	 if(!cs)
		ok=ErrMessage("ExportDots","Invalid coordinate set number.");
  }

  if(ok) {
	 rep = (RepDot*)RepDotDoNew(cs,cRepDotAreaType);
	 if(!rep) 
		ok=ErrMessage("ExportDots","Couldn't get dot representation.");
	 else {
		result=Alloc(ExportDotsObj,1);
		ErrChkPtr(result);
		result->export.fFree=(void (*)(struct Export *))ExportDotsObjFree;
		/* cannabilize the data structures */
		result->point=rep->V;
		rep->V=NULL;
		result->normal=rep->VN;
		rep->VN=NULL;
		result->type=rep->T;
		rep->T=NULL;
		result->flag=rep->F;
		rep->F=NULL;
		result->area=rep->A;
		rep->A=NULL;
		result->nPoint=rep->N;
		rep->R.fFree((Rep*)rep); /* free the remaining structures */
	 }
  }
  return(result);
}

void ExportDotsObjFree(ExportDotsObj *obj)
{
  if(obj) {
	 FreeP(obj->point);
	 FreeP(obj->normal);
	 FreeP(obj->type);
	 FreeP(obj->flag);
	 FreeP(obj->area);
  }
}

void ExportDeleteMDebug(Export *ex)
{
  if(ex) 
	 if(ex->fFree)
		ex->fFree(ex);
  FreeP(ex);
}


/*--------------------------------------------------------------------- */

