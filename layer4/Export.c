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
#include"Scene.h"

void ExportDotsObjFree(ExportDotsObj *rec);

/*--------------------------------------------------------------------- */

#ifndef _MMechIO
void *ExportMMGet(char *obj,int state)
{
  return(NULL);
}
int ExportMMUpdate(char *obj,int state,void *mmdat)
{
  return(0);
}
void ExportMMFree(void *mmdat)
{
}
#else

/* compile routines for shuttling opaque CObjects
   to and from molecular mechanics routines */

#include"mmechio.h" 

void *ExportMMGet(char *name,int state)
{
  CMMechIO *io = NULL;
  ObjectMolecule *obj;
  CoordSet *cs;
  int a,b,a0,a1,a2,at;
  CMMechIOAtom *mmat;
  int *mmbnd,*bnd;
  AtomInfoType *ai;
  float *crd;
  obj = ExecutiveFindObjectMoleculeByName(name);
  if(obj&&(state>=0)) {
    if((state<obj->NCSet)&&(!obj->DiscreteFlag)) {
      if(obj->CSet[state]) {
        cs=obj->CSet[state];
        io= (CMMechIO*)mmalloc(sizeof(CMMechIO));
        io->nAtom=cs->NIndex;
        io->atom=Alloc(CMMechIOAtom,cs->NIndex);
        io->nBond=0;
        io->bond=VLAlloc(int,cs->NIndex*9);

        mmat = io->atom;
        crd=cs->Coord;
        for(a=0;a<cs->NIndex;a++) {
          at=cs->IdxToAtm[a]; /* alway valid */
          ai=obj->AtomInfo+at;
          mmat->type = ai->customType;
          mmat->flags = ai->flags;
          mmat->coord[0]=*(crd++);
          mmat->coord[1]=*(crd++);
          mmat->coord[2]=*(crd++);
          mmat++;
        }

        bnd=obj->Bond;

        for(b=0;b<obj->NBond;b++) {
          a0 = cs->AtmToIdx[*(bnd++)];
          a1 = cs->AtmToIdx[*(bnd++)];
          a2 = *(bnd++);
          if((a0>=0)&&(a1>=0)) {
            VLACheck(io->bond,int,3*io->nBond+2);
            mmbnd = io->bond+3*io->nBond;
            mmbnd[0] = a0;
            mmbnd[1] = a1;
            mmbnd[2] = a2; /* convey valence in case we adopt a different type of ffield */
            io->nBond++;
          }
        }
      }
    }
  }
  return((void*)io);
}

int ExportMMUpdate(char *name,int state,void *mmdat)
{
  
  CMMechIO *io = mmdat;
  int result = false;
  ObjectMolecule *obj;
  CoordSet *cs;
  int a;
  CMMechIOAtom *mmat;
  AtomInfoType *ai;
  float *crd;
  obj = ExecutiveFindObjectMoleculeByName(name);
  if(io) {
    if(!obj) {
      result=ErrMessage("MMUpdate","invalid object");
    } else {
      if((state<0)||(state>=obj->NCSet)||(obj->DiscreteFlag)) {
        result=ErrMessage("MMUpdate","invalid state for object.");
      } else  {
        if(!obj->CSet[state]) {
          result=ErrMessage("MMUpdate","empty state.");
        } else {
          cs=obj->CSet[state];
          if(cs->NIndex!=io->nAtom) {
            result=ErrMessage("MMUpdate","atom count mismatch.");
            PRINTF "MMUpdate: cset %d != mmechio %d \n",cs->NIndex,io->nAtom ENDF;
          } else {
            mmat = io->atom;
            crd=cs->Coord;
            for(a=0;a<cs->NIndex;a++) {
              *(crd++) = mmat->coord[0];
              *(crd++) = mmat->coord[1];
              *(crd++) = mmat->coord[2];
              mmat++;
            }

            if(cs->fInvalidateRep)
              cs->fInvalidateRep(cs,cRepAll,cRepInvAll);
            SceneChanged();
            result = true;
          }
        }
      }
    }
  }
  return(result);
}


void ExportMMFree(void *mmdat)
{
  CMMechIO *io=(CMMechIO*)mmdat;
  if(io) {
    FreeP(io->atom);
    VLAFreeP(io->bond);
    FreeP(io);
  }
}

#endif


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
    /*	 ExecutiveSetRepVisib(name,cRepDot,1); */
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

/*---------------------------------------------------------------------*/




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

