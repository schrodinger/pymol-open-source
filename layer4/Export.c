
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
#include"os_python.h"

#include"os_predef.h"

#include"os_std.h"

#include"MemoryDebug.h"
#include"Base.h"
#include"Err.h"
#include"Export.h"
#include"RepDot.h"
#include"Executive.h"
#include"ObjectMolecule.h"
#include"Scene.h"

void ExportDotsObjFree(PyMOLGlobals * G, ExportDotsObj * rec);


/*--------------------------------------------------------------------- */


/* routines for shuttling coordinates to and from molecular mechanics routines */

ExportCoords *ExportCoordsExport(PyMOLGlobals * G, char *name, int state, int order)
{
  ExportCoords *io = NULL;
  ObjectMolecule *obj;
  CoordSet *cs;
  int a, a0;
  float *crd0, *crd1;

  obj = ExecutiveFindObjectMoleculeByName(G, name);
  if(obj && (state >= 0) && (state < obj->NCSet) && (!obj->DiscreteFlag)
     && obj->CSet[state]) {
    cs = obj->CSet[state];
    io = (ExportCoords *) mmalloc(sizeof(ExportCoords));

    if(io) {
      io->nAtom = cs->NIndex;
      io->coord = Alloc(float, cs->NIndex * 3);

      if(io->coord) {
        crd0 = cs->Coord;
        crd1 = io->coord;
        if(order) {
          /* Coordinate Set Order */
          for(a = 0; a < cs->NIndex; a++) {
            *(crd1++) = *(crd0++);
            *(crd1++) = *(crd0++);
            *(crd1++) = *(crd0++);
          }
        } else {
          /* PyMOL Atom Order */
          for(a = 0; a < obj->NAtom; a++) {
            a0 = cs->AtmToIdx[a];
            if(a0 >= 0) {
              crd0 = cs->Coord + 3 * a0;
              *(crd1++) = *(crd0++);
              *(crd1++) = *(crd0++);
              *(crd1++) = *(crd0++);
            }
          }
        }
      }
    }
  }
  return ((void *) io);
}

int ExportCoordsImport(PyMOLGlobals * G, char *name, int state, ExportCoords * io,
                       int order)
{
  int result = false;
  ObjectMolecule *obj;
  CoordSet *cs;
  int a, a0, cc;
  float *crd0, *crd1;

  obj = ExecutiveFindObjectMoleculeByName(G, name);
  if(io) {
    if(!obj) {
      result = ErrMessage(G, "ExportCoordsImport", "invalid object");
    } else {
      if((state < 0) || (state >= obj->NCSet) || (obj->DiscreteFlag)) {
        result = ErrMessage(G, "ExportCoordsImport", "invalid state for object.");
      } else {
        if(!obj->CSet[state]) {
          result = ErrMessage(G, "ExportCoordsImport", "empty state.");
        } else {
          cs = obj->CSet[state];
          if(cs->NIndex != io->nAtom) {
            result = ErrMessage(G, "ExportCoordsImport", "atom count mismatch.");
            PRINTF "ExportCoordsImport: cset %d != io %d \n", cs->NIndex,
              io->nAtom ENDF(G);
          } else {

            crd0 = cs->Coord;
            crd1 = io->coord;
            if(order) {
              /* Coordinate Set Ordering */
              for(a = 0; a < cs->NIndex; a++) {
                *(crd0++) = *(crd1++);
                *(crd0++) = *(crd1++);
                *(crd0++) = *(crd1++);
              }
            } else {
              cc = cs->NIndex;  /* array range safety */
              /* PyMOL Atom Ordering */
              for(a = 0; a < obj->NAtom; a++) {
                a0 = cs->AtmToIdx[a];
                if((a0 >= 0) && (cc--)) {
                  crd0 = cs->Coord + 3 * a0;
                  *(crd0++) = *(crd1++);
                  *(crd0++) = *(crd1++);
                  *(crd0++) = *(crd1++);
                }
              }
            }

            if(cs->fInvalidateRep)
              cs->fInvalidateRep(cs, cRepAll, cRepInvAll);
            SceneChanged(G);
            result = true;
          }
        }
      }
    }
  }
  return (result);
}

void ExportCoordsFree(ExportCoords * io)
{
  if(io) {
    FreeP(io->coord);
    FreeP(io);
  }
}

ExportDotsObj *ExportDots(PyMOLGlobals * G, char *name, int csIndex)
{
  CObject *obj;
  ObjectMolecule *objMol;
  RepDot *rep;
  CoordSet *cs = NULL;
  ExportDotsObj *result = NULL;
  int ok = true;

  obj = ExecutiveFindObjectByName(G, name);
  if(!obj)
    ok = ErrMessage(G, "ExportDots", "Not a valid object name.");
  else if(obj->type != cObjectMolecule)
    ok = ErrMessage(G, "ExportDots", "Not molecule object.");

  if(ok) {
    /*   ExecutiveSetRepVisib(G,name,cRepDot,1); */
    objMol = (ObjectMolecule *) obj;
    cs = ObjectMoleculeGetCoordSet(objMol, csIndex);
    if(!cs)
      ok = ErrMessage(G, "ExportDots", "Invalid coordinate set number.");
  }

  if(ok) {
    rep = (RepDot *) RepDotDoNew(cs, cRepDotAreaType, -1);
    if(!rep)
      ok = ErrMessage(G, "ExportDots", "Couldn't get dot representation.");
    else {
      result = Alloc(ExportDotsObj, 1);
      ErrChkPtr(G, result);
      result->export.fFree = (void (*)(struct Export *)) ExportDotsObjFree;
      /* cannabilize the data structures */
      result->point = rep->V;
      rep->V = NULL;
      result->normal = rep->VN;
      rep->VN = NULL;
      result->type = rep->T;
      rep->T = NULL;
      result->flag = rep->F;
      rep->F = NULL;
      result->area = rep->A;
      rep->A = NULL;
      result->nPoint = rep->N;
      rep->R.fFree((Rep *) rep);        /* free the remaining structures */
    }
  }
  return (result);
}


/*---------------------------------------------------------------------*/

void ExportDotsObjFree(PyMOLGlobals * G, ExportDotsObj * obj)
{
  if(obj) {
    FreeP(obj->point);
    FreeP(obj->normal);
    FreeP(obj->type);
    FreeP(obj->flag);
    FreeP(obj->area);
  }
}

void ExportDeleteMDebug(PyMOLGlobals * G, Export * ex)
{
  if(ex)
    if(ex->fFree)
      ex->fFree(ex);
  FreeP(ex);
}


/*--------------------------------------------------------------------- */
