
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
#ifndef _H_GadgetSet
#define _H_GadgetSet

#include"Rep.h"

#include <vector>

struct ObjectGadget;

struct GadgetSet {
  ~GadgetSet();

  // methods
  void update();
  void render(RenderInfo * info);
  void invalidateRep(cRep_t type, cRepInv_t level);

  PyMOLGlobals *G;
  ObjectGadget* Obj = nullptr;        /* NOT pickled -- restore manually */
  int State = 0;                      /* NOT pickled -- restore manually */
  float *Coord = nullptr;
  float *Normal = nullptr;
  float *Color = nullptr;
  int NCoord = 0;
  int NNormal = 0;
  int NColor = 0;

  CGO *PickShapeCGO = nullptr;
  CGO *PickCGO = nullptr;
  CGO *StdCGO = nullptr;
  CGO *ShapeCGO = nullptr;
  int offsetPtOP = 0;
  int offsetPtOPick = 0;
};

GadgetSet *GadgetSetNew(PyMOLGlobals * G);
PyObject *GadgetSetAsPyList(GadgetSet * I, bool incl_cgos);
int GadgetSetFromPyList(PyMOLGlobals * G, PyObject * list, GadgetSet ** cs, int version);
int GadgetSetGetExtent(GadgetSet * I, float *mn, float *mx);
int GadgetSetGetVertex(const GadgetSet * I, int index, int base, float *v);
int GadgetSetSetVertex(GadgetSet * I, int index, int base, const float *v);
std::vector<float> GadgetSetGetCoord(const GadgetSet* I);
#endif
