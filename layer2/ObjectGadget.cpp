
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectGadget.h"
#include"ObjectGadgetRamp.h"
#include"GadgetSet.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"VFont.h"

CGO *ObjectGadgetPyListFloatToCGO(PyObject * list);

int ObjectGadgetGetVertex(const ObjectGadget * I, int index, int base, float *v)
{
  GadgetSet *gs;
  int ok = false;
  if(I->CurGSet < I->NGSet) {
    gs = I->GSet[I->CurGSet];
    if(gs) {
      ok = GadgetSetGetVertex(gs, index, base, v);
    }
  }
  return (ok);
}

int ObjectGadgetSetVertex(ObjectGadget * I, int index, int base, const float *v)
{
  GadgetSet *gs;
  int ok = false;
  if(I->CurGSet < I->NGSet) {
    gs = I->GSet[I->CurGSet];
    if(gs) {
      ok = GadgetSetSetVertex(gs, index, base, v);
    }
  }
  if (index) // if 0 - xyz doesn't change, 1 - mouse position when changing colors
    I->Changed = true;
  return (ok);
}


/* in current state */
ObjectGadget *ObjectGadgetTest(PyMOLGlobals * G)
{
  ObjectGadget *I = NULL;
  GadgetSet *gs = NULL;
  CGO *cgo = NULL;
  int a;

  float coord[] = {
    0.5F, 0.5F, 0.0F,
    0.0F, 0.0F, 0.0F,
    0.3F, 0.0F, 0.0F,
    0.0F, -0.3F, 0.0F,
    0.3F, -0.3F, 0.0F,
    0.03F, -0.03F, 0.03F,
    0.27F, -0.03F, 0.03F,
    0.03F, -0.27F, 0.03F,
    0.27F, -0.27F, 0.03F,
    0.02F, -0.02F, 0.01F,
    0.28F, -0.02F, 0.01F,
    0.02F, -0.28F, 0.01F,
    0.28F, -0.28F, 0.01F,
  };

  float normal[] = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    -1.0, 0.0, 0.0,
    0.0, -1.0, 0.0,
  };

  I = new ObjectGadget(G);
  gs = GadgetSetNew(G);

  gs->NCoord = 13;
  gs->Coord = VLAlloc(float, gs->NCoord * 3);
  for(a = 0; a < gs->NCoord * 3; a++) {
    gs->Coord[a] = coord[a];
  }

  gs->NNormal = 5;
  gs->Normal = VLAlloc(float, gs->NNormal * 3);
  for(a = 0; a < gs->NNormal * 3; a++) {
    gs->Normal[a] = normal[a];
  }

  cgo = CGONewSized(G, 100);
  CGOColor(cgo, 1.0, 1.0, 1.0);

  /* top */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 2.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 5.0, 0.0);
  CGOVertex(cgo, 1.0, 6.0, 0.0);

  CGONormal(cgo, 2.0, 1.0, 0.0);
  CGOVertex(cgo, 1.0, 1.0, 0.0);
  CGOVertex(cgo, 1.0, 2.0, 0.0);
  CGOEnd(cgo);

  /* bottom */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 2.0, 4.0, 0.0);
  CGOVertex(cgo, 1.0, 3.0, 0.0);
  CGOVertex(cgo, 1.0, 4.0, 0.0);

  CGONormal(cgo, 2.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 7.0, 0.0);
  CGOVertex(cgo, 1.0, 8.0, 0.0);
  CGOEnd(cgo);

  /* left */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 2.0, 3.0, 0.0);
  CGOVertex(cgo, 1.0, 1.0, 0.0);
  CGOVertex(cgo, 1.0, 3.0, 0.0);

  CGONormal(cgo, 2.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 5.0, 0.0);
  CGOVertex(cgo, 1.0, 7.0, 0.0);
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 2.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 6.0, 0.0);
  CGOVertex(cgo, 1.0, 8.0, 0.0);

  CGONormal(cgo, 2.0, 0.0, 0.0);
  CGOVertex(cgo, 1.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 4.0, 0.0);
  CGOEnd(cgo);

  CGOColor(cgo, 1.0, 0.0, 0.0);

  /* center */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 2.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 5.0, 0.0);
  CGOVertex(cgo, 1.0, 7.0, 0.0);
  CGOVertex(cgo, 1.0, 6.0, 0.0);
  CGOVertex(cgo, 1.0, 8.0, 0.0);
  CGOEnd(cgo);

  CGOColor(cgo, 0.0, 1.0, 0.0);
  /* backr */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 2.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 9.0, 0.0);
  CGOVertex(cgo, 1.0, 10.0, 0.0);
  CGOVertex(cgo, 1.0, 11.0, 0.0);
  CGOVertex(cgo, 1.0, 12.0, 0.0);
  CGOEnd(cgo);
  CGOStop(cgo);

  gs->ShapeCGO = cgo;

  cgo = CGONewSized(G, 100);
  CGODotwidth(cgo, 5);

  CGOPickColor(cgo, 0, cPickableGadget);

  /* top */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, 1.0, 1.0, 0.0);
  CGOVertex(cgo, 1.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 5.0, 0.0);
  CGOVertex(cgo, 1.0, 6.0, 0.0);
  CGOEnd(cgo);

  /* bottom */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, 1.0, 3.0, 0.0);
  CGOVertex(cgo, 1.0, 4.0, 0.0);
  CGOVertex(cgo, 1.0, 7.0, 0.0);
  CGOVertex(cgo, 1.0, 8.0, 0.0);
  CGOEnd(cgo);

  /* left */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, 1.0, 1.0, 0.0);
  CGOVertex(cgo, 1.0, 3.0, 0.0);
  CGOVertex(cgo, 1.0, 5.0, 0.0);
  CGOVertex(cgo, 1.0, 7.0, 0.0);
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, 1.0, 6.0, 0.0);
  CGOVertex(cgo, 1.0, 8.0, 0.0);
  CGOVertex(cgo, 1.0, 2.0, 0.0);
  CGOVertex(cgo, 1.0, 4.0, 0.0);
  CGOEnd(cgo);
  CGOEnd(cgo);
  CGOStop(cgo);
  gs->PickShapeCGO = cgo;

  gs->Obj = I;
  gs->State = 0;

  I->GSet[0] = gs;
  I->NGSet = 1;
  gs->update();
  ObjectGadgetUpdateExtents(I);
  return (I);

}

void ObjectGadgetUpdateExtents(ObjectGadget * I)
{
  float maxv[3] = { FLT_MAX, FLT_MAX, FLT_MAX };
  float minv[3] = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
  int a;
  GadgetSet *ds;

  /* update extents */
  copy3f(maxv, I->ExtentMin);
  copy3f(minv, I->ExtentMax);
  I->ExtentFlag = false;
  for(a = 0; a < I->NGSet; a++) {
    ds = I->GSet[a];
    if(ds) {
      if(GadgetSetGetExtent(ds, I->ExtentMin, I->ExtentMax))
        I->ExtentFlag = true;
    }
  }
}

static PyObject *ObjectGadgetGSetAsPyList(ObjectGadget * I, bool incl_cgos)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NGSet);
  for(a = 0; a < I->NGSet; a++) {
    if(I->GSet[a]) {
      PyList_SetItem(result, a, GadgetSetAsPyList(I->GSet[a], incl_cgos));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(Py_None));
    }
  }
  return (PConvAutoNone(result));

}

static int ObjectGadgetGSetFromPyList(ObjectGadget * I, PyObject * list, int version)
{

  int ok = true;
  int a;
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    VLACheck(I->GSet, GadgetSet *, I->NGSet);
    for(a = 0; a < I->NGSet; a++) {
      if(ok){
        auto *val = PyList_GetItem(list, a);
        ok = GadgetSetFromPyList(I->G, val, &I->GSet[a], version);
      }
      if(ok && I->GSet[a]) {
        I->GSet[a]->Obj = I;
        I->GSet[a]->State = a;
      }
    }
  }
  return (ok);
}

int ObjectGadgetInitFromPyList(PyMOLGlobals * G, PyObject * list, ObjectGadget * I,
                               int version)
{
  int ok = true;
  if(ok)
    ok = (I != NULL) && (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok){
    auto *val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->GadgetType);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->NGSet);
  if(ok)
    ok = ObjectGadgetGSetFromPyList(I, PyList_GetItem(list, 3), version);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 4), &I->CurGSet);

  /*  ObjectGadgetInvalidateRep(I,cRepAll); */
  if(ok) {
    ObjectGadgetUpdateExtents(I);
  } else {
    /* cleanup? */
  }
  return (ok);
}

int ObjectGadgetNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectGadget ** result,
                              int version)
{
  int ok = true;
  ObjectGadget *I = NULL;
  int gadget_type = -1;
  PyObject *plain;
  (*result) = NULL;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);

  /* NOTE there is a serious screw-up here...ramp gadgets aren't saved right, but
     we've got to maintain backward compat...ugh */

  if(ok)
    ok = ((plain = PyList_GetItem(list, 0)) != NULL);
  if(ok)
    ok = PyList_Check(plain);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(plain, 1), &gadget_type);
  if(ok)
    switch (gadget_type) {      /* call the right routine to restore the gadget! */
    case cGadgetRamp:
      ok = ObjectGadgetRampNewFromPyList(G, list, (ObjectGadgetRamp **) result, version);
      break;
    case cGadgetPlain:
      I = new ObjectGadget(G);
      if(ok)
        ok = (I != NULL);
      if(ok)
        ok = ObjectGadgetInitFromPyList(G, list, I, version);
      if(ok)
        (*result) = I;
      break;
    default:
      ok = false;
      break;
    }
  return (ok);
}

PyObject *ObjectGadgetPlainAsPyList(ObjectGadget * I, bool incl_cgos)
{
  PyObject *result = NULL;

  /* first, dump the atoms */

  result = PyList_New(5);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->GadgetType));
  PyList_SetItem(result, 2, PyInt_FromLong(I->NGSet));
  PyList_SetItem(result, 3, ObjectGadgetGSetAsPyList(I, incl_cgos));
  PyList_SetItem(result, 4, PyInt_FromLong(I->CurGSet));
  return (PConvAutoNone(result));
}

PyObject *ObjectGadgetAsPyList(ObjectGadget * I)
{
  PyObject *result = NULL;

  /* first, dump the atoms */

  switch (I->GadgetType) {
  case cGadgetRamp:
    result = ObjectGadgetRampAsPyList((ObjectGadgetRamp *) I);
    break;
  case cGadgetPlain:
    result = ObjectGadgetPlainAsPyList(I);
    break;
  }
  return (PConvAutoNone(result));
}

ObjectGadget::~ObjectGadget()
{
  auto I = this;
  for(int a = 0; a < I->NGSet; a++)
    if(I->GSet[a]) {
      delete I->GSet[a];
      I->GSet[a] = NULL;
    }
}

void ObjectGadgetUpdateStates(ObjectGadget * I)
{
  int a;
  OrthoBusyPrime(I->G);
  for(a = 0; a < I->NGSet; a++)
    if(I->GSet[a]) {
      OrthoBusySlow(I->G, a, I->NGSet);
      /*           printf(" ObjectGadget: updating state %d of \"%s\".\n" , a+1, I->Name); */
      I->GSet[a]->update();
    }
}


/*========================================================================*/
void ObjectGadget::update()
{
  auto I = this;
  if(I->Changed) {
    ObjectGadgetUpdateStates(I);
    ObjectGadgetUpdateExtents(I);
    I->Changed = false;
  }
}


/*========================================================================*/

int ObjectGadget::getNFrame() const
{
  return NGSet;
}


/*========================================================================*/
void ObjectGadget::render(RenderInfo * info)
{
  auto I = this;
  int state = info->state;
  const RenderPass pass = info->pass;
  if(pass == RenderPass::Transparent || info->ray || info->pick) {

    ObjectPrepareContext(I, info);
    for(StateIterator iter(I->G, I->Setting.get(), state, I->NGSet);
        iter.next();) {
      GadgetSet * gs = I->GSet[iter.state];
      gs->render(info);
    }
  }
}


/*========================================================================*/
ObjectGadget::ObjectGadget(PyMOLGlobals * G) : pymol::CObject(G)
{
  type = cObjectGadget;
  GSet = pymol::vla<GadgetSet*>(10);        /* auto-zero */
}

pymol::RenderContext ObjectGadget::getRenderContext() const
{
  return pymol::RenderContext::UnitWindow;
}

