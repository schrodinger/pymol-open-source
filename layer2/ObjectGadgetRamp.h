
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
#ifndef _H_ObjectGadgetRamp
#define _H_ObjectGadgetRamp

#include"os_python.h"

#include "ObjectGadget.h"

struct ObjectMap;
struct ObjectMolecule;

#define cRampNone 0
#define cRampMap 1
#define cRampMol 2

struct ObjectGadgetRamp : public ObjectGadget {
  int RampType = 0;

  int NLevel = 0;
  pymol::vla<float> Level;
  pymol::vla<float> LevelTmp;
  pymol::vla<float> Color;
  int var_index = 0;

  /* cRampMap */

  ObjectNameType SrcName;
  int SrcState;

  int CalcMode = 0;

  /* fields below are not saved in session */
  ObjectMap *Map = nullptr;
  ObjectMolecule *Mol = nullptr;

  float border = 0.018f;
  float width = 0.9f;
  float height = 0.06f;
  float bar_height = 0.03f;
  float text_raise = 0.003f;
  float text_border = 0.004f;
  float text_scale_h = 0.04f;
  float text_scale_v = 0.02f;
  float x = (1.0f - (width + 2 * border)) / 2.0f;
  float y = 0.12f;
  ObjectGadgetRamp(PyMOLGlobals* G);
  ~ObjectGadgetRamp();

  // virtual methods
  void update() override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
};

#define cRAMP_TRADITIONAL 1
#define cRAMP_SLUDGE 2
#define cRAMP_OCEAN 3
#define cRAMP_HOT 4
#define cRAMP_GRAYABLE 5
#define cRAMP_RAINBOW 6
#define cRAMP_AFMHOT 7
#define cRAMP_GRAYSCALE 8

ObjectGadgetRamp *ObjectGadgetRampMapNewAsDefined(PyMOLGlobals * G,
                                                  ObjectGadgetRamp *I,
                                                  ObjectMap * map,
                                                  pymol::vla<float>&& level_vla,
                                                  pymol::vla<float>&& color_vla,
                                                  int map_state, float *vert_vla,
                                                  float beyond, float within, float sigma,
                                                  int zero, int calc_mode);

ObjectGadgetRamp *ObjectGadgetRampMolNewAsDefined(PyMOLGlobals * G,
                                                  ObjectGadgetRamp *I,
                                                  ObjectMolecule * mol,
                                                  pymol::vla<float>&& level_vla,
                                                  pymol::vla<float>&& color_vla,
                                                  int mol_state, int calc_mode);

int ObjectGadgetRampInterpolate(ObjectGadgetRamp * I, float level, float *color);
int ObjectGadgetRampInterVertex(ObjectGadgetRamp * I, const float *pos, float *color,
                                int state);

PyObject *ObjectGadgetRampAsPyList(ObjectGadgetRamp * I);
int ObjectGadgetRampNewFromPyList(PyMOLGlobals * G, PyObject * list,
				      ObjectGadgetRamp ** result, int version);

void ObjectGadgetRampFree(ObjectGadgetRamp * I);
#endif
