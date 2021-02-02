#pragma once

#include "Picking.h"
#include "pymol/zstring_view.h"

struct DeferredMouse;
struct PyMOLGlobals;

namespace pymol
{
struct CObject;
}

struct NamedPickContext {
  std::string name;
  int state;
};

struct NamedPicking {
  Pickable src;
  NamedPickContext context;
  NamedPicking(const Picking& pick);
};

void SceneClickObject(PyMOLGlobals* G, pymol::CObject* obj, const NamedPicking& LastPicked,
    int mode, pymol::zstring_view sel_mode_kw);
void SceneClickTransformObject(
    PyMOLGlobals* G, pymol::CObject* obj, const NamedPicking& LastPicked, int mode, bool is_single_click);
void SceneClickPickBond(PyMOLGlobals* G, int x, int y, int mode, const NamedPicking& LastPicked);
int SceneDeferredRelease(DeferredMouse* dm);
int SceneDeferredClick(DeferredMouse* dm);
void SceneClickPickNothing(PyMOLGlobals* G, int button, int mod, int mode);

