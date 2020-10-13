#pragma once

#include "Picking.h"
#include "pymol/zstring_view.h"

struct DeferredMouse;
struct PyMOLGlobals;

void SceneClickObject(PyMOLGlobals* G, CObject* obj, Picking LastPicked,
    int mode, pymol::zstring_view sel_mode_kw);
int SceneDeferredRelease(DeferredMouse* dm);
int SceneDeferredClick(DeferredMouse* dm);
