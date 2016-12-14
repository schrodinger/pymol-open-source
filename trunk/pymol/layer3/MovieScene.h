/*
 * PyMOL Scene C++ implementation
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifndef _H_MOVIESCENE
#define _H_MOVIESCENE

#include "os_python.h"
#include "Base.h"
#include "PyMOLGlobals.h"

#include <vector>
#include <string>

bool MovieSceneFunc(PyMOLGlobals * G, const char * name,
    const char * action,
    const char * message = "",
    bool store_view = true,
    bool store_color = true,
    bool store_active = true,
    bool store_rep = true,
    bool store_frame = true,
    float animate = -1.0,
    const char * new_key = "",
    bool hand = true,
    const char * sele = "all");

bool MovieSceneRecall(PyMOLGlobals * G, const char * name, float animate = -1.0,
    bool recall_view = true,
    bool recall_color = true,
    bool recall_active = true,
    bool recall_rep = true,
    bool recall_frame = true,
    const char * sele = "all");

bool MovieSceneOrder(PyMOLGlobals * G, const char * names,
    bool sort = false,
    const char * location = NULL /* "current" */);

const std::vector<std::string> & MovieSceneGetOrder(PyMOLGlobals * G);

void MovieScenesInit(PyMOLGlobals * G);
void MovieScenesFree(PyMOLGlobals * G);

PyObject * MovieScenesAsPyList(PyMOLGlobals * G);
void MovieScenesFromPyList(PyMOLGlobals * G, PyObject * o);

#endif
