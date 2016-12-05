
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

#ifndef _H_IOManager
#define _H_IOManager

#include "PyMOLGlobals.h"
#include "ObjectMolecule.h"
#include "ObjectMap.h"

const char * PlugIOManagerFindPluginByExt(PyMOLGlobals * G, const char * ext, int mask=0);

#ifdef __cplusplus
extern "C" {
#endif

int PlugIOManagerInit(PyMOLGlobals * G);
int PlugIOManagerFree(PyMOLGlobals * G);
int PlugIOManagerLoadTraj(PyMOLGlobals * G, ObjectMolecule * obj,
                          const char *fname, int frame,
                          int interval, int average, int start,
                          int stop, int max, const char *sele, int image,
                          float *shift, int quiet, const char *plugin_type);
ObjectMap *PlugIOManagerLoadVol(PyMOLGlobals * G, ObjectMap * obj,
    const char *fname, int state, int quiet, const char *plugin_type);
ObjectMolecule *PlugIOManagerLoadMol(PyMOLGlobals * G, ObjectMolecule *origObj,
    const char *fname, int state, int quiet, const char *plugin_type);
CObject * PlugIOManagerLoad(PyMOLGlobals * G, CObject ** obj_ptr,
    const char *fname, int state, int quiet, const char *plugin_type);

#ifdef __cplusplus
}
#endif

#endif
