
/* 
   A* -------------------------------------------------------------------
   B* This file contains source code for the PyMOL computer program
   C* Copyright (c) Schrodinger, LLC. 
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

#include "os_predef.h"


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#include <signal.h>
#include <tchar.h>
#endif

/* END PROPRIETARY CODE SEGMENT */

#include "os_std.h"
#include "os_gl.h"

#ifdef _DRI_WORKAROUND
#include <dlfcn.h>
#endif

#include "PyMOLGlobals.h"
#include "PyMOL.h"
#include "PyMOLOptions.h"

#include "MemoryDebug.h"
#include "Base.h"
#include "main.h"

#include "P.h"
#include "PConv.h"
#include "Ortho.h"
#include "Scene.h"
#include "Seq.h"
#include "Util.h"
#include "Control.h"
#include "Movie.h"
#include "Executive.h"
#include "Feedback.h"

int _gScaleFactor = 1;

PyObject *MainAsPyList(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return nullptr;
#else
  PyObject *result = nullptr;
  int width, height;
  result = PyList_New(2);
  width = SceneGetBlock(G)->getWidth();
  height = SceneGetBlock(G)->getHeight();
  if(SettingGetGlobal_b(G, cSetting_seq_view)
     && !SettingGetGlobal_b(G, cSetting_seq_view_overlay))
    height += SeqGetHeight(G);
  PyList_SetItem(result, 0, PyInt_FromLong(width));
  PyList_SetItem(result, 1, PyInt_FromLong(height));
  return (PConvAutoNone(result));
#endif
}

int MainFromPyList(PyMOLGlobals *G, PyObject * list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  int win_x, win_y;
  int ll = 0;
  OrthoLineType buffer;

  if(ok)
    ok = (list != nullptr);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  if(ok && (ll >= 2)) {
    if(!G->Option->presentation &&
        !G->Option->full_screen &&
        !ExecutiveIsFullScreen(G)) {
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 0), &win_x);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 1), &win_y);
      /* BlockGetSize(SceneGetBlock(G),&win_x,&win_y); * so how did this get into 0.98beta29?  */
      if(ok) {

        sprintf(buffer, "viewport %d, %d", win_x, win_y);
        PParse(G, buffer);
      }
    }
  }
  return (ok);
#endif
}
