/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warrn Lyford Delano of DeLano Scientific. 
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

#include "os_std.h"
#include "MemoryDebug.h"

#include "Base.h"

#include "ClassPyMOL.h"
#include "PyMOLGlobals.h"
#include "PyMOLOptions.h"


PyMOLGlobals *TempPyMOLGlobals;

typedef struct _ClassPyMOL {
  PyMOLGlobals *G;
} _ClassPyMOL;

const static COption Defaults = {
  true, /* pmgui */
  true, /* internal_gui*/
  true, /* show_splash */
  1,   /* internal_feedback */
  true, /* security */
  false, /* game mode */
  0, /* force_stereo */
  640, /* winX */
  480, /* winY */
  false, /* blue_line */
  0, /* winPX */
  175, /* winPY */
  true, /* external_gui */
  true, /* siginthand */
  false, /* reuse helper */
  false, /* auto reinitialize */
  false, /* keep thread alive */
  false, /* quiet */
  false, /* incentive product */
  "", /* after_load_script */
  0, /* multisample */
  1, /* window_visible */
  0, /* read_stdin */
};

COption *ClassPyMOLOptionsNew(void)
{
  COption *result = NULL;
  result = Calloc(COption,1);
  if(result)
    *result = Defaults;
  return result;
}

void ClassPyMOLOptionsFree(COption *options)
{
  FreeP(options);
}

static ClassPyMOL *_ClassPyMOLNew(void)
{
  ClassPyMOL *result = NULL;

  /* allocate global container */

  if( (result = Calloc(ClassPyMOL,1)) ) {

    if( (result->G = Calloc(PyMOLGlobals,1)) ) {
      
      /* temporary global pointer for the transition period */

      TempPyMOLGlobals=result->G;

      /* continue initialization */

    } else {
      FreeP(result);
    }
  }
  return result;
}

ClassPyMOL *ClassPyMOLNew(void)
{
  ClassPyMOL *result = _ClassPyMOLNew();
  if(result && result->G) {
    result->G->Option = Calloc(COption,1);
    if(result->G->Option)
      (*result->G->Option) = Defaults;
  }
  return result;
}

ClassPyMOL *ClassPyMOLNewWithOptions(COption *option)
{
  ClassPyMOL *result = _ClassPyMOLNew();
  if(result && result->G) {
    result->G->Option = Calloc(COption,1);
    if(result->G->Option)
      *(result->G->Option) = *option;
  }
  return result;
}

void ClassPyMOLFree(ClassPyMOL *I)
{
  /* take PyMOL down gracefully */
  FreeP(I->G);
  FreeP(I);
}

struct _PyMOLGlobals *ClassPyMOLGetGlobals(ClassPyMOL *I)
{
  return I->G;
}
