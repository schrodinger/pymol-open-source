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
#include "os_gl.h"

#include "MemoryDebug.h"

#include "Base.h"

#include "OVContext.h"

#include "MemoryDebug.h"
#include "MemoryCache.h"
#include "Err.h"
#include "Util.h"
#include "Selector.h"
#include "Color.h"
#include "Ortho.h"
#include "Scene.h"
#include "PyMOLObject.h"
#include "Executive.h"
#include "Word.h"
#include "RepMesh.h"
#include "ObjectMolecule.h"
#include "Control.h"
#include "Sphere.h"
#include "Setting.h"
#include "Ray.h"
#include "Util.h"
#include "Movie.h"
#include "P.h"
#include "Editor.h"
#include "SculptCache.h"
#include "Isosurf.h"
#include "Tetsurf.h"
#include "PConv.h"
#include "VFont.h"
#include "Wizard.h"
#include "Text.h"
#include "Character.h"
#include "Seq.h"
#include "Seeker.h"
#include "Texture.h"
#include "TestPyMOL.h"

#include "PyMOL.h"
#include "PyMOLGlobals.h"
#include "PyMOLOptions.h"

PyMOLGlobals *TempPyMOLGlobals;

typedef struct _CPyMOL {
  PyMOLGlobals *G;
  int FakeDragFlag;
  int RedisplayFlag;
  int PassiveFlag;
  int SwapFlag;
  PyMOLSwapBuffersFn *SwapFn;

} _CPyMOL;

const static CPyMOLOptions Defaults = {
  true, /* pmgui */
#ifndef _PYMOL_NOPY
  true, /* internal_gui*/
#else
  false, 
#endif
#ifndef _PYMOL_NOPY
  true, /* show_splash */
#else
  false,
#endif
#ifndef _PYMOL_NOPY
  1,   /* internal_feedback */
#else
  0, 
#endif
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

CPyMOLOptions *PyMOLOptions_New(void)
{
  CPyMOLOptions *result = NULL;
  result = Calloc(CPyMOLOptions,1);
  if(result)
    *result = Defaults;
  return result;
}

void PyMOLOptions_Free(CPyMOLOptions *options)
{
  FreeP(options);
}

static CPyMOL *_PyMOL_New(void)
{
  CPyMOL *result = NULL;

  /* allocate global container */

  if( (result = Calloc(CPyMOL,1)) ) {

    if( (result->G = Calloc(PyMOLGlobals,1)) ) {
      
      result->G->PyMOL = result; /* store the instance pointer */

      #ifndef _PYMOL_NO_PY
      /* temporary global pointer for the transition period */
      TempPyMOLGlobals=result->G;
      #endif

      /* continue initialization */

    } else {
      FreeP(result);
    }
  }
  return result;
}

static void _PyMOL_Config(CPyMOL *I)
{
    I->G->HaveGUI = I->G->Option->pmgui;
    I->G->Security = I->G->Option->security;
}

CPyMOL *PyMOL_New(void)
{
  CPyMOL *result = _PyMOL_New();
  if(result && result->G) {
    result->G->Option = Calloc(CPyMOLOptions,1);
    if(result->G->Option)
      (*result->G->Option) = Defaults;
    _PyMOL_Config(result);
  }
  return result;
}

CPyMOL *PyMOL_NewWithOptions(CPyMOLOptions *option)
{
  CPyMOL *result = _PyMOL_New();
  if(result && result->G) {
    result->G->Option = Calloc(CPyMOLOptions,1);
    if(result->G->Option)
      *(result->G->Option) = *option;
    _PyMOL_Config(result);
  }
  return result;
}

void PyMOL_Start(CPyMOL *I)
{
  PyMOLGlobals *G=I->G;

  if(G->HaveGUI) {

    /* get us into a well defined GL state */
    
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
  }

  G->Context = OVContext_New();
  MemoryCacheInit(G);
  FeedbackInit(G,G->Option->quiet);
  WordInit(G);
  UtilInit(G);
  ColorInit(G);
  CGORendererInit(G);
  SettingInitGlobal(G,true,true);  
  SettingSetGlobal_i(G,cSetting_internal_gui,G->Option->internal_gui);
  SettingSetGlobal_i(G,cSetting_internal_feedback,G->Option->internal_feedback);
  TextureInit(G);
  TextInit(G);
  CharacterInit(G);
  SphereInit(G);
  OrthoInit(G,G->Option->show_splash);
  WizardInit(G); /* must come after ortho */
  MovieInit(G);
  SceneInit(G);
  SelectorInit(G);
  SeqInit(G);
  SeekerInit(G);
  ButModeInit(G);
  ControlInit(G);
  AtomInfoInit(G);
  SculptCacheInit(G);
  VFontInit(G);
  ExecutiveInit(G);
  IsosurfInit(G);
  TetsurfInit(G);
  EditorInit(G);

  I->RedisplayFlag = true;

}

void PyMOL_Stop(CPyMOL *I)
{
  PyMOLGlobals *G=I->G;
  G->Terminating=true;
  TetsurfFree(G);
  IsosurfFree(G);
  WizardFree(G);
  SceneCleanupStereo(G);
  EditorFree(G);
  ExecutiveFree(G);
  VFontFree(G);
  SculptCacheFree(G);
  AtomInfoFree(G);
  ButModeFree(G);
  ControlFree(G);
  SeekerFree(G);
  SeqFree(G);
  SelectorFree(G);
  SceneFree(G);
  MovieFree(G);
  OrthoFree(G);
  SettingFreeGlobal(G);
  CharacterFree(G);
  TextFree(G);
  TextureFree(G);
  SphereFree(G);
  PFree();
  CGORendererFree(G);
  ColorFree(G);
  UtilFree(G);
  WordFree(G);
  FeedbackFree(G);
  MemoryCacheDone(G);
  OVContext_Del(G->Context);
   
}
void PyMOL_Free(CPyMOL *I)
{
  /* take PyMOL down gracefully */
  FreeP(I->G);
  FreeP(I);
}

struct _PyMOLGlobals *PyMOL_GetGlobals(CPyMOL *I)
{
  return I->G;
}

void PyMOL_Draw(CPyMOL *I)
{
  PyMOLGlobals *G = I->G;

  I->RedisplayFlag = false;
  OrthoBusyPrime(G);
  ExecutiveDrawNow(G);
  /* don't claim to be ready until we've drawn at least once */
  G->Ready = true; 

}

void PyMOL_Key(CPyMOL *I,unsigned char k, int x, int y, int modifiers)
{
  PyMOLGlobals *G = I->G;

  if(!WizardDoKey(G,k,x,y,modifiers))
    OrthoKey(G,k,x,y,modifiers);
}


void PyMOL_Special(CPyMOL *I,int k, int x, int y, int modifiers)
{
  PyMOLGlobals *G = I->G;

  int grabbed = false;
  char buffer[255];
  
  if(!grabbed)
    grabbed = WizardDoKey(G,(unsigned char)k,x,y,modifiers);
  
  switch(k) {
  case P_GLUT_KEY_UP:
  case P_GLUT_KEY_DOWN:
    grabbed=1;
    OrthoSpecial(G,k,x,y,modifiers);
    break;
  case P_GLUT_KEY_LEFT:
  case P_GLUT_KEY_RIGHT:      
    if(OrthoArrowsGrabbed(G)) {
      grabbed=1;
      OrthoSpecial(G,k,x,y,modifiers);
    }
    break;
  }
  
  if(!grabbed) {
    sprintf(buffer,"_special %d,%d,%d,%d",k,x,y,modifiers);
    PLog(buffer,cPLog_pml);
    PParse(buffer);
    PFlush();
  }
}

void PyMOL_Reshape(CPyMOL *I,int width, int height, int force)
{
  
  PyMOLGlobals *G = I->G;

  G->Option->winX = width;
  G->Option->winY = height;

  OrthoReshape(G,width,height,force);
}

int PyMOL_Idle(CPyMOL *I)
{
  PyMOLGlobals *G = I->G;
  int did_work = false;

  if(I->FakeDragFlag==1) {
    I->FakeDragFlag = false;
    OrthoFakeDrag(G);
    did_work = true;
  }

  if(ControlIdling(G)) {
    ExecutiveSculptIterateAll(G);
	 SceneIdle(G); 
    did_work = true;
  }

  if(SceneRovingCheckDirty(G)) {
    SceneRovingUpdate(G);
    did_work = true;
  }

  PFlush();
  return did_work;
}

void PyMOL_NeedFakeDrag(CPyMOL *I)
{
  I->FakeDragFlag = true;
}

void PyMOL_NeedRedisplay(CPyMOL *I)
{
  I->RedisplayFlag = true;
}

void PyMOL_NeedSwap(CPyMOL *I)
{
  I->SwapFlag = true;
}

void PyMOL_SetPassive(CPyMOL *I,int onOff)
{
  I->PassiveFlag = onOff;
}

int PyMOL_GetRedisplay(CPyMOL *I, int reset)
{
  PyMOLGlobals *G = I->G;
  int result = I->RedisplayFlag;

  if(result) {
    if(SettingGet_b(G,NULL,NULL,cSetting_defer_updates)) {
      result = false;
    } else {
      if(reset)
        I->RedisplayFlag = false;
    }
  }
  return result;
}

int PyMOL_GetPassive(CPyMOL *I, int reset)
{
  int result = I->PassiveFlag;
  if(reset)
    I->PassiveFlag = false;
  return result;
}

int PyMOL_GetSwap(CPyMOL *I, int reset)
{
  int result = I->SwapFlag;
  if(reset)
    I->SwapFlag = false;
  return result;
}

void PyMOL_Drag(CPyMOL *I,int x, int y, int modifiers)
{
  OrthoDrag(I->G,x,y,modifiers);
}

void PyMOL_Button(CPyMOL *I,int button, int state,int x, int y, int modifiers)
{
  OrthoButton(I->G,button,state,x,y,modifiers);
}

void PyMOL_SetSwapBuffersFn(CPyMOL *I, PyMOLSwapBuffersFn *fn)
{
  I->SwapFn = fn;
}

void PyMOL_SwapBuffers(CPyMOL *I)
{
  if(I->SwapFn) {
    I->SwapFn();
    I->SwapFlag = false;
  } else {
    I->SwapFlag = true;
  }
}

void PyMOL_RunTest(CPyMOL *I, int group, int test)
{
  TestPyMOLRun(I->G, group, test);
}
