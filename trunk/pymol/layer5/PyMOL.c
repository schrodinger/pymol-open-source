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
  
  /* dynamically mapped string constants */

  OVLexicon *Lex;
  ov_word lex_pdb, lex_c_string, lex_c_filename;

  OVOneToOne *Rep;
  ov_word lex_everything, lex_sticks, lex_spheres, lex_surface;
  ov_word lex_labels, lex_nb_spheres, lex_cartoon, lex_ribbon;
  ov_word lex_lines, lex_mesh, lex_dots, lex_dashes, lex_nonbonded;
  ov_word lex_cell, lex_cgo, lex_callback, lex_extent, lex_slice;

} _CPyMOL;

static OVstatus PyMOL_InitAPI(CPyMOL *I)
{
  OVContext *C = I->G->Context;
  OVreturn_word result;
  I->Lex = OVLexicon_New(C->heap);
  if(!I->Lex) 
    return_OVstatus_FAILURE;

  I->Rep = OVOneToOne_New(C->heap);
  if(!I->Rep)
    return_OVstatus_FAILURE;

  /* the following preprocessor macros may required GNU's cpp or VC++
     we'll see... */

#define LEX(ARG)  \
  if(!OVreturn_IS_OK( (result= OVLexicon_GetFromCString(I->Lex,#ARG))))  \
    return_OVstatus_FAILURE \
    else \
      I -> lex_ ## ARG = result.word;
  
  LEX(pdb);
  LEX(c_string);
  LEX(c_filename);

  /* string constants that are accepted on input */

#define LEX_REP(NAME,CODE) LEX(NAME) \
    if(!OVreturn_IS_OK( OVOneToOne_Set(I->Rep,I->lex_ ## NAME, CODE)))  \
      return_OVstatus_FAILURE;
  
  LEX_REP(everything,-1);
  LEX_REP(sticks,0);
  LEX_REP(spheres,1);
  LEX_REP(surface,2);
  LEX_REP(labels,3);
  LEX_REP(nb_spheres,4);
  LEX_REP(cartoon,5);
  LEX_REP(ribbon,6);
  LEX_REP(lines,7);
  LEX_REP(mesh,8);
  LEX_REP(dots,9);
  LEX_REP(dashes,10);
  LEX_REP(nonbonded,11);
  LEX_REP(cell,11);
  LEX_REP(cgo,13);
  LEX_REP(callback,14);
  LEX_REP(extent,15);
  LEX_REP(slice,16);

 return_OVstatus_SUCCESS;
}

int PyMOL_NewG3DStream(CPyMOL *I,int **array_ptr)
{
  int *return_vla = ExecutiveGetG3d(I->G);
  int result = OVstatus_FAILURE;
  if(return_vla) {
    result = VLAGetSize(return_vla)*(sizeof(G3dPrimitive)/sizeof(int));
  }
  if(array_ptr)
    *array_ptr = return_vla;
  return result;
}

int PyMOL_DelG3DStream(CPyMOL *I,int *array_ptr)
{
  VLAFreeP(array_ptr);
  return OVstatus_SUCCESS;
}


static OVstatus PyMOL_PurgeAPI(CPyMOL *I)
{
  OVOneToOne_DEL_AUTO_NULL(I->Rep);
  OVLexicon_DEL_AUTO_NULL(I->Lex);
  return_OVstatus_SUCCESS;
}

int PyMOL_Delete(CPyMOL *I,char *name)
{
  ExecutiveDelete(I->G, name);
  return OVstatus_SUCCESS; /* to do: return a real result */
}

int PyMOL_Zoom(CPyMOL *I,char *selection, float buffer,
               int state, int complete, int animate)
{
  return ExecutiveWindowZoom(I->G, selection, buffer, state, complete, animate);
}

static OVreturn_word get_rep_id(CPyMOL *I,char *representation)
{
  OVreturn_word result;

  if(!OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,representation))))
    return result;
  return OVOneToOne_GetForward(I->Rep,result.word);
}

int PyMOL_Show(CPyMOL *I,char *representation, char *selection)
{
  OrthoLineType s1;
  int ok=true;
  OVreturn_word rep;
  if(OVreturn_IS_OK( (rep= get_rep_id(I,representation)))) {
    SelectorGetTmp(I->G,selection,s1);
    ExecutiveSetRepVisib(I->G,s1,rep.word,true);
    SelectorFreeTmp(I->G,s1);
  }
  return ok;
}

int PyMOL_Hide(CPyMOL *I,char *representation, char *selection)
{
  OrthoLineType s1;
  int ok=true;
  OVreturn_word rep;
  if(OVreturn_IS_OK( (rep= get_rep_id(I,representation)))) {
    SelectorGetTmp(I->G,selection,s1);
    ExecutiveSetRepVisib(I->G,s1,rep.word,false);
    SelectorFreeTmp(I->G,s1);
  }
  return ok;
}

int PyMOL_Reinitialize(CPyMOL *I)
{
  return ExecutiveReinitialize(I->G);
}

int PyMOL_Load(CPyMOL *I,char *content, char *content_type, 
               char *content_format, char *object_name, 
               int frame, int discrete, int finish, 
               int quiet, int multiplex) /* ADD zoom */
{
  OrthoLineType buf = "";
  OVreturn_word result;
  int type_code;
  int format_code;
  int zoom = -1;

  if(!OVreturn_IS_OK( (result= OVLexicon_BorrowFromCString(I->Lex,content_type))))
    return OVstatus_FAILURE;
  else
    type_code = result.word;

  if(!OVreturn_IS_OK( (result= OVLexicon_BorrowFromCString(I->Lex,content_format))))
    return OVstatus_FAILURE;
  else
    format_code = result.word;

  if((type_code != I->lex_c_filename) &&
     (type_code != I->lex_c_string)) {
    return OVstatus_FAILURE;
  } else if(format_code != I->lex_pdb) {
    return OVstatus_FAILURE;
  }
  
  /* handling of multiplex option */

  if(multiplex==-2) /* use setting default value */
    multiplex = SettingGetGlobal_i(TempPyMOLGlobals,cSetting_multiplex);
  if(multiplex<0) /* default behavior is not to multiplex */
    multiplex = 0;

  /* handing of discete option */

  if(discrete<0) {/* use default discrete behavior for the file format 
                   * this will be the case for MOL2 and SDF */ 
    if(multiplex==1) /* if also multiplexing, then default discrete
                      * behavior is not load as discrete objects */
      discrete=0;
    else
      discrete=1; /* otherwise, allow discrete to be the default */
  }

  ExecutiveProcessPDBFile(I->G,ExecutiveGetExistingCompatible(I->G,object_name,cLoadTypePDB),content,
                          object_name,frame-1,discrete,finish
                          ,buf,NULL,quiet,
                          type_code==I->lex_c_string, multiplex,zoom);
                          
  return OVstatus_SUCCESS;
}


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
  0, /* presentation */
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

  G->Context = OVContext_New();

  PyMOL_InitAPI(I);

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
  G->Ready = true; 

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

  PyMOL_PurgeAPI(I);

  OVContext_Del(G->Context);   
}
void PyMOL_Free(CPyMOL *I)
{
  /* take PyMOL down gracefully */
  PyMOLOptions_Free(I->G->Option);
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
  if(G->HaveGUI) {

    PyMOL_PushValidContext(I);

    /* get us into a well defined GL state */

    /*glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();*/
         
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_AUTO_NORMAL);
    glDisable(GL_BLEND);
    glDisable(GL_COLOR_LOGIC_OP);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_CULL_FACE);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DITHER);
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_NORMALIZE);
    glDisable(GL_POLYGON_SMOOTH);

  } 

  I->RedisplayFlag = false;
  OrthoBusyPrime(G);
  ExecutiveDrawNow(G);

  if(G->HaveGUI) PyMOL_PopValidContext(I);
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
    did_work = true;
  }

  SceneIdle(G); 

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
  if(I->SwapFn && I->G->ValidContext) {
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

void PyMOL_PushValidContext(CPyMOL *I)
{
  if(I && I->G)
    I->G->ValidContext++;
}
void PyMOL_PopValidContext(CPyMOL *I)
{
  if(I && I->G && (I->G->ValidContext>0))
    I->G->ValidContext--;
}

void PyMOL_SetDefaultMouse(CPyMOL *I)
{
  PyMOLGlobals *G = I->G;

  ButModeSet(G,0,cButModeRotXYZ);
  ButModeSet(G,1,cButModeTransXY);
  ButModeSet(G,2,cButModeTransZ);
  ButModeSet(G,12,cButModeScaleSlab);
  ButModeSet(G,13,cButModeMoveSlab);
  ButModeSet(G,5,cButModeClipNF);
  ButModeSet(G,14,cButModeMoveSlabAndZoom);
  ButModeSet(G,15,cButModeTransZ);
  ButModeSet(G,20,cButModeCent);
  ButModeSet(G,10,cButModeOrigAt);
}

