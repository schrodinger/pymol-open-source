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
#include"os_predef.h"
#include"os_python.h"

#ifdef WIN32
#include<signal.h>
#endif

#include"os_std.h"
#include"os_gl.h"

#ifdef _PYMOL_MODULE
#ifdef _DRI_WORKAROUND
#include <dlfcn.h>
#endif
#endif

#include"MemoryDebug.h"
#include"MemoryCache.h"
#include"Err.h"
#include"Util.h"
#include"Selector.h"
#include"Color.h"
#include"Ortho.h"
#include"Scene.h"
#include"PyMOLObject.h"
#include"P.h"
#include"Executive.h"
#include"Word.h"
#include"RepMesh.h"
#include"ObjectMolecule.h"
#include"Control.h"
#include"Sphere.h"
#include"Setting.h"
#include"Ray.h"
#include"Util.h"
#include"Movie.h"
#include"main.h"
#include"Editor.h"
#include"SculptCache.h"
#include"Isosurf.h"
#include"Tetsurf.h"
#include"PConv.h"
#include"VFont.h"
#include"Wizard.h"
#include"Text.h"
#include"Character.h"
#include"Seq.h"
#include"Seeker.h"
#include"Texture.h"

#include"PyMOLGlobals.h"
#include"OVContext.h"
#include"ClassPyMOL.h"
#include"PyMOLOptions.h"

#ifdef _PYMOL_OSX
int *MacPyMOLReady = NULL;
COption *MacPyMOLOption = NULL;
#endif

void MainFree(void);
void MainTest(void);
void MainBusyIdle(void);
static void MainInit(PyMOLGlobals *G);
void MainReshape(int width, int height);
static void MainDrawLocked(void);
static void MainDrag(int x,int y);

static  ClassPyMOL *PyMOLInstance; /* eliminate */

static char **myArgv,*myArgvv[2],myArgvvv[1024];
static int myArgc;

struct _CMain {
  int DirtyFlag;
  int IdleMode;
  int SwapFlag;
  double IdleTime;
  int IdleCount;
  int ReshapeFlag;
  int DragDirtyFlag;
  int DragPassive;
  int Modifiers;
  int FinalInit;
  int TheWindow;
  int WindowIsVisible;
};

/* global options */

static void launch(COption *option);

void MainOnExit(void);

void MainDragDirty(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;
  I->DragDirtyFlag = 1;
}

static void DrawBlueLine(PyMOLGlobals *G)
{
  if(G->Option->blue_line) {
    GLint i;
    unsigned long buffer;
    GLint window_width, window_height;

    window_width=G->Option->winX;
    window_height=G->Option->winY;

    glPushAttrib(GL_ALL_ATTRIB_BITS);
 
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_BLEND);
    for(i = 0; i < 6; i++) glDisable(GL_CLIP_PLANE0 + i);
    glDisable(GL_COLOR_LOGIC_OP);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DITHER);
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_LINE_STIPPLE);
    glDisable(GL_SCISSOR_TEST);
#ifdef _PYMOL_OSX
    glDisable(GL_SHARED_TEXTURE_PALETTE_EXT);
#endif
    glDisable(GL_STENCIL_TEST);
#ifdef _PYMOL_OSX
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_3D);
    glDisable(GL_TEXTURE_CUBE_MAP);
    glDisable(GL_TEXTURE_RECTANGLE_EXT);
    glDisable(GL_VERTEX_PROGRAM_ARB);
#endif
  
    for(buffer = GL_BACK_LEFT; buffer <= GL_BACK_RIGHT; buffer++) {
      GLint matrixMode;
      GLint vp[4];
  
      glDrawBuffer(buffer);
  
      glGetIntegerv(GL_VIEWPORT, vp);
      glViewport(0, 0, window_width, window_height);
  
      glGetIntegerv(GL_MATRIX_MODE, &matrixMode);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
  
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
      glScalef(2.0f / window_width, -2.0f / window_height, 1.0f);
      glTranslatef(-window_width / 2.0f, -window_height / 2.0f, 0.0f);
 
      /* draw sync lines*/
      glColor3d(0.0f, 0.0f, 0.0f);
      glBegin(GL_LINES); /* Draw a background line*/
      glVertex3f(0.0F, window_height - 0.5F, 0.0F);
      glVertex3f((float)window_width, window_height - 0.5F, 0.0F);
      glEnd();
      glColor3d(0.0f, 0.0f, 1.0f);
      glBegin(GL_LINES); /* Draw a line of the correct length (the cross
			    over is about 40% across the screen from the left */
      glVertex3f(0.0f, window_height - 0.5f, 0.0f);
      if(buffer == GL_BACK_LEFT)
	glVertex3f(window_width * 0.30f, window_height - 0.5f, 0.0f);
      else
	glVertex3f(window_width * 0.80f, window_height - 0.5f, 0.0f);
      glEnd();
 
      glPopMatrix();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(matrixMode);
  
      glViewport(vp[0], vp[1], vp[2], vp[3]);
    } 
    glPopAttrib();
  }
}

#ifdef _PYMOL_OSX

void MainMovieCopyPrepare(int *width,int *height,int *length)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  PLockAPIAsGlut();
  MovieCopyPrepare(G,width,height,length);
  PUnlockAPIAsGlut();
}
int MainMovieCopyFrame(int frame,int width,int height,int rowbytes,void *ptr)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  int result;
  PLockAPIAsGlut();
  result = MovieCopyFrame(G,frame,width,height,rowbytes,ptr);
  PUnlockAPIAsGlut();
  return result;
}
void MainMovieCopyFinish(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  PLockAPIAsGlut();
  MovieCopyFinish(G);
  PUnlockAPIAsGlut();
}
void MainSceneGetSize(int *width,int *height)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  PLockAPIAsGlut();
  SceneGetWidthHeight(G,width,height);
  PUnlockAPIAsGlut();
}
int MainSceneCopy(int width,int height,int rowbytes,void *ptr)
{ 
  PyMOLGlobals *G = TempPyMOLGlobals;

  int result;
  PLockAPIAsGlut();
  result = SceneCopyExternal(G,width, height,rowbytes,(unsigned char *)ptr);
  PUnlockAPIAsGlut();
  return result;
}
/*========================================================================*/
void MainDoCommand(char *str1)
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  PLockAPIAsGlut();
  if(str1[0]!='_') { /* suppress internal call-backs */
    if(strncmp(str1,"cmd._",5)) {
      OrthoAddOutput(G,"PyMOL>");
      OrthoAddOutput(G,str1);
      OrthoNewLine(G,NULL,true);
    }
    PDo(str1);
  } else if(str1[1]==' ') { /* "_ command" suppresses echoing of command, but it is still logged */
    PDo(str1+2);    
  } else { 
    PDo(str1);
  }
  PUnlockAPIAsGlut();
}
/*========================================================================*/
void MainRunCommand(char *str1)
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  PLockAPIAsGlut();
  if(str1[0]!='_') { /* suppress internal call-backs */
    if(strncmp(str1,"cmd._",5)) {
      OrthoAddOutput(G,"PyMOL>");
      OrthoAddOutput(G,str1);
      OrthoNewLine(G,NULL,true);
      if(WordMatch(G,str1,"quit",true)==0) /* don't log quit */
        PLog(str1,cPLog_pml);
    }
    PParse(str1);
  } else if(str1[1]==' ') { /* "_ command" suppresses echoing of command, but it is still logged */
    if(WordMatch(G,str1+2,"quit",true)>=0) /* don't log quit */
      PLog(str1+2,cPLog_pml);
    PParse(str1+2);    
  } else { 
    PParse(str1);
  }
  PUnlockAPIAsGlut();
}

/*========================================================================*/
void MainFlushAsync(void)
{
  PLockAPIAsGlut();
  PFlush();
  PUnlockAPIAsGlut();
}
/*========================================================================*/
void MainFlush(void)
{
  PFlush();
}
/*========================================================================*/
void MainRunString(char *str)
{
  PBlock();
  PRunString(str);
  PUnblock();
}
/*========================================================================*/
PyObject *MainGetStringResult(char *str)
{
  PyObject *result;
  result = PyRun_String(str,Py_eval_input,P_globals,P_globals);
  return(result);
}

#endif

/*========================================================================*/

void MainOnExit(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
 /* 
     here we enter not knowing anything about the current state...
     so, no graceful exit is possible -- in fact under Window's we'll
     crash unless we take the drastic way out 
  */
  if(!G->Terminating) {
    G->Terminating=true;
	printf(" PyMOL: abrupt program termination.\n");
#ifdef WIN32
	TerminateProcess(GetCurrentProcess(),0); /* only way to avoid a crash */
#endif
    exit(EXIT_SUCCESS);
  }
}
/*========================================================================*/
int MainSavingUnderWhileIdle(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;
  return(I->IdleMode==2);
}
/*========================================================================*/
void MainResetIdle(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;
  I->IdleMode = 0;
}
/*========================================================================*/
void MainDirty(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  CMain *I = G->Main;
  /*  PRINTFD(G,FB_Main)
    " MainDirty: called.\n"
    ENDFD; */
  I->DirtyFlag=true;
  I->IdleMode = 0;
}
/*========================================================================*/
void MainSwapBuffers(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;
  I->SwapFlag=true;
}
/*========================================================================*/
void MainSetPassiveDrag(int onOrOff)
{  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;
  I->DragPassive = onOrOff;
}
/*========================================================================*/
void MainSetWindowVisibility(int mode)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  G->Option->window_visible = mode;
}
/*========================================================================*/
void MainTest(void)
{
}
/*========================================================================*/
static void MainButton(int button,int state,int x,int y)
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  static int glMod;  
  CMain *I = G->Main;

  glMod = p_glutGetModifiers();

  PLockAPIAsGlut();

  if(I->DragPassive) {
    I->DragPassive = false;
    MainDrag(x,y);
  } else {
    /* stay blocked here because Clicks->SexFrame->PParse */
    
    y=G->Option->winY-y;
    
    I->Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
      ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
      ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);
    
    switch(button) {
    case P_GLUT_BUTTON_SCROLL_FORWARD:
    case P_GLUT_BUTTON_SCROLL_BACKWARD:
      x=G->Option->winX/2;
      y=G->Option->winY/2; /* force into scene */
      break;
    }
    if(!OrthoButton(G,button,state,x,y,I->Modifiers))
      {
      }
  }

  PUnlockAPIAsGlut();

}
/*========================================================================*/
static void MainDrag(int x,int y)
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  CMain *I = G->Main;
  
  PLockAPIAsGlut();
  
  y=G->Option->winY-y;
  if(!OrthoDrag(G,x,y,I->Modifiers))
    {
    }
  
  if(I->DirtyFlag)
    if(G->HaveGUI) {
      p_glutPostRedisplay();
      I->DirtyFlag=false;
    }
  
  PUnlockAPIAsGlut();

}
/*========================================================================*/
static void MainPassive(int x,int y)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;

  if(I->DragPassive) { /* a harmless race condition -- we don't want
                          to slow Python down buy locking on passive
                          mouse motion */

    PLockAPIAsGlut();

    if((y<0)||(x<0)||(x>G->Option->winX)||(y>G->Option->winY)) {       /* release passive drag if mouse leaves window... */

      y=G->Option->winY-y;
     

      OrthoButton(G,P_GLUT_LEFT_BUTTON, P_GLUT_UP,x,y,I->Modifiers);
      I->DragPassive = false;

    } else {
      
      y=G->Option->winY-y;
      
      if(!OrthoDrag(G,x,y,I->Modifiers))
        {
        }
      
      if(I->DirtyFlag)
        if(G->HaveGUI) {
          p_glutPostRedisplay();
          I->DirtyFlag=false;
        }
    }
    
    PUnlockAPIAsGlut();
  }
  
}

/*========================================================================*/
static void MainDrawLocked(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  CMain *I = G->Main;
  
  
  if(I->DirtyFlag) {
    I->DirtyFlag=false;
  }
  
  OrthoBusyPrime(G);
  ExecutiveDrawNow(G);

  if(G->HaveGUI) {
    if(Feedback(G,FB_OpenGL,FB_Debugging)) {
      PyMOLCheckOpenGLErr("During Rendering");
    }
  }

  if(I->SwapFlag)
    {
      if(!(int)SettingGet(G,cSetting_suspend_updates))
        if(G->HaveGUI) {
          DrawBlueLine(G);
          p_glutSwapBuffers();
        }
      I->SwapFlag=false;
    }
  /* don't claim to be ready until we've drawn at least once */
  G->Ready = true; 
}
/*========================================================================*/
static void MainDraw(void)
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  PRINTFD(G,FB_Main)
    " MainDraw: called.\n"
    ENDFD;
  PLockAPIAsGlut();

  MainDrawLocked();
  PUnlockAPIAsGlut();
  PRINTFD(G,FB_Main)
    " MainDraw: completed.\n"
    ENDFD;
}
/*========================================================================*/
static void MainKey(unsigned char k, int x, int y)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;
  int glMod;

  glMod = p_glutGetModifiers();

  PRINTFD(G,FB_Main)
    " MainKey: %d %d %d\n",k,x,y
    ENDFD;
  PLockAPIAsGlut();

  I->Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  if(!WizardDoKey(G,k,x,y,I->Modifiers))
    OrthoKey(G,k,x,y,I->Modifiers);

  PUnlockAPIAsGlut();
  
}

/*========================================================================*/
static void MainSpecial(int k, int x, int y)
{
  PyMOLGlobals *G = TempPyMOLGlobals;
  CMain *I = G->Main;
  char buffer[255];
  int grabbed = false;
  static int glMod;  

  glMod = p_glutGetModifiers();
  PLockAPIAsGlut();

  I->Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  if(!grabbed)
    grabbed = WizardDoKey(G,(unsigned char)k,x,y,I->Modifiers);

  switch(k) {
    case P_GLUT_KEY_UP:
    case P_GLUT_KEY_DOWN:
      grabbed=1;
      OrthoSpecial(G,k,x,y,I->Modifiers);
      break;
    case P_GLUT_KEY_LEFT:
    case P_GLUT_KEY_RIGHT:      
      if(OrthoArrowsGrabbed(G)) {
        grabbed=1;
        OrthoSpecial(G,k,x,y,I->Modifiers);
      }
      break;
  }

  if(!grabbed) {
    sprintf(buffer,"_special %d,%d,%d,%d",k,x,y,I->Modifiers);
    PLog(buffer,cPLog_pml);
    PParse(buffer);
    PFlush();
  }
  PUnlockAPIAsGlut();
}

/* new window size or exposure */
/*========================================================================*/
void MainReshape(int width, int height) /* called by Glut */
{
  float h;
  PyMOLGlobals *G = TempPyMOLGlobals;

  G->Option->winX = width;
  G->Option->winY = height;

  h = ((float)height)/width;
  if(G->HaveGUI) glViewport(0, 0, (GLint) width, (GLint) height);
  
  OrthoReshape(G,width,height,false);

}
/*========================================================================*/

PyObject *MainAsPyList(void) 
{
  PyMOLGlobals *G = TempPyMOLGlobals;

#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result=NULL;
  int width,height;
  result = PyList_New(2);
  BlockGetSize(SceneGetBlock(G),&width,&height);
  if(SettingGetGlobal_b(G,cSetting_seq_view)&&!SettingGetGlobal_b(G,cSetting_seq_view_overlay))
    height+=SeqGetHeight(G);
  PyList_SetItem(result,0,PyInt_FromLong(width));
  PyList_SetItem(result,1,PyInt_FromLong(height));
  return(PConvAutoNone(result));
#endif
}

int MainFromPyList(PyObject *list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  int win_x,win_y;
  int ll=0;
  OrthoLineType buffer;
  if(ok) ok = (list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) ll=PyList_Size(list);
  if(ok&&(ll>=2)) {
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&win_x);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&win_y);
    if(ok) {
      sprintf(buffer,"viewport %d, %d",win_x,win_y);
      PParse(buffer);
    }
  }
  return(ok);
#endif
}
/*========================================================================*/
void MainDoReshape(int width, int height) /* called internally */
{
  int h,w;
  int internal_feedback;
  int force =false;
  PyMOLGlobals *G = TempPyMOLGlobals;

  if(width<0) {
    BlockGetSize(SceneGetBlock(G),&width,&h);
    if(SettingGetGlobal_b(G,cSetting_internal_gui))
      width+=SettingGetGlobal_i(G,cSetting_internal_gui_width);
    force = true;
  }

  if(height<0) { 
    BlockGetSize(SceneGetBlock(G),&w,&height);
    internal_feedback = (int)SettingGet(G,cSetting_internal_feedback);
    if(internal_feedback)
      height+=(internal_feedback-1)*cOrthoLineHeight+cOrthoBottomSceneMargin;
    if(SettingGetGlobal_b(G,cSetting_seq_view)&&!SettingGetGlobal_b(G,cSetting_seq_view_overlay))
      height+=SeqGetHeight(G);
    force = true;
  }
  if(G->HaveGUI) {
    p_glutReshapeWindow(width,height);
    glViewport(0, 0, (GLint) width, (GLint) height);
  }
  OrthoReshape(G,width,height,force);
  if(SettingGet(G,cSetting_full_screen))
    p_glutFullScreen();

}
/*========================================================================*/
static void MainInit(PyMOLGlobals *G) 
{
  /*  GLfloat one[4] = { 1,1,1,1 }; 
      GLfloat low[4] = { 0.20,0.20,0.20,1 };*/

  CMain *I = (G->Main = Calloc(CMain,1));
  
  I->DirtyFlag=true;
  I->IdleMode=2;
  I->IdleCount = 0;
  I->ReshapeFlag = false;
  I->DragDirtyFlag=0;
  I->DragPassive = false;

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

  I->IdleTime=(float)UtilGetSeconds(G);

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

}




/*========================================================================*/
void MainFree(void)
{
  PyMOLGlobals *G = ClassPyMOLGetGlobals(PyMOLInstance); /* temporary -- will change */

#ifdef WIN32
   int haveGUI = G->HaveGUI;
   int theWindow = G->Main->TheWindow;
#endif
#ifdef _PYMOL_OSX
   int haveGUI = G->HaveGUI;
   int theWindow = G->Main->TheWindow;
#endif

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
                     
  ClassPyMOLFree(PyMOLInstance);

  if(G->Option->show_splash) {
    MemoryDebugDump();
    printf(" PyMOL: normal program termination.\n");
  }
  

#ifdef WIN32
  if(haveGUI) p_glutDestroyWindow(theWindow);
  TerminateProcess(GetCurrentProcess(),0); /* only way to avoid a crash */
#endif
#ifdef _PYMOL_OSX
  if(haveGUI) {
    if(G->Option->game_mode) {
      p_glutLeaveGameMode();
      /* force a full-screen refresh to eliminate garbage on screen 
       * NOTE that we currently have to patch Apple's GLUT to make this work */
      p_glutInitWindowPosition(0,0);
      p_glutInitWindowSize(640,480);
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );            
      if(p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
        theWindow = p_glutCreateWindow("PyMOL Viewer");
        p_glutFullScreen();
        p_glutDestroyWindow(theWindow);
      }
    } else 
      p_glutDestroyWindow(theWindow);
  }
#endif

}
/*========================================================================*/
void MainRefreshNow(void) 
{ /* should only be called by the master thread, with a locked API */
  PyMOLGlobals *G = TempPyMOLGlobals;

  CMain *I = G->Main;
  if(I->SwapFlag)
    {
      if(G->HaveGUI) {
        DrawBlueLine(G);
        p_glutSwapBuffers();
      }
      I->SwapFlag=false;
    }
  if(I->DirtyFlag)
    {
      if(G->HaveGUI) 
        p_glutPostRedisplay();
      else
        MainDrawLocked();
      I->DirtyFlag=false;
    }
}
/*========================================================================*/

#ifdef _PYMOL_SHARP3D
static int Sharp3DLastWindowX = -1000;
#endif

void MainBusyIdle(void) 
{
  PyMOLGlobals *G = TempPyMOLGlobals;

  /* This is one of the few places in the program where we can be sure 
	* that we have the "glut" thread...glut doesn't seem to be completely
	* thread safe or rather thread consistent
   */

  CMain *I = G->Main;


  PRINTFD(G,FB_Main)
    " MainBusyIdle: called.\n"
    ENDFD;

#ifdef _PYMOL_SHARP3D
// keep the window on even coordinates to preserve L/R stereo...
 {
   int x,y;
   x = glutGet(GLUT_WINDOW_X);
   if(x!=Sharp3DLastWindowX) {
     Sharp3DLastWindowX=x;
     if(x&0x1) {
       y = glutGet(GLUT_WINDOW_Y);
       glutPositionWindow(x-1,y);
     }
   }
 }
#endif


  /* flush command and output queues */
  
  /*  PRINTFD(G,FB_Main)
    " MainBusyIdle: entered, IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;*/
  PLockAPIAsGlut();

  if(G->HaveGUI) {
    if(I->WindowIsVisible!=G->Option->window_visible) {
      I->WindowIsVisible = G->Option->window_visible;
      if(I->WindowIsVisible) {
        p_glutShowWindow();
        OrthoDirty(G);
      } else {
        p_glutHideWindow();
      }
    }
  }

  PRINTFD(G,FB_Main)
    " MainBusyIdle: got lock.\n"
    ENDFD;

  if(I->DragDirtyFlag==1) {
    I->DragDirtyFlag = 0;
    OrthoFakeDrag(G);
  }

  if(ControlIdling(G)) {
    ExecutiveSculptIterateAll(G);
	 SceneIdle(G); 
	 I->IdleMode=0;
  } else {
    if(!I->IdleMode) {
      I->IdleTime=UtilGetSeconds(G);
      I->IdleMode=1;
    }
  }

  if(SceneRovingCheckDirty(G)) {
    SceneRovingUpdate(G);
    I->IdleMode = 0;
  }
  PFlush();

  if(I->SwapFlag) {
    if(G->HaveGUI) {
      DrawBlueLine(G);
      p_glutSwapBuffers();
    }
    I->SwapFlag=false;
  }
  if(I->DirtyFlag&&
     ((!SettingGet_b(G,NULL,NULL,cSetting_defer_updates)||
       ControlIdling(G)))
     ) {
    if(G->HaveGUI) 
      p_glutPostRedisplay();
    else
      MainDrawLocked();
    I->DirtyFlag=false;
  }

  if(I->IdleMode) { /* avoid racing the CPU */
    if(I->IdleMode==1) {
      if(UtilGetSeconds(G)-I->IdleTime>SettingGet(G,cSetting_idle_delay)) { 
        I->IdleMode=2;
        if(G->HaveGUI)
          if(SettingGet(G,cSetting_cache_display))
             p_glutPostRedisplay(); /* trigger caching of the current scene */
      }
    }
    if(I->IdleMode==1)
      PSleep((int)SettingGet(G,cSetting_fast_idle)); /* fast idle - more responsive */
    else
      PSleep((int)SettingGet(G,cSetting_slow_idle)); /* slow idle - save CPU cycles */
  } else {
    PSleep((int)SettingGet(G,cSetting_no_idle)); /* give Tcl/Tk a chance to run */
  }
  PUnlockAPIAsGlut();

  #define FINAL_INIT_AT 10

  if(I->FinalInit<FINAL_INIT_AT)
	 {
      I->FinalInit=I->FinalInit+1;
      if(I->FinalInit==FINAL_INIT_AT) {

#ifndef _PYMOL_NOPY
        PBlock();
        /* restore working directory if asked to */
        PRunString("if os.environ.has_key('PYMOL_WD'): os.chdir(os.environ['PYMOL_WD'])");

        #ifdef _PYMOL_OSX
        PRunString("if os.getcwd()[-23:]=='.app/Contents/Resources': os.chdir('../../..')");
        #endif

        PRunString("launch_gui()");

        /*#ifndef _PYMOL_WX_GLUT
         PRunString("launch_gui()");
        #endif*/

        PRunString("adapt_to_hardware()");


        if(G->Option->incentive_product) { /* perform incentive product initialization (if any) */
          PyRun_SimpleString("try:\n   import ipymol\nexcept:\n   pass\n");
        }

        PRunString("exec_deferred()");
#ifdef _PYMOL_SHARP3D
        //PParse("load $TUT/1hpv.pdb;hide;show sticks;show surface;set surface_color,white;set transparency,0.5;stereo on");
        PParse("stereo on");
        //PParse("wizard demo,cartoon");
#endif
        PUnblock();
#endif
      }
    }
  if(I->ReshapeFlag) {
    MainDoReshape(G->Option->winX,G->Option->winY);
    I->ReshapeFlag=false;
  }
  if(!G->HaveGUI) {
    if(!OrthoCommandWaiting(G)) {
      if((!G->Option->keep_thread_alive)&&(!G->Option->read_stdin)) {
        I->IdleCount++;
        if(I->IdleCount==10) {
          PLockAPIAsGlut();
          PParse("_quit");
          PFlush();
          PUnlockAPIAsGlut();
        }
      }
    }
  }
  PRINTFD(G,FB_Main)
    " MainBusyIdle: leaving... IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;

}
void MainRepositionWindowDefault(PyMOLGlobals *G)
{
  p_glutPositionWindow(G->Option->winPX,G->Option->winPY);
  p_glutReshapeWindow(G->Option->winX,G->Option->winY);
}
/*========================================================================*/

#ifdef WIN32
BOOL WINAPI HandlerRoutine(
						     DWORD dwCtrlType   //  control signal type
)
{
	switch(dwCtrlType) {
	case CTRL_CLOSE_EVENT:
	case CTRL_BREAK_EVENT:
	case CTRL_C_EVENT:
		TerminateProcess(GetCurrentProcess(),0); /* only way to avoid a crash */
	
		break;
	}
	return 1;
}
#endif

#ifdef _PYMOL_SHARP3D
void sharp3d_prepare_context(void);
#endif

static void launch(COption *options)
{
  int multisample_mask = 0;
  int theWindow = 0;

  PyMOLGlobals *G = NULL;
  PyMOLInstance = ClassPyMOLNewWithOptions(options);
  G = ClassPyMOLGetGlobals(PyMOLInstance);

#ifdef _PYMOL_OSX
  MacPyMOLOption = G->Option;
  MacPyMOLReady = &G->Ready;
#endif
  
  G->HaveGUI = G->Option->pmgui;
  G->Security = G->Option->security;

  if(G->Option->multisample)
    multisample_mask = P_GLUT_MULTISAMPLE;
  
  if(G->Option->internal_gui&&(!G->Option->game_mode))
    G->Option->winX+=cOrthoRightSceneMargin;
  if(G->Option->internal_feedback&&(!G->Option->game_mode))
    G->Option->winY+= (G->Option->internal_feedback-1)*cOrthoLineHeight + cOrthoBottomSceneMargin;

  if(G->HaveGUI) {
    #ifndef _PYMOL_OSX
    atexit(MainOnExit); /* register callback to help prevent crashes
                                 when GLUT spontaneously kills us */
    #endif

#ifdef WIN32
SetConsoleCtrlHandler(
  HandlerRoutine,  // address of handler function
  true                          // handler to add or remove
);
 
#endif
#ifdef _PYMOL_SHARP3D
      sharp3d_prepare_context();
#endif

      p_glutInit(&myArgc, myArgv);

    switch(G->Option->force_stereo) {

    case -1: /* force mono */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask | P_GLUT_DOUBLE );
      G->StereoCapable = 0;
      break;

    case 0: /* default/autodetect (stereo on win/unix; mono on macs) */
#ifdef _PYMOL_SHARP3D
      G->Option->winX = 794+220;
      G->Option->winY = 547;

      glutInitDisplayMode( P_GLUT_RGBA| P_GLUT_DOUBLE| P_GLUT_STENCIL );
      G->StereoCapable = 1;
#else
#ifndef _PYMOL_OSX
      /* don't try stereo on OS X unless asked to do so */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask | P_GLUT_DOUBLE | P_GLUT_STEREO );
      if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
#endif
        p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask | P_GLUT_DOUBLE );
        if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
          p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );
        }
        G->StereoCapable = 0;
#ifndef _PYMOL_OSX
      } else {
        G->StereoCapable = 1;
      }
#endif
#endif
      break;

    case 1: /* force stereo (if possible) */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE | P_GLUT_STEREO );
      if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {

        G->StereoCapable = 0;
      } else {
        G->StereoCapable = 1;
      }      
      break;
    }

    if(!G->Option->game_mode) {
      if((G->Option->winPX>-1000)&&(G->Option->winPY>-1000)) {
        #ifndef _PYMOL_FINK
        p_glutInitWindowPosition(G->Option->winPX,G->Option->winPY);
        #else
        p_glutInitWindowPosition(G->Option->winPX,G->Option->winPY-22); /* somethings wrong with FinkGlut's window positioning...*/
        #endif
      }
      p_glutInitWindowSize(G->Option->winX, G->Option->winY);
 
      theWindow = p_glutCreateWindow("PyMOL Viewer");
      if(G->Option->window_visible) {
        p_glutShowWindow();
      } else {
        p_glutHideWindow();
      }
        
    } else {
      char str[255];
      sprintf(str,"%dx%d:32@120",G->Option->winX,G->Option->winY);
      p_glutGameModeString(str);
      p_glutEnterGameMode(); 
    }
  } 

#ifdef _PYMOL_SHARP3D
  	// Setup OpenGL
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
#endif

   MainInit(G);
   {
     CMain *I = G->Main;
     
     I->TheWindow = theWindow;

     PInit(G);
     
#ifdef _PYMOL_SHARP3D
     SettingSetGlobal_b(G,cSetting_overlay,1);
     SettingSetGlobal_b(G,cSetting_overlay_lines,1);
#endif
     
     if(G->HaveGUI) {
       p_glutDisplayFunc(         MainDraw );
       p_glutReshapeFunc(         MainReshape );
       p_glutKeyboardFunc(        MainKey );
       p_glutMouseFunc(           MainButton );
       p_glutMotionFunc(          MainDrag );
       p_glutPassiveMotionFunc(   MainPassive );
       p_glutSpecialFunc(         MainSpecial );
       p_glutIdleFunc(            MainBusyIdle );
       
       p_glutPostRedisplay();
     }
     
     PUnblock();
     
     if(G->HaveGUI) {
       SceneSetCardInfo(G,(char*)glGetString(GL_VENDOR),
                        (char*)glGetString(GL_RENDERER),
                        (char*)glGetString(GL_VERSION));
       if(G->Option->show_splash) {
         
         printf(" OpenGL-based graphics engine:\n");
         printf("  GL_VENDOR: %s\n",(char*)glGetString(GL_VENDOR));
         printf("  GL_RENDERER: %s\n",(char*)glGetString(GL_RENDERER));
         printf("  GL_VERSION: %s\n",(char*)glGetString(GL_VERSION));
         /*        printf("  GL_EXTENSIONS: %s\n",(char*)glGetString(GL_EXTENSIONS));*/
         if(G->StereoCapable) {
           printf(" Hardware stereo capability detected.\n");
         } else if((G->Option->force_stereo==1)&&(!G->StereoCapable)) {
           printf(" Hardware stereo not present (unable to force).\n");
         }
       } 
       if(!I->WindowIsVisible)
         MainReshape(G->Option->winX,G->Option->winY);
       p_glutMainLoop();
       PBlock(); /* if we've gotten here, then we're heading back to Python... */
     } else {
       SceneSetCardInfo(G,"none","ray trace only","none");
       if(G->Option->show_splash) 
         printf(" Command mode. No graphics front end.\n");
       MainReshape(G->Option->winX,G->Option->winY);
       MainDraw(); /* for command line processing */
       while(1) {
         MainBusyIdle();
         MainDraw();
       }
     }
  }
}

int MainCheckRedundantOpen(char *file)
{
  int result = false;
#ifndef _PYMOL_NOPY
  PBlock();
  result = PTruthCallStr(P_cmd,"check_redundant_open",file);
  PUnblock();
#endif
  return result;
}

/*========================================================================*/


#ifndef _PYMOL_MODULE
int main(int argc, char *argv[])
{
  myArgc=argc;
  myArgv=argv;

  PInitEmbedded(argc,argv);

#else
int was_main(void)
{
  myArgc=1;
  strcpy(myArgvvv,"pymol");
  myArgvv[0]=myArgvvv;
  myArgvv[1]=NULL;
  myArgv=myArgvv;
#ifdef _DRI_WORKAROUND
  dlopen("libGL.so.1",RTLD_LAZY|RTLD_GLOBAL);
#endif
  
#endif  

  { /* no matter how PyMOL was built, we always come through here... */

   COption *options = ClassPyMOLOptionsNew();

   if(options) {
     PGetOptions(options);

     /* below need to be phased out by modifying code to use
        PyMOLOption global */
     
#ifdef _PYMOL_SHARP3D
     //  InternalGUI = 0;
     options->internal_feedback = 0;
     options->show_splash = 0;
#endif

     launch(options);
     ClassPyMOLOptionsFree(options);

   }
 }
 return 0;
}






