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
-* sc
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"os_predef.h"

#include"os_std.h"
#include"os_gl.h"
#include"os_python.h"

#include"Util.h"

#include"Word.h"
#include"main.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Matrix.h"
#include"ListMacros.h"
#include"PyMOLObject.h"
#include"Scene.h"
#include"Ortho.h"
#include"Vector.h"
#include"ButMode.h"
#include"Control.h"
#include"Selector.h"
#include"Setting.h"
#include"Movie.h"
#include"MyPNG.h"
#include"P.h"
#include"Editor.h"
#include"Executive.h"
#include"Wizard.h"
#include"CGO.h"
#include"ObjectGadget.h"
#include"Seq.h"
#include"Menu.h"
#include"View.h"
#include"ObjectSlice.h"
#include"Text.h"
#include"PyMOLOptions.h"
#include"PyMOL.h"
#include"PConv.h"
#include"ScrollBar.h"

#ifdef _PYMOL_IP_EXTRAS
#include "IncentiveCopyToClipboard.h"
#endif

#ifdef _PYMOL_SHARP3D
#define cSliceMin 0.1F
#else
#define cSliceMin 1.0F
#endif

#define SceneLineHeight 127
#define SceneTopMargin 0
#define SceneBottomMargin 3
#define SceneLeftMargin 3

typedef struct ObjRec {
  CObject *obj;  
  struct ObjRec *next;
  int slot;
} ObjRec;

typedef struct {
  CDeferred deferred;
  Block *block;
  int button;
  int x;
  int y;
  int mod;
  double when;
  int mode_override;
} DeferredMouse;

typedef struct {
  CDeferred deferred;
  PyMOLGlobals *G;
  int width;
  int height;
  char *filename; /* NOTE: on heap! must free when done */
  int quiet;
  int antialias;
  float dpi;
  int entire_window;
  int format;
} DeferredImage;

typedef struct {
  CDeferred deferred;
  PyMOLGlobals *G;
  int ray_width;
  int ray_height;
  int mode;
  float angle;
  float shift;
  int quiet;
  int show_timing;
  int antialias;
} DeferredRay;

typedef struct {
  int len;
  char *name;
  int x1,y1,x2,y2,drawn;
} SceneElem;

/* allow up to 10 seconds at 30 FPS */

#define TRN_BKG 0x30
#define MAX_ANI_ELEM 300

struct _CScene {
  Block *Block;
  ObjRec *Obj;
  float RotMatrix[16]; /* WARNING: column major, as per OpenGL spec */
  float InvMatrix[16]; /* WARNING: column major, as per OpenGL spec */
  float ModMatrix[16]; /* WARNING: column major, as per OpenGL spec */
  float ProMatrix[16]; /* WARNING: column major, as per OpenGL spec */
  float PmvMatrix[16];
  float Scale;
  int Width,Height;
  int Button;
  int LastX,LastY;
  int StartX,StartY;
  int LastWinX,LastWinY;
  double LastClickTime;
  int LastButton, LastMod;
  int PossibleSingleClick;
  double LastReleaseTime;
  double SingleClickDelay;
  float ViewNormal[3],LinesNormal[3];
  float Pos[3],Origin[3];
  float H;
  float Front,Back,FrontSafe,BackSafe;
  float TextColor[3];
  double SweepTime;
  int DirtyFlag;
  int ChangedFlag;
  int CopyType,CopyNextFlag,CopyForced;
  int NFrame,HasMovie;
  ImageType *Image;
  int MovieOwnsImageFlag;
  int MovieFrameFlag;
  double LastRender,RenderTime,LastFrameTime,LastFrameAdjust;
  double LastSweep,LastSweepTime;
  float LastSweepX,LastSweepY;
  int RockFrame;
  Picking LastPicked;
  int StereoMode;
  OrthoLineType vendor,renderer,version;
  int SculptingFlag,SculptingSave;
  int RovingDirtyFlag;
  int RovingCleanupFlag;
  double RovingLastUpdate;
  int Threshold, ThresholdX, ThresholdY;
  CView *View;
  float LastPickVertex[3];
  int LastPickVertexFlag;
  int LoopFlag;
  int LoopMod;
  BlockRect LoopRect;
  CViewElem ani_elem[MAX_ANI_ELEM+1];
  int cur_ani_elem, n_ani_elem;
  int LastStateBuilt;
  int AnimationStartFlag;
  double AnimationStartTime;
  double AnimationLagTime;
  int AnimationStartFrame;
  double ApproxRenderTime;
  float VertexScale;
  float FogStart;
  float FogEnd;
  
  /* Scene Names */
  int ButtonsShown, ButtonDrag, ButtonMargin, ButtonsValid;
  int Over, Pressed, PressMode, HowFarDown, NSkip;
  int ScrollBarActive;
  int ReorderFlag;
  OrthoLineType ReorderLog;
  struct CScrollBar *ScrollBar;
  char *SceneNameVLA;
  SceneElem *SceneVLA;
  int NScene;
  CGO *AlphaCGO;

  int *SlotVLA;

  int StencilValid;
};

typedef struct {
  float unit_left,unit_right,unit_top,unit_bottom,unit_front,unit_back;
} SceneUnitContext;

static void ScenePrepareUnitContext(SceneUnitContext *context,int width,int height);

static void SceneDoRoving(PyMOLGlobals *G,float old_front,
                          float old_back,float old_origin,
                          int adjust_flag,int zoom_flag);

static float SceneGetExactScreenVertexScale(PyMOLGlobals *G,float *v1);
static void SceneRestartPerfTimer(PyMOLGlobals *G);
static void SceneRotateWithDirty(PyMOLGlobals *G,float angle,float x,float y,float z,int dirty);
static void SceneClipSetWithDirty(PyMOLGlobals *G,float front,float back,int dirty);

typedef struct {
  int n_col;
  int n_row;
  int first_slot;
  int last_slot;
  float asp_adjust;
  int active;
  int size;
  int slot;
  int mode;
  GLint cur_view[4];
  SceneUnitContext context; /* for whole-window display */
} GridInfo;

static void GridSetRayViewport(GridInfo *I, int slot, int *x, int *y, int *width, int *height)
{
  if(slot)
    I->slot = slot + I->first_slot - 1;
  else
    I->slot = slot;
  /* if we are in grid mode, then prepare the grid slot viewport */
  if(slot<0) {
    *x = I->cur_view[0];
    *y = I->cur_view[1];
    *width = I->cur_view[2];
    *height = I->cur_view[3];
  } else if(!slot) {
    int vx = 0;
    int vw = I->cur_view[2]/I->n_col;
    int vy = 0;
    int vh = I->cur_view[3]/I->n_row;
    if(I->n_col < I->n_row) {
      vw *= I->n_col;
      vh *= I->n_col;
    } else {
      vw *= I->n_row;
      vh *= I->n_row;
    }
    vx += I->cur_view[0] + (I->cur_view[2] - vw)/2;
    vy += I->cur_view[1];
    *x = vx;
    *y = vy;
    *width = vw;
    *height = vh;
  } else {
    int abs_grid_slot = slot - I->first_slot;
    int grid_col = abs_grid_slot % I->n_col;
    int grid_row = (abs_grid_slot / I->n_col);
    int vx = (grid_col*I->cur_view[2])/I->n_col;
    int vw = ((grid_col+1)*I->cur_view[2])/I->n_col - vx;
    int vy = I->cur_view[3] - ((grid_row+1)*I->cur_view[3])/I->n_row;
    int vh = (I->cur_view[3] - ((grid_row)*I->cur_view[3])/I->n_row) - vy;
    vx += I->cur_view[0];
    vy += I->cur_view[1];
    *x = vx;
    *y = vy;
    *width = vw;
    *height = vh;
  }
}
static void GridSetGLViewport(GridInfo *I, int slot)
{
  if(slot)
    I->slot = slot + I->first_slot - 1;
  else
    I->slot = slot;
  /* if we are in grid mode, then prepare the grid slot viewport */
  if(slot<0) {
    glViewport(I->cur_view[0],I->cur_view[1],I->cur_view[2],I->cur_view[3]);
  } else if(!slot) {
    int vx = 0;
    int vw = I->cur_view[2]/I->n_col;
    int vy = 0;
    int vh = I->cur_view[3]/I->n_row;
    if(I->n_col < I->n_row) {
      vw *= I->n_col;
      vh *= I->n_col;
    } else {
      vw *= I->n_row;
      vh *= I->n_row;
    }
    vx += I->cur_view[0] + (I->cur_view[2] - vw)/2;
    vy += I->cur_view[1];
    glViewport(vx,vy,vw,vh);
    ScenePrepareUnitContext(&I->context,vw,vh);
  } else {
    int abs_grid_slot = slot - I->first_slot;
    int grid_col = abs_grid_slot % I->n_col;
    int grid_row = (abs_grid_slot / I->n_col);
    int vx = (grid_col*I->cur_view[2])/I->n_col;
    int vw = ((grid_col+1)*I->cur_view[2])/I->n_col - vx;
    int vy = I->cur_view[3] - ((grid_row+1)*I->cur_view[3])/I->n_row;
    int vh = (I->cur_view[3] - ((grid_row)*I->cur_view[3])/I->n_row) - vy;
    vx += I->cur_view[0];
    vy += I->cur_view[1];
    glViewport(vx,vy,vw,vh);
    ScenePrepareUnitContext(&I->context,vw,vh);
  }
}

static void GridGetRayViewport(GridInfo *I, int width, int height)
{
  I->cur_view[0] = 0;
  I->cur_view[1] = 0;
  I->cur_view[2] = width;
  I->cur_view[3] = height;
}

static void GridGetGLViewport(GridInfo *I)
{
  glGetIntegerv(GL_VIEWPORT,I->cur_view);
}

static void GridUpdate(GridInfo *I, float asp_ratio, int mode, int size)
{
  if(mode) {
    I->size = size;
    I->mode = mode;
    {
      register int n_row = 1;
      register int n_col = 1;
      register int r_size = size;
      while((n_row * n_col) < r_size) {
        register float asp1 = asp_ratio * (n_row+1.0)/n_col;
        register float asp2 = asp_ratio * (n_row)/(n_col+1.0);
        if(asp1<1.0F) asp1 = 1.0/asp1;
        if(asp2<1.0F) asp2 = 1.0/asp2;
        if(fabs(asp1) > fabs(asp2))
          n_col++;
        else
          n_row++;
      }
      I->n_row = n_row;
      I->n_col = n_col;
    }
    if(I->size>1) {
      I->active = true;
      I->asp_adjust = (float)I->n_row / I->n_col;
      I->first_slot = 1;
      I->last_slot = I->size;
    }
  }
}

static void SceneInvalidateStencil(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  I->StencilValid = false;
}


static int SceneGetGridSize(PyMOLGlobals *G,int grid_mode)
{
  CScene *I=G->Scene;
  int slot;
  int size = 0;
    
  switch(grid_mode) {
  case 1:
    if(!I->SlotVLA)
      I->SlotVLA = VLACalloc(int,1);
    else {
      UtilZeroMem(I->SlotVLA,sizeof(int)*VLAGetSize(I->SlotVLA));
    }
    {
      int max_slot = 0;
      ObjRec *rec = NULL;
      while(ListIterate(I->Obj,rec,next)) {
        if( (slot = rec->obj->grid_slot) ) {
          slot = rec->obj->grid_slot;
          if(max_slot<slot)
            max_slot = slot;
          if(slot>0) {
            VLACheck(I->SlotVLA,int,slot);
            I->SlotVLA[slot] = 1;
          }
        }
      }
      for(slot=0;slot<=max_slot;slot++) {
        if(I->SlotVLA[slot])
          I->SlotVLA[slot] = ++size;
      }
    }
    break;
  case 2:
    if(I->SlotVLA) {
      VLAFreeP(I->SlotVLA); 
      I->SlotVLA = NULL;
    }
    {
      int max_slot = 0;
      ObjRec *rec = NULL;
      while(ListIterate(I->Obj,rec,next)) {
        if(rec->obj->fGetNFrame) {
          slot = rec->obj->fGetNFrame(rec->obj);
          if(max_slot<slot)
            max_slot = slot;
        }
      }
      size = max_slot;
    }
    break;
  }
  {
    int grid_max = SettingGetGlobal_i(G,cSetting_grid_max);
    if(grid_max>=0) 
      if(size>grid_max)
        size = grid_max;
  }
  return size;
}

int SceneHasImage(PyMOLGlobals *G)
{
  CScene *I=G->Scene;
  return(I->Image && I->Image->data);
}
int SceneMustDrawBoth(PyMOLGlobals *G)
{
  CScene *I=G->Scene;
  return (G->StereoCapable &&
          ((I->StereoMode==1) || 
           SettingGetGlobal_b(G,cSetting_stereo_double_pump_mono)));
}

static int SceneDeferClickWhen(Block *block, int button, int x, int y, double when,int mod);

static int stereo_via_adjacent_array(int stereo_mode)
{
  switch(stereo_mode) {
  case cStereo_crosseye:
  case cStereo_walleye:
  case cStereo_sidebyside:
    return true;
  }
  return false;
}

static int stereo_via_stencil(int stereo_mode)
{
  switch(stereo_mode) {
  case cStereo_stencil_by_row:
  case cStereo_stencil_by_column:
  case cStereo_stencil_checkerboard:
  case cStereo_stencil_custom:
    return true;
  }
  return false;
}

static int get_stereo_x(int x, int *last_x, int width, int *click_side)
{
  int width_2 = width/2;
  int width_3 = width/3;
  if(!last_x) {
    if(x>width_2) {
      x-=width_2;
      if(click_side) 
        *click_side = 1; /* right */
    } else {
      if(click_side)
        *click_side = -1; /* left */
    }
  } else {
    if((x-(*last_x))>width_3) {
      x-=width_2;
      if(click_side) *click_side = 1; /* right */
    } else if(((*last_x)-x)>width_3) {
      x+=width_2;
      if(click_side) *click_side = 1; /* right */
    } else {
      if(click_side) *click_side = -1; /* left */
    }
  }
  return x;
}

int SceneLoopDrag(Block *block,int x,int y,int mod);
int SceneLoopRelease(Block *block,int button,int x,int y,int mod);

int SceneLoopClick(Block *block,int button, int x,int y,int mod);

void SceneAbortAnimation(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  if(I->cur_ani_elem < I->n_ani_elem ) { /* allow user to override animation */
    I->cur_ani_elem = I->n_ani_elem;
  }
}

void ScenePrimeAnimation(PyMOLGlobals *G)
{
  if(G->HaveGUI) {
    CScene *I=G->Scene;
    UtilZeroMem(I->ani_elem,sizeof(CViewElem));
    SceneToViewElem(G,I->ani_elem,NULL);
    I->ani_elem[0].specification_level = 2;
    I->n_ani_elem = 0;
  }
}

static float SceneGetFPS(PyMOLGlobals *G)
{
  float fps = SettingGet(G,cSetting_movie_fps);
  float minTime;
  if(fps<=0.0F) {
    if(fps<0.0)
      minTime = 0.0; /* negative fps means full speed */
    else /* 0 fps means use movie_delay instead */
      minTime = SettingGet(G,cSetting_movie_delay)/1000.0;
    if(minTime>=0.0F)
      fps = 1.0F/minTime;
    else
      fps = 1000.0F;
  }
  return fps;
}

static void ScenePurgeImage(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  if(I->MovieOwnsImageFlag) {
    I->MovieOwnsImageFlag=false;
    I->Image=NULL;
  } else {
    if(I->Image) {
      FreeP(I->Image->data);
    }
    FreeP(I->Image);
  }
  I->CopyType = false;
}

void SceneInvalidateCopy(PyMOLGlobals *G,int free_buffer)
{
  register CScene *I=G->Scene;
  if(I) {
    if(I->MovieOwnsImageFlag) {
      I->MovieOwnsImageFlag=false;
      I->Image=NULL;
    } else if(free_buffer) {
      ScenePurgeImage(G);
    }
    I->CopyType=false;
  }
}

void SceneInvalidate(PyMOLGlobals *G)
{
   SceneInvalidateCopy(G,false);
   SceneDirty(G);
}

void SceneLoadAnimation(PyMOLGlobals *G, double duration,int hand)
{
  if(G->HaveGUI) {
    double now;
    int target = (int)(duration * 30);
    CScene *I=G->Scene;
    if(target<1)
      target = 1;
    if(target>MAX_ANI_ELEM)
      target = MAX_ANI_ELEM;
    UtilZeroMem(I->ani_elem+1,sizeof(CViewElem)*target);
    SceneToViewElem(G,I->ani_elem + target,NULL);
    I->ani_elem[target].specification_level = 2;
    now = UtilGetSeconds(G);
    I->ani_elem[0].timing_flag = true;
    I->ani_elem[0].timing = now + 0.01;
    I->ani_elem[target].timing_flag = true;
    I->ani_elem[target].timing = now + duration;
    ViewElemInterpolate(G,I->ani_elem, I->ani_elem + target, 
                        2.0F, 1.0F, true, 0.0F, hand, 0.0F);
    SceneFromViewElem(G,I->ani_elem,true);
    I->cur_ani_elem = 0;
    I->n_ani_elem = target;
    I->AnimationStartTime = UtilGetSeconds(G);
    I->AnimationStartFlag = true;
    I->AnimationStartFrame = SceneGetFrame(G);
    I->AnimationLagTime = 0.0;
  }
}

int SceneLoopClick(Block *block,int button, int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  I->LoopRect.left=x;
  I->LoopRect.top=y;
  I->LoopRect.right=x;
  I->LoopRect.bottom=y;
  I->LoopFlag=true;
  I->LoopMod = mod;
  OrthoSetLoopRect(G,true,&I->LoopRect);
  OrthoGrab(G,block);
  return 1;
}
int SceneLoopDrag(Block *block,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  I->LoopRect.right=x;
  I->LoopRect.bottom=y;
  OrthoSetLoopRect(G,true,&I->LoopRect);
  return 1;
}
int SceneLoopRelease(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  int tmp;
  int mode;
  mode = ButModeTranslate(G,button,I->LoopMod);

  if(I->LoopRect.top<I->LoopRect.bottom) {
    tmp=I->LoopRect.top;
    I->LoopRect.top=I->LoopRect.bottom;
    I->LoopRect.bottom=tmp;
  }
  if(I->LoopRect.right<I->LoopRect.left) {
    tmp=I->LoopRect.right;
    I->LoopRect.right=I->LoopRect.left;
    I->LoopRect.left=tmp;
  }
  OrthoSetLoopRect(G,false,&I->LoopRect);
  ExecutiveSelectRect(G,&I->LoopRect,mode);
  I->LoopFlag=false;
  OrthoUngrab(G);
  OrthoDirty(G);
  return 1;
}

static float GetFrontSafe(float front,float back)
{
  if(front>R_SMALL4) {
    if((back/front)>100.0F)
      front = back/100.0F;
  }
  if(front>back)
    front = back;
  if(front<cSliceMin)
    front = cSliceMin;
  return front;
}

static float GetBackSafe(float front_safe, float back)
{
  if((back-front_safe)<cSliceMin)
    back=front_safe+cSliceMin;
  return back;
}

#define SELE_MODE_MAX 7

static const char SelModeKW[][20]  = {
  "",
  "byresi",
  "bychain",
  "bysegi",
  "byobject",
  "bymol",
  "bca.",
};

static void SceneUpdateInvMatrix(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  register float *rm=I->RotMatrix;
  register float *im=I->InvMatrix;
  im[0] = rm[0];
  im[1] = rm[4];
  im[2] = rm[8];
  im[3] = 0.0F;
  im[4] = rm[1];
  im[5] = rm[5];
  im[6] = rm[9];
  im[7] = 0.0F;
  im[8] = rm[2];
  im[9] = rm[6];
  im[10] = rm[10];
  im[11] = 0.0F;
  im[12] = 0.0F;
  im[13] = 0.0F;
  im[14] = 0.0F;
  im[15] = 1.0F;
}

void SceneUpdateStereo(PyMOLGlobals *G)
{
  SceneSetStereo(G,SettingGetGlobal_b(G,cSetting_stereo));
}

char *SceneGetSeleModeKeyword(PyMOLGlobals *G)
{
  int sel_mode = SettingGetGlobal_i(G,cSetting_mouse_selection_mode);
  if((sel_mode>=0)&&(sel_mode<SELE_MODE_MAX))
    return (char*)SelModeKW[sel_mode];
  return (char*)SelModeKW[0];
}

static void SceneCopy(PyMOLGlobals *G,GLenum buffer,int force,int entire_window);

unsigned int SceneFindTriplet(PyMOLGlobals *G,int x,int y,GLenum gl_buffer);
unsigned int *SceneReadTriplets(PyMOLGlobals *G,int x,int y,int w,int h,GLenum gl_buffer);

void SceneDraw(Block *block);
void ScenePrepareMatrix(PyMOLGlobals *G,int mode);


#if 0
static int SceneGetObjState(PyMOLGlobals *G,CObject *obj,int state)
{
  int objState;
  if(SettingGetIfDefined_i(G,obj->Setting,cSetting_state,&objState)) {
    if(objState>0) { /* specific state */
      state=objState-1;
    } if(objState<0) { /* all states */
      state=-1;
    }
  }
  if(state>=0) { /* if all states for object is set */
    if(SettingGet_i(G,obj->Setting,NULL,cSetting_all_states))
      state=-1;
  }
  return(state);
}
#endif

float SceneGetDynamicLineWidth(RenderInfo *info, float line_width)
{
  if(info && info->dynamic_width) {
    float factor;
    if(info->vertex_scale>R_SMALL4) {
      factor = info->dynamic_width_factor / info->vertex_scale;
      if(factor > info->dynamic_width_max)
        factor = info->dynamic_width_max;
      if(factor < info->dynamic_width_min) {
        factor = info->dynamic_width_min;
      }
    } else {
      factor = info->dynamic_width_max;
    }
    return factor * line_width;
  }
  return line_width;
}

void SceneToViewElem(PyMOLGlobals *G,CViewElem *elem,char *scene_name)
{
  float *fp;
  double *dp;
  register CScene *I=G->Scene;

  /* copy rotation matrix */
  elem->matrix_flag = true;
  dp = elem->matrix;
  fp = I->RotMatrix;
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);

  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);

  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);

  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);

  /* copy position */
  elem->pre_flag = true;
  dp = elem->pre;
  fp = I->Pos;
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);
  *(dp++) = (double) *(fp++);

  /* copy origin (negative) */
  elem->post_flag = true;
  dp = elem->post;
  fp = I->Origin;
  *(dp++) = (double) (-*(fp++));
  *(dp++) = (double) (-*(fp++));
  *(dp++) = (double) (-*(fp++));

  elem->clip_flag = true;
  elem->front = I->Front;
  elem->back = I->Back;

  elem->ortho_flag = true;
  elem->ortho = SettingGet(G,cSetting_ortho) ? SettingGet(G,cSetting_field_of_view) : 
    -SettingGet(G,cSetting_field_of_view);

  {
    if(elem->scene_flag && elem->scene_name) {
      OVLexicon_DecRef(G->Lexicon,elem->scene_name);
      elem->scene_name = 0;
      elem->scene_flag = 0;
    }
  }
  {
    if(!scene_name)
      scene_name = SettingGetGlobal_s(G,cSetting_scene_current_name);
    if(scene_name && scene_name[0]) {
      OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon,scene_name);
      if(OVreturn_IS_OK(result)) {
        elem->scene_name = result.word;
        elem->scene_flag = true;
      }
    }
  }
}

void SceneFromViewElem(PyMOLGlobals *G,CViewElem *elem,int dirty)
{
  register CScene *I=G->Scene;
  float *fp;
  double *dp;
  int changed_flag = false;

  if(elem->matrix_flag) {
    dp = elem->matrix;
    fp = I->RotMatrix;

    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);

    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);

    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);

    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    changed_flag = true;
    SceneUpdateInvMatrix(G);
  }

  if(elem->pre_flag) {
    dp = elem->pre;
    fp = I->Pos;
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    *(fp++) = (float) *(dp++);
    changed_flag = true;
  }

  if(elem->post_flag) {
    dp = elem->post;
    fp = I->Origin;
    *(fp++) = (float) (-*(dp++));
    *(fp++) = (float) (-*(dp++));
    *(fp++) = (float) (-*(dp++));
    changed_flag = true;
  }   


  if(elem->clip_flag) {
    SceneClipSetWithDirty(G,elem->front,elem->back,dirty);
  }
  if(elem->ortho_flag) {
    if(elem->ortho<0.0F) {
      SettingSetGlobal_b(G,cSetting_ortho,0);
      if(elem->ortho<-(1.0F-R_SMALL4)) {
        SettingSetGlobal_f(G,cSetting_field_of_view,-elem->ortho);
      }
    } else {
      SettingSetGlobal_b(G,cSetting_ortho,(elem->ortho>0.5F));
      if(elem->ortho>(1.0F+R_SMALL4)) {
        SettingSetGlobal_f(G,cSetting_field_of_view,elem->ortho);
      }
    }
  }
  if(changed_flag) {

    SceneRestartSweepTimer(G);
    I->RockFrame = 0;
    SceneRovingDirty(G);
  }
}

void SceneCleanupStereo(PyMOLGlobals *G)
{
#ifndef _PYMOL_NOPY
  register CScene *I=G->Scene;  
  if(I->StereoMode==1)
    PSGIStereo(G,0);
#endif
}

static void ScenePrepareUnitContext(SceneUnitContext *context,int width,int height)
{
  float tw = 1.0F;
  float th = 1.0F;
  float aspRat;

  if(height) {
    aspRat = width/(float)height;
  } else {
    aspRat = 1.0F;
  }

  if(aspRat>1.0F) {
    tw = aspRat;
  } else {
    th = 1.0F/aspRat;
  }

  context->unit_left = (1.0F-tw)/2;
  context->unit_right = (tw+1.0F)/2;
  context->unit_top = (1.0F-th)/2;
  context->unit_bottom = (th+1.0F)/2;
  context->unit_front = -0.5F;
  context->unit_back = 0.5F;
  /*
  printf(
    "ScenePrepareUnitContext:%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
    context->unit_left,
    context->unit_right,
    context->unit_top, 
    context->unit_bottom,
    context->unit_front,
    context->unit_back);
  */
}
void SceneGetWidthHeight(PyMOLGlobals *G,int *width,int *height)
{
  register CScene *I=G->Scene;  
  *width = I->Width;
  *height = I->Height;
}

void SceneSetCardInfo(PyMOLGlobals *G,char *vendor,char *renderer,char *version){
  register CScene *I=G->Scene;  
  UtilNCopy(I->vendor,vendor,sizeof(OrthoLineType)-1);
  UtilNCopy(I->renderer,renderer,sizeof(OrthoLineType)-1);
  UtilNCopy(I->version,version,sizeof(OrthoLineType)-1);
}

int SceneGetStereo(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;  
  return(I->StereoMode);
}
void SceneGetCardInfo(PyMOLGlobals *G,char **vendor,char **renderer,char **version)
{
  register CScene *I=G->Scene;  
  (*vendor)=I->vendor;
  (*renderer)=I->renderer;
  (*version)=I->version;
}
void SceneSuppressMovieFrame(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;  
  I->MovieFrameFlag = false;
}
void SceneGetPos(PyMOLGlobals *G,float *pos)
{
  /* returns realspace coordinates for center of screen */

  register CScene *I=G->Scene;  
  
  PRINTFD(G,FB_Scene)
    " SceneGetPos: origin of rotation"
    ENDFD3f(I->Origin);
  /* take origin into camera coords */

  MatrixTransformC44fAs33f3f(I->RotMatrix,I->Origin,pos); 

  PRINTFD(G,FB_Scene)
    " SceneGetPos: origin in camera  "
    ENDFD3f(pos);

  /* find offset in camera coordinates */

  pos[0]=pos[0]-I->Pos[0]; 
  pos[1]=pos[1]-I->Pos[1];
  /*  pos[2]=pos[2]+I->Pos[2]; +(I->FrontSafe+I->Back)/2; */
  PRINTFD(G,FB_Scene)
    " SceneGetPos: center in camera  "
    ENDFD3f(pos);

  /* convert back to real coordinates */

  MatrixInvTransformC44fAs33f3f(I->RotMatrix,pos,pos);

  PRINTFD(G,FB_Scene)
    " SceneGetPos: center            "
    ENDFD3f(pos);

}
/*========================================================================*/
void SceneApplyRotMatrix(PyMOLGlobals *G,float *src,float *dst)
{
  register CScene *I=G->Scene;
  MatrixTransformC44fAs33f3f(I->RotMatrix,src,dst);
}
/*========================================================================*/
int SceneMultipick(PyMOLGlobals *G,Multipick *smp)
{
  register CScene *I=G->Scene;
  int click_side = 0;
  int defer_builds_mode = SettingGetGlobal_b(G,cSetting_defer_builds_mode);

  if(defer_builds_mode==5) /* force generation of a pickable version */
    SceneUpdate(G,true);

  if(OrthoGetOverlayStatus(G)||SettingGetGlobal_i(G,cSetting_text))
    SceneRender(G,NULL,0,0,NULL,0,0,0,0); /* remove overlay if present */
  SceneDontCopyNext(G);
  if(stereo_via_adjacent_array(I->StereoMode)) {
    if(smp->x > (I->Width/2))
      click_side = 1;
    else
      click_side = -1;
    smp->x = smp->x % (I->Width/2);
  }
  SceneRender(G,NULL,0,0,smp,0,0,click_side,0);
  SceneDirty(G);
  return(1);
}
/*========================================================================*/
int SceneGetNFrame(PyMOLGlobals *G,int *has_movie)
{
  register CScene *I=G->Scene;
  if(has_movie)
    *has_movie = I->HasMovie;
  return(I->NFrame);
}
/*========================================================================*/
void SceneGetView(PyMOLGlobals *G,SceneViewType view)
{
  float *p;
  int a;
  register CScene *I=G->Scene;
  p=view;
  for(a=0;a<16;a++)
    *(p++) = I->RotMatrix[a];
  *(p++) = I->Pos[0];
  *(p++) = I->Pos[1];
  *(p++) = I->Pos[2];
  *(p++) = I->Origin[0];
  *(p++) = I->Origin[1];
  *(p++) = I->Origin[2];
  *(p++) = I->Front;
  *(p++) = I->Back;
  *(p++) = SettingGet(G,cSetting_ortho) ? SettingGet(G,cSetting_field_of_view) : 
    -SettingGet(G,cSetting_field_of_view);

}
/*========================================================================*/
void SceneSetView(PyMOLGlobals *G,SceneViewType view,
                  int quiet,float animate,int hand)
{
  float *p;
  int a;
  register CScene *I=G->Scene;

  if(animate<0.0F) {
    if(SettingGetGlobal_b(G,cSetting_animation))
      animate=SettingGetGlobal_f(G,cSetting_animation_duration);
    else
      animate=0.0F;
  }
  if(animate!=0.0F)
    ScenePrimeAnimation(G);
  else {
    SceneAbortAnimation(G);
  }

  p=view;
  for(a=0;a<16;a++)
    I->RotMatrix[a] = *(p++); 
  SceneUpdateInvMatrix(G);
  I->Pos[0] = *(p++);
  I->Pos[1] = *(p++);
  I->Pos[2] = *(p++);
  I->Origin[0] = *(p++);
  I->Origin[1] = *(p++);
  I->Origin[2] = *(p++);

  I->LastSweep = 0.0F;
  I->LastSweepX = 0.0F;
  I->LastSweepY = 0.0F;
  I->SweepTime = 0.0;
  I->RockFrame = 0;

  SceneClipSet(G,p[0],p[1]);
  p+=2;
  if(p[0]<0.0F) {
    SettingSetGlobal_b(G,cSetting_ortho,0);
    if(p[0]<-(1.0F-R_SMALL4)) {
      SettingSetGlobal_f(G,cSetting_field_of_view,-p[0]);
    }
  } else {
    SettingSetGlobal_b(G,cSetting_ortho,(p[0]>0.5F));
    if(p[0]>(1.0F+R_SMALL4)) {
      SettingSetGlobal_f(G,cSetting_field_of_view,p[0]);
    }
  }
  if(!quiet) { 
    PRINTFB(G,FB_Scene,FB_Actions)
      " Scene: view updated.\n"
      ENDFB(G);
  }
  if(animate!=0.0F)
    SceneLoadAnimation(G,animate,hand);

  SceneRovingDirty(G);
}
/*========================================================================*/
void SceneDontCopyNext(PyMOLGlobals *G)
/* disables automatic copying of the image for the next rendering run */
{
  register CScene *I=G->Scene;
  I->CopyNextFlag=false;
}
/*========================================================================*/
void SceneUpdateStereoMode(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  if(I->StereoMode) {
    SceneSetStereo(G,true);
  }
}
/*========================================================================*/
void SceneSetStereo(PyMOLGlobals *G,int flag)
{
  register CScene *I=G->Scene;
  int cur_stereo = I->StereoMode;

  if(flag) {
    I->StereoMode=(int)SettingGet(G,cSetting_stereo_mode); 
  } else {
    I->StereoMode=0;
  }

  if((cur_stereo!=I->StereoMode)&&((cur_stereo==4)||(I->StereoMode==4))) {
    OrthoReshape(G,G->Option->winX,G->Option->winY,true);
#ifndef _PYMOL_NOPY
    if((cur_stereo==4)&&I->StereoMode) {
      PParse(G,"viewport");
    }
#endif
  }
  SettingSetGlobal_b(G,cSetting_stereo,flag);
  SceneInvalidateStencil(G);
  SceneInvalidate(G);
}
/*========================================================================*/
void SceneTranslate(PyMOLGlobals *G,float x,float y, float z)
{
  register CScene *I=G->Scene;
  I->Pos[0]+=x;
  I->Pos[1]+=y;
  I->Pos[2]+=z;
  I->Back-=z;
  I->Front-=z;
  if(I->Front>I->Back)
    I->Front=I->Back+cSliceMin;
  I->FrontSafe= GetFrontSafe(I->Front,I->Back);
  I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
  SceneInvalidate(G);
}

void SceneTranslateScaled(PyMOLGlobals *G,float x,float y, float z, int sdof_mode)
{
  register CScene *I=G->Scene;
  int invalidate = false;
  
  switch(sdof_mode) {
  case SDOF_NORMAL_MODE:
    if((x!=0.0F)||(y!=0.0F)) {
      float vScale = SceneGetExactScreenVertexScale(G,NULL);
      float factor = vScale * (I->Height + I->Width) / 2;
      I->Pos[0]+=x * factor;
      I->Pos[1]+=y * factor;
      invalidate = true;
    }
    if(z!=0.0F) {
      float factor = ((I->FrontSafe+I->BackSafe)/2); /* average distance within visible space */
      if(factor>0.0F) {
        factor *= z;
        I->Pos[2] += factor;
        I->Front -= factor;
        I->Back -= factor;
        I->FrontSafe = GetFrontSafe(I->Front,I->Back);
        I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
      }
      invalidate = true;
    }
    break;
  case SDOF_CLIP_MODE:
    if((x!=0.0F)||(y!=0.0F)) {
      float vScale = SceneGetExactScreenVertexScale(G,NULL);
      float factor = vScale * (I->Height + I->Width) / 2;
      I->Pos[0]+=x * factor;
      I->Pos[1]+=y * factor;
      invalidate = true;
    }
    if(z!=0.0F) {
      float factor = ((I->FrontSafe+I->BackSafe)/2); /* average distance within visible space */
      if(factor>0.0F) {
        factor *= z;
        {
          float old_front = I->Front;
          float old_back = I->Back;
          float old_origin = -I->Pos[2];
          SceneClip(G,7,factor,NULL,0);
          SceneDoRoving(G,old_front,old_back,old_origin,true,true);
        }
        invalidate = true;
      }
    }
    break;
  case SDOF_DRAG_MODE:
    {
      float v2[3];
      float scale = SettingGetGlobal_f(G,cSetting_sdof_drag_scale);
      
      { 
        /* when dragging, we treat all axes proportionately */
        float vScale = SceneGetExactScreenVertexScale(G,NULL);
        float factor = vScale * (I->Height + I->Width) / 2;
        x *= factor;
        y *= factor;
        z *= factor;
      }

      v2[0] = x * scale;
      v2[1] = y * scale;
      v2[2] = z * scale;
      
      /* transform into model coodinate space */
      MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); 
      
      EditorDrag(G,NULL,-1,cButModeMovDrag,
                 SettingGetGlobal_i(G,cSetting_state)-1,NULL,v2,NULL);
    }
    break;
  }
  if(invalidate) {
    SceneInvalidate(G);
    if(SettingGetGlobal_b(G,cSetting_roving_origin)) {
      float v2[3];
      SceneGetPos(G,v2); /* gets position of center of screen */
      SceneOriginSet(G,v2,true);
    }
    if(SettingGetGlobal_b(G,cSetting_roving_detail)) {    
      SceneRovingPostpone(G);
    }
  }
}

void SceneRotateScaled(PyMOLGlobals *G,float rx,float ry, float rz, int sdof_mode)
{
  register CScene *I=G->Scene; 
  int invalidate = false;
  float axis[3];
  switch(sdof_mode) {
  case SDOF_NORMAL_MODE:
    axis[0]=rx;
    axis[1]=ry;
    axis[2]=rz;
    {
      float angle = length3f(axis);
      normalize3f(axis);
      SceneRotate(G,60*angle,axis[0],axis[1],axis[2]);
    }
    break;
  case SDOF_CLIP_MODE:
    if( (fabs(rz)>fabs(rx)) || (fabs(rz)>fabs(rx))) {
      rx = 0.0;
      ry = 0.0;
    } else {
      rz = 0.0;
    }
    axis[0]=rx;
    axis[1]=ry;
    axis[2]=0.0;
    {
      float angle = length3f(axis);
      normalize3f(axis);
      SceneRotate(G,60*angle,axis[0],axis[1],axis[2]);
    }
    if(axis[2]!=rz) {
      SceneClip(G,5,1.0F+rz,NULL,0);
    }
    break;
  case SDOF_DRAG_MODE:
    {
      float scale = SettingGetGlobal_f(G,cSetting_sdof_drag_scale);
      float v1[3],v2[3];
      axis[0]=rx;
      axis[1]=ry;
      axis[2]=rz;

      EditorReadyDrag(G,SettingGetGlobal_i(G,cSetting_state)-1);
    
      {
        float angle = length3f(axis);
        normalize3f(axis);
        
        v1[0] = cPI*(60*angle/180.0F) * scale;
        
        /* transform into model coodinate space */
        MatrixInvTransformC44fAs33f3f(I->RotMatrix,axis,v2); 
        
        EditorDrag(G,NULL,-1,cButModeRotDrag,
                   SettingGetGlobal_i(G,cSetting_state)-1,v1,v2,NULL);
        invalidate = true;
      }
    }
    break;
  }
  if(invalidate) {
    SceneInvalidate(G);
  }
}
/*========================================================================*/

static void SceneClipSetWithDirty(PyMOLGlobals *G,float front,float back,int dirty)
{
  register CScene *I=G->Scene;
  I->Front=front;
  I->Back=back;
  if(I->Front>I->Back)
  I->Front=I->Back+cSliceMin;
  I->FrontSafe= GetFrontSafe(I->Front,I->Back);
  I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
  if(dirty) 
    SceneInvalidate(G);
  else
    SceneInvalidateCopy(G,false);
}
/*========================================================================*/
void SceneClipSet(PyMOLGlobals *G,float front,float back)
{
  SceneClipSetWithDirty(G,front,back,true);
}
/*========================================================================*/
void SceneClip(PyMOLGlobals *G,int plane,float movement,char *sele,int state) /* 0=front, 1=back*/
{
  register CScene *I=G->Scene;
  float avg;
  float mn[3],mx[3],cent[3],v0[3],offset[3],origin[3];
  switch(plane) {
  case 0: /* near */
    SceneClipSet(G,I->Front-movement,I->Back);
    break;
  case 1: /* far */
    SceneClipSet(G,I->Front,I->Back-movement);      
    break;
  case 2: /* move */
    SceneClipSet(G,I->Front-movement,I->Back-movement);    
    break;
  case 3: /* slab */
    if(sele[0]) {
      if(!ExecutiveGetExtent(G,sele,mn,mx,true,state,false))
        sele = NULL;
      else {
        average3f(mn,mx,cent); /* get center of selection */
        subtract3f(cent,I->Origin,v0); /* how far from origin? */
        MatrixTransformC44fAs33f3f(I->RotMatrix,v0,offset); /* convert to view-space */
      }
    } else {
      sele = NULL;
    }
    avg = (I->Front+I->Back)/2.0F;
    movement/=2.0F;
    if(sele) {
      avg = -I->Pos[2]-offset[2];
    }
    SceneClipSet(G,avg-movement,avg+movement);
    break;
  case 4: /* atoms */
    if(!sele) 
      sele=cKeywordAll;
    else if(!sele[0]) {
      sele=cKeywordAll;
    } 
    if(WordMatchExact(G,sele,cKeywordCenter,true)) {
      MatrixTransformC44fAs33f3f(I->RotMatrix,I->Origin,origin); /* convert to view-space */
      SceneClipSet(G,origin[2]-movement,origin[2]+movement);
    } else if(WordMatchExact(G,sele,cKeywordOrigin,true)) {
      SceneClipSet(G,-I->Pos[2]-movement,-I->Pos[2]+movement);      
    } else {
      if(!ExecutiveGetCameraExtent(G,sele,mn,mx,true,state))
        sele = NULL;
      if(sele) {
        if(sele[0]) {
          average3f(mn,mx,cent); /* get center of selection */
          MatrixTransformC44fAs33f3f(I->RotMatrix,I->Origin,origin); /* convert to view-space */
          subtract3f(mx,origin,mx); /* how far from origin? */
          subtract3f(mn,origin,mn); /* how far from origin? */
          SceneClipSet(G,-I->Pos[2]-mx[2]-movement,-I->Pos[2]-mn[2]+movement);
        } else {
          sele = NULL;
        }
      }
    }
    break;
  case 5: /* scaling */
    {
      float width = I->Front-I->Back;
      float new_width = movement * width;
      float avg = (I->Front+I->Back)/2;
      float new_front = avg+new_width/2.0F;
      float new_back = avg-new_width/2.0F;
      
      SceneClipSet(G,new_front,new_back);
    }
    break;
  case 6: /* proportional movement */
    {
      float shift = (I->Front-I->Back)*movement;
      SceneClipSet(G,I->Front+shift,I->Back+shift);    
    }
    break;
  case 7: /* linear movement */
    {
      SceneClipSet(G,I->Front+movement,I->Back+movement);    
    }
    break;
  }
}
/*========================================================================*/
void SceneSetMatrix(PyMOLGlobals *G,float *m)
{
  register CScene *I=G->Scene;
  int a;
  for(a=0;a<16;a++)
     I->RotMatrix[a]=m[a];
  SceneUpdateInvMatrix(G);
}
/*========================================================================*/
void SceneGetViewNormal(PyMOLGlobals *G,float *v)
{
  register CScene *I=G->Scene;
  copy3f(I->ViewNormal,v);
}
/*========================================================================*/
int SceneGetState(PyMOLGlobals *G)
{
  return(SettingGetGlobal_i(G,cSetting_state)-1);
}
/*========================================================================*/
float *SceneGetMatrix(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  return(I->RotMatrix);
}
/*========================================================================*/

int SceneCaptureWindow(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  int ok=true;
  
  /* check assumptions */
  
  if(ok && G->HaveGUI && G->ValidContext) {
    int draw_both = SceneMustDrawBoth(G);

    ScenePurgeImage(G);

    if(draw_both) {
      SceneCopy(G,GL_BACK_LEFT,true,true);
    } else {
      SceneCopy(G,GL_BACK,true,true);
    }
    if(!I->Image)
        ok=false;
    
    if(ok && I->Image) {
      I->DirtyFlag = false;
      I->CopyType = 2; /* suppresses display of copied image */
      if(SettingGetGlobal_b(G,cSetting_opaque_background)) 
        I->Image->needs_alpha_reset = true; 
      I->MovieOwnsImageFlag = false;
    }
  } else {
    ok=false;
  }
  return ok;
}

static int SceneMakeSizedImage(PyMOLGlobals *G,int width,
                               int height, int antialias)
{
  register CScene *I=G->Scene;
  int ok=true;
  int save_flag = false;
  int save_width=0, save_height=0;

  /* check assumptions */
    
  if( (width && height && I->Width && I->Height ) &&
      fabs(((float)(height - (width * I->Height) / (I->Width)))/height) > 0.01F)
    {
      save_width = I->Width;
      save_height = I->Height;
      save_flag = true;

      /* squish the dimensions as needed to maintain 
         aspect ratio within the current rectangle */

      if(I->Width > ( (width * I->Height) / height))
        I->Width = (width * I->Height) / height;
      else if(I->Height > ( (height * I->Width) / width))
        I->Height = ( height * I->Width) / width;
    }
  if( (!width) && (!height) ) {
    width=I->Width;
    height=I->Height;
  }
  if(width && !height) {
    height = (I->Height * width) / I->Width;
  }
  if(height && !width) {
    width = (I->Width * height) / I->Height;
  }
  if(!((width>0)&&
       (height>0)&&
       (I->Width>0)&&
       (I->Height>0))) {
    PRINTFB(G,FB_Scene,FB_Errors)
      "SceneMakeSizedImage-Error: invalid image dimensions\n"
      ENDFB(G);
    ok=false;
  }

  if(ok && G->HaveGUI && G->ValidContext) {

    int factor = 0;
    int shift = 0;
    int max_dim[2];

    glGetIntegerv(GL_MAX_VIEWPORT_DIMS,(GLint*)(void*)max_dim);

    /* clamp to what this OpenGL implementation can do */
    if(width>max_dim[0]) {
      height = (max_dim[0] * height) / width;
      width = max_dim[0];
      PRINTFB(G,FB_Scene,FB_Warnings)
        "Scene-Warning: Maximum OpenGL viewport dimension exceeded."
        ENDFB(G)
    }
    if(height>max_dim[0]) {
      width = (max_dim[1] * width ) / height;
      height = max_dim[1];
      PRINTFB(G,FB_Scene,FB_Warnings)
        "Scene-Warning: Maximum OpenGL viewport dimension exceeded."
        ENDFB(G)
      
    }

    if(antialias==1) {
      factor = 2;
      shift = 2; 
    } if(antialias>=2) {
      factor = 4;
      shift = 4;
    }

    while(factor>1) {
      if(((width*factor)<max_dim[0])&&
         ((height*factor)<max_dim[1])) {
        width = width * factor;
        height = height * factor;
        break;
      } else {
        factor = (factor>>1);
        shift = shift - 2;
        if(factor<2) {
          PRINTFB(G,FB_Scene,FB_Blather)
            "Scene-Warning: Maximum OpenGL viewport exceeded. Antialiasing disabled."
            ENDFB(G);
          break;
        }
      }
    }
    if(factor<2)
      factor = 0;

    {
      unsigned int final_buffer_size = width*height;
      unsigned int *final_image = NULL;
      int nXStep = (width/(I->Width+1)) + 1;
      int nYStep = (height/(I->Height+1)) + 1;
      int x,y;
      int draw_both = SceneMustDrawBoth(G);
      /* note here we're treating the buffer as 32-bit unsigned ints, not chars */
      
      final_image = Alloc(unsigned int,final_buffer_size);
      
      if(!final_image) { 
        ok=false;
      }
      ScenePurgeImage(G);

      if(draw_both) {
        SceneCopy(G,GL_BACK_LEFT,true,false);
      } else {
        SceneCopy(G,GL_BACK,true,false);
      }
      if(!I->Image)
        ok=false;
      
      if(ok) {
        int total_steps = nXStep * nYStep;
        
        OrthoBusyPrime(G);
        
        /* so the trick here is that we need to move the camera around
           so that we get a pixel-perfect mosaic */
        for(y=0;y<nYStep;y++) {
          int y_offset = -(I->Height*y);
          
          for(x=0;x<nXStep;x++) {
            int x_offset = -(I->Width*x);
            int a,b;
            float *v;  
            float alpha = (SettingGetGlobal_b(G,cSetting_opaque_background) ? 1.0F : 0.0F);
            unsigned int *p, *q, *qq, *pp;
            v=SettingGetfv(G,cSetting_bg_rgb);
            OrthoBusyFast(G,y*nXStep+x,total_steps);

            if(draw_both) {
              OrthoDrawBuffer(G,GL_BACK_LEFT);
            } else {
              OrthoDrawBuffer(G,GL_BACK);
            }
            
            glClearColor(v[0],v[1],v[2],alpha);
            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
            SceneInvalidateCopy(G,false);
            SceneRender(G,NULL,x_offset,y_offset,NULL,width,height,0,0);
            glClearColor(0.0,0.0,0.0,1.0);
            
            if(draw_both) {
              SceneCopy(G,GL_BACK_LEFT,true,false);
            } else {
              SceneCopy(G,GL_BACK,true,false);
            }
            
            if(I->Image) { /* the image into place */
              p = (unsigned int*)I->Image->data;
              q = final_image + (x*I->Width) + (y*I->Height)*width;
              {
                int y_limit;
                int x_limit;
                
                if(((y+1)*I->Height)>height)
                  y_limit = height - (y*I->Height);
                else
                  y_limit = I->Height;
                if(((x+1)*I->Width)>width)
                  x_limit = width - (x*I->Width);
                else
                  x_limit = I->Width;
                for(a=0;a<y_limit;a++) {
                  qq = q;
                  pp = p;
                  for(b=0;b<x_limit;b++) {
                    *(qq++) = *(pp++);
                  }
                  q += width;
                  p += I->Width;
                }
              }
            }
          }
        }


        if(!OrthoDeferredWaiting(G)) {
          if(SettingGetGlobal_i(G,cSetting_draw_mode)==-2) {
            ExecutiveSetSettingFromString(G,cSetting_draw_mode,"-1","",-1, true, true);
            SceneUpdate(G,false);
          }
        }

        SceneInvalidateCopy(G,true);
        
        if(factor) { /* are we oversampling? */
          unsigned int src_row_bytes = width * 4;
          
          width = width / factor;
          height = height / factor;
          
          {
            unsigned char *p = (unsigned char*)final_image;
            unsigned char *buffer = Alloc(unsigned char, 4 * width * height);
            unsigned char *q = buffer;
            register unsigned char *pp, *ppp, *pppp;
            register int a,b,c,d;
            register unsigned int c1,c2,c3,c4,alpha;
            unsigned int factor_col_bytes = factor * 4;
            unsigned int factor_row_bytes = factor * src_row_bytes;
            
            for(b=0;b<height;b++) { /* rows */
              pp = p;
              for(a=0;a<width;a++) { /* cols */
                c1 = c2 = c3 = c4 = 0;
                ppp = pp;
                for(d=0;d<factor;d++) { /* box rows */
                  pppp = ppp;
                  for(c=0;c<factor;c++) { /* box cols */
                    c4 += (alpha = pppp[3]);
                    c1 += *(pppp++) * alpha;
                    c2 += *(pppp++) * alpha;
                    c3 += *(pppp++) * alpha;
                    pppp++;
                  }
                  ppp += src_row_bytes;
                }
                if(c4) { /* divide out alpha channel & average */
                  c1 = c1 / c4;
                  c2 = c2 / c4;
                  c3 = c3 / c4;
                } else { /* alpha zero! so compute average RGB */
                  c1 = c2 = c3 = 0;
                  ppp = pp;
                  for(d=0;d<factor;d++) { /* box rows */
                    pppp = ppp;
                    for(c=0;c<factor;c++) { /* box cols */
                      c1 += *(pppp++);
                      c2 += *(pppp++);
                      c3 += *(pppp++);
                      pppp++;
                    }
                    ppp += src_row_bytes;
                  }
                  c1 = c1 >> shift;
                  c2 = c2 >> shift;
                  c3 = c3 >> shift;
                }
                *(q++)= c1;
                *(q++)= c2;
                *(q++)= c3;
                *(q++)= c4>>shift;
                pp += factor_col_bytes;
              }
              p+= factor_row_bytes;
            }
            
            FreeP(final_image);
            final_image = (unsigned int*)buffer;
          }
        }
        ScenePurgeImage(G);

        I->Image = Calloc(ImageType,1);
        I->Image->data = (unsigned char*)final_image;
        final_image = NULL;
        I->Image->size = final_buffer_size*4; /* in bytes, not 32-bit words */
        I->Image->width=width;
        I->Image->height=height;
        I->Image->stereo=false;
        
        I->DirtyFlag = false;
        I->CopyType = true;
        I->CopyForced = true;

        if(SettingGetGlobal_b(G,cSetting_opaque_background)) 
          I->Image->needs_alpha_reset = true;

        I->MovieOwnsImageFlag = false;
      }
      FreeP(final_image);
    }
  } else {
    ok=false;
  }

  if(save_flag) {
    I->Width = save_width;
    I->Height = save_height;
  }
  return ok;

}

/*========================================================================*/
static unsigned char *SceneImagePrepare(PyMOLGlobals *G,int prior_only)
{
  register CScene *I=G->Scene;
  unsigned char *image = NULL;
  int reset_alpha = false;
  int save_stereo = (I->StereoMode==1);
  
  if(!(I->CopyType || prior_only)) {
    if(G->HaveGUI && G->ValidContext) {
      unsigned int buffer_size;
      
      buffer_size = 4*I->Width*I->Height;
      if(save_stereo)
        image = (GLvoid*)Alloc(char,buffer_size*2);
      else
        image = (GLvoid*)Alloc(char,buffer_size);
      ErrChkPtr(G,image);
      if(SceneMustDrawBoth(G)||save_stereo) {
        glReadBuffer(GL_BACK_LEFT);
      } else {
        glReadBuffer(GL_BACK);
      }
      PyMOLReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                      GL_RGBA,GL_UNSIGNED_BYTE,(GLvoid*)(image));
      if(save_stereo) {
        glReadBuffer(GL_BACK_RIGHT);
        PyMOLReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                        GL_RGBA,GL_UNSIGNED_BYTE,(GLvoid*)(image+buffer_size));
      }
      reset_alpha = true;
      ScenePurgeImage(G);
      I->Image=Calloc(ImageType,1);
      I->Image->data = image;
      I->Image->height=I->Height;
      I->Image->width=I->Width;
      I->Image->size = buffer_size;
      if(save_stereo)
        I->Image->stereo = 1;
    }
  } else if(I->Image) {
    image = I->Image->data;
    reset_alpha = I->Image->needs_alpha_reset;
  }
  if(image) {
    int opaque_back = SettingGetGlobal_b(G,cSetting_opaque_background);
    if(opaque_back && reset_alpha) {
      unsigned char *p = (unsigned char*)image;
      int x,y;
      int width = I->Image->width;
      int height = I->Image->height;
      if(I->Image && (I->Image->data == (unsigned char*)image))
        I->Image->needs_alpha_reset = false;
      for(y=0;y<height;y++) {
        for(x=0;x<width;x++) {
          p[3]=0xFF;
          p+=4;
        }
      }
      if(save_stereo) {
        for(y=0;y<height;y++) {
          for(x=0;x<width;x++) {
            p[3]=0xFF;
            p+=4;
          }
        }
      }
    }
  }
  return (unsigned char*)image;
}

static void SceneImageFinish(PyMOLGlobals *G,char *image)
{
  register CScene *I=G->Scene;
  if(I->Image) {
    if(I->Image->data!=(unsigned char*)image) /* purge the image if this isn't the active copy */
      FreeP(image);
  } else {
    FreeP(image);
  }
}

void SceneGetImageSize(PyMOLGlobals *G,int *width,int *height)
{
  register CScene *I=G->Scene;
  GLvoid *image = SceneImagePrepare(G,false);
  if(image && I->Image) {
    *width = I->Image->width;
    *height = I->Image->height;
  } else {
    *width = I->Width;
    *height = I->Height;
  }
  SceneImageFinish(G,image); /* don't leak if(image != I->Image) */
}

int  SceneCopyExternal(PyMOLGlobals *G,int width, int height,
               int rowbytes,unsigned char *dest,int mode)
{
  GLvoid *image = SceneImagePrepare(G,false);
  register CScene *I=G->Scene;
  int result=false;
  int i,j;
  int premultiply_alpha = true;
  int red_index=0,blue_index=1,green_index=2,alpha_index=3;
  int no_alpha = (SettingGetGlobal_b(G,cSetting_opaque_background) &&
                  SettingGetGlobal_b(G,cSetting_ray_opaque_background));

  if(mode&0x1) {
    int index=0;
    while(index<4) {
      if(dest[index]=='R') red_index = index;
      if(dest[index]=='G') green_index = index;
      if(dest[index]=='B') blue_index = index;
      if(dest[index]=='A') alpha_index = index;
      index++;
    }
  }
  if(mode&0x2) {
    premultiply_alpha = false;
  }
  /*
  printf("image %p I->image %p\n");
  if(I->Image) {
    printf("%d %d %d %d\n",I->Image->width,width,I->Image->height,height);
    }*/

  if(image&&I->Image&&(I->Image->width==width)&&(I->Image->height==height)) {
    for (i=0; i< height; i++)
      {
        unsigned char *src = ((unsigned char*)image) + ((height-1)-i) * width*4;
        unsigned char *dst;
        if(mode&0x4) {
          dst = dest + (height-(i+1))*(rowbytes);
        } else {
          dst = dest + i * (rowbytes);
        }
        for (j = 0; j < width; j++) {
          if(no_alpha) {
            dst[red_index]   = src[0]; /* no alpha */
            dst[green_index] = src[1];
            dst[blue_index]  = src[2];
            dst[alpha_index] = 0xFF;
            /*            if(!(i||j)) {
              printf("no alpha\n");
              }*/
          } else if(premultiply_alpha) {
            dst[red_index]   = (((unsigned int)src[0])*src[3])/255; /* premultiply alpha */
            dst[green_index] = (((unsigned int)src[1])*src[3])/255;
            dst[blue_index]  = (((unsigned int)src[2])*src[3])/255;
            dst[alpha_index] = src[3];
            /*       if(!(i||j)) {
              printf("premult alpha\n");
              }*/
          } else {
            dst[red_index]   = src[0]; /* standard alpha */
            dst[green_index] = src[1];
            dst[blue_index]  = src[2];
            dst[alpha_index] = src[3];
            /*            if(!(i||j)) {
              printf("standard alpha\n");
              }*/
          }
          dst+=4;
          src+=4;
        }
      }
    result=true;
  } else {
    printf("image or size mismatch\n");
  }
  SceneImageFinish(G,image);  
  return(result);
}

static void interlace(unsigned int *dst,unsigned int *src,int width,int height)
{
  register int a,b;
  unsigned int *p0 = src, *p1 = src+(height*width);
  unsigned int *q = dst;
  for(a=0;a<height;a++) {
    for(b=0;b<width;b++) {
      *(q++) = *(p0++);
    }
    for(b=0;b<width;b++) {
      *(q++) = *(p1++);      
    }
  }
}
static void deinterlace(unsigned int *dst,unsigned int *src,
                        int width,int height,int swap)
{
  register int a,b;
  unsigned int *p = src;
  unsigned int *q0 = dst, *q1 = dst+(height*width);
  if(swap) {
    q0 = dst+(height*width);
    q1 = dst;
  }
    
  for(a=0;a<height;a++) {
    for(b=0;b<width;b++) {
      *(q0++) = *(p++);
    }
    for(b=0;b<width;b++) {
      *(q1++) = *(p++);      
    }
  }
}

int ScenePNG(PyMOLGlobals *G,char *png,float dpi,int quiet,
             int prior_only,int format)
{
  register CScene *I=G->Scene;
  GLvoid *image = SceneImagePrepare(G,prior_only);
  if(image && I->Image) {
    int width = I->Image->width;
    int height = I->Image->height;
    unsigned char *save_image = image;

    if((image==I->Image->data) && I->Image->stereo) {
      width = I->Image->width;
      save_image = Alloc(unsigned char, I->Image->size*2);
      interlace((unsigned int*)save_image,(unsigned int*)I->Image->data,width,height);
      width *= 2;
    }
    if(dpi<0.0F) dpi = SettingGetGlobal_f(G,cSetting_image_dots_per_inch);
    if(MyPNGWrite(G,png,save_image,width,height,dpi,format,quiet)) {
      if(!quiet) {
        PRINTFB(G,FB_Scene,FB_Actions) 
          " ScenePNG: wrote %dx%d pixel image to file \"%s\".\n",
          width,I->Image->height,png
          ENDFB(G);
      }
    } else {
      PRINTFB(G,FB_Scene,FB_Errors) 
        " ScenePNG-Error: error writing \"%s\"! Please check directory...\n",
        png
        ENDFB(G);
    }
    if(save_image && (save_image!=image))
      FreeP(save_image);
  }
  SceneImageFinish(G,image);  
  return (image!=NULL);
}
/*========================================================================*/
void ScenePerspective(PyMOLGlobals *G,int flag)
{
  float persp;
  persp=(float)(!flag);
  SettingSetfv(G,cSetting_ortho,&persp);
  SceneInvalidate(G);
}
/*========================================================================*/
int SceneGetFrame(PyMOLGlobals *G)
{
  if(MovieDefined(G))
    return(SettingGetGlobal_i(G,cSetting_frame)-1);
  else
    return(SettingGetGlobal_i(G,cSetting_state)-1);    
}
/*========================================================================*/
void SceneCountFrames(PyMOLGlobals *G) 
{
  register CScene *I=G->Scene;
  ObjRec *rec = NULL;
  int n;
  int mov_len;
  I->NFrame=0;
  while(ListIterate(I->Obj,rec,next)) {
    if(rec->obj->fGetNFrame)
      n=rec->obj->fGetNFrame(rec->obj);
    else
      n=0;
    if(n>I->NFrame)
      I->NFrame=n;
  }
  mov_len = MovieGetLength(G);
  I->HasMovie = (mov_len != 0);
  if(mov_len>0) {
    I->NFrame=mov_len;
  } else if(mov_len<0) {
    mov_len=-mov_len;
    if(I->NFrame<mov_len) /* allows you to see cached movie even w/o object */
      I->NFrame=mov_len;
  }
  PRINTFD(G,FB_Scene)
    " SceneCountFrames: leaving... I->NFrame %d\n",I->NFrame
    ENDFD
}
/*========================================================================*/
void SceneSetFrame(PyMOLGlobals *G,int mode,int frame)
{
  register CScene *I=G->Scene;
  int newFrame;
  int newState=0;
  int movieCommand = false;
  newFrame = SettingGetGlobal_i(G,cSetting_frame) -1;
  PRINTFD(G,FB_Scene)
    " SceneSetFrame: entered.\n"
    ENDFD;
  switch(mode) {
  case -1: /* movie/frame override - go to this state absolutely! */
    newState=frame;
    break;
  case 0: /* absolute frame */
    newFrame=frame; 
     break;
  case 1: /* relative frame */
    newFrame+=frame; 
     break;
  case 2: /* end */
    newFrame=I->NFrame-1; 
     break;
  case 3: /* middle with automatic movie command */
     newFrame=I->NFrame/2;
    movieCommand = true;
     break;
  case 4: /* absolute with automatic movie command */
    newFrame=frame;
    movieCommand = true;
    break;
  case 5: /* relative with automatic movie command */
     newFrame+=frame;
     movieCommand = true;
     break;
  case 6: /* end with automatic movie command */
    newFrame=I->NFrame-1; 
    movieCommand = true;
    break;
  case 7: /* absolute with forced movie command */
    newFrame=frame;
    movieCommand = true;
    break;
  case 8: /* relative with forced movie command */
     newFrame+=frame;
    movieCommand = true;
     break;
  case 9: /* end with forced movie command */
    newFrame=I->NFrame-1; 
    movieCommand = true;
    break;
  }
  SceneCountFrames(G);
  if (mode>=0) { 
    if(newFrame>=I->NFrame) newFrame=I->NFrame-1;
    if(newFrame<0) newFrame=0;
    newState = MovieFrameToIndex(G,newFrame);
    if(newFrame==0) {
      if(MovieMatrix(G,cMovieMatrixRecall)) {
        SceneAbortAnimation(G); /* if we have a programmed initial
                                   orientation, don't allow animation
                                   to override it */
      }
    }
    SettingSetGlobal_i(G,cSetting_frame,newFrame+1);
    SettingSetGlobal_i(G,cSetting_state,newState+1);
    if(movieCommand) {
      MovieDoFrameCommand(G,newFrame);
      MovieFlushCommands(G);
    }
    if(SettingGet(G,cSetting_cache_frames))
      I->MovieFrameFlag=true;
  } else {
    SettingSetGlobal_i(G,cSetting_frame,newFrame+1);
    SettingSetGlobal_i(G,cSetting_state,newState+1);
  }

  SceneInvalidate(G);
  PRINTFD(G,FB_Scene)
    " SceneSetFrame: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void SceneDirty(PyMOLGlobals *G) 
      /* This means that the current image on the screen (and/or in the buffer)
         needs to be updated */
{
  register CScene *I=G->Scene;

  PRINTFD(G,FB_Scene)
    " SceneDirty: called.\n"
    ENDFD;

  if(I) {
    if(!I->DirtyFlag) {
      I->DirtyFlag=true;
      /* SceneInvalidateCopy(G,false); */
      OrthoDirty(G);
    }
  }

}

void SceneRovingPostpone(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  float delay;
  if(SettingGet(G,cSetting_roving_detail)) {
    delay = SettingGet(G,cSetting_roving_delay);
    if(delay<0.0F) {
      I->RovingLastUpdate = UtilGetSeconds(G); /* put off delay */
    }
  }
}

void SceneRovingDirty(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;

  if(SettingGet(G,cSetting_roving_detail)) {
    SceneRovingPostpone(G);
    I->RovingDirtyFlag=true;
  }
}

/*========================================================================*/
void SceneChanged(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  I->ChangedFlag=true;
  SceneInvalidateCopy(G,false);
  SceneDirty(G);
  SeqChanged(G);
}
/*========================================================================*/
Block *SceneGetBlock(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  return(I->Block);
}
/*========================================================================*/
int SceneMakeMovieImage(PyMOLGlobals *G,int show_timing, 
                        int validate, int mode) {
  register CScene *I=G->Scene;
  float *v;
  int valid = true;
  PRINTFB(G,FB_Scene,FB_Blather)
    " Scene: Making movie image.\n"
    ENDFB(G);

  switch(mode) {
  case cSceneImage_Normal:
  case cSceneImage_Draw:
  case cSceneImage_Ray:
    break;
  default:
    if(SettingGet(G,cSetting_ray_trace_frames)) {
      mode = cSceneImage_Ray;
    } else if(SettingGet(G,cSetting_draw_frames)) {
      mode = cSceneImage_Draw;
    } else {
      mode = cSceneImage_Normal;
    }
    break;
  }

  I->DirtyFlag=false;
  switch(mode) {
  case cSceneImage_Ray:
    SceneRay(G,0,0,(int)SettingGet(G,cSetting_ray_default_renderer),
             NULL,NULL,
             0.0F,0.0F,false,NULL,show_timing,-1); 
    break;
  case cSceneImage_Draw:
    SceneMakeSizedImage(G,0,0,SettingGetGlobal_i(G,cSetting_antialias));
    break;
  case cSceneImage_Normal:
    {
      int draw_both = SceneMustDrawBoth(G);
      float alpha = (SettingGetGlobal_b(G,cSetting_opaque_background) ? 1.0F : 0.0F);
      v=SettingGetfv(G,cSetting_bg_rgb);
      if(G->HaveGUI && G->ValidContext) {
        if(draw_both) {
          OrthoDrawBuffer(G,GL_BACK_LEFT);
        } else {
          OrthoDrawBuffer(G,GL_BACK);
        }
        glClearColor(v[0],v[1],v[2],alpha);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        /* insert OpenGL context validation code here? */
        SceneRender(G,NULL,0,0,NULL,0,0,0,0);
        glClearColor(0.0,0.0,0.0,1.0);
        if(draw_both) {
          SceneCopy(G,GL_BACK_LEFT,true,false);
        } else {
          SceneCopy(G,GL_BACK,true,false);
        }
        /* insert OpenGL context validation code here? */
      }
    }
    break;
  }
  if(I->Image) {
     MovieSetImage(G,
                   MovieFrameToImage(G,SettingGetGlobal_i(G,cSetting_frame)-1),
                   I->Image);
     I->MovieOwnsImageFlag=true;
  } else {
    I->MovieOwnsImageFlag=false;
  }
  if(I->Image)
    I->CopyType=true;
  return valid;
}
/*========================================================================*/
static void SceneUpdateCameraRock(PyMOLGlobals *G,int dirty) {

  register CScene *I=G->Scene;
  float ang_cur,disp,diff;
  float sweep_angle = SettingGetGlobal_f(G,cSetting_sweep_angle);
  float sweep_speed = SettingGetGlobal_f(G,cSetting_sweep_speed);
  float sweep_phase = SettingGetGlobal_f(G,cSetting_sweep_phase);
  int sweep_mode = SettingGetGlobal_i(G,cSetting_sweep_mode);
  float shift = (float)(PI/2.0F);
  
  
  switch(sweep_mode) {
  case 0:
  case 1:
  case 2:
    if(sweep_angle<=0.0F) {
      diff = (float)((PI/180.0F)*I->RenderTime*10);
    } else {
      ang_cur = (float)(I->SweepTime*sweep_speed) + sweep_phase;
      disp = (float)(sweep_angle*(PI/180.0F)*sin(ang_cur)/2);
      diff = (float)(disp-I->LastSweep);
      I->LastSweep = disp;
    }
    switch(sweep_mode) {
    case 0:
      SceneRotateWithDirty(G,(float)(180*diff/PI),0.0F,1.0F,0.0F,dirty);
      break;
    case 1:
      SceneRotateWithDirty(G,(float)(180*diff/PI),1.0F,0.0F,0.0F,dirty);
      break;
    case 2: /* z-rotation...useless! */
      SceneRotateWithDirty(G,(float)(180*diff/PI),0.0F,0.0F,1.0F,dirty);
      break;
    }
    break;
  case 3: /* nutate */
    SceneRotateWithDirty(G,(float)(-I->LastSweepY),0.0F,1.0F,0.0F,dirty);
    SceneRotateWithDirty(G,(float)(-I->LastSweepX),1.0F,0.0F,0.0F,dirty);
    ang_cur = (float)(I->SweepTime*sweep_speed) + sweep_phase;
    
    I->LastSweepX = (float)(sweep_angle*sin(ang_cur)/2);
    I->LastSweepY = (float)(sweep_angle*sin(ang_cur+shift)/2);
    
    if(I->SweepTime*sweep_speed<PI) {
      float factor = (float)((I->SweepTime*sweep_speed)/PI);
      I->LastSweepX *= factor;
      I->LastSweepY *= factor;
    }
    SceneRotateWithDirty(G,(float)I->LastSweepX,1.0F,0.0F,0.0F,dirty);
    SceneRotateWithDirty(G,(float)I->LastSweepY,0.0F,1.0F,0.0F,dirty);
    break;
  }
}

/*========================================================================*/
void SceneIdle(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  double renderTime;
  double minTime;
  int frameFlag = false;

  if(I->PossibleSingleClick==2) {
    double now = UtilGetSeconds(G);
    double single_click_delay = I->SingleClickDelay;
    double diff = now-I->LastReleaseTime;
    if(diff>single_click_delay) {
      /* post a single click processing event */
      SceneDeferClickWhen(I->Block, 
                          I->LastButton + P_GLUT_SINGLE_LEFT,
                          I->LastWinX, I->LastWinY,
                          I->LastClickTime,
                          I->LastMod); /* push a click onto the queue */
      
      I->PossibleSingleClick = 0;
      OrthoDirty(G); /* force an update */
    }
  }
  if(!OrthoDeferredWaiting(G)) {
    if(MoviePlaying(G)) {
      renderTime = UtilGetSeconds(G) - I->LastFrameTime;
      {
        float fps = SettingGet(G,cSetting_movie_fps);
        if(fps<=0.0F) {
          if(fps<0.0)
            minTime = 0.0; /* negative fps means full speed */
          else /* 0 fps means use movie_delay instead */
            minTime = SettingGet(G,cSetting_movie_delay)/1000.0;
          if(minTime>=0)
            fps = 1.0/minTime;
          else
            fps = 1000.0F;
        } else {
          minTime = 1.0/fps;
        }
        if(renderTime >= (minTime-I->LastFrameAdjust)) {
          float adjust = (renderTime - minTime);
          if((fabs(adjust)<minTime) && (fabs(I->LastFrameAdjust)<minTime)) {
            float new_adjust = (renderTime - minTime) + I->LastFrameAdjust;
            I->LastFrameAdjust = (new_adjust + fps*I->LastFrameAdjust)/(1+fps);
          } else {
            I->LastFrameAdjust = 0.0F;
          }
          frameFlag=true;
        }
      }
    } else if(ControlRocking(G)) {
      renderTime = -I->LastSweepTime + UtilGetSeconds(G);
      minTime=SettingGet(G,cSetting_rock_delay)/1000.0;
      if(renderTime>=minTime) {
        I->LastSweepTime=UtilGetSeconds(G);
        I->SweepTime+=I->RenderTime;
        SceneUpdateCameraRock(G,true);
      }
    }
    
    if(MoviePlaying(G) && frameFlag) {
      I->LastFrameTime = UtilGetSeconds(G);
      if((SettingGetGlobal_i(G,cSetting_frame)-1)==(I->NFrame-1)) {
        if((int)SettingGet(G,cSetting_movie_loop)) {
          SceneSetFrame(G,7,0);
        } else
          MoviePlay(G,cMovieStop);
      } else 
        SceneSetFrame(G,5,1);
    }
  }
}
/*========================================================================*/
void SceneWindowSphere(PyMOLGlobals *G,float *location,float radius)
{
  register CScene *I=G->Scene;
  float v0[3];
  float dist;
  float aspRat;
  float fov;

  if(I->Height && I->Width) {
    aspRat = ((float) I->Width) / ((float) I->Height);
  } else {
    aspRat = 1.3333F;
  }

  /* find where this point is in relationship to the origin */
  subtract3f(I->Origin,location,v0); 

  dist = I->Pos[2];

  MatrixTransformC44fAs33f3f(I->RotMatrix,v0,I->Pos); /* convert to view-space */
  fov = SettingGet(G,cSetting_field_of_view);
  if(aspRat<1.0)
    fov *= aspRat;

  dist = (float)(radius/tan((fov/2.0)*cPI/180.0));

  I->Pos[2]-=dist;
  I->Front=(-I->Pos[2]-radius*1.2F);
  I->Back=(-I->Pos[2]+radius*1.2F);
  I->FrontSafe= GetFrontSafe(I->Front,I->Back);
  I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
  SceneRovingDirty(G);
}
/*========================================================================*/
void SceneRelocate(PyMOLGlobals *G,float *location)
{
  register CScene *I=G->Scene;
  float v0[3];
  float slab_width;
  float dist;

  slab_width = I->Back-I->Front;

  /* find out how far camera was from previous origin */
  dist = I->Pos[2];

  /* find where this point is in relationship to the origin */
  subtract3f(I->Origin,location,v0); 

  /*  printf("%8.3f %8.3f %8.3f\n",I->Front,I->Pos[2],I->Back);*/

  MatrixTransformC44fAs33f3f(I->RotMatrix,v0,I->Pos); /* convert to view-space */

  I->Pos[2]=dist;
  I->Front=(-I->Pos[2]-(slab_width*0.50F));
  I->Back=(-I->Pos[2]+(slab_width*0.50F));
  I->FrontSafe= GetFrontSafe(I->Front,I->Back);
  I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
  SceneRovingDirty(G);

}
/*========================================================================*/
void SceneOriginGet(PyMOLGlobals *G,float *origin)
{
  register CScene *I=G->Scene;
  copy3f(I->Origin,origin);
}
/*========================================================================*/
void SceneOriginSet(PyMOLGlobals *G,float *origin,int preserve)
{
  register CScene *I=G->Scene;
  float v0[3],v1[3];
  
  if(preserve) {/* preserve current viewing location */
    subtract3f(origin,I->Origin,v0); /* model-space translation */
    MatrixTransformC44fAs33f3f(I->RotMatrix,v0,v1); /* convert to view-space */
    add3f(I->Pos,v1,I->Pos); /* offset view to compensate */
  }
  I->Origin[0]=origin[0]; /* move origin */
  I->Origin[1]=origin[1];
  I->Origin[2]=origin[2];
  SceneInvalidate(G);
}
/*========================================================================*/
int SceneObjectAdd(PyMOLGlobals *G,CObject *obj)
{
  register CScene *I=G->Scene;
  ObjRec *rec = NULL;
  ListElemAlloc(G,rec,ObjRec);
  rec->next=NULL;
  obj->Enabled=true;
  rec->obj=obj;
  ListAppend(I->Obj,rec,next,ObjRec);
  SceneCountFrames(G);
  SceneChanged(G);
  return 1;
}
/*========================================================================*/
int SceneObjectIsActive(PyMOLGlobals *G,CObject *obj)
{
  int result=false;
  register CScene *I=G->Scene;
  ObjRec *rec = NULL;
  while(ListIterate(I->Obj,rec,next))
    if(rec->obj==obj) {
      result=true;
      break;
    }
  return result;
}
int SceneObjectDel(PyMOLGlobals *G,CObject *obj)
{
  register CScene *I=G->Scene;
  ObjRec *rec = NULL;
  int defer_builds_mode = SettingGetGlobal_b(G,cSetting_defer_builds_mode);

  if(!obj) { /* deletes all members */
    while(ListIterate(I->Obj,rec,next)) {
      if(rec) {
        if(defer_builds_mode >= 3) { 
          /* purge graphics representation when no longer used */
          if(rec->obj->fInvalidate) 
            rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvPurge,-1);
        }
        ListDetach(I->Obj,rec,next,ObjRec);
        ListElemFree(rec);
      }
    }
  } else {
    while(ListIterate(I->Obj,rec,next))
      if(rec->obj==obj)
        break;
    if(rec) {
      if(defer_builds_mode >= 3) { 
        /* purge graphics representation when no longer used */
        if(rec->obj->fInvalidate) 
          rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvPurge,-1);
      }
      rec->obj->Enabled=false;
      ListDetach(I->Obj,rec,next,ObjRec);
      ListElemFree(rec);
    }
  }
  SceneCountFrames(G);
  SceneInvalidate(G);
  return 0;
}
/*========================================================================*/
int SceneLoadPNG(PyMOLGlobals *G,char *fname,int movie_flag,int stereo,int quiet) 
{
  register CScene *I=G->Scene;
  int ok=false;
  if(I->Image) {
     if(I->MovieOwnsImageFlag) {
        I->MovieOwnsImageFlag=false;
        I->Image=NULL;
     } else {
       ScenePurgeImage(G);
     }
    I->CopyType=false;
  }
  I->Image=Calloc(ImageType,1);
  if(MyPNGRead(fname,
               (unsigned char**)&I->Image->data,
               (unsigned int*)&I->Image->width,
               (unsigned int*)&I->Image->height)) {
    I->Image->size = I->Image->width * I->Image->height * 4;
    if(!quiet) {
      PRINTFB(G,FB_Scene,FB_Details)
        " Scene: loaded image from '%s'.\n",fname
        ENDFB(G);
    }
    if((stereo>0) || ((stereo<0) &&
                    (I->Image->width == 2*I->Width)&&
                    (I->Image->height == I->Height))) {
      unsigned char *tmp = Alloc(unsigned char,I->Image->size);
      if(tmp) {
        I->Image->width /=2;
        I->Image->stereo = true;
        I->Image->size /=2;
        deinterlace((unsigned int*)tmp,
                    (unsigned int*)I->Image->data,
                    I->Image->width,I->Image->height,
                    (stereo==2));
        FreeP(I->Image->data);
        I->Image->data = tmp;
      }
    }
        
    I->CopyType=true;
    I->CopyForced=true;
    OrthoRemoveSplash(G);
    SettingSet(G,cSetting_text,0.0);
    if(movie_flag&&
       I->Image&&I->Image->data&&
       (I->Image->height==I->Height)&&
       (I->Image->width==I->Width)) {
      MovieSetImage(G,
                    MovieFrameToImage(G,
                                      SettingGetGlobal_i(G,cSetting_frame)-1)
                    ,I->Image);
      I->MovieOwnsImageFlag=true;
      I->MovieFrameFlag=true;
    } else {
      I->MovieOwnsImageFlag=false;
      I->DirtyFlag=false; /* make sure we don't overwrite image */
    }
    OrthoDirty(G);
    ok=true;
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Scene,FB_Errors)
        " Scene: unable to load image from '%s'.\n",fname
        ENDFB(G);
    }
  }
  return(ok);
}
/*static unsigned int byte_max(unsigned int value)
{
  return (value>0xFF) ? 0xFF : value;
}
*/

#define SceneClickMargin 2
#define SceneTopMargin 0
#define SceneToggleMargin 2
#define SceneRightMargin 0
#define SceneToggleWidth 17
#define SceneToggleSize 16
#define SceneToggleTextShift 4
#define SceneTextLeftMargin 1
#define SceneScrollBarMargin 1
#define SceneScrollBarWidth 13

static void draw_button(int x2,int y2, int z, int w, int h, float *light, float *dark, float *inside)
{
  glColor3fv(light);
  glBegin(GL_POLYGON);
  glVertex3i(x2,y2,z);
  glVertex3i(x2,y2+h,z);
  glVertex3i(x2+w,y2+h,z);
  glVertex3i(x2+w,y2,z);
  glEnd();
  
  glColor3fv(dark);
  glBegin(GL_POLYGON);
  glVertex3i(x2+1,y2,z);
  glVertex3i(x2+1,y2+h-1,z);
  glVertex3i(x2+w,y2+h-1,z);
  glVertex3i(x2+w,y2,z);
  glEnd();
  
  if(inside) {
    glColor3fv(inside);
    glBegin(GL_POLYGON);
    glVertex3i(x2+1,y2+1,z);
    glVertex3i(x2+1,y2+h-1,z);
    glVertex3i(x2+w-1,y2+h-1,z);
    glVertex3i(x2+w-1,y2+1,z);
    glEnd();
  } else { /* rainbow */
    glBegin(GL_POLYGON);
    glColor3f(1.0F,0.1F,0.1F);
    glVertex3i(x2+1,y2+1,z);
    glColor3f(0.1F,1.0F,0.1F);
    glVertex3i(x2+1,y2+h-1,z);
    glColor3f(1.0F,1.0F,0.1F);
    glVertex3i(x2+w-1,y2+h-1,z);
    glColor3f(0.1F,0.1F,1.0F);
    glVertex3i(x2+w-1,y2+1,z);
    glEnd();
  }

}
int SceneSetNames(PyMOLGlobals *G,PyObject *list)
{
#ifndef _PYMOL_NOPY
  register CScene *I = G->Scene;
  int ok = PConvPyListToStrVLAList(list, &I->SceneNameVLA, &I->NScene);
  if(ok) {
    VLACheck(I->SceneVLA, SceneElem, I->NScene);
    {
      int a;
      char *c = I->SceneNameVLA;
      SceneElem *elem = I->SceneVLA;
      for(a=0;a<I->NScene;a++) {
        elem->name = c;
        elem->len = strlen(c);
        elem->drawn = false;
        c += elem->len+1;
        elem++;
      }
    }
  }
  OrthoDirty(G);
  return ok;
#else
  return 0;
#endif
}

/*========================================================================*/
static void SceneDrawButtons(Block *block, int draw_for_real)
{
#ifndef _PYMOL_NOPY  
  PyMOLGlobals *G=block->G;
  register CScene *I = G->Scene;
  int x,y,xx,x2;
  char *c=NULL;
  float enabledColor[3] = { 0.5F, 0.5F, 0.5F };
  float pressedColor[3] = { 0.7F, 0.7F, 0.7F };
  float disabledColor[3] = { 0.25F, 0.25F, 0.25F };
  float lightEdge[3] = {0.6F, 0.6F, 0.6F };
  float darkEdge[3] = {0.35F, 0.35F, 0.35F };
  int charWidth = 8;
  int n_ent;
  int n_disp;
  int skip=0;
  int row = -1;
  int lineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
  int text_lift = (lineHeight/2)-5;
  int op_cnt = 1;
  
  if( ((G->HaveGUI && G->ValidContext)|| (!draw_for_real)) && 
     ((block->rect.right-block->rect.left)>6) && (I->NScene)) {
    int max_char;
    int nChar;
    I->ButtonsShown = true;

    /* do we have enough structures to warrant a scroll bar? */
    n_ent = I->NScene;
    
    n_disp = (((I->Block->rect.top-I->Block->rect.bottom)-(SceneTopMargin))/lineHeight)-1;
    if(n_disp<1) n_disp=1;
    
    {
      int i;
      for(i=0;i<I->NScene;i++)
        I->SceneVLA[i].drawn = false;
    }
    if(n_ent>n_disp) {
      int bar_maxed = ScrollBarIsMaxed(I->ScrollBar);
      if(!I->ScrollBarActive) {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        if(bar_maxed) {
          ScrollBarMaxOut(I->ScrollBar);
          I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
        } else {
          ScrollBarSetValue(I->ScrollBar,0);
          I->NSkip =0;
        }
      } else {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        if(bar_maxed)
          ScrollBarMaxOut(I->ScrollBar);
        I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
      }
      I->ScrollBarActive = 1;

    } else {
      I->ScrollBarActive = 0;
      I->NSkip =0;
    }

    max_char = (((I->Block->rect.right-I->Block->rect.left)-
		 (SceneTextLeftMargin+SceneRightMargin+4)) -
		(op_cnt*SceneToggleWidth));
    if(I->ScrollBarActive) {
      max_char -= (SceneScrollBarMargin+SceneScrollBarWidth);
    }      
    max_char/=charWidth;

    if(I->ScrollBarActive) {
      ScrollBarSetBox(I->ScrollBar,I->Block->rect.top-SceneScrollBarMargin,
                      I->Block->rect.left+SceneScrollBarMargin,
                      I->Block->rect.bottom+2,
                      I->Block->rect.left+SceneScrollBarMargin+SceneScrollBarWidth);
      if(draw_for_real) 
	ScrollBarDoDraw(I->ScrollBar);
    }
    
    skip=I->NSkip;
    x = I->Block->rect.left+SceneTextLeftMargin;
    
    /*    y = ((I->Block->rect.top-lineHeight)-SceneTopMargin)-lineHeight;*/
    
    {
      int n_vis = n_disp;
      if(n_ent<n_vis)
        n_vis = n_ent;
      y = (I->Block->rect.bottom+SceneBottomMargin)+(n_vis-1)*lineHeight;
    }

    /*    xx = I->Block->rect.right-SceneRightMargin-SceneToggleWidth*(cRepCnt+op_cnt);*/
    xx = I->Block->rect.right-SceneRightMargin-SceneToggleWidth*(op_cnt);

    if(I->ScrollBarActive) {
      x+=SceneScrollBarWidth+SceneScrollBarMargin;
    }
    {
      int i;

      for(i=0;i<n_ent;i++) {
        if(skip) {
          skip--;
        } else {
          row++;
          x2=xx;
          nChar = max_char;
          
          if((x-SceneToggleMargin)-(xx-SceneToggleMargin)>-10) {
            x2 = x+10;
          }
          {
            float toggleColor[3] = { 0.5F, 0.5F, 1.0F };
           
	    if(draw_for_real) {
	      glColor3fv(toggleColor);
	      
	      TextSetColor(G,I->Block->TextColor);
	      TextSetPos2i(G,x+2,y+text_lift);
	    }
            {
              int len;
              char *cur_name = SettingGetGlobal_s(G,cSetting_scene_current_name);
              SceneElem *elem = I->SceneVLA + i;
              int item = I->NSkip + row;
              c = elem->name;
              len = elem->len;

              x2 = xx;
              if(len>max_char)
                len = max_char;
              x2 = x + len * charWidth + 6;

              /* store rectangles for finding clicks */

              elem->drawn = true;

              elem->x1 = x;
              elem->y1 = y;
              elem->x2 = x2;
              elem->y2 = y+lineHeight;

	      if(I->ButtonMargin<x2)
		I->ButtonMargin = x2;

	      if(draw_for_real) {
	      
		if((item==I->Pressed)&&(item==I->Over)) {
		  draw_button(x,y,0,(x2-x)-1,(lineHeight-1),lightEdge,darkEdge,pressedColor);
		} else if(cur_name&&elem->name&&(!strcmp(elem->name,cur_name))) {
		  draw_button(x,y,0,(x2-x)-1,(lineHeight-1),lightEdge,darkEdge,enabledColor);
		} else {
		  draw_button(x,y,0,(x2-x)-1,(lineHeight-1),lightEdge,darkEdge,disabledColor);
		}
		
		TextSetColor(G,I->Block->TextColor);
		
		if(c) {
		  while(*c) {
		    if((nChar--)>0)
		      TextDrawChar(G,*(c++));
		    else
		      break;
		  }
		}
	      }
            }
          }
          y-=lineHeight;
          if(y<(I->Block->rect.bottom))
            break;
        }
      }
    }
    I->HowFarDown = y;
    I->ButtonsValid = true;
  }
#endif
}

static void SceneUpdateButtons(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  SceneDrawButtons(I->Block,false);
}

void SceneDraw(Block *block)
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  int overlay,text;
  
  if(G->HaveGUI && G->ValidContext) {

    I->ButtonsShown = false;

    overlay = OrthoGetOverlayStatus(G);
    text = (int)SettingGet(G,cSetting_text);

    if(((!text)||overlay) &&
       (I->CopyType == true) && 
       I->Image && I->Image->data) {
      
      int show_alpha = SettingGetGlobal_b(G,cSetting_show_alpha_checker);
      float *bg_color=SettingGetfv(G,cSetting_bg_rgb);              
      unsigned int bg_rr,bg_r = (unsigned int)(255*bg_color[0]);
      unsigned int bg_gg,bg_g = (unsigned int)(255*bg_color[1]);
      unsigned int bg_bb,bg_b = (unsigned int)(255*bg_color[2]);
      
      unsigned char *data = I->Image->data;

      int width = I->Image->width;
      int height = I->Image->height;

      
      if(I->Image->stereo) {
        int buffer;
        glGetIntegerv(GL_DRAW_BUFFER,(GLint*)&buffer);
        if(buffer==GL_BACK_RIGHT) /* hardware stereo */
          data += I->Image->size; 
        else {
          switch(OrthoGetRenderMode(G)) {
          case cStereo_geowall:
            data += I->Image->size; 
            break;
          }
        }
        /* if drawing the right buffer, then draw the right image */
      }
      
      if((height>I->Height)||
         (width>I->Width)) { /* image is oversize */
        {
          int factor = 1;
          int shift = 0;
          register int tmp_height = I->Image->height;
          register int tmp_width = I->Image->width;
          int src_row_bytes = I->Image->width * 4;
          unsigned int color_word;
          float rgba[4] = { 0.0F, 0.0F, 0.0F, 1.0F };
          
          ColorGetBkrdContColor(G,rgba,false);
          color_word = ColorGet32BitWord(G,rgba);
          
          while(tmp_height&&tmp_width&&
                ((tmp_height>(I->Height-3))||(tmp_width>(I->Width-3)))) {
            tmp_height = (tmp_height>>1);
            tmp_width = (tmp_width>>1);
            factor = (factor<<1);
            shift++;
          }
          tmp_width+=2;
          tmp_height+=2;
          
          if(tmp_height&&tmp_width) {
            unsigned int buffer_size = tmp_height * tmp_width * 4;
            unsigned char *buffer = Alloc(unsigned char, buffer_size);
            
            
            if(buffer && data) {
              unsigned char *p = data;
              unsigned char *q = buffer;
              register unsigned char *pp, *ppp, *pppp;
              register int a,b,c,d;
              register unsigned int c1,c2,c3,c4,alpha,tot,bg;
              unsigned int factor_col_bytes = factor * 4;
              unsigned int factor_row_bytes = factor * src_row_bytes;
              
              shift = shift + shift;
              for(b=0;b<tmp_height;b++) { /* rows */
                pp = p;
                if((!b)||(b==(tmp_height-1))) {
                  for(a=0;a<tmp_width;a++) { /* border */
                    *((unsigned int*)(q)) = color_word;
                    q+=4;
                  }
                } else {
                  for(a=0;a<tmp_width;a++) { /* cols */
                    ppp = pp;
                    if((!a)||(a==(tmp_width-1))) { /* border */
                      *((unsigned int*)(q)) = color_word;
                      q+=4;
                    } else {
                      c1 = c2 = c3 = c4 = tot = 0;

                      if(show_alpha&&(((a>>4)+(b>>4))&0x1)) { /* introduce checkerboard */
                        bg_rr = ((bg_r&0x80) ? bg_r - TRN_BKG : bg_r + TRN_BKG);
                        bg_gg = ((bg_g&0x80) ? bg_g - TRN_BKG : bg_g + TRN_BKG);
                        bg_bb = ((bg_b&0x80) ? bg_b - TRN_BKG : bg_b + TRN_BKG);
                      } else {
                        bg_rr = bg_r;
                        bg_gg = bg_g;
                        bg_bb = bg_b;
                      }

                      for(d=0;d<factor;d++) { /* box rows */
                        pppp = ppp;
                        for(c=0;c<factor;c++) { /* box cols */
                          alpha = pppp[3];
                          c1 += *(pppp++) * alpha;
                          c2 += *(pppp++) * alpha;
                          c3 += *(pppp++) * alpha;
                          pppp++;
                          c4 += alpha;
                          tot += 0xFF;
                        }
                        ppp += src_row_bytes;
                      }
                      if(c4) {
                        bg = tot-c4;
                        *(q++)= (c1 + bg_rr*bg)/tot;
                        *(q++)= (c2 + bg_gg*bg)/tot;
                        *(q++)= (c3 + bg_bb*bg)/tot;
                        *(q++)= 0xFF;
                      } else {
                        *(q++)= bg_rr;
                        *(q++)= bg_gg;
                        *(q++)= bg_bb;
                        *(q++)= 0xFF;
                      }
                      pp += factor_col_bytes;
                    }
                  }
                  p+= factor_row_bytes;
                }
              }

              glRasterPos3i((int)((I->Width-tmp_width)/2+I->Block->rect.left),
                            (int)((I->Height-tmp_height)/2+I->Block->rect.bottom),-10);
              PyMOLDrawPixels(tmp_width,tmp_height,GL_RGBA,GL_UNSIGNED_BYTE,buffer);
            }
            FreeP(buffer);
          }
          {
            char buffer[255];
            int text_pos = (I->Height - tmp_height) / 2 - 15;
            int x_pos, y_pos;
            if(text_pos<0) {
              text_pos = (I->Height - tmp_height) / 2 + 3;
              x_pos = (I->Width - tmp_width) / 2 + 3;
              y_pos = text_pos;
            } else {
              x_pos = (I->Width - tmp_width) / 2;
              y_pos = text_pos;
            }
                
            sprintf(buffer,"Image size = %d x %d",
                    I->Image->width,
                    I->Image->height);
                
            TextSetColor3f(G,rgba[0],rgba[1],rgba[2]);
            TextDrawStrAt(G,buffer,
                          x_pos + I->Block->rect.left,
                          y_pos + I->Block->rect.bottom);
          }
        }
      } else if(((width<I->Width)||(height<I->Height)) &&
                ((I->Width-width)>2) && ((I->Height-height)>2)) { /* but a border around image */
        
        unsigned int color_word;
        float rgba[4] = { 0.0F, 0.0F, 0.0F, 1.0F };
        register unsigned int tmp_height = height+2;
        register unsigned int tmp_width = width+2;
        unsigned int n_word = tmp_height * tmp_width;
        unsigned int *tmp_buffer = Alloc(unsigned int,n_word);
        ColorGetBkrdContColor(G,rgba,false);
        color_word = ColorGet32BitWord(G,rgba);
        
        if(tmp_buffer) {
          register unsigned int a,b;
          unsigned int *p=(unsigned int*)data;
          unsigned int *q=tmp_buffer;
            for(a=0;a<tmp_height;a++) {
              if((!a)||(a==(tmp_height-1))) {
                for(b=0;b<tmp_width;b++) 
                  *(q++) = color_word;
              } else {
                for(b=0;b<tmp_width;b++) {
                  if((!b)||(b==(tmp_width-1))) {        
                    *(q++) = color_word;
                  } else {
                    unsigned char *qq = (unsigned char*)q;
                    unsigned char *pp = (unsigned char*)p;
                    unsigned char bg;
                    if(show_alpha&&(((a>>4)+(b>>4))&0x1)) { /* introduce checkerboard */
                      bg_rr = ((bg_r&0x80) ? bg_r - TRN_BKG : bg_r + TRN_BKG);
                      bg_gg = ((bg_g&0x80) ? bg_g - TRN_BKG : bg_g + TRN_BKG);
                      bg_bb = ((bg_b&0x80) ? bg_b - TRN_BKG : bg_b + TRN_BKG);
                    } else {
                      bg_rr = bg_r;
                      bg_gg = bg_g;
                      bg_bb = bg_b;
                    }
                    if(pp[3]) {
                      bg = 0xFF-pp[3];
                      *(qq++)= (pp[0]*pp[3] + bg_rr*bg)/0xFF;
                      *(qq++)= (pp[1]*pp[3] + bg_gg*bg)/0xFF;
                      *(qq++)= (pp[2]*pp[3] + bg_bb*bg)/0xFF;
                      *(qq++)= 0xFF;
                    } else {
                      *(qq++)= bg_rr;
                      *(qq++)= bg_gg;
                      *(qq++)= bg_bb;
                      *(qq++)= 0xFF;
                    }
                    q++;
                    p++;
                  }
                }
              }
            }
            glRasterPos3i((int)((I->Width-tmp_width)/2+I->Block->rect.left),
                          (int)((I->Height-tmp_height)/2+I->Block->rect.bottom),-10);
            PyMOLDrawPixels(tmp_width,tmp_height,GL_RGBA,GL_UNSIGNED_BYTE,tmp_buffer);
            
        }
        FreeP(tmp_buffer);
      } else if(I->CopyForced) { /* near-exact fit */
          unsigned int color_word;
          float rgba[4] = { 0.0F, 0.0F, 0.0F, 1.0F };
          unsigned int n_word = height * width;
          unsigned int *tmp_buffer = Alloc(unsigned int,n_word);
          ColorGetBkrdContColor(G,rgba,false);
          color_word = ColorGet32BitWord(G,rgba);
          
          if(tmp_buffer) {
            register unsigned int a,b;
            unsigned int *p=(unsigned int*)data;
            unsigned int *q=tmp_buffer;
            for(a=0;a<(unsigned int)height;a++) {
              for(b=0;b<(unsigned int)width;b++) {
                unsigned char *qq = (unsigned char*)q;
                unsigned char *pp = (unsigned char*)p;
                unsigned char bg;
                if(show_alpha&&(((a>>4)+(b>>4))&0x1)) { /* introduce checkerboard */
                    bg_rr = ((bg_r&0x80) ? bg_r - TRN_BKG : bg_r + TRN_BKG);
                    bg_gg = ((bg_g&0x80) ? bg_g - TRN_BKG : bg_g + TRN_BKG);
                    bg_bb = ((bg_b&0x80) ? bg_b - TRN_BKG : bg_b + TRN_BKG);
                  } else {
                    bg_rr = bg_r;
                    bg_gg = bg_g;
                    bg_bb = bg_b;
                  }
                  if(pp[3]) {
                    bg = 0xFF-pp[3];
                    *(qq++)= (pp[0]*pp[3] + bg_rr*bg)/0xFF;
                    *(qq++)= (pp[1]*pp[3] + bg_gg*bg)/0xFF;
                    *(qq++)= (pp[2]*pp[3] + bg_bb*bg)/0xFF;
                    *(qq++)= 0xFF;
                  } else {
                    *(qq++)= bg_rr;
                    *(qq++)= bg_gg;
                    *(qq++)= bg_bb;
                    *(qq++)= 0xFF;
                  }
                  q++;
                  p++;
              }
            }
          }
          glRasterPos3i((int)((I->Width-width)/2+I->Block->rect.left),
                        (int)((I->Height-height)/2+I->Block->rect.bottom),-10);
          PyMOLDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,tmp_buffer);
          FreeP(tmp_buffer);
      } else { /* not a forced copy, so don't show/blend alpha */
        glRasterPos3i((int)((I->Width-width)/2+I->Block->rect.left),
                      (int)((I->Height-height)/2+I->Block->rect.bottom),-10);
        PyMOLDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,data); 
        
      }
      I->RenderTime = -I->LastRender;
      I->LastRender = UtilGetSeconds(G);
      I->RenderTime += I->LastRender;

    }
    if(SettingGetGlobal_b(G,cSetting_scene_buttons)&&
       (SettingGetGlobal_i(G,cSetting_scene_buttons_mode)==1)) {
      SceneDrawButtons(block,true);
    } else {
      I->ButtonMargin = 0;
    }
  }
}
int SceneGetButtonMargin(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  return I->ButtonMargin;
}
/*========================================================================*/

typedef unsigned char pix[4];
#define cRange 7
/*typedef pix pix_array[cRange*2+1][cRange*2+1];*/

unsigned int SceneFindTriplet(PyMOLGlobals *G,int x,int y,GLenum gl_buffer) 
{
  int result = 0;
  /*int before_check[100];
  int *int_ptr;
*/
  pix *buffer = NULL;
  pix *extra_safe_buffer = NULL;

 /*int after_check[100];*/
  /* pix_array *array_ptr;
  char *safe_place;
*/
  int a,b,d,flag;
  int h = (cRange*2+1), w = (cRange*2+1);

  int debug = false;
  unsigned char *c;
  int strict = false;
  GLint rb,gb,bb;
  int bkrd_alpha = 0xFF;
  int check_alpha = false;


  if(G->HaveGUI && G->ValidContext) { /*just in case*/
  
    glGetIntegerv(GL_RED_BITS,&rb);
    glGetIntegerv(GL_GREEN_BITS,&gb);
    glGetIntegerv(GL_BLUE_BITS,&bb);

    if((rb>=8)&&(gb>=8)&&(bb>=8))
        strict = true;

    if(Feedback(G,FB_Scene,FB_Debugging)) debug=true;
    
    glReadBuffer(gl_buffer);

    extra_safe_buffer=Alloc(pix,w*h*21); 
    buffer=extra_safe_buffer+(w*h*10);

    PyMOLReadPixels(x-cRange,y-cRange,cRange*2+1,cRange*2+1,GL_RGBA,GL_UNSIGNED_BYTE,&buffer[0][0]);

      if(debug) {
      for(a=0;a<=(cRange*2);a++)
        {
          for(b=0;b<=(cRange*2);b++)
            printf("%2x ",(buffer[a+b*w][0]+buffer[a+b*w][1]+buffer[a+b*w][2])&0xFF);
          printf("\n");
        }
      printf("\n");     
      for(a=0;a<=(cRange*2);a++)
        {
          for(b=0;b<=(cRange*2);b++)
            printf("%02x ",(buffer[a+b*w][3])&0xFF);
          printf("\n");
        }
      printf("\n");     
       for(a=0;a<=(cRange*2);a++)
        {
          for(b=0;b<=(cRange*2);b++)
            printf("%02x%02x%02x ",(buffer[a+b*w][0])&0xFF,(buffer[a+b*w][1])&0xFF,(buffer[a+b*w][2])&0xFF);
          printf("\n");
        }
       printf("\n");     
     }

     /* first, check to make sure bkrd_alpha is correct 
        (this is a bug for systems with broken alpha, such as Extreme 3D on Solaris 8 */

     flag=true;
     for(d=0;flag&&(d<cRange);d++)
       for(a=-d;flag&&(a<=d);a++)
         for(b=-d;flag&&(b<=d);b++) {
           c = &buffer[(a+cRange)+(b+cRange)*w][0];
           if(c[3]==bkrd_alpha) {
             check_alpha = true;
             flag=false;
           }
         }

     /* now find the correct pixel */

     flag=true;
     for(d=0;flag&&(d<cRange);d++)
       for(a=-d;flag&&(a<=d);a++)
         for(b=-d;flag&&(b<=d);b++) {
           c = &buffer[(a+cRange)+(b+cRange)*w][0];
           if(((c[3]==bkrd_alpha)||(!check_alpha))&&
              ((c[1]&0x8)&&
               ((!strict)||
                (((c[1]&0xF)==8)&&
                 ((c[0]&0xF)==0)&&
                 ((c[2]&0xF)==0)
                 )))) { /* only consider intact, saturated pixels */
             flag = false;
             result =  ((c[0]>>4)&0xF)+(c[1]&0xF0)+((c[2]<<4)&0xF00);
             if(debug) {
               printf("%2x %2x %2x %d\n",c[0],c[1],c[2],result);
             }
           }
         }
     FreeP(extra_safe_buffer);
  }
  return(result);
}
/*========================================================================*/
unsigned int *SceneReadTriplets(PyMOLGlobals *G,int x,int y,int w,int h,GLenum gl_buffer)
{ 
  unsigned int *result = NULL;
  pix *buffer=NULL;
  pix *extra_safe_buffer=NULL;
  int a,b;
  unsigned char *c;
  int cc = 0;
  int dim[3];
  int strict = false;
  int bkrd_alpha = 0xFF;
  int check_alpha = false;
 
  GLint rb,gb,bb;

  dim[0]=w;
  dim[1]=h;
  
  if(w<1) w=1;
  if(h<1) h=1;
  if(G->HaveGUI && G->ValidContext) { /*just in case*/
    
    glGetIntegerv(GL_RED_BITS,&rb);
    glGetIntegerv(GL_RED_BITS,&gb);
    glGetIntegerv(GL_RED_BITS,&bb);
    
    if((rb>=8)&&(gb>=8)&&(bb>=8))
        strict = true;
    
    /* create some safe RAM on either side of the read buffer -- buggy
       ReadPixels implementations tend to trash RAM surrounding the
       target block */

    extra_safe_buffer=Alloc(pix,w*h*11); 
    buffer=extra_safe_buffer+(w*h*5);
    
    result = VLAlloc(unsigned int,w*h);
    glReadBuffer(gl_buffer);
    PyMOLReadPixels(x,y,w,h,GL_RGBA,GL_UNSIGNED_BYTE,&buffer[0][0]);
    
     /* first, check to make sure bkrd_alpha is correct 
        (this is a bug for systems with broken alpha, such as Extreme 3D on Solaris 8 */

    for(a=0;a<w;a++)
      for(b=0;b<h;b++)
        {
          c = &buffer[a+b*w][0];
          if(c[3]==bkrd_alpha) {
            check_alpha = true;
          }
        }
    
    /* now read pixels */

    for(a=0;a<w;a++)
      for(b=0;b<h;b++)
        {
          c = &buffer[a+b*w][0];
          if((((c[3]==bkrd_alpha)||(!check_alpha)))&&
             ((c[1]&0x8)&&
              ((!strict)||
               (((c[1]&0xF)==8)&&
                ((c[0]&0xF)==0)&&
                ((c[2]&0xF)==0)
                )))) { /* only consider intact, saturated pixels */
            VLACheck(result,unsigned int,cc+1);
            result[cc] =  ((c[0]>>4)&0xF)+(c[1]&0xF0)+((c[2]<<4)&0xF00);
            result[cc+1] = b+a*h;
            /*printf("%2x %2x %2x %d\n",c[0],c[1],c[2],result[cc]);*/
            cc+=2;
          }
        }
    FreeP(extra_safe_buffer);
    VLASize(result,unsigned int,cc);
  }
  return(result);

}
/*========================================================================*/
static int SceneRelease(Block *block,int button,int x,int y,int mod, double when) 
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  int release_handled = false;
  if(I->ButtonsShown && I->PressMode) {
    if(I->ScrollBarActive) {
      if((x-I->Block->rect.left)<(SceneScrollBarWidth+SceneScrollBarMargin)) {
        ScrollBarDoRelease(I->ScrollBar,button,x,y,mod);
        release_handled = true;
      }
    }
    if(!release_handled) {
      int ungrab=true;
      if(I->PressMode) {
        int i;
        SceneElem *elem = I->SceneVLA;
        I->Over = -1;
        for(i=0;i<I->NScene;i++) {
          if(elem->drawn&&
             (x>=elem->x1) &&
             (y>=elem->y1) &&
             (x<elem->x2) &&
             (y<elem->y2)) {
            I->Over = i;
            break;
          }
          elem++;
        }
        if(I->Over>=0) {
          release_handled=true;
          switch(I->PressMode) {
          case 1:
            if(I->Over == I->Pressed) {
              OrthoLineType buffer;
              sprintf(buffer,"cmd.scene('''%s''')",elem->name);
              PParse(G,buffer);
              PFlush(G);
              PLog(G,buffer,cPLog_pym);
            }
            break;
          case 2:
            { 
              char *cur_name = SettingGetGlobal_s(G,cSetting_scene_current_name);
              if(cur_name && elem->name && (strcmp(cur_name, elem->name))) {
                OrthoLineType buffer;
                sprintf(buffer,"cmd.scene('''%s''')",elem->name);
                PParse(G,buffer);
                PFlush(G);
                PLog(G,buffer,cPLog_pym);
              }
            }
            break;
          case 3:
            if(I->Pressed==I->Over) {
              MenuActivate1Arg(G,I->LastWinX,I->LastWinY+20, /* scene menu */
                                 I->LastWinX,I->LastWinY,
                                 true,"scene_menu",elem->name);
              ungrab=false;
            }
            break;
          }
        }
      }
      I->LastPickVertexFlag=false;
      I->Pressed = -1;
      I->Over = -1;
      I->PressMode = 0;
      if(ungrab)
        OrthoUngrab(G);
    }
  }
  if(!release_handled) {
    ObjectMolecule *obj;
    I->LastReleaseTime = when;
    if(I->PossibleSingleClick==1) {
      double slowest_single_click = 0.25F;
      double diff = when-I->LastClickTime;
      
      slowest_single_click += I->ApproxRenderTime;
      
      if((diff<0.0)||(diff>slowest_single_click))
        I->PossibleSingleClick = 0;
      else {
        int but = -1;
        I->PossibleSingleClick = 2;
        I->SingleClickDelay = 0.15;
        
        switch(I->LastButton) {
        case P_GLUT_LEFT_BUTTON:
          but = P_GLUT_DOUBLE_LEFT;
          break;
        case P_GLUT_MIDDLE_BUTTON:
          but = P_GLUT_DOUBLE_MIDDLE;
          break;
        case P_GLUT_RIGHT_BUTTON:
          but = P_GLUT_DOUBLE_RIGHT;
          break;
        }
        if(but>0) {
          int mode=ButModeTranslate(G,but,mod);
          if(mode == cButModeNone)
            I->SingleClickDelay = 0.0; /* no double-click set? force immediate single click */
        }
      }
    }
    if(I->LoopFlag && (I->PossibleSingleClick!=2))
      return SceneLoopRelease(block,button,x,y,mod);
    OrthoUngrab(G);
    I->LoopFlag=false;
    if(I->SculptingFlag) {
      /* SettingSet(G,cSetting_sculpting,1); */
      obj=(ObjectMolecule*)I->LastPicked.context.object;
      if(obj) {
        obj->AtomInfo[I->LastPicked.src.index].protekted=I->SculptingSave;
      }
      I->SculptingFlag=0;
    }
  }
  return 1;
}
/*========================================================================*/
static void SceneDoRoving(PyMOLGlobals *G,float old_front,
                          float old_back,float old_origin,
                          int adjust_flag,int zoom_flag)
{
  EditorFavorOrigin(G,NULL);
  if((int)SettingGet(G,cSetting_roving_origin)) {

    register CScene *I=G->Scene;
    float delta_front,delta_back;
    float front_weight,back_weight,slab_width;
    float z_buffer = 3.0;
    float old_pos2 = 0.0F;
    float v2[3];

    z_buffer = SettingGet(G,cSetting_roving_origin_z_cushion);
    
    delta_front = I->Front - old_front;
    delta_back = I->Back - old_back;
    
    zero3f(v2);
    
    slab_width = I->Back - I->Front;
    
    /* first, check to make sure that the origin isn't too close to either plane */
    if((z_buffer*2)>slab_width)
      z_buffer = slab_width*0.5F;      
    
    if(old_origin<(I->Front+z_buffer)) { /* old origin behind front plane */
      front_weight = 1.0F;
      delta_front = (I->Front+z_buffer)-old_origin; /* move origin into allowed regioin */
    } else if(old_origin>(I->Back-z_buffer)) { /* old origin was behind back plane */
      front_weight = 0.0F;
      delta_back = (I->Back-z_buffer)-old_origin;
      
    } else if(slab_width>=R_SMALL4) { /* otherwise, if slab exists */
      front_weight = (old_back-old_origin)/slab_width; /* weight based on relative proximity */
    } else {
      front_weight = 0.5F;
    }
    
    back_weight = 1.0F - front_weight;
    
    if((front_weight>0.2)&&(back_weight>0.2)) { /* origin not near edge */
      if(delta_front*delta_back>0.0F) { /* planes moving in same direction */
        if(fabs(delta_front)>fabs(delta_back)) { /* so stick with whichever moves less */
          v2[2] = delta_back;
        } else {
          v2[2] = delta_front;
        }
      } else {
        /* planes moving in opposite directions (increasing slab size) */
        /* don't move origin */
      }
    } else { /* origin is near edge -- move origin with plane having highest weight */
      if(front_weight<back_weight) {
        v2[2] = delta_back;
      } else {
        v2[2] = delta_front;
      }
    }

    old_pos2 = I->Pos[2];

    MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); /* transform offset into realspace */
    subtract3f(I->Origin,v2,v2); /* calculate new origin location */
    SceneOriginSet(G,v2,true); /* move origin, preserving camera location */
    
    if(SettingGet(G,cSetting_ortho) || zoom_flag) { 
      /* we're orthoscopic, so we don't want the effective field of view 
         to change.  Thus, we have to hold Pos[2] constant, and instead
         move the planes.
      */
      float delta = old_pos2-I->Pos[2];
      I->Pos[2] += delta;
      SceneClipSet(G, I->Front - delta, I->Back - delta );
    }
    slab_width = I->Back - I->Front;
    
    /* first, check to make sure that the origin isn't too close to either plane */
    if((z_buffer*2)>slab_width)
      z_buffer = slab_width*0.5F;      
    
  }
  if((adjust_flag)&&(int)SettingGet(G,cSetting_roving_detail)) {    
    SceneRovingPostpone(G);
  }
  if(SettingGet(G,cSetting_roving_detail)) {    
    SceneRovingDirty(G);
  }
}

#define cDoubleTime 0.35

static int SceneDoXYPick(PyMOLGlobals *G, int x, int y, int click_side)
{
  CScene *I=G->Scene;
  int defer_builds_mode = SettingGetGlobal_i(G,cSetting_defer_builds_mode);

  if(defer_builds_mode==5) /* force generation of a pickable version */
    SceneUpdate(G,true);

  if(OrthoGetOverlayStatus(G)||SettingGetGlobal_i(G,cSetting_text))
    SceneRender(G,NULL,0,0,NULL,0,0,0,0); /* remove overlay if present */
  SceneDontCopyNext(G);
  
  I->LastPicked.context.object = NULL;
  SceneRender(G,&I->LastPicked,x,y,NULL,0,0,click_side,0);
  return (I->LastPicked.context.object!=NULL);
  /* did we pick something? */
}


static void SceneNoteMouseInteraction(PyMOLGlobals *G)
{
  SceneAbortAnimation(G);
  if(SettingGet_b(G,NULL,NULL,cSetting_mouse_restart_movie_delay)) {
    SceneRestartFrameTimer(G);
  }
}

/*========================================================================*/
static int SceneClick(Block *block,int button,int x,int y,
                      int mod,double when)
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  CObject *obj;
  ObjectMolecule *objMol;
  OrthoLineType buffer,buf1,buf2;
  WordType selName = "";
  register int mode = 0; /* trying to work around something... */
  int atIndex;
  char empty_string[1] = "";
  char *sel_mode_kw = empty_string;
  int is_single_click = (( button == P_GLUT_SINGLE_LEFT   ) ||
                         ( button == P_GLUT_SINGLE_MIDDLE ) ||
                         ( button == P_GLUT_SINGLE_RIGHT  ));
  int click_handled = false;
  int click_side = 0;

  if(!is_single_click) {
    int click_handled = false;

    if(I->ButtonsShown) {
      
      int i;
      SceneElem *elem = I->SceneVLA;
      
      if(I->ScrollBarActive) {
        if((x-I->Block->rect.left)<(SceneScrollBarWidth+SceneScrollBarMargin)) {
          click_handled = true;
          ScrollBarDoClick(I->ScrollBar,button,x,y,mod);      
        }
      } 
      if(!click_handled) {
        for(i=0;i<I->NScene;i++) {
          if(elem->drawn && 
             (x>=elem->x1) &&
             (y>=elem->y1) &&
             (x<elem->x2) &&
             (y<elem->y2)) {
            click_handled = true;
            break;
          }
          elem++;
        }
      }
    }

    if(!click_handled) {
      
      if( ((ButModeCheckPossibleSingleClick(G,button,mod) || (!mod))
	   &&((when-I->LastClickTime)<cDoubleTime))) {
	int dx,dy;
	dx = abs(I->LastWinX - x);
	dy = abs(I->LastWinY - y);
	if((dx<10)&&(dy<10)&&(I->LastButton==button)) {
	  switch(button) {
	  case P_GLUT_LEFT_BUTTON:
	    button = P_GLUT_DOUBLE_LEFT;
	    break;
	  case P_GLUT_MIDDLE_BUTTON:
	    button = P_GLUT_DOUBLE_MIDDLE;
	    break;
	  case P_GLUT_RIGHT_BUTTON:
	    button = P_GLUT_DOUBLE_RIGHT;
	    break;
	  }
	}
      }
    }
    
    if(ButModeCheckPossibleSingleClick(G,button,mod) || (!mod)) {
      I->PossibleSingleClick = 1;
    } else {
      char *but_mode_name = SettingGetGlobal_s(G,cSetting_button_mode_name);
      if(but_mode_name && but_mode_name[0]=='1') {
        I->PossibleSingleClick = 1;
      } else {
        I->PossibleSingleClick = 0;
      }
    }
  }

  I->LastWinX = x;
  I->LastWinY = y;
  I->LastClickTime = when;
  I->LastButton = button;
  I->LastMod = mod;
  I->Threshold = 0;
  
  if(I->ButtonsShown) {
    int i;
    SceneElem *elem = I->SceneVLA;
    
    if(I->ScrollBarActive) {
      if((x-I->Block->rect.left)<(SceneScrollBarWidth+SceneScrollBarMargin)) {
        click_handled = true;
        ScrollBarDoClick(I->ScrollBar,button,x,y,mod);      
      }
    } 
    if(!click_handled) {
      for(i=0;i<I->NScene;i++) {
        if(elem->drawn && 
           (x>=elem->x1) &&
           (y>=elem->y1) &&
           (x<elem->x2) &&
           (y<elem->y2)) {
          switch(button) {
          case P_GLUT_LEFT_BUTTON: /* normal activate (with interpolation) */
            I->Pressed = i;
            I->Over = i;
            I->PressMode = 1;
            SceneDirty(G);
            click_handled = true;
            break;
          case P_GLUT_MIDDLE_BUTTON: /* rapid browse mode */
            I->Pressed = i;
            I->PressMode = 2;
            I->Over = i;
            click_handled = true;
            {
              char *cur_name = SettingGetGlobal_s(G,cSetting_scene_current_name);
              int animate=-1;
              if(mod&cOrthoCTRL) animate=0;
              if(cur_name && elem->name && (strcmp(cur_name, elem->name))) {
                OrthoLineType buffer;
                sprintf(buffer,"cmd.scene('''%s''',animate=%d)",elem->name,animate);
                PParse(G,buffer);
                PFlush(G);
                PLog(G,buffer,cPLog_pym);
              }
            }
            break;
          case P_GLUT_RIGHT_BUTTON: /* drag or menu... */
            I->Pressed = i;
            I->PressMode = 3;
            I->Over = i;
            click_handled = true;
            break;
          }
          break;
        }
        elem++;
      }
    }
  }
  if(!click_handled) {

    mode = ButModeTranslate(G,button,mod); 

    I->Button=button;    
    I->SculptingSave = 0;
    switch(mode) {
    case cButModeScaleSlabExpand:
      SceneNoteMouseInteraction(G);
      SceneClip(G,5,1.0F+(0.2*SettingGetGlobal_f(G,cSetting_mouse_wheel_scale)),NULL,0);
      break;
    case cButModeScaleSlabShrink:
      SceneNoteMouseInteraction(G);
      SceneClip(G,5,1.0F-(0.2*SettingGetGlobal_f(G,cSetting_mouse_wheel_scale)),NULL,0);
      break;
    case cButModeMoveSlabForward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->Front;
        float old_back = I->Back;
        float old_origin = -I->Pos[2];
        SceneClip(G,6,0.1F*SettingGetGlobal_f(G,cSetting_mouse_wheel_scale),NULL,0);
        SceneDoRoving(G,old_front,old_back,old_origin,true,false);
      }
      break;
    case cButModeMoveSlabBackward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->Front;
        float old_back = I->Back;
        float old_origin = -I->Pos[2];

        SceneClip(G,6,-0.1F*SettingGetGlobal_f(G,cSetting_mouse_wheel_scale),NULL,0);
        SceneDoRoving(G,old_front,old_back,old_origin,true,false);
      }
      break;
    case cButModeZoomForward:
      SceneNoteMouseInteraction(G);
      {
        float factor = -((I->FrontSafe+I->BackSafe)/2)*0.1*
          SettingGetGlobal_f(G,cSetting_mouse_wheel_scale);
        if(factor<=0.0F) {
          I->Pos[2]+=factor;
          I->Front-=factor;
          I->Back-=factor;
          I->FrontSafe = GetFrontSafe(I->Front,I->Back);
          I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
        }
      }
      break;
    case cButModeZoomBackward:
      SceneNoteMouseInteraction(G);
      {
        float factor = ((I->FrontSafe+I->BackSafe)/2)*0.1F
          *SettingGetGlobal_f(G,cSetting_mouse_wheel_scale);
        if(factor>=0.0F) {
          I->Pos[2]+=factor;
          I->Front-=factor;
          I->Back-=factor;
          I->FrontSafe = GetFrontSafe(I->Front,I->Back);
          I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
        }
      }
      break;
    case cButModeMoveSlabAndZoomForward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->Front;
        float old_back = I->Back;
        float old_origin = -I->Pos[2];
        SceneClip(G,6,0.1F*SettingGetGlobal_f(G,cSetting_mouse_wheel_scale),NULL,0);
        SceneDoRoving(G,old_front,old_back,old_origin,true,true);
      }
      break;
    case cButModeMoveSlabAndZoomBackward:
      SceneNoteMouseInteraction(G);
      {
        float old_front = I->Front;
        float old_back = I->Back;
        float old_origin = -I->Pos[2];
        SceneClip(G,6,-0.1F*SettingGetGlobal_f(G,cSetting_mouse_wheel_scale),NULL,0);
        SceneDoRoving(G,old_front,old_back,old_origin,true,true);
      }
      break;
    case cButModeRectAdd: /* deprecated */
    case cButModeRectSub:/* deprecated */
    case cButModeRect:/* deprecated */
    case cButModeSeleAddBox:
    case cButModeSeleSetBox:
    case cButModeSeleSubBox:
      return SceneLoopClick(block,button,x,y,mod);
      break;
    case cButModeRotDrag:
    case cButModeMovDrag:
    case cButModeMovDragZ:
      SceneNoteMouseInteraction(G);
      SceneDontCopyNext(G);

      y=y-I->Block->margin.bottom;
      x=x-I->Block->margin.left;
    
      if(stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x,NULL,I->Width, NULL);

      I->LastX=x;
      I->LastY=y;     
      EditorReadyDrag(G,SettingGetGlobal_i(G,cSetting_state)-1);
      break;
    case cButModeRotXYZ:
    case cButModeTransXY:
    case cButModeTransZ:
    case cButModeClipNF:
    case cButModeClipN:    
    case cButModeClipF:    
    case cButModeRotZ:
    case cButModeInvRotZ:
      SceneNoteMouseInteraction(G);
      SceneDontCopyNext(G);

      y=y-I->Block->margin.bottom;
      x=x-I->Block->margin.left;
    
      if(stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x,NULL,I->Width, NULL);

      I->LastX=x;
      I->LastY=y;     
      break;
    case cButModePickAtom1:
    case cButModePickAtom:
    case cButModeMenu:
      if(stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x,NULL,I->Width, NULL);

      if(SceneDoXYPick(G,x,y, click_side)) {
        obj=(CObject*)I->LastPicked.context.object;
        y=y-I->Block->margin.bottom;
        x=x-I->Block->margin.left;
        I->LastX=x;
        I->LastY=y;    
        switch(obj->type) {
        case cObjectMolecule:
          switch(mode) {
          case cButModeMenu:
            {
              ObjectMolecule *objMol = (ObjectMolecule*)obj;
              int active_sele = ExecutiveGetActiveSele(G);
              if(active_sele && SelectorIsMember(G,objMol->AtomInfo[I->LastPicked.src.index].selEntry,
                                                 active_sele)) {
                ObjectNameType name;
                ExecutiveGetActiveSeleName(G,name,false,SettingGet(G,cSetting_logging));
                MenuActivate2Arg(G,I->LastWinX,I->LastWinY+20, /* selection menu */
                                 I->LastWinX,I->LastWinY,
                                 is_single_click,
                                 "pick_sele",name,name);
              } else {
                ObjectMoleculeGetAtomSele((ObjectMolecule*)obj,I->LastPicked.src.index,buffer);
                ObjectMoleculeGetAtomSeleLog((ObjectMolecule*)obj,I->LastPicked.src.index,buf1,false);
                MenuActivate2Arg(G,I->LastWinX,I->LastWinY+20,
                                 I->LastWinX,I->LastWinY,
                                 is_single_click,
                                 "pick_menu",buffer,buf1);
              }
            }
            break;
          case cButModePickAtom1:
            if(obj&&obj->type==cObjectMolecule) {
              if(Feedback(G,FB_Scene,FB_Results)) {
                if(obj->fDescribeElement)
                  obj->fDescribeElement(obj,I->LastPicked.src.index,buffer);
                PRINTF " You clicked %s -> (%s)\n",buffer,cEditorSele1 ENDF(G);
              }
              if(SettingGet(G,cSetting_logging)) {
                objMol = (ObjectMolecule*)obj;            
                ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buffer,false);
                sprintf(buf2,"cmd.edit(\"%s\",pkresi=1)",buffer);
                PLog(G,buf2,cPLog_pym);
              }
              OrthoRestorePrompt(G);
              sprintf(buffer,"%s`%d",
                      obj->Name,I->LastPicked.src.index+1);    
              EditorInactivate(G);
              SelectorCreate(G,cEditorSele1,buffer,NULL,true,NULL);
              EditorActivate(G,SettingGetGlobal_i(G,cSetting_state)-1,false);
              if(EditorActive(G)) {
                EditorDefineExtraPks(G);
              }
              WizardDoPick(G,0);
            }
            break;
          case cButModePickAtom:
            if(obj&&obj->type==cObjectMolecule){
              WordType name;
              if(obj->fDescribeElement)
                obj->fDescribeElement(obj,I->LastPicked.src.index,buffer);
              if(EditorIsBondMode(G)
                 /* &&!(EditorIsAnActiveObject(G,(ObjectMolecule*)obj))*/ ) {
                EditorInactivate(G);
                EditorLogState(G,false);
              }
              if((!EditorIsBondMode(G))&&
                 EditorDeselectIfSelected(G,
                                          (ObjectMolecule*)obj,I->LastPicked.src.index,true)) {
                PRINTF " You unpicked %s.",buffer ENDF(G);
                if(EditorActive(G)) 
                  EditorDefineExtraPks(G);
                EditorLogState(G,false);
              } else {
                if(EditorIsBondMode(G)&&
                   EditorDeselectIfSelected(G,
                                            (ObjectMolecule*)obj,I->LastPicked.src.index,false)) {
                  EditorInactivate(G);
                }
                EditorGetNextMultiatom(G,name);

                PRINTFB(G,FB_Scene,FB_Results) " You clicked %s -> (%s)\n",buffer,name ENDFB(G);
                /* TODO: logging */
              
                sprintf(buffer,"%s`%d",obj->Name,I->LastPicked.src.index+1);    
                ExecutiveDelete(G,name);
                SelectorCreate(G,name,buffer,NULL,true,NULL);
                EditorActivate(G,SettingGetGlobal_i(G,cSetting_state)-1,false);
                if(EditorActive(G)) {
                  EditorDefineExtraPks(G);
                }
                EditorLogState(G,false);
                WizardDoPick(G,0);
              }
            }
            break;
          }
          break;
        case cObjectGadget:
          break;
        default:
          EditorInactivate(G);
          break;
        }
      } else { /* no atom picked */
        switch(mode) {
        case cButModeMenu:
          MenuActivate0Arg(G,I->LastWinX,I->LastWinY,
                           I->LastWinX,I->LastWinY,
                           is_single_click,"main_menu");
          break;
        default:
          EditorInactivate(G);
          if(SettingGet(G,cSetting_logging)) {
            PLog(G,"cmd.edit()",cPLog_pym);
          }
          break;
        }
      }
      SceneDirty(G);
      break;
    case cButModePickBond:
    case cButModePkTorBnd:
      if(stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x,NULL,I->Width, &click_side);

      if(SceneDoXYPick(G,x,y,click_side)) {
        obj=(CObject*)I->LastPicked.context.object;
        y=y-I->Block->margin.bottom;
        x=x-I->Block->margin.left;
        I->LastX=x;
        I->LastY=y;    

        if(mode==cButModePkTorBnd) {
          I->Threshold = 3;
          I->ThresholdX = x;
          I->ThresholdY = y;
        }

        switch(obj->type) {
        case cObjectMolecule:

          EditorInactivate(G);
          if(Feedback(G,FB_Scene,FB_Results)) {
            if(obj->fDescribeElement)
              obj->fDescribeElement(obj,I->LastPicked.src.index,buffer);
            PRINTF " You clicked %s -> (%s)",buffer,cEditorSele1 ENDF(G);
            OrthoRestorePrompt(G);
          }

          /*        ObjectMoleculeChooseBondDir(objMol,I->LastPicked.bond,
                    &I->LastPicked.src.index,&atIndex);*/
        
          sprintf(buffer,"%s`%d",
                  obj->Name,I->LastPicked.src.index+1);    
          SelectorCreate(G,cEditorSele1,buffer,NULL,true,NULL);
          objMol = (ObjectMolecule*)obj;
          if(I->LastPicked.src.bond>=0) {
            atIndex = objMol->Bond[I->LastPicked.src.bond].index[0];
            if(atIndex == I->LastPicked.src.index)
              atIndex = objMol->Bond[I->LastPicked.src.bond].index[1];       
            if(Feedback(G,FB_Scene,FB_Results)) {
              if(obj->fDescribeElement)
                obj->fDescribeElement(obj,atIndex,buffer);
              PRINTF " You clicked %s -> (%s)",buffer,cEditorSele2 ENDF(G);
              OrthoRestorePrompt(G);
            }
          
            if(SettingGet(G,cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buf1,false);
              ObjectMoleculeGetAtomSeleLog(objMol,atIndex,buf2,false);
              sprintf(buffer,"cmd.edit(\"%s\",\"%s\")",buf1,buf2);
              PLog(G,buffer,cPLog_pym);
            }
            sprintf(buffer,"%s`%d",
                    obj->Name,atIndex+1);    
            SelectorCreate(G,cEditorSele2,buffer,NULL,true,NULL);
            EditorActivate(G,SettingGetGlobal_i(G,cSetting_state)-1,true);


            if(mode==cButModePkTorBnd) {
              /* get ready to drag */
              SceneDontCopyNext(G);
              switch(obj->type) {
              case cObjectMolecule:
                objMol = (ObjectMolecule*)obj;
                EditorPrepareDrag(G,objMol,-1,I->LastPicked.src.index,
                                  SettingGetGlobal_i(G,cSetting_state)-1, mode);
                I->SculptingFlag = 1;
                I->SculptingSave =  objMol->AtomInfo[I->LastPicked.src.index].protekted;
                objMol->AtomInfo[I->LastPicked.src.index].protekted=2;
                break;
              }
            }
            WizardDoPick(G,1);
          } else {
            WizardDoPick(G,0);
          }
          if(SettingGet(G,cSetting_auto_hide_selections))
            ExecutiveHideSelections(G);
          break;
        case cObjectGadget:
          break;
        default:
          EditorInactivate(G);
          break;
        }
      } else {
        EditorInactivate(G);
        EditorLogState(G,false); 
      }
      SceneInvalidate(G);
      break;
    case cButModeRotObj:
    case cButModeMovObj:
    case cButModeMovObjZ:
    case cButModeRotView:
    case cButModeMovView:
    case cButModeMovViewZ:
    case cButModeRotFrag:
    case cButModeMovFrag:
    case cButModeMovFragZ:
    case cButModeTorFrag:
    case cButModeMoveAtom:
    case cButModeMoveAtomZ:
      if(stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x,NULL,I->Width, &click_side);

      if(SceneDoXYPick(G,x,y,click_side)) {
        obj=(CObject*)I->LastPicked.context.object;
        y=y-I->Block->margin.bottom;
        x=x-I->Block->margin.left;
        I->LastX=x;
        I->LastY=y;    
        switch(obj->type) {
        case cObjectMolecule:

          if(I->LastPicked.src.bond==cPickableLabel) {        
            /* if user picks a label with move object/move fragment,
               then move the object/fragment, not the label */

            switch(mode) {
            case cButModeRotObj:
            case cButModeMovObj:
            case cButModeMovObjZ:
            case cButModeRotFrag:
            case cButModeMovFrag:
            case cButModeMovFragZ:
              I->LastPicked.src.bond=cPickableAtom;
              break;
            }
          }

          if(I->LastPicked.src.bond>=cPickableAtom) {
            if(Feedback(G,FB_Scene,FB_Results)) {
              if(obj->fDescribeElement) 
                obj->fDescribeElement(obj,I->LastPicked.src.index,buffer);
              PRINTF " You clicked %s",buffer ENDF(G);        
              OrthoRestorePrompt(G);
            }
          }
          objMol = (ObjectMolecule*)obj;
          EditorPrepareDrag(G,objMol,-1,I->LastPicked.src.index,
                            SettingGetGlobal_i(G,cSetting_state)-1, mode);

          if(I->LastPicked.src.bond>=cPickableAtom) {
            I->SculptingFlag = 1;
            I->SculptingSave =  objMol->AtomInfo[I->LastPicked.src.index].protekted;
            objMol->AtomInfo[I->LastPicked.src.index].protekted=2;
          }
          break;
        case cObjectSlice:
          if(ObjectSliceGetVertex((ObjectSlice*)obj,I->LastPicked.src.index,
                                  I->LastPicked.src.bond,I->LastPickVertex)) {
            I->LastPickVertexFlag=true;
          }
          break;
        case cObjectMeasurement:
          break;
        case cObjectGadget:
          break;
        default:
          EditorInactivate(G);
          break;
        }
        /*
          (int)SettingGet(G,cSetting_sculpting);
          SettingSet(G,cSetting_sculpting,0);*/
      }
      break;
 
    case cButModeSeleSet:
    case cButModeSeleToggle:
      sel_mode_kw = SceneGetSeleModeKeyword(G);
         
      /* intentional pass through */

    case cButModeLB:
    case cButModeMB:
    case cButModeRB:
    case cButModeAddToLB:
    case cButModeAddToMB:
    case cButModeAddToRB:
    case cButModeSimpleClick:
    case cButModeOrigAt:
    case cButModeCent:
    case cButModeDragMol:
    case cButModeDragObj:
      if(stereo_via_adjacent_array(I->StereoMode))
        x = get_stereo_x(x,NULL,I->Width, &click_side);

      if(SceneDoXYPick(G,x,y,click_side)) {
        obj=(CObject*)I->LastPicked.context.object;

        switch(obj->type) {
        case cObjectMolecule:
          if(Feedback(G,FB_Scene,FB_Results)) {
            if(obj->fDescribeElement) 
              obj->fDescribeElement(obj,I->LastPicked.src.index,buffer);
            PRINTF " You clicked %s",buffer ENDF(G);        
            OrthoRestorePrompt(G);
          }
          sprintf(buffer,"%s`%d",
                  obj->Name,I->LastPicked.src.index+1);
          switch(mode) {
          case cButModeLB:
          case cButModeAddToLB:
            strcpy(selName,"lb");
            break;
          case cButModeMB:
          case cButModeAddToMB:
            strcpy(selName,"mb");
            break;
          case cButModeRB:
          case cButModeAddToRB:
            strcpy(selName,"rb");
            break;
          case cButModeSeleSet:
          case cButModeSeleToggle:
            ExecutiveGetActiveSeleName(G,selName,true,SettingGet(G,cSetting_logging));
            break;
          case cButModeDragMol:
            {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buf1,false);
              sprintf(buffer,"cmd.drag(\"bymol (%s)\")",buf1);
              PParse(G,buffer);
              PLog(G,buffer,cPLog_pym);
            }
            break;
          case cButModeDragObj:
            {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buf1,false);
              sprintf(buffer,"cmd.drag(\"byobject (%s)\")",buf1);
              PParse(G,buffer);
              PLog(G,buffer,cPLog_pym);
            }
            break;
          case cButModeOrigAt:
            SceneNoteMouseInteraction(G);
            {
              float v1[3];

              if(ObjectMoleculeGetAtomTxfVertex((ObjectMolecule*)obj,
                                                SettingGetGlobal_i(G,cSetting_state)-1,
                                                I->LastPicked.src.index,v1)) {
                EditorFavorOrigin(G,v1);
                ExecutiveOrigin(G,NULL,true,NULL,v1,0);
              }
            }
            if(obj->type==cObjectMolecule) {
              if(SettingGet(G,cSetting_logging)) {
                objMol = (ObjectMolecule*)obj;            
                ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buf1,false);
                sprintf(buffer,"cmd.origin(\"%s\")",buf1);
                PLog(G,buffer,cPLog_pym);

              }
              if(Feedback(G,FB_Scene,FB_Results)) {
                if(obj->fDescribeElement) 
                  obj->fDescribeElement(obj,I->LastPicked.src.index,buffer);
                PRINTF " You clicked %s",buffer ENDF(G);        
                OrthoRestorePrompt(G);
              }
            }
            PRINTFB(G,FB_Scene,FB_Actions) 
              " Scene: Origin set.\n"
              ENDFB(G);
            break;
          case cButModeCent:
            SceneNoteMouseInteraction(G);
            {
              float v1[3];

              if(ObjectMoleculeGetAtomTxfVertex((ObjectMolecule*)obj,
                                                SettingGetGlobal_i(G,cSetting_state)-1,
                                                I->LastPicked.src.index,v1)) {
                ExecutiveCenter(G,NULL,0,true,-1,v1,true);
              }
            }
          
            if(SettingGet(G,cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buf1,false);
              sprintf(buffer,"cmd.center(\"%s\",state=-1)",buf1);
              PLog(G,buffer,cPLog_pym);
            }
            break;
          }
          switch(mode) {
          case cButModeSimpleClick:
            PyMOL_SetClickReady(G->PyMOL,obj->Name,I->LastPicked.src.index,
                                button,mod,I->LastWinX,I->Height-(I->LastWinY+1));
            break;
          case cButModeLB:
          case cButModeMB:
          case cButModeRB:
          case cButModeSeleSet:
            sprintf(buf2,"(%s(%s))",sel_mode_kw,buffer);
            SelectorCreate(G,selName,buf2,NULL,false,NULL);
            if(SettingGet(G,cSetting_auto_hide_selections))
              ExecutiveHideSelections(G);
            if(SettingGet(G,cSetting_auto_show_selections))
              ExecutiveSetObjVisib(G,selName,1,false);
            if(obj->type==cObjectMolecule) {
              if(SettingGet(G,cSetting_logging)) {
                objMol = (ObjectMolecule*)obj;            
                ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buf1,false);
                sprintf(buffer,"cmd.select('%s',\"%s(%s)\",enable=1)",selName,sel_mode_kw,buf1);
                PLog(G,buffer,cPLog_pym);
              }
            }
            WizardDoSelect(G,selName);
            break;
          case cButModeAddToLB:
          case cButModeAddToMB:
          case cButModeAddToRB:
          case cButModeSeleToggle:
            if(SelectorIndexByName(G,selName)>=0) {
              sprintf(buf2,"(((%s) or %s(%s)) and not ((%s(%s)) and %s(%s)))",
                      selName,sel_mode_kw,buffer,sel_mode_kw,buffer,sel_mode_kw,selName);
              SelectorCreate(G,selName,buf2,NULL,false,NULL);
              if(obj->type==cObjectMolecule) {
                if(SettingGet(G,cSetting_logging)) {
                  objMol = (ObjectMolecule*)obj;            
                  ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buffer,false);
                  sprintf(buf2,"(((%s) or %s(%s)) and not ((%s(%s)) and %s(%s)))",
                          selName,sel_mode_kw,buffer,sel_mode_kw,buffer,sel_mode_kw,selName);
                  sprintf(buffer,"cmd.select('%s',\"%s(%s)\",enable=1)",selName,sel_mode_kw,buf2);
                  PLog(G,buffer,cPLog_pym);
                }
              }
            } else {
              sprintf(buf2,"%s(%s)",sel_mode_kw,buffer);
              SelectorCreate(G,selName,buf2,NULL,false,NULL);
              if(obj->type==cObjectMolecule) {
                if(SettingGet(G,cSetting_logging)) {
                  objMol = (ObjectMolecule*)obj;            
                  ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.src.index,buf1,false);
                  sprintf(buffer,"cmd.select('%s',\"%s(%s)\")",selName,sel_mode_kw,buf1);
                  PLog(G,buffer,cPLog_pym);
                }
              }
            }
            if(SettingGet(G,cSetting_auto_hide_selections))
              ExecutiveHideSelections(G);
            if(SettingGet(G,cSetting_auto_show_selections))
              ExecutiveSetObjVisib(G,selName,1,false);
            WizardDoSelect(G,selName);
            break;
          }
        case cObjectGadget:
          break;
        default:
          EditorInactivate(G);
          break;
        }
      } else {
        switch(mode) {
        case cButModeSeleSet:
          {
            OrthoLineType buf2;
            ObjectNameType name;

            if(ExecutiveGetActiveSeleName(G,name, false,SettingGet(G,cSetting_logging))) {
              SelectorCreate(G,name,"none",NULL,true,NULL);
              if(SettingGet(G,cSetting_logging)) {
                sprintf(buf2,"cmd.select('%s','none')\n",name);
                PLog(G,buf2,cPLog_no_flush);
              }
              SeqDirty(G);
            }
          }
        case cButModeSeleToggle:
          {
            OrthoLineType buf2;
            ObjectNameType name;

            if(ExecutiveGetActiveSeleName(G,name, false,SettingGet(G,cSetting_logging))) {
              ExecutiveSetObjVisib(G,name,0,false);
              if(SettingGet(G,cSetting_logging)) {
                sprintf(buf2,"cmd.disable('%s')\n",name);
                PLog(G,buf2,cPLog_no_flush);
              }
            }
          }
          break;
        case cButModeSimpleClick:
          PyMOL_SetClickReady(G->PyMOL,"",-1,button,mod,I->LastWinX,I->Height-(I->LastWinY+1));
          break;
        }
        PRINTFB(G,FB_Scene,FB_Blather) 
          " SceneClick: no atom found nearby.\n"
          ENDFB(G);
        SceneInvalidate(G); /* this here to prevent display weirdness after
                               an unsuccessful picking pass... not sure it helps though*/
        OrthoRestorePrompt(G);
      }
    }
  
    I->StartX = I->LastX;
    I->StartY = I->LastY;
  }
  return(1);
}
void ScenePushRasterMatrix(PyMOLGlobals *G,float *v) 
{
  register CScene *I=G->Scene;
  float scale = SceneGetExactScreenVertexScale(G,v);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glTranslatef(v[0],v[1],v[2]); /* go to this position */
  glMultMatrixf(I->InvMatrix);
  glScalef(scale,scale,scale);
  /*  glTranslatef(-0.33F,-0.33F,0.0F); */
  
}

void ScenePopRasterMatrix(PyMOLGlobals *G)
{
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}
/*========================================================================*/
void SceneGetEyeNormal(PyMOLGlobals *G,float *v1,float *normal)
{
  register CScene *I=G->Scene;
  float p1[4],p2[4];
  float modelView[16];
  
  identity44f(modelView);
  MatrixTranslateC44f(modelView,I->Pos[0],I->Pos[1],I->Pos[2]);
  MatrixMultiplyC44f(I->RotMatrix,modelView);
  MatrixTranslateC44f(modelView,-I->Origin[0],-I->Origin[1],-I->Origin[2]);

  copy3f(v1,p1);
  p1[3] = 1.0;
#if 1
  MatrixTransformC44f4f(modelView,p1,p2); /* modelview transformation */
#else
  MatrixTransformC44f4f(I->ModMatrix,p1,p2); /* modelview transformation */
#endif
  copy3f(p2,p1);
  normalize3f(p1);
  MatrixInvTransformC44fAs33f3f(I->RotMatrix,p1,p2); 
  invert3f3f(p2,normal);
}
/*========================================================================*/
float SceneGetScreenVertexScale(PyMOLGlobals *G,float *v1) 
     /* does not require OpenGL-provided matrices */
{
  float ratio;
  register CScene *I=G->Scene;
  float vt[3];
  float fov=SettingGet(G,cSetting_field_of_view);
  float modelView[16];

  if(!v1) v1 = I->Origin;
  identity44f(modelView);
  MatrixTranslateC44f(modelView,I->Pos[0],I->Pos[1],I->Pos[2]);
  MatrixMultiplyC44f(I->RotMatrix,modelView);
  MatrixTranslateC44f(modelView,-I->Origin[0],-I->Origin[1],-I->Origin[2]);

  MatrixTransformC44f3f(modelView, v1,vt);
  if(SettingGetGlobal_i(G,cSetting_ortho)) {
    ratio = 2*(float)(fabs(I->Pos[2])*tan((fov/2.0)*cPI/180.0))/(I->Height); 
  } else {
    float front_size = 2*I->FrontSafe*((float)tan((fov/2.0F)*PI/180.0F))/(I->Height);
    ratio = front_size*(-vt[2]/I->FrontSafe);
  }
  return ratio;
}

static float SceneGetExactScreenVertexScale(PyMOLGlobals *G,float *v1)
     /* requires the OpenGL-computed matrices */
{
  /* get conversion factor from screen point to atomic coodinate */
  register CScene *I=G->Scene;
  float vl,p1[4],p2[4];
  float height_factor = I->Height/2.0F;

  return SceneGetScreenVertexScale(G,v1);
  if(!v1) v1 = I->Origin;

  /* now, scale properly given the current projection matrix */
  copy3f(v1,p1);
  p1[3] = 1.0;
  MatrixTransformC44f4f(I->ModMatrix,p1,p2); /* modelview transformation */
  copy4f(p2,p1);
  p2[1]+=1.0;
  /* projection transformation */
  MatrixTransformC44f4f(I->ProMatrix,p1,p1); 
  MatrixTransformC44f4f(I->ProMatrix,p2,p2);
  /*  dump44f(I->ProMatrix,"pro");*/
  /* perspective division */
  p1[1]=p1[1]/p1[3];
  p2[1]=p2[1]/p2[3];
  p1[1]=(p1[1]+1.0F)*(height_factor); /* viewport transformation */
  p2[1]=(p2[1]+1.0F)*(height_factor);
  vl = (float)fabs(p1[1]-p2[1]);

  if(vl<R_SMALL4)
    vl=100.0F;

  return(1.0F/vl);
}

void SceneRovingChanged(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;  
  SceneRovingDirty(G);
  I->RovingCleanupFlag=true;
}

static void SceneRovingCleanup(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;  
  char *s;  
  char buffer[OrthoLineLength];

  I->RovingCleanupFlag=false;

  s = SettingGet_s(G,NULL,NULL,cSetting_roving_selection);

  sprintf(buffer,"cmd.hide('lines','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
  sprintf(buffer,"cmd.hide('sticks','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
  sprintf(buffer,"cmd.hide('spheres','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
  sprintf(buffer,"cmd.hide('ribbon','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
  sprintf(buffer,"cmd.hide('cartoon','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
  sprintf(buffer,"cmd.hide('labels','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
  sprintf(buffer,"cmd.hide('nonbonded','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
  sprintf(buffer,"cmd.hide('nb_spheres','''%s''')",s);
  PParse(G,buffer);
  PFlush(G);
}

void SceneRovingUpdate(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  char buffer[OrthoLineLength];
  float sticks,lines,spheres,labels,ribbon,cartoon;
  float polar_contacts,polar_cutoff,nonbonded,nb_spheres;
  char byres[10] = "byres";
  char not[4] = "not";
  char empty[1] = "";
  char *p1;
  char *p2;
  char *s;
  int refresh_flag=false;
  char *name;
  float level;
  float isosurface,isomesh;

  if(I->RovingDirtyFlag&&(
                          (UtilGetSeconds(G)-I->RovingLastUpdate)>
                          fabs(SettingGet(G,cSetting_roving_delay)))) {
    
    if(I->RovingCleanupFlag)
      SceneRovingCleanup(G);
    
    s = SettingGet_s(G,NULL,NULL,cSetting_roving_selection);
    sticks = SettingGet(G,cSetting_roving_sticks);
    lines = SettingGet(G,cSetting_roving_lines);
    labels = SettingGet(G,cSetting_roving_labels);
    spheres = SettingGet(G,cSetting_roving_spheres);
    ribbon = SettingGet(G,cSetting_roving_ribbon);
    cartoon = SettingGet(G,cSetting_roving_cartoon);
    polar_contacts = SettingGet(G,cSetting_roving_polar_contacts);
    polar_cutoff = SettingGet(G,cSetting_roving_polar_cutoff);
    nonbonded = SettingGet(G,cSetting_roving_nonbonded);
    nb_spheres = SettingGet(G,cSetting_roving_nb_spheres);

    isomesh = SettingGet(G,cSetting_roving_isomesh);
    isosurface = SettingGet(G,cSetting_roving_isosurface);

    if(SettingGet(G,cSetting_roving_byres))
      p2 = byres;
    else
      p2 = empty;

    if(sticks!=0.0F) {
      if(sticks<0.0F) {
        p1=not;
        sticks=(float)fabs(sticks);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('sticks','''%s''');cmd.show('sticks','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,sticks);
      PParse(G,buffer);
      PFlush(G);
      refresh_flag=true;
    }

    if(lines!=0.0F) {
      if(lines<0.0F) {
        p1=not;
        lines=(float)fabs(lines);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('lines','''%s''');cmd.show('lines','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,lines);
      PParse(G,buffer);
      PFlush(G);
      refresh_flag=true;
    }

    if(labels!=0.0F) {
      if(labels<0.0F) {
        p1=not;
        labels=(float)fabs(labels);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('labels','''%s''');cmd.show('labels','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,labels);
      PParse(G,buffer);
      PFlush(G);
      refresh_flag=true;
    }

    if(spheres!=0.0F) {
      if(spheres<0.0F) {
        p1=not;
        spheres=(float)fabs(spheres);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('spheres','''%s''');cmd.show('spheres','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,spheres);
      PParse(G,buffer);
      PFlush(G);
      refresh_flag=true;
    }

    if(cartoon!=0.0F) {
      if(cartoon<0.0F) {
        p1=not;
        cartoon=(float)fabs(cartoon);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('cartoon','''%s''');cmd.show('cartoon','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,cartoon);
      PParse(G,buffer);
      PFlush(G);
      refresh_flag=true;
    }

    if(ribbon!=0.0F) {
      if(ribbon<0.0F) {
        p1=not;
        ribbon=(float)fabs(ribbon);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('ribbon','''%s''');cmd.show('ribbon','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,ribbon);
      PParse(G,buffer);
      PFlush(G);

      refresh_flag=true;
    }


    if(polar_contacts!=0.0F) {
      int label_flag=0;
      if(polar_contacts<0.0F) {
        p1=not;
        polar_contacts=(float)fabs(polar_contacts);
      } else {
        p1=empty;
      }
      if(polar_cutoff<0.0F) {
        label_flag=true;
        polar_cutoff=(float)fabs(polar_cutoff);
      }
      sprintf(buffer,
              "cmd.dist('rov_pc','%s & enabled & %s %s (center expand %1.3f)','same',%1.4f,mode=2,label=%d,quiet=2)",
              s,p1,p2,polar_contacts,polar_cutoff,label_flag);
      PParse(G,buffer);
      PFlush(G);

      refresh_flag=true;
    }

    if(nonbonded!=0.0F) {
      if(nonbonded<0.0F) {
        p1=not;
        nonbonded=(float)fabs(nonbonded);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('nonbonded','''%s''');cmd.show('nonbonded','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,nonbonded);
      PParse(G,buffer);
      PFlush(G);
      refresh_flag=true;
    }

    if(nb_spheres!=0.0F) {
      if(nb_spheres<0.0F) {
        p1=not;
        nb_spheres=(float)fabs(nb_spheres);
      } else {
        p1=empty;
      }
      sprintf(buffer,
              "cmd.hide('nb_spheres','''%s''');cmd.show('nb_spheres','%s & enabled & %s %s (center expand %1.3f)')",
              s,s,p1,p2,nb_spheres);
      PParse(G,buffer);
      PFlush(G);
      refresh_flag=true;
    }

    if(isomesh!=0.0F) {
      int auto_save;

      auto_save = (int)SettingGet(G,cSetting_auto_zoom);
      SettingSet(G,cSetting_auto_zoom,0);
      
      name = SettingGet_s(G,NULL,NULL,cSetting_roving_map1_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(G,name))
            {
              level = SettingGet(G,cSetting_roving_map1_level);
              sprintf(buffer,
                      "cmd.isomesh('rov_m1','%s',%8.6f,'center',%1.3f)",
                      name,level,isomesh);
              PParse(G,buffer);
              PFlush(G);
              refresh_flag=true;
            }

      name = SettingGet_s(G,NULL,NULL,cSetting_roving_map2_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(G,name))
            {
              level = SettingGet(G,cSetting_roving_map2_level);
              sprintf(buffer,
                      "cmd.isomesh('rov_m2','%s',%8.6f,'center',%1.3f)",
                      name,level,isomesh);
              PParse(G,buffer);
              PFlush(G);
              refresh_flag=true;
            }

      name = SettingGet_s(G,NULL,NULL,cSetting_roving_map3_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(G,name))
            {
              level = SettingGet(G,cSetting_roving_map3_level);
              sprintf(buffer,
                      "cmd.isomesh('rov_m3','%s',%8.6f,'center',%1.3f)",
                      name,level,isomesh);
              PParse(G,buffer);
              PFlush(G);
              refresh_flag=true;
            }
      SettingSet(G,cSetting_auto_zoom,(float)auto_save);            
    }

    if(isosurface!=0.0F) {
      int auto_save;

      auto_save = (int)SettingGet(G,cSetting_auto_zoom);
      SettingSet(G,cSetting_auto_zoom,0.0F);
      
      name = SettingGet_s(G,NULL,NULL,cSetting_roving_map1_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(G,name))
            {
              level = SettingGet(G,cSetting_roving_map1_level);
              sprintf(buffer,
                      "cmd.isosurface('rov_s1','%s',%8.6f,'center',%1.3f)",
                      name,level,isosurface);
              PParse(G,buffer);
              PFlush(G);
              refresh_flag=true;
            }

      name = SettingGet_s(G,NULL,NULL,cSetting_roving_map2_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(G,name))
            {
              level = SettingGet(G,cSetting_roving_map2_level);
              sprintf(buffer,
                      "cmd.isosurface('rov_s2','%s',%8.6f,'center',%1.3f)",
                      name,level,isosurface);
              PParse(G,buffer);
              PFlush(G);
              refresh_flag=true;
            }

      name = SettingGet_s(G,NULL,NULL,cSetting_roving_map3_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(G,name))
            {
              level = SettingGet(G,cSetting_roving_map3_level);
              sprintf(buffer,
                      "cmd.isosurface('rov_s3','%s',%8.6f,'center',%1.3f)",
                      name,level,isosurface);
              PParse(G,buffer);
              PFlush(G);
              refresh_flag=true;
            }
      SettingSet(G,cSetting_auto_zoom,(float)auto_save);            
    }


    if(refresh_flag) {
      PParse(G,"cmd.refresh()");
      PFlush(G);
    }

    I->RovingLastUpdate=UtilGetSeconds(G);
    I->RovingDirtyFlag=false;
  } 
}

/*========================================================================*/
static int SceneDrag(Block *block,int x,int y,int mod,double when)
{
  PyMOLGlobals *G = block->G;
  register CScene *I=G->Scene;
  float scale,vScale;
  float v1[3],v2[3],n1[3],n2[3],r1,r2,cp[3],v3[3];
  float dx,dy,dt;
  float axis[3],axis2[3],theta,omega;
  float old_front,old_back,old_origin;
  int mode;
  int eff_width;
  int moved_flag;
  int adjust_flag;
  int drag_handled = false;
  CObject *obj;

  if(I->PossibleSingleClick) {
    double slowest_single_click_drag = 0.15;
    if((when-I->LastClickTime)>slowest_single_click_drag) {
      I->PossibleSingleClick = 0;
    }
  }

  if(I->LoopFlag) {
    return SceneLoopDrag(block,x,y,mod);
  }
  if(I->ButtonsShown && I->PressMode) {
    if(!I->ButtonsValid) {
      SceneUpdateButtons(G);
      /* write the update code here */
    } 
    if(I->ButtonsValid) {
      SceneElem *elem = I->SceneVLA;
      int i;
      drag_handled = true;
      I->Over = -1;
      for(i=0;i<I->NScene;i++) {
        if(elem->drawn && 
           (x>=elem->x1) &&
           (y>=elem->y1) &&
           (x<elem->x2) &&
           (y<elem->y2)) {
          I->Over = i;
          OrthoDirty(G);
          break;
        }
        elem++;
      }
      switch(I->PressMode) {
      case 2:
        if(I->Over>=0) {
          if(I->Pressed!=I->Over) {
            char *cur_name = SettingGetGlobal_s(G,cSetting_scene_current_name);
            if(cur_name && elem->name && (strcmp(cur_name, elem->name))) {
              OrthoLineType buffer;
              int animate=-1;
              if(mod&cOrthoCTRL) animate=0;
              sprintf(buffer,"cmd.scene('''%s''',animate=%d)",elem->name,animate);
              PParse(G,buffer);
              PFlush(G);
              PLog(G,buffer,cPLog_pym);
            }
            I->Pressed = I->Over;
          }
        } else {
          I->Pressed = -1;
        }
      case 3:
        if((I->Over>=0)&&(I->Pressed!=I->Over))
          I->PressMode = 4; /* activate dragging */
        break;
      }

      if(I->PressMode == 4) { /* dragging */
        if((I->Over>=0)&&(I->Pressed!=I->Over)&&(I->Pressed>=0)) {

          SceneElem *pressed = I->SceneVLA + I->Pressed;
          OrthoLineType buffer;

          if(I->Over>0) { /* not over the first scene in list */
            SceneElem *first = elem-1;
            SceneElem *second = pressed;
            if(first >= pressed) {
              first = elem;
              second = pressed;
            }
            sprintf(buffer,"cmd.scene_order('''%s %s''')",
                    first->name,
                    second->name);
          } else {
            sprintf(buffer,"cmd.scene_order('''%s''',location='top')",
                    pressed->name);
          }
          PParse(G,buffer);
          PFlush(G);
          PLog(G,buffer,cPLog_pym);
          I->Pressed = I->Over;
          I->ButtonsValid = false;
        }
      }
    }
  }

  if(!drag_handled) {

    mode = ButModeTranslate(G,I->Button,mod);
  
    y=y-I->Block->margin.bottom;


    scale = (float)I->Height;
    if(scale > I->Width)
      scale = (float)I->Width;
    scale = 0.45F * scale;
    SceneInvalidateCopy(G,false);  
    SceneDontCopyNext(G);
    switch(mode) {
    case cButModePickAtom:
      obj=(CObject*)I->LastPicked.context.object;
      if(obj)
        switch(obj->type) {
        case cObjectGadget: 
          {
            ObjectGadget *gad;
            gad = (ObjectGadget*)obj;
          
            ObjectGadgetGetVertex(gad,I->LastPicked.src.index,I->LastPicked.src.bond,v1);
          
            vScale = SceneGetExactScreenVertexScale(G,v1);
            if(stereo_via_adjacent_array(I->StereoMode)) {
              x = get_stereo_x(x,&I->LastX,I->Width, NULL);
            }
          
            /* transform into model coodinate space */
            switch(obj->Context) {
            case 1:
              {
                float divisor;
                divisor = (float)I->Width;
                if(I->Height<I->Width)
                  divisor = (float)I->Height;
                v2[0] = (x-I->LastX)/divisor;
                v2[1] = (y-I->LastY)/divisor;
                v2[2] = 0;
              }
              break;
            default:
            case 0:
              v2[0] = (x-I->LastX)*vScale;
              v2[1] = (y-I->LastY)*vScale;
              v2[2] = 0;
              MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); 
              break;
            }
            add3f(v1,v2,v2);
            ObjectGadgetSetVertex(gad,I->LastPicked.src.index,I->LastPicked.src.bond,v2);
            if(gad->Obj.fUpdate)
              gad->Obj.fUpdate((CObject*)gad);
            SceneChanged(G);
          }
          break;
        }
      I->LastX=x;
      I->LastY=y;
      break;
    case cButModeRotDrag:
      eff_width = I->Width;
      if(stereo_via_adjacent_array(I->StereoMode)) {
        eff_width = I->Width/2;
        x = get_stereo_x(x,&I->LastX,I->Width, NULL);
      }
    
      if(SettingGet_b(G,NULL,NULL,cSetting_virtual_trackball)) {
        v1[0] = (float)(eff_width/2) - x;
        v1[1] = (float)(I->Height/2) - y;
      
        v2[0] = (float)(eff_width/2) - I->LastX;
        v2[1] = (float)(I->Height/2) - I->LastY;
      
      } else {
        v1[0] = (float)(I->LastX) - x;
        v1[1] = (float)(I->LastY) - y;
      
        v2[0] = 0;
        v2[1] = 0;
      }
    
      r1 = (float)sqrt1f(v1[0]*v1[0] + v1[1]*v1[1]);
      r2 = (float)sqrt1f(v2[0]*v2[0] + v2[1]*v2[1]);
    
      if(r1<scale) {
        v1[2] = (float)sqrt1f(scale*scale - r1*r1);
      } else {
        v1[2] = 0.0;
      }
    
      if(r2<scale) {
        v2[2] = (float)sqrt1f(scale*scale - r2*r2);
      } else {
        v2[2] = 0.0;
      }
      normalize23f(v1,n1);
      normalize23f(v2,n2);
      cross_product3f(n1,n2,cp);
      theta = (float)(SettingGet_f(G,NULL,NULL,cSetting_mouse_scale)*
                      2*180*asin(sqrt1f(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]))/cPI);
      dx = (v1[0]-v2[0]);
      dy = (v1[1]-v2[1]);
      dt = (float)(SettingGet_f(G,NULL,NULL,cSetting_mouse_limit)*sqrt1f(dx*dx+dy*dy)/scale);
    
      if(theta>dt) theta = dt;
    
      normalize23f(cp,axis);
    
      axis[2] = -axis[2];

      theta = theta/(1.0F+(float)fabs(axis[2]));
    
      /* transform into model coodinate space */
      MatrixInvTransformC44fAs33f3f(I->RotMatrix,axis,v2); 
      v1[0] =(float)(cPI*theta/180.0);
      EditorDrag(G,NULL,-1,mode,
                 SettingGetGlobal_i(G,cSetting_state)-1,v1,v2,v3);
      I->LastX=x;
      I->LastY=y;
      break;
    case cButModeMovDrag:
    case cButModeMovDragZ:
      if(I->Threshold) {
        if((abs(x-I->ThresholdX)>I->Threshold)||
           (abs(y-I->ThresholdY)>I->Threshold)) {
          I->Threshold = 0;
        }
      }
      if(!I->Threshold) {

        copy3f(I->Origin,v1);
        vScale = SceneGetExactScreenVertexScale(G,v1);
        if(stereo_via_adjacent_array(I->StereoMode)) {
          x = get_stereo_x(x,&I->LastX,I->Width, NULL);
        }
      
        if(mode==cButModeMovDragZ) {
          v2[0] = 0;
          v2[1] = 0;
          v2[2] = -(y-I->LastY)*vScale;
        } else {
          v2[0] = (x-I->LastX)*vScale;
          v2[1] = (y-I->LastY)*vScale;
          v2[2] = 0;
        }
      
        v3[0] = 0.0F;
        v3[1] = 0.0F;
        v3[2] = 1.0F;
      
        /* transform into model coodinate space */
        MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); 
        MatrixInvTransformC44fAs33f3f(I->RotMatrix,v3,v3); 

        EditorDrag(G,NULL,-1,mode,
                   SettingGetGlobal_i(G,cSetting_state)-1,v1,v2,v3);
      }
      I->LastX=x;
      I->LastY=y;
      break;
    case cButModeRotObj:
    case cButModeMovObj:
    case cButModeMovObjZ:
    case cButModeRotView:
    case cButModeMovView:
    case cButModeMovViewZ:
    case cButModeRotFrag:
    case cButModeMovFrag:
    case cButModeMovFragZ:
    case cButModeTorFrag:
    case cButModeMoveAtom:
    case cButModeMoveAtomZ:
    case cButModePkTorBnd:
      obj=(CObject*)I->LastPicked.context.object;
      if(obj) {
        if(I->Threshold) {
          if((abs(x-I->ThresholdX)>I->Threshold)||
             (abs(y-I->ThresholdY)>I->Threshold)) {
            I->Threshold = 0;
          }
        }
        if(!I->Threshold)
          switch(obj->type) {
          case cObjectGadget:  /* note repeated above */
            {
              ObjectGadget *gad;
              gad = (ObjectGadget*)obj;
            
              ObjectGadgetGetVertex(gad,I->LastPicked.src.index,I->LastPicked.src.bond,v1);
            
              vScale = SceneGetExactScreenVertexScale(G,v1);
              if(stereo_via_adjacent_array(I->StereoMode)) {
                x = get_stereo_x(x,&I->LastX,I->Width, NULL);
              }
            
              /* transform into model coodinate space */
              switch(obj->Context) {
              case 1:
                {
                  float divisor;
                  divisor = (float)I->Width;
                  if(I->Height<I->Width)
                    divisor = (float)I->Height;
                  v2[0] = (x-I->LastX)/divisor;
                  v2[1] = (y-I->LastY)/divisor;
                  v2[2] = 0;
                }
                break;
              default:
              case 0:
                v2[0] = (x-I->LastX)*vScale;
                v2[1] = (y-I->LastY)*vScale;
                v2[2] = 0;
                MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); 
                break;
              }
              add3f(v1,v2,v2);
              ObjectGadgetSetVertex(gad,I->LastPicked.src.index,I->LastPicked.src.bond,v2);
              if(gad->Obj.fUpdate)
                gad->Obj.fUpdate((CObject*)gad);
              SceneChanged(G);
            }
            break;
          case cObjectMolecule:
            if(ObjectMoleculeGetAtomTxfVertex((ObjectMolecule*)obj,
                                              SettingGetGlobal_i(G,cSetting_state)-1,
                                              I->LastPicked.src.index,v1)) {
              /* scale properly given the current projection matrix */
              vScale = SceneGetExactScreenVertexScale(G,v1);
              if(stereo_via_adjacent_array(I->StereoMode)) {
                x = get_stereo_x(x,&I->LastX,I->Width, NULL);
              }
            
              switch(mode) {
              case cButModeMovFragZ:
              case cButModeMovObjZ:
              case cButModeMovViewZ:
              case cButModeMoveAtomZ:
                v2[0] = 0;
                v2[1] = 0;
                v2[2] = -(y-I->LastY)*vScale;
                break;
              default:
                v2[0] = (x-I->LastX)*vScale;
                v2[1] = (y-I->LastY)*vScale;
                v2[2] = 0;
                break;
              }

              v3[0] = 0.0F;
              v3[1] = 0.0F;
              v3[2] = 1.0F;
            
              /* transform into model coodinate space */
              MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); 
              MatrixInvTransformC44fAs33f3f(I->RotMatrix,v3,v3); 

              if(I->LastPicked.src.bond >= cPickableAtom) {
                if((mode!=cButModeMoveAtom)&&(mode!=cButModeMoveAtomZ)) {
                  EditorDrag(G,(ObjectMolecule*)obj,I->LastPicked.src.index,mode,
                             SettingGetGlobal_i(G,cSetting_state)-1,v1,v2,v3);
                } else {
                  int log_trans = (int)SettingGet(G,cSetting_log_conformations);
                  ObjectMoleculeMoveAtom((ObjectMolecule*)obj,
                                         SettingGetGlobal_i(G,cSetting_state)-1,
                                         I->LastPicked.src.index,v2,1,log_trans);
                  SceneInvalidate(G);
                }
              } else {
                int log_trans = (int)SettingGet(G,cSetting_log_conformations);
                ObjectMoleculeMoveAtomLabel((ObjectMolecule*)obj,
                                            SettingGetGlobal_i(G,cSetting_state)-1,
                                            I->LastPicked.src.index,v2,1,log_trans);
                SceneInvalidate(G);
              }
            }
            break;
          case cObjectSlice:
            {
              ObjectSlice *slice = (ObjectSlice*)obj;

              if(I->LastPickVertexFlag) {
              
                copy3f(I->LastPickVertex,v1);

                vScale = SceneGetExactScreenVertexScale(G,v1);

                if(stereo_via_adjacent_array(I->StereoMode)) {
                  x = get_stereo_x(x,&I->LastX,I->Width, NULL);
                }

                v2[0] = (x-I->LastX)*vScale;
                v2[1] = (y-I->LastY)*vScale;
                v2[2] = 0;
              
                v3[0] = 0.0F;
                v3[1] = 0.0F;
                v3[2] = 1.0F;
              
                /* transform into model coodinate space */
                MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); 
                MatrixInvTransformC44fAs33f3f(I->RotMatrix,v3,v3); 

                ObjectSliceDrag(slice,SceneGetState(G),mode,v1,v2,v3);
              }
            }
            break;
          case cObjectMeasurement:
            if(ObjectDistGetLabelTxfVertex((ObjectDist*)obj,
                                           SettingGetGlobal_i(G,cSetting_state)-1,
                                           I->LastPicked.src.index,v1)) {
              /* scale properly given the current projection matrix */
              vScale = SceneGetExactScreenVertexScale(G,v1);
              if(stereo_via_adjacent_array(I->StereoMode)) {
                x = get_stereo_x(x,&I->LastX,I->Width, NULL);
              }
            
              switch(mode) {
              case cButModeMovFragZ:
              case cButModeMovObjZ:
              case cButModeMovViewZ:
              case cButModeMoveAtomZ:
                v2[0] = 0;
                v2[1] = 0;
                v2[2] = -(y-I->LastY)*vScale;
                break;
              default:
                v2[0] = (x-I->LastX)*vScale;
                v2[1] = (y-I->LastY)*vScale;
                v2[2] = 0;
                break;
              }

              v3[0] = 0.0F;
              v3[1] = 0.0F;
              v3[2] = 1.0F;
            
              /* transform into model coodinate space */
              MatrixInvTransformC44fAs33f3f(I->RotMatrix,v2,v2); 
              MatrixInvTransformC44fAs33f3f(I->RotMatrix,v3,v3); 

              if(I->LastPicked.src.bond == cPickableLabel) {
                int log_trans = (int)SettingGet(G,cSetting_log_conformations);
                ObjectDistMoveLabel((ObjectDist*)obj,
                                    SettingGetGlobal_i(G,cSetting_state)-1,
                                    I->LastPicked.src.index,v2,1,log_trans);
                SceneInvalidate(G);
              }
            }
            break;
          default:
            break;
          }
      }
      I->LastX=x;
      I->LastY=y;
      break;
    case cButModeTransXY:

      SceneNoteMouseInteraction(G);

      vScale = SceneGetExactScreenVertexScale(G,I->Origin);
      if(stereo_via_adjacent_array(I->StereoMode)) {

        x = get_stereo_x(x,&I->LastX,I->Width, NULL);
      }

      v2[0] = (x-I->LastX)*vScale;
      v2[1] = (y-I->LastY)*vScale;
      v2[2] = 0.0F;
    
      moved_flag=false;
      if(I->LastX!=x)
        {
          I->Pos[0]+=v2[0];
          I->LastX=x;
          SceneInvalidate(G);
          moved_flag=true;
        }
      if(I->LastY!=y)
        {
          I->Pos[1]+=v2[1];
          I->LastY=y;
          SceneInvalidate(G);
          moved_flag=true;
        }
    
      EditorFavorOrigin(G,NULL);
      if(moved_flag&&(int)SettingGet(G,cSetting_roving_origin)) {
        SceneGetPos(G,v2); /* gets position of center of screen */
        SceneOriginSet(G,v2,true);
      }
      if(moved_flag&&(int)SettingGet(G,cSetting_roving_detail)) {    
        SceneRovingDirty(G);
      }
      break;
    case cButModeRotXYZ:
    case cButModeRotZ:
    case cButModeInvRotZ:
    case cButModeTransZ:
    case cButModeClipNF:
    case cButModeClipN:    
    case cButModeClipF:    
    
      SceneNoteMouseInteraction(G);
      
      eff_width = I->Width;
      if(stereo_via_adjacent_array(I->StereoMode)) {
        eff_width = I->Width/2;
        x = get_stereo_x(x,&I->LastX,I->Width, NULL);
      }
     
      if(SettingGet_b(G,NULL,NULL,cSetting_virtual_trackball)) {
        v1[0] = (float)(eff_width/2) - x;
        v1[1] = (float)(I->Height/2) - y;
      
        v2[0] = (float)(eff_width/2) - I->LastX;
        v2[1] = (float)(I->Height/2) - I->LastY;
      
      } else {
        v1[0] = (float)(I->LastX) - x;
        v1[1] = (float)(I->LastY) - y;
      
        v2[0] = 0;
        v2[1] = 0;
      }


      r1 = (float)sqrt1f(v1[0]*v1[0] + v1[1]*v1[1]);
      r2 = (float)sqrt1f(v2[0]*v2[0] + v2[1]*v2[1]);
     
      if(r1<scale) {
        v1[2] = (float)sqrt1f(scale*scale - r1*r1);
      } else {
        v1[2] = 0.0;
      }

      if(r2<scale) {
        v2[2] = (float)sqrt1f(scale*scale - r2*r2);
      } else {
        v2[2] = 0.0;
      }
      normalize23f(v1,n1);
      normalize23f(v2,n2);
      cross_product3f(n1,n2,cp);
      theta = (float)(SettingGet_f(G,NULL,NULL,cSetting_mouse_scale)*
                      2*180*asin(sqrt1f(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]))/cPI);

      dx = (v1[0]-v2[0]);
      dy = (v1[1]-v2[1]);
      dt = (float)(SettingGet_f(G,NULL,NULL,cSetting_mouse_limit)*sqrt1f(dx*dx+dy*dy)/scale);
    
      if(theta>dt)
        theta = dt;

      normalize23f(cp,axis);

      theta = theta/(1.0F+(float)fabs(axis[2]));

      v1[2]=0.0;
      v2[2]=0.0;
      normalize23f(v1,n1);
      normalize23f(v2,n2);
      cross_product3f(n1,n2,cp);
      omega = (float)(2*180*asin(sqrt1f(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]))/cPI);
      normalize23f(cp,axis2);     

      old_front = I->Front;
      old_back = I->Back;
      old_origin = -I->Pos[2];

      moved_flag=false;
      adjust_flag=false;
      switch(mode) {
      case cButModeRotXYZ:
        if(I->LastX!=x)
          {
            SceneRotate(G,theta,axis[0],axis[1],-axis[2]);
            I->LastX=x;
            adjust_flag=true;
          }
        if(I->LastY!=y)
          {
            SceneRotate(G,theta,axis[0],axis[1],-axis[2]);
            I->LastY=y;
            adjust_flag=true;
          }
        break;
      case cButModeRotZ:
        if(I->LastX!=x)
          {
            SceneRotate(G,omega,axis2[0],axis2[1],-axis2[2]);
            I->LastX=x;
            adjust_flag=true;
          }
        if(I->LastY!=y)
          {
            SceneRotate(G,omega,axis2[0],axis2[1],-axis2[2]);
            I->LastY=y;
            adjust_flag=true;        
          }
        break;
      case cButModeInvRotZ:
        if(I->LastX!=x) {
	  SceneRotate(G,(I->LastX-x)/2.0F,0.0F,0.0F,1.0F);
	  I->LastX=x;
	  adjust_flag=true;        
	}
        break;
      case cButModeTransZ:
        if(I->LastY!=y) {
          float factor;
          factor = 200/((I->FrontSafe+I->BackSafe)/2);
          if(factor>=0.0F) {
            factor = (((float)y)-I->LastY)/factor;
            if(!SettingGetGlobal_b(G,cSetting_legacy_mouse_zoom))
              factor = -factor;
            I->Pos[2]+=factor;
            I->Front-=factor;
            I->Back-=factor;
            I->FrontSafe = GetFrontSafe(I->Front,I->Back);
            I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
          }
          I->LastY=y;
          SceneInvalidate(G);
          adjust_flag=true;
        }
        break;
      case cButModeClipNF:
        if(I->LastX!=x)
          {
            I->Back-=(((float)x)-I->LastX)/10;
            if(I->Back<I->Front)
              I->Back=I->Front+cSliceMin;
            I->LastX=x;
            SceneInvalidate(G);
            moved_flag=true;
          }
        if(I->LastY!=y)
          {
            I->Front-=(((float)y)-I->LastY)/10;
            if(I->Front>I->Back)
              I->Front=I->Back+cSliceMin;
            I->LastY=y;
            SceneInvalidate(G);
            moved_flag=true;
          }
        if(moved_flag) {
          I->FrontSafe = GetFrontSafe(I->Front,I->Back);
          I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
        }
        break;
      case cButModeClipN:
        if(I->LastX!=x)
          {
            I->Front-=(((float)x)-I->LastX)/10;
            if(I->Front>I->Back)
              I->Front=I->Back+cSliceMin;
            I->FrontSafe = GetFrontSafe(I->Front,I->Back);
            I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
            I->LastX=x;
            SceneInvalidate(G);
            moved_flag=true;
          }
        if(I->LastY!=y)
          {
            I->Front-=(((float)y)-I->LastY)/10;
            if(I->Front>I->Back)
              I->Front=I->Back+cSliceMin;
            I->FrontSafe = GetFrontSafe(I->Front,I->Back);
            I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
            I->LastY=y;
            SceneInvalidate(G);
            moved_flag=true;
          }
        break;
      case cButModeClipF:
        if(I->LastX!=x)
          {
            I->Back-=(((float)x)-I->LastX)/10;
            if(I->Back<I->Front)
              I->Back=I->Front+cSliceMin;
            I->LastX=x;
            SceneInvalidate(G);
            moved_flag=true;
          }
        if(I->LastY!=y)
          {
            I->Back-=(((float)y)-I->LastY)/10;
            if(I->Back<I->Front)
              I->Back=I->Front+cSliceMin;
            I->LastY=y;
            SceneInvalidate(G);
            moved_flag=true;
          }
        break;
      }
      if(moved_flag)
        SceneDoRoving(G,old_front,old_back,old_origin,adjust_flag,false);
    }
  }
  if(I->PossibleSingleClick) {
    int max_single_click_drag = 4;
    int dx = abs(I->StartX-I->LastX);
    int dy = abs(I->StartY-I->LastY);
    if((dx>max_single_click_drag)||
       (dy>max_single_click_drag)) {
      I->PossibleSingleClick = false;
    }
  }
  return(1);
}

static int SceneDeferredClick(DeferredMouse *dm)
{
  if(!SceneClick(dm->block, dm->button, dm->x, dm->y, dm->mod, dm->when)) {
  }
  return 1;
}

static int SceneDeferredImage(DeferredImage *di)
{
  PyMOLGlobals *G=di->G;
  SceneMakeSizedImage(G,di->width, di->height,di->antialias);
  if(di->filename) {
    ScenePNG(G,di->filename, di->dpi, di->quiet, false, di->format);
    FreeP(di->filename);
  } else if(G->HaveGUI &&
            SettingGetGlobal_b(G,cSetting_auto_copy_images)) {
#ifdef _PYMOL_IP_EXTRAS
    if(IncentiveCopyToClipboard(G,di->quiet)) {
    }
#else
#ifdef PYMOL_EVAL
    PRINTFB(G,FB_Scene,FB_Warnings)
      " Warning: Clipboard image transfers disabled in Evaluation builds.\n"
      ENDFB(G);
#endif
#endif    
  }
  return 1;
}

int SceneDeferImage(PyMOLGlobals *G,int width, int height, 
                    char *filename, int antialias, float dpi, 
                    int format, int quiet)
{
  DeferredImage *di = Calloc(DeferredImage,1);
  if(di) {
    DeferredInit(G,&di->deferred);
    di->G = G;
    di->width = width;
    di->height = height;
    di->antialias = antialias;
    di->deferred.fn = (DeferredFn*)SceneDeferredImage;
    di->dpi = dpi;
    di->format = format;
    di->quiet = quiet;
    if(filename) {
      int stlen = strlen(filename);
      di->filename = Alloc(char,stlen+1);
      strcpy(di->filename, filename);
    }
  }
  OrthoDefer(G,&di->deferred);
  return 1;
}

int SceneDeferClick(Block *block, int button, int x, int y, int mod)
{
  PyMOLGlobals *G=block->G;
  DeferredMouse *dm = Calloc(DeferredMouse,1);
  if(dm) {
    DeferredInit(G,&dm->deferred);
    dm->block = block;
    dm->button = button;
    dm->x = x;
    dm->y = y;
    dm->mod = mod;
    dm->when = UtilGetSeconds(G);
    dm->deferred.fn = (DeferredFn*)SceneDeferredClick;
  }
  OrthoDefer(G,&dm->deferred);
  return 1;
}


static int SceneDeferClickWhen(Block *block, int button, int x, int y, double when,int mod)
{
  PyMOLGlobals *G=block->G;
  DeferredMouse *dm = Calloc(DeferredMouse,1);
  if(dm) {
    DeferredInit(G,&dm->deferred);
    dm->block = block;
    dm->button = button;
    dm->x = x;
    dm->y = y;
    dm->when = when;
    dm->mod = mod;
    dm->deferred.fn = (DeferredFn*)SceneDeferredClick;
  }
  OrthoDefer(G,&dm->deferred);
  return 1;
}

static int SceneDeferredDrag(DeferredMouse *dm)
{
  SceneDrag(dm->block, dm->x, dm->y, dm->mod, dm->when);
  return 1;
}

int SceneDeferDrag(Block *block,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  DeferredMouse *dm = Calloc(DeferredMouse,1);
  if(dm) {
    DeferredInit(G,&dm->deferred);
    dm->block = block;
    dm->x = x;
    dm->y = y;
    dm->mod = mod;
    dm->when = UtilGetSeconds(G);
    dm->deferred.fn = (DeferredFn*)SceneDeferredDrag;
  }
  OrthoDefer(G,&dm->deferred);
  return 1;
}

static int SceneDeferredRelease(DeferredMouse *dm)
{
  SceneRelease(dm->block, dm->button, dm->x, dm->y, dm->mod, dm->when);
  return 1;
}

int SceneDeferRelease(Block *block,int button,int x,int y,int mod) 
{
  PyMOLGlobals *G=block->G;
  DeferredMouse *dm = Calloc(DeferredMouse,1);
  if(dm) {
    DeferredInit(G,&dm->deferred);
    dm->block = block;
    dm->button = button;
    dm->x = x;
    dm->y = y;
    dm->mod = mod;
    dm->when = UtilGetSeconds(G);
    dm->deferred.fn = (DeferredFn*)SceneDeferredRelease;
  }
  OrthoDefer(G,&dm->deferred);
  return 1;
}
/*========================================================================*/
void SceneFree(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  if(I->ScrollBar)
    ScrollBarFree(I->ScrollBar);
  CGOFree(I->AlphaCGO);
  VLAFreeP(I->SceneVLA);
  VLAFreeP(I->SceneNameVLA);
  VLAFreeP(I->SlotVLA);
  OrthoFreeBlock(G,I->Block);
  ListFree(I->Obj,next,ObjRec);

  ScenePurgeImage(G);
  CGOFree(G->DebugCGO);
  FreeP(G->Scene);
}
/*========================================================================*/
void SceneResetMatrix(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  identity44f(I->RotMatrix);
  SceneUpdateInvMatrix(G);
}
/*========================================================================*/
void SceneSetDefaultView(PyMOLGlobals *G) {
  register CScene *I=G->Scene;

  identity44f(I->RotMatrix);
  identity44f(I->ModMatrix);
  identity44f(I->ProMatrix);
  SceneUpdateInvMatrix(G);
  
  I->ViewNormal[0]=0.0F;
  I->ViewNormal[1]=0.0F;
  I->ViewNormal[2]=1.0F;
  
  I->Pos[0] = 0.0F;
  I->Pos[1] = 0.0F;
  I->Pos[2] = -50.0F;

  I->Origin[0] = 0.0F;
  I->Origin[1] = 0.0F;
  I->Origin[2] = 0.0F;

  I->Front=40.0F;
  I->Back=100.0F;
  I->FrontSafe= GetFrontSafe(I->Front,I->Back);
  I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
  
  I->Scale = 1.0F;
  
}
int SceneReinitialize(PyMOLGlobals *G)
{
  int ok=true;
  SceneSetDefaultView(G);
  SceneCountFrames(G);
  SceneSetFrame(G,0,0);
  SceneInvalidate(G);
  return(ok);
}
/*========================================================================*/
int  SceneInit(PyMOLGlobals *G)
{
  register CScene *I=NULL;
  if( (I=(G->Scene=Calloc(CScene,1)))) { 

    /* all defaults to zero, so only initialize non-zero elements */

    G->DebugCGO = CGONew(G);

    ListInit(I->Obj);
    
    I->LoopFlag = false;

    I->TextColor[0]=0.2F;
    I->TextColor[1]=1.0F;
    I->TextColor[2]=0.2F;
    
    I->LastClickTime = UtilGetSeconds(G);
    I->LastPickVertexFlag = false;
    
    SceneSetDefaultView(G);
    
    I->HasMovie = false;
    I->Scale = 1.0;
    I->Block = OrthoNewBlock(G,NULL);
    I->Block->fClick   = SceneDeferClick;
    I->Block->fRelease = SceneDeferRelease;
    I->Block->fDrag    = SceneDeferDrag;
    I->Block->fDraw    = SceneDraw;
    I->Block->fReshape = SceneReshape;
    I->Block->active = true;
    
    OrthoAttach(G,I->Block,cOrthoScene);
    
    I->DirtyFlag = true;

    I->LastRender = UtilGetSeconds(G);
    I->LastFrameTime = UtilGetSeconds(G);

    I->LastSweepTime = UtilGetSeconds(G);

    I->LastPicked.context.object = NULL;
    I->LastStateBuilt = -1;
    I->CopyNextFlag=true;

    SceneRestartFrameTimer(G);
    SceneRestartPerfTimer(G);

    I->Width = 640; /* standard defaults */
    I->Height = 480;

    I->VertexScale = 0.01F;

    /* scene list */

    I->ScrollBar=ScrollBarNew(G,false);
    I->Pressed = -1;
    I->Over = -1;

    I->SceneNameVLA = VLAlloc(char,10);
    I->SceneVLA = VLAlloc(SceneElem,10);

    return 1;
  } else 
    return 0;
}
/*========================================================================*/
void SceneReshape(Block *block,int width,int height)
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  int y = height;

  if(I->Block->margin.right) {
     width -= I->Block->margin.right;
     if(width<1)
        width=1;
  }
  
  if(I->Block->margin.top) {
    height -= I->Block->margin.top;
    y = height;
  }

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _PYMOL_OSX
#ifndef _MACPYMOL_XCODE
  /* workaround for broken pixel handling under OSX 
     (Who's fault: Me? Apple? NVidia?) */
  width = 8*(width/8);
  /* it's 2007, four years later, do we still need this? */
#endif
#endif
/* END PROPRIETARY CODE SEGMENT */

  I->Width = width;

  I->Height = height;

  I->Block->rect.top = I->Height;
  I->Block->rect.left = 0;
  I->Block->rect.bottom = 0;
  I->Block->rect.right = I->Width;

  if(I->Block->margin.bottom) {
     height-=I->Block->margin.bottom;
     if(height<1)
        height=1;
     I->Height=height;
     I->Block->rect.bottom=I->Block->rect.top - I->Height;
  }
  SceneDirty(G);

  if(I->CopyType&&(!I->CopyForced)) {
    SceneInvalidateCopy(G,false);
  }
  /*MovieClearImages(G);*/
  MovieSetSize(G,I->Width,I->Height);
  SceneInvalidateStencil(G);
}

/*========================================================================*/
void SceneDone(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  if(I->Block)
     OrthoFreeBlock(G,I->Block);
}
/*========================================================================*/
void SceneResetNormal(PyMOLGlobals *G,int lines)
{
  register CScene *I=G->Scene;
  if(G->HaveGUI && G->ValidContext) {
    if(lines)
      glNormal3fv(I->LinesNormal);
    else
      glNormal3fv(I->ViewNormal);
  }
}

/*========================================================================*/
static void SceneApplyImageGamma(PyMOLGlobals *G,unsigned int *buffer, int width, int height)
{
  unsigned int test;
  unsigned char *testPtr;
  int big_endian;
  float gamma = SettingGet(G,cSetting_gamma);

  if(gamma>R_SMALL4)
    gamma=1.0F/gamma;
  else
    gamma=1.0F;

  test = 0xFF000000;
  testPtr = (unsigned char*)&test;
  big_endian = (*testPtr)&&1;

  if(buffer&&height&&width) {
    register float _inv3 = 1/(255*3.0F);
    register float _1 = 1/3.0F;
    register unsigned char *p;
    register int x,y;
    register float c1,c2,c3,inp,sig;
    register unsigned int i1,i2,i3;
    p = (unsigned char*) buffer;
    for(y=0;y<height;y++) {
      for(x=0;x<width;x++) {
        c1 = p[0];
        c2 = p[1];
        c3 = p[2];
        inp = (c1+c2+c3) * _inv3;
        if(inp < R_SMALL4) 
          sig = _1;
        else
          sig = (float)(pow(inp,gamma) / inp);
        i1 = (unsigned int)(sig * c1);
        i2 = (unsigned int)(sig * c2);
        i3 = (unsigned int)(sig * c3);
        if(i1>255) i1 = 255;
        if(i2>255) i2 = 255;
        if(i3>255) i3 = 255;
        p[0] = i1;
        p[1] = i2;
        p[2] = i3;
        p+=4;
      }
    }
  }
}
/*========================================================================*/

static double accumTiming = 0.0; 

void SceneDoRay(PyMOLGlobals *G,int width,int height,int mode,
                char **headerVLA,char **charVLA,
                float angle,float shift,int quiet,
                G3dPrimitive **g3d,int show_timing,int antialias)
{
  SceneRay(G, width, height, mode,
           headerVLA, charVLA,
           angle,shift, quiet,
           g3d, show_timing, antialias);
}

static int SceneDeferredRay(DeferredRay *dr)
{
  PyMOLGlobals *G=dr->G;
  SceneRay(G, dr->ray_width, dr->ray_height, dr->mode,
           NULL, NULL, dr->angle, dr->shift, dr->quiet, 
           NULL, dr->show_timing, dr->antialias);
  if((dr->mode==0) && 
     G->HaveGUI &&
     SettingGetGlobal_b(G,cSetting_auto_copy_images)) {
#ifdef _PYMOL_IP_EXTRAS
    IncentiveCopyToClipboard(G,dr->quiet);
#else
#ifdef PYMOL_EVAL
    PRINTFB(G,FB_Scene,FB_Warnings)
      " Warning: Clipboard image transfers disabled in Evaluation Builds.\n"
      ENDFB(G);
#endif
#endif    
  }
  return 1;
}


int SceneDeferRay(PyMOLGlobals *G,
                   int ray_width,
                   int ray_height,
                   int mode,
                   float angle,
                   float shift,
                   int quiet,
                   int show_timing,
                   int antialias)
{
  DeferredRay *dr = Calloc(DeferredRay,1);
  if(dr) {
    DeferredInit(G,&dr->deferred);
    dr->G = G;
    dr->ray_width = ray_width;
    dr->ray_height = ray_height;
    dr->mode = mode;
    dr->angle = angle;
    dr->shift = shift;
    dr->quiet = quiet;
    dr->show_timing = show_timing;
    dr->antialias = antialias;
    dr->deferred.fn = (DeferredFn*)SceneDeferredRay;
  }
  OrthoDefer(G,&dr->deferred);
  return 1;
}

static void SceneUpdateAnimation(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  int rockFlag = false;
  int dirty = false;
  int movie_rock = SettingGetGlobal_b(G,cSetting_movie_rock);
  
  if(movie_rock<0) movie_rock = ControlRocking(G);

  if(MoviePlaying(G) && movie_rock) {

    if(MovieGetRealtime(G) && 
       ! SettingGetGlobal_b(G,cSetting_movie_animate_by_frame)) {
      I->SweepTime+=I->RenderTime;
      rockFlag = true;
      dirty = true; /* force a subsequent update */
    } else {
      float fps = SceneGetFPS(G); /* guaranteed to be >= 0.0F */
      if(fps>0.0F) {
        int rock_frame = SceneGetFrame(G);
        if(rock_frame!=I->RockFrame) {
          I->RockFrame = rock_frame;
          rockFlag = true;
          I->SweepTime += 1.0/fps;
        }
      } else {
        I->SweepTime += I->RenderTime;
        rockFlag = true;
      }
    }
  } else
    dirty = true;

  if(I->cur_ani_elem < I->n_ani_elem ) { /* play motion animation */
    double now;

    int cur = I->cur_ani_elem;

    if(I->AnimationStartFlag) { 
      /* allow animation timing to lag since it may take a few seconds
         to get here given geometry updates, etc. */
                                   
      I->AnimationLagTime = UtilGetSeconds(G) - I->AnimationStartTime;
      I->AnimationStartFlag = false;
    }

    if( (!MoviePlaying(G)) || 
        ((MovieGetRealtime(G) &&
          ! SettingGetGlobal_b(G,cSetting_movie_animate_by_frame)))) {
      now = UtilGetSeconds(G) - I->AnimationLagTime;
    } else {
      float fps = SceneGetFPS(G); /* guaranteed to be >= 0.0F */
      int frame = SceneGetFrame(G);
      int n_frame = 0;

      cur = 0; /* allow backwards interpolation */
      if(frame >= I->AnimationStartFrame) {
        n_frame = frame - I->AnimationStartFrame;
      } else {
        n_frame = frame + (I->NFrame - I->AnimationStartFrame);
      }
      now = I->AnimationStartTime + n_frame / fps;
    }

    while(I->ani_elem[cur].timing<now) {
      cur++;
      if(cur >= I->n_ani_elem) {
        cur = I->n_ani_elem;
        break;
      }
    }
    I->cur_ani_elem = cur;
    SceneFromViewElem(G,I->ani_elem+cur,dirty);
  }
  if(rockFlag && (I->SweepTime!=0.0)) {
    SceneUpdateCameraRock(G,dirty);
  }
}

static int SceneGetDrawFlag(GridInfo *grid, int *slot_vla, int slot)
{
  int draw_flag = false;
  if(grid && grid->active) {
    switch(grid->mode) {
    case 1: /* assigned grid slots (usually by group) */
      {
        if(((slot<0) && grid->slot) || 
           ((slot==0) && (grid->slot==0)) ||
           (slot_vla && (slot_vla[slot] == grid->slot))) {
          draw_flag = true;
        }
      }
      break;
    case 2: /* each state in a separate slot */
      draw_flag = true;
      break;
    }
  } else {
    draw_flag = true;
  }
  return draw_flag;
}

void SceneRay(PyMOLGlobals *G,
              int ray_width,int ray_height,int mode,
              char **headerVLA_ptr,
              char **charVLA_ptr,float angle,
              float shift,int quiet,
              G3dPrimitive **g3d,int show_timing,
              int antialias)
{

  register CScene *I=G->Scene;
  ObjRec *rec=NULL;
  CRay *ray =NULL;
  float height,width;
  float aspRat;
  float rayView[16];
  int curState;
  double timing;
  char *charVLA = NULL;
  char *headerVLA = NULL;
  float fov;
  int stereo_hand = 0;
  GridInfo grid;
  int grid_mode = SettingGetGlobal_i(G,cSetting_grid_mode);
  ImageType *stereo_image = NULL;
  OrthoLineType prefix = "";
 int ortho = SettingGetGlobal_i(G,cSetting_ray_orthoscopic);
  
  if(SettingGetGlobal_b(G,cSetting_defer_builds_mode) == 5) 
    SceneUpdate(G,true);
  
  if(ortho<0) ortho = SettingGetGlobal_b(G,cSetting_ortho);

  UtilZeroMem(&grid,sizeof(GridInfo));

  if(mode!=0) grid_mode = 0; /* only allow grid mode with PyMOL renderer */

  SceneUpdateAnimation(G);
  if(mode==0) 
    SceneInvalidateCopy(G,true);

  if(antialias<0) {
    antialias = (int)SettingGet(G,cSetting_antialias);
  }
  if(ray_width<0) ray_width = 0;
  if(ray_height<0) ray_height = 0;
  if((!ray_width)||(!ray_height)) {
    if(ray_width&&(!ray_height)) {
      ray_height = (ray_width*I->Height)/I->Width;
    } else if(ray_height&&(!ray_width)) {
      ray_width = (ray_height*I->Width)/I->Height;
    } else {
      ray_width=I->Width;
      ray_height=I->Height;
    }
  }

  fov=SettingGet(G,cSetting_field_of_view);

  if(SettingGet(G,cSetting_all_states)) {
    curState=-1;
  } else {
    curState=SettingGetGlobal_i(G,cSetting_state)-1;
  }

  timing = UtilGetSeconds(G); /* start timing the process */

  SceneUpdate(G,false);
    
  switch(I->StereoMode) {
  case cStereo_quadbuffer:
    stereo_hand=2;
    break;
  case cStereo_crosseye:
  case cStereo_walleye:
    ray_width = ray_width/2;
    stereo_hand=2;
    break;
  case cStereo_geowall:
  case cStereo_sidebyside:
    stereo_hand=2;
    break;
  case cStereo_stencil_by_row:
  case cStereo_stencil_by_column:
  case cStereo_stencil_checkerboard:
  case cStereo_stencil_custom:
  case cStereo_anaglyph:
    stereo_hand=2;
    break;
  }

  aspRat = ((float) ray_width) / ((float) ray_height);

  if(grid_mode) {
    int grid_size = SceneGetGridSize(G,grid_mode);
    GridUpdate(&grid, aspRat, grid_mode, grid_size);
    if(grid.active) 
      aspRat *= grid.asp_adjust;
  }

  while(1) {
    int slot;
    int tot_width = ray_width;
    int tot_height = ray_height;
    int ray_x = 0,ray_y = 0;

    if(grid.active)
      GridGetRayViewport(&grid,ray_width,ray_height);
    
    for(slot=0;slot<=grid.last_slot;slot++) {
      
      if(grid.active) { 
        GridSetRayViewport(&grid,slot,&ray_x,&ray_y,&ray_width,&ray_height);
        OrthoBusySlow(G,slot,grid.last_slot);
      }
    
      /* start afresh, looking in the negative Z direction (0,0,-1) from (0,0,0) */
      identity44f(rayView);
      
      ray = RayNew(G,antialias);
      if(!ray) break;
      
      if(stereo_hand) {
        /* stereo */
        
        float stAng,stShift;
        
        stAng = SettingGet(G,cSetting_stereo_angle);
        stShift = SettingGet(G,cSetting_stereo_shift);
        
        /* right hand */
        
        stShift = (float)(stShift*fabs(I->Pos[2])/100.0);
        stAng = (float)(stAng*atan(stShift/fabs(I->Pos[2]))*90.0/cPI);
        
        if(stereo_hand==2) { /* left hand */
          stAng=-stAng;
          stShift=-stShift;
        }
        
        angle = stAng;
        
        {
          float temp[16];
          identity44f(temp);
          MatrixRotateC44f(temp,(float)(-PI*stAng/180),0.0F,1.0F,0.0F); /* y-axis rotation */
          MatrixMultiplyC44f(temp,rayView);
        }
        
        /* move the camera to the location we are looking at */
        MatrixTranslateC44f(rayView,I->Pos[0],I->Pos[1],I->Pos[2]);
        MatrixTranslateC44f(rayView,stShift,0.0,0.0);
        
        MatrixMultiplyC44f(I->RotMatrix,rayView);
        
      } else { /* not stereo mode */
        
        /* move the camera to the location we are looking at */
        MatrixTranslateC44f(rayView,I->Pos[0],I->Pos[1],I->Pos[2]);
        
        if(shift) {
          MatrixTranslateC44f(rayView,shift,0.0F,0.0F);
        }
        /* move the camera so that we can see the origin 
         * NOTE, vector is given in the coordinates of the world's motion
         * relative to the camera */
        
        
        /* 4. rotate about the origin (the the center of rotation) */
        
        if(angle) {
          float temp[16];
          identity44f(temp);
          MatrixRotateC44f(temp,(float)(-PI*angle/180),0.0F,1.0F,0.0F);
          MatrixMultiplyC44f(I->RotMatrix,temp);
          MatrixMultiplyC44f(temp,rayView);
        } else {
          MatrixMultiplyC44f(I->RotMatrix,rayView);
        }
      }
      
      /* 5. move the origin to the center of rotation */
      MatrixTranslateC44f(rayView,-I->Origin[0],-I->Origin[1],-I->Origin[2]);
      
      if(Feedback(G,FB_Scene,FB_Debugging)) {
        fprintf(stderr,"SceneRay: %8.3f %8.3f %8.3f\n",
                I->Pos[0],I->Pos[1],I->Pos[2]);
        fprintf(stderr,"SceneRay: %8.3f %8.3f %8.3f\n",
                I->Origin[0],I->Origin[1],I->Origin[2]);
        fprintf(stderr,"SceneRay: %8.3f %8.3f %8.3f\n",
                I->RotMatrix[0],I->RotMatrix[1],I->RotMatrix[2]);
      }
      /* define the viewing volume */
      
      height  = (float)(fabs(I->Pos[2])*tan((fov/2.0)*cPI/180.0));     
      width = height*aspRat;
      
      
      OrthoBusyFast(G,0,20);
      
      {
        float pixel_scale_value = SettingGetGlobal_f(G,cSetting_ray_pixel_scale);
        float fov=SettingGet(G,cSetting_field_of_view);
        
        if(pixel_scale_value<0) pixel_scale_value = 1.0F;
        

        pixel_scale_value *= ((float)tot_height)/I->Height;
        
        if(ortho) {
          const float _1 = 1.0F;
          RayPrepare(ray,-width,width,-height,height,
                     I->FrontSafe,I->BackSafe,
                     fov, I->Pos, 
                     rayView,I->RotMatrix,aspRat,
                     ray_width, ray_height, 
                     pixel_scale_value, ortho,
                     _1, _1, /* gcc 3.2.3 blows chunks if these are 1.0F */
                     ((float)ray_height)/I->Height);
          
        } else {        
          float back_ratio;
          float back_height;
          float back_width;
          float pos;
          float fov = SettingGet(G,cSetting_field_of_view);
          pos = I->Pos[2];

          if((-pos)<I->FrontSafe) {
            pos = -I->FrontSafe;
          }

          back_ratio = -I->Back/pos;
          back_height = back_ratio*height;
          back_width = aspRat * back_height;
          RayPrepare(ray,
                     -back_width, back_width, 
                     -back_height, back_height,
                     I->FrontSafe,I->BackSafe,
                     fov, I->Pos,
                     rayView,I->RotMatrix,aspRat,
                     ray_width, ray_height, 
                     pixel_scale_value, ortho,
                     height/back_height,
                     I->FrontSafe/I->BackSafe,
                     ((float)ray_height)/I->Height);
        }
      }
      {
        int *slot_vla = I->SlotVLA;
        int state = SceneGetState(G);
        RenderInfo info;
        UtilZeroMem(&info,sizeof(RenderInfo));
        info.ray = ray;
        info.ortho = ortho;
        info.vertex_scale = SceneGetScreenVertexScale(G,NULL);

        if(SettingGetGlobal_b(G,cSetting_dynamic_width)) {
          info.dynamic_width = true;
          info.dynamic_width_factor = SettingGetGlobal_f(G,cSetting_dynamic_width_factor);
          info.dynamic_width_min = SettingGetGlobal_f(G,cSetting_dynamic_width_min);
          info.dynamic_width_max = SettingGetGlobal_f(G,cSetting_dynamic_width_max);
        }

        while(ListIterate(I->Obj,rec,next)) {
          if(rec->obj->fRender) {
            if(SceneGetDrawFlag(&grid, slot_vla, rec->obj->grid_slot)) {
              int obj_color = rec->obj->Color;
              float color[3];
              int icx;
              ColorGetEncoded(G,obj_color,color);
              RaySetContext(ray,rec->obj->Context);
              ray->fColor3fv(ray,color);
              
              if(SettingGetIfDefined_i(G,rec->obj->Setting,cSetting_ray_interior_color,&icx)) {
                float icolor[3];
                if(icx!=-1) {
                  if(icx==cColorObject) {
                    ray->fInteriorColor3fv(ray,color,false);              
                  } else {
                    ColorGetEncoded(G,icx,icolor);
                    ray->fInteriorColor3fv(ray,icolor,false);              
                  }
                } else {
                  ray->fInteriorColor3fv(ray,color,true);
                }
              } else {
                ray->fInteriorColor3fv(ray,color,true);
              }
              if((!grid.active)||(grid.mode!=2)) {
                info.state = ObjectGetCurrentState(rec->obj,false);
                rec->obj->fRender(rec->obj,&info);
              } else if(grid.slot) {
                if ( (info.state = state + grid.slot - 1) >= 0 )
                  rec->obj->fRender(rec->obj,&info);
              }
            }
          }
        }
      }

      OrthoBusyFast(G,1,20);
    
      if(mode!=2) { /* don't show pixel count for tests */
        if(!quiet) {
          PRINTFB(G,FB_Ray,FB_Blather)
            " Ray: tracing %dx%d = %d rays against %d primitives.\n",ray_width,ray_height,
            ray_width*ray_height,RayGetNPrimitives(ray)
            ENDFB(G);
        }
      }
      switch(mode) {
      case 0: /* mode 0 is built-in */
        {
          unsigned int buffer_size = 4*ray_width*ray_height;
          unsigned int *buffer=(GLvoid*)Alloc(char,buffer_size);
          unsigned int background;
          ErrChkPtr(G,buffer);
          
          RayRender(ray,buffer,timing,angle,antialias,&background);
          
          /*    RayRenderColorTable(ray,ray_width,ray_height,buffer);*/
          
          if(!grid.active) {
            I->Image=Calloc(ImageType,1);
            I->Image->data = (unsigned char*)buffer;
            I->Image->size = buffer_size;
            I->Image->width = ray_width;
            I->Image->height = ray_height;
          } else {
            if(!I->Image) { /* alloc on first pass */
              I->Image=Calloc(ImageType,1);
              if(I->Image) {
                unsigned int tot_size = 4*tot_width*tot_height;
                I->Image->data = (GLvoid*)Alloc(char,tot_size);
                I->Image->size = tot_size;
                I->Image->width = tot_width;
                I->Image->height = tot_height;
                { /* fill with background color */
                  unsigned int i;
                  unsigned int *ptr = (unsigned int*)I->Image->data;
                  for(i=0;i<tot_size;i+=4) {
                    *(ptr++) = background;
                  }
                }
              }
            }
            /* merge in the latest rendering */
            if(I->Image && I->Image->data) {
              int i,j;
              unsigned int *src = buffer;
              unsigned int *dst = (unsigned int*)I->Image->data;
              
              dst += (ray_x + ray_y*tot_width);
              
              for(i=0;i<ray_height;i++) {
                for(j=0;j<ray_width;j++) {
                  if(*src != background) 
                    *(dst) = *(src);
                  dst++;
                  src++;
                }
                dst += (tot_width - ray_width);
              }
            }
            FreeP(buffer);
          }
          I->DirtyFlag=false;
          I->CopyType = true;
          I->CopyForced = true;
          I->MovieOwnsImageFlag = false;
        }
        break;
      
      case 1: /* mode 1 is povray */
        charVLA=VLACalloc(char,100000); 
        headerVLA=VLACalloc(char,2000);
        RayRenderPOV(ray,ray_width,ray_height,&headerVLA,&charVLA,
                     I->FrontSafe,I->BackSafe,fov,angle,antialias);
        if(!(charVLA_ptr&&headerVLA_ptr)) { /* immediate mode */
          strcpy(prefix,SettingGet_s(G,NULL,NULL,cSetting_batch_prefix));
#ifndef _PYMOL_NOPY
          if(PPovrayRender(G,headerVLA,charVLA,prefix,ray_width,
                           ray_height,antialias)) {
            strcat(prefix,".png");
            SceneLoadPNG(G,prefix,false,0,false);
            I->DirtyFlag=false;
          }
#endif
          VLAFreeP(charVLA);
          VLAFreeP(headerVLA);
        } else { /* get_povray mode */
          *charVLA_ptr=charVLA;
          *headerVLA_ptr=headerVLA;
        }
        break;
      case 2: /* mode 2 is for testing of geometries */
        RayRenderTest(ray,ray_width,ray_height,I->FrontSafe,I->BackSafe,fov);
        break;
      case 3: /* mode 3 is for Jmol */
        {
          G3dPrimitive *jp = RayRenderG3d(ray,ray_width,ray_height,I->FrontSafe,I->BackSafe,fov,quiet);
          if(0) {
            int cnt = VLAGetSize(jp);
            int a;
            for(a=0;a<cnt;a++) {
              switch(jp[a].op) {
              case 1:
                printf("g3d.fillSphereCentered(gray,%d,%d,%d,%d);\n",jp[a].r,jp[a].x1,jp[a].y1,jp[a].z1);
                break;
              case 2:
                printf("triangle(%d,%d,%d,%d,%d,%d,%d,%d,%d);\n",
                       jp[a].x1,jp[a].y1,jp[a].z1,
                       jp[a].x2,jp[a].y2,jp[a].z2,
                       jp[a].x3,jp[a].y3,jp[a].z3
                       );
                break;
              case 3:
                printf("g3d.fillCylinder(gray,gray,(byte)3,%d,%d,%d,%d,%d,%d,%d);\n",
                       jp[a].r,
                       jp[a].x1,jp[a].y1,jp[a].z1,
                       jp[a].x2,jp[a].y2,jp[a].z2
                       );          
                break;
              }
            }
          }
          if(g3d) {
            *g3d = jp;
          } else {
            VLAFreeP(jp);
          }
        }
        break;
      case 4: /* VRML2 */
        {
          char *vla = VLACalloc(char,100000);
          RayRenderVRML2(ray,ray_width,ray_height,&vla,
                         I->FrontSafe,I->BackSafe,fov,angle,I->Pos[2]);
          *charVLA_ptr=vla;
        }
        break;
      case 5: /* mode 5 is OBJ MTL */
        {
          char *objVLA=VLACalloc(char,100000); 
          char *mtlVLA=VLACalloc(char,1000);
          RayRenderObjMtl(ray,ray_width,ray_height,&objVLA,&mtlVLA,
                          I->FrontSafe,I->BackSafe,fov,angle,I->Pos[2]);
          *headerVLA_ptr=objVLA;
          *charVLA_ptr=mtlVLA;
        }
        break;
      case 6: /* VRML1 -- more compatible with tools like blender */
        {
          char *vla = VLACalloc(char,100000);
          RayRenderVRML1(ray,ray_width,ray_height,&vla,
                         I->FrontSafe,I->BackSafe,fov,angle,I->Pos[2]);
          *charVLA_ptr=vla;
        }
        break;
      case cSceneRay_MODE_IDTF:
        {
          *headerVLA_ptr = VLACalloc(char,10000);
          *charVLA_ptr = VLACalloc(char,10000);
          RayRenderIDTF(ray,headerVLA_ptr,charVLA_ptr);
        }
        break;
    
      }
      RayFree(ray);
    }
    if(grid.active)
      GridSetRayViewport(&grid,-1,&ray_x,&ray_y,&ray_width,&ray_height);
    
    if((mode==0)&&I->Image&&I->Image->data) {
      SceneApplyImageGamma(G,(unsigned int*)I->Image->data,I->Image->width,I->Image->height);
    }

    stereo_hand--;
    if((I->StereoMode==0)||(stereo_hand<=0))
      break;
    else {
      stereo_image = I->Image;
      I->Image = NULL;
    }
  }

  if(stereo_image) {
    if(I->Image) {
      switch(I->StereoMode) {
      case cStereo_quadbuffer:
      case cStereo_geowall:
        /* merge the two images into one pointer */
        I->Image->data = Realloc(I->Image->data,unsigned char,I->Image->size*2);
        UtilCopyMem(I->Image->data+I->Image->size,
                    stereo_image->data,
                    I->Image->size);
        I->Image->stereo = true;
        break;
      case cStereo_crosseye:
      case cStereo_walleye:
        {
          /* merge the two images into one */
          
          unsigned char *merged_image = Alloc(unsigned char,I->Image->size*2);
          unsigned int *q=(unsigned int*)merged_image;
          unsigned int *l;
          unsigned int *r;
          register int height,width;
          register int a,b;
          
          if(I->StereoMode==2) {
            l=(unsigned int*)stereo_image->data;
            r=(unsigned int*)I->Image->data;
          } else {
            r=(unsigned int*)stereo_image->data;
            l=(unsigned int*)I->Image->data;
          }
          height = I->Image->height;
          width = I->Image->width;
          
          for(a=0;a<height;a++) {
            for(b=0;b<width;b++)
              *(q++) = *(l++);
            for(b=0;b<width;b++)
              *(q++) = *(r++);
          }
          FreeP(I->Image->data);
          I->Image->data = merged_image;
          I->Image->width*=2;
          I->Image->size*=2;
        }
        break;
      case cStereo_stencil_by_row:
      case cStereo_stencil_by_column:
      case cStereo_stencil_checkerboard:
        {
          /* merge the two images into one */
          
          unsigned char *merged_image = Alloc(unsigned char,I->Image->size);
          unsigned int *q=(unsigned int*)merged_image;
          unsigned int *l;
          unsigned int *r;
          register int height,width;
          register int a,b;
          
          l=(unsigned int*)stereo_image->data;
          r=(unsigned int*)I->Image->data;

          height = I->Image->height;
          width = I->Image->width;
          
          for(a=0;a<height;a++) {
            for(b=0;b<width;b++) {
              switch(I->StereoMode) {
              case cStereo_stencil_by_row:
                if(a&0x1) {
                  *(q++) = *(l++); r++;
                } else {
                  *(q++) = *(r++); l++;
                }
                break;
              case cStereo_stencil_by_column:
                if(b&0x1) {
                  *(q++) = *(l++); r++;
                } else {
                  *(q++) = *(r++); l++;
                }
                break;
              case cStereo_stencil_checkerboard:
                if((a+b)&0x1) {
                  *(q++) = *(l++); r++;
                } else {
                  *(q++) = *(r++); l++;
                }
                break;
              }
            }
          }
          FreeP(I->Image->data);
          I->Image->data = merged_image;
        }
        break;
      }
    }
    FreeP(stereo_image->data);
    FreeP(stereo_image);
  }
  timing = UtilGetSeconds(G)-timing;
  if(mode!=2) { /* don't show timings for tests */
    accumTiming += timing; 
    
    if(show_timing && !quiet) {
      if(!G->Interrupt) {
        PRINTFB(G,FB_Ray,FB_Details)
          " Ray: render time: %4.2f sec. = %3.1f frames/hour (%4.2f sec. accum.).\n", 
          timing,3600/timing, 
          accumTiming 
          ENDFB(G);
      } else {
         PRINTFB(G,FB_Ray,FB_Details)
          " Ray: render aborted.\n"
          ENDFB(G);
      }
    }
  }
  
  if(mode!=3) {
    OrthoDirty(G);
  }

}
/*========================================================================*/
static void SceneCopy(PyMOLGlobals *G,GLenum buffer,int force,int entire_window)
{
  register CScene *I=G->Scene;
  unsigned int buffer_size;
  
  if(force || (!(I->StereoMode||
                 SettingGet(G,cSetting_stereo_double_pump_mono)||
                 I->ButtonsShown))) {
    /* no copies while in stereo mode */
    if(force || ((!I->DirtyFlag)&&(!I->CopyType))) { 
      int x,y,w,h;
      if(entire_window) {
        x = 0;
        y = 0;
        h = OrthoGetHeight(G);
        w = OrthoGetWidth(G);
      } else {
        x = I->Block->rect.left;
        y = I->Block->rect.bottom;
        w = I->Width;
        h = I->Height;
      }
      ScenePurgeImage(G);
      buffer_size = 4*w*h;
      if(buffer_size) {
        I->Image = Calloc(ImageType,1);
        I->Image->data=(GLvoid*)Alloc(char,buffer_size);
        I->Image->size = buffer_size;
        I->Image->width = w;
        I->Image->height = h;
        if(G->HaveGUI && G->ValidContext) {
          glReadBuffer(buffer);
          PyMOLReadPixels(x,y,w,h,GL_RGBA,GL_UNSIGNED_BYTE,I->Image->data);
        }
      }
      I->CopyType = true;
      I->Image->needs_alpha_reset = true;
      I->CopyForced = force;
    }
  }
}

/*========================================================================*/
int SceneRovingCheckDirty(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  return(I->RovingDirtyFlag);
}

struct _CObjectUpdateThreadInfo {
  CObject *obj;
};

void SceneObjectUpdateThread(CObjectUpdateThreadInfo *T)
{
  if(T->obj && T->obj->fUpdate) {
    T->obj->fUpdate(T->obj);
  }
}

#ifndef _PYMOL_NOPY
static void SceneObjectUpdateSpawn(PyMOLGlobals *G,CObjectUpdateThreadInfo *Thread,int n_thread,int n_total)
{
  if(n_total==1) {
    SceneObjectUpdateThread(Thread);
  } else if(n_total){
    int blocked;
    PyObject *info_list;
    int a,n=0;
    blocked = PAutoBlock(G);
    
    PRINTFB(G,FB_Scene,FB_Blather)
      " Scene: updating objects with %d threads...\n",n_thread
      ENDFB(G);
    info_list = PyList_New(n_total);
    for(a=0;a<n_total;a++) {
      PyList_SetItem(info_list,a,PyCObject_FromVoidPtr(Thread+a,NULL));
      n++;
    }
    PXDecRef(PyObject_CallMethod(G->P_inst->cmd,"_object_update_spawn","Oi",info_list,n_thread));
    Py_DECREF(info_list);
    PAutoUnblock(G,blocked);
  }
}
#endif

/*========================================================================*/
void SceneUpdate(PyMOLGlobals *G, int force)
{
  register CScene *I=G->Scene;
  ObjRec *rec=NULL;
  int cur_state = SettingGetGlobal_i(G,cSetting_state) - 1;
  int defer_builds_mode = SettingGetGlobal_b(G,cSetting_defer_builds_mode);

  PRINTFD(G,FB_Scene)
    " SceneUpdate: entered.\n"
    ENDFD;

  OrthoBusyPrime(G);
  EditorUpdate(G);
  if(defer_builds_mode == 0) {
    if(SettingGetGlobal_i(G,cSetting_draw_mode)==-2) {
      defer_builds_mode=1;
    }
  }
  if(force || I->ChangedFlag || ((cur_state != I->LastStateBuilt) && 
                                 (defer_builds_mode>0))) {
    SceneCountFrames(G);

    if(force || (defer_builds_mode!=5)) { /* mode 5 == immediate mode */

      PyMOL_SetBusy(G->PyMOL,true); /*  race condition -- may need to be fixed */
      {
#ifndef _PYMOL_NOPY
        int n_thread  = SettingGetGlobal_i(G,cSetting_max_threads);
        int multithread = SettingGetGlobal_i(G,cSetting_async_builds);
        if(multithread && (n_thread>1)) {
          int min_start = -1;
          int max_stop = -1;
          int n_frame = SceneGetNFrame(G,NULL);
          int n_obj = 0;
          while(ListIterate(I->Obj,rec,next)) {
            int start = 0;
            int stop = n_frame;
            n_obj++;
            if(rec->obj->fGetNFrame) {
              stop = rec->obj->fGetNFrame(rec->obj);
            } 
            ObjectAdjustStateRebuildRange(rec->obj,&start,&stop);
            if(min_start<0) {
              min_start = start;
              max_stop = stop;
            } else {
              if(min_start>start)
                min_start = start;
              if(max_stop<stop)
                max_stop = stop;
            }
          }
        
          n_frame = max_stop - min_start;

          if( n_frame > n_thread ) {
            n_thread = 1;
            /* prevent n_thread * n_thread -- only multithread within
               individual object states (typically more balanced) */
          } else if( n_frame > 1 ) {
            n_thread = n_thread / n_frame;
          }
        
          if(n_thread < 1)
            n_thread = 1;
        }

        if(multithread && (n_thread>1)) {
          /* multi-threaded geometry update */
          int cnt = 0;

          rec = NULL;
          while(ListIterate(I->Obj,rec,next))
            cnt++;
        
          if(cnt) {
            CObjectUpdateThreadInfo *thread_info = Alloc(CObjectUpdateThreadInfo, cnt);
            if(thread_info) {
              cnt = 0;
              while(ListIterate(I->Obj,rec,next))
                thread_info[cnt++].obj = rec->obj;
              SceneObjectUpdateSpawn(G,thread_info,n_thread,cnt);
              FreeP(thread_info);
            }
          }
        } else 
#endif
          {
            /* single-threaded update */
            rec = NULL;
            while(ListIterate(I->Obj,rec,next))
              if(rec->obj->fUpdate) 
                rec->obj->fUpdate(rec->obj);
          }
      }
      PyMOL_SetBusy(G->PyMOL,false); /*  race condition -- may need to be fixed */
    } else { /* defer builds mode == 5 -- for now, only update non-molecular objects */
      /* single-threaded update */
      rec = NULL;
      while(ListIterate(I->Obj,rec,next)) {
        if(rec->obj->type != cObjectMolecule) {
          if(rec->obj->fUpdate) 
            rec->obj->fUpdate(rec->obj);
        }
      }
    }

    I->ChangedFlag = false;
    
    if((defer_builds_mode >= 2) && (force || (defer_builds_mode !=5)) &&
       (cur_state != I->LastStateBuilt)) { 
      /* purge graphics representation when no longer used */
      if(I->LastStateBuilt>=0) {
        while(ListIterate(I->Obj,rec,next)) {
          if(rec->obj->fInvalidate && 
             ((rec->obj->type != cObjectMolecule) || force || defer_builds_mode!=5)) {
            int static_singletons = SettingGet_b(G,rec->obj->Setting,NULL,cSetting_static_singletons);
            int async_builds = SettingGet_b(G,rec->obj->Setting,NULL,cSetting_async_builds);
            int max_threads =  SettingGet_i(G,rec->obj->Setting,NULL,cSetting_max_threads); 
            int nFrame = 0;
            if(rec->obj->fGetNFrame)
              nFrame = rec->obj->fGetNFrame(rec->obj);
            else
              nFrame = 0;
            if((nFrame>1)||(!static_singletons)) {
              int start = I->LastStateBuilt;
              int stop = start+1;
              int ste;
              if(async_builds&&(max_threads>1)) {
                if( (start/max_threads) == (cur_state/max_threads)) {
                  stop = start; /* don't purge current batch */
                } else {
                  int base = start/max_threads; /* now purge previous batch */
                  start = base * max_threads;
                  stop = (base+1) * max_threads;
                }
              }
              for(ste=start;ste<stop;ste++) {
                rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvPurge,ste);
              }
            }
          }
        }
      }
    } 
    I->LastStateBuilt = cur_state;
    WizardDoScene(G);
    if(!MovieDefined(G)) {
      if(SettingGetGlobal_i(G,cSetting_frame)!= (cur_state+1))
        SettingSetGlobal_i(G,cSetting_frame, (cur_state+1));
    }
  }

  PRINTFD(G,FB_Scene)
    " SceneUpdate: leaving...\n"
    ENDFD;
}
/*========================================================================*/
int SceneRenderCached(PyMOLGlobals *G)
{
  /* sets up a cached image buffer is one is available, or if we are
   * using cached images by default */
  register CScene *I=G->Scene;
  ImageType *image;
  int renderedFlag=false;

  PRINTFD(G,FB_Scene)
    " SceneRenderCached: entered.\n"
    ENDFD;

  if(I->DirtyFlag) {
    int moviePlaying = MoviePlaying(G);
    
    if(I->MovieFrameFlag||
       (moviePlaying&&SettingGet(G,cSetting_cache_frames))) {
      I->MovieFrameFlag=false;
      image = MovieGetImage(G,
                           MovieFrameToImage(G,
                                             SettingGetGlobal_i(G,cSetting_frame)-1));
      if(image) {
        if(I->Image && (!I->MovieOwnsImageFlag))
          ScenePurgeImage(G);
        I->MovieOwnsImageFlag = true;
        I->CopyType = true;
        I->Image = image;
        OrthoDirty(G);
        renderedFlag=true;
      } else {
        SceneMakeMovieImage(G,true,false,cSceneImage_Default);
        renderedFlag=true;
      }
    } else if(moviePlaying&&SettingGetGlobal_b(G,cSetting_ray_trace_frames)) {
      SceneRay(G,0,0,(int)SettingGet(G,cSetting_ray_default_renderer),
               NULL,NULL,0.0F,0.0F,false,NULL,true,-1); 
    }  else if(moviePlaying&&SettingGetGlobal_b(G,cSetting_draw_frames)) {
      SceneMakeSizedImage(G,0,0,SettingGetGlobal_i(G,cSetting_antialias));
    } else if(I->CopyType == true) { /* true vs. 2 */
      renderedFlag = true;
    } else {
      renderedFlag = false;
    }
    I->DirtyFlag=false;
  } else if(I->CopyType == true) { /* true vs. 2 */
    renderedFlag=true;
  }

  PRINTFD(G,FB_Scene)
    " SceneRenderCached: leaving...renderedFlag %d\n",renderedFlag
    ENDFD;

  return(renderedFlag);
}

float SceneGetSpecularValue(PyMOLGlobals *G,float spec,int limit)
{
  int n_light = SettingGetGlobal_i(G,cSetting_spec_count);  
  if(n_light<0)
    n_light = SettingGetGlobal_i(G,cSetting_light_count);  
  if(n_light>limit)
    n_light = limit;
  if(n_light>2) {
    spec = spec/pow(n_light-1,0.6F);
  }
  return spec;
}

float SceneGetReflectScaleValue(PyMOLGlobals *G,int limit)
{
  float result = 1.0F;
  register float _1 = 1.0F;
  int n_light = SettingGetGlobal_i(G,cSetting_light_count);  
  if(n_light>limit)
    n_light = limit;
  if(n_light>1) {
    float tmp[3];
    float sum = 0.0F;
    copy3f(SettingGetGlobal_3fv(G,cSetting_light),tmp);
    normalize3f(tmp);
    sum = _1 - tmp[2];
    if(n_light>2) {
      copy3f(SettingGetGlobal_3fv(G,cSetting_light2),tmp);
      normalize3f(tmp);
      sum += _1 - tmp[2];
      if(n_light>3) {
        copy3f(SettingGetGlobal_3fv(G,cSetting_light3),tmp);
        normalize3f(tmp);
        sum += _1 - tmp[2];
        if(n_light>4) {
          copy3f(SettingGetGlobal_3fv(G,cSetting_light4),tmp);
          normalize3f(tmp);
          sum += _1 - tmp[2];

          if(n_light>5) {
            copy3f(SettingGetGlobal_3fv(G,cSetting_light5),tmp);
            normalize3f(tmp);
            sum += _1 - tmp[2];

            if(n_light>6) {
              copy3f(SettingGetGlobal_3fv(G,cSetting_light6),tmp);
              normalize3f(tmp);
              sum += _1 - tmp[2];
              
              if(n_light>7) {
                copy3f(SettingGetGlobal_3fv(G,cSetting_light7),tmp);
                normalize3f(tmp);
                sum += _1 - tmp[2];
                
                if(n_light>8) {
                  copy3f(SettingGetGlobal_3fv(G,cSetting_light8),tmp);
                  normalize3f(tmp);
                  sum += _1 - tmp[2];
                }

                if(n_light>9) {
                  copy3f(SettingGetGlobal_3fv(G,cSetting_light9),tmp);
                  normalize3f(tmp);
                  sum += _1 - tmp[2];
                }
              }
            }
          }
        }
      }
    }
    sum *= 0.5;
    return result/sum;
  }
  return result;
}

static void SceneProgramLighting(PyMOLGlobals *G)
{

   /* load up the light positions relative to the camera while 
     MODELVIEW still has the identity */
  int n_light = SettingGetGlobal_i(G,cSetting_light_count);
  float direct = SettingGetGlobal_f(G,cSetting_direct);
  float f;
  float vv[4];
  float reflect = SceneGetReflectScaleValue(G,8) * SettingGetGlobal_f(G,cSetting_reflect);

  float spec_value = SettingGet(G,cSetting_specular);
  if(spec_value == 1.0F) {
    spec_value=SettingGet(G,cSetting_specular_intensity);
  }
  if(spec_value<R_SMALL4) spec_value = 0.0F;
  spec_value = SceneGetSpecularValue(G,spec_value,8);

  /* lighting */
  
  glEnable(GL_LIGHTING);
  
  vv[0] = 0.0F;
  vv[1] = 0.0F;
  vv[2] = 1.0F;
  /* workaround for flickering of specular reflections on Mac OSX 10.3.8 with nVidia hardware */
  vv[3] = 0.0F;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifndef _MACPYMOL_XCODE
#ifdef _PYMOL_OSX
  vv[3] = 0.000001F;
#endif
#endif
/* END PROPRIETARY CODE SEGMENT */

  glLightfv(GL_LIGHT0,GL_POSITION,vv);


  if(n_light>1) {

    /* workaround for flickering of specular reflections on Mac OSX 10.3.8 with nVidia hardware */
    vv[3] = 0.0F;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifndef _MACPYMOL_XCODE
#ifdef _PYMOL_OSX
    vv[3] = 0.000001F;
#endif
#endif
/* END PROPRIETARY CODE SEGMENT */

    copy3f(SettingGetGlobal_3fv(G,cSetting_light),vv);
    normalize3f(vv);
    invert3f(vv);
    
    glLightfv(GL_LIGHT1,GL_POSITION,vv);
    
    if(n_light>2) {
      copy3f(SettingGetGlobal_3fv(G,cSetting_light2),vv);
      normalize3f(vv);
      invert3f(vv);
      glLightfv(GL_LIGHT2,GL_POSITION,vv);
      
      if(n_light>3) {
        copy3f(SettingGetGlobal_3fv(G,cSetting_light3),vv);
        normalize3f(vv);
        invert3f(vv);
        glLightfv(GL_LIGHT3,GL_POSITION,vv);
        
        if(n_light>4) {
          copy3f(SettingGetGlobal_3fv(G,cSetting_light4),vv);
          normalize3f(vv);
          invert3f(vv);
          glLightfv(GL_LIGHT4,GL_POSITION,vv);

          if(n_light>5) {
            copy3f(SettingGetGlobal_3fv(G,cSetting_light5),vv);
            normalize3f(vv);
            invert3f(vv);
            glLightfv(GL_LIGHT5,GL_POSITION,vv);

            if(n_light>6) {
              copy3f(SettingGetGlobal_3fv(G,cSetting_light6),vv);
              normalize3f(vv);
              invert3f(vv);
              glLightfv(GL_LIGHT6,GL_POSITION,vv);

              if(n_light>7) {
                copy3f(SettingGetGlobal_3fv(G,cSetting_light7),vv);
                normalize3f(vv);
                invert3f(vv);
                glLightfv(GL_LIGHT7,GL_POSITION,vv);

              }
            }
          }
        }
      }
    }
  }  else {
    direct += reflect;
    if(direct>1.0F)
      direct = 1.0F;
  }

  if(SettingGet(G,cSetting_two_sided_lighting)||
     (SettingGetGlobal_i(G,cSetting_transparency_mode)==1)) {
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  } else {
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
  }
  
  /* ambient lighting */
  
  f=SettingGet(G,cSetting_ambient);
  vv[0] = f;
  vv[1] = f;
  vv[2] = f;
  vv[3] = 1.0F;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,vv);

  /* LIGHT0 is our direct light (eminating from the camera -- minus Z) */
  
  if(direct>R_SMALL4) {          

    glEnable(GL_LIGHT0);

    vv[0] = 0.0F;
    vv[1] = 0.0F;
    vv[2] = 0.0F;
    vv[3] = 1.0F;
    glLightfv(GL_LIGHT0,GL_AMBIENT,vv);
    
    vv[0] = direct;
    vv[1] = direct;
    vv[2] = direct;
    vv[3] = 1.0F;
    glLightfv(GL_LIGHT0,GL_DIFFUSE,vv);
    
    {
      float spec_direct = SettingGet(G,cSetting_spec_direct);
      float spec[4] = {0.0F,0.0F,0.0F,1.0F};
      if(spec_direct<0.0F) {
        spec[0] = spec[1] = spec[2] = spec_value;
        spec[3] = 1.0F;
      } else if(spec_direct>0.0F){
        spec[0] = spec[1] = spec[2] = spec_direct;
        spec[3] = 1.0F;
      }
      glLightfv(GL_LIGHT0,GL_SPECULAR,spec);
    }
    
  } else {
    glDisable(GL_LIGHT0);
  }
  
  /* LIGHTS1-3 are our reflected light (specular and diffuse
     reflections from a movable directional lights) */
  
  {
    float spec[4];
    if(n_light>1) { 
      float diff[4];
      float zero[4] = { 0.0F, 0.0F, 0.0F, 1.0F }; /* no ambient */
      int spec_count = SettingGetGlobal_i(G,cSetting_spec_count);  
      if(spec_count<0)
        spec_count = SettingGetGlobal_i(G,cSetting_light_count);  
      
      spec[0] = spec[1] = spec[2] = spec_value;
      spec[3] = 1.0F;
      diff[0] = diff[1] = diff[2] = reflect;
      diff[3] = 1.0F;
      glEnable(GL_LIGHT1);
      if(spec_count>=1) {
        glLightfv(GL_LIGHT1,GL_SPECULAR,spec);
      } else {
        glLightfv(GL_LIGHT1,GL_SPECULAR,zero);
      }
      glLightfv(GL_LIGHT1,GL_AMBIENT,zero);
      glLightfv(GL_LIGHT1,GL_DIFFUSE,diff);
      if(n_light>2) {
        glEnable(GL_LIGHT2);
        if(spec_count>=2) {
          glLightfv(GL_LIGHT2,GL_SPECULAR,spec);
        } else {
          glLightfv(GL_LIGHT2,GL_SPECULAR,zero);
        }
        glLightfv(GL_LIGHT2,GL_AMBIENT,zero);
        glLightfv(GL_LIGHT2,GL_DIFFUSE,diff);
        if(n_light>3) {
          glEnable(GL_LIGHT3);
          if(spec_count>=3) {
            glLightfv(GL_LIGHT3,GL_SPECULAR,spec);
          } else {
            glLightfv(GL_LIGHT3,GL_SPECULAR,zero);
          }
          glLightfv(GL_LIGHT3,GL_AMBIENT,zero);
          glLightfv(GL_LIGHT3,GL_DIFFUSE,diff);
          if(n_light>4) {
            glEnable(GL_LIGHT4);
            if(spec_count>=4) {
              glLightfv(GL_LIGHT4,GL_SPECULAR,spec);
            } else {
              glLightfv(GL_LIGHT4,GL_SPECULAR,zero);
            }
            glLightfv(GL_LIGHT4,GL_AMBIENT,zero);
            glLightfv(GL_LIGHT4,GL_DIFFUSE,diff);
            if(n_light>5) {
              glEnable(GL_LIGHT5);
              if(spec_count>=5) {
                glLightfv(GL_LIGHT5,GL_SPECULAR,spec);
              } else {
                glLightfv(GL_LIGHT5,GL_SPECULAR,zero);
              }
              glLightfv(GL_LIGHT5,GL_AMBIENT,zero);
              glLightfv(GL_LIGHT5,GL_DIFFUSE,diff);
              if(n_light>6) {
                glEnable(GL_LIGHT6);
                if(spec_count>=6) {
                  glLightfv(GL_LIGHT6,GL_SPECULAR,spec);
                } else {
                  glLightfv(GL_LIGHT6,GL_SPECULAR,zero);
                }
                glLightfv(GL_LIGHT6,GL_AMBIENT,zero);
                glLightfv(GL_LIGHT6,GL_DIFFUSE,diff);
                if(n_light>7) {
                  glEnable(GL_LIGHT7);
                  if(spec_count>=7) {
                    glLightfv(GL_LIGHT7,GL_SPECULAR,spec);
                  } else {
                    glLightfv(GL_LIGHT7,GL_SPECULAR,zero);
                  }
                  glLightfv(GL_LIGHT7,GL_AMBIENT,zero);
                  glLightfv(GL_LIGHT7,GL_DIFFUSE,diff);
                }
              }
            }
          }
        }
      }
    }
    if(n_light<2) glDisable(GL_LIGHT1);
    if(n_light<3) glDisable(GL_LIGHT2);
    if(n_light<4) glDisable(GL_LIGHT3);
    if(n_light<5) glDisable(GL_LIGHT4);
    if(n_light<6) glDisable(GL_LIGHT5);
    if(n_light<7) glDisable(GL_LIGHT6);
    if(n_light<8) glDisable(GL_LIGHT7);
  }

  {
    float ones[4] = {1.0F,1.0F,1.0F,1.0F};
    glMaterialfv(GL_FRONT,GL_SPECULAR,ones);
  }

  glMaterialf(GL_FRONT,GL_SHININESS,SettingGet(G,cSetting_shininess));

  if(0) {
    
    glColor4f(1.0,1.0,1.0,1.0);
    
    glGetLightfv(GL_LIGHT0,GL_POSITION,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_POSITION,vv)"); 
      
    glGetLightfv(GL_LIGHT1,GL_POSITION,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_POSITION,vv)");
    
    glGetFloatv(GL_LIGHT_MODEL_AMBIENT,vv);
    dump4f(vv,"glGetFloatv(GL_LIGHT_MODEL_AMBIENT,vv)");
    
    glGetFloatv(GL_LIGHT0,vv);
    printf("glGetFloatv(GL_LIGHT0) %8.3f\n",vv[0]);
    
    glGetLightfv(GL_LIGHT0,GL_AMBIENT,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_AMBIENT,vv)");
    
    glGetLightfv(GL_LIGHT0,GL_DIFFUSE,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_DIFFUSE,vv)");
    
    glGetLightfv(GL_LIGHT0,GL_SPECULAR,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_SPECULAR,vv)");
    
    glGetFloatv(GL_LIGHT1,vv);
    printf("glGetFloatv(GL_LIGHT1) %8.3f\n",vv[0]);
    
    glGetLightfv(GL_LIGHT1,GL_AMBIENT,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_AMBIENT,vv)");
    
    glGetLightfv(GL_LIGHT1,GL_DIFFUSE,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_DIFFUSE,vv)");
    
    glGetLightfv(GL_LIGHT1,GL_SPECULAR,vv);
    dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_SPECULAR,vv)");
    
    glGetFloatv(GL_LIGHT_MODEL_AMBIENT,vv);
    dump4f(vv, "glGetFloatv(GL_LIGHT_MODEL_AMBIENT,vv)");
    
    glGetMaterialfv(GL_FRONT, GL_AMBIENT, vv);
    dump4f(vv, "glGetMaterialfv(GL_FRONT,GL_AMBIENT,vv)");
    
    glGetMaterialfv(GL_FRONT, GL_DIFFUSE, vv);
    dump4f(vv, "glGetMaterialfv(GL_FRONT,GL_DIFFUSE,vv)"); 
    
    glGetMaterialfv(GL_FRONT, GL_SPECULAR, vv);
    dump4f(vv, "glGetMaterialfv(GL_FRONT,GL_SPECULAR,vv)");
    
    glGetMaterialfv(GL_FRONT, GL_SHININESS, vv);
    printf("glGetMaterialfv(GL_FRONT,GL_SHININESS, vv) %8.3f\n",vv[0]);
    
  }
}
/*========================================================================*/
static void SceneRenderAll(PyMOLGlobals *G,SceneUnitContext *context,
                           float *normal,Picking **pickVLA,
                           int pass,int fat, float width_scale,
                           GridInfo *grid, int dynamic_pass)
{
  register CScene *I=G->Scene;
  ObjRec *rec=NULL;
  float vv[4];
  int state = SceneGetState(G);
  RenderInfo info;
  UtilZeroMem(&info,sizeof(RenderInfo));
  info.pick = pickVLA;
  info.pass = pass;
  info.vertex_scale = I->VertexScale;
  info.fog_start = I->FogStart;
  info.fog_end = I->FogEnd;
  info.pmv_matrix = I->PmvMatrix;
  info.front = I->FrontSafe;
  info.sampling = 1;
  info.alpha_cgo = I->AlphaCGO;
  info.ortho = SettingGetGlobal_b(G,cSetting_ortho);
  if(I->StereoMode && dynamic_pass && (!info.pick)) {
    int stereo_mode = SettingGetGlobal_i(G,cSetting_stereo_mode);
    switch(stereo_mode) {
    case cStereo_dynamic:
    case cStereo_clone_dynamic:
      info.line_lighting = true;
      break;
    }
  }
  if(I->StereoMode) {
    float buffer;
    float stAng,stShift;
    stAng = SettingGet(G,cSetting_stereo_angle);
    stShift = SettingGet(G,cSetting_stereo_shift);
    stShift = (float)(stShift*fabs(I->Pos[2])/100.0);
    stAng = (float)(stAng*atan(stShift/fabs(I->Pos[2]))*90.0/cPI);
    buffer = fabs(I->Width*I->VertexScale*tan(cPI*stAng/180.0));
    info.stereo_front = I->FrontSafe + buffer;
  } else {
    info.stereo_front = I->FrontSafe;
  }
  info.back = I->BackSafe;
  SceneGetViewNormal(G,info.view_normal);

  if(info.alpha_cgo && (pass == 1)) {
    CGOReset(info.alpha_cgo);
    CGOSetZVector(info.alpha_cgo, I->ModMatrix[2], I->ModMatrix[6], I->ModMatrix[10]);
  }

  if(SettingGetGlobal_b(G,cSetting_dynamic_width)) {
    info.dynamic_width = true;
    info.dynamic_width_factor = SettingGetGlobal_f(G,cSetting_dynamic_width_factor);
    info.dynamic_width_min = SettingGetGlobal_f(G,cSetting_dynamic_width_min);
    info.dynamic_width_max = SettingGetGlobal_f(G,cSetting_dynamic_width_max);
  }

  if(width_scale!=0.0F) {
    info.width_scale_flag = true;
    info.width_scale = width_scale;
    info.sampling = (int)info.width_scale;
    if(info.sampling<1)
      info.sampling = 1;
  }

  {
    int *slot_vla = I->SlotVLA;
    while(ListIterate(I->Obj,rec,next)) {
      if(rec->obj->fRender) {

        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("Before fRender iteration");

        if(SceneGetDrawFlag(grid, slot_vla, rec->obj->grid_slot)) {
          glPushMatrix();
          if(fat)
            glLineWidth(3.0);
          
          switch(rec->obj->Context) {
          case 1: /* unit context */
            {
#ifndef _PYMOL_OSX
              /* workaround for MacOSX 10.4.3 */
              glPushAttrib(GL_LIGHTING_BIT); 
#endif
              
              glMatrixMode(GL_PROJECTION);
              glPushMatrix();
              glLoadIdentity();
              glMatrixMode(GL_MODELVIEW);
              glPushMatrix();
              glLoadIdentity();
              vv[0]=0.0;
              vv[1]=0.0;
              vv[2]=-1.0;
              vv[3]=0.0;
              glLightfv(GL_LIGHT0,GL_POSITION,vv);
              glLightfv(GL_LIGHT1,GL_POSITION,vv);

              if(!grid->active) {
                glOrtho(context->unit_left,
                        context->unit_right,
                        context->unit_top,
                        context->unit_bottom,
                        context->unit_front,
                        context->unit_back);
              } else { /* special unit context */
                glOrtho(grid->context.unit_left,
                        grid->context.unit_right,
                        grid->context.unit_top,
                        grid->context.unit_bottom,
                        grid->context.unit_front,
                        grid->context.unit_back);
              }
              
              glNormal3f(0.0F,0.0F,1.0F);
              info.state = ObjectGetCurrentState(rec->obj,false);
              rec->obj->fRender(rec->obj,&info);
              
              glMatrixMode(GL_PROJECTION);
              glPopMatrix();
              glMatrixMode(GL_MODELVIEW);
              glLoadIdentity();
#ifndef _PYMOL_OSX
              glPopAttrib();
#else  
              /* workaround for MacOSX 10.4.3 */
              /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
              SceneProgramLighting(G); /* an expensive workaround... */
              if(pickVLA) {
                glDisable(GL_FOG);
                glDisable(GL_COLOR_MATERIAL);
                glDisable(GL_LIGHTING);
                glDisable(GL_DITHER);
                glDisable(GL_BLEND);
                glDisable(GL_LINE_SMOOTH);
                glDisable(GL_POLYGON_SMOOTH);
                if(G->Option->multisample)    
                  glDisable(0x809D); /* GL_MULTISAMPLE_ARB */
                glShadeModel(GL_FLAT);
              }
              /* END PROPRIETARY CODE SEGMENT */
#endif
              glPopMatrix();
            }
            break;
          case 2:
            break;
          case 0: /* context/grid 0 is all slots */
          default:
        if(Feedback(G,FB_OpenGL,FB_Debugging))
            if(normal) 
              glNormal3fv(normal);
            if((!grid->active)||(grid->mode!=2)) {
              info.state = ObjectGetCurrentState(rec->obj,false);
              rec->obj->fRender(rec->obj,&info);
            } else if(grid->slot) {
              if ( (info.state = state + grid->slot - 1) >= 0 )
                rec->obj->fRender(rec->obj,&info);              
            }
            break;
          }
          glPopMatrix();
        }
        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("After fRender iteration");

      }
    }
  }

  if(info.alpha_cgo) { 
    CGOStop(info.alpha_cgo);
    /* this only works when all objects are rendered in the same frame of reference */
    if(pass == -1) {
      CGORenderGLAlpha(info.alpha_cgo, &info);
    }
  }
}

#ifdef _PYMOL_SHARP3D
void sharp3d_begin_left_stereo(void);
void sharp3d_switch_to_right_stereo(void);
void sharp3d_end_stereo(void);
#endif

/*========================================================================*/
void SceneRender(PyMOLGlobals *G,Picking *pick,int x,int y,
                 Multipick *smp,int oversize_width, int oversize_height,
                 int click_side,int force_copy)
{
  /* think in terms of the camera's world */
  register CScene *I=G->Scene;
  float fog[4];
  float *v;
  unsigned int lowBits,highBits;
  unsigned int *lowBitVLA=NULL,*highBitVLA=NULL;
  int high,low;
  float zAxis[4] = { 0.0, 0.0, 1.0, 0.0 };
  float normal[4] = { 0.0, 0.0, 1.0, 0.0 };
  float aspRat = ((float) I->Width) / ((float) I->Height);
  float height,width;
  double start_time=0.0;
  int view_save[4];
  Picking *pickVLA,*pik;
  int lastIndex=0;
  void *lastPtr=NULL;
  int index;
  int curState;
  int nPick,nHighBits,nLowBits;
  int pass;
  float fov;
  int must_render_stereo = false;
  int mono_as_quad_stereo = false;
  int stereo_using_mono_matrix = false;
  int debug_pick = 0;
  GLenum render_buffer;
  SceneUnitContext context;
  float width_scale = 0.0F;
  int stereo_mode = I->StereoMode;
  GridInfo grid;
  int grid_mode = SettingGetGlobal_i(G,cSetting_grid_mode);
  int fog_active = false;

  PRINTFD(G,FB_Scene)
    " SceneRender: entered. pick %p x %d y %d smp %p\n",
    (void*)pick,x,y,(void*)smp
    ENDFD;

  UtilZeroMem(&grid,sizeof(GridInfo));
  if(grid_mode) {
    int grid_size = SceneGetGridSize(G,grid_mode);
    GridUpdate(&grid, aspRat, grid_mode, grid_size);
    if(grid.active) 
      aspRat *= grid.asp_adjust;
  }
  
  SceneUpdateAnimation(G);

  if(SceneMustDrawBoth(G)) {
    render_buffer = GL_BACK_LEFT;
  } else {
    render_buffer = GL_BACK;
  }

  switch(stereo_mode) {
  case cStereo_walleye:
  case cStereo_crosseye:
    aspRat=aspRat/2;
    break;
  }

  fov=SettingGet(G,cSetting_field_of_view);
  if(G->HaveGUI && G->ValidContext) {
    
    if(Feedback(G,FB_OpenGL,FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender checkpoint 0");

    must_render_stereo = (stereo_mode!=0); /* are we doing stereo? */
    if(!must_render_stereo) {
      if(G->StereoCapable &&
        SettingGet_i(G,NULL,NULL,cSetting_stereo_double_pump_mono)) {
        /* force stereo rendering */
            must_render_stereo=true;
        if(stereo_mode==0) {
          mono_as_quad_stereo = true; /* rendering stereo as mono */
          stereo_using_mono_matrix = true;
        }
      } else {
        int st_mode = SettingGet_i(G,NULL,NULL,cSetting_stereo_mode);
        if(st_mode==cStereo_geowall) {
          stereo_mode = st_mode;
          must_render_stereo = true;
          stereo_using_mono_matrix = true;
        }
      }
    }

    /* if we seem to be configured for hardware stereo, 
        but can't actually do it, then fallback on mono -- 
        this would happen for instance if fullscreen is stereo-component
        and windowed is not */
    
    if(must_render_stereo && (stereo_mode<cStereo_crosseye) && !(G->StereoCapable)) {
      must_render_stereo = false;
      mono_as_quad_stereo = false;
    }

    if(must_render_stereo && stereo_via_stencil(stereo_mode)) {
      if(!I->StencilValid) {
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT,viewport);
    
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0,viewport[2],0,viewport[3],-10.0,10.0);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glTranslatef(0.33F,0.33F,0.0F); 
        
        glDisable(GL_ALPHA_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_FOG);
        glDisable(GL_NORMALIZE);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_COLOR_MATERIAL);
        glDisable(GL_LINE_SMOOTH);
        glDisable(GL_DITHER);
        glDisable(GL_BLEND);
        glShadeModel(GL_SMOOTH);
        glDisable(0x809D); /* GL_MULTISAMPLE_ARB */ 

        glDisable(GL_STENCIL_TEST);
        glClearStencil(0);
        glColorMask(false,false,false,false);
        glDepthMask(false);
        glClear(GL_STENCIL_BUFFER_BIT);

        glEnable(GL_STENCIL_TEST);
        glStencilFunc(GL_ALWAYS, 1, 1);
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

        {
          int h = viewport[3], w=viewport[2];
          glLineWidth(1.0);
          switch(stereo_mode) {
          case cStereo_stencil_by_row:
            {
              int y;
              glBegin(GL_LINES);
              for(y=0;y<h;y+=2) {
                glVertex2i(0,y);
                glVertex2i(w,y);
              }
              glEnd();
            }
            break;
          case cStereo_stencil_by_column:
            {
              int x;
              glBegin(GL_LINES);
              for(x=0;x<w;x+=2) {
                glVertex2i(x,0);
                glVertex2i(x,h);
              }
              glEnd();
            }
            break;
          case cStereo_stencil_checkerboard:
            {
              int i,m = 2* ((h>w) ? h : w);
              glBegin(GL_LINES);
              for(i=0;i<m;i+=2) {
                glVertex2i(i,0);
                glVertex2i(0,i);
              }
              glEnd();
            }
            break;
          }
        }

        glColorMask(true,true,true,true);
        glDepthMask(true);

        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();

        I->StencilValid = true;
      }
    }
        
    if(must_render_stereo) {
      if(mono_as_quad_stereo) { /* double-pumped mono */
        OrthoDrawBuffer(G,GL_BACK_LEFT);
        render_buffer = GL_BACK_LEFT;
      } else {
        switch(stereo_mode) {
        case cStereo_quadbuffer: /* hardware stereo */
        case cStereo_clone_dynamic:
          OrthoDrawBuffer(G,GL_BACK_LEFT);
          render_buffer = GL_BACK_LEFT;
          break;
        default:  /* some kind of software stereo */
          OrthoDrawBuffer(G,GL_BACK);
          render_buffer = GL_BACK;
          break;
        }
      }
    } else { /* normal mono rendering */
      OrthoDrawBuffer(G,GL_BACK);
      render_buffer = GL_BACK;
    }

    if(Feedback(G,FB_OpenGL,FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender checkpoint 1");

    glGetIntegerv(GL_VIEWPORT,(GLint*)(void*)view_save);

    if(oversize_width && oversize_height) {
      int want_view[4];
      int got_view[4];
      want_view[0] = I->Block->rect.left+x;
      want_view[1] = I->Block->rect.bottom+y;
      want_view[2] = oversize_width;
      want_view[3] = oversize_height;
      glViewport(want_view[0],want_view[1],want_view[2],want_view[3]);
      glGetIntegerv(GL_VIEWPORT,(GLint*)(void*)got_view);
      if((got_view[0]!=want_view[0])||
         (got_view[1]!=want_view[1])||
         (got_view[2]!=want_view[2])||
         (got_view[3]!=want_view[3])) {
        PRINTFB(G,FB_Scene,FB_Warnings)
          "Scene-Warning: glViewport failure.\n"
          ENDFB(G);
      }
      switch(stereo_mode) {
      case cStereo_geowall: 
        stereo_mode = 0;
        break;
      }
      stereo_using_mono_matrix = true;
      width_scale = ((float)(oversize_width))/I->Width;
    } else {
      glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height);
    }

    debug_pick = (int)SettingGet(G,cSetting_debug_pick);

    if(SettingGet(G,cSetting_line_smooth)) {
      if(!(pick||smp)) {
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
      }
    } else {
      glDisable(GL_LINE_SMOOTH);
    }
    glLineWidth(SettingGet(G,cSetting_line_width));

    glPointSize(SettingGet(G,cSetting_dot_width));

    glEnable(GL_NORMALIZE); /* get rid of this to boost performance */

    glEnable(GL_DEPTH_TEST);

    /* get matrixes for unit objects */

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    SceneProgramLighting(G); /* must be done with identity MODELVIEW */

    ScenePrepareUnitContext(&context,I->Width,I->Height);

    /* do standard 3D objects */

    /* Set up the clipping planes */
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if(SettingGet(G,cSetting_all_states)) {
      curState=-1;
    } else {
      curState=SettingGetGlobal_i(G,cSetting_state)-1;
    }

    if(!SettingGetGlobal_b(G,cSetting_ortho)) {
      gluPerspective(fov,aspRat,I->FrontSafe,I->BackSafe);
    } else {
      height  = (float)(fabs(I->Pos[2])*tan((fov/2.0)*cPI/180.0));     
      width = height*aspRat;
    
      glOrtho(-width,width,-height,height,
              I->FrontSafe,I->BackSafe);

    }

    glMatrixMode(GL_MODELVIEW);
    ScenePrepareMatrix(G,0);

    /* Save these for editing operations */

    glGetFloatv(GL_MODELVIEW_MATRIX,I->ModMatrix);
    glGetFloatv(GL_PROJECTION_MATRIX,I->ProMatrix);

    multiply44f44f44f(I->ModMatrix,I->ProMatrix,I->PmvMatrix);
    
    /* get the Z axis vector for sorting transparent objects */
    
    if(SettingGetGlobal_b(G,cSetting_transparency_global_sort) &&
       SettingGetGlobal_b(G,cSetting_transparency_mode) ) {
      if(!I->AlphaCGO)
        I->AlphaCGO = CGONew(G);
    } else if(I->AlphaCGO) {
      CGOFree(I->AlphaCGO);
      I->AlphaCGO = NULL;
    }

    /* make note of how large pixels are at the origin  */

    I->VertexScale = SceneGetScreenVertexScale(G,I->Origin);

    /* determine the direction in which we are looking relative*/

    /* 2. set the normals to reflect light back at the camera */

    MatrixInvTransformC44fAs33f3f(I->RotMatrix,zAxis,normal); 
    copy3f(normal,I->ViewNormal);
  
    if(SettingGet(G,cSetting_normal_workaround)) {
      I->LinesNormal[0]=0.0;    
      I->LinesNormal[1]=0.0;     
      I->LinesNormal[2]=1.0;
      /* for versions of GL that don't transform GL_LINES normals */
    } else {
      I->LinesNormal[0]=I->ViewNormal[0];
      I->LinesNormal[1]=I->ViewNormal[1];
      I->LinesNormal[2]=I->ViewNormal[2];
    }

    if(!(pick||smp)) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

      glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

      glEnable(GL_COLOR_MATERIAL);

      glShadeModel(GL_SMOOTH);

      glEnable(GL_DITHER);

      glAlphaFunc(GL_GREATER, 0.05F);
      glEnable(GL_ALPHA_TEST);

      if(G->Option->multisample)
        glEnable(0x809D); /* GL_MULTISAMPLE_ARB */

      I->FogStart = (I->BackSafe-I->FrontSafe)*SettingGet(G,cSetting_fog_start)+I->FrontSafe;

      glFogf(GL_FOG_MODE, GL_LINEAR);
      glFogf(GL_FOG_START, I->FogStart);
      

      {
        float fog_density = SettingGet(G,cSetting_fog);
        if((fog_density>R_SMALL8) && (fog_density!=1.0F)) {
          I->FogEnd = I->FogStart + (I->BackSafe - I->FogStart)/fog_density;
        } else {
          I->FogEnd = I->BackSafe;          
        }
        glFogf(GL_FOG_END, I->FogEnd);
        glFogf(GL_FOG_DENSITY, fog_density);
      }

      
      v=SettingGetfv(G,cSetting_bg_rgb);
      fog[0]=v[0];
      fog[1]=v[1];
      fog[2]=v[2];
      
      /* NOTE: this doesn't seem to work :( -- only raytracing can do this */
      fog[3]= (SettingGetGlobal_b(G,cSetting_opaque_background) ? 1.0F : 0.0F);
      
      glFogfv(GL_FOG_COLOR, fog);


      if(SettingGetGlobal_b(G,cSetting_depth_cue) && 
         (SettingGet(G,cSetting_fog)!=0.0F)) {
        fog_active = true;
        glEnable(GL_FOG);
      } else {
        fog_active = false;
        glDisable(GL_FOG);
      }
      glColor4ub(255,255,255,255);
      glNormal3fv(normal);
      
    } else {
      /* picking mode: we want flat, unshaded, unblended, unsmooth colors */

      glDisable(GL_FOG);
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_LIGHTING);
      glDisable(GL_DITHER);
      glDisable(GL_BLEND);
      glDisable(GL_LINE_SMOOTH);
      glDisable(GL_POLYGON_SMOOTH);
      if(G->Option->multisample)    
        glDisable(0x809D); /* GL_MULTISAMPLE_ARB */
      glShadeModel(GL_FLAT);

    }

    PRINTFD(G,FB_Scene)
    " SceneRender: matrices loaded. rendering objects...\n"
    ENDFD;
    
    /* 1. render all objects */
    if(pick||smp) {
      
      switch(stereo_mode) {
      case cStereo_crosseye:
      case cStereo_walleye:
      case cStereo_sidebyside:
        glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
        break;
      case cStereo_geowall:
        click_side = OrthoGetWrapClickSide(G);
        break;
      }
        
      glPushMatrix(); /* 1 */
      {
        if(!stereo_using_mono_matrix)
          switch(stereo_mode) {
          case cStereo_crosseye:
            ScenePrepareMatrix(G, (click_side > 0) ? 1 : 2 );
            break;
          case cStereo_walleye:
          case cStereo_geowall:
          case cStereo_sidebyside:
            ScenePrepareMatrix(G, (click_side < 0) ? 1 : 2 );
            break;
          }
      }

      if(pick) {
        /* atom picking HACK - obfuscative coding */
          
        glClearColor(0.0,0.0,0.0,0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                    
        pickVLA=VLACalloc(Picking,5000);
        pickVLA[0].src.index = 0;
        pickVLA[0].src.bond = 0;

        if(grid.active)
          GridGetGLViewport(&grid);

        {
          int slot;
          for(slot=0;slot<=grid.last_slot;slot++) {
            if(grid.active) { 
              GridSetGLViewport(&grid,slot);
            }
            SceneRenderAll(G,&context,NULL,&pickVLA,0,true,0.0F,&grid,0);
          }
          if(grid.active)
            GridSetGLViewport(&grid,-1);
        }
          
        if(debug_pick) {
          PyMOL_SwapBuffers(G->PyMOL);
          PSleep(G,1000000*debug_pick/4);
          PyMOL_SwapBuffers(G->PyMOL);
        }
        lowBits = SceneFindTriplet(G,x,y,render_buffer);
          
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
          
        pickVLA[0].src.index = 0;
        pickVLA[0].src.bond = 1;
          
        {
          int slot;
          for(slot=0;slot<=grid.last_slot;slot++) {
            if(grid.active) { 
              GridSetGLViewport(&grid,slot);
            }
            SceneRenderAll(G,&context,NULL,&pickVLA,0,true,0.0F,&grid,0);
          }
          if(grid.active)
            GridSetGLViewport(&grid,-1);
        }

          
        if(debug_pick) {
          PyMOL_SwapBuffers(G->PyMOL);
          PSleep(G,1000000*debug_pick/4);
          PyMOL_SwapBuffers(G->PyMOL);
        }
          
        highBits = SceneFindTriplet(G,x,y,render_buffer);
        index = lowBits+(highBits<<12);
          
        if(debug_pick) {
          PRINTFB(G,FB_Scene,FB_Details)
            " SceneClick-Detail: index %d < %d?\n",index,pickVLA[0].src.index
            ENDFB(G);
        }
          
        if(index&&(index<=pickVLA[0].src.index)) {
          *pick = pickVLA[index]; /* return object info */
          if(debug_pick) {
            PRINTFB(G,FB_Scene,FB_Details)
              " SceneClick-Detail: obj %p index %d bond %d\n",
              pick->context.object,
              pick->src.index,pick->src.bond
              ENDFB(G);
          }
        } else {
          pick->context.object = NULL;
        }
          
        VLAFree(pickVLA);
      } else if(smp) {
          
        /* multiple atom picking HACK - even more obfuscative coding */
        
        glClearColor(0.0,0.0,0.0,0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        pickVLA=VLACalloc(Picking,5000);
        pickVLA[0].src.index = 0;
        pickVLA[0].src.bond = 0; /* this is just a flag for first pass */
        
        if(grid.active)
          GridGetGLViewport(&grid);
        {
          int slot;
          for(slot=0;slot<=grid.last_slot;slot++) {
            if(grid.active) { 
              GridSetGLViewport(&grid,slot);
            }
            SceneRenderAll(G,&context,NULL,&pickVLA,0,true,0.0F,&grid,0);
          }
          if(grid.active)
            GridSetGLViewport(&grid,-1);
        }
        
        lowBitVLA = SceneReadTriplets(G,smp->x,smp->y,smp->w,smp->h,render_buffer);
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        pickVLA[0].src.index = 0;
        pickVLA[0].src.bond = 1; /* this is just a flag for second pass */
        
        {
          int slot;
          for(slot=0;slot<=grid.last_slot;slot++) {
            if(grid.active) { 
              GridSetGLViewport(&grid,slot);
            }
            SceneRenderAll(G,&context,NULL,&pickVLA,0,true,0.0F,&grid,0);
          }
          if(grid.active)
            GridSetGLViewport(&grid,-1);
        }
        
        highBitVLA = SceneReadTriplets(G,smp->x,smp->y,smp->w,smp->h,render_buffer);
        
        nLowBits = VLAGetSize(lowBitVLA);
/* need to scissor this */        nHighBits = VLAGetSize(highBitVLA);
        nPick=0;
        if(nLowBits&&nHighBits) {
          low = 0;
          high = 0;
          while((low<nLowBits)&&(high<nHighBits)) {
            
            if(lowBitVLA[low+1]==highBitVLA[high+1]) {
              index = lowBitVLA[low]+(highBitVLA[high]<<12);
              if(index&&(index<=pickVLA[0].src.index)) {          
                pik = pickVLA+index; /* just using as a tmp */
                if((pik->src.index!=lastIndex)||(pik->context.object!=lastPtr))
                  {
                    if(((CObject*)pik->context.object)->type==cObjectMolecule) {
                      nPick++; /* start from 1 */
                      VLACheck(smp->picked,Picking,nPick);
                      smp->picked[nPick] = *pik; /* return atom/object info -- will be redundant */
                    }
                    lastIndex=pik->src.index;                
                    lastPtr=pik->context.object;
                  }
              }
              low+=2;
              high+=2;
            } else if(lowBitVLA[low+1]<highBitVLA[high+1])
              low+=2;
            else 
              high+=2;
          }
        }
        
        smp->picked[0].src.index=nPick;
        
        VLAFree(pickVLA);
        VLAFreeP(lowBitVLA);
        VLAFreeP(highBitVLA);
      }
      glPopMatrix(); /* 1 */

    } else {
      
      /* STANDARD RENDERING */

      /* rendering for visualization */

      int times = 1;

      switch(stereo_mode) {
      case cStereo_clone_dynamic:
      case cStereo_dynamic:
        times = 2;
        break;
      }

      PRINTFD(G,FB_Scene)
        " SceneRender: I->StereoMode %d must_render_stereo %d\n    mono_as_quad_stereo %d  StereoCapable %d\n",
        stereo_mode, must_render_stereo, mono_as_quad_stereo, G->StereoCapable
        ENDFD;

      start_time = UtilGetSeconds(G);
      while(times--) {
        if(must_render_stereo) {
          /* STEREO RENDERING (real or double-pumped) */

          PRINTFD(G,FB_Scene)
            " SceneRender: left hand stereo...\n"
            ENDFD;

          if(Feedback(G,FB_OpenGL,FB_Debugging))
            PyMOLCheckOpenGLErr("before stereo glViewport 1");
        
          /* LEFT HAND STEREO */

          if(mono_as_quad_stereo) {
            OrthoDrawBuffer(G,GL_BACK_LEFT);
          } else switch(stereo_mode) {
          case cStereo_quadbuffer: /* hardware */
            OrthoDrawBuffer(G,GL_BACK_LEFT);
            break;
          case cStereo_crosseye: /* side by side, crosseye */
            if(oversize_width && oversize_height) {
              glViewport(I->Block->rect.left+oversize_width/2+x,
                         I->Block->rect.bottom+y,
                         oversize_width/2, oversize_height);
            } else {
              glViewport(I->Block->rect.left+I->Width/2,I->Block->rect.bottom,I->Width/2,I->Height);
            }
            break;
          case cStereo_walleye:
          case cStereo_sidebyside:
            if(oversize_width && oversize_height) {
              glViewport(I->Block->rect.left+x,
                         I->Block->rect.bottom+y,
                         oversize_width/2, oversize_height);
            } else {
              glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
            }
            break;
          case cStereo_geowall:
            glViewport(I->Block->rect.left,
                       I->Block->rect.bottom,I->Width,I->Height);
            break;          
          case cStereo_stencil_by_row:
          case cStereo_stencil_by_column:
          case cStereo_stencil_checkerboard:
            if(I->StencilValid) {
              glStencilFunc(GL_EQUAL, 1, 1);
              glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
              glEnable(GL_STENCIL_TEST);
            }
            break;
          case cStereo_stencil_custom:
#ifdef _PYMOL_SHARP3D
            sharp3d_begin_left_stereo();
#endif
            break;
          case cStereo_anaglyph:
            glClear(GL_ACCUM_BUFFER_BIT);
            break;
          case cStereo_clone_dynamic:
            glClear(GL_ACCUM_BUFFER_BIT);
            OrthoDrawBuffer(G,GL_BACK_LEFT);
            if(times) {
              float dynamic_strength = SettingGetGlobal_f(G,cSetting_stereo_dynamic_strength);
              float vv[4] = {0.75F, 0.75F, 0.75F, 1.0F};
              vv[0] =  dynamic_strength;
              vv[1] =  dynamic_strength;
              vv[2] =  dynamic_strength;
              glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,vv);
              glAccum(GL_ADD,0.5);
              glDisable(GL_FOG);
            }
            break;
           case cStereo_dynamic:
            if(times) {
              float dynamic_strength = SettingGetGlobal_f(G,cSetting_stereo_dynamic_strength);
              float vv[4] = {0.75F, 0.75F, 0.75F, 1.0F};
              vv[0] =  dynamic_strength;
              vv[1] =  dynamic_strength;
              vv[2] =  dynamic_strength;
              glClearAccum(0.5, 0.5, 0.5, 0.5);
              glClear(GL_ACCUM_BUFFER_BIT);
              glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,vv);
              glDisable(GL_FOG);
              glViewport(I->Block->rect.left+G->Option->winX/2,
                         I->Block->rect.bottom,I->Width,I->Height);
            } else {
              glClearAccum(0.0,0.0,0.0,0.0);
              glClear(GL_ACCUM_BUFFER_BIT);
              glViewport(I->Block->rect.left,
                         I->Block->rect.bottom,I->Width,I->Height);
            }
            break;
          }

          /* prepare the stereo transformation matrix */
        
          glPushMatrix(); /* 1 */
          ScenePrepareMatrix(G,stereo_using_mono_matrix ? 0 : 1);
        
          if(grid.active)
            GridGetGLViewport(&grid);
          {
            int slot;
            for(slot=0;slot<=grid.last_slot;slot++) {
            
              if(grid.active) { 
                GridSetGLViewport(&grid,slot);
              }
            
              /* render picked atoms */
            
              glPushMatrix(); /* 2 */
              EditorRender(G,curState);
              glPopMatrix(); /* 1 */
            
              /* render the debugging CGO */
            
              glPushMatrix();  /* 2 */
              glNormal3fv(normal);
              CGORenderGL(G->DebugCGO,NULL,NULL,NULL,NULL);
              glPopMatrix();  /* 1 */
            
              /* render all objects */
            
              glPushMatrix(); /* 2 */
            
              for(pass=1;pass>-2;pass--) { /* render opaque, then antialiased, then transparent...*/
                SceneRenderAll(G,&context,normal,NULL,pass,false,width_scale,&grid,times);
              }
              glPopMatrix(); /* 1 */
            
              /* render selections */
              glPushMatrix(); /* 2 */
              glNormal3fv(normal);
              ExecutiveRenderSelections(G,curState);
              glPopMatrix(); /* 1 */
            
            }
          }
          if(grid.active)
            GridSetGLViewport(&grid,-1);
        
          glPopMatrix(); /* 0 */

          /* RIGHT HAND STEREO */
        
          PRINTFD(G,FB_Scene)
            " SceneRender: right hand stereo...\n"
            ENDFD;
        
          if(Feedback(G,FB_OpenGL,FB_Debugging))
            PyMOLCheckOpenGLErr("before stereo glViewport 2");

          if(mono_as_quad_stereo) { /* double pumped mono */
            OrthoDrawBuffer(G,GL_BACK_RIGHT);
          } else switch(stereo_mode) {
          case cStereo_quadbuffer: /* hardware */
            OrthoDrawBuffer(G,GL_BACK_RIGHT);
            break;
          case cStereo_crosseye: /* side by side, crosseye */
            if(oversize_width && oversize_height) {
              glViewport(I->Block->rect.left+x,
                         I->Block->rect.bottom+y,
                         oversize_width/2, oversize_height);
            } else {
              glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
            }
            break;
          case cStereo_walleye: /* side by side, walleye */
          case cStereo_sidebyside:
            if(oversize_width && oversize_height) {
              glViewport(I->Block->rect.left+oversize_width/2+x,
                         I->Block->rect.bottom+y,
                         oversize_width/2, oversize_height);
            } else {
              glViewport(I->Block->rect.left+I->Width/2,I->Block->rect.bottom,I->Width/2,I->Height);
            }
            break;
          case cStereo_geowall: /* geowall */
            glViewport(I->Block->rect.left+G->Option->winX/2,
                       I->Block->rect.bottom,I->Width,I->Height);
            break;
          case cStereo_stencil_by_row:
          case cStereo_stencil_by_column:
          case cStereo_stencil_checkerboard:
            glStencilFunc(GL_EQUAL, 0, 1);
            glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
            glEnable(GL_STENCIL_TEST);
            break;
          case cStereo_stencil_custom:
#ifdef _PYMOL_SHARP3D
            sharp3d_switch_to_right_stereo();
#endif
            break;
          case cStereo_anaglyph:
            glAccum(GL_ACCUM,0.5);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            break;
          case cStereo_clone_dynamic:
            if(times) {
              glAccum(GL_ACCUM,-0.5);
            } else {
              glAccum(GL_ACCUM,0.5);
            }
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            break;
          case cStereo_dynamic:
            if(times) {
              glAccum(GL_ACCUM,-0.5);
            } else {
              glAccum(GL_ACCUM,0.5);
              glEnable(GL_SCISSOR_TEST);
            }
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
            if(!times) {
              glDisable(GL_SCISSOR_TEST);
            }
            break;
          }

          /* prepare the stereo transformation matrix */
        
          glPushMatrix(); /* 1 */
          ScenePrepareMatrix(G,stereo_using_mono_matrix ? 0 : 2);
          glClear(GL_DEPTH_BUFFER_BIT) ;

          if(grid.active)
            GridGetGLViewport(&grid);
            
          {
            int slot;

            for(slot=0;slot<=grid.last_slot;slot++) {
            
              if(grid.active) { 
                GridSetGLViewport(&grid,slot);
              }

              /* render picked atoms */
            
              glPushMatrix(); /* 2 */
              EditorRender(G,curState);
              glPopMatrix(); /* 1 */
            
              /* render the debugging CGO */
            
              glPushMatrix();  /* 2 */
              glNormal3fv(normal);
              CGORenderGL(G->DebugCGO,NULL,NULL,NULL,NULL);
              glPopMatrix();  /* 1 */
            
              /* render all objects */
            
              glPushMatrix(); /* 2 */
              for(pass=1;pass>-2;pass--) { /* render opaque, then antialiased, then transparent...*/
                SceneRenderAll(G,&context,normal,NULL,pass,false,width_scale,&grid,times);
              }        
              glPopMatrix(); /* 1 */
            
              /* render selections */
              glPushMatrix(); /* 2 */
              glNormal3fv(normal);
              ExecutiveRenderSelections(G,curState);
              glPopMatrix(); /* 1 */
            
            }
          }

          if(grid.active)
            GridSetGLViewport(&grid,-1);
          
          glPopMatrix(); /* 0 */
            
          /* restore draw buffer */

          if(mono_as_quad_stereo) { /* double pumped mono...can't draw to GL_BACK so stick with LEFT */
            OrthoDrawBuffer(G,GL_BACK_LEFT);
          } else switch(stereo_mode) {
          case cStereo_quadbuffer:
            OrthoDrawBuffer(G,GL_BACK_LEFT); /* leave us in a stereo context 
                                                (avoids problems with cards than can't handle
                                                use of mono contexts) */
            break;
          case cStereo_crosseye: 
          case cStereo_walleye:
          case cStereo_sidebyside:
            OrthoDrawBuffer(G,GL_BACK);
            break;
          case cStereo_geowall:
            break;
          case cStereo_stencil_by_row:
          case cStereo_stencil_by_column:
          case cStereo_stencil_checkerboard:
            glDisable(GL_STENCIL_TEST);
            break;
          case cStereo_stencil_custom:
#ifdef _PYMOL_SHARP3D
            sharp3d_end_stereo();
#endif
            break;
          case cStereo_anaglyph:
            glAccum(GL_ACCUM, 0.5);
            glAccum(GL_RETURN, 1.0);
            OrthoDrawBuffer(G,GL_BACK_LEFT);
            break;
          case cStereo_clone_dynamic:
            glAccum(GL_ACCUM, 0.5);
            if(times) {
              float vv[4] = {0.0F,0.0F,0.0F,0.0F};
              glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,vv);
              if(fog_active)
                glEnable(GL_FOG);
              glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
              OrthoDrawBuffer(G,GL_BACK_RIGHT);
            }
            glAccum(GL_RETURN, 1.0);
            OrthoDrawBuffer(G,GL_BACK_LEFT);
            break;
          case cStereo_dynamic:
            glAccum(GL_ACCUM, 0.5);
            if(times) {
              float vv[4] = {0.0F,0.0F,0.0F,0.0F};
              glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,vv);
              if(fog_active)
                glEnable(GL_FOG);
              glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            }
            glAccum(GL_RETURN, 1.0);
            if(times) {
              glViewport(I->Block->rect.left,
                         I->Block->rect.bottom,I->Width+2,I->Height+2);
              glScissor(I->Block->rect.left-1,
                         I->Block->rect.bottom-1,I->Width+2,I->Height+2);
              glEnable(GL_SCISSOR_TEST);
              glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
              glDisable(GL_SCISSOR_TEST);
            } else {
              glDisable(GL_SCISSOR_TEST);
            }
            break;
          }

        } else {

          /* MONOSCOPING RENDERING (not double-pumped) */
          if(grid.active)
            GridGetGLViewport(&grid);
        
          if(Feedback(G,FB_OpenGL,FB_Debugging))
            PyMOLCheckOpenGLErr("Before mono rendering");

          {
            int slot;
            for(slot=0;slot<=grid.last_slot;slot++) {
            
              if(grid.active) { 
                GridSetGLViewport(&grid,slot);
              } 

              /* mono rendering */
            
              PRINTFD(G,FB_Scene)
                " SceneRender: rendering DebugCGO...\n"
                ENDFD;
            
              glPushMatrix();
              glNormal3fv(normal);
              CGORenderGL(G->DebugCGO,NULL,NULL,NULL,NULL);
              glPopMatrix();
            
              glPushMatrix();
              PRINTFD(G,FB_Scene)
                " SceneRender: rendering picked atoms...\n"
                ENDFD;
              EditorRender(G,curState);
              glPopMatrix();
            
              PRINTFD(G,FB_Scene)
                " SceneRender: rendering opaque and antialiased...\n"
                ENDFD;
            
              for(pass=1;pass>-2;pass--) { /* render opaque then antialiased...*/
                SceneRenderAll(G,&context,normal,NULL,pass,false, width_scale, &grid,0);
              }
            
              glPushMatrix();
              PRINTFD(G,FB_Scene)
                " SceneRender: rendering selections...\n"
                ENDFD;
            
              glNormal3fv(normal);
              ExecutiveRenderSelections(G,curState);
            
              PRINTFD(G,FB_Scene)
                " SceneRender: rendering transparent objects...\n"
                ENDFD;
            
              /* render transparent */
              SceneRenderAll(G,&context,normal,NULL,-1,false, width_scale, &grid,0);
              glPopMatrix();
            }
          }
          if(grid.active) {
            GridSetGLViewport(&grid,-1);
          }
          if(Feedback(G,FB_OpenGL,FB_Debugging))
            PyMOLCheckOpenGLErr("during mono rendering");
        }
      }
    }
    
    if(!(pick||smp)) {
      glDisable(GL_FOG);
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT0);
      glDisable(GL_LIGHT1);
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_DITHER);
    } 
    glLineWidth(1.0);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_ALPHA_TEST);
    if(G->Option->multisample)    
      glDisable(0x809D); /* GL_MULTISAMPLE_ARB */
    glViewport(view_save[0],view_save[1],view_save[2],view_save[3]);

    if(Feedback(G,FB_OpenGL,FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender final checkpoint");

  }

  PRINTFD(G,FB_Scene)
    " SceneRender: rendering complete.\n"
    ENDFD;
  
  if(!(pick||smp)) { /* update frames per second field */
    I->RenderTime = -I->LastRender;
    I->LastRender = UtilGetSeconds(G);
    I->RenderTime += I->LastRender;
    I->ApproxRenderTime = I->LastRender - start_time;

    if(I->CopyNextFlag) {
      start_time = I->LastRender - start_time;
      if((start_time>0.10)||(MainSavingUnderWhileIdle()))
        if(!(ControlIdling(G)))
          if(SettingGet(G,cSetting_cache_display)) {
            if(!I->CopyType) {
              SceneCopy(G,render_buffer,false,false);
            }
          }
    } else {
      I->CopyNextFlag=true;
    }
    if(force_copy && !(I->CopyType)) {
      SceneCopy(G,render_buffer,true,false);      
      I->CopyType = 2; /* do not display force copies */
    }
  }
  PRINTFD(G,FB_Scene)
    " SceneRender: leaving...\n"
    ENDFD;
}
/*========================================================================*/
void SceneRestartFrameTimer(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  I->LastFrameTime = UtilGetSeconds(G);
}
static void SceneRestartPerfTimer(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  I->LastRender = UtilGetSeconds(G);
  I->RenderTime = 0.0;
}
void SceneRestartSweepTimer(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  I->LastSweep = 0.0F; /* continue to defer rocking until this is done */
  I->LastSweepX = 0.0F;
  I->LastSweepY = 0.0F;
  I->SweepTime = 0.0;
  SceneRestartPerfTimer(G);

}
/*========================================================================*/
void ScenePrepareMatrix(PyMOLGlobals *G,int mode)
{
  register CScene *I=G->Scene;

  float stAng,stShift;
  
  /* start afresh, looking in the negative Z direction (0,0,-1) from (0,0,0) */
  glLoadIdentity();

  if(!mode) {

    /* mono */

    /* move the camera to the location we are looking at */
    glTranslated(I->Pos[0],I->Pos[1],I->Pos[2]);
  
    /* rotate about the origin (the the center of rotation) */
    glMultMatrixf(I->RotMatrix);            
  
    /* move the origin to the center of rotation */
    glTranslatef(-I->Origin[0],-I->Origin[1],-I->Origin[2]);

  } else { 

    /* stereo */

    stAng = SettingGet(G,cSetting_stereo_angle);
    stShift = SettingGet(G,cSetting_stereo_shift);

    /* right hand */

    stShift = (float)(stShift*fabs(I->Pos[2])/100.0);
    stAng = (float)(stAng*atan(stShift/fabs(I->Pos[2]))*90.0/cPI);
    
    if(mode==2) { /* left hand */
      stAng=-stAng;
      stShift=-stShift;
    }

    PRINTFD(G,FB_Scene)
      " StereoMatrix-Debug: mode %d stAng %8.3f stShift %8.3f \n",mode,stAng,stShift
      ENDFD;
    
    glRotatef(stAng,0.0,1.0,0.0);
    glTranslatef(I->Pos[0],I->Pos[1],I->Pos[2]);
    glTranslatef(stShift,0.0,0.0);

    /* rotate about the origin (the center of rotation) */
    glMultMatrixf(I->RotMatrix);            

    /* move the origin to the center of rotation */
    glTranslatef(-I->Origin[0],-I->Origin[1],-I->Origin[2]);
  }

}
/*========================================================================*/
static void SceneRotateWithDirty(PyMOLGlobals *G,float angle,float x,float y,float z,int dirty)
{
  register CScene *I=G->Scene;
  float temp[16];
  int a;
  angle = (float)(-PI*angle/180.0);
  identity44f(temp);
  MatrixRotateC44f(temp,angle,x,y,z);
  MatrixMultiplyC44f(I->RotMatrix,temp);
  for(a=0;a<16;a++) 
    I->RotMatrix[a]=temp[a];
  SceneUpdateInvMatrix(G);
  if(dirty) {
    SceneInvalidate(G);
  } else {
    SceneInvalidateCopy(G,false);
  }

    /*  glPushMatrix();
        glLoadIdentity();
        glRotatef(angle,x,y,z);
        glMultMatrixf(I->RotMatrix);
        glGetFloatv(GL_MODELVIEW_MATRIX,I->RotMatrix);
        glPopMatrix();*/
}

void SceneRotate(PyMOLGlobals *G,float angle,float x,float y,float z)
{
  SceneRotateWithDirty(G,angle,x,y,z,true);
}
/*========================================================================*/
void SceneApplyMatrix(PyMOLGlobals *G,float *m)
{
  register CScene *I=G->Scene;
  MatrixMultiplyC44f(m,I->RotMatrix);
  SceneDirty(G);
  
  /*  glPushMatrix();
      glLoadIdentity();
      glMultMatrixf(m);
      glMultMatrixf(I->RotMatrix);
      glGetFloatv(GL_MODELVIEW_MATRIX,I->RotMatrix);
  glPopMatrix();*/
}
/*========================================================================*/
void SceneScale(PyMOLGlobals *G,float scale)
{
  register CScene *I=G->Scene;
  I->Scale*=scale;
  SceneInvalidate(G);
}










