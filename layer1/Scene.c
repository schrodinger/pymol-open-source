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
  struct CObject *obj;  
  struct ObjRec *next;
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

/* allow up to 10 seconds at 30 FPS */

#define MAX_ANI_ELEM 300

struct _CScene {
  Block *Block;
  ObjRec *Obj;
  float RotMatrix[16];
  float InvMatrix[16];
  float ModMatrix[16];
  float ProMatrix[16];
  float Scale;
  int Width,Height;
  int Button;
  int LastX,LastY;
  int StartX,StartY;
  int LastWinX,LastWinY;
  double LastClickTime;
  int LastButton;
  int PossibleSingleClick;
  double LastReleaseTime;
  double SingleClickDelay;
  float ViewNormal[3],LinesNormal[3];
  float Pos[3],Origin[3];
  float H;
  float Front,Back,FrontSafe,BackSafe;
  float TextColor[3];
  double RockTime;
  int DirtyFlag;
  int ChangedFlag;
  int CopyFlag,CopyNextFlag,CopiedFromOpenGL;
  int NFrame;
  GLvoid *ImageBuffer;
  int ImageBufferHeight,ImageBufferWidth;
  int MovieOwnsImageFlag;
  int MovieFrameFlag;
  unsigned ImageBufferSize;
  double LastRender,RenderTime,LastFrameTime;
  double LastRock,LastRockTime;
  float LastRockX,LastRockY;
  Pickable LastPicked;
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

};

typedef struct {
  float unit_left,unit_right,unit_top,unit_bottom,unit_front,unit_back;
} SceneUnitContext;

static int SceneDeferClickWhen(Block *block, int button, int x, int y, double when);
int SceneLoopDrag(Block *block,int x,int y,int mod);
int SceneLoopRelease(Block *block,int button,int x,int y,int mod);

int SceneLoopClick(Block *block,int button, int x,int y,int mod);

void ScenePrimeAnimation(PyMOLGlobals *G)
{
  if(G->HaveGUI) {
    CScene *I=G->Scene;
    UtilZeroMem(I->ani_elem,sizeof(CViewElem));
    SceneToViewElem(G,I->ani_elem);
    I->ani_elem[0].specification_level = 2;
    I->n_ani_elem = 0;
  }
}

void SceneLoadAnimation(PyMOLGlobals *G, double duration)
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
    SceneToViewElem(G,I->ani_elem + target);
    I->ani_elem[target].specification_level = 2;
    now = UtilGetSeconds(G);
    I->ani_elem[0].timing_flag = true;
    I->ani_elem[0].timing = now + 0.01;
    I->ani_elem[target].timing_flag = true;
    I->ani_elem[target].timing = now + duration;
    ViewElemInterpolate(I->ani_elem, I->ani_elem + target, 2.0F, 1.0F);
    SceneFromViewElem(G,I->ani_elem);
    I->cur_ani_elem = 0;
    I->n_ani_elem = target;
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
  UtilZeroMem(im,sizeof(float)*16);
  im[0] = rm[0];
  im[1] = rm[4];
  im[2] = rm[8];
  im[4] = rm[1];
  im[5] = rm[5];
  im[6] = rm[9];
  im[8] = rm[2];
  im[9] = rm[6];
  im[10] = rm[10];
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

void SceneCopy(PyMOLGlobals *G,GLenum buffer,int force);

unsigned int SceneFindTriplet(PyMOLGlobals *G,int x,int y,GLenum gl_buffer);
unsigned int *SceneReadTriplets(PyMOLGlobals *G,int x,int y,int w,int h,GLenum gl_buffer);

void SceneDraw(Block *block);
void ScenePrepareMatrix(PyMOLGlobals *G,int mode);

void ScenePrepareUnitContext(PyMOLGlobals *G,SceneUnitContext *context,int width,int height);

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

void SceneToViewElem(PyMOLGlobals *G,CViewElem *elem)
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
  elem->ortho = SettingGetGlobal_b(G,cSetting_ortho);
  
}

void SceneFromViewElem(PyMOLGlobals *G,CViewElem *elem)
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
    SceneClipSet(G,elem->front,elem->back);
  }
  if(elem->ortho_flag) {
    SettingSetGlobal_b(G,cSetting_ortho,elem->ortho);
  }
  if(changed_flag) {

    I->LastRock = 0.0F; /* continue to defer rocking until this is done */
    I->LastRockX = 0.0F;
    I->LastRockY = 0.0F;
    I->RockTime = 0.0;

    SceneRovingDirty(G);
  }
}

void SceneCleanupStereo(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;  
  if(I->StereoMode==1)
    PSGIStereo(0);
}

void ScenePrepareUnitContext(PyMOLGlobals *G,SceneUnitContext *context,int width,int height)
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

  PRINTFD(G,FB_Scene)
    "ScenePrepareUnitContext:%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
    context->unit_left,
    context->unit_right,
    context->unit_top, 
    context->unit_bottom,
    context->unit_front,
    context->unit_back
    ENDFD;

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

  MatrixTransform3f(I->RotMatrix,I->Origin,pos); 

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

  MatrixInvTransform3f(I->RotMatrix,pos,pos);

  PRINTFD(G,FB_Scene)
    " SceneGetPos: center            "
    ENDFD3f(pos);

}
/*========================================================================*/
void SceneApplyRotMatrix(PyMOLGlobals *G,float *src,float *dst)
{
  register CScene *I=G->Scene;
  MatrixTransform3f(I->RotMatrix,src,dst);
}
/*========================================================================*/
int SceneMultipick(PyMOLGlobals *G,Multipick *smp)
{
  register CScene *I=G->Scene;
  if(((int)SettingGet(G,cSetting_overlay))&&((int)SettingGet(G,cSetting_text)))
    SceneRender(G,NULL,0,0,NULL); /* remove overlay if present */
  SceneDontCopyNext(G);
  if(I->StereoMode>1) {
    smp->x = smp->x % (I->Width/2);
  }
  SceneRender(G,NULL,0,0,smp);
  SceneDirty(G);
  return(1);
}
/*========================================================================*/
int SceneGetNFrame(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
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
  *(p++) = SettingGet(G,cSetting_ortho);
}
/*========================================================================*/
void SceneSetView(PyMOLGlobals *G,SceneViewType view,int quiet,int animate)
{
  float *p;
  int a;
  register CScene *I=G->Scene;

  if(animate<0)
    animate=SettingGetGlobal_b(G,cSetting_animation);
  if(animate)
    ScenePrimeAnimation(G);

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

  I->LastRock = 0.0F;
  I->LastRockX = 0.0F;
  I->LastRockY = 0.0F;
  I->RockTime = 0.0;

  SceneClipSet(G,p[0],p[1]);
  p+=2;
  SettingSet(G,cSetting_ortho,*(p++));
  if(!quiet) { 
    PRINTFB(G,FB_Scene,FB_Actions)
      " Scene: view updated.\n"
      ENDFB(G);
  }
  if(animate)
    SceneLoadAnimation(G,SettingGetGlobal_f(G,cSetting_animation_duration));

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
  if(flag) 
    I->StereoMode=(int)SettingGet(G,cSetting_stereo_mode);
  else
    I->StereoMode=false;
  SettingSetGlobal_b(G,cSetting_stereo,flag);
  SceneDirty(G);
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
  SceneDirty(G);
}
/*========================================================================*/
void SceneClipSet(PyMOLGlobals *G,float front,float back)
{
  register CScene *I=G->Scene;
  I->Front=front;
  I->Back=back;
  if(I->Front>I->Back)
	 I->Front=I->Back+cSliceMin;
  I->FrontSafe= GetFrontSafe(I->Front,I->Back);
  I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
  SceneDirty(G);
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
        MatrixTransform3f(I->RotMatrix,v0,offset); /* convert to view-space */
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
    if(!ExecutiveGetCameraExtent(G,sele,mn,mx,true,state))
      sele = NULL;
    if(sele) {
      if(sele[0]) {
        average3f(mn,mx,cent); /* get center of selection */
        MatrixTransform3f(I->RotMatrix,I->Origin,origin); /* convert to view-space */
        subtract3f(mx,origin,mx); /* how far from origin? */
        subtract3f(mn,origin,mn); /* how far from origin? */
        SceneClipSet(G,-I->Pos[2]-mx[2]-movement,-I->Pos[2]-mn[2]+movement);
      } else {
        sele = NULL;
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
static unsigned char *SceneImagePrepare(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  unsigned int buffer_size;
  GLvoid *image;
  int reset_alpha = false;

  buffer_size = 4*I->Width*I->Height;
  if(!I->CopyFlag) {
	 image = (GLvoid*)Alloc(char,buffer_size);
	 ErrChkPtr(G,image);
    if(G->HaveGUI && G->ValidContext) {
      glReadBuffer(GL_BACK);
      PyMOLReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                   GL_RGBA,GL_UNSIGNED_BYTE,image);
      
      reset_alpha = true;
      I->ImageBufferHeight=I->Height;
      I->ImageBufferWidth=I->Width;
    } else {
       PRINTFB(G,FB_Scene,FB_Errors)
         " ScenePNG-WARNING: writing a blank image buffer.\n"
         ENDFB(G);
     }
  } else {
    image=I->ImageBuffer;
    reset_alpha = I->CopiedFromOpenGL;
    PRINTFB(G,FB_Scene,FB_Blather)
      " ScenePNG: writing cached image (reset_alpha=%d).\n",reset_alpha
      ENDFB(G);
  }
  if(reset_alpha&&image) {
    unsigned char *p = (unsigned char*)image;
    int x,y;
    for(y=0;y<I->Height;y++) {
      for(x=0;x<I->Width;x++) {
        p[3]=0xFF;
        p+=4;
      }
    }
  }
  return (unsigned char*)image;
}

static void SceneImageFinish(PyMOLGlobals *G,char *image)
{
  register CScene *I=G->Scene;

  if(!I->CopyFlag)
	 FreeP(image);
}

int  SceneCopyExternal(PyMOLGlobals *G,int width, int height,int rowbytes,unsigned char *dest)
{
  GLvoid *image = SceneImagePrepare(G);
  register CScene *I=G->Scene;
  int result=false;
  int i,j;
  if(image&&(I->ImageBufferWidth==width)&&(I->ImageBufferHeight==height)) {
    for (i=0; i< height; i++)
      {
        unsigned char *dst = dest + i * (rowbytes);
        unsigned char *src = ((unsigned char*)image) + ((height-1)-i) * width*4;
        for (j = 0; j < width; j++)
          {
            *dst++ = ((unsigned int)src[0]*src[3])/255; /* premultiply alpha */
            *dst++ = ((unsigned int)src[1]*src[3])/255;
            *dst++ = ((unsigned int)src[2]*src[3])/255;
            *dst++ = src[3];
            src+=4;
          }
      }
    result=true;
  }
  SceneImageFinish(G,image);  
  return(result);
}

void ScenePNG(PyMOLGlobals *G,char *png,int quiet)
{
  register CScene *I=G->Scene;
  GLvoid *image = SceneImagePrepare(G);
  if(image) {
    if(MyPNGWrite(G,png,image,I->ImageBufferWidth,I->ImageBufferHeight)) {
      if(!quiet) {
        PRINTFB(G,FB_Scene,FB_Actions) 
          " ScenePNG: wrote %dx%d pixel image to file \"%s\".\n",
          I->ImageBufferWidth,I->ImageBufferHeight,png
          ENDFB(G);
      }
    } else {
      PRINTFB(G,FB_Scene,FB_Errors) 
        " ScenePNG-Error: error writing \"%s\"! Please check directory...\n",
        png
        ENDFB(G);
    }
  }
  SceneImageFinish(G,image);  
}
/*========================================================================*/
void ScenePerspective(PyMOLGlobals *G,int flag)
{
  float persp;
  persp=(float)(!flag);
  SettingSetfv(G,cSetting_ortho,&persp);
  SceneDirty(G);
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
  while(ListIterate(I->Obj,rec,next))
	 {
      if(rec->obj->fGetNFrame)
        n=rec->obj->fGetNFrame(rec->obj);
      else
        n=0;
		if(n>I->NFrame)
		  I->NFrame=n;
	 }
  mov_len = MovieGetLength(G);
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
      MovieMatrix(G,cMovieMatrixRecall);
    }
    if(movieCommand) {
      MovieDoFrameCommand(G,newFrame);
    }
    if(SettingGet(G,cSetting_cache_frames))
      I->MovieFrameFlag=true;
  }
  SettingSetGlobal_i(G,cSetting_frame,newFrame+1);
  SettingSetGlobal_i(G,cSetting_state,newState+1);

  SceneDirty(G);
  PRINTFD(G,FB_Scene)
    " SceneSetFrame: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void ScenePurgeCopy(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  I->CopyFlag=false;
  if(I->MovieOwnsImageFlag) 
	 {
		I->MovieOwnsImageFlag=false;
		I->ImageBuffer=NULL;
	 }
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
    I->DirtyFlag=true;
    ScenePurgeCopy(G);
    OrthoDirty(G);
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
void SceneMakeMovieImage(PyMOLGlobals *G) {
  register CScene *I=G->Scene;
  float *v;

  PRINTFB(G,FB_Scene,FB_Blather)
    " Scene: Making movie image.\n"
    ENDFB(G);

  I->DirtyFlag=false;
  if(SettingGet(G,cSetting_ray_trace_frames)) {
	SceneRay(G,0,0,(int)SettingGet(G,cSetting_ray_default_renderer),NULL,NULL,
            0.0F,0.0F,false,NULL); 
  } else {
	 v=SettingGetfv(G,cSetting_bg_rgb);
    if(G->HaveGUI && G->ValidContext) {

      glDrawBuffer(GL_BACK);
      glClearColor(v[0],v[1],v[2],1.0);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glClearColor(0.0,0.0,0.0,1.0);
      SceneRender(G,NULL,0,0,NULL);
      SceneCopy(G,GL_BACK,true);
    }
  }
  if(I->ImageBuffer&&(I->ImageBufferHeight==I->Height)&&(I->ImageBufferWidth==I->Width)) {
	 MovieSetImage(G,MovieFrameToImage(G,SettingGetGlobal_i(G,cSetting_frame)-1)
                                    ,I->ImageBuffer);
    I->MovieOwnsImageFlag=true;
  } else {
    I->MovieOwnsImageFlag=false;
  }
  I->CopyFlag=true;
}

/*========================================================================*/
void SceneIdle(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  double renderTime;
  double minTime;
  int frameFlag = false;
  int rockFlag = false;

  if(I->PossibleSingleClick==2) {
    double now = UtilGetSeconds(G);
    double single_click_delay = I->SingleClickDelay;
    double diff = now-I->LastReleaseTime;
    if(diff>single_click_delay) {
      /* post a single click processing event */
      SceneDeferClickWhen(I->Block, 
                          I->LastButton + P_GLUT_SINGLE_LEFT,
                          I->LastWinX, I->LastWinY,
                          I->LastClickTime); /* push a click onto the queue */
      
      I->PossibleSingleClick = 0;
      OrthoDirty(G); /* force an update */
    }
  }
  if(MoviePlaying(G))
    {
		renderTime = -I->LastFrameTime + UtilGetSeconds(G);
		minTime=SettingGet(G,cSetting_movie_delay)/1000.0;
		if(renderTime>=minTime) {
        frameFlag=true;
        rockFlag=true;
      }
    }
  if(ControlRocking(G)&&(!rockFlag))
    {
		renderTime = -I->LastRockTime + UtilGetSeconds(G);
		minTime=SettingGet(G,cSetting_rock_delay)/1000.0;
		if(renderTime>=minTime) {
        rockFlag=true;
        I->LastRockTime=UtilGetSeconds(G);
      }
    }
  if(ControlRocking(G)&&rockFlag) {
    float ang_cur,disp,diff;
    float sweep_angle = SettingGet(G,cSetting_sweep_angle);
    float sweep_speed = SettingGet(G,cSetting_sweep_speed);
    float sweep_phase = SettingGet(G,cSetting_sweep_phase);
    int sweep_mode = SettingGet(G,cSetting_sweep_mode);
    float shift = PI/2.0F;

    I->RockTime+=I->RenderTime;

    switch(sweep_mode) {
    case 0:
    case 1:
    case 2:
      if(sweep_angle<=0.0F) {
        diff = (float)((PI/180.0F)*I->RenderTime*10);
      } else {
        ang_cur = (float)(I->RockTime*sweep_speed) + sweep_phase;
        disp = (float)(sweep_angle*(PI/180.0F)*sin(ang_cur)/2);
        diff = (float)(disp-I->LastRock);
        I->LastRock = disp;
      }
      switch(sweep_mode) {
      case 0:
        SceneRotate(G,(float)(180*diff/PI),0.0F,1.0F,0.0F);
        break;
      case 1:
        SceneRotate(G,(float)(180*diff/PI),1.0F,0.0F,0.0F);
        break;
      case 2: /* z-rotation...useless! */
        SceneRotate(G,(float)(180*diff/PI),0.0F,0.0F,1.0F);
        break;
      }
      break;
    case 3: /* nutate */
      SceneRotate(G,(float)(-I->LastRockY),0.0F,1.0F,0.0F);
      SceneRotate(G,(float)(-I->LastRockX),1.0F,0.0F,0.0F);
      ang_cur = (float)(I->RockTime*sweep_speed) + sweep_phase;
      
      I->LastRockX = (sweep_angle*sin(ang_cur)/2);
      I->LastRockY = (sweep_angle*sin(ang_cur+shift)/2);
      
      if(I->RockTime*sweep_speed<PI) {
        float factor = (I->RockTime*sweep_speed)/PI;
        I->LastRockX *= factor;
        I->LastRockY *= factor;
      }
      SceneRotate(G,(float)I->LastRockX,1.0F,0.0F,0.0F);
      SceneRotate(G,(float)I->LastRockY,0.0F,1.0F,0.0F);
      break;
    }
  }
  if(MoviePlaying(G)&&frameFlag)
	 {
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

  MatrixTransform3f(I->RotMatrix,v0,I->Pos); /* convert to view-space */
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

  MatrixTransform3f(I->RotMatrix,v0,I->Pos); /* convert to view-space */

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
  
  if(preserve) /* preserve current viewing location */
	 {
		subtract3f(origin,I->Origin,v0); /* model-space translation */
		MatrixTransform3f(I->RotMatrix,v0,v1); /* convert to view-space */
		add3f(I->Pos,v1,I->Pos); /* offset view to compensate */
	 }
  I->Origin[0]=origin[0]; /* move origin */
  I->Origin[1]=origin[1];
  I->Origin[2]=origin[2];
  SceneDirty(G);
}
/*========================================================================*/
void SceneObjectAdd(PyMOLGlobals *G,CObject *obj)
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
}
/*========================================================================*/
void SceneObjectDel(PyMOLGlobals *G,CObject *obj)
{
  register CScene *I=G->Scene;
  ObjRec *rec = NULL;

  if(!obj) {
    while(ListIterate(I->Obj,rec,next)) {
      if(rec) {
        ListDetach(I->Obj,rec,next,ObjRec);
        ListElemFree(rec);
      }
    }
  } else {
    while(ListIterate(I->Obj,rec,next))
      if(rec->obj==obj)
        break;
    if(rec) {
      rec->obj->Enabled=false;
      ListDetach(I->Obj,rec,next,ObjRec);
      ListElemFree(rec);
    }
  }
  SceneCountFrames(G);
  SceneDirty(G);
}
/*========================================================================*/
int SceneLoadPNG(PyMOLGlobals *G,char *fname,int movie_flag,int quiet) 
{
  register CScene *I=G->Scene;
  int ok=false;
  if(I->ImageBuffer) {
	 if(I->MovieOwnsImageFlag) {
		I->MovieOwnsImageFlag=false;
		I->ImageBuffer=NULL;
	 } else {
		FreeP(I->ImageBuffer);
	 }
    I->CopyFlag=false;
  }
  if(MyPNGRead(fname,(unsigned char**)&I->ImageBuffer,(unsigned int*)&I->ImageBufferWidth,(unsigned int*)&I->ImageBufferHeight)) {
    if(!quiet) {
      PRINTFB(G,FB_Scene,FB_Details)
        " Scene: loaded image from '%s'.\n",fname
        ENDFB(G);
    }
    I->CopyFlag=true;
    I->CopiedFromOpenGL = false;
    OrthoRemoveSplash(G);
    SettingSet(G,cSetting_text,0.0);
    if(movie_flag&&I->ImageBuffer&&(I->ImageBufferHeight==I->Height)&&(I->ImageBufferWidth==I->Width)) {
      MovieSetImage(G,
                    MovieFrameToImage(G,
                                      SettingGetGlobal_i(G,cSetting_frame)-1)
                    ,I->ImageBuffer);
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
/*========================================================================*/
void SceneDraw(Block *block)
{
  PyMOLGlobals *G=block->G;
  register CScene *I=G->Scene;
  int overlay,text;
  int width,height;
  int double_pump;

  if(G->HaveGUI && G->ValidContext) {
    overlay = (int)SettingGet(G,cSetting_overlay);
    text = (int)SettingGet(G,cSetting_text);
    double_pump = (int)SettingGet(G,cSetting_stereo_double_pump_mono);

    if(overlay||(!text)) 

      if(I->CopyFlag)
        {
          glReadBuffer(GL_BACK); 

          if(I->ImageBufferHeight>I->Height||I->ImageBufferWidth>I->Width) {
            TextSetColor3f(G,1.0F,0.2F,0.2F);
            TextDrawStrAt(G,"Sorry, I can't display an oversize image.",30,60);
            TextDrawStrAt(G,"To save image, use File Menu or enter \"png <filename>\".",30,40);
          } else {
            width = I->ImageBufferWidth;
            height = I->ImageBufferHeight;
            
            if((width<I->Width)||(height<I->Height)) {
              glRasterPos3i((int)((I->Width-width)/2+I->Block->rect.left),
                            (int)((I->Height-height)/2+I->Block->rect.bottom),-10);
            } else {
              glRasterPos3i(I->Block->rect.left,I->Block->rect.bottom,-10);
            }
            if(I->ImageBuffer) {
#if 1
              PyMOLDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
#else
              if(!(double_pump||(I->StereoMode==1))) {
                glDrawBuffer(GL_BACK);
                PyMOLDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
              } else {
                glDrawBuffer(GL_BACK_LEFT);
                PyMOLDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
                glDrawBuffer(GL_BACK_RIGHT);
                PyMOLDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
              }
#endif
            }

          }
          I->RenderTime = -I->LastRender;
          I->LastRender = UtilGetSeconds(G);
          I->RenderTime += I->LastRender;
          ButModeSetRate(G,(float)I->RenderTime);
        }
    
    glColor3f(1.0,1.0,1.0);
  }
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
  ObjectMolecule *obj;
  
  I->LastReleaseTime = when;
 
  if(I->PossibleSingleClick==1) {
    double slowest_single_click = 0.25;
    double diff = when-I->LastClickTime;
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
        int mode=ButModeTranslate(G,but, 0);
        if(mode == cButModeNone)
          I->SingleClickDelay = 0.0;
      }
    }
  }
  if(I->LoopFlag)
    return SceneLoopRelease(block,button,x,y,mod);
  if(I->SculptingFlag) {
    /* SettingSet(G,cSetting_sculpting,1); */
    obj=(ObjectMolecule*)I->LastPicked.ptr;
    if(obj) {
      obj->AtomInfo[I->LastPicked.index].protekted=I->SculptingSave;
    }
    I->SculptingFlag=0;
  }
  I->LastPickVertexFlag=false;
  return 1;
}
/*========================================================================*/
static void SceneDoRoving(PyMOLGlobals *G,float old_front,
                          float old_back,float old_origin,
                          int adjust_flag,int zoom_flag)
{
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

    MatrixInvTransform3f(I->RotMatrix,v2,v2); /* transform offset into realspace */
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

static int SceneDoXYPick(PyMOLGlobals *G, int x, int y)
{
  CScene *I=G->Scene;
  if(((int)SettingGet(G,cSetting_overlay))&&((int)SettingGet(G,cSetting_text)))
    SceneRender(G,NULL,0,0,NULL); /* remove overlay if present */
  SceneDontCopyNext(G);
  
  I->LastPicked.ptr = NULL;
  SceneRender(G,&I->LastPicked,x,y,NULL);
  return (I->LastPicked.ptr!=NULL);
  /* did we pick something? */
}

static void SceneNoteMouseInteraction(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  if(I->cur_ani_elem < I->n_ani_elem ) { /* allow user to override animation */
    I->cur_ani_elem = I->n_ani_elem;
  }
  if(SettingGet_b(G,NULL,NULL,cSetting_mouse_restart_movie_delay)) {
    SceneRestartTimers(G);
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

  if(!is_single_click) {
    if((!(mod&(cOrthoCTRL+cOrthoSHIFT)))&&((when-I->LastClickTime)<cDoubleTime))
      {
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
    
    if(!mod)
      I->PossibleSingleClick = 1;
    else
      I->PossibleSingleClick = 0;
  }

  I->LastWinX = x;
  I->LastWinY = y;
  I->LastClickTime = when;
  I->LastButton = button;
  I->Threshold = 0;

  mode = ButModeTranslate(G,button,mod); 

  I->Button=button;    
  I->SculptingSave = 0;
  switch(mode) {
  case cButModeScaleSlabExpand:
    SceneNoteMouseInteraction(G);
    SceneClip(G,5,1.2F,NULL,0);
    break;
  case cButModeScaleSlabShrink:
    SceneNoteMouseInteraction(G);
    SceneClip(G,5,0.8F,NULL,0);
    break;
  case cButModeMoveSlabForward:
    SceneNoteMouseInteraction(G);
    {
      float old_front = I->Front;
      float old_back = I->Back;
      float old_origin = -I->Pos[2];
      SceneClip(G,6,0.1F,NULL,0);
      SceneDoRoving(G,old_front,old_back,old_origin,true,false);
    }
    break;
  case cButModeMoveSlabBackward:
    SceneNoteMouseInteraction(G);
    {
      float old_front = I->Front;
      float old_back = I->Back;
      float old_origin = -I->Pos[2];

      SceneClip(G,6,-0.1F,NULL,0);
      SceneDoRoving(G,old_front,old_back,old_origin,true,false);
    }
    break;
  case cButModeZoomForward:
    SceneNoteMouseInteraction(G);
    {
      float factor = -((I->FrontSafe+I->BackSafe)/2)*0.1F;
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
      float factor = ((I->FrontSafe+I->BackSafe)/2)*0.1F;
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
      SceneClip(G,6,0.1F,NULL,0);
      SceneDoRoving(G,old_front,old_back,old_origin,true,true);
    }
    break;
  case cButModeMoveSlabAndZoomBackward:
    SceneNoteMouseInteraction(G);
    {
      float old_front = I->Front;
      float old_back = I->Back;
      float old_origin = -I->Pos[2];
      SceneClip(G,6,-0.1F,NULL,0);
      SceneDoRoving(G,old_front,old_back,old_origin,true,true);
    }
    break;
  case cButModeRectAdd: /* deprecated */
  case cButModeRectSub:/* deprecated */
  case cButModeRect:/* deprecated */
  case cButModeSeleAdd:
  case cButModeSeleSub:
    return SceneLoopClick(block,button,x,y,mod);
    break;
  case cButModeRotXYZ:
  case cButModeTransXY:
  case cButModeTransZ:
  case cButModeClipNF:
  case cButModeClipN:    
  case cButModeClipF:    
  case cButModeRotZ:
    SceneNoteMouseInteraction(G);
    SceneDontCopyNext(G);

    y=y-I->Block->margin.bottom;
    x=x-I->Block->margin.left;
    
    if(I->StereoMode>1)
      x = x % (I->Width/2);

    I->LastX=x;
    I->LastY=y;	 

    /*    SceneDirty(G);*/
    break;
  case cButModePickAtom1:
  case cButModePickAtom:
  case cButModeMenu:
    if(I->StereoMode>1)
      x = x % (I->Width/2);

    if(SceneDoXYPick(G,x,y)) {
      obj=(CObject*)I->LastPicked.ptr;
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
            if(active_sele && SelectorIsMember(G,objMol->AtomInfo[I->LastPicked.index].selEntry,
                                               active_sele)) {
              char name[ObjNameMax];
              ExecutiveGetActiveSeleName(G,name,false);
              MenuActivate2Arg(G,I->LastWinX,I->LastWinY+20, /* selection menu */
                               I->LastWinX,I->LastWinY,
                               is_single_click,
                               "pick_sele",name,name);
            } else {
              ObjectMoleculeGetAtomSele((ObjectMolecule*)obj,I->LastPicked.index,buffer);
              ObjectMoleculeGetAtomSeleLog((ObjectMolecule*)obj,I->LastPicked.index,buf1,false);
              MenuActivate2Arg(G,I->LastWinX,I->LastWinY+20,
                               I->LastWinX,I->LastWinY,
                               is_single_click,
                               "pick_menu",buffer,buf1);
            }
          }
          break;
        case cButModePickAtom1:
          if(obj&&obj->type==cObjectMolecule) {
            if(Feedback(G,FB_ObjectMolecule,FB_Results)) {
              if(obj->fDescribeElement)
                obj->fDescribeElement(obj,I->LastPicked.index,buffer);
              PRINTF " You clicked %s -> (%s)\n",buffer,cEditorSele1 ENDF(G);
            }
            if(SettingGet(G,cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buffer,false);
              sprintf(buf2,"cmd.edit(\"%s\",pkresi=1)",buffer);
              PLog(buf2,cPLog_pym);
            }
            OrthoRestorePrompt(G);
            sprintf(buffer,"%s`%d",
                    obj->Name,I->LastPicked.index+1);    
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
              obj->fDescribeElement(obj,I->LastPicked.index,buffer);
            if(EditorIsBondMode(G)
               /* &&!(EditorIsAnActiveObject(G,(ObjectMolecule*)obj))*/ ) {
              EditorInactivate(G);
              EditorLogState(G,false);
            }
            if((!EditorIsBondMode(G))&&
               EditorDeselectIfSelected(G,
                                        (ObjectMolecule*)obj,I->LastPicked.index,true)) {
              PRINTF " You unpicked %s.",buffer ENDF(G);
              if(EditorActive(G)) 
                EditorDefineExtraPks(G);
              EditorLogState(G,false);
            } else {
              if(EditorIsBondMode(G)&&
                 EditorDeselectIfSelected(G,
                                          (ObjectMolecule*)obj,I->LastPicked.index,false)) {
                EditorInactivate(G);
              }
              EditorGetNextMultiatom(G,name);

              PRINTF " You clicked %s -> (%s)\n",buffer,name ENDF(G);
              /* TODO: logging */
              
              sprintf(buffer,"%s`%d",obj->Name,I->LastPicked.index+1);    
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
        break;
      }
    }
    SceneDirty(G);
    break;
  case cButModePickBond:
  case cButModePkTorBnd:
    if(I->StereoMode>1)
      x = x % (I->Width/2);

    if(SceneDoXYPick(G,x,y)) {
      obj=(CObject*)I->LastPicked.ptr;
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
        if(Feedback(G,FB_ObjectMolecule,FB_Results)) {
          if(obj->fDescribeElement)
            obj->fDescribeElement(obj,I->LastPicked.index,buffer);
          PRINTF " You clicked %s -> (%s)",buffer,cEditorSele1 ENDF(G);
          OrthoRestorePrompt(G);
        }
	sprintf(buffer,"%s`%d",
		obj->Name,I->LastPicked.index+1);    
        SelectorCreate(G,cEditorSele1,buffer,NULL,true,NULL);
        objMol = (ObjectMolecule*)obj;
        if(I->LastPicked.bond>=0) {
          atIndex = objMol->Bond[I->LastPicked.bond].index[0];
          if(atIndex == I->LastPicked.index)
            atIndex = objMol->Bond[I->LastPicked.bond].index[1];              
          if(Feedback(G,FB_ObjectMolecule,FB_Results)) {
            if(obj->fDescribeElement)
              obj->fDescribeElement(obj,atIndex,buffer);
            PRINTF " You clicked %s -> (%s)",buffer,cEditorSele2 ENDF(G);
            OrthoRestorePrompt(G);
          }
          
          if(SettingGet(G,cSetting_logging)) {
            objMol = (ObjectMolecule*)obj;            
            ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1,false);
            ObjectMoleculeGetAtomSeleLog(objMol,atIndex,buf2,false);
            sprintf(buffer,"cmd.edit(\"%s\",\"%s\")",buf1,buf2);
            PLog(buffer,cPLog_pym);
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
              EditorPrepareDrag(G,objMol,I->LastPicked.index,
                                SettingGetGlobal_i(G,cSetting_state)-1);
              I->SculptingFlag = 1;
              I->SculptingSave =  objMol->AtomInfo[I->LastPicked.index].protekted;
              objMol->AtomInfo[I->LastPicked.index].protekted=2;
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
    }
    SceneDirty(G);
    break;
  case cButModeMovFrag:
  case cButModeTorFrag:
  case cButModeRotFrag:
  case cButModeMoveAtom:
    if(I->StereoMode>1)
      x = x % (I->Width/2);
    if(SceneDoXYPick(G,x,y)) {
      obj=(CObject*)I->LastPicked.ptr;
      y=y-I->Block->margin.bottom;
      x=x-I->Block->margin.left;
      I->LastX=x;
      I->LastY=y;	
      switch(obj->type) {
      case cObjectMolecule:
        
        if(Feedback(G,FB_ObjectMolecule,FB_Results)) {
          if(obj->fDescribeElement) 
            obj->fDescribeElement(obj,I->LastPicked.index,buffer);
          PRINTF " You clicked %s",buffer ENDF(G);        
          OrthoRestorePrompt(G);
        }
        objMol = (ObjectMolecule*)obj;
        EditorPrepareDrag(G,objMol,I->LastPicked.index,
                          SettingGetGlobal_i(G,cSetting_state)-1);
        I->SculptingFlag = 1;
        I->SculptingSave =  objMol->AtomInfo[I->LastPicked.index].protekted;
        objMol->AtomInfo[I->LastPicked.index].protekted=2;
        break;
      case cObjectSlice:
        
        if(ObjectSliceGetVertex((ObjectSlice*)obj,I->LastPicked.index,I->LastPicked.bond,I->LastPickVertex)) {
          I->LastPickVertexFlag=true;
        }
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

  case cButModeOrigAt:
  case cButModeCent:
    if(I->StereoMode>1)
      x = x % (I->Width/2);

    if(SceneDoXYPick(G,x,y)) {
      obj=(CObject*)I->LastPicked.ptr;

      switch(obj->type) {
      case cObjectMolecule:
        if(Feedback(G,FB_ObjectMolecule,FB_Results)) {
          if(obj->fDescribeElement) 
            obj->fDescribeElement(obj,I->LastPicked.index,buffer);
          PRINTF " You clicked %s",buffer ENDF(G);        
          OrthoRestorePrompt(G);
        }
        sprintf(buffer,"%s`%d",
                obj->Name,I->LastPicked.index+1);
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
          ExecutiveGetActiveSeleName(G,selName,true);
          break;
          
        case cButModeOrigAt:
          SceneNoteMouseInteraction(G);
#if 1
          {
            float v1[3];

            if(ObjectMoleculeGetAtomVertex((ObjectMolecule*)obj,
                                           SettingGetGlobal_i(G,cSetting_state)-1,
                                           I->LastPicked.index,v1)) {
              ExecutiveOrigin(G,NULL,true,NULL,v1,0);
            }
          }
#else
          sprintf(buf2,"origin (%s)",buffer);        
          OrthoCommandIn(G,buf2);
#endif
          if(obj->type==cObjectMolecule) {
            if(SettingGet(G,cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1,false);
              sprintf(buffer,"cmd.origin(\"%s\")",buf1);
              PLog(buffer,cPLog_pym);

            }
            if(Feedback(G,FB_ObjectMolecule,FB_Results)) {
              if(obj->fDescribeElement) 
                obj->fDescribeElement(obj,I->LastPicked.index,buffer);
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
#if 1
          {
            float v1[3];

            if(ObjectMoleculeGetAtomVertex((ObjectMolecule*)obj,
                                           SettingGetGlobal_i(G,cSetting_state)-1,
                                           I->LastPicked.index,v1)) {
              ExecutiveCenter(G,NULL,0,true,-1,v1);
            }
          }
          
#else
          sprintf(buf2,"center (%s),state=-1,animate=-1",buffer);        
          OrthoCommandIn(G,buf2);
#endif
          if(SettingGet(G,cSetting_logging)) {
            objMol = (ObjectMolecule*)obj;            
            ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1,false);
            sprintf(buffer,"cmd.center(\"%s\",state=-1)",buf1);
            PLog(buffer,cPLog_pym);
          }
          break;
        }
        switch(mode) {
        case cButModeLB:
        case cButModeMB:
        case cButModeRB:
        case cButModeSeleSet:
	  sprintf(buf2,"(%s(%s))",sel_mode_kw,buffer);
	  SelectorCreate(G,selName,buf2,NULL,false,NULL);
          if(SettingGet(G,cSetting_auto_hide_selections))
            ExecutiveHideSelections(G);
          if(SettingGet(G,cSetting_auto_show_selections))
            ExecutiveSetObjVisib(G,selName,1);
          if(obj->type==cObjectMolecule) {
            if(SettingGet(G,cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1,false);
              sprintf(buffer,"cmd.select('%s',\"%s(%s)\")",selName,sel_mode_kw,buf1);
              PLog(buffer,cPLog_pym);
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
                ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buffer,false);
                sprintf(buf2,"(((%s) or %s(%s)) and ((%s(%s)) and %s(%s)))",
                        selName,sel_mode_kw,buffer,sel_mode_kw,buffer,sel_mode_kw,selName);
                sprintf(buffer,"cmd.select('%s',\"%s(%s)\")",selName,sel_mode_kw,buf2);
                PLog(buffer,cPLog_pym);
              }
            }
          } else {
            sprintf(buf2,"%s(%s)",sel_mode_kw,buffer);
            SelectorCreate(G,selName,buf2,NULL,false,NULL);
            if(obj->type==cObjectMolecule) {
              if(SettingGet(G,cSetting_logging)) {
                objMol = (ObjectMolecule*)obj;            
                ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1,false);
                sprintf(buffer,"cmd.select('%s',\"%s(%s)\")",selName,sel_mode_kw,buf1);
                PLog(buffer,cPLog_pym);
              }
            }
          }
          if(SettingGet(G,cSetting_auto_hide_selections))
            ExecutiveHideSelections(G);
          if(SettingGet(G,cSetting_auto_show_selections))
            ExecutiveSetObjVisib(G,selName,1);
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
	  char name[ObjNameMax];
	  if(ExecutiveGetActiveSeleName(G,name, false)) {
	    SelectorCreate(G,name,"none",NULL,true,NULL);
	    if(SettingGet(G,cSetting_logging)) {
	      sprintf(buf2,"cmd.select('%s','none')\n",name);
	      PLog(buf2,cPLog_no_flush);
	    }
	    SeqDirty(G);
	  }
	}
      case cButModeSeleToggle:
	{
	  OrthoLineType buf2;
	  char name[ObjNameMax];
	  
	  if(ExecutiveGetActiveSeleName(G,name, false)) {
	    ExecutiveSetObjVisib(G,name,0);
	    if(SettingGet(G,cSetting_logging)) {
	      sprintf(buf2,"cmd.disable('%s')\n",name);
	      PLog(buf2,cPLog_no_flush);
	    }
	  }
	}
	break;
      }
      PRINTFB(G,FB_Scene,FB_Warnings) 
	" SceneClick: no atom found nearby.\n"
	ENDFB(G);
      SceneDirty(G); /* this here to prevent display weirdness after
			an unsuccessful picking pass... not sure it helps though*/
      OrthoRestorePrompt(G);
    }
  }
  
  I->StartX = I->LastX;
  I->StartY = I->LastY;

  return(1);
}
void ScenePushRasterMatrix(PyMOLGlobals *G,float *v) 
{
  register CScene *I=G->Scene;
  float scale = SceneGetScreenVertexScale(G,v);
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
float SceneGetScreenVertexScale(PyMOLGlobals *G,float *v1)
{
  /* get conversion factor from screen point to atomic coodinate */
  register CScene *I=G->Scene;
  float vl,p1[4],p2[4];
  float width_factor = I->Width/2.0F;

  if(!v1) v1 = I->Origin;

  /* now, scale properly given the current projection matrix */
  copy3f(v1,p1);
  p1[3] = 1.0;
  MatrixTransform44f4f(I->ModMatrix,p1,p2); /* modelview transformation */
  copy4f(p2,p1);
  p2[0]+=1.0;
  /* projection transformation */
  MatrixTransform44f4f(I->ProMatrix,p1,p1); 
  MatrixTransform44f4f(I->ProMatrix,p2,p2);
  /*  dump44f(I->ProMatrix,"pro");*/
  /* perspective division */
  p1[0]=p1[0]/p1[3];
  p2[0]=p2[0]/p2[3];
  p1[0]=(p1[0]+1.0F)*(width_factor); /* viewport transformation */
  p2[0]=(p2[0]+1.0F)*(width_factor);
  /*
    p1[1]=p1[1]/p1[3];
    p2[1]=p2[1]/p2[3];
    p1[2]=0.0;
    p2[2]=0.0;
    p1[1]=(p1[1]+1.0F)*(I->Height/2.0F);
    p2[1]=(p2[1]+1.0F)*(I->Height/2.0F);
    dump3f(p1,"p1");
    dump3f(p2,"p2");
    vl=(float)diff3f(p1,p2);
  */
  vl = (float)fabs(p1[0]-p2[0]);

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
  PParse(buffer);
  PFlush();
  sprintf(buffer,"cmd.hide('sticks','''%s''')",s);
  PParse(buffer);
  PFlush();
  sprintf(buffer,"cmd.hide('spheres','''%s''')",s);
  PParse(buffer);
  PFlush();
  sprintf(buffer,"cmd.hide('ribbon','''%s''')",s);
  PParse(buffer);
  PFlush();
  sprintf(buffer,"cmd.hide('cartoon','''%s''')",s);
  PParse(buffer);
  PFlush();
  sprintf(buffer,"cmd.hide('labels','''%s''')",s);
  PParse(buffer);
  PFlush();
  sprintf(buffer,"cmd.hide('nonbonded','''%s''')",s);
  PParse(buffer);
  PFlush();
  sprintf(buffer,"cmd.hide('nb_spheres','''%s''')",s);
  PParse(buffer);
  PFlush();
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
      PParse(buffer);
      PFlush();
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
      PParse(buffer);
      PFlush();
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
      PParse(buffer);
      PFlush();
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
      PParse(buffer);
      PFlush();
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
      PParse(buffer);
      PFlush();
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
      PParse(buffer);
      PFlush();

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
"cmd.dist('rov_pc','%s & (elem n+o) & enabled & %s %s (center expand %1.3f)','same',%1.4f,mode=1,labels=%d,quiet=2)",
              s,p1,p2,polar_contacts,polar_cutoff,label_flag);
      PParse(buffer);
      PFlush();

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
      PParse(buffer);
      PFlush();
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
      PParse(buffer);
      PFlush();
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
              PParse(buffer);
              PFlush();
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
              PParse(buffer);
              PFlush();
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
              PParse(buffer);
              PFlush();
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
              PParse(buffer);
              PFlush();
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
              PParse(buffer);
              PFlush();
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
              PParse(buffer);
              PFlush();
              refresh_flag=true;
            }


      SettingSet(G,cSetting_auto_zoom,(float)auto_save);            
    }


    if(refresh_flag) {
      PParse("cmd.refresh()");
      PFlush();
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
  CObject *obj;


  if(I->PossibleSingleClick) {
    double slowest_single_click_drag = 0.15;
    if((when-I->LastClickTime)>slowest_single_click_drag) {
      I->PossibleSingleClick = 0;
    }
  }

  if(I->LoopFlag)
    return SceneLoopDrag(block,x,y,mod);

  mode = ButModeTranslate(G,I->Button,mod);
  
  y=y-I->Block->margin.bottom;


  scale = (float)I->Height;
  if(scale > I->Width)
	 scale = (float)I->Width;
  scale = 0.45F * scale;

  SceneDontCopyNext(G);
  switch(mode) {
  case cButModePickAtom:
    obj=(CObject*)I->LastPicked.ptr;
    if(obj)
      switch(obj->type) {
      case cObjectGadget: {
        ObjectGadget *gad;
        
        gad = (ObjectGadget*)obj;

        ObjectGadgetGetVertex(gad,I->LastPicked.index,I->LastPicked.bond,v1);

        vScale = SceneGetScreenVertexScale(G,v1);
        if(I->StereoMode>1) {
          x = x % (I->Width/2);
          vScale*=2;
        }
        
        /* transform into model coodinate space */
        switch(obj->Context) {
        case 0:
          v2[0] = (x-I->LastX)*vScale;
          v2[1] = (y-I->LastY)*vScale;
          v2[2] = 0;
          MatrixInvTransform44fAs33f3f(I->RotMatrix,v2,v2); 
          break;
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
        }
        add3f(v1,v2,v2);
        ObjectGadgetSetVertex(gad,I->LastPicked.index,I->LastPicked.bond,v2);
        if(gad->Obj.fUpdate)
          gad->Obj.fUpdate((CObject*)gad);
        SceneChanged(G);
        /*        printf("dragging gadget\n");*/
      }
      break;
      }
    I->LastX=x;
    I->LastY=y;
    break;
  case cButModeMovFrag:
  case cButModeTorFrag:
  case cButModeRotFrag:
  case cButModeMoveAtom:
  case cButModePkTorBnd:
    obj=(CObject*)I->LastPicked.ptr;
    if(obj) {
      if(I->Threshold) {
        if((abs(x-I->ThresholdX)>I->Threshold)||
           (abs(y-I->ThresholdY)>I->Threshold)) {
          I->Threshold = 0;
        }
      }
      if(!I->Threshold)
        switch(obj->type) {
        case cObjectMolecule:
          if(ObjectMoleculeGetAtomVertex((ObjectMolecule*)obj,
                                         SettingGetGlobal_i(G,cSetting_state)-1,
                                         I->LastPicked.index,v1)) {
            /* scale properly given the current projection matrix */
            vScale = SceneGetScreenVertexScale(G,v1);
            if(I->StereoMode>1) {
              x = x % (I->Width/2);
              vScale*=2;
            }

            v2[0] = (x-I->LastX)*vScale;
            v2[1] = (y-I->LastY)*vScale;
            v2[2] = 0;
            
            v3[0] = 0.0F;
            v3[1] = 0.0F;
            v3[2] = 1.0F;
            
            /* transform into model coodinate space */
            MatrixInvTransform44fAs33f3f(I->RotMatrix,v2,v2); 
            MatrixInvTransform44fAs33f3f(I->RotMatrix,v3,v3); 

            if(mode!=cButModeMoveAtom) {
              EditorDrag(G,(ObjectMolecule*)obj,I->LastPicked.index,mode,
                         SettingGetGlobal_i(G,cSetting_state)-1,v1,v2,v3);
            } else {
              int log_trans = (int)SettingGet(G,cSetting_log_conformations);
              ObjectMoleculeMoveAtom((ObjectMolecule*)obj,SettingGetGlobal_i(G,cSetting_state)-1,
                                     I->LastPicked.index,v2,1,log_trans);
              SceneDirty(G);
            }
          }
          break;
        case cObjectSlice:
          {
            ObjectSlice *slice = (ObjectSlice*)obj;

            if(I->LastPickVertexFlag) {
              
              copy3f(I->LastPickVertex,v1);

              vScale = SceneGetScreenVertexScale(G,v1);

              if(I->StereoMode>1) {
                x = x % (I->Width/2);
                vScale*=2;
              }
              
              v2[0] = (x-I->LastX)*vScale;
              v2[1] = (y-I->LastY)*vScale;
              v2[2] = 0;
              
              v3[0] = 0.0F;
              v3[1] = 0.0F;
              v3[2] = 1.0F;
              
              /* transform into model coodinate space */
              MatrixInvTransform44fAs33f3f(I->RotMatrix,v2,v2); 
              MatrixInvTransform44fAs33f3f(I->RotMatrix,v3,v3); 

              ObjectSliceDrag(slice,SceneGetState(G),mode,v1,v2,v3);
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

    vScale = SceneGetScreenVertexScale(G,I->Origin);
    if(I->StereoMode>1) {
      x = x % (I->Width/2);
      vScale*=2;
    }

    v2[0] = (x-I->LastX)*vScale;
    v2[1] = (y-I->LastY)*vScale;
    v2[2] = 0.0F;
    
    moved_flag=false;
    if(I->LastX!=x)
      {
        I->Pos[0]+=v2[0];
        I->LastX=x;
        SceneDirty(G);
        moved_flag=true;
      }
    if(I->LastY!=y)
      {
        I->Pos[1]+=v2[1];
        I->LastY=y;
        SceneDirty(G);
        moved_flag=true;
      }
    
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
  case cButModeTransZ:
  case cButModeClipNF:
  case cButModeClipN:    
  case cButModeClipF:    
    
    SceneNoteMouseInteraction(G);

    eff_width = I->Width;
    if(I->StereoMode>1) {
      eff_width = I->Width/2;
      x = x % eff_width;
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
	 case cButModeTransZ:
		if(I->LastY!=y)
		  {
          float factor;
          factor = 200/((I->FrontSafe+I->BackSafe)/2);
          if(factor>=0.0F) {
            factor = (((float)y)-I->LastY)/factor;
            I->Pos[2]+=factor;
            I->Front-=factor;
            I->Back-=factor;
            I->FrontSafe = GetFrontSafe(I->Front,I->Back);
            I->BackSafe= GetBackSafe(I->FrontSafe,I->Back);
          }
          I->LastY=y;
          SceneDirty(G);
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
			 SceneDirty(G);
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
			 SceneDirty(G);
          moved_flag=true;
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
			 SceneDirty(G);
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
			 SceneDirty(G);
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
			 SceneDirty(G);
          moved_flag=true;
		  }
		if(I->LastY!=y)
		  {
			 I->Back-=(((float)y)-I->LastY)/10;
			 if(I->Back<I->Front)
				I->Back=I->Front+cSliceMin;
			 I->LastY=y;
			 SceneDirty(G);
          moved_flag=true;
		  }
		break;
    }
    if(moved_flag)
      SceneDoRoving(G,old_front,old_back,old_origin,adjust_flag,false);
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
  if(!SceneClick(dm->block, dm->button, dm->x, dm->y, dm->mod, dm->when))
    {
      
    }
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


static int SceneDeferClickWhen(Block *block, int button, int x, int y, double when)
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
  OrthoFreeBlock(G,I->Block);
  
  ListFree(I->Obj,next,ObjRec);
  if(!I->MovieOwnsImageFlag)
	 FreeP(I->ImageBuffer);
  
  CGOFree(G->DebugCGO);
  FreeP(G->Scene);
  
}
/*========================================================================*/
void SceneResetMatrix(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  MatrixLoadIdentity44f(I->RotMatrix);
  SceneUpdateInvMatrix(G);
}
/*========================================================================*/
void SceneSetDefaultView(PyMOLGlobals *G) {
  register CScene *I=G->Scene;

  MatrixLoadIdentity44f(I->RotMatrix);
  MatrixLoadIdentity44f(I->ModMatrix);
  MatrixLoadIdentity44f(I->ProMatrix);
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
  SceneDirty(G);
  return(ok);
}
/*========================================================================*/
int  SceneInit(PyMOLGlobals *G)
{
  register CScene *I=NULL;
  if( (I=(G->Scene=Calloc(CScene,1)))) {

    G->DebugCGO = CGONew(G);

    ListInit(I->Obj);
    
    I->LoopFlag = false;
    I->RockTime=0;
    I->TextColor[0]=0.2F;
    I->TextColor[1]=1.0F;
    I->TextColor[2]=0.2F;
    I->SculptingSave=0;
    
    I->LastClickTime = UtilGetSeconds(G);
    I->PossibleSingleClick = 0;
    I->LastReleaseTime = 0.0F;
    I->LastWinX = 0;
    I->LastWinY = 0;
    I->Threshold = 0;
    I->LastPickVertexFlag = false;
    
    SceneSetDefaultView(G);
    
    I->NFrame = 0;
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
    I->RovingDirtyFlag = false;
    I->ImageBuffer = NULL;
    I->ImageBufferWidth=0;
    I->ImageBufferHeight=0;
    I->ImageBufferSize = 0;
    I->MovieOwnsImageFlag = false;
    I->MovieFrameFlag = false;
    I->RenderTime = 0;
    I->LastRender = UtilGetSeconds(G);
    I->LastFrameTime = UtilGetSeconds(G);
    I->LastRockTime = UtilGetSeconds(G);
    I->SingleClickDelay = 0.0;
    I->LastPicked.ptr = NULL;
    
    I->CopyNextFlag=true;
    I->CopyFlag=false;
    I->CopiedFromOpenGL=false;
    I->vendor[0]=0;
    I->renderer[0]=0;
    I->version[0]=0;
    SceneRestartTimers(G);

    I->Width = 400; /* sensible defaults */
    I->Height = 300;

    I->n_ani_elem = 0;
    I->cur_ani_elem = 0;

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

#ifdef _PYMOL_OSX
  /* workaround for broken pixel handling under OSX 
     (Who's fault: Me? Apple? NVidia?) */
  width = 8*(width/8);
#endif

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

  MovieClearImages(G);
  MovieSetSize(G,I->Width,I->Height);

}

float fog_val=1.0;
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

void SceneRay(PyMOLGlobals *G,int ray_width,int ray_height,int mode,
              char **headerVLA_ptr,
              char **charVLA_ptr,float angle,float shift,int quiet,
              G3dPrimitive **g3d)
{
  register CScene *I=G->Scene;
  ObjRec *rec=NULL;
  CRay *ray;
  unsigned int buffer_size;
  float height,width;
  float aspRat;
  float white[3] = {1.0,1.0,1.0};
  unsigned int *buffer;
  float rayView[16];
  int curState;
  double timing;
  char *charVLA = NULL;
  char *headerVLA = NULL;
  float fov;
  OrthoLineType prefix = "";
  SceneUnitContext context;

  if((!ray_width)||(!ray_height)) {
    ray_width=I->Width;
    ray_height=I->Height;
  }

  fov=SettingGet(G,cSetting_field_of_view);
  aspRat = ((float) ray_width) / ((float) ray_height);

  ScenePrepareUnitContext(G,&context,ray_width,ray_height);
  if(SettingGet(G,cSetting_all_states)) {
    curState=-1;
  } else {
    curState=SettingGetGlobal_i(G,cSetting_state)-1;
  }

  ray = RayNew(G);

  SceneUpdate(G);

  timing = UtilGetSeconds(G); /* start timing the process */
  
  /* start afresh, looking in the negative Z direction (0,0,-1) from (0,0,0) */
  MatrixLoadIdentity44f(rayView);

  /* move the camera to the location we are looking at */
  MatrixTranslate44f3f(rayView,I->Pos[0],I->Pos[1],I->Pos[2]);

  if(shift) {
    MatrixTranslate44f3f(rayView,shift,0.0F,0.0F);
  }
  /* move the camera so that we can see the origin 
	* NOTE, vector is given in the coordinates of the world's motion
	* relative to the camera */

  
  /* 4. rotate about the origin (the the center of rotation) */

  if(angle) {
    float temp[16];
    MatrixLoadIdentity44f(temp);
    MatrixRotate44f3f(temp,(float)(-PI*angle/180),0.0F,1.0F,0.0F);
    MatrixMultiply44f(I->RotMatrix,temp);
    MatrixMultiply44f(temp,rayView);
  } else {
    MatrixMultiply44f(I->RotMatrix,rayView);
  }


  /* 5. move the origin to the center of rotation */
  MatrixTranslate44f3f(rayView,-I->Origin[0],-I->Origin[1],-I->Origin[2]);


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
    int ortho = SettingGetGlobal_i(G,cSetting_ray_orthoscopic);
    int ray_pixel_width;

    if(ortho<0) ortho = SettingGetGlobal_b(G,cSetting_ortho);
    
    if(SettingGetGlobal_b(G,cSetting_ray_pixel_scale_to_window)) {
      ray_pixel_width = I->Width;
    }  else {
      ray_pixel_width = ray_width;
    }
    if(ortho) {
      RayPrepare(ray,-width,width,-height,height,
                 I->FrontSafe,I->BackSafe,rayView,I->RotMatrix,aspRat,ray_pixel_width);
    } else {
      float back_ratio;
      float back_height;
      float back_width;

      back_ratio = -I->Back/I->Pos[2];
      back_height = back_ratio*height;
      back_width = aspRat * back_height;
      RayPrepare(ray,-back_width, back_width, -back_height, back_height,
                 I->FrontSafe,I->BackSafe,rayView,I->RotMatrix,aspRat,ray_pixel_width);
    }
  }

  while(ListIterate(I->Obj,rec,next))
	 {
		if(rec->obj->fRender) {
        RaySetContext(ray,rec->obj->Context);
		  ray->fColor3fv(ray,white);
		  rec->obj->fRender(rec->obj,
                          ObjectGetCurrentState(rec->obj,false),ray,NULL,0);
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
    buffer_size = 4*ray_width*ray_height;
    buffer=(GLvoid*)Alloc(char,buffer_size);
    ErrChkPtr(G,buffer);
    
    RayRender(ray,ray_width,ray_height,buffer,I->FrontSafe,I->BackSafe,timing,angle,
              fov,I->Pos);
    SceneApplyImageGamma(G,buffer,ray_width,ray_height);

    /*    RayRenderColorTable(ray,ray_width,ray_height,buffer);*/
    
    if(I->ImageBuffer) {
      if(I->MovieOwnsImageFlag) {
        I->MovieOwnsImageFlag=false;
        I->ImageBuffer=NULL;
      } else {
        FreeP(I->ImageBuffer);
      }
    }
    
    I->ImageBuffer = buffer;
    I->ImageBufferSize = buffer_size;
    I->ImageBufferWidth=ray_width;
    I->ImageBufferHeight=ray_height;
    I->DirtyFlag=false;
    I->CopyFlag = true;
    I->CopiedFromOpenGL = false;
    I->MovieOwnsImageFlag = false;
    break;

  case 1: /* mode 1 is povray */
    charVLA=VLAlloc(char,100000); 
    headerVLA=VLAlloc(char,2000);
    RayRenderPOV(ray,ray_width,ray_height,&headerVLA,&charVLA,
                 I->FrontSafe,I->BackSafe,fov,angle);
    if(!(charVLA_ptr&&headerVLA_ptr)) { /* immediate mode */
      strcpy(prefix,SettingGet_s(G,NULL,NULL,cSetting_batch_prefix));
      if(PPovrayRender(headerVLA,charVLA,prefix,ray_width,
                       ray_height,(int)SettingGet(G,cSetting_antialias))) {
        strcat(prefix,".png");
        SceneLoadPNG(G,prefix,false,false);
        I->DirtyFlag=false;
      }
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
  }
  timing = UtilGetSeconds(G)-timing;
  if(mode!=2) { /* don't show timings for tests */
	accumTiming += timing; 

	if(!quiet) {
     PRINTFB(G,FB_Ray,FB_Details)
       " Ray: total time: %4.2f sec. = %3.1f frames/hour. (%4.2f sec. accum.)\n", 
       timing,3600/timing, 
       accumTiming 
      ENDFB(G);
   }
  }
  if(mode!=3)
    OrthoDirty(G);
  RayFree(ray);
}
/*========================================================================*/
void SceneCopy(PyMOLGlobals *G,GLenum buffer,int force)
{
  register CScene *I=G->Scene;
  unsigned int buffer_size;

  if(force || (!(I->StereoMode||SettingGet(G,cSetting_stereo_double_pump_mono))))
  { /* no copies while in stereo mode */
    if((!I->DirtyFlag)&&(!I->CopyFlag)) { 
      buffer_size = 4*I->Width*I->Height;
      if(buffer_size) {
        if(I->ImageBuffer)	 {
          if(I->MovieOwnsImageFlag) {
            I->MovieOwnsImageFlag=false;
            I->ImageBuffer=NULL;
          } else if(I->ImageBufferSize!=buffer_size) {
            FreeP(I->ImageBuffer);
          }
        }
        if((I->ImageBufferWidth!=I->Width)||(I->ImageBufferHeight!=I->Height)) {
          FreeP(I->ImageBuffer);
        }
        if(!I->ImageBuffer) {
          I->ImageBuffer=(GLvoid*)Alloc(char,buffer_size);
          ErrChkPtr(G,I->ImageBuffer);
          I->ImageBufferSize = buffer_size;
          I->ImageBufferWidth=I->Width;
          I->ImageBufferHeight=I->Height;
        }
        if(G->HaveGUI && G->ValidContext) {
          glReadBuffer(buffer);
          PyMOLReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                       GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);
          I->ImageBufferWidth=I->Width;
          I->ImageBufferHeight=I->Height;
        }
      }
      I->CopyFlag = true;
      I->CopiedFromOpenGL = true;
    }
  }
}

/*========================================================================*/
int SceneRovingCheckDirty(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;

  return(I->RovingDirtyFlag);
}
/*========================================================================*/
void SceneUpdate(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  ObjRec *rec=NULL;

  PRINTFD(G,FB_Scene)
    " SceneUpdate: entered.\n"
    ENDFD;
  if(I->ChangedFlag) {
    SceneCountFrames(G);
	 while(ListIterate(I->Obj,rec,next))
      if(rec->obj->fUpdate) 
        rec->obj->fUpdate(rec->obj);
	 I->ChangedFlag=false;
    if(!MovieDefined(G)) {
      if(SettingGetGlobal_i(G,cSetting_frame)!=
         SettingGetGlobal_i(G,cSetting_state))
        SettingSetGlobal_i(G,cSetting_frame,SettingGetGlobal_i(G,cSetting_state));
    }
    WizardDoScene(G);
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
  ImageType image;
  int renderedFlag=false;

  PRINTFD(G,FB_Scene)
    " SceneRenderCached: entered.\n"
    ENDFD;

  if(I->DirtyFlag) {
	if(I->MovieFrameFlag||
	   (MoviePlaying(G)&&SettingGet(G,cSetting_cache_frames))) {
	  I->MovieFrameFlag=false;
	  image = MovieGetImage(G,
                           MovieFrameToImage(G,
                                             SettingGetGlobal_i(G,cSetting_frame)-1));
	  if(image)
		{
		  if(I->ImageBuffer)
			{
			  if(!I->MovieOwnsImageFlag) {
				mfree(I->ImageBuffer);
			  }
			}
		  I->MovieOwnsImageFlag=true;
		  I->CopyFlag=true;
		  I->ImageBuffer=image;
		  OrthoDirty(G);
		  renderedFlag=true;
		}
	  else
		{
        SceneMakeMovieImage(G);
        renderedFlag=true;
		}
	} else if(MoviePlaying(G)&&SettingGet(G,cSetting_ray_trace_frames)) {
	  SceneRay(G,0,0,(int)SettingGet(G,cSetting_ray_default_renderer),NULL,NULL,0.0F,0.0F,false,NULL); 
	} else {
	  renderedFlag=false;
	  I->CopyFlag = false;
	}
	I->DirtyFlag=false;
  } else if(I->CopyFlag) {
	renderedFlag=true;
  }
  /*  if(renderedFlag) {
	I->RenderTime = -I->LastRender;
	I->LastRender = UtilGetSeconds(G);
	I->RenderTime += I->LastRender;
	ButModeSetRate(G,I->RenderTime);
   }*/

  PRINTFD(G,FB_Scene)
    " SceneRenderCached: leaving...renderedFlag %d\n",renderedFlag
    ENDFD;

  return(renderedFlag);
}


/*========================================================================*/
static void SceneRenderAll(PyMOLGlobals *G,SceneUnitContext *context,float *normal,Pickable **pickVLA,int pass,int fat)
{
  register CScene *I=G->Scene;
  ObjRec *rec=NULL;
  float vv[4];

  while(ListIterate(I->Obj,rec,next))
    {
      glPushMatrix();
      if(fat)
        glLineWidth(3.0);
      if(rec->obj->fRender)
        switch(rec->obj->Context) {
        case 0:
          if(normal) 
            glNormal3fv(normal);
          rec->obj->fRender(rec->obj,
                            ObjectGetCurrentState(rec->obj,false),NULL,pickVLA,pass);
          break;
        case 1:
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

          glOrtho(context->unit_left,
                  context->unit_right,
                  context->unit_top,
                  context->unit_bottom,
                  context->unit_front,
                  context->unit_back);
          
          glNormal3f(0.0F,0.0F,1.0F);
          rec->obj->fRender(rec->obj,
                            ObjectGetCurrentState(rec->obj,false),NULL,pickVLA,pass);

          glMatrixMode(GL_MODELVIEW);
          glLoadIdentity();
          vv[0]=0.0;
          vv[1]=0.0;
          vv[2]=1.0;
          vv[3]=0.0;
          glLightfv(GL_LIGHT0,GL_POSITION,vv);
          glLightfv(GL_LIGHT1,GL_POSITION,vv);

          glPopMatrix();
          glMatrixMode(GL_PROJECTION);
          glPopMatrix();
          glMatrixMode(GL_MODELVIEW);
          break;
        }
      glPopMatrix();
    }
}

#ifdef _PYMOL_SHARP3D
void sharp3d_begin_left_stereo(void);
void sharp3d_switch_to_right_stereo(void);
void sharp3d_end_stereo(void);
#endif
/*========================================================================*/
void SceneRender(PyMOLGlobals *G,Pickable *pick,int x,int y,Multipick *smp)
{
  /* think in terms of the camera's world */
  register CScene *I=G->Scene;
  float fog[4];
  float *v,vv[4],f;
  unsigned int lowBits,highBits;
  unsigned int *lowBitVLA=NULL,*highBitVLA=NULL;
  int high,low;
  static float white[4] =
  {1.0, 1.0, 1.0, 1.0};
  float zero[4] = {0.0,0.0,0.0,0.0};
  float zAxis[4] = { 0.0, 0.0, 1.0, 0.0 };
  float normal[4] = { 0.0, 0.0, 1.0, 0.0 };
  float aspRat = ((float) I->Width) / ((float) I->Height);
  float height,width;
  double start_time=0.0;
  int view_save[4];
  Pickable *pickVLA,*pik;
  int lastIndex=0;
  void *lastPtr=NULL;
  int index;
  int curState;
  int nPick,nHighBits,nLowBits;
  int pass;
  float fov;
  float fog_start;
  int double_pump = false;
  int must_render_stereo = false;
  int stereo_as_mono = false;
  int debug_pick = 0;
  double now;
  GLenum render_buffer = GL_BACK;
  SceneUnitContext context;
  
  PRINTFD(G,FB_Scene)
    " SceneRender: entered. pick %p x %d y %d smp %p\n",
    (void*)pick,x,y,(void*)smp
    ENDFD;

  if(I->cur_ani_elem < I->n_ani_elem ) { /* play motion animation */
    int cur = I->cur_ani_elem;

    now = UtilGetSeconds(G);

    while(I->ani_elem[cur].timing<now) {
      cur++;
      if(cur >= I->n_ani_elem) {
        cur = I->n_ani_elem;
        break;
      }
    }
    I->cur_ani_elem = cur;
    SceneFromViewElem(G,I->ani_elem+cur);


  }

  double_pump=SettingGet_i(G,NULL,NULL,cSetting_stereo_double_pump_mono);
  
  if(I->StereoMode>1)
    aspRat=aspRat/2;

  fov=SettingGet(G,cSetting_field_of_view);
  if(G->HaveGUI && G->ValidContext) {
    
    if(Feedback(G,FB_OpenGL,FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender checkpoint 0");

    must_render_stereo = (I->StereoMode!=0);
    if(!must_render_stereo) 
      if(double_pump&&G->StereoCapable) {            /* force stereo rendering */
        must_render_stereo=true;
        stereo_as_mono=true; /* rendering stereo as mono */
      }

    if(must_render_stereo) {
      glDrawBuffer(GL_BACK_LEFT);
      render_buffer = GL_BACK_LEFT;
    } else {
      glDrawBuffer(GL_BACK);
      render_buffer = GL_BACK;
    }

    if(Feedback(G,FB_OpenGL,FB_Debugging))
      PyMOLCheckOpenGLErr("SceneRender checkpoint 1");

  
    glGetIntegerv(GL_VIEWPORT,(GLint*)view_save);
    glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height);

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

    /* load up the light positions relative to the camera while MODELVIEW still has the identity */

    vv[0] = 0.0F;
    vv[1] = 0.0F;
    vv[2] = 1.0F;
    vv[3]=0.0;
    glLightfv(GL_LIGHT0,GL_POSITION,vv);

    copy3f(SettingGetGlobal_3fv(G,cSetting_light),vv);
    normalize3f(vv);
    invert3f(vv);
    vv[3]=0.0;
    glLightfv(GL_LIGHT1,GL_POSITION,vv);

    ScenePrepareUnitContext(G,&context,I->Width,I->Height);
 
 
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
  
    /* determine the direction in which we are looking relative*/

    /* 2. set the normals to reflect light back at the camera */

    MatrixInvTransform3f(I->RotMatrix,zAxis,normal); 
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
    }  

    if(!(pick||smp)) {


      glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, white);
      glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

      vv[0]=SettingGet(G,cSetting_shininess);
      glMaterialfv(GL_FRONT,GL_SHININESS,vv);

      glShadeModel(GL_SMOOTH);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_DITHER);

      glAlphaFunc(GL_GREATER, 0.05F);
      glEnable(GL_ALPHA_TEST);

      if(G->Option->multisample)
        glEnable(0x809D); /* GL_MULTISAMPLE_ARB */

      {
        float direct = SettingGetGlobal_f(G,cSetting_direct);
        float reflect = SettingGetGlobal_f(G,cSetting_reflect) * (1.0F - direct);

        if(reflect>1.0F) reflect=1.0F;

        /* lighting */
        
        glEnable(GL_LIGHTING);
        
        /* lighting model: one or two sided? */
        
        if(SettingGet(G,cSetting_two_sided_lighting)||
           (SettingGetGlobal_i(G,cSetting_transparency_mode)==1)) {
          glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
        } else {
          glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
        }     
        
        /* add half the ambient component (perceptive kludge) */
        
        f=SettingGet(G,cSetting_ambient) * 0.5F;

        vv[0]=f;
        vv[1]=f;
        vv[2]=f;
        vv[3]=1.0;
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT,vv);
        
        /* LIGHT0 is or direct light (eminating from the camera) */
        
        if(direct>R_SMALL4) {          
          glEnable(GL_LIGHT0);
          
          vv[0]=direct;
          vv[1]=direct;
          vv[2]=direct;
          vv[3]=1.0F;
          glLightfv(GL_LIGHT0,GL_DIFFUSE,vv);
          
          f=SettingGet(G,cSetting_ambient);
          vv[0] = f;
          vv[1] = f;
          vv[2] = f;
          vv[3] = 1.0F;
          glLightfv(GL_LIGHT0,GL_AMBIENT,vv);
          
          glLightfv(GL_LIGHT0,GL_SPECULAR,zero);
          
        } else {
          glDisable(GL_LIGHT0);
        }

        
        /* LIGHT1 is our reflected light (specular and diffuse
           reflection from a movable directional light) */
        
        glEnable(GL_LIGHT1);
        
        /* load up the default ambient and diffuse lighting values */
        
        vv[0]=reflect;
        vv[1]=reflect;
        vv[2]=reflect;
        vv[3]=0.0F;
        glLightfv(GL_LIGHT1,GL_DIFFUSE,vv);

        vv[0] = 0.0F;
        vv[1] = 0.0F;
        vv[2] = 0.0F;
        vv[3] = 1.0F;
        glLightfv(GL_LIGHT1,GL_AMBIENT,vv);
        
        f = SettingGet(G,cSetting_specular);
        if(f==1.0F) {
          f=SettingGet(G,cSetting_specular_intensity);
        } 
        if(f>R_SMALL4) { /* setup specular reflections for LIGHT0 if needed */
          
          vv[0]=f;
          vv[1]=f;
          vv[2]=f;
          vv[3]=1.0;
          
          glLightfv(GL_LIGHT1,GL_SPECULAR,vv);
          glMaterialfv(GL_FRONT,GL_SPECULAR,vv);
        } else { /* no specular reflections */
          
          vv[0]=0.0;
          vv[1]=0.0;
          vv[2]=0.0;
          vv[3]=1.0;
          
          glLightfv(GL_LIGHT1, GL_SPECULAR, vv);
          glMaterialfv(GL_FRONT,GL_SPECULAR,vv); 
        }
        
        /* 

        glGetFloatv(GL_LIGHT0,vv);
        printf("glGetFloatv(GL_LIGHT0) %8.3f\n",vv[0]);
        
        glGetLightfv(GL_LIGHT0,GL_AMBIENT,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_AMBIENT,vv)");
        
        glGetLightfv(GL_LIGHT0,GL_DIFFUSE,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_DIFFUSE,vv)");
        
        glGetLightfv(GL_LIGHT0,GL_SPECULAR,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_SPECULAR,vv)");
        
        glGetLightfv(GL_LIGHT0,GL_POSITION,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT0,GL_POSITION,vv)");
        
        glGetFloatv(GL_LIGHT1,vv);
        printf("glGetFloatv(GL_LIGHT1) %8.3f\n",vv[0]);
        
        glGetLightfv(GL_LIGHT1,GL_AMBIENT,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_AMBIENT,vv)");
        
        glGetLightfv(GL_LIGHT1,GL_DIFFUSE,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_DIFFUSE,vv)");
        
        glGetLightfv(GL_LIGHT1,GL_SPECULAR,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_SPECULAR,vv)");

        glGetLightfv(GL_LIGHT1,GL_POSITION,vv);
        dump4f(vv, "glGetLightfv(GL_LIGHT1,GL_POSITION,vv)");

        glGetFloatv(GL_LIGHT_MODEL_AMBIENT,vv);
        dump4f(vv, "glGetFloatv(GL_LIGHT_MODEL_AMBIENT,vv)");

        glGetMaterialfv(GL_FRONT, GL_SPECULAR, vv);
        dump4f(vv, "glGetMaterialfv(GL_FRONT,GL_SPECULAR,vv)");
        */

      }


      if(SettingGet(G,cSetting_depth_cue)&&SettingGet(G,cSetting_fog)) {
        fog_start = (I->Back-I->FrontSafe)*SettingGet(G,cSetting_fog_start)+I->FrontSafe;
#ifdef _PYMOL_3DFX
        if(SettingGet(G,cSetting_ortho)==0.0) {
#endif
          glEnable(GL_FOG);
          glFogf(GL_FOG_MODE, GL_LINEAR);
          glHint(GL_FOG_HINT,GL_NICEST);
          glFogf(GL_FOG_START, fog_start);
#ifdef _PYMOL_3DFX
          if(I->Back>(I->FrontSafe*4.0))
            glFogf(GL_FOG_END, I->Back);
          else
            glFogf(GL_FOG_END,I->FrontSafe*4.0);
          fog_val+=0.0000001;
          if(fog_val>1.0) fog_val=0.99999;
          glFogf(GL_FOG_DENSITY, fog_val);
#else
          glFogf(GL_FOG_END, I->FrontSafe+(I->Back-I->FrontSafe)/SettingGet(G,cSetting_fog));
          glFogf(GL_FOG_DENSITY, fog_val);
#endif
          v=SettingGetfv(G,cSetting_bg_rgb);
          fog[0]=v[0];
          fog[1]=v[1];
          fog[2]=v[2];
          fog[3]=1.0;
          glFogfv(GL_FOG_COLOR, fog);
#ifdef _PYMOL_3DFX
        } else {
          glDisable(GL_FOG);
        }
#endif
      } else {
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
    if(pick) {
      /* atom picking HACK - obfuscative coding */

      switch(I->StereoMode) {
      case 2:
      case 3:
        glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
        break;
      }

      glClearColor(0.0,0.0,0.0,0.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


      pickVLA=VLAlloc(Pickable,5000);
      pickVLA[0].index=0;
      pickVLA[0].ptr=NULL;

      SceneRenderAll(G,&context,NULL,&pickVLA,0,true);
	  

      if(debug_pick) {
        PyMOL_SwapBuffers(G->PyMOL);
        PSleep(1000000*debug_pick/4);
        PyMOL_SwapBuffers(G->PyMOL);
      }
      lowBits = SceneFindTriplet(G,x,y,render_buffer);
	
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA[0].index=0;
      pickVLA[0].ptr=(void*)pick; /* this is just a flag */
	
      SceneRenderAll(G,&context,NULL,&pickVLA,0,true);

      if(debug_pick) {
        PyMOL_SwapBuffers(G->PyMOL);
        PSleep(1000000*debug_pick/4);
        PyMOL_SwapBuffers(G->PyMOL);
      }

      highBits = SceneFindTriplet(G,x,y,render_buffer);
      index = lowBits+(highBits<<12);

      if(debug_pick) {
        PRINTFB(G,FB_Scene,FB_Details)
          " SceneClick-Detail: index %d < %d?\n",index,pickVLA[0].index
          ENDFB(G);
      }
      
      if(index&&(index<=pickVLA[0].index)) {
        *pick = pickVLA[index]; /* return object info */
        if(debug_pick) {
          PRINTFB(G,FB_Scene,FB_Details)
            " SceneClick-Detail: obj %p index %d bond %d\n",pick->ptr,pick->index,pick->bond
            ENDFB(G);
        }
      } else {
        pick->ptr = NULL;
      }
      
		VLAFree(pickVLA);
		
    } else if(smp) {

      switch(I->StereoMode) {
      case 2:
      case 3:
        glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
        break;
      }

      /* multiple atom picking HACK - even more obfuscative coding */

      glClearColor(0.0,0.0,0.0,0.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA=VLAlloc(Pickable,5000);
      pickVLA[0].index=0;
      pickVLA[0].ptr=NULL;
      
      SceneRenderAll(G,&context,NULL,&pickVLA,0,true);
      
      lowBitVLA = SceneReadTriplets(G,smp->x,smp->y,smp->w,smp->h,render_buffer);

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA[0].index=0;
      pickVLA[0].ptr=(void*)smp; /* this is just a flag */
	
      SceneRenderAll(G,&context,NULL,&pickVLA,0,true);

      highBitVLA = SceneReadTriplets(G,smp->x,smp->y,smp->w,smp->h,render_buffer);
      
      nLowBits = VLAGetSize(lowBitVLA);
      nHighBits = VLAGetSize(highBitVLA);
      nPick=0;
      if(nLowBits&&nHighBits) {
		  low = 0;
		  high = 0;
		  while((low<nLowBits)&&(high<nHighBits)) {
          
          if(lowBitVLA[low+1]==highBitVLA[high+1]) {
            index = lowBitVLA[low]+(highBitVLA[high]<<12);
            if(index&&(index<=pickVLA[0].index)) {          
              pik = pickVLA+index; /* just using as a tmp */
              if((pik->index!=lastIndex)||(pik->ptr!=lastPtr))
                {
                  if(((CObject*)pik->ptr)->type==cObjectMolecule) {
                    nPick++; /* start from 1 */
                    VLACheck(smp->picked,Pickable,nPick);
                    smp->picked[nPick] = *pik; /* return atom/object info -- will be redundant */
                  }
                  lastIndex=pik->index;                
                  lastPtr=pik->ptr;
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

      smp->picked[0].index=nPick;

		VLAFree(pickVLA);
      VLAFreeP(lowBitVLA);
      VLAFreeP(highBitVLA);
    } else {

      
      /* STANDARD RENDERING */

      ButModeCaptionReset(G); /* reset the frame caption if any */
      /* rendering for visualization */


      PRINTFD(G,FB_Scene)
        " SceneRender: I->StereoMode %d must_render_stereo %d\n    stereo_as_mono %d  StereoCapable %d\n",
        I->StereoMode, must_render_stereo, stereo_as_mono, G->StereoCapable
        ENDFD;

      start_time = UtilGetSeconds(G);
      if(must_render_stereo) {
        /*stereo*/

        PRINTFD(G,FB_Scene)
          " SceneRender: left hand stereo...\n"
          ENDFD;

        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("before stereo glViewport 1");
        
        /* LEFT HAND STEREO */

        if(stereo_as_mono) {
          glDrawBuffer(GL_BACK_LEFT);
        } else switch(I->StereoMode) {
        case 1: /* hardware */
#ifdef _PYMOL_SHARP3D
          sharp3d_begin_left_stereo();
#else
          glDrawBuffer(GL_BACK_LEFT);
#endif
          break;
        case 2: /* side by side */
          glViewport(I->Block->rect.left+I->Width/2,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        case 3:
          glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        }

        /* prepare the stereo transformation matrix */

        glPushMatrix(); /* 1 */
        ScenePrepareMatrix(G,stereo_as_mono ? 0 : 1);

        /* render picked atoms */

        glPushMatrix(); /* 2 */
        EditorRender(G,curState);
        glPopMatrix(); /* 1 */

        /* render the debugging CGO */

        glPushMatrix();  /* 2 */
        glNormal3fv(normal);
        CGORenderGL(G->DebugCGO,NULL,NULL,NULL);
        glPopMatrix();  /* 1 */

        /* render all objects */

        glPushMatrix(); /* 2 */

        for(pass=1;pass>-2;pass--) { /* render opaque, then antialiased, then transparent...*/
          SceneRenderAll(G,&context,normal,NULL,pass,false);
        }
        glPopMatrix(); /* 1 */

        /* render selections */
        glPushMatrix(); /* 2 */
        glNormal3fv(normal);
        ExecutiveRenderSelections(G,curState);
        glPopMatrix(); /* 1 */
        
        glPopMatrix(); /* 0 */

        /* RIGHT HAND STEREO */

        PRINTFD(G,FB_Scene)
          " SceneRender: right hand stereo...\n"
          ENDFD;

        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("before stereo glViewport 2");

        if(stereo_as_mono) { /* double pumped mono */
          glDrawBuffer(GL_BACK_RIGHT);
        } else switch(I->StereoMode) {
        case 1: /* hardware */
#ifdef _PYMOL_SHARP3D
          sharp3d_switch_to_right_stereo();
#else
          glDrawBuffer(GL_BACK_RIGHT);
#endif
          break;
        case 2: /* side by side */
          glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        case 3:
          glViewport(I->Block->rect.left+I->Width/2,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        }

        /* prepare the stereo transformation matrix */

        glPushMatrix(); /* 1 */
        ScenePrepareMatrix(G,stereo_as_mono ? 0 : 2);
        glClear(GL_DEPTH_BUFFER_BIT) ;

        /* render picked atoms */

        glPushMatrix(); /* 2 */
        EditorRender(G,curState);
        glPopMatrix(); /* 1 */

        /* render the debugging CGO */

        glPushMatrix();  /* 2 */
        glNormal3fv(normal);
        CGORenderGL(G->DebugCGO,NULL,NULL,NULL);
        glPopMatrix();  /* 1 */

        /* render all objects */

        glPushMatrix(); /* 2 */
        for(pass=1;pass>-2;pass--) { /* render opaque, then antialiased, then transparent...*/
          SceneRenderAll(G,&context,normal,NULL,pass,false);
        }        
        glPopMatrix(); /* 1 */

        /* render selections */
        glPushMatrix(); /* 2 */
        glNormal3fv(normal);
        ExecutiveRenderSelections(G,curState);
        glPopMatrix(); /* 1 */

        glPopMatrix(); /* 0 */

        /* restore draw buffer */

        if(stereo_as_mono) { /* double pumped mono */
          glDrawBuffer(GL_BACK);
        } else switch(I->StereoMode) {
        case 1: /* hardware */
#ifdef _PYMOL_SHARP3D
          sharp3d_end_stereo();
#else
          glDrawBuffer(GL_BACK);
#endif
          break;
        case 2: /* side by side */
        case 3:
          glDrawBuffer(GL_BACK);
          break;
        }

      } else {

        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("Before mono rendering");

        /* mono rendering */

        PRINTFD(G,FB_Scene)
          " SceneRender: rendering DebugCGO...\n"
          ENDFD;
        
        glPushMatrix();
        glNormal3fv(normal);
        CGORenderGL(G->DebugCGO,NULL,NULL,NULL);
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
          SceneRenderAll(G,&context,normal,NULL,pass,false);
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
        SceneRenderAll(G,&context,normal,NULL,-1,false);
        glPopMatrix();

        if(Feedback(G,FB_OpenGL,FB_Debugging))
          PyMOLCheckOpenGLErr("during mono rendering");
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
    ButModeSetRate(G,(float)I->RenderTime);
    if(I->CopyNextFlag) {
      start_time = I->LastRender - start_time;
      if((start_time>0.10)||(MainSavingUnderWhileIdle()))
        if(!(ControlIdling(G)))
          if(SettingGet(G,cSetting_cache_display)) {
            SceneCopy(G,render_buffer,false);
          }
    } else {
      I->CopyNextFlag=true;
    }
  }
  PRINTFD(G,FB_Scene)
    " SceneRender: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void SceneRestartTimers(PyMOLGlobals *G)
{
  register CScene *I=G->Scene;
  I->LastRender = UtilGetSeconds(G);
  I->LastFrameTime = UtilGetSeconds(G);
  I->RenderTime = 0;
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

    /* rotate about the origin (the the center of rotation) */
    glMultMatrixf(I->RotMatrix);			

    /* move the origin to the center of rotation */
    glTranslatef(-I->Origin[0],-I->Origin[1],-I->Origin[2]);
  }

}
/*========================================================================*/
void SceneRotate(PyMOLGlobals *G,float angle,float x,float y,float z)
{
  register CScene *I=G->Scene;
  float temp[16];
  int a;
  angle = (float)(-PI*angle/180.0);
  MatrixLoadIdentity44f(temp);
  MatrixRotate44f3f(temp,angle,x,y,z);
  MatrixMultiply44f(I->RotMatrix,temp);
  for(a=0;a<16;a++) 
    I->RotMatrix[a]=temp[a];
  SceneUpdateInvMatrix(G);

  SceneDirty(G);

    /*  glPushMatrix();
        glLoadIdentity();
        glRotatef(angle,x,y,z);
        glMultMatrixf(I->RotMatrix);
        glGetFloatv(GL_MODELVIEW_MATRIX,I->RotMatrix);
        glPopMatrix();*/
}
/*========================================================================*/
void SceneApplyMatrix(PyMOLGlobals *G,float *m)
{
  register CScene *I=G->Scene;
  MatrixMultiply44f(m,I->RotMatrix);
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
  SceneDirty(G);
}










