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
#include"Grap.h"
#include"ObjectGadget.h"

#define cFrontMin 0.1F
#define cSliceMin 0.1F

#define SceneLineHeight 12
#define SceneTopMargin 0
#define SceneBottomMargin 3
#define SceneLeftMargin 3

typedef struct ObjRec {
  struct CObject *obj;  
  struct ObjRec *next;
} ObjRec;

float SceneGetScreenVertexScale(float *v1);

ListVarDeclare(ObjList,ObjRec);

typedef struct {
  Block *Block;
  ObjRec *Obj;
  float RotMatrix[16];
  float InvMatrix[16];
  float ModMatrix[16];
  float VewMatrix[16];
  float ProMatrix[16];
  float UnitMatrix[16];
  float Scale;
  int Width,Height;
  int Button;
  int LastX,LastY;
  float ViewNormal[3],LinesNormal[3];
  float Pos[3],Origin[3];
  float H;
  float Front,Back,FrontSafe;
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
  Pickable LastPicked;
  int StereoMode;
  OrthoLineType vendor,renderer,version;
  int SculptingFlag,SculptingSave;
  int RovingDirtyFlag;
  int RovingCleanupFlag;
  double RovingLastUpdate;
} CScene;

CScene Scene;

typedef struct {
  float unit_left,unit_right,unit_top,unit_bottom,unit_front,unit_back;
} SceneUnitContext;

void SceneCopy(GLenum buffer,int force);

unsigned int SceneFindTriplet(int x,int y,GLenum gl_buffer);
unsigned int *SceneReadTriplets(int x,int y,int w,int h,GLenum gl_buffer);

void SceneDraw(Block *block);
int SceneClick(Block *block,int button,int x,int y,int mod);
int SceneRelease(Block *block,int button,int x,int y,int mod);

int SceneDrag(Block *block,int x,int y,int mod);
void ScenePrepareMatrix(int mode);

void ScenePrepareUnitContext(SceneUnitContext *context,int width,int height);

#if 0
static int SceneGetObjState(CObject *obj,int state)
{
  int objState;
  if(SettingGetIfDefined_i(obj->Setting,cSetting_state,&objState)) {
    if(objState>0) { /* specific state */
      state=objState-1;
    } if(objState<0) { /* all states */
      state=-1;
    }
  }
  if(state>=0) { /* if all states for object is set */
    if(SettingGet_i(obj->Setting,NULL,cSetting_all_states))
      state=-1;
  }
  return(state);
}
#endif

void SceneCleanupStereo(void)
{
  CScene *I=&Scene;  
  if(I->StereoMode==1)
    PSGIStereo(0);
}

void ScenePrepareUnitContext(SceneUnitContext *context,int width,int height)
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

  PRINTFD(FB_Scene)
    "ScenePrepareUnitContext:%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
    context->unit_left,
    context->unit_right,
    context->unit_top, 
    context->unit_bottom,
    context->unit_front,
    context->unit_back
    ENDFD;

}

void SceneSetCardInfo(char *vendor,char *renderer,char *version){
  CScene *I=&Scene;  
  UtilNCopy(I->vendor,vendor,sizeof(OrthoLineType)-1);
  UtilNCopy(I->renderer,renderer,sizeof(OrthoLineType)-1);
  UtilNCopy(I->version,version,sizeof(OrthoLineType)-1);
}

int SceneGetStereo(void)
{
  CScene *I=&Scene;  
  return(I->StereoMode);
}
void SceneGetCardInfo(char **vendor,char **renderer,char **version)
{
  CScene *I=&Scene;  
  (*vendor)=I->vendor;
  (*renderer)=I->renderer;
  (*version)=I->version;
}

void SceneGetPos(float *pos)
{
  CScene *I=&Scene;  
  
  PRINTFD(FB_Scene)
    " SceneGetPos: origin of rotation"
    ENDFD3f(I->Origin);
  /* take origin into camera coords */

  MatrixTransform3f(I->RotMatrix,I->Origin,pos); 

  PRINTFD(FB_Scene)
    " SceneGetPos: origin in camera  "
    ENDFD3f(pos);

  /* find offset in camera coordinates */

  pos[0]=pos[0]-I->Pos[0]; 
  pos[1]=pos[1]-I->Pos[1];
  /*  pos[2]=pos[2]+I->Pos[2]; +(I->FrontSafe+I->Back)/2; */
  PRINTFD(FB_Scene)
    " SceneGetPos: center in camera  "
    ENDFD3f(pos);

  /* convert back to real coordinates */

  MatrixInvTransform3f(I->RotMatrix,pos,pos);

  PRINTFD(FB_Scene)
    " SceneGetPos: center            "
    ENDFD3f(pos);

}
/*========================================================================*/
void SceneApplyRotMatrix(float *src,float *dst)
{
  CScene *I=&Scene;
  MatrixTransform3f(I->RotMatrix,src,dst);
}
/*========================================================================*/
int SceneMultipick(Multipick *smp)
{
  CScene *I=&Scene;
  if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
    SceneRender(NULL,0,0,NULL); /* remove overlay if present */
  SceneDontCopyNext();
  if(I->StereoMode>1) {
    smp->x = smp->x % (I->Width/2);
  }
  SceneRender(NULL,0,0,smp);
  SceneDirty();
  return(1);
}
/*========================================================================*/
int SceneGetNFrame(void)
{
  CScene *I=&Scene;
  return(I->NFrame);
}
/*========================================================================*/
void SceneGetView(SceneViewType view)
{
  float *p;
  int a;
  CScene *I=&Scene;
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
  *(p++) = SettingGet(cSetting_ortho);
}
/*========================================================================*/
void SceneSetView(SceneViewType view,int quiet)
{
  float *p;
  int a;
  CScene *I=&Scene;
  p=view;
  for(a=0;a<16;a++)
    I->RotMatrix[a] = *(p++); 
  I->Pos[0] = *(p++);
  I->Pos[1] = *(p++);
  I->Pos[2] = *(p++);
  I->Origin[0] = *(p++);
  I->Origin[1] = *(p++);
  I->Origin[2] = *(p++);
  SceneClipSet(p[0],p[1]);
  p+=2;
  SettingSet(cSetting_ortho,*(p++));
  if(!quiet) { 
    PRINTFB(FB_Scene,FB_Actions)
      " Scene: view updated.\n"
      ENDFB;
  }
  SceneRovingDirty();
}
/*========================================================================*/
void SceneDontCopyNext(void)
/* disables automatic copying of the image for the next rendering run */
{
  CScene *I=&Scene;
  I->CopyNextFlag=false;
}
/*========================================================================*/
void SceneUpdateStereoMode(void)
{
  CScene *I=&Scene;
  if(I->StereoMode) {
    SceneSetStereo(true);
  }
}
/*========================================================================*/
void SceneSetStereo(int flag)
{
  CScene *I=&Scene;
  if(flag) 
    I->StereoMode=(int)SettingGet(cSetting_stereo_mode);
  else
    I->StereoMode=false;
  SceneDirty();
}
/*========================================================================*/
void SceneTranslate(float x,float y, float z)
{
  CScene *I=&Scene;
  I->Pos[0]+=x;
  I->Pos[1]+=y;
  I->Pos[2]+=z;
  I->Back-=z;
  I->Front-=z;
  if(I->Front>I->Back)
	 I->Front=I->Back+cSliceMin;
  if(I->Front<cFrontMin) I->Front=cFrontMin;
  I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
  SceneDirty();
}
/*========================================================================*/
void SceneClipSet(float front,float back)
{
  CScene *I=&Scene;
  I->Front=front;
  I->Back=back;
  if(I->Front>I->Back)
	 I->Front=I->Back+cSliceMin;
  if(I->Front<cFrontMin) I->Front=cFrontMin;
  I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
  SceneDirty();
}
/*========================================================================*/
void SceneClip(int plane,float movement,char *sele,int state) /* 0=front, 1=back*/
{
  CScene *I=&Scene;
  float avg;
  float mn[3],mx[3],cent[3],v0[3],offset[3],origin[3];
  switch(plane) {
  case 0: /* near */
    SceneClipSet(I->Front-movement,I->Back);
    break;
  case 1: /* far */
    SceneClipSet(I->Front,I->Back-movement);      
    break;
  case 2: /* move */
    SceneClipSet(I->Front-movement,I->Back-movement);    
    break;
  case 3: /* slab */
    if(sele[0]) {
      if(!ExecutiveGetExtent(sele,mn,mx,true,state,false))
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
    SceneClipSet(avg-movement,avg+movement);
    break;
  case 4: /* atoms */
    if(!sele) 
      sele=cKeywordAll;
    else if(!sele[0]) {
      sele=cKeywordAll;
    } 
    if(!ExecutiveGetCameraExtent(sele,mn,mx,true,state))
      sele = NULL;
    if(sele) {
      if(sele[0]) {
        average3f(mn,mx,cent); /* get center of selection */
        MatrixTransform3f(I->RotMatrix,I->Origin,origin); /* convert to view-space */
        subtract3f(mx,origin,mx); /* how far from origin? */
        subtract3f(mn,origin,mn); /* how far from origin? */
        SceneClipSet(-I->Pos[2]-mx[2]-movement,-I->Pos[2]-mn[2]+movement);
      } else {
        sele = NULL;
      }
    }
    break;
  }
}
/*========================================================================*/
void SceneSetMatrix(float *m)
{
  CScene *I=&Scene;
  int a;
  for(a=0;a<16;a++)
	 I->RotMatrix[a]=m[a];
}
/*========================================================================*/
void SceneGetViewNormal(float *v)
{
  CScene *I=&Scene;
  copy3f(I->ViewNormal,v);
}
/*========================================================================*/
int SceneGetState(void)
{
  return(SettingGetGlobal_i(cSetting_state)-1);
}
/*========================================================================*/
float *SceneGetMatrix()
{
  CScene *I=&Scene;
  return(I->RotMatrix);
}
/*========================================================================*/
void ScenePNG(char *png)
{
  CScene *I=&Scene;
  unsigned int buffer_size;
  GLvoid *image;
  int reset_alpha = false;

  buffer_size = 4*I->Width*I->Height;
  if(!I->CopyFlag) {
	 image = (GLvoid*)Alloc(char,buffer_size);
	 ErrChkPtr(image);
    if(PMGUI) {
      glReadBuffer(GL_BACK);
      glReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                   GL_RGBA,GL_UNSIGNED_BYTE,image);
      
      reset_alpha = true;
      I->ImageBufferHeight=I->Height;
      I->ImageBufferWidth=I->Width;
    } else {
       PRINTFB(FB_Scene,FB_Errors)
         " ScenePNG-WARNING: writing a blank image buffer.\n"
         ENDFB;
     }
  } else {
    PRINTFB(FB_Scene,FB_Blather)
      " ScenePNG: writing cached image.\n"
      ENDFB;
    image=I->ImageBuffer;
    reset_alpha = I->CopiedFromOpenGL;
  }
  if(reset_alpha&&image) {
    char *p = image;
    int x,y;
    for(y=0;y<I->Height;y++) {
      for(x=0;x<I->Width;x++) {
        p[3]=0xFF;
        p+=4;
      }
    }
  }
  if(MyPNGWrite(png,image,I->ImageBufferWidth,I->ImageBufferHeight)) {
    PRINTFB(FB_Scene,FB_Actions) 
      " ScenePNG: wrote %dx%d pixel image to file \"%s\".\n",
      I->ImageBufferWidth,I->ImageBufferHeight,png
      ENDFB;
  } else {
    PRINTFB(FB_Scene,FB_Errors) 
      " ScenePNG-Error: error writing \"%s\"! Please check directory...\n",
      png
      ENDFB;
  }
  
  if(!I->CopyFlag)
	 FreeP(image);
}
/*========================================================================*/
void ScenePerspective(int flag)
{
  float persp;
  persp=(float)(!flag);
  SettingSetfv(cSetting_ortho,&persp);
  SceneDirty();
}
/*========================================================================*/
int SceneGetFrame(void)
{
  return(SettingGetGlobal_i(cSetting_frame)-1);
}
/*========================================================================*/
void SceneCountFrames() 
{
  CScene *I=&Scene;
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
  mov_len = MovieGetLength();
  if(mov_len>0) {
    I->NFrame=mov_len;
  } else if(mov_len<0) {
    mov_len=-mov_len;
    if(I->NFrame<mov_len) /* allows you to see cached movie even w/o object */
      I->NFrame=mov_len;
  }
  PRINTFD(FB_Scene)
    " SceneCountFrames: leaving... I->NFrame %d\n",I->NFrame
    ENDFD
}
/*========================================================================*/
void SceneSetFrame(int mode,int frame)
{
  CScene *I=&Scene;
  int newFrame=0;
  int newState=0;
  int movieCommand = false;
  newFrame = SettingGetGlobal_i(cSetting_frame) -1;
  PRINTFD(FB_Scene)
    " SceneSetFrame: entered.\n"
    ENDFD;
  switch(mode) {
  case -1: /* movie/frame override - go to this state absolutely! */
    newState=frame;
    break;
  case 0: /* absolute */
    newFrame=frame; 
	 break;
  case 1: /* relative */
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
  SceneCountFrames();
  if (mode>=0) { 
    if(newFrame>=I->NFrame) newFrame=I->NFrame-1;
    if(newFrame<0) newFrame=0;
    newState = MovieFrameToIndex(newFrame);
    if(newFrame==0) {
      MovieMatrix(cMovieMatrixRecall);
    }
    if(movieCommand) {
      MovieDoFrameCommand(newFrame);
    }
    if(SettingGet(cSetting_cache_frames))
      I->MovieFrameFlag=true;
  }
  SettingSetGlobal_i(cSetting_frame,newFrame+1);
  SettingSetGlobal_i(cSetting_state,newState+1);
  SceneDirty();
  PRINTFD(FB_Scene)
    " SceneSetFrame: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void ScenePurgeCopy(void)
{
  CScene *I=&Scene;
  I->CopyFlag=false;
  if(I->MovieOwnsImageFlag) 
	 {
		I->MovieOwnsImageFlag=false;
		I->ImageBuffer=NULL;
	 }
}
/*========================================================================*/
void SceneDirty(void) 
	  /* This means that the current image on the screen (and/or in the buffer)
		 needs to be updated */
{
  CScene *I=&Scene;

  PRINTFD(FB_Scene)
    " SceneDirty: called.\n"
    ENDFD;

  I->DirtyFlag=true;
  ScenePurgeCopy();
  OrthoDirty();
}

void SceneRovingPostpone(void)
{
  CScene *I=&Scene;
  float delay;
  if(SettingGet(cSetting_roving_detail)) {
    delay = SettingGet(cSetting_roving_delay);
    if(delay<0.0F) {
      I->RovingLastUpdate = UtilGetSeconds(); /* put off delay */
    }
  }
}

void SceneRovingDirty(void)
{
  CScene *I=&Scene;

  if(SettingGet(cSetting_roving_detail)) {
    SceneRovingPostpone();
    I->RovingDirtyFlag=true;
  }
}

/*========================================================================*/
void SceneChanged(void)
{
  CScene *I=&Scene;
  I->ChangedFlag=true;
  SceneDirty();
}
/*========================================================================*/
Block *SceneGetBlock(void)
{
  CScene *I=&Scene;
  return(I->Block);
}
/*========================================================================*/
void SceneMakeMovieImage(void) {
  CScene *I=&Scene;
  float *v;

  I->DirtyFlag=false;
  if(SettingGet(cSetting_ray_trace_frames)) {
	SceneRay(0,0,(int)SettingGet(cSetting_ray_default_renderer),NULL,NULL,
            0.0F,0.0F,false); 
  } else {
	 v=SettingGetfv(cSetting_bg_rgb);
    if(PMGUI) {
      glDrawBuffer(GL_BACK);
      glClearColor(v[0],v[1],v[2],1.0);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glClearColor(0.0,0.0,0.0,1.0);
      SceneRender(NULL,0,0,NULL);
      SceneCopy(GL_BACK,true);
    }
  }
  if(I->ImageBuffer&&(I->ImageBufferHeight==I->Height)&&(I->ImageBufferWidth==I->Width)) {
	 MovieSetImage(MovieFrameToImage(SettingGetGlobal_i(cSetting_frame)-1)
                                    ,I->ImageBuffer);
    I->MovieOwnsImageFlag=true;
  } else {
    I->MovieOwnsImageFlag=false;
  }
  I->CopyFlag=true;
}
/*========================================================================*/
void SceneIdle(void)
{
  CScene *I=&Scene;
  double renderTime;
  double minTime;
  int frameFlag = false;
  int rockFlag = false;
  float ang_cur,disp,diff;

  if(MoviePlaying())
    {
		renderTime = -I->LastFrameTime + UtilGetSeconds();
		minTime=SettingGet(cSetting_movie_delay)/1000.0;
		if(renderTime>=minTime) {
        frameFlag=true;
        rockFlag=true;
      }
    }
  if(Control.Rocking&&(!rockFlag))
    {
		renderTime = -I->LastRockTime + UtilGetSeconds();
		minTime=SettingGet(cSetting_rock_delay)/1000.0;
		if(renderTime>=minTime) {
        rockFlag=true;
        I->LastRockTime=UtilGetSeconds();
      }
    }
  if(Control.Rocking&&rockFlag) {


	I->RockTime+=I->RenderTime;
    ang_cur = (float)(I->RockTime*SettingGet(cSetting_sweep_speed));
    
    disp = (float)(SettingGet(cSetting_sweep_angle)*(3.1415/180.0)*sin(ang_cur)/2);
    diff = (float)(disp-I->LastRock);
    SceneRotate((float)(180*diff/cPI),0.0F,1.0F,0.0F);
    I->LastRock = disp;
  }
  if(MoviePlaying()&&frameFlag)
	 {
      I->LastFrameTime = UtilGetSeconds();
      if((SettingGetGlobal_i(cSetting_frame)-1)==(I->NFrame-1)) {
        if((int)SettingGet(cSetting_movie_loop)) {
          SceneSetFrame(7,0);
        } else
          MoviePlay(cMovieStop);
      } else 
        SceneSetFrame(5,1);
	 }
}
/*========================================================================*/
void SceneWindowSphere(float *location,float radius)
{
  CScene *I=&Scene;
  float v0[3];
  float dist;
  float aspRat = ((float) I->Width) / ((float) I->Height);
  float fov;

  /* find where this point is in relationship to the origin */
  subtract3f(I->Origin,location,v0); 

  dist = I->Pos[2];
  /*  printf("%8.3f %8.3f %8.3f\n",I->Front,I->Pos[2],I->Back);*/

  MatrixTransform3f(I->RotMatrix,v0,I->Pos); /* convert to view-space */
  fov = SettingGet(cSetting_field_of_view);
  if(aspRat<1.0)
    fov *= aspRat;

  dist = (float)(radius/tan((fov/2.0)*cPI/180.0));

  I->Pos[2]-=dist;
  I->Front=(-I->Pos[2]-radius*1.2F);
  I->FrontSafe=(I->Front<cFrontMin ? cFrontMin : I->Front);  
  I->Back=(-I->Pos[2]+radius*1.55F);

  SceneRovingDirty();
  /*printf("%8.3f %8.3f %8.3f\n",I->Front,I->Pos[2],I->Back);*/
}
/*========================================================================*/
void SceneRelocate(float *location)
{
  CScene *I=&Scene;
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
  I->Front=(-I->Pos[2]-(slab_width*0.45F));
  I->FrontSafe=(I->Front<cFrontMin ? cFrontMin : I->Front);  
  I->Back=(-I->Pos[2]+(slab_width*0.55F));
  SceneRovingDirty();

}
/*========================================================================*/
void SceneOriginGet(float *origin)
{
  CScene *I=&Scene;
  copy3f(I->Origin,origin);
}
/*========================================================================*/
void SceneOriginSet(float *origin,int preserve)
{
  CScene *I=&Scene;
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
  SceneDirty();
}
/*========================================================================*/
void SceneObjectAdd(CObject *obj)
{
  CScene *I=&Scene;
  ObjRec *rec = NULL;
  ListElemAlloc(rec,ObjRec);
  rec->next=NULL;
  obj->Enabled=true;
  rec->obj=obj;
  ListAppend(I->Obj,rec,next,ObjList);
  SceneCountFrames();
  SceneChanged();
}
/*========================================================================*/
void SceneObjectDel(CObject *obj)
{
  CScene *I=&Scene;
  ObjRec *rec = NULL;

  if(!obj) {
    while(ListIterate(I->Obj,rec,next)) {
      if(rec) {
        ListDetach(I->Obj,rec,next,ObjList);
        ListElemFree(rec);
      }
    }
  } else {
    while(ListIterate(I->Obj,rec,next))
      if(rec->obj==obj)
        break;
    if(rec) {
      rec->obj->Enabled=false;
      ListDetach(I->Obj,rec,next,ObjList);
      ListElemFree(rec);
    }
  }
  SceneCountFrames();
  SceneDirty();
}
/*========================================================================*/
int SceneLoadPNG(char *fname,int movie_flag,int quiet) 
{
  CScene *I=&Scene;
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
      PRINTFB(FB_Scene,FB_Details)
        " Scene: loaded image from '%s'.\n",fname
        ENDFB;
    }
    I->CopyFlag=true;
    I->CopiedFromOpenGL = false;
    OrthoRemoveSplash();
    SettingSet(cSetting_text,0.0);
    if(movie_flag&&I->ImageBuffer&&(I->ImageBufferHeight==I->Height)&&(I->ImageBufferWidth==I->Width)) {
      MovieSetImage(MovieFrameToImage(SettingGetGlobal_i(cSetting_frame)-1)
                    ,I->ImageBuffer);
      I->MovieOwnsImageFlag=true;
      I->MovieFrameFlag=true;
    } else {
      I->MovieOwnsImageFlag=false;
      I->DirtyFlag=false; /* make sure we don't overwrite image */
    }
    OrthoDirty();
    ok=true;
  } else {
    if(!quiet) {
      PRINTFB(FB_Scene,FB_Errors)
        " Scene: unable to load image from '%s'.\n",fname
        ENDFB;
    }
  }
  return(ok);
}
/*========================================================================*/
void SceneDraw(Block *block)
{
  CScene *I=&Scene;
  int overlay,text;
  int width,height;
  int double_pump;

  if(PMGUI) {
    overlay = (int)SettingGet(cSetting_overlay);
    text = (int)SettingGet(cSetting_text);
    double_pump = (int)SettingGet(cSetting_stereo_double_pump_mono);

    if(overlay||(!text)) 

      if(I->CopyFlag)
        {
          glReadBuffer(GL_BACK); 

          if(I->ImageBufferHeight>I->Height||I->ImageBufferWidth>I->Width) {
            glColor3f(1.0F,0.2F,0.2F);
            GrapDrawStr("Sorry, I can't display an oversize image.",30,60);
            GrapDrawStr("To save image, use File Menu or enter \"png <filename>\".",30,40);
          } else {
            width = I->ImageBufferWidth;
            height = I->ImageBufferHeight;
            
            if((width<I->Width)||(height<I->Height)) {
              glRasterPos3i((int)((I->Width-width)/2+I->Block->rect.left),
                            (int)((I->Height-height)/2+I->Block->rect.bottom),0);
            } else {
              glRasterPos3i(I->Block->rect.left,I->Block->rect.bottom,0);
            }
            if(I->ImageBuffer) {
#if 1
              glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
#else
              if(!(double_pump||(I->StereoMode==1))) {
                glDrawBuffer(GL_BACK);
                glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
              } else {
                glDrawBuffer(GL_BACK_LEFT);
                glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
                glDrawBuffer(GL_BACK_RIGHT);
                glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
              }
#endif
            }

          }
          I->RenderTime = -I->LastRender;
          I->LastRender = UtilGetSeconds();
          I->RenderTime += I->LastRender;
          ButModeSetRate((float)I->RenderTime);
        }
    
    glColor3f(1.0,1.0,1.0);
  }
}
/*========================================================================*/

typedef unsigned char pix[4];
#define cRange 5
/*typedef pix pix_array[cRange*2+1][cRange*2+1];*/

unsigned int SceneFindTriplet(int x,int y,GLenum gl_buffer) 
{
  int result = 0;
  /*int before_check[100];
  int *int_ptr;
*/
  pix buffer[cRange*2+1][cRange*2+1];
 /*int after_check[100];*/
  /* pix_array *array_ptr;
  char *safe_place;
*/
  int a,b,d,flag;
  int debug = false;
  unsigned char *c;
  int strict = false;
  GLint rb,gb,bb;
  int bkrd_alpha = 0xFF;
  int check_alpha = false;

  if(PMGUI) { /*just in case*/
  
	glGetIntegerv(GL_RED_BITS,&rb);
	glGetIntegerv(GL_GREEN_BITS,&gb);
	glGetIntegerv(GL_BLUE_BITS,&bb);

	if((rb>=8)&&(gb>=8)&&(bb>=8))
		strict = true;

    if(Feedback(FB_Scene,FB_Debugging)) debug=true;
    
    glReadBuffer(gl_buffer);

/*	safe_place = (char*)malloc(1000000);
	array_ptr = (pix_array*)(safe_place+500000);*/

    /*int_ptr = before_check;
    for(a=0;a<100;a++) *(int_ptr++)=0x12345678;
    int_ptr = after_check;
    for(a=0;a<100;a++) *(int_ptr++)=0x12345678;
*/
/*printf("%d %d %d %d %p %p\n",x-cRange,y-cRange,cRange*2+1,cRange*2+1,buffer,&buffer[0][0][0]);*/
    glReadPixels(x-cRange,y-cRange,cRange*2+1,cRange*2+1,GL_RGBA,GL_UNSIGNED_BYTE,&buffer[0][0][0]);

/*	{
	for(a=0;a<=(cRange*2);a++)
       for(b=0;b<=(cRange*2);b++)
		  for(d=0;d<4;d++)
            buffer[a][b][d]=(*array_ptr)[a][b][d];
		 
         
    }
	free(safe_place);*/

      /*int_ptr = before_check;
    for(a=0;a<100;a++) {
      if(*(int_ptr++)!=0x12345678) {
        printf(" SceneFindTriplet-WARNING: OpenGL glReadPixels may have damaged stack\n");
      }
    }
    int_ptr = after_check;
    for(a=0;a<100;a++) {
      if(*(int_ptr++)!=0x12345678) {
        printf(" SceneFindTriplet-WARNING: OpenGL glReadPixels may have damaged stack\n");
      }
    }*/

	  if(debug) {
      for(a=0;a<=(cRange*2);a++)
        {
          for(b=0;b<=(cRange*2);b++)
            printf("%2x ",(buffer[a][b][0]+buffer[a][b][1]+buffer[a][b][2])&0xFF);
          printf("\n");
        }
      printf("\n");	 
      for(a=0;a<=(cRange*2);a++)
        {
          for(b=0;b<=(cRange*2);b++)
            printf("%02x ",(buffer[a][b][3])&0xFF);
          printf("\n");
        }
      printf("\n");	 
       for(a=0;a<=(cRange*2);a++)
        {
          for(b=0;b<=(cRange*2);b++)
            printf("%02x%02x%02x ",(buffer[a][b][0])&0xFF,(buffer[a][b][1])&0xFF,(buffer[a][b][2])&0xFF);
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
           c = &buffer[a+cRange][b+cRange][0];
           if(c[4]==bkrd_alpha) {
             check_alpha = true;
             flag=false;
           }
         }

     /* now find the correct pixel */

     flag=true;
     for(d=0;flag&&(d<cRange);d++)
       for(a=-d;flag&&(a<=d);a++)
         for(b=-d;flag&&(b<=d);b++) {
           c = &buffer[a+cRange][b+cRange][0];
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
  }
  return(result);
}
/*========================================================================*/
unsigned int *SceneReadTriplets(int x,int y,int w,int h,GLenum gl_buffer)
{ 
  unsigned int *result = NULL;
  pix *buffer=NULL;
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
  if(PMGUI) { /*just in case*/
    
    
    glGetIntegerv(GL_RED_BITS,&rb);
    glGetIntegerv(GL_RED_BITS,&gb);
    glGetIntegerv(GL_RED_BITS,&bb);
    
    if((rb>=8)&&(gb>=8)&&(bb>=8))
		strict = true;
    
    buffer=Alloc(pix,w*h);
    result = VLAlloc(unsigned int,w*h);
    glReadBuffer(gl_buffer);
    glReadPixels(x,y,w,h,GL_RGBA,GL_UNSIGNED_BYTE,&buffer[0][0]);
    
     /* first, check to make sure bkrd_alpha is correct 
        (this is a bug for systems with broken alpha, such as Extreme 3D on Solaris 8 */

    for(a=0;a<w;a++)
      for(b=0;b<h;b++)
        {
          c = &buffer[a+b*w][0];
          if(c[4]==bkrd_alpha) {
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
    FreeP(buffer);
    VLASize(result,unsigned int,cc);
  }
  return(result);

}
/*========================================================================*/
int SceneRelease(Block *block,int button,int x,int y,int mod) 
{
  CScene *I=&Scene;
  ObjectMolecule *obj;
  if(I->SculptingFlag) {
    /* SettingSet(cSetting_sculpting,1); */
    obj=(ObjectMolecule*)I->LastPicked.ptr;
    obj->AtomInfo[I->LastPicked.index].protekted=I->SculptingSave;
    I->SculptingFlag=0;
  }
  return(1);
}
/*========================================================================*/
int SceneClick(Block *block,int button,int x,int y,int mod)
{
  CScene *I=&Scene;
  CObject *obj;
  ObjectMolecule *objMol;
  OrthoLineType buffer,buf1,buf2;
  WordType selName = "";
  int mode;
  int atIndex;
  mode = ButModeTranslate(button,mod);
  I->Button=button;    
  I->SculptingSave = 0;

  switch(mode) {
  case cButModeRectAdd:
  case cButModeRectSub:
  case cButModeRect:
    return(0);
    break;
  case cButModeRotXYZ:
  case cButModeTransXY:
  case cButModeTransZ:
  case cButModeClipNF:
  case cButModeClipN:    
  case cButModeClipF:    
  case cButModeRotZ:
    SceneDontCopyNext();

    y=y-I->Block->margin.bottom;
    x=x-I->Block->margin.left;
    
    if(I->StereoMode>1)
      x = x % (I->Width/2);

    I->LastX=x;
    I->LastY=y;	 

    SceneDirty();
    break;
  case cButModePickAtom:
    if(I->StereoMode>1)
      x = x % (I->Width/2);

    if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
      SceneRender(NULL,0,0,NULL); /* remove overlay if present */
    SceneDontCopyNext();
    I->LastPicked.ptr = NULL;
	 SceneRender(&I->LastPicked,x,y,NULL);
	 if(I->LastPicked.ptr) { /* did we pick something? */
		obj=(CObject*)I->LastPicked.ptr;
      y=y-I->Block->margin.bottom;
      x=x-I->Block->margin.left;
      I->LastX=x;
      I->LastY=y;	
      switch(obj->type) {
      case cObjectMolecule:
        if(Feedback(FB_ObjectMolecule,FB_Results)) {
          if(obj->fDescribeElement)
            obj->fDescribeElement(obj,I->LastPicked.index,buffer);
          PRINTF " You clicked %s -> (%s)",buffer,cEditorSele1 ENDF;
          if(SettingGet(cSetting_logging)) {
            objMol = (ObjectMolecule*)obj;            
            ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buffer);
            sprintf(buf2,"cmd.edit(\"%s\",pkresi=1)",buffer);
            PLog(buf2,cPLog_pym);
          }
          OrthoRestorePrompt();
        }
        sprintf(buffer,"%s`%d",
                obj->Name,I->LastPicked.index+1);    
        SelectorCreate(cEditorSele1,buffer,NULL,true,NULL);
        ExecutiveDelete(cEditorSele2);
        EditorSetActiveObject((ObjectMolecule*)obj,
                              SettingGetGlobal_i(cSetting_state)-1);
        if(EditorActive()) {
          SelectorCreate(cEditorRes,"(byres pk1)",NULL,true,NULL);
          if(SettingGet(cSetting_auto_hide_selections))
            ExecutiveHideSelections();
        }
        
        WizardDoPick(0);
        break;
      case cObjectGadget:
        break;
      default:
        EditorSetActiveObject(NULL,0);
        break;
      }
    } else {
      EditorSetActiveObject(NULL,0);
    }
    SceneDirty();
    break;
  case cButModePickBond:
  case cButModePkTorBnd:
    if(I->StereoMode>1)
      x = x % (I->Width/2);

    if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
      SceneRender(NULL,0,0,NULL); /* remove overlay if present */
    SceneDontCopyNext();
    I->LastPicked.ptr = NULL;
	 SceneRender(&I->LastPicked,x,y,NULL);
	 if(I->LastPicked.ptr) {
		obj=(CObject*)I->LastPicked.ptr;
      y=y-I->Block->margin.bottom;
      x=x-I->Block->margin.left;
      I->LastX=x;
      I->LastY=y;	

      switch(obj->type) {
      case cObjectMolecule:
        if(Feedback(FB_ObjectMolecule,FB_Results)) {
          if(obj->fDescribeElement)
            obj->fDescribeElement(obj,I->LastPicked.index,buffer);
          PRINTF " You clicked %s -> (%s)",buffer,cEditorSele1 ENDF;
          OrthoRestorePrompt();
        }
		  sprintf(buffer,"%s`%d",
					 obj->Name,I->LastPicked.index+1);    
        SelectorCreate(cEditorSele1,buffer,NULL,true,NULL);
        objMol = (ObjectMolecule*)obj;
        if(I->LastPicked.bond>=0) {
          atIndex = objMol->Bond[I->LastPicked.bond].index[0];
          if(atIndex == I->LastPicked.index)
            atIndex = objMol->Bond[I->LastPicked.bond].index[1];              
          if(Feedback(FB_ObjectMolecule,FB_Results)) {
            if(obj->fDescribeElement)
              obj->fDescribeElement(obj,atIndex,buffer);
            PRINTF " You clicked %s -> (%s)",buffer,cEditorSele2 ENDF;
            OrthoRestorePrompt();
          }
          
          if(SettingGet(cSetting_logging)) {
            objMol = (ObjectMolecule*)obj;            
            ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1);
            ObjectMoleculeGetAtomSeleLog(objMol,atIndex,buf2);
            sprintf(buffer,"cmd.edit(\"%s\",\"%s\")",buf1,buf2);
            PLog(buffer,cPLog_pym);
          }
          sprintf(buffer,"%s`%d",
                  obj->Name,atIndex+1);    
          SelectorCreate(cEditorSele2,buffer,NULL,true,NULL);
          EditorSetActiveObject(objMol,
                                SettingGetGlobal_i(cSetting_state)-1);


          WizardDoPick(1);

          if(mode==cButModePkTorBnd) {
            /* get ready to drag */
            SceneDontCopyNext();
            switch(obj->type) {
            case cObjectMolecule:
              objMol = (ObjectMolecule*)obj;
              EditorPrepareDrag(objMol,I->LastPicked.index,
                                SettingGetGlobal_i(cSetting_state)-1);
              I->SculptingFlag = 1;
              I->SculptingSave =  objMol->AtomInfo[I->LastPicked.index].protekted;
              objMol->AtomInfo[I->LastPicked.index].protekted=2;
              break;
            }
          }
        }
        break;
      case cObjectGadget:
        break;
      default:
        EditorSetActiveObject(NULL,0);
        break;
      }
    } else {
      EditorSetActiveObject(NULL,0);
    }
    SceneDirty();
    break;
  case cButModeMovFrag:
  case cButModeTorFrag:
  case cButModeRotFrag:

    if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
      SceneRender(NULL,0,0,NULL); /* remove overlay if present */
    SceneDontCopyNext();
    if(I->StereoMode>1)
      x = x % (I->Width/2);
	 SceneRender(&I->LastPicked,x,y,NULL);
	 if(I->LastPicked.ptr) {
      obj=(CObject*)I->LastPicked.ptr;
      y=y-I->Block->margin.bottom;
      x=x-I->Block->margin.left;
      I->LastX=x;
      I->LastY=y;	
      switch(obj->type) {
      case cObjectMolecule:
        
        if(Feedback(FB_ObjectMolecule,FB_Results)) {
          if(obj->fDescribeElement) 
            obj->fDescribeElement(obj,I->LastPicked.index,buffer);
          PRINTF " You clicked %s",buffer ENDF;        
          OrthoRestorePrompt();
        }
        objMol = (ObjectMolecule*)obj;
        EditorPrepareDrag(objMol,I->LastPicked.index,
                          SettingGetGlobal_i(cSetting_state)-1);
        I->SculptingFlag = 1;
        I->SculptingSave =  objMol->AtomInfo[I->LastPicked.index].protekted;
        objMol->AtomInfo[I->LastPicked.index].protekted=2;
        break;
      case cObjectGadget:
        break;
      default:
        EditorSetActiveObject(NULL,0);
        break;
      }
      /*
        (int)SettingGet(cSetting_sculpting);
            SettingSet(cSetting_sculpting,0);*/
    }
    break;
 
  case cButModePk1:
  case cButModePk2:
  case cButModePk3:
  case cButModeAddToPk1:
  case cButModeAddToPk2:
  case cButModeAddToPk3:
  case cButModeOrigAt:
  case cButModeCent:
    if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
      SceneRender(NULL,0,0,NULL); /* remove overlay if present */
    SceneDontCopyNext();

    if(I->StereoMode>1)
      x = x % (I->Width/2);

	 SceneRender(&I->LastPicked,x,y,NULL);
	 if(I->LastPicked.ptr) {
		obj=(CObject*)I->LastPicked.ptr;

      switch(obj->type) {
      case cObjectMolecule:

        
        if(Feedback(FB_ObjectMolecule,FB_Results)) {
          if(obj->fDescribeElement) 
            obj->fDescribeElement(obj,I->LastPicked.index,buffer);
          PRINTF " You clicked %s",buffer ENDF;        
          OrthoRestorePrompt();
        }
        sprintf(buffer,"%s`%d",
                obj->Name,I->LastPicked.index+1);
        switch(mode) {
        case cButModePk1:
        case cButModeAddToPk1:
          strcpy(selName,"lb");
          break;
        case cButModePk2:
        case cButModeAddToPk2:
          strcpy(selName,"mb");
          break;
        case cButModePk3:
        case cButModeAddToPk3:
          strcpy(selName,"rb");
          break;
        case cButModeOrigAt:
          sprintf(buf2,"origin (%s)",buffer);        
          OrthoCommandIn(buf2);
          if(obj->type==cObjectMolecule) {
            if(SettingGet(cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1);
              sprintf(buffer,"cmd.origin(\"%s\")",buf1);
              PLog(buffer,cPLog_pym);
            }
          }
          PRINTFB(FB_Scene,FB_Actions) 
            " Scene: Origin set.\n"
            ENDFB;
          break;
        case cButModeCent:
          sprintf(buf2,"center (%s),state=-1",buffer);        
          OrthoCommandIn(buf2);
          if(obj->type==cObjectMolecule) {
            if(SettingGet(cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1);
              sprintf(buffer,"cmd.center(\"%s\",state=-1)",buf1);
              PLog(buffer,cPLog_pym);
            }
          }
          PRINTFB(FB_Scene,FB_Actions) 
            " Scene: Centered.\n"
            ENDFB;
          break;
        }
        switch(mode) {
        case cButModePk1:
        case cButModePk2:
        case cButModePk3:
          SelectorCreate(selName,buffer,NULL,false,NULL);
          if(SettingGet(cSetting_auto_hide_selections))
            ExecutiveHideSelections();
          if(SettingGet(cSetting_auto_show_selections))
            ExecutiveSetObjVisib(selName,1);
          if(obj->type==cObjectMolecule) {
            if(SettingGet(cSetting_logging)) {
              objMol = (ObjectMolecule*)obj;            
              ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1);
              sprintf(buffer,"cmd.select('%s',\"%s\")",selName,buf1);
              PLog(buffer,cPLog_pym);
            }
          }
          WizardDoSelect(selName);
          break;
        case cButModeAddToPk1:
        case cButModeAddToPk2:
        case cButModeAddToPk3:
          if(SelectorIndexByName(selName)>=0) {
            sprintf(buf2,"( ((%s) or (%s)) and not ((%s) in (%s)))",
                    selName,buffer,buffer,selName);
            SelectorCreate(selName,buf2,NULL,false,NULL);
            if(obj->type==cObjectMolecule) {
              if(SettingGet(cSetting_logging)) {
                objMol = (ObjectMolecule*)obj;            
                ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buffer);
                sprintf(buf2,"( ((%s) or (%s)) and not ((%s) in (%s)))",
                        selName,buffer,buffer,selName);
                sprintf(buffer,"cmd.select('%s',\"%s\")",selName,buf2);
                PLog(buffer,cPLog_pym);
              }
            }
          } else {
            SelectorCreate(selName,buffer,NULL,false,NULL);
            if(obj->type==cObjectMolecule) {
              if(SettingGet(cSetting_logging)) {
                objMol = (ObjectMolecule*)obj;            
                ObjectMoleculeGetAtomSeleLog(objMol,I->LastPicked.index,buf1);
                sprintf(buffer,"cmd.select('%s',\"%s\")",selName,buf1);
                PLog(buffer,cPLog_pym);
              }
            }
          }
          if(SettingGet(cSetting_auto_hide_selections))
            ExecutiveHideSelections();
          if(SettingGet(cSetting_auto_show_selections))
            ExecutiveSetObjVisib(selName,1);
          WizardDoSelect(selName);
          break;
        }
      case cObjectGadget:
        break;
      default:
        EditorSetActiveObject(NULL,0);
        break;
      }
	 } else {
      PRINTFB(FB_Scene,FB_Warnings) 
        " SceneClick: no atom found nearby.\n"
        ENDFB;
		OrthoRestorePrompt();
	 }
  }
  return(1);
}
/*========================================================================*/
float SceneGetScreenVertexScale(float *v1)
{
  /* get conversion factor from screen point to atomic coodinate */
  CScene *I=&Scene;
  float vl,p1[4],p2[4];
  /* now, scale properly given the current projection matrix */
  copy3f(v1,p1);
  p1[3] = 1.0;
  MatrixTransform44f4f(I->ModMatrix,p1,p2); /* modelview transformation */
  copy4f(p2,p1);
  p2[0]+=1.0;
  MatrixTransform44f4f(I->ProMatrix,p1,p1); /* projection transformation */
  MatrixTransform44f4f(I->ProMatrix,p2,p2);
  p1[0]=p1[0]/p1[3];/* perspective vision */
  p1[1]=p1[1]/p1[3];
  p1[2]=0.0;
  p2[0]=p2[0]/p2[3];
  p2[1]=p2[1]/p2[3];
  p2[2]=0.0;
  p1[0]=(p1[0]+1.0F)*(I->Width/2.0F); /* viewport transformation */
  p1[1]=(p1[1]+1.0F)*(I->Height/2.0F);
  p2[0]=(p2[0]+1.0F)*(I->Width/2.0F);
  p2[1]=(p2[1]+1.0F)*(I->Height/2.0F);
  vl=(float)diff3f(p1,p2);
  if(vl<R_SMALL4)
    vl=100.0F;
  
  return(1.0F/vl);
}

void SceneRovingChanged(void)
{
  CScene *I=&Scene;  
  SceneRovingDirty();
  I->RovingCleanupFlag=true;
}

static void SceneRovingCleanup(void)
{
  CScene *I=&Scene;  
  char *s;  
  char buffer[OrthoLineLength];

  I->RovingCleanupFlag=false;

  s = SettingGet_s(NULL,NULL,cSetting_roving_selection);

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

void SceneRovingUpdate(void)
{
  CScene *I=&Scene;
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
                          (UtilGetSeconds()-I->RovingLastUpdate)>
                          fabs(SettingGet(cSetting_roving_delay)))) {
    
    if(I->RovingCleanupFlag)
      SceneRovingCleanup();
    
    s = SettingGet_s(NULL,NULL,cSetting_roving_selection);
    sticks = SettingGet(cSetting_roving_sticks);
    lines = SettingGet(cSetting_roving_lines);
    labels = SettingGet(cSetting_roving_labels);
    spheres = SettingGet(cSetting_roving_spheres);
    ribbon = SettingGet(cSetting_roving_ribbon);
    cartoon = SettingGet(cSetting_roving_cartoon);
    polar_contacts = SettingGet(cSetting_roving_polar_contacts);
    polar_cutoff = SettingGet(cSetting_roving_polar_cutoff);
    nonbonded = SettingGet(cSetting_roving_nonbonded);
    nb_spheres = SettingGet(cSetting_roving_nb_spheres);

    isomesh = SettingGet(cSetting_roving_isomesh);
    isosurface = SettingGet(cSetting_roving_isosurface);

    if(SettingGet(cSetting_roving_byres))
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

      auto_save = (int)SettingGet(cSetting_auto_zoom);
      SettingSet(cSetting_auto_zoom,0);
      
      name = SettingGet_s(NULL,NULL,cSetting_roving_map1_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(name))
            {
              level = SettingGet(cSetting_roving_map1_level);
              sprintf(buffer,
                      "cmd.isomesh('rov_m1','%s',%8.6f,'center',%1.3f)",
                      name,level,isomesh);
              PParse(buffer);
              PFlush();
              refresh_flag=true;
            }

      name = SettingGet_s(NULL,NULL,cSetting_roving_map2_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(name))
            {
              level = SettingGet(cSetting_roving_map2_level);
              sprintf(buffer,
                      "cmd.isomesh('rov_m2','%s',%8.6f,'center',%1.3f)",
                      name,level,isomesh);
              PParse(buffer);
              PFlush();
              refresh_flag=true;
            }

      name = SettingGet_s(NULL,NULL,cSetting_roving_map3_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(name))
            {
              level = SettingGet(cSetting_roving_map3_level);
              sprintf(buffer,
                      "cmd.isomesh('rov_m3','%s',%8.6f,'center',%1.3f)",
                      name,level,isomesh);
              PParse(buffer);
              PFlush();
              refresh_flag=true;
            }


      SettingSet(cSetting_auto_zoom,(float)auto_save);            
    }

    if(isosurface!=0.0F) {
      int auto_save;

      auto_save = (int)SettingGet(cSetting_auto_zoom);
      SettingSet(cSetting_auto_zoom,0.0F);
      
      name = SettingGet_s(NULL,NULL,cSetting_roving_map1_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(name))
            {
              level = SettingGet(cSetting_roving_map1_level);
              sprintf(buffer,
                      "cmd.isosurface('rov_s1','%s',%8.6f,'center',%1.3f)",
                      name,level,isosurface);
              PParse(buffer);
              PFlush();
              refresh_flag=true;
            }

      name = SettingGet_s(NULL,NULL,cSetting_roving_map2_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(name))
            {
              level = SettingGet(cSetting_roving_map2_level);
              sprintf(buffer,
                      "cmd.isosurface('rov_s2','%s',%8.6f,'center',%1.3f)",
                      name,level,isosurface);
              PParse(buffer);
              PFlush();
              refresh_flag=true;
            }

      name = SettingGet_s(NULL,NULL,cSetting_roving_map3_name);
      if(name)
        if(name[0]) 
          if(ExecutiveFindObjectByName(name))
            {
              level = SettingGet(cSetting_roving_map3_level);
              sprintf(buffer,
                      "cmd.isosurface('rov_s3','%s',%8.6f,'center',%1.3f)",
                      name,level,isosurface);
              PParse(buffer);
              PFlush();
              refresh_flag=true;
            }


      SettingSet(cSetting_auto_zoom,(float)auto_save);            
    }


    if(refresh_flag) {
      PParse("cmd.refresh()");
      PFlush();
    }

    I->RovingLastUpdate=UtilGetSeconds();
    I->RovingDirtyFlag=false;
  } 
}
/*========================================================================*/
int SceneDrag(Block *block,int x,int y,int mod)
{
  CScene *I=&Scene;
  float scale,vScale;
  float v1[3],v2[3],n1[3],n2[3],r1,r2,cp[3],v3[3];
  float dx,dy,dt;
  float axis[3],axis2[3],theta,omega,old_z;
  int mode;
  int eff_width;
  int moved_flag;
  int adjust_flag;
  CObject *obj;

  mode = ButModeTranslate(I->Button,mod);
  
  y=y-I->Block->margin.bottom;
  scale = (float)I->Height;
  if(scale > I->Width)
	 scale = (float)I->Width;
  scale = 0.45F * scale;

  SceneDontCopyNext();
  switch(mode) {
  case cButModePickAtom:
    obj=(CObject*)I->LastPicked.ptr;
    if(obj)
      switch(obj->type) {
      case cObjectGadget: {
        ObjectGadget *gad;
        
        gad = (ObjectGadget*)obj;

        ObjectGadgetGetVertex(gad,I->LastPicked.index,I->LastPicked.bond,v1);

        vScale = SceneGetScreenVertexScale(v1);
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
        SceneChanged();
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
  case cButModePkTorBnd:
    obj=(CObject*)I->LastPicked.ptr;
    if(obj)
      switch(obj->type) {
      case cObjectMolecule:
        if(ObjectMoleculeGetAtomVertex((ObjectMolecule*)obj,
                                       SettingGetGlobal_i(cSetting_state)-1,
                                       I->LastPicked.index,v1)) {
          /* scale properly given the current projection matrix */
          vScale = SceneGetScreenVertexScale(v1);
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
          EditorDrag((ObjectMolecule*)obj,I->LastPicked.index,mode,
                     SettingGetGlobal_i(cSetting_state)-1,v1,v2,v3);
        }
        break;
      default:
        break;
      }
    I->LastX=x;
    I->LastY=y;
    break;
  case cButModeTransXY:

    vScale = SceneGetScreenVertexScale(I->Origin);
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
        SceneDirty();
        moved_flag=true;
      }
    if(I->LastY!=y)
      {
        I->Pos[1]+=v2[1];
        I->LastY=y;
        SceneDirty();
        moved_flag=true;
      }
    
    if(moved_flag&&(int)SettingGet(cSetting_roving_origin)) {
      SceneGetPos(v2);
      SceneOriginSet(v2,true);
    }
    if(moved_flag&&(int)SettingGet(cSetting_roving_detail)) {    
      SceneRovingDirty();
    }
    break;
  case cButModeRotXYZ:
  case cButModeRotZ:
  case cButModeTransZ:
  case cButModeClipNF:
  case cButModeClipN:    
  case cButModeClipF:    
    
    eff_width = I->Width;
    if(I->StereoMode>1) {
      eff_width = I->Width/2;
      x = x % eff_width;
    }

    v1[0] = (float)(eff_width/2) - x;
    v1[1] = (float)(I->Height/2) - y;
    
	 v2[0] = (float)(eff_width/2) - I->LastX;
	 v2[1] = (float)(I->Height/2) - I->LastY;
	 
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
	 theta = (float)(SettingGet_f(NULL,NULL,cSetting_mouse_scale)*
      2*180*asin(sqrt1f(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]))/cPI);

    dx = (v1[0]-v2[0]);
    dy = (v1[1]-v2[1]);
    dt = (float)(SettingGet_f(NULL,NULL,cSetting_mouse_limit)*sqrt1f(dx*dx+dy*dy)/scale);
    
    if(theta>dt)
      theta = dt;

	 normalize23f(cp,axis);

    v1[2]=0.0;
    v2[2]=0.0;
	 normalize23f(v1,n1);
	 normalize23f(v2,n2);
	 cross_product3f(n1,n2,cp);
    omega = (float)(2*180*asin(sqrt1f(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]))/cPI);
	 normalize23f(cp,axis2);	 
    old_z = (I->Back + I->FrontSafe )/ 2.0F;
    moved_flag=false;
    adjust_flag=false;
	 switch(mode) {
	 case cButModeRotXYZ:
		if(I->LastX!=x)
		  {
			 SceneRotate(theta,axis[0],axis[1],-axis[2]);
			 I->LastX=x;
          adjust_flag=true;
		  }
		if(I->LastY!=y)
		  {
			 SceneRotate(theta,axis[0],axis[1],-axis[2]);
			 I->LastY=y;
          adjust_flag=true;
		  }
		break;
	 case cButModeRotZ:
		if(I->LastX!=x)
		  {
			 SceneRotate(omega,axis2[0],axis2[1],-axis2[2]);
			 I->LastX=x;
          adjust_flag=true;
		  }
		if(I->LastY!=y)
		  {
			 SceneRotate(omega,axis2[0],axis2[1],-axis2[2]);
			 I->LastY=y;
          adjust_flag=true;		
		  }
      break;
	 case cButModeTransZ:
		if(I->LastY!=y)
		  {
          float factor;
          factor = 200/((I->Front+I->Back)/2);
			 I->Pos[2]+=(((float)y)-I->LastY)/factor;
			 I->Front-=(((float)y)-I->LastY)/factor;
			 I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
			 I->Back-=(((float)y)-I->LastY)/factor;
			 I->LastY=y;
			 SceneDirty();
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
			 SceneDirty();
          moved_flag=true;
		  }
		if(I->LastY!=y)
		  {
			 I->Front-=(((float)y)-I->LastY)/10;
			 if(I->Front>I->Back)
				I->Front=I->Back+cSliceMin;
			 if(I->Front<cFrontMin) I->Front=cFrontMin;
			 I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
			 I->LastY=y;
			 SceneDirty();
          moved_flag=true;
		  }
		break;
	 case cButModeClipN:
		if(I->LastX!=x)
		  {
			 I->Front-=(((float)x)-I->LastX)/10;
			 if(I->Front>I->Back)
				I->Front=I->Back+cSliceMin;
			 if(I->Front<cFrontMin) I->Front=cFrontMin;
			 I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
			 I->LastX=x;
			 SceneDirty();
          moved_flag=true;
		  }
		if(I->LastY!=y)
		  {
			 I->Front-=(((float)y)-I->LastY)/10;
			 if(I->Front>I->Back)
				I->Front=I->Back+cSliceMin;
			 if(I->Front<cFrontMin) I->Front=cFrontMin;
			 I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
			 I->LastY=y;
			 SceneDirty();
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
			 SceneDirty();
          moved_flag=true;
		  }
		if(I->LastY!=y)
		  {
			 I->Back-=(((float)y)-I->LastY)/10;
			 if(I->Back<I->Front)
				I->Back=I->Front+cSliceMin;
			 I->LastY=y;
			 SceneDirty();
          moved_flag=true;
		  }
		break;
    }
    if(moved_flag&&(int)SettingGet(cSetting_roving_origin)) {
      v2[0] = 0.0F;
      v2[1] = 0.0F;
      v2[2] = (I->Back + I->FrontSafe)/2.0F - old_z;

      MatrixInvTransform3f(I->RotMatrix,v2,v2);
      subtract3f(I->Origin,v2,v2);
      SceneOriginSet(v2,true);
    }
    if((adjust_flag)&&(int)SettingGet(cSetting_roving_detail)) {    
      SceneRovingPostpone();
    }
    if((moved_flag)&&(int)SettingGet(cSetting_roving_detail)) {    
      SceneRovingDirty();
    }

  }
  return(1);
}
/*========================================================================*/
void SceneFree(void)
{
  CScene *I=&Scene;
  OrthoFreeBlock(I->Block);
  
  ListFree(I->Obj,next,ObjList);
  if(!I->MovieOwnsImageFlag)
	 FreeP(I->ImageBuffer);
  

  CGOFree(DebugCGO);
}
/*========================================================================*/
void SceneResetMatrix(void)
{
  CScene *I=&Scene;
  MatrixLoadIdentity44f(I->RotMatrix);
}
/*========================================================================*/
void SceneSetDefaultView(void) {
  CScene *I=&Scene;

  MatrixLoadIdentity44f(I->RotMatrix);

  I->Pos[0] = 0.0;
  I->Pos[1] = 0.0;
  I->Pos[2] = -50.0;

  I->Origin[0] = 0.0;
  I->Origin[1] = 0.0;
  I->Origin[2] = 0.0;

  I->Front=40;
  I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
  I->Back=100;

  I->Scale = 1.0;
  
}
int SceneReinitialize(void)
{
  int ok=true;
  SceneSetDefaultView();
  SceneCountFrames();
  SceneSetFrame(0,0);
  SceneDirty();
  return(ok);
}
/*========================================================================*/
void SceneInit(void)
{
  CScene *I=&Scene;

  DebugCGO = CGONew();

  ListInit(I->Obj);

  I->RockTime=0;
  I->TextColor[0]=0.2F;
  I->TextColor[1]=1.0F;
  I->TextColor[2]=0.2F;
  I->SculptingSave=0;
  
  SceneSetDefaultView();

  I->NFrame = 0;
  I->Scale = 1.0;
  I->Block = OrthoNewBlock(NULL);
  I->Block->fClick   = SceneClick;
  I->Block->fRelease = SceneRelease;
  I->Block->fDrag    = SceneDrag;
  I->Block->fDraw    = SceneDraw;
  I->Block->fReshape = SceneReshape;
  I->Block->active = true;

  OrthoAttach(I->Block,cOrthoScene);

  I->DirtyFlag = true;
  I->RovingDirtyFlag = false;
  I->ImageBuffer = NULL;
  I->ImageBufferWidth=0;
  I->ImageBufferHeight=0;
  I->ImageBufferSize = 0;
  I->MovieOwnsImageFlag = false;
  I->MovieFrameFlag = false;
  I->RenderTime = 0;
  I->LastRender = UtilGetSeconds();
  I->LastFrameTime = UtilGetSeconds();
  I->LastRockTime = UtilGetSeconds();
  I->LastPicked.ptr = NULL;

  I->CopyNextFlag=true;
  I->CopyFlag=false;
  I->CopiedFromOpenGL=false;
  I->vendor[0]=0;
  I->renderer[0]=0;
  I->version[0]=0;
}
/*========================================================================*/
void SceneReshape(Block *block,int width,int height)
{
  CScene *I=&Scene;
  
  if(I->Block->margin.right) {
	 width -= I->Block->margin.right;
	 if(width<1)
		width=1;
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
  SceneDirty();

  MovieClearImages();
  MovieSetSize(I->Width,I->Height);

}

float fog_val=1.0;
/*========================================================================*/
void SceneDone(void)
{
  CScene *I=&Scene;
  if(I->Block)
	 OrthoFreeBlock(I->Block);
}
/*========================================================================*/
void SceneResetNormal(int lines)
{
  CScene *I=&Scene;
  if(PMGUI) {
    if(lines)
      glNormal3fv(I->LinesNormal);
    else
      glNormal3fv(I->ViewNormal);
  }
}

/*========================================================================*/

static double accumTiming = 0.0; 

void SceneRay(int ray_width,int ray_height,int mode,char **headerVLA_ptr,
              char **charVLA_ptr,float angle,float shift,int quiet)
{
  CScene *I=&Scene;
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

  fov=SettingGet(cSetting_field_of_view);
  aspRat = ((float) ray_width) / ((float) ray_height);

  ScenePrepareUnitContext(&context,ray_width,ray_height);
  if(SettingGet(cSetting_all_states)) {
    curState=-1;
  } else {
    curState=SettingGetGlobal_i(cSetting_state)-1;
  }

  ray = RayNew();

  SceneUpdate();

  timing = UtilGetSeconds(); /* start timing the process */
  
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


  if(Feedback(FB_Scene,FB_Debugging)) {
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

  OrthoBusyFast(0,20);

  RayPrepare(ray,-width,width,-height,height,
             I->FrontSafe,I->Back,rayView,aspRat,ray_width);

  while(ListIterate(I->Obj,rec,next))
	 {
		if(rec->obj->fRender) {
        RaySetContext(ray,rec->obj->Context);
		  ray->fColor3fv(ray,white);
		  rec->obj->fRender(rec->obj,
                          ObjectGetCurrentState(rec->obj,false),ray,NULL,0);
		}
	 }
  OrthoBusyFast(1,20);

  if(mode!=2) { /* don't show pixel count for tests */
    if(!quiet) {
    PRINTFB(FB_Ray,FB_Blather)
      " Ray: tracing %dx%d = %d rays against %d primitives.\n",ray_width,ray_height,
      ray_width*ray_height,RayGetNPrimitives(ray)
      ENDFB;
    }
  }
  switch(mode) {
  case 0: /* mode 0 is built-in */
    buffer_size = 4*ray_width*ray_height;
    buffer=(GLvoid*)Alloc(char,buffer_size);
    ErrChkPtr(buffer);
    
    RayRender(ray,ray_width,ray_height,buffer,I->Front,I->Back,timing,angle);

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
                 I->FrontSafe,I->Back,fov,angle);
    if(!(charVLA_ptr&&headerVLA_ptr)) { /* immediate mode */
      strcpy(prefix,SettingGet_s(NULL,NULL,cSetting_batch_prefix));
      if(PPovrayRender(headerVLA,charVLA,prefix,ray_width,
                       ray_height,(int)SettingGet(cSetting_antialias))) {
        strcat(prefix,".png");
        SceneLoadPNG(prefix,false,false);
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
    RayRenderTest(ray,ray_width,ray_height,I->FrontSafe,I->Back,fov);
    break;
  }

  timing = UtilGetSeconds()-timing;
  if(mode!=2) { /* don't show timings for tests */
	accumTiming += timing; 

	if(!quiet) {
     PRINTFB(FB_Ray,FB_Details)
       " Ray: total time: %4.2f sec. = %3.1f frames/hour. (%4.2f sec. accum.)\n", 
       timing,3600/timing, 
       accumTiming 
      ENDFB;
   }
  }
  OrthoDirty();
  RayFree(ray);
}
/*========================================================================*/
void SceneCopy(GLenum buffer,int force)
{
  CScene *I=&Scene;
  unsigned int buffer_size;

  if(force || (!(I->StereoMode||SettingGet(cSetting_stereo_double_pump_mono))))
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
          ErrChkPtr(I->ImageBuffer);
          I->ImageBufferSize = buffer_size;
          I->ImageBufferWidth=I->Width;
          I->ImageBufferHeight=I->Height;
        }
        if(PMGUI) {
          glReadBuffer(buffer);
          glReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
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
int SceneRovingCheckDirty(void)
{
  CScene *I=&Scene;

  return(I->RovingDirtyFlag);
}
/*========================================================================*/
void SceneUpdate(void)
{
  CScene *I=&Scene;
  ObjRec *rec=NULL;

  PRINTFD(FB_Scene)
    " SceneUpdate: entered.\n"
    ENDFD;
  if(I->ChangedFlag) {
    SceneCountFrames();
	 while(ListIterate(I->Obj,rec,next))
      if(rec->obj->fUpdate) 
        rec->obj->fUpdate(rec->obj);
	 I->ChangedFlag=false;
  }
  PRINTFD(FB_Scene)
    " SceneUpdate: leaving...\n"
    ENDFD;
}
/*========================================================================*/
int SceneRenderCached(void)
{
  /* sets up a cached image buffer is one is available, or if we are
   * using cached images by default */
  CScene *I=&Scene;
  ImageType image;
  int renderedFlag=false;

  PRINTFD(FB_Scene)
    " SceneRenderCached: entered.\n"
    ENDFD;

  if(I->DirtyFlag) {
	if(I->MovieFrameFlag||
	   (MoviePlaying()&&SettingGet(cSetting_cache_frames))) {
	  I->MovieFrameFlag=false;
	  image = MovieGetImage(MovieFrameToImage(SettingGetGlobal_i(cSetting_frame)-1));
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
		  OrthoDirty();
		  renderedFlag=true;
		}
	  else
		{
		  SceneMakeMovieImage();
		  renderedFlag=true;
		}
	} else if(MoviePlaying()&&SettingGet(cSetting_ray_trace_frames)) {
	  SceneRay(0,0,(int)SettingGet(cSetting_ray_default_renderer),NULL,NULL,0.0F,0.0F,false); 
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
	I->LastRender = UtilGetSeconds();
	I->RenderTime += I->LastRender;
	ButModeSetRate(I->RenderTime);
   }*/

  PRINTFD(FB_Scene)
    " SceneRenderCached: leaving...renderedFlag %d\n",renderedFlag
    ENDFD;

  return(renderedFlag);
}


/*========================================================================*/
static void SceneRenderAll(SceneUnitContext *context,float *normal,Pickable **pickVLA,int pass,int fat)
{
  CScene *I=&Scene;
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

          glPopMatrix();
          glMatrixMode(GL_PROJECTION);
          glPopMatrix();
          glMatrixMode(GL_MODELVIEW);
          break;
        }
      glPopMatrix();
    }
}

/*========================================================================*/
void SceneRender(Pickable *pick,int x,int y,Multipick *smp)
{
  /* think in terms of the camera's world */
  CScene *I=&Scene;
  ObjRec *rec=NULL;
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
  GLenum render_buffer = GL_BACK;

  SceneUnitContext context;

  PRINTFD(FB_Scene)
    " SceneRender: entered. pick %p x %d y %d smp %p\n",
    pick,x,y,smp
    ENDFD;

  double_pump=SettingGet_i(NULL,NULL,cSetting_stereo_double_pump_mono);
  
  if(I->StereoMode>1)
    aspRat=aspRat/2;

  fov=SettingGet(cSetting_field_of_view);
  if(PMGUI) {

    must_render_stereo = (I->StereoMode!=0);
    if(!must_render_stereo) 
      if(double_pump&&StereoCapable) {            /* force stereo rendering */
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
  
    glGetIntegerv(GL_VIEWPORT,(GLint*)view_save);
    glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height);
    
    debug_pick = (int)SettingGet(cSetting_debug_pick);

    if(SettingGet(cSetting_line_smooth)) {
      if(!(pick||smp)) {
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
      }
      glLineWidth(0.0);
    } else {
      glLineWidth(SettingGet(cSetting_line_width));
      glDisable(GL_LINE_SMOOTH);
    }

    glPointSize(SettingGet(cSetting_dot_width));

    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);

    /* get matrixes for unit objects */

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    ScenePrepareUnitContext(&context,I->Width,I->Height);

    /* do standard 3D objects */

    /* Set up the clipping planes */
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if(SettingGet(cSetting_all_states)) {
      curState=-1;
    } else {
      curState=SettingGetGlobal_i(cSetting_state)-1;
    }

    if(SettingGet(cSetting_ortho)==0.0) {
      gluPerspective(fov,aspRat,I->FrontSafe,I->Back);
    } else {
      height  = (float)(fabs(I->Pos[2])*tan((fov/2.0)*cPI/180.0));	 
      width = height*aspRat;
	
      glOrtho(-width,width,-height,height,
              I->FrontSafe,I->Back);

    }

    glMatrixMode(GL_MODELVIEW);
    ScenePrepareMatrix(0);

    /* Save these for editing operations */

    glGetFloatv(GL_MODELVIEW_MATRIX,I->ModMatrix);
    glGetFloatv(GL_PROJECTION_MATRIX,I->ProMatrix);
  
    /* determine the direction in which we are looking relative*/

    /* 2. set the normals to reflect light back at the camera */

    MatrixInvTransform3f(I->RotMatrix,zAxis,normal); 
    copy3f(normal,I->ViewNormal);
  
    if(SettingGet(cSetting_normal_workaround)) {
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
      glShadeModel(GL_SMOOTH);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_DITHER);

      f=SettingGet(cSetting_gl_ambient);
      vv[0]=f;
      vv[1]=f;
      vv[2]=f;
      vv[3]=1.0;

      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT,vv);
      if(SettingGet(cSetting_two_sided_lighting)) {
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
      } else {
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
      }

      f = SettingGet(cSetting_specular);
      if(f>R_SMALL4) {
        /*        glEnable(GL_LIGHT1);*/
        /*        glLightfv(GL_LIGHT1,GL_AMBIENT,zero);*/
        vv[0]=f;
        vv[1]=f;
        vv[2]=f;
        vv[3]=1.0;

        glEnable(GL_LIGHT1);
        glLightfv(GL_LIGHT0,GL_SPECULAR,zero);
        glLightfv(GL_LIGHT1,GL_SPECULAR,vv);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,zero);
        glMaterialfv(GL_FRONT,GL_SPECULAR,vv);
        vv[0]=SettingGet(cSetting_shininess);
        glMaterialfv(GL_FRONT,GL_SHININESS,vv);

        copy3f(SettingGetGlobal_3fv(cSetting_light),vv);
        normalize3f(vv);
        MatrixInvTransform44fAs33f3f(I->RotMatrix,vv,vv); 
        invert3f(vv);
        vv[3]=0.0;
        glLightfv(GL_LIGHT1,GL_POSITION,vv);

      } else {
        glMaterialfv(GL_FRONT,GL_SPECULAR,zero); 
      }
      if(SettingGet(cSetting_depth_cue)&&SettingGet(cSetting_fog)) {
        fog_start = (I->Back-I->FrontSafe)*SettingGet(cSetting_fog_start)+I->FrontSafe;
#ifdef _PYMOL_3DFX
        if(SettingGet(cSetting_ortho)==0.0) {
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
          glFogf(GL_FOG_END, I->FrontSafe+(I->Back-I->FrontSafe)/SettingGet(cSetting_fog));
          glFogf(GL_FOG_DENSITY, fog_val);
#endif
          v=SettingGetfv(cSetting_bg_rgb);
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
      glShadeModel(GL_FLAT);

    }

    PRINTFD(FB_Scene)
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

      SceneRenderAll(&context,NULL,&pickVLA,0,true);
	  

      if(debug_pick) {
        p_glutSwapBuffers();
        PSleep(1000000*debug_pick/4);
        p_glutSwapBuffers();
      }
      lowBits = SceneFindTriplet(x,y,render_buffer);
	
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA[0].index=0;
      pickVLA[0].ptr=(void*)pick; /* this is just a flag */
	
      SceneRenderAll(&context,NULL,&pickVLA,0,true);

      if(debug_pick) {
        p_glutSwapBuffers();
        PSleep(1000000*debug_pick/4);
        p_glutSwapBuffers();
      }

      highBits = SceneFindTriplet(x,y,render_buffer);
      index = lowBits+(highBits<<12);

      if(index&&(index<=pickVLA[0].index)) {
        *pick = pickVLA[index]; /* return object info */
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
      
      SceneRenderAll(&context,NULL,&pickVLA,0,true);
      
      lowBitVLA = SceneReadTriplets(smp->x,smp->y,smp->w,smp->h,render_buffer);

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA[0].index=0;
      pickVLA[0].ptr=(void*)smp; /* this is just a flag */
	
      SceneRenderAll(&context,NULL,&pickVLA,0,true);

      highBitVLA = SceneReadTriplets(smp->x,smp->y,smp->w,smp->h,render_buffer);
      
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

      ButModeCaptionReset(); /* reset the frame caption if any */
      /* rendering for visualization */


      PRINTFD(FB_Scene)
        " SceneRender: I->StereoMode %d must_render_stereo %d\n    stereo_as_mono %d  StereoCapable %d\n",
        I->StereoMode, must_render_stereo, stereo_as_mono, StereoCapable
        ENDFD;

      start_time = UtilGetSeconds();
      if(must_render_stereo) {
        /*stereo*/

        PRINTFD(FB_Scene)
          " SceneRender: left hand stereo...\n"
          ENDFD;
        
        /* LEFT HAND STEREO */

        if(stereo_as_mono) {
          glDrawBuffer(GL_BACK_LEFT);
        } else switch(I->StereoMode) {
        case 1: /* hardware */
          glDrawBuffer(GL_BACK_LEFT);
          break;
        case 2: /* side by side */
          glViewport(I->Block->rect.left+I->Width/2,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        case 3:
          glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        }
        
        glPushMatrix(); /* 1 */
        ScenePrepareMatrix(stereo_as_mono ? 0 : 1);
        for(pass=1;pass>=0;pass--) { /* render opaque then antialiased...*/
          rec=NULL;

          SceneRenderAll(&context,normal,NULL,pass,false);

        }

        glPushMatrix();  /* 2 */
        glNormal3fv(normal);
        CGORenderGL(DebugCGO,NULL,NULL,NULL);
        glPopMatrix();  /* 1 */


        glPushMatrix(); /* 2 */
        glNormal3fv(normal);
        ExecutiveRenderSelections(curState);
        EditorRender(curState);
        glPopMatrix(); /* 1 */

        /* render transparent */

        SceneRenderAll(&context,normal,NULL,-1,false);

        glPopMatrix(); /* 0 */

        /* RIGHT HAND STEREO */

        PRINTFD(FB_Scene)
          " SceneRender: right hand stereo...\n"
          ENDFD;

        if(stereo_as_mono) { /* double pumped mono */
          glDrawBuffer(GL_BACK_RIGHT);
        } else switch(I->StereoMode) {
        case 1: /* hardware */
          glDrawBuffer(GL_BACK_RIGHT);
          break;
        case 2: /* side by side */
          glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        case 3:
          glViewport(I->Block->rect.left+I->Width/2,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        }
        
        glClear(GL_DEPTH_BUFFER_BIT);        
        glPushMatrix(); /* 1 */
        ScenePrepareMatrix(stereo_as_mono ? 0 : 2);
        for(pass=1;pass>=0;pass--) { /* render opaque then antialiased...*/
          
          SceneRenderAll(&context,normal,NULL,pass,false);
          
        }
        
        glPushMatrix(); /* 2 */
        glNormal3fv(normal);
        CGORenderGL(DebugCGO,NULL,NULL,NULL);
        glPopMatrix(); /* 1 */
        
        glPushMatrix(); /* 2 */
        glNormal3fv(normal);
        ExecutiveRenderSelections(curState);
        EditorRender(curState);
        glPopMatrix(); /* 1 */
        
        /* render transparent */
        SceneRenderAll(&context,normal,NULL,-1,false);
        
        glPopMatrix(); /* 0 */

        /* restore draw buffer */

        if(must_render_stereo) { /* double pumped mono */
          glDrawBuffer(GL_BACK);
        } else switch(I->StereoMode) {
        case 1: /* hardware */
          glDrawBuffer(GL_BACK);
          break;
        case 2: /* side by side */
        case 3:
          glDrawBuffer(GL_BACK);
          break;
        }

      } else {

        /* mono rendering */

        PRINTFD(FB_Scene)
          " SceneRender: rendering opaque and antialiased...\n"
          ENDFD;
        
        for(pass=1;pass>=0;pass--) { /* render opaque then antialiased...*/
          SceneRenderAll(&context,normal,NULL,pass,false);
        }

        PRINTFD(FB_Scene)
          " SceneRender: rendering DebugCGO...\n"
          ENDFD;
        
        glPushMatrix();
        glNormal3fv(normal);
        CGORenderGL(DebugCGO,NULL,NULL,NULL);
        glPopMatrix();

        PRINTFD(FB_Scene)
          " SceneRender: rendering selections...\n"
          ENDFD;

        glPushMatrix();
        glNormal3fv(normal);
        ExecutiveRenderSelections(curState);

        PRINTFD(FB_Scene)
          " SceneRender: rendering editing...\n"
          ENDFD;

        EditorRender(curState);

        PRINTFD(FB_Scene)
          " SceneRender: rendering transparent objects...\n"
          ENDFD;

        /* render transparent */
        SceneRenderAll(&context,normal,NULL,-1,false);
        glPopMatrix();

      }
    }
    
    if(!(pick||smp)) {
      glDisable(GL_FOG);
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT1);
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_DITHER);
    }
    glLineWidth(1.0);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glViewport(view_save[0],view_save[1],view_save[2],view_save[3]);
  }

  PRINTFD(FB_Scene)
    " SceneRender: rendering complete.\n"
    ENDFD;
  
  if(!(pick||smp)) { /* update frames per second field */
    I->RenderTime = -I->LastRender;
    I->LastRender = UtilGetSeconds();
    I->RenderTime += I->LastRender;
    ButModeSetRate((float)I->RenderTime);
    if(I->CopyNextFlag) {
      start_time = I->LastRender - start_time;
      if((start_time>0.10)||(MainSavingUnderWhileIdle()))
        if(!(ControlIdling()))
          if(SettingGet(cSetting_cache_display)) {
            SceneCopy(render_buffer,false);
          }
    } else {
      I->CopyNextFlag=true;
    }
  }
  PRINTFD(FB_Scene)
    " SceneRender: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void SceneRestartTimers(void)
{
  CScene *I=&Scene;
  I->LastRender = UtilGetSeconds();
  I->LastFrameTime = UtilGetSeconds();
  I->RenderTime = 0;
}
/*========================================================================*/
void ScenePrepareMatrix(int mode)
{
  CScene *I=&Scene;
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

    stAng = SettingGet(cSetting_stereo_angle);
    stShift = SettingGet(cSetting_stereo_shift);


    stShift = (float)(stShift*fabs(I->Pos[2])/100.0);

    stAng = (float)(stAng*atan(stShift/fabs(I->Pos[2]))*90.0/cPI);


    if(mode==2) { /* left hand */
      stAng=-stAng;
      stShift=-stShift;
    }

      PRINTFD(FB_Scene)
        " StereoMatrix-Debug: mode %d stAng %8.3f stShift %8.3f \n",mode,stAng,stShift
        ENDFD;


    fflush(stdout);

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
void SceneRotate(float angle,float x,float y,float z)
{
  CScene *I=&Scene;
  float temp[16];
  int a;
  angle = (float)(-PI*angle/180.0);
  MatrixLoadIdentity44f(temp);
  MatrixRotate44f3f(temp,angle,x,y,z);
  MatrixMultiply44f(I->RotMatrix,temp);
  for(a=0;a<16;a++)
    I->RotMatrix[a]=temp[a];
  SceneDirty();

    /*  glPushMatrix();
        glLoadIdentity();
        glRotatef(angle,x,y,z);
        glMultMatrixf(I->RotMatrix);
        glGetFloatv(GL_MODELVIEW_MATRIX,I->RotMatrix);
        glPopMatrix();*/
}
/*========================================================================*/
void SceneApplyMatrix(float *m)
{
  CScene *I=&Scene;
  MatrixMultiply44f(m,I->RotMatrix);
  SceneDirty();
  
  /*  glPushMatrix();
      glLoadIdentity();
      glMultMatrixf(m);
      glMultMatrixf(I->RotMatrix);
      glGetFloatv(GL_MODELVIEW_MATRIX,I->RotMatrix);
  glPopMatrix();*/
}
/*========================================================================*/
void SceneScale(float scale)
{
  CScene *I=&Scene;
  I->Scale*=scale;
  SceneDirty();
}










