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

#define cFrontMin 0.1
#define cSliceMin 0.1

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
  float Scale;
  int Width,Height;
  int Button;
  int LastX,LastY;
  float ViewNormal[3],LinesNormal[3];
  float Pos[3],Origin[3];
  float H;
  float Front,Back,FrontSafe;
  float TextColor[3];
  float RockTime;
  int DirtyFlag;
  int ChangedFlag;
  int CopyFlag,CopyNextFlag;
  int NFrame;
  GLvoid *ImageBuffer;
  int ImageBufferHeight,ImageBufferWidth;
  int MovieOwnsImageFlag;
  int MovieFrameFlag;
  unsigned ImageBufferSize;
  double LastRender,RenderTime,LastFrameTime;
  float LastRock,LastRockTime;
  Pickable LastPicked;
  int StereoMode;
  OrthoLineType vendor,renderer,version;
  int SculptingFlag,SculptingSave;
} CScene;

CScene Scene;


unsigned int SceneFindTriplet(int x,int y);
unsigned int *SceneReadTriplets(int x,int y,int w,int h);

void SceneDraw(Block *block);
int SceneClick(Block *block,int button,int x,int y,int mod);
int SceneRelease(Block *block,int button,int x,int y,int mod);

int SceneDrag(Block *block,int x,int y,int mod);
void ScenePrepareMatrix(int mode);

static int SceneGetObjState(CObject *obj,int state)
{
  int objState;
  if(SettingGetIfDefined_i(obj->Setting,cSetting_state,&objState)) {
    if(objState>0) {
      state=objState-1;
    } if(objState<0) {
      state=-1;
    }
  }
  return(state);
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
void SceneSetView(SceneViewType view)
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
  PRINTFB(FB_Scene,FB_Actions)
    " Scene: view updated.\n"
    ENDFB;

}
/*========================================================================*/
void SceneDontCopyNext(void)
/* disables automatic copying of the image for the next rendering run */
{
  CScene *I=&Scene;
  I->CopyNextFlag=false;
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
void SceneClip(int plane,float movement) /* 0=front, 1=back*/
{
  CScene *I=&Scene;
  float avg;

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
    avg = (I->Front+I->Back)/2.0;
    SceneClipSet(avg-movement,avg+movement);
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

  buffer_size = 4*I->Width*I->Height;
  if(!I->CopyFlag) {
	 image = (GLvoid*)Alloc(char,buffer_size);
	 ErrChkPtr(image);
    if(PMGUI) {
      glReadBuffer(GL_BACK);
      glReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                   GL_RGBA,GL_UNSIGNED_BYTE,image);
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
  persp=!flag;
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
  int newFrame;
  int newState;

  newFrame = SettingGetGlobal_i(cSetting_frame) -1;
  PRINTFD(FB_Scene)
    " SceneSetFrame: entered.\n"
    ENDFD;
  switch(mode) {
  case 0:
    newFrame=frame; 
	 break;
  case 1:
    newFrame+=frame; 
	 break;
  case 2:
    newFrame=I->NFrame-1; 
	 break;
  case 3:
	 newFrame=I->NFrame/2;
	 break;
  case 4:
	 newFrame=frame;
	 break;
  case 5:
	 newFrame+=frame;
	 break;
  case 6: /* movie/frame override - go to this state absolutely! */
    newState = frame;
    break;
  }
  SceneCountFrames();
  if (mode<6) { 
    if(newFrame>=I->NFrame) newFrame=I->NFrame-1;
    if(newFrame<0) newFrame=0;
    newState = MovieFrameToIndex(newFrame);
    if(mode&4) 
      MovieDoFrameCommand(newFrame);
    if(newFrame==0)
      MovieMatrix(cMovieMatrixRecall);
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
	SceneRay(0,0,SettingGet(cSetting_ray_default_renderer),NULL,NULL); 
  } else {
	 v=SettingGetfv(cSetting_bg_rgb);
    if(PMGUI) {
      glDrawBuffer(GL_BACK);
      glClearColor(v[0],v[1],v[2],1.0);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glClearColor(0.0,0.0,0.0,1.0);
      SceneRender(NULL,0,0,NULL);
      SceneCopy(0);
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
  float renderTime;
  float minTime;
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
    ang_cur = (I->RockTime*SettingGet(cSetting_sweep_speed));
    
    disp = SettingGet(cSetting_sweep_angle)*(3.1415/180.0)*sin(ang_cur)/2;
    diff = disp-I->LastRock;
    SceneRotate(180*diff/3.1415,0.0,1.0,0.0);
    I->LastRock = disp;
  }
  if(MoviePlaying()&&frameFlag)
	 {
      I->LastFrameTime = UtilGetSeconds();
      if((SettingGetGlobal_i(cSetting_frame)-1)==(I->NFrame-1))
        SceneSetFrame(4,0);
      else
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

  dist = radius/tan((fov/2.0)*cPI/180.0);

  I->Pos[2]-=dist;
  I->Front=(-I->Pos[2]-radius*1.1);
  I->FrontSafe=(I->Front<cFrontMin ? cFrontMin : I->Front);  
  I->Back=(-I->Pos[2]+radius*1.4);
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
  MatrixTransform3f(I->RotMatrix,I->Origin,v0); 
  dist = I->Pos[2];

  /* find where this point is in relationship to the origin */
  subtract3f(I->Origin,location,v0); 

  /*  printf("%8.3f %8.3f %8.3f\n",I->Front,I->Pos[2],I->Back);*/

  MatrixTransform3f(I->RotMatrix,v0,I->Pos); /* convert to view-space */

  I->Pos[2]=dist;
  I->Front=(-I->Pos[2]-(slab_width*0.45));
  I->FrontSafe=(I->Front<cFrontMin ? cFrontMin : I->Front);  
  I->Back=(-I->Pos[2]+(slab_width*0.55));

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

  while(ListIterate(I->Obj,rec,next))
	 if(rec->obj==obj)
		break;
  if(rec) {
    rec->obj->Enabled=false;
	 ListDetach(I->Obj,rec,next,ObjList);
	 ListElemFree(rec);
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

  if(PMGUI) {
    overlay = SettingGet(cSetting_overlay);
    text = SettingGet(cSetting_text);

    if(overlay||(!text)) 

      if(I->CopyFlag)
        {
          glReadBuffer(GL_BACK); /* neither of these seems to fix the nVidia driver under windows */
          glDrawBuffer(GL_BACK); /* in full-screen mode...pixels still get written to the front buffer! */

          if(I->ImageBufferHeight>I->Height||I->ImageBufferWidth>I->Width) {
            glColor3f(1.0,0.2,0.2);
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
            glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);            
          }
          I->RenderTime = -I->LastRender;
          I->LastRender = UtilGetSeconds();
          I->RenderTime += I->LastRender;
          ButModeSetRate(I->RenderTime);
        }
    
    glColor3f(1.0,1.0,1.0);
  }
}
/*========================================================================*/

typedef unsigned char pix[4];
#define cRange 5

unsigned int SceneFindTriplet(int x,int y) 
{
  int result = 0;
  pix buffer[cRange*2+1][cRange*2+1];
  int a,b,d,flag;
  int debug = false;
  unsigned char *c;
  int strict = false;
  int rb,gb,bb;

  if(PMGUI) { /*just in case*/
  
	glGetIntegerv(GL_RED_BITS,&rb);
	glGetIntegerv(GL_RED_BITS,&gb);
	glGetIntegerv(GL_RED_BITS,&bb);

	if((rb>=8)&&(gb>=8)&&(bb>=8))
		strict = true;

    if(Feedback(FB_Scene,FB_Debugging)) debug=true;
    
    glReadBuffer(GL_BACK);
    glReadPixels(x-cRange,y-cRange,cRange*2+1,cRange*2+1,GL_RGBA,GL_UNSIGNED_BYTE,&buffer[0][0][0]);
    
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
    
    flag=true;
    for(d=0;flag&&(d<cRange);d++)
      for(a=-d;flag&&(a<=d);a++)
		  for(b=-d;flag&&(b<=d);b++) {
			c = &buffer[a+cRange][b+cRange][0];
			if((c[3]==0xFF)&&
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
unsigned int *SceneReadTriplets(int x,int y,int w,int h)
{ 
  unsigned int *result = NULL;
  pix *buffer=NULL;
  int a,b;
  unsigned char *c;
  int cc = 0;
  int dim[3];
  int strict = false;
  int rb,gb,bb;

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
    glReadBuffer(GL_BACK);
    glReadPixels(x,y,w,h,GL_RGBA,GL_UNSIGNED_BYTE,&buffer[0][0]);
    
    for(a=0;a<w;a++)
      for(b=0;b<h;b++)
        {
          c = &buffer[a+b*w][0];
          if((c[3]==0xFF)&&
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
    obj->AtomInfo[I->LastPicked.index].protected=I->SculptingSave;
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
	 if(I->LastPicked.ptr) {
		obj=(CObject*)I->LastPicked.ptr;
      if(obj->type==cObjectMolecule) {
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
      } else {
      EditorSetActiveObject(NULL,0);
      }
    } else {
      EditorSetActiveObject(NULL,0);
    }
    SceneDirty();
    break;
  case cButModePickBond:

    if(I->StereoMode>1)
      x = x % (I->Width/2);

    if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
      SceneRender(NULL,0,0,NULL); /* remove overlay if present */
    SceneDontCopyNext();
    I->LastPicked.ptr = NULL;
	 SceneRender(&I->LastPicked,x,y,NULL);
	 if(I->LastPicked.ptr) {
		obj=(CObject*)I->LastPicked.ptr;
      if(obj->type==cObjectMolecule) {
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

        }
      } else {
        EditorSetActiveObject(NULL,0);
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
      if(Feedback(FB_ObjectMolecule,FB_Results)) {
        if(obj->fDescribeElement) 
          obj->fDescribeElement(obj,I->LastPicked.index,buffer);
        PRINTF " You clicked %s",buffer ENDF;        
        OrthoRestorePrompt();
      }
      objMol = (ObjectMolecule*)obj;
      EditorPrepareDrag(objMol,I->LastPicked.index,
                        SettingGetGlobal_i(cSetting_state)-1);
      y=y-I->Block->margin.bottom;
      x=x-I->Block->margin.left;


      I->LastX=x;
      I->LastY=y;	
      I->SculptingFlag = 1;
      I->SculptingSave =  objMol->AtomInfo[I->LastPicked.index].protected;
      objMol->AtomInfo[I->LastPicked.index].protected=2;
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
    if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
      SceneRender(NULL,0,0,NULL); /* remove overlay if present */
    SceneDontCopyNext();

    if(I->StereoMode>1)
      x = x % (I->Width/2);

	 SceneRender(&I->LastPicked,x,y,NULL);
	 if(I->LastPicked.ptr) {
		obj=(CObject*)I->LastPicked.ptr;
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
  p1[0]=(p1[0]+1.0)*(I->Width/2.0); /* viewport transformation */
  p1[1]=(p1[1]+1.0)*(I->Height/2.0);
  p2[0]=(p2[0]+1.0)*(I->Width/2.0);
  p2[1]=(p2[1]+1.0)*(I->Height/2.0);
  vl=diff3f(p1,p2);
  if(vl<R_SMALL4)
    vl=100.0;
  
  return(1.0/vl);
}

/*========================================================================*/
int SceneDrag(Block *block,int x,int y,int mod)
{
  CScene *I=&Scene;
  float scale,vScale;
  float v1[3],v2[3],n1[3],n2[3],r1,r2,cp[3];
  float axis[3],axis2[3],theta,omega;
  int mode;
  int eff_width;
  CObject *obj;

  mode = ButModeTranslate(I->Button,mod);
  
  y=y-I->Block->margin.bottom;
  scale = I->Height;
  if(scale > I->Width)
	 scale = I->Width;
  scale = 0.45 * scale;

  SceneDontCopyNext();
  switch(mode) {
  case cButModeMovFrag:
  case cButModeTorFrag:
  case cButModeRotFrag:

    obj=(CObject*)I->LastPicked.ptr;
    if(obj)
      if(obj->type==cObjectMolecule) {
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
          /* transform into model coodinate space */
          MatrixInvTransform44fAs33f3f(I->RotMatrix,v2,v2); 
          EditorDrag((ObjectMolecule*)obj,I->LastPicked.index,mode,
                     SettingGetGlobal_i(cSetting_state)-1,v1,v2);
        }
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
    v2[2] = 0;
    
    if(I->LastX!=x)
      {
        I->Pos[0]+=v2[0];
        I->LastX=x;
        SceneDirty();
      }
    if(I->LastY!=y)
      {
        I->Pos[1]+=v2[1];
        I->LastY=y;
        SceneDirty();
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

    v1[0] = (eff_width/2) - x;
    v1[1] = (I->Height/2) - y;
    
	 v2[0] = (eff_width/2) - I->LastX;
	 v2[1] = (I->Height/2) - I->LastY;
	 
	 r1 = sqrt1f(v1[0]*v1[0] + v1[1]*v1[1]);
	 r2 = sqrt1f(v2[0]*v2[0] + v2[1]*v2[1]);
	 
	 if(r1<scale) {
		v1[2] = sqrt1f(scale*scale - r1*r1);
	 }
	 else {
		v1[2] = 0.0;
	 }

	 if(r2<scale) {
		v2[2] = sqrt1f(scale*scale - r2*r2);
	 } else {
		v2[2] = 0.0;
	 }
	 normalize23f(v1,n1);
	 normalize23f(v2,n2);
	 cross_product3f(n1,n2,cp);
	 theta = 2*180*asin(sqrt1f(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]))/3.14;
	 normalize23f(cp,axis);

    v1[2]=0.0;
    v2[2]=0.0;
	 normalize23f(v1,n1);
	 normalize23f(v2,n2);
	 cross_product3f(n1,n2,cp);
    omega = 2*180*asin(sqrt1f(cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]))/3.14;
	 normalize23f(cp,axis2);	 

	 switch(mode) {
	 case cButModeRotXYZ:
		if(I->LastX!=x)
		  {
			 SceneRotate(theta,axis[0],axis[1],-axis[2]);
			 I->LastX=x;
		  }
		if(I->LastY!=y)
		  {
			 SceneRotate(theta,axis[0],axis[1],-axis[2]);
			 I->LastY=y;
		  }
		break;
	 case cButModeRotZ:
		if(I->LastX!=x)
		  {
			 SceneRotate(omega,axis2[0],axis2[1],-axis2[2]);
			 I->LastX=x;
		  }
		if(I->LastY!=y)
		  {
			 SceneRotate(omega,axis2[0],axis2[1],-axis2[2]);
			 I->LastY=y;
		  }
		break;
	 case cButModeTransZ:
		if(I->LastY!=y)
		  {
			 I->Pos[2]+=(((float)y)-I->LastY)/5;
			 I->Front-=(((float)y)-I->LastY)/5;
			 I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
			 I->Back-=(((float)y)-I->LastY)/5;
			 I->LastY=y;
			 SceneDirty();
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
		  }
		if(I->LastY!=y)
		  {
			 I->Back-=(((float)y)-I->LastY)/10;
			 if(I->Back<I->Front)
				I->Back=I->Front+cSliceMin;
			 I->LastY=y;
			 SceneDirty();
		  }
		break;
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
  
  if(I->StereoMode==1) {
    PSGIStereo(0);
  }

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
/*========================================================================*/
void SceneInit(void)
{
  CScene *I=&Scene;

  DebugCGO = CGONew();

  ListInit(I->Obj);

  I->RockTime=0;
  I->TextColor[0]=0.2;
  I->TextColor[1]=1.0;
  I->TextColor[2]=0.2;
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
void SceneRay(int ray_width,int ray_height,int mode,char **headerVLA_ptr,char **charVLA_ptr)
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

  if((!ray_width)||(!ray_height)) {
    ray_width=I->Width;
    ray_height=I->Height;
  }

  fov=SettingGet(cSetting_field_of_view);
  aspRat = ((float) ray_width) / ((float) ray_height);


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

  /* move the camera so that we can see the origin 
	* NOTE, vector is given in the coordinates of the world's motion
	* relative to the camera */
  
  /* turn on depth cuing and all that jazz */
  
  /* 4. rotate about the origin (the the center of rotation) */
  MatrixMultiply44f(I->RotMatrix,rayView);
  
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

  height  = abs(I->Pos[2])*tan((fov/2.0)*cPI/180.0);	 
  width = height*aspRat;

  RayPrepare(ray,-width,width,-height,height,I->FrontSafe,I->Back,rayView);

  while(ListIterate(I->Obj,rec,next))
	 {
		if(rec->obj->fRender) {
		  ray->fColor3fv(ray,white);
		  rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),ray,NULL,0);
		}
	 }

  if(mode!=2) { /* don't show pixel count for tests */
    PRINTFB(FB_Ray,FB_Details)
      " Ray: tracing %dx%d = %d rays...\n",ray_width,ray_height,
      ray_width*ray_height
      ENDFB;
  }
  switch(mode) {
  case 0: /* mode 0 is built-in */
    buffer_size = 4*ray_width*ray_height;
    buffer=(GLvoid*)Alloc(char,buffer_size);
    ErrChkPtr(buffer);
    
    RayRender(ray,ray_width,ray_height,buffer,I->Front,I->Back,timing);
    
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
    I->MovieOwnsImageFlag = false;
    break;

  case 1: /* mode 1 is povray */
    charVLA=VLAlloc(char,100000); 
    headerVLA=VLAlloc(char,2000);
    RayRenderPOV(ray,ray_width,ray_height,&headerVLA,&charVLA,I->FrontSafe,I->Back,fov);
    if(!(charVLA_ptr&&headerVLA_ptr)) { /* immediate mode */
      strcpy(prefix,SettingGet_s(NULL,NULL,cSetting_batch_prefix));
      if(PPovrayRender(headerVLA,charVLA,prefix,ray_width,ray_height,SettingGet(cSetting_antialias))) {
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
    PRINTFB(FB_Ray,FB_Details)
      " Ray: total rendering time: %4.2f sec. = %3.1f frames per hour.\n", 
      timing,3600/timing 
      ENDFB;
  }
  OrthoDirty();
  RayFree(ray);
}
/*========================================================================*/
void SceneCopy(int buffer)
{
  CScene *I=&Scene;
  unsigned int buffer_size;

    if(!I->StereoMode) { /* no copies while in stereo mode */
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
            if(buffer)
              glReadBuffer(GL_FRONT);
            else
              glReadBuffer(GL_BACK);
            glReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                         GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);
            I->ImageBufferWidth=I->Width;
            I->ImageBufferHeight=I->Height;
          }
        }
        I->CopyFlag = true;
      }
    }

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
	  SceneRay(0,0,SettingGet(cSetting_ray_default_renderer),NULL,NULL); 
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
  float start_time=0;
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

  PRINTFD(FB_Scene)
    " SceneRender: entered. pick %p x %d y %d smp %p\n",
    pick,x,y,smp
    ENDFD;

  if(I->StereoMode>1)
    aspRat=aspRat/2;

  fov=SettingGet(cSetting_field_of_view);
  if(PMGUI) {
    glDrawBuffer(GL_BACK);
  
    glGetIntegerv(GL_VIEWPORT,(GLint*)view_save);
    glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height);
    
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
      height  = abs(I->Pos[2])*tan((fov/2.0)*cPI/180.0);	 
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
    if(SettingGet(cSetting_line_smooth)) {
      if(!(pick||smp)) {
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
      }
      glLineWidth(0.0);
    } else {
      glLineWidth(SettingGet(cSetting_line_width));
      glDisable(GL_LINE_SMOOTH);
    }

    glPointSize(SettingGet(cSetting_dot_width));

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    
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
        glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
        break;
      }
	
      glClearColor(0.0,0.0,0.0,0.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


      pickVLA=VLAlloc(Pickable,5000);
      pickVLA[0].index=0;
      pickVLA[0].ptr=NULL;
      while(ListIterate(I->Obj,rec,next))
        {
		  glLineWidth(3.0); /* insure some 100% pixels */
          glPushMatrix();
          if(rec->obj->fRender)
            rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,&pickVLA,0);
          glPopMatrix();
        }

	  
	 /* p_glutSwapBuffers();
	  PSleep(5000000);*/

      lowBits = SceneFindTriplet(x,y);
	
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA[0].index=0;
      pickVLA[0].ptr=(void*)pick; /* this is just a flag */
	
      while(ListIterate(I->Obj,rec,next))
        {
		  glLineWidth(3.0); /* insure some 100% pixels */
          glPushMatrix();
          if(rec->obj->fRender)
            rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,&pickVLA,0);
          glPopMatrix();
        }

      
	  highBits = SceneFindTriplet(x,y);
      index = lowBits+(highBits<<12);

	  /*
	  p_glutSwapBuffers();
	  PSleep(5000000);*/

      if(index&&(index<=pickVLA[0].index)) {
        *pick = pickVLA[index]; /* return object info */
      } else {
        pick->ptr = NULL;
      }
		VLAFree(pickVLA);
		
    } else if(smp) {

      switch(I->StereoMode) {
      case 2:
        glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
        break;
      }

      /* multiple atom picking HACK - even more obfuscative coding */

      glClearColor(0.0,0.0,0.0,0.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA=VLAlloc(Pickable,5000);
      pickVLA[0].index=0;
      pickVLA[0].ptr=NULL;
      while(ListIterate(I->Obj,rec,next))
        {
		  glLineWidth(3.0); /* insure some 100% pixels */
          glPushMatrix();
			 if(rec->obj->fRender)
			   rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,&pickVLA,0);
			 glPopMatrix();
        }

	
      lowBitVLA = SceneReadTriplets(smp->x,smp->y,smp->w,smp->h);

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA[0].index=0;
      pickVLA[0].ptr=(void*)smp; /* this is just a flag */
	
      while(ListIterate(I->Obj,rec,next))
        {
		  glLineWidth(3.0); /* insure some 100% pixels */
          glPushMatrix();
          if(rec->obj->fRender)
            rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,&pickVLA,0);
          glPopMatrix();
        }
      
      highBitVLA = SceneReadTriplets(smp->x,smp->y,smp->w,smp->h);
      
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
                  nPick++; /* start from 1 */
                  VLACheck(smp->picked,Pickable,nPick);
                  lastIndex=pik->index;                
                  lastPtr=pik->ptr;
                  smp->picked[nPick] = *pik; /* return atom/object info -- will be redundant */
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
      

      ButModeCaptionReset(); /* reset the frame caption if any */
      /* rendering for visualization */

      start_time = UtilGetSeconds();
      if(I->StereoMode) {
        /*stereo*/
        
        /* render left side */

        switch(I->StereoMode) {
        case 1:
          glDrawBuffer(GL_BACK_LEFT);
          break;
        case 2:
          glViewport(I->Block->rect.left+I->Width/2,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        }
          
        glPushMatrix();
        ScenePrepareMatrix(1);
        for(pass=1;pass>=0;pass--) { /* render opaque then antialiased...*/
          rec=NULL;
          while(ListIterate(I->Obj,rec,next))
            {
              glPushMatrix();
              glNormal3fv(normal);
              if(rec->obj->fRender)
                rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,NULL,pass);
              glPopMatrix();
            }
        }
        glPushMatrix();
        glNormal3fv(normal);
        CGORenderGL(DebugCGO,NULL,NULL,NULL);
        glPopMatrix();


        glPushMatrix();
        glNormal3fv(normal);
        ExecutiveRenderSelections(curState);
        EditorRender(curState);
        glPopMatrix();

        rec=NULL;
        while(ListIterate(I->Obj,rec,next)) /* render transparent */
          {
            glPushMatrix();
            glNormal3fv(normal);
            if(rec->obj->fRender)
              rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,NULL,-1);
            glPopMatrix();
          }
        glPopMatrix();
        
        /* render right side */

        switch(I->StereoMode) {
        case 1:
          glDrawBuffer(GL_BACK_RIGHT);
          break;
        case 2:
          glViewport(I->Block->rect.left,I->Block->rect.bottom,I->Width/2,I->Height);
          break;
        }

        glClear(GL_DEPTH_BUFFER_BIT);        
        glPushMatrix();
        ScenePrepareMatrix(2);
        for(pass=1;pass>=0;pass--) { /* render opaque then antialiased...*/
          rec=NULL;
          while(ListIterate(I->Obj,rec,next))
            {
              glPushMatrix();
              glNormal3fv(normal);
              if(rec->obj->fRender)
                rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,NULL,pass);
              glPopMatrix();
            }
        }

        glPushMatrix();
        glNormal3fv(normal);
        CGORenderGL(DebugCGO,NULL,NULL,NULL);
        glPopMatrix();

        glPushMatrix();
        glNormal3fv(normal);
        ExecutiveRenderSelections(curState);
        EditorRender(curState);
        glPopMatrix();

        rec=NULL;
        while(ListIterate(I->Obj,rec,next)) /* render transparent */
          {
            glPushMatrix();
            glNormal3fv(normal);
            if(rec->obj->fRender)
              rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,NULL,-1);
            glPopMatrix();
          }

        glPopMatrix();        

        switch(I->StereoMode) {
        case 1:
          glDrawBuffer(GL_BACK);
          break;
        }


      } else {

        PRINTFD(FB_Scene)
          " SceneRender: rendering opaque and antialiased...\n"
          ENDFD;
        
        /* mono */
        for(pass=1;pass>=0;pass--) { /* render opaque then antialiased...*/
          rec=NULL;
          while(ListIterate(I->Obj,rec,next))
            {
              glPushMatrix();
              glNormal3fv(normal);
              if(rec->obj->fRender)
                rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,NULL,pass);
              glPopMatrix();
            }
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

        rec=NULL;
        while(ListIterate(I->Obj,rec,next))/* render transparent */
          {
            glPushMatrix();
            glNormal3fv(normal);
            if(rec->obj->fRender)
              rec->obj->fRender(rec->obj,SceneGetObjState(rec->obj,curState),NULL,NULL,pass); 
            glPopMatrix();
          }

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
  
  if(!(pick||smp)) {
    I->RenderTime = -I->LastRender;
    I->LastRender = UtilGetSeconds();
    I->RenderTime += I->LastRender;
    ButModeSetRate(I->RenderTime);
    if(I->CopyNextFlag) {
      start_time = I->LastRender - start_time;
      if((start_time>0.10)||(MainSavingUnderWhileIdle()))
        if(!(ControlIdling()))
          if(SettingGet(cSetting_cache_display))
            SceneCopy(0);
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
    glTranslatef(I->Pos[0],I->Pos[1],I->Pos[2]);
  
    /* rotate about the origin (the the center of rotation) */
    glMultMatrixf(I->RotMatrix);			
  
    /* move the origin to the center of rotation */
    glTranslatef(-I->Origin[0],-I->Origin[1],-I->Origin[2]);

  } else {

    /* stereo */

    stAng = SettingGet(cSetting_stereo_angle);
    stShift = SettingGet(cSetting_stereo_shift);

    stShift = stShift*fabs(I->Pos[2])/100.0;

    stAng = stAng*atan(stShift/fabs(I->Pos[2]))*90.0/PI;

    if(mode==2) {
      stAng=-stAng;
      stShift=-stShift;
    }

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
  angle = -PI*angle/180;
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










