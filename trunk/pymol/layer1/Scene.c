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


#include"os_std.h"
#include"os_gl.h"

#include"Util.h"

#include"Word.h"
#include"main.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Matrix.h"
#include"ListMacros.h"
#include"Object.h"
#include"Scene.h"
#include"Ortho.h"
#include"Vector.h"
#include"ButMode.h"
#include"Control.h"
#include"Selector.h"
#include"Setting.h"
#include"Movie.h"
#include"MyPNG.h"
#include"Python.h"
#include"P.h"


#define cFrontMin 0.1
#define cSliceMin 0.1

#define SceneLineHeight 12
#define SceneTopMargin 0
#define SceneBottomMargin 3
#define SceneLeftMargin 3

#define SceneFOV 20.0

typedef struct ObjRec {
  struct Object *obj;  
  struct ObjRec *next;
} ObjRec;

ListVarDeclare(ObjList,ObjRec);

typedef struct {
  Block *Block;
  ObjRec *Obj;
  float RotMatrix[16];
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
  int StateIndex,Frame,NFrame;
  GLvoid *ImageBuffer;
  int MovieOwnsImageFlag;
  int MovieFrameFlag;
  unsigned ImageBufferSize;
  double LastRender,RenderTime,LastFrameTime;
  float LastRock,LastRockTime;
  Pickable LastPicked;
  int StereoMode;
  
} CScene;

CScene Scene;

unsigned int SceneFindTriplet(int x,int y);
void SceneDraw(Block *block);
int SceneClick(Block *block,int button,int x,int y,int mod);
int SceneDrag(Block *block,int x,int y,int mod);
void ScenePrepareMatrix(int mode);

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
  I->StereoMode=flag;
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
void SceneClip(int plane,float movement)
{
  CScene *I=&Scene;
  if(plane) {
	 I->Back-=movement;
  } else {
	 I->Front-=movement;
  }
  if(I->Front>I->Back)
	 I->Front=I->Back+cSliceMin;
  if(I->Front<cFrontMin) I->Front=cFrontMin;
  I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
  SceneDirty();

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
    }
  } else {
	 image=I->ImageBuffer;
  }
  MyPNGWrite(png,image,I->Width,I->Height);
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
  CScene *I=&Scene;
  
  return(I->Frame);
}
/*========================================================================*/
void SceneCountFrames() 
{
  CScene *I=&Scene;
  ObjRec *rec = NULL;
  int n;
  int frame;

  I->NFrame=0;
  while(ListIterate(I->Obj,rec,next,ObjList))
	 {
      if(rec->obj->fGetNFrame)
        n=rec->obj->fGetNFrame(rec->obj);
      else
        n=0;
		if(n>I->NFrame)
		  I->NFrame=n;
	 }
  if(I->NFrame<MovieGetLength())
	 I->NFrame=MovieGetLength();
  if(I->Frame>=I->NFrame) {
	 frame=I->NFrame-1;
	 if(frame<0) frame=0;
	 SceneSetFrame(0,frame);
  }
}
/*========================================================================*/
void SceneSetFrame(int mode,int frame)
{
  CScene *I=&Scene;
  switch(mode) {
  case 0:
	 I->Frame=frame;
	 break;
  case 1:
	 I->Frame+=frame;
	 break;
  case 2:
	 I->Frame=I->NFrame-1;
	 break;
  case 3:
	 I->Frame=I->NFrame/2;
	 break;
  case 4:
	 I->Frame=frame;
	 break;
  case 5:
	 I->Frame+=frame;
	 break;
  }
  if(I->Frame>=I->NFrame) I->Frame=I->NFrame-1;
  if(I->Frame<0) I->Frame=0;
  I->StateIndex = MovieFrameToIndex(I->Frame);
  if(mode&4) 
	MovieDoFrameCommand(I->Frame);
  if(I->Frame==0)
	MovieMatrix(cMovieMatrixRecall);
  if(SettingGet(cSetting_cache_frames))
	 I->MovieFrameFlag=true;
  SceneDirty();
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
  I->DirtyFlag=true;
  ScenePurgeCopy();
  OrthoDirty();
}
/*========================================================================*/
void SceneChanged(void)
{
  CScene *I=&Scene;
  I->ChangedFlag=true;
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
	SceneRay(); 
  } else {
	 v=SettingGetfv(cSetting_bg_rgb);
    if(PMGUI) {
      glDrawBuffer(GL_BACK);
      glClearColor(v[0],v[1],v[2],1.0);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glClearColor(0.0,0.0,0.0,1.0);
      SceneRender(NULL,0,0);
      SceneCopy(0);
    }
  }
  if(I->ImageBuffer)
	 MovieSetImage(MovieFrameToImage(I->Frame),I->ImageBuffer);
  I->MovieOwnsImageFlag=true;
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
	 SceneRotate(-I->LastRock,0.0,1.0,0.0);
	 I->RockTime+=I->RenderTime;
	 I->LastRock=sin(I->RockTime*SettingGet(cSetting_sweep_speed))*SettingGet(cSetting_sweep_angle);
	 SceneRotate(I->LastRock,0.0,1.0,0.0);
  }
  if(MoviePlaying()&&frameFlag)
	 {
      I->LastFrameTime = UtilGetSeconds();

      if(I->Frame==I->NFrame-1)
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
  /* find where this point is in relationship to the origin */
  subtract3f(I->Origin,location,v0); 
  /*  printf("%8.3f %8.3f %8.3f\n",I->Front,I->Pos[2],I->Back);*/

  MatrixTransform3f(v0,I->RotMatrix,I->Pos); /* convert to view-space */
  dist = radius/tan((SceneFOV/2.0)*cPI/180.0);

  I->Pos[2]-=dist;
  I->Front=(-I->Pos[2]-radius*1.5);
  I->FrontSafe=(I->Front<cFrontMin ? cFrontMin : I->Front);  
  I->Back=(-I->Pos[2]+radius*1.5);
  /*  printf("%8.3f %8.3f %8.3f\n",I->Front,I->Pos[2],I->Back);*/
}
/*========================================================================*/
void SceneOriginSet(float *origin,int preserve)
{
  CScene *I=&Scene;
  float v0[3],v1[3];
  
  if(preserve) /* preserve current viewing location */
	 {
		subtract3f(origin,I->Origin,v0); /* model-space translation */
		MatrixTransform3f(v0,I->RotMatrix,v1); /* convert to view-space */
		add3f(I->Pos,v1,I->Pos); /* offset view to compensate */
	 }
  I->Origin[0]=origin[0]; /* move origin */
  I->Origin[1]=origin[1];
  I->Origin[2]=origin[2];
  SceneDirty();
}
/*========================================================================*/
void SceneObjectAdd(Object *obj)
{
  CScene *I=&Scene;
  ObjRec *rec = NULL;
  ListElemAlloc(rec,ObjRec);
  rec->next=NULL;
  rec->obj=obj;
  ListAppend(I->Obj,rec,next,ObjList);
  SceneCountFrames();
  SceneDirty();
}
/*========================================================================*/
void SceneObjectDel(Object *obj)
{
  CScene *I=&Scene;
  ObjRec *rec = NULL;

  while(ListIterate(I->Obj,rec,next,ObjList))
	 if(rec->obj==obj)
		break;
  if(rec) {
	 ListDetach(I->Obj,rec,next,ObjList);
	 ListElemFree(rec);
  }
  SceneCountFrames();
  SceneDirty();
}
/*========================================================================*/
void SceneDraw(Block *block)
{
  CScene *I=&Scene;
  int overlay,text;

  if(PMGUI) {
    overlay = SettingGet(cSetting_overlay);
    text = SettingGet(cSetting_text);

    if(overlay||(!text)) 

      if(I->CopyFlag)
        {
          glReadBuffer(GL_BACK);
          glRasterPos3i(I->Block->rect.left,I->Block->rect.bottom,0);
          glDrawPixels(I->Width,I->Height,GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);
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
#define cRange 10

unsigned int SceneFindTriplet(int x,int y) 
{
  int result = 0;
  pix buffer[cRange*2+1][cRange*2+1];
  int a,b,d,e,flag;
  unsigned char *c;
  
  if(PMGUI) { /*just in case*/
    glReadBuffer(GL_BACK);
    glReadPixels(x-cRange,y-cRange,cRange*2+1,cRange*2+1,GL_RGBA,GL_UNSIGNED_BYTE,&buffer[0][0][0]);
    
    /*  for(a=0;a<=(cRange*2);a++)
        {
        for(b=0;b<=(cRange*2);b++)
		  printf("%2x ",(buffer[a][b][0]+buffer[a][b][1]+buffer[a][b][2])&0xFF);
        printf("\n");
        }
        printf("\n");	 
    */
    flag=true;
    for(d=0;flag&&(d<cRange);d++)
      for(a=-d;flag&&(a<=d);a++)
        for(b=-d;flag&&(b<=d);b++)
          {
            for(e=0;e<3;e++)
              if(buffer[a+cRange][b+cRange][e]) {
                flag = false;
                c=&buffer[a+cRange][b+cRange][0];
                /*  printf("%2x %2x %2x\n",c[0],c[1],c[2]);*/
                result =  ((c[0]>>4)&0xF)+(c[1]&0xF0)+((c[2]<<4)&0xF00);
                break;
              }
          }
  }
  return(result);
}
/*========================================================================*/
int SceneClick(Block *block,int button,int x,int y,int mod)
{
  CScene *I=&Scene;
  Object *obj;
  char buffer[OrthoLineLength],buf2[OrthoLineLength];
  WordType selName = "";
  int mode;

  mode = ButModeTranslate(button,mod);
  switch(mode) {
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

	 I->LastX=x;
	 I->LastY=y;	 
	 SceneDirty();
	 I->Button=button;    
    break;
  case cButModePk1:
  case cButModePk2:
  case cButModePk3:
  case cButModeAddToPk1:
  case cButModeAddToPk2:
  case cButModeAddToPk3:
    if(((int)SettingGet(cSetting_overlay))&&((int)SettingGet(cSetting_text)))
      SceneRender(NULL,0,0);
    SceneDontCopyNext();

	 SceneRender(&I->LastPicked,x,y);
	 if(I->LastPicked.ptr) {
		obj=(Object*)I->LastPicked.ptr;
      if(obj->fDescribeElement) 
        obj->fDescribeElement(obj,I->LastPicked.index);
		  sprintf(buffer,"model %s and index %i",
					 obj->Name,I->LastPicked.index+1);
		switch(mode) {
      case cButModePk1:
      case cButModeAddToPk1:
        strcpy(selName,"%pk1");
		  break;
      case cButModePk2:
      case cButModeAddToPk2:
        strcpy(selName,"%pk2");
		  break;
      case cButModePk3:
      case cButModeAddToPk3:
        strcpy(selName,"%pk3");
		  break;
      }
      switch(mode) {
      case cButModePk1:
      case cButModePk2:
      case cButModePk3:
        SelectorCreate(selName,buffer,NULL,false);
        break;
      case cButModeAddToPk1:
      case cButModeAddToPk2:
      case cButModeAddToPk3:
        if(SelectorIndexByName(selName)>=0) {
          sprintf(buf2,"( %s or (%s))",selName,buffer);
          SelectorCreate(selName,buf2,NULL,false);
        } else 
          SelectorCreate(selName,buffer,NULL,false);
        break;
      }
	 } else {
		OrthoAddOutput(" SceneClick: no atom found nearby.\n");
		OrthoNewLine(NULL);
		OrthoRestorePrompt();
	 }
  }
  return(1);
}
/*========================================================================*/
int SceneDrag(Block *block,int x,int y,int mod)
{
  CScene *I=&Scene;
  float scale;
  float v1[3],v2[3],n1[3],n2[3],r1,r2,cp[3];
  float axis[3],axis2[3],theta,omega;
  int mode;

  mode = ButModeTranslate(I->Button,mod);
  
  y=y-I->Block->margin.bottom;
  scale = I->Height;
  if(scale > I->Width)
	 scale = I->Width;
  scale = 0.38 * scale;

  SceneDontCopyNext();
  switch(mode) {
  case cButModeRotXYZ:
  case cButModeRotZ:
  case cButModeTransXY:
  case cButModeTransZ:
  case cButModeClipNF:
  case cButModeClipN:    
  case cButModeClipF:    

    v1[0] = (I->Width/2) - x;
    v1[1] = (I->Height/2) - y;
    
	 v2[0] = (I->Width/2) - I->LastX;
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
	 case cButModeTransXY:
		if(I->LastX!=x)
		  {
			 I->Pos[0]+=(((float)x)-I->LastX)/10;
			 I->LastX=x;
			 SceneDirty();
		  }
		if(I->LastY!=y)
		  {
			 I->Pos[1]+=(((float)y)-I->LastY)/10;
			 I->LastY=y;
			 SceneDirty();
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
  
  if(I->StereoMode) {
    PStereoOff();
  }
}
/*========================================================================*/
void SceneResetMatrix(void)
{
  CScene *I=&Scene;
  MatrixLoadIdentity44f(I->RotMatrix);
}
/*========================================================================*/
void SceneInit(void)
{
  CScene *I=&Scene;

  ListInit(I->Obj);

  I->RockTime=0;
  I->TextColor[0]=0.2;
  I->TextColor[1]=1.0;
  I->TextColor[2]=0.2;

  MatrixLoadIdentity44f(I->RotMatrix);

  I->Scale = 1.0;
  I->Frame=0;
  I->StateIndex=0;
  
  I->Front=40;
  I->FrontSafe= (I->Front<cFrontMin ? cFrontMin : I->Front);
  I->Back=100;
  
  I->Block = OrthoNewBlock(NULL);
  I->Block->fClick   = SceneClick;
  I->Block->fDrag    = SceneDrag;
  I->Block->fDraw    = SceneDraw;
  I->Block->fReshape = SceneReshape;
  I->Block->active = true;

  I->Pos[0] = 0.0;
  I->Pos[1] = 0.0;
  I->Pos[2] = -50.0;

  I->Origin[0] = 0.0;
  I->Origin[1] = 0.0;
  I->Origin[2] = 0.0;

  I->Scale = 1.0;

  OrthoAttach(I->Block,cOrthoScene);

  I->DirtyFlag = true;
  I->ImageBuffer = NULL;
  I->ImageBufferSize = 0;
  I->MovieOwnsImageFlag = false;
  I->MovieFrameFlag = false;
  I->RenderTime = 0;
  I->LastRender = UtilGetSeconds();
  I->LastFrameTime = UtilGetSeconds();
  I->LastRockTime = UtilGetSeconds();
  I->LastPicked.ptr = NULL;

  I->CopyNextFlag=true;
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
void SceneRay(void)
{
  CScene *I=&Scene;
  ObjRec *rec=NULL;
  CRay *ray;
  unsigned int buffer_size;
  float height,width;
  float aspRat = ((float) I->Width) / ((float) I->Height);
  float white[3] = {1.0,1.0,1.0};
  unsigned int *buffer;
  float rayView[16];
  int curState;

  if(SettingGet(cSetting_all_states)) {
    curState=-1;
  } else {
    curState=I->StateIndex;
  }

  ray = RayNew();

  SceneUpdate();

  
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

  /* define the viewing volume */

  height  = abs(I->Pos[2])*tan((SceneFOV/2.0)*cPI/180.0);	 
  width = height*aspRat;

  RayPrepare(ray,-width,width,-height,height,I->FrontSafe,I->Back,rayView);

  while(ListIterate(I->Obj,rec,next,ObjList))
	 {
		if(rec->obj->fRender) {
		  ray->fColor3fv(ray,white);
		  rec->obj->fRender(rec->obj,curState,ray,NULL);
		}
	 }

  buffer_size = 4*I->Width*I->Height;
  buffer=(GLvoid*)Alloc(char,buffer_size);
  ErrChkPtr(buffer);

  RayRender(ray,I->Width,I->Height,buffer,I->Front,I->Back);

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
  I->DirtyFlag=false;
  I->CopyFlag = true;
  I->MovieOwnsImageFlag = false;

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
      if(!I->ImageBuffer) {
        I->ImageBuffer=(GLvoid*)Alloc(char,buffer_size);
        ErrChkPtr(I->ImageBuffer);
        I->ImageBufferSize = buffer_size;
      }
      if(PMGUI) {
        if(buffer)
          glReadBuffer(GL_FRONT);
        else
          glReadBuffer(GL_BACK);
        glReadPixels(I->Block->rect.left,I->Block->rect.bottom,I->Width,I->Height,
                     GL_RGBA,GL_UNSIGNED_BYTE,I->ImageBuffer);
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

  if(I->ChangedFlag) {
    SceneCountFrames();
	 while(ListIterate(I->Obj,rec,next,ObjList))
      if(rec->obj->fUpdate) 
        rec->obj->fUpdate(rec->obj);
	 I->ChangedFlag=false;
  }
}
/*========================================================================*/
int SceneRenderCached(void)
{
  /* sets up a cached image buffer is one is available, or if we are
   * using cached images by default */
  CScene *I=&Scene;
  ImageType image;
  int renderedFlag=false;

  if(I->DirtyFlag) {
	if(I->MovieFrameFlag||
	   (MoviePlaying()&&SettingGet(cSetting_cache_frames))) {
	  I->MovieFrameFlag=false;
	  image = MovieGetImage(MovieFrameToImage(I->Frame));
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
	  SceneRay(); 
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
  return(renderedFlag);
}
/*========================================================================*/
void SceneRender(Pickable *pick,int x,int y)
{
  /* think in terms of the camera's world */
  CScene *I=&Scene;
  ObjRec *rec=NULL;
  float fog[4];
  float *v,vv[4],f;
  unsigned int lowBits,highBits;
  static float white[4] =
  {1.0, 1.0, 1.0, 1.0};
  float zAxis[4] = { 0.0, 0.0, 1.0, 0.0 };
  float normal[4] = { 0.0, 0.0, 1.0, 0.0 };
  float aspRat = ((float) I->Width) / ((float) I->Height);
  float height,width;
  float start_time=0;
  int view_save[4];
  Pickable *pickVLA;
  int index;
  int curState;

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
      curState=I->StateIndex;
    }

    if(SettingGet(cSetting_ortho)==0.0) {
      gluPerspective(SceneFOV,aspRat,I->FrontSafe,I->Back);
    } else {
      height  = abs(I->Pos[2])*tan((SceneFOV/2.0)*cPI/180.0);	 
      width = height*aspRat;
	
      glOrtho(-width,width,-height,height,
              I->FrontSafe,I->Back);
    }

    glMatrixMode(GL_MODELVIEW);
    ScenePrepareMatrix(0);
  
    /* determine the direction in which we are looking relative*/

    /* 2. set the normals to reflect light back at the camera */

    MatrixInvTransform3f(zAxis,I->RotMatrix,normal); 
    I->ViewNormal[0]=normal[0];
    I->ViewNormal[1]=normal[1];
    I->ViewNormal[2]=normal[2];	 
  
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

    if(SettingGet(cSetting_line_smooth)) {
      if(!pick) {
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
      }
      glLineWidth(0.0);
    } else {
      glLineWidth(SettingGet(cSetting_line_width));
      glDisable(GL_LINE_SMOOTH);
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    
    if(!pick) {

      glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, white);
      glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_DITHER);

      v=SettingGetfv(cSetting_ambient);
      f=SettingGet(cSetting_ambient_scale);
      vv[0]=v[0]*f;
      vv[1]=v[0]*f;
      vv[2]=v[0]*f;
      vv[3]=1.0;

      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT,vv);
      /*glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);*/
	
#ifdef _PYMOL_3DFX
      if(SettingGet(cSetting_ortho)==0.0) {
#endif

        glEnable(GL_FOG);
        glFogf(GL_FOG_MODE, GL_LINEAR);
        glHint(GL_FOG_HINT,GL_NICEST);
        glFogf(GL_FOG_START, I->FrontSafe);
        glFogf(GL_FOG_END, I->Back);
#ifdef _PYMOL_3DFX
        if(I->Back>(I->FrontSafe*4.0))
          glFogf(GL_FOG_END, I->Back);
        else
          glFogf(GL_FOG_END,I->FrontSafe*4.0);
        fog_val+=0.0000001;
        if(fog_val>1.0) fog_val=0.99999;
        glFogf(GL_FOG_DENSITY, fog_val);
#else
        glFogf(GL_FOG_END,I->Back);
        glFogf(GL_FOG_DENSITY, 1.0);
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

      glColor4ub(255,255,255,255);
      glNormal3fv(normal);
      
    } else {
      /* picking mode: we want flat, unshaded colors */

      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_LIGHTING);
      glDisable(GL_DITHER);
    }

    /* 1. render all objects */
    if(pick) {
      /* atom picking HACK - obfuscative coding */
	
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA=VLAlloc(Pickable,1000);
      pickVLA[0].index=0;
      pickVLA[0].ptr=NULL;
      while(ListIterate(I->Obj,rec,next,ObjList))
        {
          glPushMatrix();
			 if(rec->obj->fRender)
			   rec->obj->fRender(rec->obj,curState,NULL,&pickVLA);
			 glPopMatrix();
        }
	
      lowBits = SceneFindTriplet(x,y);
	
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
      pickVLA[0].index=0;
      pickVLA[0].ptr=(void*)pick; /* this is just a flag */
	
      while(ListIterate(I->Obj,rec,next,ObjList))
        {
          glPushMatrix();
          if(rec->obj->fRender)
            rec->obj->fRender(rec->obj,curState,NULL,&pickVLA);
          glPopMatrix();
        }
      highBits = SceneFindTriplet(x,y);
      index = lowBits+(highBits<<12);
	
      if(index&&(index<=pickVLA[0].index)) {
        *pick = pickVLA[index]; /* return object info */
      } else {
        pick->ptr = NULL;
      }
		VLAFree(pickVLA);
		
    } else {
      ButModeCaptionReset(); /* reset the frame caption if any */
      /* rendering for visualization */

      start_time = UtilGetSeconds();
      if(I->StereoMode) {
        /*stereo*/

        glDrawBuffer(GL_BACK_LEFT);
        glPushMatrix();
        ScenePrepareMatrix(1);
        rec=NULL;
        while(ListIterate(I->Obj,rec,next,ObjList))
          {
            glPushMatrix();
            glNormal3fv(normal);
            if(rec->obj->fRender)
              rec->obj->fRender(rec->obj,curState,NULL,NULL);
            glPopMatrix();
          }
        glPopMatrix();
        
        glDrawBuffer(GL_BACK_RIGHT);
        glClear(GL_DEPTH_BUFFER_BIT);        
        glPushMatrix();
        ScenePrepareMatrix(2);
        rec=NULL;
        while(ListIterate(I->Obj,rec,next,ObjList))
          {
            glPushMatrix();
            glNormal3fv(normal);
            if(rec->obj->fRender)
              rec->obj->fRender(rec->obj,curState,NULL,NULL);
            glPopMatrix();
          }
        glPopMatrix();        glDrawBuffer(GL_BACK);
      } else {
        
        /* mono */
        rec=NULL;
        while(ListIterate(I->Obj,rec,next,ObjList))
          {
            glPushMatrix();
            glNormal3fv(normal);
            if(rec->obj->fRender)
              rec->obj->fRender(rec->obj,curState,NULL,NULL);
            glPopMatrix();
          }
      }
    }
  
    if(!pick) {
      glDisable(GL_FOG);
      glDisable(GL_LIGHTING);
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
  if(!pick) {
    I->RenderTime = -I->LastRender;
    I->LastRender = UtilGetSeconds();
    I->RenderTime += I->LastRender;
    ButModeSetRate(I->RenderTime);
    if(I->CopyNextFlag) {
      start_time = I->LastRender - start_time;
      if((start_time>0.10)||(MainSavingUnderWhileIdle()))
        if(!(ControlIdling()))
          SceneCopy(0);
    } else {
      I->CopyNextFlag=true;
    }
  }
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










