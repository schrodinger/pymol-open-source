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

#include"os_gl.h"
#include"os_std.h"

#include"main.h"
#include"PyMOLObject.h"
#include"Color.h"
#include"Ortho.h"
#include"Scene.h"
#include"Util.h"
#include"Ray.h"

int ObjectGetNFrames(CObject *I);

void ObjectDescribeElement(struct CObject *I,int index,char *buffer);
CSetting **ObjectGetSettingHandle(struct CObject *I,int state);

/*========================================================================*/
void ObjectCombineTTT(CObject *I,float *ttt)
{
  float cpy[16];
  if(!I->TTTFlag) {
    I->TTTFlag=true;
    initializeTTT44f(cpy);
  } else {
    UtilCopyMem(cpy,I->TTT,sizeof(float)*16);
  }
  combineTTT44f44f(ttt,cpy,I->TTT);
}
/*========================================================================*/
void ObjectResetTTT(CObject *I)
{
  I->TTTFlag=false;
  SceneDirty();
}
/*========================================================================*/
void ObjectPrepareContext(CObject *I,CRay *ray)
{
  float gl[16],*ttt;
  if(PMGUI) {
    if(ray) {
      RaySetTTT(ray,I->TTTFlag,I->TTT);
    } else {
      if(I->TTTFlag) {
        /* form standard 4x4 GL matrix with TTT rotation and 2nd translation */
        ttt=I->TTT;
        gl[ 0] = ttt[ 0];
        gl[ 4] = ttt[ 1];
        gl[ 8] = ttt[ 2];
        gl[12] = ttt[ 3];
        gl[ 1] = ttt[ 4];
        gl[ 5] = ttt[ 5];
        gl[ 9] = ttt[ 6];
        gl[13] = ttt[ 7];
        gl[ 2] = ttt[ 8];
        gl[ 6] = ttt[ 9];
        gl[10] = ttt[10];
        gl[14] = ttt[11];
        gl[ 3] = 0.0;
        gl[ 7] = 0.0;
        gl[11] = 0.0;
        gl[15] = 1.0;
        /*        dump44f(gl,"ttt");
                  dump44f(gl,"gl");*/
        glMultMatrixf(gl);
        /* now add in the first translation */
        glTranslatef(ttt[12],ttt[13],ttt[14]);
      }
    }
  }
}

/*========================================================================*/
void ObjectSetTTTOrigin(CObject *I,float *origin)
{
  if(!I->TTTFlag) {
    I->TTTFlag=true;
    initializeTTT44f(I->TTT);
  }

  I->TTT[3]+=I->TTT[12]; /* remove existing origin from overall translation */
  I->TTT[7]+=I->TTT[13];
  I->TTT[11]+=I->TTT[14];

  scale3f(origin,-1.0,I->TTT+12); /* set new origin */

  I->TTT[3]+=origin[0]; /* add new origin into overall translation */
  I->TTT[7]+=origin[1];
  I->TTT[11]+=origin[2];

  SceneDirty();

}
/*========================================================================*/
CSetting **ObjectGetSettingHandle(struct CObject *I,int state)
{
  return(&I->Setting);
}
/*========================================================================*/
void ObjectDescribeElement(struct CObject *I,int index,char *buffer)
{
  buffer[0]=0;
}
/*========================================================================*/
void ObjectSetRepVis(CObject *I,int rep,int state)
{
  if((rep>=0)&&(rep<cRepCnt))
    I->RepVis[rep]=state;
}
/*========================================================================*/
void ObjectSetName(CObject *I,char *name)
{
  strcpy(I->Name,name);
}
/*========================================================================*/
void ObjectRenderUnitBox(struct CObject *this,int frame,CRay *ray,Pickable **pick,int pass);
void ObjectUpdate(struct CObject *I);

/*========================================================================*/
void ObjectUpdate(struct CObject *I)
{
  
}
/*========================================================================*/
void ObjectPurge(CObject *I)
{
  if(I) 
    SettingFreeP(I->Setting);
}
/*========================================================================*/
void ObjectFree(CObject *I)
{
  if(I)
    ObjectPurge(I);
}
/*========================================================================*/
int ObjectGetNFrames(CObject *I)
{
  return 1;
}
/*========================================================================*/
void ObjectUseColor(CObject *I)
{
  if(PMGUI) glColor3fv(ColorGet(I->Color));
}
/*========================================================================*/
static void ObjectInvalidate(CObject *this,int rep,int level,int state)
{
  
}
/*========================================================================*/
void ObjectInit(CObject *I)
{
  int a;
  I->fFree = ObjectFree;
  I->fRender = ObjectRenderUnitBox;
  I->fUpdate = ObjectUpdate;
  I->fGetNFrame = ObjectGetNFrames;
  I->fDescribeElement = ObjectDescribeElement;
  I->fGetSettingHandle = ObjectGetSettingHandle;
  I->fInvalidate = ObjectInvalidate;
  I->Name[0]=0;
  I->Color=0; /* white */
  I->ExtentFlag=false;
  I->Setting=NULL;
  I->TTTFlag=false;
  OrthoRemoveSplash();
  for(a=0;a<cRepCnt;a++) I->RepVis[a]=true;
  I->RepVis[cRepCell]=false;
  I->RepVis[cRepExtent]=false;
}
/*========================================================================*/
void ObjectRenderUnitBox(CObject *this,int frame,
                         CRay *ray,Pickable **pick,int pass)
{
  if(PMGUI) {
    glBegin(GL_LINE_LOOP);
    glVertex3i(-1,-1,-1);
    glVertex3i(-1,-1, 1);
    glVertex3i(-1, 1, 1);
    glVertex3i(-1, 1,-1);
    
    glVertex3i( 1, 1,-1);
    glVertex3i( 1, 1, 1);
    glVertex3i( 1,-1, 1);
    glVertex3i( 1,-1,-1);
    glEnd();
    
    glBegin(GL_LINES);
    glVertex3i(0,0,0);
    glVertex3i(1,0,0);
    
    glVertex3i(0,0,0);
    glVertex3i(0,3,0);
    
    glVertex3i(0,0,0);
    glVertex3i(0,0,9);

    glEnd();
  }
}




