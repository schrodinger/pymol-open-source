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
#include"Object.h"
#include"Color.h"
#include"Ortho.h"

int ObjectGetNFrames(Object *I);

void ObjectDescribeElement(struct Object *I,int index);

/*========================================================================*/
void ObjectDescribeElement(struct Object *I,int index)
{
}
/*========================================================================*/
void ObjectSetRepVis(Object *I,int rep,int state)
{
  if((rep>=0)&&(rep<cRepCnt))
    I->RepVis[rep]=state;
}
/*========================================================================*/
void ObjectSetName(Object *I,char *name)
{
  strcpy(I->Name,name);
}
/*========================================================================*/
void ObjectRenderUnitBox(struct Object *this,int frame,CRay *ray,Pickable **pick);
void ObjectUpdate(struct Object *I);

/*========================================================================*/
void ObjectUpdate(struct Object *I)
{
  
}
/*========================================================================*/
void ObjectPurge(Object *I)
{
  if(I) 
    SettingFreeP(I->Setting);
}
/*========================================================================*/
void ObjectFree(Object *I)
{
  if(I)
    ObjectPurge(I);
}
/*========================================================================*/
int ObjectGetNFrames(Object *I)
{
  return 1;
}
/*========================================================================*/
void ObjectUseColor(Object *I)
{
  if(PMGUI) glColor3fv(ColorGet(I->Color));
}
/*========================================================================*/
void ObjectInit(Object *I)
{
  int a;
  I->fFree = ObjectFree;
  I->fRender = ObjectRenderUnitBox;
  I->fUpdate = ObjectUpdate;
  I->fGetNFrame = ObjectGetNFrames;
  I->fDescribeElement = ObjectDescribeElement;
  I->Name[0]=0;
  I->Color=0;
  I->ExtentFlag=false;
  I->Setting=NULL;
  OrthoRemoveSplash();
  for(a=0;a<cRepCnt;a++) I->RepVis[a]=true;
}
/*========================================================================*/
void ObjectRenderUnitBox(Object *this,int frame,CRay *ray,Pickable **pick)
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




