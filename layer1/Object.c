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
#include <GL/glut.h>
#include<string.h>

#include"main.h"
#include"Object.h"
#include"Color.h"

int ObjectGetNFrames(Object *I);

void ObjectDescribeElement(struct Object *I,int index);

/*========================================================================*/
void ObjectDescribeElement(struct Object *I,int index)
{
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
void ObjectFree(Object *I)
{
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
  I->fFree = ObjectFree;
  I->fRender = ObjectRenderUnitBox;
  I->fUpdate = ObjectUpdate;
  I->fGetNFrame = ObjectGetNFrames;
  I->fDescribeElement = ObjectDescribeElement;
  I->Name[0]=0;
  I->Color=1;
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




