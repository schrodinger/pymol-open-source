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
#include<stdlib.h>
#include<GL/gl.h>


#include"main.h"
#include"Rep.h"
#include"MemoryDebug.h"

/*========================================================================*/

void RepRenderBox(struct Rep *this,CRay *ray,Pickable **pick);
void RepInvalidate(struct Rep *I,int level);
void RepUpdate(struct Rep *I,struct CoordSet *cs);

/*========================================================================*/
void RepUpdate(struct Rep *I,struct CoordSet *cs)
{
  if(I->MaxInvalid>=cRepInvColor) {
	 if(I->fRecolor)
		I->fRecolor(I,cs);
  }
  I->MaxInvalid=0;
}
/*========================================================================*/
void RepInvalidate(struct Rep *I,int level)
{
  if(level>I->MaxInvalid) I->MaxInvalid=level;
}
/*========================================================================*/
void RepInit(Rep *I)
{
  I->fInvalidate = RepInvalidate;
  I->fUpdate = RepUpdate;
  I->fRender = RepRenderBox;
  I->P=NULL;
  I->MaxInvalid = 0;
}
/*========================================================================*/
void RepFree(Rep *I)
{
  FreeP(I->P);
}
/*========================================================================*/
void RepRenderBox(struct Rep *this,CRay *ray,Pickable **pick)
{
  if(PMGUI) {
    glBegin(GL_LINE_LOOP);
    glVertex3i(-0.5,-0.5,-0.5);
    glVertex3i(-0.5,-0.5, 0.5);
    glVertex3i(-0.5, 0.5, 0.5);
    glVertex3i(-0.5, 0.5,-0.5);
    
    glVertex3i( 0.5, 0.5,-0.5);
    glVertex3i( 0.5, 0.5, 0.5);
    glVertex3i( 0.5,-0.5, 0.5);
    glVertex3i( 0.5,-0.5,-0.5);
    glEnd();
    
    glBegin(GL_LINES);
    glVertex3i(0,0,0);
    glVertex3i(1,0,0);
    
    glVertex3i(0,0,0);
    glVertex3i(0,2,0);
    
    glVertex3i(0,0,0);
    glVertex3i(0,0,3);
    
    glEnd();
  }

}





