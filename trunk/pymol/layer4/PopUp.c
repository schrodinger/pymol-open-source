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

#include "OOMac.h"

#include "main.h"
#include "Base.h"
#include "Pop.h"
#include "PopUp.h"
#include "Ortho.h"
#include "Util.h"
#include "P.h"

#define cPopUpLineHeight 13
#define cPopUpCharWidth 8
#define cPopUpCharMargin 2
#define cPopUpCharLift 1


typedef struct {
  Block *Block;
  int LastX,LastY;
  int Selected;
  int Width,Height;
  int NLine;
  char **Command;
  char **Text;
  int *Code;
}  CPopUp;

int PopUpRelease(Block *block,int button,int x,int y,int mod);
void PopUpDraw(Block *block);
int PopUpDrag(Block *block,int x,int y,int mod);
/*========================================================================*/
void PopUpNew(int x,int y,PyObject *list)
{
  /* assumes blocked threads (calls the Python C API) */

  int mx,a,l,cl,cmx;
  int dim[2];
  PyObject *elem;
  char *str,*c;

  OOAlloc(CPopUp);

  I->Block = OrthoNewBlock(NULL);
  I->Block->reference = (void*)I;
  I->Block->fDraw    = PopUpDraw;
  I->Block->fDrag    = PopUpDrag;
  I->Block->fRelease = PopUpRelease;
  I->Block->active = false;
  I->Block->TextColor[0]=1.0;
  I->Block->TextColor[1]=1.0;
  I->Block->TextColor[2]=1.0;

  I->Block->BackColor[0]=0.1;
  I->Block->BackColor[1]=0.1;
  I->Block->BackColor[2]=0.1;

  I->NLine=PyList_Size(list);
  I->Text = NULL;
  I->Command = NULL;
  I->Code = NULL;
  I->Selected = -1;

  mx=1;
  cmx=1;
  for(a=0;a<I->NLine;a++) {
    elem = PyList_GetItem(PyList_GetItem(list,a),1);
    l = PyString_Size(elem);
    str  = PyString_AsString(elem);
    cl=l;
    c = str;
    while(*c) {
      if(*c=='`') { /* discount the markup */
        cl-=4;
      }
      c++;
    }
    if(cl>cmx) cmx=cl;
    if(l>mx) mx=l;
  }
  I->Width = (cmx * cPopUpCharWidth) + 2 * cPopUpCharMargin;
  I->Height = (I->NLine * cPopUpLineHeight) + 2 * cPopUpCharMargin;

  dim[0]=I->NLine+1;
  dim[1]=mx+1;
  I->Text = (char**)UtilArrayMalloc((unsigned int*)dim,2,1);

  mx=1;
  for(a=0;a<I->NLine;a++) {
    l = PyString_Size(PyList_GetItem(PyList_GetItem(list,a),2));
    if(l>mx) mx=l;
  }
  dim[0]=I->NLine+1;
  dim[1]=mx+1;
  I->Command = (char**)UtilArrayMalloc((unsigned int*)dim,2,1);

  I->Code = Alloc(int,I->NLine+1);

  for(a=0;a<I->NLine;a++) {
    elem = PyList_GetItem(list,a);
    I->Code[a]=PyInt_AsLong(PyList_GetItem(elem,0));
    strcpy(I->Text[a],PyString_AsString(PyList_GetItem(elem,1)));
    strcpy(I->Command[a],PyString_AsString(PyList_GetItem(elem,2)));
  }

  I->Block->rect.top=y;
  I->Block->rect.bottom=y-I->Height;
  I->Block->rect.left=x-(I->Width)/3;
  I->Block->rect.right=x+(2*I->Width)/3;

  PopFitBlock(I->Block);

  OrthoAttach(I->Block,cOrthoTool);
  I->Block->active=true;
  OrthoGrab(I->Block);
  OrthoDirty();

}
/*========================================================================*/
int PopUpRelease(Block *block,int button,int x,int y,int mod)
{
  CPopUp *I = (CPopUp*)block->reference;
  PopUpDrag(block,x,y,mod);
  OrthoUngrab();
  OrthoDetach(I->Block);

  if(I->Selected>=0) {
    PParse(I->Command[I->Selected]);
    PFlush();
    /* PBlockAndUnlockAPI();
      PRunString(I->Command[I->Selected]);
      PLockAPIAndUnblock(); */
  }
  
  OrthoFreeBlock(I->Block);
  FreeP(I->Code);
  FreeP(I->Command);
  FreeP(I->Text);
  OOFreeP(I);
  OrthoDirty();
  return(1);
}
/*========================================================================*/
int PopUpDrag(Block *block,int x,int y,int mod)
{
  CPopUp *I = (CPopUp*)block->reference;
  
  int a;
  int was = I->Selected;
  I->LastX=x;
  I->LastY=y;

  x-=I->Block->rect.left;
  y =(I->Block->rect.top -cPopUpCharMargin) -y -1;

  if((x<0)||(x>I->Width)) 
    I->Selected=-1;
  else {
    a = y/cPopUpLineHeight;
    if(I->NLine&&(a==I->NLine))
      if((y-a*cPopUpLineHeight)<4)
        a=I->NLine-1;
    if((a<0)||(a>=I->NLine))
      I->Selected=-1;
    else if(I->Code[a]!=1)
      I->Selected=-1;
    else
      I->Selected=a;
  }
  if(was!=I->Selected)
    OrthoDirty();
  return(1);
}
/*========================================================================*/
void PopUpDraw(Block *block)
{
  CPopUp *I = (CPopUp*)block->reference;
  int x,y,a,xx;
  char *c;

  if(PMGUI) {

    glColor3fv(block->BackColor);
    BlockFill(block);
    glColor3fv(block->TextColor);

    if(I->Selected>=0) {

      x = I->Block->rect.left;
      y = I->Block->rect.top-(cPopUpLineHeight*I->Selected)-cPopUpCharMargin;
      
      glBegin(GL_POLYGON);
      glVertex2i(x,y);
      glVertex2i(x+I->Width,y);
      glVertex2i(x+I->Width,y-(cPopUpLineHeight+cPopUpCharMargin));
      glVertex2i(x,y-(cPopUpLineHeight+cPopUpCharMargin));
      glEnd();
    }

    if(I->Code[0]==2) {
        
      glColor3f(0.4,0.4,0.4);
      x = I->Block->rect.left;
      y = I->Block->rect.top;
      
      glBegin(GL_POLYGON);
      glVertex2i(x,y);
      glVertex2i(x+I->Width,y);
      glVertex2i(x+I->Width,y-(cPopUpLineHeight+cPopUpCharMargin));
      glVertex2i(x,y-(cPopUpLineHeight+cPopUpCharMargin));
      glEnd();
    }

    x = I->Block->rect.left+cPopUpCharMargin;
    y = (I->Block->rect.top-cPopUpLineHeight)-cPopUpCharMargin+1;
  
    for(a=0;a<I->NLine;a++)
      {
        if(a==I->Selected)
          glColor3fv(I->Block->BackColor);
        else
          glColor3fv(I->Block->TextColor);          
        if(I->Code[a]) {
          c=I->Text[a];
          xx=x;
          while(*c) {
            if(*c=='`') if(*(c+1)) if(*(c+2)) if(*(c+3)) {
              if(*(c+1)=='-') {
                if(a==I->Selected)
                  glColor3fv(I->Block->BackColor);
                else
                  glColor3fv(I->Block->TextColor);          
                c+=4;
              } else {
                glColor3f((*(c+1)-'0')/9.0,(*(c+2)-'0')/9.0,(*(c+3)-'0')/9.0);
                c+=4;
              }
            }
            glRasterPos4d((double)(xx),(double)(y+cPopUpCharLift),0.0,1.0);
            p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c++));
            xx = xx + 8;
          }
        } else {
          glBegin(GL_LINES);
          glVertex2i(I->Block->rect.left,y+((cPopUpLineHeight+cPopUpCharMargin)/2)-1);
          glVertex2i(I->Block->rect.right,y+((cPopUpLineHeight+cPopUpCharMargin)/2)-1);
          glEnd();
        }
        y-=cPopUpLineHeight;
        if(!a) y--;
      }
    glColor3fv(block->TextColor);
    BlockOutline(block);
  }
}
  

