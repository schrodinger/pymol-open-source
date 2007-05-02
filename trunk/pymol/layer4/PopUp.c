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

#include "OOMac.h"

#include "main.h"
#include "Base.h"
#include "Pop.h"
#include "PopUp.h"
#include "Ortho.h"
#include "Util.h"
#include "P.h"
#include "Util.h"
#include "Color.h"
#include "Text.h"
#include "PyMOL.h"

#define cPopUpLineHeight 17
#define cPopUpTitleHeight 19
#define cPopUpBarHeight 4
#define cPopUpCharWidth 8
#define cPopUpCharMargin 2
#define cPopUpCharLift 2

#define cChildDelay 0.25
#define cPassiveDelay 0.35
#define cDirtyDelay 0.05

typedef struct CPopUp {
  Block *Block;
  Block *Parent;
  Block *Child;
  int ChildLine;
  int LastX,LastY;
  int StartX,StartY;
  int Selected;
  int Width,Height;
  int NLine;
  PyObject **Sub;
  char **Command;
  char **Text;
  int *Code;
  double ChildDelay;
  double DirtyDelay;
  double PassiveDelay;
  int DirtyDelayFlag;
  int NeverDragged;
  int PlacementAffinity;
}  CPopUp;

int PopUpRelease(Block *block,int button,int x,int y,int mod);
void PopUpDraw(Block *block);
int PopUpDrag(Block *block,int x,int y,int mod);
int PopUpConvertY(CPopUp *I,int value,int mode);

/*========================================================================*/
static Block *PopUpRecursiveFind(Block *block,int x, int y)
{
  PyMOLGlobals *G=block->G;
  CPopUp *I = (CPopUp*)block->reference;
  if(I->Child) { /* favor the child */
    if(PopUpRecursiveFind(I->Child,x,y)==I->Child)
      return block;
  } 
  if(BlockRecursiveFind(block,x,y) == block) {
    OrthoGrab(G,block);
    return block;
  }
  return NULL;
}
/*========================================================================*/
Block *PopUpNew(PyMOLGlobals *G,int x,int y,int last_x,int last_y,
                int passive, PyObject *list,Block *parent)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* assumes blocked threads (calls the Python C API) */

  int mx,a,l,cl,cmx;
  int dim[2];
  PyObject *elem;
  char *str,*c;
  int blocked = PAutoBlock(G);
  OOAlloc(G,CPopUp);

  I->Block = OrthoNewBlock(G,NULL);
  I->Block->reference = (void*)I;
  I->Block->fDraw    = PopUpDraw;
  I->Block->fDrag    = PopUpDrag;
  I->Block->fRelease = PopUpRelease;
  I->Block->active = false;
  I->Block->TextColor[0]=1.0F;
  I->Block->TextColor[1]=1.0F;
  I->Block->TextColor[2]=1.0F;

  I->Block->BackColor[0]=0.1F;
  I->Block->BackColor[1]=0.1F;
  I->Block->BackColor[2]=0.1F;

  I->Parent = parent;
  I->Child = NULL;
  I->NLine=PyList_Size(list);
  I->Text = NULL;
  I->Command = NULL;
  I->Code = NULL;
  I->Selected = -1;
  I->StartX = (I->LastX = last_x);
  I->StartY = (I->LastY = last_y);
  I->ChildDelay = UtilGetSeconds(G) + cChildDelay*2.5;
  I->PassiveDelay = UtilGetSeconds(G) + cPassiveDelay;
  I->DirtyDelay = false;
  I->DirtyDelayFlag = false;
  I->NeverDragged = true;
  I->PlacementAffinity = 0;
  mx=1;
  cmx=1;
  for(a=0;a<I->NLine;a++) {
    elem = PyList_GetItem(PyList_GetItem(list,a),1);
    l = PyString_Size(elem);
    str  = PyString_AsString(elem);
    cl=l;
    c = str;
    while(*c) {
      if(*c=='\\') { /* discount the markup */
        if((((c!=str)&&(c[-1])!='\\'))||
           ((c==str)&&(c[1])&&(c[1]!='\\')))
          cl-=4;
      }
      c++;
    }
    if(cl>cmx) cmx=cl;
    if(l>mx) mx=l;
  }
  I->Width = (cmx * cPopUpCharWidth) + 2 * cPopUpCharMargin;
  
  dim[0]=I->NLine+1;
  dim[1]=mx+1;
  I->Text = (char**)UtilArrayMalloc((unsigned int*)dim,2,1);

  mx=1;
  for(a=0;a<I->NLine;a++) {
    PyObject *command = (PyList_GetItem(PyList_GetItem(list,a),2));
    if(PyString_Check(command)) {
      l = PyString_Size(command);
      if(l>mx) mx=l;
    }
  }
  dim[0]=I->NLine+1;
  dim[1]=mx+1;
  I->Command = (char**)UtilArrayMalloc((unsigned int*)dim,2,1);

  I->Code = Alloc(int,I->NLine+1);
  I->Sub = Calloc(PyObject*,I->NLine+1);
  
  for(a=0;a<I->NLine;a++) {
    PyObject *command;
    elem = PyList_GetItem(list,a);
    I->Code[a]=PyInt_AsLong(PyList_GetItem(elem,0));
    strcpy(I->Text[a],PyString_AsString(PyList_GetItem(elem,1)));
    command=PyList_GetItem(elem,2);
    if(PyString_Check(command)) {
      strcpy(I->Command[a],PyString_AsString(command));
    } else if(PyList_Check(command)) {
      Py_INCREF(command);
      I->Sub[a] = command;
    }
  }

  /* compute height */

  I->Height = 1 * cPopUpCharMargin + PopUpConvertY(I,I->NLine,true);


  I->Block->rect.top=y;
  I->Block->rect.bottom=y-I->Height;
  I->Block->rect.left=x-(I->Width)/3;
  I->Block->rect.right=x+(2*I->Width)/3;

  PopFitBlock(I->Block);

  OrthoAttach(G,I->Block,cOrthoTool);
  I->Block->active=true;
  OrthoGrab(G,I->Block);
  OrthoDirty(G);

  if(passive) 
    PyMOL_SetPassive(G->PyMOL,true);

  PAutoUnblock(G,blocked);

  return I->Block;
#endif

}
/*========================================================================*/
int PopUpConvertY(CPopUp *I,int value,int mode)
{
  int result;
  int a;
  int flag;
  if(mode) { 
    result=0;
      
    /* line to height */
    for(a=0;a<I->NLine;a++)
      {
        if(a>=value)
          break;
        switch(I->Code[a]) {
        case 0:
          result+=cPopUpBarHeight;
          break;
        case 1:
          result+=cPopUpLineHeight;
          break;
        case 2:
          result+=cPopUpTitleHeight;
          break;
        }
      }
  } else {
    flag=false;
    /* height to line */
    if(value<0) {
      result = -1;
    } else {
      result = 0;
      for(a=0;a<I->NLine;a++)
        {
          switch(I->Code[a]) {
          case 0:
            if(value<cPopUpBarHeight)
              flag=true;
            value-=cPopUpBarHeight;
            break;
          case 1:
            if(value<cPopUpLineHeight)
              flag=true;
            value-=cPopUpLineHeight;
            break;
          case 2:
            if(value<cPopUpLineHeight)
              flag=true;
            value-=cPopUpTitleHeight;
            break;
          }
          if(flag) break;
          result++;
        }
      if(!flag)
        result = -1;
      else
        if(result&&!I->Code[result])
          result--;
    }
    /* height to line */
  }
  return(result);
}
/*========================================================================*/

static void PopUpDetachRecursiveChild(Block *block)
{
  PyMOLGlobals *G=block->G;
  CPopUp *I = (CPopUp*)block->reference;

  OrthoDetach(G,block);
  if(I->Child)
    PopUpDetachRecursiveChild(I->Child);
}
/*========================================================================*/

static void PopUpForgetChild(Block *block)
{

  CPopUp *I = (CPopUp*)block->reference;
  I->Child = NULL;
}

static void PopUpRecursiveDetach(Block *block)
{
  PyMOLGlobals *G=block->G;
  CPopUp *I = (CPopUp*)block->reference;
  OrthoDetach(G,block);
  if(I->Child)
    PopUpDetachRecursiveChild(I->Child);
  if(I->Parent) {
    PopUpForgetChild(I->Parent);
    PopUpRecursiveDetach(I->Parent); 
  }
}


static void PopUpFree(Block *block)
{
  PyMOLGlobals *G=block->G;
  CPopUp *I = (CPopUp*)block->reference;

#ifndef _PYMOL_NOPY
  { /* purge code */
    int a;
    int blocked = PAutoBlock(G);
    for(a=0;a<I->NLine;a++)
      if(I->Sub[a]) {
        Py_DECREF(I->Sub[a]);
      }
    PAutoUnblock(G,blocked);
}
#endif

  OrthoDetach(G,I->Block);
  OrthoFreeBlock(G,I->Block);
  FreeP(I->Sub);
  FreeP(I->Code);
  FreeP(I->Command);
  FreeP(I->Text);
  OOFreeP(I);

}

static void PopUpRecursiveFree(Block *block)
{

  CPopUp *I = (CPopUp*)block->reference;

  if(I->Child)
    PopUpFree(I->Child);
  I->Child=NULL;
  if(I->Parent) {
    PopUpForgetChild(I->Parent);
    PopUpRecursiveFree(I->Parent);
  }
  PopUpFree(block);
}

static void PopUpFreeRecursiveChild(Block *block)
{
  CPopUp *I = (CPopUp*)block->reference;
  if(I->Child)
    PopUpFreeRecursiveChild(I->Child);
  I->Child=NULL;
  PopUpFree(block);
}


/*========================================================================*/
int PopUpRelease(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  CPopUp *I = (CPopUp*)block->reference;
  int gone_passive = false;

  if(I->NeverDragged) {
    if(I->PassiveDelay>UtilGetSeconds(G)) {    
      gone_passive = true;
      I->PassiveDelay = UtilGetSeconds(G); /* kill any further delay */
    }
  } 
  if(!gone_passive) {
    if(!I->NeverDragged)
      PopUpDrag(block,x,y,mod);

    /* go passive if we click and release on a sub-menu */

    if((I->Selected>=0)&&(I->Sub[I->Selected])) {
      if((x>=I->Block->rect.left)&&(x<=I->Block->rect.right)) {
        gone_passive = true;
      }
    }
  }  
  if(gone_passive) {
    PyMOL_SetPassive(G->PyMOL,true);
  } else {
    OrthoUngrab(G);
    PopUpRecursiveDetach(block);
    if(!I->NeverDragged) 
      if((I->Selected>=0)&&(!I->Sub[I->Selected])) {
        PLog(G,I->Command[I->Selected],cPLog_pym);
        PParse(G,I->Command[I->Selected]);
        PFlush(G);
      }
    PopUpRecursiveFree(block);
  }
  OrthoDirty(G);
  return(1);
}


/*========================================================================*/
int PopUpDrag(Block *block,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  CPopUp *I = (CPopUp*)block->reference;
  
  int a;
  int was = I->Selected;
  
  if((!I->NeverDragged)&&
     ( ((I->StartX-x)>4) ||
       ((I->StartY-y)>4)))
    I->NeverDragged=false;

  I->LastX=x;
  I->LastY=y;

  x-=I->Block->rect.left;
  y =(I->Block->rect.top -cPopUpCharMargin) -y -1;

  if((x<-2)||(x>(I->Width+2))) {
    int handled_flag = false;
    if(I->Child) {
      if(PopUpRecursiveFind(I->Child,I->LastX,I->LastY) == I->Child) {
        I->Selected = I->ChildLine;
        handled_flag = true;
      } 
    }
    if(!handled_flag) {
      if(I->Parent) { /* are we back in the parent window? */
        I->Selected=-1;
        return PopUpDrag(I->Parent,I->LastX,I->LastY,mod);
      } else if(!I->Child) {
        I->Selected=-1;
      }
    } 
  } else {
    OrthoGrab(G,block);
    a = PopUpConvertY(I,y,false);
    if(I->NLine&&(a==I->NLine))
      if((y-a*cPopUpLineHeight)<4)
        a=I->NLine-1;
    if((a<0)||(a>=I->NLine))
      I->Selected=-1;
    else {
      if(I->Code[a]==1) {
        if((I->Child)&&(I->ChildLine!=a)) {
          if(I->ChildDelay<UtilGetSeconds(G)) {
            PopUpDetachRecursiveChild(I->Child);
            PopUpFreeRecursiveChild(I->Child);
            I->Child=NULL;
            I->ChildLine=-1;
            OrthoDirty(G);
          } else {
            I->Selected = a;
          }
          PyMOL_NeedFakeDrag(G->PyMOL);
        }
      }
      
      if(I->Code[a]!=1)
        I->Selected=-1;
      else if(I->Sub[a]) {
        /* activate submenu */
        if(!I->Child) {
          I->ChildLine = a;
          if(I->ChildDelay>UtilGetSeconds(G)) {
            PyMOL_NeedFakeDrag(G->PyMOL); /* keep coming back here... */
          } else {
            I->Child = PopUpNew(G,I->LastX-300,I->LastY,I->LastX,I->LastY,
                                false,I->Sub[a],I->Block);
            {
              int target_y = block->rect.top - (PopUpConvertY(I,a,true)+ cPopUpCharMargin);
              CPopUp *child = (CPopUp*)(I->Child->reference);
              if(child->NLine)
                if(child->Code[0]!=1)
                  target_y+=cPopUpTitleHeight+2;
              child->PlacementAffinity = PopPlaceChild(I->Child,block->rect.left-5,
                                                       block->rect.right+5,target_y,I->PlacementAffinity);
            }
                          
            OrthoGrab(G,I->Block);
            I->ChildDelay = UtilGetSeconds(G) + cChildDelay; /* leave child up for a while */
          }
          PyMOL_NeedFakeDrag(G->PyMOL); /* keep coming back here... */
        } else if(I->ChildLine==a) { /* on correct line */
          I->ChildDelay = UtilGetSeconds(G) + cChildDelay; /* keep child here for a while */
        }
        I->Selected=a;
      } else 
        I->Selected=a;
    }
  }
  /* delay updates, etc. so that child menus 
     can be comfortably accessed with sloppy mousing */

  if((I->Child)&&(I->Selected!=I->ChildLine))
    PyMOL_NeedFakeDrag(G->PyMOL);

  if(was!=I->Selected) {

    I->NeverDragged = false;
    if(!I->Child) {
      /* we moved, so renew the child delay */
      I->ChildDelay = UtilGetSeconds(G) + cChildDelay;
      PyMOL_NeedFakeDrag(G->PyMOL);
    }

    if((I->Child)&&(I->Selected!=I->ChildLine)) {
      I->DirtyDelayFlag = true;
      I->DirtyDelay = UtilGetSeconds(G) + cDirtyDelay;
    }
    if(!I->DirtyDelayFlag)
      OrthoDirty(G);
  }
  if(I->DirtyDelayFlag && (I->DirtyDelay<UtilGetSeconds(G))) {
    I->DirtyDelayFlag = false;
    OrthoDirty(G);
  }
  return(1);
}
/*========================================================================*/
void PopUpDraw(Block *block)
{
  CPopUp *I = (CPopUp*)block->reference;
  PyMOLGlobals *G=block->G;
  int x,y,a,xx;
  char *c;

  
  if(G->HaveGUI && G->ValidContext ) {

    if((I->Child)&&(I->Selected!=I->ChildLine))
      PyMOL_NeedFakeDrag(G->PyMOL);
    
    /* put raised border around pop-up menu */
    
    /* bottom */
    
    glColor3f(0.2F,0.2F,0.4F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.left-2,block->rect.bottom-2);
    glVertex2i(block->rect.right+2,block->rect.bottom-2);
    glVertex2i(block->rect.right+2,block->rect.bottom+1);
    glVertex2i(block->rect.left-2,block->rect.bottom+1);
    glEnd();

    glColor3f(0.4F,0.4F,0.6F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.left-1,block->rect.bottom-1);
    glVertex2i(block->rect.right+1,block->rect.bottom-1);
    glVertex2i(block->rect.right+1,block->rect.bottom+1);
    glVertex2i(block->rect.left-1,block->rect.bottom+1);
    glEnd();

    /* right */

    glColor3f(0.2F,0.2F,0.4F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.right,block->rect.bottom-2);
    glVertex2i(block->rect.right+2,block->rect.bottom-2);
    glVertex2i(block->rect.right+2,block->rect.top);
    glVertex2i(block->rect.right,block->rect.top);
    glEnd();

    glColor3f(0.4F,0.4F,0.6F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.right,block->rect.bottom-1);
    glVertex2i(block->rect.right+1,block->rect.bottom-1);
    glVertex2i(block->rect.right+1,block->rect.top);
    glVertex2i(block->rect.right,block->rect.top);
    glEnd();

    /* top */

    glColor3f(0.5F,0.5F,0.7F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.left-2,block->rect.top+2);
    glVertex2i(block->rect.right+2,block->rect.top+2);
    glVertex2i(block->rect.right+2,block->rect.top);
    glVertex2i(block->rect.left-2,block->rect.top);
    glEnd();

    glColor3f(0.6F,0.6F,0.8F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.left-1,block->rect.top+1);
    glVertex2i(block->rect.right+1,block->rect.top+1);
    glVertex2i(block->rect.right+1,block->rect.top);
    glVertex2i(block->rect.left-1,block->rect.top);
    glEnd();

    /* left */

    glColor3f(0.5F,0.5F,0.7F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.left-2,block->rect.bottom-2);
    glVertex2i(block->rect.left,block->rect.bottom);
    glVertex2i(block->rect.left,block->rect.top);
    glVertex2i(block->rect.left-2,block->rect.top);
    glEnd();

    glColor3f(0.6F,0.6F,0.8F);
    glBegin(GL_POLYGON);
    glVertex2i(block->rect.left-1,block->rect.bottom-1);
    glVertex2i(block->rect.left,block->rect.bottom-1);
    glVertex2i(block->rect.left,block->rect.top);
    glVertex2i(block->rect.left-1,block->rect.top);
    glEnd();





    glColor3fv(block->BackColor);
    BlockFill(block);
    glColor3fv(block->TextColor);

    if(I->Selected>=0) {

      x = I->Block->rect.left;
      y = I->Block->rect.top-PopUpConvertY(I,I->Selected,true)-cPopUpCharMargin;
      
      y+=2;
      glBegin(GL_POLYGON);
      glVertex2i(x,y);
      glVertex2i(x+I->Width-1,y);
      glVertex2i(x+I->Width-1,y-(cPopUpLineHeight+3));
      glVertex2i(x,y-(cPopUpLineHeight+3));
      glEnd();
    }

    if(I->Code[0]==2) { /* menu name */
        
      glColor3f(0.3F,0.3F,0.6F);
      x = I->Block->rect.left;
      y = I->Block->rect.top;
      
      glBegin(GL_POLYGON);
      glVertex2i(x,y);
      glVertex2i(x+I->Width,y);
      glVertex2i(x+I->Width,y-(cPopUpLineHeight+cPopUpCharMargin));
      glVertex2i(x,y-(cPopUpLineHeight+cPopUpCharMargin));
      glEnd();

      glColor3f(0.2F,0.2F,0.4F);
      glBegin(GL_LINES);
      glVertex2i(x+I->Width-1,y-(cPopUpLineHeight+cPopUpCharMargin));
      glVertex2i(x,y-(cPopUpLineHeight+cPopUpCharMargin));
      glEnd();
      
    }

    x = I->Block->rect.left+cPopUpCharMargin;
    y = (I->Block->rect.top-cPopUpLineHeight)-cPopUpCharMargin+2;
  
    for(a=0;a<I->NLine;a++)
      {
        if(a==I->Selected)
          TextSetColor(G,I->Block->BackColor);
        else
          TextSetColor(G,I->Block->TextColor);          
        if(I->Code[a]) {
          c=I->Text[a];
          xx=x;
          while(*c) {
            if(*c=='\\') if(*(c+1)) if(*(c+2)) if(*(c+3)) {
              if(*(c+1)=='-') {
                if(a==I->Selected)
                  TextSetColor(G,I->Block->BackColor);
                else
                  TextSetColor(G,I->Block->TextColor);          
                c+=4;
              } else if(*(c+1)=='+') {
                c+=4;                
                TextSetColor(G,ColorGetNamed(G,c));
              }else {
                TextSetColor3f(G,(*(c+1)-'0')/9.0F,(*(c+2)-'0')/9.0F,(*(c+3)-'0')/9.0F);
                c+=4;
              }
            }

            TextSetPos2i(G,xx,y+cPopUpCharLift);
            TextDrawChar(G,*(c++));
            xx = xx + 8;
          }


          if(I->Sub[a]) {

            
            glBegin(GL_POLYGON);
            glColor3f(0.4F,0.4F,0.4F);
            glVertex2i(I->Block->rect.left-2,y+1);
            glColor3f(0.1F,0.1F,0.1F);
            glVertex2i(I->Block->rect.left,y+1);
            glVertex2i(I->Block->rect.left,y+((cPopUpLineHeight))-4);
            glColor3f(0.4F,0.4F,0.4F);
            glVertex2i(I->Block->rect.left-2,y+((cPopUpLineHeight))-4);
            glEnd();
            
            glBegin(GL_POLYGON);
            glColor3f(0.1F,0.1F,0.1F);
            glVertex2i(I->Block->rect.right,y+1);
            glColor3f(0.4F,0.4F,0.4F);
            glVertex2i(I->Block->rect.right+2,y+1);
            glVertex2i(I->Block->rect.right+2,y+((cPopUpLineHeight))-4);
            glColor3f(0.1F,0.1F,0.1F);
            glVertex2i(I->Block->rect.right,y+((cPopUpLineHeight))-4);
            glEnd();

          }
        
          y-=cPopUpLineHeight;
          if(I->Code[a]==2)
            y-=2;
        } else {
          glBegin(GL_LINES);
          glColor3f(0.3F,0.3F,0.5F);
          glVertex2i(I->Block->rect.left,y+((cPopUpLineHeight+cPopUpCharMargin)/2)+3);
          glVertex2i(I->Block->rect.right,y+((cPopUpLineHeight+cPopUpCharMargin)/2)+3);
          glColor3f(0.6F,0.6F,0.8F);
          glVertex2i(I->Block->rect.left,y+((cPopUpLineHeight+cPopUpCharMargin)/2)+4);
          glVertex2i(I->Block->rect.right,y+((cPopUpLineHeight+cPopUpCharMargin)/2)+4);
          glEnd();
          y-=cPopUpBarHeight;
        }

      }
    glColor3fv(block->TextColor);
    /*    BlockOutline(block);*/
  }
}
  

