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
#include"Err.h"
#include"main.h"

#include"MemoryDebug.h"
#include"Ortho.h"
#include"P.h"
#include"PConv.h"
#include"Grap.h"
#include"PopUp.h"

#include"Wizard.h"

#include"Executive.h"
#include"Block.h"

#define cWizBlank      0
#define cWizTypeText   1
#define cWizTypeButton 2
#define cWizTypePopUp  3

typedef struct {
  int type;
  WordType text;
  OrthoLineType code;
} WizardLine;

typedef struct {
  Block *Block;
  PyObject *Wiz;
  WizardLine *Line;
  int NLine;
  int Pressed;
}  CWizard;

CWizard Wizard;

#define cWizardLineHeight 14
#define cWizardLeftMargin 2
#define cWizardTopMargin -1
#define cWizardClickOffset 4

/*========================================================================*/
void WizardRefresh(void)
{
  CWizard *I = &Wizard;
  char *vla;
  PyObject *P_list;
  int ll;
  PyObject *i;
  int a;
  PBlock();

  /* get the current prompt */

  if(I->Wiz) {
    vla = NULL;
    if(PyObject_HasAttrString(I->Wiz,"get_prompt")) {
      P_list = PyObject_CallMethod(I->Wiz,"get_prompt","");
      if(PyErr_Occurred()) PyErr_Print();
      if(P_list) 
        PConvPyListToStringVLA(P_list,&vla);
      Py_XDECREF(P_list);
    }
    OrthoSetWizardPrompt(vla);
  }
  /* get the current panel list */

  I->NLine = 0;
  if(I->Wiz) {

    if(PyObject_HasAttrString(I->Wiz,"get_panel")) {
      P_list = PyObject_CallMethod(I->Wiz,"get_panel","");
      if(PyErr_Occurred()) PyErr_Print();
      if(P_list) {
        if(PyList_Check(P_list)) {
          ll = PyList_Size(P_list);
          VLACheck(I->Line,WizardLine,ll);
          for(a=0;a<ll;a++) {
            /* fallback defaults */
            
            I->Line[a].text[0]=0;
            I->Line[a].code[0]=0;
            I->Line[a].type = 0;
            
            i = PyList_GetItem(P_list,a);
            if(PyList_Check(i))
              if(PyList_Size(i)>2) {
                PConvPyObjectToInt(PyList_GetItem(i,0),&I->Line[a].type);
                PConvPyObjectToStrMaxLen(PyList_GetItem(i,1),I->Line[a].text,sizeof(WordType)-1);
                PConvPyObjectToStrMaxLen(PyList_GetItem(i,2),I->Line[a].code,sizeof(OrthoLineType)-1);
              }
          }
          I->NLine=ll;
        }
      }
      Py_XDECREF(P_list);
    }
  }
  if(I->NLine) {
    OrthoReshapeWizard(cWizardLineHeight*I->NLine+4);
  } else {
    OrthoReshapeWizard(0);
  }
  PUnblock();
}
/*========================================================================*/
void WizardSet(PyObject *wiz)
{
  CWizard *I = &Wizard;
  PBlock();
  if(I->Wiz) {
    if(PyObject_HasAttrString(I->Wiz,"cleanup")) {
      PXDecRef(PyObject_CallMethod(I->Wiz,"cleanup",""));
      if(PyErr_Occurred()) PyErr_Print();
    }
    Py_DECREF(I->Wiz);
    I->Wiz=NULL;
  }
  I->Wiz=wiz;
  if(I->Wiz)
    Py_INCREF(I->Wiz);
  PUnblock();
  WizardRefresh();
}
/*========================================================================*/
int WizardActive(void)
{
  CWizard *I = &Wizard;
  return(I->Wiz&&1);
}
/*========================================================================*/
Block *WizardGetBlock(void)
{
  CWizard *I=&Wizard;
  {return(I->Block);}
}
/*========================================================================*/
void WizardDoPick(int bondFlag)
{
  CWizard *I=&Wizard;
  if(I->Wiz) {
    PBlock(); 
    if(I->Wiz) {
      if(PyObject_HasAttrString(I->Wiz,"do_pick")) {
        PXDecRef(PyObject_CallMethod(I->Wiz,"do_pick","i",bondFlag));
        if(PyErr_Occurred()) PyErr_Print();
      }
    }
    PUnblock();
  }
}
/*========================================================================*/
static int WizardClick(Block *block,int button,int x,int y,int mod)
{
  CWizard *I=&Wizard;

  int a;
  PyObject *menuList=NULL;

  a=((I->Block->rect.top-(y+cWizardClickOffset))-cWizardTopMargin)/cWizardLineHeight;
  if((a>=0)&&(a<I->NLine)) {
    switch(I->Line[a].type) {
    case cWizTypeButton:
      OrthoGrab(I->Block);
      I->Pressed=a;
      OrthoDirty();
      break;
    case cWizTypePopUp:
      PBlock(); 
      if(I->Wiz) {
        if(PyObject_HasAttrString(I->Wiz,"get_menu")) {
          menuList = PyObject_CallMethod(I->Wiz,"get_menu","s",I->Line[a].code);
          if(PyErr_Occurred()) PyErr_Print();
        }
      }
     
      if(PyErr_Occurred()) PyErr_Print();
      if(menuList&&(menuList!=Py_None)) {
        y = I->Block->rect.top-(cWizardTopMargin + a*cWizardLineHeight) -2;
        
        PopUpNew(x,y,menuList);
      }
      Py_XDECREF(menuList);
      PUnblock();

      break;
    }
  }
  return(1);
}
/*========================================================================*/
static int WizardDrag(Block *block,int x,int y,int mod)
{
  CWizard *I=&Wizard;

  int a;
  a=((I->Block->rect.top-(y+cWizardClickOffset))-cWizardTopMargin)/cWizardLineHeight;

  if((x<I->Block->rect.left)||(x>I->Block->rect.right))
    a=-1;

  if(I->Pressed!=a) {
    I->Pressed=-1;
    OrthoDirty();
  }
  if((a>=0)&&(a<I->NLine)) {

    switch(I->Line[a].type) {
    case cWizTypeButton:
      if(I->Pressed!=a) {
        I->Pressed=a;
        OrthoDirty();
      }
      break;
    }
  }
  return(1);
}
/*========================================================================*/
static int WizardRelease(Block *block,int button,int x,int y,int mod)
{
  CWizard *I=&Wizard;

  int a;
  a=((I->Block->rect.top-(y+cWizardClickOffset))-cWizardTopMargin)/cWizardLineHeight;

  if(I->Pressed)
    I->Pressed=-1;
  OrthoDirty();

  OrthoUngrab();

  if((a>=0)&&(a<I->NLine)) {
    switch(I->Line[a].type) {
    case cWizTypeButton:
      if(I->Wiz) {
        PBlockAndUnlockAPI();
        PRunString(I->Line[a].code);
        PLockAPIAndUnblock();
      }
      break;
    }
  }
  return(1);
}
/*========================================================================*/
static void WizardDraw(Block *block)
{
  CWizard *I=&Wizard;
  int x,y;
  int a;

  float buttonColor[3] = { 0.5, 0.5, 1.0 };
  float menuColor[3] = { 1.0, 0.2, 0.2 };

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    
    glColor3fv(I->Block->TextColor);

    x = I->Block->rect.left+cWizardLeftMargin;
    y = (I->Block->rect.top-cWizardLineHeight)-cWizardTopMargin;

    for(a=0;a<I->NLine;a++) {
      if(I->Pressed==a) {
        glColor3fv(buttonColor);
        glBegin(GL_POLYGON);
        glVertex2i(I->Block->rect.left,y+cWizardLineHeight-1);
        glVertex2i(I->Block->rect.right,y+cWizardLineHeight-1);
        glVertex2i(I->Block->rect.right,y-3);
        glVertex2i(I->Block->rect.left,y-3);
        glEnd();
        glColor3f(0.0,0.0,0.0);
      } else {
        switch(I->Line[a].type) {
        case cWizTypeText:
          glColor3fv(I->Block->TextColor);
          break;
        case cWizTypeButton:
          glColor3fv(buttonColor);
          break;
        case cWizTypePopUp:
          glColor3fv(menuColor);
          break;
        default:
          break;
        }
      }
      GrapDrawStr(I->Line[a].text,x,y);
      y-=cWizardLineHeight;
    }
  }
}
/*========================================================================*/
PyObject *WizardGet(void)
{
  CWizard *I=&Wizard;
  return(I->Wiz);
}
/*========================================================================*/
void WizardInit(void)
{
  CWizard *I = &Wizard;

  I->Block = OrthoNewBlock(NULL);
  I->Block->fClick = WizardClick;
  I->Block->fDrag  = WizardDrag;
  I->Block->fDraw    = WizardDraw;
  I->Block->fReshape = BlockReshape;
  I->Block->fRelease = WizardRelease;
  I->Block->active = true;

  I->Block->TextColor[0]=0.2;
  I->Block->TextColor[1]=1.0;
  I->Block->TextColor[2]=0.2;

  OrthoAttach(I->Block,cOrthoTool);

  I->Line = VLAlloc(WizardLine,10);
  I->NLine = 0;
  I->Pressed = -1;

  I->Wiz = NULL;
}
/*========================================================================*/
void WizardFree(void)
{
  CWizard *I = &Wizard;
  OrthoFreeBlock(I->Block);
  VLAFreeP(I->Line);
  Py_XDECREF(I->Wiz);
}

