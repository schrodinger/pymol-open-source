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
  PyObject **Wiz;
  WizardLine *Line;
  int NLine;
  int Stack;
  int Pressed;
}  CWizard;

CWizard Wizard;

#define cWizardLineHeight 14
#define cWizardLeftMargin 2
#define cWizardTopMargin (-1)
#define cWizardClickOffset 4

void WizardPurgeStack(void)
{
  int blocked;
  int a;
  CWizard *I=&Wizard;
  blocked = PAutoBlock();
  for(a=I->Stack;a>=0;a--)
    Py_XDECREF(I->Wiz[a]);
  I->Stack = -1;
  PAutoUnblock(blocked);

}
void WizardDoSelect(char *name)
{
  OrthoLineType buf;
  CWizard *I=&Wizard;
  if(I->Stack>=0)
    if(I->Wiz[I->Stack]) {
      sprintf(buf,"cmd.get_wizard().do_select('''%s''')",name);
      PLog(buf,cPLog_pym);
      PBlock(); 
      if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_select")) {
        PXDecRef(PyObject_CallMethod(I->Wiz[I->Stack],"do_select","s",name));
        if(PyErr_Occurred()) PyErr_Print();
      }
    PUnblock();
  }
}
/*========================================================================*/
void WizardRefresh(void)
{
  CWizard *I = &Wizard;
  char *vla = NULL;
  PyObject *P_list;
  int ll;
  PyObject *i;
  int a;
  int blocked;
  blocked = PAutoBlock();

  /* get the current prompt */
  if(I->Stack>=0)
    if(I->Wiz[I->Stack]) {
      vla = NULL;
      if(PyObject_HasAttrString(I->Wiz[I->Stack],"get_prompt")) {
        P_list = PyObject_CallMethod(I->Wiz[I->Stack],"get_prompt","");
        if(PyErr_Occurred()) PyErr_Print();
        if(P_list) 
          PConvPyListToStringVLA(P_list,&vla);
        Py_XDECREF(P_list);
      }
    }

  OrthoSetWizardPrompt(vla);

  /* get the current panel list */
  
  I->NLine = 0;
  if(I->Stack>=0)
    if(I->Wiz[I->Stack]) {
      
      if(PyObject_HasAttrString(I->Wiz[I->Stack],"get_panel")) {
        P_list = PyObject_CallMethod(I->Wiz[I->Stack],"get_panel","");
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
                  PConvPyObjectToStrMaxLen(PyList_GetItem(i,1),
                                           I->Line[a].text,
                                           sizeof(WordType)-1);
                  PConvPyObjectToStrMaxLen(PyList_GetItem(i,2),
                                           I->Line[a].code,
                                           sizeof(OrthoLineType)-1);
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
  PAutoUnblock(blocked);
}
/*========================================================================*/
void WizardSet(PyObject *wiz,int replace)
{
  CWizard *I = &Wizard;
  int blocked;
  blocked = PAutoBlock();
  if(I->Wiz) {
    if((!wiz)||(wiz==Py_None)||((I->Stack>=0)&&replace)) { 
      if(I->Stack>=0) {  /* pop */
        if(I->Wiz[I->Stack]) {
          if(PyObject_HasAttrString(I->Wiz[I->Stack],"cleanup")) {
            PXDecRef(PyObject_CallMethod(I->Wiz[I->Stack],"cleanup",""));
            if(PyErr_Occurred()) PyErr_Print();
          }
          Py_DECREF(I->Wiz[I->Stack]);
          I->Wiz[I->Stack]=NULL;
          I->Stack--;
        }
      }
    }
    if(wiz&&(wiz!=Py_None)) { /* push */
      if(wiz) { 
        I->Stack++;
        VLACheck(I->Wiz,PyObject*,I->Stack);
        I->Wiz[I->Stack]=wiz;
        if(I->Wiz[I->Stack])
          Py_INCREF(I->Wiz[I->Stack]);
      }
    }
  }
  WizardRefresh();
  PAutoUnblock(blocked);
}
/*========================================================================*/
int WizardActive(void)
{
  CWizard *I = &Wizard;
  if(!I->Wiz)
    return(false);
  if(I->Stack<0)
    return false;
  return(I->Wiz[I->Stack]&&1);
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
  if(I->Stack>=0) 
    if(I->Wiz[I->Stack]) {
      if(bondFlag)
        PLog("cmd.get_wizard().do_pick(1)",cPLog_pym);
      else
        PLog("cmd.get_wizard().do_pick(0)",cPLog_pym);
      
      PBlock(); 
      if(I->Stack>=0)
        if(I->Wiz[I->Stack]) {
          if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_pick")) {
            PXDecRef(PyObject_CallMethod(I->Wiz[I->Stack],"do_pick","i",bondFlag));
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
      if(I->Stack>=0)
        if(I->Wiz[I->Stack]) {
          if(PyObject_HasAttrString(I->Wiz[I->Stack],"get_menu")) {
            menuList = PyObject_CallMethod(I->Wiz[I->Stack],"get_menu","s",I->Line[a].code);
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
  a=((I->Block->rect.top-(y+cWizardClickOffset))-
     cWizardTopMargin)/cWizardLineHeight;

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
  a=((I->Block->rect.top-(y+cWizardClickOffset))-
     cWizardTopMargin)/cWizardLineHeight;

  if(I->Pressed)
    I->Pressed=-1;
  OrthoDirty();

  OrthoUngrab();

  if((a>=0)&&(a<I->NLine)) {
    switch(I->Line[a].type) {
    case cWizTypeButton:
      if(I->Stack>=0)
        if(I->Wiz[I->Stack]) {
          PLog(I->Line[a].code,cPLog_pym);
          PParse(I->Line[a].code);
          PFlush();
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

  float buttonTextColor[3] = { 1.0,1.0,1.0 };
  float buttonActiveColor[3] = { 0.8,0.8,0.8 };
  
  float dimColor[3] = {0.3, 0.3,0.3};
  float menuBGColor[3] = {0.5, 0.5,1.0};
  float menuColor[3] = { 0.0,0.0,0.0};

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    
    glColor3fv(I->Block->TextColor);

    x = I->Block->rect.left+cWizardLeftMargin;
    y = (I->Block->rect.top-cWizardLineHeight)-cWizardTopMargin;

    for(a=0;a<I->NLine;a++) {
      if(I->Pressed==a) {
        glColor3fv(buttonActiveColor);
        glBegin(GL_POLYGON);
        glVertex2i(I->Block->rect.left,y+cWizardLineHeight-3);
        glVertex2i(I->Block->rect.right,y+cWizardLineHeight-3);
        glVertex2i(I->Block->rect.right,y-2);
        glVertex2i(I->Block->rect.left,y-2);
        glEnd();
        glColor3f(0.0,0.0,0.0);
      } else {
        switch(I->Line[a].type) {
        case cWizTypeText:

          glColor3fv(I->Block->TextColor);
          break;
        case cWizTypeButton:
          glColor3fv(dimColor);
          glBegin(GL_POLYGON);
          glVertex2i(I->Block->rect.left,y+cWizardLineHeight-3);
          glVertex2i(I->Block->rect.right,y+cWizardLineHeight-3);
          glVertex2i(I->Block->rect.right,y-2);
          glVertex2i(I->Block->rect.left,y-2);
          glEnd();
          glColor3fv(buttonTextColor);

         break;
        case cWizTypePopUp:
          glColor3fv(menuBGColor);
          glBegin(GL_POLYGON);
          glVertex2i(I->Block->rect.left,y+cWizardLineHeight-3);
          glVertex2i(I->Block->rect.right,y+cWizardLineHeight-3);
          glVertex2i(I->Block->rect.right,y-2);
          glVertex2i(I->Block->rect.left,y-2);
          glEnd();
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
  if(!I->Wiz)
    return(NULL);
  if(I->Stack<0)
    return(NULL);
  return(I->Wiz[I->Stack]);
}
/*========================================================================*/
PyObject *WizardGetStack(void)
{
  CWizard *I=&Wizard;
  int a;
  PyObject *result;

  result = PyList_New(I->Stack-(-1));
  if(I->Wiz) {
    for(a=I->Stack;a>=0;a--) {
      Py_INCREF(I->Wiz[a]);
      PyList_SetItem(result,a,I->Wiz[a]); /* steals ref */
    }
  }
  return(result);
}
/*========================================================================*/
int WizardSetStack(PyObject *list)
{
  CWizard *I=&Wizard;
  int a;
  int ok= true;
  
  if(I->Wiz) {
    WizardPurgeStack();
    if(ok) ok = (list!=NULL);
    if(ok) ok = PyList_Check(list);
    if(ok) {
      I->Stack = PyList_Size(list)-1;
      if(I->Stack>=0) {
        VLACheck(I->Wiz,PyObject*,I->Stack);
        for(a=I->Stack;a>=0;a--) {
          I->Wiz[a] = PyList_GetItem(list,a);
          Py_INCREF(I->Wiz[a]);
        }
      }
    }
    if(ok) WizardRefresh();
    if(ok) OrthoDirty();
  }
  return(ok);
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

  I->Stack = -1;
  I->Wiz = VLAlloc(PyObject*,10);
}
/*========================================================================*/
void WizardFree(void)
{
  CWizard *I = &Wizard;
  WizardPurgeStack();
  OrthoFreeBlock(I->Block);
  VLAFreeP(I->Line);
  VLAFreeP(I->Wiz)
    
}

