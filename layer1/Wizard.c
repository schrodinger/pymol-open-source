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

#define cWizardLeftMargin 2
#define cWizardTopMargin 0
#define cWizardClickOffset 2

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
    int LineHeight = SettingGetGlobal_i(cSetting_internal_gui_control_size);
    OrthoReshapeWizard(LineHeight*I->NLine+4);
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
  int LineHeight = SettingGetGlobal_i(cSetting_internal_gui_control_size);

  a=((I->Block->rect.top-(y+cWizardClickOffset))-cWizardTopMargin)/LineHeight;
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
        int my = I->Block->rect.top-(cWizardTopMargin + a*LineHeight) - 2;
        
        PopUpNew(x,my,x,y,menuList,NULL);
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
  int LineHeight = SettingGetGlobal_i(cSetting_internal_gui_control_size);

  int a;
  a=((I->Block->rect.top-(y+cWizardClickOffset))-
     cWizardTopMargin)/LineHeight;

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
    int LineHeight = SettingGetGlobal_i(cSetting_internal_gui_control_size);

  int a;
  a=((I->Block->rect.top-(y+cWizardClickOffset))-
     cWizardTopMargin)/LineHeight;

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

static void draw_button(int x2,int y2, int w, int h, float *light, float *dark, float *inside)
{
  glColor3fv(light);
  glBegin(GL_POLYGON);
  glVertex2i(x2,y2);
  glVertex2i(x2,y2+h);
  glVertex2i(x2+w,y2+h);
  glVertex2i(x2+w,y2);
  glEnd();
  
  glColor3fv(dark);
  glBegin(GL_POLYGON);
  glVertex2i(x2+1,y2);
  glVertex2i(x2+1,y2+h-1);
  glVertex2i(x2+w,y2+h-1);
  glVertex2i(x2+w,y2);
  glEnd();
  
  if(inside) {
    glColor3fv(inside);
    glBegin(GL_POLYGON);
    glVertex2i(x2+1,y2+1);
    glVertex2i(x2+1,y2+h-1);
    glVertex2i(x2+w-1,y2+h-1);
    glVertex2i(x2+w-1,y2+1);
    glEnd();
  } else { /* rainbow */
    glBegin(GL_POLYGON);
    glColor3f(1.0F,0.1F,0.1F);
    glVertex2i(x2+1,y2+1);
    glColor3f(0.1F,1.0F,0.1F);
    glVertex2i(x2+1,y2+h-1);
    glColor3f(1.0F,1.0F,0.1F);
    glVertex2i(x2+w-1,y2+h-1);
    glColor3f(0.1F,0.1F,1.0F);
    glVertex2i(x2+w-1,y2+1);
    glEnd();
  }

}

static void draw_text(char *c,int xx,int yy,float *color)
{
  glColor3fv(color);
  while(*c) {
    if(*c=='\\') if(*(c+1)) if(*(c+2)) if(*(c+3)) {
      if(*(c+1)=='-') {
        glColor3fv(color);
        c+=4;
      } else {
        glColor3f((*(c+1)-'0')/9.0F,(*(c+2)-'0')/9.0F,(*(c+3)-'0')/9.0F);
        c+=4;
      }
    }
    glRasterPos4d((double)(xx),(double)(yy),0.0,1.0);
    p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c++));
    xx = xx + 8;
  }
}

/*========================================================================*/
static void WizardDraw(Block *block)
{
  CWizard *I=&Wizard;
  int x,y;
  int a;

  float buttonTextColor[3] = { 1.0,1.0,1.0 };
  float buttonActiveColor[3] = { 0.8F,0.8F,0.8F };
  
  float dimColor[3] = {0.45F, 0.45F,0.45F};

  float dimLightEdge[3] = {0.6F, 0.6F,0.6F};
  float dimDarkEdge[3] = {0.25F, 0.25F,0.25F};

  float menuBGColor[3] = {0.5F, 0.5F,1.0};
  float menuLightEdge[3] = {0.7,0.7,0.9F};
  float menuDarkEdge[3] = {0.3,0.3,0.5F};

  float black_color[3] = {0.0F,0.0F,0.0F};
  float menuColor[3] = { 0.0,0.0,0.0};
  int LineHeight = SettingGetGlobal_i(cSetting_internal_gui_control_size);
  int text_lift = (LineHeight/2)-5;
  float *text_color ;

  text_color = menuColor;

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    
    glColor3fv(I->Block->TextColor);

    x = I->Block->rect.left+cWizardLeftMargin;
    y = (I->Block->rect.top-LineHeight)-cWizardTopMargin;

    for(a=0;a<I->NLine;a++) {
      if(I->Pressed==a) {
          draw_button(I->Block->rect.left+1,y,
                      (I->Block->rect.right-I->Block->rect.left)-1,
                      LineHeight-1,
                      dimLightEdge,
                      dimDarkEdge,
                      buttonActiveColor);
          /*        glColor3f(0.0,0.0,0.0);*/
          text_color = black_color;
      } else {
        switch(I->Line[a].type) {
        case cWizTypeText:

          glColor3fv(I->Block->TextColor);
          text_color = I->Block->TextColor;
          break;
        case cWizTypeButton:
          draw_button(I->Block->rect.left+1,y,
                      (I->Block->rect.right-I->Block->rect.left)-1,
                      LineHeight-1,
                      dimLightEdge,
                      dimDarkEdge,
                      dimColor);

          /*          glColor3fv(buttonTextColor);*/
          text_color = buttonTextColor;
         break;
        case cWizTypePopUp:
          draw_button(I->Block->rect.left+1,y,
                      (I->Block->rect.right-I->Block->rect.left)-1,
                      LineHeight-1,
                      menuLightEdge,
                      menuDarkEdge,
                      menuBGColor);
          /* glColor3fv(menuColor);*/
          text_color = menuColor;
          break;
        default:
          break;
        }
      }
      draw_text(I->Line[a].text,x+1,y+text_lift,text_color);
      /*GrapDrawStr(I->Line[a].text,x+1,y+text_lift);*/
      y-=LineHeight;
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

  I->Block->TextColor[0]=0.2F;
  I->Block->TextColor[1]=1.0F;
  I->Block->TextColor[2]=0.2F;

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

