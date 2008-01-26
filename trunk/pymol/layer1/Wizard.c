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
#include"PopUp.h"

#include"Wizard.h"

#include"Executive.h"
#include"Block.h"
#include"Text.h"

#define cWizBlank      0
#define cWizTypeText   1
#define cWizTypeButton 2
#define cWizTypePopUp  3

#define cWizEventPick    1
#define cWizEventSelect  2
#define cWizEventKey     4
#define cWizEventSpecial 8
#define cWizEventScene   16
#define cWizEventState   32
#define cWizEventFrame   64
#define cWizEventDirty  128

typedef struct {
  int type;
  WordType text;
  OrthoLineType code;
} WizardLine;

struct _CWizard {
  Block *Block;
  PyObject **Wiz;
  WizardLine *Line;
  ov_size NLine;
  ov_diff Stack;
  int Pressed;
  int EventMask;
  int Dirty;
  int LastUpdatedState;
  int LastUpdatedFrame;
  
};

#define cWizardLeftMargin 2
#define cWizardTopMargin 0
#define cWizardClickOffset 2

void WizardDirty(PyMOLGlobals *G)
{
  register CWizard *I=G->Wizard;
  I->Dirty=true;
  OrthoDirty(G);
}

int WizardUpdate(PyMOLGlobals *G)
{
  register CWizard *I=G->Wizard;
  int result = false;

  if(OrthoGetDirty(G)) {
      WizardDoDirty(G);    
  }

  {
    int frame = SettingGetGlobal_i(G,cSetting_frame);
    int state = SettingGetGlobal_i(G,cSetting_state);
    if(frame!=I->LastUpdatedFrame) {
      I->LastUpdatedFrame = frame;
      WizardDoFrame(G);
    }
    if(state!=I->LastUpdatedState) {
      I->LastUpdatedState = state;
      WizardDoState(G);
    }
  }

  if(I->Dirty) {
    WizardRefresh(G);
    I->Dirty=false;
    result = true;
  }

  return result;
}

void WizardPurgeStack(PyMOLGlobals *G)
{
#ifndef _PYMOL_NOPY
  int blocked;
  ov_diff a;
  register CWizard *I=G->Wizard;
  blocked = PAutoBlock(G);
  for(a=I->Stack;a>=0;a--)
    Py_XDECREF(I->Wiz[a]);
  I->Stack = -1;
  PAutoUnblock(G,blocked);
#endif
}
int WizardDoSelect(PyMOLGlobals *G,char *name)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  OrthoLineType buf;
  register CWizard *I=G->Wizard;
  int result = false;

  if(I->EventMask & cWizEventSelect) 
    if(I->Stack>=0)
      if(I->Wiz[I->Stack]) {
        sprintf(buf,"cmd.get_wizard().do_select('''%s''')",name);
        PLog(G,buf,cPLog_pym);
        PBlock(G); 
        if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_select")) {
          result = PTruthCallStr(I->Wiz[I->Stack],"do_select",name);
        if(PyErr_Occurred()) PyErr_Print();
        }
        PUnblock(G);
      }
  return result;
#endif
}
/*========================================================================*/
void WizardRefresh(PyMOLGlobals *G)
{

#ifndef _PYMOL_NOPY
  register CWizard *I = G->Wizard;
  char *vla = NULL;
  PyObject *P_list;
  ov_size ll;
  PyObject *i;
  ov_size a;
  int blocked;
  blocked = PAutoBlock(G);
  
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

  OrthoSetWizardPrompt(G,vla);

  /* get the current panel list */
  
  I->NLine = 0;
  if(I->Stack>=0)
    if(I->Wiz[I->Stack]) {

      I->EventMask = cWizEventPick + cWizEventSelect;
      
      if(PyObject_HasAttrString(I->Wiz[I->Stack],"get_event_mask")) {      
        i = PyObject_CallMethod(I->Wiz[I->Stack],"get_event_mask","");
        if(PyErr_Occurred()) PyErr_Print();
        if(!PConvPyIntToInt(i,&I->EventMask))
          I->EventMask = cWizEventPick + cWizEventSelect;
        Py_XDECREF(i);
      }

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
    int LineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
    OrthoReshapeWizard(G,LineHeight*I->NLine+4);
  } else {
    OrthoReshapeWizard(G,0);
  }
  PAutoUnblock(G,blocked);
#endif
}
/*========================================================================*/
void WizardSet(PyMOLGlobals *G,PyObject *wiz,int replace)
{
#ifndef _PYMOL_NOPY
  register CWizard *I = G->Wizard;
  int blocked;
  blocked = PAutoBlock(G);
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
  WizardRefresh(G);
  PAutoUnblock(G,blocked);
#endif
}
/*========================================================================*/
int WizardActive(PyMOLGlobals *G)
{
  register CWizard *I = G->Wizard;
  if(!I->Wiz)
    return(false);
  if(I->Stack<0)
    return false;
  return(I->Wiz[I->Stack]&&1);
}
/*========================================================================*/
Block *WizardGetBlock(PyMOLGlobals *G)
{
  register CWizard *I=G->Wizard;
  {return(I->Block);}
}
/*========================================================================*/
int WizardDoPick(PyMOLGlobals *G,int bondFlag)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CWizard *I=G->Wizard;
  int result=false;
  if(I->EventMask & cWizEventPick) 
    if(I->Stack>=0) 
      if(I->Wiz[I->Stack]) {
        if(bondFlag)
          PLog(G,"cmd.get_wizard().do_pick(1)",cPLog_pym);
        else
          PLog(G,"cmd.get_wizard().do_pick(0)",cPLog_pym);
        
        PBlock(G); 
        if(I->Stack>=0)
          if(I->Wiz[I->Stack]) {
            if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_pick")) {
              result = PTruthCallStr1i(I->Wiz[I->Stack],"do_pick",bondFlag);
              if(PyErr_Occurred()) PyErr_Print();
            }
          }
        PUnblock(G);
      }
  return result;
#endif
}

int WizardDoKey(PyMOLGlobals *G,unsigned char k, int x, int y, int mod)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CWizard *I=G->Wizard;
  int result=false;
  if(I->EventMask & cWizEventKey) 
    if(I->Stack>=0) 
      if(I->Wiz[I->Stack]) {
        OrthoLineType buffer;
        sprintf(buffer,"cmd.get_wizard().do_key(%d,%d,%d,%d)",k,x,y,mod);
        PLog(G,buffer,cPLog_pym);
        PBlock(G); 
        if(I->Stack>=0)
          if(I->Wiz[I->Stack]) {
            if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_key")) {
              result = PTruthCallStr4i(I->Wiz[I->Stack],"do_key",k,x,y,mod);
              if(PyErr_Occurred()) PyErr_Print();
            }
          }
        PUnblock(G);
      }
  return result;
#endif
}

int WizardDoScene(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CWizard *I=G->Wizard;
  int result=false;
  if(I->EventMask & cWizEventScene) 
    if(I->Stack>=0) 
      if(I->Wiz[I->Stack]) {
        OrthoLineType buffer;
        sprintf(buffer,"cmd.get_wizard().do_scene()");
        PLog(G,buffer,cPLog_pym);
        PBlock(G); 
        if(I->Stack>=0)
          if(I->Wiz[I->Stack]) {
            if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_scene")) {
              result = PTruthCallStr0(I->Wiz[I->Stack],"do_scene");
              if(PyErr_Occurred()) PyErr_Print();
            }
          }
        PUnblock(G);
      }
  return result;
#endif
}

int WizardDoDirty(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CWizard *I=G->Wizard;
  int result=false;
  if(I->EventMask & cWizEventDirty) 
    if(I->Stack>=0) 
      if(I->Wiz[I->Stack]) {
        OrthoLineType buffer;
        sprintf(buffer,"cmd.get_wizard().do_dirty()");
        PLog(G,buffer,cPLog_pym);
        PBlock(G); 
        if(I->Stack>=0)
          if(I->Wiz[I->Stack]) {
            if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_dirty")) {
              result = PTruthCallStr0(I->Wiz[I->Stack],"do_dirty");
              if(PyErr_Occurred()) PyErr_Print();
            }
          }
        PUnblock(G);
      }
  return result;
#endif
}


int WizardDoState(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CWizard *I=G->Wizard;
  int result=false;
  if(I->EventMask & cWizEventState) 
    if(I->Stack>=0) 
      if(I->Wiz[I->Stack]) {
        OrthoLineType buffer;
        int state = SettingGetGlobal_i(G,cSetting_state) + 1;
        sprintf(buffer,"cmd.get_wizard().do_state(%d)",state);
        PLog(G,buffer,cPLog_pym);
        PBlock(G); 
        if(I->Stack>=0)
          if(I->Wiz[I->Stack]) {
            if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_state")) {
              result = PTruthCallStr1i(I->Wiz[I->Stack],"do_state",state);
              if(PyErr_Occurred()) PyErr_Print();
            }
          }
        PUnblock(G);
      }
  return result;
#endif
}


int WizardDoFrame(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CWizard *I=G->Wizard;
  int result=false;
  if(I->EventMask & cWizEventFrame) 
    if(I->Stack>=0) 
      if(I->Wiz[I->Stack]) {
        OrthoLineType buffer;
        int frame = SettingGetGlobal_i(G,cSetting_frame) + 1;
        sprintf(buffer,"cmd.get_wizard().do_frame(%d)",frame);
        PLog(G,buffer,cPLog_pym);
        PBlock(G); 
        if(I->Stack>=0)
          if(I->Wiz[I->Stack]) {
            if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_frame")) {
              result = PTruthCallStr1i(I->Wiz[I->Stack],"do_frame",frame);
              if(PyErr_Occurred()) PyErr_Print();
            }
          }
        PUnblock(G);
      }
  return result;
#endif
}

int WizardDoSpecial(PyMOLGlobals *G,int k, int x, int y, int mod)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CWizard *I=G->Wizard;
  int result=false;

  if(I->EventMask & cWizEventSpecial) 
    if(I->Stack>=0) 
      if(I->Wiz[I->Stack]) {
        OrthoLineType buffer;
        sprintf(buffer,"cmd.get_wizard().do_special(%d,%d,%d,%d)",k,x,y,mod);
        PLog(G,buffer,cPLog_pym);
        PBlock(G); 
        if(I->Stack>=0)
          if(I->Wiz[I->Stack]) {
            if(PyObject_HasAttrString(I->Wiz[I->Stack],"do_special")) {
              result = PTruthCallStr4i(I->Wiz[I->Stack],"do_special",k,x,y,mod);
              if(PyErr_Occurred()) PyErr_Print();
            }
          }
        PUnblock(G);
      }
  return result;
#endif
}

/*========================================================================*/
static int WizardClick(Block *block,int button,int x,int y,int mod)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  PyMOLGlobals *G=block->G;
  register CWizard *I=G->Wizard;
  int a;
  PyObject *menuList=NULL;
  int LineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);

  a=((I->Block->rect.top-(y+cWizardClickOffset))-cWizardTopMargin)/LineHeight;
  if((a>=0)&&((ov_size)a<I->NLine)) {
    switch(I->Line[a].type) {
    case cWizTypeButton:
      OrthoGrab(G,I->Block);
      I->Pressed=(int)a;
      OrthoDirty(G);
      break;
    case cWizTypePopUp:
      PBlock(G); 
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
        
        PopUpNew(G,x,my,x,y,false,menuList,NULL);
      }
      Py_XDECREF(menuList);
      PUnblock(G);
      break;
    }
  }
  return(1);
#endif
}
/*========================================================================*/
static int WizardDrag(Block *block,int x,int y,int mod)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  PyMOLGlobals *G=block->G;

  register CWizard *I=G->Wizard;
  int LineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);

  int a;
  a=((I->Block->rect.top-(y+cWizardClickOffset))-
     cWizardTopMargin)/LineHeight;

  if((x<I->Block->rect.left)||(x>I->Block->rect.right))
    a=-1;

  if(I->Pressed!=a) {
    I->Pressed=-1;
    OrthoDirty(G);
  }
  if((a>=0)&&((ov_size)a<I->NLine)) {

    switch(I->Line[a].type) {
    case cWizTypeButton:
      if(I->Pressed!=a) {
        I->Pressed=a;
        OrthoDirty(G);
      }
      break;
    }
  }
  return(1);
#endif
}
/*========================================================================*/
static int WizardRelease(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;

  register CWizard *I=G->Wizard;
    int LineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);

  int a;
  a=((I->Block->rect.top-(y+cWizardClickOffset))-
     cWizardTopMargin)/LineHeight;

  if(I->Pressed)
    I->Pressed=-1;
  OrthoDirty(G);

  OrthoUngrab(G);

  if((a>=0)&&((ov_size)a<I->NLine)) {
    switch(I->Line[a].type) {
    case cWizTypeButton:
      if(I->Stack>=0)
        if(I->Wiz[I->Stack]) {
          PLog(G,I->Line[a].code,cPLog_pym);
          PParse(G,I->Line[a].code);
          PFlush(G);
        }
      break;
    }
  }
  I->Pressed = -1;
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

static void draw_text(PyMOLGlobals *G,char *c,int xx,int yy,float *color)
{
  TextSetColor(G,color);
  while(*c) {
    if(*c=='\\') if(*(c+1)) if(*(c+2)) if(*(c+3)) {
      if(*(c+1)=='-') {
        TextSetColor(G,color);
        c+=4;
      } else {
        TextSetColor3f(G,(*(c+1)-'0')/9.0F,(*(c+2)-'0')/9.0F,(*(c+3)-'0')/9.0F);
        c+=4;
      }
    }
    TextSetPos2i(G,xx,yy);
    TextDrawChar(G,*(c++));
    xx = xx + 8;
  }
}

/*========================================================================*/
static void WizardDraw(Block *block)
{
  PyMOLGlobals *G=block->G;

  register CWizard *I=G->Wizard;
  int x,y;
  int a;

  float buttonTextColor[3] = { 1.0,1.0,1.0 };
  float buttonActiveColor[3] = { 0.8F,0.8F,0.8F };
  
  float dimColor[3] = {0.45F, 0.45F,0.45F};

  float dimLightEdge[3] = {0.6F, 0.6F,0.6F};
  float dimDarkEdge[3] = {0.25F, 0.25F,0.25F};

  float menuBGColor[3] = {0.5F, 0.5F,1.0};
  float menuLightEdge[3] = {0.7F,0.7F,0.9F};
  float menuDarkEdge[3] = {0.3F,0.3F,0.5F};

  float black_color[3] = {0.0F,0.0F,0.0F};
  float menuColor[3] = { 0.0,0.0,0.0};
  int LineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
  int text_lift = (LineHeight/2)-5;
  float *text_color ;

  text_color = menuColor;

  if(G->HaveGUI && G->ValidContext && ((block->rect.right-block->rect.left)>6)) {

    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    
    glColor3fv(I->Block->TextColor);
    
    
    x = I->Block->rect.left+cWizardLeftMargin;
    y = (I->Block->rect.top-LineHeight)-cWizardTopMargin;

    for(a=0;(ov_size)a<I->NLine;a++) {
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
      draw_text(G,I->Line[a].text,x+1,y+text_lift,text_color);
      /*GrapDrawStr(I->Line[a].text,x+1,y+text_lift);*/
      y-=LineHeight;
    }
  }
}
/*========================================================================*/
PyObject *WizardGet(PyMOLGlobals *G)
{
  register CWizard *I=G->Wizard;
  if(!I->Wiz)
    return(NULL);
  if(I->Stack<0)
    return(NULL);
  return(I->Wiz[I->Stack]);
}
/*========================================================================*/
PyObject *WizardGetStack(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  register CWizard *I=G->Wizard;
  PyObject *result;

  result = PyList_New(I->Stack-(-1));
  if(I->Wiz) {
    ov_diff a;
    for(a=I->Stack;a>=0;a--) {
      Py_INCREF(I->Wiz[a]);
      PyList_SetItem(result,a,I->Wiz[a]); /* steals ref */
    }
  }
  return(result);
#endif
}
/*========================================================================*/
int WizardSetStack(PyMOLGlobals *G,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  register CWizard *I=G->Wizard;
  int ok= true;
  
  if(I->Wiz) {
    WizardPurgeStack(G);
    if(ok) ok = (list!=NULL);
    if(ok) ok = PyList_Check(list);
    if(ok) {
      I->Stack = PyList_Size(list)-1;
      if(I->Stack>=0) {
        ov_diff a;
		VLACheck(I->Wiz,PyObject*,I->Stack);
        for(a=I->Stack;a>=0;a--) {
          I->Wiz[a] = PyList_GetItem(list,a);
          Py_INCREF(I->Wiz[a]);
        }
      }
    }
    if(ok) WizardRefresh(G);
    if(ok) OrthoDirty(G);
  }
  return(ok);
#endif
}
/*========================================================================*/
int WizardInit(PyMOLGlobals *G)
{
  register CWizard *I=NULL;
  if( (I=(G->Wizard=Calloc(CWizard,1)))) {

  I->Block = OrthoNewBlock(G,NULL);
  I->Block->fClick = WizardClick;
  I->Block->fDrag  = WizardDrag;
  I->Block->fDraw    = WizardDraw;
  I->Block->fReshape = BlockReshape;
  I->Block->fRelease = WizardRelease;
  I->Block->active = true;

  I->Block->TextColor[0]=0.2F;
  I->Block->TextColor[1]=1.0F;
  I->Block->TextColor[2]=0.2F;

  I->LastUpdatedState = -1;
  I->LastUpdatedFrame = -1;

  OrthoAttach(G,I->Block,cOrthoTool);

  I->Line = VLAlloc(WizardLine,1);
  I->NLine = 0;
  I->Pressed = -1;
  I->EventMask = 0;
  I->Stack = -1;
  I->Wiz = VLAlloc(PyObject*,10);
  return 1;
  } else
    return 0;
}
/*========================================================================*/
void WizardFree(PyMOLGlobals *G)
{
  register CWizard *I = G->Wizard;
  WizardPurgeStack(G);
  OrthoFreeBlock(G,I->Block);
  VLAFreeP(I->Line);
  VLAFreeP(I->Wiz);
  FreeP(G->Wizard);
}

