
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
#include"os_python.h"

#ifndef _PYMOL_NOPY

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
#include"Scene.h"

#include"Executive.h"
#include"Block.h"
#include"Text.h"
#include"CGO.h"
#include"vla.h"
#include"pymol/type_traits.h"

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
#define cWizEventView   256
#define cWizEventPosition 512

using cWizEvent_t = int;

struct WizardLine {
  int type;
  WordType text;
  OrthoLineType code;
};

struct CWizard : public Block {
  std::vector<unique_PyObject_ptr_auto_gil> Wiz {};
  pymol::vla<WizardLine> Line {};
  ov_size NLine { 0 };
  int Pressed { -1 };
  int EventMask { 0 };
  int Dirty {};
  int LastUpdatedState { -1 };
  int LastUpdatedFrame { -1 };
  float LastUpdatedPosition[3] {};
  SceneViewType LastUpdatedView {};

  CWizard(PyMOLGlobals * G) : Block(G) {};

  int click(int button, int x, int y, int mod) override;
  int drag(int x, int y, int mod) override;
  void draw(CGO* orthoCGO) override;
  int release(int button, int x, int y, int mod) override;

  bool isEventType(cWizEvent_t eventType) const noexcept;
};

#define cWizardLeftMargin DIP2PIXEL(3)
#define cWizardTopMargin 0
#define cWizardClickOffset DIP2PIXEL(2)

void WizardDirty(PyMOLGlobals * G)
{
  CWizard *I = G->Wizard;
  I->Dirty = true;
  OrthoDirty(G);
}

int WizardUpdate(PyMOLGlobals * G)
{
  CWizard *I = G->Wizard;
  int result = false;

  if(OrthoGetDirty(G)) {
    WizardDoDirty(G);
  }

  {
    int frame = SettingGetGlobal_i(G, cSetting_frame);
    if(frame != I->LastUpdatedFrame) {
      I->LastUpdatedFrame = frame;
      WizardDoFrame(G);
    }
  }
  {
    int state = SettingGetGlobal_i(G, cSetting_state);
    if(state != I->LastUpdatedState) {
      I->LastUpdatedState = state;
      WizardDoState(G);
    }
  }
  WizardDoPosition(G, false);
  WizardDoView(G, false);
  if(I->Dirty) {
    WizardRefresh(G);
    I->Dirty = false;
    result = true;
  }
  return result;
}

void WizardPurgeStack(PyMOLGlobals * G)
{
  pymol::pautoblock block(G);
  G->Wizard->Wiz.clear();
}

bool CWizard::isEventType(cWizEvent_t eventType) const noexcept
{
  return EventMask & eventType;
}

/**
 * Calls wizard's method in python
 *
 * @param wiz wizard
 * @param funcName wizard method name
 * @param func C/C++ function pointer to call the method (see layer1/P.cpp)
 * @param fargs function arguments for func
 * @return result of func
 */

template<typename Func, typename... FuncArgs>
auto WizardCallPython(PyMOLGlobals* G, PyObject* wiz, const char* funcName, Func&& func, FuncArgs&&... fargs)
  -> pymol::result_of_t<Func(PyObject*, const char*, FuncArgs...)>
{
  using result_t = pymol::result_of_t<Func(PyObject*, const char*, FuncArgs...)>;
  result_t result{};
  assert(wiz != nullptr);
  if (PyObject_HasAttrString(wiz, funcName)) {
    result = func(wiz, funcName, std::forward<FuncArgs>(fargs)...);
    PErrPrintIfOccurred(G);
  }
  return result;
}

/**
 * Run when user selects something with the mouse, in a wizard
 */
int WizardDoSelect(PyMOLGlobals* G, const char* name, int state)
{
  CWizard *I = G->Wizard;

  /* if the event is a selection and we're listening for selections */
  if (!I->isEventType(cWizEventSelect)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  /* log if necessary */
  auto buf = pymol::string_format("cmd.get_wizard().do_select('''%s''')", name);
  PLog(G, buf, cPLog_pym);
  /* block and call (in Python) the wizard's do_select */
  pymol::pblock block(G);
  WizardCallPython(G, wiz, "do_pick_state", PTruthCallStr1i, state + 1);
  return WizardCallPython(G, wiz, "do_select", PTruthCallStr, name);
}


/*========================================================================*/
void WizardRefresh(PyMOLGlobals * G)
{

  CWizard *I = G->Wizard;
  char *vla = NULL;
  pymol::pautoblock block(G);

  /* get the current prompt */
  auto wiz = WizardGet(G);
  if (wiz) {
    auto P_list = WizardCallPython(G, wiz, "get_prompt", PyObject_CallMethod, "");
    if(P_list)
      PConvPyListToStringVLA(P_list, &vla);
    Py_XDECREF(P_list);
  }

  OrthoSetWizardPrompt(G, vla);

  /* get the current panel list */

  I->NLine = 0;
  if (wiz) {

    I->EventMask = cWizEventPick + cWizEventSelect;

    auto i = WizardCallPython(G, wiz, "get_event_mask", PyObject_CallMethod, "");
    if (i) {
      if(!PConvPyIntToInt(i, &I->EventMask))
        I->EventMask = cWizEventPick + cWizEventSelect;
      Py_XDECREF(i);
    }

    auto P_list = WizardCallPython(G, wiz, "get_panel", PyObject_CallMethod, "");
    if(P_list) {
      if(PyList_Check(P_list)) {
        std::size_t ll = PyList_Size(P_list);
        I->Line.check(ll);
        for(std::size_t a = 0u; a < ll; a++) {
          /* fallback defaults */

          I->Line[a].text[0] = 0;
          I->Line[a].code[0] = 0;
          I->Line[a].type = 0;

          i = PyList_GetItem(P_list, a);
          if(PyList_Check(i))
            if(PyList_Size(i) > 2) {
              PConvPyObjectToInt(PyList_GetItem(i, 0), &I->Line[a].type);
              PConvPyObjectToStrMaxLen(PyList_GetItem(i, 1),
                                       I->Line[a].text, sizeof(WordType) - 1);
              PConvPyObjectToStrMaxLen(PyList_GetItem(i, 2),
                                       I->Line[a].code, sizeof(OrthoLineType) - 1);
            }
        }
        I->NLine = ll;
      }
      Py_XDECREF(P_list);
    }
  }
  if(I->NLine) {
    int LineHeight = DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
    OrthoReshapeWizard(G, LineHeight * I->NLine + 4);
  } else {
    OrthoReshapeWizard(G, 0);
  }
}


/*========================================================================*/
pymol::Result<> WizardSet(PyMOLGlobals * G, PyObject * wiz, bool replace)
{
  CWizard *I = G->Wizard;
  pymol::pautoblock block(G);
  if ((!wiz) || (wiz == Py_None) || (!I->Wiz.empty() && replace)) {
    /* remove wizard from stack first */
    if (!I->Wiz.empty()) {
      auto old_wiz = std::move(I->Wiz.back());
      I->Wiz.pop_back();
      if (old_wiz) {
        /* then call cleanup, etc. */
        auto result = WizardCallPython(G, old_wiz.get(), "cleanup", PyObject_CallMethod, "");
        PXDecRef(result);
      }
    }
  }
  if (wiz && (wiz != Py_None)) {       /* push */
    I->Wiz.emplace_back(PIncRef(wiz));
  }
  WizardRefresh(G);
  return {};
}


/*========================================================================*/
int WizardActive(PyMOLGlobals * G)
{
  return !G->Wizard->Wiz.empty();
}


/*========================================================================*/
Block *WizardGetBlock(PyMOLGlobals * G)
{
  return G->Wizard;
}


/*========================================================================*/
int WizardDoPick(PyMOLGlobals * G, int bondFlag, int state)
{
  /**
   * Run when user picks something
   */
  CWizard *I = G->Wizard;

  /* process the pick if it happened and we're listening for it */
  if (!I->isEventType(cWizEventPick)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  if(bondFlag)
    PLog(G, "cmd.get_wizard().do_pick(1)", cPLog_pym);
  else
    PLog(G, "cmd.get_wizard().do_pick(0)", cPLog_pym);

  pymol::pblock block(G);
  WizardCallPython(G, wiz, "do_pick_state", PTruthCallStr1i, state + 1);
  return WizardCallPython(G, wiz, "do_pick", PTruthCallStr1i, bondFlag);
}

int WizardDoKey(PyMOLGlobals * G, unsigned char k, int x, int y, int mod)
{
  CWizard *I = G->Wizard;
  if (!I->isEventType(cWizEventKey)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  auto buffer = pymol::string_format("cmd.get_wizard().do_key(%d,%d,%d,%d)", k, x, y, mod);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_key", PTruthCallStr4i, k, x, y, mod);
}

int WizardDoPosition(PyMOLGlobals * G, int force)
{
  CWizard *I = G->Wizard;
  int result = false;
  if (!I->isEventType(cWizEventPosition)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int changed = force;
  if(!changed) {
    float pos[3];
    SceneGetCenter(G, pos);
    changed = ((fabs(pos[0] - I->LastUpdatedPosition[0]) > R_SMALL4) ||
               (fabs(pos[1] - I->LastUpdatedPosition[1]) > R_SMALL4) ||
               (fabs(pos[2] - I->LastUpdatedPosition[2]) > R_SMALL4));
  }
  if(changed) {
    SceneGetCenter(G, I->LastUpdatedPosition);
    pymol::pblock block(G);
    result = WizardCallPython(G, wiz, "do_position", PTruthCallStr0);
  }
  return result;
}

int WizardDoView(PyMOLGlobals * G, int force)
{
  CWizard *I = G->Wizard;
  if (!I->isEventType(cWizEventView)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int result = false;
  int changed = force;
  if(!changed) {
    SceneViewType view;
    SceneGetView(G, view);
    changed = !SceneViewEqual(view, I->LastUpdatedView);
  }
  if(changed) {
    SceneGetView(G, I->LastUpdatedView);
    pymol::pblock block(G);
    result = WizardCallPython(G, wiz, "do_view", PTruthCallStr0);
  }
  return result;
}

int WizardDoScene(PyMOLGlobals * G)
{
  CWizard *I = G->Wizard;
  if (!I->isEventType(cWizEventScene)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  std::string buffer = "cmd.get_wizard().do_scene()";
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_scene", PTruthCallStr0);
}

int WizardDoDirty(PyMOLGlobals * G)
{
  CWizard *I = G->Wizard;
  if (!I->isEventType(cWizEventDirty)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  std::string buffer = "cmd.get_wizard().do_dirty()";
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_dirty", PTruthCallStr0);
}

int WizardDoState(PyMOLGlobals * G)
{
  CWizard *I = G->Wizard;
  if (!I->isEventType(cWizEventState)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int state = SettingGet<int>(G, cSetting_state);
  auto buffer = pymol::string_format("cmd.get_wizard().do_state(%d)", state);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_state", PTruthCallStr1i, state);
}

int WizardDoFrame(PyMOLGlobals * G)
{
  CWizard *I = G->Wizard;
  if (!I->isEventType(cWizEventFrame)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int frame = SettingGet<int>(G, cSetting_frame) + 1;
  auto buffer = pymol::string_format("cmd.get_wizard().do_frame(%d)", frame);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_frame", PTruthCallStr1i, frame);
}

int WizardDoSpecial(PyMOLGlobals * G, int k, int x, int y, int mod)
{
  CWizard *I = G->Wizard;

  if (!I->isEventType(cWizEventSpecial)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  auto buffer = pymol::string_format("cmd.get_wizard().do_special(%d,%d,%d,%d)", k, x, y, mod);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_special", PTruthCallStr4i, k, x, y, mod);
}


/*========================================================================*/
int CWizard::click(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CWizard *I = G->Wizard;
  int LineHeight = DIP2PIXEL(SettingGet<int>(G, cSetting_internal_gui_control_size));

  int a = ((rect.top - (y + cWizardClickOffset)) - cWizardTopMargin) / LineHeight;
  if((a >= 0) && ((ov_size) a < I->NLine)) {
    switch (I->Line[a].type) {
    case cWizTypeButton:
      OrthoGrab(G, this);
      I->Pressed = (int) a;
      OrthoDirty(G);
      break;
    case cWizTypePopUp:
      {
        pymol::pblock block(G);
        auto wiz = WizardGet(G);
        PyObject* menuList = NULL;
        if (wiz) {
          menuList = WizardCallPython(G, wiz, "get_menu", PyObject_CallMethod, "s", I->Line[a].code);
        }

        if(menuList && (menuList != Py_None)) {
          int my = rect.top - (cWizardTopMargin + a * LineHeight) - 2;

          PopUpNew(G, x, my, x, y, false, menuList, NULL);
        }
        Py_XDECREF(menuList);
      }
      break;
    }
  }
  return (1);
}


/*========================================================================*/
int CWizard::drag(int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;

  CWizard *I = G->Wizard;
  int LineHeight = DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));

  int a;
  a = ((rect.top - (y + cWizardClickOffset)) - cWizardTopMargin) / LineHeight;

  if((x < rect.left) || (x > rect.right))
    a = -1;

  if(I->Pressed != a) {
    I->Pressed = -1;
    OrthoDirty(G);
  }
  if((a >= 0) && ((ov_size) a < I->NLine)) {

    switch (I->Line[a].type) {
    case cWizTypeButton:
      if(I->Pressed != a) {
        I->Pressed = a;
        OrthoDirty(G);
      }
      break;
    }
  }
  return (1);
}


/*========================================================================*/
int CWizard::release(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;

  CWizard *I = this; // TODO: Remove during Wizard Refactor
  int LineHeight = DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));

  int a;
  a = ((rect.top - (y + cWizardClickOffset)) - cWizardTopMargin) / LineHeight;

  if(I->Pressed)
    I->Pressed = -1;
  OrthoDirty(G);

  OrthoUngrab(G);

  if((a >= 0) && ((ov_size) a < I->NLine)) {
    switch (I->Line[a].type) {
    case cWizTypeButton:
      {
        if (auto wiz = WizardGet(G)) {
          PLog(G, I->Line[a].code, cPLog_pym);
          PParse(G, I->Line[a].code);
          PFlush(G);
        }
      }
      break;
    }
  }
  I->Pressed = -1;
  return (1);
}

static void draw_button(int x2, int y2, int w, int h, float *light, float *dark,
                        float *inside ORTHOCGOARG)
{
    if (orthoCGO){
      CGOColorv(orthoCGO, light);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, x2, y2, 0.f);
      CGOVertex(orthoCGO, x2, y2 + h, 0.f);
      CGOVertex(orthoCGO, x2 + w, y2, 0.f);
      CGOVertex(orthoCGO, x2 + w, y2 + h, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3fv(light);
      glBegin(GL_POLYGON);
      glVertex2i(x2, y2);
      glVertex2i(x2, y2 + h);
      glVertex2i(x2 + w, y2 + h);
      glVertex2i(x2 + w, y2);
      glEnd();
    }

    if (orthoCGO){
      CGOColorv(orthoCGO, dark);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, x2 + 1, y2, 0.f);
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
      CGOVertex(orthoCGO, x2 + w, y2, 0.f);
      CGOVertex(orthoCGO, x2 + w, y2 + h - 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3fv(dark);
      glBegin(GL_POLYGON);
      glVertex2i(x2 + 1, y2);
      glVertex2i(x2 + 1, y2 + h - 1);
      glVertex2i(x2 + w, y2 + h - 1);
      glVertex2i(x2 + w, y2);
      glEnd();
    }

  if(inside) {
    if (orthoCGO){
      CGOColorv(orthoCGO, inside);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, 0.f);
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, 0.f);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3fv(inside);
      glBegin(GL_POLYGON);
      glVertex2i(x2 + 1, y2 + 1);
      glVertex2i(x2 + 1, y2 + h - 1);
      glVertex2i(x2 + w - 1, y2 + h - 1);
      glVertex2i(x2 + w - 1, y2 + 1);
      glEnd();
    }
  } else {                      /* rainbow */
    if (orthoCGO){
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOColor(orthoCGO, 0.1F, 1.0F, 0.1F); // green
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
      CGOColor(orthoCGO, 1.0F, 1.0F, 0.1F);  // yellow
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, 0.f);
      CGOColor(orthoCGO, 1.f, 0.1f, 0.1f); // red
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, 0.f);
      CGOColor(orthoCGO, 0.1F, 0.1F, 1.0F);  // blue
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glBegin(GL_POLYGON);
      glColor3f(1.0F, 0.1F, 0.1F);
      glVertex2i(x2 + 1, y2 + 1);
      glColor3f(0.1F, 1.0F, 0.1F);
      glVertex2i(x2 + 1, y2 + h - 1);
      glColor3f(1.0F, 1.0F, 0.1F);
      glVertex2i(x2 + w - 1, y2 + h - 1);
      glColor3f(0.1F, 0.1F, 1.0F);
      glVertex2i(x2 + w - 1, y2 + 1);
      glEnd();
    }
  }

}

static void draw_text(PyMOLGlobals * G, char *c, int xx, int yy, float *color ORTHOCGOARG)
{
  TextSetColor(G, color);
  while(*c) {
    if(TextSetColorFromCode(G, c, color)) {
      c += 4;
    }
    TextSetPos2i(G, xx, yy);
    TextDrawChar(G, *(c++) ORTHOCGOARGVAR);
    xx = xx + DIP2PIXEL(8);
  }
}


/*========================================================================*/
void CWizard::draw(CGO* orthoCGO)
{
  PyMOLGlobals *G = m_G;

  CWizard *I = G->Wizard;
  int x, y;
  int a;

  float buttonTextColor[3] = { 1.0, 1.0, 1.0 };
  float buttonActiveColor[3] = { 0.8F, 0.8F, 0.8F };

  float dimColor[3] = { 0.45F, 0.45F, 0.45F };

  float dimLightEdge[3] = { 0.6F, 0.6F, 0.6F };
  float dimDarkEdge[3] = { 0.25F, 0.25F, 0.25F };

  float menuBGColor[3] = { 0.5F, 0.5F, 1.0 };
  float menuLightEdge[3] = { 0.7F, 0.7F, 0.9F };
  float menuDarkEdge[3] = { 0.3F, 0.3F, 0.5F };

  float black_color[3] = { 0.0F, 0.0F, 0.0F };
  float menuColor[3] = { 0.0, 0.0, 0.0 };
  int LineHeight = DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
  int text_lift = (LineHeight / 2) - DIP2PIXEL(5);
  float *text_color, *text_color2 = TextColor;

  text_color = menuColor;

  if(G->HaveGUI && G->ValidContext && ((rect.right - rect.left) > 6)) {

    if(SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) == InternalGUIMode::Default) {
    if (orthoCGO)
      CGOColorv(orthoCGO, BackColor);
#ifndef PURE_OPENGL_ES_2
    else
      glColor3fv(BackColor);
#endif
      fill(orthoCGO);
      drawLeftEdge(orthoCGO);
    } else {
      drawLeftEdge(orthoCGO);
    if (orthoCGO)
      CGOColor(orthoCGO, .5f, .5f, .5f);
    else
      glColor3f(0.5, 0.5, 0.5);
      drawTopEdge();
      text_color2 = OrthoGetOverlayColor(G);
    }

    if (orthoCGO)
      CGOColorv(orthoCGO, TextColor);
#ifndef PURE_OPENGL_ES_2
    else
      glColor3fv(TextColor);
#endif
    x = rect.left + cWizardLeftMargin;
    y = (rect.top - LineHeight) - cWizardTopMargin;

    for(a = 0; (ov_size) a < I->NLine; a++) {
      if(I->Pressed == a) {
        draw_button(rect.left + 1, y,
                    (rect.right - rect.left) - 1,
                    LineHeight - 1, dimLightEdge, dimDarkEdge, buttonActiveColor ORTHOCGOARGVAR);
        /*        glColor3f(0.0,0.0,0.0); */
        text_color = black_color;
      } else {
        switch (I->Line[a].type) {
        case cWizTypeText:
          text_color = text_color2;
          glColor3fv(text_color2);
          break;
        case cWizTypeButton:
          draw_button(rect.left + 1, y,
                      (rect.right - rect.left) - 1,
                      LineHeight - 1, dimLightEdge, dimDarkEdge, dimColor ORTHOCGOARGVAR);

          /*          glColor3fv(buttonTextColor); */
          text_color = buttonTextColor;
          break;
        case cWizTypePopUp:
          draw_button(rect.left + 1, y,
                      (rect.right - rect.left) - 1,
                      LineHeight - 1, menuLightEdge, menuDarkEdge, menuBGColor ORTHOCGOARGVAR);
          /* glColor3fv(menuColor); */
          text_color = menuColor;
          break;
        default:
          break;
        }
      }
      draw_text(G, I->Line[a].text, x, y + text_lift, text_color ORTHOCGOARGVAR);
      /*GrapDrawStr(I->Line[a].text,x+1,y+text_lift); */
      y -= LineHeight;
    }
  }
}


/*========================================================================*/
PyObject* WizardGet(PyMOLGlobals* G)
{
  auto I = G->Wizard;
  return I->Wiz.empty() ? nullptr : I->Wiz.back().get();
}

/*========================================================================*/
PyObject* WizardGetStack(PyMOLGlobals * G)
{
  CWizard *I = G->Wizard;
  auto result = PyList_New(I->Wiz.size());
  for (std::size_t i = 0; i < I->Wiz.size(); i++) {
    auto wiz = I->Wiz[i].get();
    Py_INCREF(wiz);
    PyList_SetItem(result, i, wiz);     /* steals ref */
  }
  return result;
}


/*========================================================================*/
pymol::Result<> WizardSetStack(PyMOLGlobals * G, PyObject * list)
{
  CWizard *I = G->Wizard;

  if (list == nullptr || !PyList_Check(list)) {
    return pymol::make_error("Invalid list.");
  }
  WizardPurgeStack(G);
  auto size = PyList_Size(list);
  pymol::pautoblock block(G);
  for (int a = 0; a < size; a++) {
    auto wiz = PyList_GetItem(list, a);
    I->Wiz.emplace_back(PIncRef(wiz));
  }
  WizardRefresh(G);
  OrthoDirty(G);
  return {};
}


/*========================================================================*/
int WizardInit(PyMOLGlobals * G)
{
  auto I = new CWizard(G);
  G->Wizard = I;
  I->active = true;
  I->TextColor[0] = 0.2F;
  I->TextColor[1] = 1.0F;
  I->TextColor[2] = 0.2F;

  OrthoAttach(G, I, cOrthoTool);

  I->Line.resize(1);
  return 1;
}


/*========================================================================*/
void WizardFree(PyMOLGlobals * G)
{
  WizardPurgeStack(G);
  DeleteP(G->Wizard);
}

std::vector<unique_PyObject_ptr_auto_gil> WizardGetWizardCopies(PyMOLGlobals* G)
{
  auto I = G->Wizard;
  std::vector<unique_PyObject_ptr_auto_gil> copies;
  copies.reserve(I->Wiz.size());
  pymol::pautoblock block(G);
  for (std::size_t i = 0u; i < I->Wiz.size(); i++) {
    copies.emplace_back(PIncRef(I->Wiz[i].get()));
  }
  return copies;
}

void WizardSetWizards(PyMOLGlobals* G, const std::vector<unique_PyObject_ptr_auto_gil>& wizs)
{
  auto I = G->Wizard;
  WizardPurgeStack(G);
  I->Wiz.reserve(wizs.size());
  pymol::pautoblock block(G);
  for (std::size_t i = 0u; i < wizs.size(); i++) {
    I->Wiz.emplace_back(PIncRef(wizs[i].get()));
  }
  WizardRefresh(G);
  WizardDirty(G);
  OrthoDirty(G);
}

#endif
