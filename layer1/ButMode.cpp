
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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"
#include"MemoryDebug.h"
#include "main.h"
#include "Base.h"
#include "ButMode.h"
#include "Scene.h"
#include "Util.h"
#include "Ortho.h"
#include "Setting.h"
#include "P.h"
#include "Text.h"
#include "Menu.h"
#include "CGO.h"
#include "Movie.h"

#define cButModeLineHeight DIP2PIXEL(12)
#define cButModeLeftMargin DIP2PIXEL(2)
#define cButModeTopMargin DIP2PIXEL(1)
#define cButModeBottomMargin DIP2PIXEL(2)

struct CButMode : public Block {
  CodeType Code[cButModeCount + 1] {};
  int NCode {};
  int Mode[cButModeInputCount] {};
  int NBut {};
  float Rate { 0.0f };
  float RateShown { 0.0f };
  float Samples { 0.0f }, Delay { 0.0f };
  float TextColor1[3] = { 0.5f, 0.5f, 1.0f };
  float TextColor2[3] = { 0.8f, 0.8f, 0.8f };
  float TextColor3[3] = { 1.0f, 0.5f, 0.5f };
  int DeferCnt { 0 };
  float DeferTime { 0.0f };

  CButMode(PyMOLGlobals * G) : Block(G) {}

  int click(int button, int x, int y, int mod) override;
  void draw(CGO *orthoCGO) override;
  bool fastDraw(CGO* orthoCGO) override;
};


/*========================================================================*/
Block *ButModeGetBlock(PyMOLGlobals * G)
{
  CButMode *I = G->ButMode;
  {
    return (I);
  }
}

int ButModeGetHeight(PyMOLGlobals * G)
{
  if(SettingGetGlobal_b(G, cSetting_mouse_grid))
    return DIP2PIXEL(124);
  else
    return DIP2PIXEL(40);
}


/*========================================================================*/
int ButModeGet(PyMOLGlobals * G, int button)
{
  CButMode *I = G->ButMode;
  if((button >= 0) && (button < I->NBut)) {
    return I->Mode[button];
  }
  return 0;
}

void ButModeSet(PyMOLGlobals * G, int button, int action)
{
  CButMode *I = G->ButMode;
  if((button >= 0) && (button < I->NBut) && (action >= 0) && (action < I->NCode)) {
    I->Mode[button] = action;
    OrthoDirty(G);
  }
}


/*========================================================================*/
void ButModeSetRate(PyMOLGlobals * G, float interval)
{
  CButMode *I = G->ButMode;

  if(interval >= 0.001F) {      /* sub-millisecond, defer... */
    if(I->DeferCnt) {
      interval = (interval + I->DeferTime) / (I->DeferCnt + 1);
      I->DeferCnt = 0;
      I->DeferTime = 0.0F;
    }
    I->Delay -= interval;
    if(interval < 1.0F) {
      I->Samples *= 0.95 * (1.0F - interval);
      I->Rate *= 0.95 * (1.0F - interval);
    } else {
      I->Samples = 0.0F;
      I->Rate = 0.0F;
    }

    I->Samples++;
    I->Rate += 1.0F / interval;
  } else {
    I->DeferCnt++;
    I->DeferTime += interval;
  }
}


/*========================================================================*/
void ButModeResetRate(PyMOLGlobals * G)
{
  CButMode *I = G->ButMode;
  I->Samples = 0.0;
  I->Rate = 0.0;
  I->RateShown = 0;
  I->Delay = 0;
}


/*========================================================================*/
void ButModeFree(PyMOLGlobals * G)
{
  DeleteP(G->ButMode);
}


/*========================================================================*/
int CButMode::click(int button, int x, int y, int mod)
{
  int dy = (y - rect.bottom) / cButModeLineHeight;
  //  int dx = (x - block->rect.left);
  int forward = 1; 
  // TAKEN OUT : BB 12/11 : Mouse position should not have
  // an effect on whether the mouse ring goes forwards or 
  // backwards.
  // forward = (dx > ((block->rect.right - block->rect.left) / 2));
  // 
  /*  register CButMode *I=block->G->ButMode; */
  if(button == P_GLUT_BUTTON_SCROLL_BACKWARD || button == P_GLUT_RIGHT_BUTTON)
    forward = !forward;
  if(mod == cOrthoSHIFT)
    forward = !forward;
  if(dy < 2) {
    if(ButModeTranslate(m_G, P_GLUT_SINGLE_LEFT, 0) != cButModePickAtom) {
      if(!forward) {
        PLog(m_G, "cmd.mouse('select_backward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse select_backward,quiet=1");
      } else {
        PLog(m_G, "cmd.mouse('select_forward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse select_forward,quiet=1");
      }
    }
  } else {
    if(button == P_GLUT_RIGHT_BUTTON) {
      MenuActivate0Arg(m_G,x,y,x,y,false,"mouse_config");
    } else {
      if(!forward) {
        PLog(m_G, "cmd.mouse('backward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse backward,quiet=1");
      } else {
        PLog(m_G, "cmd.mouse('forward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse forward,quiet=1");
      }
    }
  }
  return (1);
}

static bool ButModeDrawFastImpl(Block * block, short definitely , CGO *orthoCGO);
/*========================================================================*/
void CButMode::draw(CGO* orthoCGO)
{
  int x, y, a;
  int mode;
  const float *textColor = TextColor;
  const float *textColor2 = TextColor2;
  CButMode *I = this; // TODO: Remove during ButMode refactor
#define BLANK_STR "     "

  if(m_G->HaveGUI && m_G->ValidContext && ((rect.right - rect.left) > 6)) {
    if(SettingGet<InternalGUIMode>(m_G, cSetting_internal_gui_mode) == InternalGUIMode::Default) {
      if (orthoCGO)
	CGOColorv(orthoCGO, BackColor);
      else
	glColor3fv(BackColor);
      fill(orthoCGO);
      drawLeftEdge(orthoCGO);
    } else {
      drawLeftEdge(orthoCGO);
      if (orthoCGO)
	CGOColor(orthoCGO, .5f, .5f, .5f);
      else
	glColor3f(0.5, 0.5, 0.5);
      drawTopEdge();
      textColor2 = OrthoGetOverlayColor(m_G);
      textColor = textColor2;
    }

    x = rect.left + cButModeLeftMargin;
    y = (rect.top - cButModeLineHeight) - cButModeTopMargin;

    TextSetColor(m_G, textColor);
    TextDrawStrAt(m_G, "Mouse Mode ", x + 1, y, orthoCGO);
    TextSetColor(m_G, TextColor3);
    TextDrawStrAt(m_G, SettingGetGlobal_s(m_G, cSetting_button_mode_name), x + DIP2PIXEL(88), y, orthoCGO);
    /*    TextDrawStrAt(m_G,"2-Bttn Selecting",x+88,y); */
    y -= cButModeLineHeight;

    if(SettingGetGlobal_b(m_G, cSetting_mouse_grid)) {

      TextSetColor(m_G, TextColor3);
      TextDrawStrAt(m_G, "Buttons", x + DIP2PIXEL(6), y, orthoCGO);
      TextSetColor(m_G, TextColor1);
      /*    TextDrawStrAt(m_G,"  Left Mddl Rght Scrl",x+48,y); */
      TextDrawStrAt(m_G, "    L    M    R  Wheel", x + DIP2PIXEL(43), y, orthoCGO);

      y -= cButModeLineHeight;
      /*    glColor3fv(I->Block->TextColor);
         TextDrawStrAt(m_G,"K",x,y-4); */
      TextSetColor(m_G, TextColor3);
      TextDrawStrAt(m_G, "&", x + DIP2PIXEL(12), y, orthoCGO);
      TextDrawStrAt(m_G, "Keys", x + DIP2PIXEL(24), y, orthoCGO);
      TextSetColor(m_G, textColor2);

      TextSetPos2i(m_G, x + DIP2PIXEL(64), y);
      for(a = 0; a < 3; a++) {
        mode = Mode[a];
        if(mode < 0)
          TextDrawStr(m_G, BLANK_STR, orthoCGO);
        else
          TextDrawStr(m_G, Code[mode], orthoCGO);
      }
      mode = Mode[12];
      if(mode < 0)
        TextDrawStr(m_G, BLANK_STR, orthoCGO);
      else
        TextDrawStr(m_G, Code[mode], orthoCGO);

      y -= cButModeLineHeight;
      /*    TextSetColor(m_G,I->Block->TextColor);
         TextDrawStrAt(m_G,"e",x+5,y-1); */
      TextSetColor(m_G, TextColor1);

      TextSetColor(m_G, TextColor1);
      TextDrawStrAt(m_G, "Shft ", x + DIP2PIXEL(24), y, orthoCGO);
      TextSetColor(m_G, textColor2);
      TextSetPos2i(m_G, x + DIP2PIXEL(64), y);
      for(a = 3; a < 6; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(m_G, BLANK_STR, orthoCGO);
        else
          TextDrawStr(m_G, Code[mode], orthoCGO);
      }
      mode = Mode[13];
      if(mode < 0)
        TextDrawStr(m_G, BLANK_STR, orthoCGO);
      else
        TextDrawStr(m_G, Code[mode], orthoCGO);

      y -= cButModeLineHeight;
      /*    glColor3fv(I->Block->TextColor);
         TextDrawStrAt(m_G,"y",x+10,y+2); */
      TextSetColor(m_G, TextColor1);
      TextDrawStrAt(m_G, "Ctrl ", x + DIP2PIXEL(24), y, orthoCGO);
      TextSetColor(m_G, textColor2);
      TextSetPos2i(m_G, x + DIP2PIXEL(64), y);
      for(a = 6; a < 9; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(m_G, BLANK_STR, orthoCGO);
        else
          TextDrawStr(m_G, Code[mode], orthoCGO);
      }
      mode = I->Mode[14];
      if(mode < 0)
        TextDrawStr(m_G, BLANK_STR, orthoCGO);
      else
        TextDrawStr(m_G, Code[mode], orthoCGO);
      y -= cButModeLineHeight;

      /*    glColor3fv(I->Block->TextColor);
         TextDrawStrAt(m_G,"s",x+15,y+3); */
      TextSetColor(m_G, TextColor1);
      TextSetColor(m_G, TextColor1);
      TextDrawStrAt(m_G, "CtSh ", x + DIP2PIXEL(24), y, orthoCGO);
      TextSetColor(m_G, textColor2);
      TextSetPos2i(m_G, x + DIP2PIXEL(64), y);
      for(a = 9; a < 12; a++) {
        mode = Mode[a];
        if(mode < 0)
          TextDrawStr(m_G, BLANK_STR, orthoCGO);
        else
          TextDrawStr(m_G, Code[mode], orthoCGO);
      }
      mode = Mode[15];
      if(mode < 0)
        TextDrawStr(m_G, BLANK_STR, orthoCGO);
      else
        TextDrawStr(m_G, Code[mode], orthoCGO);

      y -= cButModeLineHeight;

      TextSetColor(m_G, TextColor);
      TextSetColor(m_G, TextColor1);
      TextDrawStrAt(m_G, " SnglClk", x - DIP2PIXEL(8), y, orthoCGO);
      TextSetColor(m_G, textColor2);
      TextSetPos2i(m_G, x + DIP2PIXEL(64), y);
      for(a = 19; a < 22; a++) {
        mode = Mode[a];
        if(mode < 0)
          TextDrawStr(m_G, BLANK_STR, orthoCGO);
        else
          TextDrawStr(m_G, Code[mode], orthoCGO);
      }
      TextSetColor(m_G, TextColor);
      y -= cButModeLineHeight;

      TextSetColor(m_G, TextColor);
      TextSetColor(m_G, TextColor1);
      TextDrawStrAt(m_G, " DblClk", x, y, orthoCGO);
      TextSetColor(m_G, textColor2);
      TextSetPos2i(m_G, x + DIP2PIXEL(64), y);
      for(a = 16; a < 19; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(m_G, BLANK_STR, orthoCGO);
        else
          TextDrawStr(m_G, Code[mode], orthoCGO);
      }
      TextSetColor(m_G, TextColor);
      y -= cButModeLineHeight;

    }

    {
      TextSetColor(m_G, textColor);
      if(ButModeTranslate(m_G, P_GLUT_SINGLE_LEFT, 0) == cButModePickAtom) {
        TextDrawStrAt(m_G, "Picking ", x, y, orthoCGO);
        TextSetColor(m_G, TextColor3);
        TextDrawStrAt(m_G, "Atoms (and Joints)", x + DIP2PIXEL(64), y, orthoCGO);
      } else {
        TextDrawStrAt(m_G, "Selecting ", x, y, orthoCGO);
        TextSetColor(m_G, TextColor3);
        switch (SettingGetGlobal_i(m_G, cSetting_mouse_selection_mode)) {
        case 0:
          TextDrawStrAt(m_G, "Atoms", x + DIP2PIXEL(80), y, orthoCGO);
          break;
        case 1:
          TextDrawStrAt(m_G, "Residues", x + DIP2PIXEL(80), y, orthoCGO);
          break;
        case 2:
          TextDrawStrAt(m_G, "Chains", x + DIP2PIXEL(80), y, orthoCGO);
          break;
        case 3:
          TextDrawStrAt(m_G, "Segments", x + DIP2PIXEL(80), y, orthoCGO);
          break;
        case 4:
          TextDrawStrAt(m_G, "Objects", x + DIP2PIXEL(80), y, orthoCGO);
          break;
        case 5:
          TextDrawStrAt(m_G, "Molecules", x + DIP2PIXEL(80), y, orthoCGO);
          break;
        case 6:
          TextDrawStrAt(m_G, "C-alphas", x + DIP2PIXEL(80), y, orthoCGO);
          break;
        }
      }
    }
  }
  if (!orthoCGO || !(SettingGetGlobal_b(m_G, cSetting_show_frame_rate) || MoviePlaying(m_G))) {
    ButModeDrawFastImpl(this, true, orthoCGO);
  }
}

bool CButMode::fastDraw(CGO* orthoCGO){
  return ButModeDrawFastImpl(this, false, orthoCGO);
}

static bool ButModeDrawFastImpl(Block * block, short definitely , CGO *orthoCGO)
{
  PyMOLGlobals *G = block->m_G;
  CButMode *I = block->m_G->ButMode;
  int x, y;
  float *textColor = I->TextColor;
  float *textColor2 = I->TextColor2;

  if (!definitely && (!(SettingGetGlobal_b(G, cSetting_show_frame_rate) || MoviePlaying(G)))) {
    return false;
  }

  x = I->rect.left + cButModeLeftMargin;
  y = I->rect.bottom + cButModeLineHeight + cButModeBottomMargin;
  
  TextSetColor(G, I->TextColor);
  y -= cButModeLineHeight;
    {
        int buffer;
        /* TODO : Why do we only do this for the back right buffer,
         for performance? */
        glGetIntegerv(GL_DRAW_BUFFER, (GLint *) & buffer);
        if(buffer != GL_BACK_RIGHT) {
            if(I->Delay <= 0.0F) {
                if(I->Samples > 0.0F)
                    I->RateShown = (I->Rate / I->Samples);
                else
                    I->RateShown = 0.0F;
                I->Delay = 0.2F;
            }
        }
    }

  {
    int has_movie = false;
    int frame_rate = SettingGetGlobal_b(G, cSetting_show_frame_rate);
    int nf;
    char rateStr[255];
    nf = SceneGetNFrame(G, &has_movie);
    if(nf == 0)
      nf = 1;
    TextSetColor(G, textColor);
    if(has_movie) {
      TextDrawStrAt(G, "Frame ", x, y, orthoCGO);
    } else {
      TextDrawStrAt(G, "State ", x, y, orthoCGO);
    }
    TextSetColor(G, textColor2);
    snprintf(rateStr, std::size(rateStr), "%4d/%4d ", SceneGetFrame(G) + 1, nf);
    TextDrawStrAt(G, rateStr, x + DIP2PIXEL(48), y, orthoCGO);
    if(frame_rate) {
      snprintf(rateStr, std::size(rateStr),"%5.1f",I->RateShown);
      TextDrawStrAt(G, rateStr, x + DIP2PIXEL(144), y, orthoCGO);
      TextSetColor(G, textColor);
      TextDrawStrAt(G, "Hz ", x + DIP2PIXEL(192), y, orthoCGO);
      TextSetColor(G, textColor2);
    } else if(has_movie) {
      TextSetColor(G, textColor);
      TextDrawStrAt(G, "State ", x + DIP2PIXEL(128), y, orthoCGO);
      TextSetColor(G, textColor2);
      snprintf(rateStr, std::size(rateStr)," %4d",SceneGetState(G)+1);
      TextDrawStrAt(G, rateStr, x + DIP2PIXEL(168), y, orthoCGO);
    } else if(frame_rate) {
    }
  }
  return true;
}


/*========================================================================*/
int ButModeInit(PyMOLGlobals * G)
{
  CButMode *I = nullptr;
  if((I = (G->ButMode = new CButMode(G)))) {

    int a;

    I->Rate = 0.0;
    I->Samples = 0.0;
    I->RateShown = 0.0;
    I->Delay = 0.0;
    I->DeferCnt = 0;
    I->DeferTime = 0.0F;
    I->NCode = cButModeCount;
    I->NBut = cButModeInputCount;

    for(a = 0; a < I->NBut; a++) {
      I->Mode[a] = -1;
    }

    UtilNCopy(I->Code[cButModeRotXYZ], "Rota ", std::size(I->Code[cButModeRotXYZ]));
    UtilNCopy(I->Code[cButModeRotZ], "RotZ ", std::size(I->Code[cButModeRotZ]));
    UtilNCopy(I->Code[cButModeTransXY], "Move ", std::size(I->Code[cButModeTransXY]));
    UtilNCopy(I->Code[cButModeTransZ], "MovZ ", std::size(I->Code[cButModeTransZ]));
    UtilNCopy(I->Code[cButModeClipNF], "Clip ", std::size(I->Code[cButModeClipNF]));
    UtilNCopy(I->Code[cButModeClipN], "ClpN ", std::size(I->Code[cButModeClipN]));
    UtilNCopy(I->Code[cButModeClipF], "ClpF ", std::size(I->Code[cButModeClipF]));
    UtilNCopy(I->Code[cButModePickAtom], "PkAt ", std::size(I->Code[cButModePickAtom]));
    UtilNCopy(I->Code[cButModePickBond], "PkBd ", std::size(I->Code[cButModePickBond]));
    UtilNCopy(I->Code[cButModeTorFrag], "TorF ", std::size(I->Code[cButModeTorFrag]));
    UtilNCopy(I->Code[cButModeRotFrag], "RotF ", std::size(I->Code[cButModeRotFrag]));
    UtilNCopy(I->Code[cButModeMovFrag], "MovF ", std::size(I->Code[cButModeMovFrag]));
    UtilNCopy(I->Code[cButModeLB], " lb  ", std::size(I->Code[cButModeLB]));
    UtilNCopy(I->Code[cButModeMB], " mb  ", std::size(I->Code[cButModeMB]));
    UtilNCopy(I->Code[cButModeRB], " rb  ", std::size(I->Code[cButModeRB]));
    UtilNCopy(I->Code[cButModeAddToLB], "+lb  ", std::size(I->Code[cButModeAddToLB]));
    UtilNCopy(I->Code[cButModeAddToMB], "+mb  ", std::size(I->Code[cButModeAddToMB]));
    UtilNCopy(I->Code[cButModeAddToRB], "+rb  ", std::size(I->Code[cButModeAddToRB]));
    UtilNCopy(I->Code[cButModeOrigAt], "Orig ", std::size(I->Code[cButModeOrigAt]));
    UtilNCopy(I->Code[cButModeRectAdd], "+lBx ", std::size(I->Code[cButModeRectAdd]));
    UtilNCopy(I->Code[cButModeRectSub], "-lBx ", std::size(I->Code[cButModeRectSub]));
    UtilNCopy(I->Code[cButModeRect], "lbBx ", std::size(I->Code[cButModeRect]));
    UtilNCopy(I->Code[cButModeNone], "  -  ", std::size(I->Code[cButModeNone]));
    UtilNCopy(I->Code[cButModeCent], "Cent ", std::size(I->Code[cButModeCent]));
    UtilNCopy(I->Code[cButModePkTorBnd], "PkTB ", std::size(I->Code[cButModePkTorBnd]));
    UtilNCopy(I->Code[cButModeScaleSlab], "Slab ", std::size(I->Code[cButModeScaleSlab]));
    UtilNCopy(I->Code[cButModeMoveSlab], "MovS ", std::size(I->Code[cButModeMoveSlab]));
    UtilNCopy(I->Code[cButModePickAtom1], "Pk1  ", std::size(I->Code[cButModePickAtom1]));
    UtilNCopy(I->Code[cButModeMoveAtom], "MovA ", std::size(I->Code[cButModeMoveAtom]));
    UtilNCopy(I->Code[cButModeMenu], "Menu ", std::size(I->Code[cButModeMenu]));
    UtilNCopy(I->Code[cButModeSeleSet], "Sele ", std::size(I->Code[cButModeSeleSet]));
    UtilNCopy(I->Code[cButModeSeleToggle], "+/-  ", std::size(I->Code[cButModeSeleToggle]));
    UtilNCopy(I->Code[cButModeSeleAddBox], "+Box ", std::size(I->Code[cButModeSeleAddBox]));
    UtilNCopy(I->Code[cButModeSeleSubBox], "-Box ", std::size(I->Code[cButModeSeleSubBox]));
    UtilNCopy(I->Code[cButModeMoveSlabAndZoom], "MvSZ ", std::size(I->Code[cButModeMoveSlabAndZoom]));
    UtilNCopy(I->Code[cButModeSimpleClick], "Clik ", std::size(I->Code[cButModeSimpleClick]));
    UtilNCopy(I->Code[cButModeRotDrag], "RotD ", std::size(I->Code[cButModeRotDrag]));
    UtilNCopy(I->Code[cButModeMovDrag], "MovD ", std::size(I->Code[cButModeMovDrag]));
    UtilNCopy(I->Code[cButModeMovDragZ], "MvDZ ", std::size(I->Code[cButModeMovDragZ]));
    UtilNCopy(I->Code[cButModeRotObj], "RotO ", std::size(I->Code[cButModeRotObj]));
    UtilNCopy(I->Code[cButModeMovObj], "MovO ", std::size(I->Code[cButModeMovObj]));
    UtilNCopy(I->Code[cButModeMovObjZ], "MvOZ ", std::size(I->Code[cButModeMovObjZ]));
    UtilNCopy(I->Code[cButModeMovFragZ], "MvFZ ", std::size(I->Code[cButModeMovFragZ]));
    UtilNCopy(I->Code[cButModeMoveAtomZ], "MvAZ ", std::size(I->Code[cButModeMoveAtomZ]));
    UtilNCopy(I->Code[cButModeDragMol], "DrgM ", std::size(I->Code[cButModeDragMol]));
    UtilNCopy(I->Code[cButModeRotView], "RotV ", std::size(I->Code[cButModeRotView]));
    UtilNCopy(I->Code[cButModeMovView], "MovV ", std::size(I->Code[cButModeMovView]));
    UtilNCopy(I->Code[cButModeMovViewZ], "MvVZ ", std::size(I->Code[cButModeMovViewZ]));
    UtilNCopy(I->Code[cButModeDragObj], "DrgO ", std::size(I->Code[cButModeDragObj]));
    UtilNCopy(I->Code[cButModeInvMoveSlabAndZoom], "IMSZ ", std::size(I->Code[cButModeInvMoveSlabAndZoom]));
    UtilNCopy(I->Code[cButModeInvTransZ], "IMvZ ", std::size(I->Code[cButModeInvTransZ]));
    UtilNCopy(I->Code[cButModeSeleSetBox], " Box ", std::size(I->Code[cButModeSeleSetBox]));
    UtilNCopy(I->Code[cButModeInvRotZ], "IRtZ ", std::size(I->Code[cButModeInvRotZ]));
    UtilNCopy(I->Code[cButModeRotL], "RotL " , std::size(I->Code[cButModeRotL]));
    UtilNCopy(I->Code[cButModeMovL], "MovL " , std::size(I->Code[cButModeMovL]));
    UtilNCopy(I->Code[cButModeMvzL], "MvzL " , std::size(I->Code[cButModeMvzL]));

    I->active = true;

    I->TextColor[0] = 0.2F;
    I->TextColor[1] = 1.0F;
    I->TextColor[2] = 0.2F;

    I->TextColor1[0] = 0.5F;
    I->TextColor1[1] = 0.5F;
    I->TextColor1[2] = 1.0F;

    I->TextColor2[0] = 0.8F;
    I->TextColor2[1] = 0.8F;
    I->TextColor2[2] = 0.8F;

    I->TextColor3[0] = 1.0F;
    I->TextColor3[1] = 0.5F;
    I->TextColor3[2] = 0.5F;

    OrthoAttach(G, I, cOrthoTool);
    return 1;
  } else
    return 0;
}


/*========================================================================*/
int ButModeCheckPossibleSingleClick(PyMOLGlobals * G, int button, int mod)
{
  int click_button = -1;
  switch (button) {
  case P_GLUT_LEFT_BUTTON:
    click_button = P_GLUT_SINGLE_LEFT;
    break;
  case P_GLUT_MIDDLE_BUTTON:
    click_button = P_GLUT_SINGLE_MIDDLE;
    break;
  case P_GLUT_RIGHT_BUTTON:
    click_button = P_GLUT_SINGLE_RIGHT;
    break;
  }
  if(click_button < 0)
    return false;
  else
    return (ButModeTranslate(G, click_button, mod) >= 0);
}

int ButModeTranslate(PyMOLGlobals * G, int button, int mod)
{
  int mode = cButModeNothing;
  CButMode *I = G->ButMode;
  switch (button) {
  case P_GLUT_LEFT_BUTTON:
    mode = 0;
    break;
  case P_GLUT_MIDDLE_BUTTON:
    mode = 1;
    break;
  case P_GLUT_RIGHT_BUTTON:
    mode = 2;
    break;
  case P_GLUT_BUTTON_SCROLL_FORWARD:
  case P_GLUT_BUTTON_SCROLL_BACKWARD:
    switch (mod) {
    case 0:
      mode = 12;
      break;
    case cOrthoSHIFT:
      mode = 13;
      break;
    case cOrthoCTRL:
      mode = 14;
      break;
    case (cOrthoCTRL + cOrthoSHIFT):
      mode = 15;
    }
    mod = 0;
    switch (I->Mode[mode]) {
    case cButModeScaleSlab:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeScaleSlabExpand;
      } else {
        return cButModeScaleSlabShrink;
      }
      break;
    case cButModeMoveSlab:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeMoveSlabForward;
      } else {
        return cButModeMoveSlabBackward;
      }
      break;
    case cButModeMoveSlabAndZoom:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeMoveSlabAndZoomForward;
      } else {
        return cButModeMoveSlabAndZoomBackward;
      }
      break;
    case cButModeInvMoveSlabAndZoom:
      if(button != P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeMoveSlabAndZoomForward;
      } else {
        return cButModeMoveSlabAndZoomBackward;
      }
      break;
    case cButModeTransZ:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeZoomForward;
      } else {
        return cButModeZoomBackward;
      }
      break;
    case cButModeInvTransZ:
      if(button != P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeZoomForward;
      } else {
        return cButModeZoomBackward;
      }
      break;
    }
    return -1;
    break;
  case P_GLUT_DOUBLE_LEFT:
  case P_GLUT_DOUBLE_MIDDLE:
  case P_GLUT_DOUBLE_RIGHT:
  case P_GLUT_SINGLE_LEFT:
  case P_GLUT_SINGLE_MIDDLE:
  case P_GLUT_SINGLE_RIGHT:
    switch (button) {
    case P_GLUT_DOUBLE_LEFT:
      mode = 16;
      break;
    case P_GLUT_DOUBLE_MIDDLE:
      mode = 17;
      break;
    case P_GLUT_DOUBLE_RIGHT:
      mode = 18;
      break;
    case P_GLUT_SINGLE_LEFT:
      mode = 19;
      break;
    case P_GLUT_SINGLE_MIDDLE:
      mode = 20;
      break;
    case P_GLUT_SINGLE_RIGHT:
      mode = 21;
      break;
    }
    switch (mod) {
    case cOrthoSHIFT:
      mode += 6;
      break;
    case cOrthoCTRL:
      mode += 12;
      break;
    case (cOrthoCTRL + cOrthoSHIFT):
      mode += 18;
      break;
    case cOrthoALT:
      mode += 24;
      break;
    case (cOrthoALT + cOrthoSHIFT):
      mode += 30;
      break;
    case (cOrthoALT + cOrthoCTRL):
      mode += 36;
      break;
    case (cOrthoALT + cOrthoCTRL + cOrthoSHIFT):
      mode += 42;
      break;
    }
    mod = 0;
    break;
  }
  switch (mod) {
  case 0:
    break;
  case cOrthoSHIFT:
    mode += 3;
    break;
  case cOrthoCTRL:
    mode += 6;
    break;
  case (cOrthoCTRL + cOrthoSHIFT):
    mode += 9;
    break;
  case cOrthoALT:
    mode += 68;
    break;
  case (cOrthoALT + cOrthoSHIFT):
    mode += 71;
    break;
  case (cOrthoALT + cOrthoCTRL):
    mode += 74;
    break;
  case (cOrthoALT + cOrthoCTRL + cOrthoSHIFT):
    mode += 77;
    break;
  }
  return (I->Mode[mode]);
}
