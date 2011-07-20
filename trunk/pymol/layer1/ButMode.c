
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

#define cButModeLineHeight 12
#define cButModeLeftMargin 2
#define cButModeTopMargin 1

struct _CButMode {
  Block *Block;
  CodeType Code[cButModeCount + 1];
  int NCode;
  int Mode[cButModeInputCount];
  int NBut;
  float Rate;
  float RateShown;
  float Samples, Delay;
  float TextColor1[3];
  float TextColor2[3];
  float TextColor3[3];
  int DeferCnt;
  float DeferTime;
};


/*========================================================================*/
Block *ButModeGetBlock(PyMOLGlobals * G)
{
  register CButMode *I = G->ButMode;
  {
    return (I->Block);
  }
}

int ButModeGetHeight(PyMOLGlobals * G)
{
  if(SettingGetGlobal_b(G, cSetting_mouse_grid))
    return 124;
  else
    return 40;
}


/*========================================================================*/
int ButModeGet(PyMOLGlobals * G, int button)
{
  register CButMode *I = G->ButMode;
  if((button >= 0) && (button < I->NBut)) {
    return I->Mode[button];
  }
  return 0;
}

void ButModeSet(PyMOLGlobals * G, int button, int action)
{
  register CButMode *I = G->ButMode;
  if((button >= 0) && (button < I->NBut) && (action >= 0) && (action < I->NCode)) {
    I->Mode[button] = action;
    OrthoDirty(G);
  }
}


/*========================================================================*/
void ButModeSetRate(PyMOLGlobals * G, float interval)
{
  register CButMode *I = G->ButMode;

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
  register CButMode *I = G->ButMode;
  I->Samples = 0.0;
  I->Rate = 0.0;
  I->RateShown = 0;
  I->Delay = 0;
}


/*========================================================================*/
void ButModeFree(PyMOLGlobals * G)
{
  register CButMode *I = G->ButMode;
  OrthoFreeBlock(G, I->Block);
  FreeP(G->ButMode);
}


/*========================================================================*/
static int ButModeClick(Block * block, int button, int x, int y, int mod)
{
  PyMOLGlobals *G = block->G;
  int dy = (y - block->rect.bottom) / cButModeLineHeight;
  int dx = (x - block->rect.left);
  int forward = (dx > ((block->rect.right - block->rect.left) / 2));
  /*  register CButMode *I=block->G->ButMode; */
  if(dy < 2) {
    if(ButModeTranslate(G, P_GLUT_SINGLE_LEFT, 0) != cButModePickAtom) {
      switch (mod) {
      case cOrthoSHIFT:
        forward = !forward;
        break;
      }
      if(!forward) {
        PLog(G, "cmd.mouse('select_backward')", cPLog_pym);
        OrthoCommandIn(G, "mouse select_backward,quiet=1");
      } else {
        PLog(G, "cmd.mouse('select_forward')", cPLog_pym);
        OrthoCommandIn(G, "mouse select_forward,quiet=1");
      }
    }
  } else {
    if(button == P_GLUT_RIGHT_BUTTON) {
      MenuActivate0Arg(G,x,y,x,y,false,"mouse_config");
    } else {
      switch (mod) {
      case cOrthoSHIFT:
        forward = !forward;
        break;
      }
      if(!forward) {
        PLog(G, "cmd.mouse('backward')", cPLog_pym);
        OrthoCommandIn(G, "mouse backward,quiet=1");
      } else {
        PLog(G, "cmd.mouse('forward')", cPLog_pym);
        OrthoCommandIn(G, "mouse forward,quiet=1");
      }
    }
  }
  return (1);
}


/*========================================================================*/
static void ButModeDraw(Block * block)
{
  PyMOLGlobals *G = block->G;
  register CButMode *I = G->ButMode;
  int x, y, a;
  char rateStr[255];
  int mode;
  int nf;
  float *textColor = I->Block->TextColor;
  float *textColor2 = I->TextColor2;

#define BLANK_STR "     "

  if(G->HaveGUI && G->ValidContext && ((block->rect.right - block->rect.left) > 6)) {
    if(SettingGetGlobal_b(G, cSetting_internal_gui_mode) == 0) {
      glColor3fv(I->Block->BackColor);
      BlockFill(I->Block);
      BlockDrawLeftEdge(I->Block);
    } else {
      BlockDrawLeftEdge(I->Block);
      glColor3f(0.5, 0.5, 0.5);
      BlockDrawTopEdge(I->Block);
      textColor2 = OrthoGetOverlayColor(G);
      textColor = textColor2;
    }

    x = I->Block->rect.left + cButModeLeftMargin;
    y = (I->Block->rect.top - cButModeLineHeight) - cButModeTopMargin;

    TextSetColor(G, textColor);
    TextDrawStrAt(G, "Mouse Mode ", x + 1, y);
    TextSetColor(G, I->TextColor3);
    TextDrawStrAt(G, SettingGetGlobal_s(G, cSetting_button_mode_name), x + 88, y);
    /*    TextDrawStrAt(G,"2-Bttn Selecting",x+88,y); */
    y -= cButModeLineHeight;

    if(SettingGetGlobal_b(G, cSetting_mouse_grid)) {

      TextSetColor(G, I->TextColor3);
      TextDrawStrAt(G, "Buttons", x + 6, y);
      TextSetColor(G, I->TextColor1);
      /*    TextDrawStrAt(G,"  Left Mddl Rght Scrl",x+48,y); */
      TextDrawStrAt(G, "    L    M    R  Wheel", x + 43, y);

      y -= cButModeLineHeight;
      /*    glColor3fv(I->Block->TextColor);
         TextDrawStrAt(G,"K",x,y-4); */
      TextSetColor(G, I->TextColor3);
      TextDrawStrAt(G, "&", x + 12, y);
      TextDrawStrAt(G, "Keys", x + 24, y);
      TextSetColor(G, textColor2);

      TextSetPos2i(G, x + 64, y);
      for(a = 0; a < 3; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(G, BLANK_STR);
        else
          TextDrawStr(G, I->Code[mode]);
      }
      mode = I->Mode[12];
      if(mode < 0)
        TextDrawStr(G, BLANK_STR);
      else
        TextDrawStr(G, I->Code[mode]);

      y -= cButModeLineHeight;
      /*    TextSetColor(G,I->Block->TextColor);
         TextDrawStrAt(G,"e",x+5,y-1); */
      TextSetColor(G, I->TextColor1);

      TextSetColor(G, I->TextColor1);
      TextDrawStrAt(G, "Shft ", x + 24, y);
      TextSetColor(G, textColor2);
      TextSetPos2i(G, x + 64, y);
      for(a = 3; a < 6; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(G, BLANK_STR);
        else
          TextDrawStr(G, I->Code[mode]);
      }
      mode = I->Mode[13];
      if(mode < 0)
        TextDrawStr(G, BLANK_STR);
      else
        TextDrawStr(G, I->Code[mode]);

      y -= cButModeLineHeight;
      /*    glColor3fv(I->Block->TextColor);
         TextDrawStrAt(G,"y",x+10,y+2); */
      TextSetColor(G, I->TextColor1);
      TextDrawStrAt(G, "Ctrl ", x + 24, y);
      TextSetColor(G, textColor2);
      TextSetPos2i(G, x + 64, y);
      for(a = 6; a < 9; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(G, BLANK_STR);
        else
          TextDrawStr(G, I->Code[mode]);
      }
      mode = I->Mode[14];
      if(mode < 0)
        TextDrawStr(G, BLANK_STR);
      else
        TextDrawStr(G, I->Code[mode]);
      y -= cButModeLineHeight;

      /*    glColor3fv(I->Block->TextColor);
         TextDrawStrAt(G,"s",x+15,y+3); */
      TextSetColor(G, I->TextColor1);
      TextSetColor(G, I->TextColor1);
      TextDrawStrAt(G, "CtSh ", x + 24, y);
      TextSetColor(G, textColor2);
      TextSetPos2i(G, x + 64, y);
      for(a = 9; a < 12; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(G, BLANK_STR);
        else
          TextDrawStr(G, I->Code[mode]);
      }
      mode = I->Mode[15];
      if(mode < 0)
        TextDrawStr(G, BLANK_STR);
      else
        TextDrawStr(G, I->Code[mode]);

      y -= cButModeLineHeight;

      TextSetColor(G, I->Block->TextColor);
      TextSetColor(G, I->TextColor1);
      TextDrawStrAt(G, " SnglClk", x - 8, y);
      TextSetColor(G, textColor2);
      TextSetPos2i(G, x + 64, y);
      for(a = 19; a < 22; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(G, BLANK_STR);
        else
          TextDrawStr(G, I->Code[mode]);
      }
      TextSetColor(G, I->Block->TextColor);
      y -= cButModeLineHeight;

      TextSetColor(G, I->Block->TextColor);
      TextSetColor(G, I->TextColor1);
      TextDrawStrAt(G, " DblClk", x, y);
      TextSetColor(G, textColor2);
      TextSetPos2i(G, x + 64, y);
      for(a = 16; a < 19; a++) {
        mode = I->Mode[a];
        if(mode < 0)
          TextDrawStr(G, BLANK_STR);
        else
          TextDrawStr(G, I->Code[mode]);
      }
      TextSetColor(G, I->Block->TextColor);
      y -= cButModeLineHeight;

    }

    {
      TextSetColor(G, textColor);
      if(ButModeTranslate(G, P_GLUT_SINGLE_LEFT, 0) == cButModePickAtom) {
        TextDrawStrAt(G, "Picking ", x, y);
        TextSetColor(G, I->TextColor3);
        TextDrawStrAt(G, "Atoms (and Joints)", x + 64, y);
      } else {
        TextDrawStrAt(G, "Selecting ", x, y);
        TextSetColor(G, I->TextColor3);
        switch (SettingGetGlobal_i(G, cSetting_mouse_selection_mode)) {
        case 0:
          TextDrawStrAt(G, "Atoms", x + 80, y);
          break;
        case 1:
          TextDrawStrAt(G, "Residues", x + 80, y);
          break;
        case 2:
          TextDrawStrAt(G, "Chains", x + 80, y);
          break;
        case 3:
          TextDrawStrAt(G, "Segments", x + 80, y);
          break;
        case 4:
          TextDrawStrAt(G, "Objects", x + 80, y);
          break;
        case 5:
          TextDrawStrAt(G, "Molecules", x + 80, y);
          break;
        case 6:
          TextDrawStrAt(G, "C-alphas", x + 80, y);
          break;
        }
      }
    }

    TextSetColor(G, I->Block->TextColor);
    y -= cButModeLineHeight;
    {
      int buffer;
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
      nf = SceneGetNFrame(G, &has_movie);
      if(nf == 0)
        nf = 1;
      TextSetColor(G, textColor);
      if(has_movie) {
        TextDrawStrAt(G, "Frame ", x, y);
      } else {
        TextDrawStrAt(G, "State ", x, y);
      }
      TextSetColor(G, textColor2);
      sprintf(rateStr, "%4d/%4d ", SceneGetFrame(G) + 1, nf);
      TextDrawStrAt(G, rateStr, x + 48, y);
      if(frame_rate) {
        sprintf(rateStr,"%5.1f",I->RateShown);
        TextDrawStrAt(G, rateStr, x + 144, y);
        TextSetColor(G, textColor);
        TextDrawStrAt(G, "Hz ", x + 192, y);
        TextSetColor(G, textColor2);
      } else if(has_movie) {
        TextSetColor(G, textColor);
        TextDrawStrAt(G, "State ", x + 128, y);
        TextSetColor(G, textColor2);
        sprintf(rateStr," %4d",SceneGetState(G)+1);
        TextDrawStrAt(G, rateStr, x + 168, y);
      } else if(frame_rate) {
      }
    }
  }
}


/*========================================================================*/
int ButModeInit(PyMOLGlobals * G)
{
  register CButMode *I = NULL;
  if((I = (G->ButMode = Calloc(CButMode, 1)))) {

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

    strcpy(I->Code[cButModeRotXYZ], "Rota ");
    strcpy(I->Code[cButModeRotZ], "RotZ ");
    strcpy(I->Code[cButModeTransXY], "Move ");
    strcpy(I->Code[cButModeTransZ], "MovZ ");
    strcpy(I->Code[cButModeClipNF], "Clip ");
    strcpy(I->Code[cButModeClipN], "ClpN ");
    strcpy(I->Code[cButModeClipF], "ClpF ");
    strcpy(I->Code[cButModePickAtom], "PkAt ");
    strcpy(I->Code[cButModePickBond], "PkBd ");
    strcpy(I->Code[cButModeTorFrag], "TorF ");
    strcpy(I->Code[cButModeRotFrag], "RotF ");
    strcpy(I->Code[cButModeMovFrag], "MovF ");
    strcpy(I->Code[cButModeLB], " lb  ");
    strcpy(I->Code[cButModeMB], " mb  ");
    strcpy(I->Code[cButModeRB], " rb  ");
    strcpy(I->Code[cButModeAddToLB], "+lb  ");
    strcpy(I->Code[cButModeAddToMB], "+mb  ");
    strcpy(I->Code[cButModeAddToRB], "+rb  ");
    strcpy(I->Code[cButModeOrigAt], "Orig ");
    strcpy(I->Code[cButModeRectAdd], "+lBx ");
    strcpy(I->Code[cButModeRectSub], "-lBx ");
    strcpy(I->Code[cButModeRect], "lbBx ");
    strcpy(I->Code[cButModeNone], "  -  ");
    strcpy(I->Code[cButModeCent], "Cent ");
    strcpy(I->Code[cButModePkTorBnd], "PkTB ");
    strcpy(I->Code[cButModeScaleSlab], "Slab ");
    strcpy(I->Code[cButModeMoveSlab], "MovS ");
    strcpy(I->Code[cButModePickAtom1], "Pk1  ");
    strcpy(I->Code[cButModeMoveAtom], "MovA ");
    strcpy(I->Code[cButModeMenu], "Menu ");
    strcpy(I->Code[cButModeSeleSet], "Sele ");
    strcpy(I->Code[cButModeSeleToggle], "+/-  ");
    strcpy(I->Code[cButModeSeleAddBox], "+Box ");
    strcpy(I->Code[cButModeSeleSubBox], "-Box ");
    strcpy(I->Code[cButModeMoveSlabAndZoom], "MvSZ ");
    strcpy(I->Code[cButModeSimpleClick], "Clik ");
    strcpy(I->Code[cButModeRotDrag], "RotD ");
    strcpy(I->Code[cButModeMovDrag], "MovD ");
    strcpy(I->Code[cButModeMovDragZ], "MvDZ ");
    strcpy(I->Code[cButModeRotObj], "RotO ");
    strcpy(I->Code[cButModeMovObj], "MovO ");
    strcpy(I->Code[cButModeMovObjZ], "MvOZ ");
    strcpy(I->Code[cButModeMovFragZ], "MvFZ ");
    strcpy(I->Code[cButModeMoveAtomZ], "MvAZ ");
    strcpy(I->Code[cButModeDragMol], "DrgM ");
    strcpy(I->Code[cButModeRotView], "RotV ");
    strcpy(I->Code[cButModeMovView], "MovV ");
    strcpy(I->Code[cButModeMovViewZ], "MvVZ ");
    strcpy(I->Code[cButModeDragObj], "DrgO ");
    strcpy(I->Code[cButModeInvMoveSlabAndZoom], "IMSZ ");
    strcpy(I->Code[cButModeInvTransZ], "IMvZ ");
    strcpy(I->Code[cButModeSeleSetBox], " Box ");
    strcpy(I->Code[cButModeInvRotZ], "IRtZ ");

    I->Block = OrthoNewBlock(G, NULL);
    I->Block->fClick = ButModeClick;
    I->Block->fDraw = ButModeDraw;
    I->Block->fReshape = BlockReshape;
    I->Block->active = true;

    I->Block->TextColor[0] = 0.2F;
    I->Block->TextColor[1] = 1.0F;
    I->Block->TextColor[2] = 0.2F;

    I->TextColor1[0] = 0.5F;
    I->TextColor1[1] = 0.5F;
    I->TextColor1[2] = 1.0F;

    I->TextColor2[0] = 0.8F;
    I->TextColor2[1] = 0.8F;
    I->TextColor2[2] = 0.8F;

    I->TextColor3[0] = 1.0F;
    I->TextColor3[1] = 0.5F;
    I->TextColor3[2] = 0.5F;

    OrthoAttach(G, I->Block, cOrthoTool);
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
  register CButMode *I = G->ButMode;
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
