
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

#include "main.h"
#include "Base.h"
#include "Control.h"
#include "Scene.h"
#include "Executive.h"
#include "Setting.h"
#include "P.h"
#include "Seq.h"
#include"Movie.h"
#include "Text.h"
#include "Util.h"
#include "Ortho.h"
#include "CGO.h"

#define cControlBoxSize DIP2PIXEL(17)
#define cControlLeftMargin DIP2PIXEL(8)
#define cControlTopMargin DIP2PIXEL(2)
#define cControlSpacing DIP2PIXEL(2)
#define cControlInnerMargin  DIP2PIXEL(4)
#define cControlSpread DIP2PIXEL(6)
#define cControlSize DIP2PIXEL(160)

#define cControlButtons 7

#define cControlMinWidth 5

#define SDOF_QUEUE_MASK 0x1F


struct CControl : public Block {
  bool DragFlag {};
  int LastPos {};
  int ExtraSpace {};
  float ButtonColor[3] { 0.5f, 0.5f, 0.5f };
  float ActiveColor[3] { 0.65f, 0.65f, 0.65f };
  int Pressed { -1 };
  int Active { -1 };
  int SaveWidth { 0 };
  double LastClickTime {};
  int SkipRelease {};
  int NButton { 9 };

  /* not saved */

  int sdofActive {};
  double sdofLastIterTime {};
  int sdofMode {};
  float sdofTrans[3] {};
  float sdofRot[3] {};
  unsigned int sdofWroteTo {}, sdofReadFrom {};       /* queue synchronization fields */
  float sdofBuffer[(SDOF_QUEUE_MASK + 1) * 6] {};

  CControl(PyMOLGlobals * G) : Block(G) {};

  int click(int button, int x, int y, int mod) override;
  void draw(CGO* ortho) override;
  int drag(int x, int y, int mod) override;
  int release(int button, int x, int y, int mod) override;
  void reshape(int width, int height) override;
};

int ControlSdofButton(PyMOLGlobals * G, int button)
{
  CControl *I = G->Control;
  if(I) {
    if(button == 1) {           /* LEFT */
      if(I->sdofMode != SDOF_DRAG_MODE) {
        I->sdofMode = SDOF_DRAG_MODE;
        OrthoAddOutput(G, " SDOF: Drag mode.\n");
      } else {
        I->sdofMode = SDOF_NORMAL_MODE;
        OrthoAddOutput(G, " SDOF: Normal mode.\n");
      }
    } else if(button == 2) {    /* RIGHT */
      if(I->sdofMode != SDOF_CLIP_MODE) {
        I->sdofMode = SDOF_CLIP_MODE;
        OrthoAddOutput(G, " SDOF: Clip mode.\n");
      } else {
        I->sdofMode = SDOF_NORMAL_MODE;
        OrthoAddOutput(G, " SDOF: Normal mode.\n");
      }
    }
    OrthoDirty(G);
  }
  return 1;
}

int ControlSdofUpdate(PyMOLGlobals * G, float tx, float ty, float tz, float rx, float ry,
                      float rz)
{
  /* may be called asynchronously anytime after CControl has been initialized */

  CControl *I = G->Control;
  if(I) {
    float tol = 0.0001;
    bool active =
      fabs(tx) >= tol || fabs(ty) >= tol || fabs(tz) >= tol ||
      fabs(rx) >= tol || fabs(ry) >= tol || fabs(rz) >= tol;

    if(active) {
      int slot = (I->sdofWroteTo + 1) & SDOF_QUEUE_MASK;
      float *buffer = I->sdofBuffer + (6 * slot);

      buffer[0] = tx;
      buffer[1] = ty;
      buffer[2] = tz;
      buffer[3] = rx;
      buffer[4] = ry;
      buffer[5] = rz;

      I->sdofWroteTo = slot;

      {
        if(!I->sdofActive) {
          I->sdofLastIterTime = UtilGetSeconds(G);
        }
      }
      /*printf("wrote %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d %d %d\n",tx,ty,tz,rx,ry,rz,I->sdofReadFrom,I->sdofWroteTo,slot); */
    }
    I->sdofActive = active;
  }
  return 1;
}

int ControlSdofIterate(PyMOLGlobals * G)
{
  CControl *I = G->Control;
  if(I->sdofWroteTo != I->sdofReadFrom && I->sdofActive) {

    int slot = I->sdofWroteTo;

    /* new data available */

    float *buffer = I->sdofBuffer + (6 * slot);

    I->sdofTrans[0] = buffer[0];
    I->sdofTrans[1] = buffer[1];
    I->sdofTrans[2] = buffer[2];
    I->sdofRot[0] = buffer[3];
    I->sdofRot[1] = buffer[4];
    I->sdofRot[2] = buffer[5];

    I->sdofReadFrom = slot;

    double now = UtilGetSeconds(G);
    double delta = now - I->sdofLastIterTime;
    I->sdofLastIterTime = now;

    {
      /* suppress small amounts of combined motion using a truncated switching function */

      float len_rot = length3f(I->sdofRot);
      float len_trans = length3f(I->sdofTrans);
      float *dom, *sub;
      if(len_rot > len_trans) {
        dom = &len_rot;
        sub = &len_trans;
      } else {
        dom = &len_trans;
        sub = &len_rot;
      }

      {
        float expo = 2.0F;
        float trnc = 0.05F;
        float sub_fxn = (*sub) / (*dom);
        if(sub_fxn < trnc) {
          sub_fxn = 0.0F;
        } else if(sub_fxn < 0.5F) {
          sub_fxn = (sub_fxn - trnc) / (0.5F - trnc);
          sub_fxn = (float) pow(sub_fxn, expo);
        } else {
          sub_fxn = 1.0F - (float) pow(1.0F - sub_fxn, expo);
        }
        *dom = 1.0F;
        *sub = sub_fxn;
        scale3f(I->sdofTrans, len_trans, I->sdofTrans);
        scale3f(I->sdofRot, len_rot, I->sdofRot);
      }
    }
    /* translate */
    SceneTranslateScaled(G,
                         delta * I->sdofTrans[0],
                         -delta * I->sdofTrans[1], -delta * I->sdofTrans[2], I->sdofMode);

    /* rotate */
    SceneRotateScaled(G,
                      2.0F * delta * I->sdofRot[0],
                      -2.0F * delta * I->sdofRot[1],
                      -2.0F * delta * I->sdofRot[2], I->sdofMode);

    SceneDirty(G);
  }
  return 1;
}

int ControlRocking(PyMOLGlobals * G)
{
  if(G->Interrupt) {
    SettingSetGlobal_b(G, cSetting_rock, false);
  }
  return SettingGetGlobal_b(G, cSetting_rock);
}

void CControl::reshape(int width, int height)
{
  PyMOLGlobals *G = m_G;
  CControl *I = G->Control;
  Block::reshape(width, height);
  /* this is a pragmatic workaround for mac X11 where the nub gets
     hidden by the window expansion tab */

  if((rect.right - rect.left) < 20) {
    rect.top = rect.top + 10;
  }

  I->ExtraSpace = ((rect.right - rect.left) - cControlSize);
  if(I->ExtraSpace < 0)
    I->ExtraSpace = 0;
}

static int which_button(CControl * I, int x, int y)
{
  int result = -1;
  x -= I->rect.left + cControlLeftMargin;
  y -= I->rect.top - cControlTopMargin;
  if(x >= 0)
    if((y <= 0) && (y > (-cControlBoxSize))) {
      int control_width =
        I->rect.right - (I->rect.left + cControlLeftMargin);
      result = (I->NButton * x) / control_width;
    }
  return result;
}

int CControl::drag(int x, int y, int mod)
{
  int delta;
  int gui_width;
  PyMOLGlobals *G = m_G;
  CControl *I = G->Control;
  if(!I->SkipRelease) {
    delta = x - I->LastPos;
    delta /= DIP2PIXEL(1);
    if(I->DragFlag) {
      if(delta) {
        gui_width = SettingGetGlobal_i(G, cSetting_internal_gui_width) - delta;
        if(gui_width < cControlMinWidth)
          gui_width = cControlMinWidth;
        delta = SettingGetGlobal_i(G, cSetting_internal_gui_width) - gui_width;
        I->LastPos = x;
        I->SaveWidth = 0;
        SettingSetGlobal_i(G, cSetting_internal_gui_width, gui_width);
        OrthoReshape(G, -1, -1, false);
      }
    } else {
      I->Active = which_button(I, x, y);
      if(I->Active != I->Pressed)
        I->Active = -1;
      OrthoInvalidateDoDraw(G);
      OrthoDirty(G);
    }
  }
  return (1);
}

int CControl::release(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CControl *I = G->Control;

  int sel = 0;

  I->LastPos = x;
  sel = which_button(I, x, y);
  if(!I->SkipRelease) {
    switch (sel) {
    case 0:
      SceneSetFrame(G, 4, 0);
      PLog(G, "cmd.rewind()", cPLog_pym);
      break;
    case 1:
      SceneSetFrame(G, 5, -1);
      PLog(G, "cmd.back()", cPLog_pym);
      break;
    case 2:
      MoviePlay(G, cMovieStop);
      if(SettingGetGlobal_b(G, cSetting_sculpting))
        SettingSetGlobal_b(G, cSetting_sculpting, 0);
      if(SettingGetGlobal_b(G, cSetting_rock))
        SettingSetGlobal_b(G, cSetting_rock, false);
      OrthoDirty(G);
      PLog(G, "cmd.mstop()", cPLog_pym);
      break;
    case 3:
      if(!MoviePlaying(G)) {
        if(mod & cOrthoCTRL) {
          PLog(G, "cmd.rewind()", cPLog_pym);
          PLog(G, "cmd.mplay()", cPLog_pym);
          SceneSetFrame(G, 4, 0);
          MoviePlay(G, cMoviePlay);
        } else {
          PLog(G, "cmd.mplay()", cPLog_pym);
          MoviePlay(G, cMoviePlay);
        }
      } else {
        MoviePlay(G, cMovieStop);
        OrthoDirty(G);
        PLog(G, "cmd.mstop()", cPLog_pym);
      }
      break;
    case 4:
      SceneSetFrame(G, 5, 1);
      PLog(G, "cmd.forward()", cPLog_pym);
      break;
    case 5:
      if(mod & cOrthoCTRL) {
        SceneSetFrame(G, 3, 0);
        PLog(G, "cmd.middle()", cPLog_pym);
      } else {
        SceneSetFrame(G, 6, 0);
        PLog(G, "cmd.ending()", cPLog_pym);
      }
      break;
    case 6:
      if(SettingGetGlobal_b(G, cSetting_seq_view)) {
        SettingSetGlobal_b(G, cSetting_seq_view, 0);
        SeqChanged(G);
        PLog(G, "cmd.set('seq_view',0)", cPLog_pym);
      } else {
        SettingSetGlobal_b(G, cSetting_seq_view, 1);
        SeqChanged(G);
        PLog(G, "cmd.set('seq_view',1)", cPLog_pym);
      }
      OrthoDirty(G);
      break;
    case 7:
      SettingSetGlobal_b(G, cSetting_rock, !SettingGetGlobal_b(G, cSetting_rock));
      if(SettingGetGlobal_b(G, cSetting_rock)) {
        SceneRestartSweepTimer(G);
        PLog(G, "cmd.rock(1)", cPLog_pym);
      } else
        PLog(G, "cmd.rock(0)", cPLog_pym);
      SceneRestartFrameTimer(G);
      OrthoDirty(G);
      break;
    case 8:
      PLog(G, "cmd.full_screen()", cPLog_pym);
#ifdef _PYMOL_NOPY
      ExecutiveFullScreen(G, -1);
#else
      PParse(G, "full_screen");
#endif
      break;
    }
    OrthoDirty(G);
    OrthoUngrab(G);
    I->LastClickTime = UtilGetSeconds(G);
    I->DragFlag = false;
    I->Active = -1;
    I->Pressed = -1;
  }
  return (1);
}

Block *ControlGetBlock(PyMOLGlobals * G)
{
  CControl *I = G->Control;
  {
    return (I);
  }
}


/*========================================================================*/
int ControlIdling(PyMOLGlobals * G)
{
  CControl *I = G->Control;
  return (I->sdofActive ||
          MoviePlaying(G) ||
          SettingGetGlobal_b(G, cSetting_rock) || SettingGetGlobal_b(G, cSetting_sculpting));
}


/*========================================================================*/
void ControlFree(PyMOLGlobals * G)
{

  DeleteP(G->Control);
}


/*========================================================================*/
pymol::Result<bool> ControlRock(PyMOLGlobals * G, int mode)
{
  switch (mode) {
  case -2:
    break;
  case -1:
    SettingSetGlobal_b(G, cSetting_rock, !SettingGetGlobal_b(G, cSetting_rock));
    if(SettingGetGlobal_b(G, cSetting_rock)) {
      SceneRestartSweepTimer(G);
    }
    break;
  case 0:
    SettingSetGlobal_b(G, cSetting_rock, false);
    break;
  case 1:
    SettingSetGlobal_b(G, cSetting_rock, true);
    SceneRestartSweepTimer(G);
    break;
  }
  if(mode != -2) {
    SceneRestartFrameTimer(G);
    OrthoDirty(G);
  }
  return SettingGetGlobal_b(G, cSetting_rock);
}


/*========================================================================*/
int CControl::click(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CControl *I = G->Control;
  I->SkipRelease = false;
  if(x < (I->rect.left + cControlLeftMargin)) {
    y -= I->rect.top - cControlTopMargin;
    if((y <= 0) && (y > (-cControlBoxSize))) {
      double now = UtilGetSeconds(m_G);
      if((now - I->LastClickTime) < 0.35) {
        if(I->SaveWidth) {
          SettingSetGlobal_i(G, cSetting_internal_gui_width, I->SaveWidth);
          OrthoReshape(G, -1, -1, false);
          I->SaveWidth = 0;
        } else {
          I->SaveWidth = SettingGetGlobal_i(G, cSetting_internal_gui_width);
          SettingSetGlobal_i(G, cSetting_internal_gui_width, cControlMinWidth);
          OrthoReshape(G, -1, -1, false);
        }
        I->SkipRelease = true;
      } else {
        I->LastPos = x;
        OrthoGrab(G, this);
        I->DragFlag = true;
        I->LastClickTime = UtilGetSeconds(G);
      }
    }
  } else {
    I->Pressed = which_button(I, x, y);
    I->Active = I->Pressed;
    if(I->Pressed)
      OrthoGrab(G, this);
    OrthoDirty(G);
  }
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
    
    CGOColorv(orthoCGO, dark);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2 + 1, y2, 0.f);
    CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2 + h - 1, 0.f);
    CGOEnd(orthoCGO);
    
    CGOColorv(orthoCGO, inside);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2 + 1, y2 + 1, 0.f);
    CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
    CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, 0.f);
    CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, 0.f);
    CGOEnd(orthoCGO);
  } else {
    glColor3fv(light);
    glBegin(GL_POLYGON);
    glVertex2i(x2, y2);
    glVertex2i(x2, y2 + h);
    glVertex2i(x2 + w, y2 + h);
    glVertex2i(x2 + w, y2);
    glEnd();
    
    glColor3fv(dark);
    glBegin(GL_POLYGON);
    glVertex2i(x2 + 1, y2);
    glVertex2i(x2 + 1, y2 + h - 1);
    glVertex2i(x2 + w, y2 + h - 1);
    glVertex2i(x2 + w, y2);
    glEnd();
    
    glColor3fv(inside);
    glBegin(GL_POLYGON);
    glVertex2i(x2 + 1, y2 + 1);
    glVertex2i(x2 + 1, y2 + h - 1);
    glVertex2i(x2 + w - 1, y2 + h - 1);
    glVertex2i(x2 + w - 1, y2 + 1);
    glEnd();
  }
}


/*========================================================================*/
void CControl::draw(CGO* orthoCGO)
{
  PyMOLGlobals *G = m_G;
  CControl *I = this; // TODO: Remove I during Control refactor
  int x, y;
  int nButton = I->NButton;
  int but_num;
  float lightEdge[3] = { 0.65F, 0.65F, 0.65F };
  float darkEdge[3] = { 0.3F, 0.3F, 0.3F };
  float pushed[3] = { 0.8F, 0.8F, 0.8F };

  if(G->HaveGUI && G->ValidContext) {

    int control_width = rect.right - (rect.left + cControlLeftMargin);

    if (orthoCGO)
      CGOColorv(orthoCGO, BackColor);
#ifndef PURE_OPENGL_ES_2
    else
      glColor3fv(BackColor);
#endif
    fill(orthoCGO);
    if (orthoCGO)
      CGOColorv(orthoCGO, TextColor);
#ifndef PURE_OPENGL_ES_2
    else
      glColor3fv(TextColor);
#endif
    {
      int top, left, bottom, right;

      left = rect.left + 1;
      bottom = rect.bottom + 1;
      top = rect.top - (cControlTopMargin - 1);
      right = left + DIP2PIXEL(5);

      /* This draws the separator on the left side of the movie control buttons */
      if (orthoCGO){
	CGOColor(orthoCGO, 0.8F, 0.8F, 0.8F);
	CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	CGOVertex(orthoCGO, right, top, 0.f);
	CGOVertex(orthoCGO, right, bottom, 0.f);
	CGOVertex(orthoCGO, left, top, 0.f);
	CGOVertex(orthoCGO, left, bottom, 0.f);
	CGOEnd(orthoCGO);
	
	CGOColor(orthoCGO, 0.3F, 0.3F, 0.3F);
	CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	CGOVertex(orthoCGO, right, top - 1, 0.f);
	CGOVertex(orthoCGO, right, bottom, 0.f);
	CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
	CGOVertex(orthoCGO, left + 1, bottom, 0.f);
	CGOEnd(orthoCGO);
	
	CGOColorv(orthoCGO, I->ButtonColor);
	
	CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	CGOVertex(orthoCGO, right - 1, top - 1, 0.f);
	CGOVertex(orthoCGO, right - 1, bottom + 1, 0.f);
	CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
	CGOVertex(orthoCGO, left + 1, bottom + 1, 0.f);
	CGOEnd(orthoCGO);
      } else {
	glColor3f(0.8F, 0.8F, 0.8F);
	glBegin(GL_POLYGON);
	glVertex2i(right, top);
	glVertex2i(right, bottom);
	glVertex2i(left, bottom);
	glVertex2i(left, top);
	glEnd();
	
	glColor3f(0.3F, 0.3F, 0.3F);
	glBegin(GL_POLYGON);
	glVertex2i(right, top - 1);
	glVertex2i(right, bottom);
	glVertex2i(left + 1, bottom);
	glVertex2i(left + 1, top - 1);
	glEnd();
	
	glColor3fv(I->ButtonColor);
	
	glBegin(GL_POLYGON);
	glVertex2i(right - 1, top - 1);
	glVertex2i(right - 1, bottom + 1);
	glVertex2i(left + 1, bottom + 1);
	glVertex2i(left + 1, top - 1);
	glEnd();
      }

    y = rect.top - cControlTopMargin;

    for(but_num = 0; but_num < nButton; but_num++) {
      int but_width;
      int but_left;
      int but_bottom;
      int but_height;

      but_left =
        rect.left + cControlLeftMargin + (but_num * control_width) / nButton;
      but_width =
        (((but_num + 1) * control_width / nButton) -
         ((but_num) * control_width / nButton)) - 1;

      but_bottom = y - (cControlBoxSize - 1);
      but_height = cControlBoxSize;

      if(but_num == I->Active) {
        draw_button(but_left, but_bottom,
                    but_width, but_height, lightEdge, darkEdge, pushed ORTHOCGOARGVAR);
      } else if(((but_num == 6) && (SettingGetGlobal_b(G, cSetting_seq_view))) ||
                ((but_num == 3) && (MoviePlaying(G))) ||
                ((but_num == 7) && (SettingGetGlobal_b(G, cSetting_rock)))){
        draw_button(but_left, but_bottom,
                    but_width, but_height, lightEdge, darkEdge, I->ActiveColor ORTHOCGOARGVAR);
      } else {
        draw_button(but_left, but_bottom,
                    but_width, but_height, lightEdge, darkEdge, I->ButtonColor ORTHOCGOARGVAR);
      }

      if(control_width > 100) {
        x = but_left + (but_width - cControlBoxSize) / 2;

	if (orthoCGO)
	  CGOColorv(orthoCGO, TextColor);
#ifndef PURE_OPENGL_ES_2
	else
	  glColor3fv(TextColor);
#endif
        switch (but_num) {
        case 0:
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLES);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		      y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - (cControlBoxSize / 2), 0.f);
	    CGOEnd(orthoCGO);


	    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin-1.f, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin-1.f, y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_TRIANGLES);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glVertex2i(x + cControlInnerMargin, y - (cControlBoxSize / 2));
	    glEnd();
	    glBegin(GL_LINES);
	    glVertex2i(x + cControlInnerMargin, y - cControlInnerMargin);
	    glVertex2i(x + cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glEnd();
	  }
          break;
        case 1:
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLES);
	    CGOVertex(orthoCGO, x + cControlBoxSize / 2 + 2, y - (cControlBoxSize / 2), 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		      y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - (cControlBoxSize / 2), 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - (cControlBoxSize / 2), 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		      y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlBoxSize / 2 + 2, y - (cControlBoxSize / 2), 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_POLYGON);
	    glVertex2i(x + cControlBoxSize / 2 + 2, y - (cControlBoxSize / 2));
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - cControlInnerMargin);
	    glVertex2i(x + cControlInnerMargin, y - (cControlBoxSize / 2));
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glEnd();
	  }
          break;
        case 2:
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_POLYGON);
	    glVertex2i(x + cControlInnerMargin, y - cControlInnerMargin);
	    glVertex2i(x + cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - cControlInnerMargin);
	    glEnd();
	  }
          break;
        case 3:
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLES);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - cControlInnerMargin + 1, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin,
		      y - (cControlBoxSize - 1) + cControlInnerMargin - 1, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize) - cControlInnerMargin,
		      y - (cControlBoxSize / 2), 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_TRIANGLES);
	    glVertex2i(x + cControlInnerMargin, y - cControlInnerMargin + 1);
	    glVertex2i(x + cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin - 1);
	    glVertex2i(x + (cControlBoxSize) - cControlInnerMargin,
		       y - (cControlBoxSize / 2));
	    glEnd();
	  }
          break;
        case 4:
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLES);
	    CGOVertex(orthoCGO, x + cControlBoxSize / 2 - 2, y - (cControlBoxSize / 2), 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize / 2), 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize / 2), 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlBoxSize / 2 - 2, y - (cControlBoxSize / 2), 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_POLYGON);
	    glVertex2i(x + cControlBoxSize / 2 - 2, y - (cControlBoxSize / 2));
	    glVertex2i(x + cControlInnerMargin, y - cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize / 2));
	    glVertex2i(x + cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glEnd();
	  }
          break;
        case 5:
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLES);
	    CGOVertex(orthoCGO, x + cControlInnerMargin, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + cControlInnerMargin,
		      y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin,
		      y - (cControlBoxSize / 2), 0.f);
	    CGOEnd(orthoCGO);

	    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin - 1.f, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin, y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize - 1) - cControlInnerMargin - 1.f, y - (cControlBoxSize - 1) + cControlInnerMargin, 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_TRIANGLES);
	    glVertex2i(x + cControlInnerMargin, y - cControlInnerMargin);
	    glVertex2i(x + cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize / 2));
	    glEnd();
	    glBegin(GL_LINES);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize - 1) - cControlInnerMargin,
		       y - (cControlBoxSize - 1) + cControlInnerMargin);
	    glEnd();
	  }
          break;
        case 6:
	  TextSetColor(G, TextColor);
          TextDrawStrAt(G, "S", x + cControlInnerMargin,
                        y - cControlBoxSize + cControlInnerMargin + 1 ORTHOCGOARGVAR);
          break;
        case 7:
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLES);
	    CGOVertex(orthoCGO, x + (cControlBoxSize / 2) + cControlSpread, y - cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize / 2),
		      y - (cControlBoxSize) + cControlInnerMargin, 0.f);
	    CGOVertex(orthoCGO, x + (cControlBoxSize / 2) - cControlSpread, y - cControlInnerMargin, 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_POLYGON);
	    glVertex2i(x + (cControlBoxSize / 2) + cControlSpread, y - cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize / 2),
		       y - (cControlBoxSize) + cControlInnerMargin);
	    glVertex2i(x + (cControlBoxSize / 2) - cControlSpread, y - cControlInnerMargin);
	    glEnd();
	  }
          break;
        case 8:
	  TextSetColor(G, TextColor);
          TextDrawStrAt(G, "F", x + cControlInnerMargin,
                        y - cControlBoxSize + cControlInnerMargin + 1 ORTHOCGOARGVAR);
          break;
        }
      }
    }
    }
  }
}


/*========================================================================*/
int ControlInit(PyMOLGlobals * G)
{
  CControl *I = NULL;
  if((I = (G->Control = new CControl(G)))) {
    I->active = true;
    I->TextColor[0] = 1.0;
    I->TextColor[1] = 0.75;
    I->TextColor[2] = 0.75;
    OrthoAttach(G, I, cOrthoTool);
    I->LastClickTime = UtilGetSeconds(G);
    return 1;
  } else
    return 0;
}
