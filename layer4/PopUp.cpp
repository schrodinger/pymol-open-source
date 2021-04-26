
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
#include "CGO.h"

#define cPopUpLineHeight DIP2PIXEL(17)
#define cPopUpTitleHeight DIP2PIXEL(19)
#define cPopUpBarHeight DIP2PIXEL(4)
#define cPopUpCharWidth DIP2PIXEL(8)
#define cPopUpCharMargin DIP2PIXEL(2)
#define cPopUpCharLift DIP2PIXEL(2)

#define cChildDelay 0.25
#define cPassiveDelay 0.45
#define cDirtyDelay 0.05

struct CPopUp : public Block {
  Block *Parent {};
  Block *Child {};
  int ChildLine {};
  int LastX {}, LastY {};
  int StartX {}, StartY {};
  int Selected {};
  int Width {}, Height {};
  int NLine {};
  PyObject **Sub {};
  char **Command {};
  char **Text {};
  int *Code {};
  double ChildDelay {};
  double DirtyDelay {};
  double PassiveDelay {};
  int DirtyDelayFlag {};
  int NeverDragged {};
  int PlacementAffinity {};

  CPopUp(PyMOLGlobals * G) : Block(G){}

  void draw(CGO *orthoCGO) override;
  int drag(int x, int y, int mod) override;
  int release(int button, int x, int y, int mod) override;
  void reshape(int width, int height) override {};
};

static
int PopUpConvertY(CPopUp * I, int value, int mode);

/*========================================================================
 * If Sub[a] is a list, return it, otherwise assume it's callable so
 * call it, assign the result to Sub[a] and return the result.
 *
 * @param Sub : array of lists or callables
 * @param a   : array index
 */
static PyObject * SubGetItem(PyMOLGlobals * G, PyObject ** Sub, const int a) {
#ifndef _PYMOL_NOPY
  PyObject *elem = Sub[a];
  ok_assert(1, elem);

  if(!PyList_Check(elem)) {
    PBlock(G);
    elem = PyObject_CallObject(elem, NULL);

    if(PyErr_Occurred())
      PyErr_Print();

    Py_DECREF(Sub[a]);
    Sub[a] = elem;

    PUnblock(G);
  }

  return elem;

ok_except1:
#endif
  return NULL;
}

/*========================================================================*/
static Block *PopUpRecursiveFind(Block * block, int x, int y)
{
  PyMOLGlobals *G = block->m_G;
  CPopUp *I = (CPopUp *) block->reference;
  if(I->Child) {                /* favor the child */
    if(PopUpRecursiveFind(I->Child, x, y) == I->Child)
      return block;
  }
  if(block->recursiveFind(x, y) == block) {
    OrthoGrab(G, block);
    return block;
  }
  return NULL;
}


/*========================================================================*/
Block *PopUpNew(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                int passive, PyObject * list, Block * parent)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* assumes blocked threads (calls the Python C API) */

  int mx, a, l, cl, cmx;
  int dim[2];
  PyObject *elem;
  const char *str, *c;
  int blocked = PAutoBlock(G);
  auto ui_light_bg = SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting);
  CPopUp *I = new CPopUp(G);

  I->reference = (void *) I;
  I->active = false;
  I->TextColor[0] = 1.0F;
  I->TextColor[1] = 1.0F;
  I->TextColor[2] = 1.0F;

  I->BackColor[0] = 0.1F;
  I->BackColor[1] = 0.1F;
  I->BackColor[2] = 0.1F;

  if(ui_light_bg != InternalGUIMode::Default) {
    I->TextColor[0] = 0.0F;
    I->TextColor[1] = 0.0F;
    I->TextColor[2] = 0.0F;
    I->BackColor[0] = 1.0F;
    I->BackColor[1] = 1.0F;
    I->BackColor[2] = 1.0F;
  }

  I->Parent = parent;
  I->Child = NULL;
  I->NLine = PyList_Size(list);
  I->Text = NULL;
  I->Command = NULL;
  I->Code = NULL;
  I->Selected = -1;
  I->StartX = (I->LastX = last_x);
  I->StartY = (I->LastY = last_y);
  I->ChildDelay = UtilGetSeconds(G) + cChildDelay * 2.5;
  I->PassiveDelay = UtilGetSeconds(G) + cPassiveDelay;
  I->DirtyDelay = false;
  I->DirtyDelayFlag = false;
  I->NeverDragged = true;
  I->PlacementAffinity = 0;
  mx = 1;
  cmx = 1;
  for(a = 0; a < I->NLine; a++) {
    PyObject* item = PyList_GetItem(list, a);
    if (PyList_Size(item) < 2) {
      continue;
    }
    elem = PyList_GetItem(item, 1);
    l = PyString_Size(elem);
    str = PyString_AsString(elem);
    cl = l;
    c = str;
    while(*c) {
      if(TextStartsWithColorCode(c)) {  /* discount the markup */
          c += 3;
          cl -= 4;
      }
      c++;
    }
    if(cl > cmx)
      cmx = cl;
    if(l > mx)
      mx = l;
  }
  I->Width = (cmx * cPopUpCharWidth) + 2 * cPopUpCharMargin;

  dim[0] = I->NLine + 1;
  dim[1] = mx + 1;
  I->Text = (char **) UtilArrayCalloc((unsigned int *) (void *) dim, 2, 1);

  mx = 1;
  for(a = 0; a < I->NLine; a++) {
    PyObject* item = PyList_GetItem(list, a);
    if (PyList_Size(item) > 2) {
      PyObject* command = PyList_GetItem(item, 2);
      if(PyString_Check(command)) {
	l = PyString_Size(command);
	if(l > mx)
	  mx = l;
      }
    }
  }
  dim[0] = I->NLine + 1;
  dim[1] = mx + 1;
  I->Command = (char **) UtilArrayCalloc((unsigned int *) (void *) dim, 2, 1);

  I->Code = pymol::malloc<int>(I->NLine + 1);
  I->Sub = pymol::calloc<PyObject *>(I->NLine + 1);

  for(a = 0; a < I->NLine; a++) {
    PyObject *command;
    elem = PyList_GetItem(list, a);
    I->Code[a] = PyInt_AsLong(PyList_GetItem(elem, 0));
    if (I->Code[a] == 0) {
      continue;
    }
    strcpy(I->Text[a], PyString_AsString(PyList_GetItem(elem, 1)));
    if (I->Code[a] == 2) {
      continue;
    }
    command = PyList_GetItem(elem, 2);
    if(command) {
      if(PyString_Check(command)) {
	strcpy(I->Command[a], PyString_AsString(command));
      } else {
	Py_INCREF(command);
	I->Sub[a] = command;
      }
    }
  }

  /* compute height */

  I->Height = 1 * cPopUpCharMargin + PopUpConvertY(I, I->NLine, true);

  I->rect.top = y;
  I->rect.bottom = y - I->Height;
  I->rect.left = x - (I->Width) / 3;
  I->rect.right = x + (2 * I->Width) / 3;

  PopFitBlock(I);

  OrthoAttach(G, I, cOrthoTool);
  I->active = true;
  OrthoGrab(G, I);
  OrthoDirty(G);

  if(passive)
    PyMOL_SetPassive(G->PyMOL, true);

  PAutoUnblock(G, blocked);

  OrthoInvalidateDoDraw(G);  
  return I;
#endif

}


/*========================================================================*/
int PopUpConvertY(CPopUp * I, int value, int mode)
{
  int result;
  int a;
  int flag;
  if(mode) {
    result = 0;

    /* line to height */
    for(a = 0; a < I->NLine; a++) {
      if(a >= value)
        break;
      switch (I->Code[a]) {
      case 0:
        result += cPopUpBarHeight;
        break;
      case 1:
        result += cPopUpLineHeight;
        break;
      case 2:
        result += cPopUpTitleHeight;
        break;
      }
    }
  } else {
    flag = false;
    /* height to line */
    if(value < 0) {
      result = -1;
    } else {
      result = 0;
      for(a = 0; a < I->NLine; a++) {
        switch (I->Code[a]) {
        case 0:
          if(value < cPopUpBarHeight)
            flag = true;
          value -= cPopUpBarHeight;
          break;
        case 1:
          if(value < cPopUpLineHeight)
            flag = true;
          value -= cPopUpLineHeight;
          break;
        case 2:
          if(value < cPopUpLineHeight)
            flag = true;
          value -= cPopUpTitleHeight;
          break;
        }
        if(flag)
          break;
        result++;
      }
      if(!flag)
        result = -1;
      else if(result && !I->Code[result])
        result--;
    }
    /* height to line */
  }
  return (result);
}


/*========================================================================*/

static void PopUpDetachRecursiveChild(Block * block)
{
  PyMOLGlobals *G = block->m_G;
  CPopUp *I = (CPopUp *) block->reference;

  OrthoDetach(G, block);
  if(I->Child)
    PopUpDetachRecursiveChild(I->Child);
}


/*========================================================================*/

static void PopUpForgetChild(Block * block)
{

  CPopUp *I = (CPopUp *) block->reference;
  I->Child = NULL;
}

static void PopUpRecursiveDetach(Block * block)
{
  PyMOLGlobals *G = block->m_G;
  CPopUp *I = (CPopUp *) block->reference;
  OrthoDetach(G, block);
  if(I->Child)
    PopUpDetachRecursiveChild(I->Child);
  if(I->Parent) {
    PopUpForgetChild(I->Parent);
    PopUpRecursiveDetach(I->Parent);
  }
}

static void PopUpFree(Block * block)
{
  PyMOLGlobals *G = block->m_G;
  CPopUp *I = (CPopUp *) block->reference;

#ifndef _PYMOL_NOPY
  {                             /* purge code */
    int a;
    int blocked = PAutoBlock(G);
    for(a = 0; a < I->NLine; a++)
      if(I->Sub[a]) {
        Py_DECREF(I->Sub[a]);
      }
    PAutoUnblock(G, blocked);
  }
#endif

  OrthoDetach(G, I);
  FreeP(I->Sub);
  FreeP(I->Code);
  FreeP(I->Command);
  FreeP(I->Text);
  DeleteP(I);

}

static void PopUpRecursiveFree(Block * block)
{

  CPopUp *I = (CPopUp *) block->reference;

  if(I->Child)
    PopUpFree(I->Child);
  I->Child = NULL;
  if(I->Parent) {
    PopUpForgetChild(I->Parent);
    PopUpRecursiveFree(I->Parent);
  }
  PopUpFree(block);
}

static void PopUpFreeRecursiveChild(Block * block)
{
  CPopUp *I = (CPopUp *) block->reference;
  if(I->Child)
    PopUpFreeRecursiveChild(I->Child);
  I->Child = NULL;
  PopUpFree(block);
}


/*========================================================================*/
int CPopUp::release(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CPopUp *I = (CPopUp *) reference;
  int gone_passive = false;

  int scroll_dy = 10;
  switch (button) {
    case PYMOL_BUTTON_SCROLL_FORWARD:
      scroll_dy *= -1;
    case PYMOL_BUTTON_SCROLL_REVERSE:
      translate(0, scroll_dy);
      return 1;
  }

  if(I->NeverDragged) {
    if(I->PassiveDelay > UtilGetSeconds(G)) {
      gone_passive = true;
      I->PassiveDelay = UtilGetSeconds(G);      /* kill any further delay */
    }
  }
  if(!gone_passive) {
    if(!I->NeverDragged)
      drag(x, y, mod);

    /* go passive if we click and release on a sub-menu */

    if((I->Selected >= 0) && (I->Sub[I->Selected])) {
      if((x >= I->rect.left) && (x <= I->rect.right)) {
        gone_passive = true;
      }
    }
  }
  if(gone_passive) {
    PyMOL_SetPassive(G->PyMOL, true);
  } else {
    OrthoUngrab(G);
    PopUpRecursiveDetach(this);
    if(!I->NeverDragged)
      if((I->Selected >= 0) && (!I->Sub[I->Selected])) {
        PLog(G, I->Command[I->Selected], cPLog_pym);
        PParse(G, I->Command[I->Selected]);
        PFlush(G);
      }
    PopUpRecursiveFree(this);
  }
  OrthoDirty(G);
  return (1);
}


/*========================================================================*/
int CPopUp::drag(int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CPopUp *I = (CPopUp *) reference;

  int a;
  int was = I->Selected;

  if((!I->NeverDragged) && (((I->StartX - x) > 4) || ((I->StartY - y) > 4)))
    I->NeverDragged = false;

  I->LastX = x;
  I->LastY = y;

  x -= I->rect.left;
  y = (I->rect.top - cPopUpCharMargin) - y - 1;

  if((x < -2) || (x > (I->Width + 2))) {
    int handled_flag = false;
    if(I->Child) {
      if(PopUpRecursiveFind(I->Child, I->LastX, I->LastY) == I->Child) {
        I->Selected = I->ChildLine;
        handled_flag = true;
      }
    }
    if(!handled_flag) {
      if(I->Parent) {           /* are we back in the parent window? */
        I->Selected = -1;
        return I->Parent->drag(I->LastX, I->LastY, mod);
      } else if(!I->Child) {
        I->Selected = -1;
      }
    }
  } else {
    OrthoGrab(G, this);
    a = PopUpConvertY(I, y, false);
    if(I->NLine && (a == I->NLine))
      if((y - a * cPopUpLineHeight) < 4)
        a = I->NLine - 1;
    if((a < 0) || (a >= I->NLine))
      I->Selected = -1;
    else {
      if(I->Code[a] == 1) {
        if((I->Child) && (I->ChildLine != a)) {
          if(I->ChildDelay < UtilGetSeconds(G)) {
            PopUpDetachRecursiveChild(I->Child);
            PopUpFreeRecursiveChild(I->Child);
            I->Child = NULL;
            I->ChildLine = -1;
            OrthoDirty(G);
	    OrthoInvalidateDoDraw(G);
          } else {
            I->Selected = a;
          }
          PyMOL_NeedFakeDrag(G->PyMOL);
        }
      }

      if(I->Code[a] != 1)
        I->Selected = -1;
      else {
        /* activate submenu */
        PyObject * sub_a = SubGetItem(G, I->Sub, a);
        if(!sub_a) {
          // do nothing
        } else if(!I->Child) {
          I->ChildLine = a;
          if(I->ChildDelay > UtilGetSeconds(G)) {
            PyMOL_NeedFakeDrag(G->PyMOL);       /* keep coming back here... */
          } else {
            I->Child = PopUpNew(G, I->LastX - 300, I->LastY, I->LastX, I->LastY,
                                false, sub_a, I);
            {
              int target_y =
                rect.top - (PopUpConvertY(I, a, true) + cPopUpCharMargin);
              CPopUp *child = (CPopUp *) (I->Child->reference);
              if(child->NLine)
                if(child->Code[0] != 1)
                  target_y += cPopUpTitleHeight + 2;
              child->PlacementAffinity = PopPlaceChild(I->Child, rect.left - 5,
                                                       rect.right + 5, target_y,
                                                       I->PlacementAffinity);
            }

            OrthoGrab(G, this);
            I->ChildDelay = UtilGetSeconds(G) + cChildDelay;    /* leave child up for a while */
          }
          PyMOL_NeedFakeDrag(G->PyMOL); /* keep coming back here... */
        } else if(I->ChildLine == a) {  /* on correct line */
          I->ChildDelay = UtilGetSeconds(G) + cChildDelay;      /* keep child here for a while */
        }
        I->Selected = a;
      }
    }
  }
  /* delay updates, etc. so that child menus 
     can be comfortably accessed with sloppy mousing */

  if((I->Child) && (I->Selected != I->ChildLine))
    PyMOL_NeedFakeDrag(G->PyMOL);

  if(was != I->Selected) {

    I->NeverDragged = false;
    if(!I->Child) {
      /* we moved, so renew the child delay */
      I->ChildDelay = UtilGetSeconds(G) + cChildDelay;
      PyMOL_NeedFakeDrag(G->PyMOL);
    }

    if((I->Child) && (I->Selected != I->ChildLine)) {
      I->DirtyDelayFlag = true;
      I->DirtyDelay = UtilGetSeconds(G) + cDirtyDelay;
    }
    if(!I->DirtyDelayFlag){
      OrthoDirty(G);
      OrthoInvalidateDoDraw(G);
    }
  }
  if(I->DirtyDelayFlag && (I->DirtyDelay < UtilGetSeconds(G))) {
    I->DirtyDelayFlag = false;
    OrthoDirty(G);
    OrthoInvalidateDoDraw(G);
  }
  return (1);
}


/*========================================================================*/
void CPopUp::draw(CGO* orthoCGO)
{
  CPopUp *I = this; // TODO: Remove I during PopUp refactor
  PyMOLGlobals *G = m_G;
  int x, y, a, xx;
  char *c;

  if(G->HaveGUI && G->ValidContext) {

    if((I->Child) && (I->Selected != I->ChildLine))
      PyMOL_NeedFakeDrag(G->PyMOL);

    /* put raised border around pop-up menu */

    if (orthoCGO){
      CGOColor(orthoCGO, 0.2F, 0.2F, 0.4F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left - 2, rect.bottom - 2, 0.f);
      CGOVertex(orthoCGO, rect.right + 2, rect.bottom - 2, 0.f);
      CGOVertex(orthoCGO, rect.left - 2, rect.bottom + 1, 0.f);
      CGOVertex(orthoCGO, rect.right + 2, rect.bottom + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.2F, 0.2F, 0.4F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.left - 2, rect.bottom - 2);
      glVertex2i(rect.right + 2, rect.bottom - 2);
      glVertex2i(rect.right + 2, rect.bottom + 1);
      glVertex2i(rect.left - 2, rect.bottom + 1);
      glEnd();
    }

    if (orthoCGO){
      CGOColor(orthoCGO, 0.4F, 0.4F, 0.6F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left - 1, rect.bottom - 1, 0.f);
      CGOVertex(orthoCGO, rect.right + 1, rect.bottom - 1, 0.f);
      CGOVertex(orthoCGO, rect.left - 1, rect.bottom + 1, 0.f);
      CGOVertex(orthoCGO, rect.right + 1, rect.bottom + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.4F, 0.4F, 0.6F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.left - 1, rect.bottom - 1);
      glVertex2i(rect.right + 1, rect.bottom - 1);
      glVertex2i(rect.right + 1, rect.bottom + 1);
      glVertex2i(rect.left - 1, rect.bottom + 1);
      glEnd();
    }
    /* right */

    if (orthoCGO){
      CGOColor(orthoCGO, 0.2F, 0.2F, 0.4F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.right, rect.bottom - 2, 0.f);
      CGOVertex(orthoCGO, rect.right + 2, rect.bottom - 2, 0.f);
      CGOVertex(orthoCGO, rect.right, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.right + 2, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.2F, 0.2F, 0.4F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.right, rect.bottom - 2);
      glVertex2i(rect.right + 2, rect.bottom - 2);
      glVertex2i(rect.right + 2, rect.top);
      glVertex2i(rect.right, rect.top);
      glEnd();
    }

    if (orthoCGO){
      CGOColor(orthoCGO, 0.4F, 0.4F, 0.6F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.right, rect.bottom - 1, 0.f);
      CGOVertex(orthoCGO, rect.right + 1, rect.bottom - 1, 0.f);
      CGOVertex(orthoCGO, rect.right, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.right + 1, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.4F, 0.4F, 0.6F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.right, rect.bottom - 1);
      glVertex2i(rect.right + 1, rect.bottom - 1);
      glVertex2i(rect.right + 1, rect.top);
      glVertex2i(rect.right, rect.top);
      glEnd();
    }
    /* top */

    if (orthoCGO){
      CGOColor(orthoCGO, 0.5F, 0.5F, 0.7F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left - 2, rect.top + 2, 0.f);
      CGOVertex(orthoCGO, rect.right + 2, rect.top + 2, 0.f);
      CGOVertex(orthoCGO, rect.left - 2, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.right + 2, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.5F, 0.5F, 0.7F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.left - 2, rect.top + 2);
      glVertex2i(rect.right + 2, rect.top + 2);
      glVertex2i(rect.right + 2, rect.top);
      glVertex2i(rect.left - 2, rect.top);
      glEnd();
    }

    if (orthoCGO){
      CGOColor(orthoCGO, 0.6F, 0.6F, 0.8F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left - 1, rect.top + 1, 0.f);
      CGOVertex(orthoCGO, rect.right + 1, rect.top + 1, 0.f);
      CGOVertex(orthoCGO, rect.left - 1, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.right + 1, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.6F, 0.6F, 0.8F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.left - 1, rect.top + 1);
      glVertex2i(rect.right + 1, rect.top + 1);
      glVertex2i(rect.right + 1, rect.top);
      glVertex2i(rect.left - 1, rect.top);
      glEnd();
    }

    /* left */
    if (orthoCGO){
      CGOColor(orthoCGO, 0.5F, 0.5F, 0.7F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left - 2, rect.bottom - 2, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.bottom, 0.f);
      CGOVertex(orthoCGO, rect.left - 2, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.5F, 0.5F, 0.7F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.left - 2, rect.bottom - 2);
      glVertex2i(rect.left, rect.bottom);
      glVertex2i(rect.left, rect.top);
      glVertex2i(rect.left - 2, rect.top);
      glEnd();
    }

    if (orthoCGO){
      CGOColor(orthoCGO, 0.6F, 0.6F, 0.8F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left - 1, rect.bottom - 1, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.bottom - 1, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.left - 1, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.6F, 0.6F, 0.8F);
      glBegin(GL_POLYGON);
      glVertex2i(rect.left - 1, rect.bottom - 1);
      glVertex2i(rect.left, rect.bottom - 1);
      glVertex2i(rect.left, rect.top);
      glVertex2i(rect.left - 1, rect.top);
      glEnd();
    }

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

    if(I->Selected >= 0) {

      x = rect.left;
      y = rect.top - PopUpConvertY(I, I->Selected, true) - cPopUpCharMargin;

      y += 2;
    if (orthoCGO){
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, x, y, 0.f);
      CGOVertex(orthoCGO, x + I->Width - 1, y, 0.f);
      CGOVertex(orthoCGO, x, y - (cPopUpLineHeight + 3), 0.f);
      CGOVertex(orthoCGO, x + I->Width - 1, y - (cPopUpLineHeight + 3), 0.f);
      CGOEnd(orthoCGO);
    } else {
      glBegin(GL_POLYGON);
      glVertex2i(x, y);
      glVertex2i(x + I->Width - 1, y);
      glVertex2i(x + I->Width - 1, y - (cPopUpLineHeight + 3));
      glVertex2i(x, y - (cPopUpLineHeight + 3));
      glEnd();
    }
    }
    if(I->Code[0] == 2) {       /* menu name */

      if (SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) ==
          InternalGUIMode::Default)
      {
        if (orthoCGO)
          CGOColor(orthoCGO, 0.3F, 0.3F, 0.6F);
        else
          glColor3f(0.3F, 0.3F, 0.6F);
      } else {
	if (orthoCGO)
	  CGOColor(orthoCGO, 1.0F, 1.0F, 1.0F);
	else
	  glColor3f(1.0F, 1.0F, 1.0F);
      }
      x = rect.left;
      y = rect.top;

      if (orthoCGO){
	CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	CGOVertex(orthoCGO, x, y, 0.f);
	CGOVertex(orthoCGO, x + I->Width, y, 0.f);
	CGOVertex(orthoCGO, x, y - (cPopUpLineHeight + cPopUpCharMargin), 0.f);
	CGOVertex(orthoCGO, x + I->Width, y - (cPopUpLineHeight + cPopUpCharMargin), 0.f);
	CGOEnd(orthoCGO);
      } else {
	glBegin(GL_POLYGON);
	glVertex2i(x, y);
	glVertex2i(x + I->Width, y);
	glVertex2i(x + I->Width, y - (cPopUpLineHeight + cPopUpCharMargin));
	glVertex2i(x, y - (cPopUpLineHeight + cPopUpCharMargin));
	glEnd();
      }

      if (orthoCGO){
	CGOColor(orthoCGO, 0.2F, 0.2F, 0.4F);
	CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	CGOVertex(orthoCGO, x + I->Width - 1, y - (cPopUpLineHeight + cPopUpCharMargin), 0.f);
	CGOVertex(orthoCGO, x + I->Width - 1, y - (cPopUpLineHeight + cPopUpCharMargin) - 1.f, 0.f);
	CGOVertex(orthoCGO, x, y - (cPopUpLineHeight + cPopUpCharMargin), 0.f);
	CGOVertex(orthoCGO, x, y - (cPopUpLineHeight + cPopUpCharMargin) - 1.f, 0.f);
	CGOEnd(orthoCGO);
      } else {
	glColor3f(0.2F, 0.2F, 0.4F);
	glBegin(GL_LINES);
	glVertex2i(x + I->Width - 1, y - (cPopUpLineHeight + cPopUpCharMargin));
	glVertex2i(x, y - (cPopUpLineHeight + cPopUpCharMargin));
	glEnd();
      }
    }

    x = rect.left + cPopUpCharMargin;
    y = (rect.top - cPopUpLineHeight) - cPopUpCharMargin + 2;

    for(a = 0; a < I->NLine; a++) {
      auto text_color = (a == I->Selected) ? BackColor : TextColor;
      TextSetColor(G, text_color);
      if(I->Code[a]) {
        c = I->Text[a];
        xx = x;
        while(*c) {
          // note: previously also supported "\\+++red", but was never used
          if(TextSetColorFromCode(G, c, text_color)) {
            c += 4;
          }

          TextSetPos2i(G, xx, y + cPopUpCharLift);
          TextDrawChar(G, *(c++) ORTHOCGOARGVAR);
          xx = xx + DIP2PIXEL(8);
        }

        if(I->Sub[a]) {

	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	    CGOColor(orthoCGO, 0.4F, 0.4F, 0.4F);
	    CGOVertex(orthoCGO, rect.left - 3, y + ((cPopUpLineHeight)) - 4, 0.f);
	    CGOColor(orthoCGO, 0.4F, 0.4F, 0.4F);
	    CGOVertex(orthoCGO, rect.left - 3, y + 1, 0.f);
	    CGOColor(orthoCGO, 0.1F, 0.1F, 0.1F);
	    CGOVertex(orthoCGO, rect.left, y + ((cPopUpLineHeight)) - 4, 0.f);
	    CGOColor(orthoCGO, 0.1F, 0.1F, 0.1F);
	    CGOVertex(orthoCGO, rect.left, y + 1, 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_POLYGON);
	    glColor3f(0.4F, 0.4F, 0.4F);
	    glVertex2i(rect.left - 3, y + 1);
	    glColor3f(0.1F, 0.1F, 0.1F);
	    glVertex2i(rect.left, y + 1);
	    glVertex2i(rect.left, y + ((cPopUpLineHeight)) - 4);
	    glColor3f(0.4F, 0.4F, 0.4F);
	    glVertex2i(rect.left - 3, y + ((cPopUpLineHeight)) - 4);
	    glEnd();
	  }
	  if (orthoCGO){
	    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	    CGOColor(orthoCGO, 0.1F, 0.2F, 0.2F);
	    CGOVertex(orthoCGO, rect.right, y + 1, 0.f);
	    CGOColor(orthoCGO, 0.4F, 0.4F, 0.4F);
	    CGOVertex(orthoCGO, rect.right + 3, y + 1, 0.f);
	    CGOColor(orthoCGO, 0.1F, 0.2F, 0.2F);
	    CGOVertex(orthoCGO, rect.right, y + ((cPopUpLineHeight)) - 4, 0.f);
	    CGOColor(orthoCGO, 0.4F, 0.4F, 0.4F);
	    CGOVertex(orthoCGO, rect.right + 3, y + ((cPopUpLineHeight)) - 4, 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_POLYGON);
	    glColor3f(0.1F, 0.2F, 0.2F);
	    glVertex2i(rect.right, y + 1);
	    glColor3f(0.4F, 0.4F, 0.4F);
	    glVertex2i(rect.right + 3, y + 1);
	    glVertex2i(rect.right + 3, y + ((cPopUpLineHeight)) - 4);
	    glColor3f(0.1F, 0.2F, 0.2F);
	    glVertex2i(rect.right, y + ((cPopUpLineHeight)) - 4);
	    glEnd();
	  }
        }

        y -= cPopUpLineHeight;
        if(I->Code[a] == 2)
          y -= 2;
      } else {
	  if (orthoCGO){
	    /* two lines between sections in the menu, one light, one dark */
	    CGOColor(orthoCGO, 0.3F, 0.3F, 0.5F);
	    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	    CGOVertex(orthoCGO, rect.right,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 4, 0.f);
	    CGOVertex(orthoCGO, rect.right,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 3, 0.f);
	    CGOVertex(orthoCGO, rect.left,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 4, 0.f);
	    CGOVertex(orthoCGO, rect.left,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 3, 0.f);
	    CGOEnd(orthoCGO);
	    CGOColor(orthoCGO, 0.6F, 0.6F, 0.8F);
	    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	    CGOVertex(orthoCGO, rect.right,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 5, 0.f);
	    CGOVertex(orthoCGO, rect.right,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 4, 0.f);
	    CGOVertex(orthoCGO, rect.left,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 5, 0.f);
	    CGOVertex(orthoCGO, rect.left,
		      y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 4, 0.f);
	    CGOEnd(orthoCGO);
	  } else {
	    glBegin(GL_LINES);
	    glColor3f(0.3F, 0.3F, 0.5F);
	    glVertex2i(rect.left,
		       y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 3);
	    glVertex2i(rect.right,
		       y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 3);
	    glColor3f(0.6F, 0.6F, 0.8F);
	    glVertex2i(rect.left,
		       y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 4);
	    glVertex2i(rect.right,
		       y + ((cPopUpLineHeight + cPopUpCharMargin) / 2) + 4);
	    glEnd();
	  }
        y -= cPopUpBarHeight;
      }
    }
    if (orthoCGO)
      CGOColorv(orthoCGO, TextColor);
#ifndef PURE_OPENGL_ES_2
    else
      glColor3fv(TextColor);
#endif
    /*    BlockOutline(block); */
  }
}
