
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

#include"main.h"
#include"Feedback.h"
#include"Rep.h"
#include"MemoryDebug.h"
#include"CoordSet.h"
#include"P.h"
#include"Util.h"


/*========================================================================*/
static
struct Rep *RepRebuild(struct Rep *I, struct CoordSet *cs, int state, int rep)
{
  Rep *tmp = NULL;

  PRINTFD(I->G, FB_Rep)
    " RepRebuild-Debug: entered: rep %d I->fNew %p\n", rep, (void *) I->fNew ENDFD;

  if(I->fNew) {
    tmp = I->fNew(cs, state);
    if(tmp) {
      tmp->fNew = I->fNew;
      I->fFree(I);
    } else {                    /* nothing returned -- visibility is zero... */
      cs->Active[rep] = false;  /* keep the old object around, but inactive */
      tmp = I;
    }
  } else
    I->fFree(I);
  return (tmp);
}


/*========================================================================*/
static
struct Rep *RepUpdate(struct Rep *I, struct CoordSet *cs, int state, int rep)
{

  PRINTFD(I->G, FB_Rep)
    " RepUpdate-Debug: entered: rep %d I->MaxInvalid %d\n", rep, I->MaxInvalid ENDFD;

  if(I->MaxInvalid) {
    if(I->MaxInvalid == cRepInvPick) {
      if((rep == cRepLine) ||
         (rep == cRepCyl) || (rep == cRepRibbon) || (rep == cRepNonbonded))
        I->MaxInvalid = cRepInvRep;
    }

    if(I->MaxInvalid < cRepInvColor) {
    } else if(I->MaxInvalid == cRepInvColor) {
      if(I->fRecolor) {
        I->fRecolor(I, cs);
      } else {
        I = I->fRebuild(I, cs, state, rep);
      }
    } else if(I->MaxInvalid <= cRepInvVisib) {
      int rebuilt = false;
      if(I->fSameVis) {
        if(!I->fSameVis(I, cs)){
          I = I->fRebuild(I, cs, state, rep);
          rebuilt = true;
        }
      }
      if(I->fSameColor) {
        if (!rebuilt){
          if(!I->fSameColor(I, cs)){
            I->fRecolor(I, cs);
          }
        }
      }
      if (!I->fSameVis && !I->fSameColor)
        I = I->fRebuild(I, cs, state, rep);
    } else if(I->MaxInvalid >= cRepInvCoord) {
      I = I->fRebuild(I, cs, state, rep);
      if(!cs->Active[rep]) {
        I->fFree(I);
        I = NULL;
      }
      /*      if(I->fNew) {
         tmp = I->fNew(cs);
         if(I->fFree) I->fFree(I);
         I=tmp;
       */
    } else {
      I = I->fRebuild(I, cs, state, rep);
    }
    if(I)
      I->MaxInvalid = 0;
  }
  return (I);
}


/*========================================================================*/
void RepInvalidate(struct Rep *I, struct CoordSet *cs, int level)
{
  if(level > I->MaxInvalid)
    I->MaxInvalid = level;
}

/*
 * Get the visRep mask according to auto_show_* settings
 */
int RepGetAutoShowMask(PyMOLGlobals * G)
{
  int mask = 0;
  if (SettingGetGlobal_b(G, cSetting_auto_show_lines))     mask |= cRepLineBit;
  if (SettingGetGlobal_b(G, cSetting_auto_show_spheres))   mask |= cRepSphereBit;
  if (SettingGetGlobal_b(G, cSetting_auto_show_nonbonded)) mask |= cRepNonbondedBit;
  return mask;
}

/*========================================================================*/
static void RepRenderBox(struct Rep *this_, RenderInfo * info)
{
  PyMOLGlobals *G = this_->G;
  if(G->HaveGUI && G->ValidContext) {
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
    glBegin(GL_LINE_LOOP);
    glVertex3f(-0.5F, -0.5F, -0.5F);
    glVertex3f(-0.5F, -0.5F, 0.5F);
    glVertex3f(-0.5F, 0.5F, 0.5F);
    glVertex3f(-0.5F, 0.5F, -0.5F);

    glVertex3f(0.5F, 0.5F, -0.5F);
    glVertex3f(0.5F, 0.5F, 0.5F);
    glVertex3f(0.5F, -0.5F, 0.5F);
    glVertex3f(0.5F, -0.5F, -0.5F);
    glEnd();

    glBegin(GL_LINES);
    glVertex3i(0, 0, 0);
    glVertex3i(1, 0, 0);

    glVertex3i(0, 0, 0);
    glVertex3i(0, 2, 0);

    glVertex3i(0, 0, 0);
    glVertex3i(0, 0, 3);

    glEnd();
#endif
  }

}


/*========================================================================*/
void RepInit(PyMOLGlobals * G, Rep * I)
{
  UtilZeroMem(I, sizeof(Rep));
  I->G = G;
  I->fInvalidate = RepInvalidate;
  I->fUpdate = RepUpdate;
  I->fRender = RepRenderBox;
  I->fRebuild = RepRebuild;
}


/*========================================================================*/
void RepPurge(Rep * I)
{
  FreeP(I->P);
}

RepIterator::RepIterator(PyMOLGlobals * G, int rep_) {
  if (rep_ < 0){
    end = cRepCnt;
    rep = -1;
  } else {
    end = rep_ + 1;
    rep = rep_ - 1;
  }
}
