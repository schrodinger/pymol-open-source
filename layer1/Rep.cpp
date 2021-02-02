
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
#include"Scene.h"

/*========================================================================*/
/**
 * Delete this instance and return a newly created instance.
 *
 * If creating the new instance fails (e.g. empty rep):
 * - Don't delete this instance
 * - Set cs->Active[rep] to false
 *
 * @todo Use move semantics instead? Return false on failure?
 */
Rep* Rep::rebuild()
{
  assert(cs);
  assert(fNew);
  Rep* tmp = fNew(cs, getState());
  if (tmp) {
    tmp->fNew = fNew;
    delete this;
    return tmp;
  }

  // nothing returned -- visibility is zero...
  // keep the old object around, but inactive (TODO why?)
  cs->Active[type()] = false;
  return this;
}

/*========================================================================*/
/**
 * Rebuild if necessary (according to invalidation status). Returns either this
 * instance, or deletes this instance and returns a new one, or NULL if the rep
 * became empty/inactive.
 */
Rep* Rep::update()
{
  assert(cs);

  if (MaxInvalid == cRepInvNone) {
    return this;
  }

  auto I = this;
  auto const rep = type();
  auto const* cs_ = cs;

  assert(cs_->Active[rep]);

  if (MaxInvalid == cRepInvPick) {
    switch (rep) {
    case cRepLine:
    case cRepCyl:
    case cRepRibbon:
    case cRepNonbonded:
      // TODO Is this needed? How about other pickable reps, like:
      // - cartoon
      // - spheres
      // - surface with pick_surface=on
      MaxInvalid = cRepInvRep;
    }
  }

  if (MaxInvalid < cRepInvColor) {
    // nothing to do
  } else if (MaxInvalid == cRepInvColor) {
    I = recolor();
  } else if (MaxInvalid > cRepInvVisib || !sameVis()) {
    I = rebuild();
  } else if (!sameColor()) {
    I = recolor();
  }

  if (!cs_->Active[rep]) {
    delete I;
    return nullptr;
  }

  if (I) {
    I->MaxInvalid = cRepInvNone;
  }

  return I;
}

/*========================================================================*/
/**
 * Request that the rep gets updated. Update happens on next scene redraw.
 *
 * @param level Invalidate to at least this level
 */
void Rep::invalidate(cRepInv_t level)
{
  auto I = this;
  SceneInvalidatePicking(I->G); // for now, if anything invalidated, then invalidate picking
  if(level > I->MaxInvalid)
    I->MaxInvalid = level;
}

/**
 * Get the visRep mask according to auto_show_* settings
 */
cRepBitmask_t RepGetAutoShowMask(PyMOLGlobals * G)
{
  cRepBitmask_t mask = 0;
  if (SettingGetGlobal_b(G, cSetting_auto_show_lines))     mask |= cRepLineBit;
  if (SettingGetGlobal_b(G, cSetting_auto_show_spheres))   mask |= cRepSphereBit;
  if (SettingGetGlobal_b(G, cSetting_auto_show_nonbonded)) mask |= cRepNonbondedBit;
  return mask;
}

/*========================================================================*/
/**
 * Derived classes should overrride this method.
 *
 * TODO: What is this default implementation useful for? Debugging?
 */
void Rep::render(RenderInfo* info)
{
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
Rep::Rep(pymol::CObject* obj_, int state)
    : G(obj_->G)
    , obj(obj_)
{
  context.object = obj_;
  context.state = state;
}

Rep::Rep(CoordSet* cs_, int state)
    : Rep(cs_->Obj, state)
{
  cs = cs_;
}

Rep::~Rep()
{
  FreeP(P);
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
