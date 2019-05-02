

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

#include <algorithm>

#include"os_predef.h"
#include"os_std.h"
#include"MemoryDebug.h"
#include"Feedback.h"
#include"Ortho.h"

CFeedback::CFeedback(PyMOLGlobals* G, int quiet) : m_G{G}
{
  if(!quiet) {
    for(auto& mask : currentLayer()) {
      mask = FB_Output | FB_Results | FB_Errors | FB_Warnings | FB_Actions |
             FB_Details;
    }

    currentMask(FB_Main) &= ~(FB_Errors); /* suppress opengl errors in main */

  }

  const char *fb_env = getenv("PYMOL_FEEDBACK");
  if(fb_env) {
    int n, sysmod, mask;
    while(sscanf(fb_env, "%i:%i%n", &sysmod, &mask, &n) > 1) {
      setMask(sysmod, mask);
      fb_env += n;
    }
  }
}

/* below we'll presume that any standard feedback on the feedback
module itself will be effected at the Python level, since feedback
levels will be changed as a matter of course inside of PyMOL in order
to quietly perform complex actions.  */

void CFeedback::push()
{
  m_stack.push_back(m_stack.back());
  PRINTFD(m_G, FB_Feedback) " Feedback: push\n" ENDFD;
}

void CFeedback::pop()
{
  if(m_stack.size() > 1) {
    m_stack.pop_back();
  }
  PRINTFD(m_G, FB_Feedback) " Feedback: pop\n" ENDFD;
}

void CFeedback::setMask(unsigned int sysmod, unsigned char mask)
{
  if((sysmod > 0) && (sysmod < FB_Total)) {
    currentMask(sysmod) = mask;
  } else if(!sysmod) {
    std::fill(currentLayer().begin(), currentLayer().end(), mask);
  }
  PRINTFD(m_G, FB_Feedback)
    " FeedbackSetMask: sysmod %d, mask 0x%02X\n", sysmod, mask ENDFD;
}

unsigned char& CFeedback::currentMask(unsigned int sysmod)
{
  return m_stack.back()[sysmod];
}

bool CFeedback::testMask(unsigned int sysmod, unsigned char mask)
{
  return currentMask(sysmod) & mask;
}

void CFeedback::disable(unsigned int sysmod, unsigned char mask)
{
  if((sysmod > 0) && (sysmod < FB_Total)) {
    auto& targetMask = currentMask(sysmod);
    targetMask &= 0xFF - mask;
  } else if(!sysmod) {
    for(auto& obj_mask : currentLayer()) {
      obj_mask &=  0xFF - mask;
    }
  }
  PRINTFD(m_G, FB_Feedback)
    " FeedbackDisable: sysmod %d, mask 0x%02X\n", sysmod, mask ENDFD;

}

void CFeedback::enable(unsigned int sysmod, unsigned char mask)
{
  if((sysmod > 0) && (sysmod < FB_Total)) {
    auto& targetMask = currentMask(sysmod);
    targetMask |= mask;
  } else if(!sysmod) {
    for(auto& obj_mask : currentLayer()) {
      obj_mask |= mask;
    }
  }
  PRINTFD(m_G, FB_Feedback)
    " FeedbackEnable: sysmod %d, mask 0x%02X\n", sysmod, mask ENDFD;

}

void CFeedback::autoAdd(unsigned int sysmod, unsigned char mask, const char *str)
{
  if(testMask(sysmod, mask))
    addColored(str, mask);
}

void CFeedback::add(const char *str)
{
  OrthoAddOutput(m_G, str);
}

void CFeedback::addColored(const char *str, unsigned char mask)
{
  add(str);
}
