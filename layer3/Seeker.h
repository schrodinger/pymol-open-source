
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
#ifndef _H_Seeker
#define _H_Seeker

#include <vector>
#include"Ortho.h"
#include"ObjectMolecule.h"
#include "Seq.h"

#define cTempSeekerSele "_seeker"
#define cTempCenterSele "_seeker_center"
#define cTempSeekerSele2 "_seeker2"


int SeekerInit(PyMOLGlobals * G);
void SeekerFree(PyMOLGlobals * G);
void SeekerUpdate(PyMOLGlobals * G);
char SeekerGetAbbr(PyMOLGlobals * G, const char *abbr, char water, char unknown);

namespace GapMode{
enum {
    NONE      = 0,
    ALL       = 1,
    SINGLE    = 2
};
}//namespace GapMode

struct SeekerDragInfo {
  int start_col;
  int last_col;
  int row;
  int dir;
  int start_toggle;
  int setting;
  int button;
};

void SeekerSetDragInfo(PyMOLGlobals* G, const SeekerDragInfo& dragInfo);
SeekerDragInfo SeekerGetDragInfo(PyMOLGlobals* G);
void SeekerSelectionCenter(PyMOLGlobals * G, int action);
void SeekerSelectionUpdateCenter(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row_num,
                                 int col_num, int start_over);

#endif
