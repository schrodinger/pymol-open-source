
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
#ifndef _H_Seq
#define _H_Seq

#include "Ortho.h"
#include "PyMOLObject.h"
#include <vector>

struct CSeqCol {
  int start;
  int stop;
  int offset;
  int atom_at;                  /* starting offset in list */
  int inverse;
  int unaligned;
  int spacer;
  int state;
  int color;
  int tag;
  int is_abbr, hint_no_space;
};

struct ObjectMolecule;
struct AtomInfoType;

struct CSeqRow {
  ov_size len, ext_len;
  int label_flag, column_label_flag;
  int color;
  pymol::vla<char> txt;
  pymol::vla<CSeqCol> col;
  pymol::vla<CSeqCol> fill;
  int nCol, cCol;
  int nFill;
  pymol::vla<int> char2col;
  pymol::vla<int> atom_lists;
  ObjectNameType name;          /* associated object */
  ObjectMolecule *obj;   /* this pointer only valid during update */
  AtomInfoType *last_ai;
  int accum, current, title_width;      /* temporary stores for aligning */
};

struct CSeqHandler {
  virtual CSeqRow* click(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button, int row, int col,
                      int mod, int x, int y) = 0;
  virtual CSeqRow* drag(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int row, int col, int mod) = 0;
  virtual CSeqRow* release(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA, int button, int row, int col,
                        int mod) = 0;
 virtual void refresh(PyMOLGlobals * G, std::vector<CSeqRow>& rowVLA) = 0;

  bool box_active{};
  int box_row{};
  int box_start_col{};
  int box_stop_col{};
};

int SeqInit(PyMOLGlobals * G);
void SeqFree(PyMOLGlobals * G);
Block *SeqGetBlock(PyMOLGlobals * G);

int SeqGetHeight(PyMOLGlobals * G);
void SeqSetHandler(PyMOLGlobals * G, CSeqHandler * handler);
void SeqSetRow(PyMOLGlobals * G, std::vector<CSeqRow>&& row, int nRow);
void SeqDirty(PyMOLGlobals * G);        /* sequence dirty -- need to update selections */
void SeqChanged(PyMOLGlobals * G);      /* sequence changed -- need to rebuild */
void SeqUpdate(PyMOLGlobals * G);

#endif
