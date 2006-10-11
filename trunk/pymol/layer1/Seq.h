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


typedef struct {
  int start;
  int stop;
  int offset; 
  int atom_at; /* starting offset in list */
  int inverse;
  int unaligned;
  int spacer;
  int state;
  int color;
  int tag;
  int is_abbr, hint_no_space;
} CSeqCol;
  
typedef struct {
  int len,ext_len;
  int label_flag,column_label_flag;
  int color;
  char *txt;
  CSeqCol *col,*fill;
  int nCol,cCol;
  int nFill;
  int *char2col;
  int *atom_lists;
  ObjectNameType name; /* associated object */
  struct ObjectMolecule *obj; /* this pointer only valid during update */
  struct AtomInfoType *last_ai;
  int accum,current,title_width; /* temporary stores for aligning */
} CSeqRow;

typedef struct {
  CSeqRow* (*fClick)   (PyMOLGlobals *G,CSeqRow* rowVLA,int button,int row,int col,int mod,int x,int y);
  CSeqRow* (*fDrag)    (PyMOLGlobals *G,CSeqRow* rowVLA,int row,int col,int mod);
  CSeqRow* (*fRelease) (PyMOLGlobals *G,CSeqRow* rowVLA,int button,int row,int col,int mod);
  void (*fRefresh) (PyMOLGlobals *G,CSeqRow* rowVLA);
  int box_active,box_row;
  int box_start_col,box_stop_col;
} CSeqHandler;

int SeqInit(PyMOLGlobals *G);
void SeqFree(PyMOLGlobals *G);
Block *SeqGetBlock(PyMOLGlobals *G);

int SeqIdling(PyMOLGlobals *G);
void SeqInterrupt(PyMOLGlobals *G);
int SeqGetHeight(PyMOLGlobals *G);
void SeqSetHandler(PyMOLGlobals *G,CSeqHandler *handler);
void SeqSetRowVLA(PyMOLGlobals *G,CSeqRow *row,int nRow);
CSeqRow *SeqGetRowVLA(void);
void SeqDirty(PyMOLGlobals *G); /* sequence dirty -- need to update selections */
void SeqChanged(PyMOLGlobals *G); /* sequence changed -- need to rebuild */
void SeqUpdate(PyMOLGlobals *G);

#endif

