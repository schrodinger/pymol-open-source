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
  int spacer;
  int state;
} CSeqCol;
  
typedef struct {
  int len,ext_len;
  int label_flag;
  int color;
  char *txt;
  CSeqCol *col;
  int nCol;
  int *char2col;
  int *atom_lists;
  char name[ObjNameMax]; /* associated object */
  struct ObjectMolecule *obj; /* this pointer only valid during update */
  int accum,current,title_width; /* temporary stores for aligning */

} CSeqRow;

typedef struct {
  CSeqRow* (*fClick)   (CSeqRow* rowVLA,int button,int row,int col,int mod,int x,int y);
  CSeqRow* (*fDrag)    (CSeqRow* rowVLA,int row,int col,int mod);
  CSeqRow* (*fRelease) (CSeqRow* rowVLA,int button,int row,int col,int mod);
  void (*fRefresh) (CSeqRow* rowVLA);
  int box_active,box_row;
  int box_start_col,box_stop_col;
} CSeqHandler;

void SeqInit(void);
void SeqFree(void);
Block *SeqGetBlock(void);

void SeqInit(void);
void SeqFree(void);
Block *SeqGetBlock(void);
int SeqIdling(void);
void SeqInterrupt(void);
int SeqGetHeight(void);
void SeqSetHandler(CSeqHandler *handler);
void SeqSetRowVLA(CSeqRow *row,int nRow);
CSeqRow *SeqGetRowVLA(void);
void SeqDirty(void); /* sequence dirty -- need to update selections */
void SeqChanged(void); /* sequence changed -- need to rebuild */
void SeqUpdate(void);

#endif

