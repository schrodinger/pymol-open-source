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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"
#include"Err.h"
#include"Util.h"
#include"Seq.h"
#include"Seeker.h"
#include"MemoryDebug.h"
#include"Executive.h"

typedef struct {
} CSeeker;

CSeeker Seeker;

CSeqHandler SeekerHandler;

static CSeqRow* SeekerClick(CSeqRow* rowVLA,int button,int row,int col,int mod)
{
  printf("%d %d %d %d\n",button,row,col,mod);
  return NULL;
}

static CSeqRow* SeekerDrag(CSeqRow* rowVLA,int row,int col,int mod)
{
  return NULL;
}

static CSeqRow* SeekerRelease(CSeqRow* rowVLA,int button,
                              int row,int col,int mod)
{
  return NULL;
}

static void SeekerReset(void)
{
}

static void SeekerUpdateSimple(void)
{
}

static char SeekerGetAbbr(char *abbr)
{
  
  switch(abbr[0]) {
  case 'A':
    switch(abbr[1]) {
    case 'L': 
      if(abbr[2]=='A')
        return 'A';
      break;
    case 'R': 
      if(abbr[2]=='G')
        return 'R';
      break;
    case 'S': 
      switch(abbr[2]) {
      case 'P':
        return 'D';
        break;
      case 'N':
        return 'N';
        break;
      }
      break;
    }
    break;
  case 'C':
    switch(abbr[1]) {
    case 'Y': 
      switch(abbr[2]) {
      case 'S':
      case 'X':
        return 'C';
        break;
      }
      break;
    }
    break;
  case 'G':
    switch(abbr[1]) {
    case 'L': 
      switch(abbr[2]) {
      case 'N':
        return 'Q';
        break;
      case 'U':
        return 'D';
        break;
      case 'Y':
        return 'G';
        break;
      }
    }
    break;
  case 'H':
    switch(abbr[1]) {
    case 'I': 
      switch(abbr[2]) {
      case 'S':
      case 'D':
      case 'E':
        return 'H';
        break;
      }
      break;
    case 'O': 
      switch(abbr[2]) {
      case 'H':
        return 'O';
        break;
      }
      break;

    }
  case 'I':
    switch(abbr[1]) {
    case 'L': 
      switch(abbr[2]) {
      case 'E':
        return 'I';
        break;
      }
    }
    break;
  case 'L':
    switch(abbr[1]) {
    case 'E': 
      switch(abbr[2]) {
      case 'U':
        return 'L';
        break;
      }
      break;
    case 'Y': 
      switch(abbr[2]) {
      case 'S':
        return 'K';
        break;
      }
      break;
    }
    break;
  case 'M':
    switch(abbr[1]) {
    case 'E': 
      switch(abbr[2]) {
      case 'T':
        return 'M';
        break;
      }
    }
    break;
  case 'P':
    switch(abbr[1]) {
    case 'H':
      switch(abbr[2]) {
      case 'E':
        return 'F';
        break;
      }
      break;     
    case 'R': 
      switch(abbr[2]) {
      case 'O':
        return 'P';
        break;
      }
      break;
    }
    break;
  case 'S':
    switch(abbr[1]) {
    case 'E': 
      switch(abbr[2]) {
      case 'R':
        return 'S';
        break;
      }
      break;
    }
    break;
  case 'T':
    switch(abbr[1]) {
    case 'H': 
      switch(abbr[2]) {
      case 'R':
        return 'R';
        break;
      }
      break;
    case 'R': 
      switch(abbr[2]) {
      case 'P':
        return 'W';
        break;
      }
      break;
    case 'Y': 
      switch(abbr[2]) {
      case 'R':
        return 'Y';
        break;
      }
      break;
    }
    break;
  case 'V':
    switch(abbr[1]) {
    case 'A': 
      switch(abbr[2]) {
      case 'L':
        return 'V';
        break;
      }
      break;
    }
    break;
  }
  return 0;
}

void SeekerUpdate(void)
{
  /*  CObject *o = NULL;
      int s;*/

  void *hidden = NULL;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  int nRow = 0;
  CSeqRow *row,*r1;
  row = VLACalloc(CSeqRow,10);
  while(ExecutiveIterateObjectMolecule(&obj,&hidden)) {
    if(obj->Obj.Enabled) {
      int a;
      int st_len = 0;
      AtomInfoType *last = NULL;
      int last_abbr = false;
      int nCol = 0;
      int nListEntries = 1; /* first list starts at 1 always... */
      int est_col = obj->NAtom/5+1;
      int est_char = obj->NAtom*4;
      int first_atom_in_label;
      VLACheck(row,CSeqRow,nRow);
      CSeqCol *c1 = NULL;/* *col */

      r1 = row+nRow;
      r1->txt = VLAlloc(char,est_char);
      r1->col = VLACalloc(CSeqCol,est_col);
      r1->atom_lists = VLACalloc(int,obj->NAtom+est_col);
      r1->char2col = VLACalloc(int,est_char);
      
      ai = obj->AtomInfo;
      for(a=0;a<obj->NAtom;a++) {
        first_atom_in_label = false;
        switch(0) {
        case 0:
          if(!AtomInfoSameResidueP(last,ai)) {
            last = ai;            
            char abbr[2] = "1";

            VLACheck(r1->col,CSeqCol,nCol);
            c1 = r1->col+nCol;
            c1->start = st_len;
            first_atom_in_label = true;

            abbr[0] = SeekerGetAbbr(ai->resn);
            
            if(!abbr[0]) {
              if(last_abbr)
                UtilConcatVLA(&r1->txt,&st_len," ");
            
              UtilConcatVLA(&r1->txt,&st_len,ai->resn);
              c1->stop = st_len;

              UtilConcatVLA(&r1->txt,&st_len," ");
            } else {
              UtilConcatVLA(&r1->txt,&st_len,abbr);
              c1->stop = st_len;
            }
            if(nCol==3)
              c1->inverse = true;
            nCol++;
            last_abbr=abbr[0];
          }
          
          break;
        case 1:
          if(!AtomInfoSameResidueP(last,ai)) {
            last = ai;

            VLACheck(r1->col,CSeqCol,nCol);
            c1 = r1->col+nCol;
            c1->start = st_len;
            first_atom_in_label = true;

            UtilConcatVLA(&r1->txt,&st_len,ai->resn);
            c1->stop = st_len;
            UtilConcatVLA(&r1->txt,&st_len," ");
            nCol++;
          }
          break;
        case 2:
          UtilConcatVLA(&r1->txt,&st_len,ai->name);
          UtilConcatVLA(&r1->txt,&st_len," ");
          break;
        }

        if(first_atom_in_label) {
          if(nCol>1)  { /* terminate current list, if any */
            VLACheck(r1->atom_lists,int,nListEntries);
            r1->atom_lists[nListEntries]=-1;
            nListEntries++;
          }
          if(c1) {
            c1->atom_at = nListEntries;
            {
              int c;
              VLACheck(r1->char2col,int,st_len);
              for(c=c1->start;c<c1->stop;c++) 
                r1->char2col[c]=nCol; /* implicitly offset by 1 for convenience! */
            }
          }
        }
        
        VLACheck(r1->atom_lists,int,nListEntries);
        r1->atom_lists[nListEntries] = a;
        nListEntries++;
        ai++;
      }

      r1->txt = VLASetSize(r1->txt,st_len+1);
      r1->len = st_len;
      r1->nCol = nCol;

      /* terminate last atom list */
      VLACheck(r1->atom_lists,int,nListEntries);
      r1->atom_lists[nListEntries] = -1;
      nListEntries++;
      nRow++;
    }
  }
  SeqSetRowVLA(row,nRow);
  SeekerHandler.fClick = SeekerClick;
  SeekerHandler.fRelease = SeekerRelease;
  SeekerHandler.fDrag = SeekerDrag;
  SeqSetHandler(&SeekerHandler);
}

static void SeekerInvalidate(void)
{
}

void SeekerInit(void)
{
  /*  CSeeker *I = &Seeker;  */
}

void SeekerFree(void)
{
}


