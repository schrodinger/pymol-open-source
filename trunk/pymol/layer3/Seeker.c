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
#include"P.h"
#include"Selector.h"
#include"Wizard.h"

#define cTempSeekerSele "_seeker"

typedef struct {
  int a;
} CSeeker;

CSeeker Seeker;

CSeqHandler SeekerHandler;


static void SeekerSelectionToggle(CSeqRow* rowVLA,int row_num,int col_num,int inc_or_excl)
{
  char selName[ObjNameMax];
  OrthoLineType buf1,buf2;
  char *buf_vla = NULL;
  int buf_size = 0;
  buf_vla = VLAlloc(char,1000);

  ExecutiveGetActiveSeleName(selName,true);
  
  {
    CSeqRow *row;
    CSeqCol *col;
    int *atom_list;
    char prefix[3]="";
    int logging = SettingGet(cSetting_logging);

    if(logging==cPLog_pml)
      strcpy(prefix,"_ ");
    row = rowVLA + row_num;
    col = row->col + col_num;

    if( ExecutiveFindObjectByName(row->name)) {
      atom_list = row->atom_lists + col->atom_at;

      /* build up a selection consisting of residue atoms */

      UtilConcatVLA(&buf_vla,&buf_size,"none");
      while((*atom_list)>=0) {
        sprintf(buf1,"|%s`%d",
                row->name,*atom_list+1);
        UtilConcatVLA(&buf_vla,&buf_size,buf1);
        atom_list++;
      }
      
      SelectorCreate(cTempSeekerSele,buf_vla,NULL,true,NULL);      
      if(logging) SelectorLogSele(cTempSeekerSele);
      
      /* selection or deselecting? */

      if(col->inverse) {
        sprintf(buf1,"((%s) or (%s))",
                selName,cTempSeekerSele);
      } else {
        sprintf(buf1,"((%s) and not (%s))",
                selName,cTempSeekerSele);
      }
      
      /* create the new active selection */

      SelectorCreate(selName,buf1,NULL,false,NULL);
      {
        sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,selName,buf1);
        PLog(buf2,cPLog_no_flush);
      }

      ExecutiveDelete(cTempSeekerSele);
      if(logging) {
        sprintf(buf2,"%scmd.delete(\"%s\")\n",prefix,cTempSeekerSele);
        PLog(buf2,cPLog_no_flush);
        PLogFlush();
      }
      
      if(SettingGet(cSetting_auto_show_selections))
        ExecutiveSetObjVisib(selName,1);
      WizardDoSelect(selName);
    }
  }
  VLAFreeP(buf_vla);
}

static CSeqRow* SeekerClick(CSeqRow* rowVLA,int button,int row_num,int col_num,int mod)
{
  CSeqRow *row;
  CSeqCol *col;
  
  row = rowVLA + row_num;
  col = row->col + col_num;
  
  if(col->inverse) {
    SeekerSelectionToggle(rowVLA,row_num,col_num,true);
  } else {
    SeekerSelectionToggle(rowVLA,row_num,col_num,false);
  }

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
  int label_flag = true;
  CSeqRow *row_vla,*row,*lab=NULL;
  row_vla = VLACalloc(CSeqRow,10);

  /* FIRST PASS: get all the residues represented properly */

  while(ExecutiveIterateObjectMolecule(&obj,&hidden)) {
    if(obj->Obj.Enabled) {
      int a;
      AtomInfoType *last = NULL,*last_segi=NULL,*last_chain = NULL;
      int last_abbr = false;
      int nCol = 0;
      int nListEntries = 1; /* first list starts at 1 always... */
      int est_col = obj->NAtom/5+1;
      int est_char = obj->NAtom*4;
      int first_atom_in_label;
      int min_pad = -1;
      CSeqCol *r1 = NULL,*l1=NULL;/* *col */

      /* allocate a row for labels, if present
         the text for the labels and the residues will line up exactly 
      */

      VLACheck(row_vla,CSeqRow,nRow);
      if(label_flag)
        {
          lab = row_vla + nRow++;
          lab->txt = VLAlloc(char,est_char);
          lab->col = VLACalloc(CSeqCol,est_col);
          lab->label_flag = true;
        }
      else
        lab = NULL;

      VLACheck(row_vla,CSeqRow,nRow);

      row = row_vla+nRow;
      row->txt = VLAlloc(char,est_char);
      row->col = VLACalloc(CSeqCol,est_col);
      row->atom_lists = VLACalloc(int,obj->NAtom+est_col);
      row->char2col = VLACalloc(int,est_char);
      row->obj = obj;
      strcpy(row->name,obj->Obj.Name);
      row->color = obj->Obj.Color;
      ai = obj->AtomInfo;

      /* copy object name onto label row */

      if(lab) {
        
        int st_len;
        /* copy label text */

        VLACheck(lab->col,CSeqCol,nCol);
        l1 = lab->col + nCol;
        l1->start = lab->len;
        strcpy(lab->txt+lab->len,"/"); lab->len++;
        strcpy(lab->txt+lab->len,obj->Obj.Name); lab->len+=strlen(obj->Obj.Name);
        strcpy(lab->txt+lab->len,"/"); lab->len++;
        strcpy(lab->txt+lab->len,ai->segi); lab->len+=strlen(ai->segi);
        strcpy(lab->txt+lab->len,"/"); lab->len++;
        strcpy(lab->txt+lab->len,ai->chain); lab->len+=strlen(ai->chain);
        strcpy(lab->txt+lab->len,"/"); lab->len++;
        l1->stop = lab->len;
        st_len = l1->stop - l1->start;

        last_segi = ai;
        last_chain = ai;
        /* blank equivalent text for sequence row */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          row->len += st_len;
          r1->stop = row->len;
          memset(row->txt+r1->start,32,st_len);

        nCol++;
      }

      for(a=0;a<obj->NAtom;a++) {
        first_atom_in_label = false;
        if(lab&&!AtomInfoSameSegmentP(last_segi,ai)) {

          int st_len;

          if(row->len<min_pad) {
            row->len = min_pad;
          }
          min_pad = -1;

          /* copy label text */
          
          VLACheck(lab->col,CSeqCol,nCol);
          l1 = lab->col + nCol;
          l1->start = lab->len;
          strcpy(lab->txt+lab->len,ai->segi); lab->len+=strlen(ai->segi);
          strcpy(lab->txt+lab->len,"/"); lab->len++;
          strcpy(lab->txt+lab->len,ai->chain); lab->len+=strlen(ai->chain);
          strcpy(lab->txt+lab->len,"/"); lab->len++;
          l1->stop = lab->len;
          st_len = l1->stop - l1->start;

          /* blank equivalent text for sequence row */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          row->len += st_len;
          r1->stop = row->len;
          memset(row->txt+r1->start,32,st_len);
          
          nCol++;
          
          last_abbr = false;

          last_segi = ai;
          last_chain = ai;

        } else if(lab&&!AtomInfoSameChainP(last_chain,ai)) {

          int st_len;

          if(row->len<min_pad) {
            row->len = min_pad;
          }
          min_pad = -1;

          /* copy label text */
          
          VLACheck(lab->col,CSeqCol,nCol);
          l1 = lab->col + nCol;
          l1->start = lab->len;
          strcpy(lab->txt+lab->len,ai->chain); lab->len+=strlen(ai->chain);
          strcpy(lab->txt+lab->len,"/"); lab->len++;
          l1->stop = lab->len;
          st_len = l1->stop - l1->start;
          
          /* blank equivalent text for sequence row */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          row->len += st_len;
          r1->stop = row->len;
          memset(row->txt+r1->start,32,st_len);
          
          nCol++;
                    
          last_abbr = false;
          last_chain = ai;
        }

        if(min_pad<0)
          min_pad = strlen(ai->resi) + row->len + 1;
        
        switch(0) {
        case 0:
          if(!AtomInfoSameResidueP(last,ai)) {
            char abbr[2] = "1";
            last = ai;            

            VLACheck(row->col,CSeqCol,nCol);
            r1 = row->col+nCol;
            r1->start = row->len;
            first_atom_in_label = true;

            abbr[0] = SeekerGetAbbr(ai->resn);
            
            if(!abbr[0]) {
              if(last_abbr) {
                UtilConcatVLA(&row->txt,&row->len," ");                
                r1->start = row->len;
              }
              
            
              UtilConcatVLA(&row->txt,&row->len,ai->resn);
              r1->stop = row->len;

              UtilConcatVLA(&row->txt,&row->len," ");
            } else {
              UtilConcatVLA(&row->txt,&row->len,abbr);
              r1->stop = row->len;
            }
            nCol++;
            last_abbr=abbr[0];
          }
          
          break;
        case 1:
          if(!AtomInfoSameResidueP(last,ai)) {
            last = ai;

            VLACheck(row->col,CSeqCol,nCol);
            r1 = row->col+nCol;
            r1->start = row->len;
            first_atom_in_label = true;

            UtilConcatVLA(&row->txt,&row->len,ai->resn);
            r1->stop = row->len;
            UtilConcatVLA(&row->txt,&row->len," ");
            nCol++;
          }
          break;
        case 2:
          UtilConcatVLA(&row->txt,&row->len,ai->name);
          UtilConcatVLA(&row->txt,&row->len," ");
          break;
        }

        if(first_atom_in_label) {
          if(nCol>1)  { /* terminate current list, if any */
            VLACheck(row->atom_lists,int,nListEntries);
            row->atom_lists[nListEntries]=-1;
            nListEntries++;
          }
          if(r1) {
            r1->atom_at = nListEntries;
            {
              int c;
              VLACheck(row->char2col,int,row->len);
              for(c=r1->start;c<r1->stop;c++) 
                row->char2col[c]=nCol; /* implicitly offset by 1 for convenience! */
            }
          }
        }
        
        VLACheck(row->atom_lists,int,nListEntries);
        row->atom_lists[nListEntries] = a;
        nListEntries++;
        ai++;
      }

      if(lab) {
        lab->txt = VLASetSize(lab->txt,row->len+1);
        lab->txt[row->len] = 0;

        VLACheck(lab->col,CSeqCol,nCol); /* make sure we've got column records for labels too */
        lab->nCol = nCol;
      }

      row->txt = VLASetSize(row->txt,row->len+1);
      row->txt[row->len] = 0;

      row->nCol = nCol;

      /* terminate last atom list */
      VLACheck(row->atom_lists,int,nListEntries);
      row->atom_lists[nListEntries] = -1;
      nListEntries++;
      nRow++;
    }
  }

  /* SECOND PASS: align columns to reflect current alignment and fixed labels */
  {
    int a,b;
    int nCol;

    for(a=0;a<nRow;a++) {
      row = row_vla + a;
      nCol = row->nCol;
      if(row->label_flag)
        lab=row;
      else {
        for(b=0;b<nCol;b++) {
          CSeqCol *r1 = row->col + b,*l1=NULL;
          r1->offset = r1->start;
          if(lab) {
            l1 = lab->col + b;
            if(l1->stop) 
              l1->offset = r1->offset;
          }
        }
        lab = NULL;
      }
    }
  }

  /* THIRD PASS: fill in labels, based on actual residue spacing */

  {
    int a,b,c;
    int nCol;

    for(a=0;a<nRow;a++) {
      lab = row_vla + a;
      if(lab->label_flag) {
        int next_open = 0;
        int *atom_list;
        int st_len;
        ObjectMolecule *obj;
        AtomInfoType *ai;
        row = lab+1;
        nCol = row->nCol;
        obj = row->obj;

        for(b=0;b<nCol;b++) {
          CSeqCol *r1 = row->col + b;
          CSeqCol *l1 = lab->col + b;


          if(l1->stop) {/* if label is already present, just line it up */
            l1->offset = r1->offset;
          } else if((r1->offset >= next_open) && (r1->atom_at)) {
            atom_list = row->atom_lists + r1->atom_at;
            
            ai = obj->AtomInfo + (*atom_list); /* get first atom in list */
            st_len = strlen(ai->resi);
            VLACheck(lab->txt,char,lab->len+st_len+1);
            strcpy(lab->txt+lab->len, ai->resi); /* copy the residue identifier */
            l1->start = lab->len;
            lab->len += st_len;
            l1->stop = lab->len;
            l1->offset = r1->offset;
            next_open = r1->offset + st_len + 1;

            /* make sure this label doesn't conflict with a fixed label */

            for(c=b+1;c<nCol;c++) { 
              CSeqCol *l2 = lab->col + c;   
              if(l2->offset&&(l2->offset<next_open)) {
                l1->start=0;
                l1->stop=0;
                break;
              }
              if((c-b)>st_len) /* only search as many columns as characters */
                break;
            }
          }
        }
      }
    }
  }
  SeqSetRowVLA(row_vla,nRow);
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


