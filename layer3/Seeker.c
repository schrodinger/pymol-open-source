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
#include"Scene.h"
#include"Menu.h"

#define cTempSeekerSele "_seeker"
#define cTempCenterSele "_seeker_center"

typedef struct {
  CSeqHandler handler; /* must be first */
  int drag_start_col, drag_last_col;
  int drag_row;
  int dragging, drag_setting;
  int drag_button;
} CSeeker;

CSeeker Seeker;


static void BuildSeleFromAtomList(char *obj_name,int *atom_list,char *sele_name,int start_fresh)
{
  char *buf_vla = NULL;
  int buf_size = 0;
  OrthoLineType buf1;
  buf_vla = VLAlloc(char,4000);

  while((*atom_list)>=0) {
    if(start_fresh) {
      UtilConcatVLA(&buf_vla,&buf_size,"none");
    } else {
      UtilConcatVLA(&buf_vla,&buf_size,"?");
      UtilConcatVLA(&buf_vla,&buf_size,sele_name);
    }
    
    while((*atom_list)>=0) {
      sprintf(buf1,"|%s`%d",obj_name,*atom_list+1);
      UtilConcatVLA(&buf_vla,&buf_size,buf1);
      atom_list++;
      if((buf_size>(sizeof(OrthoLineType)+ObjNameMax+50))) {
        SelectorCreate(sele_name,buf_vla,NULL,true,NULL);      
        buf_size = 0;
        start_fresh=false;
        break;
      }
    }
    
    if(buf_size) {
      SelectorCreate(sele_name,buf_vla,NULL,true,NULL);            
    }
  }
  VLAFreeP(buf_vla);
}


static void SeekerSelectionToggle(CSeqRow* rowVLA,int row_num,int col_num,int inc_or_excl)
{
  char selName[ObjNameMax];
  OrthoLineType buf1,buf2;

  
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

      BuildSeleFromAtomList(row->name,atom_list,cTempSeekerSele,true);
      if(logging) SelectorLogSele(cTempSeekerSele);
      
      if(!WizardDoSelect(cTempSeekerSele)) {
        
        ExecutiveGetActiveSeleName(selName,true);
        
        /* selection or deselecting? */
        
        if(inc_or_excl) {
          if(!col->spacer) {
            col->inverse = true;
            sprintf(buf1,"((%s) or (?%s))",
                    selName,cTempSeekerSele);
          }
        } else {
          if(!col->spacer) {
            col->inverse = false;
            sprintf(buf1,"((%s) and not (?%s))",
                    selName,cTempSeekerSele);
          }
        }
        
        /* create the new active selection */
        
        SelectorCreate(selName,buf1,NULL,true,NULL);
        {
          sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,selName,buf1);
          PLog(buf2,cPLog_no_flush);
        }
      }
      
      ExecutiveDelete(cTempSeekerSele);
      if(logging) {
        sprintf(buf2,"%scmd.delete(\"%s\")\n",prefix,cTempSeekerSele);
        PLog(buf2,cPLog_no_flush);
        PLogFlush();
      }
      
      if(SettingGet(cSetting_auto_show_selections))
        ExecutiveSetObjVisib(selName,1);
      SceneDirty();
    }
  }
}


static void SeekerSelectionCenter(CSeqRow* rowVLA,int row_num,int col_num,int start_over)
{
  
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
      
      BuildSeleFromAtomList(row->name,atom_list,cTempCenterSele,start_over);
      if(logging) SelectorLogSele(cTempCenterSele);
      
    }
  }
}

static void SeekerSelectionCenterGo(int action)
{
  OrthoLineType buf2;
    char prefix[3]="";
  int logging = SettingGet(cSetting_logging);
    if(logging==cPLog_pml)
      strcpy(prefix,"_ ");

    switch(action) {
    case 0: /* center cumulative*/
      ExecutiveCenter(cTempCenterSele,-1,true);
      if(logging) {
        sprintf(buf2,"%scmd.center(\"%s\")\n",prefix,cTempCenterSele);
        PLog(buf2,cPLog_no_flush);
        PLogFlush();
      }
      break;
    case 1: /* zoom */
      ExecutiveWindowZoom(cTempCenterSele,0.0,-1,false);
      if(logging) {
        sprintf(buf2,"%scmd.zoom(\"%s\")\n",prefix,cTempCenterSele);
        PLog(buf2,cPLog_no_flush);
        PLogFlush();
      }
      break;
    }
}

static CSeqRow* SeekerClick(CSeqRow* rowVLA,int button,int row_num,int col_num,int mod,int x,int y)
{
  CSeqRow *row;
  CSeqCol *col;
  /*  char selName[ObjNameMax]; */
  CSeeker *I = &Seeker;    

  row = rowVLA + row_num;
  col = row->col + col_num;
  I->dragging = false;
  I->drag_button = button;
  I->handler.box_active = true;
  I->handler.box_row = row_num;
  I->handler.box_start_col = col_num;
  I->handler.box_stop_col = col_num;

  switch(button) {
  case P_GLUT_RIGHT_BUTTON:
    {
      char name[ObjNameMax];
      ObjectMolecule *obj;
      I->handler.box_active = false;
      if(ExecutiveGetActiveSeleName(name, false) && col->inverse) {
        MenuActivate2Arg(x,y+16,x,y,"pick_option",name,name);
      } else if( (obj = ExecutiveFindObjectMoleculeByName(row->name) )) {
        OrthoLineType buffer;
        {
          int *atom_list;
          char prefix[3]="";
          int logging = SettingGet(cSetting_logging);
          
          if(logging==cPLog_pml)
            strcpy(prefix,"_ ");
          
          if( ExecutiveFindObjectByName(row->name)) {
            atom_list = row->atom_lists + col->atom_at;
            
            /* build up a selection consisting of residue atoms */
            
            if((*atom_list)>=0) {
              
              ObjectMoleculeGetAtomSele(obj,*atom_list,buffer);
              
              BuildSeleFromAtomList(row->name,atom_list,cTempSeekerSele,true);
              if(logging) SelectorLogSele(cTempSeekerSele);
              
              MenuActivate2Arg(x,y+16,x,y,"seq_option",buffer,cTempSeekerSele); 
              
            }
          }
        }
      }
    }
    break;
  case P_GLUT_MIDDLE_BUTTON:
    I->drag_start_col = col_num;
    I->drag_last_col = col_num;
    I->drag_row = row_num;
    I->dragging = true;
    SeekerSelectionCenter(rowVLA,row_num,col_num,true);
    if(mod & cOrthoCTRL) 
      SeekerSelectionCenterGo(1);
    else
      SeekerSelectionCenterGo(0);
    break;
  case P_GLUT_LEFT_BUTTON:
    if(!col->spacer) {
      I->drag_start_col = col_num;
      I->drag_last_col = col_num;
      I->drag_row = row_num;
      I->dragging = true;
      I->drag_setting = true;
      if(col->inverse) {
        SeekerSelectionToggle(rowVLA,row_num,col_num,false);
        I->drag_setting = false;
      } else {
        SeekerSelectionToggle(rowVLA,row_num,col_num,true);
      }
    }
    break;
  }

  /*  ExecutiveGetActiveSeleName(selName,true);
      WizardDoSelect(selName); */

  return NULL;
}

static void SeekerRefresh(CSeqRow *rowVLA)
{
  if(rowVLA) {
    CSeqRow *row;
    CSeqCol *col;
    int *atom_list;
    int nRow = VLAGetSize(rowVLA);
    int sele = ExecutiveGetActiveSele();
    int b;
    ObjectMolecule *obj;

    for(b=0;b<nRow;b++) {
      row = rowVLA + b;
      
      if( (obj = ExecutiveFindObjectMoleculeByName(row->name)) ) {
        register int a;
        register AtomInfoType *atInfo = obj->AtomInfo;
        register int at;
        register int selected;
        register int not_selected;
        
        if(sele<0) {
          for(a=0;a<row->nCol;a++) {
            col = row->col + a;
            col->inverse = false;
          }
        } else {
          for(a=0;a<row->nCol;a++) {
            
            col = row->col + a;
            if(!col->spacer) {
              selected = false;
              atom_list = row->atom_lists + col->atom_at;
              not_selected = true;
              
              while( (at=(*atom_list)) >=0) {
                atom_list++;
                if(SelectorIsMember(atInfo[at].selEntry,sele)) {
                  selected = true; 
                } else {
                  not_selected = true;
                }
              }
              
              if(selected)
                col->inverse = true;
              else
                col->inverse = false;
            } else 
              col->inverse = false;
          }
        }
      }
    }
  }
}

static CSeqRow* SeekerDrag(CSeqRow* rowVLA,int row,int col,int mod)
{
  CSeeker *I = &Seeker;    
  int a;

  if(I->dragging) {
    I->handler.box_stop_col = col;
    
    switch(I->drag_button) {
    case P_GLUT_LEFT_BUTTON:
      if(col != I->drag_last_col) {

        if((I->drag_last_col<I->drag_start_col) && (col>I->drag_start_col))
          {
            for(a=I->drag_last_col;a<I->drag_start_col;a++)
              SeekerSelectionToggle(rowVLA,I->drag_row,a,!I->drag_setting);                        
            I->drag_last_col = I->drag_start_col;
          }
        if((I->drag_last_col>I->drag_start_col) && (col<I->drag_start_col))
          {
            for(a=I->drag_last_col;a>I->drag_start_col;a--)
              SeekerSelectionToggle(rowVLA,I->drag_row,a,!I->drag_setting);                        
            I->drag_last_col = I->drag_start_col;
          }
        if(I->drag_start_col == I->drag_last_col) {
          if(col>I->drag_start_col) {
            I->drag_last_col = I->drag_start_col+1;
            SeekerSelectionToggle(rowVLA,I->drag_row,I->drag_last_col,I->drag_setting);          
          } else if(col<I->drag_start_col){
            I->drag_last_col = I->drag_start_col-1;          
            SeekerSelectionToggle(rowVLA,I->drag_row,I->drag_last_col,I->drag_setting);          
          }
        }
        if(I->drag_start_col < I->drag_last_col) {
          
          if( col > I->drag_last_col ) {
            for( a=I->drag_last_col+1; a<=col; a++) {
              SeekerSelectionToggle(rowVLA,I->drag_row,a,I->drag_setting);          
            }
          } else {
            for(a=I->drag_last_col; a>col ;a--) {
              SeekerSelectionToggle(rowVLA,I->drag_row,a,!I->drag_setting);          
            }
          }
        } else {
          
          if( col < I->drag_last_col) {
            for(a=I->drag_last_col-1;a>=col;a--) {
              SeekerSelectionToggle(rowVLA,I->drag_row,a,I->drag_setting);          
            }
          } else {
            for(a=I->drag_last_col; a<col ;a++) {
              SeekerSelectionToggle(rowVLA,I->drag_row,a,!I->drag_setting);          
            }
          }
        }
        I->drag_last_col = col;              
        
      }
      break;
    case P_GLUT_MIDDLE_BUTTON:
      {
        int action=0;
        int start_over = false;
        if(mod & cOrthoCTRL) {
          action = 1;
        }
        if(!(mod & cOrthoSHIFT)) {
          start_over = true;
          I->handler.box_start_col = col;
          SeekerSelectionCenter(rowVLA,I->drag_row,col,start_over);          
        } else {
          if(I->drag_start_col == I->drag_last_col) {
            if(col>I->drag_start_col) {
              I->drag_last_col = I->drag_start_col+1;
              SeekerSelectionCenter(rowVLA,I->drag_row,I->drag_last_col,start_over);          
            } else if(col<I->drag_start_col) {
              I->drag_last_col = I->drag_start_col-1;          
              SeekerSelectionCenter(rowVLA,I->drag_row,I->drag_last_col,start_over);          
            }
          }
          if(I->drag_start_col < I->drag_last_col) {
            
            if( col > I->drag_last_col ) {
              for( a=I->drag_last_col+1; a<=col; a++) {
                SeekerSelectionCenter(rowVLA,I->drag_row,a,start_over);          
              }
            }
          } else {
            
            if( col < I->drag_last_col) {
              for(a=I->drag_last_col-1;a>=col;a--) {
                SeekerSelectionCenter(rowVLA,I->drag_row,a,start_over);          
              }
            }
          }
        }
        I->drag_last_col = col;              
        
        SeekerSelectionCenterGo(action);
      }
      break;
    }
  }
  return NULL;
}
  
static CSeqRow* SeekerRelease(CSeqRow* rowVLA,int button,
                              int row,int col,int mod)
{
  CSeeker *I = &Seeker;    
  I->dragging = false;

  I->handler.box_active=false;
  return NULL;
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
  int codes;

  CSeqRow *row_vla,*row,*lab=NULL;
  row_vla = VLACalloc(CSeqRow,10);
  /* FIRST PASS: get all the residues represented properly */

  while(ExecutiveIterateObjectMolecule(&obj,&hidden)) {
    if(obj->Obj.Enabled&&(SettingGet_b(obj->Obj.Setting,NULL,cSetting_seq_view))&&
       (obj->Obj.Name[0]!='_')) {
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

      codes = SettingGet_i(obj->Obj.Setting,NULL,cSetting_seq_view_codes);
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
      row->atom_lists = VLACalloc(int,obj->NAtom+est_col+1);
      row->atom_lists[0] = -1; /* terminate the blank listQ (IMPORTANT!) */
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
        /* blank equivalent text for sequence row below the fixed label */
        VLACheck(row->col,CSeqCol,nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        row->len += st_len;
        r1->stop = row->len;
        memset(row->txt+r1->start,32,st_len);
        r1->spacer = true;
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
          r1->spacer = true;
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
          r1->spacer = true;
          nCol++;
          
          last_abbr = false;
          last_chain = ai;
        }

        if(min_pad<0)
          min_pad = strlen(ai->resi) + row->len + 1;
        
        switch(codes) {
        case 0: /* one letter residue codes */
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
        case 1: /* explicit residue codes */
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
        case 2: /* atom names */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col+nCol;
          r1->start = row->len;
          first_atom_in_label = true;
          UtilConcatVLA(&row->txt,&row->len,ai->name);
          r1->stop = row->len;
          UtilConcatVLA(&row->txt,&row->len," ");
          nCol++;
          break;
        case 3:
          if(!AtomInfoSameChainP(last,ai)) {
            last = ai;

            VLACheck(row->col,CSeqCol,nCol);
            r1 = row->col+nCol;
            r1->start = row->len;
            first_atom_in_label = true;

            UtilConcatVLA(&row->txt,&row->len,ai->chain);
            r1->stop = row->len;
            UtilConcatVLA(&row->txt,&row->len," ");
            nCol++;
          }
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
  if(nRow)
    {
    int a,b;
    int nCol;
    int maxCol = 0;
    int done_flag = false;
    /* find out the maximum number of columns */

    for(a=0;a<nRow;a++) {
      row = row_vla + a;
      nCol = row->nCol;
      row->accum = 0; /* initialize the accumulators */
      row->current = 0;
      if(maxCol<nCol)
        maxCol = nCol;
    }

    /* in the simplest mode, just start each sequence in the same column */

    switch(0) {
    case 0:
      b = 0;
      while(!done_flag) {
        int max_offset = 0;
        done_flag = true;
        for(a=0;a<nRow;a++) {
          row = row_vla + a;
          if(!row->label_flag) {
            if(b< row->nCol) {
              CSeqCol *r1 = row->col + b;
              done_flag = false;
              
              r1->offset = r1->start + row->accum;
              if(max_offset<r1->offset)
                max_offset = r1->offset;
            }
          }
        }
        for(a=0;a<nRow;a++) { 
          row = row_vla + a;
          if(!row->label_flag) {
            if(b<row->nCol) {
              CSeqCol *r1 = row->col + b;
              if(b<2) {
                if(r1->offset<max_offset) {
                  row->accum += max_offset - r1->offset;
                }
              }
              r1->offset = r1->start + row->accum;
            }
          }
        }
        b++;
      }
      break;
    }
    

    for(a=0;a<nRow;a++) {
      row = row_vla + a;
      nCol = row->nCol;
      if(row->label_flag)
        lab=row;
      else {
        for(b=0;b<nCol;b++) {
          CSeqCol *r1 = row->col + b,*l1=NULL;
          if(lab) { 
            l1 = lab->col + b; /* if a fixed label is present, get the final offset from the residue line */
            if(l1->stop) 
              l1->offset = r1->offset;
          }
        }
        lab = NULL;
      }
    }
    
  }

  /* THIRD PASS: fill in labels, based on actual residue spacing */

  if(nRow) {
    int a,b,c;
    int nCol;
    for(a=0;a<nRow;a++) {
      lab = row_vla + a;
      if(lab->label_flag) {
        int next_open = 0;
        int *atom_list;
        int st_len;
        int div,sub;
        int draw_it;
        ObjectMolecule *obj;
        AtomInfoType *ai;
        row = lab+1;
        nCol = row->nCol;
        obj = row->obj;

        div = SettingGet_i(obj->Obj.Setting,NULL,cSetting_seq_view_label_spacing);
        sub = SettingGet_i(obj->Obj.Setting,NULL,cSetting_seq_view_label_start);
        for(b=0;b<nCol;b++) {
          CSeqCol *r1 = row->col + b;
          CSeqCol *l1 = lab->col + b;


          if(l1->stop) {/* if label is already present, just line it up */
            l1->offset = r1->offset;
          } else if((r1->offset >= next_open) && (r1->atom_at)) {

            atom_list = row->atom_lists + r1->atom_at;
            ai = obj->AtomInfo + (*atom_list); /* get first atom in list */
              
            if(div>1) {
              if(! ((ai->resv-sub) % div))
                draw_it = true;
              else
                draw_it = false;
            } else {
              draw_it = true;
            }
            if(draw_it) {
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
  }

  /* FOURTH PASS: simply fill in character offsets */
  if(nRow) {
    int a,b;
    int nCol;
    int start,stop;
    for(a=0;a<nRow;a++) {
      row = row_vla + a;
      row->ext_len = 0;

      if(!row->label_flag) {
        nCol = row->nCol;

        for(b=0;b<nCol;b++) {
          CSeqCol *r1 = row->col + b;
          stop = r1->offset + (r1->stop-r1->start);
          if(row->ext_len<stop)
            row->ext_len = stop;
        }
        VLACheck(row->char2col,int,row->ext_len);
        UtilZeroMem(row->char2col,row->ext_len);
        for(b=0;b<nCol;b++) {
          CSeqCol *r1 = row->col + b;
          int c;
          start = r1->offset;
          stop = r1->offset + (r1->stop-r1->start);
          for(c=start;c<stop;c++) 
            row->char2col[c]=b+1;
        }
      }
    }
  }

  SeqSetRowVLA(row_vla,nRow);
  Seeker.handler.fClick = SeekerClick;
  Seeker.handler.fRelease = SeekerRelease;
  Seeker.handler.fDrag = SeekerDrag;
  Seeker.handler.fRefresh = SeekerRefresh;
  SeqSetHandler(&Seeker.handler);
}

void SeekerInit(void)
{
  CSeeker *I = &Seeker;  
  UtilZeroMem(I,sizeof(CSeeker));
}

void SeekerFree(void)
{
}


