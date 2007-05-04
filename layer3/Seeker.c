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
#define cTempSeekerSele2 "_seeker2"

static CSeqRow* SeekerDrag(PyMOLGlobals *G,CSeqRow* rowVLA,int row,int col,int mod);

struct _CSeeker {
  CSeqHandler handler; /* must be first */
  int drag_start_col, drag_last_col;
  int drag_row;
  int drag_dir,drag_start_toggle;
  int dragging, drag_setting;
  int drag_button;
  double LastClickTime;
};


static void SeekerBuildSeleFromAtomList(PyMOLGlobals *G,char *obj_name,int *atom_list,char *sele_name,int start_fresh)
{
  ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(G,obj_name);

  if(start_fresh) {
    SelectorCreateFromObjectIndices(G,sele_name, obj, atom_list,-1);
  } else {
    OrthoLineType buf1;

    SelectorCreateFromObjectIndices(G,cTempSeekerSele2, obj,atom_list,-1);    

    sprintf(buf1,"?%s|?%s",sele_name,cTempSeekerSele2);
    SelectorCreate(G,sele_name,buf1,NULL,true,NULL);      
    ExecutiveDelete(G,cTempSeekerSele2);    
  }
}

static void SeekerSelectionToggleRange(PyMOLGlobals *G,CSeqRow* rowVLA,int row_num,
                                  int col_first,int col_last,int inc_or_excl,
                                  int start_over)
{
  char selName[WordLength];
  OrthoLineType buf1,buf2;

  if(row_num>=0) {
    CSeqRow *row;
    CSeqCol *col;
    char prefix[3]="";
    int logging = SettingGetGlobal_i(G,cSetting_logging);
    int col_num;
    register int *atom_vla = NULL;
    register int n_at = 0;
    register int at_idx;
    register int *atom_list;

    ObjectMolecule *obj;
    if(logging==cPLog_pml)
      strcpy(prefix,"_ ");
    row = rowVLA + row_num;
    if( (obj = ExecutiveFindObjectMoleculeByName(G,row->name)) ) {
      atom_vla = VLAlloc(int,obj->NAtom/10);
      for(col_num=col_first;col_num<=col_last;col_num++) {
        col = row->col + col_num;
        if(!col->spacer) {
          if(!start_over) {
            if(inc_or_excl)
              col->inverse = true;
            else
              col->inverse = false;
          } else {
            col->inverse = true;
          }
          atom_list = row->atom_lists + col->atom_at;
          while((at_idx=(*(atom_list++)))>=0) { /* build one extra long list 
                                    so that we only call selector once*/
            VLACheck(atom_vla,int,n_at);
            atom_vla[n_at++] = at_idx;
          }
        }
      }
      VLACheck(atom_vla,int,n_at);
      atom_vla[n_at]=-1;
      SeekerBuildSeleFromAtomList(G,row->name,atom_vla,cTempSeekerSele,true);
      VLAFreeP(atom_vla);
      
      {      
        char *sele_mode_kw;
        sele_mode_kw = SceneGetSeleModeKeyword(G);
        
        if(logging) SelectorLogSele(G,cTempSeekerSele);
        
        if(!WizardDoSelect(G,cTempSeekerSele)) {
          
          ExecutiveGetActiveSeleName(G,selName,true,logging);
          
          /* selection or deselecting? */
          
          if(!start_over) {
            if(inc_or_excl) {
              sprintf(buf1,"((%s(?%s)) or %s(?%s))",
                      sele_mode_kw,selName,sele_mode_kw,cTempSeekerSele);
            } else {
              sprintf(buf1,"((%s(?%s)) and not %s(?%s))",
                      sele_mode_kw,selName,sele_mode_kw,cTempSeekerSele);
            }
          } else {
            sprintf(buf1,"%s(?%s)",sele_mode_kw,cTempSeekerSele);
          }
          
          /* create the new active selection */
          
          SelectorCreate(G,selName,buf1,NULL,true,NULL);
          {
            sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,selName,buf1);
            PLog(G,buf2,cPLog_no_flush);
          }
        }
        
        ExecutiveDelete(G,cTempSeekerSele);
        if(logging) {
          sprintf(buf2,"%scmd.delete(\"%s\")\n",prefix,cTempSeekerSele);
          PLog(G,buf2,cPLog_no_flush);
          PLogFlush(G);
        }
        
        if(SettingGet(G,cSetting_auto_show_selections))
          ExecutiveSetObjVisib(G,selName,1);
        SceneInvalidate(G);
      }
    }
  }
}

static void SeekerSelectionToggle(PyMOLGlobals *G,CSeqRow* rowVLA,int row_num,
                                  int col_num,int inc_or_excl,
                                  int start_over)
{
  char selName[WordLength];
  OrthoLineType buf1,buf2;

  if(row_num>=0) {
    CSeqRow *row;
    CSeqCol *col;
    int *atom_list;
    char prefix[3]="";
    int logging = SettingGetGlobal_i(G,cSetting_logging);

    if(logging==cPLog_pml)
      strcpy(prefix,"_ ");
    row = rowVLA + row_num;
    col = row->col + col_num;
    if(!col->spacer) 
      if( ExecutiveFindObjectByName(G,row->name)) {
        char *sele_mode_kw;
        atom_list = row->atom_lists + col->atom_at;
        
        /* build up a selection consisting of residue atoms */
        
        SeekerBuildSeleFromAtomList(G,row->name,atom_list,cTempSeekerSele,true);
        sele_mode_kw = SceneGetSeleModeKeyword(G);

        if(logging) SelectorLogSele(G,cTempSeekerSele);
        
        if(!WizardDoSelect(G,cTempSeekerSele)) {
          
          ExecutiveGetActiveSeleName(G,selName,true,logging);
          
          /* selection or deselecting? */

          if(!start_over) {
            if(inc_or_excl) {
              if(!col->spacer) {
                col->inverse = true;
                sprintf(buf1,"((%s(?%s)) or %s(%s))",
                        sele_mode_kw,selName,sele_mode_kw,cTempSeekerSele);
              }
            } else {
              if(!col->spacer) {
                col->inverse = false;
                sprintf(buf1,"((%s(?%s)) and not %s(%s))",
                        sele_mode_kw,selName,sele_mode_kw,cTempSeekerSele);
              }
            }
          } else {
            if(!col->spacer) {
              col->inverse = true;
              sprintf(buf1,"%s(%s)",sele_mode_kw,cTempSeekerSele);
            }
          }
          
          /* create the new active selection */
          
          SelectorCreate(G,selName,buf1,NULL,true,NULL);
          {
            sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,selName,buf1);
            PLog(G,buf2,cPLog_no_flush);
          }
        }

        ExecutiveDelete(G,cTempSeekerSele);
        if(logging) {
          sprintf(buf2,"%scmd.delete(\"%s\")\n",prefix,cTempSeekerSele);
          PLog(G,buf2,cPLog_no_flush);
          PLogFlush(G);
        }
        
        if(SettingGet(G,cSetting_auto_show_selections))
          ExecutiveSetObjVisib(G,selName,1);
        SceneInvalidate(G);
      }
  }
}


static void SeekerSelectionUpdateCenter(PyMOLGlobals *G,CSeqRow* rowVLA,int row_num,int col_num,int start_over)
{
  
  {
    CSeqRow *row;
    CSeqCol *col;
    CObject *obj;

    int *atom_list;
    char prefix[3]="";
    int logging = SettingGetGlobal_i(G,cSetting_logging);


    if(logging==cPLog_pml)
      strcpy(prefix,"_ ");
    if(row_num>=0) {
      row = rowVLA + row_num;
      col = row->col + col_num;
      
      if(!col->spacer)
        if( (obj = ExecutiveFindObjectByName(G,row->name))){
          
          if(col->state&& obj )
            SettingSetSmart_i(G,obj->Setting,NULL,cSetting_state,col->state);
          
          atom_list = row->atom_lists + col->atom_at;
          
          SeekerBuildSeleFromAtomList(G,row->name,atom_list,cTempCenterSele,start_over);
          if(logging) SelectorLogSele(G,cTempCenterSele);
        }
    }
  }

}

static void SeekerSelectionCenter(PyMOLGlobals *G,int action)
{
  OrthoLineType buf2;
  char prefix[3]="";
  int logging = SettingGetGlobal_i(G,cSetting_logging);
  if(logging==cPLog_pml)
    strcpy(prefix,"_ ");
  
  switch(action) {
  case 0: /* center cumulative*/
    ExecutiveCenter(G,cTempCenterSele,-1,true,-1,NULL,true);
    if(logging) {
      sprintf(buf2,"%scmd.center(\"%s\")\n",prefix,cTempCenterSele);
      PLog(G,buf2,cPLog_no_flush);
      PLogFlush(G);
    }
    break;
  case 1: /* zoom */
    ExecutiveWindowZoom(G,cTempCenterSele,0.0,-1,false,-1,true);
    if(logging) {
      sprintf(buf2,"%scmd.zoom(\"%s\")\n",prefix,cTempCenterSele);
      PLog(G,buf2,cPLog_no_flush);
      PLogFlush(G);
    }
    break;
  case 2: /* center seeker */
    {
      char selName[WordLength];
      if(ExecutiveGetActiveSeleName(G,selName,true,logging)) {
        ExecutiveCenter(G,selName,-1,true,-1,NULL,true);
        if(logging) {
          sprintf(buf2,"%scmd.center(\"%s\")\n",prefix,selName);
          PLog(G,buf2,cPLog_no_flush);
          PLogFlush(G);
        }
      }
    }
    break;
  }
}

#define cDoubleTime 0.35

static CSeqRow* SeekerClick(PyMOLGlobals *G,CSeqRow* rowVLA,int button,int row_num,int col_num,int mod,int x,int y)
{
  CSeqRow *row;
  CSeqCol *col;
  /*  char selName[WordLength]; */
  register CSeeker *I = G->Seeker;    
  int logging = SettingGetGlobal_i(G,cSetting_logging);
  int continuation = false;
  if((row_num<0)||(col_num<0)) {
    switch(button) {
    case P_GLUT_LEFT_BUTTON:
      if((UtilGetSeconds(G)-I->LastClickTime)<cDoubleTime) {
        OrthoLineType buf2;
        char name[WordLength];
        if(ExecutiveGetActiveSeleName(G,name, false,false)) {
          SelectorCreate(G,name,"none",NULL,true,NULL);
          if(SettingGet(G,cSetting_logging)) {
            sprintf(buf2,"cmd.select('%s','none')\n",name);
            PLog(G,buf2,cPLog_no_flush);
          }
          SeqDirty(G);
        }
      }
      I->LastClickTime = UtilGetSeconds(G);
      break;
    }
  } else {
    row = rowVLA + row_num;
    col = row->col + col_num;
    I->dragging = false;
    I->drag_button = button;
    I->handler.box_row = row_num;
    I->handler.box_stop_col = col_num;
    if((I->drag_row==row_num)&&
       (button==P_GLUT_LEFT_BUTTON) &&
       (mod & cOrthoSHIFT)) {
      continuation = true;
    } else {
      I->drag_row = -1; /* invalidate */
      I->handler.box_start_col = col_num;
    }
    
    switch(button) {
    case P_GLUT_RIGHT_BUTTON:
      {
        ObjectMolecule *obj;
        char name[WordLength];

        if(ExecutiveGetActiveSeleName(G,name, false,logging) && col->inverse) {
          MenuActivate2Arg(G,x,y+16,x,y,false,"pick_sele",name,name);
        } else if( (obj = ExecutiveFindObjectMoleculeByName(G,row->name) )) {
          OrthoLineType buffer;
          {
            int *atom_list;
            char prefix[3]="";
            int logging = SettingGetGlobal_i(G,cSetting_logging);
            
            if(logging==cPLog_pml)
              strcpy(prefix,"_ ");
            
            if( ExecutiveFindObjectByName(G,row->name)) {
              atom_list = row->atom_lists + col->atom_at;
              
              /* build up a selection consisting of residue atoms */
              
              if((*atom_list)>=0) {
                
                ObjectMoleculeGetAtomSele(obj,*atom_list,buffer);
                
                SeekerBuildSeleFromAtomList(G,row->name,atom_list,cTempSeekerSele,true);
                if(logging) SelectorLogSele(G,cTempSeekerSele);
                
                MenuActivate2Arg(G,x,y+16,x,y,false,"seq_option",cTempSeekerSele,buffer); 
                
              }
            }
          }
        }
      }
      break;
    case P_GLUT_MIDDLE_BUTTON:
      if(!col->spacer) {
        ObjectMolecule *obj;
        I->drag_start_col = col_num;
        I->drag_last_col = col_num;
        I->drag_row = row_num;
        I->dragging = true;
        SeekerSelectionUpdateCenter(G,rowVLA,row_num,col_num,true);
        if(mod & cOrthoCTRL) 
          SeekerSelectionCenter(G,1);
        else
          SeekerSelectionCenter(G,0);
        I->handler.box_active=true;
        if(col->state && (obj = ExecutiveFindObjectMoleculeByName(G,row->name) )) {
          SettingSetSmart_i(G,obj->Obj.Setting,NULL,cSetting_state,col->state);
          SceneChanged(G);
        }
      }
      break;
    case P_GLUT_LEFT_BUTTON:
      if(!col->spacer) {
        int start_over=false;
        int center = 0;
        ObjectMolecule *obj;
        if(mod & cOrthoCTRL) {
          center = 2;
        }
        if(!continuation) {
          I->drag_start_col = col_num;
          I->drag_last_col = col_num;
          I->drag_row = row_num;
          I->drag_dir = 0;
          I->drag_start_toggle = true;
        } else {
          int tmp;
          if(((col_num<I->drag_start_col)&&(I->drag_last_col>I->drag_start_col)) ||
             ((col_num>I->drag_start_col)&&(I->drag_last_col<I->drag_start_col))) {
              tmp = I->drag_last_col;
              I->drag_last_col=I->drag_start_col;
              I->drag_start_col = tmp;
              I->drag_dir = -I->drag_dir;
          }
        }
        I->dragging = true;


        I->handler.box_active=true;
        if(continuation) {
          SeekerDrag(G,rowVLA,row_num,col_num,mod);
        } else {
          if(col->inverse&&!start_over) {
            SeekerSelectionToggle(G,rowVLA,row_num,col_num,false,false);
            I->drag_setting = false;
          } else {
            SeekerSelectionToggle(G,rowVLA,row_num,col_num,true,start_over);
            I->drag_setting = true;
          }
        }
        if(center)
          SeekerSelectionCenter(G,2);

        if(col->state && (obj = ExecutiveFindObjectMoleculeByName(G,row->name))) {
          SettingSetSmart_i(G,obj->Obj.Setting,NULL,cSetting_state,col->state);
          SceneChanged(G);
        }
      }
      break;
    }
  }

  return NULL;
}

static void SeekerRefresh(PyMOLGlobals *G,CSeqRow *rowVLA)
{
  if(rowVLA) {
    CSeqRow *row;
    CSeqCol *col;
    int *atom_list;
    int nRow = VLAGetSize(rowVLA);
    int sele = ExecutiveGetActiveSele(G);
    int b;
    ObjectMolecule *obj;

    if(sele<0)
      sele = SelectorIndexByName(G,"_seeker_hilight");

    for(b=0;b<nRow;b++) {
      row = rowVLA + b;
      
      if( (obj = ExecutiveFindObjectMoleculeByName(G,row->name)) ) {
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
                if(SelectorIsMember(G,atInfo[at].selEntry,sele)) {
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

static CSeqRow* SeekerDrag(PyMOLGlobals *G,CSeqRow* rowVLA,int row,int col,int mod)
{
  register CSeeker *I = G->Seeker;    
  int a;

  if((row>=0)&&(col>=0)&&(I->dragging)) {
    I->handler.box_stop_col = col;
    
    switch(I->drag_button) {
    case P_GLUT_LEFT_BUTTON:
      if(col != I->drag_last_col) {

        if(I->drag_dir) {
          if(I->drag_dir>0) {
            if(col<=I->drag_start_col) {
              col = I->drag_start_col;
              if(I->drag_start_toggle) {
                SeekerSelectionToggle(G,rowVLA,I->drag_row,I->drag_start_col,!I->drag_setting,false);  
                I->drag_start_toggle = false;
              }
            } else if(col>I->drag_start_col) {
              if(!I->drag_start_toggle) {
                SeekerSelectionToggle(G,rowVLA,I->drag_row,I->drag_start_col,I->drag_setting,false);  
                I->drag_start_toggle = true;
              }
            }
          } else if(I->drag_dir<0) {
            if(col>=I->drag_start_col) {
              col = I->drag_start_col;
              if(I->drag_start_toggle) {
                SeekerSelectionToggle(G,rowVLA,I->drag_row,I->drag_start_col,!I->drag_setting,false);  
                I->drag_start_toggle = false;
              }
            } else if (col<I->drag_start_col) {
              if(!I->drag_start_toggle) {
                SeekerSelectionToggle(G,rowVLA,I->drag_row,I->drag_start_col,I->drag_setting,false);  
                I->drag_start_toggle = true;
              }
            }
          }
        }
        /*
        if(mod &cOrthoSHIFT) {
          if(I->drag_start_col == I->drag_last_col) {
            if(col>I->drag_start_col) {
              SeekerSelectionCenter(G,rowVLA,I->drag_row,I->drag_start_col+1,false);
            } else if(col<I->drag_start_col) {
              SeekerSelectionCenter(G,rowVLA,I->drag_row,I->drag_start_col-1,false);
            }
          }
          if(I->drag_start_col < I->drag_last_col) {
            if( col > I->drag_last_col ) {
              for( a=I->drag_last_col+1; a<=col; a++) {
                SeekerSelectionCenter(G,rowVLA,I->drag_row,a,false);
              }
            }
          } else {
            
            if( col < I->drag_last_col) {
              for(a=I->drag_last_col-1;a>=col;a--) {
                SeekerSelectionCenter(G,rowVLA,I->drag_row,a,false);
              }
            }
          }
          SeekerSelectionCenter(G,0);
        }
        */

        if((I->drag_last_col<I->drag_start_col) && (col>I->drag_start_col))
          {
            /*            for(a=I->drag_last_col;a<I->drag_start_col;a++)*/
            SeekerSelectionToggleRange(G,rowVLA,I->drag_row,I->drag_last_col,I->drag_start_col-1,!I->drag_setting,false);  
            I->drag_last_col = I->drag_start_col;
          }
        if((I->drag_last_col>I->drag_start_col) && (col<I->drag_start_col))
          {
            /*            for(a=I->drag_last_col;a>I->drag_start_col;a--)*/
            SeekerSelectionToggleRange(G,rowVLA,I->drag_row,I->drag_start_col+1,I->drag_last_col,!I->drag_setting,false);
            I->drag_last_col = I->drag_start_col;
          }
        if(I->drag_start_col == I->drag_last_col) {
          if(col>I->drag_start_col) {
            if(!I->drag_dir)
              I->drag_dir = 1;
            I->drag_last_col = I->drag_start_col+1;
            SeekerSelectionToggle(G,rowVLA,I->drag_row,I->drag_last_col,I->drag_setting,false);
          } else if(col<I->drag_start_col){
            if(!I->drag_dir)
              I->drag_dir = -1;
            I->drag_last_col = I->drag_start_col-1;          
            SeekerSelectionToggle(G,rowVLA,I->drag_row,I->drag_last_col,I->drag_setting,false);
          }
        }
        if(I->drag_start_col < I->drag_last_col) {
          
          if( col > I->drag_last_col ) {
            /*            for( a=I->drag_last_col+1; a<=col; a++) */
            SeekerSelectionToggleRange(G,rowVLA,I->drag_row,I->drag_last_col+1,col,I->drag_setting,false);          
          } else {
            /*            for(a=I->drag_last_col; a>col ;a--) */
            SeekerSelectionToggleRange(G,rowVLA,I->drag_row,col+1,I->drag_last_col,!I->drag_setting,false);          
          }
        } else {
          
          if( col < I->drag_last_col) {
            /*for(a=I->drag_last_col-1;a>=col;a--) */
            SeekerSelectionToggleRange(G,rowVLA,I->drag_row,col,I->drag_last_col-1,I->drag_setting,false);          
          } else {
            /*for(a=I->drag_last_col; a<col ;a++) */
            SeekerSelectionToggleRange(G,rowVLA,I->drag_row,I->drag_last_col,col-1,!I->drag_setting,false);          
          }
        }
        I->drag_last_col = col;              
       
        if(mod & cOrthoCTRL) {
          SeekerSelectionCenter(G,2);
        }
 
      }
      break;
    case P_GLUT_MIDDLE_BUTTON:
      if(col != I->drag_last_col) {
        int action=0;
        int start_over = false;

        if(mod & cOrthoCTRL) {
          action = 1;
        }
        if(!(mod & cOrthoSHIFT)) {
          start_over = true;
          I->handler.box_start_col = col;
          SeekerSelectionUpdateCenter(G,rowVLA,I->drag_row,col,start_over);          
        } else {
          if(I->drag_start_col == I->drag_last_col) {
            if(col>I->drag_start_col) {
              I->drag_last_col = I->drag_start_col+1;
              SeekerSelectionUpdateCenter(G,rowVLA,I->drag_row,I->drag_last_col,start_over);          
            } else if(col<I->drag_start_col) {
              I->drag_last_col = I->drag_start_col-1;          
              SeekerSelectionUpdateCenter(G,rowVLA,I->drag_row,I->drag_last_col,start_over);          
            }
          }
          if(I->drag_start_col < I->drag_last_col) {
            
            if( col > I->drag_last_col ) {
              for( a=I->drag_last_col+1; a<=col; a++) {
                SeekerSelectionUpdateCenter(G,rowVLA,I->drag_row,a,start_over);          
              }
            }
          } else {
            
            if( col < I->drag_last_col) {
              for(a=I->drag_last_col-1;a>=col;a--) {
                SeekerSelectionUpdateCenter(G,rowVLA,I->drag_row,a,start_over);          
              }
            }
          }
        }
        I->drag_last_col = col;              
        
        SeekerSelectionCenter(G,action);
      }
      break;
    }
  }
  return NULL;
}
  
static CSeqRow* SeekerRelease(PyMOLGlobals *G,CSeqRow* rowVLA,int button,
                              int row,int col,int mod)
{
  register CSeeker *I = G->Seeker;    
  I->dragging = false;

  I->handler.box_active=false;
  return NULL;
}

char SeekerGetAbbr(PyMOLGlobals *G,char *abbr,char water,char unknown)
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
        return 'E';
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
        return water;
        break;
      }
      break;
    case '2': 
      switch(abbr[2]) {
      case 'O':
        return water;
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
    case 'O':  /* SOL -- gromacs solvent residue */
      switch(abbr[2]) {
      case 'L':
        return water;
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
        return 'T';
        break;
      }
      break;
    case 'I': 
      switch(abbr[2]) {
      case 'P':
        return water;
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
  case 'W':
    switch(abbr[1]) {
    case 'A': 
      switch(abbr[2]) {
      case 'T':
        return water;
        break;
      }
      break;
    }
    break;

  }

  return unknown;
}

static int SeekerFindColor(PyMOLGlobals *G,AtomInfoType *ai,int n_more_plus_one)
{
  register int result = ai->color; /* default -- use first atom color */
  register AtomInfoType *ai0 =ai;
  while(1) {
    if(ai0->flags & cAtomFlag_guide) /* best use guide color */
      return ai0->color;
    if(ai0->protons == cAN_C) /* or use carbon color */
      result = ai0->color;
    n_more_plus_one--;
    if(n_more_plus_one>0) {
      ai0++;
      if(!AtomInfoSameResidueP(G,ai,ai0))
        break;
    } else 
      break;
  }
  return result;
}

static int SeekerFindTag(PyMOLGlobals *G,AtomInfoType *ai,int sele, int codes,int n_more_plus_one)
{
  register int result = 0;/* default -- no tag */
  register AtomInfoType *ai0 =ai;
  while(1) {
    int tag = SelectorIsMember(G,ai0->selEntry, sele);
    if(tag && (codes<2) && (ai0->flags & cAtomFlag_guide)) /* use guide atom if present */
      return tag;
    if(result<tag) {
      if(!result)
        result = tag;
      else if((codes<2) && (ai0->flags & cAtomFlag_guide)) /* residue based and on guide atom */
        result = tag;
    }
    n_more_plus_one--;
    if(n_more_plus_one>0) {
      int do_break = false;
      ai0++;
      switch(codes) {
      case 0:
      case 1:
        if(!AtomInfoSameResidueP(G,ai,ai0))
          do_break = true;
        break;
      case 2: /* atoms */
        do_break = true;
        break;
      case 3: /* chains */
        if(!AtomInfoSameChainP(G,ai,ai0))
          do_break = true;
        break;
      }
      if(do_break)
        break;
    } else 
      break;
  }
  return result;
}

PyObject *SeekerGetRawAlignment(PyMOLGlobals *G, int align_sele, int active_only) 
{
  PyObject *result = NULL;
#ifdef _PYMOL_NOPY
  return NULL;
#else
  int nRow = 0;
  int nCol = 0;
  CSeqRow *row_vla = NULL,*row;
  void *hidden = NULL;
  ObjectMolecule *obj;

  if(align_sele<0) {
    align_sele = ExecutiveGetActiveAlignmentSele(G);  
  }
  if(align_sele>=0) {

    row_vla = VLACalloc(CSeqRow,10);

    /* first, find out which objects are included in the alignment */

    while(ExecutiveIterateObjectMolecule(G,&obj,&hidden)) {
      if((obj->Obj.Enabled || !active_only) && (obj->Obj.Name[0]!='_')) {
        int a;
        AtomInfoType *ai = obj->AtomInfo;
        for(a=0;a<obj->NAtom;a++) {
          if(SelectorIsMember(G,ai->selEntry,align_sele)) {      
            VLACheck(row_vla,CSeqRow,nRow);
            row = row_vla + nRow;
            row->obj = obj;
            row->nCol = obj->NAtom;
            nRow++;
            break;
          }
          ai++;
        }
      }
    }

    /* next, figure out how many aligned columns exist */

    {
      int done = false;
      while(!done) {
        int a;
        int min_tag = -1;
        done = true;
        for(a=0;a<nRow;a++) {
          row = row_vla + a;
          while(row->cCol<row->nCol) { /* advance to next tag in each row & find lowest */
            AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
            int tag = SelectorIsMember(G,ai->selEntry, align_sele);     
            if(!tag) {
              row->cCol++;
            } else { /* we're at a tagged atom... */
              if(min_tag>tag)
                min_tag = tag; 
              else if(min_tag<0)
                min_tag = tag;
              done = false;
              break;
            }
          }
        }
        if(min_tag>=0) {
          nCol++;
          for(a=0;a<nRow;a++) {
            row = row_vla + a;
            if(row->cCol<row->nCol) {
              AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
              int tag = SelectorIsMember(G,ai->selEntry, align_sele);     
              if(tag == min_tag) { /* advance past this tag */
                row->cCol++;
              }
            }
          }
        }
      }
    }

    /* now populate the table */

    result = PyList_New(nCol);
  
    if(nCol) {
      int done = false;
      nCol = 0;
      { /* reset start points for our second pass */
        int a;
        for(a=0;a<nRow;a++) {
          row = row_vla + a;
          row->cCol = 0;
        }
      }
      while(!done) {
        int a;
        int min_tag = -1;
        done = true;
        for(a=0;a<nRow;a++) {
          row = row_vla + a;
          while(row->cCol<row->nCol) { /* advance to next tag in each row & find lowest */
            AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
            int tag = SelectorIsMember(G,ai->selEntry, align_sele);     
            if(!tag) {
              row->cCol++;
            } else { /* we're at a tagged atom... */
              if(min_tag>tag)
                min_tag = tag; 
              else if(min_tag<0)
                min_tag = tag;
              done = false;
              break;
            }
          }
        }
        if(min_tag>=0) {
          int n_member = 0;

          for(a=0;a<nRow;a++) {
            row = row_vla + a;
            if(row->cCol<row->nCol) {
              AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
              int tag = SelectorIsMember(G,ai->selEntry, align_sele);     
              if(tag == min_tag) { /* participates */
                n_member++;
              }
            }
          }
          {
            PyObject *column_list = PyList_New(n_member);
            n_member = 0;

            for(a=0;a<nRow;a++) {
              row = row_vla + a;
              if(row->cCol<row->nCol) {
                AtomInfoType *ai = row->obj->AtomInfo + row->cCol;
                int tag = SelectorIsMember(G,ai->selEntry, align_sele);     
                if(tag == min_tag) { /* participates */

                  PyObject *tup = PyTuple_New(2);
                  PyTuple_SetItem(tup,0,PyString_FromString(row->obj->Obj.Name));
                  PyTuple_SetItem(tup,1,PyInt_FromLong(row->cCol));
                  PyList_SetItem(column_list,n_member,tup);

                  row->cCol++; /* advance past this tag */
                  n_member++;
                }
              }
            }
            PyList_SetItem(result,nCol,column_list);
          }
          nCol++;
        }
      }
    }
  }
  VLAFreeP(row_vla);
  return result;
#endif
}

void SeekerUpdate(PyMOLGlobals *G)
{
  /*  CObject *o = NULL;
      int s;*/

  void *hidden = NULL;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  int nRow = 0;
  int label_mode = 0;
  int codes = 0;
  int max_row = 50;
  int default_color = 0;
  int align_sele = -1; /* alignment selection */
  CSeqRow *row_vla,*row,*lab=NULL;
  row_vla = VLACalloc(CSeqRow,10);
  /* FIRST PASS: get all the residues represented properly */
  label_mode = SettingGetGlobal_i(G,cSetting_seq_view_label_mode);

#if 1
  align_sele = ExecutiveGetActiveAlignmentSele(G);
#endif

  while(ExecutiveIterateObjectMolecule(G,&obj,&hidden)) {
    if(obj->Obj.Enabled&&(SettingGet_b(G,obj->Obj.Setting,NULL,cSetting_seq_view))&&
       (obj->Obj.Name[0]!='_')) {
      int a;
      AtomInfoType *last = NULL,*last_segi=NULL,*last_chain = NULL;
      CoordSet *last_disc = NULL;
      int last_state;
      int last_abbr = true;
      int last_spacer = false;
      int nCol = 0;
      int nListEntries = 1; /* first list starts at 1 always... */
      int est_col = obj->NAtom/5+1;
      int est_char = obj->NAtom*4;
      int first_atom_in_label;

      int min_pad = -1;
      CSeqCol *r1 = NULL,*l1=NULL;/* *col */

      if(nRow>=max_row)
        break;

      codes = SettingGet_i(G,obj->Obj.Setting,NULL,cSetting_seq_view_format);
      if(obj->DiscreteFlag && SettingGet_b(G,
                                           obj->Obj.Setting,
                                           NULL,
                                           cSetting_seq_view_discrete_by_state))
        codes = 4;
      default_color = SettingGet_i(G,obj->Obj.Setting,NULL,cSetting_seq_view_color);

      /* allocate a row for labels, if present
         the text for the labels and the residues will line up exactly 
      */

      VLACheck(row_vla,CSeqRow,nRow);
      if((label_mode==2)||((label_mode==1)&&(!nRow))) {
          lab = row_vla + nRow++;
          lab->txt = VLAlloc(char,est_char);
          lab->col = VLACalloc(CSeqCol,est_col);
          lab->label_flag = true;
      } else {
        lab = NULL;
      }

      VLACheck(row_vla,CSeqRow,nRow);

      row = row_vla+nRow;
      if(lab) lab = row-1; /* critical! */
      row->txt = VLAlloc(char,est_char);
      row->col = VLACalloc(CSeqCol,est_col);
      row->fill = VLACalloc(CSeqCol,est_col/8);
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
        UtilConcatVLA(&lab->txt,&lab->len,"/");
        UtilConcatVLA(&lab->txt,&lab->len,obj->Obj.Name);
        l1->stop = lab->len;
        st_len = l1->stop - l1->start;

        if(label_mode==2) {
          /* blank equivalent text for sequence row below the fixed label */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          UtilFillVLA(&row->txt,&row->len,' ',st_len);
          r1->stop = row->len;
          r1->spacer = true;
          nCol++;
        } 
      } 
      if(label_mode<2) { /* no label rows, so put object name into left-hand column */

        /* copy label text */

        VLACheck(row->col,CSeqCol,nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilConcatVLA(&row->txt,&row->len,"/");
        UtilConcatVLA(&row->txt,&row->len,obj->Obj.Name);
        r1->stop = row->len;
        r1->spacer = true;
        row->column_label_flag = true;
        row->title_width = row->len;
        nCol++;
      } else if(label_mode==3) { /* otherwise just insert a blank zero-length column */
        VLACheck(row->col,CSeqCol,nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilConcatVLA(&row->txt,&row->len,"");
        r1->stop = row->len;
        r1->spacer = true;
        nCol++;
      }

      if(lab) {
        
        int st_len;
        /* copy label text */

        VLACheck(lab->col,CSeqCol,nCol);
        l1 = lab->col + nCol;
        l1->start = lab->len;
        if(obj->NAtom) {
          UtilConcatVLA(&lab->txt,&lab->len,"/");
          UtilConcatVLA(&lab->txt,&lab->len,ai->segi);
          UtilConcatVLA(&lab->txt,&lab->len,"/");
          UtilConcatVLA(&lab->txt,&lab->len,ai->chain);
          UtilConcatVLA(&lab->txt,&lab->len,"/");
        } else {
          UtilConcatVLA(&lab->txt,&lab->len,"///");
        }

        l1->stop = lab->len;
        st_len = l1->stop - l1->start;

        last_segi = ai;
        last_chain = ai;
        /* blank equivalent text for sequence row below the fixed label */
        VLACheck(row->col,CSeqCol,nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilFillVLA(&row->txt,&row->len,' ',st_len);
        r1->stop = row->len;
        r1->spacer = true;
        nCol++;
      } else { /* if no labels, just insert a space in row below */
        VLACheck(row->col,CSeqCol,nCol);
        r1 = row->col + nCol;
        r1->start = row->len;
        UtilConcatVLA(&row->txt,&row->len," ");
        r1->stop = row->len;
        r1->spacer = true;
        nCol++;
      }


      last_state=-1;
      for(a=0;a<obj->NAtom;a++) {
        first_atom_in_label = false;
        if(lab&&!AtomInfoSameSegmentP(G,last_segi,ai)) {

          int st_len;

          if(row->len<min_pad) {
            row->len = min_pad;
          }
          min_pad = -1;

          /* copy label text */
          
          VLACheck(lab->col,CSeqCol,nCol);
          l1 = lab->col + nCol;
          l1->start = lab->len;
          UtilConcatVLA(&lab->txt,&lab->len,"/");
          UtilConcatVLA(&lab->txt,&lab->len,ai->segi);
          UtilConcatVLA(&lab->txt,&lab->len,"/");
          UtilConcatVLA(&lab->txt,&lab->len,ai->chain);
          UtilConcatVLA(&lab->txt,&lab->len,"/");
          l1->stop = lab->len;
          st_len = l1->stop - l1->start;
          
          /* blank equivalent text for sequence row */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          UtilFillVLA(&row->txt,&row->len,' ',st_len);
          r1->stop = row->len;
          r1->spacer = true;
          nCol++;
          
          last_abbr = false;
          last_spacer = true;
          last_segi = ai;
          last_chain = ai;

        } else if(lab&&!AtomInfoSameChainP(G,last_chain,ai)) {

          int st_len;

          if(row->len<min_pad) {
            row->len = min_pad;
          }
          min_pad = -1;

          /* copy label text */
          
          VLACheck(lab->col,CSeqCol,nCol);
          l1 = lab->col + nCol;
          l1->start = lab->len;
          UtilConcatVLA(&lab->txt,&lab->len,"/");
          UtilConcatVLA(&lab->txt,&lab->len,ai->chain);
          UtilConcatVLA(&lab->txt,&lab->len,"/");
          l1->stop = lab->len;
          st_len = l1->stop - l1->start;
          
          /* blank equivalent text for sequence row */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col + nCol;
          r1->start = row->len;
          UtilFillVLA(&row->txt,&row->len,' ',st_len);
          r1->stop = row->len;
          r1->spacer = true;
          nCol++;
          
          last_abbr = false;
          last_spacer = true;
          last_chain = ai;
        }

        if(min_pad<0)
          min_pad = strlen(ai->resi) + row->len + 1;
        
        switch(codes) {
        case 0: /* one letter residue codes */
          if(!AtomInfoSameResidueP(G,last,ai)) {
            char abbr[2] = "1";
            last = ai;            

            VLACheck(row->col,CSeqCol,nCol);
            r1 = row->col+nCol;
            r1->start = row->len;
            if(obj->DiscreteFlag) 
              r1->state = ai->discrete_state;

            first_atom_in_label = true;

            abbr[0] = SeekerGetAbbr(G,ai->resn,'O',0);

            r1->hint_no_space = last_abbr || last_spacer;

            if(!abbr[0]) {
              if(last_abbr) {
                UtilConcatVLA(&row->txt,&row->len," ");                
                r1->start = row->len;
              }
              
              if(ai->resn[0])
                UtilConcatVLA(&row->txt,&row->len,ai->resn);
              else
                UtilConcatVLA(&row->txt,&row->len,"''");
              
              r1->stop = row->len;
              
              UtilConcatVLA(&row->txt,&row->len," ");
            } else {
              UtilConcatVLA(&row->txt,&row->len,abbr);
              r1->is_abbr = true;
              r1->stop = row->len;
            }
            if(default_color<0)
              r1->color = SeekerFindColor(G,ai,obj->NAtom-a);
            else
              r1->color = default_color;
            if(align_sele>=0) {
              r1->tag = SeekerFindTag(G,ai,align_sele,codes,obj->NAtom-a);
            } else {
              r1->tag = 0;
            }
            nCol++;
            last_abbr=abbr[0];
          }
          
          break;
        case 1: /* explicit residue codes */
          if(!AtomInfoSameResidueP(G,last,ai)) {
            last = ai;

            VLACheck(row->col,CSeqCol,nCol);
            r1 = row->col+nCol;
            r1->start = row->len;
            if(obj->DiscreteFlag) 
              r1->state = ai->discrete_state;
            first_atom_in_label = true;

            if(ai->resn[0])
              UtilConcatVLA(&row->txt,&row->len,ai->resn);
            else
              UtilConcatVLA(&row->txt,&row->len,"''");
            r1->stop = row->len;
            if(default_color<0)
              r1->color = SeekerFindColor(G,ai,obj->NAtom-a);
            else
              r1->color = default_color;
            if(align_sele>=0) {
              r1->tag = SeekerFindTag(G,ai,align_sele,codes,obj->NAtom-a);
            } else {
              r1->tag = 0;
            }
            UtilConcatVLA(&row->txt,&row->len," ");
            nCol++;
          }
          break;
        case 2: /* atom names */
          VLACheck(row->col,CSeqCol,nCol);
          r1 = row->col+nCol;
          r1->start = row->len;
          first_atom_in_label = true;
          if(ai->name[0])
            UtilConcatVLA(&row->txt,&row->len,ai->name);
          else
            UtilConcatVLA(&row->txt,&row->len,"''");
          r1->stop = row->len;
          if(default_color<0)
            r1->color = ai->color;
          else
            r1->color = default_color;
          if(align_sele>=0) {
            r1->tag = SeekerFindTag(G,ai,align_sele,codes,obj->NAtom-a);
          } else {
            r1->tag = 0;
          }
          UtilConcatVLA(&row->txt,&row->len," ");
          nCol++;
          break;
        case 3: /* chains */
          if(!AtomInfoSameChainP(G,last,ai)) {
            last = ai;

            VLACheck(row->col,CSeqCol,nCol);
            r1 = row->col+nCol;
            r1->start = row->len;
            first_atom_in_label = true;

            if(ai->chain[0])
              UtilConcatVLA(&row->txt,&row->len,ai->chain);
            else
              UtilConcatVLA(&row->txt,&row->len,"''");
            r1->stop = row->len;
            if(default_color<0)
              r1->color = SeekerFindColor(G,ai,obj->NAtom-a);
            else
              r1->color = default_color;
            if(align_sele>=0) {
              r1->tag = SeekerFindTag(G,ai,align_sele,codes,obj->NAtom-a);
            } else {
              r1->tag = 0;
            }
            UtilConcatVLA(&row->txt,&row->len," ");
            nCol++;
          }
          break;
        case 4: /* state names */
          if(obj->DiscreteFlag) {
            CoordSet *cs;
            if((cs = obj->DiscreteCSet[a])!=last_disc) {
              last_disc = cs;
              if(cs) {
                default_color = SettingGet_i(G,cs->Setting,obj->Obj.Setting,
                                             cSetting_seq_view_color);
                VLACheck(row->col,CSeqCol,nCol);
                r1 = row->col+nCol;
                r1->start = row->len;
                r1->color = default_color;
                first_atom_in_label = true;
                
                if(cs->Name[0])
                  UtilConcatVLA(&row->txt,&row->len,cs->Name);
                else
                  UtilConcatVLA(&row->txt,&row->len,"''");
                r1->stop = row->len;
                r1->state = ai->discrete_state;
                UtilConcatVLA(&row->txt,&row->len," ");
                nCol++;
              }
            }
          } else { 
            /* non-discrete objects simply get their states enumerated
               without selections */
            
            if(last_state<0) {
              int b;
              CoordSet *cs;
              WordType buf1;
              last_state = 1;
              first_atom_in_label = true;
              for(b=0;b<obj->NCSet;b++) {
                cs = obj->CSet[b];
                if(cs) {
                  default_color = SettingGet_i(G,cs->Setting,obj->Obj.Setting,
                                               cSetting_seq_view_color);
                  
                  VLACheck(row->col,CSeqCol,nCol);
                  r1 = row->col+nCol;
                  r1->state = b+1;
                  r1->start = row->len;
                  r1->atom_at = nListEntries + 1; /* tricky & dangerous */
                  r1->color = default_color;
                  if(cs->Name[0])
                    UtilConcatVLA(&row->txt,&row->len,cs->Name);
                  else {
                    sprintf(buf1,"%d",b+1);
                    UtilConcatVLA(&row->txt,&row->len,buf1);
                  }
                  r1->stop = row->len;
                  UtilConcatVLA(&row->txt,&row->len," ");
                  nCol++;
                }
              }
            }
          }
          break;
        case 5: /* movie frames */
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
        /*        if(lab->len<row->len) {
          lab->len = row->len;
          }*/
        VLASize(lab->txt,char,lab->len+1);
        lab->txt[lab->len] = 0;
        VLACheck(lab->col,CSeqCol,nCol); /* make sure we've got column records for labels too */
        lab->nCol = nCol;

        /*if(row->len<lab->len) {
          row->len = lab->len;
          }*/
      }

      VLASize(row->txt,char,row->len+1);
      row->txt[row->len]=0;

      row->nCol = nCol;

      /* terminate last atom list */
      VLACheck(row->atom_lists,int,nListEntries);
      row->atom_lists[nListEntries] = -1;
      nListEntries++;
      nRow++;
    }
  }

  /* SECOND PASS: align columns to reflect current alignment and fixed labels */
  if(nRow) {
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

    if(align_sele<0) {
      /* in the simplest mode, just start each sequence in the same column */

      b = 0;
      while(!done_flag) {
        int max_offset = 0;
        done_flag = true;
        for(a=0;a<nRow;a++) {
          row = row_vla + a;
          if(!row->label_flag) {
            if(b < row->nCol) {
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
              if(b<3) { 
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
    } else {
      /* in alignment mode, line up the tags */
      int stagger = false;
      /* intialize current columns and get the starting character */
      int current = 0;
      int first = true;
      switch(SettingGetGlobal_i(G,cSetting_seq_view_unaligned_mode)) {
      case 0:
      case 1:
      case 2:
        stagger = false;
        break;
      default:
        stagger = true;
        break;
      }
      for(a=0;a<nRow;a++) {       
        row = row_vla + a;
        row->cCol = 0;
        if((!row->label_flag) && ((row->cCol < row->nCol))) {
          if(current < row->accum)
            current = row->accum;
        }
      }
      done_flag = false;
      while(!done_flag) {
        int hint_tagged_no_space = true;
        done_flag=true;
        {
#if 0
          if(all_spacers) {
            /* this column is only spacers, so line them up like normal */
            int max_width = 0;
            int width;
            int space_added = false;
            int rep;
            for(rep=0;rep<2;rep++)
              for(a=0;a<nRow;a++) {       
                row = row_vla + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if( (!first) && (!space_added) &&
                      (codes || (((!r1->is_abbr) && (!r1->spacer))) ||
                       (r1->is_abbr&&(!r1->hint_no_space)))) {
                    current++;
                    space_added = true;
                  }
                  done_flag = false;
                  first = false;
                  r1->offset = current;
                  width =  (r1->stop-r1->start);
                  if(max_width<width)
                    max_width = width;
                  row->cCol++;
                }
              }
            current+=max_width;
          }
#endif
          {
            /* insert untagged entries into their own columns */
            int untagged_flag = true;
            int saw_untagged_no_abbr = false;
            int hint_untagged_space = false;
            while(untagged_flag) {
              int space_added = false;
              int max_width = 0;
              untagged_flag = false;
              
              /* first get the spaces in...*/

              for(a=0;a<nRow;a++) {       
                row = row_vla + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(!r1->tag) { /* not aligned */
                    int text_len = (r1->stop-r1->start);
                    if((!first) && (!space_added) && (row->cCol>2) &&
                       (codes || (((!r1->is_abbr) && (!r1->spacer) )) ||
                        hint_untagged_space ||
                        (r1->is_abbr&&(!r1->hint_no_space)))) {
                      /* insert space */
                      current++;
                      space_added = true;
                    }
                    if(max_width<text_len)
                      max_width = text_len;
                  }
                }
              }
              
              /* then do the rest */

              for(a=0;a<nRow;a++) {       
                row = row_vla + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(!r1->tag) { /* not aligned */
                    int text_len = (r1->stop-r1->start);
                    untagged_flag = true;
                    done_flag = false;
                    saw_untagged_no_abbr |= (!r1->is_abbr)&&(!r1->spacer);

                    first = false;
                    r1->offset = current;
                    r1->unaligned = true;
                  
                    if(!r1->spacer) {
                      int aa;

                      for(aa=0;aa<nRow;aa++) { /* infill populate other rows with dashes */
                        if(aa!=a) {       
                          CSeqRow *row2 = row_vla + aa;
                          if(!row2->label_flag) { 
                            if(row2->cCol < row2->nCol) {
                              CSeqCol *r2 = row2->col + row2->cCol;
                              if(stagger||r2->tag||r2->spacer) {
                                VLACheck(row2->fill,CSeqCol,row2->nFill);
                                r2 = row2->fill + row2->nFill;
                                r2->stop = text_len;
                                r2->offset = current;
                                row2->nFill++;
                              }
                            } else {
                              CSeqCol *r2 = row2->col + row2->cCol;
                              VLACheck(row2->fill,CSeqCol,row2->nFill);
                              r2 = row2->fill + row2->nFill;
                              r2->stop = text_len;
                              r2->offset = current;
                              row2->nFill++;
                            }
                          }
                        }
                      }
                    }
                    if(stagger) 
                      current += text_len;
                    else if(max_width<text_len)
                      max_width = text_len;
                  }
                }
              }
              if(!stagger) 
                current += max_width;
              if(saw_untagged_no_abbr) {
                hint_untagged_space = true;
                hint_tagged_no_space = false;
              } else {
                hint_untagged_space = false;
                hint_tagged_no_space = true;
              }
              saw_untagged_no_abbr = false;
              for(a=0;a<nRow;a++) {       
                row = row_vla + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(!r1->tag) {
                    row->cCol++;
                  }
                }
              }
            }
          }
        }
        
        {
          /* next insert match-tagged entries into the same column */
          int min_tag = 0;
          for(a=0;a<nRow;a++) {       
            row = row_vla + a;
            if((!row->label_flag) && (row->cCol < row->nCol)) {
              CSeqCol *r1 = row->col + row->cCol;
              if(r1->tag && ((min_tag>r1->tag)||(!min_tag))) {
                min_tag = r1->tag;
              }
            }
          }
          if(min_tag) {
            int width, max_width = 0;
            int space_added = false;
            int rep;
            for(rep=0;rep<2;rep++) 
              for(a=0;a<nRow;a++) {       
                row = row_vla + a;
                if((!row->label_flag) && (row->cCol < row->nCol)) {
                  CSeqCol *r1 = row->col + row->cCol;
                  if(r1->tag == min_tag) {
                    if((!first) && (!space_added) &&
                       (codes || (((!r1->is_abbr) && (!r1->spacer))) ||
                        (r1->is_abbr&&(!(r1->hint_no_space||hint_tagged_no_space))))) {
                      /* insert space */
                      current++;
                      space_added = true;
                    }
                    done_flag = false;
                    first = false;
                    r1->offset = current;
                    width = (r1->stop-r1->start);
                    if(max_width<width)
                      max_width = width;
                    /*   row->cCol++;  */
                  }
                }
              }
            {
              int aa;
              for(aa=0;aa<nRow;aa++) { /* infill populate other rows with dashes */
                CSeqRow *row2 = row_vla + aa;
                if(!row2->label_flag) {
                  if(row2->cCol < row2->nCol) {
                    CSeqCol *r1 = row2->col + row2->cCol;
                    if( r1->tag != min_tag) {
                      CSeqCol *r2;
                      VLACheck(row2->fill,CSeqCol,row2->nFill);
                      r2 = row2->fill + row2->nFill;
                      r2->stop = max_width;
                      r2->offset = current;
                      row2->nFill++;
                    }
                  } else {
                    CSeqCol *r2;
                    VLACheck(row2->fill,CSeqCol,row2->nFill);
                    r2 = row2->fill + row2->nFill;
                    r2->stop = max_width;
                    r2->offset = current;
                    row2->nFill++;
                  }
                }
              }
            }
            for(a=0;a<nRow;a++) {       
              row = row_vla + a;
              if((!row->label_flag) && (row->cCol < row->nCol)) {
                CSeqCol *r1 = row->col + row->cCol;
                if(r1->tag == min_tag) {
                  row->cCol++;
                }
              }
            }
            current+=max_width;
          }
        }
      }
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
            l1 = lab->col + b; /* if a fixed label is present, 
                                  get the final offset from the residue line */
            if(l1->stop) 
              l1->offset = r1->offset;
          }
        }
        lab = NULL;
      }
    }
    
  }

  /* THIRD PASS: fill in labels, based on actual residue spacing */

  if(nRow&&(codes!=4)) {
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
        int n_skipped = 0;
        int last_resv = -1;
        AtomInfoType *last_ai = NULL;
        ObjectMolecule *obj;
        AtomInfoType *ai;
        row = lab+1;
        nCol = row->nCol;
        obj = row->obj;
        div = SettingGet_i(G,obj->Obj.Setting,NULL,cSetting_seq_view_label_spacing);
        sub = SettingGet_i(G,obj->Obj.Setting,NULL,cSetting_seq_view_label_start);
        for(b=0;b<nCol;b++) {
          CSeqCol *r1 = row->col + b;
          CSeqCol *l1 = lab->col + b;

          ai = NULL;
          if(r1->atom_at) {
            atom_list = row->atom_lists + r1->atom_at;
            if(*atom_list>=0)
              ai = obj->AtomInfo + (*atom_list); /* get first atom in list */
          }
          if(l1->stop) {/* if label is already present, just line it up */
            l1->offset = r1->offset;
          } else if((r1->offset >= next_open)&&ai) {
            if((div>1)&&(codes!=2)) {
              if(! ((ai->resv-sub) % div))
                draw_it = true;
              else
                draw_it = false;
            } else {
              draw_it = true;
            }
            if(ai->resv!=(last_resv+1)) /* gap in sequence?  then draw label ASAP */
              draw_it = true;
            if(n_skipped >= (div+div)) /* don't skip too many without a label! */
              draw_it = true;

            if(AtomInfoSameResidueP(G,last_ai,ai)) /* don't ever draw a residue label twice */
              draw_it = false;

            if(draw_it) {
              n_skipped = 0;
              last_ai = ai;
              l1->start = lab->len;
              if(codes==2) {
                UtilConcatVLA(&lab->txt,&lab->len,ai->resn);
                UtilConcatVLA(&lab->txt,&lab->len,"`");
              }
                
              UtilConcatVLA(&lab->txt,&lab->len,ai->resi);
              l1->stop = lab->len;
              st_len = l1->stop - l1->start + 1;
              l1->offset = r1->offset;
              next_open = r1->offset + st_len;
              
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
            } else
              n_skipped++;
          }

          if(ai)
            last_resv = ai->resv;
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
  G->Seeker->handler.fClick = SeekerClick;
  G->Seeker->handler.fRelease = SeekerRelease;
  G->Seeker->handler.fDrag = SeekerDrag;
  G->Seeker->handler.fRefresh = SeekerRefresh;
  SeqSetRowVLA(G,row_vla,nRow);
  SeqSetHandler(G,&G->Seeker->handler);
}

int SeekerInit(PyMOLGlobals *G)
{
  register CSeeker *I=NULL;
  if( (I=(G->Seeker=Calloc(CSeeker,1)))) {
    
    UtilZeroMem(I,sizeof(CSeeker));
    I->drag_row = -1;
    I->LastClickTime = UtilGetSeconds(G) - 1.0F;
    return 1;
  } else {
    return 0;
  }
}


void SeekerFree(PyMOLGlobals *G)
{
  FreeP(G->Seeker);
}


