/* 
A* -------------------------------------------------------------------
B* This fil econtains source code for the PyMOL computer program
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

#include"Version.h"
#include"main.h"
#include"Base.h"
#include"OOMac.h"
#include"Executive.h"
#include"ObjectMesh.h"
#include"ObjectDist.h"
#include"ObjectSurface.h"
#include"ObjectSlice.h"
#include"ListMacros.h"
#include"Ortho.h"
#include"Scene.h"
#include"Selector.h"
#include"Vector.h"
#include"Color.h"
#include"Setting.h"
#include"Matrix.h"
#include"P.h"
#include"PConv.h"
#include"Match.h"
#include"ObjectCGO.h"
#include"Util.h"
#include"Wizard.h"
#include"ScrollBar.h"
#include"Movie.h"
#include"ObjectGadgetRamp.h"
#include"SculptCache.h"
#include"Control.h"
#include"Menu.h"
#include"Map.h"
#include"Editor.h"
#include"RepDot.h"
#include"Seq.h"
#include"Text.h"
#include"PyMOL.h"

#define cExecObject 0
#define cExecSelection 1
#define cExecAll 2

#define cTempRectSele "_rect"
#define cLeftButSele "lb"
#define cIndicateSele "indicate"

typedef struct SpecRec {
  int type;
  WordType  name; /*only used for selections*/
  struct CObject *obj;  
  struct SpecRec *next;
  int repOn[cRepCnt];
  int visible;
  int sele_color;
  int hilight;
  int previous;
} SpecRec; /* specification record (a line in the executive window) */

struct _CExecutive {
  Block *Block;
  SpecRec *Spec;
  int Width,Height,HowFarDown;
  int ScrollBarActive;
  int NSkip;
  struct CScrollBar *ScrollBar;
  CObject *LastEdited;
  int DragMode;
  int Pressed,Over,OldVisibility,ToggleMode;
  SpecRec *LastChanged;
  int ReorderFlag;
  OrthoLineType ReorderLog;

  int oldPX,oldPY,oldWidth,oldHeight,sizeFlag;
};

static SpecRec *ExecutiveFindSpec(PyMOLGlobals *G,char *name);
static int ExecutiveDrag(Block *block,int x,int y,int mod);
static void ExecutiveSpecSetVisibility(PyMOLGlobals *G,SpecRec *rec,
                                       int new_vis,int mod);
void ExecutiveObjMolSeleOp(PyMOLGlobals *G,int sele,ObjectMoleculeOpRec *op);

static int ExecutiveCountNames(PyMOLGlobals *G)
{
  int count=0;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next))
	 count++;
  
  return(count);
}

static int ReorderOrderFn(PyMOLGlobals *G,SpecRec **rec,int l,int r)
{
  return (WordCompare(G,rec[l]->name,rec[r]->name,true)<=0);
}

int ExecutiveOrder(PyMOLGlobals *G, char *s1, int sort,int location)
{
  register CExecutive *I = G->Executive;
  int ok=true;
  CWordList *word_list = WordListNew(G,s1);
  int n_names = ExecutiveCountNames(G);
  if(n_names) {
    SpecRec **list,**subset,**sorted;
    int *index = NULL;
    int n_sel;
    int source_row = -1;
    list = Alloc(SpecRec*,n_names);
    subset = Calloc(SpecRec*,n_names);
    sorted = Calloc(SpecRec*,n_names);
    index = Alloc(int,n_names);
    if(list&&subset) {
      /* create an array of current names */
      {
        SpecRec *rec = NULL;
        int a = 0;
        /* copy all names into array */
        while(ListIterate(I->Spec,rec,next)) {
          list[a] = rec;
          a++;
        }
        for(a=0;a<n_names;a++) {
          list[a]->next = NULL;
        }
      } 
      /* transfer matching names to the subset array */
      {
        int a;
        int entry;
        int min_entry = word_list->n_word;
        for(a=n_names-1;a>0;a--) { /* skipping zeroth */
          entry = WordListMatch(G,word_list, list[a]->name, true);
          if(entry>=0) { /* append onto the new list */
            list[a]->next = subset[entry];
            subset[entry] = list[a];
            list[a] = NULL;
            if(entry<=min_entry) {
              source_row = a; /* takes the earliest first match */
              min_entry = entry;
            }
          }
        }
        if(word_list->n_word && WordMatchExact(G,word_list->start[0],cKeywordAll,true))
          location=-1; /* set to top if "all" is first in list */
      }
      /* expand the selected entries */
      {
        SpecRec *rec,*last;
        int b;
        n_sel = 0;
        for(b=0;b<word_list->n_word;b++) {
          rec=subset[b];
          while(rec) {
            sorted[n_sel++] = rec;
            last = rec;
            rec = rec->next;
            last->next = NULL;
          }
        }
      }
      /* sort the selected entries, if requested */
      if(sort) {
        UtilCopyMem(subset,sorted,sizeof(SpecRec*)*n_sel);
        {
          int a;
          UtilSortIndexGlobals(G,n_sel,subset,index,
                               (UtilOrderFnGlobals*)ReorderOrderFn);
          for(a=0;a<n_sel;a++) {
            sorted[a] = subset[index[a]];
          }
        }
      }
      /* reassemble the list using the new order */
      {
        SpecRec *spec= NULL;
        SpecRec *last= NULL;
        int a,b;
        int flag;
        for(a=0;a<n_names;a++) {
          flag=false;
          if(sorted) { /* not yet added */
            switch(location) {
            case -1: /* top */
              if(a==1) flag=true;
              break;
            case 0: 
              if(source_row>=0) {
                if(a==source_row)
                  flag=true;
              } else if(!list[a]) 
                flag=true;
              break;
            }
          }
          if(flag) {
            for(b=0;b<n_sel;b++) {
              if(sorted[b]) {
                if(last)
                  last->next = sorted[b];
                last = sorted[b];
                if(!spec)
                  spec = last;
              }
            }
            FreeP(sorted);
          }
          if(list[a]) {
            if(last) 
              last->next = list[a];
            last = list[a];
            if(!spec)
              spec=last;
          }
        }
        if(sorted) { /* still not yet readded? */
          for(b=0;b<n_sel;b++) {
            if(sorted[b]) {
              if(last)
                last->next = sorted[b];
              last = sorted[b];
              if(!spec)
                spec = last;
            }
          }
        }
        I->Spec=spec;
        OrthoDirty(G);
      }
      
      FreeP(index);
      FreeP(sorted);
      FreeP(list);
      FreeP(subset);
    }
  }
  WordListFreeP(word_list);
  return(ok);
}

ObjectMolecule **ExecutiveGetObjectMoleculeVLA(PyMOLGlobals *G,char *sele)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  ObjectMolecule **result = NULL;
  int s1=SelectorIndexByName(G,sele);  
  if(s1>=0) {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_GetObjects;
    op.obj1VLA=(ObjectMolecule**)VLAlloc(CObject*,10);
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G,s1,&op);
    result = (ObjectMolecule**)op.obj1VLA;
    VLASize(result,ObjectMolecule*,op.i1);
  }
  return result;
#endif
}

/* #define ExecLineHeight 18 */
#define ExecClickMargin 2
#define ExecTopMargin 0
#define ExecToggleMargin 2
#define ExecLeftMargin 1
#define ExecRightMargin 0
#define ExecToggleWidth 17
#define ExecToggleSize 16
#define ExecToggleTextShift 4

#define ExecOpCnt 5

typedef struct { 
  M4XAnnoType m4x;
  ObjectMolecule *obj;
} ProcPDBRec;

int ExecutivePop(PyMOLGlobals *G,char *target,char *source,int quiet)
{
  int ok = true;
  int src;
  int result = 0;

  ExecutiveDelete(G,target);
  if(ExecutiveFindObjectMoleculeByName(G,source)) {
    ok=false;
    PRINTFB(G,FB_Executive,FB_Errors)
      " Pop-Error: source selection '%s' can't be an object.\n",source
      ENDFB(G);
    
  } else {
    src = SelectorIndexByName(G,source);
    if(src<0)
      ok=false;
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        " Pop-Error: invalid source selection name '%s'\n",source
        ENDFB(G);
    } else {
      ObjectMoleculeOpRec op;
      
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_Pop;
      SelectorCreateEmpty(G,target);
      op.i1 = SelectorIndexByName(G,target);
      op.i2 = 1;
      op.i3 = 0;
      ExecutiveObjMolSeleOp(G,src,&op);
      result = op.i3;
    }
  }
  if(!result) ExecutiveDelete(G,target);
  if(!ok)
    return -1;
  else
    return result;
}

int ExecutiveGetActiveSele(PyMOLGlobals *G)
{

  char name[ObjNameMax];
  if(ExecutiveGetActiveSeleName(G,name,false))
    return SelectorIndexByName(G,name);
  else
    return -1;

}

int ExecutiveGetActiveSeleName(PyMOLGlobals *G,char *name, int create_new)
{
  /* TODO: cache/optimize to avoid table scan */

  int result=false;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection)
      if(rec->visible) {
        strcpy(name,rec->name);
        result = true;
      }
  }
  if((!result)&&create_new) {
    int sel_num = SettingGetGlobal_i(G,cSetting_sel_counter) + 1;

    SettingSetGlobal_i(G,cSetting_sel_counter,sel_num);
    sprintf(name,"sel%02d",sel_num);
    SelectorCreateEmpty(G,name);
  }
  return result;
}


int ExecutiveFixChemistry(PyMOLGlobals *G,char *s1,char *s2,int invalidate,int quiet)
{
  int sele1=SelectorIndexByName(G,s1);
  int sele2=SelectorIndexByName(G,s2);
  int ok=true;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;

  if((sele1>=0)&&(sele2>=0)) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule) {
          ObjectMoleculeFixChemistry((ObjectMolecule*)rec->obj,sele1,sele2,invalidate);
        }
    }
  }
  return ok;
}

int ExecutiveGetObjectColorIndex(PyMOLGlobals *G,char *name)
{
  int result = -1;
  CObject *obj = NULL;
  obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    result=obj->Color;
  }
  return(result);
}

int ExecutiveGetAtomVertex(PyMOLGlobals *G,char *s1,int state,int index,float *v)
{
  int ok=false;
  int sele1 = SelectorIndexByName(G,s1);

  if(sele1>=0) {
    ok=SelectorGetSingleAtomVertex(G,sele1,state,v);
  }
  return ok;
}

int ExecutiveSetName(PyMOLGlobals *G,char *old_name, char *new_name)
{
  int ok=true;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;
  int found = false;
  if(!new_name[0]) 
    ok=false;
  else
    
    while(ListIterate(I->Spec,rec,next)) {
      if(found)
        break;
      switch(rec->type) {
      case cExecObject:
        if(WordMatchExact(G,rec->obj->Name,old_name,true)) {
          ObjectSetName(rec->obj,new_name);
          UtilNCopy(rec->name,rec->obj->Name,ObjNameMax);
          if(rec->obj->type == cObjectMolecule) {
            /*
              SelectorDelete(G,old_name);
              ExecutiveUpdateObjectSelection(G,rec->obj);
            */
            SelectorSetName(G,new_name, old_name);
            SceneDirty(G);
            SeqChanged(G);
            found = true;
          }
        }
        break;
      case cExecSelection:
        if(WordMatchExact(G,rec->name,old_name,true)) {
          if(SelectorSetName(G,new_name, old_name)) {
            UtilNCopy(rec->name,new_name,ObjNameMax);
            found = true;
            OrthoDirty(G);
          }
        }
        break;
      }
    }
  if(!found)
    ok=false;
  return ok; 
}

void ExecutiveLoadMOL2(PyMOLGlobals *G,CObject *origObj,char *fname,
                       char *oname, int frame, int discrete,int finish,
                       OrthoLineType buf,int multiplex,int quiet,
                       int is_string)
{
  int ok=true;
  FILE *f;
  long size;
  char *buffer=NULL,*p;
  CObject *obj;
  char new_name[ObjNameMax] = "";
  char *next_entry = NULL;
  int repeat_flag = true;
  int n_processed = 0;

  if(is_string) {
    buffer=fname;
  } else {
    f=fopen(fname,"rb");
    
    if(!f) {
      PRINTFB(G,FB_Executive,FB_Errors)
        " ExecutiveLoadMOL2-ERROR: Unable to open file '%s'\n",fname
        ENDFB(G);
      ok=false;
    } else
      {
        PRINTFB(G,FB_Executive,FB_Blather)
          " ExecutiveLoadMOL2: Loading from %s.\n",fname
          ENDFB(G);
        
        fseek(f,0,SEEK_END);
        size=ftell(f);
        fseek(f,0,SEEK_SET);
        
        buffer=(char*)mmalloc(size+255);
        ErrChkPtr(G,buffer);
        p=buffer;
        fseek(f,0,SEEK_SET);
        fread(p,size,1,f);
        p[size]=0;
        fclose(f);
      }
  }

  while(repeat_flag&&ok) {
    char *start_at = buffer;
    int is_repeat_pass = false;
    int eff_frame = frame;

    if(next_entry) {
      start_at = next_entry;
      is_repeat_pass = true;
    }

    PRINTFD(G,FB_CCmd) " ExecutiveLoadMOL2-DEBUG: loading PDB\n" ENDFD;

    repeat_flag=false;
    next_entry = NULL;
    if(!origObj) {

      new_name[0]=0;
      obj=(CObject*)ObjectMoleculeReadMOL2Str(G,(ObjectMolecule*)origObj,
                                              start_at,eff_frame,discrete,
                                              quiet,multiplex,new_name,
                                              &next_entry);
      if(obj) {
        if(next_entry) { /* NOTE: if set, assume that multiple PDBs are present in the file */
          repeat_flag=true;
        }

        /* assign the name (if necessary) */

        if(next_entry||is_repeat_pass) {
          if(new_name[0]==0) {
            sprintf(new_name,"%s_%d",oname,n_processed+1);
          }
          ObjectSetName(obj,new_name); /* from PDB */
          ExecutiveDelete(G,new_name); /* just in case */
        } else {
          ObjectSetName(obj,oname); /* from filename/parameter */
        }

        if(obj) {
          ExecutiveManageObject(G,obj,true,true);
          if(eff_frame<0)
            eff_frame = ((ObjectMolecule*)obj)->NCSet-1;
          if(multiplex>0) {
            sprintf(buf," CmdLoad: loaded %d objects.\n",n_processed);
          } else {
            if(!is_string)
              sprintf(buf," CmdLoad: \"%s\" loaded as \"%s\".\n",
                      fname,oname);
            else
              sprintf(buf," CmdLoad: MOL2-string loaded into object \"%s\", state %d.\n",
                      oname,eff_frame+1);
          }
        }
      }
    } else {

      ObjectMoleculeReadMOL2Str(G,(ObjectMolecule*)origObj,
                                start_at,eff_frame,discrete,
                                quiet,multiplex,new_name,
                                &next_entry);

      if(finish)
        ExecutiveUpdateObjectSelection(G,origObj);
      if(eff_frame<0)
        eff_frame = ((ObjectMolecule*)origObj)->NCSet-1;
      if(!is_string) 
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,eff_frame+1);
      else
        sprintf(buf," CmdLoad: PDB-string appended into object \"%s\", state %d.\n",
                oname,eff_frame+1);
      obj = origObj;
    }

    if(obj) {
      n_processed++;
    }
  }

  if((!is_string)&&buffer) {
    mfree(buffer);
  }

}


void ExecutiveProcessPDBFile(PyMOLGlobals *G,CObject *origObj,char *fname,
                             char *oname, int frame, int discrete,int finish,
                             OrthoLineType buf,PDBInfoRec *pdb_info,int quiet,
                             int is_string)
{
  int ok=true;
  FILE *f;
  long size;
  char *buffer=NULL,*p;
  CObject *obj;
  char pdb_name[ObjNameMax] = "";
  char *next_pdb = NULL;
  int repeat_flag = true;
  ProcPDBRec *processed= NULL;
  int n_processed = 0;
  int m4x_mode = 0; /* 0 = annotate, 1 = alignment */
  ProcPDBRec *target_rec = NULL;
  char nbrhood_sele[] = "m4x_nearby";
  ProcPDBRec *current = NULL;
  PDBInfoRec pdb_info_rec;

  if(!pdb_info) {
    UtilZeroMem(&pdb_info_rec,sizeof(PDBInfoRec));
    pdb_info=&pdb_info_rec;
  }

  if(is_string) {
    buffer=fname;
  } else {
    f=fopen(fname,"rb");
    
    if(!f) {
      PRINTFB(G,FB_ObjectMolecule,FB_Errors)
        "ObjectMolecule-ERROR: Unable to open file '%s'\n",fname
        ENDFB(G);
      ok=false;
    } else
      {
        PRINTFB(G,FB_ObjectMolecule,FB_Blather)
          " ObjectMoleculeLoadPDBFile: Loading from %s.\n",fname
          ENDFB(G);
        
        fseek(f,0,SEEK_END);
        size=ftell(f);
        fseek(f,0,SEEK_SET);
        
        buffer=(char*)mmalloc(size+255);
        ErrChkPtr(G,buffer);
        p=buffer;
        fseek(f,0,SEEK_SET);
        fread(p,size,1,f);
        p[size]=0;
        fclose(f);
      }
  }

  if(ok) {
    processed = VLACalloc(ProcPDBRec,10);
  }
  while(repeat_flag&&ok) {
    char *start_at = buffer;
    int is_repeat_pass = false;
    int eff_frame = frame;
    CObject *tmpObj;

    VLACheck(processed,ProcPDBRec,n_processed);
    current = processed + n_processed;

    PRINTFD(G,FB_CCmd) " ExecutiveProcessPDBFile-DEBUG: loading PDB\n" ENDFD;

    if(next_pdb) {
      start_at = next_pdb;
      is_repeat_pass = true;
    }

    M4XAnnoInit(&current->m4x);

    repeat_flag=false;
    next_pdb = NULL;
    if(!origObj) {

      pdb_name[0]=0;
      obj=(CObject*)ObjectMoleculeReadPDBStr(G,(ObjectMolecule*)origObj,
                                             start_at,eff_frame,discrete,
                                             &current->m4x,pdb_name,
                                             &next_pdb,pdb_info,quiet);

      if(obj) {
        if(next_pdb) { /* NOTE: if set, assume that multiple PDBs are present in the file */
          repeat_flag=true;
        }

        if(current->m4x.xname_flag) { /* USER XNAME trumps the PDB Header name */                 
          ObjectSetName(obj,current->m4x.xname); /* from PDB */
          if( (tmpObj = ExecutiveFindObjectByName(G,obj->Name))) {
            if(tmpObj->type != cObjectMolecule)
              ExecutiveDelete(G,current->m4x.xname); /* just in case */
            else
              {
                if(is_repeat_pass) {  /* this is a workaround for when PLANET accidentally duplicates the target */
                  {
                    int a;
                    for(a=0; a<n_processed; a++) {
                      ProcPDBRec *cur = processed + a;
                      if(cur->obj == (ObjectMolecule*)tmpObj) {
                        cur->m4x.invisible = false;
                      }
                    }
                  }
                  ObjectMoleculeFree((ObjectMolecule*)obj);
                  obj=NULL;
                }
              }
          }
        } else if(next_pdb) {
          if(pdb_name[0]==0) {
            sprintf(pdb_name,"%s_%d",oname,n_processed+1);
          }
          ObjectSetName(obj,pdb_name); /* from PDB */
          ExecutiveDelete(G,pdb_name); /* just in case */
        } else {
          if(is_repeat_pass) {
            if(pdb_name[0]==0) {
              sprintf(pdb_name,"%s_%d",oname,n_processed+1);
            }
            ObjectSetName(obj,pdb_name); /* from PDB */
            ExecutiveDelete(G,pdb_name); /* just in case */
          } else {
            ObjectSetName(obj,oname); /* from filename/parameter */
          }
        }

        if(obj) {
          ExecutiveManageObject(G,obj,true,true);
          if(eff_frame<0)
            eff_frame = ((ObjectMolecule*)obj)->NCSet-1;
          if(buf) {
            if(!is_string)
              sprintf(buf," CmdLoad: \"%s\" loaded as \"%s\".\n",
                      fname,oname);
            else
              sprintf(buf," CmdLoad: PDB-string loaded into object \"%s\", state %d.\n",
                      oname,eff_frame+1);
          }
            
        }
      }
    } else {
      ObjectMoleculeReadPDBStr(G,(ObjectMolecule*)origObj,
                               start_at,eff_frame,discrete,&current->m4x,
                               pdb_name,&next_pdb,pdb_info,quiet);
      if(finish)
        ExecutiveUpdateObjectSelection(G,origObj);
      if(eff_frame<0)
        eff_frame = ((ObjectMolecule*)origObj)->NCSet-1;
      if(buf) {
        if(!is_string) 
          sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                  fname,oname,eff_frame+1);
        else
          sprintf(buf," CmdLoad: PDB-string appended into object \"%s\", state %d.\n",
                  oname,eff_frame+1);
      }
      obj = origObj;
    }

    if(obj&&current) {
      current->obj = (ObjectMolecule*)obj;
      n_processed++;
    }
  }

  /* BEGIN METAPHORICS ANNOTATION AND ALIGNMENT CODE */

  /* sanity check -- make sure all objects are present */
  if(ok&&n_processed) {
    int a;
    for(a=0; a<n_processed; a++) {
      ProcPDBRec *current = processed + a;
      if(!ExecutiveValidateObjectPtr(G,(CObject*)current->obj,cObjectMolecule)) {
        PRINTFB(G,FB_Executive,FB_Errors)
          " Error: Missing object! possible invalid/corrupt file.\n"
          ENDFB(G);
        ok=false;
        break;
      }
    }
  }
  
  if(ok&&n_processed) { /* first, perform any Metaphorics alignment */
    /* is there a target structure? */
    {
      int a;
      for(a=0; a<n_processed; a++) {
        ProcPDBRec *current = processed + a;
        M4XAnnoType *m4x = &current->m4x;

        if(m4x->annotated_flag && m4x->align) {
          if(WordMatchExact(G,current->obj->Obj.Name,m4x->align->target,true)) {
            target_rec = current;
            break;
          }
        }
      }
    }
    if(target_rec) { /* there is a target.. */
      
      /* first, convert all IDs to genuine atom indices */

      {
        int a;
        for(a=0; a<n_processed; a++) {
          ProcPDBRec *current = processed + a;
          M4XAnnoType *m4x = &current->m4x;
          if(m4x->align) {
            ObjectMoleculeConvertIDsToIndices(current->obj, 
                                              m4x->align->id_at_point,
                                              m4x->align->n_point);
          }
        }
      }

      /* next, peform the alignments against the target */

      {
        int a;
        char aligned_name[] = "m4x_aligned";
        char tmp_sele[ObjNameMax*3];

        SelectorCreateEmpty(G,"m4x_aligned");
                
                
        for(a=0; a<n_processed; a++) {
          ProcPDBRec *current = processed + a;
          if( current != target_rec ) {
            M4XAnnoType *m4x = &current->m4x;
            if(m4x->align) {
              ObjMolPairwise pairwise;

              m4x_mode = 1; /* performing an alignment */

              ObjMolPairwiseInit(&pairwise);
              pairwise.trg_obj = target_rec->obj;
              pairwise.mbl_obj = current->obj;
              
              /* create ordered matches */

              {
                M4XAlignType *trg_align = target_rec->m4x.align;
                M4XAlignType *mbl_align = m4x->align;

                int n_point;
                int a;
                  
                n_point = trg_align->n_point;
                if(n_point > mbl_align->n_point) n_point = mbl_align->n_point;
                
                VLACheck(pairwise.trg_vla,int,n_point);
                VLACheck(pairwise.mbl_vla,int,n_point);
                
                for(a=0;a<n_point;a++) {
                  int trg_index = trg_align->id_at_point[a];
                  int mbl_index = mbl_align->id_at_point[a];
                  if((trg_index>=0)&&(mbl_index>=0)&&
                     (mbl_align->fitness[a]>=0.0F)) {
                    pairwise.trg_vla[pairwise.n_pair] = trg_index;
                    pairwise.mbl_vla[pairwise.n_pair] = mbl_index;
                    pairwise.n_pair++;
                  }
                }

                {
                  char trg_sele[ObjNameMax],mbl_sele[ObjNameMax];
                  char align_name[ObjNameMax];
                  SelectorGetUniqueTmpName(G,trg_sele);
                  SelectorGetUniqueTmpName(G,mbl_sele);
                  
                  SelectorCreateOrderedFromObjectIndices(G,trg_sele,pairwise.trg_obj,
                                                         pairwise.trg_vla,pairwise.n_pair);
                  SelectorCreateOrderedFromObjectIndices(G,mbl_sele,pairwise.mbl_obj,
                                                         pairwise.mbl_vla,pairwise.n_pair);
                  
                  sprintf(align_name,"%s_%s_alignment",
                          pairwise.trg_obj->Obj.Name,
                          pairwise.mbl_obj->Obj.Name);

                  ExecutiveRMS(G,mbl_sele, trg_sele, 2, 0.0F, 0, 0, align_name, 0, 0, true);
                  ExecutiveColor(G,align_name,"white",0,true);
                  if(target_rec->m4x.invisible) 
                    sprintf(tmp_sele, "(%s) | (%s)",aligned_name,mbl_sele);
                  else 
                    sprintf(tmp_sele, "(%s) | (%s) | (%s)",aligned_name, trg_sele, mbl_sele);
                  SelectorCreateSimple(G,aligned_name, tmp_sele);
                  /*ExecutiveDelete(G,trg_sele);
                    ExecutiveDelete(G,mbl_sele);*/
                }

              }
              ObjMolPairwisePurge(&pairwise);
            }
          }
        }
        sprintf(tmp_sele, "bychain %s", aligned_name);
        SelectorCreateSimple(G,aligned_name, tmp_sele);
        sprintf(tmp_sele, "byres (%s expand 3.5)", aligned_name);
        SelectorCreateSimple(G,nbrhood_sele, tmp_sele);
      }
    }
  }
  
  if(ok&&n_processed) { /* next, perform any and all Metaphorics annotations */
      int a;
      int nbr_sele = SelectorIndexByName(G,nbrhood_sele);
      
      for(a=0; a<n_processed; a++) {
        ProcPDBRec *current = processed + a;
        if(current->m4x.annotated_flag) {
          char annotate_script[] = "@$PYMOL_SCRIPTS/metaphorics/annotate.pml";
          char align_script[] = "@$PYMOL_SCRIPTS/metaphorics/alignment.pml";
          char *script_file = NULL;
          
          if(a==(n_processed-1)) { /* for multi-PDB files, don't execute script until after the last file */
            switch(m4x_mode) {
            case 0:
              script_file = annotate_script;
              break;
            case 1:
              script_file = align_script;
              break;
            }
          }

          if((current!=target_rec)||(!current->m4x.invisible)) /* suppress annotations if target invisible */
            ObjectMoleculeM4XAnnotate(current->obj,&current->m4x,script_file,(m4x_mode==1),
                                      nbr_sele);
          M4XAnnoPurge(&current->m4x);
        }
      }
  }
  /* END METAPHORICS ANNOTATION AND ALIGNMENT CODE */
  
  VLAFreeP(processed);
  if((!is_string)&&buffer) {
    mfree(buffer);
  }
  
}


int  ExecutiveAssignSS(PyMOLGlobals *G,char *target,int state,char *context,int preserve,int quiet)
{
  int sele0=-1;
  int sele1=-1;
  int ok = false;
  sele0 = SelectorIndexByName(G,target);
  if(sele0>=0) {
    if(!context[0]) {
      sele1=sele0;
    } else {
      sele1 = SelectorIndexByName(G,context);
    }
    if(sele1>=0) {
      ok =  SelectorAssignSS(G,sele0,sele1,state,preserve,quiet);
    }
  }
  return(ok);
}

PyObject *ExecutiveGetVisAsPyDict(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result=NULL,*list,*repList;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int a;
  int n_vis;
  result = PyDict_New();
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->name[0]!='_') {
      list = PyList_New(4);
      PyList_SetItem(list,0,PyInt_FromLong(rec->visible));

      /* all executive entries have repOn */
      n_vis=0;
      for(a=0;a<cRepCnt;a++) {
        if(rec->repOn[a])
          n_vis++;
      }
      repList = PyList_New(n_vis);
      n_vis=0;
      for(a=0;a<cRepCnt;a++) {
        if(rec->repOn[a]) {
          PyList_SetItem(repList,n_vis,PyInt_FromLong(a));
          n_vis++;
        }
      }
      PyList_SetItem(list,1,repList);
      
      if(rec->type!=cExecObject) {
        Py_INCREF(Py_None);
        PyList_SetItem(list,2,Py_None);
        Py_INCREF(Py_None);
        PyList_SetItem(list,3,Py_None);
      } else { 
        /* objects have their own visib list too */
        n_vis=0;
        for(a=0;a<cRepCnt;a++) {
          if(rec->obj->RepVis[a])
            n_vis++;
        }
        repList = PyList_New(n_vis);
        n_vis=0;
        for(a=0;a<cRepCnt;a++) {
          if(rec->obj->RepVis[a]) {
            PyList_SetItem(repList,n_vis,PyInt_FromLong(a));
            n_vis++;
          }
        }
        PyList_SetItem(list,2,repList);
        PyList_SetItem(list,3,PyInt_FromLong(rec->obj->Color));
      }

      PyDict_SetItemString(result,rec->name,list);
      Py_DECREF(list);
    }
  }
  return(result);
#endif
}

int ExecutiveSetVisFromPyDict(PyMOLGlobals *G,PyObject *dict)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  WordType name;
  PyObject *key,*list,*col;
  PyObject *vis_list = NULL;
  int pos = 0;
  SpecRec *rec;
  int n_vis;
  int rep;
  int a;
  int ll=0;
  if(ok) ok=(dict!=NULL);
  if(ok) ok=PyDict_Check(dict);
  if(ok) {

    SceneObjectDel(G,NULL);
    while (PyDict_Next(dict, &pos, &key, &list)) {
      if(!PConvPyStrToStr(key,name,sizeof(WordType))) {
        ok=false;
      } else {
        
        rec = ExecutiveFindSpec(G,name);
        if(rec) {
          if(ok) ok = (list!=NULL);
          if(ok) ok = PyList_Check(list);
          if(ok) ll = PyList_Size(list);
          if(ok) ok = (ll>=2);
          if(ok) ok = PConvPyObjectToInt(PyList_GetItem(list,0),&rec->visible);
          if(ok) { /* rec visibility */
            vis_list = PyList_GetItem(list,1);
            if(ok) ok = (vis_list!=NULL);
            if(ok) ok = PyList_Check(vis_list);
            if(ok) {
              n_vis = PyList_Size(vis_list);
              for(a=0;a<cRepCnt;a++)
                rec->repOn[a]=false;
              for(a=0;a<n_vis;a++) {
                if(PConvPyObjectToInt(PyList_GetItem(vis_list,a),&rep)) {
                  if((rep>=0)&&(rep<cRepCnt))
                    rec->repOn[rep]=true;
                }
              }
            }
          }

          if(ok&&(rec->type==cExecObject)) { /* object properties */

            if(ll>2) { /* object visibility */
              vis_list = PyList_GetItem(list,2);
              if(ok) ok = (vis_list!=NULL);
              if(ok) if(PyList_Check(vis_list)) {
                n_vis = PyList_Size(vis_list);
                for(a=0;a<cRepCnt;a++)
                  rec->obj->RepVis[a]=false;
                for(a=0;a<n_vis;a++) {
                  if(PConvPyObjectToInt(PyList_GetItem(vis_list,a),&rep)) {
                    if((rep>=0)&&(rep<cRepCnt))
                      rec->obj->RepVis[rep]=true;
                  }
                }
              }
            }
            if(ll>3) { /* object color */
              col = PyList_GetItem(list,3);
              if(ok) ok = (col!=NULL);
              if(ok) if(PyInt_Check(col)) {
                ok = PConvPyObjectToInt(col,&rec->obj->Color);
                if(rec->obj->fInvalidate)
                  rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,-1);
              }
            }
          }
          if(rec->visible&&(rec->type==cExecObject))
            SceneObjectAdd(G,rec->obj); 
        }
      }
    }
  }
  return ok;
#endif
}

int ExecutiveIsolevel(PyMOLGlobals *G,char *name,float level,int state)
{
  int ok =true;
  CObject *obj;
  obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    switch(obj->type) {
    case cObjectMesh:
      ObjectMeshSetLevel((ObjectMesh*)obj,level,state);
        SceneChanged(G);
      break;
    case cObjectSurface:
      break;
    default:
      ok=false;
      PRINTFB(G,FB_Executive,FB_Errors)
        " Isolevel-Error: object \"%s\" is of wrong type.",name
        ENDFB(G);
      break;
    }
  }
  return(ok);

}

int ExecutiveSpectrum(PyMOLGlobals *G,char *s1,char *expr,float min,float max,int first,int last,
                      char *prefix,int digits,int byres,int quiet,
                      float *min_ret,float *max_ret)
{
  int ok=true;
  int sele1;
  int n_color,n_atom;
  ObjectMoleculeOpRec op;
  WordType buffer;
  int *color_index = NULL;
  float *value = NULL;
  int a,b;
  char pat[] = "%0Xd";
  int pref_len;
  char *at;
  float range;

  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {

    if(digits>9) digits = 9;
    pat[2]=('0'+digits);
    UtilNCopy(buffer,prefix,sizeof(WordType)-digits);
    
    pref_len = strlen(prefix);
    at = buffer+pref_len;

    n_color = abs(first-last)+1;
    if(n_color) {
      color_index = Alloc(int,n_color);
      for(a=0;a<n_color;a++) {
        b = first + ((last-first)*a)/(n_color-1);
        sprintf(at,pat,b);
        color_index[a] = ColorGetIndex(G,buffer);
      }
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_CountAtoms;
      op.i1=0;
      ExecutiveObjMolSeleOp(G,sele1,&op);
      n_atom=op.i1;
      
      if(n_atom) {
        value = Calloc(float,n_atom);
        
        if(WordMatch(G,"count",expr,true)) {
          for(a=0;a<n_atom;a++) {
            value[a]=(float)a+1;
          }
        } else if(WordMatch(G,"b",expr,true)) {
          op.code = OMOP_GetBFactors;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(G,sele1,&op);
        } else if(WordMatch(G,"q",expr,true)) {
          op.code = OMOP_GetOccupancies;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(G,sele1,&op);
        } else if(WordMatch(G,"pc",expr,true)) {
          op.code = OMOP_GetPartialCharges;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(G,sele1,&op);
        }

        if(max<min) {
          max = value[0];
          min = value[0];
          for(a=1;a<n_atom;a++) {
            if(value[a]<min) min=value[a];
            if(value[a]>max) max=value[a];
          }
        }
        range = max-min;

        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Spectrum: range (%8.5f to %8.5f).\n"
            ,min,max
            ENDFB(G);
        }
        if(range==0.0F)
          range = 1.0F;
        *min_ret = min;
        *max_ret = max;


        op.code = OMOP_Spectrum;
        op.i1 = n_color-1;
        op.i2 = n_atom;
        op.i3 = 0;
        op.i4 = byres;
        op.ii1 = color_index;
        op.ff1 = value;
        op.f1 = min;
        op.f2 = range;
        
        ExecutiveObjMolSeleOp(G,sele1,&op);

        op.code=OMOP_INVA;
        op.i1=cRepAll; 
        op.i2=cRepInvColor;
        ExecutiveObjMolSeleOp(G,sele1,&op);

      }
    }

    FreeP(color_index);
    FreeP(value);
  }
  return(ok);
}

char *ExecutiveGetChains(PyMOLGlobals *G,char *sele,int state,int *null_chain)
{
  int sele1;
  char *result = NULL;
  int chains[256];
  int a,c;
  ObjectMoleculeOpRec op;

  sele1 = SelectorIndexByName(G,sele);
  if(sele1>=0) {

    for(a=0;a<256;a++) {
      chains[a]=0;
    }
    ObjectMoleculeOpRecInit(&op);
    op.code=OMOP_GetChains;
    op.ii1 = chains;
    op.i1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op);
    c=0;
    for(a=1;a<256;a++) {
      if(chains[a]) c++;
    }
    result = Calloc(char,c+1);
    if(result) {
      c=0;
      *null_chain = chains[0];
      for(a=1;a<256;a++) {
        if(chains[a]) {
          result[c]=(char)a;
          c++;
        }
      }
    }

  } else {
    ErrMessage(G,"ExecutiveGetChains","Bad selection.");
  }
  return(result);
}

int ExecutiveValidateObjectPtr(PyMOLGlobals *G,CObject *ptr,int object_type)
{
  register CExecutive *I = G->Executive;
  int ok=false;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->obj == ptr) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==object_type) {
          ok=true;
          break;
        }
      }
    }
  }
  return(ok);
}

int ExecutiveRampMapNew(PyMOLGlobals *G,char *name,char *map_name,PyObject *range,
                        PyObject *color,int map_state,char *sele,
                        float beyond,float within,float sigma,int zero)
{
  ObjectGadgetRamp *obj = NULL;
  int ok =true;
  CObject *map_obj;
  float *vert_vla = NULL;
  map_obj = ExecutiveFindObjectByName(G,map_name);
  if(map_obj) {
    if(map_obj->type!=cObjectMap) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExecutiveRampMapNew: Error: object '%s' is not a map.\n",map_name
        ENDFB(G);
      ok=false;
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)
      "ExecutiveRampMapNew: Error: map '%s' not found.\n",map_name
      ENDFB(G);
    ok = false;
  }
  if(sele&&sele[0]) {
    vert_vla = ExecutiveGetVertexVLA(G,sele,map_state);
  }
  ok = ok && (obj=ObjectGadgetRampMapNewAsDefined(G,(ObjectMap*)map_obj,
                                                  range,color,map_state,
                                                  vert_vla,beyond,within,
                                                  sigma,zero));
  if(ok) ExecutiveDelete(G,name); 
  if(ok) ObjectSetName((CObject*)obj,name);
  if(ok) ColorRegisterExt(G,name,(void*)obj,cColorGadgetRamp);
  if(ok) ExecutiveManageObject(G,(CObject*)obj,false,false);
  VLAFreeP(vert_vla);
  return(ok);
}


#ifndef _PYMOL_NOPY
static PyObject *ExecutiveGetExecObject(PyMOLGlobals *G,SpecRec *rec)
{

  PyObject *result = NULL;
  result = PyList_New(6);
  PyList_SetItem(result,0,PyString_FromString(rec->obj->Name));
  PyList_SetItem(result,1,PyInt_FromLong(cExecObject));
  PyList_SetItem(result,2,PyInt_FromLong(rec->visible));
  PyList_SetItem(result,3,PConvIntArrayToPyList(rec->repOn,cRepCnt));
  PyList_SetItem(result,4,PyInt_FromLong(rec->obj->type));
  switch(rec->obj->type) {
  case cObjectGadget:
    PyList_SetItem(result,5,ObjectGadgetAsPyList((ObjectGadget*)rec->obj));    
    break;
  case cObjectMolecule:
    PyList_SetItem(result,5,ObjectMoleculeAsPyList((ObjectMolecule*)rec->obj));
    break;
  case cObjectDist:
    PyList_SetItem(result,5,ObjectDistAsPyList((ObjectDist*)rec->obj));
    break;
  case cObjectMap:
    PyList_SetItem(result,5,ObjectMapAsPyList((ObjectMap*)rec->obj));
    break;
  case cObjectMesh:
    PyList_SetItem(result,5,ObjectMeshAsPyList((ObjectMesh*)rec->obj));
    break;
  case cObjectSlice:
    PyList_SetItem(result,5,ObjectSliceAsPyList((ObjectSlice*)rec->obj));
    break;
  case cObjectSurface:
    PyList_SetItem(result,5,ObjectSurfaceAsPyList((ObjectSurface*)rec->obj));
    break;
  case cObjectCGO:
    PyList_SetItem(result,5,ObjectCGOAsPyList((ObjectCGO*)rec->obj));
    break;
  default: 
    PyList_SetItem(result,5,PConvAutoNone(NULL));
    break;
  }
  return(result);  
}

static int ExecutiveSetNamedEntries(PyMOLGlobals *G,PyObject *names,int version)
{
  register CExecutive *I = G->Executive;  
  int ok=true;
  int skip=false;
  int a=0,l=0;
  PyObject *cur;
  SpecRec *rec = NULL;
  int extra_int;
  int incomplete = false;

  if(ok) ok = (names!=NULL);
  if(ok) ok = PyList_Check(names);
  if(ok) l = PyList_Size(names);
  while(ok&&(a<l)) {
    cur = PyList_GetItem(names,a);
    if(cur!=Py_None) { /* skip over None w/o aborting */
      skip=false;
      rec=NULL;
      ListElemCalloc(G,rec,SpecRec); 
      rec->next=NULL;
      rec->name[0]=0;
      if(ok) ok = PyList_Check(cur);
      if(ok) ok = PConvPyStrToStr(PyList_GetItem(cur,0),rec->name,sizeof(WordType));
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,1),&rec->type);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,2),&rec->visible);
      if(ok) ok = PConvPyListToIntArrayInPlaceAutoZero(PyList_GetItem(cur,3),
                                                      rec->repOn,cRepCnt);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,4),&extra_int); 
      switch(rec->type) {
      case cExecObject:
        switch(extra_int) {
        case cObjectMolecule:
          if(ok) ok = ObjectMoleculeNewFromPyList(G,PyList_GetItem(cur,5),(ObjectMolecule**)&rec->obj);
          break;
        case cObjectDist:
          if(ok) ok = ObjectDistNewFromPyList(G,PyList_GetItem(cur,5),(ObjectDist**)&rec->obj);
          break;
        case cObjectMap:
          if(ok) ok = ObjectMapNewFromPyList(G,PyList_GetItem(cur,5),(ObjectMap**)&rec->obj);
          break;
        case cObjectMesh:
          if(ok) ok = ObjectMeshNewFromPyList(G,PyList_GetItem(cur,5),(ObjectMesh**)&rec->obj);
          break;
        case cObjectSlice:
          if(ok) ok = ObjectSliceNewFromPyList(G,PyList_GetItem(cur,5),(ObjectSlice**)&rec->obj);
          break;	  
        case cObjectSurface:
          if(ok) ok = ObjectSurfaceNewFromPyList(G,PyList_GetItem(cur,5),(ObjectSurface**)&rec->obj);
          break;
        case cObjectCGO:
          if(ok) ok = ObjectCGONewFromPyList(G,PyList_GetItem(cur,5),(ObjectCGO**)&rec->obj,version);
          break;
        case cObjectGadget:
          if(ok) ok = ObjectGadgetNewFromPyList(G,PyList_GetItem(cur,5),(ObjectGadget**)&rec->obj,version);
          break;
        default:
          PRINTFB(G,FB_Executive,FB_Errors)
            " Executive: skipping unrecognized object \"%s\" of type %d.\n",
            rec->name,rec->type
            ENDFB(G);
          skip=true;
          break;
        }
        break;
      case cExecSelection: /* on the first pass, just create an entry in the rec list */
        rec->sele_color=extra_int;
        break;
      }

      if(PyErr_Occurred()) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetNamedEntries-ERROR: after object \"%s\".\n",rec->name
          ENDFB(G);
        PyErr_Print();  
      }

      if(ok&&!skip) {
        switch(rec->type) {
        case cExecObject:        
          if(rec->visible) {
            SceneObjectAdd(G,rec->obj);
          }
          ExecutiveUpdateObjectSelection(G,rec->obj);
          break;
        }
        ListAppend(I->Spec,rec,next,SpecRec);
      } else {
        ListElemFree(rec);
      }
    }
    a++;
    if(!ok) {
      incomplete=true;
      ok=true;
    }
  }
  return(!incomplete);
}

static int ExecutiveSetSelections(PyMOLGlobals *G,PyObject *names)
{
  /* must already have objects loaded at this point... */

  int ok=true;
  int a=0,l=0;
  PyObject *cur;
  SpecRec *rec = NULL;
  int extra;
  int incomplete = false;

  if(ok) ok = (names!=NULL);
  if(ok) ok = PyList_Check(names);
  if(ok) l = PyList_Size(names);
  while(ok&&(a<l)) {
    cur = PyList_GetItem(names,a);
    if(cur!=Py_None) { /* skip over None w/o aborting */
      rec=NULL;
      ListElemCalloc(G,rec,SpecRec); 
      rec->next=NULL;

      if(ok) ok = PyList_Check(cur);
      if(ok) ok = PConvPyStrToStr(PyList_GetItem(cur,0),rec->name,sizeof(WordType));
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,1),&rec->type);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,2),&rec->visible);
      if(ok) ok = PConvPyListToIntArrayInPlaceAutoZero(PyList_GetItem(cur,3),
                                                      rec->repOn,cRepCnt);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,4),&extra);
      switch(rec->type) {
      case cExecSelection:
        ok = SelectorFromPyList(G,rec->name,PyList_GetItem(cur,5));
        break;
      }
      ListElemFree(rec);
    }
    a++;
    if(!ok) {
      incomplete=true;
      ok=true;
    }
  }
  return(!incomplete);
}

static PyObject *ExecutiveGetExecSelePyList(PyMOLGlobals *G,SpecRec *rec)
{
  PyObject *result = NULL;
  int sele;

  sele = SelectorIndexByName(G,rec->name);
  if(sele>=0) {
    result = PyList_New(6);
    PyList_SetItem(result,0,PyString_FromString(rec->name));
    PyList_SetItem(result,1,PyInt_FromLong(cExecSelection));
    PyList_SetItem(result,2,PyInt_FromLong(rec->visible));
    PyList_SetItem(result,3,PConvIntArrayToPyList(rec->repOn,cRepCnt));
    PyList_SetItem(result,4,PyInt_FromLong(-1));
    PyList_SetItem(result,5,SelectorAsPyList(G,sele));
  }
  return(PConvAutoNone(result));
}

static PyObject *ExecutiveGetNamedEntries(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;  
  PyObject *result = NULL;
  int count;
  SpecRec *rec = NULL;

  count = ExecutiveCountNames(G);
  result = PyList_New(count);

  SelectorUpdateTable(G);

  count=0;
  while(ListIterate(I->Spec,rec,next))
	 {
      switch(rec->type) {
      case cExecObject:
        PyList_SetItem(result,count,
                       ExecutiveGetExecObject(G,rec));
        break;
      case cExecSelection:
        PyList_SetItem(result,count,
                       ExecutiveGetExecSelePyList(G,rec));
        break;
      default:
        PyList_SetItem(result,count,PConvAutoNone(NULL));
        break;
      }
      count++;
    }
  return(PConvAutoNone(result));
}
#endif

int ExecutiveGetSession(PyMOLGlobals *G,PyObject *dict)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok=true;
  SceneViewType sv;
  PyObject *tmp;

  tmp = ExecutiveGetNamedEntries(G);
  PyDict_SetItemString(dict,"names",tmp);
  Py_XDECREF(tmp);

  tmp = SelectorSecretsAsPyList(G);
  PyDict_SetItemString(dict,"selector_secrets",tmp);
  Py_XDECREF(tmp);
  
  tmp = SettingGetGlobalsPyList(G);
  PyDict_SetItemString(dict,"settings",tmp);
  Py_XDECREF(tmp);

  tmp = ColorAsPyList(G);
  PyDict_SetItemString(dict,"colors",tmp);
  Py_XDECREF(tmp);

  tmp = ColorExtAsPyList(G);
  PyDict_SetItemString(dict,"color_ext",tmp);
  Py_XDECREF(tmp);

  tmp = PyInt_FromLong(_PyMOL_VERSION_int);
  PyDict_SetItemString(dict,"version",tmp);
  Py_XDECREF(tmp);

  SceneGetView(G,sv);
  tmp = PConvFloatArrayToPyList(sv,cSceneViewSize);
  PyDict_SetItemString(dict,"view",tmp);
  Py_XDECREF(tmp);

  tmp = MovieAsPyList(G);
  PyDict_SetItemString(dict,"movie",tmp);
  Py_XDECREF(tmp);

  tmp = EditorAsPyList(G);
  PyDict_SetItemString(dict,"editor",tmp);
  Py_XDECREF(tmp);

#ifndef _PYMOL_NO_MAIN
  tmp = MainAsPyList();
  PyDict_SetItemString(dict,"main",tmp);
  Py_XDECREF(tmp);
#endif

  if(Feedback(G,FB_Executive,FB_Errors)) {
    if(PyErr_Occurred()) {
      PRINTF
        " ExecutiveGetSession: a Python error occured during creation of the session object:\n"
        ENDF(G);
      PyErr_Print();
    }
  }

  return(ok);
#endif

}

#ifndef _PYMOL_NOPY
static void ExecutiveMigrateSession(PyMOLGlobals *G,int session_version)
{
  if(session_version<96) {
    SettingSetGlobal_f(G,cSetting_ray_transparency_contrast, 1.0F);
  }
  if(session_version<95) {

    { /* adjust fog to reflect current importance of seeing to the Z-slab center w/o fog */
      
      float fog_start = SettingGetGlobal_f(G,cSetting_fog_start);
      float ray_trace_fog_start = SettingGetGlobal_f(G,cSetting_ray_trace_fog_start);
      if((fog_start==0.40F)||(fog_start==0.35F)||(fog_start==0.30F)) {
        SettingSetGlobal_f(G,cSetting_fog_start,0.45F);
      }
      if((ray_trace_fog_start==0.45F)||(ray_trace_fog_start==0.40F)||(ray_trace_fog_start==0.35F)) {
        SettingSetGlobal_f(G,cSetting_ray_trace_fog_start,0.50F);
      }

    }

    { /* adjust GUI width */

      int gui_width = SettingGetGlobal_i(G,cSetting_internal_gui_width);

      if(gui_width==160) {
        SettingSetGlobal_i(G,cSetting_internal_gui_width,220);
      }
    }

    { /* enable antialiasing */

      int antialias = SettingGetGlobal_i(G,cSetting_antialias);

      if(antialias==0) {
        SettingSetGlobal_i(G,cSetting_antialias,1);
      }
      
    }
  }
}
#endif

int ExecutiveSetSession(PyMOLGlobals *G,PyObject *session)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok=true;
  int incomplete = false;
  PyObject *tmp;
  SceneViewType sv;
  int version=-1;
  int migrate_sessions = SettingGetGlobal_b(G,cSetting_session_migration);
  char active[ObjNameMax] = "";
  int  have_active = false;

  ExecutiveDelete(G,"all");
  ColorReset(G);
  if(ok) ok = PyDict_Check(session);

  if(ok) {
    tmp = PyDict_GetItemString(session,"version");
    if(tmp) {
      ok = PConvPyIntToInt(tmp,&version);
      if(ok) {
        if(version>_PyMOL_VERSION_int) {
          PRINTFB(G,FB_Executive,FB_Errors)
            "Warning: This session was created with a newer version of PyMOL (%1.2f).\n",
            version/100.
            ENDFB(G);
          if(SettingGet(G,cSetting_session_version_check)) {
            PRINTFB(G,FB_Executive,FB_Errors)
              "Error: Please update first -- see http://www.pymol.org\n"
              ENDFB(G);
            ok=false;
          } else {
            PRINTFB(G,FB_Executive,FB_Errors)
              "Warning: Some content may not load completely.\n"
              ENDFB(G);
          }
        } else {
          PRINTFB(G,FB_Executive,FB_Details)          
            " Executive: Loading version %1.2f session...\n",
            version/100.0
            ENDFB(G);
        }
      }
    }
  }

  if(ok) {
    tmp = PyDict_GetItemString(session,"colors");
    if(tmp) {
      ok = ColorFromPyList(G,tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after colors.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"color_ext");
    if(tmp) {
      ok = ColorExtFromPyList(G,tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after color_ext.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"settings");
    if(tmp) {
      ok = SettingSetGlobalsFromPyList(G,tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after settings.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"names");
    if(tmp) {
      if(ok) ok=ExecutiveSetNamedEntries(G,tmp,version);
      if(ok) ok=ExecutiveSetSelections(G,tmp);
      if(ok) have_active = ExecutiveGetActiveSeleName(G,active,false);
    }
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after names.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"selector_secrets");
    if(tmp) {
      if(ok) ok=SelectorSecretsFromPyList(G,tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after selector secrets.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }  
  if(ok) {
    tmp = PyDict_GetItemString(session,"view");
    if(tmp) {
      ok = PConvPyListToFloatArrayInPlace(tmp,sv,cSceneViewSize);
    }
    if(ok) SceneSetView(G,sv,true,0);
    
    
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after view.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  if(ok) {
    int warning;
    tmp = PyDict_GetItemString(session,"movie");
    if(tmp) {
      ok = MovieFromPyList(G,tmp,&warning);
    }
    
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after movie.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  
  if(ok) {
    tmp = PyDict_GetItemString(session,"editor");
    if(tmp) {
      ok = EditorFromPyList(G,tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after editor.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  if(ok) { /* update mouse in GUI */
    PParse("cmd.mouse(quiet=1)");
    PParse("viewport"); /* refresh window/internal_gui status */
  }
#ifndef _PYMOL_NO_MAIN
  if(ok) {
    tmp = PyDict_GetItemString(session,"main");
    if(tmp) {
      ok = MainFromPyList(tmp);
    }
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after main.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
#endif
  if(ok&&migrate_sessions) { /* migrate sessions */
    tmp = PyDict_GetItemString(session,"version");
    if(tmp) {
      ok = PConvPyIntToInt(tmp,&version);
      if(ok) {
        ExecutiveMigrateSession(G,version);
      }
    }
  }
  if(ok) {
    if(have_active)
      ExecutiveSetObjVisib(G,active,true);      
  }
  if(incomplete) {
    PRINTFB(G,FB_Executive,FB_Warnings)
      "ExectiveSetSession-Warning: restore may be incomplete.\n"
      ENDFB(G);
  }
  return(ok);
#endif
}

#define ExecScrollBarMargin 1
#define ExecScrollBarWidth 13

void ExecutiveObjMolSeleOp(PyMOLGlobals *G,int sele,ObjectMoleculeOpRec *op);

static CObject **ExecutiveSeleToObjectVLA(PyMOLGlobals *G,char *s1)
{
  /* return VLA containing list of atoms references by selection */

  CObject **result = NULL;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  CObject *obj=NULL;
  int n = 0;
  ObjectMoleculeOpRec op2;
  int sele;

  result = VLAlloc(CObject*,50);
  if(WordMatch(G,s1,cKeywordAll,true)) {
    /* all objects */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject)
          {
            VLACheck(result,CObject*,n);
            result[n]=rec->obj;
            n++;
          }
      }
  } else {
    sele = SelectorIndexByName(G,s1);
    if(sele>0) {
      ObjectMoleculeOpRecInit(&op2);
      op2.code=OMOP_GetObjects;
      op2.obj1VLA=(ObjectMolecule**)result;
      op2.i1=0;
      ExecutiveObjMolSeleOp(G,sele,&op2);
      n = op2.i1;
      result = (CObject**)op2.obj1VLA;
    } else {
      obj = ExecutiveFindObjectByName(G,s1);
      if(obj) {
        VLACheck(result,CObject*,n);
        result[n]=obj;
        n++;
      }
    }
  }
  VLASize(result,CObject*,n);
  return(result);
}

int ExecutiveGetCrystal(PyMOLGlobals *G,char *sele,float *a,float *b,float *c,
                         float *alpha,float *beta,float *gamma,
                        char *sgroup,int *defined)
{
  int ok=true;

  ObjectMolecule *objMol;
  int sele0;
  sele0 = SelectorIndexByName(G,sele);
  *defined = false;
  if(sele0<0) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Error: invalid selection.\n"
      ENDFB(G);
    ok=false;
  } else {
    objMol = SelectorGetSingleObjectMolecule(G,sele0);
    if(!objMol) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "Error: selection must refer to exactly one object.\n"
        ENDFB(G);
    ok=false;
    } else {
      if(objMol->Symmetry&&objMol->Symmetry->Crystal) {
        *a = objMol->Symmetry->Crystal->Dim[0];
        *b = objMol->Symmetry->Crystal->Dim[1];
        *c = objMol->Symmetry->Crystal->Dim[2];
        *alpha = objMol->Symmetry->Crystal->Angle[0];
        *beta = objMol->Symmetry->Crystal->Angle[1];
        *gamma = objMol->Symmetry->Crystal->Angle[2];
        UtilNCopy(sgroup, objMol->Symmetry->SpaceGroup,sizeof(WordType));
        *defined = true;
      }
    }
  }
  return(ok);
}

int ExecutiveSetCrystal(PyMOLGlobals *G,char *sele,float a,float b,float c,
                         float alpha,float beta,float gamma,char *sgroup)
{
  CObject **objVLA = NULL;
  CObject *obj;
  ObjectMolecule *objMol;
  /*  ObjectMap *objMap;*/
  int ok=true;
  CSymmetry *symmetry = NULL;
  CCrystal *crystal = NULL;
  int n_obj;
  int i;

  objVLA = ExecutiveSeleToObjectVLA(G,sele);
  n_obj = VLAGetSize(objVLA);
  if(n_obj) {
    for(i=0;i<n_obj;i++) {
      obj = objVLA[i];
      switch(obj->type) {
      case cObjectMolecule:
        if(!symmetry) {
          symmetry=SymmetryNew(G);          
          symmetry->Crystal->Dim[0]=a;
          symmetry->Crystal->Dim[1]=b;
          symmetry->Crystal->Dim[2]=c;
          symmetry->Crystal->Angle[0]=alpha;
          symmetry->Crystal->Angle[1]=beta;
          symmetry->Crystal->Angle[2]=gamma;
          UtilNCopy(symmetry->SpaceGroup,sgroup,sizeof(WordType));
          SymmetryAttemptGeneration(symmetry,false,false);
        }
        objMol = (ObjectMolecule*)obj;
        if(symmetry) {
          if(objMol->Symmetry)
            SymmetryFree(objMol->Symmetry);
          objMol->Symmetry = SymmetryCopy(symmetry);
        }
        break;
        /* not currently supported 
      case cObjectMap:
        
        if(!crystal) {
          crystal = CrystalNew(G);
          crystal->Dim[0]=a;
          crystal->Dim[1]=b;
          crystal->Dim[2]=c;
          crystal->Angle[0]=alpha;
          crystal->Angle[1]=beta;
          crystal->Angle[2]=gamma;
          CrystalUpdate(crystal);
        }
        if(crystal) {
          objMap = (ObjectMap*)obj;
          if(objMap->Crystal) {
            CrystalFree(objMap->Crystal);
          }
          objMap->Crystal=CrystalCopy(crystal);
        }
        break;
        */

      }
    }
  } else {
    ok=false;
    PRINTFB(G,FB_Executive,FB_Errors)
      " ExecutiveSetCrystal: no object selected\n"
      ENDFB(G);
  }
  if(crystal)
    CrystalFree(crystal);
  if(symmetry)
    SymmetryFree(symmetry);
  VLAFreeP(objVLA);
  return(ok);
}

int ExecutiveSmooth(PyMOLGlobals *G,char *name,int cycles,int window,int first, int last, int ends, int quiet)
{
  int sele = -1;
  ObjectMoleculeOpRec op;
  int state;
  int n_state;
  float *coord0=NULL,*coord1=NULL;
  int *flag0=NULL,*flag1=NULL;
  int a,b,c,d,st,cnt;
  float i_cnt;
  int n_atom;
  int ok=true;
  int backward;
  int forward;
  int range,offset;
  int end_skip=0;
  float *v0,*v1;
  float sum[3];
  /*  WordType all = "_all";*/

  PRINTFD(G,FB_Executive)
    " ExecutiveSmooth: entered %s,%d,%d,%d,%d,%d\n",name,cycles,first,last,window,ends
    ENDFD;

  sele=SelectorIndexByName(G,name);



  if(sele>=0) {
    if(last<0) 
      last = ExecutiveCountStates(G,name)-1;
    if(first<0)
      first = 0;
    if(last<first) {
      state=last;
      last=first;
      first=state;
    }
    n_state=last-first+1;

    backward=window/2;
    forward=window/2;

    if((forward-backward)==(window+1))
      forward--; /* even sizes window */
    
    switch(ends) {
    case 0:
      end_skip = 1;
      break;
    case 1:
      end_skip = 0;
      break;
    case 2:
      end_skip = backward;
      break;
    default:
      end_skip = 0;
      break;
    }

    if(ends) {
      range = (last-first)+1;
      offset = 0;
    } else {
      range = (last-end_skip)-(first+end_skip)+1;
      offset = end_skip;
    }
    
    PRINTFD(G,FB_Executive)
      " ExecutiveSmooth: first %d last %d n_state %d backward %d forward %d range %d\n",
      first,last,n_state,backward,forward,range
      ENDFD;

    if(n_state>=window) {
      
      /* determine storage req */
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_CountAtoms;
      op.i1=0;
      ExecutiveObjMolSeleOp(G,sele,&op);
      n_atom=op.i1;
      if(n_atom) {
        /* allocate storage */
        coord0 = Alloc(float,3*n_atom*n_state);
        coord1 = Alloc(float,3*n_atom*n_state);
        flag0 = Alloc(int,n_atom*n_state);
        flag1 = Alloc(int,n_atom*n_state);
        
        /* clear the arrays */
        
        UtilZeroMem(coord0,sizeof(float)*3*n_atom*n_state);
        UtilZeroMem(flag0,sizeof(int)*n_atom*n_state);
        
        /* get the data */
        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Smooth: copying coordinates to temporary arrays..\n"
            ENDFB(G);
        }
        op.code = OMOP_CSetIdxGetAndFlag;
        op.i1 = n_atom; 
        op.i2 = 0;
        op.cs1 = first;
        op.cs2 = last;
        op.vv1 = coord0;
        op.ii1 = flag0;
        op.nvv1 = 0;          
        ExecutiveObjMolSeleOp(G,sele,&op);    
        
        PRINTFD(G,FB_Executive)  
          " ExecutiveSmooth: got %d %d\n",op.i2,op.nvv1
          ENDFD;
        
        UtilZeroMem(coord1,sizeof(float)*3*n_atom*n_state);
        UtilZeroMem(flag1,sizeof(int)*n_atom*n_state);
        
        for(a=0;a<cycles;a++) {                
          if(!quiet) {
            PRINTFB(G,FB_Executive,FB_Actions)
              " Smooth: smoothing (pass %d)...\n",a+1
              ENDFB(G);
          }
          for(b=0;b<range;b++) {
            for(c=0;c<n_atom;c++) {
              zero3f(sum);
              cnt = 0;
              for(d=-backward;d<=forward;d++) {
                st = b + offset + d;
                if(st<0) {
                  st=0;
                } else if(st>=n_state) {
                  st=n_state-1;
                }
                /*if(c==0) printf("averaging from slot %d\n",st);*/
                cnt+=flag0[(n_atom*st)+c];
                v0 = coord0 + 3*(n_atom*st+c);
                add3f(sum,v0,sum);
              }
              if(cnt) {
                st = b + offset;
                if((st>=end_skip)&&(st<(n_state-end_skip))) {
                  /* if(c==0) printf("dumping into slot %d\n",st);*/
                  flag1[(n_atom*st)+c] = 1;
                  i_cnt = 1.0F/cnt;
                  v1 = coord1 + 3*((n_atom*st)+c);
                  scale3f(sum,i_cnt,v1);
                }
              }
            }
          }
          for(b=0;b<range;b++) {
            for(c=0;c<n_atom;c++) {
              st = b + offset;
              if(flag1[(n_atom*st)+c]) {
                v0 = coord0 + 3*((n_atom*st)+c);
                v1 = coord1 + 3*((n_atom*st)+c);
                copy3f(v1,v0);
              }
            }
          }
        }

        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Smooth: updating coordinates...\n"
            ENDFB(G);
        }
        
        /* set the new coordinates */
        
        op.code = OMOP_CSetIdxSetFlagged;
        op.i1 = n_atom; 
        op.i2 = 0;
        if(ends) {
          op.cs1 = first;
          op.cs2 = last;
          op.vv1 = coord1;
          op.ii1 = flag1;
        } else {
          op.cs1 = first+end_skip;
          op.cs2 = last-end_skip;
          op.vv1 = coord1+(end_skip*3*n_atom);
          op.ii1 = flag1+(end_skip*n_atom);
        }
        op.nvv1 = 0;
        
        ExecutiveObjMolSeleOp(G,sele,&op);      
          PRINTFD(G,FB_Executive)  
            " ExecutiveSmooth: put %d %d\n",op.i2,op.nvv1
            ENDFD;
        
        FreeP(coord0);
        FreeP(coord1);
        FreeP(flag0);
        FreeP(flag1);
      }
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)  
      " ExecutiveSmooth: selection not found\n"
      ENDFB(G);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveDebug(PyMOLGlobals *G,char *name)
{
  ObjectMolecule *obj;
  ObjectMoleculeBPRec bp;
  int a;

  obj=(ObjectMolecule*)ExecutiveFindObjectByName(G,name);
  if(obj) {
    ObjectMoleculeInitBondPath(obj,&bp);
    ObjectMoleculeGetBondPaths(obj,0,10,&bp);
    for(a=0;a<bp.n_atom;a++) {
      printf("%d %d %d\n",a,bp.list[a],bp.dist[bp.list[a]]);
    }
    
    ObjectMoleculePurgeBondPath(obj,&bp);
  }
  return(1);
}
/*========================================================================*/
int ***ExecutiveGetBondPrint(PyMOLGlobals *G,char *name,int max_bond,int max_type,int *dim)
{
  int ***result = NULL;
  CObject *obj;
  ObjectMolecule *objMol;

  obj=ExecutiveFindObjectByName(G,name);
  if(obj->type==cObjectMolecule) {
    objMol = (ObjectMolecule*)obj;
    result = ObjectMoleculeGetBondPrint(objMol,max_bond,max_type,dim);
  }
  return(result);
}
/*========================================================================*/
int ExecutiveMapNew(PyMOLGlobals *G,char *name,int type,float *grid,
                    char *sele,float buffer,
                    float *minCorner,
                    float *maxCorner,int state,int have_corners,int quiet)
{
  CObject *origObj=NULL;
  ObjectMap *objMap;
  ObjectMapState *ms = NULL;
  int a;
  float v[3];
  ObjectMapDesc _md,*md;
  int ok = true;
  int sele0 = SelectorIndexByName(G,sele);
  int isNew=true;
  int n_state;
  int valid_extent=false;
  int st;
  int st_once_flag=true;
  int n_st;

  md=&_md;

  if(state==-2) state=SceneGetState(G);

  /* remove object if it already exists */

  origObj=ExecutiveFindObjectByName(G,name);

  if(origObj) {
    if(origObj->type!=cObjectMap) {
      ExecutiveDelete(G,origObj->Name);
    } else {
      isNew=false;
    }
  }

  n_st = ExecutiveCountStates(G,NULL);

  for(st=0;st<n_st;st++) {
    if(state==-1) st_once_flag=false; /* each state, separate map, separate extent */
    if(!st_once_flag) state=st;
    
    if(strlen(sele)&&(!have_corners)) {
      valid_extent = ExecutiveGetExtent(G,sele,md->MinCorner,
                                        md->MaxCorner,true,state,false); /* TODO restrict to state */
    } else {
      valid_extent = 1;
      copy3f(minCorner,md->MinCorner);
      copy3f(maxCorner,md->MaxCorner);
    }
    copy3f(grid,md->Grid);

    subtract3f(md->MaxCorner,md->MinCorner,v);
    for(a=0;a<3;a++) { if(v[a]<0.0) swap1f(md->MaxCorner+a,md->MinCorner+a); };
    subtract3f(md->MaxCorner,md->MinCorner,v);

    if(buffer!=0.0F) {
      for(a=0;a<3;a++) {
        md->MinCorner[a]-=buffer;
        md->MaxCorner[a]+=buffer;
      }
    }
    md->mode = cObjectMap_OrthoMinMaxGrid;
    md->init_mode=-1; /* no initialization */

    /* validate grid */
    for(a=0;a<3;a++) 
      if(md->Grid[a]<=R_SMALL8) md->Grid[a]=R_SMALL8;

    if(ok) {
      if(isNew)
        objMap = ObjectMapNew(G);
      else
        objMap = (ObjectMap*)origObj;
      if(objMap) {
        int once_flag=true;
        n_state = SelectorCountStates(G,sele0);
        if(valid_extent)
          for(a=0;a<n_state;a++) {
            if(state==-3) once_flag=false; /* -2 = each state, separate map, shared extent */
            if(state==-4) state=-1; /* all states, one map */
            if(!once_flag) state=a;
            ms = ObjectMapNewStateFromDesc(G,objMap,md,state);
            if(!ms)
              ok=false;
          
            if(ok&&ms) {
            
              switch(type) {
              case 0: /* vdw */
                SelectorMapMaskVDW(G,sele0,ms,0.0F,state);
                break;
              case 1: /* coulomb */
                SelectorMapCoulomb(G,sele0,ms,0.0F,state,false,
                                   false,1.0F);
                break;
              case 2: /* gaussian */
                SelectorMapGaussian(G,sele0,ms,0.0F,state);
                break;
              case 3: /* coulomb_neutral */
                SelectorMapCoulomb(G,sele0,ms,0.0F,state,true, false,1.0F);
                break;
              case 4: /* coulomb_local */
                SelectorMapCoulomb(G,sele0,ms,
                                   SettingGetGlobal_f(G,cSetting_coulomb_cutoff),state,false,
                                   true, 2.0F);
                break;
              }
              if(!ms->Active)
                ObjectMapStatePurge(G,ms);
            }
            if(once_flag) break;
          }

        ObjectSetName((CObject*)objMap,name);
        ObjectMapUpdateExtents(objMap);
        if(isNew)
          ExecutiveManageObject(G,(CObject*)objMap,true,false);
        isNew=false;
        origObj = (CObject*)objMap;
      }
      SceneDirty(G);
    }
    if(st_once_flag)
      break;
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSculptIterateAll(PyMOLGlobals *G)
{
  int active = false;

  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;

  int state = SceneGetState(G);
  int cycles = (int)SettingGet(G,cSetting_sculpting_cycles);

  if(SettingGet(G,cSetting_sculpting)) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptIterate(objMol,state,cycles);
          active = true;
        }
      }
    }
  }
  return(active);
}
/*========================================================================*/
float ExecutiveSculptIterate(PyMOLGlobals *G,char *name,int state,int n_cycle)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  register CExecutive *I = G->Executive;
  int ok=true;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  float total_strain = 0.0F;

  if(state<0) state=SceneGetState(G);

  if(WordMatch(G,name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          total_strain+=ObjectMoleculeSculptIterate(objMol,state,n_cycle);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    total_strain=ObjectMoleculeSculptIterate((ObjectMolecule*)obj,state,n_cycle);
  }
  return(total_strain);
}
/*========================================================================*/
int ExecutiveSculptActivate(PyMOLGlobals *G,char *name,int state)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  register CExecutive *I = G->Executive;
  int ok=true;
  if(state<0) state=SceneGetState(G);

  if(WordMatch(G,name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptImprint(objMol,state);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectMoleculeSculptImprint((ObjectMolecule*)obj,state);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSculptDeactivate(PyMOLGlobals *G,char *name)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  register CExecutive *I = G->Executive;

  int ok=true;

  if(WordMatch(G,name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptClear(objMol);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectMoleculeSculptClear((ObjectMolecule*)obj);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSetGeometry(PyMOLGlobals *G,char *s1,int geom,int valence)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int ok=false;

  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_SetGeometry;
    op1.i1 = geom;
    op1.i2 = valence;
    op1.i3 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    if(op1.i3) ok=true;
  } else {
    ErrMessage(G,"SetGeometry","Invalid selection.");
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveMultiSave(PyMOLGlobals *G,char *fname,char *name,int state,int append)
{
  int result=false;
  SpecRec *tRec;
  ObjectMolecule *objMol;
  
  PRINTFD(G,FB_Executive)
    " ExecutiveMultiSave-Debug: entered %s %s.\n",fname,name
    ENDFD;
  tRec = ExecutiveFindSpec(G,name);
  if(tRec) {
    if(tRec->type==cExecObject)
      if(tRec->obj->type==cObjectMolecule) {
        objMol =(ObjectMolecule*)tRec->obj;
        result = ObjectMoleculeMultiSave(objMol,fname,state,append);
      }
  }
  return(result);
  
}
int ExecutiveMapSetBorder(PyMOLGlobals *G,char *name,float level)
{
  int result=false;
  SpecRec *tRec;
  ObjectMap *mobj;

  tRec = ExecutiveFindSpec(G,name);
  if(tRec) {
    if(tRec->type==cExecObject)
      if(tRec->obj->type==cObjectMap) {
        mobj =(ObjectMap*)tRec->obj;
        ObjectMapSetBorder(mobj,level);
        result=true;
      }
  }
  return(result);
}

int ExecutiveMapDouble(PyMOLGlobals *G,char *name,int state)
{
  int result=false;
  SpecRec *tRec;
  ObjectMap *mobj;

  tRec = ExecutiveFindSpec(G,name);
  if(tRec) {
    if(tRec->type==cExecObject)
      if(tRec->obj->type==cObjectMap) {
        mobj =(ObjectMap*)tRec->obj;
        result = ObjectMapDouble(mobj,state);
      }
  }
  return(result);
}

void ExecutiveSelectRect(PyMOLGlobals *G,BlockRect *rect,int mode)
{
  Multipick smp;
  OrthoLineType buffer,buf2;
  char selName[ObjNameMax] = cLeftButSele;
  char prefix[3]="";
  int log_box = 0;
  int logging;
  char empty_string[1] = "";
  char *sel_mode_kw = empty_string;

  logging = (int)SettingGet(G,cSetting_logging);
  if(logging)
    log_box= (int)SettingGet(G,cSetting_log_box_selections);
  /*  if(logging==cPLog_pml)
      strcpy(prefix,"_ ");*/
  smp.picked=VLAlloc(Pickable,1000);
  smp.x=rect->left;
  smp.y=rect->bottom;
  smp.w=rect->right-rect->left;
  smp.h=rect->top-rect->bottom;
  SceneMultipick(G,&smp);
  if(smp.picked[0].index) {
    SelectorCreate(G,cTempRectSele,NULL,NULL,1,&smp);
    if(log_box) SelectorLogSele(G,cTempRectSele);
    switch(mode) {
    case cButModeRect:
      if(mode==cButModeRect) {
        SelectorCreate(G,cLeftButSele,cTempRectSele,NULL,1,NULL);
        if(log_box) {
          sprintf(buf2,"%scmd.select(\"%s\",\"%s\",quiet=1)\n",prefix,cLeftButSele,cTempRectSele);
          PLog(buf2,cPLog_no_flush);
        }
      } 
      break;
    case cButModeSeleAdd:
    case cButModeSeleSub:
        ExecutiveGetActiveSeleName(G,selName,true);
        sel_mode_kw = SceneGetSeleModeKeyword(G);        
        /* intentional omission of break! */
    case cButModeRectAdd:
    case cButModeRectSub:
      if(SelectorIndexByName(G,selName)>=0) {
        if((mode==cButModeRectAdd)||(mode==cButModeSeleAdd)) {
          sprintf(buffer,"(?%s or %s(%s))",selName,sel_mode_kw,cTempRectSele);
          SelectorCreate(G,selName,buffer,NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"(%s)\")\n",prefix,selName,buffer);
            PLog(buf2,cPLog_no_flush);
          }
        } else {
          sprintf(buffer,"(%s(?%s) and not %s(%s))",sel_mode_kw,selName,sel_mode_kw,cTempRectSele);
          SelectorCreate(G,selName,buffer,NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,selName,buffer);
            PLog(buf2,cPLog_no_flush);
          }
        }
      } else {
        if((mode==cButModeRectAdd)||(mode=cButModeSeleAdd)) {
          sprintf(buffer,"%s(?%s)",sel_mode_kw,cTempRectSele);
          SelectorCreate(G,selName,buffer,NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,selName,buffer);
            PLog(buf2,cPLog_no_flush);
          }
        } else {
          SelectorCreate(G,selName,"(none)",NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"(none)\")\n",prefix,selName);
            PLog(buf2,cPLog_no_flush);
          }
        }
      }
      if(SettingGet(G,cSetting_auto_show_selections)) {
        ExecutiveSetObjVisib(G,selName,true);
      }
      break;
    }
    if(log_box) {
      sprintf(buf2,"%scmd.delete(\"%s\")\n",prefix,cTempRectSele);
      PLog(buf2,cPLog_no_flush);
      PLogFlush();
    }
    ExecutiveDelete(G,cTempRectSele);
    WizardDoSelect(G,selName);
  }
  VLAFreeP(smp.picked);
}

int ExecutiveTranslateAtom(PyMOLGlobals *G,char *sele,float *v,int state,int mode,int log)
{
  int ok=true;
  ObjectMolecule *obj0;
  int sele0 = SelectorIndexByName(G,sele);
  int i0;
  if(sele0<0) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Error: bad selection %s.\n",sele
      ENDFB(G);
    ok=false;
  } else{ 
    obj0 = SelectorGetSingleObjectMolecule(G,sele0);
    if(!obj0) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "Error: selection isn't a single atom.\n"
        ENDFB(G);
      ok=false;
    } else {
      i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
      if(i0<0) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "Error: selection isn't a single atom.\n"
          ENDFB(G);
        ok=false;
      } else {
        ObjectMoleculeMoveAtom(obj0,state,i0,v,mode,log);
      }
    }
  }
  return(ok);
}

int ExecutiveCombineObjectTTT(PyMOLGlobals *G,char *name,float *ttt)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  int ok=true;

  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectCombineTTT(obj,ttt);
    SceneDirty(G);
  }
  return(ok);
}

int ExecutiveTransformSelection(PyMOLGlobals *G,int state,char *s1,int log,float *ttt)
{
  int sele=-1;
  ObjectMolecule *obj = NULL;
  ObjectMolecule **vla = NULL;
  int nObj;
  int ok=true;
  int a;

  sele = SelectorIndexByName(G,s1);
  if(sele<0)
    ok=false;
  if(ok) {
    vla=SelectorGetObjectMoleculeVLA(G,sele);
    if(!vla) ok=false;
  }
  if(ok) {
    nObj = VLAGetSize(vla);
    for(a=0;a<nObj;a++) {
      obj=vla[a];
      ObjectMoleculeTransformSelection(obj,state,sele,ttt,log,s1);
    }
  }
  SceneDirty(G);
  VLAFreeP(vla);
  return(ok);
}

int ExecutiveTransformObjectSelection(PyMOLGlobals *G,char *name,int state,char *s1,int log,float *ttt)
{
  int sele=-1;
  ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(G,name);
  int ok=true;

  if(s1[0]) {
    sele = SelectorIndexByName(G,s1);
    if(sele<0)
      ok=false;
  }
  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
  } else if(!ok) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: selection object %s not found.\n",s1
      ENDFB(G);
  } else {
    ObjectMoleculeTransformSelection(obj,state,sele,ttt,log,s1);
  }
  SceneDirty(G);
  return(ok);
}

int ExecutiveValidName(PyMOLGlobals *G,char *name)
{
  int result=true;

  if(!ExecutiveFindSpec(G,name)) {
    if(!WordMatch(G,name,cKeywordAll,true))
      if(!WordMatch(G,name,cKeywordSame,true))
        if(!WordMatch(G,name,cKeywordCenter,true))
          if(!WordMatch(G,name,cKeywordOrigin,true))
            result=false;
  }
  return result;
}

int ExecutivePhiPsi(PyMOLGlobals *G,char *s1,ObjectMolecule ***objVLA,int **iVLA,
                    float **phiVLA,float **psiVLA,int state) 
{
  int sele1=SelectorIndexByName(G,s1);
  int result = false;
  ObjectMoleculeOpRec op1;
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.i1 = 0;
    op1.i2 = state;
    op1.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op1.i1VLA=VLAlloc(int,1000);
    op1.f1VLA=VLAlloc(float,1000);
    op1.f2VLA=VLAlloc(float,1000);
    op1.code=OMOP_PhiPsi;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    result = op1.i1;
    VLASize(op1.i1VLA,int,op1.i1);
    VLASize(op1.obj1VLA,ObjectMolecule*,op1.i1);
    VLASize(op1.f1VLA,float,op1.i1);
    VLASize(op1.f2VLA,float,op1.i1);
    *iVLA=op1.i1VLA;
    *objVLA=op1.obj1VLA;
    *phiVLA=op1.f1VLA;
    *psiVLA=op1.f2VLA;
  } else {
    *objVLA=NULL;
    *iVLA=NULL;
    *phiVLA=NULL;
    *psiVLA=NULL;
  }
  return(result);
}


float ExecutiveAlign(PyMOLGlobals *G,char *s1,char *s2,char *mat_file,float gap,float extend,int skip,
                     float cutoff,int cycles,int quiet,char *oname,
                     int state1,int state2)
{
  int sele1=SelectorIndexByName(G,s1);
  int sele2=SelectorIndexByName(G,s2);
  int *vla1=NULL;
  int *vla2=NULL;
  int na,nb;
  int c;
  float result = 0.0;
  int ok=true;
  CMatch *match = NULL;

  if((sele1>=0)&&(sele2>=0)) {
    vla1=SelectorGetResidueVLA(G,sele1);
    vla2=SelectorGetResidueVLA(G,sele2);
    if(vla1&&vla2) {
      na = VLAGetSize(vla1)/3;
      nb = VLAGetSize(vla2)/3;
      if(na&&nb) {
        match = MatchNew(G,na,nb);
        if (ok) ok = MatchResidueToCode(match,vla1,na);
        if (ok) ok = MatchResidueToCode(match,vla2,nb);
        if (ok) ok = MatchMatrixFromFile(match,mat_file,quiet);
        if (ok) ok = MatchPreScore(match,vla1,na,vla2,nb,quiet);
        result = MatchAlign(match,gap,extend,skip,quiet);
        if(match->pair) { 
          c = SelectorCreateAlignments(G,match->pair,
                                       sele1,vla1,sele2,vla2,
                                       "_align1","_align2",false);
          if(c) {
            if(!quiet) {
              PRINTFB(G,FB_Executive,FB_Actions)
                " ExecutiveAlign: %d atoms aligned.\n",c
                ENDFB(G);
            }
            result =ExecutiveRMS(G,"_align1","_align2",2,cutoff,cycles,quiet,oname,
                                 state1,state2,false);
            
          }
        }
        if(match) 
          MatchFree(match);
      }
    }
  }
  VLAFreeP(vla1);
  VLAFreeP(vla2);
  return result;
}

int ExecutivePairIndices(PyMOLGlobals *G,char *s1,char *s2,int state1,int state2,
                         int mode,float cutoff,float h_angle,
                         int **indexVLA, ObjectMolecule ***objVLA)
{
  int result = 0;
  int sele1,sele2;

  sele1 = SelectorIndexByName(G,s1);
  sele2 = SelectorIndexByName(G,s2);
  if((sele1>=0)&&(sele2>=0)) {
    result=SelectorGetPairIndices(G,sele1,state1,sele2,state2,
                                  mode,cutoff,h_angle,indexVLA,objVLA);
  } else {
    ErrMessage(G,"ExecutivePairIndices","One or more bad selections.");
  }
  return(result);
}


int ExecutiveCartoon(PyMOLGlobals *G,int type,char *s1)
{
  int sele1;
  ObjectMoleculeOpRec op1;

  sele1=SelectorIndexByName(G,s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2=0;
  if(sele1>=0) {
    op1.code=OMOP_INVA;
    op1.i1=cRepCartoon; 
    op1.i2=cRepInvRep;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    op1.code = OMOP_Cartoon;
    op1.i1 = type;
    op1.i2 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
  } else {
    ErrMessage(G,"Cartoon","Invalid selection.");
  }
  return(op1.i2);
}
/*========================================================================*/
float *ExecutiveGetVertexVLA(PyMOLGlobals *G,char *s1,int state)
{
  /* returns NULL if none found */

  float *result = NULL;
  ObjectMoleculeOpRec op1;
  int sele1;
  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.nvv1 = 0;
    op1.vv1=VLAlloc(float,1000);
    if(state>=0) {
      op1.cs1 = state;
      op1.code=OMOP_SingleStateVertices;
    } else {
      op1.code=OMOP_VERT;
    }
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    VLASize(op1.vv1,float,op1.nvv1*3);
    result = op1.vv1;
  }
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveGetSettingText(PyMOLGlobals *G,int index,char *object,int state)
{ 
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* Assumes blocked Python interpreter */
  PyObject *result = NULL;
  OrthoLineType buffer = "";
  CObject *obj = NULL;
  CSetting **handle=NULL,*set_ptr1=NULL,*set_ptr2=NULL;
  int ok=true;

  if(object)
    if(object[0]) {
      obj=ExecutiveFindObjectByName(G,object);
      if(!obj) 
        ok=false;
    } 
  if(!ok) {
    PRINTFB(G,FB_Executive,FB_Errors)
      " SettingGet-Error: object \"%s\" not found.\n",object
      ENDFB(G);
    ok=false;
  } else if(obj) {
    handle = obj->fGetSettingHandle(obj,-1);
    if(handle) set_ptr1 = *handle;
    if(state>=0) {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) 
        set_ptr2 = *handle;
      else {
        PRINTFB(G,FB_Executive,FB_Errors)
          " SettingGet-Error: object \"%s\" lacks state %d.\n",object,state+1
          ENDFB(G);
        ok=false;
      }
    }
  }
  if(ok) {
    buffer[0]=0;
  SettingGetTextValue(G,set_ptr2,set_ptr1,index,buffer);
  result=Py_BuildValue("s",buffer);
  } 
  
  return(result);
#endif

}
/*========================================================================*/
PyObject *ExecutiveGetSettingTuple(PyMOLGlobals *G,int index,char *object,int state)
{ /* Assumes blocked Python interpreter */
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;
  CSetting **handle = NULL;
  CObject *obj=NULL;
  int ok = true;
  PRINTFD(G,FB_Executive)
    " ExecutiveGetSettingTuple: object %p state %d\n",object,state
    ENDFD;

  if(object[0]==0) /* global */
    result = SettingGetTuple(G,NULL,NULL,index);
  else {

    if(strlen(object)) {
      obj=ExecutiveFindObjectByName(G,object);
      if(!obj) 
        ok=false;
    } else ok=false;
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        " Executive: object not found.\n"
        ENDFB(G);
    } else {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) 
        result = SettingGetDefinedTuple(G,*handle,index);      
    }
  }
  if(!ok) {
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return(result);
#endif
}
/*========================================================================*/
void ExecutiveSetLastObjectEdited(PyMOLGlobals *G,CObject *o)
{
  register CExecutive *I = G->Executive;
  I->LastEdited = o;
}
/*========================================================================*/
CObject *ExecutiveGetLastObjectEdited(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  return(I->LastEdited);
}
/*========================================================================*/
int ExecutiveSaveUndo(PyMOLGlobals *G,char *s1,int state)
{
  int sele1;
  ObjectMoleculeOpRec op1;

  if(state<0) state = SceneGetState(G);                
  sele1=SelectorIndexByName(G,s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2=0;
  if(sele1>=0) {
    op1.code = OMOP_SaveUndo;
    op1.i1 = state;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
  }
  return(op1.i2);
}

/*========================================================================*/
int ExecutiveSetTitle(PyMOLGlobals *G,char *name,int state,char *text)
{
  int result=false;
  ObjectMolecule *obj;
  obj =ExecutiveFindObjectMoleculeByName(G,name);
  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
  } else {
    result = ObjectMoleculeSetStateTitle(obj,state,text);
  }
  SceneDirty(G);
  return(result);
}
/*========================================================================*/
char *ExecutiveGetTitle(PyMOLGlobals *G,char *name,int state)
{
  char *result = NULL;
  ObjectMolecule *obj;
  obj =ExecutiveFindObjectMoleculeByName(G,name);
  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
  } else {
    result = ObjectMoleculeGetStateTitle(obj,state);
  }
  SceneDirty(G);
  return(result);
}
/*========================================================================*/
void ExecutiveHideSelections(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection) {
      if(rec->visible) {
        rec->visible=false;
        SceneDirty(G);
        SeqDirty(G);
      }
    }
  }
}
/*========================================================================*/
void ExecutiveRenderSelections(PyMOLGlobals *G,int curState)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  SpecRec *rec1;
  int sele;
  int no_depth;
  float min_width;
  float gl_width;
  int width;
  int max_width = (int)SettingGetGlobal_f(G,cSetting_selection_width_max);
  float width_scale = SettingGetGlobal_f(G,cSetting_selection_width_scale);

  min_width = SettingGetGlobal_f(G,cSetting_selection_width);

  if(width_scale>=0.0F) {
    width = (int)((2*SettingGetGlobal_f(G,cSetting_stick_radius)/SceneGetScreenVertexScale(G,NULL)));
  if(width<min_width)
    width = min_width;
  if(width>max_width)
    width = max_width;
  } else
    width = min_width;

  no_depth = (int)SettingGet(G,cSetting_selection_overlay);

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection) {

      if(rec->visible) {
        sele = SelectorIndexByName(G,rec->name); /* TODO: speed this up */
        if(sele>=0) {

          if(no_depth)
            glDisable(GL_DEPTH_TEST);
          glDisable(GL_FOG);

          if(rec->sele_color<0)
            glColor3f(1.0F,0.2F,0.6F);
          else
            glColor3fv(ColorGet(G,rec->sele_color));

          gl_width=(float)width;
          if(width>5) {
            if(width&0x1) {
              width++;
              gl_width = (float)width;
            }
          }
          glPointSize(gl_width);
          glBegin(GL_POINTS);
          rec1 = NULL;
          while(ListIterate(I->Spec,rec1,next)) {
            if(rec1->type==cExecObject) {
              if(rec1->obj->type==cObjectMolecule) {
                ObjectMoleculeRenderSele((ObjectMolecule*)rec1->obj,curState,sele);
              }
            }
          }
          glEnd();

          if(width>4) {
            if(width>5) 
              glPointSize(4.0F);
            else 
              glPointSize(3.0F);
            glColor3f(1.0F,1.0F,1.0F);
            
            glBegin(GL_POINTS);
            rec1 = NULL;
            while(ListIterate(I->Spec,rec1,next)) {
              if(rec1->type==cExecObject) {
                if(rec1->obj->type==cObjectMolecule) {
                  ObjectMoleculeRenderSele((ObjectMolecule*)rec1->obj,curState,sele);
                }
              }
            }
            glEnd();
          }

          if(width>2) {
            if(width>5)
              glPointSize(2.0F);
            else if(width&0x1)
              glPointSize(1.0F);
            else
              glPointSize(2.0F);              

            glColor3f(0.0F,0.0F,0.0F);
            glBegin(GL_POINTS);
            rec1 = NULL;
            while(ListIterate(I->Spec,rec1,next)) {
              if(rec1->type==cExecObject) {
                if(rec1->obj->type==cObjectMolecule) {
                  ObjectMoleculeRenderSele((ObjectMolecule*)rec1->obj,curState,sele);
                }
              }
            }
            glEnd();
          }


          if(no_depth)
            glEnable(GL_DEPTH_TEST);
          glEnable(GL_FOG);
        }
      }
    }
  }
}
/*========================================================================*/
int ExecutiveGetDistance(PyMOLGlobals *G,char *s0,char *s1,float *value,int state)
{
  Vector3f v0,v1;
  int sele0=-1,sele1=-1;
  int ok=true;

  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"GetDistance","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"GetDistance","Selection 2 invalid.");    
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"GetDistance","Selection 1 doesn't contain a single atom/vertex.");
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"GetDistance","Selection 2 doesn't contain a single atom/vertex.");
  }
  if(ok) {
    (*value)=(float)diff3f(v0,v1);
  }
  return ok;
}
/*========================================================================*/
int ExecutiveGetAngle(PyMOLGlobals *G,char *s0,char *s1,char *s2,float *value,int state)
{
  Vector3f v0,v1,v2;
  int sele0=-1,sele1=-1,sele2=-1;
  int ok=true;
  float d1[3],d2[3];
  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"GetAngle","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"GetAngle","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(G,s2))<0)
    ok = ErrMessage(G,"GetAngle","Selection 3 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"GetAngle","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"GetAngle","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele2,state,v2))
      ok = ErrMessage(G,"GetAngle","Selection 3 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    subtract3f(v0,v1,d1);
    subtract3f(v2,v1,d2);
    (*value)=rad_to_deg(get_angle3f(d1,d2));
  }
  return ok;
}
/*========================================================================*/
int ExecutiveGetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float *value,int state)
{
  Vector3f v0,v1,v2,v3;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int ok=true;
  
  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(G,s2))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 3 invalid.");
  else if((sele3 = SelectorIndexByName(G,s3))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 4 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"GetDihedral","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"GetDihedral","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele2,state,v2))
      ok = ErrMessage(G,"GetDihedral","Selection 3 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele3,state,v3))
      ok = ErrMessage(G,"GetDihedral","Selection 4 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    (*value)=rad_to_deg(get_dihedral3f(v0,v1,v2,v3));
  }
  return ok;
}
/*========================================================================*/
int ExecutiveSetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float value,int state,int quiet)
{
  Vector3f v0,v1,v2,v3;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int ok=true;
  int save_state;
  float current;
  float change;

  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(G,s2))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 3 invalid.");
  else if((sele3 = SelectorIndexByName(G,s3))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 4 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"GetDihedral","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"GetDihedral","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele2,state,v2))
      ok = ErrMessage(G,"GetDihedral","Selection 3 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele3,state,v3))
      ok = ErrMessage(G,"GetDihedral","Selection 4 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    current=rad_to_deg(get_dihedral3f(v0,v1,v2,v3));
    change=value-current;
    save_state = SceneGetState(G);                
    SceneSetFrame(G,-1,state); /* KLUDGE ALERT!
                             * necessary because the editor 
                             * can only work on the current state...this
                             * needs to be changed.*/
    EditorSelect(G,s2,s1,NULL,NULL,false,true,true);
    EditorTorsion(G,change);
    SceneSetFrame(G,-1,save_state);
    if(!quiet) {
      PRINTFB(G,FB_Editor,FB_Actions)
        " SetDihedral: adjusted to %5.3f\n",value
        ENDFB(G);
    }

  }
  return ok;
}
/*========================================================================*/
float ExecutiveGetArea(PyMOLGlobals *G,char *s0,int sta0,int load_b)
{
  ObjectMolecule *obj0;
  RepDot *rep;
  CoordSet *cs;
  float result=-1.0;
  int a,sele0;
  int known_member=-1;
  int is_member;
  int *ati;
  float *area;
  AtomInfoType *ai=NULL;
  ObjectMoleculeOpRec op;
  sele0 = SelectorIndexByName(G,s0);
  if(sele0<0) {
    ErrMessage(G,"Area","Invalid selection.");
  } else {
    obj0 = SelectorGetSingleObjectMolecule(G,sele0);
    if(!(obj0))
      ErrMessage(G,"Area","Selection must be within a single object.");
    else {
      cs = ObjectMoleculeGetCoordSet(obj0,sta0);
      if(!cs)
        ErrMessage(G,"Area","Invalid state.");
      else {
        rep = (RepDot*)RepDotDoNew(cs,cRepDotAreaType);
        if(!rep) 
          ErrMessage(G,"Area","Can't get dot representation.");
        else {

          if(load_b) {
            /* zero out B-values within selection */
            ObjectMoleculeOpRecInit(&op);
            op.code=OMOP_SetB;
            op.f1=0.0;
            op.i1=0;
            ExecutiveObjMolSeleOp(G,sele0,&op);
          }

          result=0.0;
          
          area=rep->A;
          ati=rep->Atom;
          
          is_member = false;

          for(a=0;a<rep->N;a++) {
            
            if(known_member!=(*ati)) {
              known_member=(*ati);
              ai=obj0->AtomInfo+known_member;
              is_member = SelectorIsMember(G,ai->selEntry,sele0);
            } 

            if(is_member) {
              result+=(*area);
              if(load_b)
                ai->b+=(*area);
            }
            area++;
            ati++;
          }
          
          rep->R.fFree((Rep*)rep); /* free the representation */
        }
      }
    }
  }
  return(result);
}

/*========================================================================*/
char *ExecutiveGetNames(PyMOLGlobals *G,int mode,int enabled_only)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  char *result;
  int size = 0;
  int stlen;
  result=VLAlloc(char,1000);

  while(ListIterate(I->Spec,rec,next)) {
    if(
       (rec->type==cExecObject&&((!mode)||(mode==1)||(mode==3)))||
       (rec->type==cExecSelection&&((!mode)||(mode==2)||(mode==3))))
      {
        if((mode!=3)||(rec->name[0]!='_')) {
          if((!enabled_only)||(rec->visible)) {
            stlen = strlen(rec->name);
            VLACheck(result,char,size+stlen+1);
            strcpy(result+size,rec->name);
            size+=stlen+1;
          }
        }
      }
  }
  VLASize(result,char,size);
  return(result);
}
/*========================================================================*/
int ExecutiveGetType(PyMOLGlobals *G,char *name,WordType type)
{
  SpecRec *rec = NULL;
  int ok=true;
  rec = ExecutiveFindSpec(G,name);
  if(!rec) {
    ok=false;
  } else {
    if(rec->type==cExecObject) {
      strcpy(type,"object:");
      if(rec->obj->type==cObjectMolecule)
        strcat(type,"molecule");
      else if(rec->obj->type==cObjectMap)
        strcat(type,"map");
      else if(rec->obj->type==cObjectMesh)
        strcat(type,"mesh");
      else if(rec->obj->type==cObjectSlice)
        strcat(type,"slice");
      else if(rec->obj->type==cObjectSurface)
        strcat(type,"surface");
      else if(rec->obj->type==cObjectDist)
        strcat(type,"distance");
    } else if(rec->type==cExecSelection) {
      strcpy(type,"selection");
    }
  }
  return(ok);
}

/*========================================================================*/
void ExecutiveUpdateCmd(PyMOLGlobals *G,char *s0,char *s1,int sta0,int sta1)
{
  int sele0,sele1;

  sele0 = SelectorIndexByName(G,s0);
  sele1 = SelectorIndexByName(G,s1);
  if(!(sele0&&sele1)) {
    ErrMessage(G,"Update","One or more invalid input selections.");
  } else {
    SelectorUpdateCmd(G,sele0,sele1,sta0,sta1);
  }
}
/*========================================================================*/
void ExecutiveRenameObjectAtoms(PyMOLGlobals *G,char *name,int force) 
{
  register CExecutive *I = G->Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(G,name);
    if(!os)
      ErrMessage(G," Executive","object not found.");
    else if(os->type!=cObjectMolecule) {
      ErrMessage(G," Executive","bad object type.");
      os = NULL;
    }
  }
  
  if(os||(!strlen(name))) { /* sort one or all */
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          if((!os)||(rec->obj==os)) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeRenameAtoms(obj,force);  
          }
    }
    SceneChanged(G);
  }
} 

/*========================================================================*/
int  ExecutiveInvert(PyMOLGlobals *G,int quiet)
{
  int ok=false;
  ok = EditorInvert(G,quiet);
  return(ok);
}
/*========================================================================*/
void ExecutiveFuse(PyMOLGlobals *G,char *s0,char *s1,int mode)
{
  int i0=-1;
  int i1=-1;
  int sele0,sele1,sele2;
  ObjectMolecule *obj0,*obj1;
  ObjectMoleculeOpRec op;
  
#define tmp_fuse_sele "tmp_fuse_sele"

  sele0 = SelectorIndexByName(G,s0);
  if(sele0>=0) {
    sele1 = SelectorIndexByName(G,s1);
    if(sele1>=0) {
      EditorInactivate(G);
      obj0 = SelectorGetSingleObjectMolecule(G,sele0);
      obj1 = SelectorGetSingleObjectMolecule(G,sele1);
      if(obj0)
        i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
      if(obj1)
        i1 = ObjectMoleculeGetAtomIndex(obj1,sele1);
      if(obj0&&obj1&&(i0>=0)&&(i1>=0)&&(obj0!=obj1)) {
        ObjectMoleculeVerifyChemistry(obj0);
        ObjectMoleculeVerifyChemistry(obj1);
        
        SelectorCreate(G,tmp_fuse_sele,NULL,obj0,1,NULL);
        sele2=SelectorIndexByName(G,tmp_fuse_sele);
        if(mode) {
          ObjectMoleculeOpRecInit(&op);
          op.code=OMOP_PrepareFromTemplate;
          op.ai=obj1->AtomInfo+i1;
          op.i1=mode;
          op.i2=0;
          ExecutiveObjMolSeleOp(G,sele2,&op);
        }
        SelectorDelete(G,tmp_fuse_sele);

        if((obj0->AtomInfo[i0].protons==1)&&
           (obj1->AtomInfo[i1].protons==1))
          ObjectMoleculeFuse(obj1,i1,obj0,i0,0);
        else if((obj0->AtomInfo[i0].protons!=1)&&
                (obj1->AtomInfo[i1].protons!=1))
          ObjectMoleculeFuse(obj1,i1,obj0,i0,1);
        else 
          ErrMessage(G,"Fuse","Can't fuse between a hydrogen and a non-hydrogen");
      }
    }
  }
}

/*========================================================================*/
void ExecutiveSpheroid(PyMOLGlobals *G,char *name,int average)  /* EXPERIMENTAL */
{
  register CExecutive *I = G->Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(G,name);
    if(!os)
      ErrMessage(G," Executive","object not found.");
    else if(os->type!=cObjectMolecule) {
      ErrMessage(G," Executive","bad object type.");
      os=NULL;
    }
  }
  
  if(os||(!strlen(name))) { /* sort one or all */
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          if((!os)||(rec->obj==os)) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeCreateSpheroid(obj,average);  
            ObjectMoleculeInvalidate(obj,cRepAll,cRepInvRep);
          }
    }
    SceneChanged(G);
  }
} 
/*========================================================================*/
void ExecutiveRebuildAll(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  PRINTFD(G,FB_Executive)
    " ExecutiveRebuildAll: entered.\n"
    ENDFD;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      switch(rec->obj->type) {
      case cObjectMolecule:
        ObjectMoleculeInvalidate((ObjectMolecule*)rec->obj,cRepAll,cRepInvRep);
        break;
      case cObjectDist:
        ObjectDistInvalidateRep((ObjectDist*)rec->obj,cRepAll);
        break;
      case cObjectSurface:
      case cObjectMesh:
      case cObjectSlice:
        if(rec->obj->fInvalidate) {
          rec->obj->fInvalidate((CObject*)rec->obj,cRepAll,cRepInvAll,-1);
        }
        break;
      }
    }
  }
  SeqChanged(G);
  SceneDirty(G);
}
/*========================================================================*/
void ExecutiveRebuildAllObjectDist(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      if(rec->obj->type==cObjectDist) {
        ObjectDistInvalidateRep((ObjectDist*)rec->obj,cRepAll);
      }
    }
  }
  SceneDirty(G);
}
/*========================================================================*/
void ExecutiveUndo(PyMOLGlobals *G,int dir)
{
  register CExecutive *I = G->Executive;
  CObject *o;
  ObjectMolecule *obj=NULL,*compObj;
  SpecRec *rec = NULL;

  o = ExecutiveGetLastObjectEdited(G);
  PRINTFB(G,FB_Executive,FB_Debugging)
    " ExecutiveUndo: last object %p\n",(void*)o
    ENDFB(G);
  if(o)
    if(o->type==cObjectMolecule)
      obj = (ObjectMolecule*)o;
  /* make sure this is still a real object */
  if(obj) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule) {
          compObj=(ObjectMolecule*)rec->obj;
          if(obj==compObj) {
            ObjectMoleculeUndo(obj,dir);
            break;
          }
        }
    }
  }
  
}
/*========================================================================*/
void ExecutiveSort(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;
  ObjectMoleculeOpRec op;
  int all_obj = false;
  int sele;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(G,name);
    if(!os) {
      if(!WordMatchExact(G,cKeywordAll,name,true))
        ErrMessage(G," Executive","object not found.");
      else
        all_obj=true;
    } else if(os->type!=cObjectMolecule)
      ErrMessage(G," Executive","bad object type.");
  } else {
    all_obj = true;
  }
  
  if(os||all_obj) { /* sort one or all */
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          if((rec->obj==os)||all_obj) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeSort(obj);
            sele=SelectorIndexByName(G,rec->obj->Name);
            if(sele>=0) {
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_INVA;
              op.i1=cRepAll; 
              op.i2=cRepInvRep;
              ExecutiveObjMolSeleOp(G,sele,&op);
            }
          }
    }
    SceneChanged(G);
  }
}
/*========================================================================*/
void ExecutiveRemoveAtoms(PyMOLGlobals *G,char *s1,int quiet)
{
  int sele;
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  ObjectMoleculeOpRec op;
  int flag = false;

  sele=SelectorIndexByName(G,s1);
  if(sele>=0)
	 {
		while(ListIterate(I->Spec,rec,next))
		  {
			 if(rec->type==cExecObject)
				{
				  if(rec->obj->type==cObjectMolecule)
					 {
                  ObjectMoleculeOpRecInit(&op);
                  op.code = OMOP_Remove;
                  op.i1 = 0;
						obj=(ObjectMolecule*)rec->obj;
                  ObjectMoleculeVerifyChemistry(obj); /* remember chemistry for later */
						ObjectMoleculeSeleOp(obj,sele,&op);
                  if(op.i1) {
                    if(!quiet) {
                    PRINTFD(G,FB_Editor)
                      " ExecutiveRemove-Debug: purging %i of %i atoms in %s\n",
                      op.i1,obj->NAtom,obj->Obj.Name
                      ENDFD;
                    }
                    ObjectMoleculePurge(obj);
                    if(!quiet) {
                      PRINTFB(G,FB_Editor,FB_Actions)
                        " Remove: eliminated %d atoms in model \"%s\".\n",
                        op.i1,obj->Obj.Name 
                        ENDFB(G);
                    }
                    flag=true;
                  }
					 }
				}
		  }
	 }
  /*  if(!flag) {
      ErrMessage(G,"Remove","no atoms removed.");
      }*/
}
/*========================================================================*/
void ExecutiveAddHydrogens(PyMOLGlobals *G,char *s1,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_AddHydrogens; /* 4 passes completes the job */
    ExecutiveObjMolSeleOp(G,sele1,&op);    
  }
}
/*========================================================================*/
void ExecutiveFlag(PyMOLGlobals *G,int flag,char *s1,int action,int quiet)
{
  int sele1;
  OrthoLineType buffer;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    switch(action) {
    case 0: op.code = OMOP_Flag; break;
    case 1: op.code = OMOP_FlagSet; break;
    case 2: op.code = OMOP_FlagClear; break;
    default:
      op.code = OMOP_Flag;
      break;
    }
    op.i1 = (((unsigned int)1)<<flag);
    op.i2 = ((unsigned int)0xFFFFFFFF - (((unsigned int)1)<<flag));
    op.i3 = 0;
    op.i4 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op);    
    if(Feedback(G,FB_Executive,FB_Actions)) {
      if(!quiet) {
        switch(action) {
        case 0:
          if(op.i3) {
            PRINTF " Flag: flag %d is set in %d of %d atoms.\n", flag, op.i3, op.i4 ENDF(G);
          } else {
            PRINTF " Flag: flag %d cleared on all atoms.\n", flag ENDF(G);
          }
          break;
        case 1:
          PRINTF " Flag: flag %d set on %d atoms.\n", flag, op.i3 ENDF(G);
          break;
        case 2:
          PRINTF " Flag: flag %d cleared on %d atoms.\n", flag, op.i3 ENDF(G);
          break;
        }
      }
    }
    if((int)SettingGet(G,cSetting_auto_indicate_flags)) {
      sprintf(buffer,"(flag %d)",flag);
      SelectorCreate(G,cIndicateSele,buffer,NULL,true,NULL);
      ExecutiveSetObjVisib(G,cIndicateSele,true);
      SceneDirty(G);
    }
  }

}
/*========================================================================*/
float ExecutiveOverlap(PyMOLGlobals *G,char *s1,int state1,char *s2,int state2,float adjust)
{
  int sele1,sele2;
  float result=0.0;

  if(state1<0) state1=0;
  if(state2<0) state2=0;
                 
  sele1=SelectorIndexByName(G,s1);
  sele2=SelectorIndexByName(G,s2);

  if((sele1>=0)&&(sele2>=0))
    result = SelectorSumVDWOverlap(G,sele1,state1,sele2,state2,adjust);

  return(result);
}
/*========================================================================*/
void ExecutiveProtect(PyMOLGlobals *G,char *s1,int mode,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Protect;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op);    
    if(!quiet) {
      if(Feedback(G,FB_Executive,FB_Actions)) {
        if(op.i2) {
          if(mode) {
            PRINTF " Protect: %d atoms protected from movement.\n",op.i2 ENDF(G);
          } else {
            PRINTF " Protect: %d atoms deprotected.\n", op.i2 ENDF(G);
          }
        }
      }
    }
  }
}
/*========================================================================*/
void ExecutiveMask(PyMOLGlobals *G,char *s1,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Mask;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op);
    if(Feedback(G,FB_Executive,FB_Actions)) {    
      if(op.i2) {
        if(mode) {
          PRINTF " Protect: %d atoms masked (can not be picked).\n",op.i2 ENDF(G);
        } else {
          PRINTF " Protect: %d atoms unmasked.\n", op.i2 ENDF(G);
        }
      }
    }
    op.code = OMOP_INVA; /* need to invalidate all pickable representations */
    op.i1 = cRepAll;
    op.i2 = cRepInvPick;
    ExecutiveObjMolSeleOp(G,sele1,&op);    
  }
}
/*========================================================================*/
int ExecutiveStereo(PyMOLGlobals *G,int flag)
{
  int ok=1;
  int stereo_mode;

  switch(flag) {
  case -1:
    SettingSet(G,cSetting_stereo_shift,
               -SettingGet(G,cSetting_stereo_shift));
    /* shouldn't have to swap angle -- that's implicit
       SettingSet(cSetting_stereo_angle,-SettingGet(G,cSetting_stereo_angle));*/
    break;
  default:
    
    if(G->HaveGUI) {
      stereo_mode = (int)SettingGet(G,cSetting_stereo_mode);
      
      switch(stereo_mode) {
      case 1: /* hardware stereo-in-a-window*/
        if(G->StereoCapable||SceneGetStereo(G)) {
          SceneSetStereo(G,flag);
          PSGIStereo(flag);
        } else {
          ok=false;
        }
        break;
      case 2: /* cross-eye stereo*/
      case 3:
        SceneSetStereo(G,flag);
        break;
      
      }
    }
  }
  return(ok);
}
/*========================================================================*/
void ExecutiveBond(PyMOLGlobals *G,char *s1,char *s2,int order,int add)
{
  int sele1,sele2;
  int cnt;
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  int flag = false;

  sele1=SelectorIndexByName(G,s1);
  sele2=SelectorIndexByName(G,s2);
  
  if((sele1>=0)&&(sele2>=0)) {
	 {
		while(ListIterate(I->Spec,rec,next))
		  {
			 if(rec->type==cExecObject)
				{
				  if(rec->obj->type==cObjectMolecule)
					 {
                  if(add==1) {
                    cnt = ObjectMoleculeAddBond((ObjectMolecule*)rec->obj,sele1,sele2,order);
                    if(cnt) {
                      PRINTFB(G,FB_Editor,FB_Actions)
                        " AddBond: %d bonds added to model \"%s\".\n",cnt,rec->obj->Name 
                        ENDFB(G);
                      flag=true;
                    }
                  } else if(add==2) {
                    cnt = ObjectMoleculeAdjustBonds((ObjectMolecule*)rec->obj,sele1,sele2,1,order);                    
                  }
                  else {
                    cnt = ObjectMoleculeRemoveBonds((ObjectMolecule*)rec->obj,sele1,sele2);
                    if(cnt) {
                      PRINTFB(G,FB_Editor,FB_Actions)
                        " RemoveBond: %d bonds removed from model \"%s\".\n",
                        cnt,rec->obj->Name 
                        ENDFB(G);
                      flag=true;
                    }
                  }
                }
            }
        }
      if(!flag) {
        if(add) 
          ErrMessage(G,"AddBond","no bonds added.");
        else
          ErrMessage(G,"RemoveBond","no bonds removed.");          
      }
    }
  } else if(sele1<0) {
    ErrMessage(G,"ExecutiveBond","The first selection contains no atoms.");
  } else if(sele2<0) {
    ErrMessage(G,"ExecutiveBond","The second selection contains no atoms.");
  }
}
/*========================================================================*/
float ExecutiveDist(PyMOLGlobals *G,char *nam,char *s1,char *s2,int mode,float cutoff,
                    int labels,int quiet)
{
  int sele1,sele2;
  ObjectDist *obj;
  CObject *anyObj = NULL;
  float result;
  sele1=SelectorIndexByName(G,s1);
  if(!WordMatch(G,s2,"same",true))
    sele2=SelectorIndexByName(G,s2);
  else {
    sele2 = sele1;
  }
  
  if((sele1>=0)&&(sele2>=0)) {
    anyObj = ExecutiveFindObjectByName(G,nam);
    if(anyObj)
      if(anyObj->type!=cObjectDist)
        ExecutiveDelete(G,nam);
    obj = ObjectDistNewFromSele(G,(ObjectDist*)anyObj,sele1,sele2,mode,cutoff,labels,&result);
    if(!obj) {
      ErrMessage(G,"ExecutiveDistance","No such distances found.");
    } else {
      ObjectSetName((CObject*)obj,nam);
      ExecutiveManageObject(G,(CObject*)obj,true,quiet);
      ExecutiveSetRepVisib(G,nam,cRepLine,1);
      if(!labels)
        ExecutiveSetRepVisib(G,nam,cRepLabel,0);        
    }
  } else if(sele1<0) {
    ErrMessage(G,"ExecutiveDistance","The first selection contains no atoms.");
  } else if(sele2<0) {
    ErrMessage(G,"ExecutiveDistance","The second selection contains no atoms.");
  }
  return(result);
}
/*========================================================================*/
float ExecutiveDistance(PyMOLGlobals *G,char *s1,char *s2)
{
  int sele1,sele2;
  float dist = -1.0;
  
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  
  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  sele1=SelectorIndexByName(G,s1);
  op1.i1=0;
  op2.i2=0;
  if(sele1>=0) {
    op1.code = OMOP_SUMC;
    op1.v1[0]=0.0;
    op1.v1[1]=0.0;
    op1.v1[2]=0.0;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
  } else {
    ErrMessage(G,"ExecutiveDistance","The first selection contains no atoms.");
  }
  
  sele2=SelectorIndexByName(G,s2);
  op2.i1=0;
  op2.i2=0;
  if(sele2>=0) {
    op2.code = OMOP_SUMC;
    op2.v1[0]=0.0;
    op2.v1[1]=0.0;
    op2.v1[2]=0.0;
    op2.i1=0;
    ExecutiveObjMolSeleOp(G,sele2,&op2);
  } else {
    ErrMessage(G,"ExecutiveDistance","The second selection contains no atoms.");
  }
  
  if(op1.i1&&op2.i1) {
    scale3f(op1.v1,1.0F/op1.i1,op1.v1);
    scale3f(op2.v1,1.0F/op2.i1,op2.v1);
    dist = (float)diff3f(op1.v1,op2.v1);
    PRINTFB(G,FB_Executive,FB_Results)
      " Distance: %8.3f [%i atom(s) to %i atom(s)]\n",
      dist,op1.i1,op2.i1
      ENDFB(G);
  } else {
    ErrMessage(G,"ExecutiveRMS","No atoms selected.");
  }
  return(dist);
}
/*========================================================================*/
char *ExecutiveSeleToPDBStr(PyMOLGlobals *G,char *s1,int state,int conectFlag,int mode)
{
  char *result=NULL;
  ObjectMoleculeOpRec op1;
  int sele1;
  char end_str[] = "END\n";
  int model_count = 1;
  int actual_state = 0;
  int n_state = 1;
  int a;
  char model_record[50];
  int count=0,*counter=NULL;
  PDBInfoRec pdb_info;
  ObjectMolecule *obj = NULL;


  UtilZeroMem((void*)&pdb_info,sizeof(PDBInfoRec));
  ObjectMoleculeOpRecInit(&op1);
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    obj = SelectorGetSingleObjectMolecule(G,sele1);
    if(obj)
      if(obj->DiscreteFlag) {
        counter=&count; /* discrete objects need atom counters between states */
      }
  }
  op1.i2 = 0;
  op1.charVLA=VLAlloc(char,10000);
  if(state==-2) { /* multimodel PDB */
    n_state = ExecutiveCountStates(G,s1);
  }

  if(mode==1) {
    pdb_info.is_pqr_file = true;    
  }

  for(a=0;a<n_state;a++) {
    switch(state) {
    case -2:
      sprintf(model_record,"MODEL     %4d\n",model_count++);
      UtilConcatVLA(&op1.charVLA,&op1.i2,model_record);
      actual_state = a;
      break;
    case -1:
      if(state==-1) actual_state=SceneGetState(G);
      break;
    default:
      actual_state = state;
      break;
    }
    
    if(conectFlag) {
      op1.i2=SelectorGetPDB(G,&op1.charVLA,op1.i2,sele1,
                            actual_state,conectFlag,&pdb_info,counter);
    } else {
      op1.i3 = 0; /* atIndex */
      if(sele1>=0) {
        op1.code = OMOP_PDB1;
        op1.i1 = actual_state;
        ExecutiveObjMolSeleOp(G,sele1,&op1);
      }
    }
    if(!(SettingGetGlobal_i(G,cSetting_pdb_no_end_record)))
      /* terminate with END */
      UtilConcatVLA(&op1.charVLA,&op1.i2,end_str);
    switch(state) {
    case -2:
      UtilConcatVLA(&op1.charVLA,&op1.i2,"ENDMDL\n");
      break;
    }
  }

  /* terminate (just in case) */
  VLACheck(op1.charVLA,char,op1.i2+1);
  op1.charVLA[op1.i2]=0;
  op1.i2++;
  
  result=Alloc(char,op1.i2);
  memcpy(result,op1.charVLA,op1.i2);
  VLAFreeP(op1.charVLA);
  
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveSeleToChemPyModel(PyMOLGlobals *G,char *s1,int state)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result;
  int sele1;
  sele1=SelectorIndexByName(G,s1);
  if(state<0) state=0;
  PBlock(); /*   PBlockAndUnlockAPI();*/
  result=SelectorGetChemPyModel(G,sele1,state);
  if(PyErr_Occurred()) PyErr_Print();
  PUnblock(); /* PLockAPIAndUnblock();*/
  return(result);
#endif
}
/*========================================================================*/
void ExecutiveSeleToObject(PyMOLGlobals *G,char *name,char *s1,int source,int target,int discrete)
{
  int sele1;

  sele1=SelectorIndexByName(G,s1);

  SelectorCreateObjectMolecule(G,sele1,name,target,source,discrete);
}
/*========================================================================*/
void ExecutiveCopy(PyMOLGlobals *G,char *src,char *dst)
{
  CObject *os;
  ObjectMolecule *oSrc,*oDst;
  SpecRec *rec1 = NULL,*rec2=NULL;
  int a;

  os=ExecutiveFindObjectByName(G,src);
  if(!os)
    ErrMessage(G," Executive","object not found.");
  else if(os->type!=cObjectMolecule)
    ErrMessage(G," Executive","bad object type.");
  else 
    {
      oSrc =(ObjectMolecule*)os;
      oDst = ObjectMoleculeCopy(oSrc);
      if(oDst) {
        strcpy(oDst->Obj.Name,dst);
        ExecutiveManageObject(G,(CObject*)oDst,true,false);
        rec1=ExecutiveFindSpec(G,oSrc->Obj.Name);
        rec2=ExecutiveFindSpec(G,oDst->Obj.Name);
        if(rec1&&rec2) {
          for(a=0;a<cRepCnt;a++)
            rec2->repOn[a]=rec1->repOn[a];
        }
        
        PRINTFB(G,FB_Executive,FB_Actions)
          " Executive: object %s created.\n",oDst->Obj.Name 
          ENDFB(G);
      }
    }
  SceneChanged(G);
}

/*========================================================================*/
void ExecutiveOrient(PyMOLGlobals *G,char *sele,Matrix33d mi,
                     int state,int animate)
{
  double egval[3],egvali[3];
  double evect[3][3];
  float m[4][4],mt[4][4];
  float t[3];

  int a,b;

  if(!MatrixEigensolve33d(G,(double*)mi,egval,egvali,(double*)evect)) {

	 normalize3d(evect[0]);
	 normalize3d(evect[1]);
	 normalize3d(evect[2]);

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  m[a][b]=(float)evect[b][a]; /* fill columns */
		}
	 }

    for(a=0;a<3;a++) /* expand to 4x4 */
      {
        m[3][a]=0;
        m[a][3]=0;
      }
    m[3][3]=1.0;

    normalize3f(m[0]); /* cross normalization (probably unnec.)  */
    normalize3f(m[1]);
    normalize3f(m[2]);

    for(a=0;a<3;a++) /* convert to row-major */
      for(b=0;b<3;b++)
        mt[a][b]=m[b][a];

    cross_product3f(mt[0],mt[1],t);     /* insure right-handed matrix */
    if(dot_product3f(t,mt[2])<0.0) {
      mt[2][0] = -mt[2][0];
      mt[2][1] = -mt[2][1];
      mt[2][2] = -mt[2][2];
    }

    for(a=0;a<3;a++) /* convert back to column major */
      for(b=0;b<3;b++)
        m[a][b]=mt[b][a];

    if(animate<0)
      animate=SettingGetGlobal_b(G,cSetting_animation);
    if(animate)
      ScenePrimeAnimation(G);

    SceneSetMatrix(G,m[0]); /* load matrix */

    /* there must  be a more elegant to get the PC on X and the SC
     * on Y then what is shown below, but I couldn't get it to work.
     * I tried swapping the eigen-columns around but either that is 
     * a bogus approach (?) or my code was buggy.  Hence the following...*/

    if((egval[0]<egval[2])&&(egval[2]<egval[1])) { /* X < Z < Y */
      SceneRotate(G,90,1,0,0); /*1<-->2*/
    } else if((egval[1]<egval[0])&&(egval[0]<egval[2])) { /* Y < X < Z */
      SceneRotate(G,90,0,0,1); /*0<-->1*/
    } else if((egval[1]<egval[2])&&(egval[2]<egval[0])) { /* Y < Z < X */
      SceneRotate(G,90,0,1,0); /*1<-->2*/
      SceneRotate(G,90,0,0,1); /*0<-->1*/
    } else if((egval[2]<egval[1])&&(egval[1]<egval[0])) { /* Z < Y < X */
      SceneRotate(G,90,0,1,0); /*0<-->2*/
    } else if((egval[2]<egval[0])&&(egval[0]<egval[1])) { /* Z < X < Y */
      SceneRotate(G,90,0,1,0); /*0<-->2*/
      SceneRotate(G,90,1,0,0); /*0<-->1*/
    }
    /* X < Y < Z  - do nothing - that's what we want */

    ExecutiveWindowZoom(G,sele,0.0,state,0,false);
    if(animate)
      SceneLoadAnimation(G,SettingGetGlobal_f(G,cSetting_animation_duration));

  }
}
/*========================================================================*/
void ExecutiveLabel(PyMOLGlobals *G,char *s1,char *expr,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int cnt;

  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_LABL;
    op1.s1 = expr;
    op1.i1 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    cnt = op1.i1;
    op1.code=OMOP_VISI;
    op1.i1=cRepLabel;
    op1.i2=1;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    op1.code = OMOP_INVA;
    op1.i1=cRepLabel; 
    op1.i2=cRepInvVisib;
    ExecutiveObjMolSeleOp(G,sele1,&op1);

    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Actions)
        " Label: labelled %i atoms.\n",cnt
        ENDFB(G);
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Warnings)
      " Label: no atoms selections.\n"
      ENDFB(G);
  }
}
/*========================================================================*/
int ExecutiveIterate(PyMOLGlobals *G,char *s1,char *expr,int read_only,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRecInit(&op1);
  op1.i1=0;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    op1.code = OMOP_ALTR;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = read_only;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G,FB_Executive,FB_Actions)
          " Alter: modified %i atoms.\n",op1.i1
          ENDFB(G);
      } else {
        PRINTFB(G,FB_Executive,FB_Actions)
          " Iterate: iterated over %i atoms.\n",op1.i1
          ENDFB(G);
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterate: No atoms selected.\n"
        ENDFB(G);
    }
  }
  return(op1.i1);
}
/*========================================================================*/
int ExecutiveSelectList(PyMOLGlobals *G,char *sele_name,
                        char *s1,PyObject *list,int quiet,int id_type)
{/* assumes a blocked Python interpreter */
#ifdef _PYMOL_NOPY
  return -1;
#else
  int ok=true;
  int n_eval=0;
  int sele0 = SelectorIndexByName(G,s1);
  int n_sele = 0;
  ObjectMolecule *obj = NULL;
  if(sele0>=0) obj = SelectorGetSingleObjectMolecule(G,sele0);
  if(obj) {
    int n_atom = obj->NAtom;
    int list_len = 0;
    int a;
    int index = 0;
    int *idx_list = NULL;
    if(ok) ok=PyList_Check(list);
    if(ok) {
      list_len = PyList_Size(list);
      idx_list = Alloc(int,list_len);
      ok = (idx_list!=NULL);
    }
    if(ok) {
      if(list_len) {
        for(a=0;a<list_len;a++) {
          if(ok) 
            ok = PConvPyIntToInt(PyList_GetItem(list,a),&index);
          else
            break;
          if((index<1)||(index>n_atom))
            ok=false;
          else
            idx_list[a]=index-1;
        }
        if(ok) 
          n_sele = SelectorCreateOrderedFromObjectIndices(G,sele_name,obj,idx_list,list_len);
      } else
        SelectorCreateEmpty(G,sele_name);
    }
    FreeP(idx_list);
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)
      " SelectList-Error: selection cannot span more than one object.\n"
      ENDFB(G);
  }
  if(ok) {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Actions)
        " SelectList: modified %i atoms.\n",n_eval
        ENDFB(G);
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterateList: An error occurred.\n"
        ENDFB(G);
    }
  }
  if(!ok)
    return -1;
  else
    return n_sele;
#endif
}


/*========================================================================*/
int ExecutiveIterateList(PyMOLGlobals *G,char *name,PyObject *list,int read_only,int quiet)
{
#ifdef _PYMOL_NOPY
  return -1;
#else
  int ok=true;
  int n_eval=0;
  int sele0 = SelectorIndexByName(G,name);
  PyObject *entry = NULL;
  ObjectMolecule *obj = NULL;
  if(sele0>=0) obj = SelectorGetSingleObjectMolecule(G,sele0);
  if(obj) {
    int n_atom = obj->NAtom;
    int list_len = 0;
    int a;
    int index = 0;
    char *expr = NULL;
    if(ok) ok=PyList_Check(list);
    if(ok) {
      list_len = PyList_Size(list);
      for(a=0;a<list_len;a++) {
        if(ok) entry=PyList_GetItem(list,a);
        if(ok) ok = PyList_Check(entry);
        if(ok) ok = (PyList_Size(entry)==2);
        if(ok) ok = PConvPyIntToInt(PyList_GetItem(entry,0),&index);
        if(ok) ok = PConvPyStrToStrPtr(PyList_GetItem(entry,1),&expr);
        if(ok) ok = ((index<=n_atom) && (index>0));
        if(ok) ok = PAlterAtom(obj->AtomInfo+index-1,expr,read_only,name,index-1);
        if(ok) n_eval++;
      }
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)
      " AlterList-Error: selection cannot span more than one object.\n"
      ENDFB(G);
  }
  if(ok) {
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G,FB_Executive,FB_Actions)
          " AlterList: modified %i atoms.\n",n_eval
          ENDFB(G);
      } else {
        PRINTFB(G,FB_Executive,FB_Actions)
          " IterateList: iterated over %i atoms.\n",n_eval
          ENDFB(G);
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterateList: An error occurred.\n"
        ENDFB(G);
    }
  }
  if(!ok)
    return -1;
  else
    return n_eval;
#endif
}
/*========================================================================*/
void ExecutiveIterateState(PyMOLGlobals *G,int state,char *s1,char *expr,int read_only,
                           int atomic_props,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_AlterState;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = state;
    op1.i3 = read_only;
    op1.i4 = atomic_props;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G,FB_Executive,FB_Actions)
          " AlterState: modified %i atom states.\n",op1.i1
          ENDFB(G);
      } else {
        PRINTFB(G,FB_Executive,FB_Actions)
        " IterateState: iterated over %i atom states.\n",op1.i1
          ENDFB(G);
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterateState: No atoms selected.\n"
        ENDFB(G);
    }
  }
}

typedef struct {
  int priority;
  float vertex[3];
} FitVertexRec;

static int fVertexOrdered(FitVertexRec *array,int l, int r)
{
  return(array[l].priority<=array[r].priority);
}

/*========================================================================*/
float ExecutiveRMS(PyMOLGlobals *G,char *s1,char *s2,int mode,float refine,int max_cyc,
                   int quiet,char *oname,int state1,int state2,
                   int ordered_selections)
{
  int sele1,sele2;
  float rms = -1.0;
  int a,b;
  float inv,*f,*f1,*f2;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType buffer;
  int *flag;
  int ok=true;
  int repeat;
  float v1[3],*v2;
            
  sele1=SelectorIndexByName(G,s1);

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.vv1=NULL;
  op1.vc1=NULL;
  op2.vv1=NULL;
  op2.vc1=NULL;

  
  if(sele1>=0) {
    if(state1<0) {
      op1.code = OMOP_AVRT;
    } else {
      op1.code = OMOP_StateVRT;
      op1.i1=state1;
    }
    op1.nvv1=0;
    op1.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    if(ordered_selections)
      op1.vp1=VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    for(a=0;a<op1.nvv1;a++)
      {
        inv=(float)op1.vc1[a]; /* average over coordinate sets */
        if(inv)
          {
            f=op1.vv1+(a*3);
            inv=1.0F/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }
  
  sele2=SelectorIndexByName(G,s2);
  if(sele2>=0) {

    if(state2<0) {
      op2.code = OMOP_AVRT;
    } else {
      op2.code = OMOP_StateVRT;
      op2.i1=state2;
    }
    op2.nvv1=0;
    op2.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    if(ordered_selections)
      op2.vp1=VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele2,&op2);
    for(a=0;a<op2.nvv1;a++)
      {
        inv=(float)op2.vc1[a]; /* average over coordinate sets */
        if(inv)
          {
            f=op2.vv1+(a*3);
            inv=1.0F/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }

  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1!=op2.nvv1) {
      sprintf(buffer,"Atom counts between selections don't match (%d vs %d)",
              op1.nvv1,op2.nvv1);
      ErrMessage(G,"ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(!SelectorGetSingleObjectMolecule(G,sele1)) {
        if(mode!=2) {
          PRINTFB(G,FB_Executive,FB_Warnings)
            "Executive-Warning: Mobile selection spans more than one object.\n"
            ENDFB(G);
        } else {
          PRINTFB(G,FB_Executive,FB_Errors)
            "Executive-Error: Mobile selection spans more than one object. Aborting.\n"
            ENDFB(G);
          ok=false;
        }
      }

      if(ordered_selections&&op1.vp1&&op2.vp1) {
        /* if we expected ordered selections and have priorities, 
           then we may need to sort vertices */

        int sort_flag1 = false, sort_flag2 = false;
        int well_defined1 = true, well_defined2 = true;
        
        for(a=0;a<(op1.nvv1-1);a++) {
          /*          printf("op1 vertex %d priority %d\n",a,op1.vp1[a]);
                      printf("op2 vertex %d priority %d\n",a,op2.vp1[a]);*/

          if(op1.vp1[a]>op1.vp1[a+1])
            sort_flag1 = true;
          else if(op1.vp1[a]==op1.vp1[a+1])
            well_defined1 = false;
          if(op2.vp1[a]>op2.vp1[a+1])
            sort_flag2 = true;
          else if(op2.vp1[a]==op2.vp1[a+1])
            well_defined2 = false;
        }
        
        if(sort_flag1||sort_flag2) {
          if(!(well_defined1||well_defined2)) {
            PRINTFB(G,FB_Executive,FB_Warnings) 
              "Executive-Warning: Ordering requested but not well defined.\n"
               ENDFB(G);
          } else {
            FitVertexRec *vert = Alloc(FitVertexRec,op1.nvv1);

            if(sort_flag1) {
              float *src,*dst;
              src = op1.vv1;
              for(a=0;a<op1.nvv1;a++) {              
                vert[a].priority = op1.vp1[a];
                dst=vert[a].vertex;
                copy3f(src,dst);
                src+=3;
              }
              UtilSortInPlace(G,vert,op1.nvv1,sizeof(FitVertexRec),(UtilOrderFn*)fVertexOrdered);
              dst = op1.vv1;
              for(a=0;a<op1.nvv1;a++) {              
                src=vert[a].vertex;
                copy3f(src,dst);
                dst+=3;
              }
            }

            if(sort_flag2) {
              float *src,*dst;
              src = op2.vv1;
              for(a=0;a<op2.nvv1;a++) {              
                vert[a].priority = op2.vp1[a];
                dst=vert[a].vertex;
                copy3f(src,dst);
                src+=3;
              }
              UtilSortInPlace(G,vert,op2.nvv1,sizeof(FitVertexRec),(UtilOrderFn*)fVertexOrdered);
              dst = op2.vv1;
              for(a=0;a<op2.nvv1;a++) {              
                src=vert[a].vertex;
                copy3f(src,dst);
                dst+=3;
              }
            }
            
            FreeP(vert);
          }
        }
      }
      if(mode!=0) {
        rms = MatrixFitRMS(G,op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
        repeat=true;
        b=0;
        while(repeat) {
          repeat=false;
          b++;
          if(b>max_cyc)
            break;
          if((refine>R_SMALL4)&&(rms>R_SMALL4)) {
            flag=Alloc(int,op1.nvv1);
            
            if(flag) {          
              for(a=0;a<op1.nvv1;a++) {
                MatrixApplyTTTfn3f(1,v1,op2.ttt,op1.vv1+(a*3));
                v2=op2.vv1+(a*3);
                if((diff3f(v1,v2)/rms)>refine) {
                  flag[a] = false;
                  repeat=true;
                }
                else
                  flag[a] = true;
              }
              f1 = op1.vv1;
              f2 = op2.vv1;
              for(a=0;a<op1.nvv1;a++) {
                if(!flag[a]) {
                  op2.nvv1--;
                } else {
                  copy3f(op1.vv1+(3*a),f1);
                  copy3f(op2.vv1+(3*a),f2);
                  f1+=3;
                  f2+=3;
                }
              }
              if(op2.nvv1!=op1.nvv1) {
                PRINTFB(G,FB_Executive,FB_Actions)
                  " ExecutiveRMS: %d atoms rejected during cycle %d (RMS=%0.2f).\n",op1.nvv1-op2.nvv1,b,rms
                  ENDFB(G);
              }
              op1.nvv1 = op2.nvv1;
              FreeP(flag);
              if(op1.nvv1) 
                rms = MatrixFitRMS(G,op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);            
              else
                break;
            }
          }
        }
      }
      else
        rms = MatrixGetRMS(G,op1.nvv1,op1.vv1,op2.vv1,NULL);

      if(!op1.nvv1) {
        PRINTFB(G,FB_Executive,FB_Results) 
          " Executive: Error -- no atoms left after refinement!\n"
          ENDFB(G);
        ok=false;
      }

      if(ok) {
        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Results) 
            " Executive: RMS = %8.3f (%d to %d atoms)\n", rms,op1.nvv1,op2.nvv1 
            ENDFB(G);
        }
        if(oname) 
          if(oname[0]) {
            CGO *cgo = NULL;
            ObjectCGO *ocgo;
            int auto_save;

            cgo=CGONew(G);
            /*             CGOColor(cgo,1.0,1.0,0.0); 
                           CGOLinewidth(cgo,3.0);*/
            CGOBegin(cgo,GL_LINES);
            for(a=0;a<op1.nvv1;a++) {
              CGOVertexv(cgo,op2.vv1+(a*3));
              MatrixApplyTTTfn3f(1,v1,op2.ttt,op1.vv1+(a*3));
              CGOVertexv(cgo,v1);
            }
            CGOEnd(cgo);
            CGOStop(cgo);
            ocgo = ObjectCGOFromCGO(G,NULL,cgo,0);
            ocgo->Obj.Color = ColorGetIndex(G,"yellow");
            ObjectSetName((CObject*)ocgo,oname);
            ExecutiveDelete(G,oname);
            auto_save = (int)SettingGet(G,cSetting_auto_zoom);
            SettingSet(G,cSetting_auto_zoom,0);
            ExecutiveManageObject(G,(CObject*)ocgo,true,false);
            SettingSet(G,cSetting_auto_zoom,(float)auto_save);            
            SceneDirty(G);
          }
        if(mode==2) {
          if(ok) {
            op2.code = OMOP_TTTF;
            ExecutiveObjMolSeleOp(G,sele1,&op2);
          }
        }
      }
    } else {
      ErrMessage(G,"ExecutiveRMS","No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  VLAFreeP(op1.vp1);
  VLAFreeP(op2.vp1);
  return(rms);
}
/*========================================================================*/
int *ExecutiveIdentify(PyMOLGlobals *G,char *s1,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  int *result = NULL;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_Identify;
    op2.i1=0;
    op2.i1VLA=VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    result = op2.i1VLA;
    VLASize(result,int,op2.i1);
  } 
  return(result);
}
/*========================================================================*/
int ExecutiveIdentifyObjects(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_IdentifyObjects;
    op2.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op2.i1VLA=VLAlloc(int,1000);
    op2.i1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    VLASize(op2.i1VLA,int,op2.i1);
    VLASize(op2.obj1VLA,ObjectMolecule*,op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } 
  return(op2.i1);
}
/*========================================================================*/
int ExecutiveIndex(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_Index;
    op2.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op2.i1VLA=VLAlloc(int,1000);
    op2.i1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    VLASize(op2.i1VLA,int,op2.i1);
    VLASize(op2.obj1VLA,ObjectMolecule*,op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } 
  return(op2.i1);
}
/*========================================================================*/
float *ExecutiveRMSStates(PyMOLGlobals *G,char *s1,int target,int mode,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  float *result = NULL;
  int ok=true;
  
  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.vv1=NULL;
  op2.vv1=NULL;
  sele1=SelectorIndexByName(G,s1);
  
  if(!SelectorGetSingleObjectMolecule(G,sele1)) {
    if(mode!=2) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "Executive-Warning: Mobile selection spans more than one object.\n"
        ENDFB(G);
    } else {
      PRINTFB(G,FB_Executive,FB_Errors)
        "Executive-Error: Mobile selection spans more than one object. Aborting.\n\n"
        ENDFB(G);
      ok=false;
    }
  }

  if(ok&&sele1>=0) {
    op1.code = OMOP_SVRT;
    op1.nvv1=0;
    op1.i1=target;
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    op1.i1VLA = VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele1,&op1);

    op2.vv2=op1.vv1;
    op2.nvv2=op1.nvv1;
    op2.i1VLA=op1.i1VLA;
    op2.i2=target;
    op2.i1=mode;
    op2.f1VLA=VLAlloc(float,10);
    VLASize(op2.f1VLA,float,0); /* failsafe */
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    op2.code = OMOP_SFIT;
    op2.nvv1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    result=op2.f1VLA;
    VLAFreeP(op1.vv1);
    VLAFreeP(op1.i1VLA);
    VLAFreeP(op2.vv1);
  } 
  return(result);
}
/*========================================================================*/
float ExecutiveRMSPairs(PyMOLGlobals *G,WordType *sele,int pairs,int mode)
{
  int sele1,sele2;
  int a,c;
  float rms=0.0,inv,*f;
  OrthoLineType buffer;

  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType combi,s1;

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.nvv1=0;
  op1.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
  op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1); /* auto-zero */
  op1.code = OMOP_AVRT;

  op2.nvv1=0;
  op2.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
  op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1); /* auto-zero */
  op2.code = OMOP_AVRT;

  strcpy(combi,"(");
  c=0;
  for(a=0;a<pairs;a++) {
    sele1=SelectorIndexByName(G,sele[c]);
    if(sele1>=0) ExecutiveObjMolSeleOp(G,sele1,&op1);
    strcat(combi,sele[c]);
    if(a<(pairs-1)) strcat(combi," or ");
    c++;
    sele2=SelectorIndexByName(G,sele[c]);
    if(sele2>=0) ExecutiveObjMolSeleOp(G,sele2,&op2);
    c++;
  }
  strcat(combi,")");
  for(a=0;a<op1.nvv1;a++)
    {
      inv=(float)op1.vc1[a];
      if(inv)
        {
          f=op1.vv1+(a*3);
          inv=1.0F/inv;
          *(f++)*=inv;
          *(f++)*=inv;
          *(f++)*=inv;
        }
    }
  for(a=0;a<op2.nvv1;a++)
    {
      inv=(float)op2.vc1[a];
      if(inv)
        {
          f=op2.vv1+(a*3);
          inv=1.0F/inv;
          *(f++)*=inv;
          *(f++)*=inv;
          *(f++)*=inv;
        }
    }
  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1!=op2.nvv1) {
      sprintf(buffer,"Atom counts between selection sets don't match (%d != %d).",
              op1.nvv1,op2.nvv1);
      ErrMessage(G,"ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(mode!=0)
        rms = MatrixFitRMS(G,op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
      else
        rms = MatrixGetRMS(G,op1.nvv1,op1.vv1,op2.vv1,NULL);
      PRINTFB(G,FB_Executive,FB_Results) 
        " ExecutiveRMS: RMS = %8.3f (%d to %d atoms)\n",
        rms,op1.nvv1,op2.nvv1
        ENDFB(G);
    
      op2.code = OMOP_TTTF;
      SelectorGetTmp(G,combi,s1);
      sele1=SelectorIndexByName(G,s1);
      ExecutiveObjMolSeleOp(G,sele1,&op2);
      SelectorFreeTmp(G,s1);
    } else {
      ErrMessage(G,"ExecutiveRMS","No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  return(rms);
}
/*========================================================================*/
void ExecutiveUpdateObjectSelection(PyMOLGlobals *G,struct CObject *obj)
{
  if(obj->type==cObjectMolecule) {
    SelectorUpdateObjectSele(G,(ObjectMolecule*)obj);  
  }
}
/*========================================================================*/
int ExecutiveReset(PyMOLGlobals *G,int cmd,char *name)
{
  int ok=true;
  CObject *obj;
  if(!name[0]) {
    SceneResetMatrix(G);
    ExecutiveWindowZoom(G,cKeywordAll,0.0,-1,0,0); /* reset does all states */
  } else {
    obj = ExecutiveFindObjectByName(G,name);
    if(!obj)
      ok=false;
    else
      ObjectResetTTT(obj);
  }
  return(ok);
}
/*========================================================================*/
void ExecutiveDrawNow(PyMOLGlobals *G) 
{
  PRINTFD(G,FB_Executive)
    " ExecutiveDrawNow: entered.\n"
    ENDFD;

  OrthoExecDeferred(G);

  if(!SettingGet(G,cSetting_suspend_updates)) {

    if(G->HaveGUI && G->ValidContext) {
      glMatrixMode(GL_MODELVIEW); /* why is this necessary?  is it? */
    }

    SceneUpdate(G);
    if(WizardUpdate(G))
      SceneUpdate(G);

    OrthoDoDraw(G);
    
    PyMOL_NeedSwap(G->PyMOL);
  }

  PRINTFD(G,FB_Executive)
    " ExecutiveDrawNow: leaving.\n"
    ENDFD;
}
/*========================================================================*/
int ExecutiveCountStates(PyMOLGlobals *G,char *s1)
{
  register CExecutive *I = G->Executive;
  int sele1;
  int result=0;
  int n_frame;
  SpecRec *rec = NULL;
  
  if(s1)
    if(WordMatch(G,cKeywordAll,s1,true))
      s1 = NULL;
  if(!s1) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->fGetNFrame) {
          n_frame = rec->obj->fGetNFrame(rec->obj);
          if(result<n_frame)
            result=n_frame;
        }
      }
    } 
  } else {
  sele1=SelectorIndexByName(G,s1);
    if(sele1>=0) {
      SelectorUpdateTable(G);
      result = SelectorGetSeleNCSet(G,sele1);
    }
  }
  return(result);
}
/*========================================================================*/
void ExecutiveRay(PyMOLGlobals *G,int width,int height,int mode,float angle,float shift,int quiet)
{
  SceneRay(G,width,height,mode,NULL,NULL,angle,shift,quiet);
}
/*========================================================================*/
int  ExecutiveSetSetting(PyMOLGlobals *G,int index,PyObject *tuple,char *sele,
                         int state,int quiet,int updates)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  OrthoLineType value;
  CSetting **handle=NULL;
  SettingName name;
  int nObj=0;
  int unblock;
  int ok =true;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetSetting: entered. sele \"%s\"\n",sele
    ENDFD;
  unblock = PAutoBlock();
  if(sele[0]==0) { 
    ok = SettingSetTuple(G,NULL,index,tuple);
    if(ok) {
      if(!quiet) {
        if(Feedback(G,FB_Setting,FB_Actions)) {
          SettingGetTextValue(G,NULL,NULL,index,value);
          SettingGetName(G,index,name);
          PRINTF
            " Setting: %s set to %s.\n",name,value
            ENDF(G);
        }
      }
      if(updates) 
        SettingGenerateSideEffects(G,index,sele,state);
    }
  } 
  else if(!strcmp(cKeywordAll,sele)) { /* all objects setting */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject) {
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(G,handle);
              ok = SettingSetTuple(G,*handle,index,tuple);
              nObj++;
            }
          }
        }
        if(nObj) {
          if(updates) 
            SettingGenerateSideEffects(G,index,sele,state);
        }
        if(Feedback(G,FB_Setting,FB_Actions)) {
          if(nObj&&handle) {
            SettingGetTextValue(G,*handle,NULL,index,value);
            SettingGetName(G,index,name);
            if(!quiet) {
              if(state<0) {
                PRINTF
                  " Setting: %s set to %s in %d objects.\n",name,value,nObj
                  ENDF(G);
              } else {
                PRINTF
                  " Setting: %s set to %s in %d objects, state %d.\n",
                  name,value,nObj,state+1
                  ENDF(G);
              }
            }
          }
        }
      }
  } else { /* based on a selection/object name */
    sele1=SelectorIndexByName(G,sele);
    while((ListIterate(I->Spec,rec,next)))
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule)
          {
            if(sele1>=0) {
              obj=(ObjectMolecule*)rec->obj;
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_CountAtoms;
              op.i1=0;
              ObjectMoleculeSeleOp(obj,sele1,&op);
              if(op.i1&&rec->obj->fGetSettingHandle) {
                handle = rec->obj->fGetSettingHandle(rec->obj,state);
                if(handle) {
                  SettingCheckHandle(G,handle);
                  ok = SettingSetTuple(G,*handle,index,tuple);
                  if(ok) {
                    if(updates) 
                      SettingGenerateSideEffects(G,index,sele,state);
                    if(!quiet) {
                      if(state<0) { /* object-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetTextValue(G,*handle,NULL,index,value);
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\".\n",
                            name,value,rec->obj->Name
                            ENDF(G);
                        }
                      } else { /* state-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetTextValue(G,*handle,NULL,index,value);
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\", state %d.\n",
                            name,value,rec->obj->Name,state+1
                            ENDF(G);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if(strcmp(rec->obj->Name,sele)==0) {
            if(rec->obj->fGetSettingHandle) {
              handle = rec->obj->fGetSettingHandle(rec->obj,state);
              if(handle) {
                SettingCheckHandle(G,handle);
                ok = SettingSetTuple(G,*handle,index,tuple);
                if(ok) {
                  if(updates)
                    SettingGenerateSideEffects(G,index,sele,state);
                  if(!quiet) {
                    if(state<0) { /* object-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetTextValue(G,*handle,NULL,index,value);
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\".\n",
                          name,value,rec->obj->Name
                          ENDF(G);
                      }
                    } else { /* state-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetTextValue(G,*handle,NULL,index,value);
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\", state %d.\n",
                          name,value,rec->obj->Name,state+1
                          ENDF(G);
                      }
                    }
                  }
                }
              }
            }
          }
      }
  }
  PAutoUnblock(unblock);
  return(ok);
#endif
}
/*========================================================================*/
int  ExecutiveUnsetSetting(PyMOLGlobals *G,int index,char *sele,
                         int state,int quiet,int updates)
{
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  CSetting **handle=NULL;
  SettingName name;
  int nObj=0;
  int unblock;
  int ok =true;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetSetting: entered. sele \"%s\"\n",sele
    ENDFD;
  unblock = PAutoBlock();
  if(sele[0]==0) { 
    /* do nothing */
  } 
  else if(!strcmp(cKeywordAll,sele)) { /* all objects setting */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject) {
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(G,handle);
              ok = SettingUnset(*handle,index);
              nObj++;
            }
          }
        }
        if(nObj) {
          if(updates) 
            SettingGenerateSideEffects(G,index,sele,state);
        }
        if(Feedback(G,FB_Setting,FB_Actions)) {
          if(nObj&&handle) {
            SettingGetName(G,index,name);
            if(!quiet) {
              if(state<0) {
                PRINTF
                  " Setting: %s unset in %d objects.\n",name,nObj
                  ENDF(G);
              } else {
                PRINTF
                  " Setting: %s unset in %d objects, state %d.\n",
                  name,nObj,state+1
                  ENDF(G);
              }
            }
          }
        }
      }
  } else { /* based on a selection/object name */
    sele1=SelectorIndexByName(G,sele);
    while((ListIterate(I->Spec,rec,next)))
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule)
          {
            if(sele1>=0) {
              obj=(ObjectMolecule*)rec->obj;
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_CountAtoms;
              op.i1=0;
              ObjectMoleculeSeleOp(obj,sele1,&op);
              if(op.i1&&rec->obj->fGetSettingHandle) {
                handle = rec->obj->fGetSettingHandle(rec->obj,state);
                if(handle) {
                  SettingCheckHandle(G,handle);
                  ok = SettingUnset(*handle,index);
                  if(ok) {
                    if(updates) 
                      SettingGenerateSideEffects(G,index,sele,state);
                    if(!quiet) {
                      if(state<0) { /* object-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s unset in object \"%s\".\n",
                            name,rec->obj->Name
                            ENDF(G);
                        }
                      } else { /* state-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s unset in object \"%s\", state %d.\n",
                            name,rec->obj->Name,state+1
                            ENDF(G);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if(strcmp(rec->obj->Name,sele)==0) {
            if(rec->obj->fGetSettingHandle) {
              handle = rec->obj->fGetSettingHandle(rec->obj,state);
              if(handle) {
                SettingCheckHandle(G,handle);
                ok = SettingUnset(*handle,index);
                if(ok) {
                  if(updates)
                    SettingGenerateSideEffects(G,index,sele,state);
                  if(!quiet) {
                    if(state<0) { /* object-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s unset in object \"%s\".\n",
                          name,rec->obj->Name
                          ENDF(G);
                      }
                    } else { /* state-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s unset in object \"%s\", state %d.\n",
                          name,rec->obj->Name,state+1
                          ENDF(G);
                      }
                    }
                  }
                }
              }
            }
          }
      }
  }
  PAutoUnblock(unblock);
  return(ok);
}
/*========================================================================*/
int ExecutiveColor(PyMOLGlobals *G,char *name,char *color,int flags,int quiet)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int sele;
  ObjectMoleculeOpRec op;
  int col_ind;
  int ok=false;
  int n_atm=0;
  int n_obj=0;
  char atms[]="s";
  char objs[]="s";
  char *best_match;
  col_ind = ColorGetIndex(G,color);
  if(col_ind==-1) {
    ErrMessage(G,"Color","Unknown color.");
  } else {
    best_match = ExecutiveFindBestNameMatch(G,name);
    /* per atom */
    if(!(flags&0x1)) {
      sele=SelectorIndexByName(G,name);
      if(sele>=0) {
        ok=true; 
        ObjectMoleculeOpRecInit(&op);
        op.code = OMOP_COLR;
        op.i1= col_ind;
        op.i2= 0;
        ExecutiveObjMolSeleOp(G,sele,&op);
        n_atm = op.i2;
        op.code=OMOP_INVA;
        op.i1=cRepAll; 
        op.i2=cRepInvColor;
        ExecutiveObjMolSeleOp(G,sele,&op);
      }
    }
    /* per object */
    if(strcmp(name,cKeywordAll)) {
      rec=ExecutiveFindSpec(G,name);
      if(rec) {
        if(rec->type==cExecObject) {
          rec->obj->Color=col_ind;
          if(rec->obj->fInvalidate)
            rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,cRepAll);
          n_obj++;
          ok=true;
          SceneDirty(G);
        }
      } 
    } else {
      rec=NULL;
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          rec->obj->Color=col_ind;
          if(rec->obj->fInvalidate)
            rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,cRepAll);
          n_obj++;
          ok=true;
          SceneDirty(G);
        }
      }
    }
    if(n_obj||n_atm) {
      if(n_obj<2) objs[0]=0;
      if(n_atm<2) atms[0]=0;
      if(!quiet) {
        if(n_obj&&n_atm) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Executive: Colored %d atom%s and %d object%s.\n",n_atm,atms,n_obj,objs
            ENDFB(G);
        } else if (n_obj) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Executive: Colored %d object%s.\n",n_obj,objs
            ENDFB(G);
        } else {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Executive: Colored %d atom%s.\n",n_atm,atms
            ENDFB(G);
        }
      }
    }
  }
  return(ok);
}
/*========================================================================*/
char *ExecutiveFindBestNameMatch(PyMOLGlobals *G,char *name)
{
  char *result;
  register CExecutive *I = G->Executive;
  SpecRec *rec=NULL,*best_rec = NULL;
  int best;
  int wm;

  best = 0;
  result = name;

  while(ListIterate(I->Spec,rec,next)) {
    wm = WordMatch(G,name,rec->name,true);
    if(wm<0) {
      best_rec = rec;
      best = wm;
      break;
    } else if ((best>0)&&(best<wm)) {
      best_rec=rec;
      best = wm;
    }
  }
  if(best_rec)
    result=best_rec->name;
  return(result);
}
/*========================================================================*/
static SpecRec *ExecutiveFindSpec(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next)) {
	 if(strcmp(rec->name,name)==0) 
		break;
  }
  return(rec);
}
/*========================================================================*/
void ExecutiveObjMolSeleOp(PyMOLGlobals *G,int sele,ObjectMoleculeOpRec *op) 
{
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;

  if(sele>=0)
	 {
		while(ListIterate(I->Spec,rec,next))
		  {
			 if(rec->type==cExecObject)
				{
				  if(rec->obj->type==cObjectMolecule)
					 {
						obj=(ObjectMolecule*)rec->obj;
						ObjectMoleculeSeleOp(obj,sele,op);
					 }
				}
		  }
	 }
}

/*========================================================================*/
int ExecutiveGetCameraExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,int transformed,int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int flag = false;

  if(state==-2) state=SceneGetState(G);

  PRINTFD(G,FB_Executive)
    " ExecutiveGetCameraExtent: name %s state %d\n",name,state
    ENDFD;
  
  sele=SelectorIndexByName(G,name);

  if(sele>=0) { 
    ObjectMoleculeOpRecInit(&op);
    if(state<0) {
      op.code = OMOP_CameraMinMax;
    } else {
      op.code = OMOP_CSetCameraMinMax;
      op.cs1 = state;
    }
	 op.v1[0]=FLT_MAX;
	 op.v1[1]=FLT_MAX;
	 op.v1[2]=FLT_MAX;
    op.v2[0]=FLT_MIN;
    op.v2[1]=FLT_MIN;
    op.v2[2]=FLT_MIN;
    op.i1 = 0;
    op.i2 = transformed;
    op.mat1=SceneGetMatrix(G);

    ExecutiveObjMolSeleOp(G,sele,&op);

    PRINTFD(G,FB_Executive)
      " ExecutiveGetCameraExtent: minmax over %d vertices\n",op.i1
      ENDFD;
    if(op.i1)
      flag = true;
  }
  copy3f(op.v1,mn);
  copy3f(op.v2,mx);
  
  PRINTFD(G,FB_Executive)
    " ExecutiveGetCameraExtent: returning %d\n",flag
    ENDFD;

  return(flag);  
}

/*========================================================================*/
int ExecutiveGetExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,int transformed,int state,int weighted)
{
  int sele;
  ObjectMoleculeOpRec op,op2;
  register CExecutive *I=G->Executive;
  CObject *obj;
  int flag = false;
  SpecRec *rec = NULL;
  int all_flag = false;
  float f1,f2,fmx;
  int a;

  if(WordMatch(G,cKeywordCenter,name,1)<0) {
    SceneGetPos(G,mn);
    copy3f(mn,mx);
    return 1;
  }
  if(WordMatch(G,cKeywordOrigin,name,1)<0) {
    SceneOriginGet(G,mn);
    copy3f(mn,mx);
    return 1;
  }
  if(state==-2) state=SceneGetState(G);

  PRINTFD(G,FB_Executive)
    " ExecutiveGetExtent: name %s state %d\n",name,state
    ENDFD;

  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);  
  op2.i1 = 0;
  op2.v1[0]=-1.0;
  op2.v1[1]=-1.0;
  op2.v1[2]=-1.0;
  op2.v2[0]=1.0;
  op2.v2[1]=1.0;
  op2.v2[2]=1.0;
  
  if(WordMatch(G,cKeywordAll,name,true)<0) {
    all_flag=true;
  }
  sele=SelectorIndexByName(G,name);

  if(sele>=0) { 
    if(state<0) {
      op.code = OMOP_MNMX;
    } else {
      op.code = OMOP_CSetMinMax;
      op.cs1 = state;
    }
	 op.v1[0]=FLT_MAX;
	 op.v1[1]=FLT_MAX;
	 op.v1[2]=FLT_MAX;
    op.v2[0]=FLT_MIN;
    op.v2[1]=FLT_MIN;
    op.v2[2]=FLT_MIN;
    op.i1 = 0;
    op.i2 = transformed;
    ExecutiveObjMolSeleOp(G,sele,&op);

    PRINTFD(G,FB_Executive)
      " ExecutiveGetExtent: minmax over %d vertices\n",op.i1
      ENDFD;

    if(op.i1)
      flag = true;
    if(all_flag) {
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          obj=rec->obj;
          if(obj->ExtentFlag) 
            switch(obj->type) {
            case cObjectMolecule:
              break;
            default:
              min3f(obj->ExtentMin,op.v1,op.v1);
              max3f(obj->ExtentMax,op.v2,op.v2);
              flag = true;
              break;
            }
        }
      }
    }
    if(weighted) {
      op2.i1=0;
      op2.i2=transformed;
      if(state<0) 
        op2.code = OMOP_SUMC;
      else {
        op2.code = OMOP_CSetSumVertices;
        op2.cs1 = state;
      }
      
      op2.v1[0]=0.0;
      op2.v1[1]=0.0;
      op2.v1[2]=0.0;
      ExecutiveObjMolSeleOp(G,sele,&op2);
      if(op2.i1) {
        op2.v1[0]/=op2.i1;
        op2.v1[1]/=op2.i1;
        op2.v1[2]/=op2.i1;
      }
    }
  } else {
    obj = ExecutiveFindObjectByName(G,name);
    if(obj) {
      switch(obj->type) {
      case cObjectMolecule:
        break;
      default:
        if(obj->ExtentFlag) {
          copy3f(obj->ExtentMin,op.v1);
          copy3f(obj->ExtentMax,op.v2);
          flag = true;
        } else {

          PRINTFD(G,FB_Executive)
            " ExecutiveGetExtent: no extent for object %s\n",obj->Name
            ENDFD;
          
        }
        break;
      }
    }
  }
  if(flag&&weighted) { 
    if(op2.i1) { 
      for (a=0;a<3;a++) { /* this puts origin at the weighted center */
        f1 = op2.v1[a] - op.v1[a];
        f2 = op.v2[a] - op2.v1[a];
        if(f1>f2) 
          fmx = f1;
        else
          fmx = f2;
        op.v1[a] = op2.v1[a] - fmx;
        op.v2[a] = op2.v1[a] + fmx;
      }
    }
  }
  copy3f(op.v1,mn);
  copy3f(op.v2,mx);
  
  if(all_flag) {
    rec=NULL;
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        obj=rec->obj;
        switch(rec->obj->type) {
        case cObjectMolecule:
          break;
        default:
          if(obj->ExtentFlag) {
            if(!flag) {
              copy3f(obj->ExtentMax,mx);
              copy3f(obj->ExtentMin,mn);
              flag=true;
            } else {
              max3f(obj->ExtentMax,mx,mx);
              min3f(obj->ExtentMin,mn,mn);
            }
          }
          break;
        }
      }
    }
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveGetExtent: returning %d\n",flag
    ENDFD;

  return(flag);  
}
/*========================================================================*/
static int ExecutiveGetMaxDistance(PyMOLGlobals *G,char *name,float *pos,float *dev,int transformed,int state)
{
  int sele;
  ObjectMoleculeOpRec op,op2;
  register CExecutive *I=G->Executive;
  CObject *obj;
  int flag = false;
  SpecRec *rec = NULL;
  int all_flag = false;
  float f1,fmx=0.0F;

  if(state==-2) state=SceneGetState(G);

  PRINTFD(G,FB_Executive)
    " ExecutiveGetExtent: name %s state %d\n",name,state
    ENDFD;
  
  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);

  op2.i1 = 0;
  op2.v1[0]=-1.0;
  op2.v1[1]=-1.0;
  op2.v1[2]=-1.0;
  op2.v2[0]=1.0;
  op2.v2[1]=1.0;
  op2.v2[2]=1.0;
  
  if(WordMatch(G,cKeywordAll,name,true)<0) {
      all_flag=true;
  }
  sele=SelectorIndexByName(G,name);

  if(sele>=0) {
    if(state<0) {
      op.code = OMOP_MaxDistToPt;
    } else {
      op.code = OMOP_CSetMaxDistToPt;
      op.cs1 = state;
    }
	 op.v1[0]=pos[0];
	 op.v1[1]=pos[1];
	 op.v1[2]=pos[2];
    op.i1 = 0;
    op.f1 = 0.0F;
    op.i2 = transformed;
    ExecutiveObjMolSeleOp(G,sele,&op);
    fmx = op.f1;

    if(op.i1)
      flag = true;
    if(all_flag) {
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          obj=rec->obj;
          if(obj->ExtentFlag) 
            switch(obj->type) {
            case cObjectMolecule:
              break;
            default:
              f1 = (float)diff3f(obj->ExtentMin,pos);
              if(fmx<f1) fmx = f1;
              f1 = (float)diff3f(obj->ExtentMax,pos);
              if(fmx<f1) fmx = f1;
              flag = true;
              break;
            }
        }
      }
    }
  } else {
    obj = ExecutiveFindObjectByName(G,name);
    if(obj) {
      switch(obj->type) {
      case cObjectMolecule:
        break;
      default:
        if(obj->ExtentFlag) {
          f1 = (float)diff3f(obj->ExtentMin,pos);
          if(fmx<f1) fmx = f1;
          f1 = (float)diff3f(obj->ExtentMax,pos);
          if(fmx<f1) fmx = f1;
          flag = true;
          break;
        }
      }
    } else if(all_flag) {
      rec=NULL;
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          obj=rec->obj;
          switch(rec->obj->type) {
          case cObjectMolecule:
            break;
          default:
            if(obj->ExtentFlag) {
              f1 = (float)diff3f(obj->ExtentMin,pos);
              if(fmx<f1) fmx = f1;
              f1 = (float)diff3f(obj->ExtentMax,pos);
              if(fmx<f1) fmx = f1;
            }
            break;
          }
        }
      }
    }
  }
  *dev = fmx;
  return(flag);  
}
/*========================================================================*/
int ExecutiveWindowZoom(PyMOLGlobals *G,char *name,float buffer,
                        int state,int inclusive,int animate)
{
  float center[3],radius;
  float mn[3],mx[3],df[3];
  int sele0;
  int ok=true;

  PRINTFD(G,FB_Executive)
    " ExecutiveWindowZoom-DEBUG: entered\n"
    ENDFD;
  if(ExecutiveGetExtent(G,name,mn,mx,true,state,true)) {
    if(buffer!=0.0) {
      buffer = buffer;
      mx[0]+=buffer;
      mx[1]+=buffer;
      mx[2]+=buffer;
      mn[0]-=buffer;
      mn[1]-=buffer;
      mn[2]-=buffer;
    }
    subtract3f(mx,mn,df);
    average3f(mn,mx,center);
    if(inclusive) {
      if(!ExecutiveGetMaxDistance(G,name,center,&radius,true,state))
        radius=0.0;
      radius+=buffer;
    } else {
      radius = df[0];
      if(radius<df[1]) radius=df[1];
      if(radius<df[2]) radius=df[2];
      radius=radius/2.0F;
    }
    if(radius<MAX_VDW) radius=MAX_VDW;
    PRINTFD(G,FB_Executive)
      " ExecutiveWindowZoom: zooming with radius %8.3f...state %d\n",radius,state
      ENDFD;
    PRINTFD(G,FB_Executive)
      " ExecutiveWindowZoom: on center %8.3f %8.3f %8.3f...\n",center[0],
      center[1],center[2]
      ENDFD;
    if(animate<0)
      animate=SettingGetGlobal_b(G,cSetting_animation);
    if(animate)
      ScenePrimeAnimation(G);
    SceneOriginSet(G,center,false);
    SceneWindowSphere(G,center,radius);
    if(animate)
      SceneLoadAnimation(G,SettingGetGlobal_f(G,cSetting_animation_duration));
    SceneDirty(G);
  } else {

    sele0 = SelectorIndexByName(G,name);
    if(sele0>0) { /* any valid selection except "all" */
      ErrMessage(G,"ExecutiveWindowZoom","selection doesn't specify any coordinates.");
      ok=false;
    } else if(ExecutiveValidName(G,name)) {
      PRINTFD(G,FB_Executive)
        " ExecutiveWindowZoom-DEBUG: name valid, but no extents -- using default view\n"
        ENDFD;
      SceneSetDefaultView(G);
      SceneDirty(G);
    } else {
      ErrMessage(G,"ExecutiveWindowZoom","selection or object unknown.");
      ok=false;
    }
  }

  return(ok);
}
/*========================================================================*/
int ExecutiveCenter(PyMOLGlobals *G,char *name,int state,
                    int origin,int animate)
{
  float center[3];
  float mn[3],mx[3],df[3];
  int sele0;
  int ok=true;

  if(ExecutiveGetExtent(G,name,mn,mx,true,state,true)) {
    subtract3f(mx,mn,df);
    average3f(mn,mx,center);
    PRINTFD(G,FB_Executive)
      " ExecutiveCenter: centering state %d\n",state
      ENDFD;
    PRINTFD(G,FB_Executive)
      " ExecutiveCenter: on center %8.3f %8.3f %8.3f...\n",center[0],
      center[1],center[2]
      ENDFD;

    if(animate<0)
      animate=SettingGetGlobal_b(G,cSetting_animation);

    if(animate)
      ScenePrimeAnimation(G);
    if(origin) 
      SceneOriginSet(G,center,false);
    SceneRelocate(G,center);
    SceneDirty(G);
    if(animate)
      SceneLoadAnimation(G,SettingGetGlobal_f(G,cSetting_animation_duration));
  } else {
    sele0 = SelectorIndexByName(G,name);
    if(sele0>=0) { /* any valid selection except "all" */
      ErrMessage(G,"ExecutiveCenter","selection doesn't specify any coordinates.");
      ok=false;
    } else if(ExecutiveValidName(G,name)) {
      SceneSetDefaultView(G);
      SceneDirty(G);
    } else {
      ErrMessage(G,"ExecutiveCenter","selection or object unknown.");
      ok=false;
    }
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveOrigin(PyMOLGlobals *G,char *name,int preserve,char *oname,float *pos,int state)
{
  float center[3];
  float mn[3],mx[3];
  int ok=true;
  CObject *obj = NULL;
  if(oname[0]) {
    obj = ExecutiveFindObjectByName(G,oname);
    if(!obj)
      ok=false;
  }
  if(ok) {
    if(name[0]) {
      ok = ExecutiveGetExtent(G,name,mn,mx,(oname[0]==0),state,true);
      if(ok) 
        average3f(mn,mx,center);
    } else {
      copy3f(pos,center)
        }
  }
  if(ok) {
    if(obj) {
      ObjectSetTTTOrigin(obj,center);
      PRINTFB(G,FB_Executive,FB_Blather)
        " ExecutiveCenter: origin for %s set to %8.3f %8.3f %8.3f\n",
        oname,center[0],center[1],center[2]
        ENDFB(G);
    } else {
      PRINTFB(G,FB_Executive,FB_Blather)
        " ExecutiveCenter: scene origin set to %8.3f %8.3f %8.3f\n",
        center[0],center[1],center[2]
        ENDFB(G);
      SceneOriginSet(G,center,preserve);
    }
    SceneDirty(G);
  } else
    ok=false;
  return(ok);
}
/*========================================================================*/
int ExecutiveGetMoment(PyMOLGlobals *G,char *name,Matrix33d mi,int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int a,b;
  int c=0;

  if(state==-2) state=SceneGetState(G);

  for(a=0;a<3;a++)
	 {
		for(b=0;b<3;b++)
		  mi[a][b]=0.0;
		mi[a][a]=1.0;
	 }
  
  sele=SelectorIndexByName(G,name);
  if(sele>=0) {
    ObjectMoleculeOpRecInit(&op);
    if(state<0) {
      op.code = OMOP_SUMC;
    } else {
      op.code = OMOP_CSetSumVertices;
      op.cs1=state;
    }
    
    op.v1[0]=0.0;
    op.v1[1]=0.0;
    op.v1[2]=0.0;
    op.i1=0;
    op.i2=0; /* untransformed...is this right? */
	 
	 ExecutiveObjMolSeleOp(G,sele,&op);
	 
	 if(op.i1) { /* any vertices? */
		c+=op.i1;
		scale3f(op.v1,1.0F/op.i1,op.v1); /* compute raw average */
      if(state<0) {
        op.code = OMOP_MOME;		
      } else {
        op.code = OMOP_CSetMoment;
        op.cs1=state;
      }
		for(a=0;a<3;a++)
		  for(b=0;b<3;b++)
			 op.d[a][b]=0.0;
		ExecutiveObjMolSeleOp(G,sele,&op);			 
		for(a=0;a<3;a++)
		  for(b=0;b<3;b++)
			 mi[a][b]=op.d[a][b];
	 }
  } 
  return(c);
}
/*========================================================================*/
void ExecutiveSetObjVisib(PyMOLGlobals *G,char *name,int state)
{
  register CExecutive *I = G->Executive;
  SpecRec *tRec;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetObjVisib: entered.\n"
    ENDFD;

  if(strcmp(name,cKeywordAll)==0) {
    tRec=NULL;
    while(ListIterate(I->Spec,tRec,next)) {
      if(state!=tRec->visible) {
        if(tRec->type==cExecObject) {
          if(tRec->visible)
            SceneObjectDel(G,tRec->obj);				
          else {
            SceneObjectAdd(G,tRec->obj);
          }
        }
        if((tRec->type!=cExecSelection)||(!state)) /* hide all selections, but show all */
          tRec->visible=!tRec->visible;
      }
    }
  } else {
    tRec = ExecutiveFindSpec(G,name);
    if(tRec) {
      if(tRec->type==cExecObject) {
        if(tRec->visible!=state)
          {
            if(tRec->visible)
              SceneObjectDel(G,tRec->obj);				
            else {
              SceneObjectAdd(G,tRec->obj);
            }
            tRec->visible=!tRec->visible;
          }
      }
      else if(tRec->type==cExecSelection) {
        if(tRec->visible!=state) {
          tRec->visible=!tRec->visible;
          if(tRec->visible)
            if(SettingGetGlobal_b(G,cSetting_active_selections)) {
              ExecutiveHideSelections(G);
              tRec->visible=true;
            }
          SceneDirty(G);
          SeqDirty(G);
        }
      }
    }
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveSetObjVisib: leaving...\n"
    ENDFD;

}



/*========================================================================*/
void ExecutiveFullScreen(PyMOLGlobals *G,int flag)
{

#ifndef _PYMOL_NO_GLUT
  register CExecutive *I = G->Executive;
  if(G->HaveGUI && G->ValidContext) {
    if(!SettingGet(G,cSetting_full_screen))
      {
        I->oldPX = p_glutGet(GLUT_WINDOW_X);
        I->oldPY = p_glutGet(GLUT_WINDOW_Y);
        I->oldWidth = p_glutGet(GLUT_WINDOW_WIDTH);
        I->oldHeight = p_glutGet(GLUT_WINDOW_HEIGHT);
        I->sizeFlag = true;
      }
      
    SettingSet(G,cSetting_full_screen,(float)flag);
    if(flag) {
      p_glutFullScreen();
    } else {
      if(I->sizeFlag) {
        p_glutPositionWindow(I->oldPX,I->oldPY);
        p_glutReshapeWindow(I->oldWidth,I->oldHeight);
      } else {
#ifndef _PYMOL_NO_MAIN
        MainRepositionWindowDefault(G);
#endif
      }
    }
  }
#endif

}
/*========================================================================*/
void ExecutiveSetAllVisib(PyMOLGlobals *G,int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int rep;
  int sele;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllVisib: entered.\n"
    ENDFD;


  while(ListIterate(I->Spec,rec,next)) {
	 if(rec->type==cExecObject)
		{
        switch(rec->obj->type) {
        case cObjectMolecule:
          obj=(ObjectMolecule*)rec->obj;
          sele = SelectorIndexByName(G,obj->Obj.Name);
          for(rep=0;rep<cRepCnt;rep++) 
            rec->repOn[rep]=state;
          ObjectMoleculeOpRecInit(&op);

          op.code=OMOP_VISI;
          op.i1=-1;
          op.i2=state;
          ObjectMoleculeSeleOp(obj,sele,&op);
          op.code=OMOP_INVA;
          op.i1=-1;
          op.i2=cRepInvVisib;
          ObjectMoleculeSeleOp(obj,sele,&op);				
          break;
        default:
          for(rep=0;rep<cRepCnt;rep++) {
            ObjectSetRepVis(rec->obj,rep,state);
            if(rec->obj->fInvalidate)
              rec->obj->fInvalidate(rec->obj,rep,cRepInvVisib,state);
          }
          SceneDirty(G);
          break;
        }
		}
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
int ExecutiveToggleRepVisib(PyMOLGlobals *G,char *name,int rep)
{
  int ok =true;
  int sele;
  int handled = false;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;

  PRINTFD(G,FB_Executive)
    " ExecutiveToggleRepVisib: entered.\n"
    ENDFD;

  tRec = ExecutiveFindSpec(G,name);
  if((!tRec)&&(!strcmp(name,cKeywordAll))) {
    ExecutiveToggleAllRepVisib(G,name,rep);
  }
  if(tRec) {
    if(tRec->type==cExecObject) 
      switch(tRec->obj->type) {
      case cObjectMolecule: /* do nothing -- yet */
        break;
      default: /* otherwise, toggle the representation on/off */
        if(rep>=0) {
          ObjectToggleRepVis(tRec->obj,rep);
          if(tRec->obj->fInvalidate)
            tRec->obj->fInvalidate(tRec->obj,rep,cRepInvVisib,0);
        } 
        handled = true;
        SceneChanged(G);
        break;
      }
    if(!handled)
      switch(tRec->type) {
      case cExecSelection:
      case cExecObject:
        sele=SelectorIndexByName(G,name);
        if(sele>=0) {
          ObjectMoleculeOpRecInit(&op);

          op.code=OMOP_CheckVis;
          op.i1=rep;
          op.i2=false;
          ExecutiveObjMolSeleOp(G,sele,&op);
          op.i2 = !op.i2;

          if(tRec->type==cExecObject)
            ObjectSetRepVis(tRec->obj,rep,op.i2);

          op.code=OMOP_VISI;
          op.i1=rep;
          ExecutiveObjMolSeleOp(G,sele,&op);
          op.code=OMOP_INVA;
          op.i2=cRepInvVisib;
          ExecutiveObjMolSeleOp(G,sele,&op);
        }
        break;
      }
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveToggleRepVisib: leaving...\n"
    ENDFD;
  return (ok);
}

/*========================================================================*/
void ExecutiveSetRepVisib(PyMOLGlobals *G,char *name,int rep,int state)
{
  int sele;
  int a;
  int handled = false;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetRepVisib: entered.\n"
    ENDFD;

  tRec = ExecutiveFindSpec(G,name);
  if((!tRec)&&(!strcmp(name,cKeywordAll))) {
    ExecutiveSetAllRepVisib(G,name,rep,state);
  }
  if(tRec) {
	 if(name[0]!='_') {
      /* remember visibility information for real selections */
	   if(rep>=0) {
        tRec->repOn[rep]=state;
      } else {
        for(a=0;a<cRepCnt;a++)
          tRec->repOn[a]=state; 
      }
	 }
    if(tRec->type==cExecObject) 
      switch(tRec->obj->type) {
      default:
        if(rep>=0) {
          ObjectSetRepVis(tRec->obj,rep,state);
          if(tRec->obj->fInvalidate)
            tRec->obj->fInvalidate(tRec->obj,rep,cRepInvVisib,0);
        } else {
          for(a=0;a<cRepCnt;a++) {
            tRec->repOn[a]=state; 
            ObjectSetRepVis(tRec->obj,a,state);
            if(tRec->obj->fInvalidate)
              tRec->obj->fInvalidate(tRec->obj,a,cRepInvVisib,0);
          }
        }
        SceneChanged(G);
        break;
      }
    if(!handled)
      switch(tRec->type) {
      case cExecSelection:
      case cExecObject:
        sele=SelectorIndexByName(G,name);
        if(sele>=0) {
          ObjectMoleculeOpRecInit(&op);

          op.code=OMOP_VISI;
          op.i1=rep;
          op.i2=state;
          ExecutiveObjMolSeleOp(G,sele,&op);
          op.code=OMOP_INVA;
          op.i2=cRepInvVisib;
          ExecutiveObjMolSeleOp(G,sele,&op);
        }
        break;
      }
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveSetRepVisib: leaving...\n"
    ENDFD;

}

/*========================================================================*/
int ExecutiveSetOnOffBySele(PyMOLGlobals *G,char *name,int onoff)
{
  int sele;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;

  tRec = ExecutiveFindSpec(G,name);
  if((!tRec)&&(!strcmp(name,cKeywordAll))) {
    ExecutiveSetObjVisib(G,name,onoff);
  }
  if(tRec) {
    sele=SelectorIndexByName(G,name);
    if(sele>=0) {
      ObjectMoleculeOpRecInit(&op);
      
      op.code=OMOP_OnOff;
      op.i1=onoff;
      ExecutiveObjMolSeleOp(G,sele,&op);
    }
  }
  return 1;
}


/*========================================================================*/
void ExecutiveSetAllRepVisib(PyMOLGlobals *G,char *name,int rep,int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int sele;
  int a;

  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllRepVisib: entered.\n"
    ENDFD;
  while(ListIterate(I->Spec,rec,next)) {
	 if(rec->type==cExecObject)
		{
        if(rec->name[0]!='_') {
          /* remember visibility information for real selections */
          if(rep>=0) {
            rec->repOn[rep]=state;
          } else {
            for(a=0;a<cRepCnt;a++)
              rec->repOn[a]=state; 
          }
        }   
        if(rec->type==cExecObject) {
          switch(rec->obj->type) {
          case cObjectMolecule:
            if(rep>=0) {
              rec->repOn[rep]=state;
            } else {
              for(a=0;a<cRepCnt;a++)
                rec->repOn[a]=state;
            }
            obj=(ObjectMolecule*)rec->obj;
            sele = SelectorIndexByName(G,obj->Obj.Name);
            ObjectMoleculeOpRecInit(&op);

            op.code=OMOP_VISI;
            op.i1=rep;
            op.i2=state;
            ObjectMoleculeSeleOp(obj,sele,&op);
            op.code=OMOP_INVA;
            op.i2=cRepInvVisib;
            ObjectMoleculeSeleOp(obj,sele,&op);				
            break;
          default:
            if(rep>=0) {
              ObjectSetRepVis(rec->obj,rep,state);
              if(rec->obj->fInvalidate)
                rec->obj->fInvalidate(rec->obj,rep,cRepInvVisib,state);
            } else {
              for(a=0;a<cRepCnt;a++) {
                ObjectSetRepVis(rec->obj,a,state);
                if(rec->obj->fInvalidate)
                  rec->obj->fInvalidate(rec->obj,rep,cRepInvVisib,state);
              }
            }
            SceneDirty(G);
            break;
          }
        }
		}
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllRepVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void ExecutiveToggleAllRepVisib(PyMOLGlobals *G,char *name,int rep)
{
  ObjectMoleculeOpRec op;
  int toggle_state;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  op.code=OMOP_CheckVis;
  op.i1=rep;
  op.i2=false;
  ExecutiveObjMolSeleOp(G,cSelectionAll,&op);
  toggle_state = op.i2;
  while(ListIterate(I->Spec,rec,next)) {
	 if(rec->type==cExecObject) {
      switch(rec->obj->type) {
        case cObjectMolecule:
        break;
        default:
        if(rec->repOn[rep])
          toggle_state = true;
        break;
      }
    }
  }

  ExecutiveSetAllRepVisib(G,name,rep,!toggle_state);
}
/*========================================================================*/
void ExecutiveInvalidateRep(PyMOLGlobals *G,char *name,int rep,int level)
{
  register CExecutive *I = G->Executive;
  int sele = -1;
  ObjectMoleculeOpRec op;
  int all_flag=false;
  SpecRec *rec = NULL;
  PRINTFD(G,FB_Executive)
    "ExecInvRep-Debug: %s %d %d\n",name,rep,level
    ENDFD;
  if(WordMatch(G,cKeywordAll,name,true)<0) {
    all_flag=true;
  }
  if(all_flag) {
    while(ListIterate(I->Spec,rec,next))
      if(rec->type==cExecObject) {
        if(rec->obj->fInvalidate) {
          rec->obj->fInvalidate(rec->obj,rep,cRepInvColor,cRepAll);
          SceneDirty(G);
        }
      }
  }
  sele=SelectorIndexByName(G,name);
  if(sele>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_INVA;
    op.i1=rep;
    op.i2=level;
    ExecutiveObjMolSeleOp(G,sele,&op);
  }
}


/*========================================================================*/
CObject *ExecutiveFindObjectByName(PyMOLGlobals *G,char *name)
{
  /* TODO: switch over to using a Python Dictionary to store objects 
     instead of this stupid list */

  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  CObject *obj=NULL;
  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->type==cExecObject)
		  {
			 if(strcmp(rec->obj->Name,name)==0) 
				{
				  obj=rec->obj;
				  break;
				}
		  }
	 }
  return(obj);
}
/*========================================================================*/
ObjectMap *ExecutiveFindObjectMapByName(PyMOLGlobals *G,char *name)
{
  CObject *obj;
  
  obj = ExecutiveFindObjectByName(G,name);
  if(obj)
    if(obj->type!=cObjectMap)
      obj=NULL;
  return((ObjectMap*)obj);
}

/*========================================================================*/
ObjectMolecule *ExecutiveFindObjectMoleculeByName(PyMOLGlobals *G,char *name)
{
  CObject *obj;
  
  obj = ExecutiveFindObjectByName(G,name);
  if(obj)
    if(obj->type!=cObjectMolecule)
      obj=NULL;
  return((ObjectMolecule*)obj);
}
/*========================================================================*/
Block *ExecutiveGetBlock(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  return(I->Block);
}
/*========================================================================*/
void ExecutiveSetControlsOff(PyMOLGlobals *G,char *name)
{
  SpecRec *rec;
  int a;
  rec = ExecutiveFindSpec(G,name);
  if(rec)
	 {
		for(a=0;a<cRepCnt;a++)
		  rec->repOn[a]=false;
	 }
}
/*========================================================================*/
void ExecutiveSymExp(PyMOLGlobals *G,char *name,char *oname,char *s1,float cutoff) /* TODO state */
{
  CObject *ob;
  ObjectMolecule *obj = NULL;
  ObjectMolecule *new_obj = NULL;
  ObjectMoleculeOpRec op;
  MapType *map;
  int x,y,z,a,b,c,i,j,h,k,l,n;
  CoordSet *cs,*os;
  int keepFlag,sele,tt[3];
  float *v1,*v2,m[16],tc[3],ts[3];
  OrthoLineType new_name;
  float auto_save;

  PRINTFD(G,FB_Executive)
    " ExecutiveSymExp: entered.\n"
    ENDFD;

  auto_save = SettingGet(G,cSetting_auto_zoom);
  SettingSet(G,cSetting_auto_zoom,0);
  sele=SelectorIndexByName(G,s1);
  ob = ExecutiveFindObjectByName(G,oname);
  if(ob->type==cObjectMolecule)
    obj=(ObjectMolecule*)ob;
  if(!(obj&&sele)) {
    ErrMessage(G,"ExecutiveSymExp","Invalid object");
  } else if(!obj->Symmetry) {
    ErrMessage(G,"ExecutiveSymExp","No symmetry loaded!");
  } else if(!obj->Symmetry->NSymMat) {
    ErrMessage(G,"ExecutiveSymExp","No symmetry matrices!");    
  } else {
    PRINTFB(G,FB_Executive,FB_Actions)
      " ExecutiveSymExp: Generating symmetry mates...\n"
      ENDFB(G);
    ObjectMoleculeOpRecInit(&op);
	 op.code = OMOP_SUMC;
	 op.i1 =0;
    op.i2 =0;
    op.v1[0]= 0.0;
    op.v1[1]= 0.0;
    op.v1[2]= 0.0;
    ExecutiveObjMolSeleOp(G,sele,&op);
    tc[0]=op.v1[0];
    tc[1]=op.v1[1];
    tc[2]=op.v1[2];
    if(op.i1) {
      tc[0]/=op.i1;
      tc[1]/=op.i1;
      tc[2]/=op.i1;
    }
    transform33f3f(obj->Symmetry->Crystal->RealToFrac,tc,tc);

	 op.code = OMOP_VERT;
	 op.nvv1 =0;
    op.vv1 = VLAlloc(float,10000);
    ExecutiveObjMolSeleOp(G,sele,&op);
    
    if(!op.nvv1) {
      ErrMessage(G,"ExecutiveSymExp","No atoms indicated!");          
    } else {
      map=MapNew(G,-cutoff,op.vv1,op.nvv1,NULL);
      if(map) {
        MapSetupExpress(map);  
        /* go out no more than one lattice step in each direction */
        for(x=-1;x<2;x++)
          for(y=-1;y<2;y++)
            for(z=-1;z<2;z++)
              for(a=0;a<obj->Symmetry->NSymMat;a++) {
                new_obj = ObjectMoleculeCopy(obj);
                keepFlag=false;
                for(b=0;b<new_obj->NCSet;b++) 
                  if(new_obj->CSet[b]) {
                    cs = new_obj->CSet[b];
                    os = obj->CSet[b];
                    CoordSetRealToFrac(cs,obj->Symmetry->Crystal);
                    CoordSetTransform44f(cs,obj->Symmetry->SymMatVLA+(a*16));
                    CoordSetGetAverage(cs,ts);
                    identity44f(m);
                    /* compute the effective translation resulting
                       from application of the symmetry operator so
                       that we can shift it into the cell of the
                       target selection */
                    for(c=0;c<3;c++) { /* manual rounding - rint broken */
                      ts[c]=tc[c]-ts[c];
                      if(ts[c]<0)
                        ts[c]-=0.5;
                      else
                        ts[c]+=0.5;
                      tt[c]=(int)ts[c];
                    }
                    m[3] = (float)tt[0]+x;
                    m[7] = (float)tt[1]+y;
                    m[11] = (float)tt[2]+z;
                    CoordSetTransform44f(cs,m);
                    CoordSetFracToReal(cs,obj->Symmetry->Crystal);
                    if(!keepFlag) {
                      v2 = cs->Coord;
                      n=cs->NIndex;
                      while(n--) {
                        MapLocus(map,v2,&h,&k,&l);
                        i=*(MapEStart(map,h,k,l));
                        if(i) {
                          j=map->EList[i++];
                          while(j>=0) {
                            if(within3f(op.vv1+3*j,v2,cutoff)) {
                              keepFlag=true;
                              break;
                            }
                            j=map->EList[i++];
                          }
                        }
                        v2+=3;
                        if(keepFlag) break;
                      }
                    }
                    if(keepFlag) { /* make sure that we don't aren't simply duplicating the template coordinates */
                      keepFlag = false;
                      v1 = os->Coord;
                      v2 = cs->Coord;
                      n = cs->NIndex;
                      while(n--) {
                        if(diffsq3f(v1,v2)>R_SMALL8) {
                          keepFlag = true;
                          break;
                        }
                        v1++;
                        v2++;
                      }
                    }
                  }
                if(keepFlag) { /* we need to create new object */
                  sprintf(new_name,"%s%02d%02d%02d%02d",name,a,x,y,z);
                  ObjectSetName((CObject*)new_obj,new_name);
                  ExecutiveDelete(G,new_name);
                  ExecutiveManageObject(G,(CObject*)new_obj,true,false);
                  SceneChanged(G);
                } else {
                  ((CObject*)new_obj)->fFree((CObject*)new_obj);
                }
              }
        MapFree(map);
      }
    }
    VLAFreeP(op.vv1);
  }
  PRINTFD(G,FB_Executive)
     " ExecutiveSymExp: leaving...\n"
    ENDFD;
  SettingSet(G,cSetting_auto_zoom,auto_save);
}
/*========================================================================*/
void ExecutiveDelete(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int all_flag=false;
  WordType name_copy; /* needed in case the passed string changes */

  if(WordMatch(G,name,cKeywordAll,true)<0) all_flag=true;
  strcpy(name_copy,name);
  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->type==cExecObject)
		  {
          if(I->LastEdited==cExecObject) 
            I->LastEdited=NULL;
			 if(all_flag||(WordMatch(G,name_copy,rec->obj->Name,true)<0))
				{
              if(rec->obj->type == cObjectMolecule)
                if(EditorIsAnActiveObject(G,(ObjectMolecule*)rec->obj))
                  EditorInactivate(G);
              SeqChanged(G);
              if(rec->visible) 
                SceneObjectDel(G,rec->obj);
				  SelectorDelete(G,rec->name);
				  rec->obj->fFree(rec->obj);
				  rec->obj=NULL;
				  ListDelete(I->Spec,rec,next,SpecRec);
				  rec=NULL;
				}
		  }
		else if(rec->type==cExecSelection)
		  {

			 if(all_flag||(WordMatch(G,name_copy,rec->name,true)<0))
				{
              if(all_flag||rec->visible)
                SceneDirty(G);
              SeqDirty(G);
				  SelectorDelete(G,rec->name);
				  ListDelete(I->Spec,rec,next,SpecRec);
				  rec=NULL;
				}
		  }
	 }
  if(all_flag)
    SelectorDefragment(G);
}
/*========================================================================*/
void ExecutiveDump(PyMOLGlobals *G,char *fname,char *obj)
{
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;

  SceneUpdate(G);

  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->type==cExecObject)
		  {
			 if(strcmp(rec->obj->Name,obj)==0) 
				break;
		  }
	 }
  if(rec)
	 { 
      if(rec->obj->type==cObjectMesh) {
        ObjectMeshDump((ObjectMesh*)rec->obj,fname,0);
      } else if(rec->obj->type==cObjectSurface) {
        ObjectSurfaceDump((ObjectSurface*)rec->obj,fname,0);
      } else {
        ErrMessage(G,"ExecutiveDump","Invalid object type for this operation.");
      }
	 }
  else {
    ErrMessage(G,"ExecutiveDump","Object not found.");
  }
  
}
/*========================================================================*/
void ExecutiveManageObject(PyMOLGlobals *G,CObject *obj,int allow_zoom,int quiet)
{
  int a;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;
  int exists=false;

  if(SettingGet(G,cSetting_auto_hide_selections))
    ExecutiveHideSelections(G);
  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->obj==obj) {
        exists = true;
      }
	 }
  if(!exists) {
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject)
          {
            if(strcmp(rec->obj->Name,obj->Name)==0) 
              break;
          }
      }
    if(rec) /* another object of this type already exists */
      { /* purge it */
        SceneObjectDel(G,rec->obj);
        rec->obj->fFree(rec->obj);
        rec->obj=NULL;
      }
    else 
      {
        if(!quiet)
          if(obj->Name[0]!='_') { /* suppress internal objects */
            PRINTFB(G,FB_Executive,FB_Actions)
              " Executive: object \"%s\" created.\n",obj->Name 
              ENDFB(G);
          }
      }
    if(!rec)
      ListElemCalloc(G,rec,SpecRec);

    if(WordMatch(G,cKeywordAll,obj->Name,true)<0) {
      PRINTFB(G,FB_Executive,FB_Warnings) 
        " Executive: object name \"%s\" is illegal -- renamed to 'all_'.",obj->Name
        ENDFB(G);
      strcat(obj->Name,"_"); /* don't allow object named "all" */
    }
    strcpy(rec->name,obj->Name);
    rec->type=cExecObject;
    rec->next=NULL;
    rec->obj=obj;
    if(rec->obj->type==cObjectMap) {
      rec->visible=0;
    } else {
      rec->visible=1;
      SceneObjectAdd(G,obj);
    }
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    if(rec->obj->type==cObjectMolecule)
      rec->repOn[cRepLine]=true;
    ListAppend(I->Spec,rec,next,SpecRec);
  }
  if(obj->type==cObjectMolecule) {
	 ExecutiveUpdateObjectSelection(G,obj);
  }

  if(SettingGet(G,cSetting_auto_dss)) {
    if(obj->type==cObjectMolecule) {
      ObjectMolecule *objMol = (ObjectMolecule*)obj;
      if(objMol->NCSet==1) {
        ExecutiveAssignSS(G,obj->Name,0,"",1,1);
      }
    }
  }

  if(allow_zoom)
    if(!exists) {
      switch(SettingGetGlobal_i(G,cSetting_auto_zoom)) {
      case 1: /* zoom new one */
        ExecutiveWindowZoom(G,obj->Name,0.0,-1,0,0); /* auto zoom (all states) */
        break;
      case 2: /* zoom all */
        ExecutiveWindowZoom(G,cKeywordAll,0.0,-1,0,0);
        break;
      }
    }
  SeqChanged(G);
}
/*========================================================================*/
void ExecutiveManageSelection(PyMOLGlobals *G,char *name)
{

  int a;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;
  int hide_all  = SettingGetGlobal_b(G,cSetting_active_selections);
  if(name[0]=='_')
    hide_all = false; /* hidden selections don't change active selection */
  while(ListIterate(I->Spec,rec,next))
    {
      if(rec->type==cExecSelection) {
        if(strcmp(rec->name,name)==0) 
          break;
        if(hide_all)
          rec->visible=false;
      }
    }
  if(rec&&hide_all)
    while(ListIterate(I->Spec,rec,next))
      if(rec->type==cExecSelection)
        rec->visible=false;

  if(!rec) {
    ListElemCalloc(G,rec,SpecRec);
    strcpy(rec->name,name);
    rec->type=cExecSelection;
    rec->next=NULL;
    rec->sele_color=-1;
    rec->visible=false;
    ListAppend(I->Spec,rec,next,SpecRec);
  }
  if(rec) {
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    if(name[0]!='_') {
      if(SettingGet(G,cSetting_auto_hide_selections))
        ExecutiveHideSelections(G);
      if(SettingGet(G,cSetting_auto_show_selections)) {
        rec->visible=true;
      }
    }
    if(rec->visible) SceneDirty(G);
  }
  SeqDirty(G);
}
/*========================================================================*/
static int ExecutiveClick(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;
  int n,a;
  SpecRec *rec = NULL;
  int t;
  int pass = false;
  int skip;
  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);

  if(y<I->HowFarDown) {
    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==1) 
      return SceneDeferClick(SceneGetBlock(G),button,x,y,mod);
  }
  n=((I->Block->rect.top-y)-(ExecTopMargin+ExecClickMargin))/ExecLineHeight;
  a=n;
  if(I->ScrollBarActive) {
    if((x-I->Block->rect.left)<(ExecScrollBarWidth+ExecScrollBarMargin+ExecToggleMargin)) {
      pass = 1;
      ScrollBarDoClick(I->ScrollBar,button,x,y,mod);      
    }
  } 
  skip = I->NSkip;
  if(!pass)   
    while(ListIterate(I->Spec,rec,next))
      if(rec->name[0]!='_')
        {
          if(skip) {
            skip--;
          } else {
            if(!a) {
              t = ((I->Block->rect.right-ExecRightMargin)-x)/ExecToggleWidth;
              if(t<ExecOpCnt) {
                int my = I->Block->rect.top-(ExecTopMargin + n*ExecLineHeight)-3;
                int mx = I->Block->rect.right-(ExecRightMargin + t*ExecToggleWidth);
                t = (ExecOpCnt-t)-1;
                switch(t) {
                case 0:
                  switch(rec->type) {
                  case cExecAll:
                    MenuActivate(G,mx,my,x,y,false,"all_action",rec->name);
                    break;
                  case cExecSelection:
                    MenuActivate(G,mx,my,x,y,false,"sele_action",rec->name);
                    break;
                  case cExecObject:
                    switch(rec->obj->type) {
                    case cObjectMolecule:
                      MenuActivate(G,mx,my,x,y,false,"mol_action",rec->obj->Name);
                      break;
                    case cObjectSurface:
                    case cObjectMesh:
                    case cObjectDist:
                    case cObjectMap:
                    case cObjectCGO:
                    case cObjectCallback:
                      MenuActivate(G,mx,my,x,y,false,"simple_action",rec->obj->Name);
                      break;
                    case cObjectSlice:
                      MenuActivate(G,mx,my,x,y,false,"slice_action",rec->obj->Name);
                      break;
                    case cObjectGadget:
                      MenuActivate(G,mx,my,x,y,false,"ramp_action",rec->obj->Name);
                      break;
                    }
                    break;
                  }
                  break;
                case 1:
                  switch(rec->type) {
                  case cExecAll:
                    MenuActivate(G,mx,my,x,y,false,"mol_show",cKeywordAll);
                    break;
                  case cExecSelection:
                    MenuActivate(G,mx,my,x,y,false,"mol_show",rec->name);
                    break;
                  case cExecObject:
                    switch(rec->obj->type) {
                    case cObjectMolecule:
                      MenuActivate(G,mx,my,x,y,false,"mol_show",rec->obj->Name);
                      break;
                    case cObjectCGO:
                      MenuActivate(G,mx,my,x,y,false,"cgo_show",rec->obj->Name);
                      break;
                    case cObjectDist:
                      MenuActivate(G,mx,my,x,y,false,"dist_show",rec->obj->Name);
                      break;
                    case cObjectMap:
                      MenuActivate(G,mx,my,x,y,false,"simple_show",rec->obj->Name);
                      break;
                    case cObjectSurface:
                    case cObjectMesh:
                      MenuActivate(G,mx,my,x,y,false,"mesh_show",rec->obj->Name);
                      break;
                    case cObjectSlice:
                      MenuActivate(G,mx,my,x,y,false,"slice_show",rec->obj->Name);
                      break;

                    }
                    break;
                  }
                  break;
                case 2:
                  switch(rec->type) {
                  case cExecAll:
                    MenuActivate(G,mx,my,x,y,false,"mol_hide",cKeywordAll);
                    break;
                  case cExecSelection:
                    MenuActivate(G,mx,my,x,y,false,"mol_hide",rec->name);
                    break;
                  case cExecObject:
                    switch(rec->obj->type) {
                    case cObjectMolecule:
                      MenuActivate(G,mx,my,x,y,false,"mol_hide",rec->obj->Name);
                      break;
                    case cObjectCGO:
                      MenuActivate(G,mx,my,x,y,false,"cgo_hide",rec->obj->Name);
                      break;
                    case cObjectDist:
                      MenuActivate(G,mx,my,x,y,false,"dist_hide",rec->obj->Name);
                      break;
                    case cObjectMap:
                      MenuActivate(G,mx,my,x,y,false,"simple_hide",rec->obj->Name);
                      break;
                    case cObjectSurface:
                    case cObjectMesh:
                      MenuActivate(G,mx,my,x,y,false,"mesh_hide",rec->obj->Name);
                      break;
                    case cObjectSlice:
                      MenuActivate(G,mx,my,x,y,false,"slice_hide",rec->obj->Name);
                      break;

                    }
                    break;
                  }
                  break;
                case 3:
                  switch(rec->type) {
                  case cExecAll:
                    MenuActivate(G,mx,my,x,y,false,"mol_labels","(all)");
                    break;
                  case cExecSelection:
                    MenuActivate(G,mx,my,x,y,false,"mol_labels",rec->name);
                    break;
                  case cExecObject:
                    switch(rec->obj->type) {
                    case cObjectMolecule:
                      MenuActivate(G,mx,my,x,y,false,"mol_labels",rec->obj->Name);
                      break;
                    case cObjectDist:
                      break;
                    case cObjectMap:
                    case cObjectSurface:
                    case cObjectMesh:
                    case cObjectSlice:
                      break;
                    }
                    break;
                  }
                  break;
                case 4:
                  switch(rec->type) {
                  case cExecAll:
                  case cExecSelection:
                    MenuActivate(G,mx,my,x,y,false,"mol_color",rec->name);
                    break;
                  case cExecObject:
                    switch(rec->obj->type) {
                    case cObjectMolecule:
                      MenuActivate(G,mx,my,x,y,false,"mol_color",rec->obj->Name);
                      break;
                    case cObjectDist:
                    case cObjectMap:
                    case cObjectSurface:
                    case cObjectCGO:
                    case cObjectMesh:
                      MenuActivate(G,mx,my,x,y,false,"general_color",rec->obj->Name);
                      break;
                    case cObjectSlice:
                      MenuActivate(G,mx,my,x,y,false,"slice_color",rec->obj->Name);
                      break;
                    }
                    break;
                  }
                  break;
                }
              } else {
                rec->hilight=true;
                switch(button) {
                case P_GLUT_LEFT_BUTTON:
                  
                  I->Pressed = n;
                  I->OldVisibility = rec->visible;
                  I->Over = n;
                  I->DragMode = 0;
                  I->ToggleMode=0;
                  I->LastChanged=NULL;
                  
                  if(mod&cOrthoSHIFT) {
                    ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod);
                    I->ToggleMode=1;
                  }
                  if(mod&cOrthoCTRL) {
                    I->ToggleMode=2;
                    if(!rec->visible) {
                      ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod);
                    }
                    I->LastChanged=rec;
                  } 
                  break;
                case P_GLUT_MIDDLE_BUTTON:
                  I->DragMode = 1; /* reorder */
                  I->Pressed = n;
                  I->Over = n;
                  
                  break;
                }
                OrthoGrab(G,I->Block);
                OrthoDirty(G);
              }
            }
            a--;
          }
        }
  PyMOL_NeedRedisplay(G->PyMOL);
  return(1);
}
/*========================================================================*/
static void ExecutiveSpecSetVisibility(PyMOLGlobals *G,SpecRec *rec,
                                      int new_vis,int mod)
{
  OrthoLineType buffer = "";

  if(rec->type==cExecObject)
    {
      if(rec->visible&&!new_vis) {
        if(SettingGet(G,cSetting_logging)) 
          sprintf(buffer,"cmd.disable('%s')",rec->obj->Name);
        SceneObjectDel(G,rec->obj);			
      }
      else if((!rec->visible)&&new_vis) {
        sprintf(buffer,"cmd.enable('%s')",rec->obj->Name);
        SceneObjectAdd(G,rec->obj);
      }
      SceneChanged(G);
      if(SettingGet(G,cSetting_logging)) {
        PLog(buffer,cPLog_pym);
      }
      rec->visible=new_vis;
    }
  else if(rec->type==cExecAll)
    {
      if(SettingGet(G,cSetting_logging)) {
        if(rec->visible)
          sprintf(buffer,"cmd.disable('all')");
        else
          sprintf(buffer,"cmd.enable('all')");
        PLog(buffer,cPLog_pym);
      }
      ExecutiveSetObjVisib(G,cKeywordAll,!rec->visible);
    }
  else if(rec->type==cExecSelection)
    {
      if(mod&cOrthoCTRL) {
        SettingSet(G,cSetting_selection_overlay,
                   (float)(!((int)SettingGet(G,cSetting_selection_overlay))));
        if(SettingGet(G,cSetting_logging)) {
          sprintf(buffer,"cmd.set('selection_overlay',%d)",
                  (int)SettingGet(G,cSetting_selection_overlay));
          PLog(buffer,cPLog_pym);
          sprintf(buffer,"cmd.enable('%s')",rec->name);
          PLog(buffer,cPLog_pym);
        }
        rec->visible=true; 
      } else if(mod&cOrthoSHIFT) {
        if(rec->sele_color<7)
          rec->sele_color=15;
        else {
          rec->sele_color--;
          if(rec->sele_color<7)
            rec->sele_color=15;
        }
        /* NO COMMAND EQUIVALENT FOR THIS FUNCTION YET */
        rec->visible=true;
      } else {
        
        if(rec->visible&&!new_vis) {
          if(SettingGet(G,cSetting_logging)) 
            sprintf(buffer,"cmd.disable('%s')",rec->name);
        }
        else if((!rec->visible)&&new_vis) {
          sprintf(buffer,"cmd.enable('%s')",rec->name);
        }
        if(new_vis && SettingGetGlobal_b(G,cSetting_active_selections)) {
          ExecutiveHideSelections(G);
        }
        if(SettingGet(G,cSetting_logging)) {
          PLog(buffer,cPLog_pym);
        }
        rec->visible=new_vis;
      }
      SceneChanged(G);
    }
}

static int ExecutiveRelease(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;
  int n;  
  SpecRec *rec = NULL;
  int pass = false;
  int skip;

  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);

  if(y<I->HowFarDown) {
    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==1) 
      return SceneDeferRelease(SceneGetBlock(G),button,x,y,mod);
  }

  n=((I->Block->rect.top-y)-(ExecTopMargin+ExecClickMargin))/ExecLineHeight;

  if(I->ScrollBarActive) {
    if((x-I->Block->rect.left)<(ExecScrollBarWidth+ExecScrollBarMargin+ExecToggleMargin)) {
      pass = 1;
      ScrollBarDoRelease(I->ScrollBar,button,x,y,mod);
      OrthoUngrab(G);
    }
  } 

  skip=I->NSkip;

  if(!pass)
    {
      ExecutiveDrag(block,x,y,mod); /* incorporate final changes in cursor position */
      switch(I->DragMode) {
      case 0:
        
        while(ListIterate(I->Spec,rec,next)) {
          if(rec->name[0]!='_')
            {
              if(skip) {
                skip--;
              } else if(rec->hilight) {
                if(rec->type==cExecSelection) {
                  ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,0);                    
                } else {
                  ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod);
                }
              }
            }
        }
        break;
      case 1:
        if(I->ReorderFlag) {
          I->ReorderFlag=false;
          PLog(I->ReorderLog,cPLog_no_flush);
        }
        break;
      }
    }
  
  {
    SpecRec *rec=NULL;
    while(ListIterate(I->Spec,rec,next)) {
      rec->hilight=false;
    }
  }

  I->Over = -1;
  I->Pressed = -1;
  OrthoUngrab(G);
  PyMOL_NeedRedisplay(G->PyMOL);
  return(1);
}
/*========================================================================*/
static int ExecutiveDrag(Block *block,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;
  int xx,t;
  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);

  if(y<I->HowFarDown) {
    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==1) 
      return SceneDeferDrag(SceneGetBlock(G),x,y,mod);
  }

  xx = (x-I->Block->rect.left);
  t = ((I->Block->rect.right-ExecRightMargin)-x)/ExecToggleWidth;
  if(I->ScrollBarActive) {
    xx -= (ExecScrollBarWidth+ExecScrollBarMargin);
  }
  
  {
    int row_offset;
    if((xx>=0)&&(t>=ExecOpCnt)) {
      row_offset = ((I->Block->rect.top-y)-
                    (ExecTopMargin+ExecClickMargin))/ExecLineHeight;
      I->Over = row_offset;
    } else {
      I->Over = -1;
      row_offset = -1;
      {
        SpecRec *rec=NULL;
        while(ListIterate(I->Spec,rec,next))      
          rec->hilight=false;
      }
    }
    
    if(I->Over>=0) {
      SpecRec *rec = NULL;
      int skip=I->NSkip;
      int row=0;
      switch(I->DragMode) {
      case 0:
        
        while(ListIterate(I->Spec,rec,next)) {
          if(rec->name[0]!='_')
            {
              if(skip) {
                skip--;
              } else {
                rec->hilight=false;
                if( ((row>=I->Over)&&(row<=I->Pressed))||
                    ((row>=I->Pressed)&&(row<=I->Over))) {
                  switch(I->ToggleMode) {
                  case 0:
                    if(row||(row==I->Pressed))
                      rec->hilight=true;
                    break;
                  case 1:
                    if(row)
                      ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod);
                    break;
                  case 2:
                    if((row==I->Over)&&row) {
                      if(I->LastChanged!=rec)
                        ExecutiveSpecSetVisibility(G,I->LastChanged,false,mod);
                      if(!rec->visible) {
                        ExecutiveSpecSetVisibility(G,rec,true,mod);
                        I->LastChanged=rec;
                      }
                    }
                    break;
                  }
                }
                row++;
              }
            }
        }
        break;
      case 1:
        {
          int loop_flag = (I->Over!=I->Pressed);
          while(loop_flag) {
            loop_flag=false;
            if(I->Over>I->Pressed)
              I->Over = I->Pressed+1;
            else if(I->Over<I->Pressed)
              I->Over = I->Pressed-1;
            if(I->Over!=I->Pressed) {
              SpecRec *last=NULL,*new_parent=NULL,*old_parent=NULL;
              
              while(ListIterate(I->Spec,rec,next)) {
                if(rec->name[0]!='_')
                  {
                  if(skip) {
                    skip--;
                  } else {
                    if(row==I->Pressed)
                      old_parent = last;
                    if(row==I->Over)
                      new_parent = last;
                    row++;
                    last=rec;
                  }
                }
            }
            if(new_parent&&old_parent&&(new_parent!=old_parent)) {
              SpecRec *moving = old_parent->next;
              old_parent->next = moving->next;
              if(moving!=new_parent) {
                moving->next = new_parent->next;
                new_parent->next = moving;
                if(new_parent==I->Spec) {
                  sprintf(I->ReorderLog,"cmd.order(\"%s\",location='top')\n",
                         moving->name);
                  I->ReorderFlag=true;
                } else if(new_parent&&moving) {
                  sprintf(I->ReorderLog,"cmd.order(\"%s %s\")\n",
                          new_parent->name,moving->name);
                  I->ReorderFlag=true;
                }
              } else {
                old_parent->next = moving->next;
                moving->next = old_parent->next->next;
                old_parent->next->next = moving;
                if(old_parent->next&&moving) {
                  sprintf(I->ReorderLog,"cmd.order(\"%s %s\")\n",
                          old_parent->next->name,moving->name);
                  I->ReorderFlag=true;
                }
              }
              if(I->Pressed!=I->Over)
                loop_flag=true;
              I->Pressed=I->Over;
            }
            }
          }
        }
        break;
      }
    } else if(I->LastChanged)
      ExecutiveSpecSetVisibility(G,I->LastChanged,false,mod);      
    OrthoDirty(G);
  }
  return(1);
}

static void draw_button(int x2,int y2, int w, int h, float *light, float *dark, float *inside)
{
  glColor3fv(light);
  glBegin(GL_POLYGON);
  glVertex2i(x2,y2);
  glVertex2i(x2,y2+h);
  glVertex2i(x2+w,y2+h);
  glVertex2i(x2+w,y2);
  glEnd();
  
  glColor3fv(dark);
  glBegin(GL_POLYGON);
  glVertex2i(x2+1,y2);
  glVertex2i(x2+1,y2+h-1);
  glVertex2i(x2+w,y2+h-1);
  glVertex2i(x2+w,y2);
  glEnd();
  
  if(inside) {
    glColor3fv(inside);
    glBegin(GL_POLYGON);
    glVertex2i(x2+1,y2+1);
    glVertex2i(x2+1,y2+h-1);
    glVertex2i(x2+w-1,y2+h-1);
    glVertex2i(x2+w-1,y2+1);
    glEnd();
  } else { /* rainbow */
    glBegin(GL_POLYGON);
    glColor3f(1.0F,0.1F,0.1F);
    glVertex2i(x2+1,y2+1);
    glColor3f(0.1F,1.0F,0.1F);
    glVertex2i(x2+1,y2+h-1);
    glColor3f(1.0F,1.0F,0.1F);
    glVertex2i(x2+w-1,y2+h-1);
    glColor3f(0.1F,0.1F,1.0F);
    glVertex2i(x2+w-1,y2+1);
    glEnd();
  }

}

#ifndef _PYMOL_NOPY
static void draw_button_char(PyMOLGlobals *G,int x2,int y2,char ch)
{
  TextSetColor3f(G,0.0F,0.0F,0.0F);
  TextSetPos2i(G,x2+ExecToggleTextShift,y2);
  TextDrawChar(G,ch);
}
#endif

/*========================================================================*/
static void ExecutiveDraw(Block *block)
{
  PyMOLGlobals *G=block->G;
  int x,y,xx,x2,y2;
  char *c=NULL;
  float enabledColor[3] = { 0.5F, 0.5F, 0.5F };
  float pressedColor[3] = { 0.7F, 0.7F, 0.7F };
  float disabledColor[3] = { 0.3F, 0.3F, 0.3F };
  float lightEdge[3] = {0.6F, 0.6F, 0.6F };
  float darkEdge[3] = {0.35F, 0.35F, 0.35F };
  float captionColor[3] = {0.2F, 0.8F, 0.2F };
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;
  int n_ent;
  int n_disp;
  int skip=0;
  int row = -1;
  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
  int text_lift = (ExecLineHeight/2)-5;

  if(G->HaveGUI && G->ValidContext) {
    int max_char;
    int nChar;
    /* do we have enough structures to warrant a scroll bar? */
    n_ent = 0;
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->name[0]!='_') 
        n_ent++;
    }

    n_disp = ((I->Block->rect.top-I->Block->rect.bottom)-(ExecTopMargin))/ExecLineHeight;
    if(n_disp<1) n_disp=1;
      
    if(n_ent>n_disp) {
      int bar_maxed = ScrollBarIsMaxed(I->ScrollBar);
      if(!I->ScrollBarActive) {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        if(bar_maxed) {
          ScrollBarMaxOut(I->ScrollBar);
          I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
        } else {
          ScrollBarSetValue(I->ScrollBar,0);
          I->NSkip =0;
        }
      } else {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        if(bar_maxed)
          ScrollBarMaxOut(I->ScrollBar);
        I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
      }
      I->ScrollBarActive = 1;

    } else {
      I->ScrollBarActive = 0;
      I->NSkip =0;
    }

    max_char = (((I->Block->rect.right-I->Block->rect.left)-(ExecLeftMargin+ExecRightMargin+4)) -
                     (ExecOpCnt*ExecToggleWidth));
    if(I->ScrollBarActive) {
      max_char -= (ExecScrollBarMargin+ExecScrollBarWidth);
    }      
    max_char/=8;

    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==0) {
      glColor3fv(I->Block->BackColor);
      BlockFill(I->Block);
    }

    if(I->ScrollBarActive) {
      ScrollBarSetBox(I->ScrollBar,I->Block->rect.top-ExecScrollBarMargin,
                      I->Block->rect.left+ExecScrollBarMargin,
                      I->Block->rect.bottom+2,
                      I->Block->rect.left+ExecScrollBarMargin+ExecScrollBarWidth);
      ScrollBarDoDraw(I->ScrollBar);
    }
    
    x = I->Block->rect.left+ExecLeftMargin;
    y = (I->Block->rect.top-ExecLineHeight)-ExecTopMargin;
    /*    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(cRepCnt+ExecOpCnt);*/
#ifndef _PYMOL_NOPY
    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(ExecOpCnt);
#else
    xx = I->Block->rect.right-ExecRightMargin;
#endif

    if(I->ScrollBarActive) {
      x+=ExecScrollBarWidth+ExecScrollBarMargin;
    }
    skip=I->NSkip;
    while(ListIterate(I->Spec,rec,next))
      if(rec->name[0]!='_')
        {
          if(skip) {
            skip--;
          } else {
            row++;
            x2=xx;
            y2=y;
            nChar = max_char;

            if((x-ExecToggleMargin)-(xx-ExecToggleMargin)>-10) {
              x2 = x+10;
            }
#ifndef _PYMOL_NOPY
            {
              int a;
              float toggleColor[3] = { 0.5F, 0.5F, 1.0F };
              float toggleColor2[3] = { 0.4F, 0.4F, 0.6F };
              float toggleColor3[3] = { 0.6F, 0.6F, 0.8F };
              float toggleDarkEdge[3] = { 0.3F, 0.3F, 0.5F};
              float toggleLightEdge[3] = { 0.7F, 0.7F, 0.9F};
              
              glColor3fv(toggleColor);
              for(a=0;a<ExecOpCnt;a++)
                {
                  switch(a) {
                  case 0:
                    /*
                      glColor3fv(toggleColor);
                      glBegin(GL_POLYGON);
                      glVertex2i(x2,y2+(ExecToggleSize)/2);
                      glVertex2i(x2+(ExecToggleSize)/2,y2);
                      glVertex2i(x2+ExecToggleSize,y2+(ExecToggleSize)/2);
                      glVertex2i(x2+(ExecToggleSize)/2,y2+ExecToggleSize);
                      glEnd();
                    */
                    
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor);
                    
                    draw_button_char(G,x2,y2+text_lift,'A');
                    break;
                  case 1:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor3);
                    
                    draw_button_char(G,x2,y2+text_lift,'S');
                    break;
                  case 2:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor2);
                    draw_button_char(G,x2,y2+text_lift,'H');
                    break;
                  case 3:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor);
                    draw_button_char(G,x2,y2+text_lift,'L');
                    break;
                  case 4:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                NULL);
                    draw_button_char(G,x2,y2+text_lift,'C');
                    break;
                  }
                  x2+=ExecToggleWidth;
                }
            }
#endif
            TextSetColor(G,I->Block->TextColor);
            TextSetPos2i(G,x+2,y2+text_lift);
            if((rec->type==cExecObject)||(rec->type==cExecAll)||
               (rec->type==cExecSelection))
              {
                y2=y;
                x2 = xx;
                if((x-ExecToggleMargin)-(xx-ExecToggleMargin)>-10) {
                  x2 = x+10;
                }
                if(rec->hilight||(row==I->Over)) {
                  draw_button(x,y2,(x2-x)-1,(ExecLineHeight-1),lightEdge,darkEdge,pressedColor);
                } else if(rec->visible) {
                  draw_button(x,y2,(x2-x)-1,(ExecLineHeight-1),lightEdge,darkEdge,enabledColor);
                } else {
                  draw_button(x,y2,(x2-x)-1,(ExecLineHeight-1),lightEdge,darkEdge,disabledColor);
                }

                TextSetColor(G,I->Block->TextColor);

                if(rec->type!=cExecObject)
                  c=rec->name;
                else 
                  c=rec->obj->Name;

                if(rec->type==cExecSelection)
                  if((nChar--)>0) {
                    TextDrawChar(G,'(');
                  }
              }

            if(c)
              while(*c) {
                if((nChar--)>0)
                  TextDrawChar(G,*(c++));
                else
                  break;
              }

            if(rec->type==cExecSelection)
              {
                if((nChar--)>0) {
                  TextDrawChar(G,')');
                }

                c=rec->name;
              }

            if(rec->type==cExecObject) {
              if(rec->obj->fGetCaption)
                c = rec->obj->fGetCaption(rec->obj);
              if(c && c[0] && nChar>1 && strcmp(c,rec->obj->Name)!=0)
                {
                  
                  TextSetColor(G,captionColor);
                  TextSetPos2i(G,x+2+8*(max_char-nChar),y2+text_lift);
                  if((nChar--)>0)
                    TextDrawChar(G,' ');
                  while(*c) 
                    if((nChar--)>0) 
                      TextDrawChar(G,*(c++));
                    else
                      break;
                }
            }

            y-=ExecLineHeight;
            if(y<(I->Block->rect.bottom))
              break;
          }
        }
    I->HowFarDown = y;
  }
}
/*========================================================================*/
int ExecutiveIterateObject(PyMOLGlobals *G,CObject **obj,void **hidden)
{
  int result;
  register CExecutive *I = G->Executive;
  int flag=false;
  SpecRec **rec=(SpecRec**)hidden;
  while(!flag)
	 {
		result = (int)ListIterate(I->Spec,(*rec),next);
		if(!(*rec))
		  flag=true;
		else if((*rec)->type==cExecObject)
		  flag=true;
	 }
  if(*rec)
	 (*obj)=(*rec)->obj;
  else
	 (*obj)=NULL;
  return(result);
}
/*========================================================================*/
int ExecutiveIterateObjectMolecule(PyMOLGlobals *G,ObjectMolecule **obj,void **hidden)
{
  int result;
  register CExecutive *I = G->Executive;
  int flag=false;
  SpecRec **rec=(SpecRec**)hidden;
  while(!flag)
	 {
		result = (int)ListIterate(I->Spec,(*rec),next);
		if(!(*rec))
		  flag=true;
		else if((*rec)->type==cExecObject)
        if((*rec)->obj->type==cObjectMolecule)
          flag=true;
	 }
  if(*rec)
	 (*obj)=(ObjectMolecule*)(*rec)->obj;
  else
	 (*obj)=NULL;
  return(result);
}
/*========================================================================*/
static void ExecutiveReshape(Block *block,int width,int height)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;

  BlockReshape(block,width,height);

  I->Width = block->rect.right-block->rect.left+1;
  I->Height = block->rect.top-block->rect.bottom+1;
  
}

/*========================================================================*/
int ExecutiveReinitialize(PyMOLGlobals *G)
{ 
  int ok=true;
  int blocked = false;
  /* reinitialize PyMOL */

  ExecutiveDelete(G,cKeywordAll);
  ColorReset(G);
  SettingInitGlobal(G,false,false);
  MovieReset(G);
  EditorInactivate(G);
  ControlRock(G,0);

  blocked = PAutoBlock();
  PRunString("cmd.view('*','clear')");
  PRunString("cmd.scene('*','clear')");
  WizardSet(G,NULL,false);
  PAutoUnblock(blocked);
  
  SculptCachePurge(G);
  SceneReinitialize(G);
  SelectorReinit(G);
  SeqChanged(G);

  return(ok);
}
/*========================================================================*/
int ExecutiveInit(PyMOLGlobals *G)
{
  register CExecutive *I=NULL;
  if( (I=(G->Executive=Calloc(CExecutive,1)))) {
    
  SpecRec *rec = NULL;
  int a;

  ListInit(I->Spec);
  I->Block = OrthoNewBlock(G,NULL);  
  I->Block->fRelease = ExecutiveRelease;
  I->Block->fClick   = ExecutiveClick;
  I->Block->fDrag    = ExecutiveDrag;
  I->Block->fDraw    = ExecutiveDraw;
  I->Block->fReshape = ExecutiveReshape;
  I->Block->active = true;
  I->ScrollBarActive = 0;
  I->ScrollBar=ScrollBarNew(G,false);
  OrthoAttach(G,I->Block,cOrthoTool);
  I->Pressed = -1;
  I->Over = -1;
  I->LastEdited=NULL;
  I->ReorderFlag=false;
  I->NSkip=0;
  I->HowFarDown=0;
  I->sizeFlag=false;

  ListElemCalloc(G,rec,SpecRec);
  strcpy(rec->name,"(all)");
  rec->type=cExecAll;
  rec->visible=true;
  rec->next=NULL;
  for(a=0;a<cRepCnt;a++)
	 rec->repOn[a]=false;
  ListAppend(I->Spec,rec,next,SpecRec);
  return 1;
  }
  else return 0;

}
/*========================================================================*/
void ExecutiveFree(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec=NULL;
  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->type==cExecObject)
		  rec->obj->fFree(rec->obj);
	 }
  ListFree(I->Spec,next,SpecRec);
  if(I->ScrollBar)
    ScrollBarFree(I->ScrollBar);
  OrthoFreeBlock(G,I->Block);
  I->Block=NULL;
  FreeP(G->Executive);
}


#ifdef _undefined

matrix checking code...

		double mt[3][3],mt2[3][3],pr[3][3],im[3][3],em[3][3];
	 printf("normalized matrix \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%12.3f ",evect[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 printf("tensor \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%12.3f ",mi[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	  for(a=0;a<3;a++) {
		 for(b=0;b<3;b++) {
			mt[a][b]=evect[a][b];
		 }
	  }

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  mt2[a][b]=evect[b][a];
		}
	 }

	 matrix_multiply33d33d(mt,mt2,pr);
	 printf("self product 1 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 matrix_multiply33d33d(mt,mi,im);
	 matrix_multiply33d33d(im,mt2,pr);
	 printf("diagonal product 1 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  em[a][b]=0.0;
		}
		em[a][a]=egval[a];
	 }

	 matrix_multiply33d33d(mt2,em,im);
	 matrix_multiply33d33d(im,mt,pr);
	 printf("diagonal product 4 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");
#endif


