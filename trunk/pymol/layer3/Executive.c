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

#include"Version.h"
#include"main.h"
#include"Base.h"
#include"OOMac.h"
#include"Executive.h"
#include"ObjectMesh.h"
#include"ObjectDist.h"
#include"ObjectSurface.h"
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

#include"Menu.h"
#include"Map.h"
#include"Editor.h"
#include"RepDot.h"

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
} SpecRec; /* specification record (a line in the executive window) */

ListVarDeclare(SpecList,SpecRec); 
/* NOTE: these vars are only used within calls -- not between them
 * However, they should be replaced with local vars declared inside
 * macros (assuming that's legal with all C compilers) */

typedef struct Executive {
  Block *Block;
  SpecRec *Spec;
  int Width,Height;
  int ScrollBarActive;
  int NSkip;
  struct CScrollBar *ScrollBar;
  CObject *LastEdited;
} CExecutive;

CExecutive Executive;
SpecRec *ExecutiveFindSpec(char *name);


int ExecutiveClick(Block *block,int button,int x,int y,int mod);
int ExecutiveRelease(Block *block,int button,int x,int y,int mod);
int ExecutiveDrag(Block *block,int x,int y,int mod);
void ExecutiveDraw(Block *block);
void ExecutiveReshape(Block *block,int width,int height);
int ExecutiveGetMaxDistance(char *name,float *pos,float *dev,int transformed,int state);
void ExecutiveObjMolSeleOp(int sele,ObjectMoleculeOpRec *op);

#define ExecLineHeight 14
#define ExecTopMargin 0
#define ExecToggleMargin 2
#define ExecLeftMargin 5
#define ExecRightMargin 2
#define ExecToggleWidth 14
#define ExecToggleSize 13

#define ExecOpCnt 5
#define ExecGreyVisible 0.45F
#define ExecGreyHidden 0.3F

typedef struct { 
  M4XAnnoType m4x;
  ObjectMolecule *obj;
} ProcPDBRec;

int ExecutiveGetObjectColorIndex(char *name)
{
  int result = -1;
  CObject *obj = NULL;
  obj = ExecutiveFindObjectByName(name);
  if(obj) {
    result=obj->Color;
  }
  return(result);
}

void ExecutiveProcessPDBFile(CObject *origObj,char *fname, char *oname, 
                             int frame, int discrete,int finish,OrthoLineType buf)
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


  f=fopen(fname,"rb");
  if(!f) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "ObjectMolecule-ERROR: Unable to open file '%s'\n",fname
      ENDFB;
    ok=false;
  } else
	 {
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeLoadPDBFile: Loading from %s.\n",fname
        ENDFB;
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);
    }
  if(ok) {
    processed = VLACalloc(ProcPDBRec,10);
  }
  while(repeat_flag&&ok) {
    char *start_at = buffer;
    int is_repeat_pass = false;
    ProcPDBRec *current;

    VLACheck(processed,ProcPDBRec,n_processed);
    current = processed + n_processed;

    M4XAnnoInit(&current->m4x);
    PRINTFD(FB_CCmd) " ExecutiveProcessPDBFile-DEBUG: loading PDB\n" ENDFD;

    if(next_pdb) {
      start_at = next_pdb;
      is_repeat_pass = true;
    }

    repeat_flag=false;
    next_pdb = NULL;
    if(!origObj) {

      obj=(CObject*)ObjectMoleculeReadPDBStr((ObjectMolecule*)origObj,
                                             start_at,frame,discrete,&current->m4x,pdb_name,&next_pdb);
      if(obj) {
        if(next_pdb) { /* NOTE: if set, assume that multiple PDBs are present in the file */
          repeat_flag=true;
        }
        if(current->m4x.xname_flag) { /* USER XNAME trumps the PDB Header name */                 
          ObjectSetName(obj,current->m4x.xname); /* from PDB */
          ExecutiveDelete(current->m4x.xname); /* just in case */
        } else if(next_pdb) {
          ObjectSetName(obj,pdb_name); /* from PDB */
          ExecutiveDelete(pdb_name); /* just in case */
        } else if(is_repeat_pass) {
          ObjectSetName(obj,pdb_name); /* from PDB */
          ExecutiveDelete(pdb_name); /* just in case */
        } else 
          ObjectSetName(obj,oname); /* from filename/parameter */

        ExecutiveManageObject(obj,true,false);
        if(frame<0)
          frame = ((ObjectMolecule*)obj)->NCSet-1;
        sprintf(buf," CmdLoad: \"%s\" loaded into object \"%s\", state %d.\n",
                fname,oname,frame+1);
      }
    } else {
      ObjectMoleculeReadPDBStr((ObjectMolecule*)origObj,
                               start_at,frame,discrete,&current->m4x,pdb_name,&next_pdb);
      if(finish)
        ExecutiveUpdateObjectSelection(origObj);
      if(frame<0)
        frame = ((ObjectMolecule*)origObj)->NCSet-1;
      sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
              fname,oname,frame+1);
      obj = origObj;
    }
    if(obj) {
      current->obj = (ObjectMolecule*)obj;
      n_processed++;
    }
  }

  /* BEGIN METAPHORICS ANNOTATION AND ALIGNMENT CODE */

  if(ok&&n_processed) { /* first, perform any Metaphorics alignment */
    /* is there a target structure? */
    {
      int a;
      for(a=0; a<n_processed; a++) {
        ProcPDBRec *current = processed + a;
        M4XAnnoType *m4x = &current->m4x;

        if(m4x->annotated_flag && m4x->align) {
          if(WordMatchExact(current->obj->Obj.Name,m4x->align->target,true)) {
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
            ObjectMoleculeConvertIDsToIndices(current->obj, m4x->align->id_at_point, m4x->align->n_point);
          }
        }
      }

      /* next, peform the alignments against the target */

      {
        int a;
        char aligned_name[] = "m4x_aligned";
        char tmp_sele[ObjNameMax*3];

        SelectorCreateEmpty("m4x_aligned");
                
                
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
                  SelectorGetUniqueTmpName(trg_sele);
                  SelectorGetUniqueTmpName(mbl_sele);
                  
                  SelectorCreateOrderedFromObjectIndices(trg_sele,pairwise.trg_obj,
                                                         pairwise.trg_vla,pairwise.n_pair);
                  SelectorCreateOrderedFromObjectIndices(mbl_sele,pairwise.mbl_obj,
                                                         pairwise.mbl_vla,pairwise.n_pair);
                  
                  sprintf(align_name,"%s_%s_alignment",
                          pairwise.trg_obj->Obj.Name,
                          pairwise.mbl_obj->Obj.Name);

                  ExecutiveRMS(mbl_sele, trg_sele, 2, 0.0F, 0, 0, align_name, 0, 0, true);
                  ExecutiveColor(align_name,"white",0,true);
                  if(target_rec->m4x.invisible) 
                    sprintf(tmp_sele, "(%s) | (%s)",aligned_name,mbl_sele);
                  else 
                    sprintf(tmp_sele, "(%s) | (%s) | (%s)",aligned_name, trg_sele, mbl_sele);
                  SelectorCreateSimple(aligned_name, tmp_sele);
                  /*ExecutiveDelete(trg_sele);
                    ExecutiveDelete(mbl_sele);*/
                }

              }
              ObjMolPairwisePurge(&pairwise);
            }
          }
        }
        sprintf(tmp_sele, "bychain %s", aligned_name);
        SelectorCreateSimple(aligned_name, tmp_sele);
      }
    }
  }
  
  if(ok&&n_processed) { /* next, perform any and all Metaphorics annotations */
      int a;
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
              ObjectMoleculeM4XAnnotate(current->obj,&current->m4x,script_file,(m4x_mode==1));
          M4XAnnoPurge(&current->m4x);
        }
      }
  }
  /* END METAPHORICS ANNOTATION AND ALIGNMENT CODE */
  
  VLAFreeP(processed);
  if(buffer) {
    mfree(buffer);
  }
  
}


int  ExecutiveAssignSS(char *target,int state,char *context,int preserve,int quiet)
{
  int sele0=-1;
  int sele1=-1;
  int ok = false;
  sele0 = SelectorIndexByName(target);
  if(sele0>=0) {
    if(!context[0]) {
      sele1=sele0;
    } else {
      sele1 = SelectorIndexByName(context);
    }

    if(sele1>=0) {
      ok =  SelectorAssignSS(sele0,sele1,state,preserve,quiet);
    }
  }
  return(ok);
}

PyObject *ExecutiveGetVisAsPyDict(void)
{
  PyObject *result=NULL,*list,*repList;
  CExecutive *I = &Executive;
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
}

int ExecutiveSetVisFromPyDict(PyObject *dict)
{
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

    SceneObjectDel(NULL);
    while (PyDict_Next(dict, &pos, &key, &list)) {
      if(!PConvPyStrToStr(key,name,sizeof(WordType))) {
        ok=false;
      } else {
        
        rec = ExecutiveFindSpec(name);
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
            SceneObjectAdd(rec->obj); 
        }
      }
    }
  }
  return ok;
}

int ExecutiveIsolevel(char *name,float level,int state)
{
  int ok =true;
  CObject *obj;
  obj = ExecutiveFindObjectByName(name);
  if(obj) {
    switch(obj->type) {
    case cObjectMesh:
      ObjectMeshSetLevel((ObjectMesh*)obj,level,state);
        SceneChanged();
      break;
    case cObjectSurface:
      break;
    default:
      ok=false;
      PRINTFB(FB_Executive,FB_Errors)
        " Isolevel-Error: object \"%s\" is of wrong type.",name
        ENDFB;
      break;
    }
  }
  return(ok);

}

int ExecutiveSpectrum(char *s1,char *expr,float min,float max,int first,int last,
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

  sele1 = SelectorIndexByName(s1);
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
        color_index[a] = ColorGetIndex(buffer);
      }
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_CountAtoms;
      op.i1=0;
      ExecutiveObjMolSeleOp(sele1,&op);
      n_atom=op.i1;
      
      if(n_atom) {
        value = Calloc(float,n_atom);
        
        if(WordMatch("count",expr,true)) {
          for(a=0;a<n_atom;a++) {
            value[a]=(float)a+1;
          }
        } else if(WordMatch("b",expr,true)) {
          op.code = OMOP_GetBFactors;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(sele1,&op);
        } else if(WordMatch("q",expr,true)) {
          op.code = OMOP_GetOccupancies;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(sele1,&op);
        } else if(WordMatch("pc",expr,true)) {
          op.code = OMOP_GetPartialCharges;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(sele1,&op);
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
          PRINTFB(FB_Executive,FB_Actions)
            " Spectrum: range (%8.5f to %8.5f).\n"
            ,min,max
            ENDFB;
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
        
        ExecutiveObjMolSeleOp(sele1,&op);

        op.code=OMOP_INVA;
        op.i1=cRepAll; 
        op.i2=cRepInvColor;
        ExecutiveObjMolSeleOp(sele1,&op);

      }
    }

    FreeP(color_index);
    FreeP(value);
  }
  return(ok);
}

char *ExecutiveGetChains(char *sele,int state,int *null_chain)
{
  int sele1;
  char *result = NULL;
  int chains[256];
  int a,c;
  ObjectMoleculeOpRec op;

  sele1 = SelectorIndexByName(sele);
  if(sele1>=0) {

    for(a=0;a<256;a++) {
      chains[a]=0;
    }
    ObjectMoleculeOpRecInit(&op);
    op.code=OMOP_GetChains;
    op.ii1 = chains;
    op.i1=0;
    ExecutiveObjMolSeleOp(sele1,&op);
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
    ErrMessage("ExecutiveGetChains","Bad selection.");
  }
  return(result);
}

int ExecutiveValidateObjectPtr(CObject *ptr,int object_type)
{
  CExecutive *I = &Executive;
  int ok=false;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      if(rec->obj->type==object_type) {
        ok=true;
        break;
      }
    }
  }
  return(ok);
}

int ExecutiveRampMapNew(char *name,char *map_name,PyObject *range,PyObject *color,int map_state)
{
  ObjectGadgetRamp *obj = NULL;
  int ok =true;
  CObject *map_obj;

  map_obj = ExecutiveFindObjectByName(map_name);
  if(map_obj) {
    if(map_obj->type!=cObjectMap) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExecutiveRampMapNew: Error: object '%s' is not a map.\n",map_name
        ENDFB;
      ok=false;
    }
  } else {
    PRINTFB(FB_Executive,FB_Errors)
      "ExecutiveRampMapNew: Error: map '%s' not found.\n",map_name
      ENDFB;
    ok = false;
  }
  ok = ok && (obj=ObjectGadgetRampMapNewAsDefined((ObjectMap*)map_obj,range,color,map_state));
  if(ok) ExecutiveDelete(name); 
  if(ok) ObjectSetName((CObject*)obj,name);
  if(ok) ColorRegisterExt(name,(void*)obj,cColorGadgetRamp);
  if(ok) ExecutiveManageObject((CObject*)obj,false,false);
  return(ok);
}

static int ExecutiveCountNames(void)
{
  int count=0;
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next))
	 count++;
  
  return(count);
}

static PyObject *ExecutiveGetExecObject(SpecRec *rec)
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

static int ExecutiveSetNamedEntries(PyObject *names,int version)
{
  CExecutive *I = &Executive;  
  int ok=true;
  int skip=false;
  int a=0,l=0;
  PyObject *cur;
  SpecRec *rec = NULL;
  int extra_int;

  if(ok) ok = (names!=NULL);
  if(ok) ok = PyList_Check(names);
  if(ok) l = PyList_Size(names);
  while(ok&&(a<l)) {
    cur = PyList_GetItem(names,a);
    if(cur!=Py_None) { /* skip over None w/o aborting */
      skip=false;
      rec=NULL;
      ListElemAlloc(rec,SpecRec); 
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
          if(ok) ok = ObjectMoleculeNewFromPyList(PyList_GetItem(cur,5),(ObjectMolecule**)&rec->obj);
          break;
        case cObjectDist:
          if(ok) ok = ObjectDistNewFromPyList(PyList_GetItem(cur,5),(ObjectDist**)&rec->obj);
          break;
        case cObjectMap:
          if(ok) ok = ObjectMapNewFromPyList(PyList_GetItem(cur,5),(ObjectMap**)&rec->obj);
          break;
        case cObjectMesh:
          if(ok) ok = ObjectMeshNewFromPyList(PyList_GetItem(cur,5),(ObjectMesh**)&rec->obj);
          break;
        case cObjectSurface:
          if(ok) ok = ObjectSurfaceNewFromPyList(PyList_GetItem(cur,5),(ObjectSurface**)&rec->obj);
          break;
        case cObjectCGO:
          if(ok) ok = ObjectCGONewFromPyList(PyList_GetItem(cur,5),(ObjectCGO**)&rec->obj,version);
          break;
        case cObjectGadget:
          if(ok) ok = ObjectGadgetNewFromPyList(PyList_GetItem(cur,5),(ObjectGadget**)&rec->obj,version);
          break;
        default:
          PRINTFB(FB_Executive,FB_Errors)
            " Executive:  unrecognized object \"%s\" of type %d.\n",
            rec->name,rec->type
            ENDFB;
          skip=true;
          break;
        }
        break;
      case cExecSelection: /* on the first pass, just create an entry in the rec list */
        rec->sele_color=extra_int;
        break;
      }

      if(PyErr_Occurred()) {
        PRINTFB(FB_Executive,FB_Errors)
          "ExectiveSetNamedEntries-ERROR: after object \"%s\".\n",rec->name
          ENDFB;
        PyErr_Print();  
      }

      if(ok&&!skip) {
        switch(rec->type) {
        case cExecObject:        
          if(rec->visible) {
            SceneObjectAdd(rec->obj);
          }
          ExecutiveUpdateObjectSelection(rec->obj);
          break;
        }
        ListAppend(I->Spec,rec,next,SpecList);
      } else {
        ListElemFree(rec);
      }
    }
    a++;
  }
  return(ok);
}

static int ExecutiveSetSelections(PyObject *names)
{
  /* must already have objects loaded at this point... */

  int ok=true;
  int a=0,l=0;
  PyObject *cur;
  SpecRec *rec = NULL;
  int extra;

  if(ok) ok = (names!=NULL);
  if(ok) ok = PyList_Check(names);
  if(ok) l = PyList_Size(names);
  while(ok&&(a<l)) {
    cur = PyList_GetItem(names,a);
    if(cur!=Py_None) { /* skip over None w/o aborting */
      rec=NULL;
      ListElemAlloc(rec,SpecRec); 
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
        ok = SelectorFromPyList(rec->name,PyList_GetItem(cur,5));
        break;
      }
      ListElemFree(rec);
    }
    a++;
  }
  return(ok);
}

static PyObject *ExecutiveGetExecSelePyList(SpecRec *rec)
{
  PyObject *result = NULL;
  int sele;

  sele = SelectorIndexByName(rec->name);
  if(sele>=0) {
    result = PyList_New(6);
    PyList_SetItem(result,0,PyString_FromString(rec->name));
    PyList_SetItem(result,1,PyInt_FromLong(cExecSelection));
    PyList_SetItem(result,2,PyInt_FromLong(rec->visible));
    PyList_SetItem(result,3,PConvIntArrayToPyList(rec->repOn,cRepCnt));
    PyList_SetItem(result,4,PyInt_FromLong(-1));
    PyList_SetItem(result,5,SelectorAsPyList(sele));
  }
  return(PConvAutoNone(result));
}

static PyObject *ExecutiveGetNamedEntries(void)
{
  CExecutive *I = &Executive;  
  PyObject *result = NULL;
  int count;
  SpecRec *rec = NULL;

  count = ExecutiveCountNames();
  result = PyList_New(count);

  SelectorUpdateTable();

  count=0;
  while(ListIterate(I->Spec,rec,next))
	 {
      switch(rec->type) {
      case cExecObject:
        PyList_SetItem(result,count,
                       ExecutiveGetExecObject(rec));
        break;
      case cExecSelection:
        PyList_SetItem(result,count,
                       ExecutiveGetExecSelePyList(rec));
        break;
      default:
        PyList_SetItem(result,count,PConvAutoNone(NULL));
        break;
      }
      count++;
    }
  return(PConvAutoNone(result));
}

int ExecutiveGetSession(PyObject *dict)
{
  int ok=true;
  SceneViewType sv;
  PyObject *tmp;

  tmp = ExecutiveGetNamedEntries();
  PyDict_SetItemString(dict,"names",tmp);
  Py_XDECREF(tmp);

  tmp = SelectorSecretsAsPyList();
  PyDict_SetItemString(dict,"selector_secrets",tmp);
  Py_XDECREF(tmp);
  
  tmp = SettingGetGlobalsPyList();
  PyDict_SetItemString(dict,"settings",tmp);
  Py_XDECREF(tmp);

  tmp = ColorAsPyList();
  PyDict_SetItemString(dict,"colors",tmp);
  Py_XDECREF(tmp);

  tmp = ColorExtAsPyList();
  PyDict_SetItemString(dict,"color_ext",tmp);
  Py_XDECREF(tmp);

  tmp = PyInt_FromLong(_PyMOL_VERSION_int);
  PyDict_SetItemString(dict,"version",tmp);
  Py_XDECREF(tmp);

  SceneGetView(sv);
  tmp = PConvFloatArrayToPyList(sv,cSceneViewSize);
  PyDict_SetItemString(dict,"view",tmp);
  Py_XDECREF(tmp);

  tmp = MovieAsPyList();
  PyDict_SetItemString(dict,"movie",tmp);
  Py_XDECREF(tmp);

  tmp = EditorAsPyList();
  PyDict_SetItemString(dict,"editor",tmp);
  Py_XDECREF(tmp);

  tmp = MainAsPyList();
  PyDict_SetItemString(dict,"main",tmp);
  Py_XDECREF(tmp);
  return(ok);
}

int ExecutiveSetSession(PyObject *session)
{
  int ok=true;
  PyObject *tmp;
  SceneViewType sv;
  int version=-1;

  ExecutiveDelete("all");
  ColorReset();
  if(ok) ok = PyDict_Check(session);

  if(ok&&SettingGet(cSetting_session_version_check)) {
    tmp = PyDict_GetItemString(session,"version");
    if(tmp) {
      ok = PConvPyIntToInt(tmp,&version);
      if(ok) {
        if(version>_PyMOL_VERSION_int) {
          PRINTFB(FB_Executive,FB_Errors)
            "Error: This session was created with a newer version of PyMOL (%1.2f).\n",version/100.0
            ENDFB;
          PRINTFB(FB_Executive,FB_Errors)
            "Error: Please obtain a more recent version from http://www.pymol.org\n"
            ENDFB;
          ok=false;
        } else {
          PRINTFB(FB_Executive,FB_Details)          
            " Executive: Loading version %1.2f session...\n",
            version/100.0
            ENDFB;
        }
      }
    }
  }

  if(ok) {
    tmp = PyDict_GetItemString(session,"colors");
    if(tmp) {
      ok = ColorFromPyList(tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after colors.\n"
        ENDFB;
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"color_ext");
    if(tmp) {
      ok = ColorExtFromPyList(tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after color_ext.\n"
        ENDFB;
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"settings");
    if(tmp) {
      ok = SettingSetGlobalsFromPyList(tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after settings.\n"
        ENDFB;
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"names");
    if(tmp) {
      if(ok) ok=ExecutiveSetNamedEntries(tmp,version);
      if(ok) ok=ExecutiveSetSelections(tmp);
    }
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after names.\n"
        ENDFB;
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"selector_secrets");
    if(tmp) {
      if(ok) ok=SelectorSecretsFromPyList(tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after selector secrets.\n"
        ENDFB;
    }
  }  
  if(ok) {
    tmp = PyDict_GetItemString(session,"view");
    if(tmp) {
      ok = PConvPyListToFloatArrayInPlace(tmp,sv,cSceneViewSize);
    }
    if(ok) SceneSetView(sv,true);
    
    
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after view.\n"
        ENDFB;
    }
  }
  
  if(ok) {
    int warning;
    tmp = PyDict_GetItemString(session,"movie");
    if(tmp) {
      ok = MovieFromPyList(tmp,&warning);
    }
    
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after movie.\n"
        ENDFB;
    }
  }
  
  if(ok) {
    tmp = PyDict_GetItemString(session,"editor");
    if(tmp) {
      ok = EditorFromPyList(tmp);
    }
    
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after editor.\n"
        ENDFB;
    }
  }
  if(ok) { /* update mouse in GUI */
    PParse("cmd.mouse(quiet=1)");
    PParse("viewport"); /* refresh window/internal_gui status */
  }
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
      PRINTFB(FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after main.\n"
        ENDFB;
    }
  }
  if(!ok) {
    PRINTFB(FB_Executive,FB_Warnings)
      "ExectiveSetSession-Warning: restore may be incomplete.\n"
      ENDFB;
  }
  return(ok);
}

#define ExecScrollBarMargin 2
#define ExecScrollBarWidth 13

void ExecutiveObjMolSeleOp(int sele,ObjectMoleculeOpRec *op);
CObject **ExecutiveSeleToObjectVLA(char *s1);

CObject **ExecutiveSeleToObjectVLA(char *s1)
{
  /* return VLA containing list of atoms references by selection */

  CObject **result = NULL;
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  CObject *obj=NULL;
  int n = 0;
  ObjectMoleculeOpRec op2;
  int sele;

  result = VLAlloc(CObject*,50);
  if(WordMatch(s1,cKeywordAll,true)) {
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
    sele = SelectorIndexByName(s1);
    if(sele>0) {
      ObjectMoleculeOpRecInit(&op2);
      op2.code=OMOP_GetObjects;
      op2.obj1VLA=(ObjectMolecule**)result;
      op2.i1=0;
      ExecutiveObjMolSeleOp(sele,&op2);
      n = op2.i1;
      result = (CObject**)op2.obj1VLA;
    } else {
      obj = ExecutiveFindObjectByName(s1);
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

int ExecutiveGetCrystal(char *sele,float *a,float *b,float *c,
                         float *alpha,float *beta,float *gamma,
                        char *sgroup,int *defined)
{
  int ok=true;

  ObjectMolecule *objMol;
  int sele0;
  sele0 = SelectorIndexByName(sele);
  *defined = false;
  if(sele0<0) {
    PRINTFB(FB_Executive,FB_Errors)
      "Error: invalid selection.\n"
      ENDFB;
    ok=false;
  } else {
    objMol = SelectorGetSingleObjectMolecule(sele0);
    if(!objMol) {
      PRINTFB(FB_Executive,FB_Errors)
        "Error: selection must refer to exactly one object.\n"
        ENDFB;
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

int ExecutiveSetCrystal(char *sele,float a,float b,float c,
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

  objVLA = ExecutiveSeleToObjectVLA(sele);
  n_obj = VLAGetSize(objVLA);
  if(n_obj) {
    for(i=0;i<n_obj;i++) {
      obj = objVLA[i];
      switch(obj->type) {
      case cObjectMolecule:
        if(!symmetry) {
          symmetry=SymmetryNew();          
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
          crystal = CrystalNew();
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
    PRINTFB(FB_Executive,FB_Errors)
      " ExecutiveSetCrystal: no object selected\n"
      ENDFB;
  }
  if(crystal)
    CrystalFree(crystal);
  if(symmetry)
    SymmetryFree(symmetry);
  VLAFreeP(objVLA);
  return(ok);
}

int ExecutiveSmooth(char *name,int cycles,int window,int first, int last, int ends)
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

  PRINTFD(FB_Executive)
    " ExecutiveSmooth: entered %s,%d,%d,%d,%d,%d\n",name,cycles,first,last,window,ends
    ENDFD;

  sele=SelectorIndexByName(name);



  if(sele>=0) {
    if(last<0) 
      last = ExecutiveCountStates(name)-1;
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
    
    PRINTFD(FB_Executive)
      " ExecutiveSmooth: first %d last %d n_state %d backward %d forward %d range %d\n",
      first,last,n_state,backward,forward,range
      ENDFD;

    if(n_state>=window) {
      
      /* determine storage req */
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_CountAtoms;
      op.i1=0;
      ExecutiveObjMolSeleOp(sele,&op);
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
        
        PRINTFB(FB_Executive,FB_Actions)
          " Smooth: copying coordinates to temporary arrays..\n"
          ENDFB;
        op.code = OMOP_CSetIdxGetAndFlag;
        op.i1 = n_atom; 
        op.i2 = 0;
        op.cs1 = first;
        op.cs2 = last;
        op.vv1 = coord0;
        op.ii1 = flag0;
        op.nvv1 = 0;          
        ExecutiveObjMolSeleOp(sele,&op);    
        
        PRINTFD(FB_Executive)  
          " ExecutiveSmooth: got %d %d\n",op.i2,op.nvv1
          ENDFD;
        
        UtilZeroMem(coord1,sizeof(float)*3*n_atom*n_state);
        UtilZeroMem(flag1,sizeof(int)*n_atom*n_state);
        
        for(a=0;a<cycles;a++) {                
          PRINTFB(FB_Executive,FB_Actions)
            " Smooth: smoothing (pass %d)...\n",a+1
            ENDFB;
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

        PRINTFB(FB_Executive,FB_Actions)
          " Smooth: updating coordinates...\n"
          ENDFB;
        
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
        
        ExecutiveObjMolSeleOp(sele,&op);      
        PRINTFD(FB_Executive)  
          " ExecutiveSmooth: put %d %d\n",op.i2,op.nvv1
          ENDFD;
        
        
        FreeP(coord0);
        FreeP(coord1);
        FreeP(flag0);
        FreeP(flag1);
      }
    }
  } else {
    PRINTFB(FB_Executive,FB_Errors)  
      " ExecutiveSmooth: selection not found\n"
      ENDFB;
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveDebug(char *name)
{
  ObjectMolecule *obj;
  ObjectMoleculeBPRec bp;
  int a;

  obj=(ObjectMolecule*)ExecutiveFindObjectByName(name);
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
int ***ExecutiveGetBondPrint(char *name,int max_bond,int max_type,int *dim)
{
  int ***result = NULL;
  CObject *obj;
  ObjectMolecule *objMol;

  obj=ExecutiveFindObjectByName(name);
  if(obj->type==cObjectMolecule) {
    objMol = (ObjectMolecule*)obj;
    result = ObjectMoleculeGetBondPrint(objMol,max_bond,max_type,dim);
  }
  return(result);
}
/*========================================================================*/
int ExecutiveMapNew(char *name,int type,float *grid,
                    char *sele,float buffer,
                    float *minCorner,
                    float *maxCorner,int state)
{
  CObject *origObj=NULL;
  ObjectMap *objMap;
  ObjectMapState *ms = NULL;
  int a;
  float v[3];
  ObjectMapDesc _md,*md;
  int ok = true;
  int sele0 = SelectorIndexByName(sele);
  int isNew=true;
  int n_state;
  int valid_extent=false;
  int st;
  int st_once_flag=true;
  int n_st;

  md=&_md;

  if(state==-2) state=SceneGetState();

  /* remove object if it already exists */

  origObj=ExecutiveFindObjectByName(name);

  if(origObj) {
    if(origObj->type!=cObjectMap) {
      ExecutiveDelete(origObj->Name);
    } else {
      isNew=false;
    }
  }

  n_st = ExecutiveCountStates(NULL);

  for(st=0;st<n_st;st++) {
    if(state==-1) st_once_flag=false; /* each state, separate map, separate extent */
    if(!st_once_flag) state=st;
    
    if(strlen(sele)) {
      valid_extent = ExecutiveGetExtent(sele,md->MinCorner,
                                        md->MaxCorner,true,state,false); /* TODO restrict to state */
    } else {
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
        objMap = ObjectMapNew();
      else
        objMap = (ObjectMap*)origObj;
      if(objMap) {
        int once_flag=true;
        n_state = SelectorCountStates(sele0);
        if(valid_extent)
          for(a=0;a<n_state;a++) {
            if(state==-3) once_flag=false; /* -2 = each state, separate map, shared extent */
            if(state==-4) state=-1; /* all states, one map */
            if(!once_flag) state=a;
            ms = ObjectMapNewStateFromDesc(objMap,md,state);
            if(!ms)
              ok=false;
          
            if(ok&&ms) {
            
              switch(type) {
              case 0: /* vdw */
                SelectorMapMaskVDW(sele0,ms,0.0F,state);
                break;
              case 1: /* coulomb */
                SelectorMapCoulomb(sele0,ms,50.0F,state);
                break;
              case 2: /* gaussian */
                SelectorMapGaussian(sele0,ms,0.0F,state);
                break;
              }
              if(!ms->Active)
                ObjectMapStatePurge(ms);
            }
            if(once_flag) break;
          }

        ObjectSetName((CObject*)objMap,name);
        ObjectMapUpdateExtents(objMap);
        if(isNew)
          ExecutiveManageObject((CObject*)objMap,true,false);
        isNew=false;
        origObj = (CObject*)objMap;
      }
      SceneDirty();
    }
    if(st_once_flag)
      break;
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSculptIterateAll(void)
{
  int active = false;

  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;

  int state = SceneGetState();
  int cycles = (int)SettingGet(cSetting_sculpting_cycles);

  if(SettingGet(cSetting_sculpting)) {
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
float ExecutiveSculptIterate(char *name,int state,int n_cycle)
{
  CObject *obj = ExecutiveFindObjectByName(name);
  CExecutive *I = &Executive;
  int ok=true;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  float total_strain = 0.0F;

  if(state<0) state=SceneGetState();

  if(WordMatch(name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          total_strain+=ObjectMoleculeSculptIterate(objMol,state,n_cycle);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB;
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB;
    ok=false;
  } else {
    total_strain=ObjectMoleculeSculptIterate((ObjectMolecule*)obj,state,n_cycle);
  }
  return(total_strain);
}
/*========================================================================*/
int ExecutiveSculptActivate(char *name,int state)
{
  CObject *obj = ExecutiveFindObjectByName(name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  CExecutive *I = &Executive;
  int ok=true;
  if(state<0) state=SceneGetState();

  if(WordMatch(name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptImprint(objMol,state);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB;
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB;
    ok=false;
  } else {
    ObjectMoleculeSculptImprint((ObjectMolecule*)obj,state);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSculptDeactivate(char *name)
{
  CObject *obj = ExecutiveFindObjectByName(name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  CExecutive *I = &Executive;

  int ok=true;

  if(WordMatch(name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptClear(objMol);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB;
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB;
    ok=false;
  } else {
    ObjectMoleculeSculptClear((ObjectMolecule*)obj);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSetGeometry(char *s1,int geom,int valence)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int ok=false;

  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_SetGeometry;
    op1.i1 = geom;
    op1.i2 = valence;
    op1.i3 = 0;
    ExecutiveObjMolSeleOp(sele1,&op1);
    if(op1.i3) ok=true;
  } else {
    ErrMessage("SetGeometry","Invalid selection.");
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveMultiSave(char *fname,char *name,int state,int append)
{
  int result=false;
  SpecRec *tRec;
  ObjectMolecule *objMol;
  
  PRINTFD(FB_Executive)
    " ExecutiveMultiSave-Debug: entered %s %s.\n",fname,name
    ENDFD;
  tRec = ExecutiveFindSpec(name);
  if(tRec) {
    if(tRec->type==cExecObject)
      if(tRec->obj->type==cObjectMolecule) {
        objMol =(ObjectMolecule*)tRec->obj;
        result = ObjectMoleculeMultiSave(objMol,fname,state,append);
      }
  }
  return(result);
  
}
int ExecutiveMapSetBorder(char *name,float level)
{
  int result=false;
  SpecRec *tRec;
  ObjectMap *mobj;

  tRec = ExecutiveFindSpec(name);
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

int ExecutiveMapDouble(char *name,int state)
{
  int result=false;
  SpecRec *tRec;
  ObjectMap *mobj;

  tRec = ExecutiveFindSpec(name);
  if(tRec) {
    if(tRec->type==cExecObject)
      if(tRec->obj->type==cObjectMap) {
        mobj =(ObjectMap*)tRec->obj;
        result = ObjectMapDouble(mobj,state);
      }
  }
  return(result);
}

void ExecutiveSelectRect(BlockRect *rect,int mode)
{
  Multipick smp;
  OrthoLineType buffer,buf2;
  char prefix[3]="";
  int log_box = 0;
  int logging;
  logging = (int)SettingGet(cSetting_logging);
  if(logging)
    log_box= (int)SettingGet(cSetting_log_box_selections);
  if(logging==cPLog_pml)
    strcpy(prefix,"_ ");
  smp.picked=VLAlloc(Pickable,1000);
  smp.x=rect->left;
  smp.y=rect->bottom;
  smp.w=rect->right-rect->left;
  smp.h=rect->top-rect->bottom;
  SceneMultipick(&smp);
  if(smp.picked[0].index) {
    SelectorCreate(cTempRectSele,NULL,NULL,1,&smp);
    if(log_box) SelectorLogSele(cTempRectSele);
    if(mode==cButModeRect) {
      SelectorCreate(cLeftButSele,cTempRectSele,NULL,1,NULL);
      if(log_box) {
        sprintf(buf2,"%scmd.select(\"%s\",\"%s\",quiet=1)\n",prefix,cLeftButSele,cTempRectSele);
        PLog(buf2,cPLog_no_flush);
      }
    } else if(SelectorIndexByName(cLeftButSele)>=0) {
      if(mode==cButModeRectAdd) {
        sprintf(buffer,"(%s or %s)",cLeftButSele,cTempRectSele);
        SelectorCreate(cLeftButSele,buffer,NULL,0,NULL);
        if(log_box) {
          sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,cLeftButSele,buffer);
          PLog(buf2,cPLog_no_flush);
        }
      } else {
        sprintf(buffer,"(%s and not %s)",cLeftButSele,cTempRectSele);
        SelectorCreate(cLeftButSele,buffer,NULL,0,NULL);
        if(log_box) {
          sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,cLeftButSele,buffer);
          PLog(buf2,cPLog_no_flush);
        }
      }
    } else {
      if(mode==cButModeRectAdd) {
        SelectorCreate(cLeftButSele,cTempRectSele,NULL,0,NULL);
        if(log_box) {
          sprintf(buf2,"%scmd.select(\"%s\",\"%s\")\n",prefix,cLeftButSele,cTempRectSele);
          PLog(buf2,cPLog_no_flush);
        }
      } else {
        SelectorCreate(cLeftButSele,"(none)",NULL,0,NULL);
        if(log_box) {
          sprintf(buf2,"%scmd.select(\"%s\",\"(none)\")\n",prefix,cLeftButSele);
          PLog(buf2,cPLog_no_flush);
        }
      }
    }
    if(SettingGet(cSetting_auto_show_selections)) {
      ExecutiveSetObjVisib(cLeftButSele,true);
    }
    if(log_box) {
      sprintf(buf2,"%scmd.delete(\"%s\")\n",prefix,cTempRectSele);
      PLog(buf2,cPLog_no_flush);
      PLogFlush();
    }
    ExecutiveDelete(cTempRectSele);
  }
  VLAFreeP(smp.picked);
  WizardDoSelect(cLeftButSele);
}

int ExecutiveTranslateAtom(char *sele,float *v,int state,int mode,int log)
{
  int ok=true;
  ObjectMolecule *obj0;
  int sele0 = SelectorIndexByName(sele);
  int i0;
  if(sele0<0) {
    PRINTFB(FB_Executive,FB_Errors)
      "Error: bad selection %s.\n",sele
      ENDFB;
    ok=false;
  } else{ 
    obj0 = SelectorGetSingleObjectMolecule(sele0);
    if(!obj0) {
      PRINTFB(FB_Executive,FB_Errors)
        "Error: selection isn't a single atom.\n"
        ENDFB;
      ok=false;
    } else {
      i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
      if(i0<0) {
        PRINTFB(FB_Executive,FB_Errors)
          "Error: selection isn't a single atom.\n"
          ENDFB;
        ok=false;
      } else {
        ObjectMoleculeMoveAtom(obj0,state,i0,v,mode,log);
      }
    }
  }
  return(ok);
}

int ExecutiveCombineObjectTTT(char *name,float *ttt)
{
  CObject *obj = ExecutiveFindObjectByName(name);
  int ok=true;

  if(!obj) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB;
    ok=false;
  } else {
    ObjectCombineTTT(obj,ttt);
    SceneDirty();
  }
  return(ok);
}

int ExecutiveTransformSelection(int state,char *s1,int log,float *ttt)
{
  int sele=-1;
  ObjectMolecule *obj = NULL;
  ObjectMolecule **vla = NULL;
  int nObj;
  int ok=true;
  int a;

  sele = SelectorIndexByName(s1);
  if(sele<0)
    ok=false;
  if(ok) {
    vla=SelectorGetObjectMoleculeVLA(sele);
    if(!vla) ok=false;
  }
  if(ok) {
    nObj = VLAGetSize(vla);
    for(a=0;a<nObj;a++) {
      obj=vla[a];
      ObjectMoleculeTransformSelection(obj,state,sele,ttt,log,s1);
    }
  }
  SceneDirty();
  VLAFreeP(vla);
  return(ok);
}

int ExecutiveTransformObjectSelection(char *name,int state,char *s1,int log,float *ttt)
{
  int sele=-1;
  ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(name);
  int ok=true;

  if(s1[0]) {
    sele = SelectorIndexByName(s1);
    if(sele<0)
      ok=false;
  }
  if(!obj) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB;
  } else if(!ok) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: selection object %s not found.\n",s1
      ENDFB;
  } else {
    ObjectMoleculeTransformSelection(obj,state,sele,ttt,log,s1);
  }
  SceneDirty();
  return(ok);
}

int ExecutiveValidName(char *name)
{
  int result=true;

  if(!ExecutiveFindSpec(name)) {
    if(!WordMatch(name,cKeywordAll,true))
      if(!WordMatch(name,cKeywordSame,true))
        if(!WordMatch(name,cKeywordCenter,true))
          if(!WordMatch(name,cKeywordOrigin,true))
            result=false;
  }
  return result;
}

int ExecutivePhiPsi(char *s1,ObjectMolecule ***objVLA,int **iVLA,
                    float **phiVLA,float **psiVLA,int state) 
{
  int sele1=SelectorIndexByName(s1);
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
    ExecutiveObjMolSeleOp(sele1,&op1);
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


float ExecutiveAlign(char *s1,char *s2,char *mat_file,float gap,float extend,int skip,
                     float cutoff,int cycles,int quiet,char *oname,
                     int state1,int state2)
{
  int sele1=SelectorIndexByName(s1);
  int sele2=SelectorIndexByName(s2);
  int *vla1=NULL;
  int *vla2=NULL;
  int na,nb;
  int c;
  float result = 0.0;
  int ok=true;
  CMatch *match = NULL;

  if((sele1>=0)&&(sele2>=0)) {
    vla1=SelectorGetResidueVLA(sele1);
    vla2=SelectorGetResidueVLA(sele2);
    if(vla1&&vla2) {
      na = VLAGetSize(vla1)/3;
      nb = VLAGetSize(vla2)/3;
      if(na&&nb) {
        match = MatchNew(na,nb);
        if (ok) ok = MatchResidueToCode(match,vla1,na);
        if (ok) ok = MatchResidueToCode(match,vla2,nb);
        if (ok) ok = MatchMatrixFromFile(match,mat_file);
        if (ok) ok = MatchPreScore(match,vla1,na,vla2,nb);
        result = MatchAlign(match,gap,extend,skip);
        if(match->pair) { 
          c = SelectorCreateAlignments(match->pair,
                                       sele1,vla1,sele2,vla2,
                                       "_align1","_align2",false);
          if(c) {
            PRINTFB(FB_Executive,FB_Actions)
              " ExecutiveAlign: %d atoms aligned.\n",c
              ENDFB;
            result =ExecutiveRMS("_align1","_align2",2,cutoff,cycles,quiet,oname,
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

int ExecutivePairIndices(char *s1,char *s2,int state1,int state2,
                         int mode,float cutoff,float h_angle,
                         int **indexVLA, ObjectMolecule ***objVLA)
{
  int result = 0;
  int sele1,sele2;

  sele1 = SelectorIndexByName(s1);
  sele2 = SelectorIndexByName(s2);
  if((sele1>=0)&&(sele2>=0)) {
    result=SelectorGetPairIndices(sele1,state1,sele2,state2,
                                  mode,cutoff,h_angle,indexVLA,objVLA);
  } else {
    ErrMessage("ExecutivePairIndices","One or more bad selections.");
  }
  return(result);
}

void ExecutiveFocus(void)
{ /* unfortunately, this doesn't achieve the desired effect */
  if(PMGUI) {
    p_glutPopWindow();
    p_glutShowWindow();
  }
}

int ExecutiveCartoon(int type,char *s1)
{
  int sele1;
  ObjectMoleculeOpRec op1;

  sele1=SelectorIndexByName(s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2=0;
  if(sele1>=0) {
    op1.code=OMOP_INVA;
    op1.i1=cRepCartoon; 
    op1.i2=cRepInvRep;
    ExecutiveObjMolSeleOp(sele1,&op1);
    op1.code = OMOP_Cartoon;
    op1.i1 = type;
    op1.i2 = 0;
    ExecutiveObjMolSeleOp(sele1,&op1);
  } else {
    ErrMessage("Cartoon","Invalid selection.");
  }
  return(op1.i2);
}
/*========================================================================*/
float *ExecutiveGetVertexVLA(char *s1,int state)
{
  /* returns NULL if none found */

  float *result = NULL;
  ObjectMoleculeOpRec op1;
  int sele1;
  sele1 = SelectorIndexByName(s1);
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
    ExecutiveObjMolSeleOp(sele1,&op1);
    if(op1.nvv1) {
      VLASize(op1.vv1,float,op1.nvv1*3);
      result = op1.vv1;
    } else 
      VLAFreeP(op1.vv1);
  }
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveGetSettingText(int index,char *object,int state)
{ /* Assumes blocked Python interpreter */
  PyObject *result = NULL;
  OrthoLineType buffer = "";
  CObject *obj = NULL;
  CSetting **handle=NULL,*set_ptr1=NULL,*set_ptr2=NULL;
  int ok=true;

  if(object)
    if(object[0]) {
      obj=ExecutiveFindObjectByName(object);
      if(!obj) 
        ok=false;
    } 
  if(!ok) {
    PRINTFB(FB_Executive,FB_Errors)
      " SettingGet-Error: object \"%s\" not found.\n",object
      ENDFB;
    ok=false;
  } else if(obj) {
    handle = obj->fGetSettingHandle(obj,-1);
    if(handle) set_ptr1 = *handle;
    if(state>=0) {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) 
        set_ptr2 = *handle;
      else {
        PRINTFB(FB_Executive,FB_Errors)
          " SettingGet-Error: object \"%s\" lacks state %d.\n",object,state+1
          ENDFB;
        ok=false;
      }
    }
  }
  if(ok) {
    buffer[0]=0;
  SettingGetTextValue(set_ptr2,set_ptr1,index,buffer);
  result=Py_BuildValue("s",buffer);
  } 
  
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveGetSettingTuple(int index,char *object,int state)
{ /* Assumes blocked Python interpreter */
  PyObject *result = NULL;
  CSetting **handle = NULL;
  CObject *obj=NULL;
  int ok = true;
  PRINTFD(FB_Executive)
    " ExecutiveGetSettingTuple: object %p state %d\n",object,state
    ENDFD;

  if(object[0]==0) /* global */
    result = SettingGetTuple(NULL,NULL,index);
  else {

    if(strlen(object)) {
      obj=ExecutiveFindObjectByName(object);
      if(!obj) 
        ok=false;
    } else ok=false;
    if(!ok) {
      PRINTFB(FB_Executive,FB_Errors)
        " Executive: object not found.\n"
        ENDFB;
    } else {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) 
        result = SettingGetDefinedTuple(*handle,index);      
    }
  }
  if(!ok) {
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return(result);
}
/*========================================================================*/
void ExecutiveSetLastObjectEdited(CObject *o)
{
  CExecutive *I = &Executive;
  I->LastEdited = o;
}
/*========================================================================*/
CObject *ExecutiveGetLastObjectEdited(void)
{
  CExecutive *I = &Executive;
  return(I->LastEdited);
}
/*========================================================================*/
int ExecutiveSaveUndo(char *s1,int state)
{
  int sele1;
  ObjectMoleculeOpRec op1;

  if(state<0) state = SceneGetState();                
  sele1=SelectorIndexByName(s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2=0;
  if(sele1>=0) {
    op1.code = OMOP_SaveUndo;
    op1.i1 = state;
    ExecutiveObjMolSeleOp(sele1,&op1);
  }
  return(op1.i2);
}

/*========================================================================*/
int ExecutiveSetTitle(char *name,int state,char *text)
{
  int result=false;
  ObjectMolecule *obj;
  obj =ExecutiveFindObjectMoleculeByName(name);
  if(!obj) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB;
  } else {
    result = ObjectMoleculeSetStateTitle(obj,state,text);
  }
  SceneDirty();
  return(result);
}
/*========================================================================*/
char *ExecutiveGetTitle(char *name,int state)
{
  char *result = NULL;
  ObjectMolecule *obj;
  obj =ExecutiveFindObjectMoleculeByName(name);
  if(!obj) {
    PRINTFB(FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB;
  } else {
    result = ObjectMoleculeGetStateTitle(obj,state);
  }
  SceneDirty();
  return(result);
}
/*========================================================================*/
void ExecutiveHideSelections(void)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection) {
      if(rec->visible) {
        rec->visible=false;
        SceneDirty();
      }
    }
  }
}
/*========================================================================*/
void ExecutiveRenderSelections(int curState)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  SpecRec *rec1;
  int sele;
  int no_depth;
  float width;

  no_depth = (int)SettingGet(cSetting_selection_overlay);
  width = SettingGet(cSetting_selection_width);

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection) {

      if(rec->visible) {
        sele = SelectorIndexByName(rec->name); /* TODO: speed this up */
        if(sele>=0) {
          rec1 = NULL;
          if(rec->sele_color<0)
            glColor3f(1.0F,0.2F,0.8F);
          else
            glColor3fv(ColorGet(rec->sele_color));
          glPointSize(width);
          if(no_depth)
            glDisable(GL_DEPTH_TEST);
          glBegin(GL_POINTS);
          while(ListIterate(I->Spec,rec1,next)) {
            if(rec1->type==cExecObject) {
              if(rec1->obj->type==cObjectMolecule) {
                ObjectMoleculeRenderSele((ObjectMolecule*)rec1->obj,curState,sele);
              }
            }
          }
          glEnd();
          if(no_depth)
            glEnable(GL_DEPTH_TEST);
        }
      }
    }
  }
}
/*========================================================================*/
int ExecutiveGetDihe(char *s0,char *s1,char *s2,char *s3,float *value,int state)
{
  Vector3f v0,v1,v2,v3;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int ok=true;
  
  if((sele0 = SelectorIndexByName(s0))<0)
    ok = ErrMessage("GetDihedral","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(s1))<0)
    ok = ErrMessage("GetDihedral","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(s2))<0)
    ok = ErrMessage("GetDihedral","Selection 3 invalid.");
  else if((sele3 = SelectorIndexByName(s3))<0)
    ok = ErrMessage("GetDihedral","Selection 4 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(sele0,state,v0))
      ok = ErrMessage("GetDihedral","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(sele1,state,v1))
      ok = ErrMessage("GetDihedral","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(sele2,state,v2))
      ok = ErrMessage("GetDihedral","Selection 3 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(sele3,state,v3))
      ok = ErrMessage("GetDihedral","Selection 4 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    (*value)=rad_to_deg(get_dihedral3f(v0,v1,v2,v3));
  }
  return ok;
}
/*========================================================================*/
int ExecutiveSetDihe(char *s0,char *s1,char *s2,char *s3,float value,int state)
{
  Vector3f v0,v1,v2,v3;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int ok=true;
  int save_state;
  float current;
  float change;

  if((sele0 = SelectorIndexByName(s0))<0)
    ok = ErrMessage("GetDihedral","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(s1))<0)
    ok = ErrMessage("GetDihedral","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(s2))<0)
    ok = ErrMessage("GetDihedral","Selection 3 invalid.");
  else if((sele3 = SelectorIndexByName(s3))<0)
    ok = ErrMessage("GetDihedral","Selection 4 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(sele0,state,v0))
      ok = ErrMessage("GetDihedral","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(sele1,state,v1))
      ok = ErrMessage("GetDihedral","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(sele2,state,v2))
      ok = ErrMessage("GetDihedral","Selection 3 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(sele3,state,v3))
      ok = ErrMessage("GetDihedral","Selection 4 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    current=rad_to_deg(get_dihedral3f(v0,v1,v2,v3));
    change=value-current;
    save_state = SceneGetState();                
    SceneSetFrame(-1,state); /* KLUDGE ALERT!
                             * necessary because the editor 
                             * can only work on the current state...this
                             * needs to be changed.*/
    EditorSelect(s2,s1,NULL,NULL,false,true);
    EditorTorsion(change);
    SceneSetFrame(-1,save_state);
    PRINTFB(FB_Editor,FB_Actions)
      " SetDihedral: adjusted to %5.3f\n",value
      ENDFB;

  }
  return ok;
}
/*========================================================================*/
float ExecutiveGetArea(char *s0,int sta0,int load_b)
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
  sele0 = SelectorIndexByName(s0);
  if(sele0<0) {
    ErrMessage("Area","Invalid selection.");
  } else {
    obj0 = SelectorGetSingleObjectMolecule(sele0);
    if(!(obj0))
      ErrMessage("Area","Selection must be within a single object.");
    else {
      cs = ObjectMoleculeGetCoordSet(obj0,sta0);
      if(!cs)
        ErrMessage("Area","Invalid state.");
      else {
        rep = (RepDot*)RepDotDoNew(cs,cRepDotAreaType);
        if(!rep) 
          ErrMessage("Area","Can't get dot representation.");
        else {

          if(load_b) {
            /* zero out B-values within selection */
            ObjectMoleculeOpRecInit(&op);
            op.code=OMOP_SetB;
            op.f1=0.0;
            op.i1=0;
            ExecutiveObjMolSeleOp(sele0,&op);
          }

          result=0.0;
          
          area=rep->A;
          ati=rep->Atom;
          
          is_member = false;

          for(a=0;a<rep->N;a++) {
            
            if(known_member!=(*ati)) {
              known_member=(*ati);
              ai=obj0->AtomInfo+known_member;
              is_member = SelectorIsMember(ai->selEntry,sele0);
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
char *ExecutiveGetNames(int mode,int enabled_only)
{
  CExecutive *I = &Executive;
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
int ExecutiveGetType(char *name,WordType type)
{
  SpecRec *rec = NULL;
  int ok=true;
  rec = ExecutiveFindSpec(name);
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
void ExecutiveUpdateCmd(char *s0,char *s1,int sta0,int sta1)
{
  int sele0,sele1;

  sele0 = SelectorIndexByName(s0);
  sele1 = SelectorIndexByName(s1);
  if(!(sele0&&sele1)) {
    ErrMessage("Update","One or more invalid input selections.");
  } else {
    SelectorUpdateCmd(sele0,sele1,sta0,sta1);
  }
}
/*========================================================================*/
void ExecutiveRenameObjectAtoms(char *name,int force) 
{
  CExecutive *I = &Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(name);
    if(!os)
      ErrMessage(" Executive","object not found.");
    else if(os->type!=cObjectMolecule) {
      ErrMessage(" Executive","bad object type.");
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
    SceneChanged();
  }
} 

/*========================================================================*/
int  ExecutiveInvert(char *s0,char *s1,int mode)
{
  int i0=-1;
  int i1=-1;
  int sele0,sele1;
  int ok=false;
  ObjectMolecule *obj0,*obj1;

  sele0 = SelectorIndexByName(s0);
  if(sele0<0) {
    ErrMessage("Invert","Please indicate immobile fragments with (lb) and (rb).");
  } else {
    obj0 = SelectorGetSingleObjectMolecule(sele0);
    sele1 = SelectorIndexByName(s1);
    if(sele1>=0) {
      obj1 = SelectorGetSingleObjectMolecule(sele1);
    } else {
      sele1=sele0;
      obj1=obj0;
    }
    i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
    if(obj1)
      i1 = ObjectMoleculeGetAtomIndex(obj1,sele1);
    if(!(obj0&&(obj0==obj1)&&(i0>=0)&&(i1>=0)))
      ErrMessage("Invert","Invalid immobile atoms in (lb) and (rb).");
    else {
      ok = EditorInvert(obj0,sele0,sele1,mode);
    }
  }
  return(ok);
}
/*========================================================================*/
void ExecutiveFuse(char *s0,char *s1,int mode)
{
  int i0=-1;
  int i1=-1;
  int sele0,sele1,sele2;
  ObjectMolecule *obj0,*obj1;
  ObjectMoleculeOpRec op;
  
#define tmp_fuse_sele "tmp_fuse_sele"

  sele0 = SelectorIndexByName(s0);
  if(sele0>=0) {
    sele1 = SelectorIndexByName(s1);
    if(sele1>=0) {
      EditorSetActiveObject(NULL,0);
      obj0 = SelectorGetSingleObjectMolecule(sele0);
      obj1 = SelectorGetSingleObjectMolecule(sele1);
      if(obj0)
        i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
      if(obj1)
        i1 = ObjectMoleculeGetAtomIndex(obj1,sele1);
      if(obj0&&obj1&&(i0>=0)&&(i1>=0)&&(obj0!=obj1)) {
        ObjectMoleculeVerifyChemistry(obj0);
        ObjectMoleculeVerifyChemistry(obj1);
        
        SelectorCreate(tmp_fuse_sele,NULL,obj0,1,NULL);
        sele2=SelectorIndexByName(tmp_fuse_sele);
        if(mode) {
          ObjectMoleculeOpRecInit(&op);
          op.code=OMOP_PrepareFromTemplate;
          op.ai=obj1->AtomInfo+i1;
          op.i1=mode;
          op.i2=0;
          ExecutiveObjMolSeleOp(sele2,&op);
        }
        SelectorDelete(tmp_fuse_sele);

        if((obj0->AtomInfo[i0].protons==1)&&
           (obj1->AtomInfo[i1].protons==1))
          ObjectMoleculeFuse(obj1,i1,obj0,i0,0);
        else if((obj0->AtomInfo[i0].protons!=1)&&
                (obj1->AtomInfo[i1].protons!=1))
          ObjectMoleculeFuse(obj1,i1,obj0,i0,1);
        else 
          ErrMessage("Fuse","Can't fuse between a hydrogen and a non-hydrogen");
      }
    }
  }
}

/*========================================================================*/
void ExecutiveSpheroid(char *name,int average)  /* EXPERIMENTAL */
{
  CExecutive *I = &Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(name);
    if(!os)
      ErrMessage(" Executive","object not found.");
    else if(os->type!=cObjectMolecule) {
      ErrMessage(" Executive","bad object type.");
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
    SceneChanged();
  }
} 
/*========================================================================*/
void ExecutiveRebuildAll(void)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  PRINTFD(FB_Executive)
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
        if(rec->obj->fInvalidate) {
          rec->obj->fInvalidate((CObject*)rec->obj,cRepAll,cRepInvAll,-1);
        }
        break;
      }
    }
  }
  SceneDirty();
}
/*========================================================================*/
void ExecutiveRebuildAllObjectDist(void)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      if(rec->obj->type==cObjectDist) {
        ObjectDistInvalidateRep((ObjectDist*)rec->obj,cRepAll);
      }
    }
  }
  SceneDirty();
}
/*========================================================================*/
void ExecutiveUndo(int dir)
{
  CExecutive *I = &Executive;
  CObject *o;
  ObjectMolecule *obj=NULL,*compObj;
  SpecRec *rec = NULL;

  o = ExecutiveGetLastObjectEdited();
  PRINTFB(FB_Executive,FB_Debugging)
    " ExecutiveUndo: last object %p\n",o
    ENDFB;
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
void ExecutiveSort(char *name)
{
  CExecutive *I = &Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;
  ObjectMoleculeOpRec op;
  int all_obj = false;
  int sele;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(name);
    if(!os) {
      if(!WordMatchExact(cKeywordAll,name,true))
        ErrMessage(" Executive","object not found.");
      else
        all_obj=true;
    } else if(os->type!=cObjectMolecule)
      ErrMessage(" Executive","bad object type.");
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
            sele=SelectorIndexByName(rec->obj->Name);
            if(sele>=0) {
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_INVA;
              op.i1=cRepAll; 
              op.i2=cRepInvRep;
              ExecutiveObjMolSeleOp(sele,&op);
            }
          }
    }
    SceneChanged();
  }
}
/*========================================================================*/
void ExecutiveRemoveAtoms(char *s1)
{
  int sele;
  CExecutive *I=&Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  ObjectMoleculeOpRec op;
  int flag = false;

  sele=SelectorIndexByName(s1);
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
                    PRINTFD(FB_Editor)
                      " ExecutiveRemove-Debug: purging %i of %i atoms in %s\n",
                      op.i1,obj->NAtom,obj->Obj.Name
                      ENDFD;
                    ObjectMoleculePurge(obj);
                    PRINTFB(FB_Editor,FB_Actions)
                      " Remove: eliminated %d atoms in model \"%s\".\n",
                      op.i1,obj->Obj.Name 
                      ENDFB;
                    flag=true;
                  }
					 }
				}
		  }
	 }
  /*  if(!flag) {
      ErrMessage("Remove","no atoms removed.");
      }*/
}
/*========================================================================*/
void ExecutiveAddHydrogens(char *s1)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_AddHydrogens; /* 4 passes completes the job */
    ExecutiveObjMolSeleOp(sele1,&op);    
  }
}
/*========================================================================*/
void ExecutiveFlag(int flag,char *s1,int action,int quiet)
{
  int sele1;
  OrthoLineType buffer;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(s1);
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
    ExecutiveObjMolSeleOp(sele1,&op);    
    if(Feedback(FB_Executive,FB_Actions)) {
      if(!quiet) {
        switch(action) {
        case 0:
          if(op.i3) {
            PRINTF " Flag: flag %d is set in %d of %d atoms.\n", flag, op.i3, op.i4 ENDF;
          } else {
            PRINTF " Flag: flag %d cleared on all atoms.\n", flag ENDF;
          }
          break;
        case 1:
          PRINTF " Flag: flag %d set on %d atoms.\n", flag, op.i3 ENDF;
          break;
        case 2:
          PRINTF " Flag: flag %d cleared on %d atoms.\n", flag, op.i3 ENDF;
          break;
        }
      }
    }
    if((int)SettingGet(cSetting_auto_indicate_flags)) {
      sprintf(buffer,"(flag %d)",flag);
      SelectorCreate(cIndicateSele,buffer,NULL,true,NULL);
      ExecutiveSetObjVisib(cIndicateSele,true);
      SceneDirty();
    }
  }

}
/*========================================================================*/
float ExecutiveOverlap(char *s1,int state1,char *s2,int state2,float adjust)
{
  int sele1,sele2;
  float result=0.0;

  if(state1<0) state1=0;
  if(state2<0) state2=0;
                 
  sele1=SelectorIndexByName(s1);
  sele2=SelectorIndexByName(s2);

  if((sele1>=0)&&(sele2>=0))
    result = SelectorSumVDWOverlap(sele1,state1,sele2,state2,adjust);

  return(result);
}
/*========================================================================*/
void ExecutiveProtect(char *s1,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Protect;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(sele1,&op);    
    if(Feedback(FB_Executive,FB_Actions)) {
      if(op.i2) {
        if(mode) {
          PRINTF " Protect: %d atoms protected from movement.\n",op.i2 ENDF;
        } else {
          PRINTF " Protect: %d atoms deprotected.\n", op.i2 ENDF;
        }
      }
    }
  }
}
/*========================================================================*/
void ExecutiveMask(char *s1,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Mask;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(sele1,&op);
    if(Feedback(FB_Executive,FB_Actions)) {    
      if(op.i2) {
        if(mode) {
          PRINTF " Protect: %d atoms masked (can not be picked).\n",op.i2 ENDF;
        } else {
          PRINTF " Protect: %d atoms unmasked.\n", op.i2 ENDF;
        }
      }
    }
    op.code = OMOP_INVA; /* need to invalidate all pickable representations */
    op.i1 = cRepAll;
    op.i2 = cRepInvPick;
    ExecutiveObjMolSeleOp(sele1,&op);    
  }
}
/*========================================================================*/
int ExecutiveStereo(int flag)
{
  int ok=1;
  int stereo_mode;

  switch(flag) {
  case -1:
    SettingSet(cSetting_stereo_shift,-SettingGet(cSetting_stereo_shift));
    SettingSet(cSetting_stereo_angle,-SettingGet(cSetting_stereo_angle));
    break;
  default:
    
    if(PMGUI) {
      stereo_mode = (int)SettingGet(cSetting_stereo_mode);
      
      switch(stereo_mode) {
      case 1: /* hardware stereo-in-a-window*/
        if(StereoCapable||SceneGetStereo()) {
          SceneSetStereo(flag);
          PSGIStereo(flag);
        } else {
          ok=false;
        }
        break;
      case 2: /* cross-eye stereo*/
      case 3:
        SceneSetStereo(flag);
        break;
      
      }
    }
  }
  return(ok);
}
/*========================================================================*/
void ExecutiveBond(char *s1,char *s2,int order,int add)
{
  int sele1,sele2;
  int cnt;
  CExecutive *I=&Executive;
  SpecRec *rec = NULL;
  int flag = false;

  sele1=SelectorIndexByName(s1);
  sele2=SelectorIndexByName(s2);
  
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
                      PRINTFB(FB_Editor,FB_Actions)
                        " AddBond: %d bonds added to model \"%s\".\n",cnt,rec->obj->Name 
                        ENDFB;
                      flag=true;
                    }
                  } else if(add==2) {
                    cnt = ObjectMoleculeAdjustBonds((ObjectMolecule*)rec->obj,sele1,sele2,1,order);                    
                  }
                  else {
                    cnt = ObjectMoleculeRemoveBonds((ObjectMolecule*)rec->obj,sele1,sele2);
                    if(cnt) {
                      PRINTFB(FB_Editor,FB_Actions)
                        " RemoveBond: %d bonds removed from model \"%s\".\n",
                        cnt,rec->obj->Name 
                        ENDFB;
                      flag=true;
                    }
                  }
                }
            }
        }
      if(!flag) {
        if(add) 
          ErrMessage("AddBond","no bonds added.");
        else
          ErrMessage("RemoveBond","no bonds removed.");          
      }
    }
  } else if(sele1<0) {
    ErrMessage("ExecutiveBond","The first selection contains no atoms.");
  } else if(sele2<0) {
    ErrMessage("ExecutiveBond","The second selection contains no atoms.");
  }
}
/*========================================================================*/
float ExecutiveDist(char *nam,char *s1,char *s2,int mode,float cutoff,
                    int labels,int quiet)
{
  int sele1,sele2;
  ObjectDist *obj;
  CObject *anyObj = NULL;
  float result;
  sele1=SelectorIndexByName(s1);
  if(!WordMatch(s2,"same",true))
    sele2=SelectorIndexByName(s2);
  else {
    sele2 = sele1;
  }
  
  if((sele1>=0)&&(sele2>=0)) {
    anyObj = ExecutiveFindObjectByName(nam);
    if(anyObj)
      if(anyObj->type!=cObjectDist)
        ExecutiveDelete(nam);
    obj = ObjectDistNewFromSele((ObjectDist*)anyObj,sele1,sele2,mode,cutoff,labels,&result);
    if(!obj) {
      ErrMessage("ExecutiveDistance","No such distances found.");
    } else {
      ObjectSetName((CObject*)obj,nam);
      ExecutiveManageObject((CObject*)obj,true,quiet);
      ExecutiveSetRepVisib(nam,cRepLine,1);
      if(!labels)
        ExecutiveSetRepVisib(nam,cRepLabel,0);        
    }
  } else if(sele1<0) {
    ErrMessage("ExecutiveDistance","The first selection contains no atoms.");
  } else if(sele2<0) {
    ErrMessage("ExecutiveDistance","The second selection contains no atoms.");
  }
  return(result);
}
/*========================================================================*/
float ExecutiveDistance(char *s1,char *s2)
{
  int sele1,sele2;
  float dist = -1.0;
  
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  
  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  sele1=SelectorIndexByName(s1);
  op1.i1=0;
  op2.i2=0;
  if(sele1>=0) {
    op1.code = OMOP_SUMC;
    op1.v1[0]=0.0;
    op1.v1[1]=0.0;
    op1.v1[2]=0.0;
    ExecutiveObjMolSeleOp(sele1,&op1);
  } else {
    ErrMessage("ExecutiveDistance","The first selection contains no atoms.");
  }
  
  sele2=SelectorIndexByName(s2);
  op2.i1=0;
  op2.i2=0;
  if(sele2>=0) {
    op2.code = OMOP_SUMC;
    op2.v1[0]=0.0;
    op2.v1[1]=0.0;
    op2.v1[2]=0.0;
    op2.i1=0;
    ExecutiveObjMolSeleOp(sele2,&op2);
  } else {
    ErrMessage("ExecutiveDistance","The second selection contains no atoms.");
  }
  
  if(op1.i1&&op2.i1) {
    scale3f(op1.v1,1.0F/op1.i1,op1.v1);
    scale3f(op2.v1,1.0F/op2.i1,op2.v1);
    dist = (float)diff3f(op1.v1,op2.v1);
    PRINTFB(FB_Executive,FB_Results)
      " Distance: %8.3f [%i atom(s) to %i atom(s)]\n",
      dist,op1.i1,op2.i1
      ENDFB;
  } else {
    ErrMessage("ExecutiveRMS","No atoms selected.");
  }
  return(dist);
}
/*========================================================================*/
char *ExecutiveSeleToPDBStr(char *s1,int state,int conectFlag)
{
  char *result=NULL;
  ObjectMoleculeOpRec op1;
  int sele1,l;
  char end_str[] = "END\n";

  ObjectMoleculeOpRecInit(&op1);
  sele1=SelectorIndexByName(s1);
  op1.charVLA=VLAlloc(char,10000);
  if(state<0) state=SceneGetState();
  if(conectFlag) {
    op1.i2=SelectorGetPDB(&op1.charVLA,sele1,state,conectFlag);
  } else {

    op1.i2 = 0;
    op1.i3 = 0; /* atIndex */
    if(sele1>=0) {
      op1.code = OMOP_PDB1;
      op1.i1 = state;
      ExecutiveObjMolSeleOp(sele1,&op1);
    }
  }
  if(!(int)SettingGet(cSetting_pdb_no_end_record)) { /* terminate with END */
    l=strlen(end_str);
    VLACheck(op1.charVLA,char,op1.i2+l+1);
    strcpy(op1.charVLA+op1.i2,end_str);
    op1.i2+=l+1;
  } else { /* terminate */
    VLACheck(op1.charVLA,char,op1.i2+1);
    op1.charVLA[op1.i2]=0;
    op1.i2++;
  }
  result=Alloc(char,op1.i2);
  memcpy(result,op1.charVLA,op1.i2);
  VLAFreeP(op1.charVLA);
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveSeleToChemPyModel(char *s1,int state)
{
  PyObject *result;
  int sele1;
  sele1=SelectorIndexByName(s1);
  if(state<0) state=0;
  PBlock(); /*   PBlockAndUnlockAPI();*/
  result=SelectorGetChemPyModel(sele1,state);
  if(PyErr_Occurred()) PyErr_Print();
  PUnblock(); /* PLockAPIAndUnblock();*/
  return(result);
}
/*========================================================================*/
void ExecutiveSeleToObject(char *name,char *s1,int source,int target)
{
  int sele1;

  sele1=SelectorIndexByName(s1);

  SelectorCreateObjectMolecule(sele1,name,target,source);
}
/*========================================================================*/
void ExecutiveCopy(char *src,char *dst)
{
  CObject *os;
  ObjectMolecule *oSrc,*oDst;
  SpecRec *rec1 = NULL,*rec2=NULL;
  int a;

  os=ExecutiveFindObjectByName(src);
  if(!os)
    ErrMessage(" Executive","object not found.");
  else if(os->type!=cObjectMolecule)
    ErrMessage(" Executive","bad object type.");
  else 
    {
      oSrc =(ObjectMolecule*)os;
      oDst = ObjectMoleculeCopy(oSrc);
      if(oDst) {
        strcpy(oDst->Obj.Name,dst);
        ExecutiveManageObject((CObject*)oDst,true,false);
        rec1=ExecutiveFindSpec(oSrc->Obj.Name);
        rec2=ExecutiveFindSpec(oDst->Obj.Name);
        if(rec1&&rec2) {
          for(a=0;a<cRepCnt;a++)
            rec2->repOn[a]=rec1->repOn[a];
        }
        
        PRINTFB(FB_Executive,FB_Actions)
          " Executive: object %s created.\n",oDst->Obj.Name 
          ENDFB;
      }
    }
  SceneChanged();
}

/*========================================================================*/
void ExecutiveOrient(char *sele,Matrix33d mi,int state)
{
  double egval[3],egvali[3];
  double evect[3][3];
  float m[4][4],mt[4][4];
  float t[3];

  int a,b;

  if(!MatrixEigensolve33d((double*)mi,egval,egvali,(double*)evect)) {

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

    SceneSetMatrix(m[0]); /* load matrix */

    /* there must  be a more elegant to get the PC on X and the SC
     * on Y then what is shown below, but I couldn't get it to work.
     * I tried swapping the eigen-columns around but either that is 
     * a bogus approach (?) or my code was buggy.  Hence the following...*/

    if((egval[0]<egval[2])&&(egval[2]<egval[1])) { /* X < Z < Y */
      SceneRotate(90,1,0,0); /*1<-->2*/
    } else if((egval[1]<egval[0])&&(egval[0]<egval[2])) { /* Y < X < Z */
      SceneRotate(90,0,0,1); /*0<-->1*/
    } else if((egval[1]<egval[2])&&(egval[2]<egval[0])) { /* Y < Z < X */
      SceneRotate(90,0,1,0); /*1<-->2*/
      SceneRotate(90,0,0,1); /*0<-->1*/
    } else if((egval[2]<egval[1])&&(egval[1]<egval[0])) { /* Z < Y < X */
      SceneRotate(90,0,1,0); /*0<-->2*/
    } else if((egval[2]<egval[0])&&(egval[0]<egval[1])) { /* Z < X < Y */
      SceneRotate(90,0,1,0); /*0<-->2*/
      SceneRotate(90,1,0,0); /*0<-->1*/
    }
    /* X < Y < Z  - do nothing - that's what we want */

    ExecutiveWindowZoom(sele,0.0,state,0);

  }
}
/*========================================================================*/
void ExecutiveLabel(char *s1,char *expr,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int cnt;

  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_LABL;
    op1.s1 = expr;
    op1.i1 = 0;
    ExecutiveObjMolSeleOp(sele1,&op1);
    cnt = op1.i1;
    op1.code=OMOP_VISI;
    op1.i1=cRepLabel;
    op1.i2=1;
    ExecutiveObjMolSeleOp(sele1,&op1);
    op1.code = OMOP_INVA;
    op1.i1=cRepLabel; 
    op1.i2=cRepInvVisib;
    ExecutiveObjMolSeleOp(sele1,&op1);

    if(!quiet) {
      PRINTFB(FB_Executive,FB_Actions)
        " Label: labelled %i atoms.\n",cnt
        ENDFB;
    }
  } else {
    PRINTFB(FB_Executive,FB_Warnings)
      " Label: no atoms selections.\n"
      ENDFB;
  }
}
/*========================================================================*/
int ExecutiveIterate(char *s1,char *expr,int read_only,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRecInit(&op1);
  op1.i1=0;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    op1.code = OMOP_ALTR;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = read_only;
    ExecutiveObjMolSeleOp(sele1,&op1);
    if(!quiet) {
      if(!read_only) {
        PRINTFB(FB_Executive,FB_Actions)
          " Alter: modified %i atoms.\n",op1.i1
          ENDFB;
      } else {
        PRINTFB(FB_Executive,FB_Actions)
          " Iterate: iterated over %i atoms.\n",op1.i1
          ENDFB;
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(FB_Executive,FB_Warnings)
        "ExecutiveIterate: No atoms selected.\n"
        ENDFB;
    }
  }
  return(op1.i1);
}
/*========================================================================*/
void ExecutiveIterateState(int state,char *s1,char *expr,int read_only,
                           int atomic_props,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_AlterState;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = state;
    op1.i3 = read_only;
    op1.i4 = atomic_props;
    ExecutiveObjMolSeleOp(sele1,&op1);
    if(!quiet) {
      if(!read_only) {
        PRINTFB(FB_Executive,FB_Actions)
          " AlterState: modified %i atom states.\n",op1.i1
          ENDFB;
      } else {
        PRINTFB(FB_Executive,FB_Actions)
        " IterateState: iterated over %i atom states.\n",op1.i1
          ENDFB;
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(FB_Executive,FB_Warnings)
        "ExecutiveIterateState: No atoms selected.\n"
        ENDFB;
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
float ExecutiveRMS(char *s1,char *s2,int mode,float refine,int max_cyc,
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
            
  sele1=SelectorIndexByName(s1);

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
    ExecutiveObjMolSeleOp(sele1,&op1);
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
  
  sele2=SelectorIndexByName(s2);
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
    ExecutiveObjMolSeleOp(sele2,&op2);
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
      ErrMessage("ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(!SelectorGetSingleObjectMolecule(sele1)) {
        if(mode!=2) {
          PRINTFB(FB_Executive,FB_Warnings)
            "Executive-Warning: Mobile selection spans more than one object.\n"
            ENDFB;
        } else {
          PRINTFB(FB_Executive,FB_Errors)
            "Executive-Error: Mobile selection spans more than one object. Aborting.\n"
            ENDFB;
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
            PRINTFB(FB_Executive,FB_Warnings) 
              "Executive-Warning: Ordering requested but not well defined.\n"
               ENDFB;
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
              UtilSortInPlace(vert,op1.nvv1,sizeof(FitVertexRec),(UtilOrderFn*)fVertexOrdered);
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
              UtilSortInPlace(vert,op2.nvv1,sizeof(FitVertexRec),(UtilOrderFn*)fVertexOrdered);
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
        rms = MatrixFitRMS(op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
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
                PRINTFB(FB_Executive,FB_Actions)
                  " ExecutiveRMS: %d atoms rejected during cycle %d (RMS=%0.2f).\n",op1.nvv1-op2.nvv1,b,rms
                  ENDFB;
              }
              op1.nvv1 = op2.nvv1;
              FreeP(flag);
              if(op1.nvv1) 
                rms = MatrixFitRMS(op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);            
              else
                break;
            }
          }
        }
      }
      else
        rms = MatrixGetRMS(op1.nvv1,op1.vv1,op2.vv1,NULL);

      if(!op1.nvv1) {
        PRINTFB(FB_Executive,FB_Results) 
          " Executive: Error -- no atoms left after refinement!\n"
          ENDFB;
        ok=false;
      }

      if(ok) {
        if(!quiet) {
          PRINTFB(FB_Executive,FB_Results) 
            " Executive: RMS = %8.3f (%d to %d atoms)\n", rms,op1.nvv1,op2.nvv1 
            ENDFB;
        }
        if(oname) 
          if(oname[0]) {
            CGO *cgo = NULL;
            ObjectCGO *ocgo;
            int auto_save;

            cgo=CGONew();
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
            ocgo = ObjectCGOFromCGO(NULL,cgo,0);
            ocgo->Obj.Color = ColorGetIndex("yellow");
            ObjectSetName((CObject*)ocgo,oname);
            ExecutiveDelete(oname);
            auto_save = (int)SettingGet(cSetting_auto_zoom);
            SettingSet(cSetting_auto_zoom,0);
            ExecutiveManageObject((CObject*)ocgo,true,false);
            SettingSet(cSetting_auto_zoom,(float)auto_save);            
            SceneDirty();
          }
        if(mode==2) {
          if(ok) {
            op2.code = OMOP_TTTF;
            ExecutiveObjMolSeleOp(sele1,&op2);
          }
        }
      }
    } else {
      ErrMessage("ExecutiveRMS","No atoms selected.");
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
int *ExecutiveIdentify(char *s1,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  int *result = NULL;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_Identify;
    op2.i1=0;
    op2.i1VLA=VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(sele1,&op2);
    result = op2.i1VLA;
    VLASize(result,int,op2.i1);
  } 
  return(result);
}
/*========================================================================*/
int ExecutiveIdentifyObjects(char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_IdentifyObjects;
    op2.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op2.i1VLA=VLAlloc(int,1000);
    op2.i1=0;
    ExecutiveObjMolSeleOp(sele1,&op2);
    VLASize(op2.i1VLA,int,op2.i1);
    VLASize(op2.obj1VLA,ObjectMolecule*,op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } 
  return(op2.i1);
}
/*========================================================================*/
int ExecutiveIndex(char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_Index;
    op2.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op2.i1VLA=VLAlloc(int,1000);
    op2.i1=0;
    ExecutiveObjMolSeleOp(sele1,&op2);
    VLASize(op2.i1VLA,int,op2.i1);
    VLASize(op2.obj1VLA,ObjectMolecule*,op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } 
  return(op2.i1);
}
/*========================================================================*/
float *ExecutiveRMSStates(char *s1,int target,int mode,int quiet)
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
  sele1=SelectorIndexByName(s1);
  
  if(!SelectorGetSingleObjectMolecule(sele1)) {
    if(mode!=2) {
      PRINTFB(FB_Executive,FB_Warnings)
        "Executive-Warning: Mobile selection spans more than one object.\n"
        ENDFB;
    } else {
      PRINTFB(FB_Executive,FB_Errors)
        "Executive-Error: Mobile selection spans more than one object. Aborting.\n\n"
        ENDFB;
      ok=false;
    }
  }

  if(ok&&sele1>=0) {
    op1.code = OMOP_SVRT;
    op1.nvv1=0;
    op1.i1=target;
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    op1.i1VLA = VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(sele1,&op1);

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
    ExecutiveObjMolSeleOp(sele1,&op2);
    result=op2.f1VLA;
    VLAFreeP(op1.vv1);
    VLAFreeP(op1.i1VLA);
    VLAFreeP(op2.vv1);
  } 
  return(result);
}
/*========================================================================*/
float ExecutiveRMSPairs(WordType *sele,int pairs,int mode)
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
    sele1=SelectorIndexByName(sele[c]);
    if(sele1>=0) ExecutiveObjMolSeleOp(sele1,&op1);
    strcat(combi,sele[c]);
    if(a<(pairs-1)) strcat(combi," or ");
    c++;
    sele2=SelectorIndexByName(sele[c]);
    if(sele2>=0) ExecutiveObjMolSeleOp(sele2,&op2);
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
      ErrMessage("ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(mode!=0)
        rms = MatrixFitRMS(op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
      else
        rms = MatrixGetRMS(op1.nvv1,op1.vv1,op2.vv1,NULL);
      PRINTFB(FB_Executive,FB_Results) 
        " ExecutiveRMS: RMS = %8.3f (%d to %d atoms)\n",
        rms,op1.nvv1,op2.nvv1
        ENDFB;
    
      op2.code = OMOP_TTTF;
      SelectorGetTmp(combi,s1);
      sele1=SelectorIndexByName(s1);
      ExecutiveObjMolSeleOp(sele1,&op2);
      SelectorFreeTmp(s1);
    } else {
      ErrMessage("ExecutiveRMS","No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  return(rms);
}
/*========================================================================*/
void ExecutiveUpdateObjectSelection(struct CObject *obj)
{
  if(obj->type==cObjectMolecule) {
    SelectorUpdateObjectSele((ObjectMolecule*)obj);  
  }
}
/*========================================================================*/
int ExecutiveReset(int cmd,char *name)
{
  int ok=true;
  CObject *obj;
  if(!name[0]) {
    SceneResetMatrix();
    ExecutiveWindowZoom(cKeywordAll,0.0,-1,0); /* reset does all states */
  } else {
    obj = ExecutiveFindObjectByName(name);
    if(!obj)
      ok=false;
    else
      ObjectResetTTT(obj);
  }
  return(ok);
}
/*========================================================================*/
void ExecutiveDrawNow(void) 
{
  PRINTFD(FB_Executive)
    " ExecutiveDrawNow: entered.\n"
    ENDFD;

  if(!SettingGet(cSetting_suspend_updates)) {

    if(PMGUI) {
      glMatrixMode(GL_MODELVIEW);
      /*  glClear( GL_DEPTH_BUFFER_BIT);*/
    }

    SceneUpdate();
    
    OrthoDoDraw();
    
    MainSwapBuffers();
  }

  PRINTFD(FB_Executive)
    " ExecutiveDrawNow: leaving.\n"
    ENDFD;
}
/*========================================================================*/
int ExecutiveCountStates(char *s1)
{
  CExecutive *I = &Executive;
  int sele1;
  int result=0;
  int n_frame;
  SpecRec *rec = NULL;
  
  if(s1)
    if(WordMatch(cKeywordAll,s1,true))
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
  sele1=SelectorIndexByName(s1);
    if(sele1>=0) {
      SelectorUpdateTable();
      result = SelectorGetSeleNCSet(sele1);
    }
  }
  return(result);
}
/*========================================================================*/
void ExecutiveRay(int width,int height,int mode,float angle,float shift,int quiet)
{
  SceneRay(width,height,mode,NULL,NULL,angle,shift,quiet);
}
/*========================================================================*/
int  ExecutiveSetSetting(int index,PyObject *tuple,char *sele,
                         int state,int quiet,int updates)
{
  CExecutive *I=&Executive;
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

  PRINTFD(FB_Executive)
    " ExecutiveSetSetting: entered. sele \"%s\"\n",sele
    ENDFD;
  unblock = PAutoBlock();
  if(sele[0]==0) { 
    ok = SettingSetTuple(NULL,index,tuple);
    if(ok) {
      if(!quiet) {
        if(Feedback(FB_Setting,FB_Actions)) {
          SettingGetTextValue(NULL,NULL,index,value);
          SettingGetName(index,name);
          PRINTF
            " Setting: %s set to %s.\n",name,value
            ENDF;
        }
      }
      if(updates) 
        SettingGenerateSideEffects(index,sele,state);
    }
  } 
  else if(!strcmp(cKeywordAll,sele)) { /* all objects setting */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject) {
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(handle);
              ok = SettingSetTuple(*handle,index,tuple);
              nObj++;
            }
          }
        }
        if(nObj) {
          if(updates) 
            SettingGenerateSideEffects(index,sele,state);
        }
        if(Feedback(FB_Setting,FB_Actions)) {
          if(nObj&&handle) {
            SettingGetTextValue(*handle,NULL,index,value);
            SettingGetName(index,name);
            if(!quiet) {
              if(state<0) {
                PRINTF
                  " Setting: %s set to %s in %d objects.\n",name,value,nObj
                  ENDF;
              } else {
                PRINTF
                  " Setting: %s set to %s in %d objects, state %d.\n",
                  name,value,nObj,state+1
                  ENDF;
              }
            }
          }
        }
      }
  } else { /* based on a selection/object name */
    sele1=SelectorIndexByName(sele);
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
                  SettingCheckHandle(handle);
                  ok = SettingSetTuple(*handle,index,tuple);
                  if(ok) {
                    if(updates) 
                      SettingGenerateSideEffects(index,sele,state);
                    if(!quiet) {
                      if(state<0) { /* object-specific */
                        if(Feedback(FB_Setting,FB_Actions)) {
                          SettingGetTextValue(*handle,NULL,index,value);
                          SettingGetName(index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\".\n",
                            name,value,rec->obj->Name
                            ENDF;
                        }
                      } else { /* state-specific */
                        if(Feedback(FB_Setting,FB_Actions)) {
                          SettingGetTextValue(*handle,NULL,index,value);
                          SettingGetName(index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\", state %d.\n",
                            name,value,rec->obj->Name,state+1
                            ENDF;
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
                SettingCheckHandle(handle);
                ok = SettingSetTuple(*handle,index,tuple);
                if(ok) {
                  if(updates)
                    SettingGenerateSideEffects(index,sele,state);
                  if(!quiet) {
                    if(state<0) { /* object-specific */
                      if(Feedback(FB_Setting,FB_Actions)) {
                        SettingGetTextValue(*handle,NULL,index,value);
                        SettingGetName(index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\".\n",
                          name,value,rec->obj->Name
                          ENDF;
                      }
                    } else { /* state-specific */
                      if(Feedback(FB_Setting,FB_Actions)) {
                        SettingGetTextValue(*handle,NULL,index,value);
                        SettingGetName(index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\", state %d.\n",
                          name,value,rec->obj->Name,state+1
                          ENDF;
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
int  ExecutiveUnsetSetting(int index,char *sele,
                         int state,int quiet,int updates)
{
  CExecutive *I=&Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  CSetting **handle=NULL;
  SettingName name;
  int nObj=0;
  int unblock;
  int ok =true;

  PRINTFD(FB_Executive)
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
              SettingCheckHandle(handle);
              ok = SettingUnset(*handle,index);
              nObj++;
            }
          }
        }
        if(nObj) {
          if(updates) 
            SettingGenerateSideEffects(index,sele,state);
        }
        if(Feedback(FB_Setting,FB_Actions)) {
          if(nObj&&handle) {
            SettingGetName(index,name);
            if(!quiet) {
              if(state<0) {
                PRINTF
                  " Setting: %s unset in %d objects.\n",name,nObj
                  ENDF;
              } else {
                PRINTF
                  " Setting: %s unset in %d objects, state %d.\n",
                  name,nObj,state+1
                  ENDF;
              }
            }
          }
        }
      }
  } else { /* based on a selection/object name */
    sele1=SelectorIndexByName(sele);
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
                  SettingCheckHandle(handle);
                  ok = SettingUnset(*handle,index);
                  if(ok) {
                    if(updates) 
                      SettingGenerateSideEffects(index,sele,state);
                    if(!quiet) {
                      if(state<0) { /* object-specific */
                        if(Feedback(FB_Setting,FB_Actions)) {
                          SettingGetName(index,name);
                          PRINTF
                            " Setting: %s unset in object \"%s\".\n",
                            name,rec->obj->Name
                            ENDF;
                        }
                      } else { /* state-specific */
                        if(Feedback(FB_Setting,FB_Actions)) {
                          SettingGetName(index,name);
                          PRINTF
                            " Setting: %s unset in object \"%s\", state %d.\n",
                            name,rec->obj->Name,state+1
                            ENDF;
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
                SettingCheckHandle(handle);
                ok = SettingUnset(*handle,index);
                if(ok) {
                  if(updates)
                    SettingGenerateSideEffects(index,sele,state);
                  if(!quiet) {
                    if(state<0) { /* object-specific */
                      if(Feedback(FB_Setting,FB_Actions)) {
                        SettingGetName(index,name);
                        PRINTF
                          " Setting: %s unset in object \"%s\".\n",
                          name,rec->obj->Name
                          ENDF;
                      }
                    } else { /* state-specific */
                      if(Feedback(FB_Setting,FB_Actions)) {
                        SettingGetName(index,name);
                        PRINTF
                          " Setting: %s unset in object \"%s\", state %d.\n",
                          name,rec->obj->Name,state+1
                          ENDF;
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
int ExecutiveColor(char *name,char *color,int flags,int quiet)
{
  CExecutive *I = &Executive;
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
  col_ind = ColorGetIndex(color);
  if(col_ind==-1) {
    ErrMessage("Color","Unknown color.");
  } else {
    best_match = ExecutiveFindBestNameMatch(name);
    /* per atom */
    if(!(flags&0x1)) {
      sele=SelectorIndexByName(name);
      if(sele>=0) {
        ok=true; 
        ObjectMoleculeOpRecInit(&op);
        op.code = OMOP_COLR;
        op.i1= col_ind;
        op.i2= 0;
        ExecutiveObjMolSeleOp(sele,&op);
        n_atm = op.i2;
        op.code=OMOP_INVA;
        op.i1=cRepAll; 
        op.i2=cRepInvColor;
        ExecutiveObjMolSeleOp(sele,&op);
      }
    }
    /* per object */
    if(strcmp(name,cKeywordAll)) {
      rec=ExecutiveFindSpec(name);
      if(rec) {
        if(rec->type==cExecObject) {
          rec->obj->Color=col_ind;
          n_obj++;
          ok=true;
          SceneDirty();
        }
      } 
    } else {
      rec=NULL;
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          rec->obj->Color=col_ind;
          n_obj++;
          ok=true;
          SceneDirty();
        }
      }
    }
    if(n_obj||n_atm) {
      if(n_obj<2) objs[0]=0;
      if(n_atm<2) atms[0]=0;
      if(!quiet) {
        if(n_obj&&n_atm) {
          PRINTFB(FB_Executive,FB_Actions)
            " Executive: Colored %d atom%s and %d object%s.\n",n_atm,atms,n_obj,objs
            ENDFB;
        } else if (n_obj) {
          PRINTFB(FB_Executive,FB_Actions)
            " Executive: Colored %d object%s.\n",n_obj,objs
            ENDFB;
        } else {
          PRINTFB(FB_Executive,FB_Actions)
            " Executive: Colored %d atom%s.\n",n_atm,atms
            ENDFB;
        }
      }
    }
  }
  return(ok);
}
/*========================================================================*/
char *ExecutiveFindBestNameMatch(char *name)
{
  char *result;
  CExecutive *I = &Executive;
  SpecRec *rec=NULL,*best_rec = NULL;
  int best;
  int wm;

  best = 0;
  result = name;

  while(ListIterate(I->Spec,rec,next)) {
    wm = WordMatch(name,rec->name,true);
    if(wm<0) {
      best_rec = rec;
      break;
    } else if ((wm>0)&&(best<wm)) {
      best_rec=rec;
      best = wm;
    }
  }
  if(best_rec)
    result=best_rec->name;
  return(result);
}
/*========================================================================*/
SpecRec *ExecutiveFindSpec(char *name)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next)) {
	 if(strcmp(rec->name,name)==0) 
		break;
  }
  return(rec);
}
/*========================================================================*/
void ExecutiveObjMolSeleOp(int sele,ObjectMoleculeOpRec *op) 
{
  CExecutive *I=&Executive;
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
int ExecutiveGetCameraExtent(char *name,float *mn,float *mx,int transformed,int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int flag = false;

  if(state==-2) state=SceneGetState();

  PRINTFD(FB_Executive)
    " ExecutiveGetCameraExtent: name %s state %d\n",name,state
    ENDFD;
  
  sele=SelectorIndexByName(name);

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
    op.mat1=SceneGetMatrix();

    ExecutiveObjMolSeleOp(sele,&op);

    PRINTFD(FB_Executive)
      " ExecutiveGetCameraExtent: minmax over %d vertices\n",op.i1
      ENDFD;
    if(op.i1)
      flag = true;
  }
  copy3f(op.v1,mn);
  copy3f(op.v2,mx);
  
  PRINTFD(FB_Executive)
    " ExecutiveGetCameraExtent: returning %d\n",flag
    ENDFD;

  return(flag);  
}

/*========================================================================*/
int ExecutiveGetExtent(char *name,float *mn,float *mx,int transformed,int state,int weighted)
{
  int sele;
  ObjectMoleculeOpRec op,op2;
  CExecutive *I=&Executive;
  CObject *obj;
  int flag = false;
  SpecRec *rec = NULL;
  int all_flag = false;
  float f1,f2,fmx;
  int a;

  if(WordMatch(cKeywordCenter,name,1)<0) {
    SceneGetPos(mn);
    copy3f(mn,mx);
    return 1;
  }
  if(WordMatch(cKeywordOrigin,name,1)<0) {
    SceneOriginGet(mn);
    copy3f(mn,mx);
    return 1;
  }
  if(state==-2) state=SceneGetState();

  PRINTFD(FB_Executive)
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
  
  if(WordMatch(cKeywordAll,name,true)<0) {
    all_flag=true;
  }
  sele=SelectorIndexByName(name);

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
    ExecutiveObjMolSeleOp(sele,&op);

    PRINTFD(FB_Executive)
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
      ExecutiveObjMolSeleOp(sele,&op2);
      if(op2.i1) {
        op2.v1[0]/=op2.i1;
        op2.v1[1]/=op2.i1;
        op2.v1[2]/=op2.i1;
      }
    }
  } else {
    obj = ExecutiveFindObjectByName(name);
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

          PRINTFD(FB_Executive)
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
  PRINTFD(FB_Executive)
    " ExecutiveGetExtent: returning %d\n",flag
    ENDFD;

  return(flag);  
}
/*========================================================================*/
int ExecutiveGetMaxDistance(char *name,float *pos,float *dev,int transformed,int state)
{
  int sele;
  ObjectMoleculeOpRec op,op2;
  CExecutive *I=&Executive;
  CObject *obj;
  int flag = false;
  SpecRec *rec = NULL;
  int all_flag = false;
  float f1,fmx=0.0F;

  if(state==-2) state=SceneGetState();

  PRINTFD(FB_Executive)
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
  
  if(WordMatch(cKeywordAll,name,true)<0) {
      all_flag=true;
  }
  sele=SelectorIndexByName(name);

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
    ExecutiveObjMolSeleOp(sele,&op);
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
    obj = ExecutiveFindObjectByName(name);
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
int ExecutiveWindowZoom(char *name,float buffer,int state,int inclusive)
{
  float center[3],radius;
  float mn[3],mx[3],df[3];
  int sele0;
  int ok=true;
  PRINTFD(FB_Executive)
    " ExecutiveWindowZoom-DEBUG: entered\n"
    ENDFD;
  if(ExecutiveGetExtent(name,mn,mx,true,state,true)) {
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
      if(!ExecutiveGetMaxDistance(name,center,&radius,true,state))
        radius=0.0;
      radius+=buffer;
    } else {
      radius = df[0];
      if(radius<df[1]) radius=df[1];
      if(radius<df[2]) radius=df[2];
      radius=radius/2.0F;
    }
    if(radius<MAX_VDW) radius=MAX_VDW;
    PRINTFD(FB_Executive)
      " ExecutiveWindowZoom: zooming with radius %8.3f...state %d\n",radius,state
      ENDFD;
    PRINTFD(FB_Executive)
      " ExecutiveWindowZoom: on center %8.3f %8.3f %8.3f...\n",center[0],
      center[1],center[2]
      ENDFD;
    SceneOriginSet(center,false);
    SceneWindowSphere(center,radius);
    SceneDirty();
  } else {

    sele0 = SelectorIndexByName(name);
    if(sele0>0) { /* any valid selection except "all" */
      ErrMessage("ExecutiveWindowZoom","selection doesn't specify any coordinates.");
      ok=false;
    } else if(ExecutiveValidName(name)) {
      PRINTFD(FB_Executive)
        " ExecutiveWindowZoom-DEBUG: name valid, but no extents -- using default view\n"
        ENDFD;
      SceneSetDefaultView();
      SceneDirty();
    } else {
      ErrMessage("ExecutiveWindowZoom","selection or object unknown.");
      ok=false;
    }
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveCenter(char *name,int state,int origin)
{
  float center[3];
  float mn[3],mx[3],df[3];
  int sele0;
  int ok=true;
  if(ExecutiveGetExtent(name,mn,mx,true,state,true)) {
    subtract3f(mx,mn,df);
    average3f(mn,mx,center);
    PRINTFD(FB_Executive)
      " ExecutiveCenter: centering state %d\n",state
      ENDFD;
    PRINTFD(FB_Executive)
      " ExecutiveCenter: on center %8.3f %8.3f %8.3f...\n",center[0],
      center[1],center[2]
      ENDFD;
    if(origin) 
      SceneOriginSet(center,false);
    SceneRelocate(center);
    SceneDirty();
  } else {
    sele0 = SelectorIndexByName(name);
    if(sele0>=0) { /* any valid selection except "all" */
      ErrMessage("ExecutiveCenter","selection doesn't specify any coordinates.");
      ok=false;
    } else if(ExecutiveValidName(name)) {
      SceneSetDefaultView();
      SceneDirty();
    } else {
      ErrMessage("ExecutiveCenter","selection or object unknown.");
      ok=false;
    }
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveOrigin(char *name,int preserve,char *oname,float *pos,int state)
{
  float center[3];
  float mn[3],mx[3];
  int ok=true;
  CObject *obj = NULL;
  if(oname[0]) {
    obj = ExecutiveFindObjectByName(oname);
    if(!obj)
      ok=false;
  }
  if(ok) {
    if(name[0]) {
      ok = ExecutiveGetExtent(name,mn,mx,(oname[0]==0),state,true);
      if(ok) 
        average3f(mn,mx,center);
    } else {
      copy3f(pos,center)
        }
  }
  if(ok) {
    if(obj) {
      ObjectSetTTTOrigin(obj,center);
      PRINTFB(FB_Executive,FB_Blather)
        " ExecutiveCenter: origin for %s set to %8.3f %8.3f %8.3f\n",
        oname,center[0],center[1],center[2]
        ENDFB;
    } else {
      PRINTFB(FB_Executive,FB_Blather)
        " ExecutiveCenter: scene origin set to %8.3f %8.3f %8.3f\n",
        center[0],center[1],center[2]
        ENDFB;
      SceneOriginSet(center,preserve);
    }
    SceneDirty();
  } else
    ok=false;
  return(ok);
}
/*========================================================================*/
int ExecutiveGetMoment(char *name,Matrix33d mi,int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int a,b;
  int c=0;

  if(state==-2) state=SceneGetState();

  for(a=0;a<3;a++)
	 {
		for(b=0;b<3;b++)
		  mi[a][b]=0.0;
		mi[a][a]=1.0;
	 }
  
  sele=SelectorIndexByName(name);
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
	 
	 ExecutiveObjMolSeleOp(sele,&op);
	 
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
		ExecutiveObjMolSeleOp(sele,&op);			 
		for(a=0;a<3;a++)
		  for(b=0;b<3;b++)
			 mi[a][b]=op.d[a][b];
	 }
  } 
  return(c);
}
/*========================================================================*/
void ExecutiveSetObjVisib(char *name,int state)
{
  CExecutive *I = &Executive;
  SpecRec *tRec;

  PRINTFD(FB_Executive)
    " ExecutiveSetObjVisib: entered.\n"
    ENDFD;

  if(strcmp(name,cKeywordAll)==0) {
    tRec=NULL;
    while(ListIterate(I->Spec,tRec,next)) {
      if(state!=tRec->visible) {
        if(tRec->type==cExecObject) {
          if(tRec->visible)
            SceneObjectDel(tRec->obj);				
          else {
            SceneObjectAdd(tRec->obj);
          }
        }
        if((tRec->type!=cExecSelection)||(!state)) /* hide all selections, but show all */
          tRec->visible=!tRec->visible;
      }
    }
  } else {
    tRec = ExecutiveFindSpec(name);
    if(tRec) {
      if(tRec->type==cExecObject) {
        if(tRec->visible!=state)
          {
            if(tRec->visible)
              SceneObjectDel(tRec->obj);				
            else {
              SceneObjectAdd(tRec->obj);
            }
            tRec->visible=!tRec->visible;
          }
      }
      else if(tRec->type==cExecSelection) {
        if(tRec->visible!=state) {
          tRec->visible=!tRec->visible;
          SceneChanged();
        }
      }
    }
  }
  PRINTFD(FB_Executive)
    " ExecutiveSetObjVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void ExecutiveFullScreen(int flag)
{
  if(PMGUI) {
    SettingSet(cSetting_full_screen,(float)flag);
    if(flag) {
      p_glutFullScreen();
    } else {
      p_glutReshapeWindow(640+(int)SettingGet(cSetting_internal_gui_width),
                          480+cOrthoBottomSceneMargin);
    }
  }
}
/*========================================================================*/
void ExecutiveSetAllVisib(int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int rep;
  int sele;
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;

  PRINTFD(FB_Executive)
    " ExecutiveSetAllVisib: entered.\n"
    ENDFD;


  while(ListIterate(I->Spec,rec,next)) {
	 if(rec->type==cExecObject)
		{
        switch(rec->obj->type) {
        case cObjectMolecule:
          obj=(ObjectMolecule*)rec->obj;
          sele = SelectorIndexByName(obj->Obj.Name);
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
          SceneDirty();
          break;
        }
		}
  }
  PRINTFD(FB_Executive)
    " ExecutiveSetAllVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void ExecutiveSetRepVisib(char *name,int rep,int state)
{
  int sele;
  int a;
  int handled = false;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;

  PRINTFD(FB_Executive)
    " ExecutiveSetRepVisib: entered.\n"
    ENDFD;

  tRec = ExecutiveFindSpec(name);
  if((!tRec)&&(!strcmp(name,cKeywordAll))) {
    ExecutiveSetAllRepVisib(name,rep,state);
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
            tRec->obj->fInvalidate(tRec->obj,rep,cRepInvVisib,state);
        } else {
          for(a=0;a<cRepCnt;a++) {
            tRec->repOn[a]=state; 
            ObjectSetRepVis(tRec->obj,a,state);
            if(tRec->obj->fInvalidate)
              tRec->obj->fInvalidate(tRec->obj,a,cRepInvVisib,state);
          }
        }
        SceneChanged();
        break;
      }
    if(!handled)
      switch(tRec->type) {
      case cExecSelection:
      case cExecObject:
        sele=SelectorIndexByName(name);
        if(sele>=0) {
          ObjectMoleculeOpRecInit(&op);

          op.code=OMOP_VISI;
          op.i1=rep;
          op.i2=state;
          ExecutiveObjMolSeleOp(sele,&op);
          op.code=OMOP_INVA;
          op.i2=cRepInvVisib;
          ExecutiveObjMolSeleOp(sele,&op);
        }
        break;
      }
  }
  PRINTFD(FB_Executive)
    " ExecutiveSetRepVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void ExecutiveSetAllRepVisib(char *name,int rep,int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int sele;
  int a;
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  PRINTFD(FB_Executive)
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
            sele = SelectorIndexByName(obj->Obj.Name);
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
            SceneDirty();
            break;
          }
        }
		}
  }
  PRINTFD(FB_Executive)
    " ExecutiveSetAllRepVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void ExecutiveInvalidateRep(char *name,int rep,int level)
{
  int sele = -1;
  ObjectMoleculeOpRec op;
  int all_flag=false;
  PRINTFD(FB_Executive)
    "ExecInvRep-Debug: %s %d %d\n",name,rep,level
    ENDFD;
  if(WordMatch(cKeywordAll,name,true)<0) {
    all_flag=true;
  }
  sele=SelectorIndexByName(name);
  if(sele>=0) {
    ObjectMoleculeOpRecInit(&op);
	 op.code = OMOP_INVA;
	 op.i1=rep;
	 op.i2=level;
	 ExecutiveObjMolSeleOp(sele,&op);
  }
}


/*========================================================================*/
CObject *ExecutiveFindObjectByName(char *name)
{
  CExecutive *I = &Executive;
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
ObjectMap *ExecutiveFindObjectMapByName(char *name)
{
  CObject *obj;
  
  obj = ExecutiveFindObjectByName(name);
  if(obj)
    if(obj->type!=cObjectMap)
      obj=NULL;
  return((ObjectMap*)obj);
}

/*========================================================================*/
ObjectMolecule *ExecutiveFindObjectMoleculeByName(char *name)
{
  CObject *obj;
  
  obj = ExecutiveFindObjectByName(name);
  if(obj)
    if(obj->type!=cObjectMolecule)
      obj=NULL;
  return((ObjectMolecule*)obj);
}
/*========================================================================*/
Block *ExecutiveGetBlock(void)
{
  CExecutive *I = &Executive;
  return(I->Block);
}
/*========================================================================*/
void ExecutiveSetControlsOff(char *name)
{
  SpecRec *rec;
  int a;
  rec = ExecutiveFindSpec(name);
  if(rec)
	 {
		for(a=0;a<cRepCnt;a++)
		  rec->repOn[a]=false;
	 }
}
/*========================================================================*/
void ExecutiveSymExp(char *name,char *oname,char *s1,float cutoff) /* TODO state */
{
  CObject *ob;
  ObjectMolecule *obj = NULL;
  ObjectMolecule *new_obj = NULL;
  ObjectMoleculeOpRec op;
  MapType *map;
  int x,y,z,a,b,c,i,j,h,k,l,n;
  CoordSet *cs;
  int keepFlag,sele,tt[3];
  float *v2,m[16],tc[3],ts[3];
  OrthoLineType new_name;
  float auto_save;

  PRINTFD(FB_Executive)
    " ExecutiveSymExp: entered.\n"
    ENDFD;

  auto_save = SettingGet(cSetting_auto_zoom);
  SettingSet(cSetting_auto_zoom,0);
  sele=SelectorIndexByName(s1);
  ob = ExecutiveFindObjectByName(oname);
  if(ob->type==cObjectMolecule)
    obj=(ObjectMolecule*)ob;
  if(!(obj&&sele)) {
    ErrMessage("ExecutiveSymExp","Invalid object");
  } else if(!obj->Symmetry) {
    ErrMessage("ExecutiveSymExp","No symmetry loaded!");
  } else if(!obj->Symmetry->NSymMat) {
    ErrMessage("ExecutiveSymExp","No symmetry matrices!");    
  } else {
    PRINTFB(FB_Executive,FB_Actions)
      " ExecutiveSymExp: Generating symmetry mates...\n"
      ENDFB;
    ObjectMoleculeOpRecInit(&op);
	 op.code = OMOP_SUMC;
	 op.i1 =0;
    op.i2 =0;
    op.v1[0]= 0.0;
    op.v1[1]= 0.0;
    op.v1[2]= 0.0;
    ExecutiveObjMolSeleOp(sele,&op);
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
    ExecutiveObjMolSeleOp(sele,&op);
    
    if(!op.nvv1) {
      ErrMessage("ExecutiveSymExp","No atoms indicated!");          
    } else {
      map=MapNew(-cutoff,op.vv1,op.nvv1,NULL);
      if(map) {
        MapSetupExpress(map);  

        for(x=-1;x<2;x++)
          for(y=-1;y<2;y++)
            for(z=-1;z<2;z++)
              for(a=0;a<obj->Symmetry->NSymMat;a++) {
                if(a||x||y||z) {
                  new_obj = ObjectMoleculeCopy(obj);
                  keepFlag=false;
                  for(b=0;b<new_obj->NCSet;b++) 
                    if(new_obj->CSet[b]) {
                      cs = new_obj->CSet[b];
                      CoordSetRealToFrac(cs,obj->Symmetry->Crystal);
                      CoordSetTransform44f(cs,obj->Symmetry->SymMatVLA+(a*16));
                      CoordSetGetAverage(cs,ts);
                      identity44f(m);
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
                    }
                  if(keepFlag) { /* need to create new object */
                    sprintf(new_name,"%s%02d%02d%02d%02d",name,a,x,y,z);
                    ObjectSetName((CObject*)new_obj,new_name);
                    ExecutiveDelete(new_name);
                    ExecutiveManageObject((CObject*)new_obj,true,false);
                    SceneChanged();
                  } else {
                    ((CObject*)new_obj)->fFree((CObject*)new_obj);
                  }
                }
              }
        MapFree(map);
      }
    }
    VLAFreeP(op.vv1);
  }
  PRINTFD(FB_Executive)
    " ExecutiveSymExp: leaving...\n"
    ENDFD;
  SettingSet(cSetting_auto_zoom,auto_save);
}
/*========================================================================*/
void ExecutiveDelete(char *name)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  int all_flag=false;
  WordType name_copy; /* needed in case the passed string changes */

  if(WordMatch(name,cKeywordAll,true)<0) all_flag=true;
  strcpy(name_copy,name);
  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->type==cExecObject)
		  {
          if(I->LastEdited==cExecObject) 
            I->LastEdited=NULL;
			 if(all_flag||(WordMatch(name_copy,rec->obj->Name,true)<0))
				{
              if(rec->obj==(CObject*)EditorGetActiveObject())
                EditorSetActiveObject(NULL,0);
              if(rec->visible) 
                SceneObjectDel(rec->obj);
				  SelectorDelete(rec->name);
				  rec->obj->fFree(rec->obj);
				  rec->obj=NULL;
				  ListDelete(I->Spec,rec,next,SpecList);
				  rec=NULL;
				}
		  }
		else if(rec->type==cExecSelection)
		  {

			 if(all_flag||(WordMatch(name_copy,rec->name,true)<0))
				{
              if(all_flag||rec->visible)
                SceneChanged();
				  SelectorDelete(rec->name);
				  ListDelete(I->Spec,rec,next,SpecList);
				  rec=NULL;
				}
		  }
	 }
  if(all_flag)
    SelectorDefragment();
}
/*========================================================================*/
void ExecutiveDump(char *fname,char *obj)
{
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;

  SceneUpdate();

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
        ErrMessage("ExecutiveDump","Invalid object type for this operation.");
      }
	 }
  else {
    ErrMessage("ExecutiveDump","Object not found.");
  }
  
}
/*========================================================================*/
void ExecutiveManageObject(CObject *obj,int allow_zoom,int quiet)
{
  int a;
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
  int exists=false;

  if(SettingGet(cSetting_auto_hide_selections))
    ExecutiveHideSelections();
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
        SceneObjectDel(rec->obj);
        rec->obj->fFree(rec->obj);
        rec->obj=NULL;
      }
    else 
      {
        if(!quiet)
          if(obj->Name[0]!='_') { /* suppress internal objects */
            PRINTFB(FB_Executive,FB_Actions)
              " Executive: object \"%s\" created.\n",obj->Name 
              ENDFB;
          }
      }
    if(!rec)
      ListElemAlloc(rec,SpecRec);

    if(WordMatch(cKeywordAll,obj->Name,true)<0) {
      PRINTFB(FB_Executive,FB_Warnings) 
        " Executive: object name \"%s\" is illegal -- renamed to 'all_'.",obj->Name
        ENDFB;
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
      SceneObjectAdd(obj);
    }
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    if(rec->obj->type==cObjectMolecule)
      rec->repOn[cRepLine]=true;
    ListAppend(I->Spec,rec,next,SpecList);
  }
  if(obj->type==cObjectMolecule) {
	 ExecutiveUpdateObjectSelection(obj);
  }
  if(allow_zoom)
    if(!exists) {
      switch(SettingGetGlobal_i(cSetting_auto_zoom)) {
      case 1: /* zoom new one */
        ExecutiveWindowZoom(obj->Name,0.0,-1,0); /* auto zoom (all states) */
        break;
      case 2: /* zoom all */
        ExecutiveWindowZoom(cKeywordAll,0.0,-1,0);
        break;
      }
    }

}
/*========================================================================*/
void ExecutiveManageSelection(char *name)
{

  int a;
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
  
  while(ListIterate(I->Spec,rec,next))
    {
      if(rec->type==cExecSelection)
        if(strcmp(rec->name,name)==0) 
          break;
    }
  if(!rec) {
    ListElemAlloc(rec,SpecRec);
    strcpy(rec->name,name);
    rec->type=cExecSelection;
    rec->next=NULL;
    rec->sele_color=-1;
    rec->visible=false;
    ListAppend(I->Spec,rec,next,SpecList);
  }
  if(rec) {
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    if(name[0]!='_') {
      if(SettingGet(cSetting_auto_hide_selections))
        ExecutiveHideSelections();
      if(SettingGet(cSetting_auto_show_selections)) {
        rec->visible=true;
      }
    }
    if(rec->visible) SceneDirty();
  }
}
/*========================================================================*/
int ExecutiveClick(Block *block,int button,int x,int y,int mod)
{
  CExecutive *I = &Executive;
  int n,a;
  SpecRec *rec = NULL;
  int t;
  int pass = false;
  int skip;
  n=((I->Block->rect.top-(y+2))-ExecTopMargin)/ExecLineHeight;
  a=n;
  if(I->ScrollBarActive) {
    if((x-I->Block->rect.left)<(ExecScrollBarWidth+ExecScrollBarMargin+ExecToggleMargin)) {
      pass = 1;
      ScrollBarDoClick(I->ScrollBar,button,x,y,mod);      
    }
  } 
  skip = I->NSkip;
  if(!pass)   while(ListIterate(I->Spec,rec,next))
    if(rec->name[0]!='_')
      {
        if(skip) {
          skip--;
        } else {
          if(!a) {
            t = ((I->Block->rect.right-ExecRightMargin)-x)/ExecToggleWidth;
            if(t<ExecOpCnt) {
              y = I->Block->rect.top-(ExecTopMargin + n*ExecLineHeight) - 1;
              x = I->Block->rect.right-(ExecRightMargin + t*ExecToggleWidth);
              t = (ExecOpCnt-t)-1;
              switch(t) {
              case 0:
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(x,y,"all_action",rec->name);
                  break;
                case cExecSelection:
                  MenuActivate(x,y,"sele_action",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectMolecule:
                    MenuActivate(x,y,"mol_action",rec->obj->Name);
                    break;
                  case cObjectSurface:
                  case cObjectMesh:
                  case cObjectDist:
                  case cObjectMap:
                  case cObjectCGO:
                  case cObjectCallback:
                    MenuActivate(x,y,"simple_action",rec->obj->Name);
                    break;
                  }
                  break;
                }
                break;
              case 1:
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(x,y,"mol_show",cKeywordAll);
                  break;
                case cExecSelection:
                  MenuActivate(x,y,"mol_show",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectMolecule:
                    MenuActivate(x,y,"mol_show",rec->obj->Name);
                    break;
                  case cObjectCGO:
                    MenuActivate(x,y,"cgo_show",rec->obj->Name);
                    break;
                  case cObjectDist:
                    MenuActivate(x,y,"dist_show",rec->obj->Name);
                    break;
                  case cObjectMap:
                    MenuActivate(x,y,"simple_show",rec->obj->Name);
                    break;
                  case cObjectSurface:
                  case cObjectMesh:
                    MenuActivate(x,y,"mesh_show",rec->obj->Name);
                    break;
                  }
                  break;
                }
                break;
              case 2:
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(x,y,"mol_hide",cKeywordAll);
                  break;
                case cExecSelection:
                  MenuActivate(x,y,"mol_hide",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectMolecule:
                    MenuActivate(x,y,"mol_hide",rec->obj->Name);
                    break;
                  case cObjectCGO:
                    MenuActivate(x,y,"cgo_hide",rec->obj->Name);
                    break;
                  case cObjectDist:
                    MenuActivate(x,y,"dist_hide",rec->obj->Name);
                    break;
                  case cObjectMap:
                    MenuActivate(x,y,"simple_hide",rec->obj->Name);
                    break;
                  case cObjectSurface:
                  case cObjectMesh:
                    MenuActivate(x,y,"mesh_hide",rec->obj->Name);
                    break;
                  }
                  break;
                }
                break;
              case 3:
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(x,y,"mol_labels","(all)");
                  break;
                case cExecSelection:
                  MenuActivate(x,y,"mol_labels",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectMolecule:
                    MenuActivate(x,y,"mol_labels",rec->obj->Name);
                    break;
                  case cObjectDist:
                    break;
                  case cObjectMap:
                  case cObjectSurface:
                  case cObjectMesh:
                    break;
                  }
                  break;
                }
                break;
              case 4:
                switch(rec->type) {
                case cExecAll:
                case cExecSelection:
                  MenuActivate(x,y,"mol_color",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectMolecule:
                    MenuActivate(x,y,"mol_color",rec->obj->Name);
                    break;
                  case cObjectDist:
                  case cObjectMap:
                  case cObjectSurface:
                  case cObjectCGO:
                  case cObjectMesh:
                    MenuActivate(x,y,"general_color",rec->obj->Name);
                    break;
                  }
                  break;
                }
                break;
              }
            }
          }
          a--;
        }
      }
  MainDirty();
  
  return(1);
}
/*========================================================================*/
int ExecutiveRelease(Block *block,int button,int x,int y,int mod)
{
  CExecutive *I = &Executive;
  int n;  
  SpecRec *rec = NULL;
  int t;
  OrthoLineType buffer;
  int pass = false;
  int skip;

  n=((I->Block->rect.top-(y+2))-ExecTopMargin)/ExecLineHeight;

  if(I->ScrollBarActive) {
    if((x-I->Block->rect.left)<(ExecScrollBarWidth+ExecScrollBarMargin+ExecToggleMargin)) {
      pass = 1;
      ScrollBarDoRelease(I->ScrollBar,button,x,y,mod);
      OrthoUngrab();
    }
  } 
  skip=I->NSkip;
  if(!pass) while(ListIterate(I->Spec,rec,next))
    if(rec->name[0]!='_')
      {
        if(skip) {
          skip--;
        } else 
          {
            if(!n) {
              {
                t = ((I->Block->rect.right-ExecRightMargin)-x)/ExecToggleWidth;
                if(t<ExecOpCnt) {
                  /* nothing to do anymore now that we have menus! */
                } else if(rec->type==cExecObject)
                  {
                    if(rec->visible)
                      SceneObjectDel(rec->obj);				
                    else 
                      SceneObjectAdd(rec->obj);
                    SceneChanged();
                    if(SettingGet(cSetting_logging)) {
                      if(rec->visible)
                        sprintf(buffer,"cmd.disable('%s')",rec->obj->Name);
                      else
                        sprintf(buffer,"cmd.enable('%s')",rec->obj->Name);
                      PLog(buffer,cPLog_pym);
                    }
                    rec->visible=!rec->visible;
                  }
                else if(rec->type==cExecAll)
                  {
                    if(SettingGet(cSetting_logging)) {
                      if(rec->visible)
                        sprintf(buffer,"cmd.disable('all')");
                      else
                        sprintf(buffer,"cmd.enable('all')");
                      PLog(buffer,cPLog_pym);
                    }
                    ExecutiveSetObjVisib(cKeywordAll,!rec->visible);
                  }
                else if(rec->type==cExecSelection)
                  {
                    if(mod&cOrthoCTRL) {
                      SettingSet(cSetting_selection_overlay,
                                 (float)(!((int)SettingGet(cSetting_selection_overlay))));
                      if(SettingGet(cSetting_logging)) {
                        sprintf(buffer,"cmd.set('selection_overlay',%d)",
                                (int)SettingGet(cSetting_selection_overlay));
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
                      if(SettingGet(cSetting_logging)) {
                        if(rec->visible)
                          sprintf(buffer,"cmd.disable('%s')",rec->name);
                        else
                          sprintf(buffer,"cmd.enable('%s')",rec->name);
                        PLog(buffer,cPLog_pym);
                      }
                      rec->visible=!rec->visible; 
                    }
                    SceneChanged();
                  }
              }
            }
            n--;
        }
      }
  MainDirty();
  return(1);
}
/*========================================================================*/
int ExecutiveDrag(Block *block,int x,int y,int mod)
{
  return(1);
}
/*========================================================================*/
void ExecutiveDraw(Block *block)
{
  int a,x,y,xx,x2,y2;
  char *c=NULL;
  float toggleColor[3] = { 0.5F, 0.5F, 1.0F };
  float toggleColor2[3] = { 0.3F, 0.3F, 0.6F };
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
  int n_ent;
  int n_disp;
  int skip=0;

  if(PMGUI) {

    
    /* do we have enough structures to warrant a scroll bar? */
    n_ent = 0;
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->name[0]!='_') 
        n_ent++;
    }

    n_disp = ((I->Block->rect.top-I->Block->rect.bottom)-(2+ExecTopMargin))/ExecLineHeight;
    if(n_disp<1) n_disp=1;
      
    if(n_ent>n_disp) {
      if(!I->ScrollBarActive) {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        ScrollBarSetValue(I->ScrollBar,0);
        I->NSkip =0;
      } else {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        
        I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
      }
      I->ScrollBarActive = 1;

    } else {
      I->ScrollBarActive = 0;
      I->NSkip =0;
    }
      
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);

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
    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(ExecOpCnt);
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
            x2=xx;
            y2=y-ExecToggleMargin;

            if((x-ExecToggleMargin)-(xx-ExecToggleMargin)>-10) {
              x2 = x+10;
            }
            glColor3fv(toggleColor);
            for(a=0;a<ExecOpCnt;a++)
              {
                switch(a) {
                case 0:
                  glColor3fv(toggleColor);
                  glBegin(GL_POLYGON);
                  glVertex2i(x2,y2+(ExecToggleSize)/2);
                  glVertex2i(x2+(ExecToggleSize)/2,y2);
                  glVertex2i(x2+ExecToggleSize,y2+(ExecToggleSize)/2);
                  glVertex2i(x2+(ExecToggleSize)/2,y2+ExecToggleSize);
                  glEnd();
                  break;
                case 1:
                  glColor3fv(toggleColor);
                  glBegin(GL_POLYGON);
                  glVertex2i(x2,y2);
                  glVertex2i(x2,y2+ExecToggleSize);
                  glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
                  glVertex2i(x2+ExecToggleSize,y2);
                  glEnd();
                  glColor3f(0.0,0.0,0.0);
                  glRasterPos4d((double)(x2+2),(double)(y2+2),0.0,1.0);
                  p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,'S');              
                  glColor3fv(toggleColor);
                  break;
                case 2:
                  glColor3fv(toggleColor2);
                  glBegin(GL_POLYGON);
                  glVertex2i(x2,y2);
                  glVertex2i(x2,y2+ExecToggleSize);
                  glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
                  glVertex2i(x2+ExecToggleSize,y2);
                  glEnd();
                  glColor3f(0.0,0.0,0.0);
                  glRasterPos4d((double)(x2+2),(double)(y2+2),0.0,1.0);
                  p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,'H');              
                  glColor3fv(toggleColor);
                  break;
                case 3:
                  glColor3fv(toggleColor);
                  glBegin(GL_POLYGON);
                  glVertex2i(x2,y2);
                  glVertex2i(x2,y2+ExecToggleSize);
                  glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
                  glVertex2i(x2+ExecToggleSize,y2);
                  glEnd();
                  glColor3f(0.0,0.0,0.0);
                  glRasterPos4d((double)(x2+2),(double)(y2+2),0.0,1.0);
                  p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,'L');              
                  glColor3fv(toggleColor);
                  break;
                case 4:
                  glBegin(GL_POLYGON);
                  glColor3f(1.0F,0.1F,0.1F);
                  glVertex2i(x2,y2);
                  glColor3f(0.1F,1.0F,0.1F);
                  glVertex2i(x2,y2+ExecToggleSize);
                  glColor3f(1.0F,1.0F,0.1F);
                  glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
                  glColor3f(0.1F,0.1F,1.0F);
                  glVertex2i(x2+ExecToggleSize,y2);
                  glEnd();
                  /*              glColor3f(0.0,0.0,0.0);
                                  glRasterPos4d((double)(x2+2),(double)(y2+2),0.0,1.0);
                                  p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,'C');              */
                  glColor3fv(toggleColor);
                  break;
                }
                x2+=ExecToggleWidth;
              }
        
#ifdef PYMOL_OLD_CODE 
            for(a=0;a<cRepCnt;a++)
              {
                if(rec->repOn[a]) 
                  {
                    glBegin(GL_POLYGON);
                    glVertex2i(x2,y2);
                    glVertex2i(x2,y2+ExecToggleSize);
                    glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
                    glVertex2i(x2+ExecToggleSize,y2);
                  }
                else
                  {
                    glBegin(GL_LINE_LOOP);
                    glVertex2i(x2,y2);
                    glVertex2i(x2,y2+ExecToggleSize-1);
                    glVertex2i(x2+ExecToggleSize-1,y2+ExecToggleSize-1);
                    glVertex2i(x2+ExecToggleSize-1,y2);
                  }
                glEnd();
                x2+=ExecToggleWidth;
              }
#endif

            glColor3fv(I->Block->TextColor);
            glRasterPos4d((double)(x),(double)(y),0.0,1.0);
            if((rec->type==cExecObject)||(rec->type==cExecAll)||(rec->type==cExecSelection))
              {
                y2=y-ExecToggleMargin;
                if(rec->visible)
                  glColor3f(ExecGreyVisible,ExecGreyVisible,ExecGreyVisible);
                else
                  glColor3f(ExecGreyHidden,ExecGreyHidden,ExecGreyHidden);
                x2 = xx;
                if((x-ExecToggleMargin)-(xx-ExecToggleMargin)>-10) {
                  x2 = x+10;
                }
                glBegin(GL_POLYGON);
                glVertex2i(x-ExecToggleMargin,y2);
                glVertex2i(x2-ExecToggleMargin,y2);
                glVertex2i(x2-ExecToggleMargin,y2+ExecToggleSize);
                glVertex2i(x-ExecToggleMargin,y2+ExecToggleSize);
                glEnd();
                glColor3fv(I->Block->TextColor);

                if(rec->type!=cExecObject)
                  c=rec->name;
                else 
                  c=rec->obj->Name;

                if(rec->type==cExecSelection)
                  p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,'(');
              }

            if(c)
              while(*c) 
                p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c++));

            if(rec->type==cExecSelection)
              {
                p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,')');
                c=rec->name;
              }

            y-=ExecLineHeight;
            if(y<(I->Block->rect.bottom+2))
              break;
          }
        }
  }
}
/*========================================================================*/
int ExecutiveIterateObject(CObject **obj,void **hidden)
{
  int result;
  CExecutive *I = &Executive;
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
int ExecutiveIterateObjectMolecule(ObjectMolecule **obj,void **hidden)
{
  int result;
  CExecutive *I = &Executive;
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
void ExecutiveReshape(Block *block,int width,int height)
{
  CExecutive *I = &Executive;

  BlockReshape(block,width,height);

  I->Width = block->rect.right-block->rect.left+1;
  I->Height = block->rect.top-block->rect.bottom+1;
  
}

/*========================================================================*/
int ExecutiveReinitialize(void)
{ 
  int ok=true;
  int blocked = false;
  /* reinitialize PyMOL */

  ExecutiveDelete(cKeywordAll);
  ColorReset();
  SettingInitGlobal(false,false);
  MovieReset();
  EditorInactive();

  blocked = PAutoBlock();
  PRunString("cmd.view('*','clear')");
  PRunString("cmd.scene('*','clear')");
  WizardSet(NULL,false);
  PAutoUnblock(blocked);
  
  SculptCachePurge();
  SceneReinitialize();
  SelectorReinit();

  return(ok);
}
/*========================================================================*/
void ExecutiveInit(void)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  int a;

  ListInit(I->Spec);
  I->Block = OrthoNewBlock(NULL);  
  I->Block->fRelease = ExecutiveRelease;
  I->Block->fClick   = ExecutiveClick;
  I->Block->fDrag    = ExecutiveDrag;
  I->Block->fDraw    = ExecutiveDraw;
  I->Block->fReshape = ExecutiveReshape;
  I->Block->active = true;
  I->ScrollBarActive = 0;
  I->ScrollBar=ScrollBarNew(false);
  OrthoAttach(I->Block,cOrthoTool);

  I->LastEdited=NULL;
  I->NSkip=0;
  ListElemAlloc(rec,SpecRec);
  strcpy(rec->name,"(all)");
  rec->type=cExecAll;
  rec->visible=true;
  rec->next=NULL;
  for(a=0;a<cRepCnt;a++)
	 rec->repOn[a]=false;
  ListAppend(I->Spec,rec,next,SpecList);

}
/*========================================================================*/
void ExecutiveFree(void)
{
  CExecutive *I = &Executive;
  SpecRec *rec=NULL;
  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->type==cExecObject)
		  rec->obj->fFree(rec->obj);
	 }
  ListFree(I->Spec,next,SpecList);
  if(I->ScrollBar)
    ScrollBarFree(I->ScrollBar);
  OrthoFreeBlock(I->Block);
  I->Block=NULL;
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


