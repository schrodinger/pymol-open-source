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

#include"os_std.h"
#include"os_gl.h"

#include"main.h"
#include"Base.h"
#include"OOMac.h"
#include"Executive.h"
#include"ObjectMesh.h"
#include"ObjectDist.h"
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

#include"Menu.h"
#include"Map.h"
#include"Editor.h"
#include"RepDot.h"

#define cExecObject 0
#define cExecSelection 1
#define cExecAll 2

#define cKeywordAll "all"

typedef struct SpecRec {
  int type;
  WordType  name; /*only used for selections*/
  struct Object *obj;  
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
  Object *LastEdited;
} CExecutive;

CExecutive Executive;

int ExecutiveClick(Block *block,int button,int x,int y,int mod);
int ExecutiveRelease(Block *block,int button,int x,int y,int mod);
int ExecutiveDrag(Block *block,int x,int y,int mod);
void ExecutiveDraw(Block *block);
void ExecutiveReshape(Block *block,int width,int height);

#define ExecLineHeight 14
#define ExecTopMargin 0
#define ExecToggleMargin 2
#define ExecLeftMargin 5
#define ExecRightMargin 2
#define ExecToggleWidth 14
#define ExecToggleSize 13

#define ExecOpCnt 5
#define ExecColorVisible 0.45,0.45,0.45
#define ExecColorHidden 0.3,0.3,0.3

void ExecutiveObjMolSeleOp(int sele,ObjectMoleculeOpRec *op);
SpecRec *ExecutiveFindSpec(char *name);

float ExecutiveAlign(char *s1,char *s2)
{
  int sele1=SelectorIndexByName(s1);
  int sele2=SelectorIndexByName(s2);
  int *vla1=NULL;
  int *vla2=NULL;
  int na,nb;
  int a,b,c;
  float result = 0.0;
  CMatch *match = NULL;

  if((sele1>=0)&&(sele2>=0)) {
    vla1=SelectorGetResidueVLA(sele1);
    vla2=SelectorGetResidueVLA(sele2);
    if(vla1&&vla2) {
      na = VLAGetSize(vla1)/3;
      nb = VLAGetSize(vla2)/3;
      if(na&&nb) {
        match = MatchNew(na,nb);
        for(a=0;a<na;a++) {
          for(b=0;b<nb;b++) {
            if(vla1[a*3+2] == vla2[b*3+2])
              match->mat[a][b] = 1.0;
            else
              match->mat[a][b] = -0.1;
          }
        }
        result = MatchAlign(match,-0.2,-0.1,10);
        if(match->pair) { /* alignment was successful */
          c = SelectorCreateAlignments(match->pair,
                                       sele1,vla1,sele2,vla2,
                                       "align1","align2");
          if(c) {
            PRINTFB(FB_Executive,FB_Actions)
              " ExecutiveAlign: %d atoms aligned.\n",c
              ENDFB;
              ExecutiveRMS("align1","align2",2,1.5);
          }
        }
        if(match) 
          MatchFree(match);
      }
    }
  }
  VLAFreeP(vla1);
  VLAFreeP(vla2);
  return 1.0;
}

void ExecutiveFocus(void)
{ /* unfortunately, this doesn't achieve the desired effect */
  if(PMGUI) {
    glutPopWindow();
    glutShowWindow();
  }
}

PyObject *ExecutiveGetSettingText(int index,char *object,int state)
{
  PyObject *result;
  OrthoLineType buffer;

  if(object[0]==0) /* global */ {
    buffer[0]=0;
    SettingGetTextValue(NULL,NULL,index,buffer);
    result=Py_BuildValue("s",buffer);
  } else {
    /* TODO */
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveGetSettingTuple(int index,char *object,int state)
{
  PyObject *result;
  if(object[0]==0) /* global */
    result = SettingGetTuple(NULL,NULL,index);
  else {
    /* TODO */
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return(result);
}
/*========================================================================*/
void ExecutiveSetLastObjectEdited(Object *o)
{
  CExecutive *I = &Executive;
  I->LastEdited = o;
}
/*========================================================================*/
Object *ExecutiveGetLastObjectEdited(void)
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
  OrthoDirty();
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

  no_depth = SettingGet(cSetting_selection_overlay);
  width = SettingGet(cSetting_selection_width);

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection) {

      if(rec->visible) {
        sele = SelectorIndexByName(rec->name); /* TODO: speed this up */
        if(sele>=0) {
          rec1 = NULL;
          if(rec->sele_color<0)
            glColor3f(1.0,0.2,0.8);
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
  float current;

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
    /* UNFINISHED */
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
char *ExecutiveGetNames(int mode)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  char *result;
  int size = 0;
  int stlen;
  result=VLAlloc(char,1000);

  while(ListIterate(I->Spec,rec,next)) {
    if(
       (rec->type==cExecObject&&((!mode)||(mode==1)))||
       (rec->type==cExecSelection&&((!mode)||(mode==2))))
    {
      stlen = strlen(rec->name);
      VLACheck(result,char,size+stlen+1);
      strcpy(result+size,rec->name);
      size+=stlen+1;
    }
  }
  VLASize(result,char,size);
  return(result);
}
/*========================================================================*/
void ExecutiveGetType(char *name,WordType type)
{
  SpecRec *rec = NULL;

  rec = ExecutiveFindSpec(name);
  if(!rec)
    strcpy(type,"nonexistent");
  else {
    if(rec->type==cExecObject) {
      strcpy(type,"object:");
      if(rec->obj->type==cObjectMolecule)
        strcat(type,"molecule");
      else if(rec->obj->type==cObjectMap)
        strcat(type,"map");
      else if(rec->obj->type==cObjectMesh)
        strcat(type,"mesh");
      else if(rec->obj->type==cObjectDist)
        strcat(type,"distance");
    } else if(rec->type==cExecSelection) {
      strcpy(type,"selection");
    }
  }
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
  Object *os=NULL;
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
void ExecutiveInvert(char *s0,char *s1,int mode)
{
  int i0=-1;
  int i1=-1;
  int sele0,sele1;
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
      EditorInvert(obj0,sele0,sele1,mode);
    }
  }
}
/*========================================================================*/
void ExecutiveFuse(char *s0,char *s1)
{
  int i0=-1;
  int i1=-1;
  int sele0,sele1;
  ObjectMolecule *obj0,*obj1;

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
void ExecutiveSpheroid(char *name)  /* EXPERIMENTAL */
{
  CExecutive *I = &Executive;
  Object *os=NULL;
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
            ObjectMoleculeCreateSpheroid(obj);  
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
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject)
      if(rec->obj->type==cObjectMolecule) {
        ObjectMoleculeInvalidate((ObjectMolecule*)rec->obj,cRepAll,cRepInvRep);
      }
  }
  SceneDirty();
}
/*========================================================================*/
void ExecutiveUndo(int dir)
{
  CExecutive *I = &Executive;
  Object *o;
  ObjectMolecule *obj=NULL,*compObj;
  SpecRec *rec = NULL;

  o = ExecutiveGetLastObjectEdited();
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
  Object *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;
  ObjectMoleculeOpRec op;
  int sele;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(name);
    if(!os)
      ErrMessage(" Executive","object not found.");
    else if(os->type!=cObjectMolecule)
      ErrMessage(" Executive","bad object type.");
  }
  
  if(os||(!strlen(name))) { /* sort one or all */
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          if((!os)||(rec->obj==os)) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeSort(obj);
            sele=SelectorIndexByName(rec->obj->Name);
            if(sele>=0) {
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
                  op.code = OMOP_Remove;
                  op.i1 = 0;
						obj=(ObjectMolecule*)rec->obj;
                  ObjectMoleculeVerifyChemistry(obj); /* remember chemistry for later */
						ObjectMoleculeSeleOp(obj,sele,&op);
                  if(op.i1) {
                    ObjectMoleculePurge(obj);
                    PRINTFB(FB_Editor,FB_Actions)
                      " Remove: eliminated %d atoms in model '%s'.\n",
                      op.i1,obj->Obj.Name 
                      ENDFB;
                    flag=true;
                  }
					 }
				}
		  }
	 }
  if(!flag) {
    ErrMessage("Remove","no atoms removed.");
  }
}
/*========================================================================*/
void ExecutiveAddHydrogens(char *s1)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(s1);
  if(sele1>=0) {
    op.code = OMOP_AddHydrogens;
    ExecutiveObjMolSeleOp(sele1,&op);    
  }
}
/*========================================================================*/
void ExecutiveFlag(int flag,char *s1)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(s1);
  if(sele1>=0) {
    op.code = OMOP_Flag;
    op.i1 = (((unsigned int)1)<<flag);
    op.i2 = ((unsigned int)0xFFFFFFFF - (((unsigned int)1)<<flag));
    op.i3 = 0;
    ExecutiveObjMolSeleOp(sele1,&op);    
    if(Feedback(FB_Executive,FB_Actions)) {
      if(op.i3) {
        PRINTF " Flag: flag %d set on %d atoms.\n", flag, op.i3 ENDF;
      } else {
        PRINTF " Flag: flag %d cleared on all atoms.\n", flag ENDF;
      }
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
      op.i1 = cRepLine;
      op.i2 = cRepInvPick;
      ExecutiveObjMolSeleOp(sele1,&op);    
  }
}
/*========================================================================*/
void ExecutiveStereo(int flag)
{

  if(PMGUI) {
    if(StereoCapable) {
      SceneSetStereo(flag);
    }
  }
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
                        " AddBond: %d bonds added to model '%s'.\n",cnt,rec->obj->Name 
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
                        " RemoveBond: %d bonds removed from model '%s'.\n",
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
float ExecutiveDist(char *nam,char *s1,char *s2,int mode,float cutoff)
{
  int sele1,sele2;
  ObjectDist *obj;
  float result;
  sele1=SelectorIndexByName(s1);
  sele2=SelectorIndexByName(s2);
  
  if((sele1>=0)&&(sele2>=0)) {
    obj = ObjectDistNew(sele1,sele2,mode,cutoff,&result);
    if(!obj) {
      ErrMessage("ExecutiveDistance","No such distances found.");
    } else {
      if(ExecutiveFindObjectByName(nam))
        ExecutiveDelete(nam);
      ObjectSetName((Object*)obj,nam);
      ExecutiveManageObject((Object*)obj);
      ExecutiveSetRepVisib(nam,cRepLine,1);
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
  
  sele1=SelectorIndexByName(s1);
  op1.i1=0;
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
    scale3f(op1.v1,1.0/op1.i1,op1.v1);
    scale3f(op2.v1,1.0/op2.i1,op2.v1);
    dist = diff3f(op1.v1,op2.v1);
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

  sele1=SelectorIndexByName(s1);
  op1.charVLA=VLAlloc(char,10000);
  if(conectFlag) {
    if(state<0) state=0;
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
  l=strlen(end_str);
  VLACheck(op1.charVLA,char,op1.i2+l+1);
  strcpy(op1.charVLA+op1.i2,end_str);
  op1.i2+=l+1;
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
  PBlockAndUnlockAPI();
  result=SelectorGetChemPyModel(sele1,state);
  if(PyErr_Occurred()) PyErr_Print();
  PLockAPIAndUnblock();
  return(result);
}
/*========================================================================*/
void ExecutiveSeleToObject(char *name,char *s1,int source,int target)
{
  int sele1;

  sele1=SelectorIndexByName(s1);

  if(sele1>=0) {
    SelectorCreateObjectMolecule(sele1,name,target,source);
  }
}
/*========================================================================*/
void ExecutiveCopy(char *src,char *dst)
{
  Object *os;
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
        ExecutiveManageObject((Object*)oDst);
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
void ExecutiveOrient(char *sele,Matrix33d mi)
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
		  m[a][b]=evect[b][a]; /* fill columns */
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

    ExecutiveWindowZoom(sele,0.0);

  }
}
/*========================================================================*/
void ExecutiveLabel(char *s1,char *expr)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int cnt;

  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
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

    PRINTFB(FB_Executive,FB_Actions)
      " Label: labelled %i atoms.\n",cnt
      ENDFB;
  } else {
    PRINTFB(FB_Executive,FB_Warnings)
      " Label: no atoms selections.\n"
      ENDFB;
  }
}
/*========================================================================*/
int ExecutiveIterate(char *s1,char *expr,int read_only)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  
  op1.i1=0;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    op1.code = OMOP_ALTR;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = read_only;
    ExecutiveObjMolSeleOp(sele1,&op1);
    if(!read_only) {
      PRINTFB(FB_Executive,FB_Actions)
        " Alter: modified %i atoms.\n",op1.i1
        ENDFB;
    } else {
      PRINTFB(FB_Executive,FB_Actions)
        " Iterate: iterated over %i atoms.\n",op1.i1
        ENDFB;
    }
  } else {
      PRINTFB(FB_Executive,FB_Warnings)
        "ExecutiveIterate: No atoms selected.\n"
        ENDFB;
  }
  return(op1.i1);
}
/*========================================================================*/
void ExecutiveIterateState(int state,char *s1,char *expr,int read_only)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    op1.code = OMOP_AlterState;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = state;
    op1.i3 = read_only;
    ExecutiveObjMolSeleOp(sele1,&op1);
    if(!read_only) {
      PRINTFB(FB_Executive,FB_Actions)
        " AlterState: modified %i atom states.\n",op1.i1
        ENDFB;
    } else {
      PRINTFB(FB_Executive,FB_Actions)
        " IterateState: iterated over %i atom states.\n",op1.i1
        ENDFB;
    }
  } else {
      PRINTFB(FB_Executive,FB_Warnings)
        "ExecutiveIterateState: No atoms selected.\n"
        ENDFB;
  }
}
/*========================================================================*/
float ExecutiveRMS(char *s1,char *s2,int mode,float refine)
{
  int sele1,sele2;
  float rms = -1.0;
  int a;
  float inv,*f,*f1,*f2;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType buffer;
  int *flag;

  sele1=SelectorIndexByName(s1);
  op1.vv1=NULL;
  op1.vc1=NULL;
  op2.vv1=NULL;
  op2.vc1=NULL;
  if(sele1>=0) {
    op1.code = OMOP_AVRT;
    op1.nvv1=0;
    op1.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    ExecutiveObjMolSeleOp(sele1,&op1);
    for(a=0;a<op1.nvv1;a++)
      {
        inv=op1.vc1[a]; /* average over coordinate sets */
        if(inv)
          {
            f=op1.vv1+(a*3);
            inv=1.0/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }

  sele2=SelectorIndexByName(s2);
  if(sele2>=0) {
    op2.code = OMOP_AVRT;
    op2.nvv1=0;
    op2.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    ExecutiveObjMolSeleOp(sele2,&op2);
    for(a=0;a<op2.nvv1;a++)
      {
        inv=op2.vc1[a]; /* average over coordinate sets */
        if(inv)
          {
            f=op2.vv1+(a*3);
            inv=1.0/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }
  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1==op2.nvv1) {
      if(refine>R_SMALL4) {
        flag = MatrixFilter(refine,20,2,op1.nvv1,op1.vv1,op2.vv1);

        /* eliminate unflagged vertices -- they're considered outliers */

        if(flag) {
          f1 = op1.vv1;
          f2 = op2.vv1;
          for(a=0;a<op1.nvv1;a++) {
            if(!flag[a]) {
              op2.nvv1--;
            } else {
              f1+=3;
              f2+=3;
            }
          copy3f(op1.vv1+3*a,f1);
          copy3f(op2.vv1+3*a,f2);
          }
          if(op2.nvv1!=op1.nvv1) {
            PRINTFB(FB_Executive,FB_Actions)
              " ExecutiveRMS: %d atoms thrown out..\n",op1.nvv1-op2.nvv1
              ENDFB;
          }
          op1.nvv1 = op2.nvv1;
          FreeP(flag);
        }
      }
    }
  }

  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1!=op2.nvv1) {
      sprintf(buffer,"Atom counts between selections don't match (%d vs %d)\n",
              op1.nvv1,op2.nvv1);
      ErrMessage("ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(mode!=0) 
        rms = MatrixFitRMS(op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
      else
        rms = MatrixGetRMS(op1.nvv1,op1.vv1,op2.vv1,NULL);
      PRINTFB(FB_Executive,FB_Results) 
        " Executive: RMS = %8.3f (%d to %d atoms)\n", rms,op1.nvv1,op2.nvv1 
        ENDFB
      if(mode==2) {
        op2.code = OMOP_TTTF;
        ExecutiveObjMolSeleOp(sele1,&op2);
      }
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
int *ExecutiveIdentify(char *s1,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  int *result = NULL;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
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
float *ExecutiveRMSStates(char *s1,int target,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  float *result = NULL;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    op1.code = OMOP_SVRT;
    op1.nvv1=0;
    op1.i1=target;
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    ExecutiveObjMolSeleOp(sele1,&op1);

    op2.vv2=op1.vv1;
    op2.nvv2=op1.nvv1;
    op2.i2=target;
    op2.i1=mode;
    op2.f1VLA=VLAlloc(float,10);
    VLASetSize(op2.f1VLA,0); /* failsafe */
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    op2.code = OMOP_SFIT;
    op2.nvv1=0;
    ExecutiveObjMolSeleOp(sele1,&op2);
    result=op2.f1VLA;
  } 
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
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
      inv=op1.vc1[a];
      if(inv)
        {
          f=op1.vv1+(a*3);
          inv=1.0/inv;
          *(f++)*=inv;
          *(f++)*=inv;
          *(f++)*=inv;
        }
    }
  for(a=0;a<op2.nvv1;a++)
    {
      inv=op2.vc1[a];
      if(inv)
        {
          f=op2.vv1+(a*3);
          inv=1.0/inv;
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
void ExecutiveUpdateObjectSelection(struct Object *obj)
{
  if(obj->type==cObjectMolecule) {
    SelectorUpdateObjectSele((ObjectMolecule*)obj);  
  }
}
/*========================================================================*/
void ExecutiveReset(int cmd)
{
  SceneResetMatrix();
  ExecutiveWindowZoom(cKeywordAll,0.0);
}
/*========================================================================*/
void ExecutiveDrawNow(void) 
{
  if(PMGUI) {
    glMatrixMode(GL_MODELVIEW);
  /*  glClear( GL_DEPTH_BUFFER_BIT);*/
  }

  SceneUpdate();


  OrthoDoDraw();


  MainSwapBuffers();


}
/*========================================================================*/
int ExecutiveCountStates(char *s1)
{
  int sele1;
  int result=0;
  sele1=SelectorIndexByName(s1);
  if(sele1>=0) {
    SelectorUpdateTable();
    result = SelectorGetSeleNCSet(sele1);
  }
  return(result);
}
/*========================================================================*/
void ExecutiveRay(void)
{
  SceneRay();
}
/*========================================================================*/
void ExecutiveSetSetting(int index,PyObject *tuple,char *sele,
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

  unblock = PAutoBlock();
  if(sele[0]==0) { 
    SettingSetTuple(NULL,index,tuple);
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
  else if(!strcmp(cKeywordAll,sele)) { /* all objects setting */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject) {
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
          if(handle) {
            SettingCheckHandle(handle);
            SettingSetTuple(*handle,index,tuple);
            if(updates) 
              SettingGenerateSideEffects(index,sele,state);
            nObj++;
          }
          }
        }
        if(Feedback(FB_Setting,FB_Actions)) {
          if(nObj&&handle) {
            SettingGetTextValue(*handle,NULL,index,value);
            SettingGetName(index,name);
            if(updates)
              SettingGenerateSideEffects(index,sele,state);
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
  }
 else { /* based on a selection/object name */
    sele1=SelectorIndexByName(sele);
    while((sele1>=0)&&(ListIterate(I->Spec,rec,next)))
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          {
            obj=(ObjectMolecule*)rec->obj;
            op.code=OMOP_CountAtoms;
            op.i1=0;
            ObjectMoleculeSeleOp(obj,sele1,&op);
            if(op.i1&&rec->obj->fGetSettingHandle) {
              handle = rec->obj->fGetSettingHandle(rec->obj,state);
              if(handle) {
                SettingCheckHandle(handle);
                SettingSetTuple(*handle,index,tuple);
                if(!quiet) {
                  if(state<0) { /* object-specific */
                    if(Feedback(FB_Setting,FB_Actions)) {
                      SettingGetTextValue(*handle,NULL,index,value);
                      SettingGetName(index,name);
                      PRINTF
                        " Setting: %s set to %s in object '%s'.\n",
                        name,value,rec->obj->Name
                        ENDF;
                      if(updates)
                        SettingGenerateSideEffects(index,sele,state);
                      
                    }
                  } else { /* state-specific */
                    if(Feedback(FB_Setting,FB_Actions)) {
                      SettingGetTextValue(*handle,NULL,index,value);
                      SettingGetName(index,name);
                      PRINTF
                        " Setting: %s set to %s in object '%s', state %d.\n",
                        name,value,rec->obj->Name,state+1
                        ENDF;
                      if(updates) 
                        SettingGenerateSideEffects(index,sele,state);
                      
                    }
                  }
                }
              }
            }
          }
 }
  PAutoUnblock(unblock);
}
/*========================================================================*/
void ExecutiveColor(char *name,char *color,int flags)
{
  SpecRec *rec = NULL;
  int sele;
  ObjectMoleculeOpRec op;
  int col_ind;

  col_ind = ColorGetIndex(color);
  if(col_ind<0) {
    ErrMessage("Color","unknown color.");
  } else {
    /* per atom */
    if(!(flags&0x1)) {
      sele=SelectorIndexByName(name);
      if(sele>=0) {
        op.code = OMOP_COLR;
        op.i1= col_ind;
        ExecutiveObjMolSeleOp(sele,&op);
        op.code=OMOP_INVA;
        op.i1=cRepAll; 
        op.i2=cRepInvColor;
        ExecutiveObjMolSeleOp(sele,&op);
      }
    }
    /* per object */
    rec=ExecutiveFindSpec(name);
    if(rec) {
      if(rec->type==cExecObject) {
        rec->obj->Color=col_ind;
        SceneDirty();
      }
    }
  }
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
int ExecutiveGetExtent(char *name,float *mn,float *mx)
{
  int sele;
  ObjectMoleculeOpRec op,op2;
  CExecutive *I=&Executive;
  Object *obj;
  int flag = false;
  int all_flag = false;
  SpecRec *rec = NULL;
  WordType all = "_all";
  float f1,f2,fmx;
  int a;



  op2.i1 = 0;
  op2.v1[0]=-1.0;
  op2.v1[1]=-1.0;
  op2.v1[2]=-1.0;
  op2.v2[0]=1.0;
  op2.v2[1]=1.0;
  op2.v2[2]=1.0;
  
  if(WordMatch(cKeywordAll,name,true)<0) {
    name=all;
    all_flag=true;
    SelectorCreate(all,"(all)",NULL,true);
  }
  sele=SelectorIndexByName(name);

  if(sele>=0) {
	 op.code = OMOP_MNMX;
	 op.v1[0]=FLT_MAX;
	 op.v1[1]=FLT_MAX;
	 op.v1[2]=FLT_MAX;
    op.v2[0]=FLT_MIN;
    op.v2[1]=FLT_MIN;
    op.v2[2]=FLT_MIN;
    op.i1 = 0;
    ExecutiveObjMolSeleOp(sele,&op);
    if(op.i1)
      flag = true;
    if(all_flag) {
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          obj=rec->obj;
          if(obj->ExtentFlag) 
            switch(obj->type) {
            case cObjectMesh:
            case cObjectMap:
              min3f(obj->ExtentMin,op.v1,op.v1);
              max3f(obj->ExtentMax,op.v2,op.v2);
              flag = true;
              break;
            }
        }
      }
    }
    op2.i1=0;
	 op2.code = OMOP_SUMC;
	 op2.v1[0]=0.0;
	 op2.v1[1]=0.0;
	 op2.v1[2]=0.0;
    ExecutiveObjMolSeleOp(sele,&op2);
    if(op2.i1) {
      op2.v1[0]/=op2.i1;
      op2.v1[1]/=op2.i1;
      op2.v1[2]/=op2.i1;
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
          break;
        }
      }
    }
  }
  if(all_flag) {
    ExecutiveDelete(all);
  }

  if(flag) { 
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
  return(flag);  
}
/*========================================================================*/
void ExecutiveWindowZoom(char *name,float buffer)
{
  float center[3],radius;
  float mn[3],mx[3];

  if(ExecutiveGetExtent(name,mn,mx)) {
    if(buffer!=0.0) {
      buffer = buffer;
      mx[0]+=buffer;
      mx[1]+=buffer;
      mx[2]+=buffer;
      mn[0]-=buffer;
      mn[1]-=buffer;
      mn[2]-=buffer;
    }
    radius = diff3f(mn,mx)/3.0;
    average3f(mn,mx,center);
    if(radius<MAX_VDW)
      radius=MAX_VDW;
    SceneOriginSet(center,false);
    SceneWindowSphere(center,radius);
    SceneDirty();
  }
}
/*========================================================================*/
void ExecutiveCenter(char *name,int preserve)
{
  float center[3];
  float mn[3],mx[3];

  if(ExecutiveGetExtent(name,mn,mx)) {
    average3f(mn,mx,center);
    SceneOriginSet(center,preserve);
    SceneDirty();
  }
}
/*========================================================================*/
int ExecutiveGetMoment(char *name,Matrix33d mi)
{
  int sele;
  ObjectMoleculeOpRec op;
  int a,b;
  int c=0;
  for(a=0;a<3;a++)
	 {
		for(b=0;b<3;b++)
		  mi[a][b]=0.0;
		mi[a][a]=1.0;
	 }
  
  sele=SelectorIndexByName(name);
  if(sele>=0) {
	 op.code = OMOP_SUMC;
	 op.v1[0]=0.0;
	 op.v1[1]=0.0;
	 op.v1[2]=0.0;
	 op.i1=0;
	 
	 ExecutiveObjMolSeleOp(sele,&op);
	 
	 if(op.i1) {
		c+=op.i1;
		scale3f(op.v1,1.0/op.i1,op.v1);
		op.code = OMOP_MOME;		
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
  if(strcmp(name,cKeywordAll)==0) {
    tRec=NULL;
    while(ListIterate(I->Spec,tRec,next)) {
      if(state!=tRec->visible) {
        if(tRec->type==cExecObject) {
          if(tRec->visible)
            SceneObjectDel(tRec->obj);				
          else 
            SceneObjectAdd(tRec->obj);
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
            else 
              SceneObjectAdd(tRec->obj);
            tRec->visible=!tRec->visible;
          }
      }
      else if(tRec->type==cExecSelection) {
        if(!(tRec->visible&&state)) {
          tRec->visible=!tRec->visible;
          SceneChanged();
        }
      }
    }
  }
}
/*========================================================================*/
void ExecutiveFullScreen(int flag)
{
  if(PMGUI) {
    if(flag) {
      glutFullScreen();
    } else {
      glutReshapeWindow(640+cOrthoRightSceneMargin,
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
  while(ListIterate(I->Spec,rec,next)) {
	 if(rec->type==cExecObject)
		{
		  if(rec->obj->type==cObjectMolecule)
			 {
				obj=(ObjectMolecule*)rec->obj;
				sele = SelectorIndexByName(obj->Obj.Name);
				for(rep=0;rep<cRepCnt;rep++) {
				  rec->repOn[rep]=state;
				  op.code=OMOP_VISI;
				  op.i1=rep;
				  op.i2=state;
				  ObjectMoleculeSeleOp(obj,sele,&op);
				  op.code=OMOP_INVA;
				  op.i2=cRepInvVisib;
				  ObjectMoleculeSeleOp(obj,sele,&op);				
				}
			 }
		}
  }
}
/*========================================================================*/
void ExecutiveSetRepVisib(char *name,int rep,int state)
{
  int sele;
  int a;
  int handled = false;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;

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
      if(tRec->obj->type==cObjectDist)
        {
          ObjectSetRepVis(tRec->obj,rep,state);
          SceneDirty();
        }
    if(!handled)
      switch(tRec->type) {
      case cExecSelection:
      case cExecObject:
        sele=SelectorIndexByName(name);
        if(sele>=0) {
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
          if(rec->obj->type==cObjectMolecule)
            {
              obj=(ObjectMolecule*)rec->obj;
              sele = SelectorIndexByName(obj->Obj.Name);
              rec->repOn[rep]=state;
              op.code=OMOP_VISI;
              op.i1=rep;
              op.i2=state;
              ObjectMoleculeSeleOp(obj,sele,&op);
              op.code=OMOP_INVA;
              op.i2=cRepInvVisib;
              ObjectMoleculeSeleOp(obj,sele,&op);				
            }
          else if(rec->obj->type==cObjectDist) {
            ObjectSetRepVis(rec->obj,rep,state);
            SceneDirty();
          }
        }
		}
  }
}
/*========================================================================*/
void ExecutiveInvalidateRep(char *name,int rep,int level)
{
  int sele = -1;
  ObjectMoleculeOpRec op;
  WordType all = "_all";
  int all_flag=false;
  if(WordMatch(cKeywordAll,name,true)<0) {
    name=all;
    all_flag=true;
    SelectorCreate(all,"(all)",NULL,true);
  }
  sele=SelectorIndexByName(name);
  if(sele>=0) {
	 op.code = OMOP_INVA;
	 op.i1=rep;
	 op.i2=level;
	 ExecutiveObjMolSeleOp(sele,&op);
  }
  if(all_flag) {
    ExecutiveDelete(all);
  }
}
/*========================================================================*/
Object *ExecutiveFindObjectByName(char *name)
{
  CExecutive *I = &Executive;
  SpecRec *rec = NULL;
  Object *obj=NULL;
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
ObjectMolecule *ExecutiveFindObjectMoleculeByName(char *name)
{
  Object *obj;
  
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
void ExecutiveSymExp(char *name,char *oname,char *s1,float cutoff)
{
  Object *ob;
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
      " ExecutiveSymExp: Generating symmetry mates"
      ENDFB;
	 op.code = OMOP_SUMC;
	 op.i1 =0;
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
                      m[3] = tt[0]+x;
                      m[7] = tt[1]+y;
                      m[11] = tt[2]+z;
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
                    ObjectSetName((Object*)new_obj,new_name);
                    ExecutiveManageObject((Object*)new_obj);
                    SceneChanged();
                  } else {
                    ((Object*)new_obj)->fFree((Object*)new_obj);
                  }
                }
              }
        MapFree(map);
      }
    }
    VLAFreeP(op.vv1);
  }
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
              if(rec->obj==(Object*)EditorGetActiveObject())
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
				  SelectorDelete(rec->name);
				  ListDelete(I->Spec,rec,next,SpecList);
				  rec=NULL;
				}
		  }
	 }
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
      } else {
        ErrMessage("ExecutiveDump","Invalid object type for this operation.");
      }
	 }
  else {
    ErrMessage("ExecutiveDump","Object not found.");
  }
  
}
/*========================================================================*/
void ExecutiveManageObject(Object *obj)
{
  int a;
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;
  int exists=false;

  if(SettingGet(cSetting_autohide_selections))
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
        if(obj->Name[0]!='_') { /* suppress internal objects */
          PRINTFB(FB_Executive,FB_Actions)
            " Executive: object \"%s\" created.\n",obj->Name 
            ENDFB;
        }
      }
    if(!rec)
      ListElemAlloc(rec,SpecRec);
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
  if(!exists) 
    if(SettingGet(cSetting_auto_zoom)) {
      ExecutiveWindowZoom(obj->Name,0.0);
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
    rec->visible=true;
    ListAppend(I->Spec,rec,next,SpecList);
  }
  if(rec) {
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    if(name[0]!='_') {
      if(SettingGet(cSetting_autohide_selections))
        ExecutiveHideSelections();
      if(SettingGet(cSetting_autoshow_selections)) {
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

  n=((I->Block->rect.top-(y+2))-ExecTopMargin)/ExecLineHeight;
  a=n;

  while(ListIterate(I->Spec,rec,next))
    if(rec->name[0]!='_')
	 {
		if(!a)
		  {
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
                case cObjectDist:
                case cObjectMap:
                case cObjectMesh:
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
                case cObjectDist:
                  MenuActivate(x,y,"dist_show",rec->obj->Name);
                  break;
                case cObjectMap:
                case cObjectMesh:
                  MenuActivate(x,y,"simple_show",rec->obj->Name);
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
                case cObjectDist:
                  MenuActivate(x,y,"dist_hide",rec->obj->Name);
                  break;
                case cObjectMap:
                case cObjectMesh:
                  MenuActivate(x,y,"simple_hide",rec->obj->Name);
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

  n=((I->Block->rect.top-(y+2))-ExecTopMargin)/ExecLineHeight;

  while(ListIterate(I->Spec,rec,next))
    if(rec->name[0]!='_')
      {
        if(!n)
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
                rec->visible=!rec->visible;
                SceneChanged();
              }
            else if(rec->type==cExecAll)
              {
                ExecutiveSetObjVisib(cKeywordAll,!rec->visible);
              }
            else if(rec->type==cExecSelection)
              {
                if(mod&cOrthoCTRL) {
                  SettingSet(cSetting_selection_overlay,
                             (float)(!((int)SettingGet(cSetting_selection_overlay))));
                  rec->visible=true;
                } else if(mod&cOrthoSHIFT) {
                  if(rec->sele_color<0)
                    rec->sele_color=7;
                  else {
                    rec->sele_color--;
                    if(rec->sele_color<7)
                      rec->sele_color=15;
                  }
                  rec->visible=true;
                } else {
                  rec->visible=!rec->visible; 
                }
                SceneChanged();
              }

          }
        n--;
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
  float toggleColor[3] = { 0.5, 0.5, 1.0 };
  float toggleColor2[3] = { 0.3, 0.3, 0.6 };
  SpecRec *rec = NULL;
  CExecutive *I = &Executive;

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    
    x = I->Block->rect.left+ExecLeftMargin;
    y = (I->Block->rect.top-ExecLineHeight)-ExecTopMargin;
    /*    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(cRepCnt+ExecOpCnt);*/
    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(ExecOpCnt);
    
    while(ListIterate(I->Spec,rec,next))
      if(rec->name[0]!='_')
      {
        x2=xx;
        y2=y-ExecToggleMargin;


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
              glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'S');              
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
              glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'H');              
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
              glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'L');              
              glColor3fv(toggleColor);
              break;
            case 4:
              glBegin(GL_POLYGON);
              glColor3f(1.0,0.1,0.1);
              glVertex2i(x2,y2);
              glColor3f(0.1,1.0,0.1);
              glVertex2i(x2,y2+ExecToggleSize);
              glColor3f(1.0,1.0,0.1);
              glVertex2i(x2+ExecToggleSize,y2+ExecToggleSize);
              glColor3f(0.1,0.1,1.0);
              glVertex2i(x2+ExecToggleSize,y2);
              glEnd();
              /*              glColor3f(0.0,0.0,0.0);
              glRasterPos4d((double)(x2+2),(double)(y2+2),0.0,1.0);
              glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'C');              */
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
              glColor3f(ExecColorVisible);
            else
              glColor3f(ExecColorHidden);
            glBegin(GL_POLYGON);
            glVertex2i(x-ExecToggleMargin,y2);
            glVertex2i(xx-ExecToggleMargin,y2);
            glVertex2i(xx-ExecToggleMargin,y2+ExecToggleSize);
            glVertex2i(x-ExecToggleMargin,y2+ExecToggleSize);
            glEnd();
            glColor3fv(I->Block->TextColor);

            if(rec->type!=cExecObject)
              c=rec->name;
            else 
              c=rec->obj->Name;

            if(rec->type==cExecSelection)
              glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'(');
          }

        if(c)
          while(*c) 
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));

        if(rec->type==cExecSelection)
          {
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13,')');
            c=rec->name;
          }

        y-=ExecLineHeight;
        if(y<(I->Block->rect.bottom+2))
          break;
      }
  }
}
/*========================================================================*/
int ExecutiveIterateObject(Object **obj,void **hidden)
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
void ExecutiveReshape(Block *block,int width,int height)
{
  CExecutive *I = &Executive;

  BlockReshape(block,width,height);

  I->Width = block->rect.right-block->rect.left+1;
  I->Height = block->rect.top-block->rect.bottom+1;
  
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

  OrthoAttach(I->Block,cOrthoTool);

  I->LastEdited=NULL;
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


