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

#include"PConv.h"
#include"MemoryDebug.h"
#include"Vector.h"
#include"Matrix.h"
#include"ButMode.h"
#include"Scene.h"
#include"Editor.h"
#include"Selector.h"
#include"Ortho.h"
#include"main.h"
#include"Color.h"
#include"Setting.h"
#include"Util.h"
#include"Executive.h"
#include"P.h"

typedef struct {
  ObjectMolecule *Obj;
  WordType DragSeleName;
  int ActiveState;
  int DragIndex;
  int DragSelection;
  int DragHaveAxis,DragHaveBase,DragBondFlag,DragSlowFlag;
  int PickMode; /* 1 = atom, 2 = bond, 3 = multiatom */
  int NextPickSele;
  int BondMode;
  ObjectMolecule *DragObject;
  int NFrag;
  float V0[3],V1[3],Axis[3],Center[3],DragBase[3];
  float *PosVLA;
  int ShowFrags;
}  CEditor;

CEditor Editor;

int EditorGetNFrag(void)
{
  CEditor *I = &Editor;
  if(EditorActive()) {
    return I->NFrag;
  }
  return 0;
}

void EditorDefineExtraPks(void)
{
  WordType name;
  WordType buffer;

  if(EditorGetSinglePicked(name)) {
    sprintf(buffer,"(byres %s)",name);
    SelectorCreate(cEditorRes,buffer,NULL,true,NULL);
    sprintf(buffer,"(bychain %s)",name);
    SelectorCreate(cEditorChain,buffer,NULL,true,NULL);
    sprintf(buffer,"(byobject %s)",name);
    SelectorCreate(cEditorObject,buffer,NULL,true,NULL);
    
    if(SettingGet(cSetting_auto_hide_selections))
      ExecutiveHideSelections();
  }
}

int EditorDeselectIfSelected(int index,int update)
{
  CEditor *I = &Editor;
  int result = false;
  int s,sele;
  if(I->Obj) {
    if((index>=0)&&(index<I->Obj->NAtom)) {
      s=I->Obj->AtomInfo[index].selEntry;                  
      sele = SelectorIndexByName(cEditorSele1);
      if(SelectorIsMember(s,sele)) {
        ExecutiveDelete(cEditorSele1);
        result = true;
      }
      sele = SelectorIndexByName(cEditorSele2);
      if(SelectorIsMember(s,sele)) {
        ExecutiveDelete(cEditorSele2);
        result = true;
      }
      sele = SelectorIndexByName(cEditorSele3);
      if(SelectorIsMember(s,sele)) {
        ExecutiveDelete(cEditorSele3);
        result = true;
      }
      sele = SelectorIndexByName(cEditorSele4);
      if(SelectorIsMember(s,sele)) {
        ExecutiveDelete(cEditorSele4);
        result = true;
      }
      if(result&&update)
        EditorSetActiveObject(I->Obj,I->ActiveState);
    }
  }
  
  return result;
}

int EditorIsBondMode(void)
{
  CEditor *I = &Editor;
  return(I->BondMode);
}

PyObject *EditorAsPyList(void)
{
  PyObject *result = NULL;
  CEditor *I = &Editor;

  if(!EditorActive()) {
    result = PyList_New(0); /* not editing? return null list */
  } else {
    result = PyList_New(2);
    PyList_SetItem(result,0,PyString_FromString(I->Obj->Obj.Name));
    PyList_SetItem(result,1,PyInt_FromLong(I->ActiveState));
  }
  return(PConvAutoNone(result));
}

int EditorFromPyList(PyObject *list)
{
  int ok=true;
  int active_flag = false;
  ObjectMolecule *objMol;
  int active_state;
  WordType obj_name;
  int ll;
  
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  if(ok) active_flag=(PyList_Size(list)!=0);
  if(!active_flag) {
    EditorInactive();
  } else {
    if(ok) ok=PConvPyStrToStr(PyList_GetItem(list,0),obj_name,sizeof(WordType));
    if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,1),&active_state);
    if(ok) {
      objMol=ExecutiveFindObjectMoleculeByName(obj_name);
      if(objMol) {
        EditorSetActiveObject(objMol,active_state);
        EditorDefineExtraPks();
      }
    } else {
      EditorInactive();
    }
  }
  if(!ok) {
    EditorInactive();
  }
  return(ok);
}

int EditorActive(void) {
  CEditor *I = &Editor;
  if(!I->Obj)
    return 0;
  if(I->Obj) {
    return 1;
  }

  return(I->Obj!=NULL);
}

ObjectMolecule *EditorDragObject(void)
{
  CEditor *I = &Editor;
  return(I->DragObject);
}
static void subdivide( int n, float *x, float *y);

int EditorGetSinglePicked(char *name)
{
  int cnt = 0;
  int sele;
  if((sele = SelectorIndexByName(cEditorSele1))>=0) {
    cnt++;
    if(name) strcpy(name,cEditorSele1);
  }
  if((sele = SelectorIndexByName(cEditorSele2))>=0) {
    cnt++;
    if(name) strcpy(name,cEditorSele2);
  }
  if((sele = SelectorIndexByName(cEditorSele3))>=0) {
    cnt++;
    if(name) strcpy(name,cEditorSele3);
  }
  if((sele = SelectorIndexByName(cEditorSele4))>=0) {
    cnt++;
    if(name) strcpy(name,cEditorSele4);
  }
  return (cnt==1);
}

int EditorIsObjectNotCurrent(ObjectMolecule *obj)
{
  CEditor *I = &Editor;
  if(!I->Obj) return 0;
  return (I->Obj!=obj);
}

void EditorGetNextMultiatom(char *name)
{
  CEditor *I = &Editor;
  int sele;
  sele = SelectorIndexByName(cEditorSele1);
  if(sele<0) {
    strcpy(name,cEditorSele1);
    I->NextPickSele=0;
    return;
  }
  sele = SelectorIndexByName(cEditorSele2);
  if(sele<0) {
    strcpy(name,cEditorSele2);
    I->NextPickSele=1;
    return;
  }
  sele = SelectorIndexByName(cEditorSele3);
  if(sele<0) {
    strcpy(name,cEditorSele3);
    I->NextPickSele=2;
    return;
  }
  sele = SelectorIndexByName(cEditorSele4);
  if(sele<0) {
    strcpy(name,cEditorSele4);
    I->NextPickSele=3;
    return;
  }
  I->NextPickSele = (++I->NextPickSele)&0x3;
  switch(I->NextPickSele) {
  case 0: strcpy(name,cEditorSele1); break;
  case 1: strcpy(name,cEditorSele2); break;
  case 2: strcpy(name,cEditorSele3); break;
  case 3: strcpy(name,cEditorSele4); break;
  }
  return;    
}

/*========================================================================*/
int EditorInvert(ObjectMolecule *obj,int isele0,int isele1,int mode)
{
  CEditor *I = &Editor;
  int sele0,sele1,sele2,sele3;
  int i0,frg;
  int a,s,n,a1;
  int if0=-1;
  int if1=-1;
  int ia0=-1;
  int ia1=-1;
  float v[3],v0[3],v1[3];
  float n0[3],n1[3];
  float m[16];
  int state;
  int vf,vf0,vf1;
  int ok=false;
  WordType name,bname;

  if((!I->Obj)||(I->Obj!=obj)) { 
    ErrMessage("Editor","Must pick an atom to invert.");
  } else {
    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      if(sele1>=0) {
        ErrMessage("Editor","Must edit an atom, not a bond.");        
      } else if(!(sele0>=0)) {
        ErrMessage("Editor","Must pick an atom to invert.");        
      } else {
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0);
        
        if(i0>=0) {
          for(frg=1;frg<=I->NFrag;frg++) {
            sprintf(name,"%s%1d",cEditorFragPref,frg);
            sele2=SelectorIndexByName(name);
            if(sele2>=0) {
              /* now find out which bonded atoms are immobilized */
              for(a=0;a<obj->NAtom;a++) {
                s=obj->AtomInfo[a].selEntry;                  
                if(SelectorIsMember(s,sele2)) { /* atom in fragment */
                  if(if0<0) {
                    if(SelectorIsMember(s,isele0)) { /* atom in (lb) */
                      if0=frg;
                      sprintf(bname,"%s%1d",cEditorBasePref,frg); 
                      sele3 = SelectorIndexByName(bname);
                      ia0 = ObjectMoleculeGetAtomIndex(obj,sele3);
                      /* this is the atom we want */
                    }
                  }
                  if(if1<0) {
                    if(SelectorIsMember(s,isele1)) { /* atom in (rb) */
                      if1=frg;
                      sprintf(bname,"%s%1d",cEditorBasePref,frg); 
                      sele3 = SelectorIndexByName(bname);
                      ia1 = ObjectMoleculeGetAtomIndex(obj,sele3);
                      /* this is the atom we want */
                    }
                  }
                  if((if0>=0)&&(if1>=0)) /* we got 'em */
                    break;
                }
              }
            }
          }
          /* make sure the anchor atoms are unique (for instance,
               * in the case of a cycle they won't be )*/

          if((ia0>=0)&&(ia1>=0)&&(i0>=0)&&(ia0==ia1)) {
            ObjectMoleculeUpdateNeighbors(obj);
            ia1=-1;
            sprintf(name,"%s%1d",cEditorFragPref,if0); 
            sele3 = SelectorIndexByName(name);
            n = obj->Neighbor[i0];
            n++; /* skip count */
            while(1) { /* look for another attached atom as a second base */
              a1 = obj->Neighbor[n];
              n+=2; 
              if(a1<0) break;
              if(a1!=ia0) {
                s=obj->AtomInfo[a1].selEntry;
                if(SelectorIsMember(s,sele3)) {
                  ia1 = a1;
                  break;
                }
              }
            }
          }
          if((ia0<0)||(ia1<0)||(i0<0))
            ErrMessage("Invert","couldn't find basis for inversion");
          else {

            state = SceneGetState();                
            ObjectMoleculeSaveUndo(obj,state,false);

            vf  = ObjectMoleculeGetAtomVertex(obj,state,i0,v);
            vf0 = ObjectMoleculeGetAtomVertex(obj,state,ia0,v0);
            vf1 = ObjectMoleculeGetAtomVertex(obj,state,ia1,v1);
                
            if(vf&vf0&vf1) {
              subtract3f(v,v0,n0);
              subtract3f(v,v1,n1);
              normalize3f(n0);
              normalize3f(n1);
                
              add3f(n0,n1,n0);
              normalize3f(n0);

              MatrixRotation44f(m,(float)cPI,n0[0],n0[1],n0[2]);
              m[3 ] = -v[0];
              m[7 ] = -v[1];
              m[11] = -v[2];
              m[12] =  v[0];
              m[13] =  v[1];
              m[14] =  v[2];

              for(frg=1;frg<=I->NFrag;frg++)
                switch(mode) {
                case 0 :
                  if((frg!=if0)&(frg!=if1)) {
                    sprintf(name,"%s%1d",cEditorFragPref,frg);
                    sele2=SelectorIndexByName(name);
                    ok = ObjectMoleculeTransformSelection(obj,state,sele2,m,false,NULL);
                  }
                  break;
                case 1:
                  if((frg!=if0)&(frg!=if1)) {
                    sprintf(name,"%s%1d",cEditorFragPref,frg);
                    sele2=SelectorIndexByName(name);
                    ok = ObjectMoleculeTransformSelection(obj,state,sele2,m,false,NULL);
                  }
                  break;
                }
              SceneDirty();
              I->DragIndex=-1;
              I->DragSelection=-1;
              I->DragObject=NULL;
            }
          }
        }
      }
    }
  }
  return(ok);
}
/*========================================================================*/
int EditorTorsion(float angle)
{
  CEditor *I = &Editor;
  int sele0,sele1,sele2;
  int i0,i1;
  float v0[3],v1[3];
  float d1[3],n0[3];
  float theta;
  float m[16];
  int state;
  int vf1,vf2;
  int ok=false;
  WordType sele;

  if(!I->Obj) { 
    ErrMessage("Editor","Must specify a bond first.");
  } else {
    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      strcpy(sele,cEditorFragPref);
      strcat(sele,"1");
      sele2 = SelectorIndexByName(sele);

      if(!(sele0>=0)&&(sele1>=0)&&(sele2>=0)) {
        ErrMessage("Editor","Must specify a bond first.");
      } else {
        
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0);
        i1 = ObjectMoleculeGetAtomIndex(I->Obj,sele1);
        if((i0>=0)&&(i1>=0)) {
          state = SceneGetState();
          
          vf1 = ObjectMoleculeGetAtomVertex(I->Obj,state,i0,I->V0);
          vf2 = ObjectMoleculeGetAtomVertex(I->Obj,state,i1,I->V1);
          
          if(vf1&&vf2) {
            ObjectMoleculeSaveUndo(I->Obj,SceneGetState(),false);
            
            
            subtract3f(I->V1,I->V0,I->Axis);
            average3f(I->V1,I->V0,I->Center);
            normalize3f(I->Axis);
            
            copy3f(I->V0,v1);
            copy3f(I->V1,v0);
            
            subtract3f(v1,v0,d1);
            copy3f(d1,n0);
            normalize3f(n0);
            
            theta=(float)(cPI*angle/180.0);
            MatrixRotation44f(m,theta,n0[0],n0[1],n0[2]);
            m[3 ] = -v1[0];
            m[7 ] = -v1[1];
            m[11] = -v1[2];
            m[12] =  v1[0];
            m[13] =  v1[1];
            m[14] =  v1[2];
            ok = ObjectMoleculeTransformSelection(I->Obj,state,sele2,m,false,NULL);
            SceneDirty();
            
            I->DragIndex=-1;
            I->DragSelection=-1;
            I->DragObject=NULL;
          }
        }
      }
    }
  }
  return(ok);
}

/*========================================================================*/
int EditorSelect(char *s0,char *s1,char *s2,char *s3,int pkresi,int quiet)
{
  int i0=-1;
  int i1=-1;
  int i2=-1;
  int i3=-1;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int result=false;
  int ok = true;
  ObjectMolecule *obj0=NULL,*obj1=NULL,*obj2=NULL,*obj3=NULL,*consensus_obj=NULL;

  if(s0)
    if(!*s0)
      s0=NULL;
  if(s1)
    if(!*s1)
      s1=NULL;
  if(s2)
    if(!*s2)
      s2=NULL;
  if(s3)
    if(!*s3)
      s3=NULL;

  if(s0) {
    sele0 = SelectorIndexByName(s0);
    obj0 = SelectorGetSingleObjectMolecule(sele0);
    if(obj0)
      i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
    ExecutiveDelete(cEditorSele1);
  }
 
  if(s1) {
    sele1 = SelectorIndexByName(s1);
    if(sele1>=0) {
      obj1 = SelectorGetSingleObjectMolecule(sele1);
      if(obj1)
        i1 = ObjectMoleculeGetAtomIndex(obj1,sele1);
    }
    ExecutiveDelete(cEditorSele2);
  }

  if(s2) {
    sele2 = SelectorIndexByName(s2);
    if(sele2>=0) {
      obj2 = SelectorGetSingleObjectMolecule(sele2);
      if(obj2)
        i2 = ObjectMoleculeGetAtomIndex(obj2,sele2);
    }
    ExecutiveDelete(cEditorSele3);
  }

  if(s3) {
    sele3 = SelectorIndexByName(s3);
    if(sele3>=0) {
      obj3 = SelectorGetSingleObjectMolecule(sele3);
      if(obj3)
        i3 = ObjectMoleculeGetAtomIndex(obj3,sele3);
    }
    ExecutiveDelete(cEditorSele4);
  }

  if(ok) {
    if(obj0) {
      consensus_obj = obj0;
    } 
    if(obj1) {
      if(!consensus_obj)
        consensus_obj = obj1;
      else if(consensus_obj!=obj1)
        ok=false;
    }
    if(obj2) {
      if(!consensus_obj)
        consensus_obj = obj2;
      else if(consensus_obj!=obj2)
        ok=false;
    }
    if(obj3) {
      if(!consensus_obj)
        consensus_obj = obj3;
      else if(consensus_obj!=obj3)
        ok=false;
    }
  }

  if(!consensus_obj) ok=false;
  if(ok) {
    ObjectMoleculeVerifyChemistry(consensus_obj);  
    
    if(i0>=0) SelectorCreate(cEditorSele1,s0,NULL,quiet,NULL);
    if(i1>=0) SelectorCreate(cEditorSele2,s1,NULL,quiet,NULL);
    if(i2>=0) SelectorCreate(cEditorSele3,s2,NULL,quiet,NULL);
    if(i3>=0) SelectorCreate(cEditorSele4,s3,NULL,quiet,NULL);
    
    EditorSetActiveObject(consensus_obj,SceneGetState());        

    if(pkresi)
      EditorDefineExtraPks();

    SceneDirty();
    result=true;

  } else {
    EditorSetActiveObject(NULL,0);
    ErrMessage("Editor","Invalid input.");    
  }
  return(result);
}
/*========================================================================*/
static void subdivide( int n, float *x, float *y)
{
  int a;
  if(n<3) {n=3;}
  for(a=0;a<=n;a++)
	 {
		x[a]=(float)cos(a*2*PI/n);
		y[a]=(float)sin(a*2*PI/n);
	 }
}

/*========================================================================*/

ObjectMolecule *EditorGetActiveObject(void) {
  CEditor *I = &Editor;
  return(I->Obj);
}
/*========================================================================*/
void EditorCycleValence(void)
{
  CEditor *I = &Editor;
  int sele0,sele1;
  
  if(I->Obj) {

    ObjectMoleculeVerifyChemistry(I->Obj); /* remember chemistry for later */

    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      if(sele1>=0) {
        /* bond mode */
        ObjectMoleculeAdjustBonds(I->Obj,sele0,sele1,0,0);
      }
    }
  }

}
/*========================================================================*/
void EditorAttach(char *elem,int geom,int valence,char *name)
{
  CEditor *I = &Editor;
  int i0;
  int sele0,sele1;
  int state;
  AtomInfoType *ai;
  
  ai=(AtomInfoType*)VLAMalloc(1,sizeof(AtomInfoType),1,true);
  if(I->Obj) {

    ObjectMoleculeVerifyChemistry(I->Obj); /* remember chemistry for later */

    state = SceneGetState();

    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      if(sele1>=0) {
        /* bond mode - behave like replace */
        EditorReplace(elem,geom,valence,name);
      } else {
        /* atom mode */
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0); /* slow */
        if(i0>=0) {
          UtilNCopy(ai->elem,elem,sizeof(AtomName));
          ai->geom=geom;
          ai->valence=valence;
          if(name[0])
            UtilNCopy(ai->name,name,sizeof(AtomName));
          ObjectMoleculeAttach(I->Obj,i0,ai); /* will free ai */
          ai = NULL;
        }
      }
    }
  }
  VLAFreeP(ai); /* safety */
}
/*========================================================================*/
void EditorRemove(int hydrogen)
{
  CEditor *I = &Editor;
  int sele0,sele1;
  int i0;
  int h_flag = false;
  OrthoLineType buf;

#define cEditorRemoveSele "_EditorRemove"

  if(I->Obj) {
    ObjectMoleculeVerifyChemistry(I->Obj); /* remember chemistry for later */
    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);

      if(sele1>=0) {
        /* bond mode */
        ObjectMoleculeRemoveBonds(I->Obj,sele0,sele1);
        EditorSetActiveObject(NULL,0);        
      } else {

        if(hydrogen) {
          sprintf(buf,"((neighbor %s) and hydro)",cEditorSele1);          
          h_flag = SelectorCreate(cEditorRemoveSele,buf,NULL,false,NULL);
        }

        /* atom mode */
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0); /* slow */
        if(i0>=0) {
          ExecutiveRemoveAtoms(cEditorSele1);
          EditorSetActiveObject(NULL,0);
        }

        if(h_flag) {
          ExecutiveRemoveAtoms(cEditorRemoveSele);
          SelectorDelete(cEditorRemoveSele);
        }

      }
    }
  }
#undef cEditorRemoveSele
}
/*========================================================================*/
void EditorHFill(void)
{

  CEditor *I = &Editor;
  int sele0,sele1;
  int i0;
  OrthoLineType buffer,s1;
  
  if(I->Obj) {
    ObjectMoleculeVerifyChemistry(I->Obj); /* remember chemistry for later */
    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      if(sele0>=0) {
        if(sele1>=0) 
          sprintf(buffer,"((neighbor %s) and (elem h) and not %s)",
                  cEditorSele1,cEditorSele2);
        else 
          sprintf(buffer,"((neighbor %s) and (elem h))",
                  cEditorSele1);
        SelectorGetTmp(buffer,s1);
        ExecutiveRemoveAtoms(s1);
        SelectorFreeTmp(s1);
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0); 
        I->Obj->AtomInfo[i0].chemFlag=false;
        ExecutiveAddHydrogens(cEditorSele1);
      }
      
      if(sele1>=0) {
        if(sele0>=0) 
          sprintf(buffer,"((neighbor %s) and (elem h) and not %s)",
                  cEditorSele2,cEditorSele1);
        else 
          sprintf(buffer,"((neighbor %s) and (elem h))",
                  cEditorSele2);
        SelectorGetTmp(buffer,s1);
        ExecutiveRemoveAtoms(s1);
        SelectorFreeTmp(s1);
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele1); 
        I->Obj->AtomInfo[i0].chemFlag=false;
        ExecutiveAddHydrogens(cEditorSele2);
      }
    }
  }
  
}
/*========================================================================*/
void EditorReplace(char *elem,int geom,int valence,char *name)
{
  CEditor *I = &Editor;
  int i0;
  int sele0;
  int state;
  AtomInfoType ai;
  
  UtilZeroMem(&ai,sizeof(AtomInfoType));
  if(I->Obj) {

    ObjectMoleculeVerifyChemistry(I->Obj); /* remember chemistry for later */

    state = SceneGetState();

    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele0); /* slow */
      if(i0>=0) {
        UtilNCopy(ai.elem,elem,sizeof(AtomName));
        if(name[0])
          UtilNCopy(ai.name,name,sizeof(AtomName));
        ai.geom=geom;
        ai.valence=valence;
        ObjectMoleculePrepareAtom(I->Obj,i0,&ai);
        ObjectMoleculePreposReplAtom(I->Obj,i0,&ai);
        ObjectMoleculeReplaceAtom(I->Obj,i0,&ai); /* invalidates */
        ObjectMoleculeVerifyChemistry(I->Obj);
        ObjectMoleculeFillOpenValences(I->Obj,i0);
        ObjectMoleculeSort(I->Obj);
        ObjectMoleculeUpdateIDNumbers(I->Obj);
        EditorSetActiveObject(NULL,0);
      }
    }
  }
}

static void draw_bond(float *v0,float *v1)
{
  
  float v[3],v2[3],v3[3];
  float d0[3],n0[3],n1[3],n2[3];
  float x[50],y[50];
  int nEdge;
  int c,a;
  float tube_size1=0.5F;
  float tube_size3=0.45F;

  nEdge = (int)SettingGet(cSetting_stick_quality)*2;
  if(nEdge>50)
    nEdge=50;
  
  subdivide(nEdge,x,y);


  subtract3f(v1,v0,d0);
  average3f(v1,v0,v2);
  average3f(v0,v2,v3);
  average3f(v2,v3,v2);
  copy3f(d0,n0);
  get_system1f3f(n0,n1,n2);
  
  glColor3fv(ColorGet(0));
  glBegin(GL_TRIANGLE_STRIP);
  for(a=0;a<=nEdge;a++) {
    c=a % nEdge;
    v[0] =  n1[0]*x[c] + n2[0]*y[c];
    v[1] =  n1[1]*x[c] + n2[1]*y[c];
    v[2] =  n1[2]*x[c] + n2[2]*y[c];
    normalize3f(v);
    glNormal3fv(v);
    v[0] = v2[0] + n1[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c];
    v[1] = v2[1] + n1[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c];
    v[2] = v2[2] + n1[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c];
    glVertex3fv(v);
    
    v[0] = v3[0] + n1[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c];
    v[1] = v3[1] + n1[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c];
    v[2] = v3[2] + n1[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c];
    glVertex3fv(v);
  }
  glEnd();
  
  glBegin(GL_TRIANGLE_STRIP);
  glNormal3fv(n0);
  for(a=0;a<=nEdge;a++) {
    c=a % nEdge;
    v[0] = v2[0] + n1[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c];
    v[1] = v2[1] + n1[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c];
    v[2] = v2[2] + n1[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c];
    glVertex3fv(v);
    v[0] = v2[0] + n1[0]*tube_size3*x[c] + n2[0]*tube_size3*y[c];
    v[1] = v2[1] + n1[1]*tube_size3*x[c] + n2[1]*tube_size3*y[c];
    v[2] = v2[2] + n1[2]*tube_size3*x[c] + n2[2]*tube_size3*y[c];
    glVertex3fv(v);
  }
  glEnd();
  
  glBegin(GL_TRIANGLE_STRIP);
  scale3f(n0,-1.0F,v);
  glNormal3fv(v);
  for(a=0;a<=nEdge;a++) {
    c=a % nEdge;
    v[0] = v3[0] + n1[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c];
    v[1] = v3[1] + n1[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c];
    v[2] = v3[2] + n1[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c];
    glVertex3fv(v);
    v[0] = v3[0] + n1[0]*tube_size3*x[c] + n2[0]*tube_size3*y[c];
    v[1] = v3[1] + n1[1]*tube_size3*x[c] + n2[1]*tube_size3*y[c];
    v[2] = v3[2] + n1[2]*tube_size3*x[c] + n2[2]*tube_size3*y[c];
    glVertex3fv(v);
  }
  glEnd();

}

static void draw_globe(float *v2)
{
  float v[3];
  float n0[3],n1[3],n2[3];
  float x[50],y[50];
  int nEdge;
  int a,c;
  float tube_size1=0.5F;
  float tube_size2=0.07F;

  nEdge = (int)SettingGet(cSetting_stick_quality)*2;
  if(nEdge>50)
    nEdge=50;
  
  subdivide(nEdge,x,y);
  
  
  n0[0]=1.0;
  n0[1]=0.0;
  n0[2]=0.0;
  get_system1f3f(n0,n1,n2);
  
  glColor3fv(ColorGet(0));
  glBegin(GL_TRIANGLE_STRIP);
  for(a=0;a<=nEdge;a++) {
    c=a % nEdge;
    v[0] =  n1[0]*x[c] + n2[0]*y[c];
    v[1] =  n1[1]*x[c] + n2[1]*y[c];
    v[2] =  n1[2]*x[c] + n2[2]*y[c];
    normalize3f(v);
    glNormal3fv(v);
    v[0] = v2[0] + n1[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c]+n0[0]*tube_size2;
    v[1] = v2[1] + n1[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c]+n0[1]*tube_size2;
    v[2] = v2[2] + n1[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c]+n0[2]*tube_size2;
    glVertex3fv(v);
    v[0] = v2[0] + n1[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c]-n0[0]*tube_size2;
    v[1] = v2[1] + n1[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c]-n0[1]*tube_size2;
    v[2] = v2[2] + n1[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c]-n0[2]*tube_size2;
    glVertex3fv(v);
  }
  glEnd();
  
  glBegin(GL_TRIANGLE_STRIP);
  for(a=0;a<=nEdge;a++) {
    c=a % nEdge;
    v[0] =  n0[0]*x[c] + n2[0]*y[c];
    v[1] =  n0[1]*x[c] + n2[1]*y[c];
    v[2] =  n0[2]*x[c] + n2[2]*y[c];
    normalize3f(v);
    glNormal3fv(v);
    v[0] = v2[0] + n0[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c]+n1[0]*tube_size2;
    v[1] = v2[1] + n0[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c]+n1[1]*tube_size2;
    v[2] = v2[2] + n0[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c]+n1[2]*tube_size2;
    glVertex3fv(v);
    v[0] = v2[0] + n0[0]*tube_size1*x[c] + n2[0]*tube_size1*y[c]-n1[0]*tube_size2;
    v[1] = v2[1] + n0[1]*tube_size1*x[c] + n2[1]*tube_size1*y[c]-n1[1]*tube_size2;
    v[2] = v2[2] + n0[2]*tube_size1*x[c] + n2[2]*tube_size1*y[c]-n1[2]*tube_size2;
    glVertex3fv(v);
  }
  glEnd();
  
  glBegin(GL_TRIANGLE_STRIP);
  for(a=0;a<=nEdge;a++) {
    c=a % nEdge;
    v[0] =  n0[0]*x[c] + n1[0]*y[c];
    v[1] =  n0[1]*x[c] + n1[1]*y[c];
    v[2] =  n0[2]*x[c] + n1[2]*y[c];
    normalize3f(v);
    glNormal3fv(v);
    v[0] = v2[0] + n0[0]*tube_size1*x[c] + n1[0]*tube_size1*y[c]+n2[0]*tube_size2;
    v[1] = v2[1] + n0[1]*tube_size1*x[c] + n1[1]*tube_size1*y[c]+n2[1]*tube_size2;
    v[2] = v2[2] + n0[2]*tube_size1*x[c] + n1[2]*tube_size1*y[c]+n2[2]*tube_size2;
    glVertex3fv(v);
    v[0] = v2[0] + n0[0]*tube_size1*x[c] + n1[0]*tube_size1*y[c]-n2[0]*tube_size2;
    v[1] = v2[1] + n0[1]*tube_size1*x[c] + n1[1]*tube_size1*y[c]-n2[1]*tube_size2;
    v[2] = v2[2] + n0[2]*tube_size1*x[c] + n1[2]*tube_size1*y[c]-n2[2]*tube_size2;
    glVertex3fv(v);
  }
  glEnd();
  
}


static void draw_string(float *v,char *l)
{
  glDisable(GL_DEPTH_TEST);	 
  glDisable(GL_LIGHTING);
  if(*l) {
    glColor3f(1.0,0.0,0.5);
    glRasterPos4f(v[0],v[1],v[2],1.0);
  }
  while(*l) {
    p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(l++));
  }
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);	 
}


/*========================================================================*/
void EditorRender(int state)
{
  CEditor *I = &Editor;
  int i0,i1;
  int sele1,sele2,sele3,sele4;
  float v0[3],v1[3],v2[3];

  if(I->Obj&&(state!=I->ActiveState))
    {
      EditorSetActiveObject(NULL,0);
    }
  
  if(I->Obj) {

    PRINTFD(FB_Editor)
      " EditorRender-Debug: rendering...\n"
      ENDFD;

    if(PMGUI) {
      
      sele1 = SelectorIndexByName(cEditorSele1);
      sele2 = SelectorIndexByName(cEditorSele2);
      sele3 = SelectorIndexByName(cEditorSele3);
      sele4 = SelectorIndexByName(cEditorSele4);

      if((sele1>=0)&&(sele2>=0)&&I->BondMode) {
        /* bond mode */

        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele1); /* slow */
        i1 = ObjectMoleculeGetAtomIndex(I->Obj,sele2); /* slow */
        if((i0>=0)&&(i1>=0)) {
          ObjectMoleculeGetAtomVertex(I->Obj,state,i0,v0);
          ObjectMoleculeGetAtomVertex(I->Obj,state,i1,v1);
          draw_bond(v0,v1);
        }

      } else {
        /* atom mode */
        
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele1); /* slow */
        if(i0>=0) {
          ObjectMoleculeGetAtomVertex(I->Obj,state,i0,v2);
          draw_globe(v2);
        }
        
        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele2); /* slow */
        if(i0>=0) {
          ObjectMoleculeGetAtomVertex(I->Obj,state,i0,v2);
          draw_globe(v2);
        }

        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele3); /* slow */
        if(i0>=0) {
          ObjectMoleculeGetAtomVertex(I->Obj,state,i0,v2);
          draw_globe(v2);
        }

        i0 = ObjectMoleculeGetAtomIndex(I->Obj,sele4); /* slow */
        if(i0>=0) {
          ObjectMoleculeGetAtomVertex(I->Obj,state,i0,v2);
          draw_globe(v2);
        }
      }
      if(I->ShowFrags) {
        int a;
        WordType buffer;
        for(a=0;a<I->NFrag;a++) {
          sprintf(buffer,"(%d)",a+1);
          draw_string(I->PosVLA+a*3,buffer);
        }
      }
    }
  }
}
/*========================================================================*/
void EditorInactive(void)
{
  CEditor *I = &Editor;

  PRINTFD(FB_Editor)
    " EditorInactive-Debug: callend.\n"
    ENDFD;

  I->Obj=NULL;
  I->BondMode = false;
  I->ShowFrags = false;
  SelectorDeletePrefixSet(cEditorFragPref);
  SelectorDeletePrefixSet(cEditorBasePref);
  ExecutiveDelete(cEditorSele1);      
  ExecutiveDelete(cEditorSele2);    
  ExecutiveDelete(cEditorSele3);    
  ExecutiveDelete(cEditorSele4);    
  ExecutiveDelete(cEditorSet);
  ExecutiveDelete(cEditorRes);
  ExecutiveDelete(cEditorChain);  
  ExecutiveDelete(cEditorObject);
  ExecutiveDelete(cEditorComp);
  ExecutiveDelete(cEditorLink);
  /*  if(SettingGet(cSetting_log_conformations)) PLogFlush();
      TODO: resolve this problem:
      we can't assume that Python interpreter isn't blocked
  */
  SceneDirty();
}
/*========================================================================*/
void EditorSetActiveObject(ObjectMolecule *obj,int state)
{
  int sele1,sele2,sele3,sele4;

  CEditor *I = &Editor;
  if(obj) {
    I->Obj=obj;
    sele1 = SelectorIndexByName(cEditorSele1);
    sele2 = SelectorIndexByName(cEditorSele2);
    sele3 = SelectorIndexByName(cEditorSele3);
    sele4 = SelectorIndexByName(cEditorSele4);

    if((sele1>=0)||(sele2>=0)||(sele3>=0)||(sele4>=0)) {
      
      ExecutiveDelete(cEditorComp);      
      ExecutiveDelete(cEditorRes);
      ExecutiveDelete(cEditorChain);
      ExecutiveDelete(cEditorObject);
      
      I->NFrag = SelectorSubdivideObject(cEditorFragPref,obj,
                                         sele1,sele2,
                                         sele3,sele4,
                                         cEditorBasePref,
                                         cEditorComp,
                                         &I->BondMode);
      I->ActiveState=state;
      
      if((I->NFrag>1)&&((int)SettingGet(cSetting_editor_label_fragments))) {
        SelectorComputeFragPos(obj,I->ActiveState,I->NFrag,cEditorFragPref,&I->PosVLA);
        I->ShowFrags = true;
      } else {
        I->ShowFrags = false;
      }
      if(SettingGet(cSetting_auto_hide_selections))
        ExecutiveHideSelections();
      
    } else {
      EditorInactive();
    }
  } else {
    I->NFrag = SelectorSubdivideObject(cEditorFragPref,NULL,
                                       -1,-1,-1,-1,
                                       cEditorBasePref,
                                       cEditorComp,
                                       &I->BondMode);
    EditorInactive();
  }
}
/*========================================================================*/
void EditorPrepareDrag(ObjectMolecule *obj,int index,int state)
{
  int frg;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int s;
  WordType name;
  int seleFlag= false;
  int i0,i1,i2,i3;
  CEditor *I = &Editor;
  int log_trans = (int)SettingGet(cSetting_log_conformations);
  int drag_sele = -1;
  PRINTFD(FB_Editor)
    " EditorPrepareDrag-Debug: entered. obj %p index %d",obj,index
    ENDFD;

  if(!I->Obj) { /* non-anchored */
    /* need to modify this code to move a complete covalent structure */

    I->DragObject=obj;
    I->DragIndex=index;
    I->DragSelection=-1;
  } else if(I->ActiveState!=state) {
    EditorSetActiveObject(NULL,0); /* recurses */
    return;
  } else { /* anchored */
    for(frg=1;frg<=I->NFrag;frg++) {
      sprintf(name,"%s%1d",cEditorFragPref,frg);
      drag_sele = SelectorIndexByName(name);
      if(drag_sele>=0) {
        s=obj->AtomInfo[index].selEntry;
        seleFlag=SelectorIsMember(s,drag_sele);
      }
      if(seleFlag) {
        strcpy(I->DragSeleName,name);
        break;
      }
    }
    if(seleFlag) { /* normal selection */
      
      PRINTFB(FB_Editor,FB_Actions)
        " Editor: grabbing (%s).",name
        ENDFB;

      I->DragIndex = index;
      I->DragSelection = drag_sele;
      I->DragObject = obj;
      I->DragHaveAxis = false;
      I->DragHaveBase = false;
      I->DragBondFlag = false;
      I->DragSlowFlag = false;
      
      sprintf(name,"%s%1d",cEditorBasePref,frg); /* get relevant base vertex of bond */
      sele1 = SelectorIndexByName(name);
      if(sele1>=0) {
        i1 = ObjectMoleculeGetAtomIndex(obj,sele1);
        if(i1>=0) {
          ObjectMoleculeGetAtomVertex(obj,state,i1,I->DragBase);
          I->DragHaveBase = true;
        }
      }

      /* get axis or base atom */

      {
        int cnt=0;

        if((sele0 = SelectorIndexByName(cEditorSele1))>=0) {
          if(SelectorIsAtomBondedToSele(obj,sele0,drag_sele))
            cnt++;
          else
            sele0 = -1;
        }
        if((sele1 = SelectorIndexByName(cEditorSele2))>=0) {
          if(SelectorIsAtomBondedToSele(obj,sele1,drag_sele))
            cnt++;
          else
            sele1 = -1;
        }
        if((sele2 = SelectorIndexByName(cEditorSele3))>=0) {
          if(SelectorIsAtomBondedToSele(obj,sele2,drag_sele))
            cnt++;
          else
            sele2 = -1;
        }
        if((sele3 = SelectorIndexByName(cEditorSele4))>=0) {
          if(SelectorIsAtomBondedToSele(obj,sele3,drag_sele))
            cnt++;
          else
            sele3 = -1;
        }
        
        i0 = ObjectMoleculeGetAtomIndex(obj,sele0);
        i1 = ObjectMoleculeGetAtomIndex(obj,sele1);
        i2 = ObjectMoleculeGetAtomIndex(obj,sele2);
        i3 = ObjectMoleculeGetAtomIndex(obj,sele3);
        
        if(cnt>1) { /* bond/multiatom mode */
          
          I->DragBondFlag = I->BondMode;
          
          zero3f(I->Center);
          if(i0>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i0,I->V0);          
          } else if(i1>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i1,I->V0);          
          } else if(i2>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i2,I->V0);          
          } else if(i3>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i3,I->V0);          
          }
          
          if(i0>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i0,I->V1);          
            add3f(I->V1,I->Center,I->Center);
          } 
          if(i1>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i1,I->V1);          
            add3f(I->V1,I->Center,I->Center);
          }
          if(i2>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i2,I->V1);          
            add3f(I->V1,I->Center,I->Center);
          }
          if(i3>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i3,I->V1);          
            add3f(I->V1,I->Center,I->Center);
          }
          
          {
            float div = 1.0F/cnt;
            scale3f(I->Center,div,I->Center);
          }
          
          subtract3f(I->Center,I->V0,I->Axis);
          
          normalize3f(I->Axis);
          I->DragHaveAxis=true;

        } else { /* atom mode */
          
          if(i0>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i0,I->V0);          
          } else if(i1>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i1,I->V0);          
          } else if(i2>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i2,I->V0);          
          } else if(i3>=0) {
            ObjectMoleculeGetAtomVertex(obj,state,i3,I->V0);          
          }
          if(I->DragHaveBase) {
            copy3f(I->DragBase,I->V1);
            subtract3f(I->V1,I->V0,I->Axis);
            average3f(I->V1,I->V0,I->Center);
            normalize3f(I->Axis);
            I->DragHaveAxis=true;
          }
        }
      }
    } else { /* clicked directly on an anchor atom */
      
      sele0=SelectorIndexByName(cEditorSele1);
      if(sele0<0) sele0=SelectorIndexByName(cEditorSele2);
      if(sele0<0) sele0=SelectorIndexByName(cEditorSele3);
      if(sele0<0) sele0=SelectorIndexByName(cEditorSele4);
      if(sele0>=0) {
        s=obj->AtomInfo[index].selEntry;
        seleFlag = SelectorIsMember(s,sele0);
      }

      PRINTFB(FB_Editor,FB_Actions)
        " Editor: grabbing all fragments." 
        ENDFB
      I->DragIndex = index;
      I->DragSelection = SelectorIndexByName(cEditorComp);
      strcpy(I->DragSeleName,cEditorComp);
      I->DragObject = obj;
      I->DragHaveAxis = false;
      I->DragHaveBase = false;
      I->DragBondFlag = false;

      I->DragSlowFlag = true;
      
      if(sele0>=0) { /* just provide a base vector, no valid axis exists */
        i1 = ObjectMoleculeGetAtomIndex(obj,sele0);
        if(i1>=0) {
          ObjectMoleculeGetAtomVertex(obj,state,i1,I->DragBase);
          I->DragHaveBase=true;
          I->DragBondFlag=true;
        }
      }
    } 
    if(!seleFlag) {
      I->DragIndex=-1;
      I->DragSelection=-1;
      I->DragObject=NULL;
    }
  }
  if(I->DragObject) {
    I->ShowFrags = false;
    ObjectMoleculeSaveUndo(I->DragObject,state,log_trans);
    if(SettingGet(cSetting_auto_sculpt)) {
      SettingSet(cSetting_sculpting,1);
      if(!I->DragObject->Sculpt)
        ObjectMoleculeSculptImprint(I->DragObject,state);
    }
  }
  if(log_trans) PLogFlush();

  PRINTFD(FB_Editor)
    " EditorPrepDrag-Debug: leaving Index %d Sele %d Object %p\n Axis %d Base %d BondFlag %d SlowFlag %d seleFlag %d\n",
    I->DragIndex,I->DragSelection,I->DragObject,I->DragHaveAxis,I->DragHaveBase,
    I->DragBondFlag,I->DragSlowFlag,seleFlag
    ENDFD;
}
/*========================================================================*/
void EditorDrag(ObjectMolecule *obj,int index,int mode,int state,
                float *pt,float *mov,float *z_dir)
{
  CEditor *I = &Editor;
  float v0[3],v1[3],v2[3],v3[3],v4[4],cp[3];
  float d0[3],d1[3],d2[3],n0[3],n1[3],n2[3];
  float opp,adj,theta;
  float m[16];
  int log_trans = (int)SettingGet(cSetting_log_conformations);

  PRINTFD(FB_Editor)
    " EditorDrag-Debug: entered. obj %p state %d index %d mode %d \nIndex %d Sele %d Object %p\n Axis %d Base %d BondFlag %d SlowFlag %d\n", obj,state,index,mode,
    I->DragIndex,I->DragSelection,I->DragObject,I->DragHaveAxis,I->DragHaveBase,
    I->DragBondFlag,I->DragSlowFlag
    ENDFD;

  if((index==I->DragIndex)&&(obj==I->DragObject)) {
    if(obj!=I->Obj) {
      /* non-achored actions */
      switch(mode) {
      case cButModeRotFrag:
        SceneOriginGet(v3);

        subtract3f(pt,v3,n0);
        add3f(pt,mov,n1);
        subtract3f(n1,v3,n1);
        normalize3f(n0);
        normalize3f(n1);
        cross_product3f(n0,n1,cp);
        theta = (float)asin(length3f(cp));
        normalize23f(cp,n2);        
        
        MatrixRotation44f(m,theta,n2[0],n2[1],n2[2]);
        m[3 ] = -v3[0];
        m[7 ] = -v3[1];
        m[11] = -v3[2];
        m[12] =  v3[0];
        m[13] =  v3[1];
        m[14] =  v3[2];
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m,log_trans,I->DragSeleName);
        SceneDirty();
        break;
      case cButModeTorFrag:
        ObjectMoleculeMoveAtom(obj,state,index,mov,1,log_trans);
        SceneDirty();
        break;
      case cButModeMovFrag:
        MatrixLoadIdentity44f(m);
        copy3f(mov,m+12);
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m,log_trans,I->DragSeleName);
        SceneDirty();
        break;
      }
    } else {
      switch(mode) {
      case cButModeRotFrag:
        if(I->DragHaveBase&&I->DragBondFlag) {
          copy3f(I->DragBase,v3);
        } else {
          copy3f(I->V0,v3);
        }
        if(I->DragSlowFlag) {
          SceneGetViewNormal(v4);
          scale3f(v4,-1.0F,v4);
          add3f(v3,v4,v4)
          subtract3f(pt,v4,n0);
          add3f(pt,mov,n1);
          subtract3f(n1,v4,n1);
        } else {
          subtract3f(pt,v3,n0);
          add3f(pt,mov,n1);
          subtract3f(n1,v3,n1);
        }
        normalize3f(n0);
        normalize3f(n1);
        cross_product3f(n0,n1,cp);
        theta = (float)asin(length3f(cp));
        normalize23f(cp,n2);        
        
        MatrixRotation44f(m,theta,n2[0],n2[1],n2[2]);
        m[3 ] = -v3[0];
        m[7 ] = -v3[1];
        m[11] = -v3[2];
        m[12] =  v3[0];
        m[13] =  v3[1];
        m[14] =  v3[2];
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m,log_trans,I->DragSeleName);
        SceneDirty();
        break;
      case cButModeTorFrag:
      case cButModePkTorBnd:
        if(I->DragHaveAxis) {
          subtract3f(pt,I->Center,d0);
          if(dot_product3f(d0,I->Axis)<0.0) {
            copy3f(I->V0,v1);
            copy3f(I->V1,v0);
          } else {
            copy3f(I->V0,v0);
            copy3f(I->V1,v1);
          }          
          subtract3f(v1,v0,d1);
          copy3f(d1,n0);
          normalize3f(n0);
          cross_product3f(n0,d0,n1);
          normalize3f(n1);
          
          project3f(d0,n0,v2);
          add3f(I->Center,v2,v2); /* v2 is the perpendicular point on the axis */
          subtract3f(pt,v2,d2);
          opp=(float)length3f(mov);
          adj=(float)length3f(d2);
          if(adj>R_SMALL4) {
            theta=(float)atan(opp/adj);
            if(dot_product3f(n1,mov)<0.0)
              theta=-theta;
            MatrixRotation44f(m,theta,n0[0],n0[1],n0[2]);
            m[3 ] = -v1[0];
            m[7 ] = -v1[1];
            m[11] = -v1[2];
            m[12] =  v1[0];
            m[13] =  v1[1];
            m[14] =  v1[2];
            ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m,log_trans,I->DragSeleName);
          } else {

            cross_product3f(I->Axis,z_dir,d0);
            theta = -dot_product3f(d0,mov);

            MatrixRotation44f(m,theta,n0[0],n0[1],n0[2]);
            m[3 ] = -v1[0];
            m[7 ] = -v1[1];
            m[11] = -v1[2];
            m[12] =  v1[0];
            m[13] =  v1[1];
            m[14] =  v1[2];
            ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m,log_trans,I->DragSeleName);
            
          }
        }
        SceneDirty();
        break;
      case cButModeMovFrag:
        MatrixLoadIdentity44f(m);
        copy3f(mov,m+12);
        ObjectMoleculeTransformSelection(obj,state,I->DragSelection,m,log_trans,I->DragSeleName);
        SceneDirty();
        break;
      }
    }
  }

  PRINTFD(FB_Editor)
    " EditorDrag-Debug: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void EditorInit(void)
{
  CEditor *I = &Editor;
  I->Obj = NULL;
  I->NFrag= 0;
  I->DragObject=NULL;
  I->DragIndex=-1;
  I->DragSelection=-1;
  I->NextPickSele = 0;
  I->BondMode = false;
  I->PosVLA = VLAlloc(float,30);
}
/*========================================================================*/
void EditorFree(void)
{
  CEditor *I = &Editor;
  VLAFreeP(I->PosVLA);
}


