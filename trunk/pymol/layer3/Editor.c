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
  int Active;
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

static int EditorGetEffectiveState(ObjectMolecule *obj,int state)
{
  if(!obj) obj = SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele1));
  if(!obj) obj = SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele2));
  if(!obj) obj = SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele3));
  if(!obj) obj = SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele4));

  if(obj) {
    if((obj->NCSet==1)&&(state>0))
      if(SettingGet_i(NULL,obj->Obj.Setting,cSetting_static_singletons))
        return 0;
  }
  return state;
}

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

int EditorDeselectIfSelected(ObjectMolecule *obj,int index,int update)
{
  CEditor *I = &Editor;
  int result = false;
  int s,sele;
  if(obj) {
    if((index>=0)&&(index<obj->NAtom)) {
      s=obj->AtomInfo[index].selEntry;                  
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
        EditorActivate(I->ActiveState,I->BondMode);
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
    result = PyList_New(3);
    PyList_SetItem(result,0,PyString_FromString(""));
    PyList_SetItem(result,1,PyInt_FromLong(I->ActiveState));
    PyList_SetItem(result,2,PyInt_FromLong(I->BondMode));
  }
  return(PConvAutoNone(result));
}

int EditorFromPyList(PyObject *list)
{
  int ok=true;
  int active_flag = false;
  int active_state;
  WordType obj_name;
  int ll = 0;
  int bond_mode = true;

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  if(ok) active_flag=(PyList_Size(list)!=0);
  if(!active_flag) {
    EditorInactivate();
  } else {
    if(ok) ok=PConvPyStrToStr(PyList_GetItem(list,0),obj_name,sizeof(WordType));
    if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,1),&active_state);
    if(ok&&(ll>2)) ok=PConvPyIntToInt(PyList_GetItem(list,2),&bond_mode); /* newer session files */
    if(ok) {
      EditorActivate(active_state,bond_mode);
      EditorDefineExtraPks();
    } else {
      EditorInactivate();
    }
  }
  if(!ok) {
    EditorInactivate();
  }
  return(ok);
}

int EditorActive(void) {
  CEditor *I = &Editor;
  return(I->Active);
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
  strcpy(name,cEditorSele4);
  I->NextPickSele=3;
  return;
  /*
  I->NextPickSele = (++I->NextPickSele)&0x3;
  switch(I->NextPickSele) {
  case 0: strcpy(name,cEditorSele1); break;
  case 1: strcpy(name,cEditorSele2); break;
  case 2: strcpy(name,cEditorSele3); break;
  case 3: strcpy(name,cEditorSele4); break;
  }
  return;    
  */
}

/*========================================================================*/
int EditorLogState(int pkresi)
{
  CEditor *I = &Editor;
  if(SettingGet(cSetting_logging)) {

    OrthoLineType buffer,buf1="None",buf2="None",buf3="None",buf4="None";
    int pkbond = 1;
    
    if(!EditorActive()) {
      PLog("edit",cPLog_pml);
    } else {
      int sele1,sele2,sele3,sele4;
      ObjectMolecule *obj1=NULL,*obj2=NULL,*obj3=NULL,*obj4=NULL;
      int index1,index2,index3,index4;
      
      sele1 = SelectorIndexByName(cEditorSele1);
      sele2 = SelectorIndexByName(cEditorSele2);
      sele3 = SelectorIndexByName(cEditorSele3);
      sele4 = SelectorIndexByName(cEditorSele4);

      obj1 = SelectorGetFastSingleAtomObjectIndex(sele1,&index1);
      obj2 = SelectorGetFastSingleAtomObjectIndex(sele2,&index2);
      obj3 = SelectorGetFastSingleAtomObjectIndex(sele3,&index3);
      obj4 = SelectorGetFastSingleAtomObjectIndex(sele4,&index4);

      if((sele1>=0) && (sele2>=0) && I->BondMode && obj1 && obj2) {

        /* bond mode */
        ObjectMoleculeGetAtomSeleLog(obj1,index1,buf1,true);
        ObjectMoleculeGetAtomSeleLog(obj2,index2,buf2,true);
        
      } else {

        /* atom mode */
        pkbond = 0;

        if(obj1) {
          ObjectMoleculeGetAtomSeleLog(obj1,index1,buf1,true);
        }

        if(obj2) {
          ObjectMoleculeGetAtomSeleLog(obj2,index2,buf2,true);
        }

        if(obj3) {
          ObjectMoleculeGetAtomSeleLog(obj3,index3,buf3,true);
        }

        if(obj4) {
          ObjectMoleculeGetAtomSeleLog(obj4,index4,buf4,true);
        }
      }
  
      sprintf(buffer,"cmd.edit(%s,%s,%s,%s,pkresi=%d,pkbond=%d)",
              buf1,buf2,buf3,buf4,pkresi ? 1: 0, pkbond ? 1: 0);

      PLog(buffer,cPLog_pym);

    }
  }
  return 1;
}
/*========================================================================*/

int EditorInvert(int quiet)
{
  CEditor *I = &Editor;
  int sele0,sele1,sele2;
  int i0,frg;
  int ia0=-1;
  int ia1=-1;
  float v[3],v0[3],v1[3];
  float n0[3],n1[3];
  float m[16];
  int state;
  int vf,vf0,vf1;
  int ok=false;
  int found=false;
  WordType name;
  ObjectMolecule *obj0,*obj1,*obj2;

  if(!EditorActive()) {
    ErrMessage("Editor","Must pick an atom to invert.");
  } else {
    sele0 = SelectorIndexByName(cEditorSele1);
    sele1 = SelectorIndexByName(cEditorSele2);
    sele2 = SelectorIndexByName(cEditorSele3);
    obj0 = SelectorGetFastSingleAtomObjectIndex(sele0,&i0);
    obj1 = SelectorGetFastSingleAtomObjectIndex(sele1,&ia0);
    obj2 = SelectorGetFastSingleAtomObjectIndex(sele2,&ia1);
    if(sele0<0) {
      ErrMessage("Editor","Must pick atom to invert as pk1.");        
    } else if(sele1<0) {
      ErrMessage("Editor","Must pick immobile atom in pk2.");
    } else if(sele2<0) {
      ErrMessage("Editor","Must pick immobile atom in pk3.");
    } else if(!(obj0&&(obj0==obj1)&&(obj0=obj2))) {
      ErrMessage("Editor","Must pick three atoms in the same object.");
    } else {

        state = SceneGetState();                
        ObjectMoleculeSaveUndo(obj0,state,false);
        
        vf  = ObjectMoleculeGetAtomVertex(obj0,state,i0,v);
        vf0 = ObjectMoleculeGetAtomVertex(obj0,state,ia0,v0);
        vf1 = ObjectMoleculeGetAtomVertex(obj0,state,ia1,v1);
        
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
          
          for(frg=1;frg<=I->NFrag;frg++) {
            sprintf(name,"%s%1d",cEditorFragPref,frg);
            sele2=SelectorIndexByName(name);
            
            if(ObjectMoleculeDoesAtomNeighborSele(obj0,i0,sele2) &&
               (!ObjectMoleculeDoesAtomNeighborSele(obj0,ia0,sele2)) &&
               (!ObjectMoleculeDoesAtomNeighborSele(obj0,ia1,sele2))) {
              found = true;
              ok = ObjectMoleculeTransformSelection(obj0,state,sele2,m,false,NULL);
            }
          }
          if(found) {
            if(!quiet) {
              PRINTFB(FB_Editor,FB_Actions) 
                " Editor: Inverted atom.\n"
                ENDFB;
            }
          } else {
            PRINTFB(FB_Editor,FB_Errors) 
              " Editor-Error: No free fragments found for inversion.\n"
              ENDFB;
          }

          SceneDirty();
          I->DragIndex=-1;
          I->DragSelection=-1;
          I->DragObject=NULL;
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
  ObjectMolecule *obj0=NULL,*obj1=NULL,*obj2=NULL;

  if(!EditorActive()) { 
    ErrMessage("Editor","Must specify a bond first.");
  } else {
    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      obj0 = SelectorGetFastSingleAtomObjectIndex(sele0,&i0);
      sele1 = SelectorIndexByName(cEditorSele2);
      obj1 = SelectorGetFastSingleAtomObjectIndex(sele1,&i1);
      strcpy(sele,cEditorFragPref);
      strcat(sele,"1");
      sele2 = SelectorIndexByName(sele);
      obj2 = SelectorGetFastSingleObjectMolecule(sele2);
      if(!((sele0>=0)&&(sele1>=0)&&(sele2>=0)&&(obj0==obj1))) {
        ErrMessage("Editor","Must specify a bond first.");
      } else {
        if((i0>=0)&&(i1>=0)) {
          state = SceneGetState();
          
          vf1 = ObjectMoleculeGetAtomVertex(obj0,state,i0,I->V0);
          vf2 = ObjectMoleculeGetAtomVertex(obj1,state,i1,I->V1);
          
          if(vf1&&vf2) {
            ObjectMoleculeSaveUndo(obj0,SceneGetState(),false);
            
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
            ok = ObjectMoleculeTransformSelection(obj2,state,sele2,m,false,NULL);
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
int EditorSelect(char *s0,char *s1,char *s2,char *s3,int pkresi,int pkbond,int quiet)
{
  int i0=-1;
  int i1=-1;
  int i2=-1;
  int i3=-1;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int result=false;
  int ok = true;
  ObjectMolecule *obj0=NULL,*obj1=NULL,*obj2=NULL,*obj3=NULL;

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
    obj0 = SelectorGetFastSingleAtomObjectIndex(sele0,&i0);
    ExecutiveDelete(cEditorSele1);
  }
 
  if(s1) {
    sele1 = SelectorIndexByName(s1);
    obj1 = SelectorGetFastSingleAtomObjectIndex(sele1,&i1);
    ExecutiveDelete(cEditorSele2);
  }

  if(s2) {
    sele2 = SelectorIndexByName(s2);
    obj2 = SelectorGetFastSingleAtomObjectIndex(sele2,&i2);
    ExecutiveDelete(cEditorSele3);
  }

  if(s3) {
    sele3 = SelectorIndexByName(s3);
    obj3 = SelectorGetFastSingleAtomObjectIndex(sele3,&i3);
    ExecutiveDelete(cEditorSele4);
  }

  if(!(obj0||obj1||obj2||obj3)) 
    ok = false;

  if(ok) {
    if(obj0)
      ObjectMoleculeVerifyChemistry(obj0);  
    if(obj1&&(obj1!=obj0))
      ObjectMoleculeVerifyChemistry(obj1);  
    if(obj2&&(obj2!=obj0)&&(obj2!=obj1))
      ObjectMoleculeVerifyChemistry(obj2);  
    if(obj3&&(obj3!=obj0)&&(obj3!=obj1)&&(obj3!=obj2))
      ObjectMoleculeVerifyChemistry(obj3);  

    if(i0>=0) SelectorCreate(cEditorSele1,s0,NULL,quiet,NULL);
    if(i1>=0) SelectorCreate(cEditorSele2,s1,NULL,quiet,NULL);
    if(i2>=0) SelectorCreate(cEditorSele3,s2,NULL,quiet,NULL);
    if(i3>=0) SelectorCreate(cEditorSele4,s3,NULL,quiet,NULL);
    
    EditorActivate(SceneGetState(),pkbond);        
    
    if(pkresi)
      EditorDefineExtraPks();
    
    SceneDirty();
    result=true;

  } else {
    EditorInactivate();
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
int EditorIsAnActiveObject(ObjectMolecule *obj) {
  if(EditorActive()) {
    if(obj) {
      if(obj == SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele1)))
        return true;
      if(obj == SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele2)))
        return true;
      if(obj == SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele3)))
        return true;
      if(obj == SelectorGetFastSingleObjectMolecule(SelectorIndexByName(cEditorSele4)))
        return true;
    }
  }
  return false;
}
/*========================================================================*/
void EditorCycleValence(int quiet)
{
  CEditor *I = &Editor;
  int sele0,sele1;

  if(EditorActive()) {
    ObjectMolecule *obj0,*obj1;
    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      if(sele1>=0) {
        obj0 = SelectorGetFastSingleObjectMolecule(sele0);
        obj1 = SelectorGetFastSingleObjectMolecule(sele1);
        if((obj0==obj1)&&I->BondMode) {
          ObjectMoleculeVerifyChemistry(obj0);
          ObjectMoleculeAdjustBonds(obj0,sele0,sele1,0,0);
        }
      }
    }
  }
}

/*========================================================================*/
void EditorAttach(char *elem,int geom,int valence,char *name,int quiet)
{
  int i0;
  int sele0,sele1;
  int state;
  AtomInfoType *ai;
  ObjectMolecule *obj0=NULL,*obj1=NULL;

  ai=(AtomInfoType*)VLAMalloc(1,sizeof(AtomInfoType),1,true);
  if(EditorActive()) {

    sele0 = SelectorIndexByName(cEditorSele1);
    if(sele0>=0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      obj0 = SelectorGetFastSingleObjectMolecule(sele0);
      obj1 = SelectorGetFastSingleObjectMolecule(sele1);

      if(obj0) {
        ObjectMoleculeVerifyChemistry(obj0); /* remember chemistry for later */
        state = SceneGetState();
        if(obj1) {
          if(obj0==obj1) {
            /* bond mode - behave like replace */
            EditorReplace(elem,geom,valence,name,quiet);
          }
        } else {
          /* atom mode */
          i0 = ObjectMoleculeGetAtomIndex(obj0,sele0); /* slow */
          if(i0>=0) {
            UtilNCopy(ai->elem,elem,sizeof(AtomName));
            ai->geom=geom;
            ai->valence=valence;
            if(name[0])
              UtilNCopy(ai->name,name,sizeof(AtomName));
            ObjectMoleculeAttach(obj0,i0,ai); /* will free ai */
            ai = NULL;
          }
        }
      }
    }
  }
  VLAFreeP(ai); /* safety */
}
/*========================================================================*/
void EditorRemove(int hydrogen,int quiet)
{
  int sele0,sele1;
  int i0;
  int h_flag = false;
  OrthoLineType buf;
  ObjectMolecule *obj0=NULL,*obj1=NULL;

#define cEditorRemoveSele "_EditorRemove"

  if(EditorActive()) {
    sele0 = SelectorIndexByName(cEditorSele1);
    obj0 = SelectorGetFastSingleObjectMolecule(sele0);
    ObjectMoleculeVerifyChemistry(obj0); /* remember chemistry for later */
    if((sele0>=0)&&obj0) {
      sele1 = SelectorIndexByName(cEditorSele2);
      obj1 = SelectorGetFastSingleObjectMolecule(sele1);
      if((sele1>=0)&&(obj0==obj1)) {
        /* bond mode */
        ObjectMoleculeRemoveBonds(obj0,sele0,sele1);
        EditorInactivate();
      } else {

        if(hydrogen) {
          sprintf(buf,"((neighbor %s) and hydro)",cEditorSele1);          
          h_flag = SelectorCreate(cEditorRemoveSele,buf,NULL,false,NULL);
        }

        if(SelectorGetFastSingleAtomObjectIndex(sele0,&i0)) {
          /* atom mode */
          if(i0>=0) {
            ExecutiveRemoveAtoms(cEditorSele1,quiet);
            EditorInactivate();
          }
        }

        if(h_flag) {
          ExecutiveRemoveAtoms(cEditorRemoveSele,quiet);
          SelectorDelete(cEditorRemoveSele);
        }

      }
    }
  }
#undef cEditorRemoveSele
}
/*========================================================================*/
void EditorHFill(int quiet)
{
  int sele0,sele1;
  int i0;
  OrthoLineType buffer,s1;
  ObjectMolecule *obj0=NULL,*obj1=NULL;

  if(EditorActive()) {
    sele0 = SelectorIndexByName(cEditorSele1);
    obj0 = SelectorGetFastSingleObjectMolecule(sele0);    
    ObjectMoleculeVerifyChemistry(obj0); /* remember chemistry for later */
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
        ExecutiveRemoveAtoms(s1,quiet);
        SelectorFreeTmp(s1);
        i0 = ObjectMoleculeGetAtomIndex(obj0,sele0); 
        obj0->AtomInfo[i0].chemFlag=false;
        ExecutiveAddHydrogens(cEditorSele1,quiet);
      }
      
      if(sele1>=0) {
        obj1 = SelectorGetFastSingleObjectMolecule(sele1);    
        if(sele0>=0) 
          sprintf(buffer,"((neighbor %s) and (elem h) and not %s)",
                  cEditorSele2,cEditorSele1);
        else 
          sprintf(buffer,"((neighbor %s) and (elem h))",
                  cEditorSele2);
        SelectorGetTmp(buffer,s1);
        ExecutiveRemoveAtoms(s1,quiet);
        SelectorFreeTmp(s1);
        i0 = ObjectMoleculeGetAtomIndex(obj1,sele1); 
        obj1->AtomInfo[i0].chemFlag=false;
        ExecutiveAddHydrogens(cEditorSele2,quiet);
      }
    }
  }
  
}
/*========================================================================*/
void EditorReplace(char *elem,int geom,int valence,char *name,int quiet)
{
  int i0;
  int sele0;
  int state;
  AtomInfoType ai;
  ObjectMolecule *obj0=NULL;

  UtilZeroMem(&ai,sizeof(AtomInfoType));
  if(EditorActive()) {
    sele0 = SelectorIndexByName(cEditorSele1);
    obj0 = SelectorGetFastSingleObjectMolecule(sele0);    
    ObjectMoleculeVerifyChemistry(obj0); /* remember chemistry for later */

    state = SceneGetState();

    if(sele0>=0) {
      i0 = ObjectMoleculeGetAtomIndex(obj0,sele0); /* slow */
      if(i0>=0) {
        UtilNCopy(ai.elem,elem,sizeof(AtomName));
        if(name[0])
          UtilNCopy(ai.name,name,sizeof(AtomName));
        ai.geom=geom;
        ai.valence=valence;
        ObjectMoleculePrepareAtom(obj0,i0,&ai);
        ObjectMoleculePreposReplAtom(obj0,i0,&ai);
        ObjectMoleculeReplaceAtom(obj0,i0,&ai); /* invalidates */
        ObjectMoleculeVerifyChemistry(obj0);
        ObjectMoleculeFillOpenValences(obj0,i0);
        ObjectMoleculeSort(obj0);
        ObjectMoleculeUpdateIDNumbers(obj0);
        EditorInactivate();
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

static void draw_globe(float *v2,int number)
{
  float v[3];
  float n0[3],n1[3],n2[3];
  float x[50],y[50];
  int nEdge;
  int a,c;
  float radius=0.5F;
  float width_base=0.10F;
  float width=0.0F;
  float offset = 0.0F;
  int cycle_counter;

  nEdge = (int)SettingGet(cSetting_stick_quality)*2;
  if(nEdge>50)
    nEdge=50;
  
  subdivide(nEdge,x,y);
    
  n0[0]=1.0;
  n0[1]=0.0;
  n0[2]=0.0;
  get_system1f3f(n0,n1,n2);
  
  glColor3fv(ColorGet(0));

  cycle_counter = number;
  while(cycle_counter) {
    
    switch(number) {
    case 1:
      width = width_base;
      offset = 0.0F;
      break;

    case 2:
      switch(cycle_counter) {
      case 2:
        width = width_base/2;
        offset = width_base;
        break;
      case 1:
        offset = -width_base;
        break;
      }
      break;

    case 3:
      switch(cycle_counter) {
      case 3:
        width = width_base/2.8;
        offset = 1.33*width_base;
        break;
      case 2:
        offset = 0.0F;
        break;
      case 1:
        offset = -1.33*width_base;
        break;
      }
      break;

    case 4:
      switch(cycle_counter) {
      case 4:
        width = width_base/3.2;
        offset = 2*width_base;
        break;
      case 3:
        offset = 0.66F*width_base;
        break;
      case 2:
        offset = -0.66F*width_base;
        break;
      case 1:
        offset = -2*width_base;
        break;
      }
    }

    glBegin(GL_TRIANGLE_STRIP);
    for(a=0;a<=nEdge;a++) {
      c=a % nEdge;
      v[0] =  n1[0]*x[c] + n2[0]*y[c];
      v[1] =  n1[1]*x[c] + n2[1]*y[c];
      v[2] =  n1[2]*x[c] + n2[2]*y[c];
      normalize3f(v);
      glNormal3fv(v);
      v[0] = v2[0] + n1[0]*radius*x[c] + n2[0]*radius*y[c]+n0[0]*(offset+width);
      v[1] = v2[1] + n1[1]*radius*x[c] + n2[1]*radius*y[c]+n0[1]*(offset+width);
      v[2] = v2[2] + n1[2]*radius*x[c] + n2[2]*radius*y[c]+n0[2]*(offset+width);
      glVertex3fv(v);
      v[0] = v2[0] + n1[0]*radius*x[c] + n2[0]*radius*y[c]+n0[0]*(offset-width);
      v[1] = v2[1] + n1[1]*radius*x[c] + n2[1]*radius*y[c]+n0[1]*(offset-width);
      v[2] = v2[2] + n1[2]*radius*x[c] + n2[2]*radius*y[c]+n0[2]*(offset-width);
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
      v[0] = v2[0] + n0[0]*radius*x[c] + n2[0]*radius*y[c]+n1[0]*(offset+width);
      v[1] = v2[1] + n0[1]*radius*x[c] + n2[1]*radius*y[c]+n1[1]*(offset+width);
      v[2] = v2[2] + n0[2]*radius*x[c] + n2[2]*radius*y[c]+n1[2]*(offset+width);
      glVertex3fv(v);
      v[0] = v2[0] + n0[0]*radius*x[c] + n2[0]*radius*y[c]+n1[0]*(offset-width);
      v[1] = v2[1] + n0[1]*radius*x[c] + n2[1]*radius*y[c]+n1[1]*(offset-width);
      v[2] = v2[2] + n0[2]*radius*x[c] + n2[2]*radius*y[c]+n1[2]*(offset-width);
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
      v[0] = v2[0] + n0[0]*radius*x[c] + n1[0]*radius*y[c]+n2[0]*(offset+width);
      v[1] = v2[1] + n0[1]*radius*x[c] + n1[1]*radius*y[c]+n2[1]*(offset+width);
      v[2] = v2[2] + n0[2]*radius*x[c] + n1[2]*radius*y[c]+n2[2]*(offset+width);
      glVertex3fv(v);
      v[0] = v2[0] + n0[0]*radius*x[c] + n1[0]*radius*y[c]+n2[0]*(offset-width);
      v[1] = v2[1] + n0[1]*radius*x[c] + n1[1]*radius*y[c]+n2[1]*(offset-width);
      v[2] = v2[2] + n0[2]*radius*x[c] + n1[2]*radius*y[c]+n2[2]*(offset-width);
      glVertex3fv(v);
    }
    glEnd();

    cycle_counter--;
  }
  
}

/*
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
*/

/*========================================================================*/
void EditorRender(int state)
{
  CEditor *I=&Editor;
  int sele1,sele2,sele3,sele4;
  float v0[3],v1[3],v2[3];
  ObjectMolecule *obj1=NULL,*obj2=NULL,*obj3=NULL,*obj4=NULL;
  int index1,index2,index3,index4;

  if(EditorActive()) {

    PRINTFD(FB_Editor)
      " EditorRender-Debug: rendering...\n"
      ENDFD;

    if(PMGUI) {
      
      sele1 = SelectorIndexByName(cEditorSele1);
      sele2 = SelectorIndexByName(cEditorSele2);
      sele3 = SelectorIndexByName(cEditorSele3);
      sele4 = SelectorIndexByName(cEditorSele4);

      obj1 = SelectorGetFastSingleAtomObjectIndex(sele1,&index1);
      obj2 = SelectorGetFastSingleAtomObjectIndex(sele2,&index2);
      obj3 = SelectorGetFastSingleAtomObjectIndex(sele3,&index3);
      obj4 = SelectorGetFastSingleAtomObjectIndex(sele4,&index4);

      /*      printf("%d %d %d %d\n",sele1,sele2,sele3,sele4);
      printf("%p %p %p %p\n",obj1,obj2,obj3,obj4);
      printf("%d %d %d %d\n",index1,index2,index3,index4);*/

      if((sele1>=0) && (sele2>=0) && I->BondMode && obj1 && obj2) {
        /* bond mode */

        ObjectMoleculeGetAtomVertex(obj1,state,index1,v0);
        ObjectMoleculeGetAtomVertex(obj2,state,index2,v1);
        draw_bond(v0,v1);
        
      } else {
        /* atom mode */
        
        if(obj1) {
          if(ObjectMoleculeGetAtomVertex(obj1,state,index1,v2))
            draw_globe(v2,1);
        }

        if(obj2) {
          if(ObjectMoleculeGetAtomVertex(obj2,state,index2,v2))
            draw_globe(v2,2);
        }

        if(obj3) {
          if(ObjectMoleculeGetAtomVertex(obj3,state,index3,v2))
            draw_globe(v2,3);
        }

        if(obj4) {
          if(ObjectMoleculeGetAtomVertex(obj4,state,index4,v2))
            draw_globe(v2,4);
        }
      }
      /*
        if(I->ShowFrags) {
        int a;
         WordType buffer;
        for(a=0;a<I->NFrag;a++) {
        sprintf(buffer,"(%d)",a+1);
        draw_string(I->PosVLA+a*3,buffer);
        }
      }
      */
    }
  }
}
/*========================================================================*/
void EditorInactivate(void)
{
  CEditor *I = &Editor;

  PRINTFD(FB_Editor)
    " EditorInactivate-Debug: callend.\n"
    ENDFD;

  I->BondMode = false;
  I->ShowFrags = false;
  I->NFrag = 0;
  I->Active = false;
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
void EditorActivate(int state,int enable_bond)
{
  int sele1,sele2,sele3,sele4;
  
  CEditor *I = &Editor;

  sele1 = SelectorIndexByName(cEditorSele1);
  sele2 = SelectorIndexByName(cEditorSele2);
  sele3 = SelectorIndexByName(cEditorSele3);
  sele4 = SelectorIndexByName(cEditorSele4);
  
  if((sele1>=0)||(sele2>=0)||(sele3>=0)||(sele4>=0)) {
    
    I->Active = true;
    ExecutiveDelete(cEditorComp);      
    ExecutiveDelete(cEditorRes);
    ExecutiveDelete(cEditorChain);
    ExecutiveDelete(cEditorObject);
    
    I->BondMode = enable_bond;
    I->NFrag = SelectorSubdivide(cEditorFragPref,
                                       sele1,sele2,
                                       sele3,sele4,
                                       cEditorBasePref,
                                       cEditorComp,
                                       &I->BondMode);
    
    state = EditorGetEffectiveState(NULL,state);
    I->ActiveState=state;
    
    if(0&&(I->NFrag>1)&&((int)SettingGet(cSetting_editor_label_fragments))) {
      /*      SelectorComputeFragPos(obj,I->ActiveState,I->NFrag,cEditorFragPref,&I->PosVLA);*/
      I->ShowFrags = true;
    } else {
      I->ShowFrags = false;
    }
    if(SettingGet(cSetting_auto_hide_selections))
      ExecutiveHideSelections();
    
  } else {
    EditorInactivate();
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
    " EditorPrepareDrag-Debug: entered. obj %p index %d",
    (void*)obj,index
    ENDFD;

  state = EditorGetEffectiveState(obj,state);

  if(!EditorActive()) { /* non-anchored */
    /* need to modify this code to move a complete covalent structure */

    I->DragObject=obj;
    I->DragIndex=index;
    I->DragSelection=-1;
  } else {


    /* anchored */
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
      
      PRINTFB(FB_Editor,FB_Blather)
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
    I->DragIndex,I->DragSelection,(void*)I->DragObject,
    I->DragHaveAxis,I->DragHaveBase,
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
    " EditorDrag-Debug: entered. obj %p state %d index %d mode %d \nIndex %d Sele %d Object %p\n Axis %d Base %d BondFlag %d SlowFlag %d\n", 
    (void*)obj,state,index,mode,
    I->DragIndex,I->DragSelection,(void*)I->DragObject,
    I->DragHaveAxis,I->DragHaveBase,
    I->DragBondFlag,I->DragSlowFlag
    ENDFD;

  state = EditorGetEffectiveState(obj,state);

  if((index==I->DragIndex)&&(obj==I->DragObject)) {
    if(!EditorActive()) {
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
  I->Active = false;
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


