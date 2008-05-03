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

#include"Base.h"
#include"OOMac.h"
#include"Vector.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Setting.h"
#include"Scene.h"
#include"Ray.h"
#include"ObjectDist.h"
#include"Selector.h"
#include"PConv.h"
#include"ObjectMolecule.h"
#include"Feedback.h"

void ObjectDistFree(ObjectDist *I);
void ObjectDistUpdate(ObjectDist *I);
int ObjectDistGetNFrames(ObjectDist *I);
void ObjectDistUpdateExtents(ObjectDist *I);

int ObjectDistGetLabelTxfVertex(ObjectDist *I,int state,int index,float *v)
{
  int result = 0;
  if(I->DSet) { 
    if(state<0) state=SettingGet_i(I->Obj.G,NULL,I->Obj.Setting,cSetting_state)-1;
    if(state<0) state=SceneGetState(I->Obj.G); 
    if(I->NDSet==1) state=0; /* static singletons always active here it seems */
    state = state % I->NDSet;
    {
	DistSet *ds = I->DSet[state];
    if((!ds)&&(SettingGet_b(I->Obj.G,I->Obj.Setting,NULL,cSetting_all_states))) {
      state=0;
      ds = I->DSet[state];
    }
    if(ds) {
      result = DistSetGetLabelVertex(ds,index,v);
    }}
  }
  return(result);
}

int ObjectDistMoveLabel(ObjectDist *I,int state,int index,float *v,int mode,int log)
{
  int result = 0;
  DistSet *ds;
  if(state<0) state=0;
  if(I->NDSet==1) state=0;
  state = state % I->NDSet;
  if((!I->DSet[state])&&(SettingGet_b(I->Obj.G,I->Obj.Setting,NULL,cSetting_all_states)))
    state=0;
  ds = I->DSet[state];
  if(ds) {
    result = DistSetMoveLabel(I->DSet[state],index,v,mode);
    ds->fInvalidateRep(ds,cRepLabel,cRepInvCoord);
    /*      ExecutiveUpdateCoordDepends(I->Obj.G,I);*/
  }
#if 0
  if(log) {
    OrthoLineType line,buffer;
    if(SettingGet(I->Obj.G,cSetting_logging)) {
      ObjectMoleculeGetAtomSele(I,index,buffer);
      sprintf(line,"cmd.translate_atom(\"%s\",%15.9f,%15.9f,%15.9f,%d,%d,%d)\n",
              buffer,v[0],v[1],v[2],state+1,mode,0);
      PLog(G,line,cPLog_no_flush);
    }
  }
#endif
  return(result);
}

static DistSet *ObjectDistGetDistSetFromM4XBond(PyMOLGlobals *G,
                                                ObjectMolecule *obj, 
                                                M4XBondType *hb, int n_hb,
                                                int state,int nbr_sele)
{
  int min_id,max_id,range,*lookup = NULL;
  int nv = 0;
  float *vv=NULL;
  DistSet *ds;
  ds = DistSetNew(G);
  vv = VLAlloc(float,10);
  
  /* this routine only works if IDs cover a reasonable range --
     should rewrite using a hash table */
  
  if(obj->NAtom) {
    
    /* determine range */

    {
      int a,cur_id;
      cur_id = obj->AtomInfo[0].id;
      min_id = cur_id;
      max_id = cur_id;
      for(a=1;a<obj->NAtom;a++) {
        cur_id = obj->AtomInfo[a].id;
        if(min_id>cur_id) min_id = cur_id;
        if(max_id<cur_id) max_id = cur_id;
      }
    }

    /* create cross-reference table */

    {
      int a,offset;
      
      range = max_id - min_id + 1;
      lookup = Calloc(int,range);
      for(a=0;a<obj->NAtom;a++) {
        offset = obj->AtomInfo[a].id - min_id;
        if(lookup[offset])
          lookup[offset] = -1;
        else {
          lookup[offset] = a+1;
        }
      }
    }
    
    /* iterate through IDs and get pairs */
    {
      AtomInfoType *ai1,*ai2;
      int at1,at2;
      CoordSet *cs;
      
      float *vv0,*vv1,dist;
      int idx1,idx2;
      
      int i,offset1,offset2;
      int sele_flag = false;
      
      for(i=0;i<n_hb;i++) {
        offset1 = hb[i].atom1-min_id;
        offset2 = hb[i].atom2-min_id;
        if((offset1>=0)&&(offset1<range)&&(offset2>=0)&&(offset2<range)) {
          at1 = lookup[offset1]-1;
          at2 = lookup[offset2]-1;
          if((at1>=0)&&(at2>=0)&&(at1!=at2)&&(state<obj->NCSet)) {
            cs=obj->CSet[state];

            if(cs) { 
              ai1=obj->AtomInfo+at1;
              ai2=obj->AtomInfo+at2;
              
              sele_flag =false;
              
              if(nbr_sele>=0) {
                if(SelectorIsMember(G,ai1->selEntry,nbr_sele))
                  sele_flag=true;
                if(SelectorIsMember(G,ai2->selEntry,nbr_sele))
                  sele_flag=true;
              } else {
                sele_flag =true;
              }
              
              if(sele_flag) {
                if(obj->DiscreteFlag) {
                  if(cs==obj->DiscreteCSet[at1]) {
                    idx1=obj->DiscreteAtmToIdx[at1];
                  } else {
                    idx1=-1;
                  }
                } else {
                  idx1=cs->AtmToIdx[at1];
                }
                
                if(obj->DiscreteFlag) {
                  if(cs==obj->DiscreteCSet[at2]) {
                    idx2=obj->DiscreteAtmToIdx[at2];
                  } else {
                    idx2=-1;
                  }
                  
                } else {
                  idx2=cs->AtmToIdx[at2];
                }
                
                if((idx1>=0)&&(idx2>=0)) {
                  dist=(float)diff3f(cs->Coord+3*idx1,cs->Coord+3*idx2);
                  VLACheck(vv,float,(nv*3)+5);
                  vv0 = vv+ (nv*3);
                  vv1 = cs->Coord+3*idx1;
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  vv1 = cs->Coord+3*idx2;
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  nv+=2;
                }
              }
            }
          }
        }
      }
    }
  }
  
  FreeP(lookup);
  ds->NIndex = nv;
  ds->Coord = vv;
  return(ds);
  
}

ObjectDist *ObjectDistNewFromM4XBond(PyMOLGlobals *G,ObjectDist *oldObj,
                                      struct ObjectMolecule *objMol,
                                      struct M4XBondType *hbond,
                                     int n_hbond,int nbr_sele)
{
  int a;
  ObjectDist *I;
  int n_state;
  if(!oldObj)
    I=ObjectDistNew(G);
  else {
    I=oldObj;
    for(a=0;a<I->NDSet;a++)
      if(I->DSet[a]) {
        if(I->DSet[a]->fFree)
          I->DSet[a]->fFree(I->DSet[a]);
        I->DSet[a]=NULL;
      }
    I->NDSet=0;
  }
  n_state = objMol->NCSet;
  for(a=0;a<n_state;a++) {
    VLACheck(I->DSet,DistSet*,a);
    
    I->DSet[a] = ObjectDistGetDistSetFromM4XBond(G,objMol,hbond,n_hbond,a,nbr_sele);
    
    if(I->DSet[a]) {
      I->DSet[a]->Obj = I;
        I->NDSet=a+1;
    }
  } 
  ObjectDistUpdateExtents(I);
  
  SceneChanged(G);
  return(I);
}

/*========================================================================*/

void ObjectDistUpdateExtents(ObjectDist *I)
{
  float maxv[3] = {FLT_MAX,FLT_MAX,FLT_MAX};
  float minv[3] = {-FLT_MAX,-FLT_MAX,-FLT_MAX};
  int a;
  DistSet *ds;

  /* update extents */
  copy3f(maxv,I->Obj.ExtentMin);
  copy3f(minv,I->Obj.ExtentMax);
  I->Obj.ExtentFlag=false;
  for(a=0;a<I->NDSet;a++) {
    ds = I->DSet[a];
    if(ds) {
      if(DistSetGetExtent(ds,I->Obj.ExtentMin,I->Obj.ExtentMax))
        I->Obj.ExtentFlag=true;
    }
  }
}


#ifndef _PYMOL_NOPY
static PyObject *ObjectDistDSetAsPyList(ObjectDist *I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NDSet);
  for(a=0;a<I->NDSet;a++) {
    if(I->DSet[a]) {
      PyList_SetItem(result,a,DistSetAsPyList(I->DSet[a]));
    } else {
      PyList_SetItem(result,a,PConvAutoNone(Py_None));
    }
  }
  return(PConvAutoNone(result));
}
#endif
#ifndef _PYMOL_NOPY
static int ObjectDistDSetFromPyList(ObjectDist *I,PyObject *list)
{
  int ok=true;
  int a;
  if(ok) ok=PyList_Check(list);
  if(ok) {
    VLACheck(I->DSet,DistSet*,I->NDSet);
    for(a=0;a<I->NDSet;a++) {
      if(ok) ok = DistSetFromPyList(I->Obj.G,PyList_GetItem(list,a),&I->DSet[a]);
      if(ok) I->DSet[a]->Obj = I;
    }
  }
  return(ok);
}
#endif
/*========================================================================*/
PyObject *ObjectDistAsPyList(ObjectDist *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  /* first, dump the atoms */

  result = PyList_New(4);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NDSet));
  PyList_SetItem(result,2,ObjectDistDSetAsPyList(I));
  PyList_SetItem(result,3,PyInt_FromLong(I->CurDSet));

#if 0

  CObject Obj;
  struct DistSet **DSet;
  int NDSet;
  int CurDSet;

#endif

  return(PConvAutoNone(result));  
#endif
}

int ObjectDistNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectDist **result)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  ObjectDist *I=NULL;
  (*result) = NULL;
  
  if(ok) ok=PyList_Check(list);

  I=ObjectDistNew(G);
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(G,PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NDSet);
  if(ok) ok = ObjectDistDSetFromPyList(I,PyList_GetItem(list,2));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,3),&I->CurDSet);
  
  ObjectDistInvalidateRep(I,cRepAll);
  if(ok) {
    (*result) = I;
    ObjectDistUpdateExtents(I);
  }
  else {
    /* cleanup? */
  }

  return(ok);
#endif
}


/*========================================================================*/
int ObjectDistGetNFrames(ObjectDist *I)
{
  return I->NDSet;
}
/*========================================================================*/
void ObjectDistUpdate(ObjectDist *I)
{
  int a;
  OrthoBusyPrime(I->Obj.G);
  for(a=0;a<I->NDSet;a++)
	 if(I->DSet[a]) {	
	   OrthoBusySlow(I->Obj.G,a,I->NDSet);
      /*	   printf(" ObjectDist: updating state %d of \"%s\".\n" , a+1, I->Obj.Name);*/
      if(I->DSet[a]->fUpdate)
        I->DSet[a]->fUpdate(I->DSet[a],a);
	 }
}
/*========================================================================*/
void ObjectDistInvalidateRep(ObjectDist *I,int rep)
{
  int a;
  PRINTFD(I->Obj.G,FB_ObjectDist)
    " ObjectDistInvalidateRep: entered.\n"
    ENDFD;

  for(a=0;a<I->NDSet;a++) 
	 if(I->DSet[a]) {	 
      if(I->DSet[a]->fInvalidateRep)
        I->DSet[a]->fInvalidateRep(I->DSet[a],rep,cRepInvAll);
	 }
}
/*========================================================================*/
static void ObjectDistRender(ObjectDist *I,RenderInfo *info)
{
  int a;
  int state = info->state;
  int pass = info->pass;
  CRay *ray = info->ray;

  if((pass==0)||(pass==-1)) {
    ObjectPrepareContext(&I->Obj,ray);
    if(state<0) {
      for(a=0;a<I->NDSet;a++)
        if(I->DSet[a])
          if(I->DSet[a]->fRender)
            I->DSet[a]->fRender(I->DSet[a],info);
    } else if(state<I->NDSet) {
      I->CurDSet=state % I->NDSet;
      if(I->DSet[I->CurDSet]) {
        if(I->DSet[I->CurDSet]->fRender)
          I->DSet[I->CurDSet]->fRender(I->DSet[I->CurDSet],info);
      }
    } else if(I->NDSet==1) { /* if only one coordinate set, assume static */
      if(I->DSet[0]->fRender) {
        if(SettingGet_b(I->Obj.G,I->Obj.Setting,NULL,cSetting_static_singletons)) {
          I->DSet[0]->fRender(I->DSet[0],info);
        }
      }
    }
  }
}

static CSetting **ObjectDistGetSettingHandle(ObjectDist *I,int state)
{
  if(state<0) {
    return(&I->Obj.Setting);
  } else if(state<I->NDSet) {
    if(I->DSet[state]) {
      return(&I->DSet[state]->Setting);
    } else {
      return(NULL);
    }
  } else {
    return(NULL);
  }
}

/*========================================================================*/
ObjectDist *ObjectDistNew(PyMOLGlobals *G)
{
  OOAlloc(G,ObjectDist);
  ObjectInit(G,(CObject*)I);
  I->Obj.type=cObjectMeasurement;
  I->DSet=VLAMalloc(10,sizeof(DistSet*),5,true); /* auto-zero */
  I->NDSet=0;
  I->Obj.fRender=(void (*)(CObject *, RenderInfo *info))ObjectDistRender;
  I->Obj.fFree= (void (*)(CObject *))ObjectDistFree;
  I->Obj.fUpdate= (void (*)(CObject *)) ObjectDistUpdate;
  I->Obj.fGetNFrame = (int (*)(CObject *)) ObjectDistGetNFrames;
  I->Obj.fGetSettingHandle = (CSetting **(*)(CObject *,int state))
    ObjectDistGetSettingHandle;
  I->Obj.fDescribeElement = NULL;
  I->CurDSet=0;
  I->Obj.Color=ColorGetIndex(G,"dash");
  return(I);
}

/*========================================================================*/
static void ObjectDistReset(PyMOLGlobals *G,ObjectDist *I)
{
  int a;
  for(a=0;a<I->NDSet;a++)
    if(I->DSet[a]) {
      if(I->DSet[a]->fFree)
          I->DSet[a]->fFree(I->DSet[a]);
      I->DSet[a]=NULL;
    }
  I->NDSet=0;
}

/*========================================================================*/
ObjectDist *ObjectDistNewFromSele(PyMOLGlobals *G,ObjectDist *oldObj,
                                  int sele1,int sele2,int mode,float cutoff,
                                  int labels,int reset,float *result,int state)
{
  int a,mn;
  float dist_sum=0.0,dist;
  int dist_cnt = 0;
  int n_state1,n_state2,state1,state2;
  ObjectDist *I;
  if(!oldObj)
    I=ObjectDistNew(G);
  else {
    I=oldObj;
    if(reset) 
      ObjectDistReset(G,I);
  }
  *result = 0.0;
  mn = 0;
  SelectorUpdateTable(G,state,-1);
  n_state1 = SelectorGetSeleNCSet(G,sele1);
  n_state2 = SelectorGetSeleNCSet(G,sele2);
  mn = n_state1;
  if(n_state2>mn) mn = n_state2;
  if(mn) {
    for(a=0;a<mn;a++) {
      if(state>=0) {
        if(state>=mn)
          break;
        a = state;
      }
        
      VLACheck(I->DSet,DistSet*,a);
      if(n_state1>1)
        state1=a;
      else
        state1=0;
      if(n_state2>1)
        state2=a;
      else
        state2=0;

      I->DSet[a] = SelectorGetDistSet(G,I->DSet[a],sele1,state1,sele2,
                                      state2,mode,cutoff,&dist);

      if(I->DSet[a]) {
        dist_sum+=dist;
        dist_cnt++;
        I->DSet[a]->Obj = I;
        I->NDSet=a+1;
      }

      if(state>=0)
        break;
    }
  } 
  ObjectDistUpdateExtents(I);
  ObjectDistInvalidateRep(I,cRepAll);

  if(dist_cnt)
    (*result) = dist_sum/dist_cnt;
  SceneChanged(G);
  return(I);
}

ObjectDist *ObjectDistNewFromAngleSele(PyMOLGlobals *G,ObjectDist *oldObj,
                                       int sele1, int sele2, int sele3, int mode,
                                       int labels,float *result, int reset,
                                       int state)
{
  int a,mn;
  float angle_sum=0.0;
  int angle_cnt=0;
  int n_state1,n_state2,n_state3, state1,state2,state3;
  ObjectDist *I;

  if(!oldObj) /* create object if new */
    I=ObjectDistNew(G);
  else { /* otherwise, use existing object */
    I=oldObj;
    if(reset) { /* if reseting, then clear out all existing coordinate sets */
      ObjectDistReset(G,I);
    }
  }
  *result = 0.0;

  /* count number of states in each selection */

  SelectorUpdateTable(G,state,-1);
  n_state1 = SelectorGetSeleNCSet(G,sele1);
  n_state2 = SelectorGetSeleNCSet(G,sele2);  
  n_state3 = SelectorGetSeleNCSet(G,sele3);

  /* figure out the total number of states */

  mn = n_state1;
  if(n_state2>mn) mn = n_state2;
  if(n_state3>mn) mn = n_state3;

  if(mn) {
    for(a=0;a<mn;a++) {
      if(state>=0) {
        if(state>mn)
          break;
        a = state;
      }
      /* treat selections with one state as static singletons */
      
      if(n_state1>1)
        state1=a;
      else
        state1=0;
      if(n_state2>1)
        state2=a;
      else
        state2=0;
      if(n_state3>1)
        state3=a;
      else
        state3=0;
      
      VLACheck(I->DSet,DistSet*,a);
      I->DSet[a] = SelectorGetAngleSet(G,I->DSet[a],sele1,state1,sele2,
                                       state2,sele3,state3,mode,&angle_sum,
                                       &angle_cnt);
      
      if(I->DSet[a]) {
        I->DSet[a]->Obj = I;
        if(I->NDSet<=a) I->NDSet = a+1;
      }
      if(state>=0)
        break;
    }  
  }
  /* else {
     VLAFreeP(I->DSet);
     OOFreeP(I);
     }
  */
  ObjectDistUpdateExtents(I);
  ObjectDistInvalidateRep(I,cRepAll);
  if(angle_cnt)
    (*result) = angle_sum/angle_cnt;

  SceneChanged(G);
  return(I);
}

ObjectDist *ObjectDistNewFromDihedralSele(PyMOLGlobals *G,ObjectDist *oldObj,
                                          int sele1, int sele2, int sele3, int sele4,
                                          int mode, int labels, float *result,
                                          int reset,int state)
{
  int a,mn;
  float angle_sum=0.0;
  int angle_cnt=0;
  int n_state1,n_state2,n_state3,n_state4, state1,state2,state3,state4;
  ObjectDist *I;

  if(!oldObj) /* create object if new */
    I=ObjectDistNew(G);
  else { /* otherwise, use existing object */
    I=oldObj;
    if(reset) { /* if reseting, then clear out all existing coordinate sets */
      ObjectDistReset(G,I);
    }
  }
  *result = 0.0;

  /* count number of states in each selection */

  SelectorUpdateTable(G,state,-1);

  n_state1 = SelectorGetSeleNCSet(G,sele1);
  n_state2 = SelectorGetSeleNCSet(G,sele2);  
  n_state3 = SelectorGetSeleNCSet(G,sele3);
  n_state4 = SelectorGetSeleNCSet(G,sele4);

  /* figure out the total number of states */

  mn = n_state1;
  if(n_state2>mn) mn = n_state2;
  if(n_state3>mn) mn = n_state3;
  if(n_state4>mn) mn = n_state4;

  if(mn) {
    for(a=0;a<mn;a++) {
      if(state>=0) {
        if(state>mn)
          break;
        a = state;
      }
      /* treat selections with one state as static singletons */
      
      if(n_state1>1)
        state1=a;
      else
        state1=0;
      if(n_state2>1)
        state2=a;
      else
        state2=0;
      if(n_state3>1)
        state3=a;
      else
        state3=0;
      if(n_state4>1)
        state4=a;
      else
        state4=0;
      
      VLACheck(I->DSet,DistSet*,a);
      I->DSet[a] = SelectorGetDihedralSet(G,I->DSet[a],sele1,state1,sele2,
                                          state2,sele3,state3,sele4,state4,
                                          mode,&angle_sum,
                                          &angle_cnt);
      
      if(I->DSet[a]) {
        I->DSet[a]->Obj = I;
        if(I->NDSet<=a) I->NDSet = a+1;
      }

      if(state>=0)
        break;
    }
  } 
  /* else {
     VLAFreeP(I->DSet);
     OOFreeP(I);
     }
  */
  ObjectDistUpdateExtents(I);
  ObjectDistInvalidateRep(I,cRepAll);

  if(angle_cnt)
    (*result) = angle_sum/angle_cnt;

  SceneChanged(G);
  return(I);
}

/*========================================================================*/
void ObjectDistFree(ObjectDist *I)
{
  int a;
  SceneObjectDel(I->Obj.G,(CObject*)I);
  for(a=0;a<I->NDSet;a++)
	 if(I->DSet[a]) {
      if(I->DSet[a]->fFree)
        I->DSet[a]->fFree(I->DSet[a]);
		I->DSet[a]=NULL;
	 }
  VLAFreeP(I->DSet);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
}

