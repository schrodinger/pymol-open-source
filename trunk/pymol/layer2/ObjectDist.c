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

void ObjectDistRender(ObjectDist *I,int frame,CRay *ray,Pickable **pick,int pass);
void ObjectDistFree(ObjectDist *I);
void ObjectDistUpdate(ObjectDist *I);
int ObjectDistGetNFrames(ObjectDist *I);
void ObjectDistUpdateExtents(ObjectDist *I);

static DistSet *ObjectDistGetDistSetFromM4XBond(ObjectMolecule *obj, M4XBondType *hb, int n_hb,
                                                 int state)
{
  int min_id,max_id,range,*lookup = NULL;
  int nv = 0;
  float *vv=NULL;
  DistSet *ds;
  ds = DistSetNew();
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
  
  FreeP(lookup);
  ds->NIndex = nv;
  ds->Coord = vv;
  return(ds);
  
}

ObjectDist *ObjectDistNewFromM4XBond(ObjectDist *oldObj,
                                      struct ObjectMolecule *objMol,
                                      struct M4XBondType *hbond,int n_hbond)
{
  int a;
  ObjectDist *I;
  int n_state;
  if(!oldObj)
    I=ObjectDistNew();
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
  for(a=0;a<n_state;a++)
    {
      VLACheck(I->DSet,DistSet*,a);
      
      I->DSet[a] = ObjectDistGetDistSetFromM4XBond(objMol,hbond,n_hbond,a);
      
      if(I->DSet[a]) {
        I->DSet[a]->Obj = I;
        I->NDSet=a+1;
      }
    } 
  ObjectDistUpdateExtents(I);
  
  SceneChanged();
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



static PyObject *ObjectDistDSetAsPyList(ObjectDist *I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NDSet);
  for(a=0;a<I->NDSet;a++) {
    if(I->DSet[a]) {
      PyList_SetItem(result,a,DistSetAsPyList(I->DSet[a]));
    } else {
      PyList_SetItem(result,a,Py_None);
      Py_INCREF(Py_None);
    }
  }
  return(PConvAutoNone(result));
}

static int ObjectDistDSetFromPyList(ObjectDist *I,PyObject *list)
{
  int ok=true;
  int a;
  if(ok) ok=PyList_Check(list);
  if(ok) {
    VLACheck(I->DSet,DistSet*,I->NDSet);
    for(a=0;a<I->NDSet;a++) {
      if(ok) ok = DistSetFromPyList(PyList_GetItem(list,a),&I->DSet[a]);
      if(ok) I->DSet[a]->Obj = I;
    }
  }
  return(ok);
}
/*========================================================================*/
PyObject *ObjectDistAsPyList(ObjectDist *I)
{
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
}

int ObjectDistNewFromPyList(PyObject *list,ObjectDist **result)
{
  int ok = true;
  ObjectDist *I=NULL;
  (*result) = NULL;
  
  if(ok) ok=PyList_Check(list);

  I=ObjectDistNew();
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(PyList_GetItem(list,0),&I->Obj);
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
  OrthoBusyPrime();
  for(a=0;a<I->NDSet;a++)
	 if(I->DSet[a]) {	
	   OrthoBusySlow(a,I->NDSet);
      /*	   printf(" ObjectDist: updating state %d of \"%s\".\n" , a+1, I->Obj.Name);*/
      if(I->DSet[a]->fUpdate)
        I->DSet[a]->fUpdate(I->DSet[a]);
	 }
}
/*========================================================================*/
void ObjectDistInvalidateRep(ObjectDist *I,int rep)
{
  int a;
  PRINTFD(FB_ObjectDist)
    " ObjectDistInvalidateRep: entered.\n"
    ENDFD;

  for(a=0;a<I->NDSet;a++) 
	 if(I->DSet[a]) {	 
      if(I->DSet[a]->fInvalidateRep)
        I->DSet[a]->fInvalidateRep(I->DSet[a],rep,cRepInvAll);
	 }
}
/*========================================================================*/
void ObjectDistRender(ObjectDist *I,int frame,CRay *ray,Pickable **pick,int pass)
{
  int a;
  if(!pass) {

    ObjectPrepareContext(&I->Obj,ray);
    if(frame<0) {
      for(a=0;a<I->NDSet;a++)
        if(I->DSet[a])
          if(I->DSet[a]->fRender)
            I->DSet[a]->fRender(I->DSet[a],ray,pick,pass);        
    } else if(frame<I->NDSet) {
      I->CurDSet=frame % I->NDSet;
      if(I->DSet[I->CurDSet]) {
        if(I->DSet[I->CurDSet]->fRender)
          I->DSet[I->CurDSet]->fRender(I->DSet[I->CurDSet],ray,pick,pass);
      }
    } else if(I->NDSet==1) { /* if only one coordinate set, assume static */
      if(I->DSet[0]->fRender)
        I->DSet[0]->fRender(I->DSet[0],ray,pick,pass);    
    }
  }
}

/*========================================================================*/
ObjectDist *ObjectDistNew(void)
{
  OOAlloc(ObjectDist);
  ObjectInit((CObject*)I);
  I->Obj.type=cObjectDist;
  I->DSet=VLAMalloc(10,sizeof(DistSet*),5,true); /* auto-zero */
  I->NDSet=0;
  I->Obj.fRender=(void (*)(struct CObject *, int, CRay *, Pickable **,int))ObjectDistRender;
  I->Obj.fFree= (void (*)(struct CObject *))ObjectDistFree;
  I->Obj.fUpdate= (void (*)(struct CObject *)) ObjectDistUpdate;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectDistGetNFrames;
  I->Obj.fDescribeElement = NULL;
  I->CurDSet=0;
  I->Obj.Color=ColorGetIndex("dash");
  return(I);
}

/*========================================================================*/
/*========================================================================*/
ObjectDist *ObjectDistNewFromSele(ObjectDist *oldObj,int sele1,int sele2,int mode,float cutoff,
                                  int labels,float *result)
{
  int a,mn;
  float dist_sum=0.0,dist;
  int dist_cnt = 0;
  int n_state1,n_state2,state1,state2;
  ObjectDist *I;

  if(!oldObj)
    I=ObjectDistNew();
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
  *result = 0.0;
  mn = 0;
  SelectorUpdateTable();
  n_state1 = SelectorGetSeleNCSet(sele1);
  n_state2 = SelectorGetSeleNCSet(sele2);
  mn = n_state1;
  if(n_state2>mn) mn = n_state2;
  if(mn) {
    for(a=0;a<mn;a++)
      {
        VLACheck(I->DSet,DistSet*,a);
        if(n_state1>1)
          state1=a;
        else
          state1=0;
        if(n_state2>1)
          state2=a;
        else
          state2=0;
        I->DSet[a] = SelectorGetDistSet(sele1,state1,sele2,
                                        state2,mode,cutoff,&dist);
        if(I->DSet[a]) {
          dist_sum+=dist;
          dist_cnt++;
          I->DSet[a]->Obj = I;
          I->NDSet=a+1;
        }
      }  
  } 
  /* else {
     VLAFreeP(I->DSet);
     OOFreeP(I);
     }
  */
  ObjectDistUpdateExtents(I);

  if(dist_cnt)
    (*result) = dist_sum/dist_cnt;
  SceneChanged();
  return(I);
}

/*========================================================================*/
void ObjectDistFree(ObjectDist *I)
{
  int a;
  SceneObjectDel((CObject*)I);
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

