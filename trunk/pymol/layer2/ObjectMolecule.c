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

#include"Base.h"
#include"Debug.h"
#include"Parse.h"
#include"OOMac.h"
#include"Vector.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Map.h"
#include"Selector.h"
#include"ObjectMolecule.h"
#include"Ortho.h"
#include"Util.h"
#include"Vector.h"
#include"Selector.h"
#include"Matrix.h"
#include"Scene.h"
#include"P.h"
#include"PConv.h"
#include"Executive.h"
#include"Setting.h"
#include"Sphere.h"

#define wcopy ParseWordCopy
#define nextline ParseNextLine
#define ncopy ParseNCopy
#define nskip ParseNSkip

void ObjectMoleculeRender(ObjectMolecule *I,int frame,CRay *ray,Pickable **pick);
void ObjectMoleculeCylinders(ObjectMolecule *I);
CoordSet *ObjectMoleculeMMDStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);
CoordSet *ObjectMoleculePDBStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);
CoordSet *ObjectMoleculeMOLStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);
void ObjectMoleculeAppendAtoms(ObjectMolecule *I,AtomInfoType *atInfo,CoordSet *cset);

void ObjectMoleculeFree(ObjectMolecule *I);
void ObjectMoleculeUpdate(ObjectMolecule *I);
int ObjectMoleculeGetNFrames(ObjectMolecule *I);

void ObjectMoleculeDescribeElement(ObjectMolecule *I,int index);

void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op);

int ObjectMoleculeConnect(ObjectMolecule *I,int **bond,AtomInfoType *ai,CoordSet *cs,float cutoff,int searchFlag);
void ObjectMoleculeTransformTTTf(ObjectMolecule *I,float *ttt,int state);
static int BondInOrder(int *a,int b1,int b2);
static int BondCompare(int *a,int *b);

CoordSet *ObjectMoleculeChemPyModel2CoordSet(PyObject *model,AtomInfoType **atInfoPtr);

int ObjectMoleculeGetAtomGeometry(ObjectMolecule *I,int state,int at);
void ObjectMoleculeBracketResidue(ObjectMolecule *I,AtomInfoType *ai,int *st,int *nd);

int ObjectMoleculeFindOpenValenceVector(ObjectMolecule *I,int state,int index,float *v);
void ObjectMoleculeAddSeleHydrogens(ObjectMolecule *I,int sele);

/*========================================================================*/
int ObjectMoleculeAreAtomsBonded(ObjectMolecule *I,int i0,int i1)
{
  int result=false;
  int a,*b;
  b=I->Bond;
  for (a=0;a<I->NBond;a++) {
    if(i0==*b) {
      if(i1==(*(b+1))) {
        result=true;
        break;
      }
    }
    if(i1==*b) {
      if(i0==*(b+1)) {
        result=true;
        break;
      }
    }
    b+=3;
  }
  return(result);
}
/*========================================================================*/
void ObjectMoleculeRenameAtoms(ObjectMolecule *I,int force)
{
  AtomInfoType *ai;
  int a;
  if(force) {
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++)
      (ai++)->name[0]=0;
  }
  AtomInfoUniquefyNames(NULL,0,I->AtomInfo,I->NAtom);  
}
/*========================================================================*/
void ObjectMoleculeAddSeleHydrogens(ObjectMolecule *I,int sele)
{
  int a,b;
  int n,nn;
  CoordSet *cs;
  CoordSet *tcs;
  int seleFlag=false;
  AtomInfoType *ai,*nai,fakeH;
  int repeatFlag=false;
  int nH;
  int *index;
  float v[3],v0[3];
  float d;

  UtilZeroMem(&fakeH,sizeof(AtomInfoType));
  fakeH.protons=1;
  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(SelectorIsMember(ai->selEntry,sele)) {
      seleFlag=true;
      break;
    }
    ai++;
  }
  if(seleFlag) {
    if(!ObjectMoleculeVerifyChemistry(I)) {
      ErrMessage(" AddHydrogens","missing chemical geometry information.");
    } else if(I->DiscreteFlag) {
      ErrMessage(" AddHydrogens","can't modify a discrete object.");
    } else {

      repeatFlag=true;
      while(repeatFlag) {
        repeatFlag=false;
        nH = 0;
        ObjectMoleculeUpdateNeighbors(I);
        nai = (AtomInfoType*)VLAMalloc(1000,sizeof(AtomInfoType),1,true);        
        ai=I->AtomInfo;
        for(a=0;a<I->NAtom;a++) {
          if(SelectorIsMember(ai->selEntry,sele)) {
            n = I->Neighbor[a];
            nn = I->Neighbor[n++];
            if(nn<ai->valence) {
              VLACheck(nai,AtomInfoType,nH);
              UtilNCopy((nai+nH)->elem,"H",2);
              (nai+nH)->geom=cAtomInfoSingle;
              (nai+nH)->valence=1;
              (nai+nH)->id = a; /* borrowing this field temporarily */
              ObjectMoleculePrepareAtom(I,a,nai+nH);
              nH++;
            }
          }
          ai++;
        }

        if(nH) {

          repeatFlag=true;
          cs = CoordSetNew();
          cs->Coord = VLAlloc(float,nH*3);
          cs->NIndex=nH;

          index = Alloc(int,nH);
          for(a=0;a<nH;a++) {
            index[a] = (nai+a)->id;
          }
          
          if(cs->fEnumIndices) cs->fEnumIndices(cs);

          cs->TmpLinkBond = VLAlloc(int,nH*3);
          for(a=0;a<nH;a++) {
            cs->TmpLinkBond[a*3] = (nai+a)->id;
            cs->TmpLinkBond[a*3+1] = a;
            cs->TmpLinkBond[a*3+2] = 1;
          }
          cs->NTmpLinkBond = nH;

          AtomInfoUniquefyNames(I->AtomInfo,I->NAtom,nai,nH);

          ObjectMoleculeMerge(I,nai,cs,false); /* will free nai and cs->TmpLinkBond  */
          ObjectMoleculeExtendIndices(I);
          ObjectMoleculeUpdateNeighbors(I);

          for(b=0;b<I->NCSet;b++) { /* add coordinate into the coordinate set */
            tcs = I->CSet[b];
            if(tcs) {
              for(a=0;a<nH;a++) {
                ObjectMoleculeGetAtomVertex(I,n,index[a],v0);
                ObjectMoleculeFindOpenValenceVector(I,b,index[a],v);
                d = AtomInfoGetBondLength(I->AtomInfo+index[a],&fakeH);
                scale3f(v,d,v);
                add3f(v0,v,cs->Coord+3*a);
              }
              CoordSetMerge(tcs,cs);  
            }
          }
          FreeP(index);
          if(cs->fFree)
            cs->fFree(cs);
          ObjectMoleculeSort(I);
        } else
          VLAFreeP(nai);
      }
    }
  }
}
/*========================================================================*/
void ObjectMoleculeFuse(ObjectMolecule *I,int index0,ObjectMolecule *src,int index1,int mode)
{
  int a,b;
  AtomInfoType *ai0,*ai1,*nai;
  int n,nn;
  int at0=-1;
  int at1=-1;
  int a0,a1,o1;
  int hydr1=-1;
  int anch1=-1;
  int ca0,ch0;
  int *b0,*b1;
  float *backup = NULL;
  float d,*f0,*f1;
  float va0[3],vh0[3],va1[3],vh1[3];
  float x0[3],y0[3],z0[3];
  float x1[3],y1[3],z1[3];
  float x[3],y[3],z[3];
  float t[3],t2[3];
  CoordSet *cs,*scs=NULL;
  int state1 = 0;
  CoordSet *tcs;

  ObjectMoleculeUpdateNeighbors(I);
  ObjectMoleculeUpdateNeighbors(src);

  /* make sure each link point has only one neighbor */

  ai0=I->AtomInfo;
  ai1=src->AtomInfo;
  switch(mode) {
  case 0: /* fusing by replacing hydrogens */
    
    n = I->Neighbor[index0];
    nn = I->Neighbor[n++];
    if(nn==1)
      at0 = I->Neighbor[n];
    
    n = src->Neighbor[index1];
    nn = src->Neighbor[n++];
    if(nn==1)
      at1 = src->Neighbor[n];
    
    if(src->NCSet) {
      scs = src->CSet[state1];
      anch1 = scs->AtmToIdx[at1];
      hydr1 = scs->AtmToIdx[index1];
    }
    break;
  case 1: /* fuse merely by drawing a bond */
    at0 = index0;
    at1 = index1;
    
    if(src->NCSet) {
      scs = src->CSet[state1];
      anch1 = scs->AtmToIdx[at1];
    }

    break;
  }
  
  if((at0>=0)&&(at1>=0)&&scs&&(anch1>=0)) { /* have anchors and source coordinate set */

    nai=(AtomInfoType*)VLAMalloc(src->NAtom,sizeof(AtomInfoType),1,true);
    
    /* copy atoms and atom info into a 1:1 direct mapping */

    cs = CoordSetNew();
    cs->Coord = VLAlloc(float,scs->NIndex*3);
    cs->NIndex = scs->NIndex;
    for(a=0;a<scs->NIndex;a++) {
      copy3f(scs->Coord+a*3,cs->Coord+a*3);
      a1 = scs->IdxToAtm[a];
      *(nai+a) = *(ai1+a1);
      (nai+a)->selEntry=0; /* avoid duplicating selection references -> leads to hangs ! */
    }

    /* copy internal bond information*/

    cs->TmpBond = VLAlloc(int,src->NBond*3);
    b1 = src->Bond;
    b0 = cs->TmpBond;
    cs->NTmpBond=0;
    for(a=0;a<src->NBond;a++) {
      a0 = scs->AtmToIdx[*(b1++)];
      a1 = scs->AtmToIdx[*(b1++)];
      o1 = *(b1++);
      if((a0>=0)&&(a1>=0)) {
        *(b0++) = a0;
        *(b0++) = a1;
        *(b0++) = o1;
        cs->NTmpBond++;
      }
    }

    backup = Alloc(float,cs->NIndex*3); /* make untransformed copy of coordinate set */
    for(a=0;a<cs->NIndex;a++) {
      copy3f(cs->Coord+a*3,backup+a*3);
    }
    
    switch(mode) {
    case 0:
      nai[hydr1].deleteFlag=true;
      I->AtomInfo[index0].deleteFlag=true;
      copy3f(backup+3*anch1,va1);
      copy3f(backup+3*hydr1,vh1);
      subtract3f(va1,vh1,x1); /* note reverse dir from above */
      get_system1f3f(x1,y1,z1);
      break;
    case 1:
      copy3f(backup+3*anch1,va1);
      ObjectMoleculeFindOpenValenceVector(src,state1,at1,x1);
      scale3f(x1,-1.0,x1);
      get_system1f3f(x1,y1,z1);      
      break;
    }

    /* set up the linking bond */

    cs->TmpLinkBond = VLAlloc(int,3);
    cs->NTmpLinkBond = 1;
    cs->TmpLinkBond[0]=at0;
    cs->TmpLinkBond[1]=anch1;
    cs->TmpLinkBond[2]=1;

    if(cs->fEnumIndices) cs->fEnumIndices(cs);

    d = AtomInfoGetBondLength(ai0+at0,ai1+at1);

    AtomInfoUniquefyNames(I->AtomInfo,I->NAtom,nai,cs->NIndex);
    ObjectMoleculeMerge(I,nai,cs,false); /* will free nai, cs->TmpBond and cs->TmpLinkBond  */
    ObjectMoleculeExtendIndices(I);
    ObjectMoleculeUpdateNeighbors(I);
    for(a=0;a<I->NCSet;a++) { /* add coordinate into the coordinate set */
      tcs = I->CSet[a];
      if(tcs) {
        switch(mode) {
        case 0:
          ca0 = tcs->AtmToIdx[at0]; /* anchor */
          ch0 = tcs->AtmToIdx[index0]; /* hydrogen */

          if((ca0>=0)&&(ch0>=0)) {
            copy3f(tcs->Coord+3*ca0,va0);
            copy3f(tcs->Coord+3*ch0,vh0);
            subtract3f(vh0,va0,x0);
            get_system1f3f(x0,y0,z0);

          }
          break;
        case 1:
          ca0 = tcs->AtmToIdx[at0]; /* anchor */

          if((ca0>=0)&&(ch0>=0)) {
            ObjectMoleculeFindOpenValenceVector(I,a,at0,x0);
            copy3f(tcs->Coord+3*ca0,va0);
            get_system1f3f(x0,y0,z0);
            
          }
          break;
        }
        scale3f(x0,d,t2);
        add3f(va0,t2,t2);
        
        f0=backup;
        f1=cs->Coord;
        for(b=0;b<cs->NIndex;b++) { /* brute force transformation */
          subtract3f(f0,va1,t);
          scale3f(x0,dot_product3f(t,x1),x);
          scale3f(y0,dot_product3f(t,y1),y);
          scale3f(z0,dot_product3f(t,z1),z);
          add3f(x,y,y);
          add3f(y,z,f1);
          add3f(t2,f1,f1);
          f0+=3;
          f1+=3;
        }
        CoordSetMerge(tcs,cs); 
      }
    }
    switch(mode) {
    case 0:
      ObjectMoleculePurge(I);
      break;
    }
    ObjectMoleculeSort(I);
  }
  if(cs->fFree)
    cs->fFree(cs);
  FreeP(backup);
}
/*========================================================================*/
int ObjectMoleculeVerifyChemistry(ObjectMolecule *I)
{
  int result=false;
  AtomInfoType *ai;
  int a;
  int flag;
  ai=I->AtomInfo;
  flag=true;
  for(a=0;a<I->NAtom;a++) {
    if(!ai->chemFlag) {
      flag=false;
    }
    ai++;
  }
  if(!flag) {
    if(I->CSet[0]) { /* right now this stuff is locked to state 0 */
      ObjectMoleculeInferChemFromBonds(I,0);
      ObjectMoleculeInferChemFromNeighGeom(I,0);
      /*      ObjectMoleculeInferChemForProtein(I,0);*/
    }
    flag=true;
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(!ai->chemFlag) {
        flag=false;
        break;
      }
      ai++;
    }
  }
  if(flag)
    result=true;
  return(result);
}
/*========================================================================*/
void ObjectMoleculeAttach(ObjectMolecule *I,int index,AtomInfoType *nai)
{
  int a;
  AtomInfoType *ai;
  int n,nn;
  float v[3],v0[3],d;
  CoordSet *cs;

  ObjectMoleculeUpdateNeighbors(I);
  ai=I->AtomInfo+index;
  n = I->Neighbor[index];
  nn = I->Neighbor[n++];
  
  cs = CoordSetNew();
  cs->Coord = VLAlloc(float,3);
  cs->NIndex=1;
  cs->TmpLinkBond = VLAlloc(int,3);
  cs->NTmpLinkBond = 1;
  cs->TmpLinkBond[0]=index;
  cs->TmpLinkBond[1]=0;
  cs->TmpLinkBond[2]=1;
  if(cs->fEnumIndices) cs->fEnumIndices(cs);
  ObjectMoleculePrepareAtom(I,index,nai);
  d = AtomInfoGetBondLength(ai,nai);
  ObjectMoleculeMerge(I,nai,cs,false); /* will free nai and cs->TmpLinkBond  */
  ObjectMoleculeExtendIndices(I);
  ObjectMoleculeUpdateNeighbors(I);
  for(a=0;a<I->NCSet;a++) { /* add atom to each coordinate set */
    if(I->CSet[a]) {
      ObjectMoleculeGetAtomVertex(I,a,index,v0);
      ObjectMoleculeFindOpenValenceVector(I,a,index,v);
      scale3f(v,d,v);
      add3f(v0,v,cs->Coord);
      CoordSetMerge(I->CSet[a],cs); 
    }
  }
  ObjectMoleculeSort(I);
  if(cs->fFree)
    cs->fFree(cs);
  
}
/*========================================================================*/
int ObjectMoleculeFillOpenValences(ObjectMolecule *I,int index)
{
  int a;
  AtomInfoType *ai,*nai;
  int n,nn;
  int result=0;
  int flag = true;
  float v[3],v0[3],d;
  CoordSet *cs;

  if((index>=0)&&(index<=I->NAtom)) {  
    while(1) {
      ObjectMoleculeUpdateNeighbors(I);
      ai=I->AtomInfo+index;
      n = I->Neighbor[index];
      nn = I->Neighbor[n++];
      
      if((nn>=ai->valence)||(!flag))
        break;
      flag=false;

      cs = CoordSetNew();
      cs->Coord = VLAlloc(float,3);
      cs->NIndex=1;
      cs->TmpLinkBond = VLAlloc(int,3);
      cs->NTmpLinkBond = 1;
      cs->TmpLinkBond[0]=index;
      cs->TmpLinkBond[1]=0;
      cs->TmpLinkBond[2]=1;
      if(cs->fEnumIndices) cs->fEnumIndices(cs);
      nai = (AtomInfoType*)VLAMalloc(1,sizeof(AtomInfoType),1,true);
      UtilNCopy(nai->elem,"H",2);
      nai->geom=cAtomInfoSingle;
      nai->valence=1;
      ObjectMoleculePrepareAtom(I,index,nai);
      d = AtomInfoGetBondLength(ai,nai);
      ObjectMoleculeMerge(I,nai,cs,false); /* will free nai and cs->TmpLinkBond  */
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeUpdateNeighbors(I);
      for(a=0;a<I->NCSet;a++) { /* add atom to each coordinate set */
        if(I->CSet[a]) {
          ObjectMoleculeGetAtomVertex(I,a,index,v0);
          ObjectMoleculeFindOpenValenceVector(I,a,index,v);
          scale3f(v,d,v);
          add3f(v0,v,cs->Coord);
          CoordSetMerge(I->CSet[a],cs); 
        }
      }
      if(cs->fFree)
        cs->fFree(cs);
      result++;
      flag=true;
    }
  }
  return(result);
}
/*========================================================================*/
int ObjectMoleculeFindOpenValenceVector(ObjectMolecule *I,int state,int index,float *v)
{
  #define MaxOcc 100
  CoordSet *cs;
  int nOcc = 0;
  float occ[MaxOcc*3];
  int n;
  int a1;
  float v0[3],v1[3],n0[3],t[3];
  int result = false;
  AtomInfoType *ai,*ai1;
  float y[3],z[3];

  /* default is +X */
  v[0]=1.0;
  v[1]=0.0;
  v[2]=0.0;
  
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    if((index>=0)&&(index<=I->NAtom)) {
      ai=I->AtomInfo+index;
      if(ObjectMoleculeGetAtomVertex(I,state,index,v0)) {              
        ObjectMoleculeUpdateNeighbors(I);
        n = I->Neighbor[index];
        n++; /* skip count */
        while(1) { /* look for an attached non-hydrogen as a base */
          a1 = I->Neighbor[n];
          n+=2; 
          if(a1<0) break;
          ai1=I->AtomInfo+a1;
          if(ObjectMoleculeGetAtomVertex(I,state,a1,v1)) {        
            subtract3f(v1,v0,n0);
            normalize3f(n0);
            copy3f(n0,occ+3*nOcc);
            nOcc++; 
            if(nOcc==MaxOcc) /* safety valve */
              break;
          }
        }
        if((!nOcc)||(nOcc>4)||(ai->geom==cAtomInfoNone)) {
          get_random3f(v);
        } else {
          switch(nOcc) {
          case 1:
            switch(ai->geom) {
            case cAtomInfoTetrahedral:
              get_system1f3f(occ,y,z);
              scale3f(occ,-0.334,v);
              scale3f(z,  0.943,t);
              add3f(t,v,v);              
              break;
            case cAtomInfoPlaner:
              get_system1f3f(occ,y,z);
              scale3f(occ,-0.500,v);
              scale3f(z,      0.866,t);
              add3f(t,v,v);
              break;
            case cAtomInfoLinear:
              scale3f(occ,-1.0,v);
              break;
            default:
              get_random3f(v); /* hypervalent */
              break;
            }
            break;
          case 2:
            switch(ai->geom) {
            case cAtomInfoTetrahedral:
              add3f(occ,occ+3,t);
              get_system2f3f(t,occ,z);
              scale3f(t,-1.0,v);
              scale3f(z,1.41,t);
              add3f(t,v,v);              
              break;
            case cAtomInfoPlaner:
            add3f(occ,occ+3,t);
              scale3f(t,-1.0,v);
              break;
            default:
              get_random3f(v); /* hypervalent */
              break;
            }
            break;
          case 3:
            switch(ai->geom) {
            case cAtomInfoTetrahedral:
              add3f(occ,occ+3,t);
              add3f(occ+6,t,t);
              scale3f(t,-1.0,v);
              break;
            default:
              get_random3f(v); /* hypervalent */
              break;
              
            }
            break;
          case 4:
            get_random3f(v); /* hypervalent */
            break;
          }
        }
      }
    }
  }
  normalize3f(v);
  return(result); /* return false if we're using the default */
#undef MaxOcc

}
/*========================================================================*/
void ObjectMoleculeCreateSpheroid(ObjectMolecule *I)
{
  CoordSet *cs;
  float *spheroid = NULL;
  int a,b,c,a0;
  SphereRec *sp;
  float *spl;
  float *v,*v0,*s,*f,ang,min_dist,*max_sq;
  int *i;
  float *center = NULL;
  float d0[3],n0[3],d1[3],d2[3];
  float p0[3],p1[3],p2[3];
  int t0,t1,t2,bt0,bt1,bt2;
  float dp,l,*fsum = NULL;
  float *norm = NULL;
  float spheroid_smooth;
  float spheroid_fill;
  float spheroid_ratio=0.1; /* minimum ratio of width over length */
  float spheroid_minimum = 0.02; /* minimum size - to insure valid normals */
  int row,*count=NULL,base;
  int nRow;
  sp=Sphere1;
  
  nRow = I->NAtom*sp->nDot;
  spheroid=Alloc(float,nRow);
  center=Alloc(float,I->NAtom*3);
  count=Alloc(int,I->NAtom);
  fsum=Alloc(float,nRow);
  max_sq = Alloc(float,I->NAtom);

  spl=spheroid;

  spheroid_smooth=SettingGet(cSetting_spheroid_smooth);
  spheroid_fill=SettingGet(cSetting_spheroid_fill);
  /* first compute average coordinate */

  v=center;
  i = count;
  for(a=0;a<I->NAtom;a++) {
    *(v++)=0.0;
    *(v++)=0.0;
    *(v++)=0.0;
    *(i++)=0;
  }

  for(b=0;b<I->NCSet;b++) {
    cs=I->CSet[b];
    if(cs) {
      v = cs->Coord;
      for(a=0;a<cs->NIndex;a++) {
        a0=cs->IdxToAtm[a];
        v0 = center+3*a0;
        add3f(v,v0,v0);
        (*(count+a0))++;
        v+=3;
      }
    }
  }

  i=count;
  v=center;
  for(a=0;a<I->NAtom;a++) 
    if(*i) {
      (*(v++))/=(*i);
      (*(v++))/=(*i);
      (*(v++))/=(*i++);
    } else {
      v+=3;
      i++;
    }

/* now go through and compute radial distances */

  f = fsum;
  s = spheroid;
  for(a=0;a<nRow;a++) {
    *(f++)=0.0;
    *(s++)=0.0; 
  }

  v = max_sq;
  for(a=0;a<I->NAtom;a++)
    *(v++)=0.0;

  for(b=0;b<I->NCSet;b++) {
    cs=I->CSet[b];
    if(cs) {
      v = cs->Coord;
      for(a=0;a<cs->NIndex;a++) {
        a0=cs->IdxToAtm[a];
        base = (a0*sp->nDot);
        v0 = center+(3*a0);
        subtract3f(v,v0,d0); /* subtract from average */
        l = lengthsq3f(d0);
        if(l>max_sq[a0])
          max_sq[a0]=l;
        if(l>0.0) {
          scale3f(d0,1.0/sqrt(l),n0);
          for(c=0;c<sp->nDot;c++) { /* average over spokes */
            dp=dot_product3f(sp->dot[c].v,n0);
            row = base + c;
            if(dp>=0.0) {
              ang = acos(dp);
              ang=(ang/spheroid_smooth)*(cPI/2.0); 
              if(ang>spheroid_fill)
                ang=spheroid_fill;
              /* take envelop to zero over that angle */
              if(ang<=(cPI/2.0)) {
                dp = cos(ang);
                fsum[row] += dp*dp;
                spheroid[row] += l*dp*dp*dp;
              }
            }
          }
        }
        v+=3;
      }
    }
  }

  f=fsum;
  s=spheroid;
  for(a=0;a<I->NAtom;a++) {
    min_dist = spheroid_ratio*sqrt(max_sq[a]);
    if(min_dist<spheroid_minimum)
      min_dist=spheroid_minimum;
    for(b=0;b<sp->nDot;b++) {
      if(*f>R_SMALL4) {
        (*s)=sqrt((*s)/(*(f++))); /* we put the "rm" in "rms" */
      } else {
        f++;
      }
      if(*s<min_dist)
        *s=min_dist;
      s++;
    }
  }

  /* set frame 0 coordinates to the average */

  for(b=0;b<I->NCSet;b++) {
    cs=I->CSet[b];
    if(cs) {
      v = cs->Coord;
      for(a=0;a<cs->NIndex;a++) {
        a0=cs->IdxToAtm[a];
        v0 = center+3*a0;
        copy3f(v0,v);
        v+=3;
      }
    }
  }

  /* now compute surface normals */

  norm = Alloc(float,nRow*3);
  for(a=0;a<I->NAtom;a++) {
    base = a*sp->nDot;
    for(b=0;b<sp->NTri;b++) {
      t0 = sp->Tri[b*3  ];
      t1 = sp->Tri[b*3+1];
      t2 = sp->Tri[b*3+2];
      bt0 = base + t0;
      bt1 = base + t1;
      bt2 = base + t2;
      copy3f(sp->dot[t0].v,p0);
      copy3f(sp->dot[t1].v,p1);
      copy3f(sp->dot[t2].v,p2);
      /*      scale3f(sp->dot[t0].v,spheroid[bt0],p0);
      scale3f(sp->dot[t1].v,spheroid[bt1],p1);
      scale3f(sp->dot[t2].v,spheroid[bt2],p2);*/
      subtract3f(p1,p0,d1);
      subtract3f(p2,p0,d2);
      cross_product3f(d1,d2,n0);
      normalize3f(n0);
      v = norm+bt0*3;
      add3f(n0,v,v);
      v = norm+bt1*3;
      add3f(n0,v,v);
      v = norm+bt2*3;
      add3f(n0,v,v);
    }
  }

  f=norm;
  for(a=0;a<I->NAtom;a++) {
    base = a*sp->nDot;
    for(b=0;b<sp->nDot;b++) {
      normalize3f(f);
      f+=3;
    }
  }
  
  if(I->CSet[0]) {
    I->CSet[0]->Spheroid=spheroid;
    I->CSet[0]->SpheroidNormal=norm;
    I->CSet[0]->NSpheroid=nRow;
  } else {
    FreeP(spheroid);
    FreeP(norm);
  }

  FreeP(center);
  FreeP(count);
  FreeP(fsum);
  FreeP(max_sq);

  for(b=1;b<I->NCSet;b++) { 
    cs=I->CSet[b];
    if(cs) {
      if(cs->fFree)
        cs->fFree(cs);
    }
  }
  I->NCSet=1;

  ObjectMoleculeInvalidate(I,cRepSphere,cRepInvProp);
}
/*========================================================================*/
void ObjectMoleculeReplaceAtom(ObjectMolecule *I,int index,AtomInfoType *ai)
{
  if((index>=0)&&(index<=I->NAtom)) {
    memcpy(I->AtomInfo+index,ai,sizeof(AtomInfoType));
    ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms);
    /* could we put in a refinement step here? */
  }
}
/*========================================================================*/
void ObjectMoleculePrepareAtom(ObjectMolecule *I,int index,AtomInfoType *ai)
{
  /* match existing properties of the old atom */
  int a;
  AtomInfoType *ai0;

  if((index>=0)&&(index<=I->NAtom)) {
    ai0=I->AtomInfo + index;
    ai->resv=ai0->resv;
    ai->hetatm=ai0->hetatm;
    ai->flags=ai0->flags;
    strcpy(ai->chain,ai0->chain);
    strcpy(ai->alt,ai0->alt);
    strcpy(ai->resi,ai0->resi);
    strcpy(ai->segi,ai0->segi);
    strcpy(ai->resn,ai0->resn);    
    if((ai->elem[0]==ai0->elem[0])&&(ai->elem[1]==ai0->elem[1]))
      ai->color=ai0->color;
    else
      ai->color=AtomInfoGetColor(ai);
    for(a=0;a<cRepCnt;a++)
      ai->visRep[a]=ai0->visRep[a];
    AtomInfoUniquefyNames(I->AtomInfo,I->NAtom,ai,1);
    AtomInfoAssignParameters(ai);
  }
}
/*========================================================================*/
void ObjectMoleculePreposReplAtom(ObjectMolecule *I,int index,
                                   AtomInfoType *ai)
{
  int n;
  int a1;
  AtomInfoType *ai1;
  float v0[3],v1[3],v[3];
  float d0[3],d,n0[3];
  int a;
  ObjectMoleculeUpdateNeighbors(I);
  for(a=0;a<I->NCSet;a++) {
    if(I->CSet[a]) {
      if(ObjectMoleculeGetAtomVertex(I,a,index,v0)) {
        copy3f(v0,v); /* default is direct superposition */
        n = I->Neighbor[index];
        n++; /* skip count */
        while(1) { /* look for an attached non-hydrogen as a base */
          a1 = I->Neighbor[n];
          n+=2;
          if(a1<0) break;
          ai1=I->AtomInfo+a1;
          if(ai1->protons!=1) 
            if(ObjectMoleculeGetAtomVertex(I,a,a1,v1)) {        
              d = AtomInfoGetBondLength(ai,ai1);
              subtract3f(v0,v1,n0);
              normalize3f(n0);
              scale3f(n0,d,d0);
              add3f(d0,v1,v);
              ObjectMoleculeSetAtomVertex(I,a,index,v);
              break;
            }
        }
        
      }
    }
  }
}
/*========================================================================*/
void ObjectMoleculeSaveUndo(ObjectMolecule *I,int state)
{
  CoordSet *cs;

  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter]=-1;
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    I->UndoCoord[I->UndoIter] = Alloc(float,cs->NIndex*3);
    memcpy(I->UndoCoord[I->UndoIter],cs->Coord,sizeof(float)*cs->NIndex*3);
    I->UndoState[I->UndoIter]=state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }
  I->UndoIter=cUndoMask&(I->UndoIter+1);
}
/*========================================================================*/
void ObjectMoleculeUndo(ObjectMolecule *I,int dir)
{
  CoordSet *cs;
  int state;

  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter]=-1;
  state=SceneGetState();
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    I->UndoCoord[I->UndoIter] = Alloc(float,cs->NIndex*3);
    memcpy(I->UndoCoord[I->UndoIter],cs->Coord,sizeof(float)*cs->NIndex*3);
    I->UndoState[I->UndoIter]=state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }

  I->UndoIter=cUndoMask&(I->UndoIter+dir);
  if(!I->UndoCoord[I->UndoIter])
    I->UndoIter=cUndoMask&(I->UndoIter-dir);

  if(I->UndoState[I->UndoIter]>=0) {
    state=I->UndoState[I->UndoIter];
    if(state<0) state=0;
    
    if(I->NCSet==1) state=0;
    state = state % I->NCSet;
    cs = I->CSet[state];
    if(cs) {
      if(cs->NIndex==I->UndoNIndex[I->UndoIter]) {
        memcpy(cs->Coord,I->UndoCoord[I->UndoIter],sizeof(float)*cs->NIndex*3);
        I->UndoState[I->UndoIter]=-1;
        FreeP(I->UndoCoord[I->UndoIter]);
        if(cs->fInvalidateRep)
          cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
        SceneChanged();
      }
    }
  }
}
/*========================================================================*/
int ObjectMoleculeAddBond(ObjectMolecule *I,int sele0,int sele1,int order)
{
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int s1,s2;
  int c = 0;
  int *bnd;

  ai1=I->AtomInfo;
  for(a1=0;a1<I->NAtom;a1++) {
    s1=ai1->selEntry;
    if(SelectorIsMember(s1,sele0)) {
      ai2=I->AtomInfo;
      for(a2=0;a2<I->NAtom;a2++) {
        s2=ai2->selEntry;
        if(SelectorIsMember(s2,sele1)) {
          {
            VLACheck(I->Bond,int,I->NBond*3+2);
            bnd = I->Bond+(3*I->NBond);
            bnd[0]=a1;
            bnd[1]=a2;                      
            bnd[2]=order;
            I->NBond++;
            c++;
            I->AtomInfo[a1].chemFlag=false;
            I->AtomInfo[a2].chemFlag=false;
          }
        }
        ai2++;
      }
    }
    ai1++;
  }
  if(c) {
    ObjectMoleculeInvalidate(I,cRepLine,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepCyl,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepNonbonded,cRepInvBonds);
  }
  return(c);    
}
/*========================================================================*/
int ObjectMoleculeAdjustBonds(ObjectMolecule *I,int sele0,int sele1,int mode,int order)
{
  int a0,a1;
  int offset=0;
  int *b0;
  int both;
  int s;
  int a;

  offset=0;
  b0=I->Bond;
  for(a=0;a<I->NBond;a++) {
    a0=b0[0];
    a1=b0[1];
    
    both=0;
    s=I->AtomInfo[a0].selEntry;
    if(SelectorIsMember(s,sele0))
      both++;
    s=I->AtomInfo[a1].selEntry;
    if(SelectorIsMember(s,sele1))
      both++;
    if(both<2) { /* reverse combo */
      both=0;
      s=I->AtomInfo[a1].selEntry;
      if(SelectorIsMember(s,sele0))
        both++;
      s=I->AtomInfo[a0].selEntry;
      if(SelectorIsMember(s,sele1))
        both++;
    }

    if(both==2) {
      switch(mode) {
      case 0: /* cycle */
        b0[2]++;
        if(b0[2]>3)
          b0[2]=1;
        I->AtomInfo[a0].chemFlag=false;
        I->AtomInfo[a1].chemFlag=false;
        break;
      case 1: /* set */
        b0[2]=order;
        I->AtomInfo[a0].chemFlag=false;
        I->AtomInfo[a1].chemFlag=false;
        break;
      }
      ObjectMoleculeInvalidate(I,cRepLine,cRepInvBonds);
      ObjectMoleculeInvalidate(I,cRepCyl,cRepInvBonds);
      ObjectMoleculeInvalidate(I,cRepNonbonded,cRepInvBonds);
    }
    b0+=3;
  }
  return(-offset);
}
/*========================================================================*/
int ObjectMoleculeRemoveBonds(ObjectMolecule *I,int sele0,int sele1)
{
  int a0,a1;
  int offset=0;
  int *b0,*b1;
  int both;
  int s;
  int a;

  offset=0;
  b0=I->Bond;
  b1=I->Bond;
  for(a=0;a<I->NBond;a++) {
    a0=b0[0];
    a1=b0[1];
    
    both=0;
    s=I->AtomInfo[a0].selEntry;
    if(SelectorIsMember(s,sele0))
      both++;
    s=I->AtomInfo[a1].selEntry;
    if(SelectorIsMember(s,sele1))
      both++;
    if(both<2) { /* reverse combo */
      both=0;
      s=I->AtomInfo[a1].selEntry;
      if(SelectorIsMember(s,sele0))
        both++;
      s=I->AtomInfo[a0].selEntry;
      if(SelectorIsMember(s,sele1))
        both++;
    }
    
    if(both==2) {
      offset--;
      b0+=3;
      I->AtomInfo[a0].chemFlag=false;
      I->AtomInfo[a1].chemFlag=false;
    } else if(offset) {
      *(b1++)=*(b0++); /* copy bond info */
      *(b1++)=*(b0++);
      *(b1++)=*(b0++);
    } else {
      *(b1++)=*(b0++); /* copy bond info */
      *(b1++)=*(b0++);
      *(b1++)=*(b0++);
    }
  }
  if(offset) {
    I->NBond += offset;
    VLASize(I->Bond,int,3*I->NBond);
    ObjectMoleculeInvalidate(I,cRepLine,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepCyl,cRepInvBonds);
    ObjectMoleculeInvalidate(I,cRepNonbonded,cRepInvBonds);
  }

  return(-offset);
}
/*========================================================================*/
void ObjectMoleculePurge(ObjectMolecule *I)
{
  int a,a0,a1;
  int *oldToNew = NULL;
  int offset=0;
  int *b0,*b1;
  AtomInfoType *ai0,*ai1;

  SelectorDelete(I->Obj.Name); /* remove the object selection and free up any selection entries*/

  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a]) 
      CoordSetPurge(I->CSet[a]);

  oldToNew = Alloc(int,I->NAtom);
  ai0=I->AtomInfo;
  ai1=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(ai0->deleteFlag) {
      offset--;
      ai0++;
      oldToNew[a]=-1;
    } else if(offset) {
      *(ai1++)=*(ai0++);
      oldToNew[a]=a+offset;
    } else {
      oldToNew[a]=a;
      ai0++;
      ai1++;
    }
  }
  if(offset) {
    I->NAtom += offset;
    VLASize(I->AtomInfo,AtomInfoType,I->NAtom);
    for(a=0;a<I->NCSet;a++)
      if(I->CSet[a])
        CoordSetAdjustAtmIdx(I->CSet[a],oldToNew,I->NAtom);
  }
  
  offset=0;
  b0=I->Bond;
  b1=I->Bond;
  for(a=0;a<I->NBond;a++) {
    a0=b0[0];
    a1=b0[1];
    if((oldToNew[a0]<0)||(oldToNew[a1]<0)) {
      offset--;
      b0+=3;
    } else if(offset) {
      *(b1++)=oldToNew[*(b0++)]; /* copy bond info */
      *(b1++)=oldToNew[*(b0++)];
      *(b1++)=*(b0++);
    } else {
      *(b1++)=oldToNew[*(b0++)]; /* copy bond info */
      *(b1++)=oldToNew[*(b0++)];
      *(b1++)=*(b0++);
    }
  }
  if(offset) {
    I->NBond += offset;
    VLASize(I->Bond,int,3*I->NBond);
  }
  FreeP(oldToNew);
    
  ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms);

}
/*========================================================================*/
int ObjectMoleculeGetAtomGeometry(ObjectMolecule *I,int state,int at)
{
  /* this determines hybridization from coordinates in those few cases
   * where it is unambiguous */

  int result = -1;
  int n,nn;
  float v0[3],v1[3],v2[3],v3[3];
  float d1[3],d2[3],d3[3];
  float cp1[3],cp2[3],cp3[3];
  float avg;
  float dp;
  n  = I->Neighbor[at];
  nn = I->Neighbor[n++]; /* get count */
  if(nn==4) 
    result = cAtomInfoTetrahedral; 
  else if(nn==3) {
    /* check cross products */
    ObjectMoleculeGetAtomVertex(I,state,at,v0);    
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n],v1);
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n+2],v2);
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n+4],v3);
    subtract3f(v1,v0,d1);
    subtract3f(v2,v0,d2);
    subtract3f(v3,v0,d3);
    cross_product3f(d1,d2,cp1);
    cross_product3f(d2,d3,cp2);
    cross_product3f(d3,d1,cp3);
    normalize3f(cp1);
    normalize3f(cp2);
    normalize3f(cp3);
    avg=(dot_product3f(cp1,cp2)+
         dot_product3f(cp2,cp3)+
         dot_product3f(cp3,cp1))/3.0;
    if(avg>0.75)
      result=cAtomInfoPlaner;
    else
      result=cAtomInfoTetrahedral;
  } else if(nn==2) {
    ObjectMoleculeGetAtomVertex(I,state,at,v0);    
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n],v1);
    ObjectMoleculeGetAtomVertex(I,state,I->Neighbor[n+2],v2);
    subtract3f(v1,v0,d1);
    subtract3f(v2,v0,d2);
    normalize3f(d1);
    normalize3f(d2);
    dp = dot_product3f(d1,d2);
    if(dp<-0.75)
      result=cAtomInfoLinear;
  }
  return(result);
}
/*========================================================================*/
void ObjectMoleculeInferChemForProtein(ObjectMolecule *I,int state)
{
  /* Infers chemical relations for a molecules under protein assumptions.
   * 
   * NOTE: this routine needs an all-atom model (with hydrogens!)
   * and it will make mistakes on non-protein atoms (if they haven't
   * already been assigned)
  */

  int a,n,a0,a1,nn;
  int changedFlag = true;
  
  AtomInfoType *ai,*ai0,*ai1;
  
  ObjectMoleculeUpdateNeighbors(I);

  /* first, try to find all amids and acids */
  while(changedFlag) {
    changedFlag=false;
    for(a=0;a<I->NAtom;a++) {
      ai=I->AtomInfo+a;
      if(ai->chemFlag) {
        if(ai->geom==cAtomInfoPlaner)
          if(ai->protons == cAN_C) {
            n = I->Neighbor[a];
            nn = I->Neighbor[n++];
            if(nn>1) {
              a1 = -1;
              while(1) {
                a0 = I->Neighbor[n];
                n+=2;
                if(a0<0) break;
                ai0 = I->AtomInfo+a0;
                if((ai0->protons==cAN_O)&&(!ai0->chemFlag)) {
                  a1=a0;
                  ai1=ai0; /* found candidate carbonyl */
                  break;
                }
              }
              if(a1>0) {
                n = I->Neighbor[a]+1;
                while(1) {
                  a0 = I->Neighbor[n];
                  n+=2;
                  if(a0<0) break;
                  if(a0!=a1) {
                    ai0 = I->AtomInfo+a0;
                    if(ai0->protons==cAN_O) {
                      if(!ai0->chemFlag) {
                        ai0->chemFlag=true; /* acid */
                        ai0->geom=cAtomInfoPlaner;
                        ai0->valence=1;
                        ai1->chemFlag=true;
                        ai1->geom=cAtomInfoPlaner;
                        ai1->valence=1;
                        changedFlag=true;
                        break;
                      }
                    } else if(ai0->protons==cAN_N) {
                      if(!ai0->chemFlag) { 
                        ai0->chemFlag=true; /* amide N */ 
                        ai0->geom=cAtomInfoPlaner;                            
                        ai0->valence=3;
                        ai1->chemFlag=true; /* amide =O */ 
                        ai1->geom=cAtomInfoPlaner;
                        ai1->valence=1;
                        changedFlag=true;
                        break;
                      } else if(ai0->geom==cAtomInfoPlaner) {
                        ai1->chemFlag=true; /* amide =O */
                        ai1->geom=cAtomInfoPlaner;
                        ai1->valence=1;
                        changedFlag=true;
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
      }
    }
  }
  /* then handle aldehydes and amines (partial amides - both missing a valence) */
  
  changedFlag=true;
  while(changedFlag) {
    changedFlag=false;
    for(a=0;a<I->NAtom;a++) {
      ai=I->AtomInfo+a;
      if(!ai->chemFlag) {
        if(ai->protons==cAN_C) {
          n = I->Neighbor[a];
          nn = I->Neighbor[n++];
          if(nn>1) {
            a1 = -1;
            while(1) {
              a0 = I->Neighbor[n];
              n+=2;
              if(a0<0) break;
              ai0 = I->AtomInfo+a0;
              if((ai0->protons==cAN_O)&&(!ai0->chemFlag)) { /* =O */
                ai->chemFlag=true; 
                ai->geom=cAtomInfoPlaner;
                ai->valence=1;
                ai0->chemFlag=true;
                ai0->geom=cAtomInfoPlaner;
                ai0->valence=3;
                changedFlag=true;
                break;
              }
            }
          }
        }
        else if(ai->protons==cAN_N)
          {
            if((!ai->chemFlag)||ai->geom!=cAtomInfoLinear) {
              if(ai->formalCharge==0.0) {
                ai->chemFlag=true; 
                ai->geom=cAtomInfoPlaner;
                ai->valence=3;
              }
            }
          }
      }
    }
  }

}
/*========================================================================*/
void ObjectMoleculeInferChemFromNeighGeom(ObjectMolecule *I,int state)
{
  /* infers chemical relations from neighbors and geometry 
  * NOTE: very limited in scope */

  int a,n,a0,nn;
  int changedFlag=true;
  int geom;
  int carbonVal[10];
  
  AtomInfoType *ai,*ai2;

  carbonVal[cAtomInfoTetrahedral] = 4;
  carbonVal[cAtomInfoPlaner] = 3;
  carbonVal[cAtomInfoLinear] = 2;
  
  ObjectMoleculeUpdateNeighbors(I);
  while(changedFlag) {
    changedFlag=false;
    for(a=0;a<I->NAtom;a++) {
      ai=I->AtomInfo+a;
      if(!ai->chemFlag) {
        geom=ObjectMoleculeGetAtomGeometry(I,state,a);
        switch(ai->protons) {
        case cAN_K:
          ai->chemFlag=1;
          ai->geom=cAtomInfoNone;
          ai->valence=0;
          break;
        case cAN_H:
        case cAN_F:
        case cAN_I:
        case cAN_Br:
          ai->chemFlag=1;
          ai->geom=cAtomInfoSingle;
          ai->valence=1;
          break;
        case cAN_O:
          n = I->Neighbor[a];
          nn = I->Neighbor[n++];
          if(nn!=1) { /* water, hydroxy, ether */
            ai->chemFlag=1;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=2;
          } else { /* hydroxy or carbonyl? check carbon geometry */
            a0 = I->Neighbor[n+2];
            ai2=I->AtomInfo+a0;
            if(ai2->chemFlag) {
              if((ai2->geom==cAtomInfoTetrahedral)||
                 (ai2->geom==cAtomInfoLinear)) {
                ai->chemFlag=1; /* hydroxy */
                ai->geom=cAtomInfoTetrahedral;
                ai->valence=2;
              }
            }
          }
          break;
        case cAN_C:
          if(geom>=0) {
            ai->geom = geom;
            ai->valence = carbonVal[geom];
            ai->chemFlag=true;
          } else {
            n = I->Neighbor[a];
            nn = I->Neighbor[n++];
            if(nn==1) { /* only one neighbor */
              ai2=I->AtomInfo+I->Neighbor[n];
              if(ai2->chemFlag&&(ai2->geom==cAtomInfoTetrahedral)) {
                ai->chemFlag=true; /* singleton carbon bonded to tetC must be tetC */
                ai->geom=cAtomInfoTetrahedral;
                ai->valence=4;
              }
            }
          }
          break;
        case cAN_N:
          if(geom==cAtomInfoPlaner) {
            ai->chemFlag=true;
            ai->geom=cAtomInfoPlaner;
            ai->valence=3;
          } else if(geom==cAtomInfoTetrahedral) {
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=4;
          }
          break;
        case cAN_S:
          n = I->Neighbor[a];
          nn = I->Neighbor[n++];
          if(nn==4) { /* sulfone */
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=4;
          } else if(nn==3) { /* suloxide */
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=3;
          } else if(nn==2) { /* thioether */
            ai->chemFlag=true;
            ai->geom=cAtomInfoTetrahedral;
            ai->valence=2;
          }
          break;
        case cAN_Cl:
          ai->chemFlag=1;
          if(ai->formalCharge==0.0) {
            ai->geom=cAtomInfoSingle;
            ai->valence=1;
          } else {
            ai->geom=cAtomInfoNone;
            ai->valence=0;
          }
          break;
        }
        if(ai->chemFlag)
          changedFlag=true;
      }
    }
  }
}
/*========================================================================*/
void ObjectMoleculeInferChemFromBonds(ObjectMolecule *I,int state)
{

  int a,b,*b0;
  AtomInfoType *ai,*ai0,*ai1;
  int a0,a1;
  int expect,order;
  int n,nn;
  int changedFlag;
  /* initialize accumulators on uncategorized atoms */


  ObjectMoleculeUpdateNeighbors(I);
  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(!ai->chemFlag) {
      ai->geom=0;
      ai->valence=0;
    }
    ai++;
  }
  
  /* find maximum bond order for each atom */

  b0=I->Bond;
  for(b=0;b<I->NBond;b++) {
    a0 = *(b0++);
    a1 = *(b0++);
    ai0=I->AtomInfo + a0;
    ai1=I->AtomInfo + a1;
    order = *(b0++);
    if(!ai0->chemFlag) {
      if(order>ai0->geom)
        ai0->geom=order;
      ai0->valence+=order;
    }
    if(!ai1->chemFlag) {
      if(order>ai1->geom)
        ai1->geom=order;
      ai1->valence+=order;
    }
    if(order==3) { 
      /* override existing chemistry * this is a temp fix to a pressing problem...
         we need to rethink the chemisty assignment ordering (should bond
         information come first? */
      ai0->geom = cAtomInfoLinear;
      ai1->geom = cAtomInfoLinear;
      switch(ai0->protons) {
      case cAN_C:
        ai0->valence=2;
        break;
      default:
        ai0->valence=1;
      }
      switch(ai1->protons) {
      case cAN_C:
        ai1->valence=2;
        break;
      default:
        ai1->valence=1;
      }
      ai0->chemFlag=true;
      ai1->chemFlag=true;
    }
  }

  /* now set up valences and geometries */

  ai=I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(!ai->chemFlag) {
      expect = AtomInfoGetExpectedValence(ai);
      n = I->Neighbor[a];
      nn = I->Neighbor[n++];
      if(ai->geom==3) {
        ai->geom = cAtomInfoLinear;
        switch(ai->protons) {
        case cAN_C:
          ai->valence=2;
          break;
        default:
          ai->valence=1;
        }
        ai->chemFlag=true;
      } else {
      if(expect<0) 
        expect = -expect; /* for now, just ignore this issue */
      /*      printf("%d %d %d %d\n",ai->geom,ai->valence,nn,expect);*/
      if(ai->valence==expect) { /* sum of bond orders equals valence */
        ai->chemFlag=true;
        ai->valence=nn;
        switch(ai->geom) { /* max bond order observed */
        case 0: ai->geom = cAtomInfoNone; break;
        case 2: ai->geom = cAtomInfoPlaner; break;
        case 3: ai->geom = cAtomInfoLinear; break;
        default: 
          if(expect==1) 
            ai->geom = cAtomInfoSingle;
          else
            ai->geom = cAtomInfoTetrahedral; 
          break;            
        }
      } else if(ai->valence<expect) { /* missing a bond */
        ai->chemFlag=true;
        ai->valence=nn+(expect-ai->valence); 
        switch(ai->geom) { 
        case 2: ai->geom = cAtomInfoPlaner; break;
        case 3: ai->geom = cAtomInfoLinear; break;
        default: 
          if(expect==1) 
            ai->geom = cAtomInfoSingle;
          else
            ai->geom = cAtomInfoTetrahedral; 
          break;
        }
      } else if(ai->valence>expect) {
        ai->chemFlag=true;
        ai->valence=nn;
        switch(ai->geom) { 
        case 2: ai->geom = cAtomInfoPlaner; break;
        case 3: ai->geom = cAtomInfoLinear; break;
        default: 
          if(expect==1) 
            ai->geom = cAtomInfoSingle;
          else
            ai->geom = cAtomInfoTetrahedral; 
          break;
        }
        if(nn>3)
          ai->geom = cAtomInfoTetrahedral;
      }
      }
    }
    ai++;
  }

  /* now go through and make sure conjugated amines are planer */
  changedFlag=true;
  while(changedFlag) {
    changedFlag=false;
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(ai->chemFlag) {
        if(ai->protons==cAN_N) 
          if(ai->formalCharge==0) 
            if(ai->geom==cAtomInfoTetrahedral) { 
              /* search for uncharged tetrahedral nitrogen */
              n = I->Neighbor[a]+1;
              while(1) {
                a0 = I->Neighbor[n];
                n+=2;
                if(a0<0) break;
                ai0 = I->AtomInfo+a0;
                if((ai0->chemFlag)&&(ai0->geom==cAtomInfoPlaner)&&
                   ((ai0->protons==cAN_C)||(ai1->protons==cAN_N))) {
                  ai->geom=cAtomInfoPlaner; /* found probable delocalization */
                  ai->valence=3; /* just in case...*/
                  changedFlag=true;
                  break;
                }
              }
            }
      }
      ai++; 
    }
  }

  /* now go through and make sure conjugated anions are planer */
  changedFlag=true;
  while(changedFlag) {
    changedFlag=false;
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(ai->chemFlag) {
        if(ai->protons==cAN_O) 
          if(ai->formalCharge==-1) 
            if((ai->geom==cAtomInfoTetrahedral)||
               (ai->geom==cAtomInfoSingle)) { /* search for anionic tetrahedral oxygen */
              n = I->Neighbor[a]+1;
              while(1) {
                a0 = I->Neighbor[n++];
                if(a0<0) break;
                ai0 = I->AtomInfo+a0;
                if((ai0->chemFlag)&&(ai0->geom==cAtomInfoPlaner)&&
                   ((ai0->protons==cAN_C)||(ai1->protons==cAN_N))) {
                  ai->geom=cAtomInfoPlaner; /* found probable delocalization */
                  changedFlag=true;
                  break;
                }
              }
            }
      }
      ai++;
    }
  }


}
/*========================================================================*/
void ObjectMoleculeTransformSelection(ObjectMolecule *I,int state,int sele,float *TTT) 
{
  /* if sele == -1, then the whole object state is transformed */

  int a,s;
  int flag=false;
  CoordSet *cs;
  AtomInfoType *ai;

  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    if(sele>=0) {
      ai=I->AtomInfo;
      for(a=0;a<I->NAtom;a++) {
        s=ai->selEntry;
        if(!ai->protected)
          if(SelectorIsMember(s,sele))
            {
              CoordSetTransformAtom(cs,a,TTT);
              flag=true;
            }
        ai++;
      }
    } else {
      ai=I->AtomInfo;
      for(a=0;a<I->NAtom;a++) {
        if(!ai->protected)
          CoordSetTransformAtom(cs,a,TTT);
        ai++;
      }
      flag=true;
    }
    if(flag) 
      cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
  }
}
/*========================================================================*/
int ObjectMoleculeGetAtomIndex(ObjectMolecule *I,int sele)
{
  int a,s;
  for(a=0;a<I->NAtom;a++) {
    s=I->AtomInfo[a].selEntry;
    if(SelectorIsMember(s,sele))
      return(a);
  }
  return(-1);
}
/*========================================================================*/
void ObjectMoleculeUpdateNonbonded(ObjectMolecule *I)
{
  int a,*b;
  AtomInfoType *ai;

  if(!I->DiscreteFlag) {
    ai=I->AtomInfo;
    
    for(a=0;a<I->NAtom;a++)
      (ai++)->bonded = false;
    
    b=I->Bond;
    ai=I->AtomInfo;
    for(a=0;a<I->NBond;a++)
      {
        ai[*(b++)].bonded=true;
        ai[*(b++)].bonded=true;
        b++;
      }

  }
}
/*========================================================================*/
void ObjectMoleculeUpdateNeighbors(ObjectMolecule *I)
{
  /* neighbor storage structure:
     
     0       list offset for atom 0 = n
     1       list offset for atom 1 = n + m + 1
     ...
     n-1     list offset for atom n-1

     n       count for atom 0 
     n+1     neighbor of atom 0
     n+2     bond index
     n+3     neighbor of atom 1
     n+4     bond index
     ...
     n+m     -1

     n+m+1   count for atom 1
     n+m+2   neighbor of atom 1
     n+m+3   bond index
     n+m+4   neighbor of atom 1
     n+m+5   bond index
     etc.
 */

  int size;
  int a,b,c,d,*l,l0,l1;
  if(!I->Neighbor) {
    size = (I->NAtom*3)+(I->NBond*4); 
    if(I->Neighbor) {
      VLACheck(I->Neighbor,int,size);
    } else {
      I->Neighbor=VLAlloc(int,size);
    }
    
    /* initialize */
    l = I->Neighbor;
    for(a=0;a<I->NAtom;a++)
      (*l++)=0;
    
    /* count neighbors for each atom */
    l = I->Bond;
    for(b=0;b<I->NBond;b++) {
      I->Neighbor[(*l++)]++;
      I->Neighbor[(*l++)]++;
      l++;
    }
    
    /* set up offsets and list terminators */
    c = I-> NAtom;
    for(a=0;a<I->NAtom;a++) {
      d = I->Neighbor[a];
      I->Neighbor[c]=d; /* store neighbor count */
      I->Neighbor[a]=c+d+d+1; /* set initial position to end of list, we'll fill backwards */
      I->Neighbor[I->Neighbor[a]]=-1; /* store terminator */
      c += d + d + 2;
    }
    
    /* now load neighbors in a sequential list for each atom (reverse order) */
    l = I->Bond;
    for(b=0;b<I->NBond;b++) {
      l0 = *(l++);
      l1 = *(l++);
      l++;

      I->Neighbor[l0]--;
      I->Neighbor[I->Neighbor[l0]]=b;
      I->Neighbor[l1]--;
      I->Neighbor[I->Neighbor[l1]]=b;

      I->Neighbor[l0]--;
      I->Neighbor[I->Neighbor[l0]]=l1;
      I->Neighbor[l1]--;
      I->Neighbor[I->Neighbor[l1]]=l0;
    }
    for(a=0;a<I->NAtom;a++) { /* adjust down to point to the count, not the first entry */
      if(I->Neighbor[a]>=0)
        I->Neighbor[a]--;
    }
    l=I->Neighbor;
  }
}
/*========================================================================*/
CoordSet *ObjectMoleculeChemPyModel2CoordSet(PyObject *model,AtomInfoType **atInfoPtr)
{
  int nAtom,nBond;
  int a,c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  float *f;
  int *ii,*bond=NULL;
  int ok=true;
  int autoshow_lines;
  int autoshow_nonbonded;
  int hetatm;

  PyObject *atomList = NULL;
  PyObject *bondList = NULL;
  PyObject *atom = NULL;
  PyObject *bnd = NULL;
  PyObject *index = NULL;
  PyObject *crd = NULL;
  PyObject *tmp = NULL;
  autoshow_lines = SettingGet(cSetting_autoshow_lines);
  autoshow_nonbonded = SettingGet(cSetting_autoshow_nonbonded);
  AtomInfoPrimeColors();

  nAtom=0;
  nBond=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  atomList = PyObject_GetAttrString(model,"atom");
  if(atomList) 
    nAtom = PyList_Size(atomList);
  else 
    ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get atom list");


  if(ok) {
	 coord=VLAlloc(float,3*nAtom);
	 if(atInfo)
		VLACheck(atInfo,AtomInfoType,nAtom);	 
  }

  if(ok) { 
    
	 f=coord;
	 for(a=0;a<nAtom;a++)
		{
        atom = PyList_GetItem(atomList,a);
        if(!atom) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get atom");
        crd = PyObject_GetAttrString(atom,"coord");
        if(!crd) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get coordinates");
        else {
          for(c=0;c<3;c++) {
            tmp = PyList_GetItem(crd,c);
            if (tmp) 
              ok = PConvPyObjectToFloat(tmp,f++);
            if(!ok) {
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read coordinates");
              break;
            }
          }
        }
        Py_XDECREF(crd);
        
        ai = atInfo+a;
        ai->id = a; /* chempy models are zero-based */

        if(ok) {
          tmp = PyObject_GetAttrString(atom,"name");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->name,sizeof(AtomName)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read name");
          Py_XDECREF(tmp);
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","text_type")) { 
            tmp = PyObject_GetAttrString(atom,"text_type");
            if (tmp)
              ok = PConvPyObjectToStrMaxClean(tmp,ai->textType,sizeof(TextType)-1);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read text_type");
            Py_XDECREF(tmp);
          } else {
            ai->textType[0]=0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","vdw")) { 
            tmp = PyObject_GetAttrString(atom,"vdw");
            if (tmp)
              ok = PConvPyObjectToFloat(tmp,&ai->vdw);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read vdw radius");
            Py_XDECREF(tmp);
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","numeric_type")) { 
            tmp = PyObject_GetAttrString(atom,"numeric_type");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&ai->customType);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read numeric_type");
            Py_XDECREF(tmp);
          } else {
            ai->customType = cAtomInfoNoType;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","formal_charge")) { 
            tmp = PyObject_GetAttrString(atom,"formal_charge");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&ai->formalCharge);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read formal_charge");
            Py_XDECREF(tmp);
          } else {
            ai->formalCharge = 0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","partial_charge")) { 
            tmp = PyObject_GetAttrString(atom,"partial_charge");
            if (tmp)
              ok = PConvPyObjectToFloat(tmp,&ai->partialCharge);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read partial_charge");
            Py_XDECREF(tmp);
          } else {
            ai->partialCharge = 0.0;
          }
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","flags")) {         
            tmp = PyObject_GetAttrString(atom,"flags");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,(int*)&ai->flags);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read flags");
            Py_XDECREF(tmp);
          } else {
            ai->flags = 0;
          }
        }

        if(ok) {
          tmp = PyObject_GetAttrString(atom,"resn");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->resn,sizeof(ResName)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read resn");
          Py_XDECREF(tmp);
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"resi");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->resi,sizeof(ResIdent)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read resi");
          else
            sscanf(ai->resi,"%d",&ai->resv);
          Py_XDECREF(tmp);
        }

        if(ok) {
          if(PTruthCallStr(atom,"has","resi_number")) {         
            tmp = PyObject_GetAttrString(atom,"resi_number");
            if (tmp)
              ok = PConvPyObjectToInt(tmp,&ai->resv);
            if(!ok) 
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read resi_number");
            Py_XDECREF(tmp);
          }
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"segi");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->segi,sizeof(SegIdent)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read segi");
          Py_XDECREF(tmp);
        }

		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"b");
          if (tmp)
            ok = PConvPyObjectToFloat(tmp,&ai->b);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read b value");
          Py_XDECREF(tmp);
        }

		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"q");
          if (tmp)
            ok = PConvPyObjectToFloat(tmp,&ai->q);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read occupancy");
          Py_XDECREF(tmp);
        }

        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"chain");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->chain,sizeof(Chain)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read chain");
          Py_XDECREF(tmp);
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"hetatm");
          if (tmp)
            ok = PConvPyObjectToInt(tmp,&hetatm);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read hetatm");
          else
            ai->hetatm = hetatm;
          Py_XDECREF(tmp);
        }
        
		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"alt");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->alt,sizeof(Chain)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read alternate conformation");
          Py_XDECREF(tmp);
        }

		  if(ok) {
          tmp = PyObject_GetAttrString(atom,"symbol");
          if (tmp)
            ok = PConvPyObjectToStrMaxClean(tmp,ai->elem,sizeof(AtomName)-1);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read symbol");
          Py_XDECREF(tmp);
        }
        
        for(c=0;c<cRepCnt;c++) {
          atInfo[a].visRep[c] = false;
		  }
        atInfo[a].visRep[cRepLine] = autoshow_lines; /* show lines by default */
        atInfo[a].visRep[cRepNonbonded] = autoshow_nonbonded; /* show lines by default */

		  if(ok&&atInfo) {
			 AtomInfoAssignParameters(ai);
			 atInfo[a].color=AtomInfoGetColor(ai);
		  }

		  if(!ok)
			 break;
		}
  }

  bondList = PyObject_GetAttrString(model,"bond");
  if(bondList) 
    nBond = PyList_Size(bondList);
  else
    ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get bond list");

  if(ok) {
	 bond=VLAlloc(int,3*nBond);
    ii=bond;
	 for(a=0;a<nBond;a++)
		{
        bnd = PyList_GetItem(bondList,a);
        if(!bnd) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get bond");
        index = PyObject_GetAttrString(bnd,"index");
        if(!index) 
          ok=ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't get bond indices");
        else {
          for(c=0;c<2;c++) {
            tmp = PyList_GetItem(index,c);
            if (tmp) 
              ok = PConvPyObjectToInt(tmp,ii++);
            if(!ok) {
              ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read coordinates");
              break;
            }
          }
        }
        if(ok) {
          tmp = PyObject_GetAttrString(bnd,"order");
          if (tmp)
            ok = PConvPyObjectToInt(tmp,ii++);
          if(!ok) 
            ErrMessage("ObjectMoleculeChemPyModel2CoordSet","can't read bond order");
          Py_XDECREF(tmp);
        }
        Py_XDECREF(index);
      }
  }

  Py_XDECREF(atomList);
  Py_XDECREF(bondList);

  if(ok) {
	 cset = CoordSetNew();
	 cset->NIndex=nAtom;
	 cset->Coord=coord;
	 cset->NTmpBond=nBond;
	 cset->TmpBond=bond;
  } else {
	 VLAFreeP(bond);
	 VLAFreeP(coord);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;

  if(PyErr_Occurred())
    PyErr_Print();
  return(cset);
}


/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadChemPyModel(ObjectMolecule *I,PyObject *model,int frame,int discrete)
{
  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  unsigned int nAtom = 0;
  PyObject *tmp,*mol;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
		atInfo = I->AtomInfo;
		isNew = true;
	 } else {
		atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
		isNew = false;
	 }
	 cset=ObjectMoleculeChemPyModel2CoordSet(model,&atInfo);	 

    mol = PyObject_GetAttrString(model,"molecule");
    if(mol) {
      if(PyObject_HasAttrString(mol,"title")) {
        tmp = PyObject_GetAttrString(mol,"title");
        if(tmp) {
          UtilNCopy(cset->Name,PyString_AsString(tmp),sizeof(WordType));
          Py_DECREF(tmp);
        }
      }
      Py_DECREF(mol);
    }
    if(PyObject_HasAttrString(model,"spheroid")&&
       PyObject_HasAttrString(model,"spheroid_normals"))
      {
        tmp = PyObject_GetAttrString(model,"spheroid");
        if(tmp) {
          cset->NSpheroid = PConvPyListToFloatArray(tmp,&cset->Spheroid);
          Py_DECREF(tmp);
        }
        tmp = PyObject_GetAttrString(model,"spheroid_normals");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&cset->SpheroidNormal);
          Py_DECREF(tmp);
        }
      }
    mol = PyObject_GetAttrString(model,"molecule");
    
	 nAtom=cset->NIndex;
  }
  /* include coordinate set */
  if(ok) {
    cset->Obj = I;
    cset->fEnumIndices(cset);
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
    if(isNew) {	
      I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
    } else {
      ObjectMoleculeMerge(I,atInfo,cset,true); /* NOTE: will release atInfo */
    }
    if(isNew) I->NAtom=nAtom;
    if(frame<0) frame=I->NCSet;
    VLACheck(I->CSet,CoordSet*,frame);
    if(I->NCSet<=frame) I->NCSet=frame+1;
    if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
    I->CSet[frame] = cset;
    if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,0.2,false);
    if(cset->TmpSymmetry&&(!I->Symmetry)) {
      I->Symmetry=cset->TmpSymmetry;
      cset->TmpSymmetry=NULL;
      SymmetryAttemptGeneration(I->Symmetry);
    }
    SceneCountFrames();
    ObjectMoleculeExtendIndices(I);
    ObjectMoleculeSort(I);
    ObjectMoleculeUpdateNonbonded(I);
  }
  return(I);
}


/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadCoords(ObjectMolecule *I,PyObject *coords,int frame)
{
  CoordSet *cset = NULL;
  int ok=true;
  int a,l;
  PyObject *v;
  float *f;
  a=0;
  while(a<I->NCSet) {
    if(I->CSet[a]) {
      cset=I->CSet[a];
      break;
    }
    a++;
  }
  
  if(!PyList_Check(coords)) 
    ErrMessage("LoadsCoords","passed argument is not a list");
  else {
    l = PyList_Size(coords);
    if (l==cset->NIndex) {
      cset=CoordSetCopy(cset);
      f=cset->Coord;
      for(a=0;a<l;a++) {
        v=PyList_GetItem(coords,a);
/* no error checking */
        *(f++)=PyFloat_AsDouble(PyList_GetItem(v,0)); 
        *(f++)=PyFloat_AsDouble(PyList_GetItem(v,1));
        *(f++)=PyFloat_AsDouble(PyList_GetItem(v,2));
      }
    }
  }
  /* include coordinate set */
  if(ok) {
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,cRepAll,cRepInvRep);

    if(frame<0) frame=I->NCSet;
    VLACheck(I->CSet,CoordSet*,frame);
    if(I->NCSet<=frame) I->NCSet=frame+1;
    if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
    I->CSet[frame] = cset;
    SceneCountFrames();
  }
  return(I);
}

/*========================================================================*/
static int BondInOrder(int *a,int b1,int b2)
{
  return(BondCompare(a+(b1*3),a+(b2*3))<=0);
}
/*========================================================================*/
static int BondCompare(int *a,int *b)
{
  int result;
  if(a[0]==b[0]) {
	if(a[1]==b[1]) {
	  result=0;
	} else if(a[1]>b[1]) {
	  result=1;
	} else {
	  result=-1;
	}
  } else if(a[0]>b[0]) {
	result=1;
  } else {
	result=-1;
  }
  return(result);
}
/*========================================================================*/
void ObjectMoleculeBlindSymMovie(ObjectMolecule *I)
{
  CoordSet *frac;
  int a,c;
  int x,y,z;
  float m[16];

  if(I->NCSet!=1) {
    ErrMessage("ObjectMolecule:","SymMovie only works on objects with a single state.");
  } else if(!I->Symmetry) {
    ErrMessage("ObjectMolecule:","No symmetry loaded!");
  } else if(!I->Symmetry->NSymMat) {
    ErrMessage("ObjectMolecule:","No symmetry matrices!");    
  } else if(I->CSet[0]) {
    frac = CoordSetCopy(I->CSet[0]);
    CoordSetRealToFrac(frac,I->Symmetry->Crystal);
    for(x=-1;x<2;x++)
      for(y=-1;y<2;y++)
        for(z=-1;z<2;z++)
          for(a=0;a<I->Symmetry->NSymMat;a++) {
            if(!((!a)&&(!x)&&(!y)&&(!z))) {
              c = I->NCSet;
              VLACheck(I->CSet,CoordSet*,c);
              I->CSet[c] = CoordSetCopy(frac);
              CoordSetTransform44f(I->CSet[c],I->Symmetry->SymMatVLA+(a*16));
              identity44f(m);
              m[3] = x;
              m[7] = y;
              m[11] = z;
              CoordSetTransform44f(I->CSet[c],m);
              CoordSetFracToReal(I->CSet[c],I->Symmetry->Crystal);
              I->NCSet++;
            }
          }
    frac->fFree(frac);
  }
  SceneChanged();
}

/*========================================================================*/
void ObjectMoleculeExtendIndices(ObjectMolecule *I)
{
  int a;
  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a])
      if(I->CSet[a]->fExtendIndices)
        I->CSet[a]->fExtendIndices(I->CSet[a],I->NAtom);
}
/*========================================================================*/
void ObjectMoleculeSort(ObjectMolecule *I) /* sorts atoms and bonds */
{
  int *index,*outdex;
  int a,b;
  CoordSet *cs,**dcs;
  AtomInfoType *atInfo;
  int *dAtmToIdx;

  if(!I->DiscreteFlag) {

    index=AtomInfoGetSortedIndex(I->AtomInfo,I->NAtom,&outdex);
    for(a=0;a<I->NBond;a++) { /* bonds */
      I->Bond[a*3]=outdex[I->Bond[a*3]];
      I->Bond[a*3+1]=outdex[I->Bond[a*3+1]];
    }
    
    for(a=0;a<I->NCSet;a++) { /* coordinate set mapping */
      cs=I->CSet[a];
      if(cs) {
        for(b=0;b<cs->NIndex;b++)
          cs->IdxToAtm[b]=outdex[cs->IdxToAtm[b]];
        if(cs->AtmToIdx) {
          for(b=0;b<I->NAtom;b++)
            cs->AtmToIdx[b]=-1;
          for(b=0;b<cs->NIndex;b++) {
            cs->AtmToIdx[cs->IdxToAtm[b]]=b;
          }
        }
      }
    }
    
    atInfo=(AtomInfoType*)VLAMalloc(I->NAtom,sizeof(AtomInfoType),5,true);
    /* autozero here is important */
    for(a=0;a<I->NAtom;a++)
      atInfo[a]=I->AtomInfo[index[a]];
    VLAFreeP(I->AtomInfo);
    I->AtomInfo=atInfo;
    
    if(I->DiscreteFlag) {
      dcs = VLAlloc(CoordSet*,I->NAtom);
      dAtmToIdx = VLAlloc(int,I->NAtom);
      for(a=0;a<I->NAtom;a++) {
        b=index[a];
        dcs[a] = I->DiscreteCSet[b];
        dAtmToIdx[a] = I->DiscreteAtmToIdx[b];
      }
      VLAFreeP(I->DiscreteCSet);
      VLAFreeP(I->DiscreteAtmToIdx);
      I->DiscreteCSet = dcs;
      I->DiscreteAtmToIdx = dAtmToIdx;
    }
    AtomInfoFreeSortedIndexes(index,outdex);

    UtilSortInPlace(I->Bond,I->NBond,sizeof(int)*3,(UtilOrderFn*)BondInOrder);
    /* sort...important! */
    ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms); /* important */

  }
}
/*========================================================================*/
CoordSet *ObjectMoleculeMOLStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom,nBond;
  int a,c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL;
  char cc[MAXLINELEN],resn[MAXLINELEN] = "UNK";
  float *f;
  int *ii,*bond=NULL;
  int ok=true;
  int autoshow_lines;
  int autoshow_nonbonded;
  WordType nameTmp;

  autoshow_lines = SettingGet(cSetting_autoshow_lines);
  autoshow_nonbonded = SettingGet(cSetting_autoshow_nonbonded);
  AtomInfoPrimeColors();

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  p=ParseWordCopy(nameTmp,p,sizeof(WordType)-1);
  p=nextline(p); 
  p=nextline(p);
  p=nextline(p);

  if(ok) {
	 p=ncopy(cc,p,3);
	 if(sscanf(cc,"%d",&nAtom)!=1)
		ok=ErrMessage("ReadMOLFile","bad atom count");
  }

  if(ok) {  
	 p=ncopy(cc,p,3);
	 if(sscanf(cc,"%d",&nBond)!=1)
		ok=ErrMessage("ReadMOLFile","bad bond count");
  }

  if(ok) {
	 coord=VLAlloc(float,3*nAtom);
	 if(atInfo)
		VLACheck(atInfo,AtomInfoType,nAtom);	 
  }
  
  p=nextline(p);

  /* read coordinates and atom names */

  if(ok) { 
	 f=coord;
	 for(a=0;a<nAtom;a++)
		{
		  if(ok) {
			 p=ncopy(cc,p,10);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMOLFile","bad coordinate");
		  }
		  if(ok) {
			 p=ncopy(cc,p,10);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMOLFile","bad coordinate");
		  }
		  if(ok) {
			 p=ncopy(cc,p,10);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMOLFile","bad coordinate");
		  }
		  if(ok) {
          p=nskip(p,1);
			 p=ncopy(atInfo[a].name,p,3);
			 UtilCleanStr(atInfo[a].name);
			 
          for(c=0;c<cRepCnt;c++) {
            atInfo[a].visRep[c] = false;
          }
          atInfo[a].visRep[cRepLine] = autoshow_lines; /* show lines by default */
          atInfo[a].visRep[cRepNonbonded] = autoshow_nonbonded; /* show lines by default */

		  }
		  if(ok&&atInfo) {
          atInfo[a].id = a+1;
			 strcpy(atInfo[a].resn,resn);
			 atInfo[a].hetatm=true;
			 AtomInfoAssignParameters(atInfo+a);
			 atInfo[a].color=AtomInfoGetColor(atInfo+a);
          atInfo[a].alt[0]=0;
          atInfo[a].segi[0]=0;
          atInfo[a].resi[0]=0;
		  }
		  p=nextline(p);
		  if(!ok)
			 break;
		}
  }
  if(ok) {
	 bond=VLAlloc(int,3*nBond);
	 ii=bond;
	 for(a=0;a<nBond;a++)
		{
		  if(ok) {
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",ii++)!=1)
				ok=ErrMessage("ReadMOLFile","bad bond atom");
		  }
		  
		  if(ok) {  
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",ii++)!=1)
				ok=ErrMessage("ReadMOLFile","bad bond atom");
		  }

		  if(ok) {  
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",ii++)!=1)
				ok=ErrMessage("ReadMOLFile","bad bond order");
		  }
		  if(!ok)
			 break;
		  p=nextline(p);
		}
	 ii=bond;
	 for(a=0;a<nBond;a++) {
		(*(ii++))--; /* adjust bond indexs down one */
		(*(ii++))--; 
      ii++;
	 }
  }
  if(ok) {
	 cset = CoordSetNew();
	 cset->NIndex=nAtom;
	 cset->Coord=coord;
	 cset->NTmpBond=nBond;
	 cset->TmpBond=bond;
    strcpy(cset->Name,nameTmp);
  } else {
	 VLAFreeP(bond);
	 VLAFreeP(coord);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  return(cset);
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadMOLStr(ObjectMolecule *I,char *MOLStr,int frame,int discrete)
{
  int ok = true;
  CoordSet *cset=NULL;
  AtomInfoType *atInfo;
  int isNew;
  int nAtom;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(isNew) {
    I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
    atInfo = I->AtomInfo;
    isNew = true;
  } else {
    atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
    isNew = false;
  }

  cset=ObjectMoleculeMOLStr2CoordSet(MOLStr,&atInfo);
  
  if(!cset) 
	 {
      ObjectMoleculeFree(I);
		ok=false;
	 }
  
  if(ok)
	 {
		if(frame<0)
		  frame=I->NCSet;
		if(I->NCSet<=frame)
		  I->NCSet=frame+1;
		VLACheck(I->CSet,CoordSet*,frame);
      
      nAtom=cset->NIndex;
      
      cset->Obj = I;
      cset->fEnumIndices(cset);
      if(cset->fInvalidateRep)
        cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
      if(isNew) {		
        I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
      } else {
        ObjectMoleculeMerge(I,atInfo,cset,true); /* NOTE: will release atInfo */
      }

      if(isNew) I->NAtom=nAtom;
      if(frame<0) frame=I->NCSet;
      VLACheck(I->CSet,CoordSet*,frame);
      if(I->NCSet<=frame) I->NCSet=frame+1;
      if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
      I->CSet[frame] = cset;
      
      if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,0.2,false);
      
      SceneCountFrames();
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeSort(I);
      ObjectMoleculeUpdateNonbonded(I);
	 }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadMOLFile(ObjectMolecule *obj,char *fname,int frame,int discrete)
{
  ObjectMolecule* I=NULL;
  int ok=true;
  FILE *f;
  fpos_t size;
  char *buffer,*p;

  f=fopen(fname,"r");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadMOLFile","Unable to open file!");
  else
	 {
		if(DebugState&DebugMolecule)
		  {
			 printf(" ObjectMoleculeLoadMOLFile: Loading from %s.\n",fname);
			 fflush(stdout);
		  }
		
		fseek(f,0,SEEK_END);
		fgetpos(f,&size);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);
		I=ObjectMoleculeReadMOLStr(obj,buffer,frame,discrete);
		mfree(buffer);
	 }

  return(I);
}

/*========================================================================*/
void ObjectMoleculeMerge(ObjectMolecule *I,AtomInfoType *ai,CoordSet *cs,int bondSearchFlag)
{
  int *index,*outdex,*a2i,*i2a,*bond=NULL;
  int a,b,c,lb,nb,ac,a1,a2;
  int found;
  int nAt,nBd,nBond;
  int expansionFlag = false;
  AtomInfoType *ai2;
  int oldNAtom,oldNBond;

  oldNAtom = I->NAtom;
  oldNBond = I->NBond;


  /* first, sort the coodinate set */
  
  index=AtomInfoGetSortedIndex(ai,cs->NIndex,&outdex);
  for(b=0;b<cs->NIndex;b++)
	 cs->IdxToAtm[b]=outdex[cs->IdxToAtm[b]];
  for(b=0;b<cs->NIndex;b++)
	 cs->AtmToIdx[b]=-1;
  for(b=0;b<cs->NIndex;b++)
	 cs->AtmToIdx[cs->IdxToAtm[b]]=b;
  ai2=(AtomInfoType*)VLAMalloc(cs->NIndex,sizeof(AtomInfoType),5,true); /* autozero here is important */
  for(a=0;a<cs->NIndex;a++) 
	 ai2[a]=ai[index[a]]; /* creates a sorted list of atom info records */
  VLAFreeP(ai);
  ai=ai2;

  /* now, match it up with the current object's atomic information */
	 
  for(a=0;a<cs->NIndex;a++) {
	 index[a]=-1;
	 outdex[a]=-1;
  }

  c=0;
  b=0;  
  for(a=0;a<cs->NIndex;a++) {
	 found=false;
    if(!I->DiscreteFlag) { /* don't even try matching for discrete objects */
      lb=b;
      while(b<I->NAtom) {
        ac=(AtomInfoCompare(ai+a,I->AtomInfo+b));
        if(!ac) {
          found=true;
          break;
        }
        else if(ac<0) {
          break;
        }
        b++;
      }
    }
	 if(found) {
		index[a]=b; /* store real atom index b for a in index[a] */
		b++;
	 } else {
	   index[a]=I->NAtom+c; /* otherwise, this is a new atom */
	   c++;
	   b=lb;
	 }
  }

  /* first, reassign atom info for matched atoms */

  /* allocate additional space */
  if(c)
	{
	  expansionFlag=true;
	  nAt=I->NAtom+c;
	} else {
     nAt=I->NAtom;
   }
  
  if(expansionFlag) {
	VLACheck(I->AtomInfo,AtomInfoType,nAt);
  }

  /* allocate our new x-ref tables */
  if(nAt<I->NAtom) nAt=I->NAtom;
  a2i = Alloc(int,nAt);
  i2a = Alloc(int,cs->NIndex);
  ErrChkPtr(a2i);
  ErrChkPtr(i2a);
  
  for(a=0;a<cs->NIndex;a++) /* a is in original file space */
    {
		a1=cs->IdxToAtm[a]; /* a1 is in sorted atom info space */
		a2=index[a1];
		i2a[a]=a2; /* a2 is in object space */
      if(a2<oldNAtom)
        AtomInfoCombine(I->AtomInfo+a2,ai+a1);
      else
        *(I->AtomInfo+a2)=*(ai+a1);
    }
  
  if(I->DiscreteFlag) {
    if(I->NDiscrete<nAt) {
      VLACheck(I->DiscreteAtmToIdx,int,nAt);
      VLACheck(I->DiscreteCSet,CoordSet*,nAt);    
      for(a=I->NDiscrete;a<nAt;a++) {
        I->DiscreteAtmToIdx[a]=-1;
        I->DiscreteCSet[a]=NULL;
      }
    }
    I->NDiscrete = nAt;
  }
  
  cs->NAtIndex = nAt;
  I->NAtom = nAt;
  
  FreeP(cs->AtmToIdx);
  FreeP(cs->IdxToAtm);
  cs->AtmToIdx = a2i;
  cs->IdxToAtm = i2a;

  if(I->DiscreteFlag) {
    FreeP(cs->AtmToIdx);
    for(a=0;a<cs->NIndex;a++) {
      I->DiscreteAtmToIdx[cs->IdxToAtm[a]]=a;
      I->DiscreteCSet[cs->IdxToAtm[a]] = cs;
    }
  } else {
    for(a=0;a<cs->NAtIndex;a++)
      cs->AtmToIdx[a]=-1;
    for(a=0;a<cs->NIndex;a++)
      cs->AtmToIdx[cs->IdxToAtm[a]]=a;
  }
  
  VLAFreeP(ai);
  AtomInfoFreeSortedIndexes(index,outdex);
  
  /* now find and integrate and any new bonds */
  if(expansionFlag) { /* expansion flag means we have introduced at least 1 new atom */
    nBond = ObjectMoleculeConnect(I,&bond,I->AtomInfo,cs,0.2,bondSearchFlag);
    if(nBond) {
      index=Alloc(int,nBond);
      
      c=0;
      b=0;  
      nb=0;
      for(a=0;a<nBond;a++) { /* iterate over new bonds */
        found=false;
        b=nb; /* pick up where we left off */
        while(b<I->NBond) { 
          ac=BondCompare(bond+a*3,I->Bond+b*3);
          if(!ac) { /* zero is a match */
            found=true;
            break;
          } else if(ac<0) { /* gone past position of this bond */
            break;
          }
          b++; /* no match yet, keep looking */
        }
        if(found) {
          index[a]=b; /* existing bond...*/
          nb=b+1;
        } else { /* this is a new bond, save index and increment */
          index[a]=I->NBond+c;
          c++; 
        }
      }
      /* first, reassign atom info for matched atoms */
      if(c) {
        /* allocate additional space */
        nBd=I->NBond+c;
        
        VLACheck(I->Bond,int,nBd*3);
        
        for(a=0;a<nBond;a++) /* copy the new bonds */
          {
            a2=index[a];
            if(a2 >= I->NBond) { 
              I->Bond[3*a2]=bond[3*a]; 
              I->Bond[3*a2+1]=bond[3*a+1]; 
              I->Bond[3*a2+2]=bond[3*a+2]; 
            }
          }
        I->NBond=nBd;
      }
      FreeP(index);
    }
    VLAFreeP(bond);
  }

  if(oldNAtom) {
    if(oldNAtom==I->NAtom) {
      if(oldNBond!=I->NBond) {
        ObjectMoleculeInvalidate(I,cRepAll,cRepInvBonds);
      }
    } else {
      ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms);
    }
  }

}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadPDBStr(ObjectMolecule *I,char *PDBStr,int frame,int discrete)
{
  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  unsigned int nAtom = 0;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
		atInfo = I->AtomInfo;
		isNew = true;
	 } else {
		atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
		isNew = false;
	 }
	 cset=ObjectMoleculePDBStr2CoordSet(PDBStr,&atInfo);	 
	 nAtom=cset->NIndex;
  }

  /* include coordinate set */
  if(ok) {
    cset->Obj = I;
    cset->fEnumIndices(cset);
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
    if(isNew) {		
      I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
    } else {
      ObjectMoleculeMerge(I,atInfo,cset,true); /* NOTE: will release atInfo */
    }
    if(isNew) I->NAtom=nAtom;
    if(frame<0) frame=I->NCSet;
    VLACheck(I->CSet,CoordSet*,frame);
    if(I->NCSet<=frame) I->NCSet=frame+1;
    if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
    I->CSet[frame] = cset;
    if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,0.2,true);
    if(cset->TmpSymmetry&&(!I->Symmetry)) {
      I->Symmetry=cset->TmpSymmetry;
      cset->TmpSymmetry=NULL;
      SymmetryAttemptGeneration(I->Symmetry);
    }
    SceneCountFrames();
    ObjectMoleculeExtendIndices(I);
    ObjectMoleculeSort(I);
    ObjectMoleculeUpdateNonbonded(I);
  }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadPDBFile(ObjectMolecule *obj,char *fname,int frame,int discrete)
{
  ObjectMolecule *I=NULL;
  int ok=true;
  FILE *f;
  fpos_t size;
  char *buffer,*p;

  f=fopen(fname,"r");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadPDBFile","Unable to open file!");
  else
	 {
		if(DebugState&DebugMolecule)
		  {
			printf(" ObjectMoleculeLoadPDBFile: Loading from %s.\n",fname);
		  }
		
		fseek(f,0,SEEK_END);
		fgetpos(f,&size);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);

		I=ObjectMoleculeReadPDBStr(obj,buffer,frame,discrete);

		mfree(buffer);
	 }

  return(I);
}

/*========================================================================*/
void ObjectMoleculeAppendAtoms(ObjectMolecule *I,AtomInfoType *atInfo,CoordSet *cs)
{
  int a;
  int *ii,*si;
  AtomInfoType *src,*dest;
  int nAtom,nBond;

  if(I->NAtom) {
	 nAtom = I->NAtom+cs->NIndex;
	 VLACheck(I->AtomInfo,AtomInfoType,nAtom);	 
	 dest = I->AtomInfo+I->NAtom;
	 src = atInfo;
	 for(a=0;a<cs->NIndex;a++)
		*(dest++)=*(src++);
	 I->NAtom=nAtom;
	 VLAFreeP(atInfo);
  } else {
	 if(I->AtomInfo)
		VLAFreeP(I->AtomInfo);
	 I->AtomInfo = atInfo;
	 I->NAtom=cs->NIndex;
  }
  nBond=I->NBond+cs->NTmpBond;
  if(!I->Bond)
	 I->Bond=VLAlloc(int,nBond*3);
  VLACheck(I->Bond,int,nBond*3);
  ii=I->Bond+I->NBond*3;
  si=cs->TmpBond;
  for(a=0;a<cs->NTmpBond;a++)
	 {
		*(ii++)=cs->IdxToAtm[*(si++)];
		*(ii++)=cs->IdxToAtm[*(si++)];
      *(ii++)=*(si++);
	 }
  I->NBond=nBond;
}
/*========================================================================*/
CoordSet *ObjectMoleculeGetCoordSet(ObjectMolecule *I,int setIndex)
{
  if((setIndex>=0)&&(setIndex<I->NCSet))
	 return(I->CSet[setIndex]);
  else
	 return(NULL);
}
/*========================================================================*/
void ObjectMoleculeTransformTTTf(ObjectMolecule *I,float *ttt,int frame) 
{
  int b;
  CoordSet *cs;
  for(b=0;b<I->NCSet;b++)
	{
     if((frame<0)||(frame==b)) {
       cs=I->CSet[b];
       if(cs) {
         if(cs->fInvalidateRep)
           cs->fInvalidateRep(I->CSet[b],cRepAll,cRepInvCoord);
         MatrixApplyTTTfn3f(cs->NIndex,cs->Coord,ttt,cs->Coord);
       }
     }
	}
}
/*========================================================================*/
void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op)
{
  int a,b,c,s,d;
  int a1,ind;
  float r,rms;
  float v1[3],v2,*vv1,*vv2;
  int inv_flag;
  int hit_flag = false;
  int ok = true;
  OrthoLineType buffer;
  int cnt,maxCnt;
  CoordSet *cs;
  AtomInfoType *ai;

  if(sele>=0) {
	SelectorUpdateTable();
	switch(op->code) {
	case OMOP_AddHydrogens:
     ObjectMoleculeAddSeleHydrogens(I,sele);
     break;
#ifdef _OLD_CODE
     if(!ObjectMoleculeVerifyChemistry(I)) {
       ErrMessage(" AddHydrogens","missing chemical geometry information.");
     } else {
       doneFlag=false;
       while(!doneFlag) {
         doneFlag=true;
         a=0;
         while(a<I->NAtom) {
           ai=I->AtomInfo + a;
           s=I->AtomInfo[a].selEntry;
           if(SelectorIsMember(s,sele))
             if(ObjectMoleculeFillOpenValences(I,a)) {
               hit_flag=true;
               doneFlag=false;
             }
           a++; /* realize that the atom list may have been resorted */
         }
       }
       if(hit_flag)
         ObjectMoleculeSort(I);
     }
     break;
#endif

	case OMOP_PDB1:
	  for(b=0;b<I->NCSet;b++)
       if(I->CSet[b])
		  {
			if((b==op->i1)||(op->i1<0))
			  for(a=0;a<I->NAtom;a++)
				{
				  s=I->AtomInfo[a].selEntry;
              if(SelectorIsMember(s,sele))
                {
                  if(I->DiscreteFlag) {
                    if(I->CSet[b]==I->DiscreteCSet[a])
                      ind=I->DiscreteAtmToIdx[a];
                    else
                      ind=-1;
                  } else 
                    ind=I->CSet[b]->AtmToIdx[a];
                  if(ind>=0) 
                    CoordSetAtomToPDBStrVLA(&op->charVLA,&op->i2,I->AtomInfo+a,
                                            I->CSet[b]->Coord+(3*ind),op->i3);
                  op->i3++;
                }
				}
		  }
	  break;
	case OMOP_AVRT: /* average vertex coordinate */
     cnt=op->nvv1;
     maxCnt=cnt;
	  for(b=0;b<I->NCSet;b++) {
       if(I->CSet[b])
         {
           op->nvv1=cnt;
           cnt=0;
           for(a=0;a<I->NAtom;a++)
             {
				   s=I->AtomInfo[a].selEntry;
               if(SelectorIsMember(s,sele))
                 {
                   if(I->DiscreteFlag) {
                     if(I->CSet[b]==I->DiscreteCSet[a])
                       a1=I->DiscreteAtmToIdx[a];
                     else
                       a1=-1;
                   } else 
                     a1=I->CSet[b]->AtmToIdx[a];
                   if(a1>=0) {
                     if(!cnt) op->i1++;
                     VLACheck(op->vv1,float,(op->nvv1*3)+2);
                     VLACheck(op->vc1,int,op->nvv1);
                     vv2=I->CSet[b]->Coord+(3*a1);
                     vv1=op->vv1+(op->nvv1*3);
                     *(vv1++)+=*(vv2++);
                     *(vv1++)+=*(vv2++);
                     *(vv1++)+=*(vv2++);
                     op->vc1[op->nvv1]++;
                     op->nvv1++;
                   }
                 }
             }
           if(maxCnt<op->nvv1) maxCnt=op->nvv1;
         }
     }
     op->nvv1=maxCnt;
     break;
	case OMOP_SFIT: /* state fitting within a single object */
	  for(b=0;b<I->NCSet;b++) {
       rms = -1.0;
       if(I->CSet[b]&&(b!=op->i1))
         {
           op->nvv1=0;
           for(a=0;a<I->NAtom;a++)
             {
				   s=I->AtomInfo[a].selEntry;
               if(SelectorIsMember(s,sele))
                 {
                   if(I->DiscreteFlag) {
                     if(I->CSet[b]==I->DiscreteCSet[a])
                       a1=I->DiscreteAtmToIdx[a];
                     else
                       a1=-1;
                   } else 
                     a1=I->CSet[b]->AtmToIdx[a];
                   if(a1>=0) {
                     VLACheck(op->vv1,float,(op->nvv1*3)+2);
                     vv2=I->CSet[b]->Coord+(3*a1);
                     vv1=op->vv1+(op->nvv1*3);
                     *(vv1++)=*(vv2++);
                     *(vv1++)=*(vv2++);
                     *(vv1++)=*(vv2++);
                     op->nvv1++;
                   }
                 }
             }
           if(op->nvv1!=op->nvv2) {
             sprintf(buffer,"Atom counts between selections don't match (%d vs %d)\n",
                     op->nvv1,op->nvv2);
             ErrMessage("ExecutiveFit",buffer);
             
           } else if(op->nvv1) {
             if(op->i1!=0)
               rms = MatrixFitRMS(op->nvv1,op->vv1,op->vv2,NULL,op->ttt);
             else 
               rms = MatrixGetRMS(op->nvv1,op->vv1,op->vv2,NULL);
             if(op->i1==2) 
               ObjectMoleculeTransformTTTf(I,op->ttt,b);
           } else {
             ErrMessage("ExecutiveFit","No atoms selected.");
           }
         }
       VLACheck(op->f1VLA,float,b);
       op->f1VLA[b]=rms;
     }
     VLASetSize(op->f1VLA,I->NCSet);  /* NOTE this action is object-specific! */
     break;
	case OMOP_Identify: /* identify atoms */
     for(a=0;a<I->NAtom;a++)
       {
         s=I->AtomInfo[a].selEntry;
         if(SelectorIsMember(s,sele))
           {
             VLACheck(op->i1VLA,int,op->i1);
             op->i1VLA[op->i1++]=I->AtomInfo[a].id;
           }
       }
     break;
	case OMOP_Protect: /* protect atoms from movement */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai->protected = op->i1;
             op->i2++;
           }
         ai++;
       }
     break;
	case OMOP_Mask: /* protect atoms from selection */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai->masked = op->i1;
             op->i2++;
           }
         ai++;
       }
     break;
	case OMOP_SetB: /* set B-value */
     ai=I->AtomInfo;
     for(a=0;a<I->NAtom;a++)
       {
         s=ai->selEntry;
         if(SelectorIsMember(s,sele))
           {
             ai->b = op->f1;
             op->i2++;
           }
         ai++;
       }
     break;
	case OMOP_Remove: /* flag atoms for deletion */
     ai=I->AtomInfo;
     if(I->DiscreteFlag) /* for now, can't remove atoms from discrete objects */
       ErrMessage("Remove","Can't remove atoms from discrete objects.");
     else
       for(a=0;a<I->NAtom;a++)
         {         
           ai->deleteFlag=false;
           s=ai->selEntry;
           if(SelectorIsMember(s,sele))
             {
               ai->deleteFlag=true;
               op->i1++;
             }
           ai++;
         }
     break;
	default:
	   for(a=0;a<I->NAtom;a++)
		 {
		   switch(op->code) { 
         case OMOP_Flag: 
           I->AtomInfo[a].flags &= op->i2; /* clear flag using mask */
           /* no break here - intentional  */
		   case OMOP_COLR: /* normal atom based loops */
		   case OMOP_VISI:
		   case OMOP_TTTF:
         case OMOP_ALTR:
         case OMOP_LABL:
         case OMOP_AlterState:
			 s=I->AtomInfo[a].selEntry;
          if(SelectorIsMember(s,sele))
            {
              switch(op->code) {
              case OMOP_Flag:
                I->AtomInfo[a].flags |= op->i1; /* set flag */
                op->i3++;
                break;
              case OMOP_VISI:
                if(op->i1<0)
                  for(d=0;d<cRepCnt;d++) 
                    I->AtomInfo[a].visRep[d]=op->i2;                      
                else
                  I->AtomInfo[a].visRep[op->i1]=op->i2;
                break;
              case OMOP_COLR:
                I->AtomInfo[a].color=op->i1;
                break;
              case OMOP_TTTF:
                hit_flag=true;
                break;
              case OMOP_LABL:
                if (ok) {
                  if(PLabelAtom(&I->AtomInfo[a],op->s1)) {
                    op->i1++;
                    I->AtomInfo[a].visRep[cRepLabel]=true;
                    hit_flag=true;
                  } else
                    ok=false;
                }
                break;
              case OMOP_ALTR:
                if (ok) {
                  if(PAlterAtom(&I->AtomInfo[a],op->s1,op->i2))
                    op->i1++;
                  else
                    ok=false;
                }
                break;
              case OMOP_AlterState:
                if (ok) {
                  if(op->i2<I->NCSet) {
                    cs = I->CSet[op->i2];
                    if(cs) {
                      if(I->DiscreteFlag) {
                        if(cs==I->DiscreteCSet[a])
                          a1=I->DiscreteAtmToIdx[a];
                        else
                          a1=-1;
                      } else 
                        a1=cs->AtmToIdx[a];
                      if(a1>=0) {
                        if(PAlterAtomState(cs->Coord+(a1*3),op->s1,op->i3)) {
                          op->i1++;
                          hit_flag=true;
                        } else
                          ok=false;
                      }
                    }
                  }
                }
                break;
              }
              break;
            }
			 break;
           
#ifdef PYMOL_FUTURE_CODE
         case OMOP_CSOC: /* specific coordinate set based operations */
           if(I->NCSet<op->cs1) 
             if(I->CSet[op->cs1]) {
               
               s=I->AtomInfo[a].selEntry;
               if(SelectorIsMember(s,sele))
                 {
                   switch(op->code) {
                   case OMOP_CSOC: /* object and coordinate index */
                     break;
                   }
                 }
             }
			 break;
#endif
		   default: /* coord-set based properties, iterating as all coordsets within atoms */
			 for(b=0;b<I->NCSet;b++)
			   if(I->CSet[b])
				 {
				   inv_flag=false;
				   s=I->AtomInfo[a].selEntry;
					   if(SelectorIsMember(s,sele))
						 {
						   switch(op->code) {
						   case OMOP_SUMC:
                       if(I->DiscreteFlag) {
                         if(I->CSet[b]==I->DiscreteCSet[a])
                           a1=I->DiscreteAtmToIdx[a];
                         else
                           a1=-1;
                       } else 
                         a1=I->CSet[b]->AtmToIdx[a];
							 if(a1>=0)
							   {
								 add3f(op->v1,I->CSet[b]->Coord+(3*a1),op->v1);
								 op->i1++;
							   }
							 break;
						   case OMOP_MNMX:
                       if(I->DiscreteFlag) {
                         if(I->CSet[b]==I->DiscreteCSet[a])
                           a1=I->DiscreteAtmToIdx[a];
                         else
                           a1=-1;
                       } else 
                         a1=I->CSet[b]->AtmToIdx[a];
							 if(a1>=0)
							   {
                          if(op->i1) {
                            for(c=0;c<3;c++) {
                              if(*(op->v1+c)>*(I->CSet[b]->Coord+(3*a1+c)))
                                *(op->v1+c)=*(I->CSet[b]->Coord+(3*a1+c));
                              if(*(op->v2+c)<*(I->CSet[b]->Coord+(3*a1+c)))
                                *(op->v2+c)=*(I->CSet[b]->Coord+(3*a1+c));
                            }
                          } else {
                            for(c=0;c<3;c++) {
                              *(op->v1+c)=*(I->CSet[b]->Coord+(3*a1+c));
                              *(op->v2+c)=*(I->CSet[b]->Coord+(3*a1+c));
                              op->i1=1;
                            }
                          }
							   }
							 break;
						   case OMOP_MDST: 
                       if(I->DiscreteFlag) {
                         if(I->CSet[b]==I->DiscreteCSet[a])
                           a1=I->DiscreteAtmToIdx[a];
                         else
                           a1=-1;
                       } else 
                         a1=I->CSet[b]->AtmToIdx[a];
							 if(a1>=0)
							   {
								 r=diff3f(op->v1,I->CSet[b]->Coord+(3*a1));
								 if(r>op->f1)
								   op->f1=r;
							   }
							 break;
						   case OMOP_INVA:
                       if(I->DiscreteFlag) {
                         if(I->CSet[b]==I->DiscreteCSet[a])
                           a1=I->DiscreteAtmToIdx[a];
                         else
                           a1=-1;
                       } else 
                         a1=I->CSet[b]->AtmToIdx[a]; 
                       if(a1>=0)                     /* selection touches this coordinate set */ 
                         inv_flag=true;              /* so set the invalidation flag */
                       break;
						   case OMOP_VERT: 
                       if(I->DiscreteFlag) {
                         if(I->CSet[b]==I->DiscreteCSet[a])
                           a1=I->DiscreteAtmToIdx[a];
                         else
                           a1=-1;
                       } else 
                         a1=I->CSet[b]->AtmToIdx[a];
                       if(a1>=0) {
                         VLACheck(op->vv1,float,(op->nvv1*3)+2);
                         vv2=I->CSet[b]->Coord+(3*a1);
                         vv1=op->vv1+(op->nvv1*3);
                         *(vv1++)=*(vv2++);
                         *(vv1++)=*(vv2++);
                         *(vv1++)=*(vv2++);
                         op->nvv1++;
                       }
							 break;	
						   case OMOP_SVRT:  /* gives us only vertices for a specific coordinate set */
                       if(b==op->i1) {
                         if(I->DiscreteFlag) {
                           if(I->CSet[b]==I->DiscreteCSet[a])
                             a1=I->DiscreteAtmToIdx[a];
                           else
                             a1=-1;
                         } else 
                           a1=I->CSet[b]->AtmToIdx[a];
                         if(a1>=0) {
                           VLACheck(op->vv1,float,(op->nvv1*3)+2);
                           vv2=I->CSet[b]->Coord+(3*a1);
                           vv1=op->vv1+(op->nvv1*3);
                           *(vv1++)=*(vv2++);
                           *(vv1++)=*(vv2++);
                           *(vv1++)=*(vv2++);
                           op->nvv1++;
                         }
                       }
							 break;	
 /* Moment of inertia tensor - unweighted - assumes v1 is center of molecule */
						   case OMOP_MOME: 
                       if(I->DiscreteFlag) {
                         if(I->CSet[b]==I->DiscreteCSet[a])
                           a1=I->DiscreteAtmToIdx[a];
                         else
                           a1=-1;
                       } else 
                         a1=I->CSet[b]->AtmToIdx[a];
							 if(a1>=0) {
							   subtract3f(I->CSet[b]->Coord+(3*a1),op->v1,v1);
							   v2=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]; 
							   op->d[0][0] += v2 - v1[0] * v1[0];
							   op->d[0][1] +=    - v1[0] * v1[1];
							   op->d[0][2] +=    - v1[0] * v1[2];
							   op->d[1][0] +=    - v1[1] * v1[0];
							   op->d[1][1] += v2 - v1[1] * v1[1];
							   op->d[1][2] +=    - v1[1] * v1[2];
							   op->d[2][0] +=    - v1[2] * v1[0];
							   op->d[2][1] +=    - v1[2] * v1[1];
							   op->d[2][2] += v2 - v1[2] * v1[2];
							 }
							 break;
						   }
						 }
				   switch(op->code) { /* full coord-set based */
				   case OMOP_INVA:
                 if(inv_flag) {
                   if(op->i1<0) /* invalidate all representations */
                    for(d=0;d<cRepCnt;d++) {
                      if(I->CSet[b]->fInvalidateRep)
                        I->CSet[b]->fInvalidateRep(I->CSet[b],d,op->i2);
                    }
                  else if(I->CSet[b]->fInvalidateRep) 
                    /* invalidate only that particular representation */
                    I->CSet[b]->fInvalidateRep(I->CSet[b],op->i1,op->i2);
                 }
					 break;
				   }
				 }
			 break;
		   }
		 }
	   break;
	}
	if(hit_flag) {
	  switch(op->code) {
	  case OMOP_TTTF:
       ObjectMoleculeTransformTTTf(I,op->ttt,-1);
       break;
	  case OMOP_LABL:
       ObjectMoleculeInvalidate(I,cRepLabel,cRepInvText);
       break;
     case OMOP_AlterState: /* overly coarse - doing all states, could do just 1 */
       ObjectMoleculeInvalidate(I,-1,cRepInvRep);
       SceneChanged();
       break;
	  }
	}
  }
}
/*========================================================================*/
void ObjectMoleculeDescribeElement(ObjectMolecule *I,int index) 
{
  char buffer[1024];
  AtomInfoType *ai;

  ai=I->AtomInfo+index;
  sprintf(buffer," Atom: [#%d:%s:%s:%d] %s:%s:%s:%s:%s:%s ",
			 ai->id,ai->elem,
          ai->textType,ai->customType,
          I->Obj.Name,ai->segi,ai->chain,
			 ai->resi,ai->resn,ai->name);
  OrthoAddOutput(buffer);
  OrthoNewLine(NULL);
  OrthoRestorePrompt();
}
/*========================================================================*/
int ObjectMoleculeGetNFrames(ObjectMolecule *I)
{
  return I->NCSet;
}
/*========================================================================*/
void ObjectMoleculeUpdate(ObjectMolecule *I)
{
  int a;
  OrthoBusyPrime();
  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a]) {	
	   OrthoBusySlow(a,I->NCSet);
      /*	   printf(" ObjectMolecule: updating state %d of \"%s\".\n" , a+1, I->Obj.Name);*/
      if(I->CSet[a]->fUpdate)
        I->CSet[a]->fUpdate(I->CSet[a]);
	 }
}
/*========================================================================*/
void ObjectMoleculeInvalidate(ObjectMolecule *I,int rep,int level)
{
  int a;
  if(level>=cRepInvBonds) {
    VLAFreeP(I->Neighbor);
    ObjectMoleculeUpdateNonbonded(I);
    if(level>=cRepInvAtoms) {
      SelectorUpdateObjectSele(I);
    }
  }
  for(a=0;a<I->NCSet;a++) 
	 if(I->CSet[a]) {	 
      if(I->CSet[a]->fInvalidateRep)
        I->CSet[a]->fInvalidateRep(I->CSet[a],rep,level);
	 }
}
/*========================================================================*/
int ObjectMoleculeMoveAtom(ObjectMolecule *I,int state,int index,float *v,int mode)
{
  int result = 0;
  CoordSet *cs;
  if(!I->AtomInfo[index].protected) {
    if(state<0) state=0;
    if(I->NCSet==1) state=0;
    state = state % I->NCSet;
    cs = I->CSet[state];
    if(cs) {
      result = CoordSetMoveAtom(I->CSet[state],index,v,mode);
      cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
    }
  }
  return(result);
}
/*========================================================================*/
int ObjectMoleculeGetAtomVertex(ObjectMolecule *I,int state,int index,float *v)
{
  int result = 0;
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  if(I->CSet[state]) 
    result = CoordSetGetAtomVertex(I->CSet[state],index,v);
  return(result);
}
/*========================================================================*/
int ObjectMoleculeSetAtomVertex(ObjectMolecule *I,int state,int index,float *v)
{
  int result = 0;
  if(state<0) state=0;
  if(I->NCSet==1) state=0;
  state = state % I->NCSet;
  if(I->CSet[state]) 
    result = CoordSetSetAtomVertex(I->CSet[state],index,v);
  return(result);
}
/*========================================================================*/
void ObjectMoleculeRender(ObjectMolecule *I,int state,CRay *ray,Pickable **pick)
{
  int a;
  if(state<0) {
    for(a=0;a<I->NCSet;a++)
      if(I->CSet[a])
        if(I->CSet[a]->fRender)
          I->CSet[a]->fRender(I->CSet[a],ray,pick);        
  } else if(state<I->NCSet) {
	 I->CurCSet=state % I->NCSet;
	 if(I->CSet[I->CurCSet]) {
      if(I->CSet[I->CurCSet]->fRender)
        I->CSet[I->CurCSet]->fRender(I->CSet[I->CurCSet],ray,pick);
	 }
  } else if(I->NCSet==1) { /* if only one coordinate set, assume static */
    if(I->CSet[0]->fRender)
      I->CSet[0]->fRender(I->CSet[0],ray,pick);    
  }
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeNew(int discreteFlag)
{
  int a;
  OOAlloc(ObjectMolecule);
  ObjectInit((Object*)I);
  I->Obj.type=cObjectMolecule;
  I->NAtom=0;
  I->NBond=0;
  I->CSet=VLAMalloc(10,sizeof(CoordSet*),5,true); /* auto-zero */
  I->NCSet=0;
  I->Bond=NULL;
  I->DiscreteFlag=discreteFlag;
  if(I->DiscreteFlag) { /* discrete objects don't share atoms between states */
    I->DiscreteAtmToIdx = VLAMalloc(10,sizeof(int),6,false);
    I->DiscreteCSet = VLAMalloc(10,sizeof(CoordSet*),5,false);
    I->NDiscrete=0;
  } else {
    I->DiscreteAtmToIdx = NULL;
    I->DiscreteCSet = NULL;
  }    
  I->Obj.fRender=(void (*)(struct Object *, int, CRay *, Pickable **))ObjectMoleculeRender;
  I->Obj.fFree= (void (*)(struct Object *))ObjectMoleculeFree;
  I->Obj.fUpdate=  (void (*)(struct Object *)) ObjectMoleculeUpdate;
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectMoleculeGetNFrames;
  I->Obj.fDescribeElement = (void (*)(struct Object *,int index)) ObjectMoleculeDescribeElement;
  I->AtomInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
  I->CurCSet=0;
  I->Symmetry=NULL;
  I->Neighbor=NULL;
  for(a=0;a<=cUndoMask;a++) {
    I->UndoCoord[a]=NULL;
    I->UndoState[a]=-1;
  }
  I->UndoIter=0;
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeCopy(ObjectMolecule *obj)
{
  int a;
  int *i0,*i1;
  AtomInfoType *a0,*a1;
  OOAlloc(ObjectMolecule);
  (*I)=(*obj);
  I->Symmetry=NULL; /* TODO: add  copy */

  I->CSet=VLAMalloc(I->NCSet,sizeof(CoordSet*),5,true); /* auto-zero */
  for(a=0;a<I->NCSet;a++) {
    I->CSet[a]=CoordSetCopy(obj->CSet[a]);
    I->CSet[a]->Obj=I;
  }
  I->Bond=VLAlloc(int,I->NBond*3);
  i0=I->Bond;
  i1=obj->Bond;
  for(a=0;a<I->NBond;a++) {
    *(i0++)=*(i1++);
    *(i0++)=*(i1++);
    *(i0++)=*(i1++);    
  }
  
  I->AtomInfo=VLAlloc(AtomInfoType,I->NAtom);
  a0=I->AtomInfo;
  a1=obj->AtomInfo;
  for(a=0;a<I->NAtom;a++)
    *(a0++)=*(a1++);

  for(a=0;a<I->NAtom;a++) {
    I->AtomInfo[a].selEntry=0;
  }
  
  return(I);

}

/*========================================================================*/
void ObjectMoleculeFree(ObjectMolecule *I)
{
  int a;
  SceneObjectDel((Object*)I);
  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a]) {
      if(I->CSet[a]->fFree)
        I->CSet[a]->fFree(I->CSet[a]);
		I->CSet[a]=NULL;
	 }
  if(I->Symmetry) SymmetryFree(I->Symmetry);
  VLAFreeP(I->Neighbor);
  VLAFreeP(I->DiscreteAtmToIdx);
  VLAFreeP(I->DiscreteCSet);
  VLAFreeP(I->CSet);
  VLAFreeP(I->AtomInfo);
  VLAFreeP(I->Bond);
  for(a=0;a<=cUndoMask;a++)
    FreeP(I->UndoCoord[a]);
  OOFreeP(I);
}

/*========================================================================*/
int ObjectMoleculeConnect(ObjectMolecule *I,int **bond,AtomInfoType *ai,
                          CoordSet *cs,float cutoff,int bondSearchFlag)
{
  #define cMULT 1

  int a,b,c,d,e,f,i,j;
  int a1,a2;
  float *v1,*v2,dst;
  int maxBond;
  MapType *map;
  int nBond,*ii1,*ii2;
  int flag;

  nBond = 0;
  maxBond = cs->NIndex * 8;
  (*bond) = VLAlloc(int,maxBond*3);
  if(cs->NIndex&&bondSearchFlag&&(!I->DiscreteFlag))
	 {
      map=MapNew(cutoff+MAX_VDW,cs->Coord,cs->NIndex,NULL);
      if(map)
        {
          for(i=0;i<cs->NIndex;i++)
            {
              v1=cs->Coord+(3*i);
              MapLocus(map,v1,&a,&b,&c);
              for(d=a-1;d<=a+1;d++)
                for(e=b-1;e<=b+1;e++)
                  for(f=c-1;f<=c+1;f++)
                    {
                      j = *(MapFirst(map,d,e,f));
                      while(j>=0)
                        {
                          if(i<j)
                            {
                              v2 = cs->Coord + (3*j);
                              dst = diff3f(v1,v2);										
                              
                              a1=cs->IdxToAtm[i];
                              a2=cs->IdxToAtm[j];
                              
                              dst -= ((ai[a1].vdw+ai[a2].vdw)/2);
                              
                              if( (dst <= cutoff)&&
                                  (!(ai[a1].hydrogen&&ai[a2].hydrogen))&&
                                  ((!cs->TmpBond)||(!(ai[a1].hetatm&&ai[a2].hetatm))))
                                {
                                  flag=true;
                                  if(ai[a1].alt[0]!=ai[a2].alt[0]) { /* handle alternate conformers */
                                    if(ai[a1].alt[0]&&ai[a2].alt)
                                      if(AtomInfoAltMatch(ai+a1,ai+a2))
                                        flag=false;
                                  }
                                  if(flag) {
                                    ai[a1].bonded=true;
                                    ai[a2].bonded=true;
                                    VLACheck((*bond),int,nBond*3+2);
                                    (*bond)[nBond*3  ] = a1;
                                    (*bond)[nBond*3+1] = a2;
                                    (*bond)[nBond*3+2] = 1;
                                    nBond++;
                                  }
                                }
                            }
                          j=MapNext(map,j);
                        }
                    }
            }
          MapFree(map);
        }
      if(DebugState&DebugMolecule) 
        printf("ObjectMoleculeConnect: Found %d bonds.\n",nBond);
    }

  if(cs->NTmpBond&&cs->TmpBond) {
    if(DebugState&DebugMolecule) 
      printf("ObjectMoleculeConnect: incorporating explicit bonds. %d %d\n",
             nBond,cs->NTmpBond);
    VLACheck((*bond),int,(nBond+cs->NTmpBond)*3);
    ii1=(*bond)+nBond*3;
    ii2=cs->TmpBond;
    for(a=0;a<cs->NTmpBond;a++)
      {
        a1 = cs->IdxToAtm[*(ii2++)]; /* convert bonds from index space */
        a2 = cs->IdxToAtm[*(ii2++)]; /* to atom space */
        ai[a1].bonded=true;
        ai[a2].bonded=true;
        *(ii1++)=a1;
        *(ii1++)=a2;
        *(ii1++)=*(ii2++);
      }
    nBond=nBond+cs->NTmpBond;
    VLAFreeP(cs->TmpBond);
    cs->NTmpBond=0;
  }

  if(cs->NTmpLinkBond&&cs->TmpLinkBond) {
    if(DebugState&DebugMolecule) 
      printf("ObjectMoleculeConnect: incorporating linkage bonds. %d %d\n",
             nBond,cs->NTmpLinkBond);
    VLACheck((*bond),int,(nBond+cs->NTmpLinkBond)*3);
    ii1=(*bond)+nBond*3;
    ii2=cs->TmpLinkBond;
    for(a=0;a<cs->NTmpLinkBond;a++)
      {
        a1 = *(ii2++); /* first atom is in object */
        a2 = cs->IdxToAtm[*(ii2++)]; /* second is in the cset */
        ai[a1].bonded=true;
        ai[a2].bonded=true;
        *(ii1++)=a1;
        *(ii1++)=a2;
        *(ii1++)=*(ii2++);
      }
    nBond=nBond+cs->NTmpLinkBond;
    VLAFreeP(cs->TmpLinkBond);
    cs->NTmpLinkBond=0;
  }

  if(!I->DiscreteFlag) {
    UtilSortInPlace((*bond),nBond,sizeof(int)*3,(UtilOrderFn*)BondInOrder);
    if(nBond) { /* eliminate duplicates */
      ii1=(*bond)+3;
      ii2=(*bond)+3;
      a=nBond-1;
      nBond=1;
      if(a>0) 
        while(a--) {
          if((ii2[0]!=ii1[-3])||
             (ii2[1]!=ii1[-2])) {
            *(ii1++)=*(ii2++);
            *(ii1++)=*(ii2++);
            *(ii1++)=*(ii2++);
            nBond++;
          } else {
            ii2+=3;
          }
        }
      VLASize((*bond),int,nBond*3);
    }
  }
  return(nBond);
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadMMDStr(ObjectMolecule *I,char *MMDStr,int frame,int discrete)
{
  int ok = true;
  CoordSet *cset=NULL;
  AtomInfoType *atInfo;
  int isNew;
  int nAtom;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
  cset=ObjectMoleculeMMDStr2CoordSet(MMDStr,&atInfo);  

  if(!cset) 
	 {
		VLAFreeP(atInfo);
		ok=false;
	 }
  
  if(ok)
	 {
		if(!I) 
		  I=(ObjectMolecule*)ObjectMoleculeNew(discrete);
		if(frame<0)
		  frame=I->NCSet;
		if(I->NCSet<=frame)
		  I->NCSet=frame+1;
		VLACheck(I->CSet,CoordSet*,frame);
      nAtom=cset->NIndex;
      cset->Obj = I;
      if(cset->fEnumIndices)
        cset->fEnumIndices(cset);
      if(cset->fInvalidateRep)
        cset->fInvalidateRep(cset,cRepAll,cRepInvRep);
      if(isNew) {		
        I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
        I->NAtom=nAtom;
      } else {
        ObjectMoleculeMerge(I,atInfo,cset,true); /* NOTE: will release atInfo */
      }
      if(frame<0) frame=I->NCSet;
      VLACheck(I->CSet,CoordSet*,frame);
      if(I->NCSet<=frame) I->NCSet=frame+1;
      I->CSet[frame] = cset;
      if(isNew) I->NBond = ObjectMoleculeConnect(I,&I->Bond,I->AtomInfo,cset,0.2,true);
      SceneCountFrames();
      ObjectMoleculeExtendIndices(I);
      ObjectMoleculeSort(I);
      ObjectMoleculeUpdateNonbonded(I);
	 }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadMMDFile(ObjectMolecule *obj,char *fname,
                                          int frame,char *sepPrefix,int discrete)
{
  ObjectMolecule* I=NULL;
  int ok=true;
  FILE *f;
  int oCnt=0;
  fpos_t size;
  char *buffer,*p;
  char cc[MAXLINELEN],oName[ObjNameMax];
  int nLines;
  f=fopen(fname,"r");
  if(!f)
	 ok=ErrMessage("ObjectMoleculeLoadMMDFile","Unable to open file!");
  else
	 {
		if(DebugState&DebugMolecule)
		  {
			 printf(" ObjectMoleculeLoadMMDFile: Loading from %s.\n",fname);
			 fflush(stdout);
		  }		
		fseek(f,0,SEEK_END);
		fgetpos(f,&size);
		fseek(f,0,SEEK_SET);
		buffer=(char*)mmalloc(size+255);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		p[size]=0;
		fclose(f);
      p=buffer;
      while(ok) {
        ncopy(cc,p,6);
        if(sscanf(cc,"%d",&nLines)!=1)
          break;
        if(ok) {
          if(sepPrefix) {
            I=ObjectMoleculeReadMMDStr(NULL,p,frame,discrete);
            oCnt++;
            sprintf(oName,"%s-%02d",sepPrefix,oCnt);
            ObjectSetName((Object*)I,oName);
            ExecutiveManageObject((Object*)I);
          } else {
            I=ObjectMoleculeReadMMDStr(obj,p,frame,discrete);
            obj=I;
          }
          p=nextline(p);
          while(nLines--)
            p=nextline(p);
        }
      }
		mfree(buffer);
	 }

  return(I);
}

/*========================================================================*/
CoordSet *ObjectMoleculePDBStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom;
  int a,c,llen;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  char cc[MAXLINELEN];
  int AFlag;
  int atomCount;
  int conectFlag = false;
  int *bond=NULL,*ii1,*ii2,*idx;
  int nBond=0;
  int b1,b2,nReal,maxAt;
  CSymmetry *symmetry = NULL;
  int symFlag;
  int autoshow_lines = SettingGet(cSetting_autoshow_lines);
  int autoshow_nonbonded = SettingGet(cSetting_autoshow_nonbonded);

  AtomInfoPrimeColors();

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  if(!atInfo)
    ErrFatal("PDBStr2CoordSet","need atom information record!"); /* failsafe for old version..*/

  while(*p)
	 {
		if((*p == 'A')&&(*(p+1)=='T')&&(*(p+2)=='O')&&(*(p+3)=='M'))
		  nAtom++;
		if((*p == 'H')&&(*(p+1)=='E')&&(*(p+2)=='T')&&
         (*(p+3)=='A')&&(*(p+4)=='T')&&(*(p+5)=='M'))
        nAtom++;
		if((*p == 'C')&&(*(p+1)=='O')&&(*(p+2)=='N')&&
         (*(p+3)=='E')&&(*(p+4)=='C')&&(*(p+5)=='T'))
        conectFlag=true;
      p=nextline(p);
	 }
  
  coord=VLAlloc(float,3*nAtom);

  if(atInfo)
	 VLACheck(atInfo,AtomInfoType,nAtom);

  if(conectFlag) {
    nBond=0;
    bond=VLAlloc(int,12*nAtom);  
  }
  p=buffer;
  if(DebugState & DebugMolecule) {
	 printf(" ObjectMoleculeReadPDB: Found %i atoms...\n",nAtom);
	 fflush(stdout);
  }
  a=0;
  atomCount=0;
  
  while(*p)
	 {
		AFlag=false;
		if((*p == 'A')&&(*(p+1)=='T')&&(*(p+2)=='O')&&(*(p+3)=='M'))
		  AFlag = 1;
		if((*p == 'H')&&(*(p+1)=='E')&&(*(p+2)=='T')&&
         (*(p+3)=='A')&&(*(p+4)=='T')&&(*(p+5)=='M'))
        AFlag = 2;
		if((*p == 'R')&&(*(p+1)=='E')&&(*(p+2)=='M')&&
         (*(p+3)=='A')&&(*(p+4)=='R')&&(*(p+5)=='K')&&
         (*(p+6)==' ')&&(*(p+7)=='2')&&(*(p+8)=='9')&&
         (*(p+9)=='0'))
        {
        }
		if((*p == 'C')&&(*(p+1)=='R')&&(*(p+2)=='Y')&&
         (*(p+3)=='S')&&(*(p+4)=='T')&&(*(p+5)=='1'))
        {
          if(!symmetry) symmetry=SymmetryNew();          
          if(symmetry) {
            ErrOk(" PDBStrToCoordSet","Attempting to read symmetry information");
            p=nskip(p,6);
            symFlag=true;
            p=ncopy(cc,p,9);
            if(sscanf(cc,"%f",&symmetry->Crystal->Dim[0])!=1) symFlag=false;
            p=ncopy(cc,p,9);
            if(sscanf(cc,"%f",&symmetry->Crystal->Dim[1])!=1) symFlag=false;
            p=ncopy(cc,p,9);
            if(sscanf(cc,"%f",&symmetry->Crystal->Dim[2])!=1) symFlag=false;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%f",&symmetry->Crystal->Angle[0])!=1) symFlag=false;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%f",&symmetry->Crystal->Angle[1])!=1) symFlag=false;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%f",&symmetry->Crystal->Angle[2])!=1) symFlag=false;
            p=nskip(p,1);
            p=ncopy(symmetry->SpaceGroup,p,10);
            UtilCleanStr(symmetry->SpaceGroup);
            p=ncopy(cc,p,4);
            if(sscanf(cc,"%d",&symmetry->PDBZValue)!=1) symmetry->PDBZValue=1;
            if(!symFlag) {
              ErrMessage("PDBStrToCoordSet","Error reading CRYST1 record\n");
              SymmetryFree(symmetry);
              symmetry=NULL;
            } 
          }
        }
		if((*p == 'C')&&(*(p+1)=='O')&&(*(p+2)=='N')&&
         (*(p+3)=='E')&&(*(p+4)=='C')&&(*(p+5)=='T'))
        {
          p=nskip(p,6);
          p=ncopy(cc,p,5);
          if(sscanf(cc,"%d",&b1)==1)
            while (1) {
              p=ncopy(cc,p,5);
              if(sscanf(cc,"%d",&b2)!=1)
                break;
              else {
                VLACheck(bond,int,(nBond*3)+2);
                if(b1<=b2) {
                  bond[nBond*3]=b1; /* temporarily store the atom indexes */
                  bond[nBond*3+1]=b2;
                  bond[nBond*3+2]=1;
                } else {
                  bond[nBond*3]=b2;
                  bond[nBond*3+1]=b1;
                  bond[nBond*3+2]=1;
                }
                nBond++;
              }
            }
        }
		if(AFlag)
		  {
			 llen=0;
			 while((*(p+llen))&&((*(p+llen))!=13)&&((*(p+llen))!=10)) {
				llen++;
			 }

          ai=atInfo+atomCount;

          p=nskip(p,6);
          p=ncopy(cc,p,5);
          if(!sscanf(cc,"%d",&ai->id)) ai->id=0;

          p=nskip(p,1);/* to 12 */
          p=ncopy(cc,p,4); 
          if(!sscanf(cc,"%s",ai->name)) ai->name[0]=0;
          
          p=ncopy(cc,p,1);
          if(*cc==32)
            ai->alt[0]=0;
          else {
            ai->alt[0]=*cc;
            ai->alt[1]=0;
          }

          p=ncopy(cc,p,3); 
          if(!sscanf(cc,"%s",ai->resn)) ai->resn[0]=0;

          p=nskip(p,1);
          p=ncopy(cc,p,1);
          if(*cc==' ')
            ai->chain[0]=0;
          else {
            ai->chain[0] = *cc;
            ai->chain[1] = 0;
          }

          p=ncopy(cc,p,5); /* we treat insertion records as part of the residue identifier */
          if(!sscanf(cc,"%s",ai->resi)) ai->resi[0]=0;
          if(!sscanf(cc,"%d",&ai->resv)) ai->resv=1;
          
          p=nskip(p,3);
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+a);
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+(a+1));
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+(a+2));

          p=ncopy(cc,p,6);
          if(!sscanf(cc,"%f",&ai->q))
            ai->q=1.0;
          
          p=ncopy(cc,p,6);
          if(!sscanf(cc,"%f",&ai->b))
            ai->b=0.0;

          p=nskip(p,6);
          p=ncopy(cc,p,4);
          if(!sscanf(cc,"%s",ai->segi)) ai->segi[0]=0;
          
          p=ncopy(cc,p,2);
          if(!sscanf(cc,"%s",ai->elem)) 
            ai->elem[0]=0;          
          else if(!(((ai->elem[0]>='a')&&(ai->elem[0]<='z'))||
                    ((ai->elem[0]>='A')&&(ai->elem[0]<='Z'))))
            ai->elem[0]=0;                      

          for(c=0;c<cRepCnt;c++) {
            ai->visRep[c] = false;
          }
          ai->visRep[cRepLine] = autoshow_lines; /* show lines by default */
          ai->visRep[cRepNonbonded] = autoshow_nonbonded; /* show lines by default */

          if(AFlag==1) 
            ai->hetatm=0;
          else
            ai->hetatm=1;
          
          AtomInfoAssignParameters(ai);
          ai->color=AtomInfoGetColor(ai);

          if(DebugState&DebugMolecule)
            printf("%s %s %s %s %8.3f %8.3f %8.3f %6.2f %6.2f %s\n",
                    ai->name,ai->resn,ai->resi,ai->chain,
                    *(coord+a),*(coord+a+1),*(coord+a+2),ai->b,ai->q,
                    ai->segi);

			 a+=3;
			 atomCount++;
		  }
      p=nextline(p);
	 }
  if(conectFlag) {
    UtilSortInPlace(bond,nBond,sizeof(int)*3,(UtilOrderFn*)BondInOrder);              
    if(nBond) {
      ii1=bond;
      ii2=bond+3;
      nReal=1;
      ii1[2]=1;
      a=nBond-1;
      while(a) {
        if((ii1[0]==ii2[0])&&(ii1[1]==ii2[1])) {
          ii1[2]++; /* count dup */
        } else {
          ii1+=3; /* non-dup, make copy */
          ii1[0]=ii2[0];
          ii1[1]=ii2[1];
          ii1[2]=ii2[2];
          nReal++;
        }
        ii2+=3;
        a--;
      }
      nBond=nReal;
      /* now, find atoms we're looking for */
      maxAt=nAtom;
      ii1=bond;
      for(a=0;a<nBond;a++) {
        if(ii1[0]>maxAt) maxAt=ii1[0];
        if(ii1[1]>maxAt) maxAt=ii1[1];
        ii1+=3;
      }
      for(a=0;a<nAtom;a++) 
        if(maxAt<atInfo[a].id) maxAt=atInfo[a].id;
      /* build index */
      idx = Alloc(int,maxAt+1);
      for(a=0;a<maxAt;a++) idx[a]=-1;
      for(a=0;a<nAtom;a++)
        idx[atInfo[a].id]=a;
      /* convert indices to bonds */
      ii1=bond;
      ii2=bond;
      nReal=0;
      for(a=0;a<nBond;a++) {
        ii2[0]=idx[ii1[0]];
        ii2[1]=idx[ii1[1]];
        if((ii2[0]>=0)&&(ii2[1]>=0)) {
          if(ii1[2]<=2) ii2[2]=1;
          else if(ii1[2]<=4) ii2[2]=2;
          else ii2[2]=3;
          nReal++;
          ii2+=3;
        }
        ii1+=3;
      }
      nBond=nReal;
    /* first, count and eliminate duplicates */
      FreeP(idx);
    }
  }
  if(DebugState&DebugMolecule)
    printf(" PDBStr2CoordSet: Read %d bonds from CONECT records (%p).\n",nBond,bond);
  cset = CoordSetNew();
  cset->NIndex=nAtom;
  cset->Coord=coord;
  cset->TmpBond=bond;
  cset->NTmpBond=nBond;
  if(symmetry) cset->TmpSymmetry=symmetry;
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  return(cset);
}

/*========================================================================*/
CoordSet *ObjectMoleculeMMDStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom,nBond;
  int a,c,bPart,bOrder;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  char cc[MAXLINELEN];
  float *f;
  int *ii,*bond=NULL;
  int ok=true;
  int autoshow_lines = SettingGet(cSetting_autoshow_lines);
  int autoshow_nonbonded = SettingGet(cSetting_autoshow_nonbonded);
  AtomInfoPrimeColors();

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;


  if(ok) {
	 p=ncopy(cc,p,6);
	 if(sscanf(cc,"%d",&nAtom)!=1)
		ok=ErrMessage("ReadMMDFile","bad atom count");
  }

  if(ok) {
	 coord=VLAlloc(float,3*nAtom);
	 if(atInfo)
		VLACheck(atInfo,AtomInfoType,nAtom);	 
  }

  if(!atInfo)
    ErrFatal("PDBStr2CoordSet","need atom information record!"); /* failsafe for old version..*/

  nBond=0;
  if(ok) {
	 bond=VLAlloc(int,18*nAtom);  
  }
  p=nextline(p);

  /* read coordinates and atom names */

  if(ok) { 
	 f=coord;
	 ii=bond;
	 for(a=0;a<nAtom;a++)
		{
        ai=atInfo+a;

        ai->id=a+1;
        if(ok) {
          p=ncopy(cc,p,4);
          if(sscanf(cc,"%d",&ai->customType)!=1) 
            ok=ErrMessage("ReadMMDFile","bad atom type");
        }
        if(ok) {
          if(ai->customType<=14) strcpy(ai->elem,"C");
          else if(ai->customType<=23) strcpy(ai->elem,"O");
          else if(ai->customType<=40) strcpy(ai->elem,"N");
          else if(ai->customType<=48) strcpy(ai->elem,"H");
          else if(ai->customType<=52) strcpy(ai->elem,"S");
          else if(ai->customType<=53) strcpy(ai->elem,"P");
          else if(ai->customType<=55) strcpy(ai->elem,"B");
          else if(ai->customType<=56) strcpy(ai->elem,"F");
          else if(ai->customType<=57) strcpy(ai->elem,"Cl");           
          else if(ai->customType<=58) strcpy(ai->elem,"Br");           
          else if(ai->customType<=59) strcpy(ai->elem,"I");           
          else if(ai->customType<=60) strcpy(ai->elem,"Si");           
          else if(ai->customType<=61) strcpy(ai->elem,"Du");           
          else if(ai->customType<=62) strcpy(ai->elem,"Z0");
          else if(ai->customType<=63) strcpy(ai->elem,"Lp");
          else strcpy(ai->elem,"?");
        }
        for(c=0;c<6;c++) {
          if(ok) {
            p=ncopy(cc,p,8);
            if(sscanf(cc,"%d%d",&bPart,&bOrder)!=2)
              ok=ErrMessage("ReadMMDFile","bad bond record");
            else {
              if(bPart&&bOrder&&(a<(bPart-1))) {
                nBond++;
                *(ii++)=a;
                *(ii++)=bPart-1;
                *(ii++)=bOrder;
              }
            }
          }
        }
        if(ok) {
          p=ncopy(cc,p,12);
          if(sscanf(cc,"%f",f++)!=1)
            ok=ErrMessage("ReadMMDFile","bad coordinate");
        }
        if(ok) {
          p=ncopy(cc,p,12);
          if(sscanf(cc,"%f",f++)!=1)
            ok=ErrMessage("ReadMMDFile","bad coordinate");
        }
        if(ok) {
          p=ncopy(cc,p,12);
			 if(sscanf(cc,"%f",f++)!=1)
				ok=ErrMessage("ReadMMDFile","bad coordinate");
		  }
        if(ok) {
          p=nskip(p,12);
          p=ncopy(cc,p,9);
			 if(sscanf(cc,"%f",&ai->partialCharge)!=1)
				ok=ErrMessage("ReadMMDFile","bad charge");
        }
        if(ok) {
          p=nskip(p,10);
          p=ncopy(cc,p,3);
          if(sscanf(cc,"%s",ai->resn)!=1)
            ai->resn[0]=0;
          ai->hetatm=true;
        }

        ai->segi[0]=0;
        ai->alt[0]=0;

        if(ok) {
          p=nskip(p,2);
          p=ncopy(ai->name,p,4);
          UtilCleanStr(ai->name);
          if(ai->name[0]==0) {
            strcpy(ai->name,ai->elem);
            sprintf(cc,"%02d",a+1);
            if((strlen(cc)+strlen(ai->name))>4)
              strcpy(ai->name,cc);
            else
              strcat(ai->name,cc);
          }

          for(c=0;c<cRepCnt;c++) {
            ai->visRep[c] = false;
          }
          ai->visRep[cRepLine] = autoshow_lines; /* show lines by default */
          ai->visRep[cRepNonbonded] = autoshow_nonbonded; /* show lines by default */

        }
        if(ok) {
          AtomInfoAssignParameters(ai);
          ai->color = AtomInfoGetColor(ai);
        }
        if(!ok)
          break;
        p=nextline(p);
      }
  }
  if(ok)
    bond=VLASetSize(bond,3*nBond);
  if(ok) {
	 cset = CoordSetNew();
	 cset->NIndex=nAtom;
	 cset->Coord=coord;
	 cset->NTmpBond=nBond;
	 cset->TmpBond=bond;
  } else {
	 VLAFreeP(bond);
	 VLAFreeP(coord);
  }
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  return(cset);
}


#ifdef _NO_LONGER_USED
/*========================================================================*/
static char *nextline(char *p) {
  while(*p) {
	 if(*p==0xD) { /* Mac or PC */
		if(*(p+1)==0xA) /* PC */
		  p++;
		p++;
		break;
	 }
	 if(*p==0xA) /* Unix */
		{
		  p++;
		  break;
		}
	 p++;
  }
  return p;
}
/*========================================================================*/
static char *wcopy(char *q,char *p,int n) { /* word copy */
  while(*p) {
	 if(*p<=32) 
		p++;
	 else
		break;
  }
  while(*p) {
	 if(*p<=32)
		break;
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
	 *(q++)=*(p++);
	 n--;
  }
  *q=0;
  return p;
}
/*========================================================================*/
static char *ncopy(char *q,char *p,int n) {  /* n character copy */
  while(*p) {
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
	 *(q++)=*(p++);
	 n--;
  }
  *q=0;
  return p;
}
/*========================================================================*/
static char *nskip(char *p,int n) {  /* n character skip */
  while(*p) {
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* stop at newlines */
		break;
    p++;
	 n--;
  }
  return p;
}

#endif
