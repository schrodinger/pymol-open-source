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

#include<string.h>
#include<stdio.h>
#include<ctype.h>
#include<math.h>

#include"Base.h"
#include"Map.h"
#include"Vector.h"
#include"Debug.h"
#include"Err.h"
#include"Word.h"
#include"Util.h"
#include"MemoryDebug.h"
#include"Selector.h"
#include"Executive.h"
#include"ObjectMolecule.h"
#include"CoordSet.h"
#include"DistSet.h"
#include"Word.h"
#include"Scene.h"

#define SelectorMaxDepth 100


typedef struct {
  int selection;
  int next;
} MemberType;

typedef struct {
  int level;
  int type; /* 0 = value 1 = operation 2 = pre-operation */
  unsigned int code; 
  WordType text;
  int *sele;
} EvalElem;

typedef struct {
  int model;
  int atom;
  int index;
} TableRec;

typedef struct {
  WordType *Name;
  int NSelection;
  int TmpCounter;
  MemberType *Member;
  int NMember;
  ObjectMolecule **Obj;
  TableRec *Table;
  float *Vertex;
  int *Flag1,*Flag2;
  int NAtom;
  int NModel;
  int NCSet;
  int IgnoreCase;
} SelectorType;

SelectorType Selector;

int SelectorGetInterstateVLA(int sele1,int state1,int sele2,int state2,
									  float cutoff,int **vla);
int SelectorGetArrayNCSet(int *array);

int SelectorModulate1(EvalElem *base);
int SelectorSelect0(EvalElem *base);
int SelectorSelect1(EvalElem *base);
int SelectorSelect2(EvalElem *base);
int SelectorLogic1(EvalElem *base);
int SelectorLogic2(EvalElem *base);
int *SelectorEvaluate(WordType *word);
WordType *SelectorParse(char *s);
void SelectorPurgeMembers(int sele);

#define STYP_VALU 0
#define STYP_OPR1 1
#define STYP_OPR2 2
#define STYP_SEL0 3
#define STYP_SEL1 4
#define STYP_SEL2 5
#define STYP_LIST 6
#define STYP_PRP1 7
#define STYP_SEL3 8
#define STYP_PVAL 0

#define SELE_NOT1 ( 0x0100 | STYP_OPR1 )
#define SELE_BYR1 ( 0x0200 | STYP_OPR1 )
#define SELE_AND2 ( 0x0300 | STYP_OPR2 )
#define SELE_OR_2 ( 0x0400 | STYP_OPR2 )
#define SELE_IN_2 ( 0x0500 | STYP_OPR2 )
#define SELE_ALLz ( 0x0600 | STYP_SEL0 )
#define SELE_NONz ( 0x0700 | STYP_SEL0 )
#define SELE_HETz ( 0x0800 | STYP_SEL0 )
#define SELE_HYDz ( 0x0900 | STYP_SEL0 )
#define SELE_VISz ( 0x0A00 | STYP_SEL0 )
#define SELE_ARD_ ( 0x0B00 | STYP_PRP1 )
#define SELE_EXP_ ( 0x0C00 | STYP_PRP1 )
#define SELE_NAMs ( 0x0D00 | STYP_SEL1 )
#define SELE_ELEs ( 0x0E00 | STYP_SEL1 )
#define SELE_RSIs ( 0x0F00 | STYP_SEL1 )
#define SELE_CHNs ( 0x1000 | STYP_SEL1 )
#define SELE_SEGs ( 0x1100 | STYP_SEL1 )
#define SELE_MODs ( 0x1200 | STYP_SEL1 ) 
#define SELE_IDXs ( 0x1300 | STYP_SEL1 )
#define SELE_RSNs ( 0x1400 | STYP_SEL1 )
#define SELE_SELs ( 0x1500 | STYP_SEL1 )
#define SELE_BVLx ( 0x1606 | STYP_SEL2 )
#define SELE_ALTs ( 0x1700 | STYP_SEL1 )

static WordKeyValue Keyword[] = 
{
  {  "not",      SELE_NOT1 },
  {  "!",        SELE_NOT1 },
  {  "byresi",   SELE_BYR1 },
  {  "b;",       SELE_BYR1 },
  {  "and",      SELE_AND2 },
  {  "&",        SELE_AND2 },
  {  "or",       SELE_OR_2 },
  {  "|",        SELE_OR_2 },
  {  "in",       SELE_IN_2 },
  {  "all",      SELE_ALLz }, /* 0 parameter */
  {  "+",        SELE_ALLz }, /* 0 parameter */
  {  "none",     SELE_NONz }, /* 0 parameter */
  {  "hetatm",   SELE_HETz }, /* 0 parameter */
  {  "het",      SELE_HETz }, /* 0 parameter */
  {  "hydro",    SELE_HYDz }, /* 0 parameter */
  {  "h;",       SELE_HYDz }, /* 0 parameter */
  {  "visi",     SELE_VISz }, /* 0 parameter */
  {  "v;",       SELE_VISz }, /* 0 parameter */
  {  "around",   SELE_ARD_ }, /* 1 parameter */
  {  "a;",       SELE_ARD_ }, /* 1 parameter */
  {  "expand",   SELE_EXP_ }, /* 1 parameter */
  {  "x;",       SELE_EXP_ }, /* 1 parameter */
  {  "name",     SELE_NAMs },
  {  "n;",       SELE_NAMs },
  {  "elem",     SELE_ELEs },
  {  "e;",       SELE_ELEs },
  {  "resi",     SELE_RSIs },
  {  "i;",       SELE_RSIs },
  {  "alt",      SELE_ALTs },
  {  "l;",       SELE_ALTs },
  {  "chain",    SELE_CHNs },
  {  "c;",       SELE_CHNs },
  {  "segi",     SELE_SEGs },
  {  "s;",       SELE_SEGs },
  {  "model",    SELE_MODs },
  {  "m;",       SELE_MODs },
  {  "index",    SELE_IDXs },
  {  "resn",     SELE_RSNs },
  {  "r;",       SELE_RSNs },
  {  "%",        SELE_SELs },
  {  "b",        SELE_BVLx, }, /* 2 operand selection operator */ 
  {  "", 0 }
};

#define SCMP_GTHN 0x01
#define SCMP_LTHN 0x02
#define SCMP_RANG 0x03
#define SCMP_EQAL 0x04

static WordKeyValue AtOper[] = 
{
 { ">",      SCMP_GTHN },
 { "<",      SCMP_LTHN },
 { "in",     SCMP_RANG },
 { "=",      SCMP_EQAL },
 { "", 0 }
};

static int BondInOrder(int *a,int b1,int b2);
static int BondCompare(int *a,int *b);

/*========================================================================*/
static int BondInOrder(int *a,int b1,int b2)
{
  return(BondCompare(a+b1*2,a+b2*2)<=0);
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
int SelectorGetSeleNCSet(int sele)
{
  SelectorType *I=&Selector;
  int a,s,at;
  ObjectMolecule *obj;
  int result=0;
  for(a=0;a<I->NAtom;a++) {
    obj=I->Obj[I->Table[a].model];
    at=I->Table[a].atom;
    s=obj->AtomInfo[at].selEntry;
    while(s) 
      {
        if(I->Member[s].selection==sele)
			 if(result<obj->NCSet) result=obj->NCSet;
        s=SelectorNext(s);
		}
  }
  return(result);
}
/*========================================================================*/
int SelectorGetArrayNCSet(int *array)
{
  SelectorType *I=&Selector;
  int a;
  ObjectMolecule *obj;
  int result=0;

  for(a=0;a<I->NAtom;a++) 
	 if(*(array++)) {
		obj=I->Obj[I->Table[a].model];
		if(result<obj->NCSet) result=obj->NCSet;
	 }
  return(result);
}
/*========================================================================*/
float SelectorSumVDWOverlap(int sele1,int state1,int sele2,int state2)
{
  SelectorType *I=&Selector;
  int *vla=NULL;
  int c;
  float result=0.0;
  float sumVDW,dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at1,at2;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj1,*obj2;
  int idx1,idx2;
  int a;

  SelectorUpdateTable();
  c=SelectorGetInterstateVLA(sele1,state1,sele2,state2,2*MAX_VDW,&vla);
  for(a=0;a<c;a++) {
    a1=vla[a*2];
    a2=vla[a*2+1];

    at1=I->Table[a1].atom;
    at2=I->Table[a2].atom;
    
    obj1=I->Obj[I->Table[a1].model];
    obj2=I->Obj[I->Table[a2].model];

    if((state1<obj1->NCSet)&&(state2<obj2->NCSet)) {
      cs1=obj1->CSet[state1];
      cs2=obj2->CSet[state2];
      if(cs1&&cs2) { /* should always be true */
        
        ai1=obj1->AtomInfo+at1;
        ai2=obj2->AtomInfo+at2;
        
        idx1=cs1->AtmToIdx[at1];
        idx2=cs2->AtmToIdx[at2];
        
        sumVDW=ai1->vdw+ai2->vdw;
        dist=diff3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2);

        if(dist<sumVDW) {
          result+=((sumVDW-dist)/2.0);
        }
      }
    }
  }
  VLAFreeP(vla);
  return(result);
}
/*========================================================================*/
int SelectorGetInterstateVLA(int sele1,int state1,int sele2,int state2,
										float cutoff,int **vla) /* Assumes valid tables */
{
  SelectorType *I=&Selector;
  MapType *map;
  float *v2;
  int n1,n2;
  int c,i,j,h,k,l;
  int at;
  int a,s,idx;
  ObjectMolecule *obj;
  CoordSet *cs;

  if(!(*vla))
	 (*vla)=VLAlloc(int,1000);

  c=0;
  n1=0;

  for(a=0;a<I->NAtom;a++) {
	 I->Flag1[a]=false;
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    while(s) 
      {
        if(I->Member[s].selection==sele1)
          {
            if(state1<obj->NCSet) 
              cs=obj->CSet[state1];
            else
              cs=NULL;
				if(cs) {
				  idx=cs->AtmToIdx[at];
				  if(idx>=0) {
					 copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
					 I->Flag1[a]=true;
					 n1++;
				  }
				}
			 }
        s=SelectorNext(s);
		}
  }
  /* now create and apply voxel map */
  c=0;
  if(n1) {
	 n2=0;
	 map=MapNewFlagged(-cutoff,I->Vertex,I->NAtom,NULL,I->Flag1);
	 if(map) {
		MapSetupExpress(map);
		for(a=0;a<I->NAtom;a++) {
		  at=I->Table[a].atom;
		  obj=I->Obj[I->Table[a].model];
		  s=obj->AtomInfo[at].selEntry;
		  while(s) 
			 {
				if(I->Member[s].selection==sele2)
				  {
                if(state2<obj->NCSet) 
                  cs=obj->CSet[state2];
                else
                  cs=NULL;
					 if(cs) {
						idx=cs->AtmToIdx[at];
						if(idx>=0) {
						  v2 = cs->Coord+(3*idx);
						  MapLocus(map,v2,&h,&k,&l);
						  i=*(MapEStart(map,h,k,l));
						  if(i) {
							 j=map->EList[i++];
							 while(j>=0) {
								if(within3f(I->Vertex+3*j,v2,cutoff)) {
								  VLACheck((*vla),int,c*2+1);
								  *((*vla)+c*2)=j;
								  *((*vla)+c*2+1)=a;
								  c++;
								}
								j=map->EList[i++];
							 }
						  }
						  n2++;
						}
					 }
				  }
				s=SelectorNext(s);
			 }
		}
		MapFree(map);
	 }
  }
  return(c);
}
/*========================================================================*/
int SelectorGetPDB(char **charVLA,int sele,int state,int conectFlag)
{
  SelectorType *I=&Selector;

  int a,b,b1,b2,c,d,*ii1,s,idx,at;
  int *bond=NULL;
  int nBond=0;
  int cLen =0;
  int newline;
  CoordSet *cs;
  ObjectMolecule *obj;
  AtomInfoType *atInfo;

  SelectorUpdateTable();
  c=0;
  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    I->Table[a].index=0;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    while(s) 
      {
        if(I->Member[s].selection==sele)
          {
            if(state<obj->NCSet) 
              cs=obj->CSet[state];
            else
              cs=NULL;
            if(cs) {
              idx=cs->AtmToIdx[at];
              if(idx>=0) {
                I->Table[a].index=c+1; /* NOTE marking with "1" based indexes here */
                CoordSetAtomToPDBStrVLA(charVLA,&cLen,obj->AtomInfo+at,
                                        obj->CSet[state]->Coord+(3*idx),c);
                c++;
              }
              break;
            }
          }
        s=SelectorNext(s);
      }
  }
  if(conectFlag) {
    nBond = 0;
    bond = VLAlloc(int,1000);
    for(a=0;a<I->NModel;a++) {
      obj=I->Obj[a];
      ii1=obj->Bond;
      if(state<obj->NCSet) 
        cs=obj->CSet[state];
      else
        cs=NULL;
      if(cs) {
        atInfo=obj->AtomInfo;
        for(b=0;b<obj->NBond;b++) {
          b1=ii1[0];
          b2=ii1[1];        
          if((cs->AtmToIdx[b1]>=0)&&(cs->AtmToIdx[b2]>=0)&&
             (atInfo[b1].hetatm||atInfo[b2].hetatm)) {
            b1+=obj->SeleBase;
            b2+=obj->SeleBase;
            if(I->Table[b1].index&&I->Table[b2].index) {
              VLACheck(bond,int,nBond*2+12);
              b1=I->Table[b1].index;
              b2=I->Table[b2].index;
              for(d=0;d<ii1[2];d++) {
                bond[nBond*2] = b1;
                bond[nBond*2+1] = b2;
                nBond++;
                bond[nBond*2] = b2;
                bond[nBond*2+1] = b1;
                nBond++;
              }
            }
          }
        ii1+=3;
        }
      }
    }
    UtilSortInPlace(bond,nBond,sizeof(int)*2,(UtilOrderFn*)BondInOrder);
    ii1=bond;
    b1=-1;
	 b2=-1;
    newline = false;
    for(a=0;a<nBond;a++) {
      if(a<(nBond-1)) 
        if((ii1[0]==ii1[2])&&(ii1[1]==ii1[3])) newline=true;
      if(b1!=ii1[0]||((b1==ii1[0])&&(b2==ii1[1]))||newline) {
        if(a) cLen+=sprintf((*charVLA)+cLen,"\n");
        cLen+=sprintf((*charVLA)+cLen,"CONECT%5d%5d",
                      ii1[0],ii1[1]);
        b1=ii1[0];
		  b2=ii1[1];
        newline=false;
      } else cLen+=sprintf((*charVLA)+cLen,"%5d",
                           ii1[1]);
      b2=ii1[1];
      ii1+=2;
    }
    cLen+=sprintf((*charVLA)+cLen,"\n");
    VLAFree(bond);
  }
  return(cLen);
}
/*========================================================================*/
void SelectorCreateObjectMolecule(int sele,char *name,int target,int source)
{
  SelectorType *I=&Selector;

  int a,b,a1,a2,b1,b2,c,d,*ii1,s,at;
  int *bond=NULL;
  int nBond=0;
  int nCSet,nAtom,ts;
  AtomInfoType *atInfo = NULL;
  int isNew,csFlag;
  CoordSet *cs = NULL;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj;
  Object *ob;
  ObjectMolecule *targ = NULL;

  ob=ExecutiveFindObjectByName(name);
  if(ob)
    if(ob->type==cObjectMolecule) 
      targ = (ObjectMolecule*)ob;
  if(!targ) {
    isNew=true;
    targ = ObjectMoleculeNew();
    targ->Bond = VLAlloc(int,1);
  } else {
    isNew=false;
  }

  SelectorUpdateTable();
  c=0;
  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    I->Table[a].index=-1;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    while(s) 
      {
        if(I->Member[s].selection==sele)
          {
            I->Table[a].index=c; /* Mark records  */
            c++;
            break;
          }
        s=SelectorNext(s);
      }
  }
  nAtom=c;
  if(nAtom) {
    nBond = 0;
    bond = VLAlloc(int,nAtom*3);

    for(a=0;a<I->NModel;a++) { /* find bonds wholly contained in the selection */
      obj=I->Obj[a];
      ii1=obj->Bond;
      for(b=0;b<obj->NBond;b++) {
        b1=ii1[0]+obj->SeleBase;
        b2=ii1[1]+obj->SeleBase;
        if((I->Table[b1].index>=0)&&(I->Table[b2].index>=0)) {
          VLACheck(bond,int,nBond*3+1);
          bond[nBond*3]=I->Table[b1].index; /* store what will be the new index */
          bond[nBond*3+1]=I->Table[b2].index;
          bond[nBond*3+2]=ii1[2];
          nBond++;
        }
        ii1+=3;
      }
    }
    
    atInfo = VLAlloc(AtomInfoType,nAtom); /* copy the atom info records */
    c=0;
    for(a=0;a<I->NAtom;a++) {
      if(I->Table[a].index>=0) {
        obj=I->Obj[I->Table[a].model];
        at=I->Table[a].atom;
        VLACheck(atInfo,AtomInfoType,c);
        atInfo[c] = obj->AtomInfo[at];
        atInfo[c].selEntry=0;
        c++;
      }
    }
    
    cs=CoordSetNew();  /* set up a dummy coordinate set for the merge xref */
    cs->NIndex = nAtom;
    cs->fEnumIndices(cs);
    cs->TmpBond = bond; /* load up the bonds */
    cs->NTmpBond = nBond;
    bond=NULL;

    ObjectMoleculeMerge(targ,atInfo,cs,false); /* will free atInfo */
    ObjectMoleculeExtendIndices(targ);
    
    if(!isNew) { /* recreate selection table */
      SelectorUpdateTable(); 
      
      c=0;
      for(a=0;a<I->NAtom;a++) {
        at=I->Table[a].atom;
        I->Table[a].index=-1;
        obj=I->Obj[I->Table[a].model];
        s=obj->AtomInfo[at].selEntry;
        while(s) 
          {
            if(I->Member[s].selection==sele)
              {
                I->Table[a].index=c; /* Mark records  */
                c++;
                break;
              }
            s=SelectorNext(s);
          }
      }
    }
    if(c!=nAtom) ErrFatal("SelectorCreate","inconsistent selection.");
    /* cs->IdxToAtm now has the relevant indexes for the coordinate transfer */
    
    /* get maximum state index */
    nCSet = 0;
    for(a=0;a<I->NModel;a++) { 
      if(nCSet<I->Obj[a]->NCSet)
        nCSet=I->Obj[a]->NCSet;
    }
    for(d=0;d<nCSet;d++) { /* iterate through states */
      if((source<0)||(source==d)) {
        csFlag = false;
        for(a=0;a<I->NAtom;a++)  /* any selected atoms in this state? */
          if(I->Table[a].index>=0) {
            at=I->Table[a].atom;
            obj=I->Obj[I->Table[a].model];
            if(d<obj->NCSet) {
              cs1 = obj->CSet[d];
              if(cs1) {
                if(cs1->AtmToIdx[at]>=0) {
                  csFlag=true;
                  break;
                }
              }
            }
          }
        if(csFlag) { /* copy this coordinate set */
          cs2 = CoordSetNew();
          c = 0;
          cs2->Coord = VLAlloc(float,3*nAtom);
          cs2->AtmToIdx = Alloc(int,targ->NAtom+1);
          for(a=0;a<targ->NAtom;a++) 
            cs2->AtmToIdx[a]=-1;
          cs2->NAtIndex = targ->NAtom;
          cs2->IdxToAtm = Alloc(int,nAtom);
          for(a=0;a<I->NAtom;a++)  /* any selected atoms in this state? */
            if(I->Table[a].index>=0) {
              at=I->Table[a].atom;
              obj=I->Obj[I->Table[a].model];
              if(d<obj->NCSet) {
                cs1 = obj->CSet[d];
                if(cs1) {
                  a1 = cs1->AtmToIdx[at]; /* coord index in existing object */
                  if(a1>=0) {
                    copy3f(cs1->Coord+a1*3,cs2->Coord+c*3);
                    a2 = cs->IdxToAtm[c]; /* actual merged atom index */
                    cs2->IdxToAtm[c] = a2;
                    cs2->AtmToIdx[a2] = c;
                    c++;
                  }
                }
              }
            }
          cs2->IdxToAtm=Realloc(cs2->IdxToAtm,int,c);
          VLASize(cs2->Coord,float,c*3);
          cs2->NIndex = c;
          if(target>=0)
            ts = target;
          else
            ts = d;
          VLACheck(targ->CSet,CoordSet*,ts);
          if(targ->NCSet<=ts) 
            targ->NCSet=ts+1;
          if(targ->CSet[ts])
            targ->CSet[ts]->fFree(targ->CSet[ts]);
          targ->CSet[ts]=cs2;
          cs2->Obj=targ;
        }
      }
    }               
  }
  VLAFreeP(bond); /* null-safe */
  if(cs) cs->fFree(cs);
  if(nAtom) {
    SceneCountFrames();
    ObjectMoleculeSort(targ);
    if(isNew) {
      ObjectSetName((Object*)targ,name);
      ExecutiveManageObject((Object*)targ);
    } else {
      ExecutiveUpdateObjectSelection((Object*)targ);
    }
    SceneChanged();
  } else {
    targ->Obj.fFree((Object*)targ);
  }
}


/*========================================================================*/
int SelectorMatch(int ref,int sele)
{
  SelectorType *I=&Selector;
  return(I->Member[ref].selection==sele);
}
/*========================================================================*/
int SelectorNext(int ref)
{
  SelectorType *I=&Selector;
  return(I->Member[ref].next);
}
/*========================================================================*/
int SelectorIndexByName(char *sname)
{
 OrthoLineType name;
 SelectorType *I=&Selector;
 if(sname[0]=='%')
	strcpy(name,&sname[1]);
 else
	strcpy(name,sname);		  
 return(WordIndex(I->Name,name,1,I->IgnoreCase));
}
/*========================================================================*/
void SelectorPurgeMembers(int sele) 
{
/* NOTE: this routine leaks memory in I->Member array - TODO: fix*/
  int a=0;
  int s=0;
  int l;
  Object *o = NULL;
  void *hidden = NULL;
  ObjectMolecule *obj;

  SelectorType *I=&Selector;
  if(I->Member)
	 while(ExecutiveIterateObject(&o,&hidden))
		{
		  if(o->type==cObjectMolecule)
			 {
				obj=(ObjectMolecule*)o;
				for(a=0;a<obj->NAtom;a++)
				  {
					 l=-1;
					 s=obj->AtomInfo[a].selEntry;
					 while(s)
						{
						  if(I->Member[s].selection==sele)
							 {
								if(l>=0)
								  I->Member[l].next=I->Member[s].next;
								else
								  obj->AtomInfo[a].selEntry=I->Member[s].next;
							 }
                    l=s;
						  s=I->Member[s].next;
						}
				  }
			 }
		}
}
/*========================================================================*/
void SelectorDelete(char *sele) /* should (only) be called by Executive */
{
  SelectorType *I=&Selector;
  int n;
  n=WordIndex(I->Name,sele,999,I->IgnoreCase); /* already exist? */
  if(n>=0) /* get rid of existing selection*/
	 {
		SelectorPurgeMembers(n);
		I->Name[n][0]=32; /*set to blank*/
		I->Name[n][1]=0;
	 }
}
/*========================================================================*/
void SelectorGetTmp(char *input,char *store)
{
  SelectorType *I=&Selector;
  WordType name;
  if(input[0]=='(') {
    sprintf(name,"_%d",I->TmpCounter++);
	 SelectorCreate(name,input,NULL,false);
	 strcpy(store,name);
  } else {
    strcpy(store,input);
  }
}
/*========================================================================*/
void SelectorFreeTmp(char *name)
{
  if(name[0]=='_') ExecutiveDelete(name);
}
/*========================================================================*/
void SelectorCreate(char *sname,char *sele,ObjectMolecule *obj,int quiet) 
{
  SelectorType *I=&Selector;
  int a,m,n;
  int c=0;
  int *atom=NULL;
  OrthoLineType name;
  int ok=true;
  int newFlag=false;
  int flag;
  char buffer[255];

  if(sname[0]=='%')
	 strcpy(name,&sname[1]);
  else
	 strcpy(name,sname);		  
  UtilCleanStr(name);
  if(!name[0])
	 {
		ok=ErrMessage("Select","Invalid selection name.");
		strcpy(name,sname);
		strcat(name,"<--");
		OrthoAddOutput(name);
		OrthoRestorePrompt();
	 }
  if(ok)
	 {
	   if(sele) {
		 atom=SelectorSelect(sele);
		 if(!atom) ok=false;
	   } else {
		 SelectorUpdateTable();
	   }
	 }
  if(ok)
	 {
		n=WordIndex(I->Name,name,999,I->IgnoreCase); /* already exist? */
		if(n>=0) /* get rid of existing selection*/
		  SelectorPurgeMembers(n);
		else
		  {
			 n=I->NSelection;
			 VLACheck(I->Name,WordType,n+1);
			 strcpy(I->Name[n],name);
			 strcpy(I->Name[n+1],""); /*reqd for WordIndex */
			 I->NSelection++;
			 newFlag = true;
		  }
		for(a=0;a<I->NAtom;a++)
		  {
			 flag=false;
			 if(sele) {
				if(atom[a]) flag=true;
			 } else {
				if(I->Obj[I->Table[a].model]==obj) flag=true;
			 }
			 if(flag)
				{
				  c++;
				  I->NMember++;
				  m=I->NMember;
				  VLACheck(I->Member,MemberType,m);
				  I->Member[m].selection=n;
				  I->Member[m].next = 
					 I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].selEntry;
				  I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].selEntry = m;
				}
		  }
		if(!obj) {
		  if(newFlag)
			 ExecutiveManageSelection(name);
		  else
			 ExecutiveSetControlsOff(name);
		}
	 }
  FreeP(atom);
  FreeP(I->Table);
  FreeP(I->Obj);
  I->NAtom=0;
  if(!quiet) {
    if(c) {
      if(name[0]!='_') {
        sprintf(buffer," Selector: selection \"%s\" defined with %d atoms.\n",name,c);
        OrthoAddOutput(buffer);
      }
    } else {
      sprintf(buffer," Selector: no atoms selected.\n");
      OrthoAddOutput(buffer);
    }
  }
}
/*========================================================================*/
int SelectorUpdateTable(void)
{
  int a=0;
  int c=0;
  int modelCnt;
  Object *o = NULL;
  void *hidden = NULL;
  ObjectMolecule *obj;

  SelectorType *I=&Selector;
  FreeP(I->Table);
  FreeP(I->Obj);
  FreeP(I->Vertex);
  FreeP(I->Flag1);
  FreeP(I->Flag2);
  I->NCSet = 0;
  modelCnt=0;
  while(ExecutiveIterateObject(&o,&hidden))
	 {
		if(o->type==cObjectMolecule)
		  {
			 obj=(ObjectMolecule*)o;
			 c+=obj->NAtom;
			 if(I->NCSet<obj->NCSet) I->NCSet=obj->NCSet;
		  }
		modelCnt++;
	 }
  I->Table=Alloc(TableRec,c);
  ErrChkPtr(I->Table);
  I->Obj=Alloc(ObjectMolecule*,modelCnt);
  ErrChkPtr(I->Obj);
  c=0;
  modelCnt=0;
  while(ExecutiveIterateObject(&o,&hidden))
	 {
		I->Obj[modelCnt]=NULL;
		if(o->type==cObjectMolecule)
		  {
			 obj=(ObjectMolecule*)o;
			 I->Obj[modelCnt]=obj;
          obj->SeleBase=c; /* make note of where this object starts */
			 for(a=0;a<obj->NAtom;a++)
				{
				  I->Table[c].model=modelCnt;
				  I->Table[c].atom=a;
				  c++;
				}
		  }
		modelCnt++;
	 }
  I->NModel=modelCnt;
  I->NAtom=c;
  I->Flag1=Alloc(int,c);
  ErrChkPtr(I->Flag1);
  I->Flag2=Alloc(int,c);
  ErrChkPtr(I->Flag2);
  I->Vertex=Alloc(float,c*3);
  ErrChkPtr(I->Vertex);
  return(true);
}
/*========================================================================*/
int *SelectorSelect(char *sele)
{
  WordType *parsed;
  int *result=NULL;
  SelectorUpdateTable();
  parsed=SelectorParse(sele);
  if(parsed)
	 {
		result=SelectorEvaluate(parsed);
		VLAFreeP(parsed);
	 }
  return(result);
}
/*========================================================================*/
int SelectorModulate1(EvalElem *base)
{
  SelectorType *I=&Selector;
  int a,d,e;
  int c=0;
  float dist;
  float *v2;
  CoordSet *cs;
  int ok=true;
  int nCSet;
  MapType *map;
  int i,j,h,k,l;
  int n1,at,idx;
  ObjectMolecule *obj;

  base[1].sele=base[0].sele;
  base->sele=Alloc(int,I->NAtom);
  for(a=0;a<I->NAtom;a++) base[0].sele[a]=false;
  ErrChkPtr(base->sele);
  switch(base[1].code)
	 {
	 case SELE_ARD_:
	 case SELE_EXP_:
		if(!sscanf(base[2].text,"%f",&dist))
		  ok=ErrMessage("Selector","Invalid distance.");
		if(ok)
		  {
			 for(d=0;d<I->NCSet;d++)
				{
				  n1=0;
				  for(a=0;a<I->NAtom;a++) {
					 I->Flag1[a]=false;
					 at=I->Table[a].atom;
					 obj=I->Obj[I->Table[a].model];
                if(d<obj->NCSet) 
                  cs=obj->CSet[d];
                else
                  cs=NULL;
					 if(cs) {
						idx=cs->AtmToIdx[at];
						if(idx>=0) {
						  copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
						  I->Flag1[a]=true;
						  n1++;
						}
					 }
				  }
				  if(n1) {
					 map=MapNewFlagged(-dist,I->Vertex,I->NAtom,NULL,I->Flag1);
					 if(map) {
						MapSetupExpress(map);
						nCSet=SelectorGetArrayNCSet(base[1].sele);
						for(e=0;e<nCSet;e++) {
						  for(a=0;a<I->NAtom;a++) {
							 if(base[1].sele[a])
								{
								  at=I->Table[a].atom;
								  obj=I->Obj[I->Table[a].model];
                          if(e<obj->NCSet) 
                            cs=obj->CSet[e];
                          else
                            cs=NULL;
								  if(cs) {
									 idx=cs->AtmToIdx[at];

									 if(idx>=0) {
										v2 = cs->Coord+(3*idx);
										MapLocus(map,v2,&h,&k,&l);
										i=*(MapEStart(map,h,k,l));
										if(i) {
										  j=map->EList[i++];
										  while(j>=0) {
											 if((!base[0].sele[j])&&
												 ((base[1].code==SELE_EXP_)
												  ||(!base[1].sele[j]))) /*exclude current selection */
												{
												  if(within3f(I->Vertex+3*j,v2,dist)) base[0].sele[j]=true;
												}
											 j=map->EList[i++];
										  }
										}
									 }
								  }
								}
						  }
						}
						MapFree(map);
					 }
				  }
				}
		  }
		break;
	 }
  FreeP(base[1].sele);
  if(DebugSelector&DebugState) {
	 c=0;
	 for(a=0;a<I->NAtom;a++)
		if(base[0].sele[a]) c++;
	 printf("SelectorModulate0: %d atoms selected.\n",c);
  }
  return(ok);
  
}
/*========================================================================*/
int SelectorSelect0(EvalElem *base)
{
  SelectorType *I=&Selector;
  int a,b,flag;
  int c=0;
  short int *vis;

  base->type=STYP_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);
  switch(base->code)
	 {
	 case SELE_NONz:
		for(a=0;a<I->NAtom;a++)
		  base[0].sele[a]=false;
		break;
	 case SELE_HETz:
		for(a=0;a<I->NAtom;a++)
        base[0].sele[a]=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hetatm;
		break;
	 case SELE_HYDz:
		for(a=0;a<I->NAtom;a++)
        base[0].sele[a]=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hydrogen;
		break;
	 case SELE_ALLz:
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a]=true;
			 c++;
		  }
      break;
	 case SELE_VISz:
		for(a=0;a<I->NAtom;a++)
		  {
          flag = false;
          vis = I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].visRep; 
          for(b=0;b<cRepCnt;b++) {
            if(vis[b]) {
              flag=true;
              break;
            }
          }
          base[0].sele[a]=flag;
          if(flag)
            c++;
		  }
		break;
	 }
  if(DebugSelector&DebugState)
	 printf("SelectorSelect0: %d atoms selected.\n",c);
  return(1);
}
/*========================================================================*/
int SelectorSelect1(EvalElem *base)
{
  int a,model,sele,s;
  int c=0;
  int ok=true;
  int rmin,rmax,rtest,index;
  char *p;
  
  SelectorType *I=&Selector;
  ObjectMolecule *obj;
  base->type=STYP_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);
  switch(base->code)
	 {
	 case SELE_IDXs:
		if(sscanf(base[1].text,"%i",&index)!=1)		
		  ok=ErrMessage("Selector","Invalid Range.");
		if(ok) {
		  for(a=0;a<I->NAtom;a++)
			 {
				if(I->Table[a].atom==(index-1)) {
              c++;
				  base[0].sele[a]=true;
            } else {
				  base[0].sele[a]=false;
            }
			 }
		}
		break;
	 case SELE_NAMs:
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
                       I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].name,
                       I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
		  }
		break;
	 case SELE_ELEs:
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
                       I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].elem,
                       I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
		  }
		break;
	 case SELE_SEGs:
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].segi,I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
		  }
		break;
	 case SELE_CHNs:
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
                            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].chain,
                            I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
        }
		break;
	 case SELE_ALTs:
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
                            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].alt,
                            I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
        }
		break;
	 case SELE_RSIs:
		if((p=strstr(base[1].text,":"))) /* range */
		  {
			 *p=' ';
			 if(sscanf(base[1].text,"%i%i",&rmin,&rmax)!=2)
				ok=ErrMessage("Selector","Invalid Range.");
			 if(ok)
				for(a=0;a<I->NAtom;a++)
				  {
					 if(sscanf(I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].resi,
								  "%i",&rtest))
						{
						  if((rtest>=rmin)&&(rtest<=rmax)) {
							 base[0].sele[a]=true;
                      c++;
						  } else 
							 base[0].sele[a]=false;
						}
					 else 
						base[0].sele[a]=false;
				  }
		  }
		else /* not a range */
		  for(a=0;a<I->NAtom;a++)
			 {
				if(WordMatchComma(base[1].text,
                         I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].resi,
                         I->IgnoreCase)<0)
				  {
					 base[0].sele[a]=true;
					 c++;
				  }
				else
				  base[0].sele[a]=false;
			 }
		break;
	 case SELE_RSNs:
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
                            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].resn,
                            I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
		  }
		break;
	 case SELE_SELs:
		sele=WordIndex(I->Name,base[1].text,1,I->IgnoreCase);
		if(sele>=0)
		  {
			 for(a=0;a<I->NAtom;a++)
				{
				  base[0].sele[a]=false;
				  s=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].selEntry;
				  while(s)
					 {
						if(I->Member[s].selection==sele)
						  {
							 base[0].sele[a]=true;
							 c++;
						  }
						s=I->Member[s].next;
					 }
				}
		  } else {
			 for(a=0;a<I->NAtom;a++)
            base[0].sele[a]=false;
          ok=ErrMessage("Selector","Invalid Selection Name.");          
        }
		break;
	 case SELE_MODs:
		model=0;
      obj=(ObjectMolecule*)ExecutiveFindObjectByName(base[1].text);
      if(obj)
        {
          for(a=0;a<I->NModel;a++)
            if(I->Obj[a]==obj)
              {
                model=a+1;
                break;
              }
        }
      if(!model) 
        if(sscanf(base[1].text,"%i",&model)==1)
          {
            if(model<=0)
              model=0;
            else if(model>I->NModel)
              model=0;
            else if(!I->Obj[model])
              model=0;
          }
		if(model)
		  {
			 model--;
			 for(a=0;a<I->NAtom;a++)
				{
				  if(I->Table[a].model==model)
					 {
						base[0].sele[a]=true;
						c++;
					 }
				  else
					 base[0].sele[a]=false;
				}
		  }
		else
		  ErrFatal("SelectorSelect1","Invalid Model");
		break;
	 }
  if(DebugSelector&DebugState)
	 printf("SelectorSelect1:  %d atoms selected.\n",c);
  return(ok);
}
/*========================================================================*/
int SelectorSelect2(EvalElem *base)
{
  int a;
  int c=0;
  int ok=true;
  int oper;
  float comp1;
  AtomInfoType *at1;
  SelectorType *I=&Selector;
  base->type=STYP_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);
  for(a=0;a<I->NAtom;a++) {
    base->sele[a]=0;
  }
  switch(base->code)
	 {
	 case SELE_BVLx:
      oper=WordKey(AtOper,base[1].text,4,I->IgnoreCase);
      if(!oper)
        ok=ErrMessage("Selector","Invalid Operator.");
      if(ok) {
        switch(oper) {
        case SCMP_GTHN:
        case SCMP_LTHN:
        case SCMP_EQAL:
          if (sscanf(base[2].text,"%f",&comp1)!=1) 
            ok=ErrMessage("Selector","Invalid Number");
          break;
        }
        if(ok) {
          switch(oper) {
          case SCMP_GTHN:
            switch(base->code) {
            case SELE_BVLx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->b>comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            }
            break;
          case SCMP_LTHN:
            switch(base->code) {
            case SELE_BVLx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->b<comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            }
            break;
          case SCMP_EQAL:
            switch(base->code) {
            case SELE_BVLx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->b-comp1)<0.0001) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            }
            break;
          }
          break;
        }
      }
    }
  
  if(DebugSelector&DebugState)
	 printf("SelectorSelect2: %d atoms selected.\n",c);
  return(ok);
}
/*========================================================================*/
int SelectorLogic1(EvalElem *base)
{
  int a,b;
  int c=0;
  AtomInfoType *at1,*at2;
  SelectorType *I=&Selector;
  base[0].sele=base[1].sele;
  base[1].sele=NULL;
  base[0].type=STYP_LIST;
  switch(base->code)
	 {
	 case SELE_NOT1:
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a] = ! base[0].sele[a];
			 if(base[0].sele[a])
				c++;
		  }
		break;
	 case SELE_BYR1: /* grossly inefficient */
		for(a=0;a<I->NAtom;a++)
		  {
			 if(base[0].sele[a]) 
			   {
				 at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
				 for(b=0;b<I->NAtom;b++)
				   if(!base[0].sele[b]) 
                 if(I->Table[a].model==I->Table[b].model)
                   {
                     at2=&I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                     if(at1->hetatm==at2->hetatm)
                       if(at1->chain[0]==at2->chain[0])
                         if(WordMatch(at1->resi,at2->resi,I->IgnoreCase)<0)
                           if(WordMatch(at1->segi,at2->segi,I->IgnoreCase)<0)
                             base[0].sele[b]=true;
                   }
			   }
		  }
		break;
	 }
  if(DebugSelector&DebugState)
	 printf("SelectorLogic1: %d atoms selected.\n",c);
  return(1);
}
/*========================================================================*/
int SelectorLogic2(EvalElem *base)
{
  int a,b;
  int c=0;
  SelectorType *I=&Selector;
  AtomInfoType *at1,*at2;

  switch(base[1].code)
	 {
	 case SELE_OR_2:
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a] = base[0].sele[a] || base[2].sele[a];
			 if(base[0].sele[a]) c++;
		  }
		break;
	 case SELE_AND2:
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a] = base[0].sele[a] && base[2].sele[a];
			 if(base[0].sele[a]) c++;
		  }
		break;
	 case SELE_IN_2:
		for(a=0;a<I->NAtom;a++)
		  {
			if(base[0].sele[a]) {
			  at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
			  base[0].sele[a]=0;
			  for(b=0;b<I->NAtom;b++)
				if(base[2].sele[b]) {
				  at2=&I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
				  if((tolower(at1->chain[0]))==(tolower(at2->chain[0])))
					if(WordMatch(at1->name,at2->name,I->IgnoreCase)<0)
					  if(WordMatch(at1->resi,at2->resi,I->IgnoreCase)<0)
						if(WordMatch(at1->resn,at2->resn,I->IgnoreCase)<0)
						  if(WordMatch(at1->segi,at2->segi,I->IgnoreCase)<0)
							base[0].sele[a]=1;
				}
			  if(base[0].sele[a]) c++;
			}
		  }
		break;
	 }
  FreeP(base[2].sele);
  if(DebugSelector&DebugState)
	 printf("SelectorLogic2: %d atoms selected.\n",c);
  return(1);
}
/*========================================================================*/
int *SelectorEvaluate(WordType *word)
{
  int level = 0;
  int depth = 0;
  int a,c = 0;
  int ok=true;
  unsigned int code;
  int valueFlag = 0; /* are we expecting? */
  int *result = NULL;
  int opFlag;
  char *q,*cc1,*cc2;
  OrthoLineType line;
  EvalElem Stack[SelectorMaxDepth],*e;
  SelectorType *I=&Selector;
  /* converts all keywords into code, adds them into a operation list */
  while(ok&&word[c][0])
	 {
		switch(word[c][0])
		  {
		  case '(':
			 if(valueFlag)	ok=ErrMessage("Selector","Misplaced (.");
			 if(ok) level++;
			 break;
		  case ')':
			 if(valueFlag)	ok=ErrMessage("Selector","Misplaced ).");
			 if(ok)
				{
				level--;
				if(level<0)
				  ok=ErrMessage("Selector","Syntax error.");
				}
			 if(ok&&depth) Stack[depth].level--;
			 break;
		  default:
			 if(valueFlag>0) /* standard operand */
				{
				  depth++;
				  e=Stack+depth;
				  e->level=level;
				  e->type=STYP_VALU;
              cc1 = word[c];
              cc2 = e->text;
              while(*cc1) { /* remove embedded quotes if any */
                if((*cc1!=34)&&(*cc1!=39)) {
                  *(cc2++)=*(cc1++);
                } else 
                  cc1++;
              }
              *cc2=0;
				  valueFlag--;
				}
			 else if(valueFlag<0) /* operation parameter i.e. around X<-- */
				{
				  depth++;
				  e=Stack+depth;
				  e->level=level;
				  e->type=STYP_PVAL;
				  strcpy(e->text,word[c]);
				  valueFlag++;
				} 
			 else
				{
				  code=WordKey(Keyword,word[c],4,I->IgnoreCase);
              if(DebugState&DebugSelector)
                printf("code %x\n",code);
				  if(code) 
					 { /* this is a known operation */
						depth++;
						e=Stack+depth;
						e->level=level;
						e->code=code;
                  e->type=(e->code&0xF);
						switch(e->type)
						  {
						  case STYP_SEL0:
							 valueFlag=0;
							 break;
						  case STYP_SEL1:
							 valueFlag=1;
							 break;
						  case STYP_SEL2:
							 valueFlag=2;
							 break;
						  case STYP_OPR1:
							 valueFlag=0;
							 break;
						  case STYP_OPR2:
							 valueFlag=0;
							 break;
                    case STYP_PRP1: 
							 valueFlag=-1;
							 break;
                    }
					 } else if(SelectorIndexByName(word[c])>=0) {
						depth++;
						e=Stack+depth;
						e->level=level;
						e->code=SELE_SELs;
                  e->type=STYP_SEL1;
                  valueFlag=1;
                  c--;
                } else {
                  ok=ErrMessage("Selector","Unknown keyword or selection.");
                }
            }
          break;
		  }
		if(ok)
		  do 
			 {
				opFlag=false;
				if(ok)
				  if(depth>0)
					 if((!opFlag)&&(Stack[depth].type==STYP_SEL0))
						{
						  opFlag=true;
						  ok=SelectorSelect0(&Stack[depth]);
						}
				if(ok)
				  if(depth>1)
					 if(Stack[depth].level==Stack[depth-1].level)
						{
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_SEL1)
							  &&(Stack[depth].type==STYP_VALU))
							 { /* 1 argument selection operator */
								opFlag=true;
								ok=SelectorSelect1(&Stack[depth-1]);
								depth--;
							 }
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_OPR1)
							  &&(Stack[depth].type==STYP_LIST))
							 { /* 1 argument logical operator */
								opFlag=true;
								ok=SelectorLogic1(&Stack[depth-1]);
								depth--;
							 }
						}
				if(ok)
				  if(depth>2)
					 if((Stack[depth].level==Stack[depth-1].level)&&
						 (Stack[depth].level==Stack[depth-2].level))
						{
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_OPR2)
							  &&(Stack[depth].type==STYP_LIST)
							  &&(Stack[depth-2].type==STYP_LIST))
							 { /* 2 argument logical operator */
								ok=SelectorLogic2(&Stack[depth-2]);
								depth-=2;
							 }
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_PRP1)
							  &&(Stack[depth].type==STYP_PVAL)
							  &&(Stack[depth-2].type==STYP_LIST))
							 { /* 2 argument logical operator */
								ok=SelectorModulate1(&Stack[depth-2]);
								depth-=2;
							 }
						  if(ok&&(!opFlag)&&(Stack[depth-2].type==STYP_SEL2)
							  &&(Stack[depth-1].type==STYP_VALU)
							  &&(Stack[depth].type==STYP_VALU))
							 { /* 2 argument value operator */
								ok=SelectorSelect2(&Stack[depth-2]);
								depth-=2;
							 }
						}
				if(ok)
				  if(depth>3)
					 if((Stack[depth].level==Stack[depth-1].level)&&
						 (Stack[depth].level==Stack[depth-2].level)&&
						 (Stack[depth].level==Stack[depth-3].level))
						{
						  if(ok&&(!opFlag)&&(Stack[depth-3].type==STYP_SEL3)
							  &&(Stack[depth].type==STYP_VALU)
							  &&(Stack[depth-1].type==STYP_VALU)
							  &&(Stack[depth-2].type==STYP_VALU))
							 { /* 2 argument logical operator */
                        /*								ok=SelectorSelect3(&Stack[depth-3]);*/
								depth-=3;
							 }
						}
			 }
        while(ok&&opFlag); /* part of a do while */
		if(ok) c++; /* go onto next word */
	 }
  if(ok)
	 {
      if(level>0)
		ok=ErrMessage("Selector","Malformed selection.");		
      else if(depth!=1)
        ok=ErrMessage("Selector","Malformed selection.");
      else if(Stack[depth].type!=STYP_LIST)
        ok=ErrMessage("Selector","Invalid selection.");
      else
        result=Stack[depth].sele; /* return the selection list */
	 }
  if(!ok)
	 {
		for(a=0;a<depth;a++)
		  if(Stack[a].type==STYP_LIST)
			 FreeP(Stack[a].sele);
		depth=0;
		q=line;
		*q=0;
		for(a=0;a<=c;a++)
		  {
			 if(a&&word[a][0])
				q=UtilConcat(q," ");
			 q=UtilConcat(q,word[a]);
		  }
		q=UtilConcat(q,"<--");
		OrthoAddOutput(line);
		OrthoRestorePrompt();
	 }
  return(result);
}

/*========================================================================*/
WordType *SelectorParse(char *s) {

  /* break a selection down into its constituent strings and
	  return them in a WordType VLA, null string terminated */

  WordType *r = NULL;
  int c=0;
  int w_flag=false;
  char *p = s;
  char *q = NULL;
  r=VLAlloc(WordType,100);
  while(*p) 
	 {
		if(w_flag) /* currently in a word */
		  {
			 switch(*p)
				{
				case ' ':
				  *q=0;
				  w_flag=false;
				  break;
            case ';': /* special word terminator */
				  *q++=*p;
				  *q=0;
              w_flag=false;
              break;
				case '!': /* single words */ 
				case '&': 
				case '|': 
				case '(': 
				case ')':
				case '>': 
				case '<':
            case '=':
				case '%':
				  *q=0;
				  q=r[c];
				  c++;
				  VLACheck(r,WordType,c);
				  *q++=*p;
				  *q=0;
				  w_flag=false;
				  break;
				default:
				  *q++=*p;
				  break;
				}
		  }
		else /*outside a word*/
		  {
			 switch(*p)
				{
				case '*': /* special case */
				  *q=0;
				  q=r[c];
				  c++;
				  VLACheck(r,WordType,c);
				  *q++='+';
				  *q=0;
				  w_flag=false;
              break;
				case '!': /* single words */ 
				case '&': 
				case '|': 
				case '(': 
				case ')':
				case '>': 
				case '<':
            case '=':
				case '%':
				  q=r[c];
				  *q++=(*p);
				  *q=0;
				  c++;
				  VLACheck(r,WordType,c);
				  break;
            case ' ':
              break;
				default:
				  w_flag=true;
				  q=r[c];
				  *q++=*p;
				  c++;
				  VLACheck(r,WordType,c);
				  break;
				}
		  }
		p++;
	 }
  /* end current word */
  if(w_flag) 
	 *q=0;

  /* null strings terminate the list*/
  q=r[c];
  *q=0;
  if(DebugState&DebugSelector)
	 {
		c=0;
		while(r[c][0])
		  {
			 printf("word: %s\n",r[c]);
			 c++;
		  }
	 }
  return(r);
}
/*========================================================================*/
void SelectorFree(void)
{
  SelectorType *I = &Selector;
  FreeP(I->Table);
  FreeP(I->Obj);
  VLAFreeP(I->Member);
  VLAFreeP(I->Name);
  FreeP(I->Flag1);
  FreeP(I->Flag2);
  FreeP(I->Vertex);

}
/*========================================================================*/
void SelectorInit(void)
{
  SelectorType *I = &Selector;
  I->NSelection = 0;
  I->TmpCounter = 0;
  I->Name = VLAlloc(WordType,10);
  I->Member = (MemberType*)VLAMalloc(1000,sizeof(MemberType),5,true);
  I->NMember=0;
  I->NCSet=0;
  I->Table=NULL;
  I->Obj=NULL;
  I->IgnoreCase=true;
  I->Flag1=NULL;
  I->Flag2=NULL;
  I->Vertex=NULL;
}
/*========================================================================*/


DistSet *SelectorGetDistSet(int sele1,int state1,int sele2,int state2,int mode,float cutoff)
{
  SelectorType *I=&Selector;
  int *vla=NULL;
  int c;
  float result=0.0;
  float sumVDW,dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at1,at2;
  CoordSet *cs1,*cs2;
  DistSet *ds;
  ObjectMolecule *obj1,*obj2;
  int idx1,idx2;
  int a;
  int nv = 0;
  float *vv,*vv0,*vv1;

  ds = DistSetNew();
  vv = VLAlloc(float,10000);
  SelectorUpdateTable();
  if(cutoff<0) cutoff = 1000.0;
  c=SelectorGetInterstateVLA(sele1,state1,sele2,state2,cutoff,&vla);
  for(a=0;a<c;a++) {
    a1=vla[a*2];
    a2=vla[a*2+1];

    if(a1!=a2) {
      at1=I->Table[a1].atom;
      at2=I->Table[a2].atom;
      
      obj1=I->Obj[I->Table[a1].model];
      obj2=I->Obj[I->Table[a2].model];
      
      if(state1<obj1->NCSet&&state2<obj2->NCSet) {
        cs1=obj1->CSet[state1];
        cs2=obj2->CSet[state2];
        if(cs1&&cs2) { 
          
          ai1=obj1->AtomInfo+at1;
          ai2=obj2->AtomInfo+at2;
          
          idx1=cs1->AtmToIdx[at1];
          idx2=cs2->AtmToIdx[at2];
          
          /*          sumVDW=ai1->vdw+ai2->vdw;*/

          dist=diff3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2);
          
          if(dist<cutoff) {
            VLACheck(vv,float,(nv*3)+5);
            vv0 = vv+ (nv*3);
            vv1 = cs1->Coord+3*idx1;
            *(vv0++) = *(vv1++);
            *(vv0++) = *(vv1++);
            *(vv0++) = *(vv1++);
            vv1 = cs2->Coord+3*idx2;
            *(vv0++) = *(vv1++);
            *(vv0++) = *(vv1++);
            *(vv0++) = *(vv1++);

            result+=((sumVDW-dist)/2.0);
            nv+=2;
          }
        }
      }
    }
  }
  VLAFreeP(vla);
  ds->NIndex = nv;
  ds->Coord = vv;
  return(ds);
}

/*========================================================================*/


/*

example selections
cas
backbone
(model 1 and backbone)
(name ca)
(resi 123)
(resi 200:400 and chain A)
(model a)

 */

/* In order for selections to be robust during atom insertions
   deletions, they are stored not as lists of selected atoms, but
   rather in the inverse - as atoms with selection membership
   information */

/* Each atom points into the selection heap, a series of linked
   entries where each member is selection index */

/* Definition of the selection language:

	<sele> = [(] [<not>] <prop> <val-range> [<qual> [<qual-range>] ] [)] { <SEL1> <sele> }

	Example selections:
	
   name ca
	( name ca )
	name ca around 5 {all atoms within 5 angstroms of any ca atom) }
	( resi 10 ) 
	resi 10:50
	resi 10A:50A
	resi 10,50
	chain A
	segi A
	model 1 and segi B
	model a and segi B
	not name ca

*/

/* Selection processing, left to right by default, but with parenthesis for grouping
	stack based functional language processing
	each stack level has a full selection matrix
	logical binary operators (or,and) and negation SEL1ation solely on these matrices.
*/	

/*

(not (name ca or name c) and (name s around 5 and (name c) )) around 6

0:
1: not
2: name 
2: ca

0:

1: not
2:<name ca>
not name ca around 5

force compute

0:
1: not

attrib b < 0 

*/


