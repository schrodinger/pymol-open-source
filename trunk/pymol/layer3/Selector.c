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

#define SELE_VALU 0
#define SELE_OPR1 1
#define SELE_OPR2 2
#define SELE_SEL0 3
#define SELE_SEL1 4
#define SELE_SEL2 5
#define SELE_LIST 6
#define SELE_PRP1 7
#define SELE_SEL3 8
#define SELE_PVAL 0

static WordType Keyword[] = 
{
  "not",      "NOT1",
  "!",        "NOT1",
  "byresi",   "BYR1",
  "b;",       "BYR1",
  "and",      "AND2",
  "&",        "AND2",
  "or",       "OR_2",
  "|",        "OR_2",
  "in",       "IN_2",
  "all",      "ALLz", /* 0 parameter */
  "+",        "ALLz", /* 0 parameter */
  "none",     "NONz", /* 0 parameter */
  "hetatm",   "HETz", /* 0 parameter */
  "het",      "HETz", /* 0 parameter */
  "hydro",    "HYDz", /* 0 parameter */
  "h;",       "HYDz", /* 0 parameter */
  "around",   "ARD_", /* 1 parameter */
  "a;",       "ARD_", /* 1 parameter */
  "expand",   "EXP_", /* 1 parameter */
  "x;",       "EXP_", /* 1 parameter */
  "name",     "NAMs",
  "n;",       "NAMs",
  "elem",     "ELEs",
  "e;",       "ELEs",
  "resi",     "RSIs",
  "i;",       "RSIs",
  "chain",    "CHNs",
  "c;",       "CHNs",
  "segi",     "SEGs",
  "s;",       "SEGs",
  "model",    "MODs",
  "m;",       "MODs",
  "index",    "IDXs",
  "i;",       "IDXs",
  "resn",     "RSNs",
  "r;",       "RSNs",
  "%",        "SELs",
  "b",        "BVLx", /* 2 operand selection operator */ 
  ""
};

static WordType AtOper[] = 
{
  ">",      "GTHN",
  "<",      "LTHN",
  "in",     "RANG",
  "=",      "EQAL",
  ""
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
                I->Table[a].index=c+1;
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
	 SelectorCreate(name,input,NULL);
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
void SelectorCreate(char *sname,char *sele,ObjectMolecule *obj) 
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
  if(c&&name[0]!='_') {
    sprintf(buffer," Selector: selection \"%s\" defined with %d atoms.\n",name,c);
    OrthoAddOutput(buffer);
  } else if(!c) {
    sprintf(buffer," Selector: no atoms selected.\n");
    OrthoAddOutput(buffer);
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
	 case 'ARD_':
	 case 'EXP_':
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
												 ((base[1].code=='EXP_')
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
	 printf("SelectorModulate0: %c%c%c%c : %d atoms selected.\n",
			  base[1].code>>24,base[1].code>>16,
			  base[1].code>>8,base[1].code&0xFF,c);
  }
  return(ok);
  
}
/*========================================================================*/
int SelectorSelect0(EvalElem *base)
{
  SelectorType *I=&Selector;
  int a;
  int c=0;
  base->type=SELE_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);
  switch(base->code)
	 {
	 case 'NONz':
		for(a=0;a<I->NAtom;a++)
		  base[0].sele[a]=false;
		break;
	 case 'HETz':
		for(a=0;a<I->NAtom;a++)
        base[0].sele[a]=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hetatm;
		break;
	 case 'HYDz':
		for(a=0;a<I->NAtom;a++)
        base[0].sele[a]=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hydrogen;
		break;
	 case 'ALLz':
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a]=true;
			 c++;
		  }
		break;
	 }
  if(DebugSelector&DebugState)
	 printf("SelectorSelect0: %c%c%c%c : %d atoms selected.\n",
			  base[0].code>>24,base[0].code>>16,
			  base[0].code>>8,base[0].code&0xFF,c);
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
  base->type=SELE_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);
  switch(base->code)
	 {
	 case 'IDXs':
		if(sscanf(base[1].text,"%i",&index)!=1)		
		  ok=ErrMessage("Selector","Invalid Range.");
		if(ok) {
		  for(a=0;a<I->NAtom;a++)
			 {
				if(I->Table[a].atom==(index-1))
				  base[0].sele[a]=true;
				else 
				  base[0].sele[a]=false;
				c++;
			 }
		}
		break;
	 case 'NAMs':
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
	 case 'ELEs':
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
	 case 'SEGs':
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
	 case 'CHNs':
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
	 case 'RSIs':
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
						  if((rtest>=rmin)&&(rtest<=rmax))
							 base[0].sele[a]=true;
						  else 
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
	 case 'RSNs':
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
	 case 'SELs':
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
	 case 'MODs':
		model=0;
		if(sscanf(base[1].text,"%i",&model)==1)
		  {
			 if(model<=0)
				model=0;
			 else if(model>I->NModel)
				model=0;
			 else if(!I->Obj[model])
				model=0;
		  }
		else
		  model=0;
		if(!model)
		  {
			 obj=(ObjectMolecule*)ExecutiveFindObjectByName(base[1].text);
			 if(!obj)
				model=0;
			 else
				{
				  for(a=0;a<I->NModel;a++)
					 if(I->Obj[a]==obj)
						{
						  model=a+1;
						  break;
						}
				}
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
	 printf("SelectorSelect1: %c%c%c%c %s: %d atoms selected.\n",
			  base[0].code>>24,base[0].code>>16,
			  base[0].code>>8,base[0].code&0xFF,base[1].text,c);
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
  base->type=SELE_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);
  for(a=0;a<I->NAtom;a++) {
    base->sele[a]=0;
  }
  switch(base->code)
	 {
	 case 'BVLx':
      oper=WordChoose(AtOper,base[1].text,4,I->IgnoreCase);
      if(!oper)
        ok=ErrMessage("Selector","Invalid Operator.");
      if(ok) {
        switch(oper) {
        case 'GTHN':
        case 'LTHN':
        case 'EQAL':
          if (sscanf(base[2].text,"%f",&comp1)!=1) 
            ok=ErrMessage("Selector","Invalid Number");
          break;
        }
        if(ok) {
          switch(oper) {
          case 'GTHN':
            switch(base->code) {
            case 'BVLx':
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
          case 'LTHN':
            switch(base->code) {
            case 'BVLx':
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
          case 'EQAL':
            switch(base->code) {
            case 'BVLx':
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
	 printf("SelectorSelect2: %c%c%c%c %s: %d atoms selected.\n",
			  base[0].code>>24,base[0].code>>16,
			  base[0].code>>8,base[0].code&0xFF,base[1].text,c);
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
  base[0].type=SELE_LIST;
  switch(base->code)
	 {
	 case 'NOT1':
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a] = ! base[0].sele[a];
			 if(base[0].sele[a])
				c++;
		  }
		break;
	 case 'BYR1': /* grossly inefficient */
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
	 printf("SelectorLogic1: %c%c%c%c: %d atoms selected.\n",
			  base[0].code>>24,base[0].code>>16,
			  base[0].code>>8,base[0].code&0xFF,c);
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
	 case 'OR_2':
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a] = base[0].sele[a] || base[2].sele[a];
			 if(base[0].sele[a]) c++;
		  }
		break;
	 case 'AND2':
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a] = base[0].sele[a] && base[2].sele[a];
			 if(base[0].sele[a]) c++;
		  }
		break;
	 case 'IN_2':
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
	 printf("SelectorLogic2: %c%c%c%c: %d atoms selected.\n",
			  base[1].code>>24,base[1].code>>16,
			  base[1].code>>8,base[1].code&0xFF,c);
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
  int opFlag,lt;
  char *q;
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
				  e->type=SELE_VALU;
				  strcpy(e->text,word[c]);
              if((e->text[0]==34)||(e->text[0]==39)) { /* remove surrounding quotes if any */
                strcpy(e->text,word[c]+1);
              }
              lt=strlen(e->text);
              if(lt) {
                lt--;
                if((e->text[lt]==34)||(e->text[lt]==39)) {
                  e->text[lt]=0;
                }
              }
				  valueFlag--;
				}
			 else if(valueFlag<0) /* operation parameter i.e. around X<-- */
				{
				  depth++;
				  e=Stack+depth;
				  e->level=level;
				  e->type=SELE_PVAL;
				  strcpy(e->text,word[c]);
				  valueFlag++;
				} 
			 else
				{
				  code=WordChoose(Keyword,word[c],4,I->IgnoreCase);
				  if(code) 
					 { /* this is a known operation */
						depth++;
						e=Stack+depth;
						e->level=level;
						e->code=code;
						switch((e->code)&0xFF)
						  {
						  case 'z':
							 e->type=SELE_SEL0;
							 valueFlag=0;
							 break;
						  case 's':
							 e->type=SELE_SEL1;
							 valueFlag=1;
							 break;
						  case 'x':
							 e->type=SELE_SEL2;
							 valueFlag=2;
							 break;
						  case '1':
							 e->type=SELE_OPR1;
							 break;
						  case '2':
							 e->type=SELE_OPR2;
							 break;
						  case '_':
							 e->type=SELE_PRP1;
							 valueFlag=-1;
							 break;
                    }
					 } else if(SelectorIndexByName(word[c])>=0) {
						depth++;
						e=Stack+depth;
						e->level=level;
						e->code='SELs';
                  e->type=SELE_SEL1;
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
					 if((!opFlag)&&(Stack[depth].type==SELE_SEL0))
						{
						  opFlag=true;
						  ok=SelectorSelect0(&Stack[depth]);
						}
				if(ok)
				  if(depth>1)
					 if(Stack[depth].level==Stack[depth-1].level)
						{
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==SELE_SEL1)
							  &&(Stack[depth].type==SELE_VALU))
							 { /* 1 argument selection operator */
								opFlag=true;
								ok=SelectorSelect1(&Stack[depth-1]);
								depth--;
							 }
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==SELE_OPR1)
							  &&(Stack[depth].type==SELE_LIST))
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
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==SELE_OPR2)
							  &&(Stack[depth].type==SELE_LIST)
							  &&(Stack[depth-2].type==SELE_LIST))
							 { /* 2 argument logical operator */
								ok=SelectorLogic2(&Stack[depth-2]);
								depth-=2;
							 }
						  if(ok&&(!opFlag)&&(Stack[depth-1].type==SELE_PRP1)
							  &&(Stack[depth].type==SELE_PVAL)
							  &&(Stack[depth-2].type==SELE_LIST))
							 { /* 2 argument logical operator */
								ok=SelectorModulate1(&Stack[depth-2]);
								depth-=2;
							 }
						  if(ok&&(!opFlag)&&(Stack[depth-2].type==SELE_SEL2)
							  &&(Stack[depth-1].type==SELE_VALU)
							  &&(Stack[depth].type==SELE_VALU))
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
						  if(ok&&(!opFlag)&&(Stack[depth-3].type==SELE_SEL3)
							  &&(Stack[depth].type==SELE_VALU)
							  &&(Stack[depth-1].type==SELE_VALU)
							  &&(Stack[depth-2].type==SELE_VALU))
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
      else if(Stack[depth].type!=SELE_LIST)
        ok=ErrMessage("Selector","Invalid selection.");
      else
        result=Stack[depth].sele; /* return the selection list */
	 }
  if(!ok)
	 {
		for(a=0;a<depth;a++)
		  if(Stack[a].type==SELE_LIST)
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


