
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
#include"Map.h"
#include"Vector.h"
#include"Debug.h"
#include"Err.h"
#include"Word.h"
#include"Util.h"
#include"PConv.h"
#include"P.h"

#include"MemoryDebug.h"
#include"Selector.h"
#include"Executive.h"
#include"ObjectMolecule.h"
#include"CoordSet.h"
#include"DistSet.h"
#include"Word.h"
#include"Scene.h"
#include"CGO.h"

#define SelectorMaxDepth 1000

#define cSelectorTmpPrefix "_sel_tmp_"

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
  int branch;
  float f1;
} TableRec;

typedef struct {
  WordType *Name;
  int *ID;
  int NSelection,NActive;
  int TmpCounter;
  MemberType *Member;
  int NMember;
  int FreeMember;
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
int SelectorOperator22(EvalElem *base);
int *SelectorEvaluate(WordType *word);
WordType *SelectorParse(char *s);
void SelectorPurgeMembers(int sele);
int SelectorUpdateTableSingleObject(ObjectMolecule *obj);
int  SelectorEmbedSelection(int *atom, char *name, ObjectMolecule *obj);
void SelectorClean(void);
void SelectorDeletePrefixSet(char *s);
int *SelectorGetIndexVLA(int sele);
int *SelectorApplyMultipick(Multipick *mp);

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
#define STYP_OP22 9   /* sele oper arg1 arg2 sele */

/*                  code   |   type    | priority */

#define SELE_NOT1 ( 0x0100 | STYP_OPR1 | 0x60 )
#define SELE_BYR1 ( 0x0200 | STYP_OPR1 | 0x10 )
#define SELE_AND2 ( 0x0300 | STYP_OPR2 | 0x50 )
#define SELE_OR_2 ( 0x0400 | STYP_OPR2 | 0x30 )
#define SELE_IN_2 ( 0x0500 | STYP_OPR2 | 0x30 )
#define SELE_ALLz ( 0x0600 | STYP_SEL0 | 0x80 )
#define SELE_NONz ( 0x0700 | STYP_SEL0 | 0x80 )
#define SELE_HETz ( 0x0800 | STYP_SEL0 | 0x70 )
#define SELE_HYDz ( 0x0900 | STYP_SEL0 | 0x80 )
#define SELE_VISz ( 0x0A00 | STYP_SEL0 | 0x80 )
#define SELE_ARD_ ( 0x0B00 | STYP_PRP1 | 0x20 )
#define SELE_EXP_ ( 0x0C00 | STYP_PRP1 | 0x20 )
#define SELE_NAMs ( 0x0D00 | STYP_SEL1 | 0x70 )
#define SELE_ELEs ( 0x0E00 | STYP_SEL1 | 0x70 )
#define SELE_RSIs ( 0x0F00 | STYP_SEL1 | 0x70 )
#define SELE_CHNs ( 0x1000 | STYP_SEL1 | 0x70 )
#define SELE_SEGs ( 0x1100 | STYP_SEL1 | 0x70 )
#define SELE_MODs ( 0x1200 | STYP_SEL1 | 0x70 ) 
#define SELE_IDXs ( 0x1300 | STYP_SEL1 | 0x70 )
#define SELE_RSNs ( 0x1400 | STYP_SEL1 | 0x70 )
#define SELE_SELs ( 0x1500 | STYP_SEL1 | 0x70 )
#define SELE_BVLx ( 0x1600 | STYP_SEL2 | 0x70 )
#define SELE_ALTs ( 0x1700 | STYP_SEL1 | 0x70 )
#define SELE_FLGs ( 0x1800 | STYP_SEL1 | 0x70 )
#define SELE_GAP_ ( 0x1900 | STYP_PRP1 | 0x70 )
#define SELE_TTYs ( 0x1A00 | STYP_SEL1 | 0x70 )  
#define SELE_NTYs ( 0x1B00 | STYP_SEL1 | 0x70 )
#define SELE_PCHx ( 0x1C00 | STYP_SEL2 | 0x70 )  
#define SELE_FCHx ( 0x1D00 | STYP_SEL2 | 0x70 )
#define SELE_ID_s ( 0x1E00 | STYP_SEL1 | 0x70 )
#define SELE_BNDz ( 0x1F00 | STYP_SEL0 | 0x70 )
#define SELE_LIK2 ( 0x2000 | STYP_OPR2 | 0x30 )
#define SELE_NGH1 ( 0x2100 | STYP_OPR1 | 0x10 )
#define SELE_QVLx ( 0x2200 | STYP_SEL2 | 0x70 )
#define SELE_BYO1 ( 0x2300 | STYP_OPR1 | 0x10 )
#define SELE_SSTs ( 0x2400 | STYP_SEL1 | 0x70 )
#define SELE_STAs ( 0x2500 | STYP_SEL1 | 0x70 )
#define SELE_PREz ( 0x2500 | STYP_SEL0 | 0x70 )
#define SELE_WIT_ ( 0x2600 | STYP_OP22 | 0x20 ) 

#define SEL_PREMAX 0x8

static WordKeyValue Keyword[] = 
{
  {  "not",      SELE_NOT1 },
  {  "!",        SELE_NOT1 },
  {  "neighbor", SELE_NGH1 },
  {  "nbr;",     SELE_NGH1 }, /* deprecated */
  {  "nbr.",     SELE_NGH1 },
  {  "byresi",   SELE_BYR1 },
  {  "byres",    SELE_BYR1 },
  {  "br;",      SELE_BYR1 },/* deprecated */
  {  "br.",      SELE_BYR1 },
  {  "b;",       SELE_BYR1 }, /* deprecated */
  {  "byobj",    SELE_BYO1 },
  {  "byobject", SELE_BYO1 },
  {  "bo;",      SELE_BYO1 },/* deprecated */
  {  "bo.",      SELE_BYO1 },
  {  "and",      SELE_AND2 },
  {  "&",        SELE_AND2 },
  {  "or",       SELE_OR_2 },
  {  "|",        SELE_OR_2 },
  {  "in",       SELE_IN_2 },
  {  "like",     SELE_LIK2 },
  {  "l;",       SELE_LIK2 },
  {  "l.",       SELE_LIK2 },
  {  "all",      SELE_ALLz }, /* 0 parameter */
  /*  {  "+",        SELE_ALLz },*/ /* 0 parameter */
  {  "none",     SELE_NONz }, /* 0 parameter */
  {  "hetatm",   SELE_HETz }, /* 0 parameter */
  {  "het",      SELE_HETz }, /* 0 parameter */
  {  "hydro",    SELE_HYDz }, /* 0 parameter */
  {  "h;",       SELE_HYDz }, /* deprecated */
  {  "h.",       SELE_HYDz }, /* 0 parameter */
  {  "visible",  SELE_VISz }, /* 0 parameter */
  {  "v;",       SELE_VISz }, /* 0 parameter */
  {  "v.",       SELE_VISz }, /* 0 parameter */
  {  "around",   SELE_ARD_ }, /* 1 parameter */
  {  "a;",       SELE_ARD_ }, /* deprecated */
  {  "a.",       SELE_ARD_ }, /* 1 parameter */
  {  "expand",   SELE_EXP_ }, /* 1 parameter */
  {  "x;",       SELE_EXP_ }, /* 1 parameter */
  {  "x.",       SELE_EXP_ }, /* 1 parameter */
  {  "name",     SELE_NAMs },
  {  "n;",       SELE_NAMs },/* deprecated */
  {  "n.",       SELE_NAMs },
  {  "symbol",   SELE_ELEs },
  {  "element",  SELE_ELEs },
  {  "elem",     SELE_ELEs },
  {  "e;",       SELE_ELEs },/* deprecated */
  {  "e.",       SELE_ELEs },
  {  "resi",     SELE_RSIs },
  {  "resid",    SELE_RSIs },
  {  "i;",       SELE_RSIs },/* deprecated */
  {  "i.",       SELE_RSIs },
  {  "alt",      SELE_ALTs },
  {  "flag",     SELE_FLGs },
  {  "f;",       SELE_FLGs },/* deprecated */
  {  "f.",       SELE_FLGs },
  {  "gap",      SELE_GAP_ },
  {  "partial_charge",SELE_PCHx },
  {  "pc;",      SELE_PCHx },/* deprecated */
  {  "pc.",      SELE_PCHx },
  {  "formal_charge", SELE_FCHx },
  {  "fc;",      SELE_FCHx },/* deprecated */
  {  "fc.",      SELE_FCHx },
  {  "numeric_type",SELE_NTYs },
  {  "nt;",      SELE_NTYs },/* deprecated */
  {  "nt.",      SELE_NTYs },
  {  "text_type",SELE_TTYs },
  {  "tt;",      SELE_TTYs },/* deprecated */
  {  "tt.",      SELE_TTYs },
  {  "chain",    SELE_CHNs },
  {  "c;",       SELE_CHNs },/* deprecated */
  {  "c.",       SELE_CHNs },
  {  "bonded",   SELE_BNDz },
  {  "segi",     SELE_SEGs },
  {  "s;",       SELE_SEGs },/* deprecated */
  {  "s.",       SELE_SEGs },
  {  "ss",       SELE_SSTs },
  {  "state",    SELE_STAs },
  {  "model",    SELE_MODs },
  {  "m;",       SELE_MODs },/* deprecated */
  {  "m.",       SELE_MODs },
  {  "index",    SELE_IDXs },
  {  "id",       SELE_ID_s },
  {  "resn",     SELE_RSNs },
  {  "within",   SELE_WIT_ },
  {  "present",  SELE_PREz },
  {  "w.",       SELE_WIT_ },
  {  "r;",       SELE_RSNs },/* deprecated */
  {  "r.",       SELE_RSNs },
  {  "%",        SELE_SELs },
  {  "b",        SELE_BVLx, }, /* 2 operand selection operator */ 
  {  "q",        SELE_QVLx, }, /* 2 operand selection operator */ 
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

static int BondInOrder(BondType *a,int b1,int b2);
static int BondCompare(BondType *a,BondType *b);

int SelectorWalkTree(int *atom,int *comp,int *toDo,int **stk,
                     int stkDepth,ObjectMolecule *obj,int sele1,int sele2);

/*========================================================================*/

int SelectorGetPairIndices(int sele1,int state1,int sele2,int state2,
                           int mode,float cutoff, float h_angle,
                           int **indexVLA, ObjectMolecule ***objVLA)
{
  SelectorType *I=&Selector;
  int *vla=NULL;
  int c;
  float dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at1,at2;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj1,*obj2;
  int idx1,idx2;
  int a;
  int dist_cnt = 0;
  float dir[3];
  float v1[3],v2[3];
  int flag;
  float angle_cutoff=0.0;

  if(mode==1) {
    angle_cutoff = cos(PI*h_angle/180.8);
  }
  SelectorUpdateTable();
  if(cutoff<0) cutoff = 1000.0;
  c=SelectorGetInterstateVLA(sele1,state1,sele2,state2,cutoff,&vla);
  (*indexVLA)=VLAlloc(int,1000);
  (*objVLA)=VLAlloc(ObjectMolecule*,1000);

  CGOReset(DebugCGO);

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

          if(obj1->DiscreteFlag) {
            if(cs1==obj1->DiscreteCSet[at1]) {
              idx1=obj1->DiscreteAtmToIdx[at1];
            } else {
              idx1=-1;
            }
          } else {
            idx1=cs1->AtmToIdx[at1];
          }
          
          if(obj2->DiscreteFlag) {
            if(cs2==obj2->DiscreteCSet[at2]) {
              idx2=obj2->DiscreteAtmToIdx[at2];
            } else {
              idx2=-1;
            }
            
          } else {
            idx2=cs2->AtmToIdx[at2];
          }
            
          if((idx1>=0)&&(idx2>=0)) {
            subtract3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2,dir);
            dist=length3f(dir);
            if(dist>R_SMALL4) 
              scale3f(dir,1.0/dist,dir);
            if(dist<cutoff) {
              if(mode==1) { /* coarse hydrogen bonding assessment */
                flag=false;
                if(ObjectMoleculeGetAvgHBondVector(obj1,at1,state1,v1)>0.3)
                  if(dot_product3f(v1,dir)<-angle_cutoff) 
                    flag=true;
                if(ObjectMoleculeGetAvgHBondVector(obj2,at2,state2,v2)>0.3)
                  if(dot_product3f(v2,dir)>angle_cutoff)
                    flag=true;
              } else 
                flag=true;

              if(flag) {
                VLACheck((*objVLA),ObjectMolecule*,dist_cnt+1);
                VLACheck((*indexVLA),int,dist_cnt+1);
                (*objVLA)[dist_cnt]=obj1;
                (*indexVLA)[dist_cnt]=at1;
                dist_cnt++;
                (*objVLA)[dist_cnt]=obj2;
                (*indexVLA)[dist_cnt]=at2;              
                dist_cnt++;
              }
            }
          }
        }
      }
    }
  }
  
  CGOStop(DebugCGO);
  
  VLASize((*objVLA),ObjectMolecule*,dist_cnt);
  VLASize((*indexVLA),int,dist_cnt);
  VLAFreeP(vla);
  dist_cnt = dist_cnt / 2;
  return(dist_cnt);
}

/*========================================================================*/
int  SelectorCreateAlignments(int *pair,int sele1,int *vla1,int sele2,
                              int *vla2,char *name1,char *name2,int identical)
{
  SelectorType *I=&Selector;
  int *flag1=NULL,*flag2=NULL;
  int *p;
  int a,i,np;
  int cnt; 
  int mod1,mod2; /* model indexes */
  int at1,at2,at1a,at2a; /* atoms indexes */
  int vi1,vi2; /* vla indexes */
  int index1,index2; /* indices in the selection array */
  AtomInfoType *ai1,*ai2,*ai1a,*ai2a; /* atom information pointers */
  ObjectMolecule *obj1,*obj2;
  int cmp;
  PRINTFD(FB_Selector) 
    " SelectorCreateAlignments-DEBUG: entry.\n"
    ENDFD
  cnt = 0;
  np = VLAGetSize(pair)/2;
  if(np) {
    SelectorUpdateTable();
    flag1=Alloc(int,I->NAtom);
    flag2=Alloc(int,I->NAtom);
    for(a=0;a<I->NAtom;a++)
      {
        flag1[a]=0;
        flag2[a]=0;
      }

    /* we need to create two selection arrays: for the matched 
     * atoms in the original selections */
    p = pair;
    for(i=0;i<np;i++) { /* iterate through all pairs of matched residues */
      vi1 = *(p++);
      vi2 = *(p++);

      /* find positions in the selection arrays */

      mod1 = vla1[vi1*3];
      at1 = vla1[vi1*3+1];
      
      mod2 = vla2[vi2*3];
      at2 = vla2[vi2*3+1];

      PRINTFD(FB_Selector) 
        " S.C.A.-DEBUG: mod1 %d at1 %d mod2 %d at2 %d\n",mod1,at1,mod2,at2
        ENDFD

      obj1 = I->Obj[mod1];
      obj2 = I->Obj[mod2];

      ai1 = obj1->AtomInfo+at1;
      ai2 = obj2->AtomInfo+at2;
      at1a = at1;
      at2a = at2;
      ai1a = obj1->AtomInfo+at1a;
      ai2a = obj2->AtomInfo+at2a;
      while(1) { /* match up all atoms in each residue */
        cmp = AtomInfoNameOrder(ai1a,ai2a);
        if(cmp==0) { /* atoms match */
          index1 = obj1->SeleBase + at1a;
          index2 = obj2->SeleBase + at2a;
        PRINTFD(FB_Selector) 
          " S.C.A.-DEBUG: compare %s %s %d\n",ai1a->name,ai2a->name,cmp
          ENDFD

        PRINTFD(FB_Selector) 
          " S.C.A.-DEBUG: entry %d %d\n",
          ai1a->selEntry,ai2a->selEntry
          ENDFD

          if(SelectorIsMember(ai1a->selEntry,sele1)&&
             SelectorIsMember(ai2a->selEntry,sele2)) {
            if((!identical)||(strcmp(ai1a->resn,ai2a->resn)==0)) {
              flag1[index1] = true;
              flag2[index2] = true; 
              cnt++;
            }
          }
          at1a++;
          at2a++;
        } else if(cmp<0) { /* 1 is before 2 */
          at1a++;
        } else if(cmp>0) { /* 1 is after 2 */
          at2a++;
        }
        if(at1a>=obj1->NAtom) /* make sure we're still in the same residue */
          break;
        if(at2a>=obj2->NAtom)
          break;
        ai1a = obj1->AtomInfo+at1a;
        ai2a = obj2->AtomInfo+at2a;
        if(!AtomInfoSameResidue(ai1a,ai1))
          break;
        if(!AtomInfoSameResidue(ai2a,ai2))
          break;
      }
    }
    if(cnt) {
      SelectorEmbedSelection(flag1,name1,NULL);
      SelectorEmbedSelection(flag2,name2,NULL);
    }
    FreeP(flag1);
    FreeP(flag2);
  }
  PRINTFD(FB_Selector) 
    " SelectorCreateAlignments-DEBUG: exit, cnt = %d.\n",cnt
    ENDFD

  return cnt;
}
/*========================================================================*/
int *SelectorGetResidueVLA(int sele)
{
  /* returns a VLA containing atom indices followed by residue integers
   (residue names packed as characters into integers)
   The indices are the first and last residue in the selection...
  */
  SelectorType *I=&Selector;
  int *result = NULL,*r;
  int a;
  int c;
  AtomInfoType *ai1 = NULL,*ai2;
  int at1=0,at2;
  unsigned int rcode;
  ResName rn;
  int mod1=0;
  ObjectMolecule *obj;

  SelectorUpdateTable();

  result = VLAlloc(int,I->NAtom*3);

  r = result;
  PRINTFD(FB_Selector)
    " SelectorGetResidueVLA-DEBUG: entry, sele = %d\n",sele
    ENDFD;

  if(I->NAtom) {
    for(a=0;a<I->NAtom;a++)
      {
        fflush(stdout);
        obj=I->Obj[I->Table[a].model];
        at2=I->Table[a].atom;
        if(!ai1) {
          mod1 = I->Table[a].model;
          at1 = at2;
          ai1=obj->AtomInfo+at1;
        }
        if(SelectorIsMember(obj->AtomInfo[at2].selEntry,sele)) {
          ai2=obj->AtomInfo+at2;
          if(!AtomInfoSameResidue(ai1,ai2)) {
            if(ai1) {
              *(r++)=mod1;
              *(r++)=at1;
              for(c=0;c<sizeof(ResName);c++)
                rn[c]=0;
              strcpy(rn,ai1->resn); /* store residue code as a number */
              rcode = 0;
              for(c=0;c<3;c++) {
                rcode = (rcode<<8) | rn[c];
              }
              *(r++) = rcode;
              at1 = at2;
              ai1 = ai2;
              mod1 = I->Table[a].model;
            }
          }
        }
      }
    if(ai1) { /* handle last residue */
      *(r++)=mod1;
      *(r++)=at1;
      for(c=0;c<sizeof(ResName);c++)
        rn[c]=0;
      strcpy(rn,ai1->resn); /* store residue code as a number */
      rcode = 0;
      for(c=0;c<3;c++) {
        rcode = (rcode<<8) | rn[c];
      }
      *(r++) = rcode;
    }
  }
  if(result) {
    VLASize(result,int,(r-result));
  }
  PRINTFD(FB_Selector)
    " SelectorGetResidueVLA-DEBUG: exit, result = %p, size = %d\n",result,VLAGetSize(result)
    ENDFD;
  
  return(result);
}
/*========================================================================*/
int *SelectorGetIndexVLA(int sele) /* assumes updated tables */
{
  SelectorType *I=&Selector;
  int a,c=0;
  int *result = NULL;
  ObjectMolecule *obj;
  int at1;

  result = VLAlloc(int,(I->NAtom/10)+1);
  for(a=0;a<I->NAtom;a++) {
    obj=I->Obj[I->Table[a].model];
    at1=I->Table[a].atom;
    if(SelectorIsMember(obj->AtomInfo[at1].selEntry,sele)) {
      VLACheck(result,int,c);
      result[c++]=a;
    }
  }
  VLASize(result,int,c);
  return(result);
}
/*========================================================================*/
void SelectorUpdateObjectSele(ObjectMolecule *obj)
{
  if(obj->Obj.Name[0]) {
    SelectorDelete(obj->Obj.Name);  
    SelectorCreate(obj->Obj.Name,NULL,obj,true,NULL); /* create a selection with same name */ 
  }
}

/*========================================================================*/
void SelectorLogSele(char *name)
{
  SelectorType *I=&Selector;
  int a;
  OrthoLineType line,buf1;
  int cnt=-1;
  int first = 1;
  int append=0;
  ObjectMolecule *obj;
  int at1;
  int sele;
  int logging;
  int robust;
  logging = (int)SettingGet(cSetting_logging);
  robust = (int)SettingGet(cSetting_robust_logs);
  if(logging) {
    sele = SelectorIndexByName(name);
    if(sele>=0) {
      SelectorUpdateTable();
      for(a=0;a<I->NAtom;a++)
        {
          obj=I->Obj[I->Table[a].model];
          at1=I->Table[a].atom;
          if(SelectorIsMember(obj->AtomInfo[at1].selEntry,sele)) {
            
            if(cnt<0) {
              if(first) {
                switch(logging) {
                case cPLog_pml:
                  sprintf(line,"_ select %s,(",name);
                  break;
                case cPLog_pym:
                  sprintf(line,"cmd.select(\"%s\",\"(",name);
                  break;
                }
                append=0;
                cnt=0;
                first=0;
              } else {
                switch(logging) {
                case cPLog_pml:
                  sprintf(line,"_ select %s,(%s",name,name);
                  break;
                case cPLog_pym:
                  sprintf(line,"cmd.select(\"%s\",\"(%s",name,name);
                  break;
                }
                append=1;
                cnt=0;
              }
            }
            if(append) 
              strcat(line,"|");
            if(robust) 
              ObjectMoleculeGetAtomSeleFast(obj,at1,buf1);
            else 
              sprintf(buf1,"%s`%d",obj->Obj.Name,at1+1);
            strcat(line,buf1);
            append=1;
            cnt++;
            if(strlen(line)>(sizeof(OrthoLineType)/2)) {
              switch(logging) {
              case cPLog_pml:
                strcat(line,")\n");
                break;
              case cPLog_pym:
                strcat(line,")\")\n");
                break;
              }
              PLog(line,cPLog_no_flush);
              cnt=-1;
            }
          }
        }
      if(cnt>0) {
        switch(logging) {
        case cPLog_pml:
          strcat(line,")\n");
          break;
        case cPLog_pym:
          strcat(line,")\")\n");
          break;
        }
        PLog(line,cPLog_no_flush);
        PLogFlush();
      }
    }
  }
}
/*========================================================================*/
int SelectorIsMember(int s,int sele)
{
  SelectorType *I=&Selector;
  MemberType *member,*mem;
  member=I->Member;
  while(s) 
    {
      mem = member+s;
      if(mem->selection==sele)
        break;
      s = mem->next;
    }
  return(s);
}
/*========================================================================*/
ObjectMolecule *SelectorGetSingleObjectMolecule(int sele)
{
  int a;
  ObjectMolecule *result = NULL;
  ObjectMolecule *obj;
  SelectorType *I=&Selector;
  int at1;
  SelectorUpdateTable();

  for(a=0;a<I->NAtom;a++)
    {
      obj=I->Obj[I->Table[a].model];
      at1=I->Table[a].atom;
      if(SelectorIsMember(obj->AtomInfo[at1].selEntry,sele)) {
        if(result) {
          if(result!=obj) {
            result=NULL;
            break;
          }
        } else {
          result=obj;
        }
      }
    }
  return(result);
}
/*========================================================================*/
int SelectorGetSingleAtomObjectIndex(int sele,ObjectMolecule **in_obj,int *index)
{
  int found_it = false;
  int a;
  CObject *o = NULL;
  void *hidden = NULL;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  int s;
  while(ExecutiveIterateObject(&o,&hidden))
	 {
		if(o->type==cObjectMolecule)
		  {
			 obj=(ObjectMolecule*)o;
          ai=obj->AtomInfo;
			 for(a=0;a<obj->NAtom;a++)
				{
              s=ai[a].selEntry;
              if(SelectorIsMember(s,sele)) {
                if(found_it){
                  return false; /* ADD'L EXIT POINT */
                } else {
                  found_it = true;
                  (*in_obj)=obj;
                  (*index)=a;
                }
              }
            }
        }
    }
  return(found_it);
}

/*========================================================================*/
int SelectorGetSingleAtomVertex(int sele,int state,float *v)
{
  ObjectMolecule *obj;
  int index;
  int found_it = false;
  if(SelectorGetSingleAtomObjectIndex(sele,&obj,&index))
    found_it = ObjectMoleculeGetAtomVertex(obj,state,index,v);
  return(found_it);
}
/*========================================================================*/
void SelectorDeletePrefixSet(char *pref)
{
  int a;
  SelectorType *I=&Selector;
  WordType name_copy; 

  while(1) {
    a = WordIndex(I->Name,pref,strlen(pref),false);
    if(a>=0) {
      strcpy(name_copy,I->Name[a]);
      ExecutiveDelete(name_copy); /* import to use a copy, otherwise 
                                   * you'll delete all objects  */
    } else
      break;
  }
}
/*========================================================================*/
int SelectorWalkTree(int *atom,int *comp,int *toDo,int **stk,
                     int stkDepth,ObjectMolecule *obj,int sele1,int sele2)
{
  int s;
  int c = 0;
  int a,a1;
  int seleFlag;
  AtomInfoType *ai;

  while(stkDepth) { /* this will explore a tree */
    stkDepth--;
    a=(*stk)[stkDepth];
    toDo[a]=0;
    seleFlag=false;
    ai=obj->AtomInfo+a;
    s=ai->selEntry;
    seleFlag = SelectorIsMember(s,sele1);
    if(!seleFlag)
      seleFlag = SelectorIsMember(s,sele2);      
    if(!seleFlag) {
      if(!(ai->protected==1)) { /* if not explicitly protected...*/
        atom[a]=1; /* mark this atom into the selection */
        comp[a]=1;
      }
      s=obj->Neighbor[a]; /* add neighbors onto the stack */
      s++; /* skip count */
      while(1) {
        a1 = obj->Neighbor[s];
        if(a1>=0) {
          if(toDo[a1]) {
            VLACheck((*stk),int,stkDepth);
            (*stk)[stkDepth]=a1;
            stkDepth++;
          }
        } else 
          break;
        s+=2;
      }
      c++;
    }
  }
  return (c);
}
/*========================================================================*/
int SelectorSubdivideObject(char *pref,ObjectMolecule *obj,int sele1,int sele2,
                            char *fragPref,char *compName)
{
  int a,a0,a1,a2;
  int *atom=NULL;
  int *toDo=NULL;
  int *comp=NULL;
  int nAtom;
  int nFrag = 0;
  int *p1;
  int *stk=NULL;
  int stkDepth;
  int c,s,n;
  int cycFlag=false;
  WordType name;

  PRINTFD(FB_Selector)
    " SelectorSubdivideObject: entered...\n"
    ENDFD;
  SelectorDeletePrefixSet(pref);
  SelectorDeletePrefixSet(fragPref);
  /* delete any existing matches */
  if(obj) {
    ObjectMoleculeUpdateNeighbors(obj);
    SelectorUpdateTableSingleObject(obj);
    nAtom=obj->NAtom;
    if(nAtom) {
      comp = Alloc(int,nAtom);
      p1=comp; /* first atom */
      for(a=0;a<nAtom;a++) 
        *(p1++)=0;
      atom = Alloc(int,nAtom);
      toDo = Alloc(int,nAtom);
      stk=VLAlloc(int,100);
      p1=toDo;
      for(a=0;a<nAtom;a++) 
        *(p1++)=true;
      
      if((sele1>=0)&&(sele2>=0)) { /* bond mode */
        a0 = ObjectMoleculeGetAtomIndex(obj,sele1);
        if(a0>=0) {
          stkDepth=0;
          s=obj->Neighbor[a0]; /* add neighbors onto the stack */
          s++; /* skip count */
          while(1) {
            a1 = obj->Neighbor[s];
            if(a1>=0) {
              if(toDo[a1]) {
                VLACheck(stk,int,stkDepth);
                stk[stkDepth]=a1;
                stkDepth++;
              }
            } else 
              break;
            s+=2;
          }
          p1=atom; /* first atom */
          for(a=0;a<nAtom;a++) 
            *(p1++)=0;
          atom[a0] = 1; /* create selection for this atom alone as fragment base atom */
          comp[a0] = 1;
          sprintf(name,"%s%1d",fragPref,nFrag+1);
          SelectorEmbedSelection(atom,name,NULL);
          c = SelectorWalkTree(atom,comp,toDo,&stk,stkDepth,obj,sele1,sele2) + 1;
          sprintf(name,"%s%1d",pref,nFrag+1);

          /* check for cyclic situation */
          cycFlag=false;
          a2 = ObjectMoleculeGetAtomIndex(obj,sele2);
          if(a2>=0) {
            stkDepth=0;
            s=obj->Neighbor[a2]; /* add neighbors onto the stack */
            s++; /* skip count */
            while(1) {
              a1 = obj->Neighbor[s];
              if (a1<0) break;
              if((a1>=0)&&(a1!=a0)) {
                if(!toDo[a1]) {
                  cycFlag=true; /* we have a cycle...*/
                  break;
                }
              }
              s+=2;
            }
          }
          if(cycFlag) { /* cyclic situation is a bit complex...*/

            a0 = ObjectMoleculeGetAtomIndex(obj,sele2);
            if(a0>=0) {
              stkDepth=0;
              s=obj->Neighbor[a0]; /* add neighbors onto the stack */
              s++; /* skip count */
              while(1) {
                a1 = obj->Neighbor[s];
                if(a1>=0) {
                  if(toDo[a1]) {
                    VLACheck(stk,int,stkDepth);
                    stk[stkDepth]=a1;
                    stkDepth++;
                  }
                } else 
                  break;
                s+=2;
              }
              atom[a0] = 1; 
              comp[a0] = 1;
              c = SelectorWalkTree(atom,comp,toDo,&stk,stkDepth,obj,sele1,sele2) + 1;
            }
          }
          SelectorEmbedSelection(atom,name,NULL);
          nFrag++;
        }
        
        if(!cycFlag) {
          a0 = ObjectMoleculeGetAtomIndex(obj,sele2);
          if(a0>=0) {
            stkDepth=0;
            s=obj->Neighbor[a0]; /* add neighbors onto the stack */
            s++; /* skip count */
            while(1) {
              a1 = obj->Neighbor[s];
              if(a1>=0) {
                if(toDo[a1]) {
                  VLACheck(stk,int,stkDepth);
                  stk[stkDepth]=a1;
                  stkDepth++;
                }
              } else 
                break;
              s+=2;
            }
          
            p1=atom; /* second atom */
            for(a=0;a<nAtom;a++) 
              *(p1++)=0;
            atom[a0] = 1; /* create selection for this atom alone as fragment base atom */
            comp[a0] = 1;
            sprintf(name,"%s%1d",fragPref,nFrag+1);
            SelectorEmbedSelection(atom,name,NULL);
            c = SelectorWalkTree(atom,comp,toDo,&stk,stkDepth,obj,sele1,sele2) + 1;
            sprintf(name,"%s%1d",pref,nFrag+1);
            SelectorEmbedSelection(atom,name,NULL);
            nFrag++;
          }

        }
        
      } else if(sele1>=0) { /* atom mode */
        a0 = ObjectMoleculeGetAtomIndex(obj,sele1);
        comp[a0]=1;
        n=obj->Neighbor[a0];
        n++; /* skip count */
        while(1) {
          a1 = obj->Neighbor[n];
          if(a1<0) break;
          if(toDo[a1]) {
            stkDepth=1;
            stk[0] = a1;
            p1=atom; /* first atom */
            for(a=0;a<nAtom;a++) 
              *(p1++)=0;
            atom[a1] = 1; /* create selection for this atom alone as fragment base atom */
            comp[a1] = 1;
            sprintf(name,"%s%1d",fragPref,nFrag+1);
            SelectorEmbedSelection(atom,name,NULL);
            atom[a1] = 0;
            c = SelectorWalkTree(atom,comp,toDo,&stk,stkDepth,obj,sele1,-1);
            if(c) {
              sprintf(name,"%s%1d",pref,nFrag+1);
              SelectorEmbedSelection(atom,name,NULL);
              nFrag++;
            }
          }
          n+=2;
        }
      }
      if(nFrag) {
        SelectorEmbedSelection(comp,compName,NULL);        
      }
      FreeP(toDo);
      FreeP(atom);
      FreeP(comp);
      VLAFreeP(stk);
      SelectorClean();
    }
  }
  PRINTFD(FB_Selector)
    " SelectorSubdivideObject: leaving...nFrag %d\n",nFrag
    ENDFD;

  return(nFrag);
}
/*========================================================================*/
static int BondInOrder(BondType *a,int b1,int b2)
{
  return(BondCompare(a+b1,a+b2)<=0);
}
/*========================================================================*/
static int BondCompare(BondType *a,BondType *b)
{
  int result;
  if(a->index[0]==b->index[0]) {
	if(a->index[1]==b->index[1]) {
	  result=0;
	} else if(a->index[1]>b->index[1]) {
	  result=1;
	} else {
	  result=-1;
	}
  } else if(a->index[0]>b->index[0]) {
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
    if(SelectorIsMember(s,sele))
      if(result<obj->NCSet) result=obj->NCSet;
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
float SelectorSumVDWOverlap(int sele1,int state1,int sele2,int state2,float adjust)
{
  SelectorType *I=&Selector;
  int *vla=NULL;
  int c;
  float result=0.0;
  float sumVDW=0.0,dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at1,at2;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj1,*obj2;
  int idx1,idx2;
  int a;

  if(state1<0) state1=0;
  if(state2<0) state2=0;

  SelectorUpdateTable();
  c=SelectorGetInterstateVLA(sele1,state1,sele2,state2,2*MAX_VDW+adjust,&vla);
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
       
        idx1=cs1->AtmToIdx[at1]; /* these are also pre-validated */
        idx2=cs2->AtmToIdx[at2];
        
        sumVDW=ai1->vdw+ai2->vdw+adjust;
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
    if(SelectorIsMember(s,sele1))
      {
        if(state1<obj->NCSet) 
          cs=obj->CSet[state1];
        else
          cs=NULL;
        if(cs) {
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {
            copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
            I->Flag1[a]=true;
            n1++;
          }
        }
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
        if(SelectorIsMember(s,sele2))
          {
            if(state2<obj->NCSet) 
              cs=obj->CSet[state2];
            else
              cs=NULL;
            if(cs) {
              if(obj->DiscreteFlag) {
                if(cs==obj->DiscreteCSet[at])
                  idx=obj->DiscreteAtmToIdx[at];
                else
                  idx=-1;
              } else 
                idx=cs->AtmToIdx[at];
              if(idx>=0) {
                v2 = cs->Coord+(3*idx);
                if(MapExclLocus(map,v2,&h,&k,&l)) {
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
                }
                n2++;
              }
            }
          }
		}
		MapFree(map);
	 }
  }
  return(c);
}


/*========================================================================*/
int SelectorMapMaskVDW(int sele1,ObjectMapState *oMap,float buffer)
{
  SelectorType *I=&Selector;
  MapType *map;
  float *v2;
  int n1,n2;
  int a,b,c,i,j,h,k,l;
  int at;
  int s,idx;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  CoordSet *cs;
  int state1;

  state1 = SceneGetState();
  c=0;
  n1=0;
  SelectorUpdateTable();

  for(a=0;a<I->NAtom;a++) {
	 I->Flag1[a]=false;
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(s,sele1))
      {
        if(state1<obj->NCSet) 
          cs=obj->CSet[state1];
        else
          cs=NULL;
        if(cs) {
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {
            copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
            I->Flag1[a]=true;
            n1++;
          }
        }
      }
  }
  /* now create and apply voxel map */
  c=0;
  if(n1) {
	 n2=0;
	 map=MapNewFlagged(-(buffer+MAX_VDW),I->Vertex,I->NAtom,NULL,I->Flag1);
	 if(map) {
		MapSetupExpress(map);
      
      for(a=oMap->Min[0];a<oMap->Max[0];a++) {      
        for(b=oMap->Min[1];b<oMap->Max[1];b++) {      
          for(c=oMap->Min[2];c<oMap->Max[2];c++) {      
            F3(oMap->Field->data,a,b,c)=0.0;            

            v2 = F4Ptr(oMap->Field->points,a,b,c,0);
            
            if(MapExclLocus(map,v2,&h,&k,&l)) {
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                while(j>=0) {
                  ai = I->Obj[I->Table[j].model]->AtomInfo+I->Table[j].atom;
                  if(within3f(I->Vertex+3*j,v2,ai->vdw+buffer)) {
                    F3(oMap->Field->data,a,b,c)=1.0;
                  }
                  j=map->EList[i++];
                }
              }
            }
          }
        }
		}
      oMap->Active=true;
		MapFree(map);
	 }
  }
  return(c);
}


/*========================================================================*/
int SelectorMapGaussian(int sele1,ObjectMapState *oMap,float buffer)
{
  SelectorType *I=&Selector;
  MapType *map;
  float *v2;
  int n1,n2;
  int a,b,c,i,j,h,k,l;
  int at;
  int s,idx;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  CoordSet *cs;
  int state1;
  float *point=NULL,*fp;
  int *sfidx=NULL,*ip;
  int prot;
  float d,e_val;
  double sum,sumsq;
  float mean,stdev;  
  float sf[256][9],*sfp;

  for(a=0;a<256;a++) {
    sf[256][0]=-1.0;
  }

  sf[cAN_H][0] =  0.493002;
  sf[cAN_H][1] = 10.510900;
  sf[cAN_H][2] =  0.322912;
  sf[cAN_H][3] = 26.125700;
  sf[cAN_H][4] =  0.140191;
  sf[cAN_H][5] =  3.142360;
  sf[cAN_H][6] =  0.040810;
  sf[cAN_H][7] = 57.799698;
  sf[cAN_H][8] =  0.003038;

  sf[cAN_C][0] =  2.310000;
  sf[cAN_C][1] = 20.843899;
  sf[cAN_C][2] =  1.020000;
  sf[cAN_C][3] = 10.207500;
  sf[cAN_C][4] =  1.588600;
  sf[cAN_C][5] =  0.568700;
  sf[cAN_C][6] =  0.865000;
  sf[cAN_C][7] = 51.651199;
  sf[cAN_C][8] =  0.215600;

  sf[cAN_O][0] =  3.048500;
  sf[cAN_O][1] =  13.277100;
  sf[cAN_O][2] =  2.286800;
  sf[cAN_O][3] =  5.701100;
  sf[cAN_O][4] =  1.546300;
  sf[cAN_O][5] =  0.323900;
  sf[cAN_O][6] =  0.867000;
  sf[cAN_O][7] =  32.908897;
  sf[cAN_O][8] =  0.250800;

  sf[cAN_N][0] = 12.212600;
  sf[cAN_N][1] =  0.005700;
  sf[cAN_N][2] =  3.132200;
  sf[cAN_N][3] =  9.893300;
  sf[cAN_N][4] =  2.012500;
  sf[cAN_N][5] = 28.997499;
  sf[cAN_N][6] =  1.166300;
  sf[cAN_N][7] =  0.582600;
  sf[cAN_N][8] = -11.528999;
  
  sf[cAN_S][0] =  6.905300;
  sf[cAN_S][1] =  1.467900;
  sf[cAN_S][2] =  5.203400;
  sf[cAN_S][3] = 22.215099;
  sf[cAN_S][4] =  1.437900;
  sf[cAN_S][5] =  0.253600;
  sf[cAN_S][6] =  1.586300;
  sf[cAN_S][7] = 56.172001;
  sf[cAN_S][8] =  0.866900;

  sf[cAN_Cl][0] = 11.460400;
  sf[cAN_Cl][1] =  0.010400; 
  sf[cAN_Cl][2] =  7.196400;
  sf[cAN_Cl][3] =  1.166200;
  sf[cAN_Cl][4] =  6.255600;
  sf[cAN_Cl][5] = 18.519400;
  sf[cAN_Cl][6] =  1.645500;
  sf[cAN_Cl][7] = 47.778400;
  sf[cAN_Cl][8] =  0.866900;

  sf[cAN_Br][0] = 17.178900;
  sf[cAN_Br][1] =  2.172300;
  sf[cAN_Br][2] =  5.235800;
  sf[cAN_Br][3] = 16.579599;
  sf[cAN_Br][4] =  5.637700;
  sf[cAN_Br][5] =  0.260900;
  sf[cAN_Br][6] =  3.985100;
  sf[cAN_Br][7] = 41.432800;
  sf[cAN_Br][8] =  2.955700;

  sf[cAN_I][0] = 20.147200;
  sf[cAN_I][1] = 4.347000;
  sf[cAN_I][2] = 18.994900;
  sf[cAN_I][3] = 0.381400;
  sf[cAN_I][4] = 7.513800;
  sf[cAN_I][5] = 27.765999;
  sf[cAN_I][6] = 2.273500;
  sf[cAN_I][7] = 66.877602;
  sf[cAN_I][8] = 4.071200;
  
  sf[cAN_F][0] = 3.539200;
  sf[cAN_F][1] = 10.282499;
  sf[cAN_F][2] = 2.641200;
  sf[cAN_F][3] = 4.294400;
  sf[cAN_F][4] = 1.517000;
  sf[cAN_F][5] = 0.261500;
  sf[cAN_F][6] = 1.024300;
  sf[cAN_F][7] = 26.147600;
  sf[cAN_F][8] = 0.277600;
  
  sf[cAN_K][0] =   8.218599;
  sf[cAN_K][1] =  12.794900;
  sf[cAN_K][2] =   7.439800;
  sf[cAN_K][3] =   0.774800;
  sf[cAN_K][4] =   1.051900;
  sf[cAN_K][5] = 213.186996;
  sf[cAN_K][6] =    0.865900;
  sf[cAN_K][7] =   41.684097;
  sf[cAN_K][8] =    1.422800;
  
  sf[cAN_Mg][0] = 5.420400;
  sf[cAN_Mg][1] = 2.827500;
  sf[cAN_Mg][2] = 2.173500;
  sf[cAN_Mg][3] = 79.261101;
  sf[cAN_Mg][4] =  1.226900;
  sf[cAN_Mg][5] = 0.380800;
  sf[cAN_Mg][6] =  2.307300;
  sf[cAN_Mg][7] = 7.193700;
  sf[cAN_Mg][8] = 0.858400;

  sf[cAN_Na][0] = 4.762600;
  sf[cAN_Na][1] = 3.285000;
  sf[cAN_Na][2] = 3.173600;
  sf[cAN_Na][3] = 8.842199;
  sf[cAN_Na][4] = 1.267400;
  sf[cAN_Na][5] = 0.313600;
  sf[cAN_Na][6] = 1.112800;
  sf[cAN_Na][7] = 129.423996;
  sf[cAN_Na][8] = 0.676000;

  sf[cAN_P][0] = 6.434500;
  sf[cAN_P][1] = 1.906700;
  sf[cAN_P][2] = 4.179100;
  sf[cAN_P][3] = 27.157000;
  sf[cAN_P][4] =  1.780000;
  sf[cAN_P][5] =  0.526000;
  sf[cAN_P][6] =  1.490800;
  sf[cAN_P][7] = 68.164497;
  sf[cAN_P][8] = 1.114900;
  
  sf[cAN_Zn][0] = 14.074300;
  sf[cAN_Zn][1] = 3.265500;
  sf[cAN_Zn][2] = 7.031800;
  sf[cAN_Zn][3] = 0.233300;
  sf[cAN_Zn][4] = 5.162500;
  sf[cAN_Zn][5] = 10.316299;
  sf[cAN_Zn][6] = 2.410000;
  sf[cAN_Zn][7] = 58.709702;
  sf[cAN_Zn][8] = 1.304100;
  
  sf[cAN_Ca][0] = 8.626600;
  sf[cAN_Ca][1] = 10.442100;
  sf[cAN_Ca][2] = 7.387300;
  sf[cAN_Ca][3] = 0.659900;
  sf[cAN_Ca][4] = 1.589900;
  sf[cAN_Ca][5] = 85.748398;
  sf[cAN_Ca][6] = 1.021100;
  sf[cAN_Ca][7] = 178.436996;
  sf[cAN_Ca][8] = 1.375100;
  
  sf[cAN_Cu][0] = 13.337999;
  sf[cAN_Cu][1] = 3.582800;
  sf[cAN_Cu][2] = 7.167600;
  sf[cAN_Cu][3] = 0.247000;
  sf[cAN_Cu][4] = 5.615800;;
  sf[cAN_Cu][5] = 11.396600;
  sf[cAN_Cu][6] = 1.673500;
  sf[cAN_Cu][7] = 64.812599;
  sf[cAN_Cu][8] = 1.191000;
  
  sf[cAN_Fe][0] = 11.769500;
  sf[cAN_Fe][1] = 4.761100;
  sf[cAN_Fe][2] = 7.357300;
  sf[cAN_Fe][3] = 0.307200;
  sf[cAN_Fe][4] = 3.522200;
  sf[cAN_Fe][5] = 15.353500;
  sf[cAN_Fe][6] = 2.304500;
  sf[cAN_Fe][7] = 76.880501;
  sf[cAN_Fe][8] = 1.036900;
  

  /* 
  sf[cAN_Se][0] = 17.000599;
  sf[cAN_Se][1] = 2.409800;
  sf[cAN_Se][2] = 5.819600;
  sf[cAN_Se][3] = 0.272600;
   */    

  buffer+=MAX_VDW;
  c=0;
  n1=0;
  SelectorUpdateTable();
  n1=0;
  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(s,sele1))
      {
        for(state1=0;state1<obj->NCSet;state1++) {
          if(state1<obj->NCSet) 
            cs=obj->CSet[state1];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              n1++;
            }
          }
        }
      }
  }
  point=Alloc(float,3*n1);
  sfidx=Alloc(int,n1);
  fp=point;
  ip=sfidx;
  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    ai=obj->AtomInfo + at;
    s=ai->selEntry;
    if(SelectorIsMember(s,sele1))
      {

        for(state1=0;state1<obj->NCSet;state1++) {

          if(state1<obj->NCSet) 
            cs=obj->CSet[state1];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              copy3f(cs->Coord+(3*idx),fp);
              fp+=3;
              prot = ai->protons;
              if(sf[prot][0]==-1.0)
                prot=cAN_C;
              *(ip++)=prot;
            }
          }
        }
      }
  }

  PRINTFB(FB_ObjectMap,FB_Details)
    " ObjectMap: Computing Gaussian map for %d atom positions.\n",n1
    ENDFB;
  /* now create and apply voxel map */
  c=0;
  if(n1) {
	 n2=0;
	 map=MapNew(-(buffer),point,n1,NULL);
	 if(map) {
		MapSetupExpress(map);
      sum = 0.0;
      sumsq = 0.0;
      for(a=oMap->Min[0];a<oMap->Max[0];a++) {      
        for(b=oMap->Min[1];b<oMap->Max[1];b++) {      
          for(c=oMap->Min[2];c<oMap->Max[2];c++) {      
            e_val=0.0;
            v2 = F4Ptr(oMap->Field->points,a,b,c,0);

            if(MapExclLocus(map,v2,&h,&k,&l)) {
              i=*(MapEStart(map,h,k,l));
              if(i) {

                j=map->EList[i++];
                while(j>=0) {
                  d = diff3f(point+3*j,v2);
                  if(d<buffer) {

                    d=d*d;
                    if(d<R_SMALL8) d=R_SMALL8;
                    sfp=sf[sfidx[j]];
                    e_val+=
                      (sfp[0]*exp(-sfp[1]*d))
                      +(sfp[2]*exp(-sfp[3]*d))
                      +(sfp[4]*exp(-sfp[5]*d))
                      +(sfp[6]*exp(-sfp[7]*d))
                      +sfp[8];
                  }
                  j=map->EList[i++];
                }
              }
            }
            F3(oMap->Field->data,a,b,c)=e_val;
            sum+=e_val;
            sumsq+=(e_val*e_val);
            n2++;
          }
        }
      }
      mean = sum/n2;
      stdev = sqrt1f((sumsq - (sum*sum/n2))/(n2-1));
      if((int)SettingGet(cSetting_normalize_ccp4_maps)) {

        PRINTFB(FB_ObjectMap,FB_Details)
          " ObjectMap: Normalizing: mean = %8.6f & stdev = %8.6f.\n"
          ,mean,stdev
          ENDFB;
        
        if(stdev<R_SMALL8)
          stdev=R_SMALL8;
        
        for(a=oMap->Min[0];a<oMap->Max[0];a++) {      
          for(b=oMap->Min[1];b<oMap->Max[1];b++) {      
            for(c=oMap->Min[2];c<oMap->Max[2];c++) {      
              fp = F3Ptr(oMap->Field->data,a,b,c);
              
              *fp = (*fp-mean)/stdev;
            }
          }
        }
      } else {
        PRINTFB(FB_ObjectMap,FB_Details)
          " ObjectMap: Not normalizing: mean = %8.6f and stdev = %8.6f.\n",
          mean,stdev
          ENDFB;
        
      }
      oMap->Active=true;
      MapFree(map);
    }
  }
  FreeP(point);
  FreeP(sfidx);
  return(c);
}


/*========================================================================*/
int SelectorMapCoulomb(int sele1,ObjectMapState *oMap,float cutoff)
{
  SelectorType *I=&Selector;
  MapType *map;
  float *v2;
  float distsq;
  int n1,n2;
  int a,b,c,i,j,h,k,l;
  int at;
  int s,idx;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  CoordSet *cs;
  int state1;

  state1 = SceneGetState();
  c=0;
  n1=0;
  SelectorUpdateTable();

  for(a=0;a<I->NAtom;a++) {
	 I->Flag1[a]=false;
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(s,sele1))
      {
        if(state1<obj->NCSet) 
          cs=obj->CSet[state1];
        else
          cs=NULL;
        if(cs) {
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {
            copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
            I->Flag1[a]=true;
            n1++;
          }
        }
      }
  }
  /* now create and apply voxel map */
  c=0;
  if(n1) {
	 n2=0;
	 map=MapNewFlagged(-(cutoff),I->Vertex,I->NAtom,NULL,I->Flag1);
	 if(map) {
		MapSetupExpress(map);
      
      for(a=oMap->Min[0];a<oMap->Max[0];a++) {      
        for(b=oMap->Min[1];b<oMap->Max[1];b++) {      
          for(c=oMap->Min[2];c<oMap->Max[2];c++) {      
            F3(oMap->Field->data,a,b,c)=0.0;            
            v2 = F4Ptr(oMap->Field->points,a,b,c,0);

            if(MapExclLocus(map,v2,&h,&k,&l)) {
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                while(j>=0) {
                  ai = I->Obj[I->Table[j].model]->AtomInfo+I->Table[j].atom;
                  distsq = diffsq3f(I->Vertex+3*j,v2);
                  if(distsq>R_SMALL8) {
                    F3(oMap->Field->data,a,b,c)+=
                      -ai->partialCharge/distsq;
                  }
                  j=map->EList[i++];
                }
              }
            }
          }
        }
		}
      oMap->Active=true;
		MapFree(map);
	 }
  }
  return(c);
}

/*========================================================================*/
int SelectorGetPDB(char **charVLA,int sele,int state,int conectFlag)
{
  SelectorType *I=&Selector;

  int a,b,b1,b2,c,d,s,idx,at,a1,a2;
  BondType *ii1;
  BondType *bond=NULL;
  int nBond=0;
  int cLen =0;
  int newline;
  CoordSet *cs;
  ObjectMolecule *obj;
  AtomInfoType *atInfo,*ai,*last = NULL;
  SelectorUpdateTable();
  c=0;

    /*  if(SettingGet(cSetting_save_pdb_ss)) {
  SSEntry *ss = NULL;
  int n_ss = 0;
  int ss_active = false;

      ss = VLAlloc(SSEntry,100);
      
      for(a=0;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      I->Table[a].index=0;
      obj=I->Obj[I->Table[a].model];
      s = obj->AtomInfo[at].selEntry;
      if(SelectorIsMember(s,sele)) 
        {
          if(state<obj->NCSet) 
            cs=obj->CSet[state];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              ai = obj->AtomInfo+at;

              if(ss_active) {
                
              } else {
                if((at->ss=='H')||(at->ss=='S')) {
                  VLACheck(ss,SSEntry,n_ss);
                  
                }
              }
              
              CoordSetAtomToPDBStrVLA(charVLA,&cLen,ai,
                                      obj->CSet[state]->Coord+(3*idx),c);
              last = ai;
              c++;
            }
          }
        }
    }
  }
    */

  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    I->Table[a].index=0;
    obj=I->Obj[I->Table[a].model];
    s = obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(s,sele)) 
      {
        if(state<obj->NCSet) 
          cs=obj->CSet[state];
        else
          cs=NULL;
        if(cs) {
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {
            ai = obj->AtomInfo+at;
            if(last)
              if(!last->hetatm)
                if(ai->resv!=last->resv)
                  if((abs(ai->resv-last->resv)>1)||(ai->hetatm)) {
                    CoordSetAtomToTERStrVLA(charVLA,&cLen,last,c);
                    c++;
                  }
            I->Table[a].index=c+1; /* NOTE marking with "1" based indexes here */
            CoordSetAtomToPDBStrVLA(charVLA,&cLen,ai,
                                    obj->CSet[state]->Coord+(3*idx),c);
            last = ai;
            c++;
          }
        }
      }
  }
  if(conectFlag) {
    nBond = 0;
    bond = VLAlloc(BondType,1000);
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
          b1=ii1->index[0];
          b2=ii1->index[1];   
          if(obj->DiscreteFlag) {
            if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
              a1=obj->DiscreteAtmToIdx[b1];
              a2=obj->DiscreteAtmToIdx[b2];
            } else {
              a1=-1;
              a2=-1;
            }
          } else {
            a1=cs->AtmToIdx[b1];
            a2=cs->AtmToIdx[b2];
          }
          
          if((a1>=0)&&(a2>=0)&&(atInfo[b1].hetatm||atInfo[b2].hetatm)) {
            b1+=obj->SeleBase;
            b2+=obj->SeleBase;
            if(I->Table[b1].index&&I->Table[b2].index) {
              VLACheck(bond,BondType,2*(nBond+ii1->order+2)); /* ??? */
              b1=I->Table[b1].index;
              b2=I->Table[b2].index;
              for(d=0;d<ii1->order;d++) {
                bond[nBond].index[0] = b1;
                bond[nBond].index[1] = b2;
                nBond++;
                bond[nBond].index[0] = b2;
                bond[nBond].index[1] = b1;
                nBond++;
              }
            }
          }
        ii1++;
        }
      }
    }
    UtilSortInPlace(bond,nBond,sizeof(BondType),(UtilOrderFn*)BondInOrder);
    ii1=bond;
    b1=-1;
	 b2=-1;
    newline = false;
    for(a=0;a<nBond;a++) {
      if(a<(nBond-1)) 
        if((ii1->index[0]==(ii1+1)->index[0])&&(ii1->index[1]==(ii1+1)->index[1])) newline=true;
      if((b1!=ii1->index[0])||((b1==ii1->index[0])&&(b2==ii1->index[1]))||newline) {
        VLACheck((*charVLA),char,cLen+255);
        if(a) cLen+=sprintf((*charVLA)+cLen,"\n");
        cLen+=sprintf((*charVLA)+cLen,"CONECT%5d%5d",
                      ii1->index[0],ii1->index[1]);
        b1=ii1->index[0];
		  b2=ii1->index[1];
        newline=false;
        if(a>0)
          if(((ii1-1)->index[0]==ii1->index[0])&&((ii1-1)->index[1]==ii1->index[1])) newline=true;        
      } else cLen+=sprintf((*charVLA)+cLen,"%5d",
                           ii1->index[1]);
      b2=ii1->index[1];
      ii1++;
    }
    if(cLen) {
      VLACheck((*charVLA),char,cLen+4);
      if(*((*charVLA)+cLen-1)!='\n')
        cLen+=sprintf((*charVLA)+cLen,"\n");
    }
    VLAFree(bond);
  }
  /*
    VLAFreeP(ss); */

  return(cLen);
}
/*========================================================================*/
PyObject *SelectorGetChemPyModel(int sele,int state)
{
  SelectorType *I=&Selector;
  PyObject *model=NULL,*bnd=NULL;
  PyObject *atom_list=NULL,*bond_list=NULL;
  PyObject *tmp;
  PyObject *molecule = NULL;
  int a,b,b1,b2,c,s,idx,at,a1,a2;
  BondType *ii1;
  BondType *bond=NULL;
  int nBond=0;
  int ok =true;
  CoordSet *cs;
  int single_flag = true;
  CoordSet *single_cs = NULL;
  ObjectMolecule *obj;
  AtomInfoType *atInfo,*ai;
  SelectorUpdateTable();

  model = PyObject_CallMethod(P_models,"Indexed","");
  if (!model) 
    ok = ErrMessage("CoordSetAtomToChemPyAtom","can't create model");  
  if(ok) {    
    c=0;
    for(a=0;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      I->Table[a].index=0;
      obj=I->Obj[I->Table[a].model];
      s=obj->AtomInfo[at].selEntry;
      if(SelectorIsMember(s,sele))
        {
          if(state<obj->NCSet) 
            cs=obj->CSet[state];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              I->Table[a].index=c+1; /* NOTE marking with "1" based indexes here */
              c++;
            }
          }
        }
    }
    if(c) {
      atom_list = PyList_New(c);
      PyObject_SetAttrString(model,"atom",atom_list);
      c=0;
      for(a=0;a<I->NAtom;a++) {
        if(I->Table[a].index) {
          at=I->Table[a].atom;
          obj=I->Obj[I->Table[a].model];
          cs=obj->CSet[state]; /* assuming this is valid... */
          if(obj->DiscreteFlag) {
            if(obj->CSet[state]==obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {
            if(single_flag) { /* remember whether all atoms come from a single coordinate set...*/
              if(single_cs) {
                if(single_cs!=cs)
                  single_flag=false;
              } else {
                single_cs = cs;
              }
            }
            ai = obj->AtomInfo+at;
            PyList_SetItem(atom_list,c,
                           CoordSetAtomToChemPyAtom(
                       ai,obj->CSet[state]->Coord+(3*idx),at));
            c = c + 1;
          }
        }
      }
      Py_XDECREF(atom_list);

      if(single_flag&&single_cs) { /* single coordinate set?  then set coordinate set info */
        molecule = PyObject_GetAttrString(model,"molecule");
        if(molecule) {
          if(single_cs->Name[0]) {
            PyObject_SetAttrString(molecule,"title",           /* including name/title */
                                   PyString_FromString(single_cs->Name)); 
          }
        }
        Py_XDECREF(molecule);
      }

      nBond = 0;
      bond = VLAlloc(BondType,1000);
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
            b1=ii1->index[0];
            b2=ii1->index[1];        
            if(obj->DiscreteFlag) {
              if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
                a1=obj->DiscreteAtmToIdx[b1];
                a2=obj->DiscreteAtmToIdx[b2];
              } else {
                a1=-1;
                a2=-1;
              }
            } else {
              a1=cs->AtmToIdx[b1];
              a2=cs->AtmToIdx[b2];
            }

            if((a1>=0)&&(a2>=0)) {
              b1+=obj->SeleBase;
              b2+=obj->SeleBase;
              if(I->Table[b1].index&&I->Table[b2].index) {
                VLACheck(bond,BondType,nBond);
                bond[nBond] = *ii1;
                b1=I->Table[b1].index - 1; /* counteract 1-based */
                b2=I->Table[b2].index - 1; /* indexing from above */
                bond[nBond].index[0] = b1;
                bond[nBond].index[1] = b2;
                nBond++;
              }
            }
            ii1++;
          }
        }
        if(cs) {
          if(c==cs->NIndex) { /* support for experimental spheroids - likely to change */
            if(cs->Spheroid&&cs->SpheroidNormal) {
              tmp = PConvFloatArrayToPyList(cs->Spheroid,cs->NSpheroid);
              PyObject_SetAttrString(model,"spheroid",tmp);
              Py_XDECREF(tmp);          
              tmp = PConvFloatArrayToPyList(cs->SpheroidNormal,cs->NSpheroid*3);
              PyObject_SetAttrString(model,"spheroid_normals",tmp);
              Py_XDECREF(tmp);          
            }
          }
        }

        ii1=bond;
        bond_list = PyList_New(nBond);
        PyObject_SetAttrString(model,"bond",bond_list);
        for(b=0;b<nBond;b++) {
          bnd = PyObject_CallMethod(P_chempy,"Bond","");
          if(bnd) {
            PConvInt2ToPyObjAttr(bnd,"index",ii1->index);
            PConvIntToPyObjAttr(bnd,"order",ii1->order);
            PConvIntToPyObjAttr(bnd,"id",ii1->id);
            PConvIntToPyObjAttr(bnd,"stereo",ii1->stereo);
            PyList_SetItem(bond_list,b,bnd);
          }
          ii1++;
        }
        Py_XDECREF(bond_list);
      }
      VLAFree(bond);
    }
  }
  return(model);
}
/*========================================================================*/
void SelectorUpdateCmd(int sele0,int sele1,int sta0, int sta1)
{
  SelectorType *I=&Selector;
  int a,b;
  int at0,at1;
  int *vla0=NULL;
  int *vla1=NULL;
  int c0,c1;
  int i0=0,i1;
  int cc1;
  ObjectMolecule *obj0,*obj1;
  CoordSet *cs0,*cs1;
  int matched_flag;
  int b_start;
  int ci0,ci1;
  int ccc = 0;


  PRINTFD(FB_Selector)
    " SelectorUpdateCmd-Debug: entered sta0 %d sta1 %d",sta0,sta1
    ENDFD;

  SelectorUpdateTable();

  vla0 = SelectorGetIndexVLA(sele0);
  vla1 = SelectorGetIndexVLA(sele1);

  if(!(vla0&&vla1))
    ErrMessage("Update","no coordinates updated.");
  else {
    c0 = VLAGetSize(vla0);
    c1 = VLAGetSize(vla1);

    b = 0;
    for(a=0;a<c1;a++) { /* iterate over source atoms */
      i1 = vla1[a];
      at1=I->Table[i1].atom;
      obj1=I->Obj[I->Table[i1].model];
      
      b_start = b;
      matched_flag=false;
      while(1) {
        i0 = vla0[b];
        at0=I->Table[i0].atom;
        obj0=I->Obj[I->Table[i0].model];
        if(obj0!=obj1)
          if(AtomInfoMatch(obj1->AtomInfo + at1,
                           obj0->AtomInfo + at0)) {
            matched_flag=true;
            break;
          }
        b++;
        if(b>=c0)
          b = 0;
        if(b==b_start) 
          break;
      }
      if(matched_flag) { /* atom matched, so copy coordinates */
        ccc++;
        for(cc1=0;cc1<obj1->NCSet;cc1++) { /* iterate over all source states */
          if((cc1==sta1)||(sta1<0)) {
            cs1 = obj1->CSet[cc1];
            if(cs1&&(cc1<obj0->NCSet)&&
               ((sta0<0)||(cc1==sta0)|| /* multiple or single state */
                ((sta0>=0)&&(sta1>=0)))) { /* explicit state */
              if((sta0<0)||(sta0>=obj0->NCSet))
                cs0 = obj0->CSet[cc1];
              else
                cs0 = obj0->CSet[sta0];
              if(cs0) {
                ci0 = cs0->AtmToIdx[at0];
                ci1 = cs1->AtmToIdx[at1];
                if((ci0>=0)&&(ci1>=0))
                  copy3f(cs1->Coord + 3*ci1,
                         cs0->Coord + 3*ci0);
              }
            }
          }
        }
      }
    }
    obj0=NULL;
    for(b=0;b<c1;b++) {
      obj1=I->Obj[I->Table[i0].model];
      if(obj0!=obj1) {
        obj0=obj1;
        ObjectMoleculeInvalidate(obj0,cRepAll,cRepInvCoord);
      }
    }
    PRINTFB(FB_Selector,FB_Actions)
      " Update: coordinates updated for %d atoms.\n",ccc 
      ENDFB
  }
  VLAFreeP(vla0);
  VLAFreeP(vla1);
}
/*========================================================================*/

void SelectorCreateObjectMolecule(int sele,char *name,int target,int source)
{
  SelectorType *I=&Selector;

  int a,b,a1,a2,b1,b2,c,d,s,at;
  BondType *ii1,*bond=NULL;
  int nBond=0;
  int nCSet,nAtom,ts;
  AtomInfoType *atInfo = NULL;
  int isNew,csFlag;
  CoordSet *cs = NULL;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj;
  CObject *ob;
  ObjectMolecule *targ = NULL;
  ObjectMolecule *info_src = NULL;

  ob=ExecutiveFindObjectByName(name);
  if(ob)
    if(ob->type==cObjectMolecule) 
      targ = (ObjectMolecule*)ob;
  if(!targ) {
    isNew=true;
    targ = ObjectMoleculeNew(false);
    targ->Bond = VLAlloc(BondType,1);
  } else {
    isNew=false;
  }
  
  c=0;
  SelectorUpdateTable();
  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    I->Table[a].index=-1;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(s,sele))
      {
        I->Table[a].index=c; /* Mark records as to which atom they'll be */
        c++;
        if(!info_src) info_src = obj;
      }
  }
  if(isNew&&info_src) {/* copy symmetry information, etc. */
    targ->Symmetry=SymmetryCopy(info_src->Symmetry);
  }
  nAtom=c;

  nBond = 0;
  bond = VLAlloc(BondType,nAtom*4);
  for(a=0;a<I->NModel;a++) { /* find bonds wholly contained in the selection */
    obj=I->Obj[a];
    ii1=obj->Bond;
    for(b=0;b<obj->NBond;b++) {
      b1=ii1->index[0]+obj->SeleBase;
      b2=ii1->index[1]+obj->SeleBase;
      if((I->Table[b1].index>=0)&&(I->Table[b2].index>=0)) {
        VLACheck(bond,BondType,nBond);
        bond[nBond].index[0]=I->Table[b1].index; /* store what will be the new index */
        bond[nBond].index[1]=I->Table[b2].index;
        bond[nBond].order=ii1->order;
        bond[nBond].stereo=ii1->stereo;
        nBond++;
      }
      ii1++;
    }
  }
  
  atInfo = VLAlloc(AtomInfoType,nAtom); 
  /* copy the atom info records and create new zero-based IDs */
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
  
  ObjectMoleculeMerge(targ,atInfo,cs,false,cAIC_AllMask); /* will free atInfo */
  /* cs->IdxToAtm will now have the reverse mapping from the new subset
     to the new merged molecule */

  ObjectMoleculeExtendIndices(targ);
  ObjectMoleculeUpdateNonbonded(targ);
  
  if(!isNew) { /* recreate selection table */
    SelectorUpdateTable(); 
    
    c=0;
    for(a=0;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      I->Table[a].index=-1;
      obj=I->Obj[I->Table[a].model];
      s=obj->AtomInfo[at].selEntry;
      if(SelectorIsMember(s,sele))
        {
          I->Table[a].index=c; /* Mark records  */
          c++;
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
      csFlag = true;

      /* any selected atoms in this state? */
      /*
        for(a=0;a<I->NAtom;a++)  
        if(I->Table[a].index>=0) {
        at=I->Table[a].atom;
        obj=I->Obj[I->Table[a].model];
        if(d<obj->NCSet) {
        cs1 = obj->CSet[d];
        if(cs1) {
        if(obj->DiscreteFlag) {
        if(cs1==obj->DiscreteCSet[at])
        a1=obj->DiscreteAtmToIdx[at];
        else
        a1=-1;
        } else 
        a1 = cs1->AtmToIdx[at];
        if(a1>=0) {
        csFlag=true;
        break;
        }
        }
        }
        }*/

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
                if(obj->DiscreteFlag) {
                  if(cs1==obj->DiscreteCSet[at])
                    a1=obj->DiscreteAtmToIdx[at];
                  else
                    a1=-1;
                } else 
                  a1 = cs1->AtmToIdx[at]; /* coord index in existing object */
                if(a1>=0) {
                  copy3f(cs1->Coord+a1*3,cs2->Coord+c*3);
                  a2 = cs->IdxToAtm[I->Table[a].index]; /* actual merged atom index */
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
  VLAFreeP(bond); /* null-safe */
  if(cs) cs->fFree(cs);
  if(targ->DiscreteFlag) { /* if the new object is discrete, then eliminate the AtmToIdx array */
    for(d=0;d<targ->NCSet;d++) {
      cs = targ->CSet[d];
      if(cs) {
        if(cs->AtmToIdx) {
          for(a=0;a<cs->NIndex;a++) {
            b = cs->IdxToAtm[a];
            targ->DiscreteAtmToIdx[b] = a;
            targ->DiscreteCSet[b] = cs;
          }
          FreeP(cs->AtmToIdx);
        }
      }
    }
  }
  SceneCountFrames();
  PRINTFB(FB_Selector,FB_Details)
    " Selector: found %d atoms.\n",nAtom 
    ENDFB
    ObjectMoleculeSort(targ);
  if(isNew) {
    ObjectSetName((CObject*)targ,name);
    ExecutiveManageObject((CObject*)targ,true);
  } else {
    ExecutiveUpdateObjectSelection((CObject*)targ);
  }
  SceneChanged();
}


/*========================================================================*/
int SelectorIndexByName(char *sname)
{
 OrthoLineType name;
 SelectorType *I=&Selector;
 int i=-1;
 if(sname) {
   if(sname[0]=='%')
     strcpy(name,&sname[1]);
   else
     strcpy(name,sname);		  
   i = WordIndex(I->Name,name,1,I->IgnoreCase);
   if(i>=0) i = I->ID[i];
 }
 return(i);
}
/*========================================================================*/
void SelectorPurgeMembers(int sele) 
{
  int a=0;
  int s=0;
  int l;
  int nxt;
  CObject *o = NULL;
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
                    nxt = I->Member[s].next;
						  if(I->Member[s].selection==sele)
							 {
								if(l>=0)
								  I->Member[l].next=I->Member[s].next;
								else
								  obj->AtomInfo[a].selEntry=I->Member[s].next;
                        I->Member[s].next = I->FreeMember; 
                        I->FreeMember=s;
							 }
                    l=s;
						  s=nxt;
						}
				  }
			 }
		}
}
/*========================================================================*/
void SelectorDelete(char *sele) 
     /* should (only) be called by Executive or by Selector, unless the
        named selection has never been registered with the executive 
        (i.e., temporary on-the-fly selection) */
{
  SelectorType *I=&Selector;
  int n;
  int index;
  n=WordIndex(I->Name,sele,999,I->IgnoreCase); /* already exist? */
  if(n>=0) /* get rid of existing selection*/
	 {
      index = I->ID[n];
		SelectorPurgeMembers(index);
      I->NActive--;
      strcpy(I->Name[n],I->Name[I->NActive]);
      I->ID[n]=I->ID[I->NActive];
      I->Name[I->NActive][0]=0;
	 }
}
/*========================================================================*/
int SelectorGetTmp(char *input,char *store)
{
  SelectorType *I=&Selector;
  WordType name;
  OrthoLineType buffer;
  int count = 0;
  PRINTFD(FB_Selector)
    " SelectorGetTmp-Debug: entered with '%s'.\n",input
    ENDFD;

  if(input[0]=='(') {
    sprintf(name,"%s%d",cSelectorTmpPrefix,I->TmpCounter++);
	 count = SelectorCreate(name,input,NULL,false,NULL);
	 strcpy(store,name);
  } else {
    if(ExecutiveValidName(input))
      strcpy(store,input);
    else {
      strcpy(buffer,"(");
      strcat(buffer,input);
      strcat(buffer,")");
      sprintf(name,"%s%d",cSelectorTmpPrefix,I->TmpCounter++);      
      count = SelectorCreate(name,buffer,NULL,false,NULL);
      strcpy(store,name);
    }
  }
  PRINTFD(FB_Selector)
    " SelectorGetTmp-Debug: leaving with '%s'.\n",store
    ENDFD;
  return count;
}
/*========================================================================*/
void SelectorFreeTmp(char *name)
{
  if(strncmp(name,cSelectorTmpPrefix,strlen(cSelectorTmpPrefix))==0) {
    ExecutiveDelete(name);
  }
}
/*========================================================================*/
int  SelectorEmbedSelection(int *atom, char *name, ObjectMolecule *obj)
{
  /* either atom or obj should be NULL, not both and not neither */

  SelectorType *I=&Selector;
  int flag;
  int newFlag=true;
  int n,a,m,sele;
  int c=0;
  AtomInfoType *ai;
  n=WordIndex(I->Name,name,999,I->IgnoreCase); /* already exist? */
  if(n>=0) /* get rid of existing selection*/ {
    SelectorDelete(I->Name[n]);
    newFlag = false;
  }
  n=I->NActive;
  VLACheck(I->Name,WordType,n+1);
  VLACheck(I->ID,int,n+1);
  strcpy(I->Name[n],name);
  I->Name[n+1][0]=0;
  sele = I->NSelection++;
  I->ID[n] = sele;
  I->NActive++;
  fflush(stdout);
  for(a=0;a<I->NAtom;a++)
    {
      flag=false;
      if(atom) {
        if(atom[a]) flag=true;
      } else {
        if(I->Obj[I->Table[a].model]==obj) flag=true;
      }
      if(flag)
        {
          ai = I->Obj[I->Table[a].model]->AtomInfo+I->Table[a].atom;

          c++;
          if(I->FreeMember>=0) {
            m=I->FreeMember;
            I->FreeMember=I->Member[m].next;
          } else {
            I->NMember++;
            m=I->NMember;
            VLACheck(I->Member,MemberType,m);
          }
          I->Member[m].selection=sele;
          I->Member[m].next = ai->selEntry;
          ai->selEntry = m;
        }
    }
  if(!obj) {
    if(newFlag)
      ExecutiveManageSelection(name);
    else
      ExecutiveSetControlsOff(name);
  }
  PRINTFD(FB_Selector)
    " Selector: Embedded %s, %d atoms.\n",name,c
    ENDFD;
  return(c);
}
/*========================================================================*/
int *SelectorApplyMultipick(Multipick *mp)
{
  SelectorType *I=&Selector;
  int *result;
  int a,n;
  Pickable *p;
  ObjectMolecule *obj;
  SelectorUpdateTable();  
  result = Alloc(int,I->NAtom);
  n=mp->picked[0].index;
  p=mp->picked+1;
  for(a=0;a<I->NAtom;a++) 
    result[a]=0;
  while(n--) {
    obj=(ObjectMolecule*)p->ptr;
    result[obj->SeleBase+p->index] = true;
    p++;
  }
  return(result);
}
/*========================================================================*/
int SelectorCreate(char *sname,char *sele,ObjectMolecule *obj,int quiet,Multipick *mp) 
{
  SelectorType *I=&Selector;
  int *atom=NULL;
  OrthoLineType name;
  int ok=true;
  int c=0;

  PRINTFD(FB_Selector)
    "SelectorCreate-Debug: entered...\n"
    ENDFD;

  if(sname[0]=='%')
	 strcpy(name,&sname[1]);
  else
	 strcpy(name,sname);		  
  UtilCleanStr(name);
  if(!name[0])
	 {
      PRINTFB(FB_Selector,FB_Errors)
        "Selector-Error: Invalid selection name '%s'.\n",sname
        ENDFB;
		OrthoRestorePrompt();
	 }
  if(ok)
	 {
	   if(sele) {
		 atom=SelectorSelect(sele);
		 if(!atom) ok=false;
	   } else if(obj) { /* optimized full-object selection */
        SelectorUpdateTableSingleObject(obj);
	   } else if(mp) {
        atom=SelectorApplyMultipick(mp);
      } else
        ok=false;
	 }
  if(ok)	c=SelectorEmbedSelection(atom,name,obj);
  FreeP(atom);
  SelectorClean();
  I->NAtom=0;
  if(!quiet) {
    if(name[0]!='_') {
      if(c) {
        PRINTFB(FB_Selector,FB_Results)
          " Selector: selection \"%s\" defined with %d atoms.\n",name,c
          ENDFB;
      } else {
        PRINTFB(FB_Selector,FB_Results)
          " Selector: no atoms selected.\n"
          ENDFB;
      }
    } 
  }
  PRINTFD(FB_Selector)
    " SelectorCreate: \"%s\" created with %d atoms.\n",name,c    
    ENDFD
  return(c);
}
/*========================================================================*/
void SelectorClean(void)
{
  SelectorType *I=&Selector;
  FreeP(I->Table);
  FreeP(I->Obj);
  FreeP(I->Vertex);
  FreeP(I->Flag1);
  FreeP(I->Flag2);
}
/*========================================================================*/
int SelectorUpdateTableSingleObject(ObjectMolecule *obj)
{
  int a=0;
  int c=0;
  int modelCnt;

  SelectorType *I=&Selector;

  PRINTFD(FB_Selector)
    "SelectorUpdateTableSingleObject-Debug: entered..\n"
    ENDFD;
  

  SelectorClean();

  I->NCSet = 0;
  modelCnt=0;
  c+=obj->NAtom;
  if(I->NCSet<obj->NCSet) I->NCSet=obj->NCSet;
  modelCnt++;
  I->Table=Alloc(TableRec,c);
  ErrChkPtr(I->Table);
  I->Obj=Alloc(ObjectMolecule*,modelCnt);
  ErrChkPtr(I->Obj);
  c=0;
  modelCnt=0;
  I->Obj[modelCnt]=NULL;
  I->Obj[modelCnt]=obj;
  obj->SeleBase=c; 
  for(a=0;a<obj->NAtom;a++)
    {
      I->Table[c].model=modelCnt;
      I->Table[c].atom=a;
      c++;
    }
  modelCnt++;
  I->NModel=modelCnt;
  I->NAtom=c;
  I->Flag1=Alloc(int,c);
  ErrChkPtr(I->Flag1);
  I->Flag2=Alloc(int,c);
  ErrChkPtr(I->Flag2);
  I->Vertex=Alloc(float,c*3);
  ErrChkPtr(I->Vertex);

  PRINTFD(FB_Selector)
    "SelectorUpdateTableSingleObject-Debug: leaving...\n"
    ENDFD;

  return(true);
}
/*========================================================================*/
int SelectorUpdateTable(void)
{
  int a=0;
  int c=0;
  int modelCnt;
  CObject *o = NULL;
  void *hidden = NULL;
  ObjectMolecule *obj;

  SelectorType *I=&Selector;
  SelectorClean();
  I->NCSet = 0;
  modelCnt=0;
  while(ExecutiveIterateObject(&o,&hidden))
	 {
		if(o->type==cObjectMolecule)
		  {
			 obj=(ObjectMolecule*)o;
			 c+=obj->NAtom;
			 if(I->NCSet<obj->NCSet) I->NCSet=obj->NCSet;
          modelCnt++;
		  }
	 }
  I->Table=Alloc(TableRec,c);
  ErrChkPtr(I->Table);
  I->Obj=Alloc(ObjectMolecule*,modelCnt);
  ErrChkPtr(I->Obj);
  c=0;
  modelCnt=0;
  while(ExecutiveIterateObject(&o,&hidden))
	 {
		if(o->type==cObjectMolecule)
		  {
          I->Obj[modelCnt]=NULL;
			 obj=(ObjectMolecule*)o;
			 I->Obj[modelCnt]=obj;
          obj->SeleBase=c; /* make note of where this object starts */
			 for(a=0;a<obj->NAtom;a++)
				{
				  I->Table[c].model=modelCnt;
				  I->Table[c].atom=a;
				  c++;
				}
          modelCnt++;
		  }
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
  PRINTFD(FB_Selector)
    "SelectorSelect-DEBUG: sele = '%s'\n",sele
    ENDFD;
  SelectorUpdateTable();
  parsed=SelectorParse(sele);
  if(parsed)
	 {
      if(Feedback(FB_Selector,FB_Debugging)) {
        WordType *a;
        fprintf(stderr,"SelectorSelect-DEBUG: parsed tokens:\n");
        a = parsed;
        while(1) {
          if(!a[0][0]) break;
          fprintf(stderr,"  '%s'\n",(a[0]));
          a++;
        }
        fprintf(stderr,"SelectorSelect-DEBUG: end of tokens.\n");
      }
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
                  if(obj->DiscreteFlag) {
                    if(cs==obj->DiscreteCSet[at])
                      idx=obj->DiscreteAtmToIdx[at];
                    else
                      idx=-1;
                  } else 
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
                            if(obj->DiscreteFlag) {
                              if(cs==obj->DiscreteCSet[at])
                                idx=obj->DiscreteAtmToIdx[at];
                              else
                                idx=-1;
                            } else 
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
	 case SELE_GAP_:
		if(!sscanf(base[2].text,"%f",&dist))
		  ok=ErrMessage("Selector","Invalid distance.");
		if(ok)
		  {
          for(a=0;a<I->NAtom;a++) {
            obj=I->Obj[I->Table[a].model];
            at=I->Table[a].atom;
            I->Table[a].f1 = obj->AtomInfo[at].vdw;
            base[0].sele[a]=true; /* start selected, subtract off */
            c=I->NAtom;
          }
			 for(d=0;d<I->NCSet;d++)
				{
				  n1=0;
				  for(a=0;a<I->NAtom;a++) {
                obj=I->Obj[I->Table[a].model];
                I->Flag1[a]=false;
					 at=I->Table[a].atom;
                if(d<obj->NCSet) 
                  cs=obj->CSet[d];
                else
                  cs=NULL;
					 if(cs) {
                  if(obj->DiscreteFlag) {
                    if(cs==obj->DiscreteCSet[at])
                      idx=obj->DiscreteAtmToIdx[at];
                    else
                      idx=-1;
                  } else 
                    idx=cs->AtmToIdx[at];
						if(idx>=0) {
						  copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
						  I->Flag1[a]=true;
						  n1++;
						}
					 }
				  }
				  if(n1) {
					 map=MapNewFlagged(-(dist+2*MAX_VDW),I->Vertex,I->NAtom,NULL,I->Flag1);
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
                            if(obj->DiscreteFlag) {
                              if(cs==obj->DiscreteCSet[at])
                                idx=obj->DiscreteAtmToIdx[at];
                              else
                                idx=-1;
                            } else 
                              idx=cs->AtmToIdx[at];

									 if(idx>=0) {
										v2 = cs->Coord+(3*idx);
										MapLocus(map,v2,&h,&k,&l);
										i=*(MapEStart(map,h,k,l));
										if(i) {
										  j=map->EList[i++];
										  while(j>=0) {
											 if((base[0].sele[j])&&
												 (!base[1].sele[j])) /*exclude current selection */
												{
												  if(within3f(I->Vertex+3*j,v2,dist+ /* eliminate atoms w/o gap */
                                                  I->Table[a].f1+
                                                  I->Table[j].f1)) {
                                        base[0].sele[j]=false;
                                        c--;
                                      }
												}
                                  else if (base[1].sele[j]) {
                                    base[0].sele[j]=false;
                                    c--;
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
  if(Feedback(FB_Selector,FB_Debugging)) {
	 c=0;
	 for(a=0;a<I->NAtom;a++)
		if(base[0].sele[a]) c++;
	 fprintf(stderr,"SelectorModulate0: %d atoms selected.\n",c);
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
  int state;
  int static_singletons;
  ObjectMolecule *obj,*cur_obj=NULL;
  CoordSet *cs;
  int at_idx;
  base->type=STYP_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);

  switch(base->code)
	 {
	 case SELE_NONz:
		for(a=0;a<I->NAtom;a++)
		  base[0].sele[a]=false;
		break;
	 case SELE_BNDz:
		for(a=0;a<I->NAtom;a++)
        base[0].sele[a]=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].bonded;
		break;
	 case SELE_HETz:
		for(a=0;a<I->NAtom;a++)
        base[0].sele[a]=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hetatm;
		break;
	 case SELE_HYDz:
		for(a=0;a<I->NAtom;a++)
        base[0].sele[a]=I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hydrogen;
		break;
    case SELE_PREz:
      state = SceneGetState();
      static_singletons = (int)SettingGet(cSetting_static_singletons);
      flag=false;
      cs = NULL;
      for(a=0;a<I->NAtom;a++)
        {
          base[0].sele[a]=false;
          obj = I->Obj[I->Table[a].model];
          if(obj!=cur_obj) { /* different object */
            if(state>=obj->NCSet)
              flag=false;
            else if(!obj->CSet[state]) 
              flag=false;
            else {
              cs = obj->CSet[state];
              flag=true; /* valid state */
            }
            if(!flag)
              if(I->NCSet==1)
                if(static_singletons) {
                  cs = obj->CSet[0];
                  flag=true;
                }
            cur_obj = obj;
          }
          if(flag&&cs) {
            at_idx = I->Table[a].atom;
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at_idx]) {
                if(obj->DiscreteAtmToIdx[at_idx]>=0) {
                  base[0].sele[a]=true;
                  c++;
                }
              }
            } else {
            } if(cs->AtmToIdx[at_idx]>=0) {
              base[0].sele[a]=true;
              c++;
            }
          }
        }
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
          if(I->Obj[I->Table[a].model]->Obj.Enabled) {
            vis = I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].visRep ;
            
            for(b=0;b<cRepCnt;b++) {
              if(vis[b]) {
                flag=true;
                break;
              }
            }
          }
          base[0].sele[a]=flag;
          if(flag)
            c++;
		  }
		break;
	 }
  PRINTFD(FB_Selector)
	 " SelectorSelect0: %d atoms selected.\n",c
    ENDFD;

  return(1);
}
/*========================================================================*/
int SelectorSelect1(EvalElem *base)
{
  int a,model,sele,s,at_idx;
  int c=0;
  int ok=true;
  int rmin,rmax,rtest,index,numeric,state;
  int flag;
  char *p;
  char *np;
  SelectorType *I=&Selector;
  ObjectMolecule *obj,*cur_obj = NULL;
  CoordSet *cs=NULL;

  base->type=STYP_LIST;
  base->sele=Alloc(int,I->NAtom);
  PRINTFD(FB_Selector)
    " SelectorSelect1: base: %p sele: %p\n",base,base->sele
  ENDFD;
  ErrChkPtr(base->sele);
  switch(base->code)
	 {
	 case SELE_IDXs:
		if(sscanf(base[1].text,"%i",&index)!=1)		
		  ok=ErrMessage("Selector","Invalid Index.");
		if(ok) {
        index--;
		  for(a=0;a<I->NAtom;a++)
			 {
				if(I->Table[a].atom==index) {
              c++;
				  base[0].sele[a]=true;
            } else {
				  base[0].sele[a]=false;
            }
			 }
		}
		break;
	 case SELE_ID_s:
      WordPrimeCommaMatch(base[1].text);
      for(a=0;a<I->NAtom;a++)
        {
          if(WordMatchCommaInt(base[1].text,
                               I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].id)<0)
            {
              base[0].sele[a]=true;
              c++;
            } else {
              base[0].sele[a]=false;
            }
        }
      break;
	 case SELE_NAMs:
      WordPrimeCommaMatch(base[1].text);
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchCommaExact(base[1].text,
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
	 case SELE_TTYs:
      WordPrimeCommaMatch(base[1].text);
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
                            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].textType,
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
      WordPrimeCommaMatch(base[1].text);
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
      WordPrimeCommaMatch(base[1].text);
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
              I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].segi,I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
		  }
		break;
	 case SELE_CHNs:
      WordPrimeCommaMatch(base[1].text);
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
	 case SELE_SSTs:
      WordPrimeCommaMatch(base[1].text);
		for(a=0;a<I->NAtom;a++)
		  {
			 if(WordMatchComma(base[1].text,
                            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].ssType,
                            I->IgnoreCase)<0)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
        }
		break;
	 case SELE_STAs:
      sscanf(base[1].text,"%d",&state);
      state = state - 1;
      obj = NULL;
      
      if(state<0) {
        for(a=0;a<I->NAtom;a++)
          base[0].sele[a]=false;
      } else {
        for(a=0;a<I->NAtom;a++)
          {
            base[0].sele[a]=false;
            obj = I->Obj[I->Table[a].model];
            if(obj!=cur_obj) { /* different object */
              if(state>=obj->NCSet)
                flag=false;
              else if(!obj->CSet[state]) 
                flag=false;
              else {
                cs = obj->CSet[state];
                flag=true; /* valid state */
              }
              cur_obj = obj;
            }
            if(flag&&cs) {
              at_idx = I->Table[a].atom;
              if(obj->DiscreteFlag) {
                if(cs==obj->DiscreteCSet[at_idx]) {
                  if(obj->DiscreteAtmToIdx[at_idx]>=0) {
                    base[0].sele[a]=true;
                    c++;
                  }
                }
              } else {
              } if(cs->AtmToIdx[at_idx]>=0) {
                base[0].sele[a]=true;
                c++;
              }
            }
          }
      }
		break;
	 case SELE_ALTs:
      WordPrimeCommaMatch(base[1].text);
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
	 case SELE_FLGs:
      sscanf(base[1].text,"%d",&flag);
      flag = (1<<flag);
		for(a=0;a<I->NAtom;a++)
		  {
			 if(I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags&flag)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
        }
		break;
	 case SELE_NTYs:
      if(sscanf(base[1].text,"%i",&numeric)!=1)
        ok=ErrMessage("Selector","Invalid Numeric Type.");
      if(ok)
        for(a=0;a<I->NAtom;a++)
          {
            if(I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].customType ==
               numeric) {
              base[0].sele[a]=true;
              c++;
            }
          }
      break;
	 case SELE_RSIs:
		if((p=strstr(base[1].text,":"))/* range */
         ||(p=strstr(base[1].text,"-")))/* range */
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
		else /* not a range */ {
        WordPrimeCommaMatch(base[1].text);
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
      }
		break;
	 case SELE_RSNs:
      WordPrimeCommaMatch(base[1].text);
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
          sele=I->ID[sele];
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
      index = -1;
      if((np=strstr(base[1].text,"`"))) {
        *np=0;
        if(sscanf(np+1,"%d",&index)!=1)
          index = -1;
        else
          index--;
      }
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
          if(index>=0) {
            for(a=0;a<I->NAtom;a++)
              {
                if(I->Table[a].model==model)
                  if(I->Table[a].atom==index)
                    {
                      base[0].sele[a]=true;
                      c++;
                    }
                  else {
                    base[0].sele[a]=false;                    
                  }
                else
                  base[0].sele[a]=false;
              }
          } else {
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
		  }
		else {
        PRINTFB(FB_Selector,FB_Errors)
          " Selector-Error: invalid model '%s'.\n",base[1].text
          ENDFB;
        ok=false;
      }
		break;
	 }
  PRINTFD(FB_Selector)
	 " SelectorSelect1:  %d atoms selected.\n",c
    ENDFD;
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
  int exact;
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
	 case SELE_PCHx:
	 case SELE_FCHx:
	 case SELE_BVLx:
	 case SELE_QVLx:
      oper=WordKey(AtOper,base[1].text,4,I->IgnoreCase,&exact);
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
            case SELE_QVLx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->q>comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_PCHx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->partialCharge>comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_FCHx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->formalCharge>comp1) {
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
            case SELE_QVLx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->q<comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;               
            case SELE_PCHx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->partialCharge<comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_FCHx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->formalCharge<comp1) {
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
                if(fabs(at1->b-comp1)<R_SMALL4) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_QVLx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->q-comp1)<R_SMALL4) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_PCHx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->partialCharge-comp1)<R_SMALL4) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_FCHx:
              for(a=0;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->formalCharge-comp1)<R_SMALL4) {
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
  
  PRINTFD(FB_Selector)
	 " SelectorSelect2: %d atoms selected.\n",c
    ENDFD;
  return(ok);
}
/*========================================================================*/
int SelectorLogic1(EvalElem *base)
{
  int a,b;
  int c=0;
  int flag;
  int n;
  int a0,a1,a2;
  AtomInfoType *at1,*at2;
  SelectorType *I=&Selector;
  ObjectMolecule *lastObj = NULL;
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
    case SELE_NGH1:
      base[1].sele=base[0].sele;
      base[0].sele=Alloc(int,I->NAtom);
      for(a=0;a<I->NAtom;a++) {
        base->sele[a]=0;
      }
      for(a=0;a<I->NAtom;a++) {
        if(base[1].sele[a]) {
          if(I->Obj[I->Table[a].model]!=lastObj) {
            lastObj = I->Obj[I->Table[a].model];
            ObjectMoleculeUpdateNeighbors(lastObj);
          }
          a0= I->Table[a].atom;
          n=lastObj->Neighbor[a0];
          n++;
          while(1) {
            a1=lastObj->Neighbor[n];
            if(a1<0) break;
            a2 = a1+lastObj->SeleBase;
            if(!base[1].sele[a2])
              base[0].sele[a2] =1;
            n+=2;
          }
        }
      }
      FreeP(base[1].sele);
      break;
    case SELE_BYO1: 
      base[1].sele=base[0].sele;
      base[0].sele=Alloc(int,I->NAtom);
      for(a=0;a<I->NAtom;a++) {
        base->sele[a]=0;
      }
      for(a=0;a<I->NAtom;a++) {
        if(base[1].sele[a]) {
          if(I->Obj[I->Table[a].model]!=lastObj) {
            lastObj = I->Obj[I->Table[a].model];
            b = a;
            while(b>=0) {
              if(I->Obj[I->Table[b].model]!=lastObj) 
                break;
              base[0].sele[b]=1;
              b--;
            }
            b=a+1;
            while(b<I->NAtom) {
              if(I->Obj[I->Table[b].model]!=lastObj) 
                break;
              base[0].sele[b]=1;
              b++;
            }
          }
        }
      }
      FreeP(base[1].sele);
      break;      
	 case SELE_BYR1: /* ASSUMES atoms are sorted by residue */
		for(a=0;a<I->NAtom;a++)
		  {
			 if(base[0].sele[a]) 
			   {
              at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
              b = a-1;
              while(b>=0) {
                if(!base[0].sele[b]) {
                  flag = false;
                  if(I->Table[a].model==I->Table[b].model)
                    {
                      at2=&I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                      if(at1->hetatm==at2->hetatm)
                        if(at1->chain[0]==at2->chain[0])
                          if(WordMatch(at1->resi,at2->resi,I->IgnoreCase)<0)
                            if(WordMatch(at1->segi,at2->segi,I->IgnoreCase)<0) {
                              base[0].sele[b]=true;
                              c++;
                              flag=1;
                            }
                    }
                  if(!flag)
                    break;
                }
                b--;
              }
              b = a + 1;
              while(b<I->NAtom) {
                if(!base[0].sele[b]) {
                  flag=false;
                  if(I->Table[a].model==I->Table[b].model)
                    {
                      at2=&I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                      if(at1->hetatm==at2->hetatm)
                        if(at1->chain[0]==at2->chain[0])
                          if(WordMatch(at1->resi,at2->resi,I->IgnoreCase)<0)
                            if(WordMatch(at1->segi,at2->segi,I->IgnoreCase)<0) {
                              base[0].sele[b]=true;
                              c++;
                              flag=1;
                            }
                    }
                if(!flag)
                  break;
                }
                b++;
              }
			   }
		  }
#ifdef _OLD_CODE
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
                           if(WordMatch(at1->segi,at2->segi,I->IgnoreCase)<0) {
                             base[0].sele[b]=true;
                             c++;
                           }
                   }
			   }
		  }
#endif
		break;
	 }
  PRINTFD(FB_Selector)
	 " SelectorLogic1: %d atoms selected.\n",c
    ENDFD;
  return(1);
}
/*========================================================================*/
int SelectorLogic2(EvalElem *base)
{
  int a,b;
  int c=0;
  SelectorType *I=&Selector;
  AtomInfoType *at1,*at2;
  TableRec *tr0,*tr2;
  int *s0,*s2;
  ObjectMolecule **obj;
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
              if(at1->resv==at2->resv)
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
	 case SELE_LIK2:
      s0  = base[0].sele;
      tr0  = I->Table;
      obj = I->Obj;
		for(a=0;a<I->NAtom;a++)
		  {
          if(*s0) {
            at1=&obj[tr0->model]->AtomInfo[tr0->atom];
            *s0=0;
            s2 = base[2].sele;
            tr2 = I->Table;
            for(b=0;b<I->NAtom;b++) {
              if(*s2) {
                at2=&obj[tr2->model]->AtomInfo[tr2->atom];
                if(at1->resv==at2->resv)
                  if(WordMatch(at1->name,at2->name,I->IgnoreCase)<0)
                    if(WordMatch(at1->resi,at2->resi,I->IgnoreCase)<0)
                      (*s0)=1;
              }
              s2++;
              tr2++;
            }
            if(*s0) c++;
          }
          s0++;
          tr0++;
		  }
		break;
	 }
  FreeP(base[2].sele);
  PRINTFD(FB_Selector)
	 " SelectorLogic2: %d atoms selected.\n",c
    ENDFD;
  return(1);
}
/*========================================================================*/
int SelectorOperator22(EvalElem *base)
{
  int c=0;
  int a,d,e;
  SelectorType *I=&Selector;
  ObjectMolecule *obj;

  float dist;
  float *v2;
  CoordSet *cs;
  int ok=true;
  int nCSet;
  MapType *map;
  int i,j,h,k,l;
  int n1,at,idx;

  switch(base[1].code)
	 {
	 case SELE_WIT_:
		if(!sscanf(base[2].text,"%f",&dist))
		  ok=ErrMessage("Selector","Invalid distance.");
		if(ok)
		  {
          if(dist<0.0) dist = 0.0;

          /* copy starting mask */
          for(a=0;a<I->NAtom;a++) {
            I->Flag2[a]=base[0].sele[a];
            base[0].sele[a]=false;
          }

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
                  if(obj->DiscreteFlag) {
                    if(cs==obj->DiscreteCSet[at])
                      idx=obj->DiscreteAtmToIdx[at];
                    else
                      idx=-1;
                  } else 
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
						nCSet=SelectorGetArrayNCSet(base[4].sele);
						for(e=0;e<nCSet;e++) {
						  for(a=0;a<I->NAtom;a++) {
							 if(base[4].sele[a])
								{
								  at=I->Table[a].atom;
								  obj=I->Obj[I->Table[a].model];
                          if(e<obj->NCSet) 
                            cs=obj->CSet[e];
                          else
                            cs=NULL;
								  if(cs) {
                            if(obj->DiscreteFlag) {
                              if(cs==obj->DiscreteCSet[at])
                                idx=obj->DiscreteAtmToIdx[at];
                              else
                                idx=-1;
                            } else 
                              idx=cs->AtmToIdx[at];
									 if(idx>=0) {
										v2 = cs->Coord+(3*idx);
										MapLocus(map,v2,&h,&k,&l);
										i=*(MapEStart(map,h,k,l));
										if(i) {
										  j=map->EList[i++];
										  while(j>=0) {
											 if(!base[0].sele[j])
                                    if(I->Flag2[j])
                                      if(within3f(I->Vertex+3*j,v2,dist)) 
                                        base[0].sele[j]=true;
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
          for(a=0;a<I->NAtom;a++)
            if(base[0].sele[a]) c++;
        }
      break;
    }
  FreeP(base[4].sele);
  PRINTFD(FB_Selector)
	 " SelectorOperator22: %d atoms selected.\n",c
    ENDFD;
  return(1);
}

static void remove_quotes(char *st)
{ 
  /* nasty */

  WordType store;
  char *p,*q;
  char *quote_start = NULL;
  char active_quote = 0;
  p = st;
  q = store;
  /*  printf("input [%s]\n",st);*/

  while(*p) {
    if(((*p)==34)||((*p)==39)) {
      if(quote_start&&(active_quote==*p)) { /* eliminate quotes... */
        while(quote_start<(q-1)) {
          *(quote_start)=*(quote_start+1);
          quote_start++;
        }
        q--;
        quote_start=NULL;
        p++;
      } else if(quote_start) {
        *(q++)=*(p++);
      } else {
        if(p==st) { /* at start => real quote */
          quote_start = q;
          active_quote = *p;
        } else if((*(p-1)=='+')||(*(p-1)==',')) { /* after separator => real quote */
          quote_start = q;
          active_quote = *p;
        }
        *(q++)=*(p++);
      }
    } else {
      if((*p=='+')&&(!quote_start))
        if(!((*(p+1)==0)||(*(p+1)==',')||(*(p+1)=='+')))
          *p=',';
      *(q++)=*(p++);
    }
  }
  *(q++) = 0;
  strcpy(st,store);

  /*  printf("output [%s]\n",st);*/
}

/*========================================================================*/
int *SelectorEvaluate(WordType *word)
{
  int level = 0;
  int depth = 0;
  int a,b,c = 0;
  int ok=true;
  unsigned int code;
  int valueFlag = 0; /* are we expecting? */
  int *result = NULL;
  int opFlag,opFlag2,maxLevel;
  char *q,*cc1,*cc2;
  int totDepth=0;
  int exact;
  char *np;
  OrthoLineType line;
  EvalElem *Stack=NULL,*e;
  SelectorType *I=&Selector;
  WordType tmpKW;

  Stack = Alloc(EvalElem,SelectorMaxDepth);

  Stack[0].sele=NULL;
  Stack[0].type=0;
  Stack[0].level=0;

  /* converts all keywords into code, adds them into a operation list */
  while(ok&&word[c][0])
	 {

      if(word[c][0]=='#') {
        if((!valueFlag)&&(!level)) {
          word[c][0]=0; /* terminate selection if we encounter a comment */
          word[c+1][0]=0;
          break;
        }
      }
		switch(word[c][0])
		  {
        case 0:
          break;
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
				  e->level=level<<4;
				  e->type=STYP_VALU;
              cc1 = word[c];
              cc2 = e->text;
              strcpy(e->text,cc1);
              
#if 0
              while(*cc1) { /* remove embedded quotes if any */
                if((*cc1!=34)&&(*cc1!=39)) {
                  *(cc2++)=*(cc1++);
                } else {
                  if(cc2!=e->text) 
                    if(*(cc2-1)=='+') { 
                      /* workaround for things like (alt A+''), 
                         which would fall back to atom name "A+" instead
                         of "A,", which is what this kludge fixes
                         -- this behavior necessary to accom. Na+, etc. */

                      *(cc2-1)=',';
                    }
                  cc1++;
                }
              }
              *cc2=0;
#endif

              remove_quotes(e->text);
				  valueFlag--;
				}
			 else if(valueFlag<0) /* operation parameter i.e. around X<-- */
				{
				  depth++;
				  e=Stack+depth;
				  e->level=level<<4;
				  e->type=STYP_PVAL;
				  strcpy(e->text,word[c]);
				  valueFlag++;
				} 
			 else
				{
              code=WordKey(Keyword,word[c],4,I->IgnoreCase,&exact);
              if(!code) {
                b=strlen(word[c])-1;
                if((b>2)&&(word[c][b]==';')) {
                  /* kludge to accomodate unnec. ';' usage */
                  word[c][b]=0;
                  code=WordKey(Keyword,word[c],4,I->IgnoreCase,&exact);
                }
                  
              }
              PRINTFD(FB_Selector)
                " Selector: code %x\n",code
                ENDFD;
              if((code>0)&&(!exact))  
                if(SelectorIndexByName(word[c])>=0)
                  code=0; /* favor selections over partial keyword matches */
				  if(code) 
					 { /* this is a known operation */
						depth++;
						e=Stack+depth;
						e->code=code;
						e->level=(level<<4)+((e->code&0xF0)>>4);
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
                    case STYP_OP22: 
							 valueFlag=-2;
							 break;
                    }
					 } else {
                  strcpy(tmpKW,word[c]); /* handle <object-name>[`index] syntax */
                  if((np=strstr(tmpKW,"`"))) { /* must be an object name */
                    *np=0;
                    if(SelectorIndexByName(tmpKW)>=0) {
                      depth++;
                      e=Stack+depth;
                      e->code=SELE_MODs;
                      e->level=(level<<4)+((e->code&0xF0)>>4);
                      e->type=STYP_SEL1;
                      valueFlag=1;
                      c--;
                    } else {
                      ok=ErrMessage("Selector","Unknown keyword or selection.");
                    }
                  } else { /* handle <selection-name> syntax */
                    if(SelectorIndexByName(tmpKW)>=0) {
                      depth++;
                      e=Stack+depth;
                      e->code=SELE_SELs;
                      e->level=(level<<4)+((e->code&0xF0)>>4);
                      e->type=STYP_SEL1;
                      valueFlag=1;
                      c--;
                    } else { 
                      ok=ErrMessage("Selector","Unknown keyword or selection.");
                    }
                  }
                }
            }
          break;
		  }
		if(ok) c++; /* go onto next word */
    }
  if(level>0)
    ok=ErrMessage("Selector","Malformed selection.");		
  if(ok) /* this is the main operation loop */
    {
      totDepth=depth;
      opFlag=true;
      maxLevel=-1;
      for(a=1;a<=totDepth;a++) {
        PRINTFD(FB_Selector)
          " Selector initial stack %d-%p lv: %x co: %d type: %x sele %p\n",
          a,Stack+a,Stack[a].level,Stack[a].code,Stack[a].type,Stack[a].sele
          ENDFD;
                 
        if(Stack[a].level>maxLevel) 
          maxLevel=Stack[a].level;
      }
      level=maxLevel;
      PRINTFD(FB_Selector)
        " Selector: maxLevel %d %d\n",maxLevel,totDepth
        ENDFD;
      if(level>=0) 
        while(ok) { /* loop until all ops at all levels have been tried */
          PRINTFD(FB_Selector)
            " Selector: new cycle...\n"
            ENDFD;
          depth = 1;
          opFlag=true;
          while(ok&&opFlag) { /* loop through all entries looking for ops at the current level */
            PRINTFD(FB_Selector)
              " Selector: lvl: %d de:%d-%p slv:%d co: %x typ %x sele %p td: %d\n",
              level,depth,Stack+depth,Stack[depth].level,Stack[depth].code,
              Stack[depth].type,Stack[depth].sele,totDepth
              ENDFD;
          
            opFlag=false;
            
            if(Stack[depth].level>=level) {
              Stack[depth].level=level; /* trim peaks */
            }
            if(ok) 
              if(depth>0)
                if((!opFlag)&&(Stack[depth].type==STYP_SEL0))
                  {
                    opFlag=true;
                    ok=SelectorSelect0(&Stack[depth]);
                  }
            if(ok)
              if(depth>1)
                if(Stack[depth-1].level>=Stack[depth].level)
                  {
                    if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_SEL1)
                       &&(Stack[depth].type==STYP_VALU))
                      { /* 1 argument selection operator */
                        opFlag=true;
                        ok=SelectorSelect1(&Stack[depth-1]);
                        for(a=depth+1;a<=totDepth;a++) 
                          Stack[a-1]=Stack[a];
                        totDepth--;
                      }
                    else if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_OPR1)
                       &&(Stack[depth].type==STYP_LIST))
                      { /* 1 argument logical operator */
                        opFlag=true;
                        ok=SelectorLogic1(&Stack[depth-1]);
                        for(a=depth+1;a<=totDepth;a++) 
                          Stack[a-1]=Stack[a];
                        totDepth--;
                      }
                  }
            if(ok)
              if(depth>2)
                if((Stack[depth-1].level>=Stack[depth].level)&&
                   (Stack[depth-1].level>=Stack[depth-2].level))
                  {
                    if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_OPR2)
                       &&(Stack[depth].type==STYP_LIST)
                       &&(Stack[depth-2].type==STYP_LIST))
                      { /* 2 argument logical operator */
                        ok=SelectorLogic2(&Stack[depth-2]);
                        opFlag=true;
                        for(a=depth+1;a<=totDepth;a++) 
                          Stack[a-2]=Stack[a];
                        totDepth-=2;
                      }
                    else if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_PRP1)
                       &&(Stack[depth].type==STYP_PVAL)
                       &&(Stack[depth-2].type==STYP_LIST))
                      { /* 2 argument logical operator */
                        ok=SelectorModulate1(&Stack[depth-2]);
                        opFlag=true;
                        for(a=depth+1;a<=totDepth;a++) 
                          Stack[a-2]=Stack[a];
                        totDepth-=2;
                      }
                  }
            if(ok)
              if(depth>2)
                if((Stack[depth-2].level>=Stack[depth-1].level)&&
                   (Stack[depth-2].level>=Stack[depth].level))
                  {            
                    if(ok&&(!opFlag)&&(Stack[depth-2].type==STYP_SEL2)
                       &&(Stack[depth-1].type==STYP_VALU)
                       &&(Stack[depth].type==STYP_VALU))
                      { /* 2 argument value operator */
                        ok=SelectorSelect2(&Stack[depth-2]);
                        opFlag=true;
                        for(a=depth+1;a<=totDepth;a++) 
                          Stack[a-2]=Stack[a];
                        totDepth-=2;
                      }
                  }
            if(ok)
              if(depth>3)
                if((Stack[depth-3].level>=Stack[depth].level)&&
                   (Stack[depth-3].level>=Stack[depth-1].level)&&
                   (Stack[depth-3].level>=Stack[depth-2].level))
                  {
                    if(ok&&(!opFlag)&&(Stack[depth-3].type==STYP_SEL3)
                       &&(Stack[depth].type==STYP_VALU)
                       &&(Stack[depth-1].type==STYP_VALU)
                       &&(Stack[depth-2].type==STYP_VALU))
                      { /* 2 argument logical operator */
                        /*								ok=SelectorSelect3(&Stack[depth-3]);*/
                        opFlag=true;
                        for(a=depth+1;a<=totDepth;a++) 
                          Stack[a-3]=Stack[a];
                        totDepth-=3;
                      }
                  }
            if(ok)
              if(depth>4)
                if((Stack[depth-3].level>=Stack[depth].level)&&
                   (Stack[depth-3].level>=Stack[depth-1].level)&&
                   (Stack[depth-3].level>=Stack[depth-2].level)&&
                   (Stack[depth-3].level>=Stack[depth-4].level))
                  {
                    if(ok&&(!opFlag)&&(Stack[depth-3].type==STYP_OP22)
                       &&(Stack[depth-1].type==STYP_VALU)
                       &&(Stack[depth-2].type==STYP_VALU)
                       &&(Stack[depth].type==STYP_LIST)
                       &&(Stack[depth-4].type==STYP_LIST))
                      { 
                        ok=SelectorOperator22(&Stack[depth-4]);
                        opFlag=true;
                        for(a=depth+1;a<=totDepth;a++) 
                          Stack[a-4]=Stack[a];
                        totDepth-=4;
                      }

                  }
            if(opFlag) {
              opFlag2=true; /* make note that we performed an operation */
              depth = 1; /* start back at the left hand side */
            } else {
              depth = depth + 1;
              opFlag=true;
              if(depth>totDepth)
                break;
            }
          }
          if(level)
            level--;
          else
            break;
        }
      depth=totDepth;
	 }
  if(ok)
	 {
      if(depth!=1){
        ok=ErrMessage("Selector","Malformed selection.");
      }
      else if(Stack[depth].type!=STYP_LIST)
        ok=ErrMessage("Selector","Invalid selection.");
      else
        result=Stack[totDepth].sele; /* return the selection list */
	 }
  if(!ok)
	 {
		for(a=1;a<=depth;a++) {
        PRINTFD(FB_Selector)
          " Selector: releasing %d %x %p\n",a,Stack[a].type,Stack[a].sele
          ENDFD;
		  if(Stack[a].type==STYP_LIST)
			 FreeP(Stack[a].sele);
      }
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
      PRINTFB(FB_Selector,FB_Errors) 
        "%s\n",line
        ENDFB;
		OrthoRestorePrompt();
	 }
  FreeP(Stack);
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
		if(w_flag) /* currently in a word, thus q is a valid pointer */
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
				  *q=0; /* terminate current word */
				  c++;
				  VLACheck(r,WordType,c); /* add new word */
				  q=r[c-1];
				  *q++=*p;
				  *q=0;  /* terminate current word */
				  w_flag=false;
				  break;
				default:
				  *q++=*p;
				  break;
				}
		  }
		else /*outside a word -- q is undefined */
		  {
			 switch(*p)
				{
				case '*': /* special case */
				  c++;
				  VLACheck(r,WordType,c);
				  q=r[c-1];
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
				  c++;
				  VLACheck(r,WordType,c);
				  q=r[c-1];
				  *q++=(*p);
				  *q=0;
				  break;
            case ' ':
              break;
				default:
				  w_flag=true;
				  c++;
              VLACheck(r,WordType,c);
				  q=r[c-1];
				  *q++=*p;
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
  if(Feedback(FB_Selector,FB_Debugging)) 
	 {
		c=0;
		while(r[c][0])
		  {
			 fprintf(stderr,"word: %s\n",r[c]);
			 c++;
		  }
	 }
  return(r);
}
/*========================================================================*/
void SelectorFree(void)
{
  SelectorType *I = &Selector;
  SelectorClean();
  VLAFreeP(I->Member);
  VLAFreeP(I->Name);
  VLAFreeP(I->ID);
}
/*========================================================================*/
void SelectorInit(void)
{
  SelectorType *I = &Selector;
  I->NSelection = 0;
  I->NActive=0;
  I->TmpCounter = 0;
  I->Name = VLAlloc(WordType,100);
  I->ID = VLAlloc(int,100);
  
  I->Member = (MemberType*)VLAMalloc(10000,sizeof(MemberType),5,true);
  I->NMember=0;
  I->FreeMember=-1;
  I->NCSet=0;
  I->Table=NULL;
  I->Obj=NULL;
  I->IgnoreCase=true;
  I->Flag1=NULL;
  I->Flag2=NULL;
  I->Vertex=NULL;
}
/*========================================================================*/


DistSet *SelectorGetDistSet(int sele1,int state1,int sele2,int state2,
                            int mode,float cutoff,float *result)
{
  SelectorType *I=&Selector;
  int *vla=NULL;
  int c;
  float dist;
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
  float dist_sum=0.0;
  int dist_cnt = 0;
  
  *result = 0.0;
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


          if(obj1->DiscreteFlag) {
            if(cs1==obj1->DiscreteCSet[at1]) {
              idx1=obj1->DiscreteAtmToIdx[at1];
            } else {
              idx1=-1;
            }
          } else {
            idx1=cs1->AtmToIdx[at1];
          }
          
          if(obj2->DiscreteFlag) {
            if(cs2==obj2->DiscreteCSet[at2]) {
              idx2=obj2->DiscreteAtmToIdx[at2];
            } else {
              idx2=-1;
            }
            
          } else {
            idx2=cs2->AtmToIdx[at2];
          }
            
          if((idx1>=0)&&(idx2>=0)) {
            dist=diff3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2);
            
            if(dist<cutoff) {
              dist_cnt++;
              dist_sum+=dist;
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
              
              nv+=2;
            }
          }
        }
      }
    }
  }
  if(dist_cnt)
    (*result)=dist_sum/dist_cnt;
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


