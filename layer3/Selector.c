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

#include"Base.h"
#include"Vector.h"
#include"Debug.h"
#include"Err.h"
#include"Word.h"
#include"Util.h"
#include"MemoryDebug.h"
#include"Selector.h"
#include"Executive.h"
#include"ObjectMolecule.h"

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
} TableRec;

typedef struct {
  WordType *Name;
  int NSelection;
  MemberType *Member;
  int NMember;
  ObjectMolecule **Obj;
  TableRec *Table;
  int NAtom;
  int NModel;
  int IgnoreCase;
} SelectorType;

SelectorType Selector;
int SelectorUpdateTable(void);


int SelectorModulate1(EvalElem *base);
int SelectorSelect0(EvalElem *base);
int SelectorSelect1(EvalElem *base);
int SelectorSelect3(EvalElem *base);
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
#define SELE_SEL3 5
#define SELE_LIST 6
#define SELE_PRP1 7
#define SELE_PVAL 0

static WordType Keyword[] = 
{
  "not",      "NOT1",
  "byresi",   "BYR1",
  "and",      "AND2",
  "or",       "OR_2",
  "in",       "IN_2",
  "all",      "ALLz", /* 0 parameter */
  "none",     "NONz", /* 0 parameter */
  "around",   "ARD_", /* 1 parameter */
  "expand",   "EXP_", /* 1 parameter */
  "name",     "NAMs",
  "resi",     "RSIs",
  "chain",    "CHNs",
  "segi",     "SEGs",
  "model",    "MODs",
  "index",    "IDXs",
  "resn",     "RSNs",
  "%",        "SELs",
  "attrib",   "VALx", /* 3 operand selection operator */ 
  ""
};
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
  if(c) 
    printf(" Selector: %s created with %d atoms.\n",name,c);
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
  modelCnt=0;
  while(ExecutiveIterateObject(&o,&hidden))
	 {
		if(o->type==cObjectMolecule)
		  {
			 obj=(ObjectMolecule*)o;
			 c+=obj->NAtom;
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
  int a,b,d,e;
  int c=0;
  float dist;
  float *v;
  CoordSet *cs1,*cs2;
  int ok=true;
  SelectorType *I=&Selector;
  base[1].sele=base[0].sele;
  base->sele=Alloc(int,I->NAtom);
  for(a=0;a<I->NAtom;a++)
	 base[0].sele[a]=false;
  ErrChkPtr(base->sele);
  switch(base[1].code)
	 {
	 case 'ARD_': /* Currently an c*N^2 operations - TODO fix */
	 case 'EXP_':
		if(!sscanf(base[2].text,"%f",&dist))
		  ok=ErrMessage("Selector","Invalid distance.");
		if(ok)
		  for(a=0;a<I->NAtom;a++)
			 {
				if(base[1].sele[a])
				  for(e=0;e<I->Obj[I->Table[a].model]->NCSet;e++)
                if(I->Obj[I->Table[a].model]->CSet[e])
                  {
						  cs1 = I->Obj[I->Table[a].model]->CSet[e];
                    v=cs1->Coord+(3*cs1->AtmToIdx[I->Table[a].atom]);
                    for(b=0;b<I->NAtom;b++)
                      if((!base[0].sele[b])&&((base[1].code=='EXP_')||(!base[1].sele[b]))) /*exclude current selection */
                        {
                          for(d=0;d<I->Obj[I->Table[b].model]->NCSet;d++)
                            {
										cs2 = I->Obj[I->Table[b].model]->CSet[e];
                              if(diff3f(cs2->Coord+(3*cs2->AtmToIdx[I->Table[b].atom]),v)<dist)
                                {
                                  base[0].sele[b]=true;
                                  c++;
                                  break;
                                }
                            }
                        }
                  }
			 }
		break;
	 }
  FreeP(base[1].sele);
  if(DebugSelector&DebugState)
	 printf("SelectorModulate0: %c%c%c%c : %d atoms selected.\n",
			  base[1].code>>24,base[1].code>>16,
			  base[1].code>>8,base[1].code&0xFF,c);
  return(ok);
  
}
/*========================================================================*/
int SelectorSelect0(EvalElem *base)
{
  int a;
  int c=0;
  SelectorType *I=&Selector;
  base->type=SELE_LIST;
  base->sele=Alloc(int,I->NAtom);
  ErrChkPtr(base->sele);
  switch(base->code)
	 {
	 case 'NONz':
		for(a=0;a<I->NAtom;a++)
		  base[0].sele[a]=false;
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
			 if(WordMatch(base[1].text,I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].name,I->IgnoreCase)<0)
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
			 if(WordMatch(base[1].text,I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].segi,I->IgnoreCase)<0)
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
			 if(base[1].text==
				 I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].chain)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else if(I->IgnoreCase) 
				if(tolower(base[1].text[0])==
					(tolower(I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].chain[0])))
				  {
					 base[0].sele[a]=true;
					 c++;
				  }
				else
				  base[0].sele[a]=false;
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
				if(WordMatch(base[1].text,I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].resi,I->IgnoreCase)<0)
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
			 if(WordMatch(base[1].text,I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].resn,I->IgnoreCase)<0)
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
int SelectorSelect3(EvalElem *base)
{
  base->type=SELE_LIST;
  return(1);
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
  int opFlag;
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
							 e->type=SELE_SEL3;
							 valueFlag=3;
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
						e->code=code;
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
								ok=SelectorSelect3(&Stack[depth-3]);
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
				case '(': /* single word terminators */ 
				case ')':
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
				case '(': /* single word terminators */
				case ')':
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
}
/*========================================================================*/
void SelectorInit(void)
{
  SelectorType *I = &Selector;
  I->NSelection = 0;
  I->Name = VLAlloc(WordType,10);
  I->Member = (MemberType*)VLAMalloc(1000,sizeof(MemberType),5,true);
  I->NMember=0;
  I->Table=NULL;
  I->Obj=NULL;
  I->IgnoreCase=true;
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


