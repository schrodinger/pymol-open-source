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

#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>

#include"Base.h"
#include"Debug.h"
#include"OOMac.h"
#include"Vector.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Map.h"
#include"Selector.h"
#include"Color.h"
#include"ObjectMolecule.h"
#include"Ortho.h"
#include"Util.h"
#include"Vector.h"
#include"Selector.h"
#include"Matrix.h"
#include"Scene.h"
#include"PUtils.h"

void ObjectMoleculeRender(ObjectMolecule *I,int frame,CRay *ray,Pickable **pick);
void ObjectMoleculeCylinders(ObjectMolecule *I);
CoordSet *ObjectMoleculePDBStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);
CoordSet *ObjectMoleculeMOLStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr);
void ObjectMoleculeAppendAtoms(ObjectMolecule *I,AtomInfoType *atInfo,CoordSet *cset);

void ObjectMoleculeFree(ObjectMolecule *I);
void ObjectMoleculeUpdate(ObjectMolecule *I);
int ObjectMoleculeGetNFrames(ObjectMolecule *I);

void ObjectMoleculeDescribeElement(ObjectMolecule *I,int index);

void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op);
void ObjectMoleculeExtendIndices(ObjectMolecule *I);
void ObjectMoleculeSort(ObjectMolecule *I);
void ObjectMoleculeMerge(ObjectMolecule *I,AtomInfoType *ai,CoordSet *cs);

int ObjectMoleculeConnect(int **bond,AtomInfoType *ai,CoordSet *cs,float cutoff);
void ObjectTransformTTTf(ObjectMolecule *I,float *ttt);

static int BondInOrder(int *a,int b1,int b2);
static int BondCompare(int *a,int *b);

#define MAXLINELEN 255

/*========================================================================*/
static int BondInOrder(int *a,int b1,int b2)
{
  return(BondCompare(a+b1,a+b2)<=0);
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
void ObjectMoleculeExtendIndices(ObjectMolecule *I)
{
  int a;
  for(a=0;a<I->NCSet;a++)
	 if(I->CSet[a])
		I->CSet[a]->fExtendIndices(I->CSet[a],I->NAtom);
}
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
void ObjectMoleculeSort(ObjectMolecule *I)
{
  int *index,*outdex;
  int a,b;
  CoordSet *cs;
  AtomInfoType *atInfo;
  
  index=AtomInfoGetSortedIndex(I->AtomInfo,I->NAtom,&outdex);
  
  for(a=0;a<I->NBond;a++) { /* bonds */
	I->Bond[a*2]=outdex[I->Bond[2*a]];
	I->Bond[a*2+1]=outdex[I->Bond[a*2+1]];
  }
  
  for(a=0;a<I->NCSet;a++) { /* coordinate set mapping */
	cs=I->CSet[a];
	if(cs) {
	  for(b=0;b<cs->NIndex;b++)
		cs->IdxToAtm[b]=outdex[cs->IdxToAtm[b]];
	  for(b=0;b<I->NAtom;b++)
		cs->AtmToIdx[b]=-1;
	  for(b=0;b<cs->NIndex;b++)
		cs->AtmToIdx[cs->IdxToAtm[b]]=b;
	}
  }
  
  atInfo=(AtomInfoType*)VLAMalloc(I->NAtom,sizeof(AtomInfoType),5,true);
  /* autozero here is important */
  for(a=0;a<I->NAtom;a++)
	atInfo[a]=I->AtomInfo[index[a]];
  VLAFreeP(I->AtomInfo);
  I->AtomInfo=atInfo;

  AtomInfoFreeSortedIndexes(index,outdex);
  /*
  for(a=0;a<I->NAtom;a++)
	 {
		printf("%d %s\n",I->AtomInfo[a].resv,I->AtomInfo[a].chain);
		}*/
}
/*========================================================================*/
CoordSet *ObjectMoleculeMOLStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom,nBond,nType;
  int a,c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL;
  char cc[MAXLINELEN],resn[MAXLINELEN];
  float *f;
  int *ii,*bond=NULL;
  int NColor,CColor,HColor,OColor,SColor,MColor;
  int color=0;
  int ok=true;

  NColor=ColorGetIndex("nitrogen");
  CColor=ColorGetIndex("carbon");
  HColor=ColorGetIndex("hydrogen");
  OColor=ColorGetIndex("oxygen");
  SColor=ColorGetIndex("sulfer");
  MColor=ColorGetIndex("magenta");
  
  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  p=ncopy(resn,p,sizeof(ResName)-1);
  UtilCleanStr(resn);
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
			 p++; 
			 p=ncopy(atInfo[a].name,p,3);
			 UtilCleanStr(atInfo[a].name);
			 
			 atInfo[a].visRep[0] = true; /* show sticks by default */
			 for(c=1;c<cRepCnt;c++) {
				atInfo[a].visRep[c] = false;
			 }
		  }
		  if(ok&&atInfo) {
			 strcpy(atInfo[a].resn,resn);
			 atInfo[a].hetatm=true;
			 AtomInfoAssignParameters(atInfo);
			 switch ( atInfo[a].name[0] ) /* need to move this stuff into parameters */
				{
				case 'N' : color = NColor; break;
				case 'C' : 
				  switch(atInfo[a].name[1]) {
				  case 0:
              case 32:
                color = CColor; break;
				  case 'l':
				  case 'L':
				  default:
                color = MColor; break;
				  }
				  break;
				case 'O' : color = OColor; break;
				case 'I' : color = MColor; break;
				case 'P' : color = MColor; break;
				case 'B' : color = MColor; break;
				case 'S' : color = SColor; break;
				case 'F' : color = MColor; break;
				case 'H' : color=HColor; break;
				default  : color=MColor; break;
				}
			 atInfo[a].color=color;
		  }
		  p=nextline(p);
		  if(!ok)
			 break;
		}
  }
  if(ok) {
	 bond=VLAlloc(int,2*nBond);
	 ii=bond;
	 for(a=0;a<nBond;a++)
		{
		  if(ok) {
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",ii++)!=1)
				ok=ErrMessage("ReadMOLFile","bad bond");
		  }
		  
		  if(ok) {  
			 p=ncopy(cc,p,3);
			 if(sscanf(cc,"%d",ii++)!=1)
				ok=ErrMessage("ReadMOLFile","bad bond");
		  }
		  if(!ok)
			 break;
		  p=nextline(p);
		}
	 ii=bond;
	 for(a=0;a<nBond;a++) {
		(*(ii++))--; /* adjust bond indexs down one */
		(*(ii++))--; 
	 }
  }
  if(ok&&atInfo) 
	 while(ok&&(*p)) {
		if(*p== '>') {
		  p=ncopy(cc,p,MAXLINELEN);
		  if(strstr(cc,"<LUDI.TYPES>")) {
			 p=nextline(p);
			 p=wcopy(cc,p,MAXLINELEN);
			 p=wcopy(cc,p,MAXLINELEN);
			 if(sscanf(cc,"%d",&nType)!=1)
				ok=ErrMessage("ObjectMoleculeLoadMOLFile","bad ludi record");
			 if(ok)
				if(nType!=nAtom)
				  ok=ErrMessage("ObjectMoleculeLoadMOLFile","ludi atom count mismatch");
			 if(ok) {
				for(a=0;a<nType;a++) {
				  p=nextline(p);
				  p+=16;
				  p=ncopy(cc,p,6);
				  if(sscanf(cc,"%d",&c)!=1)
					 ok=ErrMessage("ObjectMoleculeLoadMOLFile","missing ludi type count");
				  if(ok) {
					 if(!c) 
						{
						  switch(atInfo[a].name[0]) {
						  case 'S':
							 atInfo[a].ludiType=17;
							 break;
						  case 'N':
							 atInfo[a].ludiType=18;
							 break;
						  default:
							 atInfo[a].ludiType=0;
						  }
						}
					 else {
						p=ncopy(cc,p,6);
						if(sscanf(cc,"%d",&atInfo[a].ludiType)!=1) {
						  ok=ErrMessage("ObjectMoleculeLoadMOLFile","missing ludi atom type");
						} else {
						  if(atInfo[a].name[0]=='0') 
							 if(atInfo[a].ludiType==6) /* convert alcohol O into weak hbond acceptor */
								atInfo[a].ludiType=-3;
						}
					 }
				  }
				}
			 }
		  } else if(strstr(cc,"<CUSTOM.TYPES>")) {
			 p=nextline(p);
			 p=nextline(p);
			 p=wcopy(cc,p,MAXLINELEN);
			 if(sscanf(cc,"%d",&nType)!=1)
				ok=ErrMessage("ObjectMoleculeLoadMOLFile","bad custom type record");
			 if(ok)
				if(nType!=nAtom)
				  ok=ErrMessage("ObjectMoleculeLoadMOLFile","custom atom type count mismatch");
			 if(ok) {
				for(a=0;a<nType;a++) {
				  p=nextline(p);
				  p+=8;
				  p=ncopy(cc,p,3);
				  if(sscanf(cc,"%d",&atInfo[a].customType)!=1) {
					 ok=ErrMessage("ObjectMoleculeLoadMOLFile","bad custom atom type record-1");
				  } else {
					 p+=1;
					 p=ncopy(cc,p,3);
					 if(sscanf(cc,"%d",&atInfo[a].customFlag)!=1) {
						ok=ErrMessage("ObjectMoleculeLoadMOLFile","bad custom atom type record-2");
					 } else {
					   p+=1;
					   p=ncopy(cc,p,6);
					   if(sscanf(cc,"%f",&atInfo[a].vdw)!=1) {
						 ok=ErrMessage("ObjectMoleculeLoadMOLFile","bad van der waals radius");
					   } 
					 }
				  } 
				}
			 }
		  }
		}
		p=nextline(p);
	 }

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

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadMOLStr(ObjectMolecule *I,char *MOLStr,int frame)
{
  int ok = true;
  CoordSet *cset=NULL;
  AtomInfoType *atInfo;
  
  atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
  cset=ObjectMoleculeMOLStr2CoordSet(MOLStr,&atInfo);
  
  if(!cset) 
	 {
		VLAFreeP(atInfo);
		ok=false;
	 }
  
  if(ok)
	 {
		if(!I) 
		  I=(ObjectMolecule*)ObjectMoleculeNew();
		if(frame<0)
		  frame=I->NCSet;
		if(I->NCSet<=frame)
		  I->NCSet=frame+1;
		VLACheck(I->CSet,CoordSet*,frame);
		
		cset->fAppendIndices(cset,I->NAtom);
		cset->Obj=I;
		ObjectMoleculeAppendAtoms(I,atInfo,cset);
		I->CSet[frame] = cset;
		I->CSet[frame]->fInvalidateRep(I->CSet[frame],-1,0);
		SceneCountFrames();
		ObjectMoleculeExtendIndices(I);
	 }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadMOLFile(ObjectMolecule *obj,char *fname,int frame)
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
		I=ObjectMoleculeReadMOLStr(obj,buffer,frame);
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
  AtomInfoType *atInfo = NULL;
  int AFlag;
  int atomCount;
  char shstr[12];
  float dummy[3];
  char t;
  int NColor,CColor,HColor,OColor,SColor,MColor;
  int color=0;

  NColor=ColorGetIndex("nitrogen");
  CColor=ColorGetIndex("carbon");
  HColor=ColorGetIndex("hydrogen");
  OColor=ColorGetIndex("oxygen");
  SColor=ColorGetIndex("sulfer");
  MColor=ColorGetIndex("magenta");

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  while(*p)
	 {
		if(*p == 'A') if(*(p+1)=='T') if(*(p+2)=='O') if(*(p+3)=='M')
		  nAtom++;
		if( *p == 'H') if(*(p+1)=='E') if(*(p+2)=='T') if(*(p+3)=='A')
		  if(*(p+4)=='T') if(*(p+5)=='M')
			 nAtom++;
		p++;
		while(*p)
		  if(*p==0xA)
			 {
				p++;
				break;
			 }
		  else
			 p++;
	 }
  for(a=0;a<255;a++) /*to prevent hopping over end of file*/
	 *p++=0;
  
  coord=VLAlloc(float,3*nAtom);
  if(atInfo)
	 VLACheck(atInfo,AtomInfoType,nAtom);

  p=buffer;
  if(DebugState & DebugMolecule) {
	 printf(" ObjectMoleculeReadPDB: Found %i atoms...\n",nAtom);
	 fflush(stdout);
  }
  fflush(stdout);
  a=0;
  atomCount=0;
  
  while(*p)
	 {
		AFlag=false;
		if(*p == 'A') if(*(p+1)=='T') if(*(p+2)=='O') if(*(p+3)=='M')
		  AFlag = 1;
		if( *p == 'H') if(*(p+1)=='E') if(*(p+2)=='T') if(*(p+3)=='A')
		  if(*(p+4)=='T') if(*(p+5)=='M')
			 AFlag = 2;
		if(AFlag)
		  {
			 llen=0;
			 while((*(p+llen))&&((*(p+llen))!=13)&&((*(p+llen))!=10)) {
				llen++;
			 }
			 t=*(p+26);
			 *(p+26) = 0;
			 sscanf(p+22,"%s",shstr);
			 *(p+26) = t;
			 if(atInfo)
				{
				  strcpy(atInfo[atomCount].resi,shstr);
				  sscanf(shstr,"%d",&atInfo[atomCount].resv);
				}

			 t=*(p+38);
			 *(p+38) = 0;
			 sscanf(p+30,"%f",coord+a);
			 *(p+38) = t;

			 t=*(p+46);
			 *(p+46) = 0;
			 sscanf(p+38,"%f",coord+a+1);
			 *(p+46) = t;

			 t=*(p+54);
			 *(p+54) = 0;
			 sscanf(p+46,"%f",coord+a+2);
			 *(p+54) = t;
			 
			 t = *(p+60);
			 *(p+60) = 0;
			 sscanf(p+54,"%f",&dummy[0]);
			 if(atInfo)
				atInfo[atomCount].q=dummy[0];
			 *(p+60) = t;

			 t = *(p+66);
			 *(p+66) = 0;
			 sscanf(p+60,"%f",&dummy[0]);
			 if(atInfo)
				atInfo[atomCount].b=dummy[0];
			 *(p+66) = t;
			 
			 t = *(p+16);
			 *(p+16) = 0;
			 sscanf(p+12,"%s",shstr);
			 if(atInfo) 
				strcpy(atInfo[atomCount].name,shstr);
			 *(p+16) = t;
			 
			 t = *(p+20);
			 *(p+20) = 0;
			 sscanf(p+17,"%s",shstr);
			 if(atInfo) 
				strcpy(atInfo[atomCount].resn,shstr);
			 *(p+20) = t;
			 
			 if(llen>=73) {
				t = *(p+76);
				*(p+76) = 0;
				shstr[0]=0;
				sscanf(p+72,"%s",shstr);
				if(atInfo) 
				  strcpy(atInfo[atomCount].segi,shstr);
				*(p+76) = t;
				}

			 if(atInfo)
				{
				  atInfo[atomCount].chain[0]=
					 atInfo[atomCount].chain[0] = *(p+21);
				  if((atInfo[atomCount].chain[0])==32)
					 atInfo[atomCount].chain[0] = 0;
				  else
					 atInfo[atomCount].chain[1] = 0;

				  atInfo[atomCount].visRep[0] = true;

				  for(c=1;c<cRepCnt;c++) {
					 atInfo[atomCount].visRep[c] = false;
				  }
				  if(AFlag==1) 
					 atInfo[atomCount].hetatm=0;
				  else
					 atInfo[atomCount].hetatm=1;

				  AtomInfoAssignParameters(&atInfo[atomCount]);

				  switch ( atInfo[atomCount].name[0] ) /* need to move this stuff into parameters */
					{
					case 'N' : color = NColor; break;
					case 'C' : color = CColor; break;
					case 'O' : color = OColor; break;
					case 'I' : color = MColor; break;
					case 'P' : color = MColor; break;
					case 'B' : color = MColor; break;
					case 'S' : color = SColor; break;
					case 'F' : color = MColor; break;
					case 'H' : color=HColor; break;
					default  : color=MColor; break;
					}
				  atInfo[atomCount].color=color;
				}
			 a+=3;
			 atomCount++;
		  }
		p++;
		while(*p)
		  if(*p==0xA)
			 {
				p++;
				break;
			 }
		  else
			 p++;
	 }
  cset = CoordSetNew();
  cset->NIndex=nAtom;
  cset->Coord=coord;
  if(atInfoPtr)
	 *atInfoPtr = atInfo;
  return(cset);
}
/*========================================================================*/
void ObjectMoleculeMerge(ObjectMolecule *I,AtomInfoType *ai,CoordSet *cs)
{
  int *index,*outdex,*a2i,*i2a,*bond=NULL;
  int a,b,c,lb,ac,a1,a2;
  int found;
  int nAt,nBd,nBond;
  int expansionFlag = false;
  AtomInfoType *ai2;
  
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
		if(a2 >= I->NAtom) { 
		  I->AtomInfo[a2]=ai[a1]; /* copy atom info */
		}
    }
  
  cs->NAtIndex = nAt;
  I->NAtom = nAt;
  
  FreeP(cs->AtmToIdx);
  FreeP(cs->IdxToAtm);
  cs->AtmToIdx = a2i;
  cs->IdxToAtm = i2a;
  
  for(a=0;a<cs->NAtIndex;a++)
    cs->AtmToIdx[a]=-1;
  for(a=0;a<cs->NIndex;a++)
    cs->AtmToIdx[cs->IdxToAtm[a]]=a;
  
  VLAFreeP(ai);
  AtomInfoFreeSortedIndexes(index,outdex);

  /* now find and integrate and any new bonds */
  if(expansionFlag) { /* expansion flag means we have introduced at least 1 new atom */
    nBond = ObjectMoleculeConnect(&bond,I->AtomInfo,cs,0.2);
    
    if(nBond) {
      index=Alloc(int,nBond);
      
      c=0;
      b=0;  
      for(a=0;a<nBond;a++) {
        found=false;
        lb=b;
        while(b<I->NBond) {
          ac=BondCompare(bond+a,I->Bond+b);
          if(!ac) {
            found=true;
            break;
          } else if(ac<0) {
            break;
          }
          b++;
        }
        if(found) {
          index[a]=b;
          b++;
		} else {
		  index[a]=I->NBond+c;
		  c++;
		  b=lb;
		}
	  }
      /* first, reassign atom info for matched atoms */
      if(c) {
        /* allocate additional space */
        nBd=I->NBond+c;
        
        VLACheck(I->Bond,int,nBd*2);
        
        for(a=0;a<nBond;a++)
          {
            a2=index[a];
            if(a2 >= I->NBond) { 
              I->Bond[2*a2]=bond[2*a]; /* copy bond info */
              I->Bond[2*a2+1]=bond[2*a+1]; /* copy bond info */
            }
          }
        I->NBond=nBd;
      }
      VLAFreeP(bond);
      FreeP(index);
    }
  }
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadPDBStr(ObjectMolecule *I,char *PDBStr,int frame)
{
  float *cord,*coord;
  float *extent;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo;
  int ok=true;
  int isNew = true;
  int a;
  unsigned int nAtom = 0;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMolecule*)ObjectMoleculeNew();
		atInfo = I->AtomInfo;
		isNew = true;
	 } else {
		atInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
		isNew = false;
	 }
	 cset=ObjectMoleculePDBStr2CoordSet(PDBStr,&atInfo);	 
	 coord=cset->Coord;
	 nAtom=cset->NIndex;

	 if(isNew) {	 
		I->Obj.extent[0]=coord[0];
		I->Obj.extent[1]=coord[0];
		I->Obj.extent[2]=coord[1];
		I->Obj.extent[3]=coord[1];
		I->Obj.extent[4]=coord[2];
		I->Obj.extent[5]=coord[2];
	 }
	 
	 cord = coord-1;
	 for(a=0;a<nAtom;a++)
		{
		  extent = I->Obj.extent; /* a monument to obfuscation and pseudo-optimization! */
		  if(*(++cord)<*(  extent))   *extent=*cord;
		  if(*(  cord)>*(++extent))	*extent=*cord;
		  if(*(++cord)<*(++extent))	*extent=*cord;
		  if(*(  cord)>*(++extent))	*extent=*cord;
		  if(*(++cord)<*(++extent))	*extent=*cord;
		  if(*(  cord)>*(++extent))	*extent=*cord;
		}

  }

  /* prepare the cset for inclusion into the object */
  if(ok) {
	cset->fEnumIndices(cset);
	cset->fInvalidateRep(cset,-1,0);
	if(isNew) {		
	  for(a=0;a<3;a++)
		I->Obj.center[a] = (I->Obj.extent[2*a]+I->Obj.extent[2*a+1])/2.0;
	  I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
	} else {
	  ObjectMoleculeMerge(I,atInfo,cset); /* NOTE: will release atInfo */
	}
  }
  if(ok) { /* now, include it */
	 cset->Obj = I;
	 if(isNew) I->NAtom=nAtom;
	 if(frame<0) frame=I->NCSet;
	 VLACheck(I->CSet,CoordSet*,frame);
	 if(I->NCSet<=frame) I->NCSet=frame+1;
	 I->CSet[frame] = cset;
	 if(isNew) I->NBond = ObjectMoleculeConnect(&I->Bond,I->AtomInfo,cset,0.2);
	 SceneCountFrames();
	 ObjectMoleculeExtendIndices(I);
	 ObjectMoleculeSort(I);
  }
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadPDBFile(ObjectMolecule *obj,char *fname,int frame)
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

		I=ObjectMoleculeReadPDBStr(obj,buffer,frame);

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
	 I->Bond=VLAlloc(int,nBond*2);
  VLACheck(I->Bond,int,nBond*2);
  ii=I->Bond+I->NBond*2;
  si=cs->TmpBond;
  for(a=0;a<cs->NTmpBond;a++)
	 {
		*(ii++)=cs->IdxToAtm[*(si++)];
		*(ii++)=cs->IdxToAtm[*(si++)];
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
void ObjectTransformTTTf(ObjectMolecule *I,float *ttt) 
{
  int b;
  CoordSet *cs;
  for(b=0;b<I->NCSet;b++)
	{
	  cs=I->CSet[b];
	  if(cs) {
		cs->fInvalidateRep(I->CSet[b],cRepAll,0);
		MatrixApplyTTTfn3f(cs->NIndex,cs->Coord,ttt,cs->Coord);
	  }
	}
}
/*========================================================================*/
void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op)
{
  int a,b,s;
  int a1;
  float r;
  float v1[3],v2,*vv1,*vv2;
  int inv_flag;
  int hit_flag = false;
  int ok = true;

  if(sele>=0) {
	SelectorUpdateTable();
	switch(op->code) {
	case 'PDB1':
	  for(b=0;b<I->NCSet;b++)
		if(I->CSet[b])
		  {
			if((b==op->i1)||(op->i1<0))
			  for(a=0;a<I->NAtom;a++)
				{
				  s=I->AtomInfo[a].selEntry;
				  while(s) 
					{
					  if(SelectorMatch(s,sele))
						{
						  
						  CoordSetAtomToPDBStrVLA(&op->charVLA,&op->i2,I->AtomInfo+a,
												  I->CSet[b]->Coord+(3*I->CSet[b]->AtmToIdx[a]),op->i3);
						  op->i3++;
						}
					  s=SelectorNext(s);
					}
				}
		  }
	  break;
	  
	default:
	   for(a=0;a<I->NAtom;a++)
		 {
		   switch(op->code) {
		   case 'COLR':
		   case 'VISI':
		   case 'TTTF':
         case 'ALTR':
			 s=I->AtomInfo[a].selEntry;
			 while(s) 
			   {
				 if(SelectorMatch(s,sele))
				   {
					 switch(op->code) {
					 case 'VISI':
					   I->AtomInfo[a].visRep[op->i1]=op->i2;
					   break;
					 case 'COLR':
					   I->AtomInfo[a].color=op->i1;
					   break;
					 case 'TTTF':
					   hit_flag=true;
					   break;
                case 'ALTR':
                  if (ok) {
                    if(!PAlterAtom(&I->AtomInfo[a],op->s1))
                      op->i1++;
                    else
                      ok=false;
                  }
                  break;
					 }
					 break;
				   }
				 s=SelectorNext(s);
			   }
			 break;
		   default: /* coord-set based properties */
			 for(b=0;b<I->NCSet;b++)
			   if(I->CSet[b])
				 {
				   inv_flag=false;
				   s=I->AtomInfo[a].selEntry;
				   while(s) 
					 {
					   if(SelectorMatch(s,sele))
						 {
						   switch(op->code) {
						   case 'SUMC':
							 a1=I->CSet[b]->AtmToIdx[a];
							 if(a1>=0)
							   {
								 add3f(op->v1,I->CSet[b]->Coord+(3*a1),op->v1);
								 op->i1++;
							   }
							 break;
						   case 'MDST': 
							 a1=I->CSet[b]->AtmToIdx[a];
							 if(a1>=0)
							   {
								 r=diff3f(op->v1,I->CSet[b]->Coord+(3*a1));
								 if(r>op->f1)
								   op->f1=r;
							   }
							 break;
						   case 'INVA': 
							 if(I->CSet[b]->AtmToIdx[a]>=0)
							   inv_flag=true;
							 break;
						   case 'VERT': 
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
							 
							 /* Moment of inertia tensor - unweighted - assumes v1 is center of molecule */
						   case 'MOME': 
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
					   s=SelectorNext(s);
					 }
				   switch(op->code) {
				   case 'INVA':
					 if(inv_flag) 
					   I->CSet[b]->fInvalidateRep(I->CSet[b],op->i1,op->i2);
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
	  case 'TTTF':
		ObjectTransformTTTf(I,op->ttt);
		break;
	  }
	}
  }
}
/*========================================================================*/
void ObjectMoleculeDescribeElement(ObjectMolecule *I,int index) 
{
  char buffer[1024];
  sprintf(buffer,"Selected %s:%s:%s:%s:%s:%s",
			 I->Obj.Name,I->AtomInfo[index].segi,I->AtomInfo[index].chain,
			 I->AtomInfo[index].resi,I->AtomInfo[index].resn,I->AtomInfo[index].name
			 );
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
	   PRINTF " ObjectMolecule: updating state %d of \"%s\".\n" , a+1, I->Obj.Name ENDF
	   I->CSet[a]->fUpdate(I->CSet[a]);
	 }
}
/*========================================================================*/
void ObjectMoleculeInvalidateRep(ObjectMolecule *I,int rep)
{
  int a;
  for(a=0;a<I->NCSet;a++) 
	 if(I->CSet[a]) {	 
		I->CSet[a]->fInvalidateRep(I->CSet[a],rep,0);
	 }
}
/*========================================================================*/
void ObjectMoleculeRender(ObjectMolecule *I,int frame,CRay *ray,Pickable **pick)
{
  if(frame<I->NCSet) {
	 I->CurCSet=frame % I->NCSet;
	 if(I->CSet[I->CurCSet]) {
		I->CSet[I->CurCSet]->fRender(I->CSet[I->CurCSet],ray,pick);
	 }
  }
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeNew(void)
{
  OOAlloc(ObjectMolecule);
  ObjectInit((Object*)I);
  I->Obj.type=cObjectMolecule;
  I->NAtom=0;
  I->NBond=0;
  I->CSet=VLAMalloc(10,sizeof(CoordSet*),5,true); /* auto-zero */
  I->NCSet=0;
  I->Bond=NULL;
  I->Obj.fRender=(void (*)(struct Object *, int, CRay *, Pickable **))ObjectMoleculeRender;
  I->Obj.fFree= (void (*)(struct Object *))ObjectMoleculeFree;
  I->Obj.fUpdate=  (void (*)(struct Object *)) ObjectMoleculeUpdate;
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectMoleculeGetNFrames;
  I->Obj.fDescribeElement = (void (*)(struct Object *,int index)) ObjectMoleculeDescribeElement;
  I->AtomInfo=VLAMalloc(10,sizeof(AtomInfoType),2,true); /* autozero here is important */
  I->CurCSet=0;
  return(I);
}
/*========================================================================*/
ObjectMolecule *ObjectMoleculeCopy(ObjectMolecule *obj)
{
  int a,b;
  int *i0,*i1;
  AtomInfoType *a0,*a1;
  OOAlloc(ObjectMolecule);
  (*I)=(*obj);
  I->CSet=VLAMalloc(I->NCSet,sizeof(CoordSet*),5,true); /* auto-zero */
  for(a=0;a<I->NCSet;a++) {
    I->CSet[a]=CoordSetCopy(obj->CSet[a]);
    I->CSet[a]->Obj=I;
  }
  I->Bond=VLAlloc(int,I->NBond*2);
  i0=I->Bond;
  i1=obj->Bond;
  for(a=0;a<I->NBond;a++) {
    *(i0++)=*(i1++);
    *(i0++)=*(i1++);
  }
  
  I->AtomInfo=VLAlloc(AtomInfoType,I->NAtom);
  a0=I->AtomInfo;
  a1=obj->AtomInfo;
  for(a=0;a<I->NAtom;a++)
    *(a0++)=*(a1++);

  for(a=0;a<I->NAtom;a++) {
    for(b=0;b<I->NAtom;b++) {
      
    }
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
		I->CSet[a]->fFree(I->CSet[a]);
		I->CSet[a]=NULL;
	 }
  VLAFreeP(I->CSet);
  VLAFreeP(I->AtomInfo);
  VLAFreeP(I->Bond);
  OOFreeP(I);
}

/*========================================================================*/
int ObjectMoleculeConnect(int **bond,AtomInfoType *ai,CoordSet *cs,float cutoff)
{
  #define cMULT 1

  int a,b,c,d,e,f,i,j;
  int a1,a2;
  float *v1,*v2,dst;
  int maxBond;
  MapType *map;
  int nBond;

  nBond = 0;
  if(cs->NIndex)
	 {
		maxBond = cs->NIndex * 8;
		(*bond) = VLAlloc(int,maxBond*2);
		
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
										
										if( (dst <= cutoff)&&(
                                                    !((ai[a1].name[0]=='H') && 
                                                      (ai[a2].name[0]=='H'))))
										  {
											 VLACheck((*bond),int,nBond*2+1);
											 (*bond)[nBond*2  ] = a1;
											 (*bond)[nBond*2+1] = a2;
											 nBond++;
										  }
									 }
								  j=MapNext(map,j);
								}
						  }
				}
			 MapFree(map);
		  }
		(*bond) = (int*)VLASetSize((*bond),nBond*2+1);

		if(DebugState&DebugMolecule) {
		  printf("ObjectMoleculeConnect: Found %d bonds.\n",nBond);
		}
	 }
  UtilSortInPlace((*bond),nBond,sizeof(int)*2,(UtilOrderFn*)BondInOrder);
  return(nBond);
}



