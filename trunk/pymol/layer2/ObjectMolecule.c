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

#include <GL/glut.h>

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

#include"Scene.h"

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

#define MAXLINELEN 255

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
CoordSet *ObjectMoleculeMOLStr2CoordSet(char *buffer,AtomInfoType **atInfoPtr)
{
  char *p;
  int nAtom,nBond,nType;
  int a,c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL;
  char cc[MAXLINELEN],resn[MAXLINELEN];
  float vdw;
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
			 switch ( atInfo[a].name[0] ) /* need to move this stuff into parameters */
				{
				case 'N' : vdw=1.8; color=NColor; break;
				case 'C' :	
				  color=CColor;
				  switch (atInfo[a].name[1]) 
					 {
					 case 'U':
						vdw=1.35; break; /* CU */
					 case 0:
					 default:
						vdw=1.8;  break; /*incl C,CL*/
					 }
				  break;
				case 'O' :  vdw=1.5;  color=OColor; break;
				case 'I' :	vdw=2.15; color=MColor; break;
				case 'P' :	vdw=1.9; color=MColor; break;
				case 'B' :	vdw=1.9; color=MColor; break; /* incl B, BR */
				case 'S' :	vdw=1.9; color=SColor; break;
				case 'F' : 
				  switch (atInfo[a].name[1])
					 {
					 case 0:
						vdw=1.35; break;
					 case 'E': 
						vdw=0.64; break;
					 default:
						vdw=1.35; break;
					 }
				  color = MColor;
				  break;
				case 'H' :
				  color=HColor;
				  vdw=1.2;
				  break;
				default:
				  color=MColor;
				  vdw=1.8;
				  break;
				}
			 atInfo[a].vdw=vdw;
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
  float vdw;
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
		  AFlag = true;
		if( *p == 'H') if(*(p+1)=='E') if(*(p+2)=='T') if(*(p+3)=='A')
		  if(*(p+4)=='T') if(*(p+5)=='M')
			 AFlag = true;
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

				  switch ( atInfo[atomCount].name[0] )
					 {
					 case 'N' : vdw=1.8; color=NColor; break;
					 case 'C' :	
						color=CColor;
						switch (atInfo[atomCount].name[1]) 
						  {
						  case 'U':
							 vdw=1.35; break; /* CU */
						  case 0:
						  default:
							 vdw=1.8;  break; /*incl C,CL*/
						  }
						break;
					 case 'O' : vdw=1.5;  color=OColor; break;
					 case 'I' :	vdw=2.15; break;
					 case 'P' :	vdw=1.9; break;
					 case 'B' :	vdw=1.9; break; /* incl B, BR */
					 case 'S' :	vdw=1.9; color=SColor; break;
					 case 'F' : 
						switch (atInfo[atomCount].name[1])
						  {
						  case 0:
							 vdw=1.35; break;
						  case 'E': 
							 vdw=0.64; break;
						  default:
							 vdw=1.35; break;
						  }
						break;
					 case 'H' :
						color=HColor;
						vdw=1.2;
						break;
					 default:
						color=HColor;
						vdw=1.8;
						break;
					 }
				  atInfo[atomCount].vdw=vdw;
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

	 if(isNew) {		
		for(a=0;a<3;a++)
		  I->Obj.center[a] = (I->Obj.extent[2*a]+I->Obj.extent[2*a+1])/2.0;
		I->AtomInfo=atInfo; /* IMPORTANT: this is a VLA whose address may have changed! */
	 } else {
		VLAFreeP(atInfo); /* for now, we're assuming idential molecular structures in PDB frames */
	 }
  }
  
  if(ok) {
	 if(isNew) I->NAtom=nAtom;

	 if(frame<0) frame=I->NCSet;
	 VLACheck(I->CSet,CoordSet*,frame);
	 if(I->NCSet<=frame) I->NCSet=frame+1;
	 I->CSet[frame] = cset;


	 cset->Obj = I;
	 cset->fEnumIndices(cset);
	 cset->fInvalidateRep(cset,-1,0);
	 if(isNew) ObjectMoleculeConnect(I,cset,0.2);

	 SceneCountFrames();
	 ObjectMoleculeExtendIndices(I);
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
void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op)
{
  int a,b,s;
  int a1;
  float r;
  int inv_flag;

  if(sele>=0) {
	 switch(op->code) {
	 default:
		for(a=0;a<I->NAtom;a++)
		  {
			 switch(op->code) {
			 case 'COLR':
			 case 'VISI':
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
								  }
								  break;
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
		printf(" ObjectMolecule: updating \"%s\", frame %i.\n",I->Obj.Name,a+1);
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
  I->NCyl=0;
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
void ObjectMoleculeConnect(ObjectMolecule *I,CoordSet *cs,float cutoff)
{
  #define cMULT 1


  int a,b,c,d,e,f,i,j;
  int a1,a2;
  float *v1,*v2,dst;
  int maxBond;
  MapType *map;

  I->NBond = 0;
  if(cs->NIndex)
	 {
		maxBond = cs->NIndex * 8;
		I->Bond = VLAlloc(int,maxBond*2);
		
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
										
										dst -= ((I->AtomInfo[a1].vdw+I->AtomInfo[a2].vdw)/2);
										
										if( dst <= cutoff ) 
										  {
											 VLACheck(I->Bond,int,I->NBond*2+1);
											 I->Bond[I->NBond*2  ] = i;
											 I->Bond[I->NBond*2+1] = j;
											 I->NBond++;
										  }
									 }
								  j=MapNext(map,j);
								}
						  }
				}
			 MapFree(map);
		  }
		I->Bond = (int*)VLASetSize(I->Bond,I->NBond*2+1);

		if(DebugState&DebugMolecule) {
		  printf("ObjectMoleculeConnect: Found %d bonds.\n",I->NBond);
		}
	 }
}



