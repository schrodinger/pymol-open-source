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
#include"PUtils.h"
#include"Executive.h"
#include"Setting.h"


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
void ObjectMoleculeExtendIndices(ObjectMolecule *I);
void ObjectMoleculeMerge(ObjectMolecule *I,AtomInfoType *ai,CoordSet *cs);

int ObjectMoleculeConnect(int **bond,AtomInfoType *ai,CoordSet *cs,float cutoff);
void ObjectMoleculeTransformTTTf(ObjectMolecule *I,float *ttt,int state);

static int BondInOrder(int *a,int b1,int b2);
static int BondCompare(int *a,int *b);

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
    fflush(stdout);
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
void ObjectMoleculeSort(ObjectMolecule *I)
{
  int *index,*outdex;
  int a,b;
  CoordSet *cs;
  AtomInfoType *atInfo;
  
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
  int ok=true;
  int autoshow_lines;

  autoshow_lines = SettingGet(cSetting_autoshow_lines);
  AtomInfoPrimeColors();

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
          p=nskip(p,1);
			 p=ncopy(atInfo[a].name,p,3);
			 UtilCleanStr(atInfo[a].name);
			 
			 atInfo[a].visRep[0] = autoshow_lines; /* show lines by default */
			 for(c=1;c<cRepCnt;c++) {
				atInfo[a].visRep[c] = false;
			 }
		  }
		  if(ok&&atInfo) {
			 strcpy(atInfo[a].resn,resn);
			 atInfo[a].hetatm=true;
			 AtomInfoAssignParameters(atInfo+a);
			 atInfo[a].color=AtomInfoGetColor(atInfo+a);
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
				  p=nskip(p,16);
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
				  p=nskip(p,8);
				  p=ncopy(cc,p,3);
				  if(sscanf(cc,"%d",&atInfo[a].customType)!=1) {
					 ok=ErrMessage("ObjectMoleculeLoadMOLFile","bad custom atom type record-1");
				  } else {
					 p=nskip(p,1);
					 p=ncopy(cc,p,3);
					 if(sscanf(cc,"%d",&atInfo[a].customFlag)!=1) {
						ok=ErrMessage("ObjectMoleculeLoadMOLFile","bad custom atom type record-2");
					 } else {
					   p=nskip(p,1);
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
      if(cset->fAppendIndices)
        cset->fAppendIndices(cset,I->NAtom);
		cset->Obj=I;
		ObjectMoleculeAppendAtoms(I,atInfo,cset);
      if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
		I->CSet[frame] = cset;
      if(I->CSet[frame]->fInvalidateRep)
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
        
        VLACheck(I->Bond,int,nBd*3);
        
        for(a=0;a<nBond;a++)
          {
            a2=index[a];
            if(a2 >= I->NBond) { 
              I->Bond[3*a2]=bond[3*a]; /* copy bond info */
              I->Bond[3*a2+1]=bond[3*a+1]; /* copy bond info */
              I->Bond[3*a2+2]=bond[3*a+2]; /* copy bond info */
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
		I=(ObjectMolecule*)ObjectMoleculeNew();
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
    cset->fEnumIndices(cset);
    if(cset->fInvalidateRep)
      cset->fInvalidateRep(cset,-1,0);
    if(isNew) {		
      I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
    } else {
      ObjectMoleculeMerge(I,atInfo,cset); /* NOTE: will release atInfo */
    }
    cset->Obj = I;
    if(isNew) I->NAtom=nAtom;
    if(frame<0) frame=I->NCSet;
    VLACheck(I->CSet,CoordSet*,frame);
    if(I->NCSet<=frame) I->NCSet=frame+1;
    if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
    I->CSet[frame] = cset;
    if(isNew) I->NBond = ObjectMoleculeConnect(&I->Bond,I->AtomInfo,cset,0.2);
    if(cset->TmpSymmetry&&(!I->Symmetry)) {
      I->Symmetry=cset->TmpSymmetry;
      cset->TmpSymmetry=NULL;
      SymmetryAttemptGeneration(I->Symmetry);
    }
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
           cs->fInvalidateRep(I->CSet[b],cRepAll,0);
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
                    ind=I->CSet[b]->AtmToIdx[a];
						  if(ind>=0) 
                      CoordSetAtomToPDBStrVLA(&op->charVLA,&op->i2,I->AtomInfo+a,
                                              I->CSet[b]->Coord+(3*ind),op->i3);
						  op->i3++;
						}
					  s=SelectorNext(s);
					}
				}
		  }
	  break;
	case 'AVRT': /* average vertex coordinate */
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
				   while(s) 
                 {
                   if(SelectorMatch(s,sele))
                     {
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
                   s=SelectorNext(s);
                 }
             }
           if(maxCnt<op->nvv1) maxCnt=op->nvv1;
         }
     }
     op->nvv1=maxCnt;
     break;
	case 'SFIT': /* state fitting within a single object */
	  for(b=0;b<I->NCSet;b++) {
       rms = -1.0;
       if(I->CSet[b]&&(b!=op->i1))
         {
           op->nvv1=0;
           for(a=0;a<I->NAtom;a++)
             {
				   s=I->AtomInfo[a].selEntry;
				   while(s) 
					 {
					   if(SelectorMatch(s,sele))
                    {
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
                  s=SelectorNext(s);
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
             printf(" Executive: RMS = %8.3f (%d atoms)\n",
                    rms,op->nvv1);
             if(op->i1==2) 
               ObjectMoleculeTransformTTTf(I,op->ttt,b);
           } else {
             ErrMessage("ExecutiveFit","No atoms selected.");
           }
         }
       VLACheck(op->f1VLA,float,b);
       op->f1VLA[b]=rms;
     }
     VLASetSize(op->f1VLA,I->NCSet); 
     break;

	default:
	   for(a=0;a<I->NAtom;a++)
		 {
		   switch(op->code) { 
		   case 'COLR': /* atom based loops */
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
                  if(op->i1<0)
                    for(d=0;d<cRepCnt;d++) 
                      I->AtomInfo[a].visRep[d]=op->i2;                      
                  else
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
#ifdef PYMOL_FUTURE_CODE
         case 'CSOC': /* specific coordinate set based operations */
           if(I->NCSet<op->cs1) 
             if(I->CSet[op->cs1]) {
               
               s=I->AtomInfo[a].selEntry;
               while(s)
                 {
                   if(SelectorMatch(s,sele))
                     {
                       switch(op->code) {
                       case 'CSOC': /* object and coordinate index */
                         break;
                       }
                       break;
                     }
                   s=SelectorNext(s);
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
						   case 'MNMX': 
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
						   case 'SVRT':  /* gives us only vertices for a specific coordinate set */
                       if(b==op->i1) {
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
                 if(inv_flag) {
                  if(op->i1<0)
                    for(d=0;d<cRepCnt;d++) {
                      if(I->CSet[b]->fInvalidateRep)
                        I->CSet[b]->fInvalidateRep(I->CSet[b],d,op->i2);
                    }
                  else if(I->CSet[b]->fInvalidateRep)
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
	  case 'TTTF':
		ObjectMoleculeTransformTTTf(I,op->ttt,-1);
		break;
	  }
	}
  }
}
/*========================================================================*/
void ObjectMoleculeDescribeElement(ObjectMolecule *I,int index) 
{
  char buffer[1024];
  sprintf(buffer," Pick: Selected %s:%s:%s:%s:%s:%s",
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
	   printf(" ObjectMolecule: updating state %d of \"%s\".\n" , a+1, I->Obj.Name);
      if(I->CSet[a]->fUpdate)
        I->CSet[a]->fUpdate(I->CSet[a]);
	 }
}
/*========================================================================*/
void ObjectMoleculeInvalidateRep(ObjectMolecule *I,int rep)
{
  int a;
  for(a=0;a<I->NCSet;a++) 
	 if(I->CSet[a]) {	 
      if(I->CSet[a]->fInvalidateRep)
        I->CSet[a]->fInvalidateRep(I->CSet[a],rep,0);
	 }
}
/*========================================================================*/
void ObjectMoleculeRender(ObjectMolecule *I,int frame,CRay *ray,Pickable **pick)
{
  int a;
  if(frame<0) {
    for(a=0;a<I->NCSet;a++)
      if(I->CSet[a])
        if(I->CSet[a]->fRender)
          I->CSet[a]->fRender(I->CSet[a],ray,pick);        
  } else if(frame<I->NCSet) {
	 I->CurCSet=frame % I->NCSet;
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
  I->Symmetry=NULL;
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
  int nBond,*ii1,*ii2;
  nBond = 0;
  if(cs->NIndex)
	 {
		maxBond = cs->NIndex * 8;
		(*bond) = VLAlloc(int,maxBond*3);
		
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
											 VLACheck((*bond),int,nBond*3+2);
											 (*bond)[nBond*3  ] = a1;
											 (*bond)[nBond*3+1] = a2;
                                  (*bond)[nBond*3+2] = 1;
											 nBond++;
										  }
									 }
								  j=MapNext(map,j);
								}
						  }
				}
			 MapFree(map);
		  }
		(*bond) = (int*)VLASetSize((*bond),nBond*3+2);

		if(DebugState&DebugMolecule) 
		  printf("ObjectMoleculeConnect: Found %d bonds.\n",nBond);
		
	 }
  if(cs->NTmpBond&&cs->TmpBond) {
    if(DebugState&DebugMolecule) 
      printf("ObjectMoleculeConnect: incorporating CONECT bonds. %d\n",nBond);
    VLACheck((*bond),int,(nBond+cs->NTmpBond)*3);
    ii1=(*bond)+nBond*3;
    ii2=cs->TmpBond;
    for(a=0;a<cs->NTmpBond;a++)
      {
        *(ii1++)=*(ii2++);
        *(ii1++)=*(ii2++);
        *(ii1++)=*(ii2++);
      }
    nBond=nBond+cs->NTmpBond;
    VLAFreeP(cs->TmpBond);
    cs->NTmpBond=0;
  }
  
  UtilSortInPlace((*bond),nBond,sizeof(int)*3,(UtilOrderFn*)BondInOrder);
  return(nBond);
}



/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadMMDStr(ObjectMolecule *I,char *MMDStr,int frame)
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
		  I=(ObjectMolecule*)ObjectMoleculeNew();
		if(frame<0)
		  frame=I->NCSet;
		if(I->NCSet<=frame)
		  I->NCSet=frame+1;
		VLACheck(I->CSet,CoordSet*,frame);
      nAtom=cset->NIndex;

      if(cset->fEnumIndices)
        cset->fEnumIndices(cset);
      if(cset->fInvalidateRep)
        cset->fInvalidateRep(cset,-1,0);
      if(isNew) {		
        I->AtomInfo=atInfo; /* IMPORTANT to reassign: this VLA may have moved! */
        I->NAtom=nAtom;
      } else {
        ObjectMoleculeMerge(I,atInfo,cset); /* NOTE: will release atInfo */
      }
      cset->Obj = I;
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
ObjectMolecule *ObjectMoleculeLoadMMDFile(ObjectMolecule *obj,char *fname,
                                          int frame,char *sepPrefix)
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
            I=ObjectMoleculeReadMMDStr(NULL,p,frame);
            oCnt++;
            sprintf(oName,"%s-%02d",sepPrefix,oCnt);
            ObjectSetName((Object*)I,oName);
            ExecutiveManageObject((Object*)I);
          } else {
            I=ObjectMoleculeReadMMDStr(obj,p,frame);
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
  int autoshow_lines;
  CSymmetry *symmetry = NULL;
  int symFlag;
  autoshow_lines = SettingGet(cSetting_autoshow_lines);

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
  for(a=0;a<255;a++) /*to prevent hopping over end of file*/
	 *p++=0;
  
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
  fflush(stdout);
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
          if(!sscanf(cc,"%d",&ai->tmpID)) ai->tmpID=0;

          p=nskip(p,1);/* to 12 */
          p=ncopy(cc,p,4); 
          if(!sscanf(cc,"%s",ai->name)) ai->name[0]=0;
          
          p=nskip(p,1);
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

          p=ncopy(cc,p,5);
          if(!sscanf(cc,"%s",ai->resi)) ai->resi[0]=0;
          if(!sscanf(cc,"%d",&ai->resv)) ai->resv=0;
          
          p=nskip(p,3);
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+a);
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+(a+1));
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+(a+2));

          p=ncopy(cc,p,6);
          sscanf(cc,"%f",&ai->q);
          
          p=ncopy(cc,p,6);
          sscanf(cc,"%f",&ai->b);

          p=nskip(p,3);
          p=ncopy(cc,p,4);
          if(!sscanf(cc,"%s",ai->segi)) ai->segi[0]=0;
          
          ai->visRep[0] = autoshow_lines;
          for(c=1;c<cRepCnt;c++) {
            ai->visRep[c] = false;
          }

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
        if(maxAt<atInfo[a].tmpID) maxAt=atInfo[a].tmpID;
      /* build index */
      idx = Alloc(int,maxAt+1);
      for(a=0;a<maxAt;a++) idx[a]=-1;
      for(a=0;a<nAtom;a++)
        idx[atInfo[a].tmpID]=a;
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
  int autoshow_lines;

  autoshow_lines = SettingGet(cSetting_autoshow_lines);
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

        if(ok) {
          p=ncopy(cc,p,4);
          if(sscanf(cc,"%d",&ai->customType)!=1)
            ok=ErrMessage("ReadMMDFile","bad atom type");
        }

        for(c=0;c<6;c++) {
          if(ok) {
            p=ncopy(cc,p,8);
            if(sscanf(cc,"%d%d",&bPart,&bOrder)!=2)
              ok=ErrMessage("ReadMMDFile","bad bond record");
            else {
              if(bPart&&bOrder) {
                nBond++;
                *(ii++)=a;
                *(ii++)=bPart-1;
                *(ii++)=1;
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
          p=nskip(p,31);
          p=ncopy(cc,p,3);
          if(sscanf(cc,"%s",ai->resn)!=1)
            ai->resn[0]=0;
          ai->hetatm=true;
        }

        ai->segi[0]=0;

        if(ok) {
          p=nskip(p,2);
          p=ncopy(ai->name,p,4);
          UtilCleanStr(ai->name);
          if(ai->name[0]==0) {
            if(ai->customType<=14) strcpy(ai->name,"C");
            else if(ai->customType<=23) strcpy(ai->name,"O");
            else if(ai->customType<=40) strcpy(ai->name,"N");
            else if(ai->customType<=48) strcpy(ai->name,"H");
            else if(ai->customType<=52) strcpy(ai->name,"S");
            else if(ai->customType<=53) strcpy(ai->name,"P");
            else if(ai->customType<=55) strcpy(ai->name,"B");
            else if(ai->customType<=56) strcpy(ai->name,"F");
            else if(ai->customType<=57) strcpy(ai->name,"Cl");           
            else if(ai->customType<=58) strcpy(ai->name,"Br");           
            else if(ai->customType<=59) strcpy(ai->name,"I");           
            else if(ai->customType<=60) strcpy(ai->name,"Si");           
            else if(ai->customType<=61) strcpy(ai->name,"Du");           
            else if(ai->customType<=62) strcpy(ai->name,"Z0");
            else if(ai->customType<=63) strcpy(ai->name,"Lp");
            else strcpy(ai->name,"?");
            sprintf(cc,"%02d",a+1);
            if((strlen(cc)+strlen(ai->name))>4)
              strcpy(ai->name,cc);
            else
              strcat(ai->name,cc);
          }
          ai->visRep[0] = autoshow_lines; /* show lines by default */
          for(c=1;c<cRepCnt;c++) {
            ai->visRep[c] = false;
          }
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
    bond=VLASetSize(bond,3*sizeof(int)*nBond);
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
