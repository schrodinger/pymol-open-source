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

#include"AtomInfo.h"
#include"Word.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Util.h"
#include"Color.h"


static int NColor,CarbColor,HColor,OColor,SColor,MColor,IColor;

int AtomInfoInOrder(AtomInfoType *atom,int atom1,int atom2);

void AtomInfoPrimeColors(void)
{
  NColor=ColorGetIndex("nitrogen");
  CarbColor=ColorGetIndex("carbon");
  HColor=ColorGetIndex("hydrogen");
  OColor=ColorGetIndex("oxygen");
  SColor=ColorGetIndex("sulfer");
  MColor=ColorGetIndex("magenta");
  IColor=ColorGetIndex("yellow");
}

float AtomInfoGetBondLength(AtomInfoType *ai1,AtomInfoType *ai2)
{
  float result = 1.6;

  /* very simple for now just CC and CH, 
     flush this out with decent parameters later */

  if((ai1->protons>1)&&(ai2->protons>1)) {
    result=0.0;
    switch(ai1->geom) {
    case cAtomInfoLinear: result+=1.20; break;
    case cAtomInfoPlaner: result+=1.34; break;
    default: result+= 1.54; break;
    }
    switch(ai2->geom) {
    case cAtomInfoLinear: result+=1.20; break;
    case cAtomInfoPlaner: result+=1.34; break;
    default: result+= 1.54; break;
    }    
    result/=2.0;
  } else {
    result=1.05;
  }
  return(result);
}

int AtomInfoGetColor(AtomInfoType *at1)
{
  char *n=at1->elem;
  int color = 0;

  while(((*n)>='0')&&((*n)<='9')&&(*(n+1))) n++;
  switch ( (*n) )
    {
    case 'N' : color = NColor; break;
    case 'C' : 
      n++;
      switch(*(n)) {
      case 'A':
      case 'a':
        if(at1->hetatm) 
          color = IColor;
        else
          color = CarbColor;
      case 0:
      case 32:
        color = CarbColor; 
        break;
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
  return(color);
}


int *AtomInfoGetSortedIndex(AtomInfoType *rec,int n,int **outdex)
{
  int *index;
  int a;
  index = Alloc(int,n+1);
  ErrChkPtr(index);
  (*outdex)=Alloc(int,n);
  ErrChkPtr(*outdex);

  UtilSortIndex(n,rec,index,(UtilOrderFn*)AtomInfoInOrder);

  for(a=0;a<n;a++)
	(*outdex)[index[a]]=a;
  return(index);
}

void AtomInfoFreeSortedIndexes(int *index,int *outdex)
{
  FreeP(index);
  FreeP(outdex);
}

int AtomInfoCompare(AtomInfoType *at1,AtomInfoType *at2)
{
  int result;
  int wc;
  /* order by segment, chain, residue value, residue id, residue name, priority,
	 and lastly by name */

  /*  printf(":%s:%s:%d:%s:%s:%i:%s =? ",
	  at1->segi,at1->chain,at1->resv,at1->resi,at1->resn,at1->priority,at1->name);
	  printf(":%s:%s:%d:%s:%s:%i:%s\n",
	  at2->segi,at2->chain,at2->resv,at2->resi,at2->resn,at2->priority,at2->name);
  */

  wc=WordCompare(at1->segi,at2->segi,true);
  if(!wc) {
	 if(at1->chain[0]==at2->chain[0]) {
      if(at1->hetatm==at2->hetatm) {
        if(at1->resv==at2->resv) {
          wc=WordCompare(at1->resi,at2->resi,true);
          if(!wc) {
            wc=WordCompare(at1->resn,at2->resn,true);
            if(!wc) {
              if(at1->alt[0]==at2->alt[0]) {
                if(at1->priority==at2->priority) {
                  result=WordCompare(at1->name,at2->name,true);
                } else if(at1->priority<at2->priority) {
                  result=-1;
                } else {
                  result=1;
                }
              } else if((!at2->alt[0])||(at1->alt[0]&&((at1->alt[0]<at2->alt[0])))) {
                result=-1;
              } else {
                result=1;
              }
            } else {
              result=wc;
            }
          } else {
            result=wc;
          }
        } else if(at1->resv<at2->resv) {
          result=-1;
        } else {
          result=1;
        }
      } else if(at2->hetatm) {
        result=-1;
      } else {
        result=1;
      }
    } else if((!at2->chain[0])||(at1->chain[0]&&((at1->chain[0]<at2->chain[0])))) {
      result=-1;
    } else {
      result=1;
    }
  } else {
	 result=wc;
  }
  return(result);
}

int AtomInfoInOrder(AtomInfoType *atom,int atom1,int atom2)
{
  return(AtomInfoCompare(atom+atom1,atom+atom2)<=0);
}

int AtomInfoSameResidue(AtomInfoType *at1,AtomInfoType *at2)
{
  int result;
  if(at1->hetatm==at2->hetatm)
    if(at1->chain[0]==at2->chain[0])
      if(WordMatch(at1->resi,at2->resi,true)<0)
        if(WordMatch(at1->segi,at2->segi,true)<0) {
          result = true;
        } else {
          result = false;
        }
  return(result);
}

int AtomInfoMatch(AtomInfoType *at1,AtomInfoType *at2)
{
  if((tolower(at1->chain[0]))==(tolower(at2->chain[0])))
	 if(WordMatch(at1->name,at2->name,true)<0)
		if(WordMatch(at1->resi,at2->resi,true)<0)
		  if(WordMatch(at1->resn,at2->resn,true)<0)
			 if(WordMatch(at1->segi,at2->segi,true)<0)
            if((tolower(at1->alt[0]))==(tolower(at2->alt[0])))
              return 1;
  return 0;
}

int AtomInfoAltMatch(AtomInfoType *at1,AtomInfoType *at2)
{
  if((tolower(at1->chain[0]))==(tolower(at2->chain[0])))
    if(WordMatch(at1->resi,at2->resi,true)<0)
      if(WordMatch(at1->resn,at2->resn,true)<0)
        if(WordMatch(at1->segi,at2->segi,true)<0)
            return 1;
  return 0;
}

int AtomInfoGetExpectedValence(AtomInfoType *I) {
  int result=-1; /* negative indicates minimum expected valence (abs)
                but it could be higher  */

  if(I->formalCharge==0) {
    switch(I->protons) {
    case cAN_H : result= 1; break;
    case cAN_C : result= 4; break;
    case cAN_N : result= 3; break;
    case cAN_O : result= 2; break;
    case cAN_F : result= 1; break;
    case cAN_Cl: result= 1; break;
    case cAN_Br: result= 1; break;
    case cAN_I : result= 1; break;
    case cAN_Na: result= 1; break; 
    case cAN_Ca: result= 1; break; 
    case cAN_K : result= 1; break;
    case cAN_Mg: result= 2; break;
    case cAN_Zn: result=-1; break;
    case cAN_S : result=-2; break;
    case cAN_P : result=-3; break;
    }
  } else if(I->formalCharge==1) {
    switch(I->protons) {
    case cAN_N : result= 4; break;
    case cAN_O : result= 3; break;
    case cAN_Na: result= 0; break; 
    case cAN_Ca: result= 0; break; 
    case cAN_K : result= 0; break;
    case cAN_Mg: result= 1; break; 
    case cAN_Zn: result=-1; break;
    case cAN_S : result=-2; break;
    case cAN_P : result=-3; break;
    }
  } else if(I->formalCharge==-1) {
    switch(I->protons) {
    case cAN_N : result= 2; break;
    case cAN_O : result= 1; break;
    case cAN_C : result= 3; break;
    case cAN_Zn: result=-1; break;
    case cAN_S : result=-2; break;
    case cAN_P : result=-3; break;
    }
  } else if(I->formalCharge==2) {
    switch(I->protons) {
    case cAN_Mg: result= 0; break;
    case cAN_Zn: result=-1; break;
    case cAN_S : result=-2; break;
    case cAN_P : result=-3; break;
    }
  }
  return(result);
}
void AtomInfoAssignParameters(AtomInfoType *I)
{
  char *n,*e;
  int pri;
  float vdw;

  e = I->elem;
  if(!*e) {
    if(!*e) {
      n = I->name;
      while(((*n)>='0')&&((*n)<='9')&&(*(n+1))) n++;  
      while((((*n)>='A')&&((*n)<='Z'))||(((*n)>='a')&&((*n)<='z'))) {
        *(e++)=*(n++);
      }
      *e=0;
    }
    e=I->elem;
    switch ( *e ) {
    case 'C':
      if(!((*(e+1)=='l')||(*(e+1)=='L')))
        if(!(((*(e+1))=='A')&&(I->hetatm)))
          *(e+1)=0;
      break;
    default:
      *(e+1)=0;
      break;
    }
  }
  I->hydrogen=(((*I->elem)=='H')&&(!(*(I->elem+1))));
  n = I->name;
  while((*n>='0')&&(*n<='0')&&(*(n+1))) n++;
  switch ( *n )
	 {
	 case 'N' : 
	 case 'C' : 
	 case 'O' :
	 case 'S' :
		switch (*(n+1))
		  {
		  case 0: 
			 switch ( *n )
				{
				case 'N':
				  pri = 1; break;
				case 'C':
				  pri = 997; break;
				case 'O':
				  pri = 998; break;
				default:
				  pri = 1000; break;
				}
			 break;
		  case 'A': pri = 3; break;
		  case 'B': pri = 4; break;
		  case 'G': pri = 5; break;
		  case 'D': pri = 6; break;
		  case 'E': pri = 7; break;
		  case 'Z': pri = 8; break;
		  case 'H': pri = 9; break;
		  case 'I': pri = 10; break;
		  case 'J': pri = 11; break;
		  case 'K': pri = 12; break;
		  case 'L': pri = 13; break;
		  case 'M': pri = 14; break;
		  case 'N': pri = 15; break;
		  case '0':
		  case '1':
		  case '2':
		  case '3':
		  case '4':
		  case '5':
		  case '6':
		  case '7':
		  case '8':
		  case '9':
			pri=0;
			n++;
			while(*n) {
			  pri*=10;
			  pri+=(*n-'0');
			  n++;
			}
			pri+=25;
			break;
		  default: pri = 500; break;
		  }
		break;
	 default:
		pri = 1000; break;
	 }
  I->priority=pri;
  n = I->name;  
  while((*n>='0')&&(*n<='9')&&(*(n+1))) n++;
  switch ( *n )
    {
    case 'N' : vdw=1.8;  break;
    case 'C' :	
      switch (*(n+1)) 
        {
        case 'U':
          vdw=1.35; break; /* CU */
        case 0:
        default:
          vdw=1.8;  break; /*incl C,CL*/
        }
      break;
    case 'O' : vdw=1.5;  break;
    case 'I' :	vdw=2.15; break;
    case 'P' :	vdw=1.9; break;
    case 'B' :	vdw=1.9; break; /* incl B, BR */
    case 'S' :	vdw=1.9; break;
    case 'F' : 
      switch (*(n+1))
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
      vdw = 1.1;
      break;
    default:
      vdw=1.8;
      break;
    }

  e = I->elem;
  if(!e[1]) { /* single letter */
    switch(e[0]) {
    case 'H': I->protons=cAN_H; break;
    case 'C': I->protons=cAN_C; break;
    case 'N': I->protons=cAN_N; break;
    case 'O': I->protons=cAN_O; break;
    case 'F': I->protons=cAN_F; break;
    case 'S': I->protons=cAN_S; break;
    case 'P': I->protons=cAN_P; break;
    case 'K': I->protons=cAN_K; break;
    case 'I': I->protons=cAN_I; break;
    }
  } else {
    switch(e[0]) {
    case 'C':
      switch(e[1]) {
      case 'L':
      case 'l': I->protons=cAN_Cl; break;
      case 'A':
      case 'a': I->protons=cAN_Ca; break;
      }
      break;
    case 'B':
      switch(e[1]) {
      case 'R':
      case 'r': I->protons=cAN_Br; break;
      }
      break;
    case 'M':
      switch(e[1]) {
      case 'G':
      case 'g': I->protons=cAN_Mg; break;
      }
      break;
    case 'Z':
      switch(e[1]) {
      case 'N':
      case 'n': I->protons=cAN_Zn; break;
      }
      break;
    }
  }
  if(I->vdw==0.0)
    I->vdw = vdw;
}












