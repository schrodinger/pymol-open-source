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

#include<ctype.h>

#include"AtomInfo.h"
#include"Word.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Util.h"


int AtomInfoInOrder(AtomInfoType *atom,int atom1,int atom2);

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
		if(at1->resv==at2->resv) {
		  wc=WordCompare(at1->resi,at2->resi,true);
		  if(!wc) {
			 wc=WordCompare(at1->resn,at2->resn,true);
			 if(!wc) {
				if(at1->priority==at2->priority) {
				  result=WordCompare(at1->name,at2->name,true);
				} else if(at1->priority<at2->priority) {
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
	 } else if(at1->chain[0]<at2->chain[0]) {
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

int AtomInfoMatch(AtomInfoType *at1,AtomInfoType *at2)
{
  if((tolower(at1->chain[0]))==(tolower(at2->chain[0])))
	 if(WordMatch(at1->name,at2->name,true)<0)
		if(WordMatch(at1->resi,at2->resi,true)<0)
		  if(WordMatch(at1->resn,at2->resn,true)<0)
			 if(WordMatch(at1->segi,at2->segi,true)<0)
				return 1;
  return 0;
}

void AtomInfoAssignParameters(AtomInfoType *I)
{
  char *n;
  int pri;
  float vdw;
  
  n = I->name;
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
		  break;
		}
	  break;
	default:
	  vdw=1.8;
	  break;
	}
  I->vdw = vdw;
}

