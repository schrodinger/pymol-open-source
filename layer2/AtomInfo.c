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

#include"AtomInfo.h"
#include"Word.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Util.h"
#include"Color.h"
#include"PConv.h"


static int NColor,CarbColor,HColor,OColor,SColor,MColor,IColor;


/*========================================================================*/

int AtomInfoInOrder(AtomInfoType *atom,int atom1,int atom2);

/*========================================================================*/
PyObject *AtomInfoAsPyList(AtomInfoType *I)
{
  PyObject *result = NULL;

  result = PyList_New(34);
  PyList_SetItem(result, 0,PyInt_FromLong(I->resv));
  PyList_SetItem(result, 1,PyString_FromString(I->chain));
  PyList_SetItem(result, 2,PyString_FromString(I->alt));
  PyList_SetItem(result, 3,PyString_FromString(I->resi));
  PyList_SetItem(result, 4,PyString_FromString(I->segi));
  PyList_SetItem(result, 5,PyString_FromString(I->resn));
  PyList_SetItem(result, 6,PyString_FromString(I->name));
  PyList_SetItem(result, 7,PyString_FromString(I->elem));
  PyList_SetItem(result, 8,PyString_FromString(I->textType));
  PyList_SetItem(result, 9,PyString_FromString(I->label));
  PyList_SetItem(result,10,PyString_FromString(I->ssType));
  PyList_SetItem(result,11,PyInt_FromLong(I->hydrogen));  
  PyList_SetItem(result,12,PyInt_FromLong(I->customType));
  PyList_SetItem(result,13,PyInt_FromLong(I->priority));
  PyList_SetItem(result,14,PyFloat_FromDouble(I->b));
  PyList_SetItem(result,15,PyFloat_FromDouble(I->q));
  PyList_SetItem(result,16,PyFloat_FromDouble(I->vdw));
  PyList_SetItem(result,17,PyFloat_FromDouble(I->partialCharge));
  PyList_SetItem(result,18,PyInt_FromLong(I->formalCharge));
  PyList_SetItem(result,19,PyInt_FromLong((int)I->hetatm));
  PyList_SetItem(result,20,PConvSIntArrayToPyList(I->visRep,cRepCnt));
  PyList_SetItem(result,21,PyInt_FromLong(I->color));
  PyList_SetItem(result,22,PyInt_FromLong(I->id));
  PyList_SetItem(result,23,PyInt_FromLong(I->cartoon));
  PyList_SetItem(result,24,PyInt_FromLong(I->flags));
  PyList_SetItem(result,25,PyInt_FromLong((int)I->bonded));
  PyList_SetItem(result,26,PyInt_FromLong((int)I->chemFlag));
  PyList_SetItem(result,27,PyInt_FromLong((int)I->geom));
  PyList_SetItem(result,28,PyInt_FromLong((int)I->valence));
  PyList_SetItem(result,29,PyInt_FromLong((int)I->masked));
  PyList_SetItem(result,30,PyInt_FromLong((int)I->protected));
  PyList_SetItem(result,31,PyInt_FromLong((int)I->protons));
  PyList_SetItem(result,32,PyInt_FromLong(I->sculpt_id));
  PyList_SetItem(result,33,PyInt_FromLong(I->stereo));

  return(PConvAutoNone(result));
}

int AtomInfoFromPyList(AtomInfoType *I,PyObject *list)
{
  int ok=true;
  int hetatm;
  if(ok) ok = PyList_Check(list);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list, 0),&I->resv);
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 1),I->chain,sizeof(Chain));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 2),I->alt,sizeof(Chain));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 3),I->resi,sizeof(ResIdent));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 4),I->segi,sizeof(SegIdent));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 5),I->resn,sizeof(ResName));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 6),I->name,sizeof(AtomName));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 7),I->elem,sizeof(AtomName));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 8),I->textType,sizeof(TextType));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list, 9),I->label,sizeof(LabelType));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,10),I->ssType,sizeof(SSType));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,11),&I->hydrogen);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,12),&I->customType);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,13),&I->priority);
  if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,14),&I->b);
  if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,15),&I->q);
  if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,16),&I->vdw);
  if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,17),&I->partialCharge);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,18),&I->formalCharge);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,19),&hetatm);
  if(ok) I->hetatm = hetatm;
  if(ok) ok = PConvPyListToSIntArrayInPlaceAutoZero(PyList_GetItem(list,20),I->visRep,cRepCnt);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,21),&I->color);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,22),&I->id);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,23),&I->cartoon);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,24),(int*)&I->flags);
  if(ok) ok = PConvPyIntToChar(PyList_GetItem(list,25),(char*)&I->bonded);
  if(ok) ok = PConvPyIntToChar(PyList_GetItem(list,26),(char*)&I->chemFlag);
  if(ok) ok = PConvPyIntToChar(PyList_GetItem(list,27),(char*)&I->geom);
  if(ok) ok = PConvPyIntToChar(PyList_GetItem(list,28),(char*)&I->valence);
  if(ok) ok = PConvPyIntToChar(PyList_GetItem(list,29),(char*)&I->masked);
  if(ok) ok = PConvPyIntToChar(PyList_GetItem(list,30),(char*)&I->protected);
  if(ok) ok = PConvPyIntToChar(PyList_GetItem(list,31),(char*)&I->protons);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,32),&I->sculpt_id);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,33),&I->stereo);
  return(ok);
}

/*========================================================================*/
void AtomInfoCombine(AtomInfoType *dst,AtomInfoType *src,int mask)
{
  if(mask&cAIC_tt) strcpy(dst->textType,src->textType); /* use the new types */
  if(mask&cAIC_ct) dst->customType = src->customType;
  if(mask&cAIC_pc) dst->partialCharge = src->partialCharge;
  if(mask&cAIC_fc) dst->formalCharge = src->formalCharge;
  if(mask&cAIC_flags) dst->flags = src->flags;
  if(mask&cAIC_b) dst->b = src->b;
  if(mask&cAIC_q) dst->q = src->q;
  if(mask&cAIC_id) dst->id = src->id;
  dst->temp1 = src->temp1;
  dst->sculpt_id = src->sculpt_id;
  /* keep all existing names, identifiers, etc. */
  /* also keep all existing selections,
     colors, masks, and visible representations*/

}
/*========================================================================*/
void AtomInfoUniquefyNames(AtomInfoType *atInfo0,int n0,AtomInfoType *atInfo1,int n1)
{
  /* makes sure all names in atInfo1 are unique WRT 0 and 1 */

  /* tricky optimizations to avoid n^2 dependence in this operation */

  int a,b,c;
  AtomInfoType *ai0,*ai1,*lai1,*lai0;
  int st1,nd1,st0,nd0; /* starts and ends */
  int matchFlag;
  int bracketFlag;
  WordType name;

  ai1=atInfo1;
  lai0=NULL; /* last atom compared against in each object */
  lai1=NULL;
  c=1;
  /* ai1->name is the atom we're currently on */

  b=0;
  while(b<n1) {
    matchFlag=false;

    if(!ai1->name[0])
      matchFlag=true;

    if(!matchFlag) {
      /* check within object 1 */
      
      if(!lai1) bracketFlag=true;
      else if(!AtomInfoSameResidue(lai1,ai1))
        bracketFlag=true;
      else
        bracketFlag=false;
      if(bracketFlag) {
        c=1;
        AtomInfoBracketResidue(atInfo1,n1,ai1,&st1,&nd1);
        lai1=ai1;
      }
      ai0 = atInfo1 + st1;
      for(a=st1;a<=nd1;a++) {
        if(strcmp(ai1->name,ai0->name))
          ai0++;
        else if(!AtomInfoSameResidue(ai1,ai0))
          ai0++;
        else if(ai1!=ai0) {
          matchFlag=true;
          break;
        } else 
          ai0++;
      }
    }

    if(!matchFlag) {
      if(atInfo0) {
        /* check within object 2 */
        
        if(!lai0) bracketFlag=true;
        else if(!AtomInfoSameResidue(lai0,ai1))
          bracketFlag=true;
        else
          bracketFlag=false;
        if(bracketFlag) {
          AtomInfoBracketResidue(atInfo0,n0,ai1,&st0,&nd0);
          lai0=ai1;
        }
        ai0 = atInfo0 + st0;
        for(a=st0;a<=nd0;a++) {
          if(strcmp(ai1->name,ai0->name))
            ai0++;
          else if(!AtomInfoSameResidue(ai1,ai0))
            ai0++;
          else if(ai1!=ai0) {
            matchFlag=true;
            break;
          } else 
            ai0++;
        }
      }
    }

    if(matchFlag) {
      if(c<100) {
        if((c<10)&&ai1->elem[1]) /* try to keep halogens 3 or under */
          sprintf(name,"%2s%1d",ai1->elem,c);
        else 
          sprintf(name,"%1s%02d",ai1->elem,c);
      } else {
        sprintf(name,"%1d%1s%02d",c/100,ai1->elem,c%100);
      }
      name[4]=0; /* just is case we go over */
      strcpy(ai1->name,name);
      c=c+1;
    } else {
      ai1++;
      b++;
    }
  }
}

/*========================================================================*/
void AtomInfoBracketResidue(AtomInfoType *ai0,int n0,AtomInfoType *ai,int *st,int *nd)
{
  /* inefficient but reliable way to find where residue atoms are located in an object 
   * for purpose of residue-based operations */
  int a;
  AtomInfoType *ai1;
  
  *st=0;
  *nd=n0-1;
  ai1=ai0;
  for(a=0;a<n0;a++) {
    if(!AtomInfoSameResidue(ai,ai1++))
      *st=a;
    else 
      break;
  }
  ai1=ai0+n0-1;
  for(a=n0-1;a>=0;a--) {
    if(!AtomInfoSameResidue(ai,ai1--)) {
      *nd=a;
    } else 
      break;
  }
}
/*========================================================================*/
void AtomInfoBracketResidueFast(AtomInfoType *ai0,int n0,int cur,int *st,int *nd)
{
  /* efficient but unreliable way to find where residue atoms are located in an object 
   * for purpose of residue-based operations */
  int a;
  AtomInfoType *ai1;
  
  *st=cur;
  *nd=cur;
  ai0=ai0+cur;
  ai1=ai0-1;
  for(a=cur-1;a>=0;a--) {
    if(!AtomInfoSameResidue(ai0,ai1--))
      break;
    *st=a;
  }
  ai1=ai0+1;
  for(a=cur+1;a<n0;a++) {
    if(!AtomInfoSameResidue(ai0,ai1++))
      break;
    *nd=a;
  }
}

/*========================================================================*/
void AtomInfoPrimeColors(void)
{
  NColor=ColorGetIndex("nitrogen");
  if(SettingGet(cSetting_auto_color))
    CarbColor=ColorGetNext();
  else
    CarbColor=ColorGetIndex("carbon");
  HColor=ColorGetIndex("hydrogen");
  OColor=ColorGetIndex("oxygen");
  SColor=ColorGetIndex("sulfur");
  MColor=ColorGetIndex("magenta");
  IColor=ColorGetIndex("yellow");
}
/*========================================================================*/
float AtomInfoGetBondLength(AtomInfoType *ai1,AtomInfoType *ai2)
{
  float result = 1.6F;
  AtomInfoType *a1,*a2;
  
  /* very simple for now ... 
     flush this out with pattern-based parameters later on */
  
  if(ai1->protons>ai2->protons) {
    a1 = ai2;
    a2 = ai1;
  } else {
    a1 = ai1;
    a2 = ai2;
  }    
  switch(a1->protons) {
    
/* hydrogen =========================== */
  case cAN_H:
    switch(a2->protons) {
    case cAN_H: result = 0.74F; break;
    case cAN_C: result = 1.09F; break;
    case cAN_N: result = 1.01F; break;
    case cAN_S: result = 1.34F; break;
    case cAN_O: result = 0.96F; break;
    default:
      result = 1.09F; 
      break;
    }
    break;
/* carbon =========================== */

  case cAN_C: /* carbon */
    switch(a1->geom) {
      
    case cAtomInfoLinear: /* linear carbon =============*/
      switch(a2->geom) {
      case cAtomInfoLinear: 
        switch(a2->protons) {
        case cAN_C: result = 1.20F; break; /* C#^C */
        case cAN_N: result = 1.16F; break; /* C#^N */
        default: result = 1.20F; break;
        }
        break;
      case cAtomInfoPlaner: 
        switch(a2->protons) { 
        case cAN_I:  result = 2.14F; break; /* normal single bond lengths */
        case cAN_Cl: result = 1.77F; break; 
        case cAN_Br: result = 1.94F; break; 
        case cAN_F:  result = 1.35F; break; 
        case cAN_S:  result = 1.82F; break; 
        case cAN_N:  result = 1.47F; break; 
        case cAN_O:  result = 1.43F; break; 
        case cAN_P:  result = 1.84F; break;
        case cAN_C:  result = 1.54F; break;
        default: result = 1.54F; break;
        }
        break;
      default: /* ?#C-^? */
        switch(a2->protons) { 
        case cAN_I:  result = 2.14F; break; /* normal single bond lengths */
        case cAN_Cl: result = 1.77F; break; 
        case cAN_Br: result = 1.94F; break; 
        case cAN_F:  result = 1.35F; break; 
        case cAN_S:  result = 1.82F; break; 
        case cAN_N:  result = 1.47F; break; 
        case cAN_O:  result = 1.43F; break; 
        case cAN_P:  result = 1.84F; break;
        case cAN_C:  result = 1.54F; break;
        default: result = 1.54F; break;
        }
        break;
      }
      break;
      
    case cAtomInfoPlaner: /* planer carbon =============*/
      switch(a2->geom) {
      case cAtomInfoLinear: 
        switch(a2->protons) { 
        case cAN_I:  result = 2.14F; break; /* normal single bond lengths */
        case cAN_Cl: result = 1.77F; break; 
        case cAN_Br: result = 1.94F; break; 
        case cAN_F:  result = 1.35F; break; 
        case cAN_S:  result = 1.82F; break; 
        case cAN_N:  result = 1.47F; break; 
        case cAN_O:  result = 1.43F; break; 
        case cAN_P:  result = 1.84F; break;
        case cAN_C:  result = 1.54F; break;
        default: result = 1.54F; break;
        }
        break;
      case cAtomInfoPlaner: 
        switch(a2->protons) {
        case cAN_C: result = 1.34F; break; /* C=^C or ?=C-^C=? */
        case cAN_O: result = 1.20F; break; /* carbonyl */
        case cAN_N: result = 1.29F; break; /* C=^N or ?=C-^N=? */
        case cAN_S: result = 1.60F; break; /* C=^S or ?=C-^S-?=? */
        default: result = 1.34F; break;
        }
        break;
      default: /* ?#C-^? */
        switch(a2->protons) { 
        case cAN_I:  result = 2.14F; break; /* normal single bond lengths */
        case cAN_Cl: result = 1.77F; break; 
        case cAN_Br: result = 1.94F; break; 
        case cAN_F:  result = 1.35F; break; 
        case cAN_S:  result = 1.71F; break; /* mmod */
        case cAN_N:  result = 1.47F; break; 
        case cAN_O:  result = 1.43F; break; 
        case cAN_P:  result = 1.84F; break;
        case cAN_C:  result = 1.54F; break;
        default: result = 1.54F; break;
        }
        break;
      }
      break;

    default: /* tetrahedral carbon */
      switch(a2->protons) { 
      case cAN_I:  result = 2.14F; break; /* normal single bond lengths */
      case cAN_Cl: result = 1.77F; break; 
      case cAN_Br: result = 1.94F; break; 
      case cAN_F:  result = 1.35F; break; 
      case cAN_S:  result = 1.82F; break; 
      case cAN_N:  result = 1.47F; break; 
      case cAN_O:  result = 1.43F; break; 
      case cAN_P:  result = 1.84F; break;
      case cAN_C:  result = 1.54F; break;
      default: result = 1.54F; break;
      }
      break;
    }
    break;

/* nitrogen ========================= */

  case cAN_N:
    if((a1->geom==cAtomInfoPlaner)&&
       (a2->geom==cAtomInfoPlaner))
      switch(a2->protons) {
      case cAN_N:  result = 1.25F; break;
      case cAN_O:  result = 1.21F; break;
      case cAN_S:  result = 1.53F; break; /* interpolated */
      default: result = 1.25F; break;
      }
    else {
      switch(a2->protons) {
      case cAN_N:  result = 1.45F; break;
      case cAN_O:  result = 1.40F; break;
      case cAN_S:  result = 1.75F; break; /* interpolated */
      default: result = 1.45F; break;
      }
    }
    break;
    
/* oxygen =========================== */

  case cAN_O:
    if((a1->geom==cAtomInfoPlaner)&&
       (a2->geom==cAtomInfoPlaner))
      switch(a2->protons) {
      case cAN_O:  result = 1.35F; break; /* guess */
      case cAN_S:  result = 1.44F; break; /* macromodel */
      default: result = 1.35F; break;
      }
    else if(a1->geom==cAtomInfoPlaner) {
      switch(a2->protons) {
      case cAN_O:  result = 1.35F; break; /* guess */
      case cAN_S:  result = 1.44F; break; /* macromodel */
      default: result = 1.35F; break;
      }
    } else {
      switch(a2->protons) {
      case cAN_O:  result = 1.40F; break;
      case cAN_S:  result = 1.75F; break; /* interpolated */
      default: result = 1.45F; break;
      }
    }
    break;

/* sulfur =========================== */

  case cAN_S:
    switch(a2->protons) {
    case cAN_S:  result = 2.05F; break; /* interpolated */
    default: result = 1.82F; break;
    }
    break;
    
    /* fall-back to old method */
  default:

    result=0.0;
    switch(a1->geom) {
    case cAtomInfoLinear: result+=1.20F; break;
    case cAtomInfoPlaner: result+=1.34F; break;
    default: result+= 1.54F; break;
    }
    switch(a2->geom) {
    case cAtomInfoLinear: result+=1.20F; break;
    case cAtomInfoPlaner: result+=1.34F; break;
    default: result+= 1.54F; break;
    }    
    result/=2.0F;
    break;
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
        break;
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

int AtomInfoNameOrder(AtomInfoType *at1,AtomInfoType *at2) 
{
  int result;
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
  return(result);
}

int AtomInfoInOrder(AtomInfoType *atom,int atom1,int atom2)
{
  return(AtomInfoCompare(atom+atom1,atom+atom2)<=0);
}

int AtomInfoSameResidue(AtomInfoType *at1,AtomInfoType *at2)
{
  if(at1->hetatm==at2->hetatm)
    if(at1->chain[0]==at2->chain[0])
      if(at1->resv==at2->resv)
        if(WordMatch(at1->resi,at2->resi,true)<0) 
          if(WordMatch(at1->segi,at2->segi,true)<0) 
            return 1;
  return 0;
}

int AtomInfoSequential(AtomInfoType *at1,AtomInfoType *at2)
{
  char last1=0,last2=0;
  char *p;
  if(at1->hetatm==at2->hetatm)
    if(at1->chain[0]==at2->chain[0]) {
      if(WordMatch(at1->segi,at2->segi,true)<0) {
        if(at1->resv==at2->resv) {
          p=at1->resi;
          while(*p) {
            last1=(*p++);
          }
          p=at2->resi;
          while(*p) {
            last2=(*p++);
          }
          if(last1==last2)
            return 1;
          if((last1+1)==last2)
            return 1;
        } else {
          if((at1->resv+1)==at2->resv)
            return 1;
        }
      }
    }
  return 0;
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
  if(!*e) { /* try to guess atomic type from name */
    n = I->name;
    while(((*n)>='0')&&((*n)<='9')&&(*(n+1))) n++;  
    while((((*n)>='A')&&((*n)<='Z'))||(((*n)>='a')&&((*n)<='z'))) {
      *(e++)=*(n++);
    }
    *e=0;
    e=I->elem;
    switch ( *e ) {
    case 'C':
      if(*(e+1)=='A') {
        if(!(WordMatch("CA",I->resn,true)<0))
          *(e+1)=0; 
      } else if(!(
           (*(e+1)=='a')||/* CA intpreted as carbon, not calcium */
           (*(e+1)=='l')||(*(e+1)=='L')||
           (*(e+1)=='u')||(*(e+1)=='U')||
           (*(e+1)=='o')||(*(e+1)=='O')||
           (*(e+1)=='s')||(*(e+1)=='S')||
           (*(e+1)=='r')||(*(e+1)=='R')
           ))
        *(e+1)=0; 
      break;
    case 'H':
      if(!((*(e+1)=='e')
           ))
        *(e+1)=0;
      break;
    case 'N':
      if(!(
           (*(e+1)=='i')||(*(e+1)=='I')||
           (*(e+1)=='a')||(*(e+1)=='A')||
           (*(e+1)=='b')||(*(e+1)=='B')
           ))
        *(e+1)=0; 
      break;
    case 'S':
      if(!(
           (*(e+1)=='e')||(*(e+1)=='E')||
           (*(e+1)=='r')||(*(e+1)=='R')||
           (*(e+1)=='c')||(*(e+1)=='C')||
           (*(e+1)=='b')||(*(e+1)=='B')
           ))
        *(e+1)=0; 
      
      break;
    case 'O':
      if(!((*(e+1)=='s')))
        *(e+1)=0;
      break;
    default:
      break;
    }
    if(*(e+1)) *(e+1)=tolower(*(e+1));
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
  e = I->elem;  
  while((*e>='0')&&(*e<='9')&&(*(e+1))) e++;
  switch ( *e )
    {
    case 'N' : vdw=1.55F;  break;
    case 'C' :	
      switch (*(e+1)) 
        {
        case 'U':
        case 'u':
          vdw=1.40F; break; /* Cu */
        case 'L':
        case 'l':
          vdw=1.75F; break; /* Cl */
        case 0:
          vdw=1.7F;  break; /* Carbon */
        default:
          vdw=1.8F;  break; 
        }
      break;
    case 'O' : vdw=1.52F;  break;
    case 'I' :	vdw=1.98F; break;
    case 'P' :	vdw=1.80F; break;
    case 'B' :	vdw=1.85F; break; /* incl B, BR */
    case 'S' :	vdw=1.80F; break;
    case 'F' : 
      switch (*(e+1))
        {
        case 0:
          vdw=1.47F; break;
        case 'E': 
        case 'e': 
          /*          vdw=0.64F; break;*/
          vdw=1.80F; break; /* default */
        default:
          vdw=1.35F; break;
        }
      break;
    case 'H' :
      vdw = 1.2F; /* WLD */
      break;
    default:
      vdw=1.8F;
      break;
    }
  if(SettingGet(cSetting_legacy_vdw_radii)) { /* ver<0.75, old, incorrect VDW */
    if(!strcmp(e,"N")) vdw=1.8F; /* slow but compact */
    if(!strcmp(e,"C")) vdw=1.8F;
    if(!strcmp(e,"Cl")) vdw=1.8F;
    if(!strcmp(e,"O")) vdw=1.5F;
    if(!strcmp(e,"Br")) vdw=1.9F;
    if(!strcmp(e,"I")) vdw=2.15F;
    if(!strcmp(e,"S")) vdw=1.9F;
    if(!strcmp(e,"P")) vdw=1.9F;
    if(!strcmp(e,"F")) vdw=1.35F;
    if(!strcmp(e,"H")) vdw=1.1F;
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
      case 'U':
      case 'u': I->protons=cAN_Cu; break;
      }
      break;
    case 'B':
      switch(e[1]) {
      case 'R':
      case 'r': I->protons=cAN_Br; break;
      }
      break;
    case 'F':
      switch(e[1]) {
      case 'E':
      case 'e': I->protons=cAN_Fe; break;
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
    case 'S':
      switch(e[1]) {
      case 'I':
      case 'i': I->protons=cAN_Si; break;
      }
      break;
    }
  }
  if(I->vdw==0.0)
    I->vdw = vdw;
}












