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
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Scene.h"
#include"CoordSet.h"
#include"Color.h"
#include"PConv.h"
#include"P.h"
#include"ButMode.h"
#include"Matrix.h"
#include"Sphere.h"
#include"Util.h"
#include"Feedback.h"
#include"RepWireBond.h"
#include"RepCylBond.h"
#include"RepDot.h"
#include"RepMesh.h"
#include"RepSphere.h"
#include"RepRibbon.h"
#include"RepCartoon.h"
#include"RepSurface.h"
#include"RepLabel.h"
#include"RepNonbonded.h"
#include"RepNonbondedSphere.h"

#include"PyMOLGlobals.h"

static void CoordSetUpdate(CoordSet *I,int state);

void CoordSetFree(CoordSet *I);
void CoordSetRender(CoordSet *I,RenderInfo *info);
void CoordSetEnumIndices(CoordSet *I);
void CoordSetStrip(CoordSet *I);
void CoordSetInvalidateRep(CoordSet *I,int type,int level);
void CoordSetExtendIndices(CoordSet *I,int nAtom);
void CoordSetAppendIndices(CoordSet *I,int offset);

/*========================================================================*/
static  char sATOM[]="ATOM  ";
static  char sHETATM[]="HETATM";

/*========================================================================*/
int BondCompare(BondType *a,BondType *b)
{
  register int ai0=a->index[0];
  register int bi0=b->index[0];
  if(ai0==bi0) {
    register int ai1=a->index[1];
    register int bi1=b->index[1];
    if(ai1==bi1) {
      return 0;
    } else if(ai1>bi1) {
      return 1;
    } else {
      return -1;
    }
  } else if(ai0>bi0) {
    return 1;
  } else {
    return -1;
  }
}

/*========================================================================*/
int BondInOrder(BondType *a,int b1,int b2)
{
  return(BondCompare(a+b1,a+b2)<=0);
}

int CoordSetFromPyList(PyMOLGlobals *G,PyObject *list,CoordSet **cs)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  CoordSet *I = NULL;
  PyObject *tmp;
  int ok = true;
  int ll = 0;

  if(*cs) {
    CoordSetFree(*cs);
    *cs=NULL;
  }

  if(list==Py_None) { /* allow None for CSet */
    *cs = NULL;
  } else {
  
    if(ok) I=CoordSetNew(G);
    if(ok) ok = (I!=NULL);
    if(ok) ok = (list!=NULL);
    if(ok) ok = PyList_Check(list);
    if(ok) ll = PyList_Size(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->NIndex);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NAtIndex);
    if(ok) ok = PConvPyListToFloatVLA(PyList_GetItem(list,2),&I->Coord);
    if(ok) ok = PConvPyListToIntArray(PyList_GetItem(list,3),&I->IdxToAtm);
    if(ok) {
      tmp = PyList_GetItem(list,4); /* Discrete CSets don't have this */
      if(tmp!=Py_None) 
        ok = PConvPyListToIntArray(tmp,&I->AtmToIdx);
    }
    if(ok&&(ll>5)) ok = PConvPyStrToStr(PyList_GetItem(list,5),I->Name,sizeof(WordType));
    if(ok&&(ll>6)) ok = ObjectStateFromPyList(G,PyList_GetItem(list,6),&I->State);
    if(ok&&(ll>7)) I->Setting = SettingNewFromPyList(G,PyList_GetItem(list,7));
    if(ok&&(ll>8)) ok = PConvPyListToLabPosVLA(PyList_GetItem(list,8),&I->LabPos);
    if(!ok) {
      if(I)
        CoordSetFree(I);
    } else {
      *cs = I;
    }
  }
  return(ok);
#endif
}

PyObject *CoordSetAsPyList(CoordSet *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  if(I) {
    result = PyList_New(9);
    
    PyList_SetItem(result,0,PyInt_FromLong(I->NIndex));
    PyList_SetItem(result,1,PyInt_FromLong(I->NAtIndex));
    PyList_SetItem(result,2,PConvFloatArrayToPyList(I->Coord,I->NIndex*3));
    PyList_SetItem(result,3,PConvIntArrayToPyList(I->IdxToAtm,I->NIndex));
    if(I->AtmToIdx) 
      PyList_SetItem(result,4,PConvIntArrayToPyList(I->AtmToIdx,I->NAtIndex));
    else 
      PyList_SetItem(result,4,PConvAutoNone(NULL));
    PyList_SetItem(result,5,PyString_FromString(I->Name));
    PyList_SetItem(result,6,ObjectStateAsPyList(&I->State));
    PyList_SetItem(result,7,SettingAsPyList(I->Setting));
    PyList_SetItem(result,8,PConvLabPosVLAToPyList(I->LabPos,I->NIndex));
    /* TODO symmetry, spheroid, periodic box ... */
  }
  return(PConvAutoNone(result));
#endif

}

void CoordSetAdjustAtmIdx(CoordSet *I,int *lookup,int nAtom)
     /* performs second half of removal */
{
  /* NOTE: only works in a compressive mode, where lookup[a]<=a  */
  int a;
  int a0;

  PRINTFD(I->State.G,FB_CoordSet)
    " CoordSetAdjustAtmIdx-Debug: entered NAtIndex: %d NIndex %d\n I->AtmToIdx %p\n",
    I->NAtIndex,I->NIndex,(void*)I->AtmToIdx
    ENDFD;

  for(a=0;a<I->NAtIndex;a++) {
    a0=lookup[a];
    if(a0>=0) {
      I->AtmToIdx[a0] = I->AtmToIdx[a];
    }
  }
  I->NAtIndex = nAtom;
  I->AtmToIdx = Realloc(I->AtmToIdx,int,nAtom);
  for(a=0;a<I->NIndex;a++) { 
    I->IdxToAtm[a] = lookup[I->IdxToAtm[a]];
  }
  PRINTFD(I->State.G,FB_CoordSet)
    " CoordSetAdjustAtmIdx-Debug: leaving... NAtIndex: %d NIndex %d\n",
    I->NAtIndex,I->NIndex
    ENDFD;

}
/*========================================================================*/
void CoordSetMerge(CoordSet *I,CoordSet *cs) /* must be non-overlapping */
{
  int nIndex;
  int a,i0;

  nIndex = I->NIndex+cs->NIndex;
  I->IdxToAtm=Realloc(I->IdxToAtm,int,nIndex);
  VLACheck(I->Coord,float,nIndex*3);
  for(a=0;a<cs->NIndex;a++) {
    i0 = a+I->NIndex;
    I->IdxToAtm[i0] = cs->IdxToAtm[a];
    I->AtmToIdx[cs->IdxToAtm[a]] = i0;
    copy3f(cs->Coord+a*3,I->Coord+i0*3);
  }
  if(cs->LabPos) {
    if(!I->LabPos) 
      I->LabPos = VLACalloc(LabPosType,I->NIndex);
    if(I->LabPos) {
      UtilCopyMem(I->LabPos+I->NIndex,cs->LabPos,sizeof(LabPosType)*cs->NIndex);
    }
  }
  if(I->fInvalidateRep)
    I->fInvalidateRep(I,cRepAll,cRepInvAll);
  I->NIndex = nIndex;
}
/*========================================================================*/
void CoordSetPurge(CoordSet *I) 
     /* performs first half of removal  */
{
  int offset = 0;
  int a,a1,ao;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  float *c0,*c1;
  LabPosType *l0,*l1;
  obj=I->Obj;

  PRINTFD(I->State.G,FB_CoordSet)
    " CoordSetPurge-Debug: entering..."
    ENDFD;

  c0 = c1 = I->Coord;
  l0 = l1 = I->LabPos;

  for(a=0;a<I->NIndex;a++) {
    a1 = I->IdxToAtm[a];
    ai = obj->AtomInfo + a1;
    if(ai->deleteFlag) {
      offset--;
      c0+=3;
      if(l0)
        l0++;
    } else if(offset) {
        ao=a+offset;
        *(c1++)=*(c0++);
        *(c1++)=*(c0++);
        *(c1++)=*(c0++);
        if(l0) {
          *(l1++) = *(l0++);
        }
        I->AtmToIdx[a1] = ao;
        I->IdxToAtm[ao] = a1; /* no adjustment of these indexes yet...*/
    } else {
      c0+=3;
      c1+=3;
      if(l0) {
        l0++;
        l1++;
      }
    }
  }
  if(offset) {
    I->NIndex+=offset;
    VLASize(I->Coord,float,I->NIndex*3);
    if(I->LabPos)
      VLASize(I->LabPos,LabPosType,I->NIndex);
    I->IdxToAtm=Realloc(I->IdxToAtm,int,I->NIndex);
    PRINTFD(I->State.G,FB_CoordSet)
      " CoordSetPurge-Debug: I->IdxToAtm shrunk to %d\n",I->NIndex
      ENDFD;
    if(I->fInvalidateRep)
      I->fInvalidateRep(I,cRepAll,cRepInvAtoms); /* this will free Color */
  }
  PRINTFD(I->State.G,FB_CoordSet)
    " CoordSetPurge-Debug: leaving NAtIndex %d NIndex %d...\n",
    I->NAtIndex,I->NIndex
    ENDFD;

}
/*========================================================================*/
int CoordSetTransformAtomTTTf(CoordSet *I,int at,float *TTT)
{
  ObjectMolecule *obj;
  int a1 = -1;
  int result = 0;
  float *v1;

  obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(I==obj->DiscreteCSet[at])
      a1=obj->DiscreteAtmToIdx[at];
  } else 
    a1=I->AtmToIdx[at];
  
  if(a1>=0) {
    result = 1;
    v1 = I->Coord+3*a1;
    MatrixTransformTTTfN3f(1,v1,TTT,v1);
  }

  return(result);
}

/*========================================================================*/
int CoordSetTransformAtomR44f(CoordSet *I,int at,float *matrix)
{
  ObjectMolecule *obj;
  int a1 = -1;
  int result = 0;
  float *v1;

  obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(I==obj->DiscreteCSet[at])
      a1=obj->DiscreteAtmToIdx[at];
  } else 
    a1=I->AtmToIdx[at];
  
  if(a1>=0) {
    result = 1;
    v1 = I->Coord+3*a1;
    MatrixTransformR44fN3f(1,v1,matrix,v1);
  }

  return(result);
}
/*========================================================================*/
void CoordSetRecordTxfApplied(CoordSet *I,float *matrix,int homogenous)
{
  if(I->State.Matrix) { 
    double temp[16];
    if(!homogenous) {
      convertTTTfR44d(matrix,temp);
    } else {
      convert44f44d(matrix,temp);
    }
    left_multiply44d44d(temp,I->State.Matrix);
  } else {
    I->State.Matrix=Alloc(double,16);
    if ( I->State.Matrix ) {
      if(!homogenous)
        convertTTTfR44d(matrix,I->State.Matrix);
      else {
        convert44f44d(matrix,I->State.Matrix);
      }
    }
  }
  /*  dump44d(I->State.Matrix,"history");*/
}
/*========================================================================*/
int CoordSetMoveAtom(CoordSet *I,int at,float *v,int mode)
{
  ObjectMolecule *obj;
  int a1 = -1;
  int result = 0;
  float *v1;

  obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(I==obj->DiscreteCSet[at])
      a1=obj->DiscreteAtmToIdx[at];
  } else 
    a1=I->AtmToIdx[at];

  if(a1>=0) {
    result = 1;
    v1 = I->Coord+3*a1;
    if(mode) {
      add3f(v,v1,v1);
    } else {
      copy3f(v,v1);
    }
  }

  return(result);
}
/*========================================================================*/
int CoordSetMoveAtomLabel(CoordSet *I,int at,float *v,int mode)
{
  ObjectMolecule *obj;
  int a1 = -1;
  int result = 0;
  LabPosType *lp;

  obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(I==obj->DiscreteCSet[at])
      a1=obj->DiscreteAtmToIdx[at];
  } else 
    a1=I->AtmToIdx[at];
  
  if(a1>=0) {
    if(!I->LabPos) 
      I->LabPos = VLACalloc(LabPosType,I->NIndex);
    if(I->LabPos) {
      result = 1;
      lp = I->LabPos+a1;
      if(!lp->mode) {
        float *lab_pos = SettingGet_3fv(obj->Obj.G,I->Setting,obj->Obj.Setting,cSetting_label_position);
        copy3f(lab_pos,lp->pos);
      }
      lp->mode=1;
      if(mode) {
        add3f(v,lp->offset,lp->offset);
      } else { 
        copy3f(v,lp->offset);
      }
    }
  }

  return(result);
}
/*========================================================================*/
int CoordSetGetAtomVertex(CoordSet *I,int at,float *v)
{
  register ObjectMolecule *obj;
  register int a1 = -1;
  register int result = 0;

  obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(I==obj->DiscreteCSet[at])
      a1=obj->DiscreteAtmToIdx[at];
  } else 
    a1=I->AtmToIdx[at];
  
  if(a1>=0) {
    result = 1;
    copy3f(I->Coord+3*a1,v);
  }

  return(result);
}
/*========================================================================*/
int CoordSetGetAtomTxfVertex(CoordSet *I,int at,float *v)
{
  register ObjectMolecule *obj;
  register int a1 = -1;
  register int result = 0;
  
  obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(I==obj->DiscreteCSet[at])
      a1=obj->DiscreteAtmToIdx[at];
  } else 
    a1=I->AtmToIdx[at];
  
  if(a1>=0) {
    result = 1;
    copy3f(I->Coord+3*a1,v);
    if(I->State.Matrix && SettingGet_b(I->State.G,
                                       obj->Obj.Setting,I->Setting,
                                       cSetting_matrix_mode)>0) {
      /* apply state transformation */
      transform44d3f(I->State.Matrix,v,v);
    }
    if(obj->Obj.TTTFlag) { /* object transformation */
      transformTTT44f3f(obj->Obj.TTT, v,v);
    }
  }
  return(result);
}

/*========================================================================*/
int CoordSetSetAtomVertex(CoordSet *I,int at,float *v)
{
  ObjectMolecule *obj;
  int a1 = -1;
  int result = 0;

  obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(I==obj->DiscreteCSet[at])
      a1=obj->DiscreteAtmToIdx[at];
  } else 
    a1=I->AtmToIdx[at];
  
  if(a1>=0) {
    result = 1;
    copy3f(v,I->Coord+3*a1);
  }

  return(result);
}

/*========================================================================*/
void CoordSetRealToFrac(CoordSet *I,CCrystal *cryst)
{
  int a;
  float *v;
  v=I->Coord;
  for(a=0;a<I->NIndex;a++) {
    transform33f3f(cryst->RealToFrac,v,v);
    v+=3;
  }
}
/*========================================================================*/
void CoordSetTransform44f(CoordSet *I,float *mat)
{
  int a;
  float *v;
  v=I->Coord;
  for(a=0;a<I->NIndex;a++) {
    transform44f3f(mat,v,v);
    v+=3;
  }  
}
/*========================================================================*/

void CoordSetTransform33f(CoordSet *I,float *mat)
{
  int a;
  float *v;
  v=I->Coord;
  for(a=0;a<I->NIndex;a++) {
    transform33f3f(mat,v,v);
    v+=3;
  }  
}
/*========================================================================*/
void CoordSetGetAverage(CoordSet *I,float *v0)
{
  int a;
  float *v;
  double accum[3];
  if(I->NIndex) {
    v=I->Coord;
    accum[0]=*(v++);
    accum[1]=*(v++);
    accum[2]=*(v++);
    for(a=1;a<I->NIndex;a++) {
      accum[0]+=*(v++);
      accum[1]+=*(v++);
      accum[2]+=*(v++);
    }
    v0[0]=(float)(accum[0]/I->NIndex);
    v0[1]=(float)(accum[1]/I->NIndex);
    v0[2]=(float)(accum[2]/I->NIndex);
  }
}
/*========================================================================*/
void CoordSetFracToReal(CoordSet *I,CCrystal *cryst)
{
  int a;
  float *v;
  v=I->Coord;
  for(a=0;a<I->NIndex;a++) {
    transform33f3f(cryst->FracToReal,v,v);
    v+=3;
  }
}

/*========================================================================*/
void CoordSetAtomToPDBStrVLA(PyMOLGlobals *G,char **charVLA,int *c,AtomInfoType *ai,
                             float *v,int cnt,PDBInfoRec *pdb_info)
{
  char *aType;
  AtomName name;
  ResIdent resi; 
  ResName resn;
  Chain chain;
  int rl;
  int literal = (int)SettingGet(G,cSetting_pdb_literal_names);
  int reformat = (int)SettingGet(G,cSetting_pdb_reformat_names_mode);
  WordType x,y,z;

  strcpy(resn,ai->resn); 
  if(SettingGetGlobal_b(G,cSetting_pdb_truncate_residue_name)) {
    resn[3]=0; /* enforce 3-letter residue name in PDB files */
  }

  if(ai->hetatm)
	aType=sHETATM;
  else
	aType=sATOM;

  strcpy(resi,ai->resi);
  rl = strlen(resi)-1;
  if(rl>=0)
    if((resi[rl]>='0')&&(resi[rl]<='9')) {
        resi[rl+1]=' ';
        resi[rl+2]=0;
    }
  VLACheck(*charVLA,char,(*c)+1000);  
  
  if(!ai->name[0]) {
    if(!ai->elem[1])
      sprintf(name," %s",ai->elem);
    else
      sprintf(name,"%s",ai->elem);
  } else if(!literal) {
    if(strlen(ai->name)<4) { /* atom name less than length 4 */
      if(!((ai->name[0]>='0')&&(ai->name[0])<='9')) { /* doesn't start with a number */
        if((toupper(ai->elem[0])==toupper(ai->name[0]))&&
           ((!ai->elem[1]) || /* symbol len = 1 */
            (toupper(ai->elem[1])==toupper(ai->name[1])))) { /* matched len 2 */
          /* starts with corrent atomic symbol, so */
          if(strlen(ai->elem)>1) { /* if atom symbol is length 2 */
            strcpy(name,ai->name); /* then start in column 0 */
          } else { 
            switch(reformat) {
            case 2: /* amber/iupac */
              name[0]=' ';	
              strcpy(name+1,ai->name);
              break;
            case 1: /* pdb with internal pdb */
            case 3: /* pdb with internal iupac */
              if((ai->elem[0]=='H')&&(!ai->elem[1])&&(ai->name[2])) {
                AtomInfoGetPDB3LetHydroName(G,resn,ai->name,name);
              } else {
                name[0]=' ';	
                strcpy(name+1,ai->name);
              }
              break;
            case 4:
            default: /* otherwise, start in column 1 */
              name[0]=' ';	
              strcpy(name+1,ai->name);
              break;
            } 
          }
        } else { /* name doesn't start with atomic symbol */
          /* then just place it in column 1 as usual */
          name[0]=' ';	
          strcpy(name+1,ai->name);
        }
      } else { /* name starts with a number */
        switch(reformat) {
        case 2: /* make Amber compliant */
          if((ai->elem[0]==ai->name[1]) && 
             ((!ai->elem[1]) || (toupper(ai->elem[1])==toupper(ai->name[2])))) {
            /* rotate the name to place atom symbol in column 0 to comply with Amber PDB format */
            name[0]=' ';
            name[1]=ai->name[1];
            name[2]=ai->name[2];
            name[3]=ai->name[0];
            name[4]=0;
          } else {
            strcpy(name,ai->name);
          }
          break;
        default: /* otherwise, assume that number goes in column 0 */
          strcpy(name,ai->name);
          break;
        }
      } /* just stick it in column 0 and hope for the best */
    } else { /* if name is length 4 */
      if((ai->elem[0]==ai->name[0]) &&
         ((!ai->elem[1]) || /* symbol len = 1 */
          (toupper(ai->elem[1])==toupper(ai->name[1])))) { /* matched len 2 */
        /* name starts with the atomic symbol */
        if((!ai->elem[1]) && (ai->elem[0])) { /* but if element is one letter... */
          switch(reformat) {
          case 1: /* retaining PDB compliance throughout, or */
          case 3: /* saving as PDB compliant, but use IUPAC within PyMOL */
            if((ai->name[3]>='0')&&(ai->name[3]<='9')) { /* and last character is a number */
              /* rotate the name to place atom symbol in column 1 to comply with PDB format */
              name[0]=ai->name[3];
              name[1]=ai->name[0];
              name[2]=ai->name[1];
              name[3]=ai->name[2];
              name[4]=0;
            } else {
              strcpy(name,ai->name);
            }
            break;
          case 4:
          default: /* no changes */
            strcpy(name,ai->name);
            break;
          }
        } else {
          strcpy(name,ai->name);
        }
      } else { /* name does not start with the symbol... */
        if(reformat==2) { /* AMBER compliance mode */
          if((ai->name[0]>='0')&&(ai->name[0]<='9')) {
            if((ai->elem[0]==ai->name[1]) &&
               ((!(ai->elem[1])) || 
                (toupper(ai->elem[1])==toupper(ai->name[2])))) {
              /* rotate the name to place atom symbol in column 0 to comply with Amber PDB format */
              name[0]=ai->name[1];
              name[1]=ai->name[2];
              name[2]=ai->name[3];
              name[3]=ai->name[0];
              name[4]=0;
            } else {
              strcpy(name,ai->name);
            }
          } else {
            strcpy(name,ai->name);
          }
        } else {
          strcpy(name,ai->name);
        }
      }
    }
  } else { /* LITERAL mode: preserve what was in the original PDB as best PyMOL can 
            this should enable people to open and save amber pdb files without issues */
    if(strlen(ai->name)==4) {
      strcpy(name,ai->name); /* save literal contents of field */
    } else { /* under length 4? check match with atomic symbol */
      if((toupper(ai->elem[0])==toupper(ai->name[0]))&&
         ((!ai->elem[1]) || /* symbol len = 1 */
          (toupper(ai->elem[1])==toupper(ai->name[1])))) { /* matched len 2 */
        /* starts with corrent atomic symbol, so */
        if(strlen(ai->elem)>1) { /* if atom symbol is length 2 */
          strcpy(name,ai->name); /* then start in column 0 */
        } else { 
          /* otherwise, start in column 1 */
          name[0]=' ';	
          strcpy(name+1,ai->name);
        } 
      } else {
        /* otherwise, start in column 1 */
        name[0]=' ';	
        strcpy(name+1,ai->name);
      }
    }
  }
  if((int)SettingGet(G,cSetting_pdb_retain_ids)) {
    cnt = ai->id - 1;
  }
  if(cnt>99998)
    cnt=99998;

  if((!pdb_info)||(!pdb_info->is_pqr_file)) { /* relying upon short-circuit */
    sprintf(x,"%8.3f",v[0]); x[8]=0;
    sprintf(y,"%8.3f",v[1]); y[8]=0;
    sprintf(z,"%8.3f",v[2]); z[8]=0;
    (*c)+=sprintf((*charVLA)+(*c),"%6s%5i %-4s%1s%-4s%1s%5s   %s%s%s%6.2f%6.2f      %-4s%2s\n",
                  aType,cnt+1,name,ai->alt,resn,
                  ai->chain,resi,x,y,z,ai->q,ai->b,ai->segi,ai->elem);
  } else {
    if(pdb_info->is_pqr_file && pdb_info->pqr_workarounds) {
      int non_num = false;
      char *p = resi;
      while(*p) {
        if( (((*p)<='0')||((*p)>='9'))
            &&(*p!=' ')) {
          non_num=true;
          break;
        }
        p++;
      }
      if(non_num) {
        sprintf(resi,"%d",ai->resv);
        rl = strlen(resi)-1;
        if(rl>=0)
          if((resi[rl]>='0')&&(resi[rl]<='9')) {
            resi[rl+1]=' ';
            resi[rl+2]=0;
            }
      }
      chain[0]=0;
    } else {
      chain[0] = ai->chain[0];
      chain[1] = 0;
    }
    sprintf(x,"%8.3f",v[0]); x[8]=0;
    sprintf(y,"%8.3f",v[1]); y[8]=0;
    sprintf(z,"%8.3f",v[2]); z[8]=0;
      
    (*c)+=sprintf((*charVLA)+(*c),"%6s%5i %-4s%1s%-4s%1s%5s   %s%s%s %11.8f %7.3f\n",
                  aType,cnt+1,name,ai->alt,resn,
                  chain,resi,x,y,z,ai->partialCharge,ai->elec_radius);
  }
  
}
/*========================================================================*/
PyObject *CoordSetAtomToChemPyAtom(PyMOLGlobals *G,AtomInfoType *ai,float *v,int index)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *atom;
  int ok = true;

  atom = PyObject_CallMethod(P_chempy,"Atom","");
  if (!atom) 
    ok = ErrMessage(G,"CoordSetAtomToChemPyAtom","can't create atom");
  else {
    PConvFloat3ToPyObjAttr(atom,"coord",v);
    PConvStringToPyObjAttr(atom,"name",ai->name);
    PConvStringToPyObjAttr(atom,"symbol",ai->elem);
    PConvStringToPyObjAttr(atom,"resn",ai->resn);
    PConvStringToPyObjAttr(atom,"resi",ai->resi);
    PConvStringToPyObjAttr(atom,"ss",ai->ssType);
    PConvIntToPyObjAttr(atom,"resi_number",ai->resv);
    PConvIntToPyObjAttr(atom,"stereo",ai->stereo);
    PConvStringToPyObjAttr(atom,"chain",ai->chain);
    if(ai->alt[0]) 
      PConvStringToPyObjAttr(atom,"alt",ai->alt); 
    PConvStringToPyObjAttr(atom,"segi",ai->segi);
    PConvFloatToPyObjAttr(atom,"q",ai->q);
    PConvFloatToPyObjAttr(atom,"b",ai->b);
    PConvFloatToPyObjAttr(atom,"vdw",ai->vdw);
    PConvFloatToPyObjAttr(atom,"elec_radius",ai->elec_radius);
    PConvFloatToPyObjAttr(atom,"partial_charge",ai->partialCharge);
    PConvIntToPyObjAttr(atom,"formal_charge",ai->formalCharge);
    if(ai->customType!=-9999)
      PConvIntToPyObjAttr(atom,"numeric_type",ai->customType);
    if(ai->textType) {
      char null_st[1] = "";
      char *st = null_st;
      
      if(ai->textType) st = OVLexicon_FetchCString(G->Lexicon,ai->textType);
      PConvStringToPyObjAttr(atom,"text_type",st);
    }

    PConvIntToPyObjAttr(atom,"hetatm",ai->hetatm);      
    PConvIntToPyObjAttr(atom,"flags",ai->flags);      
    PConvIntToPyObjAttr(atom,"id",ai->id);  /* not necc. unique */
    PConvIntToPyObjAttr(atom,"index",index+1);  /* fragile */
  }
  if(PyErr_Occurred())
    PyErr_Print();
  return(atom);
#endif
}
/*========================================================================*/
void CoordSetAtomToTERStrVLA(PyMOLGlobals *G,char **charVLA,int *c,AtomInfoType *ai,int cnt)
{
  ResIdent resi; 
  ResName resn;
  int rl;
  int retain_ids = (int)SettingGet(G,cSetting_pdb_retain_ids);
  int ter_id;

  strcpy(resn,ai->resn); /* enforce 3-letter residue name in PDB files */
  resn[3]=0;

  strcpy(resi,ai->resi);
  rl = strlen(resi)-1;
  if(rl>=0)
    if((resi[rl]>='0')&&(resi[rl]<='9')) {
        resi[rl+1]=' ';
        resi[rl+2]=0;
    }
  VLACheck(*charVLA,char,(*c)+1000);  

  if(retain_ids) {
    ter_id = ai->id+1;
  } else {
    ter_id = cnt+1;
  }

  (*c)+=sprintf((*charVLA)+(*c),
                "%3s   %5i      %3s %1s%5s\n",
                 "TER",ter_id,resn,ai->chain,resi);
  
}

/*========================================================================*/
void CoordSetInvalidateRep(CoordSet *I,int type,int level)
{
  int a;

  /*  printf("inv %d %d \n",type,level);fflush(stdout);*/

  if(level==cRepInvVisib) {
    if(SettingGet_b(I->State.G,I->Setting,I->Obj->Obj.Setting,
                    cSetting_cartoon_side_chain_helper)) {
      if((type==cRepCyl)||(type==cRepLine)||(type==cRepSphere))
        CoordSetInvalidateRep(I,cRepCartoon,cRepInvVisib2);
      else if(type==cRepCartoon) {
        CoordSetInvalidateRep(I,cRepLine,cRepInvVisib2);
        CoordSetInvalidateRep(I,cRepCyl,cRepInvVisib2);
        CoordSetInvalidateRep(I,cRepSphere,cRepInvVisib2);
      }
    }
    if(SettingGet_b(I->State.G,I->Setting,I->Obj->Obj.Setting,
                    cSetting_ribbon_side_chain_helper)) {
      if((type==cRepCyl)||(type==cRepLine)||(type==cRepSphere))
        CoordSetInvalidateRep(I,cRepRibbon,cRepInvVisib2);
      else if(type==cRepRibbon) {
        CoordSetInvalidateRep(I,cRepLine,cRepInvVisib2);
        CoordSetInvalidateRep(I,cRepCyl,cRepInvVisib2);
        CoordSetInvalidateRep(I,cRepSphere,cRepInvVisib2);
      }
    }
    if(SettingGet_b(I->State.G,I->Setting,I->Obj->Obj.Setting,
                    cSetting_line_stick_helper)) {
      if(type==cRepCyl) 
        CoordSetInvalidateRep(I,cRepLine,cRepInvVisib2);
      else if(type==cRepLine) {
        CoordSetInvalidateRep(I,cRepCyl,cRepInvVisib2);
      }
    }
  }

  if(I->Spheroid)
    if(I->NSpheroid!=I->NAtIndex*I->SpheroidSphereSize) {
      FreeP(I->Spheroid);
      FreeP(I->SpheroidNormal);
    }
  if(level>=cRepInvColor) 
	 VLAFreeP(I->Color);
  if(type>=0) { /* representation specific */
	 if(type<I->NRep)	{
      a=type;
      if(I->Rep[a]) {
        if(I->Rep[a]->fInvalidate && (level<cRepInvPurge))
          I->Rep[a]->fInvalidate(I->Rep[a],I,level);
        else {
          I->Rep[a]->fFree(I->Rep[a]);
          I->Rep[a] = NULL;
        }
      }
      if(level>=cRepInvVisib) /* make active if visibility has changed */
        I->Active[type]=true;
	 }
  } else { /* all representations are affected */
	 for(a=0;a<I->NRep;a++)	{
      if(level>=cRepInvVisib) /* make active if visibility has changed */
        I->Active[a]=true;
		if(I->Rep[a]) {
          if(I->Rep[a]->fInvalidate && (level<cRepInvPurge))
            I->Rep[a]->fInvalidate(I->Rep[a],I,level);
          else {
            I->Rep[a]->fFree(I->Rep[a]);
            I->Rep[a] = NULL;
          }
		}
	 }
  }
  if(level>=cRepInvCoord) { /* if coordinates change, then this map becomes invalid */
    MapFree(I->Coord2Idx);
    I->Coord2Idx = NULL;
  }
  SceneChanged(I->State.G);
}
/*========================================================================*/

#define RepUpdateMacro(I,rep,new_fn,state) {\
  if(I->Active[rep]) {\
    if(!I->Rep[rep]) {\
      I->Rep[rep]=new_fn(I,state);\
      if(I->Rep[rep]) \
         I->Rep[rep]->fNew=(struct Rep *(*)(struct CoordSet *,int state))new_fn;\
    } else {\
      if(I->Rep[rep]->fUpdate)\
         I->Rep[rep] = I->Rep[rep]->fUpdate(I->Rep[rep],I,state,rep);\
    }\
  }\
OrthoBusyFast(I->State.G,rep,I->NRep);\
}
/*========================================================================*/
static void CoordSetUpdate(CoordSet *I,int state)
{
  int a;
  int i;
  ObjectMolecule *obj;
  obj=I->Obj;
  
  if(Feedback(I->Obj->Obj.G,FB_CoordSet,FB_Blather)) {
    printf(" CoordSetUpdate-Entered: object %s state %d cset %p\n",
           I->Obj->Obj.Name, state, I);
  }

  if(!I->Color) /* colors invalidated */
	 {
		I->Color=VLAlloc(int,I->NIndex);
		if(I->Color) {
        if(obj->DiscreteFlag)
          for(a=0;a<I->Obj->NAtom;a++) {
            if(obj->DiscreteCSet[a]==I) {
              i = obj->DiscreteAtmToIdx[a];
              if(i>=0) 
                I->Color[i]=obj->AtomInfo[a].color;
            }
          }
        else 
          for(a=0;a<I->Obj->NAtom;a++)
            {
              i=I->AtmToIdx[a];
              if(i>=0) 
                I->Color[i]=obj->AtomInfo[a].color;
            }
		}
	 }
  OrthoBusyFast(I->State.G,0,I->NRep);
  RepUpdateMacro(I, cRepLine,            RepWireBondNew        , state );
  RepUpdateMacro(I, cRepCyl,             RepCylBondNew         , state );
  RepUpdateMacro(I, cRepDot,             RepDotNew             , state );
  RepUpdateMacro(I, cRepMesh,            RepMeshNew            , state );
  RepUpdateMacro(I, cRepSphere,          RepSphereNew          , state );
  RepUpdateMacro(I, cRepRibbon,          RepRibbonNew          , state );
  RepUpdateMacro(I, cRepCartoon,         RepCartoonNew         , state );
  RepUpdateMacro(I, cRepSurface,         RepSurfaceNew         , state );
  RepUpdateMacro(I, cRepLabel,           RepLabelNew           , state );
  RepUpdateMacro(I, cRepNonbonded,       RepNonbondedNew       , state );
  RepUpdateMacro(I, cRepNonbondedSphere, RepNonbondedSphereNew , state );

  for(a=0;a<I->NRep;a++) 
    if(!I->Rep[a])
      I->Active[a]=false;

  SceneInvalidate(I->State.G);
  OrthoBusyFast(I->State.G,1,1);
  if(Feedback(I->Obj->Obj.G,FB_CoordSet,FB_Blather)) {
    printf(" CoordSetUpdate-Leaving: object %s state %d cset %p\n",
           I->Obj->Obj.Name, state, I);
  }

}
/*========================================================================*/
void CoordSetUpdateCoord2IdxMap(CoordSet *I, float cutoff)
{
  if(cutoff<R_SMALL4)
    cutoff = R_SMALL4;
  if(I->NIndex>10) {
    if(I->Coord2Idx) {
      if((I->Coord2IdxDiv<cutoff)||
         (((cutoff-I->Coord2IdxReq)/I->Coord2IdxReq)<-0.5F)) {
        MapFree(I->Coord2Idx);
        I->Coord2Idx = NULL;
      }
    }
    if(I->NIndex && (!I->Coord2Idx)) { /* NOTE: map based on stored coords */
      I->Coord2IdxReq = cutoff;
      I->Coord2IdxDiv = cutoff * 1.25F;
      I->Coord2Idx = MapNew(I->State.G, I->Coord2IdxDiv, I->Coord, I->NIndex, NULL);
      if(I->Coord2IdxDiv<I->Coord2Idx->Div) 
        I->Coord2IdxDiv = I->Coord2Idx->Div;
    }
  }
}
/*========================================================================*/
void CoordSetRender(CoordSet *I,RenderInfo *info)
{
  PyMOLGlobals *G = I->State.G;
  int pass = info->pass;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int a,aa;
  Rep *r;
  int float_labels = SettingGet_i(G,I->Setting,
                                  I->Obj->Obj.Setting,
                                  cSetting_float_labels);
  PRINTFD(G,FB_CoordSet)
    " CoordSetRender: entered (%p).\n",(void*)I
    ENDFD;

  if((!pass)&&I->Name[0])
    ButModeCaption(G,I->Name);
  if((!pass)&&I->SculptCGO&&(I->Obj->Obj.RepVis[cRepCGO])) {
    if(ray) {
      CGORenderRay(I->SculptCGO,ray,
                   ColorGet(G,I->Obj->Obj.Color),
                   I->Setting,I->Obj->Obj.Setting);
    } else if(G->HaveGUI && G->ValidContext) {
      if(!pick) {
        CGORenderGL(I->SculptCGO,ColorGet(G,I->Obj->Obj.Color),
                    I->Setting,I->Obj->Obj.Setting,info);
      }
    }
  }
  for(aa=0;aa<I->NRep;aa++) {
    if(aa==cRepSurface) { /* reorder */
      a=cRepCell;
    } else if(aa==cRepCell) {
      a=cRepSurface;
    } else {
      a = aa;
    }
      
    if(I->Active[a])
      if(I->Rep[a]) 
        {
          r = I->Rep[a];
          if(!ray) {
            ObjectUseColor((CObject*)I->Obj);
          } else {
            if(I->Obj) 
              ray->fWobble(ray,
                           SettingGet_i(G,I->Setting,
                                        I->Obj->Obj.Setting,
                                        cSetting_ray_texture),
                           SettingGet_3fv(G,I->Setting,
                                          I->Obj->Obj.Setting,
                                          cSetting_ray_texture_settings));
            else
              ray->fWobble(ray,
                           SettingGet_i(G,I->Setting,
                                        NULL,cSetting_ray_texture),
                           SettingGet_3fv(G,I->Setting, NULL, 
                                          cSetting_ray_texture_settings));
            ray->fColor3fv(ray,ColorGet(G,I->Obj->Obj.Color));
          }
        
          if(r->fRender) { /* do OpenGL rendering in three passes */
            if(ray||pick) {
              
              /* here we need to iterate through and apply coordinate set matrices */

              r->fRender(r,info);
            } else {

              /* here we need to iterate through and apply coordinate set matrices */

              switch(a) {
              case cRepLabel:
                if(float_labels && (pass==-1))
                  r->fRender(r,info);
                else if(pass==1)
                  r->fRender(r,info);                  
                break;
              case cRepNonbondedSphere:
              case cRepRibbon:
              case cRepDot:
              case cRepCGO:
              case cRepCallback:
                if(pass==1) r->fRender(r,info);
                break;
              case cRepLine:
              case cRepMesh:
              case cRepDash:
              case cRepNonbonded:
              case cRepCell:
              case cRepExtent:
                if(!pass) r->fRender(r,info);                
                break;
              case cRepCyl: /* render sticks differently depending on transparency */
                if(SettingGet_f(G,r->cs->Setting,
                                r->obj->Setting,
                                cSetting_stick_transparency)>0.0001) {
                  if(pass==-1)
                    r->fRender(r,info);                                
                } else if(pass==1)
                  r->fRender(r,info);
                break;

              case cRepSurface:
                /*                if(pass==-1) r->fRender(r,ray,pick);              */
                if(SettingGet_f(G,r->cs->Setting,
                                r->obj->Setting,
                                cSetting_transparency)>0.0001) {
                  if(pass==-1)
                    r->fRender(r,info);                                
                } else if(pass==1)
                  r->fRender(r,info);
                break;
              case cRepSphere: /* render spheres differently depending on transparency */
                if(SettingGet_f(G,r->cs->Setting,
                                r->obj->Setting,
                                cSetting_sphere_transparency)>0.0001) {
                  if(pass==-1)
                    r->fRender(r,info);                                
                } else if(pass==1)
                  r->fRender(r,info);
                break;
              case cRepCartoon:
                if(SettingGet_f(G,r->cs->Setting,
                                r->obj->Setting,
                                cSetting_cartoon_transparency)>0.0001) {
                  if(pass==-1)
                    r->fRender(r,info);                                
                } else if(pass==1)
                  r->fRender(r,info);
                break;
              }

            }
          }
          /*          if(ray)
                      ray->fWobble(ray,0,NULL);*/
        }
  }
  PRINTFD(G,FB_CoordSet)
    " CoordSetRender: leaving...\n"
    ENDFD;

}
/*========================================================================*/
CoordSet *CoordSetNew(PyMOLGlobals *G)
{
  int a;
  OOAlloc(G,CoordSet);
  ObjectStateInit(G,&I->State);
  I->State.G=G;
  I->fFree=CoordSetFree;
  I->fRender=CoordSetRender;
  I->fUpdate=CoordSetUpdate;
  I->fEnumIndices=CoordSetEnumIndices;
  I->fExtendIndices=CoordSetExtendIndices;
  I->fAppendIndices=CoordSetAppendIndices;
  I->fInvalidateRep=CoordSetInvalidateRep;
  I->NIndex=0;
  I->NAtIndex=0;
  I->Coord = NULL;
  I->Color = NULL;
  I->AtmToIdx = NULL;
  I->IdxToAtm = NULL;
  I->NTmpBond = 0;
  I->TmpBond = NULL;
  I->TmpLinkBond = NULL;
  I->NTmpLinkBond = 0;
  I->PeriodicBox=NULL;
  I->PeriodicBoxType=cCSet_NoPeriodicity;
  /*  I->Rep=VLAlloc(Rep*,cRepCnt);*/
  I->NRep=cRepCnt;
  I->Symmetry = NULL;
  I->Name[0]=0;
  I->Obj = NULL;
  I->Spheroid = NULL;
  I->SpheroidNormal = NULL;
  I->SpheroidSphereSize = I->State.G->Sphere->Sphere[1]->nDot; /* does this make any sense? */
  for(a=0;a<I->NRep;a++)
	 I->Rep[a] = NULL;
  I->Setting = NULL;
  I->Coord2Idx = NULL;
  I->NMatrix = 0;
  I->MatrixVLA = NULL;
  I->SculptCGO = NULL;
  I->LabPos = NULL;
  return(I);
}
/*========================================================================*/
CoordSet *CoordSetCopy(CoordSet *cs)
{
  int a;
  int nAtom;
  float *v0,*v1;
  int *i0,*i1;
  OOAlloc(cs->State.G,CoordSet);

  (*I)=(*cs);
  ObjectStateCopy(&cs->State,&I->State);
  I->Symmetry=SymmetryCopy(cs->Symmetry);
  if(I->PeriodicBox) I->PeriodicBox=CrystalCopy(I->PeriodicBox);
  I->Coord = VLAlloc(float,I->NIndex*3);
  
  v0=I->Coord;
  v1=cs->Coord;
  for(a=0;a<I->NIndex;a++) {
    *(v0++)=*(v1++);
    *(v0++)=*(v1++);
    *(v0++)=*(v1++);
  }
  if(cs->LabPos) {
    I->LabPos = VLACalloc(LabPosType,I->NIndex);
    UtilCopyMem(I->LabPos,cs->LabPos,sizeof(LabPosType)*I->NIndex);
  }
  if(I->AtmToIdx) {
    nAtom = cs->Obj->NAtom;
    I->AtmToIdx = Alloc(int,nAtom);
    i0=I->AtmToIdx;
    i1=cs->AtmToIdx;
    for(a=0;a<nAtom;a++)
      *(i0++)=*(i1++);
  }
  
  
  if(cs->MatrixVLA) { /* not used yet */
    I->MatrixVLA = VLAlloc(double,16*cs->NMatrix*sizeof(double));
    if(I->MatrixVLA) {
      UtilCopyMem(I->MatrixVLA,cs->MatrixVLA,sizeof(double)*16*cs->NMatrix);
    }
  }

  I->IdxToAtm = Alloc(int,I->NIndex);
  i0=I->IdxToAtm;
  i1=cs->IdxToAtm;
  for(a=0;a<I->NIndex;a++)
    *(i0++)=*(i1++);
  
  /*  I->Rep=VLAlloc(Rep*,I->NRep); */
  i0=I->Active;
  i1=cs->Active;
  for(a=0;a<I->NRep;a++) {
    *(i0++)=*(i1++);
    I->Rep[a] = NULL;
  }


  I->TmpBond=NULL;
  I->Color=NULL;
  I->Spheroid=NULL;
  I->SpheroidNormal=NULL;
  I->Coord2Idx = NULL;
  return(I);
}
/*========================================================================*/
void CoordSetExtendIndices(CoordSet *I,int nAtom)
{
  int a,b;
  ObjectMolecule *obj = I->Obj;
  if(obj->DiscreteFlag) {
    if(obj->NDiscrete<nAtom) {
      VLACheck(obj->DiscreteAtmToIdx,int,nAtom);
      VLACheck(obj->DiscreteCSet,CoordSet*,nAtom);    
      for(a=obj->NDiscrete;a<nAtom;a++) {
        obj->DiscreteAtmToIdx[a]=-1;
        obj->DiscreteCSet[a]=NULL;
      }
      obj->NDiscrete = nAtom;
    }

    if(I->AtmToIdx) { /* convert to discrete if necessary */
      FreeP(I->AtmToIdx);
      for(a=0;a<I->NIndex;a++) { 
        b = I->IdxToAtm[a];
        obj->DiscreteAtmToIdx[b] = a;
        obj->DiscreteCSet[b] = I;
      }
    }
  }
  if(I->NAtIndex<nAtom)
	 {
		if(I->AtmToIdx) {
        I->AtmToIdx = Realloc(I->AtmToIdx,int,nAtom);
        if(nAtom){
          ErrChkPtr(I->State.G,I->AtmToIdx);
          for(a=I->NAtIndex;a<nAtom;a++)
            I->AtmToIdx[a]=-1;
        }
        I->NAtIndex = nAtom;
      } else if(!obj->DiscreteFlag) {
        I->AtmToIdx = Alloc(int,nAtom);
        for(a=0;a<nAtom;a++)
          I->AtmToIdx[a]=-1;
        I->NAtIndex = nAtom;
		}
	 }
}
/*========================================================================*/
void CoordSetAppendIndices(CoordSet *I,int offset) 
{
  int a,b;
  ObjectMolecule *obj = I->Obj;

  I->IdxToAtm = Alloc(int,I->NIndex);
  if(I->NIndex){
    ErrChkPtr(I->State.G,I->IdxToAtm);
    for(a=0;a<I->NIndex;a++)
      I->IdxToAtm[a]=a+offset;
  }
  if(obj->DiscreteFlag) {
    VLACheck(obj->DiscreteAtmToIdx,int,I->NIndex+offset);
    VLACheck(obj->DiscreteCSet,CoordSet*,I->NIndex+offset);
    for(a=0;a<I->NIndex;a++) {
      b=a+offset;
      obj->DiscreteAtmToIdx[b] = a;
      obj->DiscreteCSet[b] = I;
    }
  } else {
    I->AtmToIdx = Alloc(int,I->NIndex+offset);
    if(I->NIndex+offset){
      ErrChkPtr(I->State.G,I->AtmToIdx);
      for(a=0;a<offset;a++)
        I->AtmToIdx[a]=-1;
      for(a=0;a<I->NIndex;a++) 
        I->AtmToIdx[a+offset]=a;
    }
  }
  I->NAtIndex = I->NIndex + offset;
}
/*========================================================================*/
void CoordSetEnumIndices(CoordSet *I)
{
  /* set up for simple case where 1 = 1, etc. */
  int a;
  I->AtmToIdx = Alloc(int,I->NIndex);
  I->IdxToAtm = Alloc(int,I->NIndex);
  if(I->NIndex) {
    ErrChkPtr(I->State.G,I->AtmToIdx);
    ErrChkPtr(I->State.G,I->IdxToAtm);
    for(a=0;a<I->NIndex;a++)
      {
        I->AtmToIdx[a]=a;
        I->IdxToAtm[a]=a;
      }
  }
  I->NAtIndex = I->NIndex;
}
/*========================================================================*/
void CoordSetStrip(CoordSet *I)
{
  int a;
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a])
		I->Rep[a]->fFree(I->Rep[a]);
  I->NRep=0;
}
/*========================================================================*/
void CoordSetFree(CoordSet *I)
{
  int a;
  ObjectMolecule *obj;
  if(I)  {
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a]) 
		I->Rep[a]->fFree(I->Rep[a]);
    obj=I->Obj;
    if(obj)
      if(obj->DiscreteFlag) /* remove references to the atoms in discrete objects */
        for(a=0;a<I->NIndex;a++) {
          obj->DiscreteAtmToIdx[I->IdxToAtm[a]]=-1;
          obj->DiscreteCSet[I->IdxToAtm[a]]=NULL;
        } 
    FreeP(I->AtmToIdx);
    FreeP(I->IdxToAtm);
    VLAFreeP(I->Color);
    MapFree(I->Coord2Idx);
    VLAFreeP(I->Coord);
    /*    VLAFreeP(I->Rep);*/
    VLAFreeP(I->TmpBond);
    if(I->Symmetry) SymmetryFree(I->Symmetry);
    if(I->PeriodicBox) CrystalFree(I->PeriodicBox);
    FreeP(I->Spheroid);
    FreeP(I->SpheroidNormal);
    SettingFreeP(I->Setting);
    ObjectStatePurge(&I->State);
    CGOFree(I->SculptCGO);
    VLAFreeP(I->LabPos);
    OOFreeP(I);
  }
}


