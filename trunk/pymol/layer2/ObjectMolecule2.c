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
#include"os_gl.h"

#include"Base.h"
#include"Debug.h"
#include"Parse.h"
#include"OOMac.h"
#include"Vector.h"
#include"PConv.h"
#include"ObjectMolecule.h"
#include"Feedback.h"
#include"Util.h"
#include"AtomInfo.h"

float ObjectMoleculeGetMaxVDW(ObjectMolecule *I)
{
  float max_vdw = 0.0F;
  int a;
  AtomInfoType *ai;
  if(I->NAtom) {
    ai=I->AtomInfo;
    for(a=0;a<I->NAtom;a++) {
      if(max_vdw<ai->vdw)
        max_vdw = ai->vdw;
      ai++;
    }
  }
  return(max_vdw);
}
/*========================================================================*/
static PyObject *ObjectMoleculeCSetAsPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NCSet);
  for(a=0;a<I->NCSet;a++) {
    if(I->CSet[a]) {
      PyList_SetItem(result,a,CoordSetAsPyList(I->CSet[a]));
    } else {
      PyList_SetItem(result,a,Py_None);
      Py_INCREF(Py_None);
    }
  }
  return(PConvAutoNone(result));
}

/*static PyObject *ObjectMoleculeDiscreteCSetAsPyList(ObjectMolecule *I)
  {
  PyObject *result = NULL;
  return(PConvAutoNone(result));
  }*/


static int ObjectMoleculeCSetFromPyList(ObjectMolecule *I,PyObject *list)
{
  int ok=true;
  int a;
  if(ok) ok=PyList_Check(list);
  if(ok) {
    VLACheck(I->CSet,CoordSet*,I->NCSet);
    for(a=0;a<I->NCSet;a++) {
      
      if(ok) ok = CoordSetFromPyList(PyList_GetItem(list,a),&I->CSet[a]);
      if(ok) 
        if(I->CSet[a]) /* WLD 030205 */
          I->CSet[a]->Obj = I;
    }
  }
  return(ok);
}

static PyObject *ObjectMoleculeBondAsPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  PyObject *bond_list;
  BondType *bond;
  int a;

  result = PyList_New(I->NBond);  
  bond = I->Bond;
  for(a=0;a<I->NBond;a++) {
    bond_list=PyList_New(5);
    PyList_SetItem(bond_list,0,PyInt_FromLong(bond->index[0]));
    PyList_SetItem(bond_list,1,PyInt_FromLong(bond->index[1]));
    PyList_SetItem(bond_list,2,PyInt_FromLong(bond->order));
    PyList_SetItem(bond_list,3,PyInt_FromLong(bond->id));
    PyList_SetItem(bond_list,4,PyInt_FromLong(bond->stereo));
    PyList_SetItem(result,a,bond_list);
    bond++;
  }

  return(PConvAutoNone(result));
}

static int ObjectMoleculeBondFromPyList(ObjectMolecule *I,PyObject *list) 
{
  int ok=true;
  int a;
  PyObject *bond_list=NULL;
  BondType *bond;
  if(ok) ok=PyList_Check(list);  
  if(ok) ok=((I->Bond=VLAlloc(BondType,I->NBond))!=NULL);
  bond = I->Bond;
  for(a=0;a<I->NBond;a++) {
    if(ok) bond_list = PyList_GetItem(list,a);
    if(ok) ok = PyList_Check(bond_list);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,0),&bond->index[0]);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,1),&bond->index[1]);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,2),&bond->order);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,3),&bond->id);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,4),&bond->stereo);
    bond++;
  }
  return(ok);
}

static PyObject *ObjectMoleculeAtomAsPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  AtomInfoType *ai;
  int a;

  result = PyList_New(I->NAtom);  
  ai = I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    PyList_SetItem(result,a,AtomInfoAsPyList(ai));
    ai++;
  }
  return(PConvAutoNone(result));
}

static int ObjectMoleculeAtomFromPyList(ObjectMolecule *I,PyObject *list) 
{
  int ok=true;
  int a;
  AtomInfoType *ai;
  if(ok) ok=PyList_Check(list);  
  VLACheck(I->AtomInfo,AtomInfoType,I->NAtom+1);
  ai = I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(ok) ok = AtomInfoFromPyList(ai,PyList_GetItem(list,a));
    ai++;
  }
  return(ok);
}

int ObjectMoleculeNewFromPyList(PyObject *list,ObjectMolecule **result)
{
  int ok = true;
  ObjectMolecule *I=NULL;
  int discrete_flag;
  int ll;
  (*result) = NULL;
  

  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,8),&discrete_flag);

  I=ObjectMoleculeNew(discrete_flag);
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NCSet);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->NBond);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,3),&I->NAtom);
  if(ok) ok = ObjectMoleculeCSetFromPyList(I,PyList_GetItem(list,4));
  if(ok) ok = CoordSetFromPyList(PyList_GetItem(list,5),&I->CSTmpl);
  if(ok) ok = ObjectMoleculeBondFromPyList(I,PyList_GetItem(list,6));
  if(ok) ok = ObjectMoleculeAtomFromPyList(I,PyList_GetItem(list,7));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,8),&I->DiscreteFlag);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,9),&I->NDiscrete);
  if(ok) I->Symmetry = SymmetryNewFromPyList(PyList_GetItem(list,10));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,11),&I->CurCSet);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,12),&I->BondCounter);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,13),&I->AtomCounter);
  if(ok&&I->DiscreteFlag) {
    int *dcs = NULL;
    int a,i;
    CoordSet *cs;
    VLACheck(I->DiscreteAtmToIdx,int,I->NDiscrete);
    VLACheck(I->DiscreteCSet,CoordSet*,I->NDiscrete);
    if(ok) ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list,14),I->DiscreteAtmToIdx,I->NDiscrete);
    if(ok) ok = PConvPyListToIntArray(PyList_GetItem(list,15),&dcs);
    if(ok) {

      for(a=0;a<I->NDiscrete;a++) {
        i = dcs[a];
        I->DiscreteCSet[a] = NULL;
        if(i>=0) {
          if(i<I->NCSet) {
            cs = I->CSet[i];
            if(cs) {
              I->DiscreteCSet[a]=cs;
            }
          }
        }
      }
    }
    FreeP(dcs);
  }
  
  ObjectMoleculeInvalidate(I,cRepAll,cRepInvAll);
  if(ok) 
    (*result) = I;
  else {
    /* cleanup? */
  }

  return(ok);
}


/*========================================================================*/
PyObject *ObjectMoleculeAsPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;


  /* first, dump the atoms */

  result = PyList_New(16);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NCSet));
  PyList_SetItem(result,2,PyInt_FromLong(I->NBond));
  PyList_SetItem(result,3,PyInt_FromLong(I->NAtom));
  PyList_SetItem(result,4,ObjectMoleculeCSetAsPyList(I));
  PyList_SetItem(result,5,CoordSetAsPyList(I->CSTmpl));
  PyList_SetItem(result,6,ObjectMoleculeBondAsPyList(I));
  PyList_SetItem(result,7,ObjectMoleculeAtomAsPyList(I));
  PyList_SetItem(result,8,PyInt_FromLong(I->DiscreteFlag));
  PyList_SetItem(result,9,PyInt_FromLong(I->NDiscrete));
  PyList_SetItem(result,10,SymmetryAsPyList(I->Symmetry));
  PyList_SetItem(result,11,PyInt_FromLong(I->CurCSet));
  PyList_SetItem(result,12,PyInt_FromLong(I->BondCounter));
  PyList_SetItem(result,13,PyInt_FromLong(I->AtomCounter));


  if(I->DiscreteFlag) {
    int *dcs;
    int a;
    CoordSet *cs;
    
    /* get coordinate set indices */
    
    for(a=0;a<I->NCSet;a++) {
      cs = I->CSet[a];
      if(cs) {
        cs->tmp_index = a;
      }
    }
    
    dcs = Alloc(int,I->NDiscrete);
    
    for(a=0;a<I->NDiscrete;a++) {
      cs = I->DiscreteCSet[a];
      if(cs) 
        dcs[a] = cs->tmp_index;
      else
        dcs[a] = -1;
    }

    PyList_SetItem(result,14,PConvIntArrayToPyList(I->DiscreteAtmToIdx,I->NDiscrete));
    PyList_SetItem(result,15,PConvIntArrayToPyList(dcs,I->NDiscrete));
    FreeP(dcs);
  } else {
    PyList_SetItem(result,14,PConvAutoNone(NULL));
    PyList_SetItem(result,15,PConvAutoNone(NULL));
  }

#if 0
  CObject Obj;
  struct CoordSet **CSet;
  int NCSet;
  struct CoordSet *CSTmpl; /* template for trajectories, etc.*/
  BondType *Bond;
  AtomInfoType *AtomInfo;
  int NAtom;
  int NBond;
  int DiscreteFlag,NDiscrete;
  int *DiscreteAtmToIdx;
  struct CoordSet **DiscreteCSet;
  int CurCSet;
  int SeleBase; /* for internal usage by  selector & only valid during selection process */
  CSymmetry *Symmetry;
  int *Neighbor;
  float *UndoCoord[cUndoMask+1];
  int UndoState[cUndoMask+1];
  int UndoNIndex[cUndoMask+1];
  int UndoIter;
  CGO *UnitCellCGO;
  int BondCounter;
  int AtomCounter;
  struct CSculpt *Sculpt;
#endif

  return(PConvAutoNone(result));  
}


/*========================================================================*/
int ObjectMoleculeConnect(ObjectMolecule *I,BondType **bond,AtomInfoType *ai,
                          struct CoordSet *cs,int bondSearchFlag)
{
  #define cMULT 1

  int a,b,c,d,e,f,i,j;
  int a1,a2;
  float *v1,*v2,dst;
  int maxBond;
  MapType *map;
  int nBond;
  BondType *ii1,*ii2;
  int flag;
  int order;
  AtomInfoType *ai1,*ai2;
  float cutoff_s;
  float cutoff_h;
  float cutoff_v;
  float cutoff;
  float max_cutoff;
  int water_flag;

  cutoff_v=SettingGet(cSetting_connect_cutoff);
  cutoff_s=cutoff_v + 0.2F;
  cutoff_h=cutoff_v - 0.2F;
  max_cutoff = cutoff_s;

  /*  FeedbackMask[FB_ObjectMolecule]=0xFF;*/
  nBond = 0;
  maxBond = cs->NIndex * 8;
  (*bond) = VLAlloc(BondType,maxBond);
  if(cs->NIndex&&bondSearchFlag) /* &&(!I->DiscreteFlag) WLD 010527 */
	 {
      switch((int)SettingGet(cSetting_connect_mode)) {
      case 0:
        /* distance-based bond location  */

      map=MapNew(max_cutoff+MAX_VDW,cs->Coord,cs->NIndex,NULL);
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
                              dst = (float)diff3f(v1,v2);										
                              
                              a1=cs->IdxToAtm[i];
                              a2=cs->IdxToAtm[j];

                              ai1=ai+a1;
                              ai2=ai+a2;
                                    
                              dst -= ((ai1->vdw+ai2->vdw)/2);

                              /* quick hack for water detection.  
                                 they don't usually don't have CONECT records 
                                 and may not be HETATMs though they are supposed to be... */
                              
                              water_flag=false;
                              if((ai1->resn[0]=='W')&&
                                 (ai1->resn[1]=='A')&&
                                 (ai1->resn[2]=='T')&&
                                 (!ai1->resn[3]))
                                water_flag=true;
                              else if((ai1->resn[0]=='H')&&
                                      (ai1->resn[1]=='O')&&
                                      (ai1->resn[2]=='H')&&
                                      (!ai1->resn[3]))
                                water_flag=true;
                              if((ai2->resn[0]=='W')&&
                                 (ai2->resn[1]=='A')&&
                                 (ai2->resn[2]=='T')&&
                                 (!ai2->resn[3]))
                                water_flag=true;
                              else if((ai2->resn[0]=='H')&&
                                      (ai2->resn[1]=='O')&&
                                      (ai2->resn[2]=='H')&&
                                      (!ai2->resn[3]))
                                water_flag=true;
                              
                              cutoff = cutoff_h;
                              
                              /* workaround for hydrogens and sulfurs... */
                              
                              if(ai1->hydrogen||ai2->hydrogen)
                                cutoff = cutoff_h;
                              else if(((ai1->elem[0]=='S')&&(!ai1->elem[1]))||
                                   ((ai2->elem[0]=='S')&&(!ai2->elem[1])))
                                cutoff = cutoff_s;
                              else
                                cutoff = cutoff_v;
                              if( (dst <= cutoff)&&
                                  (!(ai1->hydrogen&&ai2->hydrogen))&&
                                  (water_flag||(!cs->TmpBond)||(!(ai1->hetatm&&ai2->hetatm))))
                                {
                                  flag=true;
                                  if(ai1->alt[0]!=ai2->alt[0]) { /* handle alternate conformers */
                                    if(ai1->alt[0]&&ai2->alt[0])
                                        flag=false; /* don't connect atoms with different, non-NULL
                                                       alternate conformations */
                                  } else if(ai1->alt[0]&&ai2->alt[0])
                                    if(!AtomInfoSameResidue(ai1,ai2))
                                      flag=false; /* don't connect non-NULL, alt conformations in 
                                                     different residues */
                                  if(ai1->alt[0]||ai2->alt[0]) 
                                  if(water_flag) /* hack to clean up water bonds */
                                    if(!AtomInfoSameResidue(ai1,ai2))
                                      flag=false;
                                      
                                  if(flag) {
                                    ai1->bonded=true;
                                    ai2->bonded=true;
                                    VLACheck((*bond),BondType,nBond);
                                    (*bond)[nBond].index[0] = a1;
                                    (*bond)[nBond].index[1] = a2;
                                    (*bond)[nBond].stereo = 0;
                                    order = 1;
                                    if((!ai1->hetatm)&&(!ai1->resn[3])) { /* Standard PDB residue */
                                      if(AtomInfoSameResidue(ai1,ai2)) {
                                        /* nasty high-speed hack to get bond valences and formal charges 
                                           for standard residues */
                                        if(((!ai1->name[1])&&(!ai2->name[1]))&&
                                           (((ai1->name[0]=='C')&&(ai2->name[0]=='O'))||
                                            ((ai1->name[0]=='O')&&(ai2->name[0]=='C')))) {
                                          order=2;
                                        } else {
                                          switch(ai1->resn[0]) {
                                          case 'A':
                                            switch(ai1->resn[1]) {
                                            case 'R': /* ARG */
                                              if(!strcmp(ai1->name,"NH1")) 
                                                ai1->formalCharge=1;
                                              else if(!strcmp(ai2->name,"NH1")) 
                                                ai2->formalCharge=1;
                                              if(((!strcmp(ai1->name,"CZ"))&&(!strcmp(ai2->name,"NH1")))||
                                                 ((!strcmp(ai2->name,"CZ"))&&(!strcmp(ai1->name,"NH1")))) 
                                                order=2;
                                              break;
                                            case 'S': 
                                              switch(ai1->resn[2]) {
                                              case 'P': /* ASP */
                                                if(!strcmp(ai1->name,"OD2")) 
                                                  ai1->formalCharge=-1;
                                                else if(!strcmp(ai2->name,"OD2")) 
                                                  ai2->formalCharge=-1;
                                              case 'N': /* ASN or ASP */
                                                if(((!strcmp(ai1->name,"CG"))&&(!strcmp(ai2->name,"OD1")))||
                                                   ((!strcmp(ai2->name,"CG"))&&(!strcmp(ai1->name,"OD1")))) 
                                                  order=2;
                                                break;
                                              }
                                            }
                                          case 'G':
                                            switch(ai1->resn[1]) {
                                            case 'L': 
                                              switch(ai1->resn[2]) {
                                              case 'U': /* GLU */
                                                if(!strcmp(ai1->name,"OE2")) 
                                                  ai1->formalCharge=-1;
                                                else if(!strcmp(ai2->name,"OE2")) 
                                                  ai2->formalCharge=-1;
                                              case 'N': /* GLN or GLU */
                                                if(((!strcmp(ai1->name,"CD"))&&(!strcmp(ai2->name,"OE1")))||
                                                   ((!strcmp(ai2->name,"CD"))&&(!strcmp(ai1->name,"OE1")))) 
                                                  order=2;
                                                break;
                                              }
                                            }
                                            break;
                                          case 'H':
                                            switch(ai1->resn[1]) {
                                            case 'I':
                                              switch(ai1->resn[2]) {
                                              case 'P':
                                                if(!strcmp(ai1->name,"ND1")) 
                                                  ai1->formalCharge=1;
                                                else if(!strcmp(ai2->name,"ND1")) 
                                                  ai2->formalCharge=1;
                                              case 'S':
                                              case 'E':
                                                if(((!strcmp(ai1->name,"CG"))&&(!strcmp(ai2->name,"CD2")))||
                                                   ((!strcmp(ai2->name,"CG"))&&(!strcmp(ai1->name,"CD2")))) 
                                                  order=2;
                                                else if(((!strcmp(ai1->name,"CE1"))&&(!strcmp(ai2->name,"ND1")))||
                                                        ((!strcmp(ai2->name,"CE1"))&&(!strcmp(ai1->name,"ND1")))) 
                                                  order=2;
                                                break;
                                                break;
                                              case 'D':
                                                if(((!strcmp(ai1->name,"CG"))&&(!strcmp(ai2->name,"CD2")))||
                                                   ((!strcmp(ai2->name,"CG"))&&(!strcmp(ai1->name,"CD2")))) 
                                                  order=2;
                                                else if(((!strcmp(ai1->name,"CE1"))&&(!strcmp(ai2->name,"NE2")))||
                                                        ((!strcmp(ai2->name,"CE1"))&&(!strcmp(ai1->name,"NE2")))) 
                                                  order=2;
                                                break;
                                              }
                                              break;
                                            }
                                            break;
                                          case 'P':
                                            switch(ai1->resn[1]) {
                                            case 'H': /* PHE */
                                              if(ai1->resn[2]=='E') {
                                                if(((!strcmp(ai1->name,"CG"))&&(!strcmp(ai2->name,"CD1")))||
                                                   ((!strcmp(ai2->name,"CG"))&&(!strcmp(ai1->name,"CD1")))) 
                                                  order=2;
                                                else if(((!strcmp(ai1->name,"CZ"))&&(!strcmp(ai2->name,"CE1")))||
                                                        ((!strcmp(ai2->name,"CZ"))&&(!strcmp(ai1->name,"CE1")))) 
                                                  order=2;
                                                
                                                else if(((!strcmp(ai1->name,"CE2"))&&(!strcmp(ai2->name,"CD2")))||
                                                        ((!strcmp(ai2->name,"CE2"))&&(!strcmp(ai1->name,"CD2")))) 
                                                  order=2;
                                                break; 
                                              }
                                            }
                                            break;
                                          case 'L':
                                            if(!strcmp(ai1->name,"NZ")) 
                                              ai1->formalCharge=1;
                                            else if(!strcmp(ai2->name,"NZ")) 
                                              ai2->formalCharge=1;
                                            break;
                                          case 'T':
                                            switch(ai1->resn[1]) {
                                            case 'Y': /* TYR */
                                              if(ai1->resn[2]=='R') {
                                                if(((!strcmp(ai1->name,"CG"))&&(!strcmp(ai2->name,"CD1")))||
                                                   ((!strcmp(ai2->name,"CG"))&&(!strcmp(ai1->name,"CD1")))) 
                                                  order=2;
                                                else if(((!strcmp(ai1->name,"CZ"))&&(!strcmp(ai2->name,"CE1")))||
                                                        ((!strcmp(ai2->name,"CZ"))&&(!strcmp(ai1->name,"CE1")))) 
                                                  order=2;
                                                
                                                else if(((!strcmp(ai1->name,"CE2"))&&(!strcmp(ai2->name,"CD2")))||
                                                        ((!strcmp(ai2->name,"CE2"))&&(!strcmp(ai1->name,"CD2")))) 
                                                  order=2;
                                                break; 
                                              }
                                              break;
                                            case 'R':
                                              if(ai1->resn[2]=='P') {
                                                if(((!strcmp(ai1->name,"CG"))&&(!strcmp(ai2->name,"CD1")))||
                                                   ((!strcmp(ai2->name,"CG"))&&(!strcmp(ai1->name,"CD1")))) 
                                                  order=2;
                                                else if(((!strcmp(ai1->name,"CZ3"))&&(!strcmp(ai2->name,"CE3")))||
                                                        ((!strcmp(ai2->name,"CZ3"))&&(!strcmp(ai1->name,"CE3")))) 
                                                  order=2;
                                                else if(((!strcmp(ai1->name,"CZ2"))&&(!strcmp(ai2->name,"CH2")))||
                                                        ((!strcmp(ai2->name,"CZ2"))&&(!strcmp(ai1->name,"CH2")))) 
                                                  order=2;
                                                else if(((!strcmp(ai1->name,"CE2"))&&(!strcmp(ai2->name,"CD2")))||
                                                        ((!strcmp(ai2->name,"CE2"))&&(!strcmp(ai1->name,"CD2")))) 
                                                  order=2;
                                                break; 
                                              }

                                              break;
                                            }
                                            
                                          }
                                        }
                                      }
                                    }
                                    (*bond)[nBond].order = order;
                                    nBond++;
                                  }
                                }
                            }
                          j=MapNext(map,j);
                        }
                    }
            }
          MapFree(map);
        case 1: /* only use explicit connectivity from file (don't do anything) */ 
          break;
        case 2:  /* dictionary-based connectivity */
          /* TODO */
          
          break;
        }
      }
      PRINTFB(FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeConnect: Found %d bonds.\n",nBond
        ENDFB;
      if(Feedback(FB_ObjectMolecule,FB_Debugging)) {
        for(a=0;a<nBond;a++)
          printf(" ObjectMoleculeConnect: bond %d ind0 %d ind1 %d\n",
                 a,(*bond)[a].index[0],(*bond)[a].index[1]);
      }
    }

  if(cs->NTmpBond&&cs->TmpBond) {
      PRINTFB(FB_ObjectMolecule,FB_Blather) 
      " ObjectMoleculeConnect: incorporating explicit bonds. %d %d\n",
             nBond,cs->NTmpBond
        ENDFB;
    VLACheck((*bond),BondType,(nBond+cs->NTmpBond));
    ii1=(*bond)+nBond;
    ii2=cs->TmpBond;
    for(a=0;a<cs->NTmpBond;a++)
      {
        a1 = cs->IdxToAtm[ii2->index[0]]; /* convert bonds from index space */
        a2 = cs->IdxToAtm[ii2->index[1]]; /* to atom space */
        ai[a1].bonded=true;
        ai[a2].bonded=true;
        ii1->index[0]=a1;
        ii1->index[1]=a2;
        ii1->order = ii2->order;
        ii1->stereo = ii2->stereo;
        ii1++;
        ii2++;

      }
    nBond=nBond+cs->NTmpBond;
    VLAFreeP(cs->TmpBond);
    cs->NTmpBond=0;
  }

  if(cs->NTmpLinkBond&&cs->TmpLinkBond) {
    PRINTFB(FB_ObjectMolecule,FB_Blather) 
      "ObjectMoleculeConnect: incorporating linkage bonds. %d %d\n",
      nBond,cs->NTmpLinkBond
      ENDFB;
    VLACheck((*bond),BondType,(nBond+cs->NTmpLinkBond));
    ii1=(*bond)+nBond;
    ii2=cs->TmpLinkBond;
    for(a=0;a<cs->NTmpLinkBond;a++)
      {
        a1 = ii2->index[0]; /* first atom is in object */
        a2 = cs->IdxToAtm[ii2->index[1]]; /* second is in the cset */
        ai[a1].bonded=true;
        ai[a2].bonded=true;
        ii1->index[0]=a1;
        ii1->index[1]=a2;
        ii1->order = ii2->order;
        ii1->stereo = ii2->stereo;
        ii1++;
        ii2++;
      }
    nBond=nBond+cs->NTmpLinkBond;
    VLAFreeP(cs->TmpLinkBond);
    cs->NTmpLinkBond=0;
  }

  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeConnect: elminating duplicates with %d bonds...\n",nBond
    ENDFD;

  if(!I->DiscreteFlag) {
    UtilSortInPlace((*bond),nBond,sizeof(BondType),(UtilOrderFn*)BondInOrder);
    if(nBond) { /* eliminate duplicates */
      ii1=(*bond)+1;
      ii2=(*bond)+1;
      a=nBond-1;
      nBond=1;
      if(a>0) 
        while(a--) { 
          if((ii2->index[0]!=(ii1-1)->index[0])||
             (ii2->index[1]!=(ii1-1)->index[1])) {
            *(ii1++)=*(ii2++); /* copy bond */
            nBond++;
          } else {
            ii2++; /* skip bond */
          }
        }
      VLASize((*bond),BondType,nBond);
    }
  }
  PRINTFD(FB_ObjectMolecule)
    " ObjectMoleculeConnect: leaving with %d bonds...\n",nBond
    ENDFD;
  return(nBond);
}


/*========================================================================*/
void ObjectMoleculeSort(ObjectMolecule *I) /* sorts atoms and bonds */
{
  int *index,*outdex;
  int a,b;
  CoordSet *cs,**dcs;
  AtomInfoType *atInfo;
  int *dAtmToIdx;

  if(!I->DiscreteFlag) {

    index=AtomInfoGetSortedIndex(I->AtomInfo,I->NAtom,&outdex);
    for(a=0;a<I->NBond;a++) { /* bonds */
      I->Bond[a].index[0]=outdex[I->Bond[a].index[0]];
      I->Bond[a].index[1]=outdex[I->Bond[a].index[1]];
    }
    
    for(a=-1;a<I->NCSet;a++) { /* coordinate set mapping */
      if(a<0) {
        cs=I->CSTmpl;
      } else {
        cs=I->CSet[a];
      }
      
      if(cs) {
        for(b=0;b<cs->NIndex;b++)
          cs->IdxToAtm[b]=outdex[cs->IdxToAtm[b]];
        if(cs->AtmToIdx) {
          for(b=0;b<I->NAtom;b++)
            cs->AtmToIdx[b]=-1;
          for(b=0;b<cs->NIndex;b++) {
            cs->AtmToIdx[cs->IdxToAtm[b]]=b;
          }
        }
      }
    }
    
    atInfo=(AtomInfoType*)VLAMalloc(I->NAtom,sizeof(AtomInfoType),5,true);
    /* autozero here is important */
    for(a=0;a<I->NAtom;a++)
      atInfo[a]=I->AtomInfo[index[a]];
    VLAFreeP(I->AtomInfo);
    I->AtomInfo=atInfo;
    
    if(I->DiscreteFlag) {
      dcs = VLAlloc(CoordSet*,I->NAtom);
      dAtmToIdx = VLAlloc(int,I->NAtom);
      for(a=0;a<I->NAtom;a++) {
        b=index[a];
        dcs[a] = I->DiscreteCSet[b];
        dAtmToIdx[a] = I->DiscreteAtmToIdx[b];
      }
      VLAFreeP(I->DiscreteCSet);
      VLAFreeP(I->DiscreteAtmToIdx);
      I->DiscreteCSet = dcs;
      I->DiscreteAtmToIdx = dAtmToIdx;
    }
    AtomInfoFreeSortedIndexes(index,outdex);

    UtilSortInPlace(I->Bond,I->NBond,sizeof(BondType),(UtilOrderFn*)BondInOrder);
    /* sort...important! */
    ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms); /* important */

  }
}
