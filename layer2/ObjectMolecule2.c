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
#include"Selector.h"
#include"ObjectDist.h"
#include"Executive.h"
#include"P.h"
#include"ObjectCGO.h"
#include"Scene.h"

#define ntrim ParseNTrim
#define nextline ParseNextLine
#define ncopy ParseNCopy
#define nskip ParseNSkip

#define cResvMask 0x7FFF

#define cMaxOther 6

typedef struct {
  int n_cyclic_arom,  cyclic_arom[cMaxOther];
  int n_arom, arom[cMaxOther];
  int n_high_val, high_val[cMaxOther];
  int n_cyclic, cyclic[cMaxOther];
  int n_planer, planer[cMaxOther];
  int n_rest, rest[cMaxOther];
  int score;
} OtherRec;

static int populate_other(OtherRec *other,int at,AtomInfoType *ai,BondType *bd, int *neighbor)
{
  int five_cycle = false;
  int six_cycle = false;
  
  {        
    int mem[9], nbr[7];
    const int ESCAPE_MAX = 500;
    register int escape_count; 
    
    escape_count = ESCAPE_MAX; /* don't get bogged down with structures 
                                  that have unreasonable connectivity */
    mem[0] = bd->index[0];
    mem[1] = bd->index[1];
    nbr[1] = neighbor[mem[1]]+1;
    while(((mem[2] = neighbor[nbr[1]])>=0)) {
      if(mem[2]!=mem[0]) {
        nbr[2] = neighbor[mem[2]]+1;
        while(((mem[3] = neighbor[nbr[2]])>=0)) {
          if(mem[3]!=mem[1]) {
            nbr[3] = neighbor[mem[3]]+1;
            while(((mem[4] = neighbor[nbr[3]])>=0)) {
              if((mem[4]!=mem[2])&&(mem[4]!=mem[1])&&(mem[4]!=mem[0])) {            
                nbr[4] = neighbor[mem[4]]+1;              
                while(((mem[5] = neighbor[nbr[4]])>=0)) {
                  if(!(escape_count--)) goto escape;
                  if((mem[5]!=mem[3])&&(mem[5]!=mem[2])&&(mem[5]!=mem[1])) { 
                    if(mem[5]==mem[0]) { /* five-cycle */
                      five_cycle = true;
                    }
                    nbr[5] = neighbor[mem[5]]+1;              
                    while(((mem[6] = neighbor[nbr[5]])>=0)) {
                      if((mem[6]!=mem[4])&&(mem[6]!=mem[3])&&(mem[6]!=mem[2])&&(mem[6]!=mem[1])) {
                        if(mem[6]==mem[0]) {  /* six-cycle */
                          six_cycle = true;
                        }
                      }
                      nbr[5]+=2;                          
                    }
                  }
                  nbr[4]+=2;
                }
              }
              nbr[3]+=2;
            }
          }
          nbr[2]+=2;
        }
      }
      nbr[1]+=2;
    }
  }
 escape:

  if(bd->order==4) { /* aromatic */
    if((five_cycle || six_cycle) && (other->n_cyclic_arom<cMaxOther)) {
      other->cyclic_arom[other->n_cyclic_arom++]=at;
      if(five_cycle && six_cycle) 
        other->score+=34;
      else if(five_cycle)
        other->score+=33;
      else
        other->score+=32;
      return 1;
    } else if(other->n_arom<cMaxOther) {
      other->arom[other->n_arom++]=at;
      other->score+=64;
      return 1;
    }
  }
  if(bd->order>1) {
    if(other->n_high_val<cMaxOther) {
      other->high_val[other->n_high_val++]=at;
      other->score+=16;
      return 1;
    }
  }
  if(five_cycle || six_cycle) {
    if(other->n_cyclic<cMaxOther) {
      other->cyclic[other->n_cyclic++]=at;
      other->score+=8;
      return 1;
    }
  }
  if(ai->geom==cAtomInfoPlaner) {
    if(other->n_planer<cMaxOther) {
      other->planer[other->n_planer++]=at;
      other->score+=4;
      return 1;
    }
  }
  if(other->n_rest<cMaxOther) {
    other->rest[other->n_rest++]=at;
    other->score+=1;
    return 1;
  }
  return 0;
}

static int append_index(int *result, int offset, int a1, int a2, int score, int ar_count)
{
  int c;
  c=result[a1];
  while(c<offset) {
    if(result[c]==a2) { /* already entered */
      if(result[c+1]<score) {
        result[c+1]=score;
        result[c+2]=ar_count;
      }
      return offset;
    }
    c+=3;
  }
  result[offset++]=a2;
  result[offset++]=score;
  result[offset++]=ar_count;
  return offset;
}

int ObjectMoleculeAddPseudoatom(ObjectMolecule *I,int sele_index, char *name, 
                                char *resn, char *resi, char  *chain,
                                char *segi, char *elem, float vdw, 
                                int hetatm, float b, float q, char *label, 
                                float *pos, 
                                int color, int state, int mode, int quiet)
{
  PyMOLGlobals *G = I->Obj.G;
  int start_state=0, stop_state = 0;
  int nAtom = 1;
  int extant_only = false;
  int ai_merged = false;
  float pos_array[3] = { 0.0F, 0.0F, 0.0F };
  int ok=true;

  AtomInfoType *atInfo = VLACalloc(AtomInfoType,1);

  if(state>=0) { /* specific state */
    start_state = state;
    stop_state = state+1;
  } else if(state==-1) { /* current state */
    start_state = ObjectGetCurrentState(&I->Obj,true);
    stop_state = start_state+1;
  } else { /* all states */
    if(sele_index>=0) {
      start_state = 0;
      stop_state = SelectorCountStates(G,sele_index);
      if(state==-3)
        extant_only = true;
    } else {
      start_state = 0;
      stop_state = ExecutiveCountStates(G,cKeywordAll);
      if(stop_state<1)
        stop_state = 1;
    }
  }
  {
    /* match existing properties of the old atom */
    int auto_show_lines = (int)SettingGet(G,cSetting_auto_show_lines);
    int auto_show_spheres = (int)SettingGet(G,cSetting_auto_show_spheres);
    int auto_show_nonbonded = (int)SettingGet(G,cSetting_auto_show_nonbonded);
    AtomInfoType *ai = atInfo;
    ai->resv = AtomResvFromResi(resi);
    ai->hetatm=hetatm;
    ai->geom=cAtomInfoNone;
    ai->q=q;
    ai->b=b;
    strcpy(ai->chain,chain);
    strcpy(ai->resi,resi);
    strcpy(ai->segi,segi);
    strcpy(ai->resn,resn);  
    strcpy(ai->name,name);  
    strcpy(ai->elem,elem);  
    ai->visRep[cRepLine] = auto_show_lines; 
    ai->visRep[cRepNonbonded] = auto_show_nonbonded; 
    ai->visRep[cRepSphere] = auto_show_spheres; 
    ai->id=-1;
    ai->rank=-1;
    if(vdw>=0.0F) 
      ai->vdw = vdw;
    else
      ai->vdw = 1.0F;
    if(label[0]) {
      OVreturn_word ret = OVLexicon_GetFromCString(
                                                   G->Lexicon,label);
      if(OVreturn_IS_OK(ret)) {
        ai->label = ret.word;
        ai->visRep[cRepLabel] = true;
        ai->visRep[cRepLine] = false;
        ai->visRep[cRepNonbonded] = false;
        ai->visRep[cRepSphere] = false;
      }
    }
    if(color<0) {
      AtomInfoAssignColors(I->Obj.G,ai); 
      if((ai->elem[0]=='C')&&(ai->elem[1]==0)) 
        /* carbons are always colored according to the object color */
        ai->color=I->Obj.Color;
    } else {
      ai->color = color;
    }
    AtomInfoAssignParameters(I->Obj.G,ai);
    AtomInfoUniquefyNames(I->Obj.G,I->AtomInfo,I->NAtom,ai,1);
    if(!quiet) {
      PRINTFB(G,FB_ObjectMolecule,FB_Actions)
        " ObjMol: created %s/%s/%s/%s`%s/%s\n",
        I->Obj.Name,ai->segi,ai->chain,ai->resn,ai->resi,ai->name
        ENDFB(G);
    }
  }

  for(state=start_state;state<stop_state;state++) {

    CoordSet *cset = NULL;

    if((extant_only&&(state<I->NCSet)&&I->CSet[state]) ||
       !extant_only) {
      
      if(sele_index>=0) {
        ObjectMoleculeOpRec op;
        ObjectMoleculeOpRecInit(&op);
        op.code = OMOP_CSetSumVertices;
        op.cs1 = state;

        ExecutiveObjMolSeleOp(I->Obj.G,sele_index,&op);
        
        if(op.i1) {
          float factor = 1.0F/op.i1;
          scale3f(op.v1, factor, pos_array);
          pos = pos_array;

          if(vdw<0.0F) {
            switch(mode) {
            case 1:
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_CSetMaxDistToPt;
              copy3f(pos_array, op.v1);
              op.cs1 = state;
              ExecutiveObjMolSeleOp(I->Obj.G,sele_index,&op);
              vdw = op.f1;
              break;
            case 2:
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_CSetSumSqDistToPt;
              copy3f(pos_array, op.v1);
              op.cs1 = state;
              ExecutiveObjMolSeleOp(I->Obj.G,sele_index,&op);
              vdw = sqrt1f(op.d1/op.i1);
              break;
            case 0:
            default: 
              vdw = 0.5F;
              break;
            }
            if(vdw>=0.0F) 
              atInfo->vdw = vdw; /* NOTE: only uses vdw from first state selection...*/
          }
        } else {
          pos = NULL; /* skip this state */
        }
      } else if(!pos) {
        pos = pos_array;
        SceneGetPos(I->Obj.G,pos);
      }
      
      if(pos) { /* only add coordinate to state if we have position for it */

        float *coord=VLAlloc(float,3*nAtom);
        
        copy3f(pos,coord);
        
        cset = CoordSetNew(I->Obj.G);
        cset->NIndex=nAtom;
        cset->Coord=coord;
        cset->TmpBond=NULL;
        cset->NTmpBond=0;
        
        cset->Obj = I;
        if(cset->fEnumIndices) cset->fEnumIndices(cset);
        if(!ai_merged) {
          ObjectMoleculeMerge(I,atInfo,cset,false,cAIC_AllMask,true); /* NOTE: will release atInfo */
          ObjectMoleculeExtendIndices(I,-1);
          ObjectMoleculeUpdateNeighbors(I);
          ai_merged = true;
        }
        if(state>=I->NCSet) {
          VLACheck(I->CSet,CoordSet*,state);        
          I->NCSet=state+1;
        }
        if(!I->CSet[state]) {
          /* new coordinate set */
          I->CSet[state] = cset;
          cset = NULL;
        } else {
          /* merge coordinate set */
          CoordSetMerge(I->CSet[state],cset); 
          if(cset->fFree) {
            cset->fFree(cset);
            cset = NULL;
          }
        }
      }
    }
  }
  if(ai_merged) {
    ObjectMoleculeSort(I);
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
    ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms,-1);
  } else {
    VLAFreeP(atInfo);
  }
  return(ok);
}

int *ObjectMoleculeGetPrioritizedOtherIndexList(ObjectMolecule *I,CoordSet *cs)
{
  int a,b;
  int b1,b2,a1,a2,a3;
  OtherRec *o;
  OtherRec *other=Calloc(OtherRec,cs->NIndex);
  int *result = NULL;
  int offset;
  int n_alloc=0;
  BondType *bd;

  ObjectMoleculeUpdateNeighbors(I);
  bd=I->Bond;
  for(a=0;a<I->NBond;a++) {
    b1 = bd->index[0];
    b2 = bd->index[1];
    if(I->DiscreteFlag) {
      if((cs==I->DiscreteCSet[b1])&&(cs==I->DiscreteCSet[b2])) {
        a1=I->DiscreteAtmToIdx[b1];
        a2=I->DiscreteAtmToIdx[b2];
      } else {
        a1=-1;
        a2=-1;
      }
    } else {
      a1=cs->AtmToIdx[b1];
      a2=cs->AtmToIdx[b2];
    }
    if((a1>=0)&&(a2>=0)) {
      n_alloc+=populate_other(other+a1,a2,I->AtomInfo+b2,bd,I->Neighbor);
      n_alloc+=populate_other(other+a2,a1,I->AtomInfo+b1,bd,I->Neighbor);
    }
    bd++;
  }

  n_alloc = 3*(n_alloc+cs->NIndex);
  o=other;
  result = Alloc(int,n_alloc);
  for(a=0;a<cs->NIndex;a++) {
    result[a]=-1;
  }
  offset = cs->NIndex;
  bd=I->Bond;
  for(a=0;a<I->NBond;a++) {
    b1 = bd->index[0];
    b2 = bd->index[1];
    if(I->DiscreteFlag) {
      if((cs==I->DiscreteCSet[b1])&&(cs==I->DiscreteCSet[b2])) {
        a1=I->DiscreteAtmToIdx[b1];
        a2=I->DiscreteAtmToIdx[b2];
      } else {
        a1=-1;
        a2=-1;
      }
    } else {
      a1=cs->AtmToIdx[b1];
      a2=cs->AtmToIdx[b2];
    }
    if((a1>=0)&&(a2>=0))
      {
        if(result[a1]<0) {
          o = other+a1;
          result[a1]=offset;
          for(b=0;b<o->n_cyclic_arom;b++) {
            a3 = o->cyclic_arom[b];
            offset = append_index(result, offset, a1, a3, 128 + other[a3].score, 1);
          }
          for(b=0;b<o->n_arom;b++) {
            a3 = o->arom[b];
            offset = append_index(result, offset, a1, a3, 64 + other[a3].score, 1);
          }
          for(b=0;b<o->n_high_val;b++) {
            a3 = o->high_val[b];
            offset = append_index(result, offset, a1, a3, 16 + other[a3].score, 0);
          }
          for(b=0;b<o->n_cyclic;b++) {
            a3 = o->cyclic[b];
            offset = append_index(result, offset, a1, a3, 8 + other[a3].score, 0);
          }
          for(b=0;b<o->n_planer;b++) {
            a3 = o->planer[b];
            offset = append_index(result, offset, a1, a3, 2 + other[a3].score, 0);
          }
          for(b=0;b<o->n_rest;b++) {
            a3 = o->rest[b];
            offset = append_index(result, offset, a1, a3, 1 + other[a3].score, 0);
          }
          result[offset++]=-1;
        }

        if(result[a2]<0) {
          o = other+a2;
          result[a2]=offset;
          for(b=0;b<o->n_cyclic_arom;b++) {
            a3 = o->cyclic_arom[b];
            offset = append_index(result, offset, a2, a3, 128 + other[a3].score, 1);
          }
          for(b=0;b<o->n_arom;b++) {
            a3 = o->arom[b];
            offset = append_index(result, offset, a2, a3, 64 + other[a3].score, 1);
          }
          for(b=0;b<o->n_high_val;b++) {
            a3 = o->high_val[b];
            offset = append_index(result, offset, a2, a3, 16 + other[a3].score, 0);
          }
          for(b=0;b<o->n_cyclic;b++) {
            a3 = o->cyclic[b];
            offset = append_index(result, offset, a2, a3, 8 + other[a3].score, 0);
          }
          for(b=0;b<o->n_planer;b++) {
            a3 = o->planer[b];
            offset = append_index(result, offset, a2, a3, 2 + other[a3].score, 0);
          }
          for(b=0;b<o->n_rest;b++) {
            a3 = o->rest[b];
            offset = append_index(result, offset, a2, a3, 1 + other[a3].score, 0);
          }
          result[offset++]=-1;
        }

      }
    bd++;
  }
  FreeP(other);
  return result;
}

int ObjectMoleculeGetNearestBlendedColor(ObjectMolecule *I, float *point, 
                                         float cutoff, int state, float *dist,
                                         float *color, int sub_vdw)
{
  int result = -1;
  float tot_weight = 0.0F;
  float cutoff2 = cutoff * cutoff;
  register float nearest = -1.0F;


  color[0] = 0.0F;
  color[1] = 0.0F;
  color[2] = 0.0F;

  if(state<0)
    state = ObjectGetCurrentState(&I->Obj,true);
  
  if((state>=0)&&(state<I->NCSet)) {
    CoordSet *cs = I->CSet[state];
    if(cs) {
      MapType *map;
      CoordSetUpdateCoord2IdxMap(cs, cutoff);
      if(sub_vdw) {
        cutoff -= MAX_VDW;
        cutoff2 = cutoff * cutoff;
      }
      nearest = cutoff2;
      if( (map = cs->Coord2Idx)) {
        int a,b,c,d,e,f,j;
        register float test;
        register float *v;
        MapLocus(map,point,&a,&b,&c);
        for(d=a-1;d<=a+1;d++)
          for(e=b-1;e<=b+1;e++)
            for(f=c-1;f<=c+1;f++) {
              j = *(MapFirst(map,d,e,f));
              while(j>=0) {
                v = cs->Coord + (3*j);                    
                test = diffsq3f(v,point);
                if(sub_vdw) {
                  test = sqrt1f(test);
                  test -= I->AtomInfo[cs->IdxToAtm[j]].vdw;
                  if (test<0.0F)
                    test = 0.0F;
                  test = test*test;
                }
                if(test<cutoff2) {
                  float weight = cutoff - sqrt1f(test); 
                  float *at_col = ColorGet(I->Obj.G,I->AtomInfo[cs->IdxToAtm[j]].color);
                  color[0] += at_col[0] * weight;
                  color[1] += at_col[1] * weight;
                  color[2] += at_col[2] * weight;
                  tot_weight += weight;
                }
                if(test<=nearest) {
                  result = j;
                  nearest = test;
                }
                j=MapNext(map,j);
              }
            }
      } else {
        register int j;
        register float test,*v=cs->Coord;
        for(j=0;j<cs->NIndex;j++) {
          test = diffsq3f(v,point);
          if(sub_vdw) {
            test = sqrt1f(test);
            test -= I->AtomInfo[cs->IdxToAtm[j]].vdw;
            if (test<0.0F)
              test = 0.0F;
            test = test*test;
          }
          if(test<cutoff2) {
            float weight = cutoff - sqrt1f(test); 
            float *color = ColorGet(I->Obj.G,I->AtomInfo[cs->IdxToAtm[j]].color);
            color[0] += color[0] * weight;
            color[1] += color[1] * weight;
            color[2] += color[2] * weight;
            tot_weight += weight;
          }
          if(test<=nearest) {
            result = j;
            nearest = test;
          }
          v+=3;
        }
      }
      if(result>=0)
        result = cs->IdxToAtm[result];
    }
  }
  if(dist) {
    if(result>=0) {
      *dist = sqrt1f(nearest);
      if(tot_weight>0.0F) {
        color[0] /= tot_weight;
        color[1] /= tot_weight;
        color[2] /= tot_weight;
      }
    } else {
      *dist = -1.0F;
    }
  }
  return result;
}

int ObjectMoleculeGetNearestAtomIndex(ObjectMolecule *I, float *point, float cutoff, int state, float *dist)
{
  int result = -1;
  register float nearest = -1.0F;
  if(state<0)
    state = ObjectGetCurrentState(&I->Obj,true);
  if((state>=0)&&(state<I->NCSet)) {
    CoordSet *cs = I->CSet[state];
    if(cs) {
      MapType *map;
      CoordSetUpdateCoord2IdxMap(cs, cutoff);
      nearest = cutoff * cutoff;
      if( (map = cs->Coord2Idx)) {
        int a,b,c,d,e,f,j;
        register float test;
        register float *v;
        MapLocus(map,point,&a,&b,&c);
        for(d=a-1;d<=a+1;d++)
          for(e=b-1;e<=b+1;e++)
            for(f=c-1;f<=c+1;f++) {
              j = *(MapFirst(map,d,e,f));
              while(j>=0) {
                v = cs->Coord + (3*j);                    
                test = diffsq3f(v,point);
                if(test<=nearest) {
                  result = j;
                  nearest = test;
                }
                j=MapNext(map,j);
              }
            }
      } else {
        register int j;
        register float test,*v=cs->Coord;
        for(j=0;j<cs->NIndex;j++) {
          test = diffsq3f(v,point);
          if(test<=nearest) {
            result = j;
            nearest = test;
          }
          v+=3;
        }
      }
      if(result>=0)
        result = cs->IdxToAtm[result];
    }
  }
  if(dist) {
    if(result>=0) {
      *dist = sqrt1f(nearest);
    } else {
      *dist = -1.0F;
    }
  }
  return result;
}

int ObjectMoleculeGetPrioritizedOther(int *other, int a1, int a2, int *double_sided)
     
{
  int a3 = -1;
  int lvl=-1,ck,ck_lvl;
  int offset;
  int ar_count = 0;

  a3 = -1;
  lvl = -1;
  if(a1>=0) {
    offset = other[a1];
    if(offset>=0) {
      while(1) {
        ck = other[offset];
        if(ck!=a2) {
          if(ck>=0) {
            ck_lvl = other[offset+1];
            if(ck_lvl>lvl) {
              a3 = ck;
              lvl = ck_lvl;
            }
            ar_count+=other[offset+2];
          } else
            break;
        }
        offset+=3;
      }
    }
  }
  if(a2>=0) {
    offset = other[a2];
    if(offset>=0) {
      while(1) {
        ck = other[offset];
        if(ck!=a1) {
          if(ck>=0) {
            ck_lvl = other[offset+1];
            if(ck_lvl>lvl) {
              a3 = ck;
              lvl = ck_lvl;
            }
            ar_count+=other[offset+2];
          } else
            break;
        }
        offset+=3;
      }
    }
  }

  if(double_sided) {
    if(ar_count==4)
      *double_sided=true;
    else
      *double_sided=false;
  }
  return a3;
}


int ObjectMoleculeIsAtomBondedToName(ObjectMolecule *obj,int a0,char *name)
{
  int a2,s;
  int bonded =false;
  
  if(a0>=0) { 
    s=obj->Neighbor[a0]; 
    s++; /* skip count */
    while(1) {
      a2 = obj->Neighbor[s];
      if(a2<0)
        break;
      if(WordMatch(obj->Obj.G,obj->AtomInfo[a2].name,name,true)<0)
        bonded = true;
      break;
      s+=2;
    }
  }
  return bonded;
}

int ObjectMoleculeAreAtomsBonded2(ObjectMolecule *obj0,int a0, ObjectMolecule *obj1,int a1)
{
  /* assumes neighbor list is current */

  if(obj0!=obj1)
    return false;
  else {
    int a2,s;
    
    if(a0>=0) { 
      s=obj0->Neighbor[a0]; 
      s++; /* skip count */
      while(1) {
        a2 = obj0->Neighbor[s];
        if(a2<0)
          break;
        if(a1==a2)
          return true;
        s+=2;
      }
    }
  }
  return false;
}

int ObjectMoleculeDoesAtomNeighborSele(ObjectMolecule *I, int index, int sele)
{
  int result = false;
  ObjectMoleculeUpdateNeighbors(I);
  if(index<I->NAtom) {
    int a1;
    int n;
    AtomInfoType *ai;

    n = I->Neighbor[index]+1;
    while(1) { /* look for an attached non-hydrogen as a base */
      a1 = I->Neighbor[n];
      n+=2; 
      if(a1<0) break;
      ai = I->AtomInfo + a1;
      if(SelectorIsMember(I->Obj.G,ai->selEntry,sele)) {
        result=true;
        break;
      }
    }
  }
  return result;
}

static void assign_pdb_known_residue(PyMOLGlobals *G, AtomInfoType *ai1, AtomInfoType *ai2, int *bond_order)
{
  int order = *(bond_order);
  char *name1 = ai1->name;
  char *name2 = ai2->name;
  char *resn1 = ai1->resn;
    
  /* nasty high-speed hack to get bond valences and formal charges 
     for standard residues */
  if(((!name1[1])&&(!name2[1]))&&
     (((name1[0]=='C')&&(name2[0]=='O'))||
      ((name1[0]=='O')&&(name2[0]=='C')))) {
    order=2;
  } else if((!name2[1])&&(name2[0]=='C')&&(!strcmp(name1,"OXT"))) {
    ai1->formalCharge = -1;
    ai1->chemFlag = false;
  } else if((!name1[1])&&(name1[0]=='C')&&(!strcmp(name2,"OXT"))) {
    ai2->formalCharge = -1;
    ai2->chemFlag = false;
  } else {
    switch(resn1[0]) {
    case 'A':
      switch(resn1[1]) {
      case 'R': 
        switch(resn1[2]) {   
        case 'G': /* ARG... */
          switch(resn1[3]) {   
          case 0:
          case 'P': /*  ARG, ARGP */
            if(!strcmp(name1,"NH1"))  {
              ai1->formalCharge=1;
              ai1->chemFlag = false;
            } else if(!strcmp(name2,"NH1"))  {
              ai2->formalCharge=1;
              ai2->chemFlag = false;
            }
            break;
          }
          if(((!strcmp(name1,"CZ"))&&(!strcmp(name2,"NH1")))||
             ((!strcmp(name2,"CZ"))&&(!strcmp(name1,"NH1")))) 
            order=2;
          break;
        }
        break;
      case 'S': 
        switch(resn1[2]) {
        case 'P': /* ASP... */
          switch(resn1[3]) {
          case 0:
          case 'M': /* ASP, ASPM minus assumption */
            if(!strcmp(name1,"OD2")) {
              ai1->formalCharge=-1;
              ai1->chemFlag = false;
            } else if(!strcmp(name2,"OD2"))  {
              ai2->formalCharge=-1;
              ai2->chemFlag = false;
            }
            break;
          }
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"OD1")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"OD1")))) 
            order=2;
          break;
        case 'N': /* ASN  */
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"OD1")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"OD1")))) 
            order=2;
          break;
        }
        break;
      case 0:
        if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
          ai1->formalCharge=-1;
          ai1->chemFlag = false;
        } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2"))))  {
          ai2->formalCharge=-1;
          ai2->chemFlag = false;
        }
        if(((!strcmp(name1,"C8"))&&(!strcmp(name2,"N7")))||
           ((!strcmp(name2,"C8"))&&(!strcmp(name1,"N7")))) 
          order=2;
        else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"C5")))||
                ((!strcmp(name2,"C4"))&&(!strcmp(name1,"C5")))) 
          order=2;
          
        else if(((!strcmp(name1,"C6"))&&(!strcmp(name2,"N1")))||
                ((!strcmp(name2,"C6"))&&(!strcmp(name1,"N1")))) 
          order=2;
        else if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"N3")))||
                ((!strcmp(name2,"C2"))&&(!strcmp(name1,"N3")))) 
          order=2;
        else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
          order=2;
        break;
      }
      break;
    case 'C':
      if(resn1[1]==0) {
        if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
          ai1->formalCharge=-1;
          ai1->chemFlag = false;
        } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
          ai2->formalCharge=-1;
          ai2->chemFlag = false;          
        }
        if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"O2")))||
           ((!strcmp(name2,"C2"))&&(!strcmp(name1,"O2")))) 
          order=2;
        else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"N3")))||
                ((!strcmp(name2,"C4"))&&(!strcmp(name1,"N3")))) 
          order=2;
          
        else if(((!strcmp(name1,"C5"))&&(!strcmp(name2,"C6")))||
                ((!strcmp(name2,"C5"))&&(!strcmp(name1,"C6")))) 
          order=2;
        else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
          order=2;
      }
      break;
    case 'D': /* deoxy nucleic acids */
      switch(resn1[1]) {
      case 'A':
        if(resn1[2]==0) {
          if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
            ai1->formalCharge=-1;
            ai1->chemFlag = false;
          } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
            ai2->formalCharge=-1;
            ai2->chemFlag = false;
          }
          if(((!strcmp(name1,"C8"))&&(!strcmp(name2,"N7")))||
             ((!strcmp(name2,"C8"))&&(!strcmp(name1,"N7")))) 
            order=2;
          else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"C5")))||
                  ((!strcmp(name2,"C4"))&&(!strcmp(name1,"C5")))) 
            order=2;
            
          else if(((!strcmp(name1,"C6"))&&(!strcmp(name2,"N1")))||
                  ((!strcmp(name2,"C6"))&&(!strcmp(name1,"N1")))) 
            order=2;
          else if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"N3")))||
                  ((!strcmp(name2,"C2"))&&(!strcmp(name1,"N3")))) 
            order=2;
          else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                  ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
            order=2;
        }
        break;
      case 'C':
        if(resn1[2]==0) {
          if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
            ai1->formalCharge=-1;
            ai1->chemFlag = false;
          } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2"))))  {
            ai2->formalCharge=-1;
            ai2->chemFlag = false;
          }
          if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"O2")))||
             ((!strcmp(name2,"C2"))&&(!strcmp(name1,"O2")))) 
            order=2;
          else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"N3")))||
                  ((!strcmp(name2,"C4"))&&(!strcmp(name1,"N3")))) 
            order=2;
            
          else if(((!strcmp(name1,"C5"))&&(!strcmp(name2,"C6")))||
                  ((!strcmp(name2,"C5"))&&(!strcmp(name1,"C6")))) 
            order=2;
          else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                  ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
            order=2;
        }
        break;
      case 'T':
        if(resn1[2]==0) {
          if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) 
            ai1->formalCharge=-1;
          else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) 
            ai2->formalCharge=-1;
            
          if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"O2")))||
             ((!strcmp(name2,"C2"))&&(!strcmp(name1,"O2")))) 
            order=2;
          else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"O4")))||
                  ((!strcmp(name2,"C4"))&&(!strcmp(name1,"O4")))) 
            order=2;
            
          else if(((!strcmp(name1,"C5"))&&(!strcmp(name2,"C6")))||
                  ((!strcmp(name2,"C5"))&&(!strcmp(name1,"C6")))) 
            order=2;
          else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                  ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
            order=2;
        }
        break;
      case 'G':
        if(resn1[2]==0) {
          if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
            ai1->formalCharge=-1;  
            ai1->chemFlag = false;
          } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
            ai2->formalCharge=-1;
            ai2->chemFlag = false;
          }
          if(((!strcmp(name1,"C6"))&&(!strcmp(name2,"O6")))||
             ((!strcmp(name2,"C6"))&&(!strcmp(name1,"O6")))) 
            order=2;
          else if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"N3")))||
                  ((!strcmp(name2,"C2"))&&(!strcmp(name1,"N3")))) 
            order=2;
          else if(((!strcmp(name1,"C8"))&&(!strcmp(name2,"N7")))||
                  ((!strcmp(name2,"C8"))&&(!strcmp(name1,"N7")))) 
            order=2;
          else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"C5")))||
                  ((!strcmp(name2,"C4"))&&(!strcmp(name1,"C5")))) 
            order=2;
          else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                  ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
            order=2;
        }
        break;
      case 'U':
        if(resn1[2]==0) {
          if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
            ai1->formalCharge=-1;
            ai1->chemFlag = false;
          } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
            ai2->formalCharge=-1;
            ai2->chemFlag = false;
          }
                                      
          if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"O2")))||
             ((!strcmp(name2,"C2"))&&(!strcmp(name1,"O2")))) 
            order=2;
          else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"O4")))||
                  ((!strcmp(name2,"C4"))&&(!strcmp(name1,"O4")))) 
            order=2;
                                      
          else if(((!strcmp(name1,"C5"))&&(!strcmp(name2,"C6")))||
                  ((!strcmp(name2,"C5"))&&(!strcmp(name1,"C6")))) 
            order=2;
          else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                  ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
            order=2;
        }
        break;
      }
      break;
    case 'G':
      switch(resn1[1]) {
      case 'L': 
        switch(resn1[2]) {
        case 'U': /* GLU missing GLUN, GLUH, GLH handling */
          switch(resn1[3]) {                                                
          case 0:
          case 'M': /* minus */
            if(!strcmp(name1,"OE2")) {
              ai1->formalCharge=-1;
              ai1->chemFlag = false;
            } else if(!strcmp(name2,"OE2")) {
              ai2->formalCharge=-1;
              ai2->chemFlag = false;
            }
            break;
          }
          if(((!strcmp(name1,"CD"))&&(!strcmp(name2,"OE1")))||
             ((!strcmp(name2,"CD"))&&(!strcmp(name1,"OE1")))) 
            order=2;
          break;
        case 'N': /* GLN or GLU */
          if(((!strcmp(name1,"CD"))&&(!strcmp(name2,"OE1")))||
             ((!strcmp(name2,"CD"))&&(!strcmp(name1,"OE1")))) 
            order=2;
          break;
        }
        break;
      case 0:
        if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
          ai1->formalCharge=-1;
          ai1->chemFlag = false;
        } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
          ai2->formalCharge=-1;
          ai2->chemFlag = false;
        }
                                    
        if(((!strcmp(name1,"C6"))&&(!strcmp(name2,"O6")))||
           ((!strcmp(name2,"C6"))&&(!strcmp(name1,"O6")))) 
          order=2;
        else if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"N3")))||
                ((!strcmp(name2,"C2"))&&(!strcmp(name1,"N3")))) 
          order=2;
        else if(((!strcmp(name1,"C8"))&&(!strcmp(name2,"N7")))||
                ((!strcmp(name2,"C8"))&&(!strcmp(name1,"N7")))) 
          order=2;
        else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"C5")))||
                ((!strcmp(name2,"C4"))&&(!strcmp(name1,"C5")))) 
          order=2;
        else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
          order=2;
        break;
      }
      break;
    case 'H':
      switch(resn1[1]) {
      case 'I':
        switch(resn1[2]) {
        case 'P':
          if(!strcmp(name1,"ND1")) {
            ai1->formalCharge=1;
            ai1->chemFlag = false;
          } else if(!strcmp(name2,"ND1"))  {
            ai2->formalCharge=1;
            ai2->chemFlag = false;
          }
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD2")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD2")))) 
            order=2;
          else if(((!strcmp(name1,"CE1"))&&(!strcmp(name2,"ND1")))||
                  ((!strcmp(name2,"CE1"))&&(!strcmp(name1,"ND1")))) 
            order=2;
          break;
        case 'S':
          switch(resn1[3]) {
          case 'A': /* HISA Gromacs */
          case 'D': /* HISD Quanta */
            if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD2")))||
               ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD2")))) 
              order=2;
            else if(((!strcmp(name1,"CE1"))&&(!strcmp(name2,"NE2")))||
                    ((!strcmp(name2,"CE1"))&&(!strcmp(name1,"NE2")))) 
              order=2;
            break;
          case 0: /* plain HIS */
          case 'B': /* HISB Gromacs */
          case 'E': /* HISE Quanta */
            if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD2")))||
               ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD2")))) 
              order=2;
            else if(((!strcmp(name1,"CE1"))&&(!strcmp(name2,"ND1")))||
                    ((!strcmp(name2,"CE1"))&&(!strcmp(name1,"ND1")))) 
              order=2;
            break;
          case 'H': /* HISH Gromacs */
          case 'P': /* HISP Quanta */
            if(!strcmp(name1,"ND1")) {
              ai1->formalCharge=1;
              ai1->chemFlag = false;
            } else if(!strcmp(name2,"ND1")) {
              ai2->formalCharge=1;
              ai2->chemFlag = false;
            }
            if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD2")))||
               ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD2")))) 
              order=2;
            else if(((!strcmp(name1,"CE1"))&&(!strcmp(name2,"ND1")))||
                    ((!strcmp(name2,"CE1"))&&(!strcmp(name1,"ND1")))) 
              order=2;
            break;
          }
          break;
        case 'E': /* HIE */
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD2")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD2")))) 
            order=2;
          else if(((!strcmp(name1,"CE1"))&&(!strcmp(name2,"ND1")))||
                  ((!strcmp(name2,"CE1"))&&(!strcmp(name1,"ND1")))) 
            order=2;
          break;
        case 'D': /* HID */
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD2")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD2")))) 
            order=2;
          else if(((!strcmp(name1,"CE1"))&&(!strcmp(name2,"NE2")))||
                  ((!strcmp(name2,"CE1"))&&(!strcmp(name1,"NE2")))) 
            order=2;
          break;
        }
        break;
      }
      break;
    case 'I':
      if(resn1[1]==0) {
        if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
          ai1->formalCharge=-1;
          ai1->chemFlag = false;
        } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
          ai2->formalCharge=-1;
          ai2->chemFlag = false;
        }
        if(((!strcmp(name1,"C8"))&&(!strcmp(name2,"N7")))||
           ((!strcmp(name2,"C8"))&&(!strcmp(name1,"N7")))) 
          order=2;
        else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"C5")))||
                ((!strcmp(name2,"C4"))&&(!strcmp(name1,"C5")))) 
          order=2;
                                              
        else if(((!strcmp(name1,"C6"))&&(!strcmp(name2,"N1")))||
                ((!strcmp(name2,"C6"))&&(!strcmp(name1,"N1")))) 
          order=2;
        else if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"N3")))||
                ((!strcmp(name2,"C2"))&&(!strcmp(name1,"N3")))) 
          order=2;
        else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
          order=2;
      }
      break;
    case 'P':
      switch(resn1[1]) {
      case 'H': /* PHE */
        if(resn1[2]=='E') {
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD1")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD1")))) 
            order=2;
          else if(((!strcmp(name1,"CZ"))&&(!strcmp(name2,"CE1")))||
                  ((!strcmp(name2,"CZ"))&&(!strcmp(name1,"CE1")))) 
            order=2;
                                                
          else if(((!strcmp(name1,"CE2"))&&(!strcmp(name2,"CD2")))||
                  ((!strcmp(name2,"CE2"))&&(!strcmp(name1,"CD2")))) 
            order=2;
          break; 
        }
      }
      break;
    case 'L':
      switch(resn1[1]) {
      case 'Y':
        switch(resn1[2]) {
        case 'S': /* LYS. */
          switch(resn1[3]) {                                                
          case 0:
          case 'P': /* LYS, LYSP */
            if(!strcmp(name1,"NZ"))  {
              ai1->formalCharge=1;
              ai1->chemFlag = false;
            } else if(!strcmp(name2,"NZ")) {
              ai2->formalCharge=1;
              ai2->chemFlag = false;
            }
            break;
          }
          break;
        }
        break;
      }
      break;
    case 'T':
      switch(resn1[1]) {
      case 'Y': /* TYR */
        if(resn1[2]=='R') {
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD1")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD1")))) 
            order=2;
          else if(((!strcmp(name1,"CZ"))&&(!strcmp(name2,"CE1")))||
                  ((!strcmp(name2,"CZ"))&&(!strcmp(name1,"CE1")))) 
            order=2;
                                                
          else if(((!strcmp(name1,"CE2"))&&(!strcmp(name2,"CD2")))||
                  ((!strcmp(name2,"CE2"))&&(!strcmp(name1,"CD2")))) 
            order=2;
          break; 
        }
        break;
      case 'R':
        if(resn1[2]=='P') {
          if(((!strcmp(name1,"CG"))&&(!strcmp(name2,"CD1")))||
             ((!strcmp(name2,"CG"))&&(!strcmp(name1,"CD1")))) 
            order=2;
          else if(((!strcmp(name1,"CZ3"))&&(!strcmp(name2,"CE3")))||
                  ((!strcmp(name2,"CZ3"))&&(!strcmp(name1,"CE3")))) 
            order=2;
          else if(((!strcmp(name1,"CZ2"))&&(!strcmp(name2,"CH2")))||
                  ((!strcmp(name2,"CZ2"))&&(!strcmp(name1,"CH2")))) 
            order=2;
          else if(((!strcmp(name1,"CE2"))&&(!strcmp(name2,"CD2")))||
                  ((!strcmp(name2,"CE2"))&&(!strcmp(name1,"CD2")))) 
            order=2;
          break; 
        }
        break;
      case 0:
        if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
          ai1->formalCharge=-1;
          ai1->chemFlag = false;
        } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
          ai2->formalCharge=-1;
          ai2->chemFlag = false;
        }
                                    
        if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"O2")))||
           ((!strcmp(name2,"C2"))&&(!strcmp(name1,"O2")))) 
          order=2;
        else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"O4")))||
                ((!strcmp(name2,"C4"))&&(!strcmp(name1,"O4")))) 
          order=2;
                                              
        else if(((!strcmp(name1,"C5"))&&(!strcmp(name2,"C6")))||
                ((!strcmp(name2,"C5"))&&(!strcmp(name1,"C6")))) 
          order=2;
        else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
          order=2;
        break;
      }
      break;
    case 'U':
      if(resn1[1]==0) {
        if(((!strcmp(name1,"O2P"))||(!strcmp(name1,"OP2")))) {
          ai1->formalCharge=-1;
          ai1->chemFlag = false;
        } else if(((!strcmp(name2,"O2P"))||(!strcmp(name2,"OP2")))) {
          ai2->formalCharge=-1;
          ai2->chemFlag = false;
        }
        if(((!strcmp(name1,"C2"))&&(!strcmp(name2,"O2")))||
           ((!strcmp(name2,"C2"))&&(!strcmp(name1,"O2")))) 
          order=2;
        else if(((!strcmp(name1,"C4"))&&(!strcmp(name2,"O4")))||
                ((!strcmp(name2,"C4"))&&(!strcmp(name1,"O4")))) 
          order=2;
                                              
        else if(((!strcmp(name1,"C5"))&&(!strcmp(name2,"C6")))||
                ((!strcmp(name2,"C5"))&&(!strcmp(name1,"C6")))) 
          order=2;
        else if(((!strcmp(name1,"P"))&&(((!strcmp(name2,"O1P"))||(!strcmp(name2,"OP1")))))||
                ((!strcmp(name2,"P"))&&(((!strcmp(name1,"O1P"))||(!strcmp(name1,"OP1")))))) 
          order=2;
      }
      break;
    }
  }
  *(bond_order) = order;
}

void ObjectMoleculeFixChemistry(ObjectMolecule *I, int sele1, int sele2, int invalidate) 
{
  int b;
  int flag = false;
  int s1,s2;
  AtomInfoType *ai1,*ai2;
  int order;
  BondType *bond;
  bond = I->Bond;
  for(b=0;b<I->NBond;b++) {
    flag = false;
    ai1 = I->AtomInfo + bond->index[0];
    ai2 = I->AtomInfo + bond->index[1];
    s1=ai1->selEntry;
    s2=ai2->selEntry;
    
    if((SelectorIsMember(I->Obj.G,s1,sele1)&&
        SelectorIsMember(I->Obj.G,s2,sele2))||
       (SelectorIsMember(I->Obj.G,s2,sele1)&&
        SelectorIsMember(I->Obj.G,s1,sele2))) {
      order = -1;
      if(!ai1->resn[3]) { /* Standard disconnected PDB residue */
        if(AtomInfoSameResidue(I->Obj.G,ai1,ai2)) {
          assign_pdb_known_residue(I->Obj.G, ai1,ai2,&order);
        }
      }
      if(order>0) {
        bond->order = order;
        ai1->chemFlag=false;
        ai2->chemFlag=false;
        flag = true;
      } else if(invalidate) {
        ai1->chemFlag=false;
        ai2->chemFlag=false;
        flag=true;
      }
    }
    bond++;
  }
  if(flag) {
    ObjectMoleculeInvalidate(I,cRepAll,cRepInvAll,-1);
    SceneChanged(I->Obj.G);
  }
}
  


void ObjMolPairwiseInit(ObjMolPairwise *pairwise)
{
  UtilZeroMem((char*)pairwise,sizeof(ObjMolPairwise));
  pairwise->trg_vla = VLAlloc(int,10);
  pairwise->mbl_vla = VLAlloc(int,10);
}

void ObjMolPairwisePurge(ObjMolPairwise *pairwise)
{
  VLAFreeP(pairwise->trg_vla);
  VLAFreeP(pairwise->mbl_vla);
}

int ObjectMoleculeConvertIDsToIndices(ObjectMolecule *I,int *id,int n_id)
{
  /* return true if all IDs are unique, false if otherwise */

  int min_id,max_id,range,*lookup = NULL;
  int unique = true;

  /* this routine only works if IDs cover a reasonable range --
     should rewrite using a hash table */

  if(I->NAtom) {

    /* determine range */

    {
      int a,cur_id;
      cur_id = I->AtomInfo[0].id;
      min_id = cur_id;
      max_id = cur_id;
      for(a=1;a<I->NAtom;a++) {
        cur_id = I->AtomInfo[a].id;
        if(min_id>cur_id) min_id = cur_id;
        if(max_id<cur_id) max_id = cur_id;
      }
    }

    /* create cross-reference table */

    {
      int a,offset;
      
      range = max_id - min_id + 1;
      lookup = Calloc(int,range);
      for(a=0;a<I->NAtom;a++) {
        offset = I->AtomInfo[a].id - min_id;
        if(!lookup[offset])
          lookup[offset] = a+1;
        else 
          unique = false;
      }
    }
    
    /* iterate through IDs and replace with indices or -1 */

    {
      int i,offset,lkup;

      for(i=0;i<n_id;i++) {
        offset = id[i]-min_id;
        if((offset>=0)&&(offset<range)) {
          lkup = lookup[offset];
          if(lkup>0) {
            id[i] = lkup-1;
          } else {
            id[i] = -1; /* negative means no match */
          }
        } else 
          id[i] = -1;
      }
    }
  }
  
  FreeP(lookup);
  return unique;
    
}
static char *check_next_pdb_object(char *p, int skip_to_next)
{
  char *start = p;
  char *line_start;
  char cc[MAXLINELEN];  
  while(*p) {
    line_start = p;
    p=ncopy(cc,p,6);
    if((cc[0]=='H')&&
       (cc[1]=='E')&&
       (cc[2]=='A')&&
       (cc[3]=='D')&&
       (cc[4]=='E')&&
       (cc[5]=='R')) {
      if(skip_to_next) 
        return line_start;
      else
        return start;
    } else if(((cc[0]=='A')&& /* ATOM */
               (cc[1]=='T')&&
               (cc[2]=='O')&&
               (cc[3]=='M')&&
               (cc[4]==' ')&&
               (cc[5]==' '))||
              ((cc[0]=='H')&& /* HETATM */
               (cc[1]=='E')&&
               (cc[2]=='T')&&
               (cc[3]=='A')&&
               (cc[4]=='T')&&
               (cc[5]=='M'))) {
      p=nskip(p,5); 
      ParseNTrim(cc,p,14); 
      /* this is a special workaround for the bogus HETATM entry in PDB 4ZN2:
         HETATM20829              0       0.000   0.000   0.000  1.00  0.00      4ZN20773
         which screws up our PDB test case */
      if(cc[0]) {
        return start;
      }
    } else if(((cc[0]=='E')&& /* if we pass over the END of the current PDB file, then reset start */
               (cc[1]=='N')&&
               (cc[2]=='D'))) {
      if(strcmp("END",cc)==0) { /* END */
        if(skip_to_next) {
          start = line_start;
        }
      }
    }
    p=nextline(p);
  }
  return NULL;
}

static int get_multi_object_status(char *p) /* expensive -- only call
                                               this if there is no
                                               other way to determine
                                               this information */
{
  int seen_header = 0;
  char cc[MAXLINELEN];  
  while(*p) {
    p=ncopy(cc,p,6);
    if((cc[0]=='H')&&
       (cc[1]=='E')&&
       (cc[2]=='A')&&
       (cc[3]=='D')&&
       (cc[4]=='E')&&
       (cc[5]=='R')) {
      if(seen_header) {
        return 1;
      } else {
        seen_header = true;
      }
    }
    p=nextline(p);
  }
  return -1;
}

int ObjectMoleculeAutoDisableAtomNameWildcard(ObjectMolecule *I)
{
  PyMOLGlobals *G = I->Obj.G;
  char wildcard = 0;
  int found_wildcard = false;

  {
    char *tmp = SettingGet_s(G, NULL, I->Obj.Setting,cSetting_atom_name_wildcard);
    if(tmp && tmp[0]) {
      wildcard = *tmp;
    } else {
      tmp = SettingGet_s(G, NULL, I->Obj.Setting,cSetting_wildcard);      
      if(tmp) {
        wildcard = *tmp;
      }
    }
    if(wildcard==32)
      wildcard = 0;

  }

  if(wildcard) {
    register int a;
    register char *p, ch;
    register AtomInfoType *ai = I->AtomInfo;

    for(a=0;a<I->NAtom;a++) {
      p = ai->name;
      while( (ch=*(p++)) ) {
        if(ch == wildcard) {
          found_wildcard = true;
          break;
        }
      }
      ai++;
    }
    if(found_wildcard) {
      ExecutiveSetObjSettingFromString(G,cSetting_atom_name_wildcard," ",
                                       &I->Obj,-1, true, true);
    }
  }
  return found_wildcard;
}
/*========================================================================*/
#define PDB_MAX_TAGS 64

CoordSet *ObjectMoleculePDBStr2CoordSet(PyMOLGlobals *G,
                                        char *buffer,
                                        AtomInfoType **atInfoPtr,
                                        char **restart_model,
                                        char *segi_override,
                                        M4XAnnoType *m4x,
                                        char *pdb_name,
                                        char **next_pdb,
                                        PDBInfoRec *info,
                                        int quiet, int *model_number)
{

  register char *p;
  int nAtom;
  int a,b,c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL,*ai;
  int AFlag;
  char SSCode;
  int atomCount;
  int bondFlag = false;
  BondType *bond=NULL,*ii1,*ii2;
  int *idx;
  int nBond=0;
  int b1,b2,nReal,maxAt;
  CSymmetry *symmetry = NULL;
  int auto_show_lines = (int)SettingGet(G,cSetting_auto_show_lines);
  int auto_show_nonbonded = (int)SettingGet(G,cSetting_auto_show_nonbonded);
  int auto_show_spheres = (int)SettingGet(G,cSetting_auto_show_spheres);
  int reformat_names = (int)SettingGet(G,cSetting_pdb_reformat_names_mode);
  int truncate_resn = SettingGetGlobal_b(G,cSetting_pdb_truncate_residue_name);
  char *tags_in = SettingGetGlobal_s(G,cSetting_pdb_echo_tags), *tag_start[PDB_MAX_TAGS];
  int n_tags = 0;
  int foundNextModelFlag = false;
  int ssFlag = false;
  int ss_resv1=0,ss_resv2=0;
  ResIdent ss_resi1="",ss_resi2="";
  unsigned char ss_chain1=0,ss_chain2=0;
  SSEntry *ss_list = NULL;
  int n_ss = 1;
  int *(ss[256]); /* one array for each chain identifier */

  char cc[MAXLINELEN],tags[MAXLINELEN];
  char cc_saved,ctmp;
  int index;
  int ignore_pdb_segi = 0;
  int ss_valid,ss_found=false;
  SSEntry *sst;
  int ssi = 0;
  int only_read_one_model = false;
  int ignore_conect = false;
  int have_bond_order = false;
  int seen_model,in_model = false;
  int seen_conect = false;
  int have_conect = false;
  int is_end_of_object = false;
  int literal_names = SettingGetGlobal_b(G,cSetting_pdb_literal_names);
  int bogus_name_alignment = true;
  AtomName literal_name = "";

  if(tags_in&&(!quiet)&&(!*restart_model)) {
    char *p = tags;
    strcpy(tags,tags_in);
    
    while(*p) {
      while(*p==' ') /* skip spaces */
        p++;
      if(*p) {
        tag_start[n_tags] = p;
        n_tags++;
        while(*p) {
          if(*p!=',') {
            if(*p==' ')
              *p = 0;
            p++;
          } else
            break;
        }
        if(*p) { /* terminate tag */
          *p = 0;
          p++;
        }
      }
    }
  }

  if(literal_names)
    reformat_names = 0;
  
  ignore_pdb_segi = (int)SettingGet(G,cSetting_ignore_pdb_segi);

  p=buffer;
  nAtom=0;
  if(atInfoPtr)
	 atInfo = *atInfoPtr;

  if(!atInfo)
    ErrFatal(G,"PDBStr2CoordSet","need atom information record!"); /* failsafe for old version..*/

  if(buffer == *restart_model)
    only_read_one_model = true;
  else if(info && (info->multiplex>0)) {
    only_read_one_model = true;    
    if(!info->multi_object_status) { /* figure out if this is a multi-object (multi-HEADER) pdb file */
      info->multi_object_status = get_multi_object_status(p);
    }
  }
  /* PASS 1 */
  {
    int seen_end_of_atoms = false;
    
    *restart_model = NULL;
    while(*p)
      {
        if(n_tags && !quiet) {
          int skip=false;
          
          if((p[0]=='H')&& 
             (p[1]=='E')&&
             (p[2]=='A')&&
             (p[3]=='D')&&
             (p[4]=='E')&&
             (p[5]=='R')) {
            if(nAtom>0) { 
              /* don't print HEADER until next time*/
              skip=true;
            }
          }
          if(!skip) {
            /* fast unrolled string match */
            register int tc = 0;
            register char *q;
            register int same;
            while(tc<n_tags) {
              same = true;
              q = tag_start[tc];
              if(p[0] != q[0])
                same = false;
              else if(p[0]&&q[0]) {
                if((p[1] != q[1])&&!((p[1]==' ')&&!q[1]))
                  same = false;
                else if(p[1]&&q[1]) {
                  if((p[2] != q[2])&&!((p[2]==' ')&&!q[2]))
                    same = false;
                  else if(p[3]&&q[3]) {
                    if((p[3] != q[3])&&!((p[3]==' ')&&!q[3]))
                      same = false;
                    else if(p[4]&&q[4]) {
                      if((p[4] != q[4])&&!((p[4]==' ')&&!q[4]))
                        same = false;
                      else if(p[5]&&q[5]) {
                        if((p[5] != q[5])&&!((p[5]==' ')&&!q[5]))
                          same = false;
                      }
                    }
                  }
                }
              }
              if(same) {
                ParseNTrimRight(cc,p,MAXLINELEN-1);
                /*              OrthoAddOutput(G," PDB: ");*/
                OrthoAddOutput(G,cc);
                OrthoNewLine(G,NULL,true);
              }
              tc++;
            }
          }
        }
        if(((p[0]== 'A')&&(p[1]=='T')&&(p[2]=='O')&&(p[3]=='M'))|| /* ATOM */
           ((p[0]== 'H')&&(p[1]=='E')&&(p[2]=='T')&& /* HETATM */
            (p[3]=='A')&&(p[4]=='T')&&(p[5]=='M')&&(!*restart_model))) {
          if(!seen_end_of_atoms) 
            nAtom++;
          if(bogus_name_alignment) {
            ncopy(cc,nskip(p,12),4); /* copy the atom field */
            if((cc[0]==32)&&(cc[1]!=32)) { /* check to see if indentation was followed correctly */
              bogus_name_alignment = false;
            }
          }
        } else if((p[0]== 'H')&&(p[1]=='E')&&(p[2]=='L')&&(p[3]=='I')&&(p[4]=='X')) /* HELIX */
          ssFlag=true;
        else if((p[0]== 'S')&&(p[1]=='H')&&(p[2]=='E')&&(p[3]=='E')&&(p[4]=='T')) /* SHEET */
          ssFlag=true;
        else if((p[0]=='H')&& /* HEADER */
                (p[1]=='E')&&
                (p[2]=='A')&&
                (p[3]=='D')&&
                (p[4]=='E')&&
                (p[5]=='R'))
          {
            if(nAtom>0) { /* if we've already found atom records, then this must be a new pdb */
              (*next_pdb) = p;
              break;
            }
          }
        else if((p[0]== 'R')&&(p[1]=='E')&&(p[2]=='M')&& /* REMARK */
                (p[3]=='A')&&(p[4]=='R')&&(p[5]=='K')) {
          ntrim(cc,p,30);
          if(strcmp("REMARK    GENERATED BY TRJCONV",cc)==0) 
            if(info) info->ignore_header_names = true;
        } else if((p[0]== 'E')&&(p[1]=='N')&&(p[2]=='D')&& /* ENDMDL */
                  (p[3]=='M')&&(p[4]=='D')&&(p[5]=='L')&&(!*restart_model)) {
          *restart_model=nextline(p);
          seen_end_of_atoms = true;
          if(only_read_one_model) 
            break;
        } else if((p[0]== 'E')&&(p[1]=='N')&&(p[2]=='D')) { /* stop parsing after END */
          ntrim(cc,p,6);
          if(strcmp("END",cc)==0) { /* END */
            seen_end_of_atoms = true;
            if(next_pdb) {
              p = nextline(p);
              ncopy(cc,p,6);
              if(strcmp("HEADER",cc)==0) {
                (*next_pdb) = p; /* found another PDB file after this one... */
              } else if(strcmp("ENDMDL",cc)==0) {
                seen_end_of_atoms = false;
              }
            } 
            break;
          }
        } else if((p[0]== 'C')&&(p[1]=='O')&&(p[2]=='N')&&
                  (p[3]=='E')&&(p[4]=='C')&&(p[5]=='T')) { /* CONECT */
          have_conect = true;
          bondFlag=true;
        } else if((p[0]== 'U')&&(p[1]=='S')&&(p[2]=='E')&&
                (p[3]=='R')&&(!*restart_model)) {

          /* Metaphorics key 'USER     '*/
          if((p[4]==' ')&&(p[5]==' ')&&(p[6]==' ')&&
             (p[7]==' ')&&(p[8]==' ')&&m4x) {
            p = nskip(p,10);
            p = ntrim(cc,p,6);
            m4x->annotated_flag = true;
            switch(cc[0]) {
            case 'H':
              if(WordMatchExact(G,"HINT",cc,true)) {
                p = nskip(p,1);
                p = ntrim(cc,p,6); /* get context name */
                if(WordMatchExact(G,"ALIGN",cc,true)) { /* ALIGN is special */
                  if(!m4x->align) {
                    m4x->align=Calloc(M4XAlignType,1);
                    M4XAlignInit(m4x->align);
                    p = nskip(p,8);
                    p = ntrim(cc,p,6); /* get visibility of this structure */
                  }
                } else if(WordMatchExact(G,"HIDE",cc,true)) {
                  m4x->invisible = 1;
                } else {
                  if(!m4x->context) {
                    m4x->context = VLACalloc(M4XContextType,10);
                  } 
                  if(m4x->context) {
                    int cn;
                    int found=false;
                  
                    /* does context already exist ? */
                    for(cn=0;cn<m4x->n_context;cn++) {
                      if(WordMatchExact(G,m4x->context[cn].name,cc,true)) {
                        found=true;
                        break;
                      }
                    }
                  
                    /* if not, then create it */
                    if(!found) {
                      cn = m4x->n_context++;
                      VLACheck(m4x->context,M4XContextType,cn);
                      UtilNCopy(m4x->context[cn].name,cc,sizeof(WordType));
                    }
                  
                    while(*cc) {
                      p = nskip(p,1);
                      p = ntrim(cc,p,6);
                      switch(cc[0]) {
                      case 'B':
                        if(WordMatchExact(G,"BORDER",cc,true)) {
                          /* ignore PDB CONECT if BORDER present */
                          ignore_conect = true;
                          have_bond_order = true;
                          bondFlag = true;
                        }
                        break;
                      case 'S':
                        if(WordMatchExact(G,"SITE",cc,true)) {
                          if(!m4x->context[cn].site) {
                            m4x->context[cn].site=VLAlloc(int,50);
                          }
                        } 
                        break;
                      case 'L':
                        if(WordMatchExact(G,"LIGAND",cc,true)) {
                          if(!m4x->context[cn].ligand) {
                            m4x->context[cn].ligand=VLAlloc(int,50);
                          }
                        }
                        break;
                      case 'W':
                        if(WordMatchExact(G,"WATER",cc,true)) {
                          if(!m4x->context[cn].water) {
                            m4x->context[cn].water=VLAlloc(int,50);
                          }
                        } 
                        break;
                      case 'H':
                        if(WordMatchExact(G,"HBOND",cc,true)) {
                          if(!m4x->context[cn].hbond) {
                            m4x->context[cn].hbond=VLAlloc(M4XBondType,50);
                          }
                        }
                        break;
                      case 'N':
                        if(WordMatchExact(G,"NBOND",cc,true)) {
                          if(!m4x->context[cn].nbond) {
                            m4x->context[cn].nbond=VLAlloc(M4XBondType,50);
                          }
                        }
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
        }
        p=nextline(p);
      }
  }

  *restart_model=NULL;
  coord=VLAlloc(float,3*nAtom);

  if(atInfo)
	 VLACheck(atInfo,AtomInfoType,nAtom);

  if(bondFlag) {
    nBond=0;
    bond=VLACalloc(BondType,6*nAtom);  
  }
  p=buffer;
  PRINTFB(G,FB_ObjectMolecule,FB_Blather)
	 " ObjectMoleculeReadPDB: Found %i atoms...\n",nAtom
    ENDFB(G);

  if(ssFlag) {
    for(a=0;a<=255;a++) {
      ss[a]=0;
    }
    ss_list=VLAlloc(SSEntry,50);
  }

  a=0; /* WATCHOUT */
  atomCount=0;

  /* PASS 2 */
  seen_model = false;

  while(*p) {
    /*      printf("%c%c%c%c%c%c\n",p[0],p[1],p[2],p[3],p[4],p[5]);*/
    AFlag=false;
    SSCode=0;
    if((p[0]=='A')&& /* ATOM */
       (p[1]=='T')&&
       (p[2]=='O')&&
       (p[3]=='M'))
      AFlag = 1;
    else if((p[0]=='H')&& /* HETATM */
            (p[1]=='E')&&
            (p[2]=='T')&&
            (p[3]=='A')&&
            (p[4]=='T')&&
            (p[5]=='M'))
      AFlag = 2;
    else if((p[0]=='H')&& /* HEADER */
            (p[1]=='E')&&
            (p[2]=='A')&&
            (p[3]=='D')&&
            (p[4]=='E')&&
            (p[5]=='R')) {
      
      if(pdb_name) {
        if(atomCount>0) {
          /* have we already read atoms?  then this is probably a new PDB file */
          (*next_pdb) = p; /* found another PDB file after this one... */
          break;
        } else if((!info)||(!info->ignore_header_names)) { 
          /* if this isn't an MD trajectory... */
          char *pp;
          pp=nskip(p,62); /* is there a four-letter PDB code? */
          pp=ntrim(pdb_name,pp,4);
          if(pdb_name[0]==0) { /* if not, is there a plain name (for MERCK!)*/
            p=nskip(p,6);
            p=ntrim(cc,p,44);
            UtilNCopy(pdb_name,cc,WordLength);
          } else {
            p=pp;
          }
        }
      }
    } else if((p[0]=='M')&& /* MODEL */
              (p[1]=='O')&&
              (p[2]=='D')&&
              (p[3]=='E')&&
              (p[4]=='L')) {
      if(model_number) {
        int tmp;
        p = nskip(p,10);
        p = ncopy(cc,p,5);
        if(sscanf(cc,"%d",&tmp)==1)
          *model_number = tmp;
      }
      seen_model = true;
      in_model = true;
    } else if((p[0]=='H')&& /* HELIX */
              (p[1]=='E')&&
              (p[2]=='L')&&
              (p[3]=='I')&&
              (p[4]=='X')) {
      ss_valid=true;
      
      p=nskip(p,19);
      ss_chain1 = (*p);
      p=nskip(p,2);
      p=ncopy(cc,p,4);
      if(!sscanf(cc,"%s",ss_resi1)) ss_valid=false;
      if(!sscanf(cc,"%d",&ss_resv1)) ss_valid=false;
      
      p=nskip(p,6);
      ss_chain2 = (*p);
      p=nskip(p,2);
      p=ncopy(cc,p,4);
      
      if(!sscanf(cc,"%s",ss_resi2)) ss_valid=false;
      if(!sscanf(cc,"%d",&ss_resv2)) ss_valid=false;
      
      if(ss_valid) {
        PRINTFB(G,FB_ObjectMolecule,FB_Blather)
          " ObjectMolecule: read HELIX %c %s %c %s\n",
          ss_chain1,ss_resi1,ss_chain2,ss_resi2
          ENDFB(G);
        SSCode='H';
      }
      
      if(ss_chain1==' ') ss_chain1=0;
      if(ss_chain2==' ') ss_chain2=0;          
      
    } else if((p[0]=='S')&& /* SHEET */
              (p[1]=='H')&&
              (p[2]=='E')&&
              (p[3]=='E')&&
              (p[4]=='T')) {
      ss_valid=true;
      p=nskip(p,21);
      ss_chain1 = (*p);
      p=nskip(p,1);
      p=ncopy(cc,p,4);
      if(!sscanf(cc,"%s",ss_resi1)) ss_valid=false;
      if(!sscanf(cc,"%d",&ss_resv1)) ss_valid=false;
      p=nskip(p,6);
      ss_chain2 = (*p);
      p=nskip(p,1);
      p=ncopy(cc,p,4);
      if(!sscanf(cc,"%s",ss_resi2)) ss_valid=false;
      if(!sscanf(cc,"%d",&ss_resv2)) ss_valid=false;
      
      if(ss_valid) {
        PRINTFB(G,FB_ObjectMolecule,FB_Blather)
          " ObjectMolecule: read SHEET %c %s %c %s\n",
          ss_chain1,ss_resi1,ss_chain2,ss_resi2
          ENDFB(G);
        SSCode = 'S';
      }
      
      if(ss_chain1==' ') ss_chain1=0;
      if(ss_chain2==' ') ss_chain2=0;   
      
    } else if((p[0]=='E')&& /* ENDMDL */
              (p[1]=='N')&&
              (p[2]=='D')&&
              (p[3]=='M')&&
              (p[4]=='D')&&
              (p[5]=='L')) {
      if(*restart_model)
        in_model = false;
      else {
        *restart_model=nextline(p);
        in_model = false;
        if(only_read_one_model) {
          char *pp;
          pp = nextline(p); 
          if((pp[0]=='E')&& /* END we're going to be starting a new object...*/
             (pp[1]=='N')&&
             (pp[2]=='D')) {
            (*next_pdb) = check_next_pdb_object(nextline(pp),true);
            if(info && (info->multiplex == 0)) {
              /* multiplex == 0:  FORCED multimodel behavior with concatenated PDB files */
              (*restart_model) = (*next_pdb);
              (*next_pdb) = NULL;
              foundNextModelFlag = true;
              info->multi_object_status = -1;
            } else {
              is_end_of_object = true;
            }
          } else if((pp[0]=='M')&& /* not a new object...just a new state (model) */
                    (pp[1]=='O')&&
                    (pp[2]=='D')&&
                    (pp[3]=='E')&&
                    (pp[4]=='L')) {
            if(info && (info->multiplex>0)) { /* end object if we're multiplexing */
              (*next_pdb) = check_next_pdb_object(pp, true);
              (*restart_model) = NULL;
            } else 
              is_end_of_object = false;
          } else {
            if(pp[0]>32) /* more content follows... */
              (*next_pdb) = check_next_pdb_object(pp, true);
            else
              (*next_pdb) = NULL; /* at end of file */
          }
          break;
        }
      }
    } else if((p[0]=='E')&& /* END */
              (p[1]=='N')&&
              (p[2]=='D')) {
      ntrim(cc,p,6);
      if((strcmp("END",cc)==0)&&(!in_model)) { /* stop parsing here...*/
        char *pp;
        pp=nextline(p); /* ...unless we're in MODEL or next line is itself ENDMDL */
        p=ncopy(cc,p,6);
        if(!((cc[0] =='E')&&
             (cc[1] =='N')&&
             (cc[2] =='D')&&
             (cc[3] =='M')&&
             (cc[4] =='D')&&
             (cc[5] =='L'))) { /* NOTE: this test seems unnecessary given strcmp above...*/
          if(!*next_pdb) {
            (*next_pdb) = check_next_pdb_object(pp, false);
          }
          if((*next_pdb) && info && (!info->multiplex) && !(*restart_model)) {
            /* multiplex == 0:  FORCED multimodel behavior with concatenated PDB files */
            (*restart_model) = (*next_pdb);
            (*next_pdb) = NULL;
            foundNextModelFlag = true;
            info->multi_object_status = -1; 
            is_end_of_object = false;
            break;
          }
          if(*next_pdb) { /* we've found another object... */
            if(*restart_model) 
              is_end_of_object = false; /* however, if we're parsing multi-models, then we're not yet at the end */
            else
              is_end_of_object = true;
            break;
          } else if(!seen_model)
            break;
        }
      }
    } else if((p[0]=='C')&& /* CRYST1 */
              (p[1]=='R')&&
              (p[2]=='Y')&&
              (p[3]=='S')&&
              (p[4]=='T')&&
              (p[5]=='1')&& (!*restart_model))      {
      if(!symmetry) symmetry=SymmetryNew(G);          
      if(symmetry) {
        int symFlag=true;
        PRINTFB(G,FB_ObjectMolecule,FB_Blather)
          " PDBStrToCoordSet: Attempting to read symmetry information\n"
          ENDFB(G);
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
        p=ncopy(symmetry->SpaceGroup,p,11);
        UtilCleanStr(symmetry->SpaceGroup);
        p=ncopy(cc,p,3);
        if(sscanf(cc,"%d",&symmetry->PDBZValue)!=1) symmetry->PDBZValue=1;
        if(!symFlag) {
          ErrMessage(G,"PDBStrToCoordSet","Error reading CRYST1 record\n");
          SymmetryFree(symmetry);
          symmetry=NULL;
        }
      }
    } else if((p[0]=='S')&& /* SCALEn */
                (p[1]=='C')&&
                (p[2]=='A')&&
                (p[3]=='L')&&
                (p[4]=='E')&&
              (!*restart_model)&&info) {
      switch(p[5]) {
      case '1':
      case '2':
      case '3':
        {
          int flag=(p[5]-'1');
          int offset=flag*4;
          int scale_flag=true;
          p=nskip(p,10);
          p=ncopy(cc,p,10);
          if(sscanf(cc,"%f",&info->scale.matrix[offset])!=1) scale_flag=false;
          p=ncopy(cc,p,10);
          if(sscanf(cc,"%f",&info->scale.matrix[offset+1])!=1) scale_flag=false;
          p=ncopy(cc,p,10);
          if(sscanf(cc,"%f",&info->scale.matrix[offset+2])!=1) scale_flag=false;
          p=nskip(p,5);
          p=ncopy(cc,p,10);
          if(sscanf(cc,"%f",&info->scale.matrix[offset+3])!=1) scale_flag=false;                
          if(scale_flag)
            info->scale.flag[flag]=true;
          PRINTFB(G,FB_ObjectMolecule,FB_Blather)
            " PDBStrToCoordSet: SCALE%d %8.4f %8.4f %8.4f %8.4f\n",flag+1,
            info->scale.matrix[offset],
            info->scale.matrix[offset+1],
            info->scale.matrix[offset+2],
            info->scale.matrix[offset+3]
            ENDFB(G);
        }
        break;
      }
    } else if((p[0]=='C')&& /* CONECT */
              (p[1]=='O')&&
              (p[2]=='N')&&
              (p[3]=='E')&&
              (p[4]=='C')&&
              (p[5]=='T')&&
              bondFlag && (!ignore_conect) &&
              ((!*restart_model)||(!in_model))) {
      seen_conect = true;
      p=nskip(p,6);
      p=ncopy(cc,p,5);
      if(sscanf(cc,"%d",&b1)==1)
        while (1) {
          p=ncopy(cc,p,5);
          if(sscanf(cc,"%d",&b2)!=1)
            break;
          else {
            if((b1>=0)&&(b2>=0)&&(b1!=b2)) { /* IDs must be positive and distinct*/
              VLACheck(bond,BondType,nBond);
              if(b1<=b2) {
                bond[nBond].index[0]=b1; /* temporarily store the atom indexes */
                bond[nBond].index[1]=b2;
                bond[nBond].order=1;
                bond[nBond].stereo=0;
              } else {
                bond[nBond].index[0]=b2;
                bond[nBond].index[1]=b1;
                bond[nBond].order=1;
                bond[nBond].stereo=0;
              }
              nBond++;
            }
          }
        }
    } else if((p[0]=='U')&& /* USER */
              (p[1]=='S')&&
              (p[2]=='E')&&
              (p[3]=='R')&&
              (!*restart_model)) {
      /* Metaphorics key 'USER     ' */
      if((p[4]==' ')&&
         (p[5]==' ')&&
         (p[6]==' ')&&
         
         (p[7]==' ')&&
         (p[8]==' ')&&
         m4x) {
        
        int parsed_flag = false;
        
        p = nskip(p,10);
        p = ntrim(cc,p,6);
        
        /* is this a context name or a USER record? */
        switch(cc[0]) {
        case 'X':
          if(WordMatchExact(G,"XNAME",cc,true)) {  /* object name */
            p=nskip(p,1);
            p=ntrim(m4x->xname,p,10);
            if(m4x->xname[0]) {
              m4x->xname_flag = true;
              parsed_flag = true;
            }
          }
          break;
        case 'A': /* alignment information */
          if(WordMatchExact(G,"ALIGN",cc,true)) {
            if(m4x->align && m4x->align->id_at_point) {
              M4XAlignType *align = m4x->align;
              char target[11];
              int atom_id,point_id;
              float fitness;
              p=nskip(p,1);
              p=ncopy(cc,p,6);
              if(sscanf(cc,"%d",&atom_id)==1) {
                p=nskip(p,1);
                p=ntrim(target,p,10);
                if(target[0]) {
                  if(!align->target[0]) 
                    UtilNCopy(align->target,target,WordLength);
                  if(WordMatchExact(G,align->target,target,true)) { /* must match the one target allowed */
                    p=nskip(p,1);
                    p=ncopy(cc,p,6);
                    if(sscanf(cc,"%d",&point_id)==1) {                      
                      p=nskip(p,1);
                      p=ncopy(cc,p,6);
                      if(sscanf(cc,"%f",&fitness)==1) {
                        VLACheck(align->id_at_point,int,point_id);
                        VLACheck(align->fitness,float,point_id);
                          if(point_id >= align->n_point) 
                            align->n_point = point_id + 1;
                          align->id_at_point[point_id] = atom_id;
                          align->fitness[point_id] = fitness;
                          /*                        printf("read alignment atom %d to target %s point %d fitness %8.3f\n",
                                                    atom_id,target,point_id,fitness);*/
                        }
                      }
                    }
                  }
                }
              }
              parsed_flag = true;
            }
            break;
        }
          
        if((!parsed_flag)&&m4x->context) { /* named context of some sort... */
          int cn;
          int found=false;
            
          /* does context already exist ? */
          for(cn=0;cn<m4x->n_context;cn++) {
            if(WordMatchExact(G,m4x->context[cn].name,cc,true)) {
              found=true;
              break;
            }
          }
            
          if(found) {
              
            M4XContextType *cont = m4x->context+cn;
              
            p = nskip(p,1);
            p = ntrim(cc,p,6);
            switch(cc[0]) {
            case 'B':
              if(WordMatchExact(G,"BORDER",cc,true)&&bondFlag) {
                int order;
                  
                p=nskip(p,1);
                p=ncopy(cc,p,6);
                if(sscanf(cc,"%d",&b1)==1) {
                  p=nskip(p,1);
                  p=ncopy(cc,p,6);
                  if(sscanf(cc,"%d",&b2)==1) {
                    p = nskip(p,1);
                    p = ncopy(cc,p,6);
                    if(sscanf(cc,"%d",&order)==1) {      
                      if((b1>=0)&&(b2>=0)) { /* IDs must be positive */
                        VLACheck(bond,BondType,nBond);
                        if(b1<=b2) {
                          bond[nBond].index[0]=b1; /* temporarily store the atom indexes */
                          bond[nBond].index[1]=b2;
                          bond[nBond].order=order;
                          bond[nBond].stereo=0;
                            
                        } else {
                          bond[nBond].index[0]=b2;
                          bond[nBond].index[1]=b1;
                          bond[nBond].order=order;
                          bond[nBond].stereo=0;
                        }
                        nBond++;
                      }
                    }
                  }
                }
              }
              break;
            case 'S':
              if(WordMatchExact(G,"SITE",cc,true)) {
                if(cont->site) {
                  int id;
                  while(*cc) {
                    p = nskip(p,1);
                    p = ncopy(cc,p,6);
                    if(sscanf(cc,"%d",&id)==1) {
                      VLACheck(cont->site,int,cont->n_site);
                      cont->site[cont->n_site++] = id;
                    }
                  }
                }
              } 
              break;
            case 'L':
              if(WordMatchExact(G,"LIGAND",cc,true)) {
                if(cont->ligand) {
                  int id;
                  while(*cc) {
                    p = nskip(p,1);
                    p = ncopy(cc,p,6);
                    if(sscanf(cc,"%d",&id)==1) {
                      VLACheck(cont->ligand,int,cont->n_ligand);
                      cont->ligand[cont->n_ligand++] = id;
                    }
                  }
                }
              }
              break;
            case 'W':
              if(WordMatchExact(G,"WATER",cc,true)) {
                if(cont->water) {
                  int id;
                  while(*cc) {
                    p = nskip(p,1);
                    p = ncopy(cc,p,6);
                    if(sscanf(cc,"%d",&id)==1) {
                      VLACheck(cont->water,int,cont->n_water);
                      cont->water[cont->n_water++] = id;
                    }
                  }
                }
              } 
              break;
            case 'H':
              if(WordMatchExact(G,"HBOND",cc,true)) {
                if(cont->hbond) {
                  int id1,id2;
                  float strength;
                  p = nskip(p,1);
                  p = ncopy(cc,p,6);
                  if(sscanf(cc,"%d",&id1)==1) {
                    p = nskip(p,1);
                    p = ncopy(cc,p,6);
                    if(sscanf(cc,"%d",&id2)==1) {
                      p = nskip(p,1);
                      p = ncopy(cc,p,6);
                      if(sscanf(cc,"%f",&strength)==1) {                  
                        VLACheck(cont->hbond,M4XBondType,cont->n_hbond);
                        cont->hbond[cont->n_hbond].atom1 = id1;
                        cont->hbond[cont->n_hbond].atom2 = id2;
                        cont->hbond[cont->n_hbond].strength = strength;
                        cont->n_hbond++;
                      }
                    }
                  }
                }
              }
              break;
            case 'N':
              if(WordMatchExact(G,"NBOND",cc,true)) {
                if(cont->nbond) {
                  int id1,id2;
                  float strength;
                  p = nskip(p,1);
                  p = ncopy(cc,p,6);
                  if(sscanf(cc,"%d",&id1)==1) {
                    p = nskip(p,1);
                    p = ncopy(cc,p,6);
                    if(sscanf(cc,"%d",&id2)==1) {
                      p = nskip(p,1);
                      p = ncopy(cc,p,6);
                      if(sscanf(cc,"%f",&strength)==1) {                  
                        VLACheck(cont->nbond,M4XBondType,cont->n_nbond);
                        cont->nbond[cont->n_nbond].atom1 = id1;
                        cont->nbond[cont->n_nbond].atom2 = id2;
                        cont->nbond[cont->n_nbond].strength = strength;
                        cont->n_nbond++;
                      }
                    }
                  }
                }
              }
              break;
            }
          }
        }
      }
    } else if((p[0]=='A')&& /* ANISOU */
              (p[1]=='N')&&
              (p[2]=='I')&&
              (p[3]=='S')&&
              (p[4]=='O')&&
              (p[5]=='U')&&
              (!*restart_model)&&
              (atomCount)) {
      ai=atInfo + atomCount - 1;
      
      /* TODO: check atom identifier match */

      {
        int dummy; 
        p=nskip(p,6);
        p=ncopy(cc,p,5);
        if(!sscanf(cc,"%d",&dummy)) dummy = 0;
        if(dummy == ai->id) { /* ATOM ID must match */
          p=nskip(p,17);
          {
            int dummy;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%d",&dummy))
              ai->U11 = dummy/10000.0F;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%d",&dummy)) 
              ai->U22 = dummy/10000.0F;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%d",&dummy)) 
              ai->U33 = dummy/10000.0F;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%d",&dummy)) 
              ai->U12 = dummy/10000.0F;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%d",&dummy)) 
              ai->U13 = dummy/10000.0F;
            p=ncopy(cc,p,7);
            if(sscanf(cc,"%d",&dummy)) 
              ai->U23 = dummy/10000.0F;
          }
        }
      }
    }


    /* END KEYWORDS */
    
    /* Secondary structure records */
    
    
    if(SSCode) {
        
      /* pretty confusing how this works... the following efficient (i.e. array-based)
         secondary structure lookup even when there are multiple insertion codes
         and where there may be multiple SS records for the residue using different 
         insertion codes */

      if(!ss[ss_chain1]) { /* allocate new array for chain if necc. */
        ss[ss_chain1]=Calloc(int,cResvMask+1);
      }

      sst = NULL; 
      for(b=ss_resv1;b<=ss_resv2;b++) { /* iterate over all residues indicated */
        index = b & cResvMask;
        if(ss[ss_chain1][index]) sst = NULL; /* make a unique copy in the event of multiple entries for one resv */
        if(!sst) {
          VLACheck(ss_list,SSEntry,n_ss);
          ssi = n_ss++;
          sst = ss_list + ssi;
          sst->resv1 = ss_resv1;
          sst->resv2 = ss_resv2;
          sst->chain1 = ss_chain1;
          sst->chain2 = ss_chain2;
          sst->type=SSCode;
          strcpy(sst->resi1,ss_resi1);
          strcpy(sst->resi2,ss_resi2);
          ss_found=true;
        }
        sst->next = ss[ss_chain1][index];
        ss[ss_chain1][index]=ssi;
        if(sst->next) sst = NULL; /* force another unique copy */
      }
        
      if(ss_chain2!=ss_chain1) { /* handle case where chains are not the same (undefined in PDB spec?) */
        if(!ss[ss_chain2]) {
          ss[ss_chain2]=Calloc(int,cResvMask+1);
        }
        sst = NULL; 
        for(b=ss_resv1;b<=ss_resv2;b++) { /* iterate over all residues indicated */
          index = b & cResvMask;
          if(ss[ss_chain2][index]) sst = NULL; /* make a unique copy in the event of multiple entries for one resv */
          if(!sst) {
            VLACheck(ss_list,SSEntry,n_ss);
            ssi = n_ss++;
            sst = ss_list + ssi;
            sst->resv1 = ss_resv1;
            sst->resv2 = ss_resv2;
            sst->chain1 = ss_chain1;
            sst->chain2 = ss_chain2;
            sst->type=SSCode;
            strcpy(sst->resi1,ss_resi1);
            strcpy(sst->resi2,ss_resi2);
          }
          sst->next = ss[ss_chain2][index];
          ss[ss_chain2][index]=ssi;
          if(sst->next) sst = NULL; /* force another unique copy */
        }
      }
    }
    /* Atom records */

    if(AFlag&&(!*restart_model))
      {
        ai=atInfo+atomCount;
        p=nskip(p,6);
        p=ncopy(cc,p,5);
        if(!sscanf(cc,"%d",&ai->id)) ai->id=0;
        ai->rank = atomCount;

        p=nskip(p,1); /* to 12 */
        p=ncopy(literal_name,p,4); 
        if(literal_names) {
          strcpy(ai->name,literal_name);
        } else {
          ParseNTrim(ai->name,literal_name,4); 
        }

        p=ncopy(cc,p,1);
        if(*cc==32)
          ai->alt[0]=0;
        else {
          ai->alt[0]=*cc;
          ai->alt[1]=0;
        }

        p=ncopy(cc,p,4);  /* now allowing for 4-letter residues */
        if(!sscanf(cc,"%s",ai->resn)) 
          ai->resn[0]=0;
        else if(truncate_resn) /* unless specifically disabled */
          ai->resn[3]=0;

        if(ai->name[0]) {
          int name_len = strlen(ai->name);
          char name[4];
          switch(reformat_names) {
          case 1: /* pdb compliant: HH12 becomes 2HH1, etc. */
            if(name_len>3) {
              if((ai->name[0]>='A')&&((ai->name[0]<='Z'))&&
                 (ai->name[3]>='0')&&(ai->name[3]<='9')) {
                if(!(((ai->name[1]>='a')&&(ai->name[1]<='z'))||
                     ((ai->name[0]=='C')&&(ai->name[1]=='L'))|| /* try to be smart about */
                     ((ai->name[0]=='B')&&(ai->name[1]=='R'))|| /* distinguishing common atoms */
                     ((ai->name[0]=='C')&&(ai->name[1]=='A'))|| /* in all-caps from typical */
                     ((ai->name[0]=='F')&&(ai->name[1]=='E'))|| /* nonatomic abbreviations */
                     ((ai->name[0]=='C')&&(ai->name[1]=='U'))||
                     ((ai->name[0]=='N')&&(ai->name[1]=='A'))||
                     ((ai->name[0]=='N')&&(ai->name[1]=='I'))||
                     ((ai->name[0]=='M')&&(ai->name[1]=='G'))||
                     ((ai->name[0]=='M')&&(ai->name[1]=='N'))||
                     ((ai->name[0]=='H')&&(ai->name[1]=='G'))||
                     ((ai->name[0]=='S')&&(ai->name[1]=='E'))||
                     ((ai->name[0]=='S')&&(ai->name[1]=='I'))||
                     ((ai->name[0]=='Z')&&(ai->name[1]=='N'))
                     )) {
                  ctmp = ai->name[3];
                  ai->name[3]= ai->name[2];
                  ai->name[2]= ai->name[1];
                  ai->name[1]= ai->name[0];
                  ai->name[0]= ctmp;
                }
              }
            } else if(name_len==3) {
              if((ai->name[0]=='H')&&
                 (ai->name[1]>='A')&&((ai->name[1]<='Z'))&&
                 (ai->name[2]>='0')&&(ai->name[2]<='9')) {
                AtomInfoGetPDB3LetHydroName(G,ai->resn,ai->name,name);
                if(name[0]==' ')
                  strcpy(ai->name,name+1);
                else
                  strcpy(ai->name,name);
              }
            }
            break;
          case 2: /* amber compliant: 2HH1 becomes HH12 */
          case 3: /* pdb compliant, but use IUPAC within PyMOL */
            if(ai->name[0]) {
              if((ai->name[0]>='0')&&(ai->name[0]<='9')&&
                 (!((ai->name[1]>='0')&&(ai->name[1]<='9')))&&
                 (ai->name[1]!=0)) {
                switch(strlen(ai->name)) {
                default:
                  break;
                case 2:
                  ctmp = ai->name[0];
                  ai->name[0]= ai->name[1];
                  ai->name[1]= ctmp;
                  break;
                case 3:
                  ctmp = ai->name[0];
                  ai->name[0]= ai->name[1];
                  ai->name[1]= ai->name[2];
                  ai->name[2]= ctmp;
                  break;
                case 4:
                  ctmp = ai->name[0];
                  ai->name[0]= ai->name[1];
                  ai->name[1]= ai->name[2];
                  ai->name[2]= ai->name[3];
                  ai->name[3]= ctmp;
                  break;
                }
                break;
              default: /* AS IS */
                break;
              }
            }
            break;
          case 4: /* simply read trim and write back out with 3-letter names starting from the
                     second column, and four-letter names starting in the first */
            ncopy(cc,ai->name,44);
            ParseNTrim(ai->name,cc,4);             
            break;
          }
        }

        p=ncopy(cc,p,1);
        if(*cc==' ')
          ai->chain[0]=0;
        else {
          ai->chain[0] = *cc;
          ai->chain[1] = 0;
        }

        p=ncopy(cc,p,5); /* we treat insertion records as part of the residue identifier */
        if(!sscanf(cc,"%s",ai->resi)) ai->resi[0]=0;
        ai->resv = AtomResvFromResi(ai->resi);
          
        if(ssFlag) { /* get secondary structure information (if avail) */

          ss_chain1=ai->chain[0];
          index = ai->resv & cResvMask;
          if(ss[ss_chain1]) {
            ssi = ss[ss_chain1][index];
            while(ssi) {
              sst = ss_list + ssi; /* contains shared entry, or unique linked list for each residue */
              /*                printf("%d<=%d<=%d, %s<=%s<=%s ", 
                                sst->resv1,ai->resv,sst->resv2,
                                sst->resi1,ai->resi,sst->resi2);*/
              if(ai->resv>=sst->resv1)
                if(ai->resv<=sst->resv2)
                  if((ai->resv!=sst->resv1)||(WordCompare(G,ai->resi,sst->resi1,true)>=0))
                    if((ai->resv!=sst->resv2)||(WordCompare(G,ai->resi,sst->resi2,true)<=0))
                      {
                        ai->ssType[0]=sst->type;
                        /*                          printf(" Y\n");*/
                        break;
                      }
              /*                printf(" N\n");*/
              ssi = sst->next;
            }
          }
            
        } else {
          ai->cartoon = cCartoon_tube;
        }
        if((!info)||(!info->is_pqr_file)) { /* standard PDB file */

          p=nskip(p,3);
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+a);
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+(a+1));
          p=ncopy(cc,p,8);
          sscanf(cc,"%f",coord+(a+2));
            
          p=ncopy(cc,p,6);
          if(!sscanf(cc,"%f",&ai->q))
            ai->q=1.0;
            
          p=ncopy(cc,p,6);
          if(!sscanf(cc,"%f",&ai->b))
            ai->b=0.0;
            
          p=nskip(p,6);
          p=ncopy(cc,p,4);
          if(!ignore_pdb_segi) {
            if(!segi_override[0])
              {
                if(!sscanf(cc,"%s",ai->segi)) 
                  ai->segi[0]=0;
                else {
                  cc_saved=cc[3];
                  ncopy(cc,p,4); 
                  if((cc_saved=='1')&& /* atom ID overflow? (nonstandard use...)...*/
                     (cc[0]=='0')&& 
                     (cc[1]=='0')&&
                     (cc[2]=='0')&&
                     (cc[3]=='0')&&
                     atomCount) {
                    strcpy(segi_override,(ai-1)->segi);
                    strcpy(ai->segi,(ai-1)->segi);
                  }
                }
              } else {
              strcpy(ai->segi,segi_override);
            }
          } else {              
            ai->segi[0]=0;
          }
          
          p=ncopy(cc,p,2);
          if(!sscanf(cc,"%s",ai->elem)) 
            ai->elem[0]=0;          
          else if(!((((ai->elem[0]>='a')&&(ai->elem[0]<='z'))|| /* don't get confused by PDB misuse */
                     ((ai->elem[0]>='A')&&(ai->elem[0]<='Z')))&&
                    (((ai->elem[1]==0)||
                      ((ai->elem[1]>='a')&&(ai->elem[1]<='z'))|| 
                      ((ai->elem[1]>='A')&&(ai->elem[1]<='Z'))))))
            ai->elem[0]=0;                      
            
          if(!ai->elem[0]) {
            if(((literal_name[0]==' ')||((literal_name[0]>='0')&&(literal_name[0]<='9')))&&
               (literal_name[1]>='A')&&(literal_name[1]<='Z')) { /* infer element from name column */
              ai->elem[0]=literal_name[1];
              ai->elem[1]=0;
            } else if(((literal_name[0]>='A')&&(literal_name[0]<='Z'))&&
                      (((literal_name[1]>='A')&&(literal_name[1]<='Z'))||
                       ((literal_name[1]>='a')&&(literal_name[1]<='z')))) { /* infer element from name column */
                ai->elem[0]=literal_name[0];
                ai->elem[2]=0;
                if((literal_name[1]>='A')&&(literal_name[1]<='Z')) { /* second letter is capitalized */
                  if(bogus_name_alignment) { 
                    /* if other atom names aren't properly aligned */
                    ai->elem[1]=0; /* kill 2nd letter */
                  } else if(literal_name[0]=='H') { 
                    /* or if this is an ultra-bogus PDB with inconsistent 
                       indendentation, and this is likely a hydrogen */
                    ai->elem[1]=0; /* kill 2nd letter */
                  } else {
                    ai->elem[1]=tolower(literal_name[1]);
                  }
                } else
                  ai->elem[1]=literal_name[1];
            }
          }
          
          p=ncopy(cc,p,2);
          if((cc[1]=='-')||(cc[1]=='+')) {
            /* only read formal charge when sign is present */
            char ctmp = cc[0];
            cc[0] = cc[1];
            cc[1] = ctmp;
            if(!sscanf(cc,"%d",&ai->formalCharge))
              ai->formalCharge=0;
          }

          for(c=0;c<cRepCnt;c++) {
            ai->visRep[c] = false;
          }
          
          /* end normal PDB */
        } else if(info&&info->is_pqr_file) {
          /* PQR file format...not well defined, but basically PDB
             with charge and radius instead of B and Q.  Right now,
             we insist on PDB column format through the chain ID,
             and then switch over to whitespace delimited parsing 
             for the coordinates, charge, and radius */

          p=ParseWordNumberCopy(cc,p,MAXLINELEN-1); 
          sscanf(cc,"%f",coord+a);
          p=ParseWordNumberCopy(cc,p,MAXLINELEN-1);
          sscanf(cc,"%f",coord+(a+1));
          p=ParseWordNumberCopy(cc,p,MAXLINELEN-1);
          sscanf(cc,"%f",coord+(a+2));

          p=ParseWordNumberCopy(cc,p,MAXLINELEN-1);
          if(!sscanf(cc,"%f",&ai->partialCharge))
            ai->partialCharge=0.0F;

          p=ParseWordNumberCopy(cc,p,MAXLINELEN-1);            
          if(sscanf(cc,"%f",&ai->elec_radius)!=1)
            ai->elec_radius = 0.0F;
        }

        ai->visRep[cRepLine] = auto_show_lines; /* show lines by default */
        ai->visRep[cRepNonbonded] = auto_show_nonbonded; /* show lines by default */
        ai->visRep[cRepSphere] = auto_show_spheres;

        if(AFlag==1) 
          ai->hetatm=0;
        else {
          ai->hetatm=1;
          ai->flags=cAtomFlag_ignore;
        }
          
        AtomInfoAssignParameters(G,ai);
        AtomInfoAssignColors(G,ai);

        PRINTFD(G,FB_ObjectMolecule)
          "%s %s %s %s %8.3f %8.3f %8.3f %6.2f %6.2f %s\n",
          ai->name,ai->resn,ai->resi,ai->chain,
          *(coord+a),*(coord+a+1),*(coord+a+2),ai->b,ai->q,
          ai->segi
          ENDFD;

        if(atomCount<(nAtom-1)) { /* safety */
          a+=3;
          atomCount++;
        }
      }
    p=nextline(p);
  }
  
  /* END PASS 2 */
  
  if(bondFlag) {
    UtilSortInPlace(G,bond,nBond,sizeof(BondType),(UtilOrderFn*)BondInOrder);              
    if(nBond) {
      if(!have_bond_order) { /* handle PDB bond-order kludge */
        ii1=bond;
        ii2=bond+1;
        nReal=1;
        ii1->order=1;
        a=nBond-1;
        while(a) {
          if((ii1->index[0]==ii2->index[0])&&(ii1->index[1]==ii2->index[1])) {
            ii1->order++; /* count dup */
          } else {
            ii1++; /* non-dup, make copy */
            ii1->index[0]=ii2->index[0];
            ii1->index[1]=ii2->index[1];
            ii1->order=ii2->order;
            ii1->stereo=ii2->stereo;
            nReal++;
          }
          ii2++;
          a--;
        }
        nBond=nReal;
      }
      /* now, find atoms we're looking for */

      /* determine ranges */
      maxAt=nAtom;
      ii1=bond;
      for(a=0;a<nBond;a++) {
        if(ii1->index[0]>maxAt) maxAt=ii1->index[0];
        if(ii1->index[1]>maxAt) maxAt=ii1->index[1];
        ii1++;
      }
      for(a=0;a<nAtom;a++) 
        if(maxAt<atInfo[a].id) maxAt=atInfo[a].id;
      /* build index */
      maxAt++;
      idx = Alloc(int,maxAt+1);
      for(a=0;a<maxAt;a++) {
        idx[a]=-1;
      }
      for(a=0;a<nAtom;a++)
        idx[atInfo[a].id]=a;
      /* convert indices to bonds */
      ii1=bond;
      ii2=bond;
      nReal=0;
      {
        int unbond_cations = SettingGetGlobal_i(G,cSetting_pdb_unbond_cations);
        int flag;

        for(a=0;a<nBond;a++) {
          
          if((ii1->index[0]>=0)&&((ii1->index[1])>=0)) {
            if((idx[ii1->index[0]]>=0)&&(idx[ii1->index[1]]>=0)) { /* in case PDB file has bad bonds */
              ii2->index[0]=idx[ii1->index[0]];
              ii2->index[1]=idx[ii1->index[1]];
              ii2->order=ii1->order;
              if((ii2->index[0]>=0)&&(ii2->index[1]>=0)) {
                
                if(!have_bond_order) { /* handle PDB bond order kludge */
                  if(ii1->order<=2) ii2->order=1;
                  else if(ii1->order<=4) ii2->order=2;
                  else if(ii1->order<=6) ii2->order=3;
                  else ii2->order=4;
                }
                flag=true;
                if(unbond_cations) {
                  if(AtomInfoIsFreeCation(G,atInfo + ii2->index[0]))
                    flag=false;
                  else if(AtomInfoIsFreeCation(G,atInfo + ii2->index[1]))
                    flag=false;
                }
                if(flag) {
                  atInfo[ii2->index[0]].bonded=true;
                  atInfo[ii2->index[1]].bonded=true;
                  nReal++;
                  ii2++;
                }
              }
            }
          }
          ii1++;
        }
      }
      nBond=nReal;
      FreeP(idx);
    }
  }
  if(ss_found&&!quiet) {
    PRINTFB(G,FB_ObjectMolecule,FB_Details)
      " ObjectMolecule: Read secondary structure assignments.\n"
      ENDFB(G);
  }
  if(symmetry&&!quiet&&(!only_read_one_model)) {
    PRINTFB(G,FB_ObjectMolecule,FB_Details)
      " ObjectMolecule: Read crystal symmetry information.\n"
      ENDFB(G);
  }
  PRINTFB(G,FB_ObjectMolecule,FB_Blather)
   " PDBStr2CoordSet: Read %d bonds from CONECT records (%p).\n",nBond,
    (void*)bond
    ENDFB(G);
  cset = CoordSetNew(G);
  cset->NIndex=nAtom;
  cset->Coord=coord;
  cset->TmpBond=bond;
  cset->NTmpBond=nBond;
  if(symmetry) cset->Symmetry=symmetry;
  if(atInfoPtr)
	 *atInfoPtr = atInfo;

  if((*restart_model)&&(!foundNextModelFlag)) {
    /* if plan on need to reading another model into this object, 
       make sure there is another model to read...*/
    p=*restart_model;
    while(*p) {
      if((p[0]== 'H')&&(p[1]=='E')&&(p[2]=='A')&&
        (p[3]=='D')&&(p[4]=='E')&&(p[5]=='R')) {
        /* seeing HEADER means we're off the end of the existing file */
        break;
      }
      if((p[0]== 'M')&&(p[1]=='O')&&(p[2]=='D')&&
         (p[3]=='E')&&(p[4]=='L')&&(p[5]==' ')) {
        foundNextModelFlag=true;
        break;
      }
      if((p[0]== 'E')&&(p[1]=='N')&&(p[2]=='D')&&
         (p[3]=='M')&&(p[4]=='D')&&(p[5]=='L')) {
        foundNextModelFlag=true;
        break;
      }
      p=nextline(p);
    }
    if(!foundNextModelFlag) {
      *restart_model=NULL;
    }
  }

  if(ssFlag) {
    for(a=0;a<=255;a++) {
      FreeP(ss[a]);
    }
    VLAFreeP(ss_list);
  }

  if(!seen_model)
    *model_number = 1;

  if((*restart_model) && (*next_pdb)) { /* if we're mixing multistate objects and
                                           trajectories, then enforce sanity by 
                                           reading the models first... */
    if(is_end_of_object)
      (*restart_model) = NULL;
    else if((*next_pdb)<(*restart_model))
      (*next_pdb) = NULL;
  }
  return(cset);
}
/*========================================================================*/

void ObjectMoleculeM4XAnnotate(ObjectMolecule *I,M4XAnnoType *m4x,char *script_file,
                               int match_colors,int nbr_sele)
{
  int a;
  WordType name;
  M4XContextType *cont;

  if(m4x) {
    for(a=0;a<m4x->n_context;a++) {
      cont = m4x->context + a;

      if(cont->site) {
        UtilNCopy(name,I->Obj.Name,sizeof(WordType));
        UtilNConcat(name,"_",sizeof(WordType));
        UtilNConcat(name,cont->name,sizeof(WordType));
        UtilNConcat(name,"_site",sizeof(WordType));
        SelectorSelectByID(I->Obj.G,name,I,cont->site,cont->n_site);
      }
      if(cont->ligand) {
        UtilNCopy(name,I->Obj.Name,sizeof(WordType));
        UtilNConcat(name,"_",sizeof(WordType));
        UtilNConcat(name,cont->name,sizeof(WordType));
        UtilNConcat(name,"_ligand",sizeof(WordType));
        SelectorSelectByID(I->Obj.G,name,I,cont->ligand,cont->n_ligand);
      }
      if(cont->water) {
        UtilNCopy(name,I->Obj.Name,sizeof(WordType));
        UtilNConcat(name,"_",sizeof(WordType));
        UtilNConcat(name,cont->name,sizeof(WordType));
        UtilNConcat(name,"_water",sizeof(WordType));
        SelectorSelectByID(I->Obj.G,name,I,cont->water,cont->n_water);
      }
      if(cont->hbond) {
        ObjectDist *distObj;
        UtilNCopy(name,I->Obj.Name,sizeof(WordType));
        UtilNConcat(name,"_",sizeof(WordType));
        UtilNConcat(name,cont->name,sizeof(WordType));
        UtilNConcat(name,"_hbond",sizeof(WordType));
        ExecutiveDelete(I->Obj.G,name);
        distObj = ObjectDistNewFromM4XBond(I->Obj.G,NULL,
                                            I,
                                            cont->hbond,
                                           cont->n_hbond,
                                           nbr_sele);
        if(match_colors)
          distObj->Obj.Color = I->Obj.Color;
        else
          distObj->Obj.Color = ColorGetIndex(I->Obj.G,"yellow");
        ObjectSetName((CObject*)distObj,name);
        if(distObj)
          ExecutiveManageObject(I->Obj.G,(CObject*)distObj,false,true);
      }

      if(cont->nbond&&0) {
        /*        ObjectDist *distObj;*/
        UtilNCopy(name,I->Obj.Name,sizeof(WordType));
        UtilNConcat(name,"_",sizeof(WordType));
        UtilNConcat(name,cont->name,sizeof(WordType));
        UtilNConcat(name,"_nbond",sizeof(WordType));
        ExecutiveDelete(I->Obj.G,name);
        /*        distObj = ObjectDistNewFromM4XBond(I->Obj.G,NULL,
                                            I,
                                             cont->nbond,
                                             cont->n_nbond);
         if(distObj)
        ExecutiveManageObject(I->Obj.G,(CObject*)distObj,false,true); */

        {
          CGO *cgo = NULL;
          ObjectCGO *ocgo;

          
          cgo=CGONew(I->Obj.G);
          /*
            CGOBegin(cgo,GL_LINES);
            for(a=0;a<op1.nvv1;a++) {
            CGOVertexv(cgo,op2.vv1+(a*3));
            MatrixApplyTTTfn3f(1,v1,op2.ttt,op1.vv1+(a*3));
            CGOVertexv(cgo,v1);
          */
          CGOEnd(cgo);
          CGOStop(cgo);
          ocgo = ObjectCGOFromCGO(I->Obj.G,NULL,cgo,0);
          if(match_colors)
            ocgo->Obj.Color = I->Obj.Color;
          else
            ocgo->Obj.Color = ColorGetIndex(I->Obj.G,"yellow");
          ObjectSetName((CObject*)ocgo,name);
          ExecutiveDelete(I->Obj.G,name);

          ExecutiveManageObject(I->Obj.G,(CObject*)ocgo,false,true);

          SceneInvalidate(I->Obj.G);
        }

      }

    }
    if(script_file) 
      PParse(I->Obj.G,script_file);
  }
}

void ObjectMoleculeInitHBondCriteria(PyMOLGlobals *G,HBondCriteria *hbc)
{
  hbc->maxAngle = SettingGet_f(G,NULL,NULL,cSetting_h_bond_max_angle);
  hbc->maxDistAtMaxAngle = SettingGet_f(G,NULL,NULL,cSetting_h_bond_cutoff_edge);
  hbc->maxDistAtZero = SettingGet_f(G,NULL,NULL,cSetting_h_bond_cutoff_center);
  hbc->power_a = SettingGet_f(G,NULL,NULL,cSetting_h_bond_power_a);
  hbc->power_b = SettingGet_f(G,NULL,NULL,cSetting_h_bond_power_b);
  hbc->cone_dangle = (float)cos(PI*0.5*SettingGet_f(G,NULL,NULL,cSetting_h_bond_cone)/180.0F);
  if(hbc->maxDistAtMaxAngle!=0.0F) {
    hbc->factor_a = (float)(0.5/pow(hbc->maxAngle,hbc->power_a));
    hbc->factor_b = (float)(0.5/pow(hbc->maxAngle,hbc->power_b));
  }
}
    
static int ObjectMoleculeTestHBond(float *donToAcc,float *donToH, float *hToAcc, 
                          float *accPlane, HBondCriteria *hbc)
{
  float nDonToAcc[3],nDonToH[3],nAccPlane[3],nHToAcc[3];
  double angle;
  double cutoff;
  double curve;
  double dist;
  float dangle;
/* A ~~ H - D */

  normalize23f(donToAcc,nDonToAcc);
  normalize23f(hToAcc,nHToAcc);
  if(accPlane) { /* remember, plane need not exist if it's water... */
    normalize23f(accPlane,nAccPlane);
#if 0
    if(dot_product3f(nDonToAcc,nAccPlane)>(-hbc->cone_dangle)) /* don't allow D behind Acceptor plan */
      return 0;
#endif
    if(dot_product3f(nHToAcc,nAccPlane)>(-hbc->cone_dangle)) /* don't allow H behind Acceptor plane */
      return 0;
  }

  normalize23f(donToH,nDonToH);
  normalize23f(donToAcc,nDonToAcc);

         
  dangle = dot_product3f(nDonToH,nDonToAcc);
  if((dangle<1.0F)&&(dangle>0.0F))
    angle = 180.0 * acos((double)dangle) / PI;
  else if(dangle>0.0F)
    angle = 0.0;
  else
    angle = 90.0;

  if(angle > hbc->maxAngle)
    return 0;

  /* interpolate cutoff based on ADH angle */

  if(hbc->maxDistAtMaxAngle!=0.0F) {
    curve = (pow(angle, hbc->power_a) * hbc->factor_a + 
             pow(angle, hbc->power_b) * hbc->factor_b );
    
    cutoff = (hbc->maxDistAtMaxAngle * curve) + 
      (hbc->maxDistAtZero * (1.0-curve));
  } else {
    cutoff = hbc->maxDistAtZero;
  }

  /*
  printf("angle %8.3f curve %8.3f %8.3f %8.3f %8.3f\n",angle,
         curve,cutoff,hbc->maxDistAtMaxAngle,hbc->maxDistAtZero);
  */

  dist = length3f(donToAcc);  

  if(dist>cutoff) 
    return 0;
  else
    return 1;

}

/*========================================================================*/

static int ObjectMoleculeFindBestDonorH(ObjectMolecule *I,
                                        int atom,
                                        int state,
                                        float *dir,
                                        float *best,
                                        int *is_real)
{
  int result = 0;
  CoordSet *cs;
  int n,nn;
  int idx;
  int a1;
  float cand[3],cand_dir[3];
  float best_dot=0.0F,cand_dot;
  float *orig;

  ObjectMoleculeUpdateNeighbors(I);  

  if((state>=0)&&
     (state<I->NCSet)&&
     (cs=I->CSet[state])&&
     (atom<I->NAtom)) {
    
    if(I->DiscreteFlag) {
      if(cs==I->DiscreteCSet[atom]) {
        idx=I->DiscreteAtmToIdx[atom];
      } else {
        idx=-1;
      }
    } else {
      idx=cs->AtmToIdx[atom];
    }
    
    if(idx>=0) {

      orig = cs->Coord + 3*idx;
      
      /*  do we need to add any new hydrogens? */
      
      n = I->Neighbor[atom];
      nn = I->Neighbor[n++];

      /*      printf("nn %d valence %d %s\n",nn,
             I->AtomInfo[atom].valence,I->AtomInfo[atom].name);
      {
        int i;
        for(i=0;i<nn;i++) {
          printf("%d \n",I->Neighbor[n+2*i]);
        }
      }
      */

      if((nn<I->AtomInfo[atom].valence)||I->AtomInfo[atom].hb_donor) {       /* is there an implicit hydrogen? */
        if(ObjectMoleculeFindOpenValenceVector(I,state,atom,best,dir,-1)) {
          result = true;
          best_dot = dot_product3f(best,dir);
          add3f(orig,best,best);
          if(is_real) *is_real = false;
        }
      }
      /* iterate through real hydrogens looking for best match
         with desired direction */
      
      while(1) { /* look for an attached non-hydrogen as a base */
        a1 = I->Neighbor[n];
        n+=2; 
        if(a1<0) break;
        if(I->AtomInfo[a1].protons==1) { /* hydrogen */
          if(ObjectMoleculeGetAtomVertex(I,state,a1,cand)) { /* present */
            
            subtract3f(cand,orig,cand_dir);
            normalize3f(cand_dir);
            cand_dot = dot_product3f(cand_dir,dir);
            if(result) { /* improved */
              if((best_dot<cand_dot)||((is_real)&&(!*is_real))) {
                best_dot = cand_dot;
                copy3f(cand,best);
                if(is_real)  *is_real = true;
              }
            } else { /* first */
              result = true;
              copy3f(cand,best);
              best_dot = cand_dot;
              if(is_real)  *is_real = true;
            }
          }
        }
      }
    }
  }
  return result;
}

/*========================================================================*/

int ObjectMoleculeGetCheckHBond(int *h_is_real,
                                float *h_crd_ret,
                                ObjectMolecule *don_obj,
                                int don_atom,
                                int don_state,
                                ObjectMolecule *acc_obj,
                                int acc_atom,
                                int acc_state,
                                HBondCriteria *hbc)
{
  int result = 0;
  CoordSet *csD, *csA;
  int idxD,idxA;
  float *vAcc,*vDon;
  float donToAcc[3];
  float donToH[3];
  float bestH[3];
  float hToAcc[3];
  float accPlane[3],*vAccPlane = NULL;

  /* first, check for existence of coordinate sets */
  
  if((don_state>=0)&&
     (don_state<don_obj->NCSet)&&
     (csD=don_obj->CSet[don_state])&&
     (acc_state>=0)&&
     (acc_state<acc_obj->NCSet)&&
     (csA=acc_obj->CSet[acc_state])&&
     (don_atom<don_obj->NAtom)&&
     (acc_atom<acc_obj->NAtom)) {

    /* now check for coordinates of these actual atoms */
    
    if(don_obj->DiscreteFlag) {
      if(csD==don_obj->DiscreteCSet[don_atom]) {
        idxD=don_obj->DiscreteAtmToIdx[don_atom];
      } else {
        idxD=-1;
      }
    } else {
      idxD=csD->AtmToIdx[don_atom];
    }
    
    if(acc_obj->DiscreteFlag) {
      if(csA==acc_obj->DiscreteCSet[acc_atom]) {
        idxA=acc_obj->DiscreteAtmToIdx[acc_atom];
      } else {
        idxA=-1;
      }
    } else {
      idxA=csA->AtmToIdx[acc_atom];
    }
    
    if((idxA>=0)&&(idxD>=0)) {
      

      /* now get local geometries, including 
         real or virtual hydrogen atom positions */
      
      vDon = csD->Coord + 3*idxD;
      vAcc = csA->Coord + 3*idxA;
      
      subtract3f(vAcc,vDon,donToAcc);

      if(ObjectMoleculeFindBestDonorH(don_obj,
                                      don_atom,
                                      don_state,
                                      donToAcc,
                                      bestH,
                                      h_is_real))
        {

          subtract3f(bestH,vDon,donToH);
          subtract3f(vAcc,bestH,hToAcc);
          
          if(ObjectMoleculeGetAvgHBondVector(acc_obj,acc_atom,
                                             acc_state,accPlane,
                                             hToAcc)>0.1) {
            vAccPlane = &accPlane[0];
          }

#if 0
          {
            float tmp[3];
            CGOLinewidth(DebugCGO,4.0F);
            CGOBegin(DebugCGO,GL_LINES);

            CGOColor(DebugCGO,0.0F,1.0F,0.0F); /* green */
            CGOVertexv(DebugCGO,vDon);
            normalize23f(donToAcc,tmp);
            add3f(vDon,tmp,tmp);
            CGOVertexv(DebugCGO,tmp);                        

            CGOColor(DebugCGO,1.0F,0.0F,0.0F); /* red */
            CGOVertexv(DebugCGO,bestH);
            normalize23f(hToAcc,tmp);
            add3f(bestH,tmp,tmp);
            CGOVertexv(DebugCGO,tmp);            

            CGOColor(DebugCGO,0.0F,1.0F,1.0F); /* cyan */
            CGOVertexv(DebugCGO,vDon);
            normalize23f(donToH,tmp);
            add3f(vDon,tmp,tmp);
            CGOVertexv(DebugCGO,tmp);                        

            if(vAccPlane) {
              CGOColor(DebugCGO,1.0F,1.0F,0.0F); /* yellow */
              CGOVertexv(DebugCGO,vAcc);
              normalize23f(vAccPlane,tmp);
              add3f(vAcc,tmp,tmp);
              CGOVertexv(DebugCGO,tmp);                        
            }

            CGOEnd(DebugCGO);

          }
#endif
          result = ObjectMoleculeTestHBond(donToAcc,donToH,hToAcc, 
                                           vAccPlane,hbc);
          if(result && h_crd_ret && h_is_real && *h_is_real)
            copy3f(bestH,h_crd_ret);
        }
    }
  }
  
  return(result);
}
                                
/*========================================================================*/

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
#ifndef _PYMOL_NOPY
static PyObject *ObjectMoleculeCSetAsPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NCSet);
  for(a=0;a<I->NCSet;a++) {
    if(I->CSet[a]) {
      PyList_SetItem(result,a,CoordSetAsPyList(I->CSet[a]));
    } else {
      PyList_SetItem(result,a,PConvAutoNone(Py_None));
    }
  }
  return(PConvAutoNone(result));
}
#endif

/*static PyObject *ObjectMoleculeDiscreteCSetAsPyList(ObjectMolecule *I)
  {
  PyObject *result = NULL;
  return(PConvAutoNone(result));
  }*/


#ifndef _PYMOL_NOPY
static int ObjectMoleculeCSetFromPyList(ObjectMolecule *I,PyObject *list)
{
  int ok=true;
  int a;
  if(ok) ok=PyList_Check(list);
  if(ok) {
    VLACheck(I->CSet,CoordSet*,I->NCSet);
    for(a=0;a<I->NCSet;a++) {
      if(ok) ok = CoordSetFromPyList(I->Obj.G,PyList_GetItem(list,a),&I->CSet[a]);
      PRINTFB(I->Obj.G,FB_ObjectMolecule,FB_Debugging)
	     " ObjectMoleculeCSetFromPyList: ok %d after CoordSet %d\n",ok,a
	  ENDFB(I->Obj.G);

      if(ok) 
        if(I->CSet[a]) /* WLD 030205 */
          I->CSet[a]->Obj = I;
    }
  }
  return(ok);
}
#endif

#ifndef _PYMOL_NOPY
static PyObject *ObjectMoleculeBondAsPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  PyObject *bond_list;
  BondType *bond;
  int a;

  result = PyList_New(I->NBond);  
  bond = I->Bond;
  for(a=0;a<I->NBond;a++) {
    bond_list=PyList_New(7);
    PyList_SetItem(bond_list,0,PyInt_FromLong(bond->index[0]));
    PyList_SetItem(bond_list,1,PyInt_FromLong(bond->index[1]));
    PyList_SetItem(bond_list,2,PyInt_FromLong(bond->order));
    PyList_SetItem(bond_list,3,PyInt_FromLong(bond->id));
    PyList_SetItem(bond_list,4,PyInt_FromLong(bond->stereo));
    PyList_SetItem(bond_list,5,PyInt_FromLong(bond->unique_id));
    PyList_SetItem(bond_list,6,PyInt_FromLong(bond->has_setting));
    PyList_SetItem(result,a,bond_list);
    bond++;
  }

  return(PConvAutoNone(result));
}
#endif

#ifndef _PYMOL_NOPY
static int ObjectMoleculeBondFromPyList(ObjectMolecule *I,PyObject *list) 
{
  int ok=true;
  int a;
  int stereo,ll = 0;
  PyObject *bond_list=NULL;
  BondType *bond;
  if(ok) ok=PyList_Check(list);  
  if(ok) ok=((I->Bond=VLAlloc(BondType,I->NBond))!=NULL);
  bond = I->Bond;
  for(a=0;a<I->NBond;a++) {
    if(ok) bond_list = PyList_GetItem(list,a);
    if(ok) ok = PyList_Check(bond_list);
    if(ok) ll = PyList_Size(bond_list);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,0),&bond->index[0]);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,1),&bond->index[1]);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,2),&bond->order);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,3),&bond->id);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,4),&stereo);
    if(ok) bond->stereo=(short int)stereo;
    if(ok&&(ll>5)) {
      int has_setting;
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,5),&bond->unique_id);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,6),&has_setting);
      if(ok) bond->has_setting = (short int)has_setting;
      if(ok && bond->unique_id) { /* reserve existing IDs */
        bond->unique_id = SettingUniqueConvertOldSessionID(I->Obj.G,bond->unique_id);
        AtomInfoReserveUniqueID(I->Obj.G,bond->unique_id);
      }
    }
    bond++;
  }
  PRINTFB(I->Obj.G,FB_ObjectMolecule,FB_Debugging)
    " ObjectMoleculeBondFromPyList: ok %d after restore\n",ok
    ENDFB(I->Obj.G);
  
  return(ok);
}
#endif

#ifndef _PYMOL_NOPY
static PyObject *ObjectMoleculeAtomAsPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  AtomInfoType *ai;
  int a;

  result = PyList_New(I->NAtom);  
  ai = I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    PyList_SetItem(result,a,AtomInfoAsPyList(I->Obj.G,ai));
    ai++;
  }
  return(PConvAutoNone(result));
}
#endif

#ifndef _PYMOL_NOPY
static int ObjectMoleculeAtomFromPyList(ObjectMolecule *I,PyObject *list) 
{
  int ok=true;
  int a;
  AtomInfoType *ai;
  if(ok) ok=PyList_Check(list);  
  VLACheck(I->AtomInfo,AtomInfoType,I->NAtom+1);
  ai = I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(ok) ok = AtomInfoFromPyList(I->Obj.G,ai,PyList_GetItem(list,a));
    ai++;
  }
      PRINTFB(I->Obj.G,FB_ObjectMolecule,FB_Debugging)
	     " ObjectMoleculeAtomFromPyList: ok %d \n",ok
	  ENDFB(I->Obj.G);
  return(ok);
}
#endif

int ObjectMoleculeNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectMolecule **result)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
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

  I=ObjectMoleculeNew(G,discrete_flag);
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(G,PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NCSet);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->NBond);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,3),&I->NAtom);
  if(ok) ok = ObjectMoleculeCSetFromPyList(I,PyList_GetItem(list,4));
  if(ok) ok = CoordSetFromPyList(G,PyList_GetItem(list,5),&I->CSTmpl);
  if(ok) ok = ObjectMoleculeBondFromPyList(I,PyList_GetItem(list,6));
  if(ok) ok = ObjectMoleculeAtomFromPyList(I,PyList_GetItem(list,7));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,8),&I->DiscreteFlag);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,9),&I->NDiscrete);
  if(ok) I->Symmetry = SymmetryNewFromPyList(G,PyList_GetItem(list,10));
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
  
  ObjectMoleculeInvalidate(I,cRepAll,cRepInvAll,-1);
  if(ok) 
    (*result) = I;
  else {
    /* cleanup? */
  }

  return(ok);
#endif
}


/*========================================================================*/
PyObject *ObjectMoleculeAsPyList(ObjectMolecule *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
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
#endif
}


/*========================================================================*/
int ObjectMoleculeConnect(ObjectMolecule *I,BondType **bond,AtomInfoType *ai,
                          struct CoordSet *cs,int bondSearchMode,
                          int connectModeOverride)
{
  #define cMULT 1
  PyMOLGlobals *G=I->Obj.G;
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
  int repeat = true;
  int discrete_chains = SettingGetGlobal_i(G,cSetting_pdb_discrete_chains);
  int connect_bonded = SettingGetGlobal_b(G,cSetting_connect_bonded);
  int connect_mode = SettingGetGlobal_i(G,cSetting_connect_mode);
  int unbond_cations = SettingGetGlobal_i(G,cSetting_pdb_unbond_cations);
  cutoff_v=SettingGet(G,cSetting_connect_cutoff);
  cutoff_s=cutoff_v + 0.2F;
  cutoff_h=cutoff_v - 0.2F;
  max_cutoff = cutoff_s;

  if(connectModeOverride>=0) 
    connect_mode = connectModeOverride;
    
  if(connect_mode==2) { /* force use of distance-based connectivity,
                         ignoring that provided with file */
    bondSearchMode = true;
    cs->NTmpBond = 0;
    FreeP(cs->TmpBond);
  }

  /*  FeedbackMask[FB_ObjectMolecule]=0xFF;*/
  nBond = 0;
  maxBond = cs->NIndex * 8;
  (*bond) = VLACalloc(BondType,maxBond);
  while(repeat) {
    repeat = false;

    if(cs->NIndex&&bondSearchMode) { /* &&(!I->DiscreteFlag) WLD 010527 */
      
      PRINTFB(G,FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeConnect: Searching for bonds amongst %d coordinates.\n",cs->NIndex
        ENDFB(G);
      if(Feedback(G,FB_ObjectMolecule,FB_Debugging)) {
        for(a=0;a<cs->NIndex;a++)
          printf(" ObjectMoleculeConnect: coord %d %8.3f %8.3f %8.3f\n",
                 a,cs->Coord[a*3],cs->Coord[a*3+1],cs->Coord[a*3+2]);
      }
      
      switch(connect_mode) {
      case 0: /* distance-based and explicit (not HETATM to HETATM) */
      case 3: /* distance-based and explicit (even HETATM to HETATM) */
      case 2: /* distance-based only */ {
        /* distance-based bond location  */
        int violations = 0;
        int *cnt = Alloc(int,cs->NIndex);
        int valcnt;

        for(i=0;i<cs->NIndex;i++) {
          valcnt = AtomInfoGetExpectedValence(G,ai+cs->IdxToAtm[i]);
          if(valcnt<0)
            valcnt=6;
          cnt[i]=valcnt;
        }
        
        map=MapNew(G,max_cutoff+MAX_VDW,cs->Coord,cs->NIndex,NULL);
        if(map) {
          for(i=0;i<cs->NIndex;i++) {
            if(nBond>maxBond)
              break; 
            v1=cs->Coord+(3*i);
            MapLocus(map,v1,&a,&b,&c);
            for(d=a-1;d<=a+1;d++)
              for(e=b-1;e<=b+1;e++)
                for(f=c-1;f<=c+1;f++) {
                  
                  j = *(MapFirst(map,d,e,f));
                  while(j>=0) {
                    
                    if(i<j) {
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
                      if(AtomInfoKnownWaterResName(G,ai1->resn))
                        water_flag=true;
                      else if(AtomInfoKnownWaterResName(G,ai2->resn))
                        water_flag=true;
                      
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
                          (water_flag||(!cs->TmpBond)||
                           ((!(ai1->hetatm&&ai2->hetatm) || 
                             (connect_mode == 3)))) &&
                          ((discrete_chains<1) ||
                           ai1->chain[0]==ai2->chain[0]) &&
                          (connect_bonded || (!(ai1->bonded&&ai2->bonded)))) { 
                        flag=true;
                        if(water_flag)
                          if(!AtomInfoSameResidue(G,ai1,ai2))
                            flag=false;
                        
                        if(flag) {
                          if(ai1->alt[0]!=ai2->alt[0]) { /* handle alternate conformers */
                            if(ai1->alt[0]&&ai2->alt[0])
                              flag=false; /* don't connect atoms with different, non-NULL
                                             alternate conformations */
                          } else if(ai1->alt[0]&&ai2->alt[0])
                            if(!AtomInfoSameResidue(G,ai1,ai2))
                              if(ai1->alt[0]!=ai2->alt[0])
                                flag=false; /* don't connect different, non-NULL 
                                               alt conformations in 
                                               different residues */
                        }
                        
                        if(flag) {
                          if(ai1->alt[0]||ai2->alt[0]) 
                            if(water_flag) /* hack to clean up water bonds */
                              if(!AtomInfoSameResidue(G,ai1,ai2))
                                flag=false;
                        }
                        
                        if(flag && unbond_cations) {
                          if(AtomInfoIsFreeCation(G,ai1))
                            flag=false;
                          else if(AtomInfoIsFreeCation(G,ai2))
                            flag=false;
                        }

                        if(flag) {
                          VLACheck((*bond),BondType,nBond);
                          (*bond)[nBond].index[0] = a1;
                          (*bond)[nBond].index[1] = a2;
                          (*bond)[nBond].stereo = 0;
                          order = 1;
                          if(discrete_chains<0) { /* if we're allowing bonds between chains,
                                                     then make sure things don't get out of hand */
                            if(cnt[i]==-1)
                              violations++;
                            if(cnt[j]==-1)
                              violations++;
                            cnt[i]--;
                            cnt[j]--;
                            if(violations>(cs->NIndex>>3)) { 
                              /* if more than 12% of the structure has excessive #'s of bonds... */
                              PRINTFB(G,FB_ObjectMolecule,FB_Blather)
                                " ObjectMoleculeConnect: Assuming chains are discrete...\n"
                                ENDFB(G);
                              discrete_chains = 1;
                              repeat = true;
                              goto do_it_again;
                            }
                          }
                          if(ai1->hetatm&&(!ai1->resn[3])) { /* common HETATMs we should know about... */
                            switch(ai1->resn[0]) {
                            case 'M':
                              switch(ai1->resn[1]) {
                              case 'S':
                                switch(ai1->resn[2]) {
                                case 'E':
                                  if(((!ai1->name[1])&&(!ai2->name[1]))&&
                                     (((ai1->name[0]=='C')&&(ai2->name[0]=='O'))||
                                      ((ai1->name[0]=='O')&&(ai2->name[0]=='C')))) {
                                    if(AtomInfoSameResidue(G,ai1,ai2)) {
                                      order = 2;
                                    }
                                  }
                                  break;
                                }
                                break;
                              }
                              break;
                            }
                          } else if((!ai1->hetatm)) { 
                            if(AtomInfoSameResidue(I->Obj.G,ai1,ai2)) {
                              /* Standard disconnected PDB residue */
                              assign_pdb_known_residue(G,ai1,ai2,&order);
                            }
                          }
                          (*bond)[nBond].order = -order; /* store tentative valence as negative */
                          nBond++;
                        }
                      }
                    }
                    j=MapNext(map,j);
                  }
                }
          }
        do_it_again:
          MapFree(map);
          FreeP(cnt);
        }
      }
      case 1: /* only use explicit connectivity from file (don't do anything) */ 
        break;
      case 4:  /* dictionary-based connectivity */
        /* TODO */
        break;
      }
      PRINTFB(G,FB_ObjectMolecule,FB_Blather)
        " ObjectMoleculeConnect: Found %d bonds.\n",nBond
        ENDFB(G);
      if(Feedback(G,FB_ObjectMolecule,FB_Debugging)) {
        for(a=0;a<nBond;a++)
          printf(" ObjectMoleculeConnect: bond %d ind0 %d ind1 %d\n",
                 a,(*bond)[a].index[0],(*bond)[a].index[1]);
      }
    }
    if(repeat) {
      nBond = 0;
    }
  }
  if(cs->NTmpBond&&cs->TmpBond) {
    int check_conect_all = false;
    int pdb_conect_all = false;
    PRINTFB(G,FB_ObjectMolecule,FB_Blather) 
      " ObjectMoleculeConnect: incorporating explicit bonds. %d %d\n",
      nBond,cs->NTmpBond
      ENDFB(G);
    if((nBond==0) && (cs->NTmpBond>0) &&
       bondSearchMode && (connect_mode == 0) && cs->NIndex) {
      /* if we were no bonds were found, and we have explicit connectivity,
       * try to determine if we need to set pdb_conect_mode */
      for(i=0;i<cs->NIndex;i++) {
        a1=cs->IdxToAtm[i];
        ai1=ai+a1;
        if(ai1->bonded && (!ai1->hetatm)) { 
          /* apparent PDB ATOM record with explicit bonding... */
          check_conect_all = true;
          break;
        }
      }
    }

    VLACheck((*bond),BondType,(nBond+cs->NTmpBond));
    ii1=(*bond)+nBond;
    ii2=cs->TmpBond;
    {
      register int n_atom = I->NAtom;
      for(a=0;a<cs->NTmpBond;a++) {
        a1 = cs->IdxToAtm[ii2->index[0]]; /* convert bonds from index space */
        a2 = cs->IdxToAtm[ii2->index[1]]; /* to atom space */
        if((a1>=0)&&(a2>=0)&&(a1<n_atom)&&(a2<n_atom)) {
          if(check_conect_all) { 
            if((!ai[a1].hetatm)&&(!ai[a2].hetatm)) { 
              /* found bond between non-HETATMs -- so tell PyMOL to CONECT all ATOMs
               * when writing out a PDB file */
              pdb_conect_all = true;
            }
          }
          ai[a1].bonded=true;
          ai[a2].bonded=true;
          (*ii1) = (*ii2); /* note this copies owned ids and thus settings etc. */
          ii1->index[0]=a1;
          ii1->index[1]=a2;
          ii1++;
          ii2++;
          
        }
      }
    }
      
    nBond=nBond+cs->NTmpBond;
    VLAFreeP(cs->TmpBond);
    cs->NTmpBond=0;

    if(pdb_conect_all) { 
      int dummy;
      if(!SettingGetIfDefined_b(G,I->Obj.Setting,cSetting_pdb_conect_all,&dummy)) {
        CSetting **handle = NULL;
        if(I->Obj.fGetSettingHandle) {
          handle = I->Obj.fGetSettingHandle(&I->Obj,-1);
          if(handle) {
            SettingCheckHandle(G,handle);        
            SettingSet_b(*handle,cSetting_pdb_conect_all,true);
          }
        }
      }
    }
  }


  if(cs->NTmpLinkBond&&cs->TmpLinkBond) {
    PRINTFB(G,FB_ObjectMolecule,FB_Blather) 
      "ObjectMoleculeConnect: incorporating linkage bonds. %d %d\n",
      nBond,cs->NTmpLinkBond
      ENDFB(G);
    VLACheck((*bond),BondType,(nBond+cs->NTmpLinkBond));
    ii1=(*bond)+nBond;
    ii2=cs->TmpLinkBond;
    for(a=0;a<cs->NTmpLinkBond;a++)
      {
        a1 = ii2->index[0]; /* first atom is in object */
        a2 = cs->IdxToAtm[ii2->index[1]]; /* second is in the cset */
        ai[a1].bonded=true;
        ai[a2].bonded=true;
        (*ii1) = (*ii2); /* note this copies owned ids and thus settings etc. */
        ii1->index[0]=a1;
        ii1->index[1]=a2;
        ii1++;
        ii2++;
      }
    nBond=nBond+cs->NTmpLinkBond;
    VLAFreeP(cs->TmpLinkBond);
    cs->NTmpLinkBond=0;
  }

  PRINTFD(G,FB_ObjectMolecule)
    " ObjectMoleculeConnect: elminating duplicates with %d bonds...\n",nBond
    ENDFD;

  if(!I->DiscreteFlag) {
    UtilSortInPlace(G,(*bond),nBond,sizeof(BondType),(UtilOrderFn*)BondInOrder);
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
            if((ii2->order>0)&&((ii1-1)->order<0))
              (ii1-1)->order = ii2->order; /* use most certain valence */
            ii2++; /* skip bond */
          }
        }
      VLASize((*bond),BondType,nBond);
    }
  }
  /* restore bond oder positivity */

  ii1 = *bond;
  for(a=0;a<nBond;a++) {
    if(ii1->order<0)
      ii1->order = -ii1->order;
    ii1++;
  }

  PRINTFD(G,FB_ObjectMolecule)
    " ObjectMoleculeConnect: leaving with %d bonds...\n",nBond
    ENDFD;
  return(nBond);
}


/*========================================================================*/
void ObjectMoleculeSort(ObjectMolecule *I) /* sorts atoms and bonds */
{
  int *index;
  int *outdex=NULL;
  register int a,b;
  CoordSet *cs,**dcs;
  AtomInfoType *atInfo;
  int *dAtmToIdx;
  if(!I->DiscreteFlag) { /* currently, discrete objects are never sorted */
    int n_bytes = sizeof(int)*I->NAtom;
    int already_in_order = true;
    int i_NAtom = I->NAtom;
    index=AtomInfoGetSortedIndex(I->Obj.G,&I->Obj,I->AtomInfo,i_NAtom,&outdex);
    for(a=0;a<i_NAtom;a++) {
      if(index[a]!=a) {
        already_in_order = false;
        break;
      }
    }
    if(!already_in_order) { /* if we aren't already in perfect order */

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
          register int cs_NIndex = cs->NIndex;
          register int *cs_IdxToAtm = cs->IdxToAtm;
          register int *cs_AtmToIdx = cs->AtmToIdx;
          for(b=0;b<cs_NIndex;b++)
            cs_IdxToAtm[b]=outdex[cs_IdxToAtm[b]];
          if(cs_AtmToIdx) {
            memset(cs_AtmToIdx, -1, n_bytes);
            /*          for(b=0;b<i_NAtom;b++)
                        cs_AtmToIdx[b]=-1;*/
            for(b=0;b<cs_NIndex;b++) {
              cs_AtmToIdx[cs_IdxToAtm[b]]=b;
            }
          }
        }
      }
      
      atInfo=(AtomInfoType*)VLAMalloc(i_NAtom,sizeof(AtomInfoType),5,true);
      /* autozero here is important */
      for(a=0;a<i_NAtom;a++)
        atInfo[a]=I->AtomInfo[index[a]];
      VLAFreeP(I->AtomInfo);
      I->AtomInfo=atInfo;
      
      if(I->DiscreteFlag) {
        dcs = VLAlloc(CoordSet*,i_NAtom);
        dAtmToIdx = VLAlloc(int,i_NAtom);
        for(a=0;a<i_NAtom;a++) {
          b=index[a];
          dcs[a] = I->DiscreteCSet[b];
          dAtmToIdx[a] = I->DiscreteAtmToIdx[b];
        }
        VLAFreeP(I->DiscreteCSet);
        VLAFreeP(I->DiscreteAtmToIdx);
        I->DiscreteCSet = dcs;
        I->DiscreteAtmToIdx = dAtmToIdx;
      }
    }
    AtomInfoFreeSortedIndexes(I->Obj.G,index,outdex);
    
    UtilSortInPlace(I->Obj.G,I->Bond,I->NBond,sizeof(BondType),(UtilOrderFn*)BondInOrder);
    /* sort...important! */
    ObjectMoleculeInvalidate(I,cRepAll,cRepInvAtoms,-1); /* important */
    
  }
}
