/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2004 by Warren Lyford Delano of DeLano Scientific. 
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
#include"Map.h"
#include"Vector.h"
#include"Debug.h"
#include"Err.h"
#include"Word.h"
#include"Util.h"
#include"PConv.h"
#include"P.h"

#include"MemoryDebug.h"
#include"Selector.h"
#include"Executive.h"
#include"ObjectMolecule.h"
#include"CoordSet.h"
#include"DistSet.h"
#include"Word.h"
#include"Scene.h"
#include"CGO.h"
#include"Seq.h"
#include"Editor.h"
#include"Seeker.h"

#include"OVContext.h"
#include"OVLexicon.h"
#include"OVOneToAny.h"
#include"Parse.h"

#define SelectorWordLength 1024
typedef char SelectorWordType[SelectorWordLength];

#define cSelectorTmpPrefix "_sel_tmp_"
#define cSelectorTmpPattern "_sel_tmp_*"

#define cNDummyModels 2
#define cNDummyAtoms 2

#define cDummyOrigin 0
#define cDummyCenter 1

/* special selections, unknown to executive */
#define cColorectionFormat "_!c_%s_%d" 

static WordKeyValue rep_names[] = { 
  { "spheres", cRepSphere },
  { "sticks", cRepCyl },
  { "surface", cRepSurface },
  { "labels", cRepLabel },
  { "nb_spheres", cRepNonbondedSphere },
  { "cartoon", cRepCartoon },
  { "ribbon", cRepRibbon },
  { "lines", cRepLine },
  { "dots", cRepDot },
  { "mesh", cRepMesh },
  { "nonbonded", cRepNonbonded },
  { "ellipsoid", cRepEllipsoid },
  { "", },
};



typedef struct {
  int level, imp_op_level;
  int type; /* 0 = value 1 = operation 2 = pre-operation */
  unsigned int code; 
  SelectorWordType text;
  int *sele;
} EvalElem;

typedef struct {
  int model;
  int atom;
  int index;
  float f1;
} TableRec;

typedef struct {
  int ID;
  int justOneObjectFlag;
  ObjectMolecule *theOneObject;
  int justOneAtomFlag;
  int theOneAtom;
} SelectionInfoRec;

struct _CSelector {
  MemberType *Member; /* Must be first in structure, so that we can get this w/o knowing the struct */
  SelectorWordType *Name; /* this seems rather excessive, since name len < ObjNameMax ... */
  SelectionInfoRec *Info;
  int NSelection,NActive;
  int TmpCounter;
  int NMember;
  int FreeMember;
  ObjectMolecule **Obj;
  TableRec *Table;
  float *Vertex;
  int *Flag1,*Flag2;
  ov_size NAtom;
  ov_size NModel;
  int NCSet;
  int SeleBaseOffsetsValid;
  ObjectMolecule *Origin,*Center;
  OVLexicon *Lex;
  OVOneToAny *Key;
  OVOneToOne *NameOffset;

};

typedef struct {
  int depth1;
  int depth2;
  int depth3;
  int depth4;
  int sum;
  int frag;
} WalkDepthRec;

static int *SelectorSelect(PyMOLGlobals *G,char *sele,int state,int domain);
static int SelectorGetInterstateVLA(PyMOLGlobals *G,int sele1,int state1,int sele2,int state2,
									  float cutoff,int **vla);

static int SelectorGetArrayNCSet(PyMOLGlobals *G,int *array,int no_dummies);

static int SelectorModulate1(PyMOLGlobals *G,EvalElem *base,int state);
static int SelectorSelect0(PyMOLGlobals *G,EvalElem *base);
static int SelectorSelect1(PyMOLGlobals *G,EvalElem *base);
static int SelectorSelect2(PyMOLGlobals *G,EvalElem *base);
static int SelectorLogic1(PyMOLGlobals *G,EvalElem *base, int state);
static int SelectorLogic2(PyMOLGlobals *G,EvalElem *base);
static int SelectorOperator22(PyMOLGlobals *G,EvalElem *base,int state);
static int *SelectorEvaluate(PyMOLGlobals *G,SelectorWordType *word,int state);
static SelectorWordType *SelectorParse(PyMOLGlobals *G,char *s);
static void SelectorPurgeMembers(PyMOLGlobals *G,int sele);
static int  SelectorEmbedSelection(PyMOLGlobals *G,int *atom, char *name, 
                                   ObjectMolecule *obj,int no_dummies,
                                   int exec_manage);

static void SelectorClean(PyMOLGlobals *G);
static int *SelectorGetIndexVLA(PyMOLGlobals *G,int sele);
static int *SelectorApplyMultipick(PyMOLGlobals *G,Multipick *mp);
static int SelectorCheckNeighbors(PyMOLGlobals *G,int maxDepth,ObjectMolecule *obj,int at1,int at2,
                           int *zero,int *scratch);

static int *SelectorUpdateTableSingleObject(PyMOLGlobals *G,ObjectMolecule *obj,
                                            int req_state,
                                            int no_dummies,int *idx,
                                            int n_idx,int numbered_tags);

static int SelectorGetObjAtmOffset(CSelector *I,ObjectMolecule *obj,int offset)
{
  if(I->SeleBaseOffsetsValid) {
    return obj->SeleBase + offset;
  } else {
    register ov_diff stop_below = obj->SeleBase;
    register ov_diff stop_above = I->NAtom - 1;
    register int result = stop_below;
    register TableRec *i_table = I->Table;
    register ObjectMolecule **i_obj=I->Obj;
    register int step = offset;
    register int cur;
    register int proposed;
    register int prior1=-1, prior2=-1;

    /* non-linear hunt to find atom */

    result = stop_below;
    cur = i_table[result].atom;
    while(step>1) {
      if(cur < offset) {
        stop_below = result + 1;
        while(step>1) {
          proposed = result + step;
          if(proposed <= stop_above) {
            if(i_obj[i_table[proposed].model] == obj) {
              if(proposed == prior1) {
                proposed--;
                step--; /* guarantee progress (avoid flip flop) */
              }
              result = prior1 = proposed;
              break;
            } else if(stop_above>proposed) {
              stop_above = proposed - 1;
            }
          }
          step = (step>>1);
        }
      } else if(cur > offset) {
        stop_above = result - 1;
        while(step>1) {
          proposed = result - step;
          if(proposed >= stop_below) {
            if(i_obj[i_table[proposed].model] == obj) {
              if(proposed == prior2) {
                proposed++;
                step--; /* guarantee progress (avoid flip flop) */
              }
              result = prior2 = proposed;
              break;
            }
          }
          step = (step>>1);
        }
      } else 
        return result;
      cur = i_table[result].atom;
      if(cur == offset) 
        return result;
    }

    {
      /* failsafe / linear search */
      register int dir = 1;
      if(cur > offset)
        dir = -1;
      while(1) { /* TODO: optimize this search algorithm! */
        if(cur == offset)
          return result;
        if(dir>0) {
          if(result >= stop_above)
            break;
          result++;
        } else {
          if(result <= stop_below)
            break;
          result--;
        }
        if(i_obj[i_table[result].model] != obj)
          break;
        cur = i_table[result].atom;      
      }
    }
  }
  return -1;
}


#define STYP_VALU 0
#define STYP_OPR1 1
#define STYP_OPR2 2
#define STYP_SEL0 3
#define STYP_SEL1 4
#define STYP_SEL2 5
#define STYP_LIST 6
#define STYP_PRP1 7
#define STYP_SEL3 8
#define STYP_PVAL 0
#define STYP_OP22 9   /* sele oper arg1 arg2 sele */

/*                  code   |   type    | priority */

#define SELE_NOT1 ( 0x0100 | STYP_OPR1 | 0x70 )
#define SELE_BYR1 ( 0x0200 | STYP_OPR1 | 0x20 )
#define SELE_AND2 ( 0x0300 | STYP_OPR2 | 0x60 )
#define SELE_OR_2 ( 0x0400 | STYP_OPR2 | 0x40 )
#define SELE_IN_2 ( 0x0500 | STYP_OPR2 | 0x40 )
#define SELE_ALLz ( 0x0600 | STYP_SEL0 | 0x90 )
#define SELE_NONz ( 0x0700 | STYP_SEL0 | 0x90 )
#define SELE_HETz ( 0x0800 | STYP_SEL0 | 0x80 )
#define SELE_HYDz ( 0x0900 | STYP_SEL0 | 0x90 )
#define SELE_VISz ( 0x0A00 | STYP_SEL0 | 0x90 )
#define SELE_ARD_ ( 0x0B00 | STYP_PRP1 | 0x30 )
#define SELE_EXP_ ( 0x0C00 | STYP_PRP1 | 0x30 )
#define SELE_NAMs ( 0x0D00 | STYP_SEL1 | 0x80 )
#define SELE_ELEs ( 0x0E00 | STYP_SEL1 | 0x80 )
#define SELE_RSIs ( 0x0F00 | STYP_SEL1 | 0x80 )
#define SELE_CHNs ( 0x1000 | STYP_SEL1 | 0x80 )
#define SELE_SEGs ( 0x1100 | STYP_SEL1 | 0x80 )
#define SELE_MODs ( 0x1200 | STYP_SEL1 | 0x80 )
#define SELE_IDXs ( 0x1300 | STYP_SEL1 | 0x80 )
#define SELE_RSNs ( 0x1400 | STYP_SEL1 | 0x80 )
#define SELE_SELs ( 0x1500 | STYP_SEL1 | 0x80 )
#define SELE_BVLx ( 0x1600 | STYP_SEL2 | 0x80 )
#define SELE_ALTs ( 0x1700 | STYP_SEL1 | 0x80 )
#define SELE_FLGs ( 0x1800 | STYP_SEL1 | 0x80 )
#define SELE_GAP_ ( 0x1900 | STYP_PRP1 | 0x80 )
#define SELE_TTYs ( 0x1A00 | STYP_SEL1 | 0x80 )  
#define SELE_NTYs ( 0x1B00 | STYP_SEL1 | 0x80 )
#define SELE_PCHx ( 0x1C00 | STYP_SEL2 | 0x80 )  
#define SELE_FCHx ( 0x1D00 | STYP_SEL2 | 0x80 )
#define SELE_ID_s ( 0x1E00 | STYP_SEL1 | 0x80 )
#define SELE_BNDz ( 0x1F00 | STYP_SEL0 | 0x80 )
#define SELE_LIK2 ( 0x2000 | STYP_OPR2 | 0x40 )
#define SELE_NGH1 ( 0x2100 | STYP_OPR1 | 0x20 )
#define SELE_QVLx ( 0x2200 | STYP_SEL2 | 0x80 )
#define SELE_BYO1 ( 0x2300 | STYP_OPR1 | 0x20 )
#define SELE_SSTs ( 0x2400 | STYP_SEL1 | 0x80 )
#define SELE_STAs ( 0x2500 | STYP_SEL1 | 0x80 )
#define SELE_PREz ( 0x2500 | STYP_SEL0 | 0x80 )
#define SELE_WIT_ ( 0x2600 | STYP_OP22 | 0x30 ) 
#define SELE_ORIz ( 0x2700 | STYP_SEL0 | 0x90 )
#define SELE_CENz ( 0x2800 | STYP_SEL0 | 0x90 )
#define SELE_ENAz ( 0x2900 | STYP_SEL0 | 0x90 )
#define SELE_REPs ( 0x2A00 | STYP_SEL1 | 0x80 )
#define SELE_COLs ( 0x2B00 | STYP_SEL1 | 0x80 )
#define SELE_HBDs ( 0x2C00 | STYP_SEL0 | 0x80 )
#define SELE_HBAs ( 0x2D00 | STYP_SEL0 | 0x80 )
#define SELE_BYC1 ( 0x2E00 | STYP_OPR1 | 0x20 )
#define SELE_BYS1 ( 0x2F00 | STYP_OPR1 | 0x20 )
#define SELE_BYM1 ( 0x3000 | STYP_OPR1 | 0x20 )
#define SELE_BYF1 ( 0x3100 | STYP_OPR1 | 0x20 )
#define SELE_EXT_ ( 0x3200 | STYP_PRP1 | 0x30 )
#define SELE_BON1 ( 0x3300 | STYP_OPR1 | 0x50 )
#define SELE_FST1 ( 0x3400 | STYP_OPR1 | 0x30 )
#define SELE_CAS1 ( 0x3500 | STYP_OPR1 | 0x30 )
#define SELE_BEY_ ( 0x3600 | STYP_OP22 | 0x30 ) 
#define SELE_POLz ( 0x3700 | STYP_SEL0 | 0x90 )
#define SELE_SOLz ( 0x3800 | STYP_SEL0 | 0x90 )
#define SELE_ORGz ( 0x3900 | STYP_SEL0 | 0x90 )
#define SELE_INOz ( 0x3A00 | STYP_SEL0 | 0x90 )
#define SELE_GIDz ( 0x3B00 | STYP_SEL0 | 0x90 )
#define SELE_RNKs ( 0x3C00 | STYP_SEL1 | 0x80 )
#define SELE_PEPs ( 0x3D00 | STYP_SEL1 | 0x80 )
#define SELE_ACCz ( 0x3E00 | STYP_SEL0 | 0x90 )
#define SELE_DONz ( 0x3F00 | STYP_SEL0 | 0x90 )
#define SELE_LST1 ( 0x4000 | STYP_OPR1 | 0x30 )
#define SELE_NTO_ ( 0x4100 | STYP_OP22 | 0x30 ) 
#define SELE_CCLs ( 0x4200 | STYP_SEL1 | 0x80 )
#define SELE_RCLs ( 0x4300 | STYP_SEL1 | 0x80 )
#define SELE_PTDz ( 0x4400 | STYP_SEL0 | 0x90 )
#define SELE_MSKz ( 0x4500 | STYP_SEL0 | 0x90 )
#define SELE_IOR2 ( 0x4600 | STYP_OPR2 | 0x10 )
#define SELE_FXDz ( 0x4700 | STYP_SEL0 | 0x90 )
#define SELE_RSTz ( 0x4800 | STYP_SEL0 | 0x90 )
#define SELE_ANT2 ( 0x4900 | STYP_OPR2 | 0x60 )
#define SELE_BYX1 ( 0x4A00 | STYP_OPR1 | 0x20 )

#define SEL_PREMAX 0x8

static WordKeyValue Keyword[] = 
{
  {  "not",      SELE_NOT1 },
  {  "!",        SELE_NOT1 },

  {  "neighbor", SELE_NGH1 },
  {  "nbr;",     SELE_NGH1 }, /* deprecated */
  {  "nbr.",     SELE_NGH1 },

  {  "byfragment",SELE_BYF1 },
  {  "byfrag"   ,SELE_BYF1 },
  {  "bf."      ,SELE_BYF1 },

  {  "byresidue",SELE_BYR1 },
  {  "byresi",   SELE_BYR1 }, /* unofficial */
  {  "byres",    SELE_BYR1 },
  {  "br;",      SELE_BYR1 },/* deprecated */
  {  "br.",      SELE_BYR1 },
  {  "b;",       SELE_BYR1 }, /* deprecated */

  {  "bychain",  SELE_BYC1 },
  {  "bc.",      SELE_BYC1 },

  {  "byobject", SELE_BYO1 },
  {  "byobj",    SELE_BYO1 },
  {  "bo;",      SELE_BYO1 },/* deprecated */
  {  "bo.",      SELE_BYO1 },

  {  "bound_to", SELE_BON1 },
  {  "bto.",     SELE_BON1 },

  { "bymolecule",SELE_BYM1 },
  {  "bymol",    SELE_BYM1 },
  {  "bm.",      SELE_BYM1 },

  {  "bysegment",SELE_BYS1 },
  {  "byseg",    SELE_BYS1 }, 
  {  "bysegi",   SELE_BYS1 }, /* unofficial */
  {  "bs.",      SELE_BYS1 },

  {  "bycalpha", SELE_CAS1 },
  {  "bca.",     SELE_CAS1 },

  {  "first",    SELE_FST1 },
  {  "last",     SELE_LST1 },

  {  "and",      SELE_AND2 },
  {  "&",        SELE_AND2 },
  {  "or",       SELE_OR_2 },
  {  "+",        SELE_OR_2 }, /* added to mitigate damage caused by the obj1+obj2 parser bug */
  {  "-",        SELE_ANT2 }, /* added to provide natural complement to the above: an AND NOT or SUBTRACT operation*/ 
  {  "|",        SELE_OR_2 },
  {  "in",       SELE_IN_2 },

  {  "like",     SELE_LIK2 },
  {  "l;",       SELE_LIK2 },
  {  "l.",       SELE_LIK2 },

  {  cKeywordAll,      SELE_ALLz }, /* 0 parameter */

  {  cKeywordNone,     SELE_NONz }, /* 0 parameter */
  {  "hetatm",   SELE_HETz }, /* 0 parameter */
  {  "het",      SELE_HETz }, /* 0 parameter */

  {  "hydrogens",SELE_HYDz }, /* 0 parameter */
  {  "hydro",    SELE_HYDz }, /* 0 parameter */
  {  "h;",       SELE_HYDz }, /* deprecated */
  {  "h.",       SELE_HYDz }, /* 0 parameter */

  {  "hba.",     SELE_HBAs },
  {  "hbd.",     SELE_HBDs },

  {  "visible",  SELE_VISz }, /* 0 parameter */
  {  "v;",       SELE_VISz }, /* 0 parameter */
  {  "v.",       SELE_VISz }, /* 0 parameter */

  {  "around",   SELE_ARD_ }, /* 1 parameter */
  {  "a;",       SELE_ARD_ }, /* deprecated */
  {  "a.",       SELE_ARD_ }, /* 1 parameter */

  {  "expand",   SELE_EXP_ }, /* 1 parameter */
  {  "x;",       SELE_EXP_ }, /* 1 parameter */
  {  "x.",       SELE_EXP_ }, /* 1 parameter */

  {  "extend",   SELE_EXT_ }, /* 1 parameter */
  {  "xt.",      SELE_EXT_ }, /* 1 parameter */

  {  "name",     SELE_NAMs },
  {  "n;",       SELE_NAMs },/* deprecated */
  {  "n.",       SELE_NAMs },

  {  "symbol",   SELE_ELEs },
  {  "element",  SELE_ELEs },
  {  "elem",     SELE_ELEs },
  {  "e;",       SELE_ELEs },/* deprecated */
  {  "e.",       SELE_ELEs },

  {  "enabled",  SELE_ENAz },

  {  "residue",  SELE_RSIs },
  {  "resi",     SELE_RSIs },
  {  "resident", SELE_RSIs },
  {  "resid",    SELE_RSIs },
  {  "i;",       SELE_RSIs },/* deprecated */
  {  "i.",       SELE_RSIs },

  {  "rep",      SELE_REPs },

  {  "color",    SELE_COLs },
  {  "cartoon_color", SELE_CCLs },
  {  "ribbon_color",  SELE_RCLs },

  {  "altloc",   SELE_ALTs },
  {  "alt",      SELE_ALTs },

  {  "flag",     SELE_FLGs },
  {  "f;",       SELE_FLGs },/* deprecated */
  {  "f.",       SELE_FLGs },
  
  {  "gap",      SELE_GAP_ },

  {  "partial_charge",SELE_PCHx },
  {  "pc;",      SELE_PCHx },/* deprecated */
  {  "pc.",      SELE_PCHx },

  {  "masked",   SELE_MSKz },
  {  "msk.",     SELE_MSKz },

  {  "protected",SELE_PTDz },
  {  "pr.",      SELE_PTDz },

  {  "formal_charge", SELE_FCHx },
  {  "fc;",      SELE_FCHx },/* deprecated */
  {  "fc.",      SELE_FCHx },

  {  "numeric_type",SELE_NTYs },
  {  "nt;",      SELE_NTYs },/* deprecated */
  {  "nt.",      SELE_NTYs },

  {  "text_type",SELE_TTYs },
  {  "tt;",      SELE_TTYs },/* deprecated */
  {  "tt.",      SELE_TTYs },

  {  "chain",    SELE_CHNs },
  {  "c;",       SELE_CHNs },/* deprecated */
  {  "c.",       SELE_CHNs },

  {  cKeywordCenter,   SELE_CENz },
  {  "bonded",   SELE_BNDz },

  {  "segment",  SELE_SEGs },
  {  "segid",    SELE_SEGs },
  {  "segi",     SELE_SEGs },
  {  "s;",       SELE_SEGs },/* deprecated */
  {  "s.",       SELE_SEGs },

  {  "ss",       SELE_SSTs },

  {  "state",    SELE_STAs },

  {  "object",   SELE_MODs },
  {  "o.",       SELE_MODs },

  {  cKeywordOrigin,   SELE_ORIz },

  {  "model",    SELE_MODs },
  {  "m;",       SELE_MODs },/* deprecated */
  {  "m.",       SELE_MODs },

  {  "index",    SELE_IDXs },
  {  "idx.",     SELE_IDXs },

  {  "id",       SELE_ID_s },
  {  "rank",     SELE_RNKs },

  {  "within",   SELE_WIT_ },
  {  "w.",       SELE_WIT_ },

  {  "near_to",  SELE_NTO_ },
  {  "nto."   ,  SELE_NTO_ },

  {  "beyond",   SELE_BEY_ },
  {  "be.",      SELE_BEY_ },

  {  "donors",   SELE_DONz },
  {  "don.",     SELE_DONz },

  {  "acceptors",SELE_ACCz },
  {  "acc."    , SELE_ACCz },

  {  "pepseq",   SELE_PEPs },
  {  "ps.",      SELE_PEPs },
 
  /*
  {  "nucseq",  SELE_NUCs },
  {  "ns.",      SELE_NUCs },
  */

  {  "fixed",    SELE_FXDz },
  {  "fxd.",     SELE_FXDz },

  {  "restrained", SELE_RSTz },
  {  "rst.",       SELE_RSTz },

  {  "polymer",  SELE_POLz },
  {  "pol.",     SELE_POLz },

  {  "organic",  SELE_ORGz },
  {  "org.",     SELE_ORGz },

  {  "inorganic",SELE_INOz },
  {  "ino.",     SELE_INOz },

  {  "solvent",  SELE_SOLz },
  {  "sol.",     SELE_SOLz },

  {  "guide",    SELE_GIDz },

  {  "present",  SELE_PREz },
  {  "pr.",      SELE_PREz },

  {  "resname",  SELE_RSNs },
  {  "resn",     SELE_RSNs },
  {  "r;",       SELE_RSNs },/* deprecated */
  {  "r.",       SELE_RSNs },

  {  "%",        SELE_SELs },
  {  "b",        SELE_BVLx }, /* 2 operand selection operator */ 
  {  "q",        SELE_QVLx }, /* 2 operand selection operator */ 

  {  "bycell",   SELE_BYX1 },
  {  "", 0 }
};

#define SCMP_GTHN 0x01
#define SCMP_LTHN 0x02
#define SCMP_RANG 0x03
#define SCMP_EQAL 0x04

static WordKeyValue AtOper[] = 
{
 { ">",      SCMP_GTHN },
 { "<",      SCMP_LTHN },
 { "in",     SCMP_RANG },
 { "=",      SCMP_EQAL },
 { "", 0 }
};


#define cINTER_ENTRIES 11

#if 0
int SelectorResidueVLAGetDistMat(PyMOLGlobals *G, int *vla, 
                                 int n, int state,float *dist_mat, float scale)
{
  register int a,b;
  register ObjectMolecule *obj;
  register CoordSet *cs;
  float *coord = Alloc(float,3*n);
  if(scale!=0.0F)
    scale = 1.0F/scale;
  if(coord) {
    for(a=0;a<n;a++) {
      register int at_ca1;
      float *vv_ca = coord + a*3;
    
      obj = I->Obj[vla[0]];
      at_ca1 = vla[1];
    
      if(state<obj->NCSet) 
        cs=obj->CSet[state];
      else
        cs=NULL;
      if(cs) {
        register int idx_ca1 = -1;
        if(obj->DiscreteFlag) {
          if(cs==obj->DiscreteCSet[at_ca1])
            idx_ca1=obj->DiscreteAtmToIdx[at_ca1];
          else
            idx_ca1=-1;
        } else 
          idx_ca1=cs->AtmToIdx[at_ca1];
      
        if(idx_ca1>=0) {
          float *v_ca1 = cs->Coord + 3*idx_ca1;
          copy3f(v_ca1,vv_ca);
        }
      }
    }
    {
      for(a=0;a<n;a++) { /* optimize this later */
        float *vv_ca = coord + a*3;
        for(b=0;b<n;b++) {
          float *vv_cb = coord + b*3;          
          float diff = scale* diff3f(vv_ca,vv_cb);
          dist_mat[a][b] = diff;
          dist_mat[b][a] = diff;
        }
      }
    }
  }
  FreeP(coord);
  return 1;
}
#endif

int SelectorRenameObjectAtoms(PyMOLGlobals *G,ObjectMolecule *obj,int sele,int force,int update_table)
{
  int result = 0;
  int obj_nAtom = obj->NAtom;

  if(update_table) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);    
  }
  if(obj_nAtom)  {
    int *flag = Calloc(int,obj_nAtom);
    if(!flag) {
      result = -1;
    } else {
      AtomInfoType *ai=obj->AtomInfo;
      int a;
      for(a=0;a<obj_nAtom;a++) {
	if(SelectorIsMember(G,ai->selEntry,sele)) {
	  flag[a] = true;
	}
	ai++;
      }
      result=ObjectMoleculeRenameAtoms(obj,flag,force);
    }
    FreeP(flag);
  }
  return result;
}

int SelectorResidueVLAsTo3DMatchScores(PyMOLGlobals *G, CMatch *match,
                                       int *vla1,int n1,int state1,
                                       int *vla2,int n2,int state2,
                                       float seq_wt,
                                       float radius, float scale, float base,
                                       float coord_wt, float rms_exp)
{
  register CSelector *I = G->Selector;  
  int a,b,*vla;
  int n_max = (n1>n2) ? n1 : n2;
  float *inter1 = Calloc(float,cINTER_ENTRIES*n1);
  float *inter2 = Calloc(float,cINTER_ENTRIES*n2);
  float *v_ca = Calloc(float,3*n_max);
  if(inter1&&inter2&&v_ca) {
    int pass;

    for(pass=0;pass<2;pass++) {
      register ObjectMolecule *obj;
      register CoordSet *cs;
      register int *neighbor = NULL;
      register AtomInfoType *atomInfo = NULL;
      ObjectMolecule *last_obj = NULL;
      float **dist_mat;
      float *inter;
      int state;
      int n;
      if(!pass) {
        vla = vla1;
        state = state1;
        inter = inter1;
        n = n1;
        dist_mat = match->da;
      } else {
        vla = vla2;
        state = state2;
        inter = inter2;
        n = n2;
        dist_mat = match->db;
      }

      if(state<0) state = 0;
      for(a=0;a<n;a++) {
        register int at_ca1;
        float *vv_ca = v_ca + a*3;

        obj = I->Obj[vla[0]];
        at_ca1 = vla[1];
        if(obj!=last_obj) {
          ObjectMoleculeUpdateNeighbors(obj);
          last_obj = obj;
          neighbor = obj->Neighbor;          
          atomInfo = obj->AtomInfo;
        }

        if(state<obj->NCSet) 
          cs=obj->CSet[state];
        else
          cs=NULL;
        if(cs && neighbor && atomInfo) {
          register int idx_ca1 = -1;
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at_ca1])
              idx_ca1=obj->DiscreteAtmToIdx[at_ca1];
            else
              idx_ca1=-1;
          } else 
            idx_ca1=cs->AtmToIdx[at_ca1];
        
          if(idx_ca1>=0) {
            register int mem0,mem1,mem2,mem3,mem4;
            register int nbr0,nbr1,nbr2,nbr3;
            float *v_ca1 = cs->Coord + 3*idx_ca1;
            int at_cb1,idx_cb1 = -1;
            int cnt = 0;

            copy3f(v_ca1,vv_ca);
            copy3f(v_ca1,inter+8);
          
            /* find attached CB */

            mem0 = at_ca1;
            nbr0 = neighbor[mem0]+1;
            while((mem1 = neighbor[nbr0])>=0) {
              if((atomInfo[mem1].protons==cAN_C) && 
                 (strcmp(atomInfo[mem1].name,"CB")==0)) {
                at_cb1 = mem1;

                if(obj->DiscreteFlag) {
                  if(cs==obj->DiscreteCSet[at_cb1])
                    idx_cb1=obj->DiscreteAtmToIdx[at_cb1];
                  else
                    idx_cb1=-1;
                } else 
                  idx_cb1=cs->AtmToIdx[at_cb1];
                break;
              }
              nbr0+=2;
            }

            /* find remote CA, CB */
          
            if(idx_cb1>=0) {
              float *v_cb1 = cs->Coord + 3*idx_cb1;

              mem0 = at_ca1;
              nbr0 = neighbor[mem0]+1;
              while((mem1 = neighbor[nbr0])>=0) {
              
                nbr1 = neighbor[mem1]+1;
                while((mem2 = neighbor[nbr1])>=0) {
                  if(mem2!=mem0) { 
                    int at_ca2,idx_ca2 = -1;

                    nbr2 = neighbor[mem2]+1;
                    while((mem3 = neighbor[nbr2])>=0) {
                      if((mem3!=mem1)&&(mem3!=mem0)) {
                        if((atomInfo[mem3].protons==cAN_C) && 
                           (strcmp(atomInfo[mem3].name,"CA")==0)) {
                          at_ca2 = mem3;
                          if(obj->DiscreteFlag) {
                            if(cs==obj->DiscreteCSet[at_ca2])
                              idx_ca2=obj->DiscreteAtmToIdx[at_ca2];
                            else
                              idx_ca2=-1;
                          } else 
                            idx_ca2=cs->AtmToIdx[at_ca2];
                          break;
                        }
                      }
                      nbr2+=2;
                    }
                    if(idx_ca2>=0) {
                      float *v_ca2 = cs->Coord + 3*idx_ca2;
                  
                      nbr2 = neighbor[mem2]+1;
                      while((mem3 = neighbor[nbr2])>=0) {
                        if((mem3!=mem1)&&(mem3!=mem0)) {
                          int idx_cb2 = -1;
                          nbr3 = neighbor[mem3]+1;
                          while((mem4 = neighbor[nbr3])>=0) {
                            if((mem4!=mem2)&&(mem4!=mem1)&&(mem4!=mem0)) {
                              if((atomInfo[mem4].protons==cAN_C) && 
                                 (strcmp(atomInfo[mem4].name,"CB")==0)) {
                                int at_cb2 = mem4;
                                if(obj->DiscreteFlag) {
                                  if(cs==obj->DiscreteCSet[at_cb2])
                                    idx_cb2=obj->DiscreteAtmToIdx[at_cb2];
                                  else
                                    idx_cb2=-1;
                                } else 
                                  idx_cb2=cs->AtmToIdx[at_cb2];
                                break;
                              }
                            }
                            nbr3+=2;
                          }
                        
                          if(idx_cb2>=0) {
                            float *v_cb2 = NULL;
                            v_cb2 = cs->Coord + 3*idx_cb2;
                            {
                              float angle = get_dihedral3f( v_cb1, v_ca1, v_ca2, v_cb2 );
                              if(idx_cb1<idx_cb2) {
                                inter[0] = (float)cos(angle);
                                inter[1] = (float)sin(angle);
                              } else {
                                inter[2] = (float)cos(angle);
                                inter[3] = (float)sin(angle);
                              }
                            }
                            cnt++;
                          }
                        }
                        nbr2+=2;
                      }
                    }
                  }
                  nbr1+=2;
                }
                nbr0+=2;
              }
            }
          }
        }
        vla+=3;
        inter+=cINTER_ENTRIES;
      }
      if(dist_mat) {
        for(a=0;a<n;a++) { /* optimize this later */
          float *vv_ca = v_ca + a*3;
          for(b=0;b<n;b++) {
            float *vv_cb = v_ca + b*3;          
            float diff = (float)diff3f(vv_ca,vv_cb);
            dist_mat[a][b] = diff;
            dist_mat[b][a] = diff;
          }
        }
      }
      {
        MapType *map=MapNew(G,radius,v_ca,n, NULL);
        if(!pass) {
          inter = inter1;
        } else {
          inter = inter2;
        }
        if(map) {
          int i,h,k,l;
          MapSetupExpress(map);
        
          for(a=0;a<n;a++) {
            float *v_ca1 = v_ca + 3 * a;
            float *i_ca1 = inter + cINTER_ENTRIES*a;
            if(MapExclLocus(map,v_ca1,&h,&k,&l)) {
              i=*(MapEStart(map,h,k,l));
              if(i) {
                b=map->EList[i++];
                while(b>=0) {
                  float *v_ca2 = v_ca + 3 * b;
                  if(a != b) {
                    if(within3f(v_ca1,v_ca2,radius)) {
                      float *i_ca2  = inter + cINTER_ENTRIES*b;
                      i_ca1[4] += i_ca2[0]; /* add dihedral vectors head-to-tail */
                      i_ca1[5] += i_ca2[1];
                      i_ca1[6] += i_ca2[2];
                      i_ca1[7] += i_ca2[3];
                    }
                  }
                  b=map->EList[i++];
                }
              }
            }
          }
          MapFree(map);
          for(a=0;a<n;a++) {
            float nf = (float)sqrt(inter[4]*inter[4] + inter[5]*inter[5]);
            if(nf>0.0001F) {
              inter[4] = inter[4] / nf;
              inter[5] = inter[5] / nf;
            }
            nf = (float)sqrt(inter[6]*inter[6] + inter[7]*inter[7]);
            if(nf>0.0001F) {
            
              inter[6] = inter[6] / nf;
              inter[7] = inter[7] / nf;
            }
            inter+=cINTER_ENTRIES;
          }
        }
      }
    }
    {
      const float _0F = 0.0F;

      if((scale!=0.0F)||(seq_wt!=0.0F)) {
        for(a=0;a<n1;a++) {
          float *i1 = inter1+ cINTER_ENTRIES*a;
          for(b=0;b<n2;b++) {
            float *i2 = inter2 + cINTER_ENTRIES*b;
            float sm[cINTER_ENTRIES], comp1, comp2, comp3 = 1.0F;
            float score;
            int c;
            for(c=0;c<(cINTER_ENTRIES-1);c+=2) {
              if( ((i1[c] == _0F) && (i1[c+1]== _0F)) ||
                  ((i2[c] == _0F) && (i2[c+1]== _0F))) { /* handle glycine case */
                sm[c] = 1.0F;
                sm[c+1] = 1.0F;
              } else {
                sm[c] = i1[c]+i2[c];
                sm[c+1] = i1[c+1]+i2[c+1];
              }
            }
            comp1 = (float)
              ((sqrt(sm[0]*sm[0] + sm[1]*sm[1]) + 
                sqrt(sm[2]*sm[2] + sm[3]*sm[3])) * 0.25);
            comp2 = (float)
              ((sqrt(sm[4]*sm[4] + sm[5]*sm[5]) +
                sqrt(sm[6]*sm[6] + sm[7]*sm[7])) * 0.25);
            score = scale*(comp1*comp2 - base);
            if(coord_wt!=0.0) {
              float diff = (float)diff3f(i1+8,i2+8);
              comp3 = (float)-log(diff/rms_exp);
              score = (1-coord_wt)*score + coord_wt*comp3*scale;
            }
            match->mat[a][b] = seq_wt*match->mat[a][b] + score;
          }
        }
      }
    }
  }
  FreeP(inter1);
  FreeP(inter2);
  FreeP(v_ca);
  return 1;
}

int SelectorNameIsKeyword(PyMOLGlobals *G, char *name)
{
  register CSelector *I = G->Selector;
  WordType lower_name;
  OVreturn_word result;
  UtilNCopyToLower(lower_name,name,sizeof(WordType));
  if(OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,lower_name)))) {
    if(OVreturn_IS_OK( (result = OVOneToAny_GetKey(I->Key, result.word)))) {
      return 1;
    }
  }
  return 0;
}
/*========================================================================*/
static int SelectorIsSelectionDiscrete(PyMOLGlobals *G,int sele,int update_table)
{
  register CSelector *I=G->Selector;
  register ObjectMolecule **i_obj=I->Obj,*obj;
  AtomInfoType *ai;
  int result=false;
  register TableRec *i_table = I->Table, *table_a;
  ov_size a;

  if(update_table) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  }

  table_a = i_table + cNDummyAtoms;
  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    obj= i_obj[table_a->model];
    ai = obj->AtomInfo + table_a->atom;
    if(SelectorIsMember(G,ai->selEntry,sele)) {
      if(obj->DiscreteFlag) {
        result=true;
        break;
      }
    }
    table_a++;
  }
  return(result);
}

static int *SelectorUpdateTableMultiObjectIdxTag(PyMOLGlobals *G,
                                                 ObjectMolecule **obj_list,
                                                 int no_dummies,
                                                 int **idx_list,int *n_idx_list,
                                                 int n_obj)
{
  int a=0;
  int b;
  int c=0;
  int modelCnt;
  int *result = NULL;
  register CSelector *I=G->Selector;
  ObjectMolecule *obj = NULL;
  int *idx,n_idx;

  PRINTFD(G,FB_Selector)
    "SelectorUpdateTableMultiObject-Debug: entered ...\n"
    ENDFD;
 
  SelectorClean(G);

  I->SeleBaseOffsetsValid = true; /* all states -> all atoms -> offsets valid */
  I->NCSet = 0;
  if(no_dummies) {
    modelCnt = 0;
    c = 0;
  } else {
    modelCnt=cNDummyModels;
    c=cNDummyAtoms;
  }
  for(b=0;b<n_obj;b++) {
    obj = obj_list[b];
    c+=obj->NAtom;
    if(I->NCSet<obj->NCSet) I->NCSet=obj->NCSet;
    modelCnt++;
  }
  result = Calloc(int,c);
  I->Table=Calloc(TableRec,c);
  ErrChkPtr(G,I->Table);
  I->Obj=Calloc(ObjectMolecule*,modelCnt);
  ErrChkPtr(G,I->Obj);
  if(no_dummies) {
    modelCnt = 0;
    c = 0;
  } else {
    c=cNDummyAtoms;
    modelCnt=cNDummyModels;
  }
  for(b=0;b<n_obj;b++) {
    obj = obj_list[b];
    idx = idx_list[b];
    n_idx = n_idx_list[b];
    
    I->Obj[modelCnt]=obj;
    obj->SeleBase=c; 
    for(a=0;a<obj->NAtom;a++) {
      I->Table[c].model=modelCnt;
      I->Table[c].atom=a;
      c++;
    }
    if(idx&&n_idx) {
      if(n_idx>0) {
        for(a=0; a< n_idx; a++) {
          int at = idx[2*a]; /* index first */
          int pri = idx[2*a+1]; /* then priority */
          if((at>=0)&&(at<obj->NAtom)) { 
            result[obj->SeleBase + at] = pri;
          }
        }
      }
    }
    modelCnt++;
    I->NModel=modelCnt;
  }
  I->NAtom=c;
  I->Flag1=Alloc(int,c);
  ErrChkPtr(G,I->Flag1);
  I->Flag2=Alloc(int,c);
  ErrChkPtr(G,I->Flag2);
  I->Vertex=Alloc(float,c*3);
  ErrChkPtr(G,I->Vertex);
  
  PRINTFD(G,FB_Selector)
    "SelectorUpdateTableMultiObject-Debug: leaving...\n"
    ENDFD;
  
  return(result);
}

static int IntInOrder(int *list,int a,int b)
{
  return(list[a]<=list[b]);
}

static int SelectorAddName(PyMOLGlobals *G, int index)
{
  register CSelector *I = G->Selector;
  int ok=false;
  OVreturn_word result;
  OVstatus status;
  if(OVreturn_IS_OK( (result = OVLexicon_GetFromCString(I->Lex,I->Name[index])))) {
    if(OVreturn_IS_OK( (status = OVOneToOne_Set(I->NameOffset, result.word, index)))) {
      ok=true;
    }
  }
  return ok;
}

static int SelectorDelName(PyMOLGlobals *G, int index)
{
  register CSelector *I = G->Selector;
  int ok=false;
  OVreturn_word result;
  if(OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,I->Name[index])))) {
    if(OVreturn_IS_OK(OVLexicon_DecRef(I->Lex, result.word)) &&
       OVreturn_IS_OK(OVOneToOne_DelForward(I->NameOffset, result.word))) {
      ok=true;
    }
  }
  return ok;
}

int SelectorClassifyAtoms(PyMOLGlobals *G,int sele, int preserve,ObjectMolecule *only_object)
{
  register CSelector *I=G->Selector;
  ObjectMolecule *obj,*last_obj=NULL,*obj0,*obj1 = NULL;
  int a,aa,at,a0,a1;
  AtomInfoType *ai,*last_ai=NULL,*ai0,*ai1;
  unsigned int mask;
  int n_dummies = 0;

  if(only_object) {
    SelectorUpdateTableSingleObject(G,only_object,cSelectorUpdateTableAllStates,
                                    true,NULL,0,false);  
    n_dummies = 0;
  } else {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
    n_dummies = cNDummyAtoms;
  }
  a=0;
  while(a<I->NAtom) {
    obj = I->Obj[I->Table[a].model];
    at = I->Table[a].atom;
    ai = obj->AtomInfo + at; 
    
    if(SelectorIsMember(G,ai->selEntry,sele) &&
       ((!AtomInfoSameResidueP(G,ai,last_ai)))) {
      
      AtomInfoType *guide_atom = NULL;
      
      /* delimit residue */
      
      a0 = a-1;
      while(a0>=n_dummies) {
        obj0 = I->Obj[I->Table[a0].model];
        if(obj0!=obj) 
          break;
        ai0 = obj0->AtomInfo + I->Table[a0].atom;
        if(!AtomInfoSameResidue(G,ai0,ai))
          break;
        a0--;
      }
      
      a1 = a+1;
      while(a1<I->NAtom) {
        obj1 = I->Obj[I->Table[a1].model];
        if(obj1!=obj) 
          break;
        ai1 = obj1->AtomInfo + I->Table[a1].atom;
        if(!AtomInfoSameResidue(G,ai1,ai))
          break;
        a1++;
      }
      
      a0++;
      a1--;

      mask = 0;
      if(AtomInfoKnownPolymerResName(G,ai->resn))
        mask = cAtomFlag_polymer;
      else if(AtomInfoKnownWaterResName(G,ai->resn))
        mask = cAtomFlag_solvent;
      else {

        /* does this residue have a canonical atoms? */
        
        int found_ca = false;
        int found_n = false;
        int found_c = false;
        int found_o = false;
        int found_oh2 = false;
        int found_carbon = false;
        int found_cn_bond = false;
        int found_nc_bond = false;
        int found_o3_bond = false;
        int found_o3star = false;
        int found_c3star = false;
        int found_c4star = false;
        int found_c5star = false;
        int found_o5star = false;
        int found_p_bond = false;
        if(obj!=last_obj) {
          ObjectMoleculeUpdateNeighbors(obj);
          last_obj = obj;
        }

        ai0 = obj->AtomInfo + I->Table[a0].atom;
        for(aa=a0;aa<=a1;aa++) { 
          if(ai0->protons == cAN_C) {
            register char *name = ai0->name;
            found_carbon = true;
            switch(name[0]) {
            case 'C':
              switch(name[1]) {
              case 0:
                found_c = true;
                found_cn_bond = ObjectMoleculeIsAtomBondedToName(obj,I->Table[aa].atom,"N");
                break;
              case 'A':
                switch(name[2]) {
                case 0:
                  found_ca = true;
                  guide_atom = ai0;
                  break;
                }
              case '3':
                switch(name[2]) {
                case '*':
                case '\'':
                  guide_atom = ai0;
                  found_c3star = true;
                  break;
                }
                break;
              case '4':
                switch(name[2]) {
                case '*':
                case '\'':
                  found_c4star = true;
                  break;
                }
                break;
              case '5':
                switch(name[2]) {
                case '*':
                case '\'':
                  found_c5star = true;
                  break;
                }
                break;
              }
            }
          } else if(ai0->protons == cAN_N) {
            register char *name = ai0->name;
            switch(name[0]) {
            case 'N':
              switch(name[1]) {
              case 0:
                found_n = true;
                found_nc_bond = ObjectMoleculeIsAtomBondedToName(obj,I->Table[aa].atom,"C");
                break;
              }
            }
          } else if(ai0->protons == cAN_O) {
            register char *name = ai0->name;
            switch(name[0]) {
            case 'O':
              switch(name[1]) {
              case 0:
                found_o = true;
                break;
              case 'H':
                switch(name[2]) {
                case '2':
                  found_oh2 = true;
                  break;
                }
              case '3':
                switch(name[2]) {
                case '*':
                case '\'':
                  found_o3star = true;
                  found_o3_bond = ObjectMoleculeIsAtomBondedToName(obj,I->Table[aa].atom,"P");
                  break;
                }
                break;
              case '5':
                switch(name[2]) {
                case '*':
                case '\'':
                  found_o5star = true;
                  break;
                }
                break;
              }
            }
          } else if(ai0->protons == cAN_P) {

            register char *name = ai0->name;
            switch(name[0]) {
            case 'P':
              switch(name[1]) {
              case 0:
                found_p_bond = (ObjectMoleculeIsAtomBondedToName(obj,I->Table[aa].atom,"O3*") ||
                                ObjectMoleculeIsAtomBondedToName(obj,I->Table[aa].atom,"O3'"));
                break;
              }
            }
          }
          ai0++;
        }
        
        if((found_ca && found_n && found_c && found_o && (found_cn_bond||found_nc_bond))||
           (found_o3star && found_c3star && found_c4star && found_c5star && found_o5star &&
            (found_o3_bond||found_p_bond))) {
          mask = cAtomFlag_polymer;
        } else if(found_carbon)
          mask = cAtomFlag_organic;
        else if((found_o||found_oh2)&&(a1==a0))
          mask = cAtomFlag_solvent; 
        else
          mask = cAtomFlag_inorganic;
      }

      /* mark which atoms we can write to */
      
      ai0 = obj->AtomInfo + I->Table[a0].atom;
      if(preserve) {
        for(aa=a0;aa<=a1;aa++) { 
          if(SelectorIsMember(G,ai0->selEntry,sele))
            if(!(ai0->flags & cAtomFlag_class_mask))
              ai0->flags = (ai0->flags & cAtomFlag_class_mask) | mask;
          ai0++;
        }
      } else {
        for(aa=a0;aa<=a1;aa++) { 
          if(SelectorIsMember(G,ai0->selEntry,sele))
            ai0->flags = (ai0->flags & cAtomFlag_class_mask) | mask;
          ai0++;
        }
      }

      if((!guide_atom)&&(mask==cAtomFlag_polymer)) {
        ai0 = obj->AtomInfo + I->Table[a0].atom;
        for(aa=a0;aa<=a1;aa++) { 
          if(ai0->protons == cAN_C) {
            register char *name = ai0->name;
            switch(name[0]) {
            case 'C':
              switch(name[1]) {
              case 'A':
                switch(name[2]) {
                case 0:
                  guide_atom = ai0;
                  break;
                }
              case '4':
                switch(name[2]) { /* use C4* as guide atom for nucleic acids */
                case '*':
                case '\'':
                  guide_atom = ai0;
                  break;
                }
                break;
              }
            }
          }
          ai0++;
        }
      }

      if(guide_atom)
        guide_atom->flags |= cAtomFlag_guide;
      
      if(a1>(a+1))
        a = a1;
    }
    a++;
  }
  return true;
}

static void SelectionInfoInit(SelectionInfoRec *rec)
{
  rec->justOneObjectFlag = false;
  rec->justOneAtomFlag = false;
}

void SelectorComputeFragPos(PyMOLGlobals *G,ObjectMolecule *obj,int state,int n_frag, char *prefix,float **vla)
{
  register CSelector *I=G->Selector;
  WordType name;
  int *sele;
  int *cnt;
  SelectorUpdateTableSingleObject(G,obj,cSelectorUpdateTableAllStates,true,NULL,0,false);
  sele = Alloc(int,n_frag);
  cnt = Calloc(int,n_frag);
  VLACheck(*vla,float,n_frag*3+2);
  {
    int a;
    for(a=0;a<n_frag;a++) {
      sprintf(name,"%s%d",prefix,a+1);
      sele[a] = SelectorIndexByName(G,name);
      zero3f((*vla)+3*a);
    }
  }
  {
    int at,a,ati;
    AtomInfoType *ai;
    float v[3],*vp;
    int vert_flag;
    for(at=0;at<I->NAtom;at++) {
      
      ati = I->Table[at].atom;
      ai = obj->AtomInfo + ati; 
      
      vert_flag = false;
      for(a=0;a<n_frag;a++) {
        if(SelectorIsMember(G,ai->selEntry,sele[a])) {
          if(!vert_flag) {
            vert_flag=ObjectMoleculeGetAtomVertex(obj,state,ati,v);
          }
          if(vert_flag) {
            vp = (*vla)+3*a;
            add3f(v,vp,vp);
            cnt[a]++;
          }
        }
      }
    }
  }

  {
    int a;
    float div,*vp;
    for(a=0;a<n_frag;a++) {
      if(cnt[a]) {
        vp = (*vla)+3*a;        
        div = 1.0F/cnt[a];
        scale3f(vp,div,vp);
      }
    }
  }

  FreeP(sele);
  FreeP(cnt);
}

MapType *SelectorGetSpacialMapFromSeleCoord(PyMOLGlobals *G,int sele,int state,float cutoff,float **coord_vla)
{
  register CSelector *I=G->Selector;
  int *index_vla = NULL;
  float *coord = NULL;
  int n,nc=0;
  MapType *result = NULL;
  if(sele<0)
    return NULL;
  else {
    
    SelectorUpdateTable(G,state,-1);
    index_vla = SelectorGetIndexVLA(G,sele);

    if(index_vla) {
      n=VLAGetSize(index_vla);
      if(n) 
        coord = VLAlloc(float,n*3);
      if(coord) {
        int i,a;
        int st,sta;
        ObjectMolecule *obj;
        CoordSet *cs;
        int at;
        AtomInfoType *ai;
        int idx;
        float *src,*dst;
        for(i=0;i<n;i++) {
          a=index_vla[i];

          obj = I->Obj[I->Table[a].model];
          at = + I->Table[a].atom;
          ai = obj->AtomInfo + at;
          for(st=0;st<I->NCSet;st++) {

            if((state<0)||(st==state)) {

              sta = st;
              if(sta<obj->NCSet) 
                cs=obj->CSet[sta];
              else
                cs=NULL;
              if(cs) {
                if(obj->DiscreteFlag) {
                  if(cs==obj->DiscreteCSet[at])
                    idx=obj->DiscreteAtmToIdx[at];
                  else
                    idx=-1;
                } else 
                  idx=cs->AtmToIdx[at];
              } else {
                idx = -1;
              }
              if(idx>=0) {
                VLACheck(coord,float,nc*3+2);
                src = cs->Coord + 3*idx;
                dst = coord + 3*nc;
                *(dst++) = *(src++);
                *(dst++) = *(src++);
                *(dst++) = *(src++); 
                nc++;

              }
            }
          }
        }
        if(nc) {
          result = MapNew(G,cutoff,coord, nc,NULL);
        }
      }
    }
  }
  VLAFreeP(index_vla);
  if(coord)
    VLASize(coord,float,nc*3);
  *(coord_vla)=coord;
  return(result);
}

static ov_diff SelectGetNameOffset(PyMOLGlobals *G,char *name,ov_size minMatch,int ignCase)
{

  register CSelector *I=G->Selector;
  int result = -1;
  while(name[0]=='?')
    name++;
  { /* first try for perfect match using the dictionary */
    OVreturn_word res;
    if( OVreturn_IS_OK( (res = OVLexicon_BorrowFromCString(I->Lex,name)))) {
      if( OVreturn_IS_OK( (res = OVOneToOne_GetForward(I->NameOffset, res.word)))) { 
        result = res.word;
      }
    }
  }
  if(result<0) { /* not found, so try partial/ignored-case match */
    
    register int offset,wm,best_match,best_offset;
    SelectorWordType *I_Name = I->Name;
    offset=0;
    best_offset=-1;
    best_match=-1;
    while(name[0]=='?')
      name++;
    
    while(I_Name[offset][0]) {
      wm=WordMatch(G,name,I_Name[offset],ignCase);
      if(wm<0) { /* exact match is always good */
        best_offset=offset;
        best_match=wm;
        break;
      }
      if(wm>0) {
        if(best_match<wm) {
          best_match=wm;
          best_offset=offset;
        } else if(best_match==wm) { /* uh oh -- ambiguous match */
          best_offset = -1;
        }
      }
      offset++;
    }
    if((best_match<0)||(best_match>(ov_diff)minMatch))
      result=best_offset;
  }
  return(result);  
}

void SelectorSelectByID(PyMOLGlobals *G,char *name,ObjectMolecule *obj,int *id,int n_id)
{
  register CSelector *I=G->Selector;
  int min_id,max_id,range,*lookup = NULL;
  int *atom = NULL;
  /* this routine only works if IDs cover a reasonable range --
     should rewrite using a hash table */

  SelectorUpdateTableSingleObject(G,obj,cSelectorUpdateTableAllStates,true,NULL,0,false);
  atom = Calloc(int,I->NAtom);
  if(I->NAtom) {

    /* determine range */

    {
      int a,cur_id;
      cur_id = obj->AtomInfo[0].id;
      min_id = cur_id;
      max_id = cur_id;
      for(a=1;a<obj->NAtom;a++) {
        cur_id = obj->AtomInfo[a].id;
        if(min_id>cur_id) min_id = cur_id;
        if(max_id<cur_id) max_id = cur_id;
      }
    }

    /* create cross-reference table */

    {
      int a,offset;
      
      range = max_id - min_id + 1;
      lookup = Calloc(int,range);
      for(a=0;a<obj->NAtom;a++) {
        offset = obj->AtomInfo[a].id - min_id;
        if(lookup[offset])
          lookup[offset] = -1;
        else {
          lookup[offset] = a+1;
        }
      }
    }
    
    /* iterate through IDs and mark */

    {
      int i,a,offset,lkup;

      for(i=0;i<n_id;i++) {
        offset = id[i]-min_id;
        if((offset>=0)&&(offset<range)) {
          lkup = lookup[offset];
          if(lkup>0) {
            atom[lkup-1]=true;
          } else if(lkup<0) { 
            for(a=0;a<obj->NAtom;a++) {
              if(obj->AtomInfo[a].id==id[i]) 
                atom[a]=true;
            }
          }
        }
      }
    }
  }

  SelectorEmbedSelection(G,atom,name,NULL,true, -1);
  FreeP(atom);
  FreeP(lookup);
  SelectorClean(G);
}

void SelectorDefragment(PyMOLGlobals *G) 
{
  register CSelector *I=G->Selector;
  /* restore new member ordering so that CPU can continue to get good cache hit */

  int n_free = 0;
  int m;
  int *list,*l;
  int a;
  m = I->FreeMember;
  while(m) {
    n_free++;
    m = I->Member[m].next;
  }
  if(n_free) {
    list = Alloc(int,n_free);
    l=list;
    m = I->FreeMember;
    while(m) {
      *(l++) = m;
      m = I->Member[m].next;
    }
    UtilSortInPlace(G,list,n_free,sizeof(int),(UtilOrderFn*)IntInOrder);
    while(n_free>5000) { /* compact inactive members when possible */
      if(list[n_free-1]==I->NMember) {
        I->NMember--;
        n_free--;
      } else 
        break;
    }
    for(a=0;a<(n_free-1);a++) {
      I->Member[list[a]].next = list[a+1];
    }
    I->Member[list[n_free-1]].next = 0;
    I->FreeMember = list[0];
    FreeP(list);
  }
}

int SelectorWalkTree(PyMOLGlobals *G,int *atom,int *comp,int *toDo,int **stk,
                     int stkDepth,ObjectMolecule *obj,
                     int sele1,int sele2,int sele3, int sele4);

typedef struct {
  int color;
  int sele;
} ColorectionRec;

static void SelectorDeleteSeleAtOffset(PyMOLGlobals *G,int n)
{
  register CSelector *I=G->Selector;
  int id;
  id = I->Info[n].ID;
  SelectorDelName(G,n);
  SelectorPurgeMembers(G,id);
  
  I->NActive--;
  {
    OVreturn_word result;
    /* repoint the name index at the relocated entry */
    if(OVreturn_IS_OK( result = OVOneToOne_GetReverse(I->NameOffset, I->NActive))) {
      OVOneToOne_DelForward(I->NameOffset, result.word);
      OVOneToOne_Set(I->NameOffset, result.word, n);
    }
    strcpy(I->Name[n],I->Name[I->NActive]);
    I->Info[n] = I->Info[I->NActive];
    I->Name[I->NActive][0]=0;
  }
}
char *SelectorGetNameFromIndex(PyMOLGlobals *G,int index)
{
  register CSelector *I=G->Selector;
   int a;
  for(a=1;a<I->NActive;a++) {
    if(I->Info[a].ID==index) {
      return I->Name[a];
    }
  }
  return NULL;
}

#ifndef _PYMOL_NOPY
static void SelectorDeleteIndex(PyMOLGlobals *G,int index)
{
  register CSelector *I=G->Selector;
  int n=0;
  int a;
  for(a=1;a<I->NActive;a++) {
    if(I->Info[a].ID==index) {
      n=a;
      break;
    }
  }
  if(n) 
    SelectorDeleteSeleAtOffset(G,n);
}
#endif

#define cSSMaxHBond 6

#define cSSHelix3HBond          0x0001
#define cSSHelix4HBond          0x0002
#define cSSHelix5HBond          0x0004
#define cSSGotPhiPsi            0x0008
#define cSSPhiPsiHelix          0x0010
#define cSSPhiPsiNotHelix       0x0020
#define cSSPhiPsiStrand         0x0040
#define cSSPhiPsiNotStrand      0x0080
#define cSSAntiStrandSingleHB   0x0100
#define cSSAntiStrandDoubleHB   0x0200
#define cSSAntiStrandBuldgeHB   0x0400
#define cSSAntiStrandSkip       0x0800
#define cSSParaStrandSingleHB   0x1000
#define cSSParaStrandDoubleHB   0x2000
#define cSSParaStrandSkip       0x4000

#define cSSBreakSize 5
    

typedef struct {
  int real;
  int ca,n,c,o; /* indices in selection-table space */
  float phi,psi;
  char ss, ss_save;
  int flags;
  int n_acc,n_don;
  int acc[cSSMaxHBond]; /* interactions where this residue is an acceptor */
  int don[cSSMaxHBond]; /* interactions where this residue is a donor */
  ObjectMolecule *obj;
  int preserve;
  int present;
} SSResi;

int SelectorAssignSS(PyMOLGlobals *G,int target,int present,
                     int state_value,int preserve,ObjectMolecule *single_object,int quiet)
{

  /* PyMOL's secondary structure assignment algorithm: 

  General principal -- if it looks like a duck, then it's a duck:

     I. Helices
       - must have reasonably helical geometry within the helical span
       - near-ideal geometry guarantees helix assignment
       - a continuous ladder stre i+3, i+4, or i+5 hydrogen bonding
        with permissible geometry can reinforce marginal cases
       - a minimum helix is three residues with i+3 H-bond

     II. Sheets
       - Hydrogen bonding ladders are the primary guide
       - Out-of-the envelope 
       - 1-residue gaps in sheets are filled unless there
         is a turn.
  */

  register CSelector *I=G->Selector;
  SSResi *res;
  int n_res = 0;
  int state_start,state_stop,state;
  int consensus = true;
  int first_last_only = false;
  int first_pass = true;

  if(!single_object) {
    if(state_value<0) {
      switch(state_value) {
      case -2: /* api: state=-1: current global state */
      case -3: /* api: state=-2: effective object states TO DO! */
        SelectorUpdateTable(G,state_value,-1);
        break;
      default:
        SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
        break;
      }
    } else {
      SelectorUpdateTable(G,state_value,-1);
    }
  } else {
    SelectorUpdateTableSingleObject(G,single_object,state_value,false,NULL,0,0);
  }

  res = VLACalloc(SSResi,1000);

  if(state_value<0) {
    if(state_value == -4)
      consensus = false;
    if(state_value == -5)
      first_last_only = true;
    state_start = 0;
    state_stop = SelectorGetSeleNCSet(G,target);
  } else {
    state_start=state_value;
    state_stop=state_value+1;
  }
  for(state=state_start;state<state_stop;state++) {
    int a;
    ObjectMolecule *obj;
    int aa,a0,a1,at,idx;
    AtomInfoType *ai,*ai0,*ai1;
    CoordSet *cs;
    ObjectMolecule *last_obj = NULL;
    /* first, we need to count the number of residues under consideration */

    if(first_pass) {
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        
        obj = I->Obj[I->Table[a].model];
        at = + I->Table[a].atom;
        ai = obj->AtomInfo + at;
        
        /* see if CA coordinates exist...*/
        
        if(SelectorIsMember(G,ai->selEntry,present)) {
          
          if((ai->protons==cAN_C)&&
             (WordMatch(G,"CA",ai->name,true)<0)) {
            
            if(last_obj!=obj) {
              ObjectMoleculeUpdateNeighbors(obj);
              ObjectMoleculeVerifyChemistry(obj,state_value);
              last_obj=obj;
            }
            /* delimit residue */
            
            a0 = a-1;
            while(a0>=cNDummyAtoms) {
              ai0 = I->Obj[I->Table[a0].model]->AtomInfo + I->Table[a0].atom;
              if(!AtomInfoSameResidue(G,ai0,ai))
                break;
              a0--;
            }
            
            a1 = a+1;
            while(a1<I->NAtom) {
              ai1 = I->Obj[I->Table[a1].model]->AtomInfo + I->Table[a1].atom;
              if(!AtomInfoSameResidue(G,ai1,ai))
                break;
              a1++;
            }
            
            {
              int found_N = 0;
              int found_O = 0;
              int found_C = 0;
              
              /* locate key atoms */
              
              for(aa=a0+1;aa<a1;aa++) {
                ai = I->Obj[I->Table[aa].model]->AtomInfo + I->Table[aa].atom;
                if((ai->protons==cAN_C)&&
                   (WordMatch(G,"C",ai->name,true)<0)) {
                  found_C = aa;
                }
                if((ai->protons==cAN_N)&&
                   (WordMatch(G,"N",ai->name,true)<0)) {
                  found_N = aa;
                }
                if((ai->protons==cAN_O)&&
                   (WordMatch(G,"O",ai->name,true)<0)) {
                  found_O = aa;
                }
              }
              
              if((found_C)&&
                 (found_N)&&
                 (found_O)) {
                
                VLACheck(res,SSResi,n_res);
                res[n_res].n = found_N;
                res[n_res].o = found_O;
                res[n_res].c = found_C;
                res[n_res].ca = a;
                res[n_res].obj = I->Obj[I->Table[a].model];
                res[n_res].real = true;

                n_res++;
                
              } else {
                if(!quiet) {
                  PRINTFB(G,FB_Selector,FB_Warnings)
                    " AssignSS-Warning: Ignoring incomplete residue %s/%s/%s/%s/...\n",
                    obj->Obj.Name,ai->segi,ai->chain,ai->resi
                    ENDFB(G);
                }
              }
            }
          }
        }
      } /* count pass */

      if(preserve) { /* if we're in preserve mode, then mark which objects don't get changed */
        int a,b;
        char ss;
        ObjectMolecule *p_obj = NULL;
        SSResi *r,*r2;
        for(a=0;a<n_res;a++) {
          r = res+a;
          if(r->real) {
            if(p_obj!=r->obj) {
              ss = r->obj->AtomInfo[I->Table[r->ca].atom].ssType[0];
              if((ss=='S')||(ss=='H')||(ss=='s')||(ss=='h')) {
                p_obj = r->obj;
                
                b=a;
                while(b>=0) {
                  r2 = res + b;
                  if(p_obj==r2->obj) r2->preserve=true;
                  b--;
                }
                b=a+1;
                while(b<n_res) {
                  r2 = res + b;
                  if(p_obj==r2->obj) r2->preserve=true;
                  b++;
                }
              }
            }
          }
        }
      }
      /*  printf("n_res %d\n",n_res);*/

      /* now, let's repack res. into discrete chunks so that we can do easy gap & ladder analysis */
  
      {
        SSResi *res2;
        int a;
        int n_res2 = 0;
        int add_break;
        int at_ca0,at_ca1;

        res2 = VLACalloc(SSResi,n_res*2);
    
        for(a=0;a<n_res;a++) {
          add_break = false;

          if(!a) { 
            add_break=true;
          } else if(res[a].obj!=res[a-1].obj) {
            add_break=true;
          } else if(res[a].obj) {
            at_ca0 = I->Table[res[a].ca].atom;
            at_ca1 = I->Table[res[a-1].ca].atom;
            if(!ObjectMoleculeCheckBondSep(res[a].obj,at_ca0,at_ca1,3)) { /* CA->N->C->CA = 3 bonds */
              add_break=true;
            }
          }
      
          if(add_break) {
            n_res2 += cSSBreakSize;
          }
      
          VLACheck(res2,SSResi,n_res2);
          res2[n_res2] = res[a];
          n_res2++;
        }

        n_res2+=cSSBreakSize;
        VLACheck(res2,SSResi,n_res2);

        VLAFreeP(res);
        res=res2;
        n_res = n_res2;
      }
      first_pass = false;
    }

    /* okay, the rest of this loop runs for each coordinate set */
  
    {
      int b;
      for(a=0;a<n_res;a++) {
        res[a].present = res[a].real;

        if(res[a].present) {
          obj = res[a].obj;
          if(state<obj->NCSet) 
            cs=obj->CSet[state];
          else
            cs=NULL;
          for(b=0;b<4;b++) {
            if(cs) {
              switch(b) {
              case 0: at = I->Table[res[a].n].atom; break;
              case 1: at = I->Table[res[a].o].atom; break;
              case 2: at = I->Table[res[a].c].atom; break;
              default:
              case 3: at = I->Table[res[a].ca].atom; break;              
              }
              if(obj->DiscreteFlag) {
                if(cs==obj->DiscreteCSet[at])
                  idx=obj->DiscreteAtmToIdx[at];
                else
                  idx=-1;
              } else 
                idx=cs->AtmToIdx[at];
            } else 
              idx = -1;
            if(idx<0) {
              res[a].present = false;
            }
          }
        }
      }
    }
      

    /* next, we need to record hydrogen bonding relationships */

    {
    
      MapType *map;
      float *v0,*v1;
      int n1;
      int c,i,h,k,l;
      int at;
      int idx;

      int a,aa;
      int a0,a1; /* SS res space */
      int as0,as1; /* selection space */
      int at0,at1; /* object-atom space */
      int exclude;

      ObjectMolecule *obj0,*obj1;

      CoordSet *cs;
      float cutoff;
      HBondCriteria hbcRec,*hbc;
      int *zero=NULL,*scratch=NULL;


      {
        int max_n_atom = I->NAtom;
        ObjectMolecule *lastObj=NULL;
        for(a=cNDummyAtoms;a<I->NAtom;a++) {
          ObjectMolecule *obj=I->Obj[I->Table[a].model];
          if(obj!=lastObj) {
            if(max_n_atom<obj->NAtom)
              max_n_atom = obj->NAtom;
            lastObj = obj;
          }
        }
        zero=Calloc(int,max_n_atom);
        scratch=Alloc(int,max_n_atom);
      }

      for(a=0;a<n_res;a++) {
        res[a].n_acc = 0;
        res[a].n_don = 0;
      }
      hbc = &hbcRec;
      ObjectMoleculeInitHBondCriteria(G,hbc);

      /* use parameters which reflect the spirit of Kabsch and Sander
         ( i.e. long hydrogen-bonds/polar electrostatic interactions ) */

      hbc->maxAngle = 63.0F;
      hbc->maxDistAtMaxAngle = 3.2F;
      hbc->maxDistAtZero = 4.0F;
      hbc->power_a = 1.6F;
      hbc->power_b = 5.0F;
      hbc->cone_dangle = 0.0F; /* 180 deg. */
      if(hbc->maxDistAtMaxAngle!=0.0F) {
        hbc->factor_a = 0.5F/(float)pow(hbc->maxAngle,hbc->power_a);
        hbc->factor_b = 0.5F/(float)pow(hbc->maxAngle,hbc->power_b);
      }
    
      cutoff = hbc->maxDistAtMaxAngle;
      if(cutoff<hbc->maxDistAtZero) {
        cutoff = hbc->maxDistAtZero; 
      }
    
      c=0;
      n1=0;
    
      for(aa=0;aa<I->NAtom;aa++) { /* first, clear flags */
        I->Flag1[aa]=false;
        I->Flag2[aa]=false;
      }

      for(a=0;a<n_res;a++) {
        if(res[a].present) {
          obj0=res[a].obj;
      
          if(obj0) {
            /* map will contain the h-bond backbone nitrogens */
        
            aa=res[a].n;
            at=I->Table[aa].atom;    
            I->Flag2[aa] = a; /* so we can find the atom again... */
        
            if(state<obj0->NCSet) 
              cs=obj0->CSet[state];
            else
              cs=NULL;
            if(cs) {
              if(obj0->DiscreteFlag) {
                if(cs==obj0->DiscreteCSet[at])
                  idx=obj0->DiscreteAtmToIdx[at];
                else
                  idx=-1;
              } else 
                idx=cs->AtmToIdx[at];
              if(idx>=0) {
                copy3f(cs->Coord+(3*idx),I->Vertex+3*aa); /* record coordinate */
                I->Flag1[aa]=true;
#if 0
                printf(" storing donor for %s %d at %8.3f %8.3f %8.3f\n",
                       res[a].obj->AtomInfo[at].resi,idx,
                       I->Vertex[3*aa],I->Vertex[3*aa+1],I->Vertex[3*aa+2]);
#endif
                n1++;
              }
            }


            /* also copy O coordinates for usage below */
          
            aa=res[a].o;
            at=I->Table[aa].atom;    
        
            if(state<obj0->NCSet) 
              cs=obj0->CSet[state];
            else
              cs=NULL;
            if(cs) {
              if(obj0->DiscreteFlag) {
                if(cs==obj0->DiscreteCSet[at])
                  idx=obj0->DiscreteAtmToIdx[at];
                else
                  idx=-1;
              } else 
                idx=cs->AtmToIdx[at];
              if(idx>=0) {
                copy3f(cs->Coord+(3*idx),I->Vertex+3*aa); /* record coordinate */
#if 0
                printf(" storing acceptor for %s %d at %8.3f %8.3f %8.3f\n",
                       res[a].obj->AtomInfo[at].resi,idx,
                       I->Vertex[3*aa],I->Vertex[3*aa+1],I->Vertex[3*aa+2]);
#endif
              }
            }

          }
        }
      }

      if(n1) {
        map=MapNewFlagged(G,-cutoff,I->Vertex,I->NAtom,NULL,I->Flag1);
        if(map) {
          MapSetupExpress(map);
        
          for(a0=0;a0<n_res;a0++) { 
       
            if(res[a0].obj) { 

              /* now iterate through carbonyls */
              obj0 = res[a0].obj;            
              as0 = res[a0].o;
              at0 = I->Table[as0].atom;
            
              v0 = I->Vertex + 3*as0;
              if(MapExclLocus(map,v0,&h,&k,&l)) {
                i=*(MapEStart(map,h,k,l));
                if(i) {
                  as1=map->EList[i++];
                  while(as1>=0) {
                    v1 = I->Vertex+3*as1;

                    if(within3f(v0,v1,cutoff)) {
                    
                      obj1 = I->Obj[I->Table[as1].model];
                      at1 = I->Table[as1].atom;

                      if(obj0==obj1) { /* don't count hbonds between adjacent residues */
                        exclude = SelectorCheckNeighbors(G,5,obj0,at0,at1,
                                                         zero,scratch);
                      } else {
                        exclude = false;
                      }
                    
                      /*                      if(!exclude) {
                        printf("at1 %s %s vs at0 %s %s\n",
                               obj1->AtomInfo[at1].resi,
                               obj1->AtomInfo[at1].name,
                               obj0->AtomInfo[at0].resi,
                               obj0->AtomInfo[at0].name
                               );
                      }
                      */
                      if((!exclude)&&
                         ObjectMoleculeGetCheckHBond(NULL,
                                                     NULL,
                                                     obj1, /* donor first */
                                                     at1,
                                                     state,
                                                     obj0, /* then acceptor */
                                                     at0,
                                                     state,
                                                     hbc)) {
                        
                        /*                        printf(" found hbond between acceptor resi %s and donor resi %s\n",
                                                  res[a0].obj->AtomInfo[at0].resi,
                                                  res[I->Flag2[as1]].obj->AtomInfo[I->Table[as1].atom].resi);*/
                      
                        a1 = I->Flag2[as1]; /* index in SS n_res space */
                      
                        /* store acceptor link */
                      
                        n1 = res[a0].n_acc;
                        if(n1<(cSSMaxHBond-1)) {
                          res[a0].acc[n1] = a1;
                          res[a0].n_acc = n1+1;
                        }
                      
                        /* store donor link */
                      
                        n1 = res[a1].n_don;
                        if(n1<(cSSMaxHBond-1)) {
                          res[a1].don[n1] = a0;
                          res[a1].n_don = n1+1;
                        }
                      }
                    }
                    as1=map->EList[i++];
                  }
                }
              }
            }
          }
        }
        MapFree(map);
      }
      FreeP(zero);
      FreeP(scratch);
    }

    { /* compute phi, psi's */

      SSResi *r;
      int a;

      float helix_psi_delta, helix_phi_delta;
      float strand_psi_delta, strand_phi_delta;

      float helix_psi_target = SettingGet_f(G,NULL,NULL,cSetting_ss_helix_psi_target);
      float helix_psi_include = SettingGet_f(G,NULL,NULL,cSetting_ss_helix_psi_include);
      float helix_psi_exclude = SettingGet_f(G,NULL,NULL,cSetting_ss_helix_psi_exclude);

      float helix_phi_target = SettingGet_f(G,NULL,NULL,cSetting_ss_helix_phi_target);
      float helix_phi_include = SettingGet_f(G,NULL,NULL,cSetting_ss_helix_phi_include);
      float helix_phi_exclude = SettingGet_f(G,NULL,NULL,cSetting_ss_helix_phi_exclude);

      float strand_psi_target = SettingGet_f(G,NULL,NULL,cSetting_ss_strand_psi_target);
      float strand_psi_include = SettingGet_f(G,NULL,NULL,cSetting_ss_strand_psi_include);
      float strand_psi_exclude = SettingGet_f(G,NULL,NULL,cSetting_ss_strand_psi_exclude);

      float strand_phi_target = SettingGet_f(G,NULL,NULL,cSetting_ss_strand_phi_target);
      float strand_phi_include = SettingGet_f(G,NULL,NULL,cSetting_ss_strand_phi_include);
      float strand_phi_exclude = SettingGet_f(G,NULL,NULL,cSetting_ss_strand_phi_exclude);
    
      for(a=0;a<n_res;a++) {
        r = res + a;
        if(r->real&&((r-1)->real)) {
          r->flags = 0;

          if(ObjectMoleculeGetPhiPsi(r->obj,I->Table[r->ca].atom,&r->phi,&r->psi,state)) {
            r->flags |= cSSGotPhiPsi;

            helix_psi_delta = (float)fabs(r->psi - helix_psi_target); 
            strand_psi_delta  = (float)fabs(r->psi - strand_psi_target); 
            helix_phi_delta = (float)fabs(r->phi - helix_phi_target); 
            strand_phi_delta  = (float)fabs(r->phi - strand_phi_target); 
          
            if(helix_psi_delta>180.0F) helix_psi_delta = 360.0F-helix_psi_delta;
            if(strand_psi_delta>180.0F) strand_psi_delta = 360.0F-strand_psi_delta;
            if(helix_phi_delta>180.0F) helix_phi_delta = 360.0F-helix_phi_delta;
            if(strand_phi_delta>180.0F) strand_phi_delta = 360.0F-strand_phi_delta;

            /* printf("helix %d strand %d\n",helix_delta,strand_delta);*/
               
            if((helix_psi_delta>helix_psi_exclude)||
               (helix_phi_delta>helix_phi_exclude)) {
              r->flags |= cSSPhiPsiNotHelix;
            } else if((helix_psi_delta<helix_psi_include) &&
                      (helix_phi_delta<helix_phi_include)) {
              r->flags |= cSSPhiPsiHelix;
            }
          
            if((strand_psi_delta>strand_psi_exclude)||
               (strand_phi_delta>strand_phi_exclude)) {
              r->flags |= cSSPhiPsiNotStrand;
            } else if((strand_psi_delta<strand_psi_include)&&
                      (strand_phi_delta<strand_phi_include)) {
              r->flags |= cSSPhiPsiStrand;
            }
          }
        }
      }
    }

    /* by default, tentatively assign everything as loop */
    
    { 
      int a;
      for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
        if(res[a].present)
          res[a].ss = 'L';
      }
    }


    { 
      SSResi *r,*r2;
      int a,b,c;

      for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
        r = res + a;
        if(r->real) {

          /* look for tell-tale i+3,4,5 hydrogen bonds for helix  */
                
          /* is residue an acceptor for i+3,4,5 residue? */
          for(b=0;b<r->n_acc;b++) {
            r->flags |= 
              ( (r->acc[b]==(a+3)) ? cSSHelix3HBond : 0 ) |
              ( (r->acc[b]==(a+4)) ? cSSHelix4HBond : 0 ) |
              ( (r->acc[b]==(a+5)) ? cSSHelix5HBond : 0 );

          }

          /* is residue a donor for i-3,4,5 residue */
          for(b=0;b<r->n_don;b++) {
            r->flags |= 
              ( (r->don[b]==(a-3)) ? cSSHelix3HBond : 0 ) |
              ( (r->don[b]==(a-4)) ? cSSHelix4HBond : 0 ) |
              ( (r->don[b]==(a-5)) ? cSSHelix5HBond : 0 );

          }

          /*        if(r->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) {
                    printf("HelixHB %s \n",
                    r->obj->AtomInfo[I->Table[r->ca].atom].resi);
                    }
          */


          /* look for double h-bonded antiparallel beta sheet pairs:
           * 
           *  \ /\ /
           *   N  C
           *   #  O
           *   O  #
           *   C  N
           *  / \/ \
           *
           */
        
          for(b=0;b<r->n_acc;b++) { /* iterate through acceptors */
            r2 = (res + r->acc[b]);
            if(r2->real) {
              for(c=0;c<r2->n_acc;c++) {
                if(r2->acc[c] == a) { /* found a pair */
                  r->flags |= cSSAntiStrandDoubleHB;
                  r2->flags |= cSSAntiStrandDoubleHB;

                  /*                printf("anti double %s to %s\n",
                                    r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                                    r2->obj->AtomInfo[I->Table[r2->ca].atom].resi);*/
                
                }
              }
            }
          }

          /* look for antiparallel beta buldges
           * 
           *     CCNC
           *  \ / O  \ /
           *   N      C
           *   #      O
           *    O    #
           *     C  N
           *    / \/ \
           *
           */
        
          for(b=0;b<r->n_acc;b++) { /* iterate through acceptors */
            r2 = (res + r->acc[b])+1; /* go forward 1 */
            if(r2->real) {
              for(c=0;c<r2->n_acc;c++) {
                if(r2->acc[c] == a) { /* found a buldge */
                  r->flags      |= cSSAntiStrandDoubleHB;
                  r2->flags     |= cSSAntiStrandBuldgeHB;
                  (r2-1)->flags |= cSSAntiStrandBuldgeHB;

                  /*                printf("anti BULDGE %s to %s %s\n",
                                    r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                                    r2->obj->AtomInfo[I->Table[r2->ca].atom].resi,
                                    r2->obj->AtomInfo[I->Table[(r2-1)->ca].atom].resi);*/

                }
              }
            }
          }
        
          /* look for antiparallel beta sheet ladders (single or double)
           *
           *        O
           *     N  C
           *  \ / \/ \ /
           *   C      N
           *   O      #
           *   #      O
           *   N      C
           *  / \ /\ / \
           *     C  N
           *     O
           */
        
          if( (r+1)->real && (r+2)->real ) {
          
            for(b=0;b<r->n_acc;b++) { /* iterate through acceptors */
              r2 = (res + r->acc[b])-2; /* go back 2 */
              if(r2->real) {
              
                for(c=0;c<r2->n_acc;c++) {
                
                  if(r2->acc[c] == a+2) { /* found a ladder */
                  
                    (r    )->flags |= cSSAntiStrandSingleHB;
                    (r  +1)->flags |= cSSAntiStrandSkip;
                    (r  +2)->flags |= cSSAntiStrandSingleHB;
                  
                    (r2   )->flags |= cSSAntiStrandSingleHB;
                    (r2 +1)->flags |= cSSAntiStrandSkip;
                    (r2 +2)->flags |= cSSAntiStrandSingleHB;
                  
                    /*                  printf("anti ladder %s %s to %s %s\n",
                                        r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                                        r->obj->AtomInfo[I->Table[(r+2)->ca].atom].resi,
                                        r2->obj->AtomInfo[I->Table[r2->ca].atom].resi,
                                        r2->obj->AtomInfo[I->Table[(r2+2)->ca].atom].resi);*/
                  }
                }
              }
            }
          }


          /* look for parallel beta sheet ladders 
           *

           *    \ /\ /
           *     C  N
           *    O    #
           *   #      O
           *   N      C
           *  / \ /\ / \
           *     C  N
           *     O
           */
        
          if( (r+1)->real && (r+2)->real ) {
          
            for(b=0;b<r->n_acc;b++) { /* iterate through acceptors */
              r2 = (res + r->acc[b]);
              if(r2->real) {
              
                for(c=0;c<r2->n_acc;c++) {
                
                  if(r2->acc[c] == a+2) { /* found a ladder */
                  
                    (r    )->flags |= cSSParaStrandSingleHB;
                    (r  +1)->flags |= cSSParaStrandSkip;
                    (r  +2)->flags |= cSSParaStrandSingleHB;
                  
                    (r2   )->flags |= cSSParaStrandDoubleHB;
                  
                    /*                                    printf("parallel ladder %s %s to %s \n",
                                                          r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                                                          r->obj->AtomInfo[I->Table[(r+2)->ca].atom].resi,
                                                          r2->obj->AtomInfo[I->Table[r2->ca].atom].resi);*/
                  }
                }
              }
            }
          }
        }
      }
    }

  
    {
      int a;
      SSResi *r;
      /* convert flags to assignments */
    
      /* HELICES FIRST */

      for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
        r = res + a;
      
        if(r->real) {
          /* clean internal helical residues are easy to find using H-bonds */

          if(((r-1)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r  )->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r+1)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond))) {
            if(!(r->flags & (cSSPhiPsiNotHelix))) {
              r->ss = 'H';
            }
          }

          /*
            if(((r-1)->flags & (cSSHelix3HBond )) &&
            ((r  )->flags & (cSSHelix3HBond )) &&
            ((r+1)->flags & (cSSHelix3HBond ))) {
            if(!(r->flags & (cSSPhiPsiNotHelix))) {
            r->ss = 'H';
            }
            }

            if(((r-1)->flags & (cSSHelix4HBond)) &&
            ((r  )->flags & (cSSHelix4HBond)) &&
            ((r+1)->flags & (cSSHelix4HBond))) {
            if(!(r->flags & (cSSPhiPsiNotHelix))) {
            r->ss = 'H';
            }
            }

            if(((r-1)->flags & (cSSHelix5HBond)) &&
            ((r  )->flags & (cSSHelix5HBond)) &&
            ((r+1)->flags & (cSSHelix5HBond))) {
            if(!(r->flags & (cSSPhiPsiNotHelix))) {
            r->ss = 'H';
            }
            }
          */

        }
      }


      for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
        r = res + a;

        if(r->real) {

          /* occasionally they'll be one whacked out residue missing h-bonds... 
             in an otherwise good segment */

          if(((r-2)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r-1)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r-1)->flags & (cSSPhiPsiHelix                                  )) &&
             ((r  )->flags & (cSSPhiPsiHelix                                  )) &&
             ((r+1)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r+1)->flags & (cSSPhiPsiHelix                                  )) &&
             ((r+2)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) 
             ) {
            r->ss = 'h';
          }
        }
      }


      for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
        r = res + a;
        if(r->real) {
          if(r->ss=='h') {
            r->flags |= (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond);
            r->ss = 'H';
          }
        }
      }

      for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
        r = res + a;

        if(r->real) {

          /* deciding where the helix ends is trickier -- here we use helix geometry */

          if(((r  )->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r  )->flags & (cSSPhiPsiHelix                                  )) &&
             ((r+1)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r+1)->flags & (cSSPhiPsiHelix                                  )) &&
             ((r+2)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r+2)->flags & (cSSPhiPsiHelix                                  )) &&
             ((r+1)->ss=='H')
             ) {
            r->ss = 'H';
          }

          if(((r  )->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r  )->flags & (cSSPhiPsiHelix                                  )) &&
             ((r-1)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r-1)->flags & (cSSPhiPsiHelix                                  )) &&
             ((r-2)->flags & (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
             ((r-2)->flags & (cSSPhiPsiHelix                                  )) &&
             ((r-1)->ss=='H')
             ) {
            r->ss = 'H';
          }

        }
      }
    
      /* THEN SHEETS/STRANDS */
    
      for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
        r = res + a;
        if(r->real) {

          /* Antiparallel Sheets */

          if(((r  )->flags & (cSSAntiStrandDoubleHB))&&
             (!((r->flags & (cSSPhiPsiNotStrand))))) {
            (r  )->ss = 'S';
          }

          if(((r  )->flags & (cSSAntiStrandBuldgeHB))&&  /* no strand geometry filtering for buldges..*/
             ((r+1)->flags & (cSSAntiStrandBuldgeHB))) {
            (r  )->ss = 'S';
            (r+1)->ss = 'S';
          }

          if(((r-1)->flags & (cSSAntiStrandDoubleHB)) &&
             ((r  )->flags & (cSSAntiStrandSkip))     &&
             (!(((r  )->flags & (cSSPhiPsiNotStrand)))) &&
             ((r+1)->flags & (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB))) {

            (r  )->ss = 'S';
          }

          if(((r-1)->flags & (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB)) &&
             ((r  )->flags & (cSSAntiStrandSkip))     &&
             (!(((r  )->flags & (cSSPhiPsiNotStrand)))) &&
             ((r+1)->flags & (cSSAntiStrandDoubleHB))) {
            (r  )->ss = 'S';
          }

          /* include open "ladders" if PHIPSI geometry supports assignment */

          if(((r-1)->flags & (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB)) &&
             ((r-1)->flags & (cSSPhiPsiStrand)) &&
             (!(((r-1)->flags & (cSSPhiPsiNotStrand)))) &&           
             ((r  )->flags & (cSSPhiPsiStrand)) &&
             (!(((r-1)->flags & (cSSPhiPsiNotStrand)))) &&           
             ((r+1)->flags & (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB)) &&
             ((r+1)->flags & (cSSPhiPsiStrand))) {
          
            (r-1)->ss = 'S';
            (r  )->ss = 'S';
            (r+1)->ss = 'S';
          }

          /* Parallel Sheets */

          if(((r  )->flags & (cSSParaStrandDoubleHB))&&
             (!(((r  )->flags & (cSSPhiPsiNotStrand))))) {
            (r  )->ss = 'S';
          }
        
          if(((r-1)->flags & (cSSParaStrandDoubleHB)) &&
             ((r  )->flags & (cSSParaStrandSkip))     &&
             (!(((r  )->flags & (cSSPhiPsiNotStrand))))&&
             ((r+1)->flags & (cSSParaStrandSingleHB | cSSParaStrandDoubleHB))) {

            (r  )->ss = 'S';
          }

          if(((r-1)->flags & (cSSParaStrandSingleHB | cSSParaStrandDoubleHB)) &&
             ((r  )->flags & (cSSParaStrandSkip))     &&
             (!(((r  )->flags & (cSSPhiPsiNotStrand))))&&
             ((r+1)->flags & (cSSParaStrandDoubleHB))) {
            (r  )->ss = 'S';
          }

          /* include open "ladders" if PHIPSI geometry supports assignment */

          if(((r-1)->flags & (cSSParaStrandSingleHB | cSSParaStrandDoubleHB)) &&
             ((r-1)->flags & (cSSPhiPsiStrand)) &&
             ((r  )->flags & (cSSParaStrandSkip))     &&
             ((r  )->flags & (cSSPhiPsiStrand)) &&
             ((r+1)->flags & (cSSParaStrandSingleHB | cSSParaStrandDoubleHB)) &&
             ((r+1)->flags & (cSSPhiPsiStrand))) {
            
            (r-1)->ss = 'S';
            (r  )->ss = 'S';
            (r+1)->ss = 'S';
          
          }
        }
      }
    }
    
    {
      int a,b;
      SSResi *r,*r2;
      int repeat = true;
      int found;
    
      while(repeat) {
        repeat = false;
            
        for(a=cSSBreakSize;a<(n_res-cSSBreakSize);a++) {
          r = res + a;
          if(r->real) {
          
            /* make sure we don't have any 2-residue segments */
          
            if((r->ss == 'S')&&((r+1)->ss == 'S') &&
               ( ((r-1)->ss!='S') && ((r+2)->ss!='S') )) {
              r->ss = 'L';
              (r+1)->ss = 'L';
              repeat=true;
            }
            if((r->ss == 'H')&&((r+1)->ss == 'H') &&
               ( ((r-1)->ss!='H') && ((r+2)->ss!='H') )) {
              r->ss = 'L';
              (r+1)->ss = 'L';
              repeat=true;
            }

            /* make sure we don't have any 1-residue segments */
          
            if((r->ss == 'S') && ( ((r-1)->ss!='S') && ((r+1)->ss!='S') )) {
              r->ss = 'L';
              repeat=true;
            }
            if((r->ss == 'H') && ( ((r-1)->ss!='H') && ((r+1)->ss!='H') )) {
              r->ss = 'L';
              repeat=true;
            }
          
            /* double-check to make sure every terminal strand residue 
               that should have a partner has one */
          
            if((r->ss == 'S') && (((r-1)->ss!='S') || ((r+1)->ss!='S') )) {
            
              found = false;
            
              for(b=0;b<r->n_acc;b++) {
                r2 = res+r->acc[b];
                if(r2->ss == r->ss) {
                  found=true;
                  break;
                }
              }
            
              if(!found) {
                for(b=0;b<r->n_don;b++) {
                  r2 = res+r->don[b];
                  if(r2->ss == r->ss) {
                    found=true;
                    break;
                  }
                }
              }
            
              if(!found) { /* allow these strand "skip" residues to persist if a neighbor has hydrogen bonds */
                if(r->flags & (cSSAntiStrandSkip | cSSParaStrandSkip)) {
                
                  if((r+1)->ss == r->ss)
                    for(b=0;b<(r+1)->n_acc;b++) {
                      r2 = res+(r+1)->acc[b];
                      if(r2->ss == r->ss) {
                        found=true;
                        break;
                      }
                    }
                
                  if(!found) {
                    if((r-1)->ss == r->ss) {
                      for(b=0;b<(r-1)->n_don;b++) {
                        r2 = res+(r-1)->don[b];
                        if(r2->ss == r->ss) {
                          found=true;
                          break;
                        }
                      }
                    }
                  }
                }
              }
            
              if(!found) {
                r->ss = 'L';
                repeat=true;
              }
            }
          }
        }
      }
    }

    { 
      int a;
      for(a=0;a<n_res;a++) { /* now apply consensus or union behavior, if appropriate*/
        if(res[a].present) {
          if(res[a].ss_save) {
            if(res[a].ss!=res[a].ss_save) {
              if(consensus) {
                res[a].ss = res[a].ss_save = 'L';
              } else if(res[a].ss=='L')
                res[a].ss = res[a].ss_save;
            }
          }
          res[a].ss_save = res[a].ss;
        }
      }
    }

    { 
      int a,aa;
      ObjectMolecule *obj=NULL,*last_obj = NULL;
      AtomInfoType *ai;
      int changed_flag = false;

      for(a=0;a<n_res;a++) {
        if(res[a].present&&(!res[a].preserve)) {

          aa = res[a].ca;
          obj=I->Obj[I->Table[aa].model];

          if(obj!=last_obj) {
            if(changed_flag&&last_obj) {
              ObjectMoleculeInvalidate(last_obj,cRepCartoon,cRepInvRep,-1);
              SceneChanged(G);
              changed_flag=false;
            }
            last_obj=obj;
          }
          ai = obj->AtomInfo + I->Table[aa].atom;
        
          if(SelectorIsMember(G,ai->selEntry,target)) {
            ai->ssType[0] = res[a].ss;
            ai->cartoon = 0; /* switch back to auto */
            ai->ssType[1] = 0;
            changed_flag=true;
          }
        }
      }

      if(changed_flag&&last_obj) {
        ObjectMoleculeInvalidate(last_obj,cRepCartoon,cRepInvRep,-1);
        SceneChanged(G);
        changed_flag=false;
      }
    }
    if(first_last_only&&(state==state_start))
      state = state_stop-2;
  }
       
  VLAFreeP(res);
  return 1;
}

PyObject *SelectorColorectionGet(PyMOLGlobals *G,char *prefix)

{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  register CSelector *I=G->Selector;
  PyObject *result = NULL;
  int n_used=0;
  ColorectionRec *used = NULL,tmp;
  ov_size a,b,n;
  int sele;
  int found;
  int m;
  int color;
  AtomInfoType *ai;
  used=VLAlloc(ColorectionRec,1000);
  
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    ai = I->Obj[I->Table[a].model]->AtomInfo+I->Table[a].atom;
    color = ai->color;
    found = false;
    for(b=0;b<n_used;b++) {
      if(used[b].color==color) {
        tmp=used[0]; /* optimize to minimize N^2 effects */
        used[0]=used[b];
        used[b]=tmp;
        found=true;
        break;
      }
    }
    if(!found) {
      VLACheck(used,ColorectionRec,n_used);
      used[n_used]=used[0];
      used[0].color=color;
      n_used++;
    }
  }
  for(a=0;a<n_used;a++) {
    /* create selections */

    n=I->NActive;
    VLACheck(I->Name,SelectorWordType,n+1);
    VLACheck(I->Info,SelectionInfoRec,n+1);
    sele = I->NSelection++;
    used[a].sele = sele;
    sprintf(I->Name[n],cColorectionFormat,prefix,used[a].color);
    I->Name[n+1][0]=0;
    SelectorAddName(G,n);
    SelectionInfoInit(I->Info + n);
    I->Info[n].ID = sele;
    I->NActive++;
  }

  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    ai = I->Obj[I->Table[a].model]->AtomInfo+I->Table[a].atom;
    color = ai->color;
    for(b=0;b<n_used;b++) {
      if(used[b].color==color) {
        tmp=used[0]; /* optimize to minimize N^2 effects */
        used[0]=used[b];
        used[b]=tmp;

        /* add selection onto atom */
        if(I->FreeMember>0) {
          m=I->FreeMember;
          I->FreeMember=I->Member[m].next;
        } else {
          I->NMember++;
          m=I->NMember;
          VLACheck(I->Member,MemberType,m);
        }
        I->Member[m].selection = used[0].sele;
        I->Member[m].tag = 1;
        I->Member[m].next = ai->selEntry;
        ai->selEntry = m;
        break;
      }
    }
  }

  VLASize(used,ColorectionRec,n_used*2);
  result = PConvIntVLAToPyList((int*)used);
  VLAFreeP(used);
  return(result);
#endif
}

int SelectorColorectionApply(PyMOLGlobals *G,PyObject *list,char *prefix)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  register CSelector *I=G->Selector;
  int ok=true;
  ColorectionRec *used=NULL;
  ov_size n_used=0;
  int a,b;
  AtomInfoType *ai;
  ObjectMolecule *obj,*last=NULL;
  SelectorWordType name;

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) n_used = PyList_Size(list)/2;
  if(ok) ok=((used=VLAlloc(ColorectionRec,n_used))!=NULL);
  if(ok) ok=PConvPyListToIntArrayInPlace(list,(int*)used,n_used*2);
  if(ok) {

    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);

    for(b=0;b<n_used;b++) { /* update selection indices */
      sprintf(name,cColorectionFormat,prefix,used[b].color);      
      used[b].sele = SelectorIndexByName(G,name);
    }
    
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      obj = I->Obj[I->Table[a].model]; 
      ai = obj->AtomInfo+I->Table[a].atom;
      
      for(b=0;b<n_used;b++) {       
        if(SelectorIsMember(G,ai->selEntry,used[b].sele)) {
          ai->color = used[b].color;
          if(obj!=last) {
            ObjectMoleculeInvalidate(obj,cRepAll,cRepInvColor,-1);
            last = obj;
          }
          break;
        } 
      }
    }
  }
  VLAFreeP(used);
  return(ok);
#endif
}

int SelectorColorectionSetName(PyMOLGlobals *G,PyObject *list,char *prefix,char *new_prefix)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  ColorectionRec *used=NULL;
  ov_size n_used=0;
  ov_size b;
  SelectorWordType name;
  SelectorWordType new_name;

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) n_used = PyList_Size(list)/2;
  if(ok) ok=((used=VLAlloc(ColorectionRec,n_used))!=NULL);
  if(ok) ok=PConvPyListToIntArrayInPlace(list,(int*)used,n_used*2);
  if(ok) {
    for(b=0;b<n_used;b++) { /* update selection indices */
      sprintf(name,cColorectionFormat,prefix,used[b].color);      
      sprintf(new_name,cColorectionFormat,new_prefix,used[b].color);      
      SelectorSetName(G,new_name,name);
    }
  }
  VLAFreeP(used);
  return(ok);
#endif

}

int SelectorColorectionFree(PyMOLGlobals *G,PyObject *list,char *prefix)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  ColorectionRec *used=NULL;
  ov_size n_used=0;
  ov_size b;
  SelectorWordType name;

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) n_used = PyList_Size(list)/2;
  if(ok) ok=((used=VLAlloc(ColorectionRec,n_used))!=NULL);
  if(ok) ok=PConvPyListToIntArrayInPlace(list,(int*)used,n_used*2);
  if(ok) {

    for(b=0;b<n_used;b++) { /* update selection indices */
      sprintf(name,cColorectionFormat,prefix,used[b].color);      
      used[b].sele = SelectorIndexByName(G,name);
    }

    for(b=0;b<n_used;b++) {
      SelectorDeleteIndex(G,used[b].sele);
    }
  }
  VLAFreeP(used);
  return(ok);
#endif

}

PyObject *SelectorSecretsAsPyList(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  register CSelector *I=G->Selector;
  int n_secret;
  int a;
  PyObject *result,*list;
  
  n_secret=0;
  for(a=0;a<I->NActive;a++) {
    if((I->Name[a][0]=='_')&&
       (I->Name[a][1]=='!'))
      n_secret++;
  }    
  result = PyList_New(n_secret);
  n_secret=0;
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  for(a=0;a<I->NActive;a++) {
    if((I->Name[a][0]=='_')&&
       (I->Name[a][1]=='!')) {
      list = PyList_New(2);
      PyList_SetItem(list,0,PyString_FromString(I->Name[a]));
      PyList_SetItem(list,1,SelectorAsPyList(G,I->Info[a].ID));
      PyList_SetItem(result,n_secret,list);
      n_secret++;
    }
  }    
  return(result);
#endif
}

int SelectorSecretsFromPyList(PyMOLGlobals *G,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok=true;
  ov_size n_secret=0;
  ov_size a;
  PyObject *entry=NULL;
  SelectorWordType name;
  ov_size ll=0;
  if(ok) ok = (list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) n_secret = PyList_Size(list);
  if(ok) {
    for(a=0;a<n_secret;a++) {
      if(ok) entry = PyList_GetItem(list,a);
      if(ok) ok = (entry!=NULL);
      if(ok) ok = PyList_Check(entry);
      if(ok) ll = PyList_Size(entry);
      if(ok&(ll>1)) {
        if(ok) ok = PConvPyStrToStr(PyList_GetItem(entry,0),name,sizeof(SelectorWordType));
        if(ok) ok = SelectorFromPyList(G,name,PyList_GetItem(entry,1));
      }
      if(!ok) break;
    }
  }
  return(ok);
#endif
}

typedef struct {
  int atom;
  int tag;
} SelAtomTag;

PyObject *SelectorAsPyList(PyMOLGlobals *G,int sele1)
{ /* assumes SelectorUpdateTable has been called */
#ifdef _PYMOL_NOPY
  return NULL;
#else
  register CSelector *I=G->Selector;
  int a,b;
  int at;
  int s;
  SelAtomTag **vla_list = NULL;
  int n_obj = 0;
  int n_idx = 0;
  int cur = -1;
  ObjectMolecule **obj_list = NULL;
  ObjectMolecule *obj,*cur_obj = NULL;
  PyObject *result = NULL;
  PyObject *obj_pyobj;
  PyObject *idx_pyobj;
  PyObject *tag_pyobj;

  vla_list = VLAMalloc(10,sizeof(SelAtomTag*),5,true);
  obj_list = VLAlloc(ObjectMolecule*,10);

  n_idx = 0;
  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    int tag;
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if( (tag = SelectorIsMember(G,s,sele1)) ) {
      if(cur_obj!=obj) {
        if(n_idx) {
          VLASize(vla_list[cur],SelAtomTag,n_idx);
        }
        cur++;
        VLACheck(vla_list,SelAtomTag*,n_obj);
        vla_list[cur] = VLAlloc(SelAtomTag,1000);
        VLACheck(obj_list,ObjectMolecule*,n_obj);
        obj_list[cur] = obj;
        cur_obj = obj;
        n_obj++;
        n_idx=0;
      }
      VLACheck(vla_list[cur],SelAtomTag,n_idx);
      vla_list[cur][n_idx].atom = at;
      vla_list[cur][n_idx].tag = tag;
      n_idx++;
    }
  }
  if(cur_obj) {
    if(n_idx) {
      VLASize(vla_list[cur],SelAtomTag,n_idx);
    }
  }
  if(n_obj) {
    result = PyList_New(n_obj);
    for(a=0;a<n_obj;a++) {
      obj_pyobj= PyList_New(3);
      n_idx = VLAGetSize(vla_list[a]);
      idx_pyobj = PyList_New(n_idx);
      tag_pyobj = PyList_New(n_idx);
      for(b=0;b<n_idx;b++) {
        PyList_SetItem(idx_pyobj,b,PyInt_FromLong(vla_list[a][b].atom));
        PyList_SetItem(tag_pyobj,b,PyInt_FromLong(vla_list[a][b].tag));
      }
      VLAFreeP(vla_list[a]);
      PyList_SetItem(obj_pyobj,0,PyString_FromString(obj_list[a]->Obj.Name));
      PyList_SetItem(obj_pyobj,1,idx_pyobj);
      PyList_SetItem(obj_pyobj,2,tag_pyobj);
      PyList_SetItem(result,a,obj_pyobj);
    }
  } else {
    result = PyList_New(0);
  }
  VLAFreeP(vla_list);
  VLAFreeP(obj_list);
  return(result);
#endif
}


int SelectorFromPyList(PyMOLGlobals *G,char *name,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  register CSelector *I=G->Selector;
  ov_size a,b;
  ov_diff n;
  int m,sele;
  ov_size ll;
  PyObject *obj_list=NULL;
  PyObject *idx_list=NULL,*tag_list;
  ov_size n_obj=0,n_idx=0;
  int idx,tag;
  char *oname;
  ObjectMolecule *obj;
  int singleAtomFlag = true;
  int singleObjectFlag = true;
  ObjectMolecule *singleObject = NULL;
  int singleAtom = -1;
  int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);

  AtomInfoType *ai;
  if(ok) ok=PyList_Check(list);
  if(ok) n_obj = PyList_Size(list);
    
  n=SelectGetNameOffset(G,name,999,ignore_case); /* already exist? */
  if(n>=0) /* get rid of existing selection*/ {
    SelectorDelete(G,I->Name[n]);
  }

  n=I->NActive;
  VLACheck(I->Name,SelectorWordType,n+1);
  VLACheck(I->Info,SelectionInfoRec,n+1);
  strcpy(I->Name[n],name);
  I->Name[n+1][0]=0;
  SelectorAddName(G,n);
  sele = I->NSelection++;
  SelectionInfoInit(I->Info + n);
  I->Info[n].ID = sele;
  I->NActive++;
  if(ok) {
    for(a=0;a<n_obj;a++) {
      ll = 0;
      if(ok) obj_list = PyList_GetItem(list,a);
      if(ok) ok = PyList_Check(obj_list);
      if(ok) ll = PyList_Size(obj_list);
      if(ok) ok = PConvPyStrToStrPtr(PyList_GetItem(obj_list,0),&oname);
      obj=NULL;
      if(ok) obj = ExecutiveFindObjectMoleculeByName(G,oname);
      if(ok&&obj) {
        if(ok) idx_list = PyList_GetItem(obj_list,1);
        if(ll>2)
          tag_list = PyList_GetItem(obj_list,2);
        else
          tag_list = NULL;
        if(ok) ok = PyList_Check(idx_list);
        if(ok) n_idx = PyList_Size(idx_list);
        for(b=0;b<n_idx;b++) {
          if(ok) ok = PConvPyIntToInt(PyList_GetItem(idx_list,b),&idx);
          if(tag_list) 
            PConvPyIntToInt(PyList_GetItem(tag_list,b),&tag);
          else
            tag = 1;
          if(ok&&(idx<obj->NAtom)) {
            ai=obj->AtomInfo+idx;
            if(I->FreeMember>0) {
              m=I->FreeMember;
              I->FreeMember=I->Member[m].next;
            } else {
              I->NMember++;
              m=I->NMember;
              VLACheck(I->Member,MemberType,m);
            }
            I->Member[m].selection = sele;
            I->Member[m].tag = tag; 
            I->Member[m].next = ai->selEntry;
            ai->selEntry = m;

            /* take note of selections which are one atom/one object */
            if(singleObjectFlag) {
              if(singleObject) {
                if( obj!=singleObject) {
                  singleObjectFlag = false;
                }
              } else {
                singleObject = obj;
              }
            }
            
            if(singleAtomFlag) {
              if(singleAtom>=0) {
                if(idx!=singleAtom) {
                  singleAtomFlag = false;
                }
              } else {
                singleAtom = idx;
              }
            }
          }
        }
      }
    }
    { /* make note of single atom/object selections */
      SelectionInfoRec *info = I->Info + (I->NActive-1);
      if( singleObjectFlag && singleObject ) {
        info->justOneObjectFlag = true;
        info->theOneObject = singleObject;
        if( singleAtomFlag && (singleAtom>=0) ) {
          info->justOneAtomFlag = true;
          info->theOneAtom = singleAtom;
        }
      }
    }
  }

  return(ok);
#endif

}

int SelectorVdwFit(PyMOLGlobals *G,int sele1,int state1, int sele2,int state2, float buffer, int quiet)
{
  int ok=true;
  register CSelector *I=G->Selector;
  int *vla=NULL;
  int c;
  float sumVDW=0.0,dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at1,at2;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj1,*obj2;
  int idx1,idx2;
  float *adj = NULL;
  int a;

  if(state1<0) state1=0;
  if(state2<0) state2=0;

  if(state1!=state2) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,state1,-1);
  }
  
  c=SelectorGetInterstateVLA(G,sele1,state1,sele2,state2,2*MAX_VDW+buffer,&vla);
  if(c) {
    adj = Calloc(float,2*c);
    for(a=0;a<c;a++) {
      a1=vla[a*2];
      a2=vla[a*2+1];
      
      at1=I->Table[a1].atom;
      at2=I->Table[a2].atom;
      
      obj1=I->Obj[I->Table[a1].model];
      obj2=I->Obj[I->Table[a2].model];
      
      if((state1<obj1->NCSet)&&(state2<obj2->NCSet)) {
        cs1=obj1->CSet[state1];
        cs2=obj2->CSet[state2];
        if(cs1&&cs2) { /* should always be true */
          
          ai1=obj1->AtomInfo+at1;
          ai2=obj2->AtomInfo+at2;
          
          idx1=cs1->AtmToIdx[at1]; /* these are also pre-validated */
          idx2=cs2->AtmToIdx[at2];
          
          sumVDW=ai1->vdw+ai2->vdw;
          dist=(float)diff3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2);
          
          if(dist<(sumVDW+buffer)) {
            float shift = (dist-(sumVDW+buffer))/2.0F;
            adj[2*a] = ai1->vdw + shift;
            adj[2*a+1] = ai2->vdw + shift;
          } else {
            adj[2*a] = ai1->vdw;
            adj[2*a+1] = ai2->vdw;
          }
          
        }
      }
    }
    
    for(a=0;a<c;a++) {
      a1=vla[a*2];
      a2=vla[a*2+1];
      
      at1=I->Table[a1].atom;
      at2=I->Table[a2].atom;
      
      obj1=I->Obj[I->Table[a1].model];
      obj2=I->Obj[I->Table[a2].model];
      
      if((state1<obj1->NCSet)&&(state2<obj2->NCSet)) {
        cs1=obj1->CSet[state1];
        cs2=obj2->CSet[state2];
        if(cs1&&cs2) { /* should always be true */
          
          ai1=obj1->AtomInfo+at1;
          ai2=obj2->AtomInfo+at2;
          
          if(adj[2*a] < ai1->vdw) {
            ai1->vdw = adj[2*a];
          } 
          
          if(adj[2*a+1] < ai2->vdw) {
            ai2->vdw = adj[2*a+1];
          }
          
        }
      }
    }
  }
 
  VLAFreeP(vla);
  FreeP(adj);
  return ok;
}

/*========================================================================*/

int SelectorGetPairIndices(PyMOLGlobals *G,int sele1,int state1,int sele2,int state2,
                           int mode,float cutoff, float h_angle,
                           int **indexVLA, ObjectMolecule ***objVLA)
{
  register CSelector *I=G->Selector;
  int *vla=NULL;
  int c;
  float dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at1,at2;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj1,*obj2;
  int idx1,idx2;
  int a;
  int dist_cnt = 0;
  float dir[3];
  float v1[3],v2[3];
  int flag;
  float angle_cutoff=0.0;

  if(mode==1) {
    angle_cutoff = (float)cos(PI*h_angle/180.8);
  }

  if(state1<0) state1=0;
  if(state2<0) state2=0;

  if(state1!=state2) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,state1,-1);
  }
  if(cutoff<0) cutoff = 1000.0;
  c=SelectorGetInterstateVLA(G,sele1,state1,sele2,state2,cutoff,&vla);
  (*indexVLA)=VLAlloc(int,1000);
  (*objVLA)=VLAlloc(ObjectMolecule*,1000);


  for(a=0;a<c;a++) {
    a1=vla[a*2];
    a2=vla[a*2+1];

    if(a1!=a2) { 
      at1=I->Table[a1].atom;
      at2=I->Table[a2].atom;
      
      obj1=I->Obj[I->Table[a1].model];
      obj2=I->Obj[I->Table[a2].model];

      if(state1<obj1->NCSet&&state2<obj2->NCSet) {
        cs1=obj1->CSet[state1];
        cs2=obj2->CSet[state2];
        if(cs1&&cs2) { 
    
          ai1=obj1->AtomInfo+at1;
          ai2=obj2->AtomInfo+at2;

          if(obj1->DiscreteFlag) {
            if(cs1==obj1->DiscreteCSet[at1]) {
              idx1=obj1->DiscreteAtmToIdx[at1];
            } else {
              idx1=-1;
            }
          } else {
            idx1=cs1->AtmToIdx[at1];
          }
          
          if(obj2->DiscreteFlag) {
            if(cs2==obj2->DiscreteCSet[at2]) {
              idx2=obj2->DiscreteAtmToIdx[at2];
            } else {
              idx2=-1;
            }
            
          } else {
            idx2=cs2->AtmToIdx[at2];
          }
            
          if((idx1>=0)&&(idx2>=0)) {
            subtract3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2,dir);
            dist=(float)length3f(dir);
            if(dist>R_SMALL4) {
              float dist_1 = 1.0F/dist;
              scale3f(dir,dist_1,dir);
            }
            if(dist<cutoff) {
              if(mode==1) { /* coarse hydrogen bonding assessment */
                flag=false;
                if(ObjectMoleculeGetAvgHBondVector(obj1,at1,state1,v1,NULL)>0.3)
                  if(dot_product3f(v1,dir)<-angle_cutoff) 
                    flag=true;
                if(ObjectMoleculeGetAvgHBondVector(obj2,at2,state2,v2,NULL)>0.3)
                  if(dot_product3f(v2,dir)>angle_cutoff)
                    flag=true;
              } else 
                flag=true;

              if(flag) {
                VLACheck((*objVLA),ObjectMolecule*,dist_cnt+1);
                VLACheck((*indexVLA),int,dist_cnt+1);
                (*objVLA)[dist_cnt]=obj1;
                (*indexVLA)[dist_cnt]=at1;
                dist_cnt++;
                (*objVLA)[dist_cnt]=obj2;
                (*indexVLA)[dist_cnt]=at2;              
                dist_cnt++;
              }
            }
          }
        }
      }
    }
  }
  
  VLASize((*objVLA),ObjectMolecule*,dist_cnt);
  VLASize((*indexVLA),int,dist_cnt);
  VLAFreeP(vla);
  dist_cnt = dist_cnt / 2;
  return(dist_cnt);
}

/*========================================================================*/
int  SelectorCreateAlignments(PyMOLGlobals *G,
                              int *pair,int sele1,int *vla1,int sele2,
                              int *vla2,char *name1,char *name2,
                              int identical,int atomic_input)
{
  register CSelector *I=G->Selector;
  int *flag1=NULL,*flag2=NULL;
  int *p;
  int i,np;
  int cnt; 
  int mod1,mod2; /* model indexes */
  int at1,at2,at1a,at2a; /* atoms indexes */
  int vi1,vi2; /* vla indexes */
  int index1,index2; /* indices in the selection array */
  AtomInfoType *ai1,*ai2,*ai1a,*ai2a; /* atom information pointers */
  ObjectMolecule *obj1,*obj2;
  int cmp;
  PRINTFD(G,FB_Selector) 
    " SelectorCreateAlignments-DEBUG: entry.\n"
    ENDFD
  cnt = 0;
  np = VLAGetSize(pair)/2;
  if(np) {

    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1); /* unnecessary? */
    flag1=Calloc(int,I->NAtom);
    flag2=Calloc(int,I->NAtom);

    /* we need to create two selection arrays: for the matched 
     * atoms in the original selections */
    p = pair;
    for(i=0;i<np;i++) { /* iterate through all pairs of matched residues */
      vi1 = *(p++);
      vi2 = *(p++);
        
      /* find positions in the selection arrays */

      mod1 = vla1[vi1*3];
      at1 = vla1[vi1*3+1];
      
      mod2 = vla2[vi2*3];
      at2 = vla2[vi2*3+1];

      PRINTFD(G,FB_Selector) 
        " S.C.A.-DEBUG: mod1 %d at1 %d mod2 %d at2 %d\n",mod1,at1,mod2,at2
        ENDFD

      obj1 = I->Obj[mod1];
      obj2 = I->Obj[mod2];

      ai1 = obj1->AtomInfo+at1;
      ai2 = obj2->AtomInfo+at2;
      at1a = at1;
      at2a = at2;
      ai1a = ai1;
      ai2a = ai2;
      
      if(atomic_input) {
        index1 = SelectorGetObjAtmOffset(I, obj1, at1a);
        index2 = SelectorGetObjAtmOffset(I, obj2, at2a);
        flag1[index1] = true;
        flag2[index2] = true; 
        cnt++;
      } else {
        
        while(1) { /* match up all matching atom names in each residue */
          cmp = AtomInfoNameOrder(G,ai1a,ai2a);
          if(cmp==0) { /* atoms match */
            index1 = SelectorGetObjAtmOffset(I, obj1, at1a);
            index2 = SelectorGetObjAtmOffset(I, obj2, at2a);
            
            PRINTFD(G,FB_Selector) 
              " S.C.A.-DEBUG: compare %s %s %d\n",ai1a->name,ai2a->name,cmp
              ENDFD
              
              PRINTFD(G,FB_Selector) 
              " S.C.A.-DEBUG: entry %d %d\n",
              ai1a->selEntry,ai2a->selEntry
              ENDFD
              if((index1>=0) && (index2>=0)) {
                if(SelectorIsMember(G,ai1a->selEntry,sele1)&&
                   SelectorIsMember(G,ai2a->selEntry,sele2)) {
                  if((!identical)||(strcmp(ai1a->resn,ai2a->resn)==0)) {
                    flag1[index1] = true;
                    flag2[index2] = true; 
                    cnt++;
                  }
                }
              }
            at1a++;
            at2a++;
          } else if(cmp<0) { /* 1 is before 2 */
            at1a++;
          } else if(cmp>0) { /* 1 is after 2 */
            at2a++;
          }
          if(at1a>=obj1->NAtom) 
            break;
          if(at2a>=obj2->NAtom)
            break;
          ai1a = obj1->AtomInfo+at1a;
          ai2a = obj2->AtomInfo+at2a;
          /* make sure we're still in the same residue */
          if(!AtomInfoSameResidue(G,ai1a,ai1))
            break;
          if(!AtomInfoSameResidue(G,ai2a,ai2))
            break;
        }
      }
    }
    if(cnt) {
      SelectorEmbedSelection(G,flag1,name1,NULL,false, -1);
      SelectorEmbedSelection(G,flag2,name2,NULL,false, -1);
    }
    FreeP(flag1);
    FreeP(flag2);
  }
  PRINTFD(G,FB_Selector) 
    " SelectorCreateAlignments-DEBUG: exit, cnt = %d.\n",cnt
    ENDFD

  return cnt;
}
/*========================================================================*/
int SelectorCountStates(PyMOLGlobals *G,int sele)
{
  register CSelector *I=G->Selector;
  int a;
  int result=0;
  int n_frame;
  int at1;
  ObjectMolecule *last=NULL;
  ObjectMolecule *obj;
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  if(I->NAtom) {
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      obj=I->Obj[I->Table[a].model];
      if(obj!=last) {
        at1=I->Table[a].atom;
        if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele)) {
          if(obj->Obj.fGetNFrame) {
            n_frame=obj->Obj.fGetNFrame((CObject*)obj);
            if(result<n_frame)
              result=n_frame;
          }
          last=obj;
        }
      }
    }
  }
  return(result);
}
/*========================================================================*/
int SelectorCheckIntersection(PyMOLGlobals *G,int sele1,int sele2)
{
  register CSelector *I=G->Selector;
  int a;
  int at1;
  ObjectMolecule *obj;
  
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  if(I->NAtom) {
    for(a=cNDummyAtoms;a<I->NAtom;a++)
      {
        obj=I->Obj[I->Table[a].model];
        at1=I->Table[a].atom;
        if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele1) &&
           SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele2))
          return 1;
      }
  }
  return 0;
}

/*========================================================================*/
int SelectorCountAtoms(PyMOLGlobals *G,int sele,int state)
{
  register CSelector *I=G->Selector;
  int a;
  int result=0;
  int at1;
  ObjectMolecule *obj;
  
  SelectorUpdateTable(G,state,-1);
  if(I->NAtom) {
    for(a=cNDummyAtoms;a<I->NAtom;a++)
      {
        obj=I->Obj[I->Table[a].model];
        at1=I->Table[a].atom;
        if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele)) {
          result++;
        }
      }
  }
  return(result);
}


/*========================================================================*/
int *SelectorGetResidueVLA(PyMOLGlobals *G,int sele,int ca_only,ObjectMolecule *exclude)
{
  /* returns a VLA containing atom indices followed by residue integers
   (residue names packed as characters into integers)
   The indices are the first and last residue in the selection...
  */
  register CSelector *I=G->Selector;
  int *result = NULL,*r;
  int a;
  int c;
  AtomInfoType *ai1 = NULL,*ai2;
  int at1=0,at2;
  unsigned int rcode;
  ResName rn;
  int mod1=0;
  ObjectMolecule *obj;

  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);

  result = VLAlloc(int,I->NAtom*3);

  r = result;
  PRINTFD(G,FB_Selector)
    " SelectorGetResidueVLA-DEBUG: entry, sele = %d\n",sele
    ENDFD;

  if(I->NAtom) {
    if(ca_only) {
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        mod1 = I->Table[a].model;
        at1=I->Table[a].atom;
        obj=I->Obj[mod1];
        ai1 = obj->AtomInfo+at1;
        if(obj!=exclude) {
          if(SelectorIsMember(G,ai1->selEntry,sele)) {
            if(strcmp(ai1->name,"CA")==0) {
              *(r++)=mod1;
              *(r++)=at1;
              for(c=0;c<sizeof(ResName);c++)
                rn[c]=0;
              strcpy(rn,ai1->resn); /* store residue code as a number */
              rcode = 0;
              for(c=0;c<3;c++) {
                rcode = (rcode<<8) | rn[c];
              }
              *(r++) = rcode;
            }
          }
        }
      }
    } else {
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        obj=I->Obj[I->Table[a].model];
        if(obj!=exclude) {
          at2=I->Table[a].atom;
          if(SelectorIsMember(G,obj->AtomInfo[at2].selEntry,sele)) {
            if(!ai1) {
              mod1 = I->Table[a].model;
              at1 = at2;
              ai1 = obj->AtomInfo+at1;
            }
            ai2=obj->AtomInfo+at2;
            if(!AtomInfoSameResidue(G,ai1,ai2)) {
              if(ai1) {
                *(r++)=mod1;
                *(r++)=at1;
                for(c=0;c<sizeof(ResName);c++)
                  rn[c]=0;
                strcpy(rn,ai1->resn); /* store residue code as a number */
                rcode = 0;
                for(c=0;c<3;c++) {
                  rcode = (rcode<<8) | rn[c];
                }
                *(r++) = rcode;
                
                at1 = at2;
                ai1 = ai2;
                mod1 = I->Table[a].model;
              }
            }
          }
        }
      }
      if(ai1) { /* handle last residue */
        *(r++)=mod1;
        *(r++)=at1;
        for(c=0;c<sizeof(ResName);c++)
          rn[c]=0;
        strcpy(rn,ai1->resn); /* store residue code as a number */
        rcode = 0;
        for(c=0;c<3;c++) {
          rcode = (rcode<<8) | rn[c];
        }
        *(r++) = rcode;
      }
    }
  }
  if(result) {
    VLASize(result,int,(ov_size)(r-result));
  }
  PRINTFD(G,FB_Selector)
    " SelectorGetResidueVLA-DEBUG: exit, result = %p, size = %d\n",
    (void*)result,(unsigned int)VLAGetSize(result)
    ENDFD;
  
  return(result);
}
/*========================================================================*/
static int *SelectorGetIndexVLA(PyMOLGlobals *G,int sele) /* assumes updated tables */
{
  register CSelector *I=G->Selector;
  int a,c=0;
  int *result = NULL;
  ObjectMolecule *obj;
  int at1;

  result = VLAlloc(int,(I->NAtom/10)+1);
  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    obj=I->Obj[I->Table[a].model];
    at1=I->Table[a].atom;
    if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele)) {
      VLACheck(result,int,c);
      result[c++]=a;
    }
  }
  VLASize(result,int,c);
  return(result);
}
/*========================================================================*/
void SelectorUpdateObjectSele(PyMOLGlobals *G,ObjectMolecule *obj)
{
  if(obj->Obj.Name[0]) {
    SelectorDelete(G,obj->Obj.Name);  
    SelectorCreate(G,obj->Obj.Name,NULL,obj,true,NULL); /* create a selection with same name */ 
    if(SettingGetGlobal_b(G,cSetting_auto_classify_atoms))
      SelectorClassifyAtoms(G,0,false,obj);
  }
}

/*========================================================================*/
void SelectorLogSele(PyMOLGlobals *G,char *name)
{
  register CSelector *I=G->Selector;
  int a;
  OrthoLineType line,buf1;
  int cnt=-1;
  int first = 1;
  int append=0;
  ObjectMolecule *obj;
  int at1;
  int sele;
  int logging;
  int robust;
  logging = (int)SettingGet(G,cSetting_logging);
  robust = (int)SettingGet(G,cSetting_robust_logs);
  if(logging) {
    sele = SelectorIndexByName(G,name);
    if(sele>=0) {
      SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        obj=I->Obj[I->Table[a].model];
        at1=I->Table[a].atom;
        if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele)) {
          if(cnt<0) {
            if(first) {
              switch(logging) {
              case cPLog_pml:
                sprintf(line,"_ cmd.select(\"%s\",\"(",name);
                break;
              case cPLog_pym:
                sprintf(line,"cmd.select(\"%s\",\"(",name);
                break;
              }
              append=0;
              cnt=0;
              first=0;
            } else {
              switch(logging) {
              case cPLog_pml:
                sprintf(line,"_ cmd.select(\"%s\",\"(%s",name,name);
                break;
              case cPLog_pym:
                sprintf(line,"cmd.select(\"%s\",\"(%s",name,name);
                break;
              }
              append=1;
              cnt=0;
            }
          }
          if(append) 
            strcat(line,"|");
          if(robust) 
            ObjectMoleculeGetAtomSeleFast(obj,at1,buf1);
          else 
            sprintf(buf1,"%s`%d",obj->Obj.Name,at1+1);
          strcat(line,buf1);
          append=1;
          cnt++;
          if(strlen(line)>(sizeof(OrthoLineType)/2)) {
            strcat(line,")\")\n");
            PLog(G,line,cPLog_no_flush);
            cnt=-1;
          }
        }
      }
      if(cnt>0) {
        strcat(line,")\")\n");
        PLog(G,line,cPLog_no_flush);
        PLogFlush(G);
      }
    }
  }
}
/*========================================================================*/
#ifndef _PYMOL_INLINE
int SelectorIsMemberSlow(PyMOLGlobals *G,int s,int sele)
{
  register int s_reg;
  register MemberType *member = G->Selector->Member;
  if( (s_reg=s) && (sele>1) ) {
    register int sele_reg = sele;
    register MemberType *mem = member + s_reg;
    register int test_sele;
    do {
      test_sele = mem->selection;
      s_reg = mem->next;
      if(test_sele==sele_reg) {
        return mem->tag;
      }
      mem = member + s_reg;
    } while(s_reg);
    return false;
  } else if(!sele) 
    return true; /* "all" is selection number 0, unordered */
  else 
    return false; /* no atom is a member of none (1), and negative selections don't exist */
}
#endif
/*========================================================================*/
#if 0
static int SelectorPurgeMember(PyMOLGlobals *G,AtomInfoType *ai,int sele)
{/* never tested...*/
  register CSelector *I=G->Selector;
  int s=ai->selEntry;
  int l=-1;
  int result = 0;
  int nxt;
  while(s)
    {
      nxt = I->Member[s].next;
      if(I->Member[s].selection==sele)
        {
          if(l>0)
            I->Member[l].next=I->Member[s].next;
          else
            ai->selEntry=I->Member[s].next;
          I->Member[s].next = I->FreeMember; 
          I->FreeMember=s;
          result=true;
        }
      l=s;
      s=nxt;
    }
  return result;
}
#endif
/*========================================================================*/
int SelectorMoveMember(PyMOLGlobals *G,int s,int sele_old,int sele_new)
{
  register CSelector *I=G->Selector;
  int result = false;
  while(s) {
    if(I->Member[s].selection==sele_old) {
      I->Member[s].selection=sele_new;
      result = true;
    }
    s=I->Member[s].next;
  }
  return result;
}

/*========================================================================*/
static int SelectorIndexByID(PyMOLGlobals *G,int id)
{
  register CSelector *I=G->Selector;
  int i=0;
  int result = -1;
  SelectionInfoRec *info = I->Info;
  while(i<I->NActive) {
    if((info++)->ID == id) {
      result = i;
      break;
    }
    i++;
  }
  return result;
}
/*========================================================================*/
ObjectMolecule *SelectorGetFastSingleObjectMolecule(PyMOLGlobals *G,int sele)
{
  register CSelector *I=G->Selector;
  ObjectMolecule *result=NULL;
  SelectionInfoRec *info;
  int sele_idx = SelectorIndexByID(G,sele);
  if((sele_idx>=0)&&(sele_idx<I->NActive)) {
    info = I->Info + sele_idx;
    if(info->justOneObjectFlag) {
      if(ExecutiveValidateObjectPtr(G,(CObject*)info->theOneObject,cObjectMolecule))
        result = info->theOneObject;
    } else {
      result = SelectorGetSingleObjectMolecule(G,sele); /* fallback onto slow approach */
    }
  }
  return(result);
}
/*========================================================================*/
ObjectMolecule *SelectorGetFastSingleAtomObjectIndex(PyMOLGlobals *G,int sele,int *index)
{
  register CSelector *I=G->Selector;
  ObjectMolecule *result=NULL;
  SelectionInfoRec *info;
  int got_it = false;
  int sele_idx = SelectorIndexByID(G,sele);

  if((sele_idx>=0)&&(sele_idx<I->NActive)) {
    info = I->Info + sele_idx;
    if(info->justOneObjectFlag && info->justOneAtomFlag) {
      ObjectMolecule *obj = info->theOneObject;
      int at = info->theOneAtom;
      if(ExecutiveValidateObjectPtr(G,(CObject*)obj,cObjectMolecule)) {
        if((at<obj->NAtom)&&SelectorIsMember(G,obj->AtomInfo[at].selEntry,sele)) {
          result = obj;
          *index = at;
          got_it=true;
        }
      }
    }
    if(!got_it) { /* fallback onto slow approach */
      if(!SelectorGetSingleAtomObjectIndex(G,sele,&result,index))
        result = NULL;
    }
  }
  return(result);
}
/*========================================================================*/
ObjectMolecule *SelectorGetSingleObjectMolecule(PyMOLGlobals *G,int sele)
{
  /* slow way */

  int a;
  ObjectMolecule *result = NULL;
  ObjectMolecule *obj;
  register CSelector *I=G->Selector;
  int at1;
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);

  for(a=0;a<I->NAtom;a++) {
    obj=I->Obj[I->Table[a].model];
    at1=I->Table[a].atom;
    if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele)) {
      if(result) {
        if(result!=obj) {
          result=NULL;
          break;
        }
      } else {
        result=obj;
      }
    }
  }
  return(result);
}
/*========================================================================*/
ObjectMolecule *SelectorGetFirstObjectMolecule(PyMOLGlobals *G,int sele)
{
  /* slow way */

  int a;
  ObjectMolecule *result = NULL;
  ObjectMolecule *obj;
  register CSelector *I=G->Selector;
  int at1;
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);

  for(a=0;a<I->NAtom;a++) {
    obj=I->Obj[I->Table[a].model];
    at1=I->Table[a].atom;
    if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele)) {
      result = obj;
      break;
    }
  }
  return(result);
}
/*========================================================================*/
ObjectMolecule **SelectorGetObjectMoleculeVLA(PyMOLGlobals *G,int sele)
{
  int a;
  ObjectMolecule *last = NULL;
  ObjectMolecule *obj,**result = NULL;  
  register CSelector *I=G->Selector;
  int at1;
  int n=0;
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);

  result = VLAlloc(ObjectMolecule*,10);
  for(a=cNDummyAtoms;a<I->NAtom;a++)
    {
      obj=I->Obj[I->Table[a].model];
      at1=I->Table[a].atom;
      if(SelectorIsMember(G,obj->AtomInfo[at1].selEntry,sele)) {
        if(obj!=last) {
          VLACheck(result,ObjectMolecule*,n);
          result[n]=obj;
          last=obj;
          n++;
        }
      }
    }
  VLASize(result,ObjectMolecule*,n);
  return(result);
}

/*========================================================================*/
int SelectorGetSingleAtomObjectIndex(PyMOLGlobals *G,int sele,ObjectMolecule **in_obj,int *index)
{
  /* slow way */

  int found_it = false;
  int a;
  void *iterator = NULL;
  ObjectMolecule *obj = NULL;

  while(ExecutiveIterateObjectMolecule(G,&obj,&iterator)) {
    int n_atom = obj->NAtom;
    AtomInfoType *ai = obj->AtomInfo;
    for(a=0;a<n_atom;a++) {
      register int s = (ai++)->selEntry;
      if(SelectorIsMember(G,s,sele)) {
        if(found_it){
          return false; /* ADD'L EXIT POINT */
        } else {
          found_it = true;
          (*in_obj) = obj;
          (*index) = a;
        }
      }
    }
  }
  return(found_it);
}

/*========================================================================*/
int SelectorGetSingleAtomVertex(PyMOLGlobals *G,int sele,int state,float *v)
{
  ObjectMolecule *obj = NULL;
  int index = 0;
  int found_it = false;
  if(SelectorGetSingleAtomObjectIndex(G,sele,&obj,&index))
    found_it = ObjectMoleculeGetAtomTxfVertex(obj,state,index,v);
  return(found_it);
}
/*========================================================================*/
void SelectorDeletePrefixSet(PyMOLGlobals *G,char *pref)
{
  ov_diff a;
  register CSelector *I=G->Selector;
  SelectorWordType name_copy; 
  int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);

  while(1) {
    a = SelectGetNameOffset(G,pref,strlen(pref),ignore_case);
    if(a>0) {
      strcpy(name_copy,I->Name[a]);
      ExecutiveDelete(G,name_copy); /* import to use a copy, otherwise 
                                   * you'll delete all objects  */
    } else
      break;
  }
}

/*========================================================================*/
#define MAX_DEPTH 1000

static int SelectorCheckNeighbors(PyMOLGlobals *G,int maxDist,ObjectMolecule *obj,int at1,int at2,
                           int *zero,int *scratch)
{
  int s;
  int a,a1;
  int stkDepth = 0;
  int si = 0;
  int stk[MAX_DEPTH];
  int dist=0;

  zero[at1]=dist;
  scratch[si++]=at1;
  stk[stkDepth]=at1;
  stkDepth++;

  while(stkDepth) { /* this will explore a tree */
    stkDepth--;
    a=stk[stkDepth];
    dist = zero[a]+1;

    s=obj->Neighbor[a]; /* add neighbors onto the stack */
    s++; /* skip count */
    while(1) {
      a1 = obj->Neighbor[s];
      if(a1==at2) {
        while(si--) {
          zero[scratch[si]]=0;
        }
        /* EXIT POINT 1 */
        return 1;
      }
      if(a1>=0) {
        if((!zero[a1])&&(stkDepth<MAX_DEPTH)&&(dist<maxDist)) {
          zero[a1]=dist;
          scratch[si++]=a1;
          stk[stkDepth]=a1;
          stkDepth++;
        }
      } else 
        break;
      s+=2;
    }
  }
  while(si--) {
    zero[scratch[si]]=0;
  }
  /* EXIT POINT 2 */
  return 0;
}

/*========================================================================*/
int SelectorWalkTree(PyMOLGlobals *G,int *atom,int *comp,int *toDo,int **stk,
                     int stkDepth,ObjectMolecule *obj,
                     int sele1,int sele2,int sele3,int sele4)
{
  int s;
  int c = 0;
  int a,a1;
  int seleFlag;
  AtomInfoType *ai;

  while(stkDepth) { /* this will explore a tree, stopping at protected atoms */
    stkDepth--;
    a=(*stk)[stkDepth];
    toDo[a]=0;
    seleFlag=false;
    ai=obj->AtomInfo+a;
    s=ai->selEntry;
    seleFlag = SelectorIsMember(G,s,sele1);
    if(!seleFlag)
      seleFlag = SelectorIsMember(G,s,sele2);      
    if(!seleFlag)
      seleFlag = SelectorIsMember(G,s,sele3);      
    if(!seleFlag)
      seleFlag = SelectorIsMember(G,s,sele4);      
    if(!seleFlag) {
      if(!(ai->protekted==1)) { /* if not explicitly protected...*/
        atom[a]=1; /* mark this atom into the selection */
        comp[a]=1;
      }
      s=obj->Neighbor[a]; /* add neighbors onto the stack */
      s++; /* skip count */
      while(1) {
        a1 = obj->Neighbor[s];
        if(a1>=0) {
          if(toDo[a1]) {
            VLACheck((*stk),int,stkDepth);
            (*stk)[stkDepth]=a1;
            stkDepth++;
          }
        } else 
          break;
        s+=2;
      }
      c++;
    }
  }
  return (c);
}

/*========================================================================*/
static int SelectorWalkTreeDepth(PyMOLGlobals *G,int *atom,int *comp,int *toDo,int **stk,
                                 int stkDepth,ObjectMolecule *obj,
                                 int sele1,int sele2,int sele3,int sele4,
                                 int **extraStk, WalkDepthRec *wd)
{
  int s;
  int c = 0;
  int a,a1;
  int seleFlag;
  int depth;
  AtomInfoType *ai;

  wd->depth1 = -1;
  wd->depth2 = -1;
  wd->depth3 = -1;
  wd->depth4 = -1;
  VLACheck(*extraStk,int,stkDepth);
  UtilZeroMem(*extraStk,sizeof(int)*stkDepth);

  while(stkDepth) { /* this will explore a tree, stopping at protected atoms */
    stkDepth--;
    a=(*stk)[stkDepth];
    depth = ((*extraStk)[stkDepth]+1);
    seleFlag=false;
    ai=obj->AtomInfo+a;
    s=ai->selEntry;

    /* record how many cycles it take to reach each & any picked atoms */

    seleFlag = false;
    if(SelectorIsMember(G,s,sele1)) {
      if(((wd->depth1<0)||(wd->depth1>depth))) {
        wd->depth1 = depth;
      }
      seleFlag = true;
    }
    if(SelectorIsMember(G,s,sele2)) {
      if(((wd->depth2<0)||(wd->depth2>depth))) {
        wd->depth2 = depth;
      }
      seleFlag = true;
    }
    if(SelectorIsMember(G,s,sele3)) {
      if(((wd->depth3<0)||(wd->depth3>depth))) {
        wd->depth3 = depth;
      }
      seleFlag = true;
    }
    if(SelectorIsMember(G,s,sele4)) {
      if(((wd->depth4<0)||(wd->depth4>depth))) {
        wd->depth4 = depth;
      }
      seleFlag = true;
    }

    if(!seleFlag) {
      toDo[a]=0;
      if(!(ai->protekted==1)) { /* if not explicitly protected...*/
        atom[a]=1; /* mark this atom into the selection */
        comp[a]=1;
      }
      s=obj->Neighbor[a]; /* add neighbors onto the stack */
      s++; /* skip count */
      while(1) {
        a1 = obj->Neighbor[s];
        if(a1>=0) {
          if(toDo[a1]) {
            VLACheck((*stk),int,stkDepth);
            (*stk)[stkDepth]=a1;
            VLACheck((*extraStk),int,stkDepth);
            (*extraStk)[stkDepth]=depth;
            stkDepth++;
          }
        } else 
          break;
        s+=2;
      }
      c++;
    }
  }
  return (c);
}

/*========================================================================*/



int SelectorIsAtomBondedToSele(PyMOLGlobals *G,ObjectMolecule *obj,int sele1atom,int sele2)
{
  int a0,a2,s,ss;
  int bonded =false;
  ObjectMoleculeUpdateNeighbors(obj);
  
  a0 = ObjectMoleculeGetAtomIndex(obj,sele1atom);
  
  if(a0>=0) { 
    s=obj->Neighbor[a0]; 
    s++; /* skip count */
    while(1) {
      a2 = obj->Neighbor[s];
      if(a2<0)
        break;
      ss=obj->AtomInfo[a2].selEntry;
      if(SelectorIsMember(G,ss,sele2)) {
        bonded = true;
        break;
      }
      s+=2;
    }
  }
  return bonded;
}



static void update_min_walk_depth(WalkDepthRec *minWD,
                            int frag, WalkDepthRec *wd,
                            int sele1, int sele2,
                            int sele3, int sele4)
{
  /* first, does this fragment even qualify ? */
  int qualifies = true;
  int cnt=0;
  wd->sum = 0;
  if(sele1>=0) {
    if(wd->depth1<0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth1;
      cnt++;
    }
  }
  if(sele2>=0) {
    if(wd->depth2<0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth2;
      cnt++;
    }
  }
  if(sele3>=0) {
    if(wd->depth3<0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth3;
      cnt++;
    }
  }
  if(sele4>=0) {
    if(wd->depth4<0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth4;
      cnt++;
    }
  }
  if(qualifies&&(cnt>1)) {

    /* is it better than the current min? */

    if((!minWD->frag)||(wd->sum < minWD->sum)) {
      (*minWD)=(*wd);
      minWD->frag = frag;
    }
  }
}
/*========================================================================*/
int SelectorSubdivide(PyMOLGlobals *G,char *pref,int sele1,int sele2,
                      int sele3,int sele4,char *fragPref,char *compName,
                      int *bondMode)
{
  register CSelector *I=G->Selector;
  int a0=0,a1=0,a2;
  int *atom=NULL;
  int *toDo=NULL;
  int *comp=NULL;
  int *pkset=NULL;
  int set_cnt=0;
  int nFrag = 0;
  int *stk=NULL;
  int stkDepth;
  int c,s;
  int cycFlag=false;
  SelectorWordType name,link_sele="";
  ObjectMolecule *obj1=NULL,*obj2=NULL,*obj3=NULL,*obj4=NULL;
  int index1=0,index2=0,index3=0,index4=0;

  /* this is seriously getting out of hand -- need to switch over to arrays soon */

  int *atom1_base=NULL,*atom2_base=NULL,*atom3_base=NULL,*atom4_base=NULL;
  int *toDo1_base=NULL,*toDo2_base=NULL,*toDo3_base=NULL,*toDo4_base=NULL;
  int *comp1_base=NULL,*comp2_base=NULL,*comp3_base=NULL,*comp4_base=NULL;
  int *pkset1_base=NULL,*pkset2_base=NULL,*pkset3_base=NULL,*pkset4_base=NULL;

  PRINTFD(G,FB_Selector)
    " SelectorSubdivideObject: entered...\n"
    ENDFD;
  SelectorDeletePrefixSet(G,pref);
  SelectorDeletePrefixSet(G,fragPref);
  ExecutiveDelete(G,cEditorLink);
  ExecutiveDelete(G,cEditorSet);
  /* delete any existing matches */
  
  obj1 = SelectorGetFastSingleAtomObjectIndex(G,sele1,&index1);
  obj2 = SelectorGetFastSingleAtomObjectIndex(G,sele2,&index2);
  obj3 = SelectorGetFastSingleAtomObjectIndex(G,sele3,&index3);
  obj4 = SelectorGetFastSingleAtomObjectIndex(G,sele4,&index4);

  if(obj1||obj2||obj3||obj4) {
    
    if(obj1) ObjectMoleculeUpdateNeighbors(obj1);
    if(obj2) ObjectMoleculeUpdateNeighbors(obj2);
    if(obj3) ObjectMoleculeUpdateNeighbors(obj3);
    if(obj4) ObjectMoleculeUpdateNeighbors(obj4);

    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1); 
  
    comp = Calloc(int,I->NAtom);
    atom = Alloc(int,I->NAtom);
    toDo = Alloc(int,I->NAtom);
    pkset = Calloc(int,I->NAtom);

    /* NOTE: SeleBase only safe with cSelectorUpdateTableAllStates!  */

    if(obj1) {
      atom1_base  = atom + obj1->SeleBase;
      toDo1_base  = toDo + obj1->SeleBase;
      comp1_base  = comp + obj1->SeleBase;
      pkset1_base  = pkset + obj1->SeleBase;
    }

    if(obj2) {
      atom2_base  = atom + obj2->SeleBase;
      toDo2_base  = toDo + obj2->SeleBase;
      comp2_base  = comp + obj2->SeleBase;
      pkset2_base  = pkset+ obj2->SeleBase;
    }

    if(obj3) {
      atom3_base  = atom + obj3->SeleBase;
      toDo3_base  = toDo + obj3->SeleBase;
      comp3_base  = comp + obj3->SeleBase;
      pkset3_base  = pkset+ obj3->SeleBase;
    }

    if(obj4) {
      atom4_base  = atom + obj4->SeleBase;
      toDo4_base  = toDo + obj4->SeleBase;
      comp4_base  = comp + obj4->SeleBase;
      pkset4_base  = pkset+ obj4->SeleBase;
    }

    stk=VLAlloc(int,100);
    
    { 
      int a;int *p1;
      p1=toDo;
      for(a=0;a<I->NAtom;a++)
        *(p1++)=true;
    }
    
    if(*bondMode) {
      /* verify bond mode, or clear the flag */
      
      *bondMode = false;
      
      if((sele1>=0)&&(sele2>=0)&&(sele3<0)&&(sele4<0)&&
         (obj1==obj2)) { /* two selections only, in same object... */
        
        a0 = index1;
        a1 = index2;
        
        if((a0>=0)&&(a1>=0)) { 
          s=obj1->Neighbor[a0]; /* add neighbors onto the stack */
          s++; /* skip count */
          while(1) {
            a2 = obj1->Neighbor[s];
            if(a2<0)
              break;
            if(a2==a1) {
              *bondMode = true;
              break;
            }
            s+=2;
          }
        }
      }
    }
    
    /* ===== BOND MODE ===== (sele0 and sele1 only) */ 
    
    if(*bondMode) { 
      if(obj1==obj2) { /* just to be safe */

	pkset1_base[a0]=1;
	pkset1_base[a1]=1;
	SelectorEmbedSelection(G,pkset,cEditorBond,NULL,false,-1);                

        a0 = index1;
        if(a0>=0) {
          stkDepth=0;
          s=obj1->Neighbor[a0]; /* add neighbors onto the stack */
          s++; /* skip count */
          while(1) {
            a1 = obj1->Neighbor[s];
            if(a1>=0) {
              if(toDo1_base[a1]) {
                VLACheck(stk,int,stkDepth);
                stk[stkDepth]=a1;
                stkDepth++;
              }
            } else 
              break;
            s+=2;
          }
          UtilZeroMem(atom,sizeof(int)*I->NAtom);
          atom1_base[a0] = 1; /* create selection for this atom alone as fragment base atom */
          comp1_base[a0] = 1;
          sprintf(name,"%s%1d",fragPref,nFrag+1);
          SelectorEmbedSelection(G,atom,name,NULL,false,-1);
          c = SelectorWalkTree(G,atom1_base,comp1_base,toDo1_base,&stk,stkDepth,obj1,sele1,sele2,-1,-1) + 1;
          sprintf(name,"%s%1d",pref,nFrag+1);

          
          /* check for cyclic situation */
          cycFlag=false;
          a2 = index2;
          if(a2>=0) {
            stkDepth=0;
            s=obj1->Neighbor[a2]; /* add neighbors onto the stack */
            s++; /* skip count */
            while(1) {
              a1 = obj1->Neighbor[s];
              if (a1<0) break;
              if((a1>=0)&&(a1!=a0)) {
                if(!toDo1_base[a1]) {
                  cycFlag=true; /* we have a cycle...*/
                  break;
                }
              }
              s+=2;
            }
          }
          if(cycFlag) { /* cyclic situation is a bit complex...*/
            
            a0 = index2;
            if(a0>=0) {
              stkDepth=0;
              s=obj1->Neighbor[a0]; /* add neighbors onto the stack */
              s++; /* skip count */
              while(1) {
                a1 = obj1->Neighbor[s];
                if(a1>=0) {
                  if(toDo1_base[a1]) {
                    VLACheck(stk,int,stkDepth);
                    stk[stkDepth]=a1;
                    stkDepth++;
                  }
                } else 
                  break;
                s+=2;
              }
              atom1_base[a0] = 1; 
              comp1_base[a0] = 1;
              c = SelectorWalkTree(G,atom1_base,comp1_base,toDo1_base,&stk,stkDepth,obj1,sele1,sele2,-1,-1) + 1;
            }
          }
          SelectorEmbedSelection(G,atom,name,NULL,false,-1);
          nFrag++;
        }
        
        if(!cycFlag) {
          a0 = index2;
          if(a0>=0) {
            stkDepth=0;
            s=obj1->Neighbor[a0]; /* add neighbors onto the stack */
            s++; /* skip count */
            while(1) {
              a1 = obj1->Neighbor[s];
              if(a1>=0) {
                if(toDo1_base[a1]) {
                  VLACheck(stk,int,stkDepth);
                  stk[stkDepth]=a1;
                  stkDepth++;
                }
              } else 
                break;
              s+=2;
            }
            
            UtilZeroMem(atom,sizeof(int)*I->NAtom);
            atom1_base[a0] = 1; /* create selection for this atom alone as fragment base atom */
            comp1_base[a0] = 1;
            sprintf(name,"%s%1d",fragPref,nFrag+1);
            SelectorEmbedSelection(G,atom,name,NULL,false,-1);
            c = SelectorWalkTree(G,atom1_base,comp1_base,toDo1_base,&stk,stkDepth,obj1,sele1,sele2,-1,-1) + 1;
            sprintf(name,"%s%1d",pref,nFrag+1);
            SelectorEmbedSelection(G,atom,name,NULL,false,-1);
            nFrag++;
          }
        }
      }
    } else {
      /* ===== WALK MODE ===== (any combination of sele0, sele1, sele2, sele3 */
      
      int *extraStk = VLAlloc(int,50);
      WalkDepthRec curWalk,minWalk;
      minWalk.sum = 0;
      minWalk.frag = 0;
      
      if(obj1) {
        a0 = index1;
        if(a0>=0) {
          pkset1_base[a0]=1;
          set_cnt++;
          comp1_base[a0]=1;
          stkDepth=0;
          s=obj1->Neighbor[a0]; /* add neighbors onto the stack */
          s++; /* skip count */
          while(1) {
            a1 = obj1->Neighbor[s];
            if(a1<0)
              break;
            if(toDo1_base[a1]) {
              stkDepth=1;
              stk[0] = a1;
              UtilZeroMem(atom,sizeof(int)*I->NAtom);
              atom1_base[a1] = 1; /* create selection for this atom alone as fragment base atom */
              comp1_base[a1] = 1;
              sprintf(name,"%s%1d",fragPref,nFrag+1);
              SelectorEmbedSelection(G,atom,name,NULL,false,-1);
              atom1_base[a1] = 0;
              c = SelectorWalkTreeDepth(G,atom1_base,comp1_base,toDo1_base,&stk,
                                        stkDepth,obj1,sele1,sele2,sele3,sele4,
                                        &extraStk, &curWalk );
              if(c) {
                nFrag++;
                sprintf(name,"%s%1d",pref,nFrag);
                SelectorEmbedSelection(G,atom,name,NULL,false,-1);
                update_min_walk_depth(&minWalk,
                                      nFrag, &curWalk, 
                                      sele1, sele2, sele3, sele4);
              }
            }
            s+=2;
          }
        }
      }
      
      if(obj2) {
        a0 = index2;
        if(a0>=0) {
          pkset2_base[a0]=1;
          set_cnt++;
          comp2_base[a0]=1;
          stkDepth=0;
          s=obj2->Neighbor[a0]; /* add neighbors onto the stack */
          s++; /* skip count */
          while(1) {
            a1 = obj2->Neighbor[s];
            if(a1<0)
              break;
            if(toDo2_base[a1]) {
              stkDepth=1;
              stk[0] = a1;
              UtilZeroMem(atom,sizeof(int)*I->NAtom);
              atom2_base[a1] = 1; /* create selection for this atom alone as fragment base atom */
              comp2_base[a1] = 1;
              sprintf(name,"%s%1d",fragPref,nFrag+1);
              SelectorEmbedSelection(G,atom,name,NULL,false,-1);
              atom2_base[a1] = 0;
              c = SelectorWalkTreeDepth(G,atom2_base,comp2_base,toDo2_base,&stk,
                                        stkDepth,obj2,sele1,sele2,sele3,sele4,
                                        &extraStk, &curWalk );
              if(c) {
                nFrag++;
                sprintf(name,"%s%1d",pref,nFrag);
                SelectorEmbedSelection(G,atom,name,NULL,false,-1);
                update_min_walk_depth(&minWalk,
                                      nFrag, &curWalk, 
                                      sele1, sele2, sele3, sele4);
              }
            }
            s+=2;
          }
        }
      }
      
      if(obj3) {
        a0 = index3;
        if(a0>=0) {
          pkset3_base[a0]=1;
          set_cnt++;
          comp3_base[a0]=1;
          stkDepth=0;
          s=obj3->Neighbor[a0]; /* add neighbors onto the stack */
          s++; /* skip count */
          while(1) {
            a1 = obj3->Neighbor[s];
            if(a1<0)
              break;
            if(toDo3_base[a1]) {
              stkDepth=1;
              stk[0] = a1;
              UtilZeroMem(atom,sizeof(int)*I->NAtom);
              atom3_base[a1] = 1; /* create selection for this atom alone as fragment base atom */
              comp3_base[a1] = 1;
              sprintf(name,"%s%1d",fragPref,nFrag+1);
              SelectorEmbedSelection(G,atom,name,NULL,false,-1);
              atom3_base[a1] = 0;
              c = SelectorWalkTreeDepth(G,atom3_base,comp3_base,toDo3_base,&stk,
                                        stkDepth,obj3,sele1,sele2,sele3,sele4,
                                        &extraStk, &curWalk );
              if(c) {
                nFrag++;
                sprintf(name,"%s%1d",pref,nFrag);
                SelectorEmbedSelection(G,atom,name,NULL,false,-1);
                update_min_walk_depth(&minWalk,
                                      nFrag, &curWalk, 
                                      sele1, sele2, sele3, sele4);
                
              }
            }
            s+=2;
          }
        }
      }
      
      if(obj4) {
        a0 = index4;
        if(a0>=0) {
          pkset4_base[a0]=1;
          set_cnt++;
          comp4_base[a0]=1;
          stkDepth=0;
          s=obj4->Neighbor[a0]; /* add neighbors onto the stack */
          s++; /* skip count */
          while(1) {
            a1 = obj4->Neighbor[s];
            if(a1<0)
              break;
            if(toDo4_base[a1]) {
              stkDepth=1;
              stk[0] = a1;
              UtilZeroMem(atom,sizeof(int)*I->NAtom);
              atom4_base[a1] = 1; /* create selection for this atom alone as fragment base atom */
              comp4_base[a1] = 1;
              sprintf(name,"%s%1d",fragPref,nFrag+1);
              SelectorEmbedSelection(G,atom,name,NULL,false,-1);
              atom4_base[a1] = 0;
              c = SelectorWalkTreeDepth(G,atom4_base,comp4_base,toDo4_base,&stk,
                                        stkDepth,obj4,sele1,sele2,sele3,sele4,
                                        &extraStk, &curWalk);
              if(c) {
                nFrag++;
                sprintf(name,"%s%1d",pref,nFrag);
                SelectorEmbedSelection(G,atom,name,NULL,false,-1);
                update_min_walk_depth(&minWalk,
                                      nFrag, &curWalk, 
                                      sele1, sele2, sele3, sele4);
              }
            }
            s+=2;
          }
        }
      }
      
      if(minWalk.frag) { /* create the linking selection if one exists */
        sprintf(link_sele,"%s%d|?pk1|?pk2|?pk3|?pk4",pref,minWalk.frag);
      }
      VLAFreeP(extraStk);
    }
    
    if(set_cnt>1) {
      SelectorEmbedSelection(G,pkset,cEditorSet,NULL,false,-1);                
    }
    
    if(nFrag) {
      SelectorEmbedSelection(G,comp,compName,NULL,false,-1);
    }
    
    if(link_sele[0])
      SelectorCreate(G,cEditorLink,link_sele,NULL,true,NULL);
    
    FreeP(toDo);
    FreeP(atom);
    FreeP(comp);
    FreeP(pkset);
    VLAFreeP(stk);
    SelectorClean(G);
  }
  PRINTFD(G,FB_Selector)
    " SelectorSubdivideObject: leaving...nFrag %d\n",nFrag
    ENDFD;

  return(nFrag);
}

/*========================================================================*/
int SelectorGetSeleNCSet(PyMOLGlobals *G,int sele)
{
  register CSelector *I=G->Selector;

  int a,s,at = 0;
  ObjectMolecule *obj,*last_obj = NULL;
  int result=0;

  if( (obj = SelectorGetFastSingleAtomObjectIndex(G,sele,&at)) ) {
    int a = obj->NCSet;
    CoordSet *cs;
    int idx;

    while(a--) {
      cs = obj->CSet[a];
      
      if(obj->DiscreteFlag) {
        if(cs==obj->DiscreteCSet[at])
          idx=obj->DiscreteAtmToIdx[at];
        else
          idx=-1;
      } else 
        idx=cs->AtmToIdx[at];
      if(idx>=0) {
        result = a+1;
        break;
      }
    }
  } else {
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      obj=I->Obj[I->Table[a].model];
      if(obj!=last_obj) {
        at=I->Table[a].atom;
        s=obj->AtomInfo[at].selEntry;
        if(SelectorIsMember(G,s,sele)) {
          if(result<obj->NCSet) {
            result=obj->NCSet;
            last_obj = obj;
          }
        }
      }
    }
  }
  return(result);
}
/*========================================================================*/
int SelectorGetArrayNCSet(PyMOLGlobals *G,int *array,int no_dummies)
{
  register CSelector *I=G->Selector;
  int a;
  ObjectMolecule *obj;
  int result=0;
  int start= 0;
  if(no_dummies)
    start = cNDummyAtoms;
  for(a=start;a<I->NAtom;a++) {
    if(*(array++)) {
      if(a>=cNDummyAtoms) {
        obj=I->Obj[I->Table[a].model];
        if(result<obj->NCSet) result=obj->NCSet;
      } else {
        if(result<1) result=1; /* selected dummy has at least one CSet */
      }

    }
  }
  return(result);
}
/*========================================================================*/
float SelectorSumVDWOverlap(PyMOLGlobals *G,int sele1,int state1,int sele2,int state2,float adjust)
{
  register CSelector *I=G->Selector;
  int *vla=NULL;
  int c;
  float result=0.0;
  float sumVDW=0.0,dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at1,at2;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj1,*obj2;
  int idx1,idx2;
  int a;

  if(state1<0) state1=0;
  if(state2<0) state2=0;

  if(state1!=state2) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,state1,-1);
  }

  c=SelectorGetInterstateVLA(G,sele1,state1,sele2,state2,2*MAX_VDW+adjust,&vla);
  for(a=0;a<c;a++) {
    a1=vla[a*2];
    a2=vla[a*2+1];

    at1=I->Table[a1].atom;
    at2=I->Table[a2].atom;
    
    obj1=I->Obj[I->Table[a1].model];
    obj2=I->Obj[I->Table[a2].model];

    if((state1<obj1->NCSet)&&(state2<obj2->NCSet)) {
      cs1=obj1->CSet[state1];
      cs2=obj2->CSet[state2];
      if(cs1&&cs2) { /* should always be true */
        
        ai1=obj1->AtomInfo+at1;
        ai2=obj2->AtomInfo+at2;
       
        idx1=cs1->AtmToIdx[at1]; /* these are also pre-validated */
        idx2=cs2->AtmToIdx[at2];
        
        sumVDW=ai1->vdw+ai2->vdw+adjust;
        dist=(float)diff3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2);
        
        if(dist<sumVDW) {
          result+=((sumVDW-dist)/2.0F);
        }
      }
    }
  }
  VLAFreeP(vla);
  return(result);
}
/*========================================================================*/
static int SelectorGetInterstateVLA(PyMOLGlobals *G,
                                    int sele1,int state1,
                                    int sele2,int state2,
                                    float cutoff,int **vla) /* Assumes valid tables */
{
  register CSelector *I=G->Selector;
  MapType *map;
  float *v2;
  int n1,n2;
  int c,i,j,h,k,l;
  int at;
  int a,s,idx;
  ObjectMolecule *obj;
  CoordSet *cs;

  if(!(*vla))
	 (*vla)=VLAlloc(int,1000);

  c=0;
  n1=0;

  for(a=0;a<I->NAtom;a++) {
	 I->Flag1[a]=false;
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(G,s,sele1))
      {
        if(state1<obj->NCSet) 
          cs=obj->CSet[state1];
        else
          cs=NULL;
        if(cs) {
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {
            copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
            I->Flag1[a]=true;
            n1++;
          }
        }
      }
  }
  /* now create and apply voxel map */
  c=0;
  if(n1) {
	 n2=0;
	 map=MapNewFlagged(G,-cutoff,I->Vertex,I->NAtom,NULL,I->Flag1);
	 if(map) {
		MapSetupExpress(map);
		for(a=cNDummyAtoms;a<I->NAtom;a++) {
		  at=I->Table[a].atom;
		  obj=I->Obj[I->Table[a].model];
		  s=obj->AtomInfo[at].selEntry;
        if(SelectorIsMember(G,s,sele2))
          {
            if(state2<obj->NCSet) 
              cs=obj->CSet[state2];
            else
              cs=NULL;
            if(cs) {
              if(obj->DiscreteFlag) {
                if(cs==obj->DiscreteCSet[at])
                  idx=obj->DiscreteAtmToIdx[at];
                else
                  idx=-1;
              } else 
                idx=cs->AtmToIdx[at];
              if(idx>=0) {
                v2 = cs->Coord+(3*idx);
                if(MapExclLocus(map,v2,&h,&k,&l)) {
                  i=*(MapEStart(map,h,k,l));
                  if(i) {
                    j=map->EList[i++];
                    while(j>=0) {
                      if(within3f(I->Vertex+3*j,v2,cutoff)) {
                        VLACheck((*vla),int,c*2+1);
                        *((*vla)+c*2)=j;
                        *((*vla)+c*2+1)=a;
                        c++;
                      }
                      j=map->EList[i++];
                    }
                  }
                }
                n2++;
              }
            }
          }
		}
		MapFree(map);
	 }
  }
  return(c);
}
/*========================================================================*/
int SelectorMapMaskVDW(PyMOLGlobals *G,int sele1,ObjectMapState *oMap,float buffer,int state)
{
  register CSelector *I=G->Selector;
  MapType *map;
  float *v2;
  int n1,n2;
  int a,b,c,i,j,h,k,l;
  int at;
  int s,idx;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  CoordSet *cs;
  int state1,state2;
  int once_flag;

  c=0;
  n1=0;
  SelectorUpdateTable(G,state,-1);

  for(a=0;a<I->NAtom;a++) {
	 I->Flag1[a]=false;
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(G,s,sele1))
      {
        once_flag=true;
        for(state2=0;state2<obj->NCSet;state2++) {
          if(state<0) once_flag=false;
          if(!once_flag) 
            state1=state2;
          else
            state1=state;
          if(state1<obj->NCSet) 
            cs=obj->CSet[state1];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
              I->Flag1[a]=true;
              n1++;
            }
          }
          if(once_flag) break;
        }
      }
  }
  /* now create and apply voxel map */
  c=0;
  if(n1) {
	 n2=0;
	 map=MapNewFlagged(G,-(buffer+MAX_VDW),I->Vertex,I->NAtom,NULL,I->Flag1);
	 if(map) {
		MapSetupExpress(map);
      
      for(a=oMap->Min[0];a<=oMap->Max[0];a++) {      
        for(b=oMap->Min[1];b<=oMap->Max[1];b++) {      
          for(c=oMap->Min[2];c<=oMap->Max[2];c++) {      
            F3(oMap->Field->data,a,b,c)=0.0;            

            v2 = F4Ptr(oMap->Field->points,a,b,c,0);
            
            if(MapExclLocus(map,v2,&h,&k,&l)) {
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                while(j>=0) {
                  ai = I->Obj[I->Table[j].model]->AtomInfo+I->Table[j].atom;
                  if(within3f(I->Vertex+3*j,v2,ai->vdw+buffer)) {
                    F3(oMap->Field->data,a,b,c)=1.0;
                  }
                  j=map->EList[i++];
                }
              }
            }
          }
        }
		}
      oMap->Active=true;
		MapFree(map);
	 }
  }
  return(c);
}

static double max2d(double a,double b)
{
  if(a>b)
    return a;
  else
    return b;
}

static double max6d(double a,double b,double c,double d,double e,double f)
{

  if(d>a) a = d;
  if(e>b) b = e;
  if(f>c) c = f;
  if(b>a) a = b;
  if(c>a) a = c;
  return(a);
}
#define D_SMALL10 1e-10

typedef double AtomSF[11];

/*========================================================================*/
int SelectorMapGaussian(PyMOLGlobals *G,int sele1,ObjectMapState *oMap,
                        float buffer,int state,int normalize,int use_max,int quiet)
{
  register CSelector *I=G->Selector;
  MapType *map;
  float *v2;
  int n1,n2;
  int a,b,c,i,j,h,k,l;
  int at;
  int s,idx;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  CoordSet *cs;
  int state1,state2;
  float *point=NULL,*fp;
  int *sfidx=NULL,*ip;
  float *b_factor=NULL,*bf,bfact;
  float *occup=NULL,*oc;
  int prot;
  int once_flag;
  float d,e_val;
  double sum,sumsq;
  float mean,stdev;  
  double sf[256][11],*sfp;
  AtomSF *atom_sf = NULL;
  double b_adjust = (double)SettingGet(G,cSetting_gaussian_b_adjust);
  double elim = 7.0;
  double rcut2;
  float rcut;
  float max_rcut = 0.0F;
  float b_floor = SettingGet(G,cSetting_gaussian_b_floor);
  float blur_factor = 1.0F;

  {
    float resolution = SettingGet(G,cSetting_gaussian_resolution);
    if(resolution<2.0) resolution = 2.0;
    blur_factor = 2.0F / resolution; /* a gaussion_resolution of 2.0 is considered perfect */
  }

  if(b_adjust>500.0) b_adjust = 500.0; /* constrain to be somewhat reasonable */

  for(a=0;a<256;a++) {
    sf[a][0]=-1.0;
  }

  sf[cAN_H][0] =  0.493002;
  sf[cAN_H][1] = 10.510900;
  sf[cAN_H][2] =  0.322912;
  sf[cAN_H][3] = 26.125700;
  sf[cAN_H][4] =  0.140191;
  sf[cAN_H][5] =  3.142360;
  sf[cAN_H][6] =  0.040810;
  sf[cAN_H][7] = 57.799698;
  sf[cAN_H][8] =  0.003038;
  sf[cAN_H][9] = 0.0;

  /* LP currently using scattering factors of carbon 
     (Roche Pocket viewer relies upon this behavior) */

  sf[cAN_LP][0] =  2.310000;
  sf[cAN_LP][1] = 20.843899;
  sf[cAN_LP][2] =  1.020000;
  sf[cAN_LP][3] = 10.207500;
  sf[cAN_LP][4] =  1.588600;
  sf[cAN_LP][5] =  0.568700;
  sf[cAN_LP][6] =  0.865000;
  sf[cAN_LP][7] = 51.651199;
  sf[cAN_LP][8] =  0.215600;
  sf[cAN_LP][9] = 0.0;

  sf[cAN_C][0] =  2.310000;
  sf[cAN_C][1] = 20.843899;
  sf[cAN_C][2] =  1.020000;
  sf[cAN_C][3] = 10.207500;
  sf[cAN_C][4] =  1.588600;
  sf[cAN_C][5] =  0.568700;
  sf[cAN_C][6] =  0.865000;
  sf[cAN_C][7] = 51.651199;
  sf[cAN_C][8] =  0.215600;
  sf[cAN_C][9] = 0.0;

  sf[cAN_O][0] =  3.048500;
  sf[cAN_O][1] =  13.277100;
  sf[cAN_O][2] =  2.286800;
  sf[cAN_O][3] =  5.701100;
  sf[cAN_O][4] =  1.546300;
  sf[cAN_O][5] =  0.323900;
  sf[cAN_O][6] =  0.867000;
  sf[cAN_O][7] =  32.908897;
  sf[cAN_O][8] =  0.250800;
  sf[cAN_O][9] = 0.0;

  sf[cAN_N][0] = 12.212600;
  sf[cAN_N][1] =  0.005700;
  sf[cAN_N][2] =  3.132200;
  sf[cAN_N][3] =  9.893300;
  sf[cAN_N][4] =  2.012500;
  sf[cAN_N][5] = 28.997499;
  sf[cAN_N][6] =  1.166300;
  sf[cAN_N][7] =  0.582600;
  sf[cAN_N][8] = -11.528999;
  sf[cAN_N][9] = 0.0;
  
  sf[cAN_S][0] =  6.905300;
  sf[cAN_S][1] =  1.467900;
  sf[cAN_S][2] =  5.203400;
  sf[cAN_S][3] = 22.215099;
  sf[cAN_S][4] =  1.437900;
  sf[cAN_S][5] =  0.253600;
  sf[cAN_S][6] =  1.586300;
  sf[cAN_S][7] = 56.172001;
  sf[cAN_S][8] =  0.866900;
  sf[cAN_S][9] = 0.0;

  sf[cAN_Cl][0] = 11.460400;
  sf[cAN_Cl][1] =  0.010400; 
  sf[cAN_Cl][2] =  7.196400;
  sf[cAN_Cl][3] =  1.166200;
  sf[cAN_Cl][4] =  6.255600;
  sf[cAN_Cl][5] = 18.519400;
  sf[cAN_Cl][6] =  1.645500;
  sf[cAN_Cl][7] = 47.778400;
  sf[cAN_Cl][8] =  0.866900;
  sf[cAN_Cl][9] = 0.0;

  sf[cAN_Br][0] = 17.178900;
  sf[cAN_Br][1] =  2.172300;
  sf[cAN_Br][2] =  5.235800;
  sf[cAN_Br][3] = 16.579599;
  sf[cAN_Br][4] =  5.637700;
  sf[cAN_Br][5] =  0.260900;
  sf[cAN_Br][6] =  3.985100;
  sf[cAN_Br][7] = 41.432800;
  sf[cAN_Br][8] =  2.955700;
  sf[cAN_Br][9] = 0.0;

  sf[cAN_I][0] = 20.147200;
  sf[cAN_I][1] = 4.347000;
  sf[cAN_I][2] = 18.994900;
  sf[cAN_I][3] = 0.381400;
  sf[cAN_I][4] = 7.513800;
  sf[cAN_I][5] = 27.765999;
  sf[cAN_I][6] = 2.273500;
  sf[cAN_I][7] = 66.877602;
  sf[cAN_I][8] = 4.071200;
  sf[cAN_I][9] = 0.0;
  
  sf[cAN_F][0] = 3.539200;
  sf[cAN_F][1] = 10.282499;
  sf[cAN_F][2] = 2.641200;
  sf[cAN_F][3] = 4.294400;
  sf[cAN_F][4] = 1.517000;
  sf[cAN_F][5] = 0.261500;
  sf[cAN_F][6] = 1.024300;
  sf[cAN_F][7] = 26.147600;
  sf[cAN_F][8] = 0.277600;
  sf[cAN_F][9] = 0.0;
  
  sf[cAN_K][0] =   8.218599;
  sf[cAN_K][1] =  12.794900;
  sf[cAN_K][2] =   7.439800;
  sf[cAN_K][3] =   0.774800;
  sf[cAN_K][4] =   1.051900;
  sf[cAN_K][5] = 213.186996;
  sf[cAN_K][6] =    0.865900;
  sf[cAN_K][7] =   41.684097;
  sf[cAN_K][8] =    1.422800;
  sf[cAN_K][9] = 0.0;
  
  sf[cAN_Mg][0] = 5.420400;
  sf[cAN_Mg][1] = 2.827500;
  sf[cAN_Mg][2] = 2.173500;
  sf[cAN_Mg][3] = 79.261101;
  sf[cAN_Mg][4] =  1.226900;
  sf[cAN_Mg][5] = 0.380800;
  sf[cAN_Mg][6] =  2.307300;
  sf[cAN_Mg][7] = 7.193700;
  sf[cAN_Mg][8] = 0.858400;
  sf[cAN_Mg][9] = 0.0;

  sf[cAN_Na][0] = 4.762600;
  sf[cAN_Na][1] = 3.285000;
  sf[cAN_Na][2] = 3.173600;
  sf[cAN_Na][3] = 8.842199;
  sf[cAN_Na][4] = 1.267400;
  sf[cAN_Na][5] = 0.313600;
  sf[cAN_Na][6] = 1.112800;
  sf[cAN_Na][7] = 129.423996;
  sf[cAN_Na][8] = 0.676000;
  sf[cAN_Na][9] = 0.0;

  sf[cAN_P][0] = 6.434500;
  sf[cAN_P][1] = 1.906700;
  sf[cAN_P][2] = 4.179100;
  sf[cAN_P][3] = 27.157000;
  sf[cAN_P][4] =  1.780000;
  sf[cAN_P][5] =  0.526000;
  sf[cAN_P][6] =  1.490800;
  sf[cAN_P][7] = 68.164497;
  sf[cAN_P][8] = 1.114900;
  sf[cAN_P][9] = 0.0;
  
  sf[cAN_Zn][0] = 14.074300;
  sf[cAN_Zn][1] = 3.265500;
  sf[cAN_Zn][2] = 7.031800;
  sf[cAN_Zn][3] = 0.233300;
  sf[cAN_Zn][4] = 5.162500;
  sf[cAN_Zn][5] = 10.316299;
  sf[cAN_Zn][6] = 2.410000;
  sf[cAN_Zn][7] = 58.709702;
  sf[cAN_Zn][8] = 1.304100;
  sf[cAN_Zn][9] = 0.0;
  
  sf[cAN_Ca][0] = 8.626600;
  sf[cAN_Ca][1] = 10.442100;
  sf[cAN_Ca][2] = 7.387300;
  sf[cAN_Ca][3] = 0.659900;
  sf[cAN_Ca][4] = 1.589900;
  sf[cAN_Ca][5] = 85.748398;
  sf[cAN_Ca][6] = 1.021100;
  sf[cAN_Ca][7] = 178.436996;
  sf[cAN_Ca][8] = 1.375100;
  sf[cAN_Ca][9] = 0.0;

  
  sf[cAN_Cu][0] = 13.337999;
  sf[cAN_Cu][1] = 3.582800;
  sf[cAN_Cu][2] = 7.167600;
  sf[cAN_Cu][3] = 0.247000;
  sf[cAN_Cu][4] = 5.615800;
  sf[cAN_Cu][5] = 11.396600;
  sf[cAN_Cu][6] = 1.673500;
  sf[cAN_Cu][7] = 64.812599;
  sf[cAN_Cu][8] = 1.191000;
  sf[cAN_Cu][9] = 0.0;
  
  sf[cAN_Fe][0] = 11.769500;
  sf[cAN_Fe][1] = 4.761100;
  sf[cAN_Fe][2] = 7.357300;
  sf[cAN_Fe][3] = 0.307200;
  sf[cAN_Fe][4] = 3.522200;
  sf[cAN_Fe][5] = 15.353500;
  sf[cAN_Fe][6] = 2.304500;
  sf[cAN_Fe][7] = 76.880501;
  sf[cAN_Fe][8] = 1.036900;
  sf[cAN_Fe][9] = 0.0;

  
  sf[cAN_Se][0] = 17.000599;
  sf[cAN_Se][1] =  2.409800;
  sf[cAN_Se][2] = 5.819600;
  sf[cAN_Se][3] = 0.272600;
  sf[cAN_Se][4] = 3.973100 ;
  sf[cAN_Se][5] = 15.237200 ;
  sf[cAN_Se][6] = 4.354300 ;
  sf[cAN_Se][7] = 43.816299;
  sf[cAN_Se][8] = 2.840900;
  sf[cAN_Se][9] = 0.0;
      
  buffer+=MAX_VDW;
  c=0;
  n1=0;
  if(state>=cSelectorUpdateTableEffectiveStates) {
    SelectorUpdateTable(G,state,-1);
  } else {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  }
  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(G,s,sele1))
      {
        once_flag=true;
        for(state1=0;state1<obj->NCSet;state1++) {
          if(state<0) once_flag=false;
          if(!once_flag) 
            state2=state1;
          else
            state2=state;
          if(state2<obj->NCSet) 
            cs=obj->CSet[state2];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              n1++;
            }
          }
          if(once_flag) break;
        }
      }
  }
  point=Alloc(float,3*n1);
  sfidx=Alloc(int,n1);
  b_factor=Alloc(float,n1);
  occup=Alloc(float,n1);
  atom_sf=Alloc(AtomSF,n1);

  if(!quiet) {
    PRINTFB(G,FB_ObjectMap,FB_Details)
      " ObjectMap: Computing Gaussian map for %d atom positions.\n",n1
      ENDFB(G);
  }

  n1 = 0;
  fp=point;
  ip=sfidx;
  bf=b_factor;
  oc=occup;
  for(a=0;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    ai=obj->AtomInfo + at;
    s=ai->selEntry;
    if(SelectorIsMember(G,s,sele1))
      {
        once_flag=true;
        for(state1=0;state1<obj->NCSet;state1++) {
          if(state<0) once_flag=false;
          if(!once_flag) 
            state2=state1;
          else
            state2=state;
          if(state2<obj->NCSet) 
            cs=obj->CSet[state2];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              copy3f(cs->Coord+(3*idx),fp);
              prot = ai->protons;
              if(sf[prot][0]==-1.0F)
                prot=cAN_C;
              bfact = ai->b + (float)b_adjust;
              if(bfact<b_floor)
                bfact = b_floor;
              if((bfact>R_SMALL4)&&(ai->q>R_SMALL4)) {
                fp+=3;
                *(ip++)=prot;
                *(bf++)=bfact;
                *(oc++)=ai->q;
                n1++;
              }
            }
          }
          if(once_flag) break;
        }
      }
  }

  for(a=0;a<n1;a++) {
    double *src_sf;
    
    src_sf = &sf[sfidx[a]][0];
    bfact = b_factor[a];

    for(b=0;b<10;b+=2) {
      double sfa,sfb;
      sfa = src_sf[b];
      sfb = src_sf[b+1];
      
      atom_sf[a][b] = occup[a] * sfa * pow(sqrt1d(4*PI/(sfb+bfact)),3.0);
      atom_sf[a][b+1] = 4*PI*PI/(sfb+bfact);

    }
    
    rcut2 = max6d(0.0,
                  (elim + log(max2d(fabs(atom_sf[a][0]),D_SMALL10)))/atom_sf[a][1],
                  (elim + log(max2d(fabs(atom_sf[a][2]),D_SMALL10)))/atom_sf[a][3],
                  (elim + log(max2d(fabs(atom_sf[a][4]),D_SMALL10)))/atom_sf[a][5],
                  (elim + log(max2d(fabs(atom_sf[a][6]),D_SMALL10)))/atom_sf[a][7],
                  (elim + log(max2d(fabs(atom_sf[a][8]),D_SMALL10)))/atom_sf[a][9]);
    rcut = ((float)sqrt1d(rcut2)) / blur_factor;
    atom_sf[a][10] = rcut;
    if(max_rcut<rcut)
      max_rcut = rcut;
  }

  /* now create and apply voxel map */
  c=0;
  if(n1) {
	 n2=0;
	 map=MapNew(G,-max_rcut,point,n1,NULL);
	 if(map) {
		MapSetupExpress(map);
      sum = 0.0;
      sumsq = 0.0;
      for(a=oMap->Min[0];a<=oMap->Max[0];a++) {
        OrthoBusyFast(G,a-oMap->Min[0],oMap->Max[0]-oMap->Min[0]+1);
        for(b=oMap->Min[1];b<=oMap->Max[1];b++) {      
          for(c=oMap->Min[2];c<=oMap->Max[2];c++) {      
            e_val=0.0;
            v2 = F4Ptr(oMap->Field->points,a,b,c,0);
            if(MapExclLocus(map,v2,&h,&k,&l)) {
              i=*(MapEStart(map,h,k,l));
              if(i) {
                j=map->EList[i++];
                if(use_max) {
                  float e_partial;
                  while(j>=0) {
                    d = (float)diff3f(point+3*j,v2) * blur_factor; /* scale up width */
                    sfp=atom_sf[j];
                    if(d<sfp[10]) {
                      d=d*d;
                      if(d<R_SMALL8) d=R_SMALL8;
                      e_partial=(float)((sfp[0]*exp(-sfp[1]*d))
                                        +(sfp[2]*exp(-sfp[3]*d))
                                        +(sfp[4]*exp(-sfp[5]*d))
                                        +(sfp[6]*exp(-sfp[7]*d))
                                        +(sfp[8]*exp(-sfp[9]*d))) * blur_factor; /* scale down intensity */
                      if(e_partial > e_val)
                        e_val = e_partial;
                    }
                    j=map->EList[i++];
                  }
                } else {
                  while(j>=0) {
                    d = (float)diff3f(point+3*j,v2) * blur_factor; /* scale up width */
                    sfp=atom_sf[j];
                    if(d<sfp[10]) {
                      d=d*d;
                      if(d<R_SMALL8) d=R_SMALL8;
                      e_val+=(float)(
                                     (sfp[0]*exp(-sfp[1]*d))
                                     +(sfp[2]*exp(-sfp[3]*d))
                                     +(sfp[4]*exp(-sfp[5]*d))
                                     +(sfp[6]*exp(-sfp[7]*d))
                                     +(sfp[8]*exp(-sfp[9]*d))) * blur_factor; /* scale down intensity */ 
                    }
                    j=map->EList[i++];
                  }
                }
              }
            }
            F3(oMap->Field->data,a,b,c)=e_val;
            sum+=e_val;
            sumsq+=(e_val*e_val);
            n2++;
          }
        }
      }
      mean = (float)(sum/n2);
      stdev = (float)sqrt1d((sumsq - (sum*sum/n2))/(n2-1));
      if(normalize) {

        if(!quiet) {
          PRINTFB(G,FB_ObjectMap,FB_Details)
            " ObjectMap: Normalizing: mean = %8.6f & stdev = %8.6f.\n"
            ,mean,stdev
            ENDFB(G);
        }
        
        if(stdev<R_SMALL8)
          stdev=R_SMALL8;
        
        for(a=oMap->Min[0];a<=oMap->Max[0];a++) {      
          for(b=oMap->Min[1];b<=oMap->Max[1];b++) {      
            for(c=oMap->Min[2];c<=oMap->Max[2];c++) {      
              fp = F3Ptr(oMap->Field->data,a,b,c);
              
              *fp = (*fp-mean)/stdev;
            }
          }
        }
      } else {
        if(!quiet) {
          PRINTFB(G,FB_ObjectMap,FB_Details)
            " ObjectMap: Not normalizing: mean = %8.6f and stdev = %8.6f.\n",
            mean,stdev
            ENDFB(G);
        }
      }
      oMap->Active=true;
      MapFree(map);
    }
  }
  FreeP(point);
  FreeP(sfidx);
  FreeP(atom_sf);
  FreeP(b_factor);
  FreeP(occup);
  return(c);
}

/*========================================================================*/
int SelectorMapCoulomb(PyMOLGlobals *G,int sele1,ObjectMapState *oMap,
		       float cutoff,int state,
                       int neutral,int shift,float shift_power)
{
  register CSelector *I=G->Selector;
  MapType *map;
  float *v2;
  register int a,b,c,j,i;
  int h,k,l;
  int at;
  int s,idx;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  CoordSet *cs;
  int state1,state2;
  int once_flag;
  int n_at=0;
  double tot_charge = 0.0;
  float *point=NULL;
  float *charge=NULL;
  int n_point = 0;
  int n_occur;
  float *v0,*v1;
  float c_factor = 1.0F;
  float cutoff_to_power = 1.0F;
  const float _1=1.0F;

  if(shift)
    cutoff_to_power = (float)pow(cutoff,shift_power);

  c_factor=SettingGet(G,cSetting_coulomb_units_factor)/
    SettingGet(G,cSetting_coulomb_dielectric);

  c=0;
  SelectorUpdateTable(G,state,-1);

  point = VLAlloc(float,I->NAtom*3);
  charge = VLAlloc(float,I->NAtom);
    
  /* first count # of times each atom appears */

  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    ai = obj->AtomInfo + at;
    if(SelectorIsMember(G,s,sele1))
      {
        n_occur = 0;
        /* count */
        once_flag=true;
        for(state2=0;state2<obj->NCSet;state2++) {
          if(state<0) once_flag=false;
          if(!once_flag) 
            state1=state2;
          else
            state1=state;
          if(state1<obj->NCSet) 
            cs=obj->CSet[state1];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              n_occur++;
              n_at++;
            }
          }
          if(once_flag) break;
        }
        /* copy */
        if(n_occur) {
          once_flag=true;
          for(state2=0;state2<obj->NCSet;state2++) {
            if(state<0) once_flag=false;
            if(!once_flag) 
              state1=state2;
            else
              state1=state;
            if(state1<obj->NCSet) 
              cs=obj->CSet[state1];
            else
              cs=NULL;
            if(cs) {
              if(obj->DiscreteFlag) {
                if(cs==obj->DiscreteCSet[at])
                  idx=obj->DiscreteAtmToIdx[at];
                else
                  idx=-1;
              } else 
                idx=cs->AtmToIdx[at];
              if(idx>=0) {
                VLACheck(point,float,3*n_point+2);
                VLACheck(charge,float,n_point);
                v0=cs->Coord+(3*idx);
                v1=point+3*n_point;
                copy3f(v0,v1);
                charge[n_point]=ai->partialCharge*ai->q/n_occur;

                tot_charge+=charge[n_point];
                n_point++;
              }
            }
            if(once_flag) break;
          }
        }
      }
  }

  PRINTFB(G,FB_Selector,FB_Details)
    " SelectorMapCoulomb: Total charge is %0.3f for %d points (%d atoms).\n",tot_charge,n_point,n_at
    ENDFB(G);

  if(neutral&&(fabs(tot_charge)>R_SMALL4)) {
    float adjust;

    adjust = (float)(-tot_charge/n_point);

    for(a=0;a<n_point;a++) {
      charge[a]+=adjust;
    }

    PRINTFB(G,FB_Selector,FB_Details)
      " SelectorMapCoulomb: Setting net charge to zero...\n"
      ENDFB(G);
    
  }

  for(a=0;a<n_point;a++) { /* premultiply c_factor by charges */
    charge[a]*=c_factor;
  }


  /* now create and apply voxel map */
  c=0;
  if(n_point) {
    register int *min = oMap->Min;
    register int *max = oMap->Max;
    register CField *data= oMap->Field->data;
    register CField *points = oMap->Field->points;
    register float dist;

    if(cutoff>0.0F) {/* we are using a cutoff */
      if(shift) {
        PRINTFB(G,FB_Selector,FB_Details)
          " SelectorMapCoulomb: Evaluating local Coulomb potential for grid (shift=%0.2f)...\n",
          cutoff
          ENDFB(G);        
      } else {
        PRINTFB(G,FB_Selector,FB_Details)
          " SelectorMapCoulomb: Evaluating Coulomb potential for grid (cutoff=%0.2f)...\n",cutoff
          ENDFB(G);
      }

      map=MapNew(G,-(cutoff),point,n_point,NULL);
      if(map) {
        register int *elist;
        register float dx,dy,dz;
        register float cut = cutoff;
        register float cut2 = cutoff*cutoff;
        
        MapSetupExpress(map);
        elist = map->EList;      
        for(a=min[0];a<=max[0];a++) {      
          OrthoBusyFast(G,a-min[0],max[0]-min[0]+1);
          for(b=min[1];b<=max[1];b++) {      
            for(c=min[2];c<=max[2];c++) {      
              F3(data,a,b,c)=0.0F;            
              v2 = F4Ptr(points,a,b,c,0);
              
              if(MapExclLocus(map,v2,&h,&k,&l)) {
                i=*(MapEStart(map,h,k,l));
                if(i) {
                  j=elist[i++];
                  while(j>=0) {
                    v1 = point + 3*j;
                    while(1) {
                      
                      dx = v1[0]-v2[0];
                      dy = v1[1]-v2[1];
                      dx = (float)fabs(dx);
                      dy = (float)fabs(dy);
                      if(dx>cut) break;
                      dz = v1[2]-v2[2];
                      dx = dx * dx;
                      if(dy>cut) break;
                      dz = (float)fabs(dz);
                      dy = dy * dy;
                      if(dz>cut) break;
                      dx = dx + dy;
                      dz = dz * dz;
                      if(dx>cut2) break;
                      dy = dx+dz;
                      if(dy>cut2) break;
                      dist = (float)sqrt1f(dy);

                      if(dist>R_SMALL4) {
                        if(shift) {
                          if(dist<cutoff) {
                            F3(data,a,b,c) +=  (charge[j]/dist)*
                              (_1-(float)pow(dist,shift_power)/cutoff_to_power);
                          } 
                        } else {
                          F3(data,a,b,c) +=  charge[j]/dist;
                        }
                      }

                      break;
                    } 
                    j=map->EList[i++];
                  }
                }
              }
            }
          }
        }
        MapFree(map);
      } 
	 } else {
      register float *v1;
      PRINTFB(G,FB_Selector,FB_Details)
        " SelectorMapCoulomb: Evaluating Coulomb potential for grid (no cutoff)...\n"
        ENDFB(G);

      for(a=min[0];a<=max[0];a++) {      
        OrthoBusyFast(G,a-min[0],max[0]-min[0]+1);
        for(b=min[1];b<=max[1];b++) {      
          for(c=min[2];c<=max[2];c++) {  
            F3(data,a,b,c)=0.0F;            
            v1 = point;
            v2 = F4Ptr(points,a,b,c,0);
            for(j=0;j<n_point;j++) {
              dist = (float)diff3f(v1,v2);
              v1+=3;
              if(dist>R_SMALL4) {
                F3(data,a,b,c) +=  charge[j]/dist;
              }
            }
          }
        }
      }
    }
    oMap->Active=true;
  }
  VLAFreeP(point);
  VLAFreeP(charge);
  return(1);
}

/*========================================================================*/
int SelectorGetPDB(PyMOLGlobals *G,char **charVLA,int cLen,int sele,int state,
                   int conectFlag,PDBInfoRec *pdb_info,int *counter, double *ref,
                   ObjectMolecule *single_object)
{
  register CSelector *I=G->Selector;

  int a,b,b1,b2,c,d,s,idx,at,a1,a2;
  int use_ter = (int)SettingGet(G,cSetting_pdb_use_ter_records);
  int retain_ids = (int)SettingGet(G,cSetting_pdb_retain_ids);
  int conect_all = (int)SettingGet(G,cSetting_pdb_conect_all);
  double matrix[16];
  int matrix_flag = false;
  float v_tmp[3],*v_ptr;
  CoordSet *cs,*mat_cs = NULL;
  ObjectMolecule *obj,*last_obj=NULL;
  AtomInfoType *atInfo,*ai,*last = NULL;
  
  if(!single_object) {
    SelectorUpdateTable(G,state,-1);
  } else {
    SelectorUpdateTableSingleObject(G,single_object,state,
                                    false,NULL,0,false);
  }

  if(pdb_info->is_pqr_file)
    use_ter = false;
  if(counter)
    c=*counter;
  else
    c=0;
    /*  if(SettingGet(G,cSetting_save_pdb_ss)) {
  SSEntry *ss = NULL;
  int n_ss = 0;
  int ss_active = false;

      ss = VLAlloc(SSEntry,100);
      
      for(a=0;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      I->Table[a].index=0;
      obj=I->Obj[I->Table[a].model];
      s = obj->AtomInfo[at].selEntry;
      if(SelectorIsMember(G,s,sele)) 
        {
          if(state<obj->NCSet) 
            cs=obj->CSet[state];
          else
            cs=NULL;
          if(cs) {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[at])
                idx=obj->DiscreteAtmToIdx[at];
              else
                idx=-1;
            } else 
              idx=cs->AtmToIdx[at];
            if(idx>=0) {
              ai = obj->AtomInfo+at;

              if(ss_active) {
                
              } else {
                if((at->ss=='H')||(at->ss=='S')) {
                  VLACheck(ss,SSEntry,n_ss);
                  
                }
              }
              
              CoordSetAtomToPDBStrVLA(G,charVLA,&cLen,ai,
                                      obj->CSet[state]->Coord+(3*idx),c);
              last = ai;
              c++;
            }
          }
        }
    }
  }
    */

  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    I->Table[a].index=0;
    obj=I->Obj[I->Table[a].model];
    s = obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(G,s,sele)) 
      {
        if((state>=0)&&(state<obj->NCSet)) 
          cs=obj->CSet[state];
        else
          cs=NULL;
        if(cs) {
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {

            if(mat_cs!=cs) {
              /* compute the effective matrix for output coordinates */

              matrix_flag = false;
              if(ObjectGetTotalMatrix(&obj->Obj,state,false,matrix)) {
                if(ref) {
                  left_multiply44d44d(ref,matrix);
                }
                matrix_flag=true;
              } else if(ref) {
                copy44d(ref,matrix);
                matrix_flag=true;
              }
              mat_cs = cs;
            }

            ai = obj->AtomInfo+at;
            if(last)
              if(!last->hetatm)
                if(ai->resv!=last->resv)
                  if((abs(ai->resv-last->resv)>1)||(ai->hetatm)) {
                    if(use_ter) {
                      CoordSetAtomToTERStrVLA(G,charVLA,&cLen,last,c);
                      c++;
                    }
                  }
            if(retain_ids) {
              I->Table[a].index = ai->id;
            } else {
              I->Table[a].index = c+1; /* NOTE marking with "1" based indexes here */
            }
            v_ptr = cs->Coord+(3*idx);
            if(matrix_flag) {
              transform44d3f(matrix,v_ptr,v_tmp);
              v_ptr = v_tmp;
            }
            CoordSetAtomToPDBStrVLA(G,charVLA,&cLen,ai,v_ptr,c,pdb_info);
            last = ai;
            c++;

            if(!conect_all) {
              int conect_all_tmp = 0;
              if(last_obj!=obj) {
                if(SettingGetIfDefined_b(G,obj->Obj.Setting,cSetting_pdb_conect_all,&conect_all_tmp)) {
                  if(conect_all_tmp)
                    conect_all = true;
                }
              }
              last_obj = obj;
            }
          }
        }
      }
  }
  if(conectFlag&&!(pdb_info->is_pqr_file)) {
    register BondType *bond=NULL;
    register BondType *ii1;
    int nBond=0;
    int newline;

    nBond = 0;
    bond = VLACalloc(BondType,1000);
    for(a=cNDummyModels;a<I->NModel;a++) {
      obj=I->Obj[a];
      ii1=obj->Bond;
      if((state>=0)&&(state<obj->NCSet))
        cs=obj->CSet[state];
      else
        cs=NULL;
      if(cs) {
        atInfo=obj->AtomInfo;
        for(b=0;b<obj->NBond;b++) {
          b1=ii1->index[0];
          b2=ii1->index[1];   
          if(obj->DiscreteFlag) {
            if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
              a1=obj->DiscreteAtmToIdx[b1];
              a2=obj->DiscreteAtmToIdx[b2];
            } else {
              a1=-1;
              a2=-1;
            }
          } else {
            a1=cs->AtmToIdx[b1];
            a2=cs->AtmToIdx[b2];
          }
          
          if((a1>=0)&&(a2>=0)&&(conect_all||atInfo[b1].hetatm||atInfo[b2].hetatm)) {
            int i_b1 = SelectorGetObjAtmOffset(I,obj,b1);
            int i_b2 = SelectorGetObjAtmOffset(I,obj,b2);
            if((i_b1>=0)&&(i_b2>=0)) {
              if(I->Table[i_b1].index&&I->Table[i_b2].index) {
                VLACheck(bond,BondType,nBond+(2*ii1->order)+2);
                b1 = I->Table[i_b1].index;
                b2 = I->Table[i_b2].index;
                for(d=0;d<ii1->order;d++) {
                  bond[nBond].index[0] = b1;
                  bond[nBond].index[1] = b2;
                  nBond++;
                  bond[nBond].index[0] = b2;
                  bond[nBond].index[1] = b1;
                  nBond++;
                }
              }
            }
          }
        ii1++;
        }
      }
    }
    {
      register char *reg_cVLA = *charVLA;
      register int reg_cLen = cLen;
      UtilSortInPlace(G,bond,nBond,sizeof(BondType),(UtilOrderFn*)BondInOrder);
      ii1=bond;
      b1=-1;
      b2=-1;
      newline = false;
      for(a=0;a<nBond;a++) {
        if(a<(nBond-1)) {
          if((ii1->index[0]==(ii1+1)->index[0])&&(ii1->index[1]==(ii1+1)->index[1])) 
            newline=true;
        }
        if((b1!=ii1->index[0])||((b1==ii1->index[0])&&(b2==ii1->index[1]))||newline) {
          VLACheck(reg_cVLA,char,reg_cLen+255);
          if(a) reg_cLen+=sprintf(reg_cVLA+reg_cLen,"\n");
          reg_cLen+=sprintf(reg_cVLA+reg_cLen,"CONECT%5d%5d",
                        ii1->index[0],ii1->index[1]);
          b1=ii1->index[0];
          b2=ii1->index[1];
          newline=false;
          if(a>0) {
            if(((ii1-1)->index[0]==ii1->index[0])&&((ii1-1)->index[1]==ii1->index[1])) 
              newline=true;        
          }
        } else {
          VLACheck(reg_cVLA,char,reg_cLen+255);
          reg_cLen+=sprintf(reg_cVLA+reg_cLen,"%5d",
                        ii1->index[1]);
        }
        b2=ii1->index[1];
        ii1++;
      }
      if(reg_cLen) {
        VLACheck(reg_cVLA,char,reg_cLen+255);
        if(*(reg_cVLA+reg_cLen-1)!='\n')
          reg_cLen+=sprintf(reg_cVLA+reg_cLen,"\n");
      }
      (*charVLA) = reg_cVLA;
      cLen = reg_cLen;
    }
    VLAFree(bond);
  }
  /*
    VLAFreeP(ss); */
  if(counter)
    *counter=c;
  return(cLen);
}
/*========================================================================*/
PyObject *SelectorGetChemPyModel(PyMOLGlobals *G,int sele,int state,double *ref)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  register CSelector *I=G->Selector;
  PyObject *model = NULL;
  int ok = true;

  SelectorUpdateTable(G,state,-1);

  model = PyObject_CallMethod(P_models,"Indexed","");
  if (!model) 
    ok = ErrMessage(G,"CoordSetAtomToChemPyAtom","can't create model");  
  if(ok) {    
    int nAtom = 0;
    int a;
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      int at = I->Table[a].atom;
      ObjectMolecule *obj = I->Obj[I->Table[a].model];
      int s = obj->AtomInfo[at].selEntry;
      I->Table[a].index = 0;
      if(SelectorIsMember(G,s,sele)) {
        CoordSet *cs;
        if(state<obj->NCSet) 
          cs=obj->CSet[state];
        else
          cs=NULL;
        if(cs) {
          int idx;
          if(obj->DiscreteFlag) {
            if(cs == obj->DiscreteCSet[at])
              idx=obj->DiscreteAtmToIdx[at];
            else
              idx=-1;
          } else 
            idx=cs->AtmToIdx[at];
          if(idx>=0) {
            I->Table[a].index = nAtom+1; /* NOTE marking with "1" based indexes here */
            nAtom++;
          }
        }
      }
    }

    if(nAtom) {
      int single_flag = true;
      CoordSet *single_cs = NULL;
      CoordSet *mat_cs = NULL;
      {
        int c = 0;
        PyObject *atom_list = PyList_New(nAtom);
        PyObject_SetAttrString(model,"atom",atom_list);
        for(a=cNDummyAtoms;a<I->NAtom;a++) {
          if(I->Table[a].index) {
            ObjectMolecule *obj;
            CoordSet *cs;
            int idx;
            int at = I->Table[a].atom;
            obj = I->Obj[I->Table[a].model];
            cs = obj->CSet[state]; /* assuming this is valid... */
            if(obj->DiscreteFlag) {
              if(obj->CSet[state] == obj->DiscreteCSet[at])
                idx = obj->DiscreteAtmToIdx[at];
              else
                idx = -1;
            } else 
              idx = cs->AtmToIdx[at];
            if(idx>=0) {
              double matrix[16];
              int matrix_flag = false;

              if(mat_cs!=cs) {
                /* compute the effective matrix for output coordinates */
              
                matrix_flag = false;
                if(ObjectGetTotalMatrix(&obj->Obj,state,false,matrix)) {
                  if(ref) {
                    left_multiply44d44d(ref,matrix);
                  }
                  matrix_flag=true;
                } else if(ref) {
                  copy44d(ref,matrix);
                  matrix_flag=true;
                }
                mat_cs = cs;
              }

              if(single_flag) { /* remember whether all atoms come from a single coordinate set...*/
                if(single_cs) {
                  if(single_cs!=cs)
                    single_flag=false;
                } else {
                  single_cs = cs;
                }
              }

              {
                AtomInfoType *ai = obj->AtomInfo + at;
                float *v_ptr = cs->Coord+(3*idx);
                float v_tmp[3];
                if(matrix_flag) {
                  transform44d3f(matrix,v_ptr,v_tmp);
                  v_ptr = v_tmp;
                }
                {
                  RefPosType *ref_pos = cs->RefPos; /* could be NULL */
                  float *ref_ptr = NULL;
                  float ref_tmp[3];
                  if(ref_pos) {
                    ref_pos += idx;
                    if(ref_pos->specified) {
                      ref_ptr = ref_pos->coord;
                      if(matrix_flag) {
                        transform44d3f(matrix,ref_ptr,ref_tmp);
                        ref_ptr = ref_tmp;
                      }
                    }
                  }
                  if(c<nAtom) { /* safety */
                    PyObject *atom =  CoordSetAtomToChemPyAtom(G,ai,v_ptr,ref_ptr,at);
                    if(atom) {
                      PyList_SetItem(atom_list,c,atom); /* steals */
                      c = c + 1;
                    }
                  }
                }
              }
            }
          }
        }
        if(c != nAtom) {
          ok = false;
        }
        Py_XDECREF(atom_list); 
      }

      if(single_flag&&single_cs) { /* single coordinate set?  then set coordinate set info */
        PyObject *molecule = PyObject_GetAttrString(model,"molecule");
        if(molecule) {
          if(single_cs->Name[0]) {
            PyObject_SetAttrString(molecule,"title", /* including name/title */
                                   PyString_FromString(single_cs->Name)); 
          }
        } else {
          ok=false;
        }
        Py_XDECREF(molecule);
      }

      {
        BondType *bond = VLACalloc(BondType,1000);
        int nBond = 0;
        for(a=cNDummyModels;a<I->NModel;a++) {
          BondType *ii1;
          CoordSet *cs;
          ObjectMolecule *obj = I->Obj[a];
          ii1=obj->Bond;
          if(state<obj->NCSet) 
            cs=obj->CSet[state];
          else
            cs=NULL;
          if(cs) {
            int b;
            for(b=0;b<obj->NBond;b++) {
              int a1,a2;
              int b1=ii1->index[0];
              int b2=ii1->index[1];        
              if(obj->DiscreteFlag) {
                if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
                  a1=obj->DiscreteAtmToIdx[b1];
                  a2=obj->DiscreteAtmToIdx[b2];
                } else {
                  a1=-1;
                  a2=-1;
                }
              } else {
                a1=cs->AtmToIdx[b1];
                a2=cs->AtmToIdx[b2];
              }

              if((a1>=0)&&(a2>=0)) {
                int i_b1 = SelectorGetObjAtmOffset(I,obj,b1);
                int i_b2 = SelectorGetObjAtmOffset(I,obj,b2);
                if((i_b1>=0)&&(i_b2>=0)) {
                  if(I->Table[i_b1].index&&I->Table[i_b2].index) { /* selected atoms will be nonzero */
                    VLACheck(bond,BondType,nBond);
                    bond[nBond] = *ii1; /* copy all fields */
                    bond[nBond].index[0] = I->Table[i_b1].index - 1; /* counteract 1-based */
                    bond[nBond].index[1] = I->Table[i_b2].index - 1; /* indexing from above */
                    nBond++;
                  }
                }
              }
              ii1++;
            }
          }
          
          if(cs && (nAtom == cs->NIndex)) { /* support for experimental spheroids - likely to change */
            if(cs->Spheroid && cs->NSpheroid && cs->SpheroidNormal) {
              PyObject *tmp = PConvFloatArrayToPyList(cs->Spheroid,cs->NSpheroid);
              PyObject_SetAttrString(model,"spheroid",tmp);
              Py_XDECREF(tmp);          
              tmp = PConvFloatArrayToPyList(cs->SpheroidNormal,cs->NSpheroid*3);
              PyObject_SetAttrString(model,"spheroid_normals",tmp);
              Py_XDECREF(tmp);          
            }
          }

          {
            PyObject *bond_list = PyList_New(nBond);
            if(bond_list) {
              BondType *ii1 = bond;
              int b;
              PyObject_SetAttrString(model,"bond",bond_list);
              for(b=0;b<nBond;b++) {
                PyObject *bnd = PyObject_CallMethod(P_chempy,"Bond","");
                if(bnd) {
                  PConvInt2ToPyObjAttr(bnd,"index",ii1->index);
                  PConvIntToPyObjAttr(bnd,"order",ii1->order);
                  PConvIntToPyObjAttr(bnd,"id",ii1->id);
                  PConvIntToPyObjAttr(bnd,"stereo",ii1->stereo);
                  PyList_SetItem(bond_list,b,bnd); /* steals bnd reference */
                } else {
                  ok=false;
                  break;
                }
                ii1++;
              }
            } else {
              ok=false;
            }
            Py_XDECREF(bond_list);
          }
        }
        VLAFree(bond);
      }
    }
  }
  if(!ok) {
    if(model) {
      Py_DECREF(model);
    }
    model=NULL; 
  }
  return(model);
#endif
}
/*========================================================================*/
void SelectorUpdateCmd(PyMOLGlobals *G,int sele0,int sele1,int sta0, int sta1,
                       int matchmaker,int quiet)
{
  register CSelector *I=G->Selector;
  int a,b;
  int at0=0,at1;
  int *vla0=NULL;
  int *vla1=NULL;
  int c0,c1;
  int i0=0,i1;
  int cc1;
  ObjectMolecule *obj0=NULL,*obj1;
  CoordSet *cs0,*cs1;
  int matched_flag;
  int b_start;
  int ci0,ci1;
  int ccc = 0;


  PRINTFD(G,FB_Selector)
    " SelectorUpdateCmd-Debug: entered sta0 %d sta1 %d",sta0,sta1
    ENDFD;

  if((sta0<0)||(sta1<0)||(sta0!=sta1)) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,sta0,-1);
  }

  vla0 = SelectorGetIndexVLA(G,sele0);
  vla1 = SelectorGetIndexVLA(G,sele1);

  if(!(vla0&&vla1))
    ErrMessage(G,"Update","no coordinates updated.");
  else {
    c0 = VLAGetSize(vla0);
    c1 = VLAGetSize(vla1);

    b = 0;
    for(a=0;a<c1;a++) { /* iterate over source atoms */
      /* NOTE, this algorithm is N^2 and slow in the worst case...
         however the best case (N) is quite common, especially when merging 
         files written out of PyMOL */

      i1 = vla1[a];
      at1=I->Table[i1].atom;
      obj1=I->Obj[I->Table[i1].model];
      matched_flag=false;

      switch(matchmaker) {
      case 0: /* simply assume that atoms are stored in PyMOL in the identical order, one for one */
        if(b<c0) {
          i0 = vla0[b];
          at0=I->Table[i0].atom;
          obj0=I->Obj[I->Table[i0].model];
          b++;
          matched_flag=true;
        }
        break;
      case 1: /* match each pair based on atom info */
        b_start = b;
        matched_flag=false;
        while(1) {
          i0 = vla0[b];
          at0=I->Table[i0].atom;
          obj0=I->Obj[I->Table[i0].model];
          if(obj0!=obj1) {
            if(AtomInfoMatch(G,obj1->AtomInfo + at1,
                             obj0->AtomInfo + at0)) {
              matched_flag=true;
              break;
            }
	  } else if(at0 == at1) {
	    matched_flag=true;
	    break;
	  }
          b++;
          if(b>=c0)
            b = 0;
          if(b==b_start) 
            break;
        }
        break;
      case 2: /* match based on ID */
       { 
          int target = obj1->AtomInfo[at1].id;
          b_start = b;
          matched_flag=false;
          while(1) {
            i0 = vla0[b];
            at0=I->Table[i0].atom;
            obj0=I->Obj[I->Table[i0].model];
            if(obj0!=obj1) {
              if(obj0->AtomInfo[at0].id==target) {
                matched_flag=true;
                break;
              }
	    } else if(at0 == at1) {
	      matched_flag=true;
	      break;
	    }
            b++;
            if(b>=c0)
              b = 0;
            if(b==b_start) 
              break;
          }
        }
        break;
      case 3: /* match based on rank */
         { 
          int target = obj1->AtomInfo[at1].rank;
          b_start = b;
          matched_flag=false;
          while(1) {
            i0 = vla0[b];
            at0=I->Table[i0].atom;
            obj0=I->Obj[I->Table[i0].model];
            if(obj0!=obj1) {
              if(obj0->AtomInfo[at0].rank==target) {
                matched_flag=true;
                break;
              }
	    } else if(at0 == at1) {
	      matched_flag=true;
	    }
            b++;
            if(b>=c0)
              b = 0;
            if(b==b_start) 
              break;
          }
        }
        break;
      case 4: /* match based on index */
         { 
          b_start = b;
          matched_flag=false;
          while(1) {
            i0 = vla0[b];
            at0=I->Table[i0].atom;
            obj0=I->Obj[I->Table[i0].model];
            if(obj0!=obj1) {
              if(at0 == at1) {
                matched_flag=true;
                break;
              }
	    } else if(at0 == at1) {
	      matched_flag=true;
	      break;
	    }
            b++;
            if(b>=c0)
              b = 0;
            if(b==b_start) 
              break;
          }
         }
        break;
      }
      
      if(matched_flag) { /* atom matched, so copy coordinates */
        ccc++;
        for(cc1=0;cc1<obj1->NCSet;cc1++) { /* iterate over all source states */
          if((cc1==sta1)||(sta1<0)) {
            cs1 = obj1->CSet[cc1];
            if(cs1&&
               (((sta0<0)&&(cc1<obj0->NCSet))|| /* multiple states */
                (cc1==sta0)|| /* single state */
                ((sta0>=0)&&(sta1>=0)))) { /* explicit state */

              if((sta0<0)||(sta0>=obj0->NCSet)) {
                cs0 = obj0->CSet[cc1];
              } else if(sta0<obj0->NCSet) {
                cs0 = obj0->CSet[sta0];
              } else {
                cs0 = NULL;
              }

              if(cs0) {

                /* old broken code: 
                   ci0 = cs0->AtmToIdx[at0];
                   ci1 = cs1->AtmToIdx[at1];
                */

                if(obj0->DiscreteFlag) {
                  if(cs0==obj0->DiscreteCSet[at0])
                    ci0=obj0->DiscreteAtmToIdx[at0];
                  else
                    ci0=-1;
                } else 
                  ci0 = cs0->AtmToIdx[at0]; 
                
                if(obj1->DiscreteFlag) {
                  if(cs1==obj1->DiscreteCSet[at1])
                    ci1=obj1->DiscreteAtmToIdx[at1];
                  else
                    ci1=-1;
                } else 
                  ci1 = cs1->AtmToIdx[at1]; 

                if((ci0>=0)&&(ci1>=0))
                  copy3f(cs1->Coord + 3*ci1,
                         cs0->Coord + 3*ci0);
              }
            }
          }
        }
      }
    }
    obj0=NULL;
    for(b=0;b<c1;b++) {
      obj1=I->Obj[I->Table[i0].model];
      if(obj0!=obj1) {
        obj0=obj1;
        ObjectMoleculeInvalidate(obj0,cRepAll,cRepInvCoord,-1);
      }
    }
    SceneChanged(G);
    if(!quiet) {
      PRINTFB(G,FB_Selector,FB_Actions)
        " Update: coordinates updated for %d atoms.\n",ccc 
        ENDFB(G);

    }
  }
  VLAFreeP(vla0);
  VLAFreeP(vla1);
}
/*========================================================================*/

int SelectorCreateObjectMolecule(PyMOLGlobals *G,int sele,char *name,
                                 int target,int source,int discrete,
                                 int zoom,int quiet,int singletons)
{
  register CSelector *I=G->Selector;
  int ok=true;
  int a,b,a1,a2,b1,b2,c,d,s,at;
  BondType *ii1,*bond=NULL;
  int nBond=0;
  int nCSet,nAtom,ts;
  AtomInfoType *atInfo = NULL;
  int isNew;
  CoordSet *cs = NULL;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj;
  CObject *ob;
  ObjectMolecule *targ = NULL;
  ObjectMolecule *info_src = NULL;
  int static_singletons = SettingGetGlobal_b(G,cSetting_static_singletons);

  if(singletons<0)
    singletons = static_singletons;

  ob=ExecutiveFindObjectByName(G,name);
  if(ob)
    if(ob->type==cObjectMolecule) 
      targ = (ObjectMolecule*)ob;
  
  c=0;
  if(source<0) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,source,-1);
  }

  if(!targ) {
    isNew=true;
    if(discrete<0)
      discrete = SelectorIsSelectionDiscrete(G,sele,false);
    targ = ObjectMoleculeNew(G,discrete);
    targ->Bond = VLACalloc(BondType,1);
    {
      /* copy object color of previous object (if any) */
      ObjectMolecule *singleObj = NULL;
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        at=I->Table[a].atom;
        I->Table[a].index=-1;
        obj=I->Obj[I->Table[a].model];
        s=obj->AtomInfo[at].selEntry;
        if(SelectorIsMember(G,s,sele)) {
          if(!singleObj)
            singleObj = obj;
          else if(singleObj && (obj!=singleObj)) {
            singleObj = NULL;
            break;
          }
        }
      }
      if(singleObj) 
        targ->Obj.Color = singleObj->Obj.Color;
      /* should also consider copying lots of other stuff from the source object ... */
    }
  } else {
    isNew=false;
  }

  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    I->Table[a].index=-1;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(G,s,sele))
      {
        I->Table[a].index=c; /* Mark records as to which atom they'll be */
        c++;
        if(!info_src) info_src = obj;
      }
  }
  if(isNew&&info_src) {/* copy symmetry information, etc. */
    targ->Symmetry=SymmetryCopy(info_src->Symmetry);
  }
  nAtom=c;

  nBond = 0;
  bond = VLACalloc(BondType,nAtom*4);
  for(a=cNDummyModels;a<I->NModel;a++) { /* find bonds wholly contained in the selection */
    obj=I->Obj[a];
    ii1=obj->Bond;
    for(b=0;b<obj->NBond;b++) {
      b1 = SelectorGetObjAtmOffset(I,obj,ii1->index[0]);
      b2 = SelectorGetObjAtmOffset(I,obj,ii1->index[1]);
      if((b1>=0)&&(b2>=0)) {
        if((I->Table[b1].index>=0)&&(I->Table[b2].index>=0)) {
          VLACheck(bond,BondType,nBond);
          {
            BondType *dst_bond = bond+nBond;
            AtomInfoBondCopy(G, ii1, dst_bond);
            dst_bond->index[0] = I->Table[b1].index; /* store what will be the new index */
            dst_bond->index[1] = I->Table[b2].index;
            /*            printf("Selector-DEBUG %d %d\n",dst_bond->index[0],dst_bond->index[1]);*/
            nBond++;
          }
        }
      }
      ii1++;
    }
  }
  
  atInfo = VLAlloc(AtomInfoType,nAtom); 
  /* copy the atom info records and create new zero-based IDs */
  c=0;
  {
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      if(I->Table[a].index>=0) {
        obj=I->Obj[I->Table[a].model];
        at=I->Table[a].atom;
        VLACheck(atInfo,AtomInfoType,c);
        AtomInfoCopy(G, obj->AtomInfo+at, atInfo+c);
        c++;
      }
    }
  }
    
  cs=CoordSetNew(G);  /* set up a dummy coordinate set for the merge xref */
  cs->NIndex = nAtom;
  cs->fEnumIndices(cs);
  cs->TmpBond = bond; /* load up the bonds */
  cs->NTmpBond = nBond;
  bond=NULL;
  
  /*  printf("Selector-DEBUG nAtom %d\n",nAtom);*/
  ObjectMoleculeMerge(targ,atInfo,cs,false,cAIC_AllMask,true); /* will free atInfo */
  /* cs->IdxToAtm will now have the reverse mapping from the new subset
     to the new merged molecule */

  ObjectMoleculeExtendIndices(targ,-1);
  ObjectMoleculeUpdateNonbonded(targ);
  
  if(!isNew) { /* recreate selection table */

    if(source<0) {
      SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
    } else {
      SelectorUpdateTable(G,source,-1);
    }

  }

  /* get maximum state index for the selection...note that
     we'll be creating states from 1 up to the maximum required to
     capture atoms in the selection 
  */
  
  nCSet = 0;

  {
    c=0;
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      I->Table[a].index=-1;
      obj=I->Obj[I->Table[a].model];
      s=obj->AtomInfo[at].selEntry;
      if(SelectorIsMember(G,s,sele)) {
        I->Table[a].index=c; /* Mark records  */
        if(nCSet<obj->NCSet)
          nCSet=obj->NCSet;
        c++;
      }
    }
  }

  if(c!=nAtom) ErrFatal(G,"SelectorCreate","inconsistent selection.");
  /* cs->IdxToAtm now has the relevant indexes for the coordinate transfer */

  for(d=0;d<nCSet;d++) { /* iterate through states */
    if((source<0)||(source==d)) {
      cs2 = CoordSetNew(G);
      c = 0;
      cs2->Coord = VLAlloc(float,3*nAtom);
      cs2->AtmToIdx = Alloc(int,targ->NAtom+1);
      for(a=0;a<targ->NAtom;a++) 
        cs2->AtmToIdx[a]=-1;
      cs2->NAtIndex = targ->NAtom;
      cs2->IdxToAtm = Alloc(int,nAtom);
      for(a=cNDummyAtoms;a<I->NAtom;a++)  /* any selected atoms in this state? */
        if(I->Table[a].index>=0) {
          at=I->Table[a].atom;
          obj=I->Obj[I->Table[a].model];
          cs1 = NULL;
          if(d<obj->NCSet) {
            cs1 = obj->CSet[d];
          } else if(singletons&&(obj->NCSet==1)) {
            cs1 = obj->CSet[0];
          }
          if(cs1) {
            if((!cs2->Name[0])&&(cs1->Name[0])) /* copy the molecule name (if any) */
              strcpy(cs2->Name,cs1->Name);
            
            if(obj->DiscreteFlag) {
              if(cs1==obj->DiscreteCSet[at])
                a1=obj->DiscreteAtmToIdx[at];
              else
                a1=-1;
            } else 
              a1 = cs1->AtmToIdx[at]; /* coord index in existing object */
            if(a1>=0) {
              copy3f(cs1->Coord+a1*3,cs2->Coord+c*3);
              a2 = cs->IdxToAtm[I->Table[a].index]; /* actual merged atom index */
              cs2->IdxToAtm[c] = a2;
              cs2->AtmToIdx[a2] = c;
              c++;
            }
          }
        }
      cs2->IdxToAtm=Realloc(cs2->IdxToAtm,int,c);
      VLASize(cs2->Coord,float,c*3);
      cs2->NIndex = c;
      if(target>=0) {
        if(source<0)
          ts = target + d;
        else
          ts = target;
      } else {
        ts = d;
      }
      VLACheck(targ->CSet,CoordSet*,ts);
      if(targ->NCSet<=ts) 
        targ->NCSet=ts+1;
      if(targ->CSet[ts])
        targ->CSet[ts]->fFree(targ->CSet[ts]);
      targ->CSet[ts]=cs2;
      cs2->Obj=targ;
    }
  }               
  VLAFreeP(bond); /* null-safe */
  if(cs) cs->fFree(cs);
  if(targ->DiscreteFlag) { /* if the new object is discrete, then eliminate the AtmToIdx array */
    for(d=0;d<targ->NCSet;d++) {
      cs = targ->CSet[d];
      if(cs) {
        if(cs->AtmToIdx) {
          for(a=0;a<cs->NIndex;a++) {
            b = cs->IdxToAtm[a];
            targ->DiscreteAtmToIdx[b] = a;
            targ->DiscreteCSet[b] = cs;
          }
          FreeP(cs->AtmToIdx);
        }
      }
    }
  }
  SceneCountFrames(G);
  if(!quiet) {
    PRINTFB(G,FB_Selector,FB_Details)
      " Selector: found %d atoms.\n",nAtom 
      ENDFB(G);
  }
  ObjectMoleculeSort(targ);
  if(isNew) {
    ObjectSetName((CObject*)targ,name);
    ExecutiveManageObject(G,(CObject*)targ,zoom,quiet);
  } else {
    ExecutiveUpdateObjectSelection(G,(CObject*)targ);
  }
  SceneChanged(G);
  return ok;
}

/*========================================================================*/
int SelectorSetName(PyMOLGlobals *G,char *new_name, char *old_name)
{
 register CSelector *I=G->Selector;
 ov_diff i;
 int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
  
 i = SelectGetNameOffset(G,old_name,1,ignore_case);
 if(i>=0) {
   SelectorDelName(G,i);
   UtilNCopy(I->Name[i], new_name, WordLength);
   SelectorAddName(G,i);
   return true;
 } else {
   return false;
 }
}

/*========================================================================*/
int SelectorIndexByName(PyMOLGlobals *G,char *sname)
{
 OrthoLineType name;
 register CSelector *I=G->Selector;
 int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
 ov_diff i=-1;

 if(sname) {
   char *tname = sname;
   while((tname[0]=='%')||(tname[0]=='?'))
     tname++;
   strcpy(name,tname);		  
   i = SelectGetNameOffset(G,name,1,ignore_case);
   if((i>=0)&&(name[0]!='_')) { /* don't do checking on internal selections */
     char *best;
     best = ExecutiveFindBestNameMatch(G,sname); /* suppress spurious matches
                                                    of selections with non-selections */
     if((best!=sname)&&(strcmp(best,I->Name[i])))
       i=-1;
   }
   if(i>=0) i = I->Info[i].ID;
 }
 return(i);
}
/*========================================================================*/
static void SelectorPurgeMembers(PyMOLGlobals *G,int sele) 
{
  register CSelector *I=G->Selector;
  void *iterator = NULL;
  ObjectMolecule *obj = NULL;
  
  if(I->Member) {
    register MemberType *I_Member = I->Member;
    register int I_FreeMember = I->FreeMember;

    while(ExecutiveIterateObjectMolecule(G,&obj,&iterator)) {
      if(obj->Obj.type == cObjectMolecule) {
        register AtomInfoType *ai = obj->AtomInfo;
        register int a, n_atom = obj->NAtom;
        for(a=0;a<n_atom;a++) {
          register int s =  (ai++)->selEntry;
          register int l = -1;
          while( s ) {
            register MemberType *i_member_s = I_Member + s;
            register int nxt = i_member_s->next;
            if(i_member_s->selection == sele) {
              if(l>0)
                I_Member[l].next = i_member_s->next;
              else
                ai[-1].selEntry = i_member_s->next;
              i_member_s->next = I_FreeMember; 
              I_FreeMember = s;
            }
            l = s;
            s = nxt;
          }
        }
      }
    }
    I->FreeMember = I_FreeMember;
  }
}

/*========================================================================*/
int SelectorPurgeObjectMembers(PyMOLGlobals *G,ObjectMolecule *obj)
{
  int a=0;
  int s=0;
  int nxt;

  register CSelector *I=G->Selector;
  if(I->Member)
    for(a=0;a<obj->NAtom;a++)
      {
        s=obj->AtomInfo[a].selEntry;
        while(s)
          {
            nxt = I->Member[s].next;
            I->Member[s].next = I->FreeMember; 
            I->FreeMember=s;
            s=nxt;
          }
        obj->AtomInfo[a].selEntry=0;
		}
  return 1;
}

/*========================================================================*/
void SelectorDelete(PyMOLGlobals *G,char *sele) 
     /* should (only) be called by Executive or by Selector, unless the
        named selection has never been registered with the executive 
        (i.e., temporary on-the-fly selection) */
{
  ov_diff n;
  n=SelectGetNameOffset(G,sele,999,SettingGetGlobal_b(G,cSetting_ignore_case)); /* already exist? */
  if(n>=0) { /* get rid of existing selection -- but never selection 0 (all) */
    SelectorDeleteSeleAtOffset(G,n);
  }
}
/*========================================================================*/
void SelectorGetUniqueTmpName(PyMOLGlobals *G,char *name_buffer)
{
  register CSelector *I=G->Selector;
  sprintf(name_buffer,"%s%d",cSelectorTmpPrefix,I->TmpCounter++);
}
/*========================================================================*/
int SelectorGetTmp(PyMOLGlobals *G,char *input,char *store)
{
  int count = 0;
  /* ASSUMES that store is at least as big as an OrthoLineType */
  register CSelector *I=G->Selector;
  PRINTFD(G,FB_Selector)
    " SelectorGetTmp-Debug: entered with \"%s\".\n",input
    ENDFD;
  
  store[0] = 0;
  
  /* skip trivial cases */

  if( input[0] && !((input[0]=='\'')&&(input[1]=='\'')&&(!input[2]))) {
    
  /* OKAY, this routine is in flux.  Eventually this routine will...

  (1) fully parse the input recognizing selection keywords, nested
      parens, quotes, escaped strings, etc.
  
  (2) replace selection blocks with temporary selection names

  (3) return a space-separated list of names for processing

  However, right now, this routine simply handles two cases.

  A. where the input is a selection, in which case store is set to a
  temporary selection name

  B. where the input is simply a list of space-separated name patterns,
  in which case store is simply passed along as a copy of the input
  
  */
    
    int is_selection = false;
    char *p = input;
    OrthoLineType word;
    OVreturn_word result;

    while(*p) {
      p=ParseWord(word,p,sizeof(OrthoLineType));
      /* see a paren? then this must be a selection */

      if(word[0]=='(') {
        is_selection=true;
        break;
      }

      /* encounterd a selection keyword? then this must be a selection */

      if(OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,word)))) {
        if(OVreturn_IS_OK( (result = OVOneToAny_GetKey(I->Key, result.word)))) {
          if((result.word != SELE_ALLz)&&
			(result.word != SELE_ORIz)&&
			(result.word != SELE_CENz)) {
            is_selection=true;
            break;
          }
        }
      }
      
      if(!ExecutiveValidName(G,word)) { /* don't recognize the name? */
        if(!ExecutiveValidNamePattern(G,word)) { /* don't recognize this as a pattern? */
          is_selection=true; /* must be a selection */
          break;
        }
      }
    }
    if(is_selection) { /* incur the computational expense of 
                          parsing the input as an atom selection */
      WordType name;
      sprintf(name,"%s%d",cSelectorTmpPrefix,I->TmpCounter++);
      count = SelectorCreate(G,name,input,NULL,false,NULL);
      if(count>=0) {
        strcpy(store,name);
      } else {
        store[0]=0;
      }
    } else { /* otherwise, just parse the input as a space-separated list of names */
      /* not a selection */
      strcpy(store,input);
    }
  }
  PRINTFD(G,FB_Selector)
    " SelectorGetTmp-Debug: leaving with \"%s\".\n",store
    ENDFD;
  return count;

}
int SelectorCheckTmp(PyMOLGlobals *G,char *name)
{
  if(WordMatch(G,cSelectorTmpPattern,name,false)+1 == -(int)strlen(cSelectorTmpPrefix))
    return true;
  else
    return false;
}
/*========================================================================*/
void SelectorFreeTmp(PyMOLGlobals *G,char *name) /* remove temporary selections */
{
  if(name&&name[0])
    if(strncmp(name,cSelectorTmpPrefix,strlen(cSelectorTmpPrefix))==0) {
      ExecutiveDelete(G,name);
    }
}
/*========================================================================*/
static int  SelectorEmbedSelection(PyMOLGlobals *G,int *atom, char *name,
                                   ObjectMolecule *obj,int no_dummies, int exec_managed)
{
  /* either atom or obj should be NULL, not both and not neither */

  register CSelector *I=G->Selector;
  int tag;
  int newFlag=true;
  int n,a,m,sele;
  int c=0;
  int start=0;
  int singleAtomFlag = true;
  int singleObjectFlag = true;
  ObjectMolecule *singleObject = NULL,*selObj;
  int singleAtom = -1;
  int index;
  AtomInfoType *ai;

  if(exec_managed<0) {
    if(atom) /* automatic behavior: manage selections defined via atom masks */
      exec_managed=true;
    else
      exec_managed=false;
  }

  n=SelectGetNameOffset(G,name,999,SettingGetGlobal_b(G,cSetting_ignore_case)); /* already exist? */
  if(n==0) /* don't allow redefinition of "all" */
    return 0;
  if(n>0) /* get rid of existing selection*/ {
    SelectorDelete(G,I->Name[n]);
    newFlag = false;
  }

  /*  printf("I->NMember %d I->FreeMember %d\n",I->NMember,I->FreeMember);*/

  n=I->NActive;
  VLACheck(I->Name,SelectorWordType,n+1);
  VLACheck(I->Info,SelectionInfoRec,n+1);
  strcpy(I->Name[n],name);
  I->Name[n+1][0]=0;
  SelectorAddName(G,n);
  sele = I->NSelection++;
  SelectionInfoInit(I->Info + n);
  I->Info[n].ID = sele;
  I->NActive++;

  if(no_dummies) {
    start = 0;
  } else {
    start = cNDummyAtoms;
  }
  for(a=start;a<I->NAtom;a++) {
    tag = false;
    if(atom) {
      if(atom[a]) tag = atom[a];
    } else {
      if(I->Obj[I->Table[a].model]==obj) tag = 1;
    }
    if(tag) {
      selObj = I->Obj[I->Table[a].model];
      index = I->Table[a].atom;
      ai = selObj->AtomInfo+index;
      
      if(singleObjectFlag) {
        if(singleObject) {
          if(selObj!=singleObject) {
            singleObjectFlag = false;
          }
        } else {
          singleObject = selObj;
        }
      }
      
      if(singleAtomFlag) {
        if(singleAtom>=0) {
          if(index!=singleAtom) {
            singleAtomFlag = false;
          }
        } else {
          singleAtom = index;
        }
      }
      
      c++;
      if(I->FreeMember>0) {
        m=I->FreeMember;
        I->FreeMember=I->Member[m].next;
      } else {
        I->NMember++;
        m=I->NMember;
        VLACheck(I->Member,MemberType,m);
      }
      I->Member[m].selection = sele;
      I->Member[m].tag = tag;
      /* at runtime, selections can now have transient ordering --
         but these are not yet persistent through session saves & restores */
      I->Member[m].next = ai->selEntry;
      ai->selEntry = m;
    }
  }
  
  if(c) { /* if selection contains just one atom/object, then take note */
    
    SelectionInfoRec *info = I->Info + (I->NActive-1);
    if( singleObjectFlag ) {
      info->justOneObjectFlag = true;
      info->theOneObject = singleObject;
      if( singleAtomFlag ) {
        info->justOneAtomFlag = true;
        info->theOneAtom = singleAtom;
      }
    }
  }
  
  if(exec_managed) { 
    if(newFlag)
      ExecutiveManageSelection(G,name);
    else
      ExecutiveSetControlsOff(G,name);
  }
  PRINTFD(G,FB_Selector)
    " Selector: Embedded %s, %d atoms.\n",name,c
    ENDFD;
  return(c);
}
/*========================================================================*/
static int *SelectorApplySeqRowVLA(PyMOLGlobals *G,CSeqRow *rowVLA,int nRow)
{
#if 0

  register CSelector *I=G->Selector;
  int *result;
  int a,n;
  ObjectMolecule *obj;
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);  
  result = Calloc(int,I->NAtom);
  for(a=0;a<nRow;a++) {
    CSeqRow *row = rowVLA + a;
    obj = ExecutiveFindObjectMoleculeByName(G,row->name);
    if(obj) {
      for(b=0;b<row->nCol;b++)
        {
          CSeqCol *col = row->col + b;
          int *index;
          if(col->atom_at&&row->atom_lists)
            index = row->atom_lists + col->atom_at;
          while((*index)>0) {
            if(index)
              if(col->inverse)

                result[obj->SeleBase + *index] = true; 
            /* NOTE: SeleBase only safe with cSelectorUpdateTableAllStates!  */
            index++;
          }
        }
    }
  }
  return(result);
#else
  return NULL;
#endif
}
/*========================================================================*/
static int *SelectorApplyMultipick(PyMOLGlobals *G,Multipick *mp)
{
  register CSelector *I=G->Selector;
  int *result;
  int a,n;
  Picking *p;
  ObjectMolecule *obj;
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);  
  result = Alloc(int,I->NAtom);
  n=mp->picked[0].src.index;
  p=mp->picked+1;
  for(a=0;a<I->NAtom;a++) 
    result[a]=0;
  while(n--) { /* what if this object isn't a molecule object??*/
    obj=(ObjectMolecule*)p->context.object;
    /* NOTE: SeleBase only safe with cSelectorUpdateTableAllStates!  */
    result[obj->SeleBase+p->src.index] = true;
    p++;
  }
  return(result);
}
/*========================================================================*/
static int *SelectorSelectFromTagDict(PyMOLGlobals *G,OVOneToAny *id2tag)
{
  register CSelector *I=G->Selector;
  int *result=NULL;
  register int a;
  register AtomInfoType *ai;
  OVreturn_word ret;
  
  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1); /* for now, update the entire table */
  {
    register TableRec *i_table = I->Table, *table_a;
    register ObjectMolecule **i_obj = I->Obj;
    
    result = Calloc(int,I->NAtom);
    if(result) {
      table_a = i_table + cNDummyAtoms;
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        ai = i_obj[table_a->model]->AtomInfo + table_a->atom;
        if(ai->unique_id) {
          if(!OVreturn_IS_ERROR(ret = OVOneToAny_GetKey(id2tag,ai->unique_id)))
            result[a]=ret.word;
        }
        table_a++;
      }
    }
  }
  return(result);
}

/*========================================================================*/


static int _SelectorCreate(PyMOLGlobals *G,char *sname,char *sele,ObjectMolecule **obj,
                           int quiet,Multipick *mp,CSeqRow *rowVLA,
                           int nRow,int **obj_idx,int *n_idx,int n_obj,
                           OVOneToAny *id2tag, int executive_manage, int state, int domain)
{
  int *atom=NULL;
  OrthoLineType name;
  int ok=true;
  int c=0;
  int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
  ObjectMolecule *embed_obj = NULL;

  PRINTFD(G,FB_Selector)
    "SelectorCreate-Debug: entered...\n"
    ENDFD;

  if(sname[0]=='%')
	 strcpy(name,&sname[1]);
  else
	 strcpy(name,sname);		  
  if(WordMatch(G,cKeywordAll,name,ignore_case)<0) {
    name[0]=0; /* force error */
  }
  UtilCleanStr(name);
  if(!name[0]) {
    PRINTFB(G,FB_Selector,FB_Errors)
      "Selector-Error: Invalid selection name \"%s\".\n",sname
      ENDFB(G);
    OrthoRestorePrompt(G);
  }
  if(ok) {
    if(sele) {
      atom=SelectorSelect(G,sele,state,domain);
      if(!atom) ok=false;
    } else if(id2tag) {
      atom=SelectorSelectFromTagDict(G,id2tag);
    } else if(obj && obj[0]) { /* optimized full-object selection */
      if(n_obj<=0) {
        embed_obj = *obj;
        if(obj_idx && n_idx) {
          atom=SelectorUpdateTableSingleObject(G,embed_obj,cSelectorUpdateTableAllStates,
                                               false,*obj_idx,*n_idx, (n_obj==0));
        } else {
          atom=SelectorUpdateTableSingleObject(G,embed_obj,cSelectorUpdateTableAllStates,
                                               false,NULL,0, (n_obj==0));
        }
      } else {
        atom=SelectorUpdateTableMultiObjectIdxTag(G,obj,false,obj_idx,n_idx,n_obj);
      }
    } else if(mp) {
      atom=SelectorApplyMultipick(G,mp);
    } else if(rowVLA) {
      atom=SelectorApplySeqRowVLA(G,rowVLA,nRow);
    } else 
      ok=false;
  }
  if(ok) c=SelectorEmbedSelection(G,atom,name,embed_obj,false, 
                                  executive_manage);
  FreeP(atom);
  SelectorClean(G);
  if(!quiet) {
    if(name[0]!='_') {
      if(ok) {
        PRINTFB(G,FB_Selector,FB_Actions)
          " Selector: selection \"%s\" defined with %d atoms.\n",name,c
          ENDFB(G);
      }
    }
  }
  if(ok) {
    PRINTFD(G,FB_Selector)
      " SelectorCreate: \"%s\" created with %d atoms.\n",name,c    
      ENDFD;
  } else {
    PRINTFD(G,FB_Selector)
      " SelectorCreate: \"%s\" not created due to error\n",name
      ENDFD;
  }
  if(!ok)
    c=-1;
  return(c);
}

int SelectorCreateFromTagDict(PyMOLGlobals *G,char *sname, OVOneToAny *id2tag, int exec_managed)
{
  return _SelectorCreate(G,sname,NULL,NULL,true,NULL,NULL,0,NULL, NULL,0,id2tag, exec_managed, -1,-1);
}

int SelectorCreateEmpty(PyMOLGlobals *G,char *name,int exec_managed)
{
  return _SelectorCreate(G,name, "none", NULL, 1, NULL, NULL, 0, NULL, 0, 0, NULL, exec_managed, -1,-1);
}
int SelectorCreateSimple(PyMOLGlobals *G,char *name, char *sele)
{
  return _SelectorCreate(G,name, sele, NULL, 1, NULL, NULL, 0, NULL, 0, 0, NULL, -1, -1,-1);  
}
int SelectorCreateFromObjectIndices(PyMOLGlobals *G,char *sname, ObjectMolecule *obj, int *idx, int n_idx)
{
  return _SelectorCreate(G,sname,NULL,&obj,true,NULL,NULL,0,&idx,&n_idx,-1,NULL, -1, -1,-1); /* n_obj = -1 disables numbered tags */
}
int SelectorCreateOrderedFromObjectIndices(PyMOLGlobals *G,char *sname, ObjectMolecule *obj, int *idx, int n_idx)
{
  return _SelectorCreate(G,sname,NULL,&obj,true,NULL,NULL,0,&idx,&n_idx,0,NULL,-1,-1,-1); /* assigned numbered tags */
}
int SelectorCreateOrderedFromMultiObjectIdxTag(PyMOLGlobals *G,char *sname, 
                                               ObjectMolecule **obj,
                                               int **idx_tag,
                                               int *n_idx, int n_obj)
{
  return _SelectorCreate(G,sname,NULL,obj,true,NULL,NULL,0,idx_tag,n_idx,n_obj,NULL,-1,-1,-1);
}
#if 0
static int SelectorCreateFromSeqRowVLA(PyMOLGlobals *G,char *sname,CSeqRow *rowVLA,int nRow)
{
  return _SelectorCreate(G,sname,NULL,NULL,true,NULL,rowVLA,nRow,NULL,0,0,NULL,-1,-1,-1);
}
#endif

int SelectorCreate(PyMOLGlobals *G,char *sname,char *sele,ObjectMolecule *obj,int quiet,Multipick *mp)
{
  return _SelectorCreate(G,sname,sele,&obj,quiet,mp,NULL,0,NULL,0,0,NULL,-1,-1,-1);
}

int SelectorCreateWithStateDomain(PyMOLGlobals *G,char *sname,char *sele,ObjectMolecule *obj,
                                    int quiet,Multipick *mp,int state,char *domain)
{
  int domain_sele = -1;
  ObjectNameType valid_name;

  UtilNCopy(valid_name, sname, sizeof(valid_name));
  if(SettingGetGlobal_b(G,cSetting_validate_object_names)) {
    ObjectMakeValidName(valid_name);
    sname = valid_name;
  }

  if(domain && domain[0]) {
    if(!WordMatchExact(G,cKeywordAll,domain,true)) { /* allow domain=all */
      domain_sele = SelectorIndexByName(G,domain);
      if(domain_sele<0) {
        
        PRINTFB(G,FB_Selector,FB_Errors)
          "Selector-Error: Invalid domain selection name \"%s\".\n",domain
          ENDFB(G);
        return -1;
      }
    }
  }
  return _SelectorCreate(G,sname,sele,&obj,quiet,mp,NULL,0,NULL,0,0,NULL,-1,state,domain_sele);
}

/*========================================================================*/
static void SelectorClean(PyMOLGlobals *G)
{
  register CSelector *I=G->Selector;
  FreeP(I->Table);
  FreeP(I->Obj);
  FreeP(I->Vertex);
  FreeP(I->Flag1);
  FreeP(I->Flag2);
  I->NAtom = 0;
}
/*========================================================================*/
static int *SelectorUpdateTableSingleObject(PyMOLGlobals *G,ObjectMolecule *obj,
                                            int req_state,
                                            int no_dummies,int *idx,
                                            int n_idx,int numbered_tags)
{
  int a=0;
  int c=0;
  int modelCnt;
  int *result = NULL;
  int tag=true;
  int state = req_state;
  register CSelector *I=G->Selector;

  PRINTFD(G,FB_Selector)
    "SelectorUpdateTableSingleObject-Debug: entered for %s...\n",obj->Obj.Name
    ENDFD;
 
  SelectorClean(G);

  switch(req_state) {
  case cSelectorUpdateTableAllStates:
    state = req_state;
    break;
  case cSelectorUpdateTableEffectiveStates:
    state = ObjectGetCurrentState(&obj->Obj,true);
    break;
  case cSelectorUpdateTableCurrentState:
    state = SceneGetState(G);
    break;
  default: 
    if(req_state<0) 
      state = cSelectorUpdateTableAllStates; /* fail safe */
    break;
  }

  switch(req_state) {
  case cSelectorUpdateTableAllStates:
    I->SeleBaseOffsetsValid = true; /* all states -> all atoms -> offsets valid */
    break;
  default:
    I->SeleBaseOffsetsValid = false; /* not including all atoms, so atom-based offsets are invalid */
    break;
  }

  I->NCSet = 0;
  if(no_dummies) {
    modelCnt = 0;
    c = 0;
  } else {
    modelCnt=cNDummyModels;
    c=cNDummyAtoms;
  }
  c+=obj->NAtom;
  if(I->NCSet<obj->NCSet) I->NCSet=obj->NCSet;
  modelCnt++;
  I->Table=Calloc(TableRec,c);
  ErrChkPtr(G,I->Table);
  I->Obj=Calloc(ObjectMolecule*,modelCnt);
  ErrChkPtr(G,I->Obj);
  if(no_dummies) {
    modelCnt = 0;
    c = 0;
  } else {
    c=cNDummyAtoms;
    modelCnt=cNDummyModels;
  }
  I->Obj[modelCnt]=obj;

  obj->SeleBase=c; 

  if(state<0) {
    for(a=0;a<obj->NAtom;a++) {
      I->Table[c].model = modelCnt;
      I->Table[c].atom = a;
      c++;
    }
  } else if(state<obj->NCSet) {
    register TableRec *rec = I->Table + c;
    CoordSet *cs = obj->CSet[state];
    if(cs) {
      for(a=0;a<obj->NAtom;a++) {
        int ix;
        if(obj->DiscreteFlag) {
          if(cs == obj->DiscreteCSet[a])
            ix=obj->DiscreteAtmToIdx[a];
          else
            ix=-1;
        } else 
          ix=cs->AtmToIdx[a];
        if(ix>=0) {
          rec->model = modelCnt;
          rec->atom = a;
          rec++;
        }
      }
    }
    c = rec - I->Table;
  }
  
  if(idx&&n_idx) {
    result = Calloc(int,c);
    if(n_idx>0) {
      for(a=0; a< n_idx; a++) {
        int at = idx[a];
        if(numbered_tags)
          tag = a+SELECTOR_BASE_TAG;
        if((at>=0)&&(at<obj->NAtom)) { /* create an ordered selection based on the input order of the object indices */
          result[obj->SeleBase + at] = tag;
        }
      }
    } else { /* -1 terminated list */
      int *at_idx = idx;
      int at;
      a=SELECTOR_BASE_TAG+1;
      while((at=*(at_idx++))>=0) {
        if(numbered_tags) {
          tag = a++;
        }
        if((at>=0)&&(at<obj->NAtom)) { /* create an ordered selection based on the input order of the object indices */
          result[obj->SeleBase + at] = tag;
        }
      }
    }
  }
  modelCnt++;
  I->NModel=modelCnt;
  I->NAtom=c;
  I->Flag1=Alloc(int,c);
  ErrChkPtr(G,I->Flag1);
  I->Flag2=Alloc(int,c);
  ErrChkPtr(G,I->Flag2);
  I->Vertex=Alloc(float,c*3);
  ErrChkPtr(G,I->Vertex);

  PRINTFD(G,FB_Selector)
    "SelectorUpdateTableSingleObject-Debug: leaving...\n"
    ENDFD;

  return(result);
}

/*========================================================================*/
int SelectorUpdateTable(PyMOLGlobals *G,int req_state,int domain)
{
  register int a=0;
  register ov_size c=0;
  register int modelCnt;
  register int state = req_state;
  void *iterator = NULL;
  ObjectMolecule *obj = NULL;

  register CSelector *I=G->Selector;

  if(!I->Origin)
    I->Origin=ObjectMoleculeDummyNew(G,cObjectMoleculeDummyOrigin);

  if(!I->Center)
    I->Center=ObjectMoleculeDummyNew(G,cObjectMoleculeDummyCenter);

  SelectorClean(G);
  I->NCSet = 0; 

  modelCnt=cNDummyModels;
  c=cNDummyAtoms;
  while(ExecutiveIterateObjectMolecule(G,&obj,&iterator)) {
    c+=obj->NAtom;
    if(I->NCSet<obj->NCSet) I->NCSet=obj->NCSet;
    modelCnt++;
  }
  I->Table=Calloc(TableRec,c);
  ErrChkPtr(G,I->Table);
  I->Obj=Calloc(ObjectMolecule*,modelCnt);
  ErrChkPtr(G,I->Obj);

  switch(req_state) {
  case cSelectorUpdateTableAllStates:
    I->SeleBaseOffsetsValid = true; /* all states -> all atoms -> offsets valid */
    break;
  default:
    I->SeleBaseOffsetsValid = false; /* not including all atoms, so atom-based offsets are invalid */
    break;
  }
  
  c=0;
  modelCnt=0;

  obj=I->Origin;
  if(obj) {
    I->Obj[modelCnt] = I->Origin;
    obj->SeleBase=c; /* make note of where this object starts */
    for(a=0;a<obj->NAtom;a++) {
      I->Table[c].model=modelCnt;
      I->Table[c].atom=a;
      c++;
    }
    modelCnt++;
  }

  obj=I->Center;
  if(obj) {
    I->Obj[modelCnt] = I->Center;
    obj->SeleBase=c; /* make note of where this object starts */
    for(a=0;a<obj->NAtom;a++) {
      I->Table[c].model=modelCnt;
      I->Table[c].atom=a;
      c++;
    }
    modelCnt++;
  }

  if(req_state<cSelectorUpdateTableAllStates) {
    state = SceneGetState(G); /* just in case... */
  }

  while(ExecutiveIterateObjectMolecule(G,&obj,&iterator)) {
    register int skip_flag = false;
    if(req_state<0) {
      switch(req_state) {
      case cSelectorUpdateTableAllStates:
        state = -1; /* all states */
        /* proceed... */
        break;
      case cSelectorUpdateTableCurrentState:
        state = SettingGetGlobal_i(G,cSetting_state)-1;
        break;
      case cSelectorUpdateTableEffectiveStates:
        state = ObjectGetCurrentState(&obj->Obj,true);
        break;
      default: /* unknown input -- fail safe (all states)*/
        state = -1;
        break;
      }
    } else {
      if(state>=obj->NCSet) 
        skip_flag=true;
      else if(!obj->CSet[state]) 
        skip_flag=true;
    }
      
    if(!skip_flag) {
      I->Obj[modelCnt]=obj;
      { 
        register int n_atom = obj->NAtom;
        register TableRec *rec = I->Table + c;
        TableRec *start_rec = rec;
        if(state<0) { /* all states */
          if(domain<0) { /* domain=all */
            for(a=0;a<n_atom;a++) {
              rec->model = modelCnt;
              rec->atom = a;
              rec++;
            } 
          } else {
            register AtomInfoType *ai = obj->AtomInfo;
            int included_one = false;
            int excluded_one = false;
            for(a=0;a<n_atom;a++) {
              if(SelectorIsMember(G,ai->selEntry,domain)) {
                rec->model = modelCnt;
                rec->atom = a;
                rec++;
                included_one = true;
              } else {
                excluded_one = true;
              }
              ai++;
            }
            if(included_one && excluded_one)
              I->SeleBaseOffsetsValid = false; /* partial objects in domain, so
                                                  base offsets are invalid */
          }
        } else { /* specific states */
          register CoordSet *cs;
          register int idx;
          if(domain<0) {
            for(a=0;a<n_atom;a++) {
              /* does coordinate exist for this atom in the requested state? */
              if(state<obj->NCSet) 
                cs=obj->CSet[state];
              else
                cs=NULL;
              if(cs) {
                if(obj->DiscreteFlag) {
                  if(cs==obj->DiscreteCSet[a])
                    idx=obj->DiscreteAtmToIdx[a];
                  else
                    idx=-1;
                } else 
                  idx=cs->AtmToIdx[a];
                if(idx>=0) {
                  rec->model = modelCnt;
                  rec->atom = a;
                  rec++;
                }
              }
            }
          } else {
            register AtomInfoType *ai = obj->AtomInfo;
            for(a=0;a<n_atom;a++) {
              /* does coordinate exist for this atom in the requested state? */
              if(state<obj->NCSet) 
                cs=obj->CSet[state];
              else
                cs=NULL;
              if(cs) {
                if(obj->DiscreteFlag) {
                  if(cs==obj->DiscreteCSet[a])
                    idx=obj->DiscreteAtmToIdx[a];
                  else
                    idx=-1;
                } else 
                  idx=cs->AtmToIdx[a];
                if(idx>=0) {
                  if(SelectorIsMember(G,ai->selEntry,domain)) {
                    rec->model = modelCnt;
                    rec->atom = a;
                    rec++;
                  }
                }
              }
              ai++;
            }
          }
        }
        if(rec!=start_rec) { /* skip excluded models */
          modelCnt++;
          obj->SeleBase=c; /* make note of where this object starts */
          c+=(rec-start_rec);
        } else {
          obj->SeleBase=0;
        }
      }
    }
  }
  I->NModel=modelCnt;
  I->NAtom=c;
  I->Flag1=Alloc(int,c);
  ErrChkPtr(G,I->Flag1);
  I->Flag2=Alloc(int,c);
  ErrChkPtr(G,I->Flag2);
  I->Vertex=Alloc(float,c*3);
  ErrChkPtr(G,I->Vertex);
  /* printf("selector update table state=%d, natom=%d\n",req_state,c); */
  return(true);
}
/*========================================================================*/
static int *SelectorSelect(PyMOLGlobals *G,char *sele,int state,int domain)
{
  SelectorWordType *parsed;
  int *result=NULL;
  PRINTFD(G,FB_Selector)
    "SelectorSelect-DEBUG: sele = \"%s\"\n",sele
    ENDFD;
  SelectorUpdateTable(G,state,domain);
  parsed=SelectorParse(G,sele);
  if(parsed) {
    if(Feedback(G,FB_Selector,FB_Debugging)) {
      SelectorWordType *a;
      fprintf(stderr,"SelectorSelect-DEBUG: parsed tokens:\n");
      a = parsed;
      while(1) {
        if(!a[0][0]) break;
        fprintf(stderr,"  \"%s\"\n",(a[0]));
        a++;
      }
      fprintf(stderr,"SelectorSelect-DEBUG: end of tokens.\n");
    }
    result=SelectorEvaluate(G,parsed,state);
    VLAFreeP(parsed);
  }
  return(result);
}
/*========================================================================*/
static int SelectorModulate1(PyMOLGlobals *G,EvalElem *base,int state)
{
  register CSelector *I=G->Selector;
  int a,d,e;
  int c=0;
  float dist;
  int nbond;
  float *v2;
  CoordSet *cs;
  int ok=true;
  int nCSet;
  MapType *map;
  int i,j,h,k,l;
  int n1,at,idx;
  ObjectMolecule *obj;

  if(state<0) {
    switch(state) {
    case -2:
    case -3:
      state=SceneGetState(G);
      break;
    }
  }

  base[1].sele=base[0].sele; /* base1 has the mask */
  base->sele=Calloc(int,I->NAtom);
  for(a=0;a<I->NAtom;a++) base[0].sele[a]=false;
  ErrChkPtr(G,base->sele);
  switch(base[1].code) {
	 case SELE_ARD_:
	 case SELE_EXP_:
		if(!sscanf(base[2].text,"%f",&dist))
		  ok=ErrMessage(G,"Selector","Invalid distance.");
		if(ok) {
          for(d=0;d<I->NCSet;d++) {
            if((state<0)||(d==state)) {
              n1=0;
              for(a=0;a<I->NAtom;a++) {
                I->Flag1[a]=false;
                at=I->Table[a].atom;
                obj=I->Obj[I->Table[a].model];
                if(d<obj->NCSet) 
                  cs=obj->CSet[d];
                else
                  cs=NULL;
                if(cs) {
                  if(obj->DiscreteFlag) {
                    if(cs==obj->DiscreteCSet[at])
                      idx=obj->DiscreteAtmToIdx[at];
                    else
                      idx=-1;
                  } else 
                    idx=cs->AtmToIdx[at];
                  if(idx>=0) {
                    copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
                    I->Flag1[a]=true;
                    n1++;
                  }
                }
              }
              if(n1) {
                map=MapNewFlagged(G,-dist,I->Vertex,I->NAtom,NULL,I->Flag1);
                if(map) {
                  MapSetupExpress(map);
                  nCSet=SelectorGetArrayNCSet(G,base[1].sele,false);
                  for(e=0;e<nCSet;e++) {
                    if((state<0)||(e==state)) {
                      for(a=0;a<I->NAtom;a++) {
                        if(base[1].sele[a]) {
                          at=I->Table[a].atom;
                          obj=I->Obj[I->Table[a].model];
                          if(e<obj->NCSet) 
                            cs=obj->CSet[e];
                          else
                            cs=NULL;
                          if(cs) {
                            if(obj->DiscreteFlag) {
                              if(cs==obj->DiscreteCSet[at])
                                idx=obj->DiscreteAtmToIdx[at];
                              else
                                idx=-1;
                            } else 
                              idx=cs->AtmToIdx[at];
                            if(idx>=0) {
                              v2 = cs->Coord+(3*idx);
                              MapLocus(map,v2,&h,&k,&l);
                              i=*(MapEStart(map,h,k,l));
                              if(i) {
                                j=map->EList[i++];
                                while(j>=0) {
                                  if((!base[0].sele[j])&&
                                     ((base[1].code==SELE_EXP_)
                                      ||(!base[1].sele[j]))) { /*exclude current selection */
                                    if(within3f(I->Vertex+3*j,v2,dist)) base[0].sele[j]=true;
                                  }
                                  j=map->EList[i++];
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  MapFree(map);
                }
              }
            }
          }
        }
		break;

	 case SELE_EXT_:
		if(sscanf(base[2].text,"%d",&nbond)!=1)
		  ok=ErrMessage(G,"Selector","Invalid bond count.");
		if(ok)
		  {
          ObjectMolecule *lastObj = NULL;
          int a,n,a0,a1,a2;
          UtilCopyMem(base[0].sele,base[1].sele,sizeof(int)*I->NAtom);
          while((nbond--)>0) {
            int *tmp = base[1].sele;
            base[1].sele=base[0].sele;
            base[0].sele=tmp;
            for(a=cNDummyAtoms;a<I->NAtom;a++) {
              if(base[1].sele[a]) {
                if(I->Obj[I->Table[a].model]!=lastObj) {
                  lastObj = I->Obj[I->Table[a].model];
                  ObjectMoleculeUpdateNeighbors(lastObj);
                }
                a0= I->Table[a].atom;
                n=lastObj->Neighbor[a0];
                n++;
                while(1) {
                  a1=lastObj->Neighbor[n];
                  if(a1<0) break;
                  if( (a2 = SelectorGetObjAtmOffset(I,lastObj,a1)) >= 0 ) {
                    base[0].sele[a2] = 1;
                    n+=2;
                  }
                }
              }
            }
          }
          FreeP(base[1].sele);
        }
      break;

	 case SELE_GAP_:
		if(!sscanf(base[2].text,"%f",&dist))
		  ok=ErrMessage(G,"Selector","Invalid distance.");
		if(ok) {
          for(a=0;a<I->NAtom;a++) {
            obj=I->Obj[I->Table[a].model];
            at=I->Table[a].atom;
            I->Table[a].f1 = obj->AtomInfo[at].vdw;
            base[0].sele[a]=true; /* start selected, subtract off */
            c=I->NAtom;
          }
          for(d=0;d<I->NCSet;d++) {
            if((state<0)||(d==state)) {
              n1=0;
              for(a=0;a<I->NAtom;a++) {
                obj=I->Obj[I->Table[a].model];
                I->Flag1[a]=false;
                at=I->Table[a].atom;
                if(d<obj->NCSet) 
                  cs=obj->CSet[d];
                else
                  cs=NULL;
                if(cs) {
                  if(obj->DiscreteFlag) {
                    if(cs==obj->DiscreteCSet[at])
                      idx=obj->DiscreteAtmToIdx[at];
                    else
                      idx=-1;
                  } else 
                    idx=cs->AtmToIdx[at];
                  if(idx>=0) {
                    copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
                    I->Flag1[a]=true;
                    n1++;
                  }
                }
              }
              if(n1) {
                map=MapNewFlagged(G,-(dist+2*MAX_VDW),I->Vertex,I->NAtom,NULL,I->Flag1);
                if(map) {

                  MapSetupExpress(map);
                  nCSet=SelectorGetArrayNCSet(G,base[1].sele,false);
                  for(e=0;e<nCSet;e++) {
                    if((state<0)||(e==state)) {
                      for(a=0;a<I->NAtom;a++) {
                        if(base[1].sele[a]) {
                          at=I->Table[a].atom;
                          obj=I->Obj[I->Table[a].model];
                          if(e<obj->NCSet) 
                            cs=obj->CSet[e];
                          else
                            cs=NULL;
                          if(cs) {
                            if(obj->DiscreteFlag) {
                              if(cs==obj->DiscreteCSet[at])
                                idx=obj->DiscreteAtmToIdx[at];
                              else
                                idx=-1;
                            } else 
                              idx=cs->AtmToIdx[at];

                            if(idx>=0) {
                              v2 = cs->Coord+(3*idx);
                              MapLocus(map,v2,&h,&k,&l);
                              i=*(MapEStart(map,h,k,l));
                              if(i) {
                                j=map->EList[i++];
                                while(j>=0) {
                                  if((base[0].sele[j])&&
                                     (!base[1].sele[j])) /*exclude current selection */
                                    {
                                      if(within3f(I->Vertex+3*j,v2,dist+ /* eliminate atoms w/o gap */
                                                  I->Table[a].f1+
                                                  I->Table[j].f1)) {
                                        base[0].sele[j]=false;
                                        c--;
                                      }
                                    }
                                  else if (base[1].sele[j]) {
                                    base[0].sele[j]=false;
                                    c--;
                                  }
                                  j=map->EList[i++];
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  MapFree(map);
                }
              }
            }
          }
        }
		break;
  }
  FreeP(base[1].sele);
  if(Feedback(G,FB_Selector,FB_Debugging)) {
	 c=0;
	 for(a=cNDummyAtoms;a<I->NAtom;a++)
		if(base[0].sele[a]) c++;
	 fprintf(stderr,"SelectorModulate0: %d atoms selected.\n",c);
  }
  return(ok);
  
}
/*========================================================================*/
static int SelectorSelect0(PyMOLGlobals *G,EvalElem *passed_base)
{
  register CSelector *I=G->Selector;
  register TableRec *i_table = I->Table;
  register ObjectMolecule **i_obj = I->Obj;
  register int a,b,flag;
  register EvalElem *base = passed_base;
  int c=0;
  signed char *vis;
  int state;
  int static_singletons;
  ObjectMolecule *obj,*cur_obj=NULL;
  CoordSet *cs;
  int at_idx;

  base->type=STYP_LIST;
  base->sele=Calloc(int,I->NAtom);
  ErrChkPtr(G,base->sele);

  switch(base->code)
	 {
    case SELE_HBAs:
    case SELE_HBDs:
	 case SELE_DONz:
	 case SELE_ACCz:
       
      { 
        /* first, verify chemistry for all atoms... */
        ObjectMolecule *lastObj=NULL,*obj;
        int at,s;
        for(a=cNDummyAtoms;a<I->NAtom;a++) {
          at=i_table[a].atom;
          obj=i_obj[i_table[a].model];
          s=obj->AtomInfo[at].selEntry;
          if(obj!=lastObj) {
            ObjectMoleculeUpdateNeighbors(obj);
            ObjectMoleculeVerifyChemistry(obj,-1);
            lastObj = obj;
          }
        }
      }
      switch(base->code) {
      case SELE_HBAs:
	 case SELE_ACCz:
        for(a=cNDummyAtoms;a<I->NAtom;a++)
          base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].hb_acceptor;
        break;
      case SELE_HBDs:
	 case SELE_DONz:
        for(a=cNDummyAtoms;a<I->NAtom;a++)
          base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].hb_donor;
        break;

      }
      break;
	 case SELE_NONz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=false;
       break;
	 case SELE_BNDz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].bonded;
       break;
	 case SELE_HETz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].hetatm;
       break;
	 case SELE_HYDz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].hydrogen;
       break;
	 case SELE_FXDz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags & cAtomFlag_fix;
       break;
	 case SELE_RSTz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags & cAtomFlag_restrain;
       break;
	 case SELE_POLz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags & cAtomFlag_polymer;
       break;
	 case SELE_SOLz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags & cAtomFlag_solvent;
       break;
	 case SELE_PTDz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].protekted;
       break;
	 case SELE_MSKz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].masked;
       break;
	 case SELE_ORGz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags & cAtomFlag_organic;
       break;
	 case SELE_INOz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags & cAtomFlag_inorganic;
       break;
	 case SELE_GIDz:
       for(a=cNDummyAtoms;a<I->NAtom;a++)
         base[0].sele[a]=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags & cAtomFlag_guide;
       break;

    case SELE_PREz:
      state = SceneGetState(G);
      static_singletons = (int)SettingGet(G,cSetting_static_singletons);
      flag=false;
      cs = NULL;
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        base[0].sele[a]=false;
        obj = i_obj[i_table[a].model];
        if(obj!=cur_obj) { /* different object */
          if(state>=obj->NCSet)
            flag=false;
          else if(!obj->CSet[state]) 
            flag=false;
          else {
            cs = obj->CSet[state];
            flag=true; /* valid state */
          }
          if(!flag)
            if(I->NCSet==1)
              if(static_singletons) {
                cs = obj->CSet[0];
                flag=true;
              }
          cur_obj = obj;
        }
        if(flag&&cs) {
          at_idx = i_table[a].atom;
          if(obj->DiscreteFlag) {
            if(cs==obj->DiscreteCSet[at_idx]) {
              if(obj->DiscreteAtmToIdx[at_idx]>=0) {
                base[0].sele[a]=true;
                c++;
              }
            } 
          } else if(cs->AtmToIdx[at_idx]>=0) {
            base[0].sele[a]=true;
            c++;
          }
        }
      }
      break;
	 case SELE_ALLz:
       for(a=cNDummyAtoms;a<I->NAtom;a++) {
         base[0].sele[a]=true;
         c++;
       }
       break;
	 case SELE_ORIz:
       for(a=0;a<I->NAtom;a++) {
         base[0].sele[a]=false;
         c++;
       }
      if(I->Origin)
        ObjectMoleculeDummyUpdate(I->Origin,cObjectMoleculeDummyOrigin);
      base[0].sele[cDummyOrigin]=true; 
      break;
	 case SELE_CENz:
		for(a=0;a<I->NAtom;a++)
		  {
			 base[0].sele[a]=false;
			 c++;
		  }
      if(I->Center)
        ObjectMoleculeDummyUpdate(I->Center,cObjectMoleculeDummyCenter);
      base[0].sele[cDummyCenter]=true; 
      break;
	 case SELE_VISz:
      {
        ObjectMolecule *last_obj = NULL;
        AtomInfoType *ai;
        int bonded;
        for(a=cNDummyAtoms;a<I->NAtom;a++)
          {
            flag = false;
            obj = i_obj[i_table[a].model];
            if(obj->Obj.Enabled) {
              ai = obj->AtomInfo + i_table[a].atom;
              vis = ai->visRep;
              bonded = ai->bonded;
              
              if(last_obj!=obj) {
                ObjectMoleculeUpdateNeighbors(obj);
                ObjectMoleculeVerifyChemistry(obj,-1);
                last_obj=obj;
              }
              
              for(b=0;b<cRepCnt;b++) {
                switch(b) {
                case cRepCartoon:                  
                case cRepRibbon:
                case cRepLine:
                case cRepCyl:
                  if(bonded&&vis[b])
                    flag=true;
                  break;
                case cRepNonbonded:
                case cRepNonbondedSphere:
                  if((!bonded)&&(vis[b]))
                    flag = true;
                  break;
                case cRepDash: /* not current applicable to atom selections */
                case cRepCell:
                case cRepCGO:
                case cRepCallback:
                case cRepExtent:
                  break;
                default:
                  if(vis[b]) {
                    flag=true;
                    break;
                  }
                }
                if(flag)
                  break;
              }
            }
            base[0].sele[a]=flag;
            if(flag)
              c++;
          }
      }
      break;
	 case SELE_ENAz:
       for(a=cNDummyAtoms;a<I->NAtom;a++) {
         flag = (i_obj[i_table[a].model]->Obj.Enabled);
         base[0].sele[a]=flag;
         if(flag)
           c++;
       }
       break;
	 }
  PRINTFD(G,FB_Selector)
	 " SelectorSelect0: %d atoms selected.\n",c
    ENDFD;

  return(1);
}
/*========================================================================*/
static int SelectorSelect1(PyMOLGlobals *G,EvalElem *base)
{
  register CSelector *I=G->Selector;
  register CWordMatcher *matcher = NULL;
  register int a,b,c=0,hit_flag;
  register ObjectMolecule **i_obj = I->Obj, *obj, *last_obj;
  register TableRec *i_table = I->Table,*table_a;
  register int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
  register int I_NAtom = I->NAtom;
  register int *base_0_sele_a;

  int model,sele,s,at_idx,col_idx;
  int flag;
  int ok=true;
  int index,state;
#if 0 
  int rmin,rmax,rtest;
  char *p;
#endif

  char *np;
  int rep_mask[cRepCnt];
  char *wildcard = SettingGetGlobal_s(G,cSetting_wildcard);

  ObjectMolecule *cur_obj = NULL;
  CoordSet *cs=NULL;

  base->type=STYP_LIST;
  base->sele=Calloc(int,I_NAtom); /* starting with zeros */
  PRINTFD(G,FB_Selector)
    " SelectorSelect1: base: %p sele: %p\n",
    (void*)base,(void*)base->sele
  ENDFD;
  ErrChkPtr(G,base->sele);
  switch(base->code)
	 {
    case SELE_PEPs:
      if(base[1].text[0]) {
        AtomInfoType *last_ai0 = NULL, *ai0;
        for(a=cNDummyAtoms;a<I_NAtom;a++) {
          ai0 = i_obj[i_table[a].model]->AtomInfo + i_table[a].atom;
          if(!AtomInfoSameResidueP(G,ai0,last_ai0)) { /* new starting residue */
            int match_found = false;
            char *ch = base[1].text; /* sequence argument */
            AtomInfoType *ai1,*last_ai1 = NULL;
            for(b=a;b<I_NAtom;b++) {
              ai1 = i_obj[i_table[b].model]->AtomInfo + i_table[b].atom; 
              if(!AtomInfoSameResidueP(G,ai1,last_ai1)) {
                if(*ch!='-') { /* if not skipping this residue */
                  if(!((*ch=='+')||(SeekerGetAbbr(G,ai1->resn,'O',0)==*ch))) { /* if a mismatch */
                    break; 
                  }
                }
                ch++;
                if(!*ch) { /* end of sequence pattern */
                  match_found = true;
                  break;
                }
                last_ai1 = ai1;
              }
            }
            if(match_found) {
              char *ch = base[1].text; /* sequence argument */
              AtomInfoType *ai1,*last_ai1 = NULL, *ai2;
              for(b=a;b<I_NAtom;b++) {
                ai1 = i_obj[i_table[b].model]->AtomInfo + i_table[b].atom;              
                if(!AtomInfoSameResidueP(G,ai1,last_ai1)) {
                  if(*ch!='-') { /* if not skipping this residue */
                    if((*ch=='+')||(SeekerGetAbbr(G,ai1->resn,'O',0)==*ch)) { /* if matched */
                      int d;
                      for(d=b;d<I_NAtom;d++) {
                        ai2 = i_obj[i_table[d].model]->AtomInfo 
                          + i_table[d].atom; /* complete residue */            
                        if(AtomInfoSameResidue(G,ai1,ai2)) {
                          c++;
                          base[0].sele[d]=true;
                        }
                      }
                    }
                  }
                  ch++;
                  if(!*ch) { /* end of sequence pattern */
                    break;
                  }
                  last_ai1 = ai1;
                }
              }
            }
          }
        }
      }
      break;
	 case SELE_IDXs:
      {
        CWordMatchOptions options;

        WordMatchOptionsConfigInteger(&options);
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true) ) ) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchInteger(matcher, 
                                          table_a->atom + 1) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }

      }
      break;
    case SELE_ID_s:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigInteger(&options);
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchInteger(matcher,
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].id) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
      break;
	 case SELE_RNKs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigInteger(&options);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchInteger(matcher,
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].rank) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
      break;
    case SELE_NAMs:
      {
        CWordMatchOptions options;
        char *atom_name_wildcard = SettingGetGlobal_s(G,cSetting_atom_name_wildcard);

        if(!atom_name_wildcard[0]) 
          atom_name_wildcard = wildcard;

        WordMatchOptionsConfigAlphaList(&options, atom_name_wildcard[0], ignore_case);

        matcher = WordMatcherNew(G,base[1].text,&options,false);

        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        last_obj = NULL;
        for(a=cNDummyAtoms;a<I_NAtom;a++) {
          obj = i_obj[table_a->model];
          if(obj!=last_obj) {

            /* allow objects to have their own atom_name_wildcards...this is a tricky workaround
             for handling nucleic acid structures that use "*" in atom names */

            char *atom_name_wildcard = SettingGet_s(G,obj->Obj.Setting,NULL,cSetting_atom_name_wildcard);
            
            if(!atom_name_wildcard[0]) 
              atom_name_wildcard = wildcard;

            if(options.wildcard != atom_name_wildcard[0]) {
              options.wildcard = atom_name_wildcard[0];
              if(matcher)
                WordMatcherFree(matcher);
              matcher=WordMatcherNew(G,base[1].text,&options,false);
              if(!matcher)
                WordPrimeCommaMatch(G,base[1].text);
            }
            last_obj = obj;
          }
          
          if(matcher) 
            hit_flag = WordMatcherMatchAlpha(matcher,i_obj[table_a->model]->AtomInfo[table_a->atom].name);
          else
            hit_flag = (WordMatchCommaExact(G,base[1].text,
                                            obj->AtomInfo[table_a->atom].name,
                                            ignore_case)<0);

          if( ( *base_0_sele_a = hit_flag ) )
            c++;
          table_a++;
          base_0_sele_a++;
        }
        if(matcher) 
          WordMatcherFree(matcher);
      }
		break;
	 case SELE_TTYs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigAlphaList(&options,wildcard[0],ignore_case); 
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          char null_st[1] = "";
          char *st;
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            AtomInfoType *ai = i_obj[table_a->model]->AtomInfo + table_a->atom;
            st = null_st;
            if(ai->textType) {
              if(!(st = OVLexicon_FetchCString(G->Lexicon,ai->textType)))
                st = null_st;
            }
            if( ( *base_0_sele_a = 
                  WordMatcherMatchAlpha(matcher,st)));
                c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
      break;
	 case SELE_ELEs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigAlphaList(&options,wildcard[0],ignore_case);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchAlpha(matcher,
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].elem) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
		break;
	 case SELE_SEGs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigAlphaList(&options,wildcard[0],ignore_case);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchAlpha(matcher,
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].segi) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
      break;
	 case SELE_REPs:
       for(a=0;a<cRepCnt;a++)
         rep_mask[a]=false;
       WordPrimeCommaMatch(G,base[1].text);
       a=0;
       while(1) {
         if(!rep_names[a].word[0]) break;
         if(WordMatchComma(G,base[1].text,
                           rep_names[a].word,
                           ignore_case)<0)
           rep_mask[rep_names[a].value]=true;
         a++;
       }
       for(a=cNDummyAtoms;a<I_NAtom;a++)
         {
           base[0].sele[a]=false;
           for(b=0;b<cRepCnt;b++) {
             if(rep_mask[b]&&i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].visRep[b]) {
               base[0].sele[a]=true;
               c++;
               break;
             }
           }
         }
       break;
	 case SELE_COLs:
       col_idx = ColorGetIndex(G,base[1].text);
       for(a=cNDummyAtoms;a<I_NAtom;a++) {
         base[0].sele[a]=false;
         if(i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].color==col_idx) {
           base[0].sele[a]=true;
           c++;
         }
       }
       break;
	 case SELE_CCLs:
       col_idx = ColorGetIndex(G,base[1].text);
       for(a=cNDummyAtoms;a<I_NAtom;a++) {
         base[0].sele[a]=false;
         {
             AtomInfoType *ai = i_obj[i_table[a].model]->AtomInfo + i_table[a].atom;
         if(ai->has_setting) {
           int value; 
           if(SettingUniqueGet_color(G,ai->unique_id, cSetting_cartoon_color,&value)) {
             if(value == col_idx) {
               base[0].sele[a]=true;
               c++;
             }
           }
         }
         }
       }
       break;
     case SELE_RCLs:
       col_idx = ColorGetIndex(G,base[1].text);
       for(a=cNDummyAtoms;a<I_NAtom;a++) {
         base[0].sele[a]=false;
         {
             AtomInfoType *ai = i_obj[i_table[a].model]->AtomInfo + i_table[a].atom;
         if(ai->has_setting) {
           int value; 
           if(SettingUniqueGet_color(G,ai->unique_id, cSetting_ribbon_color,&value)) {
             if(value == col_idx) {
               base[0].sele[a]=true;
               c++;
             }
           }
         }
         }
       }
       break;
	 case SELE_CHNs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigAlphaList(&options,wildcard[0],ignore_case);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchAlpha(matcher,
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].chain) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
		break;
	 case SELE_SSTs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigAlphaList(&options,wildcard[0],ignore_case);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchAlpha(matcher,
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].ssType) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
		break;
	 case SELE_STAs:
      sscanf(base[1].text,"%d",&state);
      state = state - 1;
      obj = NULL;
      
      if(state<0) {
        for(a=cNDummyAtoms;a<I_NAtom;a++)
          base[0].sele[a]=false;
      } else {
        for(a=cNDummyAtoms;a<I_NAtom;a++)
          {
            base[0].sele[a]=false;
            obj = i_obj[i_table[a].model];
            if(obj!=cur_obj) { /* different object */
              if(state>=obj->NCSet)
                flag=false;
              else if(!obj->CSet[state]) 
                flag=false;
              else {
                cs = obj->CSet[state];
                flag=true; /* valid state */
              }
              cur_obj = obj;
            }
            if(flag&&cs) {
              at_idx = i_table[a].atom;
              if(obj->DiscreteFlag) {
                if(cs==obj->DiscreteCSet[at_idx]) {
                  if(obj->DiscreteAtmToIdx[at_idx]>=0) {
                    base[0].sele[a]=true;
                    c++;
                  }
                }
              } else if(cs->AtmToIdx[at_idx]>=0) {
                base[0].sele[a]=true;
                c++;
              }
            }
          }
      }
		break;
	 case SELE_ALTs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigAlphaList(&options,wildcard[0],ignore_case);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchAlpha(matcher,
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].alt) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
		break;
	 case SELE_FLGs:
      sscanf(base[1].text,"%d",&flag);
      flag = (1<<flag);
		for(a=cNDummyAtoms;a<I_NAtom;a++)
		  {
			 if(i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].flags&flag)
				{
				  base[0].sele[a]=true;
				  c++;
				}
			 else
				base[0].sele[a]=false;
        }
		break;
	 case SELE_NTYs:
      {
        CWordMatchOptions options;

        WordMatchOptionsConfigInteger(&options);
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true) ) ) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchInteger(matcher, 
                                          i_obj[table_a->model]->AtomInfo[table_a->atom].customType ) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
      break;
	 case SELE_RSIs:
      {
        CWordMatchOptions options;
        AtomInfoType *ai;

        WordMatchOptionsConfigMixed(&options,wildcard[0],ignore_case);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            ai = i_obj[table_a->model]->AtomInfo + table_a->atom;
            if( ( *base_0_sele_a = 
                  WordMatcherMatchMixed(matcher,
                                        ai->resi, ai->resv) ) )
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
      break;
#if 0
      /* delete this once you're sure the new code works! */

		if((p=strstr(base[1].text,":"))/* range */
         ||(p=strstr(base[1].text,"-")))/* range */
		  {
			 *p=' ';
			 if(sscanf(base[1].text,"%i%i",&rmin,&rmax)!=2)
				ok=ErrMessage(G,"Selector","Invalid Range.");
			 if(ok)
				for(a=cNDummyAtoms;a<I_NAtom;a++)
				  {
					 if(sscanf(i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].resi,
								  "%i",&rtest))
						{
						  if((rtest>=rmin)&&(rtest<=rmax)) {
							 base[0].sele[a]=true;
                      c++;
						  } else 
							 base[0].sele[a]=false;
						}
					 else 
						base[0].sele[a]=false;
				  }
		  }
		else /* not a range */ {
        WordPrimeCommaMatch(G,base[1].text);
		  for(a=cNDummyAtoms;a<I_NAtom;a++)
			 {
				if(WordMatchComma(G,base[1].text,
                         i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].resi,
                         ignore_case)<0)
				  {
					 base[0].sele[a]=true;
					 c++;
				  }
				else
				  base[0].sele[a]=false;
			 }
      }
		break;
#endif

	 case SELE_RSNs:
      {
        CWordMatchOptions options;
        
        WordMatchOptionsConfigAlphaList(&options,wildcard[0],ignore_case);
        
        table_a = i_table + cNDummyAtoms;
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        
        if( (matcher = WordMatcherNew(G,base[1].text,&options,true))) {
          table_a = i_table + cNDummyAtoms;
          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          
          for(a=cNDummyAtoms;a<I_NAtom;a++) {
            if( ( *base_0_sele_a = 
                  WordMatcherMatchAlpha(matcher,
                                        i_obj[table_a->model]->AtomInfo[table_a->atom].resn )))
              c++;
            table_a++;
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        }
      }
		break;
	 case SELE_SELs:
      {
        char *word = base[1].text;
        int enabled_only = false;
        CWordMatchOptions options;
        
        if(word[0]=='?') {
          word++;
          if(word[0]=='?') {
            enabled_only = true;
            word++;
          }
        }

        WordMatchOptionsConfigAlpha(&options,wildcard[0],ignore_case);
        
        if( (matcher = WordMatcherNew(G,word,&options,false))) {

          SelectorWordType *list = I->Name;
          int idx = 0;

          for(a=0;a<I_NAtom;a++) /* zero out first before iterating through selections */
            base[0].sele[a]=false;
          
          while(list[idx][0]) {
            if(WordMatcherMatchAlpha(matcher,list[idx])) {
              if((idx>=0)&&
                 ((!enabled_only)||
                  ExecutiveGetActiveSeleName(G,list[idx],false,false))) {
                register MemberType *I_Member = I->Member;
                sele=I->Info[idx].ID;
                for(a=cNDummyAtoms;a<I_NAtom;a++) {
                  s=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].selEntry;
                  while(s) {
                    if(I_Member[s].selection==sele) {
                      if(!base[0].sele[a]) {
                        base[0].sele[a] = I_Member[s].tag;
                        c++;
                      }
                    }
                    s=I_Member[s].next;
                  }
                }
              }
            }
            idx++;
          }
          WordMatcherFree(matcher);

          /* must also allow for group name pattern matches */
          
          {
            int group_list_id;
            if( (group_list_id = ExecutiveGetExpandedGroupListFromPattern(G,word))) {
              int last_was_member = false;
              last_obj = NULL;
              for(a=cNDummyAtoms;a<I_NAtom;a++) {
                if(last_obj!=i_obj[i_table[a].model]) {
                  last_obj = i_obj[i_table[a].model];
                  last_was_member = ExecutiveCheckGroupMembership(G,
                                                                  group_list_id,
                                                                  (CObject*)last_obj);
                }
                if(last_was_member && !base[0].sele[a]) {
                  base[0].sele[a] = true;
                  c++;
                }
              }
            }
            ExecutiveFreeGroupList(G,group_list_id);
          }

        } else if((!enabled_only)|| ExecutiveGetActiveSeleName(G,word,false,false)) {
          sele=SelectGetNameOffset(G,word,1,ignore_case);
          if(sele>=0) {
            register MemberType *I_Member = I->Member;
            sele=I->Info[sele].ID;
            for(a=cNDummyAtoms;a<I_NAtom;a++) {
              base[0].sele[a]=false;
              s=i_obj[i_table[a].model]->AtomInfo[i_table[a].atom].selEntry;
              while(s) {
                if(I_Member[s].selection == sele) {
                  base[0].sele[a] = I_Member[s].tag;
                  c++;
                }
                s=I_Member[s].next;
              }
            }
          } else {
            int group_list_id;
            if(  (group_list_id = ExecutiveGetExpandedGroupList(G,word)) ) {
              int last_was_member = false;
              last_obj = NULL;
              for(a=0;a<I_NAtom;a++) /* zero out first before iterating through selections */
                base[0].sele[a]=false;
              for(a=cNDummyAtoms;a<I_NAtom;a++) {
                if(last_obj!=i_obj[i_table[a].model]) {
                  last_obj = i_obj[i_table[a].model];
                  last_was_member = ExecutiveCheckGroupMembership(G,
                                                                  group_list_id,
                                                                  (CObject*)last_obj);
                 }
                if( (base[0].sele[a] = last_was_member)) c++;
              }
              ExecutiveFreeGroupList(G,group_list_id);
            } else if(base[1].text[0]=='?') { /* undefined ?sele allowed */
              for(a=cNDummyAtoms;a<I_NAtom;a++)
                base[0].sele[a]=false;
            } else {
              PRINTFB(G,FB_Selector,FB_Errors)
                "Selector-Error: Invalid selection name \"%s\".\n",word
                ENDFB(G);
              ok = false;
            }
          }
        }
      }
      break;
	 case SELE_MODs:

       /* need to change this to handle wildcarded model names */
       
       /* first, trim off and record the atom index if one exists */
       
       index = -1;
       if((np=strstr(base[1].text,"`"))) {
         *np=0;
         if(sscanf(np+1,"%d",&index)!=1)
           index = -1;
         else
           index--;
       }
       model=0;
       
       {
         CWordMatchOptions options;
         WordMatchOptionsConfigAlpha(&options,wildcard[0],ignore_case);

         if( (matcher = WordMatcherNew(G,base[1].text,&options,false))) {
           
           int obj_matches = false;
           
           for(a=0;a<I_NAtom;a++) /* zero out first before iterating through selections */
             base[0].sele[a]=false;
           
           table_a = i_table + cNDummyAtoms;
           base_0_sele_a = &base[0].sele[cNDummyAtoms];
           last_obj = NULL;
           for(a=cNDummyAtoms;a<I_NAtom;a++) {
             obj = i_obj[table_a->model];
             if(obj!=last_obj) {
               
               obj_matches = WordMatcherMatchAlpha(matcher,
                                                   i_obj[table_a->model]->Obj.Name);
               last_obj = obj;
             }
             if(obj_matches) {
               if((index<0) || (table_a->atom==index)) {
                 *base_0_sele_a = true;
                 c++;
               }
             }
             table_a++;
             base_0_sele_a++;
           }
           WordMatcherFree(matcher);
         } else {
           
           obj=(ObjectMolecule*)ExecutiveFindObjectByName(G,base[1].text);
           if(obj)
             {
               for(a=cNDummyModels;a<I->NModel;a++)
                 if(i_obj[a]==obj)
                   {
                     model=a+1;
                     break;
                   }
             }
           if(!model) 
             if(sscanf(base[1].text,"%i",&model)==1)
               {
                 if(model<=0)
                   model=0;
                 else if(model>I->NModel)
                   model=0;
                 else if(!i_obj[model])
                   model=0;
               }
           if(model)
             {
               model--;
               if(index>=0) {
                 for(a=cNDummyAtoms;a<I_NAtom;a++)
                   {
                     if(i_table[a].model==model)
                       if(i_table[a].atom==index)
                         {
                           base[0].sele[a]=true;
                           c++;
                         }
                       else {
                         base[0].sele[a]=false;                    
                       }
                     else
                       base[0].sele[a]=false;
                   }
               } else {
                 for(a=cNDummyAtoms;a<I_NAtom;a++)
                   {
                     if(i_table[a].model==model)
                       {
                         base[0].sele[a]=true;
                         c++;
                       }
                     else
                       base[0].sele[a]=false;
                   }
               }
             }
           else {
             PRINTFB(G,FB_Selector,FB_Errors)
               " Selector-Error: invalid model \"%s\".\n",base[1].text
               ENDFB(G);
             ok=false;
           }
         }
       }
       break;
	 }
  PRINTFD(G,FB_Selector)
	 " SelectorSelect1:  %d atoms selected.\n",c
    ENDFD;
  return(ok);
}
/*========================================================================*/
static int SelectorSelect2(PyMOLGlobals *G,EvalElem *base)
{
  int a;
  int c=0;
  int ok=true;
  int oper;
  float comp1;
  int exact;
  register int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);

  AtomInfoType *at1;
  register CSelector *I=G->Selector;
  base->type=STYP_LIST;
  base->sele=Calloc(int,I->NAtom);
  ErrChkPtr(G,base->sele);
  switch(base->code)
	 {
	 case SELE_PCHx:
	 case SELE_FCHx:
	 case SELE_BVLx:
	 case SELE_QVLx:
      oper=WordKey(G,AtOper,base[1].text,4,ignore_case,&exact);
      if(!oper)
        ok=ErrMessage(G,"Selector","Invalid Operator.");
      if(ok) {
        switch(oper) {
        case SCMP_GTHN:
        case SCMP_LTHN:
        case SCMP_EQAL:
          if (sscanf(base[2].text,"%f",&comp1)!=1) 
            ok=ErrMessage(G,"Selector","Invalid Number");
          break;
        }
        if(ok) {
          switch(oper) {
          case SCMP_GTHN:
            switch(base->code) {
            case SELE_BVLx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->b>comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_QVLx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->q>comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_PCHx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->partialCharge>comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_FCHx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->formalCharge>comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            }
            break;
          case SCMP_LTHN:
            switch(base->code) {
            case SELE_BVLx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->b<comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;               
            case SELE_QVLx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->q<comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;               
            case SELE_PCHx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->partialCharge<comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_FCHx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(at1->formalCharge<comp1) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            }
            break;
          case SCMP_EQAL:
            switch(base->code) {
            case SELE_BVLx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->b-comp1)<R_SMALL4) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_QVLx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->q-comp1)<R_SMALL4) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_PCHx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->partialCharge-comp1)<R_SMALL4) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            case SELE_FCHx:
              for(a=cNDummyAtoms;a<I->NAtom;a++) {
                at1=&I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom];
                if(fabs(at1->formalCharge-comp1)<R_SMALL4) {
                  base[0].sele[a]=true;
                  c++;
                } else {
                  base[0].sele[a]=false;
                }
              }
              break;
            }
            break;
          }
          break;
        }
      }
    }
  
  PRINTFD(G,FB_Selector)
	 " SelectorSelect2: %d atoms selected.\n",c
    ENDFD;
  return(ok);
}
/*========================================================================*/
static int SelectorLogic1(PyMOLGlobals *G,EvalElem *inp_base,int state)
{
  /* some cases in this function still need to be optimized
     for performance (see BYR1 for example) */

  register CSelector *I=G->Selector;
  register int a,b,tag;
  register int c=0;
  register int flag;
  register EvalElem *base = inp_base;
  register AtomInfoType *at1,*at2;
  register TableRec *i_table = I->Table, *table_a;
  register ObjectMolecule **i_obj = I->Obj;
  register int n_atom = I->NAtom;
  register int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
  int n;
  int a0,a1,a2;
  ObjectMolecule *lastObj = NULL;

  base[0].sele=base[1].sele;
  base[1].sele=NULL;
  base[0].type=STYP_LIST;
  switch(base->code) {
  case SELE_NOT1:
    {
      register int *base_0_sele_a;
      
      base_0_sele_a = base[0].sele;
      for(a=0;a<n_atom;a++) {
        if( (*base_0_sele_a = ! *base_0_sele_a) )
          c++;
        base_0_sele_a ++;
      }
    }
    break;
  case SELE_NGH1:
    base[1].sele=base[0].sele;
    base[0].sele=Calloc(int,n_atom);
    
    table_a = i_table + cNDummyAtoms;
    for(a=cNDummyAtoms;a<n_atom;a++) {
      if( (tag=base[1].sele[a]) ) {
        if(i_obj[ table_a->model]!=lastObj) {
          lastObj = i_obj[table_a->model];
          ObjectMoleculeUpdateNeighbors(lastObj);
        }
        a0 = table_a->atom;
        n = lastObj->Neighbor[a0];
        n++;
        while(1) {
          a1 = lastObj->Neighbor[n];
          if(a1<0) break;
          if( (a2 = SelectorGetObjAtmOffset(I,lastObj,a1)) >= 0 ) {
            if(!base[1].sele[a2])
              base[0].sele[a2] = tag;
          }
          n+=2;
        }
      }
      table_a++;
    }
    FreeP(base[1].sele);
    break;
  case SELE_BON1:
    base[1].sele=base[0].sele;
    base[0].sele=Calloc(int,n_atom);
    table_a = i_table + cNDummyAtoms;
    for(a=cNDummyAtoms;a<n_atom;a++) {
      if( (tag=base[1].sele[a]) ) {
        if(i_obj[table_a->model]!=lastObj) {
          lastObj = i_obj[table_a->model];
          ObjectMoleculeUpdateNeighbors(lastObj);
        }
        a0= table_a->atom;
        n=lastObj->Neighbor[a0];
        n++;
        while(1) {
          a1=lastObj->Neighbor[n];
          if(a1<0) break;
          if( (a2 = SelectorGetObjAtmOffset(I,lastObj,a1)) >= 0 ) {
            if(!base[0].sele[a2]) 
              base[0].sele[a2]=1;
            n+=2;
          }
        }
      }
      table_a++;
    }
    FreeP(base[1].sele);
    break;
  case SELE_BYO1: 
    base[1].sele=base[0].sele;
    base[0].sele=Calloc(int,n_atom);
    for(a=cNDummyAtoms;a<n_atom;a++) {
      if(base[1].sele[a]) {
        if(i_obj[i_table[a].model]!=lastObj) {
          lastObj = i_obj[i_table[a].model];
          b = a;
          while(b>=0) {
            if(i_obj[i_table[b].model]!=lastObj) 
              break;
            base[0].sele[b]=1;
            b--;
          }
          b=a+1;
          while(b<n_atom) {
            if(i_obj[i_table[b].model]!=lastObj) 
              break;
            base[0].sele[b]=1;
            b++;
          }
        }
      }
    }
    FreeP(base[1].sele);
    break;      
  case SELE_BYR1: /* ASSUMES atoms are sorted & grouped by residue */
  case SELE_CAS1: 
    {
      int *base_0_sele = base[0].sele;
      int break_atom = -1; 
      int last_tag = 0;
      table_a = i_table + cNDummyAtoms;
      for(a=cNDummyAtoms;a<n_atom;a++) {
        if( (tag = base_0_sele[a]) && ((a >= break_atom)|| (base_0_sele[a] != last_tag))) {
          at1=&i_obj[table_a->model]->AtomInfo[table_a->atom];
          b = a-1;
          while(b>=0) {
            if(!base_0_sele[b]) {
              flag = false;
              if(table_a->model==i_table[b].model) {
                at2=&i_obj[i_table[b].model]->AtomInfo[i_table[b].atom];
                if(at1->chain[0]==at2->chain[0])
                  if(at1->resv==at2->resv)
                    if(at1->discrete_state==at2->discrete_state)
                      if(WordMatch(G,at1->resi,at2->resi,ignore_case)<0)
                        if(WordMatch(G,at1->resn,at2->resn,ignore_case)<0)
                          if(WordMatch(G,at1->segi,at2->segi,ignore_case)<0) {
                            base_0_sele[b]=tag;
                            c++;
                            flag=1;
                          }
              }
              if(!flag) {
                break;
              }
            }
            b--;
          }
          b = a + 1;
          while(b<n_atom) {
            if(!base_0_sele[b]) {
              flag=false;
              if(table_a->model==i_table[b].model) {
                at2=&i_obj[i_table[b].model]->AtomInfo[i_table[b].atom];
                if(at1->chain[0]==at2->chain[0])
                  if(at1->resv==at2->resv)
                    if(at1->discrete_state==at2->discrete_state)
                      if(WordMatch(G,at1->resi,at2->resi,ignore_case)<0)
                        if(WordMatch(G,at1->resn,at2->resn,ignore_case)<0)
                          if(WordMatch(G,at1->segi,at2->segi,ignore_case)<0) {
                            base_0_sele[b]=tag;
                            c++;
                            flag=1;
                          }
              }
              if(!flag) {
                break_atom = b-1;
                last_tag = tag;
                break;
              }
            }
            b++;
          }
        }
        table_a++;
      }
      if(base->code==SELE_CAS1) {
        c=0;
        table_a = i_table + cNDummyAtoms;
        for(a=cNDummyAtoms;a<n_atom;a++) {
          if(base_0_sele[a]) {
            base_0_sele[a]=false;
            
            if(i_obj[table_a->model]->AtomInfo[table_a->atom].protons == cAN_C)
              if(WordMatchCommaExact(G,"CA",
                                     i_obj[table_a->model]->AtomInfo[table_a->atom].name,
                                     ignore_case)<0) {
                base_0_sele[a]=true;
                c++;
              }
          }
          table_a++;
        }
      }
    }
    break;
  case SELE_BYC1: /* ASSUMES atoms are sorted & grouped by chain */
    { 
      int *base_0_sele = base[0].sele;
      int break_atom_high = -1; 
      int break_atom_low = 0;
      int last_tag = 0;
      table_a = i_table + cNDummyAtoms;
      for(a=cNDummyAtoms;a<n_atom;a++) {
        if( (tag = base_0_sele[a]) && ((a >= break_atom_high) || (base_0_sele[a] != last_tag))) {
          if(tag!=last_tag)
            break_atom_low = 0;
          at1=&i_obj[table_a->model]->AtomInfo[table_a->atom];
          b = a-1;
          while(b>=break_atom_low) {
            if(!base_0_sele[b]) {
              flag = false;
              if(table_a->model==i_table[b].model)
                {
                  at2=&i_obj[i_table[b].model]->AtomInfo[i_table[b].atom];
                  if(at1->chain[0]==at2->chain[0])
                    if(WordMatch(G,at1->segi,at2->segi,ignore_case)<0) {
                      base_0_sele[b]=tag;
                      c++;
                      flag=1;
                    }
                }
              if(!flag) {
                break_atom_low = b+1;
                break;
              }
            }
            b--;
          }
          if(b<0) break_atom_low = -1;
          b = a + 1;
          while(b<n_atom) {
            if(!base_0_sele[b]) {
              flag=false;
              if(table_a->model==i_table[b].model) {
                at2=&i_obj[i_table[b].model]->AtomInfo[i_table[b].atom];
                if(at1->chain[0]==at2->chain[0])
                  if(WordMatch(G,at1->segi,at2->segi,ignore_case)<0) {
                    base_0_sele[b]=tag;
                    c++;
                    flag=1;
                  }
              }
              if(!flag) {
                break_atom_high = b-1;
                last_tag = tag;
                break;
              }
            }
            b++;
          }
        }
        table_a++;
      }
    }
    break;
  case SELE_BYS1: /* ASSUMES atoms are sorted & grouped by segi */
    { 
      int *base_0_sele = base[0].sele;
      int break_atom_high = -1; 
      int break_atom_low = 0;
      int last_tag = 0;
      table_a = i_table + cNDummyAtoms;
      for(a=cNDummyAtoms;a<n_atom;a++) {
        if( (tag = base_0_sele[a]) && ((a >= break_atom_high) || (base_0_sele[a] != last_tag))) {
          if(tag!=last_tag)
            break_atom_low = 0;
          at1=&i_obj[table_a->model]->AtomInfo[table_a->atom];
          b = a-1;
          while(b>=break_atom_low) {
            if(!base_0_sele[b]) {
              flag = false;
              if(table_a->model==i_table[b].model) {
                at2=&i_obj[i_table[b].model]->AtomInfo[i_table[b].atom];
                if(WordMatch(G,at1->segi,at2->segi,ignore_case)<0) {
                  base_0_sele[b]=tag;
                  c++;
                  flag=1;
                }
              }
              if(!flag) {
                break_atom_low = b+1;
                break;
              }
            }
            b--;
          }
          b = a + 1;
          while(b<n_atom) {
            if(!base_0_sele[b]) {
              flag=false;
              if(table_a->model==i_table[b].model) {
                at2=&i_obj[i_table[b].model]->AtomInfo[i_table[b].atom];
                if(WordMatch(G,at1->segi,at2->segi,ignore_case)<0) {
                  base_0_sele[b]=tag;
                  c++;
                  flag=1;
                }
              }
              if(!flag) {
                break_atom_high = b-1;
                last_tag = tag;
                break;
              }
            }
            b++;
          }
        }
        table_a++;
      }
    }
    break;
  case SELE_BYF1: /* first, identify all atom by fragment selection */
    {
      /* NOTE: this algorithm looks incompatible with selection
         tags...need to do some more thinking & work... */

      int n_frag = EditorGetNFrag(G);

      base[1].sele=base[0].sele;
      base[0].sele=Calloc(int,n_atom);
        
      if(n_frag) {
        int a,f,at,s;
        int *fsele;
        WordType name;
        ObjectMolecule *obj;

        fsele = Alloc(int,n_frag+1);
          
        for(f=0;f<n_frag;f++) {
          sprintf(name,"%s%1d",cEditorFragPref,f+1);
          fsele[f] = SelectorIndexByName(G,name);
        }
          
        /* mark atoms by fragment */
        for(a=0;a<n_atom;a++) {
          at=i_table[a].atom;
          obj=i_obj[i_table[a].model];
          s=obj->AtomInfo[at].selEntry;
          for(f=0;f<n_frag;f++) {
            if(SelectorIsMember(G,s,fsele[f])) {
              base[0].sele[a] = f+1;
            }
          }
        }

        /* mark fragments we keep */
        for(f=0;f<=n_frag;f++) {
          fsele[f] = 0;
        }
        for(a=0;a<n_atom;a++) { 
          int f = base[0].sele[a];
          if(base[1].sele[a]&&f)
            fsele[f] = 1;
        }
          
        /* now set flags */
        for(a=0;a<n_atom;a++) {
          c+= (base[0].sele[a] = fsele[base[0].sele[a]]);
        }

        FreeP(fsele);
      }
      FreeP(base[1].sele);
    }
    break;
  case SELE_BYM1:
    {
      int s;
      int c = 0;
      int a,at,a1,aa;
      AtomInfoType *ai;
      ObjectMolecule *obj,*lastObj=NULL;
      int *stk;
      int stkDepth = 0;
      base[1].sele=base[0].sele;
      base[0].sele=Calloc(int,n_atom);
        
      stk = VLAlloc(int,50);
        
      for(a=0;a<n_atom;a++) {
        if( (tag=base[1].sele[a]) &&(!base[0].sele[a])) {
          VLACheck(stk,int,stkDepth);
          stk[stkDepth]=a;
          stkDepth++;

          obj=i_obj[i_table[a].model];
          if(obj!=lastObj) {
            lastObj = obj;
            ObjectMoleculeUpdateNeighbors(obj);
          }
            
          while(stkDepth) { /* this will explore a tree */
            stkDepth--;
            a=stk[stkDepth];
            base[0].sele[a]=tag;
            c++;
            at=i_table[a].atom; /* start walk from this location */
            ai=obj->AtomInfo+at;

            s=obj->Neighbor[at]; /* add neighbors onto the stack */
            s++; /* skip count */              
            while(1) {
              a1 = obj->Neighbor[s];
              if(a1>=0) {
                if( (aa = SelectorGetObjAtmOffset(I,obj,a1)) >= 0 ) {
                  if(!base[0].sele[aa]) {
                    VLACheck(stk,int,stkDepth);
                    stk[stkDepth]=aa; /* add index in selector space */
                    stkDepth++;
                  }
                }
              } else 
                break;
              s+=2;
            }
          }
        }
      }
      FreeP(base[1].sele);
      VLAFreeP(stk);
    }
    break;
  case SELE_BYX1: /* by cell */
    base[1].sele=base[0].sele;
    base[0].sele=Calloc(int,n_atom);
    {		
      ObjectMolecule *obj;
      CoordSet *cs;
      int d,n1,at;
      for(d=0;d<I->NCSet;d++) {
        if((state<0)||(d==state)) {
          n1=0;
          for(a=0;a<I->NAtom;a++) {
            I->Flag1[a]=false;
            at=I->Table[a].atom;
            obj=I->Obj[I->Table[a].model];
            if(d<obj->NCSet) 
              cs=obj->CSet[d];
            else
              cs=NULL;
            if(cs) {
              CCrystal *cryst = cs->PeriodicBox;
              if((!cryst) && (obj->Symmetry))
                cryst = obj->Symmetry->Crystal;
              if(cryst) {
                int idx;
                if(obj->DiscreteFlag) {
                  if(cs == obj->DiscreteCSet[at])
                    idx = obj->DiscreteAtmToIdx[at];
                  else
                    idx = -1;
                } else 
                  idx = cs->AtmToIdx[at];
                if(idx >= 0) {
                  transform33f3f(cryst->RealToFrac, cs->Coord+(3*idx), I->Vertex+3*a);
                  I->Flag1[a]=true;
                  n1++;
                }
              }
            }
          }
          if(n1) {
            MapType *map=MapNewFlagged(G,-1.1,I->Vertex,I->NAtom,NULL,I->Flag1);
            if(map) {
              int e, nCSet;
              MapSetupExpress(map);
              nCSet=SelectorGetArrayNCSet(G,base[1].sele,false);
              for(e=0;e<nCSet;e++) {
                if((state<0)||(e==state)) {
                  for(a=0;a<I->NAtom;a++) {
                    if(base[1].sele[a]) {
                      at=I->Table[a].atom;
                      obj=I->Obj[I->Table[a].model];
                      if(e<obj->NCSet) 
                        cs=obj->CSet[e];
                      else
                        cs=NULL;
                      if(cs) {
                        CCrystal *cryst = cs->PeriodicBox;
                        if((!cryst) && (obj->Symmetry))
                          cryst = obj->Symmetry->Crystal;
                        if(cryst) {
                          int idx;
                          if(obj->DiscreteFlag) {
                            if(cs==obj->DiscreteCSet[at])
                              idx=obj->DiscreteAtmToIdx[at];
                            else
                              idx=-1;
                          } else 
                            idx=cs->AtmToIdx[at];
                          if(idx>=0) {
                            float probe[3], probe_i[3];
                            int h,i,j,k,l;
                            
                            transform33f3f(cryst->RealToFrac, cs->Coord+(3*idx), probe);
                            MapLocus(map,probe,&h,&k,&l);
                            i=*(MapEStart(map,h,k,l));
                            if(i) {
                              j=map->EList[i++];
                                  
                              probe_i[0] = (int)floor(probe[0]);
                              probe_i[1] = (int)floor(probe[1]);
                              probe_i[2] = (int)floor(probe[2]);
                              
                              while(j>=0) {
                                if( !base[0].sele[j] ) { 
                                  float *tst = I->Vertex + 3*j;
                                  base[0].sele[j] = ((probe_i[0] == (int)floor(tst[0])) &&
                                                     (probe_i[1] == (int)floor(tst[1])) &&
                                                     (probe_i[2] == (int)floor(tst[2])));
                                }
                                j=map->EList[i++];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              MapFree(map);
            }
          }
        }
      }
    }
    FreeP(base[1].sele);
    break;
  case SELE_FST1: 
    base[1].sele=base[0].sele;
    base[0].sele=Calloc(int,n_atom);
    for(a=cNDummyAtoms;a<n_atom;a++) {
      if(base[1].sele[a]) {
        base[0].sele[a] = base[1].sele[a]; /* preserve tag */
        break;
      }
    }
    FreeP(base[1].sele);
    break;      
  case SELE_LST1: 
    {
      int last = -1;
      base[1].sele=base[0].sele;
      base[0].sele=Calloc(int,n_atom);
      for(a=cNDummyAtoms;a<n_atom;a++) {
        if(base[1].sele[a]) {
          last = a;
        }
      }
      if(last>=0)
        base[0].sele[last] = base[1].sele[last]; /* preserve tag */
    }
    FreeP(base[1].sele);
    break;      	
  }
  PRINTFD(G,FB_Selector)
    " SelectorLogic1: %d atoms selected.\n",c
    ENDFD;
  return(1);
}
/*========================================================================*/
static int SelectorLogic2(PyMOLGlobals *G,EvalElem *base)
{
  register CSelector *I=G->Selector;
  register int a,b,tag;
  register int c=0;
  register int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
  register int *base_0_sele_a,*base_2_sele_a;
  register TableRec *i_table = I->Table, *table_a, *table_b;
  register ObjectMolecule **i_obj = I->Obj;
  register int n_atom = I->NAtom;

  AtomInfoType *at1,*at2;
  
  switch(base[1].code) {
    
  case SELE_OR_2:
  case SELE_IOR2:
    {
      base_0_sele_a = base[0].sele;
      base_2_sele_a = base[2].sele;
      
      for(a=0;a<n_atom;a++) {
        if( ((*base_0_sele_a)= (((*base_0_sele_a) > (*base_2_sele_a)) ?
                                (*base_0_sele_a) : (*base_2_sele_a))) ) /* use higher tag */
          c++;
        base_0_sele_a++;
        base_2_sele_a++;
      }
    }
    break;
  case SELE_AND2:
    
    base_0_sele_a = base[0].sele;
    base_2_sele_a = base[2].sele;
    
    for(a=0;a<n_atom;a++) {
      if( (*base_0_sele_a) && (*base_2_sele_a) ) {
        (*base_0_sele_a) = (((*base_0_sele_a) > (*base_2_sele_a)) ?
                            (*base_0_sele_a) : (*base_2_sele_a)); /* use higher tag */
        c++;
      } else {
        (*base_0_sele_a) = 0;
      }
      base_0_sele_a++;
      base_2_sele_a++;
    }
    break;
  case SELE_ANT2:
    base_0_sele_a = base[0].sele;
    base_2_sele_a = base[2].sele;
    
    for(a=0;a<n_atom;a++) {
      if( (*base_0_sele_a) && ! (*base_2_sele_a) ) {
        c++;
      } else {
        (*base_0_sele_a) = 0;
      }
      base_0_sele_a++;
      base_2_sele_a++;
    }
    break;
  case SELE_IN_2:
    {
      register int *base_2_sele_b;
      base_0_sele_a = &base[0].sele[cNDummyAtoms];
      table_a = i_table + cNDummyAtoms;
      for(a=cNDummyAtoms;a<n_atom;a++) {
        if( (tag=*base_0_sele_a) ) {
          at1=&i_obj[table_a->model]->AtomInfo[table_a->atom];
          *base_0_sele_a = 0;
          table_b = i_table + cNDummyAtoms;
          base_2_sele_b = &base[2].sele[cNDummyAtoms];
          for(b=cNDummyAtoms;b<n_atom;b++) {
            if(*base_2_sele_b) {
              at2=&i_obj[table_b->model]->AtomInfo[table_b->atom];
              if(at1->resv==at2->resv)
                if((tolower(at1->chain[0]))==(tolower(at2->chain[0])))
                  if(WordMatchNoWild(G,at1->name,at2->name,ignore_case)<0)
                    if(WordMatchNoWild(G,at1->resi,at2->resi,ignore_case)<0)
                      if(WordMatchNoWild(G,at1->resn,at2->resn,ignore_case)<0)
                        if(WordMatchNoWild(G,at1->segi,at2->segi,ignore_case)<0) {
                          *base_0_sele_a = tag;
                          break;
                        }
            }
            base_2_sele_b++;
            table_b++;
          }
        }
        if( *(base_0_sele_a++) )
          c++;
        table_a++;
      }
    }
    break;
  case SELE_LIK2:
    {
      register int *base_2_sele_b;
      base_0_sele_a = &base[0].sele[cNDummyAtoms];
      table_a = i_table + cNDummyAtoms;
      for(a=cNDummyAtoms;a<n_atom;a++) {
        if( (tag=*base_0_sele_a) ) {
          at1=&i_obj[table_a->model]->AtomInfo[table_a->atom];
          *base_0_sele_a = 0;
          table_b = i_table + cNDummyAtoms;
          base_2_sele_b = &base[2].sele[cNDummyAtoms];
          for(b=cNDummyAtoms;b<n_atom;b++) {
            if(*base_2_sele_b) {
              at2=&i_obj[table_b->model]->AtomInfo[table_b->atom];
              if(at1->resv==at2->resv)
                if(WordMatchNoWild(G,at1->name,at2->name,ignore_case)<0)
                  if(WordMatchNoWild(G,at1->resi,at2->resi,ignore_case)<0) {
                    *base_0_sele_a = tag;
                    break;
                  }
            }
            base_2_sele_b++;
            table_b++;
          }
        }
        if( *(base_0_sele_a++) )
          c++;
        table_a++;
      }
    }
    break;
  }
  FreeP(base[2].sele);
  PRINTFD(G,FB_Selector)
    " SelectorLogic2: %d atoms selected.\n",c
    ENDFD;
  return(1);
}
/*========================================================================*/
int SelectorOperator22(PyMOLGlobals *G,EvalElem *base,int state)
{
  int c=0;
  int a,d,e;
  register CSelector *I=G->Selector;
  ObjectMolecule *obj;

  float dist;
  float *v2;
  CoordSet *cs;
  int ok=true;
  int nCSet;
  MapType *map;
  int i,j,h,k,l;
  int n1,at,idx;
  int code = base[1].code;

  if(state<0) {
    switch(state) {
    case -2:
    case -3:
      state=SceneGetState(G);
      break;
    }
  }

  switch(code) {
  case SELE_WIT_:
  case SELE_BEY_:
  case SELE_NTO_:
    if(!sscanf(base[2].text,"%f",&dist))
      ok=ErrMessage(G,"Selector","Invalid distance.");
    if(ok) {
      if(dist<0.0) dist = 0.0;
      
      /* copy starting mask */
      for(a=0;a<I->NAtom;a++) {
        I->Flag2[a]=base[0].sele[a];
        base[0].sele[a]=false;
      }
      
      for(d=0;d<I->NCSet;d++) {
        if((state<0)||(d==state)) {
          n1=0;
          for(a=0;a<I->NAtom;a++) {
            I->Flag1[a]=false;
            at=I->Table[a].atom;
            obj=I->Obj[I->Table[a].model];
            if(d<obj->NCSet) 
              cs=obj->CSet[d];
            else
              cs=NULL;
            if(cs) {
              if(obj->DiscreteFlag) {
                if(cs==obj->DiscreteCSet[at])
                  idx=obj->DiscreteAtmToIdx[at];
                else
                  idx=-1;
              } else 
                idx=cs->AtmToIdx[at];
              if(idx>=0) {
                copy3f(cs->Coord+(3*idx),I->Vertex+3*a);
                I->Flag1[a]=true;
                n1++;
              }
            }
          }
          if(n1) {
            map=MapNewFlagged(G,-dist,I->Vertex,I->NAtom,NULL,I->Flag1);
            if(map) {
              MapSetupExpress(map);
              nCSet=SelectorGetArrayNCSet(G,base[4].sele,false);
              for(e=0;e<nCSet;e++) {
                if((state<0)||(e==state)) {
                  for(a=0;a<I->NAtom;a++) {
                    if(base[4].sele[a]) {
                      at=I->Table[a].atom;
                      obj=I->Obj[I->Table[a].model];
                      if(e<obj->NCSet) 
                        cs=obj->CSet[e];
                      else
                        cs=NULL;
                      if(cs) {
                        if(obj->DiscreteFlag) {
                          if(cs==obj->DiscreteCSet[at])
                            idx=obj->DiscreteAtmToIdx[at];
                          else
                            idx=-1;
                        } else 
                          idx=cs->AtmToIdx[at];
                        if(idx>=0) {
                          v2 = cs->Coord+(3*idx);
                          MapLocus(map,v2,&h,&k,&l);
                          i=*(MapEStart(map,h,k,l));
                          if(i) {
                            j=map->EList[i++];
                            
                            while(j>=0) {
                              if(!base[0].sele[j])
                                if(I->Flag2[j])
                                  if(within3f(I->Vertex+3*j,v2,dist)) {
                                    if((code!=SELE_NTO_)||(!base[4].sele[j]))
                                      base[0].sele[j]=true;
                                  }
                              j=map->EList[i++];
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              MapFree(map);
            }
          }
        }
      }
      if(code==SELE_BEY_) {
        for(a=0;a<I->NAtom;a++) {
          if(I->Flag2[a])
            base[0].sele[a] = !base[0].sele[a];
        }
      }
      for(a=cNDummyAtoms;a<I->NAtom;a++)
        if(base[0].sele[a]) c++;
    }
    break;
  }
  FreeP(base[4].sele);
  PRINTFD(G,FB_Selector)
    " SelectorOperator22: %d atoms selected.\n",c
    ENDFD;
  return(1);
}

static void remove_quotes(char *st)
{ 
  /* nasty */

  SelectorWordType store;
  char *p,*q;
  char *quote_start = NULL;
  char active_quote = 0;
  p = st;
  q = store;
  /*  printf("DEBUG remove_quotes: input [%s]\n",st); */

  while(*p) {
    if(((*p)==34)||((*p)==39)) {
      if(quote_start&&(active_quote==*p)) { /* eliminate quotes... */
        while(quote_start<(q-1)) {
          *(quote_start)=*(quote_start+1);
          quote_start++;
        }
        q--;
        quote_start=NULL;
        p++;
      } else if(quote_start) {
        *(q++)=*(p++);
      } else {
        if(p==st) { /* at start => real quote */
          quote_start = q;
          active_quote = *p;
        } else if((*(p-1)=='+')||(*(p-1)==',')) { /* after separator => real quote */
          quote_start = q;
          active_quote = *p;
        }
        *(q++)=*(p++);
      }
    } else {
      /* UNWORKABLE -- hopelly getting rid of this kludge will not cause major grief 
      if((*p=='+')&&(!quote_start))
        if(!((*(p+1)==0)||(*(p+1)==',')||(*(p+1)=='+')))
          *p=',';
      */
      *(q++)=*(p++);
    }
  }
  *(q++) = 0;
  strcpy(st,store);

  /*  printf("DEBUG remove_quotes: output [%s]\n",st);*/
}

/*========================================================================*/
int *SelectorEvaluate(PyMOLGlobals *G,SelectorWordType *word,int state)
{
  int level = 0, imp_op_level = 0;
  int depth = 0;
  int a,b,c = 0;
  int ok=true;
  unsigned int code;
  int valueFlag = 0; /* are we expecting? */
  int *result = NULL;
  int opFlag,opFlag2,maxLevel;
  char *q,*cc1,*cc2;
  int totDepth=0;
  int exact;
  char *np;

  register int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
  EvalElem *Stack=NULL,*e;
  SelectorWordType tmpKW;
  Stack = VLAlloc(EvalElem,100);

  UtilZeroMem(Stack,sizeof(EvalElem)); /* blank first entry */

  /* converts all keywords into code, adds them into a operation list */
  while(ok&&word[c][0]) {
    if(word[c][0]=='#') {
      if((!valueFlag)&&(!level)) {
        word[c][0]=0; /* terminate selection if we encounter a comment */
        word[c+1][0]=0;
        break;
      }
    }
    switch(word[c][0]) {
      case 0:
        break;
      case '(':
        if(valueFlag)	ok=ErrMessage(G,"Selector","Misplaced (.");
        if(ok) level++;
        break;
      case ')':
        if(valueFlag)	ok=ErrMessage(G,"Selector","Misplaced ).");
        if(ok) {
          level--;
          if(level<0)
            ok=ErrMessage(G,"Selector","Syntax error.");
          else
            imp_op_level = level;
        }
        if(ok&&depth) Stack[depth].level--;
        break;
      default:
        if(valueFlag>0) { /* standard operand */
          depth++;
          VLACheck(Stack,EvalElem,depth);
          e=Stack+depth;
          e->level=(level<<4)+1; /* each nested paren is 16 levels of priority */
          e->imp_op_level=(imp_op_level<<4)+1; 
          imp_op_level = level;
          e->type=STYP_VALU;
          cc1 = word[c];
          cc2 = e->text;
          strcpy(e->text,cc1);
          remove_quotes(e->text);
          valueFlag--;
        } else if(valueFlag<0) { /* operation parameter i.e. around X<-- */
          depth++;
          VLACheck(Stack,EvalElem,depth);
          e=Stack+depth;
          e->level=(level<<4)+1;
          e->imp_op_level=(imp_op_level<<4)+1; 
          imp_op_level = level;
          e->type=STYP_PVAL;
          strcpy(e->text,word[c]);
          valueFlag++;
        } else { /* possible keyword... */
          if((word[c][0]=='*')&&(!word[c][1]))
            code = SELE_ALLz;
          else
            code = WordKey(G,Keyword,word[c],4,ignore_case,&exact);
          if(!code) {
            b=strlen(word[c])-1;
            if((b>2)&&(word[c][b]==';')) {
              /* kludge to accomodate unnec. ';' usage */
              word[c][b]=0;
              code=WordKey(G,Keyword,word[c],4,ignore_case,&exact);
            }
          }
          PRINTFD(G,FB_Selector)
            " Selector: code %x\n",code
            ENDFD;
          if((code>0)&&(!exact))  
            if(SelectorIndexByName(G,word[c])>=0)
              code=0; /* favor selections over partial keyword matches */
          if(code) {
            /* this is a known operation */
            depth++;
            VLACheck(Stack,EvalElem,depth);
            e=Stack+depth;
            e->code=code;
            e->level=(level<<4)+((e->code&0xF0)>>4);
            e->imp_op_level=(imp_op_level<<4)+1; 
            imp_op_level = level;
            e->type=(e->code&0xF);
            switch(e->type)
              {
              case STYP_SEL0:
                valueFlag=0;
                break;
              case STYP_SEL1:
                valueFlag=1;
                break;
              case STYP_SEL2:
                valueFlag=2;
                break;
              case STYP_OPR1:
                valueFlag=0;
                break;
              case STYP_OPR2:
                valueFlag=0;
                break;
              case STYP_PRP1: 
                valueFlag=-1;
                break;
              case STYP_OP22: 
                valueFlag=-2;
                break;
              }
          } else {
            strcpy(tmpKW,word[c]); /* handle <object-name>[`index] syntax */
            if((np=strstr(tmpKW,"`"))) { /* must be an object name */
              *np=0;
              depth++;
              VLACheck(Stack,EvalElem,depth);
              e=Stack+depth;
              e->code=SELE_MODs;
              e->level=(level<<4)+((e->code&0xF0)>>4);
              e->imp_op_level=(imp_op_level<<4)+1; 
              imp_op_level = level;
              e->type=STYP_SEL1;
              valueFlag = 1;
              c--;
            } else { /* handle <selection-name> syntax */
              depth++;
              VLACheck(Stack,EvalElem,depth);
              e=Stack+depth;
              e->code=SELE_SELs;
              e->level=(level<<4)+((e->code&0xF0)>>4);
              e->imp_op_level=(imp_op_level<<4)+1; 
              imp_op_level = level;
              e->type=STYP_SEL1;
              valueFlag = 1;
              c--;
            }
          }
        }
        break;
    }
    if(ok) c++; /* go onto next word */
  }
  if(level>0)
    ok=ErrMessage(G,"Selector","Malformed selection.");		
  if(ok) {/* this is the main operation loop */
    totDepth=depth;
    opFlag=true;
    maxLevel=-1;
    for(a=1;a<=totDepth;a++) {
      PRINTFD(G,FB_Selector)
        " Selector initial stack %d-%p lv: %x co: %d type: %x sele %p\n",
        a,(void*)(Stack+a),Stack[a].level,Stack[a].code,
        Stack[a].type,(void*)Stack[a].sele
        ENDFD;
      
      if(Stack[a].level>maxLevel) 
        maxLevel=Stack[a].level;
    }
    level=maxLevel;
    PRINTFD(G,FB_Selector)
      " Selector: maxLevel %d %d\n",maxLevel,totDepth
      ENDFD;
    if(level>=0) 
      while(ok) { /* loop until all ops at all levels have been tried */

        /* order & efficiency of this algorithm could be improved... */

        PRINTFD(G,FB_Selector)
          " Selector: new cycle...\n"
          ENDFD;
        depth = 1;
        opFlag=true;
        while(ok&&opFlag) { /* loop through all entries looking for ops at the current level */
          PRINTFD(G,FB_Selector)
            " Selector: lvl: %d de:%d-%p slv:%d co: %x typ %x sele %p td: %d\n",
            level,depth,(void*)(Stack+depth),Stack[depth].level,
            Stack[depth].code,
            Stack[depth].type,(void*)Stack[depth].sele,totDepth
            ENDFD;
          
          opFlag=false;
            
          if(Stack[depth].level>=level) {
            Stack[depth].level=level; /* trim peaks */
          }
          if(ok) 
            if(depth>0)
              if((!opFlag)&&(Stack[depth].type==STYP_SEL0)) {
                opFlag=true;
                ok=SelectorSelect0(G,&Stack[depth]);
              }
          if(ok)
            if(depth>1)
              if(Stack[depth-1].level>=Stack[depth].level) {
                  if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_SEL1)
                     &&(Stack[depth].type==STYP_VALU)) {
                    /* 1 argument selection operator */
                      opFlag=true;
                      ok=SelectorSelect1(G,&Stack[depth-1]);
                      for(a=depth+1;a<=totDepth;a++) 
                        Stack[a-1]=Stack[a];
                      totDepth--;
                    } else if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_OPR1)
                          &&(Stack[depth].type==STYP_LIST)) {
                     /* 1 argument logical operator */
                      opFlag=true;
                      ok=SelectorLogic1(G,&Stack[depth-1],state);
                      for(a=depth+1;a<=totDepth;a++) 
                        Stack[a-1]=Stack[a];
                      totDepth--;
                    } else if((Stack[depth-1].type==STYP_LIST)&&
                          (Stack[depth].type==STYP_LIST)&&
                          (!((Stack[depth-1].level & 0xF) ||
                             (Stack[depth].level & 0xF)))) {
                    /* two adjacent lists at zeroth priority level
                       for the scope (lowest nibble of level is
                       zero) is an implicit OR action */
                    VLACheck(Stack,EvalElem,totDepth);
                    for(a=totDepth;a>=depth;a--)
                      Stack[a+1]=Stack[a];
                    totDepth++;
                    Stack[depth].type = STYP_OPR2;
                    Stack[depth].code = SELE_IOR2;
                    Stack[depth].level = Stack[depth].imp_op_level;
                    Stack[depth].sele = NULL;
                    Stack[depth].text[0] = 0;
                    if(level<Stack[depth].level)
                      level = Stack[depth].level;
                    opFlag=true;
                  }
                }
          if(ok)
            if(depth>2)
              if((Stack[depth-1].level>=Stack[depth].level)&&
                 (Stack[depth-1].level>=Stack[depth-2].level)) {

                  if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_OPR2)
                     &&(Stack[depth].type==STYP_LIST)
                     &&(Stack[depth-2].type==STYP_LIST)) {
                    /* 2 argument logical operator */
                      ok=SelectorLogic2(G,&Stack[depth-2]);
                      opFlag=true;
                      for(a=depth+1;a<=totDepth;a++) 
                        Stack[a-2]=Stack[a];
                      totDepth-=2;
                  } else if(ok&&(!opFlag)&&(Stack[depth-1].type==STYP_PRP1)
                            &&(Stack[depth].type==STYP_PVAL)
                            &&(Stack[depth-2].type==STYP_LIST)) {
                    /* 2 argument logical operator */
                    ok=SelectorModulate1(G,&Stack[depth-2],state);
                    opFlag=true;
                    for(a=depth+1;a<=totDepth;a++) 
                      Stack[a-2]=Stack[a];
                    totDepth-=2;
                  }
              }
          if(ok)
            if(depth>2)
              if((Stack[depth-2].level>=Stack[depth-1].level)&&
                 (Stack[depth-2].level>=Stack[depth].level)) {

                if(ok&&(!opFlag)&&(Stack[depth-2].type==STYP_SEL2)
                   &&(Stack[depth-1].type==STYP_VALU)
                   &&(Stack[depth].type==STYP_VALU)) {
                  /* 2 argument value operator */
                  ok=SelectorSelect2(G,&Stack[depth-2]);
                  opFlag=true;
                  for(a=depth+1;a<=totDepth;a++) 
                    Stack[a-2]=Stack[a];
                  totDepth-=2;
                }
              }
          if(ok)
            if(depth>3)
              if((Stack[depth-3].level>=Stack[depth].level)&&
                 (Stack[depth-3].level>=Stack[depth-1].level)&&
                 (Stack[depth-3].level>=Stack[depth-2].level)) {
                
                if(ok&&(!opFlag)&&(Stack[depth-3].type==STYP_SEL3)
                   &&(Stack[depth].type==STYP_VALU)
                   &&(Stack[depth-1].type==STYP_VALU)
                   &&(Stack[depth-2].type==STYP_VALU)) {
                  /* 2 argument logical operator */
                  /*								ok=SelectorSelect3(G,&Stack[depth-3]);*/
                  opFlag=true;
                  for(a=depth+1;a<=totDepth;a++) 
                    Stack[a-3]=Stack[a];
                  totDepth-=3;
                }
              }
          if(ok)
            if(depth>4)
              if((Stack[depth-3].level>=Stack[depth].level)&&
                 (Stack[depth-3].level>=Stack[depth-1].level)&&
                 (Stack[depth-3].level>=Stack[depth-2].level)&&
                 (Stack[depth-3].level>=Stack[depth-4].level)) {
                
                if(ok&&(!opFlag)&&(Stack[depth-3].type==STYP_OP22)
                   &&(Stack[depth-1].type==STYP_VALU)
                   &&(Stack[depth-2].type==STYP_VALU)
                   &&(Stack[depth].type==STYP_LIST)
                   &&(Stack[depth-4].type==STYP_LIST)) {
                     
                  ok=SelectorOperator22(G,&Stack[depth-4],state);
                  opFlag=true;
                  for(a=depth+1;a<=totDepth;a++) 
                    Stack[a-4]=Stack[a];
                  totDepth-=4;
                }

              }
          if(opFlag) {
            opFlag2=true; /* make note that we performed an operation */
            depth = 1; /* start back at the left hand side */
          } else {
            depth = depth + 1;
            opFlag=true;
            if(depth>totDepth)
              break;
          }
        }
        if(level)
          level--;
        else
          break;
      }
    depth=totDepth;
  }
  if(ok) {
    if(depth!=1) {
      ok=ErrMessage(G,"Selector","Malformed selection.");
    } else if(Stack[depth].type!=STYP_LIST)
      ok=ErrMessage(G,"Selector","Invalid selection.");
    else
      result=Stack[totDepth].sele; /* return the selection list */
  }
  if(!ok) {
    for(a=1;a<=depth;a++) {
      PRINTFD(G,FB_Selector)
        " Selector: releasing %d %x %p\n",a,Stack[a].type,(void*)Stack[a].sele
        ENDFD;
      if(Stack[a].type==STYP_LIST)
        FreeP(Stack[a].sele);
    }
    depth=0;
    {
      OrthoLineType line;
      for(a=0;a<=c;a++) {
        q=line;
        if(a&&word[a][0])
          q=UtilConcat(q," ");
        q=UtilConcat(q,word[a]);
        PRINTFB(G,FB_Selector,FB_Errors) 
          "%s",line
          ENDFB(G);
      }
      q=line;
      q=UtilConcat(q,"<--");
      PRINTFB(G,FB_Selector,FB_Errors) 
        "%s",line
        ENDFB(G);
      OrthoRestorePrompt(G);
    }
  }
  VLAFreeP(Stack);
  if(!ok) {
    FreeP(result);
    result = NULL;
  }
  return(result);
}

/*========================================================================*/
SelectorWordType *SelectorParse(PyMOLGlobals *G,char *s) {

  /* break a selection down into its constituent strings and
	  return them in a SelectorWordType VLA, null string terminated */

  SelectorWordType *r = NULL;
  int c=0;
  int w_flag=false;
  int quote_flag=false;
  char quote_char = '"';
  char *p = s;
  char *q = NULL, *q_base = NULL;
  r=VLAlloc(SelectorWordType,100);
  while(*p) 
	 {
		if(w_flag) /* currently in a word, thus q is a valid pointer */
		  {
          if(quote_flag) {
            if(*p!=quote_char) {
              *q++=*p;
            } else {
              quote_flag=false;
              *q++=*p;
            }
          } else switch(*p)
				{
				case ' ':
				  *q=0;
				  w_flag=false;
				  break;
            case ';': /* special word terminator */
				  *q++=*p;
				  *q=0;
              w_flag=false;
              break;
				case '!': /* single words */ 
				case '&': 
				case '|': 
				case '(': 
				case ')':
				case '>': 
				case '<':
            case '=':
				case '%':
				  *q=0; /* terminate current word */
				  c++;
				  VLACheck(r,SelectorWordType,c); /* add new word */
				  q=r[c-1];
				  *q++=*p;
				  *q=0;  /* terminate current word */
				  w_flag=false;
				  break;
            case '"':
              quote_flag = true;
              *q++=*p;
              break;
				default:
				  *q++=*p;
				  break;
				}
          if(w_flag) {
            if((q-q_base)>=sizeof(SelectorWordType)) {
              q_base[sizeof(SelectorWordType)-1]=0;
              w_flag=false;
              PRINTFB(G,FB_Selector,FB_Errors) 
                "Selector-Error: Word too long. Truncated:\nSelector-Error: %s...\n",q_base
                ENDFB(G);
            }
          }
		  }
		else /*outside a word -- q is undefined */
		  {
			 switch(*p)
				{
#if 0
				case '*': /* special case */
				  c++;
				  VLACheck(r,SelectorWordType,c);
				  q=r[c-1];
				  *q++='+';
				  *q=0;
				  w_flag=false;
              break;
#endif
				case '!': /* single words */ 
				case '&': 
				case '|': 
				case '(': 
				case ')':
				case '>': 
				case '<':
            case '=':
				case '%':
				  c++;
				  VLACheck(r,SelectorWordType,c);
				  q=r[c-1];
				  *q++=(*p);
				  *q=0;
				  break;
            case ' ':
              break;
            case '"':
              quote_flag = true;
              quote_char = *p;
              w_flag=true;
				  c++;
              VLACheck(r,SelectorWordType,c);
				  q=r[c-1];
              q_base = q;
				  *q++=*p;
              break;
				default:
				  w_flag=true;
				  c++;
              VLACheck(r,SelectorWordType,c);
				  q=r[c-1];
              q_base = q;
				  *q++=*p;
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
  if(Feedback(G,FB_Selector,FB_Debugging)) 
	 {
		c=0;
		while(r[c][0])
		  {
			 fprintf(stderr,"word: %s\n",r[c]);
			 c++;
		  }
	 }
  return(r);
}
/*========================================================================*/
void SelectorFree(PyMOLGlobals *G)
{
  register CSelector *I = G->Selector;
  SelectorClean(G);
  if(I->Origin)
    if(I->Origin->Obj.fFree)
      I->Origin->Obj.fFree((CObject*)I->Origin);
  if(I->Center)
    if(I->Center->Obj.fFree)
      I->Center->Obj.fFree((CObject*)I->Center);
  VLAFreeP(I->Member);
  VLAFreeP(I->Name);
  VLAFreeP(I->Info);
  OVLexicon_DEL_AUTO_NULL(I->Lex);
  OVOneToAny_DEL_AUTO_NULL(I->Key);
  OVOneToOne_DEL_AUTO_NULL(I->NameOffset);

  FreeP(G->Selector);
}


/*========================================================================*/

void SelectorMemoryDump(PyMOLGlobals *G)
{
  register CSelector *I = G->Selector;
  printf(" SelectorMemory: NSelection %d\n",I->NSelection);
  printf(" SelectorMemory: NActive %d\n",I->NActive);
  printf(" SelectorMemory: TmpCounter %d\n",I->TmpCounter);
  printf(" SelectorMemory: NMember %d\n",I->NMember);  
}


static void SelectorInit2(PyMOLGlobals *G)
{
  register CSelector *I = G->Selector;

  I->NSelection = 0;
  I->NActive=0;
  I->TmpCounter = 0;

  I->NMember=0;
  I->FreeMember=0;
  I->NCSet=0;
  
  I->Lex = OVLexicon_New(G->Context->heap);
  I->Key = OVOneToAny_New(G->Context->heap);
  I->NameOffset = OVOneToOne_New(G->Context->heap);

  {  /* create placeholder "all" selection, which is selection 0
      and "none" selection, which is selection 1 */
    int n;

    n=I->NActive;
    VLACheck(I->Name,SelectorWordType,n+1);
    VLACheck(I->Info,SelectionInfoRec,n+1);
    strcpy(I->Name[n],cKeywordAll); /* "all" selection = 0 */
    I->Name[n+1][0]=0;
    SelectorAddName(G,n);
    SelectionInfoInit(I->Info + n);
    I->Info[n].ID = I->NSelection++;
    I->NActive++;

    n=I->NActive;
    VLACheck(I->Name,SelectorWordType,n+1);
    VLACheck(I->Info,SelectionInfoRec,n+1);
    strcpy(I->Name[n],cKeywordNone); /* "none" selection = 1*/
    I->Name[n+1][0]=0;
    SelectorAddName(G,n);
    SelectionInfoInit(I->Info + n);
    I->Info[n].ID = I->NSelection++;
    I->NActive++;
  }


  if(I->Lex && I->Key) {
    int a=0;
    OVreturn_word result;
    while(1) {
      if(!Keyword[a].word[0]) break;
      if(OVreturn_IS_OK( (result = OVLexicon_GetFromCString(I->Lex,Keyword[a].word)))) {
        OVOneToAny_SetKey(I->Key, result.word, Keyword[a].value);
      }
      a++;
    }

  }
}

void SelectorReinit(PyMOLGlobals *G)
{
  register CSelector *I=G->Selector;
  SelectorClean(G);

  OVLexicon_DEL_AUTO_NULL(I->Lex);
  OVOneToAny_DEL_AUTO_NULL(I->Key);
  OVOneToOne_DEL_AUTO_NULL(I->NameOffset);

  SelectorInit2(G);
}

/*========================================================================*/
int SelectorInit(PyMOLGlobals *G)
{
 register CSelector *I=NULL;
  if( (I=(G->Selector=Calloc(CSelector,1)))) {
    
    I->Name = VLAlloc(SelectorWordType,10);
    I->Info = VLAlloc(SelectionInfoRec,10);
    

    I->Member = (MemberType*)VLAMalloc(100,sizeof(MemberType),5,true);
    I->Vertex=NULL;
    I->Origin=NULL;
    I->Table=NULL;
    I->Obj=NULL;
    I->Flag1=NULL;
    I->Flag2=NULL;
    
    SelectorInit2(G);
    return 1;
  } else
    return 0;
}
/*========================================================================*/


DistSet *SelectorGetDistSet(PyMOLGlobals *G,DistSet *ds,
                            int sele1,int state1,int sele2,int state2,
                            int mode,float cutoff,float *result)
{
  register CSelector *I=G->Selector;
  int *vla=NULL;
  int c;
  float dist;
  int a1,a2;
  AtomInfoType *ai1,*ai2;
  int at,at1,at2;
  CoordSet *cs1,*cs2;
  ObjectMolecule *obj,*obj1,*obj2,*lastObj;
  int idx1,idx2;
  int a;
  int nv = 0;
  float *vv=NULL,*vv0,*vv1;
  float dist_sum=0.0;
  int dist_cnt = 0;
  int s;
  int a_keeper = false;
  int *zero=NULL,*scratch=NULL,*coverage=NULL;
  HBondCriteria hbcRec,*hbc;
  int exclusion = 0;
  int bonds_only = 0;
  int from_proton = SettingGetGlobal_b(G,cSetting_h_bond_from_proton);

  switch(mode) {
  case 1:
    bonds_only = 1;
    break;
  case 2:
    exclusion = SettingGetGlobal_i(G,cSetting_h_bond_exclusion);
    break;
  case 3:
    exclusion = SettingGetGlobal_i(G,cSetting_distance_exclusion);
    break;
  }

  hbc=&hbcRec;
  *result = 0.0;
  if(!ds) {
    ds = DistSetNew(G);
  } else {
    vv = ds->Coord;
    nv = ds->NIndex;
  }
  if(!vv) {
    vv = VLAlloc(float,10);
  }

  if((state1<0)||(state2<0)||(state1!=state2)) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,state1,-1);
  }

  coverage=Calloc(int,I->NAtom);

  for(a=cNDummyAtoms;a<I->NAtom;a++) {
    at=I->Table[a].atom;
    obj=I->Obj[I->Table[a].model];
    s=obj->AtomInfo[at].selEntry;
    if(SelectorIsMember(G,s,sele1))
      coverage[a]++;
    if(SelectorIsMember(G,s,sele2))
      coverage[a]++;
  }

  if((mode==1)||(mode==2)||(mode==3)) { /* fill in all the neighbor tables */
    int max_n_atom = I->NAtom;
    lastObj=NULL;
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      obj=I->Obj[I->Table[a].model];
      s=obj->AtomInfo[at].selEntry;
      if(obj!=lastObj) {
        if(max_n_atom<obj->NAtom)
          max_n_atom = obj->NAtom;
        if(SelectorIsMember(G,s,sele1)||
           SelectorIsMember(G,s,sele2)) {
          ObjectMoleculeUpdateNeighbors(obj);
          if(mode==2)
            ObjectMoleculeVerifyChemistry(obj,-1);
          lastObj = obj;
        }
      }
    }
    zero=Calloc(int,max_n_atom);
    scratch=Alloc(int,max_n_atom);
  }

  if(mode==2) {
    ObjectMoleculeInitHBondCriteria(G,hbc);
    if(cutoff<0.0F) {
      cutoff = hbc->maxDistAtMaxAngle;
      if(cutoff<hbc->maxDistAtZero) {
        cutoff = hbc->maxDistAtZero; 
      }
    }
  }
  if(cutoff<0) cutoff = 1000.0;
  c=SelectorGetInterstateVLA(G,sele1,state1,sele2,state2,cutoff,&vla);
  for(a=0;a<c;a++) {
    a1=vla[a*2];
    a2=vla[a*2+1];

    if((a1!=a2)&&(
                  (!((coverage[a1]==2)&&(coverage[a2]==2)))|| 
                  (a1<a2))) /* eliminate reverse duplicates */
      {
      at1=I->Table[a1].atom;
      at2=I->Table[a2].atom;
      
      obj1=I->Obj[I->Table[a1].model];
      obj2=I->Obj[I->Table[a2].model];
      
      if((state1<obj1->NCSet)&&(state2<obj2->NCSet)) {
        cs1=obj1->CSet[state1];
        cs2=obj2->CSet[state2];
        if(cs1&&cs2) { 
    
          float *don_vv = NULL;
          float *acc_vv = NULL;

          ai1=obj1->AtomInfo+at1;
          ai2=obj2->AtomInfo+at2;

          if(obj1->DiscreteFlag) {
            if(cs1==obj1->DiscreteCSet[at1]) {
              idx1=obj1->DiscreteAtmToIdx[at1];
            } else {
              idx1=-1;
            }
          } else {
            idx1=cs1->AtmToIdx[at1];
          }
          
          if(obj2->DiscreteFlag) {
            if(cs2==obj2->DiscreteCSet[at2]) {
              idx2=obj2->DiscreteAtmToIdx[at2];
            } else {
              idx2=-1;
            }
            
          } else {
            idx2=cs2->AtmToIdx[at2];
          }
          
          if((idx1>=0)&&(idx2>=0)) {
            dist=(float)diff3f(cs1->Coord+3*idx1,cs2->Coord+3*idx2);
            
            if(dist<cutoff) {
              
              float h_crd[3];
              int h_real = false;

              a_keeper=true;
              if(exclusion && (obj1==obj2)) {
                a_keeper = !SelectorCheckNeighbors(G,exclusion,
                                                   obj1,at1,at2,
                                                   zero,scratch);
              } else if(bonds_only) {
                a_keeper = SelectorCheckNeighbors(G,1,
                                                  obj1,at1,at2,
                                                  zero,scratch);
              }
              if(a_keeper&&(mode==2)) {
                if(ai1->hb_donor&&ai2->hb_acceptor) {
                  a_keeper = ObjectMoleculeGetCheckHBond(&h_real,
                                                         h_crd,
                                                         obj1,
                                                         at1,
                                                         state1,
                                                         obj2,
                                                         at2,
                                                         state2,
                                                         hbc);
                  if(a_keeper) {
                    if(h_real && from_proton) 
                      don_vv = h_crd;
                    else
                      don_vv = cs1->Coord + 3*idx1;
                    acc_vv = cs2->Coord + 3*idx2;
                  }
                } else if(ai1->hb_acceptor&&ai2->hb_donor) {
                  a_keeper = ObjectMoleculeGetCheckHBond(&h_real,
                                                         h_crd,
                                                         obj2,
                                                         at2,
                                                         state2,
                                                         obj1,
                                                         at1,
                                                         state1,
                                                         hbc);
                  
                  if(a_keeper) {
                    if(h_real && from_proton) 
                      don_vv = h_crd;
                    else
                      don_vv = cs2->Coord + 3*idx2;
                    acc_vv = cs1->Coord + 3*idx1;
                  }
                } else {
                  a_keeper = false;
                }
              }
              if((sele1==sele2)&&(at1>at2))
                a_keeper = false;

              if(a_keeper) {
                
                dist_cnt++;
                dist_sum+=dist;
                VLACheck(vv,float,(nv*3)+6);
                vv0 = vv + (nv*3);

                if((mode==2)&&(don_vv)&&(acc_vv)) {
                  *(vv0++) = *(don_vv++);
                  *(vv0++) = *(don_vv++);
                  *(vv0++) = *(don_vv++);
                  *(vv0++) = *(acc_vv++);
                  *(vv0++) = *(acc_vv++);
                  *(vv0++) = *(acc_vv++);
                } else {
                  vv1 = cs1->Coord+3*idx1;
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  vv1 = cs2->Coord+3*idx2;
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                }
                
                nv+=2;
              }
            }
          }
        }
      }
    }
  }
  if(dist_cnt)
    (*result)=dist_sum/dist_cnt;
  VLAFreeP(vla);
  FreeP(zero);
  FreeP(scratch);
  FreeP(coverage);
  if(vv)
    VLASize(vv,float,(nv+1)*3);
  ds->NIndex = nv;
  ds->Coord = vv;
  return(ds);
}


DistSet *SelectorGetAngleSet(PyMOLGlobals *G, DistSet *ds,
                             int sele1,int state1,
                             int sele2,int state2,
                             int sele3,int state3,
                             int mode, float *angle_sum,
                             int *angle_cnt)
{
  register CSelector *I=G->Selector;
  float *vv = NULL;
  int nv = 0;
  int *coverage=NULL;

  if(!ds) {
    ds = DistSetNew(G);
  } else {
    vv = ds->AngleCoord;
    nv = ds->NAngleIndex;
  }
  if(!vv)
    vv = VLAlloc(float,10);

  if((state1<0)||(state2<0)||(state3<0)||(state1!=state2)||(state1!=state3)) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,state1,-1);
  }


  /* which atoms are involved? */

  {
    int a, s, at;
    ObjectMolecule *obj;

    coverage=Calloc(int,I->NAtom);
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      obj=I->Obj[I->Table[a].model];
      s=obj->AtomInfo[at].selEntry;
      if(SelectorIsMember(G,s,sele1))
        coverage[a]++;
      if(SelectorIsMember(G,s,sele2))
        coverage[a]++;
      if(SelectorIsMember(G,s,sele3))
        coverage[a]++;
    }
  }

  { /* fill in neighbor tables */
    int a, s, at;
    ObjectMolecule *obj,*lastObj = NULL;
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      obj=I->Obj[I->Table[a].model];
      s=obj->AtomInfo[at].selEntry;
      if(obj!=lastObj) {
        if(SelectorIsMember(G,s,sele1)||
           SelectorIsMember(G,s,sele2)||
           SelectorIsMember(G,s,sele3)) {
          ObjectMoleculeUpdateNeighbors(obj);
          lastObj = obj;
        }
      }
    }
  }

  {
    int a, s, at;
    ObjectMolecule *obj;
    int *list1 = VLAlloc(int,1000);
    int *list2 = VLAlloc(int,1000);
    int *list3 = VLAlloc(int,1000);
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int bonded12, bonded23;

    /* now generate three lists of atoms, one for each selection set */

    if(list1&&list2&&list3) {
      for(a=cNDummyAtoms;a<I->NAtom;a++) {
        at=I->Table[a].atom;
        obj=I->Obj[I->Table[a].model];
        s=obj->AtomInfo[at].selEntry;
        if(SelectorIsMember(G,s,sele1)) {
          VLACheck(list1,int,n1);
          list1[n1++] = a;
        }
        if(SelectorIsMember(G,s,sele2)) {
          VLACheck(list2,int,n2);
          list2[n2++] = a;
        }
        if(SelectorIsMember(G,s,sele3)) {
          VLACheck(list3,int,n3);
          list3[n3++] = a;
        }
      }
      
      /* for each set of 3 atoms in each selection... */
      
      {
        int i1,i2,i3;
        int a1,a2,a3;
        int at1,at2,at3;
        
        /*        AtomInfoType *ai1,*ai2,ai3;*/
        CoordSet *cs1,*cs2,*cs3;
        ObjectMolecule *obj1,*obj2,*obj3;
        
        int idx1,idx2,idx3;
        float angle;
        float d1[3],d2[3];
        float *v1,*v2,*v3, *vv0;

        for(i1=0;i1<n1;i1++) {
          a1 = list1[i1];
          at1=I->Table[a1].atom;
          obj1=I->Obj[I->Table[a1].model];

          if(state1<obj1->NCSet) {
            cs1=obj1->CSet[state1];

            if(cs1) {
              if(obj1->DiscreteFlag) {
                if(cs1==obj1->DiscreteCSet[at1]) {
                  idx1=obj1->DiscreteAtmToIdx[at1];
                } else {
                  idx1=-1;
                }
              } else {
                idx1=cs1->AtmToIdx[at1];
              }
              
              if(idx1>=0) {
                
                for(i2=0;i2<n2;i2++) {
                  a2 = list2[i2];
                  at2=I->Table[a2].atom;
                  obj2=I->Obj[I->Table[a2].model];
                  
                  if(state2<obj2->NCSet) {
                    
                    cs2=obj2->CSet[state2];
                    
                    if(cs2) {
                      if(obj2->DiscreteFlag) {
                        if(cs2==obj2->DiscreteCSet[at2]) {
                          idx2=obj2->DiscreteAtmToIdx[at2];
                        } else {
                          idx2=-1;
                        }
                      } else {
                        idx2=cs2->AtmToIdx[at2];
                      }
                    
                      if(idx2>=0) {
                        
                        bonded12 = ObjectMoleculeAreAtomsBonded2(obj1,at1,obj2,at2); 

                        for(i3=0;i3<n3;i3++) {
                          a3 = list3[i3];
                          
                          if((a1!=a2) && (a2!=a3) && (a1!=a3)) {
                            if((!((coverage[a1]==3)&&(coverage[a2]==3)&&(coverage[a3]==3))) 
                               ||(a1<a3)) { /* eliminate alternate-order duplicates */
                              
                              at3=I->Table[a3].atom;
                              obj3=I->Obj[I->Table[a3].model];
                              
                              if(state3<obj3->NCSet) {
                                
                                cs3=obj3->CSet[state3];
                                
                                if(cs3) { 
                                  if(obj3->DiscreteFlag) {
                                    if(cs3==obj3->DiscreteCSet[at3]) {
                                      idx3=obj3->DiscreteAtmToIdx[at3];
                                    } else {
                                      idx3=-1;
                                    }
                                  } else {
                                    idx3=cs3->AtmToIdx[at3];
                                  }
                                  
                                  if(idx3>=0) {
                                    
                                    bonded23 = ObjectMoleculeAreAtomsBonded2(obj2,at2,obj3,at3); 
                                    
                                    if(!mode || ((mode==1)&&(bonded12&&bonded23))) { /* store the 3 coordinates */
                                      
                                      v1 = cs1->Coord+3*idx1;
                                      v2 = cs2->Coord+3*idx2;
                                      v3 = cs3->Coord+3*idx3;
                                      
                                      subtract3f(v1,v2,d1);
                                      subtract3f(v3,v2,d2);
                                      
                                      angle = get_angle3f(d1,d2);
                                      
                                      (*angle_sum)+=angle;
                                      (*angle_cnt)++;
                                      
                                      VLACheck(vv,float,(nv*3)+14);
                                      vv0 = vv+ (nv*3);
                                      *(vv0++) = *(v1++);
                                      *(vv0++) = *(v1++);
                                      *(vv0++) = *(v1++);
                                      *(vv0++) = *(v2++);
                                      *(vv0++) = *(v2++);
                                      *(vv0++) = *(v2++);
                                      *(vv0++) = *(v3++);
                                      *(vv0++) = *(v3++);
                                      *(vv0++) = *(v3++);
                                      *(vv0++) = (float)!bonded12;
                                      /* show line 1 flag*/
                                      *(vv0++) = (float)!bonded23;
                                      *(vv0++) = 0.0F; 
                                      *(vv0++) = 0.0F; /* label x relative to v2 */
                                      *(vv0++) = 0.0F; /* label y relative to v2 */
                                      *(vv0++) = 0.0F; /* label z relative to v2 */
                                      nv+=5;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    VLAFreeP(list1);
    VLAFreeP(list2);
    VLAFreeP(list3);
  }

  FreeP(coverage);
  if(vv)
    VLASize(vv,float,(nv+1)*3);
  ds->NAngleIndex = nv;
  ds->AngleCoord = vv;
  return(ds);
}

DistSet *SelectorGetDihedralSet(PyMOLGlobals *G, DistSet *ds,
                                int sele1,int state1,
                                int sele2,int state2,
                                int sele3,int state3,
                                int sele4,int state4,
                                int mode, float *angle_sum,
                                int *angle_cnt)
{
  register CSelector *I=G->Selector;
  float *vv = NULL;
  int nv = 0;
  int *coverage=NULL;
  ObjectMolecule *just_one_object = NULL;
  int just_one_atom[4] = {-1, -1, -1, -1};

  if(!ds) {
    ds = DistSetNew(G);
  } else {
    vv = ds->DihedralCoord;
    nv = ds->NDihedralIndex;
  }
  if(!vv)
    vv = VLAlloc(float,10);

  if((state1<0)||(state2<0)||(state3<0)||(state4<0)||
     (state1!=state2)||(state1!=state3)||(state1!=state4)) {
    SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
  } else {
    SelectorUpdateTable(G,state1,-1);
  }

  /* which atoms are involved? */

  {
    int a, s, at;
    ObjectMolecule *obj;

    coverage=Calloc(int,I->NAtom);
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      obj=I->Obj[I->Table[a].model];
      if(!a) just_one_object = obj;
      s=obj->AtomInfo[at].selEntry;
      if(SelectorIsMember(G,s,sele1)) {
        if(obj!=just_one_object)
          just_one_object = NULL;
        else if(just_one_atom[0]==-1)
          just_one_atom[0] = a;
        else
          just_one_atom[0] = -2;
        coverage[a]++;
      }
      if(SelectorIsMember(G,s,sele2)) {
        if(obj!=just_one_object)
          just_one_object = NULL;
        else if(just_one_atom[1]==-1)
          just_one_atom[1] = a;
        else
          just_one_atom[1] = -2;
        coverage[a]++;
      }
      if(SelectorIsMember(G,s,sele3)) {
        if(obj!=just_one_object)
          just_one_object = NULL;
        else if(just_one_atom[2]==-1)
          just_one_atom[2] = a;
        else
          just_one_atom[2] = -2;
        coverage[a]++;
      }
      if(SelectorIsMember(G,s,sele4)) {
        if(obj!=just_one_object)
          just_one_object = NULL;
        else if(just_one_atom[3]==-1)
          just_one_atom[3] = a;
        else
          just_one_atom[3] = -2;
        coverage[a]++;
      }
    }
  }

  if(just_one_object) {
    ObjectMoleculeUpdateNeighbors(just_one_object);
  } else { /* fill in neighbor tables */
    int a, s, at;
    ObjectMolecule *obj,*lastObj = NULL;
    for(a=cNDummyAtoms;a<I->NAtom;a++) {
      at=I->Table[a].atom;
      obj=I->Obj[I->Table[a].model];
      s=obj->AtomInfo[at].selEntry;
      if(obj!=lastObj) {
        if(SelectorIsMember(G,s,sele1)||
           SelectorIsMember(G,s,sele2)||
           SelectorIsMember(G,s,sele3)||
           SelectorIsMember(G,s,sele4)
           ) {
          ObjectMoleculeUpdateNeighbors(obj);
          lastObj = obj;
        }
      }
    }
  }


  {
    int a, s, at;
    ObjectMolecule *obj;
    int *list1 = VLAlloc(int,1000);
    int *list2 = VLAlloc(int,1000);
    int *list3 = VLAlloc(int,1000);
    int *list4 = VLAlloc(int,1000);
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int n4 = 0;
    int bonded12, bonded23, bonded34;

    /* now generate three lists of atoms, one for each selection set */

    if(list1&&list2&&list3&&list4) {

      if(just_one_object && 
         (just_one_atom[0]>=0) &&
         (just_one_atom[1]>=0) &&
         (just_one_atom[2]>=0) &&
         (just_one_atom[3]>=0)) { /* optimal case */

        list1[0] = just_one_atom[0];
        list2[0] = just_one_atom[1];
        list3[0] = just_one_atom[2];
        list4[0] = just_one_atom[3];

        n1 = n2 = n3 = n4 = 1;

      } else {
        
        for(a=cNDummyAtoms;a<I->NAtom;a++) {
          at=I->Table[a].atom;
          obj=I->Obj[I->Table[a].model];
          s=obj->AtomInfo[at].selEntry;
          if(SelectorIsMember(G,s,sele1)) {
            VLACheck(list1,int,n1);
            list1[n1++] = a;
          }
          if(SelectorIsMember(G,s,sele2)) {
            VLACheck(list2,int,n2);
            list2[n2++] = a;
          }
          if(SelectorIsMember(G,s,sele3)) {
            VLACheck(list3,int,n3);
            list3[n3++] = a;
          }
          if(SelectorIsMember(G,s,sele4)) {
            VLACheck(list4,int,n4);
            list4[n4++] = a;
          }
        }
      }
      
      /* for each set of 3 atoms in each selection... */
      
      {
        int i1,i2,i3,i4;
        int a1,a2,a3,a4;
        int at1,at2,at3,at4;
        
        /*        AtomInfoType *ai1,*ai2,ai3;*/
        CoordSet *cs1,*cs2,*cs3,*cs4;
        ObjectMolecule *obj1,*obj2,*obj3,*obj4;
        
        int idx1,idx2,idx3,idx4;
        float angle;
        float *v1,*v2,*v3,*v4, *vv0;

        for(i1=0;i1<n1;i1++) {
          a1 = list1[i1];
          at1=I->Table[a1].atom;
          obj1=I->Obj[I->Table[a1].model];

          if(state1<obj1->NCSet) {
            cs1=obj1->CSet[state1];

            if(cs1) {
              if(obj1->DiscreteFlag) {
                if(cs1==obj1->DiscreteCSet[at1]) {
                  idx1=obj1->DiscreteAtmToIdx[at1];
                } else {
                  idx1=-1;
                }
              } else {
                idx1=cs1->AtmToIdx[at1];
              }
              
              if(idx1>=0) {
                
                for(i2=0;i2<n2;i2++) {
                  a2 = list2[i2];
                  at2=I->Table[a2].atom;
                  obj2=I->Obj[I->Table[a2].model];
                  
                  if(state2<obj2->NCSet) {
                    
                    cs2=obj2->CSet[state2];
                    
                    if(cs2) {
                      if(obj2->DiscreteFlag) {
                        if(cs2==obj2->DiscreteCSet[at2]) {
                          idx2=obj2->DiscreteAtmToIdx[at2];
                        } else {
                          idx2=-1;
                        }
                      } else {
                        idx2=cs2->AtmToIdx[at2];
                      }
                    
                      if(idx2>=0) {
                        
                        bonded12 = ObjectMoleculeAreAtomsBonded2(obj1,at1,obj2,at2); 

                        if(!mode || ((mode==1)&&bonded12)) 
                         for(i3=0;i3<n3;i3++) {
                          a3 = list3[i3];
                          at3=I->Table[a3].atom;
                          obj3=I->Obj[I->Table[a3].model];
                  
                          if(state3<obj3->NCSet) {
                    
                            cs3=obj3->CSet[state3];
                    
                            if(cs3) {
                              if(obj3->DiscreteFlag) {
                                if(cs3==obj3->DiscreteCSet[at3]) {
                                  idx3=obj3->DiscreteAtmToIdx[at3];
                                } else {
                                  idx3=-1;
                                }
                              } else {
                                idx3=cs3->AtmToIdx[at3];
                              }
                    
                              if(idx3>=0) {
                                
                                bonded23 = ObjectMoleculeAreAtomsBonded2(obj2,at2,obj3,at3); 
                                if(!mode || ((mode==1)&&bonded23)) 
                                 for(i4=0;i4<n4;i4++) {
                                  a4 = list4[i4];
                          
                                  if((a1!=a2) && (a1!=a3) && (a1!=a4) && (a2!=a3) && (a2!=a4) && (a3!=a4)) {
                                    if((!((coverage[a1]==4)&&(coverage[a2]==4)&&(coverage[a3]==4)&&(coverage[a4]==4))) 
                                       ||(a1<a4)) { /* eliminate alternate-order duplicates */
                              
                                      at4=I->Table[a4].atom;
                                      obj4=I->Obj[I->Table[a4].model];
                              
                                      if(state4<obj4->NCSet) {
                                
                                        cs4=obj4->CSet[state4];
                                
                                        if(cs4) { 
                                          if(obj4->DiscreteFlag) {
                                            if(cs4==obj4->DiscreteCSet[at4]) {
                                              idx4=obj4->DiscreteAtmToIdx[at4];
                                            } else {
                                              idx4=-1;
                                            }
                                          } else {
                                            idx4=cs3->AtmToIdx[at4];
                                          }
                                  
                                          if(idx4>=0) {
                                    
                                            bonded34 = ObjectMoleculeAreAtomsBonded2(obj3,at3,obj4,at4); 
                                    
                                            if(!mode || ((mode==1)&&bonded34)) {  /* store the 3 coordinates */
                                          
                                              v1 = cs1->Coord+3*idx1;
                                              v2 = cs2->Coord+3*idx2;
                                              v3 = cs3->Coord+3*idx3;
                                              v4 = cs4->Coord+3*idx4;
                                              
                                              angle = get_dihedral3f(v1,v2,v3,v4);
                                      
                                              (*angle_sum)+=angle;
                                              (*angle_cnt)++;
                                          
                                              VLACheck(vv,float,(nv*3)+17);
                                              vv0 = vv + (nv*3);
                                              ObjectMoleculeGetAtomTxfVertex(obj1,state1,at1,vv0);
                                              ObjectMoleculeGetAtomTxfVertex(obj2,state2,at2,vv0+3);
                                              ObjectMoleculeGetAtomTxfVertex(obj3,state3,at3,vv0+6);
                                              ObjectMoleculeGetAtomTxfVertex(obj4,state4,at4,vv0+9);
                                              vv0+=12;
                                              *(vv0++) = (float)!bonded12;
                                              *(vv0++) = (float)!bonded23;
                                              *(vv0++) = (float)!bonded34;
                                              *(vv0++) = 0.0F; /* label x relative to v2+v3/2*/
                                              *(vv0++) = 0.0F; /* label y relative to v2+v3/2 */
                                              *(vv0++) = 0.0F; /* label z relative to v2+v3/2 */
                                              nv+=6;
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    VLAFreeP(list1);
    VLAFreeP(list2);
    VLAFreeP(list3);
    VLAFreeP(list4);
  }

  FreeP(coverage);
  if(vv)
    VLASize(vv,float,(nv+1)*3);
  ds->NDihedralIndex = nv;
  ds->DihedralCoord = vv;
  return(ds);
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


