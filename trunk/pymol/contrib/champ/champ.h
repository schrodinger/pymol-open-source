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
#ifndef _H_Champ
#define _H_Champ

#include"list.h"
#include"os_python.h"

/* max bonds an atom can have */

#define MAX_BOND 12
#define MAX_RING 50

/* CUTOFF = 2*N/2*N+1 for the largest ring we spend time to detecting:
   8 = 16-17 atom cycles 
 */

#define RING_SEARCH_CUTOFF 8

/* atoms */

#define cH_Any  0xFFFFFFFF
#define cH_H    0x00000001
#define cH_C    0x00000002 
#define cH_N    0x00000004 
#define cH_O    0x00000008
#define cH_Sym  0x00000010
#define cH_S    0x00000020
#define cH_P    0x00000040
#define cH_F    0x00000080
#define cH_Cl   0x00000100
#define cH_Br   0x00000200
#define cH_I    0x00000400
#define cH_Na   0x00000800
#define cH_K    0x00001000
#define cH_Ca   0x00002000
#define cH_Mg   0x00004000
#define cH_Zn   0x00008000
#define cH_Fe   0x00010000
#define cH_Cu   0x00020000
#define cH_Se   0x00040000
#define cH_B    0x00080000
#define cH_A    0x00100000
#define cH_E    0x00200000
#define cH_G    0x00400000
#define cH_J    0x00800000
#define cH_L    0x01000000
#define cH_M    0x02000000
#define cH_Q    0x04000000
#define cH_R    0x08000000
#define cH_T    0x10000000
#define cH_X    0x20000000
#define cH_Z    0x40000000

#define cH_NotH 0xFFFFFFFE

/* charge */

#define cH_Neutral     0x00000001
#define cH_Cation      0x00000002
#define cH_Dication    0x00000004
#define cH_Anion       0x00000008
#define cH_Dianion     0x00000010
#define cH_Trication   0x00000020
#define cH_Trianion    0x00000040
#define cH_Tetcation   0x00000080
#define cH_Tetanion    0x00000100
#define cH_Pentcation  0x00000200
#define cH_Pentanion   0x00000400

/* cycles */

#define cH_Acyclic    0x00000001
#define cH_Ring3      0x00000002
#define cH_Ring4      0x00000004
#define cH_Ring5      0x00000008
#define cH_Ring6      0x00000010
#define cH_Ring7      0x00000020
#define cH_Ring8      0x00000040
#define cH_RingN      0x80000000 /* not yet implemented */

#define cH_Cyclic     0xFFFFFFFE

/* class */

#define cH_Aliphatic  0x00000001
#define cH_Aromatic   0x00000002
#define cH_Pi         0x00000004 /* non-exclusive with above */

#define cH_AnyClass   0x00000003

/* degree */

#define cH_0Bond       0x00000001
#define cH_1Bond       0x00000002
#define cH_2Bond       0x00000004
#define cH_3Bond       0x00000008
#define cH_4Bond       0x00000010
#define cH_5Bond       0x00000020
#define cH_6Bond       0x00000040
#define cH_7Bond       0x00000080
#define cH_8Bond       0x00000100

/* valence */

#define cH_0Valence       0x00000001
#define cH_1Valence       0x00000002
#define cH_2Valence       0x00000004
#define cH_3Valence       0x00000008
#define cH_4Valence       0x00000010
#define cH_5Valence       0x00000020
#define cH_6Valence       0x00000040
#define cH_7Valence       0x00000080
#define cH_8Valence       0x00000100

typedef char AtomBuffer[255]; /* maximum length for a single atom */

#define NAM_SIZE 5
#define RES_SIZE 5
#define SYM_SIZE 3

#define dTag_disable 0
#define cTag_merge   1
#define cTag_copy    2

typedef struct { 
  int link; /* memory management */
  int index;
  int bond[MAX_BOND+1]; 
  int pos_flag;
  int atom;
  int charge;
  int cycle;
  int class;
  int degree;
  int valence;
  int imp_hydro;
  int tot_hydro; 
  int hydro_flag; /* are we trying to match hydrogen counts? */
  char symbol[SYM_SIZE];
  char name[NAM_SIZE];
  char residue[RES_SIZE];
  float coord[3];
  int neg_flag;
  int not_atom;
  int not_charge;
  int not_cycle;
  int not_class;
  int not_degree;
  int not_valence;
  int comp_imp_hydro_flag; /* do we need to compute implicit hydrogens? */
  int stereo; /* based on lexical ordering: 0 = unspecified, 1 = anti-clockwise, -1 = clockwise */
  int stereo_internal; /* based on internal atom IDs */
  int mark_tmpl,mark_targ,mark_read; /* traversal */
  int first_tmpl,first_targ; /* first template stack entry */
  int first_base; 
  int ext_index; /* for maintaining atom references with PyMOL */
  unsigned int tag,not_tag;
  PyObject *chempy_atom;
} ListAtom;

/* order */

#define cH_Single      0x00000001
#define cH_Double      0x00000002
#define cH_Triple      0x00000004
#define cH_AnyOrder    0x00000007
#define cH_NoOrder     0x00000000

/* double-bond stereochem */

#define cH_Up           1
#define cH_Down        -1

/* tetrahedral stereochem */

#define cH_Anticlock    1
#define cH_Clockwise   -1

typedef struct {
  int link;     /* memory management */
  int index;
  int atom[2];  /* connected atoms  -- directionality must reflect tree */
  int pri[2]; /* forward and backward lexical priorities */
  int order;
  int class;
  int cycle; 
  int not_order;
  int not_class;
  int not_cycle;
  int direction; /* 0 = specified, 1 = up, -1 = down */
  int mark_tmpl,mark_targ,mark_read; /* traversal */
  unsigned int tag,not_tag;
  int ext_index;
  PyObject *chempy_bond;
} ListBond;

typedef struct {
  int link;
  int atom;
  int bond; /* bond index index in atom */
  int base_bond; /* absolute bond index */
  int base_atom; /* absolute atom index */
  int paren_flag;
} ListScope;

typedef struct {
  int link;
  int atom; /* root of atom list (Pat) */ 
  int bond; /* root of bond list (Bond) */
  PyObject *chempy_molecule;
  int unique_atom; /* list of unique atoms (Int) */
  int target_prep; /* has pattern been prepared as a target? */
} ListPat;

typedef struct {
  int link;
  int atom;
  int bond;
} ListMatch;

typedef struct {
  int link;
  int atom; /* root of atom list (Pat) */ 
  int bond; /* root of bond list (Bond) */
  int parent;
  int match;
  int targ_start;
  int bond_pri;
} ListTmpl;

typedef struct {
  int link;
  int atom; /* root of atom list (Pat) */ 
  int bond; /* root of bond list (Bond) */
  int bond_pri;
} ListTarg;

typedef struct {
  ListAtom  *Atom; /* lists */
  ListBond  *Bond;
  ListInt   *Int;
  ListInt2  *Int2;
  ListInt3  *Int3;
  ListTmpl  *Tmpl;
  ListTarg  *Targ;
  ListPat   *Pat;
  ListScope *Scope;
  ListMatch *Match;
  char *Str;
  int ActivePatList;
} CChamp; /* champ class */

/* prototypes */


CChamp *ChampNew(void);
void ChampFree(CChamp *I);

void ChampAtomFree(CChamp *I,int atom);
void ChampAtomFreeChain(CChamp *I,int atom);

void ChampBondFree(CChamp *I,int bond);
void ChampBondFreeChain(CChamp *I,int bond);

void ChampPatFree(CChamp *I,int index);

int ChampPatIdentical(ListAtom *p,ListAtom *a);
int ChampAtomMatch(ListAtom *p,ListAtom *a);
int ChampBondMatch(ListBond *p,ListBond *a);

int ChampSmiToPat(CChamp *I,char *c);
void ChampMemoryDump(CChamp *I);
int ChampMemoryUsage(CChamp *I);
int ChampMatch_1V1_B(CChamp *I,int pattern,int target);
int ChampMatch_1V1_Map(CChamp *I,int pattern,int target,int limit,int tag_flag);
int ChampMatch_1VN_N(CChamp *I,int pattern,int list);
int ChampMatch_1V1_N(CChamp *I,int pattern,int target,int limit,int tag_flag);

int ChampMatch_NV1_N(CChamp *I,int list,int target,int limit,int tag_flag);
int ChampExact_1VN_N(CChamp *I,int pattern,int list);

int ChampModelToPat(CChamp *I,PyObject *model);
char *ChampPatToSmiVLA(CChamp *I,int index,char *vla,int mode);
int ChampAtomToString(CChamp *I,int index,char *buf);
int ChampBondToString(CChamp *I,int index,char *buf);
void ChampPatReindex(CChamp *I,int index);
void ChampPatDump(CChamp *I,int index);
void ChampOrientBonds(CChamp *I,int index);
void ChampGeneralize(CChamp *I,int index);
void ChampDetectChirality(CChamp *I,int index);

#endif









