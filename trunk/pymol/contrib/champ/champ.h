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

#define MAX_BOND 8
#define MAX_RING 50

/* atoms */

#define cH_Any 0xFFFFFFFF
#define cH_C   0x00000001 
#define cH_N   0x00000002 
#define cH_O   0x00000004
#define cH_H   0x00000008
#define cH_S   0x00000010
#define cH_P   0x00000020
#define cH_F   0x00000040
#define cH_Cl  0x00000080
#define cH_Br  0x00000100
#define cH_I   0x00000200
#define cH_Na  0x00000400
#define cH_K   0x00000800
#define cH_Ca  0x00001000
#define cH_Mg  0x00002000
#define cH_Zn  0x00004000
#define cH_Fe  0x00008000
#define cH_Cu  0x00010000
#define cH_Se  0x00020000
#define cH_B   0x00040000
#define cH_R0  0x00080000
#define cH_R1  0x00100000
#define cH_R2  0x00200000
#define cH_R3  0x00400000
#define cH_R4  0x00800000
#define cH_R5  0x01000000
#define cH_R6  0x02000000
#define cH_R7  0x04000000
#define cH_R8  0x08000000
#define cH_R9  0x10000000
#define cH_X   0x20000000
#define cH_Y   0x40000000
#define cH_Sym 0x80000000 /* match symbol name */

/* charge */

#define cH_Neutral    0x00000001
#define cH_Cation     0x00000002
#define cH_Dication   0x00000004
#define cH_Anion      0x00000008
#define cH_Dianion    0x00000010

/* cycles */

#define cH_Acyclic    0x00000001
#define cH_Ring3      0x00000002
#define cH_Ring4      0x00000004
#define cH_Ring5      0x00000008
#define cH_Ring6      0x00000010
#define cH_Ring7      0x00000020
#define cH_Ring8      0x00000040
#define cH_RingN      0x80000000 /* not yet implemented */

/* class */

#define cH_Aliphatic  0x00000001
#define cH_Aromatic   0x00000002

/* degree */

#define cH_0Bond       0x00000001
#define cH_1Bond       0x00000002
#define cH_2Bond       0x00000004
#define cH_3Bond       0x00000008
#define cH_4Bond       0x00000010
#define cH_5Bond       0x00000020
#define cH_6Bond       0x00000040

/* total valence */

#define cH_0Valence       0x00000001
#define cH_1Valence       0x00000002
#define cH_2Valence       0x00000004
#define cH_3Valence       0x00000008
#define cH_4Valence       0x00000010
#define cH_5Valence       0x00000020
#define cH_6Valence       0x00000040

typedef char AtomBuffer[255]; /* maximum length for a single atom */

#define NAM_SIZE 4
#define RES_SIZE 4
#define SYM_SIZE 2

typedef struct { 
  int link; /* memory management */
  int bond[MAX_BOND+1]; 
  int pos_flag;
  int atom;
  int charge;
  int cycle;
  int class;
  int degree;
  int valence;
  char symbol[SYM_SIZE];
  char name[NAM_SIZE];
  char residue[RES_SIZE];
  int neg_flag;
  int not_atom;
  int not_charge;
  int not_cycle;
  int not_class;
  int not_degree;
  int not_valence;
  int tag; /* string index for tag */
  int mark_tmpl,mark_targ; /* traversal */
  int first_tmpl,first_targ; /* first template stack entry */
  PyObject *chempy_atom;
} ListAtom;

/* order */

#define cH_Single      0x00000001
#define cH_Double      0x00000002
#define cH_Triple      0x00000004

typedef struct {
  int link;     /* memory management */
  int atom[2];  /* connected atoms  -- directionality must reflect tree */
  int order;
  int class;
  int cycle; 
  int not_order;
  int not_class;
  int not_cycle;
  int tag; /* string index for tag */
  int mark_tmpl,mark_targ; /* traversal */
  PyObject *chempy_bond;
} ListBond;


typedef struct {
  int link;
  int atom;
  int bond; /* bond index index in atom */
  int base_bond; /* absolute bond index */
  int paren_flag;
} ListScope;

typedef struct {
  int link;
  int atom; /* root of atom list (Pat) */ 
  int bond; /* root of bond list (Bond) */
  PyObject *chempy_molecule;
  int unique_atom; /* list of unique atoms (Int) */
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
} ListTmpl;

typedef struct {
  int link;
  int atom; /* root of atom list (Pat) */ 
  int bond; /* root of bond list (Bond) */
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
} CChamp; /* scream class */

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
int ChampSSS_1V1_B(CChamp *I,int pattern,int target);
int ChampSSS_1VN_N(CChamp *I,int pattern,int list);
int ChampModelToPat(CChamp *I,PyObject *model);
char *ChampPatToSmiVLA(CChamp *I,int index);

#endif









