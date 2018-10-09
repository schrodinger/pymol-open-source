/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright Schrodinger, LLC.
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

#ifndef _H_AtomInfoHistory
#define _H_AtomInfoHistory

#include <map>

#include"AtomInfo.h"
#include"Util.h"
#include"Lex.h"

typedef struct AtomInfoType_1_7_6 {
  int resv;
  int customType;
  int priority;
  float b, q, vdw, partialCharge;
  int selEntry;
  int color;
  int id;                       // PDB ID
  unsigned int flags;
  int temp1;                    /* kludge fields - to remove */
  int unique_id;                /* introduced in version 0.77 */
  int discrete_state;           /* state+1 for atoms in discrete objects */
  float elec_radius;            /* radius for PB calculations */
  int rank;
  int textType;
  int custom;
  int label;
  int visRep;                   /* bitmask for all reps */

  /* be careful not to write at these as (int*) */

  signed char formalCharge;     // values typically in range -2..+2
  signed char stereo;           /* for 2D representation */
  signed char mmstereo;           /* from MMStereo */
  signed char cartoon;          /* 0 = default which is auto (use ssType) */

  // boolean flags
  signed char hetatm;
  signed char bonded;
  signed char chemFlag;         // 0,1,2
  signed char geom;             // cAtomInfo*

  signed char valence;

  // boolean flags
  signed char deleteFlag;
  signed char masked;

  signed char protekted;        // 0,1,2

  signed char protons;          /* atomic number */

  // boolean flags
  signed char hb_donor;
  signed char hb_acceptor;
  signed char has_setting;      /* setting based on unique_id */

  int chain;
  SegIdent segi;
  AtomName name;
  ElemName elem;                // redundant with "protons" ?
  ResIdent resi;
  char has_prop;
  SSType ssType;                /* blank or 'L' = turn/loop, 'H' = helix, 'S' = beta-strand/sheet */
  Chain alt;
  ResName resn;

  // replace with pointer?
  float U11, U22, U33, U12, U13, U23;

  int oldid; // for undo

  int prop_id;
  bool has_anisou() const { return U11 || U22 || U33 || U12 || U13 || U23; }
  const float * get_anisou() const { return &U11; } // for copying
  float * get_anisou() { return &U11; } // for copying

  inline void setResi(int resv_, char inscode_) {
    resv = resv_;
    AtomResiFromResv(resi, sizeof(resi), resv, inscode_);
  }

  inline char getInscode() const {
    int i = strlen(resi) - 1;
    if (i >= 0 && (resi[i] < '0' || '9' < resi[i]))
      return resi[i];
    return '\0';
  }
} AtomInfoType_1_7_6;

typedef struct AtomInfoType_1_7_7 {
  union {
    float * anisou;               // only allocate with get_anisou
    int64_t dummyanisou;
  };
  int resv;
  int customType;
  int priority;
  float b, q, vdw, partialCharge;
  int selEntry;
  int color;
  int id;                       // PDB ID
  unsigned int flags;
  int temp1;                    /* kludge fields - to remove */
  int unique_id;                /* introduced in version 0.77 */
  int discrete_state;           /* state+1 for atoms in discrete objects */
  float elec_radius;            /* radius for PB calculations */
  int rank;
  int textType;
  int custom;
  int label;
  int visRep;                   /* bitmask for all reps */
  int oldid;                    // for undo
  int prop_id;

  // boolean flags
  bool hetatm : 1;
  bool bonded : 1;
  bool deleteFlag : 1;
  bool masked : 1;
  bool hb_donor : 1;
  bool hb_acceptor : 1;
  bool has_setting : 1;      /* setting based on unique_id */
  bool has_prop : 1;

  /* be careful not to write at these as (int*) */

  signed char formalCharge;     // values typically in range -2..+2
  signed char mmstereo;           /* from MMStereo */
  signed char cartoon;          /* 0 = default which is auto (use ssType) */
  signed char geom;             // cAtomInfo*
  signed char valence;          // 0-4
  signed char protons;          /* atomic number */

  int chain;
  SegIdent segi;
  AtomName name;
  ElemName elem;                // redundant with "protons" ?
  ResIdent resi;
  SSType ssType;                /* blank or 'L' = turn/loop, 'H' = helix, 'S' = beta-strand/sheet */
  Chain alt;
  ResName resn;

  // small value optimized bitfields
  unsigned char stereo : 2;     // 0-3 Only for SDF (MOL) format in/out
  unsigned char chemFlag : 2;   // 0,1,2
  unsigned char protekted : 2;  // 0,1,2

  // no anisou support in 1.8.0
  bool has_anisou() const { return false; }
  const float * get_anisou() const { return NULL; }
  float * get_anisou() { return NULL; }

  inline void setResi(int resv_, char inscode_) {
    resv = resv_;
    AtomResiFromResv(resi, sizeof(resi), resv, inscode_);
  }

  inline char getInscode() const {
    int i = strlen(resi) - 1;
    if (i >= 0 && (resi[i] < '0' || '9' < resi[i]))
      return resi[i];
    return '\0';
  }
} AtomInfoType_1_7_7;

/*
 * This is not identical to the 1.8.2 AtomInfoType, it's missing all members
 * which are not relevant or unsupported with pse_binary_dump (anisou,
 * selEntry, temp1, oldid, prop_id, deleteFlag)
 */
typedef int lexidx_int_t;
struct AtomInfoType_1_8_1 {
  // contiguous memory instead of pointer ("short" matches PDB precision)
  short anisou[6];

  lexidx_int_t segi;
  lexidx_int_t chain;
  lexidx_int_t resn;
  lexidx_int_t name;
  lexidx_int_t textType;
  lexidx_int_t custom;
  lexidx_int_t label;

  int resv;
  int customType;
  int priority;
  float b, q, vdw, partialCharge;
  int color;
  int id;                       // PDB ID
  unsigned int flags;
  int unique_id;                /* introduced in version 0.77 */
  int discrete_state;           /* state+1 for atoms in discrete objects */
  float elec_radius;            /* radius for PB calculations */
  int rank;
  int visRep;                   /* bitmask for all reps */

  // boolean flags
  bool hetatm : 1;
  bool bonded : 1;
  bool masked : 1;
  bool hb_donor : 1;
  bool hb_acceptor : 1;
  bool has_setting : 1;      /* setting based on unique_id */

  signed char formalCharge;     // values typically in range -2..+2
  signed char cartoon;          /* 0 = default which is auto (use ssType) */
  signed char geom;             // cAtomInfo*
  signed char valence;          // 0-4
  signed char protons;          /* atomic number */

  char inscode;

  ElemName elem;               // redundant with "protons" ?
  SSType ssType;               /* blank or 'L' = turn/loop, 'H' = helix, 'S' = beta-strand/sheet */
  Chain alt;

  // small value optimized bitfields
  unsigned char stereo : 2;     // 0-3 Only for SDF (MOL) format in/out
  unsigned char chemFlag : 2;   // 0,1,2
  unsigned char protekted : 2;  // 0,1,2

  char getInscode() const { return inscode; }
  void setInscode(char c) { inscode = makeInscode(c); }

  void setResi(const char * resi) {
    if (sscanf(resi, "%d%c", &resv, &inscode) == 1 || inscode <= ' ')
      inscode = '\0';
  }

  void setResi(int resv_, char inscode_) {
    resv = resv_;
    setInscode(inscode_);
  }

  bool has_anisou() const { return anisou[0] || anisou[1] || anisou[2] || anisou[3] || anisou[4] || anisou[5]; }
  const short * get_anisou() const { return anisou; } // for copying
  short * get_anisou() { return anisou; } // for copying
};

class AtomInfoTypeConverter {
  PyMOLGlobals * G;
  int NAtom;

  template <typename D, typename S> void copy1(D *dest, const S *src);
  template <typename D, typename S> void copyN(D *dest, const S *src);
  template <typename D> D * allocCopy(const AtomInfoType *src);

public:
  AtomInfoTypeConverter(PyMOLGlobals * G_, int natom) : G(G_), NAtom(natom) {}

  std::map<lexidx_int_t, lexidx_t> lexidxmap;

  lexidx_int_t to_lexidx_int(const lexidx_t& idx) {
    return idx;
  }

  void copy(AtomInfoType * dest, const void *src, int srcversion);
  void * allocCopy(int destversion, const AtomInfoType * src);

  /*
   * For copying Lex strings
   */
  inline void copy_attr_s(lexidx_t& dest, lexidx_t src) {
    if (!lexidxmap.empty())
      src = lexidxmap[src];
    LexAssign(G, dest, src);
  }
  inline void copy_attr_s(lexidx_t& dest, const char * src) {
    LexAssign(G, dest, src);
  }
  template <size_t N> inline void copy_attr_s(char (&dest)[N], const lexidx_t& src) {
    UtilNCopy(dest, LexStr(G, src), sizeof(dest));
  }
};

#endif

