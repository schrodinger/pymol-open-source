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

#ifndef _H_BondTypeHistory
#define _H_BondTypeHistory

#include"AtomInfo.h"
#include"SymOp.h"

typedef struct BondType_1_7_6 {
  int index[2];
  int order;
  int id;
  int unique_id;
  int temp1;
  short int stereo;             /* to preserve 2D rep */
  short int has_setting;        /* setting based on unique_id */
  int oldid;
} BondType_1_7_6;

typedef struct BondType_1_7_7 {
  int index[2];
  int id;
  int unique_id;
  int oldid;
  signed char order;    // 0-4
  signed char temp1;    // bool? where used?
  signed char stereo;   // 0-6 Only for SDF (MOL) format in/out
  bool has_setting;     /* setting based on unique_id */
} BondType_1_7_7;

/*
 * This is not identical to the 1.8.2 BondType, it's missing all members
 * which are not relevant or unsupported with pse_binary_dump (oldid, temp1)
 */
struct BondType_1_8_1 {
  int index[2];
  union {
    int id;
    pymol::SymOp symop_2;
  };
  int unique_id;
  signed char order;    // 0-4
  signed char stereo;   // 0-6 Only for SDF (MOL) format in/out
  bool has_setting;     /* setting based on unique_id */
};

static_assert(sizeof(pymol::SymOp) == 4,
    "Changing SymOp layout requires a backwards "
    "compatibility strategy for pse_binary_dump");

static_assert(sizeof(BondType_1_8_1) == 20,
    "pse_binary_dump not compatible with this platform");

void Copy_Into_BondType_From_Version(const void *binstr, int bondInfo_version, BondType *Bond, int NBond);
void *Copy_To_BondType_Version(int bondInfo_version, BondType *Bond, int NBond);
#endif
