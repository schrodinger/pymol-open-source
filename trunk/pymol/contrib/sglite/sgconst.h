/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#ifndef SGCONST_H__
#define SGCONST_H__


#ifndef SG_GLOBAL
extern const int EV_100[];
extern const int EV_010[];
extern const int EV_001[];
extern const int EV_m10[];
extern const int EV_110[];
extern const int EV_111[];
#else
const int EV_100[] = {  1,  0,  0 };
const int EV_010[] = {  0,  1,  0 };
const int EV_001[] = {  0,  0,  1 };
const int EV_m10[] = { -1,  1,  0 };
const int EV_110[] = {  1,  1,  0 };
const int EV_111[] = {  1,  1,  1 };
#endif /* SG_GLOBAL */


#ifndef SG_GLOBAL
extern const int  R_1_000[];
extern const int  R_2_001[];
extern const int  R_2_1b0[];
extern const int  R_2_110[];
extern const int  R_3_001[];
extern const int  R_3_111[];
extern const int  R_3i111[];
extern const int  R_4_001[];
extern const int  R_4i001[];
extern const int  R_6_001[];
#else
const int  R_1_000[] =
{  1,  0,  0,
   0,  1,  0,
   0,  0,  1 };
const int  R_2_001[] =
{ -1,  0,  0,
   0, -1,  0,
   0,  0,  1 };
const int  R_2_1b0[] =
{  0, -1,  0,
  -1,  0,  0,
   0,  0, -1 };
const int  R_2_110[] =
{  0,  1,  0,
   1,  0,  0,
   0,  0, -1 };
const int  R_3_001[] =
{  0, -1,  0,
   1, -1,  0,
   0,  0,  1 };
const int  R_3_111[] =
{  0,  0,  1,
   1,  0,  0,
   0,  1,  0 };
const int  R_3i111[] =
{  0,  1,  0,
   0,  0,  1,
   1,  0,  0 };
const int  R_4_001[] =
{  0, -1,  0,
   1,  0,  0,
   0,  0,  1 };
const int  R_4i001[] =
{  0,  1,  0,
  -1,  0,  0,
   0,  0,  1 };
const int  R_6_001[] =
{  1, -1,  0,
   1,  0,  0,
   0,  0,  1 };
#endif /* SG_GLOBAL */


#ifndef SG_GLOBAL
extern const T_RTMx  CBMx_1_000[1];
extern const T_RTMx  CBMx_2_110[1];
extern const T_RTMx  CBMx_2_101[1];
extern const T_RTMx  CBMx_3_010[1];
extern const T_RTMx  CBMx_3_111[1];
extern const T_RTMx  CBMx_4_001[2];
extern const T_RTMx  CBMxMon_c_b[2];
extern const T_RTMx  CBMxCP[2];
extern const T_RTMx  CBMxFI[2];
extern const T_RTMx  CBMxRevObv[2];
extern const T_RTMx  CBMxHP[2];
#else
#define R(i) ((i) * (CRBF / 12))
const T_RTMx  CBMx_1_000[1] =
{{{{ R( 12), R(  0), R(  0),
     R(  0), R( 12), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMx_2_110[1] =
{{{{ R(  0), R( 12), R(  0),
     R( 12), R(  0), R(  0),
     R(  0), R(  0), R(-12) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMx_2_101[1] =
{{{{ R(  0), R(  0), R( 12),
     R(  0), R(-12), R(  0),
     R( 12), R(  0), R(  0) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMx_3_010[1] =
{{{{ R(-12), R(  0), R( 12),
     R(  0), R( 12), R(  0),
     R(-12), R(  0), R(  0) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMx_3_111[1] =
{{{{ R(  0), R(  0), R( 12),
     R( 12), R(  0), R(  0),
     R(  0), R( 12), R(  0) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMx_4_001[2] =
{{{{ R(  0), R(-12), R(  0),
     R( 12), R(  0), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }},
 {{{ R(  0), R( 12), R(  0),
     R(-12), R(  0), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMxMon_c_b[2] =
{{{{ R(  0), R( 12), R(  0),
     R(  0), R(  0), R( 12),
     R( 12), R(  0), R(  0) }, { 0, 0, 0 }
 }},
 {{{ R(  0), R(  0), R( 12),
     R( 12), R(  0), R(  0),
     R(  0), R( 12), R(  0) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMxCP[2] =
{{{{ R( 12), R( 12), R(  0),
     R( 12), R(-12), R(  0),
     R(  0), R(  0), R(-12) }, { 0, 0, 0 }
 }},
 {{{ R(  6), R(  6), R(  0),
     R(  6), R( -6), R(  0),
     R(  0), R(  0), R(-12) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMxFI[2] =
{{{{ R( 12), R( 12), R(  0),
     R(-12), R( 12), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }},
 {{{ R(  6), R( -6), R(  0),
     R(  6), R(  6), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMxRevObv[2] =
{{{{ R(-12), R(  0), R(  0),
     R(  0), R(-12), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }},
 {{{ R(-12), R(  0), R(  0),
     R(  0), R(-12), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }}
};
const T_RTMx  CBMxHP[2] =
{{{{ R( 12), R( 12), R(  0),
     R(-12), R( 24), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }},
 {{{ R(  8), R( -4), R(  0),
     R(  4), R(  4), R(  0),
     R(  0), R(  0), R( 12) }, { 0, 0, 0 }
 }}
};
#undef R
#endif /* SG_GLOBAL */


#define XS_Undefined     0
#define XS_Unknown       1
#define XS_Triclinic     2
#define XS_Monoclinic    3
#define XS_Orthorhombic  4
#define XS_Tetragonal    5
#define XS_Trigonal      6
#define XS_Hexagonal     7
#define XS_Cubic         8

#ifndef SG_GLOBAL
extern
const char *XS_Name[];
#else
const char *XS_Name[] = {
  "Undefined",
  "Unknown",
  "Triclinic",
  "Monoclinic",
  "Orthorhombic",
  "Tetragonal",
  "Trigonal",
  "Hexagonal",
  "Cubic"
};
#endif


#define Make_MGC(XS, L, P, M) (((XS * 9 + L) * 3 + P) * 42 + M)

#define MGC_Undefined Make_MGC(XS_Undefined,     0,  0,  0)
#define MGC_Unknown   Make_MGC(XS_Unknown,       0,  0,  1)
#define MGC_1         Make_MGC(XS_Triclinic,     1,  0,  2)
#define MGC_1b        Make_MGC(XS_Triclinic,     0,  0,  3)
#define MGC_2         Make_MGC(XS_Monoclinic,    2,  0,  4)
#define MGC_m         Make_MGC(XS_Monoclinic,    1,  0,  5)
#define MGC_2_m       Make_MGC(XS_Monoclinic,    0,  0,  6)
#define MGC_222       Make_MGC(XS_Orthorhombic,  2,  0,  7)
#define MGC_mm2       Make_MGC(XS_Orthorhombic,  1,  0,  8)
#define MGC_mmm       Make_MGC(XS_Orthorhombic,  0,  0,  9)
#define MGC_4         Make_MGC(XS_Tetragonal,    2,  0, 10)
#define MGC_4b        Make_MGC(XS_Tetragonal,    1,  0, 11)
#define MGC_4_m       Make_MGC(XS_Tetragonal,    0,  0, 12)
#define MGC_422       Make_MGC(XS_Tetragonal,    4,  0, 13)
#define MGC_4mm       Make_MGC(XS_Tetragonal,    3,  0, 14)
#define MGC_4b2m      Make_MGC(XS_Tetragonal,    2,  1, 15)
#define MGC_4bm2      Make_MGC(XS_Tetragonal,    1,  0, 16)
#define MGC_4_mmm     Make_MGC(XS_Tetragonal,    0,  0, 17)
#define MGC_3         Make_MGC(XS_Trigonal,      1,  0, 18)
#define MGC_3b        Make_MGC(XS_Trigonal,      0,  0, 19)
#define MGC_321       Make_MGC(XS_Trigonal,      8,  2, 20)
#define MGC_312       Make_MGC(XS_Trigonal,      7,  1, 21)
#define MGC_32        Make_MGC(XS_Trigonal,      6,  0, 22)
#define MGC_3m1       Make_MGC(XS_Trigonal,      5,  2, 23)
#define MGC_31m       Make_MGC(XS_Trigonal,      4,  1, 24)
#define MGC_3m        Make_MGC(XS_Trigonal,      3,  0, 25)
#define MGC_3bm1      Make_MGC(XS_Trigonal,      2,  2, 26)
#define MGC_3b1m      Make_MGC(XS_Trigonal,      1,  1, 27)
#define MGC_3bm       Make_MGC(XS_Trigonal,      0,  0, 28)
#define MGC_6         Make_MGC(XS_Hexagonal,     2,  0, 29)
#define MGC_6b        Make_MGC(XS_Hexagonal,     1,  0, 30)
#define MGC_6_m       Make_MGC(XS_Hexagonal,     0,  0, 31)
#define MGC_622       Make_MGC(XS_Hexagonal,     4,  0, 32)
#define MGC_6mm       Make_MGC(XS_Hexagonal,     3,  0, 33)
#define MGC_6b2m      Make_MGC(XS_Hexagonal,     2,  1, 34)
#define MGC_6bm2      Make_MGC(XS_Hexagonal,     1,  0, 35)
#define MGC_6_mmm     Make_MGC(XS_Hexagonal,     0,  0, 36)
#define MGC_23        Make_MGC(XS_Cubic,         1,  0, 37)
#define MGC_m3b       Make_MGC(XS_Cubic,         0,  0, 38)
#define MGC_432       Make_MGC(XS_Cubic,         2,  0, 39)
#define MGC_4b3m      Make_MGC(XS_Cubic,         1,  0, 40)
#define MGC_m3bm      Make_MGC(XS_Cubic,         0,  0, 41)

#define ixMG(MGC)  ((MGC) % 42)
#define ixPG(MGC) (((MGC) % 42) + (((MGC) / 42) % 3))
#define ixLG(MGC) (((MGC) % 42) + (((MGC) / (42 * 3)) % 9))
#define ixXS(MGC)  ((MGC) / (42 * 3 * 9))

#ifndef SG_GLOBAL
extern
const char *MG_Names[];
#else
const char *MG_Names[] = {
  "Undefined",
  "Unknown",
  "1",
  "-1",
  "2",
  "m",
  "2/m",
  "222",
  "mm2",
  "mmm",
  "4",
  "-4",
  "4/m",
  "422",
  "4mm",
  "-4m2",
  "-42m",
  "4/mmm",
  "3",
  "-3",
  "321",
  "312",
  "32",
  "3m1",
  "31m",
  "3m",
  "-3m1",
  "-31m",
  "-3m",
  "6",
  "-6",
  "6/m",
  "622",
  "6mm",
  "-6m2",
  "-62m",
  "6/mmm",
  "23",
  "m-3",
  "432",
  "-43m",
  "m-3m"
};
#endif /* SG_GLOBAL */


#endif /* SGCONST_H__ */
