/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#undef SG_GLOBAL
#include "sglite.h"
#include "sgconst.h"
#include "sgrefset.h"


static const char *Monoclinic_SgNumber_as_HM_List[][2] = {
         { NULL, NULL },
         { NULL, NULL },
         { NULL, NULL },
/*  3 */ { "P121",    "P112" },
/*  4 */ { "P1211",   "P1121" },
/*  5 */ { "C121",    "B112" },
/*  6 */ { "P1m1",    "P11m" },
/*  7 */ { "P1c1",    "P11b" },
/*  8 */ { "C1m1",    "B11m" },
/*  9 */ { "C1c1",    "B11b" },
/* 10 */ { "P12/m1",  "P112/m" },
/* 11 */ { "P121/m1", "P1121/m" },
/* 12 */ { "C12/m1",  "B112/m" },
/* 13 */ { "P12/c1",  "P112/b" },
/* 14 */ { "P121/c1", "P1121/b" },
/* 15 */ { "C12/c1",  "B112/b" }
};


static const char *Schoenflies_List[] = {
NULL,
"C1^1", "Ci^1", "C2^1", "C2^2", "C2^3", "Cs^1", "Cs^2", "Cs^3", "Cs^4",
"C2h^1", "C2h^2", "C2h^3", "C2h^4", "C2h^5", "C2h^6", "D2^1", "D2^2",
"D2^3", "D2^4", "D2^5", "D2^6", "D2^7", "D2^8", "D2^9", "C2v^1",
"C2v^2", "C2v^3", "C2v^4", "C2v^5", "C2v^6", "C2v^7", "C2v^8",
"C2v^9", "C2v^10", "C2v^11", "C2v^12", "C2v^13", "C2v^14", "C2v^15",
"C2v^16", "C2v^17", "C2v^18", "C2v^19", "C2v^20", "C2v^21", "C2v^22",
"D2h^1", "D2h^2", "D2h^3", "D2h^4", "D2h^5", "D2h^6", "D2h^7", "D2h^8",
"D2h^9", "D2h^10", "D2h^11", "D2h^12", "D2h^13", "D2h^14", "D2h^15",
"D2h^16", "D2h^17", "D2h^18", "D2h^19", "D2h^20", "D2h^21", "D2h^22",
"D2h^23", "D2h^24", "D2h^25", "D2h^26", "D2h^27", "D2h^28", "C4^1",
"C4^2", "C4^3", "C4^4", "C4^5", "C4^6", "S4^1", "S4^2", "C4h^1",
"C4h^2", "C4h^3", "C4h^4", "C4h^5", "C4h^6", "D4^1", "D4^2", "D4^3",
"D4^4", "D4^5", "D4^6", "D4^7", "D4^8", "D4^9", "D4^10", "C4v^1",
"C4v^2", "C4v^3", "C4v^4", "C4v^5", "C4v^6", "C4v^7", "C4v^8",
"C4v^9", "C4v^10", "C4v^11", "C4v^12", "D2d^1", "D2d^2", "D2d^3",
"D2d^4", "D2d^5", "D2d^6", "D2d^7", "D2d^8", "D2d^9", "D2d^10",
"D2d^11", "D2d^12", "D4h^1", "D4h^2", "D4h^3", "D4h^4", "D4h^5",
"D4h^6", "D4h^7", "D4h^8", "D4h^9", "D4h^10", "D4h^11", "D4h^12",
"D4h^13", "D4h^14", "D4h^15", "D4h^16", "D4h^17", "D4h^18", "D4h^19",
"D4h^20", "C3^1", "C3^2", "C3^3", "C3^4", "C3i^1", "C3i^2", "D3^1",
"D3^2", "D3^3", "D3^4", "D3^5", "D3^6", "D3^7", "C3v^1", "C3v^2",
"C3v^3", "C3v^4", "C3v^5", "C3v^6", "D3d^1", "D3d^2", "D3d^3",
"D3d^4", "D3d^5", "D3d^6", "C6^1", "C6^2", "C6^3", "C6^4", "C6^5",
"C6^6", "C3h^1", "C6h^1", "C6h^2", "D6^1", "D6^2", "D6^3", "D6^4",
"D6^5", "D6^6", "C6v^1", "C6v^2", "C6v^3", "C6v^4", "D3h^1", "D3h^2",
"D3h^3", "D3h^4", "D6h^1", "D6h^2", "D6h^3", "D6h^4", "T^1", "T^2",
"T^3", "T^4", "T^5", "Th^1", "Th^2", "Th^3", "Th^4", "Th^5", "Th^6",
"Th^7", "O^1", "O^2", "O^3", "O^4", "O^5", "O^6", "O^7", "O^8",
"Td^1", "Td^2", "Td^3", "Td^4", "Td^5", "Td^6", "Oh^1", "Oh^2",
"Oh^3", "Oh^4", "Oh^5", "Oh^6", "Oh^7", "Oh^8", "Oh^9", "Oh^10"
};


typedef struct {
  const char  *shrt;
  const char  *full;
}
T_Short_Mono_HM_Dict;

static const
T_Short_Mono_HM_Dict VolI_Short_Mono_HM_Dict[] = {
{ "P2",    "P112" },
{ "P21",   "P1121" },
{ "B2",    "B112" },
{ "C2",    "C121" },
{ "Pm",    "P11m" },
{ "Pb",    "P11b" },
{ "Pc",    "P1c1" },
{ "Bm",    "B11m" },
{ "Cm",    "C1m1" },
{ "Bb",    "B11b" },
{ "Cc",    "C1c1" },
{ "P2/m",  "P112/m" },
{ "P21/m", "P1121/m" },
{ "B2/m",  "B112/m" },
{ "C2/m",  "C12/m1" },
{ "P2/b",  "P112/b" },
{ "P2/c",  "P12/c1" },
{ "P21/b", "P1121/b" },
{ "P21/c", "P121/c1" },
{ "B2/b",  "B112/b" },
{ "C2/c",  "C12/c1" },
{ NULL, NULL }
};

static const
T_Short_Mono_HM_Dict VolA_Short_Mono_HM_Dict[] = {
{ "P2",    "P121" },
{ "P21",   "P1211" },
{ "C2",    "C121" },
{ "Pm",    "P1m1" },
{ "Pc",    "P1c1" },
{ "Cm",    "C1m1" },
{ "Cc",    "C1c1" },
{ "P2/m",  "P12/m1" },
{ "P21/m", "P121/m1" },
{ "C2/m",  "C12/m1" },
{ "P2/c",  "P12/c1" },
{ "P21/c", "P121/c1" },
{ "C2/c",  "C12/c1" },
{ NULL, NULL }
};


typedef struct {
  int         SgNumber;
  const char  *Qualif;
  const char  *HM;
  const char  *Hall;
}
T_Main_HM_Dict;

static const
T_Main_HM_Dict Main_HM_Dict[] = {
{   1, NULL,   "P 1",        " P 1\0" },
{   2, NULL,   "P -1",       "-P 1\0" },
{   3, "b",    "P 1 2 1",    " P 2y\0" },
{   3, "c",    "P 1 1 2",    " P 2\0" },
{   3, "a",    "P 2 1 1",    " P 2x\0" },
{   4, "b",    "P 1 21 1",   " P 2yb\0" },
{   4, "c",    "P 1 1 21",   " P 2c\0" },
{   4, "a",    "P 21 1 1",   " P 2xa\0" },
{   5, "b1",   "C 1 2 1",    " C 2y\0" },
{   5, "b2",   "A 1 2 1",    " A 2y\0" },
{   5, "b3",   "I 1 2 1",    " I 2y\0" },
{   5, "c1",   "A 1 1 2",    " A 2\0" },
{   5, "c2",   "B 1 1 2",    " B 2\0" },
{   5, "c3",   "I 1 1 2",    " I 2\0" },
{   5, "a1",   "B 2 1 1",    " B 2x\0" },
{   5, "a2",   "C 2 1 1",    " C 2x\0" },
{   5, "a3",   "I 2 1 1",    " I 2x\0" },
{   6, "b",    "P 1 m 1",    " P -2y\0" },
{   6, "c",    "P 1 1 m",    " P -2\0" },
{   6, "a",    "P m 1 1",    " P -2x\0" },
{   7, "b1",   "P 1 c 1",    " P -2yc\0" },
{   7, "b2",   "P 1 n 1",    " P -2yac\0" },
{   7, "b3",   "P 1 a 1",    " P -2ya\0" },
{   7, "c1",   "P 1 1 a",    " P -2a\0" },
{   7, "c2",   "P 1 1 n",    " P -2ab\0" },
{   7, "c3",   "P 1 1 b",    " P -2b\0" },
{   7, "a1",   "P b 1 1",    " P -2xb\0" },
{   7, "a2",   "P n 1 1",    " P -2xbc\0" },
{   7, "a3",   "P c 1 1",    " P -2xc\0" },
{   8, "b1",   "C 1 m 1",    " C -2y\0" },
{   8, "b2",   "A 1 m 1",    " A -2y\0" },
{   8, "b3",   "I 1 m 1",    " I -2y\0" },
{   8, "c1",   "A 1 1 m",    " A -2\0" },
{   8, "c2",   "B 1 1 m",    " B -2\0" },
{   8, "c3",   "I 1 1 m",    " I -2\0" },
{   8, "a1",   "B m 1 1",    " B -2x\0" },
{   8, "a2",   "C m 1 1",    " C -2x\0" },
{   8, "a3",   "I m 1 1",    " I -2x\0" },
{   9, "b1",   "C 1 c 1",    " C -2yc\0" },
{   9, "b2",   "A 1 n 1",    " A -2yab\0" },
{   9, "b3",   "I 1 a 1",    " I -2ya\0" },
{   9, "-b1",  "A 1 a 1",    " A -2ya\0" },
{   9, "-b2",  "C 1 n 1",    " C -2yac\0" },
{   9, "-b3",  "I 1 c 1",    " I -2yc\0" },
{   9, "c1",   "A 1 1 a",    " A -2a\0" },
{   9, "c2",   "B 1 1 n",    " B -2ab\0" },
{   9, "c3",   "I 1 1 b",    " I -2b\0" },
{   9, "-c1",  "B 1 1 b",    " B -2b\0" },
{   9, "-c2",  "A 1 1 n",    " A -2ab\0" },
{   9, "-c3",  "I 1 1 a",    " I -2a\0" },
{   9, "a1",   "B b 1 1",    " B -2xb\0" },
{   9, "a2",   "C n 1 1",    " C -2xac\0" },
{   9, "a3",   "I c 1 1",    " I -2xc\0" },
{   9, "-a1",  "C c 1 1",    " C -2xc\0" },
{   9, "-a2",  "B n 1 1",    " B -2xab\0" },
{   9, "-a3",  "I b 1 1",    " I -2xb\0" },
{  10, "b",    "P 1 2/m 1",  "-P 2y\0" },
{  10, "c",    "P 1 1 2/m",  "-P 2\0" },
{  10, "a",    "P 2/m 1 1",  "-P 2x\0" },
{  11, "b",    "P 1 21/m 1", "-P 2yb\0" },
{  11, "c",    "P 1 1 21/m", "-P 2c\0" },
{  11, "a",    "P 21/m 1 1", "-P 2xa\0" },
{  12, "b1",   "C 1 2/m 1",  "-C 2y\0" },
{  12, "b2",   "A 1 2/m 1",  "-A 2y\0" },
{  12, "b3",   "I 1 2/m 1",  "-I 2y\0" },
{  12, "c1",   "A 1 1 2/m",  "-A 2\0" },
{  12, "c2",   "B 1 1 2/m",  "-B 2\0" },
{  12, "c3",   "I 1 1 2/m",  "-I 2\0" },
{  12, "a1",   "B 2/m 1 1",  "-B 2x\0" },
{  12, "a2",   "C 2/m 1 1",  "-C 2x\0" },
{  12, "a3",   "I 2/m 1 1",  "-I 2x\0" },
{  13, "b1",   "P 1 2/c 1",  "-P 2yc\0" },
{  13, "b2",   "P 1 2/n 1",  "-P 2yac\0" },
{  13, "b3",   "P 1 2/a 1",  "-P 2ya\0" },
{  13, "c1",   "P 1 1 2/a",  "-P 2a\0" },
{  13, "c2",   "P 1 1 2/n",  "-P 2ab\0" },
{  13, "c3",   "P 1 1 2/b",  "-P 2b\0" },
{  13, "a1",   "P 2/b 1 1",  "-P 2xb\0" },
{  13, "a2",   "P 2/n 1 1",  "-P 2xbc\0" },
{  13, "a3",   "P 2/c 1 1",  "-P 2xc\0" },
{  14, "b1",   "P 1 21/c 1", "-P 2ybc\0" },
{  14, "b2",   "P 1 21/n 1", "-P 2yn\0" },
{  14, "b3",   "P 1 21/a 1", "-P 2yab\0" },
{  14, "c1",   "P 1 1 21/a", "-P 2ac\0" },
{  14, "c2",   "P 1 1 21/n", "-P 2n\0" },
{  14, "c3",   "P 1 1 21/b", "-P 2bc\0" },
{  14, "a1",   "P 21/b 1 1", "-P 2xab\0" },
{  14, "a2",   "P 21/n 1 1", "-P 2xn\0" },
{  14, "a3",   "P 21/c 1 1", "-P 2xac\0" },
{  15, "b1",   "C 1 2/c 1",  "-C 2yc\0" },
{  15, "b2",   "A 1 2/n 1",  "-A 2yab\0" },
{  15, "b3",   "I 1 2/a 1",  "-I 2ya\0" },
{  15, "-b1",  "A 1 2/a 1",  "-A 2ya\0" },
{  15, "-b2",  "C 1 2/n 1",  "-C 2yac\0" },
{  15, "-b3",  "I 1 2/c 1",  "-I 2yc\0" },
{  15, "c1",   "A 1 1 2/a",  "-A 2a\0" },
{  15, "c2",   "B 1 1 2/n",  "-B 2ab\0" },
{  15, "c3",   "I 1 1 2/b",  "-I 2b\0" },
{  15, "-c1",  "B 1 1 2/b",  "-B 2b\0" },
{  15, "-c2",  "A 1 1 2/n",  "-A 2ab\0" },
{  15, "-c3",  "I 1 1 2/a",  "-I 2a\0" },
{  15, "a1",   "B 2/b 1 1",  "-B 2xb\0" },
{  15, "a2",   "C 2/n 1 1",  "-C 2xac\0" },
{  15, "a3",   "I 2/c 1 1",  "-I 2xc\0" },
{  15, "-a1",  "C 2/c 1 1",  "-C 2xc\0" },
{  15, "-a2",  "B 2/n 1 1",  "-B 2xab\0" },
{  15, "-a3",  "I 2/b 1 1",  "-I 2xb\0" },
{  16, NULL,   "P 2 2 2",    " P 2 2\0" },
{  17, NULL,   "P 2 2 21",   " P 2c 2\0" },
{  17, "cab",  "P 21 2 2",   " P 2a 2a\0" },
{  17, "bca",  "P 2 21 2",   " P 2 2b\0" },
{  18, NULL,   "P 21 21 2",  " P 2 2ab\0" },
{  18, "cab",  "P 2 21 21",  " P 2bc 2\0" },
{  18, "bca",  "P 21 2 21",  " P 2ac 2ac\0" },
{  19, NULL,   "P 21 21 21", " P 2ac 2ab\0" },
{  20, NULL,   "C 2 2 21",   " C 2c 2\0" },
{  20, "cab",  "A 21 2 2",   " A 2a 2a\0" },
{  20, "bca",  "B 2 21 2",   " B 2 2b\0" },
{  21, NULL,   "C 2 2 2",    " C 2 2\0" },
{  21, "cab",  "A 2 2 2",    " A 2 2\0" },
{  21, "bca",  "B 2 2 2",    " B 2 2\0" },
{  22, NULL,   "F 2 2 2",    " F 2 2\0" },
{  23, NULL,   "I 2 2 2",    " I 2 2\0" },
{  24, NULL,   "I 21 21 21", " I 2b 2c\0" },
{  25, NULL,   "P m m 2",    " P 2 -2\0" },
{  25, "cab",  "P 2 m m",    " P -2 2\0" },
{  25, "bca",  "P m 2 m",    " P -2 -2\0" },
{  26, NULL,   "P m c 21",   " P 2c -2\0" },
{  26, "ba-c", "P c m 21",   " P 2c -2c\0" },
{  26, "cab",  "P 21 m a",   " P -2a 2a\0" },
{  26, "-cba", "P 21 a m",   " P -2 2a\0" },
{  26, "bca",  "P b 21 m",   " P -2 -2b\0" },
{  26, "a-cb", "P m 21 b",   " P -2b -2\0" },
{  27, NULL,   "P c c 2",    " P 2 -2c\0" },
{  27, "cab",  "P 2 a a",    " P -2a 2\0" },
{  27, "bca",  "P b 2 b",    " P -2b -2b\0" },
{  28, NULL,   "P m a 2",    " P 2 -2a\0" },
{  28, "ba-c", "P b m 2",    " P 2 -2b\0" },
{  28, "cab",  "P 2 m b",    " P -2b 2\0" },
{  28, "-cba", "P 2 c m",    " P -2c 2\0" },
{  28, "bca",  "P c 2 m",    " P -2c -2c\0" },
{  28, "a-cb", "P m 2 a",    " P -2a -2a\0" },
{  29, NULL,   "P c a 21",   " P 2c -2ac\0" },
{  29, "ba-c", "P b c 21",   " P 2c -2b\0" },
{  29, "cab",  "P 21 a b",   " P -2b 2a\0" },
{  29, "-cba", "P 21 c a",   " P -2ac 2a\0" },
{  29, "bca",  "P c 21 b",   " P -2bc -2c\0" },
{  29, "a-cb", "P b 21 a",   " P -2a -2ab\0" },
{  30, NULL,   "P n c 2",    " P 2 -2bc\0" },
{  30, "ba-c", "P c n 2",    " P 2 -2ac\0" },
{  30, "cab",  "P 2 n a",    " P -2ac 2\0" },
{  30, "-cba", "P 2 a n",    " P -2ab 2\0" },
{  30, "bca",  "P b 2 n",    " P -2ab -2ab\0" },
{  30, "a-cb", "P n 2 b",    " P -2bc -2bc\0" },
{  31, NULL,   "P m n 21",   " P 2ac -2\0" },
{  31, "ba-c", "P n m 21",   " P 2bc -2bc\0" },
{  31, "cab",  "P 21 m n",   " P -2ab 2ab\0" },
{  31, "-cba", "P 21 n m",   " P -2 2ac\0" },
{  31, "bca",  "P n 21 m",   " P -2 -2bc\0" },
{  31, "a-cb", "P m 21 n",   " P -2ab -2\0" },
{  32, NULL,   "P b a 2",    " P 2 -2ab\0" },
{  32, "cab",  "P 2 c b",    " P -2bc 2\0" },
{  32, "bca",  "P c 2 a",    " P -2ac -2ac\0" },
{  33, NULL,   "P n a 21",   " P 2c -2n\0" },
{  33, "ba-c", "P b n 21",   " P 2c -2ab\0" },
{  33, "cab",  "P 21 n b",   " P -2bc 2a\0" },
{  33, "-cba", "P 21 c n",   " P -2n 2a\0" },
{  33, "bca",  "P c 21 n",   " P -2n -2ac\0" },
{  33, "a-cb", "P n 21 a",   " P -2ac -2n\0" },
{  34, NULL,   "P n n 2",    " P 2 -2n\0" },
{  34, "cab",  "P 2 n n",    " P -2n 2\0" },
{  34, "bca",  "P n 2 n",    " P -2n -2n\0" },
{  35, NULL,   "C m m 2",    " C 2 -2\0" },
{  35, "cab",  "A 2 m m",    " A -2 2\0" },
{  35, "bca",  "B m 2 m",    " B -2 -2\0" },
{  36, NULL,   "C m c 21",   " C 2c -2\0" },
{  36, "ba-c", "C c m 21",   " C 2c -2c\0" },
{  36, "cab",  "A 21 m a",   " A -2a 2a\0" },
{  36, "-cba", "A 21 a m",   " A -2 2a\0" },
{  36, "bca",  "B b 21 m",   " B -2 -2b\0" },
{  36, "a-cb", "B m 21 b",   " B -2b -2\0" },
{  37, NULL,   "C c c 2",    " C 2 -2c\0" },
{  37, "cab",  "A 2 a a",    " A -2a 2\0" },
{  37, "bca",  "B b 2 b",    " B -2b -2b\0" },
{  38, NULL,   "A m m 2",    " A 2 -2\0" },
{  38, "ba-c", "B m m 2",    " B 2 -2\0" },
{  38, "cab",  "B 2 m m",    " B -2 2\0" },
{  38, "-cba", "C 2 m m",    " C -2 2\0" },
{  38, "bca",  "C m 2 m",    " C -2 -2\0" },
{  38, "a-cb", "A m 2 m",    " A -2 -2\0" },
{  39, NULL,   "A b m 2",    " A 2 -2b\0" },
{  39, "ba-c", "B m a 2",    " B 2 -2a\0" },
{  39, "cab",  "B 2 c m",    " B -2a 2\0" },
{  39, "-cba", "C 2 m b",    " C -2a 2\0" },
{  39, "bca",  "C m 2 a",    " C -2a -2a\0" },
{  39, "a-cb", "A c 2 m",    " A -2b -2b\0" },
{  40, NULL,   "A m a 2",    " A 2 -2a\0" },
{  40, "ba-c", "B b m 2",    " B 2 -2b\0" },
{  40, "cab",  "B 2 m b",    " B -2b 2\0" },
{  40, "-cba", "C 2 c m",    " C -2c 2\0" },
{  40, "bca",  "C c 2 m",    " C -2c -2c\0" },
{  40, "a-cb", "A m 2 a",    " A -2a -2a\0" },
{  41, NULL,   "A b a 2",    " A 2 -2ab\0" },
{  41, "ba-c", "B b a 2",    " B 2 -2ab\0" },
{  41, "cab",  "B 2 c b",    " B -2ab 2\0" },
{  41, "-cba", "C 2 c b",    " C -2ac 2\0" },
{  41, "bca",  "C c 2 a",    " C -2ac -2ac\0" },
{  41, "a-cb", "A c 2 a",    " A -2ab -2ab\0" },
{  42, NULL,   "F m m 2",    " F 2 -2\0" },
{  42, "cab",  "F 2 m m",    " F -2 2\0" },
{  42, "bca",  "F m 2 m",    " F -2 -2\0" },
{  43, NULL,   "F d d 2",    " F 2 -2d\0" },
{  43, "cab",  "F 2 d d",    " F -2d 2\0" },
{  43, "bca",  "F d 2 d",    " F -2d -2d\0" },
{  44, NULL,   "I m m 2",    " I 2 -2\0" },
{  44, "cab",  "I 2 m m",    " I -2 2\0" },
{  44, "bca",  "I m 2 m",    " I -2 -2\0" },
{  45, NULL,   "I b a 2",    " I 2 -2c\0" },
{  45, "cab",  "I 2 c b",    " I -2a 2\0" },
{  45, "bca",  "I c 2 a",    " I -2b -2b\0" },
{  46, NULL,   "I m a 2",    " I 2 -2a\0" },
{  46, "ba-c", "I b m 2",    " I 2 -2b\0" },
{  46, "cab",  "I 2 m b",    " I -2b 2\0" },
{  46, "-cba", "I 2 c m",    " I -2c 2\0" },
{  46, "bca",  "I c 2 m",    " I -2c -2c\0" },
{  46, "a-cb", "I m 2 a",    " I -2a -2a\0" },
{  47, NULL,   "P m m m",    "-P 2 2\0" },
{  48, NULL,   "P n n n",    " P 2 2 -1n\0-P 2ab 2bc\0" },
{  49, NULL,   "P c c m",    "-P 2 2c\0" },
{  49, "cab",  "P m a a",    "-P 2a 2\0" },
{  49, "bca",  "P b m b",    "-P 2b 2b\0" },
{  50, NULL,   "P b a n",    " P 2 2 -1ab\0-P 2ab 2b\0" },
{  50, "cab",  "P n c b",    " P 2 2 -1bc\0-P 2b 2bc\0" },
{  50, "bca",  "P c n a",    " P 2 2 -1ac\0-P 2a 2c\0" },
{  51, NULL,   "P m m a",    "-P 2a 2a\0" },
{  51, "ba-c", "P m m b",    "-P 2b 2\0" },
{  51, "cab",  "P b m m",    "-P 2 2b\0" },
{  51, "-cba", "P c m m",    "-P 2c 2c\0" },
{  51, "bca",  "P m c m",    "-P 2c 2\0" },
{  51, "a-cb", "P m a m",    "-P 2 2a\0" },
{  52, NULL,   "P n n a",    "-P 2a 2bc\0" },
{  52, "ba-c", "P n n b",    "-P 2b 2n\0" },
{  52, "cab",  "P b n n",    "-P 2n 2b\0" },
{  52, "-cba", "P c n n",    "-P 2ab 2c\0" },
{  52, "bca",  "P n c n",    "-P 2ab 2n\0" },
{  52, "a-cb", "P n a n",    "-P 2n 2bc\0" },
{  53, NULL,   "P m n a",    "-P 2ac 2\0" },
{  53, "ba-c", "P n m b",    "-P 2bc 2bc\0" },
{  53, "cab",  "P b m n",    "-P 2ab 2ab\0" },
{  53, "-cba", "P c n m",    "-P 2 2ac\0" },
{  53, "bca",  "P n c m",    "-P 2 2bc\0" },
{  53, "a-cb", "P m a n",    "-P 2ab 2\0" },
{  54, NULL,   "P c c a",    "-P 2a 2ac\0" },
{  54, "ba-c", "P c c b",    "-P 2b 2c\0" },
{  54, "cab",  "P b a a",    "-P 2a 2b\0" },
{  54, "-cba", "P c a a",    "-P 2ac 2c\0" },
{  54, "bca",  "P b c b",    "-P 2bc 2b\0" },
{  54, "a-cb", "P b a b",    "-P 2b 2ab\0" },
{  55, NULL,   "P b a m",    "-P 2 2ab\0" },
{  55, "cab",  "P m c b",    "-P 2bc 2\0" },
{  55, "bca",  "P c m a",    "-P 2ac 2ac\0" },
{  56, NULL,   "P c c n",    "-P 2ab 2ac\0" },
{  56, "cab",  "P n a a",    "-P 2ac 2bc\0" },
{  56, "bca",  "P b n b",    "-P 2bc 2ab\0" },
{  57, NULL,   "P b c m",    "-P 2c 2b\0" },
{  57, "ba-c", "P c a m",    "-P 2c 2ac\0" },
{  57, "cab",  "P m c a",    "-P 2ac 2a\0" },
{  57, "-cba", "P m a b",    "-P 2b 2a\0" },
{  57, "bca",  "P b m a",    "-P 2a 2ab\0" },
{  57, "a-cb", "P c m b",    "-P 2bc 2c\0" },
{  58, NULL,   "P n n m",    "-P 2 2n\0" },
{  58, "cab",  "P m n n",    "-P 2n 2\0" },
{  58, "bca",  "P n m n",    "-P 2n 2n\0" },
{  59, NULL,   "P m m n",    " P 2 2ab -1ab\0-P 2ab 2a\0" },
{  59, "cab",  "P n m m",    " P 2bc 2 -1bc\0-P 2c 2bc\0" },
{  59, "bca",  "P m n m",    " P 2ac 2ac -1ac\0-P 2c 2a\0" },
{  60, NULL,   "P b c n",    "-P 2n 2ab\0" },
{  60, "ba-c", "P c a n",    "-P 2n 2c\0" },
{  60, "cab",  "P n c a",    "-P 2a 2n\0" },
{  60, "-cba", "P n a b",    "-P 2bc 2n\0" },
{  60, "bca",  "P b n a",    "-P 2ac 2b\0" },
{  60, "a-cb", "P c n b",    "-P 2b 2ac\0" },
{  61, NULL,   "P b c a",    "-P 2ac 2ab\0" },
{  61, "ba-c", "P c a b",    "-P 2bc 2ac\0" },
{  62, NULL,   "P n m a",    "-P 2ac 2n\0" },
{  62, "ba-c", "P m n b",    "-P 2bc 2a\0" },
{  62, "cab",  "P b n m",    "-P 2c 2ab\0" },
{  62, "-cba", "P c m n",    "-P 2n 2ac\0" },
{  62, "bca",  "P m c n",    "-P 2n 2a\0" },
{  62, "a-cb", "P n a m",    "-P 2c 2n\0" },
{  63, NULL,   "C m c m",    "-C 2c 2\0" },
{  63, "ba-c", "C c m m",    "-C 2c 2c\0" },
{  63, "cab",  "A m m a",    "-A 2a 2a\0" },
{  63, "-cba", "A m a m",    "-A 2 2a\0" },
{  63, "bca",  "B b m m",    "-B 2 2b\0" },
{  63, "a-cb", "B m m b",    "-B 2b 2\0" },
{  64, NULL,   "C m c a",    "-C 2ac 2\0" },
{  64, "ba-c", "C c m b",    "-C 2ac 2ac\0" },
{  64, "cab",  "A b m a",    "-A 2ab 2ab\0" },
{  64, "-cba", "A c a m",    "-A 2 2ab\0" },
{  64, "bca",  "B b c m",    "-B 2 2ab\0" },
{  64, "a-cb", "B m a b",    "-B 2ab 2\0" },
{  65, NULL,   "C m m m",    "-C 2 2\0" },
{  65, "cab",  "A m m m",    "-A 2 2\0" },
{  65, "bca",  "B m m m",    "-B 2 2\0" },
{  66, NULL,   "C c c m",    "-C 2 2c\0" },
{  66, "cab",  "A m a a",    "-A 2a 2\0" },
{  66, "bca",  "B b m b",    "-B 2b 2b\0" },
{  67, NULL,   "C m m a",    "-C 2a 2\0" },
{  67, "ba-c", "C m m b",    "-C 2a 2a\0" },
{  67, "cab",  "A b m m",    "-A 2b 2b\0" },
{  67, "-cba", "A c m m",    "-A 2 2b\0" },
{  67, "bca",  "B m c m",    "-B 2 2a\0" },
{  67, "a-cb", "B m a m",    "-B 2a 2\0" },
{  68, NULL,   "C c c a",    " C 2 2 -1ac\0-C 2a 2ac\0" },
{  68, "ba-c", "C c c b",    " C 2 2 -1ac\0-C 2a 2c\0" },
{  68, "cab",  "A b a a",    " A 2 2 -1ab\0-A 2a 2b\0" },
{  68, "-cba", "A c a a",    " A 2 2 -1ab\0-A 2ab 2b\0" },
{  68, "bca",  "B b c b",    " B 2 2 -1ab\0-B 2ab 2b\0" },
{  68, "a-cb", "B b a b",    " B 2 2 -1ab\0-B 2b 2ab\0" },
{  69, NULL,   "F m m m",    "-F 2 2\0" },
{  70, NULL,   "F d d d",    " F 2 2 -1d\0-F 2uv 2vw\0" },
{  71, NULL,   "I m m m",    "-I 2 2\0" },
{  72, NULL,   "I b a m",    "-I 2 2c\0" },
{  72, "cab",  "I m c b",    "-I 2a 2\0" },
{  72, "bca",  "I c m a",    "-I 2b 2b\0" },
{  73, NULL,   "I b c a",    "-I 2b 2c\0" },
{  73, "ba-c", "I c a b",    "-I 2a 2b\0" },
{  74, NULL,   "I m m a",    "-I 2b 2\0" },
{  74, "ba-c", "I m m b",    "-I 2a 2a\0" },
{  74, "cab",  "I b m m",    "-I 2c 2c\0" },
{  74, "-cba", "I c m m",    "-I 2 2b\0" },
{  74, "bca",  "I m c m",    "-I 2 2a\0" },
{  74, "a-cb", "I m a m",    "-I 2c 2\0" },
{  75, NULL,   "P 4",        " P 4\0" },
{  76, NULL,   "P 41",       " P 4w\0" },
{  77, NULL,   "P 42",       " P 4c\0" },
{  78, NULL,   "P 43",       " P 4cw\0" },
{  79, NULL,   "I 4",        " I 4\0" },
{  80, NULL,   "I 41",       " I 4bw\0" },
{  81, NULL,   "P -4",       " P -4\0" },
{  82, NULL,   "I -4",       " I -4\0" },
{  83, NULL,   "P 4/m",      "-P 4\0" },
{  84, NULL,   "P 42/m",     "-P 4c\0" },
{  85, NULL,   "P 4/n",      " P 4ab -1ab\0-P 4a\0" },
{  86, NULL,   "P 42/n",     " P 4n -1n\0-P 4bc\0" },
{  87, NULL,   "I 4/m",      "-I 4\0" },
{  88, NULL,   "I 41/a",     " I 4bw -1bw\0-I 4ad\0" },
{  89, NULL,   "P 4 2 2",    " P 4 2\0" },
{  90, NULL,   "P 4 21 2",   " P 4ab 2ab\0" },
{  91, NULL,   "P 41 2 2",   " P 4w 2c\0" },
{  92, NULL,   "P 41 21 2",  " P 4abw 2nw\0" },
{  93, NULL,   "P 42 2 2",   " P 4c 2\0" },
{  94, NULL,   "P 42 21 2",  " P 4n 2n\0" },
{  95, NULL,   "P 43 2 2",   " P 4cw 2c\0" },
{  96, NULL,   "P 43 21 2",  " P 4nw 2abw\0" },
{  97, NULL,   "I 4 2 2",    " I 4 2\0" },
{  98, NULL,   "I 41 2 2",   " I 4bw 2bw\0" },
{  99, NULL,   "P 4 m m",    " P 4 -2\0" },
{ 100, NULL,   "P 4 b m",    " P 4 -2ab\0" },
{ 101, NULL,   "P 42 c m",   " P 4c -2c\0" },
{ 102, NULL,   "P 42 n m",   " P 4n -2n\0" },
{ 103, NULL,   "P 4 c c",    " P 4 -2c\0" },
{ 104, NULL,   "P 4 n c",    " P 4 -2n\0" },
{ 105, NULL,   "P 42 m c",   " P 4c -2\0" },
{ 106, NULL,   "P 42 b c",   " P 4c -2ab\0" },
{ 107, NULL,   "I 4 m m",    " I 4 -2\0" },
{ 108, NULL,   "I 4 c m",    " I 4 -2c\0" },
{ 109, NULL,   "I 41 m d",   " I 4bw -2\0" },
{ 110, NULL,   "I 41 c d",   " I 4bw -2c\0" },
{ 111, NULL,   "P -4 2 m",   " P -4 2\0" },
{ 112, NULL,   "P -4 2 c",   " P -4 2c\0" },
{ 113, NULL,   "P -4 21 m",  " P -4 2ab\0" },
{ 114, NULL,   "P -4 21 c",  " P -4 2n\0" },
{ 115, NULL,   "P -4 m 2",   " P -4 -2\0" },
{ 116, NULL,   "P -4 c 2",   " P -4 -2c\0" },
{ 117, NULL,   "P -4 b 2",   " P -4 -2ab\0" },
{ 118, NULL,   "P -4 n 2",   " P -4 -2n\0" },
{ 119, NULL,   "I -4 m 2",   " I -4 -2\0" },
{ 120, NULL,   "I -4 c 2",   " I -4 -2c\0" },
{ 121, NULL,   "I -4 2 m",   " I -4 2\0" },
{ 122, NULL,   "I -4 2 d",   " I -4 2bw\0" },
{ 123, NULL,   "P 4/m m m",  "-P 4 2\0" },
{ 124, NULL,   "P 4/m c c",  "-P 4 2c\0" },
{ 125, NULL,   "P 4/n b m",  " P 4 2 -1ab\0-P 4a 2b\0" },
{ 126, NULL,   "P 4/n n c",  " P 4 2 -1n\0-P 4a 2bc\0" },
{ 127, NULL,   "P 4/m b m",  "-P 4 2ab\0" },
{ 128, NULL,   "P 4/m n c",  "-P 4 2n\0" },
{ 129, NULL,   "P 4/n m m",  " P 4ab 2ab -1ab\0-P 4a 2a\0" },
{ 130, NULL,   "P 4/n c c",  " P 4ab 2n -1ab\0-P 4a 2ac\0" },
{ 131, NULL,   "P 42/m m c", "-P 4c 2\0" },
{ 132, NULL,   "P 42/m c m", "-P 4c 2c\0" },
{ 133, NULL,   "P 42/n b c", " P 4n 2c -1n\0-P 4ac 2b\0" },
{ 134, NULL,   "P 42/n n m", " P 4n 2 -1n\0-P 4ac 2bc\0" },
{ 135, NULL,   "P 42/m b c", "-P 4c 2ab\0" },
{ 136, NULL,   "P 42/m n m", "-P 4n 2n\0" },
{ 137, NULL,   "P 42/n m c", " P 4n 2n -1n\0-P 4ac 2a\0" },
{ 138, NULL,   "P 42/n c m", " P 4n 2ab -1n\0-P 4ac 2ac\0" },
{ 139, NULL,   "I 4/m m m",  "-I 4 2\0" },
{ 140, NULL,   "I 4/m c m",  "-I 4 2c\0" },
{ 141, NULL,   "I 41/a m d", " I 4bw 2bw -1bw\0-I 4bd 2\0" },
{ 142, NULL,   "I 41/a c d", " I 4bw 2aw -1bw\0-I 4bd 2c\0" },
{ 143, NULL,   "P 3",        " P 3\0" },
{ 144, NULL,   "P 31",       " P 31\0" },
{ 145, NULL,   "P 32",       " P 32\0" },
{ 146, NULL,   "R 3",        " R 3\0 P 3*\0" },
{ 147, NULL,   "P -3",       "-P 3\0" },
{ 148, NULL,   "R -3",       "-R 3\0-P 3*\0" },
{ 149, NULL,   "P 3 1 2",    " P 3 2\0" },
{ 150, NULL,   "P 3 2 1",    " P 3 2\"\0" },
{ 151, NULL,   "P 31 1 2",   " P 31 2 (0 0 4)\0" },
{ 152, NULL,   "P 31 2 1",   " P 31 2\"\0" },
{ 153, NULL,   "P 32 1 2",   " P 32 2 (0 0 2)\0" },
{ 154, NULL,   "P 32 2 1",   " P 32 2\"\0" },
{ 155, NULL,   "R 3 2",      " R 3 2\"\0 P 3* 2\0" },
{ 156, NULL,   "P 3 m 1",    " P 3 -2\"\0" },
{ 157, NULL,   "P 3 1 m",    " P 3 -2\0" },
{ 158, NULL,   "P 3 c 1",    " P 3 -2\"c\0" },
{ 159, NULL,   "P 3 1 c",    " P 3 -2c\0" },
{ 160, NULL,   "R 3 m",      " R 3 -2\"\0 P 3* -2\0" },
{ 161, NULL,   "R 3 c",      " R 3 -2\"c\0 P 3* -2n\0" },
{ 162, NULL,   "P -3 1 m",   "-P 3 2\0" },
{ 163, NULL,   "P -3 1 c",   "-P 3 2c\0" },
{ 164, NULL,   "P -3 m 1",   "-P 3 2\"\0" },
{ 165, NULL,   "P -3 c 1",   "-P 3 2\"c\0" },
{ 166, NULL,   "R -3 m",     "-R 3 2\"\0-P 3* 2\0" },
{ 167, NULL,   "R -3 c",     "-R 3 2\"c\0-P 3* 2n\0" },
{ 168, NULL,   "P 6",        " P 6\0" },
{ 169, NULL,   "P 61",       " P 61\0" },
{ 170, NULL,   "P 65",       " P 65\0" },
{ 171, NULL,   "P 62",       " P 62\0" },
{ 172, NULL,   "P 64",       " P 64\0" },
{ 173, NULL,   "P 63",       " P 6c\0" },
{ 174, NULL,   "P -6",       " P -6\0" },
{ 175, NULL,   "P 6/m",      "-P 6\0" },
{ 176, NULL,   "P 63/m",     "-P 6c\0" },
{ 177, NULL,   "P 6 2 2",    " P 6 2\0" },
{ 178, NULL,   "P 61 2 2",   " P 61 2 (0 0 5)\0" },
{ 179, NULL,   "P 65 2 2",   " P 65 2 (0 0 1)\0" },
{ 180, NULL,   "P 62 2 2",   " P 62 2 (0 0 4)\0" },
{ 181, NULL,   "P 64 2 2",   " P 64 2 (0 0 2)\0" },
{ 182, NULL,   "P 63 2 2",   " P 6c 2c\0" },
{ 183, NULL,   "P 6 m m",    " P 6 -2\0" },
{ 184, NULL,   "P 6 c c",    " P 6 -2c\0" },
{ 185, NULL,   "P 63 c m",   " P 6c -2\0" },
{ 186, NULL,   "P 63 m c",   " P 6c -2c\0" },
{ 187, NULL,   "P -6 m 2",   " P -6 2\0" },
{ 188, NULL,   "P -6 c 2",   " P -6c 2\0" },
{ 189, NULL,   "P -6 2 m",   " P -6 -2\0" },
{ 190, NULL,   "P -6 2 c",   " P -6c -2c\0" },
{ 191, NULL,   "P 6/m m m",  "-P 6 2\0" },
{ 192, NULL,   "P 6/m c c",  "-P 6 2c\0" },
{ 193, NULL,   "P 63/m c m", "-P 6c 2\0" },
{ 194, NULL,   "P 63/m m c", "-P 6c 2c\0" },
{ 195, NULL,   "P 2 3",      " P 2 2 3\0" },
{ 196, NULL,   "F 2 3",      " F 2 2 3\0" },
{ 197, NULL,   "I 2 3",      " I 2 2 3\0" },
{ 198, NULL,   "P 21 3",     " P 2ac 2ab 3\0" },
{ 199, NULL,   "I 21 3",     " I 2b 2c 3\0" },
{ 200, NULL,   "P m -3",     "-P 2 2 3\0" },
{ 201, NULL,   "P n -3",     " P 2 2 3 -1n\0-P 2ab 2bc 3\0" },
{ 202, NULL,   "F m -3",     "-F 2 2 3\0" },
{ 203, NULL,   "F d -3",     " F 2 2 3 -1d\0-F 2uv 2vw 3\0" },
{ 204, NULL,   "I m -3",     "-I 2 2 3\0" },
{ 205, NULL,   "P a -3",     "-P 2ac 2ab 3\0" },
{ 206, NULL,   "I a -3",     "-I 2b 2c 3\0" },
{ 207, NULL,   "P 4 3 2",    " P 4 2 3\0" },
{ 208, NULL,   "P 42 3 2",   " P 4n 2 3\0" },
{ 209, NULL,   "F 4 3 2",    " F 4 2 3\0" },
{ 210, NULL,   "F 41 3 2",   " F 4d 2 3\0" },
{ 211, NULL,   "I 4 3 2",    " I 4 2 3\0" },
{ 212, NULL,   "P 43 3 2",   " P 4acd 2ab 3\0" },
{ 213, NULL,   "P 41 3 2",   " P 4bd 2ab 3\0" },
{ 214, NULL,   "I 41 3 2",   " I 4bd 2c 3\0" },
{ 215, NULL,   "P -4 3 m",   " P -4 2 3\0" },
{ 216, NULL,   "F -4 3 m",   " F -4 2 3\0" },
{ 217, NULL,   "I -4 3 m",   " I -4 2 3\0" },
{ 218, NULL,   "P -4 3 n",   " P -4n 2 3\0" },
{ 219, NULL,   "F -4 3 c",   " F -4a 2 3\0" },
{ 220, NULL,   "I -4 3 d",   " I -4bd 2c 3\0" },
{ 221, NULL,   "P m -3 m",   "-P 4 2 3\0" },
{ 222, NULL,   "P n -3 n",   " P 4 2 3 -1n\0-P 4a 2bc 3\0" },
{ 223, NULL,   "P m -3 n",   "-P 4n 2 3\0" },
{ 224, NULL,   "P n -3 m",   " P 4n 2 3 -1n\0-P 4bc 2bc 3\0" },
{ 225, NULL,   "F m -3 m",   "-F 4 2 3\0" },
{ 226, NULL,   "F m -3 c",   "-F 4a 2 3\0" },
{ 227, NULL,   "F d -3 m",   " F 4d 2 3 -1d\0-F 4vw 2vw 3\0" },
{ 228, NULL,   "F d -3 c",   " F 4d 2 3 -1ad\0-F 4ud 2vw 3\0" },
{ 229, NULL,   "I m -3 m",   "-I 4 2 3\0" },
{ 230, NULL,   "I a -3 d",   "-I 4bd 2c 3\0" },
{ 0, NULL, NULL, NULL }
};


/* remove whitespace and underscores, map to lower case
 */
static int PreProcessSymbol(const char *RawSymbol,
                            char *WorkSymbol, int BufSize)
{
  int         l;
  const char  *r;

  l = 0;
  for (r = RawSymbol; *r; r++) {
    if (! isspace(*r) && *r != '_') {
      if (l + 2 >= BufSize) return -1;
      WorkSymbol[l] = tolower(*r);
      l++;
    }
  }
  WorkSymbol[l] = '\0';

  return 0;
}


static int StripExtension(char *Symbol)
{
  int   Ext, l;
  char  *stop, *e;

  Ext = '\0';

      stop = strrchr(Symbol, ':');
  if (stop) {
    e = stop + 1;
    if (e[0] && ! e[1]) {
      Ext = e[0];
    }
    else if (e[0] == 'o' && (e[1] == '1' || e[1] == '2') && ! e[2]) {
      Ext = e[1];
    }
  }
  else {
    l = strlen(Symbol);
    /* check if last character is S, Z, R or H */
    if (l > 0) {
      e = &Symbol[l - 1];
      switch (e[0]) {
        case 's':
        case 'z':
        case 'r':
        case 'h':
          Ext = e[0];
          stop = e;
          break;
      }
      /* check if last two characters are O1 or O2 */
      if (! Ext && l > 1) {
        e = &Symbol[l - 2];
        if (e[0] == 'o' && (e[1] == '1' || e[1] == '2')) {
          Ext = e[1];
          stop = e;
        }
      }
    }
  }

  switch (Ext) {
    case '1': break;
    case 's': Ext = '1'; break;
    case '2': break;
    case 'z': Ext = '2'; break;
    case 'r': Ext = 'R'; break;
    case 'h': Ext = 'H'; break;
    default:
      Ext = '\0';
  }

  if (Ext) *stop = '\0';

  return Ext;
}


/* remove parentheses, e.g. "P2(1)2(1)2(1)" -> "P212121"
 */
static void RemoveParentheses(char *Symbol)
{
  int   ir, is;
  char  pat[5], *m;

  const int   RotNumbers[] = { 2, 3, 4, 6 };
  const char  RotSymbols[] = "2346";
  const char  ScrSymbols[] = "012345";

  strcpy(pat, "r(s)");

  range1(ir, 4) {
    pat[0] = RotSymbols[ir];
    for (is = 1; is < RotNumbers[ir]; is++) {
      pat[2] = ScrSymbols[is];
      for (;;) {
        m = strstr(Symbol, pat);
        if (! m) break;
        *m++ = pat[0];
        *m++ = pat[2];
        for (;;) {
          m[0] = m[2];
          if (! m[0]) break;
          m++;
        }
      }
    }
  }
}


static void RemoveSpaces(const char *source, char *target)
{
  for (;; source++) {
    if (*source != ' ') {
      *target = *source;
      if (! *target) break;
      target++;
    }
  }
}


static const char *SgNumber_as_HM_from_Main_HM_Dict(int SgNumber)
{
  const T_Main_HM_Dict  *Dict;

  for (Dict = Main_HM_Dict; Dict->SgNumber; Dict++)
    if (Dict->SgNumber == SgNumber)
      return Dict->HM;

  return NULL;
}


static int SgNumber_as_HM(int TableID, int SgNumber, char *Symbol)
{
  int         i;
  const char  *HM;

  if (SgNumber < 1 || SgNumber > 230) return 0;
  if (SgNumber < 3 || SgNumber > 15) {
        HM = SgNumber_as_HM_from_Main_HM_Dict(SgNumber);
    if (HM == NULL) return IE(-1);
    RemoveSpaces(HM, Symbol);
  }
  else {
    i = 0;
    if (TableID == 'I') i = 1;
    strcpy(Symbol, Monoclinic_SgNumber_as_HM_List[SgNumber][i]);
  }

  return 1;
}


static int CmpSchoenfliesSymbols(const char *FromTable, const char *WorkSymbol)
{
  int  i;

  for (i = 0;; i++) {
      if (    FromTable[i] != WorkSymbol[i]
          && (FromTable[i] != '^' || isalpha(WorkSymbol[i])
                                  || isdigit(WorkSymbol[i])))
        return -1;
    if (! FromTable[i]) break;
  }

  return 0;
}


static int Schoenflies_as_SgNumber(const char *Symbol)
{
  int SgNumber;

  range2(SgNumber, 1, 231)
    if (CmpSchoenfliesSymbols(Schoenflies_List[SgNumber], Symbol) == 0)
      return SgNumber;

  return 0;
}


static void ShortMonoHM_as_FullMonoHM(int TableID, char *WorkSymbol)
{
  int                         i;
  const T_Short_Mono_HM_Dict  *Dict;

  if (TableID == 'I') Dict = VolI_Short_Mono_HM_Dict;
  else                Dict = VolA_Short_Mono_HM_Dict;

  for (i = 0; Dict[i].shrt; i++) {
    if (strcmp(WorkSymbol, Dict[i].shrt) == 0) {
      strcpy(WorkSymbol, Dict[i].full);
      break;
    }
  }
}


static void Reset_HM_as_Hall(T_HM_as_Hall *HM_as_Hall)
{
  HM_as_Hall->SgNumber  = 0;
  HM_as_Hall->Schoenfl  = NULL;
  HM_as_Hall->Qualif    = NULL;
  HM_as_Hall->HM        = NULL;
  HM_as_Hall->Extension = '\0';
  HM_as_Hall->Hall      = NULL;
}


static int Main_HM_Lookup(int TableID, const char *WorkSymbol, int Extension,
                          T_HM_as_Hall *HM_as_Hall)
{
  const T_Main_HM_Dict  *Dict;
  char                  HM[32], *s;
  const char            *Hall, *h;
  int                   i;

  for (Dict = Main_HM_Dict; Dict->SgNumber; Dict++) {
    RemoveSpaces(Dict->HM, HM);
    if (strcmp(HM, WorkSymbol) == 0) break;
    s = strchr(HM, '-');
    if (s) { /* reverse "-N" to "N-", e.g. "P -1" -> "P 1-" */
      s[0] = s[1];
      s[1] = '-';
      if (strcmp(HM, WorkSymbol) == 0) break;
      if (   (Dict->SgNumber >= 200 && Dict->SgNumber <= 206)
          || (Dict->SgNumber >= 221 && Dict->SgNumber <= 230)) {
        for (i = 1;; i++) { /* remove '-', e.g. "P m -3 m" -> "P m 3 m" */
          s[i] = s[i + 1];
          if (! s[i]) break;
        }
        if (strcmp(HM, WorkSymbol) == 0) break;
      }
    }
  }

  Hall = NULL;

  if (Dict->SgNumber) {
         h = strchr(Dict->Hall, '\0') + 1;
    if (*h == '\0') {
      if (Extension == '\0') Hall = Dict->Hall;
    }
    else {
      if (HM[0] == 'R') {
        if (Extension == '\0') {
          if (TableID == 'I') Extension = 'R';
          else                Extension = 'H';
        }
        if      (Extension == 'H') Hall = Dict->Hall;
        else if (Extension == 'R') Hall = h;
      }
      else {
        if (Extension == '\0') {
          if (TableID == '\0') Extension = '2';
          else                 Extension = '1';
        }
        if      (Extension == '1') Hall = Dict->Hall;
        else if (Extension == '2') Hall = h;
      }
    }
  }

  if (Hall == NULL) return 0;

  if (HM_as_Hall) {
    HM_as_Hall->SgNumber  = Dict->SgNumber;
    HM_as_Hall->Schoenfl  = Schoenflies_List[Dict->SgNumber];
    HM_as_Hall->Qualif    = (Dict->Qualif ? Dict->Qualif : "");
    HM_as_Hall->HM        = Dict->HM;
    HM_as_Hall->Extension = Extension;
    HM_as_Hall->Hall      = Hall;
  }
  return Dict->SgNumber;
}


static int HallPassThrough(const char *Symbol, T_HM_as_Hall *HM_as_Hall)
{
  int         i;
  const char  *s;

  for (s = Symbol; *s; s++) if (! isspace(*s)) break;
  for (i = 0; i < 4; i++, s++) {
    if (tolower(*s) != "hall"[i]) return 0;
  }
  if (*s == ':') s++;
  else if (! isspace(*s)) return 0;
  for (; *s; s++) if (! isspace(*s)) break;
  if (HM_as_Hall) HM_as_Hall->Hall = s;
  return 1;
}


int SgSymbolLookup(int TableID, const char *Symbol, T_HM_as_Hall *HM_as_Hall)
{
  char  WorkSymbol[64], xtrac;
  int   Extension, SgNumber, n, status;

  if (HM_as_Hall) Reset_HM_as_Hall(HM_as_Hall);

  if      (TableID == 'I' || TableID == 'i' || TableID == '1')
    TableID = 'I';
  else if (TableID == 'A' || TableID == 'a')
    TableID = 'A';
  else {
    TableID = '\0';
    if (HallPassThrough(Symbol, HM_as_Hall) != 0)
      return 0; /* Attention: HM_as_Hall->Hall is a pointer to Symbol! */
  }

  if (PreProcessSymbol(Symbol, WorkSymbol, sizeof WorkSymbol) != 0) return 0;
  Extension = StripExtension(WorkSymbol);
  WorkSymbol[0] = toupper(WorkSymbol[0]);
  RemoveParentheses(WorkSymbol);

      n = sscanf(WorkSymbol, "%d%c", &SgNumber, &xtrac);
  if (n == 1) {
        status = SgNumber_as_HM(TableID, SgNumber, WorkSymbol);
    if (status <= 0) return status;
  }
  else {
        SgNumber = Schoenflies_as_SgNumber(WorkSymbol);
    if (SgNumber != 0) {
          status = SgNumber_as_HM(TableID, SgNumber, WorkSymbol);
      if (status <= 0) return IE(-1);
    }
  }

  ShortMonoHM_as_FullMonoHM(TableID, WorkSymbol);

  return Main_HM_Lookup(TableID, WorkSymbol, Extension, HM_as_Hall);
}


int MatchTabulatedSettings(const T_SgOps *SgOps, T_HM_as_Hall *HM_as_Hall)
{
  int                   SymCType, ixPG_SgOps;
  int                   iExt, jExt;
  const char            *Hall;
  const T_Main_HM_Dict  *Dict;
  T_SgOps               TidyOps[1], TabOps[1];

  const int Extensions[2][3] = {{'\0', '1', '2'}, {'\0', 'H', 'R'}};

  if (HM_as_Hall) Reset_HM_as_Hall(HM_as_Hall);

      SymCType = GetSymCType(SgOps->nLTr, SgOps->LTr);
  if (SymCType != '\0' && SymCType != 'Q')
  {
        ixPG_SgOps = ixPG(GetPG(SgOps));
    if (ixPG_SgOps == MGC_Unknown) return -1;

    SgOpsCpy(TidyOps, SgOps);
    if (TidySgOps(TidyOps) != 0) return IE(-1);

    for (Dict = Main_HM_Dict; Dict->SgNumber; Dict++) {
      if (ixPG(RefSetMGC[Dict->SgNumber]) != ixPG_SgOps) continue;
      for (Hall = Dict->Hall; *Hall; Hall = strchr(Hall, '\0') + 1) {
        if (Hall[1] != SymCType) continue;
        ResetSgOps(TabOps);
        if (ParseHallSymbol(Hall, TabOps, PHSymOptPedantic) < 0) return IE(-1);
        if (TidySgOps(TabOps) != 0) return IE(-1);
        if (SgOpsCmp(TidyOps, TabOps) == 0) {
          if (HM_as_Hall) {
            iExt = 0;
            if (143 <= Dict->SgNumber && Dict->SgNumber < 168)
              iExt = 1;
            jExt = 0;
            if (Hall != Dict->Hall)
              jExt = 2;
            else if (*(strchr(Dict->Hall, '\0') + 1) != '\0')
              jExt = 1;
            HM_as_Hall->SgNumber  = Dict->SgNumber;
            HM_as_Hall->Schoenfl  = Schoenflies_List[Dict->SgNumber];
            HM_as_Hall->Qualif    = (Dict->Qualif ? Dict->Qualif : "");
            HM_as_Hall->HM        = Dict->HM;
            HM_as_Hall->Extension = Extensions[iExt][jExt];
            HM_as_Hall->Hall      = Hall;
          }
          return Dict->SgNumber;
        }
      }
    }
  }

  return 0;
}
