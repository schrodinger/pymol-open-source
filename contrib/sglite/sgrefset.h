/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#ifndef SGREFSET_H__
#define SGREFSET_H__


/* Hall symbols for the reference settings of the 230 crystallographic
   space groups.

   The reference settings chosen are identical to those listed in
   International Tables for Crystallography Vol. A. For the cases
   where more than one setting is given in the International Tables,
   the following choices have been made:
     - For monoclinic space groups: unique axis b and cell choice 1.
     - For space groups with two origin choices: origin choice 2.
     - Rhombohedral space groups: hexagonal axes.
*/
#ifndef SG_GLOBAL
extern
const char *RefSetHallSymbols[];
#else
const char *RefSetHallSymbols[] = {
            NULL,
  /* 001 */ " P 1",             /* P 1          */
  /* 002 */ "-P 1",             /* P -1         */
  /* 003 */ " P 2y",            /* P 1 2 1      */
  /* 004 */ " P 2yb",           /* P 1 21 1     */
  /* 005 */ " C 2y",            /* C 1 2 1      */
  /* 006 */ " P -2y",           /* P 1 m 1      */
  /* 007 */ " P -2yc",          /* P 1 c 1      */
  /* 008 */ " C -2y",           /* C 1 m 1      */
  /* 009 */ " C -2yc",          /* C 1 c 1      */
  /* 010 */ "-P 2y",            /* P 1 2/m 1    */
  /* 011 */ "-P 2yb",           /* P 1 21/m 1   */
  /* 012 */ "-C 2y",            /* C 1 2/m 1    */
  /* 013 */ "-P 2yc",           /* P 1 2/c 1    */
  /* 014 */ "-P 2ybc",          /* P 1 21/c 1   */
  /* 015 */ "-C 2yc",           /* C 1 2/c 1    */
  /* 016 */ " P 2 2",           /* P 2 2 2      */
  /* 017 */ " P 2c 2",          /* P 2 2 21     */
  /* 018 */ " P 2 2ab",         /* P 21 21 2    */
  /* 019 */ " P 2ac 2ab",       /* P 21 21 21   */
  /* 020 */ " C 2c 2",          /* C 2 2 21     */
  /* 021 */ " C 2 2",           /* C 2 2 2      */
  /* 022 */ " F 2 2",           /* F 2 2 2      */
  /* 023 */ " I 2 2",           /* I 2 2 2      */
  /* 024 */ " I 2b 2c",         /* I 21 21 21   */
  /* 025 */ " P 2 -2",          /* P m m 2      */
  /* 026 */ " P 2c -2",         /* P m c 21     */
  /* 027 */ " P 2 -2c",         /* P c c 2      */
  /* 028 */ " P 2 -2a",         /* P m a 2      */
  /* 029 */ " P 2c -2ac",       /* P c a 21     */
  /* 030 */ " P 2 -2bc",        /* P n c 2      */
  /* 031 */ " P 2ac -2",        /* P m n 21     */
  /* 032 */ " P 2 -2ab",        /* P b a 2      */
  /* 033 */ " P 2c -2n",        /* P n a 21     */
  /* 034 */ " P 2 -2n",         /* P n n 2      */
  /* 035 */ " C 2 -2",          /* C m m 2      */
  /* 036 */ " C 2c -2",         /* C m c 21     */
  /* 037 */ " C 2 -2c",         /* C c c 2      */
  /* 038 */ " A 2 -2",          /* A m m 2      */
  /* 039 */ " A 2 -2b",         /* A b m 2      */
  /* 040 */ " A 2 -2a",         /* A m a 2      */
  /* 041 */ " A 2 -2ab",        /* A b a 2      */
  /* 042 */ " F 2 -2",          /* F m m 2      */
  /* 043 */ " F 2 -2d",         /* F d d 2      */
  /* 044 */ " I 2 -2",          /* I m m 2      */
  /* 045 */ " I 2 -2c",         /* I b a 2      */
  /* 046 */ " I 2 -2a",         /* I m a 2      */
  /* 047 */ "-P 2 2",           /* P m m m      */
  /* 048 */ "-P 2ab 2bc",       /* P n n n:2    */
  /* 049 */ "-P 2 2c",          /* P c c m      */
  /* 050 */ "-P 2ab 2b",        /* P b a n:2    */
  /* 051 */ "-P 2a 2a",         /* P m m a      */
  /* 052 */ "-P 2a 2bc",        /* P n n a      */
  /* 053 */ "-P 2ac 2",         /* P m n a      */
  /* 054 */ "-P 2a 2ac",        /* P c c a      */
  /* 055 */ "-P 2 2ab",         /* P b a m      */
  /* 056 */ "-P 2ab 2ac",       /* P c c n      */
  /* 057 */ "-P 2c 2b",         /* P b c m      */
  /* 058 */ "-P 2 2n",          /* P n n m      */
  /* 059 */ "-P 2ab 2a",        /* P m m n:2    */
  /* 060 */ "-P 2n 2ab",        /* P b c n      */
  /* 061 */ "-P 2ac 2ab",       /* P b c a      */
  /* 062 */ "-P 2ac 2n",        /* P n m a      */
  /* 063 */ "-C 2c 2",          /* C m c m      */
  /* 064 */ "-C 2ac 2",         /* C m c a      */
  /* 065 */ "-C 2 2",           /* C m m m      */
  /* 066 */ "-C 2 2c",          /* C c c m      */
  /* 067 */ "-C 2a 2",          /* C m m a      */
  /* 068 */ "-C 2a 2ac",        /* C c c a:2    */
  /* 069 */ "-F 2 2",           /* F m m m      */
  /* 070 */ "-F 2uv 2vw",       /* F d d d:2    */
  /* 071 */ "-I 2 2",           /* I m m m      */
  /* 072 */ "-I 2 2c",          /* I b a m      */
  /* 073 */ "-I 2b 2c",         /* I b c a      */
  /* 074 */ "-I 2b 2",          /* I m m a      */
  /* 075 */ " P 4",             /* P 4          */
  /* 076 */ " P 4w",            /* P 41         */
  /* 077 */ " P 4c",            /* P 42         */
  /* 078 */ " P 4cw",           /* P 43         */
  /* 079 */ " I 4",             /* I 4          */
  /* 080 */ " I 4bw",           /* I 41         */
  /* 081 */ " P -4",            /* P -4         */
  /* 082 */ " I -4",            /* I -4         */
  /* 083 */ "-P 4",             /* P 4/m        */
  /* 084 */ "-P 4c",            /* P 42/m       */
  /* 085 */ "-P 4a",            /* P 4/n:2      */
  /* 086 */ "-P 4bc",           /* P 42/n:2     */
  /* 087 */ "-I 4",             /* I 4/m        */
  /* 088 */ "-I 4ad",           /* I 41/a:2     */
  /* 089 */ " P 4 2",           /* P 4 2 2      */
  /* 090 */ " P 4ab 2ab",       /* P 4 21 2     */
  /* 091 */ " P 4w 2c",         /* P 41 2 2     */
  /* 092 */ " P 4abw 2nw",      /* P 41 21 2    */
  /* 093 */ " P 4c 2",          /* P 42 2 2     */
  /* 094 */ " P 4n 2n",         /* P 42 21 2    */
  /* 095 */ " P 4cw 2c",        /* P 43 2 2     */
  /* 096 */ " P 4nw 2abw",      /* P 43 21 2    */
  /* 097 */ " I 4 2",           /* I 4 2 2      */
  /* 098 */ " I 4bw 2bw",       /* I 41 2 2     */
  /* 099 */ " P 4 -2",          /* P 4 m m      */
  /* 100 */ " P 4 -2ab",        /* P 4 b m      */
  /* 101 */ " P 4c -2c",        /* P 42 c m     */
  /* 102 */ " P 4n -2n",        /* P 42 n m     */
  /* 103 */ " P 4 -2c",         /* P 4 c c      */
  /* 104 */ " P 4 -2n",         /* P 4 n c      */
  /* 105 */ " P 4c -2",         /* P 42 m c     */
  /* 106 */ " P 4c -2ab",       /* P 42 b c     */
  /* 107 */ " I 4 -2",          /* I 4 m m      */
  /* 108 */ " I 4 -2c",         /* I 4 c m      */
  /* 109 */ " I 4bw -2",        /* I 41 m d     */
  /* 110 */ " I 4bw -2c",       /* I 41 c d     */
  /* 111 */ " P -4 2",          /* P -4 2 m     */
  /* 112 */ " P -4 2c",         /* P -4 2 c     */
  /* 113 */ " P -4 2ab",        /* P -4 21 m    */
  /* 114 */ " P -4 2n",         /* P -4 21 c    */
  /* 115 */ " P -4 -2",         /* P -4 m 2     */
  /* 116 */ " P -4 -2c",        /* P -4 c 2     */
  /* 117 */ " P -4 -2ab",       /* P -4 b 2     */
  /* 118 */ " P -4 -2n",        /* P -4 n 2     */
  /* 119 */ " I -4 -2",         /* I -4 m 2     */
  /* 120 */ " I -4 -2c",        /* I -4 c 2     */
  /* 121 */ " I -4 2",          /* I -4 2 m     */
  /* 122 */ " I -4 2bw",        /* I -4 2 d     */
  /* 123 */ "-P 4 2",           /* P 4/m m m    */
  /* 124 */ "-P 4 2c",          /* P 4/m c c    */
  /* 125 */ "-P 4a 2b",         /* P 4/n b m:2  */
  /* 126 */ "-P 4a 2bc",        /* P 4/n n c:2  */
  /* 127 */ "-P 4 2ab",         /* P 4/m b m    */
  /* 128 */ "-P 4 2n",          /* P 4/m n c    */
  /* 129 */ "-P 4a 2a",         /* P 4/n m m:2  */
  /* 130 */ "-P 4a 2ac",        /* P 4/n c c:2  */
  /* 131 */ "-P 4c 2",          /* P 42/m m c   */
  /* 132 */ "-P 4c 2c",         /* P 42/m c m   */
  /* 133 */ "-P 4ac 2b",        /* P 42/n b c:2 */
  /* 134 */ "-P 4ac 2bc",       /* P 42/n n m:2 */
  /* 135 */ "-P 4c 2ab",        /* P 42/m b c   */
  /* 136 */ "-P 4n 2n",         /* P 42/m n m   */
  /* 137 */ "-P 4ac 2a",        /* P 42/n m c:2 */
  /* 138 */ "-P 4ac 2ac",       /* P 42/n c m:2 */
  /* 139 */ "-I 4 2",           /* I 4/m m m    */
  /* 140 */ "-I 4 2c",          /* I 4/m c m    */
  /* 141 */ "-I 4bd 2",         /* I 41/a m d:2 */
  /* 142 */ "-I 4bd 2c",        /* I 41/a c d:2 */
  /* 143 */ " P 3",             /* P 3          */
  /* 144 */ " P 31",            /* P 31         */
  /* 145 */ " P 32",            /* P 32         */
  /* 146 */ " R 3",             /* R 3:h        */
  /* 147 */ "-P 3",             /* P -3         */
  /* 148 */ "-R 3",             /* R -3:h       */
  /* 149 */ " P 3 2",           /* P 3 1 2      */
  /* 150 */ " P 3 2\"",         /* P 3 2 1      */
  /* 151 */ " P 31 2 (0 0 4)",  /* P 31 1 2     */
  /* 152 */ " P 31 2\"",        /* P 31 2 1     */
  /* 153 */ " P 32 2 (0 0 2)",  /* P 32 1 2     */
  /* 154 */ " P 32 2\"",        /* P 32 2 1     */
  /* 155 */ " R 3 2\"",         /* R 3 2:h      */
  /* 156 */ " P 3 -2\"",        /* P 3 m 1      */
  /* 157 */ " P 3 -2",          /* P 3 1 m      */
  /* 158 */ " P 3 -2\"c",       /* P 3 c 1      */
  /* 159 */ " P 3 -2c",         /* P 3 1 c      */
  /* 160 */ " R 3 -2\"",        /* R 3 m:h      */
  /* 161 */ " R 3 -2\"c",       /* R 3 c:h      */
  /* 162 */ "-P 3 2",           /* P -3 1 m     */
  /* 163 */ "-P 3 2c",          /* P -3 1 c     */
  /* 164 */ "-P 3 2\"",         /* P -3 m 1     */
  /* 165 */ "-P 3 2\"c",        /* P -3 c 1     */
  /* 166 */ "-R 3 2\"",         /* R -3 m:h     */
  /* 167 */ "-R 3 2\"c",        /* R -3 c:h     */
  /* 168 */ " P 6",             /* P 6          */
  /* 169 */ " P 61",            /* P 61         */
  /* 170 */ " P 65",            /* P 65         */
  /* 171 */ " P 62",            /* P 62         */
  /* 172 */ " P 64",            /* P 64         */
  /* 173 */ " P 6c",            /* P 63         */
  /* 174 */ " P -6",            /* P -6         */
  /* 175 */ "-P 6",             /* P 6/m        */
  /* 176 */ "-P 6c",            /* P 63/m       */
  /* 177 */ " P 6 2",           /* P 6 2 2      */
  /* 178 */ " P 61 2 (0 0 5)",  /* P 61 2 2     */
  /* 179 */ " P 65 2 (0 0 1)",  /* P 65 2 2     */
  /* 180 */ " P 62 2 (0 0 4)",  /* P 62 2 2     */
  /* 181 */ " P 64 2 (0 0 2)",  /* P 64 2 2     */
  /* 182 */ " P 6c 2c",         /* P 63 2 2     */
  /* 183 */ " P 6 -2",          /* P 6 m m      */
  /* 184 */ " P 6 -2c",         /* P 6 c c      */
  /* 185 */ " P 6c -2",         /* P 63 c m     */
  /* 186 */ " P 6c -2c",        /* P 63 m c     */
  /* 187 */ " P -6 2",          /* P -6 m 2     */
  /* 188 */ " P -6c 2",         /* P -6 c 2     */
  /* 189 */ " P -6 -2",         /* P -6 2 m     */
  /* 190 */ " P -6c -2c",       /* P -6 2 c     */
  /* 191 */ "-P 6 2",           /* P 6/m m m    */
  /* 192 */ "-P 6 2c",          /* P 6/m c c    */
  /* 193 */ "-P 6c 2",          /* P 63/m c m   */
  /* 194 */ "-P 6c 2c",         /* P 63/m m c   */
  /* 195 */ " P 2 2 3",         /* P 2 3        */
  /* 196 */ " F 2 2 3",         /* F 2 3        */
  /* 197 */ " I 2 2 3",         /* I 2 3        */
  /* 198 */ " P 2ac 2ab 3",     /* P 21 3       */
  /* 199 */ " I 2b 2c 3",       /* I 21 3       */
  /* 200 */ "-P 2 2 3",         /* P m -3       */
  /* 201 */ "-P 2ab 2bc 3",     /* P n -3:2     */
  /* 202 */ "-F 2 2 3",         /* F m -3       */
  /* 203 */ "-F 2uv 2vw 3",     /* F d -3:2     */
  /* 204 */ "-I 2 2 3",         /* I m -3       */
  /* 205 */ "-P 2ac 2ab 3",     /* P a -3       */
  /* 206 */ "-I 2b 2c 3",       /* I a -3       */
  /* 207 */ " P 4 2 3",         /* P 4 3 2      */
  /* 208 */ " P 4n 2 3",        /* P 42 3 2     */
  /* 209 */ " F 4 2 3",         /* F 4 3 2      */
  /* 210 */ " F 4d 2 3",        /* F 41 3 2     */
  /* 211 */ " I 4 2 3",         /* I 4 3 2      */
  /* 212 */ " P 4acd 2ab 3",    /* P 43 3 2     */
  /* 213 */ " P 4bd 2ab 3",     /* P 41 3 2     */
  /* 214 */ " I 4bd 2c 3",      /* I 41 3 2     */
  /* 215 */ " P -4 2 3",        /* P -4 3 m     */
  /* 216 */ " F -4 2 3",        /* F -4 3 m     */
  /* 217 */ " I -4 2 3",        /* I -4 3 m     */
  /* 218 */ " P -4n 2 3",       /* P -4 3 n     */
  /* 219 */ " F -4a 2 3",       /* F -4 3 c     */
  /* 220 */ " I -4bd 2c 3",     /* I -4 3 d     */
  /* 221 */ "-P 4 2 3",         /* P m -3 m     */
  /* 222 */ "-P 4a 2bc 3",      /* P n -3 n:2   */
  /* 223 */ "-P 4n 2 3",        /* P m -3 n     */
  /* 224 */ "-P 4bc 2bc 3",     /* P n -3 m:2   */
  /* 225 */ "-F 4 2 3",         /* F m -3 m     */
  /* 226 */ "-F 4a 2 3",        /* F m -3 c     */
  /* 227 */ "-F 4vw 2vw 3",     /* F d -3 m:2   */
  /* 228 */ "-F 4ud 2vw 3",     /* F d -3 c:2   */
  /* 229 */ "-I 4 2 3",         /* I m -3 m     */
  /* 230 */ "-I 4bd 2c 3"       /* I a -3 d     */
};
#endif /* SG_GLOBAL */


/* Matrix Group Codes (Boisen & Gibbs, 1990, pp. 225-228) corresponding
   to the reference settings above.
 */
#ifndef SG_GLOBAL
extern
const int RefSetMGC[];
#else
const int RefSetMGC[] = {
  MGC_Unknown,
  MGC_1,
  MGC_1b,
  MGC_2,
  MGC_2,
  MGC_2,
  MGC_m,
  MGC_m,
  MGC_m,
  MGC_m,
  MGC_2_m,
  MGC_2_m,
  MGC_2_m,
  MGC_2_m,
  MGC_2_m,
  MGC_2_m,
  MGC_222,
  MGC_222,
  MGC_222,
  MGC_222,
  MGC_222,
  MGC_222,
  MGC_222,
  MGC_222,
  MGC_222,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mm2,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_mmm,
  MGC_4,
  MGC_4,
  MGC_4,
  MGC_4,
  MGC_4,
  MGC_4,
  MGC_4b,
  MGC_4b,
  MGC_4_m,
  MGC_4_m,
  MGC_4_m,
  MGC_4_m,
  MGC_4_m,
  MGC_4_m,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_422,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4mm,
  MGC_4b2m,
  MGC_4b2m,
  MGC_4b2m,
  MGC_4b2m,
  MGC_4bm2,
  MGC_4bm2,
  MGC_4bm2,
  MGC_4bm2,
  MGC_4bm2,
  MGC_4bm2,
  MGC_4b2m,
  MGC_4b2m,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_4_mmm,
  MGC_3,
  MGC_3,
  MGC_3,
  MGC_3,
  MGC_3b,
  MGC_3b,
  MGC_312,
  MGC_321,
  MGC_312,
  MGC_321,
  MGC_312,
  MGC_321,
  MGC_32,
  MGC_3m1,
  MGC_31m,
  MGC_3m1,
  MGC_31m,
  MGC_3m,
  MGC_3m,
  MGC_3b1m,
  MGC_3b1m,
  MGC_3bm1,
  MGC_3bm1,
  MGC_3bm,
  MGC_3bm,
  MGC_6,
  MGC_6,
  MGC_6,
  MGC_6,
  MGC_6,
  MGC_6,
  MGC_6b,
  MGC_6_m,
  MGC_6_m,
  MGC_622,
  MGC_622,
  MGC_622,
  MGC_622,
  MGC_622,
  MGC_622,
  MGC_6mm,
  MGC_6mm,
  MGC_6mm,
  MGC_6mm,
  MGC_6bm2,
  MGC_6bm2,
  MGC_6b2m,
  MGC_6b2m,
  MGC_6_mmm,
  MGC_6_mmm,
  MGC_6_mmm,
  MGC_6_mmm,
  MGC_23,
  MGC_23,
  MGC_23,
  MGC_23,
  MGC_23,
  MGC_m3b,
  MGC_m3b,
  MGC_m3b,
  MGC_m3b,
  MGC_m3b,
  MGC_m3b,
  MGC_m3b,
  MGC_432,
  MGC_432,
  MGC_432,
  MGC_432,
  MGC_432,
  MGC_432,
  MGC_432,
  MGC_432,
  MGC_4b3m,
  MGC_4b3m,
  MGC_4b3m,
  MGC_4b3m,
  MGC_4b3m,
  MGC_4b3m,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm,
  MGC_m3bm
};
#endif /* SG_GLOBAL */


/* List of 'additional generators' of the Euclidean and affine normalizers,
   corresponding to the reference settings above.
   Reference: Int. Tab. Vol. A Section 15.3
 */
typedef struct {
  const char  *K2L; /* operations which generate L from K  */
  const char  *L2N; /* operations which generate N from L  */
}
T_NormAddlG;

#ifndef SG_GLOBAL
extern
const T_NormAddlG RefSetNormAddlG[];
#else
const T_NormAddlG RefSetNormAddlG[] = {
            { NULL,   NULL },
  /*   1 */ { "-1",   NULL }, /* Note: the affine normalizers for   */
  /*   2 */ { NULL,   NULL }, /*   space groups 1 to 15 CANNOT be   */
  /*   3 */ { "-1",   NULL }, /*   represented by rotation symbols. */
  /*   4 */ { "-1",   NULL },
  /*   5 */ { "-1",   NULL },
  /*   6 */ { "-1",   NULL },
  /*   7 */ { "-1",   NULL },
  /*   8 */ { "-1",   NULL },
  /*   9 */ { "-1",   NULL },
  /*  10 */ { NULL,   NULL },
  /*  11 */ { NULL,   NULL },
  /*  12 */ { NULL,   NULL },
  /*  13 */ { NULL,   NULL },
  /*  14 */ { NULL,   NULL },
  /*  15 */ { NULL,   NULL },
  /*  16 */ { "-1",   "3* -2'" },  /* Space groups 16 to 74:      */
  /*  17 */ { "-1",   "-2'w" },    /*   the L2N column is for the */
  /*  18 */ { "-1",   "-2'" },     /*   affine normalizer only    */
  /*  19 */ { "-1",   "3* -2'd" },
  /*  20 */ { "-1",   "-2'w" },
  /*  21 */ { "-1",   "-2'" },
  /*  22 */ { "-1",   "3* -2'" },
  /*  23 */ { "-1",   "3* -2'" },
  /*  24 */ { "-1",   "3* -2'd" },
  /*  25 */ { "-1",   "-2'" },
  /*  26 */ { "-1",   NULL },
  /*  27 */ { "-1",   "-2'" },
  /*  28 */ { "-1",   NULL },
  /*  29 */ { "-1",   NULL },
  /*  30 */ { "-1",   NULL },
  /*  31 */ { "-1",   NULL },
  /*  32 */ { "-1",   "-2'" },
  /*  33 */ { "-1",   NULL },
  /*  34 */ { "-1",   "-2'" },
  /*  35 */ { "-1",   "-2'" },
  /*  36 */ { "-1",   NULL },
  /*  37 */ { "-1",   "-2'" },
  /*  38 */ { "-1",   NULL },
  /*  39 */ { "-1",   NULL },
  /*  40 */ { "-1",   NULL },
  /*  41 */ { "-1",   NULL },
  /*  42 */ { "-1",   "-2'" },
  /*  43 */ { "-1uv", "-2'" },
  /*  44 */ { "-1",   "-2'" },
  /*  45 */ { "-1",   "-2'" },
  /*  46 */ { "-1",   NULL },
  /*  47 */ { NULL,   "3* -2'" },
  /*  48 */ { NULL,   "3* -2'" },
  /*  49 */ { NULL,   "-2'" },
  /*  50 */ { NULL,   "-2'" },
  /*  51 */ { NULL,   NULL },
  /*  52 */ { NULL,   NULL },
  /*  53 */ { NULL,   NULL },
  /*  54 */ { NULL,   NULL },
  /*  55 */ { NULL,   "-2'" },
  /*  56 */ { NULL,   "-2'" },
  /*  57 */ { NULL,   NULL },
  /*  58 */ { NULL,   "-2'" },
  /*  59 */ { NULL,   "-2'" },
  /*  60 */ { NULL,   NULL },
  /*  61 */ { NULL,   "3*" },
  /*  62 */ { NULL,   NULL },
  /*  63 */ { NULL,   NULL },
  /*  64 */ { NULL,   NULL },
  /*  65 */ { NULL,   "-2'" },
  /*  66 */ { NULL,   "-2'" },
  /*  67 */ { NULL,   "-2'buv" },
  /*  68 */ { NULL,   "-2'buv" },
  /*  69 */ { NULL,   "3* -2'" },
  /*  70 */ { NULL,   "3* -2'" },
  /*  71 */ { NULL,   "3* -2'" },
  /*  72 */ { NULL,   "-2'" },
  /*  73 */ { NULL,   "3* -2'd" },
  /*  74 */ { NULL,   "-2'bd" },
  /*  75 */ { "-1",   "-2'" },
  /*  76 */ { NULL,   "2\"" },
  /*  77 */ { "-1",   "-2'" },
  /*  78 */ { NULL,   "2\"" },
  /*  79 */ { "-1",   "-2'" },
  /*  80 */ { "-1a",  "2\"" },
  /*  81 */ { "-1",   "-2'" },
  /*  82 */ { "-1",   "-2'" },
  /*  83 */ { NULL,   "-2'" },
  /*  84 */ { NULL,   "-2'" },
  /*  85 */ { NULL,   "-2'" },
  /*  86 */ { NULL,   "-2'" },
  /*  87 */ { NULL,   "-2'" },
  /*  88 */ { NULL,   "-2'd" }, /* space group 88 is the only space group */
  /*  89 */ { "-1",   NULL },   /* where the L2N column is different for  */
  /*  90 */ { "-1",   NULL },   /* the two origin choices                 */
  /*  91 */ { NULL,   NULL },
  /*  92 */ { NULL,   NULL },
  /*  93 */ { "-1",   NULL },
  /*  94 */ { "-1",   NULL },
  /*  95 */ { NULL,   NULL },
  /*  96 */ { NULL,   NULL },
  /*  97 */ { "-1",   NULL },
  /*  98 */ { "-1aw", NULL },
  /*  99 */ { "-1",   NULL },
  /* 100 */ { "-1",   NULL },
  /* 101 */ { "-1",   NULL },
  /* 102 */ { "-1",   NULL },
  /* 103 */ { "-1",   NULL },
  /* 104 */ { "-1",   NULL },
  /* 105 */ { "-1",   NULL },
  /* 106 */ { "-1",   NULL },
  /* 107 */ { "-1",   NULL },
  /* 108 */ { "-1",   NULL },
  /* 109 */ { "-1a",  NULL },
  /* 110 */ { "-1a",  NULL },
  /* 111 */ { "-1",   NULL },
  /* 112 */ { "-1",   NULL },
  /* 113 */ { "-1",   NULL },
  /* 114 */ { "-1",   NULL },
  /* 115 */ { "-1",   NULL },
  /* 116 */ { "-1",   NULL },
  /* 117 */ { "-1",   NULL },
  /* 118 */ { "-1",   NULL },
  /* 119 */ { "-1",   NULL },
  /* 120 */ { "-1",   NULL },
  /* 121 */ { "-1",   NULL },
  /* 122 */ { "-1aw", NULL },
  /* 123 */ { NULL,   NULL },
  /* 124 */ { NULL,   NULL },
  /* 125 */ { NULL,   NULL },
  /* 126 */ { NULL,   NULL },
  /* 127 */ { NULL,   NULL },
  /* 128 */ { NULL,   NULL },
  /* 129 */ { NULL,   NULL },
  /* 130 */ { NULL,   NULL },
  /* 131 */ { NULL,   NULL },
  /* 132 */ { NULL,   NULL },
  /* 133 */ { NULL,   NULL },
  /* 134 */ { NULL,   NULL },
  /* 135 */ { NULL,   NULL },
  /* 136 */ { NULL,   NULL },
  /* 137 */ { NULL,   NULL },
  /* 138 */ { NULL,   NULL },
  /* 139 */ { NULL,   NULL },
  /* 140 */ { NULL,   NULL },
  /* 141 */ { NULL,   NULL },
  /* 142 */ { NULL,   NULL },
  /* 143 */ { "-1",   "2 -2'" },
  /* 144 */ { NULL,   "2 2\"" },
  /* 145 */ { NULL,   "2 2\"" },
  /* 146 */ { "-1",   "-2\"" },
  /* 147 */ { NULL,   "2 -2'" },
  /* 148 */ { NULL,   "-2\"" },
  /* 149 */ { "-1",   "2" },
  /* 150 */ { "-1",   "2" },
  /* 151 */ { NULL,   "2" },
  /* 152 */ { NULL,   "2" },
  /* 153 */ { NULL,   "2" },
  /* 154 */ { NULL,   "2" },
  /* 155 */ { "-1",   NULL },
  /* 156 */ { "-1",   "2" },
  /* 157 */ { "-1",   "2" },
  /* 158 */ { "-1",   "2" },
  /* 159 */ { "-1",   "2" },
  /* 160 */ { "-1",   NULL },
  /* 161 */ { "-1",   NULL },
  /* 162 */ { NULL,   "2" },
  /* 163 */ { NULL,   "2" },
  /* 164 */ { NULL,   "2" },
  /* 165 */ { NULL,   "2" },
  /* 166 */ { NULL,   NULL },
  /* 167 */ { NULL,   NULL },
  /* 168 */ { "-1",   "-2'" },
  /* 169 */ { NULL,   "2\"" },
  /* 170 */ { NULL,   "2\"" },
  /* 171 */ { NULL,   "2\"" },
  /* 172 */ { NULL,   "2\"" },
  /* 173 */ { "-1",   "-2'" },
  /* 174 */ { "-1",   "-2'" },
  /* 175 */ { NULL,   "-2'" },
  /* 176 */ { NULL,   "-2'" },
  /* 177 */ { "-1",   NULL },
  /* 178 */ { NULL,   NULL },
  /* 179 */ { NULL,   NULL },
  /* 180 */ { NULL,   NULL },
  /* 181 */ { NULL,   NULL },
  /* 182 */ { "-1",   NULL },
  /* 183 */ { "-1",   NULL },
  /* 184 */ { "-1",   NULL },
  /* 185 */ { "-1",   NULL },
  /* 186 */ { "-1",   NULL },
  /* 187 */ { "-1",   NULL },
  /* 188 */ { "-1",   NULL },
  /* 189 */ { "-1",   NULL },
  /* 190 */ { "-1",   NULL },
  /* 191 */ { NULL,   NULL },
  /* 192 */ { NULL,   NULL },
  /* 193 */ { NULL,   NULL },
  /* 194 */ { NULL,   NULL },
  /* 195 */ { "-1",   "-2'" },
  /* 196 */ { "-1",   "-2'" },
  /* 197 */ { "-1",   "-2'" },
  /* 198 */ { "-1",   "-2'd" },
  /* 199 */ { "-1",   "-2'd" },
  /* 200 */ { NULL,   "-2'" },
  /* 201 */ { NULL,   "-2'" },
  /* 202 */ { NULL,   "-2'" },
  /* 203 */ { NULL,   "-2'" },
  /* 204 */ { NULL,   "-2'" },
  /* 205 */ { NULL,   NULL },
  /* 206 */ { NULL,   "-2'd" },
  /* 207 */ { "-1",   NULL },
  /* 208 */ { "-1",   NULL },
  /* 209 */ { "-1",   NULL },
  /* 210 */ { "-1d",  NULL },
  /* 211 */ { "-1",   NULL },
  /* 212 */ { NULL,   NULL },
  /* 213 */ { NULL,   NULL },
  /* 214 */ { "-1",   NULL },
  /* 215 */ { "-1",   NULL },
  /* 216 */ { "-1",   NULL },
  /* 217 */ { "-1",   NULL },
  /* 218 */ { "-1",   NULL },
  /* 219 */ { "-1",   NULL },
  /* 220 */ { "-1",   NULL },
  /* 221 */ { NULL,   NULL },
  /* 222 */ { NULL,   NULL },
  /* 223 */ { NULL,   NULL },
  /* 224 */ { NULL,   NULL },
  /* 225 */ { NULL,   NULL },
  /* 226 */ { NULL,   NULL },
  /* 227 */ { NULL,   NULL },
  /* 228 */ { NULL,   NULL },
  /* 229 */ { NULL,   NULL },
  /* 230 */ { NULL,   NULL }
};
#endif /* SG_GLOBAL */


#endif /* SGREFSET_H__ */
