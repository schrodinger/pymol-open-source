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

/* Module for internal C-level PyMOL tests...*/

#include"os_predef.h"
#include"os_std.h"
#include"Base.h"

#include"MemoryDebug.h"
#include"Feedback.h"
#include"TestPyMOL.h"

#include"ObjectCGO.h"
#include"VFont.h"
#include"ObjectGadget.h"
#include"P.h"

#include"ObjectMap.h"
#include"Executive.h"
#include"ButMode.h"
#include"Control.h"

#include"PyMOL.h"

static int TestPyMOL_00_00(PyMOLGlobals *G)
{
  ObjectMap *obj;
  ObjectMapDesc _md,*md;
  ObjectMapState *ms =NULL;

  int a;

  md=&_md;

  md->mode = cObjectMap_OrthoMinMaxGrid;

  for(a=0;a<3;a++) {
    md->Grid[a] = 0.1F;
    md->MinCorner[a] = 0.0F;
    md->MaxCorner[a] = a+1.0F;
  }
  md->init_mode = -2;
  
  obj = ObjectMapNew(G);
  if(obj) {
    ms = ObjectMapNewStateFromDesc(G,obj,md,0,true);    
    ms->Active=true;
  }
  if(obj) {
    ObjectSetName((CObject*)obj,"00_00");
    ExecutiveManageObject(G,(CObject*)obj,-1,false);
  }
  return (obj!=NULL);
}

#define STR_MAX 100

static char *get_st(const char array[][STR_MAX])
{
  char *result=NULL,*p;
  size_t c=0,l=0;
  while(array[c][0]) {
    l+=strlen(array[c]);
    c++;
  }
  result = Alloc(char,l+1);
  p=result;

  l=0;
  c=0;
  while(array[c][0]) {
    strcpy(result+l,array[c]);
    l+=strlen(array[c]);
    c++;
  }
  return result;
}

static const char pdb_01_01[][STR_MAX] = {
  "ATOM      1  N   ASP E   1       4.868 -17.809  25.188  1.00 34.37      E    N\n",
  "ATOM      2  CA  ASP E   1       3.984 -16.723  25.698  1.00 33.85      E    C\n",
  "ATOM      3  CB  ASP E   1       4.633 -16.020  26.888  1.00 35.91      E    C\n",
  "ATOM      4  CG  ASP E   1       6.016 -15.468  26.567  1.00 39.23      E    C\n",
  "ATOM      5  OD1 ASP E   1       6.340 -14.367  27.058  1.00 42.40      E    O\n",
  "ATOM      6  OD2 ASP E   1       6.787 -16.131  25.836  1.00 39.28      E    O\n",
  "ATOM      7  C   ASP E   1       3.789 -15.753  24.546  1.00 33.96      E    C\n",
  "ATOM      8  O   ASP E   1       4.456 -15.889  23.517  1.00 35.44      E    O\n",
  "ATOM      9  N   CYS E   2       2.908 -14.771  24.711  1.00 31.94      E    N\n",
  "ATOM     10  CA  CYS E   2       2.638 -13.825  23.633  1.00 29.67      E    C\n",
  "ATOM     11  CB  CYS E   2       1.198 -13.996  23.132  1.00 28.73      E    C\n",
  "ATOM     12  SG  CYS E   2       0.725 -15.681  22.638  1.00 25.85      E    S\n",
  "ATOM     13  C   CYS E   2       2.842 -12.366  24.012  1.00 30.31      E    C\n",
  "ATOM     14  O   CYS E   2       3.025 -12.039  25.188  1.00 30.74      E    O\n",
  "ATOM     15  N   ALA E   3       2.792 -11.504  22.996  1.00 30.15      E    N\n",
  "ATOM     16  CA  ALA E   3       2.923 -10.056  23.151  1.00 30.15      E    C\n",
  "ATOM     17  CB  ALA E   3       4.226  -9.566  22.552  1.00 30.82      E    C\n",
  "ATOM     18  C   ALA E   3       1.736  -9.418  22.436  1.00 29.90      E    C\n",
  "ATOM     19  O   ALA E   3       1.332  -9.867  21.362  1.00 28.16      E    O\n",
  "ATOM     20  N   TRP E   4       1.173  -8.377  23.038  1.00 31.77      E    N\n",
  "ATOM     21  CA  TRP E   4       0.007  -7.715  22.471  1.00 32.08      E    C\n",
  "ATOM     22  CB  TRP E   4      -1.233  -7.992  23.344  1.00 32.83      E    C\n",
  "ATOM     23  CG  TRP E   4      -1.500  -9.458  23.636  1.00 35.04      E    C\n",
  "ATOM     24  CD1 TRP E   4      -0.831 -10.264  24.528  1.00 35.36      E    C\n",
  "ATOM     25  CD2 TRP E   4      -2.507 -10.285  23.037  1.00 35.86      E    C\n",
  "ATOM     26  CE2 TRP E   4      -2.390 -11.576  23.610  1.00 34.80      E    C\n",
  "ATOM     27  CE3 TRP E   4      -3.500 -10.059  22.069  1.00 36.66      E    C\n",
  "ATOM     28  NE1 TRP E   4      -1.360 -11.535  24.514  1.00 33.13      E    N\n",
  "ATOM     29  CZ2 TRP E   4      -3.228 -12.638  23.249  1.00 36.81      E    C\n",
  "ATOM     30  CZ3 TRP E   4      -4.338 -11.120  21.706  1.00 37.35      E    C\n",
  "ATOM     31  CH2 TRP E   4      -4.194 -12.390  22.297  1.00 38.27      E    C\n",
  "ATOM     32  C   TRP E   4       0.231  -6.212  22.372  1.00 32.17      E    C\n",
  "ATOM     33  O   TRP E   4       0.752  -5.592  23.297  1.00 33.32      E    O\n",
  "ATOM     34  N   HIS E   5      -0.135  -5.634  21.235  1.00 32.18      E    N\n",
  "ATOM     35  CA  HIS E   5      -0.006  -4.205  21.043  1.00 32.91      E    C\n",
  "ATOM     36  CB  HIS E   5       0.791  -3.871  19.783  1.00 33.02      E    C\n",
  "ATOM     37  CG  HIS E   5       0.939  -2.396  19.549  1.00 32.86      E    C\n",
  "ATOM     38  CD2 HIS E   5       0.582  -1.619  18.499  1.00 31.48      E    C\n",
  "ATOM     39  ND1 HIS E   5       1.470  -1.542  20.495  1.00 29.87      E    N\n",
  "ATOM     40  CE1 HIS E   5       1.431  -0.303  20.038  1.00 29.68      E    C\n",
  "ATOM     41  NE2 HIS E   5       0.896  -0.322  18.831  1.00 29.87      E    N\n",
  "ATOM     42  C   HIS E   5      -1.408  -3.636  20.918  1.00 34.58      E    C\n",
  "ATOM     43  O   HIS E   5      -2.092  -3.870  19.914  1.00 35.87      E    O\n",
  "ATOM     44  N   LEU E   6      -1.838  -2.904  21.943  1.00 34.54      E    N\n",
  "ATOM     45  CA  LEU E   6      -3.165  -2.295  21.956  1.00 34.66      E    C\n",
  "ATOM     46  CB  LEU E   6      -3.266  -1.199  20.892  1.00 34.44      E    C\n",
  "ATOM     47  CG  LEU E   6      -2.302  -0.023  21.018  1.00 35.32      E    C\n",
  "ATOM     48  CD1 LEU E   6      -2.422   0.863  19.781  1.00 34.55      E    C\n",
  "ATOM     49  CD2 LEU E   6      -2.582   0.753  22.302  1.00 33.68      E    C\n",
  "ATOM     50  C   LEU E   6      -4.242  -3.339  21.698  1.00 36.23      E    C\n",
  "ATOM     51  O   LEU E   6      -5.181  -3.100  20.933  1.00 35.26      E    O\n",
  "ATOM     52  N   GLY E   7      -4.087  -4.506  22.315  1.00 37.53      E    N\n",
  "ATOM     53  CA  GLY E   7      -5.063  -5.564  22.138  1.00 38.95      E    C\n",
  "ATOM     54  C   GLY E   7      -4.781  -6.514  20.988  1.00 39.52      E    C\n",
  "ATOM     55  O   GLY E   7      -5.188  -7.669  21.049  1.00 42.08      E    O\n",
  "ATOM     56  N   GLU E   8      -4.092  -6.046  19.948  1.00 39.31      E    N\n",
  "ATOM     57  CA  GLU E   8      -3.771  -6.883  18.787  1.00 38.11      E    C\n",
  "ATOM     58  CB  GLU E   8      -3.361  -6.012  17.607  1.00 40.57      E    C\n",
  "ATOM     59  CG  GLU E   8      -4.472  -5.586  16.679  1.00 43.48      E    C\n",
  "ATOM     60  CD  GLU E   8      -3.920  -4.806  15.506  1.00 47.43      E    C\n",
  "ATOM     61  OE1 GLU E   8      -3.572  -5.421  14.467  1.00 44.53      E    O\n",
  "ATOM     62  OE2 GLU E   8      -3.800  -3.572  15.644  1.00 50.21      E    O\n",
  "ATOM     63  C   GLU E   8      -2.646  -7.877  19.066  1.00 35.96      E    C\n",
  "ATOM     64  O   GLU E   8      -1.643  -7.529  19.691  1.00 37.36      E    O\n",
  "ATOM     65  N   LEU E   9      -2.793  -9.104  18.579  1.00 32.21      E    N\n",
  "ATOM     66  CA  LEU E   9      -1.763 -10.108  18.791  1.00 27.94      E    C\n",
  "ATOM     67  CB  LEU E   9      -2.275 -11.515  18.478  1.00 28.00      E    C\n",
  "ATOM     68  CG  LEU E   9      -1.255 -12.643  18.693  1.00 27.79      E    C\n",
  "ATOM     69  CD1 LEU E   9      -0.819 -12.720  20.160  1.00 24.02      E    C\n",
  "ATOM     70  CD2 LEU E   9      -1.848 -13.963  18.233  1.00 26.33      E    C\n",
  "ATOM     71  C   LEU E   9      -0.569  -9.800  17.911  1.00 25.48      E    C\n",
  "ATOM     72  O   LEU E   9      -0.699  -9.660  16.692  1.00 25.09      E    O\n",
  "ATOM     73  N   VAL E  10       0.589  -9.697  18.547  1.00 22.25      E    N\n",
  "ATOM     74  CA  VAL E  10       1.835  -9.417  17.858  1.00 18.77      E    C\n",
  "ATOM     75  CB  VAL E  10       2.797  -8.578  18.759  1.00 17.44      E    C\n",
  "ATOM     76  CG1 VAL E  10       4.131  -8.351  18.068  1.00 17.34      E    C\n",
  "ATOM     77  CG2 VAL E  10       2.166  -7.248  19.110  1.00 17.12      E    C\n",
  "ATOM     78  C   VAL E  10       2.537 -10.717  17.473  1.00 19.09      E    C\n",
  "ATOM     79  O   VAL E  10       2.708 -11.019  16.296  1.00 18.92      E    O\n",
  "ATOM     80  N   TRP E  11       2.886 -11.525  18.465  1.00 19.41      E    N\n",
  "ATOM     81  CA  TRP E  11       3.622 -12.742  18.191  1.00 18.86      E    C\n",
  "ATOM     82  CB  TRP E  11       5.059 -12.338  17.887  1.00 14.44      E    C\n",
  "ATOM     83  CG  TRP E  11       5.840 -13.342  17.171  1.00 14.67      E    C\n",
  "ATOM     84  CD1 TRP E  11       6.700 -14.257  17.709  1.00 11.97      E    C\n",
  "ATOM     85  CD2 TRP E  11       5.886 -13.519  15.756  1.00 16.13      E    C\n",
  "ATOM     86  CE2 TRP E  11       6.804 -14.562  15.500  1.00 17.59      E    C\n",
  "ATOM     87  CE3 TRP E  11       5.250 -12.894  14.676  1.00 18.07      E    C\n",
  "ATOM     88  NE1 TRP E  11       7.284 -14.992  16.711  1.00 13.12      E    N\n",
  "ATOM     89  CZ2 TRP E  11       7.101 -14.995  14.196  1.00 16.12      E    C\n",
  "ATOM     90  CZ3 TRP E  11       5.548 -13.321  13.382  1.00 16.68      E    C\n",
  "ATOM     91  CH2 TRP E  11       6.466 -14.362  13.157  1.00 17.19      E    C\n",
  "ATOM     92  C   TRP E  11       3.637 -13.609  19.431  1.00 21.88      E    C\n",
  "ATOM     93  O   TRP E  11       3.449 -13.097  20.534  1.00 24.04      E    O\n",
  "ATOM     94  N   CYS E  12       3.878 -14.908  19.252  1.00 24.56      E    N\n",
  "ATOM     95  CA  CYS E  12       3.978 -15.845  20.367  1.00 27.16      E    C\n",
  "ATOM     96  CB  CYS E  12       2.709 -16.674  20.541  1.00 25.57      E    C\n",
  "ATOM     97  SG  CYS E  12       1.146 -15.762  20.653  1.00 29.62      E    S\n",
  "ATOM     98  C   CYS E  12       5.126 -16.810  20.107  1.00 30.17      E    C\n",
  "ATOM     99  O   CYS E  12       5.278 -17.322  18.998  1.00 29.17      E    O\n",
  "ATOM    100  N   THR E  13       5.959 -17.026  21.117  1.00 35.02      E    N\n",
  "ATOM    101  CA  THR E  13       7.053 -17.973  20.984  1.00 39.40      E    C\n",
  "ATOM    102  CB  THR E  13       8.289 -17.578  21.828  1.00 38.42      E    C\n",
  "ATOM    103  CG2 THR E  13       8.919 -16.310  21.286  1.00 37.78      E    C\n",
  "ATOM    104  OG1 THR E  13       7.908 -17.397  23.194  1.00 38.69      E    O\n",
  "ATOM    105  C   THR E  13       6.513 -19.322  21.459  1.00 42.05      E    C\n",
  "ATOM    106  O   THR E  13       5.962 -19.432  22.570  1.00 41.54      E    O\n",
  "ATOM    107  OXT THR E  13       6.606 -20.331  20.602  1.00 43.02      E    O\n",
  "END\n",
  ""};

static const char mol_01_02[][STR_MAX] = {
"MFCD02681585\n",
"  ChemPy            3D                             0\n",
"\n",
" 36 39  0  0  1  0  0  0  0  0999 V2000\n",
"   52.5122   32.1815   21.0164 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   53.1716   32.3766   20.0197 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   54.4517   31.7147   19.8169 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   54.7302   30.3758   20.0244 N   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   56.0228   30.0429   19.7805 N   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   56.5952   31.1894   19.3737 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   55.6544   32.2488   19.3610 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   58.0576   31.2220   18.9914 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   58.8985   30.2129   19.7561 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   58.1952   30.8995   17.4763 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   58.6099   32.6423   19.2418 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   53.8260   29.3227   20.4323 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   54.0916   28.8513   21.8755 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   53.2812   29.2854   22.9295 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   53.4667   28.8569   24.2265 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   52.5554   29.3273   25.3339 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   54.5278   27.9533   24.4780 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   55.3346   27.5330   23.4437 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   55.1252   27.9853   22.1551 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   52.8162   33.2369   18.9768 N   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   51.5459   33.8650   18.8063 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   50.7291   34.3547   19.8196 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   49.4781   34.9470   19.5308 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   48.6170   35.3669   20.6527 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   49.0272   35.8529   21.7164 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   47.3132   35.1061   20.3115 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   46.3612   35.5386   21.2706 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   49.0388   35.0753   18.1950 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   49.8568   34.5421   17.1857 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   51.0846   33.9564   17.4684 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   47.7773   35.7058   17.8640 N   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   46.9256   35.1276   16.8017 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   45.5301   35.6792   16.9429 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   45.5107   37.1034   16.9789 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   46.2587   37.5980   18.0839 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"   47.7345   37.1658   17.9831 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
"  1  2  2  0  0  0  0\n",
"  2  3  1  0  0  0  0\n",
"  2 20  1  0  0  0  0\n",
"  3  4  4  0  0  0  0\n",
"  3  7  4  0  0  0  0\n",
"  4  5  4  0  0  0  0\n",
"  4 12  1  0  0  0  0\n",
"  5  6  4  0  0  0  0\n",
"  6  7  4  0  0  0  0\n",
"  6  8  1  0  0  0  0\n",
"  8  9  1  0  0  0  0\n",
"  8 10  1  0  0  0  0\n",
"  8 11  1  0  0  0  0\n",
" 12 13  1  0  0  0  0\n",
" 13 14  4  0  0  0  0\n",
" 13 19  4  0  0  0  0\n",
" 14 15  4  0  0  0  0\n",
" 15 16  1  0  0  0  0\n",
" 15 17  4  0  0  0  0\n",
" 17 18  4  0  0  0  0\n",
" 18 19  4  0  0  0  0\n",
" 20 21  1  0  0  0  0\n",
" 21 22  4  0  0  0  0\n",
" 21 30  4  0  0  0  0\n",
" 22 23  4  0  0  0  0\n",
" 23 24  1  0  0  0  0\n",
" 23 28  4  0  0  0  0\n",
" 24 25  2  0  0  0  0\n",
" 24 26  1  0  0  0  0\n",
" 26 27  1  0  0  0  0\n",
" 28 29  4  0  0  0  0\n",
" 28 31  1  0  0  0  0\n",
" 29 30  4  0  0  0  0\n",
" 31 32  1  0  0  0  0\n",
" 31 36  1  0  0  0  0\n",
" 32 33  1  0  0  0  0\n",
" 33 34  1  0  0  0  0\n",
" 34 35  1  0  0  0  0\n",
" 35 36  1  0  0  0  0\n",
"M  END\n"};

int TestPyMOLRun(PyMOLGlobals *G,int group,int test)
{
  switch(group) {
  case 0: /* development tests */
    switch(test) {
    case 0:
      TestPyMOL_00_00(G); 
      break;
    case 1: 
      PBlock(G);
      VFontLoad(G,1,0,0,true); 
      PUnblock(G);
      break;
    case 2: 
      {
        CObject *obj = NULL;
        float pos[3] = {0.0,0.0,0.0};
        PBlock(G);
        obj = (CObject*)ObjectCGONewVFontTest(G,"hello",pos);
        PUnblock(G);
        if(obj) {
          ObjectSetName(obj,"hello");
          ExecutiveManageObject(G,obj,-1,false);
        }
      }
      break;
    case 3: 
      {
        CObject *obj = NULL;
        obj = (CObject*)ObjectGadgetTest(G);
        if(obj)  {
          ObjectSetName(obj,"gadget");
          ExecutiveManageObject(G,obj,-1,false);
        }
      }
      break;
    case 4:
      {
        /* try to match G3D */

        SettingSetGlobal_b(G,cSetting_ortho,1);
        SettingSet_3f(G->Setting,cSetting_light, 1.0F, -1.0F, -2.5F);;
      }
      break;
    }
    break;
  case 1: 
    /* set up for test usage as a simple viewer */
    
    PyMOL_SetDefaultMouse(G->PyMOL);

    switch(test) {
    case 1: 
      {
        char *st = get_st(pdb_01_01);

        PyMOL_CmdLoad(G->PyMOL, st, "string", "pdb", "test_01_01", 0, false, true, true, false, PYMOL_DEFAULT);
        ExecutiveSetRepVisib(G,"test_01_01",cRepCyl,1);
        ExecutiveSetRepVisib(G,"test_01_01",cRepLine,0);
        SettingSetGlobal_f(G,cSetting_sweep_speed,3.0F);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 2: 
      {
        char *st = get_st(pdb_01_01);
        PyMOL_CmdLoad(G->PyMOL,st, "string", "pdb", "test_01_02", 0, false, true, true, false, PYMOL_DEFAULT);
        ExecutiveSetRepVisib(G,"test_01_02",cRepLine,0);
        ExecutiveSetRepVisib(G,"test_01_02",cRepSurface,1);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 3: 
      {
        char *st = get_st(pdb_01_01);
        PyMOL_CmdLoad(G->PyMOL,st, "string", "pdb", "test_01_03", 0, false, true, true, false, PYMOL_DEFAULT);
        ExecutiveSetRepVisib(G,"test_01_03",cRepLine,0);
        ExecutiveSetRepVisib(G,"test_01_03",cRepCartoon,1);
        SettingSetGlobal_f(G,cSetting_sweep_speed,1.50F);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 4: 
      {
        char *st = get_st(pdb_01_01);
        PyMOL_CmdLoad(G->PyMOL,st, "string", "pdb", "test_01_04", 0, false, true, true, false, PYMOL_DEFAULT);
        ExecutiveSetRepVisib(G,"test_01_04",cRepLine,0);
        ExecutiveSetRepVisib(G,"test_01_04",cRepDot,1);
        SettingSetGlobal_f(G,cSetting_sweep_speed,1.50F);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 5: 
      {
        char *st = get_st(pdb_01_01);
        PyMOL_CmdLoad(G->PyMOL,st, "string", "pdb", "test_01_05", 0, false, true, true, false, PYMOL_DEFAULT);
        ExecutiveSetRepVisib(G,"test_01_05",cRepLine,0);
        ExecutiveSetRepVisib(G,"test_01_05",cRepSphere,1);
        SettingSetGlobal_f(G,cSetting_sweep_speed,4.50F);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 6: 
      {
        char *st = get_st(pdb_01_01);
        PyMOL_CmdLoad(G->PyMOL,st, "string", "pdb", "test_01_06", 0, false, true, true, false, PYMOL_DEFAULT);
        SettingSetGlobal_f(G,cSetting_sweep_speed,4.50F);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 7: 
      {
        char *st = get_st(mol_01_02);
        ExecutiveLoad(G,NULL,st,-1,cLoadTypeMOLStr,"test_01_07",0,-1,0,1,0,1,NULL);
        ExecutiveSetRepVisib(G,"test_01_07",cRepCyl,1);
        ExecutiveSetRepVisib(G,"test_01_07",cRepLine,0);
        SettingSetGlobal_b(G,cSetting_valence,1);
        SettingSetGlobal_f(G,cSetting_sweep_speed,0.25F);
        SettingSetGlobal_f(G,cSetting_sweep_angle,180.0F);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 8: 
      {
        char *st = get_st(mol_01_02);
        ExecutiveLoad(G,NULL,st,-1,cLoadTypeMOLStr,"test_01_08",0,-1,0,1,0,1,NULL);
        SettingSetGlobal_b(G,cSetting_valence,1);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    case 9: 
      {
        char *st = get_st(mol_01_02);
        ExecutiveLoad(G,NULL,st,-1,cLoadTypeMOLStr,"test_01_09",0,-1,0,1,0,1,NULL);
        ExecutiveSetRepVisib(G,"test_01_09",cRepMesh,1);
        ExecutiveSetRepVisib(G,"test_01_09",cRepLine,0);
        SettingSetGlobal_b(G,cSetting_valence,1);
        SettingSetGlobal_f(G,cSetting_sweep_speed,0.50F);
        SettingSetGlobal_f(G,cSetting_sweep_angle,90.0F);
        ControlRock(G,1);
        FreeP(st);
        break;
      }
      break;
    }
  }
  return 1;
}






