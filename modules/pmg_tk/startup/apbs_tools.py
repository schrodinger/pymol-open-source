#!/usr/bin/env python

# APBS TOOLS Copyright Notice
# ============================
#
# The APBS TOOLS source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# APBS TOOLS is Copyright (C) 2005 by Michael G. Lerner
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Michael G. Lerner not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# MICHAEL G. LERNER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL MICHAEL G. LERNER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

"""
Recent changes:
Switched to new APBS input format
 - changed sfrm 1 to srfm smol
 - changed bcfl to use words instead of numbers
 - changed chgm to use words instead of numbers
 - added sdens parameter (default to APBS default of 10.0)
 - tweaked 'max grid points' to give an acceptable number .. we assume nlev == 4.  see comments in the source .. it's still possible to exceed max grid points by a little bit.
"""


from __future__ import division
from __future__ import generators

import os,math
import Tkinter
from Tkinter import *
import Pmw
import distutils.spawn # used for find_executable
import traceback

#
# Global config variables
#
APBS_BINARY_LOCATION = None
APBS_PSIZE_LOCATION = None

# Python backward-compatibility...
try:
    True
except:
    True = 1
try:
    False
except:
    False = 0
#
# Cheap hack for testing purposes
#
try:
    import pymol
    REAL_PYMOL = True
except ImportError:
    REAL_PYMOL = False
    class pymol:
        class cmd:
            def load(self,name,sel=''):
                pass
            def get_names(self):
                return ['mol1','mol2','map1','map2']
            def get_type(self,thing):
                if thing.startswith('mol'):
                    return 'object:molecule'
                else:
                    return 'object:map'
            def save(self,fname,selection):
                f = file(fname,'w')
                if fname.endswith('.pqr'):
                    f.write("""HEADER   this is a header, kids
REMARK    Format:                  x       y       z        charge  radius
ATOM      1    N MET     1    0.551 11.263 11.815 -0.300 1.8240
ATOM      2   CA MET     1    0.930 12.648 11.845 0.250 1.9643
ATOM      3    C MET     1    2.423 12.834 11.663 0.500 2.1046
ATOM      4    O MET     1    3.286 11.970 11.947 -0.500 1.6612
ATOM      5   CB MET     1    0.580 13.368 13.143 -0.120 1.9643
ATOM      6   CG MET     1    1.661 12.986 14.098 0.098 1.9643
ATOM      7   SD MET     1    1.441 13.562 15.776 -0.435 1.9924
ATOM      8   CE MET     1    2.528 12.285 16.478 0.037 1.9643
ATOM      9    N ILE     2    2.656 14.070 11.163 -0.500 1.8240
ATOM     10   CA ILE     2    3.942 14.570 10.821 0.140 1.9643
ATOM     11    C ILE     2    4.414 15.563 11.777 0.500 2.1046
ATOM     12    O ILE     2    3.666 16.500 12.106 -0.500 1.6612
ATOM     13   CB ILE     2    3.780 15.262 9.428 -0.060 1.9643
ATOM     14  CG1 ILE     2    3.485 14.180 8.372 -0.120 1.9643
ATOM     15  CG2 ILE     2    4.987 16.173 9.033 -0.180 1.9643
ATOM     16   CD ILE     2    2.969 14.795 7.082 -0.180 1.9643
ATOM     17    N SER     3    5.637 15.334 12.115 -0.500 1.8240
ATOM     18   CA SER     3    6.368 16.180 13.050 0.140 1.9643
ATOM     19    C SER     3    7.703 16.646 12.468 0.500 2.1046
ATOM     20    O SER     3    8.363 15.914 11.750 -0.500 1.6612
ATOM     21   CB SER     3    6.751 15.301 14.298 0.145 1.9643
ATOM     22   OG SER     3    5.544 14.967 14.987 -0.683 1.7510
ATOM     23    N LEU     4    8.112 17.862 12.816 -0.500 1.8240
ATOM     24   CA LEU     4    9.370 18.301 12.353 0.140 1.9643
ATOM     25    C LEU     4    10.283 18.416 13.583 0.500 2.1046
ATOM     26    O LEU     4    9.765 18.812 14.627 -0.500 1.6612
ATOM     27   CB LEU     4    9.230 19.754 11.875 -0.120 1.9643
ATOM     28   CG LEU     4    8.684 19.996 10.452 -0.060 1.9643
ATOM     29  CD1 LEU     4    7.569 19.115 9.957 -0.180 1.9643
ATOM     30  CD2 LEU     4    8.477 21.499 10.173 -0.180 1.9643
ATOM     31    N ILE     5    11.555 18.104 13.495 -0.500 1.8240
ATOM     32   CA ILE     5    12.465 18.230 14.624 0.140 1.9643
ATOM     33    C ILE     5    13.636 19.013 14.159 0.500 2.1046
ATOM     34    O ILE     5    14.222 18.681 13.116 -0.500 1.6612
ATOM     35   CB ILE     5    12.811 16.919 15.273 -0.060 1.9643
ATOM     36  CG1 ILE     5    13.886 17.156 16.315 -0.120 1.9643
ATOM     37  CG2 ILE     5    13.298 15.880 14.184 -0.180 1.9643
ATOM     38   CD ILE     5    13.981 15.939 17.251 -0.180 1.9643
ATOM     39    N ALA     6    14.012 20.093 14.895 -0.500 1.8240
ATOM     40   CA ALA     6    15.137 20.934 14.463 0.140 1.9643
ATOM     41    C ALA     6    15.851 21.640 15.619 0.500 2.1046
ATOM     42    O ALA     6    15.235 21.839 16.668 -0.500 1.6612
ATOM     43   CB ALA     6    14.449 22.032 13.594 -0.180 1.9643
ATOM     44    N ALA     7    17.119 21.996 15.415 -0.500 1.8240
ATOM     45   CA ALA     7    17.962 22.700 16.346 0.140 1.9643
ATOM     46    C ALA     7    18.099 24.081 15.714 0.500 2.1046
ATOM     47    O ALA     7    18.503 24.190 14.577 -0.500 1.6612
ATOM     48   CB ALA     7    19.279 22.022 16.626 -0.180 1.9643
ATOM     49    N LEU     8    17.738 25.159 16.435 -0.500 1.8240
ATOM     50   CA LEU     8    17.767 26.510 15.918 0.140 1.9643
ATOM     51    C LEU     8    18.542 27.441 16.789 0.500 2.1046
ATOM     52    O LEU     8    18.344 27.380 17.993 -0.500 1.6612
ATOM     53   CB LEU     8    16.306 27.133 16.150 -0.120 1.9643
ATOM     54   CG LEU     8    15.245 26.745 15.191 -0.060 1.9643
ATOM     55  CD1 LEU     8    14.714 25.361 15.532 -0.180 1.9643
ATOM     56  CD2 LEU     8    14.133 27.753 15.245 -0.180 1.9643
ATOM     57    N ALA     9    19.354 28.324 16.136 -0.500 1.8240
ATOM     58   CA ALA     9    20.109 29.322 16.859 0.140 1.9643
ATOM     59    C ALA     9    19.268 30.562 16.781 0.500 2.1046
ATOM     60    O ALA     9    18.052 30.512 16.395 -0.500 1.6612
ATOM     61   CB ALA     9    21.496 29.499 16.312 -0.180 1.9643
ATOM     62    N VAL    10    19.908 31.695 17.172 -0.500 1.8240
ATOM     63   CA VAL    10    19.182 32.944 17.121 0.140 1.9643
ATOM     64    C VAL    10    18.853 33.270 15.681 0.500 2.1046
ATOM     65    O VAL    10    19.589 32.940 14.704 -0.500 1.6612
ATOM     66   CB VAL    10    19.896 34.147 17.775 -0.060 1.9643
ATOM     67  CG1 VAL    10    20.668 33.742 19.022 -0.180 1.9643
ATOM     68  CG2 VAL    10    20.692 34.984 16.815 -0.180 1.9643
ATOM     69    N ASP    11    17.742 33.900 15.528 -0.500 1.8240
ATOM     70   CA ASP    11    17.367 34.231 14.199 0.140 1.9643
ATOM     71    C ASP    11    16.970 33.019 13.371 0.500 2.1046
ATOM     72    O ASP    11    16.959 33.131 12.206 -0.500 1.6612
ATOM     73   CB ASP    11    18.378 35.124 13.467 -0.220 1.9643
ATOM     74   CG ASP    11    18.356 36.564 14.035 0.700 2.1046
ATOM     75  OD1 ASP    11    17.333 37.176 14.327 -0.800 1.6612
ATOM     76  OD2 ASP    11    19.538 37.099 14.146 -0.800 1.6612
ATOM     77    N ARG    12    16.634 31.897 13.959 -0.500 1.8240
ATOM     78   CA ARG    12    16.193 30.711 13.247 0.140 1.9643
ATOM     79    C ARG    12    17.224 30.074 12.351 0.500 2.1046
ATOM     80    O ARG    12    16.841 29.226 11.518 -0.500 1.6612
ATOM     81   CB ARG    12    14.891 30.911 12.489 -0.120 1.9643
ATOM     82   CG ARG    12    13.770 31.437 13.345 -0.050 1.9643
ATOM     83   CD ARG    12    12.435 31.161 12.721 0.190 1.9643
ATOM     84   NE ARG    12    12.462 31.677 11.381 -0.700 1.8240
ATOM     85   CZ ARG    12    12.456 32.954 11.069 0.640 1.2628
ATOM     86  NH1 ARG    12    12.421 33.916 11.940 -0.800 1.8240
ATOM     87  NH2 ARG    12    12.494 33.278 9.801 -0.800 1.8240
ATOM     88    N VAL    13    18.498 30.450 12.518 -0.500 1.8240
ATOM     89   CA VAL    13    19.531 29.853 11.712 0.140 1.9643
ATOM     90    C VAL    13    19.672 28.345 12.064 0.500 2.1046
ATOM     91    O VAL    13    19.821 27.968 13.233 -0.500 1.6612
ATOM     92   CB VAL    13    20.863 30.533 11.990 -0.060 1.9643
ATOM     93  CG1 VAL    13    22.042 29.742 11.406 -0.180 1.9643
ATOM     94  CG2 VAL    13    20.797 31.998 11.467 -0.180 1.9643
ATOM     95    N ILE    14    19.629 27.482 11.056 -0.500 1.8240
ATOM     96   CA ILE    14    19.767 26.011 11.271 0.140 1.9643
ATOM     97    C ILE    14    20.995 25.409 10.552 0.500 2.1046
ATOM     98    O ILE    14    21.399 24.278 10.747 -0.500 1.6612
ATOM     99   CB ILE    14    18.506 25.235 10.843 -0.060 1.9643
ATOM    100  CG1 ILE    14    18.151 25.494 9.365 -0.120 1.9643
ATOM    101  CG2 ILE    14    17.307 25.571 11.723 -0.180 1.9643
ATOM    102   CD ILE    14    16.862 24.785 8.834 -0.180 1.9643
ATOM    103    N GLY    15    21.623 26.185 9.656 -0.500 1.8240
ATOM    104   CA GLY    15    22.773 25.644 8.908 0.080 1.9643
ATOM    105    C GLY    15    23.579 26.679 8.113 0.500 2.1046
ATOM    106    O GLY    15    23.216 27.824 8.023 -0.500 1.6612
ATOM    107    N MET    16    24.685 26.243 7.585 -0.500 1.8240
ATOM    108   CA MET    16    25.630 27.043 6.835 0.140 1.9643
ATOM    109    C MET    16    26.562 26.064 6.140 0.500 2.1046
ATOM    110    O MET    16    26.419 24.841 6.284 -0.500 1.6612
ATOM    111   CB MET    16    26.454 27.903 7.783 -0.120 1.9643
ATOM    112   CG MET    16    27.215 27.034 8.748 0.098 1.9643
ATOM    113   SD MET    16    28.614 27.916 9.432 -0.435 1.9924
ATOM    114   CE MET    16    29.524 27.955 7.874 0.037 1.9643
ATOM    115    N GLU    17    27.509 26.591 5.380 -0.500 1.8240
ATOM    116   CA GLU    17    28.375 25.714 4.685 0.140 1.9643
ATOM    117    C GLU    17    29.197 24.820 5.561 0.500 2.1046
ATOM    118    O GLU    17    29.293 23.610 5.320 -0.500 1.6612
ATOM    119   CB GLU    17    29.332 26.532 3.738 -0.120 1.9643
ATOM    120   CG GLU    17    28.646 27.132 2.520 -0.220 1.9643
ATOM    121   CD GLU    17    28.159 26.083 1.541 0.700 2.1046
ATOM    122  OE1 GLU    17    28.336 24.880 1.823 -0.800 1.6612
ATOM    123  OE2 GLU    17    27.599 26.466 0.492 -0.800 1.6612
ATOM    124    N ASN    18    29.817 25.419 6.557 -0.500 1.8240
ATOM    125   CA ASN    18    30.671 24.623 7.428 0.140 1.9643
ATOM    126    C ASN    18    29.968 24.072 8.644 0.500 2.1046
ATOM    127    O ASN    18    28.800 24.344 8.861 -0.500 1.6612
ATOM    128   CB ASN    18    31.891 25.429 7.870 -0.120 1.9643
ATOM    129   CG ASN    18    32.673 26.009 6.672 0.500 2.1046
ATOM    130  ND2 ASN    18    32.793 25.192 5.612 -0.760 1.8240
ATOM    131  OD1 ASN    18    33.131 27.181 6.693 -0.500 1.6612
ATOM    132    N ALA    19    30.721 23.301 9.419 -0.500 1.8240
ATOM    133   CA ALA    19    30.158 22.732 10.649 0.140 1.9643
ATOM    134    C ALA    19    29.808 23.892 11.608 0.500 2.1046
ATOM    135    O ALA    19    30.553 24.855 11.718 -0.500 1.6612
ATOM    136   CB ALA    19    31.208 21.851 11.326 -0.180 1.9643
ATOM    137    N MET    20    28.661 23.767 12.282 -0.500 1.8240
ATOM    138   CA MET    20    28.173 24.753 13.256 0.140 1.9643
ATOM    139    C MET    20    29.215 24.880 14.367 0.500 2.1046
ATOM    140    O MET    20    29.858 23.879 14.759 -0.500 1.6612
ATOM    141   CB MET    20    26.859 24.207 13.808 -0.120 1.9643
ATOM    142   CG MET    20    25.790 24.287 12.739 0.098 1.9643
ATOM    143   SD MET    20    25.744 25.961 12.033 -0.435 1.9924
ATOM    144   CE MET    20    24.199 26.587 12.680 0.037 1.9643
ATOM   1257    N ARG    21    -2.365 24.238 17.120 -0.500 1.8240
ATOM   1258   CA ARG    21    -3.722 24.429 17.637 0.040 1.9643
ATOM   1259    C ARG    21    -4.816 24.158 16.587 0.700 2.1046
ATOM   1260   O1 ARG    21    -4.718 24.783 15.486 -0.800 1.6612
ATOM   1261   CB ARG    21    -3.918 25.826 18.162 -0.120 1.9643
ATOM   1262   CG ARG    21    -3.019 26.176 19.341 -0.050 1.9643
ATOM   1263   CD ARG    21    -3.363 27.548 19.971 0.190 1.9643
ATOM   1264   NE ARG    21    -2.493 27.864 21.095 -0.700 1.8240
ATOM   1265   CZ ARG    21    -2.590 28.986 21.801 0.640 1.2628
ATOM   1266  NH1 ARG    21    -3.514 29.885 21.492 -0.800 1.8240
ATOM   1267  NH2 ARG    21    -1.762 29.206 22.813 -0.800 1.8240
ATOM   1268   O2 ARG    21    -5.807 23.410 16.918 -0.800 1.6612
END

""")
                    
                elif fname.endswith('.pdb'):
                    f.write("""ATOM      1  N   MET     1       0.551  11.263  11.815  1.00 32.15           N
ATOM      2  CA  MET     1       0.930  12.648  11.845  1.00 25.55           C
ATOM      3  C   MET     1       2.423  12.834  11.663  1.00 20.72           C
ATOM      4  O   MET     1       3.286  11.970  11.947  1.00 24.07           O
ATOM      5  CB  MET     1       0.580  13.368  13.143  1.00 31.39           C
ATOM      6  CG  MET     1       1.661  12.986  14.098  1.00 46.81           C
ATOM      7  SD  MET     1       1.441  13.562  15.776  1.00 59.99           S
ATOM      8  CE  MET     1       2.528  12.285  16.478  1.00 57.80           C
ATOM      9  N   ILE     2       2.656  14.070  11.163  1.00 20.22           N
ATOM     10  CA  ILE     2       3.942  14.570  10.821  1.00 20.46           C
ATOM     11  C   ILE     2       4.414  15.563  11.777  1.00 11.61           C
ATOM     12  O   ILE     2       3.666  16.500  12.106  1.00 16.02           O
ATOM     13  CB  ILE     2       3.780  15.262   9.428  1.00 27.69           C
ATOM     14  CG1 ILE     2       3.485  14.180   8.372  1.00 27.14           C
ATOM     15  CG2 ILE     2       4.987  16.173   9.033  1.00 21.01           C
ATOM     16  CD1 ILE     2       2.969  14.795   7.082  1.00 28.78           C
ATOM     17  N   SER     3       5.637  15.334  12.115  1.00 13.62           N
ATOM     18  CA  SER     3       6.368  16.180  13.050  1.00 19.09           C
ATOM     19  C   SER     3       7.703  16.646  12.468  1.00 17.57           C
ATOM     20  O   SER     3       8.363  15.914  11.750  1.00 15.93           O
ATOM     21  CB  SER     3       6.751  15.301  14.298  1.00 14.43           C
ATOM     22  OG  SER     3       5.544  14.967  14.987  1.00 15.13           O
ATOM     23  N   LEU     4       8.112  17.862  12.816  1.00 10.24           N
ATOM     24  CA  LEU     4       9.370  18.301  12.353  1.00 11.22           C
ATOM     25  C   LEU     4      10.283  18.416  13.583  1.00 12.15           C
ATOM     26  O   LEU     4       9.765  18.812  14.627  1.00 19.68           O
ATOM     27  CB  LEU     4       9.230  19.754  11.875  1.00 18.64           C
ATOM     28  CG  LEU     4       8.684  19.996  10.452  1.00 19.41           C
ATOM     29  CD1 LEU     4       7.569  19.115   9.957  1.00 16.26           C
ATOM     30  CD2 LEU     4       8.477  21.499  10.173  1.00 17.49           C
ATOM     31  N   ILE     5      11.555  18.104  13.495  1.00 11.49           N
ATOM     32  CA  ILE     5      12.465  18.230  14.624  1.00  9.79           C
ATOM     33  C   ILE     5      13.636  19.013  14.159  1.00 15.10           C
ATOM     34  O   ILE     5      14.222  18.681  13.116  1.00 17.21           O
ATOM     35  CB  ILE     5      12.811  16.919  15.273  1.00 13.29           C
ATOM     36  CG1 ILE     5      13.886  17.156  16.315  1.00 11.26           C
ATOM     37  CG2 ILE     5      13.298  15.880  14.184  1.00 16.94           C
ATOM     38  CD1 ILE     5      13.981  15.939  17.251  1.00 11.88           C
ATOM     39  N   ALA     6      14.012  20.093  14.895  1.00 10.27           N
ATOM     40  CA  ALA     6      15.137  20.934  14.463  1.00 11.73           C
ATOM     41  C   ALA     6      15.851  21.640  15.619  1.00 17.85           C
ATOM     42  O   ALA     6      15.235  21.839  16.668  1.00 16.34           O
ATOM     43  CB  ALA     6      14.449  22.032  13.594  1.00 11.92           C
ATOM     44  N   ALA     7      17.119  21.996  15.415  1.00 12.95           N
ATOM     45  CA  ALA     7      17.962  22.700  16.346  1.00 11.97           C
ATOM     46  C   ALA     7      18.099  24.081  15.714  1.00 26.40           C
ATOM     47  O   ALA     7      18.503  24.190  14.577  1.00 16.56           O
ATOM     48  CB  ALA     7      19.279  22.022  16.626  1.00 13.38           C
ATOM     49  N   LEU     8      17.738  25.159  16.435  1.00 13.46           N
ATOM     50  CA  LEU     8      17.767  26.510  15.918  1.00 11.94           C
ATOM     51  C   LEU     8      18.542  27.441  16.789  1.00 23.46           C
ATOM     52  O   LEU     8      18.344  27.380  17.993  1.00 15.97           O
ATOM     53  CB  LEU     8      16.306  27.133  16.150  1.00 17.43           C
ATOM     54  CG  LEU     8      15.245  26.745  15.191  1.00 26.42           C
ATOM     55  CD1 LEU     8      14.714  25.361  15.532  1.00 35.76           C
ATOM     56  CD2 LEU     8      14.133  27.753  15.245  1.00 21.53           C
ATOM     57  N   ALA     9      19.354  28.324  16.136  1.00 16.75           N
ATOM     58  CA  ALA     9      20.109  29.322  16.859  1.00 11.74           C
ATOM     59  C   ALA     9      19.268  30.562  16.781  1.00 14.87           C
ATOM     60  O   ALA     9      18.052  30.512  16.395  1.00 16.15           O
ATOM     61  CB  ALA     9      21.496  29.499  16.312  1.00 18.58           C
ATOM     62  N   VAL    10      19.908  31.695  17.172  1.00 19.77           N
ATOM     63  CA  VAL    10      19.182  32.944  17.121  1.00 28.89           C
ATOM     64  C   VAL    10      18.853  33.270  15.681  1.00 24.80           C
ATOM     65  O   VAL    10      19.589  32.940  14.704  1.00 24.85           O
ATOM     66  CB  VAL    10      19.896  34.147  17.775  1.00 35.16           C
ATOM     67  CG1 VAL    10      20.668  33.742  19.022  1.00 30.41           C
ATOM     68  CG2 VAL    10      20.692  34.984  16.815  1.00 35.93           C
ATOM     69  N   ASP    11      17.742  33.900  15.528  1.00 25.72           N
ATOM     70  CA  ASP    11      17.367  34.231  14.199  1.00 31.44           C
ATOM     71  C   ASP    11      16.970  33.019  13.371  1.00 28.60           C
ATOM     72  O   ASP    11      16.959  33.131  12.206  1.00 26.55           O
ATOM     73  CB  ASP    11      18.378  35.124  13.467  1.00 29.36           C
ATOM     74  CG  ASP    11      18.356  36.564  14.035  1.00 33.87           C
ATOM     75  OD1 ASP    11      17.333  37.176  14.327  1.00 36.53           O
ATOM     76  OD2 ASP    11      19.538  37.099  14.146  1.00 34.31           O
ATOM     77  N   ARG    12      16.634  31.897  13.959  1.00 15.13           N
ATOM     78  CA  ARG    12      16.193  30.711  13.247  1.00 13.52           C
ATOM     79  C   ARG    12      17.224  30.074  12.351  1.00 15.81           C
ATOM     80  O   ARG    12      16.841  29.226  11.518  1.00 19.91           O
ATOM     81  CB  ARG    12      14.891  30.911  12.489  1.00 14.72           C
ATOM     82  CG  ARG    12      13.770  31.437  13.345  1.00 24.78           C
ATOM     83  CD  ARG    12      12.435  31.161  12.721  1.00 27.21           C
ATOM     84  NE  ARG    12      12.462  31.677  11.381  1.00 36.54           N
ATOM     85  CZ  ARG    12      12.456  32.954  11.069  1.00 37.11           C
ATOM     86  NH1 ARG    12      12.421  33.916  11.940  1.00 36.05           N
ATOM     87  NH2 ARG    12      12.494  33.278   9.801  1.00 46.75           N
ATOM     88  N   VAL    13      18.498  30.450  12.518  1.00 14.21           N
ATOM     89  CA  VAL    13      19.531  29.853  11.712  1.00 12.04           C
ATOM     90  C   VAL    13      19.672  28.345  12.064  1.00 22.54           C
ATOM     91  O   VAL    13      19.821  27.968  13.233  1.00 14.87           O
ATOM     92  CB  VAL    13      20.863  30.533  11.990  1.00 16.62           C
ATOM     93  CG1 VAL    13      22.042  29.742  11.406  1.00 14.11           C
ATOM     94  CG2 VAL    13      20.797  31.998  11.467  1.00 16.15           C
ATOM     95  N   ILE    14      19.629  27.482  11.056  1.00 14.13           N
ATOM     96  CA  ILE    14      19.767  26.011  11.271  1.00 14.28           C
ATOM     97  C   ILE    14      20.995  25.409  10.552  1.00 17.38           C
ATOM     98  O   ILE    14      21.399  24.278  10.747  1.00 16.61           O
ATOM     99  CB  ILE    14      18.506  25.235  10.843  1.00 19.53           C
ATOM    100  CG1 ILE    14      18.151  25.494   9.365  1.00 22.73           C
ATOM    101  CG2 ILE    14      17.307  25.571  11.723  1.00 14.60           C
ATOM    102  CD1 ILE    14      16.862  24.785   8.834  1.00 17.31           C
ATOM    103  N   GLY    15      21.623  26.185   9.656  1.00 21.60           N
ATOM    104  CA  GLY    15      22.773  25.644   8.908  1.00 17.53           C
ATOM    105  C   GLY    15      23.579  26.679   8.113  1.00 18.53           C
ATOM    106  O   GLY    15      23.216  27.824   8.023  1.00 16.92           O
ATOM    107  N   MET    16      24.685  26.243   7.585  1.00 22.71           N
ATOM    108  CA  MET    16      25.630  27.043   6.835  1.00 23.71           C
ATOM    109  C   MET    16      26.562  26.064   6.140  1.00 28.63           C
ATOM    110  O   MET    16      26.419  24.841   6.284  1.00 19.96           O
ATOM    111  CB  MET    16      26.454  27.903   7.783  1.00 28.96           C
ATOM    112  CG  MET    16      27.215  27.034   8.748  1.00 39.22           C
ATOM    113  SD  MET    16      28.614  27.916   9.432  1.00 47.53           S
ATOM    114  CE  MET    16      29.524  27.955   7.874  1.00 33.80           C
ATOM    115  N   GLU    17      27.509  26.591   5.380  1.00 28.41           N
ATOM    116  CA  GLU    17      28.375  25.714   4.685  1.00 29.16           C
ATOM    117  C   GLU    17      29.197  24.820   5.561  1.00 24.40           C
ATOM    118  O   GLU    17      29.293  23.610   5.320  1.00 25.76           O
ATOM    119  CB  GLU    17      29.332  26.532   3.738  0.00 20.00           C
ATOM    120  CG  GLU    17      28.646  27.132   2.520  0.00 20.00           C
ATOM    121  CD  GLU    17      28.159  26.083   1.541  0.00 20.00           C
ATOM    122  OE1 GLU    17      28.336  24.880   1.823  0.00 20.00           O
ATOM    123  OE2 GLU    17      27.599  26.466   0.492  0.00 20.00           O
ATOM    124  N   ASN    18      29.817  25.419   6.557  1.00 20.87           N
ATOM    125  CA  ASN    18      30.671  24.623   7.428  1.00 18.24           C
ATOM    126  C   ASN    18      29.968  24.072   8.644  1.00 25.72           C
ATOM    127  O   ASN    18      28.800  24.344   8.861  1.00 22.25           O
ATOM    128  CB  ASN    18      31.891  25.429   7.870  1.00 24.80           C
ATOM    129  CG  ASN    18      32.673  26.009   6.672  1.00 35.68           C
ATOM    130  ND2 ASN    18      32.793  25.192   5.612  1.00 21.02           N
ATOM    131  OD1 ASN    18      33.131  27.181   6.693  1.00 43.36           O
ATOM    132  N   ALA    19      30.721  23.301   9.419  1.00 24.77           N
ATOM    133  CA  ALA    19      30.158  22.732  10.649  1.00 23.26           C
ATOM    134  C   ALA    19      29.808  23.892  11.608  1.00 16.44           C
ATOM    135  O   ALA    19      30.553  24.855  11.718  1.00 19.16           O
ATOM    136  CB  ALA    19      31.208  21.851  11.326  1.00 17.26           C
ATOM    137  N   MET    20      28.661  23.767  12.282  1.00 22.16           N
ATOM    138  CA  MET    20      28.173  24.753  13.256  1.00 27.06           C
ATOM    139  C   MET    20      29.215  24.880  14.367  1.00 16.21           C
ATOM    140  O   MET    20      29.858  23.879  14.759  1.00 20.00           O
ATOM    141  CB  MET    20      26.859  24.207  13.808  1.00 25.41           C
ATOM    142  CG  MET    20      25.790  24.287  12.739  1.00 29.82           C
ATOM    143  SD  MET    20      25.744  25.961  12.033  1.00 41.04           S
ATOM    144  CE  MET    20      24.199  26.587  12.680  1.00 55.69           C
ATOM   1257  N   ARG    21      -2.365  24.238  17.120  1.00 26.07           N
ATOM   1258  CA  ARG    21      -3.722  24.429  17.637  1.00 38.04           C
ATOM   1259  C   ARG    21      -4.816  24.158  16.587  1.00 52.67           C
ATOM   1260  O   ARG    21      -4.718  24.783  15.486  1.00 56.66           O
ATOM   1261  CB  ARG    21      -3.918  25.826  18.162  1.00 35.05           C
ATOM   1262  CG  ARG    21      -3.019  26.176  19.341  1.00 32.34           C
ATOM   1263  CD  ARG    21      -3.363  27.548  19.971  1.00 27.62           C
ATOM   1264  NE  ARG    21      -2.493  27.864  21.095  0.00 20.00           N
ATOM   1265  CZ  ARG    21      -2.590  28.986  21.801  0.00 20.00           C
ATOM   1266  NH1 ARG    21      -3.514  29.885  21.492  0.00 20.00           N
ATOM   1267  NH2 ARG    21      -1.762  29.206  22.813  0.00 20.00           N
ATOM   1268  OXT ARG    21      -5.807  23.410  16.918  1.00 51.94           O
END
""")
                f.close()
        cmd = cmd()
    pymol = pymol()

def __init__(self):
    if 0:
        self.menuBar.addcascademenu('Plugin', 'APBSTools', 'APBS Tools',
                                    label='APBS Tools')

        self.menuBar.addmenuitem('APBSTools', 'command',
                                 'Launch APBS Tools',
                                 label='APBS Tools...',
                                 command = lambda s=self: APBSTools(s))


        self.menuBar.addmenuitem('APBSTools', 'command',
                                 'Electrostatics Wizard',
                                 label='Electrostatics Wizard',
                                 command = lambda : pymol.cmd.wizard("electrostatics"))
    else:
        self.menuBar.addmenuitem('Plugin', 'command',
                                 'Launch APBS Tools',
                                 label='APBS Tools...',
                                 command = lambda s=self: APBSTools(s))
        

   

defaults = {
    "interior_dielectric" : 2.0,# ###
    "solvent_dielectric" : 80.0,# ###
    "solvent_radius" : 1.4,# ###
    "system_temp" : 310.0,# ###
    "apbs_mode" : 'Nonlinear Poisson-Boltzmann Equation',# ###
    "ion_plus_one_conc" : 0.0,# ###
    "ion_plus_one_rad" : 2.0,# ###
    "ion_plus_two_conc" : 0.0,# ###
    "ion_plus_two_rad" : 2.0,# ###
    "ion_minus_one_conc" : 0.0,# ###
    "ion_minus_one_rad" : 2.0,# ###
    "ion_minus_two_conc" : 0.0,# ###
    "ion_minus_two_rad" : 2.0,# ###
    "max_grid_points" : 1200000,
    "potential_at_sas" : 0,
    "surface_solvent" : 0,
    #"grid_buffer" : 0,
    #"grid_buffer" : 20,
    "bcfl" : 'Single DH sphere', # Boundary condition flag for APBS ###
    "sdens": 10.0, # Specify the number of grid points per square-angstrom to use in Vacc object. Ignored when srad is 0.0 (see srad) or srfm is spl2 (see srfm). There is a direct correlation between the value used for the Vacc sphere density, the accuracy of the Vacc object, and the APBS calculation time. APBS default value is 10.0.
    "chgm" : 'Cubic b-splines', # Charge disc method for APBS ###
    }

class FileDialogButtonClassFactory:
    def get(fn,filter='*'):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class FileDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.

            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                apply(Tkinter.Button.__init__, (self, master, cnf), kw)
                self.configure(command=self.set)
            def set(self):
                fd = PmwFileDialog(self.master,filter=filter)
                fd.title('Please choose a file')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)

class VisualizationGroup(Pmw.Group):
    def __init__(self,*args,**kwargs):
        Pmw.Group.__init__(self,*args,**kwargs)
        self.refresh()
        self.show_ms = False
        self.show_pi = False
        self.show_ni = False
    def refresh(self):
        things_to_kill = 'error_label update_buttonbox mm_group ms_group pi_group ni_group'.split()
        for thing in things_to_kill:
            try:
                getattr(self,thing).destroy()
                #print "destroyed",thing
            except AttributeError:
                #print "couldn't destroy",thing

                # note: this attributeerror will also hit if getattr(self,thing) misses.
                # another note: both halves of the if/else make an update_buttonbox.
                # if you rename the one in the top half to something else, you'll cause nasty Pmw errors.
                pass

        if [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:map'] and [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:molecule']:

            self.mm_group = Pmw.Group(self.interior(),tag_text = 'Maps and Molecules')
            self.map = Pmw.OptionMenu(self.mm_group.interior(),
                                      labelpos = 'w',
                                      label_text = 'Map',
                                      items = [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:map'],
                                      )
            self.map.pack(padx=4,side=LEFT)

            self.molecule = Pmw.OptionMenu(self.mm_group.interior(),
                                           labelpos = 'w',
                                           label_text = 'Molecule',
                                           items = [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:molecule'],
                                           )
            self.molecule.pack(padx=4,side=LEFT)
            self.update_buttonbox = Pmw.ButtonBox(self.mm_group.interior(), padx=0)
            self.update_buttonbox.pack(side=LEFT)
            self.update_buttonbox.add('Update',command=self.refresh)
            self.mm_group.pack(fill = 'both', expand = 1, padx = 4, pady = 5, side=TOP)

            self.ms_group = Pmw.Group(self.interior(),tag_text='Molecular Surface')
            self.ms_buttonbox = Pmw.ButtonBox(self.ms_group.interior(), padx=0)
            self.ms_buttonbox.pack()
            self.ms_buttonbox.add('Show',command=self.showMolSurface)
            self.ms_buttonbox.add('Hide',command=self.hideMolSurface)
            self.ms_buttonbox.add('Update',command=self.updateMolSurface)
            self.ms_buttonbox.alignbuttons()
            self.surface_solvent = IntVar()
            self.surface_solvent.set(defaults['surface_solvent'])
            self.sol_checkbutton = Checkbutton(self.ms_group.interior(),
                                               text = "Solvent accesssible surface",
                                               variable = self.surface_solvent)
            self.sol_checkbutton.pack()            
            self.potential_at_sas = IntVar()
            self.potential_at_sas.set(defaults['potential_at_sas'])
            self.pot_checkbutton = Checkbutton(self.ms_group.interior(),
                                               text = "Color by potential on sol. acc. surf.",
                                               variable = self.potential_at_sas)
            self.pot_checkbutton.pack()
            self.mol_surf_low = Pmw.Counter(self.ms_group.interior(),
                                            labelpos = 'w',
                                            label_text = 'Low',
                                            orient = 'vertical',
                                            entry_width = 4,
                                            entryfield_value = -1,
                                            datatype = 'real',
                                            entryfield_validate = {'validator' : 'real'},
                                            )
            self.mol_surf_middle = Pmw.Counter(self.ms_group.interior(),
                                            labelpos = 'w',
                                            label_text = 'Middle',
                                            orient = 'vertical',
                                            entry_width = 4,
                                            entryfield_value = 0,
                                            datatype = 'real',
                                            entryfield_validate = {'validator' : 'real'}
                                            )
            self.mol_surf_high = Pmw.Counter(self.ms_group.interior(),
                                            labelpos = 'w',
                                            label_text = 'High',
                                            orient = 'vertical',
                                            entry_width = 4,
                                            entryfield_value = 1,
                                            datatype = 'real',
                                            entryfield_validate = {'validator' : 'real'}
                                            )
            bars = (self.mol_surf_low,self.mol_surf_middle,self.mol_surf_high)
            Pmw.alignlabels(bars)
            for bar in bars: bar.pack(side=LEFT)
            self.ms_group.pack(fill = 'both', expand = 1, padx = 4, pady = 5, side=LEFT)

            self.pi_group = Pmw.Group(self.interior(),tag_text='Positive Isosurface')
            self.pi_buttonbox = Pmw.ButtonBox(self.pi_group.interior(), padx=0)
            self.pi_buttonbox.pack()
            self.pi_buttonbox.add('Show',command=self.showPosSurface)
            self.pi_buttonbox.add('Hide',command=self.hidePosSurface)
            self.pi_buttonbox.add('Update',command=self.updatePosSurface)
            self.pi_buttonbox.alignbuttons()
            self.pos_surf_val = Pmw.Counter(self.pi_group.interior(),
                                            labelpos = 'w',
                                            label_text = 'Contour (kT)',
                                            orient = 'vertical',
                                            entry_width = 4,
                                            entryfield_value = 1,
                                            datatype = 'real',
                                            entryfield_validate = {'validator' : 'real', 'min':0}
                                            )
            self.pos_surf_val.pack(side=LEFT)
            self.pi_group.pack(fill = 'both', expand = 1, padx = 4, pady = 5, side=LEFT)

            self.ni_group = Pmw.Group(self.interior(),tag_text='Negative Isosurface')
            self.ni_buttonbox = Pmw.ButtonBox(self.ni_group.interior(), padx=0)
            self.ni_buttonbox.pack()
            self.ni_buttonbox.add('Show',command=self.showNegSurface)
            self.ni_buttonbox.add('Hide',command=self.hideNegSurface)
            self.ni_buttonbox.add('Update',command=self.updateNegSurface)
            self.ni_buttonbox.alignbuttons()
            self.neg_surf_val = Pmw.Counter(self.ni_group.interior(),
                                            labelpos = 'w',
                                            label_text = 'Contour (kT)',
                                            orient = 'vertical',
                                            entry_width = 4,
                                            entryfield_value = -1,
                                            datatype = 'real',
                                            entryfield_validate = {'validator' : 'real', 'max':0}
                                            )
            self.neg_surf_val.pack(side=LEFT)
            self.ni_group.pack(fill = 'both', expand = 1, padx = 4, pady = 5, side=LEFT)

        else:
            self.error_label = Tkinter.Label(self.interior(),
                                  pady = 10,
                                  justify=LEFT,
                                  text = '''You must have at least a molecule and a map loaded.
If you have a molecule and a map loaded, please click "Update"''',
                                  )
            self.error_label.pack()
            self.update_buttonbox = Pmw.ButtonBox(self.interior(), padx=0)
            self.update_buttonbox.pack()
            self.update_buttonbox.add('Update',command=self.refresh)

    def showMolSurface(self):
        self.updateMolSurface()
            
    def hideMolSurface(self):
        pymol.cmd.hide('surface',self.molecule.getvalue())

    def updateMolSurface(self):
        ramp_name = 'e_lvl'
        map_name = self.map.getvalue()
        molecule_name = self.molecule.getvalue()
        low = float(self.mol_surf_low.getvalue())
        mid = float(self.mol_surf_middle.getvalue())
        high = float(self.mol_surf_high.getvalue())
        range = [low,mid,high]
        print " APBS Tools: range is",range
        pymol.cmd.delete(ramp_name)
        pymol.cmd.ramp_new(ramp_name,map_name,range)
        pymol.cmd.set('surface_color',ramp_name,molecule_name)
        if self.surface_solvent.get()==1:
            pymol.cmd.set('surface_solvent',1,molecule_name)
            pymol.cmd.set('surface_ramp_above_mode',0,molecule_name)
        else:
            pymol.cmd.set('surface_solvent',0,molecule_name)
            pymol.cmd.set('surface_ramp_above_mode',self.potential_at_sas.get(),molecule_name)
        pymol.cmd.show('surface',molecule_name)
        pymol.cmd.refresh()
        pymol.cmd.recolor()
        
    def showPosSurface(self):
        self.updatePosSurface()
    def hidePosSurface(self):
        pymol.cmd.hide('everything','iso_pos')
    def updatePosSurface(self):
        pymol.cmd.delete('iso_pos')
        pymol.cmd.isosurface('iso_pos',self.map.getvalue(),float(self.pos_surf_val.getvalue()))
        pymol.cmd.color('blue','iso_pos')
        pymol.cmd.show('everything','iso_pos')
    def showNegSurface(self):
        self.updateNegSurface()
    def hideNegSurface(self):
        pymol.cmd.hide('everything','iso_neg')
    def updateNegSurface(self):
        pymol.cmd.delete('iso_neg')
        pymol.cmd.isosurface('iso_neg',self.map.getvalue(),float(self.neg_surf_val.getvalue()))
        pymol.cmd.color('red','iso_neg')
        pymol.cmd.show('everything','iso_neg')
        
class APBSTools:

    def setPqrFile(self,name):
        print " APBS Tools: set pqr file to",name
        self.pqr_to_use.setvalue(name)
    def getPqrFilename(self):
        if self.radiobuttons.getvalue() != 'Use another PQR':
#            print 'radiobutton said to generate it',self.radiobuttons.getvalue(),'so i am returning',self.pymol_generated_pqr_filename.getvalue()
            return self.pymol_generated_pqr_filename.getvalue()
        else:
#            print 'radiobutton told me to use another',self.radiobuttons.getvalue(),'so i am returning',self.pqr_to_use.getvalue()
            return self.pqr_to_use.getvalue()
    def setPsizeLocation(self,value):
        self.psize.setvalue(value)
    def setBinaryLocation(self,value):
        self.binary.setvalue(value)

    def __init__(self,app):
        parent = app.root
        self.parent = parent
        
        # Create the dialog.
        self.dialog = Pmw.Dialog(parent,
                                 buttons = ('Run APBS','Set grid', 'Exit APBS tools'),
                                 #defaultbutton = 'Run APBS',
                                 title = 'PyMOL APBS Tools',
                                 command = self.execute)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = Tkinter.Label(self.dialog.interior(),
                                text = 'PyMOL APBS Tools\nMichael Lerner, 2004 - www.umich.edu/~mlerner/Pymol\n(incorporates modifications by WLD)',
                                background = 'black',
                                foreground = 'white',
                                #pady = 20,
                                )
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',expand=1,padx=10,pady=10)

        # Set up the Main page
        page = self.notebook.add('Main')
        group = Pmw.Group(page,tag_text='Main options')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.selection = Pmw.EntryField(group.interior(),
                                        labelpos='w',
                                        label_text='Selection to use: ',
                                        value='(polymer)',
                                        )
        self.map = Pmw.EntryField(group.interior(),
                                  labelpos='w',
                                  label_text='What to call the resulting map: ',
                                  value='apbs_map',
                                  )
        self.radiobuttons = Pmw.RadioSelect(group.interior(),
                                            buttontype = 'radiobutton',
                                            orient = 'vertical',
                                            labelpos = 'w',
                                            )
        for text in ('Use PyMOL generated PQR and existing Hydrogens and termini',
                     'Use PyMOL generated PQR and PyMOL generated Hydrogens and termini',
                     'Use another PQR'):
            self.radiobuttons.add(text)
        self.radiobuttons.setvalue('Use PyMOL generated PQR and PyMOL generated Hydrogens and termini')

        self.pqr_to_use = Pmw.EntryField(group.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.setPqrFile,'*.pqr'),
                                         label_text='\tChoose Externally Generated PQR:',
                                         )
        

        for entry in (self.selection,self.map,self.radiobuttons,self.pqr_to_use):
        #for entry in (self.selection,self.map,self.radiobuttons,):
            entry.pack(fill='x',padx=4,pady=1) # vertical
        

        # Set up the main "Calculation" page
        page = self.notebook.add('Configuration')


        group = Pmw.Group(page,tag_text='Dielectric Constants')
        group.pack(fill = 'both', expand = 1, padx = 4, pady = 5)
        group.grid(column=0, row=0)
        self.interior_dielectric = Pmw.EntryField(group.interior(),labelpos='w',
                                   label_text = 'Protein Dielectric:',
                                   value = str(defaults['interior_dielectric']),
                                   validate = {'validator' : 'real',
                                               'min':0,}
                                   )
        self.solvent_dielectric = Pmw.EntryField(group.interior(),labelpos='w',
                                   label_text = 'Solvent Dielectric:',
                                   value = str(defaults['solvent_dielectric']),
                                   validate = {'validator' : 'real',
                                               'min':0,}
                                   )
        entries = (self.interior_dielectric,self.solvent_dielectric)
        for entry in entries:
            #entry.pack(side='left',fill='both',expand=1,padx=4) # side-by-side
            entry.pack(fill='x',expand=1,padx=4,pady=1) # vertical

        group = Pmw.Group(page,tag_text='Other')
        group.pack(fill='both',expand=1, padx=4, pady=5)
        group.grid(column=1, row=1,columnspan=4)
        self.max_grid_points = Pmw.EntryField(group.interior(),labelpos='w',
                                              label_text = 'Maximum Grid Points:',
                                              value = str(defaults['max_grid_points']),
                                              validate = {'validator' : 'real',
                                                          'min':1,}
                                              )
        self.apbs_mode = Pmw.OptionMenu(group.interior(),
                                        labelpos = 'w',
                                        label_text = 'APBS Mode',
                                        items = ('Nonlinear Poisson-Boltzmann Equation','Linear Poisson-Boltzmann Equation',),
                                        initialitem = defaults['apbs_mode'],
                                        )
        self.apbs_mode.pack(fill='x',expand=1,padx=4)
        #self.apbs_mode.grid(column=0,row=0,columnspan=3)
        self.solvent_radius = Pmw.EntryField(group.interior(),
                              labelpos = 'w',
                              label_text = 'Solvent Radius:',
                              validate = {'validator':'real','min':0},
                              value = str(defaults['solvent_radius']),
                              )
        self.system_temp = Pmw.EntryField(group.interior(),
                              labelpos = 'w',
                              label_text = 'System Temperature:',
                              validate = {'validator':'real','min':0},
                              value = str(defaults['system_temp']),
                              )
        self.sdens = Pmw.EntryField(group.interior(),
                                    labelpos = 'w',
                                    label_text = 'Vacc sphere density (grid points/A^2)',
                                    validate = {'validator':'real','min':0},
                                    value = str(defaults['sdens']),
                                    )
        self.bcfl = Pmw.OptionMenu(group.interior(),
                                   labelpos = 'w',
                                   label_text = 'Boundary Condition',
                                   items = ('Zero','Single DH sphere','Multiple DH spheres','Focusing',),
                                   initialitem = defaults['bcfl'],
                                   )
        self.chgm = Pmw.OptionMenu(group.interior(),
                                   labelpos = 'w',
                                   label_text = 'Charge disc method',
                                   items = ('Linear','Cubic b-splines',),
                                   initialitem = defaults['chgm'],
                                   )
        self.srfm = Pmw.OptionMenu(group.interior(),
                                   labelpos = 'w',
                                   label_text = 'Surface Calculation Method',
                                   items = ('Mol surf for epsilon; inflated VdW for kappa, no smoothing','Same, but with harmonic average smoothing','Cubic spline',),
                                   initialitem = 'Same, but with harmonic average smoothing',
                                   )
        #for entry in (self.apbs_mode,self.system_temp,self.solvent_radius,):
        for entry in (self.max_grid_points,self.solvent_radius,self.system_temp,self.sdens,self.apbs_mode,self.bcfl,self.chgm,self.srfm):
            entry.pack(fill='x',expand=1,padx=4,pady=1) # vertical


        group = Pmw.Group(page,tag_text='Ions')
        group.pack(fill='both',expand=1, padx=4, pady=5)
        group.grid(column=0, row=1, )
        self.ion_plus_one_conc = Pmw.EntryField(group.interior(),
                                                labelpos='w',
                                                label_text='Ion Concentration (+1):',
                                                validate={'validator':'real','min':0},
                                                value = str(defaults['ion_plus_one_conc']),
                                                )
        self.ion_plus_one_rad = Pmw.EntryField(group.interior(),
                                               labelpos='w',
                                               label_text='Ion Radius (+1):',
                                               validate={'validator':'real','min':0},
                                               value = str(defaults['ion_plus_one_rad']),
                                               )
        self.ion_minus_one_conc = Pmw.EntryField(group.interior(),
                                                 labelpos='w',
                                                 label_text='Ion Concentration (-1):',
                                                 validate={'validator':'real','min':0},
                                                 value = str(defaults['ion_minus_one_conc']),
                                                 )
        self.ion_minus_one_rad = Pmw.EntryField(group.interior(),
                                                labelpos='w',
                                                label_text='Ion Radius (-1):',
                                                validate={'validator':'real','min':0},
                                                value = str(defaults['ion_minus_one_rad']),
                                                )
        self.ion_plus_two_conc = Pmw.EntryField(group.interior(),
                                                labelpos='w',
                                                label_text='Ion Concentration (+2):',
                                                validate={'validator':'real','min':0},
                                                value = str(defaults['ion_plus_two_conc']),
                                                )
        self.ion_plus_two_rad = Pmw.EntryField(group.interior(),
                                               labelpos='w',
                                               label_text='Ion Radius (+2):',
                                               validate={'validator':'real','min':0},
                                               value = str(defaults['ion_plus_two_rad']),
                                               )
        self.ion_minus_two_conc = Pmw.EntryField(group.interior(),
                                                 labelpos='w',
                                                 label_text='Ion Concentration (-2):',
                                                 validate={'validator':'real','min':0},
                                                 value = str(defaults['ion_minus_two_conc']),
                                                 )
        self.ion_minus_two_rad = Pmw.EntryField(group.interior(),
                                                labelpos='w',
                                                label_text='Ion Radius (-2):',
                                                validate={'validator':'real','min':0},
                                                value = str(defaults['ion_minus_two_rad']),
                                                )
        entries = (self.ion_plus_one_conc,self.ion_plus_one_rad,
                      self.ion_minus_one_conc,self.ion_minus_one_rad,
                      self.ion_plus_two_conc,self.ion_plus_two_rad,
                      self.ion_minus_two_conc,self.ion_minus_two_rad,
                      )
        for entry in entries:
            entry.pack(fill='x',expand=1,padx=4)

        group = Pmw.Group(page,tag_text = 'Coarse Mesh Length')
        group.pack(fill = 'both', expand = 1, padx = 4, pady = 5)
        group.grid(column = 1, row = 0)
        for coord in 'x y z'.split():
            setattr(self,'grid_coarse_%s'%coord,Pmw.EntryField(group.interior(),
                                                               labelpos='w',
                                                               label_text=coord,
                                                               validate={'validator':'real','min':0},
                                                               value = -1,
                                                               entry_width=15,
                                                               )
                    )
            getattr(self,'grid_coarse_%s'%coord).pack(fill='x', expand=1, padx=4, pady=1)


        group = Pmw.Group(page,tag_text = 'Fine Mesh Length')
        group.pack(fill = 'both', expand = 1, padx = 4, pady = 5)
        group.grid(column = 2, row = 0)
        for coord in 'x y z'.split():
            setattr(self,'grid_fine_%s'%coord,Pmw.EntryField(group.interior(),
                                                               labelpos='w',
                                                               label_text=coord,
                                                               validate={'validator':'real','min':0},
                                                               value = -1,
                                                               entry_width=15,
                                                               )
                    )
            getattr(self,'grid_fine_%s'%coord).pack(fill='x', expand=1, padx=4, pady=1)


        group = Pmw.Group(page,tag_text = 'Grid Center')
        group.pack(fill = 'both', expand = 1, padx = 4, pady = 5)
        group.grid(column = 3, row = 0)
        for coord in 'x y z'.split():
            setattr(self,'grid_center_%s'%coord,Pmw.EntryField(group.interior(),
                                                               labelpos='w',
                                                               label_text=coord,
                                                               validate={'validator':'real'},
                                                               value = 0,
                                                               entry_width=10,
                                                               )
                    )
            getattr(self,'grid_center_%s'%coord).pack(fill='x', expand=1, padx=4, pady=1)

        group = Pmw.Group(page,tag_text = 'Grid Points')
        group.pack(fill = 'both', expand = 1, padx = 4, pady = 5)
        group.grid(column = 4, row = 0)
        for coord in 'x y z'.split():
            setattr(self,'grid_points_%s'%coord,Pmw.EntryField(group.interior(),
                                                               labelpos='w',
                                                               label_text=coord,
                                                               validate={'validator':'integer','min':0},
                                                               value = -1,
                                                               entry_width=8,
                                                               )
                    )
            getattr(self,'grid_points_%s'%coord).pack(fill='x', expand=1, padx=4, pady=1)


        page.grid_rowconfigure(2,weight=1)
        page.grid_columnconfigure(5,weight=1)
        page = self.notebook.add('APBS Location')
        group = Pmw.Group(page,tag_text='Locations')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        def quickFileValidation(s):
            if s == '': return Pmw.PARTIAL
            elif os.path.isfile(s): return Pmw.OK
            elif os.path.exists(s): return Pmw.PARTIAL
            else: return Pmw.PARTIAL

        global APBS_BINARY_LOCATION, APBS_PSIZE_LOCATION
        if APBS_BINARY_LOCATION is None:
            if 'APBS_BINARY' in os.environ:
                APBS_BINARY_LOCATION = os.environ['APBS_BINARY']
            else:
                APBS_BINARY_LOCATION = distutils.spawn.find_executable('apbs')
                if APBS_BINARY_LOCATION is None:
                    APBS_BINARY_LOCATION = ''

        if APBS_PSIZE_LOCATION is None:
            if 'APBS_PSIZE' in os.environ:
                APBS_PSIZE_LOCATION = os.environ['APBS_PSIZE']
            else:
                APBS_PSIZE_LOCATION = ''

        self.binary = Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_pyclass = FileDialogButtonClassFactory.get(self.setBinaryLocation),
                                     validate = {'validator':quickFileValidation,},
                                     value = APBS_BINARY_LOCATION,
                                     label_text = 'APBS binary location:')
        self.binary.pack(fill = 'x', padx = 20, pady = 10)
        self.psize =  Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_pyclass = FileDialogButtonClassFactory.get(self.setPsizeLocation),
                                     validate = {'validator':quickFileValidation,},
                                     #value = '/usr/local/apbs-0.3.1/tools/manip/psize.py',
                                     value = APBS_PSIZE_LOCATION,
                                     label_text = 'APBS psize.py location:',
                                     )
        self.psize.pack(fill = 'x', padx = 20, pady = 10)
        label = Tkinter.Label(group.interior(),
                              pady = 10,
                              justify=LEFT,
                              text = """You must have APBS installed on your system.

The PyMOL APBS tools can calculate proper grid dimensions and spacing (we ensure the the
fine mesh spacing is 0.5A or finer).  If you wish to use APBS's psize.py to set up the
grid, make sure that the path is set correctly above.

If PyMOL does not automatically find apbs or psize.py, you may set the environment variables
APBS_BINARY and APBS_PSIZE to point to them respectively.
""",
                              )
        label.pack()
        
        page = self.notebook.add('Temporary File Locations')
        group = Pmw.Group(page,tag_text='Locations')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.pymol_generated_pqr_filename = Pmw.EntryField(group.interior(),
                                                           labelpos = 'w',
                                                           label_text = 'Temporary PQR file: ',
                                                           value = 'pymol-generated.pqr',
                                                           )
        self.pymol_generated_pqr_filename.pack(fill = 'x', padx = 20, pady = 10)
        
        self.pymol_generated_pdb_filename = Pmw.EntryField(group.interior(),
                                                           labelpos = 'w',
                                                           label_text = 'Temporary PDB file: ',
                                                           value = 'pymol-generated.pdb',
                                                           )
        self.pymol_generated_pdb_filename.pack(fill = 'x', padx = 20, pady = 10)

        self.pymol_generated_dx_filename = Pmw.EntryField(group.interior(),
                                                          labelpos = 'w',
                                                          label_text = 'Temporary DX file: ',
                                                          value = 'pymol-generated.dx',
                                                          )
        self.pymol_generated_dx_filename.pack(fill = 'x', padx = 20, pady = 10)

        self.pymol_generated_in_filename = Pmw.EntryField(group.interior(),
                                                          labelpos = 'w',
                                                          label_text = 'APBS input file: ',
                                                          value = 'pymol-generated.in',
                                                          )
        self.pymol_generated_in_filename.pack(fill = 'x', padx = 20, pady = 10)

        # Create a visualization page
        page = self.notebook.add('Visualization')
        #group = Pmw.Group(page,tag_text='Visualization')
        #group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        group = VisualizationGroup(page,tag_text='Visualization')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        
        

        # Create a couple of other empty pages
        page = self.notebook.add('About')
        group = Pmw.Group(page, tag_text='About PyMOL APBS Tools')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        label = Tkinter.Label(group.interior(),
                              pady = 10,
                              justify=LEFT,
                              text = """
This plugin integrates PyMOL (http://pymol.org/) with APBS (http://agave.wustl.edu/apbs/).

It should be fairly self-explanatory.  In the simplest case,

1) Load a structure into PyMOL.
2) Start this plugin.
3) Make sure that the path to the APBS binary is correct on the "APBS Location" tab.
4) Click the "Set grid" button to set up the grid.
5) Click the "Run APBS" button.

Many thanks to
 - Warren DeLano for everything involving PyMOL
 - Nathan Baker and Todd Dolinsky for everything involving APBS

Created by Michael Lerner <http://www.umich.edu/~mlerner/Pymol/> mglerner@gmail.com
Carlson Group, University of Michigan <http://www.umich.edu/~carlsonh/>
""")
        label.pack()

        self.notebook.setnaturalsize()





        self.showAppModal()
         
    def showAppModal(self):
        #self.dialog.activate(geometry = 'centerscreenalways', globalMode = 'nograb')
        self.dialog.show()
        #self.dialog.activate(geometry = 'centerscreenalways')

    def execute(self, result):
        if result == 'Run APBS':
            good = self.generateApbsInputFile()
            if not good:
                return False
            if self.radiobuttons.getvalue() != 'Use another PQR':
                good = self.generatePqrFile()
                if not good:
                    return False
            if os.path.exists(self.pymol_generated_dx_filename.getvalue()):
                try:
                    os.unlink(self.pymol_generated_dx_filename.getvalue())
                except:
                    traceback.print_exc()
                    pass
            command = "%s %s" % (self.binary.getvalue(),self.pymol_generated_in_filename.getvalue())
            os.system(command)
            pymol.cmd.load(self.pymol_generated_dx_filename.getvalue(),self.map.getvalue())
        elif result == 'Set grid':
            self.runPsize()
        else:
            #
            # Doing it this way takes care of clicking on the x in the top of the
            # window, which as result set to None.
            #
            if __name__ == '__main__':
                #
                # dies with traceback, but who cares
                #
                self.parent.destroy()
            else:
                #self.dialog.deactivate(result)
                global APBS_BINARY_LOCATION, APBS_PSIZE_LOCATION
                APBS_BINARY_LOCATION = self.binary.getvalue()
                APBS_PSIZE_LOCATION = self.psize.getvalue()
                self.dialog.withdraw()

    def runPsize(self):
        class NoPsize(Exception):
            pass
        class NoPDB(Exception):
            pass
        try:
            if not self.psize.valid():
                raise NoPsize
            pdb_filename = self.pymol_generated_pdb_filename.getvalue()
            try:
                f = file(pdb_filename,'w')
                f.close()
            except:
                raise NoPDB
            #
            # Do some magic to load the psize module
            #
            import imp
            f,fname,description = imp.find_module('psize',[os.path.split(self.psize.getvalue())[0]])
            psize = imp.load_module('psize',f,fname,description)
            # WLD
            sel = "((%s) or (neighbor (%s) and hydro))"%(
                           self.selection.getvalue(), self.selection.getvalue())
            pymol.cmd.save(pdb_filename,sel)
            f.close()
            size = psize.Psize()
            size.runPsize(pdb_filename)
            coarsedim = size.getCoarseGridDims() # cglen
            finedim = size.getFineGridDims() # fglen
            # could use procgrid for multiprocessors
            finegridpoints = size.getFineGridPoints() # dime
            center = size.getCenter() # cgcent and fgcent
            
        except (NoPsize,ImportError):
            #
            # First, we need to get the dimensions of the molecule
            #
            # WLD
            sel = "((%s) or (neighbor (%s) and hydro))"%(
                           self.selection.getvalue(), self.selection.getvalue())
            model = pymol.cmd.get_model(sel)
            mins = [None,None,None]
            maxs = [None,None,None]
            for a in model.atom:
                for i in (0,1,2):
                    if mins[i] is None or a.coord[i] < mins[i]:
                        mins[i] = a.coord[i]
                    if maxs[i] is None or a.coord[i] > maxs[i]:
                        maxs[i] = a.coord[i]
            if None in mins or None in maxs:
                error_dialog = Pmw.MessageDialog(self.parent,
                                                 title = 'Error',
                                                 message_text = "No atoms were in your selection",
                                                 )
                junk = error_dialog.activate()
                return False
                
            box_length = [maxs[i] - mins[i] for i in range(3)]
            center = [(maxs[i] + mins[i])/2.0 for i in range(3)]
            #
            # psize expands the molecular dimensions by CFAC (which defaults
            # to 1.7) for the coarse grid
            #
            CFAC = 1.7
            coarsedim = [length*CFAC for length in box_length]

            #
            # psize also does something strange .. it adds a buffer FADD to
            # the box lengths to get the fine lengths.  you'd think it'd also
            # have FFAC or CADD, but we'll mimic it here.  it also has the
            # requirement that the fine grid lengths must be <= the corase
            # grid lengths.  FADD defaults to 20.
            #
            FADD = 20
            finedim = [min(coarsedim[i],box_length[i] + FADD) for i in range(3)]

            #
            # And now the hard part .. setting up the grid points.
            # From the APBS manual at http://agave.wustl.edu/apbs/doc/html/user-guide/x594.html#dime
            # we have the formula
            # n = c*2^(l+1) + 1
            # where l is the number of levels in the MG hierarchy.  The typical
            # number of levels is 4.
            #
            nlev = 4
            mult_fac = 2**(nlev + 1) # this will typically be 2^5==32
            # and c must be a non-zero integer

            # If we didn't have to be c*mult_fac + 1, this is what our grid points
            # would look like (we use the ceiling to be on the safe side .. it never
            # hurts to do too much.
            SPACE = 0.5 # default desired spacing = 0.5A
            #desired_points = [int(math.ceil(flen / SPACE)) for flen in finedim] # as integers
            desired_points = [flen / SPACE for flen in finedim] # as floats .. use int(math.ceil(..)) later

            # Now we set up our cs, taking into account mult_fac
            # (we use the ceiling to be on the safe side .. it never hurts to do
            # too much.)
            cs = [int(math.ceil(dp/mult_fac)) for dp in desired_points]

            finegridpoints = [mult_fac * c + 1 for c in cs]

            print "cs",cs
            print "finedim",finedim
            print "nlev",nlev
            print "mult_fac",mult_fac
            print "finegridpoints",finegridpoints

        except NoPDB:
            error_dialog = Pmw.MessageDialog(self.parent,
                                             title = 'Error',
                                             message_text = "Please set a temporary PDB file location",
                                             )
            junk = error_dialog.activate()
            return False

        if (finegridpoints[0]>0) and (finegridpoints[1]>0) and (finegridpoints[2]>0):
            max_grid_points = float(self.max_grid_points.getvalue())
            product = float(finegridpoints[0] * finegridpoints[1] * finegridpoints[2])
            if product>max_grid_points:
                print "Maximum number of grid points exceeded.  Old grid dimensions were",finegridpoints
                factor = pow(max_grid_points/product,0.333333333)
                finegridpoints[0] = (int(factor*finegridpoints[0]/2))*2+1
                finegridpoints[1] = (int(factor*finegridpoints[1]/2))*2+1
                finegridpoints[2] = (int(factor*finegridpoints[2]/2))*2+1
                print "Fine grid points rounded down from",finegridpoints
                #
                # Now we have to make sure that this still fits the equation n = c*2^(l+1) + 1.  Here, we'll
                # just assume nlev == 4, which means that we need to be (some constant times 32) + 1.
                #
                # This can be annoying if, e.g., you're trying to set [99, 123, 99] .. it'll get rounded to [99, 127, 99].
                # First, I'll try to round to the nearest 32*c+1.  If that doesn't work, I'll just round down.
                #
                new_gp = [0,0,0]
                for i in 0,1,2:
                    dm = divmod(finegridpoints[i] - 1,32)
                    if dm[1]>16:
                        new_gp[i] = (dm[0]+1)*32+1
                    else:
                        new_gp[i] = (dm[0])*32+1
                new_prod = new_gp[0]*new_gp[1]*new_gp[2]
                #print "tried new_prod",new_prod,"max_grid_points",max_grid_points,"small enough?",new_prod <= max_grid_points
                if new_prod <= max_grid_points:
                    #print "able to round to closest"
                    for i in 0,1,2: finegridpoints[i] = new_gp[i]
                else:
                    # darn .. have to round down.
                    # Note that this can still fail a little bit .. it can only get you back down to the next multiple <= what was in
                    # finegridpoints.  So, if finegridpoints was exactly on a multiple, like (99,129,99), you'll get rounded down to
                    # (99,127,99), which is still just a bit over the default max of 1200000.  I think that's ok.  It's the rounding error
                    # from int(factor*finegridpoints ..) above, but it'll never be a huge error.  If we needed to, we could easily fix this.
                    #
                    #print "rounding down more"
                    for i in 0,1,2:
                        #print finegridpoints[i],divmod(finegridpoints[i] - 1,32),
                        finegridpoints[i] = divmod(finegridpoints[i] - 1,32)[0]*32 + 1
                print "New grid dimensions are",finegridpoints
        print " APBS Tools: coarse grid: (%5.3f,%5.3f,%5.3f)"%tuple(coarsedim)
        self.grid_coarse_x.setvalue(coarsedim[0])
        self.grid_coarse_y.setvalue(coarsedim[1])
        self.grid_coarse_z.setvalue(coarsedim[2])
        print " APBS Tools: fine grid: (%5.3f,%5.3f,%5.3f)"%tuple(finedim)
        self.grid_fine_x.setvalue(finedim[0])
        self.grid_fine_y.setvalue(finedim[1])
        self.grid_fine_z.setvalue(finedim[2])
        print " APBS Tools: center: (%5.3f,%5.3f,%5.3f)"%tuple(center)
        self.grid_center_x.setvalue(center[0])
        self.grid_center_y.setvalue(center[1])
        self.grid_center_z.setvalue(center[2])
        print " APBS Tools: fine grid points (%d,%d,%d)"%tuple(finegridpoints)
        self.grid_points_x.setvalue(finegridpoints[0])
        self.grid_points_y.setvalue(finegridpoints[1])
        self.grid_points_z.setvalue(finegridpoints[2])
            
    def generateApbsInputFile(self):
        if self.checkInput(silent=True):
            #
            # set up our variables
            #
            pqr_filename = self.getPqrFilename()

            grid_points = [int(getattr(self,'grid_points_%s'%i).getvalue()) for i in 'x y z'.split()]
            cglen = [float(getattr(self,'grid_coarse_%s'%i).getvalue()) for i in 'x y z'.split()]
            fglen = [float(getattr(self,'grid_fine_%s'%i).getvalue()) for i in 'x y z'.split()]
            cent = [float(getattr(self,'grid_center_%s'%i).getvalue()) for i in 'x y z'.split()]

            apbs_mode = self.apbs_mode.getvalue()
            if apbs_mode == 'Nonlinear Poisson-Boltzmann Equation':
                apbs_mode = 'npbe'
            else:
                apbs_mode = 'lpbe'

            bcflmap = {'Zero': 'zero',
                       'Single DH sphere': 'sdh',
                       'Multiple DH spheres': 'mdh',
                       'Focusing': 'focus',
                       }
            bcfl = bcflmap[self.bcfl.getvalue()]

            chgmmap = {'Linear':'spl0',
                       'Cubic b-splines':'spl2',
                       }
            chgm = chgmmap[self.chgm.getvalue()]

            dx_filename = self.pymol_generated_dx_filename.getvalue()
            if dx_filename.endswith('.dx'):
                dx_filename = dx_filename[:-3]

            #
            # get the input text
            #

            apbs_input_text = _getApbsInputFile(pqr_filename,
                                                grid_points,
                                                cglen,
                                                fglen,
                                                cent,
                                                apbs_mode,
                                                bcfl,
                                                float(self.ion_plus_one_conc.getvalue()), float(self.ion_plus_one_rad.getvalue()),
                                                float(self.ion_minus_one_conc.getvalue()), float(self.ion_minus_one_rad.getvalue()),
                                                float(self.ion_plus_two_conc.getvalue()), float(self.ion_plus_two_rad.getvalue()),
                                                float(self.ion_minus_two_conc.getvalue()), float(self.ion_minus_two_rad.getvalue()),
                                                float(self.interior_dielectric.getvalue()),
                                                float(self.solvent_dielectric.getvalue()),
                                                chgm,
                                                float(self.solvent_radius.getvalue()),
                                                float(self.system_temp.getvalue()),
                                                float(self.sdens.getvalue()),
                                                dx_filename,
                                                )
                                                
                                                

            #
            # write out the input text
            #
            f = file(self.pymol_generated_in_filename.getvalue(),'w')
            f.write(apbs_input_text)
            f.close()
            return True
        else:
            self.checkInput()
            return False

    def checkInput(self,silent=False):
        """If silent is True, we'll just return a True/False value
        """
        if not silent:
            def show_error(message):
                error_dialog = Pmw.MessageDialog(self.parent,
                                                 title = 'Error',
                                                 message_text = message,
                                                 )
                junk = error_dialog.activate()
        else:
            def show_error(message):
                pass
            
        #
        # First, check to make sure we have valid locations for apbs and psize
        #
        if not self.binary.valid():
            show_error('Please set the APBS binary location')
            return False
        #
        # If the path to psize is not correct, that's fine .. we'll
        # do the calculations ourself.
        #
        
        #if not self.psize.valid():
        #    show_error("Please set APBS's psize location")
        #    return False
        
        #
        # Now check the temporary filenames
        #
        if self.radiobuttons.getvalue() != 'Use another PQR':
            if not self.pymol_generated_pqr_filename.getvalue():
                show_error('Please choose a name for the PyMOL\ngenerated PQR file')
                return False
        elif not self.pqr_to_use.valid():
            show_error('Please select a valid pqr file or tell\nPyMOL to generate one')
            return False
        if not self.pymol_generated_pdb_filename.getvalue():
            show_error('Please choose a name for the PyMOL\ngenerated PDB file')
            return False
        if not self.pymol_generated_dx_filename.getvalue():
            show_error('Please choose a name for the PyMOL\ngenerated DX file')
            return False
        if not self.pymol_generated_in_filename.getvalue():
            show_error('Please choose a name for the PyMOL\ngenerated APBS input file')
            return False
        if not self.map.getvalue():
            show_error('Please choose a name for the generated map.')
            return False
        
        
        #
        # Now, the ions
        #
        for sign in 'plus minus'.split():
            for value in 'one two'.split():
                for parm in 'conc rad'.split():
                    if not getattr(self,'ion_%s_%s_%s'%(sign,value,parm)).valid():
                        show_error('Please correct Ion concentrations and radii')
                        return False
        #
        # Now the grid
        #
        for grid_type in 'coarse fine points center'.split():
            for coord in 'x y z'.split():
                if not getattr(self,'grid_%s_%s'%(grid_type,coord)).valid():
                    show_error('Please correct grid dimensions\nby clicking on the "Set grid" button')
                    return False
        
        #
        # Now other easy things
        #
        for (message, thing) in (('solvent dielectric',self.solvent_dielectric),
                                 ('protein dielectric',self.interior_dielectric),
                                 ('solvent radius',self.solvent_radius),
                                 ('system temperature',self.system_temp),
                                 ('sdens',self.sdens),
                                 ):
            if not thing.valid():
                show_error('Please correct %s'%message)
                return False

        return True
        
        
    def generatePqrFile(self):
        """generate a pqr file

        This will also call through to champ to set the Hydrogens and charges
        if it needs to.  If it does that, it may change the value self.selection
        to take the new Hydrogens into account.

        To make it worse, APBS seems to freak out when there are chain ids.  So,
        this gets rid of the chain ids.
        """
        # WLD
        sel = "((%s) or (neighbor (%s) and hydro))"%(
            self.selection.getvalue(), self.selection.getvalue())
        
        pqr_filename = self.getPqrFilename()
        try:
            f = file(pqr_filename,'w')
            f.close()
        except:
            error_dialog = Pmw.MessageDialog(self.parent,
                                             title = 'Error',
                                             message_text = "Could not write PQR file.\nPlease check that temporary PQR filename is valid.",
                                             )
            junk = error_dialog.activate()
            return False
            
        if REAL_PYMOL:
            # PyMOL + champ == pqr
            from chempy.champ import assign
            if self.radiobuttons.getvalue() == 'Use PyMOL generated PQR and PyMOL generated Hydrogens and termini':
                pymol.cmd.remove('hydro and %s'%sel)
                assign.missing_c_termini(sel)
                assign.formal_charges(sel)
                pymol.cmd.h_add(sel)
# WLD (code now unnecessary)
#                new_hydros = '(hydro and neighbor %s)'%sel
#                sel = '%s or %s' % (sel,new_hydros)
            assign.amber99(sel)
# WLD (code now unnecessary)
#            if not self.selection.getvalue() in '(all) all'.split():
#                self.selection.setvalue(sel)
                
        #
        # Get rid of chain information
        #
# WLD -- PyMOL now does this automatically with PQR files        
#        pymol.cmd.alter(sel,'chain = ""')
        pymol.cmd.save(pqr_filename,sel)
        missed_count = pymol.cmd.count_atoms("("+sel+") and flag 23")
        if missed_count > 0:
            pymol.cmd.select("unassigned","("+sel+") and flag 23")
            error_dialog = Pmw.MessageDialog(self.parent,
                                             title = 'Error',
                                             message_text = "Unable to assign parameters for the %s atoms in selection 'unassigned'.\nPlease either remove these unassigned atoms and re-start the calculation\nor fix their parameters in the generated PQR file and run the calculation\nusing the modified PQR file (select 'Use another PQR' in 'Main')."%missed_count,
                                             )
            junk = error_dialog.activate()
            return False
        return True


#
# The classes PmwFileDialog and PmwExistingFileDialog and the _errorpop function
# are taken from the Pmw contrib directory.  The attribution given in that file
# is:
################################################################################
# Filename dialogs using Pmw
#
# (C) Rob W.W. Hooft, Nonius BV, 1998
#
# Modifications:
#
# J. Willem M. Nissink, Cambridge Crystallographic Data Centre, 8/2002
#    Added optional information pane at top of dialog; if option
#    'info' is specified, the text given will be shown (in blue).
#    Modified example to show both file and directory-type dialog
#
# No Guarantees. Distribute Freely. 
# Please send bug-fixes/patches/features to <r.hooft@euromail.com>
#
################################################################################
import os,fnmatch,time
import Tkinter,Pmw
#Pmw.setversion("0.8.5")

def _errorpop(master,text):
    d=Pmw.MessageDialog(master,
                        title="Error", 
                        message_text=text,
                        buttons=("OK",))
    d.component('message').pack(ipadx=15,ipady=15)
    d.activate()
    d.destroy()
    
class PmwFileDialog(Pmw.Dialog):
    """File Dialog using Pmw"""
    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	optiondefs = (
	    ('filter',    '*',              self.newfilter),
	    ('directory', os.getcwd(),      self.newdir),
	    ('filename',  '',               self.newfilename),
	    ('historylen',10,               None),
	    ('command',   None,             None),
            ('info',      None,             None),
	    )
	self.defineoptions(kw, optiondefs)
        # Initialise base class (after defining options).
	Pmw.Dialog.__init__(self, parent)

	self.withdraw()

        # Create the components.
	interior = self.interior()

        if self['info'] is not None:
            rowoffset=1
            dn = self.infotxt()
            dn.grid(row=0,column=0,columnspan=2,padx=3,pady=3)
        else:
            rowoffset=0

	dn = self.mkdn()
	dn.grid(row=0+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	del dn

	# Create the directory list component.
	dnb = self.mkdnb()
	dnb.grid(row=1+rowoffset,column=0,sticky='news',padx=3,pady=3)
	del dnb

	# Create the filename list component.
	fnb = self.mkfnb()
	fnb.grid(row=1+rowoffset,column=1,sticky='news',padx=3,pady=3)
	del fnb

	# Create the filter entry
	ft = self.mkft()
	ft.grid(row=2+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	del ft

	# Create the filename entry
	fn = self.mkfn()
	fn.grid(row=3+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	fn.bind('<Return>',self.okbutton)
	del fn

	# Buttonbox already exists
	bb=self.component('buttonbox')
	bb.add('OK',command=self.okbutton)
	bb.add('Cancel',command=self.cancelbutton)
	del bb

	Pmw.alignlabels([self.component('filename'),
			 self.component('filter'),
			 self.component('dirname')])

    def infotxt(self):
        """ Make information block component at the top """
        return self.createcomponent(
                'infobox',
                (), None,
                Tkinter.Label, (self.interior(),),
                width=51,
                relief='groove',
                foreground='darkblue',
                justify='left',
                text=self['info']
            )

    def mkdn(self):
        """Make directory name component"""
        return self.createcomponent(
	    'dirname',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['directory'],
	    entryfield_entry_width=40,
            entryfield_validate=self.dirvalidate,
	    selectioncommand=self.setdir,
	    labelpos='w',
	    label_text='Directory:')

    def mkdnb(self):
        """Make directory name box"""
        return self.createcomponent(
	    'dirnamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='directories',
	    labelpos='n',
	    hscrollmode='none',
	    dblclickcommand=self.selectdir)

    def mkft(self):
        """Make filter"""
        return self.createcomponent(
	    'filter',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filter'],
	    entryfield_entry_width=40,
	    selectioncommand=self.setfilter,
	    labelpos='w',
	    label_text='Filter:')

    def mkfnb(self):
        """Make filename list box"""
        return self.createcomponent(
	    'filenamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='files',
	    labelpos='n',
	    hscrollmode='none',
	    selectioncommand=self.singleselectfile,
	    dblclickcommand=self.selectfile)

    def mkfn(self):
        """Make file name entry"""
        return self.createcomponent(
	    'filename',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filename'],
	    entryfield_entry_width=40,
            entryfield_validate=self.filevalidate,
	    selectioncommand=self.setfilename,
	    labelpos='w',
	    label_text='Filename:')
    
    def dirvalidate(self,string):
        if os.path.isdir(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL
        
    def filevalidate(self,string):
        if string=='':
            return Pmw.PARTIAL
        elif os.path.isfile(string):
            return Pmw.OK
        elif os.path.exists(string):
            return Pmw.PARTIAL
        else:
            return Pmw.OK
        
    def okbutton(self):
	"""OK action: user thinks he has input valid data and wants to
           proceed. This is also called by <Return> in the filename entry"""
	fn=self.component('filename').get()
	self.setfilename(fn)
	if self.validate(fn):
	    self.canceled=0
	    self.deactivate()

    def cancelbutton(self):
	"""Cancel the operation"""
	self.canceled=1
	self.deactivate()

    def tidy(self,w,v):
	"""Insert text v into the entry and at the top of the list of 
           the combobox w, remove duplicates"""
	if not v:
	    return
	entry=w.component('entry')
	entry.delete(0,'end')
	entry.insert(0,v)
	list=w.component('scrolledlist')
	list.insert(0,v)
	index=1
	while index<list.index('end'):
	    k=list.get(index)
	    if k==v or index>self['historylen']:
		list.delete(index)
	    else:
		index=index+1
        w.checkentry()

    def setfilename(self,value):
	if not value:
	    return
	value=os.path.join(self['directory'],value)
	dir,fil=os.path.split(value)
	self.configure(directory=dir,filename=value)
        
	c=self['command']
	if callable(c):
	    c()

    def newfilename(self):
	"""Make sure a newly set filename makes it into the combobox list"""
	self.tidy(self.component('filename'),self['filename'])
	
    def setfilter(self,value):
	self.configure(filter=value)

    def newfilter(self):
	"""Make sure a newly set filter makes it into the combobox list"""
	self.tidy(self.component('filter'),self['filter'])
	self.fillit()

    def setdir(self,value):
	self.configure(directory=value)

    def newdir(self):
	"""Make sure a newly set dirname makes it into the combobox list"""
	self.tidy(self.component('dirname'),self['directory'])
	self.fillit()

    def singleselectfile(self):
	"""Single click in file listbox. Move file to "filename" combobox"""
	cs=self.component('filenamebox').curselection()
	if cs!=():
	    value=self.component('filenamebox').get(cs)
            self.setfilename(value)

    def selectfile(self):
	"""Take the selected file from the filename, normalize it, and OK"""
        self.singleselectfile()
	value=self.component('filename').get()
        self.setfilename(value)
        if value:
	    self.okbutton()

    def selectdir(self):
	"""Take selected directory from the dirnamebox into the dirname"""
	cs=self.component('dirnamebox').curselection()
	if cs!=():
	    value=self.component('dirnamebox').get(cs)
	    dir=self['directory']
	    if not dir:
		dir=os.getcwd()
	    if value:
		if value=='..':
		    dir=os.path.split(dir)[0]
		else:
		    dir=os.path.join(dir,value)
	    self.configure(directory=dir)
	    self.fillit()

    def askfilename(self,directory=None,filter=None):
	"""The actual client function. Activates the dialog, and
	   returns only after a valid filename has been entered 
           (return value is that filename) or when canceled (return 
           value is None)"""
	if directory!=None:
	    self.configure(directory=directory)
	if filter!=None:
	    self.configure(filter=filter)
	self.fillit()
        self.canceled=1 # Needed for when user kills dialog window
	self.activate()
	if self.canceled:
	    return None
	else:
	    return self.component('filename').get()

    lastdir=""
    lastfilter=None
    lasttime=0
    def fillit(self):
	"""Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
        if self.lastdir==self['directory'] and self.lastfilter==self['filter'] and self.lasttime>os.stat(self.lastdir)[8]:
            return
        self.lastdir=self['directory']
        self.lastfilter=self['filter']
        self.lasttime=time.time()
	dir=self['directory']
	if not dir:
	    dir=os.getcwd()
	dirs=['..']
	files=[]
        try:
            fl=os.listdir(dir)
            fl.sort()
        except os.error,arg:
            if arg[0] in (2,20):
                return
            raise
	for f in fl:
	    if os.path.isdir(os.path.join(dir,f)):
		dirs.append(f)
	    else:
		filter=self['filter']
		if not filter:
		    filter='*'
		if fnmatch.fnmatch(f,filter):
		    files.append(f)
	self.component('filenamebox').setlist(files)
	self.component('dirnamebox').setlist(dirs)
    
    def validate(self,filename):
	"""Validation function. Should return 1 if the filename is valid, 
           0 if invalid. May pop up dialogs to tell user why. Especially 
           suited to subclasses: i.e. only return 1 if the file does/doesn't 
           exist"""
	return 1

class PmwExistingFileDialog(PmwFileDialog):
    def filevalidate(self,string):
        if os.path.isfile(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL
        
    def validate(self,filename):
        if os.path.isfile(filename):
            return 1
        elif os.path.exists(filename):
            _errorpop(self.interior(),"This is not a plain file")
            return 0
        else:
            _errorpop(self.interior(),"Please select an existing file")
            return 0

##############################################################################
##############################################################################
###                                                                        ###
###       Now the actual APBS code                                         ###
###                                                                        ###
##############################################################################
##############################################################################

def _getApbsInputFile(pqr_filename,
                      grid_points,
                      cglen,
                      fglen,
                      cent,
                      apbs_mode,
                      bcfl,
                      ion_plus_one_conc,ion_plus_one_rad,
                      ion_minus_one_conc,ion_minus_one_rad,
                      ion_plus_two_conc,ion_plus_two_rad,
                      ion_minus_two_conc,ion_minus_two_rad,
                      interior_dielectric, solvent_dielectric,
                      chgm,
                      solvent_radius,
                      system_temp,
                      sdens,
                      dx_filename,
                      ):
    #print "system_temp",system_temp,type(system_temp)
    #print "sdens",sdens,type(sdens)
    #
    # How shall we set up the grid?  We'll use cglen, fglen, cgcent, fgcent
    # and dime.
    # This allows people to automate things (e.g. "Alanine scanning")
    #
    
    #
    # New template using mg-auto
    # See http://agave.wustl.edu/apbs/doc/api/html/#mg-auto
    #

    apbs_template = """#
# Note that most of the comments here were taken from sample
# input files that came with APBS.  You can find APBS at
# http://agave.wustl.edu/apbs/
# Note that APBS is GPL'd code.
#
read
    mol pqr %s       # read molecule 1
end
elec
    mg-auto
    dime   %d %d %d  # number of find grid points
                     # calculated by psize.py
    cglen   %f %f %f # coarse mesh lengths (A)
    fglen   %f %f %f # fine mesh lengths (A)
                     # calculated by psize.py
    cgcent %f %f %f  # (could also give (x,y,z) form psize.py) #known center
    fgcent %f %f %f  # (could also give (x,y,z) form psize.py) #known center
    %s               # solve the full nonlinear PBE with npbe
    #lpbe            # solve the linear PBE with lpbe
    bcfl %s          # Boundary condition flag
                     #  0 => Zero
                     #  1 => Single DH sphere
                     #  2 => Multiple DH spheres
                     #  4 => Focusing
                     #
    #ion 1 0.000 2.0 # Counterion declaration:
    ion  1 %f %f     # Counterion declaration:
    ion -1 %f %f     # ion <charge> <conc (M)> <radius>
    ion  2 %f %f     # ion <charge> <conc (M)> <radius>
    ion -2 %f %f     # ion <charge> <conc (M)> <radius>
    pdie %f          # Solute dielectric
    sdie %f          # Solvent dielectric
    chgm %s          # Charge disc method
                     # 0 is linear splines
                     # 1 is cubic b-splines
    mol 1            # which molecule to use
    srfm smol        # Surface calculation method
                     #  0 => Mol surface for epsilon;
                     #       inflated VdW for kappa; no
                     #       smoothing
                     #  1 => As 0 with harmoic average
                     #       smoothing
                     #  2 => Cubic spline 
    srad %f          # Solvent radius (1.4 for water)
    swin 0.3         # Surface cubic spline window .. default 0.3
    temp %f          # System temperature (298.15 default)
    sdens %f         # Specify the number of grid points per square-angstrom to use in Vacc object. Ignored when srad is 0.0 (see srad) or srfm is spl2 (see srfm). There is a direct correlation between the value used for the Vacc sphere density, the accuracy of the Vacc object, and the APBS calculation time. APBS default value is 10.0.
    gamma 0.105      # Surface tension parameter for apolar forces (in kJ/mol/A^2)
                     # only used for force calculations, so we don't care, but
                     # it's always required, and 0.105 is the default
    calcenergy no    # Energy I/O to stdout
                     #  0 => don't write out energy
                     #  1 => write out total energy
                     #  2 => write out total energy and all
                     #       components
    calcforce no     # Atomic forces I/O (to stdout)
                     #  0 => don't write out forces
                     #  1 => write out net forces on molecule
                     #  2 => write out atom-level forces
    write pot dx %s  # What to write .. this says write the potential in dx
                     # format to a file.
end
quit

"""
    return apbs_template % (pqr_filename,
                            grid_points[0], grid_points[1], grid_points[2],
                            cglen[0],cglen[1],cglen[2],
                            fglen[0],fglen[1],fglen[2],
                            cent[0],cent[1],cent[2],
                            cent[0],cent[1],cent[2],
                            apbs_mode,
                            bcfl,
                            ion_plus_one_conc,ion_plus_one_rad,
                            ion_minus_one_conc,ion_minus_one_rad,
                            ion_plus_two_conc,ion_plus_two_rad,
                            ion_minus_two_conc,ion_minus_two_rad,
                            interior_dielectric, solvent_dielectric,
                            chgm,
                            solvent_radius,
                            system_temp,
                            sdens,
                            dx_filename,
                            )                            
                            





# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')
    
    widget = APBSTools(app)
    exitButton = Tkinter.Button(app.root, text = 'Exit', command = app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
