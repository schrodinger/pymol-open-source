# -c

/print "BEGIN-LOG"

# HEADER without END

/pdbstr="""HEADER PTPBR7 \
ATOM      1  N   GLY A 252      39.368  46.327  27.875  1.00 46.28           N  \
ATOM      2  CA  GLY A 252      39.565  46.477  29.344  1.00 46.89           C  \
ATOM      3  C   GLY A 252      40.602  45.521  29.907  1.00 48.55           C  \
ATOM      4  O   GLY A 252      40.907  44.497  29.293  1.00 46.16           O  \
ATOM      5  N   SER A 253      41.144  45.862  31.078  1.00 45.72           N  \
ATOM      6  CA  SER A 253      42.152  45.043  31.757  1.00 46.42           C  \
ATOM      7  C   SER A 253      41.824  44.945  33.254  1.00 47.72           C  \
ATOM      8  O   SER A 253      41.053  45.756  33.778  1.00 46.62           O  \
ATOM      9  CB  SER A 253      43.543  45.667  31.587  1.00 45.55           C  \
HEADER STEP\
ATOM      1  N   SER E 258      55.863  54.272  24.437  1.00 37.51           N  \
ATOM      2  CA  SER E 258      56.590  54.677  23.240  1.00 36.33           C  \
ATOM      3  C   SER E 258      58.022  54.180  23.260  1.00 35.14           C  \
ATOM      4  O   SER E 258      58.647  54.090  24.306  1.00 34.89           O  \
ATOM      5  CB  SER E 258      56.566  56.199  23.083  1.00 35.95           C  \
ATOM      6  OG  SER E 258      57.340  56.607  21.975  1.00 36.25           O  \
ATOM      7  N   ARG E 259      58.546  53.842  22.088  1.00 34.24           N  """

cmd.read_pdbstr(pdbstr,"test")

print cmd.get_names()
count_atoms
reinit

# HEADER with END

cmd.read_pdbstr("""\
HEADER    OXIDOREDUCTASE                          13-JUL-94   3MDD      3MDD   2\
ATOM      1  N   GLY A  11      22.422  32.842 -15.744  1.00 24.74      3MDD 172\
ATOM      2  CA  GLY A  11      21.837  32.022 -14.672  1.00 24.52      3MDD 173\
ATOM      3  C   GLY A  11      21.575  30.609 -15.194  1.00 24.62      3MDD 174\
ATOM      4  O   GLY A  11      22.036  30.238 -16.285  1.00 24.99      3MDD 175\
ATOM      5  N   PHE A  12      20.843  29.864 -14.390  1.00 24.29      3MDD 176\
ATOM      6  CA  PHE A  12      20.481  28.469 -14.679  1.00 24.00      3MDD 177\
ATOM      7  C   PHE A  12      19.514  28.397 -15.861  1.00 23.38      3MDD 178\
ATOM   6017  NZ  LYS B 395      -4.299  -6.687  15.377  1.00 35.83      3MDD6188\
ATOM   6018  OXT LYS B 395      -0.899 -11.439  20.683  1.00 35.32      3MDD6189\
END                                                                     3MDD6517\
HEADER    ELECTRON TRANSPORT                      02-MAR-92   1MDA      1MDA   2\
ATOM      1  N   GLU H   1      15.906  81.036  72.025  1.00 50.00      1MDA 307\
ATOM      2  CA  GLU H   1      17.350  81.236  71.880  1.00 50.00      1MDA 308\
ATOM      3  C   GLU H   1      17.854  81.695  70.506  1.00 50.00      1MDA 309\
ATOM      4  O   GLU H   1      17.668  80.983  69.518  1.00 50.00      1MDA 310\
ATOM      5  CB  GLU H   1      18.222  80.229  72.680  1.00 50.00      1MDA 311\
ATOM      6  CG  GLU H   1      17.540  78.863  72.954  1.00 50.00      1MDA 312\
ATOM      7  CD  GLU H   1      17.831  77.718  72.017  1.00 43.59      1MDA 313\
ATOM      8  OE1 GLU H   1      18.953  77.453  71.614  1.00 47.42      1MDA 314\
ATOM      9  OE2 GLU H   1      16.759  77.049  71.674  1.00 37.68      1MDA 315\
ATOM   8568  CG  GLU B 105      37.556  38.617  32.271  1.00 15.99      1MDA8874\
ATOM   8569  CD  GLU B 105      36.485  37.572  32.282  1.00 36.23      1MDA8875\
ATOM   8570  OE1 GLU B 105      36.488  36.782  31.231  1.00 50.00      1MDA8876\
ATOM   8571  OE2 GLU B 105      35.647  37.508  33.152  1.00 50.00      1MDA8877\
ATOM   8572  OXT GLU B 105      41.453  38.354  30.968  1.00 50.00      1MDA8878\
END                                                                     1MDA8924\
""","test")

print cmd.get_names()
count_atoms

reinit

# just END...

cmd.read_pdbstr("""\
ATOM      1  N   GLY A  11      22.422  32.842 -15.744  1.00 24.74      3MDD 172\
ATOM      2  CA  GLY A  11      21.837  32.022 -14.672  1.00 24.52      3MDD 173\
ATOM      3  C   GLY A  11      21.575  30.609 -15.194  1.00 24.62      3MDD 174\
ATOM      4  O   GLY A  11      22.036  30.238 -16.285  1.00 24.99      3MDD 175\
ATOM      5  N   PHE A  12      20.843  29.864 -14.390  1.00 24.29      3MDD 176\
ATOM      6  CA  PHE A  12      20.481  28.469 -14.679  1.00 24.00      3MDD 177\
ATOM      7  C   PHE A  12      19.514  28.397 -15.861  1.00 23.38      3MDD 178\
ATOM   6017  NZ  LYS B 395      -4.299  -6.687  15.377  1.00 35.83      3MDD6188\
ATOM   6018  OXT LYS B 395      -0.899 -11.439  20.683  1.00 35.32      3MDD6189\
END                                                                     3MDD6517\
ATOM      1  N   GLU H   1      15.906  81.036  72.025  1.00 50.00      1MDA 307\
ATOM      2  CA  GLU H   1      17.350  81.236  71.880  1.00 50.00      1MDA 308\
ATOM      3  C   GLU H   1      17.854  81.695  70.506  1.00 50.00      1MDA 309\
ATOM      4  O   GLU H   1      17.668  80.983  69.518  1.00 50.00      1MDA 310\
ATOM      5  CB  GLU H   1      18.222  80.229  72.680  1.00 50.00      1MDA 311\
ATOM      6  CG  GLU H   1      17.540  78.863  72.954  1.00 50.00      1MDA 312\
ATOM      7  CD  GLU H   1      17.831  77.718  72.017  1.00 43.59      1MDA 313\
ATOM      8  OE1 GLU H   1      18.953  77.453  71.614  1.00 47.42      1MDA 314\
ATOM      9  OE2 GLU H   1      16.759  77.049  71.674  1.00 37.68      1MDA 315\
ATOM   8568  CG  GLU B 105      37.556  38.617  32.271  1.00 15.99      1MDA8874\
ATOM   8569  CD  GLU B 105      36.485  37.572  32.282  1.00 36.23      1MDA8875\
ATOM   8570  OE1 GLU B 105      36.488  36.782  31.231  1.00 50.00      1MDA8876\
ATOM   8571  OE2 GLU B 105      35.647  37.508  33.152  1.00 50.00      1MDA8877\
ATOM   8572  OXT GLU B 105      41.453  38.354  30.968  1.00 50.00      1MDA8878\
END                                                                     1MDA8924\
""","test")

print cmd.get_names()
count_atoms

reinit

/print "END-LOG"

