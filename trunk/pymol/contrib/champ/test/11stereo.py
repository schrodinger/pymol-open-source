from chempy.champ import Champ

ch = Champ()

# basic validation of stereochemistry handling on the simplest system

IS = [ # implicit patterns (S) 
   "[C@H](C)(N)O",
   "[C@H](N)(O)C",
   "[C@H](O)(C)N",

   "C[C@H](N)O",
   "N[C@H](O)C",
   "O[C@H](C)N",

   "[C@@H](C)(O)N",
   "[C@@H](N)(C)O",
   "[C@@H](O)(N)C",

   "C[C@@H](O)N",
   "N[C@@H](C)O",
   "O[C@@H](N)C",

   ]

IR = [ # implicit patterns (R) 
   "[C@@H](C)(N)O",
   "[C@@H](N)(O)C",
   "[C@@H](O)(C)N",

   "C[C@@H](N)O",
   "N[C@@H](O)C",
   "O[C@@H](C)N",

   "[C@H](C)(O)N",
   "[C@H](N)(C)O",
   "[C@H](O)(N)C",

   "C[C@H](O)N",
   "N[C@H](C)O",
   "O[C@H](N)C",

   ]

XS = [ # explicit patterns (S) 

   "H[C@](C)(N)O",
   "H[C@](N)(O)C",
   "H[C@](O)(C)N",

   "C[C@](H)(O)N",
   "N[C@](H)(C)O",
   "O[C@](H)(N)C",

   "O[C@](N)(C)H",
   "C[C@](O)(N)H",
   "N[C@](C)(O)H",

   "N[C@](O)(H)C",
   "O[C@](C)(H)N",
   "C[C@](N)(H)O",

   "H[C@@](C)(O)N",
   "H[C@@](N)(C)O",
   "H[C@@](O)(N)C",

   "C[C@@](H)(N)O",
   "N[C@@](H)(O)C",
   "O[C@@](H)(C)N",

   "O[C@@](N)(H)C",
   "C[C@@](O)(H)N",
   "N[C@@](C)(H)O",

   "N[C@@](O)(C)H",
   "O[C@@](C)(N)H",
   "C[C@@](N)(O)H",

   ]

XR = [ # explicit patterns (R) 

   "H[C@@](C)(N)O",
   "H[C@@](N)(O)C",
   "H[C@@](O)(C)N",

   "C[C@@](H)(O)N",
   "N[C@@](H)(C)O",
   "O[C@@](H)(N)C",

   "O[C@@](N)(C)H",
   "C[C@@](O)(N)H",
   "N[C@@](C)(O)H",

   "N[C@@](O)(H)C",
   "O[C@@](C)(H)N",
   "C[C@@](N)(H)O",

   "H[C@](C)(O)N",
   "H[C@](N)(C)O",
   "H[C@](O)(N)C",

   "C[C@](H)(N)O",
   "N[C@](H)(O)C",
   "O[C@](H)(C)N",

   "O[C@](N)(H)C",
   "C[C@](O)(H)N",
   "N[C@](C)(H)O",

   "N[C@](O)(C)H",
   "O[C@](C)(N)H",
   "C[C@](N)(O)H",

   ]

I = [ # implicit patterns (achiral)

   "[CH](C)(N)O",
   "[CH](N)(O)C",
   "[CH](O)(C)N",

   "C[CH](N)O",
   "N[CH](O)C",
   "O[CH](C)N",

   ]

X = [ # explicit patterns (achiral)

   "HC(C)(N)O",
   "HC(N)(O)C",
   "HC(O)(C)N",

   "CC(H)(O)N",
   "NC(H)(C)O",
   "OC(H)(N)C",

   "OC(N)(C)H",
   "CC(O)(N)H",
   "NC(C)(O)H",

   "NC(O)(H)C",
   "OC(C)(H)N",
   "CC(N)(H)O",

   ]

# print 

PIS = map(lambda x:ch.insert_pattern_string(x),IS)
PXS = map(lambda x:ch.insert_pattern_string(x),XS)

for a in range(0,len(PIS)):
   print IS[a],ch.pattern_get_string(PIS[a])

for a in range(0,len(PXS)):
   print XS[a],ch.pattern_get_string(PXS[a])

PIR = map(lambda x:ch.insert_pattern_string(x),IR)
PXR = map(lambda x:ch.insert_pattern_string(x),XR)

for a in range(0,len(PIR)):
   print IR[a],ch.pattern_get_string(PIR[a])

for a in range(0,len(PXR)):
   print XR[a],ch.pattern_get_string(PXR[a])


PI = map(lambda x:ch.insert_pattern_string(x),I)
PX = map(lambda x:ch.insert_pattern_string(x),X)

for a in range(0,len(PI)):
   print I[a],ch.pattern_get_string(PI[a])

for a in range(0,len(PX)):
   print X[a],ch.pattern_get_string(PX[a])


