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

# S match tests

PIS = map(lambda x:ch.insert_pattern_string(x),IS)
TIS = map(lambda x:ch.insert_pattern_string(x),IS)

PXS = map(lambda x:ch.insert_pattern_string(x),XS)
TXS = map(lambda x:ch.insert_pattern_string(x),XS)

for a in range(0,len(PIS)):
   for b in range(0,len(TIS)):
      if not (ch.match_1v1_b(PIS[a],TIS[b]) == 1):
         print "Error: PIS[%d],TIS[%d]"%(a,b)

for a in range(0,len(TIS)):
   for b in range(0,len(TIS)):
      if not (ch.match_1v1_b(TIS[a],TIS[b]) == 1):
         print "Error: TIS[%d],TIS[%d]"%(a,b)

for a in range(0,len(PXS)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(PXS[a],TXS[b]) == 1):
         print "Error: PXS[%d],TXS[%d]"%(a,b)

for a in range(0,len(TXS)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(TXS[a],TXS[b]) == 1):
         print "Error: TXS[%d],TXS[%d]"%(a,b)

for a in range(0,len(PIS)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(PIS[a],TXS[b]) == 1):
         print "Error: PIS[%d],TXS[%d]"%(a,b)

# R match tests

PIR = map(lambda x:ch.insert_pattern_string(x),IR)
TIR = map(lambda x:ch.insert_pattern_string(x),IR)

PXR = map(lambda x:ch.insert_pattern_string(x),XR)
TXR = map(lambda x:ch.insert_pattern_string(x),XR)

for a in range(0,len(PIR)):
   for b in range(0,len(TIR)):
      if not (ch.match_1v1_b(PIR[a],TIR[b]) == 1):
         print "Error: PIR[%d],TIR[%d]"%(a,b)

for a in range(0,len(TIR)):
   for b in range(0,len(TIR)):
      if not (ch.match_1v1_b(TIR[a],TIR[b]) == 1):
         print "Error: TIR[%d],TIR[%d]"%(a,b)

for a in range(0,len(PXR)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(PXR[a],TXR[b]) == 1):
         print "Error: PXR[%d],TXR[%d]"%(a,b)

for a in range(0,len(TXR)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(TXR[a],TXR[b]) == 1):
         print "Error: TXR[%d],TXR[%d]"%(a,b)

for a in range(0,len(PIR)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(PIR[a],TXR[b]) == 1):
         print "Error: PIR[%d],TXR[%d]"%(a,b)

# achiral->achiral match tests

PI = map(lambda x:ch.insert_pattern_string(x),I)
TI = map(lambda x:ch.insert_pattern_string(x),I)

PX = map(lambda x:ch.insert_pattern_string(x),X)
TX = map(lambda x:ch.insert_pattern_string(x),X)

for a in range(0,len(PI)):
   for b in range(0,len(TI)):
      if not (ch.match_1v1_b(PI[a],TI[b]) == 1):
         print "Error: PI[%d],TI[%d]"%(a,b)

for a in range(0,len(TI)):
   for b in range(0,len(TI)):
      if not (ch.match_1v1_b(TI[a],TI[b]) == 1):
         print "Error: TI[%d],TI[%d]"%(a,b)

for a in range(0,len(PX)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(PX[a],TX[b]) == 1):
         print "Error: PX[%d],TX[%d]"%(a,b)

for a in range(0,len(TX)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(TX[a],TX[b]) == 1):
         print "Error: TX[%d],TX[%d]"%(a,b)

for a in range(0,len(PI)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(PI[a],TX[b]) == 1):
         print "Error: PI[%d],TX[%d]"%(a,b)

# achiral->chiral match tests

for a in range(0,len(PI)):
   for b in range(0,len(TIS)):
      if not (ch.match_1v1_b(PI[a],TIS[b]) == 1):
         print "Error: PI[%d],TIS[%d]"%(a,b)

for a in range(0,len(TI)):
   for b in range(0,len(TIS)):
      if not (ch.match_1v1_b(TI[a],TIS[b]) == 1):
         print "Error: TI[%d],TIS[%d]"%(a,b)

for a in range(0,len(PX)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(PX[a],TXS[b]) == 1):
         print "Error: PX[%d],TXS[%d]"%(a,b)

for a in range(0,len(TX)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(TX[a],TXS[b]) == 1):
         print "Error: TX[%d],TXS[%d]"%(a,b)

for a in range(0,len(PI)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(PI[a],TXS[b]) == 1):
         print "Error: PI[%d],TXS[%d]"%(a,b)

for a in range(0,len(PI)):
   for b in range(0,len(TIR)):
      if not (ch.match_1v1_b(PI[a],TIR[b]) == 1):
         print "Error: PI[%d],TIR[%d]"%(a,b)

for a in range(0,len(TI)):
   for b in range(0,len(TIR)):
      if not (ch.match_1v1_b(TI[a],TIR[b]) == 1):
         print "Error: TI[%d],TIR[%d]"%(a,b)

for a in range(0,len(PX)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(PX[a],TXR[b]) == 1):
         print "Error: PX[%d],TXR[%d]"%(a,b)

for a in range(0,len(TX)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(TX[a],TXR[b]) == 1):
         print "Error: TX[%d],TXR[%d]"%(a,b)

for a in range(0,len(PI)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(PI[a],TXR[b]) == 1):
         print "Error: PI[%d],TXR[%d]"%(a,b)

# chiral mismatch tests

for a in range(0,len(PIS)):
   for b in range(0,len(TIR)):
      if not (ch.match_1v1_b(PIS[a],TIR[b]) == 0):
         print "Error: PIS[%d],TIR[%d]"%(a,b)

for a in range(0,len(TIS)):
   for b in range(0,len(TIR)):
      if not (ch.match_1v1_b(TIS[a],TIR[b]) == 0):
         print "Error: TIS[%d],TIR[%d]"%(a,b)

for a in range(0,len(PXS)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(PXS[a],TXR[b]) == 0):
         print "Error: PXS[%d],TXR[%d]"%(a,b)

for a in range(0,len(TXS)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(TXS[a],TXR[b]) == 0):
         print "Error: TXS[%d],TXR[%d]"%(a,b)

for a in range(0,len(PIS)):
   for b in range(0,len(TXR)):
      if not (ch.match_1v1_b(PIS[a],TXR[b]) == 0):
         print "Error: PIS[%d],TXR[%d]"%(a,b)

for a in range(0,len(PIR)):
   for b in range(0,len(TIS)):
      if not (ch.match_1v1_b(PIR[a],TIS[b]) == 0):
         print "Error: PIR[%d],TIS[%d]"%(a,b)

for a in range(0,len(TIR)):
   for b in range(0,len(TIS)):
      if not (ch.match_1v1_b(TIR[a],TIS[b]) == 0):
         print "Error: TIR[%d],TIS[%d]"%(a,b)

for a in range(0,len(PXR)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(PXR[a],TXS[b]) == 0):
         print "Error: PXR[%d],TXS[%d]"%(a,b)

for a in range(0,len(TXR)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(TXR[a],TXS[b]) == 0):
         print "Error: TXR[%d],TXS[%d]"%(a,b)

for a in range(0,len(PIR)):
   for b in range(0,len(TXS)):
      if not (ch.match_1v1_b(PIR[a],TXS[b]) == 0):
         print "Error: PIR[%d],TXS[%d]"%(a,b)

# chiral->achiral mismatch tests

for a in range(0,len(PIS)):
   for b in range(0,len(TI)):
      if not (ch.match_1v1_b(PIS[a],TI[b]) == 0):
         print "Error: PIS[%d],TI[%d]"%(a,b)

for a in range(0,len(TIS)):
   for b in range(0,len(TI)):
      if not (ch.match_1v1_b(TIS[a],TI[b]) == 0):
         print "Error: TIS[%d],TI[%d]"%(a,b)

for a in range(0,len(PXS)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(PXS[a],TX[b]) == 0):
         print "Error: PXS[%d],TX[%d]"%(a,b)

for a in range(0,len(TXS)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(TXS[a],TX[b]) == 0):
         print "Error: TXS[%d],TX[%d]"%(a,b)

for a in range(0,len(PIS)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(PIS[a],TX[b]) == 0):
         print "Error: PIS[%d],TX[%d]"%(a,b)

for a in range(0,len(PIR)):
   for b in range(0,len(TI)):
      if not (ch.match_1v1_b(PIR[a],TI[b]) == 0):
         print "Error: PIR[%d],TI[%d]"%(a,b)

for a in range(0,len(TIR)):
   for b in range(0,len(TI)):
      if not (ch.match_1v1_b(TIR[a],TI[b]) == 0):
         print "Error: TIR[%d],TI[%d]"%(a,b)

for a in range(0,len(PXR)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(PXR[a],TX[b]) == 0):
         print "Error: PXR[%d],TX[%d]"%(a,b)

for a in range(0,len(TXR)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(TXR[a],TX[b]) == 0):
         print "Error: TXR[%d],TX[%d]"%(a,b)

for a in range(0,len(PIR)):
   for b in range(0,len(TX)):
      if not (ch.match_1v1_b(PIR[a],TX[b]) == 0):
         print "Error: PIR[%d],TX[%d]"%(a,b)

print " test complete."

