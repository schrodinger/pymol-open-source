from chempy.champ import Champ

# see if we can read-in and read-out all of the NCI
# database SMILES strings without a hitch

import re
num_re = re.compile("[0-9]")
pnum_re = re.compile(r"\%[0-9][0-9]")

ch = Champ()

f = open("nci.smi")

c = 1
lst = []
while 1:
   l = f.readline()
   if not l: break
   l = string.strip(l)
   if len(l):
      if 1000*(c/1000)==c:
         print c
      c = c + 1
#      print l
      idx = ch.insert_pattern_string(l)
 #      ch.pattern_dump(idx)
      o = ch.get_pattern_string(idx)
      if len(o) != len(l): # same length?
         print "orig:  ",l
         print "champ: ",o
         ch.pattern_dump(idx)
         break
      else:
         # now, anonymize cycle identifiers and then compare...
         o = pnum_re.sub("|",o)
         o = num_re.sub("|",o)
         l = pnum_re.sub("|",l)
         l = num_re.sub("|",l)
         fail = 0
         if o!=l:
            if o!="S=1(=NC(=O)C(=C(S[C])N=1)C#N)([C])[C]":
               print "orig:  ",l
               print "champ: ",o
               ch.pattern_dump(idx)
               fail = 1
               break
         if fail:
            break
      ch.pattern_free(idx)
      
