import string
import os
import sys

# ChemPy Simple Forcefield (CSFF)
#
# A simplified force field based on Amber Parm99
# ===================================================
# This force field is  designed to be used to answer
# general questions about shape and size of small
# conformationally constrained fragments and to 
# provide a means of ramming random organic compounds
# into the Amber force-field without much concern for
# correctness.  THESE PARAMETERS WILL *NOT* PRODUCE
# REALISTIC RESULTS - SO BE JUDICIOUS WITH THEIR USE
#
#
# A                 hydrogen
# D2                bivalent carbon (nitrile, alkyne)
# D3                trivalent carbon, nonaromatic
# D4                tetravalent carbon
# DA                aromatic, non-junction carbon
# DJ                aromatic junction carbon
# J1                monovalent nitrogen
# J3                planer trivalent nitrogen
# J4                tetrahedral (tri or pentavalent) nitrogen 
# JA                planer aromatic nitrogen 
# JN                sp2 negatively charged nitrogen
# Q1                carbonyl oxygen
# Q2                ether oxygen
# QA                delocalized cationic oxygen
# QN                formal -1 oxygen or delocalized equivalent
# R1                fluorine
# R2                chlorine
# R3                bromine
# R4                iodine
# T1                thiocarbonyl
# T2                bivalent sulfer
# TA                delocalized sulfer in aromatic ring 
# T3                sulfoxide
# T4                tetravalent sulfer

def load(self,fname):
   f = open(fname)
   # skip
   l = f.readline()
   # read names & molecular weights
   self.type = []
   self.mw = {}
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
      a2 = string.strip(l[0:2])
      self.type.append(a2)
      if not self.mw.has_key(a2):
         self.mw[a2] = []
      self.mw[a2].append([l[3:]])
   # skip 1
   l = f.readline()
   # read bonds
   self.bond = {}
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
      a5 = l[0:5]
      if a5[0:2]>a5[3:5]:
         a5 = a5[3:5]+'-'+a5[0:2]
      if not self.bond.has_key(a5):
         self.bond[a5] = []
      self.bond[a5].append([l[5:]])
   # read angles
   self.angle = {}
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
      a5 = l[0:8]
      if a5[0:2]>a5[6:8]:
         a5 = a5[6:8]+'-'+a5[3:5]+'-'+a5[0:2]
      if not self.angle.has_key(a5):
         self.angle[a5] = []
      self.angle[a5].append([l[8:]])
   # read torsion 
   self.torsion = {}
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
      a5 = l[0:11]
      if not self.torsion.has_key(a5):
         self.torsion[a5] = []
      self.torsion[a5].append([l[11:]])         

   # read impropers
   self.improper = {}
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
      a5 = l[0:11]
      if not self.improper.has_key(a5):
         self.improper[a5] = []
      self.improper[a5].append([l[11:]])
   # skip
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
   # read vdw equivalents 
   self.vdw_eq = {}
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
      a4 = string.strip(l[0:4])
      l = l[4:]
      while len(l):
         self.vdw_eq[string.strip(l[0:4])] = a4
         l = l[4:]
   # skip
   l = string.strip(f.readline())
   # read vdw parameters
   self.vdw = {}
   while 1:
      l = string.strip(f.readline())
      if not len(l): break
      l = '  ' + l
      a4 = string.strip(l[0:4])
      self.vdw[a4] =  [float(l[4:20]),
                       float(l[20:37]),
                       string.strip(l[37:])]

   # read extra tinker information if present
   self.extra = {}
   while 1:
      l = f.readline()
      if not l: break
      if l[0:6] == 'TINKER':
         self.extra[string.strip(l[6:12])]  = [
            int(l[12:18]),
            int(l[18:24])]
   

class BlankObject:
   pass

os.system("sed 's/CT/D4/g;s/C[DM]/D3/g;s/C[\*AKQRVW]/DA/g;"+
          "s/-C /-DJ/g;s/C -/DJ-/g;s/C[BCN]/DJ/g;s/C[ZY]/D2/g;"+
          "s/C[FX]/DA/g' parm_simple.inp > tmp1.dat")

os.system("sed 's/N[3T]/J4/g;s/-N /-J3/g;s/N -/J3-/g;s/N\*/J3/g;"+
          "s/N[A2]/J3/g;s/N[BCX]/JA/g;s/NY/J1/g' tmp1.dat > tmp2.dat")

os.system("sed 's/H[OWSAC123P45Z]/A /g;s/-H /-A /g;s/H -/A -/g' "+
          "tmp2.dat > tmp3.dat")

os.system("sed 's/O2/QN/g;s/-O /-Q1/g;s/O -/Q1-/g;s/O[HSW]/Q2/g;' "+
          "tmp3.dat > tmp4.dat")

os.system("sed 's/SO/T4/g;s/SX/TA/g;s/SH/T2/g;s/-S /-T2/g;s/S -/T2-/g;' "+
          "tmp4.dat > tmp5.dat")

os.system("sed 's/-F /-R1/g;s/F -/R1-/g;s/-I /-R4/g;s/I -/R4-/g;"+
          "s/Cl/R2/g;s/Br/R3/g;' tmp5.dat > tmp6.dat")


tmp = BlankObject()


print "CSFF: Chemical Python Simplified Force Field by Warren L. DeLano"

load(tmp,'tmp6.dat')

kees = tmp.mw.keys()
kees.sort()
for a in kees:
   c = ''
   lst = ''
   for b in tmp.mw[a]:
      if lst!=b[0][0:8]:
         print "%s%-2s %s" %(c,a,b[0])
         c = ' '
         lst = b[0][0:8]
l = [
'QA 16.00         0.434               delocalized, cationic oxygen',
'T3 32.06         2.900               Sulfoxide',
'T1 32.06         2.900               Thiocarbonyl',
'T3 32.06         2.900               Sulfoxide',
]   
for a in l:
   print a 
print
print "A   J1  J2  J3  J4  JN Q1  Q2  QN"

kees = tmp.bond.keys()
kees.sort()
for a in kees:
   if len(tmp.bond[a])==1:
      print "%-1s%s" %(a,tmp.bond[a][0][0])
   else:
      f1 = 0.0
      f2 = 0.0
      c = 0
      for b in tmp.bond[a]:
         f1 = f1 + float(b[0][0:9])
         f2 = f2 + float(b[0][9:18])
         c = c + 1
      f1 = f1 / c
      f2 = f2 / c
      print "%-1s%7.1f%9.3f       combination of %d"%(a,f1,f2,c)

# missing bond terms
l = [
'A -QN  434.0    1.010       TEMPORARY NUSIANCE REMOVAL',
'A -JN  434.0    1.010       TEMPORARY NUSIANCE REMOVAL',
'A -JA  434.0    1.010       WLD from A -J3',
'D2-DJ  428.0    1.425       WLD from CA-D2',
'D2-D3  450.0    1.331       WLD gross estimate',
'D2-J3  428.0    1.425       WLD from CA-D2',
'D2-T2  375.0    1.750       WLD gross estimate',
'D3-JA  448.0    1.365       WLD from D3-J3',
'D3-R1  386.0    1.359       WLD from DA-R1',
'D3-R2  310.0    1.724       WLD from DJ-R2',
'D3-R3  172.0    1.890       WLD from DJ-R3',
'D3-R4  171.0    2.075       WLD from DJ-R4',
'D3-T2  367.0    1.715       WLD from DA-T2',
'D3-TA  367.0    1.715       WLD from DA-T2',
'DJ-T1  500.0    1.629       WLD gross estimate',

'D4-TA  232.0    1.810       WLD from D4-T2',
'D4-T4  232.0    1.810       WLD from D4-T2',
'D4-JA  337.0    1.463       WLD from CT-N2',
'DJ-R1  386.0    1.359       WLD from DA-R1',
'DJ-R2  310.0    1.724       WLD from DJ-R2',
'DJ-R3  172.0    1.890       WLD from DJ-R3',
'DJ-R4  171.0    2.075       WLD from DJ-R4',
'DJ-T4  427.0    1.761       WLD from DA-T4',
'DJ-J4  448.0    1.365       WLD from D3-J3',
'DJ-JN  448.0    1.365       WLD from D3-J3',
'DJ-T2  367.0    1.715       WLD from DA-T2',
'DJ-TA  327.0    1.740       WLD from DA-TA',
'J3-JA  453.7    1.360       WLD gross estimate',
'J3-J3  453.7    1.360       WLD gross estimate',
'J3-J4  453.7    1.430       WLD gross estimate',
'J3-Q2  453.7    1.400       WLD gross estimate',
'JA-JA  453.7    1.360       WLD gross estimate',
'JA-Q2  453.7    1.370       WLD gross estimate',
'JA-T4  427.0    1.680       WLD gross estimate',
'J3-QN  656.0    1.250       WLD gross estimate',
'J4-R2  310.0    1.700       WLD gross estimate',
'JN-T4  230.0    1.615       WLD from J4-T4',
'JA-T4  400.0    1.600       WLD, gross estimate',
'Q2-T4  525.0    1.430       WLD from Q1-T4',
'QA-DA  400.0    1.292       WLD, gross estimate',
]

for a in l:
   if tmp.bond.has_key(a[0:5]):
      sys.stderr.write("Duplicate bond: %s\n"%a)
   print a 
print
kees = tmp.angle.keys()
kees.sort()
for a in kees:
   if len(tmp.angle[a])==1:
      print "%-1s%s" %(a,tmp.angle[a][0][0])
   else:
      f1 = 0.0
      f2 = 0.0
      c = 0
      for b in tmp.angle[a]:
         f1 = f1 + float(b[0][0:10])
         f2 = f2 + float(b[0][10:22])
         c = c + 1
      f1 = f1 / c
      f2 = f2 / c
      print "%-1s%8.1f%12.2f    combination of %d"%(a,f1,f2,c)
# missing angle terms
l = [
'J3-QN-A     50.0      109.50    TEMPORARY NUSIANCE REMOVAL', # (atom should not be protonated...)
'A -D3-JA    50.0      119.10    WLD from A -D3-J3',
'A -D3-R1    50.0      119.10    WLD from A -D3-J3',
'A -D3-R2    50.0      119.10    WLD from A -D3-J3',
'A -D3-R3    50.0      119.10    WLD from A -D3-J3',
'A -D3-R4    50.0      119.10    WLD from A -D3-J3',

'A -D4-JA    50.0      109.50    WLD from A -D4-J3',
'A -D4-TA    50.0      109.50    WLD from A -D4-T2',

'A -DA-Q2    45.0      120.15    WLD from A -DA-JA',
'A -DA-D4    50.0      120.00    WLD from A -D3-D4',
'A -DA-TA    35.0      122.60    WLD from A -DA-T2',
'A -D4-T4    30.0      115.43    WLD from A- J4-T4',
'A -JA-D3    50.0      121.20    WLD from A -J3-D3',
'A -J3-D2    46.2      111.00    WLD gross estimate',
'A -J3-JA    45.0      120.49    WLD from A -J3-DA',
'A -J3-J3    45.0      120.49    WLD from A -J3-DA',
'A -J3-J4    50.0      118.22    WLD from A -J3-D4',
'A -J3-Q2    45.0      120.49    WLD from A -J3-DA',
'A -JA-DJ    47.9      120.70    WLD from A -J3-DJ',
'A -JA-T4    47.9      120.70    WLD from A -J3-DJ',
'A -J4-DJ    30.0      115.43    WLD from A -J4-T4',
'A -J4-J3    30.0      115.43    WLD from A -J4-T4',
'A -JN-DJ    30.0      115.43    WLD from A -J4-T4',
'A -JN-T4    30.0      115.43    WLD from A -J4-T4',
'A -Q2-D3    50.0      113.00    WLD from A -Q2-DA',
'A -T2-DJ    62.0       98.90    WLD from D4-T2-D4',
'A -T2-T2    68.0      103.70    WLD from D4-T2-T2',

'A -Q2-JA    50.0      109.50    WLD gross estimate',
'A -Q2-T4    50.0      113.00    WLD from A -Q2-DJ',
'A -Q2-J3    50.0      109.50    WLD gross estimate',
'A -DA-D3    46.2      120.76    WLD from A -DA-J3',
'A -T2-D2    43.0       96.00    WLD from A -T2-D4',
'D2-DJ-Q1    80.0      125.30    WLD from D3-DJ-Q1',
'D2-DJ-Q2    80.0      125.30    WLD from D3-DJ-Q1',
'D2-DJ-QN    80.0      125.30    WLD from D3-DJ-Q1',
'D2-D2-DJ    70.0      180.00    WLD from DJ-D2-D2',
'D2-D3-D3    63.0      120.00    WLD from D3-D3-D3',
'D2-D4-J3    40.0      111.00    WLD from D2-D4-DA',
'D2-D4-JA    40.0      111.00    WLD from D2-D4-DA',
'D2-D4-DJ    40.0      111.00    WLD from D2-D4-DA',
'D2-DJ-DA    70.0      126.55    WLD from D2-DA-DA',
'D2-DJ-JA    70.0      126.55    WLD from D2-DA-DA',
'D2-DJ-DJ    70.0      126.55    WLD from D2-DA-DA',
'D2-J3-DJ    70.0      126.55    WLD from D2-DA-DA',
'D2-T2-D4    62.0       96.00    WLD gross estimate',
'D2-T2-DJ    80.0       90.10    WLD from DA-T2-DA',

'D3-D3-R1    63.0      120.00    WLD from D3-D3-D3',
'D3-D3-R2    63.0      120.00    WLD from D3-D3-D3',
'D3-D3-R3    63.0      120.00    WLD from D3-D3-D3',
'D3-D3-R4    63.0      120.00    WLD from D3-D3-D3',
'D3-D3-T2    63.0      120.00    WLD from D3-D3-D3',

'DJ-D3-R1    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-R2    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-R3    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-R4    63.0      120.00    WLD from D3-D3-D3',
'DJ-DJ-T1    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-T2    63.0      120.00    WLD from D3-D3-D3',
'D3-D3-JA    63.0      120.00    WLD from D3-D3-D3',

'D3-D4-D3    50.0      109.50    WLD from J3-D4-J3',
'D3-D4-J3    50.0      109.50    WLD from J3-D4-J3',
'D3-D4-R1    50.0      114.00    WLD from D4-D4-DA',
'D3-D4-R2    50.0      114.00    WLD from D4-D4-DA',
'D3-D4-R3    50.0      114.00    WLD from D4-D4-DA',
'D3-D4-R4    50.0      114.00    WLD from D4-D4-DA',


'D3-DJ-DJ    63.0      120.70    WLD from D3-D3-DJ',
'D3-DJ-Q2    63.0      120.00    WLD gross estimate',
'D3-J3-J3    70.0      120.00    WLD gross estimate',
'D3-JA-Q2    70.0      117.80    WLD from DA-JA-DA',
'D3-JA-JA    70.0      117.80    WLD from DA-JA-DA',
'D3-JA-J3    70.0      117.80    WLD from DA-JA-DA',
'D3-JA-D4    70.0      117.80    WLD from DA-JA-DA',
'D3-DJ-D4    70.0      119.70    WLD from D4-D3-DJ',
'D3-DJ-D3    70.0      114.10    WLD from D3-DJ-J3',
'D3-DJ-JA    70.0      114.10    WLD from D3-DJ-J3',
'D3-DJ-QN    70.0      120.00    WLD from DA-DJ-Q2',
'D3-DJ-TA    70.0      120.00    WLD from DA-DA-T4',
'D3-D4-DJ    63.0      105.85    WLD from DJ-D4-J3',

'D3-J3-D3    70.0      121.60    WLD from D3-J3-DJ',
'D3-JA-DJ    70.0      117.80    WLD from DA-JA-DA',

'D3-TA-D4    62.0      107.10    WLD from DJ-T4-J4',

'D4-D3-Q2    80.0      125.00    WLD from D3-D3-Q2',

'D4-D4-JA    72.5      110.02    WLD from D4-D4-J3',
'D4-TA-DJ    62.0      107.10    WLD from DJ-T4-J4',
'D4-T2-D3    62.0       98.00    WLD gross estimate',
'D4-D4-T4    50.0      121.20    WLD from D4-J4-T4',
'D4-D4-TA    50.0      121.20    WLD from D4-J4-T4',
'D4-DJ-DJ    70.0      122.50    WLD from D4-DA-DA',
'D4-DJ-TA    63.0      117.00    WLD from D4-DJ-D4',
'D4-DJ-JN    63.0      117.00    WLD from D4-DJ-D4',
'D4-J3-JA    65.0      123.97    WLD from D4-J3-DA',
'D4-J3-JA    65.0      123.97    WLD from D4-J3-DA',
'D4-J3-J3    70.0      120.00    WLD gross estimate',
'D4-J4-DJ    80.0      111.20    WLD from DJ-D4-J4',
'D4-J4-J3    80.0      111.20    WLD from DJ-D4-J4',
'D4-T2-DJ    62.0       98.00    WLD gross estimate',
'D4-T4-D4    62.0      103.00    WLD gross estimate',
'D4-Q2-DA    57.5      118.22    WLD from D4-Q2-DJ',
'D4-Q2-JA    57.5      118.22    WLD from D4-Q2-DJ',
'D4-Q2-J3    57.5      118.22    WLD from D4-Q2-DJ',
'D4-Q2-T4    60.0      109.50    WLD from D4-Q2-D4',

'D4-D3-JA    70.0      119.70    WLD from D4-D3-DJ',

'D4-J3-QN    66.0      121.06    WLD from D4-J3-DJ',
'D4-J3-Q2    66.0      121.06    WLD from D4-J3-DJ',
'D4-JA-JA    70.0      120.00    WLD gross estimate',

'D4-T4-J4   100.0      106.80    WLD from J4-T4-J4',
'D4-T4-Q1   100.0      107.10    WLD from J4-T4-Q1',
'D4-T4-Q2   100.0      107.10    WLD from J4-T4-Q1',
'D4-T4-QN   100.0      107.10    WLD from J4-T4-QN',
'D4-T4-DJ    62.0      107.10    WLD from DJ-T4-J4',

'DJ-Q2-T4    57.5      118.22    WLD from D4-Q2-DJ',
'DJ-D3-D2    63.0      119.20    WLD from DJ-DJ-DJ',
'DJ-D3-Q2    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-R1    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-R2    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-R3    63.0      120.00    WLD from D3-D3-D3',
'DJ-D3-R4    63.0      120.00    WLD from D3-D3-D3',

'D4-D3-D4    63.0      117.00    WLD from D4-DJ-D4',
'DA-D3-J3    63.0      117.00    WLD from D3-D3-DA',
'DA-D3-DA    63.0      117.00    WLD from D3-D3-DA',
'DA-D3-JA    63.0      117.00    WLD from D3-D3-DA',
'DA-D4-DJ    63.0      105.85    WLD from DJ-D4-J3',

'DA-DA-D3    68.8      118.12    WLD from DA-DA-J3',
'DA-DA-QA    70.0      117.50    WLD from DA-DA-JA',

'DA-DJ-T2    70.0      120.00    WLD from DA-DA-T4',
'DA-DJ-D3    63.0      120.70    WLD from D3-D3-DJ',
'DA-DJ-J4    70.0      120.00    WLD from DA-DA-J4',
'DA-DJ-T2    70.0      109.20    WLD from DA-DA-TA',
'DA-DJ-R1    70.0      121.00    WLD from DA-DA-R1',
'DA-DJ-R2    70.0      119.40    WLD from DA-DA-R2',
'DA-DJ-R3    70.0      118.80    WLD from DA-DA-R3',
'DA-DJ-R4    70.0      118.80    WLD from DA-DA-R4',
'DA-DJ-TA    63.0      116.25    WLD from DA-DJ-DJ',
'DA-DJ-T4    70.0      120.00    WLD from DA-DA-T4',

'DA-JA-JA    70.0      117.80    WLD from DA-JA-DA',
'DA-JA-J3    70.0      117.80    WLD from DA-JA-DA',
'DA-JA-Q2    70.0      117.80    WLD from DA-JA-DA',
'DA-JA-TA    70.0      120.00    WLD from DA-DA-TA',
'DA-J3-JA    70.0      117.80    WLD from DA-JA-DA',
'DA-J3-Q2    70.0      117.80    WLD from DA-JA-DA',
'DA-J3-QN    70.0      117.80    WLD from DA-JA-DA',

'DA-Q2-J3    57.5      118.22    WLD from D4-Q2-DJ',
'DA-QA-DA    70.0      114.24    WLD from DA-JA-DA',

'DA-T4-Q2   100.0      107.80    WLD from DA-T4-QN',


'DJ-D3-DJ    63.0      119.20    WLD from DJ-DJ-DJ',
'DJ-D3-JA    70.0      117.50    WLD from DA-DA-JA',

'DJ-D4-DJ    40.0      109.50    WLD from D4-D4-D4',
'DJ-D4-R1    50.0      114.00    WLD from D4-D4-DA',
'DJ-D4-R2    50.0      114.00    WLD from D4-D4-DA',
'DJ-D4-R3    50.0      114.00    WLD from D4-D4-DA',
'DJ-D4-R4    50.0      114.00    WLD from D4-D4-DA',
'DJ-D4-T2    80.0      111.20    WLD from DJ-J4-T4',
'DJ-D4-TA    63.0      111.00    WLD gross estimate',
'DJ-D4-T4    80.0      111.20    WLD from DJ-J4-T4',

'DJ-J3-J4    70.0      121.20    WLD from D3-J3-D4',
'DJ-J3-Q2    70.0      117.80    WLD from DA-JA-DA',

'DJ-D4-JA    63.0      105.85    WLD from DJ-D4-J3',

'DJ-D3-J3    70.0      117.50    WLD from DA-DA-JA',

'DJ-DA-Q2    70.0      113.91    WLD from DA-DA-Q2',
'DJ-DA-QA    70.0      117.50    WLD from DA-DA-JA',
'DJ-DA-TA    70.0      120.00    WLD from DA-DA-TA',

'DJ-DJ-J4    70.0      120.00    WLD from D4-DJ-JA',
'DJ-DJ-TA    63.0      116.25    WLD from DA-DJ-DJ',
'DJ-DJ-T2    70.0      120.00    WLD from DA-DA-T4',
'DJ-DJ-T4    70.0      120.00    WLD from DA-DA-T4',
'DJ-DJ-QN    70.0      120.00    WLD from DA-DJ-Q2',
'DJ-DJ-R1    70.0      121.00    WLD from DA-DA-R1',
'DJ-DJ-R2    70.0      119.40    WLD from DA-DA-R2',
'DJ-DJ-R3    70.0      118.80    WLD from DA-DA-R3',
'DJ-DJ-R4    70.0      118.80    WLD from DA-DA-R4',

'DJ-J4-T4    80.0      111.20    WLD from DJ-D4-J4',

'DJ-JA-DJ    70.0      123.20    WLD from DJ-J3-DJ',
'DJ-JA-TA    70.0      110.20    WLD from DJ-DA-T2',
'DJ-JA-J3    70.0      117.80    WLD from DA-JA-DA',
'DJ-JA-JA    70.0      117.80    WLD from DA-JA-DA',
'DJ-JA-Q2    70.0      117.80    WLD from DA-JA-DA',
'DJ-J3-J3    70.0      117.80    WLD from DA-JA-DA',
'DJ-J3-JA    70.0      120.00    WLD from DA-J3-DA',
'DJ-J3-QN    70.0      120.00    WLD from DA-DJ-Q2',
'DJ-JN-T4    80.0      111.20    WLD from DJ-D4-J4',
'DJ-Q1-DJ    70.0      123.20    WLD from DJ-J3-DJ',
'DJ-Q2-JA    70.0      117.80    WLD from DA-JA-DA',
'DJ-TA-DJ    62.0       98.90    WLD from D4-T2-D4',
'DJ-TA-DA    62.0       98.90    WLD from D4-T2-D4',
'DJ-T2-DJ    80.0       90.10    WLD from DA-T2-DA',
'DJ-T4-DJ    62.0       98.90    WLD from D4-T2-D4',
'DJ-T4-J4    62.0      107.10    WLD from DA-T4-J4',
'DJ-T4-JN    62.0      107.10    WLD from DA-T4-J4',
'DJ-T4-Q1   100.0      107.80    WLD from DA-T4-Q1',
'DJ-T4-QN   100.0      107.80    WLD from DA-T4-QN',
'DJ-T4-Q2   100.0      107.80    WLD from DJ-T4-Q1',
'DJ-TA-A     62.0       98.90    WLD from DJ-T2-DJ',

'J1-D2-DJ    70.0      180.00    WLD from DJ-D2-D2',
'J1-D2-D3    70.0      180.00    WLD from DJ-D2-D2',
'J1-D2-J3    70.0      180.00    WLD from DJ-D2-D2',
'J1-D2-T2    70.0      180.00    WLD from DJ-D2-D2',

'J3-D3-JA    70.0      113.13    WLD from J3-DJ-J3',
'J3-D3-JN    70.0      113.13    WLD from J3-DJ-J3',
'J3-DJ-JN    70.0      113.13    WLD from J3-DJ-J3',
'J3-DA-TA    70.0      122.40    WLD from J3-DJ-JA',
'J3-DJ-TA    70.0      122.40    WLD from J3-DJ-JA',
'J3-DJ-T2    70.0      122.40    WLD from J3-DJ-JA',
'JN-DJ-Q1    77.5      123.10    WLD from J3-DJ-Q1',
'J3-DJ-Q2    70.0      122.50    WLD from DA-DJ-Q2',
'J3-DJ-T1    70.0      121.00    WLD gross estimate',
'J3-D3-T2    70.0      123.20    WLD from J3-DA-T2',
'J3-DJ-T2    70.0      123.20    WLD from J3-DA-T2',
'J3-D3-J3    63.0      117.00    WLD from D3-D3-DA',
'J3-D4-T2    50.0      111.65    WLD from D4-D4-T2',
'J3-J4-T4    50.0      109.30    WLD from DA-D4-T2',
'J3-JA-JA    70.0      117.80    WLD from DA-JA-DA',

'J3-J3-QN    70.0      117.80    WLD from DA-JA-DA',
'J3-DJ-T1    70.0      121.00    WLD gross estimate',
'J3-D3-TA    63.0      120.00    WLD from D3-D3-D3',

'J4-D4-Q2    50.0      109.50    WLD from D4-D4-Q2',
'J4-DJ-JA    70.0      120.00    WLD gross estimate',
'J4-DJ-Q2    70.0      120.00    WLD gross estimate',
'J4-DJ-TA    63.0      116.25    WLD from DA-DJ-DJ',

'JA-D3-T2    70.0      123.20    WLD from J3-DA-T2',
'JA-D3-JA    70.0      113.13    WLD from J3-DJ-J3',
'JA-D3-Q2    80.0      125.00    WLD from D3-D3-Q2',
'JA-DJ-JA    70.0      113.13    WLD from J3-DJ-J3',
'JA-DJ-Q2    70.0      113.13    WLD from J3-DJ-J3',
'JA-DJ-R2    70.0      119.40    WLD from DA-DA-R2',
'JA-DJ-T2    63.0      116.25    WLD from DA-DJ-DJ',
'JA-DJ-TA    63.0      116.25    WLD from DA-DJ-DJ',
'JA-TA-JA    62.0       98.90    WLD from D4-T2-D4',
'JA-TA-DJ    62.0       98.90    WLD from D4-T2-D4',
'JA-DJ-R1    70.0      121.00    WLD from DA-DA-R1',
'JA-DJ-R2    70.0      119.40    WLD from DA-DA-R2',
'JA-DJ-R3    70.0      118.80    WLD from DA-DA-R3',
'JA-DJ-R4    70.0      118.80    WLD from DA-DA-R4',
'JA-J3-JA    63.3      122.80    WLD form DA-J3-DA',
'JA-J3-Q2    63.3      122.80    WLD form DA-J3-DA',
'JN-T4-Q1   100.0      107.10    WLD from J4-T4-Q1',
'JA-T4-Q1   100.0      107.10    WLD from J4-T4-Q1',
'JA-T4-D4   100.0      107.10    WLD from J4-T4-Q1',
'JA-D3-TA    70.0      123.20    WLD from J3-DA-T2',


'R1-DJ-TA    70.0      120.30    WLD from D4-DA-TA',

'R2-DJ-TA    70.0      120.30    WLD from D4-DA-TA',
'R1-D3-R1    63.0      120.00    WLD from D3-D3-D3',
'R2-D3-R2    63.0      120.00    WLD from D3-D3-D3',

'R3-DJ-TA    70.0      120.30    WLD from D4-DA-TA',
'R3-D3-R3    63.0      120.00    WLD from D3-D3-D3',
'R4-DJ-TA    70.0      120.30    WLD from D4-DA-TA',

'R1-D4-R2    70.0      109.50    WLD gross estimate',
'R1-D4-R3    70.0      109.50    WLD gross estimate',
'R1-D4-R4    70.0      109.50    WLD gross estimate',
'R2-D4-R2    70.0      109.50    WLD gross estimate',
'R2-D4-R3    70.0      109.50    WLD gross estimate',
'R2-D4-R4    70.0      109.50    WLD gross estimate',

'R2-J4-T4    80.0      111.20    WLD from DJ-D4-J4',
'R2-J4-R2    80.0      111.20    WLD from DJ-D4-J4',

'R3-D4-R3    70.0      109.50    WLD gross estimate',
'R3-D4-R4    70.0      109.50    WLD gross estimate',
'R4-D4-R4    70.0      109.50    WLD gross estimate',



'Q1-DJ-TA    63.0      116.25    WLD from DA-DJ-DJ',
'Q1-T4-Q2   140.0      119.70    WLD from Q1-T4-Q1',
'Q1-DJ-T2    80.0      125.30    WLD from D3-DJ-Q1',
'Q1-T4-Q2    70.0      109.50    WLD gross estimate',


'Q2-D4-R1    70.0      109.50    WLD gross estimate',
'Q2-DJ-R1    75.0      120.00    WLD from DA-DJ-Q1',
'Q2-DJ-R2    75.0      120.00    WLD from DA-DJ-Q1',
'Q2-DJ-R3    75.0      120.00    WLD from DA-DJ-Q1',
'Q2-DJ-R4    75.0      120.00    WLD from DA-DJ-Q1',
'Q2-T4-QN    70.0      109.50    WLD gross estimate',
'Q2-T4-Q2    70.0      109.50    WLD gross estimate',

'QN-J3-QN    80.0      126.00    WLD from QN-DJ-QN',
'QN-T4-QN    70.0      109.50    WLD gross estimate',

'T1-D3-T2    70.0      121.00    WLD gross estimate',
'T2-D3-T2    70.0      118.00    WLD gross estimate',
'T1-DJ-T2    70.0      121.00    WLD gross estimate',
'T2-DJ-T2    70.0      118.00    WLD gross estimate',


'T2-DJ-T4    63.0      117.00    WLD from D4-DJ-D4',
'T2-D4-T2    80.0      111.20    WLD from DJ-J4-T4',

'TA-DJ-T1    70.0      121.00    WLD gross estimate',
'TA-DJ-T4    63.0      117.00    WLD from D4-DJ-D4',
'T4-J4-T4    50.0      121.20    WLD from D4-J4-T4',
]

for a in l:
   if tmp.angle.has_key(a[0:8]):
      sys.stderr.write("Duplicate angle: %s\n"%a)
   print a 
print
# missing generalized torsions (divisors/forces need to be checked...)
l = [
'X -D2-DJ-X    2    0.00          0.0             2.         WLD null',
'X -D2-D3-X    2    0.00          0.0             2.         WLD null',
'X -D2-J3-X    2    0.00          0.0             2.         WLD null',
'X -D2-T2-X    2    0.00          0.0             2.         WLD null',
'X -D3-JA-X    4    8.70        180.0             2.         WLD from X -DJ-D3-X ',
'X -D4-TA-X    6    6.00        180.0             2.         WLD from X -DA-TA-X ',
'X -D4-T4-X    6    2.40          0.0             3.         WLD from X -J4-T4-X ',
'X -D4-JA-X    6    0.00          0.0             2.         WLD from X -DA-T4-X ',
'X -DJ-T4-X    6    0.00          0.0             2.         WLD from X -DA-T4-X ',
'X -DJ-T2-X    6    0.00          0.0             2.         WLD from X -DA-T4-X ',
'X -D3-T2-X    6    0.00          0.0             2.         WLD from X -DA-T4-X ',
'X -D3-TA-X    6    0.00          0.0             2.         WLD from X -DA-T4-X ',
'X -DJ-J4-X    6    0.00          0.0             2.         WLD from X -DJ-D4-X ',
'X -DJ-TA-X    2    6.00        180.0             2.         WLD from X -DA-TA-X',
'X -J3-J3-X    2   11.60        180.0             2.         WLD from X -DA-JA-X ',
'X -DA-QA-X    2   11.60        180.0             2.         WLD from X -DA-JA-X ',
'X -J3-J3-X    2   11.60        180.0             2.         WLD from X -DA-JA-X ',
'X -J3-JA-X    2   11.60        180.0             2.         WLD from X -DA-JA-X ',
'X -J3-J4-X    2    0.00          0.0             2.         WLD null',
'X -J3-QN-X    4   11.20        180.0             2.         WLD from X -DJ-Q1-X ',
'X -JA-JA-X    2   11.60        180.0             2.         WLD from X -DA-JA-X ',
'X -JA-Q2-X    2    4.20        180.0             2.         WLD from X -DA-Q2-X ',
'X -J3-Q2-X    3    0.82          0.0             3.         WLD from X -D4-D2-X ',
'X -JA-TA-X    2    6.00        180.0             2.         WLD from X -DA-TA-X ',
'X -JA-T4-X    2    6.00        180.0             2.         WLD from X -DA-TA-X ',
'X -Q2-T4-X    4    2.40          0.0             3.         WLD from X -J4-T4-X ',
'X -JN-T4-X    4    2.40          0.0             3.         WLD from X -J4-T4-X ',
'X -DJ-JN-X    4    7.01        180.0             2.         WLD from X -DJ-J3-X ',
]
for a in l:
   if tmp.torsion.has_key(a[0:11]):
      sys.stderr.write("Duplicate torsion: %s\n"%a)
   print a 
kees = tmp.torsion.keys()
kees.sort()
kees.reverse()
for a in kees:
   if len(tmp.torsion[a])==1:
      print "%-1s%s" %(a,tmp.torsion[a][0][0])
   else:

      b = tmp.torsion[a][0]
      f1 = float(b[0][0:6])
      f2 = float(b[0][6:14])
      f3 = float(b[0][14:27])
      f4 = float(b[0][27:42])         
      c = 1
      flag = 0
      for b in tmp.torsion[a][1:]:
         if ((f1 != float(b[0][0:6])) or
             (f3 != float(b[0][14:27])) or
             (f4 != float(b[0][27:42]))):
            flag=1
            break
         f2 = f2 + float(b[0][6:14])
         c = c + 1
      if not flag:
         f2 = f2 / c
         print "%-1s%4d%8.2f%13.1f%14d.         combination of %d"%(a,f1,f2,f3,f4,c)
      else:
         flag = 0
         ck = {}
         for b in tmp.torsion[a]:
            f4 = float(b[0][27:42])
            if ck.has_key(f4):
               flag=1
               break
            ck[f4] = 1
         if not flag:
            for b in tmp.torsion[a]:
               print "%-1s%s" %(a,b[0])            
         else: # known special cases 
            if a == 'X -D4-J4-X ':
                print "%-1s%s" %(a,tmp.torsion[a][1][0])
            elif a == 'X -D4-J3-X ':
                print "%-1s%s" %(a,tmp.torsion[a][1][0])
            elif a == 'Q2-D4-D4-Q2':
                print "%-1s%s" %(a,tmp.torsion[a][1][0])
                print "%-1s%s" %(a,tmp.torsion[a][2][0])
            elif a == 'A -D4-DJ-Q1':
                print "%-1s%s" %(a,tmp.torsion[a][1][0])
                print "%-1s%s" %(a,tmp.torsion[a][2][0])
            else:
               for b in tmp.torsion[a]:
                  print " %-1s%s" %(a,b[0])
# missing specific torsions
l = [
'A -T2-T2-D4   1    3.50          0.0            -2.         WLD from D4-T2-T2-D4',
'A -T2-T2-D4   1    0.60          0.0             3.         WLD from D4-T2-T2-D4',
]
for a in l:
   if tmp.torsion.has_key(a[0:11]):
      sys.stderr.write("Duplicate torsion: %s\n"%a)
   print a 
print

kees = tmp.improper.keys()
kees.sort()
kees.reverse()
for a in kees: # no major redundancy, just print the first record
   print "%-1s%s" %(a,tmp.improper[a][0][0])
# missing specific impropers
print '''
  A   Q2  0000.     0000.                                4.  flag for fast water

J1  J1  J2  J3  J4
D2  D2  D3  D4  

MOD4      RE
  A           1.3870  0.0157             A Venstra et al JCC,8,(1992),963 
  QN          1.6612  0.2100             O  OPLS
  QA          1.6612  0.2100             O  OPLS
  Q2          1.6612  0.2100             O  OPLS
  Q1          1.6612  0.2100             O2 OPLS
  D4          1.9080  0.1094             Spellmeyer
  D3          1.9080  0.0860             Spellmeyer
  DA          1.9080  0.0860             Spellmeyer
  DJ          1.9080  0.0860             Spellmeyer
  D2          1.9080  0.0860             cp C 
  J1          1.8240  0.1700             OPLS
  J2          1.8240  0.1700             OPLS
  J3          1.8240  0.1700             OPLS
  J4          1.8240  0.1700             OPLS
  JJ          1.8240  0.1700             OPLS
  JA          1.8240  0.1700             OPLS
  JN          1.8240  0.1700             OPLS  
  R1          1.75    0.061              F  Gough et al. JCC 13,(1992),963.
  R2          1.948   0.265              Cl Fox, JPCB,102,8070,(98),flex.mdl CHCl3
  R3          2.22    0.320              Br Junmei(?)
  R4          2.35    0.40               I  JCC,7,(1986),230;  
  T4          2.0000  0.2500             S  W. Cornell CH3SH and CH3SCH3 FEP's
  TA          2.0000  0.2500             S  W. Cornell CH3SH and CH3SCH3 FEP's
  T2          2.0000  0.2500             S  W. Cornell CH3SH and CH3SCH3 FEP's
  T1          2.0000  0.2500             S  W. Cornell CH3SH and CH3SCH3 FEP's
  
END

# hydrogen      
TINKER    A      1    1 
# carbon                        
TINKER    D4     6    4 
TINKER    D3     6    3 
TINKER    DA     6    3 
TINKER    DJ     6    3 
TINKER    D2     6    2 
# fluorine      
TINKER    R1     9    1 
# chloride      
TINKER    R2    17    1  
# bromine 
TINKER    R3    35    1
# iodine        
TINKER    R4    53    1 
# nitrogen      
TINKER    J4     7    4 
TINKER    J3     7    3 
TINKER    JA     7    2
TINKER    JN     7    2
TINKER    J1     7    1
# oxygen            
TINKER    QA     8    3
TINKER    Q2     8    2 
TINKER    Q1     8    1 
TINKER    QN     8    1 
# sulfer
TINKER    T4    16    4
TINKER    T3    16    3 
TINKER    TA    16    2
TINKER    T2    16    2
TINKER    T1    16    1
'''
  
os.unlink('tmp1.dat')
os.unlink('tmp2.dat')
os.unlink('tmp3.dat')
os.unlink('tmp4.dat')
os.unlink('tmp5.dat')
os.unlink('tmp6.dat')
