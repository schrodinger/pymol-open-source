import string

f = open("parm99_wld.dat")
g = open("parm_simple.dat")
h = open("parm99_simple.dat",'w')

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

g.readline() # skip first line
while 1:
   l = g.readline()
   h.write(l)
   if not len(string.strip(l)):
      break
   
l=f.readline()
h.write(l)
g.readline() # skip this stuff for now

# BONDS

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   if not len(string.strip(l)):
      break
   h.write(l)
for l in [
'T2-S   166.0    2.038       WLD from S -S ',
'DJ-N   490.0    1.335       WLD from C -N ',
'C -J3  490.0    1.335       WLD from C -N ',
'D4-H1  340.0    1.090       WLD from CT-H1',
'C -D4  317.0    1.522       WLD from C -CT',
'D4-N   337.0    1.449       WLD from CT-N',
]:
   h.write(l+"\n")
h.write("\n")

# ANGLES

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

for l in [
'H -N -DJ    50.0      120.00    WLD from C -N -H ',   
'A -J3-C     50.0      120.00    WLD from C -N -H ',
'C -J3-D4    50.0      121.90    WLD from C -N -CT',
'CT-S -T2    68.0      103.70    WLD from CT-S -S ',
'CT-N -DJ    50.0      121.90    WLD from C -N -CT',
'CT-C -J3    70.0      116.60    WLD from CT-C -N ',
'D4-T2-S     68.0      103.70    WLD from CT-S -S ',
'D4-DJ-N     70.0      116.60    WLD from CT-C -N ',
'DJ-N -CT    50.0      121.90    WLD from C -N -CT',
'J3-C -O     80.0      122.90    WLD from N -C -O ',
'N -DJ-Q1    80.0      122.90    WLD from N -C -O ',
'N -DJ-DJ    70.0      120.00    WLD from CA-C -OH',
'C -N -D4    50.0      121.90    WLD from C -N -CT',
'D4-C -N     70.0      116.60    WLD from CT-C -N ',
'D4-C -O     80.0      120.40    WLD from CT-C -O ',
'C -D4-D4    63.0      111.10    WLD from C -CT-CT',
'C -D4-H1    50.0      109.50    WLD from C -CT-H1',
'D4-N -H     50.0      118.04    WLD from CT-N -H ',
'H1-D4-N     50.0      109.50    WLD from H1-CT-N ',
'C -D4-N     63.0      110.10    WLD from C -CT-N ',
'D4-D4-N     80.0      109.70    WLD from CT-CT-N ',
'D4-D4-H1    50.0      109.50    WLD from CT-CT-H1',
]:
   h.write(l+"\n")
h.write("\n")

# TORSIONS

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

for l in [
'X -C -J3-X    4   10.00        180.0             2.         WLD from X -C -N -X',
'X -DJ-N -X    4   10.00        180.0             2.         WLD from X -C -N -X',

'DA-DJ-DJ-Q1   4    0.00        180.0             2.         WLD on benzamide',
'DJ-DJ-DJ-Q1   4    0.00        180.0             2.         WLD on benzamide',

'D4-T2-S -CT   1    3.50          0.0            -2.         WLD from CT-S-S-CT',
'D4-T2-S -CT   1    0.60          0.0             3.         WLD from CT-S-S-CT',


'N -D4-DJ-J3   1    2.000       180.000           2.         WLD from N-CT-C -N ',
'DJ-CT-C -N    1    2.000       180.000           2.         WLD from N-CT-C -N ',

'C -J3-D4-DJ   1    0.850       180.000          -2.         WLD from C-N -CT-C ',
'C -J3-D4-DJ   1    0.800         0.000           1.         WLD from C-N -CT-C ',

'DJ-N -CT-C    1    0.850       180.000          -2.         WLD from C-N -CT-C ',
'DJ-N -CT-C    1    0.800         0.000           1.         WLD from C-N -CT-C ',

'D4-D4-J3-C    1    0.50        180.0            -4.         WLD from CT-CT-N -C',
'D4-D4-J3-C    1    0.15        180.0            -3.         WLD from CT-CT-N -C',
'D4-D4-J3-C    1    0.53          0.0             1.         WLD from CT-CT-N -C',

'CT-CT-N -DJ   1    0.50        180.0            -4.         WLD from CT-CT-N -C',
'CT-CT-N -DJ   1    0.15        180.0            -3.         WLD from CT-CT-N -C',
'CT-CT-N -DJ   1    0.53          0.0             1.         WLD from CT-CT-N -C',

'D4-D4-DJ-N    1    0.100         0.0            -4.         WLD from CT-CT-C -N',
'D4-D4-DJ-N    1    0.07          0.0             2.         WLD from CT-CT-C -N',

'CT-CT-C -J3   1    0.100         0.0            -4.         WLD from CT-CT-C -N',
'CT-CT-C -J3   1    0.07          0.0             2.         WLD from CT-CT-C -N',

'H -N -DJ-Q1   1    2.50        180.0            -2.         WLD from H -N -C -O',
'H -N -DJ-Q1   1    2.00          0.0             1.         WLD from H -N -C -O',

'A -J3-C -O    1    2.50        180.0            -2.         WLD from H -N -C -O',
'A -J3-C -O    1    2.00          0.0             1.         WLD from H -N -C -O',

]:
   h.write(l+"\n")
h.write("\n")

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

h.write("\n")

while 1:
   l = f.readline()
   h.write(l)
   if l[0:4]=='MOD4':
      break

while 1:
   l = g.readline()
   if l[0:4]=='MOD4':
      break

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

h.write("\n")

h.write(f.readline())
g.readline()
h.write(f.readline())
g.readline()
h.write(f.readline())
g.readline()

while 1:
   l = f.readline()
   if not l: break
   if l[0:6]=='TINKER':
      h.write(l)

while 1:
   l = g.readline()
   if not l: break
   if l[0:6]=='TINKER':
      h.write(l)

