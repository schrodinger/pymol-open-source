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

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   h.write(l)
   if not len(string.strip(l)):
      break

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   h.write(l)
   if not len(string.strip(l)):
      break

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   h.write(l)
   if not len(string.strip(l)):
      break

while 1:
   l = f.readline()
   if not len(string.strip(l)):
      break
   h.write(l)

while 1:
   l = g.readline()
   h.write(l)
   if not len(string.strip(l)):
      break

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
   h.write(l)
   if not len(string.strip(l)):
      break

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

