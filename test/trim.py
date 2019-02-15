import sys

f = sys.stdin
g = sys.stdout

echo = 0
while 1:
   l = f.readline()
   if not l: break
   ll=l.strip()
   if ll=='BEGIN-LOG':
      echo = 1
   elif ll=='END-LOG':
      echo = 0
   elif echo:
      l=l.replace("-0.000"," 0.000") # squish annoying negative zeros
      g.write(l)
