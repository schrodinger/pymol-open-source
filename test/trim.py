import sys
import string

f = sys.stdin
g = sys.stdout

echo = 0
while 1:
   l = f.readline()
   if not l: break
   ll=string.strip(l)
   if ll=='BEGIN-LOG':
      echo = 1
   elif ll=='END-LOG':
      echo = 0
   elif echo:
      l=string.replace(l,"-0.000"," 0.000") # squish annoying negative zeros
      g.write(l)
