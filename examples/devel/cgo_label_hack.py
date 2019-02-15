import math
from pymol.cgo import *
from pymol import cmd

axes = [
   LINEWIDTH, 3.0,
   BEGIN, LINES,
   COLOR, 0.2, 1.0, 0.2,
   
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 3.0, 0.0, 0.0,

   COLOR, 1.0, 0.2, 0.2,
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 0.0, 3.0, 0.0,

   COLOR, 0.2, 0.2, 1.0,   
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 00, 0.0, 3.0,
   END
   ]



c=0
for a in range(0,63):
   balls = [
      COLOR,  0.2, 1.0, 0.2,
      SPHERE, 1.0+math.cos(a/10.0), 1.0+math.sin(a/20.0), 1.0+math.cos(a/10.0), 0.2+math.cos(a/5.0)/10.0,

      COLOR,  1.0, 0.2, 0.2,      
      SPHERE, 2.0-math.cos(a/10.0), 1.0+math.sin(0.5+a/10.0), 1.0+math.cos(a/10.0), 0.2+math.cos(a/5.0)/10.0,      
      ]
   obj = axes + balls
   cmd.load_cgo(obj,'cgo01',c)
   c = c + 1
   
   # counter label
   
   pdb_list = [
      "HETATM%5d  C   UNK     1    %8.3f%8.3f%8.3f  1.00 10.00\n"%(c,2.0,0,2.0),
      ]
   cmd.read_pdbstr(''.join(pdb_list),'lab1',c,discrete=1)
   cmd.label("(lab1 and id %d)"%c,"'frame %d %6.3f'"%(c,math.sin(a/10.0)))


cmd.hide("nonbonded","lab1")

# axes labels


pdb_list = [
"HETATM    1  X   UNK     1    %8.3f%8.3f%8.3f  1.00 10.00\n"%(3.2,0,0),
"HETATM    2  Y   UNK     2    %8.3f%8.3f%8.3f  1.00 10.00\n"%(0,3.2,0),
"HETATM    3  Z   UNK     3    %8.3f%8.3f%8.3f  1.00 10.00\n"%(0,0,3.2),]
cmd.read_pdbstr(''.join(pdb_list),'lab2')
cmd.hide('(lab2)')
cmd.label('lab2','name')
cmd.color('white','lab2')


cmd.zoom('cgo01')
cmd.clip('far',-5)

cmd.mplay()

