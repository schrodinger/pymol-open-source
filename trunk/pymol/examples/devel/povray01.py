from pymol import cmd
import os

if not ('pept' in cmd.get_names()):
   cmd.delete('all')
   util.ray_shadows('heavy')
   cmd.do('load $PYMOL_PATH/test/dat/pept.pdb')
   cmd.do('set surface_quality=1')
   cmd.do('show surface;hide lines;')
   cmd.zoom('all',10)
   cmd.do('clip far,-40;show surface;hide lines;set smooth_color=1')
   cmd.viewport(300,300)
pov = cmd.get_povray()

f=open("tmp_pymol.pov",'w')
f.write(pov[0])
f.write("#include \"colors.inc\"\n");
f.write("#include \"stones.inc\"\n");
f.write("#include \"woods.inc\"\n");
f.write(pov[1])
f.write("plane {<0,1,0.5>, -70 texture{T_Grnt10 scale 40}}\n");
f.write("sphere {<14,4,-132>, 7.0 pigment{color Grey} finish{reflection 1.0 metallic}}\n");
for x in xrange(-30,30,5):
   f.write("sphere {<%6.4f,-12,-100>, 2.0 texture{P_WoodGrain1B scale 8}}\n"%(x))
f.write("box {<-17,17,-120>,<13,15,-160> pigment {color Red}}\n")
f.close()
cmd.refresh()
os.system("x-povray +Itmp_pymol.pov +Otmp_pymol.png +W300 +H300 +A")
cmd.load_png('tmp_pymol.png')

