from pymol.vfont import plain
from pymol.cgo import *
import string
import traceback

from pymol import cmd
from pymol import util

class Demo: # stateful class for doing effective demonstrations

   def __init__(self):
      self.last = None
      
   def __call__(self,name):
      if self.last:
         self.last(cleanup=1)
      if hasattr(self,name):
         self.last = getattr(self,name)
         self.last()
                 
   def rep_old(self,cleanup=0):
      if not cleanup:
         cmd.set("suspend_updates",1,quiet=1)
         cmd.disable()
         cmd.do("cd $PYMOL_PATH")
         cmd.delete("pept")
         cmd.delete("pept_dist")
         cmd.load("test/dat/pept.pdb")
         cmd.show("sticks","(pept and not i;5:7)")
         cmd.show("surface","(pept and i;5,6)")
         cmd.show("mesh","(pept and i;1,11,12,13)")
         cmd.show("spheres","(pept and i;2,12,9,4 and not n;c,n,o,ca)")
         cmd.show("dots","(i;8)")
         cmd.dist("pept_dist","(pept and i;1&n;OD2)","(pept and i;13&n;OG1)")
         cmd.set("dot_width","2");
         cmd.set("suspend_updates",0,quiet=1)
      else:
         cmd.delete("pept")
         cmd.delete("pept_dist")

   def rep(self,cleanup=0):
      rep_list = [ "lines","sticks","spheres","surface","mesh","dots","ribbon","cartoon" ]
      try:
         if not cleanup:
            cmd.set("suspend_updates",1,quiet=1)
            cmd.load("$PYMOL_PATH/test/dat/pept.pdb","rep1")
            cmd.alter("rep1///1-5+8-13/","ss='S'")
            cmd.cartoon("auto")
            cmd.hide("everything","rep1")
            for a in range(2,9):
               cmd.create("rep%d"%a,"rep1")
            map(lambda x,y:cmd.show(x,"rep%d"%y),
                rep_list,
                range(1,9))
            cmd.reset()
            cmd.zoom("rep1",24)
            util.cbay("rep2")
            util.cbac("rep3")
            util.cbas("rep4")
            util.cbab("rep5")
            util.cbaw("rep6")            
            util.cbay("rep8")
            
            cmd.set("suspend_updates",0,quiet=1)
            scale=0.5
            for b in range(1,20):
               cmd.set("suspend_updates",0,quiet=1)
               cmd.refresh()
               cmd.set("suspend_updates",1,quiet=1)
               xt=-3.2
               yt=1.6
               for a in range(1,5):
                  cmd.translate([xt*scale,yt*scale,0],object="rep%d"%a,camera=0)
                  xt=xt+2
               yt=-yt
               xt=-3.2
               for a in range(5,9):
                  cmd.translate([xt*scale,yt*scale,0],object="rep%d"%a,camera=0)
                  xt=xt+2
            for a in range(1,9):
               cmd.origin("rep%d"%a,object="rep%d"%a)
            cmd.mset("1")
            st = string.join(map(lambda x,y:"rotate angle=3,object=rep%d,axis=%s;"%(x,y),range(1,9),
                                 ['x','y','x','y','x','y','x','y']))
            cmd.mdo(1,st)
            cmd.set("suspend_updates",0,quiet=1)
            cmd.mplay()

            cgo = []
            axes = [[4.5,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]]

            c = 1
            for a in rep_list:
               ext = cmd.get_extent("rep%d"%c)
               pos = [(ext[0][0]+ext[1][0])/2,
                      (ext[0][1]+ext[1][1])/2+14,
                      (ext[0][2]+ext[1][2])/2]
               c = c + 1
               pos[0]=pos[0]-(measure_text(plain,a,axes)/2)
               wire_text(cgo,plain,pos,a,axes)
            cmd.set("cgo_line_width",1.5)
            cmd.set("auto_zoom",0)
            cmd.load_cgo(cgo,'reps')
            cmd.set("auto_zoom",1)
         else:
            cmd.delete("rep*")
            cmd.mset()
            cmd.mstop()
      except:
         traceback.print_exc()
         
   def raster3d(self,cleanup=0):
      if not cleanup:
         cmd.disable()
         cmd.do("cd $PYMOL_PATH")
         cmd.delete("cgo1")
         cmd.delete("cgo2")
         cmd.do("cd $PYMOL_PATH")
         cmd.load("test/dat/pept.r3d","cgo1")
         cmd.load("test/dat/3al1.r3d","cgo2")
         cmd.zoom()
      else:
         cmd.delete("cgo1")
         cmd.delete("cgo2")

   def cgo(self,cleanup=0):
      if not cleanup:
         cmd.disable()
         cmd.do("cd $PYMOL_PATH")
         cmd.do("run examples/devel/cgo03.py")
      else:
         cmd.delete("cgo03")
         cmd.mset()
         cmd.mstop()

   def anime(self,cleanup=0):
      if not cleanup:
         cmd.disable()
         cmd.delete("arg")
         cmd.fragment("arg")
         cmd.zoom("arg",2)
         cmd.show("sticks","arg")
         cmd.feedback('dis','sel','res')
         for a in xrange(1,181):
            cmd.set("suspend_updates",1,quiet=1)
            cmd.edit("(arg and n;cd)","(arg and n;cg)")
            cmd.torsion("6")
            cmd.unpick()
            cmd.edit("(arg and n;cb)","(arg and n;ca)")
            cmd.torsion("2")
            cmd.unpick()
            cmd.set("suspend_updates",0,quiet=1)         
            cmd.refresh()
         cmd.feedback('ena','sel','res')
      else:
         cmd.delete("arg")

   def cartoon(self,cleanup=0):
      if not cleanup:
         cmd.set("suspend_updates",1,quiet=1)
         cmd.disable()
         cmd.delete("1tii")      
         cmd.load("$PYMOL_PATH/test/dat/1tii.pdb")
         cmd.hide("(1tii)")
         cmd.show("cartoon","1tii")
         cmd.zoom("1tii")
         util.color_chains("1tii")
         cmd.set("suspend_updates",0,quiet=1)
         cmd.refresh()
      else:
         cmd.delete("1tii")

   def trans(self,cleanup=0):
      if not cleanup:
         cmd.set("suspend_updates",1,quiet=1)
         cmd.disable()
         cmd.delete("trans")
         cmd.load("$PYMOL_PATH/test/dat/pept.pdb","trans")
         cmd.hide("(tran)")
         cmd.show("surface","trans")
         cmd.show("sticks","trans")
         cmd.set("surface_color","white","trans")
         cmd.set("transparency",0.5,"trans")
         cmd.zoom("trans")
         cmd.set("suspend_updates",0,quiet=1)
         cmd.refresh()
      else:
         cmd.delete("trans")

   def ray(self,cleanup=0):
      if not cleanup:
         cmd.set("suspend_updates",1,quiet=1)
         cmd.disable()
         cmd.delete("ray")
         cmd.load("$PYMOL_PATH/test/dat/il2.pdb","ray")
         cmd.remove("(ray and hydro)")
         cmd.hide("lines","ray")
         cmd.show("spheres","ray")
         cmd.orient("ray")
         cmd.turn("x",90)
         util.ray_shadows('heavy')
         cmd.set("suspend_updates",0,quiet=1)
         cmd.refresh()
         cmd.do("ray")
      else:
         cmd.delete("ray")

   def sculpt(self,cleanup=0):
      if not cleanup:
         cmd.set("suspend_updates",1,quiet=1)
         cmd.disable()
         cmd.delete("sculpt")
         cmd.load("$PYMOL_PATH/test/dat/pept.pdb","sculpt")
         cmd.show("spheres","sculpt")
         cmd.set("sphere_scale","0.5","sculpt")
         cmd.set("sphere_transparency","0.5","sculpt")
         cmd.set("sphere_color","grey","sculpt")
         cmd.sculpt_activate("sculpt")
         cmd.set("auto_sculpt",1)
         cmd.set("sculpting",1)
         cmd.do("edit_mode")
         cmd.set("valence","0.05")
         cmd.set("suspend_updates",0,quiet=0)
      else:
         cmd.set("sculpting",0)
         cmd.set("auto_sculpt",0)
         cmd.delete("sculpt")
