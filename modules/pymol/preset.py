#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

import cmd
import util
import traceback

polar_contacts = "preset_polar_conts"

def simple(selection="(all)"):
   s = "("+str(selection)+")"
   util.cbc(s)
   cmd.hide("everything",s)
   cmd.show("ribbon",s)
   cmd.show("lines","(byres (("+s+" & r. CYS+CYX & n. SG) & bound_to ("+s+" & r. CYS+CYX & n. SG))) & n. CA+CB+SG")
   cmd.show("sticks","(hetatm and ("+s+"))")
   util.cnc("(hetatm and ("+s+"))")
   cmd.show("nonbonded","(hetatm and ("+s+"))")
   if polar_contacts in cmd.get_names():
      cmd.disable(polar_contacts)
   if cmd.count_atoms(s):
      cmd.center(s)

def simple_no_solv(selection="(all)"):
   simple(selection)
   s = "("+str(selection)+")"   
   cmd.hide("nonbonded","("+s+" and (hetatm|hoh+wat+h2o/))")

def ligands(selection="(all)"):
   try:
      s = "("+str(selection)+")"
      host = "_preset_host"
      water = "_preset_water"
      near_water = "_preset_water"
      lig = "_preset_lig"
      cmd.select(host,s+
        " and resn ALA+CYS+ASP+GLU+PHE+GLY+HIS+ILE+LYS+LEU+MET+MSE+ASN+PRO+GLN+ARG+SER+THR+VAL+TRP+TYR+A+C+T+G+U")
      cmd.select(water,s+" and resn WAT+HOH+H20")
      cmd.select(lig,s+" and not ("+host+"|"+water+")")
      cmd.select(near_water,s+" and ("+water+" within 5 of "+lig+")")

      util.chainbow(host)
      util.cbc(lig)
      util.cbac("(("+s+") and not elem c)")
      cmd.hide("everything",s)
      cmd.show("ribbon",host)
      cmd.show("lines","(byres ("+host+" within 5 of "+lig+"))")
      cmd.show("sticks",lig)

      if cmd.count_atoms(lig):
         cmd.dist(polar_contacts,host+"|"+near_water,lig,mode=2,quiet=1,labels=0) # hbonds
         if polar_contacts in cmd.get_names():
            cmd.enable(polar_contacts)
            cmd.hide("labels",polar_contacts)
            cmd.show("dashes",polar_contacts)            
      elif polar_contacts in cmd.get_names():
         cmd.delete(polar_contacts)
      cmd.show("nonbonded",lig+"|"+host+"|"+near_water)
      if cmd.count_atoms(lig):
         cmd.center(lig)
      cmd.delete(host)
      cmd.delete(water)
      cmd.delete(near_water)
      cmd.delete(lig)
   except:
      traceback.print_exc()
      
def technical(selection="(all)"):
   s = str(selection)
   util.chainbow(s)
   util.cbc("(hetatm and ("+s+"))")   
   util.cbac("(("+s+") and not elem c)")
   cmd.hide("everything",s)
   cmd.show("nonbonded",s)
   cmd.show("lines","((("+s+") and not hetatm) extend 1)")
   cmd.show("sticks","(hetatm and ("+s+"))")
   cmd.show("ribbon",s)
   cmd.dist(polar_contacts,s,s,mode=2,labels=0) # hbonds
   if polar_contacts in cmd.get_names():
      cmd.enable(polar_contacts)
      cmd.set("dash_width",1.5,polar_contacts)
      cmd.hide("labels",polar_contacts)
      cmd.show("dashes",polar_contacts)
   cmd.show("nonbonded","((hetatm|hoh+wat+h2o/) and ("+s+"))")

def pretty(selection="(all)"):
   s = str(selection)
   cmd.dss(s,preserve=1)
   cmd.hide("everything",s)
   cmd.show("cartoon",s)
   cmd.show("sticks","(hetatm and ("+s+"))")
   cmd.show("nb_spheres","((hetatm|hoh+wat+h2o/) and ("+s+"))")
   util.cbc("(hetatm and ("+s+"))")
   util.cbac("(hetatm and ("+s+") and not elem c)")
   cmd.spectrum("count",selection="(elem c and ("+s+") and not hetatm)")
   cmd.set("cartoon_highlight_color",-1)
   cmd.set("cartoon_fancy_helices",0)
   cmd.set("cartoon_smooth_loops",0)
   if polar_contacts in cmd.get_names():
      cmd.disable(polar_contacts)
   if cmd.count_atoms(s):
      cmd.center(s)
      
def pretty_no_solv(selection):
   pretty(selection)
   cmd.hide("nb_spheres","(hetatm|hoh+wat+h2o/)")
   
def publication(selection="(all)"):
   s = str(selection)
   pretty(s)
   cmd.set("cartoon_smooth_loops",1)
   cmd.set("cartoon_highlight_color","grey50")
   cmd.set("cartoon_fancy_helices",1)
   if polar_contacts in cmd.get_names():
      cmd.disable(polar_contacts)
   if cmd.count_atoms(s):
      cmd.center(s)

def pub_no_solv(selection="(all)"):
   publication(selection)
   s = "("+str(selection)+")"
   cmd.hide("nb_spheres","((hetatm|hoh+wat+h2o/) and "+s+")")
   
def default(selection="(all)"):
   s = "("+str(selection)+")"
   cmd.hide("everything",s)
   cmd.show("lines",s)
   cmd.show("nonbonded",s)
   color=cmd.get_object_color_index(selection)
   if color<0:
      util.cbag(selection)
   else:
      util.cnc(selection)
      cmd.color(str(color),"("+s+") and elem c")
      
