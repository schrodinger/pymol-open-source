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

polar_contacts = "p_polar_contacts"
tmp_sele = "_p_tmp"

prot_and_dna_sele = "(resn ALA+CYS+CYX+ASP+GLU+PHE+GLY+HIS+ILE+LYS+LEU+MET+MSE+ASN+PRO+GLN+ARG+SER+THR+VAL+TRP+TYR+A+C+T+G+U)"
wat_sele = "(resn WAT,H2O,HOH,TIP)"
ion_sele = "(resn CA,HG,K,NA,ZN,MG,CL)"
solv_sele = "("+wat_sele+"|"+ion_sele+")"
lig_excl = "(resn MSE)"
lig_sele = "((hetatm or not "+prot_and_dna_sele+") and not ("+solv_sele+"|"+ion_sele+"|"+lig_excl+"))"
lig_and_solv_sele = "("+lig_sele+"|"+solv_sele+")"

def simple(selection="(all)"):
   s = tmp_sele
   cmd.select(s,selection)
   util.cbc(s)
   cmd.hide("everything",s)
   cmd.show("ribbon",s)
   cmd.show("lines","(byres (("+s+" & r. CYS+CYX & n. SG) & bound_to ("+s+" & r. CYS+CYX & n. SG))) & n. CA+CB+SG")
   cmd.show("sticks","("+lig_sele+" and ("+s+"))")
   util.cnc("("+lig_and_solv_sele+" and ("+s+"))")
   cmd.show("nonbonded","("+lig_and_solv_sele+" and ("+s+"))")
   if polar_contacts in cmd.get_names():
      cmd.disable(polar_contacts)
   if cmd.count_atoms(s):
      cmd.zoom(s)
   cmd.delete(s)

def simple_no_solv(selection="(all)"):
   simple(selection)
   s = tmp_sele
   cmd.select(s,selection)
   cmd.hide("nonbonded","("+solv_sele+" and "+s+")")

def ligands(selection="(all)"):
   try:
      s = tmp_sele
      cmd.select(s,selection)
      host = "_preset_host"
      solvent = "_preset_solvent"
      near_solvent = "_preset_solvent"
      lig = "_preset_lig"
      cmd.select(host,s+" and "+prot_and_dna_sele)
      cmd.select(solvent,s+" and "+solv_sele)
      cmd.select(lig,s+" and "+lig_sele)
      cmd.select(near_solvent,s+" and ("+solvent+" within 4 of "+lig+")")

      util.chainbow(host)
      util.cbc(lig)
      util.cbac("(("+s+") and not elem c)")
      cmd.hide("everything",s)
      cmd.show("ribbon",host)
      cmd.show("lines","("+s+" and byres ("+host+" within 5 of "+lig+"))")
      cmd.show("sticks",lig)
      cmd.show("sticks",solvent+" and neighbor "+lig)
      cmd.show("lines","("+s+" and (rep lines extend 1) and "+lig+")")

      if cmd.count_atoms(lig):
         cmd.dist(polar_contacts,host+"|"+near_solvent,lig+"|"+near_solvent,mode=2,quiet=1,labels=0) # hbonds
         if polar_contacts in cmd.get_names():
            cmd.enable(polar_contacts)
            cmd.hide("labels",polar_contacts)
            cmd.show("dashes",polar_contacts)            
      elif polar_contacts in cmd.get_names():
         cmd.delete(polar_contacts)
      cmd.show("nonbonded",lig+"|"+host+"|"+near_solvent)
      if cmd.count_atoms(lig):
         cmd.zoom(lig,3)
      cmd.delete(host)
      cmd.delete(solvent)
      cmd.delete(near_solvent)
      cmd.delete(lig)
   except:
      traceback.print_exc()
      
def ligand_sites(selection="(all)"):
   try:
      s = tmp_sele
      cmd.select(s,selection)
      host = "_preset_host"
      solvent = "_preset_solvent"
      near_solvent = "_preset_solvent"
      lig = "_preset_lig"
      cmd.select(host,s+" and "+prot_and_dna_sele)
      cmd.select(solvent,s+" and "+solv_sele)
      cmd.select(lig,s+" and "+lig_sele)
      cmd.select(near_solvent,s+" and ("+solvent+" within 4 of "+lig+")")
      cmd.flag("ignore",host,"clear")
      cmd.flag("ignore",lig+"|"+solvent,"set")

      util.chainbow(host)
      util.cbc(lig)
      util.cbac("(("+s+") and not elem c)")
      cmd.hide("everything",s)
      cmd.show("ribbon",host)
      cmd.show("lines","("+s+" and byres ("+host+" within 5 of "+lig+"))")
      cmd.show("surface","("+s+" and ((rep lines expand 4) within 6 of "+lig+"))")
      cmd.set("two_sided_lighting",1) # global setting
      cmd.set("transparency",0,s)
      cmd.set("surface_quality",0,s)

      cmd.show("sticks",lig)
      cmd.show("sticks",solvent+" and neighbor "+lig)
      cmd.show("lines","("+s+" and (rep lines extend 1) and "+lig+")")

      if cmd.count_atoms(lig):
         cmd.dist(polar_contacts,host+"|"+near_solvent,lig+"|"+near_solvent,mode=2,quiet=1,labels=0) # hbonds
         if polar_contacts in cmd.get_names():
            cmd.enable(polar_contacts)
            cmd.hide("labels",polar_contacts)
            cmd.show("dashes",polar_contacts)            
      elif polar_contacts in cmd.get_names():
         cmd.delete(polar_contacts)
      cmd.show("nb_spheres",lig+"|"+host+"|"+near_solvent)
      if cmd.count_atoms(lig):
         cmd.zoom(lig,3)
      cmd.delete(host)
      cmd.delete(solvent)
      cmd.delete(near_solvent)
      cmd.delete(lig)
   except:
      traceback.print_exc()

def ligand_sites_hq(selection="(all)"):
   ligand_sites("all")
   s = tmp_sele
   cmd.select(s,selection)
   cmd.show("sticks",s+" and rep lines")
   cmd.hide("lines",s+" and rep lines")
   cmd.set("surface_quality","1",s)
   cmd.set("surface_type",0,s)

def ligand_sites_trans(selection="(all)"):
   ligand_sites("all")
   s = tmp_sele
   cmd.select(s,selection)
   cmd.show("sticks",s+" and rep lines")
   cmd.hide("lines",s+" and rep lines")
   cmd.set("transparency","0.33",s)
   cmd.set("surface_type",0,s)

def ligand_sites_mesh(selection="(all)"):
   ligand_sites("all")
   s = tmp_sele
   cmd.select(s,selection)
   cmd.show("sticks",s+" and rep lines")
   cmd.hide("lines",s+" and rep lines")
   cmd.set("surface_type","2",s)
   cmd.set("surface_quality","0",s)
   cmd.set("mesh_normals",0)

def ligand_sites_dots(selection="(all)"):
   ligand_sites("all")
   s = tmp_sele
   cmd.select(s,selection)
   cmd.show("sticks",s+" and rep lines")
   cmd.hide("lines",s+" and rep lines")
   cmd.set("surface_type","1",s)
   cmd.set("surface_quality","1",s)
   cmd.set("dot_normals",0)

def technical(selection="(all)"):
   s = tmp_sele
   cmd.select(s,selection)
   util.chainbow(s)
   util.cbc("("+lig_sele+" and ("+s+"))")   
   util.cbac("(("+s+") and not elem c)")
   cmd.hide("everything",s)
   cmd.show("nonbonded",s)
   cmd.show("lines","((("+s+") and not "+lig_sele+") extend 1)")
   cmd.show("sticks","("+lig_sele+" and ("+s+"))")
   cmd.show("ribbon",s)
   cmd.dist(polar_contacts,s,s,mode=2,labels=0) # hbonds
   if polar_contacts in cmd.get_names():
      cmd.enable(polar_contacts)
      cmd.set("dash_width",1.5,polar_contacts)
      cmd.hide("labels",polar_contacts)
      cmd.show("dashes",polar_contacts)
   cmd.show("nonbonded","(("+lig_sele+"|resn hoh+wat+h2o) and ("+s+"))")

def pretty_solv(selection="(all)"):
   s = tmp_sele
   cmd.select(s,selection)
   cmd.dss(s,preserve=1)
   cmd.hide("everything",s)
   cmd.show("cartoon",s)
   cmd.show("sticks","("+lig_sele+" and ("+s+"))")
   cmd.show("nb_spheres","(("+lig_sele+"|resn hoh+wat+h2o) and ("+s+"))")
   util.cbc("("+lig_sele+" and ("+s+"))")
   util.cbac("("+lig_sele+" and ("+s+") and not elem c)")
   cmd.spectrum("count",selection="(elem c and ("+s+") and not "+lig_sele+")")
   cmd.set("cartoon_highlight_color",-1,s)
   cmd.set("cartoon_fancy_helices",0,s)
   cmd.set("cartoon_smooth_loops",0,s)
   if polar_contacts in cmd.get_names():
      cmd.disable(polar_contacts)
   if cmd.count_atoms(s):
      cmd.zoom(s)
      
def pretty(selection):
   pretty_solv(selection)
   s = tmp_sele
   cmd.select(s,selection)
   cmd.hide("nb_spheres","("+s+" and "+lig_sele+"|resn hoh+wat+h2o)")

pretty_no_solv = pretty

def pub_solv(selection="(all)"):
   pretty_solv(selection)
   s = tmp_sele
   cmd.select(s,selection)
   cmd.set("cartoon_smooth_loops",1,s)
   cmd.set("cartoon_highlight_color","grey50",s)
   cmd.set("cartoon_fancy_helices",1,s)
   if polar_contacts in cmd.get_names():
      cmd.disable(polar_contacts)
   if cmd.count_atoms(s):
      cmd.zoom(s)

def publication(selection="(all)"):
   pub_solv(selection)
   s = tmp_sele
   cmd.select(s,selection)
   cmd.hide("nb_spheres","(("+lig_sele+"|resn hoh+wat+h2o) and "+s+")")

pub_no_solv = publication
   
def default(selection="(all)"):
   s = tmp_sele
   cmd.select(s,selection)
   cmd.hide("everything",s)
   cmd.show("lines",s)
   cmd.show("nonbonded",s)
   cmd.unset("transparency",s)
   cmd.set("two_sided_lighting",0)
   cmd.unset("dot_normals",s)
   cmd.unset("mesh_normals",s)
   cmd.unset("surface_quality",s)
   cmd.unset("surface_type",s)
   if polar_contacts in cmd.get_names():
      cmd.disable(polar_contacts)
   color=cmd.get_object_color_index(selection)
   if color<0:
      util.cbag(selection)
   else:
      util.cnc(selection)
      cmd.color(str(color),"("+s+") and elem c")
      
