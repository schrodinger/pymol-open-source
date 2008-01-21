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

polar_contacts_suffix = "_pol_conts"
default_polar_contacts = "polar_contacts"

tmp_sele = "_p_tmp"

prot_and_dna_sele = "(resn ALA+CYS+CYX+ASP+GLU+PHE+GLY+HIS+HID+HIE+HIP+HISE+HISD+HISP+ILE+LYS+LEU+MET+MSE+ASN+PRO+GLN+ARG+SER+THR+VAL+TRP+TYR+A+C+T+G+U)"
wat_sele = "(resn WAT,H2O,HOH,TIP)"
ion_sele = "(resn CA,HG,K,NA,ZN,MG,CL)"
solv_sele = "("+wat_sele+"|"+ion_sele+")"
lig_excl = "(resn MSE)"
lig_sele = "((hetatm or not "+prot_and_dna_sele+") and not ("+solv_sele+"|"+ion_sele+"|"+lig_excl+"))"
lig_and_solv_sele = "("+lig_sele+"|"+solv_sele+")"

def _get_polar_contacts_name(s,_self=cmd):
    cmd=_self
    list = cmd.get_object_list(s)
    if list!=None and (len(list)==1):
        return list[0]+polar_contacts_suffix
    else:
        return default_polar_contacts

def _prepare(selection,polar_contacts=None,_self=cmd):
    cmd=_self
    # this function should undo everything that is done by any preset function in this module
    # (except for coloring)
    s = tmp_sele
    cmd.select(s,selection)

    cmd.cartoon("auto",s)   
    cmd.hide("everything",s)
    
    cmd.set("two_sided_lighting",0) # global
    cmd.unset("transparency",s)
    cmd.unset("dot_normals",s)
    cmd.unset("mesh_normals",s)
    cmd.unset("surface_quality",s)
    cmd.unset("surface_type",selection)
    cmd.unset("sphere_scale",selection)
    cmd.unset_bond("stick_radius",s,s)
    cmd.unset_bond("stick_color",s,s)
    cmd.unset("cartoon_highlight_color",selection)
    cmd.unset("cartoon_fancy_helices",selection)
    cmd.unset("cartoon_smooth_loops",selection)
    cmd.unset("cartoon_flat_sheets",selection)
    cmd.unset("cartoon_side_chain_helper",selection)   
    cmd.unset("mesh_normals",s)
    cmd.unset("dot_normals",s)
    if polar_contacts == None:
        polar_contacts = _get_polar_contacts_name(s,_self)
        if polar_contacts in cmd.get_names('objects'):
            cmd.delete(polar_contacts)
        
def prepare(selection="(all)",_self=cmd):
    cmd=_self
    s = tmp_sele
    cmd.select(s,selection)
    _prepare(s,_self=cmd)

def simple(selection="(all)",_self=cmd):
    cmd=_self
    s = tmp_sele
    cmd.select(s,selection)
    _prepare(s,_self=cmd)
    util.cbc(s,_self=cmd)
    cmd.show("ribbon",s)
    cmd.show("lines","(byres (("+s+" & r. CYS+CYX & n. SG) & bound_to ("+s+" & r. CYS+CYX & n. SG))) & n. CA+CB+SG")
    # try to show what covalent ligands are connected to...
    cmd.show("sticks","("+lig_sele+" and ("+s+")) extend 2")
    cmd.show("sticks","byres (("+lig_sele+" and ("+s+") and not resn ACE+NAC+NME+NH2) extend 1)")
    cmd.hide("sticks","("+s+") and ((not rep sticks) extend 1)")
    cmd.show("sticks","("+lig_sele+" and ("+s+")) extend 2")
    # color by atom if lines or sticks are shown
    util.cnc("(( rep lines or rep sticks or ("+lig_and_solv_sele+")) and ("+s+"))",_self=cmd)
    cmd.show("nonbonded","("+lig_and_solv_sele+" and ("+s+"))")
    cmd.show("lines","("+lig_and_solv_sele+" and ("+s+"))")
    if cmd.count_atoms(s):
        cmd.zoom(s)
    cmd.delete(s)

def simple_no_solv(selection="(all)",_self=cmd):
    cmd=_self
    simple(selection,_self=_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.hide("nonbonded","("+solv_sele+" and "+s+")")

def ligands(selection="(all)",_self=cmd):
    cmd=_self
    try:
        s = tmp_sele
        cmd.select(s,selection)
        polar_contacts = _get_polar_contacts_name(s,_self)
        _prepare(s,polar_contacts,_self=cmd)
        host = "_preset_host"
        solvent = "_preset_solvent"
        near_solvent = "_preset_solvent"
        lig = "_preset_lig"
        cmd.select(host,s+" and "+prot_and_dna_sele)
        cmd.select(solvent,s+" and "+solv_sele)
        cmd.select(lig,s+" and "+lig_sele)
        cmd.select(near_solvent,s+" and ("+solvent+" within 4 of "+lig+")")

        util.chainbow(host,_self=cmd)
        util.cbc(lig,_self=cmd)
        util.cbac("(("+s+") and not elem c)",_self=cmd)
        cmd.hide("everything",s)
        cmd.show("ribbon",host)
        cmd.show("lines","("+s+" and byres ("+host+" within 5 of "+lig+"))")
        cmd.show("sticks",lig)
        cmd.show("sticks",solvent+" and neighbor "+lig)
        cmd.show("lines","("+s+" and (rep lines extend 1) and "+lig+")")

        if cmd.count_atoms(lig):
            cmd.dist(polar_contacts,host+"|"+near_solvent,lig+"|"+near_solvent,
                     mode=2,quiet=1,label=0,reset=1) # hbonds
            if polar_contacts in cmd.get_names():
                cmd.enable(polar_contacts)
                cmd.hide("labels",polar_contacts)
                cmd.show("dashes",polar_contacts)
        else:
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

def ball_and_stick(selection="(all)",_self=cmd):
    cmd=_self
    s = tmp_sele
    cmd.select(s,selection)
    _prepare(s,_self=cmd)
    cmd.hide("everything",s)
    cmd.set_bond("stick_color","white",s,s)
    cmd.set_bond("stick_radius","0.14",s,s)
    cmd.set("sphere_scale","0.25",s)
    cmd.show("sticks",s)
    cmd.show("spheres",s)
    
def b_factor_putty(selection="(name ca or name p)",_self=cmd):
    cmd=_self
    s = tmp_sele
    cmd.select(s,selection)
    _prepare(s,_self=cmd)
    cmd.select(s,"(name ca or name p) and ("+selection+")")
    cmd.show("cartoon",s)
    cmd.set("cartoon_flat_sheets",0,selection)
    cmd.cartoon("putty",s)
    cmd.spectrum("b",selection=s)

def ligand_cartoon(selection="(all)",_self=cmd):
    cmd=_self
    ligand_sites(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.set("cartoon_side_chain_helper",1,selection)
    cmd.show("cartoon","rep ribbon")
    cmd.hide("ribbon")
    cmd.hide("surface")
    
def ligand_sites(selection="(all)",_self=cmd):
    cmd=_self
    try:
        s = tmp_sele
        cmd.select(s,selection)
        polar_contacts = _get_polar_contacts_name(s,_self)
        _prepare(s,polar_contacts,_self=cmd)
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

        util.chainbow(host,_self=cmd)
        util.cbc(lig,_self=cmd)
        util.cbac("(("+s+") and not elem c)",_self=cmd)
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
            cmd.dist(polar_contacts,host+"|"+near_solvent,lig+"|"+near_solvent,mode=2,quiet=1,label=0,reset=1) # hbonds
            if polar_contacts in cmd.get_names():
                cmd.enable(polar_contacts)
                cmd.hide("labels",polar_contacts)
                cmd.show("dashes",polar_contacts)
        else:
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

def ligand_sites_hq(selection="(all)",_self=cmd):
    cmd=_self
    ligand_sites(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.set("surface_quality","1",selection)
    cmd.set("surface_type",0,selection)

def ligand_sites_trans(selection="(all)",_self=cmd):
    cmd=_self
    ligand_sites(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("transparency","0.33",s)
    cmd.set("surface_type",0,selection)
    cmd.set("surface_quality",0,selection)

def ligand_sites_trans_hq(selection="(all)",_self=cmd):
    cmd=_self
    ligand_sites(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("transparency","0.33",s)
    cmd.set("surface_type",0,selection)
    cmd.set("surface_quality",1,selection)
    
def ligand_sites_mesh(selection="(all)",_self=cmd):
    cmd=_self
    ligand_sites(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("surface_type","2",selection)
    cmd.set("surface_quality","0",selection)
    cmd.set("mesh_normals",0,s)
    
def ligand_sites_dots(selection="(all)",_self=cmd):
    cmd=_self
    ligand_sites(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("surface_type","1",selection)
    cmd.set("surface_quality","1",selection)
    cmd.set("dot_normals",0,s)

def technical(selection="(all)",_self=cmd):
    cmd=_self
    s = tmp_sele
    cmd.select(s,selection)
    polar_contacts = _get_polar_contacts_name(s,_self)
    _prepare(s,polar_contacts,_self=cmd)
    util.chainbow(s,_self=cmd)
    util.cbc("("+lig_sele+" and ("+s+"))",_self=cmd)   
    util.cbac("(("+s+") and not elem c)",_self=cmd)
    cmd.show("nonbonded",s)
    cmd.show("lines","((("+s+") and not "+lig_sele+") extend 1)")
    cmd.show("sticks","("+lig_sele+" and ("+s+"))")
    cmd.show("ribbon",s)
    cmd.dist(polar_contacts,s,s,mode=2,label=0,reset=1) # hbonds
    if polar_contacts in cmd.get_names():
        cmd.enable(polar_contacts)
        cmd.set("dash_width",1.5,polar_contacts)
        cmd.hide("labels",polar_contacts)
        cmd.show("dashes",polar_contacts)
    cmd.show("nonbonded","(("+lig_sele+"|resn hoh+wat+h2o) and ("+s+"))")

def pretty_solv(selection="(all)",_self=cmd):
    cmd=_self
    s = tmp_sele
    cmd.select(s,selection)
    polar_contacts = _get_polar_contacts_name(s,_self)
    _prepare(s,polar_contacts,_self=cmd)
    cmd.dss(s,preserve=1)
    cmd.cartoon("auto",s)
    cmd.show("cartoon",s)
    cmd.show("sticks","("+lig_sele+" and ("+s+"))")
    cmd.show("nb_spheres","(("+lig_sele+"|resn hoh+wat+h2o) and ("+s+"))")
    util.cbc("("+lig_sele+" and ("+s+"))",_self=cmd)
    util.cbac("("+lig_sele+" and ("+s+") and not elem c)",_self=cmd)
    cmd.spectrum("count",selection="(elem c and ("+s+") and not "+lig_sele+")")
    cmd.set("cartoon_highlight_color",-1,selection)
    cmd.set("cartoon_fancy_helices",0,selection)
    cmd.set("cartoon_smooth_loops",0,selection)
    cmd.set("cartoon_flat_sheets",1,selection)
    cmd.set("cartoon_side_chain_helper",0,selection)   
    if polar_contacts in cmd.get_names():
        cmd.disable(polar_contacts)
    if cmd.count_atoms(s):
        cmd.zoom(s)
        
def pretty(selection,_self=cmd):
    cmd=_self
    pretty_solv(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.hide("nb_spheres","("+s+" and "+lig_sele+"|resn hoh+wat+h2o)")

pretty_no_solv = pretty

def pub_solv(selection="(all)",_self=cmd):
    cmd=_self
    pretty_solv(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.set("cartoon_smooth_loops",1,selection)
    cmd.set("cartoon_highlight_color","grey50",selection)
    cmd.set("cartoon_fancy_helices",1,selection)
    cmd.set("cartoon_flat_sheets",1,selection)
    cmd.set("cartoon_side_chain_helper",0,selection)   
    if cmd.count_atoms(s):
        cmd.zoom(s)

def publication(selection="(all)",_self=cmd):
    cmd=_self
    pub_solv(selection,_self)
    s = tmp_sele
    cmd.select(s,selection)
    cmd.hide("nb_spheres","(("+lig_sele+"|resn hoh+wat+h2o) and "+s+")")

pub_no_solv = publication

def default(selection="(all)",_self=cmd):
    cmd=_self
    s = tmp_sele
    cmd.select(s,selection)
    _prepare(s,_self=cmd)
    cmd.show("lines",s)
    cmd.show("nonbonded",s)
    color=cmd.get_object_color_index(selection)
    if color<0:
        util.cbag(selection,_self=cmd)
    else:
        util.cnc(selection,_self=cmd)
        cmd.color(str(color),"("+s+") and elem c")
        
