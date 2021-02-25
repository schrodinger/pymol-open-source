#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

cmd = __import__("sys").modules["pymol.cmd"]
from . import util
import traceback

polar_contacts_suffix = "_pol_conts"
default_polar_contacts = "polar_contacts"

tmp_sele = "_p_tmp"

prot_and_dna_sele = "(resn ALA+CYS+CYX+ASP+GLU+PHE+GLY+HIS+HID+HIE+HIP+HISE+HISD+HISP+ILE+LYS+LEU+MET+MSE+ASN+PRO+GLN+ARG+SER+THR+VAL+TRP+TYR+A+C+T+G+U+DA+DC+DT+DG+DU+DI)"
wat_sele = "solvent"
ion_sele = "(resn CA,HG,K,NA,ZN,MG,CL)"
solv_sele = "("+wat_sele+"|"+ion_sele+")"
lig_excl = "(resn MSE)"
lig_sele = "((hetatm or not "+prot_and_dna_sele+") and not ("+solv_sele+"|"+ion_sele+"|"+lig_excl+"))"
lig_and_solv_sele = "("+lig_sele+"|"+solv_sele+")"

def get_sname_oname_dname(selection, _self=cmd):
    '''
    Get a named selection, object names and a distance (polar contacts) name
    '''
    _self.select(tmp_sele, selection)

    if selection not in _self.get_object_list():
        selection = ' '.join(_self.get_object_list(tmp_sele))

    if selection and ' ' not in selection:
        dname = selection + polar_contacts_suffix
    else:
        dname = default_polar_contacts

    return tmp_sele, selection, dname

def _prepare(selection,polar_contacts=None,_self=cmd):
    cmd=_self
    # this function should undo everything that is done by any preset function in this module
    # (except for coloring)

    s, selection, dname = get_sname_oname_dname(selection, _self=cmd)

    cmd.cartoon("auto",s)
    cmd.hide("everything",s)

    cmd.set("two_sided_lighting",0) # global
    cmd.unset("transparency",s)
    cmd.unset("surface_quality", selection)
    cmd.unset("surface_type",selection)
    cmd.unset("sphere_scale",s)
    cmd.unset_bond("stick_radius",s,s)
    cmd.unset_bond("stick_color",s,s)
    cmd.unset("cartoon_highlight_color",selection)
    cmd.unset("cartoon_fancy_helices",selection)
    cmd.unset("cartoon_smooth_loops",selection)
    cmd.unset("cartoon_flat_sheets",selection)
    cmd.unset("cartoon_side_chain_helper",selection)
    if polar_contacts is None:
        polar_contacts = dname
        if polar_contacts in cmd.get_names('objects'):
            cmd.delete(polar_contacts)

    return s, selection, dname

def simple(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = _prepare(selection, _self=cmd)[:2]
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
    cmd.delete(s)

def simple_no_solv(selection="(all)",_self=cmd):
    cmd=_self
    simple(selection,_self=_self)
    s, selection = get_sname_oname_dname(selection, _self=_self)[:2]
    cmd.hide("everything","("+solv_sele+" and "+s+")")
    cmd.delete(s)

def ligands(selection="(all)",_self=cmd):
    cmd=_self
    try:
        s, selection, polar_contacts = _prepare(selection, _self=cmd)
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
        util.cbac("(("+s+") and not elem C)",_self=cmd)
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
            cmd.zoom(lig,3, animate=1)
        cmd.delete(host)
        cmd.delete(solvent)
        cmd.delete(near_solvent)
        cmd.delete(lig)
    except:
        traceback.print_exc()
    cmd.delete(s)

def ball_and_stick(selection="(all)",mode=1,_self=cmd):
    cmd=_self
    s, selection = _prepare(selection, _self=cmd)[:2]
    if mode == 1:
        cmd.hide("everything",s)
        cmd.set_bond("stick_color","white",s,s)
        cmd.set_bond("stick_radius","0.14",s,s)
        cmd.set("sphere_scale","0.25",s)
        cmd.show("sticks",s)
        cmd.show("spheres",s)
    elif mode == 2:
        cmd.hide("everything",s)
        cmd.set_bond("stick_color","white",s,s)
        cmd.set_bond("stick_radius","-0.14",s,s)
        cmd.set("stick_ball","1")
        cmd.set("stick_ball_ratio",-1.0)
        cmd.set("stick_ball_color","atomic")
        cmd.show("sticks",s)
    cmd.delete(s)

def b_factor_putty(selection="(name CA+P)",_self=cmd):
    cmd=_self
    s, selection = _prepare(selection, _self=cmd)[:2]
    cmd.select(s,"(name CA+P) and ("+selection+") and present")
    cmd.show("cartoon",s)
    cmd.set("cartoon_flat_sheets",0,selection)
    cmd.cartoon("putty",s)
    cmd.spectrum("b",selection=s)
    cmd.delete(s)

def ligand_cartoon(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = ligand_sites(selection, _self)[:2]
    cmd.set("cartoon_side_chain_helper",1,selection)
    cmd.show("cartoon","rep ribbon")
    cmd.hide("ribbon")
    cmd.hide("surface")
    cmd.delete(s)

def ligand_sites(selection="(all)",_self=cmd):
    cmd=_self
    try:
        s, selection, polar_contacts = _prepare(selection, _self=cmd)
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
        util.cbac("(("+s+") and not elem C)",_self=cmd)
        cmd.hide("everything",s)
        cmd.show("ribbon",host)
        cmd.show("lines","("+s+" and byres ("+host+" within 5 of "+lig+"))")
        cmd.show("surface","("+s+" and ((rep lines expand 4) within 6 of "+lig+"))")
        cmd.set("two_sided_lighting",1) # global setting
        cmd.set("transparency",0,s)
        cmd.set("surface_quality",0, selection)

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

        # add lines because nb_spheres won't show solvent with hydrogens
        cmd.show("lines", near_solvent)

        if cmd.count_atoms(lig):
            cmd.zoom(lig,3, animate=1)
        cmd.delete(host)
        cmd.delete(solvent)
        cmd.delete(near_solvent)
        cmd.delete(lig)
    finally:
        pass

    return s, selection, polar_contacts

def ligand_sites_hq(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = ligand_sites(selection, _self)[:2]
    cmd.set("surface_quality","1",selection)
    cmd.set("surface_type",0,selection)
    cmd.delete(s)

def ligand_sites_trans(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = ligand_sites(selection, _self)[:2]
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("transparency","0.33",s)
    cmd.set("surface_type",0,selection)
    cmd.set("surface_quality",0,selection)
    cmd.delete(s)

def ligand_sites_trans_hq(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = ligand_sites(selection, _self)[:2]
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("transparency","0.33",s)
    cmd.set("surface_type",0,selection)
    cmd.set("surface_quality",1,selection)
    cmd.delete(s)

def ligand_sites_mesh(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = ligand_sites(selection, _self)[:2]
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("surface_type","2",selection)
    cmd.set("surface_quality","0",selection)
    cmd.delete(s)

def ligand_sites_dots(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = ligand_sites(selection, _self)[:2]
    cmd.show("sticks",s+" and rep lines")
    cmd.hide("lines",s+" and rep lines")
    cmd.set("surface_type","1",selection)
    cmd.set("surface_quality","1",selection)
    cmd.delete(s)

def technical(selection="(all)",_self=cmd):
    cmd=_self
    s, selection, polar_contacts = _prepare(selection, _self=cmd)
    util.chainbow(s,_self=cmd)
    util.cbc("("+lig_sele+" and ("+s+"))",_self=cmd)
    util.cbac("(("+s+") and not elem C)",_self=cmd)
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
    cmd.show("nonbonded","(("+lig_sele+"|resn HOH+WAT+H2O) and ("+s+"))")
    cmd.delete(s)

def pretty_solv(selection="(all)",_self=cmd):
    pretty(selection, solv=True, _self=_self)

def pretty(selection="(all)", *, solv=False, _self=cmd):
    cmd=_self
    s, selection = _prepare(selection, _self=cmd)[:2]
    cmd.dss(s,preserve=1)
    cmd.cartoon("auto",s)
    cmd.show("cartoon",s)
    if solv:
        cmd.show("licorice", f"({lig_sele}|{wat_sele}) and ?{s}")
    else:
        cmd.show("sticks", f"({lig_sele}) and ?{s}")
    util.cbc("("+lig_sele+" and ("+s+"))",_self=cmd)
    util.cbac("("+lig_sele+" and ("+s+") and not elem C)",_self=cmd)
    cmd.spectrum("count",selection="(elem C and ("+s+") and not "+lig_sele+")")
    cmd.set("cartoon_highlight_color",-1,selection)
    cmd.set("cartoon_fancy_helices",0,selection)
    cmd.set("cartoon_smooth_loops",0,selection)
    cmd.set("cartoon_flat_sheets",1,selection)
    cmd.set("cartoon_side_chain_helper",0,selection)
    cmd.delete(s)

pretty_no_solv = pretty

def pub_solv(selection="(all)",_self=cmd):
    publication(selection, solv=True, _self=_self)

def publication(selection="(all)", *, solv=False, _self=cmd):
    cmd=_self
    pretty(selection, solv=solv, _self=_self)
    s, selection = get_sname_oname_dname(selection, _self=_self)[:2]
    cmd.set("cartoon_smooth_loops",1,selection)
    cmd.set("cartoon_highlight_color","grey50",selection)
    cmd.set("cartoon_fancy_helices",1,selection)
    cmd.set("cartoon_flat_sheets",1,selection)
    cmd.set("cartoon_side_chain_helper",0,selection)
    cmd.delete(s)

pub_no_solv = publication

def default(selection="(all)",_self=cmd):
    cmd=_self
    s, selection = _prepare(selection, _self=cmd)[:2]
    cmd.show("lines",s)
    cmd.show("nonbonded",s)
    color=cmd.get_object_color_index(selection)
    if color<0:
        util.cbag(selection,_self=cmd)
    else:
        util.cnc(selection,_self=cmd)
        cmd.color(str(color),"("+s+") and elem C")
    cmd.delete(s)


def interface(selection='*', _self=cmd):
    '''
    Protein-Protein interface preset, mimics the BioLuminate preset
    '''
    s = _prepare(selection, _self=_self)[0]

    # temporary selection names
    s_interface = _self.get_unused_name('_iface')

    # interface atoms
    _self.select(s_interface, '?%s & (%s)' % (s, ' '.join(
        '((chain "%s") around 4.5) ' % (chain)
        for chain in _self.get_chains(s))), 0)

    # Color by chain, non-carbon by element
    util.cbc(s, _self=_self)
    _self.color('atomic', '?%s & !(elem C)' % (s))

    # Change everything to cartoons
    _self.show_as('cartoon', s)

    # interface residues as sticks
    _self.show('sticks', 'byres ?' + s_interface)
    _self.show('nb_spheres', '?' + s_interface)

    # delete temporary selections
    _self.delete(s_interface)
    _self.delete(s)


def classified(selection='*', _self=cmd):
    '''
    Equivalent of "auto_show_classified" setting. Sets representations
    according to atom classification ("auto_classify_atoms"). Does not
    change any colors or settings.
    '''
    s = _prepare(selection, _self=_self)[0]

    _self.show_as('cartoon', 'polymer & %' + s)
    _self.show_as('sticks', 'organic & %' + s)
    _self.show_as('spheres', 'inorganic & %' + s)

    _self.delete(s)
