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

import sys

cmd = __import__("sys").modules["pymol.cmd"]
import pymol
from pymol import movie
# legacy mappings, remove in PyMOL 2.0

mload = movie.load
mrock = movie.rock
mroll = movie.roll

# should match the list in layer1/Color.c:
_color_cycle = [
    26   , # /* carbon */
    5    , # /* cyan */
    154  , # /* lightmagenta */
    6    , # /* yellow */
    9    , # /* salmon */
    29   , # /* hydrogen */
    11   , # /* slate */
    13   , # /* orange */
    10   , # /* lime */
    5262 , # /* deepteal */
    12   , # /* hotpink */
    36   , # /* yelloworange */
    5271 , # /* violetpurple */
    124  , # /* grey70 */
    17   , # /* marine */
    18   , # /* olive */
    5270 , # /* smudge */
    20   , # /* teal */
    5272 , # /* dirtyviolet */
    52   , # /* wheat */
    5258 , # /* deepsalmon */
    5274 , # /* lightpink */
    5257 , # /* aquamarine */
    5256 , # /* paleyellow */
    15   , # /* limegreen */
    5277 , # /* skyblue */
    5279 , # /* warmpink */
    5276 , # /* limon */
    53   , # /* violet */
    5278 , # /* bluewhite */
    5275 , # /* greencyan */
    5269 , # /* sand */
    22   , # /* forest */
    5266 , # /* lightteal */
    5280 , # /* darksalmon */
    5267 , # /* splitpea */
    5268 , # /* raspberry */
    104  , # /* grey50 */
    23   , # /* deepblue */
    51   , # /* brown */
    ]

_color_cycle_len = len(_color_cycle)

class _verbose_cmd_proxy:
    def __init__(self, cmd):
        self.cmd = cmd
    def __getattr__(self, name):
        return getattr(self.cmd, name)
    def set(self, name, value=1):
        return self.cmd.set(name, value, quiet=0)

def color_by_area(sele, mode="molecular", state=0, palette='rainbow', _self=cmd):
    """
DESCRIPTION

    Colors molecule by surface area

ARGUMENTS

    sele = str: atom selection

    mode = str: "molecular" {default} or "solvent"
    """
    asa = 1 if mode=="solvent" else 0

    tmpObj = _self.get_unused_name("_tmp")
    tmpSel = _self.get_unused_name("_sel")
    orgSel = _self.get_unused_name("_org")

    orgN = _self.select(orgSel, sele, 0)
    _self.create(tmpObj, "byobj ?%s & ! solvent" % (orgSel), zoom=0)
    tmpN = _self.select(tmpSel, '?%s in ?%s' % (tmpObj, orgSel), 0)

    try:
        if orgN != tmpN:
            raise pymol.CmdException('color_by_area failed')

        _self.set("dot_solvent", asa, tmpObj)
        _self.set("dot_density", 3, tmpObj)

        l = []
        _self.get_area(tmpSel, load_b=1)
        _self.spectrum("b", palette, tmpSel)
        _self.iterate(tmpSel, "l_a(color)", space={'l_a': l.append})
        _self.alter(orgSel, "color=l_n()", space={'l_n': getattr(iter(l), '__next__')})

        _self.recolor(orgSel)
    finally:
        _self.delete(tmpSel)
        _self.delete(tmpObj)
        _self.delete(orgSel)

def find_surface_residues(sele, name='', _self=cmd):
    """
DESCRIPTION

    Finds those residues on the surface of a protein
    that have at least 'surface_residue_cutoff' (setting)
    exposed A**2 surface area.

    Returns the name of the selection.

ARGUMENTS

    sele = str: the object or selection in which to find exposed residues

    name = str: name of selection to create {default: exposed??}
    """
    from collections import defaultdict

    tmpObj = _self.get_unused_name("__tmp")
    selName = name or _self.get_unused_name("exposed")

    _self.select(selName, sele)
    _self.create(tmpObj, "(byobj ?%s) & ! solvent" % (selName), zoom=0)
    _self.select(selName, '?%s in ?%s' % (tmpObj, selName))

    _self.set("dot_solvent", 1, tmpObj);
    _self.get_area(selName, load_b=1)

    # threshold on what one considers an "exposed" atom (in A**2):
    surface_residue_cutoff = _self.get_setting_float("surface_residue_cutoff")

    res_area = defaultdict(int)
    _self.iterate(selName, "res_area[segi, chain, resi] += b", space=locals())

    _self.select(selName, 'none')
    _self.delete(tmpObj)

    for (k, v) in res_area.items():
        if v < surface_residue_cutoff:
            continue
        _self.select(selName, 'segi %s & chain %s & resi %s' % k, merge=1)

    return selName


def find_surface_atoms(sele, name='', cutoff=-1, _self=cmd):
    """
DESCRIPTION

    Finds those atoms on the surface of a protein
    that have at least 'cutoff' (argument) or 'surface_residue_cutoff'
    (setting) exposed A**2 surface area.

    Returns the name of the selection.

ARGUMENTS

    sele = str: the object or selection in which to find exposed residues

    name = str: name of selection to create {default: exposed??}
    """
    tmpObj = _self.get_unused_name("_tmp")
    tmpSel = _self.get_unused_name("_sel")

    _self.select(tmpSel, sele, 0)
    _self.create(tmpObj, "(byobj ?%s) & ! solvent" % (tmpSel), zoom=0)

    selName = name or _self.get_unused_name("exposed")
    _self.select(selName, '?%s in ?%s' % (tmpObj, tmpSel))

    _self.set("dot_solvent", 1, tmpObj);
    _self.get_area(selName, load_b=1)

    cutoff = float(cutoff)
    if cutoff < 0.0:
        cutoff = _self.get_setting_float("surface_residue_cutoff")

    _self.select(selName, '?%s in (?%s and b > %f)' % (tmpSel, selName, cutoff))
    _self.delete(tmpObj)
    _self.delete(tmpSel)

    return selName


def get_area(sele, state=-1, dot_solvent=0, dot_density=5, quiet=1, _self=cmd):
    '''
DESCRIPTION

    Wrapper for cmd.get_area that works on a copy of the selected object
    to set dot_solvent and dot_density.

SEE ALSO

    get_area command, dot_solvent and dot_density settings
    '''
    state, dot_density, quiet = int(state), int(dot_density), int(quiet)

    tmpSel = _self.get_unused_name("_sel")
    _self.select(tmpSel, sele, 0)

    if state < 1:
        state = _self.get_selection_state(tmpSel)

    tmpObj = _self.get_unused_name("_tmp")
    _self.create(tmpObj, "(byobj ?%s) & ! solvent" % (tmpSel), state, zoom=0)
    _self.select(tmpSel, '?%s in ?%s' % (tmpObj, tmpSel), 0)

    _self.set("dot_solvent", dot_solvent, tmpObj);
    if dot_density > -1:
        _self.set('dot_density', dot_density, tmpObj)

    r = _self.get_area(tmpSel, quiet=int(quiet))

    _self.delete(tmpSel)
    _self.delete(tmpObj)

    return r


def get_sasa(sele, state=-1, dot_density=5, quiet=1, _self=cmd):
    '''
DESCRIPTION

    Get solvent accesible surface area

SEE ALSO

    get_area command, dot_solvent and dot_density settings
    '''
    return get_area(sele, state, 1, dot_density, quiet, _self)


def mass_align(target,enabled_only=0,max_gap=50,_self=cmd):
    cmd=_self
    list = cmd.get_names("public_objects",int(enabled_only))
    [x for x in list if cmd.get_type(x)!="object:molecule"]
    if enabled_only:
        aln_object = 'aln_enabled_to'+target
    else:
        aln_object = 'aln_all_to_'+target
    cmd.delete(aln_object)
    for name in list:
        if name!=target:
            if cmd.count_atoms("(%s) and (%s)"%(target,name))==0:
                cmd.align('polymer and name CA and (%s)'%name,
                'polymer and name CA and (%s)'%target,max_gap=max_gap,quiet=0,
                          object=aln_object)

def sum_formal_charges(selection="(all)",quiet=1,_self=cmd):
    _util_sum_fc = [0]
    _self.iterate(selection, "_util_sum_fc[0] += formal_charge", space=locals())
    result = _util_sum_fc[0]
    if not int(quiet):
        print(" util.sum_formal_charges: %d" % result)
    return result

def sum_partial_charges(selection="(all)",quiet=1,_self=cmd):
    _util_sum_pc = [0.0]
    _self.iterate(selection, "_util_sum_pc[0] += partial_charge", space=locals())
    result = _util_sum_pc[0]
    if not int(quiet):
        print(" util.sum_partial_charges: sum = %0.4f"%result)
    return result

def compute_mass(selection="(all)",state=-1,implicit=False,quiet=1,_self=cmd):
    """
DESCRIPTION

    "compute_mass" calculates the atomic mass of a selection
    (in atomic mass units).
	
USAGE

    compute_mass [ selection [, state [, implicit [, quiet ]]]]

ARGUMENTS

   selection = selection, defaults to '(all)'

   state = object state, defaults to current state for each given 
           object. See notes.

   implicit = if false then only calculate masses exactly as
              in the objects; if true, then add hydrogens were
	      possible before calculating mass

EXAMPLES

  print util.compute_mass("all")

  m = util.compute_mass("organic",state=4,implicit=True)

NOTES

  If the state argument is specified and an object does not exist
  in that state, the 0 atoms will be counted for that object,
  thus resulting in a zero mass for that object.

  """
    result = 0.0
    for obj in _self.get_object_list(selection):
        if state==-1:
            state = _self.get("state",obj)
        m = _self.get_model(selection + " and " + obj,state)
        if len(m.atom)==0:
            print(" Warning: No atoms in state %d for object %s" % (state,obj))
        if implicit!=False:
            result += m.get_implicit_mass()
        else:
            result += m.get_mass()
    if not quiet:
        print(" util.compute_mass: mass = %0.4f u"%result)
    return result

def protein_assign_charges_and_radii(obj_name,_self=cmd):
    cmd=_self

    from chempy.champ import assign

    # apply a few kludges

    # convent Seleno-methionine to methionine

    cmd.alter(obj_name+"///MSE/SE","elem='S';name='SD'",quiet=1)
    cmd.alter(obj_name+"///MSE/","resn='MET'",quiet=1)
    cmd.flag("ignore",obj_name,"clear")

    # remove alternate conformers

    cmd.remove(obj_name+" and not alt ''+A")
    cmd.alter(obj_name,"alt=''")
    cmd.sort(obj_name)
    cmd.fix_chemistry(obj_name,obj_name,1)

    # make sure all atoms are included...
    cmd.alter(obj_name,"q=1.0",quiet=1)

    print(" Util: Fixing termini and assigning formal charges...")

    assign.missing_c_termini(obj_name,quiet=1,_self=_self)

    while not assign.formal_charges(obj_name,quiet=1,_self=_self):
        print(" WARNING: unrecognized or incomplete residues are being deleted:")
        cmd.iterate("(byres ("+obj_name+" and flag 23)) and flag 31",
                        'print("  "+model+"/"+segi+"/"+chain+"/"+resn+"`"+resi+"/")',quiet=1)
        cmd.remove("byres ("+obj_name+" and flag 23)") # get rid of residues that weren't assigned
        assign.missing_c_termini(obj_name,quiet=1,_self=_self)

    print(" Util: Assigning Amber 99 charges and radii...")

    cmd.h_add(obj_name)
    if not assign.amber99(obj_name,quiet=1,_self=_self):
        print(" WARNING: some unassigned atoms are being deleted:")
        cmd.iterate("byres ("+obj_name+" and flag 23)",
                        'print("  "+model+"/"+segi+"/"+chain+"/"+resn+"`"+resi+"/"+name+"? ["+elem+"]")',quiet=1)
        cmd.remove(obj_name+" and flag 23") # get rid of any atoms that weren't assigned

    # show the user what the net charges are...

    formal = sum_formal_charges(obj_name,quiet=0,_self=_self)
    partial = sum_partial_charges(obj_name,quiet=0,_self=_self)
    if round(formal)!=round(partial):
        print(" WARNING: formal and partial charge sums don't match -- there is a problem!")

def protein_vacuum_esp(selection, mode=2, border=10.0, quiet = 1, _self=cmd):
    cmd=_self

    if (selection.split() != [selection] or
         selection not in cmd.get_names('objects')):
        print(" Error: must provide an object name")
        raise cmd.QuietException
    obj_name = selection + "_e_chg"
    map_name = selection + "_e_map"
    pot_name = selection + "_e_pot"
    cmd.disable(selection)
    cmd.delete(obj_name)
    cmd.delete(map_name)
    cmd.delete(pot_name)
    cmd.create(obj_name,"((polymer and ("+selection+
               ") and (not resn A+C+T+G+U)) or ((bymol (polymer and ("+
               selection+"))) and resn NME+NHE+ACE)) and (not hydro)")
         # try to just get protein...

    protein_assign_charges_and_radii(obj_name,_self=_self)

    ext = cmd.get_extent(obj_name)
    max_length = max(abs(ext[0][0] - ext[1][0]),abs(ext[0][1] - ext[1][1]),abs(ext[0][2]-ext[1][2])) + 2*border

    # compute an grid with a maximum dimension of 50, with 10 A borders around molecule, and a 1.0 A minimum grid

    sep = max_length/50.0
    if sep<1.0: sep = 1.0
    print(" Util: Calculating electrostatic potential...")
    if mode==0: # absolute, no cutoff
        cmd.map_new(map_name,"coulomb",sep,obj_name,border)
    elif mode==1: # neutral, no cutoff
        cmd.map_new(map_name,"coulomb_neutral",sep,obj_name,border)
    else: # local, with cutoff
        cmd.map_new(map_name,"coulomb_local",sep,obj_name,border)

    cmd.ramp_new(pot_name, map_name, selection=obj_name,zero=1)
    cmd.hide("everything",obj_name)
    cmd.show("surface",obj_name)
    cmd.set("surface_color",pot_name,obj_name)
    cmd.set("surface_ramp_above_mode",1,obj_name)

def color_carbon(color,selection="(all)",_self=cmd):
    cmd=_self
    selection = str(selection)
    cmd.color(color,"(%s) and elem C"%selection)

def cbss(selection="(all)",helix_color="red",sheet_color="yellow",loop_color="green",quiet=1,_self=cmd):
    cmd=_self
    sel = str(selection)
    h = str(helix_color)
    s = str(sheet_color)
    l = str(loop_color)
    cmd.color(h,"(ss H and ("+sel+"))",quiet=quiet)
    cmd.color(s,"(ss S and ("+sel+"))",quiet=quiet)
    cmd.color(l,"((not (ss S+H)) and ("+sel+"))",quiet=quiet)

def cbag(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("carbon","(elem C and ("+s+"))",quiet=quiet)

def cbac(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("cyan","(elem C and ("+s+"))",quiet=quiet)

def cbam(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("lightmagenta","(elem C and ("+s+"))",quiet=quiet)

def cbay(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("yellow","(elem C and ("+s+"))",quiet=quiet)

def cbas(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("salmon","(elem C and ("+s+"))",quiet=quiet)

def cbaw(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("hydrogen","(elem C and ("+s+"))",quiet=quiet)

def cbab(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("slate","(elem C and ("+s+"))",quiet=quiet)

def cbao(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("brightorange","(elem C and ("+s+"))",quiet=quiet)

def cbap(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("purple","(elem C and ("+s+"))",quiet=quiet)

def cbak(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("pink","(elem C and ("+s+"))",quiet=quiet)

def cnc(selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)

def cba(color,selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color(color,"(elem C and ("+s+"))",quiet=quiet)
    cmd.color(color,s,flags=1,quiet=quiet)

def cbh(color,selection="(all)",quiet=1,_self=cmd):
    '''Wrapper around "color atomic"'''
    cmd=_self
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem H)",quiet=quiet)
    cmd.color(color,"(elem H and ("+s+"))",quiet=quiet)
    cmd.color(color,s,flags=1,quiet=quiet)

def enable_all_shaders(_self=cmd):
    """
    Turns on shaders for all representations
    """
    cmd=_self
    cmd.set("use_shaders",1)
    cmd.set("cartoon_use_shader",1)
    cmd.set("cgo_use_shader",1)
    cmd.set("dash_use_shader",1)
    cmd.set("dot_use_shader",1)
    cmd.set("line_use_shader",1)
    cmd.set("mesh_use_shader",1)
    cmd.set("nb_spheres_use_shader", 1)
    cmd.set("nonbonded_use_shader",1)
    cmd.set("ribbon_use_shader", 1)
    cmd.set("sphere_use_shader", 1)
    cmd.set("stick_use_shader",1)
    cmd.set("surface_use_shader",1)

def modernize_rendering(mode,_self=cmd):
    """
    Turns on shaders for all representations and
    updates settings to improve rendering speed
    and quality.
    """
    cmd = _verbose_cmd_proxy(_self)

    cmd.set("max_ups",0)

    enable_all_shaders(cmd)

    cmd.set("stick_ball", 0)
    cmd.set("stick_as_cylinders", 1)
    cmd.set("sphere_mode", 9)

    cmd.do("_ rebuild")

def performance(mode,_self=cmd):
    cmd = _verbose_cmd_proxy(_self)
    mode = int(mode)
    if mode==0: # maximum quality
        cmd.set('line_smooth',1)
        cmd.set('depth_cue',1)
        cmd.set('specular',1)
        cmd.set('surface_quality',1)
        cmd.set('cartoon_sampling',14)
        cmd.set('ribbon_sampling',10)
        cmd.set('transparency_mode',2)
        # new rendering
        if cmd.get_setting_int("use_shaders")==1:
            enable_all_shaders(cmd)
            # as cylinderss
            cmd.set("render_as_cylinders", 1)
            cmd.set("alignment_as_cylinders", 1)
            cmd.set("cartoon_nucleic_acid_as_cylinders", 1)
            cmd.set("dash_as_cylinders", 1)
            cmd.set("line_as_cylinders", 1)
            cmd.set("mesh_as_cylinders", 1)
            cmd.set("nonbonded_as_cylinders", 1)
            cmd.set("ribbon_as_cylinders", 1)
            cmd.set("stick_as_cylinders", 1)
            # as spheres
            cmd.set("dot_as_spheres", 1)

            # special settings
            # sticks
            cmd.set("stick_ball", 0)
            # spheres
            cmd.set("sphere_mode", 9)
            # nb_spheres
            cmd.set("nb_spheres_quality", 3)

        # old rendering
        if cmd.get_setting_int("use_shaders")==0:
            cmd.set('stick_quality',15)
            cmd.set('sphere_quality',2)
            cmd.set('stick_ball',1)
    elif mode==33:
        cmd.set('line_smooth',1)
        cmd.set('depth_cue',1)
        cmd.set('specular',1)
        cmd.set('surface_quality',0)
        cmd.set('cartoon_sampling',7)
        cmd.set('ribbon_sampling',1)
        cmd.set('transparency_mode',2)
        # new rendering
        if cmd.get_setting_int("use_shaders")==1:
            enable_all_shaders(cmd)
            # as cylinderss
            cmd.set("render_as_cylinders", 1)
            cmd.set("alignment_as_cylinders", 1)
            cmd.set("cartoon_nucleic_acid_as_cylinders", 1)
            cmd.set("dash_as_cylinders", 1)
            cmd.set("line_as_cylinders", 1)
            cmd.set("mesh_as_cylinders", 1)
            cmd.set("nonbonded_as_cylinders", 1)
            cmd.set("ribbon_as_cylinders", 1)
            cmd.set("stick_as_cylinders", 1)
            # as spheres
            cmd.set("dot_as_spheres", 1)

            # special settings
            # sticks
            cmd.set("stick_ball", 0)
            # spheres
            cmd.set("sphere_mode", 9)
            # nb_spheres
            cmd.set("nb_spheres_quality", 3)

        # old rendering
        if cmd.get_setting_int("use_shaders")==0:
            cmd.set('stick_quality',8)
            cmd.set('sphere_quality',1)
            cmd.set('stick_ball',1)

    elif mode==66: # good perfomance
        cmd.set('line_smooth',0)
        cmd.set('depth_cue',0)
        cmd.set('specular',1)
        cmd.set('surface_quality',0)
        cmd.set('cartoon_sampling',6)
        cmd.set('ribbon_sampling',1)
        cmd.set('transparency_mode',2)
        # new rendering
        if cmd.get_setting_int("use_shaders")==1:
            enable_all_shaders(cmd)
            # as cylinderss
            cmd.set("render_as_cylinders", 1)
            cmd.set("alignment_as_cylinders", 0)
            cmd.set("cartoon_nucleic_acid_as_cylinders", 0)
            cmd.set("dash_as_cylinders", 0)
            cmd.set("line_as_cylinders", 0)
            cmd.set("mesh_as_cylinders", 0)
            cmd.set("nonbonded_as_cylinders", 0)
            cmd.set("ribbon_as_cylinders", 0)
            cmd.set("stick_as_cylinders", 1)
            # as spheres
            cmd.set("dot_as_spheres", 0)

            # special settings
            # sticks
            cmd.set("stick_ball", 0)
            # spheres
            cmd.set("sphere_mode", 9)
            # nb_spheres
            cmd.set("nb_spheres_quality", 3)
        # old rendering
        if cmd.get_setting_int("use_shaders")==0:
            cmd.set('stick_quality',8)
            cmd.set('sphere_quality',1)
            cmd.set('stick_ball',0)

    else: # maximum performance
        cmd.set('line_smooth',0)
        cmd.set('depth_cue',0)
        cmd.set('specular',0)
        cmd.set('surface_quality',-1) # new
        cmd.set('stick_quality',5)
        cmd.set('sphere_quality',0)
        cmd.set('ribbon_sampling',1)
        cmd.set('cartoon_sampling',3)
        cmd.set('transparency_mode',0)
        cmd.set('max_ups',0)
        # new rendering
        if cmd.get_setting_int("use_shaders")==1:
            enable_all_shaders(cmd)
            # as cylinderss
            cmd.set("render_as_cylinders", 1)
            cmd.set("alignment_as_cylinders", 0)
            cmd.set("cartoon_nucleic_acid_as_cylinders", 0)
            cmd.set("dash_as_cylinders", 0)
            cmd.set("line_as_cylinders", 0)
            cmd.set("mesh_as_cylinders", 0)
            cmd.set("nonbonded_as_cylinders", 0)
            cmd.set("ribbon_as_cylinders", 0)
            cmd.set("stick_as_cylinders", 1)
            # as spheres
            cmd.set("dot_as_spheres", 0)
        # old rendering
        if cmd.get_setting_int("use_shaders")==0:
            cmd.set('stick_quality',5)
            cmd.set('sphere_quality',0)
            cmd.set('stick_ball',0)

    cmd.do("rebuild")


def label_chains(sele="all",_self=cmd):
    pymol=_self._pymol
    cmd=_self
    pymol.stored._cs = []
    last = None
    save = ()
    list = []
    cmd.iterate(sele,"stored._cs.append((model,chain,index))")
    for a in pymol.stored._cs:
        if (a[0:2]!=save):
            list.append(last)
            list.append(a)
            save = a[0:2]
        last = a
    if len(list):
        list.append(last)
    list = [_f for _f in list if _f]
    for a in list:
        if(a[1]==''):
            cmd.label("%s`%d"%(a[0],a[2]),'''"chain ''"''',quiet=1)
        elif(a[1]==' '):
            cmd.label("%s`%d"%(a[0],a[2]),'''"chain ' '"''',quiet=1)
        else:
            cmd.label("%s`%d"%(a[0],a[2]),"'chain '+chain",quiet=1)

def label_segments(sele="all",_self=cmd):
    pymol=_self._pymol
    cmd=_self
    pymol.stored._cs = []
    last = None
    save = ()
    list = []
    cmd.iterate(sele,"stored._cs.append((model,segi,index))")
    for a in pymol.stored._cs:
        if (a[0:2]!=save):
            list.append(last)
            list.append(a)
            save = a[0:2]
        last = a
    if len(list):
        list.append(last)
    list = [_f for _f in list if _f]
    for a in list:
        if(a[1]==''):
            cmd.label("%s`%d"%(a[0],a[2]),'''"segi ''"''',quiet=1)
        elif(a[1]==' '):
            cmd.label("%s`%d"%(a[0],a[2]),'''"segi ' '"''',quiet=1)
        else:
            cmd.label("%s`%d"%(a[0],a[2]),"'segi '+segi",quiet=1)

def cbc(selection='(all)',first_color=7,quiet=1,legacy=0,_self=cmd):
    '''
    Color all chains a different color
    '''
    legacy = int(legacy)
    for c, a in enumerate(_self.get_chains(selection)):
        if len(a.split()) != 1:
            a = '"' + a + '"'
        if legacy:
            color = c + int(first_color)
        else:
            color = _color_cycle[c % _color_cycle_len]
        if not quiet:
            print(" util.cbc: color %s, (chain %s)" % (color, a))
        _self.color(color, "(chain %s and (%s))" % (a, selection), quiet=quiet)

def color_objs(selection='all', quiet=1, _self=cmd):
    cmd=_self
    '''
    Color all objects a different color
    '''
    c = 0
    for a in cmd.get_names('public_nongroup_objects',selection=selection):
        if (selection!='all') and (selection!='(all)'):
            cmd.set_object_color(a, _color_cycle[c])
            cmd.color(_color_cycle[c],"(?%s and (%s))"%(a,selection),quiet=quiet)
        else:
            cmd.color(_color_cycle[c], a, quiet=quiet)
        c = (c + 1) % _color_cycle_len

def color_deep(color, name='all', quiet=1, _self=cmd):
    '''
    Unset all object and atom level (not global) color settings and
    apply given color.
    '''
    print(' util.color_deep: Deprecated, use cmd.color_deep() instead')
    _self.color_deep(color, name, quiet)

def chainbow(selection='(all)', palette="rainbow", quiet=1, _self=cmd):
    '''
    Color all chains in rainbow
    '''
    for model in _self.get_object_list('(' + selection + ')'):
        for a in _self.get_chains('model %s & (%s)' % (model, selection)) or ['']:
            _self.spectrum('count', palette,
                    "(chain '%s' & model %s & (%s))" % (a, model, selection),
                    byres=1, quiet=quiet)

color_chains = cbc

def ray_shadows(mode,_self=cmd):
    cmd = _verbose_cmd_proxy(_self)
    # adjustment factors for new lighting model in 0.99
    if mode=='none':
        cmd.set('ray_shadow',0)
    else:
        cmd.set('ray_shadow',1)
    if mode=='light':
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.65)
        cmd.set('reflect',0.25)
        cmd.set('shininess',40)
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('spec_count',-1)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='matte':
        cmd.set('light_count',4)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set("light2","[-0.3,-0.4,-1.0]")
        cmd.set("light3","[-0.3,-0.3,-1.0]")
        cmd.set("light4","[-0.4,-0.3,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.45)
        cmd.set('reflect',0.45)
        cmd.set('shininess',25)
        cmd.set('power',1.25)
        cmd.set('spec_count',-1)
        cmd.set('specular_intensity',0.2)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='soft':
        cmd.set('light_count',10)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set("light2","[-0.3,-0.4,-1.0]")
        cmd.set("light3","[-0.3,-0.3,-1.0]")
        cmd.set("light4","[-0.4,-0.3,-1.0]")
        cmd.set("light5","[-0.4,-0.5,-1.0]")
        cmd.set("light6","[-0.5,-0.5,-1.0]")
        cmd.set("light7","[-0.5,-0.4,-1.0]")
        cmd.set("light8","[-0.5,-0.3,-1.0]")
        cmd.set("light9","[-0.3,-0.5,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.40)
        cmd.set('reflect',0.50)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_count',1)
        cmd.set('shininess',55)
        cmd.set('power',1.0)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0.1)
        cmd.set('ray_shadow_decay_range',1.8)
    elif mode=='occlusion':
        cmd.set('light_count',9)
        cmd.set("light" ,"[-0.2,-0.2,-1.0]")
        cmd.set("light2","[-0.2, 0.0,-1.0]")
        cmd.set("light3","[-0.2, 0.2,-1.0]")
        cmd.set("light4","[ 0.0, 0.2,-1.0]")
        cmd.set("light5","[ 0.2, 0.2,-1.0]")
        cmd.set("light6","[ 0.2, 0.0,-1.0]")
        cmd.set("light7","[ 0.2,-0.2,-1.0]")
        cmd.set("light8","[ 0.0,-0.2,-1.0]")
        cmd.set('ambient',0.18)
        cmd.set('direct',0.10)
        cmd.set('reflect',0.80)
        cmd.set('shininess',10)
        cmd.set('spec_count',-1)
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0)
        cmd.set('spec_direct',0.25)
        cmd.set('ray_shadow_decay_factor',0.1)
        cmd.set('ray_shadow_decay_range',5.0)
    elif mode=="occlusion2":
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.45)
        cmd.set('reflect',0.45)
        cmd.set('shininess',55)
        cmd.set('spec_count',-1)
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
        cmd.set('ambient_occlusion_mode', 1)
    elif mode=='medium':
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.45)
        cmd.set('reflect',0.45)
        cmd.set('shininess',55)
        cmd.set('spec_count',-1)
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='heavy':
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.05)
        cmd.set('direct',0.20)
        cmd.set('reflect',0.85)
        cmd.set('spec_count',-1)
        cmd.set('shininess',90)
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='black': # best for light backgrounds
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.001)
        cmd.set('direct',0.0)
        cmd.set('reflect',1.1)
        cmd.set('spec_count',-1)
        cmd.set('power',1.0)
        cmd.set('shininess',90)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)

def ff_copy(src,dst,_self=cmd):
    pymol=_self._pymol
    cmd=_self # NOT THREAD SAFE
    pymol._rcopy = pymol.Scratch_Storage()
    pymol._rcopy.pc={}
    pymol._rcopy.tt={}
    cmd.iterate("(%s)"%src,"_rcopy.pc[name]=partial_charge")
    cmd.alter("(%s)"%dst,"partial_charge=_rcopy.pc[name]")
    cmd.iterate("(%s)"%src,"_rcopy.tt[name]=text_type")
    cmd.alter("(%s)"%dst,"text_type=_rcopy.tt[name]")
    del pymol._rcopy

def b2vdw(selection='all', _self=cmd):
    # use B values to create RMS VDW spheres
    # rms = sqrt(b/(8*(PI^2)))
    _self.alter(selection, "vdw=(b/78.9568352087)**0.5")

def phipsi(selection="(pk1)",_self=cmd):
    cmd=_self # NOT THREAD SAFE
    n_sele =   "((byres (%s)) & name N)"%selection
    c_sele =   "((byres (%s)) & name C)"%selection
    ca_sele =  "((byres (%s)) & name CA)"%selection
    cm_sele = "((neighbor (%s)) and not (byres (%s)))"%(n_sele,n_sele)
    np_sele = "((neighbor (%s)) and not (byres (%s)))"%(c_sele,c_sele)
    cmd.feedback("push")
    cmd.feedback("disable","selector","everythin")
    cm_cnt = cmd.select("_pp_cm",cm_sele)
    n_cnt = cmd.select("_pp_n",n_sele)
    c_cnt = cmd.select("_pp_c",c_sele)
    ca_cnt = cmd.select("_pp_ca",ca_sele)
    np_cnt = cmd.select("_pp_np",np_sele)
    if(cm_cnt and n_cnt and ca_cnt and c_cnt):
        phi = cmd.get_dihedral("_pp_c","_pp_ca","_pp_n","_pp_cm")
    else:
        phi = None
    if(n_cnt and ca_cnt and c_cnt and np_cnt):
        psi = cmd.get_dihedral("_pp_np","_pp_c","_pp_ca","_pp_n")
    else:
        psi = None
    cmd.feedback("pop")
    cmd.delete("_pp_cm")
    cmd.delete("_pp_n")
    cmd.delete("_pp_c")
    cmd.delete("_pp_ca")
    cmd.delete("_pp_np")
    return (phi,psi)

def rainbow(selection="(name CA and alt ''+A)",reverse=0,_self=cmd):
    '''
    Legacy spectrum coloring routine. Don't use.

    Use instead: spectrum
    '''
    print(' util.rainbow: Deprecated, use cmd.spectrum() instead')
    list = [
        (0,0,255),
        (0,0,255),
        (0,128,255),
        (0,255,255),
        (0,255,128),
        (0,255,0),
        (128,255,0),
        (255,255,0),
        (255,128,0),
        (255,0,0),
        (255,0,0)
        ]
    if reverse:
        list.reverse()
    list = ['0x%02x%02x%02x' % rgb for rgb in list]
    _self.spectrum('count', ' '.join(list), selection)


def ss(selection="(name CA and alt '',A)",state=1,_self=cmd):
    '''
    Legacy secondary structure assignment routine. Don't use.

    Use instead: dss
    '''
    print(' util.ss: Deprecated, use cmd.dss() instead')
    _self.dss(selection, state)


def colors(scheme="",_self=cmd):
    cmd=_self
    if scheme=="jmol":
        cmd.set("auto_color",0)
        cmd.set_color("hydrogen",[1.000,1.000,1.000])
        cmd.set_color("carbon",[0.567,0.567,0.567])
        cmd.set_color("nitrogen",[0.189,0.315,0.976])
        cmd.set_color("oxygen",[1.000,0.051,0.051])
        cmd.set_color("fluorine",[0.567,0.882,0.314])
        cmd.set_color("sulfur",[1.000,1.000,0.189])
        cmd.color("carbon","elem C")
        cmd.recolor()


def interchain_distances(name, selection, cutoff=None, mode=None, label=0, reset=1, _self=cmd):
    '''
DESCRIPTION

    Find distances between all chains in selection
    '''
    import itertools

    if int(reset):
        _self.distance(name, 'none', 'none', 1.0, reset=1)

    chains = _self.get_chains(selection)

    for c1, c2 in itertools.combinations(chains, 2):
        s1 = '(%s) & chain "%s"' % (selection, c1)
        s2 = '(%s) & chain "%s"' % (selection, c2)
        _self.distance(name, s1, s2, cutoff, mode, label=label)

    _self.enable(name)


def get_sasa_relative(selection='all', state=1, vis=-1, var='b', quiet=1, outfile='', _self=cmd):
    '''
DESCRIPTION

    Calculates the relative per-residue solvent accessible surface area
    and optionally labels and colors residues. The value is relative to
    full exposure of the residue, calculated by removing all other
    residues except its two next neighbors, if present.

    Loads a value beteween 0.0 (fully buried) and 1.0 (fully exposed)
    into the b-factor property, available in "iterate", "alter" and
    "label" as "b".

USAGE

    get_sasa_relative [ selection [, state [, vis [, var ]]]]

ARGUMENTS

    selection = str: atom selection {default: all}

    state = int: object state {default: 1}

    vis = 0/1: show labels and do color by exposure {default: !quiet}

    var = str: name of property to assign {default: b}

    quiet = 0/1: print results to log window

    outfile = str: filename, write to file instead of log window {default: }

EXAMPLE

    fetch 1ubq, async=0
    get_sasa_relative polymer

PYTHON API

    cmd.get_sasa_relative(...) -> dict

SEE ALSO

    get_area with "load_b=1" argument.
    '''
    import collections

    state, vis, quiet = int(state), int(vis), int(quiet)
    if vis == -1:
        vis = not quiet

    sele = _self.get_unused_name('_sele')
    tripepname = _self.get_unused_name('_tripep')

    _self.select(sele, selection, 0)

    dot_solvent = _self.get_setting_boolean('dot_solvent')
    _self.set('dot_solvent', 1, updates=0)

    try:
        for model in _self.get_object_list(sele):
            _self.get_area('?%s & ?%s' % (sele, model), state, load_b=1)

        resarea = collections.defaultdict(float)
        _self.iterate(sele, 'resarea[model,segi,chain,resi] += b', space=locals())

        for key in resarea:
            _self.create(tripepname, 'byres (/%s/%s/%s/`%s extend 1)' % key, state, 1, zoom=0)
            _self.get_area(tripepname, 1, load_b=1)

            resarea_exposed = [0.0]
            _self.iterate('/' + tripepname + '/%s/%s/`%s' % key[1:],
                    'resarea_exposed[0] += b', space=locals())
            _self.delete(tripepname)

            if resarea_exposed[0] == 0.0:
                continue

            resarea[key] /= resarea_exposed[0]

        _self.alter(sele,
                var + ' = resarea[model,segi,chain,resi]', space=locals())

        handle = open(outfile, 'w') if outfile else None

        def callback(key, area):
            w = len(key[0]) + 15
            per10 = int(10 * area)
            s = ('/%s/%s/%s/%s`%s ' % key).ljust(w)
            s += '%3.0f' % (area * 100) + '% '
            s += '|' + '=' * per10 + ' ' * (10 - per10) + '|'
            print(s, file=handle)

        if not quiet or outfile:
            _self.iterate('?%s & guide' % sele,
                    'callback((model, segi, chain, resn, resi), ' + var + ')',
                    space=locals())

        if outfile:
            handle.close()
            print(" Written results to %s" % (outfile))

        if vis:
            _self.label('?%s & guide' % (sele), '"%.1f" % ' + var)
            _self.spectrum(var, 'white blue', sele, minimum=0., maximum=1.)

    finally:
        _self.set('dot_solvent', dot_solvent, updates=0)
        _self.delete(sele)

    return resarea


def ligand_zoom(step=1, _self=cmd):
    '''
DESCRIPTION

    Zoom to the next organic molecule (by residue identifier)

    Bound to CTRL-L key.
    '''
    global _current_ligand

    s = {'ligand_set': set()}
    if _self.iterate('organic', 'ligand_set.add((model,segi,chain,resi))',
            space=s) < 1:
        return

    ligands = sorted(s["ligand_set"])

    try:
        i = ligands.index(_current_ligand)
    except (ValueError, NameError):
        i = -1

    i = (i + int(step)) % len(ligands)
    _current_ligand = ligands[i]

    # use "do" for feedback
    _self.do('zoom /%s/%s/%s & resi %s, animate=1, buffer=2' % ligands[i])
