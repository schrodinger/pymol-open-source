import math

import re
import pymol
cmd = __import__("sys").modules["pymol.cmd"]
from . import setting
from . import parsing
import threading

from .cmd import DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, \
     is_list,  is_string, is_error

QuietException = parsing.QuietException

_prefix = "_tmp_editor"
tmp_wild = _prefix + "*"
tmp_editor = _prefix + "0"
tmp_connect = _prefix + "_con"
tmp_domain = _prefix + "_dom"
tmp1 = _prefix + "1"
tmp2 = _prefix + "2"
tmp3 = _prefix + "3"
tmp4 = _prefix + "4"

# routines to assist in molecular editing

class undocontext:
    def __init__(self, cmd, sele):
        # not implemented in open-source
        pass

    def __enter__(self):
        # not implemented in open-source
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        # not implemented in open-source
        pass

def attach_fragment(selection,fragment,hydrogen,anchor,*,_self=cmd):
    '''
ARGUMENTS

    selection = str: Name of a single-atom selection. If no such named
    selection exists, then create a new object with name `fragment`.

    fragment = str: fragment name to load from fragment library

    hydrogen = int: atom ID in fragment to fuse

    anchor = int: (unused)
    '''
    remove_hydrogens = _self.get_setting_boolean("auto_remove_hydrogens")

    if selection not in _self.get_names("selections"):
        if fragment in _self.get_names("objects"):
            raise pymol.CmdException("an object with that name already exists")

        _self.fragment(fragment)

        if remove_hydrogens:
            _self.remove(f"(hydro and {fragment})")
    else:
        fragment_label = _self.get_unused_name(_prefix + "_attach_fragment")
        _self.fragment(fragment, fragment_label, origin=0)

        try:
            _self.fuse(f"{fragment_label} and id {hydrogen}", f"({selection})", 1)

            if remove_hydrogens:
                _self.remove("(hydro and pkmol)")
            elif _self.count_atoms('hydro and (neighbor pk2)'):
                _self.h_fill()
        finally:
            _self.delete(fragment_label)

def combine_fragment(selection,fragment,hydrogen,anchor,_self=cmd):
    with undocontext(_self, selection):
        _self.fragment(fragment,tmp_editor)
        try:
            if _self.get_setting_boolean("auto_remove_hydrogens"):
                _self.remove("(hydro and ?%s)" % tmp_editor)
            _self.fuse("?%s" % tmp_editor, "(%s)" % selection, 3)
        finally:
            _self.delete(tmp_editor)

def attach_amino_acid(selection,amino_acid,center=0,animate=-1,object="",hydro=-1,ss=-1,_self=cmd):
    '''
ARGUMENTS

    selection = str: named selection of single N or C atom

    amino_acid = str: fragment name to load from fragment library

    center = bool: center on new terminus (pk1)

    animate = int: animate centering

    object = str: name of new object (if selection is none)

    hydro = int (-1/0/1): keep hydrogens

    ss = int: Secondary structure 1=alpha helix, 2=antiparallel beta, 3=parallel beta, 4=flat
    '''
    r = DEFAULT_SUCCESS
    ss = int(ss)
    center = int(center)

    if hydro<0:
        hydro = not int(_self.get_setting_boolean("auto_remove_hydrogens"))
    if (selection not in _self.get_names('all')
            if selection == 'pk1' # legacy, calling functions should pass '?pk1'
            else _self.count_atoms(selection) == 0):
        if object == "":
            object = amino_acid
        # create new object
        if amino_acid in _self.get_names("objects"):
            print("Error: an object with than name already exists")
            raise QuietException
        r = _self.fragment(amino_acid,object)
        if not hydro:
            _self.remove("(hydro and %s)"%object)
        if _self.count_atoms("((%s) and name C)"%object):
            _self.edit("((%s) and name C)"%object)
        elif _self.count_atoms("((%s) and name N)"%object):
            _self.edit("((%s) and name N)"%object)
    elif _self.select(tmp_connect,"(%s) & elem N,C"%selection) != 1:
        print("Error: invalid connection point: must be one atom, name N or C.")
        _self.delete(tmp_wild)
        raise QuietException
    elif amino_acid in ["nhh","nme"] and _self.select(tmp_connect,"(%s) & elem C"%selection) != 1:
        print("Error: invalid connection point: must be C for residue '%s'"%(amino_acid))
        _self.delete(tmp_wild)
        raise QuietException
    elif amino_acid in ["ace"] and _self.select(tmp_connect,"(%s) & elem N"%selection) != 1:
        print("Error: invalid connection point: must be N for residue '%s'"%(amino_acid))
        _self.delete(tmp_wild)
        raise QuietException
    else:
        if ss<0:
            ss = _self.get_setting_int("secondary_structure")
        if ss:
            if ss==1: # helix
                phi=-57.0
                psi=-47.0
            elif ss==2: # antipara-beta
                phi=-139.0
                psi=135.0
            elif ss==3: # para-beta
                phi=-119.0
                psi=113.0
            else:
                phi=180.0
                psi=180.0
        _self.fragment(amino_acid,tmp_editor, origin=0)
        if _self.count_atoms("elem N",domain=tmp_connect):
            tmp = [ None ]
            _self.iterate(tmp_connect,"tmp[0]=resv", space={ 'tmp' : tmp })
            tmp[0] = str(tmp[0]-1) # counting down
            _self.alter(tmp_editor,"resi=tmp[0]",space={ 'tmp' : tmp})
            _self.set_geometry(tmp_connect, 3, 3) # make nitrogen planar
            try:
                _self.fuse("(%s and name C)"%(tmp_editor),tmp_connect,2)
            except pymol.CmdException as e:
                _self.delete(tmp_wild)
                raise e
            _self.select(tmp_domain, "byresi (pk1 | pk2)")

            if not hydro:
                _self.remove("(pkmol and hydro)")

            if ((_self.select(tmp1,"?pk1",domain=tmp_domain)==1) and
                (_self.select(tmp2,"?pk2",domain=tmp_domain)==1)):

                if ((_self.select(tmp3,"(name CA,CH3 & nbr. ?pk1)",domain=tmp_domain)==1) and
                    (_self.select(tmp4,"(name CA,CH3 & nbr. ?pk2)",domain=tmp_domain)==1)):
                    _self.set_dihedral(tmp4,tmp2,tmp1,tmp3,180.0)

                if hydro:
                    _self.h_fix(tmp2) # fix hydrogen position

                if ss:
                    if amino_acid[0:3]!='pro':
                        if ((_self.select(tmp4,
                                          "(!(resn PRO) & name C  & nbr. (name CA & nbr. "+tmp2+"))",
                                          domain=tmp_domain)==1) and
                            (_self.select(tmp3,
                                          "(!(resn PRO) & name CA & nbr. "+tmp2+")",
                                          domain=tmp_domain)==1)):
                            _self.set_dihedral( # PHI
                                tmp4, # C
                                tmp3, # CA
                                tmp2, # N
                                tmp1, # C
                                phi)

                    if ((_self.select(tmp4,"(name N & nbr. (name CA & nbr. "+tmp1+"))",
                                      domain=tmp_domain)==1) and
                        (_self.select(tmp3,"(name CA & nbr. "+tmp1+")",domain=tmp_domain)==1)):
                        _self.set_dihedral( # PSI (n-1)
                            tmp2, # N
                            tmp1, # C
                            tmp3, # CA
                            tmp4, # N
                            psi)

            sele = ("(name N & (byres nbr. %s) &! (byres %s))"% (tmp_connect,tmp_connect))
            if _self.select(tmp1,sele,domain=tmp_domain):
                _self.edit(tmp1)
                if center:
                    _self.center(tmp1,animate=animate)
        elif _self.count_atoms("elem C",domain=tmp_connect): # forward
            tmp = [ None ]
            _self.iterate(tmp_connect,"tmp[0]=resv", space={ 'tmp' : tmp })
            tmp[0] = str(tmp[0]+1) # counting up
            _self.alter(tmp_editor,"resi=tmp[0]",space={ 'tmp' : tmp})
            _self.set_geometry(tmp_editor + " & name N", 3, 3) # make nitrogen planar
            try:
                _self.fuse("(%s and name N)"%tmp_editor,tmp_connect,2)
            except pymol.CmdException as e:
                _self.delete(tmp_wild)
                raise e

            _self.select(tmp_domain, "byresi (pk1 | pk2)")

            if not hydro:
                _self.remove("(pkmol and hydro)")

            if (( _self.select(tmp1,"?pk1",domain=tmp_domain)==1) and
                ( _self.select(tmp2,"?pk2",domain=tmp_domain)==1)):

                if ((_self.select(tmp3,"(name CA,CH3 & nbr. ?pk1)",domain=tmp_domain)==1) and
                    (_self.select(tmp4,"(name CA,CH3 & nbr. ?pk2)",domain=tmp_domain)==1)):
                    _self.set_dihedral(tmp4,tmp2,tmp1,tmp3,180.0)
                if hydro:
                    _self.h_fix("pk1") # fix hydrogen position
                if ss:
                    if hydro and amino_acid[0:3]=='nhh': # fix amide hydrogens
                        if ((_self.select(tmp3,"(name H1 & nbr. "+tmp1+")",domain=tmp_domain)==1) and
                            (_self.select(tmp4,"(name O & nbr. "+tmp2+")",domain=tmp_domain)==1)):
                            _self.set_dihedral(
                                tmp4, # O
                                tmp2, # C
                                tmp1, # N
                                tmp3, # H1
                                180)
                    if amino_acid[0:3]!='pro':
                        if ((_self.select(tmp3,"(name CA & nbr. "+tmp1+")",domain=tmp_domain)==1) and
                            (_self.select(tmp4,"(name C & nbr. (name CA & nbr. "+tmp1+"))",domain=tmp_domain)==1)):
                            _self.set_dihedral( # PHI
                                tmp2, # C
                                tmp1, # N
                                tmp3, # CA
                                tmp4, # C
                                phi)
                    if ((_self.select(tmp3,"(name CA & nbr. "+tmp2+")",domain=tmp_domain)==1) and
                        (_self.select(tmp4,"(name N & nbr. (name CA & nbr. "+tmp2+"))",domain=tmp_domain)==1)):
                        _self.set_dihedral( # PSI (n-1)
                            tmp4, # N
                            tmp3, # CA
                            tmp2, # C
                            tmp1, # N
                            psi)
            sele = ("(name C & (byres nbr. %s) & !(byres %s))"% (tmp_connect,tmp_connect))
            if _self.select(tmp1,sele,domain=tmp_domain):
                _self.edit(tmp1)
                if center:
                    _self.center(tmp1,animate=animate)
            else:
                _self.unpick()
        elif _self.count_atoms("((%s) and elem H)"%selection):
            print("Error: please pick a nitrogen or carbonyl carbon to grow from.")
            _self.delete(tmp_wild)
            raise QuietException
        else:
            print("Error: unable to attach fragment.")
            _self.delete(tmp_wild)
            raise QuietException
    _self.delete(tmp_wild)

    return r

_aa_codes =  {
    'A' : 'ala',
    'B' : 'ace',
    'C' : 'cys',
    'D' : 'asp',
    'E' : 'glu',
    'F' : 'phe',
    'G' : 'gly',
    'H' : 'his',
    'I' : 'ile',
    'K' : 'lys',
    'L' : 'leu',
    'M' : 'met',
    'N' : 'asn',
    'P' : 'pro',
    'Q' : 'gln',
    'R' : 'arg',
    'S' : 'ser',
    'T' : 'thr',
    'V' : 'val',
    'W' : 'trp',
    'Y' : 'tyr',
    'Z' : 'nme',
    }

_fab_codes = {
    'peptide' : _aa_codes,
    }

_pure_number = re.compile("[0-9]+")

def _fab(input,name,mode,resi,chain,segi,state,dir,hydro,ss,quiet,_self=cmd):
    r = DEFAULT_ERROR
    code = _fab_codes.get(mode,None)
    quiet = int(quiet)
    resi = int(resi)
    state = int(state)
    dir = int(dir)
    hydro = int(hydro)

    if hydro < 0:
        hydro = not _self.get_setting_boolean("auto_remove_hydrogens")

    seq_len = 0
    if (mode == 'peptide') and is_string(input):
        # '123/ ADC B/234/ AFCD' to [ '123/','A','D','C','B/234/','F','C','D' ]
        frags = input.split()
        input = []
        for frag in frags:
            if '/' in frag:
                input.append(frag)
            else:
                seq_len = seq_len + len(frag)
                input.extend(list(frag))
                input.append("/") # breaks chain
    if name is None:
        name = _self.get_unused_name("obj")
    elif name in _self.get_names():
        _self.delete(name)

    if mode in [ 'peptide' ]:  # polymers
        if (seq_len>99) and not quiet:
            print(" Generating a %d residue peptide from sequence..."%seq_len)
        input.reverse()
        sequence = input
        if code is not None:
            while len(sequence):
                while len(sequence) and '/' in sequence[-1]:
                    part = sequence.pop().split('/')
                    if len(part)>1:
                        if len(part[-2]):
                            resi = int(part[-2])
                    if len(part)>2:
                        chain = part[-3]
                    if len(part)>3:
                        segi = part[-4]
                if len(sequence) and not _self.count_atoms("?pk1"): # new polymer segment
                    tmp_obj = _self.get_unused_name()
                    first = sequence.pop()
                    _self.fragment(code[first], tmp_obj)
                    if not hydro:
                        _self.remove(tmp_obj + ' and hydro')
                    _self.alter(tmp_obj,'resi="""%s""";chain="""%s""";segi="""%s"""'%(resi,chain,segi))
                    _self.create(name,tmp_obj+" or ?"+name,1,1,zoom=0)
                    tmp_sel = _self.get_unused_name()
                    if mode == 'peptide':
                        if dir>0:
                            _self.select(tmp_sel,"name C and "+tmp_obj)
                            resi = resi + 1
                        else:
                            _self.select(tmp_sel,"name N and "+tmp_obj)
                            resi = resi - 1
                    _self.edit(name+" in "+tmp_sel) # set the editor's pk1 selection
                    _self.delete(tmp_sel+" "+tmp_obj)
                if mode == 'peptide':
                    while len(sequence):
                        if '/' in sequence[-1]:
                            _self.unpick() # break chain at this point
                            break
                        if not _self.count_atoms("?pk1"):
                            break
                        else:
                            attach_amino_acid("pk1",code[sequence.pop()],animate=0,ss=ss,hydro=hydro,_self=_self)
                            if dir>0:
                                resi = resi + 1
                            else:
                                resi = resi - 1
    if not len(sequence):
        r = DEFAULT_SUCCESS

    if _self.get_setting_int('auto_zoom'):
        _self.zoom(name)

    return r

_threeNA_to_OneNA = { "atp" : "A", "ctp" : "C", "gtp" : "G", "ttp" : "T", "utp" : "U"}
_oneNA_to_threeNA = { "A" : "atp", "C" : "ctp", "G" : "gtp", "T" : "ttp", "U" : "utp"}

_base_pair = { "DNA" : {"atp" : "ttp", "ctp" : "gtp", "gtp" : "ctp", "ttp" : "atp" },
               "RNA" : {"atp" : "utp", "ctp" : "gtp", "gtp" : "ctp", "utp" : "atp" }}

_oneNA_base_pair = {"DA" : "DT", "DC" : "DG", "DG" : "DC", "DT" : "DA", "A" : "U",
                    "C" : "G", "G" : "C", "U" : "A"}

def iterate_to_list(selection: str, expression: str, *, _self=cmd):
    outlist = []
    _self.iterate(selection,f"outlist.append(({expression}))", space={"outlist":outlist})
    return outlist

def rename_three_to_one(nuc_acid, sele, nuc_type, *, _self=cmd):
    """
    Renames nucleobase from 3-letter to 1-letter representation

    :param nuc_acid: (str) 3-letter nucleic acid representation
    :param sele: (str) selection of nucleic acid to rename
    :param nuc_type: (str) "DNA" or "RNA"
    """
    new_name = _threeNA_to_OneNA[nuc_acid]
    if nuc_type == "DNA":
        new_name = "D" + new_name
    _self.alter(sele, f"resn='{new_name}'")

def fit_sugars(mobile, target, *, _self=cmd):
    """
    Fits appending base pairs to form appropriate hydrogen bond

    :param mobile: (str) selection for the sense (main) strand
    :param target: (str) selection for the antisense (opposing) strand
    """
    try:
        _self.pair_fit(f"{mobile} & name C1'",
                       f"{target} & name C1'",
                       f"{mobile} & name C2'",
                       f"{target} & name C2'",
                       f"{mobile} & name C3'",
                       f"{target} & name C3'",
                       f"{mobile} & name C4'",
                       f"{target} & name C4'",
                       f"{mobile} & name O4'",
                       f"{target} & name O4'", quiet=1)
    except:
        _self.delete(tmp_wild)
        raise pymol.CmdException("Something went wrong when fitting the new residue.")

def fit_DS_fragment(mobile_A, target_A, mobile_B, target_B, *, _self=cmd):
    """
    Fits dummy fragment to the detected structure using atoms
    on both stands for a more accurate alignment.

    :param mobile_A: (str) selection for the base being created and attached
    :param target_A: (str) selection for the base selected to build on
    :param mobile_B: (str) selection for the opposing base being created
    :param target_B: (str) selection for the detected opposing base
     """
    try:
        _self.pair_fit(f"{mobile_A} & name C1'",
                       f"{target_A} & name C1'",
                       f"{mobile_A} & name C2'",
                       f"{target_A} & name C2'",
                       f"{mobile_A} & name C5'",
                       f"{target_A} & name C5'",
                       f"{mobile_A} & name O4'",
                       f"{target_A} & name O4'",
                       f"{mobile_A} & name O3'",
                       f"{target_A} & name O3'",
                       f"{mobile_A} & name P",
                       f"{target_A} & name P",
                       f"{mobile_B} & name O3'",
                       f"{target_B} & name O3'",
                       f"{mobile_B} & name C1'",
                       f"{target_B} & name C1'",
                       f"{mobile_B} & name C2'",
                       f"{target_B} & name C2'",
                       f"{mobile_B} & name C5'",
                       f"{target_B} & name C5'",
                       f"{mobile_B} & name P",
                       f"{target_B} & name P",
                       f"{mobile_B} & name O4'",
                       f"{target_B} & name O4'", quiet=1)
    except:
        _self.delete(tmp_wild)
        raise pymol.CmdException("Something went wrong when fitting the new residue.")

def add2pO(domain, nuc_acid, resv, *, _self=cmd):
    if nuc_acid == "utp": #utp comes with O2'
        return
    c_2p = _prefix + "_c2p"
    _self.select(c_2p, "%s & resi \\%i & name %s" % (domain, resv, "C2'"))
    _self.unpick()
    _self.edit(c_2p)
    _self.attach("O", 4, 4)
    _self.unpick()
    _self.alter("(byres %s) & resi \\%i & name O01" % (domain, resv), "name=\"O2'\"")

def move_atom_in_res(atom, dummy_res, new_res, twist, rise, *, _self=cmd):
    prev_coords = _self.get_coords(f"{dummy_res} & name {atom}", state=1)
    curr_coords = _self.get_coords(f"{new_res} & name {atom}", state=1)

    if curr_coords is None or prev_coords is None:
        return
    curr_coord = curr_coords[0]
    prev_coord = prev_coords[0]

    r = math.sqrt(prev_coord[0] ** 2 +
                  prev_coord[1] ** 2)
    old_phi = math.degrees(math.atan2(prev_coord[1],prev_coord[0]))
    phi = old_phi - twist
    phi = math.radians(phi)
    new_pos = [r * math.cos(phi),
               r * math.sin(phi),
               prev_coord[2] - rise]

    trans = list(new_pos - curr_coord)
    _self.translate(trans, f"{new_res} & name {atom}", camera=0)

def move_new_res(frag_string, full_frag, old, old_oppo, double_stranded_bool=False, form="B", chain="A", antisense=False, *, _self=cmd):
    """
    Attaches new residue (or pair) onto current nucleotide chain

    :param frag_string: (str) Name of appending nucleic acid or pair 
    :param full_frag: (str) Selection of the created fragment 
    :param old: (str) Selection of previous residue
    :param old_oppo: (str) Selection of previous opposing residue
    :param double_stranded_bool: (bool) Flag represing if double helix was detected
    :param antisense (bool) Flag for antisense
    :param form: (str) DNA form ('A'/'B')
    """
    if form == 'B':
        twist = -36.0
        rise = -3.375
    elif form == 'A':
        twist = -32.7
        rise = -2.548
    else:
        raise ValueError("Form not recognized")

    rise = rise if antisense else -rise
    twist = twist if antisense else -twist

    dummy_fragment = _prefix + "_dummyfrag"
    dummy_res_A = _prefix + "_dummyresA"
    dummy_res_B = _prefix + "_dummyresB"
    new_fragment = _prefix + "_newfrag"
    new_res_A = _prefix + "_newresA"
    new_res_B = _prefix + "_newresB"

    _self.select(new_fragment, f"{full_frag}")
    _self.select(new_res_A,f"{full_frag} and chain A")
    _self.select(new_res_B,f"{full_frag} and chain B")

    _self.fragment(frag_string, dummy_fragment, origin=0)
    _self.select(dummy_res_A, f"{dummy_fragment} and chain A")
    _self.select(dummy_res_B, f"{dummy_fragment} and chain B")

    if old_oppo == "none":
        # This is the case where a single residue is being added
        _self.select(dummy_res_A, f"{dummy_fragment}")
        _self.select(new_res_A, f"{full_frag}")

    atoms_A = iterate_to_list(dummy_res_A, "name")
    atoms_B = iterate_to_list(dummy_res_B, "name")

    #A new base is created by copying the coordinates of the previous
    #base and doing a cylindrical rotation (phi degrees) and a translation
    #down the z-axis by the rise amount
    for atom in atoms_A:
        move_atom_in_res(atom, dummy_res_A, new_res_A, twist, rise)
    for atom in atoms_B:
        move_atom_in_res(atom, dummy_res_B, new_res_B, twist, rise)
    if double_stranded_bool == True:
        fit_DS_fragment(dummy_res_A, old, dummy_res_B, old_oppo)
        orient_flag = check_dummy_oriention(old,dummy_res_A)
        if orient_flag == 0:
            fit_sugars(dummy_res_A,old)
    elif double_stranded_bool == False:
        fit_sugars(dummy_res_A, old)
    else:
        raise pymol.CmdException("Double stranded bool was not provided to move_new_res")

    dummy_fragment_transform = _self.get_object_matrix(dummy_fragment)
    _self.transform_object(full_frag, dummy_fragment_transform)
    _self.delete(dummy_fragment)

def check_dummy_oriention(old, dummy_res_A, *, _self=cmd):
    dummy_orient = _prefix + "_dummy_orient"
    orient_flag = 0
    orient_flag = _self.select(dummy_orient, f"({old} & name O4') within 1.0 of ({dummy_res_A} & name O4')")

    return orient_flag

class NascentNucAcidInfo:
    def __init__(self, fragment_name, nuc_acid, nuc_type, form, dbl_helix):
        self.fragment_name = fragment_name
        self.nuc_acid = nuc_acid
        self.nuc_type = nuc_type
        self.form = form
        self.dbl_helix = dbl_helix

def attach_O5_phosphate(_self=cmd):
    if "pk1" not in _self.get_names("selections"):
        raise pymol.CmdException("Selection must be pk1 to attach O5' phosphate")

    print("This building selection has an unphosphorylated O5' end.")
    attach_fragment("pk1","phosphite",4,0)
    # Initailize selection strings
    P_center = _prefix + "_P_center"
    H_extra = _prefix + "_H_extra"
    O_one = _prefix + "_O_one"
    O_two = _prefix + "_O_two"
    O_three = _prefix + "_O_three"

    # Selection
    _self.select(P_center,"n. P01")
    _self.select(H_extra, f"h. and bound_to {P_center} or n. H02")
    _self.select(O_one, f"n. O01 and bound_to {P_center}")
    _self.select(O_two, f"n. O02 and bound_to {P_center}")
    _self.select(O_three, f"n. O03 and bound_to {P_center}")

    # Removing unnecessary atoms
    _self.remove(H_extra)
    _self.remove(O_one)

    # Fix bonding
    _self.unbond(P_center,O_three)
    _self.bond(P_center,O_three,1)
    _self.unbond(P_center,O_two)
    _self.bond(P_center,O_two,2)

    # Rename P correctly
    _self.alter(P_center,"name = 'P'")

    # Set Pk1 correctly
    _self.select("pk1",P_center)

def check_DNA_base_pair(sele_oppo_atom, selection, *, _self=cmd):
    base_pair_dist = _prefix + "_base_pair_dist"
    base_pair_bool = 0
    tmp_last_resn = iterate_to_list(selection,"resn")
    tmp_last_resn_oppo = iterate_to_list(sele_oppo_atom,"resn")

    if len(tmp_last_resn_oppo) != 0:
        last_resn = str(tmp_last_resn[0])
        last_resn_oppo = str(tmp_last_resn_oppo[0])

        if (_oneNA_base_pair[last_resn] == last_resn_oppo and
            _self.select(base_pair_dist, f"(byres {sele_oppo_atom}) within 3.5 of (byres {selection})") != 0):
                base_pair_bool = 1
        else:
            base_pair_bool = 0
    else:
        print("check_DNA_base_pair has no opposing residue to check")
    return base_pair_bool

def get_chains_oppo (chain, tmp_connect, *, _self=cmd):
    models = iterate_to_list(tmp_connect,"model")

    close_chains = []
    close_chains = _self.get_chains(f"({models[0]}) within 15.0 of {tmp_connect}")
    close_chains = [c for c in close_chains if c != chain]

    return close_chains

def get_new_chain (chain, tmp_connect, *, _self=cmd):
    models = iterate_to_list(tmp_connect,"model")
    model_chains = _self.get_chains(models[0])
    search_chain_flag = 0

    if len(model_chains) != 0:
        last_chain = f"{model_chains[-1]}"
        last_chain_front = last_chain[:-1]
        last_chain_back = last_chain[-1]
        if last_chain_back != 'Z' and last_chain_back != 'z':
            new_chain_back = chr(ord(last_chain_back)+1)
        elif last_chain_back == 'Z':
            new_chain_back = "ZA"
            print("Z chain was detected. New chain will append A")
        else:
            new_chain_back = "za"
            print("z chain was detected. New chain will append a")
    else:
        new_chain_back = "A"

    chain_oppo = last_chain_front + new_chain_back
    return chain_oppo

def check_valid_attachment(nascent, atom_selection_name, selection, resv, *, _self=cmd):
    atom_selection_name_partner = "O3'" if atom_selection_name == "P" else "P"
    atom_sele = _prefix + "atom_sele"
    _self.select(atom_sele, f"{selection}")
    bound = _prefix + "atom_sele_bound"
    if _self.count_atoms(f"(bound_to {atom_sele}) & name {atom_selection_name_partner}") != 0:
        _self.delete(tmp_wild)
        raise pymol.CmdException(f"{atom_selection_name} already bonded!")

def bond_single_stranded(tmp_editor, object, chain, resv, last_resi_sele, atom_selection_name, atom_name_oppo, *, _self=cmd):
    """
    Forms a bond between the last atoms on selected structure and the newly created fragment
    :param tmp_editor: (str) Object representing the newly created fragment
    :param object: (str) object/model name of the selected structure
    :param chain: (str) Chain ID
    :param resv: (int)
    :param last_resi_sele: (str) The selection string of the selected residue
    :param atom_selection_name: (str) Name of the atom selected
    :param atom_name_oppo: (str) Name of the corresponding opposing atom
    """
    object_fuse = _prefix + f"_{chain}_fuse"
    object_connect = _prefix + f"_{chain}_con"

    print("The program did not detect a double stranded structure, so the opposing residue will not be attached.")

    # Select and fuse
    _self.select(object_fuse, f"{last_resi_sele} & name {atom_selection_name}")
    _self.fuse(f"{tmp_editor} & chain {chain} & name {atom_name_oppo}", object_fuse, mode=3)
    _self.select(object_connect, f"{last_resi_sele} & name {atom_selection_name} & chain {chain}")

    # Target is on the new fragment
    object_bond_target = _prefix + f"_{chain}_con_target"
    if (_self.select(object_bond_target, f"{object} & resi \\{resv} & name {atom_name_oppo} & chain {chain}") == 1):
        bond_dist = _prefix + "_bond_dist"
        if (_self.select(bond_dist, f"{object_connect} within 3.0 of {object_bond_target}") != 0):
            _self.bond(object_connect, object_bond_target)
        else:
            print("Identified bond targets were too far apart, so this will not be bound")
    else:
        print("More than one bond target was identified, so this will not be bound")

def bond_double_stranded(tmp_editor, object, chain, chain_oppo, resv, resv_oppo, last_resi_sele, prev_oppo_res, atom_selection_name,
                         atom_name_oppo, *, _self=cmd):
    """
    Forms a bond between the last atoms on selected structure and the newly created fragment
    :param tmp_editor: (str) Object representing the newly created fragment
    :param object: (str) object/model name of the selected structure
    :param chain: (str) Chain ID
    :param chain_oppo: (str) Opposing chain ID
    :param resv: (int)
    :param resv_oppo: (int)
    :param last_resi_sele: (str) The selection string of the selected residue
    :param prev_oppo_res: (str) The selection string of the previous opposing residue
    :param atom_selection_name: (str) Name of the atom selected
    :param atom_name_oppo: (str) Name of the corresponding opposing atom
    """
    object_fuse = _prefix + f"_{chain}_fuse"
    object_oppo_fuse = _prefix + f"_oppo_{chain_oppo}_fuse"
    object_connect = _prefix + f"_{chain}_con"
    object_oppo_connect = _prefix + f"_oppo_{chain_oppo}_con"

    # Target is on the new fragment
    object_bond_target = _prefix + f"_{chain}_con_target"
    object_oppo_bond_target = _prefix + f"_oppo_{chain_oppo}_con_target"

    # Select and fuse
    _self.select(object_fuse, f"{last_resi_sele} & name {atom_selection_name}")
    _self.select(object_oppo_fuse, f"{prev_oppo_res} & name {atom_name_oppo}")
    _self.fuse(f"{tmp_editor} & chain {chain} & name {atom_name_oppo}", object_fuse, mode=3)

    if ((_self.select(object_bond_target, f"{object} & resi \\{resv} & name {atom_name_oppo} & chain {chain}") == 1) and
        (_self.select(object_connect, f"{last_resi_sele} & name {atom_selection_name} & chain {chain}") == 1)):
        bond_dist = _prefix + "_bond_dist"
        if (_self.select(bond_dist, f"{object_connect} within 3.0 of {object_bond_target}") != 0):
            _self.bond(object_connect, object_bond_target)
        else:
            print("Identified bond targets were too far apart, so this will not be bound")
    else:
        print("More than one bond target was found on selected chain, so this will not be bound.")

    if ((_self.select(object_oppo_bond_target, f"{object} & resi \\{resv_oppo} & name {atom_selection_name} & chain {chain_oppo}") == 1) and
        (_self.select(object_oppo_connect, f"{prev_oppo_res} & name {atom_name_oppo} & chain {chain_oppo}") == 1)):
        bond_dist = _prefix + "_bond_dist"
        if (_self.select(bond_dist, f"{object_oppo_connect} within 3.0 of {object_oppo_bond_target}") != 0):
            _self.bond(object_oppo_connect, object_oppo_bond_target)
        else:
            print("Identified bond targets were too far apart, so this will not be bound")
    else:
        print("More than one bond target was found on opposing chain, so this will not be bound.")

def attach_nuc_acid(selection, nuc_acid, nuc_type, object= "", form ="B",
                        dbl_helix=True, *, _self=cmd):
    """
    Creates new nuc acid attached to existing PDB structure 
    :param selection: (str) selection of picked nascent chain (or nothing)
    :param nuc_acid: (str) appending nucleic acid
    :param nuc_type: (str) sugar type of nucleic acid
    :param object: (str) name of appending nucleobase
    :param form: (str) DNA structure form: A, B, or Z
    :param dbl_helix: (bool) flag for double-strandedness
    """
    original_sele = _prefix + "_original_sele"
    _self.select(original_sele,selection)

    if nuc_type == "RNA" and form != 'A':
        form = 'A'
        dbl_helix = False

    nascent = NascentNucAcidInfo(nuc_acid + form, nuc_acid, nuc_type, form, dbl_helix)
    nuc_acid_partner_temp = _base_pair[nuc_type][nuc_acid].lower()
    nascent_partner =  NascentNucAcidInfo(nuc_acid_partner_temp + form, nuc_acid_partner_temp,
                                          nuc_type, form, dbl_helix)

    if _self.cmd.count_atoms(selection) == 0:
        if object == "":
            object = nuc_acid

        if dbl_helix:
            frag_string = nascent.nuc_acid + "_"+ _base_pair["DNA"][nascent.nuc_acid] + nascent.form
            _self.fragment(frag_string,object)
        elif not dbl_helix:
            _self.fragment(nascent.fragment_name, object, origin=0)
            _self.alter(object, f"segi='A';chain='A';resv=1")
            rename_three_to_one(nascent.nuc_acid, object, nascent.nuc_type)

        if nascent.nuc_type == "RNA":
            add2pO(object, nascent.nuc_acid, 1)

        if _self.count_atoms(f"{object} & segi A & name P"):
            _self.edit(f"{object} & segi A & name P")
        elif _self.count_atoms(f"{object} & segi A & name O3'"):
            _self.edit(f"{object} & segi A & name O3'")
        _self.edit(f"{object} & segi A & name O3'")
        _self.select("pk1", f"{object} & name O3' & chain A")
    elif _self.select(tmp_connect, f'{selection}') == 1:
        chain, name, object = iterate_to_list(selection,"chain,name,model")[0]

        if name == "O5'":
            attach_O5_phosphate()
            name = "P"
            # FIXME Don't use `selection` as the name, or ensure that it's
            # indeed just a name and not a selection expression.
            _self.select(selection,f"name P & bound_to {original_sele}")
            _self.delete("_pkdihe")
            _self.select(tmp_connect,f"{selection}")

        if name == "P" or name == "O3'":
            extend_nuc_acid(nascent, nascent_partner, selection, object, name, chain, _self=_self)

        _self.select("pk1","_tmp_editor_new_selection")
    else:
        _self.delete(tmp_wild)
        raise pymol.CmdException("invalid connection point: must be one atom, name O3' or P")

    _self.show("cartoon", f"byobject {selection}")
    _self.delete(tmp_wild)

def extend_nuc_acid(nascent, nascent_partner, selection,
              object, atom_selection_name, chain = 'A', *, _self=cmd):
    """
    Creates new nuc acid (or pair) or attaches to PDB chain
    :param nascent: (NascentNucAcidInfo) appending nucleic acid
    :param nascent_partner: (NascentNucAcidInfo) partner of appending nucleic acid
    :param selection: (str) selection string of selected residue
    :param object: (str) object of selected residue
    :param atom_selection_name: (str) O3' or P
    :param chain: (str) chain ID

    FIXME Eliminate `tmp_connect` pre-condition (or at least document it
    properly). Looks like `tmp_connect` must be the same as `selection`?
    """
    # Making a temporary selection
    original_sele = _prefix + "_original_sele"
    _self.select(original_sele,selection)

    # Alter the segi to match the chain selection
    chain_sele = "_chain_sele"
    _self.select(chain_sele, f"chain {chain}")
    _self.alter(chain_sele,f"segi = '{chain}'")
    _self.delete(chain_sele)

    if not nascent.dbl_helix:
        frag_string = nascent.fragment_name
        _self.fragment(frag_string, tmp_editor, origin=0)
        rename_three_to_one(nascent.nuc_acid, tmp_editor, nascent.nuc_type)
    elif nascent.dbl_helix:
        frag_string = nascent.nuc_acid + "_"+ _base_pair["DNA"][nascent.nuc_acid] + nascent.form
        _self.fragment(frag_string, tmp_editor, origin=0)
    else:
        raise pymol.CmdException("No helix state selected")

    if _self.count_atoms(f"name {atom_selection_name}", domain=tmp_connect):
        tmp_resv = iterate_to_list(tmp_connect,"resv")

        # Assign last res before adjustment is made
        last_resv = int(tmp_resv[0])

        if atom_selection_name == "O3'":
            tmp_resv[0] = str(tmp_resv[0] + 1)
        elif atom_selection_name == "P":
            tmp_resv[0] = str(tmp_resv[0] - 1)
        else:
            raise pymol.CmdException("Something went wrong with resv loop in extend_nuc")

        # Resv and resi assignment and testing
        resv = int(tmp_resv[0])
        resi = resv

        check_valid_attachment(nascent, atom_selection_name, selection, resv)

        reverse = False if atom_selection_name == "O3'" else True

        last_resi_sele = _prefix + "_last_resi"
        _self.select(last_resi_sele, f"(byobject {selection}) & chain {chain} & resi \\{last_resv}")

        if not nascent.dbl_helix:
            _self.alter(tmp_editor, f"chain='{chain}';segi='{chain}';resi=tmp_resv[0]",
                        space={'tmp_resv': tmp_resv})
            # Set parameters so that move_new_res can be used
            prev_oppo_res = "none"
            double_stranded_bool = False

            # Move and fuse
            move_new_res(frag_string, tmp_editor, last_resi_sele, prev_oppo_res, double_stranded_bool, nascent.form, antisense=reverse)
            _self.fuse(f"{tmp_editor} & name P", tmp_connect, mode=3)

            _self.select(tmp_domain, "byresi (pk1 | pk2)")

            if atom_selection_name == "O3'":
                _self.bond(f"{tmp_domain} & resi \\{last_resv} & name O3'",
                           f"{tmp_domain} & resi \\{resv} & name P")
            elif atom_selection_name == "P":
                _self.bond(f"{tmp_domain} & resi \\{resv} & name O3'",
                           f"{tmp_domain} & resi \\{last_resv} & name P")

        elif nascent.dbl_helix:
            #Initialize strings for selection
            prev_oppo_res = _prefix + "_prev_oppo_res"
            end_oppo_atom = _prefix + "_end_oppo_atom"

            #Initialize a double stranded bool for the detection of existing double strand
            double_stranded_bool = False
            base_pair_result = 0

            # Get oppo information
            atom_name_oppo = "O3'" if atom_selection_name == "P" else "P"

            # Get opposing chains and iterate looking for basepairs
            chains_oppo = get_chains_oppo(chain, tmp_connect)
            for tmp_chain_oppo in chains_oppo:
                # Create selection strings
                last_oppo_atom = _prefix + "_last_oppo_atom"
                first_oppo_atom = _prefix + "_first_oppo_atom"
                same_oppo_res =  _prefix + "_same_oppo_atom"

                _self.select(last_oppo_atom, f"last (({object} and chain {tmp_chain_oppo} and (not {tmp_editor})) and polymer.nuc)")
                _self.select(first_oppo_atom, f"first (({object} and chain {tmp_chain_oppo} and (not {tmp_editor})) and polymer.nuc)")
                base_pair_last = check_DNA_base_pair(last_oppo_atom, original_sele)
                base_pair_first = check_DNA_base_pair(first_oppo_atom, original_sele)

                same_base_flag = _self.select(same_oppo_res, f"({object} and (byres {last_oppo_atom}) and (byres {first_oppo_atom}))")

                if ((base_pair_first == 1 or base_pair_last == 1) and base_pair_result == 1 and same_base_flag == 0):
                    print("Multiple residues meet base pairing requirements. Building as if no opposing strand detected.")
                    base_pair_result = 0
                    break
                elif (base_pair_first == 1 and base_pair_last ==1):
                    chain_oppo = tmp_chain_oppo
                    base_pair_result = 1

                    check_DNA_base_pair(first_oppo_atom, original_sele)
                    atoms_first = iterate_to_list(_prefix + "_base_pair_dist","name")
                    check_DNA_base_pair(last_oppo_atom, original_sele)
                    atoms_last = iterate_to_list(_prefix + "_base_pair_dist", "name")

                    if (len(atoms_first) > len(atoms_last)):
                        _self.select(end_oppo_atom,first_oppo_atom)
                    else:
                        _self.select(end_oppo_atom,last_oppo_atom)

                elif (base_pair_last == 1):
                    base_pair_result = 1
                    chain_oppo = tmp_chain_oppo
                    _self.select(end_oppo_atom,last_oppo_atom)
                elif (base_pair_first == 1):
                    base_pair_result = 1
                    chain_oppo = tmp_chain_oppo
                    _self.select(end_oppo_atom,first_oppo_atom)
                else:
                    print("No based pair was found on chain ", tmp_chain_oppo)

            if (base_pair_result == 1):
                double_stranded_bool = True
                # Confirm base pair result and set _base_pair_dist
                check_DNA_base_pair(end_oppo_atom, original_sele)

                tmp_last_resv_oppo = iterate_to_list(end_oppo_atom,"resv")
                if len(tmp_last_resv_oppo) == 1:
                    last_resv_oppo = tmp_last_resv_oppo[0]
                else:
                    # This is arbitrary and may be changed. Using neg selected resv for now.
                    last_resv_oppo = -last_resv

                _self.select(prev_oppo_res, "byres (_tmp_editor_base_pair_dist)")
            elif (base_pair_result == 0):
                last_resv_oppo = -last_resv
                chain_oppo = get_new_chain(chain, tmp_connect)
            else:
                raise pymol.CmdException("Base pairing result is not returning 0 or 1")

            # Alter the opposing segi to match chain
            chain_oppo_sele = "_chain_oppo_sele"
            _self.select(chain_oppo_sele, f"chain {chain_oppo}")
            _self.alter(chain_oppo_sele,f"segi = '{chain_oppo}'")
            _self.delete(chain_oppo_sele)

            # Get oppo resv value based on selection
            if atom_selection_name == "O3'":
                resv_oppo =  last_resv_oppo - 1
            elif atom_selection_name == "P":
                resv_oppo = last_resv_oppo + 1

            # If picking O3', check O5' phosphate based on info found in this loop
            if (atom_selection_name == "O3'" and double_stranded_bool == True):
                tmp_phosphate_check = _prefix + "_phosphate_check"

                if(_self.select(tmp_phosphate_check, f"{prev_oppo_res} and name P ") == 0):
                    _self.select("pk1", f"{prev_oppo_res} and name O5'")
                    _self.select("pk2", f"{prev_oppo_res} and name O5'") # This is so attach frag has a pk2

                    attach_O5_phosphate()

                    _self.unpick()
                    _self.select("pk1",original_sele)
                    _self.select(tmp_connect,f"{selection}")
                    _self.select(prev_oppo_res, f"byres ({prev_oppo_res})")

                    if (_self.select(tmp_phosphate_check, f"{prev_oppo_res} and name P")==1):
                        print("Phosphate has been successfully added")

            # Move the fragment
            move_new_res(frag_string, tmp_editor, last_resi_sele, prev_oppo_res, double_stranded_bool, nascent.form, antisense=reverse)

            # Alter residues created
            _self.alter("_tmp_editor_newresA",f"chain='{chain}';segi='{chain}';resv={resv}")
            _self.alter("_tmp_editor_newresB",f"chain='{chain_oppo}';segi='{chain_oppo}';resv={resv_oppo}")

            # Connect the created residues
            if double_stranded_bool == True:
                bond_double_stranded(tmp_editor, object, chain, chain_oppo, resv, resv_oppo, last_resi_sele, prev_oppo_res,
                                     atom_selection_name, atom_name_oppo)
            elif double_stranded_bool == False:
                bond_single_stranded(tmp_editor, object, chain, resv, last_resi_sele, atom_selection_name, atom_name_oppo)
            else:
                raise pymol.CmdException("double_stranded_bool is not returning True or False")

        if nascent.nuc_type == "RNA":
            add2pO(tmp_domain, nascent.nuc_acid, resv)

    new_term = _prefix + "_new_selection"
    _self.select(new_term, f"{object} & chain {chain} & resi \\{resv} & name {atom_selection_name}")
    #_self.edit(new_term)

def fab(input,name=None,mode='peptide',resi=1,chain='',segi='',state=-1,
        dir=1,hydro=-1,ss=0,async_=0,quiet=1,_self=cmd, **kwargs):
    '''
DESCRIPTION

    Build a peptide

ARGUMENTS

    input = str: sequence in one-letter code

    name = str: name of object to create {default: }

    ss = int: Secondary structure 1=alpha helix, 2=antiparallel beta, 3=parallel beta, 4=flat

EXAMPLE

    fab ACDEFGH
    fab ACDEFGH, helix, ss=1
    '''
    async_ = int(kwargs.pop('async', async_))

    if kwargs:
        raise pymol.CmdException('unknown argument: ' + ', '.join(kwargs))

    if async_ < 1:
        r = _fab(input,name,mode,resi,chain,segi,
                 state,dir,hydro,ss,quiet,_self)
    else:
        fab_thread = threading.Thread(target=_fab, args=(input,name,mode,
                                                         resi,chain,
                                                         segi,state,dir,
                                                         hydro,ss,quiet,_self))
        fab_thread.setDaemon(1)
        fab_thread.start()
        r = DEFAULT_SUCCESS
    return r

def fnab(input, name=None, mode="DNA", form="B", dbl_helix=1, *, _self=cmd):
    """
DESCRIPTION

    Builds a nucleotide acid from sequence

    Fragments provided by:
    Lu, Xiang-Jun, Olson, Wilma K. 3DNA: a software package for the analysis,
    rebuilding and visualization of three-dimensional nucleic acid structures,
    Nucleic Acids Research, 32, W667-W675 (2004).

USAGE

    fnab input [, name [, type [, form [, dbl_helix ]]]]

ARGUMENTS

    input = str: Sequence as an array of one letter codes

    name = str: Name of the object to create {default: obj}

    mode = str: "DNA" or "RNA"

    form = str: "A" or "B"

    dbl_helix = bool (0/1): flag for using double helix in DNA

EXAMPLE

    fnab ATGCGATAC
    fnab ATGCGATAC, name=myDNA, mode=DNA, form=B, dbl_helix=1
    fnab AAUUUUCCG, mode=RNA
    """
    _self.unpick()

    if name is None:
        name = _self.get_unused_name(prefix="obj")
    elif name in _self.get_names('all'):
        name = _self.get_unused_name(prefix=name)

    dbl_helix = int(dbl_helix) > 0

    mode = mode.upper()
    if mode not in ('DNA', 'RNA'):
        raise pymol.CmdException("\"mode\" must be \"DNA\" or \"RNA\" only.")
    if mode == "RNA" and dbl_helix != 0:
        print ("Double helix RNA building is not currently supported.")
        dbl_helix = 0

    #first pass for error checking
    for oneNA in input:
        oneNA = oneNA.upper()
        threeNA = _oneNA_to_threeNA.get(oneNA)
        if threeNA is None:
            raise pymol.CmdException("\"%s\" not of %s type..." % (oneNA, mode))

        if threeNA not in _base_pair[mode]:
            raise pymol.CmdException("\"%s\" not of %s type..." % (oneNA, mode))

    for oneNA in input:
        oneNA = oneNA.upper()
        threeNA = _oneNA_to_threeNA[oneNA]
        attach_nuc_acid(selection="?pk1", nuc_acid=threeNA, nuc_type=mode,
                        object=name, form=form, dbl_helix=dbl_helix)
    _self.unpick()
    return DEFAULT_SUCCESS

def build_peptide(sequence,_self=cmd): # legacy
    for aa in sequence:
        attach_amino_acid("pk1",_aa_codes[aa])
