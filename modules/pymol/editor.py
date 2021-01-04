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
            _self.fuse("(%s and name C)"%(tmp_editor),tmp_connect,2)
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
            _self.fuse("(%s and name N)"%tmp_editor,tmp_connect,2)
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

def build_peptide(sequence,_self=cmd): # legacy
    for aa in sequence:
        attach_amino_acid("pk1",_aa_codes[aa])
