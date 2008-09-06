import pymol
import cmd
import setting
import parsing

QuietException = parsing.QuietException

tmp_editor = "_tmp_editor"
tmp_ed_save = "_tmp_ed_save"
tmp1 = "_tmp_editor1"
tmp2 = "_tmp_editor2"
tmp3 = "_tmp_editor3"
tmp4 = "_tmp_editor4"

# routines to assist in molecular editing

def attach_fragment(selection,fragment,hydrogen,anchor,_self=cmd):
    if not selection in _self.get_names("selections"):
        if fragment in _self.get_names("objects"):
            print " Error: an object with than name already exists"
            raise QuietException
        else:
            _self.fragment(fragment)
            if _self.get_setting_legacy("auto_remove_hydrogens"):
                _self.remove("(hydro and %s)"%fragment)
    else:
        _self.fragment(fragment,tmp_editor)
        if _self.count_atoms("((%s) and elem h)"%selection,quiet=1):
            _self.fuse("(%s and id %d)"%(tmp_editor,hydrogen),"(pk1)",1)
            if _self.get_setting_legacy("auto_remove_hydrogens"):
                _self.remove("(hydro and pkmol)")            
        else:
            _self.remove("(%s and id %d)"%(tmp_editor,hydrogen))
            _self.fuse("(%s and id %d)"%(tmp_editor,anchor),"(pk1)",1)
            if _self.get_setting_legacy("auto_remove_hydrogens"):
                _self.remove("(hydro and pkmol)")            
        _self.delete(tmp_editor)
        
def attach_amino_acid(selection,amino_acid,center=0,animate=-1,object="",_self=cmd):
    if (selection not in _self.get_names('all')):
        if object == "":
            object = amino_acid
        # create new object 
        if amino_acid in _self.get_names("objects"):
            print " Error: an object with than name already exists"
            raise QuietException
        _self.fragment(amino_acid,object)
        if _self.get_setting_legacy("auto_remove_hydrogens"):
            _self.remove("(hydro and %s)"%object)
        if _self.count_atoms("((%s) and name c)"%object,quiet=1):
            _self.edit("((%s) and name c)"%object)
        elif _self.count_atoms("((%s) and name n)"%object,quiet=1):
            _self.edit("((%s) and name n)"%object)
    else:
        ss = int(_self.get_setting_legacy("secondary_structure"))
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
        _self.fragment(amino_acid,tmp_editor)
        if _self.count_atoms("((%s) and elem n)"%selection,quiet=1):
            _self.select(tmp_ed_save,"(%s)"%selection)
            _self.iterate("(%s)"%selection,"stored.resv=resv")
            pymol.stored.resi = str(pymol.stored.resv-1)
            _self.alter(tmp_editor,"resi=stored.resi")
            _self.fuse("(%s and name C)"%(tmp_editor),tmp_ed_save,2)
            if _self.get_setting_legacy("auto_remove_hydrogens"):
                _self.remove("(pkmol and hydro)")
            if ((_self.count_atoms("(name ca,ch3 and neighbor ?pk1)")==1) and
                (_self.count_atoms("(name ca and neighbor ?pk2)"))):
                _self.set_dihedral("(name ca and neighbor pk2)",
                                 "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
            if (_self.select(tmp1,"?pk1")==1) and (_self.select(tmp2,"?pk2")==1):
            
                if 0:
                    if ((_self.select(tmp3,"(name ca and neighbor "+tmp2+"))")>0) and
                        (_self.select(tmp4,"(name ca and neighbor "+tmp1+")")>0)):
                        _self.set_dihedral( # PHI
                            tmp3, # CA +0
                            tmp2, # N +0 
                            tmp1, # C -1
                            tmp4, # CA -1
                            180.0) # insure that the peptide is planer

                _self.set_geometry(tmp2,3,3) # make nitrogen planer
    
                if ss:
                    if amino_acid[0:3]!='pro':
                        if ((_self.select(tmp4,
                                        "((!r;pro) and name c and neighbor (name ca and neighbor "
                                        +tmp2+"))")==1) and
                            (_self.select(tmp3,"((!r;pro) and name ca and neighbor "+tmp2+")")==1)):
                            _self.set_dihedral( # PHI
                                tmp4, # C
                                tmp3, # CA 
                                tmp2, # N
                                tmp1, # C
                                phi)

                    if ((_self.select(tmp4,"(name n and neighbor (name ca and neighbor "+tmp1+"))")==1) and
                        (_self.select(tmp3,"(name ca and neighbor "+tmp1+")")==1)):
                        _self.set_dihedral( # PSI (n-1)
                            tmp2, # N
                            tmp1, # C
                            tmp3, # CA
                            tmp4, # N
                            psi)
                _self.h_fix(tmp2) # fix hydrogen position
                        
            _self.delete(tmp1)
            _self.delete(tmp2)
            _self.delete(tmp3)
            _self.delete(tmp4)
                
            sele = ("(name N and (byres neighbor %s) and not (byres %s))"%
                      (tmp_ed_save,tmp_ed_save))
            if _self.count_atoms(sele,quiet=1):
                _self.edit(sele)
                _self.center(sele,animate=animate)
            _self.delete(tmp_ed_save)
                    
        elif _self.count_atoms("((%s) and elem c)"%selection,quiet=1):
            _self.select(tmp_ed_save,"(%s)"%selection)
            _self.iterate("(%s)"%selection,"stored.resv=resv")
            pymol.stored.resi = str(pymol.stored.resv+1)
            _self.alter(tmp_editor,"resi=stored.resi")
            _self.fuse("(%s and name N)"%(tmp_editor),tmp_ed_save,2)
            if _self.get_setting_legacy("auto_remove_hydrogens"):
                _self.remove("(pkmol and hydro)")
            if (_self.count_atoms("(name ca,ch3 and neighbor ?pk1)") and
                _self.count_atoms("(name ca,ch3 and neighbor ?pk2)")):
                _self.set_dihedral("(name ca,ch3 and neighbor pk2)",
                                 "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
            if "pk1" in _self.get_names('selections'):
                _self.set_geometry("pk1",3,3) # make nitrogen planer
                _self.h_fix("pk1") # fix hydrogen position
            if ss:
                if (_self.select(tmp1,"?pk1")==1) and (_self.select(tmp2,"?pk2")==1):
                    if amino_acid[0:3]=='nhh': # fix amide hydrogens
                        if ((_self.select(tmp3,"(name h1 and neighbor "+tmp1+")")==1) and
                            (_self.select(tmp4,"(name o and neighbor "+tmp2+")")==1)):
                            _self.set_dihedral(
                                tmp4, # O
                                tmp2, # C
                                tmp1, # N
                                tmp3, # H1
                                180)
                    if amino_acid[0:3]!='pro':
                        if ((_self.select(tmp3,"(name ca and neighbor "+tmp1+")")==1) and
                            (_self.select(tmp4,"(name c and neighbor (name ca and neighbor "+tmp1+"))")==1)):
                            _self.set_dihedral( # PHI
                                tmp2, # C
                                tmp1, # N
                                tmp3, # CA 
                                tmp4, # C
                                phi)
                    if ((_self.select(tmp3,"(name ca and neighbor "+tmp2+")")==1) and
                        (_self.select(tmp4,"(name n and neighbor (name ca and neighbor "+tmp2+"))")==1)):
                        _self.set_dihedral( # PSI (n-1)
                            tmp4, # N
                            tmp3, # CA
                            tmp2, # C
                            tmp1, # N
                            psi)
                _self.delete(tmp1)
                _self.delete(tmp2)
                _self.delete(tmp3)
                _self.delete(tmp4)                               
            sele = ("(name C and (byres neighbor %s) and not (byres %s))"%
                      (tmp_ed_save,tmp_ed_save))
            if _self.count_atoms(sele,quiet=1):
                _self.edit(sele)
                _self.center(sele,animate=animate)                
            _self.delete(tmp_ed_save)
        elif _self.count_atoms("((%s) and elem h)"%selection,quiet=1):
            print " Error: please pick a nitrogen or carbonyl carbon to grow from."
            _self.delete(tmp_editor)
            raise QuietException
        else:
            print " Error: unable to attach fragment."
    _self.delete(tmp_editor)

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

def build_peptide(sequence):
    for aa in sequence:
        attach_amino_acid("pk1",_aa_codes[aa])
        
